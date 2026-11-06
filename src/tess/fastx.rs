use std::io::{self, BufRead, BufReader, Read};
use std::num::NonZero;
mod gzip;

pub type Record = (Vec<u8>, Vec<u8>);
pub const BATCH_BASES: usize = 1 << 24;
const PREFETCH_BASES: usize = 1 << 26;
const MAX_INPUTS: usize = 16;

fn is_bgzf(header: &[u8]) -> bool {
    if header.len() < 12 || header[..3] != [0x1f, 0x8b, 8] || header[3] & 4 == 0 {
        return false;
    }
    let end = 12 + u16::from_le_bytes([header[10], header[11]]) as usize;
    let Some(extra) = header.get(12..end) else {
        return false;
    };
    let mut at = 0;
    while at + 4 <= extra.len() {
        let size = u16::from_le_bytes([extra[at + 2], extra[at + 3]]) as usize;
        if at + 4 + size > extra.len() {
            return false;
        }
        if &extra[at..at + 2] == b"BC" && size == 2 {
            return true;
        }
        at += 4 + size;
    }
    false
}

pub(crate) fn open(path: &str, workers: usize) -> io::Result<Box<dyn Read + Send>> {
    let input: Box<dyn Read + Send> = if path == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(std::fs::File::open(path)?)
    };
    let mut input = BufReader::with_capacity(1 << 22, input);
    let header = input.fill_buf()?;
    let bgzf = is_bgzf(header);
    let gzip = header.starts_with(&[0x1f, 0x8b]);
    if bgzf {
        let workers = NonZero::new(workers).unwrap();
        Ok(Box::new(
            noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, input),
        ))
    } else if gzip {
        Ok(Box::new(gzip::Reader::new(input)?))
    } else {
        Ok(Box::new(input))
    }
}

fn read(
    path: &str,
    workers: usize,
    batch_bases: usize,
    process: &mut dyn FnMut(usize, Vec<Record>) -> io::Result<()>,
) -> io::Result<()> {
    let mut reader =
        needletail::parse_fastx_reader(open(path, workers)?).map_err(io::Error::other)?;
    let (mut batch, mut bases, mut first) = (Vec::new(), 0usize, 0usize);
    while let Some(record) = reader.next() {
        let record = record.map_err(io::Error::other)?;
        let name = record
            .id()
            .split(|b| b.is_ascii_whitespace())
            .next()
            .unwrap_or_default()
            .to_vec();
        let sequence = record.seq().into_owned();
        bases += sequence.len();
        batch.push((name, sequence));
        if bases >= batch_bases {
            let batch = std::mem::take(&mut batch);
            let n = batch.len();
            process(first, batch)?;
            (first, bases) = (first + n, 0);
        }
    }
    if !batch.is_empty() {
        process(first, batch)?;
    }
    Ok(())
}

fn receive_result<T>(
    receiver: &std::sync::mpsc::Receiver<T>,
) -> Result<T, std::sync::mpsc::RecvError> {
    if rayon::current_thread_index().is_none() {
        return receiver.recv();
    }
    loop {
        match receiver.try_recv() {
            Ok(value) => return Ok(value),
            Err(std::sync::mpsc::TryRecvError::Disconnected) => {
                return Err(std::sync::mpsc::RecvError);
            }
            Err(std::sync::mpsc::TryRecvError::Empty) => {
                rayon::yield_now();
                std::thread::yield_now();
            }
        }
    }
}

pub fn parallel<P: AsRef<str> + Sync, T: Send>(
    paths: &[P],
    map: impl Fn(usize, Vec<Record>) -> T + Sync,
    process: impl FnMut(usize, T) -> io::Result<()>,
) -> io::Result<()> {
    parallel_records(paths, read, map, process)
}

pub(crate) type BatchReader<R> =
    fn(&str, usize, usize, &mut dyn FnMut(usize, Vec<R>) -> io::Result<()>) -> io::Result<()>;

pub(crate) fn parallel_records<P: AsRef<str> + Sync, R: Send, T: Send>(
    paths: &[P],
    read: BatchReader<R>,
    map: impl Fn(usize, Vec<R>) -> T + Sync,
    process: impl FnMut(usize, T) -> io::Result<()>,
) -> io::Result<()> {
    parallel_with(paths, BATCH_BASES, PREFETCH_BASES, read, map, process)
}

fn parallel_with<P: AsRef<str> + Sync, R: Send, T: Send>(
    paths: &[P],
    batch_bases: usize,
    prefetch_bases: usize,
    read: BatchReader<R>,
    map: impl Fn(usize, Vec<R>) -> T + Sync,
    mut process: impl FnMut(usize, T) -> io::Result<()>,
) -> io::Result<()> {
    use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

    if paths.is_empty() {
        return Ok(());
    }
    let count = rayon::current_num_threads();
    let prefetch = (prefetch_bases / batch_bases).max(1);
    let pending = count.saturating_sub(1).saturating_add(prefetch);
    let active = paths.len().min(count).min(MAX_INPUTS);
    let input_workers = (count / active).max(1);
    let next = AtomicUsize::new(0);
    let (send, receive) = std::sync::mpsc::sync_channel(pending);
    let map = &map;
    rayon::in_place_scope(|workers| {
        std::thread::scope(|threads| {
            let receive = receive;
            for _ in 0..active {
                let send = send.clone();
                let next = &next;
                threads.spawn(move || {
                    loop {
                        let input = next.fetch_add(1, Relaxed);
                        let Some(path) = paths.get(input).map(AsRef::as_ref) else {
                            break;
                        };
                        if let Err(error) =
                            read(path, input_workers, batch_bases, &mut |first, batch| {
                                let (result, ready) = std::sync::mpsc::sync_channel(1);
                                send.send((input, Ok(ready))).map_err(|_| {
                                    io::Error::new(
                                        io::ErrorKind::BrokenPipe,
                                        "FASTX consumer stopped",
                                    )
                                })?;
                                workers.spawn(move |_| {
                                    let _ = result.send(map(first, batch));
                                });
                                Ok(())
                            })
                        {
                            let error = io::Error::new(error.kind(), format!("{path}: {error}"));
                            let _ = send.send((input, Err(error)));
                            break;
                        }
                    }
                });
            }
            drop(send);
            while let Ok((input, ready)) = receive_result(&receive) {
                process(
                    input,
                    receive_result(&ready?).expect("FASTX worker panicked"),
                )?;
            }
            Ok(())
        })
    })
}