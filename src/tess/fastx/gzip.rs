use std::io::{self, BufRead, Cursor, Read};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::thread::JoinHandle;

const CHUNK_BYTES: usize = 4 << 20;
const CHUNKS: usize = 8;

pub(super) struct Reader {
    ready: Option<Receiver<io::Result<Vec<u8>>>>,
    recycle: Option<SyncSender<Vec<u8>>>,
    chunk: Cursor<Vec<u8>>,
    worker: Option<JoinHandle<()>>,
    eof: bool,
}

impl Reader {
    pub(super) fn new(input: impl BufRead + Send + 'static) -> io::Result<Self> {
        let (send, ready) = sync_channel(CHUNKS);
        let (recycle, free) = sync_channel(CHUNKS);
        for _ in 0..CHUNKS {
            recycle.send(vec![0; CHUNK_BYTES]).unwrap();
        }
        let worker = std::thread::Builder::new()
            .name("fastx-gzip".into())
            .spawn(move || {
                let mut decoder = flate2::bufread::MultiGzDecoder::new(input);
                while let Ok(mut chunk) = free.recv() {
                    chunk.resize(CHUNK_BYTES, 0);
                    match fill_chunk(&mut decoder, &mut chunk) {
                        Ok(0) => break,
                        Ok(size) => {
                            chunk.truncate(size);
                            if send.send(Ok(chunk)).is_err() {
                                break;
                            }
                        }
                        Err(error) => {
                            let _ = send.send(Err(error));
                            break;
                        }
                    }
                }
            })?;
        Ok(Self {
            ready: Some(ready),
            recycle: Some(recycle),
            chunk: Cursor::new(Vec::new()),
            worker: Some(worker),
            eof: false,
        })
    }
}

fn fill_chunk(reader: &mut impl Read, chunk: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < chunk.len() {
        match reader.read(&mut chunk[filled..]) {
            Ok(0) => break,
            Ok(size) => filled += size,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error),
        }
    }
    Ok(filled)
}

impl Read for Reader {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        if output.is_empty() {
            return Ok(0);
        }
        let size = self.chunk.read(output)?;
        if size > 0 || self.eof {
            return Ok(size);
        }
        if !self.chunk.get_ref().is_empty() {
            let chunk = std::mem::take(self.chunk.get_mut());
            let _ = self.recycle.as_ref().unwrap().send(chunk);
        }
        match self.ready.as_ref().unwrap().recv() {
            Ok(Ok(chunk)) => {
                self.chunk = Cursor::new(chunk);
                self.chunk.read(output)
            }
            Ok(Err(error)) => {
                self.eof = true;
                Err(error)
            }
            Err(_) => {
                self.eof = true;
                self.worker
                    .take()
                    .unwrap()
                    .join()
                    .map_err(|_| io::Error::other("gzip decoder thread panicked"))?;
                Ok(0)
            }
        }
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        self.ready.take();
        self.recycle.take();
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}