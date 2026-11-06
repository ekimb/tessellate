use crate::Unit;
use std::io::{self, Read, Write};

mod index;
pub use index::{BgzfWriter, RecordInfo, Store};

pub const MAGIC: [u8; 8] = *b"\x89TQ\x01\r\n\x1a\n";
const CHUNK_BYTES: usize = 8192;
pub type Record = (Vec<u8>, Vec<u16>);

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

pub fn stem(path: &str) -> Option<&str> {
    let lower = path.to_ascii_lowercase();
    let end = [".gz", ".bgz", ".bgzf"]
        .iter()
        .find(|suffix| lower.ends_with(*suffix))
        .map_or(lower.len(), |suffix| lower.len() - suffix.len());
    lower[..end].ends_with(".tq").then(|| &path[..end - 3])
}

pub struct Reader<R> {
    reader: R,
}
impl<R: Read> Reader<R> {
    pub fn new(mut reader: R) -> io::Result<Self> {
        let mut magic = [0; 8];
        match reader.read_exact(&mut magic) {
            Ok(()) if magic == MAGIC => Ok(Self { reader }),
            Ok(()) => Err(invalid(
                "unsupported TQ format; regenerate with tessellate fastx",
            )),
            Err(error) if error.kind() == io::ErrorKind::UnexpectedEof => Err(invalid(
                "missing or truncated binary TQ header; regenerate with tessellate fastx",
            )),
            Err(error) => Err(error),
        }
    }

    pub fn read_record(&mut self) -> io::Result<Option<Record>> {
        let mut header = [0; 12];
        match self.reader.read_exact(&mut header[..1]) {
            Ok(()) => {}
            Err(error) if error.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(error) => return Err(error),
        }
        self.reader.read_exact(&mut header[1..])?;
        let name_len = u32::from_le_bytes(header[..4].try_into().unwrap()) as usize;
        let count = usize::try_from(u64::from_le_bytes(header[4..].try_into().unwrap()))
            .ok()
            .filter(|&n| n <= isize::MAX as usize / 2)
            .ok_or_else(|| invalid("TQ unit count exceeds address space"))?;
        let mut buffer = [0; CHUNK_BYTES];
        let mut name = Vec::new();
        while name.len() < name_len {
            let n = (name_len - name.len()).min(buffer.len());
            self.reader.read_exact(&mut buffer[..n])?;
            name.extend_from_slice(&buffer[..n]);
        }
        let mut units = Vec::new();
        while units.len() < count {
            let n = (count - units.len()).min(buffer.len() / 2);
            self.reader.read_exact(&mut buffer[..2 * n])?;
            units.reserve(n);
            for bytes in buffer[..2 * n].chunks_exact(2) {
                let unit = u16::from_le_bytes([bytes[0], bytes[1]]);
                if Unit::canonical(unit).is_none() {
                    return Err(invalid("invalid canonical unit ID in TQ"));
                }
                units.push(unit);
            }
        }
        Ok(Some((name, units)))
    }
}

pub struct Writer<W> {
    writer: W,
    position: u64,
    records: Vec<u64>,
}
impl<W: Write> Writer<W> {
    pub fn new(mut writer: W) -> io::Result<Self> {
        writer.write_all(&MAGIC)?;
        Ok(Self {
            writer,
            position: MAGIC.len() as u64,
            records: Vec::new(),
        })
    }

    pub fn push(&mut self, name: &[u8], units: &[u16]) -> io::Result<()> {
        let name_len = u32::try_from(name.len()).map_err(|_| invalid("TQ record name too long"))?;
        if units.iter().any(|&unit| Unit::canonical(unit).is_none()) {
            return Err(invalid("invalid canonical unit ID in TQ"));
        }
        let next = self
            .position
            .checked_add(12 + name.len() as u64)
            .and_then(|n| n.checked_add((units.len() as u64).checked_mul(2)?))
            .ok_or_else(|| invalid("TQ length overflow"))?;
        self.writer.write_all(&name_len.to_le_bytes())?;
        self.writer.write_all(&(units.len() as u64).to_le_bytes())?;
        self.writer.write_all(name)?;
        let mut buffer = [0; CHUNK_BYTES];
        for chunk in units.chunks(buffer.len() / 2) {
            for (bytes, unit) in buffer.chunks_exact_mut(2).zip(chunk) {
                bytes.copy_from_slice(&unit.to_le_bytes());
            }
            self.writer.write_all(&buffer[..chunk.len() * 2])?;
        }
        self.records.push(self.position);
        self.position = next;
        Ok(())
    }

    pub fn into_parts(self) -> (W, Vec<u64>, u64) {
        (self.writer, self.records, self.position)
    }
}

fn read(
    path: &str,
    workers: usize,
    batch_bytes: usize,
    process: &mut dyn FnMut(usize, Vec<Record>) -> io::Result<()>,
) -> io::Result<()> {
    let mut reader = Reader::new(super::fastx::open(path, workers)?)?;
    let (mut batch, mut bytes, mut first) = (Vec::new(), 0usize, 0usize);
    while let Some(record) = reader.read_record()? {
        bytes += 12 + record.0.len() + 2 * record.1.len();
        batch.push(record);
        if bytes >= batch_bytes {
            let batch = std::mem::take(&mut batch);
            let n = batch.len();
            process(first, batch)?;
            (first, bytes) = (first + n, 0);
        }
    }
    if !batch.is_empty() {
        process(first, batch)?;
    }
    Ok(())
}

pub fn parallel<P: AsRef<str> + Sync, T: Send>(
    paths: &[P],
    map: impl Fn(usize, Vec<Record>) -> T + Sync,
    process: impl FnMut(usize, T) -> io::Result<()>,
) -> io::Result<()> {
    super::fastx::parallel_records(paths, read, map, process)
}
