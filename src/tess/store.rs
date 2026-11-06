use super::packed::{Packed, PackedView};
use super::indexed::{self, Bytes, Format};
use super::codec::invalid;
use std::io::{self, Write};
use std::num::NonZeroUsize;

pub const MAGIC: [u8; 8] = *b"\x89TSA\r\n\x1a\n";
pub struct Tess {
    bytes: indexed::Store,
    pub offsets: Vec<u64>,
}
impl Tess {
    pub fn open(stem: &str) -> Self {
        Self::try_open(stem).unwrap_or_else(|e| panic!("cannot open {stem}.tss: {e}"))
    }
    pub fn try_open(stem: &str) -> io::Result<Self> {
        let mut bytes = indexed::Store::open(&format!("{stem}.tss"), Format::Tss)?;
        let offsets = std::mem::take(&mut bytes.index.records);
        Ok(Self {bytes, offsets})
    }
    fn source_at(&self, offset: usize) -> io::Result<Bytes<'_>> {
        let offset = offset as u64;
        let rec = self.offsets.binary_search(&offset).map_err(|_| invalid("no record at byte offset"))?;
        let end = self.offsets.get(rec + 1).copied().unwrap_or(self.bytes.len());
        let mut header = [0; 16];
        self.bytes.read_at(offset, &mut header)?;
        if header[..8] != MAGIC { return Err(invalid("invalid TSS record header")); }
        let len = u64::from_le_bytes(header[8..].try_into().unwrap());
        let start = offset.checked_add(16).ok_or_else(|| invalid("record offset overflow"))?;
        if start.checked_add(len) != Some(end) { return Err(invalid("TSS record and index disagree")); }
        let len = usize::try_from(len).map_err(|_| invalid("record exceeds address space"))?;
        Ok(Bytes::Indexed {store: &self.bytes, offset: start, len})
    }
    pub fn view_at(&self, offset: usize) -> PackedView<'_> {
        self.try_view_at(offset).unwrap_or_else(|e| panic!("invalid Packed record at byte {offset}: {e}"))
    }
    pub fn try_view_at(&self, offset: usize) -> io::Result<PackedView<'_>> {
        PackedView::from_source(self.source_at(offset)?)
    }
}

pub struct Index(indexed::Index);
impl Index {
    pub fn save(&self, path: &str) -> io::Result<()> { self.0.save(path) }
}
pub struct BgzfWriter<W: Write + Send + 'static> { writer: Writer<indexed::Writer<W>> }
impl<W: Write + Send + 'static> BgzfWriter<W> {
    pub fn new(writer: W, workers: NonZeroUsize) -> Self {
        Self {writer: Writer::new(indexed::Writer::new(writer, workers, Format::Tss))}
    }
    pub fn push(&mut self, record: &Packed) -> io::Result<u64> { self.writer.push(record) }
    pub fn finish(self) -> io::Result<(W, Index)> {
        let (writer, records, decoded) = self.writer.into_parts();
        writer.finish(records, decoded).map(|(w, i)| (w, Index(i)))
    }
}
pub struct Writer<W: Write> {
    writer: W,
    position: u64,
    records: Vec<u64>,
}
impl<W: Write> Writer<W> {
    pub fn new(writer: W) -> Self {
        Self {
            writer,
            position: 0,
            records: Vec::new(),
        }
    }
    pub fn push(&mut self, record: &Packed) -> io::Result<u64> {
        self.push_bytes(record.as_bytes())
    }
    pub fn push_bytes(&mut self, bytes: &[u8]) -> io::Result<u64> {
        let offset = self.position;
        let len = bytes.len() as u64;
        let next = offset
            .checked_add(16)
            .and_then(|n| n.checked_add(len))
            .ok_or_else(|| io::Error::other("auxiliary store length overflow"))?;
        self.writer.write_all(&MAGIC)?;
        self.writer.write_all(&len.to_le_bytes())?;
        self.writer.write_all(bytes)?;
        self.records.push(offset);
        self.position = next;
        Ok(offset)
    }
    pub fn into_parts(self) -> (W, Vec<u64>, u64) {
        (self.writer, self.records, self.position)
    }
}
