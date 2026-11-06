use super::{Writer, invalid};
use crate::Unit;
use crate::tess::indexed::{self, Format};
use std::{io::{self, Write}, num::NonZeroUsize, ops::Range};

#[derive(Default, Debug)]
pub struct Index(indexed::Index);
impl Index {
    pub fn save(&self, path: &str) -> io::Result<()> { self.0.save(path) }
}

pub struct BgzfWriter<W: Write + Send + 'static> { writer: Writer<indexed::Writer<W>> }
impl<W: Write + Send + 'static> BgzfWriter<W> {
    pub fn new(writer: W, workers: NonZeroUsize) -> io::Result<Self> {
        Ok(Self {writer: Writer::new(indexed::Writer::new(writer, workers, Format::Tq))?})
    }
    pub fn push(&mut self, name: &[u8], units: &[u16]) -> io::Result<()> { self.writer.push(name, units) }
    pub fn finish(self) -> io::Result<(W, Index)> {
        let (writer, records, decoded) = self.writer.into_parts();
        writer.finish(records, decoded).map(|(w, i)| (w, Index(i)))
    }
}
#[derive(Clone, Copy, Debug)]
pub struct RecordInfo {
    namespace: u64,
    offset: u64,
    len: usize,
}
impl RecordInfo {
    pub fn len(&self) -> usize {
        self.len
    }
}

pub struct Store {
    bytes: indexed::Store,
    view_cache: crate::tess::ReferenceCache,
}
impl Store {
    pub fn open(path: &str) -> io::Result<Self> {
        Ok(Self {bytes: indexed::Store::open(path, Format::Tq)?, view_cache: Default::default()})
    }
    fn read_at(&self, at: u64, out: &mut [u8]) -> io::Result<()> { self.bytes.read_at(at, out) }
    fn header(&self, at: u64) -> io::Result<(usize, RecordInfo)> {
        let mut header = [0; 12];
        self.read_at(at, &mut header)?;
        let name_len = u32::from_le_bytes(header[..4].try_into().unwrap()) as usize;
        let len = u64::from_le_bytes(header[4..].try_into().unwrap());
        let offset = at
            .checked_add(12 + name_len as u64)
            .ok_or_else(|| invalid("TQ record overflow"))?;
        if len > isize::MAX as u64 / 2
            || offset
                .checked_add(len * 2)
                .is_none_or(|end| end > self.bytes.len())
        {
            return Err(invalid("truncated TQ record"));
        }
        Ok((
            name_len,
            RecordInfo {
                namespace: self.bytes.namespace,
                offset,
                len: len as usize,
            },
        ))
    }

    pub fn n_records(&self) -> usize {
        self.bytes.index.records.len()
    }

    pub(crate) fn view_cache(&self) -> crate::tess::ReferenceCache {
        self.view_cache
    }

    pub fn record(&self, rec: usize) -> io::Result<(Vec<u8>, RecordInfo)> {
        let &at = self
            .bytes.index
            .records
            .get(rec)
            .ok_or_else(|| invalid("TQ record out of bounds"))?;
        let (name_len, info) = self.header(at)?;
        let end = self
            .bytes.index
            .records
            .get(rec + 1)
            .copied()
            .unwrap_or(self.bytes.len());
        if info.offset + info.len as u64 * 2 != end {
            return Err(invalid("TQ record and index disagree"));
        }
        let mut name = Vec::new();
        let mut buffer = [0; 8192];
        while name.len() < name_len {
            let n = (name_len - name.len()).min(buffer.len());
            self.read_at(at + 12 + name.len() as u64, &mut buffer[..n])?;
            name.extend_from_slice(&buffer[..n]);
        }
        Ok((name, info))
    }

    pub fn units(&self, info: RecordInfo, range: Range<usize>) -> io::Result<Vec<u16>> {
        if info.namespace != self.bytes.namespace {
            return Err(invalid("TQ record belongs to another store"));
        }
        if range.start > range.end || range.end > info.len {
            return Err(invalid("TQ tile range out of bounds"));
        }
        let mut units = Vec::new();
        let mut buffer = [0; 8192];
        while units.len() < range.len() {
            let n = (range.len() - units.len()).min(buffer.len() / 2);
            self.read_at(
                info.offset + (range.start + units.len()) as u64 * 2,
                &mut buffer[..2 * n],
            )?;
            units.reserve(n);
            for bytes in buffer[..2 * n].chunks_exact(2) {
                let unit = u16::from_le_bytes([bytes[0], bytes[1]]);
                if Unit::canonical(unit).is_none() {
                    return Err(invalid("invalid canonical unit ID in TQ"));
                }
                units.push(unit);
            }
        }
        Ok(units)
    }
}
