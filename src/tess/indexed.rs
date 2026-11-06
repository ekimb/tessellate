use super::codec::invalid;
use mmap_rs::{Mmap, MmapOptions};
use std::{
    cell::RefCell,
    collections::{HashMap, VecDeque},
    fs::{self, File, Metadata, OpenOptions},
    io::{self, Read, Write},
    num::NonZeroUsize,
    ops::Range,
    rc::Rc,
    sync::atomic::{AtomicU64, Ordering},
    time::UNIX_EPOCH,
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum Format { #[default] Tq, Tss }
impl Format {
    fn magic(self) -> &'static [u8; 8] {
        match self { Self::Tq => &super::tq::MAGIC, Self::Tss => &super::store::MAGIC }
    }
    fn index_magic(self) -> &'static [u8; 8] {
        match self { Self::Tq => b"\x89TQI\x01\r\n\x1a", Self::Tss => b"\x89TSI\x01\r\n\x1a" }
    }
    pub(crate) fn path(self, path: &str) -> String {
        format!("{path}.{}", if self == Self::Tq {"tqi"} else {"tsi"})
    }
    fn first(self) -> u64 { if self == Self::Tq {8} else {0} }
    fn label(self) -> &'static str { if self == Self::Tq {"TQ"} else {"TSS"} }
}
const CACHE_BYTES: usize = 8 << 20;
static NEXT: AtomicU64 = AtomicU64::new(1);
#[derive(Default, Debug)]
pub struct Index {
    format: Format,
    bgzf: bool,
    source_len: u64,
    decoded_len: u64,
    pub(crate) records: Vec<u64>,
    blocks: Vec<(u64, u64)>,
}

fn stamp(meta: &Metadata) -> io::Result<(u64, u64)> {
    let time = meta
        .modified()?
        .duration_since(UNIX_EPOCH)
        .map_err(io::Error::other)?;
    Ok((time.as_secs(), u64::from(time.subsec_nanos())))
}

impl Index {
    pub fn save(&self, path: &str) -> io::Result<()> {
        let meta = fs::metadata(path)?;
        if meta.len() != self.source_len {
            return Err(invalid("store changed while indexing"));
        }
        let (sec, nano) = stamp(&meta)?;
        let dest = self.format.path(path);
        let temp = format!(
            "{dest}.{}.{}.tmp",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        );
        let result = (|| {
            let file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temp)?;
            let mut out = io::BufWriter::new(file);
            out.write_all(self.format.index_magic())?;
            for value in [
                u64::from(self.bgzf),
                self.source_len,
                self.decoded_len,
                sec,
                nano,
                self.records.len() as u64,
                self.blocks.len() as u64,
            ] {
                out.write_all(&value.to_le_bytes())?;
            }
            for &offset in &self.records {
                out.write_all(&offset.to_le_bytes())?;
            }
            for &(compressed, decoded) in &self.blocks {
                out.write_all(&compressed.to_le_bytes())?;
                out.write_all(&decoded.to_le_bytes())?;
            }
            out.flush()?;
            out.get_ref().sync_data()?;
            fs::rename(&temp, &dest)
        })();
        if result.is_err() {
            let _ = fs::remove_file(&temp);
        }
        result
    }

    fn load(path: &str, meta: &Metadata, format: Format) -> io::Result<Self> {
        let bytes = fs::read(format.path(path))?;
        if bytes.len() < 64 || &bytes[..8] != format.index_magic() {
            return Err(invalid("invalid store index header"));
        }
        let word = |at| u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap());
        let (sec, nano) = stamp(meta)?;
        if word(8) > 1 || word(16) != meta.len() || (word(32), word(40)) != (sec, nano) {
            return Err(invalid("stale store index"));
        }
        let (records, blocks) = (word(48), word(56));
        let expected = records
            .checked_mul(8)
            .and_then(|n| n.checked_add(blocks.checked_mul(16)?))
            .and_then(|n| n.checked_add(64));
        if expected != Some(bytes.len() as u64) {
            return Err(invalid("truncated store index"));
        }
        let start = 64 + records as usize * 8;
        let index = Self {
            format,
            bgzf: word(8) == 1,
            source_len: word(16),
            decoded_len: word(24),
            records: (64..start).step_by(8).map(word).collect(),
            blocks: (start..bytes.len())
                .step_by(16)
                .map(|at| (word(at), word(at + 8)))
                .collect(),
        };
        index.validate()?;
        Ok(index)
    }

    fn validate(&self) -> io::Result<()> {
        let records_ok = self.decoded_len >= self.format.first()
            && self.records.first().copied().unwrap_or(self.decoded_len) == self.format.first()
            && self.records.iter().all(|&p| p < self.decoded_len)
            && self.records.windows(2).all(|w| w[0] < w[1]);
        let blocks_ok = if self.bgzf && self.decoded_len == 0 {
            self.blocks.is_empty() && self.source_len > 0
        } else if self.bgzf {
            self.blocks.first() == Some(&(0, 0))
                && self
                    .blocks
                    .iter()
                    .all(|&(c, u)| c < self.source_len && u < self.decoded_len)
                && self
                    .blocks
                    .windows(2)
                    .all(|w| w[0].0 < w[1].0 && w[0].1 < w[1].1 && w[1].1 - w[0].1 <= 65536)
                && self.decoded_len - self.blocks.last().unwrap().1 <= 65536
        } else {
            self.blocks.is_empty() && self.source_len == self.decoded_len
        };
        if !records_ok || !blocks_ok {
            return Err(invalid("invalid store index offsets"));
        }
        Ok(())
    }
}
fn block_size(header: &[u8]) -> io::Result<usize> {
    if header.len() < 18
        || header[..4] != [31, 139, 8, 4]
        || header[10..16] != [6, 0, b'B', b'C', 2, 0]
    {
        return Err(invalid(
            "indexed stores require BGZF or uncompressed data; regenerate with tessellate fastx",
        ));
    }
    let size = usize::from(u16::from_le_bytes(header[16..18].try_into().unwrap())) + 1;
    if size < 26 {
        return Err(invalid("invalid BGZF block length"));
    }
    Ok(size)
}
struct BlockWriter<W> {
    writer: W,
    frame: Vec<u8>,
    index: Index,
    eof: bool,
}
impl<W: Write> Write for BlockWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let n = self.writer.write(bytes)?;
        let mut rest = &bytes[..n];
        while !rest.is_empty() {
            if self.eof {
                return Err(invalid("data after BGZF EOF block"));
            }
            let target = if self.frame.len() < 18 {
                18
            } else {
                block_size(&self.frame)?
            };
            let take = (target - self.frame.len()).min(rest.len());
            self.frame.extend_from_slice(&rest[..take]);
            rest = &rest[take..];
            if self.frame.len() >= 18 && self.frame.len() == block_size(&self.frame)? {
                let size = self.frame.len();
                let decoded = u32::from_le_bytes(self.frame[size - 4..].try_into().unwrap()) as u64;
                if decoded > 65536 {
                    return Err(invalid("invalid BGZF uncompressed size"));
                }
                if decoded > 0 {
                    self.index
                        .blocks
                        .push((self.index.source_len, self.index.decoded_len));
                    self.index.decoded_len += decoded;
                } else {
                    self.eof = true;
                }
                self.index.source_len += size as u64;
                self.frame.clear();
            }
        }
        Ok(n)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}
pub(crate) struct Writer<W: Write + Send + 'static> {
    writer: noodles_bgzf::io::MultithreadedWriter<BlockWriter<W>>,
}
impl<W: Write + Send + 'static> Writer<W> {
    pub(crate) fn new(writer: W, workers: NonZeroUsize, format: Format) -> Self {
        let blocks = BlockWriter {
            writer, frame: Vec::new(),
            index: Index { format, bgzf: true, ..Default::default() }, eof: false,
        };
        Self { writer: noodles_bgzf::io::MultithreadedWriter::with_worker_count(workers, blocks) }
    }
    pub(crate) fn finish(mut self, records: Vec<u64>, decoded: u64) -> io::Result<(W, Index)> {
        let mut blocks = self.writer.finish()?;
        blocks.flush()?;
        if !blocks.frame.is_empty() || !blocks.eof || blocks.index.decoded_len != decoded {
            return Err(invalid("incomplete BGZF output"));
        }
        blocks.index.records = records;
        blocks.index.validate()?;
        Ok((blocks.writer, blocks.index))
    }
}
impl<W: Write + Send + 'static> Write for Writer<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> { self.writer.write(bytes) }
    fn flush(&mut self) -> io::Result<()> { self.writer.flush() }
}

#[derive(Default)]
struct Cache {
    blocks: HashMap<(u64, usize), Rc<Vec<u8>>>,
    order: VecDeque<(u64, usize)>,
    bytes: usize,
}
impl Cache {
    fn insert(&mut self, key: (u64, usize), bytes: Rc<Vec<u8>>) {
        if self.blocks.contains_key(&key) {
            return;
        }
        while self.bytes + bytes.len() > CACHE_BYTES {
            let old = self.order.pop_front().unwrap();
            self.bytes -= self.blocks.remove(&old).unwrap().len();
        }
        self.bytes += bytes.len();
        self.order.push_back(key);
        self.blocks.insert(key, bytes);
    }
}
thread_local! { static CACHE: RefCell<Cache> = RefCell::new(Cache::default()); }
pub(crate) struct Store {
    map: Option<Mmap>,
    pub(crate) index: Index,
    pub(crate) namespace: u64,
}
impl Store {
    pub(crate) fn open(path: &str, format: Format) -> io::Result<Self> {
        let file = File::open(path)?;
        let meta = file.metadata()?;
        let len = usize::try_from(meta.len()).map_err(|_| invalid("store exceeds address space"))?;
        if len == 0 && format == Format::Tq {
            return Err(invalid("empty TQ file"));
        }
        let map = if len == 0 {None} else {Some(unsafe {
            MmapOptions::new(len)
                .map_err(io::Error::other)?
                .with_file(&file, 0)
                .map()
                .map_err(io::Error::other)?
        })};
        let loaded = Index::load(path, &meta, format);
        let mut store = Self {
            map,
            index: Index { format, ..Default::default() },
            namespace: NEXT.fetch_add(1, Ordering::Relaxed),
        };
        match loaded {
            Ok(index) => store.index = index,
            Err(error)
                if matches!(
                    error.kind(),
                    io::ErrorKind::NotFound | io::ErrorKind::InvalidData
                ) =>
            {
                eprintln!("[tessellate {}]\tindexing {path}", format.label());
                store.build_index()?;
                let after = file.metadata()?;
                if after.len() != meta.len() || stamp(&after)? != stamp(&meta)? {
                    return Err(invalid("store changed while indexing"));
                }
                if let Err(error) = store.index.save(path) {
                    if !matches!(error.kind(), io::ErrorKind::PermissionDenied | io::ErrorKind::ReadOnlyFilesystem) {
                        return Err(error);
                    }
                    eprintln!("[tessellate {}]\tusing an in-memory index for {path}: {error}", format.label());
                }
            }
            Err(error) => return Err(error),
        }
        store.check_magic()?;
        Ok(store)
    }

    fn bytes(&self) -> &[u8] { self.map.as_ref().map_or(&[], |map| &map[..]) }
    pub(crate) fn len(&self) -> u64 { self.index.decoded_len }
    fn check_magic(&self) -> io::Result<()> {
        if self.len() == 0 && self.index.format == Format::Tss { return Ok(()); }
        let mut magic = [0; 8];
        self.read_at(0, &mut magic)?;
        if &magic != self.index.format.magic() { return Err(invalid("unsupported store format")); }
        Ok(())
    }
    fn build_index(&mut self) -> io::Result<()> {
        self.index.source_len = self.bytes().len() as u64;
        if self.bytes().starts_with(self.index.format.magic()) || self.bytes().is_empty() {
            self.index.decoded_len = self.index.source_len;
        } else if self.bytes().starts_with(&[31, 139]) {
            self.index.bgzf = true;
            let mut at = 0;
            while at < self.bytes().len() {
                let size = block_size(&self.bytes()[at..])?;
                let end = at.checked_add(size).filter(|&n| n <= self.bytes().len())
                    .ok_or_else(|| invalid("truncated BGZF block"))?;
                let decoded = u32::from_le_bytes(self.bytes()[end - 4..end].try_into().unwrap()) as u64;
                if decoded > 65536 { return Err(invalid("invalid BGZF uncompressed size")); }
                if decoded == 0 {
                    if end != self.bytes().len() { return Err(invalid("data after BGZF EOF block")); }
                } else {
                    self.index.blocks.push((at as u64, self.index.decoded_len));
                    self.index.decoded_len += decoded;
                }
                at = end;
            }
        } else { return Err(invalid("unsupported store format")); }
        self.check_magic()?;
        let mut at = self.index.format.first();
        while at < self.len() {
            let end = match self.index.format {
                Format::Tq => {
                    let mut h = [0; 12];
                    self.read_at(at, &mut h)?;
                    let name = u32::from_le_bytes(h[..4].try_into().unwrap()) as u64;
                    let count = u64::from_le_bytes(h[4..].try_into().unwrap());
                    count.checked_mul(2).and_then(|n| n.checked_add(12 + name)).and_then(|n| at.checked_add(n))
                }
                Format::Tss => {
                    let mut h = [0; 16];
                    self.read_at(at, &mut h)?;
                    if &h[..8] != self.index.format.magic() { return Err(invalid("invalid TSS record header")); }
                    let len = u64::from_le_bytes(h[8..].try_into().unwrap());
                    at.checked_add(16).and_then(|n| n.checked_add(len))
                }
            }.filter(|&end| end <= self.len()).ok_or_else(|| invalid("truncated store record"))?;
            self.index.records.push(at);
            at = end;
        }
        self.index.validate()
    }

    fn block(&self, number: usize) -> io::Result<Rc<Vec<u8>>> {
        let key = (self.namespace, number);
        if let Some(bytes) = CACHE.with(|cache| cache.borrow().blocks.get(&key).cloned()) {
            return Ok(bytes);
        }
        let (start, decoded) = self.index.blocks[number];
        let end = self
            .index
            .blocks
            .get(number + 1)
            .map_or(self.index.decoded_len, |p| p.1);
        let start = usize::try_from(start).map_err(|_| invalid("invalid store block offset"))?;
        let size = block_size(&self.bytes()[start..])?;
        let compressed = self
            .bytes()
            .get(start..start + size)
            .ok_or_else(|| invalid("truncated BGZF block"))?;
        let mut bytes = Vec::new();
        flate2::read::GzDecoder::new(compressed)
            .take(65537)
            .read_to_end(&mut bytes)?;
        if bytes.len() as u64 != end - decoded {
            return Err(invalid("BGZF block and store index disagree"));
        }
        let bytes = Rc::new(bytes);
        CACHE.with(|cache| cache.borrow_mut().insert(key, bytes.clone()));
        Ok(bytes)
    }

    pub(crate) fn read_at(&self, mut at: u64, mut out: &mut [u8]) -> io::Result<()> {
        if at
            .checked_add(out.len() as u64)
            .is_none_or(|end| end > self.index.decoded_len)
        {
            return Err(invalid("store byte range out of bounds"));
        }
        if !self.index.bgzf {
            out.copy_from_slice(&self.bytes()[at as usize..at as usize + out.len()]);
            return Ok(());
        }
        while !out.is_empty() {
            let block = self
                .index
                .blocks
                .partition_point(|p| p.1 <= at)
                .checked_sub(1)
                .ok_or_else(|| invalid("missing BGZF block"))?;
            let bytes = self.block(block)?;
            let offset = (at - self.index.blocks[block].1) as usize;
            let n = out.len().min(
                bytes
                    .len()
                    .checked_sub(offset)
                    .filter(|&n| n > 0)
                    .ok_or_else(|| invalid("store block offset out of bounds"))?,
            );
            out[..n].copy_from_slice(&bytes[offset..offset + n]);
            at += n as u64;
            out = &mut out[n..];
        }
        Ok(())
    }

    pub(crate) fn range(&self, at: u64, len: usize) -> io::Result<std::borrow::Cow<'_, [u8]>> {
        if at.checked_add(len as u64).is_none_or(|end| end > self.len()) {
            return Err(invalid("store byte range out of bounds"));
        }
        if !self.index.bgzf {
            return Ok(std::borrow::Cow::Borrowed(&self.bytes()[at as usize..at as usize + len]));
        }
        let mut bytes = vec![0; len];
        self.read_at(at, &mut bytes)?;
        Ok(std::borrow::Cow::Owned(bytes))
    }
}
#[derive(Clone, Copy)]
pub(crate) enum Bytes<'a> {
    Slice(&'a [u8]),
    Indexed { store: &'a Store, offset: u64, len: usize },
}
impl<'a> Bytes<'a> {
    pub(crate) fn len(self) -> usize { match self { Self::Slice(b) => b.len(), Self::Indexed {len, ..} => len } }
    pub(crate) fn slice(self, range: Range<usize>) -> io::Result<Self> {
        if range.start > range.end || range.end > self.len() { return Err(invalid("byte range out of bounds")); }
        Ok(match self {
            Self::Slice(b) => Self::Slice(&b[range]),
            Self::Indexed {store, offset, ..} => Self::Indexed {store, offset: offset + range.start as u64, len: range.len()},
        })
    }
    pub(crate) fn read(self) -> io::Result<std::borrow::Cow<'a, [u8]>> {
        match self {
            Self::Slice(b) => Ok(std::borrow::Cow::Borrowed(b)),
            Self::Indexed {store, offset, len} => store.range(offset, len),
        }
    }
    pub(crate) fn read_at(self, at: usize, out: &mut [u8]) -> io::Result<()> {
        let range = at..at.checked_add(out.len()).ok_or_else(|| invalid("byte range overflow"))?;
        match self.slice(range)? {
            Self::Slice(b) => {out.copy_from_slice(b); Ok(())},
            Self::Indexed {store, offset, ..} => store.read_at(offset, out),
        }
    }
}
