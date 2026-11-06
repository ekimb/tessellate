use super::codec::{self, invalid, put_word, word};
pub use super::tile::{AMB, Tile};
use super::{EMPTY, block};
use crate::Unit;
use super::indexed::Bytes;
use std::borrow::{Borrow, Cow};
use std::io::{self, Read, Write};
use std::ops::Range;
use sucds::Serializable;

pub(crate) use super::block::Block as DecodedBlock;

const MAGIC: &[u8; 8] = b"TSSBLK2\0";
const MAGIC_V1: &[u8; 8] = b"TSSBLK1\0";
const HEADER: usize = 48;

fn check_range(range: &Range<usize>, len: usize) -> io::Result<()> {
    if range.start > range.end || range.end > len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "range out of bounds",
        ));
    }
    Ok(())
}
fn layout(header: &[u8]) -> io::Result<(usize, usize, usize, usize, usize)> {
    code_bits(header)?;
    let h = header
        .get(8..HEADER)
        .ok_or_else(|| invalid("truncated auxiliary header"))?;
    let (n, tiles, name, rate, payload) = (
        word(h, 0)?,
        word(h, 1)?,
        word(h, 2)?,
        word(h, 3)?,
        word(h, 4)?,
    );
    if rate == 0 || tiles > n || (tiles == 0) != (n == 0) {
        return Err(invalid("invalid auxiliary dimensions"));
    }
    Ok((n, tiles, name, rate, payload))
}
fn code_bits(header: &[u8]) -> io::Result<usize> {
    match header.get(..8) {
        Some(m) if m == MAGIC => Ok(block::CODE_BITS),
        Some(m) if m == MAGIC_V1 => Ok(2),
        _ => Err(invalid("unsupported auxiliary record")),
    }
}
fn record_length(header: &[u8]) -> io::Result<usize> {
    let (n, tiles, name, rate, payload) = layout(header)?;
    let count = tiles
        .div_ceil(rate)
        .checked_add(1)
        .ok_or_else(|| invalid("block count overflow"))?;
    [
        name,
        codec::boundary_bytes(n, count)?,
        codec::boundary_bytes(payload, count)?,
        payload,
    ]
    .into_iter()
    .try_fold(HEADER, |sum, n| {
        sum.checked_add(n)
            .ok_or_else(|| invalid("record length overflow"))
    })
}

pub struct Packed {
    bytes: Vec<u8>,
}
impl Packed {
    pub fn from_tiles_with_rate(
        name: &[u8],
        n: usize,
        tiles: &[Tile],
        rate: usize,
    ) -> io::Result<Self> {
        Self::encode(name, n, tiles, rate, MAGIC, block::CODE_BITS)
    }
    fn encode(
        name: &[u8],
        n: usize,
        tiles: &[Tile],
        rate: usize,
        magic: &[u8; 8],
        code_bits: usize,
    ) -> io::Result<Self> {
        if rate == 0 || (n == 0) != tiles.is_empty() {
            return Err(invalid("invalid tessellation or block rate"));
        }
        let mut starts = Vec::with_capacity(tiles.len() + 1);
        let mut start = 0usize;
        for (i, tile) in tiles.iter().enumerate() {
            starts.push(start);
            let span = tile.span as usize;
            let end = start + span;
            if span == 0
                || end > n
                || Unit::valid(tile.id.0).is_none_or(|u| u.len() != u.canon().len())
            {
                return Err(invalid("invalid tile"));
            }
            let mut previous = 0u32;
            for (offset, code) in tile.wire_ex()? {
                let terminal =
                    i + 1 == tiles.len() && end == n && offset == tile.span && code == EMPTY;
                if offset < previous || (offset >= tile.span && !terminal) || code as usize >= 1 << code_bits {
                    return Err(invalid("invalid tile exception"));
                }
                previous = offset;
            }
            start = end;
        }
        if start != n {
            return Err(invalid("tile spans do not cover the sequence"));
        }
        starts.push(n);
        let mut payload = Vec::new();
        let mut offsets = vec![0];
        let mut bases = vec![0];
        for (i, chunk) in tiles.chunks(rate).enumerate() {
            let end = starts[((i + 1) * rate).min(tiles.len())];
            payload.extend(block::encode_with(chunk, starts[i * rate], end, code_bits)?);
            offsets.push(payload.len());
            bases.push(end);
        }
        let mut bytes = magic.to_vec();
        for v in [n, tiles.len(), name.len(), rate, payload.len()] {
            put_word(&mut bytes, v);
        }
        bytes.extend(name);
        bytes.extend(codec::encode_boundaries(n, bases.len(), bases)?);
        bytes.extend(codec::encode_boundaries(
            payload.len(),
            offsets.len(),
            offsets,
        )?);
        bytes.extend(payload);
        Ok(Self { bytes })
    }
    pub fn view(bytes: &[u8]) -> io::Result<PackedView<'_>> {
        PackedView::from_bytes(bytes)
    }
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }
    pub fn as_view(&self) -> PackedView<'_> {
        Self::view(&self.bytes).expect("valid owned auxiliary record")
    }
    pub fn name(&self) -> &[u8] {
        let (_, _, name, _, _) = layout(&self.bytes).expect("valid owned auxiliary record");
        &self.bytes[HEADER..HEADER + name]
    }
    pub fn len(&self) -> usize {
        self.as_view().len()
    }
}
impl Serializable for Packed {
    fn serialize_into<W: Write>(&self, mut writer: W) -> sucds::Result<usize> {
        writer.write_all(&self.bytes)?;
        Ok(self.bytes.len())
    }
    fn deserialize_from<R: Read>(mut reader: R) -> sucds::Result<Self> {
        let mut header = [0; HEADER];
        reader.read_exact(&mut header)?;
        let length = record_length(&header)?;
        let mut bytes = Vec::new();
        bytes
            .try_reserve_exact(length)
            .map_err(|e| invalid(&e.to_string()))?;
        bytes.extend_from_slice(&header);
        bytes.resize(length, 0);
        reader.read_exact(&mut bytes[HEADER..])?;
        Self::view(&bytes)?;
        Ok(Self { bytes })
    }
    fn size_in_bytes(&self) -> usize {
        self.bytes.len()
    }
}

#[derive(Clone)]
struct Blocks<'a> {
    n: usize,
    tiles: usize,
    rate: usize,
    code_bits: usize,
    name: Cow<'a, [u8]>,
    bases: Bytes<'a>,
    offsets: Bytes<'a>,
    payload: Bytes<'a>,
}
impl Blocks<'_> {
    fn count(&self) -> usize {
        self.tiles.div_ceil(self.rate)
    }
    fn base(&self, index: usize) -> io::Result<usize> {
        codec::boundary(self.bases, self.n, self.count() + 1, index)
    }
    fn bytes(&self, index: usize) -> io::Result<(Cow<'_, [u8]>, usize)> {
        if index >= self.count() {
            return Err(invalid("tile block out of bounds"));
        }
        let count = self.count() + 1;
        let begin = codec::boundary(self.offsets, self.payload.len(), count, index)?;
        let end = codec::boundary(self.offsets, self.payload.len(), count, index + 1)?;
        let bytes = self.payload.slice(begin..end)?.read()?;
        let tiles = self.rate.min(self.tiles - index * self.rate);
        Ok((bytes, tiles))
    }
    fn decode(&self, index: usize) -> io::Result<block::Block> {
        let (bytes, tiles) = self.bytes(index)?;
        block::decode(
            &bytes,
            tiles,
            self.base(index)?,
            self.base(index + 1)?,
            index + 1 == self.count(),
            self.code_bits,
        )
    }
}
#[derive(Clone)]
pub struct PackedView<'a> {
    v: Blocks<'a>,
}
impl<'a> PackedView<'a> {
    pub fn from_bytes(bytes: &'a [u8]) -> io::Result<Self> {
        Self::from_source(Bytes::Slice(bytes))
    }
    pub(crate) fn from_source(bytes: Bytes<'a>) -> io::Result<Self> {
        let mut header = [0; HEADER];
        bytes.read_at(0, &mut header)?;
        let (n, tiles, name, rate, payload) = layout(&header)?;
        if record_length(&header)? != bytes.len() {
            return Err(invalid("auxiliary record length mismatch"));
        }
        let count = tiles.div_ceil(rate) + 1;
        let mut offset = HEADER;
        let name = bytes.slice(offset..offset + name)?.read()?;
        offset += name.len();
        let size = codec::boundary_bytes(n, count)?;
        let bases = bytes.slice(offset..offset + size)?;
        offset += size;
        let size = codec::boundary_bytes(payload, count)?;
        let offsets = bytes.slice(offset..offset + size)?;
        offset += size;
        let blocks = Blocks {
            n,
            tiles,
            rate,
            code_bits: code_bits(&header)?,
            name,
            bases,
            offsets,
            payload: bytes.slice(offset..bytes.len())?,
        };
        if blocks.base(0)? != 0
            || blocks.base(count - 1)? != n
            || codec::boundary(offsets, payload, count, 0)? != 0
            || codec::boundary(offsets, payload, count, count - 1)? != payload
        {
            return Err(invalid("invalid block directory"));
        }
        Ok(Self { v: blocks })
    }
    pub fn name(&self) -> &[u8] {
        &self.v.name
    }
    pub fn len(&self) -> usize {
        self.v.n
    }
    pub fn n_tiles(&self) -> usize {
        self.v.tiles
    }
    pub fn block_rate(&self) -> usize {
        self.v.rate
    }
    pub(crate) fn decode_block(&self, block: usize) -> io::Result<DecodedBlock> {
        self.v.decode(block)
    }
    pub fn start(&self, tile: usize) -> io::Result<usize> {
        check_range(&(tile..tile), self.n_tiles())?;
        let v = &self.v;
        if tile == v.tiles {
            Ok(v.n)
        } else if tile % v.rate == 0 {
            v.base(tile / v.rate)
        } else {
            Err(invalid("interior tile boundary requires canonical units"))
        }
    }
    pub(crate) fn block_boundaries_with_deltas(
        &self,
        block: usize,
        units: &[u16],
        deltas: &[usize],
    ) -> io::Result<Vec<usize>> {
        let v = &self.v;
        if deltas.len() != units.len() {
            return Err(invalid("unit count does not match tile block"));
        }
        let mut result = Vec::with_capacity(units.len() + 1);
        let mut start = v.base(block)?;
        result.push(start);
        for (&unit, &delta) in units.iter().zip(deltas) {
            let unit = Unit::canonical(unit).ok_or_else(|| invalid("invalid canonical unit"))?;
            start = start
                .checked_add(unit.len())
                .and_then(|x| x.checked_add(delta))
                .ok_or_else(|| invalid("tile boundary overflow"))?;
            result.push(start);
        }
        if start != v.base(block + 1)? {
            return Err(invalid(
                "canonical periods and span deltas disagree with block span",
            ));
        }
        Ok(result)
    }
    pub(crate) fn tiles_at_with<B: Borrow<DecodedBlock>>(
        &self,
        range: Range<usize>,
        units: &[u16],
        base: usize,
        mut decode: impl FnMut(usize) -> io::Result<B>,
    ) -> io::Result<Vec<Tile>> {
        check_range(&range, self.n_tiles())?;
        if range.len() != units.len() {
            return Err(invalid("unit count does not match tile range"));
        }
        if range.is_empty() {
            return Ok(Vec::new());
        }
        let v = &self.v;
        let mut result = Vec::with_capacity(units.len());
        let mut start = base;
        let mut at = range.start;
        while at < range.end {
            let block = at / v.rate;
            let first = block * v.rate;
            let stop = range.end.min(first.saturating_add(v.rate));
            let decoded = decode(block)?;
            let decoded = decoded.borrow();
            let mut exception = decoded.exceptions.partition_point(|&(position, _)| position < start);
            if at == first && start != v.base(block)? {
                return Err(invalid("incorrect block base start"));
            }
            for index in at..stop {
                let unit = Unit::canonical(units[index - range.start])
                    .ok_or_else(|| invalid("invalid canonical unit"))?;
                let phase = decoded.phases[index - first];
                let period = unit.len();
                if phase >= period {
                    return Err(invalid("phase exceeds unit period"));
                }
                let id = unit.rotate(phase);
                let span = period
                    .checked_add(decoded.deltas[index - first])
                    .ok_or_else(|| invalid("tile boundary overflow"))?;
                let end = start
                    .checked_add(span)
                    .ok_or_else(|| invalid("tile boundary overflow"))?;
                let first_exception = exception;
                let terminal = index + 1 == v.tiles;
                while exception < decoded.exceptions.len() {
                    let position = decoded.exceptions[exception].0;
                    if position > end || (position == end && !terminal) { break; }
                    exception += 1;
                }
                result.push(Tile::from_wire_iter(
                    id,
                    u32::try_from(span).map_err(|_| invalid("tile span overflow"))?,
                    decoded.exceptions[first_exception..exception].iter()
                        .map(|&(position, code)| ((position - start) as u32, code)),
                ));
                start = end;
            }
            if stop == first + decoded.deltas.len() && start != v.base(block + 1)? {
                return Err(invalid("tile spans disagree with block boundary"));
            }
            at = stop;
        }
        Ok(result)
    }
}
