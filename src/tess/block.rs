use super::codec::{invalid, put_varint, varint};
use super::{EMPTY, Tile};
use crate::Tessellable;
use std::io;

#[derive(Default)]
struct Bits {
    bytes: Vec<u8>,
    len: usize,
}
impl Bits {
    fn push(&mut self, value: usize, width: usize) {
        for bit in (0..width).rev() {
            if self.len % 8 == 0 {
                self.bytes.push(0);
            }
            self.bytes[self.len / 8] |= (((value >> bit) & 1) as u8) << (7 - self.len % 8);
            self.len += 1;
        }
    }
    fn gamma(&mut self, value: usize) {
        let width = usize::BITS as usize - value.leading_zeros() as usize;
        self.push(0, width - 1);
        self.push(value, width);
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    bit: usize,
}
impl Reader<'_> {
    fn take(&mut self, width: usize) -> io::Result<usize> {
        let mut value = 0;
        for _ in 0..width {
            let byte = *self
                .bytes
                .get(self.bit / 8)
                .ok_or_else(|| invalid("truncated tile block"))?;
            value = (value << 1) | ((byte >> (7 - self.bit % 8)) & 1) as usize;
            self.bit += 1;
        }
        Ok(value)
    }
    fn gamma(&mut self) -> io::Result<usize> {
        let mut zeros = 0;
        while self.take(1)? == 0 {
            zeros += 1;
            if zeros >= usize::BITS as usize {
                return Err(invalid("span delta overflow"));
            }
        }
        Ok((1usize << zeros) | self.take(zeros)?)
    }
}

pub(super) const CODE_BITS: usize = 3;
pub(super) fn encode_with(tiles: &[Tile], base: usize, end: usize, code_bits: usize) -> io::Result<Vec<u8>> {
    let mut lengths = Bits::default();
    let mut mask = vec![0u8; tiles.len().div_ceil(8)];
    let mut phases = Bits::default();
    let mut positions = Vec::new();
    let mut codes = Bits::default();
    let (mut previous, mut exceptions) = (base, 0);
    let mut start = base;
    for (i, tile) in tiles.iter().enumerate() {
        let span = tile.span as usize;
        let next = start + span;
        if next > end {
            return Err(invalid("unordered tiles"));
        }
        let delta = span
            .checked_sub(tile.len())
            .ok_or_else(|| invalid("tile span is shorter than its unit period"))?;
        lengths.gamma(delta + 1);
        let phase = tile.phase();
        if phase != 0 {
            mask[i / 8] |= 1 << (i % 8);
            phases.push(phase - 1, 3);
        }
        for (offset, code) in tile.wire_ex()? {
            let position = start + offset as usize;
            if position < previous
                || code as usize >= 1 << code_bits
                || position > next
                || (position == next && (next != end || code != EMPTY || i + 1 != tiles.len()))
            {
                return Err(invalid("invalid tile exception"));
            }
            put_varint(&mut positions, position - previous);
            codes.push(code as usize, code_bits);
            previous = position;
            exceptions += 1;
        }
        start = next;
    }
    let mut bytes = lengths.bytes;
    bytes.extend(mask);
    bytes.extend(phases.bytes);
    put_varint(&mut bytes, exceptions);
    bytes.extend(positions);
    bytes.extend(codes.bytes);
    Ok(bytes)
}

pub(crate) struct Block {
    pub deltas: Vec<usize>,
    pub phases: Vec<usize>,
    pub exceptions: Vec<(usize, u8)>,
}

impl Block {
    pub(crate) fn size_in_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + (self.deltas.capacity() + self.phases.capacity()) * std::mem::size_of::<usize>()
            + self.exceptions.capacity() * std::mem::size_of::<(usize, u8)>()
    }
}

pub(super) fn deltas(bytes: &[u8], count: usize) -> io::Result<(Vec<usize>, usize)> {
    let mut reader = Reader { bytes, bit: 0 };
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        values.push(reader.gamma()? - 1);
    }
    Ok((values, reader.bit.div_ceil(8)))
}

pub(super) fn decode(
    bytes: &[u8],
    count: usize,
    base: usize,
    end: usize,
    terminal: bool,
    code_bits: usize,
) -> io::Result<Block> {
    let (deltas, mut cursor) = deltas(bytes, count)?;
    let mask_len = count.div_ceil(8);
    let mask = bytes
        .get(cursor..cursor + mask_len)
        .ok_or_else(|| invalid("truncated phase bitmap"))?;
    if count % 8 != 0 && mask.last().is_some_and(|b| b >> (count % 8) != 0) {
        return Err(invalid("phase bitmap exceeds tile block"));
    }
    cursor += mask_len;
    let mut reader = Reader {
        bytes: &bytes[cursor..],
        bit: 0,
    };
    let mut phases = Vec::with_capacity(count);
    for i in 0..count {
        let phase = if mask[i / 8] >> (i % 8) & 1 != 0 {
            reader.take(3)? + 1
        } else {
            0
        };
        if phase > 6 {
            return Err(invalid("phase exceeds maximum period"));
        }
        phases.push(phase);
    }
    cursor += reader.bit.div_ceil(8);
    let n = varint(bytes, &mut cursor)?;
    if n > bytes.len().saturating_sub(cursor) {
        return Err(invalid("exception count exceeds tile block"));
    }
    let mut exceptions = Vec::with_capacity(n);
    let mut previous = base;
    for _ in 0..n {
        previous = previous
            .checked_add(varint(bytes, &mut cursor)?)
            .ok_or_else(|| invalid("exception position overflow"))?;
        if previous > end || (previous == end && !terminal) {
            return Err(invalid("exception outside tile block"));
        }
        exceptions.push((previous, 0));
    }
    let mut reader = Reader {
        bytes: &bytes[cursor..],
        bit: 0,
    };
    for (position, code) in &mut exceptions {
        *code = reader.take(code_bits)? as u8;
        if *position == end && *code != EMPTY {
            return Err(invalid("invalid terminal exception"));
        }
    }
    cursor += reader.bit.div_ceil(8);
    if cursor != bytes.len() {
        return Err(invalid("trailing tile-block bytes"));
    }
    Ok(Block {
        deltas,
        phases,
        exceptions,
    })
}
