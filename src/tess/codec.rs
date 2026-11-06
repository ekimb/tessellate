use std::io::{self, ErrorKind};
use super::indexed::Bytes;

pub(super) const SAMPLE: usize = 64;

pub(super) fn invalid(message: &str) -> io::Error {
    io::Error::new(ErrorKind::InvalidData, message)
}

pub(super) fn word(bytes: &[u8], index: usize) -> io::Result<usize> {
    let offset = index
        .checked_mul(8)
        .ok_or_else(|| invalid("word offset overflow"))?;
    let bytes = bytes
        .get(offset..)
        .and_then(|b| b.get(..8))
        .ok_or_else(|| invalid("truncated word"))?;
    usize::try_from(u64::from_le_bytes(bytes.try_into().unwrap()))
        .map_err(|_| invalid("word exceeds address space"))
}

pub(super) fn put_word(bytes: &mut Vec<u8>, value: usize) {
    bytes.extend_from_slice(&(value as u64).to_le_bytes());
}
struct BoundaryLayout {
    low_bits: usize,
    sample_bytes: usize,
    low_bytes: usize,
    high_bits: usize,
    bytes: usize,
}

impl BoundaryLayout {
    fn new(last: usize, count: usize) -> io::Result<Self> {
        if count == 0 {
            return Err(invalid("empty boundary index"));
        }
        let universe = last
            .checked_add(1)
            .ok_or_else(|| invalid("boundary universe overflow"))?;
        let ratio = (universe / count).max(1);
        let low_bits = (usize::BITS - 1 - ratio.leading_zeros()) as usize;
        let sample_bytes = count
            .div_ceil(SAMPLE)
            .checked_mul(8)
            .ok_or_else(|| invalid("boundary samples overflow"))?;
        let low_bytes = count
            .checked_mul(low_bits)
            .ok_or_else(|| invalid("boundary low bits overflow"))?
            .div_ceil(8);
        let high_bits = (last >> low_bits)
            .checked_add(count)
            .ok_or_else(|| invalid("boundary high bits overflow"))?;
        let bytes = sample_bytes
            .checked_add(low_bytes)
            .and_then(|n| n.checked_add(high_bits.div_ceil(8)))
            .ok_or_else(|| invalid("boundary length overflow"))?;
        Ok(Self {
            low_bits,
            sample_bytes,
            low_bytes,
            high_bits,
            bytes,
        })
    }
}

pub(super) fn boundary_bytes(last: usize, count: usize) -> io::Result<usize> {
    Ok(BoundaryLayout::new(last, count)?.bytes)
}

pub(super) fn encode_boundaries(
    last: usize,
    count: usize,
    positions: impl IntoIterator<Item = usize>,
) -> io::Result<Vec<u8>> {
    let layout = BoundaryLayout::new(last, count)?;
    let mut bytes = vec![0; layout.bytes];
    let mut previous = 0;
    let mut seen = 0;
    for (i, position) in positions.into_iter().enumerate() {
        if i >= count || position < previous || position > last {
            return Err(invalid("invalid boundary positions"));
        }
        let high = (position >> layout.low_bits) + i;
        if i % SAMPLE == 0 {
            let offset = i / SAMPLE * 8;
            bytes[offset..offset + 8].copy_from_slice(&(high as u64).to_le_bytes());
        }
        for bit in 0..layout.low_bits {
            let at = i * layout.low_bits + bit;
            bytes[layout.sample_bytes + at / 8] |= (((position >> bit) & 1) as u8) << (at % 8);
        }
        bytes[layout.sample_bytes + layout.low_bytes + high / 8] |= 1 << (high % 8);
        previous = position;
        seen += 1;
    }
    if seen != count {
        return Err(invalid("boundary count mismatch"));
    }
    Ok(bytes)
}

pub(super) fn boundary(bytes: Bytes<'_>, last: usize, count: usize, index: usize) -> io::Result<usize> {
    if index >= count {
        return Err(invalid("boundary index out of bounds"));
    }
    let layout = BoundaryLayout::new(last, count)?;
    if bytes.len() != layout.bytes {
        return Err(invalid("boundary length mismatch"));
    }
    let high = bytes.slice(layout.sample_bytes + layout.low_bytes..bytes.len())?;
    let mut sample = [0; 8];
    bytes.read_at(index / SAMPLE * 8, &mut sample)?;
    let sample = usize::try_from(u64::from_le_bytes(sample)).map_err(|_| invalid("boundary sample overflow"))?;
    if sample >= layout.high_bits { return Err(invalid("invalid boundary sample")); }
    let mut buffer = [0; 64];
    let mut first = sample / 8;
    let mut buffered = buffer.len().min(high.len() - first);
    high.read_at(first, &mut buffer[..buffered])?;
    if buffer[0] >> (sample % 8) & 1 == 0 {
        return Err(invalid("invalid boundary sample"));
    }
    let mut byte = sample / 8;
    let mut bits = buffer[0] & (u8::MAX << (sample % 8));
    let mut remaining = index % SAMPLE + 1;
    loop {
        let ones = bits.count_ones() as usize;
        if ones >= remaining {
            break;
        }
        remaining -= ones;
        byte += 1;
        if byte == first + buffered {
            if byte >= high.len() { return Err(invalid("truncated boundary index")); }
            first = byte;
            buffered = buffer.len().min(high.len() - first);
            high.read_at(first, &mut buffer[..buffered])?;
        }
        bits = buffer[byte - first];
    }
    for _ in 1..remaining {
        bits &= bits - 1;
    }
    let selected = byte * 8 + bits.trailing_zeros() as usize;
    if selected >= layout.high_bits || selected < index {
        return Err(invalid("invalid selected boundary"));
    }
    let mut low = 0;
    if layout.low_bits > 0 {
        let start = index * layout.low_bits;
        let first = start / 8;
        let len = (start + layout.low_bits).div_ceil(8) - first;
        let mut low_bytes = [0; 9];
        bytes.read_at(layout.sample_bytes + first, &mut low_bytes[..len])?;
        for bit in 0..layout.low_bits {
            let at = start % 8 + bit;
            low |= (((low_bytes[at / 8] >> (at % 8)) & 1) as usize) << bit;
        }
    }
    let value = (selected - index)
        .checked_shl(layout.low_bits as u32)
        .ok_or_else(|| invalid("boundary value overflow"))?
        | low;
    if value > last {
        return Err(invalid("boundary exceeds sequence"));
    }
    Ok(value)
}

pub(super) fn put_varint(bytes: &mut Vec<u8>, mut value: usize) {
    while value >= 128 {
        bytes.push(value as u8 | 128);
        value >>= 7;
    }
    bytes.push(value as u8);
}

pub(super) fn varint(bytes: &[u8], cursor: &mut usize) -> io::Result<usize> {
    let mut value = 0usize;
    for shift in (0..usize::BITS).step_by(7) {
        let b = *bytes
            .get(*cursor)
            .ok_or_else(|| invalid("truncated position"))?;
        *cursor += 1;
        if (b & 127) as usize > usize::MAX >> shift {
            return Err(invalid("position overflow"));
        }
        value |= ((b & 127) as usize) << shift;
        if b & 128 == 0 {
            return Ok(value);
        }
    }
    Err(invalid("position overflow"))
}
