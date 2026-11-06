use std::fs::File;
use std::io::{BufWriter, Write};
use sucds::Serializable;
use sucds::int_vectors::CompactVector;

pub const fn bits_for(n: usize) -> usize {
    let b = (usize::BITS - n.leading_zeros()) as usize;
    if b == 0 { 1 } else { b }
}
pub fn compact(width: usize, values: impl IntoIterator<Item = u64>) -> CompactVector {
    let it = values.into_iter();
    let mut v = CompactVector::with_capacity(it.size_hint().0, width).unwrap();
    v.extend(it).unwrap();
    v
}
pub fn write_u64s<W: Write>(
    w: &mut W,
    values: impl IntoIterator<Item = usize>,
) -> std::io::Result<()> {
    for x in values {
        w.write_all(&(x as u64).to_le_bytes())?;
    }
    Ok(())
}
pub fn u64_at(b: &[u8], i: usize) -> usize {
    u64::from_le_bytes(b[i..i + 8].try_into().unwrap()) as usize
}
pub fn read_u64(b: &mut &[u8]) -> usize {
    let x = u64_at(b, 0);
    *b = &b[8..];
    x
}
pub fn ser<T: Serializable, W: Write>(x: &T, w: &mut W) -> std::io::Result<()> {
    x.serialize_into(w)
        .map(|_| ())
        .map_err(std::io::Error::other)
}
pub fn de<T: Serializable>(b: &mut &[u8]) -> T {
    T::deserialize_from(b).unwrap()
}
pub fn create_with_capacity(path: &str, capacity: usize) -> BufWriter<File> {
    BufWriter::with_capacity(capacity, create_file(path))
}
fn create_file(path: &str) -> File {
    File::create(path).unwrap_or_else(|e| panic!("cannot write {path}: {e}"))
}
