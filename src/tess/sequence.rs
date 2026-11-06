use super::{AMB, Tile};
use crate::tess::CHUNK;
use crate::unit::{MAX_PERIOD, pow4};
use crate::Unit;
use simd_minimizers::packed_seq::wide::u64x8;
pub fn parse(sequence: &[u8], k: usize) -> Vec<Tile> {
    assert!(k <= crate::tess::budgets::MAX_ERRORS, "k must be in 0..=8");
    assert!(sequence.len() <= u32::MAX as usize);
    super::edit::parse(&Sequence::new(sequence), k)
}

#[inline(always)]
fn splat(x: u64) -> u64x8 {
    u64x8::splat(x)
}
#[inline(always)]
fn even(x: u64x8) -> u64x8 {
    let e = !(x | x >> 1u32) & splat(0x5555_5555_5555_5555);
    let e = (e | e >> 1u32) & splat(0x3333_3333_3333_3333);
    let e = (e | e >> 2u32) & splat(0x0f0f_0f0f_0f0f_0f0f);
    let e = (e | e >> 4u32) & splat(0x00ff_00ff_00ff_00ff);
    let e = (e | e >> 8u32) & splat(0x0000_ffff_0000_ffff);
    (e | e >> 16u32) & splat(0x0000_0000_ffff_ffff)
}
fn pack64(seq: &[u8], i: usize) -> (u64, u64, u64) {
    let n = seq.len();
    let (lo7, high) = (splat(0x7f7f_7f7f_7f7f_7f7f), splat(0x8080_8080_8080_8080));
    let zero = |x: u64x8| !(((x & lo7) + lo7) | x) & high;
    let acgt = (*b"ACGT").map(|c| splat(0x0101_0101_0101_0101 * c as u64));
    let b: [u8; 64] = if 64 * i + 64 <= n {
        seq[64 * i..64 * i + 64].try_into().unwrap()
    } else {
        let mut b = [b'A'; 64];
        b[..n - 64 * i].copy_from_slice(&seq[64 * i..]);
        b
    };
    let w: [u64; 8] =
        std::array::from_fn(|j| u64::from_le_bytes(b[8 * j..8 * j + 8].try_into().unwrap()));
    let u = u64x8::from(w) & splat(0xdfdf_dfdf_dfdf_dfdf);
    let ok = zero(u ^ acgt[0]) | zero(u ^ acgt[1]) | zero(u ^ acgt[2]) | zero(u ^ acgt[3]);
    let q = (u >> 1u32) & ((ok >> 6u32) | (ok >> 7u32));
    let q = (q | q >> 6u32) & splat(0x000f_000f_000f_000f);
    let q = (q | q >> 12u32) & splat(0x0000_00ff_0000_00ff);
    let q = ((q | q >> 24u32) & splat(0xffff)).to_array();
    let x = (!ok & high) >> 7u32;
    let x = (x | x >> 7u32) & splat(0x0003_0003_0003_0003);
    let x = (x | x >> 14u32) & splat(0x0000_000f_0000_000f);
    let x = ((x | x >> 28u32) & splat(0xff)).to_array();
    (
        q[0] | q[1] << 16 | q[2] << 32 | q[3] << 48,
        q[4] | q[5] << 16 | q[6] << 32 | q[7] << 48,
        (0..8).fold(0, |a, j| a | x[j] << (8 * j)),
    )
}
fn compare(cur: u64x8, before: u64x8, amb: [u64; 5]) -> [[u64; MAX_PERIOD]; 4] {
    let mut out = [[0u64; MAX_PERIOD]; 4];
    for p in 1..=MAX_PERIOD {
        let (sl, sr) = (2 * p as u32, 64 - 2 * p as u32);
        let e = even(cur ^ ((cur << sl) | (before >> sr))).to_array();
        for k in 0..4 {
            let (ac, ap) = (amb[k + 1], (amb[k] >> (64 - p)) | (amb[k + 1] << p));
            out[k][p - 1] =
                (e[2 * k] | e[2 * k + 1] << 32) & if p == 1 { !(ac ^ ap) } else { !(ac | ap) };
        }
    }
    out
}
pub(crate) struct Sequence {
    pub(crate) n: usize,
    code: Vec<u64>,
    amb: Vec<u64>,
}
impl Sequence {
    #[inline]
    pub fn amb(&self, i: usize) -> bool {
        (self.amb[1 + i / 64] >> (i % 64)) & 1 == 1
    }
    #[inline]
    pub fn sym(&self, i: usize) -> u8 {
        if self.amb(i) {
            AMB
        } else {
            ((self.code[1 + i / 32] >> (2 * (i % 32))) & 3) as u8
        }
    }
    #[inline]
    pub(crate) fn id(&self, a: usize, l: usize) -> Unit {
        Unit::from_word(l, self.kmer(a, l))
    }
    #[inline]
    pub fn kmer(&self, a: usize, l: usize) -> usize {
        let w = 1 + a / 32;
        let x = (self.code[w] as u128 | (self.code[w + 1] as u128) << 64) >> (2 * (a % 32));
        x as usize & (pow4(l) - 1)
    }
    pub(crate) fn new(seq: &[u8]) -> Self {
        let n = seq.len();
        let words = n.div_ceil(64);
        let (mut code, mut amb) = (vec![0u64; 2 * words + 10], vec![0u64; words + 2]);
        code[1..1 + 2 * words]
            .chunks_mut(CHUNK / 32)
            .zip(amb[1..1 + words].chunks_mut(CHUNK / 64))
            .enumerate()
            .for_each(|(t, (ct, at))| {
                for (j, (c, a)) in ct.chunks_mut(2).zip(at.iter_mut()).enumerate() {
                    (c[0], c[1], *a) = pack64(seq, t * (CHUNK / 64) + j);
                }
            });
        Self { code, amb, n }
    }
    pub(crate) fn masks(&self) -> Vec<[u64; MAX_PERIOD]> {
        let at = |v: &[u64], i: usize| v.get(i).copied().unwrap_or(0);
        let mut masks: Vec<[u64; MAX_PERIOD]> = (0..self.n.div_ceil(256))
            .flat_map(|b| {
                let lanes = |o: usize| {
                    u64x8::from(std::array::from_fn::<u64, 8, _>(|j| {
                        at(&self.code, 8 * b + o + j)
                    }))
                };
                compare(
                    lanes(1),
                    lanes(0),
                    std::array::from_fn(|j| at(&self.amb, 4 * b + j)),
                )
            })
            .collect();
        masks.truncate(self.n.div_ceil(64));
        masks
    }
}
