use super::alphabet;
use quadrank::{QuadRank16, QuadRanker};
use rayon::prelude::*;

const LEAF: usize = 1 << 31;
const NONE: usize = usize::MAX;
#[inline]
fn digit(code: (u64, u8), depth: usize) -> usize {
    ((code.0 >> (2 * (code.1 as usize - 1 - depth))) & 3) as usize
}
#[derive(Default)]
pub(crate) struct Qm {
    pub(crate) levels: Vec<QuadRank16>,
    pub(crate) packed: Vec<Vec<u64>>,
    pub(crate) offs: Vec<[usize; 4]>,
    pub(crate) codes: Vec<(u64, u8)>,
}
fn pack(v: &[u16], codes: &[(u64, u8)], depth: usize) -> Vec<u64> {
    v.par_chunks(32)
        .map(|c| {
            c.iter().enumerate().fold(0u64, |w, (j, &x)| {
                w | (digit(codes[x as usize], depth) as u64) << (2 * j)
            })
        })
        .collect()
}

fn encode(
    v: Vec<u16>,
    node: usize,
    depth: usize,
    codes: &[(u64, u8)],
    children: &[[usize; 4]],
) -> Vec<(usize, Vec<u64>)> {
    let packed = pack(&v, codes, depth);
    let parts: Vec<Vec<u16>> = (0..4)
        .into_par_iter()
        .map(|d| {
            if children[node][d] == NONE || children[node][d] & LEAF != 0 {
                Vec::new()
            } else {
                v.par_iter()
                    .copied()
                    .filter(|&u| digit(codes[u as usize], depth) == d)
                    .collect()
            }
        })
        .collect();
    drop(v);
    let subtrees: Vec<_> = parts
        .into_par_iter()
        .enumerate()
        .filter(|(_, part)| !part.is_empty())
        .map(|(d, part)| encode(part, children[node][d], depth + 1, codes, children))
        .collect();
    let mut out = vec![(node, packed)];
    out.extend(subtrees.into_iter().flatten());
    out
}

impl Qm {
    pub(crate) fn children(codes: &[(u64, u8)]) -> Vec<[usize; 4]> {
        let mut children = vec![[NONE; 4]];
        for (u, &code) in codes.iter().enumerate() {
            let mut node = 0;
            for depth in 0..code.1 as usize {
                let d = digit(code, depth);
                if depth + 1 == code.1 as usize {
                    assert_eq!(children[node][d], NONE, "Huffman prefix collision");
                    children[node][d] = LEAF | u;
                } else {
                    if children[node][d] == NONE {
                        children[node][d] = children.len();
                        children.push([NONE; 4]);
                    }
                    assert_eq!(children[node][d] & LEAF, 0, "Huffman prefix collision");
                    node = children[node][d];
                }
            }
        }
        children
    }
    pub(crate) fn build(v: Vec<u16>, codes: Vec<(u64, u8)>) -> Qm {
        let children = Self::children(&codes);
        let mut packed = vec![Vec::new(); children.len()];
        for (node, words) in encode(v, 0, 0, &codes, &children) {
            packed[node] = words;
        }
        Qm::with(packed, children, codes)
    }
    pub(crate) fn with(packed: Vec<Vec<u64>>, offs: Vec<[usize; 4]>, codes: Vec<(u64, u8)>) -> Qm {
        let levels = packed
            .par_iter()
            .map(|w| QuadRank16::new_packed(w))
            .collect();
        Qm {
            levels,
            packed,
            offs,
            codes,
        }
    }
    #[inline(never)]
    pub(crate) fn range(&self, mut lo: usize, mut hi: usize, v: u16) -> (usize, usize, usize) {
        let (mut node, mut less) = (0, 0);
        let code = self.codes[v as usize];
        for depth in 0..code.1 as usize {
            let d = digit(code, depth);
            let (ra, re) = unsafe {
                (
                    self.levels[node].rank4_unchecked(lo),
                    self.levels[node].rank4_unchecked(hi),
                )
            };
            for x in 0..d {
                less += (re[x] - ra[x]) as usize;
            }
            (lo, hi) = (ra[d] as usize, re[d] as usize);
            let child = self.offs[node][d];
            if child & LEAF != 0 {
                let mate = alphabet::complement(v);
                if mate != v {
                    let other = self.codes[mate as usize];
                    debug_assert_eq!((code.0 >> 2, code.1), (other.0 >> 2, other.1));
                    let md = (other.0 & 3) as usize;
                    let count = (re[md] - ra[md]) as usize;
                    if d < md { less += count } else { less -= count }
                }
                return (lo, hi - lo, less);
            }
            node = child;
        }
        unreachable!("canonical unit has no Huffman leaf")
    }
    #[inline]
    pub(crate) fn access_rank(&self, mut i: usize) -> (u16, usize) {
        let mut node = 0;
        loop {
            let (r0, r1) = unsafe {
                (
                    self.levels[node].rank4_unchecked(i),
                    self.levels[node].rank4_unchecked(i + 1),
                )
            };
            let d = (0..4).find(|&x| r1[x] != r0[x]).unwrap();
            i = r0[d] as usize;
            let child = self.offs[node][d];
            if child & LEAF != 0 {
                return ((child & !LEAF) as u16, i);
            }
            node = child;
        }
    }
}
