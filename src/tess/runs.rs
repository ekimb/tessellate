use super::packed::AMB;
use super::sequence::Sequence;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Piece {
    pub start: usize,
    pub id: Unit,
}
impl Tessellable for Piece {
    fn unit(&self) -> Unit {
        self.id
    }
}
use crate::unit::{MAX_PERIOD, rank};
use crate::{Tessellable, Unit};
use simd_minimizers::packed_seq::wide::u64x8;

const _: () = assert!(MAX_PERIOD >= 4, "gap landmarks require units of length four");

#[inline(always)]
fn tss_walk(
    m: u64,
    starts: u64,
    s: usize,
    open: &mut Option<usize>,
    emit: &mut impl FnMut(usize, usize),
) {
    let mut off = 0u32;
    if let Some(k) = *open {
        let t = m.trailing_ones();
        if t == 64 {
            return;
        }
        emit(k, s + t as usize);
        *open = None;
        off = t + 1;
    }
    while off < 64 {
        let x = starts >> off;
        if x == 0 {
            break;
        }
        let j = off + x.trailing_zeros();
        let t = (m >> j).trailing_ones();
        if j + t >= 64 {
            *open = Some(s + j as usize);
            return;
        }
        emit(s + j as usize, s + (j + t) as usize);
        off = j + t + 1;
    }
    if off < 64 {
        let lo = m.leading_ones();
        if lo > 0 && 64 - lo >= off {
            *open = Some(s + (64 - lo) as usize);
        }
    }
}
#[inline]
fn run_starts(words: &[u64], p: usize) -> [u64; 8] {
    let masks = u64x8::from(std::array::from_fn(|i| words.get(i).copied().unwrap_or(0)));
    let mut starts = masks;
    for shift in 1..p {
        starts = starts & (masks >> shift as u32);
    }
    starts.to_array()
}

impl Sequence {
    fn equal(masks: &[[u64; MAX_PERIOD]], lo: usize, hi: usize, p: usize) -> Vec<u64> {
        let mut m: Vec<u64> = (lo / 64..hi.div_ceil(64))
            .map(|i| masks[i][p - 1])
            .collect();
        m[0] &= !((1u64 << p) - 1);
        if (hi - lo) % 64 > 0 {
            *m.last_mut().unwrap() &= (1u64 << ((hi - lo) % 64)) - 1;
        }
        m
    }
    fn maximal(
        &self,
        masks: &[[u64; MAX_PERIOD]],
        lo: usize,
        hi: usize,
        p: usize,
        out: &mut Vec<(Unit, usize, usize)>,
    ) {
        if hi - lo <= p {
            return;
        }
        let mut emit = |k: usize, b: usize| {
            if b - k < p {
                return;
            }
            let a = k - p;
            let id = if p == 1 && self.amb(a) {
                Unit::N
            } else {
                let r = self.id(a, p).primitive();
                if r.len() < p {
                    return;
                }
                r
            };
            out.push((id, a, b));
        };
        let mut open: Option<usize> = None;
        for (chunk, words) in Self::equal(masks, lo, hi, p).chunks(8).enumerate() {
            let starts = run_starts(words, p);
            for (lane, &m) in words.iter().enumerate() {
                if m == 0 && open.is_none() {
                    continue;
                }
                tss_walk(
                    m,
                    starts[lane],
                    lo + 64 * (8 * chunk + lane),
                    &mut open,
                    &mut emit,
                );
            }
        }
        if let Some(k) = open {
            emit(k, hi);
        }
    }
    fn cuts(&self, lo: usize, hi: usize, cut: &mut Vec<usize>, raw: &mut Vec<Piece>) {
        let k = hi - lo;
        if k == 0 {
            return;
        }
        let r = |j: usize| -> u8 {
            let s = self.sym(lo + j);
            if s == AMB { 255 } else { rank(s) }
        };
        let push = |cut: &mut Vec<usize>, v: usize| {
            if cut.last() != Some(&v) {
                cut.push(v)
            }
        };
        cut.clear();
        cut.push(0);
        for j in 0..k {
            let y = r(j);
            if y == 255 {
                push(cut, j);
                if j + 1 < k {
                    push(cut, j + 1);
                }
                continue;
            }
            if j == 0 {
                continue;
            }
            let x = r(j - 1);
            if x == y {
                push(cut, j);
            }
            if j == k - 1 {
                continue;
            }
            let z = r(j + 1);
            if x == 255 || z == 255 {
                continue;
            }
            let c = if x < y && y > z {
                j
            } else if x > y && y < z {
                j + 1
            } else {
                continue;
            };
            push(cut, c);
        }
        push(cut, k);
        for w in cut.windows(2) {
            let (a, l) = (lo + w[0], w[1] - w[0]);
            debug_assert!((1..=4).contains(&l));
            raw.push(Piece {
                start: a,
                id: if l == 1 && self.amb(a) {
                    Unit::N
                } else {
                    self.id(a, l).primitive()
                },
            });
        }
    }
    fn candidates(&self, lo: usize, hi: usize) -> Vec<(Unit, usize, usize)> {
        let masks = self.masks();
        let mut streams = std::array::from_fn::<_, MAX_PERIOD, _>(|i| {
            let mut runs = Vec::new();
            self.maximal(&masks, lo, hi, i + 1, &mut runs);
            runs.into_iter().peekable()
        });
        let (mut rs, mut max_end, mut prev) = (Vec::new(), 0, None);
        loop {
            let next = streams.iter_mut().enumerate().filter_map(|(i, stream)| {
                stream.peek().map(|&(id, s, e)| ((s, std::cmp::Reverse(e), id.len(), id), i))
            }).min().map(|(_, i)| i);
            let Some(i) = next else { break };
            let (id, s, e) = streams[i].next().unwrap();
            if prev == Some((s, e)) || e <= max_end {
                continue;
            }
            (max_end, prev) = (e, Some((s, e)));
            rs.push((id, s, e));
        }
        rs
    }
    fn select(&self, lo: usize, hi: usize) -> Vec<(Unit, usize, usize)> {
        let rs = self.candidates(lo, hi);
        let (mut ks, mut ke): (Vec<usize>, Vec<usize>) = (
            rs.iter().map(|r| r.1).collect(),
            rs.iter().map(|r| r.2).collect(),
        );
        for i in 0..rs.len() {
            for j in i + 1..rs.len() {
                let (a, b) = (rs[j].1, rs[i].2);
                if a >= b {
                    break;
                }
                let (ci, cj) = (rs[i].0.class(), rs[j].0.class());
                if ci < cj {
                    ks[j] = ks[j].max(b);
                } else if cj < ci {
                    ke[i] = ke[i].min(a);
                } else {
                    let h = (b - a) / 2;
                    ke[i] = ke[i].min(a + h);
                    ks[j] = ks[j].max(b - h);
                }
            }
        }
        rs.iter()
            .enumerate()
            .filter(|(i, r)| ke[*i] >= ks[*i] + r.0.len())
            .map(|(i, r)| (r.0, ks[i], ke[i]))
            .collect()
    }
    pub fn raw_pieces(&self, lo: usize, hi: usize) -> Vec<Piece> {
        let mut raw: Vec<Piece> = Vec::new();
        let mut cut = Vec::new();
        let mut cur = lo;
        for (id, s, e) in self.select(lo, hi) {
            debug_assert!(s >= cur, "trimmed runs overlap");
            if s > cur {
                self.cuts(cur, s, &mut cut, &mut raw);
            }
            raw.push(Piece {
                start: s,
                id: if id.is_n() {
                    id
                } else {
                    self.id(s, id.len()).primitive()
                },
            });
            cur = e;
        }
        if cur < hi {
            self.cuts(cur, hi, &mut cut, &mut raw);
        }
        raw
    }
}