use super::chain::Chain;
use crate::Unit;

const CLIP: i32 = 5;
const ZDROP: i32 = 50;
const BAND: usize = 64;
const NEG: i32 = i32::MIN / 4;

#[derive(Default)]
pub struct Extension {
    pub pairs: Vec<(usize, usize)>,
    pub divergences: Vec<usize>,
}

enum Op {
    Diag(usize, usize, bool),
    Ins(usize),
    Del(usize),
}

fn align(read: impl Fn(usize) -> u16, n: usize, reference: impl Fn(usize) -> u16, m: usize) -> Vec<Op> {
    if n == 0 {
        return Vec::new();
    }
    let width = 2 * BAND + 1;
    let (mut prev, mut cur) = (vec![NEG; width], vec![NEG; width]);
    let mut dir: Vec<u8> = vec![0; width];
    for j in 0..=BAND.min(m) {
        prev[j + BAND] = -(j as i32);
        dir[j + BAND] = 3;
    }
    let (mut best, mut bi, mut bj) = (0, 0, 0);
    let mut rows = n;
    for i in 1..=n {
        cur.fill(NEG);
        dir.extend(std::iter::repeat_n(0, width));
        let row = i * width;
        let mut top = NEG;
        for j in i.saturating_sub(BAND)..=(i + BAND).min(m) {
            let off = j + BAND - i;
            let (mut v, mut d) = (NEG, 0u8);
            if j > 0 && prev[off] > NEG {
                let s = prev[off] + if read(i - 1) == reference(j - 1) { 1 } else { -1 };
                (v, d) = (s, 1);
            }
            if off + 1 < width && prev[off + 1] > NEG && prev[off + 1] - 1 > v {
                (v, d) = (prev[off + 1] - 1, 2);
            }
            if off > 0 && cur[off - 1] > NEG && cur[off - 1] - 1 > v {
                (v, d) = (cur[off - 1] - 1, 3);
            }
            cur[off] = v;
            dir[row + off] = d;
            if v > top {
                top = v;
            }
            if v > best {
                (best, bi, bj) = (v, i, j);
            }
        }
        std::mem::swap(&mut prev, &mut cur);
        if top < best - ZDROP {
            rows = i;
            break;
        }
    }
    let mut end = (bi, bj);
    if rows == n {
        let (mut gscore, mut gj) = (NEG, 0);
        for j in n.saturating_sub(BAND)..=(n + BAND).min(m) {
            let v = prev[j + BAND - n];
            if v > gscore {
                (gscore, gj) = (v, j);
            }
        }
        if gscore > NEG && gscore >= best - CLIP {
            end = (n, gj);
        }
    }
    let (mut i, mut j) = end;
    let mut ops = Vec::new();
    while i > 0 || j > 0 {
        match dir[i * width + j + BAND - i] {
            1 => {
                ops.push(Op::Diag(i - 1, j - 1, read(i - 1) == reference(j - 1)));
                (i, j) = (i - 1, j - 1);
            }
            2 => {
                ops.push(Op::Ins(i - 1));
                i -= 1;
            }
            3 => {
                ops.push(Op::Del(i));
                j -= 1;
            }
            _ => break,
        }
    }
    ops
}

pub fn extend(chain: &Chain, units: &[u16], reference: &[u16]) -> Extension {
    let nq = units.len();
    let oriented = |x: usize| if chain.minus { Unit(units[nq - 1 - x]).rc().0 } else { units[x] };
    let fwd = |x: usize| if chain.minus { nq - 1 - x } else { x };
    let mut ext = Extension::default();
    let first = &chain.anchors[0];
    let last = chain.anchors.iter().max_by_key(|a| a.q + a.len).unwrap();
    let (q0, r0) = (first.q as usize, first.r as usize);
    let left = align(|i| oriented(q0 - 1 - i), q0, |j| reference[r0 - 1 - j], r0.min(q0 + BAND));
    for op in left {
        match op {
            Op::Diag(i, j, same) => {
                ext.pairs.push((fwd(q0 - 1 - i), r0 - 1 - j));
                if !same {
                    ext.divergences.push(fwd(q0 - 1 - i));
                }
            }
            Op::Ins(i) => ext.divergences.push(fwd(q0 - 1 - i)),
            Op::Del(i) => ext.divergences.push(fwd(q0 - i.min(q0))),
        }
    }
    let (qe, re) = ((last.q + last.len) as usize, (last.r + last.len) as usize);
    let n = nq - qe;
    let right = align(|i| oriented(qe + i), n, |j| reference[re + j], (reference.len() - re).min(n + BAND));
    for op in right {
        match op {
            Op::Diag(i, j, same) => {
                ext.pairs.push((fwd(qe + i), re + j));
                if !same {
                    ext.divergences.push(fwd(qe + i));
                }
            }
            Op::Ins(i) => ext.divergences.push(fwd(qe + i)),
            Op::Del(i) => ext.divergences.push(fwd((qe + i).min(nq - 1))),
        }
    }
    ext
}
