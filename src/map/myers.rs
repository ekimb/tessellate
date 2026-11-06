use super::Record;
use super::chain::Chain;
use super::extend::Extension;
use std::collections::HashMap;

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub window: usize,
    pub pad: usize,
}

fn code(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

fn myers(read: &[u8], reference: &[u8]) -> u32 {
    let m = read.len();
    let mut peq = [0u64; 5];
    for (i, &b) in read.iter().enumerate() {
        peq[code(b)] |= 1 << i;
    }
    let mask = if m == 64 { !0 } else { (1u64 << m) - 1 };
    let high = 1u64 << (m - 1);
    let (mut pv, mut mv, mut score) = (mask, 0u64, m as u32);
    let mut best = score;
    for &b in reference {
        let eq = peq[code(b)];
        let xv = eq | mv;
        let xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq;
        let mut ph = mv | !(xh | pv);
        let mut mh = pv & xh;
        if ph & high != 0 {
            score += 1;
        } else if mh & high != 0 {
            score -= 1;
        }
        ph <<= 1;
        mh <<= 1;
        pv = (mh | !(xv | ph)) & mask;
        mv = ph & xv & mask;
        best = best.min(score);
    }
    best
}

fn dp(read: &[u8], reference: &[u8]) -> u32 {
    let mut prev = vec![0u32; reference.len() + 1];
    let mut cur = Vec::with_capacity(reference.len() + 1);
    for (i, &a) in read.iter().enumerate() {
        cur.clear();
        cur.push(i as u32 + 1);
        for (j, &b) in reference.iter().enumerate() {
            cur.push((prev[j] + (a != b) as u32).min(prev[j + 1] + 1).min(cur[j] + 1));
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev.into_iter().min().unwrap_or(0)
}

pub fn distance(read: &[u8], reference: &[u8]) -> u32 {
    if read.is_empty() {
        0
    } else if read.len() <= 64 {
        myers(read, reference)
    } else {
        dp(read, reference)
    }
}

fn revcomp(s: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(s.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }));
}

pub struct Placement {
    pub minus: bool,
    pub target: usize,
    lo: usize,
    hi: usize,
    map: Vec<u32>,
    divergences: Vec<usize>,
}

impl Placement {
    fn from_pairs(minus: bool, target: usize, n: usize, pairs: impl Iterator<Item = (usize, usize)>, mut divergences: Vec<usize>) -> Placement {
        let mut map = vec![u32::MAX; n];
        for (q, r) in pairs {
            if map[q] == u32::MAX {
                map[q] = r as u32;
            }
        }
        divergences.sort_unstable();
        divergences.dedup();
        let (lo, hi) = match (map.iter().position(|&r| r != u32::MAX), map.iter().rposition(|&r| r != u32::MAX)) {
            (Some(lo), Some(hi)) => (lo, hi),
            _ => (1, 0),
        };
        if lo <= hi {
            map.truncate(hi + 1);
            map.drain(..lo);
        } else {
            map.clear();
        }
        Placement { minus, target, lo, hi, map, divergences }
    }

    pub fn from_chain(chain: &Chain, read: &Record, extension: Extension) -> Placement {
        let n = read.n_tiles();
        let fwd = |q: usize| if chain.minus { n - 1 - q } else { q };
        let mut pairs: Vec<(usize, usize)> = chain.anchors.iter().flat_map(|a| (0..a.len as usize).map(move |k| (fwd(a.q as usize + k), a.r as usize + k))).collect();
        let mut divergences = Vec::new();
        for w in chain.anchors.windows(2) {
            let (j, i) = (&w[0], &w[1]);
            let (qe, q) = ((j.q + j.len) as usize, i.q as usize);
            let (re, r) = ((j.r + j.len) as usize, i.r as usize);
            if q > qe {
                let span = r.saturating_sub(re);
                pairs.extend((qe..q).map(|x| (fwd(x), re + (x - qe) * span / (q - qe))));
            }
            let shifted = i.q as i64 - j.q as i64 != i.r as i64 - j.r as i64;
            if q > qe || shifted {
                let (a, b) = (qe.min(q), qe.max(q).max(qe.min(q) + 1));
                divergences.extend((a..b).filter(|&x| x < n).map(fwd));
            }
        }
        pairs.extend(extension.pairs);
        divergences.extend(extension.divergences);
        Placement::from_pairs(chain.minus, chain.rec, n, pairs.into_iter(), divergences)
    }

    fn span(&self, wl: usize, wh: usize) -> Option<(usize, usize)> {
        let (wl, wh) = (wl.max(self.lo), wh.min(self.hi));
        if wl > wh {
            return None;
        }
        let r1 = (wl..=wh).map(|x| self.map[x - self.lo]).find(|&r| r != u32::MAX)? as usize;
        let r2 = (wl..=wh).rev().map(|x| self.map[x - self.lo]).find(|&r| r != u32::MAX)? as usize;
        Some((r1.min(r2), r1.max(r2)))
    }
}

pub fn rank(placements: &[Placement], read: &Record, refs: &[Record], opt: &Options) -> Vec<f64> {
    let n = placements.len();
    if n < 2 {
        return vec![1.0; n];
    }
    let mut anchors: Vec<usize> = placements.iter().flat_map(|p| p.divergences.iter().copied()).collect();
    anchors.sort_unstable();
    anchors.dedup();
    let (mut sum, mut count) = (vec![0f64; n], vec![0u32; n]);
    let mut spans: Vec<Option<(bool, usize, usize, usize)>> = Vec::with_capacity(n);
    let mut seen: HashMap<(bool, &[u8]), u32> = HashMap::new();
    let mut order: Vec<(u32, usize)> = Vec::with_capacity(n);
    let mut rc = Vec::new();
    let last = read.n_tiles() - 1;
    for a in anchors {
        let (wl, wh) = (a.saturating_sub(opt.window), (a + opt.window).min(last));
        spans.clear();
        spans.extend(placements.iter().map(|p| {
            let (l, h) = p.span(wl, wh)?;
            Some((p.minus, p.target, l.saturating_sub(opt.pad), (h + opt.pad).min(refs[p.target].n_tiles() - 1)))
        }));
        let covered = spans.iter().filter(|s| s.is_some()).count();
        if covered == 0 {
            continue;
        }
        if covered == n {
            let (m0, t0, l0, h0) = spans[0].unwrap();
            if spans[1..].iter().all(|s| {
                let (m, t, l, h) = s.unwrap();
                m == m0 && t == t0 && refs[t].tiles[l..=h] == refs[t0].tiles[l0..=h0]
            }) {
                continue;
            }
        }
        let fwd = read.bases(wl..=wh);
        revcomp(fwd, &mut rc);
        seen.clear();
        order.clear();
        for (k, s) in spans.iter().enumerate() {
            if let Some((minus, t, l, h)) = *s {
                let rb = refs[t].bases(l..=h);
                let d = *seen.entry((minus, rb)).or_insert_with(|| distance(if minus { &rc } else { fwd }, rb));
                order.push((d, k));
            }
        }
        order.sort_unstable();
        let mut i = 0;
        while i < order.len() {
            let mut j = i;
            while j + 1 < order.len() && order[j + 1].0 == order[i].0 {
                j += 1;
            }
            let r = (i + j) as f64 / 2.0 + 1.0;
            for &(_, k) in &order[i..=j] {
                sum[k] += r;
                count[k] += 1;
            }
            i = j + 1;
        }
        let r = (covered + 1 + n) as f64 / 2.0;
        for (k, s) in spans.iter().enumerate() {
            if s.is_none() {
                sum[k] += r;
                count[k] += 1;
            }
        }
    }
    (0..n).map(|k| if count[k] > 0 { sum[k] / count[k] as f64 } else { f64::INFINITY }).collect()
}
