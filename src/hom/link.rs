use super::{Ctx, FIX, f, put, revcomp};
use dashmap::DashMap;
use std::collections::{BTreeMap, HashMap};
use std::io::Write;

pub const HEADER: &[u8] = b"#L\tcontigA\tendA\tcontigB\tendB\tweight\treads\tmean_junction(bp_between_contig_ends;negative=overlap)\tjunction_seq\n";
pub type Key = (usize, bool, usize, bool);
pub type Val = (u64, i64, u64, u64, String);
pub type Links = DashMap<Key, (u64, i64, u64, BTreeMap<usize, String>)>;
pub type End = (usize, bool, usize);
pub type Obs = (usize, bool, usize, bool, i64, u64, String);
pub type Row = (Key, Val);

pub(crate) fn observe_with(
    start: impl Fn(usize) -> usize,
    bases: impl Fn(std::ops::Range<usize>) -> Vec<u8>,
    prev: &[End],
    prev_e: usize,
    cin: &[End],
    s: usize,
    lk: &mut Vec<Obs>,
) {
    if prev.is_empty() || cin.is_empty() {
        return;
    }
    let (gb, ge) = (start(prev_e), start(s));
    let lw = FIX / (prev.len() * cin.len()) as u64;
    for &(a, sa, oa) in prev {
        for &(b, sb, ob) in cin {
            if (a, sa) != (b, sb) {
                let (js, je) = (gb + oa, ge.saturating_sub(ob));
                let jun = je as i64 - js as i64;
                let seq: String = if jun > 0 {
                    String::from_utf8(bases(js..je)).unwrap()
                } else {
                    String::new()
                };
                lk.push(if (a, sa) <= (b, sb) {
                    (a, sa, b, sb, jun, lw, seq)
                } else {
                    (b, sb, a, sa, jun, lw, revcomp(&seq))
                });
            }
        }
    }
}
pub fn add(links: &Links, lk: Vec<Obs>) {
    for (a, sa, b, sb, gap, lw, seq) in lk {
        let mut v = links
            .entry((a, sa, b, sb))
            .or_insert((0, 0, 0, BTreeMap::new()));
        v.0 += lw;
        v.1 += gap;
        v.2 += 1;
        match v.3.get_mut(&seq.len()) {
            Some(x) => {
                if seq < *x {
                    *x = seq;
                }
            }
            None => {
                v.3.insert(seq.len(), seq);
            }
        }
    }
}
pub fn rows(links: &Links) -> Vec<Row> {
    let mut rows: Vec<Row> = links
        .iter()
        .map(|x| {
            let (w, gap, n, seqs) = x.value();
            let mean = (*gap as f64 / *n as f64).max(0.0);
            let seq = seqs
                .iter()
                .min_by(|a, b| {
                    (a.0.abs_diff(mean.round() as usize), a.0)
                        .cmp(&(b.0.abs_diff(mean.round() as usize), b.0))
                })
                .map_or(String::new(), |x| x.1.clone());
            (*x.key(), (*w, *gap, *n, 0, seq))
        })
        .collect();
    rows.sort_by(|a, b| b.1.0.cmp(&a.1.0).then(a.0.cmp(&b.0)));
    rows
}
pub fn ends(ctx: &Ctx, rows: &[Row]) -> HashMap<(usize, bool), String> {
    let mut ends: HashMap<(usize, bool), (u64, String)> = HashMap::new();
    for &((a, sa, b, sb), (w, ..)) in rows {
        for (me, side, other, oside) in [(a, sa, b, sb), (b, sb, a, sa)] {
            let e = ends.entry((me, side)).or_insert((0, String::new()));
            if w > e.0 {
                *e = (
                    w,
                    format!(
                        "{}:{}:{:.1}",
                        ctx.name(other),
                        if oside { 'R' } else { 'L' },
                        f(w)
                    ),
                );
            }
        }
    }
    ends.into_iter().map(|(k, v)| (k, v.1)).collect()
}
pub fn write(o: &mut impl Write, ctx: &Ctx, rows: &[Row]) {
    for ((a, sa, b, sb), (w, gap, n, _, seq)) in rows {
        put(
            o,
            format!(
                "L\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\n",
                ctx.name(*a),
                if *sa { 'R' } else { 'L' },
                ctx.name(*b),
                if *sb { 'R' } else { 'L' },
                f(*w),
                n,
                gap / *n as i64,
                if seq.is_empty() { "-" } else { seq.as_str() }
            )
            .as_bytes(),
        );
    }
}
