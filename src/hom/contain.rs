use super::{Ctx, Evidence, put, repeat::Cover};
use std::collections::HashMap;
use std::io::Write;

pub const HEADER: &[u8] = b"#C\tcontig\tlength\tn_tiles\tdepth\tcontained(host:covered_fraction:left_flank_bp:right_flank_bp:reciprocal|oneway:depth_ratio_to_host)\tleft_link\tright_link\n";
pub type Contained = HashMap<usize, (usize, f64, usize, usize)>;

fn covered(iv: &mut [(u32, u32)]) -> (usize, usize, usize) {
    iv.sort_unstable();
    let (mut tot, mut end) = (0usize, 0u32);
    for &(s, e) in iv.iter() {
        let s = s.max(end);
        if s < e {
            tot += (e - s) as usize;
            end = e;
        }
    }
    (tot, iv.first().map_or(0, |x| x.0 as usize), end as usize)
}
pub fn interior_gap(iv: &mut [(u32, u32)], n: usize) -> (usize, usize, usize, usize) {
    let (tot, first, last) = covered(iv);
    if tot == 0 {
        return (n, 0, 0, 0);
    }
    (last - first - tot, first, n - last, tot)
}
pub fn holds(iv: &mut [(u32, u32)], n: usize, fr: f64) -> bool {
    let (gap, l, r, tot) = interior_gap(iv, n);
    tot > 0 && l + r < tot && gap as f64 <= fr * n as f64
}

pub fn compute(fm: &crate::Fmd, cov: &mut Cover, fr: f64) -> Contained {
    let mut best: HashMap<usize, (usize, f64, usize, usize, usize)> = HashMap::new();
    for ((c, h), iv) in cov.iter_mut() {
        let n = fm.rec_tiles(*c);
        let (gap, left, right, tot) = interior_gap(iv, n);
        let e = best.entry(*c).or_insert((*h, 0.0, 0, 0, n));
        let cf = tot as f64 / n as f64;
        if cf > e.1 || (cf == e.1 && *h < e.0) {
            *e = (*h, cf, left, right, gap);
        }
    }
    best.iter()
        .filter(|(c, (_, cf, l, r, gap))| {
            ((l + r) as f64) < *cf * fm.rec_tiles(**c) as f64
                && *gap as f64 <= fr * fm.rec_tiles(**c) as f64
        })
        .map(|(c, (h, cf, l, r, _))| (*c, (*h, *cf, *l, *r)))
        .collect()
}
pub fn reciprocal(contained: &Contained, c: usize) -> bool {
    contained
        .get(&c)
        .is_some_and(|(h, ..)| contained.get(h).is_some_and(|(hh, ..)| *hh == c))
}
pub fn ratio(ctx: &Ctx, ev: &Evidence, c: usize, h: usize) -> String {
    if ev.reads {
        let d = ev.mean(ctx, h);
        if d > 0.0 {
            format!("{:.2}", ev.mean(ctx, c) / d)
        } else {
            "-".to_string()
        }
    } else {
        "-".to_string()
    }
}
pub fn write(
    o: &mut impl Write,
    ctx: &Ctx,
    ev: &Evidence,
    contained: &Contained,
    ends: &HashMap<(usize, bool), String>,
) {
    for rec in 0..ctx.fm.n_refs() {
        let n = ctx.ntiles(rec);
        let d = if ev.reads {
            format!("{:.1}", ev.mean(ctx, rec))
        } else {
            "-".to_string()
        };
        let c = contained
            .get(&rec)
            .map_or("-".to_string(), |(h, cf, l, r)| {
                let lz = ctx.get(rec);
                format!(
                    "{}:{cf:.3}:{}:{}:{}:{}",
                    ctx.name(*h),
                    lz.start(*l).unwrap(),
                    lz.len() - lz.start(n - *r).unwrap(),
                    if reciprocal(contained, rec) {
                        "reciprocal"
                    } else {
                        "oneway"
                    },
                    ratio(ctx, ev, rec, *h)
                )
            });
        let side = |s: bool| {
            ends.get(&(rec, s))
                .cloned()
                .unwrap_or_else(|| "-".to_string())
        };
        put(
            o,
            format!(
                "C\t{}\t{}\t{n}\t{d}\t{c}\t{}\t{}\n",
                ctx.name(rec),
                ctx.len(rec),
                side(false),
                side(true)
            )
            .as_bytes(),
        );
    }
}
