pub mod chain;
pub mod contain;
pub mod family;
pub mod link;
pub mod pair;
pub mod read;
pub mod refs;
pub mod repeat;
pub mod tile;
pub mod trace;

use crate::tess::View;
use crate::tess::packed::Tile;
use crate::unit::BASE;
use crate::{Fmd, Tessellable};
use dashmap::DashMap;
use refs::Refs;
use std::io::Write;
use std::sync::atomic::{AtomicU32, Ordering::Relaxed};

pub const FIX: u64 = 1 << 10;
pub const EVIDENCE_HEADER: &[u8] = b"#evidence\tper_read_max_weight_union_v1\tdepth=max_covering_weight;agree=depth-max_mismatch_weight_per_read;alternatives=nonexclusive;ambiguous_placements=any_disagreement\n";
pub type Codes = DashMap<u64, Vec<(u64, u64)>>;

pub fn f(w: u64) -> f64 {
    w as f64 / FIX as f64
}
pub fn wsum(v: Option<dashmap::mapref::one::Ref<'_, u64, Vec<(u64, u64)>>>) -> u64 {
    v.map_or(0, |x| x.iter().map(|y| y.1).sum())
}
pub fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|c| match c {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}
pub fn put(out: &mut impl Write, b: &[u8]) {
    if let Err(e) = out.write_all(b) {
        if e.kind() == std::io::ErrorKind::BrokenPipe {
            std::process::exit(0);
        }
        panic!("cannot write output: {e}");
    }
}

pub fn pack(q: &Tile, r: &Tile) -> Option<u64> {
    let e = q.ex != r.ex;
    if q.span == r.span && q.phase() == r.phase() && !e {
        return None;
    }
    let c = |x: usize, b: u32| x.min((1 << b) - 1) as u64;
    Some(
        c(q.span as usize, 16) << 48
            | c(r.span as usize, 16) << 32
            | c(q.phase(), 8) << 24
            | c(r.phase(), 8) << 16
            | c(q.ex.len(), 7) << 9
            | c(r.ex.len(), 7) << 2
            | e as u64,
    )
}
pub fn packb(b: impl Iterator<Item = u8>, n: usize) -> u64 {
    let mut v = (n.min(31) as u64) << 58;
    for (i, c) in b.take(29).enumerate() {
        v |= (c as u64 & 3) << (2 * i);
    }
    v
}
pub fn codeb(p: u64) -> String {
    let n = (p >> 58) as usize;
    let s: String = (0..n.min(29))
        .map(|i| BASE[((p >> (2 * i)) & 3) as usize] as char)
        .collect();
    if n > 29 || n == 31 {
        format!("{s}..")
    } else {
        s
    }
}
pub fn code(p: u64) -> String {
    let (sq, sr, pq, pr, eq, er) = (
        p >> 48,
        p >> 32 & 0xffff,
        p >> 24 & 0xff,
        p >> 16 & 0xff,
        p >> 9 & 0x7f,
        p >> 2 & 0x7f,
    );
    let mut c = Vec::new();
    if sq != sr {
        c.push(format!("S:{sq}:{sr}"));
    }
    if pq != pr {
        c.push(format!("P:{pq}:{pr}"));
    }
    if p & 1 != 0 {
        c.push(format!("E:{eq}:{er}"));
    }
    c.join(",")
}

pub struct Ctx<'a> {
    pub fm: &'a Fmd,
    pub refs: &'a Refs,
}
impl<'a> Ctx<'a> {
    pub fn get(&self, rec: usize) -> View<'a> {self.refs.get(rec)}
    pub fn ntiles(&self, rec: usize) -> usize {self.refs.offs[rec + 1] - self.refs.offs[rec]}
    pub fn g0(&self, rec: usize) -> u64 {self.refs.offs[rec] as u64}
    pub fn len(&self, rec: usize) -> usize {self.get(rec).len()}
    pub fn name(&self, rec: usize) -> String {String::from_utf8_lossy(self.get(rec).name()).into_owned()}
}

pub struct Evidence {
    pub reads: bool,
    pub depth: Vec<AtomicU32>,
    pub depth_p: Vec<AtomicU32>,
    pub disagree_p: Vec<AtomicU32>,
    pub diffs: Codes,
    pub brks: Codes,
    pub links: link::Links,
}
impl Evidence {
    pub fn dep(&self, g: u64) -> u64 {
        self.depth
            .get(g as usize)
            .map_or(0, |x| x.load(Relaxed) as u64)
    }
    pub fn depp(&self, g: u64) -> u64 {
        self.depth_p
            .get(g as usize)
            .map_or(0, |x| x.load(Relaxed) as u64)
    }
    pub fn disagreement(&self, g: u64) -> u64 {
        self.disagree_p
            .get(g as usize)
            .map_or(0, |x| x.load(Relaxed) as u64)
    }
    pub fn agreement(&self, g: u64) -> u64 {
        self.depp(g).saturating_sub(self.disagreement(g))
    }
    pub fn mean(&self, ctx: &Ctx, rec: usize) -> f64 {
        let (g0, n) = (ctx.g0(rec), ctx.ntiles(rec));
        f((0..n).map(|t| self.dep(g0 + t as u64)).sum::<u64>()) / n.max(1) as f64
    }
}
