use super::{
    Ctx, Evidence, code, f,
    repeat::{Diff, Site},
};
use crate::tess::View;
use std::collections::HashSet;

pub const HEADER: &[u8] = b"#D\tcontigA\ttileA\tbaseA\tspanA\tcontigB\ttileB\tbaseB\tspanB\tstrand\tcode[\tdepthA\tagreeA\taltsA(span:weight)\tdepthB\tagreeB\taltsB\tverdict]\n";

type Dist = (u64, u64, Vec<(usize, u64)>);

pub struct Pairs {
    dkeys: Vec<u64>,
}
impl Pairs {
    pub fn new(ev: &Evidence) -> Pairs {
        let mut dkeys: Vec<u64> = ev.diffs.iter().map(|x| *x.key()).collect();
        dkeys.sort_unstable();
        Pairs { dkeys }
    }
    fn dist(ev: &Evidence, g: u64) -> Dist {
        let d = ev.depp(g);
        let mut v: Vec<(usize, u64)> = Vec::new();
        if let Some(e) = ev.diffs.get(&g) {
            for &(p, w) in e.iter() {
                let sp = ((p >> 48) & 0xffff) as usize;
                match v.iter_mut().find(|x| x.0 == sp) {
                    Some(x) => x.1 += w,
                    None => v.push((sp, w)),
                }
            }
        }
        (d, ev.agreement(g), v)
    }
    fn alts(v: &[(usize, u64)]) -> String {
        let mut v = v.to_vec();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        if v.is_empty() {
            return "-".to_string();
        }
        v.iter()
            .take(4)
            .map(|y| format!("{}:{:.1}", y.0, f(y.1)))
            .collect::<Vec<_>>()
            .join(",")
    }
    fn verdict(two: bool, sa: usize, sb: usize, a: &Dist, b: &Dist) -> (&'static str, usize) {
        if !two || a.0 == 0 || b.0 == 0 {
            return ("-", 0);
        }
        let w = |d: &Dist, own: usize, x: usize| {
            if x == own {
                d.1
            } else {
                d.2.iter().find(|y| y.0 == x).map_or(0, |y| y.1)
            }
        };
        let alt = |d: &Dist, own: usize| {
            d.2.iter()
                .filter(|y| y.0 != own)
                .max_by(|x, y| x.1.cmp(&y.1).then(y.0.cmp(&x.0)))
                .map_or((own, 0), |y| (y.0, y.1))
        };
        let (ca, cb) = if sa != sb {
            (sb, sa)
        } else {
            (alt(a, sa).0, alt(b, sb).0)
        };
        let (oa, xa, ob, xb) = (w(a, sa, sa), w(a, sa, ca), w(b, sb, sb), w(b, sb, cb));
        if oa == xa || ob == xb {
            return ("-", 0);
        }
        match (xa > oa, xb > ob) {
            (false, false) => {
                if sa != sb {
                    ("het", 0)
                } else {
                    ("-", 0)
                }
            }
            (true, false) => ("hapA_error", ca),
            (false, true) => ("hapB_error", cb),
            (true, true) => {
                if sa != sb {
                    ("swap", 0)
                } else if ca == cb {
                    ("shared_error", ca)
                } else {
                    ("both_error", 0)
                }
            }
        }
    }
    fn columns(a: &Dist, b: &Dist, v: (&str, usize)) -> String {
        let vd = if v.1 > 0 {
            format!("{}:{}", v.0, v.1)
        } else {
            v.0.to_string()
        };
        format!(
            "\t{:.1}\t{:.1}\t{}\t{:.1}\t{:.1}\t{}\t{vd}",
            f(a.0),
            f(a.1),
            Self::alts(&a.2),
            f(b.0),
            f(b.1),
            Self::alts(&b.2)
        )
    }
    pub fn rows(
        &self,
        ctx: &Ctx,
        ev: &Evidence,
        out: &mut String,
        lza: &View<'_>,
        me: usize,
        s: usize,
        e: usize,
        sites: &[Site],
        diffs: &[Diff],
    ) {
        let na = String::from_utf8_lossy(lza.name()).into_owned();
        let two = sites.len() == 2;
        let ga = ctx.g0(me);
        let mut seen: HashSet<u32> = HashSet::new();
        for &(ci, qi, ri, cd) in diffs {
            let (hr, _, rc) = sites[ci as usize];
            let hr = hr as usize;
            let lzb = ctx.get(hr);
            let (qi, ri) = (qi as usize, ri as usize);
            let (sa, sb) = (
                lza.start(qi + 1).unwrap() - lza.start(qi).unwrap(),
                lzb.start(ri + 1).unwrap() - lzb.start(ri).unwrap(),
            );
            out.push_str(&format!(
                "D\t{na}\t{qi}\t{}\t{sa}\t{}\t{ri}\t{}\t{sb}\t{}\t{}",
                lza.start(qi).unwrap(),
                String::from_utf8_lossy(lzb.name()),
                lzb.start(ri).unwrap(),
                if rc { '-' } else { '+' },
                code(cd)
            ));
            if ev.reads {
                let (a, b) = (
                    Self::dist(ev, ga + qi as u64),
                    Self::dist(ev, ctx.g0(hr) + ri as u64),
                );
                let distinct = me != hr || qi != ri;
                out.push_str(&Self::columns(&a, &b, Self::verdict(two && distinct, sa, sb, &a, &b)));
            }
            out.push('\n');
            seen.insert(qi as u32);
        }
        if !ev.reads || !two {
            return;
        }
        let m = e - s;
        let (hr, hf, hc) = (sites[1].0 as usize, sites[1].1 as usize, sites[1].2);
        let gb = ctx.g0(hr);
        let lzb = ctx.get(hr);
        let nb = String::from_utf8_lossy(lzb.name()).into_owned();
        let dk = &self.dkeys;
        let mut cand: Vec<usize> = Vec::new();
        let (i0, i1) = (
            dk.partition_point(|&g| g < ga + s as u64),
            dk.partition_point(|&g| g < ga + e as u64),
        );
        for &g in &dk[i0..i1] {
            cand.push((g - ga) as usize - s);
        }
        let (j0, j1) = (
            dk.partition_point(|&g| g < gb + hf as u64),
            dk.partition_point(|&g| g < gb + (hf + m) as u64),
        );
        for &g in &dk[j0..j1] {
            let jb = (g - gb) as usize - hf;
            cand.push(if hc { m - 1 - jb } else { jb });
        }
        cand.sort_unstable();
        cand.dedup();
        for j in cand {
            let (qi, ri) = (s + j, if hc { hf + m - 1 - j } else { hf + j });
            if seen.contains(&(qi as u32)) || (me == hr && qi == ri) {
                continue;
            }
            let (sa, sb) = (
                lza.start(qi + 1).unwrap() - lza.start(qi).unwrap(),
                lzb.start(ri + 1).unwrap() - lzb.start(ri).unwrap(),
            );
            let (a, b) = (
                Self::dist(ev, ga + qi as u64),
                Self::dist(ev, gb + ri as u64),
            );
            let v = Self::verdict(true, sa, sb, &a, &b);
            if v.0 == "-" {
                continue;
            }
            out.push_str(&format!(
                "D\t{na}\t{qi}\t{}\t{sa}\t{nb}\t{ri}\t{}\t{sb}\t{}\t={}\n",
                lza.start(qi).unwrap(),
                lzb.start(ri).unwrap(),
                if hc { '-' } else { '+' },
                Self::columns(&a, &b, v)
            ));
        }
    }
}
