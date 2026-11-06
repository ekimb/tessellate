use super::{
    Ctx, Evidence,
    contain::{self, Contained},
    f, family, pack, packb,
    pair::Pairs,
};
use crate::fmd::{Fmd, Bi, Place};
use crate::tess::View;
use crate::unit::BASE;
use crate::utils::{read_u64, write_u64s};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;
pub const HEADER: &[u8] = b"#R\tcontig\ttile_start\ttile_end\tbase_start\tbase_end\tn_occ\tcopies(contig:tile:base:strand:ndiff[:depth:support:Lthrough/Lother/Lorphan:Rthrough/Rother/Rorphan])\n";
const MAGIC: &[u8; 5] = b"TREPZ";
pub type Site = (u32, u32, bool);
pub type Pair = (u32, u32, u32, u32, i64);
pub type Diff = (u32, u32, u32, u64);
pub type Cover = HashMap<(usize, usize), Vec<(u32, u32)>>;
pub type Offsets = HashMap<(usize, usize), Vec<i64>>;

#[derive(Clone, Copy, Debug, clap::Args)]
pub struct RepeatOptions {
    #[arg(short = 'k', long, default_value = "31")]
    pub k: usize,
    #[arg(short = 'O', long, default_value_t = 0)]
    pub occ: usize,
    #[arg(long, default_value_t = 0.05)]
    pub frac: f64,
}

fn sym(c: u8) -> u8 {
    BASE.iter().position(|&x| x == c).unwrap_or(0) as u8
}
fn flank(lz: &View<'_>, t: usize) -> [u64; 2] {
    let b = lz.tiles(t..t + 1, false).unwrap()[0].render();
    [
        packb(b.iter().map(|&c| sym(c)), b.len()),
        packb(b.iter().rev().map(|&c| sym(c) ^ 2), b.len()),
    ]
}

fn put_u32s(w: &mut impl Write, v: impl Iterator<Item = u32>) -> std::io::Result<()> {
    let mut buf: Vec<u8> = Vec::with_capacity(1 << 20);
    for x in v {
        buf.extend_from_slice(&x.to_le_bytes());
        if buf.len() >= 1 << 20 {
            w.write_all(&buf)?;
            buf.clear();
        }
    }
    w.write_all(&buf)
}
fn get_u32s(b: &mut &[u8], c: usize) -> Vec<u32> {
    let v = b[..4 * c]
        .chunks_exact(4)
        .map(|x| u32::from_le_bytes(x.try_into().unwrap()))
        .collect();
    *b = &b[4 * c..];
    v
}
#[derive(Clone, Copy)]
struct Candidate {
    len: usize,
    bi: Bi,
}
impl Candidate {
    fn new(len: usize, mut bi: Bi) -> Self {
        if bi.rev < bi.fwd {
            std::mem::swap(&mut bi.fwd, &mut bi.rev);
        }
        Self { len, bi }
    }
    fn key(&self) -> (usize, usize, usize) {
        (self.len, self.bi.fwd, self.bi.size)
    }
}

pub struct Table {
    pub k: usize,
    pub occ: usize,
    pub fr: f64,
    meta: (Vec<usize>, String),
    len: Vec<u32>,
    first: Vec<usize>,
    sites: Vec<Site>,
    dfirst: Vec<usize>,
    diffs: Vec<Diff>,
    pairs: Vec<Pair>,
    pub contained: Contained,
    pub fams: Vec<Vec<usize>>,
}

impl Table {
    pub fn n(&self) -> usize {self.len.len()}
    pub fn n_sites(&self) -> usize {self.sites.len()}
    pub fn len(&self, i: usize) -> usize {self.len[i] as usize}
    pub fn sites(&self, i: usize) -> &[Site] {&self.sites[self.first[i]..self.first[i + 1]]}
    pub fn n_diffs(&self) -> usize {self.diffs.len()}
    pub fn diffs(&self, i: usize) -> &[Diff] {&self.diffs[self.dfirst[i]..self.dfirst[i + 1]]}
    pub fn dt(&self, i: usize, c: usize) -> Vec<u32> {
        let mut v: Vec<u32> = self
            .diffs(i)
            .iter()
            .filter(|d| c == 0 || d.0 as usize == c)
            .map(|d| if c == 0 { d.1 } else { d.2 })
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    }
    pub fn build(fm: &Fmd, refs: &[View<'_>], repeats: RepeatOptions) -> Table {
        debug_assert_eq!(refs.len(), fm.n_refs());
        assert!(repeats.k > 0, "minimum repeat length must be positive");
        assert!(repeats.frac.is_finite() && repeats.frac >= 0.0, "repeat fraction must be finite and nonnegative");
        let discovery = Instant::now();
        let k = repeats.k;
        let occ = repeats.occ;
        let fr = repeats.frac;
        let mut t = Table {k, occ, fr, meta: (fm.lengths().to_vec(), fm.source.clone()),
            len: Vec::new(),
            first: Vec::new(),
            sites: Vec::new(),
            dfirst: Vec::new(),
            diffs: Vec::new(),
            pairs: Vec::new(),
            contained: Contained::new(),
            fams: Vec::new(),
        };
        let found: Vec<_> = refs.par_iter().enumerate().map(|(rec, reference)| {
            let query = reference.units().expect("cannot read reference units");
            let (mut candidates, mut pairs) = (Vec::new(), Vec::new());
            for (s, e, bi) in fm.smems(query.len(), |i| query[i], k, 2) {
                if occ > 0 && bi.size > occ {continue;}
                let len = e - s;
                candidates.push(Candidate::new(len, bi));
                if fm.n_refs() == 1 {continue;}
                for row in bi.rows() {
                    let Place { rec: trec, rc, tile } = fm.locate_place(row);
                    let ts = if rc { fm.rec_tiles(trec) - tile - len } else { tile };
                    if rec == trec {continue;}
                    let off = s as i64 - ts as i64;
                    pairs.push((rec as u32, trec as u32, s as u32, e as u32, if rc {i64::MIN} else {off}));
                    pairs.push((trec as u32, rec as u32, ts as u32, (ts + len) as u32, if rc {i64::MIN} else {-off}));
                }
            }
            (candidates, pairs)
        }).collect();
        let mut candidates = Vec::with_capacity(found.iter().map(|f| f.0.len()).sum());
        let mut pairs = Vec::with_capacity(found.iter().map(|f| f.1.len()).sum());
        for (cs, ps) in found { candidates.extend(cs); pairs.extend(ps); }
        let raw = candidates.len();
        let np = pairs.len();
        candidates.par_sort_unstable_by_key(Candidate::key);
        candidates.dedup_by_key(|c| c.key());
        pairs.par_sort_unstable();
        pairs.dedup();
        t.pairs = pairs;
        let nr = candidates.len();
        let ns: usize = candidates.iter().map(|c| c.bi.size).sum();
        eprintln!(
            "[tessellate repeats]\t{raw} discoveries deduplicated to {nr} repeats, {ns} sites, {np} contig pair intervals in {:.2}s",
            discovery.elapsed().as_secs_f64()
        );
        let assembly = Instant::now();
        eprintln!("[tessellate repeats]\tconstructing repeat occurrences and tile differences");

        t.first.push(0);
        t.dfirst.push(0);
        for batch in candidates.chunks(8192) {
            let repeats: Vec<_> = batch
                .par_iter()
                .map(|c| Self::materialize(fm, refs, c))
                .collect();
            for (m, cs, dt) in repeats {
                t.len.push(m);
                t.sites.extend(cs);
                t.first.push(t.sites.len());
                t.diffs.extend(dt);
                t.dfirst.push(t.diffs.len());
            }
            if t.n() % (8192 * 16) == 0 || t.n() == nr {
                eprintln!(
                    "[tessellate repeats]\trepeat table: {}/{nr} repeats constructed in {:.2}s",
                    t.n(),
                    assembly.elapsed().as_secs_f64()
                );
            }
        }
        drop(candidates);
        eprintln!(
            "[tessellate repeats]\trepeat table assembled: {nr} repeats, {ns} sites, {} distinguishing tile pairs in {:.2}s",
            t.n_diffs(),
            assembly.elapsed().as_secs_f64()
        );
        let coverage = Instant::now();
        eprintln!("[tessellate repeats]\tcomputing contig pair coverage");
        let mut cov = t.coverage();
        eprintln!(
            "[tessellate repeats]\tcoverage computed for {} contig pairs in {:.2}s",
            cov.len(),
            coverage.elapsed().as_secs_f64()
        );
        let containment = Instant::now();
        eprintln!("[tessellate repeats]\tcomputing contained contigs");
        t.contained = contain::compute(fm, &mut cov, fr);
        eprintln!(
            "[tessellate repeats]\t{} contained contigs computed in {:.2}s",
            t.contained.len(),
            containment.elapsed().as_secs_f64()
        );
        let families = Instant::now();
        eprintln!("[tessellate repeats]\tgrouping contig families");
        t.fams = family::members(fm, &mut cov, fr);
        eprintln!(
            "[tessellate repeats]\t{} contig families grouped in {:.2}s",
            t.fams.len(),
            families.elapsed().as_secs_f64()
        );
        t
    }

    fn materialize(fm: &crate::Fmd, refs: &[View<'_>], candidate: &Candidate) -> (u32, Vec<Site>, Vec<Diff>) {
        let m = candidate.len;
        let mut copies: Vec<Site> = fm.positions(&candidate.bi, candidate.bi.size).iter().map(|p| {
                (
                    p.rec as u32,
                    if p.rc {
                        fm.rec_tiles(p.rec) - p.tile - m
                    } else {
                        p.tile
                    } as u32,
                    p.rc,
                )
            })
            .collect();
        copies.sort_unstable();
        if copies[0].2 {
            for copy in &mut copies {
                copy.2 = !copy.2;
            }
            copies.sort_unstable();
        }
        let (me, s, _) = copies[0];
        let s = s as usize;
        let lza = &refs[me as usize];
        let qt = lza.tiles(s..s + m, false).unwrap();
        let mut dt = Vec::new();
        for (ci, &(rec, rf, rc)) in copies.iter().enumerate().skip(1) {
            let lz = &refs[rec as usize];
            let rf = rf as usize;
            let rt = lz.tiles(rf..rf + m, false).unwrap();
            for j in 0..m {
                let (qi, ri) = (s + j, if rc { rf + m - 1 - j } else { rf + j });
                let q = if rc { qt[j].mirror() } else { qt[j].clone() };
                if let Some(code) = pack(&q, &rt[ri - rf]) {
                    dt.push((ci as u32, qi as u32, ri as u32, code));
                }
            }
        }
        (m as u32, copies, dt)
    }

    pub fn n_contained(&self) -> usize {
        self.contained.len()
    }
    pub fn save(&self, path: &str) -> std::io::Result<()> {
        let temporary = format!("{path}.tmp.{}", std::process::id());
        let file = std::fs::OpenOptions::new().write(true).create_new(true).open(&temporary)?;
        let result = self.write(std::io::BufWriter::new(file))
            .and_then(|()| std::fs::rename(&temporary, path));
        if result.is_err() { let _ = std::fs::remove_file(&temporary); }
        result
    }
    pub fn write(&self, mut w: impl Write) -> std::io::Result<()> {
        w.write_all(MAGIC)?;
        write_u64s(
            &mut w,
            [
                self.k,
                self.occ,
                self.fr.to_bits() as usize,
                self.len.len(),
                self.sites.len(),
                self.diffs.len(),
                self.pairs.len(),
                self.contained.len(),
                self.fams.len(),
            ],
        )?;
        put_u32s(&mut w, self.len.iter().cloned())?;
        write_u64s(&mut w, self.first.iter().cloned())?;
        put_u32s(
            &mut w,
            self.sites.iter().flat_map(|&(r, s, rc)| [r, s, rc as u32]),
        )?;
        write_u64s(&mut w, self.dfirst.iter().cloned())?;
        put_u32s(
            &mut w,
            self.diffs
                .iter()
                .flat_map(|&(c, q, r, code)| [c, q, r, code as u32, (code >> 32) as u32]),
        )?;
        put_u32s(
            &mut w,
            self.pairs
                .iter()
                .flat_map(|&(c, h, s, e, o)| [c, h, s, e, o as u32, (o >> 32) as u32]),
        )?;
        let mut cs: Vec<_> = self.contained.iter().collect();
        cs.sort_unstable_by_key(|x| *x.0);
        for (&c, &(h, cf, l, r)) in cs {
            write_u64s(&mut w, [c, h, cf.to_bits() as usize, l, r])?;
        }
        write_u64s(&mut w, self.fams.iter().map(|v| v.len()))?;
        put_u32s(&mut w, self.fams.iter().flatten().map(|&r| r as u32))?;
        let (lengths, source) = &self.meta;
        write_u64s(&mut w, [lengths.len(), source.len()])?;
        write_u64s(&mut w, lengths.iter().copied())?;
        w.write_all(source.as_bytes())?;
        w.flush()
    }
    pub fn read(path: &str) -> Table {
        let bytes = std::fs::read(path)
            .unwrap_or_else(|e| panic!("cannot read {path}: {e} (build it with tessellate hom --old -i <prefix>.tbwt)"));
        assert!(bytes.starts_with(MAGIC),
            "{path}: unsupported repeat table; remove it and rerun tessellate hom --old -i <prefix>.tbwt");
        let mut b = &bytes[MAGIC.len()..];
        let (k, occ, fr) = (
            read_u64(&mut b),
            read_u64(&mut b),
            f64::from_bits(read_u64(&mut b) as u64),
        );
        let (nr, ns, nd, np, nc, nf) = (
            read_u64(&mut b),
            read_u64(&mut b),
            read_u64(&mut b),
            read_u64(&mut b),
            read_u64(&mut b),
            read_u64(&mut b),
        );
        let len = get_u32s(&mut b, nr);
        let first: Vec<usize> = (0..=nr).map(|_| read_u64(&mut b)).collect();
        let sites = get_u32s(&mut b, 3 * ns)
            .chunks_exact(3)
            .map(|x| (x[0], x[1], x[2] != 0))
            .collect();
        let dfirst: Vec<usize> = (0..=nr).map(|_| read_u64(&mut b)).collect();
        let diffs = get_u32s(&mut b, 5 * nd)
            .chunks_exact(5)
            .map(|x| (x[0], x[1], x[2], x[3] as u64 | (x[4] as u64) << 32))
            .collect();
        let pairs = get_u32s(&mut b, 6 * np)
            .chunks_exact(6)
            .map(|x| {
                (
                    x[0],
                    x[1],
                    x[2],
                    x[3],
                    (x[4] as u64 | (x[5] as u64) << 32) as i64,
                )
            })
            .collect();
        let contained = (0..nc)
            .map(|_| {
                let c = read_u64(&mut b);
                (
                    c,
                    (
                        read_u64(&mut b),
                        f64::from_bits(read_u64(&mut b) as u64),
                        read_u64(&mut b),
                        read_u64(&mut b),
                    ),
                )
            })
            .collect();
        let sizes: Vec<usize> = (0..nf).map(|_| read_u64(&mut b)).collect();
        let members = get_u32s(&mut b, sizes.iter().sum());
        let (mut fams, mut at) = (Vec::with_capacity(nf), 0usize);
        for n in sizes {
            fams.push(members[at..at + n].iter().map(|&r| r as usize).collect());
            at += n;
        }
        let (nr, source_len) = (read_u64(&mut b), read_u64(&mut b));
        let lengths: Vec<usize> = (0..nr).map(|_| read_u64(&mut b)).collect();
        let source = String::from_utf8(b[..source_len].to_vec()).unwrap();
        b = &b[source_len..];
        assert!(b.is_empty(), "trailing data in repeat table");
        Table {
            k,
            occ,
            fr,
            meta: (lengths, source),
            len,
            first,
            sites,
            dfirst,
            diffs,
            pairs,
            contained,
            fams,
        }
    }

    pub fn matches_index(&self, fm: &crate::Fmd) -> bool {
        self.meta.0 == fm.lengths() && self.meta.1 == fm.source
    }

    pub fn coverage(&self) -> Cover {
        let mut cov: Cover = HashMap::new();
        for &(c, h, s, e, _) in &self.pairs {
            cov.entry((c as usize, h as usize))
                .or_default()
                .push((s, e));
        }
        cov
    }
    pub fn offsets(&self) -> Offsets {
        let mut offs: Offsets = HashMap::new();
        for &(c, h, _, _, o) in &self.pairs {
            if o != i64::MIN {
                offs.entry((c as usize, h as usize)).or_default().push(o);
            }
        }
        offs
    }

    pub fn annotate(&self, ctx: &Ctx, ev: &Evidence) -> Vec<String> {
        let pairs = Pairs::new(ev);
        (0..self.n()).into_par_iter().map(|i| {
                let (m, cs) = (self.len(i), self.sites(i));
                let (me, s) = (cs[0].0 as usize, cs[0].1 as usize);
                let e = s + m;
                let lza = ctx.get(me);
                let mut fl: HashMap<(u32, usize), [u64; 2]> = HashMap::new();
                let mut cols: Vec<String> = Vec::new();
                for (ci, &(rec, rf, rc)) in cs.iter().enumerate() {
                    let (rec, rf) = (rec as usize, rf as usize);
                    let lz = ctx.get(rec);
                    let dt = self.dt(i, ci);
                    let mut col = format!(
                        "{}:{rf}:{}:{}:{}",
                        String::from_utf8_lossy(lz.name()),
                        lz.start(rf).unwrap(),
                        if rc { '-' } else { '+' },
                        dt.len()
                    );
                    if ev.reads {
                        let (g0, l) = (ctx.g0(rec), ctx.ntiles(rec));
                        let d: u64 = (rf..rf + m).map(|t| ev.dep(g0 + t as u64)).sum();
                        let (mut sd, mut sa) = (0u64, 0u64);
                        for &t in &dt {
                            let g = g0 + t as u64;
                            sd += ev.depp(g);
                            sa += ev.agreement(g);
                        }
                        let sup = if sd > 0 {
                            format!("{:.2}", sa as f64 / sd as f64)
                        } else {
                            "-".to_string()
                        };
                        let mut side = |left: bool| -> String {
                            let (bt, dt_) = if left {
                                ((rf > 0).then(|| (g0 + rf as u64 - 1) << 1), g0 + rf as u64)
                            } else {
                                (
                                    (rf + m < l).then(|| (g0 + (rf + m) as u64) << 1 | 1),
                                    g0 + (rf + m) as u64 - 1,
                                )
                            };
                            let Some(bt) = bt else { return "-".to_string() };
                            let d = ev.depp(dt_);
                            if d == 0 {
                                return "-".to_string();
                            }
                            let (mut other, mut orphan, mut tot) = (0u64, 0u64, 0u64);
                            if let Some(v) = ev.brks.get(&bt) {
                                for &(p, w) in v.iter() {
                                    tot += w;
                                    let hit = cs.iter().enumerate().any(|(cj, &(recj, rfj, rcj))| {
                                        if cj == ci {return false;}
                                        let same = rcj == rc;
                                        let rfj = rfj as usize;
                                        let tile = if left == same {rfj.checked_sub(1)} else {
                                            (rfj + m < ctx.ntiles(recj as usize)).then_some(rfj + m)
                                        };
                                        let Some(tile) = tile else {return false;};
                                        let codes = fl.entry((recj, tile))
                                            .or_insert_with(|| flank(&ctx.get(recj as usize), tile));
                                        codes[usize::from(!same)] == p
                                    });
                                    if hit { other += w } else { orphan += w }
                                }
                            }
                            format!(
                                "{:.2}/{:.1}/{:.1}",
                                1.0 - tot.min(d) as f64 / d as f64,
                                f(other),
                                f(orphan)
                            )
                        };
                        col.push_str(&format!(
                            ":{:.1}:{sup}:{}:{}",
                            f(d) / m as f64,
                            side(true),
                            side(false)
                        ));
                    }
                    cols.push(col);
                }
                let mut out =
                    format!(
                        "R\t{}\t{s}\t{e}\t{}\t{}\t{}\t{}\n",
                        String::from_utf8_lossy(lza.name()),
                        lza.start(s).unwrap(),
                        lza.start(e).unwrap(),
                        cs.len(),
                        cols.join(";")
                    );
                    pairs.rows(ctx, ev, &mut out, &lza, me, s, e, cs, self.diffs(i));
                out
            })
            .collect()
    }
}
