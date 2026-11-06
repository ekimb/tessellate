use super::{Ctx, read::Query, repeat::Table};
use crate::map::{Record, chain as placement, extend, myers};
use crate::tess::{Tile, store::Tess, tq};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::io::{self, Write};
use std::ops::Range;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::time::Instant;

const SCORED_TOP: usize = 5;
const SR_TRIM: usize = 2;
const SR_RESEED_OCC: usize = 10;
const SR_CANDIDATES: usize = 16;
const SR_PAD: usize = 3;
pub const HEADER: &[u8] =b"#read\trank\tstatus\treference\tstrand\tquery_tiles\tquery_bases\treference_tiles\treference_bases\tcovered_tiles\tgeometry_matches\td_sites\td_matches\tscore\tbest_gap\tbest_ties\tanchors\tmean_rank\tmapq\n";

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub gap: usize,
    pub band: usize,
    pub top: usize,
    pub min_gap: Option<usize>,
    pub vote: Option<myers::Options>,
    pub sr: bool,
}
impl Default for Options {
    fn default() -> Self { Self { gap: 1000, band: 64, top: 5, min_gap: None, vote: None, sr: false } }
}

struct Vote { options: myers::Options }

#[derive(Default)]
struct Stats {
    seed: AtomicU64, project: AtomicU64, score: AtomicU64, chain: AtomicU64, write: AtomicU64, vote: AtomicU64,
    stage1: AtomicU64, settled: AtomicU64,
    anchors: AtomicU64, projected: AtomicU64, tiles: AtomicU64, candidates: AtomicU64, voted: AtomicU64,
}
impl Stats {
    fn time(counter: &AtomicU64, started: Instant) {counter.fetch_add(started.elapsed().as_nanos() as u64, Relaxed);}
    fn report(&self) -> String {
        let s = |c: &AtomicU64| c.load(Relaxed) as f64 / 1e9;
        let n = |c: &AtomicU64| c.load(Relaxed);
        format!("stage thread-seconds: seed {:.1}, stage-1 {:.1}, project {:.1}, score {:.1}, chain {:.1}, write {:.1}, vote {:.1}; {} reads settled before projection, {} anchors located, {} after projection, {} reference tiles scored, {} candidates, {} reads voted",
            s(&self.seed), s(&self.stage1), s(&self.project), s(&self.score), s(&self.chain), s(&self.write) - s(&self.vote), s(&self.vote),
            n(&self.settled), n(&self.anchors), n(&self.projected), n(&self.tiles), n(&self.candidates), n(&self.voted))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Anchor {
    rec: usize,
    rc: bool,
    q: usize,
    r: usize,
    len: usize,
}
impl Anchor {
    pub fn placement(q: Range<usize>, rec: usize, start: usize, rc: bool, nt: usize) -> Self {
        let len = q.len();
        Self { rec, rc, q: q.start, r: if rc {nt - start - len} else {start}, len }
    }
    fn qe(self) -> usize { self.q + self.len }
    fn re(self) -> usize { self.r + self.len }
    fn diagonal(self) -> i64 { self.r as i64 - self.q as i64 }
    fn native(self, nt: usize) -> Range<usize> {
        if self.rc {nt - self.re()..nt - self.r} else {self.r..self.re()}
    }
}
fn normalize(mut anchors: Vec<Anchor>) -> Vec<Anchor> {
    anchors.sort_unstable_by_key(|a| (a.rec, a.rc, a.diagonal(), a.q, a.len));
    let mut out: Vec<Anchor> = Vec::with_capacity(anchors.len());
    for a in anchors {
        if a.len == 0 {continue;}
        if let Some(last) = out.last_mut() {
            if (last.rec, last.rc, last.diagonal()) == (a.rec, a.rc, a.diagonal()) && a.q <= last.qe() {
                last.len = last.qe().max(a.qe()) - last.q;
                continue;
            }
        }
        out.push(a);
    }
    out
}

#[derive(Clone, Copy)]
struct Occurrence { group: usize, start: usize, end: usize, rc: bool }
#[derive(Default)]
struct Intervals { sites: Vec<Occurrence>, max_end: Vec<usize> }
impl Intervals {
    fn finish(&mut self) {
        self.sites.sort_unstable_by_key(|s| (s.start, s.end, s.group, s.rc));
        let mut end = 0;
        self.max_end = self.sites.iter().map(|s| {end = end.max(s.end); end}).collect();
    }
    fn overlapping(&self, range: Range<usize>) -> impl Iterator<Item = Occurrence> + '_ {
        let end = self.sites.partition_point(|s| s.start < range.end);
        let start = self.max_end[..end].partition_point(|&e| e <= range.start);
        self.sites[start..end].iter().copied().filter(move |s| s.end > range.start)
    }
}

pub struct Index<'a> {
    table: &'a Table,
    occurrences: Vec<Intervals>,
    differing: Vec<Vec<usize>>,
    lengths: Vec<usize>,
    units: Vec<Vec<u16>>,
    records: Vec<Record>,
    vote: Option<Vote>,
    stats: Stats,
    pub options: Options,
}
impl<'a> Index<'a> {
    pub fn new(ctx: &Ctx, table: &'a Table, options: Options) -> Self {
        let mut occurrences: Vec<_> = (0..ctx.fm.n_refs()).map(|_| Intervals::default()).collect();
        let mut differing = vec![Vec::new(); ctx.fm.n_refs()];
        for group in 0..table.n() {
            let sites = table.sites(group);
            for &(rec, start, rc) in sites {
                occurrences[rec as usize].sites.push(Occurrence {
                    group, start: start as usize, end: start as usize + table.len(group), rc,
                });
            }
            for &(copy, q, r, _) in table.diffs(group) {
                differing[sites[0].0 as usize].push(q as usize);
                differing[sites[copy as usize].0 as usize].push(r as usize);
            }
        }
        for index in &mut occurrences {index.finish();}
        for sites in &mut differing {sites.sort_unstable(); sites.dedup();}
        let records: Vec<Record> = (0..ctx.fm.n_refs()).into_par_iter()
            .map(|r| Record::read(ctx.refs, r).unwrap_or_else(|e| panic!("cannot decode reference {r}: {e}"))).collect();
        let units = records.iter().map(|r| r.units()).collect();
        let vote = options.vote.map(|options| Vote { options });
        Self { table, occurrences, differing, lengths: (0..ctx.fm.n_refs()).map(|r| ctx.ntiles(r)).collect(), units, records, vote, stats: Stats::default(), options }
    }
    pub fn scan(&self, ctx: &Ctx, paths: &[String], k: usize, occ: usize, output: &mut impl Write) -> io::Result<()> {
        if k == 0 {return Err(io::Error::other("chain seed length must be positive"));}
        let fm = ctx.fm;
        for path in paths {
            let started = Instant::now();
            let prefix = tq::stem(path).ok_or_else(|| io::Error::other("chain reads require a TQ input"))?;
            let aux = Tess::try_open(prefix)?;
            let (mut reads, mut smems, mut located, mut skipped) = (0usize, 0usize, 0usize, 0usize);
            tq::parallel(
                &[path],
                |first, batch| {
                    batch.into_par_iter().enumerate().map(|(j, (name, units))| {
                        let offset = *aux.offsets.get(first + j)
                            .ok_or_else(|| io::Error::other(format!("{path}: more records than its .tss")))?;
                        let query = Query::new(&name, units, aux.try_view_at(offset as usize)?)?;
                        let (mut nsm, mut nloc, mut nskip) = (0usize, 0usize, 0usize);
                        let mut anchors = Vec::new();
                        let mut longest = (0usize, 0usize);
                        let seeding = Instant::now();
                        for (s, e, size, places) in fm.smems_no_rescan(query.units.len(), |i| query.units[i], k, 1, occ, |r, p| self.units[r][p]) {
                            nsm += 1;
                            if e - s > longest.0 {longest = (e - s, size);}
                            if occ > 0 && size > occ {nskip += 1; continue;}
                            for pl in places {
                                nloc += 1;
                                let nt = self.lengths[pl.rec];
                                let start = if pl.rc {nt - pl.tile - (e - s)} else {pl.tile};
                                anchors.push(Anchor::placement(s..e, pl.rec, start, pl.rc, nt));
                            }
                            if self.options.sr {
                                let mut more = |a: usize, b: usize| {
                                    let mut bi = fm.all();
                                    for i in (a..b).rev() {
                                        bi = fm.extend_left(bi, query.units[i]);
                                        if bi.size == 0 {return;}
                                    }
                                    if bi.size > size && (occ == 0 || bi.size <= occ) {
                                        for pl in fm.positions(&bi, bi.size) {
                                            nloc += 1;
                                            let nt = self.lengths[pl.rec];
                                            let start = if pl.rc {nt - pl.tile - (b - a)} else {pl.tile};
                                            anchors.push(Anchor::placement(a..b, pl.rec, start, pl.rc, nt));
                                        }
                                    }
                                };
                                if e - s >= k + 2 * SR_TRIM {more(s + SR_TRIM, e - SR_TRIM);}
                                if size <= SR_RESEED_OCC && e - s >= 2 * k {
                                    let mid = s + (e - s) / 2;
                                    more(s, mid);
                                    more(mid, e);
                                }
                            }
                        }
                        Stats::time(&self.stats.seed, seeding);
                        self.stats.anchors.fetch_add(anchors.len() as u64, Relaxed);
                        let candidates = self.assign(&query, anchors)?;
                        let mut rows = Vec::new();
                        let writing = Instant::now();
                        self.write(&mut rows, ctx, &query, &candidates, longest.1)?;
                        Stats::time(&self.stats.write, writing);
                        Ok((nsm, nloc, nskip, rows))
                    }).collect::<io::Result<Vec<_>>>()
                },
                |_, batch| {
                    for (nsm, nloc, nskip, rows) in batch? {
                        output.write_all(&rows)?;
                        reads += 1;
                        smems += nsm;
                        located += nloc;
                        skipped += nskip;
                    }
                    eprint!("\r[tessellate hom]\t{path}: {reads} reads chained in {:.2}s", started.elapsed().as_secs_f64());
                    Ok(())
                },
            )?;
            if reads != aux.offsets.len() {
                return Err(io::Error::other(format!("{path}: {reads} records but its .tss has {}", aux.offsets.len())));
            }
            eprintln!("\n[tessellate hom]\t{path}: {reads} records, {smems} SMEMs, {located} occurrences located, {skipped} SMEMs above --occ {occ} skipped (0 = locate all) in {:.2}s",
                started.elapsed().as_secs_f64());
            eprintln!("[tessellate hom]\t{path}: {}", self.stats.report());
        }
        Ok(())
    }
    fn project(&self, anchors: Vec<Anchor>) -> Vec<Anchor> {
        let mut anchors = normalize(anchors);
        let mut groups: HashMap<(usize, bool, i64), Vec<Range<usize>>> = HashMap::new();
        for &a in &anchors {
            let native = a.native(self.lengths[a.rec]);
            for site in self.occurrences[a.rec].overlapping(native.clone()) {
                let s = native.start.max(site.start);
                let e = native.end.min(site.end);
                let q = a.q + if a.rc {native.end - e} else {s - native.start};
                let rc = a.rc ^ site.rc;
                let cs = if site.rc {site.end - e} else {s - site.start};
                let cr = if rc {site.end - site.start - cs - (e - s)} else {cs};
                groups.entry((site.group, rc, cr as i64 - q as i64)).or_default().push(q..q + e - s);
            }
        }
        for ((group, rc, diagonal), mut ranges) in groups {
            ranges.sort_unstable_by_key(|r| (r.start, r.end));
            let mut merged: Vec<Range<usize>> = Vec::new();
            for range in ranges {
                if let Some(last) = merged.last_mut() {
                    if range.start <= last.end {last.end = last.end.max(range.end); continue;}
                }
                merged.push(range);
            }
            for q in merged {
                let start = (q.start as i64 + diagonal) as usize;
                let cs = if rc {self.table.len(group) - start - q.len()} else {start};
                for &(rec, native, reversed) in self.table.sites(group) {
                    let native = native as usize + if reversed {self.table.len(group) - cs - q.len()} else {cs};
                    anchors.push(Anchor::placement(q.clone(), rec as usize, native, rc ^ reversed, self.lengths[rec as usize]));
                }
            }
        }
        normalize(anchors)
    }

    fn score(&self, qt: &[Tile], mir: &[Tile], anchors: Vec<Anchor>) -> io::Result<Vec<Scored>> {
        let mut scored = Vec::with_capacity(anchors.len());
        for a in anchors {
            let native = a.native(self.lengths[a.rec]);
            let rt = &self.records[a.rec].tiles[native.clone()];
            let sites = &self.differing[a.rec];
            let ds = sites.partition_point(|&s| s < native.start);
            let de = sites.partition_point(|&s| s < native.end);
            let mut marks = vec![false; a.len];
            for &s in &sites[ds..de] {marks[s - native.start] = true;}
            let mut prefix = vec![Counts::default()];
            for j in 0..a.len {
                let rj = if a.rc {a.len - 1 - j} else {j};
                let tile = if a.rc {&mir[a.q + j]} else {&qt[a.q + j]};
                debug_assert_eq!(tile.id.canon(), rt[rj].id.canon());
                let matches = tile == &rt[rj];
                let prev = *prefix.last().unwrap();
                prefix.push(Counts {covered: prev.covered + 1, geometry: prev.geometry + usize::from(matches),
                    d_sites: prev.d_sites + usize::from(marks[rj]), d_matches: prev.d_matches + usize::from(marks[rj] && matches)});
            }
            self.stats.tiles.fetch_add(a.len as u64, Relaxed);
            scored.push(Scored {anchor: a, prefix});
        }
        Ok(scored)
    }

    pub fn assign(&self, query: &Query<'_>, anchors: Vec<Anchor>) -> io::Result<Vec<Candidate>> {
        let qt = query.tiles(0..query.units.len());
        let mir: Vec<_> = qt.iter().map(|t| t.mirror()).collect();
        if let Some(g) = self.options.min_gap {
            let staging = Instant::now();
            let normalized = normalize(anchors.clone());
            let first = chain(normalized.iter().copied().map(Scored::light).collect(), self.options);
            let settled = match (first.first(), first.get(1)) {
                (Some(best), Some(next)) => best.score - next.score >= g as i64,
                (Some(_), None) => true,
                (None, _) => false,
            };
            if settled {
                let top = &first[..first.len().min(SCORED_TOP)];
                let inside = |a: &Anchor| top.iter().any(|c| a.rec == c.rec && a.rc == c.rc
                    && a.q >= c.query.start && a.qe() <= c.query.end && a.r >= c.reference.start && a.re() <= c.reference.end);
                let kept: Vec<Anchor> = normalized.into_iter().filter(inside).collect();
                let scored = chain(self.score(&qt, &mir, kept)?, self.options);
                let result = if scored.is_empty() {first} else {scored};
                Stats::time(&self.stats.stage1, staging);
                self.stats.settled.fetch_add(1, Relaxed);
                self.stats.candidates.fetch_add(result.len() as u64, Relaxed);
                return Ok(result);
            }
            Stats::time(&self.stats.stage1, staging);
        }
        let projecting = Instant::now();
        let anchors = self.project(anchors);
        Stats::time(&self.stats.project, projecting);
        self.stats.projected.fetch_add(anchors.len() as u64, Relaxed);
        if anchors.is_empty() {return Ok(Vec::new());}
        let scoring = Instant::now();
        let scored = self.score(&qt, &mir, anchors)?;
        Stats::time(&self.stats.score, scoring);
        let chaining = Instant::now();
        let candidates = chain(scored, self.options);
        Stats::time(&self.stats.chain, chaining);
        self.stats.candidates.fetch_add(candidates.len() as u64, Relaxed);
        Ok(candidates)
    }

    fn rank(&self, vote: &Vote, query: &Query<'_>, selected: &[&Candidate]) -> Vec<f64> {
        let voting = Instant::now();
        self.stats.voted.fetch_add(1, Relaxed);
        let n = query.units.len();
        let read = Record::new(query.aux.name().to_vec(), query.tiles(0..n));
        let placements: Vec<myers::Placement> = selected.iter().map(|c| {
            let native = if c.rc {self.lengths[c.rec] - c.reference.end} else {c.reference.start};
            let q = if c.rc {n - c.query.end} else {c.query.start};
            let anchor = placement::Anchor {q: q as u32, r: native as u32, len: 1};
            let chain = placement::Chain {rec: c.rec, minus: c.rc, anchors: vec![anchor]};
            myers::Placement::from_chain(&chain, &read, extend::extend(&chain, &query.units, &self.units[c.rec]))
        }).collect();
        let ranks = myers::rank(&placements, &read, &self.records, &vote.options);
        Stats::time(&self.stats.vote, voting);
        ranks
    }

    fn distance(&self, read: &[u8], n: usize, c: &Candidate) -> u32 {
        let nt = self.lengths[c.rec];
        let start = c.reference.start as i64 - c.query.start as i64 - SR_PAD as i64;
        let ws = start.max(0) as usize;
        let we = ((start + (n + 2 * SR_PAD) as i64).max(0) as usize).min(nt);
        if ws >= we {return u32::MAX;}
        let (lo, hi) = if c.rc {(nt - we, nt - ws - 1)} else {(ws, we - 1)};
        let bases = self.records[c.rec].bases(lo..=hi);
        if c.rc {
            let rc: Vec<u8> = bases.iter().rev().map(|&b| match b {b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', _ => b'N'}).collect();
            myers::distance(read, &rc)
        } else {
            myers::distance(read, bases)
        }
    }

    fn write_sr(&self, out: &mut impl Write, ctx: &Ctx, query: &Query<'_>, candidates: &[Candidate], occ: usize) -> io::Result<()> {
        let name = String::from_utf8_lossy(query.aux.name());
        if candidates.is_empty() {
            return writeln!(out, "{name}\t0\tunmapped\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t*\t0\t0\t*\t0");
        }
        let n = query.units.len();
        let pool = &candidates[..candidates.len().min(SR_CANDIDATES)];
        let read: Vec<u8> = query.tiles(0..n).iter().flat_map(|t| t.render()).collect();
        let dist: Vec<u32> = if pool.len() > 1 {pool.iter().map(|c| self.distance(&read, n, c)).collect()} else {vec![0]};
        let mut order: Vec<usize> = (0..pool.len()).collect();
        order.sort_by_key(|&i| (dist[i], std::cmp::Reverse(pool[i].score), i));
        let d1 = dist[order[0]];
        let ties = order.iter().take_while(|&&i| dist[i] == d1).count();
        let best = &pool[order[0]];
        let gap = order.get(1).map(|&i| best.score.saturating_sub(pool[i].score).max(0));
        let cover = mapq(best.score, None, 1, best.anchors, Some(n));
        let cap = if occ > SR_RESEED_OCC {30} else {60};
        let top_score = pool.iter().map(|c| c.score).max().unwrap_or(0);
        let tied = pool.iter().filter(|c| c.score == top_score).count() > 1;
        let q = match order.get(1) {
            _ if tied => 0,
            None => cover.min(cap),
            Some(_) if ties > 1 => 0,
            Some(&i) => (30 * dist[i].saturating_sub(d1).saturating_sub(1)).min(cap as u32).min(cover as u32) as u8,
        };
        let top = if self.options.top == 0 {order.len()} else {self.options.top.min(order.len())}.max(ties);
        let gap = gap.map(|g| g.to_string()).unwrap_or_else(|| "*".to_owned());
        for (rank, &i) in order.iter().take(top).enumerate() {
            let c = &pool[i];
            let reference = ctx.get(c.rec);
            let r = if c.rc {self.lengths[c.rec] - c.reference.end..self.lengths[c.rec] - c.reference.start} else {c.reference.clone()};
            let strand = if c.rc {'-'} else {'+'};
            let status = if dist[i] != d1 {"secondary"} else if ties > 1 {"ambiguous"} else {"best"};
            let mq = if rank == 0 {q} else {0};
            writeln!(out, "{name}\t{}\t{status}\t{}\t{strand}\t{}..{}\t{}..{}\t{}..{}\t{}..{}\t{}\t{}\t{}\t{}\t{}\t{gap}\t{ties}\t{}\t{}\t{mq}",
                rank + 1, String::from_utf8_lossy(reference.name()), c.query.start, c.query.end, query.start(c.query.start), query.start(c.query.end),
                r.start, r.end, self.records[c.rec].tile_start(r.start), self.records[c.rec].tile_start(r.end), c.counts.covered, c.counts.geometry,
                c.counts.d_sites, c.counts.d_matches, c.score, c.anchors, dist[i])?;
        }
        Ok(())
    }

    pub fn write(&self, out: &mut impl Write, ctx: &Ctx, query: &Query<'_>, candidates: &[Candidate], occ: usize) -> io::Result<()> {
        if self.options.sr {return self.write_sr(out, ctx, query, candidates, occ);}
        let name = String::from_utf8_lossy(query.aux.name());
        let Some(best) = candidates.first() else {
            return writeln!(out, "{name}\t0\tunmapped\t*\t*\t*\t*\t*\t*\t0\t0\t0\t0\t0\t*\t0\t0\t*\t0");
        };
        let ties = candidates.iter().take_while(|c| c.score == best.score).count();
        let best_gap = candidates.get(1).map(|c| (best.score - c.score) as usize);
        let gap = best_gap.map(|g| g.to_string()).unwrap_or_else(|| "*".to_owned());
        let limit = if self.options.min_gap.zip(best_gap).is_some_and(|(g, gap)| gap >= g) {1} else {usize::MAX};
        let cutoff = if self.options.top == 0 {i64::MIN} else {
            candidates.get(self.options.top - 1).map_or(i64::MIN, |c| c.score)
        };
        let selected: Vec<&Candidate> = candidates.iter().take(limit).take_while(|c| c.score >= cutoff).collect();
        let ranks = match &self.vote {
            Some(vote) if selected.len() > 1 => Some(self.rank(vote, query, &selected)),
            _ => None,
        };
        let mut order: Vec<usize> = (0..selected.len()).collect();
        if let Some(ranks) = &ranks {
            order.sort_by(|&a, &b| ranks[a].total_cmp(&ranks[b]).then(a.cmp(&b)));
        }
        let top_ties = ranks.as_ref().map(|r| order.iter().take_while(|&&i| r[i] == r[order[0]]).count());
        for (rank, &i) in order.iter().enumerate() {
            let c = selected[i];
            let reference = ctx.get(c.rec);
            let r = if c.rc {self.lengths[c.rec] - c.reference.end..self.lengths[c.rec] - c.reference.start} else {c.reference.clone()};
            let strand = if c.rc {'-'} else {'+'};
            let status = match top_ties {
                Some(t) => if rank > 0 {"secondary"} else if t > 1 {"ambiguous"} else {"best"},
                None => if c.score == best.score {if ties > 1 {"ambiguous"} else {"best"}} else {"secondary"},
            };
            let mean = ranks.as_ref().map_or_else(|| "*".to_owned(), |r| format!("{:.3}", r[i]));
            let q = if status == "best" {mapq(c.score, best_gap, ties, c.anchors, self.options.sr.then_some(query.units.len()))} else {0};
            writeln!(out, "{name}\t{}\t{status}\t{}\t{strand}\t{}..{}\t{}..{}\t{}..{}\t{}..{}\t{}\t{}\t{}\t{}\t{}\t{gap}\t{ties}\t{}\t{mean}\t{q}",
                rank + 1, String::from_utf8_lossy(reference.name()), c.query.start, c.query.end, query.start(c.query.start), query.start(c.query.end),
                r.start, r.end, self.records[c.rec].tile_start(r.start), self.records[c.rec].tile_start(r.end), c.counts.covered, c.counts.geometry,
                c.counts.d_sites, c.counts.d_matches, c.score, c.anchors)?;
        }
        Ok(())
    }
}

fn mapq(score: i64, gap: Option<usize>, ties: usize, anchors: usize, tiles: Option<usize>) -> u8 {
    let gap = gap.unwrap_or(score.max(0) as usize);
    let s1 = score.max(1) as f64;
    let pen = match tiles {
        Some(n) => (s1 / (2.0 * n.max(1) as f64)).min(1.0),
        None => (anchors as f64 / 10.0).min(s1 / 100.0).min(1.0),
    };
    let q = pen * 40.0 * (gap as f64).min(s1) / s1 * s1.ln() - 4.343 * (ties.max(1) as f64).ln();
    q.clamp(0.0, 60.0) as u8
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts { pub covered: usize, pub geometry: usize, pub d_sites: usize, pub d_matches: usize }
impl Counts {
    fn score(self) -> i64 { (self.covered + self.geometry) as i64 }
    fn add(self, b: Self) -> Self {Self {covered: self.covered + b.covered, geometry: self.geometry + b.geometry, d_sites: self.d_sites + b.d_sites, d_matches: self.d_matches + b.d_matches}}
    fn sub(self, b: Self) -> Self {Self {covered: self.covered - b.covered, geometry: self.geometry - b.geometry, d_sites: self.d_sites - b.d_sites, d_matches: self.d_matches - b.d_matches}}
}
struct Scored { anchor: Anchor, prefix: Vec<Counts> }
impl Scored {
    fn light(anchor: Anchor) -> Self { Self {anchor, prefix: Vec::new()} }
    fn suffix(&self, trim: usize) -> Counts {
        if self.prefix.is_empty() {
            let n = self.anchor.len - trim;
            return Counts {covered: n, geometry: n, d_sites: 0, d_matches: 0};
        }
        self.prefix[self.anchor.len].sub(self.prefix[trim])
    }
}
#[derive(Debug)]
pub struct Candidate {
    pub rec: usize,
    pub rc: bool,
    pub query: Range<usize>,
    pub reference: Range<usize>,
    pub score: i64,
    pub counts: Counts,
    pub anchors: usize,
}

fn chain(mut anchors: Vec<Scored>, options: Options) -> Vec<Candidate> {
    anchors.sort_unstable_by_key(|s| {let a = s.anchor; (a.rec, a.rc, a.qe(), a.re(), a.q, a.r)});
    let n = anchors.len();
    let mut scores = vec![0i64; n];
    let mut parents = vec![None; n];
    let mut trims = vec![0; n];
    let mut history: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
    let mut record = None;
    for i in 0..n {
        let a = anchors[i].anchor;
        if record != Some((a.rec, a.rc)) {history.clear(); record = Some((a.rec, a.rc));}
        scores[i] = anchors[i].suffix(0).score();
        let d = a.diagonal();
        let band = options.band.min(i64::MAX as usize) as i64;
        for (_, previous) in history.range(d.saturating_sub(band)..=d.saturating_add(band)) {
            let first = previous.partition_point(|&j| anchors[j].anchor.qe() < a.q.saturating_sub(options.gap));
            for &j in &previous[first..] {
                let p = anchors[j].anchor;
                if p.q >= a.q || p.r >= a.r || p.qe() >= a.qe() || p.re() >= a.re() {continue;}
                let qgap = a.q.saturating_sub(p.qe());
                let rgap = a.r.saturating_sub(p.re());
                if qgap > options.gap || rgap > options.gap {continue;}
                let trim = p.qe().saturating_sub(a.q).max(p.re().saturating_sub(a.r));
                let penalty = (d - p.diagonal()).abs() + (qgap.min(rgap) / 10) as i64 + if qgap.max(rgap) > 0 {2} else {0};
                let score = scores[j] + anchors[i].suffix(trim).score() - penalty;
                if score > scores[i] {
                    scores[i] = score; parents[i] = Some(j); trims[i] = trim;
                }
            }
        }
        history.entry(d).or_default().push(i);
    }
    let mut order: Vec<_> = (0..n).collect();
    order.sort_unstable_by_key(|&i| (std::cmp::Reverse(scores[i]), i));
    let mut used = vec![false; n];
    let mut candidates: Vec<Candidate> = Vec::new();
    let keep = options.top.max(2);
    for end in order {
        if options.top > 0 && candidates.len() >= keep && scores[end] < candidates[keep - 1].score {break;}
        if used[end] {continue;}
        let mut path = vec![end];
        while let Some(parent) = parents[*path.last().unwrap()] {path.push(parent);}
        path.reverse();
        let mut counts = Counts::default();
        let mut shared = 0;
        for &i in &path {
            let extra = anchors[i].suffix(trims[i]);
            counts = counts.add(extra);
            if used[i] {shared += extra.covered;}
        }
        if shared * 5 >= counts.covered * 4 {continue;}
        for &i in &path {used[i] = true;}
        let first = anchors[path[0]].anchor;
        let last = anchors[end].anchor;
        candidates.push(Candidate {rec: first.rec, rc: first.rc, query: first.q..last.qe(), reference: first.r..last.re(),
            score: scores[end], counts, anchors: path.len()});
    }
    candidates
}