use super::{Codes, Ctx, Evidence, FIX, link, pack, packb};
use crate::Unit;
use crate::tess::{Tile, store::Tess, tq};
use std::cell::OnceCell;
use std::collections::BinaryHeap;
mod query;
use dashmap::DashMap;
pub use query::Query;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering::Relaxed};
use std::time::Instant;
fn add_depth(intervals: &mut [(u64, u64, u32)], depth: &[AtomicU32]) {
    intervals.sort_unstable();
    let Some(&(start, _, weight)) = intervals.first() else { return };
    if intervals.iter().all(|x| x.2 == weight) {
        let mut end = 0;
        for &(start, stop, weight) in intervals.iter() {
            for tile in start.max(end)..stop {
                depth[tile as usize].fetch_add(weight, Relaxed);
            }
            end = end.max(stop);
        }
        return;
    }
    let (mut at, mut tile) = (0, start);
    let mut active: BinaryHeap<(u32, u64)> = BinaryHeap::new();
    loop {
        while at < intervals.len() && intervals[at].0 <= tile {
            let (_, stop, weight) = intervals[at];
            active.push((weight, stop));
            at += 1;
        }
        while active.peek().is_some_and(|&(_, stop)| stop <= tile) {
            active.pop();
        }
        let Some(&(weight, stop)) = active.peek() else {
            if at == intervals.len() { break; }
            tile = intervals[at].0;
            continue;
        };
        let next = intervals.get(at).map_or(stop, |x| stop.min(x.0));
        for t in tile..next {
            depth[t as usize].fetch_add(weight, Relaxed);
        }
        tile = next;
    }
}
fn add_disagreement(diffs: &[(u64, u64, u32)], disagreement: &[AtomicU32]) {
    let mut at = 0;
    while at < diffs.len() {
        let (tile, _, mut weight) = diffs[at];
        at += 1;
        while at < diffs.len() && diffs[at].0 == tile {
            weight = weight.max(diffs[at].2);
            at += 1;
        }
        disagreement[tile as usize].fetch_add(weight, Relaxed);
    }
}

pub fn scan_with_trace(
    ctx: &Ctx, paths: &[String], k: usize, occ: usize, ploidy: usize, gt: usize,
    trace: Option<&super::trace::Trace>,
) -> Evidence {
    let (fm, refs) = (ctx.fm, ctx.refs);
    let reads = !paths.is_empty();
    let nt = if reads { fm.n_tiles() } else { 0 };
    let ev = Evidence {
        reads,
        depth: std::iter::repeat_with(|| AtomicU32::new(0))
            .take(nt)
            .collect(),
        depth_p: std::iter::repeat_with(|| AtomicU32::new(0))
            .take(nt)
            .collect(),
        disagree_p: std::iter::repeat_with(|| AtomicU32::new(0))
            .take(nt)
            .collect(),
        diffs: DashMap::new(),
        brks: DashMap::new(),
        links: DashMap::new(),
    };
    let (located, skipped) = (AtomicUsize::new(0), AtomicUsize::new(0));
    for path in paths {
        let t1 = Instant::now();
        let mut stats = (0usize, 0usize, 0usize);
        let aux = tq::stem(path)
            .map(Tess::open)
            .expect("could not open store");
        tq::parallel(
            &[path],
            |first, batch| {
                batch
                    .into_par_iter()
                    .enumerate()
                    .map(|(j, (name, units))| {
                        let offset = *aux
                            .offsets
                            .get(first + j)
                            .unwrap_or_else(|| panic!("{path}: more records than its .tss"));
                        let query = Query::new(&name, units, aux.view_at(offset as usize))
                            .unwrap_or_else(|e| panic!("{path}: {e}"));
                        let nq = query.units.len();
                        let mut observations = Vec::new();
                        let source_name: Option<std::sync::Arc<str>> = trace.map(|_| path.as_str().into());
                        let read_name: Option<std::sync::Arc<str>> = trace.map(|_| String::from_utf8_lossy(&name).as_ref().into());
                        let (mut iv, mut ivp, mut dif, mut ends, mut nsm): (
                            Vec<(u64, u64, u32)>,
                            Vec<(u64, u64, u32)>,
                            Vec<(u64, u64, u32)>,
                            Vec<(u64, u64, u32)>,
                            usize,
                        ) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), 0);
                        let (mut prev, mut prev_e): (Vec<link::End>, usize) = (Vec::new(), 0);
                        let mut lk: Vec<link::Obs> = Vec::new();
                        for (s, e, bi) in fm.smems(nq, |i| query.units[i], k, 1) {
                            nsm += 1;
                            if occ > 0 && bi.size > occ {
                                skipped.fetch_add(1, Relaxed);
                                continue;
                            }
                            let w = (FIX / bi.size as u64).max(1) as u32;
                            let v = if ploidy == 0 { FIX as u32 } else if bi.size <= ploidy { w } else { 0 };
                            let fwd = OnceCell::new();
                            let mir = OnceCell::new();
                            let rev_units = OnceCell::new();
                            let flanks: [OnceCell<u64>; 4] =
                                std::array::from_fn(|_| OnceCell::new());
                            let rtile = |j: usize, rc: bool| {
                                let side = usize::from(j == e);
                                *flanks[side * 2 + usize::from(rc)]
                                    .get_or_init(|| query.flank(j, rc))
                            };
                            let (mut cin, mut cout): (Vec<link::End>, Vec<link::End>) =
                                (Vec::new(), Vec::new());
                            for pl in fm.positions(&bi, bi.size) {
                                located.fetch_add(1, Relaxed);
                                let view = OnceCell::new();
                                let reference = || view.get_or_init(|| refs.get(pl.rec));
                                let (m, l) = (e - s, ctx.ntiles(pl.rec));
                                let rf = if pl.rc { l - pl.tile - m } else { pl.tile };
                                let g0 = ctx.g0(pl.rec);
                                iv.push((g0 + rf as u64, g0 + (rf + m) as u64, w));
                                if s > 0 && (if pl.rc { rf + m + gt >= l } else { rf <= gt }) {
                                    let lz = reference();
                                    cin.push((
                                        pl.rec,
                                        pl.rc,
                                        if pl.rc {
                                            lz.len() - lz.start(rf + m).unwrap()
                                        } else {
                                            lz.start(rf).unwrap()
                                        },
                                    ));
                                }
                                if e < nq && (if pl.rc { rf <= gt } else { rf + m + gt >= l }) {
                                    let lz = reference();
                                    cout.push((
                                        pl.rec,
                                        !pl.rc,
                                        if pl.rc {
                                            lz.start(rf).unwrap()
                                        } else {
                                            lz.len() - lz.start(rf + m).unwrap()
                                        },
                                    ));
                                }
                                if v == 0 {
                                    continue;
                                }
                                if let Some(trace) = trace {
                                    let mut capture = |kind, region: usize, reference: std::ops::Range<usize>, query_range: std::ops::Range<usize>| {
                                        observations.push(super::trace::Observation {
                                            region, source: source_name.as_ref().unwrap().clone(), ordinal: first + j, offset,
                                            read: read_name.as_ref().unwrap().clone(), kind, rc: pl.rc,
                                            query_bases: query.start(query_range.start)..query.start(query_range.end),
                                            query: query_range, reference, smem_query: s..e,
                                            smem_reference: rf..rf+m, occurrences: bi.size, weight: v,
                                        });
                                    };
                                    for (region, r) in trace.overlapping(pl.rec, rf..rf+m) {
                                        let q = if pl.rc {s + rf + m - r.end} else {s + r.start - rf};
                                        let len = r.len(); capture("M", region, r, q..q+len);
                                    }
                                    let (before, after) = if pl.rc {(e < nq, s > 0)} else {(s > 0, e < nq)};
                                    if before && rf > 0 {
                                        let q = if pl.rc {e} else {s-1};
                                        for (region, r) in trace.overlapping(pl.rec, rf-1..rf) {capture("B<", region, r, q..q+1);}
                                    }
                                    if after && rf+m < l {
                                        let q = if pl.rc {s-1} else {e};
                                        for (region, r) in trace.overlapping(pl.rec, rf+m..rf+m+1) {capture("B>", region, r, q..q+1);}
                                    }
                                }
                                ivp.push((g0 + rf as u64, g0 + (rf + m) as u64, v));
                                let (before, after) = if pl.rc {
                                    (e < nq, s > 0)
                                } else {
                                    (s > 0, e < nq)
                                };
                                if before && rf > 0 {
                                    ends.push((
                                        (g0 + rf as u64 - 1) << 1,
                                        rtile(if pl.rc { e } else { s - 1 }, pl.rc),
                                        v,
                                    ));
                                }
                                if after && rf + m < l {
                                    ends.push((
                                        (g0 + (rf + m) as u64) << 1 | 1,
                                        rtile(if pl.rc { s - 1 } else { e }, pl.rc),
                                        v,
                                    ));
                                }
                                let fwd = fwd.get_or_init(|| query.tiles(s..e));
                                let mir = pl.rc.then(|| {
                                    mir.get_or_init(|| {
                                        fwd.iter().map(Tile::mirror).collect::<Vec<_>>()
                                    })
                                });
                                let units: &[u16] = if pl.rc {
                                    rev_units.get_or_init(|| {
                                        query.units[s..e]
                                            .iter()
                                            .rev()
                                            .map(|&u| Unit(u).rc().0)
                                            .collect::<Vec<_>>()
                                    })
                                } else {
                                    &query.units[s..e]
                                };
                                let rt = reference().tiles_with_units(rf..rf + m, units, false).unwrap();
                                for j in 0..m {
                                    let (qt, ri) = if pl.rc {
                                        (&mir.unwrap()[j], rf + m - 1 - j)
                                    } else {
                                        (&fwd[j], rf + j)
                                    };
                                    if let Some(p) = pack(qt, &rt[ri - rf]) {
                                        dif.push((g0 + ri as u64, p, v));
                                    }
                                }
                            }
                            link::observe_with(
                                |i| query.start(i),
                                |r| query.bases(r),
                                &prev,
                                prev_e,
                                &cin,
                                s,
                                &mut lk,
                            );
                            prev = cout;
                            prev_e = e;
                        }
                        for (v, arr) in [(&mut iv, &ev.depth), (&mut ivp, &ev.depth_p)] {
                            add_depth(v, arr);
                        }
                        for v in [&mut dif, &mut ends] {
                            v.sort_unstable_by(|a, b| {
                                (a.0, a.1, std::cmp::Reverse(a.2)).cmp(&(
                                    b.0,
                                    b.1,
                                    std::cmp::Reverse(b.2),
                                ))
                            });
                            v.dedup_by_key(|x| (x.0, x.1));
                        }
                        add_disagreement(&dif, &ev.disagree_p);
                        let add = |map: &Codes, g: u64, p: u64, w: u32| {
                            let mut v = map.entry(g).or_default();
                            match v.iter_mut().find(|x| x.0 == p) {
                                Some(x) => x.1 += w as u64,
                                None => v.push((p, w as u64)),
                            }
                        };
                        for &(g, p, w) in &dif {
                            add(&ev.diffs, g, p, w);
                        }
                        for &(g, p, w) in &ends {
                            add(&ev.brks, g, p, w);
                        }
                        if let Some(trace) = trace {trace.add(observations);}
                        link::add(&ev.links, lk);
                        (nsm, dif.len())
                    })
                    .map(|(m, d)| (1usize, m, d))
                    .reduce(|| (0, 0, 0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2))
            },
            |_, counts| {
                stats.0 += counts.0;
                stats.1 += counts.1;
                stats.2 += counts.2;
                eprint!(
                    "\r[tessellate hom]\t{path}: {} reads processed in {:.2}s",
                    stats.0,
                    t1.elapsed().as_secs_f64()
                );
                Ok(())
            },
        )
        .unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
        let (n, m, d) = stats;
        assert!(
            stats.0 == aux.offsets.len(),
            "{path}: {} records but its .tss has {}",
            stats.0,
            aux.offsets.len()
        );
        eprintln!(
            "[tessellate hom]\t{path}: {n} records, {m} SMEMs, {d} read-tile differences in {:.2}s",
            t1.elapsed().as_secs_f64()
        );
    }
    if reads {
        eprintln!(
            "[tessellate hom]\t{} occurrences located, {} SMEMs above --occ {occ} skipped (0 = locate all); {}; {} contig-end links",
            located.into_inner(),
            skipped.into_inner(),
            if ploidy == 0 {"--ploidy 0: every SMEM votes with full weight".to_owned()} else {format!("votes from SMEMs with at most --ploidy {ploidy} occurrences, weighted by 1/occurrences")},
            ev.links.len()
        );
    }
    ev
}
