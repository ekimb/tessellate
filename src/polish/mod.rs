use crate::hom::refs::Refs;
use crate::map::Record;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

#[derive(Clone, Copy, Debug)]
pub struct CopyNumber {
    pub min_len: usize,
    pub max_copies: usize,
    pub max_ratio: f64,
    pub min_extra: f64,
    pub min_copy_depth: f64,
    pub min_orphan: f64,
    pub min_other: f64,
    pub max_spacing: f64,
}

#[derive(Clone, Debug)]
pub struct Options {
    pub detect_agree: f64,
    pub detect_depth: f64,
    pub detect_hifi: Option<f64>,
    pub pad: usize,
    pub min_votes: usize,
    pub min_support: usize,
    pub max_depth: f64,
    pub tie_depth: Option<f64>,
    pub tie_margin: f64,
    pub flank_votes: bool,
    pub copy_number: Option<CopyNumber>,
}

pub struct Inputs<'a> {
    pub detect: &'a str,
    pub decide: &'a str,
    pub reads: &'a str,
    pub index: &'a str,
    pub fasta: &'a str,
    pub output: &'a str,
    pub edits: Option<&'a str>,
    pub windows: Option<&'a str>,
}

struct Window {
    rec: usize,
    start: i64,
    end: i64,
}

struct Obs {
    window: usize,
    ordinal: usize,
    qa: usize,
    qb: usize,
    rc: bool,
    m0: usize,
    m1: usize,
    flank: bool,
}

struct Copy {
    base: usize,
    depth: f64,
    left: [f64; 3],
    right: [f64; 3],
}

struct Repeat {
    rec: usize,
    len: usize,
    copies: Vec<Copy>,
}

pub struct Edit {
    pub rec: usize,
    pub start: usize,
    pub end: usize,
    pub old: Vec<u8>,
    pub new: Vec<u8>,
    pub votes: usize,
    pub support: usize,
    pub kind: &'static str,
}

#[derive(Default)]
pub struct Summary {
    pub candidates: usize,
    pub windows: usize,
    pub observations: usize,
    pub flank_observations: usize,
    pub edge_dropped: usize,
    pub anchor_dropped: usize,
    pub median_depth: f64,
    pub decisions: HashMap<&'static str, usize>,
    pub edits: usize,
    pub copy_number: usize,
}

fn range(s: &str) -> (usize, usize) {
    let (a, b) = s.split_once("..").unwrap_or((s, s));
    (a.parse().unwrap_or(0), b.parse().unwrap_or(0))
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| match b {b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x}).collect()
}

fn links(s: &str) -> [f64; 3] {
    let mut out = [0.0; 3];
    if s != "-" {
        for (i, x) in s.split('/').take(3).enumerate() {out[i] = x.parse().unwrap_or(0.0);}
    }
    out
}

fn open(path: &str) -> io::Result<BufReader<File>> {
    Ok(BufReader::with_capacity(1 << 22, File::open(path)?))
}

fn detect(path: &str, names: &HashMap<String, usize>, agree: f64, depth: f64) -> io::Result<Vec<(usize, usize)>> {
    let mut out = Vec::new();
    let mut line = String::new();
    let mut reader = open(path)?;
    while reader.read_line(&mut line)? > 0 {
        if line.starts_with('X') {break;}
        if line.starts_with('T') {
            let f: Vec<&str> = line.splitn(7, '\t').collect();
            if let (Some(&rec), Ok(tile), Ok(d), Ok(a)) = (names.get(f[1]), f[2].parse::<usize>(), f[4].parse::<f64>(), f[5].parse::<f64>()) {
                if d >= depth && a / d < agree {out.push((rec, tile));}
            }
        }
        line.clear();
    }
    out.sort_unstable();
    Ok(out)
}

fn windows(candidates: &[(usize, usize)], pad: usize) -> Vec<Window> {
    let mut out: Vec<Window> = Vec::new();
    let pad = pad as i64;
    for &(rec, tile) in candidates {
        let t = tile as i64;
        match out.last_mut() {
            Some(w) if w.rec == rec && t - pad <= w.end + 1 => w.end = t + pad,
            _ => out.push(Window {rec, start: t - pad, end: t + pad}),
        }
    }
    out
}

type BreakKey = (usize, usize, usize, bool, bool);

struct Scan {
    obs: Vec<Obs>,
    breaks: HashSet<BreakKey>,
    depths: Vec<f32>,
    repeats: Vec<Repeat>,
}

fn scan(path: &str, names: &HashMap<String, usize>, wins: &[Window], cn: Option<CopyNumber>, flank: bool) -> io::Result<Scan> {
    let mut by_rec: HashMap<usize, Vec<usize>> = HashMap::new();
    for (k, w) in wins.iter().enumerate() {by_rec.entry(w.rec).or_default().push(k);}
    let starts: HashMap<usize, Vec<i64>> = by_rec.iter().map(|(&r, ks)| (r, ks.iter().map(|&k| wins[k].start).collect())).collect();
    let mut ends: HashSet<(usize, i64)> = HashSet::new();
    for w in wins {ends.insert((w.rec, w.start - 1)); ends.insert((w.rec, w.end + 1));}
    let mut region: HashMap<usize, usize> = HashMap::new();
    let (mut obs, mut breaks, mut depths, mut repeats) = (Vec::new(), HashSet::new(), Vec::new(), Vec::new());
    let mut left: HashMap<(usize, usize, bool), Vec<usize>> = HashMap::new();
    let mut right: HashMap<(usize, usize, bool), Vec<usize>> = HashMap::new();
    let mut line = String::new();
    let mut reader = open(path)?;
    while reader.read_line(&mut line)? > 0 {
        let first = line.as_bytes().first().copied();
        match first {
            Some(b'T') => {
                let f: Vec<&str> = line.splitn(6, '\t').collect();
                if let Ok(d) = f[4].parse::<f32>() {depths.push(d);}
            }
            Some(b'G') => {
                let f: Vec<&str> = line.split('\t').collect();
                if let (Ok(r), Some(&rec)) = (f[1].parse::<usize>(), names.get(f[2])) {region.insert(r, rec);}
            }
            Some(b'R') if cn.is_some() => {
                let f: Vec<&str> = line.trim_end().split('\t').collect();
                let (Some(&rec), Ok(b0), Ok(b1)) = (names.get(f[1]), f[4].parse::<usize>(), f[5].parse::<usize>()) else {line.clear(); continue};
                let cn = cn.unwrap();
                if b1 - b0 >= cn.min_len {
                    let copies: Vec<Copy> = f[7].split(';').filter_map(|c| {
                        let p: Vec<&str> = c.split(':').collect();
                        if p.len() < 9 || p[0] != f[1] {return None;}
                        Some(Copy {base: p[2].parse().ok()?, depth: p[5].parse().ok()?, left: links(p[7]), right: links(p[8])})
                    }).collect();
                    if copies.len() >= 2 && copies.len() <= cn.max_copies {repeats.push(Repeat {rec, len: b1 - b0, copies});}
                }
            }
            Some(b'X') => {
                let f: Vec<&str> = line.trim_end().split('\t').collect();
                if f.len() < 13 {line.clear(); continue;}
                let Some(&rec) = f[1].parse::<usize>().ok().and_then(|r| region.get(&r)) else {line.clear(); continue};
                let (r0, r1) = range(f[10]);
                let (m0, m1) = range(f[12]);
                let ordinal: usize = f[3].parse().unwrap_or(0);
                let rc = f[7] == "-";
                if f[6] == "M" {
                    if let Some(st) = starts.get(&rec) {
                        let ks = &by_rec[&rec];
                        let mut j = st.partition_point(|&s| s < r0 as i64);
                        let (q0, _) = range(f[8]);
                        let at = |x: usize| if rc {q0 + (r1 - 1 - x)} else {q0 + (x - r0)};
                        if flank && j > 0 {
                            let k = ks[j - 1];
                            let (a, b) = (wins[k].start, wins[k].end);
                            if a < r0 as i64 && b >= r0 as i64 && b < r1 as i64 {right.entry((k, ordinal, rc)).or_default().push(at(b as usize));}
                        }
                        while j < ks.len() && st[j] < r1 as i64 {
                            let k = ks[j];
                            let (a, b) = (wins[k].start, wins[k].end);
                            if b < r1 as i64 {
                                let (a, b) = (a as usize, b as usize);
                                let (qa, qb) = if rc {(at(b), at(a))} else {(at(a), at(b))};
                                obs.push(Obs {window: k, ordinal, qa, qb, rc, m0, m1, flank: false});
                            } else if flank {
                                left.entry((k, ordinal, rc)).or_default().push(at(a as usize));
                            }
                            j += 1;
                        }
                    }
                } else if ends.contains(&(rec, r0 as i64)) {
                    breaks.insert((ordinal, m0, m1, f[6] == "B<", rc));
                }
            }
            _ => {}
        }
        line.clear();
    }
    for ((k, ordinal, rc), ls) in &left {
        let Some(rs) = right.get(&(*k, *ordinal, *rc)) else {continue};
        let span = (wins[*k].end - wins[*k].start) as i64;
        let slack = (span / 4).max(4);
        let mut best: Option<(i64, usize, usize)> = None;
        for &l in ls {
            for &r in rs {
                let (qa, qb) = if *rc {(r, l)} else {(l, r)};
                if qb <= qa {continue;}
                let d = ((qb - qa) as i64 - span).abs();
                if d > slack {continue;}
                if best.is_none_or(|x| d < x.0) {best = Some((d, qa, qb));}
            }
        }
        if let Some((_, qa, qb)) = best {obs.push(Obs {window: *k, ordinal: *ordinal, qa, qb, rc: *rc, m0: 0, m1: 0, flank: true});}
    }
    Ok(Scan {obs, breaks, depths, repeats})
}

fn median(v: &mut [f32]) -> f64 {
    if v.is_empty() {return 0.0;}
    let n = v.len();
    let mid = n / 2;
    v.select_nth_unstable_by(mid, |a, b| a.total_cmp(b));
    if n % 2 == 1 {v[mid] as f64} else {
        let hi = v[mid];
        let lo = v[..mid].iter().copied().fold(f32::MIN, f32::max);
        (lo as f64 + hi as f64) / 2.0
    }
}

fn read_fasta(path: &str) -> io::Result<(Vec<String>, HashMap<String, Vec<u8>>)> {
    let (mut names, mut seqs, mut name, mut buf): (Vec<String>, HashMap<String, Vec<u8>>, Option<String>, Vec<u8>) = (Vec::new(), HashMap::new(), None, Vec::new());
    for line in open(path)?.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            if let Some(n) = name.take() {seqs.insert(n, std::mem::take(&mut buf));}
            let n = h.split_whitespace().next().unwrap_or_default().to_owned();
            names.push(n.clone());
            name = Some(n);
        } else {
            buf.extend(line.trim_end().bytes());
        }
    }
    if let Some(n) = name {seqs.insert(n, buf);}
    Ok((names, seqs))
}

fn copy_number(repeats: &[Repeat], median: f64, cn: CopyNumber) -> Vec<(usize, usize, usize)> {
    struct Call {rec: usize, start: usize, end: usize, len: usize, ratio: f64, cut: (usize, usize)}
    let mut calls: Vec<Call> = Vec::new();
    for r in repeats {
        let mut copies: Vec<&Copy> = r.copies.iter().collect();
        copies.sort_by_key(|c| c.base);
        let total: f64 = copies.iter().map(|c| c.depth).sum();
        let n = copies.len() as f64;
        let ratio = total / (n * median);
        let implied = total / median;
        let missing = n - implied;
        if ratio > cn.max_ratio || missing < cn.min_extra {continue;}
        if copies.iter().map(|c| c.depth).fold(f64::INFINITY, f64::min) < cn.min_copy_depth * median {continue;}
        let orphan = copies.iter().map(|c| c.left[2].max(c.right[2])).fold(0.0, f64::max);
        let other = copies.iter().map(|c| c.left[1].max(c.right[1])).fold(0.0, f64::max);
        if !(orphan >= cn.min_orphan || other >= cn.min_other) {continue;}
        let spacing = copies.windows(2).map(|w| w[1].base - w[0].base).min().unwrap_or(r.len);
        if spacing as f64 > cn.max_spacing * r.len as f64 {continue;}
        let extra = (missing.round() as usize).max(1);
        let start = copies[0].base;
        calls.push(Call {rec: r.rec, start, end: copies.last().unwrap().base + r.len, len: r.len, ratio, cut: (start, start + extra * spacing)});
    }
    calls.sort_by(|a, b| (a.rec, a.start).cmp(&(b.rec, b.start)).then(b.len.cmp(&a.len)));
    let mut events: Vec<(usize, usize, usize, f64, (usize, usize))> = Vec::new();
    for c in calls {
        match events.last_mut() {
            Some(e) if e.0 == c.rec && c.start < e.2 => {
                e.2 = e.2.max(c.end);
                if (c.len, -c.ratio) > (e.1, -e.3) {e.1 = c.len; e.3 = c.ratio; e.4 = c.cut;}
            }
            _ => events.push((c.rec, c.len, c.end, c.ratio, c.cut)),
        }
    }
    events.into_iter().map(|e| (e.0, e.4.0, e.4.1)).collect()
}

pub fn polish(inputs: &Inputs<'_>, opt: &Options) -> io::Result<Summary> {
    let started = Instant::now();
    let mut sum = Summary::default();
    let prefix = |p: &str| crate::tess::tq::stem(p).unwrap_or(p).to_owned();
    let asm = Refs::open_prefix(&prefix(inputs.index))?;
    let records: Vec<Record> = (0..asm.n_records()).into_par_iter().map(|r| Record::read(&asm, r)).collect::<io::Result<_>>()?;
    let names: HashMap<String, usize> = records.iter().enumerate().map(|(i, r)| (String::from_utf8_lossy(&r.name).into_owned(), i)).collect();
    eprintln!("[tessellate polish]\t{}: {} records decoded in {:.2}s", inputs.index, records.len(), started.elapsed().as_secs_f64());
    let t = Instant::now();
    let mut candidates = detect(inputs.detect, &names, opt.detect_agree, opt.detect_depth)?;
    if let Some(agree) = opt.detect_hifi {
        let n = candidates.len();
        candidates.extend(detect(inputs.decide, &names, agree, opt.detect_depth)?);
        candidates.sort_unstable();
        candidates.dedup();
        eprintln!("[tessellate polish]\t{}: {} more candidate tiles from --detect-hifi {agree}", inputs.decide, candidates.len() - n);
    }
    let wins = windows(&candidates, opt.pad);
    sum.candidates = candidates.len();
    sum.windows = wins.len();
    eprintln!("[tessellate polish]\t{}: {} candidate tiles in {} windows in {:.2}s", inputs.detect, sum.candidates, sum.windows, t.elapsed().as_secs_f64());
    let t = Instant::now();
    let Scan {obs, breaks, mut depths, repeats} = scan(inputs.decide, &names, &wins, opt.copy_number, opt.flank_votes)?;
    sum.median_depth = median(&mut depths);
    drop(depths);
    sum.observations = obs.len();
    sum.flank_observations = obs.iter().filter(|o| o.flank).count();
    eprintln!("[tessellate polish]\t{}: {} observations ({} from flanking SMEMs), {} break rows, {} repeat rows, median depth {:.1} in {:.2}s", inputs.decide, obs.len(), sum.flank_observations, breaks.len(), repeats.len(), sum.median_depth, t.elapsed().as_secs_f64());
    let t = Instant::now();
    let reads = Refs::open_prefix(&prefix(inputs.reads))?;
    let mut by_read: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, o) in obs.iter().enumerate() {by_read.entry(o.ordinal).or_default().push(i);}
    let groups: Vec<(usize, Vec<usize>)> = by_read.into_iter().collect();
    let in_contig = |w: &Window| w.start >= 0 && (w.end as usize) < records[w.rec].n_tiles();
    let strings: Vec<(usize, Option<Vec<u8>>, bool)> = groups.par_iter().map(|(ordinal, idx)| -> io::Result<Vec<(usize, Option<Vec<u8>>, bool)>> {
        let lo = idx.iter().map(|&i| obs[i].qa).min().unwrap();
        let hi = idx.iter().map(|&i| obs[i].qb).max().unwrap();
        let view = reads.view(*ordinal)?;
        let n = view.n_tiles();
        let tiles = if hi < n {view.tiles(lo..hi + 1, false)?} else {Vec::new()};
        Ok(idx.iter().map(|&i| {
            let o = &obs[i];
            let w = &wins[o.window];
            let (a, b) = (w.start, w.end);
            let edge = !o.flank && ((a == o.m0 as i64 && breaks.contains(&(o.ordinal, o.m0, o.m1, true, o.rc))) || (b == o.m1 as i64 - 1 && breaks.contains(&(o.ordinal, o.m0, o.m1, false, o.rc))));
            if edge {return (i, None, true);}
            if tiles.is_empty() || !in_contig(w) {return (i, None, false);}
            let rec = &records[w.rec];
            let span = &tiles[o.qa - lo..=o.qb - lo];
            let fwd: Vec<u8> = span.iter().flat_map(|t| t.render()).collect();
            let (first, last) = if o.rc {(revcomp(&span[span.len() - 1].render()), revcomp(&span[0].render()))} else {(span[0].render(), span[span.len() - 1].render())};
            let (a, b) = (a as usize, b as usize);
            if first != rec.bases(a..=a) || last != rec.bases(b..=b) {return (i, None, false);}
            (i, Some(if o.rc {revcomp(&fwd)} else {fwd}), false)
        }).collect())
    }).collect::<io::Result<Vec<_>>>()?.into_iter().flatten().collect();
    let mut vote_of: Vec<Option<Vec<u8>>> = vec![None; obs.len()];
    let mut per_window: Vec<[usize; 3]> = vec![[0; 3]; wins.len()];
    for (i, s, edge) in strings {
        let k = obs[i].window;
        per_window[k][0] += 1;
        if edge {sum.edge_dropped += 1; per_window[k][1] += 1;} else if s.is_none() {sum.anchor_dropped += 1; per_window[k][2] += 1;}
        vote_of[i] = s;
    }
    eprintln!("[tessellate polish]\t{}: read tiles rendered for {} reads; {} votes dropped at an SMEM edge next to a break, {} for anchor mismatch, in {:.2}s", inputs.reads, groups.len(), sum.edge_dropped, sum.anchor_dropped, t.elapsed().as_secs_f64());
    let t = Instant::now();
    let mut votes: Vec<Vec<(usize, usize, bool, usize)>> = vec![Vec::new(); wins.len()];
    for (i, o) in obs.iter().enumerate() {if vote_of[i].is_some() {votes[o.window].push((o.ordinal, o.qa, o.rc, i));}}
    let cap = opt.max_depth * sum.median_depth;
    let mut edits: Vec<Edit> = Vec::new();
    let mut outcome: Vec<(&'static str, usize, usize, usize, usize, String)> = vec![("off contig end", 0, 0, 0, 0, String::new()); wins.len()];
    for (k, w) in wins.iter().enumerate() {
        let mut tally = |s: &'static str| *sum.decisions.entry(s).or_default() += 1;
        if !in_contig(w) {tally("off contig end"); continue;}
        let mut seen: HashSet<(usize, usize, bool)> = HashSet::new();
        let mut order: Vec<&[u8]> = Vec::new();
        let mut count: HashMap<&[u8], usize> = HashMap::new();
        for &(ordinal, qa, rc, i) in &votes[k] {
            if !seen.insert((ordinal, qa, rc)) {continue;}
            let s = vote_of[i].as_deref().unwrap();
            let c = count.entry(s).or_insert(0);
            if *c == 0 {order.push(s);}
            *c += 1;
        }
        let n: usize = count.values().sum();
        let mut lens: Vec<(usize, usize)> = Vec::new();
        for (s, c) in &count {lens.push((s.len(), *c));}
        lens.sort_unstable();
        let mut hist: Vec<(usize, usize)> = Vec::new();
        for &(l, c) in &lens {match hist.last_mut() {Some(h) if h.0 == l => h.1 += c, _ => hist.push((l, c))}}
        let hist = hist.iter().map(|(l, c)| format!("{l}:{c}")).collect::<Vec<_>>().join(",");
        let (a0, b0) = (w.start as usize, w.end as usize);
        let clen = records[w.rec].bases(a0..=b0).len();
        outcome[k] = ("votes below minimum", n, 0, 0, clen, hist);
        if n < opt.min_votes {tally("votes below minimum"); continue;}
        outcome[k].0 = "depth above cap";
        if n as f64 > cap {tally("depth above cap"); continue;}
        let (mut acc, mut med) = (0, 0);
        for (l, c) in lens {acc += c; if 2 * acc >= n {med = l; break;}}
        let mut best: Option<(&[u8], usize)> = None;
        for s in &order {
            if s.len() != med {continue;}
            let c = count[s];
            if match best {None => true, Some((_, bc)) => c > bc} {best = Some((s, c));}
        }
        let (s, c) = best.unwrap();
        let rec = &records[w.rec];
        let (a, b) = (w.start as usize, w.end as usize);
        let contig = rec.bases(a..=b);
        let at_len = |l: usize| count.iter().filter(|(x, _)| x.len() == l).map(|(_, &v)| v).sum::<usize>();
        outcome[k].2 = med;
        outcome[k].3 = c;
        outcome[k].0 = "reads agree with contig";
        if s == contig {tally("reads agree with contig"); continue;}
        outcome[k].0 = "chosen string below support";
        if c < opt.min_support {tally("chosen string below support"); continue;}
        if let Some(f) = opt.tie_depth {
            outcome[k].0 = "high-depth tie";
            if n as f64 > f * sum.median_depth && (at_len(med) as f64) < opt.tie_margin * at_len(contig.len()) as f64 {tally("high-depth tie"); continue;}
        }
        outcome[k].0 = "edit";
        tally("edit");
        let start = rec.tile_start(a) as usize;
        edits.push(Edit {rec: w.rec, start, end: start + contig.len(), old: contig.to_vec(), new: s.to_vec(), votes: n, support: c, kind: "window"});
    }
    sum.edits = edits.len();
    eprintln!("[tessellate polish]\t{} window edits in {:.2}s", edits.len(), t.elapsed().as_secs_f64());
    if let Some(path) = inputs.windows {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(w, "#contig\tstart\tend\toutcome\tobservations\tedge_dropped\tanchor_dropped\tvotes\tmedian_len\tsupport\tcontig_len\tvote_lengths")?;
        for (k, win) in wins.iter().enumerate() {
            let rec = &records[win.rec];
            let (s, e) = if in_contig(win) {let s = rec.tile_start(win.start as usize) as usize; (s as i64, (s + rec.bases(win.start as usize..=win.end as usize).len()) as i64)} else {(-1, -1)};
            let (o, n, med, c, clen, hist) = &outcome[k];
            let [nobs, edge, anchor] = per_window[k];
            writeln!(w, "{}\t{s}\t{e}\t{o}\t{nobs}\t{edge}\t{anchor}\t{n}\t{med}\t{c}\t{clen}\t{}", String::from_utf8_lossy(&rec.name), if hist.is_empty() {"-"} else {hist})?;
        }
        w.flush()?;
    }
    if let Some(cn) = opt.copy_number {
        for (rec, s, e) in copy_number(&repeats, sum.median_depth, cn) {
            eprintln!("[tessellate polish]\tcopy number: delete {}:{}-{} ({} bp)", String::from_utf8_lossy(&records[rec].name), s, e, e - s);
            edits.retain(|x| x.rec != rec || x.end <= s || x.start >= e);
            edits.push(Edit {rec, start: s, end: e, old: Vec::new(), new: Vec::new(), votes: 0, support: 0, kind: "copy_number"});
            sum.copy_number += 1;
        }
    }
    let (order, mut seqs) = read_fasta(inputs.fasta)?;
    let mut per: HashMap<usize, Vec<&Edit>> = HashMap::new();
    for e in &edits {per.entry(e.rec).or_default().push(e);}
    let mut out = BufWriter::new(File::create(inputs.output)?);
    for name in &order {
        let mut seq = seqs.remove(name).unwrap_or_default();
        if let Some(&rec) = names.get(name) {
            let mut list = per.remove(&rec).unwrap_or_default();
            list.sort_by_key(|e| std::cmp::Reverse(e.start));
            let mut limit = seq.len();
            for e in list {
                if e.end > limit {continue;}
                if e.kind == "window" && seq[e.start..e.end] != e.old[..] {
                    return Err(io::Error::other(format!("{name}: FASTA differs from the index at {}", e.start)));
                }
                seq.splice(e.start..e.end, e.new.iter().copied());
                limit = e.start;
            }
        }
        writeln!(out, ">{name}")?;
        for chunk in seq.chunks(80) {out.write_all(chunk)?; out.write_all(b"\n")?;}
    }
    out.flush()?;
    if let Some(path) = inputs.edits {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(w, "#contig\tstart\tend\tkind\told_len\tnew\tvotes\tsupport")?;
        let mut sorted: Vec<&Edit> = edits.iter().collect();
        sorted.sort_by_key(|e| (e.rec, e.start));
        for e in sorted {
            writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", String::from_utf8_lossy(&records[e.rec].name), e.start, e.end, e.kind, e.end - e.start, String::from_utf8_lossy(&e.new), e.votes, e.support)?;
        }
        w.flush()?;
    }
    eprintln!("[tessellate polish]\t{}: {} window edits, {} copy-number deletions written in {:.2}s total", inputs.output, sum.edits, sum.copy_number, started.elapsed().as_secs_f64());
    Ok(sum)
}
