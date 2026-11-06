use super::{Ctx, FIX};
use dashmap::DashMap;
use std::{collections::HashMap, io::{self, BufRead, Write}, ops::Range, path::Path, sync::Arc};

#[derive(Clone)]
struct Region { rec: usize, start: usize, end: usize, label: String, tiles: Range<usize> }

#[derive(Clone, Eq, PartialEq)]
pub struct Observation {
    pub region: usize,
    pub source: Arc<str>,
    pub ordinal: usize,
    pub offset: u64,
    pub read: Arc<str>,
    pub kind: &'static str,
    pub rc: bool,
    pub query: Range<usize>,
    pub query_bases: Range<usize>,
    pub reference: Range<usize>,
    pub smem_query: Range<usize>,
    pub smem_reference: Range<usize>,
    pub occurrences: usize,
    pub weight: u32,
}

pub struct Trace { regions: Vec<Region>, by_record: Vec<Vec<usize>>, max_end: Vec<Vec<usize>>, rows: DashMap<usize, Vec<Observation>> }
impl Trace {
    pub fn whole(ctx: &Ctx) -> Self {
        let regions: Vec<_> = (0..ctx.fm.n_refs()).map(|rec| Region {rec, start: 0, end: ctx.len(rec),
            label: "whole_contig".to_owned(), tiles: 0..ctx.ntiles(rec)}).collect();
        let by_record = (0..regions.len()).map(|i| vec![i]).collect();
        let max_end = regions.iter().map(|r| vec![r.tiles.end]).collect();
        Self {regions, by_record, max_end, rows: DashMap::new()}
    }
    pub fn open(path: &str, ctx: &Ctx) -> io::Result<Self> {
        let mut regions = Vec::new();
        let mut by_record = vec![Vec::new(); ctx.fm.n_refs()];
        let file = std::fs::File::open(Path::new(path))?;
        for (line, row) in io::BufReader::new(file).lines().enumerate() {
            let row = row?;
            if row.is_empty() || row.starts_with('#') || row.starts_with("track ") || row.starts_with("browser ") {continue;}
            let f: Vec<_> = row.split('\t').collect();
            let invalid = || io::Error::other(format!("{path}: invalid BED line {}", line + 1));
            if f.len() < 3 {return Err(invalid());}
            let rec = (0..ctx.fm.n_refs()).find(|&i| ctx.get(i).name() == f[0].as_bytes()).ok_or_else(invalid)?;
            let view = ctx.get(rec);
            let start: usize = f[1].parse().map_err(|_| invalid())?;
            let end: usize = f[2].parse().map_err(|_| invalid())?;
            if start >= end || end > view.len() {return Err(invalid());}
            let n = ctx.ntiles(rec);
            let lower = |pos: usize| -> io::Result<usize> {
                let (mut lo, mut hi) = (0, n + 1);
                while lo < hi {let mid = lo + (hi-lo)/2; if view.start(mid)? < pos {lo = mid+1} else {hi = mid}}
                Ok(lo)
            };
            let first = lower(start.saturating_add(1))?.saturating_sub(1);
            let last = lower(end)?;
            let label = f.get(3).filter(|s| !s.is_empty()).map_or_else(|| format!("region{}", regions.len()+1), |s| (*s).to_owned());
            by_record[rec].push(regions.len());
            regions.push(Region {rec, start, end, label, tiles: first..last});
        }
        if regions.is_empty() {return Err(io::Error::other(format!("{path}: no BED regions")));}
        let mut max_end = Vec::new();
        for ids in &mut by_record {
            ids.sort_unstable_by_key(|&i| (regions[i].tiles.start, regions[i].tiles.end, i));
            let mut end = 0;
            max_end.push(ids.iter().map(|&i| {end = end.max(regions[i].tiles.end); end}).collect());
        }
        Ok(Self {regions, by_record, max_end, rows: DashMap::new()})
    }
    pub fn overlapping(&self, rec: usize, range: Range<usize>) -> impl Iterator<Item=(usize, Range<usize>)> + '_ {
        let ids = &self.by_record[rec];
        let end = ids.partition_point(|&i| self.regions[i].tiles.start < range.end);
        let start = self.max_end[rec][..end].partition_point(|&e| e <= range.start);
        ids[start..end].iter().filter_map(move |&i| {
            let t = &self.regions[i].tiles;
            let s = t.start.max(range.start); let e = t.end.min(range.end);
            (s < e).then_some((i, s..e))
        })
    }
    pub fn add(&self, rows: Vec<Observation>) {
        let mut grouped: HashMap<usize, Vec<Observation>> = HashMap::new();
        for row in rows {grouped.entry(row.region).or_default().push(row);}
        for (region, mut rows) in grouped {self.rows.entry(region).or_default().append(&mut rows);}
    }
    pub fn write(self, out: &mut impl Write, ctx: &Ctx) -> io::Result<usize> {
        writeln!(out, "#G\tregion\tcontig\tbase_start\tbase_end\tlabel\ttile_start\ttile_end")?;
        writeln!(out, "#X\tregion\tread_source\tread_ordinal\ttss_offset\tread\tevidence\tstrand\tquery_tiles\tquery_bases\treference_tiles\tsmem_query_tiles\tsmem_reference_tiles\toccurrences\tvote_weight")?;
        let mut count = 0;
        for (i, region) in self.regions.iter().enumerate() {
            writeln!(out, "G\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i+1, escape(&ctx.name(region.rec)), region.start, region.end, escape(&region.label), region.tiles.start, region.tiles.end)?;
            if let Some((_, mut rows)) = self.rows.remove(&i) {
                rows.sort_unstable_by(|a,b| ((a.source.as_ref(),a.ordinal,a.offset,a.read.as_ref(),a.kind,a.rc), (a.query.start,a.query.end,a.reference.start,a.reference.end,a.smem_query.start,a.smem_query.end,a.smem_reference.start,a.smem_reference.end,a.occurrences,a.weight)).cmp(&((b.source.as_ref(),b.ordinal,b.offset,b.read.as_ref(),b.kind,b.rc), (b.query.start,b.query.end,b.reference.start,b.reference.end,b.smem_query.start,b.smem_query.end,b.smem_reference.start,b.smem_reference.end,b.occurrences,b.weight))));
                rows.dedup();
                for x in rows {
                    writeln!(out, "X\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}..{}\t{}..{}\t{}..{}\t{}..{}\t{}..{}\t{}\t{:.3}",
                        i+1, escape(&x.source), x.ordinal, x.offset, escape(&x.read), x.kind, if x.rc {'-'} else {'+'},
                        x.query.start,x.query.end,x.query_bases.start,x.query_bases.end,x.reference.start,x.reference.end,
                        x.smem_query.start,x.smem_query.end,x.smem_reference.start,x.smem_reference.end,x.occurrences,x.weight as f64/FIX as f64)?;
                    count += 1;
                }
            }
        }
        Ok(count)
    }
}
fn escape(s: &str) -> String {
    s.chars().flat_map(|c| match c {'%' => "%25".chars().collect::<Vec<_>>(), '\t' => "%09".chars().collect(), '\n' => "%0A".chars().collect(), '\r' => "%0D".chars().collect(), _ => vec![c]}).collect()
}
