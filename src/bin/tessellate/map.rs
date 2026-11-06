use super::hom::{RepeatArgs, open_table};
use std::io::{self, BufWriter, Write};
use std::time::Instant;
use tessellate::hom::{Ctx, chain};
use tessellate::map::myers;

#[derive(clap::Args)]
pub struct Options {
    #[arg(short, long)]
    pub input: String,
    #[arg(short, long)]
    pub output: String,
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
    #[arg(long, default_value_t = 5)]
    pub top: usize,
    #[arg(long, default_value_t = 1000)]
    pub gap: usize,
    #[arg(long, default_value_t = 64)]
    pub band: usize,
    #[arg(short = 'g', long)]
    pub gt: Option<usize>,
    #[arg(long)]
    pub vote: bool,
    #[arg(long, help = "short reads (WIP)")]
    pub sr: bool,
}

pub fn run(opt: Options, repeats: RepeatArgs, reads: &[String]) -> io::Result<()> {
    if repeats.k == Some(0) {return Err(io::Error::other("-k/--k must be positive"));}
    if reads.is_empty() {return Err(io::Error::other("map requires read .tq.gz input"));}
    if opt.output == "-" {return Err(io::Error::other("map requires an output file"));}
    let options = chain::Options {top: opt.top, gap: opt.gap, band: opt.band, min_gap: opt.gt,
        vote: opt.vote.then_some(myers::Options {window: 7, pad: 3}), sr: opt.sr};
    let cmd = std::env::args().collect::<Vec<_>>().join(" ");
    let started = Instant::now();
    let (index, path, fm, refs, table) = open_table(&opt.input, repeats);
    let ctx = Ctx {fm: &fm, refs: &refs};
    eprintln!(
        "[tessellate map]\t{index}: {} records, {} tiles; {path}: {} repeats, {} sites, {} contained contigs, {} families loaded in {:.2}s",
        fm.n_refs(),
        fm.n_tiles(),
        table.n(),
        table.n_sites(),
        table.n_contained(),
        table.fams.len(),
        started.elapsed().as_secs_f64()
    );
    let (k, occ) = (repeats.k.unwrap_or(table.k), repeats.occ.unwrap_or(table.occ));
    eprintln!("[tessellate map]\tread seeding: minimum {k} tiles, maximum occurrences {occ} (0 = unlimited); k from {}, occ from {}",
        if repeats.k.is_some() {"command line"} else {".trep"},
        if repeats.occ.is_some() {"command line"} else {".trep"});
    let chains = chain::Index::new(&ctx, &table, options);
    let out = &opt.output;
    let mut writer = BufWriter::new(std::fs::File::create(out).map_err(|e| io::Error::other(format!("cannot create {out}: {e}")))?);
    writeln!(writer, "#cmd\t{cmd}")?;
    writer.write_all(chain::HEADER)?;
    chains.scan(&ctx, reads, k, occ, &mut writer)?;
    writer.flush()?;
    eprintln!("[tessellate map]\tread candidate regions written to {out}");
    Ok(())
}
