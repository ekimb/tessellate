use std::io::{self, BufWriter, Write};
use std::time::Instant;
use tessellate::Fmd;
use tessellate::hom::{Ctx, contain, family, link, pair, put, read, refs::Refs, repeat, tile};
use tessellate::hom::repeat::RepeatOptions;

pub fn tq_gz_path(path: &str) -> Result<String, String> {
    if path.to_ascii_lowercase().ends_with(".tq.gz") {Ok(path.to_owned())}
    else {Err("reads must end in .tq.gz and have a matching .tss; prepare raw reads with tessellate fastx".to_owned())}
}

#[derive(Clone, Copy, Debug, clap::Args)]
pub struct RepeatArgs {
    #[arg(short = 'k', long)]
    pub k: Option<usize>,
    #[arg(short = 'O', long)]
    pub occ: Option<usize>,
    #[arg(long, default_value_t = 0.05)]
    pub frac: f64,
}
impl RepeatArgs {
    fn discovery(self) -> RepeatOptions {
        RepeatOptions {k: self.k.unwrap_or(31), occ: self.occ.unwrap_or(0), frac: self.frac}
    }
}

#[derive(clap::Args)]
pub struct Options {
    #[arg(short, long)]
    pub input: String,
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
    #[arg(short, long)]
    pub output: Option<String>,
    #[arg(short = 'p', long, default_value_t = 2)]
    pub ploidy: usize,
    #[arg(short = 'g', long, default_value_t = 100)]
    pub gt: usize,
    #[arg(long, default_value_t = 1.0)]
    pub min: f64,
    #[arg(long, value_name = "BED", conflicts_with = "trace")]
    pub trace_bed: Option<String>,
    #[arg(long)]
    pub trace: bool,
}

pub fn hom(opt: Options, reads: &[String], repeats: RepeatArgs) -> io::Result<()> {
    if (opt.trace || opt.trace_bed.is_some()) && reads.is_empty() {return Err(io::Error::other("--trace/--trace-bed requires read .tq.gz input"));}
    if !opt.min.is_finite() || opt.min < 0.0 {
        return Err(io::Error::other("--min must be finite and nonnegative"));
    }
    let command = std::env::args().collect::<Vec<_>>().join(" ");
    report(&opt.input, reads, opt.ploidy, opt.gt, opt.min,
        opt.output.as_deref().unwrap_or("-"), &command, repeats, opt.trace, opt.trace_bed.as_deref());
    Ok(())
}

pub fn open_table(index: &str, repeats: RepeatArgs) -> (String, String, Fmd, Refs, repeat::Table) {
    let prefix = index.strip_suffix(".tbwt").unwrap_or(index);
    let index = format!("{prefix}.tbwt");
    let path = format!("{prefix}.trep");
    let fm = Fmd::read(&std::fs::read(&index).unwrap_or_else(|e| panic!("cannot read {index}: {e}")));
    let refs = Refs::open(&fm, &index);
    let table = match std::fs::metadata(&path) {
        Ok(_) => repeat::Table::read(&path),
        Err(error) if error.kind() == io::ErrorKind::NotFound => {
            let repeats = repeats.discovery();
            let views: Vec<_> = (0..fm.n_refs()).map(|rec| refs.get(rec)).collect();
            eprintln!("[tessellate hom]\t{path} missing; building repeat table from {index}: minimum {} tiles, maximum occurrences {} (0 = unlimited), fraction {}",
                repeats.k, repeats.occ, repeats.frac);
            let table = repeat::Table::build(&fm, &views, repeats);
            table.save(&path).unwrap_or_else(|e| panic!("cannot save {path}: {e}"));
            eprintln!("[tessellate hom]\t{path}: {} repeats, {} sites written", table.n(), table.n_sites());
            table
        }
        Err(error) => panic!("cannot inspect {path}: {error}"),
    };
    let index_modified = std::fs::metadata(&index).and_then(|m| m.modified()).unwrap();
    let table_modified = std::fs::metadata(&path).and_then(|m| m.modified()).unwrap();
    assert!(table.matches_index(&fm) && table_modified >= index_modified, "{path}: repeat table is obsolete or does not match {index}; remove it and rerun tessellate hom -i {index}");
    (index, path, fm, refs, table)
}

fn report(index: &str, reads: &[String], ploidy: usize, gt: usize, min: f64, out: &str, cmd: &str, repeats: RepeatArgs, trace_all: bool, trace_bed: Option<&str>) {
    assert!(min.is_finite() && min >= 0.0, "--min must be finite and nonnegative");
    let started = Instant::now();
    let (index, path, fm, refs, table) = open_table(index, repeats);
    let ctx = Ctx {
        fm: &fm,
        refs: &refs,
    };
    let nr = fm.n_refs();
    let (k, occ) = (table.k, table.occ);
    eprintln!(
        "[tessellate hom]\t{index}: {nr} records, {} tiles; {path}: {} repeats, {} sites, {} contained contigs, {} families loaded in {:.2}s",
        fm.n_tiles(),
        table.n(),
        table.n_sites(),
        table.n_contained(),
        table.fams.len(),
        started.elapsed().as_secs_f64()
    );
    if !reads.is_empty() {eprintln!("[tessellate hom]\tread matching uses saved repeat settings: minimum {k} tiles, maximum occurrences {occ} (0 = unlimited)");}
    let trace = if trace_all {Some(tessellate::hom::trace::Trace::whole(&ctx))} else {
        trace_bed.map(|p| tessellate::hom::trace::Trace::open(p, &ctx).unwrap_or_else(|e| panic!("cannot load trace BED: {e}")))
    };
    let ev = read::scan_with_trace(&ctx, reads, k, occ, ploidy, gt, trace.as_ref());
    let started = Instant::now();
    let text = table.annotate(&ctx, &ev);
    let lrows = link::rows(&ev.links);
    let ends = link::ends(&ctx, &lrows);
    let writer: Box<dyn Write + '_> = if out == "-" {
        Box::new(io::stdout().lock())
    } else {
        Box::new(std::fs::File::create(out).unwrap_or_else(|e| panic!("cannot create {out}: {e}")))
    };
    let mut o = BufWriter::with_capacity(1 << 24, writer);
    put(&mut o, format!("#cmd\t{cmd}\n").as_bytes());
    for header in [contain::HEADER, family::HEADER, repeat::HEADER, pair::HEADER] {
        put(&mut o, header);
    }
    if ev.reads {
        put(&mut o, tessellate::hom::EVIDENCE_HEADER);
        put(&mut o, link::HEADER);
        put(&mut o, tile::HEADER);
    }
    contain::write(&mut o, &ctx, &ev, &table.contained, &ends);
    if ev.reads {
        link::write(&mut o, &ctx, &lrows);
    }
    family::write(
        &mut o,
        &ctx,
        &ev,
        &table.fams,
        &table.offsets(),
        &table.contained,
        &lrows,
    );
    let mut nd = 0usize;
    for s in &text {
        nd += s.lines().filter(|l| l.starts_with("D\t")).count();
        put(&mut o, s.as_bytes());
    }
    let nt = if ev.reads {
        tile::write(&mut o, &ctx, &ev, min)
    } else {
        0
    };
    if let Some(trace) = trace {
        let n = trace.write(&mut o, &ctx).unwrap_or_else(|e| panic!("cannot write trace rows: {e}"));
        eprintln!("[tessellate hom]\t{n} read provenance rows written");
    }
    if let Err(e) = o.flush() {
        if e.kind() == io::ErrorKind::BrokenPipe {
            return;
        }
        panic!("cannot flush hom report: {e}");
    }
    eprintln!(
        "[tessellate hom]\t{}: {nr} C, {} F, {} R, {nd} D, {} L, {nt} T rows written in {:.2}s",
        out,
        table.fams.len(),
        text.len(),
        lrows.len(),
        started.elapsed().as_secs_f64()
    );
}
