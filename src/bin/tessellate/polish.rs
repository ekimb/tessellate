use std::io;
use tessellate::polish::{self, CopyNumber, Inputs};
use tessellate::tess::tq;

fn tq_gz_path(path: &str) -> Result<String, String> {
    tq::stem(path).map(|_| path.to_owned()).ok_or_else(|| "polish inputs must end in .tq.gz and have a matching .tss".to_owned())
}

#[derive(clap::Args)]
pub struct Options {
    #[arg(long, help = "qhom whose T rows select candidate tiles (e.g. from ultra-long reads)")]
    pub detect: String,
    #[arg(long, help = "qhom built with --trace whose X rows supply read votes (e.g. from HiFi reads)")]
    pub decide: String,
    #[arg(long, value_parser = tq_gz_path, help = "Reads of the --decide qhom (.tq.gz with a matching .tss)")]
    pub reads: String,
    #[arg(long, value_parser = tq_gz_path, help = "Tessellated assembly (.tq.gz with a matching .tss)")]
    pub index: String,
    #[arg(long, help = "Assembly FASTA to polish")]
    pub fasta: String,
    #[arg(short, long, help = "Polished FASTA")]
    pub output: String,
    #[arg(long, help = "Optional TSV of applied edits")]
    pub edits: Option<String>,
    #[arg(long, help = "Optional TSV with the outcome and vote counts of every window")]
    pub windows: Option<String>,
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
    #[arg(long, default_value_t = 0.5, help = "Candidate when agree/depth in --detect is below this")]
    pub detect_agree: f64,
    #[arg(long, default_value_t = 3.0, help = "Minimum --detect depth for a candidate")]
    pub detect_depth: f64,
    #[arg(long, help = "Also take candidates where agree/depth in --decide is below this")]
    pub detect_hifi: Option<f64>,
    #[arg(long, default_value_t = 2, help = "Tiles added on each side of a candidate")]
    pub pad: usize,
    #[arg(long, default_value_t = 5, help = "Minimum read votes in a window")]
    pub min_votes: usize,
    #[arg(long, default_value_t = 2, help = "Minimum votes for the chosen string")]
    pub min_support: usize,
    #[arg(long, default_value_t = 2.0, help = "Skip windows with more votes than this times the median depth")]
    pub max_depth: f64,
    #[arg(long, help = "In windows with more votes than this times the median depth, require --tie-margin")]
    pub tie_depth: Option<f64>,
    #[arg(long, default_value_t = 1.5, help = "Votes at the chosen length must be at least this times the votes at the contig length")]
    pub tie_margin: f64,
    #[arg(long, help = "Also take votes from reads whose SMEMs flank a window on both sides")]
    pub flank_votes: bool,
    #[arg(long, help = "Also delete extra repeat copies found in the --decide R rows")]
    pub copy_number: bool,
    #[arg(long, default_value_t = 1000)]
    pub cn_min_len: usize,
    #[arg(long, default_value_t = 8)]
    pub cn_max_copies: usize,
    #[arg(long, default_value_t = 0.8)]
    pub cn_max_ratio: f64,
    #[arg(long, default_value_t = 0.4)]
    pub cn_min_extra: f64,
    #[arg(long, default_value_t = 0.15)]
    pub cn_min_copy_depth: f64,
    #[arg(long, default_value_t = 30.0)]
    pub cn_min_orphan: f64,
    #[arg(long, default_value_t = 15.0)]
    pub cn_min_other: f64,
    #[arg(long, default_value_t = 4.0, help = "Keep copy-number candidates whose copy spacing is at most this times the repeat length")]
    pub cn_max_spacing: f64,
}

pub fn run(opt: Options) -> io::Result<()> {
    let copy_number = opt.copy_number.then_some(CopyNumber {
        min_len: opt.cn_min_len, max_copies: opt.cn_max_copies, max_ratio: opt.cn_max_ratio, min_extra: opt.cn_min_extra,
        min_copy_depth: opt.cn_min_copy_depth, min_orphan: opt.cn_min_orphan, min_other: opt.cn_min_other,
        max_spacing: opt.cn_max_spacing,
    });
    let options = polish::Options {
        detect_agree: opt.detect_agree, detect_depth: opt.detect_depth, detect_hifi: opt.detect_hifi, pad: opt.pad, min_votes: opt.min_votes,
        min_support: opt.min_support, max_depth: opt.max_depth, tie_depth: opt.tie_depth, tie_margin: opt.tie_margin, flank_votes: opt.flank_votes, copy_number,
    };
    let inputs = Inputs {
        detect: &opt.detect, decide: &opt.decide, reads: &opt.reads, index: &opt.index, fasta: &opt.fasta,
        output: &opt.output, edits: opt.edits.as_deref(), windows: opt.windows.as_deref(),
    };
    let summary = polish::polish(&inputs, &options)?;
    let mut decisions: Vec<_> = summary.decisions.iter().collect();
    decisions.sort();
    eprintln!("[tessellate polish]\twindows: {}", decisions.iter().map(|(k, v)| format!("{k} {v}")).collect::<Vec<_>>().join(", "));
    Ok(())
}
