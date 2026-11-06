use clap::{Parser, Subcommand};
mod bwt;
mod fastx;
mod hom;
mod map;
mod output;
mod polish;

fn stem(path: &str) -> Option<&str> {
    let lower = path.to_ascii_lowercase();
    let z = [".gz", ".bgz", ".bgzf", ".bz2", ".zst", ".xz"]
        .iter()
        .find(|z| lower.ends_with(*z))
        .map_or(0, |z| z.len());
    let rest = &lower[..lower.len() - z];
    [".fa", ".fasta", ".fna", ".fas", ".fq", ".fastq", ".tq"]
        .iter()
        .find(|e| rest.ends_with(*e))
        .map(|e| &path[..rest.len() - e.len()])
}

#[derive(Parser)]
#[command(name = "tessellate", version, about = "genomic tessellations")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(clap::Args)]
#[group(id = "CommonOptions")]
struct Options {
    #[arg(short = 'd', long, default_value_t = 1, value_parser = clap::value_parser!(u8).range(0..=8))]
    distance: u8,
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
    #[arg(short, long, required = true)]
    input: Vec<String>,
    #[arg(short, long)]
    output: Option<String>,
    #[arg(short = 'r', long, default_value_t = tessellate::tess::DEFAULT_RATE)]
    rate: usize,
}

impl Options {
    fn output(&self) -> &str {
        self.output.as_deref().unwrap_or_else(|| {
            assert!(self.input.len() == 1, "-o is required with several inputs");
            assert!(self.input[0] != "-", "-o is required with standard input");
            stem(&self.input[0]).unwrap_or(&self.input[0])
        })
    }
}

#[derive(Subcommand)]
enum Cmd {
    Fastx {
        #[command(flatten)]
        opt: Options,
    },
    Bwt {
        #[command(flatten)]
        opt: Options,
    },
    Hom {
        #[command(flatten)]
        opt: hom::Options,
        #[command(flatten)]
        repeats: hom::RepeatArgs,
        #[arg(value_parser = hom::tq_gz_path)]
        reads: Vec<String>,
    },
    Map {
        #[command(flatten)]
        opt: map::Options,
        #[command(flatten)]
        repeats: hom::RepeatArgs,
        #[arg(value_parser = hom::tq_gz_path, required = true)]
        reads: Vec<String>,
    },
    Polish {
        #[command(flatten)]
        opt: polish::Options,
    },
}

impl Cmd {
    fn threads(&self) -> usize {
        match self {
            Self::Fastx {opt} | Self::Bwt {opt} => opt.threads,
            Self::Hom {opt, ..} => opt.threads,
            Self::Map {opt, ..} => opt.threads,
            Self::Polish {opt} => opt.threads,
        }
    }
}

fn main() {
    let cli = Cli::parse();
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.cmd.threads())
        .build_global()
        .expect("failed to initialize thread pool");
    match cli.cmd {
        Cmd::Fastx { opt } => {
            fastx::fastx(&opt.input, opt.output(), opt.rate, usize::from(opt.distance));
        }
        Cmd::Bwt { opt } => {
            let [input] = opt.input.as_slice() else {
                panic!("bwt requires one .tq.gz input; combine raw inputs with tessellate fastx first");
            };
            bwt::bwt(input, opt.output(), opt.rate);
        }
        Cmd::Hom { opt, repeats, reads } => {
            if let Err(error) = hom::hom(opt, &reads, repeats) {
                eprintln!("[tessellate hom]\t{error}");
                std::process::exit(1);
            }
        }
        Cmd::Map { opt, repeats, reads } => {
            if let Err(error) = map::run(opt, repeats, &reads) {
                eprintln!("[tessellate map]\t{error}");
                std::process::exit(1);
            }
        }
        Cmd::Polish { opt } => {
            if let Err(error) = polish::run(opt) {
                eprintln!("[tessellate polish]\t{error}");
                std::process::exit(1);
            }
        }
    }
}
