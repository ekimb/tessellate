use std::path::{Path, PathBuf};
use std::time::Instant;
use crate::Fmd;
use std::io::{self, Write};
use crate::tess::{store, tq};
pub fn source_path(path: &str, output: &str) -> io::Result<String> {
    let absolute = |p: &str| -> io::Result<PathBuf> {
        let p = Path::new(p);
        let parent = p.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
        Ok(parent.canonicalize()?.join(p.file_name().ok_or_else(|| io::Error::other("missing file name"))?))
    };
    let path = absolute(path)?;
    Ok(if path.parent() == absolute(output)?.parent() {path.file_name().unwrap().to_string_lossy().into_owned()} else {path.to_string_lossy().into_owned()})
}

pub fn build_file(input: &str, output: &str, rate: usize) -> io::Result<Fmd> {
    if rate == 0 {return Err(io::Error::other("suffix-array sampling rate must be positive"));}
    if !input.ends_with(".tq.gz") {return Err(io::Error::other("bwt requires one .tq.gz input with matching .tss; run tessellate fastx first"));}
    let prefix = tq::stem(input).unwrap();
    let store = store::Tess::try_open(prefix)?;
    let (mut text, mut bounds, mut bases) = (Vec::new(), Vec::new(), 0u64);
    let loading = Instant::now();
    tq::parallel(&[input], |_, batch| batch, |_, records| {
            for (name, units) in records {
                let offset = *store.offsets.get(bounds.len()).ok_or_else(|| io::Error::other(format!("{input}: more TQ records than TSS records")))?;
                let aux = store.try_view_at(offset as usize)?;
                if name != aux.name() || units.len() != aux.n_tiles() {return Err(io::Error::other(format!("{input}: TQ and TSS records disagree")));}
                bases += aux.len() as u64;
                bounds.push(text.len());
                text.extend(units);
                text.push(0);
            }
            Ok(())
        },
    )?;
    if bounds.len() != store.offsets.len() || bounds.is_empty() {return Err(io::Error::other("TQ/TSS record counts disagree or reference is empty"));}
    eprintln!("[tessellate bwt]\t{input}: {} records loaded in {:.2}s", bounds.len(), loading.elapsed().as_secs_f64());
    eprintln!("[tessellate bwt]\t{} records, {} tiles loaded for indexing", bounds.len(), text.len() - bounds.len());
    let started = Instant::now();
    let fm = Fmd::build(text, bounds, rate, source_path(prefix, &format!("{output}.tbwt"))?);
    let index = format!("{output}.tbwt");
    let writing = Instant::now();
    eprintln!("[tessellate bwt]\twriting {index}");
    let mut out = io::BufWriter::new(std::fs::File::create(&index)?);
    fm.write(&mut out)?;
    out.flush()?;
    drop(out);
    eprintln!("[tessellate bwt]\t{index} written in {:.2}s", writing.elapsed().as_secs_f64());
    let bytes = std::fs::metadata(&index)?.len();
    let size = if bases > 0 {format!(" ({:.3} bits/base)", (bytes as f64 * 8.0) / (2.0 * bases as f64))} else {String::new()};
    eprintln!("[tessellate bwt]\t{index} {bytes} bytes{size}, including locate samples, in {:.2}s", started.elapsed().as_secs_f64());
    Ok(fm)
}
