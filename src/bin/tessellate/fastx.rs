use super::output::OutputFile;
use rayon::prelude::*;
use std::num::NonZero;
use std::time::Instant;
use tessellate::Tessellable;
use tessellate::tess::{Packed, fastx, parse, store::BgzfWriter, tq};

const OUTPUT_BUFFER_BYTES: usize = 1 << 20;
fn tessellate(batch: Vec<fastx::Record>, rate: usize, edit: usize) -> Vec<(Vec<u16>, Packed)> {
    batch
        .into_par_iter()
        .map(|(name, sequence)| {
            let tiles = parse(&sequence, edit);
            let units = tiles.iter().map(|tile| tile.canon().0).collect();
            let aux = Packed::from_tiles_with_rate(&name, sequence.len(), &tiles, rate)
                .expect("cannot encode reference auxiliary data");
            (units, aux)
        })
        .collect()
}

pub struct Output {
    prefix: String,
    units: tq::BgzfWriter<OutputFile>,
    auxiliary: BgzfWriter<OutputFile>,
    pub records: usize,
    bases: usize,
    tiles: usize,
    started: Instant,
}
impl Output {
    pub fn new(prefix: &str) -> Self {
        Self::with_capacity(prefix, OUTPUT_BUFFER_BYTES)
    }
    fn with_capacity(prefix: &str, buffer_bytes: usize) -> Self {
        let workers = NonZero::new((rayon::current_num_threads() / 2).max(1)).unwrap();
        Self {
            prefix: prefix.to_string(),
            units: tq::BgzfWriter::new(
                OutputFile::new(&format!("{prefix}.tq.gz"), buffer_bytes),
                workers,
            )
            .expect("cannot write binary TQ header"),
            auxiliary: BgzfWriter::new(OutputFile::new(&format!("{prefix}.tss"), buffer_bytes), workers),
            records: 0,
            bases: 0,
            tiles: 0,
            started: Instant::now(),
        }
    }
    pub fn push(&mut self, (units, aux): &(Vec<u16>, Packed)) -> u64 {
        self.units
            .push(aux.name(), units)
            .expect("cannot write binary TQ record");
        let offset = self
            .auxiliary
            .push(aux)
            .expect("cannot write reference auxiliary data");
        self.records += 1;
        self.bases += aux.len();
        self.tiles += units.len();
        offset
    }
    pub fn finish(self, command: &str) {
        let path = format!("{}.tq.gz", self.prefix);
        eprintln!(
            "[tessellate {command}]\tflushing and closing {path} and {}.tss",
            self.prefix
        );
        let finishing = Instant::now();
        let (compressed, aux_bytes) = std::thread::scope(|threads| {
            let auxiliary = threads.spawn(|| {
                let (auxiliary, index) = self
                    .auxiliary
                    .finish()
                    .expect("cannot finish reference auxiliary store");
                let bytes = auxiliary.bytes;
                drop(auxiliary);
                index.save(&format!("{}.tss", self.prefix)).expect("cannot write TSS index");
                eprintln!(
                    "[tessellate {command}]\t{}.tss flushed and closed in {:.2}s",
                    self.prefix,
                    finishing.elapsed().as_secs_f64()
                );
                bytes
            });
            let (units, index) = self.units.finish().expect("cannot finish BGZF unit stream");
            let compressed = units.bytes;
            drop(units);
            index.save(&path).expect("cannot write TQ index");
            eprintln!(
                "[tessellate {command}]\t{path} flushed and closed in {:.2}s",
                finishing.elapsed().as_secs_f64()
            );
            (
                compressed,
                auxiliary.join().expect("auxiliary writer panicked"),
            )
        });
        eprintln!(
            "[tessellate {command}]\t{} records, {} bp, {} tiles; {path} {compressed} bytes, {}.tss {aux_bytes} auxiliary bytes ({:.3} bits/base combined) in {:.2}s",
            self.records,
            self.bases,
            self.tiles,
            self.prefix,
            (compressed + aux_bytes) as f64 * 8.0 / self.bases.max(1) as f64,
            self.started.elapsed().as_secs_f64()
        );
    }
}

pub fn fastx(inputs: &[String], output: &str, rate: usize, edit: usize) {
    let timing = std::env::var_os("TESSELLATE_BATCH_TIMING").is_some();
    let mut writer = Output::new(output);
    let mut stats = vec![(0usize, std::time::Duration::ZERO); inputs.len()];
    fastx::parallel(
        inputs,
        |first, batch| {
            let started = Instant::now();
            let batch = tessellate(batch, rate, edit);
            (first, batch, started.elapsed(), Instant::now())
        },
        |input, (first, batch, mapped, ready)| {
            let waiting = ready.elapsed();
            let started = Instant::now();
            let records = batch.len();
            let bases = writer.bases;
            stats[input].0 += batch.len();
            for record in batch {
                writer.push(&record);
                if writer.records.is_multiple_of(100_000) {
                    eprintln!(
                        "[tessellate fastx]\ttessellated {} sequences in {:.2}s",
                        writer.records,
                        writer.started.elapsed().as_secs_f64()
                    );
                }
            }
            stats[input].1 = writer.started.elapsed();
            if timing {
                eprintln!("[tessellate fastx batch]\tinput={input}\tfirst={first}\treads={records}\tbases={}\tmap_s={:.3}\twait_s={:.3}\twrite_s={:.3}\telapsed_s={:.3}",
                    writer.bases - bases, mapped.as_secs_f64(), waiting.as_secs_f64(), started.elapsed().as_secs_f64(), writer.started.elapsed().as_secs_f64());
            }
            Ok(())
        },
    )
    .unwrap_or_else(|e| panic!("cannot read FASTX input: {e}"));
    for (input, (records, elapsed)) in inputs.iter().zip(stats) {
        eprintln!(
            "[tessellate fastx]\t{input}: {} records tessellated in {:.2}s",
            records,
            elapsed.as_secs_f64()
        );
    }
    assert!(writer.records > 0, "input contains no records");
    writer.finish("fastx");
}
