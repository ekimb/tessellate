mod alphabet;
mod samples;
pub mod build;
use samples::TextSamples;
mod wavelet;

use crate::Unit;
use crate::unit::N_CODES;
use crate::utils::{bits_for, compact, de, read_u64, ser, write_u64s};
use libsais::{SuffixArrayConstruction, ThreadCount};
use rayon::prelude::*;
use std::io::Write;
use std::ops::Range;
use std::time::Instant;
use wavelet::Qm;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Place {
    pub rec: usize,
    pub rc: bool,
    pub tile: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Bi {
    pub fwd: usize,
    pub rev: usize,
    pub size: usize,
}

impl Bi {
    #[inline]
    pub fn rows(&self) -> Range<usize> {
        self.fwd..self.fwd + self.size
    }
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }
    #[inline]
    pub(crate) fn swap(self) -> Bi {
        Bi {
            fwd: self.rev,
            rev: self.fwd,
            size: self.size,
        }
    }
}

#[derive(Clone, Copy)]
enum Cursor {
    Bi(Bi),
    At(Place),
}

const EMPTY: Bi = Bi {
    fwd: 0,
    rev: 0,
    size: 0,
};
const MAGIC: &[u8; 5] = b"TBWTS";

fn from_sa<O: Copy + Into<i64> + Sync>(
    text: &[u16],
    symbols: &[u16],
    sa: &[O],
    bounds: &[usize],
    rate: usize,
) -> (Vec<u16>, TextSamples) {
    let started = Instant::now();
    eprintln!("[tessellate bwt]\tconstructing BWT and shared text-position samples");
    let n = text.len();
    let at = |k: usize| sa[k].into() as usize;
    let bwt: Vec<u16> = (0..n)
        .into_par_iter()
        .map(|k| symbols[text[(at(k) + n - 1) % n] as usize])
        .collect();
    let samples = TextSamples::build(sa, bounds, rate);
    eprintln!(
        "[tessellate bwt]\tBWT and shared text-position samples constructed in {:.2}s",
        started.elapsed().as_secs_f64()
    );
    (bwt, samples)
}

#[derive(Default)]
pub struct Fmd {
    pub rate: usize,
    pub source: String,
    lengths: Vec<usize>,
    bwt: Qm,
    c: Vec<usize>,
    text_samples: TextSamples,
    n: usize,
}

impl Fmd {
    pub fn len(&self) -> usize {self.n}
    pub fn lengths(&self) -> &[usize] {&self.lengths}
    pub fn n_refs(&self) -> usize {self.lengths.len()}
    pub fn rec_tiles(&self, r: usize) -> usize {self.lengths[r]}
    pub fn n_tiles(&self) -> usize {(self.len() - 2 * self.n_refs()) / 2}
    pub fn all(&self) -> Bi {
        Bi {
            fwd: 0,
            rev: 0,
            size: self.len(),
        }
    }
    pub fn positions(&self, b: &Bi, max: usize) -> Vec<Place> {
        b.rows().take(max).map(|i| self.locate_place(i)).collect()
    }
    pub fn smems(
        &self,
        end: usize,
        unit: impl Fn(usize) -> u16,
        min: usize,
        lo: usize,
    ) -> Vec<(usize, usize, Bi)> {
        self.smems_while(end, unit, min, lo, || true).unwrap()
    }
    pub fn smems_while(
        &self,
        end: usize,
        unit: impl Fn(usize) -> u16,
        min: usize,
        lo: usize,
        keep_going: impl Fn() -> bool,
    ) -> Option<Vec<(usize, usize, Bi)>> {
        assert!(
            min > 0 && lo > 0,
            "SMEM length and occurrence thresholds must be positive"
        );
        let (mut x, mut out) = (0, Vec::new());
        while x + min <= end {
            if !keep_going() {
                return None;
            }
            let (mut b, mut i) = (self.all(), x + min);
            while i > x {
                if i % 1024 == 0 && !keep_going() {
                    return None;
                }
                let next = self.extend_left(b, unit(i - 1));
                if next.size < lo {
                    break;
                }
                (b, i) = (next, i - 1);
            }
            if i > x {
                x = i;
                continue;
            }
            let mut j = x + min;
            while j < end {
                if j % 1024 == 0 && !keep_going() {
                    return None;
                }
                let next = self.extend_right(b, unit(j));
                if next.size < lo {
                    break;
                }
                (b, j) = (next, j + 1);
            }
            out.push((x, j, b));
            if j == end {
                break;
            }
            let (mut b, mut s) = (self.all(), j + 1);
            while s > x + 1 {
                if s % 1024 == 0 && !keep_going() {
                    return None;
                }
                let next = self.extend_left(b, unit(s - 1));
                if next.size < lo {
                    break;
                }
                (b, s) = (next, s - 1);
            }
            x = s;
        }
        Some(out)
    }
    pub fn smems_no_rescan(
        &self,
        end: usize,
        unit: impl Fn(usize) -> u16,
        min: usize,
        lo: usize,
        occ: usize,
        text: impl Fn(usize, usize) -> u16,
    ) -> Vec<(usize, usize, usize, Vec<Place>)> {
        assert!(min > 0 && lo > 0, "SMEM length and occurrence thresholds must be positive");
        let at = |p: Place, pos: usize| -> Option<u16> {
            let n = self.lengths[p.rec];
            if pos >= n {return None;}
            Some(if p.rc {alphabet::complement(text(p.rec, n - 1 - pos))} else {text(p.rec, pos)})
        };
        let promote = |b: Bi| if lo == 1 && b.size == 1 {Cursor::At(self.locate_place(b.fwd))} else {Cursor::Bi(b)};
        let left = |c: Cursor, x: u16| -> Option<Cursor> {
            if Unit::canonical(x).is_none() {return None;}
            match c {
                Cursor::At(p) => (p.tile > 0 && at(p, p.tile - 1) == Some(x)).then_some(Cursor::At(Place {tile: p.tile - 1, ..p})),
                Cursor::Bi(b) => {let b = self.extend_left(b, x); (b.size >= lo).then_some(Cursor::Bi(b))}
            }
        };
        let right = |c: Cursor, x: u16, len: usize| -> Option<Cursor> {
            if Unit::canonical(x).is_none() {return None;}
            match c {
                Cursor::At(p) => (at(p, p.tile + len) == Some(x)).then_some(c),
                Cursor::Bi(b) => {let b = self.extend_right(b, x); (b.size >= lo).then(|| promote(b))}
            }
        };
        let (mut x, mut out, mut carry) = (0, Vec::new(), None);
        while x + min <= end {
            let (mut c, mut j) = match carry.take() {
                Some(state) => state,
                None => {
                    let (mut c, mut i) = (Cursor::Bi(self.all()), x + min);
                    while i > x {
                        match left(c, unit(i - 1)) {Some(next) => (c, i) = (next, i - 1), None => break}
                    }
                    if i > x {
                        x = i;
                        continue;
                    }
                    (c, x + min)
                }
            };
            while j < end {
                match right(c, unit(j), j - x) {Some(next) => (c, j) = (next, j + 1), None => break}
            }
            let (size, places) = match c {
                Cursor::At(p) => (1, vec![p]),
                Cursor::Bi(b) => (b.size, if occ > 0 && b.size > occ {Vec::new()} else {self.positions(&b, b.size)}),
            };
            out.push((x, j, size, places));
            if j == end {
                break;
            }
            let (mut c, mut s) = (Cursor::Bi(self.all()), j + 1);
            while s > x + 1 {
                match left(c, unit(s - 1)) {Some(next) => (c, s) = (next, s - 1), None => break}
            }
            if j + 1 - s >= min {
                carry = Some((c, j + 1));
            }
            x = s;
        }
        out
    }
    pub fn build(mut text: Vec<u16>, bounds: Vec<usize>, rate: usize, source: String) -> Fmd {
        let preparing = Instant::now();
        let m = text.len();
        eprintln!("[tessellate bwt]\tpreparing text and alphabet");
        debug_assert!(rate > 0, "suffix-array sampling rate must be positive");
        debug_assert!(!text.is_empty() && text.last() == Some(&0), "text must end with a record separator");
        debug_assert!(text.iter().all(|&u| u == 0 || Unit::canonical(u).is_some()), "index text must contain canonical unit IDs");
        debug_assert_eq!(bounds.first(), Some(&0), "first record must start at zero");
        debug_assert!(bounds.windows(2).all(|w| w[0] < w[1]) && bounds.last().unwrap() < &m, "invalid record boundaries");
        let lengths = bounds.iter().enumerate().map(|(r, &start)| bounds.get(r + 1).copied().unwrap_or(m) - start - 1).collect();
        text.extend_from_within(..m - 1);
        text[m..].reverse();
        for u in &mut text[m..] {*u = alphabet::complement(*u);}
        text.push(0);
        let n = text.len();
        let mut counts = vec![0usize; alphabet::SIZE];
        for &u in &text {counts[alphabet::index(u)] += 1;}
        let codes = alphabet::codes(&counts);
        let order = alphabet::ordered(&codes);
        let (mut c, mut total) = (vec![0usize; alphabet::SIZE], 0);
        let mut ranks = vec![0u16; N_CODES];
        for (rank, &u) in order.iter().enumerate() {
            let index = alphabet::index(u);
            ranks[u as usize] = rank as u16;
            c[index] = total;
            total += counts[index];
        }
        debug_assert_eq!(total, n);
        for u in &mut text {*u = ranks[*u as usize];}
        eprintln!("[tessellate bwt]\ttext and alphabet prepared: {n} symbols over both strands in {:.2}s", preparing.elapsed().as_secs_f64());
        let sorting = Instant::now();
        eprintln!("[tessellate bwt]\tsorting suffixes");
        let (bwt, samples) = if n < i32::MAX as usize {
            let sa = SuffixArrayConstruction::for_text(&text).in_owned_buffer32().multi_threaded(ThreadCount::openmp_default()).generalized_suffix_array().run().expect("libsais failed");
            from_sa(&text, &order, sa.suffix_array(), &bounds, rate)

        } else {
            let sa = SuffixArrayConstruction::for_text(&text).in_owned_buffer64().multi_threaded(ThreadCount::openmp_default()).generalized_suffix_array().run().expect("libsais failed");
            from_sa(&text, &order, sa.suffix_array(), &bounds, rate)
        };
        eprintln!("[tessellate bwt]\tsuffixes sorted in {:.2}s", sorting.elapsed().as_secs_f64());
        drop(text);
        drop(bounds);
        let wavelet = Instant::now();
        eprintln!("[tessellate bwt]\tconstructing wavelet and rank structures");
        let bwt = Qm::build(bwt, codes);
        eprintln!("[tessellate bwt]\twavelet and rank structures constructed in {:.2}s", wavelet.elapsed().as_secs_f64());
        let fm = Fmd {bwt, c, text_samples: samples, n, rate, source, lengths};
        fm
    }
    pub fn extend_left(&self, b: Bi, x: u16) -> Bi {
        if b.is_empty() || Unit::canonical(x).is_none() {
            return EMPTY;
        }
        let (before, size, less) = self.bwt.range(b.fwd, b.fwd + b.size, x);
        Bi {
            fwd: self.c[alphabet::index(x)] + before,
            rev: b.rev + less,
            size,
        }
    }
    pub fn extend_right(&self, b: Bi, x: u16) -> Bi {
        if Unit::canonical(x).is_none() {
            return EMPTY;
        }
        self.extend_left(b.swap(), alphabet::complement(x)).swap()
    }
    pub fn locate_place(&self, i: usize) -> Place {
        let (place, _, _) = self.locate_trace(i);
        debug_assert!(place.tile < self.rec_tiles(place.rec), "a record's end");
        place
    }
    fn locate_trace(&self, mut i: usize) -> (Place, usize, usize) {
        let (mut steps, mut periods) = (0, 0);
        loop {
            if let Some(mut place) = self.text_samples.sampled_place(i) {
                place.tile += steps;
                return (place, steps, periods);
            }
            let (v, r) = self.bwt.access_rank(i);
            debug_assert_ne!(v, 0, "record start is missing its text-position sample");
            i = self.c[alphabet::index(v)] + r;
            steps += 1;
            periods += Unit(v).len();
        }
    }
    pub fn write<W: Write>(&self, mut w: W) -> std::io::Result<()> {
        if self.source.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "index requires a TQ/TSS source prefix",
            ));
        }
        assert_eq!(
            self.bwt.packed.len(),
            self.bwt.offs.len(),
            "an index read back keeps no digits to write"
        );
        w.write_all(MAGIC)?;
        write_u64s(
            &mut w,
            [self.rate, self.n_refs(), self.n, self.source.len()],
        )?;
        w.write_all(self.source.as_bytes())?;
        write_u64s(&mut w, self.c.iter().copied())?;
        for u in alphabet::symbols() {
            let (code, digits) = self.bwt.codes[u as usize];
            w.write_all(&code.to_le_bytes())?;
            w.write_all(&[digits])?;
        }
        for p in &self.bwt.packed {
            write_u64s(&mut w, [p.len()])?;
            write_u64s(&mut w, p.iter().map(|&x| x as usize))?;
        }
        self.text_samples.write(&mut w)?;
        w.flush()
    }

    pub fn read(mut b: &[u8]) -> Self {
        assert!(
            b.starts_with(MAGIC),
            "unsupported .tbwt format; rebuild with tessellate bwt"
        );
        b = &b[MAGIC.len()..];
        let (rate, nr, n, source_len) = (
            read_u64(&mut b),
            read_u64(&mut b),
            read_u64(&mut b),
            read_u64(&mut b),
        );
        assert!(source_len > 0, "missing TQ/TSS source prefix");
        let source = String::from_utf8(b[..source_len].to_vec()).unwrap();
        b = &b[source_len..];
        let c: Vec<usize> = (0..alphabet::SIZE).map(|_| read_u64(&mut b)).collect();
        let mut codes = vec![(0, 0); N_CODES];
        for u in alphabet::symbols() {
            let code = read_u64(&mut b) as u64;
            let digits = b[0];
            b = &b[1..];
            assert!((1..=32).contains(&digits), "invalid Huffman code length");
            codes[u as usize] = (code, digits);
        }
        let offs = Qm::children(&codes);
        let packed: Vec<Vec<u64>> = (0..offs.len())
            .map(|_| {
                let m = read_u64(&mut b);
                let w = b[..8 * m]
                    .par_chunks_exact(8)
                    .map(|x| u64::from_le_bytes(x.try_into().unwrap()))
                    .collect();
                b = &b[8 * m..];
                w
            })
            .collect();
        let mut bwt = Qm::with(packed, offs, codes);
        bwt.packed = Vec::new();
        let (samples, lengths) = TextSamples::read(&mut b, nr, n, rate);
        assert!(b.is_empty(), "trailing data in .tbwt");
        Fmd {
            bwt,
            c,
            text_samples: samples,
            n,
            rate,
            source,
            lengths,
        }
    }
}
