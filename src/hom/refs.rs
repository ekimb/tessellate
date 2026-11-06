use crate::Fmd;
use crate::tess::store::Tess;
use crate::tess::{View, tq};
use std::{io, path::PathBuf};

pub struct Refs {
    ts: Tess,
    tq: tq::Store,
    pub offs: Vec<usize>,
    records: Vec<tq::RecordInfo>,
    cache: crate::tess::ReferenceCache,
}
impl Refs {
    pub fn source(fm: &Fmd, index: &str) -> PathBuf {
        std::path::Path::new(index).parent().unwrap_or(std::path::Path::new(".")).join(&fm.source)
    }
    pub fn open(fm: &Fmd, index: &str) -> Self {
        let prefix = Self::source(fm, index);
        let refs = Self::open_prefix(prefix.to_str().expect("invalid source path"))
            .unwrap_or_else(|e| panic!("cannot open {}: {e}", prefix.display()));
        refs.check_index(fm);
        refs
    }
    pub fn open_prefix(prefix: &str) -> io::Result<Self> {
        Self::try_new(Tess::try_open(prefix)?, tq::Store::open(&format!("{prefix}.tq.gz"))?)
    }
    fn check_index(&self, fm: &Fmd) {
        assert_eq!(self.records.len(), fm.n_refs(), "reference record count does not match index");
        for (r, info) in self.records.iter().enumerate() {
            assert_eq!(info.len(), fm.rec_tiles(r), "TQ tile count does not match index");
        }
    }
    fn try_new(ts: Tess, tq: tq::Store) -> io::Result<Self> {
        if ts.offsets.len() != tq.n_records() {return Err(io::Error::other("TQ and TSS record counts disagree"));}
        let mut records = Vec::with_capacity(tq.n_records());
        let mut offs = vec![0];
        for (r, &offset) in ts.offsets.iter().enumerate() {
            let aux = ts.try_view_at(offset as usize)?;
            let (name, info) = tq.record(r)?;
            if name != aux.name() || info.len() != aux.n_tiles() {return Err(io::Error::other("TQ and TSS records disagree"));}
            offs.push(offs[r] + info.len());
            records.push(info);
        }
        Ok(Self {ts, tq, offs, records, cache: Default::default()})
    }
    pub fn n_records(&self) -> usize {self.records.len()}
    pub fn view(&self, rec: usize) -> io::Result<View<'_>> {
        let info = *self.records.get(rec).ok_or_else(|| io::Error::other("reference record out of bounds"))?;
        Ok(View::indexed(self.ts.try_view_at(self.ts.offsets[rec] as usize)?, &self.tq, info, rec)?.with_cache(&self.cache))
    }
    pub fn get(&self, rec: usize) -> View<'_> {
        self.view(rec).expect("auxiliary data does not match TQ")
    }
}
