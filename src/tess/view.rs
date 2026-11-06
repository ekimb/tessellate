use super::codec::invalid;
use super::packed::DecodedBlock;
use super::reference::ReferenceCache;
use super::{PackedView, Tile, tq};
use std::borrow::Cow;
use std::io;
use std::ops::Range;
use std::rc::Rc;

#[derive(Clone)]
pub struct View<'a> {
    aux: PackedView<'a>,
    store: &'a tq::Store,
    record: tq::RecordInfo,
    rec: usize,
    cache: Option<ReferenceCache>,
}
impl<'a> View<'a> {
    pub(crate) fn with_cache(mut self, cache: &ReferenceCache) -> Self {self.cache = Some(*cache); self}
    pub fn units(&self) -> io::Result<Cow<'a, [u16]>> {self.units_range(0..self.n_tiles())}
    pub fn name(&self) -> &[u8] {self.aux.name()}
    pub fn len(&self) -> usize {self.aux.len()}
    pub fn n_tiles(&self) -> usize {self.aux.n_tiles()}
    pub fn indexed(aux: PackedView<'a>, store: &'a tq::Store, record: tq::RecordInfo, rec: usize) -> io::Result<Self> {
        if aux.n_tiles() != record.len() {return Err(invalid("TQ and auxiliary tile counts disagree"));}
        Ok(Self {aux, store, record, rec, cache: Some(store.view_cache())})
    }
    pub fn units_range(&self, range: Range<usize>) -> io::Result<Cow<'a, [u16]>> {
        if range.start > range.end || range.end > self.n_tiles() {return Err(invalid("tile range out of bounds"));}
        self.store.units(self.record, range).map(Cow::Owned)
    }
    pub fn tiles(&self, range: Range<usize>, rc: bool) -> io::Result<Vec<Tile>> {
        let units = self.units_range(range.clone())?;
        let mut tiles = self.tiles_with_units(range, &units, false)?;
        if rc {tiles.reverse(); for tile in &mut tiles {*tile = tile.mirror();}}
        Ok(tiles)
    }
    pub fn tiles_with_units(&self, range: Range<usize>, units: &[u16], rc: bool) -> io::Result<Vec<Tile>> {
        let forward: Cow<'_, [u16]> = if rc {Cow::Owned(units.iter().rev().map(|&u| crate::Unit(u).rc().0).collect())} else {Cow::Borrowed(units)};
        let base = self.start(range.start)?;
        let mut tiles = self.aux.tiles_at_with(range, &forward, base, |block| self.decoded_block(block))?;
        if rc {tiles.reverse(); for tile in &mut tiles {*tile = tile.mirror();}}
        Ok(tiles)
    }

    fn decoded_block(&self, block: usize) -> io::Result<Rc<DecodedBlock>> {
        if let Some(decoded) = self.cache.and_then(|cache| cache.decoded(self.rec, block)) {
            return Ok(decoded);
        }
        let decoded = Rc::new(self.aux.decode_block(block)?);
        if let Some(cache) = self.cache {
            cache.insert_decoded(self.rec, block, decoded.clone());
        }
        Ok(decoded)
    }
    pub fn start(&self, tile: usize) -> io::Result<usize> {
        let rate = self.aux.block_rate();
        if tile > self.n_tiles() {
            return Err(invalid("tile index out of bounds"));
        }
        if tile == self.n_tiles() || tile % rate == 0 {
            return self.aux.start(tile);
        }
        let block = tile / rate;
        if let Some(bounds) = self
            .cache
            .and_then(|cache| cache.boundaries(self.rec, block))
        {
            return Ok(bounds[tile % rate]);
        }
        let first = block * rate;
        let end = first.saturating_add(rate).min(self.n_tiles());
        let decoded = self.decoded_block(block)?;
        let bounds = self.aux.block_boundaries_with_deltas(
            block,
            &self.units_range(first..end)?,
            &decoded.deltas,
        )?;
        let start = bounds[tile % rate];
        if let Some(cache) = self.cache {
            cache.insert_boundaries(self.rec, block, bounds);
        }
        Ok(start)
    }

    pub fn bases(&self, range: Range<usize>) -> io::Result<Vec<u8>> {
        if range.start > range.end || range.end > self.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "base range out of bounds",
            ));
        }
        if range.is_empty() {
            return Ok(Vec::new());
        }
        let tile_at = |position: usize| -> io::Result<usize> {
            let (mut lo, mut hi) = (0, self.n_tiles());
            while lo + 1 < hi {
                let mid = lo + (hi - lo) / 2;
                if self.start(mid)? <= position {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            Ok(lo)
        };
        let first = tile_at(range.start)?;
        let end = tile_at(range.end - 1)? + 1;
        let tiles = self.tiles(first..end, false)?;
        let mut out = Vec::with_capacity(range.len());
        render(&tiles, self.start(first)?, range, &mut out);
        Ok(out)
    }
}

pub(crate) fn render(tiles: &[Tile], start: usize, range: Range<usize>, out: &mut Vec<u8>) {
    let mut s = start;
    for tile in tiles {
        let e = s + tile.span as usize;
        let (lo, hi) = (s.max(range.start), e.min(range.end));
        if lo >= hi {s = e; continue;}
        tile.render_range(lo - s..hi - s, out);
        s = e;
    }
}
