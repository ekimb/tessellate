use crate::tess::packed::DecodedBlock;
use crate::tess::{AMB, PackedView, Tile};
use std::cell::RefCell;
use std::collections::HashMap;
use std::io;
use std::ops::Range;
use std::rc::Rc;
pub struct Query<'a> {
    pub units: Vec<u16>,
    pub aux: PackedView<'a>,
    boundaries: RefCell<HashMap<usize, Vec<usize>>>,
    blocks: RefCell<HashMap<usize, Rc<DecodedBlock>>>,
}

impl<'a> Query<'a> {
    pub fn new(name: &[u8], units: Vec<u16>, aux: PackedView<'a>) -> io::Result<Self> {
        if aux.name() != name || aux.n_tiles() != units.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{}: TQ record does not match its .tss record",
                    String::from_utf8_lossy(name)
                ),
            ));
        }
        Ok(Self {
            units,
            aux,
            boundaries: RefCell::new(HashMap::new()),
            blocks: RefCell::new(HashMap::new()),
        })
    }

    fn decoded_block(&self, block: usize) -> Rc<DecodedBlock> {
        self.blocks
            .borrow_mut()
            .entry(block)
            .or_insert_with(|| Rc::new(self.aux.decode_block(block).unwrap()))
            .clone()
    }

    pub fn start(&self, tile: usize) -> usize {
        assert!(tile <= self.units.len(), "query tile out of bounds");
        let rate = self.aux.block_rate();
        if tile == self.units.len() || tile % rate == 0 {
            return self.aux.start(tile).unwrap();
        }
        let block = tile / rate;
        let mut boundaries = self.boundaries.borrow_mut();
        let bounds = boundaries.entry(block).or_insert_with(|| {
            let first = block * rate;
            let end = (first + rate).min(self.units.len());
            let decoded = self.decoded_block(block);
            self.aux
                .block_boundaries_with_deltas(block, &self.units[first..end], &decoded.deltas)
                .unwrap()
        });
        bounds[tile % rate]
    }

    pub fn tiles(&self, range: Range<usize>) -> Vec<Tile> {
        let base = self.start(range.start);
        self.aux
            .tiles_at_with(range.clone(), &self.units[range], base, |block| {
                Ok(self.decoded_block(block))
            })
            .unwrap()
    }

    pub fn bases(&self, range: Range<usize>) -> Vec<u8> {
        assert!(range.start <= range.end && range.end <= self.start(self.units.len()));
        if range.is_empty() {
            return Vec::new();
        }
        let tile_at = |position| {
            let (mut lo, mut hi) = (0, self.units.len());
            while lo + 1 < hi {
                let mid = lo + (hi - lo) / 2;
                if self.start(mid) <= position {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        };
        let first = tile_at(range.start);
        let end = tile_at(range.end - 1) + 1;
        let mut start = self.start(first);
        let mut out = Vec::with_capacity(range.len());
        for tile in self.tiles(first..end).iter() {
            let stop = start + tile.span as usize;
            tile.render_range(
                range.start.saturating_sub(start)..range.end.min(stop) - start,
                &mut out,
            );
            start = stop;
        }
        out
    }

    pub fn flank(&self, j: usize, rc: bool) -> u64 {
        let tiles = self.tiles(j..j + 1);
        let tile = &tiles[0];
        let n = tile.span as usize;
        let mut bases = Vec::with_capacity(n.min(29));
        tile.render_range(
            if rc {
                n.saturating_sub(29)..n
            } else {
                0..n.min(29)
            },
            &mut bases,
        );
        let code = |c| match c {
            b'A' => 0,
            b'C' => 1,
            b'T' => 2,
            b'G' => 3,
            _ => AMB,
        };
        if rc {
            super::packb(bases.iter().rev().map(|&c| code(c) ^ 2), n)
        } else {
            super::packb(bases.iter().map(|&c| code(c)), n)
        }
    }
}
