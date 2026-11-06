pub mod chain;
pub mod extend;
pub mod myers;

use crate::Tessellable;
use crate::hom::refs::Refs;
use crate::tess::Tile;
use crate::unit::BASE;
use std::io;
use std::ops::RangeInclusive;

pub struct Record {
    pub name: Vec<u8>,
    pub letters: Vec<u8>,
    pub tiles: Vec<Tile>,
    tq: Vec<u32>,
    start: Vec<u32>,
    bases: Vec<u8>,
}

impl Record {
    pub fn read(refs: &Refs, rec: usize) -> io::Result<Record> {
        let view = refs.view(rec)?;
        let tiles = view.tiles(0..view.n_tiles(), false)?;
        Ok(Record::new(view.name().to_vec(), tiles))
    }

    pub fn new(name: Vec<u8>, tiles: Vec<Tile>) -> Record {
        let name = name.split(|c| c.is_ascii_whitespace()).next().unwrap_or_default().to_vec();
        let mut record = Record {
            name,
            letters: Vec::new(),
            tiles: Vec::new(),
            tq: Vec::with_capacity(tiles.len() + 1),
            start: Vec::with_capacity(tiles.len() + 1),
            bases: Vec::new(),
        };
        for tile in &tiles {
            record.tq.push(record.letters.len() as u32);
            record.start.push(record.bases.len() as u32);
            if tile.is_n() {
                record.letters.push(b'N');
            } else {
                record.letters.extend(tile.canon().bases().map(|b| BASE[b as usize]));
            }
            record.bases.extend(tile.render());
        }
        record.tq.push(record.letters.len() as u32);
        record.start.push(record.bases.len() as u32);
        record.tiles = tiles;
        record
    }

    pub fn n_tiles(&self) -> usize {
        self.tiles.len()
    }

    pub fn bases(&self, tiles: RangeInclusive<usize>) -> &[u8] {
        &self.bases[self.start[*tiles.start()] as usize..self.start[*tiles.end() + 1] as usize]
    }

    pub fn units(&self) -> Vec<u16> {
        self.tiles.iter().map(|t| t.canon().0).collect()
    }

    pub fn tile_start(&self, tile: usize) -> u32 {
        self.start[tile]
    }

}
