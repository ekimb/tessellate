use super::*;
use sucds::bit_vectors::{Access as _, BitVector, Rank, Rank9Sel};
use sucds::int_vectors::{Access as _, CompactVector};

#[derive(Default)]
pub(super) struct TextSamples {
    marked: Rank9Sel,
    records: CompactVector,
    tiles: CompactVector,
}
impl TextSamples {
    pub(super) fn build<O: Copy + Into<i64> + Sync>(
        sa: &[O],
        bounds: &[usize],
        rate: usize,
    ) -> Self {
        let n = sa.len();
        let mut text = BitVector::from_bit(false, n);
        let (mut count, mut max_len) = (0, 0);
        for (rec, &start) in bounds.iter().enumerate() {
            let len = bounds.get(rec + 1).copied().unwrap_or(n / 2) - start - 1;
            max_len = max_len.max(len);
            count += 2 * (len.div_ceil(rate) + 1);
            for tile in (0..len).step_by(rate).chain(std::iter::once(len)) {
                text.set_bit(start + tile, true).unwrap();
                text.set_bit(n - 1 - start - tile, true).unwrap();
            }
        }
        let selected: Vec<_> = sa
            .par_iter()
            .enumerate()
            .filter_map(|(row, &pos)| {
                let pos = pos.into() as usize;
                text.access(pos).unwrap().then(|| {
                    let rc = pos >= n / 2;
                    let p = if rc { n - 1 - pos } else { pos };
                    let rec = bounds.partition_point(|&start| start <= p) - 1;
                    let tile = p - bounds[rec];
                    let tile = if rc {
                        let len = bounds.get(rec + 1).copied().unwrap_or(n / 2) - bounds[rec] - 1;
                        len - tile
                    } else {
                        tile
                    };
                    (row, Place { rec, rc, tile })
                })
            })
            .collect();
        assert_eq!(selected.len(), count);
        drop(text);
        let mut selected_rows = selected.iter().map(|&(row, _)| row).peekable();
        let marked = Rank9Sel::from_bits((0..n).map(|row| {
            if selected_rows.peek() == Some(&row) {
                selected_rows.next();
                true
            } else {
                false
            }
        }));
        Self {
            marked,
            records: compact(
                bits_for((2 * bounds.len()).saturating_sub(1)),
                selected
                    .iter()
                    .map(|&(_, p)| (2 * p.rec + usize::from(p.rc)) as u64),
            ),
            tiles: compact(
                bits_for(max_len),
                selected.iter().map(|&(_, p)| p.tile as u64),
            ),
        }
    }

    #[inline]
    pub(super) fn sampled_place(&self, row: usize) -> Option<Place> {
        if !self.marked.access(row).unwrap() {
            return None;
        }
        let slot = self.marked.rank1(row).unwrap();
        let record = self.records.access(slot).unwrap() as usize;
        Some(Place {
            rec: record / 2,
            rc: record % 2 != 0,
            tile: self.tiles.access(slot).unwrap() as usize,
        })
    }
    pub(super) fn write<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        ser(&self.marked, &mut *w)?;
        ser(&self.records, &mut *w)?;
        ser(&self.tiles, w)
    }
    pub(super) fn read(bytes: &mut &[u8], nr: usize, n: usize, rate: usize) -> (Self, Vec<usize>) {
        assert!(rate > 0, "invalid sampling rate");
        let s = Self {
            marked: de(bytes),
            records: de(bytes),
            tiles: de(bytes),
        };
        assert_eq!(
            s.marked.len(),
            n,
            "sample bitmap length disagrees with index"
        );
        let count = s.marked.rank1(n).unwrap();
        assert_eq!(
            s.records.len(),
            count,
            "sample record count disagrees with bitmap"
        );
        assert_eq!(
            s.tiles.len(),
            count,
            "sample tile count disagrees with bitmap"
        );
        let mut ends = vec![None::<usize>; 2 * nr];
        for slot in 0..count {
            let record = s.records.access(slot).unwrap() as usize;
            let tile = s.tiles.access(slot).unwrap() as usize;
            let end = ends.get_mut(record).expect("sample record out of bounds");
            *end = Some(end.unwrap_or(0).max(tile));
        }
        let lengths: Vec<_> = ends
            .chunks_exact(2)
            .map(|pair| {
                let len = pair[0].expect("missing forward record samples");
                assert_eq!(
                    pair[1],
                    Some(len),
                    "sampled strands disagree on record length"
                );
                len
            })
            .collect();
        assert_eq!(n % 2, 0, "invalid two-strand text length");
        assert_eq!(
            lengths.iter().sum::<usize>() + nr,
            n / 2,
            "sampled record lengths disagree with index"
        );
        assert_eq!(
            count,
            lengths
                .iter()
                .map(|n| 2 * (n.div_ceil(rate) + 1))
                .sum::<usize>(),
            "sample bitmap count disagrees with text blocks"
        );
        (s, lengths)
    }
}
