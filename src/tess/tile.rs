use crate::{Tessellable, Unit};
use std::{io, ops::Range};
pub const AMB: u8 = 4;
pub const EDIT_BASE: [u8; 4] = *b"ACGT";
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Tile {
    pub id: Unit,
    pub span: u32,
    pub ex: Vec<(u32, u8)>,
}
impl Tessellable for Tile {
    fn unit(&self) -> Unit {
        self.id
    }
}
impl Tile {
    pub fn template_len(&self) -> usize {
        let mut n = self.span as usize;
        for &(_, code) in &self.ex {
            if code >= 8 {
                n += (code & 7) as usize + 1;
            }
        }
        n - self.ex.iter().filter(|e| e.1 < 4).count()
    }
    pub fn edit_count(&self) -> usize {
        self.ex
            .iter()
            .map(|&(_, c)| if c >= 8 { (c & 7) as usize + 1 } else { 1 })
            .sum()
    }
    pub fn render(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.span as usize);
        self.render_range(0..self.span as usize, &mut out);
        out
    }
    pub(crate) fn render_range(&self, range: Range<usize>, out: &mut Vec<u8>) {
        assert!(range.start <= range.end && range.end <= self.span as usize);
        if self.is_n() {
            out.resize(out.len() + range.len(), b'N');
            return;
        }
        let p = self.id.len();
        let word: [u8; crate::unit::MAX_PERIOD] =
            std::array::from_fn(|i| crate::unit::BASE[self.id.base(i % p) as usize]);
        let mut t = 0;
        let mut q = 0;
        let matches = |t: usize, q: usize, n: usize, out: &mut Vec<u8>| {
            for at in q.max(range.start)..(q + n).min(range.end) {
                out.push(word[(t + at - q) % p]);
            }
        };
        for &(offset, code) in &self.ex {
            let offset = offset as usize;
            assert!(offset >= t && code < 16, "invalid tile edit");
            matches(t, q, offset - t, out);
            q += offset - t;
            t = offset;
            if code >= 8 {
                t += (code & 7) as usize + 1;
            } else {
                if range.contains(&q) {
                    out.push(EDIT_BASE[(code & 3) as usize]);
                }
                q += 1;
                t += usize::from(code >= 4);
            }
        }
        assert!(q <= self.span as usize, "edits exceed tile span");
        matches(t, q, self.span as usize - q, out);
    }
    pub fn mirror(&self) -> Tile {
        if self.is_n() {
            return self.clone();
        }
        let n = self.template_len();
        let ex = self
            .ex
            .iter()
            .rev()
            .map(|&(t, code)| {
                let consumed = if code >= 8 {
                    (code & 7) as usize + 1
                } else {
                    usize::from(code >= 4)
                };
                (
                    (n - t as usize - consumed) as u32,
                    if code < 8 { code ^ 3 } else { code },
                )
            })
            .collect();
        Tile {
            id: self.id.rotate(n % self.id.len()).rc(),
            span: self.span,
            ex,
        }
    }
    pub(crate) fn from_wire_iter(
        id: Unit,
        span: u32,
        ex: impl ExactSizeIterator<Item = (u32, u8)>,
    ) -> Self {
        let mut out = Vec::with_capacity(ex.len());
        let (mut deleted, mut inserted) = (0, 0);
        for (q, code) in ex {
            let t = q + deleted - inserted;
            if code >= 4 {
                out.push((t, code & 3));
                inserted += 1;
            } else if code == super::EMPTY {
                match out.last_mut() {
                    Some((before, run)) if *run >= 8 && *run < 15 && *before + (*run & 7) as u32 + 1 == t => *run += 1,
                    _ => out.push((t, 8)),
                }
                deleted += 1;
            } else {
                let actual =
                    crate::unit::BASE[((id.base(t as usize % id.len()) + code + 1) % 4) as usize];
                out.push((
                    t,
                    4 | EDIT_BASE.iter().position(|&b| b == actual).unwrap() as u8,
                ));
            }
        }
        Self { id, span, ex: out }
    }
    pub(crate) fn wire_ex(&self) -> io::Result<Vec<(u32, u8)>> {
        let invalid = || io::Error::new(io::ErrorKind::InvalidInput, "invalid tile edit");
        let mut out = Vec::with_capacity(self.ex.len());
        let (mut deleted, mut inserted) = (0, 0);
        for &(t, code) in &self.ex {
            let q = (t + inserted).checked_sub(deleted).ok_or_else(invalid)?;
            if code >= 16 {
                return Err(invalid());
            }
            if code >= 8 {
                let n = (code & 7) as u32 + 1;
                out.extend(std::iter::repeat_n((q, super::EMPTY), n as usize));
                deleted += n;
            } else if code < 4 {
                out.push((q, 4 | code));
                inserted += 1;
            } else {
                let actual = crate::unit::BASE
                    .iter()
                    .position(|&b| b == EDIT_BASE[(code & 3) as usize])
                    .unwrap() as u8;
                let delta = (actual + 3 - self.id.base(t as usize % self.id.len())) % 4;
                if delta == super::EMPTY {
                    return Err(invalid());
                }
                out.push((q, delta));
            }
        }
        Ok(out)
    }
}
