pub mod fmd;
pub mod hom;
pub mod map;
pub mod polish;
pub mod tess;
pub mod unit;
pub mod utils;

pub use fmd::Fmd;

use unit::{CLASS, MAX_PERIOD, N_CODE, N_CODES, OFFSET, UNIT, pow4, rank};

#[derive(Clone, Copy)]
pub(crate) struct Fields {
    pub(crate) primitive: u16,
    pub(crate) canon: u16,
    pub(crate) rc: u16,
    pub(crate) phase: u8,
    pub(crate) len: u8,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(transparent)]
pub struct Unit(pub u16);

impl Unit {
    pub const N: Unit = Unit(N_CODE);
    pub fn valid(u: u16) -> Option<Unit> {
        (u != 0 && (u as usize) < N_CODES).then_some(Unit(u))
    }
    pub fn canonical(u: u16) -> Option<Unit> {
        Self::valid(u).filter(|x| x.canon() == *x)
    }
    pub const fn from_word(k: usize, word: usize) -> Unit {
        Unit((OFFSET[k] + word) as u16)
    }
    pub const fn is_n(self) -> bool {
        self.0 == N_CODE
    }
    const fn digits(self) -> (usize, usize) {
        let (id, mut l) = (self.0 as usize, MAX_PERIOD);
        while OFFSET[l] > id {
            l -= 1;
        }
        (id - OFFSET[l], l)
    }
    pub const fn word(self) -> usize {
        self.digits().0
    }
    pub const fn base(self, i: usize) -> u8 {
        ((self.word() >> (2 * i)) & 3) as u8
    }
    pub const fn rotate(self, j: usize) -> Unit {
        if self.is_n() {
            return self;
        }
        let (v, l) = self.digits();
        Unit::from_word(l, ((v >> (2 * j)) | (v << (2 * (l - j)))) & (pow4(l) - 1))
    }
    pub const fn revcomp(self) -> Unit {
        let (v, l) = self.digits();
        let (mut w, mut i) = (0, 0);
        while i < l {
            w |= (((v >> (2 * (l - 1 - i))) & 3) ^ 2) << (2 * i);
            i += 1;
        }
        Unit::from_word(l, w)
    }
    const fn period(self) -> usize {
        let (v, l) = self.digits();
        let mut d = 1;
        while d < l {
            if l % d == 0 {
                let mut i = d;
                while i < l && (v >> (2 * i)) & 3 == (v >> (2 * (i - d))) & 3 {
                    i += 1;
                }
                if i == l {
                    return d;
                }
            }
            d += 1;
        }
        l
    }
    const fn less(self, other: Unit) -> bool {
        let ((a, l), (b, _)) = (self.digits(), other.digits());
        let mut i = 0;
        while i < l {
            let (x, y) = (
                rank(((a >> (2 * i)) & 3) as u8),
                rank(((b >> (2 * i)) & 3) as u8),
            );
            if x != y {
                return x < y;
            }
            i += 1;
        }
        false
    }
    const fn least(self) -> Unit {
        let l = self.digits().1;
        let (mut m, mut j) = (self, 1);
        while j < l {
            let w = self.rotate(j);
            if w.less(m) {
                m = w;
            }
            j += 1;
        }
        m
    }
    pub(crate) const fn representative(self) -> Unit {
        let (lw, mw) = (self.least(), self.revcomp().least());
        if lw.0 == mw.0 {
            let l = self.digits().1;
            let mut j = 0;
            while j < l {
                let w = lw.rotate(j);
                if w.0 == w.revcomp().0 {
                    return w;
                }
                j += 1;
            }
            lw
        } else if lw.less(mw) {
            lw
        } else {
            mw.revcomp()
        }
    }
    pub(crate) const fn len_form(self) -> usize {
        self.digits().1
    }
    pub(crate) const fn primitive_form(self) -> Unit {
        let (v, l) = self.digits();
        let d = self.period();
        if d < l {
            Unit::from_word(d, v & (pow4(d) - 1))
        } else {
            self
        }
    }
    fn fields(self) -> &'static Fields {
        &UNIT[self.0 as usize]
    }
    pub fn len(self) -> usize {
        self.fields().len as usize
    }
    pub fn canon(self) -> Unit {
        Unit(self.fields().canon)
    }
    pub fn rc(self) -> Unit {
        Unit(self.fields().rc)
    }
    pub fn primitive(self) -> Unit {
        Unit(self.fields().primitive)
    }
    pub fn phase(self) -> usize {
        self.fields().phase as usize
    }
    pub fn class(self) -> u16 {
        CLASS[self.0 as usize]
    }
    pub fn bases(self) -> impl Iterator<Item = u8> {
        let (v, k) = self.digits();
        (0..k).map(move |i| ((v >> (2 * i)) & 3) as u8)
    }
}

pub trait Tessellable {
    fn unit(&self) -> Unit;
    fn len(&self) -> usize {
        self.unit().len()
    }
    fn canon(&self) -> Unit {
        self.unit().canon()
    }
    fn phase(&self) -> usize {
        self.unit().phase()
    }
    fn is_n(&self) -> bool {
        self.unit().is_n()
    }
}
impl Tessellable for Unit {
    fn unit(&self) -> Unit {
        *self
    }
}
