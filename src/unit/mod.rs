use crate::{Fields, Unit};

pub const MAX_PERIOD: usize = 7;
pub const BASE: [u8; 4] = *b"ACTG";
pub const N_CODE: u16 = (OFFSET[MAX_PERIOD] + pow4(MAX_PERIOD)) as u16;
pub const N_CODES: usize = N_CODE as usize + 1;
pub const OFFSET: [usize; MAX_PERIOD + 1] = offsets();
const _: () = assert!(MAX_PERIOD <= 7 && N_CODES <= u16::MAX as usize);

const UNITS: [Fields; N_CODES] = units();
pub(crate) static UNIT: [Fields; N_CODES] = UNITS;
pub(crate) static CLASS: [u16; N_CODES] = classes();

pub(crate) const fn pow4(n: usize) -> usize {
    1 << (2 * n)
}
pub(crate) const fn rank(c: u8) -> u8 {
    c ^ (c >> 1)
}
const fn offsets() -> [usize; MAX_PERIOD + 1] {
    let mut o = [0usize; MAX_PERIOD + 1];
    let (mut l, mut at) = (1, 1);
    while l <= MAX_PERIOD {
        o[l] = at;
        at += pow4(l);
        l += 1;
    }
    o
}
const fn units() -> [Fields; N_CODES] {
    let none = Fields {
        primitive: 0,
        canon: 0,
        rc: 0,
        phase: 0,
        len: 0,
    };
    let mut t = [none; N_CODES];
    let mut id = 1;
    while id < N_CODE as usize {
        let u = Unit(id as u16);
        let p = u.primitive_form();
        if p.0 != u.0 {
            let f = t[p.0 as usize];
            t[id] = Fields {
                primitive: p.0,
                canon: f.canon,
                rc: u.revcomp().0,
                phase: f.phase,
                len: u.len_form() as u8,
            };
        } else if t[id].len == 0 {
            let (c, l) = (u.representative(), u.len_form());
            let mut k = 0;
            while k < 2 {
                let w = if k == 0 { c } else { c.revcomp() };
                let mut q = 0;
                while q < l {
                    let x = w.rotate(q);
                    t[x.0 as usize] = Fields {
                        primitive: x.0,
                        canon: w.0,
                        rc: x.revcomp().0,
                        phase: q as u8,
                        len: l as u8,
                    };
                    q += 1;
                }
                k += 1;
            }
        }
        id += 1;
    }
    t[N_CODE as usize] = Fields {
        primitive: N_CODE,
        canon: N_CODE,
        rc: N_CODE,
        phase: 0,
        len: 1,
    };
    t
}
const fn classes() -> [u16; N_CODES] {
    let mut c = [0u16; N_CODES];
    let mut next = 1u16;
    let mut id = 1;
    while id < N_CODES {
        let f = UNITS[id];
        if f.primitive as usize == id && c[id] == 0 {
            c[id] = next;
            c[f.rc as usize] = next;
            next += 1;
        }
        id += 1;
    }
    let mut id = 1;
    while id < N_CODES {
        let p = UNITS[id].primitive as usize;
        if p != id {
            c[id] = c[p];
        }
        id += 1;
    }
    let mut out = [0u16; N_CODES];
    let mut id = 1;
    while id < N_CODES {
        out[id] = c[UNITS[id].canon as usize];
        id += 1;
    }
    out
}
