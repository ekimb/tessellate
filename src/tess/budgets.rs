use crate::unit::MAX_PERIOD;
use crate::{Unit, unit::N_CODE};

pub const MAX_ERRORS: usize = 8;
pub const MAX_TEMPLATE_LEN: usize = (2 * MAX_ERRORS + 1) * (2 * MAX_PERIOD - 1) + 2 * MAX_ERRORS;

const fn offsets() -> [u32; N_CODE as usize + 1] {
    let mut out = [0; N_CODE as usize + 1];
    let mut size = 0;
    let mut id = 1;
    while id < N_CODE as usize {
        out[id] = size;
        let unit = Unit(id as u16);
        if unit.primitive_form().0 == unit.0 {
            size += (MAX_TEMPLATE_LEN + 1 - 2 * unit.len_form()) as u32;
        }
        id += 1;
    }
    out[N_CODE as usize] = size;
    out
}
pub const OFFSETS: [u32; N_CODE as usize + 1] = offsets();
pub const DATA_LEN: usize = OFFSETS[N_CODE as usize] as usize;
const DATA: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/budgets.bin"));
const _: () = assert!(DATA.len() == DATA_LEN, "budgets.bin does not match the table layout");

pub fn lookup(unit: Unit, phase: usize, length: usize, k: usize) -> usize {
    let p = unit.len();
    assert!(k <= MAX_ERRORS && !unit.is_n() && unit.primitive() == unit);
    assert!(phase < p && length >= 2 * p);
    if length > MAX_TEMPLATE_LEN {
        return k;
    }
    let oriented = unit.rotate(phase);
    let index = OFFSETS[oriented.0 as usize] as usize + length - 2 * p;
    k.min(DATA[index] as usize)
}
