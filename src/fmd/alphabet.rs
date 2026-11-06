use crate::Unit;
use crate::unit::N_CODES;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
const fn slots() -> ([u16; N_CODES], usize) {
    let mut slots = [u16::MAX; N_CODES];
    let (mut id, mut next) = (0, 0);
    while id < N_CODES {
        if crate::unit::UNIT[id].canon as usize == id {
            slots[id] = next;
            next += 1;
        }
        id += 1;
    }
    (slots, next as usize)
}
static SLOTS: [u16; N_CODES] = slots().0;
pub(super) const SIZE: usize = slots().1;
#[inline]
pub(super) fn index(unit: u16) -> usize {
    SLOTS[unit as usize] as usize
}

pub(super) fn complement(unit: u16) -> u16 {
    if unit == 0 { 0 } else { Unit(unit).rc().0 }
}

pub(super) fn symbols() -> Vec<u16> {
    let mut symbols = vec![0];
    for id in 1..N_CODES as u16 {
        let u = Unit(id);
        if u.canon() != u {
            continue;
        }
        let mate = u.rc();
        if u > mate {
            continue;
        }
        symbols.push(u.0);
        if mate != u {
            symbols.push(mate.0);
        }
    }
    symbols
}
pub(super) fn key(code: (u64, u8)) -> u64 {
    code.0 << (64 - 2 * code.1 as usize)
}

pub(super) fn ordered(codes: &[(u64, u8)]) -> Vec<u16> {
    let mut order = symbols();
    order.sort_unstable_by_key(|&u| key(codes[u as usize]));
    order
}
pub(super) fn codes(counts: &[usize]) -> Vec<(u64, u8)> {
    enum Node {
        Leaf(u16),
        Branch(usize, usize),
    }
    let mut nodes = vec![Node::Leaf(0)];
    let mut heap = BinaryHeap::from([Reverse((counts[0], 0usize))]);
    for u in symbols().into_iter().skip(1) {
        let v = complement(u);
        if u > v {
            continue;
        }
        let leaf = nodes.len();
        nodes.push(Node::Leaf(u));
        let (weight, node) = if u == v {
            (counts[index(u)], leaf)
        } else {
            nodes.push(Node::Leaf(v));
            nodes.push(Node::Branch(leaf, leaf + 1));
            (counts[index(u)] + counts[index(v)], leaf + 2)
        };
        heap.push(Reverse((weight, node)));
    }
    let mut separator = vec![false; nodes.len()];
    separator[0] = true;
    while heap.len() > 1 {
        let Reverse((a, mut left)) = heap.pop().unwrap();
        let Reverse((b, mut right)) = heap.pop().unwrap();
        if separator[right] {
            std::mem::swap(&mut left, &mut right);
        }
        let node = nodes.len();
        separator.push(separator[left] || separator[right]);
        nodes.push(Node::Branch(left, right));
        heap.push(Reverse((a + b, node)));
    }
    let mut codes = vec![(0, 0); N_CODES];
    let Reverse((_, root)) = heap.pop().unwrap();
    let mut stack = vec![(root, 0u64, 0usize)];
    while let Some((node, bits, depth)) = stack.pop() {
        match nodes[node] {
            Node::Leaf(u) => {
                assert!(depth > 0 && depth <= 64, "Huffman code exceeds 64 bits");
                codes[u as usize] = (bits << (depth % 2), depth.div_ceil(2) as u8);
            }
            Node::Branch(left, right) => {
                assert!(depth < 64, "Huffman code exceeds 64 bits");
                stack.push((right, (bits << 1) | 1, depth + 1));
                stack.push((left, bits << 1, depth + 1));
            }
        }
    }
    codes
}
