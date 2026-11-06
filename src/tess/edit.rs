use super::{Tile, runs::Piece, sequence::Sequence};
use crate::{Tessellable, Unit, tess::budgets::lookup, unit::rank};
use std::{cell::OnceCell, cmp::Ordering, ops::Range};
pub(super) struct Proposal {
    pub a: usize,
    pub b: usize,
    pub n: usize,
    pub edits: usize,
}
pub(super) fn order_proposals<T>(
    mut items: Vec<T>,
    pieces: usize,
    key: impl Fn(&T) -> (usize, usize, usize),
) -> Vec<T> {
    if items.is_empty() { return items; }
    let mut order: Vec<usize> = (0..items.len()).collect();
    let mut scratch = vec![0; items.len()];
    let mut counts = vec![0; pieces.max(crate::tess::budgets::MAX_ERRORS + 1)];
    for pass in (0..3).rev() {
        let bins = if pass == 2 { crate::tess::budgets::MAX_ERRORS + 1 } else { pieces };
        let counts = &mut counts[..bins];
        counts.fill(0);
        let digit = |i: usize| {
            let (a, b, edits) = key(&items[i]);
            [a, b, edits][pass]
        };
        for &i in &order { counts[digit(i)] += 1; }
        let mut at = 0;
        for count in counts.iter_mut() {
            let n = *count;
            *count = at;
            at += n;
        }
        for &i in &order {
            let at = &mut counts[digit(i)];
            scratch[*at] = i;
            *at += 1;
        }
        std::mem::swap(&mut order, &mut scratch);
    }
    for (to, &from) in order.iter().enumerate() { scratch[from] = to; }
    for i in 0..items.len() {
        while scratch[i] != i {
            let to = scratch[i];
            items.swap(i, to);
            scratch.swap(i, to);
        }
    }
    items.dedup_by(|a, b| {
        let (a0, a1, _) = key(a);
        let (b0, b1, _) = key(b);
        (a0, a1) == (b0, b1)
    });
    items
}
fn content_hash(seq: &Sequence, range: Range<usize>) -> u64 {
    const OFFSET: u64 = 0xcbf29ce484222325;
    const PRIME: u64 = 0x100000001b3;
    let step = |hash: u64, base: u8| (hash ^ u64::from(base)).wrapping_mul(PRIME);
    let forward = range.clone().fold(OFFSET, |h, i| step(h, seq.sym(i)));
    let reverse = range.rev().fold(OFFSET, |h, i| {
        let base = seq.sym(i);
        step(h, if base < 4 { base ^ 2 } else { base })
    });
    forward.min(reverse)
}

pub(super) fn winner_indices(
    props: &[Proposal],
    hash: impl Fn(&Proposal) -> u64,
) -> Vec<usize> {
    let rank = |q: &Proposal| (std::cmp::Reverse(q.n), q.edits);
    let hashes: Vec<OnceCell<u64>> = (0..props.len()).map(|_| OnceCell::new()).collect();
    let compare = |i: usize, j: usize| {
        rank(&props[i]).cmp(&rank(&props[j])).then_with(|| {
            hashes[i].get_or_init(|| hash(&props[i]))
                .cmp(hashes[j].get_or_init(|| hash(&props[j])))
        })
    };
    let width = props.iter().map(|q| q.b - q.a).max().unwrap_or(0);
    let mut winners = Vec::new();
    for (i, q) in props.iter().enumerate() {
        let mut from = i;
        while from > 0 && props[from - 1].a + width >= q.a {
            from -= 1;
        }
        let wins = (from..props.len())
            .take_while(|&j| props[j].a <= q.b)
            .all(|j| i == j || props[j].b < q.a || compare(j, i) == Ordering::Greater);
        if wins {
            winners.push(i);
        }
    }
    winners
}
fn gap(
    seq: &Sequence,
    range: Range<usize>,
    start: Option<Unit>,
    end: Option<Unit>,
    k: usize,
) -> Option<Tile> {
    if range.clone().any(|i| seq.amb(i)) {
        return None;
    }
    let a = range.len();
    let reverse_start = end.map(Unit::rc);
    let reverse_end = start.map(Unit::rc);
    let orientation = range.clone().map(|i| seq.sym(i))
        .cmp(range.clone().rev().map(|i| seq.sym(i) ^ 2))
        .then_with(|| (start, end).cmp(&(reverse_start, reverse_end)));
    let reverse = orientation == Ordering::Greater;
    let tile = gap_oriented(
        a,
        |j| if reverse { seq.sym(range.end - 1 - j) ^ 2 } else { seq.sym(range.start + j) },
        if reverse { reverse_start } else { start },
        if reverse { reverse_end } else { end },
        k,
    )?;
    if orientation == Ordering::Equal && tile != tile.mirror() {
        return None;
    }
    Some(if reverse { tile.mirror() } else { tile })
}

fn gap_oriented(
    a: usize,
    base: impl Fn(usize) -> u8,
    start: Option<Unit>,
    end: Option<Unit>,
    k: usize,
) -> Option<Tile> {
    let p = start.or(end)?.len();
    let mut best: Option<(usize, usize, usize, Unit)> = None;
    for n in a.saturating_sub(k)..=a + k {
        let id = start.unwrap_or_else(|| end.unwrap().rotate((p - n % p) % p));
        if end.is_some_and(|end| id.rotate(n % p) != end) {
            continue;
        }
        let mut edits = a.abs_diff(n);
        for j in 0..a.min(n) {
            if base(j) != id.base(j % p) {
                edits += 1;
                if edits > k {
                    break;
                }
            }
        }
        if edits > k {
            continue;
        }
        let key = (edits, a.abs_diff(n), n);
        if best.as_ref().is_none_or(|b| key < (b.0, b.1, b.2)) {
            best = Some((edits, a.abs_diff(n), n, id));
        }
    }
    let (edits, _, n, id) = best?;
    let mut ex = Vec::with_capacity(edits);
    for j in 0..a.min(n) {
        let actual = base(j);
        if actual != id.base(j % p) {
            ex.push((j as u32, 4 | rank(actual)));
        }
    }
    if a > n {
        for j in n..a {
            ex.push((n as u32, rank(base(j))));
        }
    } else {
        let mut at = a;
        while at < n {
            let len = (n - at).min(8);
            ex.push((at as u32, 8 | (len - 1) as u8));
            at += len;
        }
    }
    let tile = Tile {
        id,
        span: a as u32,
        ex,
    };
    debug_assert_eq!(tile.edit_count(), edits);
    Some(tile)
}

fn end(t: &[Piece], i: usize, n: usize) -> usize {
    t.get(i + 1).map_or(n, |p| p.start)
}
fn run(t: &[Piece], i: usize, n: usize) -> bool {
    !t[i].is_n() && end(t, i, n) - t[i].start > t[i].len()
}

fn witness(
    seq: &Sequence,
    t: &[Piece],
    a: usize,
    b: usize,
    l: usize,
    r: usize,
    k: usize,
) -> Option<Tile> {
    let mut tile = gap(seq, t[a].start..t[l].start, None, Some(t[l].id), k)?;
    let mut n = tile.template_len();
    let mut edits = tile.edit_count();
    let exact = |i: usize| Tile {
        id: t[i].id,
        span: (end(t, i, seq.n) - t[i].start) as u32,
        ex: Vec::new(),
    };
    append_part(&mut tile, &mut n, &mut edits, exact(l));
    if l != r {
        let part = gap(
            seq,
            end(t, l, seq.n)..t[r].start,
            Some(tile.id.rotate(n % tile.id.len())),
            Some(t[r].id),
            k - edits,
        )?;
        append_part(&mut tile, &mut n, &mut edits, part);
        append_part(&mut tile, &mut n, &mut edits, exact(r));
    }
    let part = gap(
        seq,
        end(t, r, seq.n)..end(t, b, seq.n),
        Some(tile.id.rotate(n % tile.id.len())),
        None,
        k - edits,
    )?;
    append_part(&mut tile, &mut n, &mut edits, part);
    if n > u32::MAX as usize || n < 2 * tile.id.len() || edits > lookup(tile.id, 0, n, k) {
        return None;
    }
    debug_assert_eq!(tile.template_len(), n);
    debug_assert_eq!(tile.edit_count(), edits);
    Some(tile)
}

fn append_part(tile: &mut Tile, n: &mut usize, edits: &mut usize, part: Tile) {
    let len = part.template_len();
    debug_assert_eq!(tile.id.rotate(*n % tile.id.len()), part.id);
    *edits += part.edit_count();
    tile.ex.extend(
        part.ex
            .into_iter()
            .map(|(at, op)| ((*n + at as usize) as u32, op)),
    );
    tile.span += part.span;
    *n += len;
}

pub(super) fn parse(seq: &Sequence, k: usize) -> Vec<Tile> {
    let t = seq.raw_pieces(0, seq.n);
    let mut left: Vec<usize> = (0..t.len()).collect();
    let mut right = left.clone();
    let mut candidates = Vec::new();
    let free = |i: usize| !t[i].is_n() && !run(&t, i, seq.n);
    let mut propose = |a, b, l, r| {
        let tile = witness(seq, &t, a, b, l, r, k)?;
        let meta = Proposal {
            a,
            b,
            n: tile.template_len(),
            edits: tile.edit_count(),
        };
        candidates.push((meta, tile));
        Some(())
    };
    for i in 0..t.len() {
        if !run(&t, i, seq.n) {
            continue;
        }
        let p = t[i].len();
        let mut a = i;
        while a > 0 && free(a - 1) && t[i].start - t[a - 1].start <= 4 * p {
            a -= 1;
            if propose(a, i, i, i).is_none() {
                break;
            }
            left[i] = a;
        }
        let mut b = i;
        while b + 1 < t.len() && free(b + 1) && end(&t, b + 1, seq.n) - end(&t, i, seq.n) <= 4 * p {
            b += 1;
            if propose(i, b, i, i).is_none() {
                break;
            }
            right[i] = b;
        }
        if left[i] < i && right[i] > i {
            propose(left[i], right[i], i, i);
        }
    }
    for r in 0..t.len() {
        if !run(&t, r, seq.n) {
            continue;
        }
        let p = t[r].len();
        let mut l = r;
        while l > 0 && t[r].start - end(&t, l - 1, seq.n) <= 4 * p {
            l -= 1;
            if t[l].is_n() {
                break;
            }
            if !run(&t, l, seq.n) {
                continue;
            }
            if t[l].canon() == t[r].canon() {
                for a in [l, left[l]] {
                    for b in [r, right[r]] {
                        propose(a, b, l, r);
                    }
                }
                break;
            }
            if end(&t, l, seq.n) - t[l].start > p + 1 {
                break;
            }
        }
    }
    let candidates = order_proposals(candidates, t.len(), |(p, _)| (p.a, p.b, p.edits));
    let (meta, tiles): (Vec<_>, Vec<_>) = candidates.into_iter().unzip();
    let winners = winner_indices(&meta, |p| content_hash(seq, t[p.a].start..end(&t, p.b, seq.n)));
    let mut tiles: Vec<_> = tiles.into_iter().map(Some).collect();
    let mut winners = winners.into_iter().peekable();
    let mut out = Vec::new();
    let mut i = 0;
    while i < t.len() {
        if let Some(&j) = winners.peek().filter(|&&j| meta[j].a == i) {
            out.push(tiles[j].take().unwrap());
            i = meta[j].b + 1;
            winners.next();
        } else {
            out.push(Tile {
                id: t[i].id,
                span: (end(&t, i, seq.n) - t[i].start) as u32,
                ex: Vec::new(),
            });
            i += 1;
        }
    }
    out
}