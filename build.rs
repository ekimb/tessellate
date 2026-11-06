use rayon::prelude::*;
use std::{path::PathBuf, time::Instant};
const MAX_PERIOD: usize = 7;
const MAX_ERRORS: usize = 8;
const MAX_TEMPLATE_LEN: usize = (2 * MAX_ERRORS + 1) * (2 * MAX_PERIOD - 1) + 2 * MAX_ERRORS;
const RADIUS: usize = 2 * MAX_ERRORS;
const CAP: u8 = (RADIUS + 1) as u8;
const LEN: usize = MAX_TEMPLATE_LEN + RADIUS;

const fn pow4(n: usize) -> usize {
    1 << (2 * n)
}

fn offset(l: usize) -> usize {
    1 + (1..l).map(pow4).sum::<usize>()
}

fn get(v: usize, i: usize) -> usize {
    (v >> (2 * i)) & 3
}

fn period(v: usize, l: usize) -> usize {
    (1..l)
        .find(|&d| l % d == 0 && (d..l).all(|i| get(v, i) == get(v, i - d)))
        .unwrap_or(l)
}

fn rotate(v: usize, l: usize, j: usize) -> usize {
    ((v >> (2 * j)) | (v << (2 * (l - j)))) & (pow4(l) - 1)
}

fn revcomp(v: usize, l: usize) -> usize {
    (0..l).fold(0, |w, i| w | ((get(v, l - 1 - i) ^ 2) << (2 * i)))
}

fn identity(v: usize, l: usize) -> usize {
    let r = revcomp(v, l);
    (0..l).map(|j| rotate(v, l, j).min(rotate(r, l, j))).min().unwrap()
}

struct Template {
    identity: usize,
    p: usize,
    bases: [u8; LEN],
}

fn pair(a: &Template, b: &Template, sa: &mut [u8], sb: &mut [u8]) {
    let width = 2 * RADIUS + 1;
    let (mut prev, mut cur) = (vec![CAP; width], vec![CAP; width]);
    for d in 0..width {
        let m = d as isize - RADIUS as isize;
        prev[d] = if (0..=RADIUS as isize).contains(&m) { m as u8 } else { CAP };
    }
    let record = |n: usize, m: usize, v: u8, sa: &mut [u8], sb: &mut [u8]| {
        if n >= 2 * a.p && n <= MAX_TEMPLATE_LEN && m >= 2 * b.p && v < sa[n] {
            sa[n] = v;
        }
        if m >= 2 * b.p && m <= MAX_TEMPLATE_LEN && n >= 2 * a.p && v < sb[m] {
            sb[m] = v;
        }
    };
    for d in 0..width {
        let m = d as isize - RADIUS as isize;
        if m >= 0 {
            record(0, m as usize, prev[d], sa, sb);
        }
    }
    for n in 1..=LEN {
        let mut alive = false;
        for d in 0..width {
            let m = n as isize + d as isize - RADIUS as isize;
            let v = if m < 0 || m as usize > LEN {
                CAP
            } else if m == 0 {
                (n.min(CAP as usize)) as u8
            } else {
                let m = m as usize;
                let diag = prev[d].saturating_add(u8::from(a.bases[n - 1] != b.bases[m - 1]));
                let up = if d + 1 < width { prev[d + 1].saturating_add(1) } else { CAP };
                let left = if d > 0 { cur[d - 1].saturating_add(1) } else { CAP };
                diag.min(up).min(left).min(CAP)
            };
            cur[d] = v;
            if v < CAP {
                alive = true;
                if m >= 0 {
                    record(n, m as usize, v, sa, sb);
                }
            }
        }
        if !alive {
            break;
        }
        std::mem::swap(&mut prev, &mut cur);
    }
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    let begin = Instant::now();
    let templates: Vec<Template> = (1..=MAX_PERIOD)
        .flat_map(|l| (0..pow4(l)).map(move |v| (l, v)))
        .filter(|&(l, v)| period(v, l) == l)
        .map(|(p, v)| Template {
            identity: offset(p) + identity(v, p),
            p,
            bases: std::array::from_fn(|i| get(v, i % p) as u8),
        })
        .collect();
    let data_len: usize = templates.iter().map(|t| MAX_TEMPLATE_LEN + 1 - 2 * t.p).sum();
    let merged = (0..templates.len())
        .into_par_iter()
        .fold(
            || vec![vec![CAP; MAX_TEMPLATE_LEN + 1]; templates.len()],
            |mut sep, i| {
                let mut sa = std::mem::take(&mut sep[i]);
                for j in i + 1..templates.len() {
                    if templates[j].identity == templates[i].identity {
                        continue;
                    }
                    let mut sb = std::mem::take(&mut sep[j]);
                    pair(&templates[i], &templates[j], &mut sa, &mut sb);
                    sep[j] = sb;
                }
                sep[i] = sa;
                sep
            },
        )
        .reduce(
            || vec![vec![CAP; MAX_TEMPLATE_LEN + 1]; templates.len()],
            |mut x, y| {
                for (a, b) in x.iter_mut().zip(&y) {
                    for (u, v) in a.iter_mut().zip(b) {
                        *u = (*u).min(*v);
                    }
                }
                x
            },
        );
    let mut data = Vec::with_capacity(data_len);
    for (t, sep) in templates.iter().zip(&merged) {
        for n in 2 * t.p..=MAX_TEMPLATE_LEN {
            let s = sep[n] as usize;
            assert!(s >= 1, "identical templates of different identity: n={n}");
            data.push(((s - 1) / 2).min(MAX_ERRORS) as u8);
        }
        assert_eq!(*data.last().unwrap(), MAX_ERRORS as u8, "theorem violated");
    }
    assert_eq!(data.len(), data_len);
    let out = PathBuf::from(std::env::var_os("OUT_DIR").unwrap()).join("budgets.bin");
    std::fs::write(&out, &data).unwrap();
    println!(
        "cargo:warning=generated budgets.bin ({} bytes) in {:.0}s",
        data.len(),
        begin.elapsed().as_secs_f64()
    );
}
