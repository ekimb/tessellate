use super::{
    Ctx, Evidence,
    contain::{self, Contained},
    f,
    link::Row,
    put,
    repeat::{Cover, Offsets},
};
use std::collections::{HashMap, HashSet};
use std::io::Write;

pub const HEADER: &[u8] = b"#F\tfamily\tn_members\ttotal_bp\tlink_weight(inside:self:outside)\tform\tperiod\trepresentative\trep_len\tdepth\tredundant(oneway_member:depth_ratio_to_host)\tmembers\n";

fn find(uf: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while uf[r] != r {
        r = uf[r];
    }
    let mut y = x;
    while uf[y] != r {
        let n = uf[y];
        uf[y] = r;
        y = n;
    }
    r
}

pub fn members(fm: &crate::Fmd, cov: &mut Cover, fr: f64) -> Vec<Vec<usize>> {
    let nr = fm.n_refs();
    let mut uf: Vec<usize> = (0..nr).collect();
    for ((c, h), iv) in cov.iter_mut() {
        if contain::holds(iv, fm.rec_tiles(*c), fr) {
            let (ra, rb) = (find(&mut uf, *c), find(&mut uf, *h));
            if ra != rb {
                uf[ra] = rb;
            }
        }
    }
    let mut fam: HashMap<usize, Vec<usize>> = HashMap::new();
    for r in 0..nr {
        let root = find(&mut uf, r);
        fam.entry(root).or_default().push(r);
    }
    let mut fams: Vec<Vec<usize>> = fam.into_values().filter(|v| v.len() >= 2).collect();
    for v in fams.iter_mut() {
        v.sort_unstable();
    }
    fams.sort_by(|a, b| b.len().cmp(&a.len()).then(a[0].cmp(&b[0])));
    fams
}

pub fn write(o: &mut impl Write, ctx: &Ctx, ev: &Evidence, fams: &[Vec<usize>], offs: &Offsets, contained: &Contained, rows: &[Row]) {
    let mut hosts: HashMap<usize, usize> = HashMap::new();
    for (h, _, _, _) in contained.values() {
        *hosts.entry(*h).or_default() += 1;
    }
    for (id, members) in fams.iter().enumerate() {
        let set: HashSet<usize> = members.iter().cloned().collect();
        let (mut intra, mut selfw, mut outw) = (0u64, 0u64, 0u64);
        let mut est: Vec<(f64, u64, usize, usize, String)> = Vec::new();
        for ((a, sa, b, sb), (w, gapsum, n, _, seq)) in rows {
            let (ia, ib) = (set.contains(a), set.contains(b));
            if !ia && !ib {continue;}
            if !(ia && ib) {
                outw += w;
                continue;
            }
            let gap = *gapsum as f64 / *n as f64;
            if a == b {
                if sa != sb {
                    selfw += w;
                    est.push((ctx.len(*a) as f64 + gap, *w, *a, 0, seq.clone()));
                }
                continue;
            }
            intra += w;
            let (aa, bb) = match (sa, sb) {
                (true, false) => (*a, *b),
                (false, true) => (*b, *a),
                _ => continue,
            };
            let Some(ov) = offs.get(&(aa, bb)) else {
                continue;
            };
            let na = ctx.ntiles(aa);
            let Some(t) = ov.iter().filter(|&&t| t >= 0 && (t as usize) < na).max() else {
                continue;
            };
            let lz = ctx.get(aa);
            est.push((
                (lz.start(na).unwrap() - lz.start(*t as usize).unwrap()) as f64 + gap,
                *w,
                aa,
                *t as usize,
                seq.clone(),
            ));
        }
        let windows = intra > selfw;
        let form = if !ev.reads {
            "-"
        } else if intra + selfw == 0 {
            "linear"
        } else if windows {
            "windows"
        } else {
            "loops"
        };
        let mut hist: HashMap<i64, u64> = HashMap::new();
        for (p, w, _, _, _) in &est {
            *hist.entry((*p / 50.0) as i64).or_default() += w;
        }
        let mode = hist
            .iter()
            .max_by_key(|(k, v)| (**v, std::cmp::Reverse(**k)))
            .map(|(&b, _)| b);
        let period = mode.map_or(0.0, |b| {
            let wb: u64 = est
                .iter()
                .filter(|x| (x.0 / 50.0) as i64 == b)
                .map(|x| x.1)
                .sum();
            est.iter()
                .filter(|x| (x.0 / 50.0) as i64 == b)
                .map(|x| x.0 * x.1 as f64)
                .sum::<f64>()
                / wb as f64
        });
        est.retain(|x| x.0 > 0.0);
        let hub = members
            .iter()
            .max_by_key(|&&r| (hosts.get(&r).cloned().unwrap_or(0), ctx.len(r)))
            .cloned()
            .unwrap();
        let (rep, rlen) = match form {
            "windows" => {
                let pick = est
                    .iter()
                    .filter(|x| Some((x.0 / 50.0) as i64) == mode)
                    .max_by_key(|x| x.1);
                match pick {
                    Some((p, _, a, _, g)) => {
                        let n = ctx.len(*a);
                        let body = (*p as usize).saturating_sub(g.len()).min(n);
                        (
                            format!("{}:{}-{}+gap{}", ctx.name(*a), n - body, n, g.len()),
                            *p as usize,
                        )
                    }
                    None => (format!("{}:hub", ctx.name(hub)), ctx.len(hub)),
                }
            }
            "loops" => {
                let mut v: Vec<String> = est
                    .iter()
                    .filter(|x| x.3 == 0 && set.contains(&x.2))
                    .map(|x| format!("{}:0-{}", ctx.name(x.2), ctx.len(x.2)))
                    .collect::<HashSet<_>>()
                    .into_iter()
                    .collect();
                v.sort();
                (
                    if v.is_empty() {
                        format!("{}:hub", ctx.name(hub))
                    } else {
                        v.join(",")
                    },
                    period as usize,
                )
            }
            _ => (
                format!(
                    "{}:host:{}",
                    ctx.name(hub),
                    hosts.get(&hub).cloned().unwrap_or(0)
                ),
                ctx.len(hub),
            ),
        };
        let total: usize = members.iter().map(|&r| ctx.len(r)).sum();
        let mean_depth = if ev.reads {
            let (mut sd, mut sn) = (0u64, 0usize);
            for &r in members {
                let (g0, n) = (ctx.g0(r), ctx.ntiles(r));
                sd += (0..n).map(|t| ev.dep(g0 + t as u64)).sum::<u64>();
                sn += n;
            }
            format!("{:.1}", f(sd) / sn.max(1) as f64)
        } else {
            "-".to_string()
        };
        let names: Vec<String> = members.iter().map(|&r| ctx.name(r)).collect();
        let red: Vec<String> = members
            .iter()
            .filter_map(|&r| {
                let (h, ..) = *contained.get(&r)?;
                if form == "windows" {
                    (r != hub).then(|| format!("{}:{}", ctx.name(r), contain::ratio(ctx, ev, r, h)))
                } else {
                    (!contain::reciprocal(contained, r))
                        .then(|| format!("{}:{}", ctx.name(r), contain::ratio(ctx, ev, r, h)))
                }
            })
            .collect();
        put(o, format!("F\tfamily{id}\t{}\t{total}\t{:.1}:{:.1}:{:.1}\t{form}\t{}\t{rep}\t{rlen}\t{mean_depth}\t{}\t{}\n", members.len(), f(intra), f(selfw), f(outw), if ev.reads { (period as usize).to_string() } else { "-".to_string() }, if red.is_empty() { "-".to_string() } else { red.join(",") }, names.join(",")).as_bytes());
    }
}
