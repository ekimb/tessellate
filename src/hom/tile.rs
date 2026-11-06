use super::{Ctx, Evidence, FIX, code, codeb, f, put, wsum};
use std::io::Write;

pub const HEADER: &[u8] = b"#T\tcontig\ttile\tbase\tdepth(occ<=ploidy)\tagree\tbreaks_left\tbreaks_right\tbases\tcodes(weight>=min)\n";

pub fn write(o: &mut impl Write, ctx: &Ctx, ev: &Evidence, min: f64) -> usize {
    if !ev.reads {
        return 0;
    }
    let lim = (min * FIX as f64) as u64;
    let mut keys: Vec<u64> = ev
        .diffs
        .iter()
        .map(|x| *x.key())
        .chain(ev.brks.iter().map(|x| *x.key() >> 1))
        .collect();
    keys.sort_unstable();
    keys.dedup();
    let mut ntl = 0usize;
    for g in keys {
        let rec = ctx.refs.offs.partition_point(|&x| x as u64 <= g) - 1;
        let lz = ctx.get(rec);
        let ri = (g - ctx.g0(rec)) as usize;
        let (dv, bl, br) = (
            ev.diffs.get(&g),
            ev.brks.get(&(g << 1)),
            ev.brks.get(&(g << 1 | 1)),
        );
        let (wl, wr) = (
            wsum(ev.brks.get(&(g << 1))),
            wsum(ev.brks.get(&(g << 1 | 1))),
        );
        let agree = ev.agreement(g);
        let mut codes: Vec<(u64, String)> = Vec::new();
        if let Some(v) = &dv {
            for x in v.iter().filter(|x| x.1 >= lim) {
                codes.push((x.1, format!("{}:{:.1}", code(x.0), f(x.1))));
            }
        }
        if let Some(v) = &bl {
            for x in v.iter().filter(|x| x.1 >= lim) {
                codes.push((x.1, format!("B<:{}:{:.1}", codeb(x.0), f(x.1))));
            }
        }
        if let Some(v) = &br {
            for x in v.iter().filter(|x| x.1 >= lim) {
                codes.push((x.1, format!("B>:{}:{:.1}", codeb(x.0), f(x.1))));
            }
        }
        if codes.is_empty() {
            continue;
        }
        codes.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        put(
            o,
            format!(
                "T\t{}\t{ri}\t{}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{}\t{}\n",
                String::from_utf8_lossy(lz.name()),
                lz.start(ri).unwrap(),
                f(ev.depp(g)),
                f(agree),
                f(wl),
                f(wr),
                String::from_utf8_lossy(
                    &lz.bases(lz.start(ri).unwrap()..lz.start(ri + 1).unwrap())
                        .unwrap()
                ),
                codes
                    .iter()
                    .map(|x| x.1.as_str())
                    .collect::<Vec<_>>()
                    .join(";")
            )
            .as_bytes(),
        );
        ntl += 1;
    }
    ntl
}
