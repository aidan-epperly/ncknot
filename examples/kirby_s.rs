//! The Rasmussen invariant of the knots K(n, m) left after cancelling the
//! 1-handle of the test family (tests/common/spiral.rs), next to the
//! diagram's bound w − s + 2c₊ − 1 (s = Seifert circles here) for the
//! diagram after removing kinks and bigons (`simplify_basic`).  s itself is
//! computed on the diagram from `simplify_for_local`, over 𝔽_p, p = 2³¹ − 1,
//! which is much faster than over ℚ on large members and agrees with it
//! unless p divides the torsion of Bar-Natan homology.
//!
//! `cargo run --release --example kirby_s [n_max] [m_max] [crossing limit]`, or
//! `... --example kirby_s 1,3 1,-3` for chosen members.

#[path = "../tests/common/spiral.rs"]
mod spiral;

use ncknot::kirby::KirbyDiagram;
use ncknot::local_khovanov::{rasmussen_s_mod, simplify_for_local};
use std::time::Instant;

/// One row: crossings and bound of the basic diagram, crossings of the
/// diagram s is computed on, s, and seconds; `None` past the crossing limit.
fn row(k: &KirbyDiagram, n: i64, m: i64, limit: usize) -> Option<String> {
    let raw = k.expand(&[n, m]).unwrap().diagram;
    let d = raw.simplify_basic();
    if d.crossings() > limit {
        return None;
    }
    let w = d.writhe();
    let sc = d.seifert_circles().1 as i64;
    let cp = d.seifert_graph_components(Some(1)) as i64;
    let t = Instant::now();
    let e = simplify_for_local(&raw);
    let s = rasmussen_s_mod(&e, 2147483647);
    Some(format!(
        "{n}\t{m}\t{}\t{}\t{}\t{s:?}\t{:.2}",
        d.crossings(),
        w - sc + 2 * cp - 1,
        e.crossings(),
        t.elapsed().as_secs_f64()
    ))
}

fn main() {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    let k = spiral::spiral_family().cancel(0, 2).unwrap();
    println!("n\tm\tcrossings\tbound\tsimplified\ts\tseconds");
    // explicit members: `1,3 1,-3 ...`
    if raw.iter().any(|x| x.contains(',')) {
        for pair in &raw {
            let v: Vec<i64> = pair.split(',').map(|x| x.parse().unwrap()).collect();
            println!("{}", row(&k, v[0], v[1], usize::MAX).unwrap());
        }
        return;
    }
    let a: Vec<i64> = raw.iter().map(|x| x.parse().unwrap()).collect();
    let (n_max, m_max, limit) = (
        *a.first().unwrap_or(&2),
        *a.get(1).unwrap_or(&5),
        *a.get(2).unwrap_or(&400) as usize,
    );
    for n in 1..=n_max {
        for m in (-m_max..=m_max).filter(|m| m % 2 != 0) {
            if let Some(r) = row(&k, n, m, limit) {
                println!("{r}");
            }
        }
    }
}
