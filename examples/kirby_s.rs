//! The Rasmussen invariant of the knots K(n, m) left after cancelling the
//! 1-handle of the test family (tests/common/spiral.rs), next to the
//! diagram's bound w − s + 2c₊ − 1 (s = Seifert circles here), for the
//! diagram after removing kinks and bigons.
//!
//! `cargo run --release --example kirby_s [n_max] [m_max] [crossing limit]`

#[path = "../tests/common/spiral.rs"]
mod spiral;

use std::time::Instant;

fn main() {
    let a: Vec<i64> = std::env::args().skip(1).map(|x| x.parse().unwrap()).collect();
    let (n_max, m_max, limit) = (
        *a.first().unwrap_or(&2),
        *a.get(1).unwrap_or(&5),
        *a.get(2).unwrap_or(&400),
    );
    let k = spiral::spiral_family().cancel(0, 2).unwrap();
    println!("n\tm\tcrossings\tbound\ts\tseconds");
    for n in 1..=n_max {
        for m in (-m_max..=m_max).filter(|m| m % 2 != 0) {
            let d = k.expand(&[n, m]).unwrap().diagram.simplify();
            if d.crossings() as i64 > limit {
                continue;
            }
            let w = d.writhe();
            let sc = d.seifert_circles().1 as i64;
            let cp = d.seifert_graph_components(Some(1)) as i64;
            let t = Instant::now();
            let s = d.rasmussen_s();
            println!(
                "{n}\t{m}\t{}\t{}\t{:?}\t{:.2}",
                d.crossings(),
                w - sc + 2 * cp - 1,
                s,
                t.elapsed().as_secs_f64()
            );
        }
    }
}
