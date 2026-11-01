//! The two-parameter Kirby family of `tests/common/spiral.rs`: cancel its
//! 1-handle, then tabulate the resulting knot's diagram for each `(n, m)`:
//! writhe `w`, Seifert circles `s`, components `c₊` of the Seifert graph
//! without its negative crossings, and `w − s + 2c₊ − 1`; and the same for
//! the mirror image, `−w − s + 2c₋ − 1`.
//!
//! `cargo run --release --example kirby_family [n_max] [m_max]`

#[path = "../tests/common/spiral.rs"]
mod spiral;

fn main() {
    let args: Vec<i64> = std::env::args().skip(1).map(|a| a.parse().unwrap()).collect();
    let n_max = *args.first().unwrap_or(&5);
    let m_max = *args.get(1).unwrap_or(&7);
    let k = spiral::spiral_family().cancel(0, 2).unwrap();
    println!("{}", k.describe());
    println!("n\tm\tcrossings\tn+\tn-\tw\ts\tc+\tw-s+2c+-1\tc-\t-w-s+2c--1");
    for n in 1..=n_max {
        for m in (-m_max..=m_max).filter(|m| m % 2 != 0) {
            let d = k.expand(&[n, m]).unwrap().diagram;
            let (pos, neg) = d.sign_counts();
            let w = d.writhe();
            let s = d.seifert_circles().1 as i64;
            let cp = d.seifert_graph_components(Some(1)) as i64;
            let cn = d.seifert_graph_components(Some(-1)) as i64;
            println!(
                "{n}\t{m}\t{}\t{pos}\t{neg}\t{w}\t{s}\t{cp}\t{}\t{cn}\t{}",
                d.crossings(),
                w - s + 2 * cp - 1,
                -w - s + 2 * cn - 1
            );
        }
    }
}
