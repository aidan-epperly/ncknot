//! Prints the invariants of the knots and links in `ncknot::table`.
use ncknot::{table, Ring};

fn main() {
    for &(name, _, _) in table::BRAIDS {
        let d = table::diagram(name).unwrap();
        let kh = if d.crossings() <= 10 {
            d.khovanov(Ring::Fp(2))
                .map(|k| k.total_rank().to_string())
                .unwrap_or_else(|e| e.to_string())
        } else {
            "-".into()
        };
        println!(
            "{name:10} c={:2} comps={} w={:3} V = {}\n           Δ = {}   det = {}   σ = {:?}   H1(Σ2) = {:?}   #col3 = {}   rank Kh(F2) = {kh}",
            d.crossings(),
            d.components(),
            d.writhe(),
            d.jones(),
            d.alexander(),
            d.determinant(),
            d.signature(),
            d.double_branched_cover_homology(),
            d.colorings(3).unwrap(),
        );
    }
}
