//! Timings of each invariant on torus knots and random braids.
use ncknot::{Braid, Ring};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::time::Instant;

fn time<T>(label: &str, f: impl FnOnce() -> T) -> T {
    let t = Instant::now();
    let v = f();
    println!("  {label:28} {:>9.3} s", t.elapsed().as_secs_f64());
    v
}

fn main() {
    for (p, q) in [(3, 7), (4, 5), (5, 6)] {
        let d = Braid::torus(p, q).unwrap().closure();
        println!("T({p},{q}): {} crossings", d.crossings());
        let v = time("jones", || d.jones());
        println!("    V = {v}");
        time("alexander", || d.alexander());
        let s = time("signature", || d.signature().unwrap());
        println!("    σ = {s}");
        time("S_5 homomorphisms", || {
            d.group().count_homomorphisms_to_symmetric(5, true)
        });
    }
    let mut rng = StdRng::seed_from_u64(3);
    let w: Vec<i32> = (0..40)
        .map(|_| rng.gen_range(1..5) * if rng.gen_bool(0.5) { 1 } else { -1 })
        .collect();
    let d = Braid::new(5, &w).unwrap().closure();
    println!("random 5-braid, {} crossings", d.crossings());
    time("jones", || d.jones());
    time("alexander", || d.alexander());
    time("signature", || d.signature().unwrap());
    for n in [10, 12] {
        let w: Vec<i32> = (0..n).map(|i| [1, -2, 3, -1, 2, -3][i % 6]).collect();
        let d = Braid::new(4, &w).unwrap().closure();
        println!("4-braid, {} crossings", d.crossings());
        let k = time("khovanov F_2", || d.khovanov(Ring::Fp(2)).unwrap());
        println!("    total rank {}", k.total_rank());
        let k = time("khovanov Z", || d.khovanov(Ring::Integers).unwrap());
        println!("    {}", k.poincare());
    }
}
