//! Bar-Natan's local algorithm against the cube of resolutions, and the
//! Rasmussen invariant against known values.

use ncknot::{table, Braid, Diagram, Ring};
use rand::{rngs::StdRng, Rng, SeedableRng};

fn same_ranks(d: &Diagram, p: u64, what: &str) {
    let cube = d.khovanov(Ring::Fp(p)).unwrap();
    let local = d.khovanov_local(Ring::Fp(p)).unwrap();
    assert_eq!(local.ranks, cube.ranks, "{what} over F_{p}");
}

#[test]
fn table_knots_match_the_cube() {
    for &(name, _, _) in table::BRAIDS {
        let d = table::diagram(name).unwrap();
        if d.crossings() > 10 {
            continue;
        }
        for p in [2, 3] {
            same_ranks(&d, p, name);
        }
    }
}

fn random_braid(rng: &mut StdRng, strands: usize, len: usize) -> Braid {
    let word: Vec<i32> = (0..len)
        .map(|_| {
            let g = rng.gen_range(1..strands as i32);
            if rng.gen_bool(0.5) {
                g
            } else {
                -g
            }
        })
        .collect();
    Braid::new(strands, &word).unwrap()
}

#[test]
fn random_braids_match_the_cube() {
    let mut rng = StdRng::seed_from_u64(7);
    for t in 0..40 {
        let strands = 2 + t % 3;
        let b = random_braid(&mut rng, strands, 4 + t % 7);
        let d = b.closure();
        same_ranks(&d, 2, &format!("{:?}", b.word()));
        same_ranks(&d, 3, &format!("{:?}", b.word()));
    }
    // a kink and free loops
    let kink = Diagram::parse_pd("X[1,1,2,2]").unwrap();
    same_ranks(&kink, 3, "kink");
    let b = Braid::new(4, &[1, 1, 1]).unwrap().closure(); // two idle strands
    same_ranks(&b, 3, "with free loops");
}

#[test]
fn rasmussen_known_values() {
    let tref = Braid::new(2, &[1, 1, 1]).unwrap().closure();
    assert_eq!(tref.rasmussen_s().unwrap(), 2);
    assert_eq!(tref.mirror().rasmussen_s().unwrap(), -2);
    assert_eq!(table::diagram("4_1").unwrap().rasmussen_s().unwrap(), 0);
    // positive torus knots: s = (p − 1)(q − 1)
    for (p, q) in [(2, 5), (2, 7), (3, 4), (3, 5), (3, 7), (4, 5)] {
        let d = Braid::torus(p, q).unwrap().closure();
        assert_eq!(
            d.rasmussen_s().unwrap(),
            ((p as i64 - 1) * (q as i64 - 1)),
            "T({p},{q})"
        );
    }
    // alternating knots: s = −σ (σ < 0 for positive knots here)
    for &(name, _, _) in table::BRAIDS {
        if !name.contains('_') || name.starts_with("0_") || name == "8_19" {
            continue;
        }
        let d = table::diagram(name).unwrap();
        if d.components() != 1 {
            continue;
        }
        assert_eq!(d.rasmussen_s().unwrap(), -d.signature().unwrap(), "{name}");
    }
    // over F_2 and a large prime too
    let d = Braid::torus(3, 4).unwrap().closure();
    assert_eq!(ncknot::local_khovanov::rasmussen_s_mod(&d, 2).unwrap(), 6);
    assert_eq!(
        ncknot::local_khovanov::rasmussen_s_mod(&d, 2147483647).unwrap(),
        6
    );
}

#[test]
fn rasmussen_is_a_knot_invariant() {
    let mut rng = StdRng::seed_from_u64(11);
    let mut tested = 0;
    while tested < 12 {
        let b = random_braid(&mut rng, 3, 8); // even length: the permutation can be a 3-cycle
        let d = b.closure();
        if d.components() != 1 {
            continue;
        }
        let s = d.rasmussen_s().unwrap();
        // conjugate, stabilize both ways, insert σσ⁻¹
        let w = b.word();
        let conj = Braid::new(3, &[&w[1..], &w[..1]].concat()).unwrap();
        assert_eq!(conj.closure().rasmussen_s().unwrap(), s);
        assert_eq!(b.stabilize(true).closure().rasmussen_s().unwrap(), s);
        assert_eq!(b.stabilize(false).closure().rasmussen_s().unwrap(), s);
        let mut ww = w.to_vec();
        ww.insert(2, 2);
        ww.insert(3, -2);
        assert_eq!(Braid::new(3, &ww).unwrap().closure().rasmussen_s().unwrap(), s);
        assert_eq!(d.mirror().rasmussen_s().unwrap(), -s);
        tested += 1;
    }
}
