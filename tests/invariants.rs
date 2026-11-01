//! Invariants against published values and against each other.

use ncknot::{table, BigInt, Braid, Diagram, IntLaurent, Jones, Ring};
use rand::{rngs::StdRng, Rng, SeedableRng};

fn poly(terms: &[(i64, i64)]) -> IntLaurent {
    IntLaurent::from_terms(terms.iter().map(|&(e, c)| (e, BigInt::from(c))))
}

/// `(name, Alexander coefficients from t^{-d} to t^d, determinant, signature)`.
const KNOTS: &[(&str, &[i64], i64, i64)] = &[
    ("0_1", &[1], 1, 0),
    ("3_1", &[1, -1, 1], 3, -2),
    ("4_1", &[-1, 3, -1], 5, 0),
    ("5_1", &[1, -1, 1, -1, 1], 5, -4),
    ("5_2", &[2, -3, 2], 7, -2),
    ("6_1", &[-2, 5, -2], 9, 0),
    ("6_2", &[-1, 3, -3, 3, -1], 11, -2),
    ("6_3", &[1, -3, 5, -3, 1], 13, 0),
    ("7_1", &[1, -1, 1, -1, 1, -1, 1], 7, -6),
    ("7_2", &[3, -5, 3], 11, -2),
    ("7_3", &[2, -3, 3, -3, 2], 13, -4),
    ("7_4", &[4, -7, 4], 15, -2),
    ("7_5", &[2, -4, 5, -4, 2], 17, -4),
    ("7_6", &[-1, 5, -7, 5, -1], 19, -2),
    ("7_7", &[1, -5, 9, -5, 1], 21, 0),
    ("8_19", &[1, -1, 0, 1, 0, -1, 1], 3, -6),
];

#[test]
fn alexander_determinant_signature_table() {
    for &(name, coeffs, det, sig) in KNOTS {
        let d = table::diagram(name).unwrap();
        let half = (coeffs.len() / 2) as i64;
        let want = IntLaurent::from_terms(
            coeffs
                .iter()
                .enumerate()
                .map(|(j, &c)| (j as i64 - half, BigInt::from(c))),
        );
        assert_eq!(d.alexander(), want, "{name}");
        assert_eq!(d.determinant(), BigInt::from(det), "{name}");
        // the table lists |σ| up to chirality; our braids are the chirality with σ ≤ 0
        assert_eq!(d.signature().unwrap(), sig, "{name}");
        assert_eq!(d.mirror().signature().unwrap(), -sig, "{name}");
        // H_1 of the double branched cover has order det
        let h: BigInt = d.double_branched_cover_homology().iter().product();
        assert_eq!(h, BigInt::from(det), "{name}");
    }
}

#[test]
fn jones_polynomials() {
    let jones_t = |name: &str| table::diagram(name).unwrap().jones().in_t().unwrap();
    assert_eq!(jones_t("3_1"), poly(&[(1, 1), (3, 1), (4, -1)]));
    assert_eq!(
        jones_t("4_1"),
        poly(&[(-2, 1), (-1, -1), (0, 1), (1, -1), (2, 1)])
    );
    assert_eq!(jones_t("8_19"), poly(&[(3, 1), (5, 1), (8, -1)]));
    assert_eq!(
        jones_t("borromean"),
        poly(&[(-3, -1), (-2, 3), (-1, -2), (0, 4), (1, -2), (2, 3), (3, -1)])
    );
    // the positive Hopf link and the 2-component unlink, in t^{1/2}
    let hopf = table::diagram("hopf").unwrap().jones();
    assert_eq!(hopf, Jones::from_sqrt_t(poly(&[(1, -1), (5, -1)])));
    assert_eq!(hopf.to_string(), "-t^(1/2) - t^(5/2)");
    let unlink = Diagram::unlink(2).jones();
    assert_eq!(unlink, Jones::from_sqrt_t(poly(&[(-1, -1), (1, -1)])));
    assert_eq!(Diagram::unknot().jones().in_t().unwrap(), IntLaurent::one());
}

#[test]
fn knottheory_pd_codes() {
    // KnotTheory's PD[Knot[3, 1]] is the left-handed trefoil
    let lh = Diagram::parse_pd("PD[X[1, 4, 2, 5], X[3, 6, 4, 1], X[5, 2, 6, 3]]").unwrap();
    assert_eq!(lh.writhe(), -3);
    assert_eq!(lh.jones().to_string(), "-t^-4 + t^-3 + t^-1");
    assert_eq!(lh.pd_code(), vec![[1, 4, 2, 5], [3, 6, 4, 1], [5, 2, 6, 3]]);
    assert_eq!(lh.signature().unwrap(), 2);
    assert_eq!(lh.mirror().jones(), table::diagram("3_1").unwrap().jones());
    let fig8 = Diagram::parse_pd("X[4,2,5,1], X[8,6,1,5], X[6,3,7,4], X[2,7,3,8]").unwrap();
    assert_eq!(fig8.writhe(), 0);
    assert_eq!(fig8.determinant(), BigInt::from(5));
    assert_eq!(fig8.jones(), table::diagram("4_1").unwrap().jones());
    // a label used three times, an inconsistent orientation, a virtual (non-planar) code
    assert!(Diagram::parse_pd("X[1,1,2,2], X[1,3,3,2]").is_err());
    assert!(Diagram::parse_pd("X[1,2,3,4], X[1,4,3,2]").is_err());
    assert_eq!(
        Diagram::parse_pd("X[1,3,2,4], X[3,1,4,2]").unwrap().components(),
        2
    );
    assert!(Diagram::parse_pd("X[1,2,3,4], X[3,4,1,2]").is_err());
    // a kink is fine
    assert_eq!(
        Diagram::parse_pd("X[1,1,2,2]").unwrap().jones().in_t().unwrap(),
        IntLaurent::one()
    );
}

#[test]
fn conway_and_linking() {
    let z = |name: &str| table::diagram(name).unwrap().conway().unwrap();
    assert_eq!(z("3_1"), poly(&[(0, 1), (2, 1)]));
    assert_eq!(z("4_1"), poly(&[(0, 1), (2, -1)]));
    assert_eq!(z("5_2"), poly(&[(0, 1), (2, 2)]));
    let lk = |name: &str| table::diagram(name).unwrap().linking_matrix();
    assert_eq!(lk("hopf"), vec![vec![0, 1], vec![1, 0]]);
    assert_eq!(lk("T(2,4)")[0][1], 2);
    assert!(lk("borromean").iter().flatten().all(|&x| x == 0));
}

#[test]
fn colorings_and_representations() {
    let d = |name: &str| table::diagram(name).unwrap();
    assert_eq!(d("3_1").colorings(3).unwrap(), BigInt::from(9));
    assert_eq!(d("4_1").colorings(3).unwrap(), BigInt::from(3));
    assert_eq!(d("4_1").colorings(5).unwrap(), BigInt::from(25));
    assert_eq!(Diagram::unlink(2).colorings(3).unwrap(), BigInt::from(9));
    // trefoil: 6 abelian homomorphisms to S_3 and 6 onto it; transitive are
    // the 6 onto S_3 and the 2 onto A_3
    let g = d("3_1").group();
    assert_eq!(g.count_homomorphisms_to_symmetric(3, false), BigInt::from(12));
    assert_eq!(g.count_homomorphisms_to_symmetric(3, true), BigInt::from(8));
    assert_eq!(
        d("4_1").group().count_homomorphisms_to_symmetric(3, false),
        BigInt::from(6)
    );
    // 4_1 has a homomorphism onto A_4 ⊂ S_4 (meridians to 3-cycles)
    let fig8 = d("4_1").group();
    assert!(fig8.count_homomorphisms_to_symmetric(4, true) > BigInt::from(0));
    // the unknot group is ℤ: one homomorphism per element
    assert_eq!(
        Diagram::unknot()
            .group()
            .count_homomorphisms_to_symmetric(4, false),
        BigInt::from(24)
    );
}

#[test]
fn khovanov_trefoil_over_z() {
    let kh = table::diagram("3_1").unwrap().khovanov(Ring::Integers).unwrap();
    let ranks: Vec<_> = kh.ranks.iter().map(|(&k, &r)| (k, r)).collect();
    assert_eq!(ranks, vec![((0, 1), 1), ((0, 3), 1), ((2, 5), 1), ((3, 9), 1)]);
    assert_eq!(kh.torsion.get(&(3, 7)), Some(&vec![BigInt::from(2)]));
    assert_eq!(kh.torsion.len(), 1);
    // over F_2 the Z/2 contributes in degrees 2 and 3
    let f2 = table::diagram("3_1").unwrap().khovanov(Ring::Fp(2)).unwrap();
    assert_eq!(f2.total_rank(), 6);
    assert_eq!(f2.rank(2, 7), 1);
    assert_eq!(f2.rank(3, 7), 1);
    // rationally the figure eight has rank det + 1 = 6 (alternating, thin)
    let q = table::diagram("4_1").unwrap().khovanov(Ring::Integers).unwrap();
    assert_eq!(q.total_rank(), 6);
}

/// `(q + q^{-1}) V(t)` at `t^{1/2} = −q`.
fn jones_to_kh_euler(v: &Jones) -> IntLaurent {
    let q_plus = poly(&[(1, 1), (-1, 1)]);
    v.in_sqrt_t().negate_variable().mul_ref(&q_plus)
}

fn random_braid(rng: &mut StdRng, strands: usize, len: usize) -> Braid {
    loop {
        let w: Vec<i32> = (0..len)
            .map(|_| {
                let g = rng.gen_range(1..strands as i32);
                if rng.gen_bool(0.5) {
                    g
                } else {
                    -g
                }
            })
            .collect();
        let b = Braid::new(strands, &w).unwrap();
        // every strand crosses: no free loops, connected projection
        if (1..strands as i32).all(|g| w.iter().any(|x| x.abs() == g)) {
            return b;
        }
    }
}

/// A random Markov move or braid relation: conjugation, stabilization, or
/// inserting `σ σ^{-1}`.
fn markov_move(rng: &mut StdRng, b: &Braid) -> Braid {
    let n = b.strands();
    let w = b.word();
    match rng.gen_range(0..3) {
        0 if !w.is_empty() => {
            let k = rng.gen_range(0..w.len());
            let mut v = w[k..].to_vec();
            v.extend_from_slice(&w[..k]);
            Braid::new(n, &v).unwrap()
        }
        1 => b.stabilize(rng.gen_bool(0.5)),
        _ => {
            let g = rng.gen_range(1..n as i32);
            let k = rng.gen_range(0..=w.len());
            let mut v = w[..k].to_vec();
            v.extend([g, -g]);
            v.extend_from_slice(&w[k..]);
            Braid::new(n, &v).unwrap()
        }
    }
}

#[test]
fn invariants_survive_markov_moves() {
    let mut rng = StdRng::seed_from_u64(11);
    for _ in 0..12 {
        let strands = rng.gen_range(2..=4);
        let len = rng.gen_range(3..=7);
        let b = random_braid(&mut rng, strands, len);
        let d = b.closure();
        let mut b2 = b.clone();
        for _ in 0..3 {
            b2 = markov_move(&mut rng, &b2);
        }
        let d2 = b2.closure();
        assert_eq!(d.components(), d2.components(), "{b:?} vs {b2:?}");
        assert_eq!(d.jones(), d2.jones(), "{b:?} vs {b2:?}");
        assert_eq!(d.alexander(), d2.alexander(), "{b:?} vs {b2:?}");
        assert_eq!(d.determinant(), d2.determinant(), "{b:?} vs {b2:?}");
        assert_eq!(d.signature().unwrap(), d2.signature().unwrap(), "{b:?} vs {b2:?}");
        assert_eq!(
            d.colorings(3).unwrap(),
            d2.colorings(3).unwrap(),
            "{b:?} vs {b2:?}"
        );
        assert_eq!(d.linking_matrix().len(), d2.linking_matrix().len());
        let mut s = d.double_branched_cover_homology();
        let mut s2 = d2.double_branched_cover_homology();
        s.sort();
        s2.sort();
        assert_eq!(s, s2, "{b:?} vs {b2:?}");
        if d2.crossings() <= 9 {
            let k = d.khovanov(Ring::Fp(2)).unwrap();
            let k2 = d2.khovanov(Ring::Fp(2)).unwrap();
            assert_eq!(k.ranks, k2.ranks, "{b:?} vs {b2:?}");
        }
    }
}

#[test]
fn khovanov_euler_characteristic_is_jones() {
    let mut rng = StdRng::seed_from_u64(5);
    for _ in 0..10 {
        let strands = rng.gen_range(2..=4);
        let len = rng.gen_range(2..=7);
        let d = random_braid(&mut rng, strands, len).closure();
        for ring in [Ring::Fp(2), Ring::Fp(3), Ring::Integers] {
            let kh = d.khovanov(ring).unwrap();
            assert_eq!(
                kh.euler_characteristic(),
                jones_to_kh_euler(&d.jones()),
                "{:?}",
                d.pd_code()
            );
        }
    }
    // free loops: the unlink of 2 has Kh = V ⊗ V
    let kh = Diagram::unlink(2).khovanov(Ring::Integers).unwrap();
    assert_eq!(kh.total_rank(), 4);
    let split = table::diagram("3_1").unwrap().disjoint_union(&Diagram::unknot());
    assert_eq!(split.khovanov(Ring::Fp(3)).unwrap().total_rank(), 8);
}

#[test]
fn mirror_reverse_and_split_unions() {
    let d = table::diagram("5_2").unwrap();
    let v = d.jones().in_sqrt_t().clone();
    assert_eq!(d.mirror().jones().in_sqrt_t(), &v.substitute_power(-1));
    assert_eq!(d.reverse().jones(), d.jones());
    assert_eq!(d.reverse().alexander(), d.alexander());
    // split union: Jones multiplies by −(t^{1/2} + t^{−1/2}), Alexander vanishes
    let u = d.disjoint_union(&table::diagram("3_1").unwrap());
    let factor = poly(&[(-1, -1), (1, -1)]);
    let v31 = table::diagram("3_1").unwrap().jones().in_sqrt_t().clone();
    assert_eq!(u.jones().in_sqrt_t(), &v.mul_ref(&v31).mul_ref(&factor));
    assert!(u.alexander().is_zero());
    assert_eq!(u.determinant(), BigInt::from(0));
    assert!(u.signature().is_err());
    // idle strands of a braid are crossingless loops
    let b = Braid::new(3, &[1, 1, 1]).unwrap().closure();
    assert_eq!(b.components(), 2);
    assert_eq!(b.free_loops(), 1);
    assert_eq!(b.colorings(3).unwrap(), BigInt::from(27));
}

#[test]
fn laurent_arithmetic() {
    let d = poly(&[(2, -1), (-2, -1)]);
    let p = poly(&[(-3, 2), (0, 5), (4, -1)]);
    assert_eq!(p.mul_ref(&d).div_exact(&d), Some(p.clone()));
    assert_eq!(poly(&[(0, 1)]).div_exact(&d), None);
    assert_eq!(p.to_string(), "2*t^-3 + 5 - t^4");
    assert_eq!(poly(&[(1, 1), (3, -2)]).display_as("t", 2), "t^(1/2) - 2*t^(3/2)");
}

#[test]
fn seifert_circles_and_graph() {
    // a braid closure on k strands has k Seifert circles
    let t = Braid::new(2, &[1, 1, 1]).unwrap().closure();
    assert_eq!(t.seifert_circles().1, 2);
    assert_eq!(t.seifert_graph_components(None), 1);
    assert_eq!(t.seifert_graph_components(Some(1)), 1);
    assert_eq!(t.seifert_graph_components(Some(-1)), 2);
    let f = Braid::new(3, &[1, -2, 1, -2]).unwrap().closure();
    assert_eq!(f.seifert_circles().1, 3);
    assert_eq!(f.seifert_graph_components(Some(1)), 2);
    assert_eq!(f.seifert_graph_components(Some(-1)), 2);
    // PD input agrees with the braid (KnotTheory's 3_1 is left-handed)
    let lh = ncknot::Diagram::parse_pd("PD[X[1,4,2,5], X[3,6,4,1], X[5,2,6,3]]").unwrap();
    assert_eq!(lh.seifert_circles().1, 2);
    assert_eq!(lh.seifert_graph_components(Some(-1)), 1);
}

#[test]
fn simplify_keeps_the_knot() {
    // σ1 σ2 σ2⁻¹ σ1 σ1: the bigon goes, the trefoil stays
    let d = Braid::new(3, &[1, 2, -2, 1, 1]).unwrap().closure();
    let s = d.simplify();
    assert!(s.crossings() < d.crossings());
    assert_eq!(s.jones(), d.jones());
    // a kink alone is the unknot
    let k = Diagram::parse_pd("X[1,1,2,2]").unwrap().simplify();
    assert_eq!(k.crossings(), 0);
    assert_eq!(k.components(), 1);
    // random braids with cancelling pairs inserted
    let mut rng = StdRng::seed_from_u64(5);
    for _ in 0..20 {
        let mut w: Vec<i32> = (0..6).map(|_| if rng.gen_bool(0.5) { 1 } else { -2 }).collect();
        let at = rng.gen_range(0..w.len());
        w.insert(at, 2);
        w.insert(at + 1, -2);
        let d = Braid::new(3, &w).unwrap().closure();
        let s = d.simplify();
        assert_eq!(s.jones(), d.jones(), "{w:?}");
        assert_eq!(s.components(), d.components());
    }
}

#[test]
fn reidemeister_three() {
    // every R3 move keeps the link: random braids, all their triangular faces
    let mut rng = StdRng::seed_from_u64(11);
    let mut moves = 0;
    for t in 0..40 {
        let d = random_braid(&mut rng, 3 + t % 3, 6 + t % 7).closure();
        let (v, a, sig) = (d.jones(), d.alexander(), d.signature().unwrap());
        for f in d.r3_faces() {
            let e = d.r3(&f).expect("a move at an R3 face");
            assert_eq!(e.crossings(), d.crossings());
            assert_eq!(e.jones(), v, "{:?}", d.pd_code());
            assert_eq!(e.alexander(), a);
            assert_eq!(e.signature().unwrap(), sig);
            assert_eq!(e.linking_matrix(), d.linking_matrix());
            // the move is undone at the new triangle
            assert!(e
                .r3_faces()
                .iter()
                .any(|g| e.r3(g).is_some_and(|h| h.jones() == v)));
            moves += 1;
        }
    }
    assert!(moves > 50, "only {moves} moves tested");
    // σ1σ2σ1 = σ2σ1σ2: the closure of σ1σ2σ1σ2⁻¹σ1⁻¹σ2⁻¹ is the 3-component
    // unlink, though no kink or bigon shows before an R3 move
    let d = Braid::new(3, &[1, 2, 1, -2, -1, -2]).unwrap().closure();
    let s = d.simplify();
    assert_eq!((s.crossings(), s.free_loops()), (0, 3));
    // random braids with a braid relation inserted reduce to the same size
    for _ in 0..20 {
        let b = random_braid(&mut rng, 3, 5);
        let mut w: Vec<i32> = b.word().to_vec();
        let base = b.closure().simplify().crossings();
        let at = rng.gen_range(0..=w.len());
        w.splice(at..at, [1, 2, 1, -2, -1, -2]);
        let d = Braid::new(3, &w).unwrap().closure();
        let s = d.simplify();
        assert_eq!(s.jones(), d.jones(), "{w:?}");
        assert!(s.crossings() <= base, "{w:?}: {} > {base}", s.crossings());
    }
}

#[test]
fn reidemeister_two_up_and_global_simplify() {
    let mut rng = StdRng::seed_from_u64(12);
    let mut moves = 0;
    for t in 0..30 {
        let d = random_braid(&mut rng, 3 + t % 3, 5 + t % 6).closure();
        let (v, a) = (d.jones(), d.alexander());
        for f in d.faces().iter().filter(|f| f.len() >= 2) {
            for i in 0..f.len() {
                for j in 0..f.len() {
                    for over in [true, false] {
                        let Some(e) = d.r2_up(f, i, j, over) else { continue };
                        assert_eq!(e.crossings(), d.crossings() + 2);
                        assert_eq!(e.jones(), v, "{:?} face {f:?} {i} {j} {over}", d.pd_code());
                        assert_eq!(e.alexander(), a);
                        assert_eq!(e.linking_matrix(), d.linking_matrix());
                        moves += 1;
                    }
                }
            }
        }
    }
    assert!(moves > 500, "only {moves} moves tested");
    // global simplification keeps the knot and never does worse
    for t in 0..15 {
        let d = random_braid(&mut rng, 3 + t % 2, 10).closure();
        let (s, g) = (d.simplify(), d.simplify_global(5));
        assert!(g.crossings() <= s.crossings());
        assert_eq!(g.jones(), d.jones());
        assert_eq!(g.components(), d.components());
    }
}

#[test]
fn strand_pickup() {
    let mut rng = StdRng::seed_from_u64(13);
    let mut moves = 0;
    for t in 0..200 {
        let d = random_braid(&mut rng, 3 + t % 4, 6 + t % 12).closure();
        let Some(e) = d.pickup() else { continue };
        assert!(e.crossings() < d.crossings());
        assert_eq!(e.jones(), d.jones(), "{:?}", d.pd_code());
        assert_eq!(e.alexander(), d.alexander());
        if e.is_connected_projection() {
            assert_eq!(e.signature().unwrap(), d.signature().unwrap());
        }
        // linking numbers, up to the order of the components
        let sorted = |m: Vec<Vec<i64>>| {
            let mut rows: Vec<Vec<i64>> = m
                .into_iter()
                .map(|mut r| {
                    r.sort();
                    r
                })
                .collect();
            rows.sort();
            rows
        };
        assert_eq!(sorted(e.linking_matrix()), sorted(d.linking_matrix()));
        assert_eq!(e.components(), d.components());
        moves += 1;
    }
    assert!(moves > 30, "only {moves} pickups tested");
}
