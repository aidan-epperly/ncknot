//! Kirby diagrams of families: the symbolic linking numbers, writhes and
//! framings agree with the drawn members, and the cancellation preserves the
//! boundary 3-manifold (counted by homomorphisms of its fundamental group to
//! symmetric groups).

mod common;

use common::spiral::spiral_family;
use ncknot::kirby::{HandleKind, KirbyDiagram};
use ncknot::param::Poly;
use ncknot::{BigInt, Braid, KnotGroup};

/// Homomorphisms from the boundary's fundamental group to `S_k`.
fn boundary_homs(k: &KirbyDiagram, values: &[i64], deg: usize) -> BigInt {
    let ex = k.expand(values).unwrap();
    let f = k.surgery_framings(&ex, values).unwrap();
    let rel = KnotGroup::surgery_relators(&ex.diagram, &f);
    ex.diagram.group().count_homomorphisms_with(&rel, deg, false)
}

#[test]
fn surgery_sign_convention() {
    // +1 surgery on the right-handed trefoil is the Poincaré sphere: its
    // group, the binary icosahedral group, maps onto A_5 in 120 ways; −1
    // surgery gives Σ(2,3,7), whose group has no nontrivial map to S_5.
    // +5 surgery is the lens space L(5, 4): Z/5 has 25 maps to S_5.
    let t = Braid::new(2, &[1, 1, 1]).unwrap().closure();
    let homs = |f: i64| {
        let rel = KnotGroup::surgery_relators(&t, &[f]);
        t.group().count_homomorphisms_with(&rel, 5, false)
    };
    assert_eq!(homs(1), BigInt::from(121));
    assert_eq!(homs(-1), BigInt::from(1));
    assert_eq!(homs(5), BigInt::from(25));
    // 0-surgery on the unknot is S^1 × S^2, group Z: 3! homomorphisms to S_3
    let u = ncknot::Diagram::unknot();
    let rel = KnotGroup::surgery_relators(&u, &[0]);
    assert_eq!(
        u.group().count_homomorphisms_with(&rel, 3, false),
        BigInt::from(6)
    );
}

#[test]
fn spiral_family_symbolic() {
    let k = spiral_family();
    let names = k.param_names();
    let l = k.linking();
    let (r, g, b) = (0, 1, 2);
    let show = |p: &Poly| p.display(&names);
    assert_eq!(show(&l.writhe[g]), "n");
    assert_eq!(show(&l.writhe[b]), "m");
    assert_eq!(show(&l.writhe[r]), "0");
    assert_eq!(show(&l.lk[g][b]), "0");
    assert_eq!(show(&l.lk[g][r]), "-1");
    assert_eq!(show(&l.lk[b][r]), "-1");
    for n in 1..=3 {
        for m in [-3, -1, 1, 3] {
            k.check(&[n, m]).unwrap();
        }
    }
    assert!(k.check(&[2, 2]).is_err(), "m must be odd");
    println!("{}", k.describe());
}

#[test]
fn spiral_family_cancellation() {
    let k = spiral_family();
    let c = k.cancel(0, 2).unwrap();
    println!("{}", c.describe());
    assert_eq!(c.handles.len(), 1);
    assert_eq!(c.handles[0].kind, HandleKind::Framed(Poly::zero()));
    assert_eq!(c.linking().writhe[0].display(&c.param_names()), "2*n*m + n + m");
    for n in 1..=3 {
        for m in [-3, -1, 1, 3] {
            c.check(&[n, m]).unwrap();
            let d = c.expand(&[n, m]).unwrap().diagram;
            assert_eq!(d.components(), 1);
            if d.crossings() <= 30 {
                assert_eq!(d.alexander().to_string(), "-t^-2 + 3 - t^2");
            }
        }
    }
    // the boundary is unchanged
    for (n, m) in [(1, 1), (1, -1), (2, 1), (1, 3)] {
        for deg in [3, 4] {
            let before = boundary_homs(&k, &[n, m], deg);
            let after = boundary_homs(&c, &[n, m], deg);
            println!("n = {n}, m = {m}, S_{deg}: {before} = {after}");
            assert_eq!(before, after);
        }
    }
}
