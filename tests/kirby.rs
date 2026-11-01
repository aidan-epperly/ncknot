//! Kirby diagrams of families: the symbolic linking numbers, writhes and
//! framings agree with the drawn members, and the cancellation preserves the
//! boundary 3-manifold (counted by homomorphisms of its fundamental group to
//! symmetric groups).

use ncknot::kirby::{HandleKind, KirbyDiagram, NodeKind, Units};
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

/// A two-parameter family: a dotted circle R; a 0-framed handle G winding
/// n + 1 times up and n times down through R (a spiral: an (n+1)-bundle and
/// an n-bundle joined by caps, one strand peeling off at each end); a
/// 0-framed handle B through R once, with m (odd) half-twists on two of its
/// strands, clasping the G bundles above and below R.
fn spiral_family() -> KirbyDiagram {
    let mut k = KirbyDiagram::new(&[("n", None), ("m", Some(true))]);
    let n = k.param("n");
    let m = k.param("m");
    let r = k.dotted("R");
    let g = k.framed("G", Poly::zero());
    let b = k.framed("B", Poly::zero());
    let cross = |k: &mut KirbyDiagram| k.node(NodeKind::Cross, 4);
    // crossings of R (top edge: R under; bottom edge: R over)
    let (rut, rdt, rbt, rbb, rdb, rub) = (
        cross(&mut k),
        cross(&mut k),
        cross(&mut k),
        cross(&mut k),
        cross(&mut k),
        cross(&mut k),
    );
    let split_top = k.node(NodeKind::Split, 3);
    let split_bot = k.node(NodeKind::Split, 3);
    let xa = cross(&mut k); // the peeled-off strand over the lower n-bundle
    let (bnt, bdt, bdb, bnb) = (cross(&mut k), cross(&mut k), cross(&mut k), cross(&mut k));
    let split_in = k.node(NodeKind::Split, 3);
    let bx = k.node(
        NodeKind::Twist {
            half_twists: m.clone(),
            units: Units::Ribbons(vec![1, 1]),
        },
        2,
    );
    let split_out = k.node(NodeKind::Split, 3);
    let one = Poly::constant(1);
    let rr = [(r, 1, one.clone())];
    // R: along the top edge west to east, back along the bottom
    k.edge(&rr, (rut, 2), (rdt, 0));
    k.edge(&rr, (rdt, 2), (rbt, 0));
    k.edge(&rr, (rbt, 2), (rbb, 1));
    k.edge(&rr, (rbb, 3), (rdb, 1));
    k.edge(&rr, (rdb, 3), (rub, 1));
    k.edge(&rr, (rub, 3), (rut, 0));
    // G
    let up = [(g, 1, &n + 1)];
    let nn = [(g, 1, n.clone())];
    let g1 = [(g, 1, one.clone())];
    k.edge(&up, (rub, 2), (rut, 1)); // up through R
    k.edge(&up, (rut, 3), (split_top, 0));
    k.edge(&nn, (split_top, 1), (bnt, 1)); // the n-bundle continues up
    k.edge(&nn, (bnt, 3), (bdt, 0)); // over the top cap, down
    k.edge(&nn, (bdt, 2), (rdt, 3));
    k.edge(&nn, (rdt, 1), (rdb, 2)); // down through R
    k.edge(&nn, (rdb, 0), (bdb, 0));
    k.edge(&nn, (bdb, 2), (bnb, 3)); // around the bottom cap, up
    k.edge(&nn, (bnb, 1), (xa, 0));
    k.edge(&nn, (xa, 2), (split_bot, 1));
    k.edge(&up, (split_bot, 0), (rub, 0));
    k.edge(&g1, (split_top, 2), (xa, 3)); // the peeled strand, around the left
    k.edge(&g1, (xa, 1), (split_bot, 2));
    // B
    let bb = [(b, 1, one.clone())];
    k.edge(&bb, (rbb, 2), (rbt, 1)); // up through R
    k.edge(&bb, (rbt, 3), (split_in, 2)); // the small loop, into the box
    k.edge(
        &[(b, 1, one.clone()), (b, 1, one.clone())],
        (split_in, 0),
        (bx, 0),
    );
    k.edge(
        &[(b, 1, one.clone()), (b, 1, one.clone())],
        (bx, 1),
        (split_out, 0),
    );
    k.edge(&bb, (split_out, 1), (rbb, 0)); // back round to R
    k.edge(&bb, (split_out, 2), (bdb, 3)); // the long way: west under/over G below R
    k.edge(&bb, (bdb, 1), (bnb, 0));
    k.edge(&bb, (bnb, 2), (bnt, 0)); // around the outside, east above R
    k.edge(&bb, (bnt, 2), (bdt, 1));
    k.edge(&bb, (bdt, 3), (split_in, 1));
    k.validate().unwrap();
    k
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
    assert_eq!(c.linking().writhe[0].display(&c.param_names()), "2*n*m + m + n");
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
