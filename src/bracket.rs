//! The Kauffman bracket and the Jones polynomial.
//!
//! The bracket is contracted crossing by crossing (a planar-algebra / tangle
//! contraction): after some crossings have been absorbed, the partial state
//! sum is a map from *boundary matchings* (how the open edges are joined by
//! arcs inside the absorbed region) to Laurent polynomials in `A`.  Adding a
//! crossing smooths it both ways, joins the new arcs to the old ones, and
//! multiplies by `δ = −A² − A^{−2}` for every loop that closes.  Crossings are
//! taken greedily, the next one sharing the most edges with the boundary, so
//! the cost is exponential in the boundary width, not in the crossing number.
//!
//! Smoothing convention for `X[a, b, c, d]`: the **A-smoothing** joins
//! `a–b` and `c–d`, the **B-smoothing** `a–d` and `b–c`.  With it,
//! `⟨D⟩ = Σ_s A^{#A − #B} δ^{|s| − 1}` and the Jones polynomial is
//! `V(t) = (−A³)^{−w(D)} ⟨D⟩` at `A = t^{−1/4}`; the positive (right-handed)
//! trefoil has `V = t + t³ − t⁴`.

use crate::diagram::Diagram;
use crate::laurent::Laurent;
use num_bigint::BigInt;
use std::collections::HashMap;

/// Laurent polynomials with integer coefficients.
pub type IntLaurent = Laurent<BigInt>;

/// `δ = −A² − A^{−2}`.
pub fn delta() -> IntLaurent {
    IntLaurent::from_terms([(2, BigInt::from(-1)), (-2, BigInt::from(-1))])
}

/// The order in which [`state_sum`] absorbs crossings: greedily, the
/// crossing with the most edges already on the boundary.
pub(crate) fn contraction_order(d: &Diagram) -> Vec<usize> {
    let n = d.crossings();
    let mut done = vec![false; n];
    let mut ends_done = vec![0u8; d.edges()];
    let mut order = Vec::with_capacity(n);
    for _ in 0..n {
        let c = (0..n)
            .filter(|&c| !done[c])
            .max_by_key(|&c| {
                let shared = d.crossing(c).iter().filter(|&&e| ends_done[e] > 0).count();
                (shared, std::cmp::Reverse(c))
            })
            .unwrap();
        done[c] = true;
        for e in d.crossing(c) {
            ends_done[e] += 1;
        }
        order.push(c);
    }
    order
}

/// The unnormalized bracket `Σ_s A^{#A − #B} δ^{|s|}` (one more factor of
/// `δ` than `⟨D⟩`; the crossingless unknot gives `δ`).
pub fn bracket_unnormalized(d: &Diagram) -> IntLaurent {
    let mut delta_pow = vec![IntLaurent::one()];
    let loop_factor = |k: usize, delta_pow: &mut Vec<IntLaurent>| {
        while delta_pow.len() <= k {
            let p = delta_pow.last().unwrap().mul_ref(&delta());
            delta_pow.push(p);
        }
        delta_pow[k].clone()
    };
    // key: sorted pairs of open edges joined inside the absorbed region
    let mut states: HashMap<Vec<(u32, u32)>, IntLaurent> = HashMap::new();
    states.insert(Vec::new(), IntLaurent::one());
    let mut ends_done = vec![0u8; d.edges()];
    for c in contraction_order(d) {
        let x = d.crossing(c).map(|e| e as u32);
        for &e in &x {
            ends_done[e as usize] += 1;
        }
        let smoothings = [
            ([(x[0], x[1]), (x[2], x[3])], 1i64),
            ([(x[0], x[3]), (x[1], x[2])], -1i64),
        ];
        let mut next: HashMap<Vec<(u32, u32)>, IntLaurent> = HashMap::new();
        for (key, poly) in &states {
            for (arcs, a_exp) in &smoothings {
                let (new_key, loops) = join(key, arcs, &ends_done);
                let term = poly.shift(*a_exp).mul_ref(&loop_factor(loops, &mut delta_pow));
                let slot = next.entry(new_key).or_insert_with(IntLaurent::zero);
                *slot = slot.add_ref(&term);
            }
        }
        next.retain(|_, p| !p.is_zero());
        states = next;
    }
    let total = states.remove(&Vec::new()).unwrap_or_else(IntLaurent::zero);
    total.mul_ref(&loop_factor(d.free_loops(), &mut delta_pow))
}

/// Joins the boundary matching `key` with two new arcs; edges with both ends
/// absorbed (`ends_done == 2`) are interior.  Returns the new matching of the
/// open edges and the number of loops closed.
fn join(key: &[(u32, u32)], arcs: &[(u32, u32); 2], ends_done: &[u8]) -> (Vec<(u32, u32)>, usize) {
    // adjacency of the edge-nodes through arcs (each node has degree ≤ 2)
    let mut adj: HashMap<u32, [u32; 2]> = HashMap::new();
    let mut deg: HashMap<u32, u8> = HashMap::new();
    let mut link = |u: u32, v: u32| {
        for (p, q) in [(u, v), (v, u)] {
            let k = deg.entry(p).or_insert(0);
            adj.entry(p).or_insert([u32::MAX; 2])[*k as usize] = q;
            *k += 1;
        }
    };
    for &(u, v) in key.iter().chain(arcs.iter()) {
        link(u, v);
    }
    let interior = |e: u32| ends_done[e as usize] == 2;
    let mut visited: HashMap<u32, bool> = HashMap::new();
    let mut pairs = Vec::new();
    // walk from every open node to its partner
    let mut nodes: Vec<u32> = adj.keys().copied().collect();
    nodes.sort_unstable();
    for &s in &nodes {
        if interior(s) || visited.contains_key(&s) {
            continue;
        }
        visited.insert(s, true);
        let (mut prev, mut cur) = (s, adj[&s][0]);
        while interior(cur) {
            visited.insert(cur, true);
            let nxt = if adj[&cur][0] == prev {
                adj[&cur][1]
            } else {
                adj[&cur][0]
            };
            prev = cur;
            cur = nxt;
        }
        visited.insert(cur, true);
        pairs.push((s.min(cur), s.max(cur)));
    }
    // what is left are closed loops of interior nodes
    let mut loops = 0;
    for &s in &nodes {
        if visited.contains_key(&s) {
            continue;
        }
        loops += 1;
        let (mut prev, mut cur) = (s, adj[&s][0]);
        visited.insert(s, true);
        while cur != s {
            visited.insert(cur, true);
            let nxt = if adj[&cur][0] == prev {
                adj[&cur][1]
            } else {
                adj[&cur][0]
            };
            prev = cur;
            cur = nxt;
        }
    }
    pairs.sort_unstable();
    (pairs, loops)
}

/// The Kauffman bracket `⟨D⟩` in `A`, normalized so that the crossingless
/// unknot has `⟨O⟩ = 1`.
pub fn kauffman_bracket(d: &Diagram) -> IntLaurent {
    bracket_unnormalized(d)
        .div_exact(&delta())
        .expect("the unnormalized bracket is divisible by δ")
}

/// The Jones polynomial, stored in the variable `t^{1/2}`.
#[derive(Clone, Debug, PartialEq)]
pub struct Jones {
    half: IntLaurent,
}

impl Jones {
    /// `V` as a Laurent polynomial in `s = t^{1/2}`.
    pub fn in_sqrt_t(&self) -> &IntLaurent {
        &self.half
    }

    /// `V` as a Laurent polynomial in `t`, when only integer powers occur
    /// (always for knots, and for links with an odd number of components).
    pub fn in_t(&self) -> Option<IntLaurent> {
        self.half
            .exponents_divisible_by(2)
            .then(|| self.half.compress_exponents(2))
    }

    /// `V` from a Laurent polynomial in `s = t^{1/2}`.
    pub fn from_sqrt_t(half: IntLaurent) -> Self {
        Jones { half }
    }
}

impl std::fmt::Display for Jones {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.half.display_as("t", 2))
    }
}

/// The Jones polynomial `V(t) = (−A³)^{−w} ⟨D⟩` at `A = t^{−1/4}`.
pub fn jones(d: &Diagram) -> Jones {
    let w = d.writhe();
    let mut f = kauffman_bracket(d).shift(-3 * w);
    if w.rem_euclid(2) == 1 {
        f = f.negated();
    }
    // A^k = t^{-k/4} = s^{-k/2}
    assert!(
        f.exponents_divisible_by(2),
        "bracket exponents of a link are even"
    );
    Jones {
        half: f.compress_exponents(2).substitute_power(-1),
    }
}
