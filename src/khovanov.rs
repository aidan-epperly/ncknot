//! Khovanov homology from the cube of resolutions (Bar-Natan's conventions).
//!
//! A state `s ∈ {0,1}^n` smooths crossing `c` by its 0-smoothing (the
//! A-smoothing, `a–b` and `c–d` of `X[a, b, c, d]`) or its 1-smoothing.  The
//! chain group at `s` is `V^{⊗ circles}`, `V = ⟨v₊, v₋⟩`; a generator has
//! homological degree `h = |s| − n₋` and quantum degree
//! `q = #v₊ − #v₋ + |s| + n₊ − 2n₋`.  Edges `s → s + e_c` merge two circles
//! (`m`) or split one (`Δ`), with sign `(−1)^{#1s of s before c}`.  The
//! differential preserves `q`, so the complex is split into its `q`-graded
//! pieces.  Ranks over `𝔽_p` use ncalg's sparse Howell form
//! ([`ncalg::zmod::SparseHowell`]); over `ℤ` the sparse Smith form
//! ([`ncalg::zmod::sparse_invariant_factors`]) gives ranks and torsion.
//!
//! `Σ (−1)^h q^j rank Kh^{h,j} = (q + q^{−1}) V(t)|_{t^{1/2} = −q}`.
//! The cost is `2^n` states times `2^{circles}` generators: comfortable up
//! to about 14 crossings over `𝔽_p` or `ℤ`.

use crate::bracket::IntLaurent;
use crate::diagram::Diagram;
use crate::error::{KnotError, Result};
use ncalg::zmod::{sparse_invariant_factors, SparseHowell};
use num_bigint::BigInt;
use num_traits::One;
use std::collections::{BTreeMap, HashMap};

/// The coefficient ring of the homology.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ring {
    /// The prime field `𝔽_p`.
    Fp(u64),
    /// The integers: ranks and torsion.
    Integers,
}

/// Khovanov homology, bigraded by `(h, q)`.
#[derive(Clone, Debug, PartialEq)]
pub struct Khovanov {
    /// The coefficient ring.
    pub ring: Ring,
    /// Ranks (dimensions over `𝔽_p`, free ranks over `ℤ`), zeros omitted.
    pub ranks: BTreeMap<(i64, i64), usize>,
    /// Over `ℤ`: the torsion invariant factors of each group.
    pub torsion: BTreeMap<(i64, i64), Vec<BigInt>>,
}

impl Khovanov {
    /// The rank in bidegree `(h, q)`.
    pub fn rank(&self, h: i64, q: i64) -> usize {
        self.ranks.get(&(h, q)).copied().unwrap_or(0)
    }

    /// The sum of the ranks.
    pub fn total_rank(&self) -> usize {
        self.ranks.values().sum()
    }

    /// The graded Euler characteristic `Σ (−1)^h q^j rank Kh^{h,j}`, a
    /// Laurent polynomial in `q`.
    pub fn euler_characteristic(&self) -> IntLaurent {
        IntLaurent::from_terms(self.ranks.iter().map(|(&(h, q), &r)| {
            let r = BigInt::from(r);
            (q, if h % 2 == 0 { r } else { -r })
        }))
    }

    /// The Poincaré polynomial as text, `Σ rank · t^h q^j`, with torsion
    /// summands appended as `ℤ/k·t^h q^j`.
    pub fn poincare(&self) -> String {
        let mono = |h: i64, q: i64| {
            let t = match h {
                0 => None,
                1 => Some("t".to_string()),
                _ => Some(format!("t^{h}")),
            };
            let q = match q {
                0 => None,
                1 => Some("q".to_string()),
                _ => Some(format!("q^{q}")),
            };
            match (t, q) {
                (None, None) => "1".to_string(),
                (t, q) => [t, q].into_iter().flatten().collect::<Vec<_>>().join(" "),
            }
        };
        let mut parts: Vec<String> = self
            .ranks
            .iter()
            .map(|(&(h, q), &r)| {
                if r == 1 {
                    mono(h, q)
                } else {
                    format!("{r} {}", mono(h, q))
                }
            })
            .collect();
        for (&(h, q), ts) in &self.torsion {
            for k in ts {
                parts.push(format!("Z/{k} {}", mono(h, q)));
            }
        }
        if parts.is_empty() {
            "0".into()
        } else {
            parts.join(" + ")
        }
    }
}

/// The circles of every state: `circ[s][e]` is the circle of edge `e`,
/// circles numbered by their smallest edge.
fn circles(d: &Diagram) -> Vec<(Vec<u8>, u8)> {
    let n = d.crossings();
    let m = d.edges();
    (0..1u32 << n)
        .map(|s| {
            let mut parent: Vec<usize> = (0..m).collect();
            fn find(p: &mut [usize], mut x: usize) -> usize {
                while p[x] != x {
                    p[x] = p[p[x]];
                    x = p[x];
                }
                x
            }
            let mut union = |a: usize, b: usize| {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                parent[ra] = rb;
            };
            for c in 0..n {
                let x = d.crossing(c);
                if s >> c & 1 == 0 {
                    union(x[0], x[1]);
                    union(x[2], x[3]);
                } else {
                    union(x[0], x[3]);
                    union(x[1], x[2]);
                }
            }
            let mut id = vec![u8::MAX; m];
            let mut circ = vec![0u8; m];
            let mut k = 0u8;
            for e in 0..m {
                let r = find(&mut parent, e);
                if id[r] == u8::MAX {
                    id[r] = k;
                    k += 1;
                }
                circ[e] = id[r];
            }
            (circ, k)
        })
        .collect()
}

/// Khovanov homology of the diagram over `ring`.  Declines more than 24
/// crossings or more than 31 circles in a state.
pub fn khovanov(d: &Diagram, ring: Ring) -> Result<Khovanov> {
    let n = d.crossings();
    if n > 24 {
        return Err(KnotError::Unsupported(format!(
            "{n} crossings: the cube of resolutions has 2^{n} states"
        )));
    }
    if let Ring::Fp(p) = ring {
        if !ncalg::arith::int::is_prime(p) {
            return Err(KnotError::Invalid(format!("{p} is not prime")));
        }
    }
    let free = d.free_loops();
    let circ = circles(d);
    if circ.iter().any(|(_, k)| *k as usize + free > 31) {
        return Err(KnotError::Unsupported("more than 31 circles in a state".into()));
    }
    let (np, nm) = d.sign_counts();
    let (np, nm) = (np as i64, nm as i64);
    // generators (state, bits): bit j set = v₋ on circle j (free loops after the others)
    let mut groups: HashMap<(i64, i64), Vec<(u32, u32)>> = HashMap::new();
    for (s, (_, k)) in circ.iter().enumerate() {
        let r = (s as u32).count_ones() as i64;
        let k = *k as u32 + free as u32;
        for b in 0..1u32 << k {
            let q = k as i64 - 2 * b.count_ones() as i64 + r + np - 2 * nm;
            groups.entry((r - nm, q)).or_default().push((s as u32, b));
        }
    }
    let index: HashMap<(u32, u32), usize> = groups
        .values()
        .flat_map(|g| g.iter().enumerate().map(|(i, &x)| (x, i)))
        .collect();
    let image = |s: u32, b: u32| -> Vec<((u32, u32), i64)> {
        let (cs, ks) = (&circ[s as usize].0, circ[s as usize].1 as u32);
        let mut out = Vec::new();
        for c in 0..n {
            if s >> c & 1 == 1 {
                continue;
            }
            let t = s | 1 << c;
            let (ct, kt) = (&circ[t as usize].0, circ[t as usize].1 as u32);
            let sign = if (s & ((1 << c) - 1)).count_ones() % 2 == 0 {
                1
            } else {
                -1
            };
            let x = d.crossing(c);
            let (c1, c2) = (cs[x[0]] as u32, cs[x[2]] as u32);
            let (d1, d2) = (ct[x[0]] as u32, ct[x[1]] as u32);
            // the untouched circles, by their smallest edge; then the free loops
            let mut base = 0u32;
            let mut seen = 0u64;
            for (e, &j) in cs.iter().enumerate() {
                let j = j as u32;
                if j == c1 || j == c2 || seen >> j & 1 == 1 {
                    continue;
                }
                seen |= 1 << j;
                base |= (b >> j & 1) << ct[e];
            }
            for f in 0..free as u32 {
                base |= (b >> (ks + f) & 1) << (kt + f);
            }
            let bit = |j: u32| b >> j & 1;
            if c1 != c2 {
                // merge: v₊v₊ = v₊, v₊v₋ = v₋, v₋v₋ = 0
                if bit(c1) & bit(c2) == 0 {
                    out.push(((t, base | (bit(c1) | bit(c2)) << d1), sign));
                }
            } else if bit(c1) == 0 {
                // split: v₊ ↦ v₊⊗v₋ + v₋⊗v₊
                out.push(((t, base | 1 << d2), sign));
                out.push(((t, base | 1 << d1), sign));
            } else {
                // v₋ ↦ v₋⊗v₋
                out.push(((t, base | 1 << d1 | 1 << d2), sign));
            }
        }
        out
    };
    // rank and torsion of d: C^{h,q} → C^{h+1,q}
    let mut rank: HashMap<(i64, i64), usize> = HashMap::new();
    let mut torsion: BTreeMap<(i64, i64), Vec<BigInt>> = BTreeMap::new();
    for (&(h, q), src) in &groups {
        let Some(dst) = groups.get(&(h + 1, q)) else {
            continue;
        };
        let rows: Vec<Vec<(usize, i64)>> = src
            .iter()
            .map(|&(s, b)| image(s, b).into_iter().map(|(g, v)| (index[&g], v)).collect())
            .collect();
        let r = match ring {
            Ring::Fp(p) => {
                let mut howell = SparseHowell::new(p, dst.len());
                for row in &rows {
                    let v: Vec<(usize, u64)> = row
                        .iter()
                        .map(|&(j, x)| (j, x.rem_euclid(p as i64) as u64))
                        .collect();
                    howell.insert(&v);
                }
                howell.len()
            }
            Ring::Integers => {
                let factors = sparse_invariant_factors(dst.len(), &rows);
                let tors: Vec<BigInt> = factors.iter().filter(|x| !x.is_one()).cloned().collect();
                if !tors.is_empty() {
                    torsion.insert((h + 1, q), tors);
                }
                factors.len()
            }
        };
        rank.insert((h, q), r);
    }
    let mut ranks = BTreeMap::new();
    for (&(h, q), src) in &groups {
        let out = rank.get(&(h, q)).copied().unwrap_or(0);
        let inc = rank.get(&(h - 1, q)).copied().unwrap_or(0);
        let r = src.len() - out - inc;
        if r > 0 {
            ranks.insert((h, q), r);
        }
    }
    Ok(Khovanov { ring, ranks, torsion })
}
