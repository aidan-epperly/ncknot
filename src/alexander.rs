//! Invariants of the Alexander module: the Alexander and Conway
//! polynomials, the determinant, Fox `p`-colorings and the homology of the
//! double branched cover.
//!
//! All of them come from the Alexander matrix of the Wirtinger presentation
//! ([`KnotGroup::alexander_matrix`]).  Its `(n−1)`-minors are computed over
//! `ℤ[t]` by evaluation and interpolation: each row is shifted to a
//! polynomial, the determinant is taken at `deg + 1` integer points with
//! ncalg's exact Bareiss determinant and the coefficients are recovered by an
//! exact Vandermonde solve.  Colorings use ncalg's `FqMat` ranks, branched
//! covers its Smith normal form.

use crate::bracket::IntLaurent;
use crate::diagram::Diagram;
use crate::error::{KnotError, Result};
use crate::group::KnotGroup;
use ncalg::ff::{FqMat, Gf};
use ncalg::linalg::exact;
use ncalg::zmod::IntMat;
use ncalg::BigRational;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

/// `det` of a square matrix over `ℤ[t^{±1}]`.
pub fn laurent_det(m: &[Vec<IntLaurent>]) -> IntLaurent {
    let n = m.len();
    if n == 0 {
        return IntLaurent::one();
    }
    assert!(m.iter().all(|r| r.len() == n), "square matrix");
    // shift each row into ℤ[t]
    let mut shift = 0i64;
    let mut deg = 0usize;
    let mut rows: Vec<Vec<IntLaurent>> = Vec::with_capacity(n);
    for r in m {
        let Some(lo) = r.iter().filter_map(|p| p.low_degree()).min() else {
            return IntLaurent::zero();
        };
        let hi = r.iter().filter_map(|p| p.high_degree()).max().unwrap();
        shift += lo;
        deg += (hi - lo) as usize;
        rows.push(r.iter().map(|p| p.shift(-lo)).collect());
    }
    // det at the points 0, 1, -1, 2, -2, ...
    let q = |k: i64| BigRational::from_integer(BigInt::from(k));
    let points: Vec<i64> = (0..=deg as i64)
        .map(|j| if j % 2 == 1 { (j + 1) / 2 } else { -j / 2 })
        .collect();
    let values: Vec<BigRational> = points
        .iter()
        .map(|&x| {
            let mx: Vec<Vec<BigRational>> = rows
                .iter()
                .map(|r| {
                    r.iter()
                        .map(|p| p.eval_with(&q(x), |c| BigRational::from_integer(c.clone())))
                        .collect()
                })
                .collect();
            exact::det(&mx)
        })
        .collect();
    let vandermonde: Vec<Vec<BigRational>> = points
        .iter()
        .map(|&x| (0..=deg as u32).map(|e| q(x).pow(e as i32)).collect())
        .collect();
    let coeffs = exact::solve(&vandermonde, &values).expect("distinct interpolation points");
    let ints = coeffs
        .into_iter()
        .map(|c| {
            assert!(c.is_integer(), "integer determinant");
            c.to_integer()
        })
        .collect();
    IntLaurent::from_coeffs(shift, ints)
}

/// The Alexander matrix at `t = −1` (the coloring matrix), as integers.
fn coloring_matrix(g: &KnotGroup) -> Vec<Vec<BigInt>> {
    g.alexander_matrix()
        .iter()
        .map(|r| r.iter().map(at_minus_one).collect())
        .collect()
}

/// `p(−1)`.
fn at_minus_one(p: &IntLaurent) -> BigInt {
    p.terms()
        .fold(BigInt::zero(), |s, (e, c)| if e.is_odd() { s - c } else { s + c })
}

/// The Alexander polynomial.  For a knot it is Conway-normalized:
/// `Δ(t^{−1}) = Δ(t)` and `Δ(1) = 1`.  For a link it is the one-variable
/// polynomial up to units, printed with lowest exponent 0 and a positive
/// lowest coefficient; split links (a component that never passes under
/// another, or crossingless loops) give 0.
pub fn alexander(d: &Diagram) -> IntLaurent {
    let n = d.crossings();
    if n == 0 {
        return if d.components() == 1 {
            IntLaurent::one()
        } else {
            IntLaurent::zero()
        };
    }
    let g = KnotGroup::wirtinger(d);
    let a = g.alexander_matrix();
    if g.presentation().generators != n {
        return IntLaurent::zero();
    }
    let minor: Vec<Vec<IntLaurent>> = a[..n - 1].iter().map(|r| r[1..].to_vec()).collect();
    let p = laurent_det(&minor);
    let (Some(lo), Some(hi)) = (p.low_degree(), p.high_degree()) else {
        return p;
    };
    if d.components() == 1 {
        assert!((hi - lo) % 2 == 0, "knot Alexander polynomials are symmetric");
        let p = p.shift(-(lo + hi) / 2);
        let at1 = p.terms().fold(BigInt::zero(), |s, (_, c)| s + c);
        assert!(at1.abs().is_one(), "Δ(1) = ±1 for a knot");
        if at1.is_negative() {
            p.negated()
        } else {
            p
        }
    } else {
        let p = p.shift(-lo);
        if p.coeff(0).is_negative() {
            p.negated()
        } else {
            p
        }
    }
}

/// The Conway polynomial `∇(z)` of a knot, with
/// `Δ(t) = ∇(t^{1/2} − t^{−1/2})`.
pub fn conway(d: &Diagram) -> Result<IntLaurent> {
    if d.components() != 1 {
        return Err(KnotError::Unsupported(
            "the Conway polynomial is computed for knots only".into(),
        ));
    }
    let delta = alexander(d);
    // Δ = a_0 + Σ_k a_k (t^k + t^{-k}); t^k + t^{-k} = P_k(u), u = t + t^{-1} = z² + 2
    let u = IntLaurent::from_terms([(2, BigInt::one()), (0, BigInt::from(2))]);
    let top = delta.high_degree().unwrap_or(0);
    let mut p_prev = IntLaurent::constant(BigInt::from(2)); // P_0
    let mut p_cur = u.clone(); // P_1
    let mut out = IntLaurent::constant(delta.coeff(0));
    for k in 1..=top {
        out = out.add_ref(&p_cur.scale(&delta.coeff(k)));
        let next = u.mul_ref(&p_cur).sub_ref(&p_prev);
        p_prev = std::mem::replace(&mut p_cur, next);
    }
    Ok(out)
}

/// The determinant `|Δ(−1)|`.
pub fn determinant(d: &Diagram) -> BigInt {
    at_minus_one(&alexander(d)).abs()
}

/// The number of Fox `p`-colorings (`p` prime), trivial ones included:
/// `p^{dim ker}` of the coloring matrix over `𝔽_p`.
pub fn colorings(d: &Diagram, p: u64) -> Result<BigInt> {
    let f = Gf::prime(p).map_err(KnotError::Invalid)?;
    let g = KnotGroup::wirtinger(d);
    let cols = g.presentation().generators;
    let m = coloring_matrix(&g);
    let rank = if m.is_empty() {
        0
    } else {
        let pb = BigInt::from(p);
        let rows: Vec<Vec<u64>> = m
            .iter()
            .map(|r| r.iter().map(|x| x.mod_floor(&pb).to_u64().unwrap()).collect())
            .collect();
        FqMat::from_rows(&f, &rows).rank()
    };
    Ok(BigInt::from(p).pow((cols - rank) as u32))
}

/// `H_1` of the double branched cover `Σ_2(L)`, as the orders of its cyclic
/// factors (`0` for a factor `ℤ`, factors of order 1 omitted): the cokernel
/// of the coloring matrix with one column deleted.  For a knot it is a finite
/// group of order [`determinant`].
pub fn double_branched_cover_homology(d: &Diagram) -> Vec<BigInt> {
    let g = KnotGroup::wirtinger(d);
    let cols = g.presentation().generators;
    if cols <= 1 {
        return Vec::new();
    }
    let m: Vec<Vec<BigInt>> = coloring_matrix(&g).into_iter().map(|r| r[1..].to_vec()).collect();
    let mut out = Vec::new();
    if !m.is_empty() {
        let s = IntMat::from_big_rows(&m).smith();
        out.extend(s.diag.iter().filter(|x| !x.is_one()).map(|x| x.abs()));
    }
    // generators beyond the number of relations are free
    out.extend((m.len()..cols - 1).map(|_| BigInt::zero()));
    out.sort_by(|a, b| match (a.is_zero(), b.is_zero()) {
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        _ => a.cmp(b),
    });
    out
}
