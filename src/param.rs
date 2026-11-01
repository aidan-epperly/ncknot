//! Polynomials in integer parameters, for families of diagrams.
//!
//! Widths of bundles, numbers of twists, framings and linking numbers of a
//! family are polynomials with rational coefficients in the family's
//! parameters (rational for counts like `n(n−1)/2`).  The type is ncalg's
//! [`MPoly`](ncalg::commutative::MPoly) over ℚ, with parsing
//! (`Poly::parse("n+1", &names)`), display, and evaluation at integers.

use crate::error::{KnotError, Result};
use ncalg::BigRational;
use num_traits::ToPrimitive;

/// A polynomial in the parameters.
pub type Poly = ncalg::commutative::MPoly<BigRational>;

/// Evaluation of a [`Poly`] to a machine integer.
pub trait EvalInt {
    /// The value at `values` (one per parameter), which must be an integer
    /// that fits an `i64`.
    fn eval_int(&self, values: &[i64]) -> Result<i64>;
}

impl EvalInt for Poly {
    fn eval_int(&self, values: &[i64]) -> Result<i64> {
        let v = self.eval_i64(values);
        if !v.is_integer() {
            return Err(KnotError::Invalid(format!("{v} is not an integer")));
        }
        v.to_integer()
            .to_i64()
            .ok_or_else(|| KnotError::Invalid(format!("{v} is too large")))
    }
}
