//! # ncknot — computational knot theory on ncalg
//!
//! Oriented link diagrams ([`Diagram`], from PD codes or braid closures) and
//! their invariants, computed exactly:
//!
//! * [`bracket`] — the Kauffman bracket and the Jones polynomial, by tangle
//!   contraction;
//! * [`alexander`] — Alexander and Conway polynomials, determinant, Fox
//!   `p`-colorings, `H_1` of the double branched cover;
//! * [`signature`] — Goeritz matrix and signature (Gordon–Litherland);
//! * [`group`] — Wirtinger presentations, Fox calculus, counting
//!   homomorphisms to symmetric groups;
//! * [`khovanov`] — Khovanov homology over `𝔽_p` and `ℤ`.
//! * [`kirby`] — Kirby diagrams of infinite families: bundles of `n`
//!   strands and twist boxes with `m` half-twists, with linking numbers and
//!   framings as polynomials in the parameters ([`param`]), 1-/2-handle
//!   cancellation done on whole bundles, and expansion of any member to an
//!   ordinary diagram.
//!
//! The exact algebra is ncalg's: Bareiss determinants and the inertia of
//! symmetric matrices over ℚ (`linalg::exact`), ranks over `𝔽_p` (`ff`,
//! `zmod::SparseHowell`) and Smith normal forms over ℤ (`zmod::IntMat`).
//!
//! ```
//! use ncknot::{table, Braid};
//! let trefoil = Braid::new(2, &[1, 1, 1]).unwrap().closure();
//! assert_eq!(trefoil.jones().to_string(), "t + t^3 - t^4");
//! assert_eq!(trefoil.alexander().to_string(), "t^-1 - 1 + t");
//! assert_eq!(trefoil.signature().unwrap(), -2);
//! let fig8 = table::diagram("4_1").unwrap();
//! assert_eq!(fig8.determinant(), 5.into());
//! ```

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]
#![warn(missing_docs)]
#![allow(clippy::needless_range_loop)]

pub mod alexander;
pub mod bracket;
pub mod braid;
pub mod diagram;
pub mod error;
pub mod group;
pub mod khovanov;
pub mod kirby;
pub mod laurent;
pub mod param;
pub mod signature;
pub mod table;

pub use bracket::{IntLaurent, Jones};
pub use braid::Braid;
pub use diagram::Diagram;
pub use error::{KnotError, Result};
pub use group::{KnotGroup, Presentation};
pub use khovanov::{Khovanov, Ring};
pub use laurent::Laurent;
pub use num_bigint::BigInt;

/// The invariants as methods.
impl Diagram {
    /// The Kauffman bracket in `A` ([`bracket::kauffman_bracket`]).
    pub fn kauffman_bracket(&self) -> IntLaurent {
        bracket::kauffman_bracket(self)
    }
    /// The Jones polynomial ([`bracket::jones`]).
    pub fn jones(&self) -> Jones {
        bracket::jones(self)
    }
    /// The Alexander polynomial ([`alexander::alexander`]).
    pub fn alexander(&self) -> IntLaurent {
        alexander::alexander(self)
    }
    /// The Conway polynomial of a knot ([`alexander::conway`]).
    pub fn conway(&self) -> Result<IntLaurent> {
        alexander::conway(self)
    }
    /// The determinant ([`alexander::determinant`]).
    pub fn determinant(&self) -> BigInt {
        alexander::determinant(self)
    }
    /// The number of Fox `p`-colorings ([`alexander::colorings`]).
    pub fn colorings(&self, p: u64) -> Result<BigInt> {
        alexander::colorings(self, p)
    }
    /// `H_1` of the double branched cover ([`alexander::double_branched_cover_homology`]).
    pub fn double_branched_cover_homology(&self) -> Vec<BigInt> {
        alexander::double_branched_cover_homology(self)
    }
    /// The signature ([`signature::signature`]).
    pub fn signature(&self) -> Result<i64> {
        signature::signature(self)
    }
    /// The Wirtinger presentation of the link group ([`KnotGroup::wirtinger`]).
    pub fn group(&self) -> KnotGroup {
        KnotGroup::wirtinger(self)
    }
    /// Khovanov homology over `ring` ([`khovanov::khovanov`]).
    pub fn khovanov(&self, ring: Ring) -> Result<Khovanov> {
        khovanov::khovanov(self, ring)
    }
}
