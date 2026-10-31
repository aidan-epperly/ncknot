# ncknot — computational knot theory on ncalg

A library in safe Rust (`#![forbid(unsafe_code)]`) for computing invariants
of knots and links exactly, built on the algebra of
[ncalg](https://github.com/aidan-epperly/ncalg): exact determinants and
inertia over ℚ, ranks over finite fields, and Smith normal forms over ℤ.  It
is a library only; add it to a crate with
`ncknot = { git = "https://github.com/aidan-epperly/ncknot" }` (it fetches
ncalg the same way).  A browser
front end is [ncknot-web](https://github.com/aidan-epperly/ncknot-web).

## Contents

- [Quick start](#quick-start)
- [Diagrams](#diagrams)
- [Invariants](#invariants)
- [Conventions](#conventions)
- [What ncalg provides](#what-ncalg-provides)
- [Examples, tests and timings](#examples-tests-and-timings)
- [Roadmap](#roadmap)
- [Bugs and improvements](#bugs-and-improvements)

## Quick start

```rust
use ncknot::{table, Braid, Diagram, Ring};

let trefoil = Braid::new(2, &[1, 1, 1]).unwrap().closure();   // right-handed
assert_eq!(trefoil.jones().to_string(), "t + t^3 - t^4");
assert_eq!(trefoil.alexander().to_string(), "t^-1 - 1 + t");
assert_eq!(trefoil.signature().unwrap(), -2);
assert_eq!(trefoil.colorings(3).unwrap(), 9.into());

// KnotTheory's PD code for 3_1 (the left-handed trefoil)
let lh = Diagram::parse_pd("PD[X[1,4,2,5], X[3,6,4,1], X[5,2,6,3]]").unwrap();
assert_eq!(lh.jones().to_string(), "-t^-4 + t^-3 + t^-1");

let kh = table::diagram("3_1").unwrap().khovanov(Ring::Integers).unwrap();
println!("{}", kh.poincare());  // q + q^3 + t^2 q^5 + t^3 q^9 + Z/2 t^3 q^7
```

## Diagrams

`Diagram` is an oriented link diagram stored as a planar diagram (PD) code.

| constructor | input |
|---|---|
| `Diagram::from_pd(&[[i64; 4]])`, `Diagram::parse_pd(&str)` | PD codes in KnotTheory's convention (`X[...]` lists, an optional `PD[...]` wrapper) |
| `Braid::new(n, &word).closure()` | braid words, `i` for `σ_i` and `-i` for `σ_i^{-1}` |
| `Braid::torus(p, q)` | torus links `T(p, q)` |
| `table::diagram(name)` | `0_1`, the Rolfsen knots `3_1` to `7_7`, `8_19`, `hopf`, `T(2,4)`, `borromean` |
| `Diagram::unknot()`, `Diagram::unlink(k)` | crossingless unlinks |

Construction checks the code.  Every label must occur twice, the orientations
must be consistent, and the code must be planar (Euler's count of faces), so
virtual codes are rejected.  Edges are renumbered consecutively along the
components.

Operations: `mirror`, `reverse`, `disjoint_union`, `writhe`, `sign_counts`,
`components`, `linking_matrix`, `faces`, `pd_code`.  `Braid` adds `compose`,
`inverse`, `stabilize` (Markov) and `permutation`.

## Invariants

| method | what | how |
|---|---|---|
| `kauffman_bracket` | `⟨D⟩ ∈ ℤ[A^{±1}]`, `⟨O⟩ = 1` | tangle contraction: a state sum over boundary matchings, crossings absorbed greedily |
| `jones` | `V(t)`, in powers of `t^{1/2}` | `(−A³)^{−w} ⟨D⟩` at `A = t^{−1/4}` |
| `alexander` | `Δ(t)`; Conway-normalized for knots | an `(n−1)`-minor of the Fox Jacobian of the Wirtinger presentation, by evaluation and exact interpolation |
| `conway` | `∇(z)` of a knot | from the symmetric `Δ` |
| `determinant` | `\|Δ(−1)\|` | |
| `colorings(p)` | the number of Fox `p`-colorings | `p^{dim ker}` of the coloring matrix over `𝔽_p` |
| `double_branched_cover_homology` | `H_1(Σ_2)`, as its cyclic factors | Smith form of the reduced coloring matrix |
| `signature` | `σ`, negative for positive knots | Goeritz matrix and the Gordon–Litherland correction; exact inertia over ℚ |
| `group` | the Wirtinger presentation (`KnotGroup`) | `presentation`, `alexander_matrix`, `Presentation::fox_jacobian` |
| `group().count_homomorphisms_to_symmetric(m, transitive)` | the number of homomorphisms to `S_m`, `m ≤ 10` | backtracking with relation propagation and conjugacy-class symmetry |
| `khovanov(Ring::Fp(p))`, `khovanov(Ring::Integers)` | `Kh^{h,q}`: ranks, plus torsion over ℤ | the cube of resolutions, split by `q`; sparse ranks and invariant factors |

`Khovanov` has `rank(h, q)`, `total_rank`, `euler_characteristic` (equal to
`(q + q^{−1}) V` at `t^{1/2} = −q`) and `poincare` (as text).
`Laurent<C>` is the polynomial type throughout: Laurent polynomials over any
ncalg `Coeff`, with exact division and evaluation.

## Conventions

- **PD codes.** `X[a, b, c, d]` lists edges counterclockwise from the incoming
  under-strand.  The crossing is positive if the over-strand enters at `d`.
- **Braids.** Strands run upward, and in `σ_i` strand `i` passes over strand
  `i + 1`, so `σ_i` is a positive crossing.  `(σ_1)^3` is the right-handed trefoil.
- **Kauffman bracket.** The A-smoothing of `X[a, b, c, d]` joins `a–b` and `c–d`.
- **Jones.** The right-handed trefoil has `V = t + t³ − t⁴`.  The positive
  Hopf link has `V = −t^{1/2} − t^{5/2}`.
- **Signature.** Positive knots have negative signature: `σ(3_1) = −2` for
  the right-handed trefoil.
- **Khovanov.** Bar-Natan's gradings.  The right-handed trefoil has
  `q + q³ + t²q⁵ + t³q⁹` and `ℤ/2` in `(3, 7)`.
- **Links.** The Alexander polynomial of a link is returned up to units
  (lowest exponent 0, positive lowest coefficient).  It is 0 for split links.
- **Table.** Each table entry has the chirality of its KnotInfo braid word.
  The signatures in the tests are those of that chirality.

## What ncalg provides

| ncknot step | ncalg |
|---|---|
| Alexander minors at interpolation points; the Vandermonde solve | `linalg::exact::{det, solve}` over `BigRational` |
| signature from the Goeritz matrix | `linalg::exact::inertia` (added for ncknot) |
| `p`-colorings | `ff::FqMat::rank` |
| double branched cover | `zmod::IntMat::smith` |
| Khovanov over `𝔽_p` | `zmod::SparseHowell` |
| Khovanov over ℤ | `zmod::sparse_invariant_factors` (added for ncknot) |
| primality, polynomial ring traits | `arith::int`, `Coeff` / `Field` |

## Examples, tests and timings

- `cargo run --release --example table` prints the invariants of every
  table entry.
- `cargo run --release --example bench` times each invariant.
- `cargo test --release`: the table's Alexander polynomials, determinants and
  signatures; Jones polynomials; KnotTheory PD codes; Conway polynomials and
  linking numbers; colorings and `S_m` counts; integral Khovanov homology of
  the trefoil.  Further checks:
  - every invariant is unchanged by random Markov moves (conjugation,
    stabilization, inserting `σσ^{−1}`);
  - the Khovanov Euler characteristic equals the Jones polynomial over
    `𝔽_2`, `𝔽_3` and ℤ;
  - mirror images, reversal and split unions behave as they should.

Release timings, one core, from `examples/bench.rs`:

| input | Jones | Alexander | signature | `S_5` reps | Kh over `𝔽_2` | Kh over ℤ |
|---|---|---|---|---|---|---|
| `T(3,7)`, 14 crossings | < 1 ms | 7 ms | < 1 ms | 1 ms | | |
| `T(5,6)`, 24 crossings | 2 ms | 56 ms | < 1 ms | 0.3 s | | |
| random 5-braid, 40 crossings | 2 ms | 0.55 s | < 1 ms | | | |
| 4-braid, 10 crossings | | | | | 25 ms | 31 ms |
| 4-braid, 12 crossings | | | | | 0.25 s | 0.44 s |

Khovanov homology costs `2^n` states times `2^{circles}` generators.
Comfortable sizes are up to about 14 crossings.

## Roadmap

- HOMFLY-PT and Kauffman two-variable polynomials, through the Hecke /
  BMW algebras.  The quotient algebras are ncalg's, but they need
  rational-function coefficients.
- Reidemeister simplification of diagrams.  Faster Khovanov through Bar-Natan's
  local (tangle) algorithm.
- Seifert matrices, Tristram–Levine signatures, and the multivariable
  Alexander polynomial.
- DT and Gauss codes as input; a larger knot table.

## Bugs and improvements

File one issue per problem in [`issues/`](issues/README.md) with
`python scripts/issues.py new ...`; see [CONTRIBUTING.md](CONTRIBUTING.md).
A missing feature or bug in ncalg goes in ncalg's tracker.

## License

Copyright (C) 2026 Aidan Epperly.  Licensed under the GNU General Public
License, version 3 or (at your option) any later version; see [LICENSE](LICENSE).
