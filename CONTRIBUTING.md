# Contributing to ncknot

The rules are [ncalg's](https://github.com/aidan-epperly/ncalg/blob/main/CONTRIBUTING.md); the short version:

- **Library only.**  No binaries and no result files.  A
  project that uses ncknot is its own crate with ncknot as a dependency.
- **Safe and documented.**  `#![forbid(unsafe_code)]`, `#![warn(missing_docs)]`
  (every public item has a doc comment that says what it computes), `cargo fmt`,
  `cargo clippy --all-targets` clean, `cargo doc --no-deps` without warnings.
- **One implementation per algorithm.**  Exact linear algebra (determinants,
  solves, inertia, ranks over `𝔽_p`, Smith forms, sparse invariant factors)
  is ncalg's; add a missing general routine there, with a test, rather than a
  private copy here.  Inside ncknot, a helper that two modules need lives in
  one place (`laurent`, `diagram::faces`, `group::arcs`).
- **Exact.**  Every invariant is computed exactly (`BigInt`, `BigRational`,
  `𝔽_p`); no floating point.
- **Conventions are tested.**  A new invariant is checked against published
  values for named knots *and* for invariance under random Markov moves
  (`tests/invariants.rs`); a sign or chirality convention is stated in the
  module docs and pinned by a test.

## Issues

One file per issue in `issues/`, as in ncalg:

```bash
python scripts/issues.py new bug "<what goes wrong>" --area khovanov --severity 1 --reporter "<who>"
python scripts/issues.py list
python scripts/issues.py close 3 --commit abc1234 --note "..."
```

Areas: `diagram braid laurent bracket alexander signature group khovanov
table examples docs build`.  Severity 1 is a wrong result without an error,
2 a panic or error on valid input, 3 performance or usability.

## Layout

```
src/laurent.rs     Laurent<C>: Laurent polynomials over an ncalg Coeff
src/diagram.rs     Diagram: PD codes, orientation, planarity, faces, mirror, linking
src/braid.rs       Braid: words, Markov stabilization, torus links, closure
src/table.rs       named knots and links (braid words)
src/bracket.rs     Kauffman bracket (tangle contraction), Jones
src/group.rs       Wirtinger presentation, Fox calculus, homomorphisms to S_m
src/alexander.rs   Alexander, Conway, determinant, colorings, double branched cover
src/signature.rs   Goeritz matrix, Gordon–Litherland signature
src/khovanov.rs    Khovanov homology over F_p and Z
tests/             invariants.rs
examples/          table, bench
docs/DESIGN.md     design notes
```
