# ncknot design notes

## Representation

One diagram type, the PD code. It is the input format of KnotTheory, SnapPy
and KnotInfo, and every invariant here reads it directly:
- the bracket and Khovanov homology need the smoothings of each crossing;
- the group needs the over-arcs;
- the signature needs the faces.

A braid is a constructor: `Braid::closure` writes the PD code.

Internally, edges are renumbered `0..2n` along the components, in the
direction of travel. Each edge stores its tail and head `(crossing, slot)`.
Orientation is read from the code: slot 0 is incoming and slot 2 outgoing,
and propagation along edges and strands orients the over-strands. A
component that only passes over is oriented arbitrarily.

Planarity is checked by Euler's formula, with faces traced by the rule "turn
left at every crossing". Each projection piece must have `crossings + 2`
faces. Without this check, a virtual code would give a bracket but a
meaningless signature.

Crossingless loops (`free_loops`) are kept as a count, so that unlinks and
braids with idle strands are diagrams. Each invariant handles them
explicitly:
- the bracket multiplies by `δ` per loop;
- Khovanov homology tensors with `V`;
- the group gains a free generator;
- the Alexander polynomial of a split link is 0.

## Algorithms and their costs

- **Kauffman bracket.** Tangle contraction. The partial state sum is a map
  from matchings of the open boundary edges to `ℤ[A^{±1}]`. Crossings are
  absorbed in a greedy order: the next one is the crossing sharing the most
  edges with the boundary. The cost is exponential in the boundary width,
  not in `n`; a 24-crossing torus knot takes 2 ms. The same order heuristic
  would serve a future Bar-Natan local Khovanov algorithm.
- **Alexander polynomial.** The `(n−1)`-minor of the abelianised Fox
  Jacobian lives in `ℤ[t]`, after shifting each row. It is evaluated at
  `deg + 1` small integers with ncalg's Bareiss determinant and recovered by
  an exact Vandermonde solve. The alternative, a fraction-free determinant
  over `ℤ[t]` itself, would need a polynomial Bareiss that ncalg does not
  have. Evaluation reuses `exact::det` unchanged.
- **Signature.** Gordon–Litherland with the checkerboard surface.
  - `η` and the crossing type are defined from the A-smoothing corners, so
    no picture convention is needed.
  - Of the two possible type rules, only one is invariant under Markov
    stabilization; the tests pin it.
  - The overall sign is fixed by `σ(right-handed trefoil) = −2`.
- **Homomorphisms to `S_m`.** Wirtinger relations are triples: the images of
  the over-arc and of one under-arc determine the other under-arc. The
  search has three parts:
  - propagation with a worklist and an undo trail;
  - branching on the unknown generator that is the over-arc of the most
    relations with one known under-arc;
  - restricting each component's generators to one conjugacy class, with the
    first generator only at class representatives, weighted by class size.

  On `T(3,7)` with `S_5` this cut the time from 100 s to 1 ms.
- **Khovanov homology.**
  - The cube of resolutions is built with circles numbered by their smallest
    edge, so circles away from the changed crossing keep their identity
    between states.
  - Generators are grouped by `(h, q)`.
  - Ranks over `𝔽_p` use ncalg's `SparseHowell`. Over ℤ they use
    `zmod::sparse_invariant_factors`, added to ncalg for this: Markowitz
    elimination of `±1` pivots, then a dense Smith form of the small rest.
    It took the 10-crossing integral computation from 14 s to 31 ms.

## Additions to ncalg

Both are general and tested in ncalg:
- `linalg::exact::inertia`: Sylvester inertia over an ordered field;
- `zmod::sparse_invariant_factors`.

## Candidates to move to ncalg

`Laurent<C>` is a general one-variable Laurent polynomial ring. If another
library needs one, it belongs in ncalg (as do parts of ncsim's
`phasepoly`).

## Open directions

- HOMFLY-PT through the Hecke algebra and Ocneanu's trace. ncalg's quotient
  algebras need a field. That means either rational-function coefficients
  `ℚ(a, z)` or evaluation and interpolation in two variables, the same trick
  as the Alexander polynomial.
- Bar-Natan's local Khovanov algorithm (delooping and Gaussian elimination
  in the tangle category), for 30+ crossings.
- Reidemeister and Markov simplification before the expensive invariants.
- Seifert matrices from Seifert circles: Tristram–Levine signatures and the
  Alexander module over `ℚ[t^{±1}]`.
