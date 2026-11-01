//! Kirby diagrams of infinite families, with bundles and twist boxes.
//!
//! A [`KirbyDiagram`] is a planar graph whose edges are **cables**: an edge
//! carries an ordered list of [`Group`]s, each a bundle of parallel strands
//! of one handle, all oriented the same way, whose width is a polynomial in
//! the family's parameters ([`Poly`]).  Its nodes are
//!
//! * [`NodeKind::Cross`]: one cable crosses another; every strand of the
//!   under-cable passes under every strand of the over-cable (a grid);
//! * [`NodeKind::Twist`]: half-twists `Δ^e` on a cable, with `e` a
//!   polynomial; the twisted units are single strands or flat ribbons;
//! * [`NodeKind::Split`]: a cable divides into consecutive sub-cables.
//!
//! Every handle is a 1-handle (a dotted circle) or a 2-handle with a framing
//! coefficient, also a polynomial.  Linking numbers, writhes and framings
//! are computed for all parameter values at once ([`KirbyDiagram::linking`]),
//! and moves act on whole bundles ([`KirbyDiagram::cancel`]), so their
//! results hold for every member of the family.  [`KirbyDiagram::expand`]
//! draws one member as an ordinary [`Diagram`], which is how the symbolic
//! results are checked.
//!
//! ## Conventions
//!
//! * A node lists its incident edges **counterclockwise**.  A `Cross` has
//!   four slots: the under-cable joins slots 0 and 2, the over-cable slots 1
//!   and 3.  A `Twist` has slot 0 (in) and slot 1 (out).  A `Split` has the
//!   whole cable at slot 0 and the parts at slots 1, 2, … .
//! * An edge runs from its `tail` to its `head`; its groups are listed left
//!   to right looking along that direction, and a group's `dir` is `+1` if
//!   its strands run tail to head.  At a node, the strands of an edge read
//!   counterclockwise are its strands left to right if the node is the
//!   edge's head, right to left if it is the tail.
//! * In a `Twist`, picture slot 0 at the bottom and slot 1 at the top;
//!   positions are counted left to right at the bottom, and a positive half
//!   twist is the braid `Δ` in [`crate::braid`]'s convention (the left strand
//!   crosses over), positive crossings when the strands are parallel.
//! * In a `Split`, the strands at slot 0 (counterclockwise) are joined to
//!   the strands of the other slots (counterclockwise, slot after slot) in
//!   reverse order, the only planar way.

use crate::diagram::Diagram;
use crate::error::{KnotError, Result};
use crate::param::Poly;
use std::collections::HashMap;

/// A parameter of the family, with its parity if it is fixed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Param {
    /// Its name, like `n`.
    pub name: String,
    /// `Some(true)` if it is always odd, `Some(false)` if always even.
    pub parity: Option<bool>,
}

/// A 1-handle or a 2-handle.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum HandleKind {
    /// A 1-handle, drawn as a dotted circle.
    Dotted,
    /// A 2-handle with this framing coefficient.
    Framed(Poly),
}

/// A handle: one component of the diagram.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Handle {
    /// Its name, like `G`.
    pub name: String,
    /// Dotted or framed.
    pub kind: HandleKind,
}

/// Parallel strands of one handle, oriented alike.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Group {
    /// The handle.
    pub handle: usize,
    /// `+1` if the strands run from the edge's tail to its head, `−1` if back.
    pub dir: i8,
    /// The number of strands.
    pub width: Poly,
}

/// What a node does to the cables through it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NodeKind {
    /// A crossing of two cables.
    Cross,
    /// `Δ^half_twists` on the cable from slot 0 to slot 1.
    Twist {
        /// The number of half-twists (negative: left-handed).
        half_twists: Poly,
        /// What the twist moves around each other.
        units: Units,
    },
    /// A cable dividing into sub-cables.
    Split,
}

/// The units a twist box permutes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Units {
    /// Every strand on its own: `Δ` on all strands.
    Strands,
    /// Flat ribbons: the cable's groups (counterclockwise at slot 0) taken
    /// this many at a time.  `Δ` crosses the ribbons over each other as
    /// grids without twisting them, which is what parallel copies of a
    /// twisted strand look like.
    Ribbons(Vec<usize>),
}

/// A node: its kind and the edges at its slots, counterclockwise.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Node {
    /// The kind.
    pub kind: NodeKind,
    /// The edge at each slot.
    pub slots: Vec<usize>,
}

/// An edge: a cable from one node slot to another.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Edge {
    /// The groups, left to right looking from tail to head.
    pub groups: Vec<Group>,
    /// `(node, slot)` where the edge starts.
    pub tail: (usize, usize),
    /// `(node, slot)` where the edge ends.
    pub head: (usize, usize),
}

/// A Kirby diagram of a family; see the [module docs](self).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct KirbyDiagram {
    /// The parameters.
    pub params: Vec<Param>,
    /// The handles.
    pub handles: Vec<Handle>,
    /// The nodes.
    pub nodes: Vec<Node>,
    /// The edges.
    pub edges: Vec<Edge>,
}

/// A strand group as seen from a node slot.
#[derive(Clone, Debug, PartialEq, Eq)]
struct SlotGroup {
    handle: usize,
    /// `+1` if the strands point into the node.
    into: i8,
    width: Poly,
}

/// Linking numbers and writhes of the handles, for all parameter values.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Linking {
    /// `lk[i][j]` for `i ≠ j`; the diagonal is zero.
    pub lk: Vec<Vec<Poly>>,
    /// The writhe of each handle (its blackboard framing).
    pub writhe: Vec<Poly>,
}

/// One member of a family drawn as an ordinary diagram.
#[derive(Clone, Debug)]
pub struct Expanded {
    /// The diagram.
    pub diagram: Diagram,
    /// The handle of each component of `diagram`.
    pub handle_of_component: Vec<usize>,
}

impl KirbyDiagram {
    /// An empty diagram with these parameters: `(name, parity)`, parity
    /// `Some(true)` for always odd.
    pub fn new(params: &[(&str, Option<bool>)]) -> Self {
        KirbyDiagram {
            params: params
                .iter()
                .map(|&(name, parity)| Param {
                    name: name.into(),
                    parity,
                })
                .collect(),
            handles: Vec::new(),
            nodes: Vec::new(),
            edges: Vec::new(),
        }
    }

    /// The parameter called `name`, as a polynomial.
    pub fn param(&self, name: &str) -> Poly {
        let i = self
            .params
            .iter()
            .position(|p| p.name == name)
            .unwrap_or_else(|| panic!("no parameter {name}"));
        Poly::var(i)
    }

    /// The parameter names.
    pub fn param_names(&self) -> Vec<String> {
        self.params.iter().map(|p| p.name.clone()).collect()
    }

    /// Adds a dotted circle; returns its index.
    pub fn dotted(&mut self, name: &str) -> usize {
        self.handles.push(Handle {
            name: name.into(),
            kind: HandleKind::Dotted,
        });
        self.handles.len() - 1
    }

    /// Adds a 2-handle with this framing; returns its index.
    pub fn framed(&mut self, name: &str, framing: Poly) -> usize {
        self.handles.push(Handle {
            name: name.into(),
            kind: HandleKind::Framed(framing),
        });
        self.handles.len() - 1
    }

    /// The handle called `name`.
    pub fn handle(&self, name: &str) -> Option<usize> {
        self.handles.iter().position(|h| h.name == name)
    }

    /// Adds a node with `slots` slots, to be filled by [`edge`](Self::edge).
    pub fn node(&mut self, kind: NodeKind, slots: usize) -> usize {
        self.nodes.push(Node {
            kind,
            slots: vec![usize::MAX; slots],
        });
        self.nodes.len() - 1
    }

    /// Adds an edge from `tail` to `head` (`(node, slot)` each) carrying
    /// `groups` (`(handle, dir, width)`, left to right); returns its index.
    pub fn edge(
        &mut self,
        groups: &[(usize, i8, Poly)],
        tail: (usize, usize),
        head: (usize, usize),
    ) -> usize {
        let e = self.edges.len();
        for &(n, s) in &[tail, head] {
            assert_eq!(
                self.nodes[n].slots[s],
                usize::MAX,
                "slot {s} of node {n} is taken"
            );
            self.nodes[n].slots[s] = e;
        }
        self.edges.push(Edge {
            groups: groups
                .iter()
                .map(|(h, d, w)| Group {
                    handle: *h,
                    dir: *d,
                    width: w.clone(),
                })
                .collect(),
            tail,
            head,
        });
        e
    }

    /// Checks that every slot is filled and the slot counts fit the kinds.
    pub fn validate(&self) -> Result<()> {
        for (i, n) in self.nodes.iter().enumerate() {
            let ok = match n.kind {
                NodeKind::Cross => n.slots.len() == 4,
                NodeKind::Twist { .. } => n.slots.len() == 2,
                NodeKind::Split => !n.slots.is_empty(),
            };
            if !ok {
                return Err(KnotError::Invalid(format!(
                    "node {i} has {} slots",
                    n.slots.len()
                )));
            }
            if let Some(s) = n.slots.iter().position(|&e| e == usize::MAX) {
                return Err(KnotError::Invalid(format!("slot {s} of node {i} is empty")));
            }
        }
        for (e, edge) in self.edges.iter().enumerate() {
            for &(n, s) in &[edge.tail, edge.head] {
                if self.nodes[n].slots[s] != e {
                    return Err(KnotError::Invalid(format!("edge {e} and node {n} disagree")));
                }
            }
            for g in &edge.groups {
                if g.handle >= self.handles.len() || g.dir.abs() != 1 {
                    return Err(KnotError::Invalid(format!("bad group on edge {e}")));
                }
            }
        }
        Ok(())
    }

    /// The groups of the edge at `(node, slot)`, counterclockwise.
    fn slot_groups(&self, node: usize, slot: usize) -> Vec<SlotGroup> {
        let e = &self.edges[self.nodes[node].slots[slot]];
        let at_head = e.head == (node, slot);
        let gs = e.groups.iter().map(|g| SlotGroup {
            handle: g.handle,
            into: if at_head { g.dir } else { -g.dir },
            width: g.width.clone(),
        });
        if at_head {
            gs.collect()
        } else {
            gs.rev().collect()
        }
    }

    /// Linking numbers and writhes, as polynomials in the parameters.
    pub fn linking(&self) -> Linking {
        let h = self.handles.len();
        let mut lk = vec![vec![Poly::zero(); h]; h];
        let mut writhe = vec![Poly::zero(); h];
        let half = Poly::rational(num_rational::BigRational::new(1.into(), 2.into()));
        let mut add = |a: usize, b: usize, c: Poly| {
            if a == b {
                writhe[a] = &writhe[a] + &c;
            } else {
                let c = &c * &half;
                lk[a][b] = &lk[a][b] + &c;
                lk[b][a] = &lk[b][a] + &c;
            }
        };
        for (i, n) in self.nodes.iter().enumerate() {
            match &n.kind {
                NodeKind::Cross => {
                    // under from slot 0 to 2, over from 1 to 3: negative for
                    // strands along those directions
                    for g in self.slot_groups(i, 0) {
                        for o in self.slot_groups(i, 1) {
                            let s = -(g.into as i64) * o.into as i64;
                            add(g.handle, o.handle, &(&g.width * &o.width) * s);
                        }
                    }
                }
                NodeKind::Twist { half_twists, units } => {
                    let gs = self.slot_groups(i, 0);
                    let unit = unit_of_groups(units, gs.len());
                    for a in 0..gs.len() {
                        if matches!(units, Units::Strands) {
                            let w = &gs[a].width;
                            let pairs = &(w * &(w - 1)) * &half;
                            add(gs[a].handle, gs[a].handle, half_twists * &pairs);
                        }
                        for b in a + 1..gs.len() {
                            if unit[a] == unit[b] && !matches!(units, Units::Strands) {
                                continue;
                            }
                            let s = gs[a].into as i64 * gs[b].into as i64;
                            let c = &(&gs[a].width * &gs[b].width) * half_twists;
                            add(gs[a].handle, gs[b].handle, c * s);
                        }
                    }
                }
                NodeKind::Split => {}
            }
        }
        Linking { lk, writhe }
    }

    /// The framing coefficient of a 2-handle (`None` for a dotted circle).
    pub fn framing(&self, h: usize) -> Option<&Poly> {
        match &self.handles[h].kind {
            HandleKind::Framed(f) => Some(f),
            HandleKind::Dotted => None,
        }
    }

    /// The linking matrix with the framings on the diagonal and `0` for the
    /// dotted circles: the surgery description of the boundary.
    pub fn linking_matrix(&self) -> Vec<Vec<Poly>> {
        let mut m = self.linking().lk;
        for (i, row) in m.iter_mut().enumerate() {
            row[i] = self.framing(i).cloned().unwrap_or_default();
        }
        m
    }

    /// Checks the parameter values against the declared parities.
    fn check_values(&self, values: &[i64]) -> Result<()> {
        if values.len() != self.params.len() {
            return Err(KnotError::Invalid(format!(
                "{} values for {} parameters",
                values.len(),
                self.params.len()
            )));
        }
        for (p, &v) in self.params.iter().zip(values) {
            if let Some(odd) = p.parity {
                if (v.rem_euclid(2) == 1) != odd {
                    return Err(KnotError::Invalid(format!(
                        "{} must be {}",
                        p.name,
                        if odd { "odd" } else { "even" }
                    )));
                }
            }
        }
        Ok(())
    }

    /// The member of the family at these parameter values, as an ordinary
    /// diagram (every bundle drawn as parallel strands).
    pub fn expand(&self, values: &[i64]) -> Result<Expanded> {
        self.validate()?;
        self.check_values(values)?;
        let mut ex = Expansion::default();
        // label of strand k (left to right) of each edge, its handle and dir
        let mut base = Vec::with_capacity(self.edges.len());
        let mut strands: Vec<Vec<(usize, i8)>> = Vec::with_capacity(self.edges.len());
        for (e, edge) in self.edges.iter().enumerate() {
            let mut s = Vec::new();
            for g in &edge.groups {
                let w = g.width.eval_int(values)?;
                if w < 0 {
                    return Err(KnotError::Invalid(format!("negative width {w} on edge {e}")));
                }
                s.extend(std::iter::repeat_n((g.handle, g.dir), w as usize));
            }
            base.push(ex.labels.len());
            for &(h, _) in &s {
                ex.fresh(h);
            }
            strands.push(s);
        }
        // the strands at a slot, counterclockwise: (label, handle, points into the node)
        let at = |node: usize, slot: usize| -> Vec<(usize, usize, bool)> {
            let e = self.nodes[node].slots[slot];
            let at_head = self.edges[e].head == (node, slot);
            let s = &strands[e];
            let one = |k: usize| {
                let (h, d) = s[k];
                (base[e] + k, h, (d > 0) == at_head)
            };
            if at_head {
                (0..s.len()).map(one).collect()
            } else {
                (0..s.len()).rev().map(one).collect()
            }
        };
        for (i, n) in self.nodes.iter().enumerate() {
            match &n.kind {
                NodeKind::Cross => {
                    let (s0, s1, s2, s3) = (at(i, 0), at(i, 1), at(i, 2), at(i, 3));
                    if s0.len() != s2.len() || s1.len() != s3.len() {
                        return Err(KnotError::Invalid(format!(
                            "crossing {i}: the cables change width"
                        )));
                    }
                    let (w0, w1) = (s0.len(), s1.len());
                    // under strand i runs from s0[i] up to s2[w0-1-i], over strand j
                    // from s1[j] (east) to s3[w1-1-j] (west)
                    let mut u = Vec::with_capacity(w0);
                    for k in 0..w0 {
                        u.push(ex.path(s0[k], s2[w0 - 1 - k], w1, i)?);
                    }
                    let mut o = Vec::with_capacity(w1);
                    for j in 0..w1 {
                        o.push(ex.path(s1[j], s3[w1 - 1 - j], w0, i)?);
                    }
                    for k in 0..w0 {
                        for j in 0..w1 {
                            let ccw = [u[k][j], o[j][w0 - 1 - k], u[k][j + 1], o[j][w0 - k]];
                            ex.crossing(ccw, s0[k].2);
                        }
                    }
                }
                NodeKind::Twist { half_twists, units } => {
                    let (s0, s1) = (at(i, 0), at(i, 1));
                    if s0.len() != s1.len() {
                        return Err(KnotError::Invalid(format!("twist {i}: the cable changes width")));
                    }
                    let e = half_twists.eval_int(values)?;
                    let widths: Vec<usize> = match units {
                        Units::Strands => vec![1; s0.len()],
                        Units::Ribbons(counts) => {
                            let e0 = self.nodes[i].slots[0];
                            let gw: Vec<usize> = {
                                let edge = &self.edges[e0];
                                let mut v: Vec<usize> = edge
                                    .groups
                                    .iter()
                                    .map(|g| g.width.eval_int(values).map(|w| w as usize))
                                    .collect::<Result<_>>()?;
                                if edge.head != (i, 0) {
                                    v.reverse();
                                }
                                v
                            };
                            if counts.iter().sum::<usize>() != gw.len() {
                                return Err(KnotError::Invalid(format!(
                                    "twist {i}: units do not cover the cable"
                                )));
                            }
                            let mut out = Vec::new();
                            let mut k = 0;
                            for &c in counts {
                                out.push(gw[k..k + c].iter().sum());
                                k += c;
                            }
                            out
                        }
                    };
                    let word = delta_power(&widths, e);
                    let mut cur: Vec<(usize, usize, bool)> = s0.clone();
                    for g in word {
                        let l = g.unsigned_abs() as usize - 1;
                        let (a, b) = (cur[l], cur[l + 1]);
                        // the strand from the right ends at l, the one from the left at l+1
                        let (nl, nr) = (ex.fresh(b.1), ex.fresh(a.1));
                        if g > 0 {
                            // left strand over, to l+1; right strand under, to l
                            ex.crossing([b.0, nr, nl, a.0], b.2);
                            cur[l] = (nl, b.1, b.2);
                            cur[l + 1] = (nr, a.1, a.2);
                        } else {
                            // right strand over, to l; left strand under, to l+1
                            ex.crossing([a.0, b.0, nr, nl], a.2);
                            cur[l] = (nl, b.1, b.2);
                            cur[l + 1] = (nr, a.1, a.2);
                        }
                    }
                    let w = s0.len();
                    for p in 0..w {
                        // cur[p] runs up into slot 1's strand, which must leave
                        ex.join(cur[p], s1[w - 1 - p], i)?;
                    }
                }
                NodeKind::Split => {
                    let s0 = at(i, 0);
                    let rest: Vec<_> = (1..n.slots.len()).flat_map(|s| at(i, s)).collect();
                    if s0.len() != rest.len() {
                        return Err(KnotError::Invalid(format!(
                            "split {i}: {} strands become {}",
                            s0.len(),
                            rest.len()
                        )));
                    }
                    let w = s0.len();
                    for k in 0..w {
                        ex.join(s0[k], rest[w - 1 - k], i)?;
                    }
                }
            }
        }
        ex.finish(self.handles.len())
    }

    /// Checks the symbolic linking numbers and writhes against the member at
    /// `values`, and that every handle is one component there.
    pub fn check(&self, values: &[i64]) -> Result<()> {
        let ex = self.expand(values)?;
        let d = &ex.diagram;
        let h = self.handles.len();
        let mut comp_of = vec![None; h];
        for (c, &hh) in ex.handle_of_component.iter().enumerate() {
            if comp_of[hh].replace(c).is_some() {
                return Err(KnotError::Invalid(format!(
                    "handle {} is more than one component at {values:?}",
                    self.handles[hh].name
                )));
            }
        }
        let comp_of: Vec<usize> = comp_of
            .into_iter()
            .enumerate()
            .map(|(i, c)| {
                c.ok_or_else(|| KnotError::Invalid(format!("handle {} is empty", self.handles[i].name)))
            })
            .collect::<Result<_>>()?;
        let lkm = d.linking_matrix();
        let mut writhe = vec![0i64; d.components()];
        for c in 0..d.crossings() {
            let x = d.crossing(c);
            let (a, b) = (d.edge_component(x[0]), d.edge_component(x[1]));
            if a == b {
                writhe[a] += d.sign(c);
            }
        }
        let sym = self.linking();
        for i in 0..h {
            let w = sym.writhe[i].eval_int(values)?;
            if w != writhe[comp_of[i]] {
                return Err(KnotError::Invalid(format!(
                    "writhe of {}: {} predicted, {} drawn at {values:?}",
                    self.handles[i].name, w, writhe[comp_of[i]]
                )));
            }
            for j in 0..h {
                if i != j {
                    let l = sym.lk[i][j].eval_int(values)?;
                    if l != lkm[comp_of[i]][comp_of[j]] {
                        return Err(KnotError::Invalid(format!(
                            "lk({}, {}): {} predicted, {} drawn at {values:?}",
                            self.handles[i].name, self.handles[j].name, l, lkm[comp_of[i]][comp_of[j]]
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    /// The framings at `values`, one per component of [`expand`](Self::expand)'s
    /// diagram, with `0` for dotted circles: the surgery description of the
    /// boundary 3-manifold.
    pub fn surgery_framings(&self, ex: &Expanded, values: &[i64]) -> Result<Vec<i64>> {
        ex.handle_of_component
            .iter()
            .map(|&h| match &self.handles[h].kind {
                HandleKind::Dotted => Ok(0),
                HandleKind::Framed(f) => f.eval_int(values),
            })
            .collect()
    }

    /// Cancels the dotted circle `dotted` against the 2-handle `handle`,
    /// which passes through it once.  First every other strand through the
    /// dotted circle slides over `handle` (each becomes a parallel copy of
    /// `handle`'s arc, with `handle`'s framing), then both are erased.
    ///
    /// The dotted circle must bound a clean disk: inside it, strands only go
    /// straight through (over one side of the circle, under the other) and
    /// side by side, as dotted circles are usually drawn.  `handle` must have
    /// width 1 and not run through twist boxes with strand units of
    /// unknown width.  The new framings follow from homology,
    /// `f' = f − 2k·lk(X, handle) + k²·f_handle` with `k = lk(X, dotted) /
    /// lk(handle, dotted)`, and the new linking numbers are checked against
    /// the same formula.
    pub fn cancel(&self, dotted: usize, handle: usize) -> Result<KirbyDiagram> {
        self.validate()?;
        if self.handles[dotted].kind != HandleKind::Dotted {
            return Err(KnotError::Invalid(format!(
                "{} is not a dotted circle",
                self.handles[dotted].name
            )));
        }
        let f_b = self
            .framing(handle)
            .cloned()
            .ok_or_else(|| KnotError::Invalid(format!("{} is not a 2-handle", self.handles[handle].name)))?;
        let one = Poly::constant(1);
        let before = self.linking();
        let lk_br = before.lk[handle][dotted]
            .as_integer()
            .filter(|v| v.abs() == 1)
            .ok_or_else(|| {
                KnotError::Invalid(format!(
                    "{} links {} {} times, not once",
                    self.handles[handle].name,
                    self.handles[dotted].name,
                    before.lk[handle][dotted].display(&self.param_names())
                ))
            })?;

        // walk around the dotted circle: (node, slot it enters by)
        let r_edges: Vec<usize> = (0..self.edges.len())
            .filter(|&e| self.edges[e].groups.iter().any(|g| g.handle == dotted))
            .collect();
        for &e in &r_edges {
            let gs = &self.edges[e].groups;
            if gs.len() != 1 || gs[0].width != one {
                return Err(KnotError::Unsupported(
                    "the dotted circle must be a single strand".into(),
                ));
            }
        }
        let mut walk: Vec<(usize, usize)> = Vec::new();
        let e0 = r_edges[0];
        let mut e = e0;
        loop {
            let edge = &self.edges[e];
            let (node, slot) = if edge.groups[0].dir > 0 {
                edge.head
            } else {
                edge.tail
            };
            if self.nodes[node].kind != NodeKind::Cross {
                return Err(KnotError::Unsupported(
                    "the dotted circle may only cross other strands".into(),
                ));
            }
            walk.push((node, slot));
            e = self.nodes[node].slots[(slot + 2) % 4];
            if e == e0 {
                break;
            }
            if walk.len() > r_edges.len() {
                return Err(KnotError::Invalid("the dotted circle does not close up".into()));
            }
        }
        let k2 = walk.len();
        let pos: HashMap<usize, usize> = walk.iter().enumerate().map(|(i, &(n, _))| (n, i)).collect();
        if pos.len() != k2 {
            return Err(KnotError::Unsupported("the dotted circle crosses itself".into()));
        }
        // the interior side: +1 right of the walk, −1 left; the strand there
        // must run straight to another crossing with the circle
        let inner_slot = |i: usize, side: i32| (walk[i].1 + if side > 0 { 1 } else { 3 }) % 4;
        let mut partner = None;
        for side in [1, -1] {
            let mut p = vec![usize::MAX; k2];
            let mut ok = true;
            for i in 0..k2 {
                let (node, _) = walk[i];
                let s = inner_slot(i, side);
                let edge = &self.edges[self.nodes[node].slots[s]];
                let other = if edge.tail == (node, s) {
                    edge.head
                } else {
                    edge.tail
                };
                match pos.get(&other.0) {
                    Some(&j) if j != i && inner_slot(j, side) == other.1 => p[i] = j,
                    _ => {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                partner = Some((side, p));
                break;
            }
        }
        let (side, partner) = partner.ok_or_else(|| {
            KnotError::Unsupported("the strands inside the dotted circle must go straight through it".into())
        })?;
        // over-ends: crossings where the strand passes over the circle
        let over_end = |i: usize| walk[i].1 % 2 == 0; // the circle is under
        for i in 0..k2 {
            if over_end(i) == over_end(partner[i]) {
                return Err(KnotError::Unsupported(
                    "a strand crosses the dotted circle's disk without passing through it".into(),
                ));
            }
        }
        // the over-ends must be consecutive, the arcs parallel
        let start = (0..k2)
            .find(|&i| over_end(i) && !over_end((i + k2 - 1) % k2))
            .ok_or_else(|| KnotError::Invalid("no strand passes through the dotted circle".into()))?;
        let k = k2 / 2;
        let order: Vec<usize> = (0..k2).map(|t| (start + t) % k2).collect();
        for t in 0..k {
            if !over_end(order[t]) || partner[order[t]] != order[k2 - 1 - t] {
                return Err(KnotError::Unsupported(
                    "the strands through the dotted circle must run side by side".into(),
                ));
            }
        }
        // arcs left to right, looking out of the disk on the over side
        let mut arcs: Vec<(usize, usize)> = (0..k).map(|t| (order[t], order[k2 - 1 - t])).collect();
        if side < 0 {
            arcs.reverse();
        }
        let arc_edge = |i: usize| self.nodes[walk[i].0].slots[inner_slot(i, side)];
        let outer = |i: usize| (walk[i].0, (inner_slot(i, side) + 2) % 4);
        let p = arcs
            .iter()
            .position(|&(o, _)| self.edges[arc_edge(o)].groups.iter().any(|g| g.handle == handle))
            .ok_or_else(|| {
                KnotError::Invalid(format!(
                    "{} does not pass through the dotted circle",
                    self.handles[handle].name
                ))
            })?;
        let b_arc = &self.edges[arc_edge(arcs[p].0)];
        if b_arc.groups.len() != 1 || b_arc.groups[0].width != one {
            return Err(KnotError::Unsupported(format!(
                "{} must pass through the dotted circle alone, as one strand",
                self.handles[handle].name
            )));
        }
        for (j, &(o, _)) in arcs.iter().enumerate() {
            if j != p && self.edges[arc_edge(o)].groups.iter().any(|g| g.handle == handle) {
                return Err(KnotError::Unsupported(format!(
                    "{} passes through the dotted circle more than once",
                    self.handles[handle].name
                )));
            }
        }
        // slots of the two new junctions, after the copy cable at slot 0
        let over_order: Vec<usize> = (0..p).rev().chain((p + 1..k).rev()).collect();
        let under_order: Vec<usize> = (p + 1..k).chain(0..p).collect();
        let groups_at = |end: (usize, usize)| -> Vec<SlotGroup> { self.slot_groups(end.0, end.1) };
        let copies = |ends: &[(usize, usize)]| -> Vec<SlotGroup> {
            let rest: Vec<SlotGroup> = ends.iter().flat_map(|&e| groups_at(e)).collect();
            rest.into_iter()
                .rev()
                .map(|g| SlotGroup { into: -g.into, ..g })
                .collect()
        };
        let over_ends: Vec<(usize, usize)> = over_order.iter().map(|&j| outer(arcs[j].0)).collect();
        let under_ends: Vec<(usize, usize)> = under_order.iter().map(|&j| outer(arcs[j].1)).collect();
        let c_ccw = copies(&over_ends);
        // the copy cable on handle's outer edge on the over side, then along
        // handle's direction of travel: (handle, dir relative to it, width)
        let b_over = outer(arcs[p].0);
        let b_edge = self.nodes[b_over.0].slots[b_over.1];
        let at_head = self.edges[b_edge].head == b_over;
        let mut c_lr: Vec<Group> = c_ccw
            .iter()
            .map(|g| Group {
                handle: g.handle,
                dir: if at_head { g.into } else { -g.into },
                width: g.width.clone(),
            })
            .collect();
        if !at_head {
            c_lr.reverse();
        }
        let bg = &self.edges[b_edge].groups;
        if bg.len() != 1 || bg[0].handle != handle {
            return Err(KnotError::Unsupported(format!(
                "{} must leave the dotted circle alone",
                self.handles[handle].name
            )));
        }
        let c_travel = orient(&c_lr, bg[0].dir);

        // the new diagram: drop the circle, its crossings and the arcs inside
        let w_b = before.writhe[handle].clone();
        let r_nodes: Vec<usize> = walk.iter().map(|w| w.0).collect();
        let dead_edges: Vec<usize> = r_edges
            .iter()
            .copied()
            .chain(arcs.iter().map(|&(o, _)| arc_edge(o)))
            .collect();
        let mut out = KirbyDiagram {
            params: self.params.clone(),
            handles: self.handles.clone(),
            nodes: Vec::new(),
            edges: Vec::new(),
        };
        let mut node_map = vec![usize::MAX; self.nodes.len()];
        for (i, n) in self.nodes.iter().enumerate() {
            if r_nodes.contains(&i) {
                continue;
            }
            node_map[i] = out.nodes.len();
            let kind = match &n.kind {
                NodeKind::Twist { half_twists, units } => {
                    let gs = self.slot_groups(i, 0);
                    let has_b = gs.iter().any(|g| g.handle == handle);
                    let units = match units {
                        _ if !has_b => units.clone(),
                        Units::Strands => {
                            if gs.iter().any(|g| g.width != one) {
                                return Err(KnotError::Unsupported(format!(
                                    "{} runs through a twist box of strands with other bundles",
                                    self.handles[handle].name
                                )));
                            }
                            Units::Ribbons(
                                gs.iter()
                                    .map(|g| if g.handle == handle { c_travel.len() } else { 1 })
                                    .collect(),
                            )
                        }
                        Units::Ribbons(counts) => {
                            let mut v = Vec::new();
                            let mut j = 0;
                            for &c in counts {
                                v.push(
                                    gs[j..j + c]
                                        .iter()
                                        .map(|g| if g.handle == handle { c_travel.len() } else { 1 })
                                        .sum(),
                                );
                                j += c;
                            }
                            Units::Ribbons(v)
                        }
                    };
                    NodeKind::Twist {
                        half_twists: half_twists.clone(),
                        units,
                    }
                }
                kind => kind.clone(),
            };
            out.nodes.push(Node {
                kind,
                slots: vec![usize::MAX; n.slots.len()],
            });
        }
        let junction_over = out.node(NodeKind::Split, k);
        let junction_under = out.node(NodeKind::Split, k);
        let mut new_end = HashMap::new();
        new_end.insert(b_over, (junction_over, 0));
        new_end.insert(outer(arcs[p].1), (junction_under, 0));
        for (s, &end) in over_ends.iter().enumerate() {
            new_end.insert(end, (junction_over, s + 1));
        }
        for (s, &end) in under_ends.iter().enumerate() {
            new_end.insert(end, (junction_under, s + 1));
        }
        let map_end = |end: (usize, usize)| new_end.get(&end).copied().unwrap_or((node_map[end.0], end.1));
        for (e, edge) in self.edges.iter().enumerate() {
            if dead_edges.contains(&e) {
                continue;
            }
            let mut groups = Vec::new();
            for g in &edge.groups {
                if g.handle == handle {
                    if g.width != one {
                        return Err(KnotError::Unsupported(format!(
                            "{} must have width 1",
                            self.handles[handle].name
                        )));
                    }
                    groups.extend(orient(&c_travel, g.dir));
                } else {
                    groups.push(g.clone());
                }
            }
            let g: Vec<(usize, i8, Poly)> = groups.into_iter().map(|g| (g.handle, g.dir, g.width)).collect();
            out.edge(&g, map_end(edge.tail), map_end(edge.head));
        }
        // the copies must fit the under side too
        let under_rest: Vec<SlotGroup> = (1..=k - 1)
            .flat_map(|s| out.slot_groups(junction_under, s))
            .collect();
        let want: Vec<SlotGroup> = under_rest
            .into_iter()
            .rev()
            .map(|g| SlotGroup { into: -g.into, ..g })
            .collect();
        if merge(&want) != merge(&out.slot_groups(junction_under, 0)) {
            return Err(KnotError::Invalid(
                "internal error: the copies do not fit on the under side".into(),
            ));
        }
        // handle's framing: f_b − w_b full twists on the copies
        let twists = &(&f_b - &w_b) * 2;
        if !twists.is_zero() {
            let c_edge = out.nodes[junction_over].slots[0];
            let edge = out.edges[c_edge].clone();
            let t = out.node(
                NodeKind::Twist {
                    half_twists: twists,
                    units: Units::Strands,
                },
                2,
            );
            // c_edge now ends at the twist; a new edge carries on
            let old_head = edge.head;
            out.edges[c_edge].head = (t, 0);
            out.nodes[t].slots[0] = c_edge;
            out.nodes[old_head.0].slots[old_head.1] = usize::MAX;
            let g: Vec<(usize, i8, Poly)> = edge
                .groups
                .iter()
                .map(|g| (g.handle, g.dir, g.width.clone()))
                .collect();
            out.edge(&g, (t, 1), old_head);
        }
        // new framings, then drop the two handles
        let kx: Vec<Poly> = (0..self.handles.len())
            .map(|x| &before.lk[x][dotted] * lk_br)
            .collect();
        for x in 0..self.handles.len() {
            if let HandleKind::Framed(f) = &self.handles[x].kind {
                if x != handle {
                    let nf = &(f - &(&(&kx[x] * &before.lk[x][handle]) * 2)) + &(&(&kx[x] * &kx[x]) * &f_b);
                    out.handles[x].kind = HandleKind::Framed(nf);
                }
            }
        }
        let keep: Vec<usize> = (0..self.handles.len())
            .filter(|&x| x != dotted && x != handle)
            .collect();
        let mut hmap = vec![usize::MAX; self.handles.len()];
        for (i, &x) in keep.iter().enumerate() {
            hmap[x] = i;
        }
        out.handles = keep.iter().map(|&x| out.handles[x].clone()).collect();
        for edge in &mut out.edges {
            for g in &mut edge.groups {
                g.handle = hmap[g.handle];
            }
        }
        out.validate()?;
        // the linking numbers must follow lk(X − k_X B, Y − k_Y B)
        let after = out.linking();
        for (i, &x) in keep.iter().enumerate() {
            for (j, &y) in keep.iter().enumerate() {
                if i == j {
                    continue;
                }
                let want = &(&(&before.lk[x][y] - &(&kx[y] * &before.lk[x][handle]))
                    - &(&kx[x] * &before.lk[y][handle]))
                    + &(&(&kx[x] * &kx[y]) * &f_b);
                if after.lk[i][j] != want {
                    return Err(KnotError::Invalid(format!(
                        "internal error: lk({}, {}) is {}, homology says {}",
                        out.handles[i].name,
                        out.handles[j].name,
                        after.lk[i][j].display(&self.param_names()),
                        want.display(&self.param_names())
                    )));
                }
            }
        }
        Ok(out)
    }

    /// A summary: the handles with their framings, and the linking numbers
    /// and writhes, as polynomials in the parameters.
    pub fn describe(&self) -> String {
        let names = self.param_names();
        let l = self.linking();
        let mut s = String::new();
        for (i, h) in self.handles.iter().enumerate() {
            let kind = match &h.kind {
                HandleKind::Dotted => "dotted circle".to_string(),
                HandleKind::Framed(f) => format!("framing {}", f.display(&names)),
            };
            s += &format!("{}: {kind}, writhe {}\n", h.name, l.writhe[i].display(&names));
        }
        for i in 0..self.handles.len() {
            for j in i + 1..self.handles.len() {
                s += &format!(
                    "lk({}, {}) = {}\n",
                    self.handles[i].name,
                    self.handles[j].name,
                    l.lk[i][j].display(&names)
                );
            }
        }
        s
    }
}

/// A cable list for a handle running with (`dir = 1`) or against (`−1`) the
/// direction it was given in.
fn orient(gs: &[Group], dir: i8) -> Vec<Group> {
    if dir > 0 {
        gs.to_vec()
    } else {
        gs.iter()
            .rev()
            .map(|g| Group {
                dir: -g.dir,
                ..g.clone()
            })
            .collect()
    }
}

/// Adjacent groups of one handle and direction combined.
fn merge(gs: &[SlotGroup]) -> Vec<SlotGroup> {
    let mut out: Vec<SlotGroup> = Vec::new();
    for g in gs {
        if g.width.is_zero() {
            continue;
        }
        match out.last_mut() {
            Some(l) if l.handle == g.handle && l.into == g.into => l.width = &l.width + &g.width,
            _ => out.push(g.clone()),
        }
    }
    out
}

/// The unit of each group of a twist's cable (counterclockwise at slot 0).
fn unit_of_groups(units: &Units, n: usize) -> Vec<usize> {
    match units {
        Units::Strands => (0..n).collect(),
        Units::Ribbons(counts) => counts
            .iter()
            .enumerate()
            .flat_map(|(u, &c)| std::iter::repeat_n(u, c))
            .collect(),
    }
}

/// `Δ^e` on units of these widths, as a braid word on all the strands
/// (`i` is `σ_i`, 1-based; the left strand crosses over).
fn delta_power(widths: &[usize], e: i64) -> Vec<i32> {
    let mut out = Vec::new();
    let mut cur = widths.to_vec();
    for _ in 0..e.unsigned_abs() {
        if e > 0 {
            out.extend(delta_once(&cur));
        } else {
            // Δ^{-1} here is the inverse of Δ on the reversed widths
            let rev: Vec<usize> = cur.iter().rev().copied().collect();
            out.extend(delta_once(&rev).iter().rev().map(|g| -g));
        }
        cur.reverse();
    }
    out
}

/// One `Δ` on units of these widths:
/// `(U_1 ⋯ U_{k−1})(U_1 ⋯ U_{k−2}) ⋯ U_1`, `U_i` crossing unit `i` over unit `i+1`.
fn delta_once(widths: &[usize]) -> Vec<i32> {
    let mut word = Vec::new();
    let mut w = widths.to_vec();
    let k = w.len();
    for top in (1..k).rev() {
        for i in 0..top {
            let p: usize = w[..i].iter().sum();
            let (a, b) = (w[i], w[i + 1]);
            for s in (0..a).rev() {
                for t in 0..b {
                    word.push((p + s + t + 1) as i32);
                }
            }
            w.swap(i, i + 1);
        }
    }
    word
}

/// Labels and crossings of an expansion under construction.
#[derive(Default)]
struct Expansion {
    /// Union-find parent and handle of each label.
    labels: Vec<(usize, usize)>,
    /// Crossings, counterclockwise from the incoming under-strand.
    code: Vec<[usize; 4]>,
}

impl Expansion {
    fn fresh(&mut self, handle: usize) -> usize {
        self.labels.push((self.labels.len(), handle));
        self.labels.len() - 1
    }

    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.labels[r].0 != r {
            r = self.labels[r].0;
        }
        let mut y = x;
        while self.labels[y].0 != r {
            let n = self.labels[y].0;
            self.labels[y].0 = r;
            y = n;
        }
        r
    }

    /// Joins the strand ending at `a` to the one starting at `b` (one must
    /// point into the node, the other out).
    fn join(&mut self, a: (usize, usize, bool), b: (usize, usize, bool), node: usize) -> Result<()> {
        if a.1 != b.1 {
            return Err(KnotError::Invalid(format!(
                "node {node} joins strands of two handles"
            )));
        }
        if a.2 == b.2 {
            return Err(KnotError::Invalid(format!(
                "node {node}: the orientations do not match"
            )));
        }
        let (ra, rb) = (self.find(a.0), self.find(b.0));
        self.labels[ra].0 = rb;
        Ok(())
    }

    /// The labels along a strand from boundary `a` to boundary `b` through
    /// `cuts` crossings.
    fn path(
        &mut self,
        a: (usize, usize, bool),
        b: (usize, usize, bool),
        cuts: usize,
        node: usize,
    ) -> Result<Vec<usize>> {
        if a.1 != b.1 {
            return Err(KnotError::Invalid(format!(
                "node {node} joins strands of two handles"
            )));
        }
        if a.2 == b.2 {
            return Err(KnotError::Invalid(format!(
                "node {node}: the orientations do not match"
            )));
        }
        if cuts == 0 {
            self.join(a, b, node)?;
            return Ok(vec![a.0]);
        }
        let mut v = vec![a.0];
        for _ in 1..cuts {
            v.push(self.fresh(a.1));
        }
        v.push(b.0);
        Ok(v)
    }

    /// A crossing from its labels counterclockwise starting with an
    /// under-strand end; `forward`: the under-strand leaves through `ccw[2]`.
    fn crossing(&mut self, ccw: [usize; 4], forward: bool) {
        self.code.push(if forward {
            ccw
        } else {
            [ccw[2], ccw[3], ccw[0], ccw[1]]
        });
    }

    fn finish(mut self, n_handles: usize) -> Result<Expanded> {
        let code: Vec<[usize; 4]> = self
            .code
            .clone()
            .iter()
            .map(|x| x.map(|l| self.find(l)))
            .collect();
        let mut used = vec![false; self.labels.len()];
        for x in &code {
            for &l in x {
                used[l] = true;
            }
        }
        // closed strands that never cross: free loops, by handle
        let mut free: Vec<usize> = Vec::new();
        let mut seen = vec![false; self.labels.len()];
        for l in 0..self.labels.len() {
            let r = self.find(l);
            if !used[r] && !seen[r] {
                seen[r] = true;
                free.push(self.labels[r].1);
            }
        }
        let pd: Vec<[i64; 4]> = code.iter().map(|x| x.map(|l| l as i64)).collect();
        let diagram = Diagram::build(&pd, free.len())?;
        let mut handle_of_component = vec![usize::MAX; diagram.components()];
        for (c, x) in code.iter().enumerate() {
            let comp = diagram.edge_component(diagram.crossing(c)[0]);
            let h = self.labels[x[0]].1;
            if handle_of_component[comp] != usize::MAX && handle_of_component[comp] != h {
                return Err(KnotError::Invalid("a component carries two handles".into()));
            }
            handle_of_component[comp] = h;
        }
        let n_cross = diagram.components() - diagram.free_loops();
        for (j, h) in free.into_iter().enumerate() {
            handle_of_component[n_cross + j] = h;
        }
        debug_assert!(handle_of_component.iter().all(|&h| h < n_handles));
        Ok(Expanded {
            diagram,
            handle_of_component,
        })
    }
}
