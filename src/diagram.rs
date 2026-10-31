//! Oriented link diagrams as planar diagram (PD) codes.
//!
//! A crossing is `X[a, b, c, d]`: the four edge labels met going
//! **counterclockwise** around it, starting from the **incoming
//! under-strand** (KnotTheory's convention).  So the under-strand runs
//! `a → c` and the over-strand joins `b` and `d`; the crossing is positive
//! when the over-strand enters at `d` (`d → b`) and negative when it enters
//! at `b`.  Every label occurs exactly twice.  The orientation of the
//! over-strands is read off by following each component; a component that
//! only ever passes over gets the orientation of its first label's first
//! occurrence as incoming.
//!
//! [`Diagram::from_pd`] renumbers the edges `0, 1, …` consecutively along
//! each component in the direction of travel; [`Diagram::pd_code`] prints
//! them 1-based, so a code in KnotTheory's numbering comes back unchanged.
//! Crossingless unknotted components (`free_loops`) are carried along so that
//! braid closures with idle strands and unlinks are diagrams too.

use crate::error::{KnotError, Result};
use std::collections::HashMap;

/// A corner of the diagram: crossing `c`, between slot `k` and slot `k + 1`
/// (mod 4) counterclockwise.
pub type Corner = (usize, u8);

/// An oriented link diagram; see the [module docs](self).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Diagram {
    x: Vec<[usize; 4]>,
    over_in: Vec<u8>,
    /// `ends[e] = [tail, head]`: the slot the edge leaves and the slot it enters.
    ends: Vec<[(usize, u8); 2]>,
    edge_comp: Vec<usize>,
    n_comp: usize,
    free_loops: usize,
}

impl Diagram {
    /// The crossingless unknot.
    pub fn unknot() -> Self {
        Self::unlink(1)
    }

    /// The crossingless unlink of `k` components.
    pub fn unlink(k: usize) -> Self {
        Diagram {
            x: Vec::new(),
            over_in: Vec::new(),
            ends: Vec::new(),
            edge_comp: Vec::new(),
            n_comp: 0,
            free_loops: k,
        }
    }

    /// The diagram with these crossings; see the [module docs](self) for the
    /// convention.  Labels are arbitrary integers, each used exactly twice.
    pub fn from_pd(code: &[[i64; 4]]) -> Result<Self> {
        Self::build(code, 0)
    }

    /// Parses `X[1,4,2,5], X[3,6,4,1], X[5,2,6,3]` (an optional `PD[...]`
    /// wrapper is ignored) or a bare list of 4-tuples like `[[1,4,2,5], ...]`.
    pub fn parse_pd(s: &str) -> Result<Self> {
        let nums: Vec<i64> = s
            .split(|ch: char| !(ch.is_ascii_digit() || ch == '-'))
            .filter(|t| !t.is_empty())
            .map(|t| {
                t.parse::<i64>()
                    .map_err(|e| KnotError::Invalid(format!("{t}: {e}")))
            })
            .collect::<Result<_>>()?;
        if nums.len() % 4 != 0 {
            return Err(KnotError::Invalid(format!(
                "{} labels is not a multiple of 4",
                nums.len()
            )));
        }
        let code: Vec<[i64; 4]> = nums.chunks(4).map(|c| [c[0], c[1], c[2], c[3]]).collect();
        Self::from_pd(&code)
    }

    pub(crate) fn build(code: &[[i64; 4]], free_loops: usize) -> Result<Self> {
        let n = code.len();
        if n == 0 {
            return Ok(Self::unlink(free_loops));
        }
        // labels -> provisional ids, in order of first appearance
        let mut id: HashMap<i64, usize> = HashMap::new();
        let mut order: Vec<i64> = Vec::new();
        let mut halves: Vec<Vec<(usize, u8)>> = Vec::new();
        for (c, x) in code.iter().enumerate() {
            for (s, &l) in x.iter().enumerate() {
                let e = *id.entry(l).or_insert_with(|| {
                    order.push(l);
                    halves.push(Vec::new());
                    halves.len() - 1
                });
                halves[e].push((c, s as u8));
            }
        }
        if let Some(e) = halves.iter().position(|h| h.len() != 2) {
            return Err(KnotError::Invalid(format!(
                "label {} occurs {} times, not twice",
                order[e],
                halves[e].len()
            )));
        }
        let hid = |(c, s): (usize, u8)| 4 * c + s as usize;
        let mut partner = vec![0usize; 4 * n];
        for h in &halves {
            partner[hid(h[0])] = hid(h[1]);
            partner[hid(h[1])] = hid(h[0]);
        }
        // orientation: incoming[h]; slot 0 in, slot 2 out; an edge has one end
        // in and one out; a strand enters at s and leaves at s + 2
        let mut incoming: Vec<Option<bool>> = vec![None; 4 * n];
        let mut stack: Vec<(usize, bool)> = Vec::new();
        for c in 0..n {
            stack.push((4 * c, true));
            stack.push((4 * c + 2, false));
        }
        let through = |h: usize| 4 * (h / 4) + (h % 4 + 2) % 4;
        let mut seeds = (0..4 * n).filter(|h| h % 2 == 1);
        loop {
            while let Some((h, v)) = stack.pop() {
                match incoming[h] {
                    Some(w) if w == v => continue,
                    Some(_) => {
                        return Err(KnotError::Invalid(format!(
                            "inconsistent orientation at crossing {} (slot {}): \
                             the first label of each crossing must be the incoming under-strand",
                            h / 4,
                            h % 4
                        )))
                    }
                    None => {
                        incoming[h] = Some(v);
                        stack.push((partner[h], !v));
                        stack.push((through(h), !v));
                    }
                }
            }
            // a component that only passes over: orient it by its first slot
            match seeds.by_ref().find(|&h| incoming[h].is_none()) {
                Some(h) => stack.push((h, true)),
                None => break,
            }
        }
        // renumber edges along components, starting from labels in first-appearance order
        let mut new_id: Vec<Option<usize>> = vec![None; halves.len()];
        let mut ends: Vec<[(usize, u8); 2]> = Vec::new();
        let mut edge_comp = Vec::new();
        let mut n_comp = 0;
        let half_of = |h: usize| (h / 4, (h % 4) as u8);
        let edge_at = |h: usize| id[&code[h / 4][h % 4]];
        for start in 0..halves.len() {
            if new_id[start].is_some() {
                continue;
            }
            let mut e = start;
            while new_id[e].is_none() {
                let [a, b] = [hid(halves[e][0]), hid(halves[e][1])];
                let (tail, head) = if incoming[a] == Some(true) { (b, a) } else { (a, b) };
                new_id[e] = Some(ends.len());
                ends.push([half_of(tail), half_of(head)]);
                edge_comp.push(n_comp);
                e = edge_at(through(head));
            }
            n_comp += 1;
        }
        let x: Vec<[usize; 4]> = code.iter().map(|x| x.map(|l| new_id[id[&l]].unwrap())).collect();
        let over_in = (0..n)
            .map(|c| if incoming[4 * c + 1] == Some(true) { 1 } else { 3 })
            .collect();
        let d = Diagram {
            x,
            over_in,
            ends,
            edge_comp,
            n_comp,
            free_loops,
        };
        // Euler: each piece of a planar projection, its faces traced on its own
        // sphere, has crossings + 2 faces
        if d.faces().len() != n + 2 * d.projection_pieces() {
            return Err(KnotError::Invalid(
                "the code is not planar (a virtual diagram?)".into(),
            ));
        }
        Ok(d)
    }

    /// The PD code with edges numbered `1, 2, …` along the components.
    pub fn pd_code(&self) -> Vec<[usize; 4]> {
        self.x.iter().map(|x| x.map(|e| e + 1)).collect()
    }

    /// Number of crossings.
    pub fn crossings(&self) -> usize {
        self.x.len()
    }

    /// Number of edges (`2 ×` crossings).
    pub fn edges(&self) -> usize {
        self.ends.len()
    }

    /// Number of link components, crossingless loops included.
    pub fn components(&self) -> usize {
        self.n_comp + self.free_loops
    }

    /// Number of crossingless unknotted components.
    pub fn free_loops(&self) -> usize {
        self.free_loops
    }

    /// The edges of crossing `c`, counterclockwise from the incoming under-strand.
    pub fn crossing(&self, c: usize) -> [usize; 4] {
        self.x[c]
    }

    /// The slot (1 or 3) at which the over-strand enters crossing `c`.
    pub fn over_in_slot(&self, c: usize) -> u8 {
        self.over_in[c]
    }

    /// The sign `±1` of crossing `c`.
    pub fn sign(&self, c: usize) -> i64 {
        if self.over_in[c] == 3 {
            1
        } else {
            -1
        }
    }

    /// Numbers of positive and negative crossings `(n₊, n₋)`.
    pub fn sign_counts(&self) -> (usize, usize) {
        let p = (0..self.crossings()).filter(|&c| self.sign(c) > 0).count();
        (p, self.crossings() - p)
    }

    /// The writhe: the sum of the crossing signs.
    pub fn writhe(&self) -> i64 {
        (0..self.crossings()).map(|c| self.sign(c)).sum()
    }

    /// The component (`0..components − free_loops`) that edge `e` belongs to.
    pub fn edge_component(&self, e: usize) -> usize {
        self.edge_comp[e]
    }

    /// `[tail, head]` of edge `e`: the (crossing, slot) it leaves and enters.
    pub fn edge_ends(&self, e: usize) -> [(usize, u8); 2] {
        self.ends[e]
    }

    /// Linking numbers `lk(i, j)` of the components (free loops last, with
    /// zero linking); the diagonal is zero.
    pub fn linking_matrix(&self) -> Vec<Vec<i64>> {
        let k = self.components();
        let mut m = vec![vec![0i64; k]; k];
        for c in 0..self.crossings() {
            let (i, j) = (self.edge_comp[self.x[c][0]], self.edge_comp[self.x[c][1]]);
            if i != j {
                m[i][j] += self.sign(c);
                m[j][i] += self.sign(c);
            }
        }
        for row in &mut m {
            for v in row.iter_mut() {
                debug_assert_eq!(*v % 2, 0);
                *v /= 2;
            }
        }
        m
    }

    fn code_i64(&self, f: impl Fn(&[usize; 4], u8) -> [usize; 4]) -> Vec<[i64; 4]> {
        (0..self.crossings())
            .map(|c| f(&self.x[c], self.over_in[c]).map(|e| e as i64))
            .collect()
    }

    /// The mirror image (every crossing switched).
    pub fn mirror(&self) -> Self {
        // the over-strand becomes the under-strand; start at its incoming slot
        let code = self.code_i64(|x, o| {
            let o = o as usize;
            [x[o], x[(o + 1) % 4], x[(o + 2) % 4], x[(o + 3) % 4]]
        });
        Self::build(&code, self.free_loops).expect("mirror of a valid diagram")
    }

    /// The same diagram with every component's orientation reversed.
    pub fn reverse(&self) -> Self {
        let code = self.code_i64(|x, _| [x[2], x[3], x[0], x[1]]);
        Self::build(&code, self.free_loops).expect("reverse of a valid diagram")
    }

    /// The split union of two diagrams.
    pub fn disjoint_union(&self, other: &Self) -> Self {
        let off = self.edges();
        let mut code = self.code_i64(|x, _| *x);
        code.extend(other.code_i64(|x, _| x.map(|e| e + off)));
        Self::build(&code, self.free_loops + other.free_loops).expect("union of valid diagrams")
    }

    /// The (crossing, slot) at the other end of the edge leaving `(c, s)`.
    pub(crate) fn partner(&self, c: usize, s: u8) -> (usize, u8) {
        let [t, h] = self.ends[self.x[c][s as usize]];
        if t == (c, s) {
            h
        } else {
            t
        }
    }

    /// The faces of the projection (crossingless loops excluded), each as
    /// its cycle of corners.  A face keeps the corner `(c, k)` on its left
    /// when it enters crossing `c` through slot `k + 1` and leaves through
    /// slot `k`.
    pub fn faces(&self) -> Vec<Vec<Corner>> {
        let n = self.crossings();
        let mut seen = vec![false; 4 * n];
        let mut faces = Vec::new();
        for start in 0..4 * n {
            if seen[start] {
                continue;
            }
            let mut face = Vec::new();
            let mut h = start;
            while !seen[h] {
                seen[h] = true;
                let (c, k) = (h / 4, (h % 4) as u8);
                face.push((c, k));
                // leave c through slot k, arrive at c2 through slot s2, turn left
                let (c2, s2) = self.partner(c, k);
                h = 4 * c2 + ((s2 + 3) % 4) as usize;
            }
            faces.push(face);
        }
        faces
    }

    /// The number of connected pieces of the projection, crossingless loops
    /// not counted.
    pub fn projection_pieces(&self) -> usize {
        let n = self.crossings();
        let mut seen = vec![false; n];
        let mut pieces = 0;
        for start in 0..n {
            if seen[start] {
                continue;
            }
            pieces += 1;
            seen[start] = true;
            let mut stack = vec![start];
            while let Some(c) = stack.pop() {
                for s in 0..4 {
                    let (c2, _) = self.partner(c, s);
                    if !seen[c2] {
                        seen[c2] = true;
                        stack.push(c2);
                    }
                }
            }
        }
        pieces
    }

    /// Is the projection (without its crossingless loops) connected?
    pub fn is_connected_projection(&self) -> bool {
        self.projection_pieces() <= 1
    }
}
