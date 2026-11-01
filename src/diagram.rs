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

    /// The Seifert circle of each edge, numbered `0, 1, …`, and the number
    /// of circles (crossingless loops count, after the others).  The
    /// oriented smoothing of `X[a, b, c, d]` joins each incoming end to the
    /// outgoing end next to it: `a–b`, `d–c` if positive, `a–d`, `b–c` if
    /// negative.
    pub fn seifert_circles(&self) -> (Vec<usize>, usize) {
        let m = self.edges();
        let mut parent: Vec<usize> = (0..m).collect();
        fn find(p: &mut [usize], x: usize) -> usize {
            let mut r = x;
            while p[r] != r {
                r = p[r];
            }
            p[x] = r;
            r
        }
        for c in 0..self.crossings() {
            let [a, b, cc, d] = self.x[c];
            let pairs = if self.sign(c) > 0 {
                [(a, b), (d, cc)]
            } else {
                [(a, d), (b, cc)]
            };
            for (u, v) in pairs {
                let (ru, rv) = (find(&mut parent, u), find(&mut parent, v));
                parent[ru] = rv;
            }
        }
        let mut id = vec![usize::MAX; m];
        let mut circle = vec![0; m];
        let mut k = 0;
        for e in 0..m {
            let r = find(&mut parent, e);
            if id[r] == usize::MAX {
                id[r] = k;
                k += 1;
            }
            circle[e] = id[r];
        }
        (circle, k + self.free_loops)
    }

    /// The number of connected components of the Seifert graph (a vertex per
    /// Seifert circle, an edge per crossing) keeping only the crossings of
    /// sign `sign` (`None`: all of them).  Isolated circles count.
    pub fn seifert_graph_components(&self, sign: Option<i64>) -> usize {
        let (circle, k) = self.seifert_circles();
        let mut parent: Vec<usize> = (0..k).collect();
        fn find(p: &mut [usize], x: usize) -> usize {
            let mut r = x;
            while p[r] != r {
                r = p[r];
            }
            p[x] = r;
            r
        }
        let mut comps = k;
        for c in 0..self.crossings() {
            if sign.is_some_and(|s| s != self.sign(c)) {
                continue;
            }
            // the two circles meeting at a crossing hold its two incoming ends
            let x = self.x[c];
            let other_in = if self.sign(c) > 0 { x[3] } else { x[1] };
            let (u, v) = (
                find(&mut parent, circle[x[0]]),
                find(&mut parent, circle[other_in]),
            );
            if u != v {
                parent[u] = v;
                comps -= 1;
            }
        }
        comps
    }

    /// A diagram of the same link with fewer crossings.  Removes kinks
    /// (Reidemeister I) and bigons where one strand passes over the other
    /// twice (Reidemeister II), and reroutes strands that pass only over or
    /// only under along shorter paths ([`Self::pickup`]); then searches for
    /// more of these through Reidemeister III moves, made at random from a
    /// fixed seed (so the result is reproducible), keeping the smallest
    /// diagram met.  The search
    /// stops after `20 n` moves (`n` crossings, at least 10) without a new
    /// smallest diagram.  Components that lose all their crossings become
    /// crossingless loops.
    pub fn simplify(&self) -> Diagram {
        let mut cur = self.reduce();
        let mut best = cur.clone();
        let mut rng: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut stale = 0;
        while stale < 20 * best.crossings().max(10) {
            let moves = cur.r3_faces();
            if moves.is_empty() {
                break;
            }
            // xorshift64
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            stale += 1;
            let Some(next) = cur.r3(&moves[(rng % moves.len() as u64) as usize]) else {
                continue;
            };
            cur = next.reduce();
            if cur.crossings() < best.crossings() {
                best = cur.clone();
                stale = 0;
            }
        }
        best
    }

    /// The triangular faces at which a Reidemeister III move can be made:
    /// three distinct crossings, and strands not cyclically over one another.
    pub fn r3_faces(&self) -> Vec<Vec<Corner>> {
        self.faces()
            .into_iter()
            .filter(|f| {
                f.len() == 3 && f[0].0 != f[1].0 && f[1].0 != f[2].0 && f[0].0 != f[2].0 && !self.r3_cyclic(f)
            })
            .collect()
    }

    /// The Reidemeister III move at a triangular face (from [`Self::faces`]),
    /// or `None` if the face is not one.
    ///
    /// Number the face's sides `e₁, e₂, e₃` in the order the face lists its
    /// corners, and the strands along them 1, 2, 3; crossing `cᵢ` joins
    /// sides `eᵢ₋₁` and `eᵢ`.  Strand `i` runs from its outer end `Pᵢ` at
    /// `cᵢ` to its outer end `Qᵢ` at `cᵢ₊₁`.  Counterclockwise around the
    /// triangle the outer ends are `P₁ Q₃ | P₂ Q₁ | P₃ Q₂` (two at each
    /// crossing).  The move slides one strand across the opposite crossing:
    /// afterwards strands `i` and `j` cross next to the two outer ends that
    /// neither of them had at their old crossing, over and under as before,
    /// so the order of the crossings along each strand is reversed.
    pub fn r3(&self, face: &[Corner]) -> Option<Diagram> {
        let &[(c1, k1), (c2, k2), (c3, k3)] = face else {
            return None;
        };
        if c1 == c2 || c2 == c3 || c1 == c3 {
            return None;
        }
        let slot = |c: usize, k: u8, d: u8| self.x[c][((k + d) % 4) as usize] as i64;
        let (e1, e2, e3) = (slot(c1, k1, 0), slot(c2, k2, 0), slot(c3, k3, 0));
        let (p1, q3) = (slot(c1, k1, 2), slot(c1, k1, 3));
        let (p2, q1) = (slot(c2, k2, 2), slot(c2, k2, 3));
        let (p3, q2) = (slot(c3, k3, 2), slot(c3, k3, 3));
        // over/under: a strand is over where it sits in odd slots; side eᵢ
        // leaves cᵢ by slot kᵢ and enters cᵢ₊₁ by slot kᵢ₊₁ + 1
        let o13 = k1 % 2 == 1; // strand 1 over strand 3, at c1
        let o12 = (k2 + 1) % 2 == 1; // strand 1 over strand 2, at c2
        let o23 = (k3 + 1) % 2 == 1; // strand 2 over strand 3, at c3
        if self.r3_cyclic(face) {
            return None;
        }
        // orientation: strand i runs P → Q when side eᵢ leaves cᵢ
        let fwd = |e: i64, c: usize, k: u8| self.ends[e as usize][0] == (c, k);
        let (f1, f2, f3) = (fwd(e1, c1, k1), fwd(e2, c2, k2), fwd(e3, c3, k3));
        // the new crossings: arms counterclockwise as (label, strand, incoming)
        let n13 = [(q1, 1, !f1), (p3, 3, f3), (e1, 1, f1), (e3, 3, !f3)];
        let n12 = [(q2, 2, !f2), (p1, 1, f1), (e2, 2, f2), (e1, 1, !f1)];
        let n23 = [(q3, 3, !f3), (p2, 2, f2), (e3, 3, f3), (e2, 2, !f2)];
        // KnotTheory's order: from the incoming under-strand
        let pd = |arms: [(i64, u8, bool); 4], under: u8| -> [i64; 4] {
            let s = (0..4).find(|&s| arms[s].1 == under && arms[s].2).unwrap();
            [0, 1, 2, 3].map(|i| arms[(s + i) % 4].0)
        };
        let mut code = self.code_i64(|x, _| *x);
        code[c1] = pd(n13, if o13 { 3 } else { 1 });
        code[c2] = pd(n12, if o12 { 2 } else { 1 });
        code[c3] = pd(n23, if o23 { 3 } else { 2 });
        Diagram::build(&code, self.free_loops).ok()
    }

    /// [`Self::simplify`], with backtracking to get out of local minima: from
    /// the smallest diagram so far, a few random moves that add crossings
    /// ([`Self::r2_up`]), then `simplify` again; the result is kept when it is
    /// smaller.  Stops after `rounds` rounds without a smaller diagram (10 is
    /// a good default).  Reproducible, like `simplify`.
    pub fn simplify_global(&self, rounds: usize) -> Diagram {
        let mut best = self.simplify();
        let mut rng: u64 = 0x2545_f491_4f6c_dd1d;
        let mut next = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };
        let mut stale = 0;
        while stale < rounds && best.crossings() > 0 {
            stale += 1;
            let mut d = best.clone();
            for _ in 0..1 + next() % 10 {
                let faces: Vec<Vec<Corner>> = d.faces().into_iter().filter(|f| f.len() >= 2).collect();
                let f = &faces[(next() % faces.len() as u64) as usize];
                let (i, j) = (
                    (next() % f.len() as u64) as usize,
                    (next() % f.len() as u64) as usize,
                );
                if let Some(e) = d.r2_up(f, i, j, next() % 2 == 0) {
                    d = e;
                }
            }
            let d = d.simplify();
            if d.crossings() < best.crossings() {
                best = d;
                stale = 0;
            }
        }
        best
    }

    /// The Reidemeister II move that adds two crossings: in a face (from
    /// [`Self::faces`]), a finger of side `i` is pushed across side `j`
    /// (`i ≠ j`), over it if `over` and under it otherwise.  Side `i` is the
    /// edge leaving the face's corner `i`.  `None` if `i = j` or the result
    /// is not a valid diagram.
    ///
    /// With the face on the left, side `e` runs `A → B` and side `f` runs
    /// `C → D`.  The finger leaves `e`, crosses `f` at `X`, turns and crosses
    /// back at `Y`; along `e` the pieces are `e_a, e_mid, e_b`, and along `f`
    /// (which meets `Y` first) `f_a, f_mid, f_b`.  Counterclockwise, `X` has
    /// `e_a, f_mid, e_mid, f_b` and `Y` has `e_b, f_a, e_mid, f_mid`.
    pub fn r2_up(&self, face: &[Corner], i: usize, j: usize, over: bool) -> Option<Diagram> {
        let n = face.len();
        if i == j || i >= n || j >= n {
            return None;
        }
        let (a, ka) = face[i];
        let (b, kb) = face[(i + 1) % n];
        let (c, kc) = face[j];
        let (d, kd) = face[(j + 1) % n];
        let e = self.x[a][ka as usize];
        let f = self.x[c][kc as usize];
        if e == f {
            return None;
        }
        // orientation: forward when the edge leaves its corner's crossing
        let ef = self.ends[e][0] == (a, ka);
        let ff = self.ends[f][0] == (c, kc);
        let mut code = self.code_i64(|x, _| *x);
        let top = 2 * self.crossings() as i64;
        let (e_a, e_mid, e_b) = (e as i64, top, top + 1);
        let (f_a, f_mid, f_b) = (f as i64, top + 2, top + 3);
        // the far ends of e and f take their new labels
        code[b][((kb + 1) % 4) as usize] = e_b;
        code[d][((kd + 1) % 4) as usize] = f_b;
        // arms counterclockwise as (label, on e, incoming)
        let xs = [
            (e_a, true, ef),
            (f_mid, false, ff),
            (e_mid, true, !ef),
            (f_b, false, !ff),
        ];
        let ys = [
            (e_b, true, !ef),
            (f_a, false, ff),
            (e_mid, true, ef),
            (f_mid, false, !ff),
        ];
        let pd = |arms: [(i64, bool, bool); 4]| -> [i64; 4] {
            // the under-strand is f when e passes over
            let s = (0..4).find(|&s| arms[s].1 != over && arms[s].2).unwrap();
            [0, 1, 2, 3].map(|t| arms[(s + t) % 4].0)
        };
        code.push(pd(xs));
        code.push(pd(ys));
        Diagram::build(&code, self.free_loops).ok()
    }

    /// Strand pickup: the over-arc (a run of a component between two
    /// under-passes, passing over everything in between) or under-arc whose
    /// rerouting saves the most crossings, lifted off and laid back along a
    /// shortest path through the faces of the rest; `None` if no arc can be
    /// shortened.  An arc that passes over (under) everything it meets can go
    /// anywhere, so the link is unchanged.
    pub fn pickup(&self) -> Option<Diagram> {
        let over = self.pickup_over();
        let under = self.mirror().pickup_over().map(|d| d.mirror());
        match (over, under) {
            (Some(a), Some(b)) => Some(if b.crossings() < a.crossings() { b } else { a }),
            (a, b) => a.or(b),
        }
    }

    /// [`Self::pickup`] for over-arcs only.
    fn pickup_over(&self) -> Option<Diagram> {
        let n = self.crossings();
        let m = self.edges();
        if n == 0 {
            return None;
        }
        // the faces on the left and right of each edge (in its direction)
        let faces = self.faces();
        let (mut left, mut right) = (vec![0; m], vec![0; m]);
        for (fi, f) in faces.iter().enumerate() {
            for &(c, k) in f {
                let l = self.x[c][k as usize];
                if self.ends[l][0] == (c, k) {
                    left[l] = fi;
                } else {
                    right[l] = fi;
                }
            }
        }
        let next_edge = |e: usize| {
            let (c, s) = self.ends[e][1];
            self.x[c][((s + 2) % 4) as usize]
        };
        let is_under_head = |e: usize| self.ends[e][1].1 % 2 == 0;
        // the over-arcs: edges a₀ … a_k, a₀ leaving an under-pass and a_k
        // entering the next one, over k crossings between
        let mut arcs: Vec<Vec<usize>> = Vec::new();
        let mut seen = vec![false; m];
        for e0 in 0..m {
            if seen[e0] || !is_under_head(e0) {
                continue;
            }
            // e0 ends at an under-pass: walk the component from there
            let mut e = e0;
            loop {
                seen[e] = true;
                let mut arc = vec![next_edge(e)];
                while !is_under_head(*arc.last().unwrap()) {
                    let l = next_edge(*arc.last().unwrap());
                    arc.push(l);
                }
                e = *arc.last().unwrap();
                if arc.len() >= 2 {
                    arcs.push(arc);
                }
                if e == e0 {
                    break;
                }
            }
        }
        let mut best: Option<(usize, Diagram)> = None;
        for arc in &arcs {
            let k = arc.len() - 1;
            if best.as_ref().is_some_and(|(saved, _)| *saved >= k) {
                continue;
            }
            if let Some((len, d)) = self.reroute(arc, &left, &right, faces.len()) {
                if len < k && best.as_ref().is_none_or(|(saved, _)| k - len > *saved) {
                    best = Some((k - len, d));
                }
            }
        }
        best.map(|(_, d)| d)
    }

    /// Lifts the over-arc `arc` (edges `a₀ … a_k`) and lays it back along a
    /// shortest path in the dual graph of the rest: the new length and the
    /// diagram.  `None` if that is no shorter, or lifting it leaves a
    /// crossingless loop.
    fn reroute(
        &self,
        arc: &[usize],
        left: &[usize],
        right: &[usize],
        nfaces: usize,
    ) -> Option<(usize, Diagram)> {
        fn find(p: &mut [usize], x: usize) -> usize {
            let mut r = x;
            while p[r] != r {
                r = p[r];
            }
            p[x] = r;
            r
        }
        let k = arc.len() - 1;
        let lifted: Vec<usize> = arc[..k].iter().map(|&a| self.ends[a][1].0).collect();
        let mut removed = vec![false; self.crossings()];
        let ends_at = [self.ends[arc[0]][0].0, self.ends[arc[k]][1].0];
        for &c in &lifted {
            // the arc passes over one crossing twice, or over its own end
            if removed[c] || ends_at.contains(&c) {
                return None;
            }
            removed[c] = true;
        }
        let in_arc: std::collections::HashSet<usize> = arc.iter().copied().collect();
        // faces of the rest: merged across the arc's edges
        let mut fp: Vec<usize> = (0..nfaces).collect();
        for &a in arc {
            let (u, v) = (find(&mut fp, left[a]), find(&mut fp, right[a]));
            fp[u] = v;
        }
        // the head of an edge once the lifted crossings are gone (following
        // the under-strand through them); None if it closes into a loop
        let head = |mut l: usize| -> Option<(usize, u8)> {
            for _ in 0..=self.crossings() {
                let (c, s) = self.ends[l][1];
                if !removed[c] {
                    return Some((c, s));
                }
                l = self.x[c][((s + 2) % 4) as usize];
            }
            None
        };
        let tail_label = |mut l: usize| -> Option<usize> {
            for _ in 0..=self.crossings() {
                let (c, s) = self.ends[l][0];
                if !removed[c] {
                    return Some(l);
                }
                l = self.x[c][((s + 2) % 4) as usize];
            }
            None
        };
        // breadth-first search over the merged faces, crossing edges of the rest
        let start = find(&mut fp, left[arc[0]]);
        let goal = find(&mut fp, left[arc[k]]);
        let mut adj: HashMap<usize, Vec<(usize, usize, bool)>> = HashMap::new(); // face -> (face, edge, left→right)
        for l in 0..self.edges() {
            if in_arc.contains(&l) {
                continue;
            }
            let (u, v) = (find(&mut fp, left[l]), find(&mut fp, right[l]));
            if u != v {
                adj.entry(u).or_default().push((v, l, true));
                adj.entry(v).or_default().push((u, l, false));
            }
        }
        let mut prev: HashMap<usize, (usize, usize, bool)> = HashMap::new();
        let mut queue = std::collections::VecDeque::from([start]);
        let mut reached = start == goal;
        while let Some(u) = queue.pop_front() {
            if reached {
                break;
            }
            for &(v, l, lr) in adj.get(&u).map(|v| v.as_slice()).unwrap_or(&[]) {
                if v != start && !prev.contains_key(&v) {
                    prev.insert(v, (u, l, lr));
                    if v == goal {
                        reached = true;
                        break;
                    }
                    queue.push_back(v);
                }
            }
        }
        if !reached {
            return None;
        }
        let mut path: Vec<(usize, bool)> = Vec::new(); // (edge crossed, left→right)
        let mut v = goal;
        while v != start {
            let (u, l, lr) = prev[&v];
            path.push((l, lr));
            v = u;
        }
        path.reverse();
        let len = path.len();
        if len >= k {
            return None;
        }
        // the new code: lifted crossings gone, their under-edges merged
        // (labelled by the piece leaving a kept crossing)
        let mut code: Vec<[i64; 4]> = Vec::new();
        let mut label = vec![0i64; self.edges()];
        for l in 0..self.edges() {
            label[l] = if in_arc.contains(&l) {
                l as i64
            } else {
                tail_label(l)? as i64
            };
        }
        let mut pos: HashMap<(usize, u8), (usize, usize)> = HashMap::new(); // old slot -> new slot
        for c in 0..self.crossings() {
            if removed[c] {
                continue;
            }
            for s in 0..4u8 {
                pos.insert((c, s), (code.len(), s as usize));
            }
            code.push(self.x[c].map(|l| label[l]));
        }
        let mut fresh = self.edges() as i64;
        let (a0, ak) = (arc[0] as i64, arc[k] as i64);
        if len == 0 {
            let (c, s) = self.ends[arc[k]][1];
            let (i, j) = pos[&(c, s)];
            code[i][j] = a0;
        }
        let mut piece = a0;
        for (j, &(l, lr)) in path.iter().enumerate() {
            let g = label[l];
            let (hc, hs) = head(l)?;
            let (i, t) = pos[&(hc, hs)];
            let g2 = fresh;
            fresh += 1;
            code[i][t] = g2;
            let out = if j + 1 == len { ak } else { fresh };
            if j + 1 < len {
                fresh += 1;
            }
            code.push(if lr {
                [g, out, g2, piece]
            } else {
                [g, piece, g2, out]
            });
            piece = out;
        }
        Diagram::build(&code, self.free_loops).ok().map(|d| (len, d))
    }

    /// Whether the strands of a triangular face pass cyclically over one
    /// another (strand 1 over 2 over 3 over 1, or the reverse), when no
    /// Reidemeister III move is possible.
    fn r3_cyclic(&self, face: &[Corner]) -> bool {
        let (k1, k2, k3) = (face[0].1, face[1].1, face[2].1);
        let (o13, o12, o23) = (k1 % 2 == 1, (k2 + 1) % 2 == 1, (k3 + 1) % 2 == 1);
        (o12 && o23 && !o13) || (!o12 && !o23 && o13)
    }

    /// Reidemeister I and II moves and strand pickups until none is left.
    fn reduce(&self) -> Diagram {
        let mut d = self.simplify_basic();
        while let Some(e) = d.pickup() {
            d = e.simplify_basic();
        }
        d
    }

    /// Removes kinks and bigons where one strand passes over the other twice
    /// (Reidemeister I and II) until there are none; the first step of
    /// [`Self::simplify`].
    pub fn simplify_basic(&self) -> Diagram {
        let mut code: Vec<[i64; 4]> = self.code_i64(|x, _| *x);
        let mut free = self.free_loops;
        loop {
            // Reidemeister I: an edge from a slot to the next one
            let kink = code.iter().position(|x| (0..4).any(|i| x[i] == x[(i + 1) % 4]));
            if let Some(c) = kink {
                let x = code.remove(c);
                let i = (0..4).find(|&i| x[i] == x[(i + 1) % 4]).unwrap();
                let (y, z) = (x[(i + 2) % 4], x[(i + 3) % 4]);
                if y == z {
                    free += 1;
                } else {
                    merge(&mut code, z, y);
                }
                continue;
            }
            let Ok(d) = Diagram::build(&code, free) else { break };
            // Reidemeister II: a bigon with one edge over at both crossings
            let mut done = false;
            for face in d.faces() {
                if face.len() != 2 || face[0].0 == face[1].0 {
                    continue;
                }
                let (c1, k1) = face[0];
                let c2 = face[1].0;
                let x1 = code[c1];
                let x2 = code[c2];
                let edges = [x1[k1 as usize], x1[(k1 as usize + 1) % 4]];
                let slot = |x: &[i64; 4], e: i64| (0..4).find(|&s| x[s] == e).unwrap();
                let over = edges
                    .iter()
                    .copied()
                    .find(|&e| slot(&x1, e) % 2 == 1 && slot(&x2, e) % 2 == 1);
                let Some(e_over) = over else { continue };
                let e_under = if edges[0] == e_over { edges[1] } else { edges[0] };
                if slot(&x1, e_under) % 2 != 0 || slot(&x2, e_under) % 2 != 0 {
                    continue;
                }
                // the far ends of the two strands
                let a = x1[(slot(&x1, e_over) + 2) % 4];
                let b = x2[(slot(&x2, e_over) + 2) % 4];
                let u1 = x1[(slot(&x1, e_under) + 2) % 4];
                let u2 = x2[(slot(&x2, e_under) + 2) % 4];
                let (hi, lo) = (c1.max(c2), c1.min(c2));
                code.remove(hi);
                code.remove(lo);
                for (p, q) in [(a, b), (u1, u2)] {
                    if p == q {
                        free += 1;
                    } else {
                        merge(&mut code, q, p);
                    }
                }
                done = true;
                break;
            }
            if !done {
                break;
            }
        }
        Diagram::build(&code, free).unwrap_or_else(|_| self.clone())
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

/// Renames edge label `from` to `to` (two edges joined into one).
fn merge(code: &mut [[i64; 4]], from: i64, to: i64) {
    for x in code.iter_mut() {
        for l in x.iter_mut() {
            if *l == from {
                *l = to;
            }
        }
    }
}
