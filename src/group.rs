//! Knot groups: the Wirtinger presentation, Fox calculus and
//! homomorphisms to symmetric groups.
//!
//! The Wirtinger generators are the **arcs** of the diagram (maximal
//! over-strands, from one undercrossing to the next) plus one generator per
//! crossingless loop.  A crossing with over-arc `k`, incoming under-arc `i`,
//! outgoing under-arc `o` and sign `ε` gives the relator
//! `x_k^ε x_i x_k^{−ε} x_o^{−1}`.

use crate::bracket::IntLaurent;
use crate::diagram::Diagram;
use num_bigint::BigInt;

/// A finitely presented group.  A letter `g > 0` is generator `g − 1`, a
/// letter `−g` its inverse.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Presentation {
    /// Number of generators.
    pub generators: usize,
    /// The relators, as words in the letters `±1, …, ±generators`.
    pub relators: Vec<Vec<i32>>,
}

impl Presentation {
    /// The Fox Jacobian `(∂r_i/∂x_j)` pushed to `ℤ[t^{±1}]` by `x_j ↦ t^{degree[j]}`.
    pub fn fox_jacobian(&self, degree: &[i64]) -> Vec<Vec<IntLaurent>> {
        assert_eq!(degree.len(), self.generators);
        self.relators
            .iter()
            .map(|r| {
                let mut row = vec![IntLaurent::zero(); self.generators];
                let mut prefix = 0i64; // the image of the prefix is t^prefix
                for &l in r {
                    let g = l.unsigned_abs() as usize - 1;
                    if l > 0 {
                        // ∂(u x)/∂x = ∂u/∂x + u
                        row[g] = row[g].add_ref(&IntLaurent::t_pow(prefix));
                        prefix += degree[g];
                    } else {
                        // ∂(u x^{-1})/∂x = ∂u/∂x − u x^{-1}
                        prefix -= degree[g];
                        row[g] = row[g].sub_ref(&IntLaurent::t_pow(prefix));
                    }
                }
                row
            })
            .collect()
    }
}

/// The group of a link diagram with its Wirtinger structure.
#[derive(Clone, Debug)]
pub struct KnotGroup {
    presentation: Presentation,
    /// `(over, in, out, sign)` per crossing.
    wirtinger: Vec<(usize, usize, usize, i64)>,
    /// The link component of each generator.
    component: Vec<usize>,
}

/// The arc of every edge: edges `b` and `d` of a crossing lie on one arc.
pub(crate) fn arcs(d: &Diagram) -> (Vec<usize>, usize) {
    let m = d.edges();
    let mut parent: Vec<usize> = (0..m).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut y = x;
        while p[y] != r {
            let nxt = p[y];
            p[y] = r;
            y = nxt;
        }
        r
    }
    for c in 0..d.crossings() {
        let x = d.crossing(c);
        let (a, b) = (find(&mut parent, x[1]), find(&mut parent, x[3]));
        parent[a] = b;
    }
    let mut id = vec![usize::MAX; m];
    let mut arc = vec![0; m];
    let mut n = 0;
    for e in 0..m {
        let r = find(&mut parent, e);
        if id[r] == usize::MAX {
            id[r] = n;
            n += 1;
        }
        arc[e] = id[r];
    }
    (arc, n)
}

impl KnotGroup {
    /// The Wirtinger presentation of the diagram's link group.
    pub fn wirtinger(d: &Diagram) -> Self {
        let (arc, n_arcs) = arcs(d);
        let mut component = vec![0; n_arcs];
        for e in 0..d.edges() {
            component[arc[e]] = d.edge_component(e);
        }
        let n_crossing_comps = d.components() - d.free_loops();
        component.extend((0..d.free_loops()).map(|j| n_crossing_comps + j));
        let wirtinger: Vec<_> = (0..d.crossings())
            .map(|c| {
                let x = d.crossing(c);
                (arc[x[1]], arc[x[0]], arc[x[2]], d.sign(c))
            })
            .collect();
        let letter = |g: usize, e: i64| (g as i32 + 1) * e as i32;
        let relators = wirtinger
            .iter()
            .map(|&(k, i, o, s)| vec![letter(k, s), letter(i, 1), letter(k, -s), letter(o, -1)])
            .collect();
        KnotGroup {
            presentation: Presentation {
                generators: n_arcs + d.free_loops(),
                relators,
            },
            wirtinger,
            component,
        }
    }

    /// The presentation.
    pub fn presentation(&self) -> &Presentation {
        &self.presentation
    }

    /// The Alexander matrix: the Fox Jacobian with every meridian `↦ t`.
    pub fn alexander_matrix(&self) -> Vec<Vec<IntLaurent>> {
        self.presentation
            .fox_jacobian(&vec![1; self.presentation.generators])
    }

    /// The number of homomorphisms to the symmetric group `S_m`, `m ≤ 10`
    /// (`transitive_only`: only those whose image acts transitively on
    /// `{0, …, m−1}`).  A backtracking search: the images of a crossing's
    /// over-arc and one under-arc determine the other under-arc, the search
    /// branches on the generator that unlocks most relations, generators of
    /// one link component run over one conjugacy class, and the first
    /// generator only over class representatives.
    pub fn count_homomorphisms_to_symmetric(&self, m: usize, transitive_only: bool) -> BigInt {
        self.count_homomorphisms_with(&[], m, transitive_only)
    }

    /// Relators that turn the link group into the fundamental group of the
    /// 3-manifold from surgery on the link with these framings (one per
    /// component, free loops last): `λ_c μ_c^{f_c}` for each component `c`,
    /// with `λ_c` the longitude of linking number 0 with `c`.  Going along
    /// `c` from the arc `μ_c`, passing under arcs `x_{k_1}, …, x_{k_N}` at
    /// crossings of signs `ε_j`, gives the peripheral element
    /// `x_{k_N}^{ε_N} ⋯ x_{k_1}^{ε_1}`; it equals `λ_c μ_c^{w_c}`, `w_c` the
    /// writhe of `c`.
    pub fn surgery_relators(d: &Diagram, framings: &[i64]) -> Vec<Vec<i32>> {
        assert_eq!(framings.len(), d.components());
        let (arc, n_arcs) = arcs(d);
        let n_cross = d.components() - d.free_loops();
        let mut out = Vec::new();
        let power = |g: usize, k: i64| -> Vec<i32> {
            let l = g as i32 + 1;
            std::iter::repeat_n(if k > 0 { l } else { -l }, k.unsigned_abs() as usize).collect()
        };
        for c in 0..n_cross {
            let edges: Vec<usize> = (0..d.edges()).filter(|&e| d.edge_component(e) == c).collect();
            let mut word = Vec::new();
            for &e in &edges {
                let (x, slot) = d.edge_ends(e)[1];
                if slot == 0 {
                    let over = arc[d.crossing(x)[1]] as i32 + 1;
                    word.push(over * d.sign(x) as i32);
                }
            }
            word.reverse();
            let writhe: i64 = (0..d.crossings())
                .filter(|&x| {
                    let k = d.crossing(x);
                    d.edge_component(k[0]) == c && d.edge_component(k[1]) == c
                })
                .map(|x| d.sign(x))
                .sum();
            word.extend(power(arc[edges[0]], framings[c] - writhe));
            out.push(word);
        }
        for j in 0..d.free_loops() {
            out.push(power(n_arcs + j, framings[n_cross + j]));
        }
        out
    }

    /// The number of homomorphisms to `S_m` that also kill the words in
    /// `extra` (letters as in [`Presentation`]); see
    /// [`count_homomorphisms_to_symmetric`](Self::count_homomorphisms_to_symmetric).
    /// With [`surgery_relators`](Self::surgery_relators) this counts
    /// homomorphisms from the fundamental group of a surgered 3-manifold.
    pub fn count_homomorphisms_with(&self, extra: &[Vec<i32>], m: usize, transitive_only: bool) -> BigInt {
        assert!((1..=MAX_DEGREE).contains(&m), "S_m with 1 ≤ m ≤ {MAX_DEGREE}");
        let n = self.presentation.generators;
        if n == 0 {
            return BigInt::from(1);
        }
        let mut touching = vec![Vec::new(); n];
        for (t, &(k, i, o, _)) in self.wirtinger.iter().enumerate() {
            for g in [k, i, o] {
                if !touching[g].contains(&t) {
                    touching[g].push(t);
                }
            }
        }
        let all = all_permutations(m);
        let mut search = Search {
            group: self,
            touching,
            img: vec![None; n],
            trail: Vec::new(),
            m,
            transitive_only,
            extra,
            classes: all.iter().map(|p| cycle_type(p, m)).collect(),
            all,
        };
        let mut total = BigInt::from(0);
        for (rep, size) in class_representatives(m) {
            let count = if search.assign(0, rep) { search.run() } else { 0 };
            search.undo(0);
            total += BigInt::from(count) * size;
        }
        total
    }
}

/// The largest `m` for [`KnotGroup::count_homomorphisms_to_symmetric`].
pub const MAX_DEGREE: usize = 10;

/// Backtracking state: the images so far and the order they were set in.
struct Search<'a> {
    group: &'a KnotGroup,
    touching: Vec<Vec<usize>>,
    img: Vec<Option<Perm>>,
    trail: Vec<usize>,
    m: usize,
    transitive_only: bool,
    extra: &'a [Vec<i32>],
    all: Vec<Perm>,
    classes: Vec<Vec<usize>>,
}

impl Search<'_> {
    /// Sets generator `g` to `p` and propagates; `false` on a contradiction
    /// (the trail still records what was set, for [`undo`](Self::undo)).
    fn assign(&mut self, g: usize, p: Perm) -> bool {
        self.img[g] = Some(p);
        self.trail.push(g);
        let mut queue = vec![g];
        while let Some(g) = queue.pop() {
            for ti in 0..self.touching[g].len() {
                let (k, i, o, s) = self.group.wirtinger[self.touching[g][ti]];
                let Some(pk) = self.img[k] else { continue };
                // x_o = x_k^s x_i x_k^{-s}
                match (self.img[i], self.img[o]) {
                    (Some(pi), Some(po)) => {
                        if conjugate(&pk, &pi, s, self.m) != po {
                            return false;
                        }
                    }
                    (Some(pi), None) => {
                        self.img[o] = Some(conjugate(&pk, &pi, s, self.m));
                        self.trail.push(o);
                        queue.push(o);
                    }
                    (None, Some(po)) => {
                        self.img[i] = Some(conjugate(&pk, &po, -s, self.m));
                        self.trail.push(i);
                        queue.push(i);
                    }
                    (None, None) => {}
                }
            }
        }
        true
    }

    fn word_is_identity(&self, w: &[i32]) -> bool {
        let mut p = identity();
        for &l in w {
            let x = self.img[l.unsigned_abs() as usize - 1].unwrap();
            let x = if l > 0 { x } else { inverse(&x, self.m) };
            // p ∘ x
            let mut q = identity();
            for i in 0..self.m {
                q[i] = p[x[i] as usize];
            }
            p = q;
        }
        p == identity()
    }

    fn undo(&mut self, mark: usize) {
        for g in self.trail.drain(mark..) {
            self.img[g] = None;
        }
    }

    /// The number of completions of the current partial assignment.
    fn run(&mut self) -> u64 {
        let n = self.img.len();
        // branch on the unknown generator that is the over-arc of the most
        // relations with exactly one known under-arc
        let mut best: Option<(usize, usize)> = None;
        for g in (0..n).filter(|&g| self.img[g].is_none()) {
            let score = self.touching[g]
                .iter()
                .filter(|&&t| {
                    let (k, i, o, _) = self.group.wirtinger[t];
                    k == g && (self.img[i].is_some() != self.img[o].is_some())
                })
                .count();
            if best.is_none_or(|(_, b)| score > b) {
                best = Some((g, score));
            }
        }
        let Some((g, _)) = best else {
            let ok = (!self.transitive_only || is_transitive(&self.img, self.m))
                && self.extra.iter().all(|w| self.word_is_identity(w));
            return u64::from(ok);
        };
        // generators of one component are conjugate: same cycle type
        let class = (0..n)
            .find(|&h| self.group.component[h] == self.group.component[g] && self.img[h].is_some())
            .map(|h| cycle_type(&self.img[h].unwrap(), self.m));
        let mut count = 0;
        for j in 0..self.all.len() {
            if class.as_ref().is_some_and(|c| *c != self.classes[j]) {
                continue;
            }
            let mark = self.trail.len();
            if self.assign(g, self.all[j]) {
                count += self.run();
            }
            self.undo(mark);
        }
        count
    }
}

/// A permutation of `0..m` (entries from `m` on are fixed points): `p[i]`
/// is the image of `i`.
type Perm = [u8; MAX_DEGREE];

fn identity() -> Perm {
    std::array::from_fn(|i| i as u8)
}

/// `k x k^{-1}` (`s = 1`) or `k^{-1} x k` (`s = −1`), maps composed right to left.
fn conjugate(k: &Perm, x: &Perm, s: i64, m: usize) -> Perm {
    let mut out = identity();
    if s > 0 {
        // k x k^{-1} sends k(i) to k(x(i))
        for i in 0..m {
            out[k[i] as usize] = k[x[i] as usize];
        }
    } else {
        // k^{-1} x k sends i to k^{-1}(x(k(i)))
        let mut kinv = identity();
        for j in 0..m {
            kinv[k[j] as usize] = j as u8;
        }
        for i in 0..m {
            out[i] = kinv[x[k[i] as usize] as usize];
        }
    }
    out
}

fn inverse(p: &Perm, m: usize) -> Perm {
    let mut q = identity();
    for i in 0..m {
        q[p[i] as usize] = i as u8;
    }
    q
}

fn cycle_type(p: &Perm, m: usize) -> Vec<usize> {
    let mut seen = [false; MAX_DEGREE];
    let mut t = Vec::new();
    for s in 0..m {
        if seen[s] {
            continue;
        }
        let (mut len, mut i) = (0, s);
        while !seen[i] {
            seen[i] = true;
            i = p[i] as usize;
            len += 1;
        }
        t.push(len);
    }
    t.sort_unstable();
    t
}

fn all_permutations(m: usize) -> Vec<Perm> {
    fn rec(p: &mut Perm, k: usize, m: usize, out: &mut Vec<Perm>) {
        if k == m {
            out.push(*p);
            return;
        }
        for j in k..m {
            p.swap(k, j);
            rec(p, k + 1, m, out);
            p.swap(k, j);
        }
    }
    let mut out = Vec::new();
    rec(&mut identity(), 0, m, &mut out);
    out
}

/// One permutation per cycle type, with the size of its conjugacy class.
fn class_representatives(m: usize) -> Vec<(Perm, BigInt)> {
    fn partitions(n: usize, max: usize) -> Vec<Vec<usize>> {
        if n == 0 {
            return vec![vec![]];
        }
        (1..=max.min(n))
            .rev()
            .flat_map(|k| {
                partitions(n - k, k).into_iter().map(move |mut p| {
                    p.insert(0, k);
                    p
                })
            })
            .collect()
    }
    let fact = |n: usize| (1..=n).fold(BigInt::from(1), |a, k| a * k);
    partitions(m, m)
        .into_iter()
        .map(|parts| {
            let mut p = identity();
            let mut start = 0u8;
            for &k in &parts {
                for j in 0..k as u8 {
                    p[(start + j) as usize] = start + (j + 1) % k as u8;
                }
                start += k as u8;
            }
            // |class| = m! / Π_k k^{a_k} a_k!
            let mut denom = BigInt::from(1);
            for k in 1..=m {
                let a = parts.iter().filter(|&&x| x == k).count();
                denom *= BigInt::from(k).pow(a as u32) * fact(a);
            }
            (p, fact(m) / denom)
        })
        .collect()
}

fn is_transitive(img: &[Option<Perm>], m: usize) -> bool {
    let mut seen = vec![false; m];
    let mut stack = vec![0usize];
    seen[0] = true;
    while let Some(i) = stack.pop() {
        for p in img.iter().flatten() {
            let j = p[i] as usize;
            if !seen[j] {
                seen[j] = true;
                stack.push(j);
            }
        }
    }
    seen.iter().all(|&s| s)
}
