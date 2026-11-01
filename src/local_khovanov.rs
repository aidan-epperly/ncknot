//! Khovanov homology and the Rasmussen invariant by Bar-Natan's local
//! algorithm ("Fast Khovanov homology computations", arXiv:math/0606318).
//!
//! The crossings are added one at a time.  After each one the complex is a
//! complex of *tangles*: its objects are crossingless matchings of the open
//! ends with a quantum shift, its maps linear combinations of dotted
//! cobordisms.  Closed circles are removed at once (*delooping*:
//! `○ ≅ ∅{+1} ⊕ ∅{−1}`) and every isomorphism in the differential is
//! cancelled (*Gaussian elimination*), so the complex stays small: its size
//! depends on how many ends are open, not on how many crossings have been
//! added.  The crossings are added in an order that keeps few ends open.
//!
//! The theory is Bar-Natan's over `R = F[h]` (`F` a field, `deg h = −2`),
//! with the Frobenius algebra `A = R[X]/(X² − hX)`: a cobordism with dots
//! reduces by `X² = hX` and neck-cutting.  A map between matchings `P → Q`
//! is an element of `A^{⊗ c}`, `c` the number of circles of `P ∪ Q`: a
//! decoration `1` or `X` on the disk each circle bounds, times a power of
//! `h`.  Composition and the planar gluing that adds a crossing glue such
//! disks along intervals; a glued component with `k` boundary circles,
//! genus `g` and `d` dots is `Δ^{(k)}(X^d (2X − h)^g)`.
//!
//! At the end the objects are empty matchings and the differential has
//! entries `c·h^k`, `k ≥ 1`.
//! * Setting `h = 0` gives Khovanov homology over `F`: the generators that
//!   are left ([`khovanov`]).
//! * Over `F[h]`, the homology of a knot in degree 0 has two free
//!   generators, in quantum degrees `s ± 1`, with `s` the Rasmussen
//!   invariant over `F` ([`rasmussen_s`]).  It comes from a graded Smith
//!   form, which for monomial entries is ordinary elimination with the pivot
//!   of least `h`-power first.
//!
//! Conventions are those of [`crate::khovanov`] (Bar-Natan's gradings):
//! the right-handed trefoil has `s = 2`.

use crate::diagram::Diagram;
use crate::error::{KnotError, Result};
use crate::khovanov::{Khovanov, Ring};
use ncalg::{BigRational, Field, Fp};
use smallvec::{smallvec, SmallVec};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

/// Coefficient fields: `Send + Sync` so large cancellations can run on
/// several threads.
pub(crate) trait Coef: Field + Send + Sync {}
impl<C: Field + Send + Sync> Coef for C {}

// ------------------------------------------------------------------ algebra

/// A polynomial in `h`, coefficients from `h^0` up.
type HPoly<C> = Vec<C>;

fn pmul<C: Coef>(a: &HPoly<C>, b: &HPoly<C>) -> HPoly<C> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![C::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j].add_assign_ref(&x.mul_ref(y));
        }
    }
    out
}

/// One term of a map: decoration (bit `i` set: circle `i` carries `X`),
/// power of `h`, coefficient.
#[derive(Clone, Debug)]
struct Term<C> {
    mask: u64,
    hpow: u32,
    c: C,
}

/// A map between two matchings, as a sum of terms (sorted, merged,
/// nonzero).  Most maps have one term, which is kept inline.
type Hom<C> = SmallVec<[Term<C>; 1]>;

fn normalize<C: Coef>(mut v: Vec<Term<C>>) -> Hom<C> {
    v.sort_by_key(|t| (t.mask, t.hpow));
    let mut out: Vec<Term<C>> = Vec::with_capacity(v.len());
    for t in v {
        match out.last_mut() {
            Some(l) if l.mask == t.mask && l.hpow == t.hpow => l.c.add_assign_ref(&t.c),
            _ => out.push(t),
        }
    }
    out.into_iter().filter(|t| !t.c.is_zero()).collect()
}

fn scale_hom<C: Coef>(h: &Hom<C>, s: &C) -> Hom<C> {
    h.iter()
        .map(|t| Term {
            mask: t.mask,
            hpow: t.hpow,
            c: t.c.mul_ref(s),
        })
        .collect()
}

/// The shape of a surface made of disks glued along boundary intervals:
/// its connected components, each with its genus and the final boundary
/// circles on it.  Maps with the same shape differ only in their dots.
struct Topo {
    /// The component of each piece.
    comp: Vec<usize>,
    /// Per component: genus and boundary circles.
    comps: Vec<(u32, Vec<usize>)>,
    circles: usize,
}

fn find(p: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while p[r] != r {
        r = p[r];
    }
    let mut y = x;
    while p[y] != r {
        let n = p[y];
        p[y] = r;
        y = n;
    }
    r
}

impl Topo {
    /// `pieces` disks; `glue` the identified pairs of boundary intervals;
    /// `circles[i]` the piece on which final boundary circle `i` lies.
    fn new(pieces: usize, glue: &[(usize, usize)], circles: &[usize]) -> Self {
        assert!(circles.len() <= 64, "more than 64 circles in one map");
        let mut parent: Vec<usize> = (0..pieces).collect();
        for &(a, b) in glue {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            parent[ra] = rb;
        }
        let mut index = vec![usize::MAX; pieces];
        let mut comp = vec![0; pieces];
        let mut stats: Vec<(i64, i64, Vec<usize>)> = Vec::new(); // pieces, gluings, circles
        for p in 0..pieces {
            let r = find(&mut parent, p);
            if index[r] == usize::MAX {
                index[r] = stats.len();
                stats.push((0, 0, Vec::new()));
            }
            comp[p] = index[r];
            stats[index[r]].0 += 1;
        }
        for &(a, _) in glue {
            stats[comp[a]].1 += 1;
        }
        for (i, &p) in circles.iter().enumerate() {
            stats[comp[p]].2.push(i);
        }
        let comps = stats
            .into_iter()
            .map(|(pc, gl, cs)| {
                let twice_g = 2 - (pc - gl) - cs.len() as i64;
                assert!(
                    twice_g >= 0 && twice_g % 2 == 0,
                    "a surface of negative or odd genus"
                );
                ((twice_g / 2) as u32, cs)
            })
            .collect();
        Topo {
            comp,
            comps,
            circles: circles.len(),
        }
    }

    /// The map with these dots on the pieces, times `c·h^hpow`, as terms
    /// over the final circles (not normalized).
    fn evaluate<C: Coef>(&self, dots: impl Fn(usize) -> u32, c: &C, hpow: u32, out: &mut Vec<Term<C>>) {
        let mut per = vec![0u32; self.comps.len()];
        for p in 0..self.comp.len() {
            per[self.comp[p]] += dots(p);
        }
        // the product over the components, as (decoration, polynomial in h)
        let mut acc: Vec<(u64, HPoly<C>)> = vec![(0, vec![c.clone()])];
        for (ci, (g, circles)) in self.comps.iter().enumerate() {
            let local = component_value::<C>(per[ci], *g, circles.len());
            if local.is_empty() {
                return;
            }
            let mut next = Vec::with_capacity(acc.len() * local.len());
            for (m, p) in &acc {
                for (lm, lp) in &local {
                    let mut mm = *m;
                    for (bit, &k) in circles.iter().enumerate() {
                        if lm >> bit & 1 == 1 {
                            mm |= 1 << k;
                        }
                    }
                    next.push((mm, pmul(p, lp)));
                }
            }
            acc = next;
        }
        for (m, p) in acc {
            for (j, x) in p.into_iter().enumerate() {
                if !x.is_zero() {
                    out.push(Term {
                        mask: m,
                        hpow: hpow + j as u32,
                        c: x,
                    });
                }
            }
        }
        let _ = self.circles;
    }
}

/// A connected surface with `k` boundary circles, genus `g` and `d` dots,
/// as an element of `A^{⊗k}`: `Δ^{(k)}(X^d (2X − h)^g)` (for `k = 0`, its
/// value under the counit).  Small cases are the common ones; they are
/// computed once and cached.
fn component_value<C: Coef>(d: u32, g: u32, k: usize) -> Vec<(u64, HPoly<C>)> {
    // e = X^d (2X − h)^g = a1 + aX·X, with integer coefficients
    let (mut a1, mut ax): (Vec<i64>, Vec<i64>) = if d == 0 {
        (vec![1], Vec::new())
    } else {
        let mut v = vec![0; d as usize];
        v[d as usize - 1] = 1;
        (Vec::new(), v)
    };
    for _ in 0..g {
        // (a1 + aX·X)(2X − h) = −h·a1 + X(2·a1 + h·aX)
        let mut n1 = vec![0];
        n1.extend(a1.iter().map(|x| -x));
        let mut nx: Vec<i64> = a1.iter().map(|x| 2 * x).collect();
        let mut hax = vec![0];
        hax.extend(ax.iter().copied());
        if nx.len() < hax.len() {
            nx.resize(hax.len(), 0);
        }
        for (x, y) in nx.iter_mut().zip(&hax) {
            *x += y;
        }
        a1 = n1;
        ax = nx;
    }
    let conv = |v: &[i64]| -> HPoly<C> { v.iter().map(|&x| C::from_i64(x)).collect() };
    let nonzero = |v: &[i64]| v.iter().any(|&x| x != 0);
    if k == 0 {
        // the counit picks out X
        return if nonzero(&ax) {
            vec![(0, conv(&ax))]
        } else {
            Vec::new()
        };
    }
    let mut out = Vec::new();
    if nonzero(&a1) {
        let a1c = conv(&a1);
        for (m, p) in delta_one(k) {
            let pc: HPoly<C> = p.iter().map(|&x| C::from_i64(x)).collect();
            out.push((m, pmul(&pc, &a1c)));
        }
    }
    if nonzero(&ax) {
        out.push(((1u64 << k) - 1, conv(&ax)));
    }
    out
}

/// Terms (decoration, integer polynomial in `h`) of an element of `A^{⊗k}`.
type IntTensor = Vec<(u64, Vec<i64>)>;

/// `Δ^{(k)}(1)` as (decoration of `k` circles, integer polynomial in `h`),
/// cached.
fn delta_one(k: usize) -> IntTensor {
    thread_local! {
        static CACHE: std::cell::RefCell<HashMap<usize, IntTensor>> = std::cell::RefCell::new(HashMap::new());
    }
    if let Some(v) = CACHE.with(|c| c.borrow().get(&k).cloned()) {
        return v;
    }
    let mut v: Vec<(u64, Vec<i64>)> = vec![(0, vec![1])];
    for j in 1..k {
        let mut next = Vec::with_capacity(v.len() * 3);
        for (m, p) in v {
            if m >> (j - 1) & 1 == 0 {
                // Δ(1) = 1⊗X + X⊗1 − h·1⊗1
                next.push((m | 1 << j, p.clone()));
                next.push((m | 1 << (j - 1), p.clone()));
                let mut hp = vec![0];
                hp.extend(p.iter().map(|x| -x));
                next.push((m, hp));
            } else {
                // Δ(X) = X⊗X
                next.push((m | 1 << j, p));
            }
        }
        v = next;
    }
    CACHE.with(|c| c.borrow_mut().insert(k, v.clone()));
    v
}

// ------------------------------------------------------------------ pictures

/// How many maps are built before they are added to a complex.
const CHUNK: usize = 1 << 16;

/// A crossingless matching of boundary labels: arcs `(a, b)`, `a < b`, sorted.
type Matching = Vec<(u32, u32)>;

fn arc_of(m: &Matching) -> HashMap<u32, usize> {
    let mut map = HashMap::with_capacity(2 * m.len());
    for (i, &(a, b)) in m.iter().enumerate() {
        map.insert(a, i);
        map.insert(b, i);
    }
    map
}

fn other(arc: (u32, u32), x: u32) -> u32 {
    if arc.0 == x {
        arc.1
    } else {
        arc.0
    }
}

/// The circles of `P ∪ Q` (two matchings of the same labels), each as its
/// arcs of `P` and of `Q`, in a canonical order.
fn circles(p: &Matching, q: &Matching) -> Vec<(Vec<usize>, Vec<usize>)> {
    let find_arc = |m: &Matching, x: u32| m.iter().position(|&(a, b)| a == x || b == x).unwrap();
    let mut seen = vec![false; p.len()];
    let mut out = Vec::new();
    for i in 0..p.len() {
        if seen[i] {
            continue;
        }
        let (mut ps, mut qs) = (Vec::new(), Vec::new());
        let start = p[i].0;
        let mut x = start;
        loop {
            let pi = find_arc(p, x);
            seen[pi] = true;
            ps.push(pi);
            let y = other(p[pi], x);
            let qi = find_arc(q, y);
            qs.push(qi);
            x = other(q[qi], y);
            if x == start {
                break;
            }
        }
        out.push((ps, qs));
    }
    out
}

/// An elementary arc of a glued picture: an arc of the old matching, or one
/// of the two arcs of the crossing's smoothing.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Elem {
    Old(usize),
    New(usize),
}

/// A crossing being added: its four labels counterclockwise from the
/// incoming under-strand, and what is shared with the tangle so far.
struct Step {
    labels: [u32; 4],
    /// Labels on both sides: the tangle's boundary and the crossing.
    shared: Vec<u32>,
    /// Pairs of slots of the crossing joined by an edge (a kink).
    self_pairs: Vec<(usize, usize)>,
}

/// The smoothings of `X[a, b, c, d]`: 0 joins slots `0–1`, `2–3`; 1 joins
/// `0–3`, `1–2`.
const SMOOTHING: [[(usize, usize); 2]; 2] = [[(0, 1), (2, 3)], [(0, 3), (1, 2)]];

/// An old matching glued to a smoothing of the crossing: the new matching
/// (canonical), the elementary arcs of each new arc and of each loop.
struct Glued {
    m: Matching,
    arcs: Vec<Vec<Elem>>,
    loops: Vec<Vec<Elem>>,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum End {
    Old(u32),
    Slot(usize),
}

fn glue(old: &Matching, sm: usize, step: &Step) -> Glued {
    let sigma = SMOOTHING[sm];
    // the two ends of each elementary arc
    let ends = |e: Elem| -> (End, End) {
        match e {
            Elem::Old(i) => (End::Old(old[i].0), End::Old(old[i].1)),
            Elem::New(j) => (End::Slot(sigma[j].0), End::Slot(sigma[j].1)),
        }
    };
    let mut at: HashMap<End, Elem> = HashMap::new();
    for i in 0..old.len() {
        let (a, b) = ends(Elem::Old(i));
        at.insert(a, Elem::Old(i));
        at.insert(b, Elem::Old(i));
    }
    for j in 0..2 {
        let (a, b) = ends(Elem::New(j));
        at.insert(a, Elem::New(j));
        at.insert(b, Elem::New(j));
    }
    // where an end continues to (shared labels, kinks); None: a new boundary point
    let mut link: HashMap<End, End> = HashMap::new();
    for &s in &step.shared {
        let slot = step.labels.iter().position(|&l| l == s).unwrap();
        link.insert(End::Old(s), End::Slot(slot));
        link.insert(End::Slot(slot), End::Old(s));
    }
    for &(i, j) in &step.self_pairs {
        link.insert(End::Slot(i), End::Slot(j));
        link.insert(End::Slot(j), End::Slot(i));
    }
    let label = |e: End| match e {
        End::Old(l) => l,
        End::Slot(s) => step.labels[s],
    };
    let mut used: HashSet<Elem> = HashSet::new();
    let mut arcs: Vec<((u32, u32), Vec<Elem>)> = Vec::new();
    // walk from every free end
    let free: Vec<End> = at.keys().copied().filter(|e| !link.contains_key(e)).collect();
    let mut free_sorted = free;
    free_sorted.sort_by_key(|&e| label(e));
    for start in free_sorted {
        let e0 = at[&start];
        if used.contains(&e0) {
            continue;
        }
        let mut chain = Vec::new();
        let mut end = start;
        loop {
            let e = at[&end];
            used.insert(e);
            chain.push(e);
            let (a, b) = ends(e);
            let far = if a == end { b } else { a };
            match link.get(&far) {
                None => {
                    let (x, y) = (label(start), label(far));
                    arcs.push(((x.min(y), x.max(y)), chain));
                    break;
                }
                Some(&next) => end = next,
            }
        }
    }
    // what is left closes up into loops
    let mut loops = Vec::new();
    let all: Vec<Elem> = (0..old.len())
        .map(Elem::Old)
        .chain((0..2).map(Elem::New))
        .collect();
    for &e0 in &all {
        if used.contains(&e0) {
            continue;
        }
        let mut chain = Vec::new();
        let (start, _) = ends(e0);
        let mut end = start;
        loop {
            let e = at[&end];
            if used.contains(&e) {
                break;
            }
            used.insert(e);
            chain.push(e);
            let (a, b) = ends(e);
            let far = if a == end { b } else { a };
            end = link[&far];
        }
        loops.push(chain);
    }
    arcs.sort_by_key(|a| a.0);
    Glued {
        m: arcs.iter().map(|a| a.0).collect(),
        arcs: arcs.into_iter().map(|a| a.1).collect(),
        loops,
    }
}

// ------------------------------------------------------------------ complexes

/// The surface of a composition `X → M → Y` of maps between matchings.
fn comp_surface(x: &Matching, mid: &Matching, y: &Matching) -> (usize, Topo) {
    let c1 = circles(x, mid);
    let c2 = circles(mid, y);
    let n1 = c1.len();
    // the circle of X ∪ M and of M ∪ Y through each arc of M
    let mut via1 = vec![0; mid.len()];
    for (i, (_, ms)) in c1.iter().enumerate() {
        for &a in ms {
            via1[a] = i;
        }
    }
    let mut via2 = vec![0; mid.len()];
    for (i, (ms, _)) in c2.iter().enumerate() {
        for &a in ms {
            via2[a] = n1 + i;
        }
    }
    let mut x_circle = vec![0; x.len()];
    for (i, (xs, _)) in c1.iter().enumerate() {
        for &a in xs {
            x_circle[a] = i;
        }
    }
    let fin = circles(x, y);
    let glue: Vec<(usize, usize)> = (0..mid.len()).map(|a| (via1[a], via2[a])).collect();
    let circ: Vec<usize> = fin.iter().map(|(xs, _)| x_circle[xs[0]]).collect();
    (n1, Topo::new(n1 + c2.len(), &glue, &circ))
}

struct Obj {
    /// The matching, as an index into [`Complex::mats`].
    m: u32,
    deg: i32,
    q: i32,
    alive: bool,
}

/// The surface of a composition `X → M → Y` (pieces: the circles of
/// `X ∪ M`, then those of `M ∪ Y`), and the number of the first.
type CompTopo = Arc<(usize, Topo)>;

/// A complex of matchings: objects, and the differential's entries.  The
/// matchings are interned, and the surface of each composition of three of
/// them is computed once.
struct Complex<C> {
    objs: Vec<Obj>,
    out: Vec<HashMap<usize, Hom<C>>>,
    inc: Vec<HashSet<usize>>,
    mats: Vec<Matching>,
    mat_id: HashMap<Matching, u32>,
    topos: HashMap<(u32, u32, u32), CompTopo>,
}

impl<C: Coef> Complex<C> {
    fn new() -> Self {
        Complex {
            objs: Vec::new(),
            out: Vec::new(),
            inc: Vec::new(),
            mats: Vec::new(),
            mat_id: HashMap::new(),
            topos: HashMap::new(),
        }
    }

    fn mat(&self, x: usize) -> &Matching {
        &self.mats[self.objs[x].m as usize]
    }

    fn add_obj(&mut self, m: Matching, deg: i32, q: i32) -> usize {
        let m = match self.mat_id.get(&m) {
            Some(&id) => id,
            None => {
                let id = self.mats.len() as u32;
                self.mats.push(m.clone());
                self.mat_id.insert(m, id);
                id
            }
        };
        self.objs.push(Obj {
            m,
            deg,
            q,
            alive: true,
        });
        self.out.push(HashMap::new());
        self.inc.push(HashSet::new());
        self.objs.len() - 1
    }

    fn add_entry(&mut self, a: usize, b: usize, h: Hom<C>) {
        if h.is_empty() {
            return;
        }
        let e = self.out[a].entry(b).or_default();
        if e.is_empty() {
            *e = h;
            self.inc[b].insert(a);
            return;
        }
        let mut v: Vec<Term<C>> = std::mem::take(e).into_vec();
        v.extend(h);
        let v = normalize(v);
        if v.is_empty() {
            self.out[a].remove(&b);
            self.inc[b].remove(&a);
        } else {
            *self.out[a].get_mut(&b).unwrap() = v;
            self.inc[b].insert(a);
        }
    }

    fn remove(&mut self, x: usize) {
        self.objs[x].alive = false;
        for y in std::mem::take(&mut self.out[x]).into_keys() {
            self.inc[y].remove(&x);
        }
        for w in std::mem::take(&mut self.inc[x]) {
            self.out[w].remove(&x);
        }
    }

    /// The surface of `X → M → Y`, cached.
    fn comp_topo(&mut self, x: u32, mid: u32, y: u32) -> CompTopo {
        if let Some(t) = self.topos.get(&(x, mid, y)) {
            return t.clone();
        }
        let (xm, mm, ym) = (
            &self.mats[x as usize],
            &self.mats[mid as usize],
            &self.mats[y as usize],
        );
        let t = Arc::new(comp_surface(xm, mm, ym));
        self.topos.insert((x, mid, y), t.clone());
        t
    }

    /// `g ∘ f` for `f: X → M`, `g: M → Y`, given the surface of `X → M → Y`.
    fn compose(topo: &(usize, Topo), f: &Hom<C>, g: &Hom<C>) -> Hom<C> {
        let (n1, ref topo) = *topo;
        let mut out = Vec::new();
        for s in f {
            for t in g {
                let dots = |p: usize| -> u32 {
                    if p < n1 {
                        (s.mask >> p & 1) as u32
                    } else {
                        (t.mask >> (p - n1) & 1) as u32
                    }
                };
                topo.evaluate(dots, &s.c.mul_ref(&t.c), s.hpow + t.hpow, &mut out);
            }
        }
        normalize(out)
    }

    /// The Markowitz cost of cancelling `b1 → b2`: the number of new maps
    /// it can create, `(incoming of b2 − 1)·(outgoing of b1 − 1)`.
    fn cost(&self, b1: usize, b2: usize) -> usize {
        (self.inc[b2].len() - 1) * (self.out[b1].len() - 1)
    }

    /// Cancels every isomorphism in the differential, cheapest first.  The
    /// costs change as cancellations proceed, so they are kept in a heap and
    /// rechecked when they come up; isomorphisms that cancellations create
    /// join the heap.
    fn eliminate(&mut self) {
        use std::cmp::Reverse;
        use std::collections::BinaryHeap;
        let mut heap = BinaryHeap::new();
        for b1 in 0..self.objs.len() {
            if !self.objs[b1].alive {
                continue;
            }
            for (&b2, h) in &self.out[b1] {
                if self.is_iso(b1, b2, h).is_some() {
                    heap.push(Reverse((self.cost(b1, b2), b1, b2)));
                }
            }
        }
        while let Some(Reverse((cost, b1, b2))) = heap.pop() {
            if !self.objs[b1].alive || !self.objs[b2].alive {
                continue;
            }
            let Some(h) = self.out[b1].get(&b2) else { continue };
            let Some(c) = self.is_iso(b1, b2, h) else { continue };
            let now = self.cost(b1, b2);
            if now > cost && heap.peek().is_some_and(|Reverse((next, _, _))| *next < now) {
                heap.push(Reverse((now, b1, b2)));
                continue;
            }
            for (x, y) in self.cancel(b1, b2, c) {
                if let Some(h) = self.out[x].get(&y) {
                    if self.is_iso(x, y, h).is_some() {
                        heap.push(Reverse((self.cost(x, y), x, y)));
                    }
                }
            }
        }
    }

    /// The scalar `c` if the map `b1 → b2` is `c·id`.
    fn is_iso(&self, b1: usize, b2: usize, h: &Hom<C>) -> Option<C> {
        let (o1, o2) = (&self.objs[b1], &self.objs[b2]);
        if o1.q != o2.q || o1.m != o2.m {
            return None;
        }
        h.iter().find(|t| t.mask == 0 && t.hpow == 0).map(|t| t.c.clone())
    }

    /// Gaussian elimination of the isomorphism `c·id: b1 → b2`; returns the
    /// entries it changed.
    fn cancel(&mut self, b1: usize, b2: usize, c: C) -> Vec<(usize, usize)> {
        let inv = c.inverse().negated();
        let xs: Vec<usize> = self.inc[b2].iter().copied().filter(|&x| x != b1).collect();
        let ys: Vec<usize> = self.out[b1].keys().copied().filter(|&y| y != b2).collect();
        let mid = self.objs[b1].m;
        let mut touched = Vec::new();
        // a chunk of sources at a time, which bounds the memory the new maps
        // take before they are added (adding x → y changes neither x → b2
        // nor b1 → y, since x ≠ b1 and y ≠ b2)
        for xc in xs.chunks((CHUNK / ys.len().max(1)).max(1)) {
            let mut jobs = Vec::with_capacity(xc.len() * ys.len());
            for &x in xc {
                for &y in &ys {
                    let t = self.comp_topo(self.objs[x].m, mid, self.objs[y].m);
                    jobs.push((x, y, t));
                }
            }
            let this = &*self;
            let run = |&(x, y, ref t): &(usize, usize, CompTopo)| {
                let comp = Self::compose(t, &this.out[x][&b2], &this.out[b1][&y]);
                (!comp.is_empty()).then(|| (x, y, scale_hom(&comp, &inv)))
            };
            #[cfg(feature = "parallel")]
            let upds: Vec<(usize, usize, Hom<C>)> = if jobs.len() >= 64 {
                use rayon::prelude::*;
                jobs.par_iter().filter_map(run).collect()
            } else {
                jobs.iter().filter_map(run).collect()
            };
            #[cfg(not(feature = "parallel"))]
            let upds: Vec<(usize, usize, Hom<C>)> = jobs.iter().filter_map(run).collect();
            drop(jobs);
            for (x, y, upd) in upds {
                touched.push((x, y));
                self.add_entry(x, y, upd);
            }
        }
        self.remove(b1);
        self.remove(b2);
        touched
    }
}

/// The order to add crossings in: a beam search over orders that grow the
/// tangle along its boundary, scored by the total of `2^{open ends}` over
/// the steps.  The size of the complexes grows with the number of matchings,
/// which grows like that; on cabled diagrams this keeps the largest complex
/// 20 times smaller than adding the crossing that shares the most edges.
fn crossing_order(d: &Diagram) -> Vec<usize> {
    beam_order(d, 256)
}

/// A diagram of the same link to run the local algorithm on: of
/// [`Diagram::simplify`] and [`Diagram::simplify_global`] (10 rounds), the
/// one with the cheaper crossing order ([`order_cost`]).  Fewer crossings do
/// not always mean a faster computation.
pub fn simplify_for_local(d: &Diagram) -> Diagram {
    let (a, b) = (d.simplify(), d.simplify_global(10));
    if order_cost(&b).0 < order_cost(&a).0 {
        b
    } else {
        a
    }
}

/// The cost of the crossing order the local algorithm uses on `d`: the
/// total of `2^{open ends}` over the steps, and the most ends open at once.
/// Of two diagrams of a knot, the one with the smaller total is usually the
/// faster to compute with.
pub fn order_cost(d: &Diagram) -> (f64, usize) {
    let mut done = vec![false; d.crossings()];
    let (mut total, mut most) = (0.0, 0);
    for c in crossing_order(d) {
        done[c] = true;
        let k = open_ends(d, &done);
        total += (2f64).powi(k as i32);
        most = most.max(k);
    }
    (total, most)
}

/// Open ends after adding the crossings in `done` (edges seen once).
fn open_ends(d: &Diagram, done: &[bool]) -> usize {
    let mut seen: HashMap<usize, u32> = HashMap::new();
    for (c, &x) in done.iter().enumerate() {
        if x {
            for &e in &d.crossing(c) {
                *seen.entry(e).or_insert(0) += 1;
            }
        }
    }
    seen.values().filter(|&&k| k == 1).count()
}

/// Beam search, keeping the `width` best partial orders (one per set of
/// crossings added).
fn beam_order(d: &Diagram, width: usize) -> Vec<usize> {
    let n = d.crossings();
    // (cost, order, done)
    let mut beam: Vec<(f64, Vec<usize>, Vec<bool>)> = vec![(0.0, Vec::new(), vec![false; n])];
    for _ in 0..n {
        let mut next: Vec<(f64, Vec<usize>, Vec<bool>)> = Vec::new();
        for (cost, order, done) in &beam {
            // candidates: crossings touching an open edge (any, at the start)
            let mut open: HashSet<usize> = HashSet::new();
            let mut cnt: HashMap<usize, u32> = HashMap::new();
            for (c, &x) in done.iter().enumerate() {
                if x {
                    for &e in &d.crossing(c) {
                        *cnt.entry(e).or_insert(0) += 1;
                    }
                }
            }
            for (&e, &k) in &cnt {
                if k == 1 {
                    open.insert(e);
                }
            }
            let cands: Vec<usize> = (0..n)
                .filter(|&c| !done[c] && (open.is_empty() || d.crossing(c).iter().any(|e| open.contains(e))))
                .collect();
            for c in cands {
                let mut nd = done.clone();
                nd[c] = true;
                let k = open_ends(d, &nd);
                let mut no = order.clone();
                no.push(c);
                next.push((cost + (2f64).powi(k as i32), no, nd));
            }
        }
        // the cheapest order to each set of crossings, then the cheapest sets
        next.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut seen: HashSet<Vec<bool>> = HashSet::new();
        next.retain(|x| seen.insert(x.2.clone()));
        next.truncate(width);
        beam = next;
    }
    beam.swap_remove(0).1
}

/// The reduced complex of the whole diagram: objects are empty matchings;
/// returns them with their (homological, quantum) degrees after the global
/// shift, and the differential.
fn reduced_complex<C: Coef>(d: &Diagram) -> Complex<C> {
    let mut cx: Complex<C> = Complex::new();
    // free loops: start from their delooped copies
    let start = cx.add_obj(Vec::new(), 0, 0);
    let _ = start;
    let mut boundary: HashSet<u32> = HashSet::new();
    for c in crossing_order(d) {
        let x = d.crossing(c);
        let labels = x.map(|e| e as u32);
        let mut shared = Vec::new();
        let mut self_pairs = Vec::new();
        for i in 0..4 {
            for j in i + 1..4 {
                if labels[i] == labels[j] {
                    self_pairs.push((i, j));
                }
            }
        }
        for &l in &labels {
            if boundary.contains(&l) && !shared.contains(&l) {
                shared.push(l);
            }
        }
        let step = Step {
            labels,
            shared: shared.clone(),
            self_pairs: self_pairs.clone(),
        };
        // the new boundary
        for &l in &labels {
            let twice = labels.iter().filter(|&&m| m == l).count() == 2;
            if boundary.contains(&l) {
                boundary.remove(&l);
            } else if !twice {
                boundary.insert(l);
            }
        }
        let t0 = std::time::Instant::now();
        cx = add_crossing(&cx, &step);
        let t_add = t0.elapsed().as_secs_f64();
        let before = cx.objs.len();
        let built: usize = cx.out.iter().map(|r| r.len()).sum();
        cx.eliminate();
        if std::env::var_os("NCKNOT_TRACE").is_some() {
            let alive = cx.objs.iter().filter(|o| o.alive).count();
            let entries: usize = cx.out.iter().map(|r| r.len()).sum();
            let terms: usize = cx.out.iter().flat_map(|r| r.values()).map(|h| h.len()).sum();
            eprintln!(
                "crossing {c}: {} open ends, {before} objects ({built} maps) -> {alive} ({entries} maps, {terms} terms), {:.3}s (adding {t_add:.3}s)",
                boundary.len(),
                t0.elapsed().as_secs_f64()
            );
        }
    }
    // free loops: each doubles the complex (∅{+1} ⊕ ∅{−1})
    for _ in 0..d.free_loops() {
        let mut next: Complex<C> = Complex::new();
        let mut map = HashMap::new();
        for (i, o) in cx.objs.iter().enumerate() {
            if !o.alive {
                continue;
            }
            for s in [1, -1] {
                map.insert((i, s), next.add_obj(Vec::new(), o.deg, o.q + s));
            }
        }
        for (i, row) in cx.out.iter().enumerate() {
            for (&j, h) in row {
                for s in [1, -1] {
                    next.add_entry(map[&(i, s)], map[&(j, s)], h.clone());
                }
            }
        }
        cx = next;
    }
    // the global shift [−n₋]{n₊ − 2n₋}
    let (np, nn) = d.sign_counts();
    for o in &mut cx.objs {
        o.deg -= nn as i32;
        o.q += np as i32 - 2 * nn as i32;
    }
    cx
}

/// A map between two pre-objects, apart from its terms: the surface that
/// glues on the crossing's part (pieces: the first `np` carry the old map's
/// dots, the crossing's come after), and how the result splits over the
/// delooped copies (`main` circles of the new matchings, `ls` and `lt`
/// loops).
struct Plan {
    topo: Topo,
    np: usize,
    ls: usize,
    lt: usize,
    main: usize,
}

/// (source pre-object, target pre-object, plan, terms, sign).
type Job<'a, C> = ((usize, usize), (usize, usize), Arc<Plan>, &'a Hom<C>, C);

/// Adds a crossing to the complex: the tensor product with
/// `smoothing 0 → smoothing 1{1}`, then delooping.  The surfaces depend only
/// on the matchings, so each is computed once; the maps are then built in
/// parallel.
fn add_crossing<C: Coef>(cx: &Complex<C>, step: &Step) -> Complex<C> {
    let mut next: Complex<C> = Complex::new();
    // per (matching, smoothing): the glued picture; per pre-object
    // (old object, smoothing): its delooped copies
    let mut glued: HashMap<(u32, usize), Glued> = HashMap::new();
    let mut ids: HashMap<(usize, usize, u64), usize> = HashMap::new();
    for (i, o) in cx.objs.iter().enumerate() {
        if !o.alive {
            continue;
        }
        for sm in 0..2 {
            let g = glued
                .entry((o.m, sm))
                .or_insert_with(|| glue(cx.mat(i), sm, step));
            let l = g.loops.len();
            assert!(l < 63, "too many loops");
            for choice in 0..(1u64 << l) {
                let ones = choice.count_ones() as i32;
                let q = o.q + sm as i32 + (l as i32 - ones) - ones;
                ids.insert((i, sm, choice), next.add_obj(g.m.clone(), o.deg + sm as i32, q));
            }
        }
    }
    let slot_of = |l: u32| step.labels.iter().position(|&m| m == l).unwrap();
    let id_hom: Hom<C> = smallvec![Term {
        mask: 0,
        hpow: 0,
        c: C::one(),
    }];
    let mut plans: HashMap<(u32, u32, usize), Arc<Plan>> = HashMap::new();
    let run = |(src, tgt, plan, f, sign): &Job<C>| -> Vec<(usize, usize, Hom<C>)> {
        let mut terms = Vec::new();
        for t in f.iter() {
            let dots = |p: usize| if p < plan.np { (t.mask >> p & 1) as u32 } else { 0 };
            plan.topo.evaluate(dots, &t.c, t.hpow, &mut terms);
        }
        let h = normalize(terms);
        deloop(&h, plan.ls, plan.lt, plan.main, sign)
            .into_iter()
            .map(|(cs, ct, v)| (ids[&(src.0, src.1, cs)], ids[&(tgt.0, tgt.1, ct)], v))
            .collect()
    };
    // the maps of a chunk of old objects at a time, which bounds the memory
    // the maps take before they are added
    let alive: Vec<usize> = (0..cx.objs.len()).filter(|&i| cx.objs[i].alive).collect();
    let mut k = 0;
    while k < alive.len() {
        let mut jobs: Vec<Job<C>> = Vec::new();
        while k < alive.len() && jobs.len() < CHUNK {
            let i = alive[k];
            k += 1;
            let o = &cx.objs[i];
            // (a) the old differential, with the identity on the smoothing
            for (&j, f) in &cx.out[i] {
                for sm in 0..2 {
                    let key = (o.m, cx.objs[j].m, sm);
                    let plan = plans
                        .entry(key)
                        .or_insert_with(|| {
                            let (gs, gt) = (&glued[&(key.0, sm)], &glued[&(key.1, sm)]);
                            Arc::new(differential_plan(
                                cx.mat(i),
                                cx.mat(j),
                                gs,
                                gt,
                                sm,
                                step,
                                &slot_of,
                            ))
                        })
                        .clone();
                    jobs.push(((i, sm), (j, sm), plan, f, C::one()));
                }
            }
            // (b) the saddle, with the identity on the old object; sign (−1)^deg
            let plan = plans
                .entry((o.m, o.m, 2))
                .or_insert_with(|| {
                    let (gs, gt) = (&glued[&(o.m, 0)], &glued[&(o.m, 1)]);
                    Arc::new(saddle_plan(cx.mat(i), gs, gt, step))
                })
                .clone();
            let sign = if o.deg % 2 == 0 {
                C::one()
            } else {
                C::one().negated()
            };
            jobs.push(((i, 0), (i, 1), plan, &id_hom, sign));
        }
        #[cfg(feature = "parallel")]
        let maps: Vec<Vec<(usize, usize, Hom<C>)>> = {
            use rayon::prelude::*;
            jobs.par_iter().with_min_len(64).map(run).collect()
        };
        #[cfg(not(feature = "parallel"))]
        let maps: Vec<Vec<(usize, usize, Hom<C>)>> = jobs.iter().map(run).collect();
        for (a, b, h) in maps.into_iter().flatten() {
            next.add_entry(a, b, h);
        }
    }
    next
}

/// The plan of a map `P → Q` of the old differential with the identity on
/// smoothing `sm`.  Pieces: the circles of `P ∪ Q`, then the smoothing's two
/// strips.
fn differential_plan(
    p: &Matching,
    q: &Matching,
    gs: &Glued,
    gt: &Glued,
    sm: usize,
    step: &Step,
    slot_of: &impl Fn(u32) -> usize,
) -> Plan {
    let pc = circles(p, q);
    let np = pc.len();
    let mut p_circle = vec![0; p.len()];
    let mut q_circle = vec![0; q.len()];
    for (k, (ps, qs)) in pc.iter().enumerate() {
        for &a in ps {
            p_circle[a] = k;
        }
        for &a in qs {
            q_circle[a] = k;
        }
    }
    let pa = arc_of(p);
    let strip = |slot: usize| -> usize {
        np + SMOOTHING[sm]
            .iter()
            .position(|&(a, b)| a == slot || b == slot)
            .unwrap()
    };
    let mut gl = Vec::new();
    for &s in &step.shared {
        gl.push((p_circle[pa[&s]], strip(slot_of(s))));
    }
    for &(a, b) in &step.self_pairs {
        gl.push((strip(a), strip(b)));
    }
    let piece_src = |e: Elem| match e {
        Elem::Old(a) => p_circle[a],
        Elem::New(k) => np + k,
    };
    let piece_tgt = |e: Elem| match e {
        Elem::Old(a) => q_circle[a],
        Elem::New(k) => np + k,
    };
    let (circ, ls, lt) = final_circles(gs, gt, piece_src, piece_tgt);
    Plan {
        topo: Topo::new(np + 2, &gl, &circ),
        np,
        ls,
        lt,
        main: circ_count(gs, gt),
    }
}

/// The plan of the saddle from smoothing 0 to smoothing 1 on `P`.  Pieces:
/// the arcs of `P` (strips), then the saddle; no dots come from an old map.
fn saddle_plan(p: &Matching, gs: &Glued, gt: &Glued, step: &Step) -> Plan {
    let pa = arc_of(p);
    let saddle = p.len();
    let mut gl = Vec::new();
    for &s in &step.shared {
        gl.push((pa[&s], saddle));
    }
    for _ in &step.self_pairs {
        gl.push((saddle, saddle));
    }
    let piece = |e: Elem| match e {
        Elem::Old(a) => a,
        Elem::New(_) => saddle,
    };
    let (circ, ls, lt) = final_circles(gs, gt, piece, piece);
    Plan {
        topo: Topo::new(saddle + 1, &gl, &circ),
        np: 0,
        ls,
        lt,
        main: circ_count(gs, gt),
    }
}

fn circ_count(gs: &Glued, gt: &Glued) -> usize {
    circles(&gs.m, &gt.m).len()
}

/// The final boundary circles of a map between two glued pictures: the
/// circles of the new matchings, then the source's loops, then the
/// target's; each as the piece it lies on.  Also the numbers of loops.
fn final_circles(
    gs: &Glued,
    gt: &Glued,
    piece_src: impl Fn(Elem) -> usize,
    piece_tgt: impl Fn(Elem) -> usize,
) -> (Vec<usize>, usize, usize) {
    let mut circ: Vec<usize> = circles(&gs.m, &gt.m)
        .iter()
        .map(|(ps, _)| piece_src(gs.arcs[ps[0]][0]))
        .collect();
    circ.extend(gs.loops.iter().map(|l| piece_src(l[0])));
    circ.extend(gt.loops.iter().map(|l| piece_tgt(l[0])));
    (circ, gs.loops.len(), gt.loops.len())
}

/// Splits a map between pre-objects over their delooped copies: (source
/// copy, target copy, map).  A source loop decorated `u` meets the copy
/// decorated `w` in `ε(u·w)`; a target loop decorated `u` lands in the copy
/// `u`.
fn deloop<C: Coef>(h: &Hom<C>, ls: usize, lt: usize, main: usize, sign: &C) -> Vec<(u64, u64, Hom<C>)> {
    let low = if main == 64 { u64::MAX } else { (1u64 << main) - 1 };
    let mut acc: HashMap<(u64, u64), Vec<Term<C>>> = HashMap::new();
    for t in h {
        let base = t.mask & low;
        let srcbits = (t.mask >> main) & ((1u64 << ls) - 1);
        let tgtbits = (t.mask >> (main + ls)) & ((1u64 << lt) - 1);
        // source loops: u = 1 → copy X only; u = X → copy 1, and copy X times h
        let mut opts: Vec<(u64, u32)> = vec![(0, 0)];
        for k in 0..ls {
            let u = srcbits >> k & 1;
            let mut nx = Vec::with_capacity(opts.len() * 2);
            for (c, e) in opts {
                if u == 0 {
                    nx.push((c | 1 << k, e));
                } else {
                    nx.push((c, e));
                    nx.push((c | 1 << k, e + 1));
                }
            }
            opts = nx;
        }
        for (c, e) in opts {
            acc.entry((c, tgtbits)).or_default().push(Term {
                mask: base,
                hpow: t.hpow + e,
                c: t.c.mul_ref(sign),
            });
        }
    }
    acc.into_iter()
        .map(|((cs, ct), v)| (cs, ct, normalize(v)))
        .collect()
}

// ------------------------------------------------------------------ results

fn kh_ranks<C: Coef>(d: &Diagram) -> BTreeMap<(i64, i64), usize> {
    let cx: Complex<C> = reduced_complex(d);
    let mut ranks = BTreeMap::new();
    for o in cx.objs.iter().filter(|o| o.alive) {
        *ranks.entry((o.deg as i64, o.q as i64)).or_insert(0) += 1;
    }
    ranks
}

/// Khovanov homology over `𝔽_p` by the local algorithm (`p` one of 2, 3, 5,
/// 7, 11, 13 or 2³¹ − 1).  The ranks agree with
/// [`crate::khovanov::khovanov`]'s; diagrams with many more crossings are
/// within reach.
pub fn khovanov(d: &Diagram, ring: Ring) -> Result<Khovanov> {
    let ranks = match ring {
        Ring::Fp(2) => kh_ranks::<Fp<2>>(d),
        Ring::Fp(3) => kh_ranks::<Fp<3>>(d),
        Ring::Fp(5) => kh_ranks::<Fp<5>>(d),
        Ring::Fp(7) => kh_ranks::<Fp<7>>(d),
        Ring::Fp(11) => kh_ranks::<Fp<11>>(d),
        Ring::Fp(13) => kh_ranks::<Fp<13>>(d),
        Ring::Fp(2147483647) => kh_ranks::<Fp<2147483647>>(d),
        _ => {
            return Err(KnotError::Unsupported(format!(
                "the local algorithm works over 𝔽_p for p in 2, 3, 5, 7, 11, 13, 2³¹ − 1, not {ring:?}"
            )))
        }
    };
    Ok(Khovanov {
        ring,
        ranks,
        torsion: BTreeMap::new(),
    })
}

/// The Rasmussen invariant over `F` from the reduced complex.
fn s_invariant<C: Coef>(d: &Diagram) -> Result<i64> {
    if d.components() != 1 {
        return Err(KnotError::Unsupported("the s-invariant is for knots".into()));
    }
    let cx: Complex<C> = reduced_complex(d);
    let alive: Vec<usize> = (0..cx.objs.len()).filter(|&i| cx.objs[i].alive).collect();
    let in_deg = |k: i32| -> Vec<usize> { alive.iter().copied().filter(|&i| cx.objs[i].deg == k).collect() };
    let (cm, c0, c1) = (in_deg(-1), in_deg(0), in_deg(1));
    let q = |i: usize| cx.objs[i].q as i64;
    // an entry c·h^k of the map `from → to`, as the scalar c (k is implied by the degrees)
    let entry = |a: usize, b: usize| -> C {
        match cx.out[a].get(&b) {
            None => C::zero(),
            Some(h) => {
                assert!(h.len() == 1 && h[0].mask == 0, "a non-monomial entry");
                assert_eq!(2 * h[0].hpow as i64, q(b) - q(a), "an entry of the wrong degree");
                h[0].c.clone()
            }
        }
    };
    // d0: rows c1, columns c0;  dm: rows c0, columns cm
    let mut d0: Vec<Vec<C>> = c1
        .iter()
        .map(|&b| c0.iter().map(|&a| entry(a, b)).collect())
        .collect();
    let mut dm: Vec<Vec<C>> = c0
        .iter()
        .map(|&b| cm.iter().map(|&a| entry(a, b)).collect())
        .collect();
    let (r1, r0): (Vec<i64>, Vec<i64>) = (
        c1.iter().map(|&i| q(i)).collect(),
        c0.iter().map(|&i| q(i)).collect(),
    );
    // 1. the kernel of d0, by elimination with the least power of h first;
    //    column operations on d0 change the basis of C^0, so dm's rows follow
    let mut col_done = vec![false; c0.len()];
    let mut row_done = vec![false; c1.len()];
    loop {
        let mut best: Option<(i64, usize, usize)> = None;
        for (j, row) in d0.iter().enumerate() {
            if row_done[j] {
                continue;
            }
            for (i, x) in row.iter().enumerate() {
                if !col_done[i] && !x.is_zero() {
                    let k = r1[j] - r0[i];
                    if best.is_none_or(|(b, _, _)| k < b) {
                        best = Some((k, j, i));
                    }
                }
            }
        }
        let Some((_, j, i)) = best else { break };
        let piv_inv = d0[j][i].inverse();
        // clear row j: col_l −= β col_i, and in dm row_i += β row_l
        for l in 0..c0.len() {
            if l == i || d0[j][l].is_zero() {
                continue;
            }
            let beta = d0[j][l].mul_ref(&piv_inv);
            for t in 0..c1.len() {
                let v = d0[t][i].mul_ref(&beta);
                d0[t][l].sub_assign_ref(&v);
            }
            for s in 0..cm.len() {
                let v = dm[l][s].mul_ref(&beta);
                dm[i][s].add_assign_ref(&v);
            }
        }
        // clear column i: row operations on C^1 only
        for t in 0..c1.len() {
            if t == j || d0[t][i].is_zero() {
                continue;
            }
            let beta = d0[t][i].mul_ref(&piv_inv);
            for l in 0..c0.len() {
                let v = d0[j][l].mul_ref(&beta);
                d0[t][l].sub_assign_ref(&v);
            }
        }
        row_done[j] = true;
        col_done[i] = true;
    }
    // the image of d_{−1} lies in the kernel
    for (i, &done) in col_done.iter().enumerate() {
        if done && dm[i].iter().any(|x| !x.is_zero()) {
            return Err(KnotError::Invalid("internal: d0 ∘ d−1 ≠ 0".into()));
        }
    }
    // 2. the cokernel of d_{−1} into the kernel: rows without a pivot are free
    let z: Vec<usize> = (0..c0.len()).filter(|&i| !col_done[i]).collect();
    let rm: Vec<i64> = cm.iter().map(|&i| q(i)).collect();
    let mut zrow_done = vec![false; z.len()];
    let mut mcol_done = vec![false; cm.len()];
    loop {
        let mut best: Option<(i64, usize, usize)> = None;
        for (zi, &row) in z.iter().enumerate() {
            if zrow_done[zi] {
                continue;
            }
            for s in 0..cm.len() {
                if !mcol_done[s] && !dm[row][s].is_zero() {
                    let k = r0[row] - rm[s];
                    if best.is_none_or(|(b, _, _)| k < b) {
                        best = Some((k, zi, s));
                    }
                }
            }
        }
        let Some((_, zi, s)) = best else { break };
        let row = z[zi];
        let piv_inv = dm[row][s].inverse();
        for &other in &z {
            if other == row || dm[other][s].is_zero() {
                continue;
            }
            let beta = dm[other][s].mul_ref(&piv_inv);
            for t in 0..cm.len() {
                let v = dm[row][t].mul_ref(&beta);
                dm[other][t].sub_assign_ref(&v);
            }
        }
        for t in 0..cm.len() {
            if t == s || dm[row][t].is_zero() {
                continue;
            }
            let beta = dm[row][t].mul_ref(&piv_inv);
            for &other in &z {
                let v = dm[other][s].mul_ref(&beta);
                dm[other][t].sub_assign_ref(&v);
            }
        }
        zrow_done[zi] = true;
        mcol_done[s] = true;
    }
    let free: Vec<i64> = z
        .iter()
        .enumerate()
        .filter(|&(zi, _)| !zrow_done[zi])
        .map(|(_, &row)| r0[row])
        .collect();
    if free.len() != 2 || (free[0] - free[1]).abs() != 2 {
        return Err(KnotError::Invalid(format!(
            "internal: the free part of degree 0 has generators in degrees {free:?}"
        )));
    }
    Ok((free[0] + free[1]) / 2)
}

/// The Rasmussen invariant `s` of a knot over ℚ, from Bar-Natan homology by
/// the local algorithm.  `s` bounds the slice genus: `|s| ≤ 2 g₄`.
pub fn rasmussen_s(d: &Diagram) -> Result<i64> {
    s_invariant::<BigRational>(d)
}

/// The Rasmussen invariant over `𝔽_p` (`p` one of 2, 3, 5, 7, 2³¹ − 1); it
/// can differ from the rational one for some knots and `p = 2`.
pub fn rasmussen_s_mod(d: &Diagram, p: u64) -> Result<i64> {
    match p {
        2 => s_invariant::<Fp<2>>(d),
        3 => s_invariant::<Fp<3>>(d),
        5 => s_invariant::<Fp<5>>(d),
        7 => s_invariant::<Fp<7>>(d),
        2147483647 => s_invariant::<Fp<2147483647>>(d),
        _ => Err(KnotError::Unsupported(format!("p = {p}"))),
    }
}
