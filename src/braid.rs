//! Braids and their closures.
//!
//! A braid on `n` strands is a word in the Artin generators `σ_1, …,
//! σ_{n-1}`, written as nonzero integers: `i` is `σ_i`, `-i` is `σ_i^{-1}`.
//! Strands run upward; in `σ_i` the strand at position `i` crosses **over**
//! the one at `i + 1`, which makes `σ_i` a positive crossing in the closure.

use crate::diagram::Diagram;
use crate::error::{KnotError, Result};

/// A braid word on `strands` strands.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Braid {
    strands: usize,
    word: Vec<i32>,
}

impl Braid {
    /// The braid `σ_{|w_1|}^{±1} σ_{|w_2|}^{±1} ⋯` on `strands ≥ 1` strands.
    pub fn new(strands: usize, word: &[i32]) -> Result<Self> {
        if strands == 0 {
            return Err(KnotError::Invalid("a braid needs at least one strand".into()));
        }
        if let Some(&g) = word
            .iter()
            .find(|&&g| g == 0 || g.unsigned_abs() as usize >= strands)
        {
            return Err(KnotError::Invalid(format!(
                "generator {g} is not σ_i^±1 with 1 ≤ i < {strands}"
            )));
        }
        Ok(Braid {
            strands,
            word: word.to_vec(),
        })
    }

    /// Number of strands.
    pub fn strands(&self) -> usize {
        self.strands
    }

    /// The word.
    pub fn word(&self) -> &[i32] {
        &self.word
    }

    /// The product `self · other` (on the same number of strands).
    pub fn compose(&self, other: &Braid) -> Result<Braid> {
        if self.strands != other.strands {
            return Err(KnotError::Invalid(
                "braids on different numbers of strands".into(),
            ));
        }
        let mut w = self.word.clone();
        w.extend_from_slice(&other.word);
        Braid::new(self.strands, &w)
    }

    /// The inverse braid.
    pub fn inverse(&self) -> Braid {
        Braid {
            strands: self.strands,
            word: self.word.iter().rev().map(|g| -g).collect(),
        }
    }

    /// The Markov stabilization `β σ_n^{±1}` on `n + 1` strands (same closure).
    pub fn stabilize(&self, positive: bool) -> Braid {
        let mut w = self.word.clone();
        let g = self.strands as i32;
        w.push(if positive { g } else { -g });
        Braid {
            strands: self.strands + 1,
            word: w,
        }
    }

    /// The torus link `T(p, q) = (σ_1 ⋯ σ_{p-1})^q` on `p` strands (`q < 0`
    /// gives the mirror image).
    pub fn torus(p: usize, q: i32) -> Result<Braid> {
        let s = q.signum();
        let w: Vec<i32> = (0..q.unsigned_abs())
            .flat_map(|_| (1..p as i32).map(move |i| s * i))
            .collect();
        Braid::new(p, &w)
    }

    /// The permutation of the strands: position `k` at the bottom ends at `perm[k]`.
    pub fn permutation(&self) -> Vec<usize> {
        let mut at: Vec<usize> = (0..self.strands).collect(); // at[pos] = strand
        for &g in &self.word {
            let i = g.unsigned_abs() as usize - 1;
            at.swap(i, i + 1);
        }
        let mut perm = vec![0; self.strands];
        for (pos, &k) in at.iter().enumerate() {
            perm[k] = pos;
        }
        perm
    }

    /// The closure as a PD diagram.  Strands that never cross become
    /// crossingless loops.
    pub fn closure(&self) -> Diagram {
        let n = self.strands;
        // cur[k]: the label of the open edge at position k; labels 0..n are the bottom edges
        let mut cur: Vec<i64> = (0..n as i64).collect();
        let mut next = n as i64;
        let mut code: Vec<[i64; 4]> = Vec::new();
        for &g in &self.word {
            let l = g.unsigned_abs() as usize - 1;
            let r = l + 1;
            let (nl, nr) = (next, next + 1);
            next += 2;
            // corners counterclockwise from the incoming under-strand
            code.push(if g > 0 {
                [cur[r], nr, nl, cur[l]] // under BR→TL, over BL→TR
            } else {
                [cur[l], cur[r], nr, nl] // under BL→TR, over BR→TL
            });
            cur[l] = nl;
            cur[r] = nr;
        }
        // close up: the top edge at position k is the bottom edge k
        let mut free = 0;
        let mut rename: Vec<(i64, i64)> = Vec::new();
        for (k, &top) in cur.iter().enumerate() {
            if top == k as i64 {
                free += 1;
            } else {
                rename.push((top, k as i64));
            }
        }
        for x in &mut code {
            for l in x.iter_mut() {
                if let Some(&(_, k)) = rename.iter().find(|&&(t, _)| t == *l) {
                    *l = k;
                }
            }
        }
        Diagram::build(&code, free).expect("braid closures are valid diagrams")
    }
}
