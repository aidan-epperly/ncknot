//! The Goeritz matrix and the signature (Gordon–Litherland).
//!
//! Checkerboard-color the faces of the projection.  With the unshaded faces
//! `X_0, …, X_m`, every crossing joins two unshaded corners and has an index
//! `η(c) = ±1`: `−1` when its unshaded corners are the regions that the
//! A-smoothing merges (the corners between `b, c` and between `d, a` of
//! `X[a, b, c, d]`), `+1` otherwise.  The pre-Goeritz matrix has
//! `G_ij = −Σ η(c)` over crossings joining `X_i ≠ X_j` and zero row sums; the
//! Goeritz matrix deletes row and column 0.  Gordon–Litherland:
//! `σ(L) = sign(G) − μ`, where `μ = Σ η(c)` over the crossings at which the
//! oriented (Seifert) smoothing merges the shaded corners.  With this sign
//! convention positive knots have negative signature (the right-handed
//! trefoil has `σ = −2`).  The inertia of `G` is exact (ncalg
//! `linalg::exact::inertia` over ℚ).

use crate::diagram::Diagram;
use crate::error::{KnotError, Result};
use ncalg::linalg::exact;
use ncalg::BigRational;
use num_bigint::BigInt;

/// The Goeritz matrix of the checkerboard coloring in which the face
/// containing corner 0 of crossing 0 is unshaded, and the correction `μ`.
pub fn goeritz(d: &Diagram) -> Result<(Vec<Vec<i64>>, i64)> {
    let n = d.crossings();
    if n == 0 {
        return Ok((Vec::new(), 0));
    }
    if !d.is_connected_projection() {
        return Err(KnotError::Unsupported(
            "Goeritz matrix of a diagram whose projection is disconnected".into(),
        ));
    }
    let faces = d.faces();
    let mut face_of = vec![0usize; 4 * n];
    for (f, face) in faces.iter().enumerate() {
        for &(c, k) in face {
            face_of[4 * c + k as usize] = f;
        }
    }
    // 2-color the faces: corners k and k + 1 of a crossing differ
    let mut color: Vec<Option<bool>> = vec![None; faces.len()];
    let mut stack = vec![(face_of[0], false)];
    while let Some((f, col)) = stack.pop() {
        match color[f] {
            Some(c) if c == col => continue,
            Some(_) => return Err(KnotError::Invalid("faces are not 2-colorable".into())),
            None => color[f] = Some(col),
        }
        for &(c, k) in &faces[f] {
            for k2 in [(k + 1) % 4, (k + 3) % 4] {
                stack.push((face_of[4 * c + k2 as usize], !col));
            }
        }
    }
    let white: Vec<usize> = (0..faces.len()).filter(|&f| color[f] == Some(false)).collect();
    let index = |f: usize| white.iter().position(|&w| w == f).unwrap();
    let m = white.len();
    let mut pre = vec![vec![0i64; m]; m];
    let mut mu = 0;
    for c in 0..n {
        // the unshaded corners: {0, 2} or {1, 3}; A merges {1, 3}
        let white_is_a = color[face_of[4 * c + 1]] == Some(false);
        let eta = if white_is_a { -1 } else { 1 };
        let (f1, f2) = if white_is_a {
            (face_of[4 * c + 1], face_of[4 * c + 3])
        } else {
            (face_of[4 * c], face_of[4 * c + 2])
        };
        let (i, j) = (index(f1), index(f2));
        if i != j {
            pre[i][j] -= eta;
            pre[j][i] -= eta;
            pre[i][i] += eta;
            pre[j][j] += eta;
        }
        // the Seifert smoothing is the A-smoothing at positive crossings
        let seifert_is_a = d.sign(c) > 0;
        if seifert_is_a != white_is_a {
            mu += eta;
        }
    }
    let g = pre[1..].iter().map(|r| r[1..].to_vec()).collect();
    Ok((g, mu))
}

/// The signature of the link (see the [module docs](self) for the sign).
/// Crossingless loops contribute nothing; a disconnected projection is
/// declined.
pub fn signature(d: &Diagram) -> Result<i64> {
    let (g, mu) = goeritz(d)?;
    let q: Vec<Vec<BigRational>> = g
        .iter()
        .map(|r| {
            r.iter()
                .map(|&x| BigRational::from_integer(BigInt::from(x)))
                .collect()
        })
        .collect();
    let (pos, neg, _) = exact::inertia(&q);
    Ok(pos as i64 - neg as i64 - mu)
}
