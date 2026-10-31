//! Small knots and links by name, as braid closures (Rolfsen names; braid
//! words from KnotInfo).  The chirality is that of the braid word: `3_1` is
//! the right-handed (positive) trefoil.  Each entry is checked against its
//! Alexander polynomial and determinant in the tests.

use crate::braid::Braid;
use crate::diagram::Diagram;

/// `(name, strands, braid word)`.
pub const BRAIDS: &[(&str, usize, &[i32])] = &[
    ("3_1", 2, &[1, 1, 1]),
    ("4_1", 3, &[1, -2, 1, -2]),
    ("5_1", 2, &[1, 1, 1, 1, 1]),
    ("5_2", 3, &[1, 1, 1, 2, -1, 2]),
    ("6_1", 4, &[1, 1, 2, -1, -3, 2, -3]),
    ("6_2", 3, &[1, 1, 1, -2, 1, -2]),
    ("6_3", 3, &[1, 1, -2, 1, -2, -2]),
    ("7_1", 2, &[1, 1, 1, 1, 1, 1, 1]),
    ("7_2", 4, &[1, 1, 1, 2, -1, 2, 3, -2, 3]),
    ("7_3", 3, &[1, 1, 1, 1, 1, 2, -1, 2]),
    ("7_4", 4, &[1, 1, 2, -1, 2, 2, 3, -2, 3]),
    ("7_5", 3, &[1, 1, 1, 1, 2, -1, 2, 2]),
    ("7_6", 4, &[1, 1, -2, 1, 3, -2, 3]),
    ("7_7", 4, &[1, -2, 1, -2, 3, -2, 3]),
    ("8_19", 3, &[1, 2, 1, 2, 1, 2, 1, 2]),
    ("hopf", 2, &[1, 1]),
    ("T(2,4)", 2, &[1, 1, 1, 1]),
    ("borromean", 3, &[1, -2, 1, -2, 1, -2]),
];

/// The braid of a named knot or link from [`BRAIDS`].
pub fn braid(name: &str) -> Option<Braid> {
    BRAIDS
        .iter()
        .find(|(n, _, _)| *n == name)
        .map(|&(_, s, w)| Braid::new(s, w).expect("table braids are valid"))
}

/// The diagram of a named knot or link: `"0_1"` (the unknot) or a name from
/// [`BRAIDS`].
pub fn diagram(name: &str) -> Option<Diagram> {
    if name == "0_1" {
        return Some(Diagram::unknot());
    }
    braid(name).map(|b| b.closure())
}
