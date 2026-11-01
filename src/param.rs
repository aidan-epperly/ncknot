//! Polynomials in integer parameters, for families of diagrams.
//!
//! A [`Poly`] is a polynomial with rational coefficients in the parameters
//! `x_0, x_1, …` of a family (named `n`, `m`, … by the family).  Widths of
//! bundles, numbers of twists, framings and linking numbers of a family are
//! all such polynomials; rational coefficients are needed for counts like
//! `n(n−1)/2`.

use crate::error::{KnotError, Result};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::collections::BTreeMap;
use std::ops::{Add, Mul, Neg, Sub};

/// A polynomial in the parameters; see the [module docs](self).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Poly {
    /// Exponent vector (no trailing zeros) → nonzero coefficient.
    terms: BTreeMap<Vec<u32>, BigRational>,
}

fn trim(mut e: Vec<u32>) -> Vec<u32> {
    while e.last() == Some(&0) {
        e.pop();
    }
    e
}

impl Poly {
    /// The zero polynomial.
    pub fn zero() -> Self {
        Poly::default()
    }

    /// The constant `c`.
    pub fn constant(c: i64) -> Self {
        Self::rational(BigRational::from_integer(c.into()))
    }

    /// The rational constant `c`.
    pub fn rational(c: BigRational) -> Self {
        let mut p = Poly::zero();
        if !c.is_zero() {
            p.terms.insert(Vec::new(), c);
        }
        p
    }

    /// The parameter `x_i`.
    pub fn var(i: usize) -> Self {
        let mut e = vec![0; i + 1];
        e[i] = 1;
        let mut p = Poly::zero();
        p.terms.insert(e, BigRational::one());
        p
    }

    /// Is this the zero polynomial?
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// The value, if the polynomial is constant.
    pub fn as_constant(&self) -> Option<BigRational> {
        match self.terms.len() {
            0 => Some(BigRational::zero()),
            1 => self.terms.get(&Vec::new()).cloned(),
            _ => None,
        }
    }

    /// The value, if the polynomial is an integer constant.
    pub fn as_integer(&self) -> Option<i64> {
        self.as_constant()
            .filter(|c| c.is_integer())
            .and_then(|c| c.to_integer().to_i64())
    }

    /// The value at `values` (one per parameter).
    pub fn eval(&self, values: &[i64]) -> BigRational {
        let mut s = BigRational::zero();
        for (e, c) in &self.terms {
            let mut t = c.clone();
            for (i, &k) in e.iter().enumerate() {
                let v = BigInt::from(*values.get(i).unwrap_or(&0));
                t *= BigRational::from_integer(num_traits::pow(v, k as usize));
            }
            s += t;
        }
        s
    }

    /// The value at `values`, which must be an integer that fits an `i64`.
    pub fn eval_int(&self, values: &[i64]) -> Result<i64> {
        let v = self.eval(values);
        if !v.is_integer() {
            return Err(KnotError::Invalid(format!("{v} is not an integer")));
        }
        v.to_integer()
            .to_i64()
            .ok_or_else(|| KnotError::Invalid(format!("{v} is too large")))
    }

    /// The parity (`Some(true)` odd) for every value of the parameters with
    /// the given parities, if it is determined: the coefficients must be
    /// integers, and every parameter in a term with an odd coefficient must
    /// have a known parity.
    pub fn parity(&self, parities: &[Option<bool>]) -> Option<bool> {
        let mut odd = false;
        for (e, c) in &self.terms {
            if !c.is_integer() {
                return None;
            }
            if c.to_integer().is_even() {
                continue;
            }
            let mut t = true;
            for (i, &k) in e.iter().enumerate() {
                if k > 0 {
                    t &= (*parities.get(i)?)?;
                }
            }
            odd ^= t;
        }
        Some(odd)
    }

    /// The polynomial written with these parameter names, e.g. `n^2 + 2*n*m - 1`.
    pub fn display(&self, names: &[String]) -> String {
        if self.is_zero() {
            return "0".into();
        }
        let deg = |e: &Vec<u32>| e.iter().sum::<u32>();
        let mut terms: Vec<_> = self.terms.iter().collect();
        terms.sort_by(|a, b| deg(b.0).cmp(&deg(a.0)).then(b.0.cmp(a.0)));
        let mut s = String::new();
        for (j, (e, c)) in terms.into_iter().enumerate() {
            let mono: Vec<String> = e
                .iter()
                .enumerate()
                .filter(|(_, &k)| k > 0)
                .map(|(i, &k)| {
                    let name = names.get(i).cloned().unwrap_or_else(|| format!("x{i}"));
                    if k == 1 {
                        name
                    } else {
                        format!("{name}^{k}")
                    }
                })
                .collect();
            let neg = c.is_negative();
            let a = c.abs();
            let coeff = if a.is_integer() {
                a.to_integer().to_string()
            } else {
                format!("{}/{}", a.numer(), a.denom())
            };
            let body = match (mono.is_empty(), a.is_one()) {
                (true, _) => coeff,
                (false, true) => mono.join("*"),
                (false, false) => format!("{coeff}*{}", mono.join("*")),
            };
            match (j, neg) {
                (0, false) => s += &body,
                (0, true) => s += &format!("-{body}"),
                (_, false) => s += &format!(" + {body}"),
                (_, true) => s += &format!(" - {body}"),
            }
        }
        s
    }

    /// Parses an expression in the named parameters: integers, names,
    /// `+ - * /` (division by a constant only), `^` with a natural exponent,
    /// parentheses, and implicit products like `2n` or `n(n+1)`.
    pub fn parse(text: &str, names: &[String]) -> Result<Poly> {
        let mut p = Parser {
            s: text.chars().filter(|c| !c.is_whitespace()).collect(),
            i: 0,
            names,
        };
        if p.s.is_empty() {
            return Err(KnotError::Invalid("empty expression".into()));
        }
        let v = p.sum()?;
        if p.i < p.s.len() {
            return Err(KnotError::Invalid(format!("unexpected '{}' in {text}", p.s[p.i])));
        }
        Ok(v)
    }

    fn add_term(&mut self, e: Vec<u32>, c: BigRational) {
        if c.is_zero() {
            return;
        }
        let e = trim(e);
        let entry = self.terms.entry(e.clone()).or_insert_with(BigRational::zero);
        *entry += c;
        if entry.is_zero() {
            self.terms.remove(&e);
        }
    }
}

struct Parser<'a> {
    s: Vec<char>,
    i: usize,
    names: &'a [String],
}

impl Parser<'_> {
    fn peek(&self) -> Option<char> {
        self.s.get(self.i).copied()
    }

    fn sum(&mut self) -> Result<Poly> {
        let mut v = self.product()?;
        while let Some(c @ ('+' | '-')) = self.peek() {
            self.i += 1;
            let t = self.product()?;
            v = if c == '+' { &v + &t } else { &v - &t };
        }
        Ok(v)
    }

    fn product(&mut self) -> Result<Poly> {
        let mut v = self.power()?;
        loop {
            match self.peek() {
                Some('*') => {
                    self.i += 1;
                    v = &v * &self.power()?;
                }
                Some('/') => {
                    self.i += 1;
                    let d = self.power()?;
                    let c = d
                        .as_constant()
                        .filter(|c| !c.is_zero())
                        .ok_or_else(|| KnotError::Invalid("division only by a nonzero number".into()))?;
                    v = &v * &Poly::rational(c.recip());
                }
                // implicit product: 2n, n(n+1), (n+1)(n+2)
                Some(c) if c == '(' || c.is_alphabetic() => v = &v * &self.power()?,
                _ => return Ok(v),
            }
        }
    }

    fn power(&mut self) -> Result<Poly> {
        if self.peek() == Some('-') {
            self.i += 1;
            return Ok(-self.power()?);
        }
        let base = self.atom()?;
        if self.peek() == Some('^') {
            self.i += 1;
            let start = self.i;
            while self.peek().is_some_and(|c| c.is_ascii_digit()) {
                self.i += 1;
            }
            let k: u32 = self.s[start..self.i]
                .iter()
                .collect::<String>()
                .parse()
                .map_err(|_| KnotError::Invalid("^ needs a natural number".into()))?;
            let mut v = Poly::constant(1);
            for _ in 0..k {
                v = &v * &base;
            }
            return Ok(v);
        }
        Ok(base)
    }

    fn atom(&mut self) -> Result<Poly> {
        match self.peek() {
            Some('(') => {
                self.i += 1;
                let v = self.sum()?;
                if self.peek() != Some(')') {
                    return Err(KnotError::Invalid("missing )".into()));
                }
                self.i += 1;
                Ok(v)
            }
            Some(c) if c.is_ascii_digit() => {
                let start = self.i;
                while self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    self.i += 1;
                }
                let t: String = self.s[start..self.i].iter().collect();
                Ok(Poly::rational(BigRational::from_integer(
                    t.parse::<BigInt>().unwrap(),
                )))
            }
            Some(c) if c.is_alphabetic() || c == '_' => {
                let start = self.i;
                while self.peek().is_some_and(|c| c.is_alphanumeric() || c == '_') {
                    self.i += 1;
                }
                let name: String = self.s[start..self.i].iter().collect();
                // a name is a parameter, or a product of one-letter parameters (nm)
                if let Some(k) = self.names.iter().position(|x| *x == name) {
                    return Ok(Poly::var(k));
                }
                let mut v = Poly::constant(1);
                for ch in name.chars() {
                    let k = self
                        .names
                        .iter()
                        .position(|x| *x == ch.to_string())
                        .ok_or_else(|| KnotError::Invalid(format!("unknown parameter {name}")))?;
                    v = &v * &Poly::var(k);
                }
                Ok(v)
            }
            Some(c) => Err(KnotError::Invalid(format!("unexpected '{c}'"))),
            None => Err(KnotError::Invalid("the expression ends too early".into())),
        }
    }
}

impl From<i64> for Poly {
    fn from(c: i64) -> Self {
        Poly::constant(c)
    }
}

impl Add for &Poly {
    type Output = Poly;
    fn add(self, o: &Poly) -> Poly {
        let mut p = self.clone();
        for (e, c) in &o.terms {
            p.add_term(e.clone(), c.clone());
        }
        p
    }
}

impl Sub for &Poly {
    type Output = Poly;
    fn sub(self, o: &Poly) -> Poly {
        self + &(-o)
    }
}

impl Neg for &Poly {
    type Output = Poly;
    fn neg(self) -> Poly {
        Poly {
            terms: self.terms.iter().map(|(e, c)| (e.clone(), -c)).collect(),
        }
    }
}

impl Mul for &Poly {
    type Output = Poly;
    fn mul(self, o: &Poly) -> Poly {
        let mut p = Poly::zero();
        for (e1, c1) in &self.terms {
            for (e2, c2) in &o.terms {
                let len = e1.len().max(e2.len());
                let e = (0..len)
                    .map(|i| [e1, e2].iter().map(|e| e.get(i).unwrap_or(&0)).sum())
                    .collect();
                p.add_term(e, c1 * c2);
            }
        }
        p
    }
}

macro_rules! by_value {
    ($tr:ident, $f:ident) => {
        impl $tr for Poly {
            type Output = Poly;
            fn $f(self, o: Poly) -> Poly {
                (&self).$f(&o)
            }
        }
        impl $tr<&Poly> for Poly {
            type Output = Poly;
            fn $f(self, o: &Poly) -> Poly {
                (&self).$f(o)
            }
        }
        impl $tr<i64> for Poly {
            type Output = Poly;
            fn $f(self, o: i64) -> Poly {
                (&self).$f(&Poly::constant(o))
            }
        }
        impl $tr<i64> for &Poly {
            type Output = Poly;
            fn $f(self, o: i64) -> Poly {
                self.$f(&Poly::constant(o))
            }
        }
    };
}
by_value!(Add, add);
by_value!(Sub, sub);
by_value!(Mul, mul);

impl Neg for Poly {
    type Output = Poly;
    fn neg(self) -> Poly {
        -&self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic_and_parity() {
        let n = Poly::var(0);
        let m = Poly::var(1);
        let names = vec!["n".to_string(), "m".to_string()];
        let p = &(&n + 1) * &n - &m * 2;
        assert_eq!(p.display(&names), "n^2 + n - 2*m");
        assert_eq!(p.eval_int(&[3, 5]).unwrap(), 2);
        // n(n+1) is always even
        assert_eq!((&(&n + 1) * &n).parity(&[None, None]), None);
        assert_eq!((&m * 2 + 1).parity(&[None, None]), Some(true));
        assert_eq!(m.parity(&[None, Some(true)]), Some(true));
        assert!((&p - &p).is_zero());
        let parse = |t: &str| Poly::parse(t, &names).unwrap().display(&names);
        assert_eq!(parse("n^2 + n - 2*m"), "n^2 + n - 2*m");
        assert_eq!(parse("2n(n+1)/2 - nm"), "n^2 - n*m + n");
        assert_eq!(parse("-(m - 1)"), "-m + 1");
        assert!(Poly::parse("k + 1", &names).is_err());
        assert!(Poly::parse("n/m", &names).is_err());
    }
}
