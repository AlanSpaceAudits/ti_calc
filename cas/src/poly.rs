//! Univariate polynomials over the rationals: arithmetic, gcd, factoring, roots.
use crate::expr::*;
use crate::simplify::{add_v, expand, factor_int, mul, mul_v, pow};
use num_bigint::BigInt;
use num_complex::Complex64 as C;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

#[derive(Clone, Debug, PartialEq)]
pub struct Poly(pub Vec<Rat>);

fn ri(n: i64) -> Rat {
    Rat::from_integer(BigInt::from(n))
}

impl Poly {
    pub fn new(mut v: Vec<Rat>) -> Poly {
        while v.last().is_some_and(|c| c.is_zero()) {
            v.pop();
        }
        Poly(v)
    }
    pub fn zero() -> Poly {
        Poly(vec![])
    }
    pub fn constant(r: Rat) -> Poly {
        Poly::new(vec![r])
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
    pub fn deg(&self) -> usize {
        self.0.len().saturating_sub(1)
    }
    pub fn lead(&self) -> Rat {
        self.0.last().cloned().unwrap_or_else(Rat::zero)
    }
    pub fn coef(&self, i: usize) -> Rat {
        self.0.get(i).cloned().unwrap_or_else(Rat::zero)
    }
    pub fn add(&self, o: &Poly) -> Poly {
        Poly::new((0..self.0.len().max(o.0.len())).map(|i| self.coef(i) + o.coef(i)).collect())
    }
    pub fn sub(&self, o: &Poly) -> Poly {
        Poly::new((0..self.0.len().max(o.0.len())).map(|i| self.coef(i) - o.coef(i)).collect())
    }
    pub fn scale(&self, k: &Rat) -> Poly {
        Poly::new(self.0.iter().map(|c| c * k).collect())
    }
    pub fn mul(&self, o: &Poly) -> Poly {
        if self.is_zero() || o.is_zero() {
            return Poly::zero();
        }
        let mut v = vec![Rat::zero(); self.0.len() + o.0.len() - 1];
        for (i, a) in self.0.iter().enumerate() {
            for (j, b) in o.0.iter().enumerate() {
                v[i + j] += a * b;
            }
        }
        Poly::new(v)
    }
    pub fn powi(&self, n: u32) -> Poly {
        (0..n).fold(Poly::constant(Rat::one()), |acc, _| acc.mul(self))
    }
    pub fn divrem(&self, d: &Poly) -> (Poly, Poly) {
        let mut r = self.0.clone();
        if d.is_zero() || r.len() < d.0.len() {
            return (Poly::zero(), self.clone());
        }
        let mut q = vec![Rat::zero(); r.len() - d.0.len() + 1];
        let dl = d.lead();
        for i in (0..q.len()).rev() {
            let c = &r[i + d.0.len() - 1] / &dl;
            if !c.is_zero() {
                for (j, dc) in d.0.iter().enumerate() {
                    r[i + j] -= &c * dc;
                }
            }
            q[i] = c;
        }
        (Poly::new(q), Poly::new(r))
    }
    pub fn monic(&self) -> Poly {
        if self.is_zero() { self.clone() } else { self.scale(&self.lead().recip()) }
    }
    pub fn gcd(&self, o: &Poly) -> Poly {
        let (mut a, mut b) = (self.clone(), o.clone());
        while !b.is_zero() {
            let r = a.divrem(&b).1;
            a = b;
            b = r;
        }
        a.monic()
    }
    pub fn deriv(&self) -> Poly {
        Poly::new(self.0.iter().enumerate().skip(1).map(|(i, c)| c * ri(i as i64)).collect())
    }
    pub fn eval(&self, x: &Rat) -> Rat {
        self.0.iter().rev().fold(Rat::zero(), |acc, c| acc * x + c)
    }
    pub fn evalc(&self, x: C) -> C {
        self.0.iter().rev().fold(C::new(0.0, 0.0), |acc, c| acc * x + rat_f64(c))
    }
    /// (content, primitive part with coprime integer coefficients and positive leading coefficient)
    pub fn primitive(&self) -> (Rat, Poly) {
        if self.is_zero() {
            return (Rat::zero(), self.clone());
        }
        let mut num = BigInt::zero();
        let mut den = BigInt::one();
        for c in &self.0 {
            num = num.gcd(c.numer());
            den = den.lcm(c.denom());
        }
        let mut content = Rat::new(num, den);
        if self.lead().is_negative() {
            content = -content;
        }
        (content.clone(), self.scale(&content.recip()))
    }
    pub fn to_expr(&self, x: &Expr) -> Expr {
        add_v(self.0.iter().enumerate().filter(|(_, c)| !c.is_zero()).map(|(i, c)| mul(Expr::Num(c.clone()), pow(x.clone(), int(i as i64)))).collect())
    }
    pub fn compose_linear_shift(&self, s: &Rat) -> Poly {
        // p(x + s)
        let xs = Poly::new(vec![s.clone(), Rat::one()]);
        self.0.iter().rev().fold(Poly::zero(), |acc, c| acc.mul(&xs).add(&Poly::constant(c.clone())))
    }
}

/// Coefficients of `e` as a polynomial in `x`; each coefficient is free of `x`.
pub fn coeffs_expr(e: &Expr, x: &str) -> Option<Vec<Expr>> {
    let e = expand(e);
    let terms: Vec<Expr> = if let Expr::Add(v) = &e { v.clone() } else { vec![e.clone()] };
    let mut out: Vec<Vec<Expr>> = Vec::new();
    for t in terms {
        let factors: Vec<Expr> = if let Expr::Mul(v) = &t { v.clone() } else { vec![t.clone()] };
        let mut deg = 0usize;
        let mut coef = Vec::new();
        for f in factors {
            match &f {
                Expr::Sym(s) if s == x => deg += 1,
                Expr::Pow(b, n) if matches!(&**b, Expr::Sym(s) if s == x) => {
                    let k = n.as_i64().filter(|k| (0..=2000).contains(k))?;
                    deg += k as usize;
                }
                f if !f.has_sym(x) => coef.push(f.clone()),
                _ => return None,
            }
        }
        if out.len() <= deg {
            out.resize(deg + 1, Vec::new());
        }
        out[deg].push(mul_v(coef));
    }
    Some(out.into_iter().map(add_v).collect())
}

pub fn from_expr(e: &Expr, x: &str) -> Option<Poly> {
    let c = coeffs_expr(e, x)?;
    Some(Poly::new(c.iter().map(|k| k.as_rat().cloned()).collect::<Option<Vec<Rat>>>()?))
}

pub fn divisors(n: &BigInt) -> Option<Vec<BigInt>> {
    let mut ds = vec![BigInt::one()];
    for (p, k) in factor_int(n) {
        let base = ds.clone();
        let mut pk = BigInt::one();
        for _ in 0..k {
            pk *= &p;
            ds.extend(base.iter().map(|d| d * &pk));
        }
        if ds.len() > 4000 {
            return None;
        }
    }
    Some(ds)
}

/// Yun's square-free decomposition of a primitive polynomial.
pub fn squarefree(p: &Poly) -> Vec<(Poly, u32)> {
    let mut out = Vec::new();
    if p.deg() == 0 {
        return out;
    }
    let dp = p.deriv();
    let g = p.gcd(&dp);
    let mut c = p.divrem(&g).0;
    let mut d = dp.divrem(&g).0.sub(&c.deriv());
    let mut i = 1;
    while c.deg() > 0 {
        let a = c.gcd(&d);
        if a.deg() > 0 {
            out.push((a.primitive().1, i));
        }
        c = c.divrem(&a).0;
        d = d.divrem(&a).0.sub(&c.deriv());
        i += 1;
        if i > 10_000 {
            break;
        }
    }
    out
}

fn rational_roots(p: &Poly) -> Vec<Rat> {
    let mut roots = Vec::new();
    let (Some(a), Some(b)) = (divisors(p.0[0].numer()), divisors(p.lead().numer())) else { return roots };
    if a.len() * b.len() > 40_000 {
        return roots;
    }
    for q in &b {
        for n in &a {
            if !n.gcd(q).is_one() {
                continue;
            }
            for s in [1, -1] {
                let r = Rat::new(n * s, q.clone());
                if p.eval(&r).is_zero() && !roots.contains(&r) {
                    roots.push(r);
                }
            }
        }
    }
    roots
}

/// Kronecker's method: find an integer factor of exact degree `d`.
fn kronecker(g: &Poly, d: usize) -> Option<Poly> {
    let pts: Vec<Rat> = (0..=d as i64).map(|i| ri(if i % 2 == 1 { (i + 1) / 2 } else { -i / 2 })).collect();
    let mut divs: Vec<Vec<BigInt>> = Vec::new();
    let mut combos = 1usize;
    for (i, x) in pts.iter().enumerate() {
        let v = g.eval(x);
        let mut ds = divisors(v.numer())?;
        if i > 0 {
            let negs: Vec<BigInt> = ds.iter().map(|d| -d).collect();
            ds.extend(negs);
        }
        combos = combos.checked_mul(ds.len())?;
        // ponytail: search capped; factors needing a bigger search are left unsplit
        if combos > 120_000 {
            return None;
        }
        divs.push(ds);
    }
    // Lagrange basis
    let basis: Vec<Poly> = (0..=d)
        .map(|i| {
            let mut b = Poly::constant(Rat::one());
            for j in 0..=d {
                if j != i {
                    b = b.mul(&Poly::new(vec![-pts[j].clone(), Rat::one()])).scale(&(&pts[i] - &pts[j]).recip());
                }
            }
            b
        })
        .collect();
    let mut idx = vec![0usize; d + 1];
    loop {
        let mut cand = Poly::zero();
        for i in 0..=d {
            cand = cand.add(&basis[i].scale(&Rat::from_integer(divs[i][idx[i]].clone())));
        }
        if cand.deg() == d && cand.0.iter().all(|c| c.is_integer()) {
            let (q, r) = g.divrem(&cand);
            if r.is_zero() && q.0.iter().all(|c| c.is_integer()) {
                return Some(cand.primitive().1);
            }
        }
        let mut k = 0;
        loop {
            idx[k] += 1;
            if idx[k] < divs[k].len() {
                break;
            }
            idx[k] = 0;
            k += 1;
            if k > d {
                return None;
            }
        }
    }
}

fn factor_sf(g: &Poly, out: &mut Vec<Poly>) {
    let mut g = g.primitive().1;
    for r in rational_roots(&g) {
        let lin = Poly::new(vec![-Rat::from_integer(r.numer().clone()), Rat::from_integer(r.denom().clone())]);
        g = g.divrem(&lin).0;
        out.push(lin);
    }
    split(&g.primitive().1, out);
}

fn split(g: &Poly, out: &mut Vec<Poly>) {
    if g.deg() == 0 {
        return;
    }
    if g.deg() >= 4 {
        for d in 2..=(g.deg() / 2).min(4) {
            if let Some(f) = kronecker(g, d) {
                let q = g.divrem(&f).0.primitive().1;
                split(&f, out);
                split(&q, out);
                return;
            }
        }
    }
    out.push(g.clone());
}

/// Factor over the rationals: content and primitive integer factors with multiplicity.
pub fn factor(p: &Poly) -> (Rat, Vec<(Poly, u32)>) {
    let (content, prim) = p.primitive();
    let mut out: Vec<(Poly, u32)> = Vec::new();
    if prim.is_zero() {
        return (content, out);
    }
    let k = prim.0.iter().take_while(|c| c.is_zero()).count();
    if k > 0 {
        out.push((Poly::new(vec![Rat::zero(), Rat::one()]), k as u32));
    }
    let rest = Poly::new(prim.0[k..].to_vec());
    for (part, mult) in squarefree(&rest) {
        let mut fs = Vec::new();
        factor_sf(&part, &mut fs);
        out.extend(fs.into_iter().map(|f| (f, mult)));
    }
    out.sort_by(|a, b| a.0.deg().cmp(&b.0.deg()).then_with(|| b.0.0.cmp(&a.0.0)));
    (content, out)
}

/// All complex roots (Durand-Kerner with Newton polish). Input should be square-free.
pub fn roots_numeric(p: &Poly) -> Vec<C> {
    let n = p.deg();
    if n == 0 {
        return vec![];
    }
    let m = p.monic();
    let dp = m.deriv();
    let bound = 1.0 + m.0[..n].iter().map(|c| rat_f64(c).abs()).fold(0.0, f64::max);
    let mut z: Vec<C> = (0..n).map(|k| C::from_polar(bound.min(1e6).max(0.5) * 0.9, 0.7 + 2.0 * std::f64::consts::PI * k as f64 / n as f64)).collect();
    for _ in 0..2000 {
        let mut delta: f64 = 0.0;
        for i in 0..n {
            let mut den = C::new(1.0, 0.0);
            for j in 0..n {
                if j != i {
                    den *= z[i] - z[j];
                }
            }
            let step = m.evalc(z[i]) / den;
            z[i] -= step;
            delta = delta.max(step.norm());
        }
        if delta < 1e-15 * bound {
            break;
        }
    }
    for r in z.iter_mut() {
        for _ in 0..4 {
            let d = dp.evalc(*r);
            if d.norm() > 0.0 {
                let s = m.evalc(*r) / d;
                if s.norm().is_finite() {
                    *r -= s;
                }
            }
        }
        if r.im.abs() < 1e-11 * (1.0 + r.re.abs()) {
            r.im = 0.0;
        }
        if r.re.abs() < 1e-11 * (1.0 + r.im.abs()) {
            r.re = 0.0;
        }
    }
    z.sort_by(|a, b| a.re.total_cmp(&b.re).then(a.im.total_cmp(&b.im)));
    z
}

/// Solve a square rational linear system.
pub fn solve_rat(mut a: Vec<Vec<Rat>>, mut b: Vec<Rat>) -> Option<Vec<Rat>> {
    let n = b.len();
    for col in 0..n {
        let p = (col..n).find(|&i| !a[i][col].is_zero())?;
        a.swap(p, col);
        b.swap(p, col);
        let piv = a[col][col].clone();
        for j in col..n {
            a[col][j] = &a[col][j] / &piv;
        }
        b[col] = &b[col] / &piv;
        for i in 0..n {
            if i != col && !a[i][col].is_zero() {
                let f = a[i][col].clone();
                for j in col..n {
                    let t = &f * &a[col][j];
                    a[i][j] -= t;
                }
                let t = &f * &b[col];
                b[i] -= t;
            }
        }
    }
    Some(b)
}

pub struct PartialFractions {
    pub quotient: Poly,
    /// (numerator, irreducible factor, power)
    pub terms: Vec<(Poly, Poly, u32)>,
}

pub fn partial_fractions(num: &Poly, den: &Poly) -> Option<PartialFractions> {
    let (quotient, rem) = num.divrem(den);
    let (content, factors) = factor(den);
    let rem = rem.scale(&content.recip());
    let n = factors.iter().map(|(f, m)| f.deg() * *m as usize).sum::<usize>();
    if n == 0 {
        return Some(PartialFractions { quotient, terms: vec![] });
    }
    let full = factors.iter().fold(Poly::constant(Rat::one()), |acc, (f, m)| acc.mul(&f.powi(*m)));
    // unknown numerator coefficients, one basis polynomial each
    let mut basis: Vec<Poly> = Vec::new();
    let mut slots: Vec<(usize, u32, usize)> = Vec::new();
    for (fi, (f, m)) in factors.iter().enumerate() {
        for j in 1..=*m {
            let cof = full.divrem(&f.powi(j)).0;
            for l in 0..f.deg() {
                let mut xl = vec![Rat::zero(); l + 1];
                xl[l] = Rat::one();
                basis.push(cof.mul(&Poly::new(xl)));
                slots.push((fi, j, l));
            }
        }
    }
    let a: Vec<Vec<Rat>> = (0..n).map(|row| basis.iter().map(|b| b.coef(row)).collect()).collect();
    let rhs: Vec<Rat> = (0..n).map(|row| rem.coef(row)).collect();
    let sol = solve_rat(a, rhs)?;
    let mut terms: Vec<(Poly, Poly, u32)> = Vec::new();
    for (fi, (f, m)) in factors.iter().enumerate() {
        for j in 1..=*m {
            let mut c = vec![Rat::zero(); f.deg()];
            for (s, v) in slots.iter().zip(&sol) {
                if s.0 == fi && s.1 == j {
                    c[s.2] = v.clone();
                }
            }
            let np = Poly::new(c);
            if !np.is_zero() {
                terms.push((np, f.clone(), j));
            }
        }
    }
    Some(PartialFractions { quotient, terms })
}

