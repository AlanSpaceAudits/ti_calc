//! Rational-function algebra: common denominators, cancellation, simplify, factor, partial fractions.
use crate::expr::*;
use crate::poly::{self};
use crate::simplify::*;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed};

/// Split into numerator and denominator, combining sums over a common denominator.
pub fn together(e: &Expr) -> (Expr, Expr) {
    match e {
        Expr::Num(r) => (big(r.numer().clone()), big(r.denom().clone())),
        Expr::Add(v) => {
            let parts: Vec<(Expr, Expr)> = v.iter().map(together).collect();
            // common denominator: each distinct base to its highest integer power
            let mut bases: Vec<(Expr, i64)> = Vec::new();
            let mut num_lcm = BigInt::one();
            for (_, d) in &parts {
                for f in factors_of(d) {
                    if let Some(n) = f.as_int() {
                        num_lcm = num_lcm.lcm(n);
                        continue;
                    }
                    let (b, k) = base_exp(&f);
                    match bases.iter_mut().find(|(bb, _)| *bb == b) {
                        Some(slot) => slot.1 = slot.1.max(k),
                        None => bases.push((b, k)),
                    }
                }
            }
            let mut dfac: Vec<Expr> = vec![big(num_lcm)];
            dfac.extend(bases.into_iter().map(|(b, k)| pow(b, int(k))));
            let den = mul_v(dfac);
            let num = add_v(parts.into_iter().map(|(n, d)| mul_v(vec![n, den.clone(), pow(d, int(-1))])).collect());
            (num, den)
        }
        Expr::Mul(v) => {
            let (mut n, mut d) = (Vec::new(), Vec::new());
            for f in v {
                let (a, b) = together(f);
                n.push(a);
                d.push(b);
            }
            (mul_v(n), mul_v(d))
        }
        Expr::Pow(b, x) => {
            if let Some(k) = x.as_i64() {
                let (n, d) = together(b);
                return if k > 0 { (pow(n, int(k)), pow(d, int(k))) } else { (pow(d, int(-k)), pow(n, int(-k))) };
            }
            if x.has_neg_coeff() { (int(1), pow((**b).clone(), neg((**x).clone()))) } else { (e.clone(), int(1)) }
        }
        _ => (e.clone(), int(1)),
    }
}

fn factors_of(e: &Expr) -> Vec<Expr> {
    match e {
        Expr::Mul(v) => v.clone(),
        e => vec![e.clone()],
    }
}

fn base_exp(f: &Expr) -> (Expr, i64) {
    if let Expr::Pow(b, k) = f {
        if let Some(k) = k.as_i64() {
            return ((**b).clone(), k);
        }
    }
    (f.clone(), 1)
}

/// Cancel a common polynomial factor when numerator and denominator are univariate over Q.
pub fn cancel(n: &Expr, d: &Expr) -> (Expr, Expr) {
    let (n, d) = (expand(n), expand(d));
    if d.is_one() {
        return (n, d);
    }
    if let Some(r) = d.as_rat() {
        return (expand(&mul(n, Expr::Num(r.recip()))), int(1));
    }
    let mut syms = n.free_syms();
    syms.extend(d.free_syms());
    if syms.len() == 1 {
        let x = syms.iter().next().unwrap();
        if let (Some(pn), Some(pd)) = (poly::from_expr(&n, x), poly::from_expr(&d, x)) {
            if pn.is_zero() {
                return (int(0), int(1));
            }
            let g = pn.gcd(&pd);
            let (qn, qd) = (pn.divrem(&g).0, pd.divrem(&g).0);
            let (c, qd) = qd.primitive();
            let xs = sym(x);
            return (qn.scale(&c.recip()).to_expr(&xs), qd.to_expr(&xs));
        }
    }
    // multivariate: divide when the denominator is a single term or equals the numerator's content
    if n == d {
        return (int(1), int(1));
    }
    if let Some(q) = div_exact(&n, &d) {
        return (q, int(1));
    }
    if !matches!(d, Expr::Add(_)) {
        let r = expand(&div(n.clone(), d.clone()));
        if !r.has(&|x| matches!(x, Expr::Pow(_, k) if k.has_neg_coeff())) {
            return (r, int(1));
        }
    }
    (n, d)
}

/// Polynomial long division in one variable with symbolic coefficients; Some(q) when n = q*d exactly.
fn div_exact(n: &Expr, d: &Expr) -> Option<Expr> {
    let neg_pow = |e: &Expr| e.has(&|x| matches!(x, Expr::Pow(_, k) if k.has_neg_coeff()));
    for var in d.free_syms() {
        let (Some(mut rem), Some(dc)) = (poly::coeffs_expr(n, &var), poly::coeffs_expr(d, &var)) else { continue };
        let lead = dc.last()?.clone();
        if dc.len() < 2 || matches!(lead, Expr::Add(_)) || rem.len() > 60 {
            continue;
        }
        let mut quo = vec![int(0); rem.len().saturating_sub(dc.len()) + 1];
        let mut ok = true;
        while rem.len() >= dc.len() {
            let k = rem.len() - dc.len();
            let q = expand(&div(rem.last().unwrap().clone(), lead.clone()));
            if neg_pow(&q) {
                ok = false;
                break;
            }
            for (i, c) in dc.iter().enumerate() {
                rem[k + i] = expand(&sub(rem[k + i].clone(), mul(q.clone(), c.clone())));
            }
            quo[k] = q;
            rem.pop();
        }
        if ok && rem.iter().all(|c| c.is_zero()) {
            let x = sym(&var);
            return Some(expand(&add_v(quo.into_iter().enumerate().map(|(i, c)| mul(c, pow(x.clone(), int(i as i64)))).collect())));
        }
    }
    None
}

pub fn ratsimp(e: &Expr) -> Expr {
    match e {
        Expr::Rel(r, a, b) => return Expr::Rel(*r, Box::new(ratsimp(a)), Box::new(ratsimp(b))),
        Expr::List(_) | Expr::Matrix(_) => return e.map(&mut |c| ratsimp(c)),
        _ => {}
    }
    let (n, d) = together(e);
    let (n, d) = cancel(&n, &d);
    // a factored denominator is usually the shorter form
    let fd = if matches!(d, Expr::Add(_)) { factor_poly_expr(&d) } else { d.clone() };
    div(n, if fd.size() < d.size() { fd } else { d })
}

/// Clear a two-term radical denominator by multiplying with its conjugate.
pub fn rationalize(e: &Expr) -> Expr {
    let (n, d) = together(e);
    let Expr::Add(terms) = &d else { return e.clone() };
    let is_rad = |t: &Expr| t.has(&|x| matches!(x, Expr::Pow(_, k) if k.as_rat().is_some_and(|r| *r.denom() == BigInt::from(2))));
    if terms.len() != 2 || !d.is_constant() || terms.iter().filter(|t| is_rad(t)).count() == 0 {
        return e.clone();
    }
    let flip = if is_rad(&terms[1]) { 1 } else { 0 };
    let conj = add(terms[1 - flip].clone(), neg(terms[flip].clone()));
    let nd = expand(&mul(d.clone(), conj.clone()));
    if nd.as_rat().is_none() || nd.is_zero() {
        return e.clone();
    }
    expand(&div(mul(n, conj), nd))
}

fn trig_pythag(e: &Expr, from: &str, to: &str) -> Expr {
    // from(u)^2 -> 1 - to(u)^2, including higher even powers
    match e {
        Expr::Pow(b, k) => {
            if let (Expr::Func(n, a), Some(k)) = (&**b, k.as_i64()) {
                if n == from && k >= 2 {
                    let base = sub(int(1), pow(f1(to, a[0].clone()), int(2)));
                    return mul(pow(base, int(k / 2)), pow((**b).clone(), int(k % 2)));
                }
            }
            pow(trig_pythag(b, from, to), trig_pythag(k, from, to))
        }
        Expr::Add(v) => add_v(v.iter().map(|c| trig_pythag(c, from, to)).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(|c| trig_pythag(c, from, to)).collect()),
        Expr::Func(n, v) => func(n, v.iter().map(|c| trig_pythag(c, from, to)).collect()),
        e => e.clone(),
    }
}

fn tan_to_sincos(e: &Expr) -> Expr {
    match e {
        Expr::Func(n, a) if n == "tan" => {
            let u = tan_to_sincos(&a[0]);
            div(f1("sin", u.clone()), f1("cos", u))
        }
        Expr::Add(v) => add_v(v.iter().map(tan_to_sincos).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(tan_to_sincos).collect()),
        Expr::Pow(b, x) => pow(tan_to_sincos(b), tan_to_sincos(x)),
        Expr::Func(n, v) => func(n, v.iter().map(tan_to_sincos).collect()),
        e => e.clone(),
    }
}

fn score(e: &Expr) -> usize {
    // prefer fewer nodes; penalize unexpanded nesting slightly less than raw size
    e.size()
}

pub fn simplify(e: &Expr) -> Expr {
    match e {
        Expr::Rel(r, a, b) => return Expr::Rel(*r, Box::new(simplify(a)), Box::new(simplify(b))),
        Expr::List(_) | Expr::Matrix(_) => return e.map(&mut |c| simplify(c)),
        _ => {}
    }
    // simplify inside function arguments first
    let e = match e {
        Expr::Func(n, v) => func(n, v.iter().map(simplify).collect()),
        Expr::Add(v) => add_v(v.iter().map(simplify_args).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(simplify_args).collect()),
        Expr::Pow(b, x) => pow(simplify_args(b), simplify_args(x)),
        e => e.clone(),
    };
    let mut cands = vec![e.clone(), expand(&e), ratsimp(&e), rationalize(&e)];
    let has_trig = e.has(&|x| matches!(x, Expr::Func(n, _) if n == "sin" || n == "cos" || n == "tan"));
    if has_trig {
        let t = tan_to_sincos(&e);
        for base in [e.clone(), t] {
            for (a, b) in [("sin", "cos"), ("cos", "sin")] {
                let r = trig_pythag(&expand(&base), a, b);
                cands.push(expand(&r));
                cands.push(ratsimp(&r));
            }
            cands.push(ratsimp(&base));
        }
    }
    if e.has(&|x| matches!(x, Expr::Func(n, _) if n == "ln")) {
        cands.push(combine_logs(&expand(&e)));
    }
    cands.into_iter().filter(|c| !c.is_undef() || e.is_undef()).min_by_key(score).unwrap_or(e)
}

fn simplify_args(e: &Expr) -> Expr {
    match e {
        Expr::Func(n, v) => func(n, v.iter().map(simplify).collect()),
        Expr::Add(v) => add_v(v.iter().map(simplify_args).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(simplify_args).collect()),
        Expr::Pow(b, x) => pow(simplify_args(b), simplify_args(x)),
        e => e.clone(),
    }
}

/// a*ln(x) + b*ln(y) -> ln(x^a * y^b) for rational a, b.
pub fn combine_logs(e: &Expr) -> Expr {
    let Expr::Add(terms) = e else { return e.clone() };
    let mut inside = Vec::new();
    let mut rest = Vec::new();
    for t in terms {
        let (c, r) = match t {
            Expr::Mul(v) if v.len() == 2 && v[0].as_rat().is_some() => (v[0].clone(), v[1].clone()),
            t => (int(1), t.clone()),
        };
        match &r {
            Expr::Func(n, a) if n == "ln" => inside.push(pow(a[0].clone(), c)),
            _ => rest.push(t.clone()),
        }
    }
    if inside.len() < 2 {
        return e.clone();
    }
    rest.push(f1("ln", ratsimp(&mul_v(inside))));
    add_v(rest)
}

fn int_factor_string(n: &BigInt) -> String {
    let parts: Vec<String> = factor_int(n).into_iter().map(|(p, k)| if k == 1 { p.to_string() } else { format!("{p}^{k}") }).collect();
    let body = if parts.is_empty() { n.abs().to_string() } else { parts.join("*") };
    if n.is_negative() { format!("-{body}") } else { body }
}

fn factor_poly_expr(e: &Expr) -> Expr {
    let syms = e.free_syms();
    if syms.len() == 1 {
        let x = syms.iter().next().unwrap();
        if let Some(p) = poly::from_expr(e, x) {
            let (c, fs) = poly::factor(&p);
            let xs = sym(x);
            let mut out = vec![Expr::Num(c)];
            out.extend(fs.into_iter().map(|(f, m)| pow(f.to_expr(&xs), int(m as i64))));
            return mul_v(out);
        }
    }
    factor_common(e)
}

/// Multivariate fallback: pull out the common monomial, then difference of squares.
fn factor_common(e: &Expr) -> Expr {
    let Expr::Add(terms) = e else { return e.clone() };
    let mut g: Option<Rat> = None;
    let mut common: Option<Vec<(Expr, i64)>> = None;
    for t in terms {
        let fs = factors_of(t);
        let c = fs.iter().find_map(|f| f.as_rat().cloned()).unwrap_or_else(Rat::one);
        g = Some(match g {
            None => c.abs(),
            Some(g) => Rat::new(g.numer().gcd(c.numer()), g.denom().lcm(c.denom())),
        });
        let mine: Vec<(Expr, i64)> = fs.iter().filter(|f| !f.is_number()).map(base_exp).filter(|(_, k)| *k > 0).collect();
        common = Some(match common {
            None => mine,
            Some(cm) => cm.into_iter().filter_map(|(b, k)| mine.iter().find(|(bb, _)| *bb == b).map(|(_, k2)| (b, k.min(*k2)))).collect(),
        });
    }
    let mut outer = vec![Expr::Num(g.unwrap_or_else(Rat::one))];
    outer.extend(common.unwrap_or_default().into_iter().map(|(b, k)| pow(b, int(k))));
    let outer = mul_v(outer);
    let inner = expand(&div(e.clone(), outer.clone()));
    let inner = diff_squares(&inner).unwrap_or(inner);
    mul(outer, inner)
}

fn perfect_square(t: &Expr) -> Option<Expr> {
    let mut out = Vec::new();
    for f in factors_of(t) {
        match &f {
            Expr::Num(r) if r.is_positive() => {
                let s = sqrt(f.clone());
                s.as_rat()?;
                out.push(s);
            }
            Expr::Pow(b, k) if k.as_i64().is_some_and(|k| k % 2 == 0) => out.push(pow((**b).clone(), int(k.as_i64().unwrap() / 2))),
            _ => return None,
        }
    }
    Some(mul_v(out))
}

fn diff_squares(e: &Expr) -> Option<Expr> {
    let Expr::Add(v) = e else { return None };
    if v.len() != 2 {
        return None;
    }
    for (p, n) in [(&v[0], &v[1]), (&v[1], &v[0])] {
        if n.has_neg_coeff() && !p.has_neg_coeff() {
            let (a, b) = (perfect_square(p)?, perfect_square(&neg(n.clone()))?);
            let (fa, fb) = (sub(a.clone(), b.clone()), add(a, b));
            return Some(mul(diff_squares(&fa).unwrap_or(fa), fb));
        }
    }
    None
}

pub fn factor(e: &Expr) -> Expr {
    match e {
        Expr::Num(r) if r.is_integer() => return Expr::Str(int_factor_string(r.numer())),
        Expr::Num(r) => return Expr::Str(format!("{}/({})", int_factor_string(r.numer()), int_factor_string(r.denom()))),
        Expr::Rel(r, a, b) => return Expr::Rel(*r, Box::new(factor(a)), Box::new(factor(b))),
        Expr::List(_) | Expr::Matrix(_) => return e.map(&mut |c| factor(c)),
        _ => {}
    }
    let (n, d) = together(e);
    let (n, d) = cancel(&n, &d);
    div(factor_poly_expr(&n), factor_poly_expr(&d))
}

pub fn partfrac(e: &Expr, x: &str) -> Option<Expr> {
    let (n, d) = together(e);
    let (pn, pd) = (poly::from_expr(&n, x)?, poly::from_expr(&d, x)?);
    let pf = poly::partial_fractions(&pn, &pd)?;
    let xs = sym(x);
    let mut terms = vec![pf.quotient.to_expr(&xs)];
    for (num, f, k) in pf.terms {
        terms.push(div(num.to_expr(&xs), pow(f.to_expr(&xs), int(k as i64))));
    }
    Some(add_v(terms))
}

/// Rewrite products and powers of sin/cos as sums of sin/cos of combined angles.
pub fn trig_linearize(e: &Expr) -> Expr {
    let e = expand(e);
    let terms: Vec<Expr> = if let Expr::Add(v) = &e { v.clone() } else { vec![e.clone()] };
    add_v(terms.iter().map(linearize_term).collect())
}

fn linearize_term(t: &Expr) -> Expr {
    // (coefficient, optional (is_sin, angle))
    let mut acc: Vec<(Expr, Option<(bool, Expr)>)> = vec![(int(1), None)];
    let mut other = Vec::new();
    for f in factors_of(t) {
        let (b, k) = base_exp(&f);
        let trig = match &b {
            Expr::Func(n, a) if (n == "sin" || n == "cos") && (1..=16).contains(&k) => Some((n == "sin", a[0].clone())),
            _ => None,
        };
        let Some((is_sin, ang)) = trig else {
            other.push(f.clone());
            continue;
        };
        for _ in 0..k {
            let mut next = Vec::new();
            for (c, cur) in acc {
                match cur {
                    None => next.push((c, Some((is_sin, ang.clone())))),
                    Some((s2, a2)) => {
                        let h = mul(c, rat(1, 2));
                        let (p, m) = (add(a2.clone(), ang.clone()), sub(a2.clone(), ang.clone()));
                        match (s2, is_sin) {
                            (true, true) => {
                                next.push((h.clone(), Some((false, m))));
                                next.push((neg(h), Some((false, p))));
                            }
                            (false, false) => {
                                next.push((h.clone(), Some((false, m))));
                                next.push((h, Some((false, p))));
                            }
                            (true, false) => {
                                next.push((h.clone(), Some((true, p))));
                                next.push((h, Some((true, m))));
                            }
                            (false, true) => {
                                next.push((h.clone(), Some((true, p))));
                                next.push((neg(h), Some((true, m))));
                            }
                        }
                    }
                }
            }
            acc = next;
        }
    }
    let rest = mul_v(other);
    add_v(
        acc.into_iter()
            .map(|(c, tr)| {
                let f = match tr {
                    None => int(1),
                    Some((true, a)) => f1("sin", a),
                    Some((false, a)) => f1("cos", a),
                };
                mul_v(vec![c, f, rest.clone()])
            })
            .collect(),
    )
}

pub fn is_zero_expr(e: &Expr) -> bool {
    e.is_zero() || ratsimp(e).is_zero() || is_zero_num(e)
}

