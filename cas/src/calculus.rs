//! Differentiation, integration, limits, series, sums.
use crate::algebra::{ratsimp, together, trig_linearize};
use crate::expr::*;
use crate::numeric::{self, evalc};
use crate::poly::{self, Poly, coeffs_expr};
use crate::simplify::*;
use num_bigint::BigInt;
use num_complex::Complex64 as C;
use num_traits::{One, Signed, Zero};

pub fn diff(e: &Expr, x: &str) -> Expr {
    if !e.has_sym(x) {
        return match e {
            Expr::List(_) | Expr::Matrix(_) => e.map(&mut |c| diff(c, x)),
            _ => int(0),
        };
    }
    match e {
        Expr::Sym(_) => int(1),
        Expr::Add(v) => add_v(v.iter().map(|t| diff(t, x)).collect()),
        Expr::Mul(v) => add_v(
            (0..v.len())
                .map(|i| {
                    let mut f = v.clone();
                    f[i] = diff(&v[i], x);
                    mul_v(f)
                })
                .collect(),
        ),
        Expr::Pow(b, n) => {
            let (b, n) = ((**b).clone(), (**n).clone());
            if !n.has_sym(x) {
                mul_v(vec![n.clone(), pow(b.clone(), sub(n, int(1))), diff(&b, x)])
            } else if !b.has_sym(x) {
                mul_v(vec![e.clone(), f1("ln", b), diff(&n, x)])
            } else {
                mul(e.clone(), add(mul(diff(&n, x), f1("ln", b.clone())), mul_v(vec![n, diff(&b, x), pow(b, int(-1))])))
            }
        }
        Expr::Func(name, a) if a.len() == 1 => {
            let u = a[0].clone();
            let du = diff(&u, x);
            let u2 = || pow(u.clone(), int(2));
            let d = match name.as_str() {
                "sin" => f1("cos", u.clone()),
                "cos" => neg(f1("sin", u.clone())),
                "tan" => pow(f1("cos", u.clone()), int(-2)),
                "asin" => pow(sub(int(1), u2()), rat(-1, 2)),
                "acos" => neg(pow(sub(int(1), u2()), rat(-1, 2))),
                "atan" => pow(add(int(1), u2()), int(-1)),
                "sinh" => f1("cosh", u.clone()),
                "cosh" => f1("sinh", u.clone()),
                "tanh" => pow(f1("cosh", u.clone()), int(-2)),
                "asinh" => pow(add(u2(), int(1)), rat(-1, 2)),
                "acosh" => pow(sub(u2(), int(1)), rat(-1, 2)),
                "atanh" => pow(sub(int(1), u2()), int(-1)),
                "ln" => pow(u.clone(), int(-1)),
                "abs" => f1("sign", u.clone()),
                "erf" => mul_v(vec![int(2), pow(PI, rat(-1, 2)), pow(E, neg(u2()))]),
                "sign" | "floor" | "ceil" | "round" | "ipart" => int(0),
                _ => Expr::Func(format!("{name}'"), vec![u.clone()]),
            };
            mul(d, du)
        }
        Expr::Rel(r, a, b) => Expr::Rel(*r, Box::new(diff(a, x)), Box::new(diff(b, x))),
        Expr::List(_) | Expr::Matrix(_) => e.map(&mut |c| diff(c, x)),
        _ => Expr::Func("diff".into(), vec![e.clone(), sym(x)]),
    }
}

/// e = m*x + b with m, b free of x and m nonzero.
fn linear(e: &Expr, x: &str) -> Option<(Expr, Expr)> {
    let c = coeffs_expr(e, x)?;
    (c.len() == 2).then(|| (c[1].clone(), c[0].clone()))
}

fn split_const(e: &Expr, x: &str) -> (Expr, Expr) {
    let fs: Vec<Expr> = if let Expr::Mul(v) = e { v.clone() } else { vec![e.clone()] };
    let (k, w): (Vec<Expr>, Vec<Expr>) = fs.into_iter().partition(|f| !f.has_sym(x));
    (mul_v(k), mul_v(w))
}

fn is_trig(e: &Expr, name: &str) -> Option<Expr> {
    match e {
        Expr::Func(n, a) if n == name => Some(a[0].clone()),
        _ => None,
    }
}

fn table(e: &Expr, x: &str) -> Option<Expr> {
    let xs = sym(x);
    match e {
        Expr::Sym(_) => Some(mul(rat(1, 2), pow(xs, int(2)))),
        Expr::Pow(b, n) if !n.has_sym(x) => {
            let n = (**n).clone();
            if let Some((m, _)) = linear(b, x) {
                // (m x + b)^n
                if n == int(-1) {
                    return Some(div(f1("ln", f1("abs", (**b).clone())), m));
                }
                let n1 = add(n, int(1));
                return Some(div(pow((**b).clone(), n1.clone()), mul(n1, m)));
            }
            if let Expr::Func(name, a) = &**b {
                let (m, _) = linear(&a[0], x)?;
                let u = a[0].clone();
                let k = n.as_i64()?;
                let (s, c) = (f1("sin", u.clone()), f1("cos", u.clone()));
                return match (name.as_str(), k) {
                    ("cos", -2) => Some(div(f1("tan", u), m)),
                    ("sin", -2) => Some(neg(div(c, mul(s, m)))),
                    ("cos", -1) => Some(div(f1("ln", f1("abs", add(pow(c, int(-1)), f1("tan", u)))), m)),
                    ("sin", -1) => Some(div(f1("ln", f1("abs", f1("tan", mul(rat(1, 2), u)))), m)),
                    ("tan", k) if k >= 2 => {
                        let lower = integrate(&pow((**b).clone(), int(k - 2)), x, 4)?;
                        Some(sub(div(pow((**b).clone(), int(k - 1)), mul(int(k - 1), m)), lower))
                    }
                    ("cosh", -2) => Some(div(f1("tanh", u), m)),
                    _ => None,
                };
            }
            // sqrt(c - a x^2), sqrt(a x^2 + c)
            if n == rat(1, 2) {
                let c = coeffs_expr(b, x)?;
                if c.len() == 3 && c[1].is_zero() {
                    let (a2, c0) = (numeric::real_const(&c[2])?, numeric::real_const(&c[0])?);
                    let first = mul_v(vec![rat(1, 2), xs.clone(), e.clone()]);
                    if a2 < 0.0 && c0 > 0.0 {
                        let a = neg(c[2].clone());
                        return Some(add(first, mul_v(vec![rat(1, 2), c[0].clone(), pow(sqrt(a.clone()), int(-1)), f1("asin", mul(xs, sqrt(div(a, c[0].clone()))))])));
                    }
                    if a2 > 0.0 {
                        let a = sqrt(c[2].clone());
                        return Some(add(first, mul_v(vec![rat(1, 2), c[0].clone(), pow(a.clone(), int(-1)), f1("ln", f1("abs", add(mul(a, xs), e.clone())))])));
                    }
                }
            }
            // 1/sqrt(c - a x^2), 1/sqrt(a x^2 + c)
            if n == rat(-1, 2) {
                let c = coeffs_expr(b, x)?;
                if c.len() == 3 && c[1].is_zero() {
                    let (a2, c0) = (numeric::real_const(&c[2])?, numeric::real_const(&c[0])?);
                    if a2 < 0.0 && c0 > 0.0 {
                        let a = neg(c[2].clone());
                        return Some(div(f1("asin", mul(xs, sqrt(div(a.clone(), c[0].clone())))), sqrt(a)));
                    }
                    if a2 > 0.0 {
                        let a = sqrt(c[2].clone());
                        return Some(div(f1("ln", f1("abs", add(mul(a.clone(), xs), sqrt((**b).clone())))), a));
                    }
                }
            }
            None
        }
        Expr::Pow(a, u) if !a.has_sym(x) => {
            if let Some((m, _)) = linear(u, x) {
                return Some(div(e.clone(), mul(f1("ln", (**a).clone()), m)));
            }
            // e^(c0 - a x^2) -> erf
            let c = coeffs_expr(u, x)?;
            if **a == E && c.len() == 3 && c[1].is_zero() && numeric::real_const(&c[2])? < 0.0 {
                let ra = sqrt(neg(c[2].clone()));
                return Some(mul_v(vec![pow(E, c[0].clone()), sqrt(PI), rat(1, 2), pow(ra.clone(), int(-1)), f1("erf", mul(ra, xs))]));
            }
            None
        }
        Expr::Func(name, a) if a.len() == 1 => {
            let u = a[0].clone();
            let (m, _) = linear(&u, x)?;
            let u2 = pow(u.clone(), int(2));
            let r = match name.as_str() {
                "sin" => neg(f1("cos", u)),
                "cos" => f1("sin", u),
                "tan" => neg(f1("ln", f1("abs", f1("cos", u)))),
                "sinh" => f1("cosh", u),
                "cosh" => f1("sinh", u),
                "tanh" => f1("ln", f1("cosh", u)),
                "ln" => sub(mul(u.clone(), f1("ln", u.clone())), u),
                "asin" => add(mul(u.clone(), f1("asin", u)), sqrt(sub(int(1), u2))),
                "acos" => sub(mul(u.clone(), f1("acos", u)), sqrt(sub(int(1), u2))),
                "atan" => sub(mul(u.clone(), f1("atan", u)), mul(rat(1, 2), f1("ln", add(int(1), u2)))),
                "abs" => mul_v(vec![rat(1, 2), u.clone(), f1("abs", u)]),
                _ => return None,
            };
            Some(div(r, m))
        }
        _ => None,
    }
}

fn rational(e: &Expr, x: &str) -> Option<Expr> {
    let (n, d) = together(e);
    if !d.has_sym(x) {
        return None;
    }
    let (pn, pd) = (poly::from_expr(&n, x)?, poly::from_expr(&d, x)?);
    let pf = poly::partial_fractions(&pn, &pd)?;
    let xs = sym(x);
    let r = |v: i64| Rat::from_integer(BigInt::from(v));
    let mut out = Vec::new();
    // quotient polynomial
    let q = Poly::new(std::iter::once(Rat::zero()).chain(pf.quotient.0.iter().enumerate().map(|(i, c)| c / r(i as i64 + 1))).collect());
    out.push(q.to_expr(&xs));
    for (num, f, k) in pf.terms {
        let fe = f.to_expr(&xs);
        match (f.deg(), k) {
            (1, 1) => out.push(mul(Expr::Num(num.coef(0) / f.coef(1)), f1("ln", f1("abs", fe)))),
            (1, k) => out.push(mul(Expr::Num(num.coef(0) / (f.coef(1) * r(1 - k as i64))), pow(fe, int(1 - k as i64)))),
            (2, 1) => {
                // (B x + C) / (a x^2 + b x + c), irreducible
                let (a, b, c) = (f.coef(2), f.coef(1), f.coef(0));
                let (bb, cc) = (num.coef(1), num.coef(0));
                let k1 = &bb / (r(2) * &a);
                let k2 = &cc - &k1 * &b;
                out.push(mul(Expr::Num(k1), f1("ln", fe)));
                let disc = r(4) * &a * &c - &b * &b;
                let sq = sqrt(Expr::Num(disc));
                let lin = Poly::new(vec![b.clone(), r(2) * &a]).to_expr(&xs);
                out.push(mul_v(vec![Expr::Num(r(2) * k2), pow(sq.clone(), int(-1)), f1("atan", ratsimp(&div(lin, sq.clone())))]));
            }
            _ => return None,
        }
    }
    Some(add_v(out))
}

fn subexprs(e: &Expr, x: &str, out: &mut Vec<Expr>) {
    if !e.has_sym(x) || matches!(e, Expr::Sym(_)) {
        return;
    }
    if !out.contains(e) {
        out.push(e.clone());
    }
    for c in e.children() {
        subexprs(c, x, out);
    }
}

fn substitution(e: &Expr, x: &str, depth: usize) -> Option<Expr> {
    let mut cands = Vec::new();
    for c in e.children() {
        subexprs(c, x, &mut cands);
    }
    // also inner arguments one level down from the whole expression
    let t = sym("@u");
    // linear inner expression under a root or function: x = (t - b)/m
    for u in cands.iter().filter(|u| matches!(u, Expr::Add(_))) {
        if let Some((m, b)) = linear(u, x) {
            let r = expand(&div(subst(e, x, &div(sub(t.clone(), b), m.clone())), m));
            if r.size() < 200 {
                if let Some(f) = integrate(&r, "@u", depth - 1) {
                    return Some(subst(&f, "@u", u));
                }
            }
        }
    }
    for u in cands {
        if linear(&u, x).is_some() && !matches!(u, Expr::Add(_)) {
            continue;
        }
        let du = diff(&u, x);
        if du.is_zero() {
            continue;
        }
        let q = div(e.clone(), du);
        let mut r = replace_gen(&q, &u, &t);
        if r.has_sym(x) {
            r = replace_gen(&ratsimp_light(&q), &u, &t);
        }
        if r.has_sym(x) {
            continue;
        }
        if let Some(f) = integrate(&r, "@u", depth - 1) {
            return Some(subst(&f, "@u", &u));
        }
    }
    None
}

fn ratsimp_light(e: &Expr) -> Expr {
    // structural cancellation is automatic; a rational pass catches polynomial factors
    let r = ratsimp(e);
    if r.size() < e.size() { r } else { e.clone() }
}

fn by_parts(e: &Expr, x: &str, depth: usize) -> Option<Expr> {
    let fs: Vec<Expr> = if let Expr::Mul(v) = e { v.clone() } else { vec![e.clone()] };
    let is_log_like = |f: &Expr| match f {
        Expr::Func(n, _) => ["ln", "asin", "acos", "atan", "asinh", "acosh", "atanh"].contains(&n.as_str()),
        Expr::Pow(b, k) => matches!(&**b, Expr::Func(n, _) if n == "ln") && k.as_i64().is_some_and(|k| k > 0),
        _ => false,
    };
    // u = log-like factor, dv = the rest
    if let Some(i) = fs.iter().position(is_log_like) {
        let u = fs[i].clone();
        let dv = mul_v(fs.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, f)| f.clone()).collect());
        if poly_like(&dv, x) {
            let v = integrate(&dv, x, depth - 1)?;
            let rest = integrate(&ratsimp_light(&mul(v.clone(), diff(&u, x))), x, depth - 1)?;
            return Some(sub(mul(u, v), rest));
        }
    }
    // u = polynomial, dv = exp/trig of a linear argument
    let is_osc = |f: &Expr| match f {
        Expr::Func(n, a) => ["sin", "cos", "sinh", "cosh"].contains(&n.as_str()) && linear(&a[0], x).is_some(),
        Expr::Pow(b, u) => !b.has_sym(x) && linear(u, x).is_some(),
        _ => false,
    };
    let (osc, rest): (Vec<Expr>, Vec<Expr>) = fs.into_iter().partition(is_osc);
    if osc.len() == 1 && !rest.is_empty() {
        let u = mul_v(rest.clone());
        if coeffs_expr(&u, x).is_some_and(|c| c.len() >= 2) {
            let v = table(&osc[0], x)?;
            let rest = integrate(&mul(v.clone(), diff(&u, x)), x, depth - 1)?;
            return Some(sub(mul(u, v), rest));
        }
    }
    // e^(a x) * sin/cos(b x)
    if osc.len() == 2 && rest.is_empty() {
        for (p, q) in [(&osc[0], &osc[1]), (&osc[1], &osc[0])] {
            if let (Expr::Pow(_, _), Expr::Func(n, arg)) = (p, q) {
                if n == "sin" || n == "cos" {
                    let ex = table(p, x)?; // p / a
                    let a = ratsimp(&div(p.clone(), ex));
                    let (b, _) = linear(&arg[0], x)?;
                    let den = add(pow(a.clone(), int(2)), pow(b.clone(), int(2)));
                    let (s, c) = (f1("sin", arg[0].clone()), f1("cos", arg[0].clone()));
                    let inner = if n == "sin" { sub(mul(a, s), mul(b, c)) } else { add(mul(a, c), mul(b, s)) };
                    return Some(div(mul(p.clone(), inner), den));
                }
            }
        }
    }
    None
}

fn poly_like(e: &Expr, x: &str) -> bool {
    // x^n for any constant n, or a polynomial
    coeffs_expr(e, x).is_some() || matches!(e, Expr::Pow(b, n) if matches!(&**b, Expr::Sym(s) if s == x) && !n.has_sym(x))
}

fn has_trig_product(e: &Expr) -> bool {
    let trig = |f: &Expr| is_trig(f, "sin").is_some() || is_trig(f, "cos").is_some();
    e.has(&|n| match n {
        Expr::Pow(b, k) => trig(b) && k.as_i64().is_some_and(|k| k >= 2),
        Expr::Mul(v) => v.iter().filter(|f| trig(f) || matches!(f, Expr::Pow(b, _) if trig(b))).count() >= 2,
        _ => false,
    })
}

/// Antiderivative with respect to x, or None when no rule applies.
pub fn integrate(e: &Expr, x: &str, depth: usize) -> Option<Expr> {
    if !e.has_sym(x) {
        return Some(mul(e.clone(), sym(x)));
    }
    if depth == 0 {
        return None;
    }
    if let Expr::Add(v) = e {
        return Some(add_v(v.iter().map(|t| integrate(t, x, depth)).collect::<Option<Vec<Expr>>>()?));
    }
    let (k, w) = split_const(e, x);
    if !k.is_one() {
        return Some(mul(k, integrate(&w, x, depth)?));
    }
    if let Some(r) = table(e, x) {
        return Some(r);
    }
    if let Some(r) = rational(e, x) {
        return Some(r);
    }
    if let Some(r) = by_parts(e, x, depth) {
        return Some(r);
    }
    if has_trig_product(e) {
        let lin = trig_linearize(e);
        if lin != *e {
            if let Some(r) = integrate(&lin, x, depth - 1) {
                return Some(r);
            }
        }
    }
    if let Some(r) = substitution(e, x, depth) {
        return Some(r);
    }
    let ex = expand(e);
    if ex != *e {
        return integrate(&ex, x, depth - 1);
    }
    None
}

/// Check F' = f numerically at a few points (other symbols get fixed sample values).
fn verify(f: &Expr, big_f: &Expr, x: &str) -> bool {
    let d = diff(big_f, x);
    let mut syms = f.free_syms();
    syms.extend(d.free_syms());
    let mut good = 0;
    for (i, xv) in [0.37, 1.13, 2.71, -0.83, 5.2, 0.05].iter().enumerate() {
        let env: Vec<(&str, C)> = syms.iter().enumerate().map(|(j, s)| (s.as_str(), C::new(if s == x { *xv } else { 0.7 + 0.31 * j as f64 + 0.1 * i as f64 }, 0.0))).collect();
        let (Some(a), Some(b)) = (evalc(f, &env), evalc(&d, &env)) else { continue };
        if !a.re.is_finite() || !b.re.is_finite() || a.im.abs() > 1e-9 || b.im.abs() > 1e-9 || !a.im.is_finite() || !b.im.is_finite() {
            continue;
        }
        if (a - b).norm() > 1e-6 * (1.0 + a.norm()) {
            return false;
        }
        good += 1;
    }
    good > 0
}

pub fn antiderivative(e: &Expr, x: &str) -> Option<Expr> {
    let r = integrate(e, x, 6)?;
    verify(e, &r, x).then_some(r)
}

fn at(f: &Expr, x: &str, a: &Expr) -> Expr {
    let v = subst(f, x, a);
    let bad = v.is_undef() || v.has(&|n| matches!(n, Expr::Const(Const::Inf))) || (v.is_constant() && numeric::eval_const(&v).is_none_or(|z| !z.re.is_finite()));
    if bad || *a == INF || *a == neg(INF) { limit(f, x, a, 0).unwrap_or(v) } else { v }
}

/// Definite integral. Returns (value, is_approximate).
pub fn definite(e: &Expr, x: &str, a: &Expr, b: &Expr) -> Result<(Expr, bool), String> {
    let bound = |v: &Expr| -> Option<f64> {
        if *v == INF {
            Some(f64::INFINITY)
        } else if *v == neg(INF) {
            Some(f64::NEG_INFINITY)
        } else if v.is_constant() {
            numeric::real_const(v)
        } else {
            None
        }
    };
    let numeric_val = match (bound(a), bound(b)) {
        (Some(lo), Some(hi)) if e.free_syms().iter().all(|s| s == x) => {
            // split at interior poles so a divergent piece cannot cancel against another
            let (n, d) = together(e);
            let (l, h) = (lo.min(hi), lo.max(hi));
            let mut cuts = vec![l];
            if d.has_sym(x) {
                let (sl, sh) = (l.max(-1e4), h.min(1e4));
                for p in numeric::find_roots(&|t| numeric::evalf(&d, x, t), sl, sh, 4000) {
                    if p > l + 1e-12 && p < h - 1e-12 && numeric::evalf(&n, x, p).abs() > 1e-9 {
                        cuts.push(p);
                    }
                }
            }
            cuts.push(h);
            let (mut v, mut err) = (0.0, 0.0);
            for w in cuts.windows(2) {
                let (pv, pe) = numeric::quad(&|t| numeric::evalf(e, x, t), w[0], w[1]);
                v += pv;
                let ok = pe <= 1e-7 * (1.0 + pv.abs());
                err += if ok { pe } else { f64::INFINITY };
                if !ok && cuts.len() > 2 {
                    // a piece next to a pole that does not settle: divergent
                    v = f64::NAN;
                }
            }
            Some((if lo > hi { -v } else { v }, err))
        }
        _ => None,
    };
    let converged = |(v, err): (f64, f64)| v.is_finite() && err <= 1e-7 * (1.0 + v.abs());
    if let Some(big_f) = antiderivative(e, x) {
        let exact = sub(at(&big_f, x, b), at(&big_f, x, a));
        let exact = if exact.is_constant() && exact.size() < 80 { expand(&exact) } else { exact };
        let Some(nv) = numeric_val else { return Ok((exact, false)) };
        if let Some(ev) = numeric::real_const(&exact).filter(|v| v.is_finite()) {
            // agreement guards against a discontinuous antiderivative; a kink only slows the quadrature
            let tol = if converged(nv) { 1e-6 } else { 1e-3 };
            if nv.0.is_finite() && (ev - nv.0).abs() <= tol * (1.0 + nv.0.abs()) {
                return Ok((exact, false));
            }
        }
        if converged(nv) {
            return Ok((Expr::Float(nv.0), true));
        }
        if exact.is_undef() || exact.has(&|n| matches!(n, Expr::Const(Const::Inf))) {
            return Ok((exact, false));
        }
        return Err("integral does not converge on this interval".into());
    }
    match numeric_val {
        Some(nv) if converged(nv) => Ok((Expr::Float(nv.0), true)),
        Some(_) => Err("integral does not converge".into()),
        None => Err("no closed form found; give numeric bounds for a numeric result".into()),
    }
}

fn finite_value(v: &Expr) -> bool {
    if v.is_undef() || v.has(&|n| matches!(n, Expr::Const(Const::Inf))) {
        return false;
    }
    !v.is_constant() || numeric::eval_const(v).is_some_and(|z| z.re.is_finite() && z.im.is_finite())
}

fn numeric_limit(e: &Expr, x: &str, a: f64, dir: i32) -> Option<Expr> {
    if !e.free_syms().iter().all(|s| s == x) {
        return None;
    }
    let side = |s: f64| -> Option<f64> {
        let mut vals: Vec<f64> = Vec::new();
        for k in 2..=9 {
            let xv = if a.is_infinite() { a.signum() * 10f64.powi(k + 1) } else { a + s * 10f64.powi(-k) };
            let v = numeric::evalf(e, x, xv);
            if v.is_nan() {
                break;
            }
            vals.push(v);
        }
        let n = vals.len();
        if n < 3 {
            return None;
        }
        let last = vals[n - 1];
        if last.is_infinite() {
            return Some(last);
        }
        // steadily growing magnitude: diverges
        if last.abs() > 1e6 && vals.windows(2).all(|w| w[1].abs() >= w[0].abs() && w[1].signum() == w[0].signum()) {
            return Some(f64::INFINITY * last.signum());
        }
        // steadily shrinking magnitude: tends to zero
        if last.abs() < 1e-4 && n >= 5 && vals[n - 5..].windows(2).all(|w| w[1].abs() < w[0].abs()) {
            return Some(0.0);
        }
        ((last - vals[n - 2]).abs() <= 1e-6 * (1.0 + last.abs())).then_some(last)
    };
    let (l, r) = (side(-1.0), side(1.0));
    let v = match dir {
        d if d < 0 => l?,
        d if d > 0 => r?,
        _ if a.is_infinite() => r?,
        _ => {
            let (l, r) = (l?, r?);
            if l.is_infinite() && r.is_infinite() && l == r {
                l
            } else if l.is_finite() && r.is_finite() && (l - r).abs() <= 1e-5 * (1.0 + l.abs()) {
                0.5 * (l + r)
            } else {
                return Some(UNDEF);
            }
        }
    };
    Some(if v.is_infinite() {
        if v > 0.0 { INF } else { neg(INF) }
    } else if v.abs() < 1e-9 {
        Expr::Float(0.0)
    } else {
        Expr::Float(v)
    })
}

fn cancelled_changed(n: &Expr, d: &Expr, reduced: &Expr) -> bool {
    div(n.clone(), d.clone()) != *reduced
}

/// Limit of e as x -> a. dir: -1 from the left, 1 from the right, 0 two-sided.
pub fn limit(e: &Expr, x: &str, a: &Expr, dir: i32) -> Option<Expr> {
    limit_d(e, x, a, dir, 0)
}

fn limit_d(e: &Expr, x: &str, a: &Expr, dir: i32, depth: usize) -> Option<Expr> {
    if !e.has_sym(x) {
        return Some(e.clone());
    }
    let neg_inf = neg(INF);
    let av = if *a == INF {
        f64::INFINITY
    } else if *a == neg_inf {
        f64::NEG_INFINITY
    } else {
        numeric::real_const(a).unwrap_or(f64::NAN)
    };
    let direct = subst(e, x, a);
    if direct == INF || direct == neg_inf {
        return Some(direct);
    }
    if finite_value(&direct) {
        // a removable form such as 0^0 can still slip through; cross-check numerically when possible
        return Some(direct);
    }
    if depth > 6 {
        return numeric_limit(e, x, av, dir);
    }
    if let Expr::Add(v) = e {
        let parts: Option<Vec<Expr>> = v.iter().map(|t| limit_d(t, x, a, dir, depth + 1)).collect();
        if let Some(p) = parts {
            let s = add_v(p);
            if !s.is_undef() {
                return Some(s);
            }
        }
    }
    // u^v with x in the exponent: exp(lim v ln u)
    if let Expr::Pow(b, ex) = e {
        if ex.has_sym(x) && b.has_sym(x) {
            let l = limit_d(&mul((**ex).clone(), f1("ln", (**b).clone())), x, a, dir, depth + 1)?;
            let r = pow(E, l);
            if !r.is_undef() {
                return Some(r);
            }
        }
    }
    // quotient: L'Hopital on 0/0 and inf/inf
    let (n, d) = together(e);
    if d.has_sym(x) && poly::from_expr(&n, x).is_some() && poly::from_expr(&d, x).is_some() {
        let (cn, cd) = crate::algebra::cancel(&n, &d);
        let reduced = div(cn, cd);
        if reduced != *e && reduced.size() < e.size() + 4 && (cancelled_changed(&n, &d, &reduced)) {
            return limit_d(&reduced, x, a, dir, depth + 1);
        }
    }
    if d.has_sym(x) {
        let (ln, ld) = (limit_d(&n, x, a, dir, depth + 1), limit_d(&d, x, a, dir, depth + 1));
        if let (Some(ln), Some(ld)) = (&ln, &ld) {
            let is_inf = |v: &Expr| *v == INF || *v == neg_inf;
            if (crate::simplify::is_zero_num(ln) && crate::simplify::is_zero_num(ld)) || (is_inf(ln) && is_inf(ld)) {
                let (dn, dd) = (diff(&n, x), diff(&d, x));
                let q = ratsimp_light(&div(dn, dd));
                if q.size() < 400 {
                    if let Some(r) = limit_d(&q, x, a, dir, depth + 1) {
                        if !matches!(r, Expr::Float(_)) {
                            return Some(r);
                        }
                    }
                }
            } else if crate::simplify::is_zero_num(ld) && finite_value(ln) && !crate::simplify::is_zero_num(ln) {
                return numeric_limit(e, x, av, dir);
            } else if finite_value(ld) && !crate::simplify::is_zero_num(ld) {
                let r = div(ln.clone(), ld.clone());
                if !r.is_undef() {
                    return Some(r);
                }
            } else if is_inf(ld) && finite_value(ln) {
                return Some(int(0));
            }
        }
    } else if let Expr::Mul(v) = &n {
        // 0 * inf: move a vanishing factor into the denominator
        let lims: Vec<Option<Expr>> = v.iter().map(|f| limit_d(f, x, a, dir, depth + 1)).collect();
        if let Some(i) = lims.iter().position(|l| l.as_ref().is_some_and(crate::simplify::is_zero_num)) {
            if lims.iter().any(|l| l.as_ref().is_some_and(|l| *l == INF || *l == neg_inf)) && depth < 4 {
                let _ = i;
                let is_fn = |f: &Expr| f.has(&|n| matches!(n, Expr::Func(..)));
                let (top, alg): (Vec<Expr>, Vec<Expr>) = v.iter().cloned().partition(is_fn);
                let (rest, zero_f) = if !top.is_empty() && !alg.is_empty() {
                    (mul_v(top), mul_v(alg))
                } else {
                    (mul_v(v.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, f)| f.clone()).collect()), v[i].clone())
                };
                let q = ratsimp_light(&div(diff(&rest, x), diff(&pow(zero_f, int(-1)), x)));
                if let Some(r) = limit_d(&div(q, d.clone()), x, a, dir, depth + 2) {
                    if !matches!(r, Expr::Float(_)) {
                        return Some(r);
                    }
                }
            }
        }
    }
    numeric_limit(e, x, av, dir)
}

pub fn taylor(e: &Expr, x: &str, a: &Expr, n: usize) -> Option<Expr> {
    let mut d = e.clone();
    let mut terms = Vec::new();
    let mut fact = BigInt::one();
    let dx = sub(sym(x), a.clone());
    for k in 0..=n {
        if k > 0 {
            d = diff(&d, x);
            fact *= k;
            if d.size() > 20_000 {
                d = ratsimp(&d);
            }
        }
        let c = {
            let v = subst(&d, x, a);
            if finite_value(&v) { v } else { limit(&d, x, a, 0)? }
        };
        let c = if c.is_constant() && c.size() < 60 { expand(&c) } else { c };
        terms.push(mul_v(vec![c, Expr::Num(Rat::new(BigInt::one(), fact.clone())), pow(dx.clone(), int(k as i64))]));
    }
    Some(add_v(terms))
}

fn bernoulli_plus(n: usize) -> Vec<Rat> {
    // B_0..B_n with B_1 = +1/2
    let r = |v: i64| Rat::from_integer(BigInt::from(v));
    let mut b: Vec<Rat> = Vec::new();
    for m in 0..=n {
        let mut s = Rat::zero();
        let mut binom = Rat::one();
        for (k, bk) in b.iter().enumerate() {
            s += &binom * bk;
            binom = binom * r(m as i64 + 1 - k as i64) / r(k as i64 + 1);
        }
        b.push(if m == 0 { Rat::one() } else { -s / r(m as i64 + 1) });
    }
    if n >= 1 {
        b[1] = Rat::new(BigInt::one(), BigInt::from(2));
    }
    b
}

/// Sum_{k=1}^{n} k^p as a polynomial in n (Faulhaber).
fn power_sum(p: usize, n: &Expr) -> Expr {
    let r = |v: i64| Rat::from_integer(BigInt::from(v));
    let b = bernoulli_plus(p);
    let mut binom = Rat::one();
    let mut terms = Vec::new();
    for j in 0..=p {
        terms.push(mul(Expr::Num(&binom * &b[j] / r(p as i64 + 1)), pow(n.clone(), int((p + 1 - j) as i64))));
        binom = binom * r(p as i64 + 1 - j as i64) / r(j as i64 + 1);
    }
    add_v(terms)
}

pub fn sum(e: &Expr, k: &str, a: &Expr, b: &Expr, product: bool) -> Result<Expr, String> {
    if let (Some(lo), Some(hi)) = (a.as_i64(), b.as_i64()) {
        if hi - lo > 1_000_000 {
            return Err("too many terms".into());
        }
        let mut items = Vec::with_capacity((hi - lo + 1).max(0) as usize);
        for i in lo..=hi {
            items.push(subst(e, k, &int(i)));
        }
        return Ok(if product { mul_v(items) } else { add_v(items) });
    }
    let unevaluated = || Expr::Func(if product { "product" } else { "sum" }.into(), vec![e.clone(), sym(k), a.clone(), b.clone()]);
    if product {
        return Ok(unevaluated());
    }
    if !e.has_sym(k) {
        return Ok(mul(e.clone(), add_v(vec![b.clone(), neg(a.clone()), int(1)])));
    }
    // polynomial summand
    if *b != INF {
        if let Some(c) = coeffs_expr(e, k) {
            if c.len() <= 30 {
                let am1 = sub(a.clone(), int(1));
                let total = add_v(
                    c.iter()
                        .enumerate()
                        .map(|(p, coef)| {
                            let s = if p == 0 { sub(b.clone(), am1.clone()) } else { sub(power_sum(p, b), power_sum(p, &am1)) };
                            mul(coef.clone(), s)
                        })
                        .collect(),
                );
                return Ok(crate::algebra::factor(&expand(&total)));
            }
        }
    }
    // geometric: c * r^k
    let (c, w) = split_const(e, k);
    if let Expr::Pow(base, ex) = &w {
        if !base.has_sym(k) {
            if let Some((m, off)) = linear(ex, k) {
                let r = pow((**base).clone(), m);
                let lead = mul(c, pow((**base).clone(), off));
                let first = pow(r.clone(), a.clone());
                if *b == INF {
                    if numeric::real_const(&r).is_some_and(|v| v.abs() < 1.0) {
                        return Ok(ratsimp(&mul(lead, div(first, sub(int(1), r)))));
                    }
                    return Err("series diverges or ratio is not a known constant".into());
                }
                let count = add_v(vec![b.clone(), neg(a.clone()), int(1)]);
                let s = mul(lead, div(mul(first, sub(pow(r.clone(), count), int(1))), sub(r.clone(), int(1))));
                return Ok(if r.is_number() { expand(&s) } else { s });
            }
        }
    }
    // zeta at even integers: sum of c/k^(2m) from 1
    if *b == INF && a.is_one() {
        if let Expr::Pow(base, ex) = &w {
            if matches!(&**base, Expr::Sym(s) if s == k) {
                if let Some(p) = ex.as_i64().filter(|p| *p <= -2 && p % 2 == 0 && *p >= -60) {
                    let p = (-p) as usize;
                    let bern = bernoulli_plus(p)[p].clone();
                    let fact = (1..=p).fold(BigInt::one(), |f, j| f * j);
                    let coef = bern.abs() * Rat::from_integer(BigInt::from(2).pow(p as u32 - 1)) / Rat::from_integer(fact);
                    return Ok(mul_v(vec![c, Expr::Num(coef), pow(PI, int(p as i64))]));
                }
            }
        }
    }
    // numeric infinite sum: partial sum plus an Euler-Maclaurin tail
    if *b == INF && e.free_syms().iter().all(|s| s == k) {
        if let Some(lo) = a.as_i64() {
            let f = |t: f64| numeric::evalf(e, k, t);
            let n0 = lo.max(1) + 2000;
            let (tail, terr) = numeric::quad(&f, n0 as f64, f64::INFINITY);
            if tail.is_finite() && terr <= 1e-10 * (1.0 + tail.abs()) && f(n0 as f64).abs() < 1.0 {
                let head: f64 = (lo..n0).map(|i| f(i as f64)).sum();
                let d1 = numeric::nderiv(&f, n0 as f64);
                return Ok(Expr::Float(head + tail + 0.5 * f(n0 as f64) - d1 / 12.0));
            }
            let term = |i: i64| numeric::evalf(e, k, i as f64);
            let partial = |n: i64| (lo..lo + n).map(term).sum::<f64>();
            // Richardson on partial sums at N, 2N, 4N assuming a power-law tail
            let (s1, s2, s3) = (partial(20_000), partial(40_000), partial(80_000));
            let (d1, d2) = (s2 - s1, s3 - s2);
            if d2.abs() < 1e-14 * (1.0 + s3.abs()) {
                return Ok(Expr::Float(s3));
            }
            let ratio = d2 / d1;
            if ratio.abs() < 0.9 && ratio.is_finite() {
                return Ok(Expr::Float(s3 + d2 * ratio / (1.0 - ratio)));
            }
            return Err("series does not converge".into());
        }
    }
    Ok(unevaluated())
}

/// Classic RK4 for y' = f(x, y).
pub fn ode_rk4(f: &Expr, x: &str, y: &str, x0: f64, y0: f64, x1: f64, steps: usize) -> Vec<(f64, f64)> {
    let g = |xv: f64, yv: f64| evalc(f, &[(x, C::new(xv, 0.0)), (y, C::new(yv, 0.0))]).map_or(f64::NAN, |z| z.re);
    let h = (x1 - x0) / steps as f64;
    let (mut xv, mut yv) = (x0, y0);
    let mut out = vec![(xv, yv)];
    for _ in 0..steps {
        let k1 = g(xv, yv);
        let k2 = g(xv + h / 2.0, yv + h * k1 / 2.0);
        let k3 = g(xv + h / 2.0, yv + h * k2 / 2.0);
        let k4 = g(xv + h, yv + h * k3);
        yv += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        xv += h;
        out.push((xv, yv));
    }
    out
}

