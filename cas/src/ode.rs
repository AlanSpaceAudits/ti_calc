//! Symbolic ordinary differential equations: first order (linear, separable, Bernoulli, exact)
//! and second order linear with constant coefficients.
use crate::algebra::{self, ratsimp, together};
use crate::calculus::{antiderivative, diff};
use crate::expr::*;
use crate::numeric::evalc;
use crate::poly::coeffs_expr;
use crate::simplify::*;
use crate::solve;
use num_complex::Complex64 as C;

fn strip_abs(e: &Expr) -> Expr {
    match e {
        Expr::Func(n, a) if n == "abs" => strip_abs(&a[0]),
        Expr::Add(v) => add_v(v.iter().map(strip_abs).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(strip_abs).collect()),
        Expr::Pow(b, x) => pow(strip_abs(b), strip_abs(x)),
        Expr::Func(n, v) => func(n, v.iter().map(strip_abs).collect()),
        e => e.clone(),
    }
}

fn int(e: &Expr, x: &str) -> Result<Expr, String> {
    antiderivative(e, x).or_else(|| antiderivative(&algebra::simplify(e), x)).ok_or_else(|| "an integral in the solution has no closed form".to_string())
}

fn tidy(e: Expr) -> Expr {
    let ex = expand(&e);
    let s = algebra::simplify(&ex);
    if s.size() <= ex.size() { s } else { ex }
}

/// y' + p*y = q
fn linear_first(p: &Expr, q: &Expr, x: &str, c: &Expr) -> Result<Expr, String> {
    let mu = pow(E, strip_abs(&int(p, x)?));
    let rhs = int(&tidy(mul(mu.clone(), q.clone())), x)?;
    Ok(tidy(div(add(rhs, c.clone()), mu)))
}

/// Terms of f as coefficient(x) * y^k.
fn powers_of(f: &Expr, y: &str) -> Option<Vec<(Expr, Expr)>> {
    let f = expand(f);
    let terms: Vec<Expr> = if let Expr::Add(v) = &f { v.clone() } else { vec![f.clone()] };
    let mut out: Vec<(Expr, Expr)> = Vec::new();
    for t in terms {
        let fs: Vec<Expr> = if let Expr::Mul(v) = &t { v.clone() } else { vec![t.clone()] };
        let (mut k, mut coef) = (int_e(0), Vec::new());
        for fac in fs {
            match &fac {
                Expr::Sym(s) if s == y => k = add(k, int_e(1)),
                Expr::Pow(b, n) if matches!(&**b, Expr::Sym(s) if s == y) && n.as_rat().is_some() => k = add(k, (**n).clone()),
                other if !other.has_sym(y) => coef.push(other.clone()),
                _ => return None,
            }
        }
        match out.iter_mut().find(|(kk, _)| *kk == k) {
            Some(slot) => slot.1 = add(slot.1.clone(), mul_v(coef)),
            None => out.push((k, mul_v(coef))),
        }
    }
    Some(out)
}

fn int_e(n: i64) -> Expr {
    crate::expr::int(n)
}

fn explicit(implicit: &Expr, y: &str) -> Option<Vec<Expr>> {
    match solve::solve(implicit, y, false) {
        Ok((solve::Solution::Set(v), false)) if !v.is_empty() && v.iter().all(|s| s.size() < 120) => Some(v),
        _ => None,
    }
}

enum Sol {
    Explicit(Vec<Expr>),
    /// psi(x, y) = 0
    Implicit(Expr),
}

fn first_order(f: &Expr, x: &str, y: &str) -> Result<Sol, String> {
    let c1 = sym("C1");
    // linear
    if let Some(c) = coeffs_expr(f, y) {
        if c.len() <= 2 {
            let p = neg(c.get(1).cloned().unwrap_or(int_e(0)));
            return Ok(Sol::Explicit(vec![linear_first(&p, &c[0], x, &c1)?]));
        }
    }
    // separable: every factor depends on x only or on y only
    let fac = algebra::factor(f);
    let fs: Vec<Expr> = if let Expr::Mul(v) = &fac { v.clone() } else { vec![fac.clone()] };
    if fs.iter().all(|t| !t.has_sym(x) || !t.has_sym(y)) {
        let (h, g): (Vec<Expr>, Vec<Expr>) = fs.into_iter().partition(|t| t.has_sym(y));
        let lhs = int(&pow(mul_v(h), int_e(-1)), y)?;
        let rhs = int(&mul_v(g), x)?;
        let psi = sub(sub(lhs, rhs), c1.clone());
        return Ok(match explicit(&psi, y) {
            Some(v) => Sol::Explicit(v.into_iter().map(tidy).collect()),
            None => Sol::Implicit(psi),
        });
    }
    // Bernoulli: y' = a(x)*y + b(x)*y^n
    if let Some(pw) = powers_of(f, y) {
        if pw.len() == 2 {
            let one = pw.iter().find(|(k, _)| k.is_one());
            let other = pw.iter().find(|(k, _)| !k.is_one() && !k.is_zero());
            if let (Some((_, a)), Some((n, b))) = (one, other) {
                let m = sub(int_e(1), n.clone());
                let v = linear_first(&mul(neg(m.clone()), a.clone()), &mul(m.clone(), b.clone()), x, &c1)?;
                return Ok(Sol::Explicit(vec![pow(v, pow(m, int_e(-1)))]));
            }
        }
    }
    // exact: M dx + N dy = 0
    let (n, d) = together(f);
    let (m_, n_) = (neg(n), d);
    if ratsimp(&sub(diff(&m_, y), diff(&n_, x))).is_zero() {
        let mut psi = int(&m_, x)?;
        let rest = ratsimp(&sub(n_, diff(&psi, y)));
        if !rest.has_sym(x) {
            psi = add(psi, int(&rest, y)?);
            let psi = sub(psi, c1);
            return Ok(match explicit(&psi, y) {
                Some(v) => Sol::Explicit(v.into_iter().map(tidy).collect()),
                None => Sol::Implicit(psi),
            });
        }
    }
    Err("no method for this first-order equation (tried linear, separable, Bernoulli, exact)".into())
}

/// Numeric value with every free symbol set to a positive sample.
fn sample(e: &Expr) -> Option<f64> {
    let syms = e.free_syms();
    let env: Vec<(&str, C)> = syms.iter().enumerate().map(|(i, s)| (s.as_str(), C::new(0.83 + 0.37 * i as f64, 0.0))).collect();
    evalc(e, &env).filter(|z| z.im.abs() < 1e-9).map(|z| z.re)
}

fn second_order(a: &Expr, b: &Expr, c: &Expr, g: &Expr, x: &str) -> Result<Expr, String> {
    let xs = sym(x);
    let disc = expand(&sub(pow(b.clone(), int_e(2)), mul_v(vec![int_e(4), a.clone(), c.clone()])));
    // symbolic parameters are taken as positive when choosing the root type
    let dv = sample(&disc).ok_or("coefficients must be real")?;
    let two_a = mul(int_e(2), a.clone());
    let alpha = ratsimp(&div(neg(b.clone()), two_a.clone()));
    let (y1, y2) = if disc.is_zero() || dv.abs() < 1e-14 {
        let e1 = pow(E, mul(alpha.clone(), xs.clone()));
        (e1.clone(), mul(xs.clone(), e1))
    } else if dv > 0.0 {
        let s = div(sqrt(disc.clone()), two_a.clone());
        let (r1, r2) = (expand(&sub(alpha.clone(), s.clone())), expand(&add(alpha.clone(), s)));
        (pow(E, mul(r1, xs.clone())), pow(E, mul(r2, xs.clone())))
    } else {
        let beta = algebra::simplify(&div(sqrt(expand(&neg(disc.clone()))), two_a));
        let env = pow(E, mul(alpha.clone(), xs.clone()));
        (mul(env.clone(), f1("cos", mul(beta.clone(), xs.clone()))), mul(env, f1("sin", mul(beta, xs.clone()))))
    };
    let mut y = add(mul(sym("C1"), y1.clone()), mul(sym("C2"), y2.clone()));
    if !g.is_zero() {
        // variation of parameters
        let w = algebra::simplify(&sub(mul(y1.clone(), diff(&y2, x)), mul(diff(&y1, x), y2.clone())));
        let gn = div(g.clone(), a.clone());
        let u1 = int(&tidy(div(mul(y2.clone(), gn.clone()), w.clone())), x)?;
        let u2 = int(&tidy(div(mul(y1.clone(), gn), w)), x)?;
        let part = sub(mul(y2, u2), mul(y1, u1));
        let lin = algebra::trig_linearize(&part);
        y = add(y, tidy(if lin.size() < part.size() { lin } else { part }));
    }
    Ok(y)
}

fn verify(f: &Expr, x: &str, y: &str, sol: &Expr) -> bool {
    let (d1, d2) = (diff(sol, x), diff(&diff(sol, x), x));
    let mut r = subst(f, &format!("{y}''"), &d2);
    r = subst(&r, &format!("{y}'"), &d1);
    r = subst(&r, y, sol);
    let syms = r.free_syms();
    let mut good = 0;
    for k in 0..5 {
        let env: Vec<(&str, C)> = syms.iter().enumerate().map(|(i, s)| (s.as_str(), C::new(0.61 + 0.29 * i as f64 + 0.17 * k as f64, 0.0))).collect();
        let Some(v) = evalc(&r, &env) else { continue };
        if !v.re.is_finite() || !v.im.is_finite() {
            continue;
        }
        let terms: Vec<&Expr> = if let Expr::Add(t) = &r { t.iter().collect() } else { vec![&r] };
        let scale: f64 = terms.iter().filter_map(|t| evalc(t, &env)).map(|t| t.norm()).filter(|t| t.is_finite()).sum();
        if v.norm() > 1e-6 * (1.0 + scale) {
            return false;
        }
        good += 1;
    }
    good > 0
}

/// Solve the ODE `eq` for y(x). `ics` is empty, [x0, y0] or [x0, y0, y'(x0)].
pub fn desolve(eq: &Expr, x: &str, y: &str, ics: &[Expr]) -> Result<Expr, String> {
    let f = match eq {
        Expr::Rel(Rel::Eq, l, r) => sub((**l).clone(), (**r).clone()),
        Expr::Rel(..) => return Err("desolve needs an equation".into()),
        e => e.clone(),
    };
    let (yp, ypp) = (format!("{y}'"), format!("{y}''"));
    let ys = sym(y);
    let free = |e: &Expr| !e.has_sym(y) && !e.has_sym(&yp) && !e.has_sym(&ypp);
    let sol = if f.has_sym(&ypp) {
        let c2 = coeffs_expr(&f, &ypp).filter(|c| c.len() == 2).ok_or("equation must be linear in the second derivative")?;
        let c1 = coeffs_expr(&c2[0], &yp).filter(|c| c.len() <= 2).ok_or("equation must be linear in y'")?;
        let c0 = coeffs_expr(&c1[0], y).filter(|c| c.len() <= 2).ok_or("equation must be linear in y")?;
        let (a, b, c) = (c2[1].clone(), c1.get(1).cloned().unwrap_or(int_e(0)), c0.get(1).cloned().unwrap_or(int_e(0)));
        let g = neg(c0[0].clone());
        if ![&a, &b, &c].iter().all(|k| free(k) && !k.has_sym(x)) || !free(&g) {
            return Err("second-order equations need constant coefficients: a*y'' + b*y' + c*y = g(x)".into());
        }
        Sol::Explicit(vec![second_order(&a, &b, &c, &g, x)?])
    } else if f.has_sym(&yp) {
        let c = coeffs_expr(&f, &yp).filter(|c| c.len() == 2).ok_or(format!("equation must be linear in {yp}"))?;
        first_order(&ratsimp(&div(neg(c[0].clone()), c[1].clone())), x, y)?
    } else {
        return Err(format!("no derivative found: write {y}' and {y}''"));
    };
    let consts = ["C1", "C2"];
    let apply_ics = |e: &Expr, implicit: bool| -> Result<Expr, String> {
        if ics.is_empty() {
            return Ok(e.clone());
        }
        let mut eqs = Vec::new();
        if implicit {
            eqs.push(subst(&subst(e, x, &ics[0]), y, &ics[1]));
        } else {
            eqs.push(sub(subst(e, x, &ics[0]), ics[1].clone()));
            if let Some(v) = ics.get(2) {
                eqs.push(sub(subst(&diff(e, x), x, &ics[0]), v.clone()));
            }
        }
        let vars: Vec<String> = consts.iter().filter(|c| e.has_sym(c)).take(eqs.len()).map(|c| c.to_string()).collect();
        if vars.len() != eqs.len() {
            return Err("number of initial conditions does not match the order".into());
        }
        let (sols, _) = solve::solve_system(&eqs, &vars, false)?;
        let row = sols.first().ok_or("initial conditions cannot be met")?;
        let mut out = e.clone();
        for (v, val) in vars.iter().zip(row) {
            out = subst(&out, v, val);
        }
        Ok(tidy(out))
    };
    match sol {
        Sol::Explicit(v) => {
            let mut out = Vec::new();
            for s in v {
                if !verify(&f, x, y, &s) {
                    continue;
                }
                match apply_ics(&s, false) {
                    Ok(e) => out.push(Expr::Rel(Rel::Eq, Box::new(ys.clone()), Box::new(e))),
                    Err(e) if ics.is_empty() => return Err(e),
                    Err(_) => {}
                }
            }
            if out.is_empty() {
                return Err(if ics.is_empty() { "solution failed verification".into() } else { "initial conditions cannot be met".into() });
            }
            Ok(func("or", out))
        }
        Sol::Implicit(psi) => Ok(Expr::Rel(Rel::Eq, Box::new(apply_ics(&psi, true)?), Box::new(int_e(0)))),
    }
}
