//! Equation solving: exact polynomial roots, isolation, substitution, numeric fallback, inequalities, systems.
use crate::algebra::{ratsimp, together};
use crate::expr::*;
use crate::matrix;
use crate::numeric::{self, evalc};
use crate::poly::{self, Poly, coeffs_expr};
use crate::simplify::*;
use num_bigint::BigInt;
use num_complex::Complex64 as C;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

pub struct Ctx {
    pub complex: bool,
    pub approx: bool,
}

fn tidy(e: Expr) -> Expr {
    if e.is_constant() && e.size() < 60 { expand(&e) } else { e }
}

fn rn(r: &Rat) -> Expr {
    Expr::Num(r.clone())
}

/// All n-th roots of a rational number, exact.
fn nth_roots(r: &Rat, n: usize) -> Vec<Expr> {
    let mag = pow(Expr::Num(r.abs()), rat(1, n as i64));
    (0..n)
        .map(|k| {
            let num = if r.is_negative() { 2 * k as i64 + 1 } else { 2 * k as i64 };
            let th = mul(rat(num, n as i64), PI);
            tidy(mul(mag.clone(), add(f1("cos", th.clone()), mul(I, f1("sin", th)))))
        })
        .collect()
}

fn quadratic(a: &Expr, b: &Expr, c: &Expr) -> Vec<Expr> {
    let disc = expand(&sub(pow(b.clone(), int(2)), mul_v(vec![int(4), a.clone(), c.clone()])));
    let two_a = mul(int(2), a.clone());
    if disc.is_zero() {
        return vec![tidy(div(neg(b.clone()), two_a))];
    }
    let s = sqrt(disc);
    vec![tidy(div(sub(neg(b.clone()), s.clone()), two_a.clone())), tidy(div(add(neg(b.clone()), s), two_a))]
}

fn cubic(f: &Poly) -> Vec<Expr> {
    let (a, b, c, d) = (f.coef(3), f.coef(2), f.coef(1), f.coef(0));
    let r = |n: i64| Rat::from_integer(BigInt::from(n));
    let p = (r(3) * &a * &c - &b * &b) / (r(3) * &a * &a);
    let q = (r(2) * &b * &b * &b - r(9) * &a * &b * &c + r(27) * &a * &a * &d) / (r(27) * &a * &a * &a);
    let shift = rn(&(-&b / (r(3) * &a)));
    let ts: Vec<Expr> = if p.is_zero() {
        nth_roots(&-q, 3)
    } else {
        let disc = &q * &q / r(4) + &p * &p * &p / r(27);
        if disc.is_positive() {
            let sq = sqrt(rn(&disc));
            let half = rn(&(-&q / r(2)));
            let real_cbrt = |v: Expr| {
                if numeric::real_const(&v).unwrap_or(0.0) >= 0.0 { pow(v, rat(1, 3)) } else { neg(pow(expand(&neg(v)), rat(1, 3))) }
            };
            let u = real_cbrt(add(half.clone(), sq.clone()));
            let v = real_cbrt(sub(half, sq));
            let s = add(u.clone(), v.clone());
            let w = mul_v(vec![rat(1, 2), sqrt(int(3)), I, sub(u, v)]);
            vec![s.clone(), add(mul(rat(-1, 2), s.clone()), w.clone()), sub(mul(rat(-1, 2), s), w)]
        } else {
            let m = mul(int(2), sqrt(rn(&(-&p / r(3)))));
            let arg = mul(rn(&(r(3) * &q / (r(2) * &p))), sqrt(rn(&(-r(3) / &p))));
            let th = mul(rat(1, 3), f1("acos", arg));
            (0..3).map(|k| mul(m.clone(), f1("cos", sub(th.clone(), mul(rat(2 * k, 3), PI))))).collect()
        }
    };
    ts.into_iter().map(|t| add(t, shift.clone())).collect()
}

fn quartic(f: &Poly, cx: &mut Ctx) -> Option<Vec<Expr>> {
    let r = |n: i64| Rat::from_integer(BigInt::from(n));
    let a4 = f.coef(4);
    let (b, c, d, e) = (f.coef(3) / &a4, f.coef(2) / &a4, f.coef(1) / &a4, f.coef(0) / &a4);
    let p = &c - r(3) * &b * &b / r(8);
    let q = &d - &b * &c / r(2) + &b * &b * &b / r(8);
    let rr = &e - &b * &d / r(4) + &b * &b * &c / r(16) - r(3) * &b * &b * &b * &b / r(256);
    let shift = rn(&(-&b / r(4)));
    let mut ys = Vec::new();
    if q.is_zero() {
        for z in quadratic(&int(1), &rn(&p), &rn(&rr)) {
            let s = sqrt(z);
            ys.push(s.clone());
            ys.push(neg(s));
        }
    } else {
        let res = Poly::new(vec![-(&q * &q), r(2) * &p * &p - r(8) * &rr, r(8) * &p, r(8)]);
        let cands = solve_poly_exact(&res, cx);
        let m = cands.into_iter().find(|m| numeric::eval_const(m).is_some_and(|z| z.im.abs() < 1e-9 && z.re > 1e-12))?;
        let s = sqrt(mul(int(2), m.clone()));
        for s1 in [1, -1] {
            let inner = neg(add_v(vec![rn(&(r(2) * &p)), mul(int(2), m.clone()), mul_v(vec![int(s1), rn(&(r(2) * &q)), pow(s.clone(), int(-1))])]));
            let root = sqrt(tidy(inner));
            for s2 in [1, -1] {
                ys.push(mul(rat(1, 2), add(mul(int(s1), s.clone()), mul(int(s2), root.clone()))));
            }
        }
    }
    Some(ys.into_iter().map(|y| tidy(add(y, shift.clone()))).collect())
}

fn check_roots(f: &Poly, roots: &[Expr]) -> bool {
    let scale = f.0.iter().map(|c| rat_f64(c).abs()).fold(0.0, f64::max);
    roots.len() == f.deg()
        && roots.iter().all(|r| match numeric::eval_const(r) {
            Some(z) => f.evalc(z).norm() <= 1e-7 * scale * (1.0 + z.norm().powi(f.deg() as i32)),
            None => false,
        })
}

/// All complex roots of a rational polynomial; exact where a closed form exists.
pub fn solve_poly_exact(p: &Poly, cx: &mut Ctx) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    for (f, _) in poly::factor(p).1 {
        let n = f.deg();
        let exact: Option<Vec<Expr>> = match n {
            0 => Some(vec![]),
            1 => Some(vec![rn(&(-f.coef(0) / f.coef(1)))]),
            2 => Some(quadratic(&rn(&f.coef(2)), &rn(&f.coef(1)), &rn(&f.coef(0)))),
            _ if f.0[1..n].iter().all(|c| c.is_zero()) => Some(nth_roots(&(-f.coef(0) / f.coef(n)), n)),
            3 => Some(cubic(&f)),
            4 => quartic(&f, cx),
            _ => {
                // polynomial in x^k
                let k = (1..=n).filter(|&i| !f.coef(i).is_zero()).fold(0, num_integer::gcd);
                if k > 1 {
                    let reduced = Poly::new((0..=n / k).map(|i| f.coef(i * k)).collect());
                    let mut v = Vec::new();
                    for y in solve_poly_exact(&reduced, cx) {
                        match y.as_rat() {
                            Some(r) => v.extend(nth_roots(r, k)),
                            None => {
                                let base = pow(y, rat(1, k as i64));
                                for j in 0..k {
                                    let th = mul(rat(2 * j as i64, k as i64), PI);
                                    v.push(tidy(mul(base.clone(), add(f1("cos", th.clone()), mul(I, f1("sin", th))))));
                                }
                            }
                        }
                    }
                    Some(v)
                } else {
                    None
                }
            }
        };
        match exact.filter(|r| check_roots(&f, r)) {
            Some(r) => out.extend(r),
            None => {
                cx.approx = true;
                out.extend(poly::roots_numeric(&f).into_iter().map(from_c));
            }
        }
    }
    out
}

fn solve_poly_coeffs(mut c: Vec<Expr>, cx: &mut Ctx) -> Option<Vec<Expr>> {
    while c.last().is_some_and(|k| k.is_zero()) {
        c.pop();
    }
    if c.len() <= 1 {
        return Some(vec![]);
    }
    if let Some(rc) = c.iter().map(|k| k.as_rat().cloned()).collect::<Option<Vec<Rat>>>() {
        return Some(solve_poly_exact(&Poly::new(rc), cx));
    }
    // float coefficients: numeric roots
    if c.iter().all(|k| k.is_constant()) && c.iter().any(|k| k.has_float()) && c.len() > 3 {
        let vals: Option<Vec<C>> = c.iter().map(numeric::eval_const).collect();
        if let Some(v) = vals.filter(|v| v.iter().all(|z| z.im == 0.0)) {
            let p = Poly::new(v.iter().map(|z| Rat::from_float(z.re)).collect::<Option<Vec<Rat>>>()?);
            let g = p.gcd(&p.deriv());
            cx.approx = true;
            return Some(poly::roots_numeric(&p.divrem(&g).0).into_iter().map(from_c).collect());
        }
    }
    let n = c.len() - 1;
    match n {
        1 if c[1].is_number() => Some(vec![expand(&div(neg(c[0].clone()), c[1].clone()))]),
        1 => Some(vec![ratsimp(&div(neg(c[0].clone()), c[1].clone()))]),
        2 => Some(quadratic(&c[2], &c[1], &c[0])),
        _ if c[1..n].iter().all(|k| k.is_zero()) => {
            let r = ratsimp(&div(neg(c[0].clone()), c[n].clone()));
            let root = pow(r, rat(1, n as i64));
            Some(if n % 2 == 0 { vec![neg(root.clone()), root] } else { vec![root] })
        }
        _ => None,
    }
}

fn periodic(k: i64) -> Expr {
    mul_v(vec![int(k), sym("@n"), PI])
}

/// Solutions u of name(u) = r.
fn invert(name: &str, r: &Expr) -> Option<Vec<Expr>> {
    let a = |n: &str| f1(n, r.clone());
    Some(match name {
        "sin" => vec![add(a("asin"), periodic(2)), add_v(vec![PI, neg(a("asin")), periodic(2)])],
        "cos" => vec![add(a("acos"), periodic(2)), add(neg(a("acos")), periodic(2))],
        "tan" => vec![add(a("atan"), periodic(1))],
        "asin" => vec![a("sin")],
        "acos" => vec![a("cos")],
        "atan" => vec![a("tan")],
        "sinh" => vec![a("asinh")],
        "cosh" => vec![a("acosh"), neg(a("acosh"))],
        "tanh" => vec![a("atanh")],
        "asinh" => vec![a("sinh")],
        "acosh" => vec![a("cosh")],
        "atanh" => vec![a("tanh")],
        "ln" => vec![pow(E, r.clone())],
        "abs" => vec![r.clone(), neg(r.clone())],
        _ => return None,
    })
}

fn union(into: &mut Vec<Expr>, more: Vec<Expr>) {
    for m in more {
        if !into.contains(&m) {
            into.push(m);
        }
    }
}

/// Solve lhs = rhs where only `lhs` contains x, by peeling the outermost operation.
fn isolate(lhs: &Expr, rhs: Expr, x: &str, cx: &mut Ctx, depth: usize) -> Option<Vec<Expr>> {
    match lhs {
        Expr::Sym(s) if s == x => Some(vec![rhs]),
        Expr::Mul(v) => {
            let (w, k): (Vec<Expr>, Vec<Expr>) = v.iter().cloned().partition(|f| f.has_sym(x));
            if w.len() != 1 {
                return None;
            }
            isolate(&w[0], div(rhs, mul_v(k)), x, cx, depth)
        }
        Expr::Pow(b, e) => {
            let (bx, ex) = (b.has_sym(x), e.has_sym(x));
            let mut out = Vec::new();
            if bx && !ex {
                let mut targets = Vec::new();
                if let Some(n) = e.as_i64() {
                    if n < 0 {
                        return isolate(&pow((**b).clone(), int(-n)), pow(rhs, int(-1)), x, cx, depth);
                    }
                    let root = pow(rhs.clone(), rat(1, n));
                    if cx.complex {
                        for j in 0..n.min(64) {
                            let th = mul(rat(2 * j, n), PI);
                            targets.push(tidy(mul(root.clone(), add(f1("cos", th.clone()), mul(I, f1("sin", th))))));
                        }
                    } else {
                        targets.push(root.clone());
                        if n % 2 == 0 {
                            targets.push(neg(root));
                        }
                    }
                } else {
                    let t = pow(rhs, pow((**e).clone(), int(-1)));
                    // even numerator: u^(p/q) is even in u under the real-root convention
                    if e.as_rat().is_some_and(|r| r.numer().is_even()) {
                        targets.push(neg(t.clone()));
                    }
                    targets.push(t);
                }
                for t in targets {
                    union(&mut out, roots(&sub((**b).clone(), t), x, cx, depth + 1)?);
                }
            } else if ex && !bx {
                let t = div(f1("ln", rhs), f1("ln", (**b).clone()));
                union(&mut out, roots(&sub((**e).clone(), t), x, cx, depth + 1)?);
            } else {
                return None;
            }
            Some(out)
        }
        Expr::Func(name, a) if a.len() == 1 => {
            let mut out = Vec::new();
            for t in invert(name, &rhs)? {
                union(&mut out, roots(&sub(a[0].clone(), t), x, cx, depth + 1)?);
            }
            Some(out)
        }
        _ => None,
    }
}

fn generators(e: &Expr, x: &str, out: &mut Vec<Expr>) {
    let push = |g: Expr, out: &mut Vec<Expr>| {
        if !out.contains(&g) {
            out.push(g);
        }
    };
    match e {
        Expr::Func(_, _) if e.has_sym(x) => push(e.clone(), out),
        Expr::Pow(b, ex) if ex.has_sym(x) && !b.has_sym(x) => {
            // unit exponential: strip the numeric coefficient of the exponent
            let unit = match &**ex {
                Expr::Mul(v) if v[0].is_number() => mul_v(v[1..].to_vec()),
                other => other.clone(),
            };
            push(pow((**b).clone(), unit), out);
            push(e.clone(), out);
        }
        Expr::Pow(b, ex) if b.has_sym(x) && ex.as_rat().is_some_and(|r| !r.is_integer()) => push(pow((**b).clone(), Expr::Num(Rat::new(BigInt::one(), ex.as_rat().unwrap().denom().clone()))), out),
        _ => {}
    }
    for c in e.children() {
        generators(c, x, out);
    }
}

fn residual_ok(f: &Expr, x: &str, s: &Expr) -> bool {
    if !s.is_constant() {
        return true;
    }
    let Some(z) = numeric::eval_const(s) else { return true };
    if !z.re.is_finite() {
        return false;
    }
    let env = [(x, z)];
    let Some(v) = evalc(f, &env) else { return true };
    let terms: Vec<&Expr> = if let Expr::Add(t) = f { t.iter().collect() } else { vec![f] };
    let scale: f64 = terms.iter().filter_map(|t| evalc(t, &env)).map(|t| t.norm()).filter(|t| t.is_finite()).sum();
    v.norm() <= 1e-7 * (1.0 + scale)
}

fn a_sin_b_cos(n: &Expr, x: &str, cx: &mut Ctx, depth: usize) -> Option<Vec<Expr>> {
    // a*sin(u) + b*cos(u) + c = 0  ->  R*sin(u + phi) = -c
    let Expr::Add(terms) = n else { return None };
    let (mut a, mut b, mut c, mut u) = (int(0), int(0), Vec::new(), None::<Expr>);
    for t in terms {
        if !t.has_sym(x) {
            c.push(t.clone());
            continue;
        }
        let (k, fs): (Vec<Expr>, Vec<Expr>) = (if let Expr::Mul(v) = t { v.clone() } else { vec![t.clone()] }).into_iter().partition(|f| !f.has_sym(x));
        let [Expr::Func(name, arg)] = &fs[..] else { return None };
        if u.as_ref().is_some_and(|u| *u != arg[0]) {
            return None;
        }
        u = Some(arg[0].clone());
        match name.as_str() {
            "sin" => a = add(a, mul_v(k)),
            "cos" => b = add(b, mul_v(k)),
            _ => return None,
        }
    }
    if a.is_zero() || b.is_zero() {
        return None;
    }
    let r = sqrt(add(pow(a.clone(), int(2)), pow(b.clone(), int(2))));
    let (av, bv) = (numeric::real_const(&a)?, numeric::real_const(&b)?);
    let mut phi = f1("atan", div(b, a));
    if av < 0.0 {
        phi = if bv >= 0.0 { add(phi, PI) } else { sub(phi, PI) };
    }
    let target = tidy(div(neg(add_v(c)), r));
    let mut out = Vec::new();
    for t in invert("sin", &target)? {
        union(&mut out, roots(&add_v(vec![u.clone()?, phi.clone(), neg(t)]), x, cx, depth + 1)?);
    }
    Some(out)
}

/// Solutions of f = 0 for x, or None when no symbolic method applies.
fn roots(f: &Expr, x: &str, cx: &mut Ctx, depth: usize) -> Option<Vec<Expr>> {
    if depth > 10 {
        return None;
    }
    let (n, _) = together(f);
    if !n.has_sym(x) {
        return Some(vec![]);
    }
    if let Expr::Mul(v) = &n {
        let mut out = Vec::new();
        let mut ok = true;
        for fac in v.iter().filter(|f| f.has_sym(x)) {
            match roots(fac, x, cx, depth + 1) {
                Some(r) => union(&mut out, r),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            return Some(out);
        }
    }
    if let Expr::Pow(b, k) = &n {
        if !k.has_sym(x) && numeric::real_const(k).is_some_and(|v| v > 0.0) {
            return roots(b, x, cx, depth + 1);
        }
    }
    if let Some(c) = coeffs_expr(&n, x) {
        if let Some(r) = solve_poly_coeffs(c, cx) {
            return Some(r);
        }
    }
    // single occurrence: peel
    let (with_x, without): (Vec<Expr>, Vec<Expr>) = (if let Expr::Add(v) = &n { v.clone() } else { vec![n.clone()] }).into_iter().partition(|t| t.has_sym(x));
    if with_x.len() == 1 {
        if let Some(r) = isolate(&with_x[0], neg(add_v(without.clone())), x, cx, depth) {
            return Some(r);
        }
    }
    if let Some(r) = a_sin_b_cos(&n, x, cx, depth) {
        return Some(r);
    }
    // several logarithms: merge into one
    if with_x.iter().filter(|t| t.has(&|e| matches!(e, Expr::Func(f, _) if f == "ln"))).count() >= 2 {
        let merged = crate::algebra::combine_logs(&n);
        if merged != n {
            if let Some(r) = roots(&merged, x, cx, depth + 1) {
                return Some(r.into_iter().filter(|s| residual_ok(&n, x, s)).collect());
            }
        }
    }
    // f(u) = f(v) for one-to-one f
    if with_x.len() == 2 && without.is_empty() {
        if let (Expr::Func(n1, a1), Expr::Mul(m)) = (&with_x[0], &with_x[1]) {
            if let [c, Expr::Func(n2, a2)] = &m[..] {
                if *c == int(-1) && n1 == n2 && ["ln", "asin", "atan", "sinh", "asinh", "tanh"].contains(&n1.as_str()) {
                    return roots(&sub(a1[0].clone(), a2[0].clone()), x, cx, depth + 1);
                }
            }
        }
    }
    // polynomial in a generator
    let mut gens = Vec::new();
    generators(&n, x, &mut gens);
    let tname = format!("@t{depth}");
    let t = sym(&tname);
    for u in gens.iter().filter(|u| **u != n) {
        let r = replace_gen(&n, u, &t);
        if r.has_sym(x) || !r.has_sym(&tname) {
            continue;
        }
        let Some(ts) = roots(&r, &tname, cx, depth + 1) else { continue };
        let mut out = Vec::new();
        let mut ok = true;
        for tv in ts {
            match isolate(u, tv, x, cx, depth + 1) {
                Some(r) => union(&mut out, r),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            return Some(out);
        }
    }
    // one square-root term: isolate it and square
    let rad: Vec<&Expr> = with_x.iter().filter(|t| t.has(&|e| matches!(e, Expr::Pow(b, k) if b.has_sym(x) && k.as_rat().is_some_and(|r| *r.denom() == BigInt::from(2))))).collect();
    if rad.len() == 1 && depth < 4 {
        let t = rad[0].clone();
        let rest = sub(n.clone(), t.clone());
        let sq = expand(&sub(pow(t, int(2)), pow(rest, int(2))));
        if let Some(r) = roots(&sq, x, cx, depth + 1) {
            return Some(r.into_iter().filter(|s| residual_ok(&n, x, s)).collect());
        }
    }
    None
}

fn numeric_roots(f: &Expr, x: &str) -> Option<Vec<Expr>> {
    let syms = f.free_syms();
    if syms.len() != 1 || !syms.contains(x) {
        return None;
    }
    let g = |v: f64| numeric::evalf(f, x, v);
    // ponytail: numeric search window is [-100, 100]; use nsolve(eq, x, guess) outside it
    let mut r = numeric::find_roots(&g, -100.0, 100.0, 20_000);
    r.truncate(64);
    Some(r.into_iter().map(Expr::Float).collect())
}

pub enum Solution {
    All,
    Set(Vec<Expr>),
}

fn sort_solutions(v: &mut [Expr]) {
    v.sort_by(|a, b| {
        let key = |e: &Expr| numeric::eval_const(e).map(|z| (z.im.abs() > 1e-12, z.re, z.im)).unwrap_or((true, f64::INFINITY, 0.0));
        let (ka, kb) = (key(a), key(b));
        ka.0.cmp(&kb.0).then(ka.1.total_cmp(&kb.1)).then(ka.2.total_cmp(&kb.2))
    });
}

/// Solve f = 0 for x.
pub fn solve(f: &Expr, x: &str, complex: bool) -> Result<(Solution, bool), String> {
    let f = simp(f);
    if !f.has_sym(x) {
        return Ok((if crate::algebra::is_zero_expr(&f) { Solution::All } else { Solution::Set(vec![]) }, false));
    }
    if ratsimp(&f).is_zero() {
        return Ok((Solution::All, false));
    }
    let mut cx = Ctx { complex, approx: false };
    let mut sols = match roots(&f, x, &mut cx, 0) {
        Some(s) => s,
        None => {
            cx.approx = true;
            numeric_roots(&f, x).ok_or_else(|| format!("cannot solve for {x}"))?
        }
    };
    sols.retain(|s| residual_ok(&f, x, s));
    if !complex {
        sols.retain(|s| !s.is_constant() || numeric::eval_const(s).is_none_or(|z| z.im.abs() <= 1e-9 * (1.0 + z.re.abs())));
    }
    // dedupe numerically equal constants
    let mut out: Vec<Expr> = Vec::new();
    for s in sols {
        let dup = out.iter().any(|o| {
            *o == s
                || match (numeric::eval_const(o), numeric::eval_const(&s)) {
                    (Some(a), Some(b)) if o.is_constant() && s.is_constant() => (a - b).norm() <= 1e-10 * (1.0 + a.norm()),
                    _ => false,
                }
        });
        if !dup {
            out.push(s);
        }
    }
    sort_solutions(&mut out);
    Ok((Solution::Set(out), cx.approx))
}

pub fn solution_expr(x: &str, s: Solution) -> Expr {
    match s {
        Solution::All => Expr::Bool(true),
        Solution::Set(v) => func("or", v.into_iter().map(|r| Expr::Rel(Rel::Eq, Box::new(sym(x)), Box::new(r))).collect()),
    }
}

/// Solve f REL 0 over the reals.
pub fn solve_inequality(f: &Expr, rel: Rel, x: &str) -> Result<Expr, String> {
    let (n, d) = together(f);
    let mut pts: Vec<Expr> = Vec::new();
    for part in [&n, &d] {
        if part.has_sym(x) {
            match solve(part, x, false)?.0 {
                Solution::All => return Err("inequality is degenerate".into()),
                Solution::Set(v) => union(&mut pts, v),
            }
        }
    }
    let mut cp: Vec<(f64, Expr)> = Vec::new();
    for p in pts {
        let v = if p.is_constant() { numeric::real_const(&p) } else { None }.ok_or("inequality needs numeric critical points")?;
        if !cp.iter().any(|(w, _)| (w - v).abs() < 1e-12) {
            cp.push((v, p));
        }
    }
    cp.sort_by(|a, b| a.0.total_cmp(&b.0));
    let holds = |v: f64| {
        let y = numeric::evalf(f, x, v);
        !y.is_nan() && rel.holds(if y.abs() < 1e-12 { 0.0 } else { y }, 0.0)
    };
    // alternating sequence: interval, point, interval, ...
    let k = cp.len();
    let mut seq = Vec::with_capacity(2 * k + 1);
    for i in 0..=k {
        let probe = if k == 0 {
            0.0
        } else if i == 0 {
            cp[0].0 - 1.0
        } else if i == k {
            cp[k - 1].0 + 1.0
        } else {
            0.5 * (cp[i - 1].0 + cp[i].0)
        };
        seq.push(holds(probe));
        if i < k {
            // exact check at the critical point: numerator root satisfies non-strict relations only
            let at = numeric::evalf(&d, x, cp[i].0).abs() > 1e-12 && numeric::evalf(&n, x, cp[i].0).abs() < 1e-9 && matches!(rel, Rel::Le | Rel::Ge | Rel::Eq);
            let other = numeric::evalf(&d, x, cp[i].0).abs() > 1e-12 && numeric::evalf(&n, x, cp[i].0).abs() >= 1e-9 && holds(cp[i].0);
            seq.push(at || other);
        }
    }
    let xs = sym(x);
    let mut pieces = Vec::new();
    let mut i = 0;
    while i < seq.len() {
        if !seq[i] {
            i += 1;
            continue;
        }
        let s = i;
        while i + 1 < seq.len() && seq[i + 1] {
            i += 1;
        }
        let e = i;
        i += 1;
        if s == 0 && e == seq.len() - 1 {
            return Ok(Expr::Bool(true));
        }
        if s == e && s % 2 == 1 {
            pieces.push(Expr::Rel(Rel::Eq, Box::new(xs.clone()), Box::new(cp[s / 2].1.clone())));
            continue;
        }
        let lower = if s == 0 { None } else if s % 2 == 1 { Some((Rel::Ge, cp[s / 2].1.clone())) } else { Some((Rel::Gt, cp[s / 2 - 1].1.clone())) };
        let upper = if e == seq.len() - 1 { None } else if e % 2 == 1 { Some((Rel::Le, cp[e / 2].1.clone())) } else { Some((Rel::Lt, cp[e / 2].1.clone())) };
        let mut conds = Vec::new();
        if let (Some((r, v)), true) = (&lower, upper.is_some()) {
            // both bounds: a < x, x < b prints as a < x < b
            conds.push(Expr::Rel(r.flip(), Box::new(v.clone()), Box::new(xs.clone())));
        } else if let Some((r, v)) = lower {
            conds.push(Expr::Rel(r, Box::new(xs.clone()), Box::new(v)));
        }
        if let Some((r, v)) = upper {
            conds.push(Expr::Rel(r, Box::new(xs.clone()), Box::new(v)));
        }
        pieces.push(func("and", conds));
    }
    Ok(func("or", pieces))
}

fn linear_system(eqs: &[Expr], vars: &[String]) -> Option<Result<Vec<Vec<Expr>>, String>> {
    let mut aug: matrix::Mat = Vec::new();
    for eq in eqs {
        let mut row = Vec::new();
        let mut rest = expand(eq);
        for v in vars {
            let c = coeffs_expr(&rest, v)?;
            if c.len() > 2 {
                return None;
            }
            let k = c.get(1).cloned().unwrap_or(int(0));
            if vars.iter().any(|w| k.has_sym(w)) {
                return None;
            }
            row.push(k);
            rest = c[0].clone();
        }
        if vars.iter().any(|w| rest.has_sym(w)) {
            return None;
        }
        row.push(neg(rest));
        aug.push(row);
    }
    let (m, piv, _) = matrix::rref(&aug);
    let n = vars.len();
    if piv.contains(&n) {
        return Some(Ok(vec![]));
    }
    if piv.len() < n {
        return Some(Err("system has infinitely many solutions".into()));
    }
    Some(Ok(vec![(0..n).map(|i| m[i][n].clone()).collect()]))
}

fn eliminate(eqs: &[Expr], vars: &[String], cx: &mut Ctx, depth: usize) -> Option<Vec<Vec<(String, Expr)>>> {
    if vars.is_empty() {
        return Some(if eqs.iter().all(crate::algebra::is_zero_expr) { vec![vec![]] } else { vec![] });
    }
    if eqs.is_empty() || depth > 6 {
        return None;
    }
    // prefer an equation that is linear in some variable
    let mut order: Vec<(usize, usize, usize)> = Vec::new();
    for (i, eq) in eqs.iter().enumerate() {
        for (j, v) in vars.iter().enumerate() {
            if eq.has_sym(v) {
                let deg = coeffs_expr(&together(eq).0, v).map_or(99, |c| c.len() - 1);
                order.push((deg, i, j));
            }
        }
    }
    order.sort();
    for (_, i, j) in order.into_iter().take(4) {
        let v = &vars[j];
        let mut sub_cx = Ctx { complex: cx.complex, approx: false };
        let Some(sols) = roots(&eqs[i], v, &mut sub_cx, 0) else { continue };
        if sub_cx.approx {
            continue;
        }
        let rest_vars: Vec<String> = vars.iter().filter(|w| *w != v).cloned().collect();
        let mut out = Vec::new();
        let mut ok = true;
        for s in sols {
            let rest_eqs: Vec<Expr> = eqs.iter().enumerate().filter(|(k, _)| *k != i).map(|(_, e)| together(&subst(e, v, &s)).0).map(|e| if e.size() < 400 { expand(&e) } else { e }).filter(|e| !e.is_zero()).collect();
            match eliminate(&rest_eqs, &rest_vars, cx, depth + 1) {
                Some(sub_sols) => {
                    for mut ss in sub_sols {
                        let mut val = s.clone();
                        for (w, wv) in &ss {
                            val = subst(&val, w, wv);
                        }
                        ss.push((v.clone(), tidy(val)));
                        out.push(ss);
                    }
                }
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            return Some(out);
        }
    }
    None
}

fn newton_system(eqs: &[Expr], vars: &[String]) -> Vec<Vec<Expr>> {
    let n = vars.len();
    let eval = |p: &[f64]| -> Vec<f64> {
        let env: Vec<(&str, C)> = vars.iter().zip(p).map(|(v, x)| (v.as_str(), C::new(*x, 0.0))).collect();
        eqs.iter().map(|e| evalc(e, &env).filter(|z| z.im.abs() < 1e-9).map_or(f64::NAN, |z| z.re)).collect()
    };
    let mut found: Vec<Vec<f64>> = Vec::new();
    let seeds = [0.5, -0.5, 1.5, -1.5, 3.0, -3.0, 7.0, -7.0];
    let total = seeds.len().pow(n.min(3) as u32);
    for s in 0..total {
        let mut p: Vec<f64> = (0..n).map(|i| seeds[(s / seeds.len().pow(i.min(2) as u32)) % seeds.len()] * (1.0 + 0.1 * i as f64)).collect();
        let mut ok = false;
        for _ in 0..60 {
            let f = eval(&p);
            if f.iter().any(|v| v.is_nan()) {
                break;
            }
            if f.iter().map(|v| v.abs()).fold(0.0, f64::max) < 1e-12 {
                ok = true;
                break;
            }
            let mut jac = vec![vec![0.0; n + 1]; n];
            for j in 0..n {
                let h = 1e-6 * (1.0 + p[j].abs());
                let mut q = p.clone();
                q[j] += h;
                let fp = eval(&q);
                q[j] -= 2.0 * h;
                let fm = eval(&q);
                for i in 0..n {
                    jac[i][j] = (fp[i] - fm[i]) / (2.0 * h);
                }
            }
            for i in 0..n {
                jac[i][n] = -f[i];
            }
            // Gaussian elimination with partial pivoting
            let mut singular = false;
            for c in 0..n {
                let piv = (c..n).max_by(|a, b| jac[*a][c].abs().total_cmp(&jac[*b][c].abs())).unwrap();
                if jac[piv][c].abs() < 1e-14 {
                    singular = true;
                    break;
                }
                jac.swap(c, piv);
                for r in 0..n {
                    if r != c {
                        let k = jac[r][c] / jac[c][c];
                        for cc in c..=n {
                            jac[r][cc] -= k * jac[c][cc];
                        }
                    }
                }
            }
            if singular {
                break;
            }
            for i in 0..n {
                p[i] += jac[i][n] / jac[i][i];
            }
        }
        if ok && !found.iter().any(|q| q.iter().zip(&p).all(|(a, b)| (a - b).abs() < 1e-8 * (1.0 + a.abs()))) {
            found.push(p);
        }
    }
    found.sort_by(|a, b| a[0].total_cmp(&b[0]));
    found.into_iter().map(|p| p.into_iter().map(Expr::Float).collect()).collect()
}

/// Solve a system of equations (each given as f = 0). Returns one value vector per solution.
pub fn solve_system(eqs: &[Expr], vars: &[String], complex: bool) -> Result<(Vec<Vec<Expr>>, bool), String> {
    if let Some(r) = linear_system(eqs, vars) {
        return r.map(|v| (v, false));
    }
    let mut cx = Ctx { complex, approx: false };
    let nums: Vec<Expr> = eqs.iter().map(|e| together(e).0).collect();
    if let Some(sols) = eliminate(&nums, vars, &mut cx, 0) {
        let mut out: Vec<Vec<Expr>> = Vec::new();
        for s in sols {
            let row: Vec<Expr> = vars.iter().map(|v| s.iter().find(|(w, _)| w == v).map(|(_, e)| e.clone()).unwrap_or_else(|| sym(v))).collect();
            let vals: Option<Vec<C>> = row.iter().map(|e| if e.is_constant() { numeric::eval_const(e) } else { None }).collect();
            if let Some(vals) = &vals {
                if !complex && vals.iter().any(|z| z.im.abs() > 1e-9) {
                    continue;
                }
                let env: Vec<(&str, C)> = vars.iter().map(|v| v.as_str()).zip(vals.iter().cloned()).collect();
                if eqs.iter().any(|e| evalc(e, &env).is_none_or(|z| !(z.norm() < 1e-6))) {
                    continue;
                }
            }
            if !out.contains(&row) {
                out.push(row);
            }
        }
        return Ok((out, cx.approx));
    }
    if eqs.len() == vars.len() {
        let all: std::collections::BTreeSet<String> = eqs.iter().flat_map(|e| e.free_syms()).collect();
        if all.iter().all(|s| vars.contains(s)) {
            return Ok((newton_system(eqs, vars), true));
        }
    }
    Err("cannot solve this system".into())
}

