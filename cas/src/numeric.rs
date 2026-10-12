//! Floating-point evaluation, special functions, root finding and quadrature.
use crate::expr::{Const, Expr, Rat, rat_f64};
use num_complex::Complex64 as C;
use num_integer::Integer;
use num_traits::{Signed, ToPrimitive};
use std::f64::consts::PI;

pub fn c(re: f64) -> C {
    C::new(re, 0.0)
}

/// Rational power of a real base. Odd denominators give the real root of a negative base.
pub fn rat_pow(b: C, e: &Rat) -> C {
    if b.im == 0.0 && b.re < 0.0 && !e.is_integer() && e.denom().is_odd() {
        let m = (-b.re).powf(rat_f64(e));
        return c(if e.numer().is_odd() { -m } else { m });
    }
    if e.is_integer() {
        if let Some(n) = e.numer().to_i32() {
            return if b.im == 0.0 { c(b.re.powi(n)) } else { b.powi(n) };
        }
    }
    if b.im == 0.0 && b.re >= 0.0 {
        return c(b.re.powf(rat_f64(e)));
    }
    b.powc(c(rat_f64(e)))
}

pub fn evalc(e: &Expr, env: &[(&str, C)]) -> Option<C> {
    Some(match e {
        Expr::Num(r) => c(rat_f64(r)),
        Expr::Float(f) => c(*f),
        Expr::Sym(s) => env.iter().find(|(n, _)| n == s)?.1,
        Expr::Const(k) => match k {
            Const::Pi => c(PI),
            Const::E => c(std::f64::consts::E),
            Const::I => C::new(0.0, 1.0),
            Const::Inf => c(f64::INFINITY),
            Const::Undef => c(f64::NAN),
        },
        Expr::Add(v) => {
            let mut s = c(0.0);
            for t in v {
                s += evalc(t, env)?;
            }
            s
        }
        Expr::Mul(v) => {
            let mut s = c(1.0);
            for t in v {
                s *= evalc(t, env)?;
            }
            s
        }
        Expr::Pow(b, x) => {
            let bv = evalc(b, env)?;
            if let Expr::Num(r) = &**x {
                return Some(rat_pow(bv, r));
            }
            let xv = evalc(x, env)?;
            if bv.im == 0.0 && xv.im == 0.0 && (bv.re > 0.0 || (bv.re == 0.0 && xv.re > 0.0)) {
                c(bv.re.powf(xv.re))
            } else if bv.im == 0.0 && bv.re == 0.0 {
                c(if xv.re == 0.0 { 1.0 } else { f64::NAN })
            } else {
                bv.powc(xv)
            }
        }
        Expr::Func(name, args) => {
            let mut a = Vec::with_capacity(args.len());
            for x in args {
                a.push(evalc(x, env)?);
            }
            call(name, &a)?
        }
        _ => return None,
    })
}

/// Real value of `e` at `var = x`; NaN when undefined or non-real.
pub fn evalf(e: &Expr, var: &str, x: f64) -> f64 {
    match evalc(e, &[(var, c(x))]) {
        Some(v) if v.im.abs() <= 1e-12 * (1.0 + v.re.abs()) => v.re,
        _ => f64::NAN,
    }
}

pub fn eval_const(e: &Expr) -> Option<C> {
    evalc(e, &[])
}

/// Real numeric value of a constant expression.
pub fn real_const(e: &Expr) -> Option<f64> {
    let v = eval_const(e)?;
    (v.im.abs() <= 1e-12 * (1.0 + v.re.abs()) && !v.re.is_nan()).then_some(v.re)
}

fn real_fn(a: C, dom: impl Fn(f64) -> bool, f: impl Fn(f64) -> f64, g: impl Fn(C) -> C) -> C {
    if a.im == 0.0 && dom(a.re) { c(f(a.re)) } else { g(a) }
}

pub fn call(name: &str, a: &[C]) -> Option<C> {
    let x = *a.first()?;
    let r = |i: usize| a.get(i).map(|v| v.re);
    Some(match (name, a.len()) {
        ("sin", 1) => real_fn(x, |_| true, f64::sin, |z| z.sin()),
        ("cos", 1) => real_fn(x, |_| true, f64::cos, |z| z.cos()),
        ("tan", 1) => real_fn(x, |_| true, f64::tan, |z| z.tan()),
        ("asin", 1) => real_fn(x, |v| v.abs() <= 1.0, f64::asin, |z| z.asin()),
        ("acos", 1) => real_fn(x, |v| v.abs() <= 1.0, f64::acos, |z| z.acos()),
        ("atan", 1) => real_fn(x, |_| true, f64::atan, |z| z.atan()),
        ("sinh", 1) => real_fn(x, |_| true, f64::sinh, |z| z.sinh()),
        ("cosh", 1) => real_fn(x, |_| true, f64::cosh, |z| z.cosh()),
        ("tanh", 1) => real_fn(x, |_| true, f64::tanh, |z| z.tanh()),
        ("asinh", 1) => real_fn(x, |_| true, f64::asinh, |z| z.asinh()),
        ("acosh", 1) => real_fn(x, |v| v >= 1.0, f64::acosh, |z| z.acosh()),
        ("atanh", 1) => real_fn(x, |v| v.abs() < 1.0, f64::atanh, |z| z.atanh()),
        ("ln", 1) => real_fn(x, |v| v > 0.0, f64::ln, |z| if z.norm() == 0.0 { c(f64::NEG_INFINITY) } else { z.ln() }),
        ("exp", 1) => x.exp(),
        ("sqrt", 1) => real_fn(x, |v| v >= 0.0, f64::sqrt, |z| z.sqrt()),
        ("abs", 1) => c(x.norm()),
        ("sign", 1) => {
            if x.im == 0.0 {
                c(if x.re > 0.0 { 1.0 } else if x.re < 0.0 { -1.0 } else { 0.0 })
            } else {
                x / x.norm()
            }
        }
        ("floor", 1) => c(x.re.floor()),
        ("ceil", 1) => c(x.re.ceil()),
        ("round", 1) => c(x.re.round()),
        ("round", 2) => {
            let p = 10f64.powi(r(1)? as i32);
            c((x.re * p).round() / p)
        }
        ("ipart", 1) => c(x.re.trunc()),
        ("fpart", 1) => c(x.re.fract()),
        ("factorial", 1) => gamma(x + 1.0),
        ("gamma", 1) => gamma(x),
        ("erf", 1) => c(erf(x.re)),
        ("erfc", 1) => c(1.0 - erf(x.re)),
        ("re", 1) => c(x.re),
        ("im", 1) => c(x.im),
        ("conj", 1) => x.conj(),
        ("arg", 1) => c(x.arg()),
        ("atan2", 2) => c(x.re.atan2(r(1)?)),
        ("mod", 2) => {
            let m = r(1)?;
            c(x.re - m * (x.re / m).floor())
        }
        ("ncr", 2) => {
            let k = r(1)?;
            (gamma(x + 1.0) / (gamma(c(k + 1.0)) * gamma(x - k + 1.0))).re.round().into()
        }
        ("npr", 2) => {
            let k = r(1)?;
            (gamma(x + 1.0) / gamma(x - k + 1.0)).re.round().into()
        }
        ("min", _) => c(a.iter().map(|v| v.re).fold(f64::INFINITY, f64::min)),
        ("max", _) => c(a.iter().map(|v| v.re).fold(f64::NEG_INFINITY, f64::max)),
        _ => return None,
    })
}

pub fn gamma(z: C) -> C {
    const G: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    if z.im == 0.0 && z.re > 0.0 && z.re <= 171.0 && z.re.fract() == 0.0 {
        return c((1..z.re as u64).fold(1.0, |p, k| p * k as f64));
    }
    if z.im == 0.0 && z.re > 171.7 {
        return c(f64::INFINITY);
    }
    if z.re < 0.5 {
        return PI / ((PI * z).sin() * gamma(1.0 - z));
    }
    let z = z - 1.0;
    let mut x = c(G[0]);
    for (i, g) in G.iter().enumerate().skip(1) {
        x += g / (z + i as f64);
    }
    let t = z + 7.5;
    (2.0 * PI).sqrt() * t.powc(z + 0.5) * (-t).exp() * x
}

pub fn lgamma(mut x: f64) -> f64 {
    let mut shift = 0.0;
    while x < 15.0 {
        shift -= x.ln();
        x += 1.0;
    }
    let x2 = x * x;
    shift + (x - 0.5) * x.ln() - x + 0.5 * (2.0 * PI).ln() + 1.0 / (12.0 * x) - 1.0 / (360.0 * x * x2)
        + 1.0 / (1260.0 * x2 * x2 * x)
        - 1.0 / (1680.0 * x2 * x2 * x2 * x)
}

/// Regularized lower incomplete gamma P(a, x).
pub fn gamma_p(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        let (mut sum, mut term, mut n) = (1.0 / a, 1.0 / a, a);
        for _ in 0..1000 {
            n += 1.0;
            term *= x / n;
            sum += term;
            if term.abs() < sum.abs() * 1e-16 {
                break;
            }
        }
        sum * (-x + a * x.ln() - lgamma(a)).exp()
    } else {
        let tiny = 1e-300;
        let mut b = x + 1.0 - a;
        let mut cc = 1.0 / tiny;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..1000 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < tiny {
                d = tiny;
            }
            cc = b + an / cc;
            if cc.abs() < tiny {
                cc = tiny;
            }
            d = 1.0 / d;
            let del = d * cc;
            h *= del;
            if (del - 1.0).abs() < 1e-16 {
                break;
            }
        }
        1.0 - (-x + a * x.ln() - lgamma(a)).exp() * h
    }
}

pub fn erf(x: f64) -> f64 {
    let p = gamma_p(0.5, x * x);
    if x < 0.0 { -p } else { p }
}

/// Regularized incomplete beta I_x(a, b).
pub fn beta_i(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (lgamma(a + b) - lgamma(a) - lgamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) { bt * betacf(a, b, x) / a } else { 1.0 - bt * betacf(b, a, 1.0 - x) / b }
}

fn betacf(a: f64, b: f64, x: f64) -> f64 {
    let tiny = 1e-300;
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let mut cc = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < tiny {
        d = tiny;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..2000 {
        let m = m as f64;
        let m2 = 2.0 * m;
        for aa in [m * (b - m) * x / ((qam + m2) * (a + m2)), -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))] {
            d = 1.0 + aa * d;
            if d.abs() < tiny {
                d = tiny;
            }
            cc = 1.0 + aa / cc;
            if cc.abs() < tiny {
                cc = tiny;
            }
            d = 1.0 / d;
            h *= d * cc;
        }
        if (d * cc - 1.0).abs() < 1e-16 {
            break;
        }
    }
    h
}

pub fn norm_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2))
}

pub fn inv_norm(p: f64) -> f64 {
    if !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }
    // Acklam's rational approximation, then one Halley step against the exact CDF.
    const A: [f64; 6] = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
    const B: [f64; 5] = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
    const CC: [f64; 6] = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
    const D: [f64; 4] = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
    let x = if p < 0.02425 {
        let q = (-2.0 * p.ln()).sqrt();
        (((((CC[0] * q + CC[1]) * q + CC[2]) * q + CC[3]) * q + CC[4]) * q + CC[5]) / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p > 1.0 - 0.02425 {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((CC[0] * q + CC[1]) * q + CC[2]) * q + CC[3]) * q + CC[4]) * q + CC[5]) / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    };
    if !x.is_finite() {
        return x;
    }
    let e = norm_cdf(x) - p;
    let u = e * (2.0 * PI).sqrt() * (x * x / 2.0).exp();
    x - u / (1.0 + x * u / 2.0)
}

pub fn t_cdf(t: f64, df: f64) -> f64 {
    let p = 0.5 * beta_i(df / 2.0, 0.5, df / (df + t * t));
    if t > 0.0 { 1.0 - p } else { p }
}

/// Brent's method on a bracketing interval.
pub fn brent(f: &dyn Fn(f64) -> f64, mut a: f64, mut b: f64) -> Option<f64> {
    let (mut fa, mut fb) = (f(a), f(b));
    if fa == 0.0 {
        return Some(a);
    }
    if fb == 0.0 {
        return Some(b);
    }
    if !(fa * fb < 0.0) {
        return None;
    }
    let (mut cx, mut fc) = (a, fa);
    let mut d = b - a;
    let mut e = d;
    for _ in 0..200 {
        if fb * fc > 0.0 {
            cx = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = cx;
            cx = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol = 2.0 * f64::EPSILON * b.abs() + 1e-300;
        let m = 0.5 * (cx - b);
        if m.abs() <= tol || fb == 0.0 {
            return Some(b);
        }
        if e.abs() >= tol && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (mut p, mut q);
            if a == cx {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                let qq = fa / fc;
                let r = fb / fc;
                p = s * (2.0 * m * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            if 2.0 * p < (3.0 * m * q - (tol * q).abs()).min((e * q).abs()) {
                e = d;
                d = p / q;
            } else {
                d = m;
                e = d;
            }
        } else {
            d = m;
            e = d;
        }
        a = b;
        fa = fb;
        b += if d.abs() > tol { d } else { tol.copysign(m) };
        fb = f(b);
        if fb.is_nan() {
            return None;
        }
    }
    Some(b)
}

/// Real roots of `f` on [lo, hi]: sign changes on a sample grid plus tangent roots at local minima of |f|.
pub fn find_roots(f: &dyn Fn(f64) -> f64, lo: f64, hi: f64, n: usize) -> Vec<f64> {
    let mut roots: Vec<f64> = Vec::new();
    let h = (hi - lo) / n as f64;
    let xs: Vec<f64> = (0..=n).map(|i| lo + h * i as f64).collect();
    let ys: Vec<f64> = xs.iter().map(|&x| f(x)).collect();
    let scale = ys.iter().filter(|y| y.is_finite()).fold(0.0f64, |m, y| m.max(y.abs())).max(1e-300);
    let push = |r: f64, roots: &mut Vec<f64>| {
        // Reject poles: the function must actually be small at the root.
        if f(r).abs() <= 1e-7 * scale.min(1e6).max(1.0) && !roots.iter().any(|q| (q - r).abs() <= 1e-9 * (1.0 + r.abs())) {
            roots.push(r);
        }
    };
    for i in 0..n {
        let (a, b, fa, fb) = (xs[i], xs[i + 1], ys[i], ys[i + 1]);
        if fa == 0.0 {
            push(a, &mut roots);
        } else if fa.is_finite() && fb.is_finite() && fa * fb < 0.0 {
            if let Some(r) = brent(f, a, b) {
                push(r, &mut roots);
            }
        } else if i > 0 && fa.is_finite() && fb.is_finite() && ys[i - 1].is_finite() && fa.abs() < ys[i - 1].abs() && fa.abs() <= fb.abs() {
            // local minimum of |f| at xs[i]: golden-section, accept if it reaches zero
            let (mut l, mut r) = (xs[i - 1], b);
            let g = 0.618_033_988_749_894_9;
            for _ in 0..100 {
                let (m1, m2) = (r - g * (r - l), l + g * (r - l));
                if f(m1).abs() < f(m2).abs() { r = m2 } else { l = m1 }
            }
            let m = 0.5 * (l + r);
            if f(m).abs() <= 1e-10 * scale.min(1.0) {
                push(m, &mut roots);
            }
        }
    }
    if ys[n] == 0.0 {
        push(xs[n], &mut roots);
    }
    roots.sort_by(f64::total_cmp);
    roots
}

fn tanh_sinh(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> (f64, f64) {
    let (mid, h) = (0.5 * (a + b), 0.5 * (b - a));
    let fin = |v: f64| if v.is_finite() { v } else { 0.0 };
    let mut sum = fin(f(mid)) * PI / 2.0;
    let mut step = 1.0;
    let (mut prev, mut err) = (f64::NAN, f64::INFINITY);
    for level in 0..11 {
        let mut k = 1;
        loop {
            let t = k as f64 * step;
            let u = PI / 2.0 * t.sinh();
            let delta = 1.0 / (u.exp() * u.cosh()); // 1 - tanh(u), computed without cancellation
            let w = PI / 2.0 * t.cosh() / (u.cosh() * u.cosh());
            let term = w * (fin(f(a + h * delta)) + fin(f(b - h * delta)));
            sum += term;
            if t > 6.5 || h * delta == 0.0 || (t > 3.0 && term.abs() <= 1e-18 * sum.abs()) {
                break;
            }
            k += if level == 0 { 1 } else { 2 };
        }
        let est = h * step * sum;
        err = (est - prev).abs();
        if level > 3 && err <= 1e-13 * est.abs().max(1e-300) {
            return (est, err);
        }
        prev = est;
        step /= 2.0;
    }
    (prev, err)
}

/// Definite integral with error estimate; handles infinite bounds and endpoint singularities.
pub fn quad(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> (f64, f64) {
    if a == b {
        return (0.0, 0.0);
    }
    if a > b {
        let (v, e) = quad(f, b, a);
        return (-v, e);
    }
    match (a.is_infinite(), b.is_infinite()) {
        (false, false) => tanh_sinh(f, a, b),
        (false, true) => tanh_sinh(&|t| f(a + t / (1.0 - t)) / ((1.0 - t) * (1.0 - t)), 0.0, 1.0),
        (true, false) => tanh_sinh(&|t| f(b - t / (1.0 - t)) / ((1.0 - t) * (1.0 - t)), 0.0, 1.0),
        (true, true) => {
            let (v1, e1) = quad(f, 0.0, f64::INFINITY);
            let (v2, e2) = quad(f, f64::NEG_INFINITY, 0.0);
            (v1 + v2, e1 + e2)
        }
    }
}

pub fn nderiv(f: &dyn Fn(f64) -> f64, x: f64) -> f64 {
    let h = 1e-4 * (1.0 + x.abs());
    (f(x - 2.0 * h) - 8.0 * f(x - h) + 8.0 * f(x + h) - f(x + 2.0 * h)) / (12.0 * h)
}

pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % p == 0 {
            return n == p;
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    let mulmod = |a: u64, b: u64| ((a as u128 * b as u128) % n as u128) as u64;
    'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let (mut x, mut e, mut base) = (1u64, d, a % n);
        while e > 0 {
            if e & 1 == 1 {
                x = mulmod(x, base);
            }
            base = mulmod(base, base);
            e >>= 1;
        }
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..s - 1 {
            x = mulmod(x, x);
            if x == n - 1 {
                continue 'w;
            }
        }
        return false;
    }
    true
}

