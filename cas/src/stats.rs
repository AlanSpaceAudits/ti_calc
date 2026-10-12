//! List statistics, regressions and probability distributions.
use crate::expr::*;
use crate::numeric::{self as nm, beta_i, gamma_p, lgamma, norm_cdf};
use crate::simplify::*;
use std::f64::consts::PI;

fn list(e: &Expr) -> Result<&Vec<Expr>, String> {
    match e {
        Expr::List(v) if !v.is_empty() => Ok(v),
        Expr::Matrix(m) if m.len() == 1 && !m[0].is_empty() => Ok(&m[0]),
        _ => Err("expected a non-empty list like {1, 2, 3}".into()),
    }
}

fn floats(e: &Expr) -> Result<Vec<f64>, String> {
    list(e)?.iter().map(|x| if x.is_constant() { nm::real_const(x) } else { None }.ok_or_else(|| "list must contain real numbers".to_string())).collect()
}

fn num(e: &Expr) -> Result<f64, String> {
    let inf = *e == INF;
    let ninf = *e == neg(INF);
    if inf {
        return Ok(f64::INFINITY);
    }
    if ninf {
        return Ok(f64::NEG_INFINITY);
    }
    if e.is_constant() { nm::real_const(e) } else { None }.ok_or_else(|| "expected a real number".to_string())
}

fn eq(name: &str, v: Expr) -> Expr {
    Expr::Rel(Rel::Eq, Box::new(sym(name)), Box::new(v))
}

fn mean(v: &[Expr]) -> Expr {
    div(add_v(v.to_vec()), int(v.len() as i64))
}

fn variance(v: &[Expr], sample: bool) -> Expr {
    let m = mean(v);
    let ss = add_v(v.iter().map(|x| pow(sub(x.clone(), m.clone()), int(2))).collect());
    expand(&div(ss, int(v.len() as i64 - sample as i64)))
}

fn median(v: &[Expr]) -> Result<Expr, String> {
    let mut s: Vec<(f64, Expr)> = v.iter().map(|x| num(x).map(|f| (f, x.clone()))).collect::<Result<_, _>>()?;
    s.sort_by(|a, b| a.0.total_cmp(&b.0));
    let n = s.len();
    Ok(if n % 2 == 1 { s[n / 2].1.clone() } else { div(add(s[n / 2 - 1].1.clone(), s[n / 2].1.clone()), int(2)) })
}

fn quantile(sorted: &[f64], lo: usize, hi: usize) -> f64 {
    // median of sorted[lo..hi] (TI quartile method: medians of the lower and upper halves)
    let n = hi - lo;
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 { sorted[lo + n / 2] } else { 0.5 * (sorted[lo + n / 2 - 1] + sorted[lo + n / 2]) }
}

fn lin_fit(x: &[f64], y: &[f64]) -> (f64, f64, f64) {
    let n = x.len() as f64;
    let (sx, sy) = (x.iter().sum::<f64>(), y.iter().sum::<f64>());
    let sxx: f64 = x.iter().map(|v| v * v).sum();
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
    let syy: f64 = y.iter().map(|v| v * v).sum();
    let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    let icpt = (sy - slope * sx) / n;
    let r = (n * sxy - sx * sy) / ((n * sxx - sx * sx) * (n * syy - sy * sy)).sqrt();
    (slope, icpt, r)
}

fn xy(args: &[Expr]) -> Result<(Vec<f64>, Vec<f64>), String> {
    if args.len() != 2 {
        return Err("expected two lists: xs, ys".into());
    }
    let (x, y) = (floats(&args[0])?, floats(&args[1])?);
    if x.len() != y.len() || x.len() < 2 {
        return Err("lists must have the same length (at least 2)".into());
    }
    Ok((x, y))
}

fn fit_result(model: Expr, r: f64) -> Expr {
    Expr::List(vec![eq("y", model), eq("r", Expr::Float(r)), eq("r2", Expr::Float(r * r))])
}

fn tail(args: &[Expr], i: usize, default: f64) -> Result<f64, String> {
    args.get(i).map_or(Ok(default), num)
}

fn ln_choose(n: f64, k: f64) -> f64 {
    lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma(n - k + 1.0)
}

fn binom_pdf(n: f64, p: f64, k: f64) -> f64 {
    if k < 0.0 || k > n || k.fract() != 0.0 {
        return 0.0;
    }
    if p == 0.0 {
        return (k == 0.0) as i32 as f64;
    }
    if p == 1.0 {
        return (k == n) as i32 as f64;
    }
    (ln_choose(n, k) + k * p.ln() + (n - k) * (1.0 - p).ln()).exp()
}

fn poisson_pdf(l: f64, k: f64) -> f64 {
    if k < 0.0 || k.fract() != 0.0 { 0.0 } else { (k * l.ln() - l - lgamma(k + 1.0)).exp() }
}

fn inv_cdf(cdf: &dyn Fn(f64) -> f64, p: f64) -> f64 {
    let (mut lo, mut hi) = (-1.0, 1.0);
    while cdf(lo) > p && lo > -1e12 {
        lo *= 2.0;
    }
    while cdf(hi) < p && hi < 1e12 {
        hi *= 2.0;
    }
    nm::brent(&|x| cdf(x) - p, lo, hi).unwrap_or(f64::NAN)
}

/// Returns None when `name` is not a statistics function.
pub fn call(name: &str, args: &[Expr]) -> Option<Result<Expr, String>> {
    let a = |i: usize| args.get(i).ok_or_else(|| format!("{name}: missing argument"));
    let f = |v: f64| Ok(Expr::Float(v));
    Some((|| -> Result<Expr, String> {
        match name {
            "mean" => Ok(mean(list(a(0)?)?)),
            "median" => median(list(a(0)?)?),
            "variance" => Ok(variance(list(a(0)?)?, true)),
            "stdev" => Ok(sqrt(variance(list(a(0)?)?, true))),
            "stdevp" => Ok(sqrt(variance(list(a(0)?)?, false))),
            "total" => Ok(add_v(list(a(0)?)?.clone())),
            "sort" => {
                let mut s: Vec<(f64, Expr)> = list(a(0)?)?.iter().map(|x| num(x).map(|f| (f, x.clone()))).collect::<Result<_, _>>()?;
                s.sort_by(|p, q| p.0.total_cmp(&q.0));
                Ok(Expr::List(s.into_iter().map(|p| p.1).collect()))
            }
            "onevar" => {
                let v = list(a(0)?)?;
                let mut s = floats(a(0)?)?;
                s.sort_by(f64::total_cmp);
                let n = s.len();
                Ok(Expr::List(vec![
                    eq("mean", mean(v)),
                    eq("sum", add_v(v.clone())),
                    eq("sumsq", add_v(v.iter().map(|x| pow(x.clone(), int(2))).collect())),
                    eq("sx", if n > 1 { sqrt(variance(v, true)) } else { UNDEF }),
                    eq("sigmax", sqrt(variance(v, false))),
                    eq("n", int(n as i64)),
                    eq("min", func("min", v.clone())),
                    eq("q1", Expr::Float(quantile(&s, 0, n / 2))),
                    eq("med", median(v)?),
                    eq("q3", Expr::Float(quantile(&s, n.div_ceil(2), n))),
                    eq("max", func("max", v.clone())),
                ]))
            }
            "linreg" => {
                let (x, y) = xy(args)?;
                let (m, b, r) = lin_fit(&x, &y);
                Ok(fit_result(add(mul(Expr::Float(m), sym("x")), Expr::Float(b)), r))
            }
            "lnreg" => {
                let (x, y) = xy(args)?;
                let lx: Vec<f64> = x.iter().map(|v| v.ln()).collect();
                let (m, b, r) = lin_fit(&lx, &y);
                Ok(fit_result(add(Expr::Float(b), mul(Expr::Float(m), f1("ln", sym("x")))), r))
            }
            "expreg" => {
                let (x, y) = xy(args)?;
                let ly: Vec<f64> = y.iter().map(|v| v.ln()).collect();
                let (m, b, r) = lin_fit(&x, &ly);
                Ok(fit_result(mul(Expr::Float(b.exp()), pow(Expr::Float(m.exp()), sym("x"))), r))
            }
            "pwrreg" => {
                let (x, y) = xy(args)?;
                let lx: Vec<f64> = x.iter().map(|v| v.ln()).collect();
                let ly: Vec<f64> = y.iter().map(|v| v.ln()).collect();
                let (m, b, r) = lin_fit(&lx, &ly);
                Ok(fit_result(mul(Expr::Float(b.exp()), pow(sym("x"), Expr::Float(m))), r))
            }
            "quadreg" => {
                let (x, y) = xy(args)?;
                if x.len() < 3 {
                    return Err("quadreg needs at least 3 points".into());
                }
                // normal equations for a x^2 + b x + c
                let s = |p: i32| x.iter().map(|v| v.powi(p)).sum::<f64>();
                let t = |p: i32| x.iter().zip(&y).map(|(v, w)| v.powi(p) * w).sum::<f64>();
                let mut m = [[s(4), s(3), s(2), t(2)], [s(3), s(2), s(1), t(1)], [s(2), s(1), s(0), t(0)]];
                for c in 0..3 {
                    let p = (c..3).max_by(|i, j| m[*i][c].abs().total_cmp(&m[*j][c].abs())).unwrap();
                    m.swap(c, p);
                    for r in 0..3 {
                        if r != c {
                            let k = m[r][c] / m[c][c];
                            for cc in c..4 {
                                m[r][cc] -= k * m[c][cc];
                            }
                        }
                    }
                }
                let (qa, qb, qc) = (m[0][3] / m[0][0], m[1][3] / m[1][1], m[2][3] / m[2][2]);
                let ybar = y.iter().sum::<f64>() / y.len() as f64;
                let ss_res: f64 = x.iter().zip(&y).map(|(v, w)| (w - (qa * v * v + qb * v + qc)).powi(2)).sum();
                let ss_tot: f64 = y.iter().map(|w| (w - ybar).powi(2)).sum();
                let xs = sym("x");
                let model = add_v(vec![mul(Expr::Float(qa), pow(xs.clone(), int(2))), mul(Expr::Float(qb), xs), Expr::Float(qc)]);
                Ok(Expr::List(vec![eq("y", model), eq("r2", Expr::Float(1.0 - ss_res / ss_tot))]))
            }
            "normalpdf" => {
                let (x, mu, sd) = (num(a(0)?)?, tail(args, 1, 0.0)?, tail(args, 2, 1.0)?);
                f((-0.5 * ((x - mu) / sd).powi(2)).exp() / (sd * (2.0 * PI).sqrt()))
            }
            "normalcdf" => {
                let (lo, hi, mu, sd) = (num(a(0)?)?, num(a(1)?)?, tail(args, 2, 0.0)?, tail(args, 3, 1.0)?);
                f(norm_cdf((hi - mu) / sd) - norm_cdf((lo - mu) / sd))
            }
            "invnorm" => {
                let (p, mu, sd) = (num(a(0)?)?, tail(args, 1, 0.0)?, tail(args, 2, 1.0)?);
                f(mu + sd * nm::inv_norm(p))
            }
            "tpdf" => {
                let (x, v) = (num(a(0)?)?, num(a(1)?)?);
                f((lgamma((v + 1.0) / 2.0) - lgamma(v / 2.0)).exp() / (v * PI).sqrt() * (1.0 + x * x / v).powf(-(v + 1.0) / 2.0))
            }
            "tcdf" => {
                let (lo, hi, v) = (num(a(0)?)?, num(a(1)?)?, num(a(2)?)?);
                f(nm::t_cdf(hi, v) - nm::t_cdf(lo, v))
            }
            "invt" => {
                let (p, v) = (num(a(0)?)?, num(a(1)?)?);
                f(inv_cdf(&|x| nm::t_cdf(x, v), p))
            }
            "chi2pdf" => {
                let (x, k) = (num(a(0)?)?, num(a(1)?)?);
                f(if x < 0.0 { 0.0 } else { ((k / 2.0 - 1.0) * x.ln() - x / 2.0 - (k / 2.0) * 2f64.ln() - lgamma(k / 2.0)).exp() })
            }
            "chi2cdf" => {
                let (lo, hi, k) = (num(a(0)?)?, num(a(1)?)?, num(a(2)?)?);
                let c = |x: f64| if x <= 0.0 { 0.0 } else if x.is_infinite() { 1.0 } else { gamma_p(k / 2.0, x / 2.0) };
                f(c(hi) - c(lo))
            }
            "fcdf" => {
                let (lo, hi, d1, d2) = (num(a(0)?)?, num(a(1)?)?, num(a(2)?)?, num(a(3)?)?);
                let c = |x: f64| if x <= 0.0 { 0.0 } else if x.is_infinite() { 1.0 } else { beta_i(d1 / 2.0, d2 / 2.0, d1 * x / (d1 * x + d2)) };
                f(c(hi) - c(lo))
            }
            "binompdf" => f(binom_pdf(num(a(0)?)?, num(a(1)?)?, num(a(2)?)?)),
            "binomcdf" => {
                let (n, p, k) = (num(a(0)?)?, num(a(1)?)?, num(a(2)?)?);
                f((0..=k.floor().min(n) as i64).map(|j| binom_pdf(n, p, j as f64)).sum())
            }
            "poissonpdf" => f(poisson_pdf(num(a(0)?)?, num(a(1)?)?)),
            "poissoncdf" => {
                let (l, k) = (num(a(0)?)?, num(a(1)?)?);
                f((0..=k.floor() as i64).map(|j| poisson_pdf(l, j as f64)).sum())
            }
            "geometpdf" => {
                let (p, k) = (num(a(0)?)?, num(a(1)?)?);
                f(if k < 1.0 { 0.0 } else { p * (1.0 - p).powf(k - 1.0) })
            }
            "geometcdf" => {
                let (p, k) = (num(a(0)?)?, num(a(1)?)?);
                f(1.0 - (1.0 - p).powf(k.floor()))
            }
            _ => Err(String::new()),
        }
    })())
    .filter(|r| !matches!(r, Err(e) if e.is_empty()))
}

pub const NAMES: &[&str] = &[
    "mean", "median", "variance", "stdev", "stdevp", "total", "sort", "onevar", "linreg", "lnreg", "expreg", "pwrreg", "quadreg", "normalpdf", "normalcdf", "invnorm",
    "tpdf", "tcdf", "invt", "chi2pdf", "chi2cdf", "fcdf", "binompdf", "binomcdf", "poissonpdf", "poissoncdf", "geometpdf", "geometcdf",
];
