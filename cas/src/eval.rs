//! The calculator session: variables, modes, and command dispatch.
use crate::algebra;
use crate::calculus;
use crate::expr::*;
use crate::matrix;
use crate::numeric::{self, evalc};
use crate::parse::parse;
use crate::print::{Fmt, show};
use crate::simplify::*;
use crate::solve;
use crate::stats;
use crate::units;
use num_bigint::BigInt;
use num_complex::Complex64 as C;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Angle {
    Rad,
    Deg,
    Grad,
}

#[derive(Clone, Debug)]
pub struct Settings {
    pub angle: Angle,
    pub fmt: Fmt,
    /// evaluate everything to decimals
    pub approx: bool,
    /// allow non-real results and return complex roots from solve
    pub complex: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Settings { angle: Angle::Rad, fmt: Fmt::default(), approx: false, complex: false }
    }
}

pub struct Output {
    pub text: String,
    pub approx: Option<String>,
    pub note: Option<String>,
    pub value: Expr,
}

#[derive(Default)]
pub struct Calc {
    pub settings: Settings,
    pub vars: HashMap<String, Expr>,
    pub funcs: HashMap<String, (Vec<String>, Expr)>,
    pub ans: Option<Expr>,
    note: std::cell::RefCell<Option<String>>,
}

const COMMANDS: &[&str] = &[
    "simplify", "expand", "factor", "ratsimp", "partfrac", "solve", "csolve", "nsolve", "zeros", "czeros", "diff", "integrate", "nint", "nderiv", "limit", "taylor",
    "sum", "product", "det", "inv", "transpose", "rref", "eigenvalues", "eigenvectors", "identity", "dot", "cross", "norm", "trace", "rank", "dim", "convert", "approx",
    "exact", "subst", "numer", "denom", "isprime", "hex", "bin", "oct", "del", "odesolve", "fmin", "fmax", "seq", "len", "polar", "rect", "root", "log2", "atan2", "mod",
    "ncr", "npr", "gcd", "lcm", "min", "max", "int", "frac", "and", "or", "tangent", "arclen", "charpoly", "nullspace", "units", "constants", "nextprime", "degree",
    "coeffs", "trigexpand", "tolist", "augment", "unit",
];

/// (syntax, description) pairs for the function catalog.
pub const CATALOG: &[(&str, &str)] = &[
    ("solve(eq, x)", "solve an equation or inequality for x; exact when possible"),
    ("solve({eq1, eq2}, {x, y})", "solve a system of equations"),
    ("csolve(eq, x)", "solve including complex roots"),
    ("nsolve(eq, x, guess)", "numeric root near a starting guess"),
    ("zeros(expr, x)", "list of roots of an expression"),
    ("simplify(expr)", "simplify algebraic and trig expressions"),
    ("expand(expr)", "multiply out products and powers"),
    ("factor(expr)", "factor a polynomial or an integer"),
    ("ratsimp(expr)", "combine over a common denominator and cancel"),
    ("partfrac(expr, x)", "partial fraction decomposition"),
    ("diff(expr, x)", "derivative; diff(expr, x, n) for the n-th"),
    ("integrate(expr, x)", "antiderivative"),
    ("integrate(expr, x, a, b)", "definite integral; numeric when no closed form"),
    ("nint(expr, x, a, b)", "numeric integral"),
    ("nderiv(expr, x, a)", "numeric derivative at a point"),
    ("limit(expr, x, a)", "limit; limit(expr, x, a, 1) from the right, -1 from the left"),
    ("taylor(expr, x, n)", "Taylor polynomial of order n; taylor(expr, x, n, a) about a"),
    ("sum(expr, k, a, b)", "summation; closed form for polynomial and geometric terms"),
    ("product(expr, k, a, b)", "product over an integer range"),
    ("tangent(expr, x, a)", "tangent line at x = a"),
    ("arclen(expr, x, a, b)", "arc length of y = expr"),
    ("fmin(expr, x, a, b)", "minimum on an interval (fmax for maximum)"),
    ("odesolve(f, x, y, x0, y0, x1)", "y(x1) for y' = f(x, y), y(x0) = y0 (RK4)"),
    ("expr | x = value", "evaluate with a substitution"),
    ("x := value", "store a variable; f(x) := expr defines a function"),
    ("del(x)", "delete a variable or function"),
    ("[[1, 2], [3, 4]]", "matrix literal; also [1, 2; 3, 4]"),
    ("det(m)", "determinant"),
    ("inv(m)", "inverse; also m^-1"),
    ("transpose(m)", "transpose"),
    ("rref(m)", "reduced row echelon form"),
    ("eigenvalues(m)", "eigenvalues"),
    ("eigenvectors(m)", "eigenvectors, one row per eigenvalue"),
    ("charpoly(m, x)", "characteristic polynomial"),
    ("nullspace(m)", "basis of the null space"),
    ("rank(m)", "rank"),
    ("trace(m)", "trace"),
    ("identity(n)", "n by n identity matrix"),
    ("dot(u, v)", "dot product"),
    ("cross(u, v)", "cross product"),
    ("norm(v)", "Euclidean norm"),
    ("{1, 2, 3}", "list literal; arithmetic applies element by element"),
    ("seq(expr, k, a, b)", "list generated from a formula"),
    ("mean(list)", "mean; also median, stdev, stdevp, variance, total, sort"),
    ("onevar(list)", "one-variable statistics summary"),
    ("linreg(xs, ys)", "linear regression; also quadreg, expreg, lnreg, pwrreg"),
    ("normalcdf(lo, hi, mu, sigma)", "normal probability; normalpdf, invnorm"),
    ("tcdf(lo, hi, df)", "Student t probability; tpdf, invt"),
    ("chi2cdf(lo, hi, df)", "chi-square probability; fcdf(lo, hi, d1, d2)"),
    ("binompdf(n, p, k)", "binomial probability; binomcdf"),
    ("poissonpdf(lambda, k)", "Poisson probability; poissoncdf, geometpdf, geometcdf"),
    ("ncr(n, k)", "combinations; npr(n, k) permutations; n! factorial"),
    ("gcd(a, b)", "greatest common divisor; lcm(a, b); mod(a, b)"),
    ("isprime(n)", "primality test; nextprime(n); factor(n)"),
    ("convert(value, from, to)", "unit conversion, e.g. convert(5, mi, km); units() lists them"),
    ("constants()", "physical constants such as _c, _g, _h, _k, _Na"),
    ("approx(expr)", "decimal value; exact(x) converts a decimal to a fraction"),
    ("hex(n)", "base conversion; bin(n), oct(n); literals 0xFF, 0b101"),
    ("polar(z)", "{modulus, argument} of a complex number; rect(r, theta)"),
    ("re(z)", "real part; im(z), conj(z), abs(z), arg(z)"),
    ("sin(x)", "trig: sin cos tan sec csc cot, asin acos atan, sinh cosh tanh and inverses"),
    ("ln(x)", "natural log; log(x) base 10; log(x, b) base b; exp(x)"),
    ("sqrt(x)", "square root; cbrt(x); root(x, n); abs(x); sign(x)"),
    ("floor(x)", "floor; ceil, round(x, digits), ipart, fpart, min, max"),
    ("gamma(x)", "gamma function; erf(x)"),
];

pub fn is_builtin(name: &str) -> bool {
    COMMANDS.contains(&name) || UNARY_MATH.contains(&name) || stats::NAMES.contains(&name)
}

fn err<T>(s: impl Into<String>) -> Result<T, String> {
    Err(s.into())
}

fn raw_subst(e: &Expr, name: &str, val: &Expr) -> Expr {
    match e {
        Expr::Sym(s) if s == name => val.clone(),
        e => e.map(&mut |c| raw_subst(c, name, val)),
    }
}

fn sym_name(e: &Expr) -> Result<String, String> {
    match e {
        Expr::Sym(s) => Ok(s.clone()),
        _ => err("expected a variable name"),
    }
}

fn as_matrix(e: &Expr) -> Result<matrix::Mat, String> {
    match e {
        Expr::Matrix(m) => Ok(m.clone()),
        Expr::List(v) if !v.is_empty() => Ok(vec![v.clone()]),
        _ => err("expected a matrix like [[1, 2], [3, 4]]"),
    }
}

fn as_vector(e: &Expr) -> Result<Vec<Expr>, String> {
    match e {
        Expr::List(v) => Ok(v.clone()),
        Expr::Matrix(m) if m.len() == 1 => Ok(m[0].clone()),
        Expr::Matrix(m) if m.iter().all(|r| r.len() == 1) => Ok(m.iter().map(|r| r[0].clone()).collect()),
        _ => err("expected a vector like [1, 2, 3]"),
    }
}

fn real(e: &Expr) -> Result<f64, String> {
    if *e == INF {
        return Ok(f64::INFINITY);
    }
    if *e == neg(INF) {
        return Ok(f64::NEG_INFINITY);
    }
    if e.is_constant() { numeric::real_const(e) } else { None }.ok_or_else(|| "expected a real number".to_string())
}

fn to_exact(x: f64) -> Expr {
    // continued fraction with a tolerance near double precision
    if !x.is_finite() {
        return Expr::Float(x);
    }
    let (mut h0, mut h1, mut k0, mut k1) = (BigInt::zero(), BigInt::one(), BigInt::one(), BigInt::zero());
    let mut v = x;
    for _ in 0..40 {
        let a = v.floor();
        let ai = BigInt::from(a as i64);
        let (h2, k2) = (&ai * &h1 + &h0, &ai * &k1 + &k0);
        h0 = h1;
        h1 = h2;
        k0 = k1;
        k1 = k2;
        let approx = h1.to_f64().unwrap() / k1.to_f64().unwrap();
        if (approx - x).abs() <= 1e-13 * x.abs().max(1e-300) || (v - a).abs() < 1e-12 {
            break;
        }
        v = 1.0 / (v - a);
        if v.abs() > 1e15 {
            break;
        }
    }
    Expr::Num(Rat::new(h1, k1))
}

fn exactify(e: &Expr) -> Expr {
    match e {
        Expr::Float(x) => to_exact(*x),
        Expr::Add(v) => add_v(v.iter().map(exactify).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(exactify).collect()),
        Expr::Pow(b, x) => pow(exactify(b), exactify(x)),
        Expr::Func(n, v) => func(n, v.iter().map(exactify).collect()),
        e => e.map(&mut |c| exactify(c)),
    }
}

fn radix_str(n: &BigInt, radix: u32, prefix: &str) -> String {
    format!("{}{}{}", if n.is_negative() { "-" } else { "" }, prefix, n.abs().to_str_radix(radix).to_uppercase())
}

impl Calc {
    pub fn new() -> Calc {
        Calc::default()
    }

    fn knows(&self, name: &str) -> bool {
        is_builtin(name) || self.funcs.contains_key(name)
    }

    fn set_note(&self, s: &str) {
        *self.note.borrow_mut() = Some(s.to_string());
    }

    fn angle_factor(&self) -> Option<Expr> {
        match self.settings.angle {
            Angle::Rad => None,
            Angle::Deg => Some(mul(rat(1, 180), PI)),
            Angle::Grad => Some(mul(rat(1, 200), PI)),
        }
    }

    fn lookup(&self, name: &str, shadow: &[String]) -> Expr {
        if shadow.iter().any(|s| s == name) {
            return sym(name);
        }
        match name {
            "pi" => return PI,
            "e" => return E,
            "i" => return I,
            "inf" | "infinity" => return INF,
            "undef" => return UNDEF,
            "true" => return Expr::Bool(true),
            "false" => return Expr::Bool(false),
            "ans" => return self.ans.clone().unwrap_or(int(0)),
            _ => {}
        }
        if let Some(v) = self.vars.get(name) {
            return v.clone();
        }
        if let Some((_, v, _)) = units::CONSTANTS.iter().find(|(n, _, _)| *n == name) {
            return if v.fract() == 0.0 && v.abs() < 1e15 { int(*v as i64) } else { Expr::Float(*v) };
        }
        sym(name)
    }

    /// Evaluate a raw tree. `shadow` lists names that must stay symbolic.
    fn ev(&self, e: &Expr, shadow: &[String]) -> Result<Expr, String> {
        Ok(match e {
            Expr::Sym(s) => self.lookup(s, shadow),
            Expr::Add(v) => add_v(v.iter().map(|x| self.ev(x, shadow)).collect::<Result<_, _>>()?),
            Expr::Mul(v) => mul_v(v.iter().map(|x| self.ev(x, shadow)).collect::<Result<_, _>>()?),
            Expr::Pow(b, x) => pow(self.ev(b, shadow)?, self.ev(x, shadow)?),
            Expr::Rel(r, a, b) => {
                let (a, b) = (self.ev(a, shadow)?, self.ev(b, shadow)?);
                if a.is_constant() && b.is_constant() {
                    if let (Some(x), Some(y)) = (numeric::real_const(&a), numeric::real_const(&b)) {
                        let d = sub(a.clone(), b.clone());
                        let (x, y) = if is_zero_num(&expand(&d)) { (0.0, 0.0) } else { (x, y) };
                        return Ok(Expr::Bool(r.holds(x, y)));
                    }
                }
                Expr::Rel(*r, Box::new(a), Box::new(b))
            }
            Expr::List(v) => Expr::List(v.iter().map(|x| self.ev(x, shadow)).collect::<Result<_, _>>()?),
            Expr::Matrix(m) => Expr::Matrix(m.iter().map(|r| r.iter().map(|x| self.ev(x, shadow)).collect::<Result<_, _>>()).collect::<Result<_, _>>()?),
            Expr::Func(name, args) => self.call(name, args, shadow)?,
            other => other.clone(),
        })
    }

    fn call(&self, name: &str, raw: &[Expr], shadow: &[String]) -> Result<Expr, String> {
        // forms that bind a variable: evaluate with that name kept symbolic
        const BINDERS: &[&str] = &[
            "diff", "integrate", "nint", "nderiv", "limit", "taylor", "sum", "product", "solve", "csolve", "nsolve", "zeros", "czeros", "partfrac", "fmin", "fmax", "seq",
            "tangent", "arclen", "charpoly", "degree", "coeffs",
        ];
        let need = |n: usize| if raw.len() < n { err(format!("{name} needs at least {n} argument(s)")) } else { Ok(()) };
        if name == "|" {
            need(2)?;
            let conds: Vec<Expr> = match &raw[1] {
                Expr::Func(n, v) if n == "and" => v.clone(),
                c => vec![c.clone()],
            };
            let mut names = shadow.to_vec();
            let mut pairs = Vec::new();
            for c in &conds {
                let Expr::Rel(Rel::Eq, l, r) = c else { return err("use: expr | x = value") };
                let n = sym_name(l)?;
                pairs.push((n.clone(), self.ev(r, shadow)?));
                names.push(n);
            }
            let mut v = self.ev(&raw[0], &names)?;
            for (n, val) in &pairs {
                v = subst(&v, n, val);
            }
            return Ok(v);
        }
        if name == "convert" {
            need(3)?;
            let v = self.ev(&raw[0], shadow)?;
            let (from, to) = (sym_name(&raw[1])?, sym_name(&raw[2])?);
            return Ok(Expr::Float(units::convert(real(&v)?, &from, &to)?));
        }
        if name == "del" {
            return err("del works at the top level: del(x)");
        }
        let mut names = shadow.to_vec();
        let mut bound: Vec<String> = Vec::new();
        if BINDERS.contains(&name) && raw.len() >= 2 {
            match &raw[1] {
                Expr::Sym(s) => bound.push(s.clone()),
                Expr::List(v) => {
                    for x in v {
                        bound.push(sym_name(x)?);
                    }
                }
                _ => {}
            }
            names.extend(bound.iter().cloned());
        }
        let a: Vec<Expr> = raw.iter().map(|x| self.ev(x, &names)).collect::<Result<_, _>>()?;

        // user-defined functions
        if let Some((params, body)) = self.funcs.get(name) {
            if params.len() != a.len() {
                return err(format!("{name} takes {} argument(s)", params.len()));
            }
            let mut b = body.clone();
            // two-phase rename avoids capture when an argument mentions another parameter
            for (i, p) in params.iter().enumerate() {
                b = raw_subst(&b, p, &sym(&format!("@p{i}")));
            }
            let mut inner = names.clone();
            for i in 0..params.len() {
                inner.push(format!("@p{i}"));
            }
            let mut v = self.ev(&b, &inner)?;
            for (i, val) in a.iter().enumerate() {
                v = subst(&v, &format!("@p{i}"), val);
            }
            return Ok(v);
        }

        // angle mode
        if let (Some(k), 1) = (self.angle_factor(), a.len()) {
            match name {
                "sin" | "cos" | "tan" | "sec" | "csc" | "cot" => return Ok(func(name, vec![mul(a[0].clone(), k)])),
                "asin" | "acos" | "atan" | "asec" | "acsc" | "acot" | "arg" => return Ok(div(func(name, a), k)),
                _ => {}
            }
        }

        if let Some(r) = stats::call(name, &a) {
            return r;
        }

        let var = || -> Result<String, String> {
            match bound.first() {
                Some(v) => Ok(v.clone()),
                None => {
                    let s = a[0].free_syms();
                    let s: Vec<&String> = s.iter().filter(|n| !n.starts_with('@')).collect();
                    if s.len() == 1 { Ok(s[0].clone()) } else { err(format!("{name}: say which variable, e.g. {name}(expr, x)")) }
                }
            }
        };
        let cx = self.settings.complex;
        Ok(match name {
            "simplify" => {
                need(1)?;
                algebra::simplify(&a[0])
            }
            "expand" => {
                need(1)?;
                expand(&a[0])
            }
            "trigexpand" => {
                need(1)?;
                algebra::trig_linearize(&a[0])
            }
            "factor" => {
                need(1)?;
                algebra::factor(&a[0])
            }
            "ratsimp" => {
                need(1)?;
                algebra::ratsimp(&a[0])
            }
            "numer" | "denom" => {
                need(1)?;
                let (n, d) = algebra::together(&a[0]);
                let (n, d) = algebra::cancel(&n, &d);
                if name == "numer" { n } else { d }
            }
            "partfrac" => {
                need(1)?;
                algebra::partfrac(&a[0], &var()?).ok_or("not a rational function with rational coefficients")?
            }
            "degree" | "coeffs" => {
                need(1)?;
                let c = crate::poly::coeffs_expr(&a[0], &var()?).ok_or("not a polynomial")?;
                if name == "degree" { int(c.len() as i64 - 1) } else { Expr::List(c.into_iter().rev().collect()) }
            }
            "subst" => {
                need(3)?;
                replace(&a[0], &a[1], &a[2])
            }
            "approx" => {
                need(1)?;
                to_float(&a[0])
            }
            "exact" => {
                need(1)?;
                exactify(&a[0])
            }
            "solve" | "csolve" | "zeros" | "czeros" => {
                need(1)?;
                let complex = cx || name.starts_with('c');
                let list = name.ends_with("zeros");
                // systems
                let eqs: Option<Vec<Expr>> = match &a[0] {
                    Expr::List(v) => Some(v.clone()),
                    Expr::Func(n, v) if n == "and" && bound.len() > 1 => Some(v.clone()),
                    _ => None,
                };
                if let Some(eqs) = eqs {
                    let vars: Vec<String> = if bound.is_empty() { eqs.iter().flat_map(|e| e.free_syms()).collect::<std::collections::BTreeSet<_>>().into_iter().collect() } else { bound.clone() };
                    let fs: Vec<Expr> = eqs
                        .iter()
                        .map(|e| match e {
                            Expr::Rel(Rel::Eq, l, r) => Ok(sub((**l).clone(), (**r).clone())),
                            Expr::Rel(..) => err("systems of inequalities are not supported"),
                            Expr::Bool(true) => Ok(int(0)),
                            Expr::Bool(false) => Ok(int(1)),
                            e => Ok(e.clone()),
                        })
                        .collect::<Result<_, _>>()?;
                    let (sols, approx) = solve::solve_system(&fs, &vars, complex)?;
                    if approx {
                        self.set_note("numeric solution");
                    }
                    if list {
                        return Ok(Expr::Matrix(sols).pipe_nonempty());
                    }
                    return Ok(func("or", sols.into_iter().map(|row| func("and", vars.iter().zip(row).map(|(v, e)| Expr::Rel(Rel::Eq, Box::new(sym(v)), Box::new(e))).collect())).collect()));
                }
                let x = var()?;
                let (f, rel) = match &a[0] {
                    Expr::Rel(r, l, rr) => (sub((**l).clone(), (**rr).clone()), *r),
                    Expr::Bool(b) => return Ok(Expr::Bool(*b)),
                    e => (e.clone(), Rel::Eq),
                };
                if rel != Rel::Eq && rel != Rel::Ne {
                    return solve::solve_inequality(&f, rel, &x);
                }
                let (sol, approx) = solve::solve(&f, &x, complex)?;
                if approx {
                    self.set_note("numeric solution");
                }
                if rel == Rel::Ne {
                    return Ok(match sol {
                        solve::Solution::All => Expr::Bool(false),
                        solve::Solution::Set(v) => func("and", v.into_iter().map(|r| Expr::Rel(Rel::Ne, Box::new(sym(&x)), Box::new(r))).collect()),
                    });
                }
                if list {
                    match sol {
                        solve::Solution::All => return err("every value is a zero"),
                        solve::Solution::Set(v) => Expr::List(v),
                    }
                } else {
                    solve::solution_expr(&x, sol)
                }
            }
            "nsolve" => {
                need(2)?;
                let x = var()?;
                let f = match &a[0] {
                    Expr::Rel(Rel::Eq, l, r) => sub((**l).clone(), (**r).clone()),
                    e => e.clone(),
                };
                let g = |v: f64| numeric::evalf(&f, &x, v);
                let guess = a.get(2).map_or(Ok(0.0), real)?;
                let root = if let (Some(lo), Some(hi)) = (a.get(2), a.get(3)) {
                    numeric::find_roots(&g, real(lo)?, real(hi)?, 2000).first().copied()
                } else {
                    // Newton from the guess, then a widening bracket search
                    let mut xv = guess;
                    let mut ok = false;
                    for _ in 0..100 {
                        let (fv, dv) = (g(xv), numeric::nderiv(&g, xv));
                        if fv.abs() < 1e-13 {
                            ok = true;
                            break;
                        }
                        if !fv.is_finite() || !dv.is_finite() || dv == 0.0 {
                            break;
                        }
                        xv -= fv / dv;
                    }
                    if ok && g(xv).abs() < 1e-9 {
                        Some(xv)
                    } else {
                        (0..8).find_map(|k| {
                            let w = 10f64 * 4f64.powi(k);
                            let r = numeric::find_roots(&g, guess - w, guess + w, 4000);
                            r.into_iter().min_by(|p, q| (p - guess).abs().total_cmp(&(q - guess).abs()))
                        })
                    }
                };
                match root {
                    Some(r) => Expr::Rel(Rel::Eq, Box::new(sym(&x)), Box::new(Expr::Float(r))),
                    None => return err("no root found"),
                }
            }
            "diff" => {
                need(1)?;
                let x = var()?;
                let n = a.get(2).map_or(Some(1), |n| n.as_i64()).filter(|n| (0..=50).contains(n)).ok_or("order must be a small non-negative integer")?;
                let mut d = a[0].clone();
                for _ in 0..n {
                    d = calculus::diff(&d, &x);
                }
                d
            }
            "nderiv" => {
                need(3)?;
                let x = var()?;
                Expr::Float(numeric::nderiv(&|v| numeric::evalf(&a[0], &x, v), real(&a[2])?))
            }
            "integrate" => {
                need(1)?;
                let x = var()?;
                if a.len() >= 4 {
                    let (v, approx) = calculus::definite(&a[0], &x, &a[2], &a[3])?;
                    if approx {
                        self.set_note("numeric integral");
                    }
                    v
                } else {
                    calculus::antiderivative(&a[0], &x).ok_or("no closed-form antiderivative found; use integrate(expr, x, a, b) for a numeric value")?
                }
            }
            "nint" => {
                need(4)?;
                let x = var()?;
                let (v, e) = numeric::quad(&|t| numeric::evalf(&a[0], &x, t), real(&a[2])?, real(&a[3])?);
                if !(e <= 1e-6 * (1.0 + v.abs())) {
                    return err("integral does not converge");
                }
                Expr::Float(v)
            }
            "arclen" => {
                need(4)?;
                let x = var()?;
                let d = calculus::diff(&a[0], &x);
                let integrand = sqrt(add(int(1), pow(d, int(2))));
                let (v, approx) = calculus::definite(&algebra::simplify(&integrand), &x, &a[2], &a[3])?;
                if approx {
                    self.set_note("numeric integral");
                }
                v
            }
            "tangent" => {
                need(3)?;
                let x = var()?;
                let slope = subst(&calculus::diff(&a[0], &x), &x, &a[2]);
                let y0 = subst(&a[0], &x, &a[2]);
                expand(&add(mul(slope, sub(sym(&x), a[2].clone())), y0))
            }
            "limit" => {
                need(3)?;
                let x = var()?;
                let dir = a.get(3).map_or(Ok(0.0), real)? as i32;
                let r = calculus::limit(&a[0], &x, &a[2], dir).ok_or("limit not found")?;
                if matches!(r, Expr::Float(_)) {
                    self.set_note("numeric limit");
                }
                r
            }
            "taylor" => {
                need(3)?;
                let x = var()?;
                let n = a[2].as_i64().filter(|n| (0..=40).contains(n)).ok_or("order must be an integer from 0 to 40")?;
                let at = a.get(3).cloned().unwrap_or(int(0));
                calculus::taylor(&a[0], &x, &at, n as usize).ok_or("series not found")?
            }
            "sum" | "product" => {
                if a.len() == 1 {
                    let v = as_vector(&a[0])?;
                    return Ok(if name == "sum" { add_v(v) } else { mul_v(v) });
                }
                need(4)?;
                let r = calculus::sum(&a[0], &var()?, &a[2], &a[3], name == "product")?;
                if matches!(r, Expr::Float(_)) && !a[0].has_float() {
                    self.set_note("numeric sum");
                }
                r
            }
            "seq" => {
                need(4)?;
                let k = var()?;
                let (lo, hi) = (a[2].as_i64().ok_or("seq bounds must be integers")?, a[3].as_i64().ok_or("seq bounds must be integers")?);
                if hi - lo > 100_000 {
                    return err("too many terms");
                }
                Expr::List((lo..=hi).map(|i| subst(&a[0], &k, &int(i))).collect())
            }
            "fmin" | "fmax" => {
                need(4)?;
                let x = var()?;
                let sgn = if name == "fmin" { 1.0 } else { -1.0 };
                let g = |v: f64| sgn * numeric::evalf(&a[0], &x, v);
                let (lo, hi) = (real(&a[2])?, real(&a[3])?);
                // coarse scan, then golden-section refinement around the best sample
                let n = 2000;
                let best = (0..=n).map(|i| lo + (hi - lo) * i as f64 / n as f64).filter(|v| g(*v).is_finite()).min_by(|p, q| g(*p).total_cmp(&g(*q))).ok_or("function undefined on the interval")?;
                let h = (hi - lo) / n as f64;
                let (mut l, mut r) = ((best - h).max(lo), (best + h).min(hi));
                let gr = 0.618_033_988_749_894_9;
                for _ in 0..200 {
                    let (m1, m2) = (r - gr * (r - l), l + gr * (r - l));
                    if g(m1) < g(m2) { r = m2 } else { l = m1 }
                }
                let mut xm = 0.5 * (l + r);
                // Newton polish on the derivative (golden section stalls near sqrt(eps))
                for _ in 0..3 {
                    let h = 1e-4 * (1.0 + xm.abs());
                    let (d1, d2) = ((g(xm + h) - g(xm - h)) / (2.0 * h), (g(xm + h) - 2.0 * g(xm) + g(xm - h)) / (h * h));
                    let step = d1 / d2;
                    if step.is_finite() && step.abs() < h && xm - step > lo && xm - step < hi {
                        xm -= step;
                    }
                }
                Expr::List(vec![Expr::Rel(Rel::Eq, Box::new(sym(&x)), Box::new(Expr::Float(xm))), Expr::Rel(Rel::Eq, Box::new(sym("y")), Box::new(Expr::Float(sgn * g(xm))))])
            }
            "odesolve" => {
                need(6)?;
                let (x, y) = (sym_name(&raw[1])?, sym_name(&raw[2])?);
                let f = self.ev(&raw[0], &[names.clone(), vec![x.clone(), y.clone()]].concat())?;
                let pts = calculus::ode_rk4(&f, &x, &y, real(&a[3])?, real(&a[4])?, real(&a[5])?, 2000);
                Expr::Float(pts.last().unwrap().1)
            }
            "det" => {
                need(1)?;
                matrix::det(&as_matrix(&a[0])?).ok_or("det needs a square matrix")?
            }
            "inv" => {
                need(1)?;
                Expr::Matrix(matrix::inverse(&as_matrix(&a[0])?).ok_or("matrix is singular or not square")?)
            }
            "transpose" => {
                need(1)?;
                Expr::Matrix(matrix::transpose(&as_matrix(&a[0])?))
            }
            "rref" => {
                need(1)?;
                Expr::Matrix(matrix::rref(&as_matrix(&a[0])?).0)
            }
            "rank" => {
                need(1)?;
                int(matrix::rref(&as_matrix(&a[0])?).1.len() as i64)
            }
            "trace" => {
                need(1)?;
                matrix::trace(&as_matrix(&a[0])?).ok_or("trace needs a square matrix")?
            }
            "dim" | "len" => {
                need(1)?;
                match &a[0] {
                    Expr::List(v) => int(v.len() as i64),
                    Expr::Matrix(m) => Expr::List(vec![int(m.len() as i64), int(m[0].len() as i64)]),
                    _ => return err("expected a list or matrix"),
                }
            }
            "identity" => {
                need(1)?;
                let n = a[0].as_i64().filter(|n| (1..=64).contains(n)).ok_or("identity(n) needs n from 1 to 64")?;
                Expr::Matrix(matrix::identity(n as usize))
            }
            "nullspace" => {
                need(1)?;
                let v = matrix::nullspace(&as_matrix(&a[0])?);
                if v.is_empty() { Expr::List(vec![]) } else { Expr::Matrix(v) }
            }
            "tolist" => {
                need(1)?;
                Expr::List(as_matrix(&a[0])?.into_iter().flatten().collect())
            }
            "augment" => {
                need(2)?;
                let (p, q) = (as_matrix(&a[0])?, as_matrix(&a[1])?);
                if p.len() != q.len() {
                    return err("augment needs matrices with the same number of rows");
                }
                Expr::Matrix(p.into_iter().zip(q).map(|(mut r, s)| {
                    r.extend(s);
                    r
                }).collect())
            }
            "charpoly" | "eigenvalues" | "eigenvectors" => {
                need(1)?;
                let m = as_matrix(&a[0])?;
                let n = m.len();
                if m.iter().any(|r| r.len() != n) {
                    return err("needs a square matrix");
                }
                let lam = if name == "charpoly" { bound.first().cloned().unwrap_or("x".into()) } else { "@lam".to_string() };
                let shifted = |l: &Expr| -> matrix::Mat { (0..n).map(|i| (0..n).map(|j| if i == j { sub(m[i][j].clone(), l.clone()) } else { m[i][j].clone() }).collect()).collect() };
                // det(A - lambda I) by cofactors needs symbolic entries; lambda makes it symbolic
                let cp = expand(&matrix::det(&shifted(&sym(&lam))).unwrap());
                if name == "charpoly" {
                    return Ok(if n % 2 == 1 { expand(&neg(cp)) } else { cp });
                }
                let (sol, approx) = solve::solve(&cp, &lam, true)?;
                let solve::Solution::Set(vals) = sol else { return err("eigenvalues not found") };
                if approx {
                    self.set_note("numeric eigenvalues");
                }
                if name == "eigenvalues" {
                    Expr::List(vals)
                } else {
                    let mut rows = Vec::new();
                    for l in &vals {
                        let ns = matrix::nullspace(&shifted(l));
                        if ns.is_empty() && approx {
                            // numeric eigenvalue: relax by rounding the shifted matrix's smallest pivot
                            return err("eigenvectors need exact eigenvalues; try a matrix with exact entries");
                        }
                        rows.extend(ns);
                    }
                    if rows.is_empty() {
                        return err("eigenvectors not found");
                    }
                    Expr::Matrix(rows)
                }
            }
            "dot" => {
                need(2)?;
                let (u, v) = (as_vector(&a[0])?, as_vector(&a[1])?);
                if u.len() != v.len() {
                    return err("dot needs vectors of equal length");
                }
                add_v(u.into_iter().zip(v).map(|(p, q)| mul(p, q)).collect())
            }
            "cross" => {
                need(2)?;
                let (u, v) = (as_vector(&a[0])?, as_vector(&a[1])?);
                if u.len() != 3 || v.len() != 3 {
                    return err("cross needs two 3-vectors");
                }
                let c = |i: usize, j: usize| sub(mul(u[i].clone(), v[j].clone()), mul(u[j].clone(), v[i].clone()));
                Expr::Matrix(vec![vec![c(1, 2), c(2, 0), c(0, 1)]])
            }
            "norm" => {
                need(1)?;
                let v: Vec<Expr> = match &a[0] {
                    Expr::Matrix(m) => m.iter().flatten().cloned().collect(),
                    other => as_vector(other)?,
                };
                sqrt(add_v(v.into_iter().map(|p| pow(f1("abs", p), int(2))).collect()))
            }
            "unit" => {
                need(1)?;
                let v = as_vector(&a[0])?;
                let n = sqrt(add_v(v.iter().map(|p| pow(p.clone(), int(2))).collect()));
                Expr::Matrix(vec![v.into_iter().map(|p| div(p, n.clone())).collect()])
            }
            "isprime" => {
                need(1)?;
                let n = a[0].as_int().and_then(|n| n.to_u64()).ok_or("isprime needs a non-negative integer below 2^64")?;
                Expr::Bool(numeric::is_prime_u64(n))
            }
            "nextprime" => {
                need(1)?;
                let mut n = a[0].as_int().and_then(|n| n.to_u64()).ok_or("nextprime needs a non-negative integer below 2^64")? + 1;
                while !numeric::is_prime_u64(n) {
                    n += 1;
                }
                big(BigInt::from(n))
            }
            "hex" | "bin" | "oct" => {
                need(1)?;
                let n = a[0].as_int().ok_or("base conversion needs an integer")?;
                Expr::Str(match name {
                    "hex" => radix_str(n, 16, "0x"),
                    "bin" => radix_str(n, 2, "0b"),
                    _ => radix_str(n, 8, "0o"),
                })
            }
            "polar" => {
                need(1)?;
                let k = self.angle_factor();
                let ang = f1("arg", a[0].clone());
                Expr::List(vec![f1("abs", a[0].clone()), match k {
                    Some(k) => div(ang, k),
                    None => ang,
                }])
            }
            "rect" => {
                need(2)?;
                let th = match self.angle_factor() {
                    Some(k) => mul(a[1].clone(), k),
                    None => a[1].clone(),
                };
                expand(&mul(a[0].clone(), add(f1("cos", th.clone()), mul(I, f1("sin", th)))))
            }
            "atan2" => {
                need(2)?;
                let r = match (real(&a[0]), real(&a[1])) {
                    (Ok(y), Ok(x)) if a.iter().any(|v| v.has_float()) => Expr::Float(y.atan2(x)),
                    _ => f1("arg", add(a[1].clone(), mul(I, a[0].clone()))),
                };
                match self.angle_factor() {
                    Some(k) => div(r, k),
                    None => r,
                }
            }
            "units" => Expr::Str(units::unit_list()),
            "constants" => Expr::Str(units::CONSTANTS.iter().map(|(n, v, d)| format!("{n} = {v:e}  ({d})")).collect::<Vec<_>>().join("\n")),
            _ => func(name, a),
        })
    }

    fn check_real(&self, raw: &Expr, v: &Expr) -> Result<(), String> {
        if self.settings.complex {
            return Ok(());
        }
        let is_i = |e: &Expr| matches!(e, Expr::Const(Const::I)) || matches!(e, Expr::Sym(s) if s == "i") || matches!(e, Expr::Func(n, _) if n == "csolve" || n == "czeros" || n == "eigenvalues" || n == "eigenvectors" || n == "rect");
        if v.has(&|e| matches!(e, Expr::Const(Const::I))) && !raw.has(&is_i) {
            return err("non-real result (switch to Complex mode to see it)");
        }
        Ok(())
    }

    /// Parse and evaluate without touching session state. `free` names stay symbolic.
    pub fn eval_str(&self, src: &str, free: &[&str]) -> Result<Expr, String> {
        let raw = parse(src, &|n| self.knows(n))?;
        let shadow: Vec<String> = free.iter().map(|s| s.to_string()).collect();
        let v = self.ev(&raw, &shadow)?;
        Ok(if self.settings.approx { to_float(&v) } else { v })
    }

    pub fn show(&self, e: &Expr) -> String {
        show(e, &self.settings.fmt)
    }

    /// Evaluate one input line, updating variables and `ans`.
    pub fn eval_line(&mut self, src: &str) -> Result<Output, String> {
        *self.note.borrow_mut() = None;
        let raw = parse(src, &|n| self.knows(n))?;
        if let Expr::Func(n, args) = &raw {
            if n == "del" {
                let mut removed = Vec::new();
                for a in args {
                    let name = sym_name(a)?;
                    if self.vars.remove(&name).is_some() | self.funcs.remove(&name).is_some() {
                        removed.push(name);
                    }
                }
                let text = if removed.is_empty() { "nothing to delete".to_string() } else { format!("deleted {}", removed.join(", ")) };
                return Ok(Output { text: text.clone(), approx: None, note: None, value: Expr::Str(text) });
            }
            if n == ":=" {
                match &args[0] {
                    Expr::Sym(name) => {
                        if ["pi", "e", "i", "inf", "undef", "ans", "true", "false"].contains(&name.as_str()) || is_builtin(name) {
                            return err(format!("'{name}' is reserved"));
                        }
                        let v = self.ev(&args[1], &[])?;
                        let v = if self.settings.approx { to_float(&v) } else { v };
                        if v.has_sym(name) {
                            return err("a variable cannot be defined in terms of itself");
                        }
                        self.vars.insert(name.clone(), v.clone());
                        self.funcs.remove(name);
                        return Ok(self.output(v));
                    }
                    Expr::Func(name, params) => {
                        if is_builtin(name) {
                            return err(format!("'{name}' is a built-in function"));
                        }
                        let params: Vec<String> = params.iter().map(sym_name).collect::<Result<_, _>>()?;
                        self.funcs.insert(name.clone(), (params.clone(), args[1].clone()));
                        self.vars.remove(name);
                        // show the definition with the body evaluated symbolically where possible
                        let body = self.ev(&args[1], &params).unwrap_or_else(|_| args[1].clone());
                        let text = format!("{}({}) := {}", name, params.join(", "), self.show(&body));
                        return Ok(Output { text, approx: None, note: None, value: body });
                    }
                    _ => return err("bad assignment"),
                }
            }
        }
        let v = tidy_constant(self.ev(&raw, &[])?);
        self.check_real(&raw, &v)?;
        let v = if self.settings.approx { to_float(&v) } else { v };
        self.ans = Some(v.clone());
        Ok(self.output(v))
    }

    fn output(&self, v: Expr) -> Output {
        let text = self.show(&v);
        let approx = if matches!(v, Expr::Str(_) | Expr::Bool(_)) {
            None
        } else {
            fn wants(e: &Expr) -> bool {
                match e {
                    Expr::Rel(_, l, r) if matches!(**l, Expr::Sym(_)) => wants(r),
                    Expr::Rel(_, l, r) => wants(l) || wants(r),
                    Expr::List(_) | Expr::Matrix(_) => e.children().iter().any(|c| wants(c)),
                    Expr::Func(n, a) if n == "and" || n == "or" => a.iter().any(wants),
                    Expr::Num(r) => !r.is_integer(),
                    Expr::Float(_) => false,
                    e => e.is_constant(),
                }
            }
            let f = to_float(&v);
            let s = self.show(&f);
            // "2." and "2" carry the same information
            let strip = |t: &str| {
                let c: Vec<char> = t.chars().collect();
                (0..c.len()).filter(|&i| !(c[i] == '.' && !c.get(i + 1).is_some_and(|d| d.is_ascii_digit()))).map(|i| c[i]).collect::<String>()
            };
            (wants(&v) && strip(&s) != strip(&text) && !f.is_undef()).then_some(s)
        };
        Output { text, approx, note: self.note.borrow().clone(), value: v }
    }

    /// Numeric function of one variable for plotting.
    pub fn plot_fn(&self, src: &str, var: &str) -> Result<Expr, String> {
        let e = self.eval_str(src, &[var])?;
        let extra: Vec<String> = e.free_syms().into_iter().filter(|s| s != var).collect();
        if !extra.is_empty() {
            return err(format!("undefined: {}", extra.join(", ")));
        }
        evalc(&e, &[(var, C::new(0.123, 0.0))]).ok_or("cannot be plotted")?;
        Ok(e)
    }
}

/// Constants: multiply out products of sums and rationalize a + b*sqrt(c) denominators when that is shorter.
fn tidy_constant(v: Expr) -> Expr {
    match &v {
        Expr::List(_) | Expr::Matrix(_) => return v.map(&mut |c| tidy_constant(c.clone())),
        Expr::Rel(r, a, b) => return Expr::Rel(*r, a.clone(), Box::new(tidy_constant((**b).clone()))),
        Expr::Func(n, a) if n == "and" || n == "or" => return Expr::Func(n.clone(), a.iter().map(|c| tidy_constant(c.clone())).collect()),
        _ => {}
    }
    if !v.is_constant() || v.has_float() || v.size() > 60 || !v.has(&|e| matches!(e, Expr::Add(_))) {
        return v;
    }
    let mut best = v.clone();
    let ex = expand(&v);
    if ex.size() < best.size() {
        best = ex;
    }
    let r = algebra::rationalize(&v);
    if r.size() < best.size() {
        best = r;
    }
    best
}

trait PipeNonEmpty {
    fn pipe_nonempty(self) -> Expr;
}
impl PipeNonEmpty for Expr {
    fn pipe_nonempty(self) -> Expr {
        match self {
            Expr::Matrix(m) if m.is_empty() => Expr::List(vec![]),
            e => e,
        }
    }
}

