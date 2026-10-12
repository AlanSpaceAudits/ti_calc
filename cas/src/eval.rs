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
    pub figure: Option<crate::figure::Figure>,
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
    "coeffs", "trigexpand", "tolist", "augment", "unit", "dms", "todms", "todd", "torad", "desolve", "lagrange", "hamiltonian", "hamilton", "fieldeq", "fieldsol", "fieldcheck", "orbperiod", "visviva", "vescape", "parallaxdist", "distmod", "absmag",
];

const GEOMETRY: &[&str] = &[
    "inscribedangle", "centralangle", "chordlen", "tangentchord", "sagitta", "arclength", "sectorarea", "segmentarea", "chordsangle", "secantsangle", "powerpoint", "vecadd",
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
    ("desolve(eq, x, y)", "solve an ODE written with y' and y''; desolve(eq, x, y, x0, y0, yp0) with initial conditions"),
    ("lagrange(L, q, t)", "Euler-Lagrange equation of motion; write the velocity as q'; lagrange(L, {x, y}, t) for several coordinates"),
    ("hamiltonian(L, q)", "Hamiltonian H(q, p) from a Lagrangian; momenta are p, or p_x, p_y for several coordinates"),
    ("hamilton(H, q, p)", "Hamilton's equations q' = dH/dp, p' = -dH/dq"),
    ("fieldeq()", "list the standard field equations; fieldeq(wave) shows one with its frame and solutions"),
    ("fieldsol(wave, 1)", "closed-form solution number n of a field equation"),
    ("fieldcheck(wave, u)", "substitute u into a field equation: true, or the residual that is left"),
    ("galilean(expr, v)", "Galilean boost x -> x - v*t; lorentzt(expr, v) is the Lorentz boost"),
    ("lorentz(v)", "Lorentz factor; beta(v), rapidity(v); c is the symbol c, or _c when v is a number"),
    ("timedil(t, v)", "time dilation; lencon(L, v) length contraction; veladd(u, v) velocity addition"),
    ("relmom(m, v)", "relativistic momentum; relenergy(m, v), relke(m, v), restenergy(m)"),
    ("reldoppler(f, v)", "relativistic Doppler shift for a source approaching at speed v"),
    ("schwarzschild(M)", "Schwarzschild radius; gravdil(t, M, r) gravitational time dilation"),
    ("inscribedangle(central)", "inscribed angle theorem: half the central angle on the same arc"),
    ("centralangle(chord, r)", "central angle that a chord subtends; chordlen(r, angle) is the reverse"),
    ("tangentchord(arc)", "tangent-chord theorem: half the intercepted arc"),
    ("chordsangle(arc1, arc2)", "angle between two chords that cross inside a circle; secantsangle for outside"),
    ("sagitta(r, chord)", "height of an arc; arclength(r, angle), sectorarea(r, angle), segmentarea(r, angle)"),
    ("powerpoint(d, r)", "power of a point at distance d from the centre of a circle of radius r"),
    ("vecadd(u, v)", "vector sum, drawn as a parallelogram for 2D vectors"),
    ("jd(y, m, d, h, min, s)", "Julian date (UT); caldate(jd) converts back"),
    ("gmst(jd)", "Greenwich mean sidereal time as an angle; lst(jd, east_longitude)"),
    ("hms(h, m, s)", "angle from hours, minutes, seconds of time; tohms(angle) shows one"),
    ("eq2hor(ha, dec, lat)", "altitude and azimuth from hour angle and declination; hor2eq(alt, az, lat)"),
    ("eq2ecl(ra, dec)", "equatorial to ecliptic coordinates; ecl2eq(lon, lat); obliquity(jd)"),
    ("angsep(ra1, dec1, ra2, dec2)", "angular separation of two directions"),
    ("sunpos(jd)", "low-precision Sun: right ascension, declination, ecliptic longitude, distance"),
    ("riseset(dec, lat)", "hour angle and azimuth at rising and setting, hours above the horizon"),
    ("refraction(alt)", "atmospheric refraction at an apparent altitude; airmass(alt)"),
    ("kepler(M, e)", "eccentric anomaly from mean anomaly; trueanomaly(E, e)"),
    ("orbperiod(a, mu)", "orbital period; visviva(mu, r, a) orbital speed; vescape(mu, r)"),
    ("parallaxdist(p)", "distance in parsecs from parallax in arcseconds; distmod(d), absmag(m, d)"),
    ("sphsss(a, b, c)", "spherical triangle from three sides; also sphsas(a, C, b), sphasa(A, c, B), sphaaa(A, B, C)"),
    ("sphexcess(A, B, C)", "spherical excess; spharea(A, B, C, R) triangle area"),
    ("sphdist(lat1, lon1, lat2, lon2)", "great-circle central angle; add a radius as fifth argument for distance"),
    ("sphbearing(lat1, lon1, lat2, lon2)", "initial great-circle bearing from north"),
    ("navtri(lat, dec, lha)", "navigational triangle: altitude and azimuth of a body"),
    ("todms(x)", "show an angle as degrees, minutes, seconds; todd(x) decimal degrees; torad(x) radians"),
    ("dms(d, m, s)", "enter an angle in degrees, minutes, seconds; also 12°34'56\""),
    ("polar(z)", "{modulus, argument} of a complex number; rect(r, theta)"),
    ("re(z)", "real part; im(z), conj(z), abs(z), arg(z)"),
    ("sin(x)", "trig: sin cos tan sec csc cot, asin acos atan, sinh cosh tanh and inverses"),
    ("ln(x)", "natural log; log(x) base 10; log(x, b) base b; exp(x)"),
    ("sqrt(x)", "square root; cbrt(x); root(x, n); abs(x); sign(x)"),
    ("floor(x)", "floor; ceil, round(x, digits), ipart, fpart, min, max"),
    ("gamma(x)", "gamma function; erf(x)"),
];

pub fn is_builtin(name: &str) -> bool {
    COMMANDS.contains(&name) || crate::astro::NAMES.contains(&name) || GEOMETRY.contains(&name) || crate::physics::REL_NAMES.contains(&name) || crate::physics::SPH_NAMES.contains(&name) || UNARY_MATH.contains(&name) || stats::NAMES.contains(&name)
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

    /// Degrees per unit of the current angle mode.
    pub fn degrees_per_unit(&self) -> f64 {
        match self.settings.angle {
            Angle::Rad => 180.0 / std::f64::consts::PI,
            Angle::Deg => 1.0,
            Angle::Grad => 0.9,
        }
    }

    /// Evaluate a raw tree with the given names kept symbolic.
    pub(crate) fn ev_free(&self, raw: &Expr, free: &[&str]) -> Result<Expr, String> {
        let shadow: Vec<String> = free.iter().map(|s| s.to_string()).collect();
        self.ev(raw, &shadow)
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
        // mechanics and ODE forms: coordinates and their primed derivatives stay symbolic
        if matches!(name, "desolve" | "lagrange" | "hamiltonian" | "hamilton") {
            need(2)?;
            let names_of = |e: &Expr| -> Result<Vec<String>, String> {
                match e {
                    Expr::List(v) => v.iter().map(sym_name).collect(),
                    e => Ok(vec![sym_name(e)?]),
                }
            };
            let mut hold = shadow.to_vec();
            let primed = |ns: &[String], hold: &mut Vec<String>| {
                for n in ns {
                    hold.extend([n.clone(), format!("{n}'"), format!("{n}''")]);
                }
            };
            match name {
                "desolve" => {
                    need(3)?;
                    let (x, y) = (sym_name(&raw[1])?, sym_name(&raw[2])?);
                    hold.push(x.clone());
                    primed(std::slice::from_ref(&y), &mut hold);
                    let eq = self.ev(&raw[0], &hold)?;
                    let ics: Vec<Expr> = raw[3..].iter().map(|e| self.ev(e, shadow)).collect::<Result<_, _>>()?;
                    if ics.len() == 1 {
                        return err("initial conditions: desolve(eq, x, y, x0, y0) or desolve(eq, x, y, x0, y0, yp0)");
                    }
                    return crate::ode::desolve(&eq, &x, &y, &ics);
                }
                "lagrange" => {
                    need(3)?;
                    let (qs, t) = (names_of(&raw[1])?, sym_name(&raw[2])?);
                    hold.push(t.clone());
                    primed(&qs, &mut hold);
                    let l = self.ev(&raw[0], &hold)?;
                    let mut eqs = crate::physics::lagrange(&l, &qs, &t);
                    return Ok(if eqs.len() == 1 { eqs.pop().unwrap() } else { Expr::List(eqs) });
                }
                "hamiltonian" => {
                    let qs = names_of(&raw[1])?;
                    primed(&qs, &mut hold);
                    for q in &qs {
                        hold.push(crate::physics::momentum_name(q, qs.len() == 1));
                    }
                    let l = self.ev(&raw[0], &hold)?;
                    return crate::physics::hamiltonian(&l, &qs);
                }
                _ => {
                    need(3)?;
                    let (qs, ps) = (names_of(&raw[1])?, names_of(&raw[2])?);
                    if qs.len() != ps.len() {
                        return err("hamilton needs one momentum per coordinate");
                    }
                    primed(&qs, &mut hold);
                    primed(&ps, &mut hold);
                    let h = self.ev(&raw[0], &hold)?;
                    return Ok(Expr::List(crate::physics::hamilton(&h, &qs, &ps)));
                }
            }
        }
        if matches!(name, "fieldeq" | "fieldsol" | "fieldcheck") {
            if raw.is_empty() {
                return Ok(Expr::Str(crate::physics::field_list()));
            }
            let f = crate::physics::field(&sym_name(&raw[0])?)?;
            return Ok(match name {
                "fieldeq" => Expr::Str(crate::physics::field_info(f, &|e| self.show(e))),
                "fieldsol" => {
                    let sols = (f.solutions)();
                    let n = match raw.get(1) {
                        Some(e) => self.ev(e, shadow)?.as_i64().ok_or("solution number must be an integer")?,
                        None => 1,
                    };
                    sols.get((n - 1).max(0) as usize).filter(|_| n >= 1).map(|s| s.1.clone()).ok_or(format!("{} has {} solution(s)", f.name, sols.len()))?
                }
                _ => {
                    need(2)?;
                    // the field variables stay symbolic even if the session has values for them
                    let mut hold = shadow.to_vec();
                    hold.extend(["x", "y", "z", "t", "theta"].map(String::from));
                    match crate::physics::field_check(f, &self.ev(&raw[1], &hold)?) {
                        None => Expr::Bool(true),
                        Some(r) => Expr::Rel(Rel::Eq, Box::new(sym("residual")), Box::new(r)),
                    }
                }
            });
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
        if let Some(r) = crate::astro::call(name, &a, self.degrees_per_unit()) {
            return r;
        }
        if GEOMETRY.contains(&name) {
            // circle theorems; angles and arcs are in the current angle mode
            let k = self.angle_factor().unwrap_or(int(1));
            let g = |i: usize| a.get(i).cloned().ok_or_else(|| format!("{name}: missing argument"));
            let half = rat(1, 2);
            return Ok(match name {
                "inscribedangle" | "tangentchord" => mul(half, g(0)?),
                "chordsangle" => mul(half, add(g(0)?, g(1)?)),
                "secantsangle" => mul(half, sub(g(0)?, g(1)?)),
                "centralangle" => div(mul(int(2), f1("asin", div(g(0)?, mul(int(2), g(1)?)))), k),
                "chordlen" => mul_v(vec![int(2), g(0)?, f1("sin", mul_v(vec![half, g(1)?, k]))]),
                "sagitta" => sub(g(0)?, sqrt(sub(pow(g(0)?, int(2)), mul(rat(1, 4), pow(g(1)?, int(2)))))),
                "arclength" => mul_v(vec![g(0)?, g(1)?, k]),
                "sectorarea" => mul_v(vec![half, pow(g(0)?, int(2)), g(1)?, k]),
                "segmentarea" => mul_v(vec![half, pow(g(0)?, int(2)), sub(mul(g(1)?, k.clone()), f1("sin", mul(g(1)?, k)))]),
                "powerpoint" => sub(pow(g(0)?, int(2)), pow(g(1)?, int(2))),
                _ => add(g(0)?, g(1)?),
            });
        }
        if crate::physics::SPH_NAMES.contains(&name) {
            let k = self.angle_factor().unwrap_or(int(1));
            if let Some(r) = crate::physics::spherical(name, &a, &k) {
                return r;
            }
        }
        if crate::physics::REL_NAMES.contains(&name) {
            // speed of light: the symbol c for symbolic input, the SI value for numeric input
            // the speed arguments decide: a numeric speed means SI units
            let speed_args: &[usize] = match name {
                "lorentz" | "beta" | "rapidity" => &[0],
                "veladd" => &[0, 1],
                "restenergy" | "schwarzschild" => &[0],
                "gravdil" => &[1, 2],
                _ => &[1],
            };
            let numeric = speed_args.iter().all(|i| a.get(*i).is_some_and(|e| e.is_constant()));
            let c = match self.vars.get("c") {
                Some(v) => v.clone(),
                None if numeric => int(299_792_458),
                None => sym("c"),
            };
            let g = match self.vars.get("G") {
                Some(v) => v.clone(),
                None if numeric => Expr::Float(6.674_30e-11),
                None => sym("G"),
            };
            if let Some(r) = crate::physics::relativity(name, &a, &c, &g) {
                return r;
            }
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
            "dms" => {
                need(1)?;
                // sign of the degrees applies to the whole angle
                let d = a[0].clone();
                let frac = add(div(a.get(1).cloned().unwrap_or(int(0)), int(60)), div(a.get(2).cloned().unwrap_or(int(0)), int(3600)));
                let deg = if d.has_neg_coeff() { sub(d, frac) } else { add(d, frac) };
                match self.settings.angle {
                    Angle::Deg => deg,
                    Angle::Rad => mul(deg, mul(rat(1, 180), PI)),
                    Angle::Grad => mul(deg, rat(10, 9)),
                }
            }
            "todd" | "torad" | "todms" => {
                need(1)?;
                // the argument is an angle in the current angle mode
                let rad = match self.angle_factor() {
                    Some(k) => mul(a[0].clone(), k),
                    None => a[0].clone(),
                };
                if name == "torad" {
                    return Ok(rad);
                }
                let deg = mul(rad, div(int(180), PI));
                if name == "todd" {
                    return Ok(if deg.as_rat().is_some_and(|r| r.is_integer()) { deg } else { to_float(&deg) });
                }
                let v = real(&deg)?;
                // work in hundred-thousandths of an arcsecond so rounding carries correctly
                let total = (v.abs() * 3600.0 * 1e5).round() as u128;
                let (d, m, s) = (total / 360_000_000, total / 6_000_000 % 60, (total % 6_000_000) as f64 / 1e5);
                let secs = format!("{s:.5}");
                let secs = secs.trim_end_matches('0').trim_end_matches('.');
                Expr::Str(format!("{}{}°{}'{}\"", if v < 0.0 { "-" } else { "" }, d, m, secs))
            }
            "orbperiod" => {
                // orbperiod(semi-major axis, gravitational parameter mu = G*M)
                need(2)?;
                mul_v(vec![int(2), PI, sqrt(div(pow(a[0].clone(), int(3)), a[1].clone()))])
            }
            "visviva" => {
                // visviva(mu, r, a): orbital speed
                need(3)?;
                sqrt(mul(a[0].clone(), sub(div(int(2), a[1].clone()), pow(a[2].clone(), int(-1)))))
            }
            "vescape" => {
                need(2)?;
                sqrt(div(mul(int(2), a[0].clone()), a[1].clone()))
            }
            "parallaxdist" => {
                // parallax in arcseconds -> distance in parsecs
                need(1)?;
                pow(a[0].clone(), int(-1))
            }
            "distmod" => {
                need(1)?;
                sub(mul(int(5), func("log", vec![a[0].clone()])), int(5))
            }
            "absmag" => {
                // absmag(apparent magnitude, distance in parsecs)
                need(2)?;
                add(sub(a[0].clone(), mul(int(5), func("log", vec![a[1].clone()]))), int(5))
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
        let is_i = |e: &Expr| matches!(e, Expr::Const(Const::I)) || matches!(e, Expr::Sym(s) if s == "i") || matches!(e, Expr::Func(n, _) if ["csolve", "czeros", "eigenvalues", "eigenvectors", "rect", "fieldsol", "fieldcheck", "fieldeq", "hamiltonian", "hamilton", "lagrange", "desolve"].contains(&n.as_str()));
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
                return Ok(Output { figure: None, text: text.clone(), approx: None, note: None, value: Expr::Str(text) });
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
                        let figure = crate::figure::build(self, &args[1], &v);
                        let mut out = self.output(v);
                        out.figure = figure;
                        return Ok(out);
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
                        return Ok(Output { figure: None, text, approx: None, note: None, value: body });
                    }
                    _ => return err("bad assignment"),
                }
            }
        }
        let v = tidy_constant(self.ev(&raw, &[])?);
        self.check_real(&raw, &v)?;
        let v = if self.settings.approx { to_float(&v) } else { v };
        self.ans = Some(v.clone());
        let figure = crate::figure::build(self, &raw, &v);
        let mut out = self.output(v);
        out.figure = figure;
        Ok(out)
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
        Output { figure: None, text, approx, note: self.note.borrow().clone(), value: v }
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

