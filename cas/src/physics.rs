//! Relativity formulas, Lagrangian and Hamiltonian mechanics, and a library of standard field equations.
//! The first seven field equations and their frames follow Murgueitio Ramirez (2024).
use crate::algebra;
use crate::calculus::diff;
use crate::expr::*;
use crate::numeric::evalc;
use crate::poly::coeffs_expr;
use crate::simplify::*;
use crate::solve;
use num_complex::Complex64 as C;

fn s(n: &str) -> Expr {
    sym(n)
}
fn sq(e: Expr) -> Expr {
    pow(e, int(2))
}

/// Lorentz factor 1/sqrt(1 - v^2/c^2).
pub fn lorentz(v: &Expr, c: &Expr) -> Expr {
    pow(sub(int(1), div(sq(v.clone()), sq(c.clone()))), rat(-1, 2))
}

pub const REL_NAMES: &[&str] = &[
    "lorentz", "beta", "rapidity", "timedil", "lencon", "veladd", "relmom", "relenergy", "relke", "restenergy", "reldoppler", "lorentzt", "galilean", "schwarzschild", "gravdil",
];

/// Relativity functions. `c` is the speed of light to use; `g` the gravitational constant.
pub fn relativity(name: &str, a: &[Expr], c: &Expr, g: &Expr) -> Option<Result<Expr, String>> {
    let need = |n: usize| if a.len() < n { Err(format!("{name} needs {n} argument(s)")) } else { Ok(()) };
    Some((|| -> Result<Expr, String> {
        Ok(match name {
            "lorentz" => {
                need(1)?;
                lorentz(&a[0], c)
            }
            "beta" => {
                need(1)?;
                div(a[0].clone(), c.clone())
            }
            "rapidity" => {
                need(1)?;
                f1("atanh", div(a[0].clone(), c.clone()))
            }
            "timedil" => {
                need(2)?;
                mul(a[0].clone(), lorentz(&a[1], c))
            }
            "lencon" => {
                need(2)?;
                div(a[0].clone(), lorentz(&a[1], c))
            }
            "veladd" => {
                need(2)?;
                algebra::ratsimp(&div(add(a[0].clone(), a[1].clone()), add(int(1), div(mul(a[0].clone(), a[1].clone()), sq(c.clone())))))
            }
            "relmom" => {
                need(2)?;
                mul_v(vec![lorentz(&a[1], c), a[0].clone(), a[1].clone()])
            }
            "relenergy" => {
                need(2)?;
                mul_v(vec![lorentz(&a[1], c), a[0].clone(), sq(c.clone())])
            }
            "relke" => {
                need(2)?;
                mul_v(vec![sub(lorentz(&a[1], c), int(1)), a[0].clone(), sq(c.clone())])
            }
            "restenergy" => {
                need(1)?;
                mul(a[0].clone(), sq(c.clone()))
            }
            "reldoppler" => {
                // observed frequency for a source approaching at speed v (negative v: receding)
                need(2)?;
                let b = div(a[1].clone(), c.clone());
                mul(a[0].clone(), sqrt(div(add(int(1), b.clone()), sub(int(1), b))))
            }
            "galilean" => {
                // x -> x - v*t
                need(2)?;
                subst(&a[0], "x", &sub(s("x"), mul(a[1].clone(), s("t"))))
            }
            "lorentzt" => {
                // x -> gamma*(x - v*t), t -> gamma*(t - v*x/c^2), applied together
                need(2)?;
                let gam = lorentz(&a[1], c);
                let nx = mul(gam.clone(), sub(s("x"), mul(a[1].clone(), s("t"))));
                let nt = mul(gam, sub(s("t"), div(mul(a[1].clone(), s("x")), sq(c.clone()))));
                let tmp = subst(&subst(&a[0], "x", &s("@x")), "t", &s("@t"));
                subst(&subst(&tmp, "@x", &nx), "@t", &nt)
            }
            "schwarzschild" => {
                need(1)?;
                div(mul_v(vec![int(2), g.clone(), a[0].clone()]), sq(c.clone()))
            }
            "gravdil" => {
                // proper time at radius r from mass M, for coordinate time t
                need(3)?;
                mul(a[0].clone(), sqrt(sub(int(1), div(mul_v(vec![int(2), g.clone(), a[1].clone()]), mul(a[2].clone(), sq(c.clone()))))))
            }
            _ => return Err(String::new()),
        })
    })())
    .filter(|r| !matches!(r, Err(e) if e.is_empty()))
}

/// Total time derivative of f(q, q', t): coordinates and velocities depend on t.
fn ddt(f: &Expr, qs: &[String], t: &str) -> Expr {
    let mut terms = vec![diff(f, t)];
    for q in qs {
        let (qd, qdd) = (format!("{q}'"), format!("{q}''"));
        terms.push(mul(diff(f, q), s(&qd)));
        terms.push(mul(diff(f, &qd), s(&qdd)));
    }
    add_v(terms)
}

fn clean(e: Expr) -> Expr {
    let ex = expand(&e);
    let sm = algebra::simplify(&ex);
    if sm.size() <= ex.size() { sm } else { ex }
}

/// Euler-Lagrange equations d/dt(dL/dq') - dL/dq = 0, one per coordinate, solved for q'' when linear in it.
pub fn lagrange(l: &Expr, qs: &[String], t: &str) -> Vec<Expr> {
    qs.iter()
        .map(|q| {
            let (qd, qdd) = (format!("{q}'"), format!("{q}''"));
            let eom = clean(sub(ddt(&diff(l, &qd), qs, t), diff(l, q)));
            let others = qs.iter().filter(|o| *o != q).any(|o| eom.has_sym(&format!("{o}''")));
            match coeffs_expr(&eom, &qdd) {
                Some(c) if c.len() == 2 && !others => Expr::Rel(Rel::Eq, Box::new(s(&qdd)), Box::new(clean(div(neg(c[0].clone()), c[1].clone())))),
                _ => Expr::Rel(Rel::Eq, Box::new(eom), Box::new(int(0))),
            }
        })
        .collect()
}

pub fn momentum_name(q: &str, single: bool) -> String {
    if single { "p".to_string() } else { format!("p_{q}") }
}

/// Legendre transform: H(q, p) = sum p*q' - L with the velocities eliminated.
pub fn hamiltonian(l: &Expr, qs: &[String]) -> Result<Expr, String> {
    let single = qs.len() == 1;
    let vels: Vec<String> = qs.iter().map(|q| format!("{q}'")).collect();
    let eqs: Vec<Expr> = qs.iter().zip(&vels).map(|(q, v)| sub(s(&momentum_name(q, single)), diff(l, v))).collect();
    let (sols, _) = solve::solve_system(&eqs, &vels, false).map_err(|_| "cannot solve the momenta for the velocities".to_string())?;
    let row = sols.first().ok_or("the Lagrangian is singular: momenta do not determine the velocities")?;
    let mut h = sub(add_v(qs.iter().zip(&vels).map(|(q, v)| mul(s(&momentum_name(q, single)), s(v))).collect()), l.clone());
    for (v, val) in vels.iter().zip(row) {
        h = subst(&h, v, val);
    }
    Ok(clean(h))
}

/// Hamilton's equations: q' = dH/dp, p' = -dH/dq.
pub fn hamilton(h: &Expr, qs: &[String], ps: &[String]) -> Vec<Expr> {
    let mut out = Vec::new();
    for (q, p) in qs.iter().zip(ps) {
        out.push(Expr::Rel(Rel::Eq, Box::new(s(&format!("{q}'"))), Box::new(clean(diff(h, p)))));
        out.push(Expr::Rel(Rel::Eq, Box::new(s(&format!("{p}'"))), Box::new(clean(neg(diff(h, q))))));
    }
    out
}

pub struct FieldEq {
    pub name: &'static str,
    pub title: &'static str,
    pub equation: &'static str,
    pub frame: &'static str,
    pub symbols: &'static str,
    pub residual: fn(&Expr) -> Expr,
    pub solutions: fn() -> Vec<(&'static str, Expr)>,
}

fn d2(u: &Expr, x: &str) -> Expr {
    diff(&diff(u, x), x)
}
fn lap3(u: &Expr) -> Expr {
    add_v(vec![d2(u, "x"), d2(u, "y"), d2(u, "z")])
}
fn radius() -> Expr {
    sqrt(add_v(vec![sq(s("x")), sq(s("y")), sq(s("z"))]))
}
fn uf(name: &str, arg: Expr) -> Expr {
    Expr::Func(name.to_string(), vec![arg])
}
fn ct() -> Expr {
    mul(s("c"), s("t"))
}
fn dalembert() -> Expr {
    add(uf("F", sub(s("x"), ct())), uf("G", add(s("x"), ct())))
}
fn mode_arg() -> Expr {
    div(mul_v(vec![s("n"), PI, s("x")]), s("L"))
}

pub const FIELD_EQS: &[FieldEq] = &[
    FieldEq {
        name: "wave",
        title: "Classical wave equation (vibrating string)",
        equation: "d2u/dt2 = c^2 * d2u/dx2, with c^2 = T/mu",
        frame: "the frame in which the string is at rest",
        symbols: "u(x, t); c wave speed; modes use A, n, L",
        residual: |u| sub(d2(u, "t"), mul(sq(s("c")), d2(u, "x"))),
        solutions: || vec![("d'Alembert travelling waves", dalembert()), ("standing mode of a string of length L", mul_v(vec![s("A"), f1("sin", mode_arg()), f1("cos", div(mul_v(vec![s("n"), PI, s("c"), s("t")]), s("L")))]))],
    },
    FieldEq {
        name: "hooke",
        title: "Ideal spring (Hooke's law)",
        equation: "m * d2u/dt2 = -k*(u - x0)",
        frame: "the frame at rest with respect to the equilibrium position",
        symbols: "u(t) position; m mass; k spring constant; x0 equilibrium; A, phi",
        residual: |u| add(mul(s("m"), d2(u, "t")), mul(s("k"), sub(u.clone(), s("x0")))),
        solutions: || vec![("harmonic oscillation", add(s("x0"), mul(s("A"), f1("cos", add(mul(sqrt(div(s("k"), s("m"))), s("t")), s("phi"))))))],
    },
    FieldEq {
        name: "kepler",
        title: "Kepler problem (inverse-square orbit)",
        equation: "d2r/dt2 + (mu/r^3)*r = 0; orbit form: d2(1/r)/dtheta2 + 1/r = mu/h^2",
        frame: "the frame in which the source of the potential is at rest",
        symbols: "u = r(theta); mu gravitational parameter; h specific angular momentum; ecc eccentricity",
        residual: |u| {
            let w = pow(u.clone(), int(-1));
            sub(add(d2(&w, "theta"), w), div(s("mu"), sq(s("h"))))
        },
        solutions: || vec![("conic section", div(div(sq(s("h")), s("mu")), add(int(1), mul(s("ecc"), f1("cos", s("theta"))))))],
    },
    FieldEq {
        name: "sound",
        title: "Acoustic wave equation",
        equation: "d2u/dt2 = c^2 * (d2u/dx2 + d2u/dy2 + d2u/dz2)",
        frame: "the frame in which the medium (the air) is at rest",
        symbols: "u(x, y, z, t) pressure; c speed of sound; A, k",
        residual: |u| sub(d2(u, "t"), mul(sq(s("c")), lap3(u))),
        solutions: || {
            vec![
                ("plane wave", mul(s("A"), f1("cos", mul(s("k"), sub(s("x"), ct()))))),
                ("spherical wave", div(mul(s("A"), f1("sin", mul(s("k"), sub(radius(), ct())))), radius())),
            ]
        },
    },
    FieldEq {
        name: "heat",
        title: "Heat equation",
        equation: "du/dt = alpha * d2u/dx2",
        frame: "the frame in which the heated material is at rest",
        symbols: "u(x, t) temperature; alpha diffusivity; k",
        residual: |u| sub(diff(u, "t"), mul(s("alpha"), d2(u, "x"))),
        solutions: || {
            vec![
                ("decaying mode", mul(pow(E, neg(mul_v(vec![s("alpha"), sq(s("k")), s("t")]))), f1("sin", mul(s("k"), s("x"))))),
                ("heat kernel", div(pow(E, neg(div(sq(s("x")), mul_v(vec![int(4), s("alpha"), s("t")])))), sqrt(mul_v(vec![int(4), PI, s("alpha"), s("t")])))),
            ]
        },
    },
    FieldEq {
        name: "schrodinger",
        title: "Schrodinger equation (free particle, particle in a box)",
        equation: "i*hbar * du/dt = -(hbar^2/(2*m)) * d2u/dx2",
        frame: "the frame in which the standing waves are stationary (the box is at rest)",
        symbols: "u(x, t) wave function; hbar; m mass; n, L, k",
        residual: |u| add(mul_v(vec![I, s("hbar"), diff(u, "t")]), mul(div(sq(s("hbar")), mul(int(2), s("m"))), d2(u, "x"))),
        solutions: || {
            let energy_t = div(mul_v(vec![sq(s("n")), sq(PI), s("hbar"), s("t")]), mul_v(vec![int(2), s("m"), sq(s("L"))]));
            vec![
                ("stationary state of a box of width L", mul_v(vec![sqrt(div(int(2), s("L"))), f1("sin", mode_arg()), pow(E, neg(mul(I, energy_t)))])),
                ("free plane wave", pow(E, mul(I, sub(mul(s("k"), s("x")), div(mul_v(vec![s("hbar"), sq(s("k")), s("t")]), mul(int(2), s("m"))))))),
            ]
        },
    },
    FieldEq {
        name: "em",
        title: "Electromagnetic wave equation (light in vacuum)",
        equation: "d2u/dt2 = c^2 * d2u/dx2",
        frame: "Lorentz invariant, not Galilean invariant (historically: the frame in which the ether is at rest)",
        symbols: "u(x, t) field component; c speed of light; E0, k",
        residual: |u| sub(d2(u, "t"), mul(sq(s("c")), d2(u, "x"))),
        solutions: || vec![("plane wave", mul(s("E0"), f1("cos", mul(s("k"), sub(s("x"), ct()))))), ("d'Alembert travelling waves", dalembert())],
    },
    FieldEq {
        name: "laplace",
        title: "Laplace equation",
        equation: "d2u/dx2 + d2u/dy2 + d2u/dz2 = 0",
        frame: "static field",
        symbols: "u(x, y, z); k",
        residual: lap3,
        solutions: || vec![("point potential 1/r", pow(radius(), int(-1))), ("separated solution", mul(pow(E, mul(s("k"), s("x"))), f1("sin", mul(s("k"), s("y"))))), ("harmonic polynomial", sub(sq(s("x")), sq(s("y"))))],
    },
    FieldEq {
        name: "poisson",
        title: "Poisson equation (uniform source)",
        equation: "d2u/dx2 + d2u/dy2 + d2u/dz2 = -rho/eps0",
        frame: "static field",
        symbols: "u(x, y, z) potential; rho charge density; eps0",
        residual: |u| add(lap3(u), div(s("rho"), s("eps0"))),
        solutions: || vec![("uniform charge density", neg(div(mul(s("rho"), add_v(vec![sq(s("x")), sq(s("y")), sq(s("z"))])), mul(int(6), s("eps0")))))],
    },
    FieldEq {
        name: "diffusion",
        title: "Diffusion equation (three dimensions)",
        equation: "du/dt = D * (d2u/dx2 + d2u/dy2 + d2u/dz2)",
        frame: "the frame in which the medium is at rest",
        symbols: "u(x, y, z, t) concentration; D diffusion constant",
        residual: |u| sub(diff(u, "t"), mul(s("D"), lap3(u))),
        solutions: || {
            let dt4 = mul_v(vec![int(4), s("D"), s("t")]);
            vec![("point-source kernel", mul(pow(mul(PI, dt4.clone()), rat(-3, 2)), pow(E, neg(div(add_v(vec![sq(s("x")), sq(s("y")), sq(s("z"))]), dt4)))))]
        },
    },
    FieldEq {
        name: "kleingordon",
        title: "Klein-Gordon equation",
        equation: "(1/c^2) * d2u/dt2 - d2u/dx2 + (m*c/hbar)^2 * u = 0",
        frame: "Lorentz invariant",
        symbols: "u(x, t); m mass; c; hbar; k",
        residual: |u| add_v(vec![div(d2(u, "t"), sq(s("c"))), neg(d2(u, "x")), mul(sq(div(mul(s("m"), s("c")), s("hbar"))), u.clone())]),
        solutions: || {
            let omega = mul(s("c"), sqrt(add(sq(s("k")), sq(div(mul(s("m"), s("c")), s("hbar"))))));
            vec![("plane wave", pow(E, mul(I, sub(mul(s("k"), s("x")), mul(omega, s("t"))))))]
        },
    },
];

pub fn field(name: &str) -> Result<&'static FieldEq, String> {
    FIELD_EQS.iter().find(|f| f.name == name).ok_or_else(|| format!("unknown field equation '{name}'; fieldeq() lists them"))
}

pub fn field_list() -> String {
    FIELD_EQS.iter().map(|f| format!("{:12} {}: {}", f.name, f.title, f.equation)).collect::<Vec<_>>().join("\n")
}

pub fn field_info(f: &FieldEq, show: &dyn Fn(&Expr) -> String) -> String {
    let mut out = format!("{}\n{}\nholds in: {}\nsymbols: {}", f.title, f.equation, f.frame, f.symbols);
    for (i, (label, e)) in (f.solutions)().iter().enumerate() {
        out.push_str(&format!("\nsolution {} ({}): {}", i + 1, label, show(e)));
    }
    out
}

/// Substitute u into the equation. Ok(None) means it is satisfied; Ok(Some(r)) returns the nonzero residual.
pub fn field_check(f: &FieldEq, u: &Expr) -> Option<Expr> {
    let r = (f.residual)(u);
    let r = expand(&r);
    if r.is_zero() {
        return None;
    }
    let sm = algebra::simplify(&r);
    if sm.is_zero() {
        return None;
    }
    // numeric test at sample points (parameters positive)
    let syms = sm.free_syms();
    let mut good = 0;
    for k in 0..5 {
        let env: Vec<(&str, C)> = syms.iter().enumerate().map(|(i, n)| (n.as_str(), C::new(0.57 + 0.23 * i as f64 + 0.19 * k as f64, 0.0))).collect();
        let Some(v) = evalc(&r, &env) else { return Some(sm) };
        let terms: Vec<&Expr> = if let Expr::Add(t) = &r { t.iter().collect() } else { vec![&r] };
        let scale: f64 = terms.iter().filter_map(|t| evalc(t, &env)).map(|t| t.norm()).filter(|t| t.is_finite()).sum();
        if !v.norm().is_finite() {
            continue;
        }
        if v.norm() > 1e-8 * (1.0 + scale) {
            return Some(sm);
        }
        good += 1;
    }
    if good > 0 { None } else { Some(sm) }
}

pub const SPH_NAMES: &[&str] = &["sphsss", "sphsas", "sphasa", "sphaaa", "sphexcess", "spharea", "sphdist", "sphbearing", "navtri"];

fn named(pairs: Vec<(&str, Expr)>) -> Expr {
    Expr::List(pairs.into_iter().map(|(n, v)| Expr::Rel(Rel::Eq, Box::new(s(n)), Box::new(v))).collect())
}

/// Spherical trigonometry. Angles and arc sides are in the current angle mode; `k` converts them to radians.
pub fn spherical(name: &str, a: &[Expr], k: &Expr) -> Option<Result<Expr, String>> {
    let need = |n: usize| if a.len() < n { Err(format!("{name} needs {n} argument(s)")) } else { Ok(()) };
    let rad = |i: usize| mul(a[i].clone(), k.clone());
    let out = |e: Expr| {
        let v = div(e, k.clone());
        if !v.is_constant() {
            return v;
        }
        // exact only when it comes out short (table angles); otherwise a decimal
        let ex = expand(&v);
        if ex.size() <= 8 && !ex.has(&|n| matches!(n, Expr::Func(..) | Expr::Const(Const::I))) { ex } else { to_float(&v) }
    };
    let (sin, cos) = (|e: Expr| f1("sin", e), |e: Expr| f1("cos", e));
    let acos = |e: Expr| f1("acos", e);
    // angle opposite side p in a triangle with sides p, q, r (law of cosines for sides)
    let angle = |p: &Expr, q: &Expr, r: &Expr| acos(div(sub(cos(p.clone()), mul(cos(q.clone()), cos(r.clone()))), mul(sin(q.clone()), sin(r.clone()))));
    // side opposite angle P in a triangle with angles P, Q, R (law of cosines for angles)
    let side = |p: &Expr, q: &Expr, r: &Expr| acos(div(add(cos(p.clone()), mul(cos(q.clone()), cos(r.clone()))), mul(sin(q.clone()), sin(r.clone()))));
    Some((|| -> Result<Expr, String> {
        Ok(match name {
            "sphsss" => {
                need(3)?;
                let (x, y, z) = (rad(0), rad(1), rad(2));
                named(vec![("A", out(angle(&x, &y, &z))), ("B", out(angle(&y, &x, &z))), ("C", out(angle(&z, &x, &y)))])
            }
            "sphsas" => {
                // side a, included angle C, side b
                need(3)?;
                let (x, cc, y) = (rad(0), rad(1), rad(2));
                let z = acos(add(mul(cos(x.clone()), cos(y.clone())), mul_v(vec![sin(x.clone()), sin(y.clone()), cos(cc)])));
                named(vec![("c", out(z.clone())), ("A", out(angle(&x, &y, &z))), ("B", out(angle(&y, &x, &z)))])
            }
            "sphasa" => {
                // angle A, included side c, angle B
                need(3)?;
                let (p, z, q) = (rad(0), rad(1), rad(2));
                let r = acos(add(neg(mul(cos(p.clone()), cos(q.clone()))), mul_v(vec![sin(p.clone()), sin(q.clone()), cos(z)])));
                named(vec![("C", out(r.clone())), ("a", out(side(&p, &q, &r))), ("b", out(side(&q, &p, &r)))])
            }
            "sphaaa" => {
                need(3)?;
                let (p, q, r) = (rad(0), rad(1), rad(2));
                named(vec![("a", out(side(&p, &q, &r))), ("b", out(side(&q, &p, &r))), ("c", out(side(&r, &p, &q)))])
            }
            "sphexcess" => {
                need(3)?;
                out(sub(add_v(vec![rad(0), rad(1), rad(2)]), PI))
            }
            "spharea" => {
                // angles A, B, C and sphere radius R (default 1)
                need(3)?;
                let r = a.get(3).cloned().unwrap_or(int(1));
                mul(sub(add_v(vec![rad(0), rad(1), rad(2)]), PI), sq(r))
            }
            "sphdist" => {
                // lat1, lon1, lat2, lon2 [, R]: central angle, or arc length when R is given
                need(4)?;
                let (p1, l1, p2, l2) = (rad(0), rad(1), rad(2), rad(3));
                let h = add(sq(sin(mul(rat(1, 2), sub(p2.clone(), p1.clone())))), mul_v(vec![cos(p1), cos(p2), sq(sin(mul(rat(1, 2), sub(l2, l1))))]));
                let ang = mul(int(2), f1("asin", sqrt(h)));
                match a.get(4) {
                    Some(r) => mul(ang, r.clone()),
                    None => out(ang),
                }
            }
            "sphbearing" | "navtri" => {
                let (north, east, extra) = if name == "sphbearing" {
                    need(4)?;
                    let (p1, l1, p2, l2) = (rad(0), rad(1), rad(2), rad(3));
                    let dl = sub(l2, l1);
                    (sub(mul(cos(p1.clone()), sin(p2.clone())), mul_v(vec![sin(p1), cos(p2.clone()), cos(dl.clone())])), mul(sin(dl), cos(p2)), None)
                } else {
                    // latitude, declination, local hour angle -> altitude and azimuth (from north, eastward)
                    need(3)?;
                    let (lat, dec, lha) = (rad(0), rad(1), rad(2));
                    let alt = f1("asin", add(mul(sin(lat.clone()), sin(dec.clone())), mul_v(vec![cos(lat.clone()), cos(dec.clone()), cos(lha.clone())])));
                    (sub(mul(sin(dec.clone()), cos(lat.clone())), mul_v(vec![cos(dec.clone()), sin(lat), cos(lha.clone())])), neg(mul(sin(lha), cos(dec))), Some(alt))
                };
                let mut az = f1("arg", add(north, mul(I, east)));
                if az.is_constant() && crate::numeric::real_const(&az).is_some_and(|v| v < 0.0) {
                    az = add(az, mul(int(2), PI));
                }
                match extra {
                    Some(alt) => named(vec![("alt", out(alt)), ("az", out(az))]),
                    None => out(az),
                }
            }
            _ => return Err(String::new()),
        })
    })())
    .filter(|r| !matches!(r, Err(e) if e.is_empty()))
}
