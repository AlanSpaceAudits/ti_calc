//! Geometric figures for results: curves, shaded areas, vector and circle diagrams, spheres, phase portraits.
//! A figure is a list of plain shapes; a front end only has to draw them.
use crate::eval::Calc;
use crate::expr::*;
use crate::numeric::evalc;
use crate::simplify::*;
use num_complex::Complex64 as C;
use std::collections::HashMap;
use std::f64::consts::{PI as PI_F, TAU};

pub enum Shape {
    Line { pts: Vec<[f64; 2]>, color: u8, dashed: bool, name: String },
    Arrow { from: [f64; 2], to: [f64; 2], color: u8, name: String },
    Fill { pts: Vec<[f64; 2]>, color: u8 },
    Dot { at: [f64; 2], color: u8, name: String },
    Text { at: [f64; 2], text: String },
}

#[derive(Default)]
pub struct Figure {
    pub title: String,
    pub caption: String,
    pub shapes: Vec<Shape>,
    /// equal scale on both axes (geometry)
    pub equal: bool,
    /// square frame with ordinary axes (phase portraits, curves in the plane)
    pub square: bool,
    /// xmin, xmax, ymin, ymax
    pub bounds: Option<[f64; 4]>,
    /// expression in x that the Graph tab can take over
    pub graph_src: Option<String>,
}

const GREY: u8 = 6;

impl Figure {
    fn new(title: impl Into<String>) -> Figure {
        Figure { title: title.into(), ..Default::default() }
    }
    fn line(&mut self, pts: Vec<[f64; 2]>, color: u8, name: &str) {
        self.shapes.push(Shape::Line { pts, color, dashed: false, name: name.to_string() });
    }
    fn dashed(&mut self, pts: Vec<[f64; 2]>, color: u8, name: &str) {
        self.shapes.push(Shape::Line { pts, color, dashed: true, name: name.to_string() });
    }
    fn arrow(&mut self, from: [f64; 2], to: [f64; 2], color: u8, name: &str) {
        self.shapes.push(Shape::Arrow { from, to, color, name: name.to_string() });
    }
    fn dot(&mut self, at: [f64; 2], color: u8, name: &str) {
        self.shapes.push(Shape::Dot { at, color, name: name.to_string() });
    }
    fn text(&mut self, at: [f64; 2], text: impl Into<String>) {
        self.shapes.push(Shape::Text { at, text: text.into() });
    }
    fn grow(&mut self, b: [f64; 4]) {
        self.bounds = Some(match self.bounds {
            None => b,
            Some(o) => [o[0].min(b[0]), o[1].max(b[1]), o[2].min(b[2]), o[3].max(b[3])],
        });
    }
    /// Sampled curve y = f(x) with breaks at undefined points and across asymptotes.
    fn curve(&mut self, f: &dyn Fn(f64) -> f64, lo: f64, hi: f64, color: u8, name: &str, dashed: bool) {
        let n = 500;
        let pts: Vec<[f64; 2]> = (0..=n).map(|i| lo + (hi - lo) * i as f64 / n as f64).map(|x| [x, f(x)]).collect();
        let mut ys: Vec<f64> = pts.iter().map(|p| p[1]).filter(|y| y.is_finite()).collect();
        if ys.is_empty() {
            return;
        }
        ys.sort_by(f64::total_cmp);
        let (p_lo, p_hi) = (ys[ys.len() * 3 / 100], ys[(ys.len() * 97 / 100).min(ys.len() - 1)]);
        let span = (p_hi - p_lo).max(1e-9);
        let (y0, y1) = (p_lo - 0.25 * span, p_hi + 0.25 * span);
        let mut segs: Vec<Vec<[f64; 2]>> = vec![vec![]];
        for p in pts {
            let bad = !p[1].is_finite();
            let jump = segs.last().unwrap().last().is_some_and(|q: &[f64; 2]| (q[1] - p[1]).abs() > 3.0 * span);
            if bad || jump {
                segs.push(vec![]);
            }
            if !bad {
                segs.last_mut().unwrap().push(p);
            }
        }
        for s in segs.into_iter().filter(|s| s.len() > 1) {
            self.shapes.push(Shape::Line { pts: s, color, dashed, name: name.to_string() });
        }
        self.grow([lo, hi, y0.min(0.0 - 0.05 * span).max(y0 - span), y1]);
    }
}

fn visible_syms(e: &Expr) -> Vec<String> {
    e.free_syms().into_iter().filter(|s| !s.starts_with('@')).collect()
}

/// Replace every symbol except `keep` by a sample value so the expression can be drawn.
fn with_defaults(e: &Expr, keep: &[&str]) -> (Expr, String) {
    let mut out = e.clone();
    let mut used = Vec::new();
    for s in visible_syms(e) {
        if keep.contains(&s.as_str()) {
            continue;
        }
        let v = match s.as_str() {
            "ecc" => 0.5,
            "L" => 3.0,
            _ => 1.0,
        };
        out = subst(&out, &s, &Expr::Float(v));
        used.push(format!("{s} = {v}"));
    }
    (out, if used.is_empty() { String::new() } else { format!("drawn with {}", used.join(", ")) })
}

fn real_fn<'a>(e: &'a Expr, var: &'a str) -> impl Fn(f64) -> f64 + 'a {
    move |x| match evalc(e, &[(var, C::new(x, 0.0))]) {
        Some(z) if z.im.abs() <= 1e-9 * (1.0 + z.re.abs()) => z.re,
        _ => f64::NAN,
    }
}

fn part_fn<'a>(e: &'a Expr, var: &'a str, imag: bool) -> impl Fn(f64) -> f64 + 'a {
    move |x| evalc(e, &[(var, C::new(x, 0.0))]).map_or(f64::NAN, |z| if imag { z.im } else { z.re })
}

fn circle_pts(c: [f64; 2], r: f64, a0: f64, a1: f64) -> Vec<[f64; 2]> {
    let n = 96;
    (0..=n).map(|i| a0 + (a1 - a0) * i as f64 / n as f64).map(|t| [c[0] + r * t.cos(), c[1] + r * t.sin()]).collect()
}

fn on_circle(r: f64, t: f64) -> [f64; 2] {
    [r * t.cos(), r * t.sin()]
}

fn fmt(v: f64) -> String {
    let s = format!("{v:.4}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" { "0".into() } else { s.to_string() }
}

/// Zero-level contour of f on a grid, stitched into polylines.
fn contour(f: &dyn Fn(f64, f64) -> f64, xr: [f64; 2], yr: [f64; 2], n: usize) -> Vec<Vec<[f64; 2]>> {
    let gx = |i: usize| xr[0] + (xr[1] - xr[0]) * i as f64 / n as f64;
    let gy = |j: usize| yr[0] + (yr[1] - yr[0]) * j as f64 / n as f64;
    let vals: Vec<Vec<f64>> = (0..=n).map(|i| (0..=n).map(|j| f(gx(i), gy(j))).collect()).collect();
    // crossing point on a grid edge, keyed by (i, j, direction)
    let cross = |i: usize, j: usize, horiz: bool| -> Option<[f64; 2]> {
        let (a, b) = if horiz { (vals[i][j], vals[i + 1][j]) } else { (vals[i][j], vals[i][j + 1]) };
        if !a.is_finite() || !b.is_finite() || (a < 0.0) == (b < 0.0) || (a - b).abs() > 1e6 {
            return None;
        }
        let t = a / (a - b);
        Some(if horiz { [gx(i) + t * (gx(i + 1) - gx(i)), gy(j)] } else { [gx(i), gy(j) + t * (gy(j + 1) - gy(j))] })
    };
    type Key = (usize, usize, bool);
    let mut segs: Vec<(Key, Key)> = Vec::new();
    let mut pts: HashMap<Key, [f64; 2]> = HashMap::new();
    for i in 0..n {
        for j in 0..n {
            let edges: [Key; 4] = [(i, j, true), (i, j + 1, true), (i, j, false), (i + 1, j, false)];
            let hit: Vec<Key> = edges.into_iter().filter(|k| cross(k.0, k.1, k.2).map(|p| pts.insert(*k, p)).is_some()).collect();
            match hit.len() {
                2 => segs.push((hit[0], hit[1])),
                4 => {
                    segs.push((hit[0], hit[2]));
                    segs.push((hit[1], hit[3]));
                }
                _ => {}
            }
        }
    }
    let mut adj: HashMap<Key, Vec<usize>> = HashMap::new();
    for (i, (a, b)) in segs.iter().enumerate() {
        adj.entry(*a).or_default().push(i);
        adj.entry(*b).or_default().push(i);
    }
    let mut used = vec![false; segs.len()];
    let mut out = Vec::new();
    for start in 0..segs.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut chain = std::collections::VecDeque::from([segs[start].0, segs[start].1]);
        for front in [false, true] {
            loop {
                let end = if front { *chain.front().unwrap() } else { *chain.back().unwrap() };
                let Some(&next) = adj[&end].iter().find(|s| !used[**s]) else { break };
                used[next] = true;
                let other = if segs[next].0 == end { segs[next].1 } else { segs[next].0 };
                if front { chain.push_front(other) } else { chain.push_back(other) }
            }
        }
        out.push(chain.into_iter().map(|k| pts[&k]).collect());
    }
    out
}

type V3 = [f64; 3];
fn dot3(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross3(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn unit(a: V3) -> V3 {
    let n = dot3(a, a).sqrt().max(1e-12);
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Orthographic view of the unit sphere, looking along -w.
struct Sphere {
    u: V3,
    v: V3,
    w: V3,
}

impl Sphere {
    fn facing(toward: V3) -> Sphere {
        let w = unit(toward);
        let up = if w[2].abs() > 0.95 { [1.0, 0.0, 0.0] } else { [0.0, 0.0, 1.0] };
        let u = unit(cross3(up, w));
        Sphere { u, v: cross3(w, u), w }
    }
    fn proj(&self, p: V3) -> ([f64; 2], bool) {
        ([dot3(p, self.u), dot3(p, self.v)], dot3(p, self.w) >= -1e-9)
    }
    /// Draw a path on the sphere: solid on the near side, dashed on the far side.
    fn path(&self, fig: &mut Figure, pts: &[V3], color: u8, name: &str) {
        self.path_opt(fig, pts, color, name, true)
    }
    fn path_opt(&self, fig: &mut Figure, pts: &[V3], color: u8, name: &str, back: bool) {
        let mut run: Vec<[f64; 2]> = Vec::new();
        let mut run_front = true;
        for p in pts {
            let (q, front) = self.proj(*p);
            if front != run_front && !run.is_empty() {
                run.push(q);
                let done = std::mem::take(&mut run);
                if done.len() > 1 {
                    if run_front {
                        fig.line(done, color, name)
                    } else if back {
                        fig.dashed(done, color, name)
                    }
                }
            }
            run_front = front;
            run.push(q);
        }
        if run.len() > 1 {
            if run_front {
                fig.line(run, color, name)
            } else if back {
                fig.dashed(run, color, name)
            }
        }
    }
    fn arc(&self, fig: &mut Figure, a: V3, b: V3, color: u8, name: &str) {
        let omega = dot3(a, b).clamp(-1.0, 1.0).acos();
        let pts: Vec<V3> = (0..=64)
            .map(|i| {
                let t = i as f64 / 64.0;
                if omega < 1e-9 {
                    return a;
                }
                let (ka, kb) = (((1.0 - t) * omega).sin() / omega.sin(), (t * omega).sin() / omega.sin());
                [ka * a[0] + kb * b[0], ka * a[1] + kb * b[1], ka * a[2] + kb * b[2]]
            })
            .collect();
        self.path(fig, &pts, color, name);
    }
    /// Circle of points at angular distance `rad` from the axis.
    fn ring(&self, fig: &mut Figure, axis: V3, rad: f64, color: u8, name: &str) {
        self.ring_opt(fig, axis, rad, color, name, true)
    }
    /// Graticule lines: near side only.
    fn grid_ring(&self, fig: &mut Figure, axis: V3, rad: f64) {
        self.ring_opt(fig, axis, rad, GREY, "", false)
    }
    fn ring_opt(&self, fig: &mut Figure, axis: V3, rad: f64, color: u8, name: &str, back: bool) {
        let ax = unit(axis);
        let e1 = unit(cross3(if ax[2].abs() > 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 0.0, 1.0] }, ax));
        let e2 = cross3(ax, e1);
        let pts: Vec<V3> = (0..=120)
            .map(|i| i as f64 / 120.0 * TAU)
            .map(|t| {
                let (c, s) = (rad.sin() * t.cos(), rad.sin() * t.sin());
                [rad.cos() * ax[0] + c * e1[0] + s * e2[0], rad.cos() * ax[1] + c * e1[1] + s * e2[1], rad.cos() * ax[2] + c * e1[2] + s * e2[2]]
            })
            .collect();
        self.path_opt(fig, &pts, color, name, back);
    }
    fn outline(&self, fig: &mut Figure) {
        fig.line(circle_pts([0.0, 0.0], 1.0, 0.0, TAU), GREY, "");
        fig.equal = true;
        fig.bounds = Some([-1.25, 1.25, -1.25, 1.25]);
    }
    fn label(&self, fig: &mut Figure, p: V3, text: &str, color: u8) {
        let (q, _) = self.proj(p);
        fig.dot(q, color, "");
        fig.text([q[0] * 1.12 + 0.02, q[1] * 1.12 + 0.02], text);
    }
}

/// Unit vector from latitude and longitude in degrees.
fn latlon(lat: f64, lon: f64) -> V3 {
    let (p, l) = (lat.to_radians(), lon.to_radians());
    [p.cos() * l.cos(), p.cos() * l.sin(), p.sin()]
}

fn list_value(v: &Expr, name: &str) -> Option<f64> {
    let Expr::List(items) = v else { return None };
    items.iter().find_map(|it| match it {
        Expr::Rel(Rel::Eq, l, r) if **l == sym(name) => crate::numeric::real_const(&to_float(r)),
        _ => None,
    })
}

fn solutions_of(v: &Expr, var: &str) -> Vec<f64> {
    let items: Vec<&Expr> = match v {
        Expr::Func(n, a) if n == "or" => a.iter().collect(),
        Expr::List(a) => a.iter().collect(),
        e => vec![e],
    };
    items
        .into_iter()
        .filter_map(|it| match it {
            Expr::Rel(_, l, r) if **l == sym(var) && r.is_constant() => crate::numeric::real_const(r),
            e if e.is_constant() => crate::numeric::real_const(e),
            _ => None,
        })
        .collect()
}

fn span_around(points: &[f64]) -> (f64, f64) {
    if points.is_empty() {
        return (-10.0, 10.0);
    }
    let (lo, hi) = (points.iter().cloned().fold(f64::INFINITY, f64::min), points.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    let pad = ((hi - lo) * 0.6).max(3.0);
    (lo - pad, hi + pad)
}

fn pulse(e: &Expr) -> Expr {
    // stand-ins for the arbitrary functions F and G of a travelling-wave solution
    match e {
        Expr::Func(n, a) if (n == "F" || n == "G") && a.len() == 1 => {
            let p = pow(E, neg(pow(pulse(&a[0]), int(2))));
            if n == "F" { p } else { mul(rat(1, 2), p) }
        }
        Expr::Add(v) => add_v(v.iter().map(pulse).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(pulse).collect()),
        Expr::Pow(b, x) => pow(pulse(b), pulse(x)),
        Expr::Func(n, v) => func(n, v.iter().map(pulse).collect()),
        e => e.clone(),
    }
}

fn parallelogram(u: [f64; 2], v: [f64; 2], labels: [&str; 3], title: &str) -> Figure {
    let mut fig = Figure::new(title);
    let s = [u[0] + v[0], u[1] + v[1]];
    fig.dashed(vec![u, s], GREY, "");
    fig.dashed(vec![v, s], GREY, "");
    fig.arrow([0.0, 0.0], u, 0, labels[0]);
    fig.arrow([0.0, 0.0], v, 1, labels[1]);
    fig.arrow([0.0, 0.0], s, 2, labels[2]);
    let xs = [0.0, u[0], v[0], s[0]];
    let ys = [0.0, u[1], v[1], s[1]];
    let (x0, x1) = (xs.iter().cloned().fold(f64::INFINITY, f64::min), xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    let (y0, y1) = (ys.iter().cloned().fold(f64::INFINITY, f64::min), ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    let pad = 0.15 * (x1 - x0).max(y1 - y0).max(1e-9);
    fig.bounds = Some([x0 - pad, x1 + pad, y0 - pad, y1 + pad]);
    fig.equal = true;
    fig
}

fn vec2(e: &Expr) -> Option<[f64; 2]> {
    let items: Vec<&Expr> = match e {
        Expr::List(v) => v.iter().collect(),
        Expr::Matrix(m) => m.iter().flatten().collect(),
        _ => return None,
    };
    if items.len() != 2 {
        return None;
    }
    Some([crate::numeric::real_const(&to_float(items[0]))?, crate::numeric::real_const(&to_float(items[1]))?])
}

fn spacetime(beta: f64, lorentz: bool) -> Figure {
    let mut fig = Figure::new(if lorentz { "Lorentz boost: both axes tilt toward the light line" } else { "Galilean boost: only the time axis tilts" });
    let r = 4.0;
    fig.line(vec![[-r, 0.0], [r, 0.0]], GREY, "x");
    fig.line(vec![[0.0, -r], [0.0, r]], GREY, "c*t");
    fig.dashed(vec![[-r, -r], [r, r]], 3, "light: x = c*t");
    fig.dashed(vec![[-r, r], [r, -r]], 3, "light: x = -c*t");
    // t' axis: the worldline x = v*t of the moving origin
    fig.line(vec![[-r * beta, -r], [r * beta, r]], 0, "t' axis");
    if lorentz {
        fig.line(vec![[-r, -r * beta], [r, r * beta]], 1, "x' axis");
    } else {
        fig.line(vec![[-r, 0.02], [r, 0.02]], 1, "x' axis (same as x)");
    }
    fig.text([r * 0.92, 0.25], "x");
    fig.text([0.2, r * 0.92], "c*t");
    fig.caption = format!("v/c = {}", fmt(beta));
    fig.bounds = Some([-r, r, -r, r]);
    fig.equal = true;
    fig
}

fn circle_figure(kind: &str, ang_deg: f64) -> Figure {
    let a = ang_deg.to_radians();
    let mut fig = Figure::new(match kind {
        "inscribed" => "Inscribed angle = half the central angle on the same arc",
        "tangentchord" => "Tangent-chord angle = half the intercepted arc",
        _ => "Central angle and its chord",
    });
    fig.line(circle_pts([0.0, 0.0], 1.0, 0.0, TAU), GREY, "");
    let o = [0.0, 0.0];
    // the arc AB sits at the bottom, symmetric about the vertical axis
    let (ta, tb) = (-PI_F / 2.0 - a / 2.0, -PI_F / 2.0 + a / 2.0);
    let (pa, pb) = (on_circle(1.0, ta), on_circle(1.0, tb));
    match kind {
        "tangentchord" => {
            // tangent at A, chord AB; the arc AB measures ang_deg
            let dir = [-(ta.sin()), ta.cos()];
            fig.line(vec![[pa[0] - 1.4 * dir[0], pa[1] - 1.4 * dir[1]], [pa[0] + 1.4 * dir[0], pa[1] + 1.4 * dir[1]]], 1, "tangent at A");
            fig.line(vec![pa, pb], 0, "chord AB");
            fig.line(circle_pts(o, 1.0, ta, tb), 2, "intercepted arc");
            fig.dashed(vec![o, pa], GREY, "");
            fig.dashed(vec![o, pb], GREY, "");
            let t0 = dir[1].atan2(dir[0]);
            fig.line(circle_pts(pa, 0.25, t0, t0 + a / 2.0), 3, "");
            fig.caption = format!("arc AB = {}, angle between tangent and chord = {}", fmt(ang_deg), fmt(ang_deg / 2.0));
        }
        "inscribed" => {
            let pp = on_circle(1.0, PI_F / 2.0 + 0.5);
            fig.line(vec![pa, o, pb], 0, "central angle AOB");
            fig.line(vec![pa, pp, pb], 1, "inscribed angle APB");
            fig.line(circle_pts(o, 1.0, ta, tb), 2, "arc AB");
            fig.line(circle_pts(o, 0.22, ta, tb), 0, "");
            let (d1, d2) = ((pa[1] - pp[1]).atan2(pa[0] - pp[0]), (pb[1] - pp[1]).atan2(pb[0] - pp[0]));
            fig.line(circle_pts(pp, 0.3, d1, d2), 1, "");
            fig.dot(pp, 1, "");
            fig.text([pp[0] - 0.05, pp[1] + 0.1], "P");
            fig.caption = format!("central angle = {}, inscribed angle = {}", fmt(ang_deg), fmt(ang_deg / 2.0));
        }
        _ => {
            fig.line(vec![pa, o, pb], 0, "radii");
            fig.line(vec![pa, pb], 1, "chord AB");
            fig.line(circle_pts(o, 0.22, ta, tb), 0, "");
            let mid = [0.0, pa[1]];
            fig.dashed(vec![mid, [0.0, -1.0]], 3, "sagitta");
            fig.caption = format!("central angle = {}; chord = 2*r*sin(angle/2)", fmt(ang_deg));
        }
    }
    fig.dot(o, GREY, "");
    fig.text([0.05, 0.06], "O");
    fig.dot(pa, 0, "");
    fig.dot(pb, 0, "");
    fig.text([pa[0] - 0.1, pa[1] - 0.1], "A");
    fig.text([pb[0] + 0.06, pb[1] - 0.1], "B");
    fig.bounds = Some([-1.5, 1.5, -1.5, 1.5]);
    fig.equal = true;
    fig
}

fn sphere_triangle(a: f64, b: f64, cc: f64) -> Figure {
    // vertex C at the pole; sides b = CA and a = CB with angle C between them (degrees)
    let vc: V3 = [0.0, 0.0, 1.0];
    let va = latlon(90.0 - b, 0.0);
    let vb = latlon(90.0 - a, cc);
    let view = Sphere::facing(unit([va[0] + vb[0] + vc[0] + 0.6, va[1] + vb[1] + vc[1] - 0.2, va[2] + vb[2] + vc[2]]));
    let mut fig = Figure::new("Spherical triangle: sides are arcs of great circles");
    view.outline(&mut fig);
    view.ring(&mut fig, vc, PI_F / 2.0, GREY, "");
    view.arc(&mut fig, vb, vc, 0, "side a");
    view.arc(&mut fig, vc, va, 1, "side b");
    view.arc(&mut fig, va, vb, 2, "side c");
    view.label(&mut fig, va, "A", 0);
    view.label(&mut fig, vb, "B", 1);
    view.label(&mut fig, vc, "C", 2);
    fig
}

/// Build a figure for one evaluated input line, when it has a natural picture.
pub fn build(calc: &Calc, raw: &Expr, value: &Expr) -> Option<Figure> {
    let kdeg = calc.degrees_per_unit();
    if let Expr::Func(name, args) = raw {
        let ev = |i: usize, free: &[&str]| -> Option<Expr> { calc.ev_free(args.get(i)?, free).ok() };
        let num = |i: usize| -> Option<f64> { ev(i, &[]).and_then(|e| crate::numeric::real_const(&to_float(&e))) };
        let var = |i: usize| -> Option<String> {
            match args.get(i)? {
                Expr::Sym(s) => Some(s.clone()),
                _ => None,
            }
        };
        match name.as_str() {
            "solve" | "csolve" | "zeros" | "nsolve" => {
                // system in two unknowns: the curves and their intersections
                if let (Some(Expr::List(eqs)), Some(Expr::List(vs))) = (args.first(), args.get(1)) {
                    if eqs.len() == 2 && vs.len() == 2 {
                        let (vx, vy) = (var_name(&vs[0])?, var_name(&vs[1])?);
                        let mut fig = Figure::new("Each equation is a curve; solutions are the crossing points");
                        let mut pts: Vec<[f64; 2]> = Vec::new();
                        let rows: Vec<&Expr> = match value {
                            Expr::Func(n, a) if n == "or" => a.iter().collect(),
                            e => vec![e],
                        };
                        for row in rows {
                            let items: Vec<&Expr> = match row {
                                Expr::Func(n, a) if n == "and" => a.iter().collect(),
                                e => vec![e],
                            };
                            let get = |nm: &str| items.iter().find_map(|it| match it {
                                Expr::Rel(Rel::Eq, l, r) if **l == sym(nm) => crate::numeric::real_const(&to_float(r)),
                                _ => None,
                            });
                            if let (Some(px), Some(py)) = (get(&vx), get(&vy)) {
                                pts.push([px, py]);
                            }
                        }
                        let (x0, x1) = span_around(&pts.iter().map(|p| p[0]).collect::<Vec<_>>());
                        let (y0, y1) = span_around(&pts.iter().map(|p| p[1]).collect::<Vec<_>>());
                        for (i, eq) in eqs.iter().enumerate() {
                            let e = calc.ev_free(eq, &[&vx, &vy]).ok()?;
                            let f = match &e {
                                Expr::Rel(Rel::Eq, l, r) => sub((**l).clone(), (**r).clone()),
                                e => e.clone(),
                            };
                            if visible_syms(&f).iter().any(|s| *s != vx && *s != vy) {
                                return None;
                            }
                            let g = |x: f64, y: f64| evalc(&f, &[(vx.as_str(), C::new(x, 0.0)), (vy.as_str(), C::new(y, 0.0))]).filter(|z| z.im.abs() < 1e-9).map_or(f64::NAN, |z| z.re);
                            for seg in contour(&g, [x0, x1], [y0, y1], 140) {
                                fig.line(seg, i as u8, &calc.show(&e));
                            }
                        }
                        for p in pts {
                            fig.dot(p, 2, &format!("({}, {})", fmt(p[0]), fmt(p[1])));
                        }
                        fig.bounds = Some([x0, x1, y0, y1]);
                        fig.square = true;
                        return Some(fig);
                    }
                }
                let x = var(1).or_else(|| visible_syms(&ev(0, &[])?).first().cloned())?;
                let e = ev(0, &[&x])?;
                let (l, r) = match &e {
                    Expr::Rel(_, l, r) => ((**l).clone(), (**r).clone()),
                    e => (e.clone(), int(0)),
                };
                if visible_syms(&l).iter().chain(&visible_syms(&r)).any(|s| *s != x) {
                    return None;
                }
                let roots = solutions_of(value, &x);
                let (lo, hi) = span_around(&roots);
                let mut fig = Figure::new("Both sides of the equation; solutions are where they meet");
                fig.curve(&real_fn(&l, &x), lo, hi, 0, &calc.show(&l), false);
                fig.curve(&real_fn(&r, &x), lo, hi, 1, &calc.show(&r), false);
                for rt in roots {
                    let y = real_fn(&r, &x)(rt);
                    if y.is_finite() {
                        fig.dot([rt, y], 2, &format!("{x} = {}", fmt(rt)));
                    }
                }
                return Some(fig);
            }
            "diff" => {
                let x = var(1)?;
                let f = ev(0, &[&x])?;
                if visible_syms(&f).iter().any(|s| *s != x) || visible_syms(value).iter().any(|s| *s != x) {
                    return None;
                }
                let mut fig = Figure::new("A function and its derivative (the slope of the first curve)");
                fig.curve(&real_fn(&f, &x), -6.0, 6.0, 0, &calc.show(&f), false);
                fig.curve(&real_fn(value, &x), -6.0, 6.0, 1, &format!("derivative: {}", calc.show(value)), false);
                return Some(fig);
            }
            "integrate" | "nint" => {
                let x = var(1)?;
                let f = ev(0, &[&x])?;
                if visible_syms(&f).iter().any(|s| *s != x) {
                    return None;
                }
                let g = real_fn(&f, &x);
                if let (Some(a), Some(b)) = (num(2), num(3)) {
                    let (a, b) = (a.max(-1e3), b.min(1e3));
                    let pad = ((b - a) * 0.25).max(0.5);
                    let mut fig = Figure::new("The definite integral is the signed area under the curve");
                    let n = 200;
                    let mut poly = vec![[a, 0.0]];
                    poly.extend((0..=n).map(|i| a + (b - a) * i as f64 / n as f64).map(|t| [t, g(t)]).filter(|p| p[1].is_finite()));
                    poly.push([b, 0.0]);
                    fig.shapes.push(Shape::Fill { pts: poly, color: 0 });
                    fig.curve(&g, a - pad, b + pad, 0, &calc.show(&f), false);
                    fig.caption = format!("area from {} to {} = {}", fmt(a), fmt(b), calc.show(value));
                    return Some(fig);
                }
                if visible_syms(value).iter().any(|s| *s != x) {
                    return None;
                }
                let mut fig = Figure::new("A function and an antiderivative (the area collected so far)");
                fig.curve(&g, -6.0, 6.0, 0, &calc.show(&f), false);
                fig.curve(&real_fn(value, &x), -6.0, 6.0, 1, &format!("antiderivative: {}", calc.show(value)), false);
                return Some(fig);
            }
            "limit" => {
                let x = var(1)?;
                let f = ev(0, &[&x])?;
                if visible_syms(&f).iter().any(|s| *s != x) {
                    return None;
                }
                let a = num(2)?;
                let (lo, hi) = if a.is_infinite() { if a > 0.0 { (0.0, 40.0) } else { (-40.0, 0.0) } } else { (a - 4.0, a + 4.0) };
                let mut fig = Figure::new("The curve near the limit point");
                fig.curve(&real_fn(&f, &x), lo, hi, 0, &calc.show(&f), false);
                if let Some(l) = crate::numeric::real_const(&to_float(value)).filter(|v| v.is_finite()) {
                    if a.is_finite() {
                        fig.dot([a, l], 2, &format!("limit = {}", fmt(l)));
                    } else {
                        fig.dashed(vec![[lo, l], [hi, l]], 2, &format!("limit = {}", fmt(l)));
                    }
                }
                return Some(fig);
            }
            "taylor" => {
                let x = var(1)?;
                let f = ev(0, &[&x])?;
                if visible_syms(&f).iter().any(|s| *s != x) {
                    return None;
                }
                let a = num(3).unwrap_or(0.0);
                let mut fig = Figure::new("The Taylor polynomial hugs the function near the expansion point");
                fig.curve(&real_fn(&f, &x), a - 5.0, a + 5.0, 0, &calc.show(&f), false);
                let keep = fig.bounds;
                fig.curve(&real_fn(value, &x), a - 5.0, a + 5.0, 1, "Taylor polynomial", true);
                fig.bounds = keep;
                return Some(fig);
            }
            "desolve" => {
                let (x, y) = (var(1)?, var(2)?);
                let first = match value {
                    Expr::Func(n, a) if n == "or" => a.first()?.clone(),
                    e => e.clone(),
                };
                let Expr::Rel(Rel::Eq, l, r) = &first else { return None };
                let mut fig = Figure::new("Solution curves of the differential equation");
                if **l == sym(&y) {
                    let (e, cap) = with_defaults(r, &[&x, "C1", "C2"]);
                    let family: Vec<(f64, f64)> = match (e.has_sym("C1"), e.has_sym("C2")) {
                        (true, true) => vec![(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (-1.0, 1.0), (2.0, -1.0)],
                        (true, false) => vec![(-2.0, 0.0), (-1.0, 0.0), (0.5, 0.0), (1.0, 0.0), (2.0, 0.0)],
                        _ => vec![(0.0, 0.0)],
                    };
                    for (i, (c1, c2)) in family.iter().enumerate() {
                        let ei = subst(&subst(&e, "C1", &Expr::Float(*c1)), "C2", &Expr::Float(*c2));
                        let label = if family.len() == 1 { calc.show(r) } else if e.has_sym("C2") { format!("C1 = {c1}, C2 = {c2}") } else { format!("C1 = {c1}") };
                        fig.curve(&real_fn(&ei, &x), -4.0, 4.0, i as u8, &label, false);
                    }
                    fig.caption = cap;
                } else {
                    // implicit solution psi(x, y) = 0: level curves for several constants
                    let (e, cap) = with_defaults(l, &[&x, &y, "C1"]);
                    for (i, c1) in [-2.0, -1.0, 0.0, 1.0, 2.0].iter().enumerate() {
                        let ei = subst(&e, "C1", &Expr::Float(*c1));
                        let g = |px: f64, py: f64| evalc(&ei, &[(x.as_str(), C::new(px, 0.0)), (y.as_str(), C::new(py, 0.0))]).filter(|z| z.im.abs() < 1e-9).map_or(f64::NAN, |z| z.re);
                        for seg in contour(&g, [-4.0, 4.0], [-4.0, 4.0], 120) {
                            fig.line(seg, i as u8, &format!("C1 = {c1}"));
                        }
                    }
                    fig.bounds = Some([-4.0, 4.0, -4.0, 4.0]);
                    fig.square = true;
                    fig.caption = cap;
                }
                return Some(fig);
            }
            "hamilton" | "hamiltonian" => {
                let (h, q, p) = if name == "hamilton" { (ev(0, &[&var(1)?, &var(2)?])?, var(1)?, var(2)?) } else { (value.clone(), var(1)?, "p".to_string()) };
                let (e, cap) = with_defaults(&h, &[&q, &p]);
                let g = |a: f64, b: f64| evalc(&e, &[(q.as_str(), C::new(a, 0.0)), (p.as_str(), C::new(b, 0.0))]).filter(|z| z.im.abs() < 1e-9).map_or(f64::NAN, |z| z.re);
                let mut fig = Figure::new("Phase portrait: each closed curve is a motion with constant energy H");
                let levels: Vec<f64> = [0.3, 0.9, 1.6, 2.4, 3.3].iter().map(|r| g(*r, 0.0)).filter(|v| v.is_finite()).collect();
                for (i, lv) in levels.iter().enumerate() {
                    for seg in contour(&|a, b| g(a, b) - lv, [-4.0, 4.0], [-4.0, 4.0], 130) {
                        fig.line(seg, i as u8, &format!("H = {}", fmt(*lv)));
                    }
                }
                fig.square = true;
                fig.text([3.4, 0.2], q.clone());
                fig.text([0.15, 3.6], p.clone());
                fig.bounds = Some([-4.0, 4.0, -4.0, 4.0]);
                fig.caption = cap;
                return Some(fig);
            }
            "fieldsol" | "fieldcheck" => {
                let field = var(0)?;
                let u = if name == "fieldsol" { value.clone() } else { ev(1, &["x", "y", "z", "t", "theta"])? };
                let u = pulse(&u);
                if field == "kepler" {
                    let (e, cap) = with_defaults(&u, &["theta"]);
                    let r = real_fn(&e, "theta");
                    let mut fig = Figure::new("Orbit r(theta) with the source at the origin");
                    fig.line((0..=400).map(|i| i as f64 / 400.0 * TAU).map(|t| [r(t) * t.cos(), r(t) * t.sin()]).filter(|p| p[0].is_finite() && p[1].is_finite()).collect(), 0, "orbit");
                    fig.dot([0.0, 0.0], 1, "source (focus)");
                    fig.equal = true;
                    fig.caption = cap;
                    return Some(fig);
                }
                let main = if u.has_sym("x") { "x" } else if u.has_sym("t") { "t" } else { return None };
                let (e, cap) = with_defaults(&u, &[main, "t"]);
                let e = if main == "x" { e } else { e.clone() };
                let mut fig = Figure::new(if main == "x" { "Snapshots of the field along x at successive times" } else { "The solution against time" });
                let (lo, hi) = if main == "x" { (-6.0, 6.0) } else { (0.0, 12.0) };
                if main == "x" && e.has_sym("t") {
                    // start slightly after zero when the solution is singular at t = 0
                    let t0 = if part_fn(&subst(&e, "t", &Expr::Float(0.0)), "x", false)(0.37).is_finite() { 0.0 } else { 0.1 };
                    for (i, dt) in [0.0, 0.4, 0.8, 1.2].iter().enumerate() {
                        let et = subst(&e, "t", &Expr::Float(t0 + dt));
                        fig.curve(&part_fn(&et, "x", false), lo, hi, i as u8, &format!("t = {}", fmt(t0 + dt)), false);
                    }
                    let e0 = subst(&e, "t", &Expr::Float(t0 + 0.4));
                    if (-3..=3).any(|k| part_fn(&e0, "x", true)(k as f64 * 0.7).abs() > 1e-9) {
                        fig.curve(&part_fn(&e0, "x", true), lo, hi, 5, &format!("imaginary part, t = {}", fmt(t0 + 0.4)), true);
                    }
                } else {
                    fig.curve(&part_fn(&e, main, false), lo, hi, 0, &calc.show(&u), false);
                }
                fig.caption = cap;
                return Some(fig);
            }
            "galilean" | "lorentzt" => {
                let beta = num(1).map(|v| v / 299_792_458.0).filter(|b| b.abs() < 1.0 && b.abs() > 1e-6).unwrap_or(0.5);
                return Some(spacetime(beta, name == "lorentzt"));
            }
            "lorentz" | "timedil" | "lencon" | "relmom" | "relenergy" | "relke" | "reldoppler" => {
                let vi = if name == "lorentz" { 0 } else { 1 };
                let beta = num(vi).map(|v| v / 299_792_458.0).filter(|b| b.abs() < 1.0);
                let inverse = name == "lencon";
                let mut fig = Figure::new(if inverse { "Length contraction factor against v/c" } else { "Lorentz factor against v/c" });
                let g = move |b: f64| if inverse { (1.0 - b * b).sqrt() } else { 1.0 / (1.0 - b * b).sqrt() };
                fig.curve(&g, 0.0, 0.995, 0, if inverse { "sqrt(1 - v^2/c^2)" } else { "1/sqrt(1 - v^2/c^2)" }, false);
                fig.bounds = Some([0.0, 1.0, 0.0, if inverse { 1.1 } else { 8.0 }]);
                if let Some(b) = beta {
                    fig.dot([b.abs(), g(b)], 2, &format!("v/c = {}, factor = {}", fmt(b.abs()), fmt(g(b))));
                }
                return Some(fig);
            }
            "veladd" => {
                let bv = num(1).map(|v| v / 299_792_458.0).filter(|b| b.abs() < 1.0).unwrap_or(0.5);
                let mut fig = Figure::new("Velocity addition: relativistic curve and the straight Galilean sum");
                fig.curve(&|u: f64| (u + bv) / (1.0 + u * bv), 0.0, 1.0, 0, "(u + v)/(1 + u*v/c^2)", false);
                fig.curve(&|u: f64| u + bv, 0.0, 1.0, 1, "u + v", true);
                fig.dashed(vec![[0.0, 1.0], [1.0, 1.0]], 3, "c");
                if let Some(bu) = num(0).map(|v| v / 299_792_458.0).filter(|b| b.abs() <= 1.0) {
                    fig.dot([bu, (bu + bv) / (1.0 + bu * bv)], 2, "this case");
                }
                fig.caption = format!("axes in units of c; v/c = {}", fmt(bv));
                fig.bounds = Some([0.0, 1.0, 0.0, 1.6]);
                return Some(fig);
            }
            "vecadd" => {
                let (u, v) = (vec2(&ev(0, &[])?)?, vec2(&ev(1, &[])?)?);
                return Some(parallelogram(u, v, ["u", "v", "u + v"], "Vector addition: the diagonal of the parallelogram"));
            }
            "det" => {
                let Expr::Matrix(m) = ev(0, &[])? else { return None };
                if m.len() != 2 || m[0].len() != 2 {
                    return None;
                }
                let g = |e: &Expr| crate::numeric::real_const(&to_float(e));
                let (u, v) = ([g(&m[0][0])?, g(&m[1][0])?], [g(&m[0][1])?, g(&m[1][1])?]);
                let mut fig = parallelogram(u, v, ["column 1", "column 2", "sum"], "The determinant is the signed area of this parallelogram");
                fig.shapes.insert(0, Shape::Fill { pts: vec![[0.0, 0.0], u, [u[0] + v[0], u[1] + v[1]], v], color: 0 });
                fig.caption = format!("area = |det| = {}", fmt((u[0] * v[1] - u[1] * v[0]).abs()));
                return Some(fig);
            }
            "inscribedangle" | "tangentchord" => return Some(circle_figure(if name == "inscribedangle" { "inscribed" } else { "tangentchord" }, num(0)? * kdeg)),
            "chordlen" | "arclength" | "sectorarea" | "segmentarea" => return Some(circle_figure("central", num(1)? * kdeg)),
            "centralangle" | "sagitta" => {
                let ang = if name == "centralangle" { crate::numeric::real_const(&to_float(value))? * kdeg } else { 2.0 * (num(1)? / (2.0 * num(0)?)).asin().to_degrees() };
                return Some(circle_figure("central", ang));
            }
            "sphsss" => return Some(sphere_triangle(num(0)? * kdeg, num(1)? * kdeg, list_value(value, "C")? * kdeg)),
            "sphsas" => return Some(sphere_triangle(num(0)? * kdeg, num(2)? * kdeg, num(1)? * kdeg)),
            "sphasa" => return Some(sphere_triangle(list_value(value, "a")? * kdeg, list_value(value, "b")? * kdeg, list_value(value, "C")? * kdeg)),
            "sphaaa" => return Some(sphere_triangle(list_value(value, "a")? * kdeg, list_value(value, "b")? * kdeg, num(2)? * kdeg)),
            "sphdist" | "sphbearing" | "angsep" => {
                // both argument orders reduce to (latitude-like, longitude-like) pairs
                let geo = name != "angsep";
                let (lat1, lon1, lat2, lon2) = if geo { (num(0)?, num(1)?, num(2)?, num(3)?) } else { (num(1)?, num(0)?, num(3)?, num(2)?) };
                let (p1, p2) = (latlon(lat1 * kdeg, lon1 * kdeg), latlon(lat2 * kdeg, lon2 * kdeg));
                let view = Sphere::facing([p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2] + 0.35]);
                let mut fig = Figure::new("Great-circle arc between two points on the sphere");
                view.outline(&mut fig);
                let pole: V3 = [0.0, 0.0, 1.0];
                for lat in [-60.0f64, -30.0, 0.0, 30.0, 60.0] {
                    view.grid_ring(&mut fig, pole, (90.0 - lat).to_radians());
                }
                for lon in [0.0f64, 30.0, 60.0, 90.0, 120.0, 150.0] {
                    view.grid_ring(&mut fig, latlon(0.0, lon + 90.0), PI_F / 2.0);
                }
                view.arc(&mut fig, p1, p2, 0, "great-circle arc");
                view.label(&mut fig, p1, "1", 1);
                view.label(&mut fig, p2, "2", 2);
                view.label(&mut fig, pole, if geo { "N" } else { "pole" }, GREY);
                return Some(fig);
            }
            "navtri" | "eq2hor" => {
                // horizon frame: x north, y east, z up
                let (lat, dec, ha) = if name == "navtri" { (num(0)? * kdeg, num(1)? * kdeg, num(2)? * kdeg) } else { (num(2)? * kdeg, num(1)? * kdeg, num(0)? * kdeg) };
                let (alt, az) = crate::astro::eq_to_hor(ha, dec, lat);
                let hz = |alt: f64, az: f64| -> V3 { [alt.to_radians().cos() * az.to_radians().cos(), alt.to_radians().cos() * az.to_radians().sin(), alt.to_radians().sin()] };
                let (z, p, x) = ([0.0, 0.0, 1.0], hz(lat, 0.0), hz(alt, az));
                let view = Sphere::facing([-0.35, if az > 180.0 { 0.75 } else { -0.75 }, 0.55]);
                let mut fig = Figure::new("Navigational triangle on the celestial sphere: pole P, zenith Z, body X");
                view.outline(&mut fig);
                view.ring(&mut fig, z, PI_F / 2.0, GREY, "horizon");
                view.ring(&mut fig, p, PI_F / 2.0, 3, "celestial equator");
                view.ring(&mut fig, p, (90.0 - dec).to_radians(), 4, "daily path of the body");
                view.arc(&mut fig, p, z, 0, "co-latitude PZ");
                view.arc(&mut fig, z, x, 1, "zenith distance ZX");
                view.arc(&mut fig, p, x, 2, "polar distance PX");
                view.label(&mut fig, p, "P", 0);
                view.label(&mut fig, z, "Z", 1);
                view.label(&mut fig, x, "X", 2);
                view.label(&mut fig, hz(0.0, 0.0), "N", GREY);
                view.label(&mut fig, hz(0.0, 90.0), "E", GREY);
                view.label(&mut fig, hz(0.0, 180.0), "S", GREY);
                view.label(&mut fig, hz(0.0, 270.0), "W", GREY);
                fig.caption = format!("altitude {} deg, azimuth {} deg", fmt(alt), fmt(az));
                return Some(fig);
            }
            "kepler" => {
                let (m, e) = (num(0)? * kdeg, num(1)?);
                let ea = crate::astro::kepler_e(m.to_radians(), e);
                let b = (1.0 - e * e).sqrt();
                let mut fig = Figure::new("Kepler's equation: eccentric anomaly E on the circle, the body on the ellipse");
                fig.line(circle_pts([0.0, 0.0], 1.0, 0.0, TAU), GREY, "auxiliary circle");
                fig.line((0..=200).map(|i| i as f64 / 200.0 * TAU).map(|t| [t.cos(), b * t.sin()]).collect(), 0, "orbit");
                let (qc, qe, focus) = ([ea.cos(), ea.sin()], [ea.cos(), b * ea.sin()], [e, 0.0]);
                fig.line(vec![[0.0, 0.0], qc], 1, "eccentric anomaly E");
                fig.dashed(vec![qc, [ea.cos(), 0.0]], GREY, "");
                fig.line(vec![focus, qe], 2, "radius from the focus");
                fig.dot(focus, 2, "focus");
                fig.dot(qe, 0, "body");
                fig.dot(qc, 1, "");
                fig.caption = format!("M = {} deg, e = {}, E = {} deg", fmt(m), fmt(e), fmt(ea.to_degrees()));
                fig.bounds = Some([-1.3, 1.3, -1.3, 1.3]);
                fig.equal = true;
                return Some(fig);
            }
            "linreg" | "quadreg" | "expreg" | "lnreg" | "pwrreg" => {
                let (Expr::List(xs), Expr::List(ys)) = (ev(0, &[])?, ev(1, &[])?) else { return None };
                let g = |e: &Expr| crate::numeric::real_const(&to_float(e));
                let pts: Vec<[f64; 2]> = xs.iter().zip(&ys).filter_map(|(a, b)| Some([g(a)?, g(b)?])).collect();
                let Expr::List(items) = value else { return None };
                let Expr::Rel(_, _, model) = items.first()? else { return None };
                let (x0, x1) = (pts.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min), pts.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max));
                let pad = ((x1 - x0) * 0.15).max(0.5);
                let mut fig = Figure::new("Data points and the fitted curve");
                fig.curve(&real_fn(model, "x"), x0 - pad, x1 + pad, 0, &calc.show(model), false);
                for p in pts {
                    fig.dot(p, 1, "");
                }
                return Some(fig);
            }
            "normalcdf" | "normalpdf" => {
                let cdf = name == "normalcdf";
                let (mu, sd) = (num(if cdf { 2 } else { 1 }).unwrap_or(0.0), num(if cdf { 3 } else { 2 }).unwrap_or(1.0));
                let pdf = move |x: f64| (-0.5 * ((x - mu) / sd).powi(2)).exp() / (sd * TAU.sqrt());
                let (lo, hi) = (mu - 4.0 * sd, mu + 4.0 * sd);
                let mut fig = Figure::new(if cdf { "The probability is the area under the normal curve" } else { "Normal density" });
                if cdf {
                    let (a, b) = (num(0)?.max(lo), num(1)?.min(hi));
                    let mut poly = vec![[a, 0.0]];
                    poly.extend((0..=120).map(|i| a + (b - a) * i as f64 / 120.0).map(|t| [t, pdf(t)]));
                    poly.push([b, 0.0]);
                    fig.shapes.push(Shape::Fill { pts: poly, color: 0 });
                    fig.caption = format!("shaded area = {}", calc.show(value));
                } else {
                    let x = num(0)?;
                    fig.dot([x, pdf(x)], 2, "");
                }
                fig.curve(&pdf, lo, hi, 0, "normal density", false);
                return Some(fig);
            }
            "polar" | "rect" => {
                let z = if name == "polar" { ev(0, &[])? } else { value.clone() };
                let z = crate::numeric::eval_const(&z)?;
                let mut fig = parallelogram([z.re, 0.0], [0.0, z.im], ["real part", "imaginary part", "z"], "A complex number as a point in the plane");
                let r = z.norm();
                fig.line(circle_pts([0.0, 0.0], r * 0.3, 0.0, z.arg()), 3, "argument");
                fig.caption = format!("modulus {}, argument {} rad", fmt(r), fmt(z.arg()));
                return Some(fig);
            }
            _ => {}
        }
    }
    // any result with exactly one free variable: draw it
    let body = match value {
        Expr::Rel(Rel::Eq, l, r) if matches!(**l, Expr::Sym(_)) => (**r).clone(),
        Expr::Add(_) | Expr::Mul(_) | Expr::Pow(..) | Expr::Func(..) => value.clone(),
        _ => return None,
    };
    let syms = visible_syms(&body);
    if syms.len() != 1 || body.has(&|e| matches!(e, Expr::List(_) | Expr::Matrix(_) | Expr::Rel(..) | Expr::Str(_))) {
        return None;
    }
    let x = &syms[0];
    let f = real_fn(&body, x);
    if !(-40..=40).any(|k| f(k as f64 * 0.251 + 0.013).is_finite()) {
        return None;
    }
    let mut fig = Figure::new(format!("Graph against {x}"));
    fig.curve(&f, -10.0, 10.0, 0, &calc.show(&body), false);
    if x == "x" {
        fig.graph_src = Some(calc.show(&body));
    }
    Some(fig)
}

fn var_name(e: &Expr) -> Option<String> {
    match e {
        Expr::Sym(s) => Some(s.clone()),
        _ => None,
    }
}
