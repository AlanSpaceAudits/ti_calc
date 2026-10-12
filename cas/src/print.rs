//! Infix output.
use crate::expr::*;
use crate::simplify::{neg, pow};
use num_traits::{One, Signed};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Notation {
    Normal,
    Sci,
    Eng,
}

#[derive(Clone, Copy, Debug)]
pub struct Fmt {
    pub notation: Notation,
    /// significant digits in float output
    pub digits: usize,
    /// fixed number of decimals, if set
    pub fix: Option<usize>,
}

impl Default for Fmt {
    fn default() -> Self {
        Fmt { notation: Notation::Normal, digits: 12, fix: None }
    }
}

fn trim_zeros(s: &str) -> String {
    if s.contains('.') { s.trim_end_matches('0').to_string() } else { s.to_string() }
}

pub fn fmt_float(x: f64, f: &Fmt) -> String {
    if x.is_nan() {
        return "undef".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "∞".into() } else { "-∞".into() };
    }
    if x == 0.0 {
        return match f.fix {
            Some(n) if n > 0 => format!("{:.*}", n, 0.0),
            _ => "0.".into(),
        };
    }
    let digits = f.digits.clamp(1, 17);
    let sci = |x: f64, step: i32| -> String {
        let s = format!("{:.*e}", digits - 1, x);
        let (m, e) = s.split_once('e').unwrap();
        let mut exp: i32 = e.parse().unwrap();
        let mut mant: f64 = m.parse().unwrap();
        if step == 3 {
            let shift = exp.rem_euclid(3);
            mant *= 10f64.powi(shift);
            exp -= shift;
        }
        let dec = match f.fix {
            Some(n) => n,
            None => (digits as i32 - 1 - if step == 3 { exp.rem_euclid(3).max(0) } else { 0 }).max(0) as usize,
        };
        let ms = format!("{:.*}", dec, mant);
        let ms = if f.fix.is_some() { ms } else { trim_zeros(&ms) };
        format!("{ms}E{exp}")
    };
    match f.notation {
        Notation::Sci => return sci(x, 1),
        Notation::Eng => return sci(x, 3),
        Notation::Normal => {}
    }
    let exp10 = x.abs().log10().floor() as i32;
    if exp10 < -5 || exp10 >= digits.max(10) as i32 {
        return sci(x, 1);
    }
    if let Some(n) = f.fix {
        return format!("{:.*}", n, x);
    }
    let dec = (digits as i32 - 1 - exp10).max(0) as usize;
    let s = trim_zeros(&format!("{:.*}", dec, x));
    if s.contains('.') { s } else { format!("{s}.") }
}

fn term_degree(t: &Expr) -> f64 {
    let fs: Vec<&Expr> = if let Expr::Mul(v) = t { v.iter().collect() } else { vec![t] };
    fs.iter()
        .map(|f| match f {
            Expr::Sym(_) => 1.0,
            Expr::Pow(b, n) if matches!(**b, Expr::Sym(_)) => n.as_f64().unwrap_or(0.0),
            _ => 0.0,
        })
        .sum()
}

const P_OR: u8 = 0;
const P_AND: u8 = 1;
const P_REL: u8 = 2;
const P_ADD: u8 = 3;
const P_MUL: u8 = 4;
const P_NEG: u8 = 5;
const P_POW: u8 = 6;
const P_ATOM: u8 = 7;

fn wrap(s: (String, u8), min: u8) -> String {
    if s.1 < min { format!("({})", s.0) } else { s.0 }
}

fn is_neg_exp(e: &Expr) -> bool {
    e.has_neg_coeff()
}

fn show_mul(v: &[Expr], f: &Fmt) -> (String, u8) {
    let mut negative = false;
    let mut num: Vec<Expr> = Vec::new();
    let mut den: Vec<Expr> = Vec::new();
    for (i, fac) in v.iter().enumerate() {
        match fac {
            Expr::Num(r) if i == 0 => {
                negative = r.is_negative();
                let a = r.abs();
                if !a.numer().is_one() {
                    num.push(big(a.numer().clone()));
                }
                if !a.denom().is_one() {
                    den.push(big(a.denom().clone()));
                }
            }
            Expr::Float(x) if i == 0 && (*x < 0.0 || *x == 1.0) && v.len() > 1 => {
                negative = *x < 0.0;
                if x.abs() != 1.0 {
                    num.push(Expr::Float(-x));
                }
            }
            Expr::Pow(b, e) if is_neg_exp(e) => den.push(pow((**b).clone(), neg((**e).clone()))),
            other => num.push(other.clone()),
        }
    }
    // merge radicals with equal exponents: sqrt(2)*sqrt(3) prints as sqrt(6)
    let merge = |fs: Vec<Expr>| -> Vec<Expr> {
        let mut out: Vec<Expr> = Vec::new();
        for fac in fs {
            if let Expr::Pow(b, e) = &fac {
                if let (Some(br), Some(er)) = (b.as_rat(), e.as_rat()) {
                    if !er.is_integer() && br.is_integer() {
                        if let Some(slot) = out.iter_mut().find(|o| matches!(o, Expr::Pow(ob, oe) if oe.as_rat() == Some(er) && ob.as_int().is_some())) {
                            if let Expr::Pow(ob, _) = slot {
                                **ob = Expr::Num(ob.as_rat().unwrap() * br);
                            }
                            continue;
                        }
                    }
                }
            }
            out.push(fac);
        }
        out
    };
    let (mut num, den) = (merge(num), merge(den));
    // exponentials and functions print after plain variables: x*e^x, not e^x*x
    let late = |x: &Expr| match x {
        Expr::Func(..) => 2,
        Expr::Pow(b, _) if matches!(**b, Expr::Const(Const::E) | Expr::Func(..)) => 2,
        Expr::Const(Const::I) => 1,
        _ => 0,
    };
    num.sort_by_key(late);
    let join = |fs: &[Expr]| fs.iter().map(|x| wrap(show_p(x, f), P_NEG + 1)).collect::<Vec<_>>().join("*");
    let mut s = if num.is_empty() { "1".to_string() } else { join(&num) };
    if !den.is_empty() {
        if den.len() == 1 {
            s = format!("{}/{}", s, wrap(show_p(&den[0], f), P_POW));
        } else {
            s = format!("{}/({})", s, join(&den));
        }
    }
    if negative { (format!("-{s}"), P_NEG) } else { (s, P_MUL) }
}

fn show_p(e: &Expr, f: &Fmt) -> (String, u8) {
    match e {
        Expr::Num(r) => {
            let p = if r.is_negative() {
                P_NEG
            } else if r.is_integer() {
                P_ATOM
            } else {
                P_MUL
            };
            (if r.is_integer() { r.numer().to_string() } else { format!("{}/{}", r.numer(), r.denom()) }, p)
        }
        Expr::Float(x) => (fmt_float(*x, f), if *x < 0.0 { P_NEG } else { P_ATOM }),
        Expr::Sym(s) => (s.clone(), P_ATOM),
        Expr::Const(c) => (
            match c {
                Const::Pi => "π",
                Const::E => "e",
                Const::I => "i",
                Const::Inf => "∞",
                Const::Undef => "undef",
            }
            .to_string(),
            P_ATOM,
        ),
        Expr::Add(v) => {
            let mut terms: Vec<&Expr> = v.iter().collect();
            let has_i = |t: &Expr| t.has(&|x| matches!(x, Expr::Const(Const::I)));
            let syms: Vec<String> = e.free_syms().into_iter().collect();
            let exps = |t: &Expr| -> Vec<f64> {
                let fs: Vec<&Expr> = if let Expr::Mul(v) = t { v.iter().collect() } else { vec![t] };
                syms.iter()
                    .map(|s| {
                        fs.iter()
                            .map(|f| match f {
                                Expr::Sym(n) if n == s => 1.0,
                                Expr::Pow(b, k) if matches!(&**b, Expr::Sym(n) if n == s) => k.as_f64().unwrap_or(0.0),
                                _ => 0.0,
                            })
                            .sum()
                    })
                    .collect()
            };
            terms.sort_by(|a, b| {
                has_i(a)
                    .cmp(&has_i(b))
                    .then(term_degree(b).total_cmp(&term_degree(a)))
                    .then_with(|| exps(b).partial_cmp(&exps(a)).unwrap_or(std::cmp::Ordering::Equal))
                    .then(a.is_number().cmp(&b.is_number()))
            });
            if terms[0].has_neg_coeff() {
                let d0 = term_degree(terms[0]);
                if let Some(k) = terms.iter().position(|t| !t.has_neg_coeff() && term_degree(t) == d0 && !has_i(t)) {
                    let t = terms.remove(k);
                    terms.insert(0, t);
                }
            }
            let mut s = String::new();
            for (i, t) in terms.iter().enumerate() {
                if i == 0 {
                    s.push_str(&wrap(show_p(t, f), P_ADD + 1));
                } else if t.has_neg_coeff() {
                    s.push_str(" - ");
                    s.push_str(&wrap(show_p(&neg((*t).clone()), f), P_MUL));
                } else {
                    s.push_str(" + ");
                    s.push_str(&wrap(show_p(t, f), P_MUL));
                }
            }
            (s, P_ADD)
        }
        Expr::Mul(v) => show_mul(v, f),
        Expr::Pow(b, x) => {
            if **x == rat(1, 2) {
                return (format!("sqrt({})", show(b, f)), P_ATOM);
            }
            if is_neg_exp(x) {
                return show_mul(std::slice::from_ref(e), f);
            }
            let xs = show_p(x, f);
            let simple = matches!(**x, Expr::Sym(_) | Expr::Const(_)) || x.as_int().is_some_and(|n| !n.is_negative());
            let xs = if simple { xs.0 } else { format!("({})", xs.0) };
            (format!("{}^{}", wrap(show_p(b, f), P_ATOM), xs), P_POW)
        }
        Expr::Func(name, a) => match name.as_str() {
            "and" if a.len() == 2 && matches!((&a[0], &a[1]), (Expr::Rel(_, _, m1), Expr::Rel(_, m2, _)) if m1 == m2 && matches!(**m1, Expr::Sym(_))) => {
                let (Expr::Rel(r1, lo, mid), Expr::Rel(r2, _, hi)) = (&a[0], &a[1]) else { unreachable!() };
                let op = |r: &Rel| match r {
                    Rel::Lt => "<",
                    Rel::Le => "<=",
                    Rel::Gt => ">",
                    Rel::Ge => ">=",
                    Rel::Eq => "=",
                    Rel::Ne => "!=",
                };
                (format!("{} {} {} {} {}", wrap(show_p(lo, f), P_ADD), op(r1), show(mid, f), op(r2), wrap(show_p(hi, f), P_ADD)), P_REL)
            }
            "and" | "or" => {
                let (p, sep) = if name == "and" { (P_AND, " and ") } else { (P_OR, " or ") };
                (a.iter().map(|x| wrap(show_p(x, f), p + 1)).collect::<Vec<_>>().join(sep), p)
            }
            "factorial" if a.len() == 1 => (format!("{}!", wrap(show_p(&a[0], f), P_ATOM)), P_POW),
            _ => (format!("{}({})", name, a.iter().map(|x| show(x, f)).collect::<Vec<_>>().join(", ")), P_ATOM),
        },
        Expr::Rel(r, a, b) => {
            let op = match r {
                Rel::Eq => "=",
                Rel::Ne => "!=",
                Rel::Lt => "<",
                Rel::Le => "<=",
                Rel::Gt => ">",
                Rel::Ge => ">=",
            };
            (format!("{} {} {}", wrap(show_p(a, f), P_ADD), op, wrap(show_p(b, f), P_ADD)), P_REL)
        }
        Expr::List(v) => (format!("{{{}}}", v.iter().map(|x| show(x, f)).collect::<Vec<_>>().join(", ")), P_ATOM),
        Expr::Matrix(m) => (format!("[{}]", m.iter().map(|r| format!("[{}]", r.iter().map(|x| show(x, f)).collect::<Vec<_>>().join(", "))).collect::<Vec<_>>().join(", ")), P_ATOM),
        Expr::Bool(b) => (b.to_string(), P_ATOM),
        Expr::Str(s) => (s.clone(), P_ATOM),
    }
}

pub fn show(e: &Expr, f: &Fmt) -> String {
    show_p(e, f).0
}

