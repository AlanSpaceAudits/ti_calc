//! Automatic simplification: every constructor here returns a canonical expression.
use crate::expr::*;
use crate::matrix;
use crate::numeric::{self, is_prime_u64};
use num_bigint::BigInt;
use num_complex::Complex64 as C;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

pub fn nadd(a: &Expr, b: &Expr) -> Expr {
    match (a, b) {
        (Expr::Num(x), Expr::Num(y)) => Expr::Num(x + y),
        _ => Expr::Float(a.as_f64().unwrap_or(f64::NAN) + b.as_f64().unwrap_or(f64::NAN)),
    }
}
pub fn nmul(a: &Expr, b: &Expr) -> Expr {
    match (a, b) {
        (Expr::Num(x), Expr::Num(y)) => Expr::Num(x * y),
        _ => Expr::Float(a.as_f64().unwrap_or(f64::NAN) * b.as_f64().unwrap_or(f64::NAN)),
    }
}
fn nneg(a: &Expr) -> Expr {
    nmul(a, &int(-1))
}

pub fn add(a: Expr, b: Expr) -> Expr {
    add_v(vec![a, b])
}
pub fn sub(a: Expr, b: Expr) -> Expr {
    add_v(vec![a, neg(b)])
}
pub fn mul(a: Expr, b: Expr) -> Expr {
    mul_v(vec![a, b])
}
pub fn div(a: Expr, b: Expr) -> Expr {
    mul_v(vec![a, pow(b, int(-1))])
}
pub fn neg(a: Expr) -> Expr {
    mul_v(vec![int(-1), a])
}
pub fn sqrt(a: Expr) -> Expr {
    pow(a, rat(1, 2))
}
pub fn f1(name: &str, a: Expr) -> Expr {
    func(name, vec![a])
}

pub fn from_c(z: C) -> Expr {
    if z.im == 0.0 || z.im.abs() < 1e-15 * z.re.abs() {
        Expr::Float(z.re)
    } else if z.re.abs() < 1e-15 * z.im.abs() {
        gauss(int(0), Expr::Float(z.im))
    } else {
        gauss(Expr::Float(z.re), Expr::Float(z.im))
    }
}

/// Recognize a literal complex number: returns numeric (re, im).
pub fn as_gauss(e: &Expr) -> Option<(Expr, Expr)> {
    let im_part = |t: &Expr| -> Option<Expr> {
        match t {
            Expr::Const(Const::I) => Some(int(1)),
            Expr::Mul(v) if v.len() == 2 && v[0].is_number() && v[1] == I => Some(v[0].clone()),
            _ => None,
        }
    };
    if e.is_number() {
        return Some((e.clone(), int(0)));
    }
    if let Some(b) = im_part(e) {
        return Some((int(0), b));
    }
    if let Expr::Add(v) = e {
        if v.len() == 2 {
            for (x, y) in [(&v[0], &v[1]), (&v[1], &v[0])] {
                if x.is_number() {
                    if let Some(b) = im_part(y) {
                        return Some((x.clone(), b));
                    }
                }
            }
        }
    }
    None
}

fn gauss(re: Expr, im: Expr) -> Expr {
    if im.is_zero() {
        return re;
    }
    let ip = if im.is_one() { I } else { Expr::Mul(vec![im, I]) };
    if re.is_zero() { ip } else { Expr::Add(vec![re, ip]) }
}

/// Prime factorization by trial division. A large leftover cofactor is returned as-is.
pub fn factor_int(n: &BigInt) -> Vec<(BigInt, u32)> {
    let mut out = Vec::new();
    let mut n = n.abs();
    if n <= BigInt::one() {
        return out;
    }
    let mut p = 2u64;
    // ponytail: trial division to 1e6; cofactors with two larger prime factors stay unfactored
    while p <= 1_000_000 {
        let bp = BigInt::from(p);
        if &bp * &bp > n {
            break;
        }
        let mut k = 0;
        while (&n % &bp).is_zero() {
            n /= &bp;
            k += 1;
        }
        if k > 0 {
            out.push((bp, k));
        }
        p += if p == 2 { 1 } else { 2 };
    }
    if n > BigInt::one() {
        // perfect-power check on the leftover
        let mut done = false;
        if n.to_u64().is_none_or(|v| !is_prime_u64(v)) {
            for k in (2..=n.bits().min(64) as u32).rev() {
                let r = n.nth_root(k);
                if r.pow(k) == n {
                    out.push((r, k));
                    done = true;
                    break;
                }
            }
        }
        if !done {
            out.push((n, 1));
        }
    }
    out
}

fn rat_root_pow(b: &Rat, e: &Rat) -> Expr {
    let raw = || Expr::Pow(Box::new(Expr::Num(b.clone())), Box::new(Expr::Num(e.clone())));
    if b.is_negative() {
        let nb = -b;
        let pos = if nb.is_one() { int(1) } else { rat_root_pow(&nb, e) };
        let (p, q) = (e.numer(), e.denom());
        return if q.is_odd() {
            if p.is_odd() { neg(pos) } else { pos }
        } else if *q == BigInt::from(2) {
            mul_v(vec![pow(I, big(p.clone())), pos])
        } else {
            let m1 = Expr::Pow(Box::new(int(-1)), Box::new(Expr::Num(e.clone())));
            if pos.is_one() { m1 } else { mul_v(vec![m1, pos]) }
        };
    }
    if e.numer().bits() > 12 {
        return raw();
    }
    let mut coef = Rat::one();
    let mut rads: Vec<Expr> = Vec::new();
    let mut parts: Vec<(BigInt, Rat)> = Vec::new();
    for (p, k) in factor_int(b.numer()) {
        parts.push((p, Rat::from_integer(BigInt::from(k)) * e));
    }
    for (p, k) in factor_int(b.denom()) {
        parts.push((p, -Rat::from_integer(BigInt::from(k)) * e));
    }
    parts.sort();
    for (p, ex) in parts {
        let w = ex.floor();
        let fr = &ex - &w;
        let Some(wi) = w.to_integer().to_i32() else { return raw() };
        coef *= Rat::from_integer(p.clone()).pow(wi);
        if !fr.is_zero() {
            rads.push(Expr::Pow(Box::new(big(p)), Box::new(Expr::Num(fr))));
        }
    }
    if rads.is_empty() {
        return Expr::Num(coef);
    }
    if !coef.is_one() {
        rads.insert(0, Expr::Num(coef));
    }
    if rads.len() == 1 { rads.pop().unwrap() } else { Expr::Mul(rads) }
}

fn is_struct(e: &Expr) -> bool {
    matches!(e, Expr::Matrix(_) | Expr::List(_))
}

fn zip_list(a: &Expr, b: &Expr, f: &dyn Fn(Expr, Expr) -> Expr) -> Expr {
    match (a, b) {
        (Expr::List(x), Expr::List(y)) => {
            if x.len() != y.len() {
                return UNDEF;
            }
            Expr::List(x.iter().zip(y).map(|(p, q)| f(p.clone(), q.clone())).collect())
        }
        (Expr::List(x), s) => Expr::List(x.iter().map(|p| f(p.clone(), s.clone())).collect()),
        (s, Expr::List(y)) => Expr::List(y.iter().map(|q| f(s.clone(), q.clone())).collect()),
        _ => UNDEF,
    }
}

fn add_struct(flat: Vec<Expr>) -> Expr {
    let mut it = flat.into_iter();
    let mut acc = it.next().unwrap();
    for b in it {
        acc = match (&acc, &b) {
            (Expr::Matrix(x), Expr::Matrix(y)) => matrix::add(x, y).map(Expr::Matrix).unwrap_or(UNDEF),
            (Expr::Matrix(_), _) | (_, Expr::Matrix(_)) => UNDEF,
            (Expr::List(_), _) | (_, Expr::List(_)) => zip_list(&acc, &b, &add),
            _ => add(acc, b),
        };
    }
    acc
}

fn mul_struct(flat: Vec<Expr>) -> Expr {
    let mut it = flat.into_iter();
    let mut acc = it.next().unwrap();
    for b in it {
        acc = match (&acc, &b) {
            (Expr::Matrix(x), Expr::Matrix(y)) => matrix::mul(x, y).map(Expr::Matrix).unwrap_or(UNDEF),
            (Expr::Matrix(_), Expr::List(_)) | (Expr::List(_), Expr::Matrix(_)) => UNDEF,
            (Expr::Matrix(x), s) | (s, Expr::Matrix(x)) => Expr::Matrix(matrix::scale(x, s)),
            (Expr::List(_), _) | (_, Expr::List(_)) => zip_list(&acc, &b, &mul),
            _ => mul(acc, b),
        };
    }
    acc
}

fn split_coeff(t: Expr) -> (Expr, Expr) {
    match t {
        Expr::Mul(mut v) if v[0].is_number() => {
            let c = v.remove(0);
            (c, if v.len() == 1 { v.pop().unwrap() } else { Expr::Mul(v) })
        }
        t => (int(1), t),
    }
}

fn with_coeff(c: Expr, rest: Expr) -> Expr {
    if c.is_one() {
        return rest;
    }
    match rest {
        Expr::Mul(mut v) => {
            v.insert(0, c);
            Expr::Mul(v)
        }
        r => Expr::Mul(vec![c, r]),
    }
}

pub fn add_v(args: Vec<Expr>) -> Expr {
    let mut flat = Vec::with_capacity(args.len());
    for a in args {
        match a {
            Expr::Add(v) => flat.extend(v),
            // c*(a + b) inside a sum distributes so like terms can combine
            Expr::Mul(v) if v.len() == 2 && v[0].is_number() && matches!(v[1], Expr::Add(_)) => {
                if let Expr::Add(ts) = &v[1] {
                    for t in ts {
                        match mul_v(vec![v[0].clone(), t.clone()]) {
                            Expr::Add(w) => flat.extend(w),
                            x => flat.push(x),
                        }
                    }
                }
            }
            x => flat.push(x),
        }
    }
    if flat.iter().any(|e| e.is_undef()) {
        return UNDEF;
    }
    if flat.iter().any(is_struct) {
        return add_struct(flat);
    }
    let (mut acc, mut acc_im) = (int(0), int(0));
    let mut terms: Vec<(Expr, Expr)> = Vec::new();
    let (mut pinf, mut ninf) = (false, false);
    for t in flat {
        if t.is_number() {
            acc = nadd(&acc, &t);
            continue;
        }
        let (c, r) = split_coeff(t);
        if r == INF {
            if c.is_neg_number() { ninf = true } else { pinf = true }
        } else if r == I {
            acc_im = nadd(&acc_im, &c);
        } else if let Some(slot) = terms.iter_mut().find(|(rr, _)| *rr == r) {
            slot.1 = nadd(&slot.1, &c);
        } else {
            terms.push((r, c));
        }
    }
    if pinf && ninf {
        return UNDEF;
    }
    if pinf || ninf {
        return if pinf { INF } else { Expr::Mul(vec![int(-1), INF]) };
    }
    if matches!(acc, Expr::Float(_)) || matches!(acc_im, Expr::Float(_)) {
        let mut z = C::new(acc.as_f64().unwrap(), acc_im.as_f64().unwrap());
        let mut exact_im = !matches!(acc_im, Expr::Float(_));
        terms.retain(|(r, c)| {
            if r.is_constant() {
                if let Some(v) = numeric::eval_const(r) {
                    let t = v * c.as_f64().unwrap();
                    exact_im &= t.im == 0.0;
                    z += t;
                    return false;
                }
            }
            true
        });
        acc = Expr::Float(z.re);
        if !exact_im || matches!(acc_im, Expr::Float(_)) {
            acc_im = Expr::Float(z.im);
        }
    }
    let mut out: Vec<Expr> = terms.into_iter().filter(|(_, c)| !c.is_zero()).map(|(r, c)| with_coeff(c, r)).collect();
    if !acc_im.is_zero() {
        out.push(gauss(int(0), acc_im));
    }
    let float_zero = matches!(acc, Expr::Float(_));
    if !acc.is_zero() {
        out.push(acc);
    }
    out.sort_by(cmp);
    match out.len() {
        0 => {
            if float_zero {
                Expr::Float(0.0)
            } else {
                int(0)
            }
        }
        1 => out.pop().unwrap(),
        _ => Expr::Add(out),
    }
}

pub fn mul_v(args: Vec<Expr>) -> Expr {
    let mut flat = Vec::with_capacity(args.len());
    for a in args {
        match a {
            Expr::Mul(v) => flat.extend(v),
            x => flat.push(x),
        }
    }
    if flat.iter().any(|e| e.is_undef()) {
        return UNDEF;
    }
    if flat.iter().any(is_struct) {
        return mul_struct(flat);
    }
    let (mut re, mut im) = (int(1), int(0));
    let mut bases: Vec<(Expr, Vec<Expr>, Option<Expr>)> = Vec::new();
    let mut n_inf = 0;
    for f in flat {
        if let Some((a, b)) = as_gauss(&f) {
            let nre = nadd(&nmul(&re, &a), &nneg(&nmul(&im, &b)));
            im = nadd(&nmul(&re, &b), &nmul(&im, &a));
            re = nre;
            continue;
        }
        if f == INF {
            n_inf += 1;
            continue;
        }
        let (b, e) = match &f {
            Expr::Pow(b, e) => ((**b).clone(), (**e).clone()),
            _ => (f.clone(), int(1)),
        };
        if let Some(slot) = bases.iter_mut().find(|(bb, _, _)| *bb == b) {
            slot.1.push(e);
            slot.2 = None;
        } else {
            bases.push((b, vec![e], Some(f)));
        }
    }
    if re.is_zero() && im.is_zero() {
        if n_inf > 0 {
            return UNDEF;
        }
        return if re.as_rat().is_some() && im.as_rat().is_some() && !bases.iter().any(|b| b.0.has_float()) { int(0) } else { Expr::Float(0.0) };
    }
    let mut out: Vec<Expr> = Vec::new();
    let mut again = false;
    for (b, es, orig) in bases {
        let p = match orig {
            Some(f) => f,
            None => {
                let p = pow(b, add_v(es));
                again |= as_gauss(&p).is_some() || matches!(p, Expr::Mul(_)) || p == INF;
                p
            }
        };
        if !p.is_one() {
            out.push(p);
        }
    }
    if again {
        out.push(gauss(re, im));
        out.extend(std::iter::repeat_n(INF, n_inf));
        return mul_v(out);
    }
    if n_inf > 0 {
        if out.is_empty() && im.is_zero() {
            return if re.is_neg_number() { Expr::Mul(vec![int(-1), INF]) } else { INF };
        }
        // sign of a constant cofactor decides the direction
        if im.is_zero() && out.iter().all(|f| f.is_constant()) {
            if let Some(v) = numeric::real_const(&Expr::Mul(out.clone())) {
                let s = v * re.as_f64().unwrap();
                return if s < 0.0 { Expr::Mul(vec![int(-1), INF]) } else { INF };
            }
        }
        out.push(INF);
    }
    if matches!(re, Expr::Float(_)) || matches!(im, Expr::Float(_)) {
        let mut z = C::new(re.as_f64().unwrap(), im.as_f64().unwrap());
        let before = out.len();
        out.retain(|f| {
            if f.is_constant() {
                if let Some(v) = numeric::eval_const(f) {
                    z *= v;
                    return false;
                }
            }
            true
        });
        if out.len() != before {
            let g = from_c(z);
            let (a, b) = as_gauss(&g).unwrap();
            re = a;
            im = b;
        }
    }
    if !im.is_zero() {
        if re.is_zero() {
            out.push(I);
            re = im;
        } else {
            out.push(gauss(re, im));
            re = int(1);
        }
    }
    out.sort_by(cmp);
    if out.is_empty() {
        return re;
    }
    if re.is_one() {
        return if out.len() == 1 { out.pop().unwrap() } else { Expr::Mul(out) };
    }
    out.insert(0, re);
    Expr::Mul(out)
}

fn is_neg_inf(e: &Expr) -> bool {
    matches!(e, Expr::Mul(v) if v.len() == 2 && v[0] == int(-1) && v[1] == INF)
}

pub fn pow(b: Expr, e: Expr) -> Expr {
    if b.is_undef() || e.is_undef() {
        return UNDEF;
    }
    match (&b, &e) {
        (Expr::Matrix(m), _) => {
            return match e.as_i64() {
                Some(n) => matrix::powi(m, n).map(Expr::Matrix).unwrap_or(UNDEF),
                None => UNDEF,
            };
        }
        (Expr::List(_), _) | (_, Expr::List(_)) => return zip_list(&b, &e, &pow),
        (_, Expr::Matrix(_)) => return UNDEF,
        _ => {}
    }
    let inf_like = |v: &Expr| *v == INF || is_neg_inf(v);
    if (inf_like(&e) && b.is_one()) || (inf_like(&b) && e.is_zero()) {
        return UNDEF;
    }
    if e.is_one() {
        return b;
    }
    if e.is_zero() {
        return if matches!(e, Expr::Float(_)) || matches!(b, Expr::Float(_)) { Expr::Float(1.0) } else { int(1) };
    }
    if b.is_one() {
        return b;
    }
    // infinities
    if b == INF {
        if let Some(v) = numeric::real_const(&e) {
            return if v > 0.0 { INF } else { int(0) };
        }
    }
    if e == INF || is_neg_inf(&e) {
        if let Some(v) = numeric::real_const(&b) {
            let grows = (v.abs() > 1.0) == (e == INF);
            if v.abs() == 1.0 || (v < 0.0 && grows) {
                return UNDEF;
            }
            return if grows { INF } else { int(0) };
        }
    }
    if b.is_zero() {
        if let Some(v) = numeric::real_const(&e) {
            return if v > 0.0 { b } else { UNDEF };
        }
        return Expr::Pow(Box::new(b), Box::new(e));
    }
    // floats
    if (matches!(b, Expr::Float(_)) || matches!(e, Expr::Float(_)) || (b.has_float() && as_gauss(&b).is_some()) || (e.has_float() && as_gauss(&e).is_some()))
        && b.is_constant()
        && e.is_constant()
    {
        if let Some(bv) = numeric::eval_const(&b) {
            if let Expr::Num(r) = &e {
                return from_c(numeric::rat_pow(bv, r));
            }
            if let Some(v) = numeric::eval_const(&e) {
                let z = if bv.im == 0.0 && v.im == 0.0 && bv.re > 0.0 { C::new(bv.re.powf(v.re), 0.0) } else { bv.powc(v) };
                return from_c(z);
            }
        }
    }
    match (&b, &e) {
        (Expr::Num(x), Expr::Num(y)) => {
            if y.is_integer() {
                let Some(n) = y.numer().to_i32() else { return Expr::Float(rat_f64(x).powf(rat_f64(y))) };
                if (x.numer().bits() + x.denom().bits()) * n.unsigned_abs() as u64 > 4_000_000 {
                    return Expr::Float(rat_f64(x).powf(n as f64));
                }
                return Expr::Num(x.pow(n));
            }
            return rat_root_pow(x, y);
        }
        (Expr::Const(Const::I), Expr::Num(y)) if y.is_integer() => {
            let k = (y.numer().mod_floor(&BigInt::from(4))).to_i32().unwrap();
            return [int(1), I, int(-1), Expr::Mul(vec![int(-1), I])][k as usize].clone();
        }
        (Expr::Const(Const::E), Expr::Func(n, a)) if n == "ln" => return a[0].clone(),
        // e^(r*pi*i) with a table angle
        (Expr::Const(Const::E), Expr::Mul(v)) if v.contains(&I) && v.contains(&PI) && v.iter().all(|f| *f == I || *f == PI || f.as_rat().is_some()) => {
            let r = v.iter().find_map(|f| f.as_rat().cloned()).unwrap_or_else(Rat::one);
            let half = Rat::new(BigInt::one(), BigInt::from(2));
            if let (Some(s), Some(c)) = (sin_pi(&r), sin_pi(&(&r + &half))) {
                return add(c, mul(I, s));
            }
        }
        (Expr::Const(Const::E), Expr::Mul(v)) if v.len() == 2 && v[0].is_number() && matches!(&v[1], Expr::Func(n, _) if n == "ln") => {
            if let Expr::Func(_, a) = &v[1] {
                return pow(a[0].clone(), v[0].clone());
            }
        }
        (Expr::Pow(b2, e2), _) => {
            if e.as_int().is_some() {
                return pow((**b2).clone(), mul((**e2).clone(), e));
            }
            // (u^n)^(1/n) for even n is |u|
            if let (Some(n), Some(r)) = (e2.as_int(), e.as_rat()) {
                if n.is_even() && (Rat::from_integer(n.clone()) * r).is_one() {
                    return f1("abs", (**b2).clone());
                }
            }
            // (u^a)^b with a,b rational and u a positive number keeps combining
            if b2.as_rat().is_some_and(|r| r.is_positive()) && e2.as_rat().is_some() && e.as_rat().is_some() {
                return pow((**b2).clone(), mul((**e2).clone(), e));
            }
        }
        (Expr::Mul(v), _) => {
            if e.as_int().is_some() {
                return mul_v(v.iter().map(|f| pow(f.clone(), e.clone())).collect());
            }
            // pull positive numeric factors and positive radicals out of the base
            let positive = |f: &Expr| match f {
                Expr::Num(r) => r.is_positive(),
                Expr::Float(x) => *x > 0.0,
                Expr::Pow(bb, _) => bb.as_rat().is_some_and(|r| r.is_positive()),
                Expr::Const(Const::Pi) | Expr::Const(Const::E) => true,
                _ => false,
            };
            if v.iter().any(positive) {
                let (pos, rest): (Vec<Expr>, Vec<Expr>) = v.iter().cloned().partition(positive);
                let mut parts: Vec<Expr> = pos.into_iter().map(|f| pow(f, e.clone())).collect();
                if !rest.is_empty() {
                    let r = if rest.len() == 1 { rest[0].clone() } else { Expr::Mul(rest) };
                    parts.push(match &r {
                        Expr::Num(_) | Expr::Float(_) => pow(r, e.clone()),
                        _ => Expr::Pow(Box::new(r), Box::new(e.clone())),
                    });
                }
                return mul_v(parts);
            }
        }
        _ => {}
    }
    // complex literal to an integer power
    if let (Some((re, im)), Some(n)) = (as_gauss(&b), e.as_i64()) {
        if !im.is_zero() && n.abs() <= 4096 {
            let mut base = b.clone();
            if n < 0 {
                let d = nadd(&nmul(&re, &re), &nmul(&im, &im));
                let dinv = match &d {
                    Expr::Num(r) => Expr::Num(r.recip()),
                    _ => Expr::Float(1.0 / d.as_f64().unwrap()),
                };
                base = gauss(nmul(&re, &dinv), nmul(&nneg(&im), &dinv));
            }
            let (mut acc, mut k) = (int(1), n.unsigned_abs());
            while k > 0 {
                if k & 1 == 1 {
                    acc = mul(acc, base.clone());
                }
                base = mul(base.clone(), base);
                k >>= 1;
            }
            return acc;
        }
    }
    Expr::Pow(Box::new(b), Box::new(e))
}

fn pi_multiple(e: &Expr) -> Option<Rat> {
    match e {
        Expr::Const(Const::Pi) => Some(Rat::one()),
        Expr::Num(r) if r.is_zero() => Some(Rat::zero()),
        Expr::Mul(v) if v.len() == 2 && v[1] == PI => v[0].as_rat().cloned(),
        _ => None,
    }
}

fn sin_pi(r: &Rat) -> Option<Expr> {
    let two = Rat::from_integer(BigInt::from(2));
    let one = Rat::one();
    let half = Rat::new(BigInt::one(), BigInt::from(2));
    let mut r = r - (r / &two).floor() * &two;
    if r >= one {
        return sin_pi(&(r - one)).map(neg);
    }
    if r > half {
        r = &one - r;
    }
    let (n, d) = (r.numer().to_i64()?, r.denom().to_i64()?);
    Some(match (n, d) {
        (0, _) => int(0),
        (1, 2) => int(1),
        (1, 6) => rat(1, 2),
        (1, 4) => div(sqrt(int(2)), int(2)),
        (1, 3) => div(sqrt(int(3)), int(2)),
        (1, 12) => div(sub(sqrt(int(6)), sqrt(int(2))), int(4)),
        (5, 12) => div(add(sqrt(int(6)), sqrt(int(2))), int(4)),
        (1, 10) => div(sub(sqrt(int(5)), int(1)), int(4)),
        (3, 10) => div(add(sqrt(int(5)), int(1)), int(4)),
        (1, 5) => div(sqrt(sub(int(10), mul(int(2), sqrt(int(5))))), int(4)),
        (2, 5) => div(sqrt(add(int(10), mul(int(2), sqrt(int(5))))), int(4)),
        _ => return None,
    })
}

fn half_pi() -> Expr {
    mul(rat(1, 2), PI)
}

/// asin of an exact value, when it is a table entry.
fn asin_exact(x: &Expr) -> Option<Expr> {
    if x.has_neg_coeff() {
        return asin_exact(&neg(x.clone())).map(neg);
    }
    for (n, d) in [(0, 1), (1, 12), (1, 10), (1, 6), (1, 5), (1, 4), (3, 10), (1, 3), (2, 5), (5, 12), (1, 2)] {
        let r = Rat::new(BigInt::from(n), BigInt::from(d));
        if sin_pi(&r).as_ref() == Some(x) {
            return Some(mul(Expr::Num(r), PI));
        }
    }
    None
}

fn perfect_power(n: &BigInt) -> Option<(BigInt, u32)> {
    if n.bits() > 4096 {
        return None;
    }
    for k in (2..=n.bits() as u32).rev() {
        let r = n.nth_root(k);
        if r > BigInt::one() && r.pow(k) == *n {
            return Some((r, k));
        }
    }
    None
}

fn factorial(n: u64) -> BigInt {
    (1..=n).fold(BigInt::one(), |p, k| p * k)
}

pub const UNARY_MATH: &[&str] = &[
    "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh", "ln", "abs", "sign", "floor", "ceil", "ipart", "fpart",
    "factorial", "gamma", "erf", "re", "im", "conj", "arg", "sqrt", "exp", "sec", "csc", "cot", "log", "cbrt", "asec", "acsc", "acot", "round", "erfc",
];

pub fn func(name: &str, args: Vec<Expr>) -> Expr {
    if args.iter().any(|a| a.is_undef()) {
        return UNDEF;
    }
    if args.len() == 1 && UNARY_MATH.contains(&name) {
        match &args[0] {
            Expr::List(v) => return Expr::List(v.iter().map(|x| f1(name, x.clone())).collect()),
            Expr::Matrix(m) => return Expr::Matrix(m.iter().map(|r| r.iter().map(|x| f1(name, x.clone())).collect()).collect()),
            _ => {}
        }
    }
    let x = args.first().cloned().unwrap_or(int(0));
    // rewrites into canonical primitives
    match (name, args.len()) {
        ("sqrt", 1) => return sqrt(x),
        ("cbrt", 1) => return pow(x, rat(1, 3)),
        ("root", 2) => return pow(x, pow(args[1].clone(), int(-1))),
        ("exp", 1) => return pow(E, x),
        ("log", 1) => return div(f1("ln", x), f1("ln", int(10))),
        ("log2", 1) => return div(f1("ln", x), f1("ln", int(2))),
        ("log", 2) => return div(f1("ln", x), f1("ln", args[1].clone())),
        ("sec", 1) => return pow(f1("cos", x), int(-1)),
        ("csc", 1) => return pow(f1("sin", x), int(-1)),
        ("cot", 1) => return div(f1("cos", x.clone()), f1("sin", x)),
        ("asec", 1) => return f1("acos", pow(x, int(-1))),
        ("acsc", 1) => return f1("asin", pow(x, int(-1))),
        ("acot", 1) => return f1("atan", pow(x, int(-1))),
        ("erfc", 1) => return sub(int(1), f1("erf", x)),
        ("int", 1) => return f1("floor", x),
        ("frac", 1) => return f1("fpart", x),
        _ => {}
    }
    // floating-point evaluation
    if args.iter().any(|a| a.has_float()) && args.iter().all(|a| a.is_constant()) {
        let vals: Option<Vec<C>> = args.iter().map(numeric::eval_const).collect();
        if let Some(z) = vals.and_then(|v| numeric::call(name, &v)) {
            return from_c(z);
        }
    }
    let keep = || Expr::Func(name.to_string(), args.clone());
    let odd = |x: &Expr| neg(f1(name, neg(x.clone())));
    match (name, args.len()) {
        ("sin" | "cos" | "tan", 1) => {
            if x == INF || is_neg_inf(&x) {
                return UNDEF;
            }
            if let Some(r) = pi_multiple(&x) {
                let half = Rat::new(BigInt::one(), BigInt::from(2));
                let (s, c) = (sin_pi(&r), sin_pi(&(&r + &half)));
                match name {
                    "sin" => {
                        if let Some(s) = s {
                            return s;
                        }
                    }
                    "cos" => {
                        if let Some(c) = c {
                            return c;
                        }
                    }
                    _ => {
                        if let (Some(s), Some(c)) = (s, c) {
                            if [1, 2, 3, 4, 6].contains(&r.denom().to_i64().unwrap_or(0)) {
                                return div(s, c);
                            }
                        }
                    }
                }
            }
            if let Expr::Func(inner, a) = &x {
                if inner.len() == 4 && inner[1..] == name[..] && inner.starts_with('a') {
                    return a[0].clone();
                }
            }
            if x.has_neg_coeff() {
                return if name == "cos" { f1("cos", neg(x)) } else { odd(&x) };
            }
        }
        ("asin", 1) => {
            if let Some(v) = asin_exact(&x) {
                return v;
            }
            if x.has_neg_coeff() {
                return odd(&x);
            }
        }
        ("acos", 1) => {
            if let Some(v) = asin_exact(&x) {
                return sub(half_pi(), v);
            }
        }
        ("atan", 1) => {
            if x == INF {
                return half_pi();
            }
            if is_neg_inf(&x) {
                return neg(half_pi());
            }
            if x.has_neg_coeff() {
                return odd(&x);
            }
            for (v, n, d) in [(int(0), 0, 1), (div(sqrt(int(3)), int(3)), 1, 6), (int(1), 1, 4), (sqrt(int(3)), 1, 3)] {
                if x == v {
                    return mul(rat(n, d), PI);
                }
            }
        }
        ("sinh" | "tanh" | "asinh" | "atanh", 1) => {
            if x.is_zero() {
                return int(0);
            }
            if x.has_neg_coeff() {
                return odd(&x);
            }
        }
        ("cosh", 1) => {
            if x.is_zero() {
                return int(1);
            }
            if x.has_neg_coeff() {
                return f1("cosh", neg(x));
            }
        }
        ("acosh", 1) => {
            if x.is_one() {
                return int(0);
            }
        }
        ("ln", 1) => {
            if x == E {
                return int(1);
            }
            if x == INF {
                return INF;
            }
            if x == I {
                return mul(rat(1, 2), mul(I, PI));
            }
            match &x {
                Expr::Num(r) if r.is_zero() => return neg(INF),
                Expr::Num(r) if r.is_one() => return int(0),
                Expr::Num(r) if r.is_negative() => return add(f1("ln", Expr::Num(-r)), mul(I, PI)),
                Expr::Num(r) if r.numer().is_one() => return neg(f1("ln", big(r.denom().clone()))),
                Expr::Num(r) if r.is_integer() => {
                    if let Some((root, k)) = perfect_power(r.numer()) {
                        return mul(int(k as i64), f1("ln", big(root)));
                    }
                }
                Expr::Pow(b, e) if **b == E => return (**e).clone(),
                Expr::Pow(b, e) if b.as_rat().is_some_and(|r| r.is_positive()) && e.as_rat().is_some() => {
                    return mul((**e).clone(), f1("ln", (**b).clone()));
                }
                _ => {}
            }
        }
        ("abs", 1) => {
            if let Some((re, im)) = as_gauss(&x) {
                return if im.is_zero() {
                    if re.is_neg_number() { neg(re) } else { re }
                } else {
                    sqrt(nadd(&nmul(&re, &re), &nmul(&im, &im)))
                };
            }
            if x == INF || is_neg_inf(&x) {
                return INF;
            }
            if x.is_constant() {
                if let Some(v) = numeric::real_const(&x) {
                    return if v < 0.0 { neg(x) } else { x };
                }
            }
            if matches!(&x, Expr::Func(n, _) if n == "abs") {
                return x;
            }
            if matches!(&x, Expr::Pow(_, e) if e.as_int().is_some_and(|n| n.is_even())) {
                return x;
            }
            // a positive constant to a (real) power is positive
            if matches!(&x, Expr::Pow(b, _) if matches!(**b, Expr::Const(Const::E) | Expr::Const(Const::Pi)) || b.as_rat().is_some_and(|r| r.is_positive())) {
                return x;
            }
            if x.has_neg_coeff() {
                return f1("abs", neg(x));
            }
        }
        ("sign", 1) => {
            if x.is_constant() {
                if let Some(v) = numeric::real_const(&x) {
                    return int(if v > 0.0 { 1 } else if v < 0.0 { -1 } else { 0 });
                }
            }
        }
        ("floor" | "ceil" | "ipart" | "fpart" | "round", 1) => {
            if let Expr::Num(r) = &x {
                return Expr::Num(match name {
                    "floor" => r.floor(),
                    "ceil" => r.ceil(),
                    "ipart" => r.trunc(),
                    "fpart" => r.fract(),
                    _ => r.round(),
                });
            }
            if x.is_constant() {
                if let Some(v) = numeric::real_const(&x) {
                    if v.abs() < 9e15 {
                        let w = match name {
                            "floor" => v.floor(),
                            "ceil" => v.ceil(),
                            "ipart" => v.trunc(),
                            "round" => v.round(),
                            _ => return sub(x.clone(), f1("ipart", x)),
                        };
                        return int(w as i64);
                    }
                }
            }
        }
        ("round", 2) => {
            if let (Expr::Num(r), Some(d)) = (&x, args[1].as_i64()) {
                let p = Rat::from_integer(BigInt::from(10)).pow(d as i32);
                return Expr::Num((r * &p).round() / p);
            }
        }
        ("factorial", 1) => {
            if let Some(n) = x.as_i64() {
                if (0..=20000).contains(&n) {
                    return big(factorial(n as u64));
                }
                if n < 0 {
                    return UNDEF;
                }
            }
        }
        ("gamma", 1) => {
            if let Some(n) = x.as_i64() {
                if (1..=20000).contains(&n) {
                    return big(factorial(n as u64 - 1));
                }
                return UNDEF;
            }
            // half-integers: gamma(n + 1/2) is a rational multiple of sqrt(pi)
            if let Expr::Num(r) = &x {
                if *r.denom() == BigInt::from(2) {
                    if let Some(n2) = r.numer().to_i64().filter(|n| n.abs() < 400) {
                        let mut coef = Rat::one();
                        let mut k = 1;
                        while k < n2 {
                            coef *= Rat::new(BigInt::from(k), BigInt::from(2));
                            k += 2;
                        }
                        let mut k = n2;
                        while k < 1 {
                            coef /= Rat::new(BigInt::from(k), BigInt::from(2));
                            k += 2;
                        }
                        return mul(Expr::Num(coef), sqrt(PI));
                    }
                }
            }
        }
        ("erf", 1) => {
            if x.is_zero() {
                return int(0);
            }
            if x == INF {
                return int(1);
            }
            if x.has_neg_coeff() {
                return odd(&x);
            }
        }
        ("ncr" | "npr", 2) => {
            if let (Some(n), Some(k)) = (x.as_i64(), args[1].as_i64()) {
                if n >= 0 && k >= 0 && n <= 100000 {
                    if k > n {
                        return int(0);
                    }
                    let k2 = if name == "ncr" { k.min(n - k) } else { k };
                    let mut p = BigInt::one();
                    for j in 0..k2 {
                        p *= n - j;
                    }
                    if name == "ncr" {
                        p /= factorial(k2 as u64);
                    }
                    return big(p);
                }
            }
        }
        ("gcd" | "lcm", 2) => {
            if let (Some(a), Some(b)) = (x.as_int(), args[1].as_int()) {
                return big(if name == "gcd" { a.gcd(b) } else { a.lcm(b) });
            }
        }
        ("mod", 2) => {
            if let (Expr::Num(a), Expr::Num(b)) = (&x, &args[1]) {
                if b.is_zero() {
                    return x;
                }
                return Expr::Num(a - b * (a / b).floor());
            }
        }
        ("min" | "max", _) => {
            let items: Vec<Expr> = if let (1, Expr::List(v)) = (args.len(), &x) { v.clone() } else { args.clone() };
            let vals: Option<Vec<f64>> = items.iter().map(|a| if a.is_constant() { numeric::real_const(a) } else { None }).collect();
            if let Some(vals) = vals {
                if !vals.is_empty() {
                    let mut best = 0;
                    for (i, v) in vals.iter().enumerate() {
                        if (name == "min" && *v < vals[best]) || (name == "max" && *v > vals[best]) {
                            best = i;
                        }
                    }
                    return items[best].clone();
                }
            }
        }
        ("re" | "im" | "conj" | "arg", 1) => {
            let g = as_gauss(&x).or_else(|| if x.is_constant() { numeric::eval_const(&x).filter(|z| z.im == 0.0).map(|_| (x.clone(), int(0))) } else { None });
            if let Some((re, im)) = g {
                match name {
                    "re" => return re,
                    "im" => return im,
                    "conj" => return add(re, mul(neg(im), I)),
                    _ => {
                        let z = numeric::eval_const(&x).unwrap_or(C::new(f64::NAN, 0.0));
                        if z.norm() == 0.0 {
                            return UNDEF;
                        }
                        if im.is_zero() {
                            return if z.re > 0.0 { int(0) } else { PI };
                        }
                        if re.is_zero() {
                            return if z.im > 0.0 { half_pi() } else { neg(half_pi()) };
                        }
                        let base = f1("atan", div(im, re));
                        return if z.re > 0.0 {
                            base
                        } else if z.im > 0.0 {
                            add(base, PI)
                        } else {
                            sub(base, PI)
                        };
                    }
                }
            }
        }
        ("and" | "or", _) => {
            let mut flat: Vec<Expr> = Vec::new();
            for a in &args {
                match a {
                    Expr::Func(n, v) if n == name => flat.extend(v.iter().cloned()),
                    Expr::Bool(b) => {
                        if *b != (name == "and") {
                            return Expr::Bool(*b);
                        }
                    }
                    a => {
                        if !flat.contains(a) {
                            flat.push(a.clone())
                        }
                    }
                }
            }
            return match flat.len() {
                0 => Expr::Bool(name == "and"),
                1 => flat.pop().unwrap(),
                _ => Expr::Func(name.to_string(), flat),
            };
        }
        _ => {}
    }
    keep()
}

/// Canonicalize a raw (parser-built) tree.
pub fn simp(e: &Expr) -> Expr {
    match e {
        Expr::Add(v) => add_v(v.iter().map(simp).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(simp).collect()),
        Expr::Pow(b, x) => pow(simp(b), simp(x)),
        Expr::Func(n, v) => func(n, v.iter().map(simp).collect()),
        other => other.map(&mut |c| simp(c)),
    }
}

/// Structural replacement of every occurrence of `from` by `to`, then re-simplified.
pub fn replace(e: &Expr, from: &Expr, to: &Expr) -> Expr {
    if e == from {
        return to.clone();
    }
    match e {
        Expr::Add(v) => add_v(v.iter().map(|c| replace(c, from, to)).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(|c| replace(c, from, to)).collect()),
        Expr::Pow(b, x) => pow(replace(b, from, to), replace(x, from, to)),
        Expr::Func(n, v) => func(n, v.iter().map(|c| replace(c, from, to)).collect()),
        other => other.map(&mut |c| replace(c, from, to)),
    }
}

pub fn subst(e: &Expr, name: &str, val: &Expr) -> Expr {
    if !e.has_sym(name) {
        return e.clone();
    }
    replace(e, &sym(name), val)
}

fn expand_mul(a: &Expr, b: &Expr) -> Expr {
    let ta: Vec<Expr> = if let Expr::Add(v) = a { v.clone() } else { vec![a.clone()] };
    let tb: Vec<Expr> = if let Expr::Add(v) = b { v.clone() } else { vec![b.clone()] };
    let mut out = Vec::with_capacity(ta.len() * tb.len());
    for x in &ta {
        for y in &tb {
            out.push(mul(x.clone(), y.clone()));
        }
    }
    add_v(out)
}

pub fn expand(e: &Expr) -> Expr {
    match e {
        Expr::Add(v) => add_v(v.iter().map(expand).collect()),
        Expr::Mul(v) => {
            let mut acc = int(1);
            for f in v {
                let f = expand(f);
                if acc.size() * f.size() > 2_000_000 {
                    return e.clone();
                }
                acc = expand_mul(&acc, &f);
            }
            acc
        }
        Expr::Pow(b, x) => {
            let (b, x) = (expand(b), expand(x));
            if let (Expr::Add(_), Some(n)) = (&b, x.as_i64()) {
                if (2..=200).contains(&n.abs()) {
                    let mut acc = b.clone();
                    for _ in 1..n.abs() {
                        acc = expand_mul(&acc, &b);
                    }
                    return if n < 0 { pow(acc, int(-1)) } else { acc };
                }
            }
            // split a sum in the exponent of a positive integer power: handled by pow itself
            let p = pow(b, x);
            if matches!(p, Expr::Mul(_)) { expand_products(&p) } else { p }
        }
        Expr::Func(n, v) => func(n, v.iter().map(expand).collect()),
        other => other.map(&mut |c| expand(c)),
    }
}

fn expand_products(e: &Expr) -> Expr {
    if let Expr::Mul(v) = e {
        let mut acc = int(1);
        for f in v {
            acc = expand_mul(&acc, f);
        }
        acc
    } else {
        e.clone()
    }
}

/// Replace exact numbers and constants by floats.
pub fn to_float(e: &Expr) -> Expr {
    if e.is_constant() {
        if let Some(z) = numeric::eval_const(e) {
            return from_c(z);
        }
    }
    match e {
        Expr::Num(r) => Expr::Float(rat_f64(r)),
        Expr::Const(Const::Pi) => Expr::Float(std::f64::consts::PI),
        Expr::Const(Const::E) => Expr::Float(std::f64::consts::E),
        Expr::Add(v) => add_v(v.iter().map(to_float).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(to_float).collect()),
        Expr::Pow(b, x) => {
            let fb = to_float(b);
            if x.as_rat().is_some() && !fb.is_constant() {
                // keep exact exponents on symbolic bases: x^2 stays x^2
                if x.as_int().is_some() {
                    return pow(fb, (**x).clone());
                }
                return pow(fb, to_float(x));
            }
            if let (Some(z), Expr::Num(r)) = (if fb.is_constant() { numeric::eval_const(&fb) } else { None }, &**x) {
                return from_c(numeric::rat_pow(z, r));
            }
            pow(fb, to_float(x))
        }
        Expr::Func(n, v) => func(n, v.iter().map(to_float).collect()),
        Expr::Rel(r, a, b) => {
            // keep the solved-for symbol side untouched
            Expr::Rel(*r, Box::new(to_float(a)), Box::new(to_float(b)))
        }
        other => other.map(&mut |c| to_float(c)),
    }
}

pub fn is_zero_num(e: &Expr) -> bool {
    if e.is_zero() {
        return true;
    }
    if e.is_constant() {
        if let Some(z) = numeric::eval_const(e) {
            return z.norm() < 1e-12;
        }
    }
    false
}


/// Replace a generator `u` by `t`, also rewriting b^(k*ex) as t^k when u = b^ex.
pub fn replace_gen(e: &Expr, u: &Expr, t: &Expr) -> Expr {
    if e == u {
        return t.clone();
    }
    if let (Expr::Pow(b, ex), Expr::Pow(ub, uex)) = (e, u) {
        if b == ub {
            if let Some(k) = crate::algebra::ratsimp(&div((**ex).clone(), (**uex).clone())).as_rat() {
                if k.is_integer() {
                    return pow(t.clone(), Expr::Num(k.clone()));
                }
            }
        }
    }
    match e {
        Expr::Add(v) => add_v(v.iter().map(|c| replace_gen(c, u, t)).collect()),
        Expr::Mul(v) => mul_v(v.iter().map(|c| replace_gen(c, u, t)).collect()),
        Expr::Pow(b, x) => pow(replace_gen(b, u, t), replace_gen(x, u, t)),
        Expr::Func(n, v) => func(n, v.iter().map(|c| replace_gen(c, u, t)).collect()),
        e => e.clone(),
    }
}
