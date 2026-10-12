use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::cmp::Ordering;
use std::collections::BTreeSet;

pub type Rat = BigRational;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Const {
    Pi,
    E,
    I,
    Inf,
    Undef,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Rel {
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
}

impl Rel {
    pub fn flip(self) -> Rel {
        match self {
            Rel::Lt => Rel::Gt,
            Rel::Le => Rel::Ge,
            Rel::Gt => Rel::Lt,
            Rel::Ge => Rel::Le,
            r => r,
        }
    }
    pub fn text(self) -> &'static str {
        match self {
            Rel::Eq => "=",
            Rel::Ne => "≠",
            Rel::Lt => "<",
            Rel::Le => "≤",
            Rel::Gt => ">",
            Rel::Ge => "≥",
        }
    }
    pub fn holds(self, a: f64, b: f64) -> bool {
        match self {
            Rel::Eq => a == b,
            Rel::Ne => a != b,
            Rel::Lt => a < b,
            Rel::Le => a <= b,
            Rel::Gt => a > b,
            Rel::Ge => a >= b,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    Num(Rat),
    Float(f64),
    Sym(String),
    Const(Const),
    Add(Vec<Expr>),
    Mul(Vec<Expr>),
    Pow(Box<Expr>, Box<Expr>),
    Func(String, Vec<Expr>),
    Rel(Rel, Box<Expr>, Box<Expr>),
    List(Vec<Expr>),
    Matrix(Vec<Vec<Expr>>),
    Bool(bool),
    Str(String),
}

pub const PI: Expr = Expr::Const(Const::Pi);
pub const E: Expr = Expr::Const(Const::E);
pub const I: Expr = Expr::Const(Const::I);
pub const INF: Expr = Expr::Const(Const::Inf);
pub const UNDEF: Expr = Expr::Const(Const::Undef);

pub fn int(n: i64) -> Expr {
    Expr::Num(Rat::from_integer(BigInt::from(n)))
}
pub fn big(n: BigInt) -> Expr {
    Expr::Num(Rat::from_integer(n))
}
pub fn rat(n: i64, d: i64) -> Expr {
    Expr::Num(Rat::new(BigInt::from(n), BigInt::from(d)))
}
pub fn sym(s: &str) -> Expr {
    Expr::Sym(s.to_string())
}
pub fn rat_f64(r: &Rat) -> f64 {
    r.to_f64().unwrap_or_else(|| {
        // ponytail: huge ratios lose precision here; fine for display-only approximations
        let n = r.numer().to_f64().unwrap_or(f64::NAN);
        let d = r.denom().to_f64().unwrap_or(f64::NAN);
        n / d
    })
}

impl Expr {
    pub fn is_zero(&self) -> bool {
        match self {
            Expr::Num(r) => r.is_zero(),
            Expr::Float(f) => *f == 0.0,
            _ => false,
        }
    }
    pub fn is_one(&self) -> bool {
        matches!(self, Expr::Num(r) if r.is_one())
    }
    pub fn is_number(&self) -> bool {
        matches!(self, Expr::Num(_) | Expr::Float(_))
    }
    pub fn is_undef(&self) -> bool {
        matches!(self, Expr::Const(Const::Undef))
    }
    pub fn as_rat(&self) -> Option<&Rat> {
        match self {
            Expr::Num(r) => Some(r),
            _ => None,
        }
    }
    pub fn as_int(&self) -> Option<&BigInt> {
        match self {
            Expr::Num(r) if r.is_integer() => Some(r.numer()),
            _ => None,
        }
    }
    pub fn as_i64(&self) -> Option<i64> {
        self.as_int().and_then(|n| n.to_i64())
    }
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Expr::Num(r) => Some(rat_f64(r)),
            Expr::Float(f) => Some(*f),
            _ => None,
        }
    }
    pub fn is_neg_number(&self) -> bool {
        match self {
            Expr::Num(r) => r.is_negative(),
            Expr::Float(f) => *f < 0.0,
            _ => false,
        }
    }
    /// True when the leading numeric coefficient is negative (used for odd/even rules and printing).
    pub fn has_neg_coeff(&self) -> bool {
        match self {
            Expr::Mul(v) => v[0].is_neg_number(),
            e => e.is_neg_number(),
        }
    }
    pub fn has_sym(&self, name: &str) -> bool {
        match self {
            Expr::Sym(s) => s == name,
            _ => self.children().iter().any(|c| c.has_sym(name)),
        }
    }
    pub fn has_float(&self) -> bool {
        match self {
            Expr::Float(_) => true,
            _ => self.children().iter().any(|c| c.has_float()),
        }
    }
    pub fn has(&self, pred: &dyn Fn(&Expr) -> bool) -> bool {
        pred(self) || self.children().iter().any(|c| c.has(pred))
    }
    pub fn children(&self) -> Vec<&Expr> {
        match self {
            Expr::Add(v) | Expr::Mul(v) | Expr::Func(_, v) | Expr::List(v) => v.iter().collect(),
            Expr::Pow(a, b) | Expr::Rel(_, a, b) => vec![a, b],
            Expr::Matrix(m) => m.iter().flatten().collect(),
            _ => vec![],
        }
    }
    pub fn free_syms(&self) -> BTreeSet<String> {
        let mut s = BTreeSet::new();
        self.collect_syms(&mut s);
        s
    }
    fn collect_syms(&self, out: &mut BTreeSet<String>) {
        if let Expr::Sym(s) = self {
            out.insert(s.clone());
        }
        for c in self.children() {
            c.collect_syms(out);
        }
    }
    /// No free symbols and no infinities/undef: can be evaluated to a number.
    pub fn is_constant(&self) -> bool {
        !self.has(&|e| {
            matches!(
                e,
                Expr::Sym(_)
                    | Expr::Const(Const::Inf)
                    | Expr::Const(Const::Undef)
                    | Expr::List(_)
                    | Expr::Matrix(_)
                    | Expr::Rel(..)
                    | Expr::Bool(_)
                    | Expr::Str(_)
            )
        })
    }
    /// Rebuild with `f` applied to each child (children only, not self).
    pub fn map(&self, f: &mut dyn FnMut(&Expr) -> Expr) -> Expr {
        match self {
            Expr::Add(v) => Expr::Add(v.iter().map(|e| f(e)).collect()),
            Expr::Mul(v) => Expr::Mul(v.iter().map(|e| f(e)).collect()),
            Expr::Func(n, v) => Expr::Func(n.clone(), v.iter().map(|e| f(e)).collect()),
            Expr::List(v) => Expr::List(v.iter().map(|e| f(e)).collect()),
            Expr::Pow(a, b) => Expr::Pow(Box::new(f(a)), Box::new(f(b))),
            Expr::Rel(r, a, b) => Expr::Rel(*r, Box::new(f(a)), Box::new(f(b))),
            Expr::Matrix(m) => Expr::Matrix(m.iter().map(|r| r.iter().map(|e| f(e)).collect()).collect()),
            e => e.clone(),
        }
    }
    pub fn size(&self) -> usize {
        1 + self.children().iter().map(|c| c.size()).sum::<usize>()
    }
}

fn rank(e: &Expr) -> u8 {
    match e {
        Expr::Num(_) | Expr::Float(_) => 0,
        Expr::Const(_) => 1,
        Expr::Sym(_) => 2,
        Expr::Func(..) => 3,
        Expr::Add(_) => 4,
        Expr::Mul(_) => 5,
        Expr::Pow(..) => 6,
        Expr::Rel(..) => 7,
        Expr::List(_) => 8,
        Expr::Matrix(_) => 9,
        Expr::Bool(_) => 10,
        Expr::Str(_) => 11,
    }
}

fn cmp_slice(a: &[Expr], b: &[Expr]) -> Ordering {
    for (x, y) in a.iter().zip(b) {
        let o = cmp(x, y);
        if o != Ordering::Equal {
            return o;
        }
    }
    a.len().cmp(&b.len())
}

fn cmp_raw(a: &Expr, b: &Expr) -> Ordering {
    let (ra, rb) = (rank(a), rank(b));
    if ra != rb {
        return ra.cmp(&rb);
    }
    match (a, b) {
        (Expr::Num(x), Expr::Num(y)) => x.cmp(y),
        (Expr::Float(x), Expr::Float(y)) => x.total_cmp(y),
        (Expr::Num(_), Expr::Float(_)) => Ordering::Less,
        (Expr::Float(_), Expr::Num(_)) => Ordering::Greater,
        (Expr::Const(x), Expr::Const(y)) => x.cmp(y),
        (Expr::Sym(x), Expr::Sym(y)) => x.cmp(y),
        (Expr::Func(n, x), Expr::Func(m, y)) => n.cmp(m).then_with(|| cmp_slice(x, y)),
        (Expr::Add(x), Expr::Add(y)) | (Expr::Mul(x), Expr::Mul(y)) | (Expr::List(x), Expr::List(y)) => cmp_slice(x, y),
        (Expr::Pow(a1, a2), Expr::Pow(b1, b2)) => cmp(a1, b1).then_with(|| cmp(a2, b2)),
        (Expr::Rel(r, a1, a2), Expr::Rel(s, b1, b2)) => r.cmp(s).then_with(|| cmp(a1, b1)).then_with(|| cmp(a2, b2)),
        (Expr::Matrix(x), Expr::Matrix(y)) => {
            let fx: Vec<Expr> = x.iter().flatten().cloned().collect();
            let fy: Vec<Expr> = y.iter().flatten().cloned().collect();
            cmp_slice(&fx, &fy)
        }
        (Expr::Bool(x), Expr::Bool(y)) => x.cmp(y),
        (Expr::Str(x), Expr::Str(y)) => x.cmp(y),
        _ => Ordering::Equal,
    }
}

/// Canonical total order. Powers sort next to their base so like bases are adjacent.
pub fn cmp(a: &Expr, b: &Expr) -> Ordering {
    let one = int(1);
    let (ba, ea) = match a {
        Expr::Pow(b, e) => (&**b, &**e),
        _ => (a, &one),
    };
    let (bb, eb) = match b {
        Expr::Pow(b, e) => (&**b, &**e),
        _ => (b, &one),
    };
    cmp_raw(ba, bb).then_with(|| cmp_raw(ea, eb))
}
