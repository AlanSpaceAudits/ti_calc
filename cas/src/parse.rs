//! Tokenizer and recursive-descent parser. Produces a raw (unsimplified) expression tree.
use crate::expr::*;
use num_bigint::BigInt;
use num_traits::Num;

#[derive(Clone, Debug, PartialEq)]
enum Tok {
    Num(Expr),
    Id(String),
    Op(&'static str),
}

const OPS: &[&str] = &[":=", "->", "<=", ">=", "!=", "==", "+", "-", "*", "/", "^", "(", ")", "[", "]", "{", "}", ",", ";", "=", "<", ">", "!", "|", "√", "°", "²", "³", "'", "\""];

fn lex(src: &str) -> Result<Vec<Tok>, String> {
    let s: Vec<char> = src.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < s.len() {
        let ch = s[i];
        if ch.is_whitespace() {
            i += 1;
            continue;
        }
        // based literals
        if ch == '0' && i + 2 < s.len() && matches!(s[i + 1], 'x' | 'b' | 'o') && s[i + 2].is_ascii_hexdigit() {
            let radix = match s[i + 1] {
                'x' => 16,
                'b' => 2,
                _ => 8,
            };
            let start = i + 2;
            let mut j = start;
            while j < s.len() && s[j].is_digit(radix) {
                j += 1;
            }
            if j > start {
                let txt: String = s[start..j].iter().collect();
                out.push(Tok::Num(big(BigInt::from_str_radix(&txt, radix).map_err(|e| e.to_string())?)));
                i = j;
                continue;
            }
        }
        if ch.is_ascii_digit() || (ch == '.' && s.get(i + 1).is_some_and(|c| c.is_ascii_digit())) {
            let start = i;
            let mut float = false;
            while i < s.len() && s[i].is_ascii_digit() {
                i += 1;
            }
            if i < s.len() && s[i] == '.' {
                float = true;
                i += 1;
                while i < s.len() && s[i].is_ascii_digit() {
                    i += 1;
                }
            }
            if i < s.len() && matches!(s[i], 'e' | 'E') {
                let mut j = i + 1;
                if j < s.len() && matches!(s[j], '+' | '-') {
                    j += 1;
                }
                if j < s.len() && s[j].is_ascii_digit() {
                    while j < s.len() && s[j].is_ascii_digit() {
                        j += 1;
                    }
                    float = true;
                    i = j;
                }
            }
            let txt: String = s[start..i].iter().collect();
            out.push(Tok::Num(if float {
                Expr::Float(txt.parse::<f64>().map_err(|_| format!("bad number '{txt}'"))?)
            } else {
                big(txt.parse::<BigInt>().map_err(|_| format!("bad number '{txt}'"))?)
            }));
            continue;
        }
        if ch.is_alphabetic() || ch == '_' || ch == '@' {
            if ch == 'π' {
                out.push(Tok::Id("pi".into()));
                i += 1;
                continue;
            }
            let start = i;
            while i < s.len() && (s[i].is_alphanumeric() || s[i] == '_' || s[i] == '@') && s[i] != 'π' {
                i += 1;
            }
            let id: String = s[start..i].iter().collect();
            out.push(Tok::Id(match id.as_str() {
                "Σ" => "sum".into(),
                "Π" => "product".into(),
                _ => id,
            }));
            continue;
        }
        let mapped: Option<&'static str> = match ch {
            '×' | '·' | '⋅' => Some("*"),
            '÷' => Some("/"),
            '−' | '–' => Some("-"),
            '≤' => Some("<="),
            '≥' => Some(">="),
            '≠' => Some("!="),
            '→' => Some("->"),
            _ => None,
        };
        if let Some(m) = mapped {
            out.push(Tok::Op(m));
            i += 1;
            continue;
        }
        match ch {
            '∞' => {
                out.push(Tok::Id("inf".into()));
                i += 1;
                continue;
            }
            '∫' => {
                out.push(Tok::Id("integrate".into()));
                i += 1;
                continue;
            }
            'Σ' => {
                out.push(Tok::Id("sum".into()));
                i += 1;
                continue;
            }
            _ => {}
        }
        let rest: String = s[i..(i + 2).min(s.len())].iter().collect();
        match OPS.iter().find(|op| rest.starts_with(**op)) {
            Some(op) => {
                out.push(Tok::Op(op));
                i += op.chars().count();
            }
            None => return Err(format!("unexpected character '{ch}'")),
        }
    }
    Ok(out)
}

fn neg_raw(e: Expr) -> Expr {
    match e {
        Expr::Num(r) => Expr::Num(-r),
        Expr::Float(f) => Expr::Float(-f),
        e => Expr::Mul(vec![int(-1), e]),
    }
}

fn alias(name: &str) -> String {
    let l = name.to_lowercase();
    match l.as_str() {
        "arcsin" => "asin",
        "arccos" => "acos",
        "arctan" => "atan",
        "arcsinh" => "asinh",
        "arccosh" => "acosh",
        "arctanh" => "atanh",
        "d" | "derivative" => "diff",
        "integral" => "integrate",
        "lim" => "limit",
        "series" => "taylor",
        "eigvl" => "eigenvalues",
        "eigvc" => "eigenvectors",
        "comdenom" | "together" | "cancel" => "ratsimp",
        "zeroes" => "zeros",
        "inverse" => "inv",
        "ceiling" => "ceil",
        "sgn" => "sign",
        "binomial" | "choose" => "ncr",
        "var" => "variance",
        "stddev" | "std" => "stdev",
        "prod" => "product",
        other => other,
    }
    .to_string()
}

struct Parser<'a> {
    t: Vec<Tok>,
    pos: usize,
    is_func: &'a dyn Fn(&str) -> bool,
}

impl Parser<'_> {
    fn peek(&self) -> Option<&Tok> {
        self.t.get(self.pos)
    }
    fn is_op(&self, op: &str) -> bool {
        matches!(self.peek(), Some(Tok::Op(o)) if *o == op)
    }
    fn is_word(&self, w: &str) -> bool {
        matches!(self.peek(), Some(Tok::Id(s)) if s == w)
    }
    fn eat(&mut self, op: &str) -> bool {
        let hit = self.is_op(op);
        if hit {
            self.pos += 1;
        }
        hit
    }
    fn expect(&mut self, op: &str) -> Result<(), String> {
        if self.eat(op) { Ok(()) } else { Err(format!("expected '{op}'")) }
    }

    fn expr(&mut self) -> Result<Expr, String> {
        let e = self.or_expr()?;
        if self.eat("|") {
            let cond = self.or_expr()?;
            return Ok(Expr::Func("|".into(), vec![e, cond]));
        }
        Ok(e)
    }
    fn or_expr(&mut self) -> Result<Expr, String> {
        let mut v = vec![self.and_expr()?];
        while self.is_word("or") {
            self.pos += 1;
            v.push(self.and_expr()?);
        }
        Ok(if v.len() == 1 { v.pop().unwrap() } else { Expr::Func("or".into(), v) })
    }
    fn and_expr(&mut self) -> Result<Expr, String> {
        let mut v = vec![self.rel()?];
        while self.is_word("and") {
            self.pos += 1;
            v.push(self.rel()?);
        }
        Ok(if v.len() == 1 { v.pop().unwrap() } else { Expr::Func("and".into(), v) })
    }
    fn rel_op(&mut self) -> Option<Rel> {
        let r = match self.peek() {
            Some(Tok::Op("=")) | Some(Tok::Op("==")) => Rel::Eq,
            Some(Tok::Op("!=")) => Rel::Ne,
            Some(Tok::Op("<")) => Rel::Lt,
            Some(Tok::Op("<=")) => Rel::Le,
            Some(Tok::Op(">")) => Rel::Gt,
            Some(Tok::Op(">=")) => Rel::Ge,
            _ => return None,
        };
        self.pos += 1;
        Some(r)
    }
    fn rel(&mut self) -> Result<Expr, String> {
        let first = self.sum()?;
        let mut parts = Vec::new();
        let mut left = first.clone();
        while let Some(r) = self.rel_op() {
            let right = self.sum()?;
            parts.push(Expr::Rel(r, Box::new(left), Box::new(right.clone())));
            left = right;
        }
        Ok(match parts.len() {
            0 => first,
            1 => parts.pop().unwrap(),
            _ => Expr::Func("and".into(), parts),
        })
    }
    fn sum(&mut self) -> Result<Expr, String> {
        let mut v = vec![self.term()?];
        loop {
            if self.eat("+") {
                v.push(self.term()?);
            } else if self.eat("-") {
                let t = self.term()?;
                v.push(neg_raw(t));
            } else {
                break;
            }
        }
        Ok(if v.len() == 1 { v.pop().unwrap() } else { Expr::Add(v) })
    }
    fn starts_atom(&self) -> bool {
        match self.peek() {
            Some(Tok::Num(_)) => true,
            Some(Tok::Id(s)) => s != "and" && s != "or",
            Some(Tok::Op(o)) => matches!(*o, "(" | "[" | "{" | "√"),
            None => false,
        }
    }
    fn term(&mut self) -> Result<Expr, String> {
        let mut v = vec![self.unary()?];
        loop {
            if self.eat("*") {
                v.push(self.unary()?);
            } else if self.eat("/") {
                let d = self.unary()?;
                v.push(Expr::Pow(Box::new(d), Box::new(int(-1))));
            } else if self.starts_atom() {
                // implicit multiplication binds like explicit multiplication
                v.push(self.power()?);
            } else {
                break;
            }
        }
        Ok(if v.len() == 1 { v.pop().unwrap() } else { Expr::Mul(v) })
    }
    fn unary(&mut self) -> Result<Expr, String> {
        if self.eat("-") {
            return Ok(neg_raw(self.unary()?));
        }
        if self.eat("+") {
            return self.unary();
        }
        self.power()
    }
    fn power(&mut self) -> Result<Expr, String> {
        let base = self.postfix()?;
        if self.eat("^") {
            let ex = self.unary()?;
            return Ok(Expr::Pow(Box::new(base), Box::new(ex)));
        }
        Ok(base)
    }
    fn postfix(&mut self) -> Result<Expr, String> {
        let mut e = self.atom()?;
        loop {
            if self.is_op("!") {
                self.pos += 1;
                e = Expr::Func("factorial".into(), vec![e]);
            } else if self.eat("²") {
                e = Expr::Pow(Box::new(e), Box::new(int(2)));
            } else if self.eat("³") {
                e = Expr::Pow(Box::new(e), Box::new(int(3)));
            } else if self.eat("°") {
                // 12°34'56.7": minutes and seconds are optional
                let mut parts = vec![e, int(0), int(0)];
                for (slot, mark) in [(1, "'"), (2, "\"")] {
                    if let (Some(Tok::Num(n)), Some(Tok::Op(m))) = (self.t.get(self.pos), self.t.get(self.pos + 1)) {
                        if *m == mark {
                            parts[slot] = n.clone();
                            self.pos += 2;
                        }
                    }
                }
                e = Expr::Func("dms".into(), parts);
            } else {
                break;
            }
        }
        Ok(e)
    }
    fn args(&mut self, close: &str) -> Result<Vec<Expr>, String> {
        let mut v = Vec::new();
        if self.eat(close) {
            return Ok(v);
        }
        loop {
            v.push(self.expr()?);
            if self.eat(close) {
                return Ok(v);
            }
            self.expect(",")?;
        }
    }
    fn atom(&mut self) -> Result<Expr, String> {
        let Some(tok) = self.peek().cloned() else { return Err("unexpected end of input".into()) };
        self.pos += 1;
        match tok {
            Tok::Num(n) => Ok(n),
            Tok::Id(name) => {
                if self.is_op("(") {
                    let fname = alias(&name);
                    if (self.is_func)(&fname) {
                        self.pos += 1;
                        return Ok(Expr::Func(fname, self.args(")")?));
                    }
                    if (self.is_func)(&name) {
                        self.pos += 1;
                        return Ok(Expr::Func(name, self.args(")")?));
                    }
                }
                Ok(Expr::Sym(name))
            }
            Tok::Op("(") => {
                let e = self.expr()?;
                self.expect(")")?;
                Ok(e)
            }
            Tok::Op("√") => {
                let e = self.postfix()?;
                Ok(Expr::Pow(Box::new(e), Box::new(rat(1, 2))))
            }
            Tok::Op("{") => Ok(Expr::List(self.args("}")?)),
            Tok::Op("[") => {
                // [[a, b], [c, d]]  or  [a, b; c, d]  or a row vector [a, b]
                let mut rows: Vec<Vec<Expr>> = vec![vec![]];
                let mut nested = false;
                if self.eat("]") {
                    return Err("empty matrix".into());
                }
                loop {
                    if self.is_op("[") {
                        nested = true;
                        self.pos += 1;
                        let r = self.args("]")?;
                        rows.last_mut().unwrap().extend(r);
                        if self.eat("]") {
                            break;
                        }
                        self.expect(",")?;
                        rows.push(vec![]);
                        continue;
                    }
                    rows.last_mut().unwrap().push(self.expr()?);
                    if self.eat("]") {
                        break;
                    }
                    if self.eat(";") {
                        rows.push(vec![]);
                    } else {
                        self.expect(",")?;
                    }
                }
                let _ = nested;
                let w = rows[0].len();
                if w == 0 || rows.iter().any(|r| r.len() != w) {
                    return Err("matrix rows must have the same length".into());
                }
                Ok(Expr::Matrix(rows))
            }
            Tok::Op(o) => Err(format!("unexpected '{o}'")),
        }
    }
}

/// Parse one input line. Assignments come back as Func(":=", [target, value]).
pub fn parse(src: &str, is_func: &dyn Fn(&str) -> bool) -> Result<Expr, String> {
    let t = lex(src)?;
    if t.is_empty() {
        return Err("empty input".into());
    }
    // locate a top-level := or ->
    let mut depth = 0;
    let mut split = None;
    for (i, tok) in t.iter().enumerate() {
        match tok {
            Tok::Op("(") | Tok::Op("[") | Tok::Op("{") => depth += 1,
            Tok::Op(")") | Tok::Op("]") | Tok::Op("}") => depth -= 1,
            Tok::Op(o @ (":=" | "->")) if depth == 0 => {
                split = Some((i, *o));
                break;
            }
            _ => {}
        }
    }
    let parse_all = |toks: &[Tok], defining: Option<&str>| -> Result<Expr, String> {
        let allow = |n: &str| is_func(n) || defining == Some(n);
        let mut p = Parser { t: toks.to_vec(), pos: 0, is_func: &allow };
        let e = p.expr()?;
        if p.pos < p.t.len() {
            return Err(match &p.t[p.pos] {
                Tok::Op(o) => format!("unexpected '{o}'"),
                Tok::Id(s) => format!("unexpected '{s}'"),
                Tok::Num(_) => "unexpected number".into(),
            });
        }
        Ok(e)
    };
    let target = |toks: &[Tok]| -> Result<Expr, String> {
        match toks {
            [Tok::Id(n)] => Ok(Expr::Sym(n.clone())),
            [Tok::Id(n), Tok::Op("("), mid @ .., Tok::Op(")")] => {
                let mut params = Vec::new();
                for (k, m) in mid.iter().enumerate() {
                    match (k % 2, m) {
                        (0, Tok::Id(p)) => params.push(Expr::Sym(p.clone())),
                        (1, Tok::Op(",")) => {}
                        _ => return Err("function parameters must be names".into()),
                    }
                }
                Ok(Expr::Func(n.clone(), params))
            }
            _ => Err("assignment target must be a name or f(x)".into()),
        }
    };
    match split {
        Some((i, ":=")) => {
            let lhs = target(&t[..i])?;
            let fname = if let Expr::Func(n, _) = &lhs { Some(n.as_str()) } else { None };
            Ok(Expr::Func(":=".into(), vec![lhs.clone(), parse_all(&t[i + 1..], fname)?]))
        }
        Some((i, _)) => {
            let lhs = target(&t[i + 1..])?;
            Ok(Expr::Func(":=".into(), vec![lhs, parse_all(&t[..i], None)?]))
        }
        None => parse_all(&t, None),
    }
}
