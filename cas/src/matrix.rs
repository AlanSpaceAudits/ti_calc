//! Matrix arithmetic over symbolic entries.
use crate::algebra::ratsimp;
use crate::expr::*;
use crate::simplify::*;

pub type Mat = Vec<Vec<Expr>>;

fn norm(e: Expr) -> Expr {
    if e.is_number() {
        e
    } else if e.is_constant() {
        expand(&e)
    } else {
        ratsimp(&e)
    }
}

pub fn dims(m: &Mat) -> (usize, usize) {
    (m.len(), m.first().map_or(0, |r| r.len()))
}

pub fn add(a: &Mat, b: &Mat) -> Option<Mat> {
    (dims(a) == dims(b)).then(|| a.iter().zip(b).map(|(r, s)| r.iter().zip(s).map(|(x, y)| crate::simplify::add(x.clone(), y.clone())).collect()).collect())
}

pub fn scale(a: &Mat, s: &Expr) -> Mat {
    a.iter().map(|r| r.iter().map(|x| crate::simplify::mul(s.clone(), x.clone())).collect()).collect()
}

pub fn mul(a: &Mat, b: &Mat) -> Option<Mat> {
    let ((n, k), (k2, m)) = (dims(a), dims(b));
    if k != k2 {
        return None;
    }
    Some((0..n).map(|i| (0..m).map(|j| norm(add_v((0..k).map(|t| crate::simplify::mul(a[i][t].clone(), b[t][j].clone())).collect()))).collect()).collect())
}

pub fn identity(n: usize) -> Mat {
    (0..n).map(|i| (0..n).map(|j| int((i == j) as i64)).collect()).collect()
}

pub fn transpose(a: &Mat) -> Mat {
    let (n, m) = dims(a);
    (0..m).map(|j| (0..n).map(|i| a[i][j].clone()).collect()).collect()
}

pub fn powi(a: &Mat, n: i64) -> Option<Mat> {
    let (r, c) = dims(a);
    if r != c {
        return None;
    }
    let mut base = if n < 0 { inverse(a)? } else { a.clone() };
    let mut acc = identity(r);
    let mut k = n.unsigned_abs();
    while k > 0 {
        if k & 1 == 1 {
            acc = mul(&acc, &base)?;
        }
        base = mul(&base, &base)?;
        k >>= 1;
    }
    Some(acc)
}

/// Reduced row echelon form. Returns (matrix, pivot columns, determinant sign/scale factor).
pub fn rref(a: &Mat) -> (Mat, Vec<usize>, Expr) {
    let mut m: Mat = a.iter().map(|r| r.iter().cloned().map(norm).collect()).collect();
    let (n, c) = dims(&m);
    let mut pivots = Vec::new();
    let mut det = int(1);
    let mut row = 0;
    for col in 0..c {
        if row >= n {
            break;
        }
        let Some(p) = (row..n).find(|&i| !is_zero_num(&m[i][col])) else { continue };
        if p != row {
            m.swap(p, row);
            det = neg(det);
        }
        let piv = m[row][col].clone();
        det = crate::simplify::mul(det, piv.clone());
        for j in 0..c {
            m[row][j] = norm(div(m[row][j].clone(), piv.clone()));
        }
        for i in 0..n {
            if i != row && !is_zero_num(&m[i][col]) {
                let f = m[i][col].clone();
                for j in 0..c {
                    m[i][j] = norm(sub(m[i][j].clone(), crate::simplify::mul(f.clone(), m[row][j].clone())));
                }
            }
        }
        pivots.push(col);
        row += 1;
    }
    (m, pivots, det)
}

pub fn det(a: &Mat) -> Option<Expr> {
    let (n, c) = dims(a);
    if n != c {
        return None;
    }
    // cofactor expansion keeps symbolic results polynomial for small matrices
    if n <= 4 && a.iter().flatten().any(|e| !e.is_constant()) {
        return Some(expand(&cofactor_det(a)));
    }
    let (_, piv, d) = rref(a);
    Some(if piv.len() < n { int(0) } else { norm(d) })
}

fn cofactor_det(a: &Mat) -> Expr {
    let n = a.len();
    if n == 0 {
        return int(1);
    }
    if n == 1 {
        return a[0][0].clone();
    }
    let mut terms = Vec::new();
    for j in 0..n {
        if a[0][j].is_zero() {
            continue;
        }
        let minor: Mat = a[1..].iter().map(|r| r.iter().enumerate().filter(|(k, _)| *k != j).map(|(_, e)| e.clone()).collect()).collect();
        let t = crate::simplify::mul(a[0][j].clone(), cofactor_det(&minor));
        terms.push(if j % 2 == 1 { neg(t) } else { t });
    }
    add_v(terms)
}

pub fn inverse(a: &Mat) -> Option<Mat> {
    let (n, c) = dims(a);
    if n != c {
        return None;
    }
    let id = identity(n);
    let aug: Mat = a.iter().zip(&id).map(|(r, i)| r.iter().chain(i).cloned().collect()).collect();
    let (m, piv, _) = rref(&aug);
    if piv.len() < n || piv[n - 1] != n - 1 {
        return None;
    }
    Some(m.into_iter().map(|r| r[n..].to_vec()).collect())
}

/// Basis of the null space, one vector per free column.
pub fn nullspace(a: &Mat) -> Vec<Vec<Expr>> {
    let (m, piv, _) = rref(a);
    let (_, c) = dims(a);
    let mut out = Vec::new();
    for free in (0..c).filter(|j| !piv.contains(j)) {
        let mut v = vec![int(0); c];
        v[free] = int(1);
        for (r, &pc) in piv.iter().enumerate() {
            v[pc] = norm(neg(m[r][free].clone()));
        }
        out.push(v);
    }
    out
}

pub fn trace(a: &Mat) -> Option<Expr> {
    let (n, c) = dims(a);
    (n == c).then(|| add_v((0..n).map(|i| a[i][i].clone()).collect()))
}
