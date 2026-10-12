#!/usr/bin/env python3
"""Check ticalc-cli results against sympy. Run: cargo build -p cas && python3 tests/crosscheck.py"""
import random
import subprocess
import sys
from pathlib import Path

import sympy as sp
from sympy.parsing.sympy_parser import convert_xor, parse_expr, standard_transformations

CLI = Path(__file__).resolve().parent.parent / "target" / "debug" / "ticalc-cli"
x, y, a, b, n, k, t = sp.symbols("x y a b n k t", real=True)
LOCALS = {"e": sp.E, "i": sp.I, "π": sp.pi, "ln": sp.log, "abs": sp.Abs, "sign": sp.sign, "∞": sp.oo,
          "x": x, "y": y, "a": a, "b": b, "n": n, "k": k, "t": t}
TR = standard_transformations + (convert_xor,)


def parse(s):
    return parse_expr(s.replace("π", "pi").replace("∞", "oo"), local_dict=LOCALS | {"pi": sp.pi, "oo": sp.oo}, transformations=TR)


def run(lines):
    out = subprocess.run([CLI], input="\n".join(lines) + "\n", capture_output=True, text=True, timeout=600).stdout.splitlines()
    return [o.split("  [figure")[0].split("  ≈ ")[0].split("  (")[0].strip() for o in out]


def close(p, q, tol=1e-8):
    p, q = complex(p), complex(q)
    return abs(p - q) <= tol * (1 + abs(q))


def same(mine, ref):
    """Numeric equality at random real points."""
    syms = sorted((mine.free_symbols | ref.free_symbols), key=str)
    random.seed(1)
    hits = 0
    for _ in range(8):
        sub = {s: sp.Float(random.uniform(0.3, 1.7)) for s in syms}
        try:
            p, q = sp.N(mine.subs(sub)), sp.N(ref.subs(sub))
            if not (p.is_number and q.is_number):
                return False
            if not close(p, q):
                return False
            hits += 1
        except (TypeError, ValueError, ZeroDivisionError):
            continue
    return hits > 0


# (calculator input, kind, sympy reference)
T = [
    ("2^100", "val", "2**100"), ("1/3+1/6", "val", "Rational(1,2)"), ("sqrt(8)", "val", "sqrt(8)"),
    ("(1+i)^10", "val", "(1+I)**10"), ("(2+3i)/(4-i)", "val", "(2+3*I)/(4-I)"), ("sin(pi/12)", "val", "sin(pi/12)"),
    ("cos(3pi/10)", "val", "cos(3*pi/10)"), ("tan(pi/6)", "val", "tan(pi/6)"), ("acos(-1/2)", "val", "acos(Rational(-1,2))"),
    ("ln(e^5)", "val", "5"), ("log(81,3)", "val", "4"), ("25!", "val", "factorial(25)"), ("ncr(40,17)", "val", "binomial(40,17)"),
    ("gamma(7/2)", "val", "gamma(Rational(7,2))"), ("(-8)^(1/3)", "val", "-2"), ("16^(3/4)", "val", "8"),
    ("sqrt(12)+sqrt(27)", "val", "5*sqrt(3)"), ("1/(sqrt(3)-sqrt(2))", "val", "sqrt(3)+sqrt(2)"), ("e^(i*pi/3)", "val", "exp(I*pi/3)"),
    ("abs(3-4i)", "val", "5"), ("mod(-7,3)", "val", "2"), ("gcd(1071,462)", "val", "21"), ("erf(0.5)", "val", "erf(0.5)"),
    ("gamma(4.5)", "val", "gamma(4.5)"), ("sinh(1.2)", "val", "sinh(1.2)"), ("atan(2.)", "val", "atan(2.0)"),
    ("expand((x+2)^5)", "expr", "expand((x+2)**5)"), ("expand((a+b)^3*(a-b))", "expr", "expand((a+b)**3*(a-b))"),
    ("factor(x^4-16)", "expr", "x**4-16"), ("factor(x^6-1)", "expr", "x**6-1"), ("factor(2x^3+3x^2-8x+3)", "expr", "2*x**3+3*x**2-8*x+3"),
    ("simplify((x^3-1)/(x-1))", "expr", "x**2+x+1"), ("simplify(sin(x)^4-cos(x)^4)", "expr", "sin(x)**4-cos(x)**4"),
    ("ratsimp(1/x+1/(x+1)+1/(x+2))", "expr", "1/x+1/(x+1)+1/(x+2)"), ("partfrac((3x+5)/(x^3-x^2-x+1),x)", "expr", "(3*x+5)/(x**3-x**2-x+1)"),
    ("diff(x^3*e^(2x)*sin(x),x)", "expr", "diff(x**3*exp(2*x)*sin(x),x)"), ("diff(atan(x^2)/ln(x),x)", "expr", "diff(atan(x**2)/log(x),x)"),
    ("diff(x^x,x)", "expr", "diff(x**x,x)"), ("diff(asin(x)*sqrt(1-x^2),x,2)", "expr", "diff(asin(x)*sqrt(1-x**2),x,2)"),
    ("diff(tanh(x)*cosh(2x),x)", "expr", "diff(tanh(x)*cosh(2*x),x)"),
    ("taylor(e^x*cos(x),x,6)", "expr", "series(exp(x)*cos(x),x,0,7).removeO()"), ("taylor(ln(x),x,4,1)", "expr", "series(log(x),x,1,5).removeO()"),
    ("taylor(1/(1+x^2),x,6)", "expr", "series(1/(1+x**2),x,0,7).removeO()"),
    ("sum(k^4,k,1,n)", "expr", "summation(k**4,(k,1,n))"), ("sum(k*(k+1),k,1,n)", "expr", "summation(k*(k+1),(k,1,n))"),
    ("sum(3^k,k,0,n)", "expr", "summation(3**k,(k,0,n))"), ("sum(1/k^6,k,1,inf)", "val", "zeta(6)"), ("sum((2/3)^k,k,2,inf)", "val", "Rational(4,3)"),
    ("sum(1/k^3,k,1,inf)", "val", "zeta(3)"), ("sum(1/(k^2+1),k,0,inf)", "val", "(1+pi*coth(pi))/2"),
    ("limit(sin(x)/x,x,0)", "val", "1"), ("limit((1+3/x)^x,x,inf)", "val", "exp(3)"), ("limit((x^3-8)/(x-2),x,2)", "val", "12"),
    ("limit((1-cos(x))/x^2,x,0)", "val", "Rational(1,2)"), ("limit(x*ln(x),x,0,1)", "val", "0"), ("limit((5x^3+x)/(2x^3-7),x,inf)", "val", "Rational(5,2)"),
    ("limit((e^x-1-x)/x^2,x,0)", "val", "Rational(1,2)"), ("limit(x^(1/x),x,inf)", "val", "1"), ("limit((tan(x)-x)/x^3,x,0)", "val", "Rational(1,3)"),
    ("integrate(x^3*e^x,x)", "int", "x**3*exp(x)"), ("integrate(x^2*cos(3x),x)", "int", "x**2*cos(3*x)"), ("integrate(1/(x^2-5x+6),x)", "int", "1/(x**2-5*x+6)"),
    ("integrate((x^3+1)/(x^2+4),x)", "int", "(x**3+1)/(x**2+4)"), ("integrate(sin(x)^4,x)", "int", "sin(x)**4"), ("integrate(sin(3x)*sin(5x),x)", "int", "sin(3*x)*sin(5*x)"),
    ("integrate(x*ln(x)^2,x)", "int", "x*log(x)**2"), ("integrate(e^(3x)*sin(2x),x)", "int", "exp(3*x)*sin(2*x)"), ("integrate(x/sqrt(x^2+9),x)", "int", "x/sqrt(x**2+9)"),
    ("integrate(1/sqrt(4-x^2),x)", "int", "1/sqrt(4-x**2)"), ("integrate(tan(x)^3,x)", "int", "tan(x)**3"), ("integrate(x*atan(x),x)", "int", "x*atan(x)"),
    ("integrate(cos(x)*e^(sin(x)),x)", "int", "cos(x)*exp(sin(x))"), ("integrate(1/(x*(x+1)^2),x)", "int", "1/(x*(x+1)**2)"), ("integrate(sqrt(3x+2),x)", "int", "sqrt(3*x+2)"),
    ("integrate(x^2*sqrt(x-1),x)", "int", "x**2*sqrt(x-1)"), ("integrate(1/(e^x+e^(-x)),x)", "int", "1/(exp(x)+exp(-x))"), ("integrate(ln(x)/x^2,x)", "int", "log(x)/x**2"),
    ("integrate(sqrt(9-x^2),x)", "int", "sqrt(9-x**2)"), ("integrate(2^x*x,x)", "int", "2**x*x"), ("integrate(1/(x^3-1),x)", "int", "1/(x**3-1)"),
    ("integrate(x^2,x,0,3)", "val", "9"), ("integrate(sin(x)^2,x,0,pi)", "val", "pi/2"), ("integrate(e^(-x^2),x,0,inf)", "val", "sqrt(pi)/2"),
    ("integrate(1/(1+x^2),x,0,1)", "val", "pi/4"), ("integrate(x*e^(-x),x,0,inf)", "val", "1"), ("integrate(sin(x)/x,x,0,pi)", "val", "Si(pi)"),
    ("integrate(ln(x),x,1,e)", "val", "1"), ("integrate(1/sqrt(x),x,0,4)", "val", "4"), ("integrate(sqrt(tan(x)),x,0,1)", "val", "Integral(sqrt(tan(x)),(x,0,1))"),
    ("integrate(e^(-x)*cos(x),x,0,inf)", "val", "Rational(1,2)"),
    ("solve(x^2-7x+10=0,x)", "roots", "x**2-7*x+10"), ("solve(x^3-2x^2-5x+6=0,x)", "roots", "x**3-2*x**2-5*x+6"), ("solve(x^3-4x+1=0,x)", "roots", "x**3-4*x+1"),
    ("solve(x^3+2x-7=0,x)", "roots", "x**3+2*x-7"), ("solve(x^4-3x^3+x-1=0,x)", "roots", "x**4-3*x**3+x-1"), ("solve(2x^4-5x^2+1=0,x)", "roots", "2*x**4-5*x**2+1"),
    ("solve(x^5-3x+1=0,x)", "roots", "x**5-3*x+1"), ("csolve(x^4+x+1=0,x)", "croots", "x**4+x+1"), ("csolve(x^3-2x+5=0,x)", "croots", "x**3-2*x+5"),
    ("csolve(x^6+2x^3+5=0,x)", "croots", "x**6+2*x**3+5"), ("csolve(x^4+3x^3-x^2+2x+7=0,x)", "croots", "x**4+3*x**3-x**2+2*x+7"),
    ("solve(e^(2x)-5e^x+6=0,x)", "sols", "[log(2),log(3)]"), ("solve(2^(x+1)=3^x,x)", "sols", "[log(2)/log(Rational(3,2))]"), ("solve(sqrt(2x+3)=x,x)", "sols", "[3]"),
    ("solve(ln(x)^2-ln(x)-2=0,x)", "sols", "[exp(-1),exp(2)]"), ("solve(1/(x-1)+1/(x+1)=1,x)", "sols", "[1-sqrt(2),1+sqrt(2)]"), ("solve(x*e^x=2,x)", "sols", "[LambertW(2)]"),
    ("solve(abs(3x-2)=7,x)", "sols", "[Rational(-5,3),3]"), ("solve(x^(2/3)=4,x)", "sols", "[-8,8]"), ("solve(cos(x)=x^2,x)", "sols", "[-0.824132312302522,0.824132312302522]"),
    ("det([[2,1,3],[0,-1,4],[5,2,1]])", "val", "Matrix([[2,1,3],[0,-1,4],[5,2,1]]).det()"),
    ("det([[1,2,3,4],[2,0,1,5],[3,1,0,2],[4,5,2,1]])", "val", "Matrix([[1,2,3,4],[2,0,1,5],[3,1,0,2],[4,5,2,1]]).det()"),
    ("normalcdf(-1.96,1.96)", "val", "erf(1.96/sqrt(2))"), ("invnorm(0.9)", "val", "sqrt(2)*erfinv(0.8)"), ("binompdf(20,0.3,7)", "val", "binomial(20,7)*0.3**7*0.7**13"),
    ("desolve(y' + 2y = e^x, x, y)", "ode", "Derivative(Y, x) + 2*Y - exp(x)"), ("desolve(y' = x*y, x, y)", "ode", "Derivative(Y, x) - x*Y"),
    ("desolve(y' + y/x = x^2, x, y)", "ode", "Derivative(Y, x) + Y/x - x**2"), ("desolve(y' = y^2*x, x, y)", "ode", "Derivative(Y, x) - Y**2*x"),
    ("desolve(y'' - 5y' + 6y = e^x, x, y)", "ode", "Derivative(Y, x, 2) - 5*Derivative(Y, x) + 6*Y - exp(x)"),
    ("desolve(y'' + 2y' + 5y = 0, x, y)", "ode", "Derivative(Y, x, 2) + 2*Derivative(Y, x) + 5*Y"),
    ("desolve(y'' + 9y = cos(2x), x, y)", "ode", "Derivative(Y, x, 2) + 9*Y - cos(2*x)"), ("desolve(y'' + 4y' + 4y = x, x, y)", "ode", "Derivative(Y, x, 2) + 4*Derivative(Y, x) + 4*Y - x"),
    ("lorentz(12/13*_c)", "val", "Rational(13,5)"), ("veladd(0.6*_c, 0.7*_c)/_c", "val", "1.3/1.42"), ("jd(1987, 4, 10, 19, 21)", "val", "2446896.30625"),
    ("gmst(2446896.30625)*180/pi", "val", "Mod(280.46061837 + 360.98564736629*(2446896.30625-2451545) + 0.000387933*((2446896.30625-2451545)/36525)**2, 360)"),
    ("sphdist(0.5, 0.2, 1.1, 2.0)", "val", "acos(sin(0.5)*sin(1.1) + cos(0.5)*cos(1.1)*cos(1.8))"),
    ("kepler(1.0, 0.3)", "val", "nsolve(Symbol('E') - 0.3*sin(Symbol('E')) - 1.0, 1.0)"), ("chordlen(3, 1.2)", "val", "6*sin(0.6)"),
    ("tcdf(-1.5,1.5,7)", "val", "1-betainc_regularized(3.5,0.5,0,7/(7+2.25))"), ("chi2cdf(0,3,4)", "val", "lowergamma(2,1.5)/gamma(2)"), ("poissoncdf(2.5,3)", "val", "exp(-2.5)*(1+2.5+2.5**2/2+2.5**3/6)"),
]


def main():
    outs = run([t[0] for t in T])
    fails = 0
    for (inp, kind, ref), out in zip(T, outs):
        try:
            if kind in ("val", "expr"):
                r = parse_expr(ref, local_dict={**vars(sp), "x": x, "y": y, "a": a, "b": b, "n": n, "k": k}) if not isinstance(ref, sp.Basic) else ref
                ok = same(parse(out), sp.sympify(r).doit())
            elif kind == "ode":
                C1, C2 = sp.symbols("C1 C2", real=True)
                sol = parse_expr(out.split(" or ")[0].split("=", 1)[1].replace("^", "**").replace("e**", "E**"), local_dict=LOCALS | {"C1": C1, "C2": C2}, transformations=TR)
                Y = sp.Symbol("Y")
                res = parse_expr(ref, local_dict={**vars(sp), "x": x, "Y": Y})
                res = res.replace(lambda e: isinstance(e, sp.Derivative), lambda e: sp.diff(sol, *e.variables)).subs(Y, sol)
                ok = same(res, sp.Integer(0) * x)
            elif kind == "int":
                f = parse_expr(ref, local_dict={**vars(sp), "x": x})
                ok = same(sp.diff(parse(out), x), f)
            else:
                mine = [] if out == "false" else [sp.N(parse(p.split("=", 1)[1])) for p in out.split(" or ")]
                if kind == "sols":
                    want = [sp.N(v) for v in eval(ref, {**vars(sp)})]
                else:
                    want = [sp.N(v) for v in sp.Poly(parse_expr(ref, local_dict={"x": x}), x).nroots(n=20)]
                    if kind == "roots":
                        want = [v for v in want if abs(sp.im(v)) < 1e-12]
                ok = len(mine) == len(want) and all(any(close(m, w, 1e-7) for w in want) for m in mine)
        except Exception as ex:  # a parse failure is a test failure
            ok = False
            out = f"{out}   [{type(ex).__name__}: {ex}]"
        if not ok:
            fails += 1
            print(f"FAIL {inp}\n     got {out}\n     ref {ref}")
    print(f"{len(T) - fails}/{len(T)} passed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
