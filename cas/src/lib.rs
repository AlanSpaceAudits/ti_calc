//! Computer algebra engine for the calculator.
pub mod algebra;
pub mod astro;
pub mod calculus;
pub mod demos;
pub mod eval;
pub mod expr;
pub mod figure;
pub mod matrix;
pub mod numeric;
pub mod ode;
pub mod parse;
pub mod physics;
pub mod poly;
pub mod print;
pub mod simplify;
pub mod solve;
pub mod stats;
pub mod units;

pub use eval::{Angle, CATALOG, Calc, Output, Settings};
pub use expr::Expr;
pub use figure::{Figure, Shape};
pub use print::{Fmt, Notation};

#[cfg(test)]
mod tests {
    use super::*;

    fn run(src: &str) -> String {
        Calc::new().eval_line(src).map(|o| o.text).unwrap_or_else(|e| format!("error: {e}"))
    }

    #[test]
    fn core_results() {
        for (input, want) in [
            ("1/3+1/6", "1/2"),
            ("sqrt(8)", "2*sqrt(2)"),
            ("(1+i)^2", "2*i"),
            ("sin(pi/6)", "1/2"),
            ("expand((x+1)^3)", "x^3 + 3*x^2 + 3*x + 1"),
            ("factor(x^4+4)", "(x^2 - 2*x + 2)*(x^2 + 2*x + 2)"),
            ("solve(x^2-5x+6=0,x)", "x = 2 or x = 3"),
            ("solve(x^2<4,x)", "-2 < x < 2"),
            ("solve({x+y=3,x-y=1},{x,y})", "x = 2 and y = 1"),
            ("diff(x^3*sin(x),x)", "x^3*cos(x) + 3*x^2*sin(x)"),
            ("integrate(x*e^x,x)", "x*e^x - e^x"),
            ("integrate(e^(-x^2),x,-inf,inf)", "sqrt(π)"),
            ("limit((1+1/x)^x,x,inf)", "e"),
            ("sum(k^2,k,1,n)", "n*(n + 1)*(2*n + 1)/6"),
            ("det([[1,2],[3,4]])", "-2"),
            ("inv([[1,2],[3,4]])", "[[-2, 1], [3/2, -1/2]]"),
            ("mean({1,2,3,4})", "5/2"),
            ("0.1+0.2", "0.3"),
            ("integrate(1/x,x,-1,1)", "error: integral does not converge on this interval"),
        ] {
            assert_eq!(run(input), want, "input: {input}");
        }
    }

    #[test]
    fn demos_run_clean() {
        for d in demos::DEMOS {
            for s in demos::run(d) {
                assert!(s.output.is_ok(), "{} / {}: {:?}", d.title, s.input, s.output);
            }
        }
    }

    #[test]
    fn field_solutions_satisfy_their_equations() {
        for f in physics::FIELD_EQS {
            for (label, u) in (f.solutions)() {
                assert!(physics::field_check(f, &u).is_none(), "{}: {}", f.name, label);
            }
        }
    }

    #[test]
    fn ode_and_mechanics() {
        for (input, want) in [
            ("desolve(y' = y, x, y, 0, 2)", "y = 2*e^x"),
            ("desolve(y'' + y = 0, x, y, 0, 1, 0)", "y = cos(x)"),
            ("desolve(y'' + 4y = sin(x), x, y)", "y = C1*cos(2*x) + C2*sin(2*x) + sin(x)/3"),
            ("lagrange(m*x'^2/2 - k*x^2/2, x, t)", "x'' = -k*x/m"),
            ("hamilton(p^2/(2m) + k*x^2/2, x, p)", "{x' = p/m, p' = -k*x}"),
            ("lorentz(3/5*_c)", "5/4"),
            ("simplify(lorentzt(x^2 - c^2*t^2, v))", "-c^2*t^2 + x^2"),
        ] {
            assert_eq!(run(input), want, "input: {input}");
        }
    }

    #[test]
    fn modes() {
        let mut c = Calc::new();
        c.settings.angle = Angle::Deg;
        assert_eq!(c.eval_line("sin(30)").unwrap().text, "1/2");
        assert_eq!(c.eval_line("atan(1)").unwrap().text, "45");
        c.settings.fmt.notation = Notation::Sci;
        assert_eq!(c.eval_line("123456.789").unwrap().text, "1.23456789E5");
        c.settings.complex = true;
        assert_eq!(c.eval_line("sqrt(-4)").unwrap().text, "2*i");
    }
}
