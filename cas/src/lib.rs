//! Computer algebra engine for the calculator.
pub mod algebra;
pub mod calculus;
pub mod eval;
pub mod expr;
pub mod matrix;
pub mod numeric;
pub mod parse;
pub mod poly;
pub mod print;
pub mod simplify;
pub mod solve;
pub mod stats;
pub mod units;

pub use eval::{Angle, CATALOG, Calc, Output, Settings};
pub use expr::Expr;
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
