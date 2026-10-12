//! Line-oriented front end: one expression per line on stdin, result on stdout.
use cas::{Angle, Calc, Notation};
use std::io::{BufRead, IsTerminal, Write};

fn main() {
    let mut calc = Calc::new();
    let tty = std::io::stdin().is_terminal();
    let stdin = std::io::stdin();
    loop {
        if tty {
            print!("> ");
            std::io::stdout().flush().ok();
        }
        let mut line = String::new();
        if stdin.lock().read_line(&mut line).unwrap_or(0) == 0 {
            break;
        }
        let line = line.trim();
        match line {
            "" => continue,
            "quit" | "exit" => break,
            "help" => {
                for (s, d) in cas::CATALOG {
                    println!("{s:32} {d}");
                }
                println!("modes: deg rad grad | exact approx | real complex | normal sci eng | fix N | float");
                continue;
            }
            "deg" => calc.settings.angle = Angle::Deg,
            "rad" => calc.settings.angle = Angle::Rad,
            "grad" => calc.settings.angle = Angle::Grad,
            "exact" => calc.settings.approx = false,
            "approx" => calc.settings.approx = true,
            "real" => calc.settings.complex = false,
            "complex" => calc.settings.complex = true,
            "normal" => calc.settings.fmt.notation = Notation::Normal,
            "sci" => calc.settings.fmt.notation = Notation::Sci,
            "eng" => calc.settings.fmt.notation = Notation::Eng,
            "float" => calc.settings.fmt.fix = None,
            _ if line.starts_with("fix ") => calc.settings.fmt.fix = line[4..].trim().parse().ok(),
            _ => {
                match calc.eval_line(line) {
                    Ok(o) => {
                        let mut s = o.text;
                        if let Some(a) = o.approx {
                            s.push_str(&format!("  ≈ {a}"));
                        }
                        if let Some(n) = o.note {
                            s.push_str(&format!("  ({n})"));
                        }
                        println!("{s}");
                    }
                    Err(e) => println!("error: {e}"),
                }
                continue;
            }
        }
        if tty {
            println!("ok");
        }
    }
}
