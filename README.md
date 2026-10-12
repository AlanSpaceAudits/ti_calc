# TI-Calc CAS

A computer algebra calculator in Rust. Target: TI-89 / Nspire CAS feature set.

- `cas/` is the math engine (library plus the `ticalc-cli` terminal front end).
- `app/` is the desktop program (`ticalc`): history tape, keypad, graphing, function catalog.

## Build and run

```sh
cargo build --release
./target/release/ticalc          # desktop program
./target/release/ticalc-cli      # terminal; type `help` for the function list
```

Release builds for Linux, Windows and the web go to `dist/`:

```sh
./scripts/build-dist.sh
```

This needs `mingw-w64-gcc`, `trunk`, and a rustup toolchain with the `x86_64-pc-windows-gnu` and `wasm32-unknown-unknown` targets. `dist/web/` is a static site; serve it with any web server.

Start-up arguments for the desktop program: `--eval "expr"`, `--plot "f(x)"`, `--tab graph|catalog`.

## Input

| Input | Meaning |
| --- | --- |
| `2x^2 + 3x`, `2(x+1)`, `2pi` | implicit multiplication |
| `1.5`, `1E-7` | a decimal point gives a decimal result; without one the result is exact |
| `x := 5`, `5 -> x`, `f(x) := x^2 + 1` | store a variable, define a function |
| `x^2 + 1 \| x = 3` | evaluate with a substitution |
| `[[1, 2], [3, 4]]`, `[1, 2; 3, 4]` | matrix |
| `{1, 2, 3}` | list |
| `pi`, `e`, `i`, `inf`, `ans` | constants and the last answer |
| `0xFF`, `0b1010` | hexadecimal and binary integers |
| `12°34'56"`, `dms(12, 34, 56)` | angle in degrees, minutes, seconds |
| `todms(x)`, `todd(x)`, `torad(x)` | show an angle as DMS, decimal degrees, or radians |

Modes: angle (RAD / DEG / GRAD), Exact / Approx, Real / Complex, Normal / Sci / Eng notation, Fix digits.

The Catalog tab (or `help` in the terminal) lists every function with its syntax.

## Tests

```sh
cargo test -p cas                                  # engine tests
cargo build -p cas && python3 tests/crosscheck.py  # 118 results compared with sympy
```

## Limits

- Symbolic integration uses rules (table, substitution, parts, partial fractions, trig reduction). An integral outside these rules gives a numeric value when the bounds are numbers.
- Factoring is complete for one variable over the rationals up to moderate size. For several variables it finds common factors and differences of squares.
- Numeric equation search covers -100 to 100. Use `nsolve(eq, x, guess)` outside that range.
- In DEG and GRAD modes, symbolic results show the conversion factor (for example `sin(π*x/180)`).
- Floats are 64-bit. Integers and fractions have no size limit.
