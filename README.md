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

## Figures and demos

A result with a natural picture gets a figure under it in the history tape:

| Input | Figure |
| --- | --- |
| an expression in one variable | its graph |
| `solve(eq, x)`, `solve({eq1, eq2}, {x, y})` | both sides or both curves, with the solutions marked |
| `diff`, `integrate`, `limit`, `taylor` | function with derivative, antiderivative, shaded area, limit point, or Taylor polynomial |
| `desolve` | family of solution curves |
| `hamilton(H, q, p)` | phase portrait |
| `fieldsol`, `fieldcheck` | snapshots of the field at successive times |
| `vecadd`, `det` of a 2 by 2 matrix | parallelogram |
| `inscribedangle`, `tangentchord`, `centralangle`, `chordlen` | circle diagram |
| `sphsss`, `sphsas`, `sphdist`, `navtri`, `eq2hor` | triangle or arc on the sphere |
| `lorentz`, `veladd`, `lorentzt`, `galilean` | factor curve or spacetime diagram |
| `kepler(M, e)` | ellipse with the eccentric anomaly |
| `linreg`, `normalcdf` | data with fit, shaded probability |

The Demos tab (or `demo` in the terminal) has 24 worked examples. Each step shows the explanation, the input to type, the result and the figure.

## Function groups

- Differential equations: `desolve(eq, x, y)` with `y'` and `y''`; first order linear, separable, Bernoulli, exact; second order linear with constant coefficients; initial conditions.
- Mechanics: `lagrange(L, q, t)`, `hamiltonian(L, q)`, `hamilton(H, q, p)`. Write a velocity as `q'`.
- Field equations: `fieldeq()`, `fieldeq(name)`, `fieldsol(name, n)`, `fieldcheck(name, u)` for `wave`, `hooke`, `kepler`, `sound`, `heat`, `schrodinger`, `em`, `laplace`, `poisson`, `diffusion`, `kleingordon`. Boosts: `galilean(expr, v)`, `lorentzt(expr, v)`.
- Relativity: `lorentz`, `beta`, `rapidity`, `timedil`, `lencon`, `veladd`, `relmom`, `relenergy`, `relke`, `restenergy`, `reldoppler`, `schwarzschild`, `gravdil`.
- Circle geometry: `inscribedangle`, `centralangle`, `chordlen`, `tangentchord`, `chordsangle`, `secantsangle`, `sagitta`, `arclength`, `sectorarea`, `segmentarea`, `powerpoint`, `vecadd`.
- Spherical trigonometry: `sphsss`, `sphsas`, `sphasa`, `sphaaa`, `sphexcess`, `spharea`, `sphdist`, `sphbearing`, `navtri`.
- Astronomy: `jd`, `caldate`, `gmst`, `lst`, `hms`, `tohms`, `eq2hor`, `hor2eq`, `eq2ecl`, `ecl2eq`, `obliquity`, `angsep`, `sunpos`, `riseset`, `refraction`, `airmass`, `kepler`, `trueanomaly`, `orbperiod`, `visviva`, `vescape`, `parallaxdist`, `distmod`, `absmag`.

Angles, arcs and right ascension use the current angle mode.

## Sources

- Bennett, G. G. (1982). The calculation of astronomical refraction in marine navigation. *Journal of Navigation, 35*(2), 255-259.
- Kasten, F., & Young, A. T. (1989). Revised optical air mass tables and approximation formula. *Applied Optics, 28*(22), 4735-4738.
- Meeus, J. (1998). *Astronomical algorithms* (2nd ed.). Willmann-Bell.
- Murgueitio Ramirez, S. (2024). Galileo's ship and the relativity principle. *Nous* (forthcoming). arXiv:2408.06439.
- U.S. Naval Observatory & HM Nautical Almanac Office. *The Astronomical Almanac*, Section C: low precision formulas for the Sun.

## Tests

```sh
cargo test -p cas                                  # engine tests
cargo build -p cas && python3 tests/crosscheck.py  # 133 results compared with sympy
```

## Limits

- `desolve` covers the equation types listed above. Other equations return an error; `odesolve` gives a numeric value.
- Astronomy has the Sun only (low precision, about 0.01 degree). There are no Moon, planet or star catalog positions.
- Figures use sample values (usually 1) for symbols other than the plotted variable; the caption lists them.
- Symbolic integration uses rules (table, substitution, parts, partial fractions, trig reduction). An integral outside these rules gives a numeric value when the bounds are numbers.
- Factoring is complete for one variable over the rationals up to moderate size. For several variables it finds common factors and differences of squares.
- Numeric equation search covers -100 to 100. Use `nsolve(eq, x, guess)` outside that range.
- In DEG and GRAD modes, symbolic results show the conversion factor (for example `sin(π*x/180)`).
- Floats are 64-bit. Integers and fractions have no size limit.
