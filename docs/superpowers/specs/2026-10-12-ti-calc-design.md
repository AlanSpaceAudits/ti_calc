# TI-Calc CAS design

Date: 2026-10-12

## Goal

A calculator and equation solver at the TI-89 / Nspire CAS level, written in Rust, as one program that runs on Linux and Windows and can be compiled for the web.

## Structure

Cargo workspace with two crates.

- `cas`: engine library with no GUI dependency.
  - `expr`: expression tree (exact rationals, floats, symbols, constants, sums, products, powers, functions, relations, lists, matrices).
  - `simplify`: automatic simplification into a canonical form, expand, substitution, float conversion.
  - `poly`: polynomials over the rationals (gcd, square-free factoring, rational roots, Kronecker factoring, numeric roots, partial fractions).
  - `algebra`: common denominator, cancel, simplify, factor, rationalize, trig linearization.
  - `solve`: exact polynomial roots to degree 4, isolation, substitution, numeric fallback, inequalities, linear and nonlinear systems.
  - `calculus`: derivative, antiderivative, definite integral, limit, Taylor series, sums.
  - `matrix`, `stats`, `units`: linear algebra, statistics and distributions, unit conversion and constants.
  - `numeric`: complex float evaluation, special functions, root finding, quadrature.
  - `ode`: symbolic ODE solver. `physics`: relativity, Lagrangian and Hamiltonian mechanics, field-equation library, spherical trigonometry. `astro`: positional astronomy.
  - `figure`: builds a list of plain shapes (lines, arrows, fills, dots, text) for a result; front ends only draw them.
  - `demos`: worked examples as (explanation, input) steps, run in a fresh session.
  - `parse`, `print`: input grammar and output format.
  - `eval`: session state (variables, functions, modes) and command dispatch.
- `app`: eframe/egui desktop program. Tabs: Calculator, Graph, Demos, Catalog.

## Decisions

- Own CAS core in pure Rust. No C libraries, so cross-compilation needs no extra setup.
- Dependencies: `num-bigint`, `num-rational`, `num-complex`, `eframe`, `egui_plot`.
- Exact arithmetic by default. A decimal point in the input or Approx mode gives decimal results.
- Every antiderivative is checked numerically against the integrand before it is returned.
- Every exact polynomial root is checked numerically; a failed check falls back to numeric roots.
- Definite integrals compare the exact value with numeric quadrature and split the interval at poles.

## Verification

- `cargo test -p cas`: engine and mode tests.
- `tests/crosscheck.py`: 133 results compared with sympy (values, expressions, antiderivatives, roots).

## Not built

- Risch integration, multivariate polynomial factoring, partial differential equation solving beyond the library of standard solutions.
- Moon, planet and star positions.
- 2D pretty-printed output.
