//! Worked examples: each step pairs a short explanation with the exact input that produces it.
use crate::eval::{Angle, Calc};

pub struct Demo {
    pub group: &'static str,
    pub title: &'static str,
    pub about: &'static str,
    /// (explanation, input). The inputs "deg" and "rad" switch the angle mode.
    pub steps: &'static [(&'static str, &'static str)],
}

pub const DEMOS: &[Demo] = &[
    Demo {
        group: "Algebra",
        title: "Quadratic equation",
        about: "Solve 2x^2 - 4x - 3 = 0 with the quadratic formula, then in one step.",
        steps: &[
            ("Store the coefficients a, b, c.", "a := 2"),
            ("", "b := -4"),
            ("", "c := -3"),
            ("The discriminant b^2 - 4ac decides the number of real roots.", "disc := b^2 - 4*a*c"),
            ("First root: (-b - sqrt(disc)) / (2a).", "(-b - sqrt(disc))/(2*a)"),
            ("Second root: (-b + sqrt(disc)) / (2a).", "(-b + sqrt(disc))/(2*a)"),
            ("solve does the same directly. Write the equation, then the unknown.", "solve(2x^2 - 4x - 3 = 0, x)"),
        ],
    },
    Demo {
        group: "Algebra",
        title: "System of equations",
        about: "Solve x + y = 7 and x^2 + y^2 = 25 by substitution, then in one step.",
        steps: &[
            ("Solve the linear equation for y.", "solve(x + y = 7, y)"),
            ("Substitute y = 7 - x into the second equation and solve for x.", "solve(x^2 + (7 - x)^2 = 25, x)"),
            ("Each x gives a y.", "7 - x | x = 3"),
            ("", "7 - x | x = 4"),
            ("One step: a list of equations and a list of unknowns.", "solve({x + y = 7, x^2 + y^2 = 25}, {x, y})"),
        ],
    },
    Demo {
        group: "Calculus",
        title: "Integration by parts",
        about: "Integrate x*e^x using the rule: integral of u dv = u*v - integral of v du.",
        steps: &[
            ("Choose u = x. Its derivative du:", "du := diff(x, x)"),
            ("Choose dv = e^x dx. Its antiderivative v:", "v := integrate(e^x, x)"),
            ("Apply the rule: u*v minus the integral of v*du.", "x*v - integrate(v*du, x)"),
            ("integrate gives the same result.", "integrate(x*e^x, x)"),
            ("Check: the derivative returns the integrand.", "simplify(diff(x*e^x - e^x, x))"),
        ],
    },
    Demo {
        group: "Calculus",
        title: "Limit, derivative, Taylor series",
        about: "Three views of sin(x) near 0.",
        steps: &[
            ("The limit of sin(x)/x as x goes to 0.", "limit(sin(x)/x, x, 0)"),
            ("The derivative from its definition: limit of (f(x+h) - f(x))/h.", "limit((sin(x + h) - sin(x))/h, h, 0)"),
            ("diff gives it directly.", "diff(sin(x), x)"),
            ("Taylor polynomial of order 5 about 0.", "taylor(sin(x), x, 5)"),
        ],
    },
    Demo {
        group: "Calculus",
        title: "Area under a curve",
        about: "The definite integral as an area, and the antiderivative that measures it.",
        steps: &[
            ("The area under x^2 from 0 to 2.", "integrate(x^2, x, 0, 2)"),
            ("The antiderivative F(x).", "integrate(x^2, x)"),
            ("The area is F(2) - F(0).", "x^3/3 | x = 2"),
            ("A probability is an area too: the normal curve between -1 and 2.", "normalcdf(-1, 2)"),
        ],
    },
    Demo {
        group: "Differential equations",
        title: "First-order linear ODE",
        about: "Solve y' + 2y = e^x with an integrating factor. Standard form: y' + P*y = Q.",
        steps: &[
            ("Here P = 2. The integrating factor is mu = e^(integral of P).", "mu := e^(integrate(2, x))"),
            ("Multiply Q = e^x by mu and integrate.", "w := integrate(mu*e^x, x)"),
            ("The general solution is y = (w + C1)/mu.", "ys := expand((w + C1)/mu)"),
            ("Check by substitution: y' + 2y - e^x must be 0.", "simplify(diff(ys, x) + 2*ys - e^x)"),
            ("desolve does all of it. Write the derivative as y'.", "desolve(y' + 2y = e^x, x, y)"),
            ("With an initial condition y(0) = 1.", "desolve(y' + 2y = e^x, x, y, 0, 1)"),
        ],
    },
    Demo {
        group: "Differential equations",
        title: "Separable ODE",
        about: "Solve y' = x*y. Separate: dy/y = x dx, then integrate both sides.",
        steps: &[
            ("Left side: integral of 1/y dy.", "integrate(1/y, y)"),
            ("Right side: integral of x dx.", "integrate(x, x)"),
            ("So ln|y| = x^2/2 + C. Solve for y.", "solve(ln(y) = x^2/2 + C, y)"),
            ("desolve gives the same family (C1 replaces e^C).", "desolve(y' = x*y, x, y)"),
        ],
    },
    Demo {
        group: "Differential equations",
        title: "Second-order: damped oscillator",
        about: "Solve y'' + 2y' + 5y = 0. Try y = e^(r*x); the equation becomes r^2 + 2r + 5 = 0.",
        steps: &[
            ("Roots of the characteristic equation.", "csolve(r^2 + 2r + 5 = 0, r)"),
            ("Roots -1 +/- 2i give e^(-x) times cos(2x) and sin(2x).", "desolve(y'' + 2y' + 5y = 0, x, y)"),
            ("Initial conditions y(0) = 1, y'(0) = 0.", "desolve(y'' + 2y' + 5y = 0, x, y, 0, 1, 0)"),
            ("With a forcing term the solver adds a particular solution.", "desolve(y'' + 4y = sin(x), x, y)"),
        ],
    },
    Demo {
        group: "Mechanics",
        title: "Euler-Lagrange: mass on a spring",
        about: "L = T - V = m*x'^2/2 - k*x^2/2. The equation of motion is d/dt(dL/dx') - dL/dx = 0.",
        steps: &[
            ("Momentum term dL/dx'. Write the velocity as x'.", "diff(m*x'^2/2 - k*x^2/2, x')"),
            ("Force term dL/dx.", "diff(m*x'^2/2 - k*x^2/2, x)"),
            ("lagrange builds d/dt(dL/dx') - dL/dx = 0 and solves it for x''.", "lagrange(m*x'^2/2 - k*x^2/2, x, t)"),
            ("Solve that equation of motion for x(t).", "desolve(m*x'' = -k*x, t, x)"),
            ("Release from rest at x = A.", "desolve(m*x'' = -k*x, t, x, 0, A, 0)"),
        ],
    },
    Demo {
        group: "Mechanics",
        title: "Euler-Lagrange: pendulum and central force",
        about: "The same recipe for other coordinates.",
        steps: &[
            ("Pendulum of length l, angle q: L = m*l^2*q'^2/2 + m*g*l*cos(q).", "lagrange(m*l^2*q'^2/2 + m*g*l*cos(q), q, t)"),
            ("Projectile: two coordinates in a list.", "lagrange(m*(x'^2 + y'^2)/2 - m*g*y, {x, y}, t)"),
            ("Central force in polar coordinates r, th.", "lagrange(m*(r'^2 + r^2*th'^2)/2 + G*M*m/r, {r, th}, t)"),
        ],
    },
    Demo {
        group: "Mechanics",
        title: "Hamiltonian of the spring",
        about: "Legendre transform: p = dL/dx', H = p*x' - L written with p.",
        steps: &[
            ("The momentum p = dL/dx'.", "diff(m*x'^2/2 - k*x^2/2, x')"),
            ("hamiltonian eliminates the velocity and returns H(x, p).", "H := hamiltonian(m*x'^2/2 - k*x^2/2, x)"),
            ("Hamilton's equations: x' = dH/dp, p' = -dH/dx.", "hamilton(H, x, p)"),
        ],
    },
    Demo {
        group: "Field equations",
        title: "Wave equation on a string",
        about: "d2u/dt2 = c^2 d2u/dx2 holds in the frame in which the string is at rest (Murgueitio Ramirez, 2024).",
        steps: &[
            ("The library entry: equation, frame, symbols, solutions.", "fieldeq(wave)"),
            ("A travelling wave. Write your own trial solution in x and t.", "u1 := sin(x - c*t)"),
            ("Substitute it. Both sides are second derivatives:", "diff(u1, t, 2) - c^2*diff(u1, x, 2)"),
            ("fieldcheck does that substitution for any trial solution.", "fieldcheck(wave, sin(x - c*t))"),
            ("The standing mode of a string of length L.", "fieldsol(wave, 2)"),
            ("", "fieldcheck(wave, fieldsol(wave, 2))"),
            ("The same wave in a frame moving at speed v: x -> x - v*t.", "galilean(sin(x - c*t), v)"),
            ("It no longer satisfies the equation. The residual carries the terms in v.", "fieldcheck(wave, galilean(sin(x - c*t), v))"),
        ],
    },
    Demo {
        group: "Field equations",
        title: "Heat equation",
        about: "du/dt = alpha d2u/dx2 holds in the frame in which the heated material is at rest (Murgueitio Ramirez, 2024).",
        steps: &[
            ("", "fieldeq(heat)"),
            ("Separate variables: try u = e^(-s*t)*sin(k*x) and find s.", "solve(diff(e^(-s*t)*sin(k*x), t) = alpha*diff(e^(-s*t)*sin(k*x), x, 2), s)"),
            ("So each mode decays at the rate alpha*k^2.", "fieldcheck(heat, e^(-alpha*k^2*t)*sin(k*x))"),
            ("A trial that fails returns what is left over.", "fieldcheck(heat, x^2)"),
            ("The heat kernel (a point of heat spreading out).", "fieldcheck(heat, fieldsol(heat, 2))"),
        ],
    },
    Demo {
        group: "Field equations",
        title: "Schrodinger equation: particle in a box",
        about: "i*hbar du/dt = -(hbar^2/2m) d2u/dx2 holds in the frame in which the box is at rest (Murgueitio Ramirez, 2024).",
        steps: &[
            ("", "fieldeq(schrodinger)"),
            ("A stationary state sin(n*pi*x/L)*e^(-i*En*t/hbar). Find the energy En that makes it a solution.", "solve(i*hbar*diff(sin(n*pi*x/L)*e^(-i*En*t/hbar), t) = -hbar^2/(2m)*diff(sin(n*pi*x/L)*e^(-i*En*t/hbar), x, 2), En)"),
            ("The library solution uses that energy.", "fieldcheck(schrodinger, fieldsol(schrodinger, 1))"),
            ("Normalization: the integral of the squared amplitude over the box is 1 (here n = 1).", "integrate((sqrt(2/L)*sin(pi*x/L))^2, x, 0, L)"),
        ],
    },
    Demo {
        group: "Field equations",
        title: "Light: Galilean and Lorentz boosts",
        about: "The electromagnetic wave equation keeps its form under a Lorentz boost and not under a Galilean boost (Murgueitio Ramirez, 2024).",
        steps: &[
            ("", "fieldeq(em)"),
            ("A plane wave satisfies it.", "fieldcheck(em, sin(x - c*t))"),
            ("Galilean boost x -> x - v*t: a residual in v is left.", "fieldcheck(em, galilean(sin(x - c*t), v))"),
            ("Lorentz boost: the equation is still satisfied.", "fieldcheck(em, lorentzt(sin(x - c*t), v))"),
            ("The interval x^2 - c^2 t^2 is unchanged by the Lorentz boost.", "simplify(lorentzt(x^2 - c^2*t^2, v))"),
        ],
    },
    Demo {
        group: "Field equations",
        title: "Kepler orbit, spring, sound",
        about: "The other lab-frame systems in Murgueitio Ramirez (2024). fieldeq() lists every equation in the library.",
        steps: &[
            ("", "fieldeq()"),
            ("Kepler: the conic section solves the orbit equation.", "fieldcheck(kepler, fieldsol(kepler))"),
            ("Hooke's law: oscillation about the equilibrium x0.", "fieldcheck(hooke, x0 + A*cos(sqrt(k/m)*t + phi))"),
            ("Sound: a spherical wave from a point source.", "fieldcheck(sound, fieldsol(sound, 2))"),
        ],
    },
    Demo {
        group: "Relativity",
        title: "Lorentz factor and time dilation",
        about: "With a symbolic speed the result keeps c as a symbol. With a number it uses c = 299792458 m/s (the constant _c).",
        steps: &[
            ("The Lorentz factor.", "lorentz(v)"),
            ("At 60 % of c. Write speeds as a fraction times _c.", "lorentz(3/5*_c)"),
            ("A 10 s interval on a clock moving at 0.8c.", "timedil(10, 0.8*_c)"),
            ("Length contraction.", "lencon(L0, v)"),
            ("Velocity addition: 0.9c and 0.9c, as a fraction of c.", "veladd(0.9*_c, 0.9*_c)/_c"),
            ("Rest energy of 1 kg in joules.", "restenergy(1)"),
            ("Kinetic energy: the low-speed limit is m*v^2/2.", "taylor(relke(m, v), v, 4)"),
        ],
    },
    Demo {
        group: "Geometry",
        title: "Circle theorems",
        about: "Angles in a circle, each with its figure. Angles are in degrees here.",
        steps: &[
            ("Work in degrees.", "deg"),
            ("Inscribed angle theorem: an angle on the circle is half the central angle on the same arc.", "inscribedangle(80)"),
            ("Thales: an angle on a diameter (central angle 180) is a right angle.", "inscribedangle(180)"),
            ("Central angle of a chord of length 1 in a circle of radius 1.", "centralangle(1, 1)"),
            ("The chord from the angle: 2*r*sin(angle/2).", "chordlen(2, 60)"),
            ("Tangent-chord theorem: the angle between a tangent and a chord is half the intercepted arc.", "tangentchord(100)"),
            ("Two chords that cross inside the circle: half the sum of the two arcs.", "chordsangle(80, 40)"),
            ("Two secants that meet outside: half the difference of the arcs.", "secantsangle(80, 40)"),
            ("Sagitta (arc height) of a 6-unit chord in a circle of radius 5.", "sagitta(5, 6)"),
        ],
    },
    Demo {
        group: "Geometry",
        title: "Vectors, parallelograms, determinants",
        about: "Vector addition is the diagonal of a parallelogram; a 2 by 2 determinant is its area.",
        steps: &[
            ("Add two vectors. The figure shows the parallelogram.", "vecadd([3, 1], [1, 2])"),
            ("Length of the sum.", "norm([4, 3])"),
            ("The determinant of the matrix with those columns is the area between them.", "det([[3, 1], [1, 2]])"),
            ("Dot product and the angle between the vectors.", "acos(dot([3, 1], [1, 2])/(norm([3, 1])*norm([1, 2])))"),
            ("Two straight lines as a system: the solution is the crossing point.", "solve({3x + y = 9, x + 2y = 8}, {x, y})"),
        ],
    },
    Demo {
        group: "Astronomy",
        title: "Date, sidereal time, Sun position",
        about: "From a calendar date to the direction of the Sun. Angles in degrees; right ascension is an angle (15 degrees per hour).",
        steps: &[
            ("Work in degrees.", "deg"),
            ("Julian date of 2024 June 20, 20:51 UT.", "d := jd(2024, 6, 20, 20, 51)"),
            ("Greenwich mean sidereal time as an angle.", "gmst(d)"),
            ("The same as hours, minutes, seconds.", "tohms(gmst(d))"),
            ("Local sidereal time at longitude 74 W (east is positive).", "lst(d, -74)"),
            ("The Sun: right ascension, declination, ecliptic longitude, distance.", "sunpos(d)"),
            ("Obliquity of the ecliptic on that date.", "obliquity(d)"),
            ("Ecliptic longitude 90, latitude 0 in equatorial coordinates.", "ecl2eq(90, 0)"),
        ],
    },
    Demo {
        group: "Astronomy",
        title: "Altitude, azimuth, rising and setting",
        about: "The navigational triangle joins the pole P, the zenith Z and the body X.",
        steps: &[
            ("Work in degrees.", "deg"),
            ("Hour angle 30, declination 20, latitude 40: altitude and azimuth.", "eq2hor(30, 20, 40)"),
            ("The altitude from the law of cosines for the side ZX.", "asin(sin(40)*sin(20) + cos(40)*cos(20)*cos(30))"),
            ("Back from altitude and azimuth to hour angle and declination.", "hor2eq(57.4850799244, 240.938807385, 40)"),
            ("Rising and setting of a body at declination 23.44 seen from latitude 40.", "riseset(23.44, 40)"),
            ("Refraction at the horizon, as degrees, minutes, seconds.", "todms(refraction(0))"),
            ("Angular separation of two directions (ra, dec pairs).", "angsep(10, 20, 50, 40)"),
        ],
    },
    Demo {
        group: "Astronomy",
        title: "Kepler's equation and orbits",
        about: "Mean anomaly M grows evenly with time. Kepler's equation E - e*sin(E) = M gives the eccentric anomaly E.",
        steps: &[
            ("Work in degrees.", "deg"),
            ("Eccentric anomaly for M = 60, eccentricity 0.6.", "E1 := kepler(60, 0.6)"),
            ("Check: E - e*sin(E) in degrees returns M.", "E1 - 0.6*sin(E1)*180/pi"),
            ("True anomaly (the angle at the focus).", "trueanomaly(E1, 0.6)"),
            ("Orbital period from the semi-major axis a and mu = G*M.", "orbperiod(a, mu)"),
            ("Kepler's third law: the period of Earth (a = 1 au, mu of the Sun) in days.", "orbperiod(1.495978707E11, 1.32712440018E20)/86400"),
            ("Orbital speed at distance r (vis-viva).", "visviva(mu, r, a)"),
        ],
    },
    Demo {
        group: "Spherical trigonometry",
        title: "Great-circle distance and bearing",
        about: "New York (40.7128 N, 74.0060 W) to London (51.5074 N, 0.1278 W). West longitudes are negative.",
        steps: &[
            ("Work in degrees.", "deg"),
            ("Central angle between the two points.", "sphdist(40.7128, -74.0060, 51.5074, -0.1278)"),
            ("The same angle in degrees, minutes, seconds.", "todms(ans)"),
            ("Give a radius as the fifth argument for the distance (6371 km).", "sphdist(40.7128, -74.0060, 51.5074, -0.1278, 6371)"),
            ("Initial bearing, measured from north toward east.", "sphbearing(40.7128, -74.0060, 51.5074, -0.1278)"),
        ],
    },
    Demo {
        group: "Spherical trigonometry",
        title: "Spherical triangles",
        about: "Sides are arcs measured as angles. Capital letters are the angles opposite the sides.",
        steps: &[
            ("Work in degrees.", "deg"),
            ("Three sides of 90 degrees: one octant of the sphere.", "sphsss(90, 90, 90)"),
            ("Its spherical excess (angle sum minus 180).", "sphexcess(90, 90, 90)"),
            ("Three sides known: law of cosines for sides.", "sphsss(60, 70, 80)"),
            ("Two sides and the included angle.", "sphsas(60, 50, 70)"),
            ("The law of cosines written out for that third side.", "acos(cos(60)*cos(70) + sin(60)*sin(70)*cos(50))"),
            ("Navigational triangle: latitude 40, declination 20, local hour angle 30.", "navtri(40, 20, 30)"),
        ],
    },
];

pub struct StepResult {
    pub note: &'static str,
    pub input: &'static str,
    pub output: Result<String, String>,
    pub figure: Option<crate::figure::Figure>,
}

/// Run a demo in a fresh session.
pub fn run(d: &Demo) -> Vec<StepResult> {
    let mut calc = Calc::new();
    d.steps
        .iter()
        .map(|(note, input)| {
            let mut figure = None;
            let output = match *input {
                "deg" => {
                    calc.settings.angle = Angle::Deg;
                    Ok("angle mode: degrees".to_string())
                }
                "rad" => {
                    calc.settings.angle = Angle::Rad;
                    Ok("angle mode: radians".to_string())
                }
                src => calc.eval_line(src).map(|mut o| {
                    figure = o.figure.take();
                    match o.approx {
                        Some(a) => format!("{}  ≈ {}", o.text, a),
                        None => o.text,
                    }
                }),
            };
            StepResult { note, input, output, figure }
        })
        .collect()
}
