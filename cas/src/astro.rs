//! Positional astronomy and orbital mechanics. Angles use the current angle mode; `kdeg` is degrees per mode unit.
//! Formulas: Meeus (1998) for dates, sidereal time, obliquity and coordinate rotations; the Astronomical Almanac
//! low-precision Sun; Bennett (1982) refraction; Kasten and Young (1989) air mass.
use crate::expr::*;
use crate::numeric;
use crate::simplify::neg;

pub const NAMES: &[&str] = &[
    "jd", "caldate", "gmst", "lst", "hms", "tohms", "eq2hor", "hor2eq", "eq2ecl", "ecl2eq", "obliquity", "angsep", "refraction", "airmass", "sunpos", "riseset", "kepler",
    "trueanomaly", "orbperiod", "visviva", "vescape", "parallaxdist", "distmod", "absmag",
];

fn num(e: &Expr) -> Result<f64, String> {
    if *e == INF {
        return Ok(f64::INFINITY);
    }
    if *e == neg(INF) {
        return Ok(f64::NEG_INFINITY);
    }
    if e.is_constant() { numeric::real_const(e) } else { None }.ok_or_else(|| "expected a real number".to_string())
}

fn named(pairs: Vec<(&str, f64)>) -> Expr {
    Expr::List(pairs.into_iter().map(|(n, v)| Expr::Rel(Rel::Eq, Box::new(sym(n)), Box::new(Expr::Float(v)))).collect())
}

pub fn julian_date(y: f64, m: f64, d: f64) -> f64 {
    let (mut y, mut m) = (y, m);
    if m <= 2.0 {
        y -= 1.0;
        m += 12.0;
    }
    // Gregorian calendar from 1582-10-15
    let gregorian = (y, m, d) >= (1582.0, 10.0, 15.0);
    let a = (y / 100.0).floor();
    let b = if gregorian { 2.0 - a + (a / 4.0).floor() } else { 0.0 };
    (365.25 * (y + 4716.0)).floor() + (30.6001 * (m + 1.0)).floor() + d + b - 1524.5
}

pub fn calendar_date(jd: f64) -> (i64, i64, f64) {
    let z = (jd + 0.5).floor();
    let f = jd + 0.5 - z;
    let a = if z < 2_299_161.0 {
        z
    } else {
        let alpha = ((z - 1_867_216.25) / 36_524.25).floor();
        z + 1.0 + alpha - (alpha / 4.0).floor()
    };
    let b = a + 1524.0;
    let c = ((b - 122.1) / 365.25).floor();
    let d = (365.25 * c).floor();
    let e = ((b - d) / 30.6001).floor();
    let day = b - d - (30.6001 * e).floor() + f;
    let month = if e < 14.0 { e - 1.0 } else { e - 13.0 };
    let year = if month > 2.0 { c - 4716.0 } else { c - 4715.0 };
    (year as i64, month as i64, day)
}

/// Greenwich mean sidereal time in degrees.
pub fn gmst_deg(jd: f64) -> f64 {
    let t = (jd - 2_451_545.0) / 36_525.0;
    (280.460_618_37 + 360.985_647_366_29 * (jd - 2_451_545.0) + 0.000_387_933 * t * t - t * t * t / 38_710_000.0).rem_euclid(360.0)
}

pub fn obliquity_deg(jd: f64) -> f64 {
    let t = (jd - 2_451_545.0) / 36_525.0;
    23.439_291_111 - 0.013_004_167 * t - 1.639e-7 * t * t + 5.036e-7 * t * t * t
}

/// Altitude and azimuth (from north, eastward) from hour angle, declination, latitude; all degrees.
pub fn eq_to_hor(ha: f64, dec: f64, lat: f64) -> (f64, f64) {
    let (h, d, p) = (ha.to_radians(), dec.to_radians(), lat.to_radians());
    let alt = (p.sin() * d.sin() + p.cos() * d.cos() * h.cos()).asin();
    let az = (-h.sin() * d.cos()).atan2(d.sin() * p.cos() - d.cos() * p.sin() * h.cos());
    (alt.to_degrees(), az.to_degrees().rem_euclid(360.0))
}

/// Eccentric anomaly from mean anomaly (radians) by Newton's method.
pub fn kepler_e(m: f64, e: f64) -> f64 {
    let mut ea = if e < 0.8 { m } else { std::f64::consts::PI };
    for _ in 0..100 {
        let step = (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
        ea -= step;
        if step.abs() < 1e-14 {
            break;
        }
    }
    ea
}

fn sexagesimal(v: f64, marks: [&str; 3]) -> String {
    let total = (v.abs() * 3600.0 * 1e4).round() as u128;
    let (a, b, c) = (total / 36_000_000, total / 600_000 % 60, (total % 600_000) as f64 / 1e4);
    let secs = format!("{c:.4}");
    let secs = secs.trim_end_matches('0').trim_end_matches('.');
    format!("{}{}{}{}{}{}{}", if v < 0.0 { "-" } else { "" }, a, marks[0], b, marks[1], secs, marks[2])
}

pub fn call(name: &str, a: &[Expr], kdeg: f64) -> Option<Result<Expr, String>> {
    if !NAMES.contains(&name) {
        return None;
    }
    let need = |n: usize| if a.len() < n { Err(format!("{name} needs {n} argument(s)")) } else { Ok(()) };
    let x = |i: usize| -> Result<f64, String> { num(a.get(i).ok_or_else(|| format!("{name}: missing argument"))?) };
    let deg = |i: usize| -> Result<f64, String> { Ok(x(i)? * kdeg) };
    let opt = |i: usize, d: f64| -> Result<f64, String> { a.get(i).map_or(Ok(d), num) };
    let ang = |d: f64| Expr::Float(d / kdeg);
    Some((|| -> Result<Expr, String> {
        Ok(match name {
            "jd" => {
                // jd(year, month, day [, hour, minute, second]) in UT
                need(3)?;
                let day = x(2)? + (opt(3, 0.0)? + opt(4, 0.0)? / 60.0 + opt(5, 0.0)? / 3600.0) / 24.0;
                Expr::Float(julian_date(x(0)?, x(1)?, day))
            }
            "caldate" => {
                need(1)?;
                let (y, m, d) = calendar_date(x(0)?);
                let hours = d.fract() * 24.0;
                Expr::Str(format!("{y}-{m:02}-{:02} {} UT", d.floor() as i64, sexagesimal(hours, ["h", "m", "s"])))
            }
            "gmst" => {
                need(1)?;
                ang(gmst_deg(x(0)?))
            }
            "lst" => {
                // lst(jd, east longitude)
                need(2)?;
                ang((gmst_deg(x(0)?) + deg(1)?).rem_euclid(360.0))
            }
            "hms" => {
                // an angle entered as hours, minutes, seconds of time
                need(1)?;
                let h = x(0)?;
                let frac = opt(1, 0.0)? / 60.0 + opt(2, 0.0)? / 3600.0;
                ang(15.0 * if h < 0.0 { h - frac } else { h + frac })
            }
            "tohms" => {
                need(1)?;
                Expr::Str(sexagesimal(deg(0)? / 15.0, ["h", "m", "s"]))
            }
            "eq2hor" => {
                // eq2hor(hour angle, declination, latitude)
                need(3)?;
                let (alt, az) = eq_to_hor(deg(0)?, deg(1)?, deg(2)?);
                named(vec![("alt", alt / kdeg), ("az", az / kdeg)])
            }
            "hor2eq" => {
                // hor2eq(altitude, azimuth, latitude)
                need(3)?;
                let (h, az, p) = (deg(0)?.to_radians(), deg(1)?.to_radians(), deg(2)?.to_radians());
                let dec = (p.sin() * h.sin() + p.cos() * h.cos() * az.cos()).asin();
                let ha = (-az.sin() * h.cos()).atan2((h.sin() - p.sin() * dec.sin()) / p.cos());
                named(vec![("ha", ha.to_degrees().rem_euclid(360.0) / kdeg), ("dec", dec.to_degrees() / kdeg)])
            }
            "eq2ecl" | "ecl2eq" => {
                // (ra, dec [, obliquity]) or (lon, lat [, obliquity]); default obliquity is J2000
                need(2)?;
                let (l, b) = (deg(0)?.to_radians(), deg(1)?.to_radians());
                let eps = match a.get(2) {
                    Some(e) => num(e)? * kdeg,
                    None => obliquity_deg(2_451_545.0),
                }
                .to_radians();
                let s = if name == "eq2ecl" { -1.0 } else { 1.0 };
                let lat = (b.sin() * eps.cos() + s * b.cos() * eps.sin() * l.sin()).asin();
                let lon = (l.sin() * eps.cos() - s * b.tan() * eps.sin()).atan2(l.cos());
                let (n1, n2) = if name == "eq2ecl" { ("lon", "lat") } else { ("ra", "dec") };
                named(vec![(n1, lon.to_degrees().rem_euclid(360.0) / kdeg), (n2, lat.to_degrees() / kdeg)])
            }
            "obliquity" => {
                need(1)?;
                ang(obliquity_deg(x(0)?))
            }
            "angsep" => {
                // angsep(ra1, dec1, ra2, dec2)
                need(4)?;
                let (a1, d1, a2, d2) = (deg(0)?.to_radians(), deg(1)?.to_radians(), deg(2)?.to_radians(), deg(3)?.to_radians());
                let h = ((d2 - d1) / 2.0).sin().powi(2) + d1.cos() * d2.cos() * ((a2 - a1) / 2.0).sin().powi(2);
                ang((2.0 * h.sqrt().asin()).to_degrees())
            }
            "refraction" => {
                // Bennett: apparent altitude -> refraction
                need(1)?;
                let h = deg(0)?;
                ang(1.0 / ((h + 7.31 / (h + 4.4)).to_radians().tan()) / 60.0)
            }
            "airmass" => {
                need(1)?;
                let h = deg(0)?;
                Expr::Float(1.0 / (h.to_radians().sin() + 0.505_72 * (h + 6.079_95).powf(-1.6364)))
            }
            "sunpos" => {
                // low-precision apparent Sun for a Julian date
                need(1)?;
                let n = x(0)? - 2_451_545.0;
                let l = (280.460 + 0.985_647_4 * n).rem_euclid(360.0);
                let g = (357.528 + 0.985_600_3 * n).rem_euclid(360.0).to_radians();
                let lam = (l + 1.915 * g.sin() + 0.020 * (2.0 * g).sin()).to_radians();
                let eps = (23.439 - 0.000_000_4 * n).to_radians();
                let ra = (eps.cos() * lam.sin()).atan2(lam.cos()).to_degrees().rem_euclid(360.0);
                let dec = (eps.sin() * lam.sin()).asin().to_degrees();
                let dist = 1.000_14 - 0.016_71 * g.cos() - 0.000_14 * (2.0 * g).cos();
                Expr::List(vec![
                    Expr::Rel(Rel::Eq, Box::new(sym("ra")), Box::new(ang(ra))),
                    Expr::Rel(Rel::Eq, Box::new(sym("dec")), Box::new(ang(dec))),
                    Expr::Rel(Rel::Eq, Box::new(sym("lon")), Box::new(ang(lam.to_degrees().rem_euclid(360.0)))),
                    Expr::Rel(Rel::Eq, Box::new(sym("dist_au")), Box::new(Expr::Float(dist))),
                ])
            }
            "riseset" => {
                // riseset(declination, latitude [, horizon altitude]): hour angle at rising/setting
                need(2)?;
                let (d, p, h0) = (deg(0)?.to_radians(), deg(1)?.to_radians(), (opt(2, 0.0)? * kdeg).to_radians());
                let c = (h0.sin() - p.sin() * d.sin()) / (p.cos() * d.cos());
                if c > 1.0 {
                    return Err("the body never rises at this latitude".into());
                }
                if c < -1.0 {
                    return Err("the body never sets at this latitude".into());
                }
                let ha = c.acos().to_degrees();
                let az_rise = (d.sin() / p.cos()).clamp(-1.0, 1.0).acos().to_degrees();
                named(vec![("ha", ha / kdeg), ("az_rise", az_rise / kdeg), ("az_set", (360.0 - az_rise) / kdeg), ("hours_up", 2.0 * ha / 15.0)])
            }
            "kepler" => {
                // kepler(mean anomaly, eccentricity) -> eccentric anomaly
                need(2)?;
                let e = x(1)?;
                if !(0.0..1.0).contains(&e) {
                    return Err("eccentricity must be from 0 up to (not including) 1".into());
                }
                ang(kepler_e(deg(0)?.to_radians(), e).to_degrees())
            }
            "trueanomaly" => {
                // trueanomaly(eccentric anomaly, eccentricity)
                need(2)?;
                let (ea, e) = (deg(0)?.to_radians(), x(1)?);
                ang((2.0 * (((1.0 + e) / (1.0 - e)).sqrt() * (ea / 2.0).tan()).atan()).to_degrees().rem_euclid(360.0))
            }
            _ => return Err(String::new()),
        })
    })())
    .filter(|r| !matches!(r, Err(e) if e.is_empty()))
}
