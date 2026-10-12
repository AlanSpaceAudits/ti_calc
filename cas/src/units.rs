//! Unit conversion and physical constants (SI 2019 values).

/// (name, dimension, factor to the SI base of that dimension)
const UNITS: &[(&str, &str, f64)] = &[
    // length (m)
    ("m", "length", 1.0),
    ("km", "length", 1e3),
    ("cm", "length", 1e-2),
    ("mm", "length", 1e-3),
    ("um", "length", 1e-6),
    ("nm", "length", 1e-9),
    ("in", "length", 0.0254),
    ("ft", "length", 0.3048),
    ("yd", "length", 0.9144),
    ("mi", "length", 1609.344),
    ("nmi", "length", 1852.0),
    ("au", "length", 149_597_870_700.0),
    ("ly", "length", 9.460_730_472_580_8e15),
    ("pc", "length", 3.085_677_581_491_367e16),
    ("angstrom", "length", 1e-10),
    // mass (kg)
    ("kg", "mass", 1.0),
    ("g", "mass", 1e-3),
    ("mg", "mass", 1e-6),
    ("tonne", "mass", 1e3),
    ("lb", "mass", 0.453_592_37),
    ("oz", "mass", 0.028_349_523_125),
    ("slug", "mass", 14.593_902_937_206_36),
    ("ton", "mass", 907.184_74),
    ("amu", "mass", 1.660_539_066_60e-27),
    // time (s)
    ("s", "time", 1.0),
    ("ms", "time", 1e-3),
    ("us", "time", 1e-6),
    ("ns", "time", 1e-9),
    ("min", "time", 60.0),
    ("hr", "time", 3600.0),
    ("day", "time", 86400.0),
    ("week", "time", 604_800.0),
    ("yr", "time", 31_557_600.0),
    // area (m^2)
    ("m2", "area", 1.0),
    ("km2", "area", 1e6),
    ("cm2", "area", 1e-4),
    ("ft2", "area", 0.092_903_04),
    ("in2", "area", 6.4516e-4),
    ("mi2", "area", 2_589_988.110_336),
    ("acre", "area", 4046.856_422_4),
    ("ha", "area", 1e4),
    // volume (m^3)
    ("m3", "volume", 1.0),
    ("l", "volume", 1e-3),
    ("ml", "volume", 1e-6),
    ("cm3", "volume", 1e-6),
    ("ft3", "volume", 0.028_316_846_592),
    ("in3", "volume", 1.638_706_4e-5),
    ("gal", "volume", 3.785_411_784e-3),
    ("qt", "volume", 9.463_529_46e-4),
    ("pt", "volume", 4.731_764_73e-4),
    ("cup", "volume", 2.365_882_365e-4),
    ("floz", "volume", 2.957_352_956_25e-5),
    ("tbsp", "volume", 1.478_676_478_125e-5),
    ("tsp", "volume", 4.928_921_593_75e-6),
    // speed (m/s)
    ("mps", "speed", 1.0),
    ("kph", "speed", 1.0 / 3.6),
    ("mph", "speed", 0.447_04),
    ("fps", "speed", 0.3048),
    ("knot", "speed", 1852.0 / 3600.0),
    // force (N)
    ("N", "force", 1.0),
    ("kN", "force", 1e3),
    ("lbf", "force", 4.448_221_615_260_5),
    ("dyn", "force", 1e-5),
    ("kgf", "force", 9.806_65),
    // energy (J)
    ("J", "energy", 1.0),
    ("kJ", "energy", 1e3),
    ("cal", "energy", 4.184),
    ("kcal", "energy", 4184.0),
    ("eV", "energy", 1.602_176_634e-19),
    ("kWh", "energy", 3.6e6),
    ("Btu", "energy", 1055.055_852_62),
    ("ftlb", "energy", 1.355_817_948_331_400_4),
    ("erg", "energy", 1e-7),
    // power (W)
    ("W", "power", 1.0),
    ("kW", "power", 1e3),
    ("MW", "power", 1e6),
    ("hp", "power", 745.699_871_582_270_2),
    // pressure (Pa)
    ("Pa", "pressure", 1.0),
    ("kPa", "pressure", 1e3),
    ("MPa", "pressure", 1e6),
    ("bar", "pressure", 1e5),
    ("atm", "pressure", 101_325.0),
    ("psi", "pressure", 6894.757_293_168_361),
    ("torr", "pressure", 101_325.0 / 760.0),
    ("mmHg", "pressure", 133.322_387_415),
    // angle (rad)
    ("rad", "angle", 1.0),
    ("deg", "angle", std::f64::consts::PI / 180.0),
    ("grad", "angle", std::f64::consts::PI / 200.0),
    ("arcmin", "angle", std::f64::consts::PI / 10_800.0),
    ("arcsec", "angle", std::f64::consts::PI / 648_000.0),
    ("rev", "angle", 2.0 * std::f64::consts::PI),
    // data (bytes)
    ("B", "data", 1.0),
    ("bit", "data", 0.125),
    ("kB", "data", 1e3),
    ("MB", "data", 1e6),
    ("GB", "data", 1e9),
    ("TB", "data", 1e12),
    ("KiB", "data", 1024.0),
    ("MiB", "data", 1_048_576.0),
    ("GiB", "data", 1_073_741_824.0),
];

fn temp_to_kelvin(v: f64, u: &str) -> Option<f64> {
    Some(match u {
        "K" => v,
        "degC" | "C" => v + 273.15,
        "degF" | "F" => (v - 32.0) * 5.0 / 9.0 + 273.15,
        "degR" | "R" => v * 5.0 / 9.0,
        _ => return None,
    })
}

fn kelvin_to(v: f64, u: &str) -> Option<f64> {
    Some(match u {
        "K" => v,
        "degC" | "C" => v - 273.15,
        "degF" | "F" => (v - 273.15) * 9.0 / 5.0 + 32.0,
        "degR" | "R" => v * 9.0 / 5.0,
        _ => return None,
    })
}

pub fn convert(v: f64, from: &str, to: &str) -> Result<f64, String> {
    if let (Some(k), true) = (temp_to_kelvin(v, from), kelvin_to(0.0, to).is_some()) {
        return Ok(kelvin_to(k, to).unwrap());
    }
    let find = |u: &str| UNITS.iter().find(|(n, _, _)| *n == u).ok_or_else(|| format!("unknown unit '{u}'"));
    let (a, b) = (find(from)?, find(to)?);
    if a.1 != b.1 {
        return Err(format!("cannot convert {} ({}) to {} ({})", from, a.1, to, b.1));
    }
    Ok(v * a.2 / b.2)
}

pub fn unit_list() -> String {
    let mut out = String::new();
    let mut last = "";
    for (n, d, _) in UNITS {
        if *d != last {
            out.push_str(&format!("\n{d}: "));
            last = d;
        }
        out.push_str(n);
        out.push(' ');
    }
    out.push_str("\ntemperature: K degC degF degR");
    out.trim().to_string()
}

/// (name, value in SI units, description)
pub const CONSTANTS: &[(&str, f64, &str)] = &[
    ("_c", 299_792_458.0, "speed of light, m/s"),
    ("_g", 9.806_65, "standard gravity, m/s^2"),
    ("_G", 6.674_30e-11, "gravitational constant, m^3/(kg s^2)"),
    ("_h", 6.626_070_15e-34, "Planck constant, J s"),
    ("_hbar", 1.054_571_817e-34, "reduced Planck constant, J s"),
    ("_k", 1.380_649e-23, "Boltzmann constant, J/K"),
    ("_Na", 6.022_140_76e23, "Avogadro constant, 1/mol"),
    ("_R", 8.314_462_618, "gas constant, J/(mol K)"),
    ("_q", 1.602_176_634e-19, "elementary charge, C"),
    ("_me", 9.109_383_7015e-31, "electron mass, kg"),
    ("_mp", 1.672_621_923_69e-27, "proton mass, kg"),
    ("_mn", 1.674_927_498_04e-27, "neutron mass, kg"),
    ("_eps0", 8.854_187_8128e-12, "vacuum permittivity, F/m"),
    ("_mu0", 1.256_637_062_12e-6, "vacuum permeability, N/A^2"),
    ("_sigma", 5.670_374_419e-8, "Stefan-Boltzmann constant, W/(m^2 K^4)"),
    ("_F", 96_485.332_12, "Faraday constant, C/mol"),
    ("_Rinf", 10_973_731.568_160, "Rydberg constant, 1/m"),
    ("_a0", 5.291_772_109_03e-11, "Bohr radius, m"),
    ("_atm", 101_325.0, "standard atmosphere, Pa"),
];
