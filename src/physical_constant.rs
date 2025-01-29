//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*! # Physical Constants

This chapter describes macros for the values of physical constants,
such as the speed of light, c, and gravitational constant, G. The
values are available in different unit systems, including the standard
MKSA system (meters, kilograms, seconds, amperes) and the CGSM system
(centimeters, grams, seconds, gauss), which is commonly used in
Astronomy.

## References and Further Reading

The authoritative sources for physical constants are the 2006 CODATA
recommended values, published in the article below. Further
information on the values of physical constants is also available from
the NIST website.

P.J. Mohr, B.N. Taylor, D.B. Newell, “CODATA Recommended Values of the
Fundamental Physical Constants: 2006”, Reviews of Modern Physics,
80(2), pp. 633–730 (2008).

<http://www.physics.nist.gov/cuu/Constants/index.html>
<http://physics.nist.gov/Pubs/SP811/appenB9.html>
*/

pub mod mksa {
    // Fundamental Constants
    /// The speed of light in vacuum, c. m / s
    pub const SPEED_OF_LIGHT: f64 = 2.99792458e8;
    /// The permeability of free space, \mu_0. This constant is defined in
    /// the MKSA system only.  `kg m / A^2 s^2`
    pub const VACUUM_PERMEABILITY: f64 = 1.25663706144e-6;
    /// The permittivity of free space, \epsilon_0. This constant is
    /// defined in the MKSA system only.  `A^2 s^4 / kg m^3`
    pub const VACUUM_PERMITTIVITY: f64 = 8.854187817e-12;
    /// Planck’s constant, h. kg m^2 / s
    pub const PLANCKS_CONSTANT_H: f64 = 6.62606896e-34;
    /// Planck’s constant divided by 2\pi, \hbar. kg m^2 / s
    pub const PLANCKS_CONSTANT_HBAR: f64 = 1.05457162825e-34;
    /// Avogadro’s number, N_a. 1 / mol
    pub static NUM_AVOGADRO: f64 = 6.02214199e23;
    /// The molar charge of 1 Faraday. A s / mol
    pub const FARADAY: f64 = 9.64853429775e4;
    /// The Boltzmann constant, k. kg m^2 / K s^2
    pub const BOLTZMANN: f64 = 1.3806504e-23;
    /// The molar gas constant, R_0. kg m^2 / K mol s^2
    pub const MOLAR_GAS: f64 = 8.314472e0;
    /// The standard gas volume, V_0. m^3 / mol
    pub const STANDARD_GAS_VOLUME: f64 = 2.2710981e-2;
    /// The Stefan-Boltzmann radiation constant, \sigma. kg / K^4 s^3
    pub const STEFAN_BOLTZMANN_CONSTANT: f64 = 5.67040047374e-8;
    /// The magnetic field of 1 Gauss. kg / A s^2
    pub const GAUSS: f64 = 1e-4;

    // Astronomy and Astrophysics
    /// The length of 1 astronomical unit (mean earth-sun distance), au. m
    pub const ASTRONOMICAL_UNIT: f64 = 1.49597870691e11;
    /// The gravitational constant, G. m^3 / kg s^2
    pub const GRAVITATIONAL_CONSTANT: f64 = 6.673e-11;
    /// The distance of 1 light-year, ly. m
    pub const LIGHT_YEAR: f64 = 9.46053620707e15;
    /// The distance of 1 parsec, pc. m
    pub const PARSEC: f64 = 3.08567758135e16;
    /// The standard gravitational acceleration on Earth, g. m / s^2
    pub const GRAV_ACCEL: f64 = 9.80665e0;
    /// The mass of the Sun. kg
    pub const SOLAR_MASS: f64 = 1.98892e30;

    // Atomic and Nuclear Physics
    /// The charge of the electron, e. A s
    pub const ELECTRON_CHARGE: f64 = 1.602176487e-19;
    /// The energy of 1 electron volt, eV. kg m^2 / s^2
    pub const ELECTRON_VOLT: f64 = 1.602176487e-19;
    /// The unified atomic mass, amu. kg
    pub const UNIFIED_ATOMIC_MASS: f64 = 1.660538782e-27;
    /// The mass of the electron, m_e. kg
    pub const MASS_ELECTRON: f64 = 9.10938188e-31;
    /// The mass of the muon, m_\mu. kg
    pub const MASS_MUON: f64 = 1.88353109e-28;
    /// The mass of the proton, m_p. kg
    pub const MASS_PROTON: f64 = 1.67262158e-27;
    /// The mass of the neutron, m_n. kg
    pub const MASS_NEUTRON: f64 = 1.67492716e-27;
    /// The electromagnetic fine structure constant \alpha. 1
    pub static NUM_FINE_STRUCTURE: f64 = 7.297352533e-3;
    /// The Rydberg constant, Ry, in units of energy. This is related
    /// to the Rydberg inverse wavelength `R_\infty by Ry = h c
    /// R_\infty. kg m^2 / s^2`
    pub const RYDBERG: f64 = 2.17987196968e-18;
    /// The Bohr radius, a_0. m
    pub const BOHR_RADIUS: f64 = 5.291772083e-11;
    /// The length of 1 angstrom. m
    pub const ANGSTROM: f64 = 1e-10;
    /// The area of 1 barn. m^2
    pub const BARN: f64 = 1e-28;
    /// The Bohr Magneton, \mu_B. A m^2
    pub const BOHR_MAGNETON: f64 = 9.27400899e-24;
    /// The Nuclear Magneton, \mu_N. A m^2
    pub const NUCLEAR_MAGNETON: f64 = 5.05078317e-27;
    /// The absolute value of the magnetic moment of the electron,
    /// \mu_e. The physical magnetic moment of the electron is
    /// negative. A m^2
    pub const ELECTRON_MAGNETIC_MOMENT: f64 = 9.28476362e-24;
    /// The magnetic moment of the proton, \mu_p. A m^2
    pub const PROTON_MAGNETIC_MOMENT: f64 = 1.410606633e-26;
    /// The Thomson cross section, \sigma_T. m^2
    pub const THOMSON_CROSS_SECTION: f64 = 6.65245893699e-29;
    /// The electric dipole moment of 1 Debye, D. A s^2 / m^2
    pub const DEBYE: f64 = 3.33564095198e-30;

    // Measurement of Time
    /// The number of seconds in 1 minute. s
    pub const MINUTE: f64 = 6e1f64;
    /// The number of seconds in 1 hour. s
    pub const HOUR: f64 = 3.6e3f64;
    /// The number of seconds in 1 day. s
    pub const DAY: f64 = 8.64e4f64;
    /// The number of seconds in 1 week. s
    pub const WEEK: f64 = 6.048e5f64;

    // Imperial Units
    /// The length of 1 inch. m
    pub const INCH: f64 = 2.54e-2;
    /// The length of 1 foot. m
    pub const FOOT: f64 = 3.048e-1;
    /// The length of 1 yard. m
    pub const YARD: f64 = 9.144e-1;
    /// The length of 1 mile. m
    pub const MILE: f64 = 1.609344e3;
    /// The length of 1 mil (1/1000th of an inch). m
    pub const MIL: f64 = 2.54e-5;

    // Speed and Nautical Units
    /// The speed of 1 kilometer per hour. m / s
    pub const KILOMETERS_PER_HOUR: f64 = 2.77777777778e-1;
    /// The speed of 1 mile per hour. m / s
    pub const MILES_PER_HOUR: f64 = 4.4704e-1;
    /// The length of 1 nautical mile. m
    pub const NAUTICAL_MILE: f64 = 1.852e3;
    /// The length of 1 fathom. m
    pub const FATHOM: f64 = 1.8288e0;
    /// The speed of 1 knot. m / s
    pub const KNOT: f64 = 5.14444444444e-1;

    // Printers Units
    /// The length of 1 printer’s point (1/72 inch). m
    pub const POINT: f64 = 3.52777777778e-4;
    /// The length of 1 TeX point (1/72.27 inch). m
    pub const TEXPOINT: f64 = 3.51459803515e-4;

    // Volume, Area and Length
    /// The length of 1 micron. m
    pub const MICRON: f64 = 1e-6;
    /// The area of 1 hectare. m^2
    pub const HECTARE: f64 = 1e4;
    /// The area of 1 acre. m^2
    pub const ACRE: f64 = 4.04685642241e3;
    /// The volume of 1 liter. m^3
    pub const LITER: f64 = 1e-3;
    /// The volume of 1 US gallon. m^3
    pub const US_GALLON: f64 = 3.78541178402e-3;
    /// The volume of 1 Canadian gallon. m^3
    pub const CANADIAN_GALLON: f64 = 4.54609e-3;
    /// The volume of 1 UK gallon. m^3
    pub const UK_GALLON: f64 = 4.546092e-3;
    /// The volume of 1 quart. m^3
    pub const QUART: f64 = 9.46352946004e-4;
    /// The volume of 1 pint. m^3
    pub const PINT: f64 = 4.73176473002e-4;
    /// m^3
    pub const CUP: f64 = 2.36588236501e-4;

    // Mass and Weight
    /// The mass of 1 pound. kg
    pub const POUND_MASS: f64 = 4.5359237e-1;
    /// The mass of 1 ounce. kg
    pub const OUNCE_MASS: f64 = 2.8349523125e-2;
    /// The mass of 1 ton. kg
    pub const TON: f64 = 9.0718474e2;
    /// The mass of 1 metric ton (1000 kg). kg
    pub const METRIC_TON: f64 = 1e3;
    /// The mass of 1 UK ton. kg
    pub const UK_TON: f64 = 1.0160469088e3;
    /// The mass of 1 troy ounce. kg
    pub const TROY_OUNCE: f64 = 3.1103475e-2;
    /// The mass of 1 carat. kg
    pub const CARAT: f64 = 2e-4;
    /// The force of 1 gram weight. kg m / s^2
    pub const GRAM_FORCE: f64 = 9.80665e-3;
    /// The force of 1 pound weight. kg m / s^2
    pub const POUND_FORCE: f64 = 4.44822161526e0;
    /// The force of 1 kilopound weight. kg m / s^2
    pub const KILOPOUND_FORCE: f64 = 4.44822161526e3;
    /// The force of 1 poundal. kg m / s^2
    pub const POUNDAL: f64 = 1.38255e-1;

    // Thermal Energy and Power
    /// The energy of 1 calorie. kg m^2 / s^2
    pub const CALORIE: f64 = 4.1868e0;
    /// The energy of 1 British Thermal Unit, btu. kg m^2 / s^2
    pub const BTU: f64 = 1.05505585262e3;
    /// The energy of 1 Therm. kg m^2 / s^2
    pub const THERM: f64 = 1.05506e8;
    /// The power of 1 horsepower. kg m^2 / s^3
    pub const HORSEPOWER: f64 = 7.457e2;

    // Pressure
    /// The pressure of 1 bar. kg / m s^2
    pub const BAR: f64 = 1e5;
    /// The pressure of 1 standard atmosphere. kg / m s^2
    pub const STD_ATMOSPHERE: f64 = 1.01325e5;
    /// The pressure of 1 torr. kg / m s^2
    pub const TORR: f64 = 1.33322368421e2;
    /// The pressure of 1 meter of mercury. kg / m s^2
    pub const METER_OF_MERCURY: f64 = 1.33322368421e5;
    /// The pressure of 1 inch of mercury. kg / m s^2
    pub const INCH_OF_MERCURY: f64 = 3.38638815789e3;
    /// The pressure of 1 inch of water. kg / m s^2
    pub const INCH_OF_WATER: f64 = 2.490889e2;
    /// The pressure of 1 pound per square inch. kg / m s^2
    pub const PSI: f64 = 6.89475729317e3;

    // Viscosity
    /// The dynamic viscosity of 1 poise. kg m^-1 s^-1
    pub const POISE: f64 = 1e-1;
    /// The kinematic viscosity of 1 stokes. m^2 / s
    pub const STOKES: f64 = 1e-4;

    // Light and Illumination
    /// The luminance of 1 stilb. cd / m^2
    pub const STILB: f64 = 1e4;
    /// The luminous flux of 1 lumen. cd sr
    pub const LUMEN: f64 = 1e0;
    /// The illuminance of 1 lux. cd sr / m^2
    pub const LUX: f64 = 1e0;
    /// The illuminance of 1 phot. cd sr / m^2
    pub const PHOT: f64 = 1e4;
    /// The illuminance of 1 footcandle. cd sr / m^2
    pub const FOOTCANDLE: f64 = 1.076e1;
    /// The luminance of 1 lambert. cd sr / m^2
    pub const LAMBERT: f64 = 1e4;
    /// The luminance of 1 footlambert. cd sr / m^2
    pub const FOOTLAMBERT: f64 = 1.07639104e1;

    // Radioactivity
    /// The activity of 1 curie. 1 / s
    pub const CURIE: f64 = 3.7e10;
    /// The exposure of 1 roentgen. A s / kg
    pub const ROENTGEN: f64 = 2.58e-4;
    /// The absorbed dose of 1 rad. m^2 / s^2
    pub const RAD: f64 = 1e-2;

    // Force and Energy
    /// The SI unit of force, 1 Newton. kg m / s^2
    pub const NEWTON: f64 = 1e0;
    /// The force of 1 Dyne = 10^-5 Newton. kg m / s^2
    pub const DYNE: f64 = 1e-5;
    /// The SI unit of energy, 1 Joule. kg m^2 / s^2
    pub const JOULE: f64 = 1e0;
    /// The energy 1 erg = 10^-7 Joule. kg m^2 / s^2
    pub const ERG: f64 = 1e-7;
}

pub mod num {
    // Prefixes : These constants are dimensionless scaling factors.
    /// 10^24
    pub const YOTTA: f64 = 1e24;
    /// 10^21
    pub const ZETTA: f64 = 1e21;
    /// 10^18
    pub const EXA: f64 = 1e18;
    /// 10^15
    pub const PETA: f64 = 1e15;
    /// 10^12
    pub const TERA: f64 = 1e12;
    /// 10^9
    pub const GIGA: f64 = 1e9;
    /// 10^6
    pub const MEGA: f64 = 1e6;
    /// 10^3
    pub const KILO: f64 = 1e3;
    /// 10^-3
    pub const MILLI: f64 = 1e-3;
    /// 10^-6
    pub const MICRO: f64 = 1e-6;
    /// 10^-9
    pub const NANO: f64 = 1e-9;
    /// 10^-12
    pub const PICO: f64 = 1e-12;
    /// 10^-15
    pub const FEMTO: f64 = 1e-15;
    /// 10^-18
    pub const ATTO: f64 = 1e-18;
    /// 10^-21
    pub const ZEPTO: f64 = 1e-21;
    /// 10^-24
    pub const YOCTO: f64 = 1e-24;
}
