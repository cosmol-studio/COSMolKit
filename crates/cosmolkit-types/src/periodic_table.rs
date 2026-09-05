//! Source-aligned periodic-table data shared by COSMolKit value types.

#[derive(Debug, Clone, Copy)]
struct Entry {
    symbol: &'static str,
    row: u8,
    outer_electrons: i32,
    valences: &'static [i32],
}

const TABLE: [Entry; 119] = [
    Entry {
        symbol: "*",
        row: 0,
        outer_electrons: 0,
        valences: &[-1],
    },
    Entry {
        symbol: "H",
        row: 1,
        outer_electrons: 1,
        valences: &[1],
    },
    Entry {
        symbol: "He",
        row: 1,
        outer_electrons: 2,
        valences: &[0],
    },
    Entry {
        symbol: "Li",
        row: 2,
        outer_electrons: 1,
        valences: &[1, -1],
    },
    Entry {
        symbol: "Be",
        row: 2,
        outer_electrons: 2,
        valences: &[2],
    },
    Entry {
        symbol: "B",
        row: 2,
        outer_electrons: 3,
        valences: &[3],
    },
    Entry {
        symbol: "C",
        row: 2,
        outer_electrons: 4,
        valences: &[4],
    },
    Entry {
        symbol: "N",
        row: 2,
        outer_electrons: 5,
        valences: &[3],
    },
    Entry {
        symbol: "O",
        row: 2,
        outer_electrons: 6,
        valences: &[2],
    },
    Entry {
        symbol: "F",
        row: 2,
        outer_electrons: 7,
        valences: &[1],
    },
    Entry {
        symbol: "Ne",
        row: 2,
        outer_electrons: 8,
        valences: &[0],
    },
    Entry {
        symbol: "Na",
        row: 3,
        outer_electrons: 1,
        valences: &[1, -1],
    },
    Entry {
        symbol: "Mg",
        row: 3,
        outer_electrons: 2,
        valences: &[2, -1],
    },
    Entry {
        symbol: "Al",
        row: 3,
        outer_electrons: 3,
        valences: &[3],
    },
    Entry {
        symbol: "Si",
        row: 3,
        outer_electrons: 4,
        valences: &[4],
    },
    Entry {
        symbol: "P",
        row: 3,
        outer_electrons: 5,
        valences: &[3, 5],
    },
    Entry {
        symbol: "S",
        row: 3,
        outer_electrons: 6,
        valences: &[2, 4, 6],
    },
    Entry {
        symbol: "Cl",
        row: 3,
        outer_electrons: 7,
        valences: &[1],
    },
    Entry {
        symbol: "Ar",
        row: 3,
        outer_electrons: 8,
        valences: &[0],
    },
    Entry {
        symbol: "K",
        row: 4,
        outer_electrons: 1,
        valences: &[1, -1],
    },
    Entry {
        symbol: "Ca",
        row: 4,
        outer_electrons: 2,
        valences: &[2, -1],
    },
    Entry {
        symbol: "Sc",
        row: 4,
        outer_electrons: 3,
        valences: &[-1],
    },
    Entry {
        symbol: "Ti",
        row: 4,
        outer_electrons: 4,
        valences: &[-1],
    },
    Entry {
        symbol: "V",
        row: 4,
        outer_electrons: 5,
        valences: &[-1],
    },
    Entry {
        symbol: "Cr",
        row: 4,
        outer_electrons: 6,
        valences: &[-1],
    },
    Entry {
        symbol: "Mn",
        row: 4,
        outer_electrons: 7,
        valences: &[-1],
    },
    Entry {
        symbol: "Fe",
        row: 4,
        outer_electrons: 8,
        valences: &[-1],
    },
    Entry {
        symbol: "Co",
        row: 4,
        outer_electrons: 9,
        valences: &[-1],
    },
    Entry {
        symbol: "Ni",
        row: 4,
        outer_electrons: 10,
        valences: &[-1],
    },
    Entry {
        symbol: "Cu",
        row: 4,
        outer_electrons: 11,
        valences: &[-1],
    },
    Entry {
        symbol: "Zn",
        row: 4,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Ga",
        row: 4,
        outer_electrons: 3,
        valences: &[3],
    },
    Entry {
        symbol: "Ge",
        row: 4,
        outer_electrons: 4,
        valences: &[4],
    },
    Entry {
        symbol: "As",
        row: 4,
        outer_electrons: 5,
        valences: &[3, 5],
    },
    Entry {
        symbol: "Se",
        row: 4,
        outer_electrons: 6,
        valences: &[2, 4, 6],
    },
    Entry {
        symbol: "Br",
        row: 4,
        outer_electrons: 7,
        valences: &[1],
    },
    Entry {
        symbol: "Kr",
        row: 4,
        outer_electrons: 8,
        valences: &[0],
    },
    Entry {
        symbol: "Rb",
        row: 5,
        outer_electrons: 1,
        valences: &[1, -1],
    },
    Entry {
        symbol: "Sr",
        row: 5,
        outer_electrons: 2,
        valences: &[2, -1],
    },
    Entry {
        symbol: "Y",
        row: 5,
        outer_electrons: 3,
        valences: &[-1],
    },
    Entry {
        symbol: "Zr",
        row: 5,
        outer_electrons: 4,
        valences: &[-1],
    },
    Entry {
        symbol: "Nb",
        row: 5,
        outer_electrons: 5,
        valences: &[-1],
    },
    Entry {
        symbol: "Mo",
        row: 5,
        outer_electrons: 6,
        valences: &[-1],
    },
    Entry {
        symbol: "Tc",
        row: 5,
        outer_electrons: 7,
        valences: &[-1],
    },
    Entry {
        symbol: "Ru",
        row: 5,
        outer_electrons: 8,
        valences: &[-1],
    },
    Entry {
        symbol: "Rh",
        row: 5,
        outer_electrons: 9,
        valences: &[-1],
    },
    Entry {
        symbol: "Pd",
        row: 5,
        outer_electrons: 10,
        valences: &[-1],
    },
    Entry {
        symbol: "Ag",
        row: 5,
        outer_electrons: 11,
        valences: &[-1],
    },
    Entry {
        symbol: "Cd",
        row: 5,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "In",
        row: 5,
        outer_electrons: 3,
        valences: &[3],
    },
    Entry {
        symbol: "Sn",
        row: 5,
        outer_electrons: 4,
        valences: &[2, 4],
    },
    Entry {
        symbol: "Sb",
        row: 5,
        outer_electrons: 5,
        valences: &[3, 5],
    },
    Entry {
        symbol: "Te",
        row: 5,
        outer_electrons: 6,
        valences: &[2, 4, 6],
    },
    Entry {
        symbol: "I",
        row: 5,
        outer_electrons: 7,
        valences: &[1, 3, 5],
    },
    Entry {
        symbol: "Xe",
        row: 5,
        outer_electrons: 8,
        valences: &[0, 2, 4, 6],
    },
    Entry {
        symbol: "Cs",
        row: 6,
        outer_electrons: 1,
        valences: &[1],
    },
    Entry {
        symbol: "Ba",
        row: 6,
        outer_electrons: 2,
        valences: &[2, -1],
    },
    Entry {
        symbol: "La",
        row: 6,
        outer_electrons: 3,
        valences: &[-1],
    },
    Entry {
        symbol: "Ce",
        row: 6,
        outer_electrons: 4,
        valences: &[-1],
    },
    Entry {
        symbol: "Pr",
        row: 6,
        outer_electrons: 3,
        valences: &[-1],
    },
    Entry {
        symbol: "Nd",
        row: 6,
        outer_electrons: 4,
        valences: &[-1],
    },
    Entry {
        symbol: "Pm",
        row: 6,
        outer_electrons: 5,
        valences: &[-1],
    },
    Entry {
        symbol: "Sm",
        row: 6,
        outer_electrons: 6,
        valences: &[-1],
    },
    Entry {
        symbol: "Eu",
        row: 6,
        outer_electrons: 7,
        valences: &[-1],
    },
    Entry {
        symbol: "Gd",
        row: 6,
        outer_electrons: 8,
        valences: &[-1],
    },
    Entry {
        symbol: "Tb",
        row: 6,
        outer_electrons: 9,
        valences: &[-1],
    },
    Entry {
        symbol: "Dy",
        row: 6,
        outer_electrons: 10,
        valences: &[-1],
    },
    Entry {
        symbol: "Ho",
        row: 6,
        outer_electrons: 11,
        valences: &[-1],
    },
    Entry {
        symbol: "Er",
        row: 6,
        outer_electrons: 12,
        valences: &[-1],
    },
    Entry {
        symbol: "Tm",
        row: 6,
        outer_electrons: 13,
        valences: &[-1],
    },
    Entry {
        symbol: "Yb",
        row: 6,
        outer_electrons: 14,
        valences: &[-1],
    },
    Entry {
        symbol: "Lu",
        row: 6,
        outer_electrons: 15,
        valences: &[-1],
    },
    Entry {
        symbol: "Hf",
        row: 6,
        outer_electrons: 4,
        valences: &[-1],
    },
    Entry {
        symbol: "Ta",
        row: 6,
        outer_electrons: 5,
        valences: &[-1],
    },
    Entry {
        symbol: "W",
        row: 6,
        outer_electrons: 6,
        valences: &[-1],
    },
    Entry {
        symbol: "Re",
        row: 6,
        outer_electrons: 7,
        valences: &[-1],
    },
    Entry {
        symbol: "Os",
        row: 6,
        outer_electrons: 8,
        valences: &[-1],
    },
    Entry {
        symbol: "Ir",
        row: 6,
        outer_electrons: 9,
        valences: &[-1],
    },
    Entry {
        symbol: "Pt",
        row: 6,
        outer_electrons: 10,
        valences: &[-1],
    },
    Entry {
        symbol: "Au",
        row: 6,
        outer_electrons: 11,
        valences: &[-1],
    },
    Entry {
        symbol: "Hg",
        row: 6,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Tl",
        row: 6,
        outer_electrons: 3,
        valences: &[-1],
    },
    Entry {
        symbol: "Pb",
        row: 6,
        outer_electrons: 4,
        valences: &[2, 4],
    },
    Entry {
        symbol: "Bi",
        row: 6,
        outer_electrons: 5,
        valences: &[3, 5],
    },
    Entry {
        symbol: "Po",
        row: 6,
        outer_electrons: 6,
        valences: &[2, 4, 6],
    },
    Entry {
        symbol: "At",
        row: 6,
        outer_electrons: 7,
        valences: &[1, 3, 5],
    },
    Entry {
        symbol: "Rn",
        row: 6,
        outer_electrons: 8,
        valences: &[0],
    },
    Entry {
        symbol: "Fr",
        row: 7,
        outer_electrons: 1,
        valences: &[1],
    },
    Entry {
        symbol: "Ra",
        row: 7,
        outer_electrons: 2,
        valences: &[2, -1],
    },
    Entry {
        symbol: "Ac",
        row: 7,
        outer_electrons: 3,
        valences: &[-1],
    },
    Entry {
        symbol: "Th",
        row: 7,
        outer_electrons: 4,
        valences: &[-1],
    },
    Entry {
        symbol: "Pa",
        row: 7,
        outer_electrons: 3,
        valences: &[-1],
    },
    Entry {
        symbol: "U",
        row: 7,
        outer_electrons: 4,
        valences: &[-1],
    },
    Entry {
        symbol: "Np",
        row: 7,
        outer_electrons: 5,
        valences: &[-1],
    },
    Entry {
        symbol: "Pu",
        row: 7,
        outer_electrons: 6,
        valences: &[-1],
    },
    Entry {
        symbol: "Am",
        row: 7,
        outer_electrons: 7,
        valences: &[-1],
    },
    Entry {
        symbol: "Cm",
        row: 7,
        outer_electrons: 8,
        valences: &[-1],
    },
    Entry {
        symbol: "Bk",
        row: 7,
        outer_electrons: 9,
        valences: &[-1],
    },
    Entry {
        symbol: "Cf",
        row: 7,
        outer_electrons: 10,
        valences: &[-1],
    },
    Entry {
        symbol: "Es",
        row: 7,
        outer_electrons: 11,
        valences: &[-1],
    },
    Entry {
        symbol: "Fm",
        row: 7,
        outer_electrons: 12,
        valences: &[-1],
    },
    Entry {
        symbol: "Md",
        row: 7,
        outer_electrons: 13,
        valences: &[-1],
    },
    Entry {
        symbol: "No",
        row: 7,
        outer_electrons: 14,
        valences: &[-1],
    },
    Entry {
        symbol: "Lr",
        row: 7,
        outer_electrons: 15,
        valences: &[-1],
    },
    Entry {
        symbol: "Rf",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Db",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Sg",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Bh",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Hs",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Mt",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Ds",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Rg",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Cn",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Nh",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Fl",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Mc",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Lv",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Ts",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
    Entry {
        symbol: "Og",
        row: 7,
        outer_electrons: 2,
        valences: &[-1],
    },
];

// BEGIN RDKIT CPP DATA atomicData::Mass / PeriodicTable::getAtomicWeight
// RDKit✔️✔️: atomicData stores the source-provided average atomic mass for each
// RDKit✔️✔️: atomic number and returns it through getAtomicWeight().
const ATOMIC_WEIGHTS: [f64; 119] = [
    0.0, 1.008, 4.003, 6.941, 9.012, 10.812, 12.011, 14.007, 15.999, 18.998, 20.18, 22.99, 24.305, 26.982, 28.086,
    30.974, 32.067, 35.453, 39.948, 39.098, 40.078, 44.956, 47.867, 50.944, 51.996, 54.938, 55.845, 58.933, 58.693,
    63.546, 65.39, 69.723, 72.61, 74.922, 78.96, 79.904, 83.8, 85.468, 87.62, 88.906, 91.224, 92.906, 95.94, 98.0,
    101.07, 102.906, 106.42, 107.868, 112.412, 114.818, 118.711, 121.76, 127.6, 126.904, 131.29, 132.905, 137.328,
    138.906, 140.116, 140.908, 144.24, 145.0, 150.36, 151.964, 157.25, 158.925, 162.5, 164.93, 167.26, 168.934, 173.04,
    174.967, 178.49, 180.948, 183.84, 186.207, 190.23, 192.217, 195.078, 196.967, 200.59, 204.383, 207.2, 208.98,
    209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.038, 231.036, 238.029, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0,
    252.0, 257.0, 258.0, 259.0, 262.0, 267.0, 268.0, 269.0, 270.0, 269.0, 278.0, 281.0, 281.0, 285.0, 284.179, 289.190,
    288.193, 293.204, 292.207, 294.214,
];
// END RDKIT CPP DATA atomicData::Mass / PeriodicTable::getAtomicWeight

// BEGIN RDKIT CPP FUNCTION PeriodicTable::getRb0 / atomicData::Rb0
// RDKit✔️✔️: double Rb0() const { return rB0; }
// RDKit✔️✔️: double getRb0(UINT atomicNumber) const {
// RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
// RDKit✔️✔️:   return byanum[atomicNumber].Rb0();
// RDKit✔️✔️: }
// RDKit✔️✔️: //  rB0
// RDKit✔️✔️: //  ...
// RDKit✔️✔️: istr >> rB0;
const RB0: [f64; 119] = [
    0.0, 0.33, 0.7, 1.23, 0.9, 0.82, 0.77, 0.7, 0.66, 0.611, 0.7, 1.54, 1.36, 1.18, 0.937, 0.89, 1.04, 0.997, 1.74,
    2.03, 1.74, 1.44, 1.32, 1.22, 1.18, 1.17, 1.17, 1.16, 1.15, 1.17, 1.25, 1.26, 1.188, 1.2, 1.17, 1.167, 1.91, 2.16,
    1.91, 1.62, 1.45, 1.34, 1.3, 1.27, 1.25, 1.25, 1.28, 1.34, 1.48, 1.44, 1.385, 1.4, 1.378, 1.387, 1.98, 2.35, 1.98,
    1.69, 1.83, 1.82, 1.81, 1.8, 1.8, 1.99, 1.79, 1.76, 1.75, 1.74, 1.73, 1.72, 1.94, 1.72, 1.44, 1.34, 1.3, 1.28,
    1.26, 1.27, 1.3, 1.34, 1.49, 1.48, 1.48, 1.45, 1.46, 1.45, 2.4, 2.0, 1.9, 1.88, 1.79, 1.61, 1.58, 1.55, 1.53, 1.07,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
];
// END RDKIT CPP FUNCTION PeriodicTable::getRb0 / atomicData::Rb0

#[must_use]
pub fn atomic_number_from_symbol(symbol: &str) -> Option<u8> {
    match symbol {
        "*" => Some(0),
        "H" => Some(1),
        "He" => Some(2),
        "Li" => Some(3),
        "Be" => Some(4),
        "B" => Some(5),
        "C" => Some(6),
        "N" => Some(7),
        "O" => Some(8),
        "F" => Some(9),
        "Ne" => Some(10),
        "Na" => Some(11),
        "Mg" => Some(12),
        "Al" => Some(13),
        "Si" => Some(14),
        "P" => Some(15),
        "S" => Some(16),
        "Cl" => Some(17),
        "Ar" => Some(18),
        "K" => Some(19),
        "Ca" => Some(20),
        "Sc" => Some(21),
        "Ti" => Some(22),
        "V" => Some(23),
        "Cr" => Some(24),
        "Mn" => Some(25),
        "Fe" => Some(26),
        "Co" => Some(27),
        "Ni" => Some(28),
        "Cu" => Some(29),
        "Zn" => Some(30),
        "Ga" => Some(31),
        "Ge" => Some(32),
        "As" => Some(33),
        "Se" => Some(34),
        "Br" => Some(35),
        "Kr" => Some(36),
        "Rb" => Some(37),
        "Sr" => Some(38),
        "Y" => Some(39),
        "Zr" => Some(40),
        "Nb" => Some(41),
        "Mo" => Some(42),
        "Tc" => Some(43),
        "Ru" => Some(44),
        "Rh" => Some(45),
        "Pd" => Some(46),
        "Ag" => Some(47),
        "Cd" => Some(48),
        "In" => Some(49),
        "Sn" => Some(50),
        "Sb" => Some(51),
        "Te" => Some(52),
        "I" => Some(53),
        "Xe" => Some(54),
        "Cs" => Some(55),
        "Ba" => Some(56),
        "La" => Some(57),
        "Ce" => Some(58),
        "Pr" => Some(59),
        "Nd" => Some(60),
        "Pm" => Some(61),
        "Sm" => Some(62),
        "Eu" => Some(63),
        "Gd" => Some(64),
        "Tb" => Some(65),
        "Dy" => Some(66),
        "Ho" => Some(67),
        "Er" => Some(68),
        "Tm" => Some(69),
        "Yb" => Some(70),
        "Lu" => Some(71),
        "Hf" => Some(72),
        "Ta" => Some(73),
        "W" => Some(74),
        "Re" => Some(75),
        "Os" => Some(76),
        "Ir" => Some(77),
        "Pt" => Some(78),
        "Au" => Some(79),
        "Hg" => Some(80),
        "Tl" => Some(81),
        "Pb" => Some(82),
        "Bi" => Some(83),
        "Po" => Some(84),
        "At" => Some(85),
        "Rn" => Some(86),
        "Fr" => Some(87),
        "Ra" => Some(88),
        "Ac" => Some(89),
        "Th" => Some(90),
        "Pa" => Some(91),
        "U" => Some(92),
        "Np" => Some(93),
        "Pu" => Some(94),
        "Am" => Some(95),
        "Cm" => Some(96),
        "Bk" => Some(97),
        "Cf" => Some(98),
        "Es" => Some(99),
        "Fm" => Some(100),
        "Md" => Some(101),
        "No" => Some(102),
        "Lr" => Some(103),
        "Rf" => Some(104),
        "Db" => Some(105),
        "Sg" => Some(106),
        "Bh" => Some(107),
        "Hs" => Some(108),
        "Mt" => Some(109),
        "Ds" => Some(110),
        "Rg" => Some(111),
        "Cn" => Some(112),
        "Nh" | "Uut" => Some(113),
        "Fl" => Some(114),
        "Mc" | "Uup" => Some(115),
        "Lv" => Some(116),
        "Ts" => Some(117),
        "Og" => Some(118),
        _ => None,
    }
}

#[must_use]
pub fn symbol(atomic_number: u8) -> Option<&'static str> {
    TABLE.get(usize::from(atomic_number)).map(|entry| entry.symbol)
}

#[must_use]
pub fn period(atomic_number: u8) -> Option<u8> {
    TABLE.get(usize::from(atomic_number)).map(|entry| entry.row)
}

#[must_use]
pub fn outer_electrons(atomic_number: u8) -> Option<i32> {
    TABLE.get(usize::from(atomic_number)).map(|entry| entry.outer_electrons)
}

#[must_use]
pub fn valences(atomic_number: u8) -> Option<&'static [i32]> {
    TABLE.get(usize::from(atomic_number)).map(|entry| entry.valences)
}

#[must_use]
pub fn rb0(atomic_number: u8) -> f64 {
    RB0.get(usize::from(atomic_number)).copied().unwrap_or(0.0)
}

#[must_use]
pub fn atomic_weight(atomic_number: u8) -> f64 {
    ATOMIC_WEIGHTS.get(usize::from(atomic_number)).copied().unwrap_or(0.0)
}
