macro_rules! element_symbols {
    ($($num:literal => $symbol:literal),+ $(,)?) => {
        pub(crate) const ELEMENT_SYMBOLS: [&str; 119] = [$($symbol),+];

        pub(crate) fn atomic_number(symbol: &str) -> Option<u8> {
            match symbol {
                $($symbol => Some($num),)+
                _ => None,
            }
        }
    };
}

element_symbols! {
    0 => "*", 1 => "H", 2 => "He", 3 => "Li", 4 => "Be", 5 => "B", 6 => "C", 7 => "N",
    8 => "O", 9 => "F", 10 => "Ne", 11 => "Na", 12 => "Mg", 13 => "Al", 14 => "Si",
    15 => "P", 16 => "S", 17 => "Cl", 18 => "Ar", 19 => "K", 20 => "Ca", 21 => "Sc",
    22 => "Ti", 23 => "V", 24 => "Cr", 25 => "Mn", 26 => "Fe", 27 => "Co", 28 => "Ni",
    29 => "Cu", 30 => "Zn", 31 => "Ga", 32 => "Ge", 33 => "As", 34 => "Se", 35 => "Br",
    36 => "Kr", 37 => "Rb", 38 => "Sr", 39 => "Y", 40 => "Zr", 41 => "Nb", 42 => "Mo",
    43 => "Tc", 44 => "Ru", 45 => "Rh", 46 => "Pd", 47 => "Ag", 48 => "Cd", 49 => "In",
    50 => "Sn", 51 => "Sb", 52 => "Te", 53 => "I", 54 => "Xe", 55 => "Cs", 56 => "Ba",
    57 => "La", 58 => "Ce", 59 => "Pr", 60 => "Nd", 61 => "Pm", 62 => "Sm", 63 => "Eu",
    64 => "Gd", 65 => "Tb", 66 => "Dy", 67 => "Ho", 68 => "Er", 69 => "Tm", 70 => "Yb",
    71 => "Lu", 72 => "Hf", 73 => "Ta", 74 => "W", 75 => "Re", 76 => "Os", 77 => "Ir",
    78 => "Pt", 79 => "Au", 80 => "Hg", 81 => "Tl", 82 => "Pb", 83 => "Bi", 84 => "Po",
    85 => "At", 86 => "Rn", 87 => "Fr", 88 => "Ra", 89 => "Ac", 90 => "Th", 91 => "Pa",
    92 => "U", 93 => "Np", 94 => "Pu", 95 => "Am", 96 => "Cm", 97 => "Bk", 98 => "Cf",
    99 => "Es", 100 => "Fm", 101 => "Md", 102 => "No", 103 => "Lr", 104 => "Rf",
    105 => "Db", 106 => "Sg", 107 => "Bh", 108 => "Hs", 109 => "Mt", 110 => "Ds",
    111 => "Rg", 112 => "Cn", 113 => "Nh", 114 => "Fl", 115 => "Mc", 116 => "Lv",
    117 => "Ts", 118 => "Og",
}

pub(crate) fn element_symbol(atomic_num: u8) -> Option<&'static str> {
    ELEMENT_SYMBOLS.get(usize::from(atomic_num)).copied()
}

pub(crate) const OUTER_ELECTRONS: [i32; 119] = [
    0, 1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2,
    3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3, 4, 5,
    6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4,
    3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
];

pub(crate) fn n_outer_electrons(atomic_num: u8) -> Option<i32> {
    OUTER_ELECTRONS.get(usize::from(atomic_num)).copied()
}

pub(crate) const MOST_COMMON_ISOTOPES: [i32; 119] = [
    0, 1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 40, 39, 40, 45, 48, 51, 52,
    55, 56, 59, 58, 63, 64, 69, 74, 75, 80, 79, 84, 85, 88, 89, 90, 93, 98, 97, 102, 103, 106, 107,
    114, 115, 120, 121, 130, 127, 132, 133, 138, 139, 140, 141, 142, 145, 152, 153, 158, 159, 164,
    165, 166, 169, 174, 175, 180, 181, 184, 187, 192, 193, 195, 197, 202, 205, 208, 209, 209, 210,
    222, 223, 226, 227, 232, 231, 238, 236, 238, 241, 243, 247, 249, 252, 257, 258, 259, 262, 267,
    268, 271, 270, 269, 278, 281, 281, 285, 284, 289, 288, 293, 292, 294,
];

pub(crate) fn most_common_isotope(atomic_num: u8) -> Option<i32> {
    MOST_COMMON_ISOTOPES.get(usize::from(atomic_num)).copied()
}

pub(crate) const RDKIT_VALENCE_LISTS: [&[i32]; 119] = [
    &[-1],
    &[1],
    &[0],
    &[1],
    &[2],
    &[3],
    &[4],
    &[3],
    &[2],
    &[1],
    &[0],
    &[1],
    &[2, -1],
    &[3, 6],
    &[4, 6],
    &[3, 5, 7],
    &[2, 4, 6],
    &[1],
    &[0],
    &[1],
    &[2, -1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[3],
    &[4],
    &[3, 5, 7],
    &[2, 4, 6],
    &[1],
    &[0],
    &[1],
    &[2, -1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[3],
    &[2, 4],
    &[3, 5, 7],
    &[2, 4, 6],
    &[1, 3, 5],
    &[0, 2, 4, 6],
    &[1],
    &[2, -1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[3],
    &[2, 4],
    &[3, 5, 7],
    &[2, 4, 6],
    &[1, 3, 5],
    &[0],
    &[1],
    &[2, -1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
    &[-1],
];

pub(crate) fn valence_list(atomic_num: u8) -> Option<&'static [i32]> {
    RDKIT_VALENCE_LISTS.get(usize::from(atomic_num)).copied()
}

pub(crate) fn exact_isotope_mass(atomic_num: u8, isotope: u16) -> Option<f64> {
    match (atomic_num, isotope) {
        (1, 2) => Some(2.014_101_778),
        (1, 3) => Some(3.016_049_278),
        (6, 13) => Some(13.003_354_835),
        (6, 14) => Some(14.003_241_989),
        (7, 15) => Some(15.000_108_899),
        (8, 17) => Some(16.999_131_757),
        (8, 18) => Some(17.999_159_613),
        _ => None,
    }
}

pub(crate) fn average_atomic_weight(atomic_num: u8) -> Option<f64> {
    match atomic_num {
        1 => Some(1.008),
        2 => Some(4.003),
        3 => Some(6.941),
        4 => Some(9.012),
        5 => Some(10.812),
        6 => Some(12.011),
        7 => Some(14.007),
        8 => Some(15.999),
        9 => Some(18.998),
        10 => Some(20.18),
        11 => Some(22.99),
        12 => Some(24.305),
        13 => Some(26.982),
        14 => Some(28.086),
        15 => Some(30.974),
        16 => Some(32.067),
        17 => Some(35.453),
        18 => Some(39.948),
        19 => Some(39.098),
        20 => Some(40.078),
        35 => Some(79.904),
        53 => Some(126.904),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn element_symbol_table_covers_rdkit_periodic_table() {
        assert_eq!(super::ELEMENT_SYMBOLS.len(), 119);
        assert_eq!(super::element_symbol(0), Some("*"));
        assert_eq!(super::element_symbol(19), Some("K"));
        assert_eq!(super::element_symbol(26), Some("Fe"));
        assert_eq!(super::element_symbol(29), Some("Cu"));
        assert_eq!(super::element_symbol(118), Some("Og"));
        assert_eq!(super::atomic_number("Og"), Some(118));
        assert_eq!(super::atomic_number("Xx"), None);
    }

    #[test]
    fn rdkit_periodic_tables_cover_expected_range() {
        assert_eq!(super::OUTER_ELECTRONS.len(), 119);
        assert_eq!(super::MOST_COMMON_ISOTOPES.len(), 119);
        assert_eq!(super::RDKIT_VALENCE_LISTS.len(), 119);
        assert_eq!(super::n_outer_electrons(118), Some(2));
        assert_eq!(super::most_common_isotope(6), Some(12));
        assert_eq!(super::most_common_isotope(88), Some(226));
        assert_eq!(super::valence_list(15), Some(&[3, 5, 7][..]));
        assert_eq!(super::valence_list(200), None);
    }
}
