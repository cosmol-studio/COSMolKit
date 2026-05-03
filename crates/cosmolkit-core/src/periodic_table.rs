pub(crate) const ELEMENT_SYMBOLS: [&str; 119] = [
    "*", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

pub(crate) fn element_symbol(atomic_num: u8) -> Option<&'static str> {
    ELEMENT_SYMBOLS.get(atomic_num as usize).copied()
}

pub(crate) fn atomic_number(symbol: &str) -> Option<u8> {
    ELEMENT_SYMBOLS
        .iter()
        .position(|candidate| *candidate == symbol)
        .and_then(|idx| u8::try_from(idx).ok())
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
}
