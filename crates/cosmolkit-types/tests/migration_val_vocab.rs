use std::{collections::HashSet, str::FromStr};

use cosmolkit_types::{
    BondDirection, BondOrder, BondStereo, ChiralTag, ELEMENTS, ELEMENTS_WITH_DUMMY, Element,
    Hybridization,
};

const SYMBOLS: [&str; 119] = [
    "*", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

const BOND_ORDERS: [(BondOrder, i64, &str); 22] = [
    (BondOrder::Unspecified, 0, "UNSPECIFIED"),
    (BondOrder::Single, 1, "SINGLE"),
    (BondOrder::Double, 2, "DOUBLE"),
    (BondOrder::Triple, 3, "TRIPLE"),
    (BondOrder::Quadruple, 4, "QUADRUPLE"),
    (BondOrder::Quintuple, 5, "QUINTUPLE"),
    (BondOrder::Hextuple, 6, "HEXTUPLE"),
    (BondOrder::OneAndHalf, 7, "ONEANDAHALF"),
    (BondOrder::TwoAndHalf, 8, "TWOANDAHALF"),
    (BondOrder::ThreeAndHalf, 9, "THREEANDAHALF"),
    (BondOrder::FourAndHalf, 10, "FOURANDAHALF"),
    (BondOrder::FiveAndHalf, 11, "FIVEANDAHALF"),
    (BondOrder::Aromatic, 12, "AROMATIC"),
    (BondOrder::Ionic, 13, "IONIC"),
    (BondOrder::Hydrogen, 14, "HYDROGEN"),
    (BondOrder::ThreeCenter, 15, "THREECENTER"),
    (BondOrder::DativeOne, 16, "DATIVEONE"),
    (BondOrder::Dative, 17, "DATIVE"),
    (BondOrder::DativeLeft, 18, "DATIVEL"),
    (BondOrder::DativeRight, 19, "DATIVER"),
    (BondOrder::Other, 20, "OTHER"),
    (BondOrder::Zero, 21, "ZERO"),
];

const CHIRAL_TAGS: [(ChiralTag, i64, &str); 9] = [
    (ChiralTag::Unspecified, 0, "CHI_UNSPECIFIED"),
    (ChiralTag::TetrahedralCw, 1, "CHI_TETRAHEDRAL_CW"),
    (ChiralTag::TetrahedralCcw, 2, "CHI_TETRAHEDRAL_CCW"),
    (ChiralTag::Other, 3, "CHI_OTHER"),
    (ChiralTag::Tetrahedral, 4, "CHI_TETRAHEDRAL"),
    (ChiralTag::Allene, 5, "CHI_ALLENE"),
    (ChiralTag::SquarePlanar, 6, "CHI_SQUAREPLANAR"),
    (ChiralTag::TrigonalBipyramidal, 7, "CHI_TRIGONALBIPYRAMIDAL"),
    (ChiralTag::Octahedral, 8, "CHI_OCTAHEDRAL"),
];

const HYBRIDIZATIONS: [(Hybridization, i64, &str); 9] = [
    (Hybridization::Unspecified, 0, "UNSPECIFIED"),
    (Hybridization::S, 1, "S"),
    (Hybridization::Sp, 2, "SP"),
    (Hybridization::Sp2, 3, "SP2"),
    (Hybridization::Sp3, 4, "SP3"),
    (Hybridization::Sp2d, 5, "SP2D"),
    (Hybridization::Sp3d, 6, "SP3D"),
    (Hybridization::Sp3d2, 7, "SP3D2"),
    (Hybridization::Other, 8, "OTHER"),
];

const BOND_DIRECTIONS: [(BondDirection, i64, &str); 7] = [
    (BondDirection::None, 0, "NONE"),
    (BondDirection::BeginWedge, 1, "BEGINWEDGE"),
    (BondDirection::BeginDash, 2, "BEGINDASH"),
    (BondDirection::EndDownRight, 3, "ENDDOWNRIGHT"),
    (BondDirection::EndUpRight, 4, "ENDUPRIGHT"),
    (BondDirection::EitherDouble, 5, "EITHERDOUBLE"),
    (BondDirection::Unknown, 6, "UNKNOWN"),
];

const BOND_STEREOS: [(BondStereo, i64, &str); 8] = [
    (BondStereo::None, 0, "STEREONONE"),
    (BondStereo::Any, 1, "STEREOANY"),
    (BondStereo::Z, 2, "STEREOZ"),
    (BondStereo::E, 3, "STEREOE"),
    (BondStereo::Cis, 4, "STEREOCIS"),
    (BondStereo::Trans, 5, "STEREOTRANS"),
    (BondStereo::AtropCw, 6, "STEREOATROPCW"),
    (BondStereo::AtropCcw, 7, "STEREOATROPCCW"),
];

macro_rules! assert_enum_conversions {
    ($type:ty, $cases:expr) => {
        for &(value, code, name) in $cases {
            assert_eq!(value.rdkit_code(), code);
            assert_eq!(value.rdkit_name(), name);
            assert_eq!(<$type>::from_rdkit_code(code), Some(value));
            assert_eq!(<$type>::from_rdkit_name(name), Some(value));
            assert_eq!(value.to_string(), name);
        }
    };
}

#[test]
fn element_number_and_symbol_domains_match_the_source_table() {
    for (atomic_number, expected_symbol) in SYMBOLS.iter().enumerate() {
        let element = Element::from_atomic_number(atomic_number as u8).unwrap();
        assert_eq!(element.atomic_number(), atomic_number as u8);
        assert_eq!(element.symbol(), *expected_symbol);
        assert_eq!(Element::from_symbol(expected_symbol), Some(element));
    }
    assert_eq!(Element::from_atomic_number(119), None);
    assert_eq!(Element::from_atomic_number(u8::MAX), None);
}

#[test]
fn element_symbol_parsing_is_case_sensitive_and_supports_source_aliases() {
    assert_eq!(Element::from_symbol("Uut"), Some(Element::NH));
    assert_eq!(Element::from_symbol("Uup"), Some(Element::MC));
    assert_eq!(Element::from_symbol("Nh").unwrap().symbol(), "Nh");
    assert_eq!(Element::from_symbol("Mc").unwrap().symbol(), "Mc");
    for invalid in ["", "Xx", "uut", "uup", "nh", "mc", "CL"] {
        assert_eq!(Element::from_symbol(invalid), None);
        let error = Element::from_str(invalid).unwrap_err();
        assert_eq!(error.input(), invalid);
    }
}

#[test]
fn element_iteration_has_exact_length_and_atomic_number_order() {
    assert_eq!(Element::iter().len(), 118);
    assert_eq!(Element::iter_with_dummy().len(), 119);
    assert_eq!(Element::iter().collect::<Vec<_>>(), ELEMENTS);
    assert_eq!(
        Element::iter_with_dummy().collect::<Vec<_>>(),
        ELEMENTS_WITH_DUMMY
    );
    assert_eq!(Element::iter().next(), Some(Element::H));
    assert_eq!(Element::iter().next_back(), Some(Element::OG));
    assert_eq!(Element::iter_with_dummy().next(), Some(Element::DUMMY));
}

#[test]
fn element_display_parse_and_serde_use_canonical_symbols() {
    for &element in &ELEMENTS_WITH_DUMMY {
        assert_eq!(element.to_string(), element.symbol());
        assert_eq!(Element::from_str(element.symbol()), Ok(element));
        let encoded = serde_json::to_string(&element).unwrap();
        assert_eq!(encoded, format!("\"{}\"", element.symbol()));
        assert_eq!(serde_json::from_str::<Element>(&encoded).unwrap(), element);
    }
    assert_eq!(
        serde_json::from_str::<Element>("\"Uut\"").unwrap(),
        Element::NH
    );
    assert!(serde_json::from_str::<Element>("\"Xx\"").is_err());
}

#[test]
fn every_enum_member_round_trips_through_source_code_and_name() {
    assert_enum_conversions!(BondOrder, &BOND_ORDERS);
    assert_enum_conversions!(ChiralTag, &CHIRAL_TAGS);
    assert_enum_conversions!(Hybridization, &HYBRIDIZATIONS);
    assert_enum_conversions!(BondDirection, &BOND_DIRECTIONS);
    assert_enum_conversions!(BondStereo, &BOND_STEREOS);
}

#[test]
fn enum_conversions_reject_unknown_codes_and_names() {
    for code in [-1, 22, i64::MAX] {
        assert_eq!(BondOrder::from_rdkit_code(code), None);
    }
    for code in [-1, 9, i64::MAX] {
        assert_eq!(ChiralTag::from_rdkit_code(code), None);
        assert_eq!(Hybridization::from_rdkit_code(code), None);
    }
    for code in [-1, 7, i64::MAX] {
        assert_eq!(BondDirection::from_rdkit_code(code), None);
    }
    for code in [-1, 8, i64::MAX] {
        assert_eq!(BondStereo::from_rdkit_code(code), None);
    }
    assert_eq!(BondOrder::from_rdkit_name("single"), None);
    assert_eq!(ChiralTag::from_rdkit_name("TETRAHEDRAL_CW"), None);
    assert_eq!(Hybridization::from_rdkit_name("sp3"), None);
    assert_eq!(BondDirection::from_rdkit_name("ENDUPRIGHT_"), None);
    assert_eq!(BondStereo::from_rdkit_name("NONE"), None);
}

#[test]
fn corrected_source_ordering_regressions_are_fixed() {
    assert_eq!(BondOrder::Hydrogen.rdkit_code(), 14);
    assert_eq!(BondOrder::ThreeCenter.rdkit_code(), 15);
    assert_eq!(BondOrder::DativeOne.rdkit_code(), 16);
    assert_eq!(BondOrder::Dative.rdkit_code(), 17);
    assert_eq!(BondOrder::DativeLeft.rdkit_code(), 18);
    assert_eq!(BondOrder::DativeRight.rdkit_code(), 19);
    assert_eq!(BondDirection::EndDownRight.rdkit_code(), 3);
    assert_eq!(BondDirection::EndUpRight.rdkit_code(), 4);
}

#[test]
fn explicit_member_lists_have_exact_cardinality_and_unique_values() {
    assert_eq!(BOND_ORDERS.len(), 22);
    assert_eq!(CHIRAL_TAGS.len(), 9);
    assert_eq!(HYBRIDIZATIONS.len(), 9);
    assert_eq!(BOND_DIRECTIONS.len(), 7);
    assert_eq!(BOND_STEREOS.len(), 8);

    for codes in [
        BOND_ORDERS.iter().map(|entry| entry.1).collect::<Vec<_>>(),
        CHIRAL_TAGS.iter().map(|entry| entry.1).collect(),
        HYBRIDIZATIONS.iter().map(|entry| entry.1).collect(),
        BOND_DIRECTIONS.iter().map(|entry| entry.1).collect(),
        BOND_STEREOS.iter().map(|entry| entry.1).collect(),
    ] {
        assert_eq!(
            codes.iter().copied().collect::<HashSet<_>>().len(),
            codes.len()
        );
    }
}
