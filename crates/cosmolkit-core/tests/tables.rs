use cosmolkit_core::{
    PeriodicTableError, atomic_mass, element_info, isotope_abundance, isotope_mass,
    most_common_isotope, most_common_isotope_mass,
};
use cosmolkit_types::Element;

#[test]
fn element_info_covers_dummy_and_every_canonical_source_row_in_order() {
    for element in Element::iter_with_dummy() {
        let info = element_info(element);
        assert_eq!(info.element, element);
        assert_eq!(info.atomic_number, element.atomic_number());
        assert_eq!(info.symbol, element.symbol());
        assert!(info.period <= 7);
        assert!(!info.valences.is_empty());
        assert!(info.rb0 >= 0.0);
        assert!(info.atomic_weight >= 0.0);
    }

    let dummy = element_info(Element::DUMMY);
    assert_eq!(dummy.valences, &[-1]);
    assert_eq!(dummy.atomic_weight, 0.0);

    let carbon = element_info(Element::C);
    assert_eq!(carbon.period, 2);
    assert_eq!(carbon.outer_electrons, 4);
    assert_eq!(carbon.valences, &[4]);
    assert_eq!(carbon.rb0, 0.77);
    assert_eq!(carbon.atomic_weight, 12.011);

    let sulfur = element_info(Element::S);
    assert_eq!(sulfur.valences, &[2, 4, 6]);
    let xenon = element_info(Element::XE);
    assert_eq!(xenon.valences, &[0, 2, 4, 6]);
    let oganesson = element_info(Element::OG);
    assert_eq!(oganesson.period, 7);
    assert_eq!(oganesson.atomic_weight, 294.0);
}

#[test]
fn legacy_symbols_select_the_canonical_numeric_row() {
    let nihonium = Element::from_symbol("Uut").expect("legacy Nh alias");
    let moscovium = Element::from_symbol("Uup").expect("legacy Mc alias");
    assert_eq!(nihonium, Element::NH);
    assert_eq!(moscovium, Element::MC);
    assert_eq!(element_info(nihonium).symbol, "Nh");
    assert_eq!(element_info(moscovium).symbol, "Mc");
    assert!(Element::from_symbol("Xx").is_none());
    assert!(Element::from_atomic_number(119).is_none());
}

#[test]
fn atomic_mass_distinguishes_average_and_exact_isotope_values() {
    assert_eq!(atomic_mass(Element::C, None), Ok(12.011));
    assert_eq!(atomic_mass(Element::C, Some(12)), Ok(12.0));
    assert_eq!(atomic_mass(Element::C, Some(13)), Ok(13.00335484));
    assert_eq!(atomic_mass(Element::O, Some(18)), Ok(17.999161));
    assert_eq!(atomic_mass(Element::U, Some(238)), Ok(238.0507882));
}

#[test]
fn absent_isotopes_return_exact_structured_errors_without_mass_number_fallback() {
    for (element, isotope) in [
        (Element::DUMMY, 0),
        (Element::C, 1),
        (Element::C, 99),
        (Element::OG, u16::MAX),
    ] {
        let expected = PeriodicTableError::UnknownIsotope { element, isotope };
        assert_eq!(isotope_mass(element, isotope), Err(expected));
        assert_eq!(isotope_abundance(element, isotope), Err(expected));
        assert_eq!(atomic_mass(element, Some(isotope)), Err(expected));
    }
}

#[test]
fn present_zero_abundance_isotopes_are_not_confused_with_missing_rows() {
    assert_eq!(isotope_mass(Element::H, 3), Ok(3.016049278));
    assert_eq!(isotope_abundance(Element::H, 3), Ok(0.0));
    assert_eq!(isotope_mass(Element::OG, 294), Ok(294.21392));
    assert_eq!(isotope_abundance(Element::OG, 294), Ok(0.0));
}

#[test]
fn most_common_isotope_columns_preserve_source_equal_different_and_missing_rows() {
    const DIFFERENT_MASS_ROWS: [(Element, u16, f64, f64); 4] = [
        (Element::RN, 222, 222.0175706, 222.0175777),
        (Element::RA, 226, 226.0254026, 226.0254098),
        (Element::RG, 281, 281.16537, 281.16636),
        (Element::CN, 285, 285.17411, 285.17712),
    ];
    const MISSING_ISOTOPE_ROWS: [(Element, u16, f64); 5] = [
        (Element::DB, 268, 268.12545),
        (Element::SG, 271, 271.13347),
        (Element::BH, 270, 270.13362),
        (Element::HS, 269, 269.13406),
        (Element::MT, 278, 278.15481),
    ];

    let mut equal_rows = 0;
    let mut different_rows = 0;
    let mut missing_rows = 0;
    for element in Element::iter() {
        let isotope = most_common_isotope(element);
        assert_ne!(isotope, 0, "{} must have a common isotope", element);
        let periodic_mass = most_common_isotope_mass(element);
        match isotope_mass(element, isotope) {
            Ok(isotope_mass) if isotope_mass == periodic_mass => equal_rows += 1,
            Ok(isotope_mass) => {
                let expected = DIFFERENT_MASS_ROWS
                    .iter()
                    .find(|(expected_element, ..)| *expected_element == element)
                    .expect("every unequal source pair must be inventoried");
                assert_eq!(isotope, expected.1, "{} mass number", element);
                assert_eq!(periodic_mass, expected.2, "{} periodic-row mass", element);
                assert_eq!(isotope_mass, expected.3, "{} isotope-row mass", element);
                different_rows += 1;
            }
            Err(error) => {
                let expected = MISSING_ISOTOPE_ROWS
                    .iter()
                    .find(|(expected_element, ..)| *expected_element == element)
                    .expect("every missing source isotope row must be inventoried");
                assert_eq!(isotope, expected.1, "{} mass number", element);
                assert_eq!(periodic_mass, expected.2, "{} periodic-row mass", element);
                assert_eq!(
                    error,
                    PeriodicTableError::UnknownIsotope { element, isotope },
                    "{} missing isotope-row error",
                    element
                );
                missing_rows += 1;
            }
        }
    }
    assert_eq!(equal_rows, 109);
    assert_eq!(different_rows, DIFFERENT_MASS_ROWS.len());
    assert_eq!(missing_rows, MISSING_ISOTOPE_ROWS.len());
}

#[test]
fn pinned_radium_and_radon_regressions_match_rdkit_source_tests() {
    assert_eq!(most_common_isotope(Element::RA), 226);
    assert!((isotope_mass(Element::RA, 226).unwrap() - 226.02540).abs() < 1.0e-4);
    assert_eq!(most_common_isotope(Element::RN), 222);
    assert!((isotope_mass(Element::RN, 222).unwrap() - 222.01757).abs() < 1.0e-4);
}

#[test]
fn immutable_table_lookups_are_deterministic_across_threads() {
    let handles = (0..8)
        .map(|_| {
            std::thread::spawn(|| {
                (
                    element_info(Element::BR),
                    isotope_mass(Element::BR, 79),
                    isotope_abundance(Element::BR, 79),
                )
            })
        })
        .collect::<Vec<_>>();

    for handle in handles {
        let (info, mass, abundance) = handle.join().expect("lookup thread");
        assert_eq!(info.atomic_number, 35);
        assert_eq!(info.symbol, "Br");
        assert_eq!(mass, Ok(78.9183371));
        assert_eq!(abundance, Ok(50.69));
    }
}
