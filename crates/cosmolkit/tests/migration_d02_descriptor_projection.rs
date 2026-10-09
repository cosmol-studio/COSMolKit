//! Fixed facade/projection regressions; no reference program or corpus setup.
//! Proposed conditions only; independent p1/ROOT acceptance remains required.
use cosmolkit::{
    BINDING_CONTRACT, BindingDefault, BindingItem, BindingKind, BindingOwner, BindingReceiver,
    DescriptorError, DescriptorReadError, FunctionStatus, Molecule, SmilesParseParams, StateModel,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() < tolerance,
        "{actual:?} vs source {expected:?}"
    );
}

fn fixture() -> serde_json::Value {
    serde_json::from_str(include_str!(
        "../../../testdata/descriptors/d02_canonical_fixed.json"
    ))
    .unwrap()
}

#[test]
fn d02_fixed_source_scalar_queries() {
    let data = fixture();
    let cases = data["scalar_projection_cases"].as_array().unwrap();
    assert_eq!(cases.len(), 17);
    for case in cases {
        let name = case["canonical"].as_str().unwrap();
        let molecule = Molecule::from_smiles(case["smiles"].as_str().unwrap()).unwrap();
        let before = molecule.clone();
        let actual = match name {
            "hall_kier_alpha" => molecule.hall_kier_alpha(),
            "kappa_1" => molecule.kappa_1(),
            "kappa_2" => molecule.kappa_2(),
            "kappa_3" => molecule.kappa_3(),
            "phi" => molecule.phi(),
            "chi_0_v" => molecule.chi_0_v(),
            "chi_1_v" => molecule.chi_1_v(),
            "chi_2_v" => molecule.chi_2_v(),
            "chi_3_v" => molecule.chi_3_v(),
            "chi_4_v" => molecule.chi_4_v(),
            "chi_n_v" => molecule.chi_n_v(case["order"].as_u64().unwrap().try_into().unwrap()),
            "chi_0_n" => molecule.chi_0_n(),
            "chi_1_n" => molecule.chi_1_n(),
            "chi_2_n" => molecule.chi_2_n(),
            "chi_3_n" => molecule.chi_3_n(),
            "chi_4_n" => molecule.chi_4_n(),
            "chi_n_n" => molecule.chi_n_n(case["order"].as_u64().unwrap().try_into().unwrap()),
            _ => panic!("unrecognized frozen query {name}"),
        }
        .unwrap();
        close(
            actual,
            case["expected"].as_f64().unwrap(),
            case["tolerance_literal"].as_str().unwrap().parse().unwrap(),
        );
        assert_eq!(molecule, before, "query preserves receiver {name}");
    }
}

#[test]
fn d02_hall_kier_contributions_are_owned_atom_rows() {
    let data = fixture();
    for case in data["hall_contributions"].as_array().unwrap() {
        let molecule = Molecule::from_smiles(case["smiles"].as_str().unwrap()).unwrap();
        let (alpha, mut rows) = molecule.hall_kier_alpha_with_contributions().unwrap();
        close(alpha, case["alpha"].as_f64().unwrap(), 1e-12);
        assert_eq!(rows.len(), molecule.num_atoms());
        for (actual, expected) in rows.iter().zip(case["rows"].as_array().unwrap()) {
            close(*actual, expected.as_f64().unwrap(), 1e-12);
        }
        rows[0] = 123.0;
        let (_, fresh) = molecule.hall_kier_alpha_with_contributions().unwrap();
        close(fresh[0], case["rows"][0].as_f64().unwrap(), 1e-12);
    }
}

#[test]
fn d02_mqn_full42_order_force_and_owned_values() {
    let data = fixture();
    for case in data["mqn_cases"].as_array().unwrap() {
        let molecule = Molecule::from_smiles(case["smiles"].as_str().unwrap()).unwrap();
        let expected = case["values"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| u32::try_from(x.as_u64().unwrap()).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(expected.len(), 42);
        let before = molecule.clone();
        let mut result = molecule.mqns(false).unwrap();
        assert_eq!(result, expected);
        assert_eq!(molecule.mqns(true).unwrap(), expected);
        result[0] = u32::MAX;
        assert_eq!(molecule.mqns(false).unwrap(), expected);
        assert_eq!(molecule, before);
    }
}

#[test]
fn d02_generic_chi_zero_is_not_fixed_zero_and_u32_wrap_is_preserved() {
    let molecule = Molecule::from_smiles("C").unwrap();
    assert_eq!(molecule.chi_0_v().unwrap(), 0.0);
    assert_eq!(molecule.chi_0_n().unwrap(), 0.0);
    assert_eq!(molecule.chi_n_v(0).unwrap(), 1.0);
    assert_eq!(molecule.chi_n_n(0).unwrap(), 1.0);
    assert_eq!(molecule.chi_n_v(u32::MAX).unwrap(), 0.0);
    assert_eq!(molecule.chi_n_n(u32::MAX).unwrap(), 0.0);
}

#[test]
fn d02_raw_molecule_missing_valence_has_no_invented_cause() {
    let molecule = Molecule::from_smiles_with_params(
        "CCC",
        &SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    let before = molecule.clone();
    for result in [
        molecule.mqns(false).map(|_| ()),
        molecule.chi_0_v().map(|_| ()),
        molecule.chi_1_v().map(|_| ()),
        molecule.chi_2_v().map(|_| ()),
        molecule.chi_3_v().map(|_| ()),
        molecule.chi_4_v().map(|_| ()),
        molecule.chi_n_v(2).map(|_| ()),
        molecule.chi_0_n().map(|_| ()),
        molecule.chi_1_n().map(|_| ()),
        molecule.chi_2_n().map(|_| ()),
        molecule.chi_3_n().map(|_| ()),
        molecule.chi_4_n().map(|_| ()),
        molecule.chi_n_n(2).map(|_| ()),
    ] {
        let error = result.unwrap_err();
        assert!(matches!(error, DescriptorReadError::MissingPreparedValence));
        assert!(std::error::Error::source(&error).is_none());
    }
    assert_eq!(molecule.hall_kier_alpha().unwrap(), 0.0);
    assert_eq!(
        molecule.hall_kier_alpha_with_contributions().unwrap(),
        (0.0, vec![0.0; 3])
    );
    assert_eq!(molecule.kappa_1().unwrap(), 3.0);
    assert_eq!(molecule.kappa_2().unwrap(), 2.0);
    assert_eq!(molecule.kappa_3().unwrap(), 0.0);
    assert_eq!(molecule.phi().unwrap(), 2.0);
    assert_eq!(molecule, before);
}

#[test]
fn d02_owned_domain_error_is_structural() {
    // Transport construction, not a claimed live-query failure.
    let error = DescriptorReadError::Algorithm {
        source: DescriptorError::InvalidHallKierContributionRows {
            actual: 1,
            minimum: 2,
        },
    };
    assert!(matches!(
        std::error::Error::source(&error)
            .unwrap()
            .downcast_ref::<DescriptorError>(),
        Some(DescriptorError::InvalidHallKierContributionRows {
            actual: 1,
            minimum: 2
        })
    ));
}

#[test]
fn d02_binding_contracts_cover_exact_query_signatures_and_defaults() {
    for (name, javascript, output, parameter) in [
        ("hall_kier_alpha", "hallKierAlpha", "f64", ""),
        (
            "hall_kier_alpha_with_contributions",
            "hallKierAlphaWithContributions",
            "(f64,Vec<f64>)",
            "",
        ),
        ("kappa_1", "kappa1", "f64", ""),
        ("kappa_2", "kappa2", "f64", ""),
        ("kappa_3", "kappa3", "f64", ""),
        ("phi", "phi", "f64", ""),
        ("mqns", "mqns", "Vec<u32>", "force"),
        ("chi_0_v", "chi0V", "f64", ""),
        ("chi_1_v", "chi1V", "f64", ""),
        ("chi_2_v", "chi2V", "f64", ""),
        ("chi_3_v", "chi3V", "f64", ""),
        ("chi_4_v", "chi4V", "f64", ""),
        ("chi_n_v", "chiNV", "f64", "order"),
        ("chi_0_n", "chi0N", "f64", ""),
        ("chi_1_n", "chi1N", "f64", ""),
        ("chi_2_n", "chi2N", "f64", ""),
        ("chi_3_n", "chi3N", "f64", ""),
        ("chi_4_n", "chi4N", "f64", ""),
        ("chi_n_n", "chiNN", "f64", "order"),
    ] {
        let id = format!("Molecule.{name}");
        let entries = BINDING_CONTRACT
            .iter()
            .filter(|x| x.semantic_id == id)
            .collect::<Vec<_>>();
        assert_eq!(entries.len(), 1);
        let entry = entries[0];
        assert_eq!(entry.item, BindingItem::Callable);
        assert_eq!(entry.owner, BindingOwner::Molecule);
        assert_eq!(
            entry.rust_path.replace(' ', ""),
            format!("crate::Molecule::{name}")
        );
        assert_eq!(entry.python_name, name);
        assert_eq!(entry.javascript_name, javascript);
        assert_eq!(entry.feature, "cap-descriptors");
        assert_eq!(entry.status, FunctionStatus::Experimental);
        let callable = entry.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(callable.receiver, Some(BindingReceiver::Shared));
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
        assert_eq!(
            callable.error_type.unwrap().replace(' ', ""),
            "crate::DescriptorReadError"
        );
        assert_eq!(callable.output_type.replace(' ', ""), output);
        if parameter.is_empty() {
            assert!(callable.parameters.is_empty());
        } else {
            assert_eq!(callable.parameters.len(), 1);
            let p = callable.parameters[0];
            assert_eq!(p.name, parameter);
            assert_eq!(
                p.type_name,
                if parameter == "force" { "bool" } else { "u32" }
            );
            assert_eq!(
                p.default,
                if parameter == "force" {
                    BindingDefault::Value("false")
                } else {
                    BindingDefault::Required
                }
            );
        }
    }
}
