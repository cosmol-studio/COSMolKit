//! The original nine single-transform conditions at their detached owner boundary.
use super::*;
use crate::TautomerCatalog;
fn fixture(text: &str) -> Result<TautomerRecord, Box<dyn std::error::Error>> {
    super::stereo_tests::fixture_from_smiles(text)
}
fn builtin_transform(name: &str) -> TautomerTransform {
    crate::TautomerCatalog::current()
        .unwrap()
        .transforms()
        .iter()
        .find(|t| t.name().as_bytes() == name.as_bytes())
        .unwrap()
        .clone()
}
fn query(text: &str) -> cosmolkit_search::QueryGraph {
    cosmolkit_search::parse_smarts(text, &Default::default()).unwrap()
}
fn first_transform_match(
    candidate: &TautomerRecord,
    transform: &TautomerTransform,
) -> SubstructMatchResult {
    transform_matches(candidate, &CoordinateBlock::default(), transform)
        .unwrap()
        .into_iter()
        .next()
        .unwrap_or_else(|| panic!("{:?} must match focused fixture", transform.name()))
}
// Only adapt the original test's known-key set to the production borrowed lookup.
fn apply_single_transform(
    source: &TautomerRecord,
    candidate: &TautomerRecord,
    transform: &TautomerTransform,
    matched: &SubstructMatchResult,
    atoms: &BTreeSet<AtomId>,
    bonds: &BTreeSet<BondId>,
    keys: &BTreeSet<cosmolkit_model::PropertyText>,
    params: TautomerParams,
) -> Result<TautomerExpansionAttempt<Arc<TautomerRecord>>, TautomerRunError> {
    apply_tautomer_transform_match(
        source,
        candidate,
        &CoordinateBlock::default(),
        transform,
        matched,
        atoms,
        bonds,
        &|key| keys.contains(key),
        params,
    )
}
#[test]
fn single_transform_application_reproduces_every_builtin_bond_edit_shape() {
    let cases = [
        ("CC=O", "1,3 (thio)keto/enol f", "C=CO", ""),
        ("CC=C=O", "keten/ynol f", "CC#CO", "#-"),
        ("CC#CO", "keten/ynol r", "CC=C=O", "=="),
        (
            "NC(N)=S(=O)=O",
            "formamidinesulfinic acid f",
            "N=C(N)S(=O)O",
            "=--",
        ),
        (
            "NC(=N)S(=O)O",
            "formamidinesulfinic acid r",
            "NC(N)=S(=O)=O",
            "==-",
        ),
        ("C#N", "isocyanide f", "[C-]#[NH+]", "#"),
        ("OP(O)O", "phosphonic acid f", "O=[PH](O)O", "="),
        ("[PH](=O)(O)(O)", "phosphonic acid r", "OP(O)O", "-"),
    ];

    for (input, transform_name, expected_smiles, expected_bonds) in cases {
        let source = fixture(input).expect("parse edit-shape fixture");
        let source_before = source.clone();
        let candidate = kekulized(&source).expect("kekulize candidate");
        let candidate_before = candidate.clone();
        let transform = builtin_transform(transform_name);
        assert_eq!(
            transform.bond_types(),
            TautomerCatalog::from_data(&[("edit shape", "[C]", expected_bonds, "")])
                .unwrap()
                .transforms()[0]
                .bond_types(),
            "transform {transform_name}"
        );
        let matched = first_transform_match(&candidate, &transform);

        let attempt = apply_single_transform(
            &source,
            &candidate,
            &transform,
            &matched,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &BTreeSet::new(),
            TautomerParams::default().with_reassign_stereo(false),
        )
        .unwrap_or_else(|error| panic!("transform {transform_name}: {error}"));
        let TautomerExpansionAttempt::Product(product) = attempt else {
            panic!("transform {transform_name} did not produce a unique tautomer");
        };
        let molecule = product.tautomer.as_ref().clone();

        assert_eq!(
            product.canonical_smiles.as_bytes(),
            expected_smiles.as_bytes(),
            "{transform_name}"
        );
        assert_eq!(
            canonical_smiles(molecule.view(&CoordinateBlock::default()))
                .expect("write product canonical SMILES")
                .as_bytes(),
            expected_smiles.as_bytes(),
            "{transform_name}"
        );
        assert_eq!(source, source_before, "source changed for {transform_name}");
        assert_eq!(
            candidate, candidate_before,
            "candidate changed for {transform_name}"
        );
    }
}

#[test]
fn single_transform_application_uses_match_endpoint_order_and_unions_modified_sets() {
    let source = fixture("CC=O").expect("parse acetaldehyde");
    let candidate = kekulized(&source).expect("kekulize candidate");
    let transform = builtin_transform("1,3 (thio)keto/enol f");
    let matched = first_transform_match(&candidate, &transform);
    let donor = AtomId::new(matched.atom_mapping[0]);
    let acceptor = AtomId::new(*matched.atom_mapping.last().expect("acceptor"));
    let retained_atom = AtomId::new(matched.atom_mapping[1]);
    let retained_bond = BondId::new(matched.bond_mapping[0]);

    let attempt = apply_single_transform(
        &source,
        &candidate,
        &transform,
        &matched,
        &BTreeSet::from([retained_atom]),
        &BTreeSet::from([retained_bond]),
        &BTreeSet::new(),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .expect("apply ordered endpoint transform");
    let TautomerExpansionAttempt::Product(product) = attempt else {
        panic!("ordered endpoint transform must produce a product");
    };

    assert_eq!(
        product.tautomer.topology.atoms[donor.index()].explicit_hydrogens(),
        2
    );
    assert_eq!(
        product.tautomer.topology.atoms[acceptor.index()].explicit_hydrogens(),
        1
    );
    assert!(product.tautomer.topology.atoms[donor.index()].no_implicit());
    assert!(product.tautomer.topology.atoms[acceptor.index()].no_implicit());
    assert_eq!(
        product.modified_atoms,
        BTreeSet::from([donor, retained_atom, acceptor])
    );
    assert!(product.modified_bonds.contains(&retained_bond));
    assert!(
        matched
            .bond_mapping
            .iter()
            .all(|&bond| product.modified_bonds.contains(&BondId::new(bond)))
    );
}

#[test]
fn single_transform_application_handles_explicit_implicit_and_isotopic_hydrogens() {
    for input in ["CC=O", "[CH3]C=O"] {
        let source = fixture(input).expect("parse hydrogen fixture");
        let candidate = kekulized(&source).expect("kekulize candidate");
        let transform = builtin_transform("1,3 (thio)keto/enol f");
        let matched = first_transform_match(&candidate, &transform);
        let attempt = apply_single_transform(
            &source,
            &candidate,
            &transform,
            &matched,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &BTreeSet::new(),
            TautomerParams::default().with_reassign_stereo(false),
        )
        .expect("apply implicit/explicit hydrogen transform");
        let TautomerExpansionAttempt::Product(product) = attempt else {
            panic!("hydrogen transform must produce a product");
        };
        assert_eq!(product.canonical_smiles.as_bytes(), b"C=CO");
    }

    let mut source = fixture("CC=O").expect("parse isotopic-H fixture");
    source.topology.atoms[0].set_tracked_isotopic_hydrogens(vec![2, 3]);
    let candidate = kekulized(&source).expect("kekulize candidate");
    let transform = builtin_transform("1,3 (thio)keto/enol f");
    let matched = first_transform_match(&candidate, &transform);
    for (remove, expected) in [(true, &[][..]), (false, &[2, 3][..])] {
        let attempt = apply_single_transform(
            &source,
            &candidate,
            &transform,
            &matched,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &BTreeSet::new(),
            TautomerParams::default()
                .with_remove_isotopic_hydrogens(remove)
                .with_reassign_stereo(false),
        )
        .expect("apply isotopic-H transform");
        let TautomerExpansionAttempt::Product(product) = attempt else {
            panic!("isotopic-H transform must produce a product");
        };
        assert_eq!(
            product.tautomer.topology.atoms[0].tracked_isotopic_hydrogens(),
            expected
        );
    }
}

#[test]
fn single_transform_application_applies_source_ordered_charge_deltas() {
    let source = fixture("C#N").expect("parse hydrogen cyanide");
    let candidate = kekulized(&source).expect("kekulize candidate");
    let transform = builtin_transform("isocyanide f");
    let matched = first_transform_match(&candidate, &transform);
    assert_eq!(transform.charges(), &[-1, 1]);

    let attempt = apply_single_transform(
        &source,
        &candidate,
        &transform,
        &matched,
        &BTreeSet::new(),
        &BTreeSet::new(),
        &BTreeSet::new(),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .expect("apply isocyanide charge transform");
    let TautomerExpansionAttempt::Product(product) = attempt else {
        panic!("isocyanide transform must produce a product");
    };
    assert_eq!(
        product.tautomer.topology.atoms[matched.atom_mapping[0]].formal_charge(),
        -1
    );
    assert_eq!(
        product.tautomer.topology.atoms[matched.atom_mapping[1]].formal_charge(),
        1
    );
    assert_eq!(product.canonical_smiles.as_bytes(), b"[C-]#[NH+]");
}

#[test]
fn single_transform_application_records_bonds_changed_only_by_sanitization() {
    let source = fixture("Cc1nc2ccccc2[nH]1").expect("parse source-comment sanitization fixture");
    let source_before = source.clone();
    let candidate = kekulized(&source).expect("kekulize candidate");
    let catalog = TautomerCatalog::current().expect("load catalog");
    let mut witnessed = None;

    'transforms: for transform in catalog.transforms() {
        for matched in
            transform_matches(&candidate, &CoordinateBlock::default(), transform).unwrap()
        {
            let directly_edited = matched
                .bond_mapping
                .iter()
                .copied()
                .map(BondId::new)
                .collect::<BTreeSet<_>>();
            let attempt = apply_single_transform(
                &source,
                &candidate,
                transform,
                &matched,
                &BTreeSet::new(),
                &BTreeSet::new(),
                &BTreeSet::new(),
                TautomerParams::default().with_reassign_stereo(false),
            )
            .unwrap_or_else(|error| panic!("{:?}: {error}", transform.name()));
            if let TautomerExpansionAttempt::Product(product) = attempt {
                let sanitize_only = product
                    .modified_bonds
                    .difference(&directly_edited)
                    .copied()
                    .collect::<BTreeSet<_>>();
                if !sanitize_only.is_empty() {
                    witnessed = Some((product, sanitize_only));
                    break 'transforms;
                }
            }
        }
    }

    let (product, sanitize_only) = witnessed
        .expect("source-comment fixture must expose a sanitization-only bond-order change");
    assert!(sanitize_only.iter().all(|bond| {
        source.topology.bonds.as_slice()[bond.index()].order()
            != product.tautomer.topology.bonds[bond.index()].order()
    }));
    assert_eq!(
        canonical_smiles(product.tautomer.view(&CoordinateBlock::default())).unwrap(),
        product.canonical_smiles
    );
    assert_eq!(source, source_before);
}

#[test]
fn single_transform_application_reports_duplicate_after_stereo_and_before_product_kekulization() {
    let source = fixture("CC=O").expect("parse duplicate fixture");
    let candidate = kekulized(&source).expect("kekulize candidate");
    let transform = builtin_transform("1,3 (thio)keto/enol f");
    let matched = first_transform_match(&candidate, &transform);
    let first = apply_single_transform(
        &source,
        &candidate,
        &transform,
        &matched,
        &BTreeSet::new(),
        &BTreeSet::new(),
        &BTreeSet::new(),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .expect("compute duplicate key");
    let TautomerExpansionAttempt::Product(first) = first else {
        panic!("first transform must be unique");
    };

    let duplicate = apply_single_transform(
        &source,
        &candidate,
        &transform,
        &matched,
        &BTreeSet::new(),
        &BTreeSet::new(),
        &BTreeSet::from([first.canonical_smiles.clone()]),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .expect("duplicate is a non-error branch");
    let TautomerExpansionAttempt::Duplicate {
        canonical_smiles,
        modified_atoms,
        modified_bonds,
    } = duplicate
    else {
        panic!("existing canonical key must select duplicate branch");
    };
    assert_eq!(canonical_smiles, first.canonical_smiles);
    assert!(!modified_atoms.is_empty());
    assert!(!modified_bonds.is_empty());
}

#[test]
fn single_transform_application_catches_only_the_kekulize_failure_branch() {
    let parsed = cosmolkit_smiles::parse_smiles(
        "c1cccc1.CC",
        &cosmolkit_smiles::SmilesParseParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .expect("parse odd aromatic ring fixture");
    let source = prepared(TautomerRecordView {
        topology: &parsed.topology,
        coordinates: &parsed.coordinates,
        properties: &parsed.properties,
        valence: None,
        rings: None,
    })
    .expect("prepare non-strict valence cache");
    let candidate = source.clone();
    let transform = TautomerTransform::new(
        "focused kekulize failure",
        query("[C]-[C]"),
        Vec::new(),
        Vec::new(),
    )
    .expect("construct focused transform");
    let matched = first_transform_match(&candidate, &transform);
    let source_before = source.clone();
    let candidate_before = candidate.clone();

    let attempt = apply_single_transform(
        &source,
        &candidate,
        &transform,
        &matched,
        &BTreeSet::new(),
        &BTreeSet::new(),
        &BTreeSet::new(),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .expect("Kekulize failure is the source-defined recoverable branch");
    let TautomerExpansionAttempt::RecoverableKekulizeFailure {
        modified_atoms,
        modified_bonds,
    } = attempt
    else {
        panic!("odd aromatic ring must fail in partial-sanitize Kekulize stage");
    };
    assert_eq!(modified_atoms.len(), 2);
    assert_eq!(modified_bonds.len(), 1);
    assert_eq!(source, source_before);
    assert_eq!(candidate, candidate_before);
}

#[test]
fn source_property_failure_precedes_transform_sanitize_and_preserves_inputs() {
    let source = fixture("CC=O").unwrap();
    let mut candidate = kekulized(&source).unwrap();
    let transform = builtin_transform("1,3 (thio)keto/enol f");
    let matched = first_transform_match(&candidate, &transform);
    candidate
        .properties
        .set_prop("__computedProps", cosmolkit_model::PropertyValue::Int(7))
        .unwrap();
    candidate.topology.atoms[0]
        .set_prop("__computedProps", cosmolkit_model::PropertyValue::Int(8))
        .unwrap();
    let source_before = source.clone();
    let candidate_before = candidate.clone();
    let run = |candidate: &TautomerRecord| {
        apply_single_transform(
            &source,
            candidate,
            &transform,
            &matched,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &BTreeSet::new(),
            TautomerParams::default().with_reassign_stereo(false),
        )
    };
    let error = run(&candidate).unwrap_err();
    assert!(matches!(
        error,
        TautomerRunError::MoleculeProperty(
            cosmolkit_model::MoleculePropertyError::ComputedListKind(_)
        )
    ));
    assert_eq!(source, source_before);
    assert_eq!(candidate, candidate_before);

    // Fix only the molecule field: the original competing atom error remains
    // observable at the next source stage, rather than being swallowed.
    candidate
        .properties
        .set_prop(
            "__computedProps",
            cosmolkit_model::PropertyValue::StringVector(Vec::new()),
        )
        .unwrap();
    let before_atom_failure = candidate.clone();
    let error = run(&candidate).unwrap_err();
    assert!(matches!(
        error,
        TautomerRunError::Sanitize(SanitizeError::AtomProperty(
            cosmolkit_model::AtomPropertyError::ComputedListKind(_)
        ))
    ));
    assert_eq!(candidate, before_atom_failure);

    candidate.topology.atoms[0]
        .set_prop(
            "__computedProps",
            cosmolkit_model::PropertyValue::StringVector(Vec::new()),
        )
        .unwrap();
    let before_success = candidate.clone();
    let TautomerExpansionAttempt::Product(product) = run(&candidate).unwrap() else {
        panic!("valid source property state must preserve the original keto/enol product");
    };
    assert_eq!(product.canonical_smiles.as_bytes(), b"C=CO");
    assert_eq!(candidate, before_success);
    assert_eq!(source, source_before);
}

#[test]
fn single_transform_application_propagates_non_kekulize_sanitize_failures() {
    let source = fixture("CCC").expect("parse propagation fixture");
    let mut candidate = source.clone();
    candidate.topology.bonds[1].set_order(BondOrder::Other);
    let transform = TautomerTransform::new(
        "focused sanitize propagation",
        query("[C]-[C]"),
        Vec::new(),
        Vec::new(),
    )
    .expect("construct focused transform");
    let matched = SubstructMatchResult {
        atom_mapping: vec![0, 1],
        bond_mapping: vec![0],
    };
    let source_before = source.clone();
    let candidate_before = candidate.clone();

    let error = apply_single_transform(
        &source,
        &candidate,
        &transform,
        &matched,
        &BTreeSet::new(),
        &BTreeSet::new(),
        &BTreeSet::new(),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .expect_err("bad unrelated bond type must propagate from property-cache sanitization");
    assert!(matches!(
        error,
        TautomerRunError::Sanitize(cosmolkit_core::SanitizeError::Properties {
            stage: cosmolkit_core::SanitizeStage::Properties,
            ..
        })
    ));
    assert_eq!(source, source_before);
    assert_eq!(candidate, candidate_before);
}

#[test]
fn single_transform_application_rejects_invalid_mappings_without_mutating_inputs() {
    let source = fixture("CC=O").expect("parse mapping fixture");
    let candidate = kekulized(&source).expect("kekulize candidate");
    let transform = builtin_transform("1,3 (thio)keto/enol f");
    let source_before = source.clone();
    let candidate_before = candidate.clone();
    let error = apply_single_transform(
        &source,
        &candidate,
        &transform,
        &SubstructMatchResult {
            atom_mapping: vec![0, 1],
            bond_mapping: vec![0, 1],
        },
        &BTreeSet::new(),
        &BTreeSet::new(),
        &BTreeSet::new(),
        TautomerParams::default(),
    )
    .expect_err("short atom mapping must return a structured error");

    assert!(matches!(
        error,
        TautomerRunError::AtomMappingCount {
            expected: 3,
            actual: 2
        }
    ));
    assert_eq!(source, source_before);
    assert_eq!(candidate, candidate_before);
}
