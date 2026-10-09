//! Original configuration and factory conditions at the canonical value boundary.
use super::*;
#[test]
fn enumerator_configuration_defaults_match_every_source_default() {
    let p = TautomerParams::default();
    assert_eq!(p.max_tautomers(), 1000);
    assert_eq!(p.max_transforms(), 1000);
    assert!(p.remove_sp3_stereo());
    assert!(p.remove_bond_stereo());
    assert!(p.remove_isotopic_hydrogens());
    assert!(p.reassign_stereo());
    assert_eq!(p.transform_count(), 37);
    assert!(p.callback().is_none());
    assert!(p.scorer().is_none());
}
#[test]
fn enumerator_configuration_setters_cover_zero_maximum_and_every_boolean() {
    let mut p = TautomerParams::default();
    p.set_max_tautomers(0);
    p.set_max_transforms(u32::MAX);
    p.set_remove_sp3_stereo(false);
    p.set_remove_bond_stereo(false);
    p.set_remove_isotopic_hydrogens(false);
    p.set_reassign_stereo(false);
    assert_eq!(p.max_tautomers(), 0);
    assert_eq!(p.max_transforms(), u32::MAX);
    assert!(!p.remove_sp3_stereo());
    assert!(!p.remove_bond_stereo());
    assert!(!p.remove_isotopic_hydrogens());
    assert!(!p.reassign_stereo());
    let p = TautomerParams::default()
        .with_max_tautomers(u32::MAX)
        .with_max_transforms(0)
        .with_remove_sp3_stereo(false)
        .with_remove_bond_stereo(false)
        .with_remove_isotopic_hydrogens(false)
        .with_reassign_stereo(false);
    assert_eq!(p.max_tautomers(), u32::MAX);
    assert_eq!(p.max_transforms(), 0);
    assert!(!p.remove_sp3_stereo());
    assert!(!p.remove_bond_stereo());
    assert!(!p.remove_isotopic_hydrogens());
    assert!(!p.reassign_stereo());
}
#[test]
fn enumerator_configuration_clone_shares_catalog_but_not_option_state() {
    let original = TautomerParams::from_transform_data(&[("custom", "[O]-[C]", "=", "")])
        .unwrap()
        .with_max_tautomers(17);
    let mut clone = original.clone();
    assert!(Arc::ptr_eq(&original.catalog, &clone.catalog));
    assert_eq!(original.catalog.transforms().len(), 1);
    clone.set_max_tautomers(2);
    clone.set_remove_sp3_stereo(false);
    assert_eq!(original.max_tautomers(), 17);
    assert!(original.remove_sp3_stereo());
    assert_eq!(clone.max_tautomers(), 2);
    assert!(!clone.remove_sp3_stereo());
}
struct Decision(bool);
impl TautomerEnumerationCallback for Decision {
    fn should_continue(
        &self,
        _: TautomerMoleculeView<'_>,
        _: TautomerProgress<'_>,
    ) -> Result<bool, TautomerRunError> {
        Ok(self.0)
    }
}
#[test]
fn enumerator_configuration_clone_from_replaces_catalog_options_and_callback() {
    let callback: Arc<dyn TautomerEnumerationCallback> = Arc::new(Decision(true));
    let mut source = TautomerParams::default().with_max_transforms(9);
    source.set_callback(Some(callback.clone()));
    let mut target = TautomerParams::v1().unwrap();
    target.clone_from(&source);
    assert!(Arc::ptr_eq(&source.catalog, &target.catalog));
    assert_eq!(target.max_transforms(), 9);
    assert!(Arc::ptr_eq(target.callback().unwrap(), &callback));
}
#[test]
fn enumerator_configuration_callback_is_borrowed_and_replaceable() {
    // Canonical values own shared callbacks; the read accessor borrows that Arc.
    let first: Arc<dyn TautomerEnumerationCallback> = Arc::new(Decision(true));
    let second: Arc<dyn TautomerEnumerationCallback> = Arc::new(Decision(false));
    let mut p = TautomerParams::default();
    p.set_callback(Some(first.clone()));
    assert!(Arc::ptr_eq(p.callback().unwrap(), &first));
    p.set_callback(Some(second.clone()));
    assert!(Arc::ptr_eq(p.callback().unwrap(), &second));
    p.set_callback(None);
    assert!(p.callback().is_none());
}
#[test]
fn canonicalization_and_factories_share_current_and_v1_catalog_paths() {
    let mut source = Molecule::from_smiles("CC(C)=O").unwrap();
    let default = TautomerParams::default();
    let current = TautomerParams::from_transform_file("").unwrap();
    let v1 = TautomerParams::v1().unwrap();
    assert_eq!(default.transform_count(), 37);
    assert_eq!(current.transform_count(), 37);
    assert_eq!(v1.transform_count(), 36);
    assert_eq!(
        default.policy,
        cosmolkit_tautomer::TautomerParams::default()
    );
    assert_eq!(v1.policy, cosmolkit_tautomer::TautomerParams::default());
    let a = source.canonical_tautomer_with_params(&default).unwrap();
    let b = source.canonical_tautomer_with_params(&current).unwrap();
    let c = source.canonical_tautomer_with_params(&v1).unwrap();
    assert_eq!(a, b);
    assert_eq!(a.to_smiles().unwrap(), c.to_smiles().unwrap());
}
#[test]
fn canonicalization_and_factories_honor_options_without_mutating_configuration_or_source() {
    let mut source = Molecule::from_smiles("CC(C)=O")
        .unwrap()
        .to_builder()
        .with_property("source_id".into(), "canonical-options".into())
        .unwrap()
        .build()
        .unwrap();
    let before = source.clone();
    let p = TautomerParams::default()
        .with_max_transforms(0)
        .with_remove_sp3_stereo(false)
        .with_remove_bond_stereo(false)
        .with_remove_isotopic_hydrogens(false)
        .with_reassign_stereo(true);
    let options = p.policy;
    let selected = source.canonical_tautomer_with_params(&p).unwrap();
    assert_eq!(selected.to_smiles().unwrap().as_bytes(), b"CC(C)=O");
    assert_eq!(
        selected.property("source_id"),
        Some(&cosmolkit_model::PropertyValue::String(
            "canonical-options".into()
        ))
    );
    assert_eq!(p.policy, options);
    assert!(p.reassign_stereo());
    assert_eq!(source, before);
}
#[test]
fn canonicalization_and_factories_custom_scorer_selects_from_the_sole_enumeration() {
    struct Scorer(std::sync::Mutex<Vec<String>>);
    impl TautomerScorer for Scorer {
        fn score(&self, m: TautomerMoleculeView<'_>) -> Result<i32, TautomerRunError> {
            let key = m.to_smiles()?;
            self.0.lock().unwrap().push(fixed_key_text(&key).to_owned());
            Ok(if key.as_bytes() == b"C=C(C)O" {
                100
            } else {
                -100
            })
        }
    }
    let mut source = Molecule::from_smiles("CC(C)=O").unwrap();
    let before = source.clone();
    let scorer = Arc::new(Scorer(Default::default()));
    let mut p = TautomerParams::default();
    p.set_scorer(Some(scorer.clone()));
    let selected = source.canonical_tautomer_with_params(&p).unwrap();
    let mut calls = scorer.0.lock().unwrap().clone();
    calls.sort();
    assert_eq!(calls, ["C=C(C)O", "CC(C)=O"]);
    assert_eq!(selected.to_smiles().unwrap().as_bytes(), b"C=C(C)O");
    assert_eq!(
        selected.property("_StereochemDone"),
        Some(&cosmolkit_model::PropertyValue::Int(1))
    );
    assert_eq!(source, before);
}
#[test]
fn canonicalization_and_factories_are_independent_of_input_tautomer_and_atom_order() {
    let endpoints = ["CC(C)=O", "C=C(C)O", "O=C(C)C"].map(|s| {
        Molecule::from_smiles(s)
            .unwrap()
            .canonical_tautomer()
            .unwrap()
            .to_smiles()
            .unwrap()
    });
    assert!(endpoints.iter().all(|s| s == &endpoints[0]));
    assert_eq!(endpoints[0].as_bytes(), b"CC(C)=O");
}
#[test]
fn canonicalization_and_factories_canonical_selection_paths_follow_quickcopy_and_restoration() {
    // Canonical API replacement of the retired private in-place correspondence.
    // Exercise independently computed public selection paths on an enol that
    // actually changes bonds, retaining the original outer-state conditions.
    let molecule = Molecule::from_smiles("C=C(C)O").unwrap();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![crate::Conformer2D::new(
            0,
            (0..molecule.num_atoms()).map(|i| [i as f64, 0.5]).collect(),
        )],
        conformers_3d: vec![crate::Conformer3D::new(
            0,
            (0..molecule.num_atoms())
                .map(|i| [i as f64, i as f64 + 0.5, -(i as f64)])
                .collect(),
            true,
        )],
        source_coordinate_dim: None,
        source_conformer_order: Some(vec![
            cosmolkit_model::CoordinateDimension::TwoD,
            cosmolkit_model::CoordinateDimension::ThreeD,
        ]),
    };
    let props = molecule
        .properties()
        .clone()
        .with_name("in-place-source")
        .with_prop("source_id", "in-place-correspondence")
        .unwrap()
        .with_sdf_data_field("dataset", "tautomer");
    let mut source = crate::MoleculeBuilder::from_parts(
        molecule.topology().clone(),
        coordinates.clone(),
        props.clone(),
    )
    .build()
    .unwrap();
    let before = source.clone();
    let original_coordinates = source.coordinate_block_runtime().clone();
    let original_properties = source.properties().clone();
    let selected = source.canonical_tautomer().unwrap();
    assert_eq!(selected.to_smiles().unwrap().as_bytes(), b"CC(C)=O");
    assert_ne!(selected.bonds(), source.bonds());
    let mut enumeration = source.enumerate_tautomers().unwrap();
    assert_eq!(enumeration.len(), 2);
    let enumeration_before = enumeration.clone();
    let from_enumeration = enumeration.canonical_tautomer().unwrap();
    let mut candidates = enumeration.iter().cloned().collect::<Vec<_>>();
    let candidates_before = candidates.clone();
    let from_iterable = canonical_tautomer_from_molecules(&mut candidates).unwrap();
    // .6 canonicalize restores the complete original molecule dictionary and
    // appends source conformers with fresh global IDs. pickCanonical on the
    // enumeration/iterable has no original-input restoration step: the chosen
    // generated keto product retains quickCopy's empty conformer state and only
    // the computed properties produced by sanitization and stereo assignment.
    let expected_coordinates = CoordinateBlock {
        conformers_2d: vec![original_coordinates.conformers_2d[0].clone().with_id(0)],
        conformers_3d: vec![original_coordinates.conformers_3d[0].clone().with_id(1)],
        source_coordinate_dim: None,
        source_conformer_order: original_coordinates.source_conformer_order.clone(),
    };
    assert_eq!(selected.coordinate_block_runtime(), &expected_coordinates);
    assert_eq!(selected.properties(), &original_properties);
    let mut expected_product_properties = crate::MoleculeProperties::default();
    expected_product_properties
        .set_prop(
            "__computedProps",
            cosmolkit_model::PropertyValue::StringVector(Vec::new()),
        )
        .unwrap();
    expected_product_properties
        .set_computed_prop("numArom", 0_i32)
        .unwrap();
    expected_product_properties
        .set_computed_prop("_StereochemDone", 1_i32)
        .unwrap();
    for picked in [&from_enumeration, &from_iterable] {
        assert_eq!(
            picked.coordinate_block_runtime(),
            &CoordinateBlock::default()
        );
        assert_eq!(picked.properties(), &expected_product_properties);
    }
    for replacement in [&selected, &from_enumeration, &from_iterable] {
        assert_eq!(replacement.atoms(), selected.atoms());
        assert_eq!(replacement.bonds(), selected.bonds());
        assert!(
            replacement
                .derived_cache_runtime()
                .valence_assignment()
                .is_some()
        );
        assert_eq!(
            replacement.derived_cache_runtime().valence_assignment(),
            selected.derived_cache_runtime().valence_assignment()
        );
        assert!(!Arc::ptr_eq(
            &replacement.coordinates_arc_runtime(),
            &source.coordinates_arc_runtime()
        ));
    }
    assert_eq!(enumeration, enumeration_before);
    assert_eq!(candidates, candidates_before);
    assert_eq!(source, before);
    // The candidate operation writes its property block; unchanged input blocks
    // retain sharing with the original observer, independently of that output.
    assert!(Arc::ptr_eq(
        &source.properties_arc_runtime(),
        &before.properties_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.coordinates_arc_runtime(),
        &before.coordinates_arc_runtime()
    ));
}
#[test]
fn canonicalization_and_factories_do_not_expose_a_duplicate_in_place_api() {
    let source = include_str!("tautomer.rs");
    assert!(!source.contains("pub fn canonicalize_in_place"));
    assert!(!source.contains("fn canonicalize_in_place_compat_with"));
}
#[test]
fn enumeration_matches_pcs_fused_ring_max_transform_boundary() {
    let mut source = Molecule::from_smiles("Cc1nc2c(nc1C)C(=O)C1=C(C2=O)C2C=CC1CC2").unwrap();
    let result = source.enumerate_tautomers().unwrap();
    // Target .6 source keeps lazy Kekulize, preserved ring topology during
    // product sanitization, and sparse stereo cleanup; the .1 count was 272.
    // Keep the max-transform boundary and exact target count as separate checks.
    assert_eq!(result.len(), 259);
    assert_eq!(
        result.status(),
        TautomerEnumerationStatus::MaxTransformsReached
    );
}

fn fixed_key_text(key: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(key.as_bytes()).expect("original fixed ASCII test observation")
}

#[test]
fn search04_registered_enumeration_propagates_failed_lazy_attempt_without_mutating_input() {
    let mut source = Molecule::from_smiles_with_params(
        "c1cccc1",
        &crate::SmilesParseParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let before = source.clone();
    let topology = source.topology().clone();
    let properties = source.properties().clone();
    let coordinates = source.coordinate_block_runtime().clone();
    let cache = source.derived_cache_runtime().clone();
    let error = source.enumerate_tautomers().unwrap_err();
    assert!(matches!(
        error,
        OperationError::Tautomer(TautomerRunError::Kekulize(
            cosmolkit_core::KekulizeError::NotKekulizable { .. }
        ))
    ));
    assert_eq!(source, before);
    assert_eq!(source.topology(), &topology);
    assert_eq!(source.properties(), &properties);
    assert_eq!(source.coordinate_block_runtime(), &coordinates);
    assert_eq!(source.derived_cache_runtime(), &cache);
}
