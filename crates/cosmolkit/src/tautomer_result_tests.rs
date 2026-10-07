//! Six original rich-result value conditions through the canonical facade.
use super::*;
fn enumeration_result_fixture(status: TautomerEnumerationStatus) -> TautomerEnumeration {
    let mut entries = [("z", "CCC"), ("a", "C"), ("m", "CC")]
        .map(|(key, text)| {
            (
                cosmolkit_model::PropertyText::from(key),
                Molecule::from_smiles(text).unwrap(),
            )
        })
        .to_vec();
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    TautomerEnumeration {
        entries,
        status,
        modified_atoms: BTreeSet::from([AtomId::new(0), AtomId::new(2)]),
        modified_bonds: BTreeSet::from([BondId::new(1)]),
    }
}
#[test]
fn enumeration_result_default_is_empty_and_completed() {
    let result = TautomerEnumeration::default();

    assert!(result.is_empty());
    assert_eq!(result.len(), 0);
    assert_eq!(result.status(), TautomerEnumerationStatus::Completed);
    assert!(result.modified_atoms().is_empty());
    assert!(result.modified_bonds().is_empty());
    assert!(result.iter().next().is_none());
}

#[test]
fn enumeration_result_single_entry_supports_all_access_paths() {
    let molecule = Molecule::from_smiles("CCO").expect("parse sole tautomer");
    let result = TautomerEnumeration {
        entries: vec![("CCO".into(), molecule)],
        ..Default::default()
    };

    assert_eq!(result.len(), 1);
    assert!(!result.is_empty());
    assert_eq!(
        result
            .canonical_smiles()
            .into_iter()
            .map(fixed_key_text)
            .collect::<Vec<_>>(),
        ["CCO"]
    );
    assert_eq!(result.get(0).unwrap().num_atoms(), 3);
}

#[test]
fn enumeration_result_preserves_canonical_smiles_order_in_every_projection() {
    let result = enumeration_result_fixture(TautomerEnumerationStatus::Completed);

    assert_eq!(
        result
            .canonical_smiles()
            .into_iter()
            .map(fixed_key_text)
            .collect::<Vec<_>>(),
        ["a", "m", "z"]
    );
    assert_eq!(
        result
            .entries()
            .map(|(smiles, molecule)| (fixed_key_text(smiles), molecule.num_atoms()))
            .collect::<Vec<_>>(),
        [("a", 1), ("m", 2), ("z", 3)]
    );
    assert_eq!(
        result.iter().map(Molecule::num_atoms).collect::<Vec<_>>(),
        [1, 2, 3]
    );
    assert_eq!(
        result
            .iter()
            .rev()
            .map(Molecule::num_atoms)
            .collect::<Vec<_>>(),
        [3, 2, 1]
    );
    assert_eq!(
        (&result)
            .into_iter()
            .map(Molecule::num_atoms)
            .collect::<Vec<_>>(),
        [1, 2, 3]
    );
}

#[test]
fn enumeration_result_random_access_and_projections_cover_bounds() {
    let result = enumeration_result_fixture(TautomerEnumerationStatus::Completed);

    assert_eq!(result.get(0).unwrap().num_atoms(), 1);
    assert_eq!(result[1].num_atoms(), 2);
    assert_eq!(result.get(2).map(Molecule::num_atoms), Some(3));
    assert!(result.get(3).is_none());
    assert_eq!(result.get(3), None);
    assert_eq!(
        result.iter().map(Molecule::num_atoms).collect::<Vec<_>>(),
        [1, 2, 3]
    );
}

#[test]
fn enumeration_result_clone_has_independent_entries_and_typed_modified_sets() {
    let original = enumeration_result_fixture(TautomerEnumerationStatus::Completed);
    let mut cloned = original.clone();
    cloned.entries.remove(0);
    cloned.modified_atoms.insert(AtomId::new(9));
    cloned.modified_bonds.clear();

    assert_eq!(original.len(), 3);
    assert_eq!(cloned.len(), 2);
    assert_eq!(
        original.modified_atoms(),
        &BTreeSet::from([AtomId::new(0), AtomId::new(2)])
    );
    assert_eq!(original.modified_bonds(), &BTreeSet::from([BondId::new(1)]));
    assert!(cloned.modified_atoms().contains(&AtomId::new(9)));
    assert!(cloned.modified_bonds().is_empty());
}

#[test]
fn enumeration_result_exposes_every_source_status_without_reinterpretation() {
    for status in [
        TautomerEnumerationStatus::Completed,
        TautomerEnumerationStatus::MaxTautomersReached,
        TautomerEnumerationStatus::MaxTransformsReached,
        TautomerEnumerationStatus::Canceled,
    ] {
        assert_eq!(enumeration_result_fixture(status).status(), status);
    }
}

fn fixed_key_text(key: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(key.as_bytes()).expect("original fixed ASCII test observation")
}
