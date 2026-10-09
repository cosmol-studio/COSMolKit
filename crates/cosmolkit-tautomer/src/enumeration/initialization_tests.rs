//! Original seven initialization conditions, adapted to explicit detached data.
use super::*;
use cosmolkit_core::{RingFindType, RingInfo};
use cosmolkit_types::BondOrder;
fn fixture(text: &str) -> TautomerRecord {
    crate::engine::stereo_tests::fixture_from_smiles(text).unwrap()
}
fn initialize(
    record: &TautomerRecord,
) -> (TautomerRecord, TautomerExpansionState<Arc<TautomerRecord>>) {
    let coordinates = CoordinateBlock::default();
    let view = record.view(&coordinates);
    initialize_candidates_with_key(view, canonical_smiles(view).unwrap()).unwrap()
}
#[test]
fn enumeration_initialization_updates_only_missing_valence_and_symm_sssr_caches() {
    let parsed = cosmolkit_smiles::parse_smiles(
        "c1ccccc1",
        &cosmolkit_smiles::SmilesParseParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let view = TautomerRecordView {
        topology: &parsed.topology,
        coordinates: &parsed.coordinates,
        properties: &parsed.properties,
        valence: None,
        rings: None,
    };
    let (uncached, _) = initialize_candidates_with_key(view, "c1ccccc1".into()).unwrap();
    assert_eq!(uncached.valence.explicit_valence.len(), 6);
    assert!(uncached.rings.is_symm_sssr());
    let cached = fixture("c1ccccc1");
    let (retained, _) = initialize(&cached);
    assert_eq!(retained.valence, cached.valence);
    assert_eq!(retained.rings, cached.rings);
    let mut fast = cached.clone();
    fast.rings = RingInfo::new(
        RingFindType::Fast,
        fast.topology.atoms.len(),
        fast.topology.bonds.len(),
    );
    assert!(initialize(&fast).0.rings.is_symm_sssr());
}
#[test]
fn enumeration_initialization_canonical_kekulizes_aromatic_input_without_clearing_flags() {
    let (_, mut state) = initialize(&fixture("c1ccccc1"));
    assert!(state.candidates[b"c1ccccc1".as_slice()].kekulized.is_none());
    get_cached_kekulized(state.candidates.get_mut(b"c1ccccc1".as_slice()).unwrap()).unwrap();
    assert_eq!(
        state
            .candidates
            .keys()
            .map(fixed_key_text)
            .collect::<Vec<_>>(),
        ["c1ccccc1"]
    );
    let topology = &state.candidates[b"c1ccccc1".as_slice()]
        .kekulized
        .as_ref()
        .unwrap()
        .topology;
    assert!(
        topology
            .atoms
            .iter()
            .all(cosmolkit_model::Atom::is_aromatic)
    );
    assert!(
        topology
            .bonds
            .iter()
            .all(cosmolkit_model::Bond::is_aromatic)
    );
    assert_eq!(
        topology
            .bonds
            .iter()
            .filter(|b| b.order() == BondOrder::Double)
            .count(),
        3
    );
    assert_eq!(
        topology
            .bonds
            .iter()
            .filter(|b| b.order() == BondOrder::Single)
            .count(),
        3
    );
}
#[test]
fn enumeration_initialization_handles_empty_and_disconnected_inputs() {
    let empty = cosmolkit_model::TopologyBlock::default();
    let coordinates = CoordinateBlock::default();
    let properties = cosmolkit_model::MoleculeProperties::default();
    let (empty, state) = initialize_candidates_with_key(
        TautomerRecordView {
            topology: &empty,
            coordinates: &coordinates,
            properties: &properties,
            valence: None,
            rings: None,
        },
        cosmolkit_model::PropertyText::new(),
    )
    .unwrap();
    assert_eq!(empty.topology.atoms.len(), 0);
    assert_eq!(empty.topology.bonds.len(), 0);
    assert!(state.candidates.contains_key(b"".as_slice()));
    let (disconnected, state) = initialize(&fixture("O.CC"));
    assert!(state.candidates.contains_key(b"CC.O".as_slice()));
    assert_eq!(disconnected.topology.atoms.len(), 3);
    assert_eq!(disconnected.topology.bonds.len(), 1);
}
#[test]
fn enumeration_initialization_is_independent_of_input_atom_order() {
    let first = initialize(&fixture("CC(=O)O")).1;
    let second = initialize(&fixture("OC(=O)C")).1;
    assert_eq!(
        first
            .candidates
            .keys()
            .map(fixed_key_text)
            .collect::<Vec<_>>(),
        ["CC(=O)O"]
    );
    assert_eq!(
        first.candidates.keys().collect::<Vec<_>>(),
        second.candidates.keys().collect::<Vec<_>>()
    );
}
#[test]
fn enumeration_initialization_inserts_one_unfinished_candidate_in_key_order() {
    let (_, state) = initialize(&fixture("OC(=O)C"));
    assert_eq!(state.candidates.len(), 1);
    let value = &state.candidates[b"CC(=O)O".as_slice()];
    assert!(value.tautomer.is_some());
    assert!(value.kekulized.is_none());
    assert_eq!(value.num_modified_atoms, 0);
    assert_eq!(value.num_modified_bonds, 0);
    assert!(!value.done);
    assert!(state.modified_atoms.is_empty());
    assert!(state.modified_bonds.is_empty());
    assert_eq!(state.num_transforms, 0);
    assert_eq!(state.status, TautomerEnumerationStatus::Completed);
}
#[test]
fn enumeration_initialization_reports_invalid_aromatic_input_as_a_structured_error() {
    let parsed = cosmolkit_smiles::parse_smiles(
        "c1cccc1",
        &cosmolkit_smiles::SmilesParseParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let (_, mut state) = initialize_candidates_with_key(
        TautomerRecordView {
            topology: &parsed.topology,
            coordinates: &parsed.coordinates,
            properties: &parsed.properties,
            valence: None,
            rings: None,
        },
        "c1cccc1".into(),
    )
    .unwrap();
    let error =
        get_cached_kekulized(state.candidates.get_mut(b"c1cccc1".as_slice()).unwrap()).unwrap_err();
    assert!(matches!(error, TautomerRunError::Kekulize(_)));
}
#[test]
fn enumeration_initialization_never_mutates_the_input() {
    for text in ["c1ccccc1", "CC(=O)O", "O.CC"] {
        let source = fixture(text);
        let before = source.clone();
        let _ = initialize(&source);
        assert_eq!(source, before, "{text}");
    }
}

fn fixed_key_text(key: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(key.as_bytes()).expect("original fixed ASCII test key")
}
