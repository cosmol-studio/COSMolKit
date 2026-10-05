#![cfg(feature = "cap-fingerprints")]
use cosmolkit::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, ChiralTag, CoordinateBlock, Element,
    Molecule, MoleculeProperties, OperationError, PropertyValue, TopologyBlock,
};

fn tetrahedron(tag: ChiralTag, computed: Option<bool>) -> Molecule {
    let atoms = [Element::C, Element::F, Element::CL, Element::BR, Element::I]
        .into_iter()
        .enumerate()
        .map(|(i, element)| {
            let spec = AtomSpec::new(element);
            Atom::from_spec(
                AtomId::new(i),
                if i == 0 {
                    spec.with_chiral_tag(tag)
                } else {
                    spec
                },
            )
        })
        .collect();
    let bonds = (1..5)
        .map(|i| {
            Bond::from_spec(
                BondId::new(i - 1),
                BondSpec::new(AtomId::new(0), AtomId::new(i), BondOrder::Single),
            )
        })
        .collect();
    let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
    let mut properties = MoleculeProperties::default()
        .with_prop("user", "preserved")
        .unwrap();
    if let Some(value) = computed {
        properties
            .set_prop("_CIPComputed", if value { "1" } else { "0" })
            .unwrap();
    }
    Molecule::from_parts(topology, CoordinateBlock::default(), properties).unwrap()
}

fn shared(left: &Molecule, right: &Molecule) {
    assert!(std::ptr::eq(left.topology(), right.topology()));
    assert!(std::ptr::eq(left.properties(), right.properties()));
}

#[test]
fn legacy_and_no_chirality_leave_cip_absent_and_share_blocks() {
    let source = tetrahedron(ChiralTag::TetrahedralCcw, None)
        .with_assigned_valence()
        .unwrap();
    for (chiral, legacy) in [(false, true), (false, false), (true, true)] {
        let output = source
            .with_atom_pair_atom_code(AtomId::new(0), 0, chiral, legacy)
            .unwrap();
        assert_eq!(output.code, 36);
        assert!(output.molecule.properties().prop("_CIPComputed").is_none());
        shared(&source, &output.molecule);
    }
}

#[test]
fn modern_missing_marker_installs_owner_assignment_atomically() {
    let source = tetrahedron(ChiralTag::TetrahedralCcw, None)
        .with_assigned_valence()
        .unwrap();
    let peer = source.clone();
    let output = source
        .with_atom_pair_atom_code(AtomId::new(0), 0, true, false)
        .unwrap();
    // Pinned CIP/source getAtomCode: tetrahedral Ccw center here is S, code36|1024.
    assert_eq!(output.code, 1060);
    assert_eq!(
        output
            .molecule
            .atom(AtomId::new(0))
            .unwrap()
            .prop("_CIPCode"),
        Some(&PropertyValue::from("S"))
    );
    assert_eq!(output.molecule.properties().prop("_CIPComputed"), Some("1"));
    assert!(
        output
            .molecule
            .properties()
            .is_prop_computed("_CIPComputed")
    );
    assert_eq!(output.molecule.property("user"), Some("preserved"));
    assert!(!std::ptr::eq(source.topology(), output.molecule.topology()));
    assert!(!std::ptr::eq(
        source.properties(),
        output.molecule.properties()
    ));
    shared(&source, &peer);
    assert_eq!(source, peer);
    assert!(source.properties().prop("_CIPComputed").is_none());
    let repeated = output
        .molecule
        .with_atom_pair_atom_code(AtomId::new(0), 0, true, false)
        .unwrap();
    assert_eq!(repeated.code, 1060);
    shared(&output.molecule, &repeated.molecule);
}

#[test]
fn present_false_and_true_use_presence_not_truth_and_keep_ordinary_marker() {
    for value in [false, true] {
        let source = tetrahedron(ChiralTag::TetrahedralCcw, Some(value))
            .with_assigned_valence()
            .unwrap();
        let result = source
            .with_atom_pair_atom_code(AtomId::new(0), 0, true, false)
            .unwrap();
        assert_eq!(result.code, 36);
        assert_eq!(
            result.molecule.properties().prop("_CIPComputed"),
            Some(if value { "1" } else { "0" })
        );
        assert!(
            !result
                .molecule
                .properties()
                .is_prop_computed("_CIPComputed")
        );
        shared(&source, &result.molecule);
    }
}

#[test]
fn modern_unspecified_atom_does_not_trigger_molecule_cip_assignment() {
    let source = tetrahedron(ChiralTag::Unspecified, None)
        .with_assigned_valence()
        .unwrap();
    let result = source
        .with_atom_pair_atom_code(AtomId::new(0), 0, true, false)
        .unwrap();
    assert_eq!(result.code, 36);
    assert!(result.molecule.properties().prop("_CIPComputed").is_none());
    shared(&source, &result.molecule);
}

#[test]
fn num_pi_cache_error_precedes_modern_cip_and_failure_keeps_input() {
    let source = tetrahedron(ChiralTag::TetrahedralCcw, None);
    let peer = source.clone();
    let error = source
        .with_atom_pair_atom_code(AtomId::new(0), u32::MAX, true, false)
        .unwrap_err();
    assert!(matches!(
        error,
        OperationError::AtomCode(cosmolkit_fingerprints::AtomCodeError::Valence(
            cosmolkit_core::ValenceError::PiElectronExplicitValenceCacheNotInitialized { .. }
        ))
    ));
    assert_eq!(source, peer);
    shared(&source, &peer);
    assert!(source.properties().prop("_CIPComputed").is_none());
    assert!(
        source
            .atom(AtomId::new(0))
            .unwrap()
            .prop("_CIPCode")
            .is_none()
    );
}

#[test]
fn branch_subtract_unsigned_extreme_and_atom_bounds_preserve_source_order() {
    let source = tetrahedron(ChiralTag::TetrahedralCcw, None)
        .with_assigned_valence()
        .unwrap();
    assert_eq!(
        source
            .with_atom_pair_atom_code(AtomId::new(0), u32::MAX, true, false)
            .unwrap()
            .code,
        1056
    );
    let error = source
        .with_atom_pair_atom_code(AtomId::new(5), 0, true, false)
        .unwrap_err();
    assert!(matches!(
        error,
        OperationError::AtomCode(cosmolkit_fingerprints::AtomCodeError::Valence(
            cosmolkit_core::ValenceError::AtomOutOfRange { .. }
        ))
    ));
    assert!(source.properties().prop("_CIPComputed").is_none());
}
