use cosmolkit_core::{BondDirectionStereoError, assign_chiral_types_from_bond_dirs};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer3D, SourceAtomValenceFacts,
    TopologyBlock,
};
use cosmolkit_types::{BondDirection, BondOrder, Element};

fn graph(facts: SourceAtomValenceFacts) -> TopologyBlock {
    let mut atoms = [Element::C, Element::C, Element::O, Element::CL]
        .into_iter()
        .enumerate()
        .map(|(i, element)| Atom::from_spec(AtomId::new(i), AtomSpec::new(element)))
        .collect::<Vec<_>>();
    atoms[1].set_source_valence_facts(facts);
    let bonds = [(1, 0), (1, 2), (1, 3)]
        .into_iter()
        .enumerate()
        .map(|(i, (a, b))| {
            let mut bond = Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
            );
            if i == 0 {
                bond.set_direction(BondDirection::BeginWedge);
            }
            bond
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
}
fn conformer() -> Conformer3D {
    Conformer3D::new(
        0,
        vec![
            [-3.9163, 5.4767, 0.0],
            [-3.9163, 3.9367, 0.0],
            [-2.5826, 3.1667, 0.0],
            [-5.25, 3.1667, 0.0],
        ],
        false,
    )
}

#[test]
fn source_wedge_updates_only_processed_atom_cache_and_recalculates_after_h_conversion() {
    let mut topology = graph(SourceAtomValenceFacts::UNINITIALIZED);
    assign_chiral_types_from_bond_dirs(&mut topology, &conformer(), false).unwrap();
    assert_eq!(topology.atoms[1].explicit_hydrogens(), 1);
    assert_eq!(
        topology.atoms[1].source_valence_facts(),
        SourceAtomValenceFacts {
            explicit_valence: 4,
            implicit_valence: 0
        }
    );
    for i in [0, 2, 3] {
        assert_eq!(
            topology.atoms[i].source_valence_facts(),
            SourceAtomValenceFacts::UNINITIALIZED
        );
    }
}

#[test]
fn source_wedge_reuses_initialized_count_instead_of_recomputing_chemistry() {
    let facts = SourceAtomValenceFacts {
        explicit_valence: 3,
        implicit_valence: 0,
    };
    let mut topology = graph(facts);
    assign_chiral_types_from_bond_dirs(&mut topology, &conformer(), false).unwrap();
    assert_eq!(topology.atoms[1].explicit_hydrogens(), 0);
    assert_eq!(topology.atoms[1].source_valence_facts(), facts);
    for i in [0, 2, 3] {
        assert_eq!(
            topology.atoms[i].source_valence_facts(),
            SourceAtomValenceFacts::UNINITIALIZED
        );
    }
}

#[test]
fn source_tiny_reference_vector_throws_after_cache_update_before_chiral_and_h_changes() {
    let mut topology = graph(SourceAtomValenceFacts::UNINITIALIZED);
    let coords = Conformer3D::new(
        0,
        vec![
            [0.0, 1e-100, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ],
        false,
    );
    let result = assign_chiral_types_from_bond_dirs(&mut topology, &coords, false);
    assert!(
        matches!(result, Err(BondDirectionStereoError::InvalidState(_))),
        "{result:?}"
    );
    assert_eq!(
        topology.atoms[1].source_valence_facts(),
        SourceAtomValenceFacts {
            explicit_valence: 3,
            implicit_valence: 1
        }
    );
    assert_eq!(topology.atoms[1].explicit_hydrogens(), 0);
    assert_eq!(
        topology.atoms[1].chiral_tag(),
        cosmolkit_types::ChiralTag::Unspecified
    );
    for i in [0, 2, 3] {
        assert_eq!(
            topology.atoms[i].source_valence_facts(),
            SourceAtomValenceFacts::UNINITIALIZED
        );
    }
}

#[test]
fn source_no_implicit_bypasses_negative_implicit_cache_without_recalculation() {
    let facts = SourceAtomValenceFacts {
        explicit_valence: 3,
        implicit_valence: -1,
    };
    let mut topology = graph(facts);
    topology.atoms[1].set_no_implicit(true);
    assign_chiral_types_from_bond_dirs(&mut topology, &conformer(), false).unwrap();
    assert_eq!(topology.atoms[1].source_valence_facts(), facts);
    assert_eq!(topology.atoms[1].explicit_hydrogens(), 0);
}

#[test]
fn source_existing_tag_skips_cache_and_unrelated_atoms_are_not_calculated() {
    let mut topology = graph(SourceAtomValenceFacts::UNINITIALIZED);
    topology.atoms[1].set_chiral_tag(cosmolkit_types::ChiralTag::TetrahedralCw);
    for i in [0, 2, 3] {
        topology.atoms[i].set_explicit_hydrogens(255);
    }
    let before = topology.clone();
    assign_chiral_types_from_bond_dirs(&mut topology, &conformer(), false).unwrap();
    assert_eq!(topology, before);
}

#[test]
fn source_h_conversion_strict_failure_retains_earlier_h_and_chiral_effects() {
    let facts = SourceAtomValenceFacts {
        explicit_valence: 3,
        implicit_valence: 1,
    };
    let mut topology = graph(facts);
    topology.bonds[1].set_order(BondOrder::Double);
    // C with a multiple bond yields source CHI_UNSPECIFIED; H conversion is
    // unconditional on this result. Adding H then causes strict valence 5.
    topology.atoms[1].set_chiral_tag(cosmolkit_types::ChiralTag::TetrahedralCw);
    let result = assign_chiral_types_from_bond_dirs(&mut topology, &conformer(), true);
    assert!(
        matches!(result, Err(BondDirectionStereoError::Valence(_))),
        "{result:?}"
    );
    assert_eq!(
        topology.atoms[1].chiral_tag(),
        cosmolkit_types::ChiralTag::Unspecified
    );
    assert_eq!(topology.atoms[1].explicit_hydrogens(), 1);
    assert_eq!(topology.atoms[1].source_valence_facts(), facts);
    for i in [0, 2, 3] {
        assert_eq!(
            topology.atoms[i].source_valence_facts(),
            SourceAtomValenceFacts::UNINITIALIZED
        );
    }
}
