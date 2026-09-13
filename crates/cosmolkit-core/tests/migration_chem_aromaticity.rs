use cosmolkit_core::{
    AromaticityError, AromaticityModel, AromaticityParams, RingFindType, RingInfo,
    RingSearchParams, assign_aromaticity, symmetrized_sssr,
};
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
use cosmolkit_types::{BondOrder, Element, Hybridization};

fn topology(specs: Vec<AtomSpec>, edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
    let atoms = specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
        .collect();
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(index, &(begin, end, order))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn cycle(specs: Vec<AtomSpec>) -> TopologyBlock {
    let count = specs.len();
    let edges = (0..count)
        .map(|index| {
            (
                index,
                (index + 1) % count,
                if index % 2 == 0 {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                },
            )
        })
        .collect::<Vec<_>>();
    topology(specs, &edges)
}

fn carbon_cycle(size: usize) -> TopologyBlock {
    cycle((0..size).map(|_| AtomSpec::new(Element::C)).collect())
}

fn sp2_carbon_cycle(size: usize) -> TopologyBlock {
    cycle(
        (0..size)
            .map(|_| AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2))
            .collect(),
    )
}

fn rings(topology: &TopologyBlock) -> RingInfo {
    symmetrized_sssr(topology, &RingSearchParams::default()).unwrap()
}

fn assignment(
    topology: &TopologyBlock,
    model: AromaticityModel,
) -> cosmolkit_core::AromaticityAssignment {
    assign_aromaticity(topology, &rings(topology), &AromaticityParams { model }).unwrap()
}

fn aromatic_atom_count(topology: &TopologyBlock) -> usize {
    topology
        .atoms
        .iter()
        .filter(|atom| atom.is_aromatic())
        .count()
}

fn aromatic_bond_count(topology: &TopologyBlock) -> usize {
    topology
        .bonds
        .iter()
        .filter(|bond| bond.is_aromatic())
        .count()
}

#[test]
fn defaults_empty_acyclic_custom_and_dimension_errors_are_exact_and_atomic() {
    assert_eq!(
        AromaticityParams::default(),
        AromaticityParams {
            model: AromaticityModel::Rdkit,
        }
    );

    let empty = TopologyBlock::default();
    let empty_rings = rings(&empty);
    let empty_result =
        assign_aromaticity(&empty, &empty_rings, &AromaticityParams::default()).unwrap();
    assert_eq!(empty_result.topology, empty);
    assert_eq!(empty_result.aromatic_ring_count, 0);

    let acyclic = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::O)],
        &[(0, 1, BondOrder::Single)],
    );
    let acyclic_snapshot = acyclic.clone();
    let acyclic_result = assignment(&acyclic, AromaticityModel::Rdkit);
    assert_eq!(acyclic_result.topology, acyclic_snapshot);
    assert_eq!(acyclic_result.aromatic_ring_count, 0);
    assert_eq!(acyclic, acyclic_snapshot);

    let custom_snapshot = acyclic.clone();
    assert!(matches!(
        assign_aromaticity(
            &acyclic,
            &rings(&acyclic),
            &AromaticityParams {
                model: AromaticityModel::Custom,
            },
        ),
        Err(AromaticityError::UnsupportedModel {
            model: AromaticityModel::Custom,
            ..
        })
    ));
    assert_eq!(acyclic, custom_snapshot);

    let wrong_dimensions = RingInfo::new(RingFindType::Sssr, 3, 1);
    assert!(matches!(
        assign_aromaticity(&acyclic, &wrong_dimensions, &AromaticityParams::default()),
        Err(AromaticityError::RingInfoDimensionMismatch {
            ring_atom_count: 3,
            ring_bond_count: 1,
            topology_atom_count: 2,
            topology_bond_count: 1,
        })
    ));

    let mut malformed = acyclic.clone();
    malformed.atoms[0] = Atom::from_spec(AtomId::new(7), AtomSpec::new(Element::C));
    assert!(matches!(
        assign_aromaticity(
            &malformed,
            &RingInfo::new(RingFindType::Sssr, 2, 1),
            &AromaticityParams::default(),
        ),
        Err(AromaticityError::InvalidTopology(_))
    ));
}

#[test]
fn simple_model_covers_benzene_pyrrole_and_four_six_eight_member_boundaries() {
    let benzene = carbon_cycle(6);
    let benzene_result = assignment(&benzene, AromaticityModel::Simple);
    assert_eq!(benzene_result.aromatic_ring_count, 1);
    assert_eq!(aromatic_atom_count(&benzene_result.topology), 6);
    assert_eq!(aromatic_bond_count(&benzene_result.topology), 6);

    let pyrrole = topology(
        vec![
            AtomSpec::new(Element::N).with_explicit_hydrogens(1),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Double),
            (4, 0, BondOrder::Single),
        ],
    );
    let pyrrole_result = assignment(&pyrrole, AromaticityModel::Simple);
    assert_eq!(pyrrole_result.aromatic_ring_count, 1);
    assert_eq!(aromatic_atom_count(&pyrrole_result.topology), 5);

    for size in [4, 8] {
        let graph = carbon_cycle(size);
        let result = assignment(&graph, AromaticityModel::Simple);
        assert_eq!(result.aromatic_ring_count, 0, "ring size {size}");
        assert_eq!(aromatic_atom_count(&result.topology), 0, "ring size {size}");
        assert_eq!(aromatic_bond_count(&result.topology), 0, "ring size {size}");
    }
}

#[test]
fn azulene_distinguishes_rdkit_fused_simple_and_mmff94_policies() {
    let graph = topology(
        vec![AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2); 10],
        &[
            (0, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Double),
            (4, 1, BondOrder::Single),
            (1, 0, BondOrder::Single),
            (0, 5, BondOrder::Single),
            (5, 6, BondOrder::Double),
            (6, 7, BondOrder::Single),
            (7, 8, BondOrder::Double),
            (8, 9, BondOrder::Single),
            (9, 1, BondOrder::Double),
        ],
    );
    let ring_info = rings(&graph);
    assert_eq!(ring_info.num_rings(), 2);

    let rdkit = assign_aromaticity(
        &graph,
        &ring_info,
        &AromaticityParams {
            model: AromaticityModel::Rdkit,
        },
    )
    .unwrap();
    assert_eq!(rdkit.aromatic_ring_count, 2);
    assert_eq!(aromatic_bond_count(&rdkit.topology), 10);
    assert!(!rdkit.topology.bonds[4].is_aromatic());

    let simple = assign_aromaticity(
        &graph,
        &ring_info,
        &AromaticityParams {
            model: AromaticityModel::Simple,
        },
    )
    .unwrap();
    assert_eq!(simple.aromatic_ring_count, 0);
    assert_eq!(aromatic_bond_count(&simple.topology), 0);

    let mmff = assign_aromaticity(
        &graph,
        &ring_info,
        &AromaticityParams {
            model: AromaticityModel::Mmff94,
        },
    )
    .unwrap();
    assert_eq!(aromatic_bond_count(&mmff.topology), 0);
}

#[test]
fn mdl_positive_negative_and_mixed_rows_match_the_pinned_matrix() {
    let benzene = carbon_cycle(6);
    let positive = assignment(&benzene, AromaticityModel::Mdl);
    assert_eq!(positive.aromatic_ring_count, 1);
    assert_eq!(aromatic_atom_count(&positive.topology), 6);

    let pyrrole = topology(
        vec![
            AtomSpec::new(Element::N).with_explicit_hydrogens(1),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Double),
            (4, 0, BondOrder::Single),
        ],
    );
    let negative = assignment(&pyrrole, AromaticityModel::Mdl);
    assert_eq!(negative.aromatic_ring_count, 0);
    assert_eq!(aromatic_atom_count(&negative.topology), 0);

    let mixed = topology(
        vec![AtomSpec::new(Element::C); 8],
        &[
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Double),
            (5, 0, BondOrder::Single),
            (2, 6, BondOrder::Single),
            (6, 7, BondOrder::Double),
            (7, 3, BondOrder::Single),
        ],
    );
    let mixed_result = assignment(&mixed, AromaticityModel::Mdl);
    assert_eq!(mixed_result.aromatic_ring_count, 1);
    assert!(mixed_result.topology.atoms[0].is_aromatic());
    assert!(!mixed_result.topology.atoms[6].is_aromatic());
}

#[test]
fn mmff94_covers_positive_negative_hetero_hydrogen_and_input_immutability() {
    let benzene = sp2_carbon_cycle(6);
    let snapshot = benzene.clone();
    let positive = assignment(&benzene, AromaticityModel::Mmff94);
    assert_eq!(positive.aromatic_ring_count, 2);
    assert_eq!(aromatic_atom_count(&positive.topology), 6);
    assert_eq!(benzene, snapshot);

    let mut non_sp2 = benzene.clone();
    non_sp2.atoms[3].set_hybridization(Hybridization::Sp3);
    let negative = assignment(&non_sp2, AromaticityModel::Mmff94);
    assert_eq!(aromatic_atom_count(&negative.topology), 0);

    let cationic_n = topology(
        vec![
            AtomSpec::new(Element::N)
                .with_formal_charge(1)
                .with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
        ],
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Double),
            (4, 0, BondOrder::Single),
        ],
    );
    let adjusted = assignment(&cationic_n, AromaticityModel::Mmff94);
    assert_eq!(adjusted.topology.atoms[0].explicit_hydrogens(), 1);
    assert_eq!(cationic_n.atoms[0].explicit_hydrogens(), 0);
    assert_eq!(aromatic_atom_count(&adjusted.topology), 5);
}

#[test]
fn issue_1730_preexisting_aromatic_component_does_not_suppress_assignment() {
    let mut graph = topology(
        vec![AtomSpec::new(Element::C); 12],
        &[
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Double),
            (5, 0, BondOrder::Single),
            (5, 6, BondOrder::Single),
            (6, 7, BondOrder::Aromatic),
            (7, 8, BondOrder::Aromatic),
            (8, 9, BondOrder::Aromatic),
            (9, 10, BondOrder::Aromatic),
            (10, 11, BondOrder::Aromatic),
            (11, 6, BondOrder::Aromatic),
        ],
    );
    for atom in &mut graph.atoms[6..] {
        atom.set_aromatic(true);
    }
    for bond in &mut graph.bonds[7..] {
        bond.set_aromatic(true);
    }
    let snapshot = graph.clone();
    let result = assignment(&graph, AromaticityModel::Rdkit);
    assert_eq!(aromatic_atom_count(&result.topology), 12);
    assert_eq!(aromatic_bond_count(&result.topology), 12);
    assert!(!result.topology.bonds[6].is_aromatic());
    assert_eq!(graph, snapshot);
}

#[test]
fn issue_1703_zero_and_dative_exocyclic_bonds_do_not_hide_ring_aromaticity() {
    for external_order in [BondOrder::Zero, BondOrder::Dative] {
        let graph = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::N),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::N),
                AtomSpec::new(Element::FE),
            ],
            &[
                (0, 1, BondOrder::Double),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Double),
                (3, 4, BondOrder::Single),
                (4, 5, BondOrder::Double),
                (5, 0, BondOrder::Single),
                (5, 6, external_order),
            ],
        );
        let result = assignment(&graph, AromaticityModel::Rdkit);
        assert_eq!(
            aromatic_atom_count(&result.topology),
            6,
            "{external_order:?}"
        );
        assert_eq!(
            aromatic_bond_count(&result.topology),
            6,
            "{external_order:?}"
        );
        assert!(!result.topology.bonds[6].is_aromatic());
    }
}

#[test]
fn issue_1936_radical_carbocation_ring_remains_nonaromatic() {
    let mut specs = vec![AtomSpec::new(Element::C); 7];
    specs[4] = specs[4]
        .clone()
        .with_formal_charge(1)
        .with_radical_electrons(1);
    let graph = topology(
        specs,
        &[
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Single),
            (5, 6, BondOrder::Double),
            (6, 0, BondOrder::Single),
        ],
    );
    let snapshot = graph.clone();
    let result = assignment(&graph, AromaticityModel::Rdkit);
    assert_eq!(aromatic_atom_count(&result.topology), 0);
    assert_eq!(aromatic_bond_count(&result.topology), 0);
    assert_eq!(graph, snapshot);
}

#[test]
fn all_successful_models_preserve_row_identity_validate_and_are_deterministic() {
    let mut graph = sp2_carbon_cycle(6);
    graph.atoms[2].set_prop("atom-note", "kept").unwrap();
    graph.bonds[4].set_prop("bond-note", "kept").unwrap();
    let ring_info = rings(&graph);

    for model in [
        AromaticityModel::Rdkit,
        AromaticityModel::Simple,
        AromaticityModel::Mdl,
        AromaticityModel::Mmff94,
    ] {
        let params = AromaticityParams { model };
        let first = assign_aromaticity(&graph, &ring_info, &params).unwrap();
        let second = assign_aromaticity(&graph, &ring_info, &params).unwrap();
        assert_eq!(first, second, "{model:?}");
        first.topology.validate().unwrap();
        assert_eq!(first.topology.atoms.len(), graph.atoms.len());
        assert_eq!(first.topology.bonds.len(), graph.bonds.len());
        assert!(
            first
                .topology
                .atoms
                .iter()
                .enumerate()
                .all(|(row, atom)| atom.id() == AtomId::new(row))
        );
        assert!(
            first
                .topology
                .bonds
                .iter()
                .enumerate()
                .all(|(row, bond)| bond.id() == BondId::new(row))
        );
        assert_eq!(first.topology.atoms[2].prop("atom-note"), Some("kept"));
        assert_eq!(first.topology.bonds[4].prop("bond-note"), Some("kept"));
    }
}
