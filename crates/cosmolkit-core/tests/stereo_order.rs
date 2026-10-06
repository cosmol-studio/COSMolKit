use cosmolkit_core::{
    StereoOrderError, TetrahedralLigand, atom_nonzero_degree, bond_affects_atom_chirality,
    count_swaps_to_interconvert, incident_tetrahedral_bond_order, invert_tetrahedral_tag,
    remap_tetrahedral_center, tetrahedral_tag_after_order_change,
};
use cosmolkit_model::{
    Atom, AtomId, AtomMapping, AtomSpec, Bond, BondId, BondMapping, BondSpec,
    MappingValidationError, TopologyBlock, TopologyMapping, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, ChiralTag, Element};

fn atom(id: usize) -> Atom {
    Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C))
}

fn bond(id: usize, begin: usize, end: usize, order: BondOrder) -> Bond {
    Bond::from_spec(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
    )
}

fn topology(atom_count: usize, bonds: Vec<Bond>) -> TopologyBlock {
    TopologyBlock::try_from_parts((0..atom_count).map(atom).collect(), bonds, vec![], vec![])
        .unwrap()
}

fn ligand(id: usize) -> TetrahedralLigand {
    TetrahedralLigand::Bond(BondId::new(id))
}

fn mapping(
    atom_old_to_new: Vec<Option<usize>>,
    atom_new_to_old: Vec<Option<usize>>,
    bond_old_to_new: Vec<Option<usize>>,
    bond_new_to_old: Vec<Option<usize>>,
) -> TopologyMapping {
    TopologyMapping {
        atoms: AtomMapping {
            old_to_new: atom_old_to_new
                .into_iter()
                .map(|id| id.map(AtomId::new))
                .collect(),
            new_to_old: atom_new_to_old
                .into_iter()
                .map(|id| id.map(AtomId::new))
                .collect(),
        },
        bonds: BondMapping {
            old_to_new: bond_old_to_new
                .into_iter()
                .map(|id| id.map(BondId::new))
                .collect(),
            new_to_old: bond_new_to_old
                .into_iter()
                .map(|id| id.map(BondId::new))
                .collect(),
        },
    }
}

#[test]
fn perturbation_order_matches_all_three_rdkit_atom_examples() {
    let reference = [0, 1, 2, 3];
    assert_eq!(
        count_swaps_to_interconvert(&reference, &[1, 0, 2, 3]),
        Ok(1)
    );
    assert_eq!(
        count_swaps_to_interconvert(&reference, &[1, 2, 3, 0]),
        Ok(3)
    );
    assert_eq!(
        count_swaps_to_interconvert(&reference, &[1, 2, 0, 3]),
        Ok(2)
    );
}

#[test]
fn permutation_scan_is_deterministic_for_duplicates_and_reports_exact_failures() {
    assert_eq!(count_swaps_to_interconvert(&[1, 1, 2], &[2, 1, 1]), Ok(2));
    assert_eq!(
        count_swaps_to_interconvert(&[1, 2], &[1]),
        Err(StereoOrderError::PermutationLength {
            reference: 2,
            probe: 1,
        })
    );
    assert_eq!(
        count_swaps_to_interconvert(&[1, 1, 2], &[1, 2, 2]),
        Err(StereoOrderError::MissingProbeValue {
            reference_position: 1,
        })
    );
}

#[test]
fn tetrahedral_inversion_accepts_only_explicit_cw_and_ccw_tags() {
    assert_eq!(
        invert_tetrahedral_tag(ChiralTag::TetrahedralCw),
        Ok(ChiralTag::TetrahedralCcw)
    );
    assert_eq!(
        invert_tetrahedral_tag(ChiralTag::TetrahedralCcw),
        Ok(ChiralTag::TetrahedralCw)
    );
    for tag in [
        ChiralTag::Unspecified,
        ChiralTag::Tetrahedral,
        ChiralTag::Allene,
        ChiralTag::SquarePlanar,
        ChiralTag::TrigonalBipyramidal,
        ChiralTag::Octahedral,
        ChiralTag::Other,
    ] {
        assert_eq!(
            invert_tetrahedral_tag(tag),
            Err(StereoOrderError::UnsupportedTetrahedralTag { tag })
        );
    }
}

#[test]
fn bond_chirality_effect_matches_zero_unspecified_and_directional_dative_source_branches() {
    let center = AtomId::new(0);
    assert_eq!(
        bond_affects_atom_chirality(&bond(0, 0, 1, BondOrder::Unspecified), center),
        Ok(false)
    );
    assert_eq!(
        bond_affects_atom_chirality(&bond(0, 0, 1, BondOrder::Zero), center),
        Ok(false)
    );
    let dative = bond(0, 0, 1, BondOrder::Dative);
    assert_eq!(bond_affects_atom_chirality(&dative, center), Ok(false));
    assert_eq!(
        bond_affects_atom_chirality(&dative, AtomId::new(1)),
        Ok(true)
    );
    assert_eq!(
        bond_affects_atom_chirality(&bond(0, 0, 1, BondOrder::DativeOne), center),
        Ok(true)
    );
    assert_eq!(
        bond_affects_atom_chirality(&bond(0, 0, 1, BondOrder::Single), AtomId::new(2)),
        Err(StereoOrderError::BondNotIncident {
            bond: BondId::new(0),
            center: AtomId::new(2),
        })
    );
}

#[test]
fn incident_order_preserves_bond_table_order_and_nonzero_degree_filters_exactly() {
    let topology = topology(
        5,
        vec![
            bond(0, 0, 2, BondOrder::Single),
            bond(1, 1, 0, BondOrder::Dative),
            bond(2, 0, 3, BondOrder::Zero),
            bond(3, 4, 0, BondOrder::Single),
        ],
    );
    assert_eq!(
        incident_tetrahedral_bond_order(&topology, AtomId::new(0)),
        Ok(vec![ligand(0), ligand(1), ligand(3)])
    );
    assert_eq!(atom_nonzero_degree(&topology, AtomId::new(0)), Ok(3));
    assert_eq!(atom_nonzero_degree(&topology, AtomId::new(3)), Ok(0));
}

#[test]
fn topology_lookup_errors_preserve_center_and_validation_details() {
    let topology = topology(2, vec![bond(0, 0, 1, BondOrder::Single)]);
    assert_eq!(
        atom_nonzero_degree(&topology, AtomId::new(2)),
        Err(StereoOrderError::CenterOutOfRange {
            center: AtomId::new(2),
            atom_count: 2,
        })
    );

    let mut invalid = topology.clone();
    invalid.atoms[0] = invalid.atoms[0].clone().with_id(AtomId::new(1));
    assert_eq!(
        incident_tetrahedral_bond_order(&invalid, AtomId::new(0)),
        Err(StereoOrderError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            }
        ))
    );
}

#[test]
fn order_change_preserves_even_and_inverts_odd_permutations_with_implicit_ligand() {
    let reference = [ligand(0), ligand(1), ligand(2), TetrahedralLigand::Implicit];
    assert_eq!(
        tetrahedral_tag_after_order_change(ChiralTag::TetrahedralCw, &reference, &reference,),
        Ok(ChiralTag::TetrahedralCw)
    );
    assert_eq!(
        tetrahedral_tag_after_order_change(
            ChiralTag::TetrahedralCw,
            &reference,
            &[ligand(1), ligand(2), ligand(0), TetrahedralLigand::Implicit],
        ),
        Ok(ChiralTag::TetrahedralCw)
    );
    assert_eq!(
        tetrahedral_tag_after_order_change(
            ChiralTag::TetrahedralCw,
            &reference,
            &[ligand(1), ligand(0), ligand(2), TetrahedralLigand::Implicit],
        ),
        Ok(ChiralTag::TetrahedralCcw)
    );
}

#[test]
fn ligand_order_validation_rejects_capacity_duplicates_and_multiple_implicit_slots() {
    assert_eq!(
        tetrahedral_tag_after_order_change(
            ChiralTag::TetrahedralCw,
            &[ligand(0), ligand(1), ligand(2), ligand(3), ligand(4)],
            &[],
        ),
        Err(StereoOrderError::TooManyLigands {
            actual: 5,
            maximum: 4,
        })
    );
    assert_eq!(
        tetrahedral_tag_after_order_change(
            ChiralTag::TetrahedralCw,
            &[ligand(0), ligand(0)],
            &[ligand(0), ligand(0)],
        ),
        Err(StereoOrderError::DuplicateBondLigand {
            position: 1,
            bond: BondId::new(0),
        })
    );
    assert_eq!(
        tetrahedral_tag_after_order_change(
            ChiralTag::TetrahedralCw,
            &[TetrahedralLigand::Implicit, TetrahedralLigand::Implicit],
            &[TetrahedralLigand::Implicit, TetrahedralLigand::Implicit],
        ),
        Err(StereoOrderError::MultipleImplicitLigands {
            first: 0,
            second: 1,
        })
    );
}

#[test]
fn identity_remap_preserves_center_order_tag_and_all_inputs() {
    let mapping = TopologyMapping::identity(2, 3);
    let mapping_snapshot = mapping.clone();
    let source = vec![ligand(0), ligand(1), ligand(2)];
    let source_snapshot = source.clone();
    let target = source.clone();
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(1),
            &source,
            ChiralTag::TetrahedralCcw,
            &mapping,
            2,
            2,
            3,
            3,
            None,
            &target,
        )
        .unwrap(),
        cosmolkit_core::TetrahedralRemap {
            center: AtomId::new(1),
            ligands: target,
            tag: ChiralTag::TetrahedralCcw,
        }
    );
    assert_eq!(source, source_snapshot);
    assert_eq!(mapping, mapping_snapshot);
}

#[test]
fn reorder_remap_uses_both_atom_and_bond_directions_and_adjusts_parity() {
    let mapping = mapping(
        vec![Some(1), Some(0)],
        vec![Some(1), Some(0)],
        vec![Some(2), Some(0), Some(1)],
        vec![Some(1), Some(2), Some(0)],
    );
    let result = remap_tetrahedral_center(
        AtomId::new(1),
        &[ligand(0), ligand(1), ligand(2)],
        ChiralTag::TetrahedralCw,
        &mapping,
        2,
        2,
        3,
        3,
        None,
        &[ligand(0), ligand(2), ligand(1)],
    )
    .unwrap();
    assert_eq!(result.center, AtomId::new(0));
    assert_eq!(result.ligands, vec![ligand(0), ligand(2), ligand(1)]);
    assert_eq!(result.tag, ChiralTag::TetrahedralCcw);
}

#[test]
fn removed_explicit_ligand_can_be_named_as_the_single_implicit_replacement() {
    let mapping = mapping(
        vec![Some(0), Some(1)],
        vec![Some(0), Some(1)],
        vec![None, Some(0), Some(1), Some(2)],
        vec![Some(1), Some(2), Some(3)],
    );
    let source = [ligand(1), ligand(2), ligand(3), ligand(0)];
    let target = [ligand(0), ligand(1), ligand(2), TetrahedralLigand::Implicit];
    let result = remap_tetrahedral_center(
        AtomId::new(0),
        &source,
        ChiralTag::TetrahedralCw,
        &mapping,
        2,
        2,
        4,
        3,
        Some(BondId::new(0)),
        &target,
    )
    .unwrap();
    assert_eq!(result.ligands, target);
    assert_eq!(result.tag, ChiralTag::TetrahedralCw);
}

#[test]
fn removed_ligand_without_explicit_policy_fails_closed() {
    let mapping = mapping(
        vec![Some(0)],
        vec![Some(0)],
        vec![None, Some(0)],
        vec![Some(1)],
    );
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(0),
            &[ligand(0), ligand(1)],
            ChiralTag::TetrahedralCw,
            &mapping,
            1,
            1,
            2,
            1,
            None,
            &[TetrahedralLigand::Implicit, ligand(0)],
        ),
        Err(StereoOrderError::RemovedLigand {
            bond: BondId::new(0),
        })
    );
}

#[test]
fn implicit_replacement_policy_distinguishes_range_presence_and_retention_errors() {
    let identity = TopologyMapping::identity(1, 2);
    let source = [ligand(0)];
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(0),
            &source,
            ChiralTag::TetrahedralCw,
            &identity,
            1,
            1,
            2,
            2,
            Some(BondId::new(2)),
            &source,
        ),
        Err(StereoOrderError::ImplicitReplacementOutOfRange {
            bond: BondId::new(2),
            old_bond_count: 2,
        })
    );
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(0),
            &source,
            ChiralTag::TetrahedralCw,
            &identity,
            1,
            1,
            2,
            2,
            Some(BondId::new(1)),
            &source,
        ),
        Err(StereoOrderError::ImplicitReplacementNotInSource {
            bond: BondId::new(1),
        })
    );
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(0),
            &source,
            ChiralTag::TetrahedralCw,
            &identity,
            1,
            1,
            2,
            2,
            Some(BondId::new(0)),
            &source,
        ),
        Err(StereoOrderError::ImplicitReplacementRetained {
            bond: BondId::new(0),
            mapped: BondId::new(0),
        })
    );
}

#[test]
fn remap_validates_mapping_lengths_before_any_indexing() {
    let invalid = mapping(vec![], vec![Some(0)], vec![], vec![]);
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(99),
            &[],
            ChiralTag::TetrahedralCw,
            &invalid,
            1,
            1,
            0,
            0,
            None,
            &[],
        ),
        Err(StereoOrderError::InvalidMapping(
            MappingValidationError::Length {
                entity: "atom",
                direction: "old-to-new",
                actual: 0,
                expected: 1,
            }
        ))
    );
}

#[test]
fn remap_preserves_mapping_out_of_range_and_inverse_mismatch_categories() {
    let out_of_range = mapping(vec![Some(1)], vec![Some(0)], vec![], vec![]);
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(0),
            &[],
            ChiralTag::TetrahedralCw,
            &out_of_range,
            1,
            1,
            0,
            0,
            None,
            &[],
        ),
        Err(StereoOrderError::InvalidMapping(
            MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "old-to-new",
                row: 0,
                mapped: 1,
                target_count: 1,
            }
        ))
    );

    let inverse = mapping(
        vec![Some(0), Some(1)],
        vec![Some(1), Some(0)],
        vec![],
        vec![],
    );
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(0),
            &[],
            ChiralTag::TetrahedralCw,
            &inverse,
            2,
            2,
            0,
            0,
            None,
            &[],
        ),
        Err(StereoOrderError::InvalidMapping(
            MappingValidationError::InverseMismatch {
                entity: "atom",
                direction: "old-to-new",
                row: 0,
                mapped: 0,
            }
        ))
    );
}

#[test]
fn remap_distinguishes_out_of_range_removed_center_and_ligand_rows() {
    let identity = TopologyMapping::identity(1, 1);
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(1),
            &[],
            ChiralTag::TetrahedralCw,
            &identity,
            1,
            1,
            1,
            1,
            None,
            &[],
        ),
        Err(StereoOrderError::CenterOutOfRange {
            center: AtomId::new(1),
            atom_count: 1,
        })
    );
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(0),
            &[ligand(1)],
            ChiralTag::TetrahedralCw,
            &identity,
            1,
            1,
            1,
            1,
            None,
            &[],
        ),
        Err(StereoOrderError::BondOutOfRange {
            bond: BondId::new(1),
            bond_count: 1,
        })
    );

    let removed_center = mapping(vec![None, Some(0)], vec![Some(1)], vec![], vec![]);
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(0),
            &[],
            ChiralTag::TetrahedralCw,
            &removed_center,
            2,
            1,
            0,
            0,
            None,
            &[],
        ),
        Err(StereoOrderError::RemovedCenter {
            center: AtomId::new(0),
        })
    );
}

#[test]
fn append_mapping_does_not_fabricate_a_source_identity_for_new_bonds() {
    let mapping = TopologyMapping::with_appended(1, 1, 1, 1);
    assert_eq!(
        remap_tetrahedral_center(
            AtomId::new(0),
            &[ligand(0)],
            ChiralTag::TetrahedralCw,
            &mapping,
            1,
            2,
            1,
            2,
            None,
            &[ligand(1)],
        ),
        Err(StereoOrderError::MissingProbeValue {
            reference_position: 0,
        })
    );
}
