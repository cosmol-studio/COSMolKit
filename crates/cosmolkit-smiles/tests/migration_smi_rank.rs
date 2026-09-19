#[path = "../src/canonical_rank.rs"]
mod canonical_rank;

use canonical_rank::{CanonicalRankPolicy, rank_component_atoms, rank_component_atoms_with_policy};
use cosmolkit_core::{
    CanonicalRankError, CanonicalRankParams, rank_fragment_atoms, rank_mol_atoms_with_params,
};
use cosmolkit_smiles::{parse_smiles, write_smiles};

fn topology(smiles: &str) -> cosmolkit_model::TopologyBlock {
    parse_smiles(smiles, &Default::default())
        .unwrap_or_else(|error| panic!("failed to parse {smiles}: {error}"))
        .topology
}

fn untied_policy() -> CanonicalRankPolicy {
    CanonicalRankPolicy {
        break_ties: false,
        ..CanonicalRankPolicy::default()
    }
}

#[test]
fn source_rank_vectors_cover_chains_branches_rings_and_labels() {
    for (smiles, expected) in [
        ("OCC", vec![1, 2, 0]),
        ("OC(C)C", vec![2, 3, 0, 1]),
        ("c1ccncc1", vec![0, 1, 3, 5, 4, 2]),
        ("C12(CCCCC1)CCCCC2", vec![10, 6, 2, 0, 3, 7, 8, 4, 1, 5, 9]),
        ("C1C2C3C1C2C3", vec![0, 2, 4, 3, 5, 1]),
        ("[13CH3]CO", vec![0, 2, 1]),
        ("[CH3:7]CO", vec![2, 1, 0]),
    ] {
        let topology = topology(smiles);
        let component = (0..topology.atoms.len()).collect::<Vec<_>>();
        assert_eq!(
            rank_component_atoms(&topology, &component).unwrap(),
            expected,
            "{smiles}"
        );
    }
}

#[test]
fn component_projection_is_deterministic_and_rejects_bad_indices() {
    let topology = topology("CC.O");
    assert_eq!(
        rank_component_atoms(&topology, &[2]).unwrap(),
        vec![usize::MAX, usize::MAX, 0]
    );
    assert_eq!(
        rank_component_atoms(&topology, &[]).unwrap(),
        vec![usize::MAX; 3]
    );
    assert_eq!(
        rank_component_atoms(&topology, &[0, 0]),
        Err(CanonicalRankError::DuplicateComponentAtom { atom_index: 0 })
    );
    assert_eq!(
        rank_component_atoms(&topology, &[3]),
        Err(CanonicalRankError::ComponentAtomOutOfRange {
            atom_index: 3,
            atom_count: 3,
        })
    );

    let first = rank_component_atoms(&topology, &[0, 1]).unwrap();
    let second = rank_component_atoms(&topology, &[0, 1]).unwrap();
    assert_eq!(first, second);
}

#[test]
fn tie_isotope_and_atom_map_options_remain_independent() {
    let ethane = topology("CC");
    assert_eq!(
        rank_component_atoms_with_policy(&ethane, &[0, 1], untied_policy()).unwrap(),
        vec![0, 0]
    );
    assert_eq!(rank_component_atoms(&ethane, &[0, 1]).unwrap(), vec![0, 1]);

    let isotopes = topology("[13CH4].[12CH4]");
    let mut policy = untied_policy();
    policy.include_isotopes = false;
    assert_eq!(
        rank_component_atoms_with_policy(&isotopes, &[0, 1], policy).unwrap(),
        vec![0, 0]
    );
    policy.include_isotopes = true;
    assert_ne!(
        rank_component_atoms_with_policy(&isotopes, &[0, 1], policy).unwrap()[0],
        rank_component_atoms_with_policy(&isotopes, &[0, 1], policy).unwrap()[1]
    );

    let maps = topology("[CH4:7].[CH4:3]");
    let mut policy = untied_policy();
    policy.include_atom_maps = false;
    assert_eq!(
        rank_component_atoms_with_policy(&maps, &[0, 1], policy).unwrap(),
        vec![0, 0]
    );
    policy.include_atom_maps = true;
    let ranked = rank_component_atoms_with_policy(&maps, &[0, 1], policy).unwrap();
    assert_ne!(ranked[0], ranked[1]);
}

#[test]
fn every_source_option_axis_reaches_the_unique_core_engine() {
    let topology = topology("F/C=C/[C@H](Cl)Br.C1CC[C@H](C)C1 |o1:3,9|");
    for mutate in [
        (|params: &mut CanonicalRankParams| params.include_chirality = false)
            as fn(&mut CanonicalRankParams),
        |params: &mut CanonicalRankParams| params.include_chiral_presence = true,
        |params: &mut CanonicalRankParams| params.include_stereo_groups = false,
        |params: &mut CanonicalRankParams| params.use_non_stereo_ranks = true,
        |params: &mut CanonicalRankParams| params.include_ring_stereo = false,
    ] {
        let mut params = CanonicalRankParams::default();
        mutate(&mut params);
        let first = rank_mol_atoms_with_params(&topology, &params).unwrap();
        let second = rank_mol_atoms_with_params(&topology, &params).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.len(), topology.atoms.len());
    }
}

#[test]
fn mask_errors_are_typed_and_empty_topology_is_stable() {
    assert_eq!(
        rank_mol_atoms_with_params(
            &cosmolkit_model::TopologyBlock::default(),
            &CanonicalRankParams::default(),
        )
        .unwrap(),
        Vec::<usize>::new()
    );

    let topology = topology("CC");
    assert_eq!(
        rank_fragment_atoms(&topology, &[true], &[true]),
        Err(CanonicalRankError::AtomMaskLength {
            expected: 2,
            actual: 1,
        })
    );
    assert_eq!(
        rank_fragment_atoms(&topology, &[true, true], &[]),
        Err(CanonicalRankError::BondMaskLength {
            expected: 1,
            actual: 0,
        })
    );
}

#[test]
fn canonical_writer_uses_component_ranks_for_real_disconnected_inputs() {
    for (left, right, expected) in [
        ("OCC", "CCO", "CCO"),
        ("OC(C)C", "CC(C)O", "CC(C)O"),
        ("N.CCO", "OCC.N", "CCO.N"),
    ] {
        let left = parse_smiles(left, &Default::default()).unwrap();
        let right = parse_smiles(right, &Default::default()).unwrap();
        assert_eq!(write_smiles(&left).unwrap(), expected);
        assert_eq!(write_smiles(&right).unwrap(), expected);
    }
}
