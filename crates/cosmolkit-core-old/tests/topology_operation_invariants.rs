//! Topology-changing operation invariant contract tests.
//!
//! These tests intentionally describe the target internal testing API, not the
//! smaller public surface that happens to exist today. Every topology-changing
//! operation is expected to expose a `_with_report()` test/internal entry point
//! returning the transformed molecule plus a complete topology edit report.
//!
//! The point of this suite is to make operation correctness explicit:
//! 1. The input molecule remains unchanged.
//! 2. The output molecule is globally self-consistent.
//! 3. The output follows the operation's atom/bond/coordinate/property/stereo
//!    remapping contract.

#![cfg(feature = "topology-contract-tests")]

use cosmolkit_core::{
    Atom, Bond, BondOrder, Molecule,
    testing::topology::{
        AtomMapping, BondMapping, CachePolicy, CoordinatePolicy, CowPolicy, InvariantProfile,
        KnownFailureSet, OperationSpec, OperationStatus, PropPolicy, StereoPolicy, TopologyClass,
        TopologyEditOutcome, assert_adjacency_equals_fresh_recompute,
        assert_atom_mapping_is_total_for_survivors, assert_atom_props_follow_mapping,
        assert_atom_state_invariants, assert_bond_endpoint_mapping,
        assert_bond_props_follow_mapping, assert_bond_state_invariants, assert_cache_consistency,
        assert_conformer_invariants, assert_coords_follow_mapping, assert_expected_failure_shape,
        assert_graph_invariants, assert_io_roundtrip_invariants, assert_molecule_invariants,
        assert_original_unchanged, assert_props_preserved_or_remapped,
        assert_rdkit_operation_parity, assert_stereo_indices_valid, assert_stereo_policy_followed,
        assert_topology_counts_follow_spec, assert_xfail_record_is_still_relevant,
        decorate_with_test_conformers, decorate_with_test_properties, force_all_caches,
        fresh_topology_signature, snapshot_molecule, topology_signature,
    },
};

const TOPOLOGY_CORPUS: &str = include_str!("../../../tests/corpus/topology/core.csv");
const COW_CORPUS: &str = include_str!("../../../tests/corpus/topology/cow_small.csv");
const KNOWN_FAILURES: &str =
    include_str!("../../../tests/known_failures/topology_invariants.jsonl");

const COORD_EPSILON: f64 = 1.0e-9;

#[derive(Clone, Debug)]
struct TopologyCase {
    id: String,
    smiles: String,
    tags: Vec<String>,
}

impl TopologyCase {
    fn has_tag(&self, tag: &str) -> bool {
        self.tags.iter().any(|item| item == tag)
    }
}

fn parse_cases(csv: &str) -> Vec<TopologyCase> {
    csv.lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let mut cols = line.splitn(3, ',');
            let id = cols.next().expect("case_id column").to_owned();
            let smiles = cols.next().expect("smiles column").to_owned();
            let tags = cols
                .next()
                .unwrap_or_default()
                .split_whitespace()
                .map(str::to_owned)
                .collect();
            TopologyCase { id, smiles, tags }
        })
        .collect()
}

fn parse_molecule(case: &TopologyCase) -> Molecule {
    Molecule::from_smiles(&case.smiles).unwrap_or_else(|err| {
        panic!(
            "failed to parse topology corpus case {} ({}): {err}",
            case.id, case.smiles
        )
    })
}

fn prepared_molecule(case: &TopologyCase) -> Molecule {
    let mut mol = parse_molecule(case);
    decorate_with_test_properties(&mut mol, &case.id);
    decorate_with_test_conformers(&mut mol, &case.id);
    force_all_caches(&mut mol);
    mol
}

fn run_topology_suite<F>(
    spec: OperationSpec,
    cases: &[TopologyCase],
    profile: InvariantProfile,
    op: F,
) where
    F: Fn(&Molecule, &TopologyCase) -> TopologyEditOutcome,
{
    let known_failures = KnownFailureSet::from_jsonl(KNOWN_FAILURES)
        .expect("topology known-failure records should parse");

    for case in cases {
        if !spec.applies_to_tags(&case.tags) {
            continue;
        }

        let before = prepared_molecule(case);
        let before_snapshot = snapshot_molecule(&before);
        let before_fresh_signature = fresh_topology_signature(&before);

        let outcome = op(&before, case);

        match outcome.status {
            OperationStatus::Success {
                molecule: after,
                report,
            } => {
                if let Some(xfail) = known_failures.find(&case.id, spec.name) {
                    assert_xfail_record_is_still_relevant(xfail, &before, &after, &report, &spec);
                }

                assert_original_unchanged(&before, &before_snapshot, &spec);
                assert_molecule_invariants(&after, profile.core());
                assert_topology_operation_contract(
                    case,
                    &spec,
                    &before,
                    &after,
                    &report.atom_mapping,
                    &report.bond_mapping,
                    profile,
                );

                if profile.check_io_roundtrip {
                    assert_io_roundtrip_invariants(&after, &spec, case);
                }
                if profile.check_rdkit_parity {
                    assert_rdkit_operation_parity(&before, &after, &report, &spec, case);
                }

                assert_eq!(
                    topology_signature(&before),
                    before_fresh_signature,
                    "operation {} on case {} left the input molecule with stale topology-derived state",
                    spec.name,
                    case.id
                );
            }
            OperationStatus::Error { error } => {
                let xfail = known_failures
                    .find(&case.id, spec.name)
                    .unwrap_or_else(|| {
                        panic!(
                            "unexpected topology operation error: operation={}, case={} ({}), error={error:?}",
                            spec.name, case.id, case.smiles
                        )
                    });
                assert_expected_failure_shape(xfail, &error);
                assert_original_unchanged(&before, &before_snapshot, &spec);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn assert_topology_operation_contract(
    case: &TopologyCase,
    spec: &OperationSpec,
    before: &Molecule,
    after: &Molecule,
    atom_mapping: &AtomMapping,
    bond_mapping: &BondMapping,
    profile: InvariantProfile,
) {
    assert_topology_counts_follow_spec(before, after, atom_mapping, bond_mapping, spec, case);

    assert_atom_mapping_is_total_for_survivors(before, after, atom_mapping, spec, case);
    assert_bond_endpoint_mapping(before, after, atom_mapping, bond_mapping, spec, case);

    assert_graph_invariants(after);
    assert_atom_state_invariants(after);
    assert_bond_state_invariants(after);
    assert_conformer_invariants(after);
    assert_stereo_indices_valid(after);
    assert_props_preserved_or_remapped(before, after, atom_mapping, bond_mapping, spec, case);

    match spec.coordinate_policy {
        CoordinatePolicy::PreserveRows
        | CoordinatePolicy::RemapRows
        | CoordinatePolicy::ExtendRows => {
            assert_coords_follow_mapping(before, after, atom_mapping, COORD_EPSILON, spec, case);
        }
        CoordinatePolicy::InvalidateAll => {
            assert!(
                after.coords_2d().is_none() && after.conformers_3d().is_empty(),
                "operation {} on case {} must invalidate conformers",
                spec.name,
                case.id
            );
        }
        CoordinatePolicy::ExplicitlyDocumentedDrop => {}
    }

    match spec.prop_policy {
        PropPolicy::PreserveAll | PropPolicy::RemapAtomAndBondProps => {
            assert_atom_props_follow_mapping(before, after, atom_mapping, spec, case);
            assert_bond_props_follow_mapping(before, after, bond_mapping, spec, case);
        }
        PropPolicy::DropMappedPropsWithReason => {}
    }

    assert_stereo_policy_followed(before, after, atom_mapping, bond_mapping, spec, case);

    match spec.cache_policy {
        CachePolicy::MustInvalidateOrRefresh => {
            assert_adjacency_equals_fresh_recompute(after);
            assert_cache_consistency(after, profile.cache());
        }
    }
}

fn full_operation_specs() -> Vec<OperationSpec> {
    vec![
        OperationSpec {
            name: "with_hydrogens",
            topology_class: TopologyClass::Strong,
            coordinate_policy: CoordinatePolicy::ExtendRows,
            prop_policy: PropPolicy::PreserveAll,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["hydrogens"],
            rdkit_parity: true,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "without_hydrogens",
            topology_class: TopologyClass::Strong,
            coordinate_policy: CoordinatePolicy::RemapRows,
            prop_policy: PropPolicy::RemapAtomAndBondProps,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["hydrogens", "explicit_h"],
            rdkit_parity: true,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "kekulize",
            topology_class: TopologyClass::Weak,
            coordinate_policy: CoordinatePolicy::PreserveRows,
            prop_policy: PropPolicy::PreserveAll,
            stereo_policy: StereoPolicy::PreserveValidReferences,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["aromatic", "kekulize"],
            rdkit_parity: true,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "sanitize",
            topology_class: TopologyClass::Weak,
            coordinate_policy: CoordinatePolicy::PreserveRows,
            prop_policy: PropPolicy::PreserveAll,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["sanitize"],
            rdkit_parity: true,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "add_atom",
            topology_class: TopologyClass::Strong,
            coordinate_policy: CoordinatePolicy::ExtendRows,
            prop_policy: PropPolicy::PreserveAll,
            stereo_policy: StereoPolicy::PreserveValidReferences,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["basic", "edit_atom"],
            rdkit_parity: false,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "remove_atom",
            topology_class: TopologyClass::Strong,
            coordinate_policy: CoordinatePolicy::RemapRows,
            prop_policy: PropPolicy::RemapAtomAndBondProps,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["remove_atom", "basic", "stereo"],
            rdkit_parity: false,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "replace_atom",
            topology_class: TopologyClass::Weak,
            coordinate_policy: CoordinatePolicy::PreserveRows,
            prop_policy: PropPolicy::RemapAtomAndBondProps,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["edit_atom", "basic"],
            rdkit_parity: false,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "add_bond",
            topology_class: TopologyClass::Strong,
            coordinate_policy: CoordinatePolicy::PreserveRows,
            prop_policy: PropPolicy::PreserveAll,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["basic", "disconnected"],
            rdkit_parity: false,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "remove_bond",
            topology_class: TopologyClass::Strong,
            coordinate_policy: CoordinatePolicy::PreserveRows,
            prop_policy: PropPolicy::RemapAtomAndBondProps,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["basic", "ring", "double_bond"],
            rdkit_parity: false,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "set_bond_type",
            topology_class: TopologyClass::Weak,
            coordinate_policy: CoordinatePolicy::PreserveRows,
            prop_policy: PropPolicy::PreserveAll,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["double_bond", "basic"],
            rdkit_parity: false,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "set_aromatic",
            topology_class: TopologyClass::Weak,
            coordinate_policy: CoordinatePolicy::PreserveRows,
            prop_policy: PropPolicy::PreserveAll,
            stereo_policy: StereoPolicy::PreserveValidReferences,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["aromatic"],
            rdkit_parity: false,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "set_formal_charge",
            topology_class: TopologyClass::Weak,
            coordinate_policy: CoordinatePolicy::PreserveRows,
            prop_policy: PropPolicy::PreserveAll,
            stereo_policy: StereoPolicy::PreserveValidReferences,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["charged", "sanitize"],
            rdkit_parity: false,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "update_implicit_hydrogens",
            topology_class: TopologyClass::Weak,
            coordinate_policy: CoordinatePolicy::PreserveRows,
            prop_policy: PropPolicy::PreserveAll,
            stereo_policy: StereoPolicy::PreserveValidReferences,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["hydrogens", "charged", "sanitize"],
            rdkit_parity: true,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "update_stereo",
            topology_class: TopologyClass::Weak,
            coordinate_policy: CoordinatePolicy::PreserveRows,
            prop_policy: PropPolicy::PreserveAll,
            stereo_policy: StereoPolicy::RecomputeFromCurrentTopology,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["stereo"],
            rdkit_parity: true,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "renumber_atoms",
            topology_class: TopologyClass::Strong,
            coordinate_policy: CoordinatePolicy::RemapRows,
            prop_policy: PropPolicy::RemapAtomAndBondProps,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["renumber", "stereo", "ring"],
            rdkit_parity: true,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "fragment",
            topology_class: TopologyClass::Strong,
            coordinate_policy: CoordinatePolicy::RemapRows,
            prop_policy: PropPolicy::RemapAtomAndBondProps,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["fragment", "disconnected", "ring"],
            rdkit_parity: false,
            io_roundtrip: true,
        },
        OperationSpec {
            name: "combine",
            topology_class: TopologyClass::Strong,
            coordinate_policy: CoordinatePolicy::ExtendRows,
            prop_policy: PropPolicy::RemapAtomAndBondProps,
            stereo_policy: StereoPolicy::RemapOrDropInvalid,
            cache_policy: CachePolicy::MustInvalidateOrRefresh,
            cow_policy: CowPolicy::MustNotMutateOriginal,
            required_tags: &["combine", "disconnected"],
            rdkit_parity: false,
            io_roundtrip: true,
        },
    ]
}

#[test]
fn topology_operations_preserve_global_molecule_invariants() {
    let cases = parse_cases(TOPOLOGY_CORPUS);

    for spec in full_operation_specs() {
        let profile = if spec.rdkit_parity {
            InvariantProfile::core_with_rdkit()
        } else if spec.io_roundtrip {
            InvariantProfile::core_with_io()
        } else {
            InvariantProfile::core()
        };

        run_topology_suite(spec.clone(), &cases, profile, |mol, case| {
            apply_operation(spec.name, mol, case)
        });
    }
}

#[test]
fn cow_isolation_uses_small_topology_subset_not_the_full_parity_corpus() {
    let cases = parse_cases(COW_CORPUS);
    let profile = InvariantProfile::cow_only();

    for spec in full_operation_specs()
        .into_iter()
        .filter(|spec| spec.cow_policy == CowPolicy::MustNotMutateOriginal)
    {
        run_topology_suite(spec.clone(), &cases, profile, |mol, case| {
            apply_operation(spec.name, mol, case)
        });
    }
}

#[test]
#[should_panic(expected = "invalid bond endpoint")]
fn invariant_helpers_reject_invalid_bond_endpoints() {
    let cases = parse_cases(TOPOLOGY_CORPUS);
    let base = prepared_molecule(&cases[1]);

    let mut invalid_bond_endpoint = base.clone();
    let invalid_index = invalid_bond_endpoint.atoms().len();
    invalid_bond_endpoint
        .bonds_mut()
        .first_mut()
        .expect("case should have a bond")
        .end_atom = invalid_index;
    assert_graph_invariants(&invalid_bond_endpoint);
}

#[test]
#[should_panic(expected = "stale adjacency")]
fn invariant_helpers_reject_stale_adjacency() {
    let cases = parse_cases(TOPOLOGY_CORPUS);
    let mut stale_adjacency = prepared_molecule(&cases[1]);
    force_all_caches(&mut stale_adjacency);
    let replacement_endpoint = stale_adjacency.atoms().len() - 1;
    stale_adjacency
        .bonds_mut()
        .first_mut()
        .expect("case should have a bond")
        .end_atom = replacement_endpoint;
    assert_adjacency_equals_fresh_recompute(&stale_adjacency);
}

#[test]
#[should_panic(expected = "coordinate row count")]
fn invariant_helpers_reject_misaligned_conformer_rows() {
    let cases = parse_cases(TOPOLOGY_CORPUS);
    let base = prepared_molecule(&cases[1]);
    let mut invalid_coords = base.clone();
    invalid_coords
        .conformers_3d_mut()
        .first_mut()
        .expect("test conformer should exist")
        .pop();
    assert_conformer_invariants(&invalid_coords);
}

#[test]
#[should_panic(expected = "stereo")]
fn invariant_helpers_reject_invalid_stereo_references() {
    let cases = parse_cases(TOPOLOGY_CORPUS);
    let base = prepared_molecule(&cases[1]);
    let mut invalid_stereo = base.clone();
    let invalid_index = invalid_stereo.atoms().len();
    invalid_stereo
        .bonds_mut()
        .first_mut()
        .expect("case should have a bond")
        .stereo_atoms
        .push(invalid_index);
    assert_stereo_indices_valid(&invalid_stereo);
}

fn apply_operation(
    operation: &'static str,
    mol: &Molecule,
    case: &TopologyCase,
) -> TopologyEditOutcome {
    match operation {
        "with_hydrogens" => mol.with_hydrogens_with_report(),
        "without_hydrogens" => mol.without_hydrogens_with_report(),
        "kekulize" => mol.with_kekulized_bonds_with_report(false),
        "sanitize" => mol.sanitize_with_report(),
        "add_atom" => mol.add_atom_with_report(test_atom_for_case(case)),
        "remove_atom" => mol.remove_atom_with_report(target_atom_for_case(mol, case)),
        "replace_atom" => {
            mol.replace_atom_with_report(target_atom_for_case(mol, case), test_atom_for_case(case))
        }
        "add_bond" => mol.add_bond_with_report(test_bond_for_case(mol, case)),
        "remove_bond" => mol.remove_bond_with_report(target_bond_for_case(mol, case)),
        "set_bond_type" => {
            mol.set_bond_type_with_report(target_bond_for_case(mol, case), BondOrder::Single)
        }
        "set_aromatic" => mol.set_aromatic_with_report(target_atom_for_case(mol, case), true),
        "set_formal_charge" => {
            mol.set_formal_charge_with_report(target_atom_for_case(mol, case), 1)
        }
        "update_implicit_hydrogens" => mol.update_implicit_hydrogens_with_report(),
        "update_stereo" => mol.update_stereo_with_report(),
        "renumber_atoms" => mol.renumber_atoms_with_report(reverse_order(mol)),
        "fragment" => mol.fragment_with_report(fragment_atoms_for_case(mol, case)),
        "combine" => mol.combine_with_report(&combine_partner_for_case(case)),
        other => panic!("unhandled topology operation spec {other}"),
    }
}

fn test_atom_for_case(case: &TopologyCase) -> Atom {
    Atom::test_builder()
        .atomic_num(if case.has_tag("charged") { 7 } else { 6 })
        .formal_charge(if case.has_tag("charged") { 1 } else { 0 })
        .prop("topology_test_atom", &case.id)
        .build()
}

fn test_bond_for_case(mol: &Molecule, case: &TopologyCase) -> Bond {
    let (begin_atom, end_atom) = if case.has_tag("disconnected") {
        first_two_fragment_bridge_atoms(mol)
    } else {
        (0, mol.atoms().len() - 1)
    };

    Bond::test_builder()
        .begin_atom(begin_atom)
        .end_atom(end_atom)
        .order(BondOrder::Single)
        .prop("topology_test_bond", &case.id)
        .build()
}

fn target_atom_for_case(mol: &Molecule, case: &TopologyCase) -> usize {
    if case.has_tag("explicit_h") {
        mol.atoms()
            .iter()
            .find(|atom| atom.atomic_num == 1)
            .map(|atom| atom.index)
            .unwrap_or(0)
    } else if case.has_tag("stereo") {
        mol.atoms()
            .iter()
            .find(|atom| !matches!(atom.chiral_tag, cosmolkit_core::ChiralTag::Unspecified))
            .map(|atom| atom.index)
            .unwrap_or(0)
    } else {
        mol.atoms().len() / 2
    }
}

fn target_bond_for_case(mol: &Molecule, case: &TopologyCase) -> usize {
    if case.has_tag("double_bond") {
        mol.bonds()
            .iter()
            .find(|bond| bond.order == BondOrder::Double)
            .map(|bond| bond.index)
            .unwrap_or(0)
    } else {
        mol.bonds().len() / 2
    }
}

fn reverse_order(mol: &Molecule) -> Vec<usize> {
    (0..mol.atoms().len()).rev().collect()
}

fn fragment_atoms_for_case(mol: &Molecule, case: &TopologyCase) -> Vec<usize> {
    if case.has_tag("disconnected") {
        mol.connected_components()
            .first()
            .expect("disconnected case should have at least one component")
            .clone()
    } else {
        (0..(mol.atoms().len() + 1) / 2).collect()
    }
}

fn combine_partner_for_case(case: &TopologyCase) -> Molecule {
    let partner_smiles = if case.has_tag("charged") {
        "[NH4+]"
    } else {
        "O"
    };
    let mut partner =
        Molecule::from_smiles(partner_smiles).expect("topology combine partner should parse");
    decorate_with_test_properties(&mut partner, &format!("{}_partner", case.id));
    decorate_with_test_conformers(&mut partner, &format!("{}_partner", case.id));
    partner
}

fn first_two_fragment_bridge_atoms(mol: &Molecule) -> (usize, usize) {
    let components = mol.connected_components();
    assert!(
        components.len() >= 2,
        "disconnected add_bond test requires at least two fragments"
    );
    (components[0][0], components[1][0])
}
