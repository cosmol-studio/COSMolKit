use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_core::{
    Molecule, SmilesWriteParams, ValenceModel, assign_valence, cached_valence_assignment,
    io::{
        molfile::read_mol_record_from_str_with_params,
        sdf::{SdfCoordinateMode, SdfReadParams},
    },
};
use serde::Deserialize;

#[derive(Debug, Deserialize, PartialEq, Eq)]
struct AtomState(u8, i8, String, Option<u16>, bool, u8, usize, i32, i32, usize, i32, u8);

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
struct BondState(usize, usize, String, String, String, Vec<usize>, bool);

#[derive(Debug, Deserialize)]
struct StereoDirectionRecord {
    case_id: String,
    source_row: usize,
    dimension: String,
    molblock: String,
    canonical_smiles: String,
    atoms: Vec<AtomState>,
    bonds: Vec<BondState>,
    positions: Vec<[f64; 3]>,
}

fn fixture_records() -> Vec<StereoDirectionRecord> {
    include_str!("../../../testdata/molblock/fixtures/molblock_stereo_direction_components.jsonl")
        .lines()
        .enumerate()
        .map(|(line_idx, line)| {
            serde_json::from_str(line).unwrap_or_else(|error| {
                panic!(
                    "invalid MolBlock stereo-direction fixture line {}: {error}",
                    line_idx + 1
                )
            })
        })
        .collect()
}

fn atom_states(molecule: &Molecule) -> Vec<AtomState> {
    let computed_valence;
    let valence = if let Some(cached) = cached_valence_assignment(molecule) {
        cached
    } else {
        computed_valence =
            assign_valence(molecule, ValenceModel::RdkitLike).expect("stress fixture should have computable valence");
        &computed_valence
    };
    let mut degrees = vec![0usize; molecule.num_atoms()];
    for bond in molecule.bonds() {
        degrees[bond.begin().index()] += 1;
        degrees[bond.end().index()] += 1;
    }

    molecule
        .atoms()
        .iter()
        .map(|atom| {
            let idx = atom.id().index();
            let explicit_valence = valence.explicit_valence[idx];
            let implicit_hydrogens = valence.implicit_hydrogens[idx];
            AtomState(
                atom.atomic_number(),
                atom.formal_charge(),
                atom.chiral_tag().rdkit_name().to_owned(),
                atom.isotope(),
                atom.is_aromatic(),
                atom.explicit_hydrogens(),
                degrees[idx],
                explicit_valence,
                implicit_hydrogens,
                usize::from(atom.explicit_hydrogens()) + usize::try_from(implicit_hydrogens).unwrap(),
                explicit_valence + implicit_hydrogens,
                atom.radical_electrons(),
            )
        })
        .collect()
}

fn bond_states(molecule: &Molecule) -> Vec<BondState> {
    molecule
        .bonds()
        .iter()
        .map(|bond| {
            BondState(
                bond.begin().index(),
                bond.end().index(),
                bond.order().rdkit_name().to_owned(),
                bond.direction().rdkit_name().to_owned(),
                bond.stereo().rdkit_name().to_owned(),
                bond.stereo_atoms()
                    .map(|atoms| atoms.into_iter().map(|atom| atom.index()).collect())
                    .unwrap_or_default(),
                bond.is_aromatic(),
            )
        })
        .collect()
}

fn find(parent: &mut [usize], mut index: usize) -> usize {
    while parent[index] != index {
        parent[index] = parent[parent[index]];
        index = parent[index];
    }
    index
}

fn union(parent: &mut [usize], left: usize, right: usize) {
    let left_root = find(parent, left);
    let right_root = find(parent, right);
    if left_root != right_root {
        parent[right_root] = left_root;
    }
}

fn stereo_direction_states_equivalent(expected: &[BondState], actual: &[BondState]) -> bool {
    if expected.len() != actual.len() {
        return false;
    }

    let is_directional = |direction: &str| matches!(direction, "ENDUPRIGHT" | "ENDDOWNRIGHT");
    for (expected_bond, actual_bond) in expected.iter().zip(actual) {
        if expected_bond.0 != actual_bond.0
            || expected_bond.1 != actual_bond.1
            || expected_bond.2 != actual_bond.2
            || expected_bond.4 != actual_bond.4
            || expected_bond.5 != actual_bond.5
            || expected_bond.6 != actual_bond.6
        {
            return false;
        }
        if expected_bond.3 != actual_bond.3 && !(is_directional(&expected_bond.3) && is_directional(&actual_bond.3)) {
            return false;
        }
    }

    let direction_bonds = expected
        .iter()
        .enumerate()
        .filter_map(|(index, bond)| is_directional(&bond.3).then_some(index))
        .collect::<Vec<_>>();
    let mut parent = (0..expected.len()).collect::<Vec<_>>();
    let mut constrained = BTreeSet::new();

    for bond in expected {
        if bond.2 != "DOUBLE" || !matches!(bond.4.as_str(), "E" | "Z" | "CIS" | "TRANS") {
            continue;
        }
        let adjacent = direction_bonds
            .iter()
            .copied()
            .filter(|index| {
                let directional = &expected[*index];
                directional.0 == bond.0 || directional.0 == bond.1 || directional.1 == bond.0 || directional.1 == bond.1
            })
            .collect::<Vec<_>>();
        constrained.extend(adjacent.iter().copied());
        if let Some((&first, rest)) = adjacent.split_first() {
            for &other in rest {
                union(&mut parent, first, other);
            }
        }
    }

    let mut component_inversion = BTreeMap::new();
    for index in direction_bonds {
        let inverted = expected[index].3 != actual[index].3;
        if inverted && !constrained.contains(&index) {
            return false;
        }
        let root = find(&mut parent, index);
        if component_inversion
            .insert(root, inverted)
            .is_some_and(|previous| previous != inverted)
        {
            return false;
        }
    }
    true
}

fn assert_coordinates(record: &StereoDirectionRecord, molecule: &Molecule) {
    if record.dimension == "2d" {
        let actual = molecule
            .coordinates_2d()
            .unwrap_or_else(|| panic!("{} should contain 2D coordinates", record.case_id));
        assert_eq!(actual.len(), record.positions.len());
        for (atom_idx, (actual, expected)) in actual.iter().zip(&record.positions).enumerate() {
            assert!(
                (actual[0] - expected[0]).abs() <= 1.0e-8
                    && (actual[1] - expected[1]).abs() <= 1.0e-8
                    && expected[2].abs() <= 1.0e-8,
                "{} coordinate mismatch at atom {atom_idx}",
                record.case_id
            );
        }
    } else {
        let actual = molecule
            .conformers_3d()
            .first()
            .unwrap_or_else(|| panic!("{} should contain 3D coordinates", record.case_id))
            .coordinates();
        assert_eq!(actual.len(), record.positions.len());
        for (atom_idx, (actual, expected)) in actual.iter().zip(&record.positions).enumerate() {
            assert!(
                actual
                    .iter()
                    .zip(expected)
                    .all(|(actual, expected)| (actual - expected).abs() <= 1.0e-8),
                "{} coordinate mismatch at atom {atom_idx}",
                record.case_id
            );
        }
    }
}

#[test]
fn parsed_stereo_direction_components_match_stable_rdkit_semantics() {
    let records = fixture_records();
    let mut source_rows = BTreeSet::new();

    for record in &records {
        let coordinate_mode = match record.dimension.as_str() {
            "2d" => SdfCoordinateMode::Require2D,
            "3d" => SdfCoordinateMode::Require3D,
            dimension => panic!("invalid fixture dimension {dimension:?}"),
        };
        let molecule = read_mol_record_from_str_with_params(
            &record.molblock,
            SdfReadParams {
                sanitize: true,
                remove_hs: false,
                strict_parsing: true,
                process_property_lists: false,
                coordinate_mode,
                ..Default::default()
            },
        )
        .unwrap_or_else(|error| panic!("failed to parse {}: {error}", record.case_id))
        .molecule;

        assert_eq!(
            molecule.to_smiles_with_params(&SmilesWriteParams::default()).unwrap(),
            record.canonical_smiles,
            "canonical SMILES mismatch for {}",
            record.case_id
        );
        assert_eq!(
            atom_states(&molecule),
            record.atoms,
            "atom-state mismatch for {}",
            record.case_id
        );
        assert!(
            stereo_direction_states_equivalent(&record.bonds, &bond_states(&molecule)),
            "bond-state mismatch beyond a uniform stereo-component inversion for {}",
            record.case_id
        );
        assert_coordinates(record, &molecule);
        source_rows.insert(record.source_row);
    }

    assert_eq!(records.len(), 14, "fixture must retain both dimensions");
    assert_eq!(source_rows.len(), 7, "fixture must retain every source row");
}

#[test]
fn stereo_direction_equivalence_rejects_local_and_unconstrained_inversions() {
    let mut expected = vec![
        BondState(0, 1, "SINGLE".into(), "ENDUPRIGHT".into(), "NONE".into(), vec![], false),
        BondState(1, 2, "DOUBLE".into(), "NONE".into(), "E".into(), vec![0, 3], false),
        BondState(
            2,
            3,
            "SINGLE".into(),
            "ENDDOWNRIGHT".into(),
            "NONE".into(),
            vec![],
            false,
        ),
        BondState(4, 5, "SINGLE".into(), "ENDUPRIGHT".into(), "NONE".into(), vec![], false),
    ];
    let mut uniform = expected.clone();
    uniform[0].3 = "ENDDOWNRIGHT".into();
    uniform[2].3 = "ENDUPRIGHT".into();
    assert!(stereo_direction_states_equivalent(&expected, &uniform));

    let mut local = uniform.clone();
    local[2].3 = "ENDDOWNRIGHT".into();
    assert!(!stereo_direction_states_equivalent(&expected, &local));

    expected[3].3 = "ENDDOWNRIGHT".into();
    assert!(!stereo_direction_states_equivalent(&expected, &uniform));
}
