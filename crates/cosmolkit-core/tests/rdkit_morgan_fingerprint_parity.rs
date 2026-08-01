use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit_core::fingerprint::{
    AdditionalOutput, FingerprintFuncArguments, MorganFingerprintGenerator,
    getMorganGeneratorWithParams,
};
use cosmolkit_core::{
    Molecule, MorganAtomInvariantsGenerator, MorganBondInvariantsGenerator,
    MorganFingerprintParams, morgan_get_fingerprint, morgan_get_fingerprint_as_bit_vect,
    morgan_get_hashed_fingerprint,
};
use serde::Deserialize;

mod common;
use common::parity_data;

#[derive(Debug, Deserialize)]
struct MorganBranch {
    ok: bool,
    on_bits: Option<Vec<usize>>,
    num_on_bits: Option<usize>,
    tanimoto_to_previous: Option<f64>,
    error: Option<String>,
    sparse_count: MorganSparseCountGolden,
    sparse_bit: MorganSparseBitGolden,
    hashed_count: MorganSparseCountGolden,
    explicit_bit: MorganExplicitBitGolden,
    additional_output: Option<MorganAdditionalOutputGolden>,
}

#[derive(Debug, Deserialize)]
struct MorganSparseCountGolden {
    ok: bool,
    nonzero_elements: Option<BTreeMap<String, i32>>,
    error: Option<String>,
    additional_output: Option<MorganAdditionalOutputGolden>,
}

#[derive(Debug, Deserialize)]
struct MorganSparseBitGolden {
    ok: bool,
    on_bits: Option<Vec<u64>>,
    num_on_bits: Option<usize>,
    error: Option<String>,
    additional_output: Option<MorganAdditionalOutputGolden>,
}

#[derive(Debug, Deserialize)]
struct MorganExplicitBitGolden {
    ok: bool,
    on_bits: Option<Vec<usize>>,
    num_on_bits: Option<usize>,
    error: Option<String>,
    additional_output: Option<MorganAdditionalOutputGolden>,
}

#[derive(Debug, Deserialize)]
struct MorganAdditionalOutputGolden {
    atom_counts: Vec<u32>,
    atom_to_bits: Vec<Vec<u64>>,
    bit_info_map: BTreeMap<String, Vec<[usize; 2]>>,
    atoms_per_bit: BTreeMap<String, Vec<Vec<usize>>>,
}

#[derive(Debug, Deserialize)]
struct MorganRecord {
    smiles: String,
    rdkit_ok: bool,
    branches: BTreeMap<String, MorganBranch>,
    error: Option<String>,
}

fn load_golden() -> Vec<MorganRecord> {
    let path = parity_data::golden_path("morgan_fingerprint.jsonl");
    assert!(
        path.exists(),
        "missing RDKit Morgan fingerprint golden: {}. Regenerate all goldens with \
         .venv/bin/python tools/testdata/rdkit/generate_all.py --python .venv/bin/python --profile {} --suite fingerprint --clean --jobs 4",
        parity_data::profile_name(),
        path.display()
    );
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!("failed to open {}: {err}", path.display());
    });
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(idx, line)| {
            let line = line.unwrap_or_else(|err| {
                panic!("failed to read {} line {}: {err}", path.display(), idx + 1)
            });
            serde_json::from_str(&line).unwrap_or_else(|err| {
                panic!("failed to parse {} line {}: {err}", path.display(), idx + 1)
            })
        })
        .collect()
}

fn branch_func_args(
    name: &str,
    mol: &Molecule,
    collect_additional_output: bool,
) -> FingerprintFuncArguments {
    let mut args = FingerprintFuncArguments {
        from_atoms: match name {
            "r2_n2048_from_atom_0" if !mol.atoms().is_empty() => Some(vec![0]),
            _ => None,
        },
        ignore_atoms: match name {
            "r2_n2048_ignore_atom_0" if !mol.atoms().is_empty() => Some(vec![0]),
            _ => None,
        },
        custom_atom_invariants: match name {
            "r2_n2048_custom_atom_invariants" => {
                Some((0..mol.atoms().len()).map(|idx| idx as u32 + 1).collect())
            }
            "r2_n2048_only_nonzero_custom_atom" => Some(
                (0..mol.atoms().len())
                    .map(|idx| if idx % 2 == 0 { 0 } else { idx as u32 + 1 })
                    .collect(),
            ),
            _ => None,
        },
        custom_bond_invariants: match name {
            "r2_n2048_custom_bond_invariants" => {
                Some((0..mol.bonds().len()).map(|idx| idx as u32 + 7).collect())
            }
            _ => None,
        },
        ..Default::default()
    };
    if collect_additional_output {
        let mut output = AdditionalOutput::new();
        output.allocate_atom_counts();
        output.allocate_atom_to_bits();
        output.allocate_bit_info_map();
        output.allocate_atoms_per_bit();
        args.additional_output = Some(output);
    }
    args
}

fn branch_params(name: &str, mol: &Molecule) -> MorganFingerprintParams {
    match name {
        "r2_n2048_chiral_false_bonds_true" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            use_chirality: false,
            use_bond_types: true,
            ..Default::default()
        },
        "r2_n2048_chiral_true_bonds_true" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            use_chirality: true,
            use_bond_types: true,
            ..Default::default()
        },
        "r3_n2048_chiral_false_bonds_true" => MorganFingerprintParams {
            radius: 3,
            n_bits: 2048,
            use_chirality: false,
            use_bond_types: true,
            ..Default::default()
        },
        "r2_n1024_chiral_false_bonds_false" => MorganFingerprintParams {
            radius: 2,
            n_bits: 1024,
            use_chirality: false,
            use_bond_types: false,
            ..Default::default()
        },
        "r2_n2048_count_sim_bounds_1_2_4_8" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            count_simulation: true,
            count_bounds: vec![1, 2, 4, 8],
            ..Default::default()
        },
        "r2_n2048_redundant_true" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            include_redundant_environments: true,
            ..Default::default()
        },
        "r2_n2048_ring_membership_false" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership: false,
            },
            ..Default::default()
        },
        "r2_n2048_from_atom_0" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            from_atoms: (!mol.atoms().is_empty()).then(|| vec![0]),
            ..Default::default()
        },
        "r2_n2048_ignore_atom_0" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            ignore_atoms: (!mol.atoms().is_empty()).then(|| vec![0]),
            ..Default::default()
        },
        "r2_n2048_custom_atom_invariants" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            custom_atom_invariants: Some(
                (0..mol.atoms().len()).map(|idx| idx as u32 + 1).collect(),
            ),
            ..Default::default()
        },
        "r2_n2048_custom_bond_invariants" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            custom_bond_invariants: Some(
                (0..mol.bonds().len()).map(|idx| idx as u32 + 7).collect(),
            ),
            ..Default::default()
        },
        "r2_n2048_only_nonzero_custom_atom" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            only_nonzero_invariants: true,
            custom_atom_invariants: Some(
                (0..mol.atoms().len())
                    .map(|idx| if idx % 2 == 0 { 0 } else { idx as u32 + 1 })
                    .collect(),
            ),
            ..Default::default()
        },
        "r2_n2048_atom_inv_gen_ring_false" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership: false,
            },
            ..Default::default()
        },
        "r2_n2048_feature_atom_inv_gen" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Feature,
            ..Default::default()
        },
        "r2_n2048_bond_inv_gen_no_bond_types" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            bond_invariants_generator: Some(MorganBondInvariantsGenerator {
                use_bond_types: false,
                use_chirality: false,
            }),
            ..Default::default()
        },
        "r2_n2048_num_bits_per_feature_2" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            num_bits_per_feature: 2,
            ..Default::default()
        },
        "r2_n2048_additional_output" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            collect_additional_output: true,
            ..Default::default()
        },
        other => panic!("unknown Morgan fingerprint branch '{other}'"),
    }
}

fn branch_generator(
    name: &str,
    mol: &Molecule,
) -> cosmolkit_core::fingerprint::MorganFingerprintGenerator {
    let params = branch_params(name, mol);
    let atom_invariants_generator = Some(params.atom_invariants_generator.clone());
    let bond_invariants_generator =
        Some(params.bond_invariants_generator.clone().unwrap_or_else(|| {
            MorganBondInvariantsGenerator {
                use_bond_types: params.use_bond_types,
                use_chirality: params.use_chirality,
            }
        }));
    let mut generator = getMorganGeneratorWithParams(
        params.radius,
        params.count_simulation,
        params.use_chirality,
        params.use_bond_types,
        params.only_nonzero_invariants,
        params.include_redundant_environments,
        atom_invariants_generator,
        bond_invariants_generator,
        params.n_bits as u32,
        params.count_bounds.clone(),
        true,
        true,
    )
    .unwrap_or_else(|err| {
        panic!("failed to construct Morgan generator for branch {name} with {params:?}: {err}")
    });
    generator
        .fingerprint_arguments
        .fingerprint_arguments
        .d_num_bits_per_feature = params.num_bits_per_feature;
    generator
}

fn bimap(entries: &[(usize, &[(usize, u32)])]) -> BTreeMap<usize, Vec<(usize, u32)>> {
    entries
        .iter()
        .map(|(bit, values)| (*bit, values.to_vec()))
        .collect()
}

fn normalized_nonzero_elements(golden: &BTreeMap<String, i32>) -> BTreeMap<u64, i32> {
    golden
        .iter()
        .map(|(bit, count)| (bit.parse::<u64>().unwrap(), *count))
        .collect()
}

fn normalized_bit_info_map(
    golden: &BTreeMap<String, Vec<[usize; 2]>>,
) -> BTreeMap<usize, Vec<(usize, u32)>> {
    golden
        .iter()
        .map(|(bit, entries)| {
            (
                bit.parse::<usize>().unwrap(),
                entries
                    .iter()
                    .map(|entry| (entry[0], entry[1] as u32))
                    .collect(),
            )
        })
        .collect()
}

fn normalized_atoms_per_bit(
    golden: &BTreeMap<String, Vec<Vec<usize>>>,
) -> BTreeMap<usize, Vec<Vec<usize>>> {
    golden
        .iter()
        .map(|(bit, atoms_per_bit)| (bit.parse::<usize>().unwrap(), atoms_per_bit.clone()))
        .collect()
}

fn assert_additional_output_matches(
    context: &str,
    actual: &AdditionalOutput,
    expected: &MorganAdditionalOutputGolden,
) {
    let actual_atom_counts = actual
        .atom_counts
        .as_ref()
        .unwrap_or_else(|| panic!("{context}: missing atom_counts"));
    let actual_atom_to_bits = actual
        .atom_to_bits
        .as_ref()
        .unwrap_or_else(|| panic!("{context}: missing atom_to_bits"));
    let actual_bit_info_map = actual
        .bit_info_map
        .as_ref()
        .unwrap_or_else(|| panic!("{context}: missing bit_info_map"));
    let actual_atoms_per_bit = actual
        .atoms_per_bit
        .as_ref()
        .unwrap_or_else(|| panic!("{context}: missing atoms_per_bit"));
    assert_eq!(
        actual_atom_counts, &expected.atom_counts,
        "{context}: atom_counts mismatch"
    );
    assert_eq!(
        actual_atom_to_bits, &expected.atom_to_bits,
        "{context}: atom_to_bits mismatch"
    );
    let actual_bit_info_map = actual_bit_info_map
        .iter()
        .map(|(&bit, entries)| {
            (
                bit as usize,
                entries
                    .iter()
                    .map(|&(atom, layer)| (atom as usize, layer))
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    assert_eq!(
        actual_bit_info_map,
        normalized_bit_info_map(&expected.bit_info_map),
        "{context}: bit_info_map mismatch"
    );
    let actual_atoms_per_bit = actual_atoms_per_bit
        .iter()
        .map(|(&bit, atoms_per_bit)| (bit as usize, atoms_per_bit.clone()))
        .collect::<BTreeMap<_, _>>();
    assert_eq!(
        actual_atoms_per_bit,
        normalized_atoms_per_bit(&expected.atoms_per_bit),
        "{context}: atoms_per_bit mismatch"
    );
}

fn assert_optional_additional_output_matches(
    context: &str,
    args: FingerprintFuncArguments,
    expected: &Option<MorganAdditionalOutputGolden>,
) {
    if let Some(expected) = expected {
        let actual = args
            .additional_output
            .as_ref()
            .unwrap_or_else(|| panic!("{context}: COSMolKit did not return additional output"));
        assert_additional_output_matches(context, &actual, expected);
    } else {
        assert!(
            args.additional_output.is_none(),
            "{context}: unexpected additional output in COSMolKit"
        );
    }
}

fn assert_sparse_count_output_matches(
    context: &str,
    branch_name: &str,
    generator: &MorganFingerprintGenerator,
    mol: &Molecule,
    expected: &MorganSparseCountGolden,
) {
    assert!(
        expected.ok,
        "{context}: RDKit sparse-count golden recorded error {:?}",
        expected.error
    );
    let mut args = branch_func_args(branch_name, mol, expected.additional_output.is_some());
    let actual = generator
        .getSparseCountFingerprint(mol, &mut args)
        .unwrap_or_else(|err| panic!("{context}: COSMolKit sparse-count failed: {err}"));
    let expected_nonzero = expected
        .nonzero_elements
        .as_ref()
        .unwrap_or_else(|| panic!("{context}: RDKit sparse-count missing nonzero_elements"));
    assert_eq!(
        actual.nonzero_elements(),
        &normalized_nonzero_elements(expected_nonzero),
        "{context}: sparse-count nonzero elements mismatch"
    );
    assert_optional_additional_output_matches(
        &format!("{context}: sparse-count additional-output"),
        args,
        &expected.additional_output,
    );
}

fn assert_sparse_bit_output_matches(
    context: &str,
    branch_name: &str,
    generator: &MorganFingerprintGenerator,
    mol: &Molecule,
    expected: &MorganSparseBitGolden,
) {
    assert!(
        expected.ok,
        "{context}: RDKit sparse-bit golden recorded error {:?}",
        expected.error
    );
    let mut args = branch_func_args(branch_name, mol, expected.additional_output.is_some());
    let actual = generator
        .getSparseFingerprint(mol, &mut args)
        .unwrap_or_else(|err| panic!("{context}: COSMolKit sparse-bit failed: {err}"));
    let expected_on_bits = expected
        .on_bits
        .as_ref()
        .unwrap_or_else(|| panic!("{context}: RDKit sparse-bit missing on_bits"));
    assert_eq!(
        actual.on_bits().iter().copied().collect::<Vec<_>>(),
        *expected_on_bits,
        "{context}: sparse-bit on bits mismatch"
    );
    assert_eq!(
        actual.on_bits().len(),
        expected
            .num_on_bits
            .expect("RDKit ok sparse-bit should include num_on_bits"),
        "{context}: sparse-bit on-bit count mismatch"
    );
    assert_optional_additional_output_matches(
        &format!("{context}: sparse-bit additional-output"),
        args,
        &expected.additional_output,
    );
}

fn assert_hashed_count_output_matches(
    context: &str,
    branch_name: &str,
    generator: &MorganFingerprintGenerator,
    mol: &Molecule,
    expected: &MorganSparseCountGolden,
) {
    assert!(
        expected.ok,
        "{context}: RDKit hashed-count golden recorded error {:?}",
        expected.error
    );
    let mut args = branch_func_args(branch_name, mol, expected.additional_output.is_some());
    let actual = generator
        .getCountFingerprint(mol, &mut args)
        .unwrap_or_else(|err| panic!("{context}: COSMolKit hashed-count failed: {err}"));
    let expected_nonzero = expected
        .nonzero_elements
        .as_ref()
        .unwrap_or_else(|| panic!("{context}: RDKit hashed-count missing nonzero_elements"));
    assert_eq!(
        actual.nonzero_elements(),
        &normalized_nonzero_elements(expected_nonzero),
        "{context}: hashed-count nonzero elements mismatch"
    );
    assert_optional_additional_output_matches(
        &format!("{context}: hashed-count additional-output"),
        args,
        &expected.additional_output,
    );
}

fn assert_explicit_bit_output_matches(
    context: &str,
    branch_name: &str,
    generator: &MorganFingerprintGenerator,
    mol: &Molecule,
    expected: &MorganExplicitBitGolden,
) -> cosmolkit_core::Fingerprint {
    assert!(
        expected.ok,
        "{context}: RDKit explicit-bit golden recorded error {:?}",
        expected.error
    );
    let mut args = branch_func_args(branch_name, mol, expected.additional_output.is_some());
    let actual = generator
        .getFingerprint(mol, &mut args)
        .unwrap_or_else(|err| panic!("{context}: COSMolKit explicit-bit failed: {err}"));
    let expected_on_bits = expected
        .on_bits
        .as_ref()
        .unwrap_or_else(|| panic!("{context}: RDKit explicit-bit missing on_bits"));
    assert_eq!(
        actual.on_bits(),
        *expected_on_bits,
        "{context}: explicit-bit on bits mismatch"
    );
    assert_eq!(
        actual.on_bits().len(),
        expected
            .num_on_bits
            .expect("RDKit ok explicit-bit should include num_on_bits"),
        "{context}: explicit-bit on-bit count mismatch"
    );
    assert_optional_additional_output_matches(
        &format!("{context}: explicit-bit additional-output"),
        args,
        &expected.additional_output,
    );
    actual
}

#[test]
fn public_morgan_sparse_count_wrapper_matches_rdkit_golden() {
    let mol = Molecule::from_smiles("CC").unwrap();

    let output =
        morgan_get_fingerprint(&mol, 1, None, None, false, true, true, false, true, false).unwrap();

    assert_eq!(
        output.fingerprint.nonzero_elements(),
        &BTreeMap::from([(2246728737, 2), (3545175291, 1)])
    );
    assert_eq!(
        output.atoms_setting_bits,
        Some(bimap(&[
            (2246728737, &[(0, 0), (1, 0)]),
            (3545175291, &[(0, 1)])
        ]))
    );
}

#[test]
fn public_morgan_sparse_bit_wrapper_matches_rdkit_golden() {
    let mol = Molecule::from_smiles("CC").unwrap();

    let output =
        morgan_get_fingerprint(&mol, 1, None, None, false, true, false, false, true, false)
            .unwrap();

    assert_eq!(
        output.fingerprint.nonzero_elements(),
        &BTreeMap::from([(2246728737, 1), (3545175291, 1)])
    );
    assert_eq!(
        output.atoms_setting_bits,
        Some(bimap(&[
            (2246728737, &[(0, 0), (1, 0)]),
            (3545175291, &[(0, 1)])
        ]))
    );
}

#[test]
fn public_morgan_hashed_count_wrapper_matches_rdkit_golden() {
    let mol = Molecule::from_smiles("CC").unwrap();

    let output =
        morgan_get_hashed_fingerprint(&mol, 1, 2048, None, None, false, true, false, true, false)
            .unwrap();

    assert_eq!(
        output.fingerprint.nonzero_elements(),
        &BTreeMap::from([(1057, 2), (1275, 1)])
    );
    assert_eq!(
        output.atoms_setting_bits,
        Some(bimap(&[(1057, &[(0, 0), (1, 0)]), (1275, &[(0, 1)])]))
    );
}

#[test]
fn public_morgan_explicit_bit_wrapper_matches_rdkit_golden() {
    let mol = Molecule::from_smiles("CC").unwrap();

    let output = morgan_get_fingerprint_as_bit_vect(
        &mol, 1, 2048, None, None, false, true, false, true, false,
    )
    .unwrap();

    assert_eq!(output.fingerprint.on_bits(), vec![1057, 1275]);
    assert_eq!(
        output.atoms_setting_bits,
        Some(bimap(&[(1057, &[(0, 0), (1, 0)]), (1275, &[(0, 1)])]))
    );
}

#[test]
fn public_morgan_wrappers_match_custom_invariants_and_from_atoms_golden() {
    let mol = Molecule::from_smiles("CC").unwrap();

    let sparse = morgan_get_fingerprint(
        &mol,
        1,
        Some(vec![10, 20]),
        Some(vec![1]),
        false,
        true,
        true,
        false,
        true,
        true,
    )
    .unwrap();
    assert_eq!(
        sparse.fingerprint.nonzero_elements(),
        &BTreeMap::from([(20, 1), (3205493690, 1)])
    );
    assert_eq!(
        sparse.atoms_setting_bits,
        Some(bimap(&[(20, &[(1, 0)]), (3205493690, &[(1, 1)])]))
    );

    let hashed = morgan_get_hashed_fingerprint(
        &mol,
        1,
        2048,
        Some(vec![10, 20]),
        Some(vec![1]),
        false,
        true,
        false,
        true,
        true,
    )
    .unwrap();
    assert_eq!(
        hashed.fingerprint.nonzero_elements(),
        &BTreeMap::from([(20, 1), (954, 1)])
    );

    let explicit = morgan_get_fingerprint_as_bit_vect(
        &mol,
        1,
        2048,
        Some(vec![10, 20]),
        Some(vec![1]),
        false,
        true,
        false,
        true,
        true,
    )
    .unwrap();
    assert_eq!(explicit.fingerprint.on_bits(), vec![20, 954]);
}

#[test]
fn public_morgan_explicit_bit_wrapper_distinguishes_chirality_like_rdkit() {
    let r_mol = Molecule::from_smiles("C[C@H](F)Cl").unwrap();
    let s_mol = Molecule::from_smiles("C[C@@H](F)Cl").unwrap();

    let r_output = morgan_get_fingerprint_as_bit_vect(
        &r_mol, 1, 2048, None, None, true, true, false, true, false,
    )
    .unwrap();
    let s_output = morgan_get_fingerprint_as_bit_vect(
        &s_mol, 1, 2048, None, None, true, true, false, true, false,
    )
    .unwrap();

    assert_eq!(
        r_output.fingerprint.on_bits(),
        vec![1, 145, 283, 914, 1050, 1057, 1683, 1928]
    );
    assert_eq!(
        s_output.fingerprint.on_bits(),
        vec![1, 144, 283, 914, 1050, 1057, 1683, 1928]
    );
    assert_eq!(
        r_output.atoms_setting_bits,
        Some(bimap(&[
            (1, &[(1, 0)]),
            (145, &[(1, 1)]),
            (283, &[(0, 1)]),
            (914, &[(3, 1)]),
            (1050, &[(2, 1)]),
            (1057, &[(0, 0)]),
            (1683, &[(3, 0)]),
            (1928, &[(2, 0)]),
        ]))
    );
}

#[test]
fn public_morgan_count_simulation_matches_rdkit_golden() {
    let mol = Molecule::from_smiles("CC").unwrap();
    let params = MorganFingerprintParams {
        radius: 1,
        n_bits: 2048,
        count_simulation: true,
        count_bounds: vec![1, 2, 4, 8],
        include_redundant_environments: true,
        collect_additional_output: true,
        ..Default::default()
    };
    let actual = mol.morgan_fingerprint_with_output(&params).unwrap();
    let actual_additional_output = actual
        .additional_output
        .expect("COSMolKit should return additional output");

    assert_eq!(actual.fingerprint.on_bits(), vec![132, 133, 1004, 1005]);
    assert_eq!(actual_additional_output.atom_counts, vec![2, 2]);
    assert_eq!(
        actual_additional_output.atom_to_bits,
        vec![vec![132, 133, 1004, 1005], vec![132, 133, 1004, 1005]]
    );
    assert_eq!(
        actual_additional_output.bit_info_map,
        bimap(&[
            (132, &[(0, 0), (1, 0)]),
            (133, &[(0, 0), (1, 0)]),
            (1004, &[(0, 1), (1, 1)]),
            (1005, &[(0, 1), (1, 1)])
        ])
    );
    assert_eq!(
        actual_additional_output.atoms_per_bit,
        BTreeMap::from([
            (132, vec![vec![0], vec![1]]),
            (133, vec![vec![0], vec![1]]),
            (1004, vec![vec![0, 1], vec![1, 0]]),
            (1005, vec![vec![0, 1], vec![1, 0]])
        ])
    );
}

#[test]
fn morgan_fingerprint_golden_has_one_record_per_smiles() {
    let expected = parity_data::count_smiles_rows();
    let records = load_golden();
    assert_eq!(
        records.len(),
        expected,
        "Morgan fingerprint golden row count must match testdata/smiles/corpus/smiles_small.smi"
    );
    assert!(
        records.iter().any(|record| {
            record.rdkit_ok
                && record
                    .branches
                    .values()
                    .any(|branch| branch.tanimoto_to_previous.is_some())
        }),
        "Morgan fingerprint golden should include adjacent-row Tanimoto baselines"
    );
}

#[test]
fn morgan_fingerprint_matches_rdkit_golden_across_param_branches() {
    let records = load_golden();
    let mut previous_by_branch: BTreeMap<String, cosmolkit_core::Fingerprint> = BTreeMap::new();

    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}) is RDKit-not-ok but has no error",
                row_idx + 1,
                record.smiles
            );
            previous_by_branch.clear();
            continue;
        }

        let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
            panic!(
                "COSMolKit failed to parse row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });

        for (branch_name, expected_branch) in &record.branches {
            assert!(
                expected_branch.ok,
                "row {} ({}) branch {}: RDKit golden recorded error {:?}",
                row_idx + 1,
                record.smiles,
                branch_name,
                expected_branch.error
            );
            let context = format!(
                "row {} ({}) branch {}",
                row_idx + 1,
                record.smiles,
                branch_name
            );
            let generator = branch_generator(branch_name, &mol);
            assert_sparse_count_output_matches(
                &context,
                branch_name,
                &generator,
                &mol,
                &expected_branch.sparse_count,
            );
            assert_sparse_bit_output_matches(
                &context,
                branch_name,
                &generator,
                &mol,
                &expected_branch.sparse_bit,
            );
            assert_hashed_count_output_matches(
                &context,
                branch_name,
                &generator,
                &mol,
                &expected_branch.hashed_count,
            );
            let actual = assert_explicit_bit_output_matches(
                &context,
                branch_name,
                &generator,
                &mol,
                &expected_branch.explicit_bit,
            );
            let expected_on_bits = expected_branch.on_bits.as_ref().unwrap_or_else(|| {
                panic!(
                    "row {} ({}) branch {}: RDKit ok branch missing on_bits",
                    row_idx + 1,
                    record.smiles,
                    branch_name
                )
            });
            assert_eq!(
                actual.on_bits(),
                *expected_on_bits,
                "Morgan fingerprint must be RDKit bit-identical at row {} ({}) branch {}",
                row_idx + 1,
                record.smiles,
                branch_name
            );
            assert_eq!(
                actual.on_bits().len(),
                expected_branch
                    .num_on_bits
                    .expect("RDKit ok branch should include num_on_bits"),
                "Morgan fingerprint on-bit count mismatch at row {} ({}) branch {}",
                row_idx + 1,
                record.smiles,
                branch_name
            );
            assert_eq!(
                expected_branch.on_bits,
                expected_branch.explicit_bit.on_bits,
                "Morgan golden top-level on_bits must mirror explicit_bit at row {} ({}) branch {}",
                row_idx + 1,
                record.smiles,
                branch_name
            );
            assert_eq!(
                expected_branch.num_on_bits,
                expected_branch.explicit_bit.num_on_bits,
                "Morgan golden top-level num_on_bits must mirror explicit_bit at row {} ({}) branch {}",
                row_idx + 1,
                record.smiles,
                branch_name
            );
            if let Some(top_level_additional_output) = &expected_branch.additional_output {
                let explicit_additional_output = expected_branch
                    .explicit_bit
                    .additional_output
                    .as_ref()
                    .unwrap_or_else(|| {
                        panic!(
                            "Morgan golden top-level additional_output exists without explicit_bit additional_output at row {} ({}) branch {}",
                            row_idx + 1,
                            record.smiles,
                            branch_name
                        )
                    });
                assert_eq!(
                    top_level_additional_output.atom_counts,
                    explicit_additional_output.atom_counts,
                    "Morgan golden top-level additional_output must mirror explicit_bit atom_counts at row {} ({}) branch {}",
                    row_idx + 1,
                    record.smiles,
                    branch_name
                );
            }

            if let Some(expected_tanimoto) = expected_branch.tanimoto_to_previous {
                let previous = previous_by_branch.get(branch_name).unwrap_or_else(|| {
                    panic!(
                        "row {} ({}) branch {} has Tanimoto golden but no previous fingerprint",
                        row_idx + 1,
                        record.smiles,
                        branch_name
                    )
                });
                let actual_tanimoto = actual
                    .tanimoto(previous)
                    .expect("same branch fingerprints should have same bit length");
                assert!(
                    (actual_tanimoto - expected_tanimoto).abs() <= 1e-12,
                    "Morgan Tanimoto must match RDKit at row {} ({}) branch {}: actual {}, expected {}",
                    row_idx + 1,
                    record.smiles,
                    branch_name,
                    actual_tanimoto,
                    expected_tanimoto
                );
            }

            previous_by_branch.insert(branch_name.clone(), actual);
        }
    }
}
