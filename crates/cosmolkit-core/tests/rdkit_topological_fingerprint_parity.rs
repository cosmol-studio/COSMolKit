use std::collections::BTreeMap;

use cosmolkit_core::{
    Molecule, TopologicalFingerprintOutputRequest, TopologicalFingerprintParams, topological_fingerprint_with_output,
};

fn bits(smiles: &str, params: TopologicalFingerprintParams) -> Vec<usize> {
    let molecule = Molecule::from_smiles(smiles).unwrap_or_else(|err| panic!("failed to parse {smiles}: {err}"));
    topological_fingerprint_with_output(&molecule, &params, TopologicalFingerprintOutputRequest::default())
        .unwrap_or_else(|err| panic!("RDKFingerprint failed for {smiles}: {err}"))
        .fingerprint
        .on_bits()
}

fn params_for_branch(branch: &str, atom_count: usize) -> TopologicalFingerprintParams {
    let mut params = TopologicalFingerprintParams::default();
    match branch {
        "fp64" => params.fp_size = 64,
        "n1" => params.num_bits_per_feature = 1,
        "linear" => params.branched_paths = false,
        "no_bond_order" => params.use_bond_order = false,
        "from_atom_1" => params.from_atoms = Some(vec![1]),
        "custom_invariants" => {
            params.atom_invariants = Some((1..=atom_count as u32).collect());
        }
        "density" => {
            params.fp_size = 64;
            params.num_bits_per_feature = 3;
            params.target_density = 0.2;
            params.min_size = 16;
        }
        _ => {}
    }
    params
}

#[test]
fn rdkit_topological_fingerprint_matches_pinned_exact_bit_matrix() {
    let cases: &[(&str, &[(&str, Vec<usize>)])] = &[
        (
            "CCO",
            &[
                ("default", vec![562, 1183, 1308, 1339, 1728, 1772]),
                ("fp64", vec![0, 28, 31, 44, 50, 59]),
                ("n1", vec![1308, 1339, 1728]),
                ("linear", vec![562, 1183, 1308, 1339, 1728, 1772]),
                ("no_bond_order", vec![562, 1183, 1308, 1339, 1728, 1772]),
                ("from_atom_1", vec![562, 1183, 1308, 1339, 1728, 1772]),
                ("custom_invariants", vec![367, 375, 1265, 1296, 1441, 1468]),
                ("density", vec![0, 1, 8, 12, 18, 24, 27, 28, 31]),
            ],
        ),
        (
            "CC(C)C",
            &[
                ("default", vec![173, 709, 1308, 1405, 1772, 1813]),
                ("fp64", vec![5, 21, 28, 44, 45, 61]),
                ("n1", vec![1308, 1405, 1813]),
                ("linear", vec![709, 1308, 1772, 1813]),
                ("no_bond_order", vec![173, 709, 1308, 1405, 1772, 1813]),
                ("from_atom_1", vec![173, 709, 1308, 1405, 1772, 1813]),
                (
                    "custom_invariants",
                    vec![
                        283, 367, 375, 437, 661, 1265, 1296, 1441, 1468, 1524, 1685, 1784, 1941, 1986,
                    ],
                ),
                ("density", vec![5, 7, 12, 13, 21, 24, 28, 29]),
            ],
        ),
        (
            "C1CCCCC1",
            &[
                (
                    "default",
                    vec![148, 433, 709, 786, 803, 875, 1308, 1772, 1813, 1817, 1869, 1927],
                ),
                ("fp64", vec![5, 7, 13, 18, 20, 21, 25, 28, 35, 43, 44, 49]),
                ("n1", vec![148, 433, 875, 1308, 1813, 1869]),
                ("linear", vec![148, 709, 803, 875, 1308, 1772, 1813, 1817, 1869, 1927]),
                (
                    "no_bond_order",
                    vec![148, 433, 709, 786, 803, 875, 1308, 1772, 1813, 1817, 1869, 1927],
                ),
                (
                    "from_atom_1",
                    vec![148, 433, 709, 786, 803, 875, 1308, 1772, 1813, 1817, 1869, 1927],
                ),
            ],
        ),
        (
            "c1ccccc1",
            &[
                (
                    "default",
                    vec![103, 161, 194, 294, 330, 792, 842, 1026, 1287, 1784, 1889, 1907],
                ),
                ("fp64", vec![2, 7, 10, 24, 33, 38, 39, 51, 56]),
                ("n1", vec![161, 194, 294, 842, 1287, 1889]),
                ("linear", vec![103, 161, 194, 294, 330, 792, 842, 1026, 1784, 1889]),
                (
                    "no_bond_order",
                    vec![172, 402, 528, 899, 1127, 1411, 1494, 1510, 1571, 1859, 1906, 2030],
                ),
                (
                    "from_atom_1",
                    vec![103, 161, 194, 294, 330, 792, 842, 1026, 1287, 1784, 1889, 1907],
                ),
            ],
        ),
        (
            "CC.CC",
            &[
                ("default", vec![1308, 1772]),
                ("fp64", vec![28, 44]),
                ("n1", vec![1308]),
            ],
        ),
        (
            "C=C",
            &[
                ("default", vec![308, 1194]),
                ("fp64", vec![42, 52]),
                ("n1", vec![308]),
                ("no_bond_order", vec![1308, 1772]),
            ],
        ),
    ];

    for &(smiles, branches) in cases {
        let atom_count = Molecule::from_smiles(smiles).unwrap().num_atoms();
        for (branch, expected) in branches {
            assert_eq!(
                bits(smiles, params_for_branch(branch, atom_count)),
                *expected,
                "{smiles} / {branch}"
            );
        }
    }
}

#[test]
fn rdkit_topological_fingerprint_preserves_source_provenance_before_folding() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let mut params = TopologicalFingerprintParams::default();
    params.fp_size = 64;
    params.num_bits_per_feature = 1;
    let output = topological_fingerprint_with_output(
        &molecule,
        &params,
        TopologicalFingerprintOutputRequest {
            atom_bits: true,
            bit_info: true,
        },
    )
    .unwrap();
    assert_eq!(output.fingerprint.on_bits(), vec![0, 28, 59]);
    assert_eq!(
        output.output.atom_bits,
        Some(vec![vec![28, 0], vec![28, 59, 0], vec![59, 0]])
    );
    assert_eq!(
        output.output.bit_info,
        Some(BTreeMap::from([
            (0, vec![vec![0, 1]]),
            (28, vec![vec![0]]),
            (59, vec![vec![1]]),
        ]))
    );

    params.target_density = 0.2;
    params.min_size = 16;
    let folded = topological_fingerprint_with_output(
        &molecule,
        &params,
        TopologicalFingerprintOutputRequest {
            atom_bits: true,
            bit_info: true,
        },
    )
    .unwrap();
    assert_eq!(folded.fingerprint.n_bits(), 16);
    assert_eq!(folded.output.bit_info, output.output.bit_info);
}

#[test]
fn rdkit_topological_fingerprint_rejects_source_precondition_ranges() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    for (field, expected) in [
        ("min_path", "minPath==0"),
        ("max_path", "maxPath<minPath"),
        ("fp_size", "fpSize==0"),
        ("num_bits_per_feature", "nBitsPerHash==0"),
    ] {
        let mut params = TopologicalFingerprintParams::default();
        match field {
            "min_path" => params.min_path = 0,
            "max_path" => params.max_path = 0,
            "fp_size" => params.fp_size = 0,
            _ => params.num_bits_per_feature = 0,
        }
        let err = topological_fingerprint_with_output(&molecule, &params, Default::default()).unwrap_err();
        assert!(err.to_string().contains(expected), "{field}: {err}");
    }
}
