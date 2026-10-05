// Original d892 complete2 corpus conditions, both complete original datasets.
use std::collections::BTreeMap;

use cosmolkit::{Molecule, TopologicalFingerprintParams};
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::sync::OnceLock;

#[derive(Debug, Deserialize)]
struct GoldenRecord {
    row: usize,
    smiles: String,
    rdkit_ok: bool,
    branches: BTreeMap<String, GoldenBranch>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct GoldenBranch {
    parameters: GoldenParameters,
    ok: bool,
    on_bits: Option<Vec<usize>>,
    num_bits: Option<usize>,
    num_on_bits: Option<usize>,
    error: Option<String>,
}

#[derive(Debug, Deserialize, serde::Serialize)]
struct GoldenParameters {
    name: String,
    #[serde(rename = "minPath")]
    min_path: u32,
    #[serde(rename = "maxPath")]
    max_path: u32,
    #[serde(rename = "fpSize")]
    fp_size: u32,
    #[serde(rename = "nBitsPerHash")]
    num_bits_per_feature: u32,
    #[serde(rename = "useHs")]
    use_hs: bool,
    #[serde(rename = "tgtDensity")]
    target_density: f64,
    #[serde(rename = "minSize")]
    min_size: u32,
    #[serde(rename = "branchedPaths")]
    branched_paths: bool,
    #[serde(rename = "useBondOrder")]
    use_bond_order: bool,
    #[serde(rename = "fromAtoms")]
    from_atoms: Option<String>,
    #[serde(rename = "atomInvariants")]
    atom_invariants: Option<String>,
}

struct Profile {
    name: &'static str,
    corpus: Vec<String>,
    golden: Vec<GoldenRecord>,
}
fn sha(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}
fn selected_profiles() -> &'static [Profile] {
    static ALL: OnceLock<Vec<Profile>> = OnceLock::new();
    ALL.get_or_init(|| {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let profile_bytes = std::fs::read(
            root.join("tools/testdata/rdkit/rdkit_topological_fingerprint_profile.json"),
        )
        .unwrap();
        assert_eq!(
            sha(&profile_bytes),
            "227d51b15c6387885fd38502ed69b8a48691700d343f2fd3103ae9e0738f34ce",
            "original14 profile bytes"
        );
        let profile: serde_json::Value = serde_json::from_slice(&profile_bytes).unwrap();
        let branches = profile["branches"].as_array().unwrap();
        assert_eq!(branches.len(), 14);
        [
            (
                "small152",
                "COSMOLKIT_RDK_ORIGINAL_SMALL_GOLDEN",
                "1a1e1fd9bcb60e702814d984434c68a1bb5726c045f7bb1bf2b7af7ed3f11fb1",
                "testdata/smiles/corpus/smiles_small.smi",
                "47380e477dc2ab4b3c2b7cd62754e52b718bb9f3f0b977d47eb97db45074442e",
                152,
            ),
            (
                "original5000",
                "COSMOLKIT_RDK_ORIGINAL5000_GOLDEN",
                "0d443196f862ad5d2fa0c9be3c16b8fc3464800d0ea64f87f894b2eba1f4cc69",
                "testdata/smiles/corpus/smiles_5000.smi",
                "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849",
                5000,
            ),
        ]
        .into_iter()
        .map(
            |(name, variable, golden_sha, corpus_file, corpus_sha, count)| {
                let path = std::env::var(variable).unwrap_or_else(|_| {
                    panic!("prepare pinned original14 references and set {variable}")
                });
                let bytes = std::fs::read(path).unwrap();
                assert_eq!(sha(&bytes), golden_sha, "{name}: immutable golden bytes");
                let golden: Vec<GoldenRecord> = bytes
                    .split(|&b| b == b'\n')
                    .filter(|l| !l.is_empty())
                    .map(|l| serde_json::from_slice(l).unwrap())
                    .collect();
                let corpus_bytes = std::fs::read(root.join(corpus_file)).unwrap();
                assert_eq!(
                    sha(&corpus_bytes),
                    corpus_sha,
                    "{name}: original input bytes"
                );
                let corpus: Vec<String> = std::str::from_utf8(&corpus_bytes)
                    .unwrap()
                    .lines()
                    .map(str::trim)
                    .filter(|l| !l.is_empty() && !l.starts_with('#'))
                    .map(str::to_owned)
                    .collect();
                assert_eq!(corpus.len(), count);
                assert_eq!(golden.len(), count);
                for (row, (record, smiles)) in golden.iter().zip(&corpus).enumerate() {
                    assert_eq!(record.row, row);
                    assert_eq!(&record.smiles, smiles);
                    if record.rdkit_ok {
                        assert!(record.error.is_none());
                        assert_eq!(record.branches.len(), 14);
                        for parameter in branches {
                            let branch_name = parameter["name"].as_str().unwrap();
                            let branch = record.branches.get(branch_name).unwrap();
                            assert_eq!(branch.parameters.name, branch_name);
                            // Validate every original parameter field before any CK call.
                            let actual_parameter =
                                serde_json::to_value(&branch.parameters).unwrap();
                            for (k, v) in parameter.as_object().unwrap() {
                                assert_eq!(
                                    &actual_parameter[k], v,
                                    "{name} row {row} {branch_name} {k}"
                                );
                            }
                            assert!(branch.ok);
                            assert!(branch.error.is_none());
                            assert_eq!(branch.num_on_bits, branch.on_bits.as_ref().map(Vec::len));
                            assert!(branch.num_bits.is_some());
                        }
                    } else {
                        assert!(record.error.is_some());
                        assert!(record.branches.is_empty());
                    }
                }
                Profile {
                    name,
                    corpus,
                    golden,
                }
            },
        )
        .collect()
    })
}

fn params_from_golden(
    parameters: &GoldenParameters,
    atom_count: usize,
) -> TopologicalFingerprintParams {
    let mut params = TopologicalFingerprintParams {
        min_path: parameters.min_path,
        max_path: parameters.max_path,
        fp_size: parameters.fp_size,
        num_bits_per_feature: parameters.num_bits_per_feature,
        use_hs: parameters.use_hs,
        target_density: parameters.target_density,
        min_size: parameters.min_size,
        branched_paths: parameters.branched_paths,
        use_bond_order: parameters.use_bond_order,
        ..Default::default()
    };
    if parameters.from_atoms.as_deref() == Some("first") && atom_count > 0 {
        params.from_atoms = Some(vec![0]);
    }
    if parameters.atom_invariants.as_deref() == Some("index_plus_one") {
        params.atom_invariants = Some((1..=atom_count as u32).collect());
    }
    params
}

#[test]
fn rdkit_topological_fingerprint_golden_has_one_record_per_active_corpus_row() {
    for profile in selected_profiles() {
        let corpus = &profile.corpus;
        let golden = &profile.golden;
        eprintln!(
            "complete original14 profile {} / {} rows",
            profile.name,
            corpus.len()
        );
        assert_eq!(
            golden.len(),
            corpus.len(),
            "RDKFingerprint golden must have one row per active corpus input"
        );
        for (row, (record, smiles)) in golden.iter().zip(corpus).enumerate() {
            assert_eq!(record.row, row, "golden row index changed at row {row}");
            assert_eq!(record.smiles, *smiles, "golden SMILES changed at row {row}");
            if record.rdkit_ok {
                assert!(
                    record.error.is_none(),
                    "RDKit-success row {row} has an error"
                );
                assert!(
                    !record.branches.is_empty(),
                    "RDKit-success row {row} has no branches"
                );
            } else {
                assert!(
                    record.branches.is_empty(),
                    "RDKit-failed row {row} has branches"
                );
                assert!(
                    record.error.is_some(),
                    "RDKit-failed row {row} has no error"
                );
            }
        }
    }
}

#[test]
fn rdkit_topological_fingerprint_matches_every_active_golden_profile_exactly() {
    for profile in selected_profiles() {
        let corpus = &profile.corpus;
        let golden = &profile.golden;
        eprintln!(
            "complete original14 profile {} / {} rows",
            profile.name,
            corpus.len()
        );
        assert_eq!(golden.len(), corpus.len());

        let expected_branch_names = golden
            .first()
            .map(|record| record.branches.keys().cloned().collect::<Vec<_>>())
            .unwrap_or_default();
        assert!(
            !expected_branch_names.is_empty(),
            "RDKFingerprint profile is empty"
        );

        golden
        .iter()
        .zip(corpus.iter())
        .enumerate()
        .for_each(|(row, (record, smiles))| {
        if !record.rdkit_ok {
            let parse_result = Molecule::from_smiles(smiles);
            assert!(
                parse_result.is_err(),
                "row {row} RDKit rejected but COSMolKit parsed"
            );
            return;
        }

        let molecule = Molecule::from_smiles(smiles)
            .unwrap_or_else(|error| panic!("row {row} ({smiles}) failed to parse: {error}"));
        assert_eq!(
            record.branches.keys().cloned().collect::<Vec<_>>(),
            expected_branch_names,
            "RDKFingerprint branch set changed at row {row}"
        );

        for branch_name in &expected_branch_names {
            let branch = record
                .branches
                .get(branch_name)
                .unwrap_or_else(|| panic!("row {row} missing branch {branch_name}"));
            assert_eq!(
                branch.parameters.name, *branch_name,
                "row {row} branch metadata mismatch"
            );
            let params = params_from_golden(&branch.parameters, molecule.num_atoms());
            let actual = molecule.topological_fingerprint_with_params(&params);

            if !branch.ok {
                assert!(
                    actual.is_err(),
                    "row {row} ({smiles}) branch {branch_name} unexpectedly succeeded; RDKit error: {:?}",
                    branch.error
                );
                continue;
            }

            let actual = actual.unwrap_or_else(|error| {
                panic!("row {row} ({smiles}) branch {branch_name} failed: {error}")
            });
            assert_eq!(
                actual.n_bits() as usize,
                branch
                    .num_bits
                    .unwrap_or_else(|| panic!("row {row} missing num_bits")),
                "row {row} ({smiles}) branch {branch_name} fingerprint size mismatch"
            );
            assert_eq!(
                actual.on_bits().into_iter().map(|b| b as usize).collect::<Vec<_>>(),
                branch
                    .on_bits
                    .clone()
                    .unwrap_or_else(|| panic!("row {row} missing on_bits")),
                "row {row} ({smiles}) branch {branch_name} exact bits mismatch"
            );
            assert_eq!(
                actual.on_bits().len(),
                branch
                    .num_on_bits
                    .unwrap_or_else(|| panic!("row {row} missing num_on_bits")),
                "row {row} ({smiles}) branch {branch_name} on-bit count mismatch"
            );
        }
    });
    }
}
