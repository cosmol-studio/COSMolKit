//! Explicit local tests; preparing references compiles official InChI C.
use cosmolkit::{InchiError, InchiErrorKind, Molecule, SdfReadParams, inchi_to_key};
use serde_json::{Value, json};
use std::io::{BufRead, BufReader, Write};

fn generation_error(error: InchiError) -> Value {
    if error.kind == InchiErrorKind::Toolkit {
        if let Some(detail) = error.detail.strip_prefix("ToolkitError: ") {
            if detail.starts_with("Explicit valence for atom # ") {
                return json!({"stage":"inchi_error", "kind":"valence", "detail":detail});
            }
            if let Some(atoms) = detail
                .strip_prefix("could not kekulize molecule; remaining atoms: [")
                .and_then(|s| s.strip_suffix(']'))
            {
                let parsed: Option<Vec<usize>> = atoms
                    .split(", ")
                    .map(|atom| {
                        atom.strip_prefix("AtomId(")?
                            .strip_suffix(')')?
                            .parse()
                            .ok()
                    })
                    .collect();
                if let Some(atoms) = parsed {
                    return json!({"stage":"inchi_error", "kind":"kekulize", "atoms":atoms});
                }
            }
        }
    }
    json!({"stage":"unexpected_error", "detail":error.to_string()})
}

fn compare(permutations: bool) {
    let path = std::env::var("CK_INCHI_OFFICIAL_REFERENCE")
        .expect("prepare tests/support/prepare_official.py and set CK_INCHI_OFFICIAL_REFERENCE");
    let input = std::fs::File::open(&path).expect("open prepared reference");
    let report_path = format!(
        "{path}.{}.report.jsonl",
        if permutations { "permutations" } else { "sdf" }
    );
    let mut report = std::fs::File::create(&report_path).unwrap();
    let (
        mut matched,
        mut mismatches,
        mut rejected,
        mut panics,
        mut unprepared,
        mut upstream_unstable,
    ) = (0, 0, 0, 0, 0, 0);
    let (mut expected_path_difference, mut matching_errors) = (0, 0);
    let mut baseline: Option<(Value, Value)> = None;
    for line in BufReader::new(input).lines() {
        let row: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let variant = row["variant"].as_u64().unwrap();
        if !permutations && variant != 0 {
            continue;
        }
        // The official consumer compares only the reserialized permutations.
        // Comparing variant 0 here would conflate SDF reader/writer changes
        // with atom-order invariance in the InChI engine.
        if variant == 1 {
            baseline = Some((row["inchi"].clone(), row["key"].clone()));
        }
        let Some(block) = row["molblock"].as_str() else {
            unprepared += 1;
            writeln!(report, "{}", row).unwrap();
            continue;
        };
        assert_eq!(
            row["rdkit_version"], "2026.09.1",
            "prepare the pinned RDKit companion reference"
        );
        let rdkit = &row["rdkit"];
        let expected = (row["inchi"].as_str().unwrap(), row["key"].as_str().unwrap());
        let actual = std::panic::catch_unwind(|| -> Value {
            // Preserve explicit hydrogens and avoid adding a separate sanitization policy.
            let params = SdfReadParams {
                sanitize: false,
                remove_hs: false,
                strict_parsing: false,
                ..Default::default()
            };
            let mol = match Molecule::from_sdf_with_params(block, &params) {
                Ok(mol) => mol,
                Err(error) => return json!({"stage":"sdf_error", "detail":error.to_string()}),
            };
            let result = match mol.to_inchi_result() {
                Ok(result) => result,
                Err(error) => return generation_error(error),
            };
            let key = if result.inchi().is_empty() {
                String::new()
            } else {
                match inchi_to_key(result.inchi()) {
                    Ok(key) => key,
                    Err(error) => return json!({"stage":"key_error", "detail":error.to_string()}),
                }
            };
            json!({"stage":"output", "inchi":result.inchi(), "key":key})
        });
        if permutations
            && variant > 1
            && !expected.0.is_empty()
            && baseline
                .as_ref()
                .is_some_and(|(i, k)| i != &row["inchi"] || k != &row["key"])
        {
            upstream_unstable += 1;
            writeln!(
                report,
                "{}",
                json!({"category":"upstream_noninvariance", "reference":row})
            )
            .unwrap();
        }
        let category = match &actual {
            Err(_) => {
                panics += 1;
                "native_panic"
            }
            Ok(_)
                if rdkit["stage"] == "sdf_rejected"
                    || (rdkit["stage"] == "output" && rdkit["inchi"] == "") =>
            {
                rejected += 1;
                "source_rejected_no_output_oracle"
            }
            Ok(value) if value == rdkit && rdkit["stage"] == "inchi_error" => {
                matching_errors += 1;
                "matching_source_error"
            }
            Ok(value)
                if value == rdkit && value["inchi"] == expected.0 && value["key"] == expected.1 =>
            {
                matched += 1;
                "matched"
            }
            Ok(value) if value == rdkit => {
                expected_path_difference += 1;
                "rdkit_matched_direct_sdf_difference"
            }
            _ => {
                mismatches += 1;
                "mismatch"
            }
        };
        if category != "matched" {
            writeln!(
                report,
                "{}",
                json!({"category":category,"reference":row,"actual":format!("{actual:?}")})
            )
            .unwrap();
        }
    }
    eprintln!(
        "matched={matched} rdkit_matched_direct_sdf_difference={expected_path_difference} matching_source_errors={matching_errors} mismatches={mismatches} source_rejected={rejected} native_panics={panics} unprepared={unprepared} upstream_noninvariance={upstream_unstable}; report={report_path}"
    );
    assert!(matched > 0, "no completed output comparisons");
    assert_eq!(mismatches + panics + unprepared, 0, "see retained report");
}

#[test]
#[ignore = "local only: requires explicitly prepared official C oracle"]
fn official_sdf_outputs() {
    compare(false);
}

#[test]
#[ignore = "local only: requires explicitly prepared official C oracle and permutations"]
fn official_atom_permutation_invariance() {
    compare(true);
}
