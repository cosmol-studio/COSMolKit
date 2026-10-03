use super::*;
use registry::Pair;
use std::cell::Cell;

fn fingerprint_tasks() -> Result<Vec<&'static Task>> {
    Ok(registry::TASKS
        .iter()
        .filter(|t| !matches!(t.operation, registry::Operation::Molecular(_)))
        .collect())
}
fn fingerprint_corpus() -> Corpus {
    Corpus {
        fingerprints: registry::builtin(),
        molecules: vec![],
    }
}
fn fingerprint(input: &Input) -> &registry::FingerprintInput {
    let Input::Fingerprint(input) = input else {
        panic!("fingerprint expected")
    };
    input
}

#[test]
fn public_input_adapter_preserves_extreme_counts_and_stored_zero() {
    let entries = vec![(1, i32::MIN), (2, i32::MAX), (3, 0)];
    let cases = Corpus {
        fingerprints: vec![Pair {
            id: "adapter_boundaries".into(),
            length: 4,
            left: entries.clone(),
            right: entries.clone(),
        }],
        molecules: vec![],
    };
    let tasks = fingerprint_tasks().unwrap();
    registry::validate(&cases, &tasks).unwrap();
    for task in tasks {
        for input in registry::expand(&cases, task) {
            assert_eq!(
                execute::run(&input).unwrap().output,
                registry::Value::Fingerprint(registry::FingerprintValue {
                    length: 4,
                    entries: entries.clone()
                })
            );
        }
    }
}

// Framework tests intentionally use synthetic reference values. They test
// validation/comparison machinery, NOT RDKit parity; the CLI uses the oracle.
fn fixture(data: &Path, task: &Task, cases: &Corpus) -> PathBuf {
    let inputs = registry::expand(cases, task);
    let records = synthetic(&inputs).unwrap();
    publish(data, task, &inputs, &records).unwrap();
    generation(data, task, &encode(&inputs).unwrap())
}

// Only framework tests compose these two explicit stages, using synthetic
// values. Production tests expose no preparation/generator entrypoint.
fn prepare_then_compare(
    tasks: &[&Task],
    cases: &Corpus,
    data: &Path,
    mut generate: impl FnMut(&[Input]) -> Result<Vec<Record>>,
    executor: impl FnMut(&Input) -> Result<Record>,
) -> Result<Vec<Comparison>> {
    prepare_with(tasks, cases, data, |_, inputs| generate(inputs))?;
    run(tasks, cases, data, executor)
}

#[test]
fn missing_last_task_stops_before_first_operation() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    fixture(temp.path(), tasks[0], &cases);
    let calls = Cell::new(0);
    let result = prepare_then_compare(
        &tasks,
        &cases,
        temp.path(),
        |_| Err("oracle unavailable".into()),
        |i| {
            calls.set(calls.get() + 1);
            execute::run(i)
        },
    );
    assert!(result.unwrap_err().contains("fuzzy_or"));
    assert_eq!(calls.get(), 0);
}

fn synthetic(inputs: &[Input]) -> Result<Vec<Record>> {
    Ok(inputs
        .iter()
        .map(|input| Record {
            input: input.clone(),
            output: registry::Value::Fingerprint(registry::FingerprintValue {
                length: fingerprint(input).case.length,
                entries: vec![],
            }),
        })
        .collect())
}

#[test]
fn explicit_preparation_finishes_before_comparison_then_reuses_cache() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    let generated = Cell::new(0);
    let executed = Cell::new(0);
    let report = prepare_then_compare(
        &tasks,
        &cases,
        temp.path(),
        |inputs| {
            generated.set(generated.get() + 1);
            synthetic(inputs)
        },
        |input| {
            assert_eq!(generated.get(), tasks.len());
            executed.set(executed.get() + 1);
            Ok(synthetic(std::slice::from_ref(input))?.remove(0))
        },
    )
    .unwrap();
    assert_eq!(executed.get(), cases.fingerprints.len() * 4);
    assert!(report.iter().all(|row| row.matches));
    let (_, preparation) = prepare_with(&tasks, &cases, temp.path(), |_, _| {
        panic!("valid cache must not invoke oracle")
    })
    .unwrap();
    assert_eq!(
        preparation,
        Preparation {
            reused_tasks: 2,
            generated_tasks: 0,
            rows: report.len()
        }
    );
}

#[test]
fn repairs_only_corrupt_task_and_preserves_evidence() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    fixture(temp.path(), tasks[0], &cases);
    let broken = fixture(temp.path(), tasks[1], &cases);
    fs::write(broken.join("reference.json"), b"[]").unwrap();
    let (_, preparation) = prepare_with(&tasks, &cases, temp.path(), |_, inputs| {
        assert_eq!(
            fingerprint(&inputs[0]).operation,
            registry::Operation::FuzzyOr
        );
        synthetic(inputs)
    })
    .unwrap();
    assert_eq!(preparation.generated_tasks, 1);
    assert_eq!(preparation.reused_tasks, 1);
    let backup = fs::read_dir(temp.path())
        .unwrap()
        .map(|e| e.unwrap().path())
        .find(|p| {
            p.file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with(".invalid-")
        })
        .unwrap();
    assert_eq!(read(&backup.join("reference.json")).unwrap(), b"[]");
    assert_eq!(
        preflight(&tasks, &cases, temp.path()).unwrap().len(),
        cases.fingerprints.len() * 4
    );
}

#[test]
fn malformed_oracle_never_publishes_or_executes() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    let result = prepare_then_compare(
        &tasks,
        &cases,
        temp.path(),
        |_| Ok(vec![]),
        |_| panic!("must not execute"),
    );
    assert!(result.unwrap_err().contains("reference row count"));
    assert!(
        !generation(
            temp.path(),
            tasks[0],
            &encode(&registry::expand(&cases, tasks[0])).unwrap()
        )
        .exists()
    );
}

#[test]
fn failed_final_generation_does_not_execute_prepared_first_task() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    let result = prepare_then_compare(
        &tasks,
        &cases,
        temp.path(),
        |inputs| {
            if fingerprint(&inputs[0]).operation == registry::Operation::FuzzyOr {
                Err("oracle failure".into())
            } else {
                synthetic(inputs)
            }
        },
        |_| panic!("global barrier must prevent execution"),
    );
    assert!(result.unwrap_err().contains("0 Rust operation calls"));
    assert!(preflight(&tasks[..1], &cases, temp.path()).is_ok());
}

#[test]
fn stale_manifest_is_repaired_and_unselected_tasks_are_not_generated() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(Some("fuzzy_and")).unwrap();
    let cases = fingerprint_corpus();
    let directory = fixture(temp.path(), tasks[0], &cases);
    let mut manifest: Manifest =
        serde_json::from_slice(&read(&directory.join("manifest.json")).unwrap()).unwrap();
    manifest.rdkit_version = "stale".into();
    fs::write(directory.join("manifest.json"), encode(&manifest).unwrap()).unwrap();
    let (_, preparation) = prepare_with(&tasks, &cases, temp.path(), |_, inputs| {
        assert_eq!(
            fingerprint(&inputs[0]).operation,
            registry::Operation::FuzzyAnd
        );
        synthetic(inputs)
    })
    .unwrap();
    assert_eq!(preparation.generated_tasks, 1);
    assert_eq!(preparation.rows, cases.fingerprints.len() * 2);
    assert!(preflight(&fingerprint_tasks().unwrap(), &cases, temp.path()).is_err());
}

#[test]
fn comparison_failures_never_rewrite_reference_data() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(Some("fuzzy_and")).unwrap();
    let cases = fingerprint_corpus();
    let directory = fixture(temp.path(), tasks[0], &cases);
    let before = read(&directory.join("reference.json")).unwrap();
    let report = prepare_then_compare(
        &tasks,
        &cases,
        temp.path(),
        |_| panic!("no oracle for valid cache"),
        |_| Err("CK failed".into()),
    )
    .unwrap();
    assert!(report.iter().all(|row| !row.matches));
    assert_eq!(read(&directory.join("reference.json")).unwrap(), before);
}

#[test]
fn corrupted_reference_stops_before_execution() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    fixture(temp.path(), tasks[0], &cases);
    let last = fixture(temp.path(), tasks[1], &cases);
    fs::write(last.join("reference.json"), b"[]").unwrap();
    assert!(
        preflight(&tasks, &cases, temp.path())
            .err()
            .unwrap()
            .contains("corrupted")
    );
}

#[test]
fn selected_task_does_not_require_unselected_task() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(Some("fuzzy_and")).unwrap();
    let cases = fingerprint_corpus();
    fixture(temp.path(), tasks[0], &cases);
    assert_eq!(
        preflight(&tasks, &cases, temp.path()).unwrap().len(),
        cases.fingerprints.len() * 2
    );
}

#[test]
fn wrong_case_identity_is_rejected_even_with_updated_checksum() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(Some("fuzzy_or")).unwrap();
    let cases = fingerprint_corpus();
    let directory = fixture(temp.path(), tasks[0], &cases);
    let mut records: Vec<LabeledRecord> =
        serde_json::from_slice(&read(&directory.join("reference.json")).unwrap()).unwrap();
    if let Input::Fingerprint(input) = &mut records[0].record.input {
        input.width = registry::Width::U64;
    }
    let reference = encode(&records).unwrap();
    let input = read(&directory.join("input.json")).unwrap();
    fs::write(directory.join("reference.json"), &reference).unwrap();
    fs::write(
        directory.join("manifest.json"),
        encode(&identity(tasks[0], &input, &reference, records.len())).unwrap(),
    )
    .unwrap();
    assert!(
        preflight(&tasks, &cases, temp.path())
            .err()
            .unwrap()
            .contains("case/parameter")
    );
}

#[test]
fn comparison_checks_values_and_does_not_drop_failures() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    for task in &tasks {
        fixture(temp.path(), task, &cases);
    }
    let ready = preflight(&tasks, &cases, temp.path()).unwrap();
    let expected = ready.len();
    let report = compare(ready, |input| {
        Ok(Record {
            input: input.clone(),
            output: registry::Value::Fingerprint(registry::FingerprintValue {
                length: fingerprint(input).case.length,
                entries: vec![(1, 99)],
            }),
        })
    });
    assert_eq!(report.len(), expected);
    assert!(report.iter().all(|row| !row.matches));
}

#[test]
fn registry_rejects_unknown_task_and_invalid_corpus() {
    assert!(registry::select(Some("fuzzy_annd")).is_err());
    let tasks = fingerprint_tasks().unwrap();
    assert!(registry::validate(&Corpus::default(), &tasks).is_err());
    let mut cases = fingerprint_corpus();
    cases.fingerprints.push(cases.fingerprints[0].clone());
    assert!(registry::validate(&cases, &tasks).is_err());
    let mut cases = fingerprint_corpus();
    cases.fingerprints[0].left = vec![(16, 3)];
    assert!(registry::validate(&cases, &tasks).is_err());
}

#[test]
fn molecular_missing_final_reference_blocks_all_rust_calls() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(None).unwrap();
    let mut cases = fingerprint_corpus();
    cases.molecules.push(registry::SmilesCase {
        id: "ethanol".into(),
        smiles: "CCO".into(),
    });
    let result = prepare_then_compare(
        &tasks,
        &cases,
        temp.path(),
        |inputs| {
            if let Input::Molecular {
                profile: registry::molecule_plan::Profile::RemoveHydrogens { .. },
                ..
            } = &inputs[0]
            {
                return Err("last reference missing".into());
            }
            inputs
                .iter()
                .map(|input| {
                    let output = match input {
                        Input::Fingerprint(input) => {
                            registry::Value::Fingerprint(registry::FingerprintValue {
                                length: input.case.length,
                                entries: vec![],
                            })
                        }
                        Input::Molecular { .. } => {
                            registry::Value::Molecular(molecular::Outcome::Error {
                                stage: molecular::Stage::Parse,
                                detail: "synthetic diagnostic".into(),
                            })
                        }
                    };
                    Ok(Record {
                        input: input.clone(),
                        output,
                    })
                })
                .collect()
        },
        |_| panic!("global molecular barrier bypassed"),
    );
    assert!(result.unwrap_err().contains("last reference missing"));
}

#[test]
fn molecular_schema_and_comparison_do_not_accept_wrong_types_or_equal_errors() {
    use registry::molecule_plan::Profile;
    assert!(
        molecular::validate_output(
            &Profile::MolecularWeight { only_heavy: false },
            &molecular::Outcome::Text("18".into())
        )
        .is_err()
    );
    let input = Input::Molecular {
        case: registry::SmilesCase {
            id: "invalid".into(),
            smiles: "C(".into(),
        },
        profile: Profile::SanitizeAll,
    };
    let record = Record {
        input,
        output: registry::Value::Molecular(molecular::Outcome::Error {
            stage: molecular::Stage::Parse,
            detail: "same text".into(),
        }),
    };
    let ready = Ready {
        records: label_records(
            registry::select(Some("sanitize")).unwrap()[0],
            std::slice::from_ref(&record),
        )
        .unwrap(),
    };
    assert!(!compare(ready, |_| Ok(record.clone()))[0].matches);
}

#[test]
fn descriptor_query_parity_heavy_registry_matrix_expands_both_policies_in_order() {
    use registry::molecule_plan::{Category, Profile, TaskId};
    let task = registry::select(Some("num_heavy_atoms")).unwrap();
    assert_eq!(task.len(), 1);
    assert_eq!(task[0].operation.name(), "num_heavy_atoms");
    assert_eq!(TaskId::NumHeavyAtoms.name(), "num_heavy_atoms");
    assert_eq!(TaskId::NumHeavyAtoms.category(), Category::Descriptors);
    // Profile expansion order is false then true (INPUT preparation
    // axis; constructor sanitize=true lives in the runner, not here).
    let profiles = TaskId::NumHeavyAtoms.profiles();
    assert_eq!(
        profiles,
        vec![
            Profile::NumHeavyAtoms {
                remove_hydrogens: false
            },
            Profile::NumHeavyAtoms {
                remove_hydrogens: true
            },
        ]
    );
    // Exact ten-profile global expansion is checked per selected task:
    // two cases x two policies.
    let mut cases = fingerprint_corpus();
    cases.molecules.push(registry::SmilesCase {
        id: "ethanol".into(),
        smiles: "CCO".into(),
    });
    cases.molecules.push(registry::SmilesCase {
        id: "methane".into(),
        smiles: "C".into(),
    });
    let inputs = registry::expand(&cases, task[0]);
    assert_eq!(inputs.len(), 4);
    for input in &inputs {
        assert_eq!(input.task_name(), "num_heavy_atoms");
    }
    assert_eq!(task[0].count(&cases), 4);
}

/// RING-LIVE-PUBLIC T1 synthetic coverage: eleven registrations x two
/// profiles, keys/generators/outcome kinds, actual adapter routing on
/// fixed benzene/pyridine/cyclohexane/dummy/H/D inputs, missing-data
/// global barrier, wrong-kind rejection, and every original registered
/// task retained.
#[test]
fn ring_descriptor_framework_registers_eleven_tasks_with_two_profiles_each() {
    use registry::molecule_plan::{Category, Profile, TaskId};
    const RING_TASKS: [(&str, TaskId); 11] = [
        ("num_rings", TaskId::NumRings),
        ("num_heterocycles", TaskId::NumHeterocycles),
        ("num_aromatic_rings", TaskId::NumAromaticRings),
        ("num_saturated_rings", TaskId::NumSaturatedRings),
        ("num_aliphatic_rings", TaskId::NumAliphaticRings),
        ("num_aromatic_heterocycles", TaskId::NumAromaticHeterocycles),
        ("num_aromatic_carbocycles", TaskId::NumAromaticCarbocycles),
        (
            "num_aliphatic_heterocycles",
            TaskId::NumAliphaticHeterocycles,
        ),
        ("num_aliphatic_carbocycles", TaskId::NumAliphaticCarbocycles),
        (
            "num_saturated_heterocycles",
            TaskId::NumSaturatedHeterocycles,
        ),
        ("num_saturated_carbocycles", TaskId::NumSaturatedCarbocycles),
    ];
    const GENERATORS: [&str; 11] = [
        "generate_num_rings",
        "generate_num_heterocycles",
        "generate_num_aromatic_rings",
        "generate_num_saturated_rings",
        "generate_num_aliphatic_rings",
        "generate_num_aromatic_heterocycles",
        "generate_num_aromatic_carbocycles",
        "generate_num_aliphatic_heterocycles",
        "generate_num_aliphatic_carbocycles",
        "generate_num_saturated_heterocycles",
        "generate_num_saturated_carbocycles",
    ];
    // 11 tasks x 2 profiles = 22 registrations; every original task
    // (fingerprint, notation, chemistry, descriptor, stereo, depiction)
    // is still present.
    let total_ring = RING_TASKS
        .iter()
        .map(|(name, _)| registry::select(Some(name)).unwrap().len())
        .sum::<usize>();
    assert_eq!(total_ring, 11);
    for (index, (name, id)) in RING_TASKS.iter().enumerate() {
        let task = &registry::select(Some(name)).unwrap()[0];
        assert_eq!(task.operation.name(), *name);
        assert_eq!(task.generator, GENERATORS[index]);
        assert!(matches!(task.corpus_type, registry::CorpusType::Smiles));
        assert_eq!(id.name(), *name);
        assert_eq!(id.category(), Category::Descriptors);
        // Two explicit profiles: remove_hydrogens false then true,
        // sanitize=true handled by the public runner.
        assert_eq!(id.profiles().len(), 2);
    }
    // Original registered tasks retained.
    assert!(
        registry::TASKS
            .iter()
            .any(|t| t.operation.name() == "num_heavy_atoms")
    );
    assert!(
        registry::TASKS
            .iter()
            .any(|t| t.operation.name() == "sanitize")
    );
    assert!(
        registry::TASKS
            .iter()
            .any(|t| t.operation.name() == "fuzzy_and")
    );
}

#[test]
fn ring_descriptor_framework_execution_routes_public_queries_and_rejects_wrong_kind() {
    use registry::molecule_plan::Profile;
    // Fixed inputs: benzene (aromatic carbocycle), pyridine (aromatic
    // heterocycle), cyclohexane (saturated carbocycle), dummy ring,
    // explicit-H and deuterated cycles.
    let inputs = [
        (
            "benzene",
            "c1ccccc1",
            Profile::NumRings {
                remove_hydrogens: false,
            },
            1,
        ),
        (
            "pyridine",
            "n1ccccc1",
            Profile::NumHeterocycles {
                remove_hydrogens: false,
            },
            1,
        ),
        (
            "cyclohexane",
            "C1CCCCC1",
            Profile::NumSaturatedRings {
                remove_hydrogens: false,
            },
            1,
        ),
        (
            "benzene-aromatic",
            "c1ccccc1",
            Profile::NumAromaticCarbocycles {
                remove_hydrogens: false,
            },
            1,
        ),
        (
            "pyridine-aromatic-hetero",
            "n1ccccc1",
            Profile::NumAromaticHeterocycles {
                remove_hydrogens: false,
            },
            1,
        ),
        (
            "cyclohexane-aliphatic-carbo",
            "C1CCCCC1",
            Profile::NumAliphaticCarbocycles {
                remove_hydrogens: false,
            },
            1,
        ),
        (
            "dummy-hetero",
            "*1CCCC1",
            Profile::NumSaturatedHeterocycles {
                remove_hydrogens: false,
            },
            1,
        ),
        (
            "dummy-aliphatic-hetero",
            "*1CCCC1",
            Profile::NumAliphaticHeterocycles {
                remove_hydrogens: false,
            },
            1,
        ),
        (
            "h-cycle-rings",
            "[H]C1CCCCC1",
            Profile::NumRings {
                remove_hydrogens: false,
            },
            1,
        ),
        (
            "h-cycle-rings-removed",
            "[H]C1CCCCC1",
            Profile::NumRings {
                remove_hydrogens: true,
            },
            1,
        ),
        (
            "d-cycle-saturated-carbo",
            "[2H]C1CCCCC1",
            Profile::NumSaturatedCarbocycles {
                remove_hydrogens: false,
            },
            1,
        ),
        (
            "benzene-aliphatic",
            "c1ccccc1",
            Profile::NumAliphaticRings {
                remove_hydrogens: false,
            },
            0,
        ),
    ];
    for (label, smiles, profile, expected) in inputs {
        let input = Input::Molecular {
            case: registry::SmilesCase {
                id: label.into(),
                smiles: smiles.into(),
            },
            profile,
        };
        let record = molecular::run(&input).unwrap();
        let registry::Value::Molecular(molecular::Outcome::Unsigned(actual)) = record.output else {
            panic!("{label}: expected Unsigned outcome");
        };
        assert_eq!(actual, expected, "{label}");
        // Wrong-kind outputs stay rejected for every ring profile shape.
        assert!(
            molecular::validate_output(&profile, &molecular::Outcome::Float64Bits(0)).is_err(),
            "{label}: wrong kind must be rejected"
        );
        assert!(
            molecular::validate_output(&profile, &molecular::Outcome::Unsigned(expected)).is_ok(),
            "{label}: unsigned accepted"
        );
    }
}

#[test]
fn descriptor_query_parity_heavy_execution_uses_public_query_and_typed_unsigned() {
    use registry::molecule_plan::Profile;
    let input = Input::Molecular {
        case: registry::SmilesCase {
            id: "ethanol".into(),
            smiles: "CCO".into(),
        },
        profile: Profile::NumHeavyAtoms {
            remove_hydrogens: true,
        },
    };
    let record = molecular::run(&input).unwrap();
    let registry::Value::Molecular(molecular::Outcome::Unsigned(3)) = record.output else {
        panic!("expected Unsigned(3) heavy atoms for CCO");
    };
    // Wrong result kind is rejected by the schema validator.
    assert!(
        molecular::validate_output(
            &Profile::NumHeavyAtoms {
                remove_hydrogens: false
            },
            &molecular::Outcome::Float64Bits(0)
        )
        .is_err()
    );
    assert!(
        molecular::validate_output(
            &Profile::NumHeavyAtoms {
                remove_hydrogens: false
            },
            &molecular::Outcome::Unsigned(3)
        )
        .is_ok()
    );
    // The retained-hydrogens policy changes the CONSTRUCTOR input, not
    // the query: an explicit-H input keeps its explicit H rows.
    let retained = Input::Molecular {
        case: registry::SmilesCase {
            id: "ammonia-explicit".into(),
            smiles: "[H]N([H])[H]".into(),
        },
        profile: Profile::NumHeavyAtoms {
            remove_hydrogens: false,
        },
    };
    let record = molecular::run(&retained).unwrap();
    let registry::Value::Molecular(molecular::Outcome::Unsigned(1)) = record.output else {
        panic!("expected Unsigned(1) heavy atom for explicit-H ammonia");
    };
}

#[test]
fn descriptor_query_parity_total_registry_and_execution_use_calc_num_atoms_branch() {
    use registry::molecule_plan::{Category, Profile, TaskId};
    let task = registry::select(Some("total_atom_count")).unwrap();
    assert_eq!(task.len(), 1);
    assert_eq!(task[0].operation.name(), "total_atom_count");
    assert_eq!(TaskId::TotalAtomCount.name(), "total_atom_count");
    assert_eq!(TaskId::TotalAtomCount.category(), Category::Descriptors);
    assert_eq!(
        TaskId::TotalAtomCount.profiles(),
        vec![
            Profile::TotalAtomCount {
                remove_hydrogens: false
            },
            Profile::TotalAtomCount {
                remove_hydrogens: true
            },
        ]
    );
    // CalcNumAtoms branch: rows + attached hydrogens; methane total is 5
    // and the total is NOT the raw atom-table row count.
    let input = Input::Molecular {
        case: registry::SmilesCase {
            id: "methane".into(),
            smiles: "C".into(),
        },
        profile: Profile::TotalAtomCount {
            remove_hydrogens: true,
        },
    };
    let record = molecular::run(&input).unwrap();
    let registry::Value::Molecular(molecular::Outcome::Unsigned(5)) = record.output else {
        panic!("expected Unsigned(5) total atoms for methane");
    };
    // Wrong-kind rejection stays typed.
    assert!(
        molecular::validate_output(
            &Profile::TotalAtomCount {
                remove_hydrogens: false
            },
            &molecular::Outcome::Text("5".into())
        )
        .is_err()
    );
}

#[test]
fn descriptor_query_parity_hba_registry_and_execution_use_direct_no_branch() {
    use registry::molecule_plan::{Category, Profile, TaskId};
    let task = registry::select(Some("lipinski_hba")).unwrap();
    assert_eq!(task.len(), 1);
    assert_eq!(TaskId::LipinskiHBA.name(), "lipinski_hba");
    assert_eq!(TaskId::LipinskiHBA.category(), Category::Descriptors);
    assert_eq!(
        TaskId::LipinskiHBA.profiles(),
        vec![
            Profile::LipinskiHBA {
                remove_hydrogens: false
            },
            Profile::LipinskiHBA {
                remove_hydrogens: true
            },
        ]
    );
    // DIRECT N/O count discriminators: urea O=C(N)N has 3 and quaternary
    // [NH4+] still counts 1 (general recursive NumHBA would exclude it).
    for (smiles, expected) in [("O=C(N)N", 3_u32), ("[NH4+]", 1)] {
        let input = Input::Molecular {
            case: registry::SmilesCase {
                id: format!("case:{smiles}"),
                smiles: smiles.into(),
            },
            profile: Profile::LipinskiHBA {
                remove_hydrogens: true,
            },
        };
        let record = molecular::run(&input).unwrap();
        let registry::Value::Molecular(molecular::Outcome::Unsigned(actual)) = record.output else {
            panic!("expected Unsigned for {smiles}");
        };
        assert_eq!(actual, expected);
    }
    assert!(
        molecular::validate_output(
            &Profile::LipinskiHBA {
                remove_hydrogens: false
            },
            &molecular::Outcome::Unsigned(0)
        )
        .is_ok()
    );
}

#[test]
fn descriptor_query_parity_hbd_registry_and_execution_use_hydrogen_sum_branch() {
    use registry::molecule_plan::{Category, Profile, TaskId};
    let task = registry::select(Some("lipinski_hbd")).unwrap();
    assert_eq!(task.len(), 1);
    assert_eq!(TaskId::LipinskiHBD.name(), "lipinski_hbd");
    assert_eq!(TaskId::LipinskiHBD.category(), Category::Descriptors);
    assert_eq!(
        TaskId::LipinskiHBD.profiles(),
        vec![
            Profile::LipinskiHBD {
                remove_hydrogens: false
            },
            Profile::LipinskiHBD {
                remove_hydrogens: true
            },
        ]
    );
    // Hydrogen-SUM discriminators: ammonium 4 and deuterated water 2
    // (isotopic H neighbors count); the donor-ATOM count would differ.
    for (smiles, expected) in [("[NH4+]", 4_u32), ("[2H]O[2H]", 2)] {
        let input = Input::Molecular {
            case: registry::SmilesCase {
                id: format!("case:{smiles}"),
                smiles: smiles.into(),
            },
            profile: Profile::LipinskiHBD {
                remove_hydrogens: true,
            },
        };
        let record = molecular::run(&input).unwrap();
        let registry::Value::Molecular(molecular::Outcome::Unsigned(actual)) = record.output else {
            panic!("expected Unsigned for {smiles}");
        };
        assert_eq!(actual, expected);
    }
}

#[test]
fn descriptor_query_parity_fraction_registry_and_execution_use_exact_bits() {
    use registry::molecule_plan::{Category, Profile, TaskId};
    let task = registry::select(Some("fraction_csp3")).unwrap();
    assert_eq!(task.len(), 1);
    assert_eq!(TaskId::FractionCSP3.name(), "fraction_csp3");
    assert_eq!(TaskId::FractionCSP3.category(), Category::Descriptors);
    assert_eq!(
        TaskId::FractionCSP3.profiles(),
        vec![
            Profile::FractionCSP3 {
                remove_hydrogens: false
            },
            Profile::FractionCSP3 {
                remove_hydrogens: true
            },
        ]
    );
    // Exact Float64Bits discriminators: cyclohexane is exactly 1.0 and
    // zero-carbon ammonia (retained) is exactly +0.0; both inside the
    // finite [0,1] validation band. Out-of-band bits are rejected.
    for (smiles, remove_hydrogens, expected_bits) in [
        ("C1CCCCC1", true, 0x3ff0_0000_0000_0000_u64),
        ("[H]N([H])[H]", false, 0x0000_0000_0000_0000_u64),
    ] {
        let input = Input::Molecular {
            case: registry::SmilesCase {
                id: format!("case:{smiles}"),
                smiles: smiles.into(),
            },
            profile: Profile::FractionCSP3 { remove_hydrogens },
        };
        let record = molecular::run(&input).unwrap();
        let registry::Value::Molecular(molecular::Outcome::Float64Bits(bits)) = record.output
        else {
            panic!("expected Float64Bits for {smiles}");
        };
        assert_eq!(bits, expected_bits, "{smiles}");
    }
    assert!(
        molecular::validate_output(
            &Profile::FractionCSP3 {
                remove_hydrogens: true
            },
            &molecular::Outcome::Float64Bits(2.0_f64.to_bits())
        )
        .is_err(),
        "out-of-band CSP3 bits must be rejected"
    );
    assert!(
        molecular::validate_output(
            &Profile::FractionCSP3 {
                remove_hydrogens: true
            },
            &molecular::Outcome::Unsigned(1)
        )
        .is_err(),
        "wrong-kind Unsigned must be rejected for CSP3"
    );
}

#[test]
fn descriptor_query_global_preflight_missing_final_new_task_reference_blocks_rust_calls() {
    // A missing FINAL reference for one of the NEW query tasks stops all
    // selected Rust calls before any execution (global preflight barrier).
    let temp = tempfile::tempdir().unwrap();
    let tasks = registry::select(Some("num_heavy_atoms")).unwrap();
    let mut cases = fingerprint_corpus();
    cases.molecules.push(registry::SmilesCase {
        id: "ethanol".into(),
        smiles: "CCO".into(),
    });
    let result = prepare_then_compare(
        &tasks,
        &cases,
        temp.path(),
        |inputs| {
            if let Input::Molecular {
                profile: registry::molecule_plan::Profile::NumHeavyAtoms { .. },
                ..
            } = &inputs[inputs.len() - 1]
            {
                return Err("final new-task reference missing".into());
            }
            inputs
                .iter()
                .map(|input| {
                    let output = match input {
                        Input::Molecular { .. } => {
                            registry::Value::Molecular(molecular::Outcome::Unsigned(0))
                        }
                        Input::Fingerprint(_) => unreachable!("query-only selection"),
                    };
                    Ok(Record {
                        input: input.clone(),
                        output,
                    })
                })
                .collect()
        },
        |_| panic!("global preflight barrier bypassed for new query task"),
    );
    assert!(
        result
            .unwrap_err()
            .contains("final new-task reference missing")
    );
}

#[test]
fn descriptor_query_global_preflight_five_tasks_expand_exactly_ten_profiles() {
    // The five query tasks contribute exactly TWO profiles each (the
    // remove_hydrogens input-preparation axis, false then true): ten
    // profiles total across the five families.
    use registry::molecule_plan::TaskId;
    let mut total = 0usize;
    for id in [
        TaskId::NumHeavyAtoms,
        TaskId::TotalAtomCount,
        TaskId::LipinskiHBA,
        TaskId::LipinskiHBD,
        TaskId::FractionCSP3,
    ] {
        let profiles = id.profiles();
        assert_eq!(
            profiles.len(),
            2,
            "{} must have exactly two profiles",
            id.name()
        );
        total += profiles.len();
    }
    assert_eq!(total, 10, "exact ten-profile expansion");
    // Wrong-result-kind rejection stays enforced for every new family.
    for id in [
        TaskId::NumHeavyAtoms,
        TaskId::TotalAtomCount,
        TaskId::LipinskiHBA,
        TaskId::LipinskiHBD,
    ] {
        let profile = id.profiles()[0];
        assert!(
            molecular::validate_output(&profile, &molecular::Outcome::Text("x".into())).is_err(),
            "{} must reject Text results",
            id.name()
        );
    }
    let fraction = TaskId::FractionCSP3.profiles()[0];
    assert!(molecular::validate_output(&fraction, &molecular::Outcome::Unsigned(1)).is_err());
}

#[test]
fn molecular_corpus_preserves_blank_records_cx_whitespace_and_duplicates() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("input.smi");
    fs::write(&path, "CCO\n\nC |$label$| name\nCCO\n").unwrap();
    let cases = molecular::read_corpus(&path).unwrap();
    assert_eq!(cases.len(), 4);
    assert_eq!(cases[1].smiles, "");
    assert_eq!(cases[2].smiles, "C |$label$| name");
    assert_ne!(cases[0].id, cases[3].id);
}

#[test]
fn selected_input_family_does_not_load_an_unselected_corpus() {
    let tasks = registry::select(Some("fuzzy_and")).unwrap();
    let cases = corpus(&[], &tasks).unwrap();
    assert!(cases.molecules.is_empty());
    assert_eq!(cases.fingerprints.len(), 5000);
    assert!(
        corpus(
            &[CorpusSource {
                corpus_type: CorpusType::Smiles,
                path: "missing.smi".into()
            }],
            &tasks
        )
        .unwrap_err()
        .contains("input family")
    );
}

#[test]
fn fingerprint_matrix_is_deterministic_complete_and_valid_for_both_widths() {
    use registry::fingerprint_corpus::*;
    let pairs = generate();
    assert_eq!(PAIRS, 5000);
    assert_eq!(pairs, generate());
    assert_eq!(pairs.len(), PAIRS);
    let cases = Corpus {
        fingerprints: pairs,
        molecules: vec![],
    };
    let tasks = fingerprint_tasks().unwrap();
    registry::validate(&cases, &tasks).unwrap();
    assert_eq!(tasks.iter().map(|t| t.count(&cases)).sum::<usize>(), 20_000);
    for (shape, left_mask, right_mask) in SHAPES {
        for length in LENGTHS {
            for (counts, palette) in COUNTS {
                let prefix = format!("fp5000-v1/{shape}/{length}/{counts}/");
                let cell: Vec<_> = cases
                    .fingerprints
                    .iter()
                    .filter(|p| p.id.starts_with(&prefix))
                    .collect();
                assert_eq!(cell.len(), SAMPLES_PER_CELL);
                for pair in cell {
                    assert_eq!(pair.length, length);
                    assert_eq!(pair.left.len(), left_mask.count_ones() as usize);
                    assert_eq!(pair.right.len(), right_mask.count_ones() as usize);
                    for entries in [&pair.left, &pair.right] {
                        assert!(entries.windows(2).all(|w| w[0].0 < w[1].0));
                        assert!(
                            entries
                                .iter()
                                .all(|&(key, value)| key < length && palette.contains(&value))
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn coordinate_comparison_checks_shape_nonfinite_topology_and_declared_tolerance() {
    let sample = |x: f64, count| molecular::Outcome::Coordinates2d {
        topology: molecular::Topology {
            atoms: vec![],
            bonds: vec![],
        },
        xy_bits: vec![[x.to_bits(), 0]; count],
    };
    assert!(molecular::matches(&sample(0.0, 1), &sample(0.5e-8, 1)));
    assert!(!molecular::matches(&sample(0.0, 1), &sample(2e-8, 1)));
    assert!(!molecular::matches(&sample(0.0, 1), &sample(0.0, 0)));
    assert!(!molecular::matches(
        &sample(f64::NAN, 1),
        &sample(f64::NAN, 1)
    ));
    assert!(!molecular::matches(
        &sample(f64::INFINITY, 1),
        &sample(f64::INFINITY, 1)
    ));
}

#[test]
fn distance_matrix_profiles_and_schema_preserve_all_entries_and_float_bits() {
    use registry::molecule_plan::{Profile, TaskId};
    let profiles = TaskId::DistanceMatrix.profiles();
    assert_eq!(profiles.len(), 4);
    let profile = Profile::DistanceMatrix {
        use_bond_order: true,
        use_atom_weights: true,
    };
    assert!(profiles.contains(&profile));
    let matrix = molecular::Outcome::Matrix {
        dimension: 1,
        values_bits: vec![f64::INFINITY.to_bits()],
    };
    // Atomic number zero has an infinite weighted diagonal in the source.
    assert!(molecular::validate_output(&profile, &matrix).is_ok());
    assert!(molecular::matches(&matrix, &matrix));
    let wrong_shape = molecular::Outcome::Matrix {
        dimension: 2,
        values_bits: vec![0],
    };
    assert!(molecular::validate_output(&profile, &wrong_shape).is_err());
    let positive_zero = molecular::Outcome::Matrix {
        dimension: 1,
        values_bits: vec![0],
    };
    let negative_zero = molecular::Outcome::Matrix {
        dimension: 1,
        values_bits: vec![(-0.0f64).to_bits()],
    };
    assert!(!molecular::matches(&positive_zero, &negative_zero));
}

#[test]
fn test_stage_missing_references_never_prepares_or_executes() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    let calls = Cell::new(0);
    assert!(
        run(&tasks, &cases, temp.path(), |_| {
            calls.set(calls.get() + 1);
            Err("must not execute".into())
        })
        .unwrap_err()
        .contains("Prepare references first")
    );
    assert_eq!(calls.get(), 0);
    assert_eq!(fs::read_dir(temp.path()).unwrap().count(), 0);
}

#[test]
fn labels_are_checked_even_with_updated_reference_checksums() {
    let tasks = registry::select(Some("fuzzy_and_fingerprint_pairs")).unwrap();
    let cases = fingerprint_corpus();
    for field in ["test", "corpus_type", "case_id", "parameters"] {
        let temp = tempfile::tempdir().unwrap();
        let directory = fixture(temp.path(), tasks[0], &cases);
        let mut rows: Vec<LabeledRecord> =
            serde_json::from_slice(&read(&directory.join("reference.json")).unwrap()).unwrap();
        match field {
            "test" => rows[0].label.test = "fuzzy_or_fingerprint_pairs".into(),
            "corpus_type" => rows[0].label.corpus_type = CorpusType::Sdf,
            "case_id" => rows[0].label.case_id = "wrong-case".into(),
            "parameters" => rows[0].label.parameters = serde_json::json!({"width":"U64"}),
            _ => unreachable!(),
        }
        let reference = encode(&rows).unwrap();
        let input = read(&directory.join("input.json")).unwrap();
        fs::write(directory.join("reference.json"), &reference).unwrap();
        fs::write(
            directory.join("manifest.json"),
            encode(&identity(tasks[0], &input, &reference, rows.len())).unwrap(),
        )
        .unwrap();
        assert!(
            preflight(&tasks, &cases, temp.path())
                .err()
                .unwrap()
                .contains("label mismatch"),
            "{field}"
        );
    }
}

#[test]
fn preparation_is_sequential_and_preserves_full_parameter_order() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = fingerprint_tasks().unwrap();
    let cases = fingerprint_corpus();
    let mut names = Vec::new();
    let (_, preparation) = prepare_with(&tasks, &cases, temp.path(), |task, inputs| {
        names.push(task.key());
        assert_eq!(inputs, registry::expand(&cases, task));
        synthetic(inputs)
    })
    .unwrap();
    assert_eq!(
        names,
        ["fuzzy_and_fingerprint_pairs", "fuzzy_or_fingerprint_pairs"]
    );
    assert_eq!(preparation.rows, 24);
}

#[test]
fn registry_keys_are_unique_function_and_corpus_pairs() {
    let mut keys = std::collections::BTreeSet::new();
    for task in registry::TASKS {
        assert!(keys.insert(task.key()));
        assert_eq!(registry::select(Some(&task.key())).unwrap().len(), 1);
        assert_eq!(
            task.generator,
            format!("generate_{}", task.operation.name())
        );
    }
    let tasks = fingerprint_tasks().unwrap();
    for kind in [
        CorpusType::Smiles,
        CorpusType::Pdb,
        CorpusType::Cif,
        CorpusType::Mmcif,
        CorpusType::Sdf,
    ] {
        assert!(
            corpus(
                &[CorpusSource {
                    corpus_type: kind,
                    path: "not-a-real-file".into()
                }],
                &tasks
            )
            .unwrap_err()
            .contains("input family")
        );
    }
}
