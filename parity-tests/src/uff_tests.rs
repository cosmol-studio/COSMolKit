//! Pipeline regressions inject data; ordinary cargo tests invoke no oracle.
use super::*;
use std::cell::Cell;

fn fake_coordinate_row(conformer_id: usize) -> uff::CoordinateRow {
    let offset = conformer_id as f64;
    uff::CoordinateRow {
        conformer_id,
        xyz_bits: vec![
            [offset.to_bits(), 1.0_f64.to_bits(), 0.0_f64.to_bits()],
            [
                10.0_f64.to_bits(),
                (offset + 2.0).to_bits(),
                3.0_f64.to_bits(),
            ],
        ],
    }
}

fn fake_records(inputs: &[Input]) -> Result<Vec<Record>> {
    inputs
        .iter()
        .map(|input| {
            let Input::Uff(row) = input else {
                return Err("UFF recipe expected".into());
            };
            let mut prepared = row.clone();
            let output = match row.profile {
                uff::Profile::Coverage { .. } => uff::Observation::Coverage(true),
                uff::Profile::Optimization { .. } => {
                    prepared.preparation = Some(uff::GeometryPreparation::Ready(uff::Geometry {
                        molblock: "synthetic pipeline-only M  END".into(),
                        atom_count: 2,
                        coordinate_rows: Vec::new(),
                    }));
                    uff::Observation::Optimized {
                        status: 1,
                        energy_bits: 2.0_f64.to_bits(),
                        xyz_bits: vec![[0; 3]; 2],
                    }
                }
                uff::Profile::ConformerOptimization {
                    conformer_count, ..
                } => {
                    prepared.preparation = Some(uff::GeometryPreparation::Ready(uff::Geometry {
                        molblock: "synthetic pipeline-only M  END".into(),
                        atom_count: 2,
                        coordinate_rows: [7, 3, 11]
                            .into_iter()
                            .take(conformer_count)
                            .map(fake_coordinate_row)
                            .collect(),
                    }));
                    uff::Observation::OptimizedConformers {
                        conformers: [7, 3, 11]
                            .into_iter()
                            .take(conformer_count)
                            .map(|conformer_id| uff::ConformerObservation {
                                conformer_id,
                                status: 1,
                                energy_bits: 2.0_f64.to_bits(),
                                xyz_bits: vec![[0; 3]; 2],
                            })
                            .collect(),
                    }
                }
            };
            Ok(Record {
                input: Input::Uff(prepared),
                output: registry::Value::Uff(output),
            })
        })
        .collect()
}

fn tasks() -> [&'static Task; 2] {
    [
        registry::select(Some("uff_has_all_molecule_params")).unwrap()[0],
        registry::select(Some("uff_optimize")).unwrap()[0],
    ]
}
fn cases() -> Corpus {
    Corpus {
        molecules: vec![
            registry::SmilesCase {
                id: "one".into(),
                smiles: "CC".into(),
            },
            registry::SmilesCase {
                id: "two".into(),
                smiles: "CCO".into(),
            },
        ],
        ..Default::default()
    }
}

fn prepare_then_run_with(
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
fn uff_pipeline_complete_profiles_and_common_geometry_identity() {
    let tasks = tasks();
    let cases = cases();
    assert_eq!(tasks[0].count(&cases), 4);
    assert_eq!(tasks[1].count(&cases), 2);
    assert_eq!(
        uff::profiles(registry::Operation::UffCoverage),
        vec![
            uff::Profile::Coverage {
                add_hydrogens: false,
            },
            uff::Profile::Coverage {
                add_hydrogens: true,
            },
        ]
    );
    let profiles = uff::profiles(registry::Operation::UffOptimization);
    assert_eq!(
        profiles,
        vec![uff::Profile::Optimization {
            add_hydrogens: true,
            max_iterations: 1,
            vdw_threshold: 100,
            ignore_interfragment_interactions: true,
            conformer_id: None,
        }]
    );
    for (i, profile) in profiles.iter().enumerate() {
        assert!(!profiles[..i].contains(profile));
    }
    let temp = tempfile::tempdir().unwrap();
    let (ready, preparation) = prepare_with(&tasks, &cases, temp.path(), |_, inputs| {
        fake_records(inputs)
    })
    .unwrap();
    assert_eq!(preparation.rows, 6);
    assert_eq!(preparation.generated_tasks, 2);
    let comparisons = compare(ready, |input| {
        let Input::Uff(row) = input else {
            panic!("UFF expected")
        };
        if matches!(row.profile, uff::Profile::Optimization { .. }) {
            assert!(row.preparation.is_some());
        }
        Ok(fake_records(std::slice::from_ref(input))?.remove(0))
    });
    assert!(comparisons.iter().all(|row| row.matches));
    let (_, reused) = prepare_with(&tasks, &cases, temp.path(), |_, _| {
        panic!("valid common geometry cache regenerated")
    })
    .unwrap();
    assert_eq!(reused.generated_tasks, 0);
    assert_eq!(reused.reused_tasks, 2);
    let inputs = registry::expand(&cases, tasks[1]);
    let mut records = fake_records(&inputs).unwrap();
    let Input::Uff(row) = &mut records[0].input else {
        unreachable!()
    };
    row.case.smiles = "N".into();
    assert!(check_records(tasks[1], &inputs, &records).is_err());
}

#[test]
fn uff_pipeline_global_preflight_prevents_calls_on_missing_geometry() {
    let tasks = tasks();
    let cases = cases();
    let temp = tempfile::tempdir().unwrap();
    let result = prepare_then_run_with(
        &tasks,
        &cases,
        temp.path(),
        |inputs| {
            if matches!(
                &inputs[0],
                Input::Uff(uff::UffInput {
                    profile: uff::Profile::Optimization { .. },
                    ..
                })
            ) {
                return Err("common embedding failed".into());
            }
            fake_records(inputs)
        },
        |_| panic!("UFF operation ran before all common geometries were ready"),
    );
    assert!(result.unwrap_err().contains("0 Rust operation calls"));
    let inputs = registry::expand(&cases, tasks[1]);
    let mut records = fake_records(&inputs).unwrap();
    let Input::Uff(row) = &mut records[0].input else {
        unreachable!()
    };
    row.preparation = None;
    assert!(check_records(tasks[1], &inputs, &records).is_err());
}

#[test]
fn uff_pipeline_bitwise_observation_and_errors_are_not_false_passes() {
    let input = registry::expand(&cases(), tasks()[1]).remove(0);
    let expected = uff::Observation::Optimized {
        status: 0,
        energy_bits: 1.0_f64.to_bits(),
        xyz_bits: vec![[0; 3]],
    };
    assert!(uff::matches(&input, &expected, &expected));
    let changed = uff::Observation::Optimized {
        status: 0,
        energy_bits: 1.0_f64.to_bits() + 1,
        xyz_bits: vec![[0; 3]],
    };
    assert!(!uff::matches(&input, &expected, &changed));
    let error = uff::Observation::Error {
        stage: molecular::Stage::Operation,
        detail: "source error".into(),
        reason: None,
    };
    assert!(!uff::matches(&input, &error, &error));
}

#[test]
fn uff_pipeline_recorded_source_rejections_are_strict_error_behavior_passes() {
    let task = tasks()[1];
    let cases = cases();
    let inputs = registry::expand(&cases, task);
    let mut records = fake_records(&inputs).unwrap();
    let Input::Uff(row) = &mut records[0].input else {
        unreachable!()
    };
    row.preparation = Some(uff::GeometryPreparation::Rejected {
        stage: molecular::Stage::Preparation,
        detail: "ValueError: UFF common geometry embedding failed: one".into(),
    });
    records[0].output = registry::Value::Uff(uff::Observation::Error {
        stage: molecular::Stage::Preparation,
        detail: "ValueError: UFF common geometry embedding failed: one".into(),
        reason: None,
    });
    check_records(task, &inputs, &records).unwrap();
    let actual = uff::run(&records[0].input).unwrap();
    let registry::Value::Uff(actual_value) = actual.output else {
        unreachable!()
    };
    let registry::Value::Uff(expected_value) = &records[0].output else {
        unreachable!()
    };
    assert!(matches!(actual_value, uff::Observation::Error { .. }));
    assert!(uff::matches(
        &records[0].input,
        expected_value,
        &actual_value
    ));
    let uff::Observation::Error {
        stage,
        detail,
        reason,
    } = &actual_value
    else {
        unreachable!()
    };
    assert_eq!(*stage, molecular::Stage::Preparation);
    assert_eq!(
        detail,
        "ValueError: UFF common geometry embedding failed: one"
    );
    assert_eq!(
        *reason,
        Some(uff::ExpectedErrorReason::EmbeddingRejected {
            case_id: "one".into()
        })
    );
    for wrong in [
        uff::Observation::Error {
            stage: molecular::Stage::Operation,
            detail: detail.clone(),
            reason: reason.clone(),
        },
        uff::Observation::Error {
            stage: stage.clone(),
            detail: "another rejection".into(),
            reason: reason.clone(),
        },
        uff::Observation::Error {
            stage: stage.clone(),
            detail: detail.clone(),
            reason: Some(uff::ExpectedErrorReason::EmbeddingRejected {
                case_id: "other".into(),
            }),
        },
        uff::Observation::Error {
            stage: stage.clone(),
            detail: detail.clone(),
            reason: None,
        },
        uff::Observation::Coverage(true),
    ] {
        assert!(!uff::matches(&records[0].input, expected_value, &wrong));
    }
    let mut wrong_input = records[0].input.clone();
    let Input::Uff(wrong_row) = &mut wrong_input else {
        unreachable!()
    };
    wrong_row.preparation = None;
    assert!(!uff::matches(&wrong_input, expected_value, &actual_value));
    assert!(!uff::matches(
        &records[0].input,
        &uff::Observation::Coverage(true),
        &actual_value
    ));
    let temp = tempfile::tempdir().unwrap();
    let report = prepare_then_run_with(
        &[task],
        &cases,
        temp.path(),
        |_| Ok(records.clone()),
        |input| {
            if input == &records[0].input {
                return uff::run(input);
            }
            Ok(records
                .iter()
                .find(|record| &record.input == input)
                .unwrap()
                .clone())
        },
    )
    .unwrap();
    assert_eq!(report.len(), 2);
    assert_eq!(report.iter().filter(|row| !row.matches).count(), 0);
}

#[test]
fn uff_all_pipeline_profiles_and_shared_coordinate_rows_are_exact() {
    use std::collections::BTreeMap;

    let task = registry::select(Some("uff_optimize_conformers")).unwrap()[0];
    let cases = cases();
    let expected_profiles = vec![uff::Profile::ConformerOptimization {
        add_hydrogens: true,
        max_iterations: 1,
        vdw_threshold: 100,
        ignore_interfragment_interactions: true,
        conformer_count: 2,
    }];
    let profiles = uff::profiles(registry::Operation::UffConformerOptimization);
    assert_eq!(profiles, expected_profiles);
    assert_eq!(profiles.len(), 1);
    for (index, profile) in profiles.iter().enumerate() {
        assert!(!profiles[..index].contains(profile));
    }
    assert_eq!(task.count(&cases), 2);

    let temp = tempfile::tempdir().unwrap();
    let (ready, preparation) = prepare_with(&[task], &cases, temp.path(), |_, inputs| {
        fake_records(inputs)
    })
    .unwrap();
    assert_eq!(preparation.rows, 2);
    assert_eq!(preparation.generated_tasks, 1);
    let mut base_geometry = BTreeMap::<String, (String, usize)>::new();
    let mut rows_by_count = BTreeMap::<(String, usize), Vec<uff::CoordinateRow>>::new();
    let mut observed = 0;
    let comparisons = compare(ready, |input| {
        let Input::Uff(row) = input else {
            panic!("UFF expected")
        };
        let case_index = observed / 1;
        let profile_index = observed % 1;
        assert_eq!(&row.case, &cases.molecules[case_index]);
        assert_eq!(row.profile, expected_profiles[profile_index]);
        let uff::Profile::ConformerOptimization {
            conformer_count, ..
        } = row.profile
        else {
            panic!("all-conformer profile expected")
        };
        let Some(uff::GeometryPreparation::Ready(geometry)) = row.preparation.as_ref() else {
            panic!("prepared common geometry expected")
        };
        assert_eq!(geometry.atom_count, 2);
        assert_eq!(geometry.coordinate_rows.len(), conformer_count);
        assert_eq!(
            geometry
                .coordinate_rows
                .iter()
                .map(|coordinate_row| coordinate_row.conformer_id)
                .collect::<Vec<_>>(),
            [7, 3, 11]
                .into_iter()
                .take(conformer_count)
                .collect::<Vec<_>>()
        );
        assert!(geometry.coordinate_rows.iter().all(|coordinate_row| {
            coordinate_row.xyz_bits.len() == geometry.atom_count
                && coordinate_row
                    .xyz_bits
                    .iter()
                    .flatten()
                    .all(|bits| f64::from_bits(*bits).is_finite())
        }));
        let previous = base_geometry.insert(
            row.case.id.clone(),
            (geometry.molblock.clone(), geometry.atom_count),
        );
        if let Some(previous) = previous {
            assert_eq!(previous.0, geometry.molblock);
            assert_eq!(previous.1, geometry.atom_count);
        }
        let key = (row.case.id.clone(), conformer_count);
        let previous_rows = rows_by_count
            .entry(key)
            .or_insert_with(|| geometry.coordinate_rows.clone());
        assert_eq!(
            previous_rows.as_slice(),
            geometry.coordinate_rows.as_slice()
        );
        observed += 1;
        Ok(fake_records(std::slice::from_ref(input))?.remove(0))
    });
    assert_eq!(observed, 2);
    assert!(comparisons.iter().all(|comparison| comparison.matches));
    for case in &cases.molecules {
        let two = &rows_by_count[&(case.id.clone(), 2)];
        assert_eq!(two.len(), 2);
    }
}

#[test]
fn uff_all_pipeline_missing_last_geometry_blocks_every_selected_call() {
    let task = registry::select(Some("uff_optimize_conformers")).unwrap()[0];
    let cases = cases();
    let inputs = registry::expand(&cases, task);
    assert_eq!(inputs.len(), 2);
    let calls = Cell::new(0);
    let temp = tempfile::tempdir().unwrap();
    let error = prepare_then_run_with(
        &[task],
        &cases,
        temp.path(),
        |inputs| {
            let mut records = fake_records(inputs)?;
            let Some(last) = records.last_mut() else {
                return Err("selected task unexpectedly empty".into());
            };
            let Input::Uff(row) = &mut last.input else {
                return Err("last selected input was not UFF".into());
            };
            assert_eq!(&row.case, &cases.molecules[1]);
            assert_eq!(row.profile, uff::profiles(task.operation)[0]);
            row.preparation = None;
            Ok(records)
        },
        |_| {
            calls.set(calls.get() + 1);
            Err("executor must not run before all geometry preflights".into())
        },
    )
    .unwrap_err();
    assert!(error.contains("0 Rust operation calls"));
    assert_eq!(calls.get(), 0);
}

#[test]
fn uff_new_reference_labels_ignore_prepared_geometry_and_reject_wrong_identity() {
    let task = registry::select(Some("uff_optimize")).unwrap()[0];
    let cases = cases();
    let recipe = registry::expand(&cases, task)[0].clone();
    let record = fake_records(std::slice::from_ref(&recipe))
        .unwrap()
        .remove(0);

    task.validate_reference(&recipe, &record.input, &record.output)
        .unwrap();
    let recipe_label = reference_label(task, &recipe).unwrap();
    let prepared_label = reference_label(task, &record.input).unwrap();
    assert_eq!(recipe_label, prepared_label);
    assert_eq!(prepared_label.test, task.key());
    assert_eq!(prepared_label.corpus_type, registry::CorpusType::Smiles);
    let Input::Uff(recipe_row) = &recipe else {
        unreachable!()
    };
    assert_eq!(prepared_label.case_id, recipe_row.case.id);
    assert_eq!(
        prepared_label.parameters,
        serde_json::to_value(recipe_row.profile).unwrap()
    );

    let mut wrong_case = record.input.clone();
    let Input::Uff(row) = &mut wrong_case else {
        unreachable!()
    };
    row.case.id = "wrong-case".into();
    assert!(
        task.validate_reference(&recipe, &wrong_case, &record.output)
            .is_err()
    );

    let mut wrong_profile = record.input.clone();
    let Input::Uff(row) = &mut wrong_profile else {
        unreachable!()
    };
    let uff::Profile::Optimization {
        add_hydrogens,
        vdw_threshold,
        ignore_interfragment_interactions,
        conformer_id,
        ..
    } = row.profile
    else {
        unreachable!()
    };
    row.profile = uff::Profile::Optimization {
        add_hydrogens,
        max_iterations: 2,
        vdw_threshold,
        ignore_interfragment_interactions,
        conformer_id,
    };
    assert!(
        task.validate_reference(&recipe, &wrong_profile, &record.output)
            .is_err()
    );

    assert!(
        task.validate_reference(
            &recipe,
            &record.input,
            &registry::Value::Fingerprint(registry::FingerprintValue {
                length: 0,
                entries: Vec::new(),
            }),
        )
        .is_err()
    );
}

#[test]
fn uff_new_reference_reuses_valid_cached_preparation() {
    let task = registry::select(Some("uff_optimize")).unwrap()[0];
    let cases = cases();
    let temp = tempfile::tempdir().unwrap();
    let (ready, generated) = prepare_with(&[task], &cases, temp.path(), |_, inputs| {
        fake_records(inputs)
    })
    .unwrap();
    assert_eq!(ready.len(), 2);
    assert_eq!(generated.rows, 2);
    assert_eq!(generated.generated_tasks, 1);

    let (reused_ready, reused) = prepare_with(&[task], &cases, temp.path(), |_, _| {
        panic!("verified UFF geometry cache was regenerated")
    })
    .unwrap();
    assert_eq!(reused_ready.len(), 2);
    assert_eq!(reused.rows, 2);
    assert_eq!(reused.generated_tasks, 0);
    assert_eq!(reused.reused_tasks, 1);
}

#[test]
fn uff_new_reference_missing_last_geometry_blocks_every_executor_call() {
    let task = registry::select(Some("uff_optimize_conformers")).unwrap()[0];
    let cases = cases();
    let inputs = registry::expand(&cases, task);
    assert_eq!(inputs.len(), 2);
    let calls = Cell::new(0);
    let temp = tempfile::tempdir().unwrap();
    let error = prepare_then_run_with(
        &[task],
        &cases,
        temp.path(),
        |inputs| {
            let mut records = fake_records(inputs)?;
            let Some(last) = records.last_mut() else {
                return Err("selected task unexpectedly empty".into());
            };
            let Input::Uff(row) = &mut last.input else {
                return Err("last selected input was not UFF".into());
            };
            assert_eq!(&row.case, &cases.molecules[1]);
            assert_eq!(row.profile, uff::profiles(task.operation)[0]);
            row.preparation = None;
            Ok(records)
        },
        |_| {
            calls.set(calls.get() + 1);
            Err("executor ran before complete reference preflight".into())
        },
    )
    .unwrap_err();
    assert!(error.contains("0 Rust operation calls"));
    assert_eq!(calls.get(), 0);
}

#[test]
fn uff_original_center_parameter_errors_use_precise_source_cause() {
    // Unchanged original prepared inputs and raw reference errors from the
    // two 5000-case reports; this fixed regression never invokes an oracle.
    let records: Vec<Record> = serde_json::from_str(include_str!(
        "../../testdata/uff/expected_error_line320.json"
    ))
    .unwrap();
    assert_eq!(records.len(), 2);
    for record in records {
        let registry::Value::Uff(expected) = &record.output else {
            unreachable!()
        };
        let actual = uff::run(&record.input).unwrap();
        let registry::Value::Uff(actual) = &actual.output else {
            unreachable!()
        };
        assert!(uff::matches(&record.input, expected, actual), "{actual:?}");
        let uff::Observation::Error {
            stage,
            detail,
            reason,
        } = actual
        else {
            unreachable!()
        };
        assert_eq!(*stage, molecular::Stage::Operation);
        assert_eq!(
            *reason,
            Some(uff::ExpectedErrorReason::SourceTbpCenterParamsMissing {
                center_atom_index: 1
            })
        );
        // A different index, missing structured cause, wrong stage, generic
        // construction failure or unexpected success cannot pass.
        for wrong in [
            uff::Observation::Error {
                stage: molecular::Stage::Preparation,
                detail: detail.clone(),
                reason: reason.clone(),
            },
            uff::Observation::Error {
                stage: stage.clone(),
                detail: detail.clone(),
                reason: Some(uff::ExpectedErrorReason::SourceTbpCenterParamsMissing {
                    center_atom_index: 0,
                }),
            },
            uff::Observation::Error {
                stage: stage.clone(),
                detail: "Construction(other)".into(),
                reason: None,
            },
            uff::Observation::Error {
                stage: stage.clone(),
                detail: String::new(),
                reason: reason.clone(),
            },
            uff::Observation::Optimized {
                status: 0,
                energy_bits: 0,
                xyz_bits: vec![],
            },
        ] {
            assert!(!uff::matches(&record.input, expected, &wrong));
        }
        let uff::Observation::Error {
            stage: expected_stage,
            detail: expected_detail,
            ..
        } = expected
        else {
            unreachable!()
        };
        for wrong_detail in [
            expected_detail.replace("at2Params", "at1Params"),
            expected_detail.replace("line 79", "line 78"),
            expected_detail.replace("AngleBend.cpp", "BondStretch.cpp"),
            "RuntimeError: unexpected Construction failure".into(),
        ] {
            let wrong = uff::Observation::Error {
                stage: expected_stage.clone(),
                detail: wrong_detail,
                reason: None,
            };
            assert!(!uff::matches(&record.input, &wrong, actual));
        }
        let mut wrong_input = record.input.clone();
        let Input::Uff(row) = &mut wrong_input else {
            unreachable!()
        };
        row.case.id = "line:321".into();
        assert!(!uff::matches(&wrong_input, expected, actual));
        assert!(!uff::matches(
            &record.input,
            &uff::Observation::Coverage(true),
            actual
        ));
    }
}
