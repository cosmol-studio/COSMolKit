use super::*;
use registry::Pair;
use std::cell::Cell;

#[test]
fn svg_schema_rejects_wrong_kind_missing_payload_and_parameter_expansion() {
    use molecular::Outcome;
    use registry::molecule_plan::{Profile, TaskId};
    assert_eq!(TaskId::Svg.profiles(), vec![Profile::SvgDefault]);
    assert_eq!(
        serde_json::to_value(Profile::SvgDefault).unwrap(),
        "SvgDefault"
    );
    assert!(serde_json::from_str::<Profile>(r#"{"SvgDefault":{"width":301}}"#).is_err());
    for output in [
        Outcome::Unsigned(1),
        Outcome::Text(String::new()),
        Outcome::Text("<svg>".into()),
    ] {
        assert!(molecular::validate_output(&Profile::SvgDefault, &output).is_err());
    }
    assert!(serde_json::from_str::<Outcome>(r#"{"Text":null}"#).is_err());
    assert!(serde_json::from_str::<Outcome>(r#"{"Text":{}}"#).is_err());
    assert!(
        molecular::validate_output(&Profile::SvgDefault, &Outcome::Text("<svg></svg>".into()))
            .is_ok()
    );
    let error = Outcome::Error {
        stage: molecular::Stage::Parse,
        detail: "retained invalid input".into(),
    };
    assert!(molecular::validate_output(&Profile::SvgDefault, &error).is_ok());
    assert!(!molecular::svg_matches(&error, &error));
}

#[test]
fn svg_typed_comparison_normalizes_only_the_four_literal_branding_substitutions() {
    use molecular::Outcome;
    use registry::molecule_plan::Profile;
    let expected = "<svg xmlns:rdkit='http://www.rdkit.org/xml'><rdkit:mol/><path d='M 1.0,2.0 L 3.0,4.0'/></svg>";
    let actual = "<svg xmlns:cosmolkit='https://www.cosmol.org'><cosmolkit:mol/><path d='M 1.0,2.0 L 3.0,4.0'/></svg>";
    let cases = Corpus {
        bio_cases: vec![],
        fingerprints: vec![],
        molecules: vec![registry::SmilesCase {
            id: "svg-synthetic".into(),
            smiles: "CCO".into(),
        }],
    };
    let tasks = registry::select(Some("svg_smiles")).unwrap();
    let data = tempfile::tempdir().unwrap();
    let inputs = registry::expand(&cases, tasks[0]);
    assert_eq!(inputs.len(), 1);
    publish(
        data.path(),
        tasks[0],
        &inputs,
        &[Record {
            input: inputs[0].clone(),
            output: registry::Value::Molecular(Outcome::Text(expected.into())),
        }],
    )
    .unwrap();
    let compare = |text: &str| {
        run(&tasks, &cases, data.path(), |input| {
            Ok(Record {
                input: input.clone(),
                output: registry::Value::Molecular(Outcome::Text(text.into())),
            })
        })
        .unwrap()[0]
            .matches
    };
    assert!(compare(actual));
    for changed in [
        actual.replace("1.0", "1.1"),
        actual.replace("M ", "L "),
        actual.replace("<path", " <path"),
        actual.replace("1.0", "1.00"),
        actual.replace("https://www.cosmol.org", "https://different.invalid"),
    ] {
        assert!(
            !compare(&changed),
            "must preserve exact SVG content: {changed}"
        );
    }
    // A formula Text observation keeps its existing exact comparison.
    assert!(!molecular::matches(
        &Outcome::Text(expected.into()),
        &Outcome::Text(actual.into())
    ));
    assert_eq!(inputs[0].task_name(), "svg");
    assert!(matches!(
        &inputs[0],
        Input::Molecular {
            profile: Profile::SvgDefault,
            ..
        }
    ));
}

#[test]
fn svg_current_canonical_identity_preserves_uri_path_and_glyph_bytes() {
    use molecular::Outcome;
    let native = "<svg xmlns:rdkit='http://www.rdkit.org/xml'><rdkit:mol/><path d='M 1.0,2.0 L 3.0,4.0'/><text>N</text></svg>";
    let canonical = native
        .replace(
            "xmlns:rdkit='http://www.rdkit.org/xml'",
            "xmlns:ck='https://kit.cosmol.org/'",
        )
        .replace("rdkit:", "ck:");
    let expected = Outcome::Text(native.into());
    assert!(molecular::svg_matches(
        &expected,
        &Outcome::Text(canonical.clone())
    ));
    for changed in [
        canonical.replace("https://kit.cosmol.org/", "https://kit.cosmol.org/wrong"),
        canonical.replace("https://kit.cosmol.org/", "https://kit.cosmol.org"),
        canonical.replace("xmlns:ck", "xmlns:other"),
        canonical.replace("1.0", "1.1"),
        canonical.replace("<text>N</text>", "<text>O</text>"),
        canonical.replace("<path", " <path"),
    ] {
        assert!(
            !molecular::svg_matches(&expected, &Outcome::Text(changed.clone())),
            "must reject unapproved metadata or changed drawing bytes: {changed}"
        );
    }
    assert!(!molecular::matches(&expected, &Outcome::Text(canonical)));
    for (native, canonical) in [
        ("<text>rdkit:A</text>", "<text>ck:A</text>"),
        (
            "<text label='rdkit:A'>N</text>",
            "<text label='ck:A'>N</text>",
        ),
        ("<!-- rdkit:A -->", "<!-- ck:A -->"),
        ("<![CDATA[rdkit:A]]>", "<![CDATA[ck:A]]>"),
        (
            "<text>xmlns:rdkit='http://www.rdkit.org/xml'</text>",
            "<text>xmlns:ck='https://kit.cosmol.org/'</text>",
        ),
    ] {
        let native = format!("<svg xmlns:rdkit='http://www.rdkit.org/xml'>{native}</svg>");
        let canonical = format!("<svg xmlns:ck='https://kit.cosmol.org/'>{canonical}</svg>");
        assert!(
            !molecular::svg_matches(&Outcome::Text(native), &Outcome::Text(canonical)),
            "ordinary content must remain exact"
        );
    }
    let native = Outcome::Text(
        "<svg xmlns:rdkit='http://www.rdkit.org/xml'><rdkit:mol rdkit:numAtoms='2'/></svg>".into(),
    );
    let canonical = Outcome::Text(
        "<svg xmlns:ck='https://kit.cosmol.org/'><ck:mol ck:numAtoms='2'/></svg>".into(),
    );
    assert!(molecular::svg_matches(&native, &canonical));
    assert!(!molecular::svg_matches(
        &Outcome::Text("<svg><rdkit:mol/></svg>".into()),
        &Outcome::Text("<svg><ck:mol/></svg>".into()),
    ));
}

#[test]
fn svg_registration_preserves_complete_cargo_task_key_census() {
    let keys: Vec<_> = registry::TASKS.iter().map(|task| task.key()).collect();
    assert_eq!(
        keys,
        [
            "bio_pdb_output_pdb",
            "bio_pdb_output_cif",
            "fuzzy_and_fingerprint_pairs",
            "fuzzy_or_fingerprint_pairs",
            "smiles_read_smiles",
            "sanitize_smiles",
            "kekulize_smiles",
            "molecular_weight_smiles",
            "exact_molecular_weight_smiles",
            "molecular_formula_smiles",
            "num_heavy_atoms_smiles",
            "total_atom_count_smiles",
            "lipinski_hba_smiles",
            "lipinski_hbd_smiles",
            "fraction_csp3_smiles",
            "num_heteroatoms_smiles",
            "num_hba_smiles",
            "num_hbd_smiles",
            "add_hydrogens_smiles",
            "remove_hydrogens_smiles",
            "coordinates_2d_smiles",
            "svg_smiles",
            "distance_matrix_smiles",
            "num_rings_smiles",
            "num_heterocycles_smiles",
            "num_aromatic_rings_smiles",
            "num_saturated_rings_smiles",
            "num_aliphatic_rings_smiles",
            "num_aromatic_heterocycles_smiles",
            "num_aromatic_carbocycles_smiles",
            "num_aliphatic_heterocycles_smiles",
            "num_aliphatic_carbocycles_smiles",
            "num_saturated_heterocycles_smiles",
            "num_saturated_carbocycles_smiles",
            "uff_has_all_molecule_params_smiles",
            "uff_optimize_smiles",
            "uff_optimize_conformers_smiles",
            "morgan_fingerprint_smiles",
            "morgan_sparse_fingerprint_smiles",
            "morgan_count_fingerprint_smiles",
            "morgan_sparse_count_fingerprint_smiles",
            "chi_0_smiles",
            "chi_1_smiles",
            "hall_kier_alpha_smiles",
            "hall_kier_alpha_with_contributions_smiles",
            "kappa_1_smiles",
            "kappa_2_smiles",
            "kappa_3_smiles",
            "phi_smiles",
            "mqns_smiles",
            "chi_0_v_smiles",
            "chi_1_v_smiles",
            "chi_2_v_smiles",
            "chi_3_v_smiles",
            "chi_4_v_smiles",
            "chi_0_n_smiles",
            "chi_1_n_smiles",
            "chi_2_n_smiles",
            "chi_3_n_smiles",
            "chi_4_n_smiles",
            "chi_n_v_smiles",
            "chi_n_n_smiles",
            "substructure_match_smiles",
            "tautomer_enumeration_smiles",
            "tautomer_canonicalization_smiles",
        ]
    );
    let task = registry::select(Some("svg_smiles")).unwrap()[0];
    assert_eq!(task.generator, "generate_svg");
    assert_eq!(task.corpus_type, CorpusType::Smiles);
    let source = include_str!("../tests/reference_parity.rs");
    let registered: Vec<_> = source
        .split("tests!(")
        .nth(1)
        .unwrap()
        .split(");")
        .next()
        .unwrap()
        .split(',')
        .map(str::trim)
        .filter(|key| !key.is_empty())
        .collect();
    assert_eq!(
        keys.iter().map(String::as_str).collect::<Vec<_>>(),
        registered
    );
}

fn fingerprint_tasks() -> Result<Vec<&'static Task>> {
    Ok(registry::TASKS
        .iter()
        .filter(|t| {
            matches!(
                t.operation,
                registry::Operation::FuzzyAnd | registry::Operation::FuzzyOr
            )
        })
        .collect())
}

/// HETERO-PUBLIC synthetic registration: one new SMILES task with the two
/// remove-H profiles, the named generator, the public Unsigned adapter
/// routing, and every original task retained.
#[test]
fn heteroatoms_descriptor_task_registration_and_adapter_routing() {
    use registry::molecule_plan::{Category, Profile, TaskId};
    let task = registry::select(Some("num_heteroatoms")).unwrap();
    assert_eq!(task.len(), 1);
    assert_eq!(task[0].operation.name(), "num_heteroatoms");
    assert_eq!(task[0].generator, "generate_num_heteroatoms");
    assert!(matches!(task[0].corpus_type, registry::CorpusType::Smiles));
    assert_eq!(TaskId::NumHeteroatoms.name(), "num_heteroatoms");
    assert_eq!(TaskId::NumHeteroatoms.category(), Category::Descriptors);
    assert_eq!(
        TaskId::NumHeteroatoms.profiles(),
        vec![
            Profile::NumHeteroatoms {
                remove_hydrogens: false
            },
            Profile::NumHeteroatoms {
                remove_hydrogens: true
            },
        ]
    );
    let input = Input::Molecular {
        case: registry::SmilesCase {
            id: "cco".into(),
            smiles: "CCO".into(),
        },
        profile: Profile::NumHeteroatoms {
            remove_hydrogens: false,
        },
    };
    let record = molecular::run(&input).unwrap();
    let registry::Value::Molecular(molecular::Outcome::Unsigned(1)) = record.output else {
        panic!("expected Unsigned(1) for CCO");
    };
    let dummy = Input::Molecular {
        case: registry::SmilesCase {
            id: "dummy".into(),
            smiles: "*".into(),
        },
        profile: Profile::NumHeteroatoms {
            remove_hydrogens: false,
        },
    };
    let record = molecular::run(&dummy).unwrap();
    let registry::Value::Molecular(molecular::Outcome::Unsigned(1)) = record.output else {
        panic!("expected Unsigned(1) for dummy (Z0 counts)");
    };
    assert!(
        molecular::validate_output(
            &Profile::NumHeteroatoms {
                remove_hydrogens: false
            },
            &molecular::Outcome::Float64Bits(0)
        )
        .is_err()
    );
    assert!(
        registry::TASKS
            .iter()
            .any(|t| t.operation.name() == "fraction_csp3")
    );
    assert!(
        registry::TASKS
            .iter()
            .any(|t| t.operation.name() == "fuzzy_and")
    );
}

/// HBA-PUBLIC synthetic registration: the ONE new SMILES task with the
/// two remove-H profiles, the named generator, the public Unsigned
/// adapter routing, wrong-kind comparison rejection, the
/// general-vs-Lipinski literal distinctions on real molecules, and every
/// original task retained.
#[test]
fn hba_descriptor_task_registration_and_adapter_routing() {
    use registry::molecule_plan::{Category, Profile, TaskId};
    let task = registry::select(Some("num_hba")).unwrap();
    assert_eq!(task.len(), 1);
    assert_eq!(task[0].operation.name(), "num_hba");
    assert_eq!(task[0].generator, "generate_num_hba");
    assert!(matches!(task[0].corpus_type, registry::CorpusType::Smiles));
    assert_eq!(TaskId::NumHba.name(), "num_hba");
    assert_eq!(TaskId::NumHba.category(), Category::Descriptors);
    assert_eq!(
        TaskId::NumHba.profiles(),
        vec![
            Profile::NumHba {
                remove_hydrogens: false
            },
            Profile::NumHba {
                remove_hydrogens: true
            },
        ]
    );
    // Acid literal: general 1 vs direct Lipinski 2 on the SAME molecule.
    let acid = Input::Molecular {
        case: registry::SmilesCase {
            id: "acid".into(),
            smiles: "CC(=O)O".into(),
        },
        profile: Profile::NumHba {
            remove_hydrogens: false,
        },
    };
    let record = molecular::run(&acid).unwrap();
    let registry::Value::Molecular(molecular::Outcome::Unsigned(1)) = record.output else {
        panic!("expected Unsigned(1) for the acid general count");
    };
    let lipinski = Input::Molecular {
        case: registry::SmilesCase {
            id: "acid".into(),
            smiles: "CC(=O)O".into(),
        },
        profile: Profile::LipinskiHBA {
            remove_hydrogens: false,
        },
    };
    let record = molecular::run(&lipinski).unwrap();
    let registry::Value::Molecular(molecular::Outcome::Unsigned(2)) = record.output else {
        panic!("expected Unsigned(2) for the acid direct Lipinski count");
    };
    // Thiophene literal: general 1 vs direct Lipinski 0.
    let thiophene = Input::Molecular {
        case: registry::SmilesCase {
            id: "thiophene".into(),
            smiles: "c1ccsc1".into(),
        },
        profile: Profile::NumHba {
            remove_hydrogens: false,
        },
    };
    let record = molecular::run(&thiophene).unwrap();
    let registry::Value::Molecular(molecular::Outcome::Unsigned(1)) = record.output else {
        panic!("expected Unsigned(1) for the thiophene general count");
    };
    let lipinski_thiophene = Input::Molecular {
        case: registry::SmilesCase {
            id: "thiophene".into(),
            smiles: "c1ccsc1".into(),
        },
        profile: Profile::LipinskiHBA {
            remove_hydrogens: false,
        },
    };
    let record = molecular::run(&lipinski_thiophene).unwrap();
    let registry::Value::Molecular(molecular::Outcome::Unsigned(0)) = record.output else {
        panic!("expected Unsigned(0) for the thiophene direct Lipinski count");
    };
    // Wrong comparison kind is rejected.
    assert!(
        molecular::validate_output(
            &Profile::NumHba {
                remove_hydrogens: false
            },
            &molecular::Outcome::Float64Bits(0)
        )
        .is_err()
    );
    // Every original task is retained.
    assert!(
        registry::TASKS
            .iter()
            .any(|t| t.operation.name() == "num_heteroatoms")
    );
    assert!(
        registry::TASKS
            .iter()
            .any(|t| t.operation.name() == "lipinski_hba")
    );
}

/// HBD-PUBLIC synthetic registration: the ONE new SMILES task with the
/// two remove-H profiles, the named generator, the public Unsigned
/// adapter routing, wrong-kind comparison rejection, and every original
/// task retained.
#[test]
fn hbd_descriptor_task_registration_and_adapter_routing() {
    use registry::molecule_plan::{Category, Profile, TaskId};
    let task = registry::select(Some("num_hbd")).unwrap();
    assert_eq!(task.len(), 1);
    assert_eq!(task[0].operation.name(), "num_hbd");
    assert_eq!(task[0].generator, "generate_num_hbd");
    assert!(matches!(task[0].corpus_type, registry::CorpusType::Smiles));
    assert_eq!(TaskId::NumHbd.name(), "num_hbd");
    assert_eq!(TaskId::NumHbd.category(), Category::Descriptors);
    assert_eq!(
        TaskId::NumHbd.profiles(),
        vec![
            Profile::NumHbd {
                remove_hydrogens: false
            },
            Profile::NumHbd {
                remove_hydrogens: true
            },
        ]
    );
    // The public adapter routes ONLY through the public num_hbd query:
    // CCO has exactly one donor atom under the keep policy.
    let input = Input::Molecular {
        case: registry::SmilesCase {
            id: "cco".into(),
            smiles: "CCO".into(),
        },
        profile: Profile::NumHbd {
            remove_hydrogens: false,
        },
    };
    let record = molecular::run(&input).unwrap();
    let registry::Value::Molecular(molecular::Outcome::Unsigned(1)) = record.output else {
        panic!("expected Unsigned(1) for CCO");
    };
    // Wrong comparison kind is rejected.
    assert!(
        molecular::validate_output(
            &Profile::NumHbd {
                remove_hydrogens: false
            },
            &molecular::Outcome::Float64Bits(0)
        )
        .is_err()
    );
    // Every original task is retained.
    assert!(
        registry::TASKS
            .iter()
            .any(|t| t.operation.name() == "num_hba")
    );
    assert!(
        registry::TASKS
            .iter()
            .any(|t| t.operation.name() == "lipinski_hbd")
    );
}
fn fingerprint_corpus() -> Corpus {
    Corpus {
        fingerprints: registry::builtin(),
        molecules: vec![],
        bio_cases: vec![],
    }
}

fn bio_corpus() -> Corpus {
    Corpus {
        fingerprints: vec![],
        molecules: vec![],
        bio_cases: vec![
            registry::BioPdbCase {
                id: "bio_c02".into(),
                text: "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nTER       2      ALA A   1                                                      \nEND".into(),
                format: registry::BioPdbCorpusFormat::Pdb,
            },
            registry::BioPdbCase {
                id: "bio_c03".into(),
                text: "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nHETATM    2  C1  LIG B   2       4.000   5.000   6.000  1.00 20.00           C  \nEND".into(),
                format: registry::BioPdbCorpusFormat::Pdb,
            },
        ],
    }
}

#[cfg(test)]
mod bio_pdb_output_tests {
    use crate::registry::{
        self, BioPdbCase, BioPdbOutputProfile, CorpusType, Input, Operation, Task,
    };

    /// bio_pdb_output_pdb: task is registered with the correct key,
    /// corpus type, and generator name.
    #[test]
    fn bio_pdb_output_pdb_task_registration() {
        let tasks = registry::select(Some("bio_pdb_output_pdb")).unwrap();
        assert_eq!(tasks.len(), 1);
        let task = tasks[0];
        assert_eq!(task.operation, Operation::BioPdbOutput);
        assert_eq!(task.corpus_type, CorpusType::Pdb);
        assert_eq!(task.generator, "generate_bio_pdb_output_pdb");
        assert_eq!(task.key(), "bio_pdb_output_pdb");
    }

    /// bio_pdb_output_cif: task is registered with the correct key,
    /// corpus type, and generator name.
    #[test]
    fn bio_pdb_output_cif_task_registration() {
        let tasks = registry::select(Some("bio_pdb_output_cif")).unwrap();
        assert_eq!(tasks.len(), 1);
        let task = tasks[0];
        assert_eq!(task.operation, Operation::BioPdbOutput);
        assert_eq!(task.corpus_type, CorpusType::Cif);
        assert_eq!(task.generator, "generate_bio_pdb_output_cif");
        assert_eq!(task.key(), "bio_pdb_output_cif");
    }

    /// Executor: a simple PDB case produces nonempty output text via the
    /// public BioStructure facade.
    #[test]
    fn bio_pdb_output_executor_produces_text() {
        let input = Input::BioPdbOutput {
            case: BioPdbCase {
                id: "c02".to_string(),
                text: "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nTER       2      ALA A   1                                                      \nEND".to_string(),
                format: crate::registry::BioPdbCorpusFormat::Pdb,
            },
            profile: BioPdbOutputProfile {
                ter_records: true,
                numbered_ter: true,
                ter_ignores_type: false,
                preserve_serial: false,
                end_record: true,
            },
        };
        let record = crate::execute::run(&input).expect("executor produces output");
        match record.output {
            crate::registry::Value::BioPdbOutput(value) => {
                assert!(!value.text.is_empty(), "output text is nonempty");
                assert!(value.text.contains("ATOM"), "contains ATOM");
                assert!(value.text.contains("TER"), "contains TER");
            }
            other => panic!("unexpected Value variant: {other:?}"),
        }
    }

    /// Validate: task.validate_reference accepts nonempty output.
    #[test]
    fn bio_pdb_output_validate_accepts_nonempty() {
        let tasks = registry::select(Some("bio_pdb_output_pdb")).unwrap();
        let task = tasks[0];
        let input = Input::BioPdbOutput {
            case: BioPdbCase {
                id: "test".to_string(),
                text: "ATOM".to_string(),
                format: crate::registry::BioPdbCorpusFormat::Pdb,
            },
            profile: BioPdbOutputProfile::ALL[0],
        };
        let value = crate::registry::Value::BioPdbOutput(crate::registry::BioPdbOutputValue {
            text: "SOME OUTPUT".to_string(),
            error: None,
        });
        assert!(task.validate_reference(&input, &input, &value).is_ok());
    }

    /// All 32 profiles are present and distinct.
    #[test]
    fn bio_pdb_output_profile_all32() {
        assert_eq!(BioPdbOutputProfile::ALL.len(), 32);
        // All profiles are pairwise distinct.
        for i in 0..32 {
            for j in (i + 1)..32 {
                assert_ne!(
                    BioPdbOutputProfile::ALL[i],
                    BioPdbOutputProfile::ALL[j],
                    "profiles {i} and {j} are identical"
                );
            }
        }
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
        bio_cases: vec![],
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
    let (ready, first_preparation) = prepare_with(&tasks, &cases, temp.path(), |_, inputs| {
        generated.set(generated.get() + 1);
        synthetic(inputs)
    })
    .unwrap();
    assert_eq!(generated.get(), tasks.len());
    assert_eq!(ready.len(), 24);
    assert_eq!(
        first_preparation,
        Preparation {
            reused_tasks: 0,
            generated_tasks: 2,
            rows: 24,
        }
    );
    for task in &tasks {
        let inputs = registry::expand(&cases, task);
        let input = encode(&inputs).unwrap();
        let directory = generation(temp.path(), task, &input);
        let manifest: Manifest =
            serde_json::from_slice(&read(&directory.join("manifest.json")).unwrap()).unwrap();
        assert_eq!(manifest.schema, 3);
        let records: Vec<LabeledRecord> =
            serde_json::from_slice(&read(&directory.join("reference.json")).unwrap()).unwrap();
        assert_eq!(records.len(), inputs.len());
        for (input, row) in inputs.iter().zip(records) {
            assert_eq!(row.label, reference_label(task, input).unwrap());
            assert_eq!(row.record.input, *input);
        }
    }

    // The producer is complete before the separate read-only comparison stage.
    let report = run(&tasks, &cases, temp.path(), |input| {
        assert_eq!(generated.get(), tasks.len());
        executed.set(executed.get() + 1);
        Ok(synthetic(std::slice::from_ref(input))?.remove(0))
    })
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
    let tasks: Vec<_> = registry::select(None)
        .unwrap()
        .into_iter()
        .filter(|task| {
            matches!(
                task.operation,
                registry::Operation::FuzzyAnd
                    | registry::Operation::FuzzyOr
                    | registry::Operation::Molecular(_)
            )
        })
        .collect();
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
                        Input::Search(_) => {
                            registry::Value::Search(crate::search::Outcome::Error {
                                stage: "Parse".into(),
                                kind: "SmilesParse".into(),
                            })
                        }
                        Input::Uff(_) => {
                            return Err("UFF reference intentionally unavailable".into());
                        }
                        Input::Fingerprint(input) => {
                            registry::Value::Fingerprint(registry::FingerprintValue {
                                length: input.case.length,
                                entries: vec![],
                            })
                        }
                        Input::BioPdbOutput { .. } => {
                            registry::Value::BioPdbOutput(registry::BioPdbOutputValue {
                                text: String::new(),
                                error: None,
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
                        Input::Search(_) => {
                            registry::Value::Search(crate::search::Outcome::Error {
                                stage: "Parse".into(),
                                kind: "SmilesParse".into(),
                            })
                        }
                        Input::Molecular { .. } => {
                            registry::Value::Molecular(molecular::Outcome::Unsigned(0))
                        }
                        Input::Fingerprint(_) => unreachable!("query-only selection"),
                        Input::BioPdbOutput { .. } => {
                            registry::Value::BioPdbOutput(registry::BioPdbOutputValue {
                                text: String::new(),
                                error: None,
                            })
                        }
                        Input::Uff(_) => unreachable!("query-only selection"),
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
        bio_cases: vec![],
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
        // BioPdbOutput tasks have corpus-type-specific generators.
        if task.operation == registry::Operation::BioPdbOutput {
            assert!(task.generator.starts_with("generate_bio_pdb_output"));
        } else {
            assert_eq!(
                task.generator,
                format!("generate_{}", task.operation.name())
            );
        }
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

#[test]
fn uff_new_registration_order_profiles_and_canonical_list_rows() {
    let tasks = registry::select(None).unwrap();
    let expected = [
        (
            "bio_pdb_output_pdb",
            "generate_bio_pdb_output_pdb",
            CorpusType::Pdb,
        ),
        (
            "bio_pdb_output_cif",
            "generate_bio_pdb_output_cif",
            CorpusType::Cif,
        ),
        (
            "fuzzy_and_fingerprint_pairs",
            "generate_fuzzy_and",
            CorpusType::FingerprintPairs,
        ),
        (
            "fuzzy_or_fingerprint_pairs",
            "generate_fuzzy_or",
            CorpusType::FingerprintPairs,
        ),
        (
            "smiles_read_smiles",
            "generate_smiles_read",
            CorpusType::Smiles,
        ),
        ("sanitize_smiles", "generate_sanitize", CorpusType::Smiles),
        ("kekulize_smiles", "generate_kekulize", CorpusType::Smiles),
        (
            "molecular_weight_smiles",
            "generate_molecular_weight",
            CorpusType::Smiles,
        ),
        (
            "exact_molecular_weight_smiles",
            "generate_exact_molecular_weight",
            CorpusType::Smiles,
        ),
        (
            "molecular_formula_smiles",
            "generate_molecular_formula",
            CorpusType::Smiles,
        ),
        (
            "num_heavy_atoms_smiles",
            "generate_num_heavy_atoms",
            CorpusType::Smiles,
        ),
        (
            "total_atom_count_smiles",
            "generate_total_atom_count",
            CorpusType::Smiles,
        ),
        (
            "lipinski_hba_smiles",
            "generate_lipinski_hba",
            CorpusType::Smiles,
        ),
        (
            "lipinski_hbd_smiles",
            "generate_lipinski_hbd",
            CorpusType::Smiles,
        ),
        (
            "fraction_csp3_smiles",
            "generate_fraction_csp3",
            CorpusType::Smiles,
        ),
        (
            "num_heteroatoms_smiles",
            "generate_num_heteroatoms",
            CorpusType::Smiles,
        ),
        ("num_hba_smiles", "generate_num_hba", CorpusType::Smiles),
        ("num_hbd_smiles", "generate_num_hbd", CorpusType::Smiles),
        (
            "add_hydrogens_smiles",
            "generate_add_hydrogens",
            CorpusType::Smiles,
        ),
        (
            "remove_hydrogens_smiles",
            "generate_remove_hydrogens",
            CorpusType::Smiles,
        ),
        (
            "coordinates_2d_smiles",
            "generate_coordinates_2d",
            CorpusType::Smiles,
        ),
        ("svg_smiles", "generate_svg", CorpusType::Smiles),
        (
            "distance_matrix_smiles",
            "generate_distance_matrix",
            CorpusType::Smiles,
        ),
        ("num_rings_smiles", "generate_num_rings", CorpusType::Smiles),
        (
            "num_heterocycles_smiles",
            "generate_num_heterocycles",
            CorpusType::Smiles,
        ),
        (
            "num_aromatic_rings_smiles",
            "generate_num_aromatic_rings",
            CorpusType::Smiles,
        ),
        (
            "num_saturated_rings_smiles",
            "generate_num_saturated_rings",
            CorpusType::Smiles,
        ),
        (
            "num_aliphatic_rings_smiles",
            "generate_num_aliphatic_rings",
            CorpusType::Smiles,
        ),
        (
            "num_aromatic_heterocycles_smiles",
            "generate_num_aromatic_heterocycles",
            CorpusType::Smiles,
        ),
        (
            "num_aromatic_carbocycles_smiles",
            "generate_num_aromatic_carbocycles",
            CorpusType::Smiles,
        ),
        (
            "num_aliphatic_heterocycles_smiles",
            "generate_num_aliphatic_heterocycles",
            CorpusType::Smiles,
        ),
        (
            "num_aliphatic_carbocycles_smiles",
            "generate_num_aliphatic_carbocycles",
            CorpusType::Smiles,
        ),
        (
            "num_saturated_heterocycles_smiles",
            "generate_num_saturated_heterocycles",
            CorpusType::Smiles,
        ),
        (
            "num_saturated_carbocycles_smiles",
            "generate_num_saturated_carbocycles",
            CorpusType::Smiles,
        ),
        (
            "uff_has_all_molecule_params_smiles",
            "generate_uff_has_all_molecule_params",
            CorpusType::Smiles,
        ),
        (
            "uff_optimize_smiles",
            "generate_uff_optimize",
            CorpusType::Smiles,
        ),
        (
            "uff_optimize_conformers_smiles",
            "generate_uff_optimize_conformers",
            CorpusType::Smiles,
        ),
        (
            "morgan_fingerprint_smiles",
            "generate_morgan_fingerprint",
            CorpusType::Smiles,
        ),
        (
            "morgan_sparse_fingerprint_smiles",
            "generate_morgan_sparse_fingerprint",
            CorpusType::Smiles,
        ),
        (
            "morgan_count_fingerprint_smiles",
            "generate_morgan_count_fingerprint",
            CorpusType::Smiles,
        ),
        (
            "morgan_sparse_count_fingerprint_smiles",
            "generate_morgan_sparse_count_fingerprint",
            CorpusType::Smiles,
        ),
        ("chi_0_smiles", "generate_chi_0", CorpusType::Smiles),
        ("chi_1_smiles", "generate_chi_1", CorpusType::Smiles),
        (
            "hall_kier_alpha_smiles",
            "generate_hall_kier_alpha",
            CorpusType::Smiles,
        ),
        (
            "hall_kier_alpha_with_contributions_smiles",
            "generate_hall_kier_alpha_with_contributions",
            CorpusType::Smiles,
        ),
        ("kappa_1_smiles", "generate_kappa_1", CorpusType::Smiles),
        ("kappa_2_smiles", "generate_kappa_2", CorpusType::Smiles),
        ("kappa_3_smiles", "generate_kappa_3", CorpusType::Smiles),
        ("phi_smiles", "generate_phi", CorpusType::Smiles),
        ("mqns_smiles", "generate_mqns", CorpusType::Smiles),
        ("chi_0_v_smiles", "generate_chi_0_v", CorpusType::Smiles),
        ("chi_1_v_smiles", "generate_chi_1_v", CorpusType::Smiles),
        ("chi_2_v_smiles", "generate_chi_2_v", CorpusType::Smiles),
        ("chi_3_v_smiles", "generate_chi_3_v", CorpusType::Smiles),
        ("chi_4_v_smiles", "generate_chi_4_v", CorpusType::Smiles),
        ("chi_0_n_smiles", "generate_chi_0_n", CorpusType::Smiles),
        ("chi_1_n_smiles", "generate_chi_1_n", CorpusType::Smiles),
        ("chi_2_n_smiles", "generate_chi_2_n", CorpusType::Smiles),
        ("chi_3_n_smiles", "generate_chi_3_n", CorpusType::Smiles),
        ("chi_4_n_smiles", "generate_chi_4_n", CorpusType::Smiles),
        ("chi_n_v_smiles", "generate_chi_n_v", CorpusType::Smiles),
        ("chi_n_n_smiles", "generate_chi_n_n", CorpusType::Smiles),
        (
            "substructure_match_smiles",
            "generate_substructure_match",
            CorpusType::Smiles,
        ),
        (
            "tautomer_enumeration_smiles",
            "generate_tautomer_enumeration",
            CorpusType::Smiles,
        ),
        (
            "tautomer_canonicalization_smiles",
            "generate_tautomer_canonicalization",
            CorpusType::Smiles,
        ),
    ];
    assert_eq!(tasks.len(), expected.len());
    assert_eq!(
        tasks
            .iter()
            .map(|task| (task.key(), task.generator, task.corpus_type))
            .collect::<Vec<_>>(),
        expected
            .iter()
            .map(|(key, generator, corpus_type)| (key.to_string(), *generator, *corpus_type))
            .collect::<Vec<_>>()
    );

    let uff_operations = [
        (registry::Operation::UffCoverage, 2),
        (registry::Operation::UffOptimization, 1),
        (registry::Operation::UffConformerOptimization, 1),
    ];
    for (operation, expected_profiles) in uff_operations {
        assert_eq!(uff::profiles(operation).len(), expected_profiles);
    }

    // The canonical list consists only of registry identity and profile counts;
    // this test never invokes the CK executor or the Python oracle.
    let cases = Corpus {
        bio_cases: vec![],
        fingerprints: registry::builtin(),
        molecules: vec![registry::SmilesCase {
            id: "list-row".into(),
            smiles: "C".into(),
        }],
    };
    let list_rows = tasks
        .iter()
        .map(|task| {
            format!(
                "{}: {} cases; {}",
                task.key(),
                task.count(&cases),
                task.generator
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        list_rows,
        [
            "bio_pdb_output_pdb: 0 cases; generate_bio_pdb_output_pdb",
            "bio_pdb_output_cif: 0 cases; generate_bio_pdb_output_cif",
            "fuzzy_and_fingerprint_pairs: 12 cases; generate_fuzzy_and",
            "fuzzy_or_fingerprint_pairs: 12 cases; generate_fuzzy_or",
            "smiles_read_smiles: 4 cases; generate_smiles_read",
            "sanitize_smiles: 1 cases; generate_sanitize",
            "kekulize_smiles: 2 cases; generate_kekulize",
            "molecular_weight_smiles: 2 cases; generate_molecular_weight",
            "exact_molecular_weight_smiles: 2 cases; generate_exact_molecular_weight",
            "molecular_formula_smiles: 4 cases; generate_molecular_formula",
            "num_heavy_atoms_smiles: 2 cases; generate_num_heavy_atoms",
            "total_atom_count_smiles: 2 cases; generate_total_atom_count",
            "lipinski_hba_smiles: 2 cases; generate_lipinski_hba",
            "lipinski_hbd_smiles: 2 cases; generate_lipinski_hbd",
            "fraction_csp3_smiles: 2 cases; generate_fraction_csp3",
            "num_heteroatoms_smiles: 2 cases; generate_num_heteroatoms",
            "num_hba_smiles: 2 cases; generate_num_hba",
            "num_hbd_smiles: 2 cases; generate_num_hbd",
            "add_hydrogens_smiles: 2 cases; generate_add_hydrogens",
            "remove_hydrogens_smiles: 2 cases; generate_remove_hydrogens",
            "coordinates_2d_smiles: 1 cases; generate_coordinates_2d",
            "svg_smiles: 1 cases; generate_svg",
            "distance_matrix_smiles: 4 cases; generate_distance_matrix",
            "num_rings_smiles: 2 cases; generate_num_rings",
            "num_heterocycles_smiles: 2 cases; generate_num_heterocycles",
            "num_aromatic_rings_smiles: 2 cases; generate_num_aromatic_rings",
            "num_saturated_rings_smiles: 2 cases; generate_num_saturated_rings",
            "num_aliphatic_rings_smiles: 2 cases; generate_num_aliphatic_rings",
            "num_aromatic_heterocycles_smiles: 2 cases; generate_num_aromatic_heterocycles",
            "num_aromatic_carbocycles_smiles: 2 cases; generate_num_aromatic_carbocycles",
            "num_aliphatic_heterocycles_smiles: 2 cases; generate_num_aliphatic_heterocycles",
            "num_aliphatic_carbocycles_smiles: 2 cases; generate_num_aliphatic_carbocycles",
            "num_saturated_heterocycles_smiles: 2 cases; generate_num_saturated_heterocycles",
            "num_saturated_carbocycles_smiles: 2 cases; generate_num_saturated_carbocycles",
            "uff_has_all_molecule_params_smiles: 2 cases; generate_uff_has_all_molecule_params",
            "uff_optimize_smiles: 1 cases; generate_uff_optimize",
            "uff_optimize_conformers_smiles: 1 cases; generate_uff_optimize_conformers",
            "morgan_fingerprint_smiles: 16 cases; generate_morgan_fingerprint",
            "morgan_sparse_fingerprint_smiles: 16 cases; generate_morgan_sparse_fingerprint",
            "morgan_count_fingerprint_smiles: 16 cases; generate_morgan_count_fingerprint",
            "morgan_sparse_count_fingerprint_smiles: 16 cases; generate_morgan_sparse_count_fingerprint",
            "chi_0_smiles: 1 cases; generate_chi_0",
            "chi_1_smiles: 1 cases; generate_chi_1",
            "hall_kier_alpha_smiles: 1 cases; generate_hall_kier_alpha",
            "hall_kier_alpha_with_contributions_smiles: 1 cases; generate_hall_kier_alpha_with_contributions",
            "kappa_1_smiles: 1 cases; generate_kappa_1",
            "kappa_2_smiles: 1 cases; generate_kappa_2",
            "kappa_3_smiles: 1 cases; generate_kappa_3",
            "phi_smiles: 1 cases; generate_phi",
            "mqns_smiles: 2 cases; generate_mqns",
            "chi_0_v_smiles: 1 cases; generate_chi_0_v",
            "chi_1_v_smiles: 1 cases; generate_chi_1_v",
            "chi_2_v_smiles: 1 cases; generate_chi_2_v",
            "chi_3_v_smiles: 1 cases; generate_chi_3_v",
            "chi_4_v_smiles: 1 cases; generate_chi_4_v",
            "chi_0_n_smiles: 1 cases; generate_chi_0_n",
            "chi_1_n_smiles: 1 cases; generate_chi_1_n",
            "chi_2_n_smiles: 1 cases; generate_chi_2_n",
            "chi_3_n_smiles: 1 cases; generate_chi_3_n",
            "chi_4_n_smiles: 1 cases; generate_chi_4_n",
            "chi_n_v_smiles: 7 cases; generate_chi_n_v",
            "chi_n_n_smiles: 7 cases; generate_chi_n_n",
            "substructure_match_smiles: 18 cases; generate_substructure_match",
            "tautomer_enumeration_smiles: 2 cases; generate_tautomer_enumeration",
            "tautomer_canonicalization_smiles: 2 cases; generate_tautomer_canonicalization",
        ]
        .map(str::to_string)
    );
}

fn expected_morgan_profiles(
    output: registry::molecule_plan::MorganOutputKind,
) -> Vec<registry::molecule_plan::Profile> {
    use registry::molecule_plan::{MorganInvariantKind, Profile};

    let mut profiles = Vec::with_capacity(16);
    for radius in [2, 3] {
        for include_chirality in [false, true] {
            for invariants in [
                MorganInvariantKind::Connectivity,
                MorganInvariantKind::Features,
            ] {
                for count_simulation in [false, true] {
                    profiles.push(Profile::Morgan {
                        output,
                        radius,
                        include_chirality,
                        invariants,
                        count_simulation,
                    });
                }
            }
        }
    }
    profiles
}

#[test]
fn parity_morgan_registry_profiles_names_counts_defaults_and_selection() {
    use registry::molecule_plan::{
        Comparison, InputState, MorganOutputKind, Prerequisite, Profile, TaskId,
    };

    let families = [
        (
            TaskId::MorganFingerprint,
            MorganOutputKind::DenseBits,
            "morgan_fingerprint",
        ),
        (
            TaskId::MorganSparseFingerprint,
            MorganOutputKind::SparseBits,
            "morgan_sparse_fingerprint",
        ),
        (
            TaskId::MorganCountFingerprint,
            MorganOutputKind::HashedCounts,
            "morgan_count_fingerprint",
        ),
        (
            TaskId::MorganSparseCountFingerprint,
            MorganOutputKind::SparseCounts,
            "morgan_sparse_count_fingerprint",
        ),
    ];
    let case = registry::SmilesCase {
        id: "registry-ethanol".into(),
        smiles: "CCO".into(),
    };
    let cases = Corpus {
        bio_cases: vec![],
        fingerprints: vec![],
        molecules: vec![case.clone()],
    };
    let mut all_profiles = Vec::with_capacity(64);

    for (id, output, task_name) in families {
        assert_eq!(id.name(), task_name);
        assert_eq!(
            id.category(),
            registry::molecule_plan::Category::FingerprintGeneration
        );

        let matrix = registry::molecule_plan::TASKS
            .iter()
            .find(|task| task.id == id)
            .unwrap();
        assert_eq!(matrix.input, InputState::SanitizedHydrogensRemoved);
        assert_eq!(
            matrix.comparison,
            Comparison::MorganFingerprintAndAdditionalOutput
        );
        assert_eq!(
            matrix.prerequisite,
            Prerequisite::PublicValenceReadoutAndMolecularPipeline
        );

        let profiles = id.profiles();
        assert_eq!(profiles, expected_morgan_profiles(output));
        assert_eq!(profiles.len(), 16);
        let default_profile = Profile::Morgan {
            output,
            radius: 3,
            include_chirality: false,
            invariants: registry::molecule_plan::MorganInvariantKind::Connectivity,
            count_simulation: false,
        };
        assert_eq!(
            profiles
                .iter()
                .filter(|profile| **profile == default_profile)
                .count(),
            1
        );
        assert!(profiles.iter().all(|profile| matches!(
            profile,
            Profile::Morgan { output: profile_output, .. } if *profile_output == output
        )));

        let selected = registry::select(Some(task_name)).unwrap();
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].operation, registry::Operation::Molecular(id));
        assert_eq!(selected[0].count(&cases), 16);
        let expanded = registry::expand(&cases, selected[0]);
        assert_eq!(expanded.len(), 16);
        for input in &expanded {
            assert_eq!(input.task_name(), task_name);
            let Input::Molecular {
                case: expanded_case,
                profile,
            } = input
            else {
                panic!("Morgan registry expansion must be molecular")
            };
            assert_eq!(expanded_case, &case);
            assert!(profiles.contains(profile));
        }
        all_profiles.extend(profiles);
    }

    assert_eq!(all_profiles.len(), 64);
    for (index, profile) in all_profiles.iter().enumerate() {
        assert!(!all_profiles[..index].contains(profile));
    }

    let defaults = cosmolkit::MorganFingerprintParams::default();
    assert!(defaults.from_atoms.is_none());
    assert!(defaults.ignore_atoms.is_none());
    assert!(defaults.custom_atom_invariants.is_none());
    assert!(defaults.custom_bond_invariants.is_none());
    assert_eq!(defaults.conformer_id, -1);
    assert!(matches!(
        defaults.invariants,
        cosmolkit::MorganInvariants::Connectivity
    ));
    assert_eq!(defaults.generator, cosmolkit::MorganParams::default());
    assert_eq!(defaults.generator.radius, 3);
    assert!(!defaults.generator.include_chirality);
    assert!(defaults.generator.use_bond_types);
    assert!(defaults.generator.include_ring_membership);
    assert!(!defaults.generator.only_nonzero_invariants);
    assert!(!defaults.generator.include_redundant_environments);
    assert_eq!(defaults.generator.fp_size, 2048);
    assert!(!defaults.generator.count_simulation);
    assert_eq!(defaults.generator.count_bounds, [1, 2, 4, 8]);
    assert_eq!(defaults.generator.bits_per_feature, 1);

    let selected_all = registry::select(None).unwrap();
    assert_eq!(selected_all.len(), registry::TASKS.len());
    assert_eq!(selected_all.len(), 65);
    assert_eq!(
        selected_all[62].operation,
        crate::registry::Operation::SubstructureMatch
    );
    assert!(
        selected_all
            .iter()
            .zip(registry::TASKS)
            .all(|(selected, registered)| std::ptr::eq(*selected, registered))
    );
    assert_eq!(
        selected_all
            .iter()
            .map(|task| task.operation.name())
            .collect::<Vec<_>>(),
        [
            "bio_pdb_output",
            "bio_pdb_output",
            "fuzzy_and",
            "fuzzy_or",
            "smiles_read",
            "sanitize",
            "kekulize",
            "molecular_weight",
            "exact_molecular_weight",
            "molecular_formula",
            "num_heavy_atoms",
            "total_atom_count",
            "lipinski_hba",
            "lipinski_hbd",
            "fraction_csp3",
            "num_heteroatoms",
            "num_hba",
            "num_hbd",
            "add_hydrogens",
            "remove_hydrogens",
            "coordinates_2d",
            "svg",
            "distance_matrix",
            "num_rings",
            "num_heterocycles",
            "num_aromatic_rings",
            "num_saturated_rings",
            "num_aliphatic_rings",
            "num_aromatic_heterocycles",
            "num_aromatic_carbocycles",
            "num_aliphatic_heterocycles",
            "num_aliphatic_carbocycles",
            "num_saturated_heterocycles",
            "num_saturated_carbocycles",
            "uff_has_all_molecule_params",
            "uff_optimize",
            "uff_optimize_conformers",
            "morgan_fingerprint",
            "morgan_sparse_fingerprint",
            "morgan_count_fingerprint",
            "morgan_sparse_count_fingerprint",
            "chi_0",
            "chi_1",
            "hall_kier_alpha",
            "hall_kier_alpha_with_contributions",
            "kappa_1",
            "kappa_2",
            "kappa_3",
            "phi",
            "mqns",
            "chi_0_v",
            "chi_1_v",
            "chi_2_v",
            "chi_3_v",
            "chi_4_v",
            "chi_0_n",
            "chi_1_n",
            "chi_2_n",
            "chi_3_n",
            "chi_4_n",
            "chi_n_v",
            "chi_n_n",
            "substructure_match",
            "tautomer_enumeration",
            "tautomer_canonicalization",
        ]
    );
    assert_eq!(
        selected_all
            .iter()
            .map(|task| task.key())
            .collect::<Vec<_>>(),
        [
            "bio_pdb_output_pdb",
            "bio_pdb_output_cif",
            "fuzzy_and_fingerprint_pairs",
            "fuzzy_or_fingerprint_pairs",
            "smiles_read_smiles",
            "sanitize_smiles",
            "kekulize_smiles",
            "molecular_weight_smiles",
            "exact_molecular_weight_smiles",
            "molecular_formula_smiles",
            "num_heavy_atoms_smiles",
            "total_atom_count_smiles",
            "lipinski_hba_smiles",
            "lipinski_hbd_smiles",
            "fraction_csp3_smiles",
            "num_heteroatoms_smiles",
            "num_hba_smiles",
            "num_hbd_smiles",
            "add_hydrogens_smiles",
            "remove_hydrogens_smiles",
            "coordinates_2d_smiles",
            "svg_smiles",
            "distance_matrix_smiles",
            "num_rings_smiles",
            "num_heterocycles_smiles",
            "num_aromatic_rings_smiles",
            "num_saturated_rings_smiles",
            "num_aliphatic_rings_smiles",
            "num_aromatic_heterocycles_smiles",
            "num_aromatic_carbocycles_smiles",
            "num_aliphatic_heterocycles_smiles",
            "num_aliphatic_carbocycles_smiles",
            "num_saturated_heterocycles_smiles",
            "num_saturated_carbocycles_smiles",
            "uff_has_all_molecule_params_smiles",
            "uff_optimize_smiles",
            "uff_optimize_conformers_smiles",
            "morgan_fingerprint_smiles",
            "morgan_sparse_fingerprint_smiles",
            "morgan_count_fingerprint_smiles",
            "morgan_sparse_count_fingerprint_smiles",
            "chi_0_smiles",
            "chi_1_smiles",
            "hall_kier_alpha_smiles",
            "hall_kier_alpha_with_contributions_smiles",
            "kappa_1_smiles",
            "kappa_2_smiles",
            "kappa_3_smiles",
            "phi_smiles",
            "mqns_smiles",
            "chi_0_v_smiles",
            "chi_1_v_smiles",
            "chi_2_v_smiles",
            "chi_3_v_smiles",
            "chi_4_v_smiles",
            "chi_0_n_smiles",
            "chi_1_n_smiles",
            "chi_2_n_smiles",
            "chi_3_n_smiles",
            "chi_4_n_smiles",
            "chi_n_v_smiles",
            "chi_n_n_smiles",
            "substructure_match_smiles",
            "tautomer_enumeration_smiles",
            "tautomer_canonicalization_smiles",
        ]
    );

    assert_eq!(registry::molecule_plan::TASKS.len(), 60);
    assert_eq!(
        registry::molecule_plan::TASKS
            .iter()
            .map(|task| task.id.profiles().len())
            .collect::<Vec<_>>(),
        [
            4, 4, 1, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4,
            8, 1, 1, 2, 16, 16, 16, 16, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7,
            7, 2, 2
        ]
    );
}

#[test]
fn parity_morgan_registry_smiles_tasks_bind_literal_generators_and_function_keys() {
    use registry::molecule_plan::TaskId;

    let cases = [
        (
            TaskId::MorganFingerprint,
            "morgan_fingerprint",
            "morgan_fingerprint_smiles",
            "generate_morgan_fingerprint",
        ),
        (
            TaskId::MorganSparseFingerprint,
            "morgan_sparse_fingerprint",
            "morgan_sparse_fingerprint_smiles",
            "generate_morgan_sparse_fingerprint",
        ),
        (
            TaskId::MorganCountFingerprint,
            "morgan_count_fingerprint",
            "morgan_count_fingerprint_smiles",
            "generate_morgan_count_fingerprint",
        ),
        (
            TaskId::MorganSparseCountFingerprint,
            "morgan_sparse_count_fingerprint",
            "morgan_sparse_count_fingerprint_smiles",
            "generate_morgan_sparse_count_fingerprint",
        ),
    ];

    for (id, function, key, generator) in cases {
        let task = registry::select(Some(function)).unwrap();
        assert_eq!(task.len(), 1);
        assert_eq!(task[0].operation, registry::Operation::Molecular(id));
        assert_eq!(task[0].operation.name(), function);
        assert_eq!(task[0].corpus_type, CorpusType::Smiles);
        assert_eq!(task[0].key(), key);
        assert_eq!(task[0].key(), format!("{function}_smiles"));
        assert_eq!(task[0].generator, generator);
    }
}

#[test]
fn parity_morgan_schema_output_tags_lengths_and_signed_indices() {
    use molecular::{
        MorganAdditionalOutput, MorganDenseBitsOutput, MorganHashedCountsOutput,
        MorganSparseBitsOutput, MorganSparseCountsOutput, Outcome,
    };
    use registry::molecule_plan::{MorganInvariantKind, MorganOutputKind, Profile};

    let profile = |output| Profile::Morgan {
        output,
        radius: 2,
        include_chirality: false,
        invariants: MorganInvariantKind::Connectivity,
        count_simulation: false,
    };
    let empty_output = || MorganAdditionalOutput {
        atom_counts: None,
        atom_to_bits: None,
        bit_info_map: None,
        bit_paths: None,
        atoms_per_bit: None,
    };

    let dense = Outcome::MorganDenseBits(MorganDenseBitsOutput {
        length: 2048,
        on_bits: vec![0, 2047],
        additional_output: empty_output(),
    });
    assert!(molecular::validate_output(&profile(MorganOutputKind::DenseBits), &dense).is_ok());
    assert!(molecular::validate_output(&profile(MorganOutputKind::SparseBits), &dense).is_err());
    assert!(
        molecular::validate_output(
            &profile(MorganOutputKind::DenseBits),
            &Outcome::MorganDenseBits(MorganDenseBitsOutput {
                length: 2047,
                on_bits: vec![0, 2046],
                additional_output: empty_output(),
            })
        )
        .is_err()
    );

    let sparse_bits = Outcome::MorganSparseBits(MorganSparseBitsOutput {
        length: u32::MAX,
        on_bits: vec![i32::MIN, -1, 0, i32::MAX],
        additional_output: empty_output(),
    });
    assert!(
        molecular::validate_output(&profile(MorganOutputKind::SparseBits), &sparse_bits).is_ok()
    );
    assert!(
        molecular::validate_output(
            &profile(MorganOutputKind::SparseBits),
            &Outcome::MorganSparseBits(MorganSparseBitsOutput {
                length: u32::MAX - 1,
                on_bits: vec![],
                additional_output: empty_output(),
            })
        )
        .is_err()
    );
    let sparse_bits_json = serde_json::to_value(&sparse_bits).unwrap();
    assert_eq!(
        sparse_bits_json["MorganSparseBits"]["on_bits"],
        serde_json::json!([i32::MIN, -1, 0, i32::MAX])
    );

    let hashed_counts = Outcome::MorganHashedCounts(MorganHashedCountsOutput {
        length: 2048,
        entries: vec![(0, 2), (2047, 1)],
        additional_output: empty_output(),
    });
    assert!(
        molecular::validate_output(&profile(MorganOutputKind::HashedCounts), &hashed_counts)
            .is_ok()
    );
    assert!(
        molecular::validate_output(
            &profile(MorganOutputKind::HashedCounts),
            &Outcome::MorganHashedCounts(MorganHashedCountsOutput {
                length: 2048,
                entries: vec![(2048, 1)],
                additional_output: empty_output(),
            })
        )
        .is_err()
    );

    let sparse_counts = Outcome::MorganSparseCounts(MorganSparseCountsOutput {
        length: u64::MAX,
        entries: vec![(0, -1), (u64::MAX, 1)],
        additional_output: empty_output(),
    });
    assert!(
        molecular::validate_output(&profile(MorganOutputKind::SparseCounts), &sparse_counts)
            .is_ok()
    );
}

#[test]
fn parity_morgan_schema_requires_every_additional_output_field_and_keeps_options() {
    use molecular::MorganAdditionalOutput;

    let wire = serde_json::to_value(MorganAdditionalOutput {
        atom_counts: None,
        atom_to_bits: Some(vec![]),
        bit_info_map: None,
        bit_paths: None,
        atoms_per_bit: None,
    })
    .unwrap();
    assert_eq!(wire["atom_counts"], serde_json::Value::Null);
    assert_eq!(wire["atom_to_bits"], serde_json::json!([]));
    let decoded: MorganAdditionalOutput = serde_json::from_value(wire.clone()).unwrap();
    assert_eq!(decoded.atom_counts, None);
    assert_eq!(decoded.atom_to_bits, Some(vec![]));
    assert_eq!(decoded, serde_json::from_value(wire.clone()).unwrap());

    for field in [
        "atom_counts",
        "atom_to_bits",
        "bit_info_map",
        "bit_paths",
        "atoms_per_bit",
    ] {
        let mut omitted = wire.clone();
        omitted.as_object_mut().unwrap().remove(field);
        assert!(
            serde_json::from_value::<MorganAdditionalOutput>(omitted).is_err(),
            "missing FingerprintAdditionalOutput field {field} must be rejected"
        );
    }

    let mut malformed = wire;
    malformed["atom_counts"] = serde_json::json!([4294967296u64]);
    assert!(serde_json::from_value::<MorganAdditionalOutput>(malformed).is_err());
}

#[test]
fn parity_morgan_schema_preserves_nested_order_and_rejects_malformed_metadata() {
    use molecular::{MorganAdditionalOutput, MorganDenseBitsOutput, Outcome};
    use registry::molecule_plan::{MorganInvariantKind, MorganOutputKind, Profile};

    let additional_output = MorganAdditionalOutput {
        atom_counts: Some(vec![3, 0]),
        atom_to_bits: Some(vec![vec![9, 9, 4], vec![]]),
        bit_info_map: Some(vec![(4, vec![(1, 2), (1, 2), (0, 1)]), (9, vec![(2, 3)])]),
        bit_paths: Some(vec![(4, vec![vec![2, 1, 1], vec![2, 1, 1]]), (9, vec![])]),
        atoms_per_bit: Some(vec![(4, vec![vec![2, 2, 0, 0]]), (9, vec![])]),
    };
    let outcome = Outcome::MorganDenseBits(MorganDenseBitsOutput {
        length: 2048,
        on_bits: vec![4, 9],
        additional_output: additional_output.clone(),
    });
    let profile = Profile::Morgan {
        output: MorganOutputKind::DenseBits,
        radius: 2,
        include_chirality: false,
        invariants: MorganInvariantKind::Connectivity,
        count_simulation: false,
    };
    assert!(molecular::validate_output(&profile, &outcome).is_ok());
    let roundtrip: Outcome =
        serde_json::from_value(serde_json::to_value(&outcome).unwrap()).unwrap();
    assert_eq!(roundtrip, outcome);

    let mut unsorted_map = additional_output.clone();
    unsorted_map.bit_info_map = Some(vec![(9, vec![]), (4, vec![])]);
    assert!(
        molecular::validate_output(
            &profile,
            &Outcome::MorganDenseBits(MorganDenseBitsOutput {
                length: 2048,
                on_bits: vec![4, 9],
                additional_output: unsorted_map,
            })
        )
        .is_err()
    );

    let mut mismatched_per_atom = additional_output;
    mismatched_per_atom.atom_to_bits = Some(vec![vec![]]);
    assert!(
        molecular::validate_output(
            &profile,
            &Outcome::MorganDenseBits(MorganDenseBitsOutput {
                length: 2048,
                on_bits: vec![4, 9],
                additional_output: mismatched_per_atom,
            })
        )
        .is_err()
    );
}

#[test]
fn parity_morgan_executor_dispatches_all_public_families_and_channels() {
    use molecular::{
        MorganAdditionalOutput, MorganDenseBitsOutput, MorganHashedCountsOutput,
        MorganSparseBitsOutput, MorganSparseCountsOutput, Outcome,
    };
    use registry::molecule_plan::{MorganInvariantKind, MorganOutputKind, Profile};

    let profile = |output| Profile::Morgan {
        output,
        radius: 1,
        include_chirality: false,
        invariants: MorganInvariantKind::Connectivity,
        count_simulation: false,
    };
    let expected_hashed = vec![
        (1, 1),
        (80, 3),
        (222, 1),
        (294, 2),
        (482, 1),
        (807, 1),
        (1_057, 2),
        (1_420, 1),
        (1_544, 2),
    ];
    let expected_sparse_counts = vec![
        (864_662_311, 1),
        (1_506_563_592, 2),
        (1_535_166_686, 1),
        (2_245_273_601, 1),
        (2_245_384_272, 3),
        (2_246_728_737, 2),
        (3_098_934_668, 1),
        (3_542_456_614, 2),
        (4_022_716_898, 1),
    ];
    let expected_sparse_bits = [
        -2_049_693_695,
        -2_049_583_024,
        -2_048_238_559,
        -1_196_032_628,
        -752_510_682,
        -272_250_398,
        864_662_311,
        1_506_563_592,
        1_535_166_686,
    ];
    let expected_dense_bits = [1, 80, 222, 294, 482, 807, 1_057, 1_420, 1_544];
    let expected_folded_metadata = MorganAdditionalOutput {
        atom_counts: Some(vec![2; 7]),
        atom_to_bits: Some(vec![
            vec![1_057, 294],
            vec![80, 1_544],
            vec![1, 1_420],
            vec![80, 1_544],
            vec![1_057, 294],
            vec![80, 482],
            vec![807, 222],
        ]),
        bit_info_map: Some(vec![
            (1, vec![(2, 0)]),
            (80, vec![(1, 0), (3, 0), (5, 0)]),
            (222, vec![(6, 1)]),
            (294, vec![(0, 1), (4, 1)]),
            (482, vec![(5, 1)]),
            (807, vec![(6, 0)]),
            (1_057, vec![(0, 0), (4, 0)]),
            (1_420, vec![(2, 1)]),
            (1_544, vec![(1, 1), (3, 1)]),
        ]),
        bit_paths: Some(vec![]),
        atoms_per_bit: Some(vec![
            (1, vec![vec![2]]),
            (80, vec![vec![1], vec![3], vec![5]]),
            (222, vec![vec![6, 5]]),
            (294, vec![vec![0, 1], vec![4, 3]]),
            (482, vec![vec![5, 2, 6]]),
            (807, vec![vec![6]]),
            (1_057, vec![vec![0], vec![4]]),
            (1_420, vec![vec![2, 1, 3, 5]]),
            (1_544, vec![vec![1, 0, 2], vec![3, 2, 4]]),
        ]),
    };
    let expected_unfolded_metadata = MorganAdditionalOutput {
        atom_counts: Some(vec![2; 7]),
        atom_to_bits: Some(vec![
            vec![2_246_728_737, 3_542_456_614],
            vec![2_245_384_272, 1_506_563_592],
            vec![2_245_273_601, 3_098_934_668],
            vec![2_245_384_272, 1_506_563_592],
            vec![2_246_728_737, 3_542_456_614],
            vec![2_245_384_272, 4_022_716_898],
            vec![864_662_311, 1_535_166_686],
        ]),
        bit_info_map: Some(vec![
            (864_662_311, vec![(6, 0)]),
            (1_506_563_592, vec![(1, 1), (3, 1)]),
            (1_535_166_686, vec![(6, 1)]),
            (2_245_273_601, vec![(2, 0)]),
            (2_245_384_272, vec![(1, 0), (3, 0), (5, 0)]),
            (2_246_728_737, vec![(0, 0), (4, 0)]),
            (3_098_934_668, vec![(2, 1)]),
            (3_542_456_614, vec![(0, 1), (4, 1)]),
            (4_022_716_898, vec![(5, 1)]),
        ]),
        bit_paths: Some(vec![]),
        atoms_per_bit: Some(vec![
            (864_662_311, vec![vec![6]]),
            (1_506_563_592, vec![vec![1, 0, 2], vec![3, 2, 4]]),
            (1_535_166_686, vec![vec![6, 5]]),
            (2_245_273_601, vec![vec![2]]),
            (2_245_384_272, vec![vec![1], vec![3], vec![5]]),
            (2_246_728_737, vec![vec![0], vec![4]]),
            (3_098_934_668, vec![vec![2, 1, 3, 5]]),
            (3_542_456_614, vec![vec![0, 1], vec![4, 3]]),
            (4_022_716_898, vec![vec![5, 2, 6]]),
        ]),
    };

    for output_kind in [
        MorganOutputKind::DenseBits,
        MorganOutputKind::SparseBits,
        MorganOutputKind::HashedCounts,
        MorganOutputKind::SparseCounts,
    ] {
        let input = Input::Molecular {
            case: registry::SmilesCase {
                id: "fixed_rdkit_morgan_bit_info_executor".into(),
                smiles: "CCC(CC)CO".into(),
            },
            profile: profile(output_kind),
        };
        let record = molecular::run(&input).expect("the public executor returns a record");
        assert_eq!(record.input, input, "{output_kind:?}");
        let registry::Value::Molecular(outcome) = record.output else {
            panic!("Morgan result must remain a typed molecular outcome");
        };
        assert!(molecular::validate_output(&profile(output_kind), &outcome).is_ok());

        let additional_output = match &outcome {
            Outcome::MorganDenseBits(MorganDenseBitsOutput {
                additional_output, ..
            })
            | Outcome::MorganSparseBits(MorganSparseBitsOutput {
                additional_output, ..
            })
            | Outcome::MorganHashedCounts(MorganHashedCountsOutput {
                additional_output, ..
            })
            | Outcome::MorganSparseCounts(MorganSparseCountsOutput {
                additional_output, ..
            }) => additional_output,
            other => panic!("unexpected Morgan result family: {other:?}"),
        };
        let expected_metadata = match output_kind {
            MorganOutputKind::DenseBits | MorganOutputKind::HashedCounts => {
                &expected_folded_metadata
            }
            MorganOutputKind::SparseBits | MorganOutputKind::SparseCounts => {
                &expected_unfolded_metadata
            }
        };
        assert_eq!(additional_output, expected_metadata, "{output_kind:?}");

        match (output_kind, outcome) {
            (MorganOutputKind::DenseBits, Outcome::MorganDenseBits(actual)) => {
                assert_eq!(actual.length, 2048);
                assert_eq!(actual.on_bits, expected_dense_bits);
            }
            (MorganOutputKind::SparseBits, Outcome::MorganSparseBits(actual)) => {
                assert_eq!(actual.length, u32::MAX);
                assert_eq!(actual.on_bits, expected_sparse_bits);
            }
            (MorganOutputKind::HashedCounts, Outcome::MorganHashedCounts(actual)) => {
                assert_eq!(actual.length, 2048);
                assert_eq!(actual.entries, expected_hashed);
            }
            (MorganOutputKind::SparseCounts, Outcome::MorganSparseCounts(actual)) => {
                assert_eq!(actual.length, u64::MAX);
                assert_eq!(actual.entries, expected_sparse_counts);
            }
            (expected, actual) => panic!("{expected:?} produced {actual:?}"),
        }
    }
}

#[test]
fn parity_morgan_executor_retains_invalid_input_as_a_parse_diagnostic() {
    use molecular::{Outcome, Stage};
    use registry::molecule_plan::{MorganInvariantKind, MorganOutputKind, Profile};

    let input = Input::Molecular {
        case: registry::SmilesCase {
            id: "fixed_unclosed_ring_morgan_diagnostic".into(),
            smiles: "C1CC".into(),
        },
        profile: Profile::Morgan {
            output: MorganOutputKind::DenseBits,
            radius: 2,
            include_chirality: false,
            invariants: MorganInvariantKind::Connectivity,
            count_simulation: false,
        },
    };
    let record = molecular::run(&input).expect("chemistry errors remain diagnostic records");
    assert_eq!(record.input, input);
    let registry::Value::Molecular(Outcome::Error { stage, detail }) = record.output else {
        panic!("an invalid SMILES must remain a molecular error diagnostic");
    };
    assert_eq!(stage, Stage::Parse);
    assert!(!detail.is_empty());
}

fn morgan_preflight_tasks() -> Vec<&'static Task> {
    [
        "morgan_fingerprint",
        "morgan_sparse_fingerprint",
        "morgan_count_fingerprint",
        "morgan_sparse_count_fingerprint",
    ]
    .into_iter()
    .flat_map(|name| registry::select(Some(name)).expect("Morgan task is registered"))
    .collect()
}

fn morgan_preflight_corpus() -> Corpus {
    Corpus {
        bio_cases: vec![],
        fingerprints: vec![],
        molecules: vec![registry::SmilesCase {
            id: "empty-morgan-preflight".into(),
            smiles: String::new(),
        }],
    }
}

// These synthetic rows exercise the existing typed cache boundary for the
// empty molecule; they are not reference chemistry expectations.
fn morgan_preflight_record(input: &Input) -> Record {
    use molecular::{
        MorganAdditionalOutput, MorganDenseBitsOutput, MorganHashedCountsOutput,
        MorganSparseBitsOutput, MorganSparseCountsOutput, Outcome,
    };
    use registry::molecule_plan::{MorganOutputKind, Profile};

    let Profile::Morgan { output, .. } = (match input {
        Input::Molecular { profile, .. } => profile,
        _ => panic!("Morgan preflight input must be molecular"),
    }) else {
        panic!("Morgan preflight profile must use the Morgan variant")
    };
    let empty_additional_output = || MorganAdditionalOutput {
        atom_counts: Some(vec![]),
        atom_to_bits: Some(vec![]),
        bit_info_map: Some(vec![]),
        bit_paths: Some(vec![]),
        atoms_per_bit: Some(vec![]),
    };
    let output = match output {
        MorganOutputKind::DenseBits => Outcome::MorganDenseBits(MorganDenseBitsOutput {
            length: 2048,
            on_bits: vec![],
            additional_output: empty_additional_output(),
        }),
        MorganOutputKind::SparseBits => Outcome::MorganSparseBits(MorganSparseBitsOutput {
            length: u32::MAX,
            on_bits: vec![],
            additional_output: empty_additional_output(),
        }),
        MorganOutputKind::HashedCounts => Outcome::MorganHashedCounts(MorganHashedCountsOutput {
            length: 2048,
            entries: vec![],
            additional_output: empty_additional_output(),
        }),
        MorganOutputKind::SparseCounts => Outcome::MorganSparseCounts(MorganSparseCountsOutput {
            length: u64::MAX,
            entries: vec![],
            additional_output: empty_additional_output(),
        }),
    };
    Record {
        input: input.clone(),
        output: registry::Value::Molecular(output),
    }
}

fn morgan_preflight_fixture(data: &Path, task: &Task, cases: &Corpus) -> PathBuf {
    let inputs = registry::expand(cases, task);
    let records: Vec<_> = inputs.iter().map(morgan_preflight_record).collect();
    publish(data, task, &inputs, &records).unwrap();
    generation(data, task, &encode(&inputs).unwrap())
}

#[test]
fn parity_morgan_preflight_missing_last_task_blocks_every_rust_call() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = morgan_preflight_tasks();
    let cases = morgan_preflight_corpus();
    assert_eq!(tasks.len(), 4);
    for task in &tasks[..3] {
        assert_eq!(registry::expand(&cases, task).len(), 16);
        morgan_preflight_fixture(temp.path(), task, &cases);
    }

    let generated = Cell::new(0);
    let executed = Cell::new(0);
    let error = prepare_with(&tasks, &cases, temp.path(), |_, inputs| {
        generated.set(generated.get() + 1);
        assert_eq!(inputs.len(), 16);
        assert!(
            inputs
                .iter()
                .all(|input| input.task_name() == "morgan_sparse_count_fingerprint")
        );
        Err("last selected Morgan reference missing".into())
    })
    .err()
    .unwrap();
    assert!(error.contains("morgan_sparse_count_fingerprint"));
    assert!(error.contains("0 Rust operation calls"));
    assert!(error.contains("last selected Morgan reference missing"));
    assert_eq!(generated.get(), 1);
    let error = run(&tasks, &cases, temp.path(), |_| {
        executed.set(executed.get() + 1);
        Err("executor ran before every selected task was ready".into())
    })
    .unwrap_err();
    assert!(error.contains("morgan_sparse_count_fingerprint"));
    assert!(error.contains("0 Rust operation calls"));
    assert_eq!(executed.get(), 0);

    assert_eq!(
        preflight(&tasks[..3], &cases, temp.path()).unwrap().len(),
        48
    );
    let error = preflight(&tasks, &cases, temp.path()).err().unwrap();
    assert!(error.contains("morgan_sparse_count_fingerprint"));
    assert!(error.contains("0 Rust operation calls"));
}

#[test]
fn parity_morgan_preflight_rejects_missing_metadata_fields_with_current_identity() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = morgan_preflight_tasks();
    let task = tasks[0];
    let cases = morgan_preflight_corpus();
    let directory = morgan_preflight_fixture(temp.path(), task, &cases);
    let input = read(&directory.join("input.json")).unwrap();
    let inputs = registry::expand(&cases, task);
    let mut reference: serde_json::Value =
        serde_json::from_slice(&read(&directory.join("reference.json")).unwrap()).unwrap();
    let additional_output = &mut reference.as_array_mut().unwrap()[0]["record"]["output"]["Molecular"]
        ["MorganDenseBits"]["additional_output"];
    assert!(
        additional_output
            .as_object_mut()
            .unwrap()
            .remove("atom_counts")
            .is_some()
    );
    let reference = serde_json::to_vec_pretty(&reference).unwrap();
    fs::write(directory.join("reference.json"), &reference).unwrap();
    fs::write(
        directory.join("manifest.json"),
        encode(&identity(task, &input, &reference, inputs.len())).unwrap(),
    )
    .unwrap();

    let error = preflight(&[task], &cases, temp.path()).err().unwrap();
    assert!(error.contains("atom_counts"), "{error}");
}

#[test]
fn parity_morgan_preflight_rejects_per_atom_count_mismatch() {
    use molecular::{MorganDenseBitsOutput, Outcome};

    let temp = tempfile::tempdir().unwrap();
    let tasks = morgan_preflight_tasks();
    let task = tasks[0];
    let cases = morgan_preflight_corpus();
    let directory = morgan_preflight_fixture(temp.path(), task, &cases);
    let input = read(&directory.join("input.json")).unwrap();
    let inputs = registry::expand(&cases, task);
    let mut records: Vec<LabeledRecord> =
        serde_json::from_slice(&read(&directory.join("reference.json")).unwrap()).unwrap();
    let registry::Value::Molecular(Outcome::MorganDenseBits(MorganDenseBitsOutput {
        additional_output,
        ..
    })) = &mut records[0].record.output
    else {
        panic!("dense Morgan cache fixture must contain a dense output")
    };
    additional_output.atom_counts = Some(vec![0]);
    let reference = encode(&records).unwrap();
    fs::write(directory.join("reference.json"), &reference).unwrap();
    fs::write(
        directory.join("manifest.json"),
        encode(&identity(task, &input, &reference, inputs.len())).unwrap(),
    )
    .unwrap();

    let error = preflight(&[task], &cases, temp.path()).err().unwrap();
    assert!(error.contains("malformed or wrong-kind molecular reference result"));
}

#[test]
fn parity_morgan_preflight_rejects_reference_row_count_mismatch() {
    let temp = tempfile::tempdir().unwrap();
    let tasks = morgan_preflight_tasks();
    let task = tasks[0];
    let cases = morgan_preflight_corpus();
    let directory = morgan_preflight_fixture(temp.path(), task, &cases);
    let input = read(&directory.join("input.json")).unwrap();
    let inputs = registry::expand(&cases, task);
    let mut records: Vec<LabeledRecord> =
        serde_json::from_slice(&read(&directory.join("reference.json")).unwrap()).unwrap();
    assert_eq!(records.len(), 16);
    records.pop();
    let reference = encode(&records).unwrap();
    fs::write(directory.join("reference.json"), &reference).unwrap();
    fs::write(
        directory.join("manifest.json"),
        encode(&identity(task, &input, &reference, inputs.len())).unwrap(),
    )
    .unwrap();

    let error = preflight(&[task], &cases, temp.path()).err().unwrap();
    assert!(error.contains("reference row count mismatch"));
}

#[test]
fn parity_morgan_preflight_rejects_stale_profile_and_reference_identities() {
    use registry::molecule_plan::Profile;

    let task = morgan_preflight_tasks()[0];
    let cases = morgan_preflight_corpus();
    let temp = tempfile::tempdir().unwrap();
    let directory = morgan_preflight_fixture(temp.path(), task, &cases);
    let expected_inputs = registry::expand(&cases, task);
    let expected_input = read(&directory.join("input.json")).unwrap();
    assert_eq!(expected_input, encode(&expected_inputs).unwrap());
    let mut records: Vec<LabeledRecord> =
        serde_json::from_slice(&read(&directory.join("reference.json")).unwrap()).unwrap();
    let Input::Molecular {
        profile: Profile::Morgan { radius, .. },
        ..
    } = &mut records[0].record.input
    else {
        panic!("Morgan input must retain its serialized profile")
    };
    *radius = 3;
    records[0].label = reference_label(task, &expected_inputs[0]).unwrap();
    assert_eq!(
        records[0].label,
        reference_label(task, &expected_inputs[0]).unwrap()
    );
    let stale_reference = encode(&records).unwrap();
    fs::write(directory.join("reference.json"), &stale_reference).unwrap();
    fs::write(
        directory.join("manifest.json"),
        encode(&identity(
            task,
            &expected_input,
            &stale_reference,
            records.len(),
        ))
        .unwrap(),
    )
    .unwrap();
    let error = preflight(&[task], &cases, temp.path()).err().unwrap();
    assert!(
        error.contains("reference case/parameter mismatch"),
        "{error}"
    );

    for stale_field in ["schema", "registry", "oracle", "reference pin"] {
        let temp = tempfile::tempdir().unwrap();
        let directory = morgan_preflight_fixture(temp.path(), task, &cases);
        let mut manifest: Manifest =
            serde_json::from_slice(&read(&directory.join("manifest.json")).unwrap()).unwrap();
        match stale_field {
            "schema" => manifest.schema = 2,
            "registry" => manifest.registry_sha256 = "stale-registry".into(),
            "oracle" => manifest.oracle_sha256 = "stale-oracle".into(),
            "reference pin" => manifest.reference_pin_sha256 = "stale-pin".into(),
            _ => unreachable!(),
        }
        fs::write(directory.join("manifest.json"), encode(&manifest).unwrap()).unwrap();
        let error = preflight(&[task], &cases, temp.path()).err().unwrap();
        assert!(
            error.contains("stale or corrupted manifest/input/reference"),
            "{stale_field}: {error}"
        );
    }
}
