//! Author delivery proposal, not independent acceptance or a replacement checker.
//! Transports all original 131 profiles using the ROOT-frozen random200
//! selection and 112 workers; original5000 reference/matcher bytes are retained.
//! through the public batch interface. No reference implementation is executed.
use cosmolkit as ck;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    error::Error,
    fs,
    io::{BufWriter, Write},
    path::PathBuf,
};

fn bytes(path: &str, expected: &str) -> Vec<u8> {
    let value = fs::read(path).unwrap_or_else(|e| panic!("required original input {path}: {e}"));
    assert_eq!(
        Sha256::digest(&value)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>(),
        expected,
        "original input checksum: {path}"
    );
    value
}
fn u(v: &Value, key: &str) -> u32 {
    u32::try_from(v[key].as_u64().unwrap()).unwrap()
}
fn b(v: &Value, key: &str) -> bool {
    v[key].as_bool().unwrap()
}
fn vector(v: &Value) -> Option<Vec<u32>> {
    if v.is_null() {
        None
    } else {
        Some(
            v.as_array()
                .unwrap()
                .iter()
                .map(|x| u32::try_from(x.as_u64().unwrap()).unwrap())
                .collect(),
        )
    }
}
fn selected(v: &Value, key: &str, n: usize) -> Option<Vec<u32>> {
    (v[key].as_str() == Some("first") && n != 0).then(|| vec![0])
}
fn ap(v: &Value, molecule: &ck::Molecule) -> ck::AtomPairFingerprintParams {
    ck::AtomPairFingerprintParams {
        generator: ck::AtomPairParams {
            min_distance: u(v, "minDistance"),
            max_distance: u(v, "maxDistance"),
            include_chirality: b(v, "includeChirality"),
            use_2d: b(v, "use2D"),
            count_simulation: b(v, "countSimulation"),
            fp_size: u(v, "fpSize"),
            count_bounds: vector(&v["countBounds"]).unwrap(),
            bits_per_feature: u(v, "numBitsPerFeature"),
        },
        from_atoms: selected(v, "fromAtoms", molecule.num_atoms()),
        ignore_atoms: selected(v, "ignoreAtoms", molecule.num_atoms()),
        custom_atom_invariants: (v["customAtomInvariants"].as_str() == Some("index_plus_11")).then(
            || {
                (0..molecule.num_atoms())
                    .map(|i| u32::try_from(i + 11).unwrap())
                    .collect()
            },
        ),
        ..Default::default()
    }
}
fn morgan(v: &Value, molecule: &ck::Molecule) -> ck::MorganFingerprintParams {
    let mut generator = ck::MorganParams::default();
    generator.radius = u(v, "radius");
    generator.fp_size = u(v, "nBits");
    generator.include_chirality = b(v, "includeChirality");
    generator.use_bond_types =
        if v["bondInvariantsGenerator"].as_str() == Some("morgan_no_bond_types") {
            false
        } else {
            b(v, "useBondTypes")
        };
    generator.include_ring_membership =
        if v["atomInvariantsGenerator"].as_str() == Some("morgan_ring_false") {
            false
        } else {
            v["includeRingMembership"].as_bool().unwrap_or(true)
        };
    generator.count_simulation = v["countSimulation"].as_bool().unwrap_or(false);
    generator.count_bounds =
        vector(&v["countBounds"]).unwrap_or_else(|| ck::MorganParams::default().count_bounds);
    generator.only_nonzero_invariants = v["onlyNonzeroInvariants"].as_bool().unwrap_or(false);
    generator.include_redundant_environments =
        v["includeRedundantEnvironments"].as_bool().unwrap_or(false);
    generator.bits_per_feature = v["numBitsPerFeature"]
        .as_u64()
        .map(|x| u32::try_from(x).unwrap())
        .unwrap_or(1);
    let custom_atom_invariants = match v["customAtomInvariants"].as_str() {
        Some("index_plus_one") => Some(
            (0..molecule.num_atoms())
                .map(|i| u32::try_from(i + 1).unwrap())
                .collect(),
        ),
        Some("even_atoms_zero_odd_index_plus_one") => Some(
            (0..molecule.num_atoms())
                .map(|i| {
                    if i % 2 == 0 {
                        0
                    } else {
                        u32::try_from(i + 1).unwrap()
                    }
                })
                .collect(),
        ),
        None => None,
        other => panic!("unmodeled original invariant option {other:?}"),
    };
    ck::MorganFingerprintParams {
        generator,
        from_atoms: selected(v, "fromAtoms", molecule.num_atoms()),
        ignore_atoms: selected(v, "ignoreAtoms", molecule.num_atoms()),
        custom_atom_invariants,
        custom_bond_invariants: (v["customBondInvariants"].as_str() == Some("index_plus_seven"))
            .then(|| {
                (0..molecule.num_bonds())
                    .map(|i| u32::try_from(i + 7).unwrap())
                    .collect()
            }),
        invariants: if v["atomInvariantsGenerator"].as_str() == Some("morgan_feature_default") {
            ck::MorganInvariants::Features
        } else {
            ck::MorganInvariants::Connectivity
        },
        ..Default::default()
    }
}
fn layered(v: &Value, reference: &Value) -> ck::LayeredFingerprintParams {
    let a = &reference["resolved_arguments"];
    ck::LayeredFingerprintParams {
        layers: ck::LayeredFingerprintLayers::from_bits_retain(u(a, "layerFlags")),
        min_path: u(v, "minPath"),
        max_path: u(v, "maxPath"),
        fp_size: u(v, "fpSize"),
        branched_paths: b(v, "branchedPaths"),
        from_atoms: vector(&a["fromAtoms"]),
        atom_counts: vector(&a["atomCounts"]),
        set_only_bits: vector(&a["setOnlyBits"])
            .map(|bits| ck::Fingerprint::from_on_bits(u(v, "fpSize"), bits).unwrap()),
    }
}
fn bits(value: &ck::Fingerprint, family: &str) -> Value {
    let mut v = json!({"ok":true,"error":null,"on_bits":value.on_bits()});
    v[match family {
        "atom_pair" => "length",
        "morgan" => "num_on_bits",
        "layered" => "num_bits",
        "pattern" => "n_bits",
        _ => unreachable!(),
    }] = json!(if family == "morgan" {
        value.on_bits().len() as u64
    } else {
        u64::from(value.n_bits())
    });
    v
}
fn output(value: &ck::BatchFingerprintOutput, family: &str, collect: bool) -> Value {
    let mut v = bits(value.fingerprint(), family);
    if collect {
        let a = value.additional_output().unwrap();
        assert!(a.bit_paths().is_none());
        v["additional_output"] = json!({"atom_counts":a.atom_counts(),"atom_to_bits":a.atom_to_bits(),"bit_info_map":a.bit_info_map(),"atoms_per_bit":a.atoms_per_bit()});
    } else {
        assert!(value.additional_output().is_err());
    }
    v
}
fn error(e: &ck::BatchValidationError) -> Value {
    let mut records = Vec::new();
    for r in &e.record_errors {
        let mut chain = Vec::new();
        let mut current: Option<&(dyn Error + 'static)> = Some(r);
        while let Some(c) = current {
            let mut v = json!({"message":c.to_string(),"debug":format!("{c:?}")});
            if let Some(ck::FingerprintError::SparseIndexOutOfRange { index, size }) =
                c.downcast_ref::<ck::FingerprintError>()
            {
                v["kind"] = json!("SparseIndexOutOfRange");
                v["index"] = json!(index);
                v["size"] = json!(size);
            }
            chain.push(v);
            current = c.source();
        }
        records.push(json!({"index":r.index,"operation":r.operation,"message":r.message,"cause_chain":chain}));
    }
    json!({"ok":false,"error":e.to_string(),"error_count":e.errors,"record_errors":records})
}
fn invoke(
    batch: &ck::MoleculeBatch,
    molecule: &ck::Molecule,
    family: &str,
    form: &str,
    branch: &Value,
    reference: &Value,
    jobs: usize,
) -> Result<Vec<Value>, ck::BatchValidationError> {
    let query = ck::BatchQueryParams {
        n_jobs: Some(jobs),
        progress_bar: Some(false),
        ..Default::default()
    };
    let collect = branch["additionalOutput"].as_bool().unwrap_or(false);
    macro_rules! project {
        ($expr:expr,$map:expr) => {
            $expr.map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        ($map)(value.expect("original parsed valid input must produce a value"))
                    })
                    .collect()
            })
        };
    }
    match (family, form) {
        ("atom_pair", "explicit_bit") => project!(
            batch.fingerprint_atom_pair_list_with_params(&ap(branch, molecule), &query),
            |x: ck::Fingerprint| bits(&x, family)
        ),
        ("atom_pair", "sparse_bit") => project!(
            batch.fingerprint_atom_pair_sparse_bits_list_with_params(&ap(branch, molecule), &query),
            |x: ck::SparseBitFingerprint| json!({"ok":true,"error":null,"length":x.n_bits(),"on_bits":x.on_bits()})
        ),
        ("atom_pair", "sparse_count") => project!(
            batch
                .fingerprint_atom_pair_sparse_count_list_with_params(&ap(branch, molecule), &query),
            |x: ck::SparseCountFingerprint| json!({"ok":true,"error":null,"length":x.length(),"nonzero_elements":x.nonzero_elements()})
        ),
        ("atom_pair", "count") => project!(
            batch.fingerprint_atom_pair_count_list_with_params(&ap(branch, molecule), &query),
            |x: ck::SparseCountFingerprint32| json!({"ok":true,"error":null,"length":x.length(),"nonzero_elements":x.nonzero_elements()})
        ),
        ("atom_pair", "with_output") => project!(
            batch.fingerprint_atom_pair_with_output_list_with_params(
                &ap(branch, molecule),
                collect,
                &query
            ),
            |x: ck::BatchFingerprintOutput| output(&x, family, collect)
        ),
        ("morgan", "explicit_bit") => project!(
            batch.fingerprint_morgan_list_with_params(&morgan(branch, molecule), &query),
            |x: ck::Fingerprint| {
                assert_eq!(x.n_bits(), u(branch, "nBits"));
                bits(&x, family)
            }
        ),
        ("morgan", "with_output") => project!(
            batch.fingerprint_morgan_with_output_list_with_params(
                &morgan(branch, molecule),
                collect,
                &query
            ),
            |x: ck::BatchFingerprintOutput| {
                assert_eq!(x.fingerprint().n_bits(), u(branch, "nBits"));
                output(&x, family, collect)
            }
        ),
        ("layered", "explicit_bit") => project!(
            batch.fingerprint_layered_list_with_params(&layered(branch, reference), &query),
            |x: ck::Fingerprint| bits(&x, family)
        ),
        ("layered", "with_output") => project!(
            batch.fingerprint_layered_with_output_list_with_params(
                &layered(branch, reference),
                &query
            ),
            |x: ck::LayeredFingerprintResult| {
                let mut v = bits(x.fingerprint(), family);
                v["atom_counts"] = json!(x.atom_counts());
                v
            }
        ),
        ("pattern", "explicit_bit") => project!(
            batch.pattern_fingerprint_list_with_params(
                &ck::PatternFingerprintParams {
                    n_bits: usize::try_from(u(branch, "fpSize")).unwrap(),
                    tautomeric: b(branch, "tautomericFingerprint")
                },
                &query
            ),
            |x: ck::Fingerprint| bits(&x, family)
        ),
        _ => panic!("unmodeled frozen original profile {family}/{form}"),
    }
}
fn expected(family: &str, form: &str, source: &Value) -> Value {
    if family == "atom_pair" || family == "morgan" {
        let mut v = source[if form == "with_output" {
            "explicit_bit"
        } else {
            form
        }]
        .clone();
        if form != "with_output" {
            v.as_object_mut().unwrap().remove("additional_output");
        }
        v
    } else {
        let mut v = json!({});
        for k in [
            "ok",
            "error",
            "on_bits",
            if family == "layered" {
                "num_bits"
            } else {
                "n_bits"
            },
        ] {
            v[k] = source[k].clone();
        }
        if family == "layered" && form == "with_output" {
            v["atom_counts"] = source["atom_counts"].clone();
        }
        v
    }
}
fn matches(actual: &Value, target: &Value, family: &str, branch: &Value) -> bool {
    if target["ok"] == true {
        return actual == target;
    }
    if actual["ok"] == true {
        return false;
    }
    if family == "atom_pair" {
        if let Some(raw) = target["error"]
            .as_str()
            .and_then(|s| s.strip_prefix("IndexError: "))
        {
            let index = raw.parse::<u64>().unwrap();
            let size = 1_u64
                << if b(branch, "includeChirality") {
                    27
                } else {
                    23
                };
            return actual["record_errors"]
                .as_array()
                .unwrap()
                .iter()
                .flat_map(|r| r["cause_chain"].as_array().unwrap())
                .any(|c| {
                    c["kind"] == "SparseIndexOutOfRange" && c["index"] == index && c["size"] == size
                });
        }
    }
    // Retain the original literal branch as a delivery proposal. Canonical
    // exception spelling alone is not independently accepted chemistry parity.
    actual["error"] == target["error"]
}
fn molecule(batch: &ck::MoleculeBatch, i: usize) -> &ck::Molecule {
    match batch.get(i).unwrap() {
        ck::BatchRecord::Molecule(m) => m,
        _ => panic!("original row {i} failed parsing"),
    }
}

// This adapter retains two distinct calls, including the raw invalid-mask
// error. It never claims the public API received that mask. Full pinned update
// body and p1 input review are hash-bound before any CK call; no source parsing
// or oracle execution produces expected values.
fn preflight_pattern_input_boundary(rows: &[Value], branches: &Value) {
    // RDKit❗✔️:   PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
    // RDKit❗✔️:                "bad setOnlyBits size");
    // Pin 351f8f378f8ad6bbd517980c38896e66bf907af8, full
    // PatternFingerprints.cpp SHA7618dcfa67271c9750780e5a742bdef6bce76b2d1c4a498f7205396ce0e23571.
    // This is input-identity validation (linear in rows), not a production port
    // or an independently accepted assertion about public invalid-mask parity.
    let proof = PathBuf::from(std::env::var("COSMOLKIT_BATCH_TRANSPORT_PROOF").unwrap());
    for (name, checksum) in [
        (
            "PatternFingerprints.cpp",
            "7618dcfa67271c9750780e5a742bdef6bce76b2d1c4a498f7205396ce0e23571",
        ),
        (
            "pattern-golden-generator.py",
            "0766fab72e344e6f53a9f7a14a1fb9d267067a75febbd483f5c543f2eb116dc7",
        ),
        (
            "p1-readonly-condition-proposal.json",
            "d8bf4c6cb2e23dd187e459329cf7764b78031afb7859a116a78caccfb5399be5",
        ),
    ] {
        bytes(proof.join(name).to_str().unwrap(), checksum);
    }
    let raw_error = "RuntimeError: Pre-condition Violation\n\tbad setOnlyBits size\n\tViolation occurred on line 346 in file Code/GraphMol/Fingerprints/PatternFingerprints.cpp\n\tFailed Expression: !setOnlyBits || setOnlyBits->getNumBits() == fpSize\n\tRDKIT: 2026.03.1\n\tBOOST: \n";
    for row in rows {
        for branch in branches.as_array().unwrap() {
            assert_eq!(
                row["branches"][branch["name"].as_str().unwrap()]["parameters"],
                *branch
            );
        }
        let raw = &row["branches"]["set_only_bits_wrong_width"];
        assert_eq!(
            raw["parameters"],
            json!({"atomCounts":"none","fpSize":257,"name":"set_only_bits_wrong_width","setOnlyBits":"wrong_width","tautomericFingerprint":false})
        );
        assert_eq!(raw["ok"], false);
        assert_eq!(raw["error"], raw_error);
        assert!(raw["n_bits"].is_null() && raw["on_bits"].is_null());
        assert!(raw["atom_counts_before"].is_null() && raw["atom_counts_after"].is_null());
        assert_eq!(raw["set_only_bits"], json!([]));
        let valid = &row["branches"]["set_only_bits_zero"];
        assert_eq!(
            valid["parameters"],
            json!({"atomCounts":"none","fpSize":257,"name":"set_only_bits_zero","setOnlyBits":"zero","tautomericFingerprint":false})
        );
        assert_eq!(valid["ok"], true);
        assert!(valid["error"].is_null());
        assert_eq!(valid["n_bits"], 257);
        assert_eq!(valid["set_only_bits"], json!([]));
        assert!(
            valid["on_bits"]
                .as_array()
                .unwrap()
                .iter()
                .all(|x| x.as_u64().unwrap() < 257)
        );
    }
}

fn pattern_raw_arguments(source: &Value) -> Value {
    let p = &source["parameters"];
    let width = match p["setOnlyBits"].as_str().unwrap() {
        "none" => Value::Null,
        "wrong_width" => json!(u(p, "fpSize") + 1),
        _ => p["fpSize"].clone(),
    };
    json!({"fpSize":p["fpSize"],"tautomericFingerprint":p["tautomericFingerprint"],"atomCounts":source["atom_counts_before"],"setOnlyBits":{"mode":p["setOnlyBits"],"width":width,"on_bits":source["set_only_bits"]},"original_parameters":p})
}

fn current_user_conditions() -> (Vec<usize>, usize) {
    let condition_path = std::env::var("COSMOLKIT_BATCH_CURRENT_USER_CONDITIONS")
        .expect("ROOT-frozen current user BATCH conditions required");
    let condition: Value = serde_json::from_slice(&bytes(
        &condition_path,
        "5e6dac1276f736f9cfa1a3955835de47f7603d79ad891870f9c2f10e4b06f6e0",
    ))
    .unwrap();
    let selection_path = std::env::var("COSMOLKIT_BATCH_RANDOM200_SELECTION")
        .expect("ROOT-frozen same random200 selection required");
    let selection: Value = serde_json::from_slice(&bytes(
        &selection_path,
        "bf26ef80887f11dca72ea0e3d267429a874ba5fa5051a1e6839d91425eb46c5e",
    ))
    .unwrap();
    assert_eq!(selection["population_size"], 5000);
    assert_eq!(selection["sample_size"], 200);
    assert_eq!(selection["without_replacement"], true);
    assert_eq!(selection["seed"], condition["seed"]);
    assert_eq!(condition["sample_size"], 200);
    assert_eq!(condition["expected_batch_corpus_comparisons"], 26200);
    assert_eq!(
        condition["random_sample_sha256"],
        "bf26ef80887f11dca72ea0e3d267429a874ba5fa5051a1e6839d91425eb46c5e"
    );
    let indices: Vec<usize> = selection["indices_zero_based"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| usize::try_from(v.as_u64().unwrap()).unwrap())
        .collect();
    assert_eq!(indices.len(), 200);
    assert!(indices.windows(2).all(|v| v[0] < v[1]));
    assert!(indices.iter().all(|i| *i < 5000));
    let workers = usize::try_from(condition["batch_workers"].as_u64().unwrap()).unwrap();
    assert_eq!(workers, 112);
    assert_eq!(
        workers,
        usize::try_from(condition["available_cpus"].as_u64().unwrap()).unwrap() / 2
    );
    (indices, workers)
}

#[test]
fn all_original131_profiles_fixed_random200_workers112() {
    let (original_indices, worker_count) = current_user_conditions();
    let sample_len = original_indices.len();
    let selected_indices: std::collections::BTreeSet<usize> =
        original_indices.iter().copied().collect();
    let freeze_path = std::env::var("COSMOLKIT_BATCH_REFERENCE_FREEZE").expect(
        "required pinned original reference-freeze.json; preparation never runs an oracle here",
    );
    let freeze_bytes = bytes(
        &freeze_path,
        "0fb54512f712f13f88e6649638d27968596b07391475caaa109da27513a060be",
    );
    let freeze: Value = serde_json::from_slice(&freeze_bytes).unwrap();
    assert_eq!(freeze["acceptance"], false);
    assert_eq!(freeze["jobs"], json!([1, 4]));
    assert_eq!(freeze["filters"], json!([]));
    assert_eq!(freeze["pin"], "351f8f378f8ad6bbd517980c38896e66bf907af8");
    let corpus_path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../testdata/smiles/corpus/smiles_5000.smi");
    let corpus_bytes = bytes(
        corpus_path.to_str().unwrap(),
        freeze["corpus"]["sha256"].as_str().unwrap(),
    );
    let all_smiles: Vec<String> = std::str::from_utf8(&corpus_bytes)
        .unwrap()
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_owned)
        .collect();
    assert_eq!(all_smiles.len(), 5000);
    let smiles: Vec<String> = original_indices
        .iter()
        .map(|i| all_smiles[*i].clone())
        .collect();
    // Preflight and own ALL four original reference snapshots before ANY CK call.
    // These branches are the frozen original proposal's input, not a second
    // production registry or any modification to the independent checker.
    let mut references = BTreeMap::new();
    let mut profile_count = 0;
    for family in ["atom_pair", "morgan", "layered", "pattern"] {
        let spec = &freeze["families"][family];
        let snapshots = PathBuf::from(
            std::env::var("COSMOLKIT_BATCH_REFERENCE_SNAPSHOTS")
                .expect("prospectively frozen byte-identical original snapshots required"),
        );
        let data = bytes(
            snapshots.join(format!("{family}.jsonl")).to_str().unwrap(),
            spec["golden_sha256"].as_str().unwrap(),
        );
        assert_eq!(data.len() as u64, spec["golden_bytes"].as_u64().unwrap());
        bytes(
            snapshots
                .join(format!("{family}-profile.json"))
                .to_str()
                .unwrap(),
            spec["profile_sha256"].as_str().unwrap(),
        );
        let branches: Vec<&str> = spec["branches"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v["name"].as_str().unwrap())
            .collect();
        assert_eq!(
            branches.len(),
            match family {
                "atom_pair" => 10,
                "morgan" => 17,
                "layered" => 18,
                "pattern" => 11,
                _ => unreachable!(),
            }
        );
        profile_count += branches.len()
            * match family {
                "atom_pair" => 5,
                "pattern" => 1,
                _ => 2,
            };
        let mut rows = Vec::new();
        let mut offset = 0;
        let mut population_count = 0;
        for (i, line) in data.split_inclusive(|x| *x == b'\n').enumerate() {
            assert_eq!(spec["offsets"][i], offset);
            offset += line.len() as u64;
            population_count += 1;
            if !selected_indices.contains(&i) {
                continue;
            }
            let row: Value = serde_json::from_slice(line).unwrap();
            assert_eq!(row["smiles"], all_smiles[i]);
            let keys: std::collections::BTreeSet<&str> = row["branches"]
                .as_object()
                .unwrap()
                .keys()
                .map(String::as_str)
                .collect();
            assert_eq!(keys, branches.iter().copied().collect());
            rows.push(row);
        }
        assert_eq!(population_count, 5000);
        assert_eq!(rows.len(), sample_len);
        references.insert(family, rows);
    }
    assert_eq!(profile_count, 131);
    preflight_pattern_input_boundary(
        &references["pattern"],
        &freeze["families"]["pattern"]["branches"],
    );
    let output_dir = PathBuf::from(
        std::env::var("COSMOLKIT_BATCH_TRANSPORT_OUTPUT")
            .expect("explicit ignored evidence output directory required"),
    );
    fs::create_dir_all(&output_dir).unwrap();
    let create = |name: &str| {
        BufWriter::new(
            fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(output_dir.join(name))
                .unwrap(),
        )
    };
    let mut observations = create("full-ordered-observations.jsonl");
    let mut failures = create("all-failures.jsonl");
    let mut pattern_witness = create("Pattern-reference-transport-witness.jsonl");
    let mut raw_boundary_observations = 0;
    let mut single_peer_recovery_calls = 0;
    let batch = ck::MoleculeBatch::from_smiles_list_with_params(
        &smiles,
        &ck::SmilesParseParams::default(),
        &ck::BatchParams {
            errors: Some(ck::BatchErrorMode::KeepErrors),
            n_jobs: Some(1),
            progress_bar: Some(false),
        },
    )
    .unwrap();
    assert_eq!(batch.valid_mask(), vec![true; sample_len]);
    assert!(batch.errors().is_empty());
    assert_eq!(batch.len(), sample_len);
    let mut failed_profiles = 0;
    let mut comparisons = 0;
    let mut mismatches = 0;
    let mut executed_profiles = 0;
    for family in ["atom_pair", "morgan", "layered", "pattern"] {
        let forms: &[&str] = match family {
            "atom_pair" => &[
                "explicit_bit",
                "sparse_count",
                "count",
                "sparse_bit",
                "with_output",
            ],
            "pattern" => &["explicit_bit"],
            _ => &["explicit_bit", "with_output"],
        };
        for branch in freeze["families"][family]["branches"].as_array().unwrap() {
            let name = branch["name"].as_str().unwrap();
            let refs: Vec<&Value> = references[family]
                .iter()
                .map(|r| &r["branches"][name])
                .collect();
            // Keep original raw rows and input grouping. Only the explicitly
            // different public wrong-width call uses the same-row valid zero
            // mask reference; matcher and raw Error bytes remain unchanged.
            let public_reference_branch =
                if family == "pattern" && name == "set_only_bits_wrong_width" {
                    "set_only_bits_zero"
                } else {
                    name
                };
            let public_refs: Vec<&Value> = references[family]
                .iter()
                .map(|r| &r["branches"][public_reference_branch])
                .collect();
            let mut keys = BTreeMap::<String, usize>::new();
            let mut groups = Vec::<Vec<usize>>::new();
            for (i, reference) in refs.iter().enumerate() {
                let m = molecule(&batch, i);
                let key = if family == "layered" {
                    serde_json::to_string(&reference["resolved_arguments"]).unwrap()
                } else {
                    serde_json::to_string(&json!([
                        if !branch["customAtomInvariants"].is_null()
                            || !branch["fromAtoms"].is_null()
                            || !branch["ignoreAtoms"].is_null()
                        {
                            Some(m.num_atoms())
                        } else {
                            None
                        },
                        if !branch["customBondInvariants"].is_null() {
                            Some(m.num_bonds())
                        } else {
                            None
                        }
                    ]))
                    .unwrap()
                };
                let group = *keys.entry(key).or_insert_with(|| {
                    groups.push(Vec::new());
                    groups.len() - 1
                });
                groups[group].push(i);
            }
            let mut coverage: Vec<usize> = groups.iter().flatten().copied().collect();
            coverage.sort_unstable();
            assert_eq!(coverage, (0..sample_len).collect::<Vec<_>>());
            for form in forms {
                let start_mismatches = mismatches;
                let mut serial = None;
                for jobs in [worker_count] {
                    let mut values = vec![Value::Null; sample_len];
                    let mut aggregates = Vec::new();
                    for ids in &groups {
                        let subset = ck::MoleculeBatch::from_records(
                            ids.iter().map(|i| batch.get(*i).unwrap().clone()).collect(),
                            ck::BatchErrorMode::KeepErrors,
                        )
                        .unwrap();
                        match invoke(
                            &subset,
                            molecule(&batch, ids[0]),
                            family,
                            form,
                            branch,
                            refs[ids[0]],
                            jobs,
                        ) {
                            Ok(outputs) => {
                                assert_eq!(outputs.len(), ids.len());
                                for (&i, v) in ids.iter().zip(outputs) {
                                    values[i] = v;
                                }
                            }
                            Err(e) => {
                                let indices: Vec<usize> =
                                    e.record_errors.iter().map(|r| ids[r.index]).collect();
                                let expected_indices: Vec<usize> = ids
                                    .iter()
                                    .copied()
                                    .filter(|i| {
                                        expected(family, form, public_refs[*i])["ok"] != true
                                    })
                                    .collect();
                                if indices != expected_indices {
                                    mismatches += 1;
                                    writeln!(failures,"{}",json!({"family":family,"branch":name,"form":form,"jobs":jobs,"kind":"aggregate_indices","actual":indices,"expected":expected_indices})).unwrap();
                                }
                                aggregates.push(json!({"global_indices":ids,"error":error(&e)}));
                                // Execute EVERY original peer after the aggregate failed,
                                // preserving both the bulk error and successful peer values.
                                for &i in ids {
                                    single_peer_recovery_calls += 1;
                                    let one = ck::MoleculeBatch::from_records(
                                        vec![batch.get(i).unwrap().clone()],
                                        ck::BatchErrorMode::KeepErrors,
                                    )
                                    .unwrap();
                                    values[i] = match invoke(
                                        &one,
                                        molecule(&batch, i),
                                        family,
                                        form,
                                        branch,
                                        refs[i],
                                        1,
                                    ) {
                                        Ok(mut v) => {
                                            assert_eq!(v.len(), 1);
                                            v.remove(0)
                                        }
                                        Err(e) => error(&e),
                                    };
                                }
                            }
                        }
                    }
                    assert!(values.iter().all(|v| !v.is_null()));
                    for i in 0..sample_len {
                        let target = expected(family, form, public_refs[i]);
                        comparisons += 1;
                        if !matches(&values[i], &target, family, branch) {
                            mismatches += 1;
                            writeln!(failures,"{}",json!({"family":family,"branch":name,"form":form,"jobs":jobs,"row":i,"kind":"original","expected":target,"actual":values[i]})).unwrap();
                        }
                    }
                    if family == "pattern" {
                        let witnesses: Vec<Value> = (0..sample_len).map(|i| json!({
                            "row":original_indices[i],"selected_position":i,"smiles":smiles[i],
                            "raw_source_arguments":pattern_raw_arguments(refs[i]),
                            "raw_source_result":refs[i],
                            "actual_public_arguments":{"n_bits":branch["fpSize"],"tautomeric":branch["tautomericFingerprint"],"atomCounts":null,"setOnlyBits":null,"n_jobs":jobs},
                            "public_reference_branch":public_reference_branch,
                            "public_equivalent_reference_arguments":pattern_raw_arguments(public_refs[i]),
                            "public_equivalent_reference_result":public_refs[i],
                            "public_expected_projection":expected(family,form,public_refs[i]),
                            "actual_public_result":values[i],
                            "raw_invalid_mask_was_passed_to_public":false
                        })).collect();
                        writeln!(
                            pattern_witness,
                            "{}",
                            json!({"branch":name,"n_jobs":jobs,"case_count":sample_len,"original_indices":original_indices,"rows":witnesses})
                        )
                        .unwrap();
                        pattern_witness.flush().unwrap();
                        if name == "set_only_bits_wrong_width" {
                            raw_boundary_observations += sample_len;
                        }
                    }
                    writeln!(observations,"{}",json!({"family":family,"branch":name,"form":form,"n_jobs":jobs,"case_count":sample_len,"original_indices":original_indices,"partition_global_indices":groups,"aggregate_errors":aggregates,"values":values})).unwrap();
                    observations.flush().unwrap();
                    if let Some(previous) = &serial {
                        if previous != &values {
                            mismatches += 1;
                            writeln!(failures,"{}",json!({"family":family,"branch":name,"form":form,"kind":"serial_parallel_values_differ"})).unwrap();
                        }
                    } else {
                        serial = Some(values);
                    }
                }
                executed_profiles += 1;
                if start_mismatches != mismatches {
                    failed_profiles += 1;
                }
                failures.flush().unwrap();
                println!(
                    "PROFILE {executed_profiles}/131 {family}/{name}/{form}: comparisons=200 workers=112 mismatches={}",
                    mismatches - start_mismatches
                );
            }
        }
    }
    assert_eq!(executed_profiles, 131);
    assert_eq!(comparisons, 26_200);
    assert_eq!(raw_boundary_observations, 200);
    fs::write(output_dir.join("summary.json"),serde_json::to_vec_pretty(&json!({"profiles":executed_profiles,"passed_profiles":executed_profiles-failed_profiles,"failed_profiles":failed_profiles,"comparisons":comparisons,"comparison_scope":"public-equivalent-input; original raw wrong-width Error retained separately","raw_wrong_width_boundary_observations":raw_boundary_observations,"raw_invalid_mask_public_parity_claimed":false,"mismatches":mismatches,"original_population_cases":5000,"sample_cases":sample_len,"original_indices":original_indices,"selection_sha256":"bf26ef80887f11dca72ea0e3d267429a874ba5fa5051a1e6839d91425eb46c5e","conditions_sha256":"5e6dac1276f736f9cfa1a3955835de47f7603d79ad891870f9c2f10e4b06f6e0","jobs":[worker_count],"additional_single_peer_recovery_calls":single_peer_recovery_calls,"additional_single_peer_recovery_workers":1,"filters":[],"proposal_only":true,"acceptance":false,"Python_tests_native_probes":"USER_SKIP"})).unwrap()).unwrap();
    assert_eq!(
        mismatches,
        0,
        "full original proposal mismatches remain; all rows/profile results retained in {}",
        output_dir.display()
    );
}
