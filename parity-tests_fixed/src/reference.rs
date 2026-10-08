use crate::{
    Corpus, Input, Record, Result, descriptors, digest, directory, encode, read, registry, root,
    workflow::Spec,
};
use serde_json::{Value, json};
use std::{
    io::Write,
    process::{Command, Stdio},
};

pub(crate) fn uses_gemmi(spec: &Spec) -> bool {
    matches!(*spec, Spec::Corpus(t) if t.operation == registry::Operation::BioPdbOutput)
        || matches!(*spec, Spec::Special(s) if matches!(s.schema, registry::SpecialRegressionSchema::BioMmcifSwitches))
}

pub(crate) fn source_digest(spec: &Spec) -> Result<String> {
    // JSON decoding is part of preparation: float_roundtrip must preserve the
    // native f64 inputs, including rank-deficient alignment coordinates.
    let mut paths = vec![
        directory().join("tools/reference.py"),
        directory().join("Cargo.toml"),
    ];
    if matches!(*spec, Spec::Corpus(t) if matches!(t.operation, registry::Operation::Fingerprint(_)))
    {
        paths.push(directory().join("tools/fingerprints.py"));
    }
    if uses_gemmi(spec) {
        paths.push(directory().join("testdata/reference/gemmi.json"));
    } else {
        paths.push(directory().join("src/descriptors.rs"));
        paths.extend(
            [
                "tools/oracles/rdkit/fingerprint_values_pilot.py",
                "tools/testdata/rdkit/_generate_molecular_descriptors_golden.py",
                "tools/testdata/rdkit/_generate_smiles_writer_golden.py",
                "tools/testdata/rdkit/_tautomer_oracle.py",
                "tools/testdata/rdkit/tautomer_profile.json",
                "tools/testdata/rdkit/_generate_tetrahedral_stereo_geometry.py",
                "tools/testdata/rdkit/_generate_molalign_golden.py",
            ]
            .map(|p| root().join(p)),
        );
    }
    paths.sort();
    let identities: Vec<_> = paths
        .iter()
        .map(|p| {
            Ok((
                p.strip_prefix(root())
                    .map_err(|e| e.to_string())?
                    .to_string_lossy()
                    .into_owned(),
                digest(&read(p)?),
            ))
        })
        .collect::<Result<_>>()?;
    Ok(digest(&encode(&identities)?))
}

fn invoke(request: &Value) -> Result<Value> {
    let mut child = Command::new(root().join(".venv/bin/python"))
        .arg(directory().join("tools/reference.py"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .map_err(|e| format!("reference Python: {e}"))?;
    let mut stdin = child.stdin.take().ok_or("reference stdin unavailable")?;
    let bytes = encode(request)?;
    let writer = std::thread::spawn(move || stdin.write_all(&bytes));
    let output = child.wait_with_output().map_err(|e| e.to_string())?;
    let sent = writer.join().map_err(|_| "reference writer panicked")?;
    if !output.status.success() {
        return Err(format!(
            "reference generator {} (see diagnostics above)",
            output.status
        ));
    }
    sent.map_err(|e| e.to_string())?;
    serde_json::from_slice(&output.stdout).map_err(|e| format!("reference output: {e}"))
}

pub(crate) fn generate(
    spec: &Spec,
    cases: &Corpus,
    inputs: &Value,
    threads: usize,
    descriptor_rows: &mut Option<Vec<Value>>,
) -> Result<Vec<Value>> {
    if let Spec::Corpus(task) = *spec {
        if descriptors::handles(task) {
            if descriptor_rows.is_none() {
                *descriptor_rows = Some(
                    serde_json::from_value(invoke(
                        &json!({"kind":"descriptors", "task":spec.key(), "corpus":cases.molecules,"threads":threads}),
                    )?)
                    .map_err(|e| e.to_string())?,
                );
            }
            let rows = descriptor_rows.as_ref().unwrap();
            if rows.len() != cases.molecules.len() {
                return Err("descriptor reference case count mismatch".into());
            }
            let mut records = Vec::new();
            for (case, row) in cases.molecules.iter().zip(rows) {
                if row["smiles"] != case.smiles {
                    return Err("descriptor reference changed input/order".into());
                }
                let registry::Operation::Molecular(id) = task.operation else {
                    return Err("descriptor operation expected".into());
                };
                for profile in id.profiles() {
                    let output = if row["rdkit_ok"] == true {
                        descriptors::observation(profile, row)?
                    } else {
                        crate::molecular::Outcome::Error {
                            stage: crate::molecular::Stage::Parse,
                            detail: "MolFromSmiles returned None".into(),
                        }
                    };
                    records.push(
                        serde_json::to_value(Record {
                            input: Input::Molecular {
                                case: case.clone(),
                                profile,
                            },
                            output: registry::Value::Molecular(output),
                        })
                        .map_err(|e| e.to_string())?,
                    );
                }
            }
            return Ok(records);
        }
    }
    let mut request = match *spec {
        Spec::Special(s) => json!({"kind":s.key,"input":inputs,"threads":threads}),
        Spec::Batch => json!({"kind":"batch_smiles","input":inputs,"threads":threads}),
        Spec::Corpus(task) => {
            let parameters = match task.operation {
                registry::Operation::MolAlign => Ok(Value::Null),
                registry::Operation::PersistentMmff | registry::Operation::PersistentUff => {
                    Ok(Value::Null)
                }
                registry::Operation::Fingerprint(_) => Ok(Value::Null),
                registry::Operation::SmilesWrite => {
                    serde_json::to_value(crate::smiles_write::profiles())
                }
                registry::Operation::Molecular(id) => serde_json::to_value(id.profiles()),
                registry::Operation::SubstructureMatch => {
                    serde_json::to_value(crate::search::profiles())
                }
                op @ (registry::Operation::UffCoverage
                | registry::Operation::UffOptimization
                | registry::Operation::UffConformerOptimization) => {
                    serde_json::to_value(crate::uff::profiles(op))
                }
                registry::Operation::BioPdbOutput => Ok(Value::Null),
                op @ (registry::Operation::MmffCoverage
                | registry::Operation::MmffOptimization
                | registry::Operation::MmffConformerOptimization) => {
                    serde_json::to_value(crate::mmff::profiles(op))
                }
            }
            .map_err(|e| e.to_string())?;
            json!({"kind":"corpus", "generator":task.generator, "input":inputs,
                   "corpus":cases.molecules,"parameters":parameters,"threads":threads})
        }
    };
    request["task"] = json!(spec.key());
    serde_json::from_value(invoke(&request)?).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    #[test]
    fn reference_transport_preserves_native_coordinate_bits() {
        // Actual native MolAlign coordinates from the retained failing rows.
        // Reading them in normal `prepare` builds must match test builds;
        // enabling float_roundtrip only through a dev-dependency is too late.
        for expected in [
            0.9389657496748851_f64,
            0.012999633836427429,
            0.9999155011900349,
            0.43526753128646156,
            1.0457859682280435,
            0.9857534564380883,
            -0.0,
        ] {
            let encoded = serde_json::to_vec(&expected).unwrap();
            let value: serde_json::Value = serde_json::from_slice(&encoded).unwrap();
            assert_eq!(value.as_f64().unwrap().to_bits(), expected.to_bits());
            let persisted = serde_json::to_vec(&value).unwrap();
            let decoded: f64 = serde_json::from_slice(&persisted).unwrap();
            assert_eq!(decoded.to_bits(), expected.to_bits());
        }
    }
}
