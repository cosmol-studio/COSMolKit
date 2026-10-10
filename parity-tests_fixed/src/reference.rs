use crate::{
    Corpus, Input, Record, Result, descriptors, digest, directory, encode, read, registry, root,
    workflow::Spec,
};
use serde_json::{Value, json};
use std::{
    io::{BufRead, BufReader, BufWriter, Seek, SeekFrom, Write},
    process::{Command, Stdio},
};

pub(crate) fn uses_gemmi(spec: &Spec) -> bool {
    matches!(*spec, Spec::Corpus(t) if t.operation == registry::Operation::BioPdbOutput)
        || matches!(*spec, Spec::Special(s) if matches!(s.schema, registry::SpecialRegressionSchema::BioMmcifSwitches))
}

pub(crate) fn source_digest(spec: &Spec) -> Result<String> {
    source_digest_at(spec, &root())
}

pub(crate) fn source_digest_at(spec: &Spec, checkout: &std::path::Path) -> Result<String> {
    let directory = || checkout.join("parity-tests_fixed");
    let root = || checkout.to_path_buf();
    // JSON decoding is part of preparation: float_roundtrip must preserve the
    // native f64 inputs, including rank-deficient alignment coordinates.
    let mut paths = vec![
        directory().join("tools/reference.py"),
        directory().join("tools/forcefield_preparation.py"),
        directory().join("Cargo.toml"),
    ];
    if let Spec::Special(s) = spec {
        let generators: &[&str] = match s.schema {
            registry::SpecialRegressionSchema::ConformerFixed19 => {
                &["tools/testdata/rdkit/_generate_conformer_generation_golden.py"]
            }
            registry::SpecialRegressionSchema::ConformerLibrary => {
                &["tools/testdata/rdkit/_generate_conformer_generation_library_golden.py"]
            }
            registry::SpecialRegressionSchema::ForcefieldProperties => &[
                "tools/testdata/rdkit/_generate_forcefield_coverage_golden.py",
                "tools/testdata/rdkit/_generate_forcefield_params_golden.py",
            ],
            _ => &[],
        };
        paths.extend(generators.iter().map(|p| root().join(p)));
    }
    if matches!(*spec, Spec::Special(s) if matches!(s.schema, registry::SpecialRegressionSchema::Mcs))
    {
        paths.push(directory().join("tools/mcs.py"));
    }
    if matches!(*spec, Spec::Corpus(t) if matches!(t.operation, registry::Operation::Fingerprint(_)))
    {
        paths.push(directory().join("tools/fingerprints.py"));
    }
    if uses_gemmi(spec) {
        paths.push(directory().join("testdata/reference/gemmi.json"));
    } else {
        if matches!(*spec, Spec::Special(s) if matches!(s.schema,
            registry::SpecialRegressionSchema::ForcefieldOptimizers | registry::SpecialRegressionSchema::MmffBuiltin))
        {
            paths.extend([
                directory().join("tools/forcefield_regression.py"),
                root().join("tools/testdata/rdkit/_generate_forcefield_params_golden.py"),
                root().join("tools/testdata/rdkit/_generate_mmff_builtin_golden.py"),
            ]);
        }
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

#[derive(serde::Serialize)]
struct NativeRequest<'a> {
    #[serde(skip_serializing_if = "Option::is_none")]
    corpus: Option<&'a [registry::SmilesCase]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    generator: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    input: Option<&'a Value>,
    kind: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    parameters: Option<Value>,
    task: &'static str,
    threads: usize,
}

fn invoke_output<T: serde::Serialize + Sync>(request: &T) -> Result<tempfile::NamedTempFile> {
    // The native command/inputs are unchanged. A private file replaces the
    // whole stdout Vec so native failures still precede JSON decoding errors.
    let artifacts = directory().join("reports");
    std::fs::create_dir_all(&artifacts).map_err(|e| e.to_string())?;
    let output = tempfile::NamedTempFile::new_in(artifacts).map_err(|e| e.to_string())?;
    let mut child = Command::new(root().join(".venv/bin/python"))
        .arg(directory().join("tools/reference.py"))
        .stdin(Stdio::piped())
        .stdout(Stdio::from(output.reopen().map_err(|e| e.to_string())?))
        .stderr(Stdio::inherit())
        .spawn()
        .map_err(|e| format!("reference Python: {e}"))?;
    let stdin = child.stdin.take().ok_or("reference stdin unavailable")?;
    std::thread::scope(|scope| -> Result<()> {
        // Borrow the complete input rather than cloning it into a request
        // Value and allocating a second complete encoded request Vec.
        let writer = scope.spawn(move || -> Result<()> {
            let mut stdin = BufWriter::new(stdin);
            serde_json::to_writer(&mut stdin, request).map_err(|e| e.to_string())?;
            stdin.flush().map_err(|e| e.to_string())
        });
        let status = child.wait().map_err(|e| e.to_string())?;
        let sent = writer.join().map_err(|_| "reference writer panicked")?;
        if !status.success() {
            return Err(format!(
                "reference generator {status} (see diagnostics above)"
            ));
        }
        sent
    })?;
    Ok(output)
}

fn invoke<T: serde::Serialize + Sync>(request: &T) -> Result<Value> {
    let output = invoke_output(request)?;
    serde_json::from_reader(BufReader::new(output.reopen().map_err(|e| e.to_string())?))
        .map_err(|e| format!("reference output: {e}"))
}

/// StreamDeserializer<Value> would deserialize the entire top-level array.
/// SeqAccess reads and releases each complete original array element instead.
fn decode_rows<R: std::io::Read + Seek, F: FnMut(Value) -> Result<()>>(
    mut reader: R,
    emit: &mut F,
) -> Result<usize> {
    let start = reader.stream_position().map_err(|e| e.to_string())?;
    let is_array = {
        let mut peek = BufReader::new(&mut reader);
        loop {
            let bytes = peek.fill_buf().map_err(|e| e.to_string())?;
            let whitespace = bytes
                .iter()
                .take_while(|byte| matches!(**byte, b' ' | b'\t' | b'\n' | b'\r'))
                .count();
            if let Some(first) = bytes.get(whitespace) {
                break *first == b'[';
            }
            let count = bytes.len();
            if count == 0 {
                break false;
            }
            peek.consume(count);
        }
    };
    reader
        .seek(SeekFrom::Start(start))
        .map_err(|e| e.to_string())?;
    if !is_array {
        // Preserve original syntax-before-shape errors for invalid top-level
        // responses. Successful native responses always use the streamed array.
        let value: Value =
            serde_json::from_reader(reader).map_err(|e| format!("reference output: {e}"))?;
        return serde_json::from_value::<Vec<Value>>(value)
            .map(|rows| rows.len())
            .map_err(|e| e.to_string());
    }
    struct RowsVisitor<'a, F>(&'a mut F);
    impl<'de, F: FnMut(Value) -> Result<()>> serde::de::Visitor<'de> for RowsVisitor<'_, F> {
        type Value = usize;
        fn expecting(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_str("a sequence")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut sequence: A,
        ) -> std::result::Result<usize, A::Error> {
            let mut count = 0;
            while let Some(row) = sequence.next_element::<Value>()? {
                (self.0)(row).map_err(serde::de::Error::custom)?;
                count += 1;
            }
            Ok(count)
        }
    }
    let mut parser = serde_json::Deserializer::from_reader(reader);
    let count = serde::de::Deserializer::deserialize_seq(&mut parser, RowsVisitor(emit))
        .map_err(|e| format!("reference output: {e}"))?;
    parser.end().map_err(|e| format!("reference output: {e}"))?;
    Ok(count)
}

fn invoke_rows<T: serde::Serialize + Sync, F: FnMut(Value) -> Result<()>>(
    request: &T,
    emit: &mut F,
) -> Result<usize> {
    let output = invoke_output(request)?;
    decode_rows(
        BufReader::new(output.reopen().map_err(|e| e.to_string())?),
        emit,
    )
}

pub(crate) fn generate<F: FnMut(Value) -> Result<()>>(
    spec: &Spec,
    cases: &Corpus,
    inputs: &Value,
    threads: usize,
    descriptor_rows: &mut Option<Vec<Value>>,
    mut emit: F,
) -> Result<usize> {
    if let Spec::Corpus(task) = *spec {
        if descriptors::handles(task) {
            if descriptor_rows.is_none() {
                *descriptor_rows = Some(
                    serde_json::from_value(invoke(&NativeRequest {
                        corpus: Some(&cases.molecules),
                        generator: None,
                        input: None,
                        kind: "descriptors",
                        parameters: None,
                        task: spec.key(),
                        threads,
                    })?)
                    .map_err(|e| e.to_string())?,
                );
            }
            let rows = descriptor_rows.as_ref().unwrap();
            if rows.len() != cases.molecules.len() {
                return Err("descriptor reference case count mismatch".into());
            }
            let mut count = 0;
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
                    emit(
                        serde_json::to_value(Record {
                            input: Input::Molecular {
                                case: case.clone(),
                                profile,
                            },
                            output: registry::Value::Molecular(output),
                        })
                        .map_err(|e| e.to_string())?,
                    )?;
                    count += 1;
                }
            }
            return Ok(count);
        }
    }
    let request = match *spec {
        Spec::Special(s) => NativeRequest {
            corpus: None,
            generator: None,
            input: Some(inputs),
            kind: s.key,
            parameters: None,
            task: spec.key(),
            threads,
        },
        Spec::Batch => NativeRequest {
            corpus: None,
            generator: None,
            input: Some(inputs),
            kind: "batch_smiles",
            parameters: None,
            task: spec.key(),
            threads,
        },
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
            NativeRequest {
                corpus: Some(&cases.molecules),
                generator: Some(task.generator),
                input: Some(inputs),
                kind: "corpus",
                parameters: Some(parameters),
                task: spec.key(),
                threads,
            }
        }
    };
    invoke_rows(&request, &mut emit)
}

#[cfg(test)]
mod tests {
    use super::*;
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
    #[test]
    fn native_array_stream_preserves_rows_float_bits_and_complete_payloads() {
        let values = [
            0.9389657496748851_f64,
            0.012999633836427429,
            0.9999155011900349,
            0.43526753128646156,
            1.0457859682280435,
            0.9857534564380883,
            -0.0,
        ];
        let input: Vec<Value> = values
            .iter()
            .enumerate()
            .map(|(index, value)| {
                json!({"index":index,"native":value,"bits":u64::MAX,
                   "payload":{"masks":[true,false],"text":"full\nrecord"}})
            })
            .collect();
        let mut observed = Vec::new();
        let count = decode_rows(std::io::Cursor::new(encode(&input).unwrap()), &mut |row| {
            observed.push(row);
            Ok(())
        })
        .unwrap();
        assert_eq!(count, input.len());
        assert_eq!(observed, input);
        for (row, value) in observed.iter().zip(values) {
            assert_eq!(row["native"].as_f64().unwrap().to_bits(), value.to_bits());
            assert_eq!(row["bits"].as_u64().unwrap(), u64::MAX);
        }
    }

    #[test]
    fn non_array_native_errors_keep_original_syntax_then_shape_precedence() {
        for raw in [
            b"{\"x\":[}".as_slice(),
            b"{\"x\":[1]}".as_slice(),
            b"null".as_slice(),
            b"true".as_slice(),
            b"".as_slice(),
        ] {
            let original = match serde_json::from_slice::<Value>(raw) {
                Err(error) => format!("reference output: {error}"),
                Ok(value) => serde_json::from_value::<Vec<Value>>(value)
                    .unwrap_err()
                    .to_string(),
            };
            let actual = decode_rows(std::io::Cursor::new(raw), &mut |_| Ok(())).unwrap_err();
            assert_eq!(actual, original);
        }
    }

    #[test]
    fn native_array_stream_rejects_late_syntax_trailing_data_and_sink_failures() {
        for raw in [b"[1,2,]".as_slice(), b"[1,2] false".as_slice()] {
            assert!(decode_rows(std::io::Cursor::new(raw), &mut |_| Ok(())).is_err());
        }
        let error = decode_rows(std::io::Cursor::new(b"[1,2]"), &mut |_| {
            Err("private reference spool write failed".into())
        })
        .unwrap_err();
        assert!(error.contains("private reference spool write failed"));
    }

    #[test]
    fn borrowed_native_request_serializes_exact_original_request_bytes() {
        let cases = vec![registry::SmilesCase {
            id: "one".into(),
            smiles: "CCO".into(),
        }];
        let input = json!([{"native":0.9389657496748851,"bits":u64::MAX}]);
        let request = NativeRequest {
            corpus: Some(&cases),
            generator: Some("generate_test"),
            input: Some(&input),
            kind: "corpus",
            parameters: Some(Value::Null),
            task: "test_task",
            threads: 7,
        };
        assert_eq!(
            encode(&request).unwrap(),
            encode(&json!({
                "kind":"corpus","generator":"generate_test","input":input,
                "corpus":cases,"parameters":null,"task":"test_task","threads":7,
            }))
            .unwrap()
        );
        let request = NativeRequest {
            corpus: None,
            generator: None,
            input: Some(&input),
            kind: "batch_smiles",
            parameters: None,
            task: "batch_smiles",
            threads: 4,
        };
        assert_eq!(
            encode(&request).unwrap(),
            encode(&json!({
                "kind":"batch_smiles","input":input,"task":"batch_smiles","threads":4,
            }))
            .unwrap()
        );
    }
}
