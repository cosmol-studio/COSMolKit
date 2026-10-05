use super::registry::{FingerprintValue, Input, Operation, Record, Value, Width};
use cosmolkit::{SparseCountFingerprint, SparseCountFingerprint32};

pub fn run(input: &Input) -> Result<Record, String> {
    if let Input::Search(row) = input {
        return crate::search::run(row);
    }
    if matches!(input, Input::Uff(_)) {
        return crate::uff::run(input);
    }
    let original = input;
    let Input::Fingerprint(input) = input else {
        if let Input::BioPdbOutput { case, profile } = original {
            let structure = match case.format {
                crate::registry::BioPdbCorpusFormat::Cif => {
                    cosmolkit::BioStructure::from_mmcif(&case.text)
                        .map_err(|error| error.to_string())
                }
                crate::registry::BioPdbCorpusFormat::Pdb => {
                    cosmolkit::BioStructure::from_pdb(&case.text).map_err(|error| error.to_string())
                }
            };
            let structure = match structure {
                Ok(structure) => structure,
                Err(error) => {
                    return Ok(Record {
                        input: original.clone(),
                        output: Value::BioPdbOutput(crate::registry::BioPdbOutputValue {
                            text: String::new(),
                            error: Some(crate::registry::BioPdbOutputError::Parse {
                                format: case.format,
                                message: error.to_string(),
                            }),
                        }),
                    });
                }
            };
            let params = cosmolkit::BioPdbWriteParams {
                ter_records: profile.ter_records,
                numbered_ter: profile.numbered_ter,
                ter_ignores_type: profile.ter_ignores_type,
                preserve_serial: profile.preserve_serial,
                end_record: profile.end_record,
            };
            let (text, error) = match structure.to_pdb_with_params(&params) {
                Ok(t) => (t, None),
                Err(e) => (
                    String::new(),
                    Some(crate::registry::BioPdbOutputError::Write {
                        message: e.to_string(),
                    }),
                ),
            };
            return Ok(Record {
                input: original.clone(),
                output: Value::BioPdbOutput(super::registry::BioPdbOutputValue { text, error }),
            });
        }
        return crate::molecular::run(input);
    };
    macro_rules! execute {
        ($ty:ty, $key:ty) => {{
            let build = |entries: &[(u64, i32)]| -> Result<$ty, String> {
                let mut v = <$ty>::new(input.case.length as $key);
                for &(key, count) in entries {
                    if count == 0 {
                        v.set_value(key as $key, 1).map_err(|e| e.to_string())?;
                    }
                }
                // Materialize stored zeros first; never offset nonzero counts,
                // so i32 extrema need no arbitrary corpus restriction.
                v = v.with_added_scalar(-1).map_err(|e| e.to_string())?;
                for &(key, count) in entries {
                    if count != 0 {
                        v.set_value(key as $key, count).map_err(|e| e.to_string())?;
                    }
                }
                Ok(v)
            };
            let left = build(&input.case.left)?;
            let right = build(&input.case.right)?;
            let before = (left.clone(), right.clone());
            let result = match input.operation {
                Operation::FuzzyAnd => left.fuzzy_and(&right),
                Operation::FuzzyOr => left.fuzzy_or(&right),
                Operation::BioPdbOutput => {
                    return Err("bio_pdb_output in fingerprint input".into());
                }
                Operation::SubstructureMatch
                | Operation::Molecular(_)
                | Operation::UffCoverage
                | Operation::UffOptimization
                | Operation::UffConformerOptimization => {
                    return Err("molecular operation in fingerprint input".into());
                }
            }
            .map_err(|e| e.to_string())?;
            if (left, right) != before {
                return Err("operands changed".into());
            }
            FingerprintValue {
                length: result.length() as u64,
                entries: result
                    .nonzero_elements()
                    .iter()
                    .map(|(&k, &v)| (k as u64, v))
                    .collect(),
            }
        }};
    }
    let output = match input.width {
        Width::U32 => execute!(SparseCountFingerprint32, u32),
        Width::U64 => execute!(SparseCountFingerprint, u64),
    };
    Ok(Record {
        input: original.clone(),
        output: Value::Fingerprint(output),
    })
}
