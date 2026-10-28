//! Public InChI identifiers, keys and SMILES -> InChI -> molecule roundtrips.
use crate::{Result, chemistry, value::Value};
use cosmolkit::{InchiReadParams, InchiWriteParams, Molecule};
use serde::{Deserialize, Serialize};
use serde_json::{Value as Json, json};

/// Exact engine rejection, distinct from read/sanitization or toolkit errors.
#[derive(Debug)]
pub(crate) struct GenerationRejected(pub i32);

impl std::fmt::Display for GenerationRejected {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "InChI generation rejected with engine status {}", self.0)
    }
}
impl std::error::Error for GenerationRejected {}

pub(crate) fn generation_rejected(identifier: &str, code: i32) -> bool {
    // InChI inchi_api.h: inchi_Ret_ERROR = 2, inchi_Ret_FATAL = 3.
    // No warning, timeout, unknown status or nonempty output is normalized.
    identifier.is_empty() && matches!(code, 2 | 3)
}

fn identifier(mol: &Molecule, params: &InchiWriteParams) -> Result<String> {
    match mol.to_inchi_with_params(params) {
        Ok(identifier) => Ok(identifier),
        Err(error)
            if error.operation == "to_inchi"
                && error.kind == cosmolkit::InchiErrorKind::InvalidInput =>
        {
            // Preserve the scalar API check; recover engine status only on its
            // explicit empty-identifier rejection, never from error text.
            let result = mol.to_inchi_result_with_params(params)?;
            if generation_rejected(result.inchi(), result.return_code()) {
                Err(Box::new(GenerationRejected(result.return_code())))
            } else {
                Err(Box::new(error))
            }
        }
        Err(error) => Err(Box::new(error)),
    }
}

pub const TASKS: &[(&str, &str)] = &[
    ("inchi_write", "inchi"),
    ("inchi_key", "inchi"),
    ("inchi_roundtrip", "inchi"),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generation_rejection_requires_empty_output_and_exact_engine_error_status() {
        for code in [-1, 0, 1, 4, 5] {
            assert!(!generation_rejected("", code));
        }
        for code in [2, 3] {
            assert!(generation_rejected("", code));
            assert!(!generation_rejected("InChI=1S/CH4/h1H4", code));
        }
    }

    #[test]
    fn generation_rejections_compare_by_stage_and_code() {
        let error = crate::tasks::failed_error(&GenerationRejected(2));
        assert_eq!(
            error,
            crate::tasks::outcome(
                1,
                Value::List(vec![
                    Value::Text("InchiGenerationRejected".into()),
                    Value::Int(2),
                ])
            )
        );
        assert_ne!(error, crate::tasks::failed_error(&GenerationRejected(3)));
        assert_ne!(error, crate::tasks::outcome(0, Value::Text(String::new())));
        assert_ne!(
            error,
            crate::tasks::failed_error(&cosmolkit::InchiError {
                operation: "to_inchi",
                kind: cosmolkit::InchiErrorKind::InvalidInput,
                detail: "the InChI engine returned no identifier".into(),
            })
        );
        assert_ne!(
            error,
            crate::tasks::failed_error(&cosmolkit::InchiError {
                operation: "mol_from_inchi",
                kind: cosmolkit::InchiErrorKind::SanitizeFailed,
                detail: "invalid valence".into(),
            })
        );
    }
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct WriteProfile {
    options: String,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct ReadProfile {
    read_sanitize: bool,
    read_remove_hs: bool,
}

pub fn profiles(task: &str) -> Result<Vec<Json>> {
    Ok(match task {
        "inchi_write" | "inchi_key" => ["", "/FixedH", "/SNon"]
            .map(|options| json!({"options": options}))
            .into(),
        "inchi_roundtrip" => [false, true]
            .into_iter()
            .flat_map(|sanitize| {
                [false, true]
                .map(move |remove_hs| json!({"read_sanitize": sanitize, "read_remove_hs": remove_hs}))
            })
            .collect(),
        _ => return Err(format!("unknown InChI task {task}").into()),
    })
}

pub fn normalize(task: &str, profile: Json) -> Result<Json> {
    Ok(match task {
        "inchi_write" | "inchi_key" => {
            serde_json::to_value(serde_json::from_value::<WriteProfile>(profile)?)?
        }
        "inchi_roundtrip" => serde_json::to_value(serde_json::from_value::<ReadProfile>(profile)?)?,
        _ => return Err(format!("unknown InChI task {task}").into()),
    })
}

pub fn scalar(mol: &Molecule, task: &str, profile: &Json) -> Result<Value> {
    Ok(match task {
        "inchi_write" => {
            let p: WriteProfile = serde_json::from_value(profile.clone())?;
            Value::Text(identifier(mol, &InchiWriteParams::new(p.options))?)
        }
        "inchi_key" => {
            let p: WriteProfile = serde_json::from_value(profile.clone())?;
            let params = InchiWriteParams::new(p.options);
            let identifier = identifier(mol, &params)?;
            Value::List(vec![
                Value::Text(mol.to_inchi_key_with_params(&params)?),
                Value::Text(cosmolkit::inchi_to_key(&identifier)?),
            ])
        }
        "inchi_roundtrip" => {
            let p: ReadProfile = serde_json::from_value(profile.clone())?;
            let identifier = identifier(mol, &InchiWriteParams::default())?;
            let parsed = Molecule::from_inchi_with_params(
                &identifier,
                &InchiReadParams::new(p.read_sanitize, p.read_remove_hs),
            )?;
            Value::List(vec![Value::Text(identifier), chemistry::state(&parsed)])
        }
        _ => return Err(format!("unknown InChI task {task}").into()),
    })
}
