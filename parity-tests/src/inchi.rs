//! Public InChI identifiers, keys and SMILES -> InChI -> molecule roundtrips.
use crate::{Result, chemistry, value::Value};
use cosmolkit::{InchiReadParams, InchiWriteParams, Molecule};
use serde::{Deserialize, Serialize};
use serde_json::{Value as Json, json};

pub const TASKS: &[(&str, &str)] = &[
    ("inchi_write", "inchi"),
    ("inchi_key", "inchi"),
    ("inchi_roundtrip", "inchi"),
];

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
            Value::Text(mol.to_inchi_with_params(&InchiWriteParams::new(p.options))?)
        }
        "inchi_key" => {
            let p: WriteProfile = serde_json::from_value(profile.clone())?;
            let params = InchiWriteParams::new(p.options);
            let identifier = mol.to_inchi_with_params(&params)?;
            Value::List(vec![
                Value::Text(mol.to_inchi_key_with_params(&params)?),
                Value::Text(cosmolkit::inchi_to_key(&identifier)?),
            ])
        }
        "inchi_roundtrip" => {
            let p: ReadProfile = serde_json::from_value(profile.clone())?;
            let identifier = mol.to_inchi()?;
            let parsed = Molecule::from_inchi_with_params(
                &identifier,
                &InchiReadParams::new(p.read_sanitize, p.read_remove_hs),
            )?;
            Value::List(vec![Value::Text(identifier), chemistry::state(&parsed)])
        }
        _ => return Err(format!("unknown InChI task {task}").into()),
    })
}
