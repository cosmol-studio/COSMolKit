//! Coordinate PDB projection using the public BIO facade and pinned Gemmi.
use crate::{Result, tasks, value::Value};
use cosmolkit::{BioPdbWriteParams, BioStructure};
use rayon::prelude::*;
use serde::Deserialize;
use serde_json::{Value as Json, json};

pub const TASKS: &[(&str, &str)] = &[("bio_pdb_output", "bio")];
pub const GEMMI: &str = "0.7.5";
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    id: String,
    format: String,
    text: String,
}

pub fn profiles() -> Vec<Json> {
    (0..32)
        .map(|mask| {
            json!({
                "gemmi_version":GEMMI,"ter_records":mask&1!=0,"numbered_ter":mask&2!=0,
                "ter_ignores_type":mask&4!=0,"preserve_serial":mask&8!=0,"end_record":mask&16!=0
            })
        })
        .collect()
}
pub fn normalize(p: Json) -> Result<Json> {
    let fields = p.as_object().ok_or("BIO profile must be an object")?;
    if fields.len() != 6
        || p["gemmi_version"] != GEMMI
        || [
            "ter_records",
            "numbered_ter",
            "ter_ignores_type",
            "preserve_serial",
            "end_record",
        ]
        .iter()
        .any(|k| !p[*k].is_boolean())
    {
        return Err("invalid BIO profile or Gemmi version".into());
    }
    Ok(p)
}
pub fn evaluate(
    profiles: &[Json],
    rows: &[String],
    pool: &rayon::ThreadPool,
) -> Result<Vec<Vec<Vec<u8>>>> {
    pool.install(|| {
        rows.par_iter()
            .map(|text| -> Result<Vec<Vec<u8>>> {
                let row: Input = serde_json::from_str(text)?;
                if row.id.is_empty() {
                    return Err("BIO input needs an id".into());
                }
                let structure = match row.format.as_str() {
                    "pdb" => BioStructure::from_pdb(&row.text)
                        .map_err(|e| Box::new(e) as Box<dyn std::error::Error>),
                    "cif" => BioStructure::from_mmcif(&row.text)
                        .map_err(|e| Box::new(e) as Box<dyn std::error::Error>),
                    _ => return Err("BIO format must be pdb or cif".into()),
                };
                Ok(profiles
                    .iter()
                    .map(|p| match &structure {
                        Err(error) => tasks::failed_error(error.as_ref()),
                        Ok(s) => s
                            .to_pdb_with_params(&BioPdbWriteParams {
                                ter_records: p["ter_records"].as_bool().unwrap(),
                                numbered_ter: p["numbered_ter"].as_bool().unwrap(),
                                ter_ignores_type: p["ter_ignores_type"].as_bool().unwrap(),
                                preserve_serial: p["preserve_serial"].as_bool().unwrap(),
                                end_record: p["end_record"].as_bool().unwrap(),
                            })
                            .map(|text| tasks::success(Value::Bytes(text.into_bytes())))
                            .unwrap_or_else(|e| tasks::failed_error(&e)),
                    })
                    .collect())
            })
            .collect()
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bio_fields_and_all_five_switches_are_explicit() {
        let p = profiles();
        assert_eq!(p.len(), 32);
        for profile in &p {
            assert_eq!(normalize(profile.clone()).unwrap(), *profile);
        }
        let rows = include_str!("../testdata/bio/small.jsonl")
            .lines()
            .map(str::to_owned)
            .collect::<Vec<_>>();
        assert_eq!(rows.len(), 2);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        let output = evaluate(&p, &rows, &pool).unwrap();
        assert_eq!(output.len(), 2);
        assert!(output.iter().flatten().all(|v| v[0] == 0));
    }
}
