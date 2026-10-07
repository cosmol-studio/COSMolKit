use crate::{Result, registry};
use serde_json::Value;
use std::collections::BTreeSet;
pub struct Snapshot {
    pub fixture: Value,
    pub rows: Vec<Value>,
}
pub(crate) fn validate(
    input: &[u8],
    reference: &[u8],
    count: usize,
    schema: registry::SpecialRegressionSchema,
) -> Result<Snapshot> {
    let fixture: Value = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    // These are distinct pinned spellings, not a version-normalization rule:
    // rdBase reports 2026.03.1; Python distribution/fixture reports 2026.3.1.
    let pin: Value = serde_json::from_str(include_str!("../testdata/reference/rdkit.json"))
        .map_err(|e| e.to_string())?;
    if fixture["schema_version"] != 1
        || fixture["reference"]["version"] != pin["python_distribution_version"]
        || (matches!(schema, registry::SpecialRegressionSchema::TautomerBranches)
            && fixture["reference"]["source_revision"] != pin["source_revision"])
    {
        return Err("special regression fixture schema/reference mismatch".into());
    }
    let mut ids = Vec::new();
    let tables: &[&str] = match schema {
        registry::SpecialRegressionSchema::StructureTags => &["cases", "octahedral_switch_cases"],
        registry::SpecialRegressionSchema::TautomerBranches => &["cases"],
    };
    for &field in tables {
        let cases = fixture[field]
            .as_array()
            .ok_or_else(|| format!("missing fixture table: {field}"))?;
        for case in cases {
            ids.push(
                case["case_id"]
                    .as_str()
                    .ok_or("missing fixture case id")?
                    .to_owned(),
            );
        }
    }
    if ids.len() != count || ids.iter().collect::<BTreeSet<_>>().len() != count {
        return Err("special regression case census/uniqueness mismatch".into());
    }
    let reference = std::str::from_utf8(reference).map_err(|e| e.to_string())?;
    let rows: Vec<Value> = reference
        .lines()
        .map(|line| serde_json::from_str(line).map_err(|e| e.to_string()))
        .collect::<Result<_>>()?;
    if rows.len() != count {
        return Err("special regression reference row count mismatch".into());
    }
    for (id, row) in ids.iter().zip(&rows) {
        if matches!(schema, registry::SpecialRegressionSchema::TautomerBranches) {
            validate_tautomer_row(&fixture, id, row)?;
            continue;
        }
        if row["case_id"].as_str() != Some(id)
            || !matches!(row["status"].as_str(), Some("ok" | "error"))
            || !row["before"].is_object()
            || !row["after"].is_object()
            || !row["environment"].is_object()
            || !row["conf_id"].is_i64()
            || !row["replace_existing_tags"].is_boolean()
        {
            return Err(format!(
                "special regression reference identity/schema: {id}"
            ));
        }
    }
    Ok(Snapshot { fixture, rows })
}

fn validate_tautomer_row(fixture: &Value, id: &str, row: &Value) -> Result<()> {
    let fail = || format!("special regression tautomer identity/schema: {id}");
    let case = fixture["cases"]
        .as_array()
        .ok_or_else(fail)?
        .iter()
        .find(|case| case["case_id"] == id)
        .ok_or_else(fail)?;
    if row["case_id"] != id
        || row["schema_version"] != 1
        || !case["smiles"].is_string()
        || !case["row"].is_u64()
        || !case["sanitize"].is_boolean()
        || !case["remove_hs"].is_boolean()
        || ["row", "smiles", "sanitize", "remove_hs", "source"]
            .iter()
            .any(|key| row[key] != case[key])
        || row["parse"]["ok"] != true
        || !row["parse"]["error"].is_null()
    {
        return Err(fail());
    }
    let branches = row["branches"].as_object().ok_or_else(fail)?;
    let parameters = fixture["branches"].as_array().ok_or_else(fail)?;
    if parameters.len() != 2
        || branches.len() != 2
        || parameters[0]["name"] != "default"
        || parameters[1]["name"] != "v1"
    {
        return Err(fail());
    }
    for params in parameters {
        let name = params["name"].as_str().ok_or_else(fail)?;
        let branch = branches.get(name).ok_or_else(fail)?;
        let smiles = branch["ordered_smiles"].as_array().ok_or_else(fail)?;
        let states = branch["molecule_states"].as_array().ok_or_else(fail)?;
        let scores = branch["scores"].as_array().ok_or_else(fail)?;
        if branch["parameters"] != *params
            || branch["ok"] != true
            || !branch["error"].is_null()
            || !branch["status"].is_string()
            || smiles.is_empty()
            || smiles.iter().any(|v| !v.is_string())
            || states.len() != smiles.len()
            || scores.len() != smiles.len()
            || !branch["modified_atoms"].is_array()
            || !branch["modified_bonds"].is_array()
            || !branch["canonical_smiles"].is_string()
            || branch["canonical_state"]["isomeric_smiles"] != branch["canonical_smiles"]
        {
            return Err(fail());
        }
        for (state, expected_smiles) in states.iter().zip(smiles).chain(std::iter::once((
            &branch["canonical_state"],
            &branch["canonical_smiles"],
        ))) {
            if state["isomeric_smiles"] != *expected_smiles
                || !state["atoms"].is_array()
                || !state["bonds"].is_array()
            {
                return Err(fail());
            }
            for atom in state["atoms"].as_array().unwrap() {
                if [
                    "atomic_number",
                    "formal_charge",
                    "explicit_hydrogens",
                    "isotope",
                    "radical_electrons",
                ]
                .iter()
                .any(|key| !atom[key].is_i64())
                    || ["no_implicit", "aromatic"]
                        .iter()
                        .any(|key| !atom[key].is_boolean())
                    || ["chiral_tag", "hybridization"]
                        .iter()
                        .any(|key| !atom[key].is_string())
                    || !(atom["cip_code"].is_string()
                        || (atom.get("cip_code").is_some() && atom["cip_code"].is_null()))
                {
                    return Err(fail());
                }
            }
            for bond in state["bonds"].as_array().unwrap() {
                if ["begin", "end"].iter().any(|key| !bond[key].is_u64())
                    || ["aromatic", "conjugated"]
                        .iter()
                        .any(|key| !bond[key].is_boolean())
                    || ["bond_type", "direction", "stereo"]
                        .iter()
                        .any(|key| !bond[key].is_string())
                    || !bond["stereo_atoms"]
                        .as_array()
                        .is_some_and(|ids| ids.iter().all(Value::is_u64))
                {
                    return Err(fail());
                }
            }
        }
        if scores.iter().any(|score| {
            ["ring", "substructure", "hetero_hydrogen", "total"]
                .iter()
                .any(|key| !score[key].is_i64())
        }) {
            return Err(fail());
        }
    }
    Ok(())
}

pub fn preflight(key: &str, _data: &std::path::Path) -> Result<Snapshot> {
    let snapshot = crate::testing::special_snapshot(key)?;
    let spec = registry::SPECIAL_REGRESSIONS
        .iter()
        .find(|s| s.key == key)
        .ok_or("unknown special regression")?;
    validate(
        &crate::encode(&snapshot.inputs)?,
        &crate::workflow::jsonl(&snapshot.rows)?,
        spec.rows,
        spec.schema,
    )
}
