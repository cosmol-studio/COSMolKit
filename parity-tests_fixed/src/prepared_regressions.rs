//! Preflight for relocated reference-dependent regressions, before CK calls.
use crate::{Result, registry::SpecialRegressionSchema};
use serde_json::Value;

fn xyz(value: &Value) -> bool {
    value.as_array().is_some_and(|rows| {
        rows.iter().all(|row| {
            row.as_array().is_some_and(|point| {
                point.len() == 3 && point.iter().all(|v| v.as_f64().is_some_and(f64::is_finite))
            })
        })
    })
}

fn optional(value: &Value, predicate: impl FnOnce(&Value) -> bool) -> bool {
    value.is_null() || predicate(value)
}

pub(crate) fn validate(
    fixture: &Value,
    rows: &[Value],
    count: usize,
    schema: SpecialRegressionSchema,
) -> Result<()> {
    let pin: Value = serde_json::from_str(include_str!("../testdata/reference/rdkit.json"))
        .map_err(|e| e.to_string())?;
    let cases = fixture["cases"]
        .as_array()
        .ok_or("missing regression cases")?;
    if fixture["schema_version"] != 1
        || fixture["reference"] != pin
        || cases.len() != count
        || rows.len() != count
    {
        return Err("prepared regression reference identity/census mismatch".into());
    }
    for (index, (case, row)) in cases.iter().zip(rows).enumerate() {
        let valid = match schema {
            SpecialRegressionSchema::ConformerFixed19 => {
                ["case_id", "mode", "source_kind", "source"]
                    .iter()
                    .all(|key| case[*key].is_string() && row[*key] == case[*key])
                    && row["preset"] == case["preset_name"]
                    && row["attrs"] == case["attrs"]
                    && row["rdkit_ok"].is_boolean()
                    && optional(&row["status"], Value::is_i64)
                    && optional(&row["ids"], |v| {
                        v.as_array().is_some_and(|v| v.iter().all(Value::is_i64))
                    })
                    && optional(&row["failure_counts"], |v| {
                        v.as_array().is_some_and(|v| v.iter().all(Value::is_u64))
                    })
                    && optional(&row["conformers"], |v| {
                        v.as_array().is_some_and(|v| v.iter().all(xyz))
                    })
                    && optional(&row["error"], Value::is_string)
                    && (case["source_kind"] != "fixture_mol" || case["mol_block"].is_string())
            }
            SpecialRegressionSchema::ConformerLibrary => {
                case.is_string()
                    && row["smiles"] == *case
                    && row["seed"] == 61453
                    && row["preset"] == "ETKDGv3"
                    && row["max_iterations"] == 3
                    && row["timeout"] == 0
                    && ["rdkit_parse_ok", "rdkit_add_hs_ok", "rdkit_embed_ok"]
                        .iter()
                        .all(|key| row[*key].is_boolean())
                    && optional(&row["status"], Value::is_i64)
                    && optional(&row["coords"], xyz)
                    && optional(&row["error_stage"], Value::is_string)
                    && optional(&row["error"], Value::is_string)
                    && (row["rdkit_embed_ok"] != true
                        || (row["status"] == 0 && row["coords"].is_array()))
            }
            SpecialRegressionSchema::ForcefieldProperties => {
                case.is_string()
                    && row["smiles"] == *case
                    && row["rdkit_ok"].is_boolean()
                    && optional(&row["error"], Value::is_string)
                    && ["uff", "mmff", "uff_explicit_h", "mmff_explicit_h"]
                        .iter()
                        .all(|key| {
                            let result = &row[*key];
                            result["ok"].is_boolean()
                                && optional(&result["has_all"], Value::is_boolean)
                                && optional(&result["error"], Value::is_string)
                                && (result["ok"] != true || result["has_all"].is_boolean())
                                && (!key.starts_with("mmff") || mmff_fields(result))
                        })
            }
            _ => return Err("wrong prepared regression schema".into()),
        };
        if !valid {
            return Err(format!(
                "prepared regression input/result schema mismatch at row {}",
                index + 1
            ));
        }
    }
    Ok(())
}

fn mmff_fields(result: &Value) -> bool {
    let types = &result["atom_types"];
    if types.is_null() {
        return result["formal_charges"].is_null()
            && result["partial_charges"].is_null()
            && result["has_all"] != true;
    }
    let Some(types) = types.as_array() else {
        return false;
    };
    types
        .iter()
        .all(|v| v.as_u64().is_some_and(|v| v <= u8::MAX as u64))
        && ["formal_charges", "partial_charges"].iter().all(|key| {
            result[*key].as_array().is_some_and(|charges| {
                charges.len() == types.len()
                    && charges
                        .iter()
                        .all(|v| v.as_f64().is_some_and(f64::is_finite))
            })
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn mmff_reference_requires_every_type_and_charge_row() {
        let mut result = json!({"has_all":true,"atom_types":[1],
            "formal_charges":[0.0],"partial_charges":[0.0]});
        assert!(mmff_fields(&result));
        result["partial_charges"] = json!([]);
        assert!(!mmff_fields(&result));
        result["atom_types"] = Value::Null;
        assert!(!mmff_fields(&result));
    }

    #[test]
    fn coordinates_require_finite_three_component_rows() {
        assert!(xyz(&json!([[1.0, 2.0, 3.0]])));
        assert!(!xyz(&json!([[1.0, 2.0]])));
        assert!(!xyz(&json!([[null, 2.0, 3.0]])));
    }
}
