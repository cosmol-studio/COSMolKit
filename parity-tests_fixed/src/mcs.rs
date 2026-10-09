//! Designated pinned FMCS regressions. Chemistry stays in the public facade.
use crate::{Result, special_regression::Snapshot};
use cosmolkit::{McsAtomComparator, McsBondComparator, McsParameters, Molecule};
use serde_json::{Value, json};
use std::collections::BTreeSet;

fn array<'a>(value: &'a Value, key: &str) -> Result<&'a [Value]> {
    value[key]
        .as_array()
        .map(Vec::as_slice)
        .ok_or_else(|| format!("missing MCS table: {key}"))
}

fn parameters(value: &Value) -> Result<McsParameters> {
    let mut params = McsParameters::default();
    macro_rules! fields {
        ($target:expr, $values:expr, $($field:ident),+ $(,)?) => {
            for (name, value) in $values.as_object().ok_or("MCS options must be an object")? {
                match name.as_str() {
                    $(stringify!($field) => $target.$field = serde_json::from_value(value.clone())
                        .map_err(|e| format!("invalid MCS option {name}: {e}"))?,)+
                    _ => return Err(format!("unknown MCS option: {name}")),
                }
            }
        };
    }
    for (name, value) in value.as_object().ok_or("missing MCS parameters")? {
        match name.as_str() {
            "atom_comparator" => {
                params.atom_comparator = match value.as_str() {
                    Some("any") => McsAtomComparator::AtomCompareAny,
                    Some("elements") => McsAtomComparator::AtomCompareElements,
                    Some("isotopes") => McsAtomComparator::AtomCompareIsotopes,
                    Some("any_heavy_atom") => McsAtomComparator::AtomCompareAnyHeavyAtom,
                    _ => return Err("invalid MCS atom comparator".into()),
                }
            }
            "bond_comparator" => {
                params.bond_comparator = match value.as_str() {
                    Some("any") => McsBondComparator::BondCompareAny,
                    Some("order") => McsBondComparator::BondCompareOrder,
                    Some("order_exact") => McsBondComparator::BondCompareOrderExact,
                    _ => return Err("invalid MCS bond comparator".into()),
                }
            }
            "atom_compare_parameters" => fields!(
                params.atom_compare_parameters,
                value,
                match_valences,
                match_chiral_tag,
                match_formal_charge,
                ring_matches_ring_only,
                complete_rings_only,
                match_isotope,
                max_distance
            ),
            "bond_compare_parameters" => fields!(
                params.bond_compare_parameters,
                value,
                ring_matches_ring_only,
                complete_rings_only,
                match_fused_rings,
                match_fused_rings_strict,
                match_stereo
            ),
            _ => fields!(
                params,
                json!({name: value}),
                store_all,
                maximize_bonds,
                threshold,
                timeout,
                verbose,
                initial_seed
            ),
        }
    }
    if params.timeout != 30 {
        return Err("designated MCS calls require the recorded 30-second timeout".into());
    }
    Ok(params)
}

pub(crate) fn validate(fixture: &Value, rows: &[Value], count: usize) -> Result<()> {
    let pin: Value = serde_json::from_str(include_str!("../testdata/reference/rdkit.json"))
        .map_err(|e| e.to_string())?;
    if fixture["schema_version"] != 1 || fixture["reference"] != pin {
        return Err("MCS fixture schema/reference mismatch".into());
    }
    let molecules = array(fixture, "molecules")?;
    for recipe in molecules {
        if !matches!(recipe["format"].as_str(), Some("smiles" | "mol"))
            || !recipe["text"].is_string()
            || !recipe["sanitize"].is_boolean()
            || !recipe["remove_hs"].is_boolean()
        {
            return Err("invalid frozen MCS molecule".into());
        }
    }
    let cases = array(fixture, "cases")?;
    if cases.len() != count || rows.len() != count {
        return Err("MCS reference/fixture case count mismatch".into());
    }
    let mut ids = BTreeSet::new();
    let mut pairs = BTreeSet::new();
    for (case, row) in cases.iter().zip(rows) {
        let id = case["case_id"].as_str().ok_or("missing MCS case ID")?;
        if !ids.insert(id) {
            return Err(format!("duplicate MCS case: {id}"));
        }
        let inputs = array(case, "inputs")?;
        if inputs.len() < 2
            || inputs
                .iter()
                .any(|i| i.as_u64().is_none_or(|i| i >= molecules.len() as u64))
        {
            return Err(format!("invalid MCS input indices: {id}"));
        }
        parameters(&case["parameters"])?;
        if ["case_id", "inputs", "parameters", "source"]
            .iter()
            .any(|key| row[key] != case[key])
        {
            return Err(format!("MCS reference identity mismatch: {id}"));
        }
        match row["status"].as_str() {
            Some("ok") => {
                let result = &row["result"];
                if !result["atom_count"].is_u64()
                    || !result["bond_count"].is_u64()
                    || !result["completed"].is_boolean()
                    || !result["smarts"].is_string()
                    || array(result, "degenerate")?.iter().any(|s| !s.is_string())
                {
                    return Err(format!("invalid MCS reference result: {id}"));
                }
                if !result["query"].is_null() {
                    let query = &result["query"];
                    if query["atom_count"] != result["atom_count"]
                        || query["bond_count"] != result["bond_count"]
                        || !query["smarts"].is_string()
                        || array(result, "query_matches")?.len() != inputs.len()
                        || array(result, "query_matches")?
                            .iter()
                            .any(|v| !v.is_boolean())
                    {
                        return Err(format!("invalid MCS query reference: {id}"));
                    }
                } else if !result["query_matches"].is_null() {
                    return Err(format!("MCS matches without a query: {id}"));
                }
            }
            Some("error")
                if row["error"]["type"].is_string() && row["error"]["message"].is_string() =>
            {
                ()
            }
            _ => return Err(format!("invalid MCS reference status: {id}")),
        }
        if count == 210 {
            if inputs.len() != 2 {
                return Err("JNK1 case must be a pair".into());
            }
            pairs.insert((inputs[0].as_u64().unwrap(), inputs[1].as_u64().unwrap()));
        }
    }
    if count == 210
        && (molecules.len() != 21
            || pairs
                != (0..21)
                    .flat_map(|a| (a + 1..21).map(move |b| (a, b)))
                    .collect())
    {
        return Err("JNK1 must retain every one of the 210 ordered-index pairs".into());
    }
    Ok(())
}

fn molecule(recipe: &Value) -> Result<Molecule> {
    let text = recipe["text"]
        .as_str()
        .ok_or("missing frozen molecule text")?;
    let sanitize = recipe["sanitize"].as_bool().ok_or("missing sanitize")?;
    let remove_hs = recipe["remove_hs"].as_bool().ok_or("missing remove_hs")?;
    match recipe["format"].as_str() {
        Some("smiles") => Molecule::from_smiles_with_params(
            text,
            &cosmolkit::SmilesParseParams {
                sanitize,
                remove_hs,
                ..Default::default()
            },
        )
        .map_err(|e| e.to_string()),
        Some("mol") => Molecule::from_mol_with_params(
            text,
            &cosmolkit::SdfReadParams {
                sanitize,
                remove_hs,
                ..Default::default()
            },
        )
        .map_err(|e| e.to_string()),
        _ => Err("unknown frozen molecule format".into()),
    }
}

fn text(value: &cosmolkit::PropertyText) -> Result<&str> {
    std::str::from_utf8(value.as_bytes()).map_err(|e| e.to_string())
}

fn observe(fixture: &Value, case: &Value) -> Result<Value> {
    let recipes = array(fixture, "molecules")?;
    let molecules = array(case, "inputs")?
        .iter()
        .map(|index| molecule(&recipes[index.as_u64().unwrap() as usize]))
        .collect::<Result<Vec<_>>>()?;
    let before = molecules
        .iter()
        .map(|mol| mol.to_binary().map_err(|e| e.to_string()))
        .collect::<Result<Vec<_>>>()?;
    let params = parameters(&case["parameters"])?;
    let refs = molecules.iter().collect::<Vec<_>>();
    let result = cosmolkit::maximum_common_substructure_with_params(&refs, &params);
    let after = molecules
        .iter()
        .map(|mol| mol.to_binary().map_err(|e| e.to_string()))
        .collect::<Result<Vec<_>>>()?;
    if before != after {
        return Err("MCS modified an input molecule".into());
    }
    let result = result.map_err(|e| e.to_string())?;
    let (query, matches) = if let Some(query) = &result.query {
        let smarts = cosmolkit::write_smarts(query, &cosmolkit::SmartsWriteParams::default())
            .map_err(|e| e.to_string())?;
        (
            json!({"atom_count":query.num_atoms(), "bond_count":query.num_bonds(), "smarts":text(&smarts)?}),
            json!(
                molecules
                    .iter()
                    .map(|mol| mol.has_substruct_match(query).map_err(|e| e.to_string()))
                    .collect::<Result<Vec<_>>>()?
            ),
        )
    } else {
        (Value::Null, Value::Null)
    };
    let mut degenerate = result
        .degenerate
        .keys()
        .map(|key| text(key).map(str::to_owned))
        .collect::<Result<Vec<_>>>()?;
    degenerate.sort(); // Map-key order only; never rewrite SMARTS or alternatives.
    Ok(
        json!({"atom_count":result.atom_count, "bond_count":result.bond_count,
        "completed":result.completed, "smarts":text(&result.smarts)?,
        "degenerate":degenerate, "query":query, "query_matches":matches}),
    )
}

pub fn compare_case(key: &str, snapshot: &Snapshot, id: &str) -> Result<()> {
    if !matches!(key, "mcs_upstream" | "mcs_jnk1") {
        return Err("unknown designated MCS regression".into());
    }
    let cases = array(&snapshot.fixture, "cases")?;
    let index = cases
        .iter()
        .position(|case| case["case_id"] == id)
        .ok_or("unknown MCS case ID")?;
    let case = &cases[index];
    let expected = &snapshot.rows[index];
    let actual = match observe(&snapshot.fixture, case) {
        Ok(result) => json!({"status":"ok", "result":result}),
        Err(error) => json!({"status":"error", "error":error}),
    };
    let passed = expected["status"] == "ok"
        && actual["status"] == "ok"
        && expected["result"]["completed"] == true
        && actual["result"]["completed"] == true
        && expected["result"] == actual["result"];
    let folder = crate::directory().join("reports").join(key);
    std::fs::create_dir_all(&folder).map_err(|e| e.to_string())?;
    let path = folder.join(format!("{id}.json"));
    let report = json!({"case":case, "expected":expected, "actual":actual, "passed":passed});
    std::fs::write(
        &path,
        serde_json::to_vec_pretty(&report).map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    if passed {
        Ok(())
    } else {
        Err(format!(
            "{key}/{id}: exact MCS comparison failed; {}",
            path.display()
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mcs_parameter_projection_and_rejections() {
        let params = parameters(&json!({"timeout":30,"store_all":true,"atom_comparator":"isotopes",
            "bond_comparator":"order_exact","atom_compare_parameters":{"max_distance":3.0,"match_formal_charge":true},
            "bond_compare_parameters":{"match_stereo":true}})).unwrap();
        assert!(params.store_all && params.atom_compare_parameters.match_formal_charge);
        assert_eq!(
            params.atom_comparator,
            McsAtomComparator::AtomCompareIsotopes
        );
        assert_eq!(
            params.bond_comparator,
            McsBondComparator::BondCompareOrderExact
        );
        assert_eq!(params.atom_compare_parameters.max_distance, 3.0);
        assert!(params.bond_compare_parameters.match_stereo);
        for value in [
            json!({"timeout":30,"unknown":true}),
            json!({"timeout":30,"atom_comparator":"typo"}),
            json!({"timeout":30,"atom_compare_parameters":{"max_distance":true}}),
            json!({"timeout":0}),
        ] {
            assert!(parameters(&value).is_err(), "{value}");
        }
    }

    #[test]
    fn mcs_preflight_rejects_missing_or_reidentified_rows() {
        let fixture: Value =
            serde_json::from_str(include_str!("../testdata/special/mcs_upstream.json")).unwrap();
        let rows: Vec<_> = fixture["cases"]
            .as_array()
            .unwrap()
            .iter()
            .map(|case| {
                let mut row = case.clone();
                row["status"] = json!("error");
                row["error"] = json!({"type":"ExampleError","message":"retained, never a match"});
                row
            })
            .collect();
        validate(&fixture, &rows, 44).unwrap();
        assert!(validate(&fixture, &rows[..43], 44).is_err());
        let mut changed = rows;
        changed[0]["parameters"]["timeout"] = json!(31);
        assert!(validate(&fixture, &changed, 44).is_err());
    }

    #[test]
    fn mcs_jnk1_preflight_requires_the_complete_pair_matrix() {
        let mut fixture: Value =
            serde_json::from_str(include_str!("../testdata/special/mcs_jnk1.json")).unwrap();
        let rows: Vec<_> = fixture["cases"]
            .as_array()
            .unwrap()
            .iter()
            .map(|case| {
                let mut row = case.clone();
                row["status"] = json!("error");
                row["error"] = json!({"type":"ExampleError","message":"not a passing comparison"});
                row
            })
            .collect();
        validate(&fixture, &rows, 210).unwrap();
        fixture["cases"][0]["inputs"] = fixture["cases"][1]["inputs"].clone();
        let mut changed = rows;
        changed[0]["inputs"] = fixture["cases"][0]["inputs"].clone();
        assert!(
            validate(&fixture, &changed, 210)
                .unwrap_err()
                .contains("210")
        );
    }
}
