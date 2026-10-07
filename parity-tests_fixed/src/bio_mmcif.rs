//! Fixed mmCIF output-switch matrix, using the existing BIO corpus inputs.
use crate::{Result, special_regression::Snapshot};
use cosmolkit::{BioMmcifWriteParams, BioStructure};
use serde_json::Value;

// This is parameter projection, not a second implementation of the writer.
// The original native matrix toggles one group from MmcifOutputGroups(all).
macro_rules! switches {
    ($($field:ident),* $(,)?) => {
        pub(crate) const FLAGS: &[&str] = &[$(stringify!($field),)* "auth_all"];
        fn params(all: bool, flag: &str, value: bool) -> Result<BioMmcifWriteParams> {
            if !FLAGS.contains(&flag) {
                return Err(format!("unknown mmCIF output switch: {flag}"));
            }
            Ok(BioMmcifWriteParams {
                $($field: if flag == stringify!($field) { value } else { all },)*
                auth_all: flag == "auth_all" && value,
                // Python Gemmi cannot toggle modres; retain the global default.
                modres: all,
                ..BioMmcifWriteParams::default()
            })
        }
    };
}
switches!(
    atoms,
    block_name,
    entry,
    database_status,
    author,
    cell,
    symmetry,
    entity,
    entity_poly,
    struct_ref,
    chem_comp,
    exptl,
    diffrn,
    reflns,
    refine,
    title_keywords,
    ncs,
    struct_asym,
    origx,
    struct_conf,
    struct_sheet,
    struct_biol,
    assembly,
    conn,
    cis,
    scale,
    atom_type,
    entity_poly_seq,
    tls,
    software,
    group_pdb,
);

pub(crate) fn validate(fixture: &Value, rows: &[Value], count: usize) -> Result<()> {
    let pin: Value = serde_json::from_str(include_str!("../testdata/reference/gemmi.json"))
        .map_err(|e| e.to_string())?;
    let cases = fixture["cases"].as_array().ok_or("missing BIO cases")?;
    if fixture["schema_version"] != 1
        || fixture["reference"]["version"] != pin["version"]
        || fixture["reference"]["source_revision"] != pin["source_revision"]
        || fixture["flags"] != serde_json::json!(FLAGS)
        || cases.len() != 2
        || cases[0]["case_id"] == cases[1]["case_id"]
        || rows.len() != count
        || count != cases.len() * 2 * FLAGS.len()
    {
        return Err("BIO mmCIF fixture/reference identity or census mismatch".into());
    }
    let mut index = 0;
    for case in cases {
        if !case["case_id"].is_string()
            || !case["text"].as_str().is_some_and(|text| !text.is_empty())
            || !matches!(case["format"].as_str(), Some("pdb" | "mmcif"))
        {
            return Err("BIO mmCIF input schema mismatch".into());
        }
        for all in [false, true] {
            for &flag in FLAGS {
                let row = &rows[index];
                if row["case_id"] != case["case_id"]
                    || row["all_groups"] != all
                    || row["flag"] != flag
                    || row["value"] != (flag == "auth_all" || !all)
                    || !row["text"].is_string()
                {
                    return Err(format!(
                        "BIO mmCIF reference row {index} identity/schema mismatch"
                    ));
                }
                index += 1;
            }
        }
    }
    Ok(())
}

pub fn compare(snapshot: &Snapshot) -> Result<()> {
    let mut results = Vec::new();
    let cases = snapshot.fixture["cases"]
        .as_array()
        .ok_or("missing BIO cases")?;
    for row in &snapshot.rows {
        let case = cases
            .iter()
            .find(|case| case["case_id"] == row["case_id"])
            .ok_or("unknown BIO reference case")?;
        let input = case["text"].as_str().ok_or("missing BIO input text")?;
        // Retain the original file-read API and identical filename context.
        // Use the preflight snapshot, not a potentially changed source file.
        let folder = tempfile::tempdir().map_err(|e| e.to_string())?;
        let filename =
            std::path::Path::new(case["input"].as_str().ok_or("missing BIO input path")?)
                .file_name()
                .ok_or("missing BIO input filename")?;
        let path = folder.path().join(filename);
        std::fs::write(&path, input.as_bytes()).map_err(|e| e.to_string())?;
        let structure = BioStructure::read(&path).map_err(|e| e.to_string());
        let params = params(
            row["all_groups"].as_bool().ok_or("missing all_groups")?,
            row["flag"].as_str().ok_or("missing output switch")?,
            row["value"].as_bool().ok_or("missing switch value")?,
        )?;
        let actual = structure.and_then(|structure| {
            structure
                .to_mmcif_with_params(&params)
                .map_err(|e| e.to_string())
        });
        let matched = actual
            .as_ref()
            .is_ok_and(|text| Some(text.as_str()) == row["text"].as_str());
        results.push(serde_json::json!({
            "case_id": row["case_id"], "all_groups": row["all_groups"],
            "flag": row["flag"], "value": row["value"], "matched": matched,
            "expected": row["text"], "actual": actual,
        }));
    }
    let failed = results.iter().filter(|row| row["matched"] != true).count();
    let directory = crate::directory().join("reports");
    std::fs::create_dir_all(&directory).map_err(|e| e.to_string())?;
    let report = tempfile::Builder::new()
        .prefix("bio-mmcif-switches-")
        .tempdir_in(directory)
        .map_err(|e| e.to_string())?
        .keep()
        .join("comparison.json");
    std::fs::write(&report, crate::encode(&serde_json::json!({
        "task": "bio_mmcif_switches", "compared": results.len(), "failed": failed, "rows": results,
    }))?).map_err(|e| e.to_string())?;
    println!(
        "bio_mmcif_switches: {}/{} matched; {}",
        snapshot.rows.len() - failed,
        snapshot.rows.len(),
        report.display()
    );
    if failed != 0 {
        return Err(format!(
            "BIO mmCIF switches: {failed} mismatches; {}",
            report.display()
        ));
    }
    Ok(())
}
