//! Fixed MMFF/UFF reference matrices. Preparation is separate from comparison.
use crate::{Result, registry::SpecialRegressionSchema, special_regression::Snapshot};
use cosmolkit::{
    MmffConformerOptimizationParams, MmffEvaluationParams, MmffOptimizationParams,
    MmffPropertiesParams, Molecule, UffConformerOptimizationParams, UffEvaluationParams,
    UffOptimizationParams,
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::collections::BTreeSet;

#[derive(Deserialize)]
struct Numerical {
    ok: bool,
    needs_more: Option<i32>,
    energy: Option<f64>,
    gradient: Option<Vec<f64>>,
    coords: Option<Vec<[f64; 3]>>,
    error: Option<String>,
}
#[derive(Deserialize)]
struct Multiple {
    ok: bool,
    initial_coords: Option<Vec<Vec<[f64; 3]>>>,
    conformer_results: Option<Vec<Numerical>>,
    error: Option<String>,
}
#[derive(Deserialize)]
struct Optimizer {
    cxsmiles: String,
    coords: Vec<[f64; 3]>,
    has_all: bool,
    initial: Numerical,
    single: Numerical,
    multi: Multiple,
}
#[derive(Deserialize)]
struct Builtin {
    rdkit_ok: bool,
    num_atoms: Option<usize>,
    has_all: Option<bool>,
    props_ok: bool,
    atom_types: Option<Vec<u8>>,
    error: Option<String>,
}

pub(crate) fn validate(
    fixture: &Value,
    rows: &[Value],
    count: usize,
    schema: SpecialRegressionSchema,
) -> Result<()> {
    let pin: Value = serde_json::from_str(include_str!("../testdata/reference/rdkit.json"))
        .map_err(|e| e.to_string())?;
    if fixture["schema_version"] != 1
        || fixture["reference"]["version"] != pin["python_distribution_version"]
        || fixture["reference"]["source_revision"] != pin["source_revision"]
    {
        return Err("force-field fixture/reference identity mismatch".into());
    }
    let cases = fixture["cases"]
        .as_array()
        .ok_or("missing force-field cases")?;
    if cases.len() != count || rows.len() != count {
        return Err("force-field fixed case census mismatch".into());
    }
    let mut ids = BTreeSet::new();
    for (case, row) in cases.iter().zip(rows) {
        let id = case["case_id"]
            .as_str()
            .ok_or("missing force-field case ID")?;
        if !ids.insert(id) || row["case_id"] != id || row["input"] != *case {
            return Err(format!("force-field case identity/order mismatch: {id}"));
        }
        if matches!(schema, SpecialRegressionSchema::ForcefieldOptimizers) {
            if fixture["comparison"] != "bitexact"
                || fixture["parameters"]
                    != json!({
                        "max_iterations":2,"seed":61453,"nonbonded_threshold":100,
                        "ignore_interfragment_interactions":true,"num_threads":1,
                        "multi_conformers":2,"mmff_variant":"MMFF94"
                    })
                || !matches!(case["forcefield"].as_str(), Some("mmff" | "uff"))
            {
                return Err("force-field optimizer parameters mismatch".into());
            }
            let output: Optimizer =
                serde_json::from_value(row["output"].clone()).map_err(|e| format!("{id}: {e}"))?;
            let n = output.coords.len();
            let valid = |r: &Numerical| {
                r.ok && r.error.is_none() && r.needs_more.is_some() && r.energy.is_some()
            };
            if n == 0
                || (case["forcefield"] == "mmff" && !output.has_all)
                || !valid(&output.initial)
                || !valid(&output.single)
                || output.initial.gradient.as_ref().map(Vec::len) != Some(3 * n)
                || output.single.coords.as_ref().map(Vec::len) != Some(n)
                || !output.multi.ok
                || output.multi.error.is_some()
                || !output
                    .multi
                    .initial_coords
                    .as_ref()
                    .is_some_and(|c| c.len() == 2 && c.iter().all(|c| c.len() == n))
                || !output.multi.conformer_results.as_ref().is_some_and(|r| {
                    r.len() == 2
                        && r.iter()
                            .all(|r| valid(r) && r.coords.as_ref().map(Vec::len) == Some(n))
                })
            {
                return Err(format!("incomplete native force-field result: {id}"));
            }
        } else {
            let output: Builtin =
                serde_json::from_value(row["output"].clone()).map_err(|e| format!("{id}: {e}"))?;
            if (!output.rdkit_ok && output.error.is_none())
                || (output.props_ok && output.atom_types.as_ref().map(Vec::len) != output.num_atoms)
                || (output.rdkit_ok && output.has_all.is_none())
            {
                return Err(format!("incomplete native MMFF type result: {id}"));
            }
        }
    }
    Ok(())
}

fn exact(
    label: &str,
    actual: impl IntoIterator<Item = f64>,
    expected: impl IntoIterator<Item = f64>,
) -> Result<()> {
    let a: Vec<_> = actual.into_iter().map(f64::to_bits).collect();
    let b: Vec<_> = expected.into_iter().map(f64::to_bits).collect();
    if a != b {
        let first = a.iter().zip(&b).position(|(a, b)| a != b);
        return Err(format!(
            "{label} bit mismatch at {first:?}: lengths {}/{}; bits {:?}/{:?}",
            a.len(),
            b.len(),
            first.map(|i| a[i]),
            first.map(|i| b[i])
        ));
    }
    Ok(())
}
fn result(
    label: &str,
    status: i32,
    energy: f64,
    coords: &[[f64; 3]],
    expected: &Numerical,
) -> Result<()> {
    if Some(status) != expected.needs_more {
        return Err(format!(
            "{label} status: {status} != {:?}",
            expected.needs_more
        ));
    }
    exact(
        &format!("{label} energy"),
        [energy],
        [expected.energy.ok_or("missing energy")?],
    )?;
    exact(
        &format!("{label} coordinates"),
        coords.iter().flatten().copied(),
        expected
            .coords
            .as_ref()
            .ok_or("missing coordinates")?
            .iter()
            .flatten()
            .copied(),
    )
}

fn optimizer(row: &Value) -> Result<()> {
    let expected: Optimizer =
        serde_json::from_value(row["output"].clone()).map_err(|e| e.to_string())?;
    let mol = Molecule::from_smiles(&expected.cxsmiles).map_err(|e| e.to_string())?;
    if mol.conformers_3d().len() != 1 {
        return Err("expected one initial conformer".into());
    }
    exact(
        "initial coordinates",
        mol.conformers_3d()[0]
            .coordinates()
            .iter()
            .flatten()
            .copied(),
        expected.coords.iter().flatten().copied(),
    )?;
    let multi = Molecule::from_parts(
        mol.topology().clone(),
        cosmolkit::CoordinateBlock {
            conformers_3d: expected
                .multi
                .initial_coords
                .as_ref()
                .ok_or("missing multi coordinates")?
                .iter()
                .enumerate()
                .map(|(id, coords)| cosmolkit::Conformer3D::new(id, coords.clone(), true))
                .collect(),
            ..Default::default()
        },
        mol.properties().clone(),
    )
    .map_err(|e| e.to_string())?;
    let (
        has_all,
        initial_status,
        initial_energy,
        gradient,
        single_status,
        single_energy,
        single,
        multiple,
    ): (_, _, _, _, _, _, _, Vec<(i32, f64, Vec<[f64; 3]>)>) =
        if row["input"]["forcefield"] == "mmff" {
            let has_all = mol
                .mmff_has_all_molecule_params()
                .map_err(|e| e.to_string())?;
            let params = MmffConformerOptimizationParams {
                num_threads: 1,
                max_iterations: 0,
                mmff_variant: "MMFF94".into(),
                non_bonded_threshold: 100.0,
                ignore_interfragment_interactions: true,
            };
            let initial = mol
                .with_mmff_optimized_conformers_with_params(&params)
                .map_err(|e| e.to_string())?;
            if initial.conformer_results.len() != 1 {
                return Err("initial MMFF conformer count".into());
            }
            let evaluation = mol
                .mmff_energy_gradient_with_params(&MmffEvaluationParams {
                    mmff_variant: "MMFF94".into(),
                    non_bonded_threshold: 100.0,
                    conformer_id: None,
                    ignore_interfragment_interactions: true,
                })
                .map_err(|e| e.to_string())?
                .ok_or("MMFF unavailable")?;
            let single = mol
                .with_mmff_optimized_with_params(&MmffOptimizationParams {
                    mmff_variant: "MMFF94".into(),
                    max_iterations: 2,
                    non_bonded_threshold: 100.0,
                    conformer_id: None,
                    ignore_interfragment_interactions: true,
                })
                .map_err(|e| e.to_string())?;
            let single_energy = mol
                .with_mmff_optimized_conformers_with_params(&MmffConformerOptimizationParams {
                    max_iterations: 2,
                    ..params.clone()
                })
                .map_err(|e| e.to_string())?;
            if single_energy.conformer_results.len() != 1 {
                return Err("single MMFF conformer count".into());
            }
            let multiple = multi
                .with_mmff_optimized_conformers_with_params(&MmffConformerOptimizationParams {
                    max_iterations: 2,
                    ..params
                })
                .map_err(|e| e.to_string())?;
            if multiple.conformer_results.len() != multiple.molecule.conformers_3d().len() {
                return Err("MMFF result/coordinate count".into());
            }
            (
                has_all,
                initial.conformer_results[0].needs_more,
                initial.conformer_results[0].energy,
                evaluation.gradient().to_vec(),
                single.needs_more,
                single_energy.conformer_results[0].energy,
                single.molecule,
                multiple
                    .conformer_results
                    .iter()
                    .zip(multiple.molecule.conformers_3d())
                    .map(|(r, c)| (r.needs_more, r.energy, c.coordinates().to_vec()))
                    .collect(),
            )
        } else {
            let mol = mol.with_assigned_valence().map_err(|e| e.to_string())?;
            let multi = multi.with_assigned_valence().map_err(|e| e.to_string())?;
            let has_all = mol
                .uff_has_all_molecule_params()
                .map_err(|e| e.to_string())?;
            let params = UffConformerOptimizationParams {
                num_threads: 1,
                max_iterations: 0,
                vdw_threshold: 100.0,
                ignore_interfragment_interactions: true,
            };
            let initial = mol
                .with_uff_optimized_conformers_with_params(&params)
                .map_err(|e| e.to_string())?;
            if initial.conformers.len() != 1 {
                return Err("initial UFF conformer count".into());
            }
            let evaluation = mol
                .uff_energy_gradient_with_params(&UffEvaluationParams {
                    vdw_threshold: 100.0,
                    conformer_id: None,
                    ignore_interfragment_interactions: true,
                })
                .map_err(|e| e.to_string())?;
            let single = mol
                .with_uff_optimized_with_params(&UffOptimizationParams {
                    max_iterations: 2,
                    vdw_threshold: 100.0,
                    conformer_id: None,
                    ignore_interfragment_interactions: true,
                })
                .map_err(|e| e.to_string())?;
            let multiple = multi
                .with_uff_optimized_conformers_with_params(&UffConformerOptimizationParams {
                    max_iterations: 2,
                    ..params
                })
                .map_err(|e| e.to_string())?;
            if multiple.conformers.len() != multiple.molecule.conformers_3d().len() {
                return Err("UFF result/coordinate count".into());
            }
            (
                has_all,
                initial.conformers[0].status,
                initial.conformers[0].energy,
                evaluation.gradient().to_vec(),
                single.status,
                single.energy,
                single.molecule,
                multiple
                    .conformers
                    .iter()
                    .zip(multiple.molecule.conformers_3d())
                    .map(|(r, c)| (r.status, r.energy, c.coordinates().to_vec()))
                    .collect(),
            )
        };
    if has_all != expected.has_all || Some(initial_status) != expected.initial.needs_more {
        return Err("availability/initial status mismatch".into());
    }
    let mut checks = vec![
        exact(
            "initial energy",
            [initial_energy],
            [expected.initial.energy.ok_or("missing initial energy")?],
        ),
        exact(
            "initial gradient",
            gradient,
            expected
                .initial
                .gradient
                .ok_or("missing initial gradient")?,
        ),
    ];
    if single.conformers_3d().len() != 1 {
        return Err("single coordinate count".into());
    }
    checks.push(result(
        "single",
        single_status,
        single_energy,
        single.conformers_3d()[0].coordinates(),
        &expected.single,
    ));
    let results = expected
        .multi
        .conformer_results
        .ok_or("missing multi results")?;
    if multiple.len() != results.len() {
        return Err("multi result count mismatch".into());
    }
    for (i, ((status, energy, coords), expected)) in multiple.iter().zip(&results).enumerate() {
        checks.push(result(
            &format!("multi {i}"),
            *status,
            *energy,
            coords,
            expected,
        ));
    }
    let failures: Vec<_> = checks
        .into_iter()
        .filter_map(std::result::Result::err)
        .collect();
    if failures.is_empty() {
        Ok(())
    } else {
        Err(failures.join("; "))
    }
}

fn builtin(row: &Value) -> Result<()> {
    let expected: Builtin =
        serde_json::from_value(row["output"].clone()).map_err(|e| e.to_string())?;
    if !expected.rdkit_ok {
        return Ok(());
    } // Original native parse-error branch, validated in preflight.
    let mol = Molecule::from_smiles(row["input"]["smiles"].as_str().ok_or("missing SMILES")?)
        .map_err(|e| e.to_string())?;
    if Some(
        mol.mmff_has_all_molecule_params()
            .map_err(|e| e.to_string())?,
    ) != expected.has_all
    {
        return Err("MMFF availability mismatch".into());
    }
    if expected.props_ok {
        let props = mol
            .mmff_properties_with_params(&MmffPropertiesParams {
                mmff_variant: row["input"]["variant"]
                    .as_str()
                    .ok_or("missing variant")?
                    .into(),
            })
            .map_err(|e| e.to_string())?;
        let types: Vec<_> = (0..mol.num_atoms())
            .map(|i| props.atom_type(i).map_err(|e| e.to_string()))
            .collect::<Result<_>>()?;
        if Some(mol.num_atoms()) != expected.num_atoms || Some(types) != expected.atom_types {
            return Err("MMFF atom count/types mismatch".into());
        }
    }
    Ok(())
}

fn compare(snapshot: &Snapshot, key: &str, call: fn(&Value) -> Result<()>) {
    let results: Vec<_> = snapshot
        .rows
        .iter()
        .map(|row| {
            let outcome = call(row);
            json!({"case_id":row["case_id"],"matches":outcome.is_ok(),"error":outcome.err()})
        })
        .collect();
    let failures: Vec<_> = results.iter().filter(|r| r["matches"] != true).collect();
    let reports = crate::directory().join("reports");
    std::fs::create_dir_all(&reports).unwrap();
    let report = reports.join(format!("{key}.json"));
    std::fs::write(
        &report,
        serde_json::to_vec_pretty(&json!({"total":results.len(),
        "failures":failures.len(),"results":results}))
        .unwrap(),
    )
    .unwrap();
    assert!(
        failures.is_empty(),
        "{}/{} fixed cases failed; report {}: {:?}",
        failures.len(),
        snapshot.rows.len(),
        report.display(),
        failures.iter().take(5).collect::<Vec<_>>()
    );
}
pub fn compare_optimizers(snapshot: &Snapshot) {
    compare(snapshot, "forcefield_optimizers", optimizer);
}
pub fn compare_builtin(snapshot: &Snapshot) {
    compare(snapshot, "mmff_builtin", builtin);
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn preflight_preserves_false_uff_coverage_and_rejects_missing_gradient() {
        let mut fixture: Value = serde_json::from_str(include_str!(
            "../testdata/special/forcefield_optimizers.json"
        ))
        .unwrap();
        let case = fixture["cases"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["case_id"] == "uff_row_47")
            .unwrap()
            .clone();
        fixture["cases"] = json!([case]);
        // This is a schema probe, not a chemistry expectation. Native UFF can
        // return a field even when HasAllMoleculeParams is false (Se2+2).
        let final_result = json!({"ok":true,"needs_more":1,"energy":1.0,
            "coords":[[0.0,0.0,0.0]],"error":null});
        let mut rows = vec![json!({"case_id":case["case_id"],"input":case,"output":{
            "cxsmiles":"C","coords":[[0.0,0.0,0.0]],"has_all":false,
            "initial":{"ok":true,"needs_more":1,"energy":1.0,"gradient":[0.0,0.0,0.0],"error":null},
            "single":final_result,"multi":{"ok":true,"error":null,
                "initial_coords":[[[0.0,0.0,0.0]],[[0.0,0.0,0.0]]],
                "conformer_results":[final_result,final_result]}
        }})];
        assert!(
            validate(
                &fixture,
                &rows,
                1,
                SpecialRegressionSchema::ForcefieldOptimizers
            )
            .is_ok()
        );
        rows[0]["output"]["initial"]["gradient"] = Value::Null;
        assert!(
            validate(
                &fixture,
                &rows,
                1,
                SpecialRegressionSchema::ForcefieldOptimizers
            )
            .is_err()
        );
        assert!(
            validate(
                &fixture,
                &rows,
                2,
                SpecialRegressionSchema::ForcefieldOptimizers
            )
            .is_err()
        );
    }
    #[test]
    fn comparison_rejects_one_bit_and_signed_zero_differences() {
        assert!(exact("one bit", [1.0], [f64::from_bits(1.0f64.to_bits() + 1)]).is_err());
        assert!(exact("zero", [0.0], [-0.0]).is_err());
        assert!(exact("length", [1.0], []).is_err());
        assert!(exact("equal", [1.0, -0.0], [1.0, -0.0]).is_ok());
    }
}
