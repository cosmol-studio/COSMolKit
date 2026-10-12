//! Public geometry/numerical calls. Shared coordinates use exact dyadic values.
use crate::{Result, tasks, value::Value};
use cosmolkit::*;
use rayon::prelude::*;
use serde_json::{Value as Json, json};

pub const TASKS: &[(&str, &str)] = &[
    ("coordinates_2d", "depiction"),
    ("svg", "depiction"),
    ("bounds_matrix", "conformer"),
    ("embedding", "conformer"),
    ("alignment", "alignment"),
    ("uff_available", "forcefields"),
    ("mmff_available", "forcefields"),
    ("uff_optimize", "forcefields"),
    ("mmff_optimize", "forcefields"),
    ("uff_optimize_conformers", "forcefields"),
    ("mmff_optimize_conformers", "forcefields"),
    ("uff_force_field", "forcefields"),
    ("mmff_force_field", "forcefields"),
];

pub fn profiles(task: &str) -> Result<Vec<Json>> {
    let mut rows:Vec<Json>=match task {
        "coordinates_2d" => [false,true].map(|canonical|json!({"canonical":canonical,"samples":0,"flips":0,"seed":0,"permute":false,"ring_templates":false})).into(),
        "svg" => vec![json!({"width":300,"height":300,"no_freetype":true})],
        "bounds_matrix" => vec![json!({"add_hs":false}),json!({"add_hs":true})],
        "embedding" => ["dg","etdg","etkdg","etkdg_v2","etkdg_v3"].into_iter().flat_map(|method|[1,2].map(move|count|json!({"method":method,"count":count,"seed":20260306,"add_hs":true,"random_coords":false,"max_iterations":100,"prune_rms":-1.0,"timeout":10}))).collect(),
        "alignment" => ["transform","best","coordinate","align","all","conformers"].into_iter().flat_map(|operation|[false,true].map(move|weighted|json!({"operation":operation,"weighted":weighted,"reflect":false,"explicit_map":true,"symmetrize":true,"max_iterations":50,"max_matches":1000,"seed":17}))).collect(),
        "uff_available"|"mmff_available" => [false,true].map(|add_hs|json!({"add_hs":add_hs})).into(),
        "uff_optimize"|"uff_optimize_conformers"|"uff_force_field" => vec![json!({"variant":"UFF","iterations":2,"threshold":100.0,"ignore_interfragment":true,"add_hs":true,"seed":17,"fixed_first":false,"force_tolerance":1e-4,"energy_tolerance":1e-6})],
        "mmff_optimize"|"mmff_optimize_conformers"|"mmff_force_field" => ["MMFF94","MMFF94s"].map(|variant|json!({"variant":variant,"iterations":2,"threshold":100.0,"ignore_interfragment":true,"add_hs":true,"seed":17,"fixed_first":false,"force_tolerance":1e-4,"energy_tolerance":1e-6})).into(),
        _ => return Err(format!("unknown geometry task {task}").into()),
    };
    if task.ends_with("force_field") {
        for p in &mut rows {
            p["coordinate_update"] = json!("initialize_distances");
        }
    }
    if matches!(task, "mmff_optimize" | "mmff_optimize_conformers") {
        for p in &mut rows {
            p["optimizer_observation"] = json!("isolated_queries_v1");
        }
    }
    Ok(rows)
}
pub fn normalize(task: &str, p: Json) -> Result<Json> {
    let defaults = profiles(task)?;
    if task == "svg" && p["no_freetype"] != true {
        return Err("CK's SVG backend uses builtin fonts (no_freetype=true)".into());
    }
    if task.ends_with("force_field") && p["coordinate_update"] != "initialize_distances" {
        return Err("invalid native coordinate-update policy".into());
    }
    if matches!(task, "mmff_optimize" | "mmff_optimize_conformers")
        && p["optimizer_observation"] != "isolated_queries_v1"
    {
        return Err("invalid MMFF optimizer observation policy".into());
    }
    let object = p.as_object().ok_or("geometry profile must be an object")?;
    let schema = defaults[0].as_object().unwrap();
    if object.len() != schema.len()
        || object.iter().any(|(k, v)| {
            schema.get(k).is_none_or(|expected| {
                if expected.is_boolean() {
                    !v.is_boolean()
                } else if expected.is_string() {
                    !v.is_string()
                } else {
                    !v.is_number()
                }
            })
        })
    {
        return Err(format!("invalid parameters for {task}").into());
    }
    if task == "embedding"
        && (!matches!(
            p["method"].as_str(),
            Some("dg" | "etdg" | "etkdg" | "etkdg_v2" | "etkdg_v3")
        ) || p["count"].as_u64().is_none_or(|n| n == 0 || n > 64))
    {
        return Err("invalid embedding method/count".into());
    }
    if task == "alignment"
        && !matches!(
            p["operation"].as_str(),
            Some("transform" | "best" | "coordinate" | "align" | "all" | "conformers")
        )
    {
        return Err("invalid alignment operation".into());
    }
    if task.starts_with("mmff_")
        && task != "mmff_available"
        && !matches!(p["variant"].as_str(), Some("MMFF94" | "MMFF94s"))
    {
        return Err("invalid MMFF variant".into());
    }
    for key in [
        "width",
        "height",
        "iterations",
        "max_iterations",
        "max_matches",
        "samples",
        "flips",
        "timeout",
    ] {
        if let Some(v) = p.get(key) {
            if v.as_u64().is_none_or(|n| n > i32::MAX as u64) {
                return Err(format!("invalid nonnegative integer {key}").into());
            }
        }
    }
    if p.get("seed")
        .is_some_and(|v| v.as_u64().is_none_or(|n| n > i32::MAX as u64))
    {
        return Err("invalid seed".into());
    }
    Ok(p)
}
fn boolean(p: &Json, key: &str) -> bool {
    p[key].as_bool().expect("normalized boolean")
}
fn integer(p: &Json, key: &str) -> i32 {
    p[key].as_i64().expect("normalized integer") as i32
}
fn number(p: &Json, key: &str) -> f64 {
    p[key].as_f64().expect("normalized number")
}
pub fn xyz(atom_count: usize, seed: u64) -> Vec<[f64; 3]> {
    let mut state = seed;
    (0..atom_count)
        .map(|_| {
            std::array::from_fn(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (((state >> 32) % 257) as f64 - 128.0) / 16.0
            })
        })
        .collect()
}
fn with_coordinates(m: &Molecule, seed: u64, count: usize) -> Result<Molecule> {
    let mut result = m.with_cleared_3d_conformers()?;
    for index in 0..count {
        result = result.with_added_3d_conformer(
            xyz(m.num_atoms(), seed + index as u64)
                .into_iter()
                .map(Vec::from)
                .collect(),
        )?;
    }
    Ok(result)
}
fn points<const N: usize>(rows: &[[f64; N]]) -> Value {
    Value::List(
        rows.iter()
            .map(|row| Value::List(row.iter().copied().map(Value::float).collect()))
            .collect(),
    )
}
pub(crate) fn coordinates(m: &Molecule) -> Value {
    Value::List(
        m.conformers_3d()
            .iter()
            .map(|c| {
                Value::List(vec![
                    Value::Int(c.id() as i64),
                    Value::Bool(c.is_3d()),
                    points(c.coordinates()),
                ])
            })
            .collect(),
    )
}
fn two_d_params(p: &Json) -> Coordinate2DParams {
    Coordinate2DParams {
        canonical_orientation: boolean(p, "canonical"),
        samples: integer(p, "samples") as u32,
        flips_per_sample: integer(p, "flips") as u32,
        sample_seed: integer(p, "seed"),
        permute_degree_four: boolean(p, "permute"),
        use_ring_templates: boolean(p, "ring_templates"),
        ..Default::default()
    }
}
fn two_d_value(m: &Molecule) -> Result<Value> {
    Ok(Value::List(vec![
        crate::chemistry::state(m),
        points(
            m.coordinates_2d()
                .ok_or("missing generated 2D coordinates")?,
        ),
    ]))
}
fn namespace(svg: String) -> String {
    // Approved product identity only; geometry, styles, text and numbers stay exact.
    svg.replace(
        "xmlns:rdkit='http://www.rdkit.org/xml'",
        "xmlns:tool='__tool_namespace__'",
    )
    .replace(
        "xmlns:ck='https://kit.cosmol.org/'",
        "xmlns:tool='__tool_namespace__'",
    )
    .replace("rdkit:", "tool:")
    .replace("ck:", "tool:")
}
fn alignment(m: &Molecule, p: &Json) -> Result<Value> {
    let probe = with_coordinates(m, integer(p, "seed") as u64, 3)?;
    let reference = with_coordinates(m, integer(p, "seed") as u64 + 1, 1)?;
    let weights = boolean(p, "weighted").then(|| {
        (0..m.num_atoms())
            .map(|i| 1.0 + (i % 3) as f64 / 2.0)
            .collect::<Vec<_>>()
    });
    let map = boolean(p, "explicit_map").then(|| {
        (0..m.num_atoms())
            .map(|i| AlignmentAtomMap {
                probe_atom: i,
                reference_atom: i,
            })
            .collect::<Vec<_>>()
    });
    let params = AlignmentParameters {
        atom_map: map.clone(),
        weights: weights.clone(),
        reflect: boolean(p, "reflect"),
        max_iterations: integer(p, "max_iterations") as u32,
        ..Default::default()
    };
    let best = BestAlignmentParameters {
        atom_maps: map.clone().into_iter().collect(),
        weights: weights.clone(),
        reflect: params.reflect,
        max_iterations: params.max_iterations,
        max_matches: integer(p, "max_matches"),
        symmetrize_conjugated_terminal_groups: boolean(p, "symmetrize"),
        ignore_hydrogens: false,
        ..Default::default()
    };
    let result_value = |r: AlignmentResult| {
        Value::List(vec![
            Value::float(r.rmsd),
            points(&r.transform.matrix),
            Value::List(
                r.atom_map
                    .iter()
                    .map(|a| {
                        Value::List(vec![
                            Value::Int(a.probe_atom as i64),
                            Value::Int(a.reference_atom as i64),
                        ])
                    })
                    .collect(),
            ),
        ])
    };
    let mut result_coordinates = coordinates(&probe);
    let result = match p["operation"].as_str().unwrap() {
        "transform" => result_value(probe.alignment_transform_to_with_params(&reference, &params)?),
        "best" => {
            let r = probe.best_alignment_to_with_params(&reference, &best)?;
            let rmsd = probe.best_rmsd_to_with_params(&reference, &best)?;
            Value::List(vec![result_value(r), Value::float(rmsd)])
        }
        "coordinate" => Value::float(probe.coordinate_rmsd_to_with_params(
            &reference,
            &CoordinateRmsdParameters {
                atom_maps: best.atom_maps.clone(),
                weights: weights.clone(),
                max_matches: best.max_matches,
                symmetrize_conjugated_terminal_groups: best.symmetrize_conjugated_terminal_groups,
                ..Default::default()
            },
        )?),
        "align" => {
            let (aligned, r) = probe.with_alignment_to_with_params(&reference, &params)?;
            result_coordinates = coordinates(&aligned);
            result_value(r)
        }
        "all" => {
            let rows = probe.all_conformer_best_rmsds_with_params(&AllConformerRmsdParameters {
                atom_maps: best.atom_maps.clone(),
                weights: weights.clone(),
                max_matches: best.max_matches,
                symmetrize_conjugated_terminal_groups: best.symmetrize_conjugated_terminal_groups,
                ignore_hydrogens: false,
                ..Default::default()
            })?;
            Value::List(
                rows.iter()
                    .map(|r| {
                        Value::List(vec![
                            Value::Int(r.probe_conformer_id as i64),
                            Value::Int(r.reference_conformer_id as i64),
                            Value::float(r.rmsd),
                        ])
                    })
                    .collect(),
            )
        }
        "conformers" => {
            let (aligned, r) =
                probe.with_aligned_conformers_with_params(&ConformerAlignmentParameters {
                    weights,
                    reflect: params.reflect,
                    max_iterations: params.max_iterations,
                    ..Default::default()
                })?;
            result_coordinates = coordinates(&aligned);
            Value::List(r.rmsds.into_iter().map(Value::float).collect())
        }
        _ => return Err("invalid alignment operation".into()),
    };
    Ok(Value::List(vec![
        result,
        result_coordinates,
        coordinates(&probe),
        coordinates(&reference),
    ]))
}
fn forcefield(m: &Molecule, task: &str, p: &Json) -> Result<Value> {
    let m = if boolean(p, "add_hs") {
        m.with_hydrogens()?
    } else {
        m.clone()
    };
    let mmff = task.starts_with("mmff_");
    let available = if mmff {
        m.mmff_has_all_molecule_params()?
    } else {
        m.uff_has_all_molecule_params()?
    };
    if task.ends_with("available") {
        return Ok(Value::Bool(available));
    }
    // Absence is part of the factory's source-defined parameterization result,
    // not a relabeling of a failing call for an available parameterization.
    if !available {
        return Ok(Value::List(vec![Value::Bool(false)]));
    }
    let mut m = with_coordinates(
        &m,
        integer(p, "seed") as u64,
        if task.ends_with("conformers") { 2 } else { 1 },
    )?;
    let original = coordinates(&m);
    let threshold = number(p, "threshold");
    let ignore = boolean(p, "ignore_interfragment");
    let iterations = integer(p, "iterations");
    let variant = p["variant"].as_str().unwrap();
    let make_field = |mol: &Molecule| -> Result<MolecularForceField> {
        Ok(if mmff {
            mol.mmff_force_field_with_params(&MmffForceFieldParams::new(
                None,
                variant.into(),
                threshold,
                ignore,
            ))?
        } else {
            mol.uff_force_field_with_params(&UffForceFieldParams::new(None, threshold, ignore))?
        })
    };
    let mut outcomes = Vec::new();
    if task.ends_with("force_field") {
        let mut field = make_field(&m)?;
        if boolean(p, "fixed_first") && m.num_atoms() > 0 {
            field.set_fixed_atoms_(&[AtomId::new(0)])?;
        }
        let snapshot = |f: &MolecularForceField| -> Result<Value> {
            Ok(Value::List(vec![
                Value::float(f.energy()?),
                points(&f.gradient()?),
                points(&f.positions()),
            ]))
        };
        outcomes.push(snapshot(&field)?);
        let moved = xyz(m.num_atoms(), integer(p, "seed") as u64 + 9);
        field.set_positions_(&moved)?;
        outcomes.push(snapshot(&field)?);
        let minimized = field.minimize_with_params_(&ForceFieldMinimizeParams::new(
            iterations as u32,
            number(p, "force_tolerance"),
            number(p, "energy_tolerance"),
        ))?;
        outcomes.push(Value::Bool(minimized.converged()));
        outcomes.push(snapshot(&field)?);
        outcomes.push(Value::List(
            field
                .fixed_atoms()
                .iter()
                .map(|a| Value::Int(a.index() as i64))
                .collect(),
        ));
    } else {
        // Two separate public optimizer calls, each limited to two iterations.
        for _ in 0..2 {
            let result = if task.ends_with("conformers") {
                if mmff {
                    let r = m.with_mmff_optimized_conformers_with_params(
                        &MmffConformerOptimizationParams {
                            max_iterations: iterations,
                            mmff_variant: variant.into(),
                            non_bonded_threshold: threshold,
                            ignore_interfragment_interactions: ignore,
                            num_threads: 1,
                        },
                    )?;
                    m = r.molecule;
                    Value::List(
                        r.conformer_results
                            .iter()
                            .map(|r| {
                                Value::List(vec![
                                    Value::Int(r.status_code() as i64),
                                    Value::float(r.energy()),
                                ])
                            })
                            .collect(),
                    )
                } else {
                    let r = m.with_uff_optimized_conformers_with_params(
                        &UffConformerOptimizationParams {
                            max_iterations: iterations,
                            vdw_threshold: threshold,
                            ignore_interfragment_interactions: ignore,
                            num_threads: 1,
                        },
                    )?;
                    m = r.molecule;
                    Value::List(
                        r.conformers
                            .iter()
                            .map(|r| {
                                Value::List(vec![
                                    Value::Int(r.status as i64),
                                    Value::float(r.energy),
                                ])
                            })
                            .collect(),
                    )
                }
            } else if mmff {
                let r = m.with_mmff_optimized_with_params(&MmffOptimizationParams {
                    max_iterations: iterations,
                    mmff_variant: variant.into(),
                    non_bonded_threshold: threshold,
                    ignore_interfragment_interactions: ignore,
                    ..Default::default()
                })?;
                let status = r.status_code();
                m = r.molecule;
                Value::List(vec![
                    Value::Int(status as i64),
                    Value::float(make_field(&m)?.energy()?),
                ])
            } else {
                let r = m.with_uff_optimized_with_params(&UffOptimizationParams {
                    max_iterations: iterations,
                    vdw_threshold: threshold,
                    ignore_interfragment_interactions: ignore,
                    ..Default::default()
                })?;
                m = r.molecule;
                Value::List(vec![Value::Int(r.status as i64), Value::float(r.energy)])
            };
            outcomes.push(Value::List(vec![result, coordinates(&m)]));
        }
    }
    Ok(Value::List(vec![
        Value::Bool(true),
        Value::List(outcomes),
        original,
        coordinates(&m),
    ]))
}
fn scalar(m: &Molecule, task: &str, p: &Json) -> Result<Value> {
    Ok(match task {
        "coordinates_2d" => two_d_value(&m.add_2d_coordinates_with_params(&two_d_params(p))?)?,
        "svg" => Value::Bytes(
            namespace(m.to_svg(integer(p, "width") as u32, integer(p, "height") as u32)?)
                .into_bytes(),
        ),
        "bounds_matrix" => {
            let m = if boolean(p, "add_hs") {
                m.with_hydrogens()?
            } else {
                m.clone()
            };
            Value::List(
                m.dg_bounds_matrix()?
                    .into_iter()
                    .map(|r| Value::List(r.into_iter().map(Value::float).collect()))
                    .collect(),
            )
        }
        "embedding" => {
            let m = if boolean(p, "add_hs") {
                m.with_hydrogens()?
            } else {
                m.clone()
            };
            let mut params = match p["method"].as_str().unwrap() {
                "dg" => EmbedParams::dg(),
                "etdg" => EmbedParams::etdg(),
                "etkdg" => EmbedParams::etkdg(),
                "etkdg_v2" => EmbedParams::etkdg_v2(),
                "etkdg_v3" => EmbedParams::etkdg_v3(),
                _ => return Err("invalid embedding method".into()),
            };
            params.random_seed = integer(p, "seed");
            params.num_threads = 1;
            params.max_iterations = integer(p, "max_iterations") as u32;
            params.use_random_coords = boolean(p, "random_coords");
            params.prune_rms_thresh = number(p, "prune_rms");
            params.timeout = integer(p, "timeout") as u32;
            params.track_failures = true;
            let count = integer(p, "count");
            if count == 1 {
                let r = m.with_3d_conformer_result_with_params(&params)?;
                Value::List(vec![
                    Value::List(vec![Value::Int(r.conf_id as i64)]),
                    coordinates(&r.molecule),
                    Value::List(
                        r.params
                            .failures
                            .iter()
                            .map(|n| Value::Int(*n as i64))
                            .collect(),
                    ),
                ])
            } else {
                let r = m.with_3d_conformers_result_with_params(count as u32, &params)?;
                Value::List(vec![
                    Value::List(r.conf_ids.iter().map(|n| Value::Int(*n as i64)).collect()),
                    coordinates(&r.molecule),
                    Value::List(
                        r.params
                            .failures
                            .iter()
                            .map(|n| Value::Int(*n as i64))
                            .collect(),
                    ),
                ])
            }
        }
        "alignment" => alignment(m, p)?,
        _ => forcefield(m, task, p)?,
    })
}
pub fn evaluate(
    task: &str,
    profiles: &[Json],
    rows: &[String],
    first: usize,
    threads: usize,
    pool: &rayon::ThreadPool,
) -> Result<Vec<Vec<Vec<u8>>>> {
    if matches!(task, "svg" | "bounds_matrix") {
        return query_batch(task, profiles, rows, first, threads, pool);
    }
    if task != "coordinates_2d" {
        return Ok(pool.install(|| {
            rows.par_iter()
                .map(|text| match Molecule::from_smiles(text) {
                    Err(e) => vec![tasks::parse_failure(&e); profiles.len()],
                    Ok(m) => profiles
                        .iter()
                        .map(|p| {
                            scalar(&m, task, p)
                                .map(tasks::success)
                                .unwrap_or_else(|e| tasks::failed_error(e.as_ref()))
                        })
                        .collect(),
                })
                .collect()
        }));
    }
    let split = 10usize.saturating_sub(first).min(rows.len());
    let mut output = vec![Vec::new(); rows.len()];
    let params = BatchParams {
        n_jobs: Some(threads),
        errors: Some(BatchErrorMode::KeepErrors),
        progress_bar: Some(false),
    };
    let batch = MoleculeBatch::from_smiles_list_with_params(
        &rows[split..],
        &SmilesParseParams::default(),
        &params,
    )?;
    for p in profiles {
        for (slot, text) in output[..split].iter_mut().zip(&rows[..split]) {
            slot.push(match Molecule::from_smiles(text) {
                Err(e) => tasks::parse_failure(&e),
                Ok(m) => scalar(&m, task, p)
                    .map(tasks::success)
                    .unwrap_or_else(|e| tasks::failed_error(e.as_ref())),
            });
        }
        let transformed = batch.add_2d_coordinates_with_params(&two_d_params(p), &params)?;
        for (slot, (before, after)) in output[split..]
            .iter_mut()
            .zip(batch.records().iter().zip(transformed.records()))
        {
            slot.push(match (before, after) {
                (BatchRecord::Error(_), _) => {
                    tasks::outcome(1, Value::Text("parse_rejected".into()))
                }
                (BatchRecord::Molecule(_), BatchRecord::Molecule(m)) => {
                    tasks::success(two_d_value(m)?)
                }
                (BatchRecord::Molecule(_), BatchRecord::Error(e)) => {
                    tasks::outcome(5, Value::Text(e.to_string()))
                }
            });
        }
    }
    Ok(output)
}

fn query_batch(
    task: &str,
    profiles: &[Json],
    rows: &[String],
    first: usize,
    threads: usize,
    pool: &rayon::ThreadPool,
) -> Result<Vec<Vec<Vec<u8>>>> {
    let molecules = pool.install(|| {
        rows.par_iter()
            .map(|text| Molecule::from_smiles(text))
            .collect::<Vec<_>>()
    });
    let split = 10usize.saturating_sub(first).min(rows.len());
    let mut output = molecules
        .iter()
        .map(|mol| {
            vec![
                mol.as_ref()
                    .err()
                    .map(tasks::parse_failure)
                    .unwrap_or_default();
                profiles.len()
            ]
        })
        .collect::<Vec<_>>();
    let query = BatchQueryParams {
        n_jobs: Some(threads),
        progress_bar: Some(false),
        ..Default::default()
    };
    let transform = BatchParams {
        errors: Some(BatchErrorMode::KeepErrors),
        n_jobs: Some(threads),
        progress_bar: Some(false),
    };
    for (pi, p) in profiles.iter().enumerate() {
        for (i, mol) in molecules[..split].iter().enumerate() {
            if let Ok(mol) = mol {
                output[i][pi] = scalar(mol, task, p)
                    .map(tasks::success)
                    .unwrap_or_else(|e| tasks::failed_error(e.as_ref()));
            }
        }
        let mut indices = (split..rows.len())
            .filter(|i| molecules[*i].is_ok())
            .collect::<Vec<_>>();
        while !indices.is_empty() {
            let mut batch = MoleculeBatch::from_records(
                indices
                    .iter()
                    .map(|i| BatchRecord::Molecule(molecules[*i].as_ref().unwrap().clone()))
                    .collect(),
                BatchErrorMode::KeepErrors,
            )?;
            if task == "bounds_matrix" && boolean(p, "add_hs") {
                batch = batch.with_hydrogens_with_params(&AddHsParams::default(), &transform)?;
                for (i, record) in batch.records().iter().enumerate() {
                    if let BatchRecord::Error(e) = record {
                        output[indices[i]][pi] = tasks::failed_error(e);
                    }
                }
            }
            let values = if task == "svg" {
                batch
                    .to_svg_list_with_params(
                        integer(p, "width") as u32,
                        integer(p, "height") as u32,
                        &query,
                    )
                    .map(|rows| {
                        rows.into_iter()
                            .map(|v| v.map(|s| Value::Bytes(namespace(s).into_bytes())))
                            .collect::<Vec<_>>()
                    })
            } else {
                batch.dg_bounds_matrix_list_with_params(&query).map(|rows| {
                    rows.into_iter()
                        .map(|v| {
                            v.map(|matrix| {
                                Value::List(
                                    matrix
                                        .into_iter()
                                        .map(|r| {
                                            Value::List(r.into_iter().map(Value::float).collect())
                                        })
                                        .collect(),
                                )
                            })
                        })
                        .collect()
                })
            };
            match values {
                Ok(values) => {
                    if values.len() != indices.len() {
                        return Err("batch query changed dimensions".into());
                    }
                    for (i, v) in indices.iter().zip(values) {
                        if let Some(v) = v {
                            output[*i][pi] = tasks::success(v);
                        } else if output[*i][pi].is_empty() {
                            return Err("batch omitted a valid result".into());
                        }
                    }
                    break;
                }
                Err(error) => {
                    let mut failed = Vec::new();
                    for e in error.record_errors {
                        let i = *indices
                            .get(e.index)
                            .ok_or("batch error index out of bounds")?;
                        output[i][pi] = tasks::failed_error(&e);
                        failed.push(i);
                    }
                    if failed.is_empty() {
                        return Err("batch failed without indexed errors".into());
                    }
                    indices.retain(|i| !failed.contains(i));
                }
            }
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mmff_optimizer_cache_requires_nonmutating_observations() {
        for task in ["mmff_optimize", "mmff_optimize_conformers"] {
            let p = profiles(task).unwrap().remove(0);
            let mut old = p.clone();
            old.as_object_mut().unwrap().remove("optimizer_observation");
            assert!(normalize(task, old).is_err());
            let mut wrong = p.clone();
            wrong["optimizer_observation"] = json!("prepare_live_input");
            assert!(normalize(task, wrong).is_err());
            assert_eq!(normalize(task, p.clone()).unwrap(), p);
        }
    }
    #[test]
    fn geometry_profiles_and_coordinate_recipe_are_explicit() {
        assert_eq!(xyz(2, 17), xyz(2, 17));
        assert_ne!(xyz(2, 17), xyz(2, 18));
        for (task, _) in TASKS {
            for p in profiles(task).unwrap() {
                assert_eq!(normalize(task, p.clone()).unwrap(), p);
            }
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        for task in [
            "coordinates_2d",
            "svg",
            "bounds_matrix",
            "uff_force_field",
            "mmff_force_field",
        ] {
            let p = profiles(task).unwrap();
            let results = evaluate(task, &p, &vec!["CCO".into(); 12], 0, 2, &pool).unwrap();
            assert_eq!(results.len(), 12);
            for row in &results {
                for result in row {
                    assert_eq!(
                        result[0],
                        0,
                        "{task}: {:?}",
                        crate::value::decode(&result[1..])
                    );
                }
            }
        }
    }
}
