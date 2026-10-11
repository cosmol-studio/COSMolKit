//! Fixed numerical regressions; every float is encoded with its original bits.
use crate::{Result, tasks, value::Value};
use cosmolkit::*;
use serde_json::{Value as Json, json};

fn embed_params(case: &Json) -> Result<EmbedParams> {
    let mut p = match case["preset_name"]
        .as_str()
        .ok_or("missing embedding preset")?
    {
        "DG" => EmbedParams::dg(),
        "KDG" => EmbedParams::kdg(),
        "ETDG" => EmbedParams::etdg(),
        "ETDGv2" => EmbedParams::etdg_v2(),
        "ETKDG" => EmbedParams::etkdg(),
        "ETKDGv2" => EmbedParams::etkdg_v2(),
        "ETKDGv3" => EmbedParams::etkdg_v3(),
        "srETKDGv3" => EmbedParams::sr_etkdg_v3(),
        _ => return Err("unknown embedding preset".into()),
    };
    for (name, value) in case["attrs"]
        .as_object()
        .ok_or("missing embedding attributes")?
    {
        match name.as_str() {
            "randomSeed" => p.random_seed = serde_json::from_value(value.clone())?,
            "numThreads" => {
                if value != 1 {
                    return Err(
                        "per-case embedding is single threaded; outer pool owns concurrency".into(),
                    );
                }
                p.num_threads = 1;
            }
            "maxIterations" => p.max_iterations = serde_json::from_value(value.clone())?,
            "timeout" => p.timeout = serde_json::from_value(value.clone())?,
            "trackFailures" => p.track_failures = serde_json::from_value(value.clone())?,
            "useRandomCoords" => p.use_random_coords = serde_json::from_value(value.clone())?,
            "embedFragmentsSeparately" => {
                p.embed_fragments_separately = serde_json::from_value(value.clone())?
            }
            "pruneRmsThresh" => p.prune_rms_thresh = serde_json::from_value(value.clone())?,
            "useSymmetryForPruning" => {
                p.use_symmetry_for_pruning = serde_json::from_value(value.clone())?
            }
            "enableSequentialRandomSeeds" => {
                p.enable_sequential_random_seeds = serde_json::from_value(value.clone())?
            }
            "forceTransAmides" => p.force_trans_amides = serde_json::from_value(value.clone())?,
            "useExpTorsionAnglePrefs" => {
                p.use_exp_torsion_angle_prefs = serde_json::from_value(value.clone())?
            }
            "useBasicKnowledge" => p.use_basic_knowledge = serde_json::from_value(value.clone())?,
            _ => return Err(format!("unknown embedding attribute {name}").into()),
        }
    }
    if let Some(v) = case.get("coord_map") {
        p.coord_map = Some(serde_json::from_value(v.clone())?);
    }
    if let Some(v) = case.get("cpci") {
        p.cpci = Some(
            v.as_object()
                .ok_or("invalid CPCI map")?
                .iter()
                .map(|(k, v)| {
                    let (a, b) = k.split_once('-').ok_or("invalid CPCI atom pair")?;
                    Ok((
                        (a.parse()?, b.parse()?),
                        v.as_f64().ok_or("invalid CPCI value")?,
                    ))
                })
                .collect::<Result<_>>()?,
        );
    }
    Ok(p)
}
fn embedding(task: &str, row: &Json) -> Result<Vec<u8>> {
    let c = if task == "conformer_library" {
        json!({"mode":"single","source_kind":"smiles_with_hydrogens","source":row["case"],"preset_name":"ETKDGv3","attrs":{"maxIterations":3,"randomSeed":61453,"numThreads":1,"timeout":0}})
    } else {
        row["case"].clone()
    };
    let m = if c["source_kind"] == "fixture_mol" {
        Molecule::from_mol_with_params(
            c["mol_block"].as_str().ok_or("missing MOL input")?,
            &SdfReadParams {
                sanitize: true,
                remove_hs: false,
                ..Default::default()
            },
        )?
    } else {
        let m = match Molecule::from_smiles(c["source"].as_str().ok_or("missing embedding SMILES")?)
        {
            Ok(m) => m,
            Err(e) => return Ok(tasks::parse_failure(&e)),
        };
        if c["source_kind"] == "smiles_with_hydrogens" {
            m.with_hydrogens()?
        } else {
            m
        }
    };
    let p = embed_params(&c)?;
    let before = (
        crate::chemistry::state(&m),
        crate::geometry::coordinates(&m),
    );
    let (output, status, failures) = if c["mode"] == "single" {
        let r = m.with_3d_conformer_result_with_params(&p)?;
        (r.molecule, Value::Int(r.conf_id as i64), r.params.failures)
    } else {
        let r = m.with_3d_conformers_result_with_params(
            c["num_confs"].as_u64().ok_or("missing conformer count")? as u32,
            &p,
        )?;
        (
            r.molecule,
            Value::List(
                r.conf_ids
                    .into_iter()
                    .map(|v| Value::Int(v as i64))
                    .collect(),
            ),
            r.params.failures,
        )
    };
    if before
        != (
            crate::chemistry::state(&m),
            crate::geometry::coordinates(&m),
        )
    {
        return Err("embedding modified the source".into());
    }
    Ok(tasks::success(Value::List(vec![
        status,
        crate::geometry::coordinates(&output),
        if p.track_failures {
            Value::List(failures.into_iter().map(|v| Value::Int(v as i64)).collect())
        } else {
            Value::Null
        },
    ])))
}
pub fn observe(task: &str, row: &Json) -> Result<Vec<u8>> {
    match task {
        "conformer_fixed19" | "conformer_library" => embedding(task, row),
        "molalign_focused" => Ok(tasks::success(focused_alignment(&row["case"])?)),
        "forcefield_optimizers" => Ok(tasks::success(optimizer(row)?)),
        _ => Err("unknown numerical special regression".into()),
    }
}

fn optimizer(row: &Json) -> Result<Value> {
    let m = Molecule::from_smiles(
        row["single_cxsmiles"]
            .as_str()
            .ok_or("missing shared optimizer geometry")?,
    )?;
    let mmff = row["case"]["forcefield"] == "mmff";
    let field = if mmff {
        m.mmff_force_field_with_params(&MmffForceFieldParams::new(
            Some(0),
            "MMFF94".into(),
            100.0,
            true,
        ))?
    } else {
        m.uff_force_field_with_params(&UffForceFieldParams::new(Some(0), 100.0, true))?
    };
    let initial = Value::List(vec![
        Value::float(field.energy()?),
        Value::List(
            field
                .gradient()?
                .iter()
                .flatten()
                .copied()
                .map(Value::float)
                .collect(),
        ),
        crate::geometry::coordinates(&m),
    ]);
    let single = if mmff {
        let r = m.with_mmff_optimized_with_params(&MmffOptimizationParams {
            max_iterations: 2,
            mmff_variant: "MMFF94".into(),
            non_bonded_threshold: 100.0,
            ignore_interfragment_interactions: true,
            ..Default::default()
        })?;
        Value::List(vec![
            Value::Int(r.status_code() as i64),
            Value::float(r.molecule.mmff_force_field()?.energy()?),
            crate::geometry::coordinates(&r.molecule),
        ])
    } else {
        let r = m.with_uff_optimized_with_params(&UffOptimizationParams {
            max_iterations: 2,
            vdw_threshold: 100.0,
            ignore_interfragment_interactions: true,
            ..Default::default()
        })?;
        Value::List(vec![
            Value::Int(r.status as i64),
            Value::float(r.energy),
            crate::geometry::coordinates(&r.molecule),
        ])
    };
    let mut multi = m.with_cleared_3d_conformers()?;
    for c in row["multi_conformers"]
        .as_array()
        .ok_or("missing shared multiconformer geometry")?
    {
        multi = multi.with_added_3d_conformer(serde_json::from_value(c[2].clone())?)?;
    }
    let before = crate::geometry::coordinates(&multi);
    let multiple = if mmff {
        let r =
            multi.with_mmff_optimized_conformers_with_params(&MmffConformerOptimizationParams {
                max_iterations: 2,
                mmff_variant: "MMFF94".into(),
                non_bonded_threshold: 100.0,
                ignore_interfragment_interactions: true,
                num_threads: 1,
            })?;
        Value::List(vec![
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
            ),
            crate::geometry::coordinates(&r.molecule),
        ])
    } else {
        let r =
            multi.with_uff_optimized_conformers_with_params(&UffConformerOptimizationParams {
                max_iterations: 2,
                vdw_threshold: 100.0,
                ignore_interfragment_interactions: true,
                num_threads: 1,
            })?;
        Value::List(vec![
            Value::List(
                r.conformers
                    .iter()
                    .map(|r| Value::List(vec![Value::Int(r.status as i64), Value::float(r.energy)]))
                    .collect(),
            ),
            crate::geometry::coordinates(&r.molecule),
        ])
    };
    if before != crate::geometry::coordinates(&multi) {
        return Err("forcefield optimizer modified the source".into());
    }
    let available = if mmff {
        m.mmff_has_all_molecule_params()?
    } else {
        m.uff_has_all_molecule_params()?
    };
    Ok(Value::List(vec![
        Value::Bool(available),
        initial,
        single,
        before,
        multiple,
        crate::geometry::coordinates(&m),
    ]))
}

fn alignment_molecule(row: &Json) -> Result<Molecule> {
    let m = Molecule::from_smiles(row["smiles"].as_str().ok_or("missing alignment SMILES")?)?;
    let conformers = row["conformers"]
        .as_array()
        .ok_or("missing alignment conformers")?
        .iter()
        .map(|c| {
            Ok(Conformer3D::new(
                c["id"].as_u64().ok_or("missing conformer id")? as usize,
                serde_json::from_value(c["coordinates"].clone())?,
                true,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(Molecule::from_parts(
        m.topology().clone(),
        CoordinateBlock {
            conformers_3d: conformers,
            ..Default::default()
        },
        m.properties().clone(),
    )?)
}
fn alignment_value(r: AlignmentResult) -> Value {
    Value::List(vec![
        Value::float(r.rmsd),
        Value::List(
            r.transform
                .matrix
                .into_iter()
                .map(|row| Value::List(row.into_iter().map(Value::float).collect()))
                .collect(),
        ),
        Value::List(
            r.atom_map
                .into_iter()
                .map(|a| {
                    Value::List(vec![
                        Value::Int(a.probe_atom as i64),
                        Value::Int(a.reference_atom as i64),
                    ])
                })
                .collect(),
        ),
    ])
}
fn alignment_error(e: AlignmentError) -> Result<Value> {
    Ok(match e {
        AlignmentError::ConformerNotFound { id } => Value::List(vec![
            Value::Text("conformer_not_found".into()),
            Value::Int(id as i64),
        ]),
        AlignmentError::WeightCountMismatch {
            map_len,
            weight_len,
        } => Value::List(vec![
            Value::Text("weight_count_mismatch".into()),
            Value::Int(map_len as i64),
            Value::Int(weight_len as i64),
        ]),
        AlignmentError::NoSubstructureMatch => {
            Value::List(vec![Value::Text("no_substructure_match".into())])
        }
        e => return Err(e.into()),
    })
}
fn focused_alignment(c: &Json) -> Result<Value> {
    let probe = alignment_molecule(c.get("probe").unwrap_or(&c["molecule"]))?;
    let reference = c.get("reference").map(alignment_molecule).transpose()?;
    let before = crate::geometry::coordinates(&probe);
    let before_ref = reference
        .as_ref()
        .map(crate::geometry::coordinates)
        .unwrap_or(Value::Null);
    let mut out = Vec::new();
    for p in c["calls"].as_array().ok_or("missing alignment calls")? {
        let mut after = before.clone();
        let maps: Vec<Vec<[usize; 2]>> =
            serde_json::from_value(p.get("atom_maps").cloned().unwrap_or(json!([])))?;
        let maps = maps
            .into_iter()
            .map(|m| {
                m.into_iter()
                    .map(|a| AlignmentAtomMap {
                        probe_atom: a[0],
                        reference_atom: a[1],
                    })
                    .collect()
            })
            .collect::<Vec<_>>();
        let map = p
            .get("atom_map")
            .map(|v| serde_json::from_value::<Vec<[usize; 2]>>(v.clone()))
            .transpose()?
            .map(|m| {
                m.into_iter()
                    .map(|a| AlignmentAtomMap {
                        probe_atom: a[0],
                        reference_atom: a[1],
                    })
                    .collect()
            });
        let weights = p
            .get("weights")
            .map(|v| serde_json::from_value(v.clone()))
            .transpose()?;
        let n = |key: &str, default: i64| p[key].as_i64().unwrap_or(default);
        let b = |key: &str, default: bool| p[key].as_bool().unwrap_or(default);
        let basic = AlignmentParameters {
            probe_conformer_id: n("probe_conformer_id", -1) as i32,
            reference_conformer_id: n("reference_conformer_id", -1) as i32,
            atom_map: map,
            weights: weights.clone(),
            reflect: b("reflect", false),
            max_iterations: n("max_iterations", 50) as u32,
        };
        let best = BestAlignmentParameters {
            probe_conformer_id: basic.probe_conformer_id,
            reference_conformer_id: basic.reference_conformer_id,
            atom_maps: maps.clone(),
            weights: weights.clone(),
            reflect: basic.reflect,
            max_iterations: basic.max_iterations,
            max_matches: n("max_matches", 1_000_000) as i32,
            symmetrize_conjugated_terminal_groups: b("symmetrize_conjugated_terminal_groups", true),
            ignore_hydrogens: b("ignore_hydrogens", true),
            num_threads: n("num_threads", 1) as i32,
        };
        let result: std::result::Result<Value, AlignmentError> = match p["operation"]
            .as_str()
            .ok_or("missing alignment operation")?
        {
            "alignment_transform" => probe
                .alignment_transform_to_with_params(
                    reference.as_ref().ok_or("missing reference")?,
                    &basic,
                )
                .map(alignment_value),
            "best_alignment" => probe
                .best_alignment_to_with_params(
                    reference.as_ref().ok_or("missing reference")?,
                    &best,
                )
                .map(alignment_value),
            "coordinate_rmsd" => probe
                .coordinate_rmsd_to_with_params(
                    reference.as_ref().ok_or("missing reference")?,
                    &CoordinateRmsdParameters {
                        probe_conformer_id: best.probe_conformer_id,
                        reference_conformer_id: best.reference_conformer_id,
                        atom_maps: maps,
                        weights: weights.clone(),
                        max_matches: best.max_matches,
                        symmetrize_conjugated_terminal_groups: best
                            .symmetrize_conjugated_terminal_groups,
                    },
                )
                .map(Value::float),
            "all_conformer_best_rms" => probe
                .all_conformer_best_rmsds_with_params(&AllConformerRmsdParameters {
                    atom_maps: maps,
                    weights: weights.clone(),
                    max_matches: best.max_matches,
                    symmetrize_conjugated_terminal_groups: best
                        .symmetrize_conjugated_terminal_groups,
                    ignore_hydrogens: best.ignore_hydrogens,
                    num_threads: best.num_threads,
                })
                .map(|rows| {
                    Value::List(
                        rows.into_iter()
                            .map(|r| {
                                Value::List(vec![
                                    Value::Int(r.probe_conformer_id as i64),
                                    Value::Int(r.reference_conformer_id as i64),
                                    Value::float(r.rmsd),
                                ])
                            })
                            .collect(),
                    )
                }),
            "align_conformers" => {
                let params = ConformerAlignmentParameters {
                    atom_indices: p
                        .get("atom_indices")
                        .map(|v| serde_json::from_value(v.clone()))
                        .transpose()?,
                    conformer_ids: p
                        .get("conformer_ids")
                        .map(|v| serde_json::from_value(v.clone()))
                        .transpose()?,
                    weights: weights.clone(),
                    reflect: basic.reflect,
                    max_iterations: basic.max_iterations,
                };
                match probe.with_aligned_conformers_with_params(&params) {
                    Ok((m, r)) => {
                        after = crate::geometry::coordinates(&m);
                        Ok(Value::List(r.rmsds.into_iter().map(Value::float).collect()))
                    }
                    Err(OperationError::Alignment(e)) => Err(e),
                    Err(e) => return Err(e.into()),
                }
            }
            _ => return Err("unknown alignment operation".into()),
        };
        let result = match result {
            Ok(v) => Value::List(vec![Value::Bool(true), v]),
            Err(e) => Value::List(vec![Value::Bool(false), alignment_error(e)?]),
        };
        out.push(Value::List(vec![
            result,
            after,
            before.clone(),
            before_ref.clone(),
        ]));
    }
    if before != crate::geometry::coordinates(&probe)
        || before_ref
            != reference
                .as_ref()
                .map(crate::geometry::coordinates)
                .unwrap_or(Value::Null)
    {
        return Err("alignment modified a source".into());
    }
    Ok(Value::List(out))
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn every_fixed_embedding_attribute_is_projected() {
        for row in include_str!("../testdata/special/conformer_fixed19.jsonl").lines() {
            let row: Json = serde_json::from_str(row).unwrap();
            let p = embed_params(&row["case"]).unwrap();
            assert_eq!(p.num_threads, 1);
        }
        assert!(embed_params(&json!({"preset_name":"DG","attrs":{"bogus":1}})).is_err());
    }
}
