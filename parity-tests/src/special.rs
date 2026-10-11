//! Fixed regressions use explicit self-contained rows, not historical runners.
use crate::{Result, tasks, value::Value};
use cosmolkit::*;
use rayon::prelude::*;
use serde_json::{Value as Json, json};

pub const TASKS: &[(&str, &str)] = &[
    ("structure_tags", "special"),
    ("forcefield_optimizers", "special"),
    ("molalign_focused", "special"),
    ("conformer_fixed19", "special"),
    ("conformer_library", "special"),
    ("tautomer_focused", "special"),
    ("tautomer_long_conjugated", "special"),
    ("mmff_builtin", "special"),
    ("forcefield_properties", "special"),
    ("mcs_upstream", "special"),
    ("mcs_jnk1", "special"),
    ("bio_mmcif_switches", "special"),
];
pub fn input(task: &str) -> Option<&'static str> {
    match task {
        "structure_tags" => Some("testdata/special/structure_tags.jsonl"),
        "forcefield_optimizers" => Some("testdata/special/forcefield_optimizers.jsonl"),
        "molalign_focused" => Some("testdata/special/molalign_focused.jsonl"),
        "conformer_fixed19" => Some("testdata/special/conformer_fixed19.jsonl"),
        "conformer_library" => Some("testdata/special/conformer_library.jsonl"),
        "tautomer_focused" => Some("testdata/special/tautomer_focused.jsonl"),
        "tautomer_long_conjugated" => Some("testdata/special/tautomer_long_conjugated.jsonl"),
        "mmff_builtin" => Some("testdata/special/mmff_builtin.jsonl"),
        "forcefield_properties" => Some("testdata/special/forcefield_properties.jsonl"),
        "mcs_upstream" => Some("testdata/special/mcs_upstream.jsonl"),
        "mcs_jnk1" => Some("testdata/special/mcs_jnk1.jsonl"),
        "bio_mmcif_switches" => Some("testdata/special/bio_mmcif_switches.jsonl"),
        _ => None,
    }
}
pub fn profiles(task: &str) -> Vec<Json> {
    vec![if task == "structure_tags" {
        json!({"result":"typed_structure_state_v1","error_result":"no_return_value","source_preservation":"cow_atomic"})
    } else if task == "bio_mmcif_switches" {
        json!({"result":"complete_mmcif_switch_matrix_v1","gemmi_version":crate::bio::GEMMI})
    } else {
        json!({"result":"complete_fixed_case_v1"})
    }]
}
pub fn normalize(task: &str, p: Json) -> Result<Json> {
    if !profiles(task).contains(&p) {
        return Err("fixed regression parameters belong in its input rows".into());
    }
    Ok(p)
}
pub(crate) fn molecule(row: &Json) -> Result<Molecule> {
    let text = row["text"].as_str().ok_or("missing molecule text")?;
    let sanitize = row["sanitize"].as_bool().ok_or("missing sanitize")?;
    let remove_hs = row["remove_hs"].as_bool().ok_or("missing remove_hs")?;
    match row["format"].as_str() {
        Some("smiles") => Ok(Molecule::from_smiles_with_params(
            text,
            &SmilesParseParams {
                sanitize,
                remove_hs,
                ..Default::default()
            },
        )?),
        Some("mol") => Ok(Molecule::from_mol_with_params(
            text,
            &SdfReadParams {
                sanitize,
                remove_hs,
                ..Default::default()
            },
        )?),
        _ => Err("unknown fixed input format".into()),
    }
}
fn mmff(m: &Molecule, variant: &str) -> Result<Value> {
    let available = m.mmff_has_all_molecule_params()?;
    let p = m.mmff_properties_with_params(&MmffPropertiesParams {
        mmff_variant: variant.into(),
    })?;
    let atoms = if p.is_valid() {
        Value::List(
            (0..m.num_atoms())
                .map(|i| {
                    Ok(Value::List(vec![
                        Value::Int(p.atom_type(i)? as i64),
                        Value::float(p.formal_charge(i)?),
                        Value::float(p.partial_charge(i)?),
                    ]))
                })
                .collect::<Result<_>>()?,
        )
    } else {
        Value::Null
    };
    Ok(Value::List(vec![
        Value::Int(m.num_atoms() as i64),
        Value::Bool(available),
        Value::Bool(p.is_valid()),
        atoms,
    ]))
}

fn mcs_params(values: &Json) -> Result<McsParameters> {
    let mut p = McsParameters::default();
    macro_rules! fields {($target:expr, $value:expr, $($name:ident),*) => {
        for (key,value) in $value.as_object().ok_or("MCS parameters must be an object")? {
            match key.as_str() {$(stringify!($name)=>$target.$name=serde_json::from_value(value.clone())?,)* _=>return Err(format!("unknown MCS field {key}").into())}
        }
    };}
    for (key, v) in values.as_object().ok_or("missing MCS parameters")? {
        match key.as_str() {
            "atom_comparator" => {
                p.atom_comparator = match v.as_str() {
                    Some("any") => McsAtomComparator::AtomCompareAny,
                    Some("elements") => McsAtomComparator::AtomCompareElements,
                    Some("isotopes") => McsAtomComparator::AtomCompareIsotopes,
                    Some("any_heavy_atom") => McsAtomComparator::AtomCompareAnyHeavyAtom,
                    _ => return Err("invalid MCS atom comparator".into()),
                }
            }
            "bond_comparator" => {
                p.bond_comparator = match v.as_str() {
                    Some("any") => McsBondComparator::BondCompareAny,
                    Some("order") => McsBondComparator::BondCompareOrder,
                    Some("order_exact") => McsBondComparator::BondCompareOrderExact,
                    _ => return Err("invalid MCS bond comparator".into()),
                }
            }
            "atom_compare_parameters" => fields!(
                p.atom_compare_parameters,
                v,
                match_valences,
                match_chiral_tag,
                match_formal_charge,
                ring_matches_ring_only,
                complete_rings_only,
                match_isotope,
                max_distance
            ),
            "bond_compare_parameters" => fields!(
                p.bond_compare_parameters,
                v,
                ring_matches_ring_only,
                complete_rings_only,
                match_fused_rings,
                match_fused_rings_strict,
                match_stereo
            ),
            _ => fields!(
                p,
                json!({key:v}),
                store_all,
                maximize_bonds,
                threshold,
                timeout,
                verbose,
                initial_seed
            ),
        }
    }
    Ok(p)
}
fn mcs(row: &Json) -> Result<Value> {
    let mols = row["molecules"]
        .as_array()
        .ok_or("missing MCS molecules")?
        .iter()
        .map(molecule)
        .collect::<Result<Vec<_>>>()?;
    let before = mols
        .iter()
        .map(Molecule::to_binary)
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let result = maximum_common_substructure_with_params(
        &mols.iter().collect::<Vec<_>>(),
        &mcs_params(&row["case"]["parameters"])?,
    )?;
    if !result.completed {
        return Ok(Value::List(vec![Value::Bool(false)]));
    }
    let query = if let Some(q) = &result.query {
        Value::List(vec![
            Value::Int(q.num_atoms() as i64),
            Value::Int(q.num_bonds() as i64),
            Value::Bytes(write_smarts(q, &SmartsWriteParams::default())?.into_bytes()),
            Value::List(
                mols.iter()
                    .map(|m| m.has_substruct_match(q).map(Value::Bool))
                    .collect::<std::result::Result<_, _>>()?,
            ),
        ])
    } else {
        Value::Null
    };
    let mut alternatives = result
        .degenerate
        .keys()
        .map(|s| s.as_bytes().to_vec())
        .collect::<Vec<_>>();
    alternatives.sort(); // Unordered map keys only, not SMARTS normalization.
    if before
        != mols
            .iter()
            .map(Molecule::to_binary)
            .collect::<std::result::Result<Vec<_>, _>>()?
    {
        return Err("MCS modified an input".into());
    }
    Ok(Value::List(vec![
        Value::Bool(true),
        Value::Int(result.atom_count as i64),
        Value::Int(result.bond_count as i64),
        Value::Bytes(result.smarts.into_bytes()),
        Value::List(alternatives.into_iter().map(Value::Bytes).collect()),
        query,
    ]))
}

macro_rules! cif_switches {($($field:ident),*)=>{
    pub const CIF_FLAGS:&[&str]=&[$(stringify!($field),)*"auth_all"];
    fn cif_params(all:bool,flag:&str)->BioMmcifWriteParams { BioMmcifWriteParams {
        $($field:if flag==stringify!($field){!all}else{all},)*auth_all:flag=="auth_all",modres:all,..Default::default()
    }}
};}
cif_switches!(
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
    group_pdb
);
fn cif(row: &Json) -> Result<Value> {
    let folder = tempfile::tempdir()?;
    let filename = std::path::Path::new(
        row["case"]["input"]
            .as_str()
            .ok_or("missing BIO filename")?,
    )
    .file_name()
    .ok_or("missing filename")?;
    let path = folder.path().join(filename);
    std::fs::write(&path, row["text"].as_str().ok_or("missing BIO text")?)?;
    let s = BioStructure::read(&path).map_err(|e| format!("BIO fixed input read failed: {e}"))?;
    Ok(Value::List(
        [false, true]
            .into_iter()
            .flat_map(|all| CIF_FLAGS.iter().map(move |flag| (all, *flag)))
            .map(|(all, flag)| {
                s.to_mmcif_with_params(&cif_params(all, flag))
                    .map(|s| Value::Bytes(s.into_bytes()))
            })
            .collect::<std::result::Result<_, _>>()?,
    ))
}
fn observe(task: &str, row: &Json) -> Result<Vec<u8>> {
    if task == "structure_tags" {
        return crate::structure_tags::observe(row);
    }
    if matches!(
        task,
        "conformer_fixed19" | "conformer_library" | "molalign_focused" | "forcefield_optimizers"
    ) {
        return crate::special_geometry::observe(task, row);
    }
    let c = &row["case"];
    if task.starts_with("tautomer_") {
        let mol = Molecule::from_smiles_with_params(
            c["smiles"].as_str().ok_or("missing tautomer SMILES")?,
            &SmilesParseParams {
                sanitize: c["sanitize"].as_bool().ok_or("missing sanitize")?,
                remove_hs: c["remove_hs"].as_bool().ok_or("missing remove_hs")?,
                ..Default::default()
            },
        );
        let mol = match mol {
            Ok(m) => m,
            Err(e) => return Ok(tasks::parse_failure(&e)),
        };
        let values = row["branches"]
            .as_array()
            .ok_or("missing branches")?
            .iter()
            .map(|branch| -> Result<Value> {
                let mut branch = branch.clone();
                branch.as_object_mut().unwrap().remove("name");
                let p: crate::search_tautomer::TautomerProfile = serde_json::from_value(branch)?;
                crate::search_tautomer::tautomer(&mol, &p.params()?)
            })
            .collect::<Result<_>>()?;
        return Ok(tasks::success(Value::List(values)));
    }
    let value = match task {
        "mmff_builtin" => {
            let m = match Molecule::from_smiles(c["smiles"].as_str().ok_or("missing MMFF SMILES")?)
            {
                Ok(m) => m,
                Err(e) => return Ok(tasks::parse_failure(&e)),
            };
            mmff(&m, c["variant"].as_str().ok_or("missing MMFF variant")?)?
        }
        "forcefield_properties" => {
            let m = match Molecule::from_smiles(c.as_str().ok_or("missing properties SMILES")?) {
                Ok(m) => m,
                Err(e) => return Ok(tasks::parse_failure(&e)),
            };
            let h = m.with_hydrogens()?;
            Value::List(
                [&m, &h]
                    .into_iter()
                    .map(|m| {
                        Ok(Value::List(vec![
                            Value::Bool(m.uff_has_all_molecule_params()?),
                            mmff(m, "MMFF94")?,
                        ]))
                    })
                    .collect::<Result<_>>()?,
            )
        }
        "mcs_upstream" | "mcs_jnk1" => mcs(row)?,
        "bio_mmcif_switches" => cif(row)?,
        _ => return Err("unknown special regression".into()),
    };
    Ok(tasks::success(value))
}
pub fn evaluate(
    task: &str,
    rows: &[String],
    pool: &rayon::ThreadPool,
) -> Result<Vec<Vec<Vec<u8>>>> {
    pool.install(|| {
        rows.par_iter()
            .map(|text| {
                let row: Json = serde_json::from_str(text)?;
                Ok(vec![
                    observe(task, &row).unwrap_or_else(|e| tasks::failed_error(e.as_ref())),
                ])
            })
            .collect()
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn special_inputs_are_self_contained_and_have_complete_censuses() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        for (task, count) in [
            ("structure_tags", 77),
            ("forcefield_optimizers", 24),
            ("molalign_focused", 8),
            ("conformer_fixed19", 19),
            ("conformer_library", 152),
            ("tautomer_focused", 18),
            ("tautomer_long_conjugated", 1),
            ("mmff_builtin", 2052),
            ("forcefield_properties", 152),
            ("mcs_upstream", 44),
            ("mcs_jnk1", 210),
            ("bio_mmcif_switches", 2),
        ] {
            let input = std::fs::read_to_string(root.join(input(task).unwrap())).unwrap();
            assert_eq!(input.lines().count(), count);
            for line in input.lines() {
                let row: Json = serde_json::from_str(line).unwrap();
                assert!(
                    row.get(if task == "structure_tags" {
                        "case_id"
                    } else {
                        "case"
                    })
                    .is_some()
                );
            }
        }
        assert_eq!(CIF_FLAGS.len(), 32);
        assert!(mcs_params(&json!({"made_up":true})).is_err());
    }
}
