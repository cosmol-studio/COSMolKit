use cosmolkit_tautomer::*;
use serde_json::{Value, json};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};
fn molecule_state(
    molecule: &TautomerRecord,
    coordinates: &cosmolkit_model::CoordinateBlock,
) -> Result<Value, Box<dyn std::error::Error>> {
    let atoms = molecule
        .topology
        .atoms
        .iter()
        .map(|atom| {
            Ok::<Value, Box<dyn std::error::Error>>(json!({
                "atomic_number": atom.atomic_number(),
                "formal_charge": atom.formal_charge(),
                "explicit_hydrogens": atom.explicit_hydrogens(),
                "no_implicit": atom.no_implicit(),
                "isotope": atom.isotope().unwrap_or(0),
                "radical_electrons": atom.radical_electrons(),
                "aromatic": atom.is_aromatic(),
                "chiral_tag": atom.chiral_tag().rdkit_name(),
                "hybridization": format!("{:?}", atom.hybridization()).to_ascii_uppercase(),
                "cip_code": atom.prop("_CIPCode").map(|v| v.as_string()).transpose()?,
            }))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bonds = molecule
        .topology.bonds
        .iter()
        .map(|bond| {
            let stereo_atoms = bond
                .stereo_atoms()
                .map(|atoms| vec![atoms[0].index(), atoms[1].index()])
                .unwrap_or_default();
            json!({
                "begin": bond.begin().index(),
                "end": bond.end().index(),
                "bond_type": bond.order().rdkit_name(),
                "aromatic": bond.is_aromatic(),
                "conjugated": bond.is_conjugated(),
                "direction": bond.direction().rdkit_name(),
                "stereo": bond.stereo().rdkit_name().strip_prefix("STEREO").expect("source enum name"),
                "stereo_atoms": stereo_atoms,
            })
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "isomeric_smiles": cosmolkit_smiles::write_smiles(cosmolkit_smiles::SmilesRecordView {topology:&molecule.topology,coordinates,properties:&molecule.properties})?,
        "atoms": atoms,
        "bonds": bonds,
    }))
}
fn options(branch: &Value) -> TautomerParams {
    TautomerParams::default()
        .with_max_tautomers(branch["max_tautomers"].as_u64().unwrap() as u32)
        .with_max_transforms(branch["max_transforms"].as_u64().unwrap() as u32)
        .with_remove_sp3_stereo(branch["remove_sp3_stereo"].as_bool().unwrap())
        .with_remove_bond_stereo(branch["remove_bond_stereo"].as_bool().unwrap())
        .with_remove_isotopic_hydrogens(branch["remove_isotopic_hydrogens"].as_bool().unwrap())
        .with_reassign_stereo(branch["reassign_stereo"].as_bool().unwrap())
}

pub fn validate_oracle(
    relative_path: &str,
    expected_rows: usize,
    expected_branches: usize,
    label: &str,
) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(relative_path);
    let rows = BufReader::new(File::open(path).unwrap())
        .lines()
        .map(|line| serde_json::from_str(&line.unwrap()).unwrap())
        .collect::<Vec<Value>>();
    assert_eq!(rows.len(), expected_rows);
    compare_rows(
        &rows,
        expected_branches,
        &Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target"),
        label,
    );
}

/// Compare owned, preflighted reference snapshots without rereading/generating data.
pub fn compare_rows(
    reference_rows: &[Value],
    expected_branches: usize,
    evidence_dir: &Path,
    label: &str,
) {
    let current = TautomerCatalog::current().unwrap();
    let v1 = TautomerCatalog::v1().unwrap();
    let mut failures = Vec::new();
    let evidence_path = evidence_dir.join(format!("TAU-{label}-detached-fields.jsonl"));
    let mut fields_file = std::io::BufWriter::new(File::create(evidence_path).unwrap());
    let mut checked = 0;
    let mut rows = 0;
    for row in reference_rows {
        rows += 1;
        let parsed = cosmolkit_smiles::parse_smiles(
            row["smiles"].as_str().unwrap(),
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: row["sanitize"].as_bool().unwrap(),
                remove_hydrogens: row["remove_hs"].as_bool().unwrap(),
                ..Default::default()
            },
        );
        if !row["parse"]["ok"].as_bool().unwrap() {
            assert!(parsed.is_err());
            continue;
        }
        let mut parsed = parsed.unwrap();
        let parse_params = cosmolkit_smiles::SmilesParseParams {
            sanitize: row["sanitize"].as_bool().unwrap(),
            remove_hydrogens: row["remove_hs"].as_bool().unwrap(),
            ..Default::default()
        };
        // Reuse the canonical facade's existing detached parse stage chain.
        let mut final_valence = None;
        let mut final_rings = None;
        if parse_params.remove_hydrogens {
            let assignment = cosmolkit_core::remove_hydrogens_with_params(
                parsed.topology,
                parsed.coordinates,
                parsed.properties,
                &cosmolkit_core::RemoveHsParams {
                    update_explicit_count: true,
                    sanitize: parse_params.sanitize,
                    ..Default::default()
                },
            )
            .unwrap();
            parsed = cosmolkit_smiles::SmilesRecord {
                topology: assignment.topology,
                coordinates: assignment.coordinates,
                properties: assignment.properties,
            };
            final_valence = assignment.final_valence;
            final_rings = assignment.final_rings;
        } else if parse_params.sanitize {
            let assignment =
                cosmolkit_core::sanitize_topology(&parsed.topology, &Default::default()).unwrap();
            parsed.topology = assignment.topology;
            final_valence = assignment.final_valence;
            final_rings = assignment.final_rings;
        }
        let parsed = cosmolkit_smiles::finalize_smiles_stereo(
            parsed,
            &parse_params,
            &mut final_valence,
            &mut final_rings,
        )
        .unwrap();
        let before = parsed.clone();
        for (branch, expected) in row["branches"].as_object().unwrap() {
            checked += 1;
            let cfg = &expected["parameters"];
            let params = options(cfg);
            let catalog = if cfg["catalog"] == "v1" {
                &v1
            } else {
                &current
            };
            let view = TautomerRecordView {
                topology: &parsed.topology,
                coordinates: &parsed.coordinates,
                properties: &parsed.properties,
                valence: if parse_params.sanitize {
                    final_valence.as_ref()
                } else {
                    None
                },
                rings: final_rings.as_ref(),
            };
            let actual = (|| -> Result<Value, Box<dyn std::error::Error>> {
                let result = enumerate_with_catalog(view, catalog, params, None)?;
                let canonical = pick_canonical_with(&result, &parsed.coordinates, |view| {
                    Ok(score_tautomer(view)?.total())
                })?;
                let mut scores = Vec::new();
                let mut states = Vec::new();
                for (_, candidate) in &result.entries {
                    let score = score_tautomer(candidate.view(&parsed.coordinates))?;
                    scores.push(json!({"ring":score.ring(),"substructure":score.substructure(),"hetero_hydrogen":score.hetero_hydrogen(),"total":score.total()}));
                    states.push(molecule_state(candidate, &parsed.coordinates)?);
                }
                Ok(
                    json!({"ordered_smiles":result.entries.iter().map(|(key,_)|key).collect::<Vec<_>>(),
                    "status":format!("{:?}",result.status),"modified_atoms":result.modified_atoms.iter().map(|id|id.index()).collect::<Vec<_>>(),
                    "modified_bonds":result.modified_bonds.iter().map(|id|id.index()).collect::<Vec<_>>(),"scores":scores,"molecule_states":states,
                    "canonical_smiles":cosmolkit_smiles::write_smiles(cosmolkit_smiles::SmilesRecordView {topology:&canonical.topology,coordinates:&parsed.coordinates,properties:&canonical.properties})?,
                    "canonical_state":molecule_state(&canonical,&parsed.coordinates)?}),
                )
            })();
            match actual {
                Ok(value) if expected["ok"]==true=>{
                    for (field,got) in value.as_object().unwrap() {
                        if got!=&expected[field] {
                            use std::io::Write;
                            serde_json::to_writer(&mut fields_file,&json!({"case":row["case_id"],"branch":branch,"field":field,"expected":expected[field],"actual":got})).unwrap();
                            fields_file.write_all(b"\n").unwrap();
                            failures.push(json!({"case":row["case_id"],"branch":branch,"field":field}));
                        }
                    }
                },
                Ok(_)=>failures.push(json!({"case":row["case_id"],"branch":branch,"failure":"expected structural error"})),
                Err(error)=>failures.push(json!({"case":row["case_id"],"branch":branch,"error":error.to_string()})),
            }
            assert_eq!(parsed, before, "source changed");
        }
    }
    std::io::Write::flush(&mut fields_file).unwrap();
    assert_eq!(rows, reference_rows.len());
    assert_eq!(checked, expected_branches);
    let out = evidence_dir.join(format!("TAU-{label}-detached-failures.json"));
    std::fs::write(
        &out,
        serde_json::to_vec_pretty(&json!({"rows":rows,"branches":checked,"failures":failures}))
            .unwrap(),
    )
    .unwrap();
    assert!(
        failures.is_empty(),
        "{} differences across {checked} source branches; full evidence in {}; first: {:?}",
        failures.len(),
        out.display(),
        failures
            .first()
            .map(|f| (&f["case"], &f["branch"], &f["field"], &f["error"]))
    );
}
