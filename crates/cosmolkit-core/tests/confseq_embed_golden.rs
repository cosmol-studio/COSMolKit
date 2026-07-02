use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{
    ChiralTag, ConfSeqDecodeOptions, ConfSeqTemplateBackend, EmbedParameters, Molecule,
    decode_confseq_with_options,
    io::sdf::{SdfCoordinateMode, SdfReadParams, read_sdf_from_str_with_params},
};
use serde::{Deserialize, Serialize};
use serde_json::Value;

mod common;
use common::parity_data;

const DISTANCE_MATRIX_TOLERANCE: f64 = 1.0e-6;
const COLLISION_THRESHOLD: f64 = 0.4;

#[derive(Debug, Deserialize)]
struct ConfSeqEmbedGoldenRecord {
    case_id: String,
    in_smiles: String,
    td_smiles: String,
    preset: String,
    random_seed: i32,
    optimize_with_uff: bool,
    rdkit_ok: bool,
    status: Option<i32>,
    coords: Option<Vec<[f64; 3]>>,
    sdf: Option<String>,
    pre_add_hs_summary: Option<Value>,
    with_h_summary: Option<Value>,
    error: Option<String>,
}

fn load_golden() -> Vec<ConfSeqEmbedGoldenRecord> {
    let explicit_path = std::env::var("COSMOLKIT_CONFSEQ_EMBED_GOLDEN").ok();
    let path = explicit_path
        .as_deref()
        .map(PathBuf::from)
        .unwrap_or_else(|| parity_data::golden_path("confseq_embed_template.jsonl"));
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate RDKit goldens with `{}` or set COSMOLKIT_CONFSEQ_EMBED_GOLDEN: {err}",
            path.display(),
            parity_data::regenerate_command()
        )
    });
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(idx, line)| {
            let line = line.unwrap_or_else(|err| {
                panic!("failed to read {} line {}: {err}", path.display(), idx + 1)
            });
            serde_json::from_str(&line).unwrap_or_else(|err| {
                panic!("failed to parse {} line {}: {err}", path.display(), idx + 1)
            })
        })
        .collect()
}

fn first_conformer_coords(molecule: &Molecule) -> &[[f64; 3]] {
    molecule
        .conformers_3d()
        .first()
        .expect("molecule should have a 3D conformer")
        .coordinates()
}

fn min_nonbonded_distance(molecule: &Molecule) -> Option<(usize, usize, f64)> {
    let coords = first_conformer_coords(molecule);
    let mut bonded = std::collections::HashSet::new();
    for bond in molecule.bonds() {
        let begin = bond.begin().index();
        let end = bond.end().index();
        bonded.insert((begin.min(end), begin.max(end)));
    }

    let mut best: Option<(usize, usize, f64)> = None;
    for i in 0..molecule.num_atoms() {
        for j in (i + 1)..molecule.num_atoms() {
            if bonded.contains(&(i, j)) {
                continue;
            }
            let d = distance(coords[i], coords[j]);
            if best.is_none_or(|(_, _, current)| d < current) {
                best = Some((i, j, d));
            }
        }
    }
    best
}

fn distance(left: [f64; 3], right: [f64; 3]) -> f64 {
    let dx = left[0] - right[0];
    let dy = left[1] - right[1];
    let dz = left[2] - right[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn distance_matrix_rmsd(left: &[[f64; 3]], right: &[[f64; 3]]) -> f64 {
    assert_eq!(left.len(), right.len());
    let mut sum_sq = 0.0;
    let mut count = 0_usize;
    for i in 0..left.len() {
        for j in (i + 1)..left.len() {
            let delta = distance(left[i], left[j]) - distance(right[i], right[j]);
            sum_sq += delta * delta;
            count += 1;
        }
    }
    (sum_sq / count as f64).sqrt()
}

#[derive(Debug, Clone, Serialize)]
struct MoleculeSummary {
    num_atoms: usize,
    num_bonds: usize,
    atoms: Vec<AtomSummary>,
    bonds: Vec<BondSummary>,
}

#[derive(Debug, Clone, Serialize)]
struct AtomSummary {
    idx: usize,
    atomic_num: u8,
    symbol: String,
    formal_charge: i8,
    explicit_hs: u8,
    isotope: u16,
    chiral_tag: &'static str,
    hybridization: String,
    is_aromatic: bool,
    no_implicit: bool,
}

#[derive(Debug, Clone, Serialize)]
struct BondSummary {
    idx: usize,
    begin: usize,
    end: usize,
    order: String,
    is_aromatic: bool,
}

fn molecule_summary(molecule: &Molecule) -> Value {
    let atoms = molecule
        .atoms()
        .iter()
        .map(|atom| AtomSummary {
            idx: atom.id().index(),
            atomic_num: atom.atomic_number(),
            symbol: element_symbol(atom.atomic_number()).to_string(),
            formal_charge: atom.formal_charge(),
            explicit_hs: atom.explicit_hydrogens(),
            isotope: atom.isotope().unwrap_or(0),
            chiral_tag: match atom.chiral_tag() {
                ChiralTag::Unspecified => "CHI_UNSPECIFIED",
                ChiralTag::TetrahedralCw => "CHI_TETRAHEDRAL_CW",
                ChiralTag::TetrahedralCcw => "CHI_TETRAHEDRAL_CCW",
                ChiralTag::Other => "CHI_OTHER",
                ChiralTag::Tetrahedral => "CHI_TETRAHEDRAL",
                ChiralTag::Allene => "CHI_ALLENE",
                ChiralTag::SquarePlanar => "CHI_SQUAREPLANAR",
                ChiralTag::TrigonalBipyramidal => "CHI_TRIGONALBIPYRAMIDAL",
                ChiralTag::Octahedral => "CHI_OCTAHEDRAL",
            },
            hybridization: format!("{:?}", atom.hybridization()).to_uppercase(),
            is_aromatic: atom.is_aromatic(),
            no_implicit: atom.no_implicit(),
        })
        .collect();
    let bonds = molecule
        .bonds()
        .iter()
        .map(|bond| BondSummary {
            idx: bond.id().index(),
            begin: bond.begin().index(),
            end: bond.end().index(),
            order: format!("{:?}", bond.order()).to_uppercase(),
            is_aromatic: bond.is_aromatic(),
        })
        .collect();
    serde_json::to_value(MoleculeSummary {
        num_atoms: molecule.num_atoms(),
        num_bonds: molecule.num_bonds(),
        atoms,
        bonds,
    })
    .expect("summary should serialize")
}

fn element_symbol(atomic_number: u8) -> &'static str {
    match atomic_number {
        0 => "*",
        1 => "H",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        16 => "S",
        17 => "Cl",
        35 => "Br",
        _ => "?",
    }
}

fn normalize_rdkit_summary(value: &Value) -> Value {
    let mut value = value.clone();
    if let Some(atoms) = value.get_mut("atoms").and_then(Value::as_array_mut) {
        for atom in atoms {
            if let Some(tag) = atom.get("chiral_tag").and_then(Value::as_str) {
                *atom.get_mut("chiral_tag").expect("tag exists") =
                    Value::String(tag.strip_prefix("CHI_").unwrap_or(tag).to_string());
            }
            if let Some(hybridization) = atom.get("hybridization").and_then(Value::as_str) {
                *atom.get_mut("hybridization").expect("hybridization exists") = Value::String(
                    hybridization
                        .strip_prefix("SP")
                        .map_or_else(|| hybridization.to_string(), |value| format!("SP{value}")),
                );
            }
        }
    }
    value
}

fn normalized_cosmolkit_summary(molecule: &Molecule) -> Value {
    let mut value = molecule_summary(molecule);
    if let Some(atoms) = value.get_mut("atoms").and_then(Value::as_array_mut) {
        for atom in atoms {
            if let Some(tag) = atom.get("chiral_tag").and_then(Value::as_str) {
                *atom.get_mut("chiral_tag").expect("tag exists") =
                    Value::String(tag.strip_prefix("CHI_").unwrap_or(tag).to_string());
            }
        }
    }
    value
}

fn first_json_diff(path: &str, left: &Value, right: &Value) -> Option<String> {
    match (left, right) {
        (Value::Object(left), Value::Object(right)) => {
            let mut keys = left.keys().chain(right.keys()).collect::<Vec<_>>();
            keys.sort();
            keys.dedup();
            for key in keys {
                let next = if path.is_empty() {
                    key.to_string()
                } else {
                    format!("{path}.{key}")
                };
                match (left.get(key), right.get(key)) {
                    (Some(left), Some(right)) => {
                        if let Some(diff) = first_json_diff(&next, left, right) {
                            return Some(diff);
                        }
                    }
                    _ => {
                        return Some(format!(
                            "{next}: left={:?} right={:?}",
                            left.get(key),
                            right.get(key)
                        ));
                    }
                }
            }
            None
        }
        (Value::Array(left), Value::Array(right)) => {
            let max_len = left.len().max(right.len());
            for idx in 0..max_len {
                let next = format!("{path}[{idx}]");
                match (left.get(idx), right.get(idx)) {
                    (Some(left), Some(right)) => {
                        if let Some(diff) = first_json_diff(&next, left, right) {
                            return Some(diff);
                        }
                    }
                    _ => {
                        return Some(format!(
                            "{next}: left={:?} right={:?}",
                            left.get(idx),
                            right.get(idx)
                        ));
                    }
                }
            }
            None
        }
        _ if left == right => None,
        _ => Some(format!("{path}: left={left:?} right={right:?}")),
    }
}

fn decode_stage(
    record: &ConfSeqEmbedGoldenRecord,
    apply_angles: bool,
    apply_dihedrals: bool,
) -> Molecule {
    let mut embed_params = EmbedParameters::etkdg();
    embed_params.random_seed = 0;
    let options = ConfSeqDecodeOptions {
        embed_params,
        optimize_with_uff: false,
        apply_angles,
        apply_dihedrals,
        template_backend: ConfSeqTemplateBackend::DistanceGeometry,
        ..ConfSeqDecodeOptions::default()
    };
    decode_confseq_with_options(&record.in_smiles, &record.td_smiles, &options).unwrap_or_else(
        |err| {
            let diag =
                cosmolkit_core::confseq::diagnostics::distance_geometry_dihedral_pair_diagnostics(
                    &record.in_smiles,
                    &record.td_smiles,
                    &options,
                )
                .ok();
            panic!(
                "{} stage decode failed: {err}; dihedral_pair_diag={diag:#?}",
                record.case_id
            )
        },
    )
}

fn stage_report(stage: &str, molecule: &Molecule, expected_coords: &[[f64; 3]]) -> String {
    let actual_coords = first_conformer_coords(molecule);
    let dmat_rmsd = distance_matrix_rmsd(actual_coords, expected_coords);
    let min_nonbond = min_nonbonded_distance(molecule);
    format!("{stage}: dmat_rmsd={dmat_rmsd:.6} min_nonbond={min_nonbond:?}")
}

#[test]
fn confseq_embedding_input_state_matches_rdkit_confseq_golden() {
    for (row_idx, record) in load_golden().iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let parsed = cosmolkit_core::confseq::diagnostics::parse_confseq_for_diagnostics(
            &record.in_smiles,
            &record.td_smiles,
        )
        .unwrap_or_else(|err| {
            panic!(
                "row {} ({}) diagnostic parse failed: {err}",
                row_idx + 1,
                record.case_id
            )
        });
        let pre_add_hs = cosmolkit_core::confseq::diagnostics::prepare_p_chiral_for_diagnostics(
            &parsed.stripped_smiles,
            &parsed.chiral_tags_by_atom,
        )
        .unwrap_or_else(|err| {
            panic!(
                "row {} ({}) diagnostic prepare failed: {err}",
                row_idx + 1,
                record.case_id
            )
        });
        let with_h = pre_add_hs.with_hydrogens().unwrap_or_else(|err| {
            panic!(
                "row {} ({}) COSMolKit AddHs failed: {err}",
                row_idx + 1,
                record.case_id
            )
        });

        let actual_pre = normalized_cosmolkit_summary(&pre_add_hs);
        let expected_pre =
            normalize_rdkit_summary(record.pre_add_hs_summary.as_ref().unwrap_or_else(|| {
                panic!(
                    "row {} ({}) missing RDKit pre-AddHs summary",
                    row_idx + 1,
                    record.case_id
                )
            }));
        assert!(
            actual_pre == expected_pre,
            "row {} ({}) pre-AddHs state differs from RDKit/ConfSeq: {}",
            row_idx + 1,
            record.case_id,
            first_json_diff("", &actual_pre, &expected_pre)
                .unwrap_or_else(|| "unknown difference".to_string())
        );
        let actual_with_h = normalized_cosmolkit_summary(&with_h);
        let expected_with_h =
            normalize_rdkit_summary(record.with_h_summary.as_ref().unwrap_or_else(|| {
                panic!(
                    "row {} ({}) missing RDKit AddHs summary",
                    row_idx + 1,
                    record.case_id
                )
            }));
        assert!(
            actual_with_h == expected_with_h,
            "row {} ({}) AddHs state differs from RDKit/ConfSeq: {}",
            row_idx + 1,
            record.case_id,
            first_json_diff("", &actual_with_h, &expected_with_h)
                .unwrap_or_else(|| "unknown difference".to_string())
        );
    }
}

#[test]
fn rdkit_confseq_embed_golden_template_is_structurally_sane() {
    for (row_idx, record) in load_golden().iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let Some(sdf) = &record.sdf else {
            panic!(
                "row {} ({}) RDKit success row missing SDF",
                row_idx + 1,
                record.case_id
            );
        };
        let expected = read_sdf_from_str_with_params(
            sdf,
            SdfReadParams {
                sanitize: true,
                remove_hs: false,
                coordinate_mode: SdfCoordinateMode::Require3D,
                ..Default::default()
            },
        )
        .unwrap_or_else(|err| {
            panic!(
                "row {} ({}) generated RDKit SDF should parse: {err}",
                row_idx + 1,
                record.case_id
            )
        })
        .molecule;
        if let Some((left, right, dist)) = min_nonbonded_distance(&expected) {
            assert!(
                dist > COLLISION_THRESHOLD,
                "row {} ({}) RDKit ConfSeq embed golden has nonbonded collision {left}-{right}: {dist:.6} A",
                row_idx + 1,
                record.case_id
            );
        }
    }
}

#[test]
#[ignore = "source-level ConfSeq ETKDG parity diagnostic; run explicitly when changing DistanceGeometry template decode"]
fn confseq_distance_geometry_template_matches_rdkit_confseq_embed_golden() {
    for (row_idx, record) in load_golden().iter().enumerate() {
        assert_eq!(record.preset, "ETKDG");
        assert_eq!(record.random_seed, 0);
        assert!(!record.optimize_with_uff);

        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}) failed RDKit golden row must include error",
                row_idx + 1,
                record.case_id
            );
            continue;
        }
        assert_eq!(
            record.status,
            Some(0),
            "row {} ({}) RDKit EmbedMolecule should succeed for this ConfSeq golden",
            row_idx + 1,
            record.case_id
        );

        let mut embed_params = EmbedParameters::etkdg();
        embed_params.random_seed = 0;
        let options = ConfSeqDecodeOptions {
            embed_params,
            optimize_with_uff: false,
            apply_angles: false,
            apply_dihedrals: false,
            template_backend: ConfSeqTemplateBackend::DistanceGeometry,
            ..ConfSeqDecodeOptions::default()
        };
        let actual = decode_confseq_with_options(&record.in_smiles, &record.td_smiles, &options)
            .unwrap_or_else(|err| {
                panic!(
                    "row {} ({}) COSMolKit ConfSeq template decode failed: {err}",
                    row_idx + 1,
                    record.case_id
                )
            });

        let expected_coords = record.coords.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) RDKit success row missing coords",
                row_idx + 1,
                record.case_id
            )
        });
        let actual_coords = first_conformer_coords(&actual);
        assert_eq!(
            actual_coords.len(),
            expected_coords.len(),
            "row {} ({}) atom count mismatch",
            row_idx + 1,
            record.case_id
        );
        let Some(sdf) = &record.sdf else {
            panic!(
                "row {} ({}) RDKit success row missing SDF",
                row_idx + 1,
                record.case_id
            );
        };
        let expected = read_sdf_from_str_with_params(
            sdf,
            SdfReadParams {
                sanitize: true,
                remove_hs: false,
                coordinate_mode: SdfCoordinateMode::Require3D,
                ..Default::default()
            },
        )
        .unwrap_or_else(|err| {
            panic!(
                "row {} ({}) generated RDKit SDF should parse: {err}",
                row_idx + 1,
                record.case_id
            )
        })
        .molecule;
        if let Some((left, right, dist)) = min_nonbonded_distance(&expected) {
            assert!(
                dist > COLLISION_THRESHOLD,
                "row {} ({}) RDKit ConfSeq embed golden has nonbonded collision {left}-{right}: {dist:.6} A",
                row_idx + 1,
                record.case_id
            );
        }
        let expected_min = min_nonbonded_distance(&expected);
        let actual_min = min_nonbonded_distance(&actual);
        let dmat_rmsd = distance_matrix_rmsd(actual_coords, expected_coords);
        let angle_only = decode_stage(record, true, false);
        let full = decode_stage(record, true, true);
        let stage_report = [
            stage_report("template", &actual, expected_coords),
            stage_report("angle_only", &angle_only, expected_coords),
            stage_report("full", &full, expected_coords),
        ]
        .join("; ");
        assert!(
            dmat_rmsd <= DISTANCE_MATRIX_TOLERANCE,
            "row {} ({}) ConfSeq template internal-distance RMSD mismatch: {dmat_rmsd:.6} A; expected_min_nonbond={expected_min:?}; actual_min_nonbond={actual_min:?}; stages={stage_report}",
            row_idx + 1,
            record.case_id
        );
    }
}
