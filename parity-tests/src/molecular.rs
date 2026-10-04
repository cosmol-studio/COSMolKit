//! Public Rust execution and typed molecular observations. No chemistry algorithms.
use crate::registry::{Input, Record, SmilesCase, Value, molecule_plan::Profile};
use cosmolkit::{
    AddHsParams, Coordinate2DParams, KekulizeParams, Molecule, RemoveHsParams, SanitizeOperations,
    SanitizeParams, SmilesParseParams,
};
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtomRow {
    pub id: usize,
    pub atomic_number: u8,
    pub isotope: u16,
    pub charge: i8,
    pub explicit_h: u8,
    pub no_implicit: bool,
    pub radicals: u8,
    pub aromatic: bool,
    pub hybridization: i64,
    pub chiral: i64,
    pub permutation: u32,
    pub atom_map: u32,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BondRow {
    pub id: usize,
    pub begin: usize,
    pub end: usize,
    pub order: i64,
    pub aromatic: bool,
    pub conjugated: bool,
    pub direction: i64,
    pub stereo: i64,
    pub stereo_atoms: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Topology {
    pub atoms: Vec<AtomRow>,
    pub bonds: Vec<BondRow>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Stage {
    Parse,
    Preparation,
    Operation,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Outcome {
    Matrix {
        dimension: usize,
        values_bits: Vec<u64>,
    },
    Float64Bits(u64),
    Unsigned(u32),
    Text(String),
    Topology(Topology),
    Coordinates2d {
        topology: Topology,
        xy_bits: Vec<[u64; 2]>,
    },
    // Retained diagnostics, never counted as parity passes until source-backed
    // cross-language error categories are registered. Equal strings are not proof.
    Error {
        stage: Stage,
        detail: String,
    },
}

pub fn read_corpus(path: &Path) -> Result<Vec<SmilesCase>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    // One complete SMILES/CXSMILES record per line. Do not split whitespace,
    // trim labels, deduplicate molecules, or filter invalid chemistry.
    Ok(text
        .lines()
        .enumerate()
        .map(|(index, smiles)| SmilesCase {
            id: format!("line:{}", index + 1),
            smiles: smiles.into(),
        })
        .collect())
}

pub fn topology(mol: &Molecule) -> Topology {
    Topology {
        atoms: mol
            .atoms()
            .iter()
            .map(|a| AtomRow {
                id: a.id().index(),
                atomic_number: a.atomic_number(),
                isotope: a.isotope().unwrap_or(0),
                charge: a.formal_charge(),
                explicit_h: a.explicit_hydrogens(),
                no_implicit: a.no_implicit(),
                radicals: a.radical_electrons(),
                aromatic: a.is_aromatic(),
                hybridization: a.hybridization().rdkit_code(),
                chiral: a.chiral_tag().rdkit_code(),
                permutation: a.chiral_permutation().unwrap_or(0),
                atom_map: a.atom_map().unwrap_or(0),
            })
            .collect(),
        bonds: mol
            .bonds()
            .iter()
            .map(|b| BondRow {
                id: b.id().index(),
                begin: b.begin().index(),
                end: b.end().index(),
                order: b.order().rdkit_code(),
                aromatic: b.is_aromatic(),
                conjugated: b.is_conjugated(),
                direction: b.direction().rdkit_code(),
                stereo: b.stereo().rdkit_code(),
                stereo_atoms: b
                    .stereo_atoms()
                    .map(|ids| ids.map(|id| id.index()).to_vec())
                    .unwrap_or_default(),
            })
            .collect(),
    }
}

pub fn validate_output(profile: &Profile, output: &Outcome) -> Result<(), String> {
    use Profile::*;
    let valid = match (profile, output) {
        (
            DistanceMatrix { .. },
            Outcome::Matrix {
                dimension,
                values_bits,
            },
        ) => dimension.checked_mul(*dimension) == Some(values_bits.len()),
        (_, Outcome::Error { detail, .. }) => !detail.is_empty(),
        (MolecularWeight { .. } | ExactMolecularWeight { .. }, Outcome::Float64Bits(bits)) => {
            f64::from_bits(*bits).is_finite()
        }
        (MolecularFormula { .. }, Outcome::Text(_)) => true,
        (SvgDefault, Outcome::Text(svg)) => svg.contains("<svg") && svg.contains("</svg>"),
        (NumHeavyAtoms { .. }, Outcome::Unsigned(_)) => true,
        (TotalAtomCount { .. }, Outcome::Unsigned(_)) => true,
        (LipinskiHBA { .. } | LipinskiHBD { .. }, Outcome::Unsigned(_)) => true,
        (
            NumRings { .. }
            | NumHeterocycles { .. }
            | NumAromaticRings { .. }
            | NumSaturatedRings { .. }
            | NumAliphaticRings { .. }
            | NumAromaticHeterocycles { .. }
            | NumAromaticCarbocycles { .. }
            | NumAliphaticHeterocycles { .. }
            | NumAliphaticCarbocycles { .. }
            | NumSaturatedHeterocycles { .. }
            | NumSaturatedCarbocycles { .. },
            Outcome::Unsigned(_),
        ) => true,
        // CSP3 uses the existing exact-bits observation with a finite
        // [0,1] validation band.
        (FractionCSP3 { .. }, Outcome::Float64Bits(bits)) => {
            let value = f64::from_bits(*bits);
            (0.0..=1.0).contains(&value)
        }
        (Coordinates2dDefault, Outcome::Coordinates2d { topology, xy_bits }) => {
            valid_topology(topology)
                && topology.atoms.len() == xy_bits.len()
                && xy_bits
                    .iter()
                    .flatten()
                    .all(|bits| f64::from_bits(*bits).is_finite())
        }
        (
            SmilesRead { .. }
            | SanitizeAll
            | Kekulize { .. }
            | AddHydrogens { .. }
            | RemoveHydrogens { .. },
            Outcome::Topology(t),
        ) => valid_topology(t),
        _ => false,
    };
    if valid {
        Ok(())
    } else {
        Err("malformed or wrong-kind molecular reference result".into())
    }
}

fn valid_topology(t: &Topology) -> bool {
    t.atoms.iter().enumerate().all(|(i, a)| a.id == i)
        && t.bonds.iter().enumerate().all(|(i, b)| {
            b.id == i
                && b.begin < t.atoms.len()
                && b.end < t.atoms.len()
                && (b.stereo_atoms.is_empty() || b.stereo_atoms.len() == 2)
                && b.stereo_atoms.iter().all(|&a| a < t.atoms.len())
        })
}

pub fn run(input: &Input) -> Result<Record, String> {
    let Input::Molecular { case, profile } = input else {
        return Err("expected molecular input".into());
    };
    let mut stage = Stage::Parse;
    let result = (|| -> Result<Outcome, String> {
        let (sanitize, remove_hydrogens) = match profile {
            Profile::SmilesRead {
                sanitize,
                remove_hydrogens,
            } => (*sanitize, *remove_hydrogens),
            Profile::SanitizeAll => (false, false),
            Profile::NumHeavyAtoms { remove_hydrogens } => (true, *remove_hydrogens),
            Profile::TotalAtomCount { remove_hydrogens } => (true, *remove_hydrogens),
            Profile::LipinskiHBA { remove_hydrogens }
            | Profile::LipinskiHBD { remove_hydrogens }
            | Profile::FractionCSP3 { remove_hydrogens }
            | Profile::NumRings { remove_hydrogens }
            | Profile::NumHeterocycles { remove_hydrogens }
            | Profile::NumAromaticRings { remove_hydrogens }
            | Profile::NumSaturatedRings { remove_hydrogens }
            | Profile::NumAliphaticRings { remove_hydrogens }
            | Profile::NumAromaticHeterocycles { remove_hydrogens }
            | Profile::NumAromaticCarbocycles { remove_hydrogens }
            | Profile::NumAliphaticHeterocycles { remove_hydrogens }
            | Profile::NumAliphaticCarbocycles { remove_hydrogens }
            | Profile::NumSaturatedHeterocycles { remove_hydrogens }
            | Profile::NumSaturatedCarbocycles { remove_hydrogens } => (true, *remove_hydrogens),
            _ => (true, true),
        };
        let mol = Molecule::from_smiles_with_params(
            &case.smiles,
            &SmilesParseParams {
                sanitize,
                remove_hydrogens,
                allow_cxsmiles: true,
                strict_cxsmiles: true,
                parse_name: true,
                skip_cleanup: false,
                debug_parse: false,
                replacements: Default::default(),
            },
        )
        .map_err(|e| e.to_string())?;
        stage = Stage::Operation;
        use Profile::*;
        let transformed = match profile {
            SvgDefault => {
                return mol
                    .to_svg(300, 300)
                    .map(Outcome::Text)
                    .map_err(|e| e.to_string());
            }
            DistanceMatrix {
                use_bond_order,
                use_atom_weights,
            } => {
                let matrix = mol
                    .distance_matrix_with_params(&cosmolkit::DistanceMatrixParams {
                        use_bond_order: *use_bond_order,
                        use_atom_weights: *use_atom_weights,
                    })
                    .map_err(|e| e.to_string())?;
                return Ok(Outcome::Matrix {
                    dimension: matrix.dimension(),
                    values_bits: matrix.values().iter().map(|v| v.to_bits()).collect(),
                });
            }
            MolecularWeight { only_heavy } => {
                return mol
                    .molecular_weight_with_options(*only_heavy)
                    .map(|n| Outcome::Float64Bits(n.to_bits()))
                    .map_err(|e| e.to_string());
            }
            ExactMolecularWeight { only_heavy } => {
                return mol
                    .exact_molecular_weight_with_options(*only_heavy)
                    .map(|n| Outcome::Float64Bits(n.to_bits()))
                    .map_err(|e| e.to_string());
            }
            MolecularFormula {
                separate_isotopes,
                abbreviate_h_isotopes,
            } => {
                return mol
                    .molecular_formula_with_options(*separate_isotopes, *abbreviate_h_isotopes)
                    .map(Outcome::Text)
                    .map_err(|e| e.to_string());
            }
            NumHeavyAtoms { .. } => {
                return mol
                    .num_heavy_atoms()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            TotalAtomCount { .. } => {
                return mol
                    .total_atom_count()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            LipinskiHBA { .. } => {
                return mol
                    .lipinski_hba()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            LipinskiHBD { .. } => {
                return mol
                    .lipinski_hbd()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            FractionCSP3 { .. } => {
                return mol
                    .fraction_csp3()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            NumRings { .. } => {
                return mol
                    .num_rings()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumHeterocycles { .. } => {
                return mol
                    .num_heterocycles()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumAromaticRings { .. } => {
                return mol
                    .num_aromatic_rings()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumSaturatedRings { .. } => {
                return mol
                    .num_saturated_rings()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumAliphaticRings { .. } => {
                return mol
                    .num_aliphatic_rings()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumAromaticHeterocycles { .. } => {
                return mol
                    .num_aromatic_heterocycles()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumAromaticCarbocycles { .. } => {
                return mol
                    .num_aromatic_carbocycles()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumAliphaticHeterocycles { .. } => {
                return mol
                    .num_aliphatic_heterocycles()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumAliphaticCarbocycles { .. } => {
                return mol
                    .num_aliphatic_carbocycles()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumSaturatedHeterocycles { .. } => {
                return mol
                    .num_saturated_heterocycles()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumSaturatedCarbocycles { .. } => {
                return mol
                    .num_saturated_carbocycles()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            SmilesRead { .. } => mol,
            SanitizeAll => mol
                .sanitize_with_params(&SanitizeParams {
                    operations: SanitizeOperations::ALL,
                })
                .map_err(|e| e.to_string())?,
            Kekulize {
                clear_aromatic_flags,
            } => mol
                .with_kekulized_bonds_with_params(&KekulizeParams {
                    mark_atoms_bonds: *clear_aromatic_flags,
                    canonical: true,
                    max_backtracks: 100,
                })
                .map_err(|e| e.to_string())?,
            AddHydrogens { explicit_only } => mol
                .with_hydrogens_with_params(&add_hydrogens_params(*explicit_only))
                .map_err(|e| e.to_string())?,
            RemoveHydrogens { sanitize } => {
                stage = Stage::Preparation;
                let mol = mol
                    .with_hydrogens_with_params(&add_hydrogens_params(false))
                    .map_err(|e| e.to_string())?;
                stage = Stage::Operation;
                mol.without_hydrogens_with_params(&RemoveHsParams {
                    sanitize: *sanitize,
                    remove_degree_zero: false,
                    remove_higher_degrees: false,
                    remove_only_h_neighbors: false,
                    remove_isotopes: false,
                    remove_and_track_isotopes: false,
                    remove_dummy_neighbors: false,
                    remove_defining_bond_stereo: false,
                    remove_with_wedged_bond: true,
                    remove_with_query: false,
                    remove_mapped: true,
                    remove_in_sgroups: true,
                    show_warnings: true,
                    remove_nonimplicit: true,
                    update_explicit_count: false,
                    remove_hydrides: false,
                    remove_nontetrahedral_neighbors: false,
                })
                .map_err(|e| e.to_string())?
            }
            Coordinates2dDefault => {
                let result = mol
                    .with_2d_coordinates_with_params(&Coordinate2DParams {
                        coordinate_map: Default::default(),
                        canonical_orientation: false,
                        clear_existing_2d: true,
                        flips_per_sample: 0,
                        samples: 0,
                        sample_seed: 0,
                        permute_degree_four: false,
                        force_rdkit: false,
                        use_ring_templates: false,
                    })
                    .map_err(|e| e.to_string())?;
                let xy = result
                    .coordinates_2d()
                    .ok_or("2D result has no coordinate block")?;
                return Ok(Outcome::Coordinates2d {
                    topology: topology(&result),
                    xy_bits: xy.iter().map(|xy| xy.map(f64::to_bits)).collect(),
                });
            }
            _ => return Err("profile is not executable".into()),
        };
        Ok(Outcome::Topology(topology(&transformed)))
    })();
    Ok(Record {
        input: input.clone(),
        output: Value::Molecular(result.unwrap_or_else(|detail| Outcome::Error { stage, detail })),
    })
}

fn add_hydrogens_params(explicit_only: bool) -> AddHsParams {
    AddHsParams {
        explicit_only,
        add_coords: false,
        add_residue_info: false,
        skip_queries: false,
        only_on_atoms: None,
    }
}

pub fn matches(expected: &Outcome, actual: &Outcome) -> bool {
    match (expected, actual) {
        (Outcome::Error { .. }, _) | (_, Outcome::Error { .. }) => false,
        (
            Outcome::Coordinates2d {
                topology: a,
                xy_bits: ax,
            },
            Outcome::Coordinates2d {
                topology: b,
                xy_bits: bx,
            },
        ) => {
            a == b
                && ax.len() == bx.len()
                && ax.iter().flatten().zip(bx.iter().flatten()).all(|(a, b)| {
                    let (a, b) = (f64::from_bits(*a), f64::from_bits(*b));
                    a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-8
                })
        }
        _ => expected == actual,
    }
}

/// Only the four literal tool-identifier substitutions from the old SVG test.
pub fn svg_matches(expected: &Outcome, actual: &Outcome) -> bool {
    let (Outcome::Text(expected), Outcome::Text(actual)) = (expected, actual) else {
        return false;
    };
    fn normalize(svg: &str) -> String {
        svg.replace(
            "xmlns:rdkit='http://www.rdkit.org/xml'",
            "xmlns:tool='__tool_namespace__'",
        )
        .replace(
            "xmlns:cosmolkit='https://www.cosmol.org'",
            "xmlns:tool='__tool_namespace__'",
        )
        .replace("rdkit:", "tool:")
        .replace("cosmolkit:", "tool:")
    }
    normalize(expected) == normalize(actual)
}
