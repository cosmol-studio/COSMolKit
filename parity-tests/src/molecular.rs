//! Public Rust execution and typed molecular observations. No chemistry algorithms.
use crate::registry::{
    Input, Record, SmilesCase, Value,
    molecule_plan::{MorganInvariantKind, MorganOutputKind, Profile},
};
use cosmolkit::{
    AddHsParams, Coordinate2DParams, FingerprintAdditionalOutput, KekulizeParams, Molecule,
    MorganFingerprintParams, MorganInvariants, MorganParams, RemoveHsParams, SanitizeOperations,
    SanitizeParams, SmilesParseParams,
};
use serde::de::{self, MapAccess, SeqAccess, Visitor};
use serde::{Deserialize, Serialize};
use std::fmt;
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

/// The five optional FingerprintAdditionalOutput values emitted by the pinned Morgan
/// generator. Map entries are represented as key-sorted rows; each nested
/// source vector remains in its original order and may contain duplicates.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct MorganAdditionalOutput {
    pub atom_counts: Option<Vec<u32>>,
    pub atom_to_bits: Option<Vec<Vec<u64>>>,
    pub bit_info_map: Option<Vec<(u64, Vec<(u32, u32)>)>>,
    pub bit_paths: Option<Vec<(u64, Vec<Vec<i32>>)>>,
    pub atoms_per_bit: Option<Vec<(u64, Vec<Vec<i32>>)>>,
}

const MORGAN_ADDITIONAL_OUTPUT_FIELDS: &[&str] = &[
    "atom_counts",
    "atom_to_bits",
    "bit_info_map",
    "bit_paths",
    "atoms_per_bit",
];

struct MorganAdditionalOutputVisitor;

impl<'de> Visitor<'de> for MorganAdditionalOutputVisitor {
    type Value = MorganAdditionalOutput;

    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("a Morgan FingerprintAdditionalOutput object with all five fields")
    }

    fn visit_map<A>(self, mut map: A) -> Result<Self::Value, A::Error>
    where
        A: MapAccess<'de>,
    {
        let mut atom_counts = None;
        let mut atom_to_bits = None;
        let mut bit_info_map = None;
        let mut bit_paths = None;
        let mut atoms_per_bit = None;

        while let Some(field) = map.next_key::<String>()? {
            match field.as_str() {
                "atom_counts" => {
                    if atom_counts.is_some() {
                        return Err(de::Error::duplicate_field("atom_counts"));
                    }
                    atom_counts = Some(map.next_value()?);
                }
                "atom_to_bits" => {
                    if atom_to_bits.is_some() {
                        return Err(de::Error::duplicate_field("atom_to_bits"));
                    }
                    atom_to_bits = Some(map.next_value()?);
                }
                "bit_info_map" => {
                    if bit_info_map.is_some() {
                        return Err(de::Error::duplicate_field("bit_info_map"));
                    }
                    bit_info_map = Some(map.next_value()?);
                }
                "bit_paths" => {
                    if bit_paths.is_some() {
                        return Err(de::Error::duplicate_field("bit_paths"));
                    }
                    bit_paths = Some(map.next_value()?);
                }
                "atoms_per_bit" => {
                    if atoms_per_bit.is_some() {
                        return Err(de::Error::duplicate_field("atoms_per_bit"));
                    }
                    atoms_per_bit = Some(map.next_value()?);
                }
                _ => {
                    return Err(de::Error::unknown_field(
                        &field,
                        MORGAN_ADDITIONAL_OUTPUT_FIELDS,
                    ));
                }
            }
        }

        Ok(MorganAdditionalOutput {
            atom_counts: atom_counts.ok_or_else(|| de::Error::missing_field("atom_counts"))?,
            atom_to_bits: atom_to_bits.ok_or_else(|| de::Error::missing_field("atom_to_bits"))?,
            bit_info_map: bit_info_map.ok_or_else(|| de::Error::missing_field("bit_info_map"))?,
            bit_paths: bit_paths.ok_or_else(|| de::Error::missing_field("bit_paths"))?,
            atoms_per_bit: atoms_per_bit
                .ok_or_else(|| de::Error::missing_field("atoms_per_bit"))?,
        })
    }

    fn visit_seq<A>(self, mut sequence: A) -> Result<Self::Value, A::Error>
    where
        A: SeqAccess<'de>,
    {
        let atom_counts: Option<Vec<u32>> = sequence
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(0, &self))?;
        let atom_to_bits: Option<Vec<Vec<u64>>> = sequence
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(1, &self))?;
        let bit_info_map: Option<Vec<(u64, Vec<(u32, u32)>)>> = sequence
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(2, &self))?;
        let bit_paths: Option<Vec<(u64, Vec<Vec<i32>>)>> = sequence
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(3, &self))?;
        let atoms_per_bit: Option<Vec<(u64, Vec<Vec<i32>>)>> = sequence
            .next_element()?
            .ok_or_else(|| de::Error::invalid_length(4, &self))?;
        if sequence.next_element::<de::IgnoredAny>()?.is_some() {
            return Err(de::Error::invalid_length(6, &self));
        }

        Ok(MorganAdditionalOutput {
            atom_counts,
            atom_to_bits,
            bit_info_map,
            bit_paths,
            atoms_per_bit,
        })
    }
}

impl<'de> Deserialize<'de> for MorganAdditionalOutput {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        deserializer.deserialize_struct(
            "MorganAdditionalOutput",
            MORGAN_ADDITIONAL_OUTPUT_FIELDS,
            MorganAdditionalOutputVisitor,
        )
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MorganDenseBitsOutput {
    pub length: u32,
    pub on_bits: Vec<u32>,
    pub additional_output: MorganAdditionalOutput,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MorganSparseBitsOutput {
    pub length: u32,
    /// Raw signed `int` observations from SparseBitVect's set<int> storage.
    pub on_bits: Vec<i32>,
    pub additional_output: MorganAdditionalOutput,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MorganHashedCountsOutput {
    pub length: u32,
    pub entries: Vec<(u32, i32)>,
    pub additional_output: MorganAdditionalOutput,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MorganSparseCountsOutput {
    pub length: u64,
    pub entries: Vec<(u64, i32)>,
    pub additional_output: MorganAdditionalOutput,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Outcome {
    TautomerEnumeration(crate::tautomer::Enumeration),
    TautomerCanonicalization(crate::tautomer::Canonicalization),
    Matrix {
        dimension: usize,
        values_bits: Vec<u64>,
    },
    MorganDenseBits(MorganDenseBitsOutput),
    MorganSparseBits(MorganSparseBitsOutput),
    MorganHashedCounts(MorganHashedCountsOutput),
    MorganSparseCounts(MorganSparseCountsOutput),
    Float64Bits(u64),
    Float64PairBits {
        first: u64,
        second: u64,
    },
    Float64VectorBits(Vec<u64>),
    Float64VectorsBits(Vec<Vec<u64>>),
    LabuteAsaContributionsBits {
        asa: u64,
        atom_contributions: Vec<u64>,
        hydrogen_contribution: u64,
    },
    Float64ContributionsBits {
        value: u64,
        atom_contributions: Vec<u64>,
    },
    UnsignedVector(Vec<u32>),
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
        (TautomerEnumeration { .. }, Outcome::TautomerEnumeration(value)) => {
            crate::tautomer::valid_enumeration(value)
        }
        (TautomerCanonicalization { .. }, Outcome::TautomerCanonicalization(value)) => {
            crate::tautomer::valid_canonicalization(value)
        }
        (
            DistanceMatrix { .. },
            Outcome::Matrix {
                dimension,
                values_bits,
            },
        ) => dimension.checked_mul(*dimension) == Some(values_bits.len()),
        (
            Morgan {
                output: crate::registry::molecule_plan::MorganOutputKind::DenseBits,
                ..
            },
            Outcome::MorganDenseBits(value),
        ) => {
            value.length == MORGAN_CONFIGURED_SIZE
                && value.on_bits.iter().all(|&bit| bit < value.length)
                && strictly_increasing(&value.on_bits)
                && valid_morgan_additional_output(&value.additional_output)
        }
        (
            Morgan {
                output: crate::registry::molecule_plan::MorganOutputKind::SparseBits,
                ..
            },
            Outcome::MorganSparseBits(value),
        ) => {
            // The source length is the complete u32 domain. Its set<int>
            // observation is signed-ascending, so negative rows are retained.
            value.length == MORGAN_SPARSE_BITS_LENGTH
                && strictly_increasing(&value.on_bits)
                && valid_morgan_additional_output(&value.additional_output)
        }
        (
            Morgan {
                output: crate::registry::molecule_plan::MorganOutputKind::HashedCounts,
                ..
            },
            Outcome::MorganHashedCounts(value),
        ) => {
            value.length == MORGAN_CONFIGURED_SIZE
                && value.entries.iter().all(|&(index, _)| index < value.length)
                && strictly_increasing_keys(&value.entries)
                && valid_morgan_additional_output(&value.additional_output)
        }
        (
            Morgan {
                output: crate::registry::molecule_plan::MorganOutputKind::SparseCounts,
                ..
            },
            Outcome::MorganSparseCounts(value),
        ) => {
            // SparseIntVect permits the terminal index when length is the
            // maximum u64 value; every u64 index is therefore source-valid.
            value.length == MORGAN_SPARSE_COUNTS_LENGTH
                && strictly_increasing_keys(&value.entries)
                && valid_morgan_additional_output(&value.additional_output)
        }
        (_, Outcome::Error { detail, .. }) => !detail.is_empty(),
        (MolecularWeight { .. } | ExactMolecularWeight { .. }, Outcome::Float64Bits(bits)) => {
            f64::from_bits(*bits).is_finite()
        }
        (
            Chi0
            | Chi1
            | HallKierAlpha
            | Kappa1
            | Kappa2
            | Kappa3
            | Phi
            | Chi0V
            | Chi1V
            | Chi2V
            | Chi3V
            | Chi4V
            | Chi0N
            | Chi1N
            | Chi2N
            | Chi3N
            | Chi4N
            | ChiNV { .. }
            | ChiNN { .. },
            Outcome::Float64Bits(bits),
        ) => f64::from_bits(*bits).is_finite(),
        (
            HallKierAlphaWithContributions,
            Outcome::Float64ContributionsBits {
                value,
                atom_contributions,
            },
        ) => {
            f64::from_bits(*value).is_finite()
                && atom_contributions
                    .iter()
                    .all(|b| f64::from_bits(*b).is_finite())
        }
        (
            NumAmideBonds
            | NumSpiroAtoms
            | NumBridgeheadAtoms
            | NumAtomStereoCenters
            | NumUnspecifiedAtomStereoCenters
            | NumRotatableBonds { .. },
            Outcome::Unsigned(_),
        ) => true,
        (
            LabuteAsa { .. }
            | Tpsa { .. }
            | SlogpVsa1
            | SlogpVsa2
            | SlogpVsa3
            | SlogpVsa4
            | SlogpVsa5
            | SlogpVsa6
            | SlogpVsa7
            | SlogpVsa8
            | SlogpVsa9
            | SlogpVsa10
            | SlogpVsa11
            | SlogpVsa12
            | SmrVsa1
            | SmrVsa2
            | SmrVsa3
            | SmrVsa4
            | SmrVsa5
            | SmrVsa6
            | SmrVsa7
            | SmrVsa8
            | SmrVsa9
            | SmrVsa10
            | Qed
            | Chi0VWithParams { .. }
            | Chi1VWithParams { .. }
            | Chi2VWithParams { .. }
            | Chi3VWithParams { .. }
            | Chi4VWithParams { .. }
            | Chi0NWithParams { .. }
            | Chi1NWithParams { .. }
            | Chi2NWithParams { .. }
            | Chi3NWithParams { .. }
            | Chi4NWithParams { .. }
            | ChiNVWithParams { .. }
            | ChiNNWithParams { .. },
            Outcome::Float64Bits(bits),
        ) => f64::from_bits(*bits).is_finite(),
        (CrippenDescriptors { .. }, Outcome::Float64PairBits { first, second }) => [first, second]
            .into_iter()
            .all(|b| f64::from_bits(*b).is_finite()),
        (
            LabuteAsaContributions { .. },
            Outcome::LabuteAsaContributionsBits {
                asa,
                atom_contributions,
                hydrogen_contribution,
            },
        ) => [asa, hydrogen_contribution]
            .into_iter()
            .chain(atom_contributions)
            .all(|b| f64::from_bits(*b).is_finite()),
        (SlogpVsa { bins, force } | SmrVsa { bins, force }, Outcome::Float64VectorBits(values)) => {
            let expected = match bins {
                crate::registry::molecule_plan::VsaBins::Default => {
                    if matches!(profile, SlogpVsa { .. }) {
                        12
                    } else {
                        10
                    }
                }
                crate::registry::molecule_plan::VsaBins::CustomDuplicates => 6,
            };
            (force.is_some() || matches!(bins, crate::registry::molecule_plan::VsaBins::Default))
                && values.len() == expected
                && values.iter().all(|b| f64::from_bits(*b).is_finite())
        }
        (LabuteAsaCacheSequence, Outcome::Float64VectorBits(values)) => {
            values.len() == 4 && values.iter().all(|b| f64::from_bits(*b).is_finite())
        }
        (ChiNVCacheSequence | ChiNNCacheSequence, Outcome::Float64VectorsBits(values)) => {
            values.len() == 3
                && values
                    .iter()
                    .all(|v| v.len() == 7 && v.iter().all(|b| f64::from_bits(*b).is_finite()))
        }
        (SlogpVsaCacheSequence, Outcome::Float64VectorsBits(values)) => {
            values.iter().map(Vec::len).collect::<Vec<_>>() == [12, 12, 6]
                && values
                    .iter()
                    .flatten()
                    .all(|b| f64::from_bits(*b).is_finite())
        }
        (SmrVsaCacheSequence, Outcome::Float64VectorsBits(values)) => {
            values.iter().map(Vec::len).collect::<Vec<_>>() == [10, 6]
                && values
                    .iter()
                    .flatten()
                    .all(|b| f64::from_bits(*b).is_finite())
        }
        (Mqns { .. }, Outcome::UnsignedVector(values)) => values.len() == 42,
        (MolecularFormula { .. }, Outcome::Text(_)) => true,
        (SvgDefault, Outcome::Text(svg)) => svg.contains("<svg") && svg.contains("</svg>"),
        (NumHeavyAtoms { .. }, Outcome::Unsigned(_)) => true,
        (TotalAtomCount { .. }, Outcome::Unsigned(_)) => true,
        (LipinskiHBA { .. } | LipinskiHBD { .. }, Outcome::Unsigned(_)) => true,
        (NumHeteroatoms { .. }, Outcome::Unsigned(_)) => true,
        (NumHba { .. }, Outcome::Unsigned(_)) => true,
        (NumHbd { .. }, Outcome::Unsigned(_)) => true,
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

const MORGAN_CONFIGURED_SIZE: u32 = 2048;
const MORGAN_SPARSE_BITS_LENGTH: u32 = u32::MAX;
const MORGAN_SPARSE_COUNTS_LENGTH: u64 = u64::MAX;

fn strictly_increasing<T: Ord>(values: &[T]) -> bool {
    values.windows(2).all(|pair| pair[0] < pair[1])
}

fn strictly_increasing_keys<K: Ord, V>(rows: &[(K, V)]) -> bool {
    rows.windows(2).all(|pair| {
        let (left, _) = &pair[0];
        let (right, _) = &pair[1];
        left < right
    })
}

fn valid_morgan_additional_output(output: &MorganAdditionalOutput) -> bool {
    let per_atom_lengths_match = match (&output.atom_counts, &output.atom_to_bits) {
        (Some(counts), Some(bits)) => counts.len() == bits.len(),
        _ => true,
    };
    per_atom_lengths_match
        && output
            .bit_info_map
            .as_ref()
            .is_none_or(|rows| strictly_increasing_keys(rows))
        && output
            .bit_paths
            .as_ref()
            .is_none_or(|rows| strictly_increasing_keys(rows))
        && output
            .atoms_per_bit
            .as_ref()
            .is_none_or(|rows| strictly_increasing_keys(rows))
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
            | Profile::NumHeteroatoms { remove_hydrogens }
            | Profile::NumHba { remove_hydrogens }
            | Profile::NumHbd { remove_hydrogens }
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
            TautomerEnumeration { parameters } => {
                return crate::tautomer::enumerate(&mol, *parameters);
            }
            TautomerCanonicalization { parameters } => {
                return crate::tautomer::canonicalize(&mol, *parameters);
            }
            Chi0 => {
                return mol
                    .chi_0()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi1 => {
                return mol
                    .chi_1()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            HallKierAlpha => {
                return mol
                    .hall_kier_alpha()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            HallKierAlphaWithContributions => {
                return mol
                    .hall_kier_alpha_with_contributions()
                    .map(
                        |(value, atom_contributions)| Outcome::Float64ContributionsBits {
                            value: value.to_bits(),
                            atom_contributions: atom_contributions
                                .into_iter()
                                .map(f64::to_bits)
                                .collect(),
                        },
                    )
                    .map_err(|e| e.to_string());
            }
            Kappa1 => {
                return mol
                    .kappa_1()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Kappa2 => {
                return mol
                    .kappa_2()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Kappa3 => {
                return mol
                    .kappa_3()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Phi => {
                return mol
                    .phi()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Mqns { force } => {
                return mol
                    .mqns(*force)
                    .map(Outcome::UnsignedVector)
                    .map_err(|e| e.to_string());
            }
            Chi0V => {
                return mol
                    .chi_0_v()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi1V => {
                return mol
                    .chi_1_v()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi2V => {
                return mol
                    .chi_2_v()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi3V => {
                return mol
                    .chi_3_v()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi4V => {
                return mol
                    .chi_4_v()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi0N => {
                return mol
                    .chi_0_n()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi1N => {
                return mol
                    .chi_1_n()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi2N => {
                return mol
                    .chi_2_n()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi3N => {
                return mol
                    .chi_3_n()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi4N => {
                return mol
                    .chi_4_n()
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            ChiNV { order } => {
                return mol
                    .chi_n_v(*order)
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }
            ChiNN { order } => {
                return mol
                    .chi_n_n(*order)
                    .map(|value| Outcome::Float64Bits(value.to_bits()))
                    .map_err(|e| e.to_string());
            }

            NumAmideBonds => {
                return mol
                    .num_amide_bonds()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumSpiroAtoms => {
                return mol
                    .num_spiro_atoms()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumBridgeheadAtoms => {
                return mol
                    .num_bridgehead_atoms()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumAtomStereoCenters => {
                return mol
                    .num_atom_stereo_centers()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumUnspecifiedAtomStereoCenters => {
                return mol
                    .num_unspecified_atom_stereo_centers()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumRotatableBonds { mode } => {
                return match mode {
                    None => mol.num_rotatable_bonds(),
                    Some(mode) => {
                        use crate::registry::molecule_plan::RotatableBondMode as M;
                        let params = match mode {
                            M::Default => cosmolkit::RotatableBondsOptions::Default,
                            M::NonStrict => cosmolkit::RotatableBondsOptions::NonStrict,
                            M::Strict => cosmolkit::RotatableBondsOptions::Strict,
                            M::StrictLinkages => cosmolkit::RotatableBondsOptions::StrictLinkages,
                        };
                        mol.num_rotatable_bonds_with_params(&params)
                    }
                }
                .map(Outcome::Unsigned)
                .map_err(|e| e.to_string());
            }
            CrippenDescriptors {
                include_hydrogens,
                force,
            } => {
                return (match include_hydrogens {
                    None => mol.crippen_descriptors(),
                    Some(flag) => mol.crippen_descriptors_with_params(*flag, *force),
                })
                .map(|v| Outcome::Float64PairBits {
                    first: v.logp.to_bits(),
                    second: v.molar_refractivity.to_bits(),
                })
                .map_err(|e| e.to_string());
            }
            LabuteAsa {
                include_hydrogens,
                force,
            } => {
                return (match include_hydrogens {
                    None => mol.labute_asa(),
                    Some(flag) => mol.labute_asa_with_params(*flag, *force),
                })
                .map(|v| Outcome::Float64Bits(v.to_bits()))
                .map_err(|e| e.to_string());
            }
            LabuteAsaContributions {
                include_hydrogens,
                force,
            } => {
                return (match include_hydrogens {
                    None => mol.labute_asa_contributions(),
                    Some(flag) => mol.labute_asa_contributions_with_params(*flag, *force),
                })
                .map(|v| Outcome::LabuteAsaContributionsBits {
                    asa: v.asa.to_bits(),
                    atom_contributions: v.atom_contributions.iter().map(|x| x.to_bits()).collect(),
                    hydrogen_contribution: v.hydrogen_contribution.to_bits(),
                })
                .map_err(|e| e.to_string());
            }
            Tpsa {
                include_sulfur_phosphorus,
                force,
            } => {
                return (match include_sulfur_phosphorus {
                    None => mol.tpsa(),
                    Some(flag) => mol.tpsa_with_params(*flag, *force),
                })
                .map(|v| Outcome::Float64Bits(v.to_bits()))
                .map_err(|e| e.to_string());
            }
            SlogpVsa { bins, force } => {
                return (match force {
                    None => mol.slogp_vsa(),
                    Some(flag) => mol.slogp_vsa_with_params(bins.values(), *flag),
                })
                .map(|v| Outcome::Float64VectorBits(v.iter().map(|x| x.to_bits()).collect()))
                .map_err(|e| e.to_string());
            }
            SmrVsa { bins, force } => {
                return (match force {
                    None => mol.smr_vsa(),
                    Some(flag) => mol.smr_vsa_with_params(bins.values(), *flag),
                })
                .map(|v| Outcome::Float64VectorBits(v.iter().map(|x| x.to_bits()).collect()))
                .map_err(|e| e.to_string());
            }
            SlogpVsa1 => {
                return mol
                    .slogp_vsa_1()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa2 => {
                return mol
                    .slogp_vsa_2()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa3 => {
                return mol
                    .slogp_vsa_3()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa4 => {
                return mol
                    .slogp_vsa_4()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa5 => {
                return mol
                    .slogp_vsa_5()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa6 => {
                return mol
                    .slogp_vsa_6()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa7 => {
                return mol
                    .slogp_vsa_7()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa8 => {
                return mol
                    .slogp_vsa_8()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa9 => {
                return mol
                    .slogp_vsa_9()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa10 => {
                return mol
                    .slogp_vsa_10()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa11 => {
                return mol
                    .slogp_vsa_11()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SlogpVsa12 => {
                return mol
                    .slogp_vsa_12()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SmrVsa1 => {
                return mol
                    .smr_vsa_1()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SmrVsa2 => {
                return mol
                    .smr_vsa_2()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SmrVsa3 => {
                return mol
                    .smr_vsa_3()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SmrVsa4 => {
                return mol
                    .smr_vsa_4()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SmrVsa5 => {
                return mol
                    .smr_vsa_5()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SmrVsa6 => {
                return mol
                    .smr_vsa_6()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SmrVsa7 => {
                return mol
                    .smr_vsa_7()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SmrVsa8 => {
                return mol
                    .smr_vsa_8()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SmrVsa9 => {
                return mol
                    .smr_vsa_9()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            SmrVsa10 => {
                return mol
                    .smr_vsa_10()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Qed => {
                return mol
                    .qed()
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi0VWithParams { force } => {
                return mol
                    .chi_0_v_with_params(*force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi1VWithParams { force } => {
                return mol
                    .chi_1_v_with_params(*force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi2VWithParams { force } => {
                return mol
                    .chi_2_v_with_params(*force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi3VWithParams { force } => {
                return mol
                    .chi_3_v_with_params(*force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi4VWithParams { force } => {
                return mol
                    .chi_4_v_with_params(*force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi0NWithParams { force } => {
                return mol
                    .chi_0_n_with_params(*force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi1NWithParams { force } => {
                return mol
                    .chi_1_n_with_params(*force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi2NWithParams { force } => {
                return mol
                    .chi_2_n_with_params(*force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi3NWithParams { force } => {
                return mol
                    .chi_3_n_with_params(*force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            Chi4NWithParams { force } => {
                return mol
                    .chi_4_n_with_params(*force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            ChiNVWithParams { order, force } => {
                return mol
                    .chi_n_v_with_params(*order, *force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            ChiNNWithParams { order, force } => {
                return mol
                    .chi_n_n_with_params(*order, *force)
                    .map(|v| Outcome::Float64Bits(v.to_bits()))
                    .map_err(|e| e.to_string());
            }
            LabuteAsaCacheSequence => {
                let values = [(false, false), (true, false), (true, true), (false, false)]
                    .into_iter()
                    .map(|(include, force)| {
                        mol.labute_asa_with_params(include, force)
                            .map(f64::to_bits)
                            .map_err(|e| e.to_string())
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                return Ok(Outcome::Float64VectorBits(values));
            }
            SlogpVsaCacheSequence => {
                let cold = mol.slogp_vsa().map_err(|e| e.to_string())?;
                let warm = mol.slogp_vsa().map_err(|e| e.to_string())?;
                let forced = mol
                    .slogp_vsa_with_params(Some(&[-0.2, 0.0, 0.25, 0.25, 0.8]), true)
                    .map_err(|e| e.to_string())?;
                return Ok(Outcome::Float64VectorsBits(
                    [cold, warm, forced]
                        .into_iter()
                        .map(|v| v.into_iter().map(f64::to_bits).collect())
                        .collect(),
                ));
            }
            SmrVsaCacheSequence => {
                mol.slogp_vsa().map_err(|e| e.to_string())?;
                mol.slogp_vsa_with_params(Some(&[-0.2, 0.0, 0.25, 0.25, 0.8]), true)
                    .map_err(|e| e.to_string())?;
                let warm = mol.smr_vsa().map_err(|e| e.to_string())?;
                let forced = mol
                    .smr_vsa_with_params(Some(&[-0.2, 0.0, 0.25, 0.25, 0.8]), true)
                    .map_err(|e| e.to_string())?;
                return Ok(Outcome::Float64VectorsBits(
                    [warm, forced]
                        .into_iter()
                        .map(|v| v.into_iter().map(f64::to_bits).collect())
                        .collect(),
                ));
            }
            ChiNVCacheSequence | ChiNNCacheSequence => {
                let mut values = Vec::new();
                for force in [false, false, true] {
                    let row = (0..=6)
                        .map(|order| {
                            if matches!(profile, ChiNVCacheSequence) {
                                mol.chi_n_v_with_params(order, force)
                            } else {
                                mol.chi_n_n_with_params(order, force)
                            }
                            .map(f64::to_bits)
                            .map_err(|e| e.to_string())
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    values.push(row);
                }
                return Ok(Outcome::Float64VectorsBits(values));
            }
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
                    .molecular_weight_with_params(*only_heavy)
                    .map(|n| Outcome::Float64Bits(n.to_bits()))
                    .map_err(|e| e.to_string());
            }
            ExactMolecularWeight { only_heavy } => {
                return mol
                    .exact_molecular_weight_with_params(*only_heavy)
                    .map(|n| Outcome::Float64Bits(n.to_bits()))
                    .map_err(|e| e.to_string());
            }
            MolecularFormula {
                separate_isotopes,
                abbreviate_h_isotopes,
            } => {
                return mol
                    .molecular_formula_with_params(*separate_isotopes, *abbreviate_h_isotopes)
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
            NumHeteroatoms { .. } => {
                return mol
                    .num_heteroatoms()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumHba { .. } => {
                return mol
                    .num_hba()
                    .map(Outcome::Unsigned)
                    .map_err(|e| e.to_string());
            }
            NumHbd { .. } => {
                return mol
                    .num_hbd()
                    .map(Outcome::Unsigned)
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
            Morgan {
                output,
                radius,
                include_chirality,
                invariants,
                count_simulation,
            } => {
                return run_morgan(
                    &mol,
                    *output,
                    *radius,
                    *include_chirality,
                    *invariants,
                    *count_simulation,
                );
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

fn run_morgan(
    molecule: &Molecule,
    output_kind: MorganOutputKind,
    radius: u32,
    include_chirality: bool,
    invariant_kind: MorganInvariantKind,
    count_simulation: bool,
) -> Result<Outcome, String> {
    let params = MorganFingerprintParams {
        generator: MorganParams {
            radius,
            include_chirality,
            count_simulation,
            ..MorganParams::default()
        },
        invariants: match invariant_kind {
            MorganInvariantKind::Connectivity => MorganInvariants::Connectivity,
            MorganInvariantKind::Features => MorganInvariants::Features,
        },
        ..MorganFingerprintParams::default()
    };

    let mut additional_output = FingerprintAdditionalOutput::default();
    additional_output.allocate_atom_counts();
    additional_output.allocate_atom_to_bits();
    additional_output.allocate_bit_info_map();
    additional_output.allocate_bit_paths();
    additional_output.allocate_atoms_per_bit();

    match output_kind {
        MorganOutputKind::DenseBits => {
            let fingerprint = molecule
                .morgan_fingerprint_with_params(&params, Some(&mut additional_output))
                .map_err(|error| error.to_string())?;
            Ok(Outcome::MorganDenseBits(MorganDenseBitsOutput {
                length: fingerprint.n_bits(),
                on_bits: fingerprint.on_bits(),
                additional_output: observed_morgan_additional_output(&additional_output),
            }))
        }
        MorganOutputKind::SparseBits => {
            let fingerprint = molecule
                .morgan_sparse_fingerprint_with_params(&params, Some(&mut additional_output))
                .map_err(|error| error.to_string())?;
            Ok(Outcome::MorganSparseBits(MorganSparseBitsOutput {
                length: fingerprint.n_bits(),
                on_bits: fingerprint.on_bits(),
                additional_output: observed_morgan_additional_output(&additional_output),
            }))
        }
        MorganOutputKind::HashedCounts => {
            let fingerprint = molecule
                .morgan_count_fingerprint_with_params(&params, Some(&mut additional_output))
                .map_err(|error| error.to_string())?;
            Ok(Outcome::MorganHashedCounts(MorganHashedCountsOutput {
                length: fingerprint.length(),
                entries: fingerprint
                    .nonzero_elements()
                    .iter()
                    .map(|(&index, &count)| (index, count))
                    .collect(),
                additional_output: observed_morgan_additional_output(&additional_output),
            }))
        }
        MorganOutputKind::SparseCounts => {
            let fingerprint = molecule
                .morgan_sparse_count_fingerprint_with_params(&params, Some(&mut additional_output))
                .map_err(|error| error.to_string())?;
            Ok(Outcome::MorganSparseCounts(MorganSparseCountsOutput {
                length: fingerprint.length(),
                entries: fingerprint
                    .nonzero_elements()
                    .iter()
                    .map(|(&index, &count)| (index, count))
                    .collect(),
                additional_output: observed_morgan_additional_output(&additional_output),
            }))
        }
    }
}

fn observed_morgan_additional_output(
    output: &FingerprintAdditionalOutput,
) -> MorganAdditionalOutput {
    MorganAdditionalOutput {
        atom_counts: output.atom_counts().map(<[u32]>::to_vec),
        atom_to_bits: output.atom_to_bits().map(<[Vec<u64>]>::to_vec),
        bit_info_map: output.bit_info_map().map(|map| {
            map.iter()
                .map(|(&bit, centers)| (bit, centers.clone()))
                .collect()
        }),
        bit_paths: output.bit_paths().map(|map| {
            map.iter()
                .map(|(&bit, paths)| (bit, paths.clone()))
                .collect()
        }),
        atoms_per_bit: output.atoms_per_bit().map(|map| {
            map.iter()
                .map(|(&bit, atoms)| (bit, atoms.clone()))
                .collect()
        }),
    }
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

/// Preserve the four original substitutions and the approved canonical identity.
pub fn svg_matches(expected: &Outcome, actual: &Outcome) -> bool {
    let (Outcome::Text(expected), Outcome::Text(actual)) = (expected, actual) else {
        return false;
    };
    // ROOT-SVG-CANONICAL-METADATA-20261005: new CK handling is limited to
    // namespace declarations and QNames bound to the exact canonical URI.
    // Text, comments, CDATA and ordinary attribute values remain byte-exact.
    fn canonical_identity(svg: &str) -> String {
        let bytes = svg.as_bytes();
        let mut replacements = Vec::new();
        let mut scopes = vec![false];
        let mut cursor = 0;
        while cursor < bytes.len() {
            if bytes[cursor] != b'<' {
                cursor += 1;
                continue;
            }
            let suffix = &svg[cursor..];
            let opaque_end = if suffix.starts_with("<!--") {
                Some("-->")
            } else if suffix.starts_with("<![CDATA[") {
                Some("]]>")
            } else if suffix.starts_with("<?") {
                Some("?>")
            } else {
                None
            };
            if let Some(end) = opaque_end {
                cursor += suffix.find(end).map_or(suffix.len(), |n| n + end.len());
                continue;
            }
            let closing = bytes.get(cursor + 1) == Some(&b'/');
            let name_start = cursor + 1 + usize::from(closing);
            if !bytes
                .get(name_start)
                .is_some_and(|b| b.is_ascii_alphabetic() || *b == b'_')
            {
                cursor += 1;
                continue;
            }
            let mut end = name_start;
            let mut quote = None;
            while end < bytes.len() {
                match (quote, bytes[end]) {
                    (Some(q), b) if q == b => quote = None,
                    (None, b'\'' | b'"') => quote = Some(bytes[end]),
                    (None, b'>') => break,
                    _ => {}
                }
                end += 1;
            }
            if end == bytes.len() {
                break;
            }
            let name_end = (name_start..end)
                .find(|&n| bytes[n].is_ascii_whitespace() || bytes[n] == b'/')
                .unwrap_or(end);
            let mut names = vec![(name_start, name_end)];
            let mut active = *scopes.last().unwrap();
            let mut pos = name_end;
            while !closing && pos < end {
                while pos < end && bytes[pos].is_ascii_whitespace() {
                    pos += 1;
                }
                if pos == end || bytes[pos] == b'/' {
                    break;
                }
                let attr_start = pos;
                while pos < end
                    && !bytes[pos].is_ascii_whitespace()
                    && !matches!(bytes[pos], b'=' | b'/')
                {
                    pos += 1;
                }
                let attr_end = pos;
                if attr_start == attr_end {
                    break;
                }
                names.push((attr_start, attr_end));
                while pos < end && bytes[pos].is_ascii_whitespace() {
                    pos += 1;
                }
                if bytes.get(pos) != Some(&b'=') {
                    break;
                }
                pos += 1;
                while pos < end && bytes[pos].is_ascii_whitespace() {
                    pos += 1;
                }
                let Some(&delimiter @ (b'\'' | b'"')) = bytes.get(pos) else {
                    break;
                };
                pos += 1;
                while pos < end && bytes[pos] != delimiter {
                    pos += 1;
                }
                if pos < end {
                    pos += 1;
                }
                if &svg[attr_start..attr_end] == "xmlns:ck" {
                    active = &svg[attr_start..pos] == "xmlns:ck='https://kit.cosmol.org/'";
                    if active {
                        replacements.push((attr_start, pos, "xmlns:tool='__tool_namespace__'"));
                    }
                }
            }
            if active {
                for (start, stop) in names {
                    if svg[start..stop].starts_with("ck:") {
                        replacements.push((start, start + 3, "tool:"));
                    }
                }
            }
            if closing {
                if scopes.len() > 1 {
                    scopes.pop();
                }
            } else if bytes[end - 1] != b'/' {
                scopes.push(active);
            }
            cursor = end + 1;
        }
        replacements.sort_unstable_by_key(|&(start, _, _)| start);
        let mut result = String::with_capacity(svg.len());
        let mut copied = 0;
        for (start, end, replacement) in replacements {
            result.push_str(&svg[copied..start]);
            result.push_str(replacement);
            copied = end;
        }
        result.push_str(&svg[copied..]);
        result
    }
    fn normalize(svg: &str) -> String {
        let original = svg
            .replace(
                "xmlns:rdkit='http://www.rdkit.org/xml'",
                "xmlns:tool='__tool_namespace__'",
            )
            .replace(
                "xmlns:cosmolkit='https://www.cosmol.org'",
                "xmlns:tool='__tool_namespace__'",
            )
            .replace("rdkit:", "tool:")
            .replace("cosmolkit:", "tool:");
        canonical_identity(&original)
    }
    normalize(expected) == normalize(actual)
}
