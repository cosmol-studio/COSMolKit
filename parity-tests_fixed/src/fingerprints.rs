//! One seeded parameter combination per molecule and task, shared with RDKit.
use crate::registry::{Input, Record, SmilesCase, Value};
use cosmolkit::{
    Fingerprint, LayeredFingerprintLayers, LayeredFingerprintParams, Molecule,
    MorganFingerprintParams, MorganParams, PatternFingerprintParams, SparseCountFingerprint,
    SparseCountFingerprint32, TopologicalFingerprintParams,
};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

pub const SEED: u64 = 0x434b_4650_2026_1007;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Kind {
    Maccs,
    Topological,
    Layered,
    Pattern,
    FuzzyAnd,
    FuzzyOr,
}
impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Self::Maccs => "fingerprint_maccs",
            Self::Topological => "fingerprint_topological",
            Self::Layered => "fingerprint_layered",
            Self::Pattern => "fingerprint_pattern",
            Self::FuzzyAnd => "fuzzy_and",
            Self::FuzzyOr => "fuzzy_or",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Params {
    Maccs,
    Topological {
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        bits_per_feature: u32,
        use_hs: bool,
        density_milli: u32,
        min_size: u32,
        branched: bool,
        bond_order: bool,
        custom_invariants: bool,
        roots: Roots,
    },
    Layered {
        layers: u32,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        branched: bool,
        roots: Roots,
        counts: Counts,
        mask: Mask,
    },
    Pattern {
        fp_size: u32,
        tautomeric: bool,
    },
    Fuzzy {
        union: bool,
        wide: bool,
        signed: bool,
        fp_size: u32,
        radius: u32,
    },
}
impl Params {
    fn kind(&self) -> Kind {
        match self {
            Self::Maccs => Kind::Maccs,
            Self::Topological { .. } => Kind::Topological,
            Self::Layered { .. } => Kind::Layered,
            Self::Pattern { .. } => Kind::Pattern,
            Self::Fuzzy { union: false, .. } => Kind::FuzzyAnd,
            Self::Fuzzy { union: true, .. } => Kind::FuzzyOr,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Roots {
    All,
    First,
    Terminals,
}
impl Roots {
    fn resolve(self, n: usize) -> Option<Vec<u32>> {
        match self {
            Self::All => None,
            Self::First => Some((0..n.min(1) as u32).collect()),
            Self::Terminals => Some(match n {
                0 => vec![],
                1 => vec![0],
                _ => vec![0, n as u32 - 1],
            }),
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Counts {
    Absent,
    Zero,
    Seeded,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Mask {
    Absent,
    Even,
    ModThree,
    Empty,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FingerprintInput {
    pub case: SmilesCase,
    pub right: Option<SmilesCase>,
    pub seed: u64,
    pub params: Params,
}
impl FingerprintInput {
    pub fn task_name(&self) -> &'static str {
        self.params.kind().name()
    }
}

// Local, specified PRNG: reproducible across platforms and dependency updates.
struct Random(u64);
impl Random {
    fn pick<T: Copy>(&mut self, values: &[T]) -> T {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        values[((z ^ (z >> 31)) % values.len() as u64) as usize]
    }
    fn boolean(&mut self) -> bool {
        self.pick(&[false, true])
    }
}

pub fn inputs(cases: &[SmilesCase], kind: Kind) -> Vec<Input> {
    cases
        .iter()
        .enumerate()
        .map(|(index, case)| {
            let mut hash = Sha256::new();
            for part in [
                SEED.to_le_bytes().as_slice(),
                kind.name().as_bytes(),
                case.id.as_bytes(),
                case.smiles.as_bytes(),
            ] {
                hash.update((part.len() as u64).to_le_bytes());
                hash.update(part);
            }
            let seed = u64::from_le_bytes(hash.finalize()[..8].try_into().unwrap());
            let mut rng = Random(seed);
            let params = match kind {
                Kind::Maccs => Params::Maccs,
                Kind::Topological => {
                    let min_path = rng.pick(&[1, 2]);
                    Params::Topological {
                        min_path,
                        max_path: rng.pick(&[4, 7]),
                        fp_size: rng.pick(&[128, 256, 1024, 2048, 4096]),
                        bits_per_feature: rng.pick(&[1, 2, 4]),
                        use_hs: rng.boolean(),
                        density_milli: rng.pick(&[0, 200, 350]),
                        min_size: rng.pick(&[64, 128]),
                        branched: rng.boolean(),
                        bond_order: rng.boolean(),
                        custom_invariants: rng.boolean(),
                        roots: rng.pick(&[Roots::All, Roots::First, Roots::Terminals]),
                    }
                }
                Kind::Layered => Params::Layered {
                    layers: rng.pick(&[u32::MAX, 0x3f, 1, 2, 4, 8, 16, 32, 7, 0xffffffc0, 0]),
                    min_path: rng.pick(&[1, 2]),
                    max_path: rng.pick(&[4, 7]),
                    fp_size: rng.pick(&[64, 256, 257, 2048, 4096]),
                    branched: rng.boolean(),
                    roots: rng.pick(&[Roots::All, Roots::First, Roots::Terminals]),
                    counts: rng.pick(&[Counts::Absent, Counts::Zero, Counts::Seeded]),
                    mask: rng.pick(&[Mask::Absent, Mask::Even, Mask::ModThree, Mask::Empty]),
                },
                Kind::Pattern => Params::Pattern {
                    fp_size: rng.pick(&[64, 256, 257, 2048, 4096]),
                    tautomeric: rng.boolean(),
                },
                Kind::FuzzyAnd | Kind::FuzzyOr => Params::Fuzzy {
                    union: kind == Kind::FuzzyOr,
                    wide: rng.boolean(),
                    signed: rng.boolean(),
                    fp_size: rng.pick(&[64, 257, 2048]),
                    radius: rng.pick(&[0, 1, 2, 3]),
                },
            };
            Input::Fingerprint(FingerprintInput {
                case: case.clone(),
                seed,
                params,
                right: matches!(kind, Kind::FuzzyAnd | Kind::FuzzyOr)
                    .then(|| cases[(index + 1) % cases.len()].clone()),
            })
        })
        .collect()
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Bits {
    pub length: u32,
    pub on_bits: Vec<u32>,
}
impl From<Fingerprint> for Bits {
    fn from(fp: Fingerprint) -> Self {
        Self {
            length: fp.n_bits(),
            on_bits: fp.on_bits(),
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CountVector {
    pub length: u64,
    pub entries: Vec<(u64, i32)>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Observation {
    /// Actual isolated reference process failure; never a chemical value or match.
    ReferenceProcessFailure {
        exit_code: i64,
        process_id: u32,
        stderr: String,
    },
    Maccs {
        raw: Bits,
        public: Bits,
    },
    Bits {
        fingerprint: Bits,
        atom_counts: Option<Vec<u32>>,
    },
    Counts {
        left: CountVector,
        right: CountVector,
        result: CountVector,
    },
}

pub fn validate_output(row: &FingerprintInput, output: &Observation) -> Result<(), String> {
    let valid_bits = |fp: &Bits| {
        fp.length > 0
            && fp.on_bits.iter().all(|&i| i < fp.length)
            && fp.on_bits.windows(2).all(|pair| pair[0] < pair[1])
    };
    let valid_counts = |fp: &CountVector, length| {
        fp.length == length
            && fp.entries.iter().all(|&(i, _)| i < length)
            && fp.entries.windows(2).all(|pair| pair[0].0 < pair[1].0)
    };
    let valid = match (&row.params, output) {
        (
            Params::Layered {
                branched: false,
                roots: Roots::All,
                ..
            },
            Observation::ReferenceProcessFailure {
                exit_code,
                process_id,
                ..
            },
        ) => (*exit_code < 0 || *exit_code >= 0x8000_0000) && *process_id != 0,
        (Params::Maccs, Observation::Maccs { raw, public }) => {
            raw.length == 167
                && public.length == 166
                && !raw.on_bits.contains(&0)
                && valid_bits(raw)
                && valid_bits(public)
        }
        (
            Params::Topological { fp_size, .. },
            Observation::Bits {
                fingerprint,
                atom_counts,
            },
        ) => valid_bits(fingerprint) && fingerprint.length <= *fp_size && atom_counts.is_none(),
        (
            Params::Layered {
                fp_size, counts, ..
            },
            Observation::Bits {
                fingerprint,
                atom_counts,
            },
        ) => {
            valid_bits(fingerprint)
                && fingerprint.length == *fp_size
                && atom_counts.is_none() == (*counts == Counts::Absent)
        }
        (
            Params::Pattern { fp_size, .. },
            Observation::Bits {
                fingerprint,
                atom_counts,
            },
        ) => valid_bits(fingerprint) && fingerprint.length == *fp_size && atom_counts.is_none(),
        (
            Params::Fuzzy { fp_size, wide, .. },
            Observation::Counts {
                left,
                right,
                result,
            },
        ) => {
            let length = u64::from(*fp_size) + if *wide { 1u64 << 32 } else { 0 };
            row.right.is_some()
                && [left, right, result]
                    .iter()
                    .all(|fp| valid_counts(fp, length))
        }
        _ => false,
    };
    if valid {
        Ok(())
    } else {
        Err("fingerprint reference kind/shape mismatch".into())
    }
}

pub fn matches(expected: &Observation, actual: &Observation) -> bool {
    !matches!(expected, Observation::ReferenceProcessFailure { .. }) && expected == actual
}

pub fn run(row: &FingerprintInput) -> Result<Record, String> {
    let mol = Molecule::from_smiles(&row.case.smiles).map_err(|e| e.to_string())?;
    let n = mol.num_atoms();
    let observed = match row.params {
        Params::Maccs => Observation::Maccs {
            raw: mol
                .fingerprint_maccs_raw()
                .map_err(|e| e.to_string())?
                .into(),
            public: mol.fingerprint_maccs().map_err(|e| e.to_string())?.into(),
        },
        Params::Topological {
            min_path,
            max_path,
            fp_size,
            bits_per_feature,
            use_hs,
            density_milli,
            min_size,
            branched,
            bond_order,
            custom_invariants,
            roots,
        } => {
            let params = TopologicalFingerprintParams {
                min_path,
                max_path,
                fp_size,
                num_bits_per_feature: bits_per_feature,
                use_hs,
                target_density: f64::from(density_milli) / 1000.0,
                min_size,
                branched_paths: branched,
                use_bond_order: bond_order,
                atom_invariants: custom_invariants.then(|| (1..=n as u32).collect()),
                from_atoms: roots.resolve(n),
                ignore_atoms: None,
            };
            Observation::Bits {
                fingerprint: mol
                    .fingerprint_topological_with_params(&params)
                    .map_err(|e| e.to_string())?
                    .into(),
                atom_counts: None,
            }
        }
        Params::Layered {
            layers,
            min_path,
            max_path,
            fp_size,
            branched,
            roots,
            counts,
            mask,
        } => {
            let atom_counts = match counts {
                Counts::Absent => None,
                Counts::Zero => Some(vec![0; n]),
                Counts::Seeded => Some((0..n as u32).map(|i| i + 10).collect()),
            };
            let set_only_bits = match mask {
                Mask::Absent => None,
                _ => Some(
                    Fingerprint::from_on_bits(
                        fp_size,
                        (0..fp_size).filter(|i| match mask {
                            Mask::Even => i % 2 == 0,
                            Mask::ModThree => i % 3 == 0,
                            _ => false,
                        }),
                    )
                    .map_err(|e| e.to_string())?,
                ),
            };
            let params = LayeredFingerprintParams {
                layers: LayeredFingerprintLayers::from_bits_retain(layers),
                min_path,
                max_path,
                fp_size,
                branched_paths: branched,
                from_atoms: roots.resolve(n),
                atom_counts,
                set_only_bits,
            };
            let before = params.clone();
            let result = mol
                .fingerprint_layered_with_output_with_params(&params)
                .map_err(|e| e.to_string())?;
            if params != before {
                return Err("Layered mutated input parameters".into());
            }
            Observation::Bits {
                fingerprint: result.fingerprint.into(),
                atom_counts: result.atom_counts,
            }
        }
        Params::Pattern {
            fp_size,
            tautomeric,
        } => Observation::Bits {
            fingerprint: mol
                .fingerprint_pattern_with_params(&PatternFingerprintParams {
                    n_bits: fp_size as usize,
                    tautomeric,
                })
                .map_err(|e| e.to_string())?
                .into(),
            atom_counts: None,
        },
        Params::Fuzzy {
            union,
            wide,
            signed,
            fp_size,
            radius,
        } => {
            let right = row.right.as_ref().ok_or("missing fuzzy right operand")?;
            let other = Molecule::from_smiles(&right.smiles).map_err(|e| e.to_string())?;
            let params = MorganFingerprintParams {
                generator: MorganParams {
                    fp_size,
                    radius,
                    ..Default::default()
                },
                ..Default::default()
            };
            let left = mol
                .fingerprint_morgan_count_with_params(&params, None)
                .map_err(|e| e.to_string())?;
            let right = other
                .fingerprint_morgan_count_with_params(&params, None)
                .map_err(|e| e.to_string())?;
            fuzzy(left, right, union, wide, signed)?
        }
    };
    Ok(Record {
        input: Input::Fingerprint(row.clone()),
        output: Value::Fingerprint(observed),
    })
}

fn fuzzy(
    left: SparseCountFingerprint32,
    right: SparseCountFingerprint32,
    union: bool,
    wide: bool,
    signed: bool,
) -> Result<Observation, String> {
    // Exercise both source index widths, including keys above u32::MAX.
    let offset = if wide { 1u64 << 32 } else { 0 };
    let length = u64::from(left.length()) + offset;
    let entries = |fp: &SparseCountFingerprint32| {
        fp.nonzero_elements()
            .iter()
            .map(|(&i, &v)| {
                (
                    u64::from(i) + offset,
                    if signed && i % 2 == 1 { -v } else { v },
                )
            })
            .collect()
    };
    let l = CountVector {
        length,
        entries: entries(&left),
    };
    let r = CountVector {
        length,
        entries: entries(&right),
    };
    macro_rules! evaluate {
        ($ty:ty, $width:ty) => {{
            let make = |v: &CountVector| -> Result<$ty, String> {
                let mut fp = <$ty>::new(v.length as $width);
                for &(i, value) in &v.entries {
                    fp.set_value(i as $width, value)
                        .map_err(|e| e.to_string())?;
                }
                Ok(fp)
            };
            let a = make(&l)?;
            let b = make(&r)?;
            let before = (a.clone(), b.clone());
            let result = if union {
                a.fuzzy_or(&b)
            } else {
                a.fuzzy_and(&b)
            }
            .map_err(|e| e.to_string())?;
            if (a, b) != before {
                return Err("fuzzy operation mutated an operand".into());
            }
            CountVector {
                length: result.length() as u64,
                entries: result
                    .nonzero_elements()
                    .iter()
                    .map(|(&i, &v)| (i as u64, v))
                    .collect(),
            }
        }};
    }
    let result = if wide {
        evaluate!(SparseCountFingerprint, u64)
    } else {
        evaluate!(SparseCountFingerprint32, u32)
    };
    Ok(Observation::Counts {
        left: l,
        right: r,
        result,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reference_process_failure_requires_the_isolated_branch_and_real_exit() {
        let cases = vec![SmilesCase {
            id: "fixed-CCO".into(),
            smiles: "CCO".into(),
        }];
        let Input::Fingerprint(mut row) = inputs(&cases, Kind::Layered).remove(0) else {
            unreachable!()
        };
        row.params = Params::Layered {
            layers: 63,
            min_path: 1,
            max_path: 7,
            fp_size: 2048,
            branched: false,
            roots: Roots::All,
            counts: Counts::Seeded,
            mask: Mask::Absent,
        };
        for exit_code in [-11, -9, 0xc000_0005] {
            let failure = Observation::ReferenceProcessFailure {
                exit_code,
                process_id: 123,
                stderr: "native diagnostic".into(),
            };
            validate_output(&row, &failure).unwrap();
            assert!(!matches(&failure, &failure));
        }
        for (exit_code, process_id) in [(0, 123), (1, 123), (-11, 0)] {
            assert!(
                validate_output(
                    &row,
                    &Observation::ReferenceProcessFailure {
                        exit_code,
                        process_id,
                        stderr: String::new(),
                    }
                )
                .is_err()
            );
        }
        let failure = Observation::ReferenceProcessFailure {
            exit_code: -11,
            process_id: 123,
            stderr: String::new(),
        };
        if let Params::Layered { roots, .. } = &mut row.params {
            *roots = Roots::First;
        }
        assert!(validate_output(&row, &failure).is_err());
        row.params = Params::Maccs;
        assert!(validate_output(&row, &failure).is_err());
    }

    #[test]
    fn fingerprint_seed_is_stable_and_one_profile_per_case() {
        let cases: Vec<_> = (0..1000)
            .map(|i| SmilesCase {
                id: i.to_string(),
                smiles: "CCO".into(),
            })
            .collect();
        for kind in [
            Kind::Maccs,
            Kind::Topological,
            Kind::Layered,
            Kind::Pattern,
            Kind::FuzzyAnd,
            Kind::FuzzyOr,
        ] {
            let rows = inputs(&cases, kind);
            assert_eq!(rows, inputs(&cases, kind));
            assert_eq!(rows.len(), cases.len());
            assert!(rows.iter().all(|row| row.task_name() == kind.name()));
            let distinct: std::collections::BTreeSet<_> = rows
                .iter()
                .map(|row| match row {
                    Input::Fingerprint(row) => serde_json::to_string(&row.params).unwrap(),
                    _ => unreachable!(),
                })
                .collect();
            assert!(distinct.len() >= if kind == Kind::Maccs { 1 } else { 10 });
        }
    }
}
