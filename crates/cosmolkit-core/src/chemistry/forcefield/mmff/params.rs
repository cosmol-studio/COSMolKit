//! Source-backed RDKit MMFF parameter structures.

use std::borrow::Cow;
use std::collections::BTreeMap;
use std::fmt;
use std::num::{ParseFloatError, ParseIntError};

include!(concat!(env!("OUT_DIR"), "/mmff_defaults_generated.rs"));

/// MMFF atom type equivalence levels.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MmffDef {
    // RDKit✔️✔️: std::uint8_t eqLevel[4];
    pub eq_level: [u8; 4],
}

/// MMFF atom-type properties.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MmffProp {
    // RDKit✔️✔️: std::uint8_t atno;
    pub atno: u8,
    // RDKit✔️✔️: std::uint8_t crd;
    pub crd: u8,
    // RDKit✔️✔️: std::uint8_t val;
    pub val: u8,
    // RDKit✔️✔️: std::uint8_t pilp;
    pub pilp: u8,
    // RDKit✔️✔️: std::uint8_t mltb;
    pub mltb: u8,
    // RDKit✔️✔️: std::uint8_t arom;
    pub arom: u8,
    // RDKit✔️✔️: std::uint8_t linh;
    pub linh: u8,
    // RDKit✔️✔️: std::uint8_t sbmb;
    pub sbmb: u8,
}

/// MMFF partial bond charge increments.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffPbci {
    // RDKit✔️✔️: double pbci;
    pub pbci: f64,
    // RDKit✔️✔️: double fcadj;
    pub fcadj: f64,
}

/// MMFF bond-charge-increment parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffChg {
    // RDKit✔️✔️: double bci;
    pub bci: f64,
}

/// MMFF bond-stretching parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffBond {
    // RDKit✔️✔️: double kb;
    pub kb: f64,
    // RDKit✔️✔️: double r0;
    pub r0: f64,
}

/// Parameters for Herschbach-Laurie's version of Badger's rule.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffHerschbachLaurie {
    // RDKit❗✔️: class RDKIT_FORCEFIELD_EXPORT MMFFHerschbachLaurie {
    // RDKit❗✔️:  public:
    // RDKit❗✔️:   double a_ij;
    pub a_ij: f64,
    // RDKit❗✔️:   double d_ij;
    pub d_ij: f64,
    // RDKit❗✔️:   double dp_ij;
    pub dp_ij: f64,
    // RDKit❗✔️: };
}

/// Covalent-radius and Pauling-electronegativity values for the MMFF bond rule.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffCovRadPauEle {
    // RDKit❗✔️: class RDKIT_FORCEFIELD_EXPORT MMFFCovRadPauEle {
    // RDKit❗✔️:  public:
    // RDKit❗✔️:   double r0;
    pub r0: f64,
    // RDKit❗✔️:   double chi;
    pub chi: f64,
    // RDKit❗✔️: };
}

/// MMFF angle-bending parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffAngle {
    // RDKit✔️✔️: double ka;
    pub ka: f64,
    // RDKit✔️✔️: double theta0;
    pub theta0: f64,
}

/// MMFF stretch-bend parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffStbn {
    // RDKit✔️✔️: double kbaIJK;
    pub kba_ijk: f64,
    // RDKit✔️✔️: double kbaKJI;
    pub kba_kji: f64,
}

/// MMFF out-of-plane bending parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffOop {
    // RDKit✔️✔️: double koop;
    pub koop: f64,
}

/// MMFF torsion parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffTor {
    // RDKit✔️✔️: double V1;
    pub v1: f64,
    // RDKit✔️✔️: double V2;
    pub v2: f64,
    // RDKit✔️✔️: double V3;
    pub v3: f64,
}

/// MMFF non-bonded Van der Waals parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffVdw {
    // RDKit✔️✔️: double alpha_i;
    pub alpha_i: f64,
    // RDKit✔️✔️: double N_i;
    pub n_i: f64,
    // RDKit✔️✔️: double A_i;
    pub a_i: f64,
    // RDKit✔️✔️: double G_i;
    pub g_i: f64,
    // RDKit✔️✔️: double R_star;
    pub r_star: f64,
    // RDKit✔️✔️: std::uint8_t DA;
    pub da: u8,
}

/// MMFF combined Van der Waals pair constants.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffVdwRijstarEps {
    // RDKit✔️✔️: double R_ij_starUnscaled;
    pub r_ij_star_unscaled: f64,
    // RDKit✔️✔️: double epsilonUnscaled;
    pub epsilon_unscaled: f64,
    // RDKit✔️✔️: double R_ij_star;
    pub r_ij_star: f64,
    // RDKit✔️✔️: double epsilon;
    pub epsilon: f64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MmffParamError {
    MissingDefaultData {
        symbol: &'static str,
    },
    MalformedDefaultData {
        symbol: &'static str,
        line: String,
    },
    MalformedMmffDefLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffPropLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffPbciLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffChgLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffBondLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffAngleLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffStbnLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffDfsbLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffOopLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffTorLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffVdwConstantsLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffVdwLine {
        line_number: usize,
        column_count: usize,
    },
    MalformedMmffVdwDa {
        line_number: usize,
    },
    ParseInt {
        line_number: usize,
        column_name: &'static str,
        value: String,
    },
    ParseFloat {
        line_number: usize,
        column_name: &'static str,
        value: String,
    },
}

impl fmt::Display for MmffParamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingDefaultData { symbol } => {
                write!(f, "missing RDKit MMFF default data symbol {symbol}")
            }
            Self::MalformedDefaultData { symbol, line } => {
                write!(
                    f,
                    "malformed RDKit MMFF default data string line for {symbol}: {line}"
                )
            }
            Self::MalformedMmffDefLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF definition line {line_number}: expected at least 6 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffPropLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF property line {line_number}: expected at least 9 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffPbciLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF PBCI line {line_number}: expected at least 4 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffChgLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF charge line {line_number}: expected at least 4 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffBondLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF bond line {line_number}: expected at least 5 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffAngleLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF angle line {line_number}: expected at least 6 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffStbnLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF stretch-bend line {line_number}: expected at least 6 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffDfsbLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF default stretch-bend line {line_number}: expected at least 5 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffOopLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF out-of-plane line {line_number}: expected at least 5 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffTorLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF torsion line {line_number}: expected at least 8 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffVdwConstantsLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF Van der Waals constants line {line_number}: expected at least 5 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffVdwLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed MMFF Van der Waals line {line_number}: expected at least 6 tab-separated columns, found {column_count}"
            ),
            Self::MalformedMmffVdwDa { line_number } => write!(
                f,
                "malformed MMFF Van der Waals DA token at line {line_number}: expected a non-empty string"
            ),
            Self::ParseInt {
                line_number,
                column_name,
                value,
            } => write!(
                f,
                "invalid MMFF definition integer at line {line_number}, column {column_name}: {value}"
            ),
            Self::ParseFloat {
                line_number,
                column_name,
                value,
            } => write!(
                f,
                "invalid MMFF definition float at line {line_number}, column {column_name}: {value}"
            ),
        }
    }
}

impl std::error::Error for MmffParamError {}

/// MMFF atom type equivalence definition collection.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MmffDefCollection {
    source_mmff_def: Cow<'static, str>,
    params: MmffDefStorage,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum MmffDefStorage {
    Static(&'static [(u8, MmffDef)]),
    Owned(Vec<MmffDef>),
}

impl MmffDefCollection {
    pub fn new(mmff_def: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFDefCollection::MMFFDefCollection(std::string mmffDef) {
        // RDKit✔️✔️:   if (mmffDef.empty()) {
        // RDKit✔️✔️:     mmffDef = defaultMMFFDef;
        // RDKit✔️✔️:   }
        let source_mmff_def = if mmff_def.is_empty() {
            Cow::Borrowed(default_mmff_def()?)
        } else {
            Cow::Owned(mmff_def.to_owned())
        };
        let params = if mmff_def.is_empty() {
            MmffDefStorage::Static(DEFAULT_MMFF_DEF_ROWS)
        } else {
            MmffDefStorage::Owned(parse_mmff_def(source_mmff_def.as_ref())?)
        };
        Ok(Self {
            source_mmff_def,
            params,
        })
    }

    pub fn source_mmff_def(&self) -> &str {
        self.source_mmff_def.as_ref()
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffDefStorage::Static(params) => params.len(),
            MmffDefStorage::Owned(params) => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffDefStorage::Static(params) => params.is_empty(),
            MmffDefStorage::Owned(params) => params.is_empty(),
        }
    }

    pub fn get(&self, atom_type: u32) -> Option<&MmffDef> {
        // RDKit✔️✔️:   const MMFFDef *operator()(const unsigned int atomType) const {
        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res = d_params.find(atomType);
        // RDKit❌❌:     return ((res != d_params.end()) ? &((*res).second) : NULL);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     return ((atomType && (atomType <= d_params.size()))
        // RDKit✔️✔️:                 ? &d_params[atomType - 1]
        // RDKit✔️✔️:                 : nullptr);
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:   }
        match &self.params {
            MmffDefStorage::Static(params) => {
                if atom_type == 0 {
                    None
                } else {
                    params
                        .binary_search_by_key(&(atom_type as u8), |(key, _)| *key)
                        .ok()
                        .map(|idx| &params[idx].1)
                }
            }
            MmffDefStorage::Owned(params) => {
                if atom_type != 0 && (atom_type as usize) <= params.len() {
                    params.get(atom_type as usize - 1)
                } else {
                    None
                }
            }
        }
    }

    fn torsion_equivalence_level(&self, atom_type: u32, level: usize) -> Option<u8> {
        let definition = self.get(atom_type)?;
        if let Some(&equivalent_type) = definition.eq_level.get(level) {
            return Some(equivalent_type);
        }

        // RDKit✔️✔️:       unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iWildCard];
        // RDKit✔️✔️:       unsigned int canLAtomType = (*mmffDef)(lAtomType)->eqLevel[lWildCard];
        //
        // The source fallback loop has one additional `iter == 4` pass for a
        // type-5 torsion with a secondary type, although `eqLevel` contains
        // four entries. The pinned RDKit build observes zero for both terminal
        // types on that pass, selecting the source table's `0-...-0` wildcard
        // row when present. Model that observable source behavior explicitly;
        // indexing beyond the Rust array would abort, while reproducing the
        // C++ out-of-bounds access would be undefined behavior.
        (level == definition.eq_level.len()).then_some(0)
    }
}

/// MMFF atom-type property collection.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MmffPropCollection {
    source_mmff_prop: Cow<'static, str>,
    params: MmffPropStorage,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum MmffPropStorage {
    Static(&'static [(u8, MmffProp)]),
    Owned {
        params: Vec<MmffProp>,
        i_atom_type: Vec<u8>,
    },
}

impl MmffPropCollection {
    pub fn new(mmff_prop: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFPropCollection::MMFFPropCollection(std::string mmffProp) {
        // RDKit✔️✔️:   if (mmffProp.empty()) {
        // RDKit✔️✔️:     mmffProp = defaultMMFFProp;
        // RDKit✔️✔️:   }
        let source_mmff_prop = if mmff_prop.is_empty() {
            Cow::Borrowed(default_mmff_prop()?)
        } else {
            Cow::Owned(mmff_prop.to_owned())
        };
        let params = if mmff_prop.is_empty() {
            MmffPropStorage::Static(DEFAULT_MMFF_PROP_ROWS)
        } else {
            let (params, i_atom_type) = parse_mmff_prop(source_mmff_prop.as_ref())?;
            MmffPropStorage::Owned {
                params,
                i_atom_type,
            }
        };
        Ok(Self {
            source_mmff_prop,
            params,
        })
    }

    pub fn source_mmff_prop(&self) -> &str {
        self.source_mmff_prop.as_ref()
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffPropStorage::Static(params) => params.len(),
            MmffPropStorage::Owned { params, .. } => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffPropStorage::Static(params) => params.is_empty(),
            MmffPropStorage::Owned { params, .. } => params.is_empty(),
        }
    }

    pub fn get(&self, atom_type: u32) -> Option<&MmffProp> {
        // RDKit✔️✔️:   const MMFFProp *operator()(const unsigned int atomType) const {
        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res = d_params.find(atomType);
        // RDKit❌❌:     return ((res != d_params.end()) ? &((*res).second) : NULL);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto bounds =
        // RDKit✔️✔️:         std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), atomType);
        // RDKit✔️✔️:     return ((bounds.first != bounds.second)
        // RDKit✔️✔️:                 ? &d_params[bounds.first - d_iAtomType.begin()]
        // RDKit✔️✔️:                 : nullptr);
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:   }
        match &self.params {
            MmffPropStorage::Static(params) => params
                .binary_search_by_key(&(atom_type as u8), |(key, _)| *key)
                .ok()
                .map(|idx| &params[idx].1),
            MmffPropStorage::Owned {
                params,
                i_atom_type,
            } => {
                let lower = i_atom_type.partition_point(|&probe| u32::from(probe) < atom_type);
                let upper = i_atom_type.partition_point(|&probe| u32::from(probe) <= atom_type);
                if lower != upper {
                    params.get(lower)
                } else {
                    None
                }
            }
        }
    }
}

/// MMFF partial bond charge increment collection.
#[derive(Clone, Debug, PartialEq)]
pub struct MmffPbciCollection {
    source_mmff_pbci: Cow<'static, str>,
    params: MmffPbciStorage,
}

#[derive(Clone, Debug, PartialEq)]
enum MmffPbciStorage {
    Static(&'static [(u8, MmffPbci)]),
    Owned(Vec<MmffPbci>),
}

impl MmffPbciCollection {
    pub fn new(mmff_pbci: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFPBCICollection::MMFFPBCICollection(std::string mmffPBCI) {
        // RDKit✔️✔️:   if (mmffPBCI.empty()) {
        // RDKit✔️✔️:     mmffPBCI = defaultMMFFPBCI;
        // RDKit✔️✔️:   }
        let source_mmff_pbci = if mmff_pbci.is_empty() {
            Cow::Borrowed(default_mmff_pbci()?)
        } else {
            Cow::Owned(mmff_pbci.to_owned())
        };
        let params = if mmff_pbci.is_empty() {
            MmffPbciStorage::Static(DEFAULT_MMFF_PBCI_ROWS)
        } else {
            MmffPbciStorage::Owned(parse_mmff_pbci(source_mmff_pbci.as_ref())?)
        };
        Ok(Self {
            source_mmff_pbci,
            params,
        })
    }

    pub fn source_mmff_pbci(&self) -> &str {
        self.source_mmff_pbci.as_ref()
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffPbciStorage::Static(params) => params.len(),
            MmffPbciStorage::Owned(params) => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffPbciStorage::Static(params) => params.is_empty(),
            MmffPbciStorage::Owned(params) => params.is_empty(),
        }
    }

    pub fn get(&self, atom_type: u32) -> Option<&MmffPbci> {
        // RDKit✔️✔️:   const MMFFPBCI *operator()(const unsigned int atomType) const {
        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res = d_params.find(atomType);
        // RDKit❌❌:     return ((res != d_params.end()) ? &((*res).second) : NULL);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     return ((atomType && (atomType <= d_params.size()))
        // RDKit✔️✔️:                 ? &d_params[atomType - 1]
        // RDKit✔️✔️:                 : nullptr);
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:   }
        match &self.params {
            MmffPbciStorage::Static(params) => {
                if atom_type != 0 && (atom_type as usize) <= params.len() {
                    params.get(atom_type as usize - 1).map(|(_, param)| param)
                } else {
                    None
                }
            }
            MmffPbciStorage::Owned(params) => {
                if atom_type != 0 && (atom_type as usize) <= params.len() {
                    params.get(atom_type as usize - 1)
                } else {
                    None
                }
            }
        }
    }
}

/// MMFF bond-charge-increment collection.
#[derive(Clone, Debug, PartialEq)]
pub struct MmffChgCollection {
    source_mmff_chg: Cow<'static, str>,
    params: MmffChgStorage,
}

#[derive(Clone, Debug, PartialEq)]
enum MmffChgStorage {
    Static(&'static [(u8, u8, u8, MmffChg)]),
    Owned {
        params: Vec<MmffChg>,
        bond_type: Vec<u8>,
        i_atom_type: Vec<u8>,
        j_atom_type: Vec<u8>,
    },
}

impl MmffChgCollection {
    pub fn new(mmff_chg: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFChgCollection::MMFFChgCollection(std::string mmffChg) {
        // RDKit✔️✔️:   if (mmffChg.empty()) {
        // RDKit✔️✔️:     mmffChg = defaultMMFFChg;
        // RDKit✔️✔️:   }
        let source_mmff_chg = if mmff_chg.is_empty() {
            Cow::Borrowed(default_mmff_chg()?)
        } else {
            Cow::Owned(mmff_chg.to_owned())
        };
        let params = if mmff_chg.is_empty() {
            MmffChgStorage::Static(DEFAULT_MMFF_CHG_ROWS)
        } else {
            let parsed = parse_mmff_chg(source_mmff_chg.as_ref())?;
            MmffChgStorage::Owned {
                params: parsed.params,
                bond_type: parsed.bond_type,
                i_atom_type: parsed.i_atom_type,
                j_atom_type: parsed.j_atom_type,
            }
        };
        Ok(Self {
            source_mmff_chg,
            params,
        })
    }

    pub fn source_mmff_chg(&self) -> &str {
        self.source_mmff_chg.as_ref()
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffChgStorage::Static(params) => params.len(),
            MmffChgStorage::Owned { params, .. } => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffChgStorage::Static(params) => params.is_empty(),
            MmffChgStorage::Owned { params, .. } => params.is_empty(),
        }
    }

    pub fn get_mmff_chg_params(
        &self,
        bond_type: u32,
        i_atom_type: u32,
        j_atom_type: u32,
    ) -> (i32, Option<&MmffChg>) {
        // RDKit✔️✔️:   const std::pair<int, const MMFFChg *> getMMFFChgParams(
        // RDKit✔️✔️:       const unsigned int bondType, const unsigned int iAtomType,
        // RDKit✔️✔️:       const unsigned int jAtomType) const {
        // RDKit✔️✔️:     int sign = -1;
        // RDKit✔️✔️:     const MMFFChg *mmffChgParams = nullptr;
        // RDKit✔️✔️:     unsigned int canIAtomType = iAtomType;
        // RDKit✔️✔️:     unsigned int canJAtomType = jAtomType;
        // RDKit✔️✔️:     if (iAtomType > jAtomType) {
        // RDKit✔️✔️:       canIAtomType = jAtomType;
        // RDKit✔️✔️:       canJAtomType = iAtomType;
        // RDKit✔️✔️:       sign = 1;
        // RDKit✔️✔️:     }
        let mut sign = -1;
        let mut can_i_atom_type = i_atom_type;
        let mut can_j_atom_type = j_atom_type;
        if i_atom_type > j_atom_type {
            can_i_atom_type = j_atom_type;
            can_j_atom_type = i_atom_type;
            sign = 1;
        }

        if let MmffChgStorage::Static(params) = &self.params {
            let key = (
                bond_type as u8,
                can_i_atom_type as u8,
                can_j_atom_type as u8,
            );
            let param = params
                .binary_search_by_key(&key, |(bond, i_atom, j_atom, _)| (*bond, *i_atom, *j_atom))
                .ok()
                .map(|idx| &params[idx].3);
            return (sign, param);
        }

        let MmffChgStorage::Owned {
            params,
            bond_type: bond_types,
            i_atom_type: i_atom_types,
            j_atom_type: j_atom_types,
        } = &self.params
        else {
            unreachable!("static MMFFChg storage returned above")
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res1 = d_params[bondType].find(canIAtomType);
        // RDKit❌❌:     if (res1 != d_params[bondType].end()) {
        // RDKit❌❌:       const auto res2 = ((*res1).second).find(canJAtomType);
        // RDKit❌❌:       if (res2 != ((*res1).second).end()) {
        // RDKit❌❌:         mmffChgParams = &((*res2).second);
        // RDKit❌❌:       }
        // RDKit❌❌:     }
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto bounds =
        // RDKit✔️✔️:         std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), canIAtomType);
        let Some((i_begin, i_end)) =
            equal_range_u8(i_atom_types, 0, i_atom_types.len(), can_i_atom_type)
        else {
            // RDKit✔️✔️:     return std::make_pair(sign, mmffChgParams);
            return (sign, None);
        };

        // RDKit✔️✔️:     if (bounds.first != bounds.second) {
        // RDKit✔️✔️:       bounds = std::equal_range(
        // RDKit✔️✔️:           d_jAtomType.begin() + (bounds.first - d_iAtomType.begin()),
        // RDKit✔️✔️:           d_jAtomType.begin() + (bounds.second - d_iAtomType.begin()),
        // RDKit✔️✔️:           canJAtomType);
        let Some((j_begin, j_end)) = equal_range_u8(j_atom_types, i_begin, i_end, can_j_atom_type)
        else {
            return (sign, None);
        };

        // RDKit✔️✔️:       if (bounds.first != bounds.second) {
        // RDKit✔️✔️:         bounds = std::equal_range(
        // RDKit✔️✔️:             d_bondType.begin() + (bounds.first - d_jAtomType.begin()),
        // RDKit✔️✔️:             d_bondType.begin() + (bounds.second - d_jAtomType.begin()),
        // RDKit✔️✔️:             bondType);
        let Some((bond_begin, _bond_end)) = equal_range_u8(bond_types, j_begin, j_end, bond_type)
        else {
            return (sign, None);
        };

        // RDKit✔️✔️:         if (bounds.first != bounds.second) {
        // RDKit✔️✔️:           mmffChgParams = &d_params[bounds.first - d_bondType.begin()];
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:
        // RDKit✔️✔️:     return std::make_pair(sign, mmffChgParams);
        // RDKit✔️✔️:   }
        (sign, params.get(bond_begin))
    }
}

/// MMFF bond-stretching parameter collection.
#[derive(Clone, Debug, PartialEq)]
pub struct MmffBondCollection {
    source_mmff_bond: Cow<'static, str>,
    params: MmffBondStorage,
}

#[derive(Clone, Debug, PartialEq)]
enum MmffBondStorage {
    Static(&'static [(u8, u8, u8, MmffBond)]),
    Owned {
        params: Vec<MmffBond>,
        bond_type: Vec<u8>,
        i_atom_type: Vec<u8>,
        j_atom_type: Vec<u8>,
    },
}

impl MmffBondCollection {
    pub fn new(mmff_bond: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFBondCollection::MMFFBondCollection(std::string mmffBond) {
        // RDKit✔️✔️:   if (mmffBond.empty()) {
        // RDKit✔️✔️:     mmffBond = defaultMMFFBond;
        // RDKit✔️✔️:   }
        let source_mmff_bond = if mmff_bond.is_empty() {
            Cow::Borrowed(default_mmff_bond()?)
        } else {
            Cow::Owned(mmff_bond.to_owned())
        };
        let params = if mmff_bond.is_empty() {
            MmffBondStorage::Static(DEFAULT_MMFF_BOND_ROWS)
        } else {
            let parsed = parse_mmff_bond(source_mmff_bond.as_ref())?;
            MmffBondStorage::Owned {
                params: parsed.params,
                bond_type: parsed.bond_type,
                i_atom_type: parsed.i_atom_type,
                j_atom_type: parsed.j_atom_type,
            }
        };
        Ok(Self {
            source_mmff_bond,
            params,
        })
    }

    pub fn source_mmff_bond(&self) -> &str {
        self.source_mmff_bond.as_ref()
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffBondStorage::Static(params) => params.len(),
            MmffBondStorage::Owned { params, .. } => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffBondStorage::Static(params) => params.is_empty(),
            MmffBondStorage::Owned { params, .. } => params.is_empty(),
        }
    }

    pub fn get(&self, bond_type: u32, atom_type: u32, nbr_atom_type: u32) -> Option<&MmffBond> {
        // RDKit✔️✔️:   const MMFFBond *operator()(const unsigned int bondType,
        // RDKit✔️✔️:                              const unsigned int atomType,
        // RDKit✔️✔️:                              const unsigned int nbrAtomType) const {
        // RDKit✔️✔️:     const MMFFBond *mmffBondParams = nullptr;
        // RDKit✔️✔️:     unsigned int canAtomType = atomType;
        // RDKit✔️✔️:     unsigned int canNbrAtomType = nbrAtomType;
        // RDKit✔️✔️:     if (atomType > nbrAtomType) {
        // RDKit✔️✔️:       canAtomType = nbrAtomType;
        // RDKit✔️✔️:       canNbrAtomType = atomType;
        // RDKit✔️✔️:     }
        let mut can_atom_type = atom_type;
        let mut can_nbr_atom_type = nbr_atom_type;
        if atom_type > nbr_atom_type {
            can_atom_type = nbr_atom_type;
            can_nbr_atom_type = atom_type;
        }

        if let MmffBondStorage::Static(params) = &self.params {
            let key = (
                bond_type as u8,
                can_atom_type as u8,
                can_nbr_atom_type as u8,
            );
            return params
                .binary_search_by_key(&key, |(bond, atom, nbr_atom, _)| (*bond, *atom, *nbr_atom))
                .ok()
                .map(|idx| &params[idx].3);
        }

        let MmffBondStorage::Owned {
            params,
            bond_type: bond_types,
            i_atom_type: i_atom_types,
            j_atom_type: j_atom_types,
        } = &self.params
        else {
            unreachable!("static MMFFBond storage returned above")
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res1 = d_params.find(bondType);
        // RDKit❌❌:     std::map<const unsigned int,
        // RDKit❌❌:              std::map<const unsigned int, MMFFBond>>::const_iterator res2;
        // RDKit❌❌:     std::map<const unsigned int, MMFFBond>::const_iterator res3;
        // RDKit❌❌:     if (res1 != d_params.end()) {
        // RDKit❌❌:       res2 = ((*res1).second).find(canAtomType);
        // RDKit❌❌:       if (res2 != ((*res1).second).end()) {
        // RDKit❌❌:         res3 = ((*res2).second).find(canNbrAtomType);
        // RDKit❌❌:         if (res3 != ((*res2).second).end()) {
        // RDKit❌❌:           mmffBondParams = &((*res3).second);
        // RDKit❌❌:         }
        // RDKit❌❌:       }
        // RDKit❌❌:     }
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto bounds =
        // RDKit✔️✔️:         std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), canAtomType);
        let (i_begin, i_end) = equal_range_u8(i_atom_types, 0, i_atom_types.len(), can_atom_type)?;

        // RDKit✔️✔️:     if (bounds.first != bounds.second) {
        // RDKit✔️✔️:       bounds = std::equal_range(
        // RDKit✔️✔️:           d_jAtomType.begin() + (bounds.first - d_iAtomType.begin()),
        // RDKit✔️✔️:           d_jAtomType.begin() + (bounds.second - d_iAtomType.begin()),
        // RDKit✔️✔️:           canNbrAtomType);
        let (j_begin, j_end) = equal_range_u8(j_atom_types, i_begin, i_end, can_nbr_atom_type)?;

        // RDKit✔️✔️:       if (bounds.first != bounds.second) {
        // RDKit✔️✔️:         bounds = std::equal_range(
        // RDKit✔️✔️:             d_bondType.begin() + (bounds.first - d_jAtomType.begin()),
        // RDKit✔️✔️:             d_bondType.begin() + (bounds.second - d_jAtomType.begin()),
        // RDKit✔️✔️:             bondType);
        let (bond_begin, _bond_end) = equal_range_u8(bond_types, j_begin, j_end, bond_type)?;

        // RDKit✔️✔️:         if (bounds.first != bounds.second) {
        // RDKit✔️✔️:           mmffBondParams = &d_params[bounds.first - d_bondType.begin()];
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:
        // RDKit✔️✔️:     return mmffBondParams;
        // RDKit✔️✔️:   }
        params.get(bond_begin)
    }
}

/// MMFF angle-bending parameter collection.
#[derive(Clone, Debug, PartialEq)]
pub struct MmffAngleCollection {
    source_mmff_angle: Cow<'static, str>,
    params: MmffAngleStorage,
}

#[derive(Clone, Debug, PartialEq)]
enum MmffAngleStorage {
    Static(&'static [(u8, u8, u8, u8, MmffAngle)]),
    Owned {
        params: Vec<MmffAngle>,
        angle_type: Vec<u8>,
        i_atom_type: Vec<u8>,
        j_atom_type: Vec<u8>,
        k_atom_type: Vec<u8>,
    },
}

impl MmffAngleCollection {
    pub fn new(mmff_angle: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFAngleCollection::MMFFAngleCollection(std::string mmffAngle) {
        // RDKit✔️✔️:   if (mmffAngle.empty()) {
        // RDKit✔️✔️:     unsigned int i = 0;
        // RDKit✔️✔️:     while (defaultMMFFAngleData[i] != "EOS") {
        // RDKit✔️✔️:       mmffAngle += defaultMMFFAngleData[i];
        // RDKit✔️✔️:       ++i;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        let source_mmff_angle = if mmff_angle.is_empty() {
            Cow::Borrowed(default_mmff_angle()?)
        } else {
            Cow::Owned(mmff_angle.to_owned())
        };
        let params = if mmff_angle.is_empty() {
            MmffAngleStorage::Static(DEFAULT_MMFF_ANGLE_ROWS)
        } else {
            let parsed = parse_mmff_angle(source_mmff_angle.as_ref())?;
            MmffAngleStorage::Owned {
                params: parsed.params,
                angle_type: parsed.angle_type,
                i_atom_type: parsed.i_atom_type,
                j_atom_type: parsed.j_atom_type,
                k_atom_type: parsed.k_atom_type,
            }
        };
        Ok(Self {
            source_mmff_angle,
            params,
        })
    }

    pub fn source_mmff_angle(&self) -> &str {
        self.source_mmff_angle.as_ref()
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffAngleStorage::Static(params) => params.len(),
            MmffAngleStorage::Owned { params, .. } => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffAngleStorage::Static(params) => params.is_empty(),
            MmffAngleStorage::Owned { params, .. } => params.is_empty(),
        }
    }

    pub fn get(
        &self,
        mmff_def: &MmffDefCollection,
        angle_type: u32,
        i_atom_type: u32,
        j_atom_type: u32,
        k_atom_type: u32,
    ) -> Option<&MmffAngle> {
        // RDKit✔️✔️:   const MMFFAngle *operator()(const MMFFDefCollection *mmffDef,
        // RDKit✔️✔️:                               const unsigned int angleType,
        // RDKit✔️✔️:                               const unsigned int iAtomType,
        // RDKit✔️✔️:                               const unsigned int jAtomType,
        // RDKit✔️✔️:                               const unsigned int kAtomType) const {
        // RDKit✔️✔️:     const MMFFAngle *mmffAngleParams = nullptr;
        // RDKit✔️✔️:     unsigned int iter = 0;
        let mut iter = 0_usize;

        if let MmffAngleStorage::Static(params) = &self.params {
            let i_def = mmff_def.get(i_atom_type)?;
            let k_def = mmff_def.get(k_atom_type)?;
            while iter < 4 {
                let mut can_i_atom_type = u32::from(i_def.eq_level[iter]);
                let mut can_k_atom_type = u32::from(k_def.eq_level[iter]);
                if can_i_atom_type > can_k_atom_type {
                    std::mem::swap(&mut can_i_atom_type, &mut can_k_atom_type);
                }
                let key = (
                    angle_type as u8,
                    can_i_atom_type as u8,
                    j_atom_type as u8,
                    can_k_atom_type as u8,
                );
                if let Ok(idx) =
                    params.binary_search_by_key(&key, |(angle, i, j, k, _)| (*angle, *i, *j, *k))
                {
                    return Some(&params[idx].4);
                }
                iter += 1;
            }
            return None;
        }

        let MmffAngleStorage::Owned {
            params,
            angle_type: angle_types,
            i_atom_type: i_atom_types,
            j_atom_type: j_atom_types,
            k_atom_type: k_atom_types,
        } = &self.params
        else {
            unreachable!("static MMFFAngle storage returned above")
        };

        // RDKit✔️✔️: // For bending of the i-j-k angle, a five-stage process based
        // RDKit✔️✔️: // in the level combinations 1-1-1,2-2-2,3-2-3,4-2-4, and
        // RDKit✔️✔️: // 5-2-5 is used. (MMFF.I, note 68, page 519)
        // RDKit✔️✔️: // We skip 1-1-1 since Level 2 === Level 1
        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     while ((iter < 4) && (!mmffAngleParams)) {
        // RDKit❌❌:       unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iter];
        // RDKit❌❌:       unsigned int canKAtomType = (*mmffDef)(kAtomType)->eqLevel[iter];
        // RDKit❌❌:       if (canIAtomType > canKAtomType) {
        // RDKit❌❌:         std::swap(canIAtomType, canKAtomType);
        // RDKit❌❌:       }
        // RDKit❌❌:       const auto res1 = d_params.find(angleType);
        // RDKit❌❌:       if (res1 != d_params.end()) {
        // RDKit❌❌:         const auto res2 = ((*res1).second).find(canIAtomType);
        // RDKit❌❌:         if (res2 != ((*res1).second).end()) {
        // RDKit❌❌:           const auto res3 = ((*res2).second).find(jAtomType);
        // RDKit❌❌:           if (res3 != ((*res2).second).end()) {
        // RDKit❌❌:             const auto res4 = ((*res3).second).find(canKAtomType);
        // RDKit❌❌:             if (res4 != ((*res3).second).end()) {
        // RDKit❌❌:               mmffAngleParams = &((*res4).second);
        // RDKit❌❌:             }
        // RDKit❌❌:           }
        // RDKit❌❌:         }
        // RDKit❌❌:       }
        // RDKit❌❌:       ++iter;
        // RDKit❌❌:     }
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto jBounds =
        // RDKit✔️✔️:         std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
        let (j_begin, j_end) = equal_range_u8(j_atom_types, 0, j_atom_types.len(), j_atom_type)?;

        let i_def = mmff_def.get(i_atom_type)?;
        let k_def = mmff_def.get(k_atom_type)?;

        // RDKit✔️✔️:     if (jBounds.first != jBounds.second) {
        // RDKit✔️✔️:       while ((iter < 4) && (!mmffAngleParams)) {
        while iter < 4 {
            // RDKit✔️✔️:         unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iter];
            // RDKit✔️✔️:         unsigned int canKAtomType = (*mmffDef)(kAtomType)->eqLevel[iter];
            // RDKit✔️✔️:         if (canIAtomType > canKAtomType) {
            // RDKit✔️✔️:           std::swap(canIAtomType, canKAtomType);
            // RDKit✔️✔️:         }
            let mut can_i_atom_type = u32::from(i_def.eq_level[iter]);
            let mut can_k_atom_type = u32::from(k_def.eq_level[iter]);
            if can_i_atom_type > can_k_atom_type {
                std::mem::swap(&mut can_i_atom_type, &mut can_k_atom_type);
            }

            // RDKit✔️✔️:         auto bounds = std::equal_range(
            // RDKit✔️✔️:             d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
            // RDKit✔️✔️:             d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
            // RDKit✔️✔️:             canIAtomType);
            if let Some((i_begin, i_end)) =
                equal_range_u8(i_atom_types, j_begin, j_end, can_i_atom_type)
            {
                // RDKit✔️✔️:         if (bounds.first != bounds.second) {
                // RDKit✔️✔️:           bounds = std::equal_range(
                // RDKit✔️✔️:               d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
                // RDKit✔️✔️:               d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
                // RDKit✔️✔️:               canKAtomType);
                if let Some((k_begin, k_end)) =
                    equal_range_u8(k_atom_types, i_begin, i_end, can_k_atom_type)
                {
                    // RDKit✔️✔️:           if (bounds.first != bounds.second) {
                    // RDKit✔️✔️:             bounds = std::equal_range(
                    // RDKit✔️✔️:                 d_angleType.begin() + (bounds.first - d_kAtomType.begin()),
                    // RDKit✔️✔️:                 d_angleType.begin() + (bounds.second - d_kAtomType.begin()),
                    // RDKit✔️✔️:                 angleType);
                    if let Some((angle_begin, _angle_end)) =
                        equal_range_u8(angle_types, k_begin, k_end, angle_type)
                    {
                        // RDKit✔️✔️:             if (bounds.first != bounds.second) {
                        // RDKit✔️✔️:               mmffAngleParams = &d_params[bounds.first - d_angleType.begin()];
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:           }
                        // RDKit✔️✔️:         }
                        // RDKit✔️✔️:         ++iter;
                        // RDKit✔️✔️:       }
                        // RDKit✔️✔️:     }
                        // RDKit✔️✔️: #endif
                        // RDKit✔️✔️:
                        // RDKit✔️✔️:     return mmffAngleParams;
                        // RDKit✔️✔️:   }
                        return params.get(angle_begin);
                    }
                }
            }
            iter += 1;
        }
        None
    }
}

/// MMFF stretch-bend parameter collection.
#[derive(Clone, Debug, PartialEq)]
pub struct MmffStbnCollection {
    source_mmff_stbn: Cow<'static, str>,
    params: MmffStbnStorage,
}

#[derive(Clone, Debug, PartialEq)]
enum MmffStbnStorage {
    Static(&'static [(u8, u8, u8, u8, MmffStbn)]),
    Owned {
        params: Vec<MmffStbn>,
        stretch_bend_type: Vec<u8>,
        i_atom_type: Vec<u8>,
        j_atom_type: Vec<u8>,
        k_atom_type: Vec<u8>,
    },
}

impl MmffStbnCollection {
    pub fn new(mmff_stbn: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFStbnCollection::MMFFStbnCollection(std::string mmffStbn) {
        // RDKit✔️✔️:   if (mmffStbn.empty()) {
        // RDKit✔️✔️:     mmffStbn = defaultMMFFStbn;
        // RDKit✔️✔️:   }
        let source_mmff_stbn = if mmff_stbn.is_empty() {
            Cow::Borrowed(default_mmff_stbn()?)
        } else {
            Cow::Owned(mmff_stbn.to_owned())
        };
        let params = if mmff_stbn.is_empty() {
            MmffStbnStorage::Static(DEFAULT_MMFF_STBN_ROWS)
        } else {
            let parsed = parse_mmff_stbn(source_mmff_stbn.as_ref())?;
            MmffStbnStorage::Owned {
                params: parsed.params,
                stretch_bend_type: parsed.stretch_bend_type,
                i_atom_type: parsed.i_atom_type,
                j_atom_type: parsed.j_atom_type,
                k_atom_type: parsed.k_atom_type,
            }
        };
        Ok(Self {
            source_mmff_stbn,
            params,
        })
    }

    pub fn source_mmff_stbn(&self) -> &str {
        self.source_mmff_stbn.as_ref()
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffStbnStorage::Static(params) => params.len(),
            MmffStbnStorage::Owned { params, .. } => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffStbnStorage::Static(params) => params.is_empty(),
            MmffStbnStorage::Owned { params, .. } => params.is_empty(),
        }
    }

    pub fn get_mmff_stbn_params(
        &self,
        stretch_bend_type: u32,
        bond_type1: u32,
        bond_type2: u32,
        i_atom_type: u32,
        j_atom_type: u32,
        k_atom_type: u32,
    ) -> (bool, Option<&MmffStbn>) {
        // RDKit✔️✔️:   const std::pair<bool, const MMFFStbn *> getMMFFStbnParams(
        // RDKit✔️✔️:       const unsigned int stretchBendType, const unsigned int bondType1,
        // RDKit✔️✔️:       const unsigned int bondType2, const unsigned int iAtomType,
        // RDKit✔️✔️:       const unsigned int jAtomType, const unsigned int kAtomType) const {
        // RDKit✔️✔️:     const MMFFStbn *mmffStbnParams = nullptr;
        // RDKit✔️✔️:     bool swap = false;
        // RDKit✔️✔️:     unsigned int canIAtomType = iAtomType;
        // RDKit✔️✔️:     unsigned int canKAtomType = kAtomType;
        // RDKit✔️✔️:     unsigned int canStretchBendType = stretchBendType;
        let mut swap = false;
        let mut can_i_atom_type = i_atom_type;
        let mut can_k_atom_type = k_atom_type;
        let can_stretch_bend_type = stretch_bend_type;

        // RDKit✔️✔️:     if (iAtomType > kAtomType) {
        // RDKit✔️✔️:       canIAtomType = kAtomType;
        // RDKit✔️✔️:       canKAtomType = iAtomType;
        // RDKit✔️✔️:       swap = true;
        // RDKit✔️✔️:     } else if (iAtomType == kAtomType) {
        // RDKit✔️✔️:       swap = (bondType1 < bondType2);
        // RDKit✔️✔️:     }
        if i_atom_type > k_atom_type {
            can_i_atom_type = k_atom_type;
            can_k_atom_type = i_atom_type;
            swap = true;
        } else if i_atom_type == k_atom_type {
            swap = bond_type1 < bond_type2;
        }

        if let MmffStbnStorage::Static(params) = &self.params {
            let key = (
                can_stretch_bend_type as u8,
                can_i_atom_type as u8,
                j_atom_type as u8,
                can_k_atom_type as u8,
            );
            let param = params
                .binary_search_by_key(&key, |(stbn, i, j, k, _)| (*stbn, *i, *j, *k))
                .ok()
                .map(|idx| &params[idx].4);
            return (swap, param);
        }

        let MmffStbnStorage::Owned {
            params,
            stretch_bend_type: stretch_bend_types,
            i_atom_type: i_atom_types,
            j_atom_type: j_atom_types,
            k_atom_type: k_atom_types,
        } = &self.params
        else {
            unreachable!("static MMFFStbn storage returned above")
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res1 = d_params.find(canStretchBendType);
        // RDKit❌❌:     if (res1 != d_params.end()) {
        // RDKit❌❌:       const auto res2 = ((*res1).second).find(canIAtomType);
        // RDKit❌❌:       if (res2 != ((*res1).second).end()) {
        // RDKit❌❌:         const auto res3 = ((*res2).second).find(jAtomType);
        // RDKit❌❌:         if (res3 != ((*res2).second).end()) {
        // RDKit❌❌:           const auto res4 = ((*res3).second).find(canKAtomType);
        // RDKit❌❌:           if (res4 != ((*res3).second).end()) {
        // RDKit❌❌:             mmffStbnParams = &((*res4).second);
        // RDKit❌❌:           }
        // RDKit❌❌:         }
        // RDKit❌❌:       }
        // RDKit❌❌:     }
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto jBounds =
        // RDKit✔️✔️:         std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
        let Some((j_begin, j_end)) =
            equal_range_u8(j_atom_types, 0, j_atom_types.len(), j_atom_type)
        else {
            // RDKit✔️✔️:     return std::make_pair(swap, mmffStbnParams);
            return (swap, None);
        };

        // RDKit✔️✔️:     if (jBounds.first != jBounds.second) {
        // RDKit✔️✔️:       auto bounds = std::equal_range(
        // RDKit✔️✔️:           d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
        // RDKit✔️✔️:           d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
        // RDKit✔️✔️:           canIAtomType);
        let Some((i_begin, i_end)) = equal_range_u8(i_atom_types, j_begin, j_end, can_i_atom_type)
        else {
            return (swap, None);
        };

        // RDKit✔️✔️:       if (bounds.first != bounds.second) {
        // RDKit✔️✔️:         bounds = std::equal_range(
        // RDKit✔️✔️:             d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
        // RDKit✔️✔️:             d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
        // RDKit✔️✔️:             canKAtomType);
        let Some((k_begin, k_end)) = equal_range_u8(k_atom_types, i_begin, i_end, can_k_atom_type)
        else {
            return (swap, None);
        };

        // RDKit✔️✔️:         if (bounds.first != bounds.second) {
        // RDKit✔️✔️:           bounds = std::equal_range(
        // RDKit✔️✔️:               d_stretchBendType.begin() + (bounds.first - d_kAtomType.begin()),
        // RDKit✔️✔️:               d_stretchBendType.begin() + (bounds.second - d_kAtomType.begin()),
        // RDKit✔️✔️:               canStretchBendType);
        let Some((stbn_begin, _stbn_end)) =
            equal_range_u8(stretch_bend_types, k_begin, k_end, can_stretch_bend_type)
        else {
            return (swap, None);
        };

        // RDKit✔️✔️:           if (bounds.first != bounds.second) {
        // RDKit✔️✔️:             mmffStbnParams =
        // RDKit✔️✔️:                 &d_params[bounds.first - d_stretchBendType.begin()];
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:
        // RDKit✔️✔️:     return std::make_pair(swap, mmffStbnParams);
        // RDKit✔️✔️:   }
        (swap, params.get(stbn_begin))
    }
}

/// MMFF default stretch-bend parameter collection.
#[derive(Clone, Debug, PartialEq)]
pub struct MmffDfsbCollection {
    source_mmff_dfsb: Cow<'static, str>,
    params: MmffDfsbStorage,
}

#[derive(Clone, Debug, PartialEq)]
enum MmffDfsbStorage {
    Static(&'static [(u32, u32, u32, MmffStbn)]),
    Owned(BTreeMap<u32, BTreeMap<u32, BTreeMap<u32, MmffStbn>>>),
}

impl MmffDfsbCollection {
    pub fn new(mmff_dfsb: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFDfsbCollection::MMFFDfsbCollection(std::string mmffDfsb) {
        // RDKit✔️✔️:   if (mmffDfsb.empty()) {
        // RDKit✔️✔️:     mmffDfsb = defaultMMFFDfsb;
        // RDKit✔️✔️:   }
        let source_mmff_dfsb = if mmff_dfsb.is_empty() {
            Cow::Borrowed(default_mmff_dfsb()?)
        } else {
            Cow::Owned(mmff_dfsb.to_owned())
        };
        let params = if mmff_dfsb.is_empty() {
            MmffDfsbStorage::Static(DEFAULT_MMFF_DFSB_ROWS)
        } else {
            MmffDfsbStorage::Owned(parse_mmff_dfsb(source_mmff_dfsb.as_ref())?)
        };
        Ok(Self {
            source_mmff_dfsb,
            params,
        })
    }

    pub fn source_mmff_dfsb(&self) -> &str {
        self.source_mmff_dfsb.as_ref()
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffDfsbStorage::Static(params) => params.len(),
            MmffDfsbStorage::Owned(params) => params
                .values()
                .map(|row2| row2.values().map(BTreeMap::len).sum::<usize>())
                .sum(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffDfsbStorage::Static(params) => params.is_empty(),
            MmffDfsbStorage::Owned(params) => params.is_empty(),
        }
    }

    pub fn get_mmff_dfsb_params(
        &self,
        periodic_table_row1: u32,
        periodic_table_row2: u32,
        periodic_table_row3: u32,
    ) -> (bool, Option<&MmffStbn>) {
        // RDKit✔️✔️:   const std::pair<bool, const MMFFStbn *> getMMFFDfsbParams(
        // RDKit✔️✔️:       const unsigned int periodicTableRow1,
        // RDKit✔️✔️:       const unsigned int periodicTableRow2,
        // RDKit✔️✔️:       const unsigned int periodicTableRow3) const {
        // RDKit✔️✔️:     const MMFFStbn *mmffDfsbParams = nullptr;
        // RDKit✔️✔️:     bool swap = false;
        // RDKit✔️✔️:     unsigned int canPeriodicTableRow1 = periodicTableRow1;
        // RDKit✔️✔️:     unsigned int canPeriodicTableRow3 = periodicTableRow3;
        let mut swap = false;
        let mut can_periodic_table_row1 = periodic_table_row1;
        let mut can_periodic_table_row3 = periodic_table_row3;

        // RDKit✔️✔️:     if (periodicTableRow1 > periodicTableRow3) {
        // RDKit✔️✔️:       canPeriodicTableRow1 = periodicTableRow3;
        // RDKit✔️✔️:       canPeriodicTableRow3 = periodicTableRow1;
        // RDKit✔️✔️:       swap = true;
        // RDKit✔️✔️:     }
        if periodic_table_row1 > periodic_table_row3 {
            can_periodic_table_row1 = periodic_table_row3;
            can_periodic_table_row3 = periodic_table_row1;
            swap = true;
        }

        // RDKit✔️✔️:     const auto res1 = d_params.find(canPeriodicTableRow1);
        // RDKit✔️✔️:     if (res1 != d_params.end()) {
        // RDKit✔️✔️:       const auto res2 = ((*res1).second).find(periodicTableRow2);
        // RDKit✔️✔️:       if (res2 != ((*res1).second).end()) {
        // RDKit✔️✔️:         const auto res3 = ((*res2).second).find(canPeriodicTableRow3);
        // RDKit✔️✔️:         if (res3 != ((*res2).second).end()) {
        // RDKit✔️✔️:           mmffDfsbParams = &((*res3).second);
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     return std::make_pair(swap, mmffDfsbParams);
        // RDKit✔️✔️:   }
        let params = match &self.params {
            MmffDfsbStorage::Static(params) => params
                .binary_search_by_key(
                    &(
                        can_periodic_table_row1,
                        periodic_table_row2,
                        can_periodic_table_row3,
                    ),
                    |(row1, row2, row3, _)| (*row1, *row2, *row3),
                )
                .ok()
                .map(|idx| &params[idx].3),
            MmffDfsbStorage::Owned(params) => params
                .get(&can_periodic_table_row1)
                .and_then(|row2| row2.get(&periodic_table_row2))
                .and_then(|row3| row3.get(&can_periodic_table_row3)),
        };
        (swap, params)
    }
}

/// MMFF out-of-plane bending parameter collection.
#[derive(Clone, Debug, PartialEq)]
pub struct MmffOopCollection {
    source_mmff_oop: Cow<'static, str>,
    params: MmffOopStorage,
}

#[derive(Clone, Debug, PartialEq)]
enum MmffOopStorage {
    Static(&'static [(u8, u8, u8, u8, MmffOop)]),
    Owned {
        params: Vec<MmffOop>,
        i_atom_type: Vec<u8>,
        j_atom_type: Vec<u8>,
        k_atom_type: Vec<u8>,
        l_atom_type: Vec<u8>,
    },
}

impl MmffOopCollection {
    pub fn new(is_mmffs: bool, mmff_oop: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFOopCollection::MMFFOopCollection(const bool isMMFFs, std::string mmffOop) {
        // RDKit✔️✔️:   if (mmffOop.empty()) {
        // RDKit✔️✔️:     mmffOop = (isMMFFs ? defaultMMFFsOop : defaultMMFFOop);
        // RDKit✔️✔️:   }
        let source_mmff_oop = if mmff_oop.is_empty() {
            if is_mmffs {
                Cow::Borrowed(default_mmffs_oop()?)
            } else {
                Cow::Borrowed(default_mmff_oop()?)
            }
        } else {
            Cow::Owned(mmff_oop.to_owned())
        };
        let params = if mmff_oop.is_empty() {
            if is_mmffs {
                MmffOopStorage::Static(DEFAULT_MMFFS_OOP_ROWS)
            } else {
                MmffOopStorage::Static(DEFAULT_MMFF_OOP_ROWS)
            }
        } else {
            let parsed = parse_mmff_oop(source_mmff_oop.as_ref())?;
            MmffOopStorage::Owned {
                params: parsed.params,
                i_atom_type: parsed.i_atom_type,
                j_atom_type: parsed.j_atom_type,
                k_atom_type: parsed.k_atom_type,
                l_atom_type: parsed.l_atom_type,
            }
        };
        Ok(Self {
            source_mmff_oop,
            params,
        })
    }

    pub fn source_mmff_oop(&self) -> &str {
        self.source_mmff_oop.as_ref()
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffOopStorage::Static(params) => params.len(),
            MmffOopStorage::Owned { params, .. } => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffOopStorage::Static(params) => params.is_empty(),
            MmffOopStorage::Owned { params, .. } => params.is_empty(),
        }
    }

    pub fn get_mmff_oop_params(
        &self,
        mmff_def: &MmffDefCollection,
        i_atom_type: u32,
        j_atom_type: u32,
        k_atom_type: u32,
        l_atom_type: u32,
    ) -> Option<&MmffOop> {
        // RDKit✔️✔️:   const MMFFOop *operator()(const MMFFDefCollection *mmffDef,
        // RDKit✔️✔️:                             const unsigned int iAtomType,
        // RDKit✔️✔️:                             const unsigned int jAtomType,
        // RDKit✔️✔️:                             const unsigned int kAtomType,
        // RDKit✔️✔️:                             const unsigned int lAtomType) const {
        // RDKit✔️✔️:     const MMFFOop *mmffOopParams = nullptr;
        // RDKit✔️✔️:     unsigned int iter = 0;
        // RDKit✔️✔️:     std::vector<unsigned int> canIKLAtomType(3);
        let mut iter = 0_usize;

        // RDKit✔️✔️: // For out-of-plane bending ijk; I , where j is the central
        // RDKit✔️✔️: // atom [cf. eq. (511, the five-stage protocol 1-1-1; 1, 2-2-2; 2,
        // RDKit✔️✔️: // 3-2-3;3, 4-2-4;4, 5-2-5;5 is used. The final stage provides
        // RDKit✔️✔️: // wild-card defaults for all except the central atom.
        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     while ((iter < 4) && (!mmffOopParams)) {
        // RDKit❌❌:       canIKLAtomType[0] = (*mmffDef)(iAtomType)->eqLevel[iter];
        // RDKit❌❌:       unsigned int canJAtomType = jAtomType;
        // RDKit❌❌:       canIKLAtomType[1] = (*mmffDef)(kAtomType)->eqLevel[iter];
        // RDKit❌❌:       canIKLAtomType[2] = (*mmffDef)(lAtomType)->eqLevel[iter];
        // RDKit❌❌:       std::sort(canIKLAtomType.begin(), canIKLAtomType.end());
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto jBounds =
        // RDKit✔️✔️:         std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
        if let MmffOopStorage::Static(params) = &self.params {
            let mut iter = 0_usize;
            while iter < 4 {
                let mut can_ikl_atom_type = [
                    u32::from(mmff_def.get(i_atom_type)?.eq_level[iter]),
                    u32::from(mmff_def.get(k_atom_type)?.eq_level[iter]),
                    u32::from(mmff_def.get(l_atom_type)?.eq_level[iter]),
                ];
                can_ikl_atom_type.sort_unstable();
                let key = (
                    can_ikl_atom_type[0] as u8,
                    j_atom_type as u8,
                    can_ikl_atom_type[1] as u8,
                    can_ikl_atom_type[2] as u8,
                );
                if let Ok(idx) =
                    params.binary_search_by_key(&key, |(i, j, k, l, _)| (*i, *j, *k, *l))
                {
                    return Some(&params[idx].4);
                }
                iter += 1;
            }
            return None;
        }

        let MmffOopStorage::Owned {
            params,
            i_atom_type: i_atom_types,
            j_atom_type: j_atom_types,
            k_atom_type: k_atom_types,
            l_atom_type: l_atom_types,
        } = &self.params
        else {
            unreachable!("static MMFFOop storage returned above")
        };

        let (j_begin, j_end) = equal_range_u8(j_atom_types, 0, j_atom_types.len(), j_atom_type)?;

        // RDKit✔️✔️:     if (jBounds.first != jBounds.second) {
        // RDKit✔️✔️:       while ((iter < 4) && (!mmffOopParams)) {
        while iter < 4 {
            // RDKit✔️✔️:         canIKLAtomType[0] = (*mmffDef)(iAtomType)->eqLevel[iter];
            // RDKit✔️✔️:         canIKLAtomType[1] = (*mmffDef)(kAtomType)->eqLevel[iter];
            // RDKit✔️✔️:         canIKLAtomType[2] = (*mmffDef)(lAtomType)->eqLevel[iter];
            let mut can_ikl_atom_type = [
                u32::from(mmff_def.get(i_atom_type)?.eq_level[iter]),
                u32::from(mmff_def.get(k_atom_type)?.eq_level[iter]),
                u32::from(mmff_def.get(l_atom_type)?.eq_level[iter]),
            ];

            // RDKit✔️✔️:         std::sort(canIKLAtomType.begin(), canIKLAtomType.end());
            can_ikl_atom_type.sort_unstable();

            // RDKit✔️✔️:         auto bounds = std::equal_range(
            // RDKit✔️✔️:             d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
            // RDKit✔️✔️:             d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
            // RDKit✔️✔️:             canIKLAtomType[0]);
            if let Some((i_begin, i_end)) =
                equal_range_u8(i_atom_types, j_begin, j_end, can_ikl_atom_type[0])
            {
                // RDKit✔️✔️:         if (bounds.first != bounds.second) {
                // RDKit✔️✔️:           bounds = std::equal_range(
                // RDKit✔️✔️:               d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
                // RDKit✔️✔️:               d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
                // RDKit✔️✔️:               canIKLAtomType[1]);
                if let Some((k_begin, k_end)) =
                    equal_range_u8(k_atom_types, i_begin, i_end, can_ikl_atom_type[1])
                {
                    // RDKit✔️✔️:           if (bounds.first != bounds.second) {
                    // RDKit✔️✔️:             bounds = std::equal_range(
                    // RDKit✔️✔️:                 d_lAtomType.begin() + (bounds.first - d_kAtomType.begin()),
                    // RDKit✔️✔️:                 d_lAtomType.begin() + (bounds.second - d_kAtomType.begin()),
                    // RDKit✔️✔️:                 canIKLAtomType[2]);
                    if let Some((l_begin, _l_end)) =
                        equal_range_u8(l_atom_types, k_begin, k_end, can_ikl_atom_type[2])
                    {
                        // RDKit✔️✔️:             if (bounds.first != bounds.second) {
                        // RDKit✔️✔️:               mmffOopParams = &d_params[bounds.first - d_lAtomType.begin()];
                        // RDKit✔️✔️:             }
                        return params.get(l_begin);
                    }
                    // RDKit✔️✔️:           }
                }
                // RDKit✔️✔️:         }
            }

            // RDKit✔️✔️:         ++iter;
            iter += 1;
        }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:
        // RDKit✔️✔️:     return mmffOopParams;
        // RDKit✔️✔️:   }
        None
    }
}

/// MMFF torsion parameter collection.
#[derive(Clone, Debug, PartialEq)]
pub struct MmffTorCollection {
    source_mmff_tor: Cow<'static, str>,
    params: MmffTorStorage,
}

#[derive(Clone, Debug, PartialEq)]
enum MmffTorStorage {
    Static(&'static [(u8, u8, u8, u8, u8, MmffTor)]),
    Owned {
        params: Vec<MmffTor>,
        tor_type: Vec<u8>,
        i_atom_type: Vec<u8>,
        j_atom_type: Vec<u8>,
        k_atom_type: Vec<u8>,
        l_atom_type: Vec<u8>,
    },
}

impl MmffTorCollection {
    pub fn new(is_mmffs: bool, mmff_tor: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFTorCollection::MMFFTorCollection(const bool isMMFFs, std::string mmffTor) {
        // RDKit✔️✔️:   if (mmffTor.empty()) {
        // RDKit✔️✔️:     mmffTor = (isMMFFs ? defaultMMFFsTor : defaultMMFFTor);
        // RDKit✔️✔️:   }
        let source_mmff_tor = if mmff_tor.is_empty() {
            if is_mmffs {
                Cow::Borrowed(default_mmffs_tor()?)
            } else {
                Cow::Borrowed(default_mmff_tor()?)
            }
        } else {
            Cow::Owned(mmff_tor.to_owned())
        };
        let params = if mmff_tor.is_empty() {
            if is_mmffs {
                MmffTorStorage::Static(DEFAULT_MMFFS_TOR_ROWS)
            } else {
                MmffTorStorage::Static(DEFAULT_MMFF_TOR_ROWS)
            }
        } else {
            let parsed = parse_mmff_tor(source_mmff_tor.as_ref())?;
            MmffTorStorage::Owned {
                params: parsed.params,
                tor_type: parsed.tor_type,
                i_atom_type: parsed.i_atom_type,
                j_atom_type: parsed.j_atom_type,
                k_atom_type: parsed.k_atom_type,
                l_atom_type: parsed.l_atom_type,
            }
        };
        Ok(Self {
            source_mmff_tor,
            params,
        })
    }

    pub fn source_mmff_tor(&self) -> &str {
        self.source_mmff_tor.as_ref()
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffTorStorage::Static(params) => params.len(),
            MmffTorStorage::Owned { params, .. } => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffTorStorage::Static(params) => params.is_empty(),
            MmffTorStorage::Owned { params, .. } => params.is_empty(),
        }
    }

    pub fn get_mmff_tor_params(
        &self,
        mmff_def: &MmffDefCollection,
        tor_type: (u32, u32),
        i_atom_type: u32,
        j_atom_type: u32,
        k_atom_type: u32,
        l_atom_type: u32,
    ) -> (u32, Option<&MmffTor>) {
        // RDKit✔️✔️:   const std::pair<const unsigned int, const MMFFTor *> getMMFFTorParams(
        // RDKit✔️✔️:       const MMFFDefCollection *mmffDef,
        // RDKit✔️✔️:       const std::pair<unsigned int, unsigned int> torType,
        // RDKit✔️✔️:       const unsigned int iAtomType, const unsigned int jAtomType,
        // RDKit✔️✔️:       const unsigned int kAtomType, const unsigned int lAtomType) const {
        // RDKit✔️✔️:     const MMFFTor *mmffTorParams = nullptr;
        // RDKit✔️✔️:     unsigned int iter = 0;
        // RDKit✔️✔️:     unsigned int iWildCard = 0;
        // RDKit✔️✔️:     unsigned int lWildCard = 0;
        // RDKit✔️✔️:     unsigned int canTorType = torType.first;
        // RDKit✔️✔️:     unsigned int maxIter = 5;
        if let MmffTorStorage::Static(params) = &self.params {
            let mut mmff_tor_params = None;
            let mut iter = 0_usize;
            let mut can_tor_type = tor_type.0;
            let mut max_iter = 5_usize;
            while ((iter < max_iter) && (mmff_tor_params.is_none() || max_iter == 4))
                || ((iter == 4) && (tor_type.0 == 5) && (tor_type.1 != 0))
            {
                if (max_iter == 5) && (iter == 4) {
                    max_iter = 4;
                    iter = 0;
                    can_tor_type = tor_type.1;
                }
                let mut i_wild_card = iter;
                let mut l_wild_card = iter;
                if iter == 1 {
                    i_wild_card = 1;
                    l_wild_card = 3;
                } else if iter == 2 {
                    i_wild_card = 3;
                    l_wild_card = 1;
                }
                let Some(can_i_atom_type) =
                    mmff_def.torsion_equivalence_level(i_atom_type, i_wild_card)
                else {
                    return (can_tor_type, None);
                };
                let Some(can_l_atom_type) =
                    mmff_def.torsion_equivalence_level(l_atom_type, l_wild_card)
                else {
                    return (can_tor_type, None);
                };
                let mut can_i_atom_type = u32::from(can_i_atom_type);
                let mut can_j_atom_type = j_atom_type;
                let mut can_k_atom_type = k_atom_type;
                let mut can_l_atom_type = u32::from(can_l_atom_type);
                if can_j_atom_type > can_k_atom_type {
                    std::mem::swap(&mut can_j_atom_type, &mut can_k_atom_type);
                    std::mem::swap(&mut can_i_atom_type, &mut can_l_atom_type);
                } else if (can_j_atom_type == can_k_atom_type)
                    && (can_i_atom_type > can_l_atom_type)
                {
                    std::mem::swap(&mut can_i_atom_type, &mut can_l_atom_type);
                }
                let key = (
                    can_tor_type as u8,
                    can_i_atom_type as u8,
                    can_j_atom_type as u8,
                    can_k_atom_type as u8,
                    can_l_atom_type as u8,
                );
                if let Ok(idx) =
                    params.binary_search_by_key(&key, |(tor, i, j, k, l, _)| (*tor, *i, *j, *k, *l))
                {
                    mmff_tor_params = Some(&params[idx].5);
                    if max_iter == 4 {
                        break;
                    }
                }
                iter += 1;
            }
            return (can_tor_type, mmff_tor_params);
        }

        let MmffTorStorage::Owned {
            params,
            tor_type: tor_types,
            i_atom_type: i_atom_types,
            j_atom_type: j_atom_types,
            k_atom_type: k_atom_types,
            l_atom_type: l_atom_types,
        } = &self.params
        else {
            unreachable!("static MMFFTor storage returned above")
        };

        let mut mmff_tor_params = None;
        let mut iter = 0_usize;
        let mut can_tor_type = tor_type.0;
        let mut max_iter = 5_usize;

        // RDKit✔️✔️: // For i-j-k-2 torsion interactions, a five-stage
        // RDKit✔️✔️: // process based on level combinations 1-1-1-1, 2-2-2-2,
        // RDKit✔️✔️: // 3-2-2-5, 5-2-2-3, and 5-2-2-5 is used, where stages 3
        // RDKit✔️✔️: // and 4 correspond to "half-default" or "half-wild-card" entries.
        // RDKit✔️✔️:     while (((iter < maxIter) && ((!mmffTorParams) || (maxIter == 4))) ||
        // RDKit✔️✔️:            ((iter == 4) && (torType.first == 5) && torType.second)) {
        while ((iter < max_iter) && (mmff_tor_params.is_none() || max_iter == 4))
            || ((iter == 4) && (tor_type.0 == 5) && (tor_type.1 != 0))
        {
            // RDKit✔️✔️:       // The rule of setting the torsion type to the value it had
            // RDKit✔️✔️:       // before being set to 5 as a last resort in case parameters
            // RDKit✔️✔️:       // could not be found is not mentioned in MMFF.IV; it was
            // RDKit✔️✔️:       // empirically discovered due to a number of tests in the
            // RDKit✔️✔️:       // MMFF validation suite otherwise failing
            // RDKit✔️✔️:       if ((maxIter == 5) && (iter == 4)) {
            // RDKit✔️✔️:         maxIter = 4;
            // RDKit✔️✔️:         iter = 0;
            // RDKit✔️✔️:         canTorType = torType.second;
            // RDKit✔️✔️:       }
            if (max_iter == 5) && (iter == 4) {
                max_iter = 4;
                iter = 0;
                can_tor_type = tor_type.1;
            }

            // RDKit✔️✔️:       iWildCard = iter;
            // RDKit✔️✔️:       lWildCard = iter;
            // RDKit✔️✔️:       if (iter == 1) {
            // RDKit✔️✔️:         iWildCard = 1;
            // RDKit✔️✔️:         lWildCard = 3;
            // RDKit✔️✔️:       } else if (iter == 2) {
            // RDKit✔️✔️:         iWildCard = 3;
            // RDKit✔️✔️:         lWildCard = 1;
            // RDKit✔️✔️:       }
            let mut i_wild_card = iter;
            let mut l_wild_card = iter;
            if iter == 1 {
                i_wild_card = 1;
                l_wild_card = 3;
            } else if iter == 2 {
                i_wild_card = 3;
                l_wild_card = 1;
            }

            // RDKit✔️✔️:       unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iWildCard];
            // RDKit✔️✔️:       unsigned int canJAtomType = jAtomType;
            // RDKit✔️✔️:       unsigned int canKAtomType = kAtomType;
            // RDKit✔️✔️:       unsigned int canLAtomType = (*mmffDef)(lAtomType)->eqLevel[lWildCard];
            let Some(can_i_atom_type) =
                mmff_def.torsion_equivalence_level(i_atom_type, i_wild_card)
            else {
                return (can_tor_type, None);
            };
            let Some(can_l_atom_type) =
                mmff_def.torsion_equivalence_level(l_atom_type, l_wild_card)
            else {
                return (can_tor_type, None);
            };
            let mut can_i_atom_type = u32::from(can_i_atom_type);
            let mut can_j_atom_type = j_atom_type;
            let mut can_k_atom_type = k_atom_type;
            let mut can_l_atom_type = u32::from(can_l_atom_type);

            // RDKit✔️✔️:       if (canJAtomType > canKAtomType) {
            // RDKit✔️✔️:         unsigned int temp = canKAtomType;
            // RDKit✔️✔️:         canKAtomType = canJAtomType;
            // RDKit✔️✔️:         canJAtomType = temp;
            // RDKit✔️✔️:         temp = canLAtomType;
            // RDKit✔️✔️:         canLAtomType = canIAtomType;
            // RDKit✔️✔️:         canIAtomType = temp;
            // RDKit✔️✔️:       } else if ((canJAtomType == canKAtomType) &&
            // RDKit✔️✔️:                  (canIAtomType > canLAtomType)) {
            // RDKit✔️✔️:         unsigned int temp = canLAtomType;
            // RDKit✔️✔️:         canLAtomType = canIAtomType;
            // RDKit✔️✔️:         canIAtomType = temp;
            // RDKit✔️✔️:       }
            if can_j_atom_type > can_k_atom_type {
                std::mem::swap(&mut can_j_atom_type, &mut can_k_atom_type);
                std::mem::swap(&mut can_i_atom_type, &mut can_l_atom_type);
            } else if (can_j_atom_type == can_k_atom_type) && (can_i_atom_type > can_l_atom_type) {
                std::mem::swap(&mut can_i_atom_type, &mut can_l_atom_type);
            }

            // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
            // RDKit❌❌:       const auto res1 = d_params.find(canTorType);
            // RDKit❌❌:       if (res1 != d_params.end()) {
            // RDKit❌❌:         const auto res2 = ((*res1).second).find(canIAtomType);
            // RDKit❌❌:         if (res2 != ((*res1).second).end()) {
            // RDKit❌❌:           const auto res3 = ((*res2).second).find(canJAtomType);
            // RDKit❌❌:           if (res3 != ((*res2).second).end()) {
            // RDKit❌❌:             const auto res4 = ((*res3).second).find(canKAtomType);
            // RDKit❌❌:             if (res4 != ((*res3).second).end()) {
            // RDKit❌❌:               const auto res5 = ((*res4).second).find(canLAtomType);
            // RDKit❌❌:               if (res5 != ((*res4).second).end()) {
            // RDKit❌❌:                 mmffTorParams = &((*res5).second);
            // RDKit❌❌:                 if (maxIter == 4) {
            // RDKit❌❌:                   break;
            // RDKit❌❌:                 }
            // RDKit❌❌:               }
            // RDKit❌❌:             }
            // RDKit❌❌:           }
            // RDKit❌❌:         }
            // RDKit❌❌:       }
            // RDKit✔️✔️: #else
            // RDKit✔️✔️:       auto jBounds = std::equal_range(d_jAtomType.begin(), d_jAtomType.end(),
            // RDKit✔️✔️:                                       canJAtomType);
            if let Some((j_begin, j_end)) =
                equal_range_u8(j_atom_types, 0, j_atom_types.len(), can_j_atom_type)
            {
                // RDKit✔️✔️:       if (jBounds.first != jBounds.second) {
                // RDKit✔️✔️:         auto bounds = std::equal_range(
                // RDKit✔️✔️:             d_kAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
                // RDKit✔️✔️:             d_kAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
                // RDKit✔️✔️:             canKAtomType);
                if let Some((k_begin, k_end)) =
                    equal_range_u8(k_atom_types, j_begin, j_end, can_k_atom_type)
                {
                    // RDKit✔️✔️:         if (bounds.first != bounds.second) {
                    // RDKit✔️✔️:           bounds = std::equal_range(
                    // RDKit✔️✔️:               d_iAtomType.begin() + (bounds.first - d_kAtomType.begin()),
                    // RDKit✔️✔️:               d_iAtomType.begin() + (bounds.second - d_kAtomType.begin()),
                    // RDKit✔️✔️:               canIAtomType);
                    if let Some((i_begin, i_end)) =
                        equal_range_u8(i_atom_types, k_begin, k_end, can_i_atom_type)
                    {
                        // RDKit✔️✔️:           if (bounds.first != bounds.second) {
                        // RDKit✔️✔️:             bounds = std::equal_range(
                        // RDKit✔️✔️:                 d_lAtomType.begin() + (bounds.first - d_iAtomType.begin()),
                        // RDKit✔️✔️:                 d_lAtomType.begin() + (bounds.second - d_iAtomType.begin()),
                        // RDKit✔️✔️:                 canLAtomType);
                        if let Some((l_begin, l_end)) =
                            equal_range_u8(l_atom_types, i_begin, i_end, can_l_atom_type)
                        {
                            // RDKit✔️✔️:             if (bounds.first != bounds.second) {
                            // RDKit✔️✔️:               bounds = std::equal_range(
                            // RDKit✔️✔️:                   d_torType.begin() + (bounds.first - d_lAtomType.begin()),
                            // RDKit✔️✔️:                   d_torType.begin() + (bounds.second - d_lAtomType.begin()),
                            // RDKit✔️✔️:                   canTorType);
                            if let Some((tor_begin, _tor_end)) =
                                equal_range_u8(tor_types, l_begin, l_end, can_tor_type)
                            {
                                // RDKit✔️✔️:               if (bounds.first != bounds.second) {
                                // RDKit✔️✔️:                 mmffTorParams = &d_params[bounds.first - d_torType.begin()];
                                mmff_tor_params = params.get(tor_begin);
                                // RDKit✔️✔️:                 if (maxIter == 4) {
                                // RDKit✔️✔️:                   break;
                                // RDKit✔️✔️:                 }
                                if max_iter == 4 {
                                    break;
                                }
                                // RDKit✔️✔️:               }
                            }
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️: #endif
            // RDKit✔️✔️:       ++iter;
            iter += 1;
        }

        // RDKit✔️✔️:     return std::make_pair(canTorType, mmffTorParams);
        // RDKit✔️✔️:   }
        (can_tor_type, mmff_tor_params)
    }
}

/// MMFF non-bonded Van der Waals parameter collection.
#[derive(Clone, Debug, PartialEq)]
pub struct MmffVdwCollection {
    source_mmff_vdw: Cow<'static, str>,
    power: f64,
    b: f64,
    beta: f64,
    darad: f64,
    daeps: f64,
    params: MmffVdwStorage,
}

#[derive(Clone, Debug, PartialEq)]
enum MmffVdwStorage {
    Static(&'static [(u8, MmffVdw)]),
    Owned {
        params: Vec<MmffVdw>,
        atom_type: Vec<u8>,
    },
}

impl MmffVdwCollection {
    pub fn new(mmff_vdw: &str) -> Result<Self, MmffParamError> {
        // RDKit✔️✔️: MMFFVdWCollection::MMFFVdWCollection(std::string mmffVdW) {
        // RDKit✔️✔️:   if (mmffVdW.empty()) {
        // RDKit✔️✔️:     mmffVdW = defaultMMFFVdW;
        // RDKit✔️✔️:   }
        let source_mmff_vdw = if mmff_vdw.is_empty() {
            Cow::Borrowed(default_mmff_vdw()?)
        } else {
            Cow::Owned(mmff_vdw.to_owned())
        };
        if mmff_vdw.is_empty() {
            Ok(Self {
                source_mmff_vdw,
                power: DEFAULT_MMFF_VDW_POWER,
                b: DEFAULT_MMFF_VDW_B,
                beta: DEFAULT_MMFF_VDW_BETA,
                darad: DEFAULT_MMFF_VDW_DARAD,
                daeps: DEFAULT_MMFF_VDW_DAEPS,
                params: MmffVdwStorage::Static(DEFAULT_MMFF_VDW_ROWS),
            })
        } else {
            let parsed = parse_mmff_vdw(source_mmff_vdw.as_ref())?;
            Ok(Self {
                source_mmff_vdw,
                power: parsed.power,
                b: parsed.b,
                beta: parsed.beta,
                darad: parsed.darad,
                daeps: parsed.daeps,
                params: MmffVdwStorage::Owned {
                    params: parsed.params,
                    atom_type: parsed.atom_type,
                },
            })
        }
    }

    pub fn source_mmff_vdw(&self) -> &str {
        self.source_mmff_vdw.as_ref()
    }

    pub fn power(&self) -> f64 {
        self.power
    }

    pub fn b(&self) -> f64 {
        self.b
    }

    pub fn beta(&self) -> f64 {
        self.beta
    }

    pub fn darad(&self) -> f64 {
        self.darad
    }

    pub fn daeps(&self) -> f64 {
        self.daeps
    }

    pub fn len(&self) -> usize {
        match &self.params {
            MmffVdwStorage::Static(params) => params.len(),
            MmffVdwStorage::Owned { params, .. } => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            MmffVdwStorage::Static(params) => params.is_empty(),
            MmffVdwStorage::Owned { params, .. } => params.is_empty(),
        }
    }

    pub fn get(&self, atom_type: u32) -> Option<&MmffVdw> {
        // RDKit✔️✔️:   const MMFFVdW *operator()(const unsigned int atomType) const {
        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res = d_params.find(atomType);
        // RDKit❌❌:     return (res != d_params.end() ? &((*res).second) : NULL);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto bounds =
        // RDKit✔️✔️:         std::equal_range(d_atomType.begin(), d_atomType.end(), atomType);
        match &self.params {
            MmffVdwStorage::Static(params) => params
                .binary_search_by_key(&(atom_type as u8), |(key, _)| *key)
                .ok()
                .map(|idx| &params[idx].1),
            MmffVdwStorage::Owned {
                params,
                atom_type: atom_types,
            } => {
                let (atom_begin, _atom_end) =
                    equal_range_u8(atom_types, 0, atom_types.len(), atom_type)?;

                // RDKit✔️✔️:     return ((bounds.first != bounds.second)
                // RDKit✔️✔️:                 ? &d_params[bounds.first - d_atomType.begin()]
                // RDKit✔️✔️:                 : nullptr);
                // RDKit✔️✔️: #endif
                // RDKit✔️✔️:   }
                params.get(atom_begin)
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct ParsedMmffChg {
    params: Vec<MmffChg>,
    bond_type: Vec<u8>,
    i_atom_type: Vec<u8>,
    j_atom_type: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq)]
struct ParsedMmffBond {
    params: Vec<MmffBond>,
    bond_type: Vec<u8>,
    i_atom_type: Vec<u8>,
    j_atom_type: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq)]
struct ParsedMmffAngle {
    params: Vec<MmffAngle>,
    angle_type: Vec<u8>,
    i_atom_type: Vec<u8>,
    j_atom_type: Vec<u8>,
    k_atom_type: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq)]
struct ParsedMmffStbn {
    params: Vec<MmffStbn>,
    stretch_bend_type: Vec<u8>,
    i_atom_type: Vec<u8>,
    j_atom_type: Vec<u8>,
    k_atom_type: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq)]
struct ParsedMmffOop {
    params: Vec<MmffOop>,
    i_atom_type: Vec<u8>,
    j_atom_type: Vec<u8>,
    k_atom_type: Vec<u8>,
    l_atom_type: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq)]
struct ParsedMmffTor {
    params: Vec<MmffTor>,
    tor_type: Vec<u8>,
    i_atom_type: Vec<u8>,
    j_atom_type: Vec<u8>,
    k_atom_type: Vec<u8>,
    l_atom_type: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq)]
struct ParsedMmffVdw {
    power: f64,
    b: f64,
    beta: f64,
    darad: f64,
    daeps: f64,
    params: Vec<MmffVdw>,
    atom_type: Vec<u8>,
}

pub fn default_mmff_def() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_DEF)
}

pub fn default_mmff_prop() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_PROP)
}

pub fn default_mmff_pbci() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_PBCI)
}

pub fn default_mmff_chg() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_CHG)
}

pub fn default_mmff_bond() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_BOND)
}

pub fn default_mmff_bndk() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_BNDK)
}

pub fn default_mmff_herschbach_laurie() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_HERSCHBACH_LAURIE)
}

pub fn default_mmff_cov_rad_pau_ele() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_COV_RAD_PAU_ELE)
}

pub(crate) fn default_mmff_bndk_params(
    atomic_num: u32,
    nbr_atomic_num: u32,
) -> Option<&'static MmffBond> {
    // RDKit❗✔️: const MMFFBond *operator()(const int atomicNum,
    // RDKit❗✔️:                              const int nbrAtomicNum) const {
    // RDKit❗✔️:   const MMFFBond *mmffBndkParams = nullptr;
    // RDKit❗✔️:   unsigned int canAtomicNum = atomicNum;
    // RDKit❗✔️:   unsigned int canNbrAtomicNum = nbrAtomicNum;
    // RDKit❗✔️:   if (atomicNum > nbrAtomicNum) {
    // RDKit❗✔️:     canAtomicNum = nbrAtomicNum;
    // RDKit❗✔️:     canNbrAtomicNum = atomicNum;
    // RDKit❗✔️:   }
    let key = if atomic_num <= nbr_atomic_num {
        (atomic_num as u8, nbr_atomic_num as u8)
    } else {
        (nbr_atomic_num as u8, atomic_num as u8)
    };
    // RDKit❗✔️:   auto bounds = std::equal_range(d_iAtomicNum.begin(), d_iAtomicNum.end(),
    // RDKit❗✔️:                                  canAtomicNum);
    // RDKit❗✔️:   if (bounds.first != bounds.second) {
    // RDKit❗✔️:     bounds = std::equal_range(
    // RDKit❗✔️:         d_jAtomicNum.begin() + (bounds.first - d_iAtomicNum.begin()),
    // RDKit❗✔️:         d_jAtomicNum.begin() + (bounds.second - d_iAtomicNum.begin()),
    // RDKit❗✔️:         canNbrAtomicNum);
    // RDKit❗✔️:     if (bounds.first != bounds.second) {
    // RDKit❗✔️:       mmffBndkParams = &d_params[bounds.first - d_jAtomicNum.begin()];
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return mmffBndkParams;
    // RDKit❗✔️: }
    DEFAULT_MMFF_BNDK_ROWS
        .binary_search_by_key(&key, |(i, j, _)| (*i, *j))
        .ok()
        .map(|idx| &DEFAULT_MMFF_BNDK_ROWS[idx].2)
}

pub(crate) fn default_mmff_herschbach_laurie_params(
    i_row: u32,
    j_row: u32,
) -> Option<&'static MmffHerschbachLaurie> {
    // RDKit❗✔️: const MMFFHerschbachLaurie *operator()(const int iRow, const int jRow) const {
    // RDKit❗✔️:   const MMFFHerschbachLaurie *mmffHerschbachLaurieParams = nullptr;
    // RDKit❗✔️:   unsigned int canIRow = iRow;
    // RDKit❗✔️:   unsigned int canJRow = jRow;
    // RDKit❗✔️:   if (iRow > jRow) {
    // RDKit❗✔️:     canIRow = jRow;
    // RDKit❗✔️:     canJRow = iRow;
    // RDKit❗✔️:   }
    let key = if i_row <= j_row {
        (i_row as u8, j_row as u8)
    } else {
        (j_row as u8, i_row as u8)
    };
    // RDKit❗✔️:   auto bounds = std::equal_range(d_iRow.begin(), d_iRow.end(), canIRow);
    // RDKit❗✔️:   if (bounds.first != bounds.second) {
    // RDKit❗✔️:     bounds = std::equal_range(
    // RDKit❗✔️:         d_jRow.begin() + (bounds.first - d_iRow.begin()),
    // RDKit❗✔️:         d_jRow.begin() + (bounds.second - d_iRow.begin()), canJRow);
    // RDKit❗✔️:     if (bounds.first != bounds.second) {
    // RDKit❗✔️:       mmffHerschbachLaurieParams = &d_params[bounds.first - d_jRow.begin()];
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return mmffHerschbachLaurieParams;
    // RDKit❗✔️: }
    DEFAULT_MMFF_HERSCHBACH_LAURIE_ROWS
        .binary_search_by_key(&key, |(i, j, _)| (*i, *j))
        .ok()
        .map(|idx| &DEFAULT_MMFF_HERSCHBACH_LAURIE_ROWS[idx].2)
}

pub(crate) fn default_mmff_cov_rad_pau_ele_params(
    atomic_num: u32,
) -> Option<&'static MmffCovRadPauEle> {
    // RDKit❗✔️: const MMFFCovRadPauEle *operator()(const unsigned int atomicNum) const {
    // RDKit❗✔️:   auto bounds =
    // RDKit❗✔️:       std::equal_range(d_atomicNum.begin(), d_atomicNum.end(), atomicNum);
    // RDKit❗✔️:   return ((bounds.first != bounds.second)
    // RDKit❗✔️:               ? &d_params[bounds.first - d_atomicNum.begin()]
    // RDKit❗✔️:               : nullptr);
    // RDKit❗✔️: }
    DEFAULT_MMFF_COV_RAD_PAU_ELE_ROWS
        .binary_search_by_key(&(atomic_num as u8), |(key, _)| *key)
        .ok()
        .map(|idx| &DEFAULT_MMFF_COV_RAD_PAU_ELE_ROWS[idx].1)
}

pub fn default_mmff_angle() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_ANGLE)
}

pub fn default_mmff_stbn() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_STBN)
}

pub fn default_mmff_dfsb() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_DFSB)
}

pub fn default_mmff_oop() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_OOP)
}

pub fn default_mmffs_oop() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFFS_OOP)
}

pub fn default_mmff_tor() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_TOR)
}

pub fn default_mmffs_tor() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFFS_TOR)
}

pub fn default_mmff_vdw() -> Result<&'static str, MmffParamError> {
    Ok(DEFAULT_MMFF_VDW)
}

fn parse_mmff_def(mmff_def: &str) -> Result<Vec<MmffDef>, MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffDef);
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   unsigned int oldAtomType = 0;
    // RDKit✔️✔️:   unsigned int atomType;
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut params = Vec::new();
    let mut old_atom_type = 0_u8;
    for (line_idx, line) in mmff_def.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 6 {
            return Err(MmffParamError::MalformedMmffDefLine {
                line_number,
                column_count: columns.len(),
            });
        }

        // RDKit✔️✔️:       // skip first token
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       atomType = (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        let atom_type = parse_u8(line_number, "TYPE", columns[1])?;

        // RDKit✔️✔️:       // Level 2 (currently = Level 1, see MMFF.I page 513)
        // RDKit✔️✔️:       mmffDefObj.eqLevel[0] =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       // Level 3
        // RDKit✔️✔️:       mmffDefObj.eqLevel[1] =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       // Level 4
        // RDKit✔️✔️:       mmffDefObj.eqLevel[2] =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       // Level 5
        // RDKit✔️✔️:       mmffDefObj.eqLevel[3] =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        let mmff_def_obj = MmffDef {
            eq_level: [
                parse_u8(line_number, "Level 2", columns[2])?,
                parse_u8(line_number, "Level 3", columns[3])?,
                parse_u8(line_number, "Level 4", columns[4])?,
                parse_u8(line_number, "Level 5", columns[5])?,
            ],
        };

        // RDKit✔️✔️:       if (atomType != oldAtomType) {
        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:         d_params[atomType] = mmffDefObj;
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:         d_params.push_back(mmffDefObj);
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:         oldAtomType = atomType;
        // RDKit✔️✔️:       }
        if atom_type != old_atom_type {
            params.push(mmff_def_obj);
            old_atom_type = atom_type;
        }
    }
    Ok(params)
}

fn parse_mmff_prop(mmff_prop: &str) -> Result<(Vec<MmffProp>, Vec<u8>), MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffProp);
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut params = Vec::new();
    let mut i_atom_type = Vec::new();
    for (line_idx, line) in mmff_prop.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 9 {
            return Err(MmffParamError::MalformedMmffPropLine {
                line_number,
                column_count: columns.len(),
            });
        }

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int atomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_iAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        i_atom_type.push(parse_u8(line_number, "atype", columns[0])?);

        // RDKit✔️✔️:       mmffPropObj.atno =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       mmffPropObj.crd =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       mmffPropObj.val =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       mmffPropObj.pilp =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       mmffPropObj.mltb =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       mmffPropObj.arom =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       mmffPropObj.linh =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       mmffPropObj.sbmb =
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️✔️:       ++token;
        let mmff_prop_obj = MmffProp {
            atno: parse_u8(line_number, "aspec", columns[1])?,
            crd: parse_u8(line_number, "crd", columns[2])?,
            val: parse_u8(line_number, "val", columns[3])?,
            pilp: parse_u8(line_number, "pilp", columns[4])?,
            mltb: parse_u8(line_number, "mltb", columns[5])?,
            arom: parse_u8(line_number, "arom", columns[6])?,
            linh: parse_u8(line_number, "lin", columns[7])?,
            sbmb: parse_u8(line_number, "sbmb", columns[8])?,
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[atomType] = mmffPropObj;
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_params.push_back(mmffPropObj);
        // RDKit✔️✔️: #endif
        params.push(mmff_prop_obj);
    }
    Ok((params, i_atom_type))
}

fn parse_mmff_chg(mmff_chg: &str) -> Result<ParsedMmffChg, MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffChg);
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut parsed = ParsedMmffChg {
        params: Vec::new(),
        bond_type: Vec::new(),
        i_atom_type: Vec::new(),
        j_atom_type: Vec::new(),
    };
    for (line_idx, line) in mmff_chg.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 4 {
            return Err(MmffParamError::MalformedMmffChgLine {
                line_number,
                column_count: columns.len(),
            });
        }

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int bondType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_bondType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .bond_type
            .push(parse_u8(line_number, "bondType", columns[0])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_iAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .i_atom_type
            .push(parse_u8(line_number, "iAtomType", columns[1])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_jAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .j_atom_type
            .push(parse_u8(line_number, "jAtomType", columns[2])?);

        let mmff_chg_obj = MmffChg {
            // RDKit✔️✔️:       mmffChgObj.bci = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            bci: parse_f64(line_number, "bci", columns[3])?,
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[bondType][iAtomType][jAtomType] = mmffChgObj;
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_params.push_back(mmffChgObj);
        // RDKit✔️✔️: #endif
        parsed.params.push(mmff_chg_obj);
    }
    Ok(parsed)
}

fn parse_mmff_bond(mmff_bond: &str) -> Result<ParsedMmffBond, MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffBond);
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut parsed = ParsedMmffBond {
        params: Vec::new(),
        bond_type: Vec::new(),
        i_atom_type: Vec::new(),
        j_atom_type: Vec::new(),
    };
    for (line_idx, line) in mmff_bond.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 5 {
            return Err(MmffParamError::MalformedMmffBondLine {
                line_number,
                column_count: columns.len(),
            });
        }

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int bondType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_bondType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .bond_type
            .push(parse_u8(line_number, "bondType", columns[0])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int atomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_iAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .i_atom_type
            .push(parse_u8(line_number, "atomType", columns[1])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_jAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .j_atom_type
            .push(parse_u8(line_number, "jAtomType", columns[2])?);

        let mmff_bond_obj = MmffBond {
            // RDKit✔️✔️:       mmffBondObj.kb = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            kb: parse_f64(line_number, "kb", columns[3])?,
            // RDKit✔️✔️:       mmffBondObj.r0 = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            r0: parse_f64(line_number, "r0", columns[4])?,
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[bondType][atomType][jAtomType] = mmffBondObj;
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_params.push_back(mmffBondObj);
        // RDKit✔️✔️: #endif
        parsed.params.push(mmff_bond_obj);
    }
    Ok(parsed)
}

fn parse_mmff_angle(mmff_angle: &str) -> Result<ParsedMmffAngle, MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffAngle);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut parsed = ParsedMmffAngle {
        params: Vec::new(),
        angle_type: Vec::new(),
        i_atom_type: Vec::new(),
        j_atom_type: Vec::new(),
        k_atom_type: Vec::new(),
    };
    for (line_idx, line) in mmff_angle.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 6 {
            return Err(MmffParamError::MalformedMmffAngleLine {
                line_number,
                column_count: columns.len(),
            });
        }

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int angleType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_angleType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .angle_type
            .push(parse_u8(line_number, "angleType", columns[0])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_iAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .i_atom_type
            .push(parse_u8(line_number, "iAtomType", columns[1])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_jAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .j_atom_type
            .push(parse_u8(line_number, "jAtomType", columns[2])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int kAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_kAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .k_atom_type
            .push(parse_u8(line_number, "kAtomType", columns[3])?);

        let mmff_angle_obj = MmffAngle {
            // RDKit✔️✔️:       mmffAngleObj.ka = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            ka: parse_f64(line_number, "ka", columns[4])?,
            // RDKit✔️✔️:       mmffAngleObj.theta0 = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            theta0: parse_f64(line_number, "theta0", columns[5])?,
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[angleType][iAtomType][jAtomType][kAtomType] = mmffAngleObj;
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_params.push_back(mmffAngleObj);
        // RDKit✔️✔️: #endif
        parsed.params.push(mmff_angle_obj);
    }
    Ok(parsed)
}

fn parse_mmff_stbn(mmff_stbn: &str) -> Result<ParsedMmffStbn, MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffStbn);
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut parsed = ParsedMmffStbn {
        params: Vec::new(),
        stretch_bend_type: Vec::new(),
        i_atom_type: Vec::new(),
        j_atom_type: Vec::new(),
        k_atom_type: Vec::new(),
    };
    for (line_idx, line) in mmff_stbn.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 6 {
            return Err(MmffParamError::MalformedMmffStbnLine {
                line_number,
                column_count: columns.len(),
            });
        }

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int stretchBendType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_stretchBendType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .stretch_bend_type
            .push(parse_u8(line_number, "stretchBendType", columns[0])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_iAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .i_atom_type
            .push(parse_u8(line_number, "iAtomType", columns[1])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_jAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .j_atom_type
            .push(parse_u8(line_number, "jAtomType", columns[2])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int kAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_kAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .k_atom_type
            .push(parse_u8(line_number, "kAtomType", columns[3])?);

        let mmff_stbn_obj = MmffStbn {
            // RDKit✔️✔️:       mmffStbnObj.kbaIJK = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            kba_ijk: parse_f64(line_number, "kbaIJK", columns[4])?,
            // RDKit✔️✔️:       mmffStbnObj.kbaKJI = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            kba_kji: parse_f64(line_number, "kbaKJI", columns[5])?,
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[stretchBendType][iAtomType][jAtomType][kAtomType] = mmffStbnObj;
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_params.push_back(mmffStbnObj);
        // RDKit✔️✔️: #endif
        parsed.params.push(mmff_stbn_obj);
    }
    Ok(parsed)
}

fn parse_mmff_dfsb(
    mmff_dfsb: &str,
) -> Result<BTreeMap<u32, BTreeMap<u32, BTreeMap<u32, MmffStbn>>>, MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffDfsb);
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut params = BTreeMap::new();
    for (line_idx, line) in mmff_dfsb.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 5 {
            return Err(MmffParamError::MalformedMmffDfsbLine {
                line_number,
                column_count: columns.len(),
            });
        }

        // RDKit✔️✔️:       auto iAtomicNum = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️:       ++token;
        let i_atomic_num = parse_u32(line_number, "iAtomicNum", columns[0])?;

        // RDKit✔️✔️:       auto jAtomicNum = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️:       ++token;
        let j_atomic_num = parse_u32(line_number, "jAtomicNum", columns[1])?;

        // RDKit✔️✔️:       auto kAtomicNum = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️:       ++token;
        let k_atomic_num = parse_u32(line_number, "kAtomicNum", columns[2])?;

        let mmff_stbn_obj = MmffStbn {
            // RDKit✔️✔️:       mmffStbnObj.kbaIJK = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            kba_ijk: parse_f64(line_number, "kbaIJK", columns[3])?,
            // RDKit✔️✔️:       mmffStbnObj.kbaKJI = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            kba_kji: parse_f64(line_number, "kbaKJI", columns[4])?,
        };

        // RDKit✔️✔️:       d_params[iAtomicNum][jAtomicNum][kAtomicNum] = mmffStbnObj;
        params
            .entry(i_atomic_num)
            .or_insert_with(BTreeMap::new)
            .entry(j_atomic_num)
            .or_insert_with(BTreeMap::new)
            .insert(k_atomic_num, mmff_stbn_obj);
    }
    Ok(params)
}

fn parse_mmff_oop(mmff_oop: &str) -> Result<ParsedMmffOop, MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffOop);
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut parsed = ParsedMmffOop {
        params: Vec::new(),
        i_atom_type: Vec::new(),
        j_atom_type: Vec::new(),
        k_atom_type: Vec::new(),
        l_atom_type: Vec::new(),
    };
    for (line_idx, line) in mmff_oop.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 5 {
            return Err(MmffParamError::MalformedMmffOopLine {
                line_number,
                column_count: columns.len(),
            });
        }

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_iAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .i_atom_type
            .push(parse_u8(line_number, "iAtomType", columns[0])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_jAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .j_atom_type
            .push(parse_u8(line_number, "jAtomType", columns[1])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int kAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_kAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .k_atom_type
            .push(parse_u8(line_number, "kAtomType", columns[2])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int lAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_lAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .l_atom_type
            .push(parse_u8(line_number, "lAtomType", columns[3])?);

        let mmff_oop_obj = MmffOop {
            // RDKit✔️✔️:       mmffOopObj.koop = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            koop: parse_f64(line_number, "koop", columns[4])?,
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[iAtomType][jAtomType][kAtomType][lAtomType] = mmffOopObj;
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_params.push_back(mmffOopObj);
        // RDKit✔️✔️: #endif
        parsed.params.push(mmff_oop_obj);
    }
    Ok(parsed)
}

fn parse_mmff_tor(mmff_tor: &str) -> Result<ParsedMmffTor, MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffTor);
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut parsed = ParsedMmffTor {
        params: Vec::new(),
        tor_type: Vec::new(),
        i_atom_type: Vec::new(),
        j_atom_type: Vec::new(),
        k_atom_type: Vec::new(),
        l_atom_type: Vec::new(),
    };
    for (line_idx, line) in mmff_tor.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 8 {
            return Err(MmffParamError::MalformedMmffTorLine {
                line_number,
                column_count: columns.len(),
            });
        }

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int torType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_torType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .tor_type
            .push(parse_u8(line_number, "torType", columns[0])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_iAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .i_atom_type
            .push(parse_u8(line_number, "iAtomType", columns[1])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_jAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .j_atom_type
            .push(parse_u8(line_number, "jAtomType", columns[2])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int kAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_kAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .k_atom_type
            .push(parse_u8(line_number, "kAtomType", columns[3])?);

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int lAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_lAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        parsed
            .l_atom_type
            .push(parse_u8(line_number, "lAtomType", columns[4])?);

        let mmff_tor_obj = MmffTor {
            // RDKit✔️✔️:       mmffTorObj.V1 = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            v1: parse_f64(line_number, "V1", columns[5])?,
            // RDKit✔️✔️:       mmffTorObj.V2 = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            v2: parse_f64(line_number, "V2", columns[6])?,
            // RDKit✔️✔️:       mmffTorObj.V3 = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            v3: parse_f64(line_number, "V3", columns[7])?,
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[torType][iAtomType][jAtomType][kAtomType][lAtomType] =
        // RDKit❌❌:           mmffTorObj;
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_params.push_back(mmffTorObj);
        // RDKit✔️✔️: #endif
        parsed.params.push(mmff_tor_obj);
    }
    Ok(parsed)
}

fn parse_mmff_vdw(mmff_vdw: &str) -> Result<ParsedMmffVdw, MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffVdW);
    // RDKit✔️✔️:   bool firstLine = true;
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut first_line = true;
    let mut parsed = ParsedMmffVdw {
        power: 0.0,
        b: 0.0,
        beta: 0.0,
        darad: 0.0,
        daeps: 0.0,
        params: Vec::new(),
        atom_type: Vec::new(),
    };

    for (line_idx, line) in mmff_vdw.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();

        // RDKit✔️✔️:       if (firstLine) {
        // RDKit✔️✔️:         firstLine = false;
        if first_line {
            first_line = false;
            if columns.len() < 5 {
                return Err(MmffParamError::MalformedMmffVdwConstantsLine {
                    line_number,
                    column_count: columns.len(),
                });
            }

            // RDKit✔️✔️:         this->power = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:         ++token;
            parsed.power = parse_f64(line_number, "power", columns[0])?;
            // RDKit✔️✔️:         this->B = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:         ++token;
            parsed.b = parse_f64(line_number, "B", columns[1])?;
            // RDKit✔️✔️:         this->Beta = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:         ++token;
            parsed.beta = parse_f64(line_number, "Beta", columns[2])?;
            // RDKit✔️✔️:         this->DARAD = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:         ++token;
            parsed.darad = parse_f64(line_number, "DARAD", columns[3])?;
            // RDKit✔️✔️:         this->DAEPS = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:         ++token;
            parsed.daeps = parse_f64(line_number, "DAEPS", columns[4])?;
        } else {
            // RDKit✔️✔️:       } else {
            // RDKit✔️✔️:         MMFFVdW mmffVdWObj;
            if columns.len() < 6 {
                return Err(MmffParamError::MalformedMmffVdwLine {
                    line_number,
                    column_count: columns.len(),
                });
            }

            // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
            // RDKit❌❌:         unsigned int atomType = boost::lexical_cast<unsigned int>(*token);
            // RDKit✔️✔️: #else
            // RDKit✔️✔️:         d_atomType.push_back(
            // RDKit✔️✔️:             (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
            // RDKit✔️✔️: #endif
            // RDKit✔️✔️:         ++token;
            parsed
                .atom_type
                .push(parse_u8(line_number, "atomType", columns[0])?);

            // RDKit✔️✔️:         mmffVdWObj.alpha_i = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:         ++token;
            let alpha_i = parse_f64(line_number, "alpha_i", columns[1])?;
            // RDKit✔️✔️:         mmffVdWObj.N_i = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:         ++token;
            let n_i = parse_f64(line_number, "N_i", columns[2])?;
            // RDKit✔️✔️:         mmffVdWObj.A_i = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:         ++token;
            let a_i = parse_f64(line_number, "A_i", columns[3])?;
            // RDKit✔️✔️:         mmffVdWObj.G_i = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:         ++token;
            let g_i = parse_f64(line_number, "G_i", columns[4])?;
            // RDKit✔️✔️:         mmffVdWObj.DA = (boost::lexical_cast<std::string>(*token)).at(0);
            // RDKit✔️✔️:         ++token;
            let da = *columns[5]
                .as_bytes()
                .first()
                .ok_or(MmffParamError::MalformedMmffVdwDa { line_number })?;

            // RDKit✔️✔️:         mmffVdWObj.R_star =
            // RDKit✔️✔️:             mmffVdWObj.A_i * pow(mmffVdWObj.alpha_i, this->power);
            let mmff_vdw_obj = MmffVdw {
                alpha_i,
                n_i,
                a_i,
                g_i,
                r_star: a_i * alpha_i.powf(parsed.power),
                da,
            };

            // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
            // RDKit❌❌:         d_params[atomType] = mmffVdWObj;
            // RDKit✔️✔️: #else
            // RDKit✔️✔️:         d_params.push_back(mmffVdWObj);
            // RDKit✔️✔️: #endif
            parsed.params.push(mmff_vdw_obj);
        }
    }

    Ok(parsed)
}

fn parse_mmff_pbci(mmff_pbci: &str) -> Result<Vec<MmffPbci>, MmffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(mmffPBCI);
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!(inStream.eof())) {
    // RDKit✔️✔️:     if (inLine[0] != '*') {
    let mut params = Vec::new();
    for (line_idx, line) in mmff_pbci.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 4 {
            return Err(MmffParamError::MalformedMmffPbciLine {
                line_number,
                column_count: columns.len(),
            });
        }

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       // IMPORTANT: skip the first field
        // RDKit❌❌:       ++token;
        // RDKit❌❌:       unsigned int atomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       // IMPORTANT: skip the first two fields
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        //
        // The vector branch ignores both leading fields and relies on row
        // order for lookup by atomType - 1.
        let mmff_pbci_obj = MmffPbci {
            // RDKit✔️✔️:       mmffPBCIObj.pbci = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            pbci: parse_f64(line_number, "pbci", columns[2])?,
            // RDKit✔️✔️:       mmffPBCIObj.fcadj = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       ++token;
            fcadj: parse_f64(line_number, "fcadj", columns[3])?,
        };

        // RDKit✔️✔️: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[atomType] = mmffPBCIObj;
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_params.push_back(mmffPBCIObj);
        // RDKit✔️✔️: #endif
        params.push(mmff_pbci_obj);
    }
    Ok(params)
}

fn equal_range_u8(values: &[u8], start: usize, end: usize, target: u32) -> Option<(usize, usize)> {
    let lower = start + values[start..end].partition_point(|&probe| u32::from(probe) < target);
    let upper = start + values[start..end].partition_point(|&probe| u32::from(probe) <= target);
    if lower != upper {
        Some((lower, upper))
    } else {
        None
    }
}

fn parse_u8(
    line_number: usize,
    column_name: &'static str,
    value: &str,
) -> Result<u8, MmffParamError> {
    let parsed = value
        .parse::<u32>()
        .map_err(|_err: ParseIntError| MmffParamError::ParseInt {
            line_number,
            column_name,
            value: value.to_owned(),
        })?;
    // RDKit✔️✔️:       atomType = (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
    // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
    Ok(parsed as u8)
}

fn parse_u32(
    line_number: usize,
    column_name: &'static str,
    value: &str,
) -> Result<u32, MmffParamError> {
    // RDKit✔️✔️:       auto iAtomicNum = boost::lexical_cast<unsigned int>(*token);
    value
        .parse::<u32>()
        .map_err(|_err: ParseIntError| MmffParamError::ParseInt {
            line_number,
            column_name,
            value: value.to_owned(),
        })
}

fn parse_f64(
    line_number: usize,
    column_name: &'static str,
    value: &str,
) -> Result<f64, MmffParamError> {
    value
        .parse::<f64>()
        .map_err(|_err: ParseFloatError| MmffParamError::ParseFloat {
            line_number,
            column_name,
            value: value.to_owned(),
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mmff_generated_bndk_lookup_matches_source_and_canonicalizes_order() {
        let forward = default_mmff_bndk_params(6, 8).expect("C-O Bndk row exists");
        let reversed = default_mmff_bndk_params(8, 6).expect("O-C Bndk row exists");

        assert_eq!(forward, &MmffBond { kb: 5.4, r0: 1.393 });
        assert_eq!(reversed, forward);
    }

    #[test]
    fn mmff_generated_herschbach_laurie_lookup_matches_source_and_canonicalizes_order() {
        let forward = default_mmff_herschbach_laurie_params(2, 4).expect("period-pair row exists");
        let reversed =
            default_mmff_herschbach_laurie_params(4, 2).expect("reversed period-pair row exists");

        assert_eq!(
            forward,
            &MmffHerschbachLaurie {
                a_ij: 2.61,
                d_ij: 1.28,
                dp_ij: 1.4,
            }
        );
        assert_eq!(reversed, forward);
    }

    #[test]
    fn mmff_generated_cov_rad_pau_ele_lookup_matches_source() {
        assert_eq!(
            default_mmff_cov_rad_pau_ele_params(8),
            Some(&MmffCovRadPauEle { r0: 0.72, chi: 3.5 })
        );
        assert_eq!(default_mmff_cov_rad_pau_ele_params(2), None);
    }

    #[test]
    fn mmff_mmffdefcollection_loads_default_definition_table() {
        let collection = MmffDefCollection::new("").expect("default MMFFDef parses");

        assert!(collection.source_mmff_def().contains("CR\t1\t1\t1\t1\t0"));
        assert_eq!(
            collection.get(1),
            Some(&MmffDef {
                eq_level: [1, 1, 1, 0]
            })
        );
        assert_eq!(
            collection.get(37),
            Some(&MmffDef {
                eq_level: [37, 2, 1, 0],
            })
        );
        assert_eq!(
            collection.get(82),
            Some(&MmffDef {
                eq_level: [82, 9, 8, 0],
            })
        );
    }

    #[test]
    fn mmff_mmffdefcollection_uses_source_vector_lookup_bounds() {
        let collection =
            MmffDefCollection::new("A\t3\t7\t8\t9\t10\n").expect("custom MMFFDef parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get(1),
            Some(&MmffDef {
                eq_level: [7, 8, 9, 10],
            })
        );
        assert_eq!(collection.get(0), None);
        assert_eq!(collection.get(2), None);
        assert_eq!(collection.get(3), None);
    }

    #[test]
    fn mmff_mmffdefcollection_skips_comment_alias_and_consecutive_duplicate_rows() {
        let collection = MmffDefCollection::new(
            "*\n\
             *\tALIAS\t1\t1\t1\t1\t0\n\
             A\t1\t1\t2\t3\t4\n\
             B\t1\t5\t6\t7\t8\n\
             C\t2\t9\t10\t11\t12\n",
        )
        .expect("custom MMFFDef parses");

        assert_eq!(collection.len(), 2);
        assert_eq!(
            collection.get(1),
            Some(&MmffDef {
                eq_level: [1, 2, 3, 4],
            })
        );
        assert_eq!(
            collection.get(2),
            Some(&MmffDef {
                eq_level: [9, 10, 11, 12],
            })
        );
    }

    #[test]
    fn mmff_mmffdefcollection_casts_unsigned_values_to_uint8_like_rdkit() {
        let collection =
            MmffDefCollection::new("A\t1\t300\t301\t302\t303\n").expect("custom MMFFDef parses");

        assert_eq!(
            collection.get(1),
            Some(&MmffDef {
                eq_level: [44, 45, 46, 47],
            })
        );
    }

    #[test]
    fn mmff_mmffdefcollection_rejects_malformed_line_with_too_few_columns() {
        let err = MmffDefCollection::new("A\t1\t2\t3\t4\n").expect_err("line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffDefLine {
                line_number: 1,
                column_count: 5,
            }
        );
    }

    #[test]
    fn mmff_mmffdefcollection_rejects_non_unsigned_integer_tokens() {
        let err = MmffDefCollection::new("A\tbad\t2\t3\t4\t5\n").expect_err("atom type is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseInt {
                line_number: 1,
                column_name: "TYPE",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffpropcollection_loads_default_property_table() {
        let collection = MmffPropCollection::new("").expect("default MMFFProp parses");

        assert!(
            collection
                .source_mmff_prop()
                .contains("* MMFFPROP - MMFF atom-type properties")
        );
        assert_eq!(
            collection.get(1),
            Some(&MmffProp {
                atno: 6,
                crd: 4,
                val: 4,
                pilp: 0,
                mltb: 0,
                arom: 0,
                linh: 0,
                sbmb: 0,
            })
        );
        assert_eq!(
            collection.get(37),
            Some(&MmffProp {
                atno: 6,
                crd: 3,
                val: 4,
                pilp: 0,
                mltb: 2,
                arom: 1,
                linh: 0,
                sbmb: 1,
            })
        );
        assert_eq!(
            collection.get(99),
            Some(&MmffProp {
                atno: 12,
                crd: 0,
                val: 0,
                pilp: 0,
                mltb: 0,
                arom: 0,
                linh: 0,
                sbmb: 0,
            })
        );
    }

    #[test]
    fn mmff_mmffpropcollection_uses_source_equal_range_lookup_bounds() {
        let collection = MmffPropCollection::new(
            "3\t6\t3\t4\t0\t2\t0\t0\t1\n\
             3\t7\t2\t3\t0\t2\t1\t0\t0\n\
             8\t8\t1\t2\t0\t2\t0\t0\t0\n",
        )
        .expect("custom MMFFProp parses");

        assert_eq!(collection.len(), 3);
        assert_eq!(
            collection.get(3),
            Some(&MmffProp {
                atno: 6,
                crd: 3,
                val: 4,
                pilp: 0,
                mltb: 2,
                arom: 0,
                linh: 0,
                sbmb: 1,
            })
        );
        assert_eq!(collection.get(0), None);
        assert_eq!(collection.get(1), None);
        assert_eq!(collection.get(9), None);
        assert_eq!(collection.get(300), None);
    }

    #[test]
    fn mmff_mmffpropcollection_skips_comment_rows() {
        let collection = MmffPropCollection::new(
            "*\n\
             * atype aspec crd val pilp mltb arom lin sbmb\n\
             1\t6\t4\t4\t0\t0\t0\t0\t0\n",
        )
        .expect("custom MMFFProp parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get(1),
            Some(&MmffProp {
                atno: 6,
                crd: 4,
                val: 4,
                pilp: 0,
                mltb: 0,
                arom: 0,
                linh: 0,
                sbmb: 0,
            })
        );
    }

    #[test]
    fn mmff_mmffpropcollection_casts_unsigned_values_to_uint8_like_rdkit() {
        let collection = MmffPropCollection::new("300\t301\t302\t303\t304\t305\t306\t307\t308\n")
            .expect("custom MMFFProp parses");

        assert_eq!(
            collection.get(44),
            Some(&MmffProp {
                atno: 45,
                crd: 46,
                val: 47,
                pilp: 48,
                mltb: 49,
                arom: 50,
                linh: 51,
                sbmb: 52,
            })
        );
        assert_eq!(collection.get(300), None);
    }

    #[test]
    fn mmff_mmffpropcollection_rejects_malformed_line_with_too_few_columns() {
        let err =
            MmffPropCollection::new("1\t6\t4\t4\t0\t0\t0\t0\n").expect_err("line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffPropLine {
                line_number: 1,
                column_count: 8,
            }
        );
    }

    #[test]
    fn mmff_mmffpropcollection_rejects_non_unsigned_integer_tokens() {
        let err = MmffPropCollection::new("1\tbad\t4\t4\t0\t0\t0\t0\t0\n")
            .expect_err("property value is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseInt {
                line_number: 1,
                column_name: "aspec",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffpbcicollection_loads_default_pbci_table() {
        let collection = MmffPbciCollection::new("").expect("default MMFFPBCI parses");

        assert!(
            collection
                .source_mmff_pbci()
                .contains("* MMFF Partial Bond Charge Incs")
        );
        assert_eq!(
            collection.get(1),
            Some(&MmffPbci {
                pbci: 0.0,
                fcadj: 0.0,
            })
        );
        assert_eq!(
            collection.get(32),
            Some(&MmffPbci {
                pbci: -0.732,
                fcadj: 0.5,
            })
        );
        assert_eq!(
            collection.get(99),
            Some(&MmffPbci {
                pbci: 2.0,
                fcadj: 0.0,
            })
        );
    }

    #[test]
    fn mmff_mmffpbcicollection_uses_source_vector_lookup_bounds() {
        let collection = MmffPbciCollection::new("9\t5\t1.250\t0.125\tignored\n")
            .expect("custom MMFFPBCI parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get(1),
            Some(&MmffPbci {
                pbci: 1.25,
                fcadj: 0.125,
            })
        );
        assert_eq!(collection.get(0), None);
        assert_eq!(collection.get(2), None);
        assert_eq!(collection.get(5), None);
    }

    #[test]
    fn mmff_mmffpbcicollection_skips_comment_rows() {
        let collection = MmffPbciCollection::new(
            "*\n\
             * type pbci fcadj\n\
             0\t1\t-0.135\t0.000\tE94\n",
        )
        .expect("custom MMFFPBCI parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get(1),
            Some(&MmffPbci {
                pbci: -0.135,
                fcadj: 0.0,
            })
        );
    }

    #[test]
    fn mmff_mmffpbcicollection_rejects_malformed_line_with_too_few_columns() {
        let err = MmffPbciCollection::new("0\t1\t-0.135\n").expect_err("line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffPbciLine {
                line_number: 1,
                column_count: 3,
            }
        );
    }

    #[test]
    fn mmff_mmffpbcicollection_rejects_non_float_tokens() {
        let err = MmffPbciCollection::new("0\t1\tbad\t0.000\n").expect_err("pbci is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseFloat {
                line_number: 1,
                column_name: "pbci",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffchgcollection_loads_default_charge_table() {
        let collection = MmffChgCollection::new("").expect("default MMFFChg parses");

        assert!(
            collection
                .source_mmff_chg()
                .contains("* MMFF BOND-CHARGE INCREMENTS")
        );
        let (sign, params) = collection.get_mmff_chg_params(0, 1, 2);
        assert_eq!(sign, -1);
        assert_eq!(params, Some(&MmffChg { bci: -0.1382 }));

        let (sign, params) = collection.get_mmff_chg_params(1, 2, 3);
        assert_eq!(sign, -1);
        assert_eq!(params, Some(&MmffChg { bci: -0.0144 }));
    }

    #[test]
    fn mmff_mmffchgcollection_canonicalizes_atom_order_and_flips_sign() {
        let collection = MmffChgCollection::new("").expect("default MMFFChg parses");

        let (sign, params) = collection.get_mmff_chg_params(0, 2, 1);

        assert_eq!(sign, 1);
        assert_eq!(params, Some(&MmffChg { bci: -0.1382 }));
    }

    #[test]
    fn mmff_mmffchgcollection_uses_source_nested_equal_range_bounds() {
        let collection = MmffChgCollection::new(
            "0\t1\t2\t-0.1250\tA\n\
             1\t1\t2\t-0.2500\tB\n\
             0\t1\t3\t-0.3750\tC\n\
             0\t2\t3\t-0.5000\tD\n",
        )
        .expect("custom MMFFChg parses");

        assert_eq!(collection.len(), 4);
        assert_eq!(
            collection.get_mmff_chg_params(1, 1, 2),
            (-1, Some(&MmffChg { bci: -0.25 }))
        );
        assert_eq!(collection.get_mmff_chg_params(2, 1, 2), (-1, None));
        assert_eq!(collection.get_mmff_chg_params(0, 1, 4), (-1, None));
        assert_eq!(collection.get_mmff_chg_params(0, 4, 1), (1, None));
    }

    #[test]
    fn mmff_mmffchgcollection_skips_comment_rows() {
        let collection = MmffChgCollection::new(
            "*\n\
             * types bci\n\
             0\t1\t2\t-0.1250\tA\n",
        )
        .expect("custom MMFFChg parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get_mmff_chg_params(0, 1, 2),
            (-1, Some(&MmffChg { bci: -0.125 }))
        );
    }

    #[test]
    fn mmff_mmffchgcollection_casts_unsigned_keys_to_uint8_like_rdkit() {
        let collection =
            MmffChgCollection::new("300\t301\t302\t-0.1250\n").expect("custom MMFFChg parses");

        assert_eq!(
            collection.get_mmff_chg_params(44, 45, 46),
            (-1, Some(&MmffChg { bci: -0.125 }))
        );
        assert_eq!(collection.get_mmff_chg_params(300, 301, 302), (-1, None));
    }

    #[test]
    fn mmff_mmffchgcollection_rejects_malformed_line_with_too_few_columns() {
        let err = MmffChgCollection::new("0\t1\t2\n").expect_err("line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffChgLine {
                line_number: 1,
                column_count: 3,
            }
        );
    }

    #[test]
    fn mmff_mmffchgcollection_rejects_non_unsigned_integer_tokens() {
        let err = MmffChgCollection::new("bad\t1\t2\t-0.1250\n").expect_err("bond type is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseInt {
                line_number: 1,
                column_name: "bondType",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffchgcollection_rejects_non_float_tokens() {
        let err = MmffChgCollection::new("0\t1\t2\tbad\n").expect_err("bci is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseFloat {
                line_number: 1,
                column_name: "bci",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffbondcollection_loads_default_bond_table() {
        let collection = MmffBondCollection::new("").expect("default MMFFBond parses");

        assert!(
            collection
                .source_mmff_bond()
                .contains("* MMFF BOND PARAMETERS")
        );
        assert_eq!(
            collection.get(0, 1, 1),
            Some(&MmffBond {
                kb: 4.258,
                r0: 1.508
            })
        );
        assert_eq!(
            collection.get(0, 1, 2),
            Some(&MmffBond {
                kb: 4.539,
                r0: 1.482
            })
        );
        assert_eq!(
            collection.get(1, 2, 3),
            Some(&MmffBond {
                kb: 4.565,
                r0: 1.468
            })
        );
    }

    #[test]
    fn mmff_mmffbondcollection_canonicalizes_atom_order_without_sign_change() {
        let collection = MmffBondCollection::new("").expect("default MMFFBond parses");

        assert_eq!(
            collection.get(0, 2, 1),
            Some(&MmffBond {
                kb: 4.539,
                r0: 1.482
            })
        );
    }

    #[test]
    fn mmff_mmffbondcollection_uses_source_nested_equal_range_bounds() {
        let collection = MmffBondCollection::new(
            "0\t1\t2\t1.250\t1.125\tA\n\
             1\t1\t2\t2.250\t2.125\tB\n\
             0\t1\t3\t3.250\t3.125\tC\n\
             0\t2\t3\t4.250\t4.125\tD\n",
        )
        .expect("custom MMFFBond parses");

        assert_eq!(collection.len(), 4);
        assert_eq!(
            collection.get(1, 1, 2),
            Some(&MmffBond {
                kb: 2.25,
                r0: 2.125
            })
        );
        assert_eq!(collection.get(2, 1, 2), None);
        assert_eq!(collection.get(0, 1, 4), None);
        assert_eq!(collection.get(0, 4, 1), None);
    }

    #[test]
    fn mmff_mmffbondcollection_skips_comment_rows() {
        let collection = MmffBondCollection::new(
            "*\n\
             * types kb r0 Source\n\
             0\t1\t2\t1.250\t1.125\tA\n",
        )
        .expect("custom MMFFBond parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get(0, 1, 2),
            Some(&MmffBond {
                kb: 1.25,
                r0: 1.125
            })
        );
    }

    #[test]
    fn mmff_mmffbondcollection_casts_unsigned_keys_to_uint8_like_rdkit() {
        let collection = MmffBondCollection::new("300\t301\t302\t1.250\t1.125\n")
            .expect("custom MMFFBond parses");

        assert_eq!(
            collection.get(44, 45, 46),
            Some(&MmffBond {
                kb: 1.25,
                r0: 1.125
            })
        );
        assert_eq!(collection.get(300, 301, 302), None);
    }

    #[test]
    fn mmff_mmffbondcollection_rejects_malformed_line_with_too_few_columns() {
        let err = MmffBondCollection::new("0\t1\t2\t1.250\n").expect_err("line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffBondLine {
                line_number: 1,
                column_count: 4,
            }
        );
    }

    #[test]
    fn mmff_mmffbondcollection_rejects_non_unsigned_integer_tokens() {
        let err =
            MmffBondCollection::new("bad\t1\t2\t1.250\t1.125\n").expect_err("bond type is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseInt {
                line_number: 1,
                column_name: "bondType",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffbondcollection_rejects_non_float_tokens() {
        let err = MmffBondCollection::new("0\t1\t2\tbad\t1.125\n").expect_err("kb is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseFloat {
                line_number: 1,
                column_name: "kb",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffanglecollection_loads_default_angle_array_table() {
        let defs = MmffDefCollection::new("").expect("default MMFFDef parses");
        let collection = MmffAngleCollection::new("").expect("default MMFFAngle parses");

        assert!(
            collection
                .source_mmff_angle()
                .contains("* MMFF ANGLE PARAMETERS")
        );
        assert_eq!(
            collection.get(&defs, 0, 1, 1, 1),
            Some(&MmffAngle {
                ka: 0.851,
                theta0: 109.608,
            })
        );
        assert_eq!(
            collection.get(&defs, 1, 1, 2, 3),
            Some(&MmffAngle {
                ka: 0.698,
                theta0: 116.104,
            })
        );
    }

    #[test]
    fn mmff_mmffanglecollection_canonicalizes_terminal_atom_types() {
        let defs = MmffDefCollection::new(
            "A\t1\t1\t5\t6\t7\n\
             B\t2\t2\t2\t2\t2\n\
             C\t3\t3\t8\t9\t10\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffAngleCollection::new("0\t1\t2\t3\t1.250\t120.500\tA\n")
            .expect("custom MMFFAngle parses");

        assert_eq!(
            collection.get(&defs, 0, 3, 2, 1),
            Some(&MmffAngle {
                ka: 1.25,
                theta0: 120.5,
            })
        );
    }

    #[test]
    fn mmff_mmffanglecollection_uses_source_level_fallback_and_nested_bounds() {
        let defs = MmffDefCollection::new(
            "A\t1\t1\t5\t6\t7\n\
             B\t2\t2\t2\t2\t2\n\
             C\t3\t3\t8\t9\t10\n\
             D\t4\t4\t4\t4\t4\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffAngleCollection::new(
            "1\t1\t2\t3\t2.250\t121.500\tB\n\
             0\t5\t2\t8\t3.250\t122.500\tC\n\
             0\t1\t4\t3\t4.250\t123.500\tD\n",
        )
        .expect("custom MMFFAngle parses");

        assert_eq!(collection.len(), 3);
        assert_eq!(
            collection.get(&defs, 1, 1, 2, 3),
            Some(&MmffAngle {
                ka: 2.25,
                theta0: 121.5,
            })
        );
        assert_eq!(
            collection.get(&defs, 0, 1, 2, 3),
            Some(&MmffAngle {
                ka: 3.25,
                theta0: 122.5,
            })
        );
        assert_eq!(collection.get(&defs, 2, 1, 2, 3), None);
        assert_eq!(collection.get(&defs, 0, 1, 9, 3), None);
        assert_eq!(collection.get(&defs, 0, 1, 2, 99), None);
    }

    #[test]
    fn mmff_mmffanglecollection_skips_comment_rows() {
        let defs = MmffDefCollection::new(
            "A\t1\t1\t1\t1\t1\n\
             B\t2\t2\t2\t2\t2\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffAngleCollection::new(
            "*\n\
             * angle types ka theta0\n\
             0\t1\t2\t1\t1.250\t120.500\tA\n",
        )
        .expect("custom MMFFAngle parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get(&defs, 0, 1, 2, 1),
            Some(&MmffAngle {
                ka: 1.25,
                theta0: 120.5,
            })
        );
    }

    #[test]
    fn mmff_mmffanglecollection_casts_unsigned_keys_to_uint8_like_rdkit() {
        let defs = MmffDefCollection::new(
            "A\t1\t300\t300\t300\t300\n\
             B\t2\t301\t301\t301\t301\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffAngleCollection::new("300\t44\t45\t44\t1.250\t120.500\n")
            .expect("custom MMFFAngle parses");

        assert_eq!(
            collection.get(&defs, 44, 1, 45, 1),
            Some(&MmffAngle {
                ka: 1.25,
                theta0: 120.5,
            })
        );
        assert_eq!(collection.get(&defs, 300, 1, 301, 1), None);
    }

    #[test]
    fn mmff_mmffanglecollection_rejects_malformed_line_with_too_few_columns() {
        let err = MmffAngleCollection::new("0\t1\t2\t3\t1.250\n").expect_err("line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffAngleLine {
                line_number: 1,
                column_count: 5,
            }
        );
    }

    #[test]
    fn mmff_mmffanglecollection_rejects_non_unsigned_integer_tokens() {
        let err = MmffAngleCollection::new("bad\t1\t2\t3\t1.250\t120.500\n")
            .expect_err("angle type is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseInt {
                line_number: 1,
                column_name: "angleType",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffanglecollection_rejects_non_float_tokens() {
        let err =
            MmffAngleCollection::new("0\t1\t2\t3\tbad\t120.500\n").expect_err("ka is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseFloat {
                line_number: 1,
                column_name: "ka",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffstbncollection_loads_default_stretch_bend_table() {
        let collection = MmffStbnCollection::new("").expect("default MMFFStbn parses");

        assert!(
            collection
                .source_mmff_stbn()
                .contains("* MMFF STRETCH-BEND PARAMETERS")
        );
        assert_eq!(
            collection.get_mmff_stbn_params(0, 0, 0, 1, 1, 1),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 0.206,
                    kba_kji: 0.206,
                }),
            )
        );
        assert_eq!(
            collection.get_mmff_stbn_params(0, 0, 0, 1, 1, 2),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 0.136,
                    kba_kji: 0.197,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmffstbncollection_canonicalizes_terminal_atom_types_and_reports_swap() {
        let collection = MmffStbnCollection::new("0\t1\t2\t3\t1.250\t1.125\tA\n")
            .expect("custom MMFFStbn parses");

        assert_eq!(
            collection.get_mmff_stbn_params(0, 0, 0, 3, 2, 1),
            (
                true,
                Some(&MmffStbn {
                    kba_ijk: 1.25,
                    kba_kji: 1.125,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmffstbncollection_sets_swap_for_equal_terminals_from_bond_order() {
        let collection = MmffStbnCollection::new("0\t1\t2\t1\t1.250\t1.125\tA\n")
            .expect("custom MMFFStbn parses");

        assert_eq!(
            collection.get_mmff_stbn_params(0, 1, 2, 1, 2, 1),
            (
                true,
                Some(&MmffStbn {
                    kba_ijk: 1.25,
                    kba_kji: 1.125,
                }),
            )
        );
        assert_eq!(
            collection.get_mmff_stbn_params(0, 2, 1, 1, 2, 1),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 1.25,
                    kba_kji: 1.125,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmffstbncollection_uses_source_nested_equal_range_bounds() {
        let collection = MmffStbnCollection::new(
            "0\t1\t2\t3\t1.250\t1.125\tA\n\
             1\t1\t2\t3\t2.250\t2.125\tB\n\
             0\t1\t2\t4\t3.250\t3.125\tC\n\
             0\t1\t5\t3\t4.250\t4.125\tD\n",
        )
        .expect("custom MMFFStbn parses");

        assert_eq!(collection.len(), 4);
        assert_eq!(
            collection.get_mmff_stbn_params(1, 0, 0, 1, 2, 3),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 2.25,
                    kba_kji: 2.125,
                }),
            )
        );
        assert_eq!(
            collection.get_mmff_stbn_params(2, 0, 0, 1, 2, 3),
            (false, None)
        );
        assert_eq!(
            collection.get_mmff_stbn_params(0, 0, 0, 1, 2, 9),
            (false, None)
        );
        assert_eq!(
            collection.get_mmff_stbn_params(0, 0, 0, 1, 9, 3),
            (false, None)
        );
        assert_eq!(
            collection.get_mmff_stbn_params(0, 0, 0, 9, 2, 1),
            (true, None)
        );
    }

    #[test]
    fn mmff_mmffstbncollection_skips_comment_rows() {
        let collection = MmffStbnCollection::new(
            "*\n\
             * types I J K kbaIJK kbaKJI\n\
             0\t1\t2\t3\t1.250\t1.125\tA\n",
        )
        .expect("custom MMFFStbn parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get_mmff_stbn_params(0, 0, 0, 1, 2, 3),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 1.25,
                    kba_kji: 1.125,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmffstbncollection_casts_unsigned_keys_to_uint8_like_rdkit() {
        let collection = MmffStbnCollection::new("300\t301\t302\t303\t1.250\t1.125\n")
            .expect("custom MMFFStbn parses");

        assert_eq!(
            collection.get_mmff_stbn_params(44, 0, 0, 45, 46, 47),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 1.25,
                    kba_kji: 1.125,
                }),
            )
        );
        assert_eq!(
            collection.get_mmff_stbn_params(300, 0, 0, 301, 302, 303),
            (false, None)
        );
    }

    #[test]
    fn mmff_mmffstbncollection_rejects_malformed_line_with_too_few_columns() {
        let err = MmffStbnCollection::new("0\t1\t2\t3\t1.250\n").expect_err("line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffStbnLine {
                line_number: 1,
                column_count: 5,
            }
        );
    }

    #[test]
    fn mmff_mmffstbncollection_rejects_non_unsigned_integer_tokens() {
        let err = MmffStbnCollection::new("bad\t1\t2\t3\t1.250\t1.125\n")
            .expect_err("stretch-bend type is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseInt {
                line_number: 1,
                column_name: "stretchBendType",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffstbncollection_rejects_non_float_tokens() {
        let err =
            MmffStbnCollection::new("0\t1\t2\t3\tbad\t1.125\n").expect_err("kbaIJK is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseFloat {
                line_number: 1,
                column_name: "kbaIJK",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffdfsbcollection_loads_default_stretch_bend_table() {
        let collection = MmffDfsbCollection::new("").expect("default MMFFDfsb parses");

        assert!(
            collection
                .source_mmff_dfsb()
                .contains("* MMFF Default Stretch-Bend Parameters")
        );
        assert_eq!(
            collection.get_mmff_dfsb_params(0, 1, 0),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 0.15,
                    kba_kji: 0.15,
                }),
            )
        );
        assert_eq!(
            collection.get_mmff_dfsb_params(4, 2, 4),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 0.25,
                    kba_kji: 0.25,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmffdfsbcollection_canonicalizes_terminal_rows_and_reports_swap() {
        let collection =
            MmffDfsbCollection::new("1\t2\t4\t1.250\t1.125\n").expect("custom MMFFDfsb parses");

        assert_eq!(
            collection.get_mmff_dfsb_params(4, 2, 1),
            (
                true,
                Some(&MmffStbn {
                    kba_ijk: 1.25,
                    kba_kji: 1.125,
                }),
            )
        );
        assert_eq!(
            collection.get_mmff_dfsb_params(1, 2, 4),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 1.25,
                    kba_kji: 1.125,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmffdfsbcollection_returns_none_for_missing_nested_map_branches() {
        let collection =
            MmffDfsbCollection::new("1\t2\t4\t1.250\t1.125\n").expect("custom MMFFDfsb parses");

        assert_eq!(collection.get_mmff_dfsb_params(0, 2, 4), (false, None));
        assert_eq!(collection.get_mmff_dfsb_params(1, 9, 4), (false, None));
        assert_eq!(collection.get_mmff_dfsb_params(1, 2, 9), (false, None));
        assert_eq!(collection.get_mmff_dfsb_params(9, 2, 1), (true, None));
    }

    #[test]
    fn mmff_mmffdfsbcollection_skips_comment_rows() {
        let collection = MmffDfsbCollection::new(
            "*\n\
             * rows I J K kbaIJK kbaKJI\n\
             1\t2\t4\t1.250\t1.125\n",
        )
        .expect("custom MMFFDfsb parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get_mmff_dfsb_params(1, 2, 4),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 1.25,
                    kba_kji: 1.125,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmffdfsbcollection_duplicate_keys_overwrite_like_std_map_assignment() {
        let collection = MmffDfsbCollection::new(
            "1\t2\t4\t1.250\t1.125\n\
             1\t2\t4\t2.250\t2.125\n",
        )
        .expect("custom MMFFDfsb parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get_mmff_dfsb_params(1, 2, 4),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 2.25,
                    kba_kji: 2.125,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmffdfsbcollection_preserves_unsigned_int_keys_without_uint8_cast() {
        let collection = MmffDfsbCollection::new("300\t301\t302\t1.250\t1.125\n")
            .expect("custom MMFFDfsb parses");

        assert_eq!(
            collection.get_mmff_dfsb_params(300, 301, 302),
            (
                false,
                Some(&MmffStbn {
                    kba_ijk: 1.25,
                    kba_kji: 1.125,
                }),
            )
        );
        assert_eq!(collection.get_mmff_dfsb_params(44, 301, 46), (false, None));
    }

    #[test]
    fn mmff_mmffdfsbcollection_rejects_malformed_line_with_too_few_columns() {
        let err = MmffDfsbCollection::new("1\t2\t4\t1.250\n").expect_err("line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffDfsbLine {
                line_number: 1,
                column_count: 4,
            }
        );
    }

    #[test]
    fn mmff_mmffdfsbcollection_rejects_non_unsigned_integer_tokens() {
        let err = MmffDfsbCollection::new("bad\t2\t4\t1.250\t1.125\n")
            .expect_err("iAtomicNum is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseInt {
                line_number: 1,
                column_name: "iAtomicNum",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffdfsbcollection_rejects_non_float_tokens() {
        let err = MmffDfsbCollection::new("1\t2\t4\tbad\t1.125\n").expect_err("kbaIJK is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseFloat {
                line_number: 1,
                column_name: "kbaIJK",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffoopcollection_loads_default_mmff_and_mmffs_tables() {
        let defs = MmffDefCollection::new("").expect("default MMFFDef parses");
        let mmff = MmffOopCollection::new(false, "").expect("default MMFFOop parses");
        let mmffs = MmffOopCollection::new(true, "").expect("default MMFFsOop parses");

        assert!(
            mmff.source_mmff_oop()
                .contains("* MMFF OUT-OF-PLANE PARAMETERS")
        );
        assert!(
            mmffs
                .source_mmff_oop()
                .contains("* MMFF94s OUT-OF-PLANE PARAMETERS")
        );
        assert_eq!(
            mmff.get_mmff_oop_params(&defs, 1, 10, 1, 1),
            Some(&MmffOop { koop: -0.02 })
        );
        assert_eq!(
            mmffs.get_mmff_oop_params(&defs, 1, 10, 1, 1),
            Some(&MmffOop { koop: 0.015 })
        );
    }

    #[test]
    fn mmff_mmffoopcollection_sorts_outer_atom_types_before_lookup() {
        let defs = MmffDefCollection::new(
            "A\t1\t1\t1\t1\t1\n\
             B\t2\t2\t2\t2\t2\n\
             C\t3\t3\t3\t3\t3\n\
             D\t4\t4\t4\t4\t4\n",
        )
        .expect("custom MMFFDef parses");
        let collection =
            MmffOopCollection::new(false, "1\t2\t3\t4\t1.250\n").expect("custom MMFFOop parses");

        assert_eq!(
            collection.get_mmff_oop_params(&defs, 4, 2, 1, 3),
            Some(&MmffOop { koop: 1.25 })
        );
    }

    #[test]
    fn mmff_mmffoopcollection_uses_source_level_fallback_and_nested_bounds() {
        let defs = MmffDefCollection::new(
            "A\t1\t9\t1\t1\t1\n\
             B\t2\t2\t2\t2\t2\n\
             C\t3\t8\t3\t3\t3\n\
             D\t4\t7\t4\t4\t4\n\
             E\t5\t6\t6\t6\t6\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffOopCollection::new(
            false,
            "1\t2\t3\t4\t1.250\n\
             5\t2\t3\t4\t2.250\n\
             1\t5\t3\t4\t3.250\n\
             1\t2\t5\t4\t4.250\n",
        )
        .expect("custom MMFFOop parses");

        assert_eq!(collection.len(), 4);
        assert_eq!(
            collection.get_mmff_oop_params(&defs, 1, 2, 3, 4),
            Some(&MmffOop { koop: 1.25 })
        );
        assert_eq!(collection.get_mmff_oop_params(&defs, 1, 9, 3, 4), None);
        assert_eq!(collection.get_mmff_oop_params(&defs, 4, 2, 3, 2), None);
        assert_eq!(collection.get_mmff_oop_params(&defs, 1, 2, 4, 5), None);
        assert_eq!(collection.get_mmff_oop_params(&defs, 1, 2, 3, 99), None);
    }

    #[test]
    fn mmff_mmffoopcollection_skips_comment_rows() {
        let defs = MmffDefCollection::new(
            "A\t1\t1\t1\t1\t1\n\
             B\t2\t2\t2\t2\t2\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffOopCollection::new(
            false,
            "*\n\
             * oop types koop\n\
             1\t2\t1\t1\t1.250\tA\n",
        )
        .expect("custom MMFFOop parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get_mmff_oop_params(&defs, 1, 2, 1, 1),
            Some(&MmffOop { koop: 1.25 })
        );
    }

    #[test]
    fn mmff_mmffoopcollection_casts_unsigned_keys_to_uint8_like_rdkit() {
        let defs = MmffDefCollection::new(
            "A\t1\t300\t300\t300\t300\n\
             B\t2\t301\t301\t301\t301\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffOopCollection::new(false, "44\t45\t44\t44\t1.250\n")
            .expect("custom MMFFOop parses");

        assert_eq!(
            collection.get_mmff_oop_params(&defs, 1, 45, 1, 1),
            Some(&MmffOop { koop: 1.25 })
        );
        assert_eq!(collection.get_mmff_oop_params(&defs, 1, 301, 1, 1), None);
    }

    #[test]
    fn mmff_mmffoopcollection_duplicate_vector_keys_return_first_match() {
        let defs = MmffDefCollection::new(
            "A\t1\t1\t1\t1\t1\n\
             B\t2\t2\t2\t2\t2\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffOopCollection::new(
            false,
            "1\t2\t1\t1\t1.250\n\
             1\t2\t1\t1\t2.250\n",
        )
        .expect("custom MMFFOop parses");

        assert_eq!(collection.len(), 2);
        assert_eq!(
            collection.get_mmff_oop_params(&defs, 1, 2, 1, 1),
            Some(&MmffOop { koop: 1.25 })
        );
    }

    #[test]
    fn mmff_mmffoopcollection_rejects_malformed_line_with_too_few_columns() {
        let err = MmffOopCollection::new(false, "1\t2\t1\t1\n").expect_err("line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffOopLine {
                line_number: 1,
                column_count: 4,
            }
        );
    }

    #[test]
    fn mmff_mmffoopcollection_rejects_non_unsigned_integer_tokens() {
        let err = MmffOopCollection::new(false, "bad\t2\t1\t1\t1.250\n")
            .expect_err("iAtomType is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseInt {
                line_number: 1,
                column_name: "iAtomType",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffoopcollection_rejects_non_float_tokens() {
        let err = MmffOopCollection::new(false, "1\t2\t1\t1\tbad\n").expect_err("koop is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseFloat {
                line_number: 1,
                column_name: "koop",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmfftorcollection_loads_default_mmff_and_mmffs_tables() {
        let defs = MmffDefCollection::new("").expect("default MMFFDef parses");
        let mmff = MmffTorCollection::new(false, "").expect("default MMFFTor parses");
        let mmffs = MmffTorCollection::new(true, "").expect("default MMFFsTor parses");

        assert!(mmff.source_mmff_tor().contains("* MMFF TORSION PARAMETERS"));
        assert!(
            mmffs
                .source_mmff_tor()
                .contains("* MMFF94s TORSION PARAMETERS")
        );
        assert_eq!(
            mmff.get_mmff_tor_params(&defs, (0, 0), 1, 1, 1, 1),
            (
                0,
                Some(&MmffTor {
                    v1: 0.103,
                    v2: 0.681,
                    v3: 0.332,
                }),
            )
        );
        assert_eq!(
            mmff.get_mmff_tor_params(&defs, (0, 0), 5, 1, 1, 10),
            (
                0,
                Some(&MmffTor {
                    v1: 0.0,
                    v2: 0.0,
                    v3: 0.427,
                }),
            )
        );
        assert_eq!(
            mmffs.get_mmff_tor_params(&defs, (0, 0), 5, 1, 1, 10),
            (
                0,
                Some(&MmffTor {
                    v1: 0.0,
                    v2: 0.0,
                    v3: 0.418,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmfftorcollection_canonicalizes_central_and_terminal_atom_types() {
        let defs = MmffDefCollection::new(
            "A\t1\t1\t1\t1\t1\n\
             B\t2\t2\t2\t2\t2\n\
             C\t3\t3\t3\t3\t3\n\
             D\t4\t4\t4\t4\t4\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffTorCollection::new(false, "0\t4\t2\t3\t1\t1.250\t1.125\t1.000\n")
            .expect("custom MMFFTor parses");

        assert_eq!(
            collection.get_mmff_tor_params(&defs, (0, 0), 1, 3, 2, 4),
            (
                0,
                Some(&MmffTor {
                    v1: 1.25,
                    v2: 1.125,
                    v3: 1.0,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmfftorcollection_canonicalizes_terminal_order_when_central_types_equal() {
        let defs = MmffDefCollection::new(
            "A\t1\t4\t4\t4\t4\n\
             B\t2\t2\t2\t2\t2\n\
             C\t3\t1\t1\t1\t1\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffTorCollection::new(false, "0\t1\t2\t2\t4\t1.250\t1.125\t1.000\n")
            .expect("custom MMFFTor parses");

        assert_eq!(
            collection.get_mmff_tor_params(&defs, (0, 0), 1, 2, 2, 3),
            (
                0,
                Some(&MmffTor {
                    v1: 1.25,
                    v2: 1.125,
                    v3: 1.0,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmfftorcollection_uses_source_fallback_stages_and_nested_bounds() {
        let defs = MmffDefCollection::new(
            "A\t1\t5\t6\t7\t8\n\
             B\t2\t2\t2\t2\t2\n\
             C\t3\t3\t3\t3\t3\n\
             D\t4\t9\t10\t11\t12\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffTorCollection::new(
            false,
            "0\t5\t2\t3\t9\t1.250\t1.125\t1.000\n\
             1\t5\t2\t3\t9\t2.250\t2.125\t2.000\n\
             0\t6\t2\t3\t10\t3.250\t3.125\t3.000\n\
             0\t5\t2\t4\t9\t4.250\t4.125\t4.000\n",
        )
        .expect("custom MMFFTor parses");

        assert_eq!(collection.len(), 4);
        assert_eq!(
            collection.get_mmff_tor_params(&defs, (0, 0), 1, 2, 3, 4),
            (
                0,
                Some(&MmffTor {
                    v1: 1.25,
                    v2: 1.125,
                    v3: 1.0,
                }),
            )
        );
        assert_eq!(
            collection.get_mmff_tor_params(&defs, (1, 0), 1, 2, 3, 4),
            (
                1,
                Some(&MmffTor {
                    v1: 2.25,
                    v2: 2.125,
                    v3: 2.0,
                }),
            )
        );
        assert_eq!(
            collection.get_mmff_tor_params(&defs, (0, 0), 1, 2, 3, 99),
            (0, None)
        );
        assert_eq!(
            collection.get_mmff_tor_params(&defs, (0, 0), 1, 2, 99, 4),
            (0, None)
        );
        assert_eq!(
            collection.get_mmff_tor_params(&defs, (0, 0), 99, 2, 3, 4),
            (0, None)
        );
        assert_eq!(
            collection.get_mmff_tor_params(&defs, (0, 0), 1, 99, 3, 4),
            (0, None)
        );
        assert_eq!(
            collection.get_mmff_tor_params(&defs, (2, 3), 1, 2, 3, 4),
            (3, None)
        );
    }

    #[test]
    fn mmff_mmfftorcollection_resets_type_five_search_to_secondary_type() {
        let defs = MmffDefCollection::new(
            "A\t1\t5\t6\t7\t8\n\
             B\t2\t2\t2\t2\t2\n\
             C\t3\t3\t3\t3\t3\n\
             D\t4\t9\t10\t11\t12\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffTorCollection::new(false, "2\t5\t2\t3\t9\t1.250\t1.125\t1.000\n")
            .expect("custom MMFFTor parses");

        assert_eq!(
            collection.get_mmff_tor_params(&defs, (5, 2), 1, 2, 3, 4),
            (
                2,
                Some(&MmffTor {
                    v1: 1.25,
                    v2: 1.125,
                    v3: 1.0,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmfftorcollection_uses_source_final_wildcard_stage_after_type_five_fallback() {
        let defs = MmffDefCollection::new(
            "A\t1\t5\t6\t7\t8\n\
             B\t2\t2\t2\t2\t2\n\
             C\t3\t3\t3\t3\t3\n\
             D\t4\t9\t10\t11\t12\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffTorCollection::new(false, "2\t0\t2\t3\t0\t1.250\t1.125\t1.000\n")
            .expect("custom MMFFTor parses");

        assert_eq!(
            collection.get_mmff_tor_params(&defs, (5, 2), 1, 2, 3, 4),
            (
                2,
                Some(&MmffTor {
                    v1: 1.25,
                    v2: 1.125,
                    v3: 1.0,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmfftorcollection_final_wildcard_stage_matches_pinned_rdkit_parameters() {
        let defs = MmffDefCollection::new("").expect("default MMFFDef parses");
        let collection = MmffTorCollection::new(false, "").expect("default MMFFTor parses");

        assert_eq!(
            collection.get_mmff_tor_params(&defs, (5, 2), 2, 3, 17, 1),
            (
                2,
                Some(&MmffTor {
                    v1: 0.0,
                    v2: 1.423,
                    v3: 0.0,
                }),
            )
        );
        assert_eq!(
            collection.get_mmff_tor_params(&defs, (5, 1), 64, 54, 2, 1),
            (1, None)
        );
    }

    #[test]
    fn mmff_mmfftorcollection_skips_comment_rows() {
        let defs = MmffDefCollection::new(
            "A\t1\t1\t1\t1\t1\n\
             B\t2\t2\t2\t2\t2\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffTorCollection::new(
            false,
            "*\n\
             * torsion types V1 V2 V3\n\
             0\t1\t2\t2\t1\t1.250\t1.125\t1.000\tA\n",
        )
        .expect("custom MMFFTor parses");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get_mmff_tor_params(&defs, (0, 0), 1, 2, 2, 1),
            (
                0,
                Some(&MmffTor {
                    v1: 1.25,
                    v2: 1.125,
                    v3: 1.0,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmfftorcollection_casts_unsigned_keys_to_uint8_like_rdkit() {
        let defs = MmffDefCollection::new(
            "A\t1\t300\t300\t300\t300\n\
             B\t2\t301\t301\t301\t301\n",
        )
        .expect("custom MMFFDef parses");
        let collection =
            MmffTorCollection::new(false, "300\t44\t45\t45\t44\t1.250\t1.125\t1.000\n")
                .expect("custom MMFFTor parses");

        assert_eq!(
            collection.get_mmff_tor_params(&defs, (44, 0), 1, 45, 45, 1),
            (
                44,
                Some(&MmffTor {
                    v1: 1.25,
                    v2: 1.125,
                    v3: 1.0,
                }),
            )
        );
        assert_eq!(
            collection.get_mmff_tor_params(&defs, (300, 0), 1, 301, 301, 1),
            (0, None)
        );
    }

    #[test]
    fn mmff_mmfftorcollection_duplicate_vector_keys_return_first_match() {
        let defs = MmffDefCollection::new(
            "A\t1\t1\t1\t1\t1\n\
             B\t2\t2\t2\t2\t2\n",
        )
        .expect("custom MMFFDef parses");
        let collection = MmffTorCollection::new(
            false,
            "0\t1\t2\t2\t1\t1.250\t1.125\t1.000\n\
             0\t1\t2\t2\t1\t2.250\t2.125\t2.000\n",
        )
        .expect("custom MMFFTor parses");

        assert_eq!(
            collection.get_mmff_tor_params(&defs, (0, 0), 1, 2, 2, 1),
            (
                0,
                Some(&MmffTor {
                    v1: 1.25,
                    v2: 1.125,
                    v3: 1.0,
                }),
            )
        );
    }

    #[test]
    fn mmff_mmfftorcollection_rejects_malformed_line_with_too_few_columns() {
        let err = MmffTorCollection::new(false, "0\t1\t2\t2\t1\t1.250\t1.125\n")
            .expect_err("line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffTorLine {
                line_number: 1,
                column_count: 7,
            }
        );
    }

    #[test]
    fn mmff_mmfftorcollection_rejects_non_unsigned_integer_tokens() {
        let err = MmffTorCollection::new(false, "bad\t1\t2\t2\t1\t1.250\t1.125\t1.000\n")
            .expect_err("torsion type is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseInt {
                line_number: 1,
                column_name: "torType",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmfftorcollection_rejects_non_float_tokens() {
        let err = MmffTorCollection::new(false, "0\t1\t2\t2\t1\tbad\t1.125\t1.000\n")
            .expect_err("V1 is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseFloat {
                line_number: 1,
                column_name: "V1",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffvdwcollection_loads_default_table_and_constants() {
        let collection = MmffVdwCollection::new("").expect("default MMFFVdW parses");

        assert!(
            collection
                .source_mmff_vdw()
                .contains("*  power      B       Beta     DARAD      DAEPS")
        );
        assert_eq!(collection.power(), 0.25);
        assert_eq!(collection.b(), 0.2);
        assert_eq!(collection.beta(), 12.0);
        assert_eq!(collection.darad(), 0.8);
        assert_eq!(collection.daeps(), 0.5);
        assert_eq!(
            collection.get(1),
            Some(&MmffVdw {
                alpha_i: 1.05,
                n_i: 2.49,
                a_i: 3.89,
                g_i: 1.282,
                r_star: 3.89_f64 * 1.05_f64.powf(0.25),
                da: b'-',
            })
        );
        assert_eq!(
            collection.get(6),
            Some(&MmffVdw {
                alpha_i: 0.70,
                n_i: 3.15,
                a_i: 3.89,
                g_i: 1.282,
                r_star: 3.89_f64 * 0.70_f64.powf(0.25),
                da: b'A',
            })
        );
        assert_eq!(
            collection.get(21),
            Some(&MmffVdw {
                alpha_i: 0.150,
                n_i: 0.800,
                a_i: 4.200,
                g_i: 1.209,
                r_star: 4.200_f64 * 0.150_f64.powf(0.25),
                da: b'D',
            })
        );
    }

    #[test]
    fn mmff_mmffvdwcollection_uses_source_equal_range_lookup() {
        let collection = MmffVdwCollection::new(
            "0.5\t0.2\t12.0\t0.8\t0.5\n\
             1\t4.000\t2.000\t3.000\t1.000\t-\tA\n\
             2\t9.000\t3.000\t4.000\t1.500\tA\tB\n\
             4\t16.000\t4.000\t5.000\t2.000\tD\tC\n",
        )
        .expect("custom MMFFVdW parses");

        assert_eq!(collection.len(), 3);
        assert_eq!(
            collection.get(2),
            Some(&MmffVdw {
                alpha_i: 9.0,
                n_i: 3.0,
                a_i: 4.0,
                g_i: 1.5,
                r_star: 4.0_f64 * 9.0_f64.powf(0.5),
                da: b'A',
            })
        );
        assert_eq!(collection.get(3), None);
        assert_eq!(collection.get(0), None);
    }

    #[test]
    fn mmff_mmffvdwcollection_skips_comment_rows_before_constants_and_params() {
        let collection = MmffVdwCollection::new(
            "*\n\
             * constants\n\
             0.5\t0.2\t12.0\t0.8\t0.5\n\
             * params\n\
             1\t4.000\t2.000\t3.000\t1.000\t-\tA\n",
        )
        .expect("custom MMFFVdW parses");

        assert_eq!(collection.power(), 0.5);
        assert_eq!(collection.len(), 1);
        assert_eq!(
            collection.get(1),
            Some(&MmffVdw {
                alpha_i: 4.0,
                n_i: 2.0,
                a_i: 3.0,
                g_i: 1.0,
                r_star: 3.0_f64 * 4.0_f64.powf(0.5),
                da: b'-',
            })
        );
    }

    #[test]
    fn mmff_mmffvdwcollection_casts_unsigned_keys_to_uint8_like_rdkit() {
        let collection = MmffVdwCollection::new(
            "0.5\t0.2\t12.0\t0.8\t0.5\n\
             300\t4.000\t2.000\t3.000\t1.000\t-\tA\n",
        )
        .expect("custom MMFFVdW parses");

        assert_eq!(
            collection.get(44),
            Some(&MmffVdw {
                alpha_i: 4.0,
                n_i: 2.0,
                a_i: 3.0,
                g_i: 1.0,
                r_star: 3.0_f64 * 4.0_f64.powf(0.5),
                da: b'-',
            })
        );
        assert_eq!(collection.get(300), None);
    }

    #[test]
    fn mmff_mmffvdwcollection_duplicate_vector_keys_return_first_match() {
        let collection = MmffVdwCollection::new(
            "0.5\t0.2\t12.0\t0.8\t0.5\n\
             1\t4.000\t2.000\t3.000\t1.000\t-\tA\n\
             1\t9.000\t3.000\t4.000\t1.500\tA\tB\n",
        )
        .expect("custom MMFFVdW parses");

        assert_eq!(
            collection.get(1),
            Some(&MmffVdw {
                alpha_i: 4.0,
                n_i: 2.0,
                a_i: 3.0,
                g_i: 1.0,
                r_star: 3.0_f64 * 4.0_f64.powf(0.5),
                da: b'-',
            })
        );
    }

    #[test]
    fn mmff_mmffvdwcollection_rejects_malformed_constants_line() {
        let err = MmffVdwCollection::new("0.5\t0.2\t12.0\t0.8\n")
            .expect_err("constants line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffVdwConstantsLine {
                line_number: 1,
                column_count: 4,
            }
        );
    }

    #[test]
    fn mmff_mmffvdwcollection_rejects_malformed_parameter_line() {
        let err = MmffVdwCollection::new(
            "0.5\t0.2\t12.0\t0.8\t0.5\n\
             1\t4.000\t2.000\t3.000\t1.000\n",
        )
        .expect_err("parameter line is malformed");

        assert_eq!(
            err,
            MmffParamError::MalformedMmffVdwLine {
                line_number: 2,
                column_count: 5,
            }
        );
    }

    #[test]
    fn mmff_mmffvdwcollection_rejects_empty_da_token() {
        let err = MmffVdwCollection::new(
            "0.5\t0.2\t12.0\t0.8\t0.5\n\
             1\t4.000\t2.000\t3.000\t1.000\t\n",
        )
        .expect_err("DA token is malformed");

        assert_eq!(err, MmffParamError::MalformedMmffVdwDa { line_number: 2 });
    }

    #[test]
    fn mmff_mmffvdwcollection_rejects_non_unsigned_integer_tokens() {
        let err = MmffVdwCollection::new(
            "0.5\t0.2\t12.0\t0.8\t0.5\n\
             bad\t4.000\t2.000\t3.000\t1.000\t-\n",
        )
        .expect_err("atom type is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseInt {
                line_number: 2,
                column_name: "atomType",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffvdwcollection_rejects_non_float_constants() {
        let err =
            MmffVdwCollection::new("bad\t0.2\t12.0\t0.8\t0.5\n").expect_err("power is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseFloat {
                line_number: 1,
                column_name: "power",
                value: "bad".to_owned(),
            }
        );
    }

    #[test]
    fn mmff_mmffvdwcollection_rejects_non_float_parameters() {
        let err = MmffVdwCollection::new(
            "0.5\t0.2\t12.0\t0.8\t0.5\n\
             1\tbad\t2.000\t3.000\t1.000\t-\n",
        )
        .expect_err("alpha_i is invalid");

        assert_eq!(
            err,
            MmffParamError::ParseFloat {
                line_number: 2,
                column_name: "alpha_i",
                value: "bad".to_owned(),
            }
        );
    }
}
