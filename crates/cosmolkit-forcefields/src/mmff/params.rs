use super::params_text::{
    MmffParamParseCause, MmffParamParseError, MmffParamTable, parse_mmff_f64, parse_mmff_u32,
    source_lines, tokenize_mmff_line,
};
use std::collections::BTreeMap;
use std::sync::OnceLock;

#[cfg(test)]
use std::sync::atomic::{AtomicUsize, Ordering};

#[cfg(test)]
static DEFAULT_MMFF_DEF_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_PROP_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_PBCI_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_ANGLE_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_OOP_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_OOP_S_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_TOR_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_TOR_S_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_VDW_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_CHG_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_STBN_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_DFSB_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_HERSCHBACH_LAURIE_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);
#[cfg(test)]
static DEFAULT_MMFF_COV_RAD_PAU_ELE_CONSTRUCTIONS: AtomicUsize = AtomicUsize::new(0);

const DEFAULT_MMFF_AROMATIC_TYPES: [u8; 17] = [
    37, 38, 39, 44, 58, 59, 63, 64, 65, 66, 69, 76, 78, 79, 80, 81, 82,
];

struct MMFFAromCollection {
    d_params: Vec<u8>,
}

impl MMFFAromCollection {
    fn new(mmff_arom: Option<&[u8]>) -> Self {
        // RDKit source, Params.cpp:23-25 and 29-41:
        // RDKit✔️✔️: const std::vector<std::uint8_t> defaultMMFFArom = {
        // RDKit✔️✔️:     37, 38, 39, 44, 58, 59, 63, 64, 65, 66, 69, 76, 78, 79, 80, 81, 82};
        // RDKit✔️✔️: MMFFAromCollection::MMFFAromCollection(
        // RDKit✔️✔️:     const std::vector<std::uint8_t> *mmffArom) {
        // RDKit✔️✔️:   if (!mmffArom) {
        // RDKit✔️✔️:     mmffArom = &defaultMMFFArom;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   d_params.clear();
        // RDKit✔️✔️:   d_params.resize(mmffArom->size());
        // RDKit✔️✔️:   std::copy(mmffArom->begin(), mmffArom->end(), d_params.begin());
        // RDKit✔️✔️: }

        let source = mmff_arom.unwrap_or(&DEFAULT_MMFF_AROMATIC_TYPES);
        let mut d_params = Vec::new();
        d_params.clear();
        d_params.resize(source.len(), 0);
        d_params.copy_from_slice(source);
        Self { d_params }
    }

    fn is_mmff_aromatic(&self, atom_type: u32) -> bool {
        // RDKit source, Params.h:156-159:
        // RDKit✔️✔️:   bool isMMFFAromatic(const unsigned int atomType) const {
        // RDKit✔️✔️:     return std::find(d_params.begin(), d_params.end(), atomType) !=
        // RDKit✔️✔️:            d_params.end();
        // RDKit✔️✔️:   }
        self.d_params
            .iter()
            .any(|&stored_type| u32::from(stored_type) == atom_type)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct MmffDef {
    pub(super) eq_level: [u8; 4],
}

#[derive(Debug)]
pub(super) struct MmffDefCollection {
    d_params: Vec<MmffDef>,
}

impl MmffDefCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:45-94:
        // RDKit✔️🔝: MMFFDefCollection::MMFFDefCollection(std::string mmffDef) {
        // RDKit✔️🔝:   if (mmffDef.empty()) {
        // RDKit✔️🔝:     mmffDef = defaultMMFFDef;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   std::istringstream inStream(mmffDef);
        // RDKit✔️🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   unsigned int oldAtomType = 0;
        // RDKit✔️🔝:   unsigned int atomType;
        // RDKit✔️🔝:   while (!(inStream.eof())) {
        // RDKit❗✔️:     if (inLine[0] != '*') {
        // RDKit✔️🔝:       MMFFDef mmffDefObj;
        // RDKit✔️🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit✔️🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit✔️🔝:       tokenizer::iterator token = tokens.begin();
        // RDKit✔️🔝:
        // RDKit✔️🔝:       // skip first token
        // RDKit✔️🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int atomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️🔝: #else
        // RDKit✔️🔝:       atomType = (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝: #endif
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       // Level 2 (currently = Level 1, see MMFF.I page 513)
        // RDKit✔️🔝:       mmffDefObj.eqLevel[0] =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       // Level 3
        // RDKit✔️🔝:       mmffDefObj.eqLevel[1] =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       // Level 4
        // RDKit✔️🔝:       mmffDefObj.eqLevel[2] =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       // Level 5
        // RDKit✔️🔝:       mmffDefObj.eqLevel[3] =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       if (atomType != oldAtomType) {
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:         d_params[atomType] = mmffDefObj;
        // RDKit✔️🔝: #else
        // RDKit✔️🔝:         d_params.push_back(mmffDefObj);
        // RDKit✔️🔝: #endif
        // RDKit✔️🔝:         oldAtomType = atomType;
        // RDKit✔️🔝:       }
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // Behavior review: for the frozen default vector build, choose the
        // embedded table only for empty input, ignore byte-zero-star comments,
        // parse the label/type/four levels in source order, narrow each value
        // to u8, and skip only a type equal to the last appended type. Missing
        // or invalid tokens return the packet's typed safety error because the
        // source indexes/dereferences missing tokens with undefined behavior.
        // Complexity review: borrow the selected source, scan each line and
        // token once, and append only retained fixed-size rows to one Vec. This
        // avoids the source's owned stream buffer, per-line strings and token
        // strings while preserving O(input bytes) parsing and vector order.
        let source = if text.is_empty() {
            include_str!("default_def.tsv")
        } else {
            text
        };

        let mut old_atom_type = 0_u32;
        let mut d_params = Vec::new();
        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }

            let mut tokens = tokenize_mmff_line(line);
            if tokens.next().is_none() {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Def,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            }

            let mut column = 1;
            let mut next_u8 = || -> Result<u8, MmffParamParseError> {
                let current_column = column;
                column += 1;
                let Some(cell) = tokens.next() else {
                    return Err(MmffParamParseError {
                        table: MmffParamTable::Def,
                        line: physical_line,
                        column: current_column,
                        cause: MmffParamParseCause::MissingToken,
                    });
                };
                let Some(value) = parse_mmff_u32(cell) else {
                    return Err(MmffParamParseError {
                        table: MmffParamTable::Def,
                        line: physical_line,
                        column: current_column,
                        cause: MmffParamParseCause::InvalidUnsigned {
                            cell: cell.to_owned(),
                        },
                    });
                };
                Ok(value as u8)
            };

            let atom_type = u32::from(next_u8()?);
            let level_2 = next_u8()?;
            let level_3 = next_u8()?;
            let level_4 = next_u8()?;
            let level_5 = next_u8()?;
            let eq_level = [level_2, level_3, level_4, level_5];
            if atom_type != old_atom_type {
                d_params.push(MmffDef { eq_level });
                old_atom_type = atom_type;
            }
        }

        Ok(Self { d_params })
    }

    pub(super) fn get(&self, atom_type: u32) -> Option<&MmffDef> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:171-179:
        // RDKit✔️✔️:   const MMFFDef *operator()(const unsigned int atomType) const {
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res = d_params.find(atomType);
        // RDKit❌❌:     return ((res != d_params.end()) ? &((*res).second) : NULL);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     return ((atomType && (atomType <= d_params.size()))
        // RDKit✔️✔️:                 ? &d_params[atomType - 1]
        // RDKit✔️✔️:                 : nullptr);
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:   }
        // The frozen vector build tests the one-based position, not the
        // stored atom-type key; a positive in-range value returns that slot.
        // usize conversion follows the source comparison against vector size.
        // The check plus direct indexing is constant-time and allocation-free.
        let atom_type = usize::try_from(atom_type).ok()?;
        if atom_type > 0 && atom_type <= self.d_params.len() {
            Some(&self.d_params[atom_type - 1])
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MmffAngle {
    pub(super) ka: f64,
    pub(super) theta0: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MmffAngleCollection {
    d_i_atom_type: Vec<u8>,
    d_j_atom_type: Vec<u8>,
    d_k_atom_type: Vec<u8>,
    d_angle_type: Vec<u8>,
    d_params: Vec<MmffAngle>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum MmffAngleLookupError {
    MissingDefinition { atom_type: u32 },
}

impl std::fmt::Display for MmffAngleLookupError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingDefinition { atom_type } => {
                write!(
                    formatter,
                    "missing MMFF definition for atom type {atom_type}"
                )
            }
        }
    }
}

impl std::error::Error for MmffAngleLookupError {}

impl MmffAngleCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:2095-2153:
        // RDKit❗🔝: MMFFAngleCollection::MMFFAngleCollection(std::string mmffAngle) {
        // RDKit❗🔝:   if (mmffAngle.empty()) {
        // RDKit❗🔝:     unsigned int i = 0;
        // RDKit❗🔝:     while (defaultMMFFAngleData[i] != "EOS") {
        // RDKit❗🔝:       mmffAngle += defaultMMFFAngleData[i];
        // RDKit❗🔝:       ++i;
        // RDKit❗🔝:     }
        // RDKit❗🔝:   }
        // RDKit❗🔝:   std::istringstream inStream(mmffAngle);
        // RDKit❗🔝:
        // RDKit❗🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   while (!(inStream.eof())) {
        // RDKit❗🔝:     if (inLine[0] != '*') {
        // RDKit❗🔝:       MMFFAngle mmffAngleObj;
        // RDKit❗🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit❗🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit❗🔝:       tokenizer::iterator token = tokens.begin();
        // RDKit❗🔝:
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int angleType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_angleType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_iAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_jAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int kAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_kAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffAngleObj.ka = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffAngleObj.theta0 = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[angleType][iAtomType][jAtomType][kAtomType] = mmffAngleObj;
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_params.push_back(mmffAngleObj);
        // RDKit❗🔝: #endif
        // RDKit❗🔝:     }
        // RDKit❗🔝:     inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   }
        // RDKit❗🔝: }
        // Behavior review: for the pinned vector build, empty input selects the
        // exact asset, consumed key cells parse as u32 then narrow modulo 256,
        // and ka/theta0 are appended in source row order without sorting or
        // duplicate removal. Unused suffix cells are ignored. Existing shared
        // line/token helpers preserve LF-only rows, one trailing CR removal,
        // byte-zero comments, and dropped empty TAB fields. A blank row or a
        // missing/invalid consumed token returns the existing typed safety
        // error instead of claiming parity for source undefined behavior; the
        // alternate STD_MAP branch remains unsupported.
        // Complexity review: the source concatenates owned lines and stages
        // owned stream/token text, then appends one row at a time. This parser
        // walks borrowed lines/tokens once and writes five source-parallel
        // vectors, preserving O(input bytes) parsing without staging copies,
        // sorting, or building an index.
        let source = if text.is_empty() {
            include_str!("default_angle.tsv")
        } else {
            text
        };

        let mut d_i_atom_type = Vec::new();
        let mut d_j_atom_type = Vec::new();
        let mut d_k_atom_type = Vec::new();
        let mut d_angle_type = Vec::new();
        let mut d_params = Vec::new();

        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(angle_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            };
            let Some(angle_type) = parse_mmff_u32(angle_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: angle_type_cell.to_owned(),
                    },
                });
            };

            let Some(i_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(i_atom_type) = parse_mmff_u32(i_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: i_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(j_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(j_atom_type) = parse_mmff_u32(j_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: j_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(k_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(k_atom_type) = parse_mmff_u32(k_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: k_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(ka_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(ka) = parse_mmff_f64(ka_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: ka_cell.to_owned(),
                    },
                });
            };

            let Some(theta0_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 5,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(theta0) = parse_mmff_f64(theta0_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: physical_line,
                    column: 5,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: theta0_cell.to_owned(),
                    },
                });
            };

            d_i_atom_type.push(i_atom_type as u8);
            d_j_atom_type.push(j_atom_type as u8);
            d_k_atom_type.push(k_atom_type as u8);
            d_angle_type.push(angle_type as u8);
            d_params.push(MmffAngle { ka, theta0 });
        }

        Ok(Self {
            d_i_atom_type,
            d_j_atom_type,
            d_k_atom_type,
            d_angle_type,
            d_params,
        })
    }

    pub(super) fn get(
        &self,
        definitions: &MmffDefCollection,
        angle_type: u32,
        i_atom_type: u32,
        j_atom_type: u32,
        k_atom_type: u32,
    ) -> Result<Option<&MmffAngle>, MmffAngleLookupError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:507-599:
        // RDKit✔️✔️:   const MMFFAngle *operator()(const MMFFDefCollection *mmffDef,
        // RDKit✔️✔️:                               const unsigned int angleType,
        // RDKit✔️✔️:                               const unsigned int iAtomType,
        // RDKit✔️✔️:                               const unsigned int jAtomType,
        // RDKit✔️✔️:                               const unsigned int kAtomType) const {
        // RDKit✔️✔️:     const MMFFAngle *mmffAngleParams = nullptr;
        // RDKit✔️✔️:     unsigned int iter = 0;
        //
        // RDKit✔️✔️: // For bending of the i-j-k angle, a five-stage process based
        // RDKit✔️✔️: // in the level combinations 1-1-1,2-2-2,3-2-3,4-2-4, and
        // RDKit✔️✔️: // 5-2-5 is used. (MMFF.I, note 68, page 519)
        // RDKit✔️✔️: // We skip 1-1-1 since Level 2 === Level 1
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
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
        // RDKit❌❌: #else
        // RDKit✔️✔️:     auto jBounds =
        // RDKit✔️✔️:         std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
        // RDKit✔️✔️:     if (jBounds.first != jBounds.second) {
        // RDKit✔️✔️:       while ((iter < 4) && (!mmffAngleParams)) {
        // RDKit❗✔️:         unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iter];
        // RDKit❗✔️:         unsigned int canKAtomType = (*mmffDef)(kAtomType)->eqLevel[iter];
        // RDKit✔️✔️:         if (canIAtomType > canKAtomType) {
        // RDKit✔️✔️:           std::swap(canIAtomType, canKAtomType);
        // RDKit✔️✔️:         }
        //
        // RDKit✔️✔️:         auto bounds = std::equal_range(
        // RDKit✔️✔️:             d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
        // RDKit✔️✔️:             d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
        // RDKit✔️✔️:             canIAtomType);
        // RDKit✔️✔️:         if (bounds.first != bounds.second) {
        // RDKit✔️✔️:           bounds = std::equal_range(
        // RDKit✔️✔️:               d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
        // RDKit✔️✔️:               d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
        // RDKit✔️✔️:               canKAtomType);
        // RDKit✔️✔️:           if (bounds.first != bounds.second) {
        // RDKit✔️✔️:             bounds = std::equal_range(
        // RDKit✔️✔️:                 d_angleType.begin() + (bounds.first - d_kAtomType.begin()),
        // RDKit✔️✔️:                 d_angleType.begin() + (bounds.second - d_kAtomType.begin()),
        // RDKit✔️✔️:                 angleType);
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
        // Behavior review: preserve the source's full-width central-j guard
        // before either positional Def lookup, then resolve i before k. Try
        // exactly four equivalence levels, canonicalize only the i/k pair,
        // and return the first angle row from the earliest matching stage.
        // The typed MissingDefinition error replaces only source null-pointer
        // dereference on missing i/k Def positions; valid misses stay None.
        // Complexity review: one central-j binary range and at most three
        // nested binary ranges for each of four stages, with constant-time
        // positional Def access. The closure allocates nothing and returns a
        // reference into this collection; no collection or Def rows are cloned.
        let equal_range = |values: &[u8], start: usize, end: usize, key: u32| {
            let values = &values[start..end];
            let lower = values.partition_point(|stored| u32::from(*stored) < key);
            let upper = values.partition_point(|stored| u32::from(*stored) <= key);
            (start + lower, start + upper)
        };

        let (j_start, j_end) = equal_range(
            &self.d_j_atom_type,
            0,
            self.d_j_atom_type.len(),
            j_atom_type,
        );
        if j_start == j_end {
            return Ok(None);
        }

        let i_definition =
            definitions
                .get(i_atom_type)
                .ok_or(MmffAngleLookupError::MissingDefinition {
                    atom_type: i_atom_type,
                })?;
        let k_definition =
            definitions
                .get(k_atom_type)
                .ok_or(MmffAngleLookupError::MissingDefinition {
                    atom_type: k_atom_type,
                })?;

        for stage in 0..4 {
            let mut canonical_i = i_definition.eq_level[stage];
            let mut canonical_k = k_definition.eq_level[stage];
            if canonical_i > canonical_k {
                std::mem::swap(&mut canonical_i, &mut canonical_k);
            }

            let (i_start, i_end) =
                equal_range(&self.d_i_atom_type, j_start, j_end, u32::from(canonical_i));
            if i_start == i_end {
                continue;
            }

            let (k_start, k_end) =
                equal_range(&self.d_k_atom_type, i_start, i_end, u32::from(canonical_k));
            if k_start == k_end {
                continue;
            }

            let (angle_start, angle_end) =
                equal_range(&self.d_angle_type, k_start, k_end, angle_type);
            if angle_start != angle_end {
                return Ok(Some(&self.d_params[angle_start]));
            }
        }

        Ok(None)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MmffOop {
    pub(super) koop: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MmffOopCollection {
    d_i_atom_type: Vec<u8>,
    d_j_atom_type: Vec<u8>,
    d_k_atom_type: Vec<u8>,
    d_l_atom_type: Vec<u8>,
    d_params: Vec<MmffOop>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum MmffOopLookupError {
    MissingDefinition { atom_type: u32 },
}

impl std::fmt::Display for MmffOopLookupError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingDefinition { atom_type } => {
                write!(
                    formatter,
                    "missing MMFF definition for atom type {atom_type}"
                )
            }
        }
    }
}

impl std::error::Error for MmffOopLookupError {}

impl MmffOopCollection {
    pub(super) fn from_text(is_mmff_s: bool, text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:4937-4988:
        // RDKit❗🔝: MMFFOopCollection::MMFFOopCollection(const bool isMMFFs, std::string mmffOop) {
        // RDKit❗🔝:   if (mmffOop.empty()) {
        // RDKit❗🔝:     mmffOop = (isMMFFs ? defaultMMFFsOop : defaultMMFFOop);
        // RDKit❗🔝:   }
        // RDKit❗🔝:   std::istringstream inStream(mmffOop);
        // RDKit❗🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   while (!(inStream.eof())) {
        // RDKit❗🔝:     if (inLine[0] != '*') {
        // RDKit❗🔝:       MMFFOop mmffOopObj;
        // RDKit❗🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit❗🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit❗🔝:       tokenizer::iterator token = tokens.begin();
        // RDKit❗🔝:
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_iAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_jAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int kAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_kAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int lAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_lAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffOopObj.koop = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[iAtomType][jAtomType][kAtomType][lAtomType] = mmffOopObj;
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_params.push_back(mmffOopObj);
        // RDKit❗🔝: #endif
        // RDKit❗🔝:     }
        // RDKit❗🔝:     inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   }
        // RDKit❗🔝: }
        // Behavior review: the source-selected vector build reads the matching
        // default only for empty input, then processes borrowed LF-terminated
        // lines and dropped-empty TAB tokens through the existing lexical
        // owners. It narrows each consumed u32 key with Rust's modulo-256 u8
        // cast, appends all five values in source order, keeps duplicates, and
        // ignores suffix tokens. Typed Oop parse errors replace source
        // undefined token dereferences or conversion exceptions for malformed
        // rows; this is a structural safety boundary, not native error parity.
        // The alternate RDKIT_MMFF_PARAMS_USE_STD_MAP branch is unmodeled.
        // Complexity review: this is one O(input bytes) borrowed pass with
        // amortized vector appends and no per-line or per-token owned copies.
        // It avoids the source stream/string staging while preserving row
        // order; the lookup's separate fixed-stack sort is not involved here.
        let source = if text.is_empty() {
            if is_mmff_s {
                include_str!("default_oop_s.tsv")
            } else {
                include_str!("default_oop.tsv")
            }
        } else {
            text
        };

        let mut d_i_atom_type = Vec::new();
        let mut d_j_atom_type = Vec::new();
        let mut d_k_atom_type = Vec::new();
        let mut d_l_atom_type = Vec::new();
        let mut d_params = Vec::new();

        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(i_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Oop,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            };
            let Some(i_atom_type) = parse_mmff_u32(i_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Oop,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: i_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(j_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Oop,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(j_atom_type) = parse_mmff_u32(j_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Oop,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: j_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(k_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Oop,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(k_atom_type) = parse_mmff_u32(k_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Oop,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: k_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(l_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Oop,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(l_atom_type) = parse_mmff_u32(l_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Oop,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: l_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(koop_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Oop,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(koop) = parse_mmff_f64(koop_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Oop,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: koop_cell.to_owned(),
                    },
                });
            };

            d_i_atom_type.push(i_atom_type as u8);
            d_j_atom_type.push(j_atom_type as u8);
            d_k_atom_type.push(k_atom_type as u8);
            d_l_atom_type.push(l_atom_type as u8);
            d_params.push(MmffOop { koop });
        }

        Ok(Self {
            d_i_atom_type,
            d_j_atom_type,
            d_k_atom_type,
            d_l_atom_type,
            d_params,
        })
    }

    pub(super) fn get(
        &self,
        definitions: &MmffDefCollection,
        i_atom_type: u32,
        j_atom_type: u32,
        k_atom_type: u32,
        l_atom_type: u32,
    ) -> Result<Option<&MmffOop>, MmffOopLookupError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:729-797:
        // RDKit❗🔝:   const MMFFOop *operator()(const MMFFDefCollection *mmffDef,
        // RDKit❗🔝:                             const unsigned int iAtomType,
        // RDKit❗🔝:                             const unsigned int jAtomType,
        // RDKit❗🔝:                             const unsigned int kAtomType,
        // RDKit❗🔝:                             const unsigned int lAtomType) const {
        // RDKit❗🔝:     const MMFFOop *mmffOopParams = nullptr;
        // RDKit❗🔝:     unsigned int iter = 0;
        // RDKit❗🔝:     std::vector<unsigned int> canIKLAtomType(3);
        // RDKit❗🔝: // For out-of-plane bending ijk; I , where j is the central
        // RDKit❗🔝: // atom [cf. eq. (511, the five-stage protocol 1-1-1; 1, 2-2-2; 2,
        // RDKit❗🔝: // 3-2-3;3, 4-2-4;4, 5-2-5;5 is used. The final stage provides
        // RDKit❗🔝: // wild-card defaults for all except the central atom.
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     while ((iter < 4) && (!mmffOopParams)) {
        // RDKit❌❌:       canIKLAtomType[0] = (*mmffDef)(iAtomType)->eqLevel[iter];
        // RDKit❌❌:       unsigned int canJAtomType = jAtomType;
        // RDKit❌❌:       canIKLAtomType[1] = (*mmffDef)(kAtomType)->eqLevel[iter];
        // RDKit❌❌:       canIKLAtomType[2] = (*mmffDef)(lAtomType)->eqLevel[iter];
        // RDKit❌❌:       std::sort(canIKLAtomType.begin(), canIKLAtomType.end());
        // RDKit❌❌:       const auto res1 = d_params.find(canIKLAtomType[0]);
        // RDKit❌❌:       if (res1 != d_params.end()) {
        // RDKit❌❌:         const auto res2 = ((*res1).second).find(canJAtomType);
        // RDKit❌❌:         if (res2 != ((*res1).second).end()) {
        // RDKit❌❌:           const auto res3 = ((*res2).second).find(canIKLAtomType[1]);
        // RDKit❌❌:           if (res3 != ((*res2).second).end()) {
        // RDKit❌❌:             const auto res4 = ((*res3).second).find(canIKLAtomType[2]);
        // RDKit❌❌:             if (res4 != ((*res3).second).end()) {
        // RDKit❌❌:               mmffOopParams = &((*res4).second);
        // RDKit❌❌:             }
        // RDKit❌❌:           }
        // RDKit❌❌:         }
        // RDKit❌❌:       }
        // RDKit❌❌:       ++iter;
        // RDKit❌❌:     }
        // RDKit❌❌: #else
        // RDKit❗🔝:     auto jBounds =
        // RDKit❗🔝:         std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
        // RDKit❗🔝:     if (jBounds.first != jBounds.second) {
        // RDKit❗🔝:       while ((iter < 4) && (!mmffOopParams)) {
        // RDKit❗🔝:         canIKLAtomType[0] = (*mmffDef)(iAtomType)->eqLevel[iter];
        // RDKit❗🔝:         canIKLAtomType[1] = (*mmffDef)(kAtomType)->eqLevel[iter];
        // RDKit❗🔝:         canIKLAtomType[2] = (*mmffDef)(lAtomType)->eqLevel[iter];
        // RDKit❗🔝:         std::sort(canIKLAtomType.begin(), canIKLAtomType.end());
        // RDKit❗🔝:         auto bounds = std::equal_range(
        // RDKit❗🔝:             d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
        // RDKit❗🔝:             d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
        // RDKit❗🔝:             canIKLAtomType[0]);
        // RDKit❗🔝:         if (bounds.first != bounds.second) {
        // RDKit❗🔝:           bounds = std::equal_range(
        // RDKit❗🔝:               d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
        // RDKit❗🔝:               d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
        // RDKit❗🔝:               canIKLAtomType[1]);
        // RDKit❗🔝:           if (bounds.first != bounds.second) {
        // RDKit❗🔝:             bounds = std::equal_range(
        // RDKit❗🔝:                 d_lAtomType.begin() + (bounds.first - d_kAtomType.begin()),
        // RDKit❗🔝:                 d_lAtomType.begin() + (bounds.second - d_kAtomType.begin()),
        // RDKit❗🔝:                 canIKLAtomType[2]);
        // RDKit❗🔝:             if (bounds.first != bounds.second) {
        // RDKit❗🔝:               mmffOopParams = &d_params[bounds.first - d_lAtomType.begin()];
        // RDKit❗🔝:             }
        // RDKit❗🔝:           }
        // RDKit❗🔝:         }
        // RDKit❗🔝:         ++iter;
        // RDKit❗🔝:       }
        // RDKit❗🔝:     }
        // RDKit❗🔝: #endif
        // RDKit❗🔝:
        // RDKit❗🔝:     return mmffOopParams;
        // RDKit❗🔝:   }
        // Behavior review: after the full-width central-j miss guard, read
        // positional Def rows in source i/k/l order. For each of four stages,
        // sort only their three equivalence values and search the j-bounded
        // rows by i, k, then l. Return the first source row at the first
        // matching stage. Missing Def positions replace only source null
        // dereferences with the packet's typed safety error; a valid miss
        // remains None. The alternate map build is not modeled.
        // Complexity review: the selected vector path does one central-j
        // binary range and at most three nested binary ranges in each of four
        // stages. The fixed stack array below avoids the source's temporary
        // three-element vector allocation without changing the sort order.
        let equal_range = |values: &[u8], start: usize, end: usize, key: u32| {
            let values = &values[start..end];
            let lower = values.partition_point(|stored| u32::from(*stored) < key);
            let upper = values.partition_point(|stored| u32::from(*stored) <= key);
            (start + lower, start + upper)
        };

        let (j_start, j_end) = equal_range(
            &self.d_j_atom_type,
            0,
            self.d_j_atom_type.len(),
            j_atom_type,
        );
        if j_start == j_end {
            return Ok(None);
        }

        let i_definition =
            definitions
                .get(i_atom_type)
                .ok_or(MmffOopLookupError::MissingDefinition {
                    atom_type: i_atom_type,
                })?;
        let k_definition =
            definitions
                .get(k_atom_type)
                .ok_or(MmffOopLookupError::MissingDefinition {
                    atom_type: k_atom_type,
                })?;
        let l_definition =
            definitions
                .get(l_atom_type)
                .ok_or(MmffOopLookupError::MissingDefinition {
                    atom_type: l_atom_type,
                })?;

        for stage in 0..4 {
            let mut canonical_ikl = [
                u32::from(i_definition.eq_level[stage]),
                u32::from(k_definition.eq_level[stage]),
                u32::from(l_definition.eq_level[stage]),
            ];
            canonical_ikl.sort_unstable();

            let (i_start, i_end) =
                equal_range(&self.d_i_atom_type, j_start, j_end, canonical_ikl[0]);
            if i_start == i_end {
                continue;
            }

            let (k_start, k_end) =
                equal_range(&self.d_k_atom_type, i_start, i_end, canonical_ikl[1]);
            if k_start == k_end {
                continue;
            }

            let (l_start, l_end) =
                equal_range(&self.d_l_atom_type, k_start, k_end, canonical_ikl[2]);
            if l_start != l_end {
                return Ok(Some(&self.d_params[l_start]));
            }
        }

        Ok(None)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MmffTor {
    pub(super) v1: f64,
    pub(super) v2: f64,
    pub(super) v3: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MmffTorCollection {
    d_tor_type: Vec<u8>,
    d_i_atom_type: Vec<u8>,
    d_j_atom_type: Vec<u8>,
    d_k_atom_type: Vec<u8>,
    d_l_atom_type: Vec<u8>,
    d_params: Vec<MmffTor>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum MmffTorLookupError {
    MissingDefinition { atom_type: u32 },
    InvalidEquivalentLevel { level: usize },
}

impl std::fmt::Display for MmffTorLookupError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingDefinition { atom_type } => {
                write!(
                    formatter,
                    "missing MMFF definition for atom type {atom_type}"
                )
            }
            Self::InvalidEquivalentLevel { level } => {
                write!(
                    formatter,
                    "MMFF torsion equivalence level {level} is outside 0..4"
                )
            }
        }
    }
}

impl std::error::Error for MmffTorLookupError {}

impl MmffTorCollection {
    pub(super) fn from_text(is_mmff_s: bool, text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:5249-5312:
        // RDKit❗🔝: MMFFTorCollection::MMFFTorCollection(const bool isMMFFs, std::string mmffTor) {
        // RDKit❗🔝:   if (mmffTor.empty()) {
        // RDKit❗🔝:     mmffTor = (isMMFFs ? defaultMMFFsTor : defaultMMFFTor);
        // RDKit❗🔝:   }
        // RDKit❗🔝:   std::istringstream inStream(mmffTor);
        // RDKit❗🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   while (!(inStream.eof())) {
        // RDKit❗🔝:     if (inLine[0] != '*') {
        // RDKit❗🔝:       MMFFTor mmffTorObj;
        // RDKit❗🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit❗🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit❗🔝:       tokenizer::iterator token = tokens.begin();
        // RDKit❗🔝:
        // RDKit❗🔝: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❗🔝:       unsigned int torType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_torType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❗🔝:       unsigned int iAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_iAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❗🔝:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_jAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❗🔝:       unsigned int kAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_kAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❗🔝:       unsigned int lAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_lAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffTorObj.V1 = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffTorObj.V2 = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffTorObj.V3 = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❗🔝:       d_params[torType][iAtomType][jAtomType][kAtomType][lAtomType] =
        // RDKit❗🔝:           mmffTorObj;
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_params.push_back(mmffTorObj);
        // RDKit❗🔝: #endif
        // RDKit❗🔝:     }
        // RDKit❗🔝:     inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   }
        // RDKit❗🔝: }
        // Behavior review: the pinned vector build selects the matching
        // default only for empty input; custom text ignores the variant flag.
        // It narrows each consumed u32 key modulo 256, parses V1/V2/V3 in
        // source order, keeps duplicate rows and ignores unused suffixes.
        // Shared line/token helpers retain LF row splitting, one trailing CR
        // removal, byte-zero comments, and dropped-empty TAB fields. Invalid
        // or missing consumed cells return typed safety errors in place of
        // source conversion exceptions/token dereferences; this does not claim
        // native malformed-input parity. The alternate STD_MAP branch remains
        // unmodeled.
        // Complexity review: one borrowed O(input bytes) pass appends each
        // consumed cell to five key vectors and one value vector. It does not
        // stage owned lines/tokens, sort, deduplicate, or construct an index.
        let source = if text.is_empty() {
            if is_mmff_s {
                include_str!("default_tor_s.tsv")
            } else {
                include_str!("default_tor.tsv")
            }
        } else {
            text
        };

        let mut d_tor_type = Vec::new();
        let mut d_i_atom_type = Vec::new();
        let mut d_j_atom_type = Vec::new();
        let mut d_k_atom_type = Vec::new();
        let mut d_l_atom_type = Vec::new();
        let mut d_params = Vec::new();

        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(tor_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            };
            let Some(tor_type) = parse_mmff_u32(tor_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: tor_type_cell.to_owned(),
                    },
                });
            };

            let Some(i_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(i_atom_type) = parse_mmff_u32(i_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: i_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(j_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(j_atom_type) = parse_mmff_u32(j_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: j_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(k_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(k_atom_type) = parse_mmff_u32(k_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: k_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(l_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(l_atom_type) = parse_mmff_u32(l_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: l_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(v1_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 5,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(v1) = parse_mmff_f64(v1_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 5,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: v1_cell.to_owned(),
                    },
                });
            };

            let Some(v2_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 6,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(v2) = parse_mmff_f64(v2_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 6,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: v2_cell.to_owned(),
                    },
                });
            };

            let Some(v3_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 7,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(v3) = parse_mmff_f64(v3_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Tor,
                    line: physical_line,
                    column: 7,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: v3_cell.to_owned(),
                    },
                });
            };

            d_tor_type.push(tor_type as u8);
            d_i_atom_type.push(i_atom_type as u8);
            d_j_atom_type.push(j_atom_type as u8);
            d_k_atom_type.push(k_atom_type as u8);
            d_l_atom_type.push(l_atom_type as u8);
            d_params.push(MmffTor { v1, v2, v3 });
        }

        Ok(Self {
            d_tor_type,
            d_i_atom_type,
            d_j_atom_type,
            d_k_atom_type,
            d_l_atom_type,
            d_params,
        })
    }

    pub(super) fn get<'a>(
        &'a self,
        definitions: &MmffDefCollection,
        tor_type: (u32, u32),
        i_atom_type: u32,
        j_atom_type: u32,
        k_atom_type: u32,
        l_atom_type: u32,
    ) -> Result<(u32, Option<&'a MmffTor>), MmffTorLookupError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:822-937:
        // RDKit❗✔️:   const std::pair<const unsigned int, const MMFFTor *> getMMFFTorParams(
        // RDKit❗✔️:       const MMFFDefCollection *mmffDef,
        // RDKit❗✔️:       const std::pair<unsigned int, unsigned int> torType,
        // RDKit❗✔️:       const unsigned int iAtomType, const unsigned int jAtomType,
        // RDKit❗✔️:       const unsigned int kAtomType, const unsigned int lAtomType) const {
        // RDKit❗✔️:     const MMFFTor *mmffTorParams = nullptr;
        // RDKit❗✔️:     unsigned int iter = 0;
        // RDKit❗✔️:     unsigned int iWildCard = 0;
        // RDKit❗✔️:     unsigned int lWildCard = 0;
        // RDKit❗✔️:     unsigned int canTorType = torType.first;
        // RDKit❗✔️:     unsigned int maxIter = 5;
        // RDKit❗✔️: // For i-j-k-2 torsion interactions, a five-stage
        // RDKit❗✔️: // process based on level combinations 1-1-1-1, 2-2-2-2,
        // RDKit❗✔️: // 3-2-2-5, 5-2-2-3, and 5-2-2-5 is used, where stages 3
        // RDKit❗✔️: // and 4 correspond to "half-default" or "half-wild-card" entries.
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌: #else
        // RDKit❌❌: #endif
        // RDKit❗✔️:
        // RDKit❗✔️:     while (((iter < maxIter) && ((!mmffTorParams) || (maxIter == 4))) ||
        // RDKit❗✔️:            ((iter == 4) && (torType.first == 5) && torType.second)) {
        // RDKit❗✔️:       // The rule of setting the torsion type to the value it had
        // RDKit❗✔️:       // before being set to 5 as a last resort in case parameters
        // RDKit❗✔️:       // could not be found is not mentioned in MMFF.IV; it was
        // RDKit❗✔️:       // empirically discovered due to a number of tests in the
        // RDKit❗✔️:       // MMFF validation suite otherwise failing
        // RDKit❗✔️:       if ((maxIter == 5) && (iter == 4)) {
        // RDKit❗✔️:         maxIter = 4;
        // RDKit❗✔️:         iter = 0;
        // RDKit❗✔️:         canTorType = torType.second;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       iWildCard = iter;
        // RDKit❗✔️:       lWildCard = iter;
        // RDKit❗✔️:       if (iter == 1) {
        // RDKit❗✔️:         iWildCard = 1;
        // RDKit❗✔️:         lWildCard = 3;
        // RDKit❗✔️:       } else if (iter == 2) {
        // RDKit❗✔️:         iWildCard = 3;
        // RDKit❗✔️:         lWildCard = 1;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iWildCard];
        // RDKit❗✔️:       unsigned int canJAtomType = jAtomType;
        // RDKit❗✔️:       unsigned int canKAtomType = kAtomType;
        // RDKit❗✔️:       unsigned int canLAtomType = (*mmffDef)(lAtomType)->eqLevel[lWildCard];
        // RDKit❗✔️:       if (canJAtomType > canKAtomType) {
        // RDKit❗✔️:         unsigned int temp = canKAtomType;
        // RDKit❗✔️:         canKAtomType = canJAtomType;
        // RDKit❗✔️:         canJAtomType = temp;
        // RDKit❗✔️:         temp = canLAtomType;
        // RDKit❗✔️:         canLAtomType = canIAtomType;
        // RDKit❗✔️:         canIAtomType = temp;
        // RDKit❗✔️:       } else if ((canJAtomType == canKAtomType) &&
        // RDKit❗✔️:                  (canIAtomType > canLAtomType)) {
        // RDKit❗✔️:         unsigned int temp = canLAtomType;
        // RDKit❗✔️:         canLAtomType = canIAtomType;
        // RDKit❗✔️:         canIAtomType = temp;
        // RDKit❗✔️:       }
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
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
        // RDKit❌❌: #else
        // RDKit❗✔️:       auto jBounds = std::equal_range(d_jAtomType.begin(), d_jAtomType.end(),
        // RDKit❗✔️:                                       canJAtomType);
        // RDKit❗✔️:       if (jBounds.first != jBounds.second) {
        // RDKit❗✔️:         auto bounds = std::equal_range(
        // RDKit❗✔️:             d_kAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
        // RDKit❗✔️:             d_kAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
        // RDKit❗✔️:             canKAtomType);
        // RDKit❗✔️:         if (bounds.first != bounds.second) {
        // RDKit❗✔️:           bounds = std::equal_range(
        // RDKit❗✔️:               d_iAtomType.begin() + (bounds.first - d_kAtomType.begin()),
        // RDKit❗✔️:               d_iAtomType.begin() + (bounds.second - d_kAtomType.begin()),
        // RDKit❗✔️:               canIAtomType);
        // RDKit❗✔️:           if (bounds.first != bounds.second) {
        // RDKit❗✔️:             bounds = std::equal_range(
        // RDKit❗✔️:                 d_lAtomType.begin() + (bounds.first - d_iAtomType.begin()),
        // RDKit❗✔️:                 d_lAtomType.begin() + (bounds.second - d_iAtomType.begin()),
        // RDKit❗✔️:                 canLAtomType);
        // RDKit❗✔️:             if (bounds.first != bounds.second) {
        // RDKit❗✔️:               bounds = std::equal_range(
        // RDKit❗✔️:                   d_torType.begin() + (bounds.first - d_lAtomType.begin()),
        // RDKit❗✔️:                   d_torType.begin() + (bounds.second - d_lAtomType.begin()),
        // RDKit❗✔️:                   canTorType);
        // RDKit❗✔️:               if (bounds.first != bounds.second) {
        // RDKit❗✔️:                 mmffTorParams = &d_params[bounds.first - d_torType.begin()];
        // RDKit❗✔️:                 if (maxIter == 4) {
        // RDKit❗✔️:                   break;
        // RDKit❗✔️:                 }
        // RDKit❗✔️:               }
        // RDKit❗✔️:             }
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❌❌: #endif
        // RDKit❗✔️:       ++iter;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     return std::make_pair(canTorType, mmffTorParams);
        // RDKit❗✔️:   }
        // Behavior review: preserve the source's five-iteration initial loop:
        // it performs up to four wildcard lookups, then transitions at iter
        // four by resetting the index and changing the lookup type to
        // tor_type.1 before trying the fallback stages.
        // Preserve stage order, asymmetric half-wildcards, definition lookup
        // order, J/K then I/L canonicalization, first equal-range row, and an
        // earlier initial hit unless the source's fifth-stage condition enters
        // fallback. Return the active lookup type on both hit and miss. Missing
        // definitions become typed errors in source dereference order; the
        // source's iter-four access outside the four-element Def array becomes
        // InvalidEquivalentLevel only after the I definition lookup.
        // Complexity review: at most four initial and four fallback table
        // searches, each with five nested binary range searches over existing
        // sorted vectors; a final source-undefined level access exits before
        // any table search. The local range closure and scalar state allocate
        // nothing, compare full-width queries against promoted u8 keys, and
        // return a borrowed row without cloning or scanning the table.
        let equal_range = |values: &[u8], start: usize, end: usize, key: u32| {
            let values = &values[start..end];
            let lower = values.partition_point(|stored| u32::from(*stored) < key);
            let upper = values.partition_point(|stored| u32::from(*stored) <= key);
            (start + lower, start + upper)
        };

        let mut mmff_tor_params: Option<&'a MmffTor> = None;
        let mut iter = 0_usize;
        let mut can_tor_type = tor_type.0;
        let mut max_iter = 5_usize;

        while ((iter < max_iter) && (mmff_tor_params.is_none() || max_iter == 4))
            || ((iter == 4) && (tor_type.0 == 5) && tor_type.1 != 0)
        {
            if max_iter == 5 && iter == 4 {
                max_iter = 4;
                iter = 0;
                can_tor_type = tor_type.1;
            }

            let (i_wildcard, l_wildcard) = match iter {
                1 => (1, 3),
                2 => (3, 1),
                other => (other, other),
            };

            let i_definition =
                definitions
                    .get(i_atom_type)
                    .ok_or(MmffTorLookupError::MissingDefinition {
                        atom_type: i_atom_type,
                    })?;
            let Some(&can_i_atom_type) = i_definition.eq_level.get(i_wildcard) else {
                return Err(MmffTorLookupError::InvalidEquivalentLevel { level: i_wildcard });
            };

            let mut can_j_atom_type = j_atom_type;
            let mut can_k_atom_type = k_atom_type;

            let l_definition =
                definitions
                    .get(l_atom_type)
                    .ok_or(MmffTorLookupError::MissingDefinition {
                        atom_type: l_atom_type,
                    })?;
            let Some(&can_l_atom_type) = l_definition.eq_level.get(l_wildcard) else {
                return Err(MmffTorLookupError::InvalidEquivalentLevel { level: l_wildcard });
            };

            let mut can_i_atom_type = u32::from(can_i_atom_type);
            let mut can_l_atom_type = u32::from(can_l_atom_type);
            if can_j_atom_type > can_k_atom_type {
                std::mem::swap(&mut can_j_atom_type, &mut can_k_atom_type);
                std::mem::swap(&mut can_i_atom_type, &mut can_l_atom_type);
            } else if can_j_atom_type == can_k_atom_type && can_i_atom_type > can_l_atom_type {
                std::mem::swap(&mut can_i_atom_type, &mut can_l_atom_type);
            }

            let (j_start, j_end) = equal_range(
                &self.d_j_atom_type,
                0,
                self.d_j_atom_type.len(),
                can_j_atom_type,
            );
            if j_start != j_end {
                let (k_start, k_end) =
                    equal_range(&self.d_k_atom_type, j_start, j_end, can_k_atom_type);
                if k_start != k_end {
                    let (i_start, i_end) =
                        equal_range(&self.d_i_atom_type, k_start, k_end, can_i_atom_type);
                    if i_start != i_end {
                        let (l_start, l_end) =
                            equal_range(&self.d_l_atom_type, i_start, i_end, can_l_atom_type);
                        if l_start != l_end {
                            let (tor_start, tor_end) =
                                equal_range(&self.d_tor_type, l_start, l_end, can_tor_type);
                            if tor_start != tor_end {
                                mmff_tor_params = Some(&self.d_params[tor_start]);
                                if max_iter == 4 {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            iter += 1;
        }

        Ok((can_tor_type, mmff_tor_params))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MmffPbci {
    pub(super) pbci: f64,
    pub(super) fcadj: f64,
}

#[derive(Debug)]
pub(super) struct MmffPbciCollection {
    d_params: Vec<MmffPbci>,
}

impl MmffPbciCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:586-620:
        // RDKit❗🔝: MMFFPBCICollection::MMFFPBCICollection(std::string mmffPBCI) {
        // RDKit✔️🔝:   if (mmffPBCI.empty()) {
        // RDKit✔️🔝:     mmffPBCI = defaultMMFFPBCI;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   std::istringstream inStream(mmffPBCI);
        // RDKit✔️🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   while (!(inStream.eof())) {
        // RDKit❗🔝:     if (inLine[0] != '*') {
        // RDKit✔️🔝:       MMFFPBCI mmffPBCIObj;
        // RDKit✔️🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit✔️🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit✔️🔝:       tokenizer::iterator token = tokens.begin();
        //
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       // IMPORTANT: skip the first field
        // RDKit❌❌:       ++token;
        // RDKit❌❌:       unsigned int atomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit✔️🔝:       // IMPORTANT: skip the first two fields
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝: #endif
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       mmffPBCIObj.pbci = boost::lexical_cast<double>(*token);
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       mmffPBCIObj.fcadj = boost::lexical_cast<double>(*token);
        // RDKit✔️🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[atomType] = mmffPBCIObj;
        // RDKit✔️🔝: #else
        // RDKit✔️🔝:       d_params.push_back(mmffPBCIObj);
        // RDKit✔️🔝: #endif
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // Behavior review: empty input selects the embedded default; the
        // vector build skips the first two tokens, parses PBCI then FCAdj,
        // ignores later fields, and appends every row in source order. A
        // zero-byte record and absent required tokens are represented by the
        // frozen typed CK safety errors because the source indexes/dereferences
        // those malformed records. Numeric conversion errors retain the
        // physical line, nonempty-token column, and only the failing cell.
        // Complexity review: one pass over borrowed lines and tab tokens,
        // appending each fixed-size row once. This keeps O(input bytes) work
        // while avoiding the source-owned stream/line/token staging copies.
        let source = if text.is_empty() {
            include_str!("default_pbci.tsv")
        } else {
            text
        };

        let mut d_params = Vec::new();
        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }
            if line.is_empty() {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Pbci,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            }

            let mut tokens = tokenize_mmff_line(line);
            for column in 0..2 {
                if tokens.next().is_none() {
                    return Err(MmffParamParseError {
                        table: MmffParamTable::Pbci,
                        line: physical_line,
                        column,
                        cause: MmffParamParseCause::MissingToken,
                    });
                }
            }

            let Some(pbci_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Pbci,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(pbci) = parse_mmff_f64(pbci_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Pbci,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: pbci_cell.to_owned(),
                    },
                });
            };

            let Some(fcadj_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Pbci,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(fcadj) = parse_mmff_f64(fcadj_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Pbci,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: fcadj_cell.to_owned(),
                    },
                });
            };

            d_params.push(MmffPbci { pbci, fcadj });
        }

        Ok(Self { d_params })
    }

    pub(super) fn get(&self, atom_type: u32) -> Option<&MmffPbci> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:225-235:
        // RDKit✔️✔️:   const MMFFPBCI *operator()(const unsigned int atomType) const {
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res = d_params.find(atomType);
        // RDKit❌❌:     return ((res != d_params.end()) ? &((*res).second) : NULL);
        // RDKit✔️🔝: #else
        // RDKit✔️🔝:     return ((atomType && (atomType <= d_params.size()))
        // RDKit✔️🔝:                 ? &d_params[atomType - 1]
        // RDKit✔️🔝:                 : nullptr);
        // RDKit❌❌: #endif
        // RDKit✔️✔️:   }
        //
        // The vector contract is one-based positional lookup; the row's
        // skipped label is not consulted. The guarded direct index is O(1),
        // borrows the stored row, and allocates nothing.
        let atom_type = usize::try_from(atom_type).ok()?;
        if atom_type > 0 && atom_type <= self.d_params.len() {
            Some(&self.d_params[atom_type - 1])
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MmffChg {
    pub(super) bci: f64,
}

#[derive(Debug)]
pub(super) struct MmffChgCollection {
    d_params: Vec<MmffChg>,
    d_i_atom_type: Vec<u8>,
    d_j_atom_type: Vec<u8>,
    d_bond_type: Vec<u8>,
}

impl MmffChgCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:732-776:
        // RDKit❗🔝: MMFFChgCollection::MMFFChgCollection(std::string mmffChg) {
        // RDKit❗🔝:   if (mmffChg.empty()) {
        // RDKit❗🔝:     mmffChg = defaultMMFFChg;
        // RDKit❗🔝:   }
        // RDKit❗🔝:   std::istringstream inStream(mmffChg);
        // RDKit❗🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   while (!(inStream.eof())) {
        // RDKit❗🔝:     if (inLine[0] != '*') {
        // RDKit❗🔝:       MMFFChg mmffChgObj;
        // RDKit❗🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit❗🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit❗🔝:       tokenizer::iterator token = tokens.begin();
        //
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int bondType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_bondType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_iAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_jAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffChgObj.bci = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[bondType][iAtomType][jAtomType] = mmffChgObj;
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_params.push_back(mmffChgObj);
        // RDKit❌❌: #endif
        // RDKit❗🔝:     }
        // RDKit❗🔝:     inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   }
        // RDKit❗🔝: }
        // Behavior review: the empty source selects the committed default
        // asset; source_lines preserves RDKit's final unterminated-record
        // behavior and removes one CR. Only byte-zero-star lines are comments.
        // Numeric conversion order is bond, i atom, j atom, then BCI. The
        // source's uint8_t casts keep the low eight bits, including wrapped
        // negative spellings accepted by parse_mmff_u32. Successful rows are
        // appended to four parallel vectors in source order; duplicates remain
        // and no sort or canonicalization is performed. Missing cells and
        // conversion failures retain their first consumed token as a typed
        // structural error instead of dereferencing an absent source token.
        // Complexity review: one borrowed line/token pass and four amortized
        // vector appends per row, O(input bytes) overall. This avoids source
        // stream/line/token staging copies; no secondary index is built.
        let source = if text.is_empty() {
            include_str!("default_chg.tsv")
        } else {
            text
        };

        let mut d_params = Vec::new();
        let mut d_i_atom_type = Vec::new();
        let mut d_j_atom_type = Vec::new();
        let mut d_bond_type = Vec::new();

        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }
            if line.is_empty() {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Chg,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(bond_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Chg,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(bond_type) = parse_mmff_u32(bond_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Chg,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: bond_type_cell.to_owned(),
                    },
                });
            };

            let Some(i_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Chg,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(i_atom_type) = parse_mmff_u32(i_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Chg,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: i_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(j_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Chg,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(j_atom_type) = parse_mmff_u32(j_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Chg,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: j_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(bci_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Chg,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(bci) = parse_mmff_f64(bci_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Chg,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: bci_cell.to_owned(),
                    },
                });
            };

            d_bond_type.push(bond_type as u8);
            d_i_atom_type.push(i_atom_type as u8);
            d_j_atom_type.push(j_atom_type as u8);
            d_params.push(MmffChg { bci });
        }

        Ok(Self {
            d_params,
            d_i_atom_type,
            d_j_atom_type,
            d_bond_type,
        })
    }

    pub(super) fn get(
        &self,
        bond_type: u32,
        i_atom_type: u32,
        j_atom_type: u32,
    ) -> (i32, Option<&MmffChg>) {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:251-292:
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
        // RDKit❗✔️:     }
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res1 = d_params[bondType].find(canIAtomType);
        // RDKit❌❌:     if (res1 != d_params[bondType].end()) {
        // RDKit❌❌:       const auto res2 = ((*res1).second).find(canJAtomType);
        // RDKit❌❌:       if (res2 != ((*res1).second).end()) {
        // RDKit❌❌:         mmffChgParams = &((*res2).second);
        // RDKit❌❌:       }
        // RDKit❌❌:     }
        // RDKit❌❌: #else
        // RDKit✔️✔️:     auto bounds =
        // RDKit✔️✔️:         std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), canIAtomType);
        // RDKit✔️✔️:     if (bounds.first != bounds.second) {
        // RDKit✔️✔️:       bounds = std::equal_range(
        // RDKit✔️✔️:           d_jAtomType.begin() + (bounds.first - d_iAtomType.begin()),
        // RDKit✔️✔️:           d_jAtomType.begin() + (bounds.second - d_iAtomType.begin()),
        // RDKit✔️✔️:           canJAtomType);
        // RDKit✔️✔️:       if (bounds.first != bounds.second) {
        // RDKit✔️✔️:         bounds = std::equal_range(
        // RDKit✔️✔️:             d_bondType.begin() + (bounds.first - d_jAtomType.begin()),
        // RDKit✔️✔️:             d_bondType.begin() + (bounds.second - d_jAtomType.begin()),
        // RDKit✔️✔️:             bondType);
        // RDKit✔️✔️:         if (bounds.first != bounds.second) {
        // RDKit✔️✔️:           mmffChgParams = &d_params[bounds.first - d_bondType.begin()];
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit❌❌: #endif
        //
        // RDKit✔️✔️:     return std::make_pair(sign, mmffChgParams);
        // RDKit✔️✔️:   }
        // Behavior review: swap only when the original i query is greater
        // than j; equal queries keep sign -1. Each stored u8 is widened for
        // comparison with the full-u32 query. The nested ranges are the vector
        // equal_range sequence: i over all rows, j within i, then bond within
        // the matching j range. The final lower bound returns the first
        // duplicate; misses preserve the previously computed sign. Input must
        // retain the source-required partition order; this code does not sort.
        // Complexity review: each source equal_range is represented by lower
        // and upper partition points over random-access slices, giving three
        // nested O(log N) bound searches and no allocation, clone or scan.
        let (canonical_i_atom_type, canonical_j_atom_type, sign) = if i_atom_type > j_atom_type {
            (j_atom_type, i_atom_type, 1)
        } else {
            (i_atom_type, j_atom_type, -1)
        };

        let i_start = self
            .d_i_atom_type
            .partition_point(|&stored| u32::from(stored) < canonical_i_atom_type);
        let i_end = self
            .d_i_atom_type
            .partition_point(|&stored| u32::from(stored) <= canonical_i_atom_type);
        if i_start == i_end {
            return (sign, None);
        }

        let j_range = &self.d_j_atom_type[i_start..i_end];
        let j_start_in_range =
            j_range.partition_point(|&stored| u32::from(stored) < canonical_j_atom_type);
        let j_end_in_range =
            j_range.partition_point(|&stored| u32::from(stored) <= canonical_j_atom_type);
        if j_start_in_range == j_end_in_range {
            return (sign, None);
        }
        let j_start = i_start + j_start_in_range;
        let j_end = i_start + j_end_in_range;

        let bond_range = &self.d_bond_type[j_start..j_end];
        let bond_start_in_range =
            bond_range.partition_point(|&stored| u32::from(stored) < bond_type);
        let bond_end_in_range =
            bond_range.partition_point(|&stored| u32::from(stored) <= bond_type);
        if bond_start_in_range == bond_end_in_range {
            return (sign, None);
        }

        let row_index = j_start + bond_start_in_range;
        (sign, Some(&self.d_params[row_index]))
    }
}

const DEFAULT_HERSCHBACH_LAURIE_TEXT: &str = include_str!("default_herschbach_laurie.tsv");

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MmffHerschbachLaurie {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Params.h:88-93:
    // RDKit✔️✔️: class RDKIT_FORCEFIELD_EXPORT MMFFHerschbachLaurie {
    // RDKit✔️✔️:  public:
    // RDKit✔️✔️:   double a_ij;
    // RDKit✔️✔️:   double d_ij;
    // RDKit✔️✔️:   double dp_ij;
    // RDKit✔️✔️: };
    pub(super) a_ij: f64,
    pub(super) d_ij: f64,
    pub(super) dp_ij: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MmffHerschbachLaurieCollection {
    d_i_row: Vec<u8>,
    d_j_row: Vec<u8>,
    d_params: Vec<MmffHerschbachLaurie>,
}

impl MmffHerschbachLaurieCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:1962-2003:
        // RDKit✔️🔝: MMFFHerschbachLaurieCollection::MMFFHerschbachLaurieCollection(
        // RDKit✔️🔝:     std::string mmffHerschbachLaurie) {
        // RDKit✔️🔝:   if (mmffHerschbachLaurie.empty()) {
        // RDKit✔️🔝:     mmffHerschbachLaurie = defaultMMFFHerschbachLaurie;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   std::istringstream inStream(mmffHerschbachLaurie);
        // RDKit✔️🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   while (!(inStream.eof())) {
        // RDKit❗✔️:     if (inLine[0] != '*') {
        // RDKit✔️🔝:       MMFFHerschbachLaurie mmffHerschbachLaurieObj;
        // RDKit✔️🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit✔️🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit✔️🔝:       tokenizer::iterator token = tokens.begin();
        //
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iRow = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️🔝: #else
        // RDKit❗🔝:       d_iRow.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jRow = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️🔝: #else
        // RDKit❗🔝:       d_jRow.push_back((std::uint8_t)boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffHerschbachLaurieObj.a_ij = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffHerschbachLaurieObj.d_ij = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffHerschbachLaurieObj.dp_ij = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[iRow][jRow] = mmffHerschbachLaurieObj;
        // RDKit✔️🔝: #else
        // RDKit✔️🔝:       d_params.push_back(mmffHerschbachLaurieObj);
        // RDKit✔️🔝: #endif
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // Behavior review: empty input selects the exact embedded source table;
        // only byte-zero-star lines are comments. The five consumed fields are
        // parsed in source order; key casts retain the low eight bits, later
        // tokens are ignored, and rows including duplicates remain in input
        // order. Missing or malformed cells and an empty processed line return
        // the existing typed safety errors instead of source undefined behavior.
        // Complexity review: one pass over borrowed newline records and tab
        // tokens, with three amortized vector appends per row. This preserves
        // row order and the two aligned keys while avoiding source stream,
        // line-string, and token-string copies: O(input bytes), O(rows) storage.
        let source = if text.is_empty() {
            DEFAULT_HERSCHBACH_LAURIE_TEXT
        } else {
            text
        };

        let mut d_i_row = Vec::new();
        let mut d_j_row = Vec::new();
        let mut d_params = Vec::new();

        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }
            if line.is_empty() {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(i_row_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(i_row) = parse_mmff_u32(i_row_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: i_row_cell.to_owned(),
                    },
                });
            };

            let Some(j_row_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(j_row) = parse_mmff_u32(j_row_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: j_row_cell.to_owned(),
                    },
                });
            };

            let Some(a_ij_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(a_ij) = parse_mmff_f64(a_ij_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: a_ij_cell.to_owned(),
                    },
                });
            };

            let Some(d_ij_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(d_ij) = parse_mmff_f64(d_ij_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: d_ij_cell.to_owned(),
                    },
                });
            };

            let Some(dp_ij_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(dp_ij) = parse_mmff_f64(dp_ij_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::HerschbachLaurie,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: dp_ij_cell.to_owned(),
                    },
                });
            };

            d_i_row.push(i_row as u8);
            d_j_row.push(j_row as u8);
            d_params.push(MmffHerschbachLaurie { a_ij, d_ij, dp_ij });
        }

        Ok(Self {
            d_i_row,
            d_j_row,
            d_params,
        })
    }

    pub(super) fn get(&self, i_row: i32, j_row: i32) -> Option<&MmffHerschbachLaurie> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:433-463:
        // RDKit✔️✔️:   const MMFFHerschbachLaurie *operator()(const int iRow, const int jRow) const {
        // RDKit✔️✔️:     const MMFFHerschbachLaurie *mmffHerschbachLaurieParams = nullptr;
        // RDKit✔️✔️:     unsigned int canIRow = iRow;
        // RDKit✔️✔️:     unsigned int canJRow = jRow;
        // RDKit✔️✔️:     if (iRow > jRow) {
        // RDKit✔️✔️:       canIRow = jRow;
        // RDKit✔️✔️:       canJRow = iRow;
        // RDKit✔️✔️:     }
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res1 = d_params.find(canIRow);
        // RDKit❌❌:     std::map<const unsigned int, MMFFHerschbachLaurie>::const_iterator res2;
        // RDKit❌❌:     if (res1 != d_params.end()) {
        // RDKit❌❌:       res2 = ((*res1).second).find(canJRow);
        // RDKit❌❌:       if (res2 != ((*res1).second).end()) {
        // RDKit❌❌:         mmffHerschbachLaurieParams = &((*res2).second);
        // RDKit❌❌:       }
        // RDKit❌❌:     }
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto bounds = std::equal_range(d_iRow.begin(), d_iRow.end(), canIRow);
        // RDKit✔️✔️:     if (bounds.first != bounds.second) {
        // RDKit✔️✔️:       bounds = std::equal_range(
        // RDKit✔️✔️:           d_jRow.begin() + (bounds.first - d_iRow.begin()),
        // RDKit✔️✔️:           d_jRow.begin() + (bounds.second - d_iRow.begin()), canJRow);
        // RDKit✔️✔️:       if (bounds.first != bounds.second) {
        // RDKit✔️✔️:         mmffHerschbachLaurieParams = &d_params[bounds.first - d_jRow.begin()];
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit❌❌: #endif
        //
        // RDKit✔️✔️:     return mmffHerschbachLaurieParams;
        // RDKit✔️✔️:   }
        // Behavior review: copy each original signed argument to its u32 source
        // key before comparing the original signed pair. Only `i_row > j_row`
        // swaps the converted keys. Two nested lower/upper-bound searches over
        // the source-partitioned key vectors return the first duplicate row;
        // full-width queries are never narrowed and stored rows are borrowed.
        // Complexity review: two nested O(log N) range searches, no scan,
        // allocation, sorting, or row clone. This matches equal_range's
        // comparison complexity for source-sorted rows; no material lookup
        // improvement over the source search is established.
        let mut can_i_row = i_row as u32;
        let mut can_j_row = j_row as u32;
        if i_row > j_row {
            std::mem::swap(&mut can_i_row, &mut can_j_row);
        }

        let i_start = self
            .d_i_row
            .partition_point(|&stored| u32::from(stored) < can_i_row);
        let i_end = self
            .d_i_row
            .partition_point(|&stored| u32::from(stored) <= can_i_row);
        if i_start == i_end {
            return None;
        }

        let j_range = &self.d_j_row[i_start..i_end];
        let j_start_in_range = j_range.partition_point(|&stored| u32::from(stored) < can_j_row);
        let j_end_in_range = j_range.partition_point(|&stored| u32::from(stored) <= can_j_row);
        if j_start_in_range == j_end_in_range {
            return None;
        }

        Some(&self.d_params[i_start + j_start_in_range])
    }
}

const DEFAULT_COV_RAD_PAU_ELE_TEXT: &str = include_str!("default_cov_rad_pau_ele.tsv");

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MmffCovRadPauEle {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Params.h:94-102:
    // RDKit✔️✔️: //! class to store covalent radius and Pauling electronegativity
    // RDKit✔️✔️: //! values for MMFF bond stretching empirical rule
    // RDKit✔️✔️: class RDKIT_FORCEFIELD_EXPORT MMFFCovRadPauEle {
    // RDKit✔️✔️:  public:
    // RDKit✔️✔️:   double r0;
    // RDKit✔️✔️:   double chi;
    // RDKit✔️✔️: };
    pub(super) r0: f64,
    pub(super) chi: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MmffCovRadPauEleCollection {
    d_atomic_num: Vec<u8>,
    d_params: Vec<MmffCovRadPauEle>,
}

impl MmffCovRadPauEleCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:2037-2070:
        // RDKit✔️🔝: MMFFCovRadPauEleCollection::MMFFCovRadPauEleCollection(
        // RDKit✔️🔝:     std::string mmffCovRadPauEle) {
        // RDKit✔️🔝:   if (mmffCovRadPauEle.empty()) {
        // RDKit✔️🔝:     mmffCovRadPauEle = defaultMMFFCovRadPauEle;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   std::istringstream inStream(mmffCovRadPauEle);
        // RDKit✔️🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   while (!(inStream.eof())) {
        // RDKit❗✔️:     if (inLine[0] != '*') {
        // RDKit✔️🔝:       MMFFCovRadPauEle mmffCovRadPauEleObj;
        // RDKit✔️🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit✔️🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit✔️🔝:       tokenizer::iterator token = tokens.begin();
        // RDKit❌❌:
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int atomicNum = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️🔝: #else
        // RDKit❗🔝:       d_atomicNum.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffCovRadPauEleObj.r0 = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffCovRadPauEleObj.chi = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[atomicNum] = mmffCovRadPauEleObj;
        // RDKit✔️🔝: #else
        // RDKit✔️🔝:       d_params.push_back(mmffCovRadPauEleObj);
        // RDKit✔️🔝: #endif
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // Behavior review: empty input selects the pinned CovRadPauEle text;
        // byte-zero-star records are comments. Each processed record consumes
        // the atomic number, r0, and chi in source order, narrows only the
        // stored key to u8, and ignores later tokens. Rows and duplicates stay
        // in input order. Missing/malformed consumed tokens and blank rows use
        // the existing typed input-safety errors in place of source UB.
        // Complexity review: scan borrowed records and tokens once, then append
        // one key and one fixed-size value per row. This keeps O(input bytes)
        // parsing and O(rows) storage while avoiding source stream, line, and
        // token string copies.
        let source = if text.is_empty() {
            DEFAULT_COV_RAD_PAU_ELE_TEXT
        } else {
            text
        };

        let mut d_atomic_num = Vec::new();
        let mut d_params = Vec::new();
        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }
            if line.is_empty() {
                return Err(MmffParamParseError {
                    table: MmffParamTable::CovRadPauEle,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(atomic_num_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::CovRadPauEle,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(atomic_num) = parse_mmff_u32(atomic_num_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::CovRadPauEle,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: atomic_num_cell.to_owned(),
                    },
                });
            };

            let Some(r0_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::CovRadPauEle,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(r0) = parse_mmff_f64(r0_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::CovRadPauEle,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: r0_cell.to_owned(),
                    },
                });
            };

            let Some(chi_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::CovRadPauEle,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(chi) = parse_mmff_f64(chi_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::CovRadPauEle,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: chi_cell.to_owned(),
                    },
                });
            };

            d_atomic_num.push(atomic_num as u8);
            d_params.push(MmffCovRadPauEle { r0, chi });
        }

        Ok(Self {
            d_atomic_num,
            d_params,
        })
    }

    pub(super) fn get(&self, atomic_num: u32) -> Option<&MmffCovRadPauEle> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:484-495:
        // RDKit✔️✔️:   const MMFFCovRadPauEle *operator()(const unsigned int atomicNum) const {
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res = d_params.find(atomicNum);
        // RDKit❌❌:     return ((res != d_params.end()) ? &((*res).second) : NULL);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto bounds =
        // RDKit✔️✔️:         std::equal_range(d_atomicNum.begin(), d_atomicNum.end(), atomicNum);
        // RDKit✔️✔️:     return ((bounds.first != bounds.second)
        // RDKit✔️✔️:                 ? &d_params[bounds.first - d_atomicNum.begin()]
        // RDKit✔️✔️:                 : nullptr);
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:   }
        // Behavior review: compare the full u32 query against widened stored
        // u8 keys, returning the first equal-range row by borrow; queries are
        // never narrowed and duplicate source rows remain observable as the
        // first matching row. The vector source assumes sorted keys; this
        // lookup does not sort or otherwise rewrite input rows.
        // Complexity review: two logarithmic partition points implement the
        // source equal_range without allocation, cloning, or a linear scan.
        let start = self
            .d_atomic_num
            .partition_point(|&stored| u32::from(stored) < atomic_num);
        let end = self
            .d_atomic_num
            .partition_point(|&stored| u32::from(stored) <= atomic_num);
        (start != end).then(|| &self.d_params[start])
    }
}

const DEFAULT_MMFF_STBN_TEXT: &str = include_str!("default_stbn.tsv");

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MmffStbn {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Params.h:110-115:
    // RDKit✔️✔️: //! class to store MMFF parameters for stretch-bending
    // RDKit✔️✔️: class RDKIT_FORCEFIELD_EXPORT MMFFStbn {
    // RDKit✔️✔️:  public:
    // RDKit✔️✔️:   double kbaIJK;
    // RDKit✔️✔️:   double kbaKJI;
    // RDKit✔️✔️: };
    pub(super) kba_ijk: f64,
    pub(super) kba_kji: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MmffStbnCollection {
    d_i_atom_type: Vec<u8>,
    d_j_atom_type: Vec<u8>,
    d_k_atom_type: Vec<u8>,
    d_stretch_bend_type: Vec<u8>,
    d_params: Vec<MmffStbn>,
}

impl MmffStbnCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:4516-4570:
        // RDKit✔️✔️: MMFFStbnCollection::MMFFStbnCollection(std::string mmffStbn) {
        // RDKit✔️✔️:   if (mmffStbn.empty()) {
        // RDKit✔️✔️:     mmffStbn = defaultMMFFStbn;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   std::istringstream inStream(mmffStbn);
        // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
        // RDKit✔️✔️:   while (!(inStream.eof())) {
        // RDKit❗✔️:     if (inLine[0] != '*') {
        // RDKit✔️✔️:       MMFFStbn mmffStbnObj;
        // RDKit✔️✔️:       boost::char_separator<char> tabSep("\t");
        // RDKit✔️✔️:       tokenizer tokens(inLine, tabSep);
        // RDKit✔️✔️:       tokenizer::iterator token = tokens.begin();
        // RDKit✔️✔️:
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int stretchBendType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_stretchBendType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_iAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_jAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int kAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_kAtomType.push_back(
        // RDKit✔️✔️:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       mmffStbnObj.kbaIJK = boost::lexical_cast<double>(*token);
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       mmffStbnObj.kbaKJI = boost::lexical_cast<double>(*token);
        // RDKit✔️✔️:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[stretchBendType][iAtomType][jAtomType][kAtomType] = mmffStbnObj;
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:       d_params.push_back(mmffStbnObj);
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     inLine = RDKit::getLine(inStream);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // Behavior review: empty text selects the committed source default;
        // each LF-terminated, non-comment row consumes stretch-bend type, i/j/k,
        // then kbaIJK/kbaKJI in source order. The vector mode narrows only the
        // four stored keys to u8, retains row order and duplicates, and ignores
        // unused suffix tokens. The source map alternative remains unmodeled.
        // Existing borrowed line/token and numeric helpers preserve the shared
        // MMFF table grammar. For source-undefined empty/missing/bad tokens,
        // CK returns the frozen typed safety errors instead of invoking UB or
        // exposing a partial collection.
        // Complexity review: one borrowed pass over records and tokens with
        // four key appends and one fixed-value append per row; O(input bytes)
        // parsing, O(rows) storage, and no per-line/token owned string copies.
        let source = if text.is_empty() {
            DEFAULT_MMFF_STBN_TEXT
        } else {
            text
        };

        let mut d_i_atom_type = Vec::new();
        let mut d_j_atom_type = Vec::new();
        let mut d_k_atom_type = Vec::new();
        let mut d_stretch_bend_type = Vec::new();
        let mut d_params = Vec::new();

        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }
            if line.is_empty() {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(stretch_bend_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(stretch_bend_type) = parse_mmff_u32(stretch_bend_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: stretch_bend_type_cell.to_owned(),
                    },
                });
            };

            let Some(i_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(i_atom_type) = parse_mmff_u32(i_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: i_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(j_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(j_atom_type) = parse_mmff_u32(j_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: j_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(k_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(k_atom_type) = parse_mmff_u32(k_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: k_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(kba_ijk_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(kba_ijk) = parse_mmff_f64(kba_ijk_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: kba_ijk_cell.to_owned(),
                    },
                });
            };

            let Some(kba_kji_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 5,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(kba_kji) = parse_mmff_f64(kba_kji_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Stbn,
                    line: physical_line,
                    column: 5,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: kba_kji_cell.to_owned(),
                    },
                });
            };

            d_stretch_bend_type.push(stretch_bend_type as u8);
            d_i_atom_type.push(i_atom_type as u8);
            d_j_atom_type.push(j_atom_type as u8);
            d_k_atom_type.push(k_atom_type as u8);
            d_params.push(MmffStbn { kba_ijk, kba_kji });
        }

        Ok(Self {
            d_i_atom_type,
            d_j_atom_type,
            d_k_atom_type,
            d_stretch_bend_type,
            d_params,
        })
    }

    pub(super) fn get(
        &self,
        stretch_bend_type: u32,
        bond_type1: u32,
        bond_type2: u32,
        i_atom_type: u32,
        j_atom_type: u32,
        k_atom_type: u32,
    ) -> (bool, Option<&MmffStbn>) {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:607-665:
        // RDKit✔️✔️:   const std::pair<bool, const MMFFStbn *> getMMFFStbnParams(
        // RDKit✔️✔️:       const unsigned int stretchBendType, const unsigned int bondType1,
        // RDKit✔️✔️:       const unsigned int bondType2, const unsigned int iAtomType,
        // RDKit✔️✔️:       const unsigned int jAtomType, const unsigned int kAtomType) const {
        // RDKit✔️✔️:     const MMFFStbn *mmffStbnParams = nullptr;
        // RDKit✔️✔️:     bool swap = false;
        // RDKit✔️✔️:     unsigned int canIAtomType = iAtomType;
        // RDKit✔️✔️:     unsigned int canKAtomType = kAtomType;
        // RDKit✔️✔️:     unsigned int canStretchBendType = stretchBendType;
        // RDKit✔️✔️:     if (iAtomType > kAtomType) {
        // RDKit✔️✔️:       canIAtomType = kAtomType;
        // RDKit✔️✔️:       canKAtomType = iAtomType;
        // RDKit✔️✔️:       swap = true;
        // RDKit✔️✔️:     } else if (iAtomType == kAtomType) {
        // RDKit✔️✔️:       swap = (bondType1 < bondType2);
        // RDKit✔️✔️:     }
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
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
        // RDKit✔️✔️:     if (jBounds.first != jBounds.second) {
        // RDKit✔️✔️:       auto bounds = std::equal_range(
        // RDKit✔️✔️:           d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
        // RDKit✔️✔️:           d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
        // RDKit✔️✔️:           canIAtomType);
        // RDKit✔️✔️:       if (bounds.first != bounds.second) {
        // RDKit✔️✔️:         bounds = std::equal_range(
        // RDKit✔️✔️:             d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
        // RDKit✔️✔️:             d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
        // RDKit✔️✔️:             canKAtomType);
        // RDKit✔️✔️:         if (bounds.first != bounds.second) {
        // RDKit✔️✔️:           bounds = std::equal_range(
        // RDKit✔️✔️:               d_stretchBendType.begin() + (bounds.first - d_kAtomType.begin()),
        // RDKit✔️✔️:               d_stretchBendType.begin() + (bounds.second - d_kAtomType.begin()),
        // RDKit✔️✔️:               canStretchBendType);
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
        // Behavior review: compare i/k using the original full-width u32 inputs;
        // reverse terminal atom keys only when i>k, and use bond-type order for
        // the equal-terminal swap flag. The vector-mode four-level equal_range
        // order is j, i, k, stretch-bend type; queries are never narrowed and
        // the first matching stored value is borrowed. The computed swap is
        // returned on misses too. The alternate source map branch is unmodeled.
        // Complexity review: four nested logarithmic partition searches over
        // source-partitioned u8 key slices, with no scan, sort, allocation, or
        // value clone; this retains the source equal_range lookup shape.
        let mut swap = false;
        let mut can_i_atom_type = i_atom_type;
        let mut can_k_atom_type = k_atom_type;
        let can_stretch_bend_type = stretch_bend_type;
        if i_atom_type > k_atom_type {
            can_i_atom_type = k_atom_type;
            can_k_atom_type = i_atom_type;
            swap = true;
        } else if i_atom_type == k_atom_type {
            swap = bond_type1 < bond_type2;
        }

        let j_start = self
            .d_j_atom_type
            .partition_point(|&stored| u32::from(stored) < j_atom_type);
        let j_end = self
            .d_j_atom_type
            .partition_point(|&stored| u32::from(stored) <= j_atom_type);
        if j_start == j_end {
            return (swap, None);
        }

        let i_keys = &self.d_i_atom_type[j_start..j_end];
        let i_start_in_range =
            i_keys.partition_point(|&stored| u32::from(stored) < can_i_atom_type);
        let i_end_in_range = i_keys.partition_point(|&stored| u32::from(stored) <= can_i_atom_type);
        if i_start_in_range == i_end_in_range {
            return (swap, None);
        }
        let i_start = j_start + i_start_in_range;
        let i_end = j_start + i_end_in_range;

        let k_keys = &self.d_k_atom_type[i_start..i_end];
        let k_start_in_range =
            k_keys.partition_point(|&stored| u32::from(stored) < can_k_atom_type);
        let k_end_in_range = k_keys.partition_point(|&stored| u32::from(stored) <= can_k_atom_type);
        if k_start_in_range == k_end_in_range {
            return (swap, None);
        }
        let k_start = i_start + k_start_in_range;
        let k_end = i_start + k_end_in_range;

        let stretch_bend_keys = &self.d_stretch_bend_type[k_start..k_end];
        let stretch_start =
            stretch_bend_keys.partition_point(|&stored| u32::from(stored) < can_stretch_bend_type);
        let stretch_end =
            stretch_bend_keys.partition_point(|&stored| u32::from(stored) <= can_stretch_bend_type);
        if stretch_start == stretch_end {
            return (swap, None);
        }

        let row_index = k_start + stretch_start;
        (swap, Some(&self.d_params[row_index]))
    }
}

const DEFAULT_MMFF_DFSB_TEXT: &str = include_str!("default_dfsb.tsv");

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MmffDfsbCollection {
    d_params: BTreeMap<u32, BTreeMap<u32, BTreeMap<u32, MmffStbn>>>,
}

impl MmffDfsbCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:4866-4896:
        // RDKit❗✔️: MMFFDfsbCollection::MMFFDfsbCollection(std::string mmffDfsb) {
        // RDKit❗✔️:   if (mmffDfsb.empty()) {
        // RDKit❗✔️:     mmffDfsb = defaultMMFFDfsb;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   std::istringstream inStream(mmffDfsb);
        // RDKit❗✔️:   std::string inLine = RDKit::getLine(inStream);
        // RDKit❗✔️:   while (!(inStream.eof())) {
        // RDKit❗✔️:     if (inLine[0] != '*') {
        // RDKit❗✔️:       MMFFStbn mmffStbnObj;
        // RDKit❗✔️:       boost::char_separator<char> tabSep("\t");
        // RDKit❗✔️:       tokenizer tokens(inLine, tabSep);
        // RDKit❗✔️:       tokenizer::iterator token = tokens.begin();
        // RDKit❗✔️:
        // RDKit❗✔️:       auto iAtomicNum = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗✔️:       ++token;
        // RDKit❗✔️:       auto jAtomicNum = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗✔️:       ++token;
        // RDKit❗✔️:       auto kAtomicNum = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗✔️:       ++token;
        // RDKit❗✔️:       mmffStbnObj.kbaIJK = boost::lexical_cast<double>(*token);
        // RDKit❗✔️:       ++token;
        // RDKit❗✔️:       mmffStbnObj.kbaKJI = boost::lexical_cast<double>(*token);
        // RDKit❗✔️:       ++token;
        // RDKit❗✔️:       d_params[iAtomicNum][jAtomicNum][kAtomicNum] = mmffStbnObj;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     inLine = RDKit::getLine(inStream);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior review: empty text selects this collection's committed
        // default asset; each processed row consumes full-width i/j/k followed
        // by kbaIJK/kbaKJI. Nested ordered-map assignment preserves the exact
        // source keys and makes the last duplicate complete key win. Reverse
        // endpoint canonicalization occurs only in get, never during parsing.
        // Source unchecked token dereferences become the existing typed
        // table/line/column errors for malformed rows and blank records.
        // Complexity review: one borrowed row/token pass; each unique nested
        // key uses ordered-map lookup/insertion, matching the source std::map
        // logarithmic levels. No key narrowing, sorting, flat scans, or value
        // copies are added; duplicate assignment replaces the existing leaf.
        let source = if text.is_empty() {
            DEFAULT_MMFF_DFSB_TEXT
        } else {
            text
        };

        let mut d_params = BTreeMap::<u32, BTreeMap<u32, BTreeMap<u32, MmffStbn>>>::new();

        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }
            if line.is_empty() {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(i_atomic_num_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(i_atomic_num) = parse_mmff_u32(i_atomic_num_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: i_atomic_num_cell.to_owned(),
                    },
                });
            };

            let Some(j_atomic_num_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(j_atomic_num) = parse_mmff_u32(j_atomic_num_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: j_atomic_num_cell.to_owned(),
                    },
                });
            };

            let Some(k_atomic_num_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(k_atomic_num) = parse_mmff_u32(k_atomic_num_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: k_atomic_num_cell.to_owned(),
                    },
                });
            };

            let Some(kba_ijk_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(kba_ijk) = parse_mmff_f64(kba_ijk_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: kba_ijk_cell.to_owned(),
                    },
                });
            };

            let Some(kba_kji_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(kba_kji) = parse_mmff_f64(kba_kji_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Dfsb,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: kba_kji_cell.to_owned(),
                    },
                });
            };

            let _ = d_params
                .entry(i_atomic_num)
                .or_default()
                .entry(j_atomic_num)
                .or_default()
                .insert(k_atomic_num, MmffStbn { kba_ijk, kba_kji });
        }

        Ok(Self { d_params })
    }

    pub(super) fn get(
        &self,
        periodic_table_row1: u32,
        periodic_table_row2: u32,
        periodic_table_row3: u32,
    ) -> (bool, Option<&MmffStbn>) {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:690-715:
        // RDKit❗✔️:   const std::pair<bool, const MMFFStbn *> getMMFFDfsbParams(
        // RDKit❗✔️:       const unsigned int periodicTableRow1,
        // RDKit❗✔️:       const unsigned int periodicTableRow2,
        // RDKit❗✔️:       const unsigned int periodicTableRow3) const {
        // RDKit❗✔️:     const MMFFStbn *mmffDfsbParams = nullptr;
        // RDKit❗✔️:     bool swap = false;
        // RDKit❗✔️:     unsigned int canPeriodicTableRow1 = periodicTableRow1;
        // RDKit❗✔️:     unsigned int canPeriodicTableRow3 = periodicTableRow3;
        // RDKit❗✔️:     if (periodicTableRow1 > periodicTableRow3) {
        // RDKit❗✔️:       canPeriodicTableRow1 = periodicTableRow3;
        // RDKit❗✔️:       canPeriodicTableRow3 = periodicTableRow1;
        // RDKit❗✔️:       swap = true;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     const auto res1 = d_params.find(canPeriodicTableRow1);
        // RDKit❗✔️:     if (res1 != d_params.end()) {
        // RDKit❗✔️:       const auto res2 = ((*res1).second).find(periodicTableRow2);
        // RDKit❗✔️:       if (res2 != ((*res1).second).end()) {
        // RDKit❗✔️:         const auto res3 = ((*res2).second).find(canPeriodicTableRow3);
        // RDKit❗✔️:         if (res3 != ((*res2).second).end()) {
        // RDKit❗✔️:           mmffDfsbParams = &((*res3).second);
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     return std::make_pair(swap, mmffDfsbParams);
        // RDKit❗✔️:   }
        // Behavior review: canonicalize only the first/third full-width query
        // rows when row1>row3, preserve the middle row and return the computed
        // swap flag even if any nested lookup misses. Constructor storage is
        // deliberately not canonicalized, so noncanonical input keys remain
        // observable only at their stored triple. Hits borrow the map leaf.
        // Complexity review: three ordered-map lookups are O(log n) at each
        // nested level, with no scan, allocation, value clone, or key rewrite;
        // this matches the source's three nested std::map searches.
        let (canonical_row1, canonical_row3, swap) = if periodic_table_row1 > periodic_table_row3 {
            (periodic_table_row3, periodic_table_row1, true)
        } else {
            (periodic_table_row1, periodic_table_row3, false)
        };

        let params = self
            .d_params
            .get(&canonical_row1)
            .and_then(|row2_map| row2_map.get(&periodic_table_row2))
            .and_then(|row3_map| row3_map.get(&canonical_row3));

        (swap, params)
    }
}

pub(super) fn default_mmff_stbn()
-> Result<&'static MmffStbnCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:73-76:
    // RDKit✔️✔️: const MMFFStbnCollection *getMMFFStbn() {
    // RDKit✔️✔️:   static MMFFStbnCollection ds_instance;
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // Behavior review: retain the Stbn default collection or its typed parse
    // error in a process-lifetime cell; the existing empty-text constructor
    // branch selects the exact committed source asset. Return a borrow without
    // allowing custom constructor input to initialize or replace this default.
    // Complexity review: source constructs its default once; this OnceLock
    // parses the asset once and warm calls borrow the stored result in O(1),
    // with no warm reparse or collection clone.
    static DEFAULT_MMFF_STBN: OnceLock<Result<MmffStbnCollection, MmffParamParseError>> =
        OnceLock::new();

    let result = DEFAULT_MMFF_STBN.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_STBN_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffStbnCollection::from_text("")
    });
    result.as_ref().map_err(|error| error)
}

pub(super) fn default_mmff_dfsb()
-> Result<&'static MmffDfsbCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:78-81:
    // RDKit✔️✔️: const MMFFDfsbCollection *getMMFFDfsb() {
    // RDKit✔️✔️:   static MMFFDfsbCollection ds_instance;
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // Behavior review: retain the Dfsb default map or its typed parse error in
    // its own process-lifetime cell; the existing empty-text constructor
    // branch selects the exact committed source asset. Return a borrow and do
    // not let another collection initialize or replace this default.
    // Complexity review: source constructs its default once; this OnceLock
    // parses the asset once and warm calls borrow the stored result in O(1),
    // with no warm reparse or collection clone.
    static DEFAULT_MMFF_DFSB: OnceLock<Result<MmffDfsbCollection, MmffParamParseError>> =
        OnceLock::new();

    let result = DEFAULT_MMFF_DFSB.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_DFSB_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffDfsbCollection::from_text("")
    });
    result.as_ref().map_err(|error| error)
}

pub(super) fn default_mmff_herschbach_laurie()
-> Result<&'static MmffHerschbachLaurieCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:58-61:
    // RDKit✔️✔️: const MMFFHerschbachLaurieCollection *getMMFFHerschbachLaurie() {
    // RDKit✔️✔️:   static MMFFHerschbachLaurieCollection ds_instance;
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // Retain this default's parsed collection or typed parse error independently
    // and return a borrow, matching the source function-local static lifetime.
    // The first call parses O(asset bytes); warm calls only access the OnceLock
    // and borrow the stored outcome without reparsing or cloning vectors.
    static DEFAULT_MMFF_HERSCHBACH_LAURIE: OnceLock<
        Result<MmffHerschbachLaurieCollection, MmffParamParseError>,
    > = OnceLock::new();

    let result = DEFAULT_MMFF_HERSCHBACH_LAURIE.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_HERSCHBACH_LAURIE_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffHerschbachLaurieCollection::from_text("")
    });
    result.as_ref().map_err(|error| error)
}

pub(super) fn default_mmff_cov_rad_pau_ele()
-> Result<&'static MmffCovRadPauEleCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:83-86:
    // RDKit✔️✔️: const MMFFCovRadPauEleCollection *getMMFFCovRadPauEle() {
    // RDKit✔️✔️:   static MMFFCovRadPauEleCollection ds_instance;
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // Retain this default's parsed collection or typed parse error independently
    // and return a borrow, matching the source function-local static lifetime.
    // The first call parses O(asset bytes); warm calls only access the OnceLock
    // and borrow the stored outcome without reparsing or cloning vectors.
    static DEFAULT_MMFF_COV_RAD_PAU_ELE: OnceLock<
        Result<MmffCovRadPauEleCollection, MmffParamParseError>,
    > = OnceLock::new();

    let result = DEFAULT_MMFF_COV_RAD_PAU_ELE.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_COV_RAD_PAU_ELE_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffCovRadPauEleCollection::from_text("")
    });
    result.as_ref().map_err(|error| error)
}

pub(super) fn default_mmff_pbci()
-> Result<&'static MmffPbciCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:63-66:
    // RDKit✔️✔️: const MMFFPBCICollection *getMMFFPBCI() {
    // RDKit✔️✔️:   static MMFFPBCICollection ds_instance;
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // Params.h:236 gives the collection constructor an empty default source;
    // from_text("") selects the exact embedded PBCI asset. Retain either the
    // parsed collection or its typed parse error and return a borrow to it.
    // One independent OnceLock initialization parses O(asset bytes); warm
    // calls do O(1) state checks without reparsing or cloning the row vector.
    static DEFAULT_MMFF_PBCI: OnceLock<Result<MmffPbciCollection, MmffParamParseError>> =
        OnceLock::new();

    let result = DEFAULT_MMFF_PBCI.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_PBCI_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffPbciCollection::from_text("")
    });
    result.as_ref().map_err(|error| error)
}

pub(super) fn default_mmff_angle()
-> Result<&'static MmffAngleCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:68-71:
    // RDKit✔️✔️: const MMFFAngleCollection *getMMFFAngle() {
    // RDKit✔️✔️:   static MMFFAngleCollection ds_instance;
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // Behavior review: retain this independent default's parsed collection or
    // typed parse error in one process-lifetime cell. The existing empty-text
    // constructor selects the exact Angle asset; custom input cannot initialize
    // or replace this default. Return a borrow of the stored outcome.
    // Complexity review: parse O(asset bytes) once; warm calls perform an O(1)
    // OnceLock lookup and borrow the value/error without cloning or reparsing.
    static DEFAULT_MMFF_ANGLE: OnceLock<Result<MmffAngleCollection, MmffParamParseError>> =
        OnceLock::new();

    let result = DEFAULT_MMFF_ANGLE.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_ANGLE_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffAngleCollection::from_text("")
    });
    result.as_ref().map_err(|error| error)
}

pub(super) fn default_mmff_oop(
    is_mmff_s: bool,
) -> Result<&'static MmffOopCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:94-98:
    // RDKit❗✔️: const MMFFOopCollection *getMMFFOop(const bool isMMFFs) {
    // RDKit❗✔️:   static MMFFOopCollection MMFF94(false, "");
    // RDKit❗✔️:   static MMFFOopCollection MMFF94s(true, "");
    // RDKit❗✔️:   return (isMMFFs) ? &MMFF94s : &MMFF94;
    // RDKit❗✔️: }
    // Behavior review: initialize and retain the regular collection outcome
    // first, then the MMFFs collection outcome, before selecting either by
    // flag. Each cell preserves its typed parse error for a borrowed return;
    // custom Oop construction remains independent of these process defaults.
    // Complexity review: the first call parses both Oop assets once; warm
    // calls perform two O(1) OnceLock lookups and return a borrow without
    // cloning or reparsing either collection.
    static DEFAULT_MMFF_OOP: OnceLock<Result<MmffOopCollection, MmffParamParseError>> =
        OnceLock::new();
    static DEFAULT_MMFF_OOP_S: OnceLock<Result<MmffOopCollection, MmffParamParseError>> =
        OnceLock::new();

    let regular = DEFAULT_MMFF_OOP.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_OOP_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffOopCollection::from_text(false, "")
    });
    let mmff_s = DEFAULT_MMFF_OOP_S.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_OOP_S_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffOopCollection::from_text(true, "")
    });

    let result = if is_mmff_s { mmff_s } else { regular };
    result.as_ref().map_err(|error| error)
}

pub(super) fn default_mmff_tor(
    is_mmff_s: bool,
) -> Result<&'static MmffTorCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:88-92:
    // RDKit❗✔️: const MMFFTorCollection *getMMFFTor(const bool isMMFFs) {
    // RDKit❗✔️:   static MMFFTorCollection MMFF94(false, "");
    // RDKit❗✔️:   static MMFFTorCollection MMFF94s(true, "");
    // RDKit❗✔️:   return (isMMFFs) ? &MMFF94s : &MMFF94;
    // RDKit❗✔️: }
    // Behavior review: construct and retain the regular collection first and
    // the MMFFs collection second before selecting either by the flag. Keep
    // each typed parse error in its own OnceLock and return a borrow to that
    // stored outcome. Custom Tor constructors do not interact with these
    // process defaults.
    // Complexity review: the first call parses both fixed assets once; warm
    // calls perform two O(1) OnceLock lookups and return a borrow without
    // cloning or reparsing either collection.
    static DEFAULT_MMFF_TOR: OnceLock<Result<MmffTorCollection, MmffParamParseError>> =
        OnceLock::new();
    static DEFAULT_MMFF_TOR_S: OnceLock<Result<MmffTorCollection, MmffParamParseError>> =
        OnceLock::new();

    let regular = DEFAULT_MMFF_TOR.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_TOR_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffTorCollection::from_text(false, "")
    });
    let mmff_s = DEFAULT_MMFF_TOR_S.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_TOR_S_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffTorCollection::from_text(true, "")
    });

    let result = if is_mmff_s { mmff_s } else { regular };
    result.as_ref().map_err(|error| error)
}

pub(super) fn default_mmff_chg() -> Result<&'static MmffChgCollection, &'static MmffParamParseError>
{
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:48-51:
    // RDKit✔️✔️: const MMFFChgCollection *getMMFFChg() {
    // RDKit✔️✔️:   static MMFFChgCollection ds_instance;
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // Params.h:305 gives the collection constructor an empty default source;
    // from_text("") selects this collection's exact committed Chg asset.
    // Retain either the parsed collection or its typed parse error and return
    // a borrow to that stored outcome.
    // One independent OnceLock initialization parses O(asset bytes); warm
    // calls do O(1) state checks without reparsing or cloning the row vectors.
    static DEFAULT_MMFF_CHG: OnceLock<Result<MmffChgCollection, MmffParamParseError>> =
        OnceLock::new();

    let result = DEFAULT_MMFF_CHG.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_CHG_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffChgCollection::from_text("")
    });
    result.as_ref().map_err(|error| error)
}

pub(super) fn default_mmff_def() -> Result<&'static MmffDefCollection, &'static MmffParamParseError>
{
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:53-56:
    // RDKit✔️✔️: const MMFFDefCollection *getMMFFDef() {
    // RDKit✔️✔️:   static MMFFDefCollection ds_instance;
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // The default argument in Params.h:182 is empty; the existing
    // from_text constructor selects the exact embedded default asset.
    // Behavior review: retain either the one owned default collection or its
    // typed construction error and return a borrow to that stored outcome.
    // Complexity review: OnceLock initializes once and subsequent calls do
    // O(1) state checks without reparsing or cloning the row vector.
    static DEFAULT_MMFF_DEF: OnceLock<Result<MmffDefCollection, MmffParamParseError>> =
        OnceLock::new();

    let result = DEFAULT_MMFF_DEF.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_DEF_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffDefCollection::from_text("")
    });
    result.as_ref().map_err(|error| error)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct MmffProp {
    pub(super) atno: u8,
    pub(super) crd: u8,
    pub(super) val: u8,
    pub(super) pilp: u8,
    pub(super) mltb: u8,
    pub(super) arom: u8,
    pub(super) linh: u8,
    pub(super) sbmb: u8,
}

#[derive(Debug)]
pub(super) struct MmffPropCollection {
    d_i_atom_type: Vec<u8>,
    d_params: Vec<MmffProp>,
}

impl MmffPropCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:427-479:
        // RDKit✔️🔝: MMFFPropCollection::MMFFPropCollection(std::string mmffProp) {
        // RDKit✔️🔝:   if (mmffProp.empty()) {
        // RDKit✔️🔝:     mmffProp = defaultMMFFProp;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   std::istringstream inStream(mmffProp);
        // RDKit✔️🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   while (!(inStream.eof())) {
        // RDKit❗✔️:     if (inLine[0] != '*') {
        // RDKit✔️🔝:       MMFFProp mmffPropObj;
        // RDKit✔️🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit✔️🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit✔️🔝:       tokenizer::iterator token = tokens.begin();
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int atomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit✔️🔝: #else
        // RDKit✔️🔝:       d_iAtomType.push_back(
        // RDKit✔️🔝:           (std::uint8_t)boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝: #endif
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       mmffPropObj.atno =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       mmffPropObj.crd =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       mmffPropObj.val =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       mmffPropObj.pilp =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       mmffPropObj.mltb =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       mmffPropObj.arom =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       mmffPropObj.linh =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit✔️🔝:       mmffPropObj.sbmb =
        // RDKit✔️🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token));
        // RDKit✔️🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[atomType] = mmffPropObj;
        // RDKit✔️🔝: #else
        // RDKit✔️🔝:       d_params.push_back(mmffPropObj);
        // RDKit✔️🔝: #endif
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     inLine = RDKit::getLine(inStream);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // Behavior review: the vector build selects the embedded default only
        // for empty input, skips byte-zero-star records, narrows the key and
        // each of the eight fields to u8, and appends every row without sorting
        // or duplicate removal. Missing/invalid tokens are returned as the
        // packet's typed safety errors where C++ dereferences/lexes undefined
        // input; this does not claim a successful source result for those rows.
        // Complexity review: both output Vecs grow in source order while lines
        // and tab fields remain borrowed; shared helpers scan each input byte
        // once with O(1) parser state and no per-row token staging allocation.
        let source = if text.is_empty() {
            include_str!("default_prop.tsv")
        } else {
            text
        };

        let mut d_i_atom_type = Vec::new();
        let mut d_params = Vec::new();
        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Prop,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            };
            let Some(atom_type_value) = parse_mmff_u32(atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Prop,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: atom_type_cell.to_owned(),
                    },
                });
            };
            d_i_atom_type.push(atom_type_value as u8);

            let mut column = 1;
            let mut next_u8 = || -> Result<u8, MmffParamParseError> {
                let current_column = column;
                column += 1;
                let Some(cell) = tokens.next() else {
                    return Err(MmffParamParseError {
                        table: MmffParamTable::Prop,
                        line: physical_line,
                        column: current_column,
                        cause: MmffParamParseCause::MissingToken,
                    });
                };
                let Some(value) = parse_mmff_u32(cell) else {
                    return Err(MmffParamParseError {
                        table: MmffParamTable::Prop,
                        line: physical_line,
                        column: current_column,
                        cause: MmffParamParseCause::InvalidUnsigned {
                            cell: cell.to_owned(),
                        },
                    });
                };
                Ok(value as u8)
            };

            let atno = next_u8()?;
            let crd = next_u8()?;
            let val = next_u8()?;
            let pilp = next_u8()?;
            let mltb = next_u8()?;
            let arom = next_u8()?;
            let linh = next_u8()?;
            let sbmb = next_u8()?;
            d_params.push(MmffProp {
                atno,
                crd,
                val,
                pilp,
                mltb,
                arom,
                linh,
                sbmb,
            });
        }

        Ok(Self {
            d_i_atom_type,
            d_params,
        })
    }

    pub(super) fn get(&self, atom_type: u32) -> Option<&MmffProp> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:197-208:
        // RDKit✔️✔️:   const MMFFProp *operator()(const unsigned int atomType) const {
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res = d_params.find(atomType);
        // RDKit❌❌:     return ((res != d_params.end()) ? &((*res).second) : NULL);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto bounds =
        // RDKit✔️✔️:         std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), atomType);
        // RDKit✔️✔️:     return ((bounds.first != bounds.second)
        // RDKit✔️✔️:                 ? &d_params[bounds.first - d_iAtomType.begin()]
        // RDKit✔️✔️:                 : nullptr);
        // RDKit❌❌: #endif
        // RDKit✔️✔️:   }
        // Behavior review: lower-bound search over the preserved byte keys
        // selects the same first duplicate as `bounds.first`; each key is
        // widened to u32 before comparison, so the query is never narrowed.
        // As in source `equal_range`, lookup requires keys partitioned for the
        // query and does not sort or repair custom input.
        // Complexity review: one random-access binary lower-bound search is
        // O(log N), followed by fixed comparisons and one borrowed row index;
        // no temporary collection, allocation, or linear scan is added.
        let first = self
            .d_i_atom_type
            .partition_point(|&stored| u32::from(stored) < atom_type);
        let stored = self.d_i_atom_type.get(first)?;
        if u32::from(*stored) != atom_type {
            return None;
        }

        Some(&self.d_params[first])
    }
}

pub(super) fn default_mmff_prop()
-> Result<&'static MmffPropCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:28-31:
    // RDKit✔️✔️: const MMFFPropCollection *getMMFFProp() {
    // RDKit✔️✔️:   static MMFFPropCollection ds_instance("");
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // Behavior review: the explicit empty source selects the Prop asset in
    // the existing constructor and its independently stored Result is borrowed.
    // Complexity review: one initialization parses O(asset bytes) into the
    // existing parallel vectors; later calls perform O(1) state checks and
    // never allocate, reparse or copy either vector.
    static DEFAULT_MMFF_PROP: OnceLock<Result<MmffPropCollection, MmffParamParseError>> =
        OnceLock::new();

    let result = DEFAULT_MMFF_PROP.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_PROP_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffPropCollection::from_text("")
    });
    result.as_ref().map_err(|error| error)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MmffBond {
    pub(super) kb: f64,
    pub(super) r0: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MmffBondCollection {
    d_bond_type: Vec<u8>,
    d_i_atom_type: Vec<u8>,
    d_j_atom_type: Vec<u8>,
    d_params: Vec<MmffBond>,
}

impl MmffBondCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:1292-1338:
        // RDKit❗🔝: MMFFBondCollection::MMFFBondCollection(std::string mmffBond) {
        // RDKit❗🔝:   if (mmffBond.empty()) {
        // RDKit❗🔝:     mmffBond = defaultMMFFBond;
        // RDKit❗🔝:   }
        // RDKit❗🔝:   std::istringstream inStream(mmffBond);
        // RDKit❗🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   while (!(inStream.eof())) {
        // RDKit❗✔️:     if (inLine[0] != '*') {
        // RDKit❗🔝:       MMFFBond mmffBondObj;
        // RDKit❗🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit❗🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit❗🔝:       tokenizer::iterator token = tokens.begin();
        //
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int bondType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_bondType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int atomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_iAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_jAtomType.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffBondObj.kb = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffBondObj.r0 = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[bondType][atomType][jAtomType] = mmffBondObj;
        // RDKit❗🔝: #else
        // RDKit❗🔝:       d_params.push_back(mmffBondObj);
        // RDKit❗🔝: #endif
        // RDKit❗🔝:     }
        // RDKit❗🔝:     inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   }
        // RDKit❗🔝: }
        // The source vector branch preserves row order and narrows all three
        // stored keys after unsigned parsing. Existing borrowed line/token
        // and numeric owners keep valid cells allocation-free; this also
        // avoids the source stream's owned line/token copies. Empty or
        // malformed consumed rows return the established typed safety errors
        // instead of claiming parity for source undefined behavior.
        let source = if text.is_empty() {
            include_str!("default_bond.tsv")
        } else {
            text
        };

        let mut d_bond_type = Vec::new();
        let mut d_i_atom_type = Vec::new();
        let mut d_j_atom_type = Vec::new();
        let mut d_params = Vec::new();

        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(bond_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bond,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            };
            let Some(bond_type) = parse_mmff_u32(bond_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bond,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: bond_type_cell.to_owned(),
                    },
                });
            };

            let Some(i_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bond,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(i_atom_type) = parse_mmff_u32(i_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bond,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: i_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(j_atom_type_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bond,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(j_atom_type) = parse_mmff_u32(j_atom_type_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bond,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: j_atom_type_cell.to_owned(),
                    },
                });
            };

            let Some(kb_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bond,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(kb) = parse_mmff_f64(kb_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bond,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: kb_cell.to_owned(),
                    },
                });
            };

            let Some(r0_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bond,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(r0) = parse_mmff_f64(r0_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bond,
                    line: physical_line,
                    column: 4,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: r0_cell.to_owned(),
                    },
                });
            };

            d_bond_type.push(bond_type as u8);
            d_i_atom_type.push(i_atom_type as u8);
            d_j_atom_type.push(j_atom_type as u8);
            d_params.push(MmffBond { kb, r0 });
        }

        Ok(Self {
            d_bond_type,
            d_i_atom_type,
            d_j_atom_type,
            d_params,
        })
    }

    pub(super) fn get(
        &self,
        bond_type: u32,
        atom_type: u32,
        nbr_atom_type: u32,
    ) -> Option<&MmffBond> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:315-360:
        // RDKit❗✔️:   const MMFFBond *operator()(const unsigned int bondType,
        // RDKit❗✔️:                              const unsigned int atomType,
        // RDKit❗✔️:                              const unsigned int nbrAtomType) const {
        // RDKit❗✔️:     const MMFFBond *mmffBondParams = nullptr;
        // RDKit❗✔️:     unsigned int canAtomType = atomType;
        // RDKit❗✔️:     unsigned int canNbrAtomType = nbrAtomType;
        // RDKit❗✔️:     if (atomType > nbrAtomType) {
        // RDKit❗✔️:       canAtomType = nbrAtomType;
        // RDKit❗✔️:       canNbrAtomType = atomType;
        // RDKit❗✔️:     }
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
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
        // RDKit❗✔️: #else
        // RDKit❗✔️:     auto bounds =
        // RDKit❗✔️:         std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), canAtomType);
        // RDKit❗✔️:     if (bounds.first != bounds.second) {
        // RDKit❗✔️:       bounds = std::equal_range(
        // RDKit❗✔️:           d_jAtomType.begin() + (bounds.first - d_iAtomType.begin()),
        // RDKit❗✔️:           d_jAtomType.begin() + (bounds.second - d_iAtomType.begin()),
        // RDKit❗✔️:           canNbrAtomType);
        // RDKit❗✔️:       if (bounds.first != bounds.second) {
        // RDKit❗✔️:         bounds = std::equal_range(
        // RDKit❗✔️:             d_bondType.begin() + (bounds.first - d_jAtomType.begin()),
        // RDKit❗✔️:             d_bondType.begin() + (bounds.second - d_jAtomType.begin()),
        // RDKit❗✔️:             bondType);
        // RDKit❗✔️:         if (bounds.first != bounds.second) {
        // RDKit❗✔️:           mmffBondParams = &d_params[bounds.first - d_bondType.begin()];
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❌❌: #endif
        //
        // RDKit❗✔️:     return mmffBondParams;
        // RDKit❗✔️:   }
        // Behavior review: preserve the full-width query, swap only the
        // endpoint pair, and search the nested source-sorted i/j/bondType
        // slices in order. `partition_point` supplies the same lower and upper
        // bounds as each vector `equal_range`, so the first duplicate row is
        // returned. Complexity review: three O(log N) random-access searches
        // within progressively narrower ranges; no lookup allocation, key
        // narrowing, sort, or row clone is introduced.
        let (canonical_i_atom_type, canonical_j_atom_type) = if atom_type > nbr_atom_type {
            (nbr_atom_type, atom_type)
        } else {
            (atom_type, nbr_atom_type)
        };

        let i_start = self
            .d_i_atom_type
            .partition_point(|&stored| u32::from(stored) < canonical_i_atom_type);
        let i_end = self
            .d_i_atom_type
            .partition_point(|&stored| u32::from(stored) <= canonical_i_atom_type);
        if i_start == i_end {
            return None;
        }

        let j_range = &self.d_j_atom_type[i_start..i_end];
        let j_start_in_range =
            j_range.partition_point(|&stored| u32::from(stored) < canonical_j_atom_type);
        let j_end_in_range =
            j_range.partition_point(|&stored| u32::from(stored) <= canonical_j_atom_type);
        if j_start_in_range == j_end_in_range {
            return None;
        }
        let j_start = i_start + j_start_in_range;
        let j_end = i_start + j_end_in_range;

        let bond_range = &self.d_bond_type[j_start..j_end];
        let bond_start_in_range =
            bond_range.partition_point(|&stored| u32::from(stored) < bond_type);
        let bond_end_in_range =
            bond_range.partition_point(|&stored| u32::from(stored) <= bond_type);
        if bond_start_in_range == bond_end_in_range {
            return None;
        }

        Some(&self.d_params[j_start + bond_start_in_range])
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MmffBndkCollection {
    d_i_atomic_num: Vec<u8>,
    d_j_atomic_num: Vec<u8>,
    d_params: Vec<MmffBond>,
}

impl MmffBndkCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffParamParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:1850-1889:
        // RDKit❗🔝: MMFFBndkCollection::MMFFBndkCollection(std::string mmffBndk) {
        // RDKit❗🔝:   if (mmffBndk.empty()) {
        // RDKit❗🔝:     mmffBndk = defaultMMFFBndk;
        // RDKit❗🔝:   }
        // RDKit❗🔝:   std::istringstream inStream(mmffBndk);
        // RDKit❗🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   while (!(inStream.eof())) {
        // RDKit❗✔️:     if (inLine[0] != '*') {
        // RDKit❗🔝:       MMFFBond mmffBondObj;
        // RDKit❗🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit❗🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit❗🔝:       tokenizer::iterator token = tokens.begin();
        //
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int iAtomicNum = boost::lexical_cast<unsigned int>(*token);
        // RDKit❌❌: #else
        // RDKit❗🔝:       d_iAtomicNum.push_back(
        // RDKit❗🔝:           (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❌❌: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       unsigned int jAtomicNum = boost::lexical_cast<unsigned int>(*token);
        // RDKit❌❌: #else
        // RDKit❗🔝:       d_jAtomicNum.push_back(
        // RDKit❗🔝:           (std::uint8_t)boost::lexical_cast<unsigned int>(*token));
        // RDKit❌❌: #endif
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffBondObj.r0 = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❗🔝:       mmffBondObj.kb = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:       ++token;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:       d_params[iAtomicNum][jAtomicNum] = mmffBondObj;
        // RDKit❌❌: #else
        // RDKit❗🔝:       d_params.push_back(mmffBondObj);
        // RDKit❌❌: #endif
        // RDKit❗🔝:     }
        // RDKit❗🔝:     inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   }
        // RDKit❗🔝: }
        // Behavior review: preserve source row order and u8 key casts. The
        // source field order is i/j/r0/kb, so the shared value stores r0 and
        // kb in their named fields without changing either parsed float.
        // Malformed consumed cells return the existing typed input-safety
        // errors instead of claiming parity for unchecked token dereferences.
        // Complexity review: one forward borrowed-record/token pass and three
        // row vectors; there is no sorting, deduplication, or token Vec, and
        // only invalid-token errors allocate cell text.
        let source = if text.is_empty() {
            include_str!("default_bndk.tsv")
        } else {
            text
        };

        let mut d_i_atomic_num = Vec::new();
        let mut d_j_atomic_num = Vec::new();
        let mut d_params = Vec::new();

        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }

            let mut tokens = tokenize_mmff_line(line);
            let Some(i_atomic_num_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bndk,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                });
            };
            let Some(i_atomic_num) = parse_mmff_u32(i_atomic_num_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bndk,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: i_atomic_num_cell.to_owned(),
                    },
                });
            };

            let Some(j_atomic_num_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bndk,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(j_atomic_num) = parse_mmff_u32(j_atomic_num_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bndk,
                    line: physical_line,
                    column: 1,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: j_atomic_num_cell.to_owned(),
                    },
                });
            };

            let Some(r0_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bndk,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(r0) = parse_mmff_f64(r0_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bndk,
                    line: physical_line,
                    column: 2,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: r0_cell.to_owned(),
                    },
                });
            };

            let Some(kb_cell) = tokens.next() else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bndk,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::MissingToken,
                });
            };
            let Some(kb) = parse_mmff_f64(kb_cell) else {
                return Err(MmffParamParseError {
                    table: MmffParamTable::Bndk,
                    line: physical_line,
                    column: 3,
                    cause: MmffParamParseCause::InvalidFloat {
                        cell: kb_cell.to_owned(),
                    },
                });
            };

            d_i_atomic_num.push(i_atomic_num as u8);
            d_j_atomic_num.push(j_atomic_num as u8);
            d_params.push(MmffBond { kb, r0 });
        }

        Ok(Self {
            d_i_atomic_num,
            d_j_atomic_num,
            d_params,
        })
    }

    pub(super) fn get(&self, atomic_num: i32, nbr_atomic_num: i32) -> Option<&MmffBond> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:381-414:
        // RDKit❗✔️:   const MMFFBond *operator()(const int atomicNum,
        // RDKit❗✔️:                              const int nbrAtomicNum) const {
        // RDKit❗✔️:     const MMFFBond *mmffBndkParams = nullptr;
        // RDKit❗✔️:     unsigned int canAtomicNum = atomicNum;
        // RDKit❗✔️:     unsigned int canNbrAtomicNum = nbrAtomicNum;
        // RDKit❗✔️:     if (atomicNum > nbrAtomicNum) {
        // RDKit❗✔️:       canAtomicNum = nbrAtomicNum;
        // RDKit❗✔️:       canNbrAtomicNum = atomicNum;
        // RDKit❗✔️:     }
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res1 = d_params.find(canAtomicNum);
        // RDKit❌❌:     std::map<const unsigned int, MMFFBond>::const_iterator res2;
        // RDKit❌❌:     if (res1 != d_params.end()) {
        // RDKit❌❌:       res2 = ((*res1).second).find(canNbrAtomicNum);
        // RDKit❌❌:       if (res2 != ((*res1).second).end()) {
        // RDKit❌❌:         mmffBndkParams = &((*res2).second);
        // RDKit❌❌:       }
        // RDKit❌❌:     }
        // RDKit❌❌: #else
        // RDKit❗✔️:     auto bounds = std::equal_range(d_iAtomicNum.begin(), d_iAtomicNum.end(),
        // RDKit❗✔️:                                    canAtomicNum);
        // RDKit❗✔️:     if (bounds.first != bounds.second) {
        // RDKit❗✔️:       bounds = std::equal_range(
        // RDKit❗✔️:           d_jAtomicNum.begin() + (bounds.first - d_iAtomicNum.begin()),
        // RDKit❗✔️:           d_jAtomicNum.begin() + (bounds.second - d_iAtomicNum.begin()),
        // RDKit❗✔️:           canNbrAtomicNum);
        // RDKit❗✔️:       if (bounds.first != bounds.second) {
        // RDKit❗✔️:         mmffBndkParams = &d_params[bounds.first - d_jAtomicNum.begin()];
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❌❌: #endif
        //
        // RDKit❗✔️:     return mmffBndkParams;
        // RDKit❗✔️:   }
        // Behavior review: make modulo-2^32 copies before testing the signed
        // original arguments, then compare the unswapped or swapped full-width
        // copies with the stored u8 keys widened to u32. Nested source-sorted
        // ranges return the first duplicate. Complexity review: two O(log N)
        // searches on borrowed slices; the lookup allocates and clones nothing.
        let mut can_atomic_num = atomic_num as u32;
        let mut can_nbr_atomic_num = nbr_atomic_num as u32;
        if atomic_num > nbr_atomic_num {
            std::mem::swap(&mut can_atomic_num, &mut can_nbr_atomic_num);
        }

        let i_start = self
            .d_i_atomic_num
            .partition_point(|&stored| u32::from(stored) < can_atomic_num);
        let i_end = self
            .d_i_atomic_num
            .partition_point(|&stored| u32::from(stored) <= can_atomic_num);
        if i_start == i_end {
            return None;
        }

        let j_range = &self.d_j_atomic_num[i_start..i_end];
        let j_start_in_range =
            j_range.partition_point(|&stored| u32::from(stored) < can_nbr_atomic_num);
        let j_end_in_range =
            j_range.partition_point(|&stored| u32::from(stored) <= can_nbr_atomic_num);
        if j_start_in_range == j_end_in_range {
            return None;
        }

        Some(&self.d_params[i_start + j_start_in_range])
    }
}

pub(super) fn default_mmff_bond()
-> Result<&'static MmffBondCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:43-46:
    // RDKit✔️✔️: const MMFFBondCollection *getMMFFBond() {
    // RDKit✔️✔️:   static MMFFBondCollection ds_instance;
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // Params.h:362 supplies the empty constructor default; from_text("") uses
    // the exact committed Bond asset. Behavior review: preserve one independent
    // cached construction result and return a borrow of its value or error.
    // Complexity review: first use parses O(asset bytes); warm calls perform
    // O(1) lock access and borrow without allocating, reparsing, or cloning.
    static DEFAULT_MMFF_BOND: OnceLock<Result<MmffBondCollection, MmffParamParseError>> =
        OnceLock::new();

    let result = DEFAULT_MMFF_BOND.get_or_init(|| MmffBondCollection::from_text(""));
    result.as_ref().map_err(|error| error)
}

pub(super) fn default_mmff_bndk()
-> Result<&'static MmffBndkCollection, &'static MmffParamParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:38-41:
    // RDKit✔️✔️: const MMFFBndkCollection *getMMFFBndk() {
    // RDKit✔️✔️:   static MMFFBndkCollection ds_instance;
    // RDKit✔️✔️:   return &ds_instance;
    // RDKit✔️✔️: }
    // Params.h:416 supplies the empty constructor default; from_text("") uses
    // the exact committed Bndk asset. Behavior review: preserve a second,
    // independent cached result and return a borrow of its value or error.
    // Complexity review: first use parses O(asset bytes); warm calls perform
    // O(1) lock access and borrow without allocating, reparsing, or cloning.
    static DEFAULT_MMFF_BNDK: OnceLock<Result<MmffBndkCollection, MmffParamParseError>> =
        OnceLock::new();

    let result = DEFAULT_MMFF_BNDK.get_or_init(|| MmffBndkCollection::from_text(""));
    result.as_ref().map_err(|error| error)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct MmffVdw {
    pub(super) alpha_i: f64,
    pub(super) n_i: f64,
    pub(super) a_i: f64,
    pub(super) g_i: f64,
    pub(super) r_star: f64,
    pub(super) da: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MmffVdwCollection {
    power: f64,
    pub(super) b: f64,
    pub(super) beta: f64,
    pub(super) darad: f64,
    pub(super) daeps: f64,
    d_atom_type: Vec<u8>,
    d_params: Vec<MmffVdw>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum MmffVdwParseError {
    Table(MmffParamParseError),
    MissingHeader,
}

impl std::fmt::Display for MmffVdwParseError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Table(error) => std::fmt::Display::fmt(error, formatter),
            Self::MissingHeader => {
                formatter.write_str("MMFF VdW source has no processed constants record")
            }
        }
    }
}

impl std::error::Error for MmffVdwParseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Table(error) => Some(error),
            Self::MissingHeader => None,
        }
    }
}

impl MmffVdwCollection {
    pub(super) fn from_text(text: &str) -> Result<Self, MmffVdwParseError> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.cpp:8246-8300:
        // RDKit❗🔝: MMFFVdWCollection::MMFFVdWCollection(std::string mmffVdW) {
        // RDKit❗🔝:   if (mmffVdW.empty()) {
        // RDKit❗🔝:     mmffVdW = defaultMMFFVdW;
        // RDKit❗🔝:   }
        // RDKit❗🔝:   std::istringstream inStream(mmffVdW);
        // RDKit❗🔝:   bool firstLine = true;
        // RDKit❗🔝:   std::string inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   while (!(inStream.eof())) {
        // RDKit❗🔝:     if (inLine[0] != '*') {
        // RDKit❗🔝:       boost::char_separator<char> tabSep("\t");
        // RDKit❗🔝:       tokenizer tokens(inLine, tabSep);
        // RDKit❗🔝:       tokenizer::iterator token = tokens.begin();
        // RDKit❗🔝:       if (firstLine) {
        // RDKit❗🔝:         firstLine = false;
        // RDKit❗🔝:         this->power = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:         this->B = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:         this->Beta = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:         this->DARAD = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:         this->DAEPS = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:       } else {
        // RDKit❗✔️:         MMFFVdW mmffVdWObj;
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:         unsigned int atomType = boost::lexical_cast<unsigned int>(*token);
        // RDKit❗🔝: #else
        // RDKit❗🔝:         d_atomType.push_back(
        // RDKit❗🔝:             (std::uint8_t)(boost::lexical_cast<unsigned int>(*token)));
        // RDKit❗🔝: #endif
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:         mmffVdWObj.alpha_i = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:         mmffVdWObj.N_i = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:         mmffVdWObj.A_i = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:         mmffVdWObj.G_i = boost::lexical_cast<double>(*token);
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:         mmffVdWObj.DA = (boost::lexical_cast<std::string>(*token)).at(0);
        // RDKit❗🔝:         ++token;
        // RDKit❗🔝:         mmffVdWObj.R_star =
        // RDKit❗🔝:             mmffVdWObj.A_i * pow(mmffVdWObj.alpha_i, this->power);
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:         d_params[atomType] = mmffVdWObj;
        // RDKit❗🔝: #else
        // RDKit❗🔝:         d_params.push_back(mmffVdWObj);
        // RDKit❗🔝: #endif
        // RDKit❗🔝:       }
        // RDKit❗🔝:     }
        // RDKit❗🔝:     inLine = RDKit::getLine(inStream);
        // RDKit❗🔝:   }
        // RDKit❗🔝: }
        // Behavior review: borrow the source text and reuse the existing line,
        // tab, unsigned and binary64 helpers; consume the header once, then
        // append each narrowed key and row in source order. The explicit
        // MissingHeader and typed errors replace source-undefined access for
        // absent or malformed fields. Complexity review: the borrowed scan
        // avoids source stream/string/token copies while retaining O(input)
        // construction and amortized Vec appends; no row sorting or map path
        // is introduced. The source-vector path is the only modeled build.
        let source = if text.is_empty() {
            include_str!("default_vdw.tsv")
        } else {
            text
        };

        let mut header: Option<(f64, f64, f64, f64, f64)> = None;
        let mut d_atom_type = Vec::new();
        let mut d_params = Vec::new();

        for (line_index, line) in source_lines(source).enumerate() {
            let physical_line = line_index + 1;
            if line.as_bytes().first() == Some(&b'*') {
                continue;
            }
            if line.is_empty() {
                return Err(MmffVdwParseError::Table(MmffParamParseError {
                    table: MmffParamTable::Vdw,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                }));
            }

            let mut tokens = tokenize_mmff_line(line);
            if header.is_none() {
                let (power, b, beta, darad, daeps) = {
                    let mut column = 0;
                    let mut next_float = || -> Result<f64, MmffVdwParseError> {
                        let current_column = column;
                        column += 1;
                        let Some(cell) = tokens.next() else {
                            return Err(MmffVdwParseError::Table(MmffParamParseError {
                                table: MmffParamTable::Vdw,
                                line: physical_line,
                                column: current_column,
                                cause: MmffParamParseCause::MissingToken,
                            }));
                        };
                        let Some(value) = parse_mmff_f64(cell) else {
                            return Err(MmffVdwParseError::Table(MmffParamParseError {
                                table: MmffParamTable::Vdw,
                                line: physical_line,
                                column: current_column,
                                cause: MmffParamParseCause::InvalidFloat {
                                    cell: cell.to_owned(),
                                },
                            }));
                        };
                        Ok(value)
                    };
                    (
                        next_float()?,
                        next_float()?,
                        next_float()?,
                        next_float()?,
                        next_float()?,
                    )
                };
                header = Some((power, b, beta, darad, daeps));
                continue;
            }

            let Some(atom_type_cell) = tokens.next() else {
                return Err(MmffVdwParseError::Table(MmffParamParseError {
                    table: MmffParamTable::Vdw,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::MissingToken,
                }));
            };
            let Some(atom_type) = parse_mmff_u32(atom_type_cell) else {
                return Err(MmffVdwParseError::Table(MmffParamParseError {
                    table: MmffParamTable::Vdw,
                    line: physical_line,
                    column: 0,
                    cause: MmffParamParseCause::InvalidUnsigned {
                        cell: atom_type_cell.to_owned(),
                    },
                }));
            };
            let atom_type = atom_type as u8;
            d_atom_type.push(atom_type);

            let (alpha_i, n_i, a_i, g_i) = {
                let mut column = 1;
                let mut next_float = || -> Result<f64, MmffVdwParseError> {
                    let current_column = column;
                    column += 1;
                    let Some(cell) = tokens.next() else {
                        return Err(MmffVdwParseError::Table(MmffParamParseError {
                            table: MmffParamTable::Vdw,
                            line: physical_line,
                            column: current_column,
                            cause: MmffParamParseCause::MissingToken,
                        }));
                    };
                    let Some(value) = parse_mmff_f64(cell) else {
                        return Err(MmffVdwParseError::Table(MmffParamParseError {
                            table: MmffParamTable::Vdw,
                            line: physical_line,
                            column: current_column,
                            cause: MmffParamParseCause::InvalidFloat {
                                cell: cell.to_owned(),
                            },
                        }));
                    };
                    Ok(value)
                };
                (next_float()?, next_float()?, next_float()?, next_float()?)
            };
            let Some(da_cell) = tokens.next() else {
                return Err(MmffVdwParseError::Table(MmffParamParseError {
                    table: MmffParamTable::Vdw,
                    line: physical_line,
                    column: 5,
                    cause: MmffParamParseCause::MissingToken,
                }));
            };
            let Some(&da) = da_cell.as_bytes().first() else {
                return Err(MmffVdwParseError::Table(MmffParamParseError {
                    table: MmffParamTable::Vdw,
                    line: physical_line,
                    column: 5,
                    cause: MmffParamParseCause::MissingToken,
                }));
            };
            let Some((power, _, _, _, _)) = header else {
                return Err(MmffVdwParseError::MissingHeader);
            };
            let r_star = a_i * alpha_i.powf(power);
            d_params.push(MmffVdw {
                alpha_i,
                n_i,
                a_i,
                g_i,
                r_star,
                da,
            });
        }

        let Some((power, b, beta, darad, daeps)) = header else {
            return Err(MmffVdwParseError::MissingHeader);
        };
        Ok(Self {
            power,
            b,
            beta,
            darad,
            daeps,
            d_atom_type,
            d_params,
        })
    }

    pub(super) fn get(&self, atom_type: u32) -> Option<&MmffVdw> {
        // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
        // Code/ForceField/MMFF/Params.h:977-990:
        // RDKit✔️✔️:   const MMFFVdW *operator()(const unsigned int atomType) const {
        // RDKit❌❌: #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        // RDKit❌❌:     const auto res = d_params.find(atomType);
        // RDKit❌❌:     return (res != d_params.end() ? &((*res).second) : NULL);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:     auto bounds =
        // RDKit✔️✔️:         std::equal_range(d_atomType.begin(), d_atomType.end(), atomType);
        // RDKit✔️✔️:     return ((bounds.first != bounds.second)
        // RDKit✔️✔️:                 ? &d_params[bounds.first - d_atomType.begin()]
        // RDKit✔️✔️:                 : nullptr);
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:   }
        // The pinned vector build uses equal_range and returns its first match.
        // The key vector's source-sorted precondition is retained; the full
        // u32 target is compared without narrowing. partition_point gives the
        // same lower bound in O(log n) time without allocation or row scans.
        let index = self
            .d_atom_type
            .partition_point(|&stored| u32::from(stored) < atom_type);
        let stored = self.d_atom_type.get(index)?;
        if u32::from(*stored) != atom_type {
            return None;
        }
        self.d_params.get(index)
    }
}

pub(super) fn default_mmff_vdw() -> Result<&'static MmffVdwCollection, &'static MmffVdwParseError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:100-103:
    // RDKit❗✔️: const MMFFVdWCollection *getMMFFVdW() {
    // RDKit❗✔️:   static MMFFVdWCollection ds_instance;
    // RDKit❗✔️:   return &ds_instance;
    // RDKit❗✔️: }
    // Behavior review: retain the one parsed embedded collection or its typed
    // safety error and return a borrow to the stored result. Complexity
    // review: OnceLock initializes once, then warm calls perform O(1) state
    // checks without reparsing or cloning either vector.
    static DEFAULT_MMFF_VDW: OnceLock<Result<MmffVdwCollection, MmffVdwParseError>> =
        OnceLock::new();

    let result = DEFAULT_MMFF_VDW.get_or_init(|| {
        #[cfg(test)]
        DEFAULT_MMFF_VDW_CONSTRUCTIONS.fetch_add(1, Ordering::Relaxed);
        MmffVdwCollection::from_text("")
    });
    result.as_ref().map_err(|error| error)
}

#[cfg(test)]
mod tests {
    use super::{
        DEFAULT_MMFF_ANGLE_CONSTRUCTIONS, DEFAULT_MMFF_CHG_CONSTRUCTIONS,
        DEFAULT_MMFF_COV_RAD_PAU_ELE_CONSTRUCTIONS, DEFAULT_MMFF_DEF_CONSTRUCTIONS,
        DEFAULT_MMFF_DFSB_CONSTRUCTIONS, DEFAULT_MMFF_HERSCHBACH_LAURIE_CONSTRUCTIONS,
        DEFAULT_MMFF_OOP_CONSTRUCTIONS, DEFAULT_MMFF_OOP_S_CONSTRUCTIONS,
        DEFAULT_MMFF_PBCI_CONSTRUCTIONS, DEFAULT_MMFF_PROP_CONSTRUCTIONS,
        DEFAULT_MMFF_STBN_CONSTRUCTIONS, DEFAULT_MMFF_TOR_CONSTRUCTIONS,
        DEFAULT_MMFF_TOR_S_CONSTRUCTIONS, DEFAULT_MMFF_VDW_CONSTRUCTIONS, MMFFAromCollection,
        MmffAngle, MmffAngleCollection, MmffAngleLookupError, MmffBndkCollection,
        MmffBondCollection, MmffChgCollection, MmffCovRadPauEle, MmffCovRadPauEleCollection,
        MmffDef, MmffDefCollection, MmffDfsbCollection, MmffHerschbachLaurie,
        MmffHerschbachLaurieCollection, MmffOop, MmffOopCollection, MmffOopLookupError,
        MmffParamParseCause, MmffParamParseError, MmffParamTable, MmffPbciCollection, MmffProp,
        MmffPropCollection, MmffStbn, MmffStbnCollection, MmffTor, MmffTorCollection,
        MmffTorLookupError, MmffVdw, MmffVdwCollection, MmffVdwParseError, default_mmff_angle,
        default_mmff_bndk, default_mmff_bond, default_mmff_chg, default_mmff_cov_rad_pau_ele,
        default_mmff_def, default_mmff_dfsb, default_mmff_herschbach_laurie, default_mmff_oop,
        default_mmff_pbci, default_mmff_prop, default_mmff_stbn, default_mmff_tor,
        default_mmff_vdw,
    };
    use std::sync::atomic::Ordering;

    fn fnv1a64(bytes: &[u8]) -> u64 {
        let mut hash = 0xcbf29ce484222325_u64;
        for &byte in bytes {
            hash = (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3);
        }
        hash
    }

    fn sha256_for_fixed_asset_test(bytes: &[u8]) -> [u8; 32] {
        const INITIAL: [u32; 8] = [
            0x6a09_e667,
            0xbb67_ae85,
            0x3c6e_f372,
            0xa54f_f53a,
            0x510e_527f,
            0x9b05_688c,
            0x1f83_d9ab,
            0x5be0_cd19,
        ];
        const ROUND: [u32; 64] = [
            0x428a_2f98,
            0x7137_4491,
            0xb5c0_fbcf,
            0xe9b5_dba5,
            0x3956_c25b,
            0x59f1_11f1,
            0x923f_82a4,
            0xab1c_5ed5,
            0xd807_aa98,
            0x1283_5b01,
            0x2431_85be,
            0x550c_7dc3,
            0x72be_5d74,
            0x80de_b1fe,
            0x9bdc_06a7,
            0xc19b_f174,
            0xe49b_69c1,
            0xefbe_4786,
            0x0fc1_9dc6,
            0x240c_a1cc,
            0x2de9_2c6f,
            0x4a74_84aa,
            0x5cb0_a9dc,
            0x76f9_88da,
            0x983e_5152,
            0xa831_c66d,
            0xb003_27c8,
            0xbf59_7fc7,
            0xc6e0_0bf3,
            0xd5a7_9147,
            0x06ca_6351,
            0x1429_2967,
            0x27b7_0a85,
            0x2e1b_2138,
            0x4d2c_6dfc,
            0x5338_0d13,
            0x650a_7354,
            0x766a_0abb,
            0x81c2_c92e,
            0x9272_2c85,
            0xa2bf_e8a1,
            0xa81a_664b,
            0xc24b_8b70,
            0xc76c_51a3,
            0xd192_e819,
            0xd699_0624,
            0xf40e_3585,
            0x106a_a070,
            0x19a4_c116,
            0x1e37_6c08,
            0x2748_774c,
            0x34b0_bcb5,
            0x391c_0cb3,
            0x4ed8_aa4a,
            0x5b9c_ca4f,
            0x682e_6ff3,
            0x748f_82ee,
            0x78a5_636f,
            0x84c8_7814,
            0x8cc7_0208,
            0x90be_fffa,
            0xa450_6ceb,
            0xbef9_a3f7,
            0xc671_78f2,
        ];

        let bit_len = (bytes.len() as u64).wrapping_mul(8);
        let mut padded = bytes.to_vec();
        padded.push(0x80);
        while padded.len() % 64 != 56 {
            padded.push(0);
        }
        padded.extend_from_slice(&bit_len.to_be_bytes());

        let mut state = INITIAL;
        for block in padded.chunks_exact(64) {
            let mut words = [0_u32; 64];
            for (index, chunk) in block.chunks_exact(4).enumerate() {
                words[index] = u32::from_be_bytes(chunk.try_into().expect("four-byte word"));
            }
            for index in 16..64 {
                let x = words[index - 15];
                let sigma0 = x.rotate_right(7) ^ x.rotate_right(18) ^ (x >> 3);
                let y = words[index - 2];
                let sigma1 = y.rotate_right(17) ^ y.rotate_right(19) ^ (y >> 10);
                words[index] = words[index - 16]
                    .wrapping_add(sigma0)
                    .wrapping_add(words[index - 7])
                    .wrapping_add(sigma1);
            }

            let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = state;
            for index in 0..64 {
                let sum1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
                let choose = (e & f) ^ (!e & g);
                let temp1 = h
                    .wrapping_add(sum1)
                    .wrapping_add(choose)
                    .wrapping_add(ROUND[index])
                    .wrapping_add(words[index]);
                let sum0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
                let majority = (a & b) ^ (a & c) ^ (b & c);
                let temp2 = sum0.wrapping_add(majority);

                h = g;
                g = f;
                f = e;
                e = d.wrapping_add(temp1);
                d = c;
                c = b;
                b = a;
                a = temp1.wrapping_add(temp2);
            }

            for (word, value) in state.iter_mut().zip([a, b, c, d, e, f, g, h]) {
                *word = word.wrapping_add(value);
            }
        }

        let mut digest = [0_u8; 32];
        for (chunk, word) in digest.chunks_exact_mut(4).zip(state) {
            chunk.copy_from_slice(&word.to_be_bytes());
        }
        digest
    }

    const EXPECTED_DEFAULT_MMFF_AROMATIC_TYPES: [u8; 17] = [
        37, 38, 39, 44, 58, 59, 63, 64, 65, 66, 69, 76, 78, 79, 80, 81, 82,
    ];

    #[test]
    fn mmff_arom_constructor_and_unsigned_lookup_match_source_rows() {
        let default = MMFFAromCollection::new(None);
        let empty = MMFFAromCollection::new(Some(&[]));
        let mut custom_source = vec![0, 37, 37, 255];
        let custom = MMFFAromCollection::new(Some(&custom_source));
        custom_source[0] = 1;
        custom_source.push(42);

        assert_eq!(default.d_params, EXPECTED_DEFAULT_MMFF_AROMATIC_TYPES);
        assert!(empty.d_params.is_empty());
        assert_eq!(custom.d_params, [0, 37, 37, 255]);
        assert_eq!(custom_source, [1, 37, 37, 255, 42]);

        for atom_type in 0_u32..=u8::MAX.into() {
            let expected_default = EXPECTED_DEFAULT_MMFF_AROMATIC_TYPES
                .iter()
                .any(|&stored| u32::from(stored) == atom_type);
            let expected_custom = [0_u8, 37, 37, 255]
                .iter()
                .any(|&stored| u32::from(stored) == atom_type);

            assert_eq!(default.is_mmff_aromatic(atom_type), expected_default);
            assert!(!empty.is_mmff_aromatic(atom_type));
            assert_eq!(custom.is_mmff_aromatic(atom_type), expected_custom);
        }

        for (atom_type, expected_custom) in [(255_u32, true), (256, false), (u32::MAX, false)] {
            let expected_default = EXPECTED_DEFAULT_MMFF_AROMATIC_TYPES
                .iter()
                .any(|&stored| u32::from(stored) == atom_type);

            assert_eq!(default.is_mmff_aromatic(atom_type), expected_default);
            assert!(!empty.is_mmff_aromatic(atom_type));
            assert_eq!(custom.is_mmff_aromatic(atom_type), expected_custom);
        }
    }

    #[test]
    fn mmff_dp_m01_assets_match_frozen_bytes_hashes_and_notices() {
        const DEFAULT_MMFF_DEF: &[u8] = include_bytes!("default_def.tsv");
        const DEFAULT_MMFF_PROP: &[u8] = include_bytes!("default_prop.tsv");
        const PARAMETER_PROVENANCE: &str = include_str!("parameter-provenance.md");

        const MERCK_NOTICE: &[u8] =
            b"*          Copyright (c) Merck and Co., Inc., 1994, 1995, 1996\n*                         All Rights Reserved\n";

        assert_eq!(DEFAULT_MMFF_DEF.len(), 6775);
        assert_eq!(fnv1a64(DEFAULT_MMFF_DEF), 0xea4ce6723b17849c);
        assert!(
            DEFAULT_MMFF_DEF
                .windows(MERCK_NOTICE.len())
                .any(|window| window == MERCK_NOTICE)
        );

        assert_eq!(DEFAULT_MMFF_PROP.len(), 2027);
        assert_eq!(fnv1a64(DEFAULT_MMFF_PROP), 0xadcc42a1ca3b802b);
        assert!(
            DEFAULT_MMFF_PROP
                .windows(MERCK_NOTICE.len())
                .any(|window| window == MERCK_NOTICE)
        );

        assert!(PARAMETER_PROVENANCE.contains("third_party/rdkit/license.txt"));
        assert!(PARAMETER_PROVENANCE.contains("Copyright Kevlin Henney"));
        assert!(PARAMETER_PROVENANCE.contains("Copyright Alexander Nasonov, 2006-2010."));
        assert!(PARAMETER_PROVENANCE.contains("Copyright Antony Polukhin, 2011-2024."));
        assert!(PARAMETER_PROVENANCE.contains("Copyright John R. Bandela 2001."));
        assert!(PARAMETER_PROVENANCE.contains("Boost Software License - Version 1.0"));
    }

    fn format_def_row(label: &str, fields: [&str; 5], separator: &str) -> String {
        format!(
            "{label}{separator}{}{separator}{}{separator}{}{separator}{}{separator}{}",
            fields[0], fields[1], fields[2], fields[3], fields[4]
        )
    }

    fn format_prop_row(key: &str, fields: [&str; 8], separator: &str) -> String {
        format!(
            "{key}{separator}{}{separator}{}{separator}{}{separator}{}{separator}{}{separator}{}{separator}{}{separator}{}",
            fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], fields[7]
        )
    }

    fn prop(fields: [u8; 8]) -> MmffProp {
        MmffProp {
            atno: fields[0],
            crd: fields[1],
            val: fields[2],
            pilp: fields[3],
            mltb: fields[4],
            arom: fields[5],
            linh: fields[6],
            sbmb: fields[7],
        }
    }

    #[test]
    fn mmff_dp_m05_def_constructor_matches_source_rows_and_typed_edges() {
        let source_rows: [(&str, u32, [u8; 4]); 6] = [
            ("ZERO", 0, [0, 0, 0, 0]),
            ("TYPE1_FIRST", 1, [1, 2, 3, 4]),
            ("TYPE1_DUPLICATE", 1, [5, 6, 7, 8]),
            ("TYPE2", 2, [9, 10, 11, 12]),
            ("TYPE1_LATER", 1, [13, 14, 15, 16]),
            ("TYPE257", 257, [17, 18, 19, 20]),
        ];
        let expected = [
            MmffDef {
                eq_level: [1, 2, 3, 4],
            },
            MmffDef {
                eq_level: [9, 10, 11, 12],
            },
            MmffDef {
                eq_level: [13, 14, 15, 16],
            },
        ];
        let mut profile_calls = 0;
        for newline_profile in 0..3 {
            for separator in ["\t", "\t\t"] {
                for has_comment in [false, true] {
                    let mut lines =
                        Vec::with_capacity(source_rows.len() + if has_comment { 1 } else { 0 });
                    if has_comment {
                        lines.push("* fixed source comment".to_owned());
                    }
                    for &(label, atom_type, levels) in &source_rows {
                        lines.push(format!(
                            "{label}{separator}{atom_type}{separator}{}{separator}{}{separator}{}{separator}{}",
                            levels[0], levels[1], levels[2], levels[3]
                        ));
                    }
                    let source = match newline_profile {
                        0 => format!("{}\n", lines.join("\n")),
                        1 => format!("{}\r\n", lines.join("\r\n")),
                        _ => lines.join("\n"),
                    };

                    assert!(profile_calls < 12, "unexpected extra profile call");
                    let actual = MmffDefCollection::from_text(&source)
                        .unwrap_or_else(|error| panic!("profile failed: {error}"));
                    assert_eq!(actual.d_params.as_slice(), expected.as_slice());
                    profile_calls += 1;
                }
            }
        }
        assert_eq!(profile_calls, 12);

        let prefix = "ZERO\t0\t0\t0\t0\t0\nTYPE1\t1\t1\t2\t3\t4\nTYPE2\t2\t9\t10\t11\t12\n";
        let terminated = format!("{prefix}TYPE3\t3\t21\t22\t23\t24\n");
        let unterminated = format!("{prefix}TYPE3\t3\t21\t22\t23\t24");
        let eof_expected = [
            MmffDef {
                eq_level: [1, 2, 3, 4],
            },
            MmffDef {
                eq_level: [9, 10, 11, 12],
            },
            MmffDef {
                eq_level: [21, 22, 23, 24],
            },
        ];
        let mut eof_calls = 0;
        assert!(eof_calls < 2);
        let actual = MmffDefCollection::from_text(&terminated).unwrap();
        assert_eq!(actual.d_params.as_slice(), eof_expected.as_slice());
        eof_calls += 1;
        assert!(eof_calls < 2);
        let actual = MmffDefCollection::from_text(&unterminated).unwrap();
        assert_eq!(actual.d_params.as_slice(), &eof_expected[..2]);
        eof_calls += 1;
        assert_eq!(eof_calls, 2);

        let roles: [(&str, [&str; 5], &str, usize); 3] = [
            ("ordinary", ["2", "1", "2", "3", "4"], "", 1),
            ("initial-zero", ["0", "0", "0", "0", "0"], "", 1),
            (
                "duplicate-type-one",
                ["1", "5", "6", "7", "8"],
                "FIRST\t1\t1\t2\t3\t4\n",
                2,
            ),
        ];
        let mut invalid_calls = 0;
        for (role, base_fields, prefix, error_line) in roles {
            for column in 1..=5 {
                let mut fields = base_fields;
                fields[column - 1] = "1x";
                let row = format_def_row(role, fields, "\t");
                let source = format!("{prefix}{row}\n");

                assert!(invalid_calls < 15, "unexpected extra invalid-row call");
                let error = MmffDefCollection::from_text(&source).unwrap_err();
                assert_eq!(
                    error,
                    MmffParamParseError {
                        table: MmffParamTable::Def,
                        line: error_line,
                        column,
                        cause: MmffParamParseCause::InvalidUnsigned {
                            cell: "1x".to_owned(),
                        },
                    },
                    "role {role}, numeric column {column}"
                );
                invalid_calls += 1;
            }
        }
        assert_eq!(invalid_calls, 15);

        let mut empty_calls = 0;
        assert!(empty_calls < 1);
        let error = MmffDefCollection::from_text("\n").unwrap_err();
        assert_eq!(
            error,
            MmffParamParseError {
                table: MmffParamTable::Def,
                line: 1,
                column: 0,
                cause: MmffParamParseCause::EmptyProcessedLine,
            }
        );
        empty_calls += 1;
        assert_eq!(empty_calls, 1);

        let missing_cases: [(&str, usize); 6] = [
            ("\t\n", 0),
            ("ROW\n", 1),
            ("ROW\t1\n", 2),
            ("ROW\t1\t2\n", 3),
            ("ROW\t1\t2\t3\n", 4),
            ("ROW\t1\t2\t3\t4\n", 5),
        ];
        let mut missing_calls = 0;
        for (source, present_tokens) in missing_cases {
            assert!(missing_calls < 6, "unexpected extra missing-token call");
            let error = MmffDefCollection::from_text(source).unwrap_err();
            let (column, cause) = if present_tokens == 0 {
                (0, MmffParamParseCause::EmptyProcessedLine)
            } else {
                (present_tokens, MmffParamParseCause::MissingToken)
            };
            assert_eq!(
                error,
                MmffParamParseError {
                    table: MmffParamTable::Def,
                    line: 1,
                    column,
                    cause,
                },
                "present nonempty tokens {present_tokens}"
            );
            missing_calls += 1;
        }
        assert_eq!(missing_calls, 6);

        assert_eq!(
            profile_calls + eof_calls + invalid_calls + empty_calls + missing_calls,
            36
        );
    }

    #[test]
    fn mmff_dp_m06_def_positional_lookup_matches_source_boundaries() {
        const CUSTOM_ROWS: [MmffDef; 3] = [
            MmffDef {
                eq_level: [1, 2, 3, 4],
            },
            MmffDef {
                eq_level: [9, 10, 11, 12],
            },
            MmffDef {
                eq_level: [13, 14, 15, 16],
            },
        ];
        const GAP_ROWS: [MmffDef; 2] =
            [MmffDef { eq_level: [8; 4] }, MmffDef { eq_level: [12; 4] }];

        let comment_only_input = "* fixed comment-only collection\n".to_owned();
        let custom_input = concat!(
            "ZERO\t0\t0\t0\t0\t0\n",
            "TYPE1_FIRST\t1\t1\t2\t3\t4\n",
            "TYPE1_DUPLICATE\t1\t5\t6\t7\t8\n",
            "TYPE2\t2\t9\t10\t11\t12\n",
            "TYPE1_LATER\t1\t13\t14\t15\t16\n",
            "TYPE257\t257\t17\t18\t19\t20\n",
        )
        .to_owned();
        let gap_input = concat!("GAP8\t8\t8\t8\t8\t8\n", "GAP12\t12\t12\t12\t12\t12\n",).to_owned();
        let inputs: [&str; 4] = [
            "",
            comment_only_input.as_str(),
            custom_input.as_str(),
            gap_input.as_str(),
        ];
        let input_snapshots: [String; 4] = std::array::from_fn(|index| inputs[index].to_owned());
        let input_addresses = inputs.map(|input| input.as_ptr());

        let default = MmffDefCollection::from_text(inputs[0]).unwrap();
        let comment_only = MmffDefCollection::from_text(inputs[1]).unwrap();
        let custom = MmffDefCollection::from_text(inputs[2]).unwrap();
        let gaps = MmffDefCollection::from_text(inputs[3]).unwrap();
        let collections = [&default, &comment_only, &custom, &gaps];
        let snapshots: [Vec<MmffDef>; 4] =
            std::array::from_fn(|index| collections[index].d_params.clone());

        assert_eq!(default.d_params.len(), 95);
        assert!(comment_only.d_params.is_empty());
        assert_eq!(custom.d_params.as_slice(), CUSTOM_ROWS.as_slice());
        assert_eq!(gaps.d_params.as_slice(), GAP_ROWS.as_slice());
        assert_eq!(default.d_params[82].eq_level, [87; 4]);
        assert_eq!(default.d_params[94].eq_level, [99; 4]);

        let default_field_hash = |rows: &[MmffDef]| {
            let fields: Vec<u8> = rows.iter().flat_map(|row| row.eq_level).collect();
            assert_eq!(fields.len(), 380);
            fnv1a64(&fields)
        };
        assert_eq!(default_field_hash(&default.d_params), 0x3e49_ed6b_8e84_98d5);

        let mut matrix_calls = 0;
        for atom_type in (0_u32..=256).chain([u32::MAX]) {
            for collection_index in 0..collections.len() {
                assert!(matrix_calls < 1032, "unexpected extra Def lookup");
                let collection = collections[collection_index];
                assert_eq!(collection.d_params.as_slice(), snapshots[collection_index]);
                assert_eq!(
                    inputs[collection_index].as_ptr(),
                    input_addresses[collection_index]
                );
                assert_eq!(inputs[collection_index], input_snapshots[collection_index]);

                let expected_index = usize::try_from(atom_type)
                    .ok()
                    .and_then(|position| position.checked_sub(1))
                    .filter(|&index| index < snapshots[collection_index].len());
                let actual = collection.get(atom_type);
                matrix_calls += 1;

                match (actual, expected_index) {
                    (Some(actual), Some(index)) => {
                        assert_eq!(actual, &snapshots[collection_index][index]);
                        assert!(std::ptr::eq(actual, &collection.d_params[index]));
                    }
                    (None, None) => {}
                    (actual, expected) => panic!(
                        "collection {collection_index}, query {atom_type}: \
                         actual {actual:?}, expected position {expected:?}"
                    ),
                }

                if collection_index == 0 && atom_type == 83 {
                    assert_eq!(actual.map(|row| row.eq_level), Some([87; 4]));
                }
                if collection_index == 0 && atom_type == 95 {
                    assert_eq!(actual.map(|row| row.eq_level), Some([99; 4]));
                }
                if collection_index == 2 && (1..=3).contains(&atom_type) {
                    assert_eq!(
                        actual.map(|row| row.eq_level),
                        Some(CUSTOM_ROWS[atom_type as usize - 1].eq_level)
                    );
                }
                if collection_index == 3 && (1..=2).contains(&atom_type) {
                    assert_eq!(
                        actual.map(|row| row.eq_level),
                        Some(GAP_ROWS[atom_type as usize - 1].eq_level)
                    );
                }

                assert_eq!(collection.d_params.as_slice(), snapshots[collection_index]);
                assert_eq!(
                    inputs[collection_index].as_ptr(),
                    input_addresses[collection_index]
                );
                assert_eq!(inputs[collection_index], input_snapshots[collection_index]);
            }
        }
        assert_eq!(matrix_calls, 1032);
        assert_eq!(default_field_hash(&default.d_params), 0x3e49_ed6b_8e84_98d5);

        let identity_cases: [(usize, u32, Option<[u8; 4]>); 9] = [
            (0, 83, Some([87; 4])),
            (0, 95, Some([99; 4])),
            (1, 1, None),
            (2, 1, Some([1, 2, 3, 4])),
            (2, 2, Some([9, 10, 11, 12])),
            (2, 3, Some([13, 14, 15, 16])),
            (3, 1, Some([8; 4])),
            (3, 2, Some([12; 4])),
            (3, 12, None),
        ];
        let mut repeated_calls = 0;
        for (collection_index, atom_type, expected_fields) in identity_cases {
            let collection = collections[collection_index];
            assert_eq!(collection.d_params.as_slice(), snapshots[collection_index]);
            assert_eq!(
                inputs[collection_index].as_ptr(),
                input_addresses[collection_index]
            );
            assert_eq!(inputs[collection_index], input_snapshots[collection_index]);

            assert!(repeated_calls < 18);
            let first = collection.get(atom_type);
            repeated_calls += 1;
            assert!(repeated_calls < 18);
            let second = collection.get(atom_type);
            repeated_calls += 1;
            match (first, second, expected_fields) {
                (Some(first), Some(second), Some(expected)) => {
                    assert_eq!(first.eq_level, expected);
                    assert!(std::ptr::eq(first, second));
                }
                (None, None, None) => {}
                (first, second, expected) => {
                    panic!("repeated query {atom_type}: {first:?}, {second:?}, {expected:?}")
                }
            }

            assert_eq!(collection.d_params.as_slice(), snapshots[collection_index]);
            assert_eq!(
                inputs[collection_index].as_ptr(),
                input_addresses[collection_index]
            );
            assert_eq!(inputs[collection_index], input_snapshots[collection_index]);
        }
        assert_eq!(repeated_calls, 18);
    }

    #[test]
    fn mmff_dp_m07_prop_constructor_preserves_parallel_rows_and_typed_edges() {
        const KEYS: [u8; 4] = [0, 1, 1, 255];
        const EXPECTED_PROPS: [MmffProp; 4] = [
            MmffProp {
                atno: 6,
                crd: 4,
                val: 4,
                pilp: 0,
                mltb: 0,
                arom: 0,
                linh: 0,
                sbmb: 2,
            },
            MmffProp {
                atno: 8,
                crd: 1,
                val: 12,
                pilp: 1,
                mltb: 1,
                arom: 0,
                linh: 0,
                sbmb: 2,
            },
            MmffProp {
                atno: 7,
                crd: 3,
                val: 34,
                pilp: 0,
                mltb: 1,
                arom: 0,
                linh: 0,
                sbmb: 2,
            },
            MmffProp {
                atno: 26,
                crd: 0,
                val: 0,
                pilp: 0,
                mltb: 0,
                arom: 0,
                linh: 0,
                sbmb: 2,
            },
        ];
        let profile_rows: [(&str, [&str; 8]); 4] = [
            ("0", ["6", "4", "4", "0", "0", "0", "0", "2"]),
            ("1", ["8", "1", "12", "1", "1", "0", "0", "2"]),
            ("1", ["7", "3", "34", "0", "1", "0", "0", "2"]),
            ("255", ["26", "0", "0", "0", "0", "0", "0", "2"]),
        ];

        let mut profile_calls = 0;
        for newline_profile in 0..3 {
            for separator in ["\t", "\t\t"] {
                for has_comment in [false, true] {
                    let mut lines = Vec::with_capacity(5);
                    if has_comment {
                        lines.push("* fixed byte-zero comment".to_owned());
                    }
                    for &(key, fields) in &profile_rows {
                        lines.push(format_prop_row(key, fields, separator));
                    }
                    let source = match newline_profile {
                        0 => format!("{}\n", lines.join("\n")),
                        1 => format!("{}\r\n", lines.join("\r\n")),
                        _ => lines.join("\n"),
                    };
                    let source_snapshot = source.clone();
                    let source_address = source.as_ptr();
                    let expected_len = if newline_profile == 2 { 3 } else { 4 };

                    assert!(profile_calls < 12, "unexpected extra Prop profile");
                    let actual = MmffPropCollection::from_text(&source)
                        .unwrap_or_else(|error| panic!("Prop profile failed: {error}"));
                    profile_calls += 1;

                    assert_eq!(actual.d_i_atom_type.as_slice(), &KEYS[..expected_len]);
                    assert_eq!(actual.d_params.as_slice(), &EXPECTED_PROPS[..expected_len]);
                    assert_eq!(actual.d_i_atom_type.len(), actual.d_params.len());
                    assert!(actual.d_params.iter().all(|row| row.sbmb == 2));
                    assert_eq!(source.as_ptr(), source_address);
                    assert_eq!(source, source_snapshot);
                }
            }
        }
        assert_eq!(profile_calls, 12);

        const BASE_FIELDS: [u32; 9] = [1, 10, 20, 30, 40, 50, 60, 70, 80];
        let mut conversion_calls = 0;
        let mut valid_conversion_calls = 0;
        let mut invalid_conversion_calls = 0;
        for field_index in 0..9 {
            for replacement in ["256", "-1", "1x"] {
                let mut cells = BASE_FIELDS.map(|value| value.to_string());
                cells[field_index] = replacement.to_owned();
                let mut source = cells.join("\t");
                source.push('\n');
                let source_snapshot = source.clone();
                let source_address = source.as_ptr();

                assert!(conversion_calls < 27, "unexpected extra Prop field case");
                let actual = MmffPropCollection::from_text(&source);
                conversion_calls += 1;
                if replacement == "1x" {
                    invalid_conversion_calls += 1;
                    assert_eq!(
                        actual.unwrap_err(),
                        MmffParamParseError {
                            table: MmffParamTable::Prop,
                            line: 1,
                            column: field_index,
                            cause: MmffParamParseCause::InvalidUnsigned {
                                cell: "1x".to_owned(),
                            },
                        },
                        "numeric token column {field_index}"
                    );
                } else {
                    valid_conversion_calls += 1;
                    let actual = actual.unwrap();
                    let narrowed = if replacement == "256" { 0 } else { 255 };
                    let expected_key = if field_index == 0 {
                        narrowed
                    } else {
                        BASE_FIELDS[0] as u8
                    };
                    let mut expected_fields = [10_u8, 20, 30, 40, 50, 60, 70, 80];
                    if field_index != 0 {
                        expected_fields[field_index - 1] = narrowed;
                    }
                    assert_eq!(actual.d_i_atom_type.as_slice(), &[expected_key]);
                    assert_eq!(actual.d_params.as_slice(), &[prop(expected_fields)]);
                    assert_eq!(actual.d_i_atom_type.len(), actual.d_params.len());
                }
                assert_eq!(source.as_ptr(), source_address);
                assert_eq!(source, source_snapshot);
            }
        }
        assert_eq!(conversion_calls, 27);
        assert_eq!(valid_conversion_calls, 18);
        assert_eq!(invalid_conversion_calls, 9);

        let mut missing_calls = 0;
        for present_tokens in 0..=8 {
            let mut source = if present_tokens == 0 {
                "\t".to_owned()
            } else {
                (1..=present_tokens)
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>()
                    .join("\t")
            };
            source.push('\n');
            let source_snapshot = source.clone();
            let source_address = source.as_ptr();

            assert!(
                missing_calls < 9,
                "unexpected extra Prop missing-token case"
            );
            let error = MmffPropCollection::from_text(&source).unwrap_err();
            missing_calls += 1;
            let (column, cause) = if present_tokens == 0 {
                (0, MmffParamParseCause::EmptyProcessedLine)
            } else {
                (present_tokens, MmffParamParseCause::MissingToken)
            };
            assert_eq!(
                error,
                MmffParamParseError {
                    table: MmffParamTable::Prop,
                    line: 1,
                    column,
                    cause,
                },
                "present nonempty numeric tokens {present_tokens}"
            );
            assert_eq!(source.as_ptr(), source_address);
            assert_eq!(source, source_snapshot);
        }
        assert_eq!(missing_calls, 9);

        let blankline = "\n".to_owned();
        let blankline_snapshot = blankline.clone();
        let blankline_address = blankline.as_ptr();
        assert!(missing_calls + 1 <= 10);
        assert_eq!(
            MmffPropCollection::from_text(&blankline).unwrap_err(),
            MmffParamParseError {
                table: MmffParamTable::Prop,
                line: 1,
                column: 0,
                cause: MmffParamParseCause::EmptyProcessedLine,
            }
        );
        assert_eq!(blankline.as_ptr(), blankline_address);
        assert_eq!(blankline, blankline_snapshot);

        let unsorted_source = format!(
            "{}\n{}\n",
            format_prop_row("2", ["1", "2", "3", "4", "5", "6", "7", "8"], "\t"),
            format_prop_row("1", ["9", "10", "11", "12", "13", "14", "15", "16"], "\t")
        );
        let unsorted_snapshot = unsorted_source.clone();
        let unsorted_address = unsorted_source.as_ptr();
        let mut unsorted_calls = 0;
        assert!(unsorted_calls < 1);
        let unsorted = MmffPropCollection::from_text(&unsorted_source).unwrap();
        unsorted_calls += 1;
        assert_eq!(unsorted.d_i_atom_type, [2, 1]);
        assert_eq!(
            unsorted.d_params,
            [
                prop([1, 2, 3, 4, 5, 6, 7, 8]),
                prop([9, 10, 11, 12, 13, 14, 15, 16]),
            ]
        );
        assert_eq!(unsorted.d_i_atom_type.len(), unsorted.d_params.len());
        assert_eq!(unsorted_source.as_ptr(), unsorted_address);
        assert_eq!(unsorted_source, unsorted_snapshot);
        assert_eq!(unsorted_calls, 1);

        assert_eq!(
            profile_calls + conversion_calls + missing_calls + 1 + unsorted_calls,
            50
        );
    }

    #[test]
    fn mmff_dp_m08_prop_lookup_matches_source_rows_and_boundaries() {
        let custom_rows = [
            prop([6, 4, 4, 0, 0, 0, 0, 2]),
            prop([8, 1, 12, 1, 1, 0, 0, 2]),
            prop([7, 3, 34, 0, 1, 0, 0, 2]),
            prop([26, 0, 0, 0, 0, 0, 0, 2]),
        ];
        assert_ne!(custom_rows[1], custom_rows[2]);
        let custom_source = format!(
            "{}\n{}\n{}\n{}\n",
            format_prop_row("0", ["6", "4", "4", "0", "0", "0", "0", "2"], "\t"),
            format_prop_row("1", ["8", "1", "12", "1", "1", "0", "0", "2"], "\t"),
            format_prop_row("1", ["7", "3", "34", "0", "1", "0", "0", "2"], "\t"),
            format_prop_row("255", ["26", "0", "0", "0", "0", "0", "0", "2"], "\t"),
        );
        let gap_rows = [
            prop([3, 5, 7, 9, 11, 13, 15, 17]),
            prop([18, 16, 14, 12, 10, 8, 6, 4]),
        ];
        let gap_source = format!(
            "{}\n{}\n",
            format_prop_row("2", ["3", "5", "7", "9", "11", "13", "15", "17"], "\t"),
            format_prop_row("8", ["18", "16", "14", "12", "10", "8", "6", "4"], "\t"),
        );

        let collections = [
            MmffPropCollection::from_text("").unwrap(),
            MmffPropCollection::from_text("* fixed comment-only table\n").unwrap(),
            MmffPropCollection::from_text(&custom_source).unwrap(),
            MmffPropCollection::from_text(&gap_source).unwrap(),
        ];
        let snapshots: Vec<(Vec<u8>, Vec<MmffProp>)> = collections
            .iter()
            .map(|collection| {
                (
                    collection.d_i_atom_type.clone(),
                    collection.d_params.clone(),
                )
            })
            .collect();

        let default = &collections[0];
        assert_eq!(default.d_i_atom_type.len(), 95);
        assert_eq!(default.d_params.len(), 95);
        for (index, key) in (1_u8..=82).enumerate() {
            assert_eq!(default.d_i_atom_type[index], key);
        }
        for (index, key) in (87_u8..=99).enumerate() {
            assert_eq!(default.d_i_atom_type[82 + index], key);
        }
        assert_eq!(default.d_params[0], prop([6, 4, 4, 0, 0, 0, 0, 0]));
        assert_eq!(default.d_params[31], prop([8, 1, 12, 1, 1, 0, 0, 0]));
        assert_eq!(default.d_params[54], prop([7, 3, 34, 0, 1, 0, 0, 0]));
        assert_eq!(default.d_params[82], prop([26, 0, 0, 0, 0, 0, 0, 0]));
        assert_eq!(default.d_params[94], prop([12, 0, 0, 0, 0, 0, 0, 0]));
        let mut flattened_default = Vec::with_capacity(855);
        for (&key, row) in default.d_i_atom_type.iter().zip(&default.d_params) {
            flattened_default.push(key);
            flattened_default.extend_from_slice(&[
                row.atno, row.crd, row.val, row.pilp, row.mltb, row.arom, row.linh, row.sbmb,
            ]);
        }
        assert_eq!(flattened_default.len(), 855);
        assert_eq!(fnv1a64(&flattened_default), 0xee7f_6305_eede_e4a8);

        assert!(collections[1].d_i_atom_type.is_empty());
        assert!(collections[1].d_params.is_empty());
        assert_eq!(collections[2].d_i_atom_type.as_slice(), &[0, 1, 1, 255]);
        assert_eq!(collections[2].d_params.as_slice(), &custom_rows);
        assert_eq!(collections[3].d_i_atom_type.as_slice(), &[2, 8]);
        assert_eq!(collections[3].d_params.as_slice(), &gap_rows);

        let mut actual_calls = 0;
        for (collection_index, collection) in collections.iter().enumerate() {
            for query in (0_u32..=256).chain(std::iter::once(u32::MAX)) {
                assert!(actual_calls < 1032, "unexpected extra Prop lookup");
                assert_eq!(collection.d_i_atom_type.len(), collection.d_params.len());
                assert_eq!(
                    collection.d_i_atom_type.as_slice(),
                    snapshots[collection_index].0.as_slice()
                );
                assert_eq!(
                    collection.d_params.as_slice(),
                    snapshots[collection_index].1.as_slice()
                );
                assert!(
                    collection
                        .d_i_atom_type
                        .windows(2)
                        .all(|pair| pair[0] <= pair[1])
                );

                let expected_index = match collection_index {
                    0 => match query {
                        1..=82 => Some((query - 1) as usize),
                        87..=99 => Some((query - 5) as usize),
                        _ => None,
                    },
                    1 => None,
                    2 => match query {
                        0 => Some(0),
                        1 => Some(1),
                        255 => Some(3),
                        _ => None,
                    },
                    3 => match query {
                        2 => Some(0),
                        8 => Some(1),
                        _ => None,
                    },
                    _ => panic!("unexpected Prop collection {collection_index}"),
                };

                let actual = collection.get(query);
                actual_calls += 1;
                match (actual, expected_index) {
                    (Some(found), Some(index)) => {
                        assert_eq!(u32::from(collection.d_i_atom_type[index]), query);
                        assert!(std::ptr::eq(found, &collection.d_params[index]));
                        match collection_index {
                            2 => assert_eq!(*found, custom_rows[index]),
                            3 => assert_eq!(*found, gap_rows[index]),
                            _ => {}
                        }
                    }
                    (None, None) => {}
                    (actual, expected) => panic!(
                        "Prop query {query} in collection {collection_index}: {actual:?}, expected index {expected:?}"
                    ),
                }

                if collection_index == 0 {
                    match query {
                        83..=86 | 256 | u32::MAX => assert!(actual.is_none()),
                        87 => assert_eq!(actual.copied(), Some(prop([26, 0, 0, 0, 0, 0, 0, 0]))),
                        99 => assert_eq!(actual.copied(), Some(prop([12, 0, 0, 0, 0, 0, 0, 0]))),
                        _ => {}
                    }
                }
                if collection_index == 2 && query == 1 {
                    assert_eq!(actual.copied(), Some(custom_rows[1]));
                    assert_ne!(actual.copied(), Some(custom_rows[2]));
                }

                assert_eq!(
                    collection.d_i_atom_type.as_slice(),
                    snapshots[collection_index].0.as_slice()
                );
                assert_eq!(
                    collection.d_params.as_slice(),
                    snapshots[collection_index].1.as_slice()
                );
            }
        }
        assert_eq!(actual_calls, 1032);
    }

    #[test]
    fn mmff_dp_m09_default_accessors_are_stable_across_threads() {
        fn def_field_hash(collection: &MmffDefCollection) -> u64 {
            let mut fields = Vec::with_capacity(380);
            for row in &collection.d_params {
                fields.extend_from_slice(&row.eq_level);
            }
            assert_eq!(fields.len(), 380);
            fnv1a64(&fields)
        }

        fn prop_field_hash(collection: &MmffPropCollection) -> u64 {
            assert_eq!(collection.d_i_atom_type.len(), collection.d_params.len());
            let mut fields = Vec::with_capacity(855);
            for (&key, row) in collection.d_i_atom_type.iter().zip(&collection.d_params) {
                fields.extend_from_slice(&[
                    key, row.atno, row.crd, row.val, row.pilp, row.mltb, row.arom, row.linh,
                    row.sbmb,
                ]);
            }
            assert_eq!(fields.len(), 855);
            fnv1a64(&fields)
        }

        let mut acquisition_calls = 0;
        assert!(acquisition_calls < 1, "unexpected extra Def warm-up call");
        let default_def = default_mmff_def().expect("frozen Def asset must parse");
        acquisition_calls += 1;
        assert_eq!(
            DEFAULT_MMFF_DEF_CONSTRUCTIONS.load(Ordering::Relaxed),
            1,
            "the Def accessor retains one actual constructor result"
        );

        assert!(acquisition_calls < 2, "unexpected extra Prop warm-up call");
        let default_prop = default_mmff_prop().expect("frozen Prop asset must parse");
        acquisition_calls += 1;
        assert_eq!(acquisition_calls, 2);
        let initialized_constructor_counts = (
            DEFAULT_MMFF_DEF_CONSTRUCTIONS.load(Ordering::Relaxed),
            DEFAULT_MMFF_PROP_CONSTRUCTIONS.load(Ordering::Relaxed),
        );
        assert_eq!(initialized_constructor_counts, (1, 1));

        assert_eq!(default_def.d_params.len(), 95);
        assert_eq!(default_def.d_params[0].eq_level, [1, 1, 1, 0]);
        assert_eq!(default_def.d_params[19].eq_level, [20, 1, 1, 0]);
        assert_eq!(default_def.d_params[69].eq_level, [70, 70, 70, 70]);
        assert_eq!(default_def.d_params[81].eq_level, [82, 9, 8, 0]);
        assert_eq!(default_def.d_params[82].eq_level, [87; 4]);
        assert_eq!(default_def.d_params[94].eq_level, [99; 4]);
        assert_eq!(def_field_hash(default_def), 0x3e49_ed6b_8e84_98d5);

        assert_eq!(default_prop.d_i_atom_type.len(), 95);
        assert_eq!(default_prop.d_params.len(), 95);
        assert_eq!(default_prop.d_i_atom_type[0], 1);
        assert_eq!(default_prop.d_params[0], prop([6, 4, 4, 0, 0, 0, 0, 0]));
        assert_eq!(default_prop.d_i_atom_type[31], 32);
        assert_eq!(default_prop.d_params[31], prop([8, 1, 12, 1, 1, 0, 0, 0]));
        assert_eq!(default_prop.d_i_atom_type[54], 55);
        assert_eq!(default_prop.d_params[54], prop([7, 3, 34, 0, 1, 0, 0, 0]));
        assert_eq!(default_prop.d_i_atom_type[82], 87);
        assert_eq!(default_prop.d_params[82], prop([26, 0, 0, 0, 0, 0, 0, 0]));
        assert_eq!(default_prop.d_i_atom_type[94], 99);
        assert_eq!(default_prop.d_params[94], prop([12, 0, 0, 0, 0, 0, 0, 0]));
        assert_eq!(prop_field_hash(default_prop), 0xee7f_6305_eede_e4a8);

        let def_snapshot = default_def.d_params.clone();
        let prop_key_snapshot = default_prop.d_i_atom_type.clone();
        let prop_row_snapshot = default_prop.d_params.clone();

        let per_thread_accessor_calls = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(4);
            for _worker_index in 0..4 {
                let def_anchor = default_def;
                let prop_anchor = default_prop;
                let def_rows = &def_snapshot;
                let prop_keys = &prop_key_snapshot;
                let prop_rows = &prop_row_snapshot;
                handles.push(scope.spawn(move || {
                    let mut accessor_calls = 0;
                    for _repeat in 0..3 {
                        assert!(accessor_calls < 6, "unexpected extra accessor request");
                        assert_eq!(
                            (
                                DEFAULT_MMFF_DEF_CONSTRUCTIONS.load(Ordering::Relaxed),
                                DEFAULT_MMFF_PROP_CONSTRUCTIONS.load(Ordering::Relaxed),
                            ),
                            initialized_constructor_counts
                        );
                        let actual_def = default_mmff_def().expect("cached Def result");
                        accessor_calls += 1;
                        assert!(std::ptr::eq(actual_def, def_anchor));
                        assert_eq!(actual_def.d_params.as_slice(), def_rows.as_slice());
                        assert_eq!(def_field_hash(actual_def), 0x3e49_ed6b_8e84_98d5);
                        assert_eq!(
                            DEFAULT_MMFF_DEF_CONSTRUCTIONS.load(Ordering::Relaxed),
                            initialized_constructor_counts.0
                        );

                        assert!(accessor_calls < 6, "unexpected extra accessor request");
                        let actual_prop = default_mmff_prop().expect("cached Prop result");
                        accessor_calls += 1;
                        assert!(std::ptr::eq(actual_prop, prop_anchor));
                        assert_eq!(actual_prop.d_i_atom_type.as_slice(), prop_keys.as_slice());
                        assert_eq!(actual_prop.d_params.as_slice(), prop_rows.as_slice());
                        assert_eq!(prop_field_hash(actual_prop), 0xee7f_6305_eede_e4a8);
                        assert_eq!(
                            DEFAULT_MMFF_PROP_CONSTRUCTIONS.load(Ordering::Relaxed),
                            initialized_constructor_counts.1
                        );
                    }
                    accessor_calls
                }));
            }
            handles
                .into_iter()
                .map(|handle| handle.join().expect("scoped accessor worker panicked"))
                .collect::<Vec<_>>()
        });
        assert_eq!(per_thread_accessor_calls, [6; 4]);
        assert_eq!(per_thread_accessor_calls.iter().sum::<usize>(), 24);

        let mut custom_constructor_calls = 0;
        assert!(
            custom_constructor_calls < 1,
            "unexpected extra custom Def row"
        );
        let custom_def = MmffDefCollection::from_text("custom\t250\t2\t3\t4\t5\n")
            .expect("valid independent custom Def row");
        custom_constructor_calls += 1;
        assert_eq!(
            custom_def.d_params.as_slice(),
            &[MmffDef {
                eq_level: [2, 3, 4, 5]
            }]
        );
        assert!(!std::ptr::eq(&custom_def, default_def));
        assert_eq!(default_def.d_params.as_slice(), def_snapshot.as_slice());
        assert_eq!(
            default_prop.d_i_atom_type.as_slice(),
            prop_key_snapshot.as_slice()
        );
        assert_eq!(
            default_prop.d_params.as_slice(),
            prop_row_snapshot.as_slice()
        );
        assert_eq!(
            (
                DEFAULT_MMFF_DEF_CONSTRUCTIONS.load(Ordering::Relaxed),
                DEFAULT_MMFF_PROP_CONSTRUCTIONS.load(Ordering::Relaxed),
            ),
            initialized_constructor_counts
        );

        assert!(
            custom_constructor_calls < 2,
            "unexpected extra custom Prop row"
        );
        let custom_prop = MmffPropCollection::from_text("250\t6\t4\t3\t2\t1\t0\t1\t2\n")
            .expect("valid independent custom Prop row");
        custom_constructor_calls += 1;
        assert_eq!(custom_prop.d_i_atom_type, [250]);
        assert_eq!(custom_prop.d_params, [prop([6, 4, 3, 2, 1, 0, 1, 2])]);
        assert!(!std::ptr::eq(&custom_prop, default_prop));
        assert_eq!(default_def.d_params.as_slice(), def_snapshot.as_slice());
        assert_eq!(
            default_prop.d_i_atom_type.as_slice(),
            prop_key_snapshot.as_slice()
        );
        assert_eq!(
            default_prop.d_params.as_slice(),
            prop_row_snapshot.as_slice()
        );
        assert_eq!(custom_constructor_calls, 2);
        assert_eq!(
            (
                DEFAULT_MMFF_DEF_CONSTRUCTIONS.load(Ordering::Relaxed),
                DEFAULT_MMFF_PROP_CONSTRUCTIONS.load(Ordering::Relaxed),
            ),
            initialized_constructor_counts
        );

        let mut post_custom_accessor_calls = 0;
        assert!(post_custom_accessor_calls < 1);
        let def_after_custom = default_mmff_def().expect("retained Def after custom construction");
        post_custom_accessor_calls += 1;
        assert!(std::ptr::eq(def_after_custom, default_def));
        assert_eq!(
            def_after_custom.d_params.as_slice(),
            def_snapshot.as_slice()
        );

        assert!(post_custom_accessor_calls < 2);
        let prop_after_custom =
            default_mmff_prop().expect("retained Prop after custom construction");
        post_custom_accessor_calls += 1;
        assert!(std::ptr::eq(prop_after_custom, default_prop));
        assert_eq!(
            prop_after_custom.d_i_atom_type.as_slice(),
            prop_key_snapshot.as_slice()
        );
        assert_eq!(
            prop_after_custom.d_params.as_slice(),
            prop_row_snapshot.as_slice()
        );
        assert_eq!(post_custom_accessor_calls, 2);
        assert_eq!(
            (
                DEFAULT_MMFF_DEF_CONSTRUCTIONS.load(Ordering::Relaxed),
                DEFAULT_MMFF_PROP_CONSTRUCTIONS.load(Ordering::Relaxed),
            ),
            initialized_constructor_counts
        );
    }

    #[test]
    fn mmff_pc_pbci_constructor_formats_preserve_rows_and_lookup_products() {
        const LF_SOURCE: &str = "discard\tNOT_NUMERIC\t1.25\t-0\textra\nnotnumber\talsoignored\t-2.5\t0.5\nsame\tsame\t3\t4\n";
        const CRLF_SOURCE: &str = "discard\tNOT_NUMERIC\t1.25\t-0\textra\r\nnotnumber\talsoignored\t-2.5\t0.5\r\nsame\tsame\t3\t4\r\n";
        const COMMENT_SOURCE: &str = "*comment\ndiscard\tNOT_NUMERIC\t1.25\t-0\textra\nnotnumber\talsoignored\t-2.5\t0.5\nsame\tsame\t3\t4\n";
        const UNTERMINATED_SOURCE: &str = "discard\tNOT_NUMERIC\t1.25\t-0\textra\nnotnumber\talsoignored\t-2.5\t0.5\nsame\tsame\t3\t4";

        const FULL_ROWS: [[u64; 2]; 3] = [
            [0x3ff4_0000_0000_0000, 0x8000_0000_0000_0000],
            [0xc004_0000_0000_0000, 0x3fe0_0000_0000_0000],
            [0x4008_0000_0000_0000, 0x4010_0000_0000_0000],
        ];
        const TWO_ROWS: [[u64; 2]; 2] = [FULL_ROWS[0], FULL_ROWS[1]];

        const FULL_LOOKUPS: [(u32, Option<(u64, u64)>); 6] = [
            (0, None),
            (1, Some((0x3ff4_0000_0000_0000, 0x8000_0000_0000_0000))),
            (2, Some((0xc004_0000_0000_0000, 0x3fe0_0000_0000_0000))),
            (3, Some((0x4008_0000_0000_0000, 0x4010_0000_0000_0000))),
            (4, None),
            (u32::MAX, None),
        ];
        const TWO_ROW_LOOKUPS: [(u32, Option<(u64, u64)>); 6] = [
            (0, None),
            (1, Some((0x3ff4_0000_0000_0000, 0x8000_0000_0000_0000))),
            (2, Some((0xc004_0000_0000_0000, 0x3fe0_0000_0000_0000))),
            (3, None),
            (4, None),
            (u32::MAX, None),
        ];

        let formats: [(&str, &[[u64; 2]], &[(u32, Option<(u64, u64)>)]); 4] = [
            (LF_SOURCE, &FULL_ROWS, &FULL_LOOKUPS),
            (CRLF_SOURCE, &FULL_ROWS, &FULL_LOOKUPS),
            (COMMENT_SOURCE, &FULL_ROWS, &FULL_LOOKUPS),
            (UNTERMINATED_SOURCE, &TWO_ROWS, &TWO_ROW_LOOKUPS),
        ];

        let mut actual_formats = 0;
        let mut actual_get_calls = 0;
        for (source, expected_rows, expected_queries) in formats {
            let input = source.to_owned();
            let input_before = input.as_bytes().to_vec();
            let input_address = input.as_ptr();
            let collection = MmffPbciCollection::from_text(&input)
                .expect("frozen PBCI constructor fixture must parse");

            assert_eq!(input.as_ptr(), input_address);
            assert_eq!(input.as_bytes(), input_before);
            assert_eq!(collection.d_params.len(), expected_rows.len());
            for (row_index, (actual, expected_bits)) in
                collection.d_params.iter().zip(expected_rows).enumerate()
            {
                assert_eq!(
                    (actual.pbci.to_bits(), actual.fcadj.to_bits()),
                    (expected_bits[0], expected_bits[1]),
                    "format {actual_formats}, row {row_index}"
                );
            }

            assert_eq!(expected_queries.len(), 6);
            for &(atom_type, expected_bits) in expected_queries {
                assert!(actual_get_calls < 24, "unexpected extra PBCI get call");
                let actual = collection
                    .get(atom_type)
                    .map(|row| (row.pbci.to_bits(), row.fcadj.to_bits()));
                assert_eq!(
                    actual, expected_bits,
                    "format {actual_formats}, atom type {atom_type}"
                );
                actual_get_calls += 1;
            }
            actual_formats += 1;
        }

        assert_eq!(actual_formats, 4);
        assert_eq!(actual_get_calls, 24);
    }

    #[test]
    fn mmff_pc_pbci_constructor_reports_fixed_typed_errors() {
        let cases: [(&str, usize, usize, MmffParamParseCause); 7] = [
            (
                "*header\nx\ty\tbroken\t0\n",
                2,
                2,
                MmffParamParseCause::InvalidFloat {
                    cell: "broken".to_owned(),
                },
            ),
            (
                "x\ty\t1\tbroken\n",
                1,
                3,
                MmffParamParseCause::InvalidFloat {
                    cell: "broken".to_owned(),
                },
            ),
            ("\t\n", 1, 0, MmffParamParseCause::MissingToken),
            ("x\t\n", 1, 1, MmffParamParseCause::MissingToken),
            ("x\ty\t\n", 1, 2, MmffParamParseCause::MissingToken),
            ("x\ty\t1\t\n", 1, 3, MmffParamParseCause::MissingToken),
            ("*header\n\n", 2, 0, MmffParamParseCause::EmptyProcessedLine),
        ];

        let mut actual_error_calls = 0;
        for (source, line, column, cause) in cases {
            assert!(actual_error_calls < 7, "unexpected extra PBCI error case");
            let actual = MmffPbciCollection::from_text(source)
                .expect_err("frozen malformed PBCI record must return its typed error");
            assert_eq!(
                actual,
                MmffParamParseError {
                    table: MmffParamTable::Pbci,
                    line,
                    column,
                    cause,
                }
            );
            actual_error_calls += 1;
        }
        assert_eq!(actual_error_calls, 7);
    }

    #[test]
    fn mmff_pc_pbci_default_asset_sentinels_and_shared_borrows() {
        const SENTINEL_TYPES: [u32; 4] = [1, 83, 87, 99];
        const EXPECTED_BITS: [(u64, u64); 4] = [
            (0x0000_0000_0000_0000, 0x0000_0000_0000_0000),
            (0x0000_0000_0000_0000, 0x0000_0000_0000_0000),
            (0x4000_0000_0000_0000, 0x0000_0000_0000_0000),
            (0x4000_0000_0000_0000, 0x0000_0000_0000_0000),
        ];

        fn sentinel_bits(collection: &MmffPbciCollection) -> [(u64, u64); 4] {
            SENTINEL_TYPES.map(|atom_type| {
                let row = collection
                    .get(atom_type)
                    .expect("frozen default PBCI sentinel must be present");
                (row.pbci.to_bits(), row.fcadj.to_bits())
            })
        }

        assert_eq!(DEFAULT_MMFF_PBCI_CONSTRUCTIONS.load(Ordering::Relaxed), 0);
        let collection = default_mmff_pbci().expect("frozen default PBCI asset must parse");
        assert_eq!(DEFAULT_MMFF_PBCI_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(collection.d_params.len(), 99);
        assert_eq!(sentinel_bits(collection), EXPECTED_BITS);

        for atom_type in [0, 100, u32::MAX] {
            assert!(collection.get(atom_type).is_none());
        }
        assert_eq!(collection.get(1).map(|row| row.pbci.to_bits()), Some(0));
        assert_eq!(
            collection.get(99).map(|row| row.pbci.to_bits()),
            Some(0x4000_0000_0000_0000)
        );

        let collection_address = collection as *const MmffPbciCollection as usize;
        let workers = (0..8)
            .map(|_| {
                std::thread::spawn(|| {
                    (0..32)
                        .map(|_| {
                            let collection = default_mmff_pbci()
                                .expect("cached default PBCI collection must remain available");
                            (
                                collection as *const MmffPbciCollection as usize,
                                sentinel_bits(collection),
                            )
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect::<Vec<_>>();

        let mut actual_shared_accesses = 0;
        for worker in workers {
            for (actual_address, actual_bits) in worker
                .join()
                .expect("default PBCI access worker must complete")
            {
                assert!(actual_shared_accesses < 8 * 32);
                assert_eq!(actual_address, collection_address);
                assert_eq!(actual_bits, EXPECTED_BITS);
                actual_shared_accesses += 1;
            }
        }
        assert_eq!(actual_shared_accesses, 8 * 32);
        assert_eq!(DEFAULT_MMFF_PBCI_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn mmff_pc_chg_default_asset_sentinels_and_shared_borrows() {
        const FIRST_BCI_BITS: u64 = 0x0000_0000_0000_0000;
        const LAST_BCI_BITS: u64 = 0xbfd9_9999_9999_999a;

        fn sentinel_observation(
            collection: &MmffChgCollection,
        ) -> (i32, u64, i32, u64, i32, usize, usize, usize) {
            let (first_sign, first) = collection.get(0, 1, 1);
            let first = first.expect("source-first default Chg row must be present");
            let (last_sign, last) = collection.get(0, 80, 81);
            let last = last.expect("source-last default Chg row must be present");
            let (reverse_sign, reverse) = collection.get(0, 81, 80);
            let reverse = reverse.expect("reverse source-last default Chg row must be present");
            assert!(std::ptr::eq(last, reverse));

            (
                first_sign,
                first.bci.to_bits(),
                last_sign,
                last.bci.to_bits(),
                reverse_sign,
                first as *const _ as usize,
                last as *const _ as usize,
                reverse as *const _ as usize,
            )
        }

        assert_eq!(DEFAULT_MMFF_CHG_CONSTRUCTIONS.load(Ordering::Relaxed), 0);
        let collection = default_mmff_chg().expect("frozen default Chg asset must parse");
        assert_eq!(DEFAULT_MMFF_CHG_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(collection.d_params.len(), 498);
        assert_eq!(collection.d_i_atom_type.len(), 498);
        assert_eq!(collection.d_j_atom_type.len(), 498);
        assert_eq!(collection.d_bond_type.len(), 498);

        let expected_observation = (-1, FIRST_BCI_BITS, -1, LAST_BCI_BITS, 1);
        let initial = sentinel_observation(collection);
        assert_eq!(initial.0, expected_observation.0);
        assert_eq!(initial.1, expected_observation.1);
        assert_eq!(initial.2, expected_observation.2);
        assert_eq!(initial.3, expected_observation.3);
        assert_eq!(initial.4, expected_observation.4);
        assert_eq!(initial.6, initial.7);

        let collection_address = collection as *const MmffChgCollection as usize;
        let bond_snapshot = collection.d_bond_type.clone();
        let i_snapshot = collection.d_i_atom_type.clone();
        let j_snapshot = collection.d_j_atom_type.clone();
        let bci_snapshot = collection
            .d_params
            .iter()
            .map(|row| row.bci.to_bits())
            .collect::<Vec<_>>();
        let workers = (0..8)
            .map(|_| {
                std::thread::spawn(|| {
                    (0..32)
                        .map(|_| {
                            let collection = default_mmff_chg()
                                .expect("cached default Chg collection must remain available");
                            (
                                collection as *const MmffChgCollection as usize,
                                sentinel_observation(collection),
                            )
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect::<Vec<_>>();

        let mut actual_shared_accesses = 0;
        for worker in workers {
            for (actual_address, actual_observation) in worker
                .join()
                .expect("default Chg access worker must complete")
            {
                assert!(actual_shared_accesses < 8 * 32);
                assert_eq!(actual_address, collection_address);
                assert_eq!(actual_observation.0, expected_observation.0);
                assert_eq!(actual_observation.1, expected_observation.1);
                assert_eq!(actual_observation.2, expected_observation.2);
                assert_eq!(actual_observation.3, expected_observation.3);
                assert_eq!(actual_observation.4, expected_observation.4);
                assert_eq!(actual_observation.5, initial.5);
                assert_eq!(actual_observation.6, initial.6);
                assert_eq!(actual_observation.7, initial.7);
                actual_shared_accesses += 1;
            }
        }
        assert_eq!(actual_shared_accesses, 8 * 32);
        assert_eq!(collection.d_bond_type.as_slice(), bond_snapshot.as_slice());
        assert_eq!(collection.d_i_atom_type.as_slice(), i_snapshot.as_slice());
        assert_eq!(collection.d_j_atom_type.as_slice(), j_snapshot.as_slice());
        assert_eq!(
            collection
                .d_params
                .iter()
                .map(|row| row.bci.to_bits())
                .collect::<Vec<_>>(),
            bci_snapshot
        );
        assert_eq!(DEFAULT_MMFF_CHG_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn mmff_pc_chg_ctor_formats_casts_and_typed_errors() {
        const LF_SOURCE: &str = "0\t1\t1\t0.125\tignored\n0\t1\t2\t0.25\tignored\n0\t1\t2\t0.5\tignored\n1\t1\t2\t-0.75\tignored\n0\t2\t3\t1.25\tignored\n1\t2\t3\t-1.5\tignored\n";
        const CRLF_SOURCE: &str = "0\t1\t1\t0.125\tignored\r\n0\t1\t2\t0.25\tignored\r\n0\t1\t2\t0.5\tignored\r\n1\t1\t2\t-0.75\tignored\r\n0\t2\t3\t1.25\tignored\r\n1\t2\t3\t-1.5\tignored\r\n";
        const COMMENT_SOURCE: &str = "*comment\n0\t1\t1\t0.125\tignored\n0\t1\t2\t0.25\tignored\n0\t1\t2\t0.5\tignored\n1\t1\t2\t-0.75\tignored\n0\t2\t3\t1.25\tignored\n1\t2\t3\t-1.5\tignored\n";
        const UNTERMINATED_SOURCE: &str = "0\t1\t1\t0.125\tignored\n0\t1\t2\t0.25\tignored\n0\t1\t2\t0.5\tignored\n1\t1\t2\t-0.75\tignored\n0\t2\t3\t1.25\tignored\n1\t2\t3\t-1.5\tignored";
        const FULL_BOND_TYPES: [u8; 6] = [0, 0, 0, 1, 0, 1];
        const FULL_I_ATOM_TYPES: [u8; 6] = [1, 1, 1, 1, 2, 2];
        const FULL_J_ATOM_TYPES: [u8; 6] = [1, 2, 2, 2, 3, 3];
        const FULL_BCI_BITS: [u64; 6] = [
            0x3fc0_0000_0000_0000,
            0x3fd0_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0xbfe8_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0xbff8_0000_0000_0000,
        ];
        const FORMATS: [(&str, usize); 4] = [
            (LF_SOURCE, 6),
            (CRLF_SOURCE, 6),
            (COMMENT_SOURCE, 6),
            (UNTERMINATED_SOURCE, 5),
        ];

        let mut actual_formats = 0;
        for (source, expected_len) in FORMATS {
            let input = source.to_owned();
            let input_before = input.as_bytes().to_vec();
            let input_address = input.as_ptr();
            let collection = MmffChgCollection::from_text(&input)
                .expect("frozen Chg constructor fixture must parse");

            assert_eq!(input.as_ptr(), input_address);
            assert_eq!(input.as_bytes(), input_before);
            assert_eq!(collection.d_params.len(), expected_len);
            assert_eq!(collection.d_bond_type.len(), expected_len);
            assert_eq!(collection.d_i_atom_type.len(), expected_len);
            assert_eq!(collection.d_j_atom_type.len(), expected_len);
            assert_eq!(
                collection.d_bond_type.as_slice(),
                &FULL_BOND_TYPES[..expected_len]
            );
            assert_eq!(
                collection.d_i_atom_type.as_slice(),
                &FULL_I_ATOM_TYPES[..expected_len]
            );
            assert_eq!(
                collection.d_j_atom_type.as_slice(),
                &FULL_J_ATOM_TYPES[..expected_len]
            );
            for (row_index, (actual, expected_bits)) in collection
                .d_params
                .iter()
                .zip(&FULL_BCI_BITS[..expected_len])
                .enumerate()
            {
                assert_eq!(
                    actual.bci.to_bits(),
                    *expected_bits,
                    "format {actual_formats}, row {row_index}"
                );
            }
            actual_formats += 1;
        }
        assert_eq!(actual_formats, 4);

        const CAST_SOURCE: &str = "256\t-1\t256\t1.25\n-1\t256\t-1\t-0\n";
        let cast_input = CAST_SOURCE.to_owned();
        let cast_input_before = cast_input.as_bytes().to_vec();
        let cast_collection = MmffChgCollection::from_text(&cast_input)
            .expect("frozen narrowing controls must parse");
        assert_eq!(cast_input.as_bytes(), cast_input_before);
        assert_eq!(cast_collection.d_bond_type.as_slice(), [0, 255]);
        assert_eq!(cast_collection.d_i_atom_type.as_slice(), [255, 0]);
        assert_eq!(cast_collection.d_j_atom_type.as_slice(), [0, 255]);
        assert_eq!(
            cast_collection
                .d_params
                .iter()
                .map(|row| row.bci.to_bits())
                .collect::<Vec<_>>(),
            [0x3ff4_0000_0000_0000, 0x8000_0000_0000_0000]
        );

        let cases: [(&str, usize, usize, MmffParamParseCause); 9] = [
            (
                "bad\t1\t2\t0.125\n",
                1,
                0,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "bad".to_owned(),
                },
            ),
            (
                "0\tbad\t2\t0.125\n",
                1,
                1,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "bad".to_owned(),
                },
            ),
            (
                "0\t1\tbad\t0.125\n",
                1,
                2,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "bad".to_owned(),
                },
            ),
            (
                "0\t1\t2\tbad\n",
                1,
                3,
                MmffParamParseCause::InvalidFloat {
                    cell: "bad".to_owned(),
                },
            ),
            ("\t\n", 1, 0, MmffParamParseCause::MissingToken),
            ("0\t\n", 1, 1, MmffParamParseCause::MissingToken),
            ("0\t1\t\n", 1, 2, MmffParamParseCause::MissingToken),
            ("0\t1\t2\t\n", 1, 3, MmffParamParseCause::MissingToken),
            ("*header\n\n", 2, 0, MmffParamParseCause::EmptyProcessedLine),
        ];

        let mut actual_error_calls = 0;
        for (source, line, column, cause) in cases {
            assert!(actual_error_calls < 9, "unexpected extra Chg error case");
            let actual = MmffChgCollection::from_text(source)
                .expect_err("frozen malformed Chg record must return its typed error");
            assert_eq!(
                actual,
                MmffParamParseError {
                    table: MmffParamTable::Chg,
                    line,
                    column,
                    cause,
                }
            );
            actual_error_calls += 1;
        }
        assert_eq!(actual_error_calls, 9);
    }

    #[test]
    fn mmff_pc_chg_lookup_source_matrix_and_full_width_boundaries() {
        const SOURCE: &str = "0\t1\t1\t0.125\n0\t1\t2\t0.25\n0\t1\t2\t0.5\n1\t1\t2\t-0.75\n0\t2\t3\t1.25\n1\t2\t3\t-1.5\n";
        const MATRIX: [(u32, u32, u32, i32, Option<(usize, u64)>); 18] = [
            (0, 1, 1, -1, Some((0, 0x3fc0_0000_0000_0000))),
            (1, 1, 1, -1, None),
            (2, 1, 1, -1, None),
            (0, 1, 2, -1, Some((1, 0x3fd0_0000_0000_0000))),
            (1, 1, 2, -1, Some((3, 0xbfe8_0000_0000_0000))),
            (2, 1, 2, -1, None),
            (0, 2, 1, 1, Some((1, 0x3fd0_0000_0000_0000))),
            (1, 2, 1, 1, Some((3, 0xbfe8_0000_0000_0000))),
            (2, 2, 1, 1, None),
            (0, 2, 3, -1, Some((4, 0x3ff4_0000_0000_0000))),
            (1, 2, 3, -1, Some((5, 0xbff8_0000_0000_0000))),
            (2, 2, 3, -1, None),
            (0, 3, 2, 1, Some((4, 0x3ff4_0000_0000_0000))),
            (1, 3, 2, 1, Some((5, 0xbff8_0000_0000_0000))),
            (2, 3, 2, 1, None),
            (0, 9, 9, -1, None),
            (1, 9, 9, -1, None),
            (2, 9, 9, -1, None),
        ];
        const BOUNDARIES: [(u32, u32, u32, i32); 7] = [
            (256, 1, 2, -1),
            (u32::MAX, 1, 2, -1),
            (257, 1, 2, -1),
            (0, 257, 258, -1),
            (0, 256, 2, 1),
            (0, u32::MAX, 2, 1),
            (0, 2, u32::MAX, -1),
        ];

        let input = SOURCE.to_owned();
        let input_before = input.as_bytes().to_vec();
        let input_address = input.as_ptr();
        let collection =
            MmffChgCollection::from_text(&input).expect("sorted source Chg fixture must parse");
        assert_eq!(input.as_ptr(), input_address);
        assert_eq!(input.as_bytes(), input_before);
        assert_eq!(collection.d_params.len(), 6);
        assert_eq!(collection.d_i_atom_type.as_slice(), [1, 1, 1, 1, 2, 2]);
        assert_eq!(collection.d_j_atom_type.as_slice(), [1, 2, 2, 2, 3, 3]);
        assert_eq!(collection.d_bond_type.as_slice(), [0, 0, 0, 1, 0, 1]);

        let input_snapshot = input.as_bytes().to_vec();
        let bond_snapshot = collection.d_bond_type.clone();
        let i_snapshot = collection.d_i_atom_type.clone();
        let j_snapshot = collection.d_j_atom_type.clone();
        let bci_snapshot = collection
            .d_params
            .iter()
            .map(|row| row.bci.to_bits())
            .collect::<Vec<_>>();

        let mut actual_matrix_calls = 0;
        for (bond_type, i_atom_type, j_atom_type, expected_sign, expected_row) in MATRIX {
            assert!(actual_matrix_calls < 18, "unexpected extra Chg lookup call");
            let (actual_sign, actual_row) = collection.get(bond_type, i_atom_type, j_atom_type);
            assert_eq!(actual_sign, expected_sign);
            match expected_row {
                Some((expected_index, expected_bits)) => {
                    let actual_row = actual_row.expect("frozen Chg lookup row must exist");
                    assert_eq!(actual_row.bci.to_bits(), expected_bits);
                    assert!(std::ptr::eq(
                        actual_row,
                        &collection.d_params[expected_index]
                    ));
                }
                None => assert!(actual_row.is_none()),
            }
            assert_eq!(input.as_bytes(), input_snapshot.as_slice());
            assert_eq!(collection.d_bond_type.as_slice(), bond_snapshot.as_slice());
            assert_eq!(collection.d_i_atom_type.as_slice(), i_snapshot.as_slice());
            assert_eq!(collection.d_j_atom_type.as_slice(), j_snapshot.as_slice());
            assert_eq!(
                collection
                    .d_params
                    .iter()
                    .map(|row| row.bci.to_bits())
                    .collect::<Vec<_>>(),
                bci_snapshot
            );
            actual_matrix_calls += 1;
        }
        assert_eq!(actual_matrix_calls, 18);

        let mut actual_boundary_calls = 0;
        for (bond_type, i_atom_type, j_atom_type, expected_sign) in BOUNDARIES {
            assert!(
                actual_boundary_calls < 7,
                "unexpected extra Chg boundary call"
            );
            let (actual_sign, actual_row) = collection.get(bond_type, i_atom_type, j_atom_type);
            assert_eq!(actual_sign, expected_sign);
            assert!(actual_row.is_none());
            assert_eq!(input.as_bytes(), input_snapshot.as_slice());
            assert_eq!(collection.d_bond_type.as_slice(), bond_snapshot.as_slice());
            assert_eq!(collection.d_i_atom_type.as_slice(), i_snapshot.as_slice());
            assert_eq!(collection.d_j_atom_type.as_slice(), j_snapshot.as_slice());
            assert_eq!(
                collection
                    .d_params
                    .iter()
                    .map(|row| row.bci.to_bits())
                    .collect::<Vec<_>>(),
                bci_snapshot
            );
            actual_boundary_calls += 1;
        }
        assert_eq!(actual_boundary_calls, 7);
    }

    #[test]
    fn mmff_bond_collection_formats_lookup_narrowing_and_consumed_errors() {
        const SOURCE: &str = concat!(
            "0\t1\t1\t1.25\t1.5\n",
            "1\t1\t1\t2.5\t1.75\n",
            "0\t1\t2\t3.75\t2.0\n",
            "0\t1\t2\t4.0\t2.25\n",
            "1\t1\t2\t5.0\t2.5\n",
            "0\t2\t2\t6.25\t2.75\n",
        );
        const EXPECTED_KEYS: [(u8, u8, u8); 6] = [
            (0, 1, 1),
            (1, 1, 1),
            (0, 1, 2),
            (0, 1, 2),
            (1, 1, 2),
            (0, 2, 2),
        ];
        const EXPECTED_VALUE_BITS: [(u64, u64); 6] = [
            (0x3ff4_0000_0000_0000, 0x3ff8_0000_0000_0000),
            (0x4004_0000_0000_0000, 0x3ffc_0000_0000_0000),
            (0x400e_0000_0000_0000, 0x4000_0000_0000_0000),
            (0x4010_0000_0000_0000, 0x4002_0000_0000_0000),
            (0x4014_0000_0000_0000, 0x4004_0000_0000_0000),
            (0x4019_0000_0000_0000, 0x4006_0000_0000_0000),
        ];
        const QUERIES: [(u32, u32, u32, Option<usize>); 9] = [
            (0, 1, 1, Some(0)),
            (1, 1, 1, Some(1)),
            (0, 1, 2, Some(2)),
            (0, 2, 1, Some(2)),
            (1, 2, 1, Some(4)),
            (0, 2, 2, Some(5)),
            (2, 1, 2, None),
            (0, 257, 1, None),
            (256, 1, 1, None),
        ];

        let sources = [
            SOURCE.to_owned(),
            SOURCE.replace('\n', "\r\n"),
            SOURCE.replace('\t', "\t\t"),
            SOURCE.trim_end_matches('\n').to_owned(),
        ];
        let mut actual_get_calls = 0;
        for (format_index, source) in sources.iter().enumerate() {
            let source_before = source.as_bytes().to_vec();
            let source_address = source.as_ptr();
            let collection = MmffBondCollection::from_text(source)
                .expect("each frozen Bond table format must parse");
            assert_eq!(source.as_ptr(), source_address);
            assert_eq!(source.as_bytes(), source_before.as_slice());

            let expected_len = if format_index == 3 { 5 } else { 6 };
            assert_eq!(collection.d_params.len(), expected_len);
            assert_eq!(collection.d_bond_type.len(), expected_len);
            assert_eq!(collection.d_i_atom_type.len(), expected_len);
            assert_eq!(collection.d_j_atom_type.len(), expected_len);
            for row_index in 0..expected_len {
                assert_eq!(
                    (
                        collection.d_bond_type[row_index],
                        collection.d_i_atom_type[row_index],
                        collection.d_j_atom_type[row_index],
                    ),
                    EXPECTED_KEYS[row_index]
                );
                assert_eq!(
                    collection.d_params[row_index].kb.to_bits(),
                    EXPECTED_VALUE_BITS[row_index].0
                );
                assert_eq!(
                    collection.d_params[row_index].r0.to_bits(),
                    EXPECTED_VALUE_BITS[row_index].1
                );
            }

            for (query_index, &(bond_type, atom_type, nbr_atom_type, expected_index)) in
                QUERIES.iter().enumerate()
            {
                assert!(actual_get_calls < 36, "unexpected extra Bond lookup call");
                let before = collection.clone();
                let actual = collection.get(bond_type, atom_type, nbr_atom_type);
                assert_eq!(
                    collection, before,
                    "lookup mutated {format_index}/{query_index}"
                );
                assert_eq!(source.as_bytes(), source_before.as_slice());

                let expected_index = if format_index == 3 && expected_index == Some(5) {
                    None
                } else {
                    expected_index
                };
                match expected_index {
                    Some(row_index) => {
                        let actual = actual.expect("frozen Bond row must be present");
                        let expected = &collection.d_params[row_index];
                        assert!(std::ptr::eq(actual, expected));
                        assert_eq!(actual.kb.to_bits(), EXPECTED_VALUE_BITS[row_index].0);
                        assert_eq!(actual.r0.to_bits(), EXPECTED_VALUE_BITS[row_index].1);
                    }
                    None => assert!(actual.is_none()),
                }
                actual_get_calls += 1;
            }
        }
        assert_eq!(actual_get_calls, 36);

        let cast = MmffBondCollection::from_text("256\t257\t258\t1.25\t1.5\n")
            .expect("unsigned stored keys narrow to their source u8 carriers");
        assert_eq!(cast.d_bond_type, [0]);
        assert_eq!(cast.d_i_atom_type, [1]);
        assert_eq!(cast.d_j_atom_type, [2]);
        let stored_cast = cast.get(0, 1, 2).expect("narrowed stored keys must match");
        assert_eq!(stored_cast.kb.to_bits(), 0x3ff4_0000_0000_0000);
        assert_eq!(stored_cast.r0.to_bits(), 0x3ff8_0000_0000_0000);
        assert!(cast.get(256, 257, 258).is_none());

        assert_eq!(
            MmffBondCollection::from_text("\n").unwrap_err(),
            MmffParamParseError {
                table: MmffParamTable::Bond,
                line: 1,
                column: 0,
                cause: MmffParamParseCause::EmptyProcessedLine,
            }
        );
        assert_eq!(
            MmffBondCollection::from_text("0\t1\t1\t1.25\n").unwrap_err(),
            MmffParamParseError {
                table: MmffParamTable::Bond,
                line: 1,
                column: 4,
                cause: MmffParamParseCause::MissingToken,
            }
        );
        assert_eq!(
            MmffBondCollection::from_text("0\tX\t1\t1.25\t1.5\n").unwrap_err(),
            MmffParamParseError {
                table: MmffParamTable::Bond,
                line: 1,
                column: 1,
                cause: MmffParamParseCause::InvalidUnsigned { cell: "X".into() },
            }
        );
        assert_eq!(
            MmffBondCollection::from_text("0\t1\t1\tbad\t1.5\n").unwrap_err(),
            MmffParamParseError {
                table: MmffParamTable::Bond,
                line: 1,
                column: 3,
                cause: MmffParamParseCause::InvalidFloat { cell: "bad".into() },
            }
        );
        let trailing = MmffBondCollection::from_text("0\t1\t1\t1.25\t1.5\tbad\n")
            .expect("source constructor does not consume trailing columns");
        assert_eq!(trailing.d_params.len(), 1);
        assert_eq!(
            trailing.get(0, 1, 1).unwrap().kb.to_bits(),
            0x3ff4_0000_0000_0000
        );
    }

    #[test]
    fn mmff_bndk_collection_formats_signed_lookup_narrowing_and_errors() {
        const SOURCE: &str = concat!(
            "1\t1\t1.5\t1.25\n",
            "1\t2\t1.75\t2.5\n",
            "1\t2\t2.0\t3.75\n",
            "2\t2\t2.25\t4.0\n",
        );
        const EXPECTED_KEYS: [(u8, u8); 4] = [(1, 1), (1, 2), (1, 2), (2, 2)];
        const EXPECTED_VALUE_BITS: [(u64, u64); 4] = [
            (0x3ff8_0000_0000_0000, 0x3ff4_0000_0000_0000),
            (0x3ffc_0000_0000_0000, 0x4004_0000_0000_0000),
            (0x4000_0000_0000_0000, 0x400e_0000_0000_0000),
            (0x4002_0000_0000_0000, 0x4010_0000_0000_0000),
        ];
        const QUERIES: [(i32, i32, Option<usize>); 8] = [
            (1, 1, Some(0)),
            (1, 2, Some(1)),
            (2, 1, Some(1)),
            (2, 2, Some(3)),
            (0, 1, None),
            (257, 1, None),
            (-1, 1, None),
            (i32::MIN, 1, None),
        ];

        let sources = [
            SOURCE.to_owned(),
            SOURCE.replace('\n', "\r\n"),
            SOURCE.replace('\t', "\t\t"),
            SOURCE.trim_end_matches('\n').to_owned(),
        ];
        let mut actual_get_calls = 0;
        for (format_index, source) in sources.iter().enumerate() {
            let source_before = source.as_bytes().to_vec();
            let source_address = source.as_ptr();
            let collection = MmffBndkCollection::from_text(source)
                .expect("each frozen Bndk table format must parse");
            assert_eq!(source.as_ptr(), source_address);
            assert_eq!(source.as_bytes(), source_before.as_slice());

            let expected_len = if format_index == 3 { 3 } else { 4 };
            assert_eq!(collection.d_params.len(), expected_len);
            assert_eq!(collection.d_i_atomic_num.len(), expected_len);
            assert_eq!(collection.d_j_atomic_num.len(), expected_len);
            for row_index in 0..expected_len {
                assert_eq!(
                    (
                        collection.d_i_atomic_num[row_index],
                        collection.d_j_atomic_num[row_index],
                    ),
                    EXPECTED_KEYS[row_index]
                );
                assert_eq!(
                    collection.d_params[row_index].r0.to_bits(),
                    EXPECTED_VALUE_BITS[row_index].0
                );
                assert_eq!(
                    collection.d_params[row_index].kb.to_bits(),
                    EXPECTED_VALUE_BITS[row_index].1
                );
            }

            for (query_index, &(atomic_num, nbr_atomic_num, expected_index)) in
                QUERIES.iter().enumerate()
            {
                assert!(actual_get_calls < 32, "unexpected extra Bndk lookup call");
                let before = collection.clone();
                let actual = collection.get(atomic_num, nbr_atomic_num);
                assert_eq!(
                    collection, before,
                    "lookup mutated {format_index}/{query_index}"
                );
                assert_eq!(source.as_bytes(), source_before.as_slice());

                let expected_index = if format_index == 3 && expected_index == Some(3) {
                    None
                } else {
                    expected_index
                };
                match expected_index {
                    Some(row_index) => {
                        let actual = actual.expect("frozen Bndk row must be present");
                        let expected = &collection.d_params[row_index];
                        assert!(std::ptr::eq(actual, expected));
                        assert_eq!(actual.r0.to_bits(), EXPECTED_VALUE_BITS[row_index].0);
                        assert_eq!(actual.kb.to_bits(), EXPECTED_VALUE_BITS[row_index].1);
                    }
                    None => assert!(actual.is_none()),
                }
                actual_get_calls += 1;
            }
        }
        assert_eq!(actual_get_calls, 32);

        let cast = MmffBndkCollection::from_text("257\t258\t1.5\t1.25\n")
            .expect("unsigned stored keys narrow to their source u8 carriers");
        assert_eq!(cast.d_i_atomic_num, [1]);
        assert_eq!(cast.d_j_atomic_num, [2]);
        let stored_cast = cast.get(1, 2).expect("narrowed stored keys must match");
        assert_eq!(stored_cast.r0.to_bits(), 0x3ff8_0000_0000_0000);
        assert_eq!(stored_cast.kb.to_bits(), 0x3ff4_0000_0000_0000);
        assert!(cast.get(257, 258).is_none());

        assert_eq!(
            MmffBndkCollection::from_text("\n").unwrap_err(),
            MmffParamParseError {
                table: MmffParamTable::Bndk,
                line: 1,
                column: 0,
                cause: MmffParamParseCause::EmptyProcessedLine,
            }
        );
        assert_eq!(
            MmffBndkCollection::from_text("1\t2\t1.5\n").unwrap_err(),
            MmffParamParseError {
                table: MmffParamTable::Bndk,
                line: 1,
                column: 3,
                cause: MmffParamParseCause::MissingToken,
            }
        );
        assert_eq!(
            MmffBndkCollection::from_text("X\t2\t1.5\t1.25\n").unwrap_err(),
            MmffParamParseError {
                table: MmffParamTable::Bndk,
                line: 1,
                column: 0,
                cause: MmffParamParseCause::InvalidUnsigned { cell: "X".into() },
            }
        );
        assert_eq!(
            MmffBndkCollection::from_text("1\t2\tbad\t1.25\n").unwrap_err(),
            MmffParamParseError {
                table: MmffParamTable::Bndk,
                line: 1,
                column: 2,
                cause: MmffParamParseCause::InvalidFloat { cell: "bad".into() },
            }
        );
    }

    #[test]
    fn mmff_bond_defaults_assets_sentinels_and_shared_borrows() {
        const DEFAULT_BOND_TEXT: &str = include_str!("default_bond.tsv");
        const DEFAULT_BNDK_TEXT: &str = include_str!("default_bndk.tsv");

        assert_eq!(
            sha256_for_fixed_asset_test(b"abc"),
            [
                0xba, 0x78, 0x16, 0xbf, 0x8f, 0x01, 0xcf, 0xea, 0x41, 0x41, 0x40, 0xde, 0x5d, 0xae,
                0x22, 0x23, 0xb0, 0x03, 0x61, 0xa3, 0x96, 0x17, 0x7a, 0x9c, 0xb4, 0x10, 0xff, 0x61,
                0xf2, 0x00, 0x15, 0xad,
            ]
        );
        assert_eq!(DEFAULT_BOND_TEXT.len(), 12_172);
        assert_eq!(
            sha256_for_fixed_asset_test(DEFAULT_BOND_TEXT.as_bytes()),
            [
                0x2d, 0x6a, 0xaa, 0x67, 0x0c, 0x42, 0x5f, 0x98, 0xec, 0x9a, 0x8b, 0x18, 0xa3, 0xcd,
                0x6e, 0x5a, 0xef, 0x50, 0xf1, 0x4a, 0x31, 0x2a, 0xe8, 0xd9, 0x2e, 0x0a, 0xd7, 0xfa,
                0x79, 0xca, 0xa2, 0xbb,
            ]
        );
        assert_eq!(DEFAULT_BNDK_TEXT.len(), 1_422);
        assert_eq!(
            sha256_for_fixed_asset_test(DEFAULT_BNDK_TEXT.as_bytes()),
            [
                0xb8, 0x76, 0x7e, 0xf0, 0x84, 0xb9, 0x0c, 0x0c, 0x99, 0xc6, 0x70, 0xff, 0x2b, 0x21,
                0x00, 0x76, 0x57, 0x13, 0x52, 0x29, 0x51, 0xfc, 0xd8, 0x2f, 0x61, 0xc4, 0x2b, 0x57,
                0xd7, 0xbc, 0x59, 0x55,
            ]
        );
        assert_eq!(
            DEFAULT_BOND_TEXT
                .lines()
                .filter(|line| !line.starts_with('*'))
                .count(),
            493
        );
        assert_eq!(
            DEFAULT_BNDK_TEXT
                .lines()
                .filter(|line| !line.starts_with('*'))
                .count(),
            58
        );

        let default_bond_bytes = DEFAULT_BOND_TEXT.as_bytes().to_vec();
        let default_bndk_bytes = DEFAULT_BNDK_TEXT.as_bytes().to_vec();

        let custom_bond = MmffBondCollection::from_text("0\t1\t1\t9\t10\n")
            .expect("independent Bond construction must parse");
        let custom_bndk = MmffBndkCollection::from_text("1\t6\t9\t10\n")
            .expect("independent Bndk construction must parse");
        assert_eq!(custom_bond.d_params[0].kb.to_bits(), 9.0_f64.to_bits());
        assert_eq!(custom_bndk.d_params[0].r0.to_bits(), 9.0_f64.to_bits());

        let bond = default_mmff_bond().expect("fixed Bond asset must parse");
        let bndk = default_mmff_bndk().expect("fixed Bndk asset must parse");
        assert_eq!(bond.d_params.len(), 493);
        assert_eq!(bndk.d_params.len(), 58);

        assert_eq!(
            (
                bond.d_bond_type[0],
                bond.d_i_atom_type[0],
                bond.d_j_atom_type[0],
            ),
            (0, 1, 1)
        );
        assert_eq!(bond.d_params[0].kb.to_bits(), 0x4011_0831_26e9_78d5);
        assert_eq!(bond.d_params[0].r0.to_bits(), 0x3ff8_20c4_9ba5_e354);
        assert_eq!(
            (
                bond.d_bond_type[246],
                bond.d_i_atom_type[246],
                bond.d_j_atom_type[246],
            ),
            (0, 10, 37)
        );
        assert_eq!(bond.d_params[246].kb.to_bits(), 0x4015_ed91_6872_b021);
        assert_eq!(bond.d_params[246].r0.to_bits(), 0x3ff6_51eb_851e_b852);
        assert_eq!(
            (
                bond.d_bond_type[492],
                bond.d_i_atom_type[492],
                bond.d_j_atom_type[492],
            ),
            (0, 80, 81)
        );
        assert_eq!(bond.d_params[492].kb.to_bits(), 0x4020_7958_1062_4dd3);
        assert_eq!(bond.d_params[492].r0.to_bits(), 0x3ff5_5c28_f5c2_8f5c);

        assert_eq!((bndk.d_i_atomic_num[0], bndk.d_j_atomic_num[0]), (1, 6));
        assert_eq!(bndk.d_params[0].r0.to_bits(), 0x3ff1_5810_624d_d2f2);
        assert_eq!(bndk.d_params[0].kb.to_bits(), 0x4014_9999_9999_999a);
        assert_eq!((bndk.d_i_atomic_num[29], bndk.d_j_atomic_num[29]), (8, 8));
        assert_eq!(bndk.d_params[29].r0.to_bits(), 0x3ff7_ae14_7ae1_47ae);
        assert_eq!(bndk.d_params[29].kb.to_bits(), 0x400c_cccc_cccc_cccd);
        assert_eq!((bndk.d_i_atomic_num[57], bndk.d_j_atomic_num[57]), (53, 53));
        assert_eq!(bndk.d_params[57].r0.to_bits(), 0x4005_5c28_f5c2_8f5c);
        assert_eq!(bndk.d_params[57].kb.to_bits(), 0x3ff9_9999_9999_999a);

        let bond_row = bond.get(0, 37, 10).expect("reversed Bond sentinel exists");
        let bndk_row = bndk.get(6, 1).expect("reversed Bndk sentinel exists");
        assert!(std::ptr::eq(bond_row, &bond.d_params[246]));
        assert!(std::ptr::eq(bndk_row, &bndk.d_params[0]));
        let bond_before = bond.clone();
        let bndk_before = bndk.clone();

        let (bond_calls, bndk_calls) = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..8)
                .map(|_| {
                    scope.spawn(move || {
                        let mut bond_calls = 0;
                        let mut bndk_calls = 0;
                        for _ in 0..16 {
                            let cached_bond =
                                default_mmff_bond().expect("cached Bond result remains valid");
                            assert!(std::ptr::eq(cached_bond, bond));
                            assert!(std::ptr::eq(
                                cached_bond
                                    .get(0, 37, 10)
                                    .expect("Bond sentinel remains present"),
                                bond_row
                            ));
                            bond_calls += 1;

                            let cached_bndk =
                                default_mmff_bndk().expect("cached Bndk result remains valid");
                            assert!(std::ptr::eq(cached_bndk, bndk));
                            assert!(std::ptr::eq(
                                cached_bndk
                                    .get(6, 1)
                                    .expect("Bndk sentinel remains present"),
                                bndk_row
                            ));
                            bndk_calls += 1;
                        }
                        (bond_calls, bndk_calls)
                    })
                })
                .collect();

            handles
                .into_iter()
                .fold((0, 0), |(bond_total, bndk_total), handle| {
                    let (bond_calls, bndk_calls) = handle.join().expect("scoped borrow worker");
                    (bond_total + bond_calls, bndk_total + bndk_calls)
                })
        });
        assert_eq!((bond_calls, bndk_calls), (128, 128));
        assert_eq!(*bond, bond_before);
        assert_eq!(*bndk, bndk_before);
        assert_eq!(DEFAULT_BOND_TEXT.as_bytes(), default_bond_bytes.as_slice());
        assert_eq!(DEFAULT_BNDK_TEXT.as_bytes(), default_bndk_bytes.as_slice());
    }

    fn mmff_sb_stbn_value_bits(value: &MmffStbn) -> [u64; 2] {
        [value.kba_ijk.to_bits(), value.kba_kji.to_bits()]
    }

    fn mmff_sb_stbn_storage_addresses(collection: &MmffStbnCollection) -> [usize; 5] {
        [
            collection.d_i_atom_type.as_ptr() as usize,
            collection.d_j_atom_type.as_ptr() as usize,
            collection.d_k_atom_type.as_ptr() as usize,
            collection.d_stretch_bend_type.as_ptr() as usize,
            collection.d_params.as_ptr() as usize,
        ]
    }

    fn assert_mmff_sb_stbn_state_unchanged(
        source: &str,
        input_bytes_before: &[u8],
        collection: &MmffStbnCollection,
        keys_before: &[Vec<u8>; 4],
        values_before: &[[u64; 2]],
        addresses_before: [usize; 5],
    ) {
        assert_eq!(source.as_bytes(), input_bytes_before);
        assert_eq!(
            collection.d_i_atom_type.as_slice(),
            keys_before[0].as_slice()
        );
        assert_eq!(
            collection.d_j_atom_type.as_slice(),
            keys_before[1].as_slice()
        );
        assert_eq!(
            collection.d_k_atom_type.as_slice(),
            keys_before[2].as_slice()
        );
        assert_eq!(
            collection.d_stretch_bend_type.as_slice(),
            keys_before[3].as_slice()
        );
        assert_eq!(
            collection
                .d_params
                .iter()
                .map(mmff_sb_stbn_value_bits)
                .collect::<Vec<_>>(),
            values_before
        );
        assert_eq!(mmff_sb_stbn_storage_addresses(collection), addresses_before);
    }

    #[test]
    fn mmff_sb_stbn_lookup_matrix_has_exact_48_calls() {
        const SOURCE: &str = concat!(
            "*fixed\n",
            "0\t1\t2\t1\t1.25\t2.5\n",
            "1\t1\t2\t1\t3.75\t4\n",
            "1\t1\t2\t1\t9\t10\n",
            "0\t1\t2\t3\t0.5\t0.75\n",
            "2\t1\t2\t3\t0.125\t0.25\n",
            "0\t3\t2\t3\t5\t6\n",
            "0\t1\t4\t3\t7\t8\n",
        );
        const EXPECTED_VALUE_BITS: [[u64; 2]; 7] = [
            [0x3ff4_0000_0000_0000, 0x4004_0000_0000_0000],
            [0x400e_0000_0000_0000, 0x4010_0000_0000_0000],
            [0x4022_0000_0000_0000, 0x4024_0000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x3fe8_0000_0000_0000],
            [0x3fc0_0000_0000_0000, 0x3fd0_0000_0000_0000],
            [0x4014_0000_0000_0000, 0x4018_0000_0000_0000],
            [0x401c_0000_0000_0000, 0x4020_0000_0000_0000],
        ];
        const QUERIES: [(u32, u32, u32, u32, u32, u32, bool, Option<usize>); 12] = [
            (0, 0, 0, 1, 2, 1, false, Some(0)),
            (0, 0, 1, 1, 2, 1, true, Some(0)),
            (1, 2, 1, 1, 2, 1, false, Some(1)),
            (1, 1, 2, 1, 2, 1, true, Some(1)),
            (0, 0, 0, 1, 2, 3, false, Some(3)),
            (0, 0, 0, 3, 2, 1, true, Some(3)),
            (2, 0, 0, 1, 2, 3, false, Some(4)),
            (0, 1, 2, 3, 2, 3, true, Some(5)),
            (0, 0, 0, 1, 4, 3, false, Some(6)),
            (1, 0, 0, 1, 2, 3, false, None),
            (0, 0, 0, 1, 9, 3, false, None),
            (0, 0, 0, 257, 2, 3, true, None),
        ];
        const EXPECTED_I: [u8; 7] = [1, 1, 1, 1, 1, 3, 1];
        const EXPECTED_J: [u8; 7] = [2, 2, 2, 2, 2, 2, 4];
        const EXPECTED_K: [u8; 7] = [1, 1, 1, 3, 3, 3, 3];
        const EXPECTED_STRETCH: [u8; 7] = [0, 1, 1, 0, 2, 0, 0];

        let crlf = SOURCE.replace('\n', "\r\n");
        let doubled_tabs = SOURCE.replace('\t', "\t\t");
        let unterminated = SOURCE
            .strip_suffix('\n')
            .expect("fixed Stbn table ends in LF")
            .to_owned();
        let formats = [
            (SOURCE.to_owned(), false),
            (crlf, false),
            (doubled_tabs, false),
            (unterminated, true),
        ];

        let mut total_get_calls = 0;
        for (source, drops_final_row) in &formats {
            let collection =
                MmffStbnCollection::from_text(source).expect("fixed Stbn source format parses");
            let expected_rows = if *drops_final_row { 6 } else { 7 };
            assert_eq!(collection.d_i_atom_type.len(), expected_rows);
            assert_eq!(collection.d_j_atom_type.len(), expected_rows);
            assert_eq!(collection.d_k_atom_type.len(), expected_rows);
            assert_eq!(collection.d_stretch_bend_type.len(), expected_rows);
            assert_eq!(collection.d_params.len(), expected_rows);
            assert_eq!(
                collection.d_i_atom_type.as_slice(),
                &EXPECTED_I[..expected_rows]
            );
            assert_eq!(
                collection.d_j_atom_type.as_slice(),
                &EXPECTED_J[..expected_rows]
            );
            assert_eq!(
                collection.d_k_atom_type.as_slice(),
                &EXPECTED_K[..expected_rows]
            );
            assert_eq!(
                collection.d_stretch_bend_type.as_slice(),
                &EXPECTED_STRETCH[..expected_rows]
            );
            assert_eq!(
                collection
                    .d_params
                    .iter()
                    .map(mmff_sb_stbn_value_bits)
                    .collect::<Vec<_>>(),
                &EXPECTED_VALUE_BITS[..expected_rows]
            );

            let mut format_get_calls = 0;
            for (stretch, bond1, bond2, i, j, k, expected_swap, expected_row) in QUERIES {
                let expected_row = if *drops_final_row && expected_row == Some(6) {
                    None
                } else {
                    expected_row
                };
                let input_bytes_before = source.as_bytes().to_vec();
                let keys_before = [
                    collection.d_i_atom_type.clone(),
                    collection.d_j_atom_type.clone(),
                    collection.d_k_atom_type.clone(),
                    collection.d_stretch_bend_type.clone(),
                ];
                let values_before: Vec<_> = collection
                    .d_params
                    .iter()
                    .map(mmff_sb_stbn_value_bits)
                    .collect();
                let addresses_before = mmff_sb_stbn_storage_addresses(&collection);
                let query_before = (stretch, bond1, bond2, i, j, k);

                let (actual_swap, actual) = collection.get(stretch, bond1, bond2, i, j, k);
                format_get_calls += 1;
                total_get_calls += 1;

                assert_eq!((stretch, bond1, bond2, i, j, k), query_before);
                assert_mmff_sb_stbn_state_unchanged(
                    source,
                    &input_bytes_before,
                    &collection,
                    &keys_before,
                    &values_before,
                    addresses_before,
                );
                assert_eq!(actual_swap, expected_swap);
                if let Some(row_index) = expected_row {
                    let actual = actual.expect("fixed Stbn lookup hit");
                    assert!(
                        std::ptr::eq(actual, &collection.d_params[row_index]),
                        "query {query_before:?} should borrow literal row {row_index}"
                    );
                    assert_eq!(
                        mmff_sb_stbn_value_bits(actual),
                        EXPECTED_VALUE_BITS[row_index]
                    );
                } else {
                    assert!(actual.is_none(), "unexpected Stbn query hit");
                }
            }

            assert_eq!(format_get_calls, 12);
        }
        assert_eq!(total_get_calls, 48);
    }

    #[test]
    fn mmff_sb_stbn_cast_controls_preserve_full_width_queries_and_signed_zero() {
        let source = "257\t258\t259\t260\t-0\t1.25\n";
        let collection = MmffStbnCollection::from_text(source).expect("fixed cast row parses");
        assert_eq!(collection.d_stretch_bend_type, [1]);
        assert_eq!(collection.d_i_atom_type, [2]);
        assert_eq!(collection.d_j_atom_type, [3]);
        assert_eq!(collection.d_k_atom_type, [4]);
        const EXPECTED_BITS: [u64; 2] = [0x8000_0000_0000_0000, 0x3ff4_0000_0000_0000];
        assert_eq!(
            mmff_sb_stbn_value_bits(&collection.d_params[0]),
            EXPECTED_BITS
        );
        const QUERIES: [(u32, u32, u32, u32, u32, u32, bool, Option<usize>); 5] = [
            (1, 0, 0, 2, 3, 4, false, Some(0)),
            (257, 0, 0, 2, 3, 4, false, None),
            (1, 0, 0, 258, 3, 4, true, None),
            (1, 0, 0, 2, 259, 4, false, None),
            (1, 0, 0, 2, 3, 260, false, None),
        ];

        let mut get_calls = 0;
        for (stretch, bond1, bond2, i, j, k, expected_swap, expected_row) in QUERIES {
            let input_bytes_before = source.as_bytes().to_vec();
            let keys_before = [
                collection.d_i_atom_type.clone(),
                collection.d_j_atom_type.clone(),
                collection.d_k_atom_type.clone(),
                collection.d_stretch_bend_type.clone(),
            ];
            let values_before: Vec<_> = collection
                .d_params
                .iter()
                .map(mmff_sb_stbn_value_bits)
                .collect();
            let addresses_before = mmff_sb_stbn_storage_addresses(&collection);
            let query_before = (stretch, bond1, bond2, i, j, k);

            let (actual_swap, actual) = collection.get(stretch, bond1, bond2, i, j, k);
            get_calls += 1;

            assert_eq!((stretch, bond1, bond2, i, j, k), query_before);
            assert_mmff_sb_stbn_state_unchanged(
                source,
                &input_bytes_before,
                &collection,
                &keys_before,
                &values_before,
                addresses_before,
            );
            assert_eq!(actual_swap, expected_swap);
            if let Some(row_index) = expected_row {
                let actual = actual.expect("small query observes cast stored keys");
                assert!(std::ptr::eq(actual, &collection.d_params[row_index]));
                assert_eq!(mmff_sb_stbn_value_bits(actual), EXPECTED_BITS);
            } else {
                assert!(actual.is_none(), "full-width query keys are not narrowed");
            }
        }
        assert_eq!(get_calls, 5);
    }

    #[test]
    fn mmff_sb_stbn_constructor_reports_frozen_typed_errors() {
        let valid_cells = ["0", "1", "2", "1", "1.25", "2.5"];
        let mut error_constructors = 0;

        for column in 0..valid_cells.len() {
            let mut cells = valid_cells;
            cells[column] = "X";
            let source = format!("*fixed\n{}\n", cells.join("\t"));
            let input_bytes_before = source.as_bytes().to_vec();
            let result = MmffStbnCollection::from_text(&source);
            error_constructors += 1;
            assert_eq!(source.as_bytes(), input_bytes_before.as_slice());

            let error = result.expect_err("each consumed invalid Stbn cell is rejected");
            assert_eq!(error.table, MmffParamTable::Stbn);
            assert_eq!(error.line, 2);
            assert_eq!(error.column, column);
            let expected_cause = if column < 4 {
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                }
            } else {
                MmffParamParseCause::InvalidFloat {
                    cell: "X".to_owned(),
                }
            };
            assert_eq!(&error.cause, &expected_cause);
            if column == 0 {
                assert_eq!(
                    error.to_string(),
                    "MMFF Stbn parse error at physical line 2 token 0: invalid unsigned integer \"X\""
                );
            }
        }

        let missing_final_cell = "*fixed\n0\t1\t2\t1\t1.25\n";
        let input_bytes_before = missing_final_cell.as_bytes().to_vec();
        let result = MmffStbnCollection::from_text(missing_final_cell);
        error_constructors += 1;
        assert_eq!(missing_final_cell.as_bytes(), input_bytes_before.as_slice());
        let error = result.expect_err("missing kbaKJI is reported");
        assert_eq!(error.table, MmffParamTable::Stbn);
        assert_eq!(error.line, 2);
        assert_eq!(error.column, 5);
        assert_eq!(error.cause, MmffParamParseCause::MissingToken);
        assert_eq!(
            error.to_string(),
            "MMFF Stbn parse error at physical line 2 token 5: missing token"
        );

        let empty_processed_line = "*fixed\n\n";
        let input_bytes_before = empty_processed_line.as_bytes().to_vec();
        let result = MmffStbnCollection::from_text(empty_processed_line);
        error_constructors += 1;
        assert_eq!(
            empty_processed_line.as_bytes(),
            input_bytes_before.as_slice()
        );
        let error = result.expect_err("empty processed line is source-undefined");
        assert_eq!(error.table, MmffParamTable::Stbn);
        assert_eq!(error.line, 2);
        assert_eq!(error.column, 0);
        assert_eq!(error.cause, MmffParamParseCause::EmptyProcessedLine);
        assert_eq!(
            error.to_string(),
            "MMFF Stbn parse error at physical line 2 token 0: empty processed line"
        );
        assert_eq!(error_constructors, 8);

        let unused_trailing_token = "*fixed\n0\t1\t2\t1\t1.25\t2.5\tX\n";
        let input_bytes_before = unused_trailing_token.as_bytes().to_vec();
        let collection = MmffStbnCollection::from_text(unused_trailing_token)
            .expect("source ignores the unused trailing token");
        assert_eq!(
            unused_trailing_token.as_bytes(),
            input_bytes_before.as_slice()
        );
        assert_eq!(collection.d_stretch_bend_type, [0]);
        assert_eq!(collection.d_i_atom_type, [1]);
        assert_eq!(collection.d_j_atom_type, [2]);
        assert_eq!(collection.d_k_atom_type, [1]);
        assert_eq!(
            mmff_sb_stbn_value_bits(&collection.d_params[0]),
            [0x3ff4_0000_0000_0000, 0x4004_0000_0000_0000]
        );
    }

    type MmffSbDfsbSnapshot = Vec<(u32, u32, u32, [u64; 2], usize)>;

    fn mmff_sb_dfsb_snapshot(collection: &MmffDfsbCollection) -> MmffSbDfsbSnapshot {
        collection
            .d_params
            .iter()
            .flat_map(|(&row1, row2_map)| {
                row2_map.iter().flat_map(move |(&row2, row3_map)| {
                    row3_map.iter().map(move |(&row3, value)| {
                        (
                            row1,
                            row2,
                            row3,
                            mmff_sb_stbn_value_bits(value),
                            value as *const MmffStbn as usize,
                        )
                    })
                })
            })
            .collect()
    }

    fn mmff_sb_dfsb_map_address(collection: &MmffDfsbCollection) -> usize {
        &collection.d_params as *const _ as usize
    }

    fn assert_mmff_sb_dfsb_state_unchanged(
        source: &str,
        input_bytes_before: &[u8],
        collection: &MmffDfsbCollection,
        map_address_before: usize,
        snapshot_before: &MmffSbDfsbSnapshot,
    ) {
        assert_eq!(source.as_bytes(), input_bytes_before);
        assert_eq!(mmff_sb_dfsb_map_address(collection), map_address_before);
        assert_eq!(
            mmff_sb_dfsb_snapshot(collection).as_slice(),
            snapshot_before.as_slice()
        );
    }

    #[test]
    fn mmff_sb_dfsb_lookup_matrix_has_exact_48_calls() {
        const SOURCE: &str = concat!(
            "*fixed\n",
            "1\t2\t1\t1.25\t2.5\n",
            "1\t2\t1\t9\t10\n",
            "1\t2\t3\t0.5\t0.75\n",
            "3\t2\t3\t5\t6\n",
            "1\t4\t3\t7\t8\n",
            "257\t258\t259\t-0\t1.25\n",
        );
        const EXPECTED_LEAVES: [((u32, u32, u32), [u64; 2]); 5] = [
            ((1, 2, 1), [0x4022_0000_0000_0000, 0x4024_0000_0000_0000]),
            ((1, 2, 3), [0x3fe0_0000_0000_0000, 0x3fe8_0000_0000_0000]),
            ((1, 4, 3), [0x401c_0000_0000_0000, 0x4020_0000_0000_0000]),
            ((3, 2, 3), [0x4014_0000_0000_0000, 0x4018_0000_0000_0000]),
            (
                (257, 258, 259),
                [0x8000_0000_0000_0000, 0x3ff4_0000_0000_0000],
            ),
        ];
        const QUERIES: [(u32, u32, u32, bool, Option<((u32, u32, u32), [u64; 2])>); 12] = [
            (
                1,
                2,
                1,
                false,
                Some(((1, 2, 1), [0x4022_0000_0000_0000, 0x4024_0000_0000_0000])),
            ),
            (
                1,
                2,
                3,
                false,
                Some(((1, 2, 3), [0x3fe0_0000_0000_0000, 0x3fe8_0000_0000_0000])),
            ),
            (
                3,
                2,
                1,
                true,
                Some(((1, 2, 3), [0x3fe0_0000_0000_0000, 0x3fe8_0000_0000_0000])),
            ),
            (
                3,
                2,
                3,
                false,
                Some(((3, 2, 3), [0x4014_0000_0000_0000, 0x4018_0000_0000_0000])),
            ),
            (
                1,
                4,
                3,
                false,
                Some(((1, 4, 3), [0x401c_0000_0000_0000, 0x4020_0000_0000_0000])),
            ),
            (
                3,
                4,
                1,
                true,
                Some(((1, 4, 3), [0x401c_0000_0000_0000, 0x4020_0000_0000_0000])),
            ),
            (
                257,
                258,
                259,
                false,
                Some((
                    (257, 258, 259),
                    [0x8000_0000_0000_0000, 0x3ff4_0000_0000_0000],
                )),
            ),
            (
                259,
                258,
                257,
                true,
                Some((
                    (257, 258, 259),
                    [0x8000_0000_0000_0000, 0x3ff4_0000_0000_0000],
                )),
            ),
            (1, 258, 3, false, None),
            (0, 2, 3, false, None),
            (9, 2, 3, true, None),
            (u32::MAX, 2, 1, true, None),
        ];

        let crlf = SOURCE.replace('\n', "\r\n");
        let doubled_tabs = SOURCE.replace('\t', "\t\t");
        let unterminated = SOURCE
            .strip_suffix('\n')
            .expect("fixed Dfsb table ends in LF")
            .to_owned();
        let formats = [
            (SOURCE.to_owned(), false),
            (crlf, false),
            (doubled_tabs, false),
            (unterminated, true),
        ];

        let mut total_get_calls = 0;
        for (source, drops_final_row) in &formats {
            let collection =
                MmffDfsbCollection::from_text(source).expect("fixed Dfsb source format parses");
            let expected_leaf_count = if *drops_final_row { 4 } else { 5 };
            let initial_snapshot = mmff_sb_dfsb_snapshot(&collection);
            assert_eq!(initial_snapshot.len(), expected_leaf_count);
            let expected_rows: Vec<_> = EXPECTED_LEAVES[..expected_leaf_count]
                .iter()
                .map(|((row1, row2, row3), bits)| (*row1, *row2, *row3, *bits))
                .collect();
            let actual_rows: Vec<_> = initial_snapshot
                .iter()
                .map(|(row1, row2, row3, bits, _)| (*row1, *row2, *row3, *bits))
                .collect();
            assert_eq!(actual_rows, expected_rows);

            let mut format_get_calls = 0;
            for (row1, row2, row3, expected_swap, expected_leaf) in QUERIES {
                let expected_leaf = if *drops_final_row
                    && expected_leaf.is_some_and(|(key, _)| key == (257, 258, 259))
                {
                    None
                } else {
                    expected_leaf
                };
                let input_bytes_before = source.as_bytes().to_vec();
                let map_address_before = mmff_sb_dfsb_map_address(&collection);
                let snapshot_before = mmff_sb_dfsb_snapshot(&collection);
                let query_before = (row1, row2, row3);

                let (actual_swap, actual) = collection.get(row1, row2, row3);
                format_get_calls += 1;
                total_get_calls += 1;

                assert_eq!((row1, row2, row3), query_before);
                assert_mmff_sb_dfsb_state_unchanged(
                    source,
                    &input_bytes_before,
                    &collection,
                    map_address_before,
                    &snapshot_before,
                );
                assert_eq!(actual_swap, expected_swap);
                if let Some((key, expected_bits)) = expected_leaf {
                    let actual = actual.expect("fixed Dfsb lookup hit");
                    let expected_value = collection
                        .d_params
                        .get(&key.0)
                        .and_then(|row2_map| row2_map.get(&key.1))
                        .and_then(|row3_map| row3_map.get(&key.2))
                        .expect("fixed expected Dfsb key exists");
                    assert!(std::ptr::eq(actual, expected_value));
                    assert_eq!(mmff_sb_stbn_value_bits(actual), expected_bits);
                } else {
                    assert!(actual.is_none(), "unexpected Dfsb query hit");
                }
            }

            assert_eq!(format_get_calls, 12);
        }
        assert_eq!(total_get_calls, 48);
    }

    #[test]
    fn mmff_sb_dfsb_noncanonical_storage_stays_unreachable() {
        let source = "3\t2\t1\t5\t6\n";
        let collection = MmffDfsbCollection::from_text(source)
            .expect("fixed noncanonical Dfsb source row parses");
        let expected_snapshot = mmff_sb_dfsb_snapshot(&collection);
        assert_eq!(expected_snapshot.len(), 1);
        assert_eq!(
            expected_snapshot
                .iter()
                .map(|(row1, row2, row3, bits, _)| (*row1, *row2, *row3, *bits))
                .collect::<Vec<_>>(),
            [(3, 2, 1, [0x4014_0000_0000_0000, 0x4018_0000_0000_0000])]
        );

        const QUERIES: [(u32, u32, u32, bool); 3] =
            [(3, 2, 1, true), (1, 2, 3, false), (1, 2, 1, false)];
        let mut get_calls = 0;
        for (row1, row2, row3, expected_swap) in QUERIES {
            let input_bytes_before = source.as_bytes().to_vec();
            let map_address_before = mmff_sb_dfsb_map_address(&collection);
            let snapshot_before = mmff_sb_dfsb_snapshot(&collection);
            let query_before = (row1, row2, row3);

            let (actual_swap, actual) = collection.get(row1, row2, row3);
            get_calls += 1;

            assert_eq!((row1, row2, row3), query_before);
            assert_mmff_sb_dfsb_state_unchanged(
                source,
                &input_bytes_before,
                &collection,
                map_address_before,
                &snapshot_before,
            );
            assert_eq!(actual_swap, expected_swap);
            assert!(actual.is_none(), "noncanonical stored key is not rewritten");
        }
        assert_eq!(get_calls, 3);
    }

    #[test]
    fn mmff_sb_dfsb_constructor_reports_frozen_typed_errors() {
        let valid_cells = ["1", "2", "1", "1.25", "2.5"];
        let mut error_constructors = 0;

        for column in 0..valid_cells.len() {
            let mut cells = valid_cells;
            cells[column] = "X";
            let source = format!("*fixed\n{}\n", cells.join("\t"));
            let input_bytes_before = source.as_bytes().to_vec();
            let result = MmffDfsbCollection::from_text(&source);
            error_constructors += 1;
            assert_eq!(source.as_bytes(), input_bytes_before.as_slice());

            let error = result.expect_err("each consumed invalid Dfsb cell is rejected");
            assert_eq!(error.table, MmffParamTable::Dfsb);
            assert_eq!(error.line, 2);
            assert_eq!(error.column, column);
            let expected_cause = if column < 3 {
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                }
            } else {
                MmffParamParseCause::InvalidFloat {
                    cell: "X".to_owned(),
                }
            };
            assert_eq!(&error.cause, &expected_cause);
            if column == 0 {
                assert_eq!(
                    error.to_string(),
                    "MMFF Dfsb parse error at physical line 2 token 0: invalid unsigned integer \"X\""
                );
            }
        }

        let missing_final_cell = "*fixed\n1\t2\t3\t5\n";
        let input_bytes_before = missing_final_cell.as_bytes().to_vec();
        let result = MmffDfsbCollection::from_text(missing_final_cell);
        error_constructors += 1;
        assert_eq!(missing_final_cell.as_bytes(), input_bytes_before.as_slice());
        let error = result.expect_err("missing kbaKJI is reported");
        assert_eq!(error.table, MmffParamTable::Dfsb);
        assert_eq!(error.line, 2);
        assert_eq!(error.column, 4);
        assert_eq!(error.cause, MmffParamParseCause::MissingToken);
        assert_eq!(
            error.to_string(),
            "MMFF Dfsb parse error at physical line 2 token 4: missing token"
        );

        let empty_processed_line = "*fixed\n\n";
        let input_bytes_before = empty_processed_line.as_bytes().to_vec();
        let result = MmffDfsbCollection::from_text(empty_processed_line);
        error_constructors += 1;
        assert_eq!(
            empty_processed_line.as_bytes(),
            input_bytes_before.as_slice()
        );
        let error = result.expect_err("empty processed line is source-undefined");
        assert_eq!(error.table, MmffParamTable::Dfsb);
        assert_eq!(error.line, 2);
        assert_eq!(error.column, 0);
        assert_eq!(error.cause, MmffParamParseCause::EmptyProcessedLine);
        assert_eq!(
            error.to_string(),
            "MMFF Dfsb parse error at physical line 2 token 0: empty processed line"
        );
        assert_eq!(error_constructors, 7);

        let unused_trailing_token = "*fixed\n1\t2\t3\t5\t6\tX\n";
        let input_bytes_before = unused_trailing_token.as_bytes().to_vec();
        let collection = MmffDfsbCollection::from_text(unused_trailing_token)
            .expect("source ignores the unused trailing token");
        assert_eq!(
            unused_trailing_token.as_bytes(),
            input_bytes_before.as_slice()
        );
        assert_eq!(
            mmff_sb_dfsb_snapshot(&collection)
                .iter()
                .map(|(row1, row2, row3, bits, _)| (*row1, *row2, *row3, *bits))
                .collect::<Vec<_>>(),
            [(1, 2, 3, [0x4014_0000_0000_0000, 0x4018_0000_0000_0000])]
        );
    }

    fn mmff_emp_herschbach_value_bits(value: &MmffHerschbachLaurie) -> [u64; 3] {
        [
            value.a_ij.to_bits(),
            value.d_ij.to_bits(),
            value.dp_ij.to_bits(),
        ]
    }

    #[test]
    fn mmff_emp_herschbach_lookup_matrix_has_exact_32_calls() {
        const SOURCE: &str = concat!(
            "*fixed\n",
            "1\t1\t1.25\t0.5\t2\n",
            "1\t1\t9\t9.5\t10\n",
            "1\t2\t2.5\t1.5\t3\n",
            "2\t2\t3.75\t2.5\t4\n",
        );
        const EXPECTED_VALUE_BITS: [[u64; 3]; 4] = [
            [
                0x3ff4_0000_0000_0000,
                0x3fe0_0000_0000_0000,
                0x4000_0000_0000_0000,
            ],
            [
                0x4022_0000_0000_0000,
                0x4023_0000_0000_0000,
                0x4024_0000_0000_0000,
            ],
            [
                0x4004_0000_0000_0000,
                0x3ff8_0000_0000_0000,
                0x4008_0000_0000_0000,
            ],
            [
                0x400e_0000_0000_0000,
                0x4004_0000_0000_0000,
                0x4010_0000_0000_0000,
            ],
        ];
        const QUERIES: [(i32, i32, Option<usize>); 8] = [
            (1, 1, Some(0)),
            (1, 2, Some(2)),
            (2, 1, Some(2)),
            (2, 2, Some(3)),
            (0, 1, None),
            (257, 1, None),
            (-1, 1, None),
            (i32::MIN, 1, None),
        ];

        let crlf = SOURCE.replace('\n', "\r\n");
        let doubled_tabs = SOURCE.replace('\t', "\t\t");
        let unterminated = SOURCE
            .strip_suffix('\n')
            .expect("fixed table ends in LF")
            .to_owned();
        let formats = [
            (SOURCE.to_owned(), false),
            (crlf, false),
            (doubled_tabs, false),
            (unterminated, true),
        ];

        let mut total_get_calls = 0;
        for (source, drops_final_row) in &formats {
            let collection = MmffHerschbachLaurieCollection::from_text(source)
                .expect("fixed Herschbach format parses");
            let expected_rows = if *drops_final_row { 3 } else { 4 };
            assert_eq!(collection.d_i_row.len(), expected_rows);
            assert_eq!(collection.d_j_row.len(), expected_rows);
            assert_eq!(collection.d_params.len(), expected_rows);

            let mut expected_rows_by_query = QUERIES.map(|(_, _, row)| row);
            if *drops_final_row {
                expected_rows_by_query[3] = None;
            }

            let mut format_get_calls = 0;
            for ((i_row, j_row, _), expected_row) in QUERIES.into_iter().zip(expected_rows_by_query)
            {
                let input_bytes_before = source.as_bytes().to_vec();
                let i_keys_before = collection.d_i_row.clone();
                let j_keys_before = collection.d_j_row.clone();
                let values_before: Vec<_> = collection
                    .d_params
                    .iter()
                    .map(mmff_emp_herschbach_value_bits)
                    .collect();
                let i_keys_address = collection.d_i_row.as_ptr();
                let j_keys_address = collection.d_j_row.as_ptr();
                let values_address = collection.d_params.as_ptr();
                let query_before = (i_row, j_row);

                let actual = collection.get(i_row, j_row);
                format_get_calls += 1;
                total_get_calls += 1;

                assert_eq!((i_row, j_row), query_before);
                assert_eq!(source.as_bytes(), input_bytes_before.as_slice());
                assert_eq!(collection.d_i_row.as_slice(), i_keys_before.as_slice());
                assert_eq!(collection.d_j_row.as_slice(), j_keys_before.as_slice());
                assert_eq!(
                    collection
                        .d_params
                        .iter()
                        .map(mmff_emp_herschbach_value_bits)
                        .collect::<Vec<_>>(),
                    values_before
                );
                assert_eq!(collection.d_i_row.as_ptr(), i_keys_address);
                assert_eq!(collection.d_j_row.as_ptr(), j_keys_address);
                assert_eq!(collection.d_params.as_ptr(), values_address);

                if let Some(row_index) = expected_row {
                    let actual = actual.expect("fixed Herschbach lookup hit");
                    let fixed_row = &collection.d_params[row_index];
                    assert!(std::ptr::eq(actual, fixed_row));
                    assert_eq!(
                        mmff_emp_herschbach_value_bits(actual),
                        EXPECTED_VALUE_BITS[row_index]
                    );
                } else {
                    assert!(actual.is_none(), "unexpected Herschbach query hit");
                }
            }

            assert_eq!(format_get_calls, 8);
        }
        assert_eq!(total_get_calls, 32);
    }

    #[test]
    fn mmff_emp_herschbach_cast_controls_preserve_full_width_queries_and_signed_zero() {
        let source = "257\t258\t1.25\t-0\t2\n";
        let collection =
            MmffHerschbachLaurieCollection::from_text(source).expect("fixed key-cast row parses");
        assert_eq!(collection.d_i_row, [1]);
        assert_eq!(collection.d_j_row, [2]);

        let mut get_calls = 0;
        let input_bytes_before = source.as_bytes().to_vec();
        let i_keys_before = collection.d_i_row.clone();
        let j_keys_before = collection.d_j_row.clone();
        let values_before: Vec<_> = collection
            .d_params
            .iter()
            .map(mmff_emp_herschbach_value_bits)
            .collect();
        let i_keys_address = collection.d_i_row.as_ptr();
        let j_keys_address = collection.d_j_row.as_ptr();
        let values_address = collection.d_params.as_ptr();
        let actual = collection.get(1, 2);
        get_calls += 1;
        assert_eq!(source.as_bytes(), input_bytes_before.as_slice());
        assert_eq!(collection.d_i_row.as_slice(), i_keys_before.as_slice());
        assert_eq!(collection.d_j_row.as_slice(), j_keys_before.as_slice());
        assert_eq!(
            collection
                .d_params
                .iter()
                .map(mmff_emp_herschbach_value_bits)
                .collect::<Vec<_>>(),
            values_before
        );
        assert_eq!(collection.d_i_row.as_ptr(), i_keys_address);
        assert_eq!(collection.d_j_row.as_ptr(), j_keys_address);
        assert_eq!(collection.d_params.as_ptr(), values_address);
        let actual = actual.expect("narrowed stored keys match the small query");
        assert!(std::ptr::eq(actual, &collection.d_params[0]));
        assert_eq!(
            mmff_emp_herschbach_value_bits(actual),
            [
                0x3ff4_0000_0000_0000,
                0x8000_0000_0000_0000,
                0x4000_0000_0000_0000
            ]
        );

        let input_bytes_before = source.as_bytes().to_vec();
        let i_keys_before = collection.d_i_row.clone();
        let j_keys_before = collection.d_j_row.clone();
        let values_before: Vec<_> = collection
            .d_params
            .iter()
            .map(mmff_emp_herschbach_value_bits)
            .collect();
        let i_keys_address = collection.d_i_row.as_ptr();
        let j_keys_address = collection.d_j_row.as_ptr();
        let values_address = collection.d_params.as_ptr();
        let actual = collection.get(257, 258);
        get_calls += 1;
        assert_eq!(source.as_bytes(), input_bytes_before.as_slice());
        assert_eq!(collection.d_i_row.as_slice(), i_keys_before.as_slice());
        assert_eq!(collection.d_j_row.as_slice(), j_keys_before.as_slice());
        assert_eq!(
            collection
                .d_params
                .iter()
                .map(mmff_emp_herschbach_value_bits)
                .collect::<Vec<_>>(),
            values_before
        );
        assert_eq!(collection.d_i_row.as_ptr(), i_keys_address);
        assert_eq!(collection.d_j_row.as_ptr(), j_keys_address);
        assert_eq!(collection.d_params.as_ptr(), values_address);
        assert!(actual.is_none(), "full-width query keys are not narrowed");
        assert_eq!(get_calls, 2);
    }

    #[test]
    fn mmff_emp_herschbach_constructor_reports_frozen_typed_errors() {
        let valid_cells = ["1", "2", "3.0", "4.0", "5.0"];
        let mut error_constructors = 0;

        for column in 0..valid_cells.len() {
            let mut cells = valid_cells;
            cells[column] = "X";
            let source = format!("*fixed\n{}\n", cells.join("\t"));
            let input_bytes_before = source.as_bytes().to_vec();
            let result = MmffHerschbachLaurieCollection::from_text(&source);
            error_constructors += 1;
            assert_eq!(source.as_bytes(), input_bytes_before.as_slice());

            let error = result.expect_err("each consumed invalid cell is rejected");
            assert_eq!(error.table, MmffParamTable::HerschbachLaurie);
            assert_eq!(error.line, 2);
            assert_eq!(error.column, column);
            let expected_cause = if column < 2 {
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                }
            } else {
                MmffParamParseCause::InvalidFloat {
                    cell: "X".to_owned(),
                }
            };
            assert_eq!(&error.cause, &expected_cause);
            if column == 0 {
                assert_eq!(
                    error.to_string(),
                    "MMFF HerschbachLaurie parse error at physical line 2 token 0: invalid unsigned integer \"X\""
                );
            }
        }

        let missing_final_cell = "*fixed\n1\t2\t3.0\t4.0\n";
        let input_bytes_before = missing_final_cell.as_bytes().to_vec();
        let result = MmffHerschbachLaurieCollection::from_text(missing_final_cell);
        error_constructors += 1;
        assert_eq!(missing_final_cell.as_bytes(), input_bytes_before.as_slice());
        let error = result.expect_err("missing dp_ij is reported");
        assert_eq!(error.table, MmffParamTable::HerschbachLaurie);
        assert_eq!(error.line, 2);
        assert_eq!(error.column, 4);
        assert_eq!(error.cause, MmffParamParseCause::MissingToken);
        assert_eq!(
            error.to_string(),
            "MMFF HerschbachLaurie parse error at physical line 2 token 4: missing token"
        );

        let empty_processed_line = "*fixed\n\n";
        let input_bytes_before = empty_processed_line.as_bytes().to_vec();
        let result = MmffHerschbachLaurieCollection::from_text(empty_processed_line);
        error_constructors += 1;
        assert_eq!(
            empty_processed_line.as_bytes(),
            input_bytes_before.as_slice()
        );
        let error = result.expect_err("empty processed line is source-undefined");
        assert_eq!(error.table, MmffParamTable::HerschbachLaurie);
        assert_eq!(error.line, 2);
        assert_eq!(error.column, 0);
        assert_eq!(error.cause, MmffParamParseCause::EmptyProcessedLine);
        assert_eq!(
            error.to_string(),
            "MMFF HerschbachLaurie parse error at physical line 2 token 0: empty processed line"
        );
        assert_eq!(error_constructors, 7);

        let unused_trailing_token = "*fixed\n1\t2\t3.0\t4.0\t5.0\tX\n";
        let input_bytes_before = unused_trailing_token.as_bytes().to_vec();
        let collection = MmffHerschbachLaurieCollection::from_text(unused_trailing_token)
            .expect("source ignores the unused trailing token");
        assert_eq!(
            unused_trailing_token.as_bytes(),
            input_bytes_before.as_slice()
        );
        assert_eq!(collection.d_i_row, [1]);
        assert_eq!(collection.d_j_row, [2]);
        assert_eq!(
            mmff_emp_herschbach_value_bits(&collection.d_params[0]),
            [
                0x4008_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x4014_0000_0000_0000
            ]
        );
    }

    fn mmff_emp_cov_value_bits(value: &MmffCovRadPauEle) -> [u64; 2] {
        [value.r0.to_bits(), value.chi.to_bits()]
    }

    #[test]
    fn mmff_emp_cov_lookup_matrix_has_exact_32_calls() {
        const SOURCE: &str = concat!(
            "*fixed\n",
            "1\t0.5\t1.25\n",
            "1\t9\t10\n",
            "6\t0.75\t2.5\n",
            "8\t1\t3.5\n",
        );
        const EXPECTED_VALUE_BITS: [[u64; 2]; 4] = [
            [0x3fe0_0000_0000_0000, 0x3ff4_0000_0000_0000],
            [0x4022_0000_0000_0000, 0x4024_0000_0000_0000],
            [0x3fe8_0000_0000_0000, 0x4004_0000_0000_0000],
            [0x3ff0_0000_0000_0000, 0x400c_0000_0000_0000],
        ];
        const QUERIES: [(u32, Option<usize>); 8] = [
            (1, Some(0)),
            (6, Some(2)),
            (8, Some(3)),
            (0, None),
            (3, None),
            (257, None),
            (u32::MAX, None),
            (2, None),
        ];

        let crlf = SOURCE.replace('\n', "\r\n");
        let doubled_tabs = SOURCE.replace('\t', "\t\t");
        let unterminated = SOURCE
            .strip_suffix('\n')
            .expect("fixed table ends in LF")
            .to_owned();
        let formats = [
            (SOURCE.to_owned(), false),
            (crlf, false),
            (doubled_tabs, false),
            (unterminated, true),
        ];

        let mut total_get_calls = 0;
        for (source, drops_final_row) in &formats {
            let collection = MmffCovRadPauEleCollection::from_text(source)
                .expect("fixed CovRadPauEle format parses");
            let expected_rows = if *drops_final_row { 3 } else { 4 };
            assert_eq!(collection.d_atomic_num.len(), expected_rows);
            assert_eq!(collection.d_params.len(), expected_rows);
            let expected_key_rows: &[u8] = if *drops_final_row {
                &[1, 1, 6]
            } else {
                &[1, 1, 6, 8]
            };
            assert_eq!(collection.d_atomic_num.as_slice(), expected_key_rows);
            let parsed_value_bits: Vec<_> = collection
                .d_params
                .iter()
                .map(mmff_emp_cov_value_bits)
                .collect();
            assert_eq!(
                parsed_value_bits.as_slice(),
                &EXPECTED_VALUE_BITS[..expected_rows]
            );

            let mut expected_rows_by_query = QUERIES.map(|(_, row)| row);
            if *drops_final_row {
                expected_rows_by_query[2] = None;
            }

            let mut format_get_calls = 0;
            for ((atomic_num, _), expected_row) in QUERIES.into_iter().zip(expected_rows_by_query) {
                let input_bytes_before = source.as_bytes().to_vec();
                let keys_before = collection.d_atomic_num.clone();
                let values_before: Vec<_> = collection
                    .d_params
                    .iter()
                    .map(mmff_emp_cov_value_bits)
                    .collect();
                let keys_address = collection.d_atomic_num.as_ptr();
                let values_address = collection.d_params.as_ptr();
                let query_before = atomic_num;

                let actual = collection.get(atomic_num);
                format_get_calls += 1;
                total_get_calls += 1;

                assert_eq!(atomic_num, query_before);
                assert_eq!(source.as_bytes(), input_bytes_before.as_slice());
                assert_eq!(collection.d_atomic_num.as_slice(), keys_before.as_slice());
                assert_eq!(
                    collection
                        .d_params
                        .iter()
                        .map(mmff_emp_cov_value_bits)
                        .collect::<Vec<_>>(),
                    values_before
                );
                assert_eq!(collection.d_atomic_num.as_ptr(), keys_address);
                assert_eq!(collection.d_params.as_ptr(), values_address);

                if let Some(row_index) = expected_row {
                    let actual = actual.expect("fixed CovRadPauEle lookup hit");
                    let fixed_row = &collection.d_params[row_index];
                    assert!(std::ptr::eq(actual, fixed_row));
                    assert_eq!(
                        mmff_emp_cov_value_bits(actual),
                        EXPECTED_VALUE_BITS[row_index]
                    );
                } else {
                    assert!(actual.is_none(), "unexpected CovRadPauEle query hit");
                }
            }

            assert_eq!(format_get_calls, 8);
        }
        assert_eq!(total_get_calls, 32);
    }

    #[test]
    fn mmff_emp_cov_cast_controls_preserve_full_width_queries_and_signed_zero() {
        let source = "257\t-0\t1.25\n";
        let collection = MmffCovRadPauEleCollection::from_text(source)
            .expect("fixed atomic-number cast row parses");
        assert_eq!(collection.d_atomic_num, [1]);

        let mut get_calls = 0;
        let input_bytes_before = source.as_bytes().to_vec();
        let keys_before = collection.d_atomic_num.clone();
        let values_before: Vec<_> = collection
            .d_params
            .iter()
            .map(mmff_emp_cov_value_bits)
            .collect();
        let keys_address = collection.d_atomic_num.as_ptr();
        let values_address = collection.d_params.as_ptr();
        let actual = collection.get(1);
        get_calls += 1;
        assert_eq!(source.as_bytes(), input_bytes_before.as_slice());
        assert_eq!(collection.d_atomic_num.as_slice(), keys_before.as_slice());
        assert_eq!(
            collection
                .d_params
                .iter()
                .map(mmff_emp_cov_value_bits)
                .collect::<Vec<_>>(),
            values_before
        );
        assert_eq!(collection.d_atomic_num.as_ptr(), keys_address);
        assert_eq!(collection.d_params.as_ptr(), values_address);
        let actual = actual.expect("narrowed stored key matches the small query");
        assert!(std::ptr::eq(actual, &collection.d_params[0]));
        assert_eq!(
            mmff_emp_cov_value_bits(actual),
            [0x8000_0000_0000_0000, 0x3ff4_0000_0000_0000]
        );

        let input_bytes_before = source.as_bytes().to_vec();
        let keys_before = collection.d_atomic_num.clone();
        let values_before: Vec<_> = collection
            .d_params
            .iter()
            .map(mmff_emp_cov_value_bits)
            .collect();
        let keys_address = collection.d_atomic_num.as_ptr();
        let values_address = collection.d_params.as_ptr();
        let actual = collection.get(257);
        get_calls += 1;
        assert_eq!(source.as_bytes(), input_bytes_before.as_slice());
        assert_eq!(collection.d_atomic_num.as_slice(), keys_before.as_slice());
        assert_eq!(
            collection
                .d_params
                .iter()
                .map(mmff_emp_cov_value_bits)
                .collect::<Vec<_>>(),
            values_before
        );
        assert_eq!(collection.d_atomic_num.as_ptr(), keys_address);
        assert_eq!(collection.d_params.as_ptr(), values_address);
        assert!(actual.is_none(), "full-width query keys are not narrowed");
        assert_eq!(get_calls, 2);
    }

    #[test]
    fn mmff_emp_cov_constructor_reports_frozen_typed_errors() {
        let valid_cells = ["1", "3.0", "4.0"];
        let mut error_constructors = 0;

        for column in 0..valid_cells.len() {
            let mut cells = valid_cells;
            cells[column] = "X";
            let source = format!("*fixed\n{}\n", cells.join("\t"));
            let input_bytes_before = source.as_bytes().to_vec();
            let result = MmffCovRadPauEleCollection::from_text(&source);
            error_constructors += 1;
            assert_eq!(source.as_bytes(), input_bytes_before.as_slice());

            let error = result.expect_err("each consumed invalid cell is rejected");
            assert_eq!(error.table, MmffParamTable::CovRadPauEle);
            assert_eq!(error.line, 2);
            assert_eq!(error.column, column);
            let expected_cause = if column == 0 {
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                }
            } else {
                MmffParamParseCause::InvalidFloat {
                    cell: "X".to_owned(),
                }
            };
            assert_eq!(&error.cause, &expected_cause);
            if column == 0 {
                assert_eq!(
                    error.to_string(),
                    "MMFF CovRadPauEle parse error at physical line 2 token 0: invalid unsigned integer \"X\""
                );
            }
        }

        let missing_final_cell = "*fixed\n1\t3.0\n";
        let input_bytes_before = missing_final_cell.as_bytes().to_vec();
        let result = MmffCovRadPauEleCollection::from_text(missing_final_cell);
        error_constructors += 1;
        assert_eq!(missing_final_cell.as_bytes(), input_bytes_before.as_slice());
        let error = result.expect_err("missing chi is reported");
        assert_eq!(error.table, MmffParamTable::CovRadPauEle);
        assert_eq!(error.line, 2);
        assert_eq!(error.column, 2);
        assert_eq!(error.cause, MmffParamParseCause::MissingToken);
        assert_eq!(
            error.to_string(),
            "MMFF CovRadPauEle parse error at physical line 2 token 2: missing token"
        );

        let empty_processed_line = "*fixed\n\n";
        let input_bytes_before = empty_processed_line.as_bytes().to_vec();
        let result = MmffCovRadPauEleCollection::from_text(empty_processed_line);
        error_constructors += 1;
        assert_eq!(
            empty_processed_line.as_bytes(),
            input_bytes_before.as_slice()
        );
        let error = result.expect_err("empty processed line is source-undefined");
        assert_eq!(error.table, MmffParamTable::CovRadPauEle);
        assert_eq!(error.line, 2);
        assert_eq!(error.column, 0);
        assert_eq!(error.cause, MmffParamParseCause::EmptyProcessedLine);
        assert_eq!(
            error.to_string(),
            "MMFF CovRadPauEle parse error at physical line 2 token 0: empty processed line"
        );
        assert_eq!(error_constructors, 5);

        let unused_trailing_token = "*fixed\n1\t3.0\t4.0\tX\n";
        let input_bytes_before = unused_trailing_token.as_bytes().to_vec();
        let collection = MmffCovRadPauEleCollection::from_text(unused_trailing_token)
            .expect("source ignores the unused trailing token");
        assert_eq!(
            unused_trailing_token.as_bytes(),
            input_bytes_before.as_slice()
        );
        assert_eq!(collection.d_atomic_num, [1]);
        assert_eq!(
            mmff_emp_cov_value_bits(&collection.d_params[0]),
            [0x4008_0000_0000_0000, 0x4010_0000_0000_0000]
        );
    }

    #[test]
    fn mmff_emp_defaults_assets_sentinels_and_128_warm_borrows() {
        const HERSCHBACH_TEXT: &str = include_str!("default_herschbach_laurie.tsv");
        const COV_RAD_PAU_ELE_TEXT: &str = include_str!("default_cov_rad_pau_ele.tsv");

        assert_eq!(HERSCHBACH_TEXT.len(), 546);
        assert_eq!(
            sha256_for_fixed_asset_test(HERSCHBACH_TEXT.as_bytes()),
            [
                0x41, 0xdc, 0xb5, 0xb0, 0xd5, 0x8e, 0x14, 0xe9, 0xaa, 0x1c, 0xf8, 0xd0, 0x43, 0x40,
                0xa0, 0xa4, 0x00, 0x62, 0x73, 0xbd, 0x27, 0x3c, 0x7f, 0xc6, 0xed, 0x7d, 0x67, 0x4d,
                0x07, 0xb4, 0x57, 0xd3,
            ]
        );
        assert_eq!(
            HERSCHBACH_TEXT
                .lines()
                .filter(|line| !line.starts_with('*'))
                .count(),
            25
        );
        assert_eq!(COV_RAD_PAU_ELE_TEXT.len(), 254);
        assert_eq!(
            sha256_for_fixed_asset_test(COV_RAD_PAU_ELE_TEXT.as_bytes()),
            [
                0xad, 0x2f, 0x75, 0x7d, 0xd0, 0x38, 0xbb, 0x1b, 0x51, 0x61, 0xdb, 0x00, 0x3d, 0xeb,
                0x2e, 0xdc, 0x6f, 0xa1, 0xa5, 0x05, 0x59, 0x91, 0xd4, 0xd8, 0x9f, 0x35, 0xe3, 0xfd,
                0xce, 0xa6, 0xfd, 0x48,
            ]
        );
        assert_eq!(
            COV_RAD_PAU_ELE_TEXT
                .lines()
                .filter(|line| !line.starts_with('*'))
                .count(),
            18
        );

        let h_constructions_before =
            DEFAULT_MMFF_HERSCHBACH_LAURIE_CONSTRUCTIONS.load(Ordering::Relaxed);
        let c_constructions_before =
            DEFAULT_MMFF_COV_RAD_PAU_ELE_CONSTRUCTIONS.load(Ordering::Relaxed);
        let custom_herschbach = MmffHerschbachLaurieCollection::from_text("1\t1\t9\t8\t7\n")
            .expect("custom Herschbach rows remain independent");
        let custom_cov = MmffCovRadPauEleCollection::from_text("15\t9\t8\n")
            .expect("custom CovRadPauEle rows remain independent");
        assert_eq!(
            mmff_emp_herschbach_value_bits(&custom_herschbach.d_params[0]),
            [
                0x4022_0000_0000_0000,
                0x4020_0000_0000_0000,
                0x401c_0000_0000_0000
            ]
        );
        assert_eq!(
            mmff_emp_cov_value_bits(&custom_cov.d_params[0]),
            [0x4022_0000_0000_0000, 0x4020_0000_0000_0000]
        );
        assert_eq!(
            DEFAULT_MMFF_HERSCHBACH_LAURIE_CONSTRUCTIONS.load(Ordering::Relaxed),
            h_constructions_before
        );
        assert_eq!(
            DEFAULT_MMFF_COV_RAD_PAU_ELE_CONSTRUCTIONS.load(Ordering::Relaxed),
            c_constructions_before
        );

        let mut baseline_accessor_calls = (0, 0);
        let herschbach =
            default_mmff_herschbach_laurie().expect("fixed Herschbach default asset must parse");
        baseline_accessor_calls.0 += 1;
        let cov =
            default_mmff_cov_rad_pau_ele().expect("fixed CovRadPauEle default asset must parse");
        baseline_accessor_calls.1 += 1;
        assert_eq!(baseline_accessor_calls, (1, 1));

        assert_eq!(herschbach.d_i_row.len(), 25);
        assert_eq!(herschbach.d_j_row.len(), 25);
        assert_eq!(herschbach.d_params.len(), 25);
        assert_eq!(cov.d_atomic_num.len(), 18);
        assert_eq!(cov.d_params.len(), 18);

        assert_eq!((herschbach.d_i_row[0], herschbach.d_j_row[0]), (0, 0));
        assert_eq!(
            mmff_emp_herschbach_value_bits(&herschbach.d_params[0]),
            [
                0x3ff4_28f5_c28f_5c29,
                0x3f99_9999_9999_999a,
                0x3f99_9999_9999_999a
            ]
        );
        assert_eq!((herschbach.d_i_row[12], herschbach.d_j_row[12]), (1, 4));
        assert_eq!(
            mmff_emp_herschbach_value_bits(&herschbach.d_params[12]),
            [
                0x4002_a3d7_0a3d_70a4,
                0x3fe5_c28f_5c28_f5c3,
                0x3ff1_eb85_1eb8_51ec
            ]
        );
        assert_eq!((herschbach.d_i_row[24], herschbach.d_j_row[24]), (4, 5));
        assert_eq!(
            mmff_emp_herschbach_value_bits(&herschbach.d_params[24]),
            [
                0x4006_147a_e147_ae14,
                0x3ff4_0000_0000_0000,
                0x3ff8_28f5_c28f_5c29
            ]
        );
        assert_eq!(
            herschbach
                .get(4, 1)
                .map(|row| std::ptr::eq(row, &herschbach.d_params[12])),
            Some(true)
        );

        assert_eq!(cov.d_atomic_num[0], 1);
        assert_eq!(
            mmff_emp_cov_value_bits(&cov.d_params[0]),
            [0x3fd5_1eb8_51eb_851f, 0x4001_9999_9999_999a]
        );
        assert_eq!(cov.d_atomic_num[9], 15);
        assert_eq!(
            mmff_emp_cov_value_bits(&cov.d_params[9]),
            [0x3ff1_70a3_d70a_3d71, 0x4000_7ae1_47ae_147b]
        );
        assert_eq!(cov.d_atomic_num[17], 53);
        assert_eq!(
            mmff_emp_cov_value_bits(&cov.d_params[17]),
            [0x3ff5_47ae_147a_e148, 0x4001_ae14_7ae1_47ae]
        );
        assert_eq!(
            cov.get(15).map(|row| std::ptr::eq(row, &cov.d_params[9])),
            Some(true)
        );

        let herschbach_i_rows_before = herschbach.d_i_row.clone();
        let herschbach_j_rows_before = herschbach.d_j_row.clone();
        let herschbach_values_before: Vec<_> = herschbach
            .d_params
            .iter()
            .map(mmff_emp_herschbach_value_bits)
            .collect();
        let herschbach_i_rows_address = herschbach.d_i_row.as_ptr() as usize;
        let herschbach_j_rows_address = herschbach.d_j_row.as_ptr() as usize;
        let herschbach_values_address = herschbach.d_params.as_ptr() as usize;
        let herschbach_row_12 = herschbach
            .get(4, 1)
            .expect("fixed Herschbach default warm row exists");
        let herschbach_collection_address =
            herschbach as *const MmffHerschbachLaurieCollection as usize;

        let cov_keys_before = cov.d_atomic_num.clone();
        let cov_values_before: Vec<_> = cov.d_params.iter().map(mmff_emp_cov_value_bits).collect();
        let cov_keys_address = cov.d_atomic_num.as_ptr() as usize;
        let cov_values_address = cov.d_params.as_ptr() as usize;
        let cov_row_9 = cov
            .get(15)
            .expect("fixed CovRadPauEle default warm row exists");
        let cov_collection_address = cov as *const MmffCovRadPauEleCollection as usize;

        let (herschbach_warm_calls, cov_warm_calls) = std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for _ in 0..8 {
                let herschbach_i_rows_before = &herschbach_i_rows_before;
                let herschbach_j_rows_before = &herschbach_j_rows_before;
                let herschbach_values_before = &herschbach_values_before;
                let cov_keys_before = &cov_keys_before;
                let cov_values_before = &cov_values_before;
                handles.push(scope.spawn(move || {
                    let mut herschbach_calls = 0;
                    let mut cov_calls = 0;
                    for _ in 0..16 {
                        let cached_herschbach = default_mmff_herschbach_laurie()
                            .expect("cached Herschbach default remains valid");
                        herschbach_calls += 1;
                        assert_eq!(
                            cached_herschbach as *const MmffHerschbachLaurieCollection as usize,
                            herschbach_collection_address
                        );
                        let cached_herschbach_row = cached_herschbach
                            .get(4, 1)
                            .expect("Herschbach warm sentinel remains present");
                        assert!(std::ptr::eq(cached_herschbach_row, herschbach_row_12));
                        assert_eq!(
                            mmff_emp_herschbach_value_bits(cached_herschbach_row),
                            [
                                0x4002_a3d7_0a3d_70a4,
                                0x3fe5_c28f_5c28_f5c3,
                                0x3ff1_eb85_1eb8_51ec,
                            ]
                        );
                        assert_eq!(
                            cached_herschbach.d_i_row.as_slice(),
                            herschbach_i_rows_before.as_slice()
                        );
                        assert_eq!(
                            cached_herschbach.d_j_row.as_slice(),
                            herschbach_j_rows_before.as_slice()
                        );
                        assert_eq!(
                            cached_herschbach
                                .d_params
                                .iter()
                                .map(mmff_emp_herschbach_value_bits)
                                .collect::<Vec<_>>()
                                .as_slice(),
                            herschbach_values_before.as_slice()
                        );
                        assert_eq!(
                            cached_herschbach.d_i_row.as_ptr() as usize,
                            herschbach_i_rows_address
                        );
                        assert_eq!(
                            cached_herschbach.d_j_row.as_ptr() as usize,
                            herschbach_j_rows_address
                        );
                        assert_eq!(
                            cached_herschbach.d_params.as_ptr() as usize,
                            herschbach_values_address
                        );

                        let cached_cov = default_mmff_cov_rad_pau_ele()
                            .expect("cached CovRadPauEle default remains valid");
                        cov_calls += 1;
                        assert_eq!(
                            cached_cov as *const MmffCovRadPauEleCollection as usize,
                            cov_collection_address
                        );
                        let cached_cov_row = cached_cov
                            .get(15)
                            .expect("CovRadPauEle warm sentinel remains present");
                        assert!(std::ptr::eq(cached_cov_row, cov_row_9));
                        assert_eq!(
                            mmff_emp_cov_value_bits(cached_cov_row),
                            [0x3ff1_70a3_d70a_3d71, 0x4000_7ae1_47ae_147b]
                        );
                        assert_eq!(
                            cached_cov.d_atomic_num.as_slice(),
                            cov_keys_before.as_slice()
                        );
                        assert_eq!(
                            cached_cov
                                .d_params
                                .iter()
                                .map(mmff_emp_cov_value_bits)
                                .collect::<Vec<_>>()
                                .as_slice(),
                            cov_values_before.as_slice()
                        );
                        assert_eq!(cached_cov.d_atomic_num.as_ptr() as usize, cov_keys_address);
                        assert_eq!(cached_cov.d_params.as_ptr() as usize, cov_values_address);
                    }
                    (herschbach_calls, cov_calls)
                }));
            }

            handles
                .into_iter()
                .fold((0, 0), |(herschbach_total, cov_total), handle| {
                    let (herschbach_calls, cov_calls) =
                        handle.join().expect("scoped default accessor worker");
                    (herschbach_total + herschbach_calls, cov_total + cov_calls)
                })
        });
        assert_eq!((herschbach_warm_calls, cov_warm_calls), (128, 128));
        assert_eq!(
            DEFAULT_MMFF_HERSCHBACH_LAURIE_CONSTRUCTIONS.load(Ordering::Relaxed),
            1
        );
        assert_eq!(
            DEFAULT_MMFF_COV_RAD_PAU_ELE_CONSTRUCTIONS.load(Ordering::Relaxed),
            1
        );
    }

    #[test]
    fn mmff_sb_defaults_assets_sentinels_and_128_warm_borrows_each() {
        const STBN_TEXT: &str = include_str!("default_stbn.tsv");
        const DFSB_TEXT: &str = include_str!("default_dfsb.tsv");

        assert_eq!(STBN_TEXT.len(), 7621);
        assert_eq!(
            sha256_for_fixed_asset_test(STBN_TEXT.as_bytes()),
            [
                0x3d, 0x54, 0xfe, 0x7c, 0x56, 0x10, 0x8c, 0xcb, 0xbe, 0x7c, 0x06, 0xdc, 0x93, 0xad,
                0xb6, 0xb9, 0x11, 0xac, 0xb4, 0x49, 0x54, 0x16, 0xf7, 0xff, 0x4f, 0xaa, 0x4b, 0xc4,
                0xa7, 0x3f, 0xc5, 0xb4,
            ]
        );
        assert_eq!(
            STBN_TEXT
                .lines()
                .filter(|line| !line.starts_with('*'))
                .count(),
            282
        );
        assert_eq!(DFSB_TEXT.len(), 711);
        assert_eq!(
            sha256_for_fixed_asset_test(DFSB_TEXT.as_bytes()),
            [
                0xe4, 0x15, 0x0c, 0x47, 0x19, 0xc1, 0x5b, 0x56, 0x9e, 0x77, 0x93, 0x06, 0xc7, 0xd4,
                0x64, 0x7f, 0x9f, 0x14, 0xf9, 0x49, 0xc4, 0x71, 0x2a, 0xf3, 0x9e, 0xa6, 0x7a, 0x1e,
                0x45, 0x8d, 0xdf, 0x55,
            ]
        );
        assert_eq!(
            DFSB_TEXT
                .lines()
                .filter(|line| !line.starts_with('*'))
                .count(),
            30
        );

        let stbn_constructions_before = DEFAULT_MMFF_STBN_CONSTRUCTIONS.load(Ordering::Relaxed);
        let dfsb_constructions_before = DEFAULT_MMFF_DFSB_CONSTRUCTIONS.load(Ordering::Relaxed);
        let custom_stbn = MmffStbnCollection::from_text("1\t1\t2\t3\t1.25\t2.5\n")
            .expect("custom Stbn rows remain independent");
        let custom_dfsb = MmffDfsbCollection::from_text("1\t2\t3\t1.25\t2.5\n")
            .expect("custom Dfsb rows remain independent");
        assert_eq!(
            mmff_sb_stbn_value_bits(&custom_stbn.d_params[0]),
            [0x3ff4_0000_0000_0000, 0x4004_0000_0000_0000]
        );
        assert_eq!(
            mmff_sb_stbn_value_bits(
                custom_dfsb
                    .d_params
                    .get(&1)
                    .and_then(|row2| row2.get(&2))
                    .and_then(|row3| row3.get(&3))
                    .expect("custom Dfsb full-width key is stored")
            ),
            [0x3ff4_0000_0000_0000, 0x4004_0000_0000_0000]
        );
        assert_eq!(
            DEFAULT_MMFF_STBN_CONSTRUCTIONS.load(Ordering::Relaxed),
            stbn_constructions_before
        );
        assert_eq!(
            DEFAULT_MMFF_DFSB_CONSTRUCTIONS.load(Ordering::Relaxed),
            dfsb_constructions_before
        );

        let stbn = default_mmff_stbn().expect("fixed Stbn default asset must parse");
        let dfsb = default_mmff_dfsb().expect("fixed Dfsb default asset must parse");
        assert_eq!(DEFAULT_MMFF_STBN_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(DEFAULT_MMFF_DFSB_CONSTRUCTIONS.load(Ordering::Relaxed), 1);

        assert_eq!(stbn.d_i_atom_type.len(), 282);
        assert_eq!(stbn.d_j_atom_type.len(), 282);
        assert_eq!(stbn.d_k_atom_type.len(), 282);
        assert_eq!(stbn.d_stretch_bend_type.len(), 282);
        assert_eq!(stbn.d_params.len(), 282);
        assert_eq!(
            mmff_sb_stbn_storage_addresses(stbn)[..4],
            [
                stbn.d_i_atom_type.as_ptr() as usize,
                stbn.d_j_atom_type.as_ptr() as usize,
                stbn.d_k_atom_type.as_ptr() as usize,
                stbn.d_stretch_bend_type.as_ptr() as usize,
            ]
        );

        assert_eq!(
            (
                stbn.d_stretch_bend_type[0],
                stbn.d_i_atom_type[0],
                stbn.d_j_atom_type[0],
                stbn.d_k_atom_type[0],
            ),
            (0, 1, 1, 1)
        );
        assert_eq!(
            mmff_sb_stbn_value_bits(&stbn.d_params[0]),
            [0x3fca_5e35_3f7c_ed91, 0x3fca_5e35_3f7c_ed91]
        );
        assert_eq!(
            (
                stbn.d_stretch_bend_type[141],
                stbn.d_i_atom_type[141],
                stbn.d_j_atom_type[141],
                stbn.d_k_atom_type[141],
            ),
            (0, 1, 18, 1)
        );
        assert_eq!(
            mmff_sb_stbn_value_bits(&stbn.d_params[141]),
            [0x3f97_8d4f_df3b_645a, 0x3f97_8d4f_df3b_645a]
        );
        assert_eq!(
            (
                stbn.d_stretch_bend_type[281],
                stbn.d_i_atom_type[281],
                stbn.d_j_atom_type[281],
                stbn.d_k_atom_type[281],
            ),
            (0, 78, 81, 80)
        );
        assert_eq!(
            mmff_sb_stbn_value_bits(&stbn.d_params[281]),
            [0x3fd7_6c8b_4395_8106, 0x3fda_d0e5_6041_8937]
        );
        let (stbn_swap, stbn_sentinel) = stbn.get(0, 0, 0, 78, 81, 80);
        assert!(!stbn_swap);
        let stbn_sentinel = stbn_sentinel.expect("fixed Stbn warm sentinel exists");
        assert!(std::ptr::eq(stbn_sentinel, &stbn.d_params[281]));
        assert_eq!(
            mmff_sb_stbn_value_bits(stbn_sentinel),
            [0x3fd7_6c8b_4395_8106, 0x3fda_d0e5_6041_8937]
        );

        let dfsb_snapshot_before = mmff_sb_dfsb_snapshot(dfsb);
        assert_eq!(dfsb_snapshot_before.len(), 30);
        assert_eq!(
            (0, 1, 0),
            (
                dfsb_snapshot_before[0].0,
                dfsb_snapshot_before[0].1,
                dfsb_snapshot_before[0].2
            )
        );
        assert_eq!(
            dfsb_snapshot_before[0].3,
            [0x3fc3_3333_3333_3333, 0x3fc3_3333_3333_3333]
        );
        assert_eq!(
            (2, 1, 3),
            (
                dfsb_snapshot_before[19].0,
                dfsb_snapshot_before[19].1,
                dfsb_snapshot_before[19].2
            )
        );
        assert_eq!(
            dfsb_snapshot_before[19].3,
            [0x3fe0_0000_0000_0000, 0x3fe0_0000_0000_0000]
        );
        assert_eq!(
            (4, 2, 4),
            (
                dfsb_snapshot_before[29].0,
                dfsb_snapshot_before[29].1,
                dfsb_snapshot_before[29].2
            )
        );
        assert_eq!(
            dfsb_snapshot_before[29].3,
            [0x3fd0_0000_0000_0000, 0x3fd0_0000_0000_0000]
        );
        let dfsb_row_0 = dfsb
            .d_params
            .get(&0)
            .and_then(|row2| row2.get(&1))
            .and_then(|row3| row3.get(&0))
            .expect("fixed Dfsb row 0 exists");
        let dfsb_row_2_1_3 = dfsb
            .d_params
            .get(&2)
            .and_then(|row2| row2.get(&1))
            .and_then(|row3| row3.get(&3))
            .expect("fixed Dfsb row (2, 1, 3) exists");
        let dfsb_row_29 = dfsb
            .d_params
            .get(&4)
            .and_then(|row2| row2.get(&2))
            .and_then(|row3| row3.get(&4))
            .expect("fixed Dfsb row 29 exists");
        assert_eq!(
            dfsb_snapshot_before[0].4,
            dfsb_row_0 as *const MmffStbn as usize
        );
        assert_eq!(
            dfsb_snapshot_before[19].4,
            dfsb_row_2_1_3 as *const MmffStbn as usize
        );
        assert_eq!(
            dfsb_snapshot_before[29].4,
            dfsb_row_29 as *const MmffStbn as usize
        );
        let (dfsb_swap_0, dfsb_hit_0) = dfsb.get(0, 1, 0);
        assert!(!dfsb_swap_0);
        assert!(std::ptr::eq(
            dfsb_hit_0.expect("fixed Dfsb row 0 lookup"),
            dfsb_row_0
        ));
        assert_eq!(
            mmff_sb_stbn_value_bits(dfsb_row_0),
            [0x3fc3_3333_3333_3333, 0x3fc3_3333_3333_3333]
        );
        let (dfsb_swap_2_1_3, dfsb_hit_2_1_3) = dfsb.get(3, 1, 2);
        assert!(dfsb_swap_2_1_3);
        assert!(std::ptr::eq(
            dfsb_hit_2_1_3.expect("fixed Dfsb (2, 1, 3) reverse lookup"),
            dfsb_row_2_1_3
        ));
        assert_eq!(
            mmff_sb_stbn_value_bits(dfsb_row_2_1_3),
            [0x3fe0_0000_0000_0000, 0x3fe0_0000_0000_0000]
        );
        let (dfsb_swap_29, dfsb_hit_29) = dfsb.get(4, 2, 4);
        assert!(!dfsb_swap_29);
        assert!(std::ptr::eq(
            dfsb_hit_29.expect("fixed Dfsb row 29 lookup"),
            dfsb_row_29
        ));
        assert_eq!(
            mmff_sb_stbn_value_bits(dfsb_row_29),
            [0x3fd0_0000_0000_0000, 0x3fd0_0000_0000_0000]
        );

        let stbn_keys_before = [
            stbn.d_i_atom_type.clone(),
            stbn.d_j_atom_type.clone(),
            stbn.d_k_atom_type.clone(),
            stbn.d_stretch_bend_type.clone(),
        ];
        let stbn_values_before: Vec<_> =
            stbn.d_params.iter().map(mmff_sb_stbn_value_bits).collect();
        let stbn_addresses_before = mmff_sb_stbn_storage_addresses(stbn);
        let stbn_collection_address = stbn as *const MmffStbnCollection as usize;
        let stbn_input_bytes_before = STBN_TEXT.as_bytes().to_vec();
        let dfsb_map_address_before = mmff_sb_dfsb_map_address(dfsb);
        let dfsb_collection_address = dfsb as *const MmffDfsbCollection as usize;
        let dfsb_input_bytes_before = DFSB_TEXT.as_bytes().to_vec();

        let (stbn_warm_calls, dfsb_warm_calls) = std::thread::scope(|scope| {
            let stbn_keys_before = &stbn_keys_before;
            let stbn_values_before = &stbn_values_before;
            let stbn_input_bytes_before = &stbn_input_bytes_before;
            let dfsb_snapshot_before = &dfsb_snapshot_before;
            let dfsb_input_bytes_before = &dfsb_input_bytes_before;
            let mut handles = Vec::with_capacity(8);
            for _ in 0..8 {
                let stbn_keys_before = stbn_keys_before;
                let stbn_values_before = stbn_values_before;
                let stbn_input_bytes_before = stbn_input_bytes_before;
                let dfsb_snapshot_before = dfsb_snapshot_before;
                let dfsb_input_bytes_before = dfsb_input_bytes_before;
                handles.push(scope.spawn(move || {
                    let mut stbn_calls = 0;
                    let mut dfsb_calls = 0;
                    for _ in 0..16 {
                        let cached_stbn =
                            default_mmff_stbn().expect("cached Stbn default remains valid");
                        stbn_calls += 1;
                        assert_eq!(
                            cached_stbn as *const MmffStbnCollection as usize,
                            stbn_collection_address
                        );
                        let (swap, row) = cached_stbn.get(0, 0, 0, 78, 81, 80);
                        assert!(!swap);
                        let row = row.expect("Stbn warm sentinel remains present");
                        assert!(std::ptr::eq(row, stbn_sentinel));
                        assert_eq!(
                            mmff_sb_stbn_value_bits(row),
                            [0x3fd7_6c8b_4395_8106, 0x3fda_d0e5_6041_8937]
                        );
                        assert_mmff_sb_stbn_state_unchanged(
                            STBN_TEXT,
                            stbn_input_bytes_before,
                            cached_stbn,
                            stbn_keys_before,
                            stbn_values_before,
                            stbn_addresses_before,
                        );

                        let cached_dfsb =
                            default_mmff_dfsb().expect("cached Dfsb default remains valid");
                        dfsb_calls += 1;
                        assert_eq!(
                            cached_dfsb as *const MmffDfsbCollection as usize,
                            dfsb_collection_address
                        );
                        let (swap, row) = cached_dfsb.get(3, 1, 2);
                        assert!(swap);
                        let row = row.expect("Dfsb warm reverse sentinel remains present");
                        assert!(std::ptr::eq(row, dfsb_row_2_1_3));
                        assert_eq!(
                            mmff_sb_stbn_value_bits(row),
                            [0x3fe0_0000_0000_0000, 0x3fe0_0000_0000_0000]
                        );
                        assert_mmff_sb_dfsb_state_unchanged(
                            DFSB_TEXT,
                            dfsb_input_bytes_before,
                            cached_dfsb,
                            dfsb_map_address_before,
                            dfsb_snapshot_before,
                        );
                    }
                    (stbn_calls, dfsb_calls)
                }));
            }

            handles
                .into_iter()
                .fold((0, 0), |(stbn_total, dfsb_total), handle| {
                    let (stbn_calls, dfsb_calls) =
                        handle.join().expect("scoped MMFF-SB default worker");
                    (stbn_total + stbn_calls, dfsb_total + dfsb_calls)
                })
        });
        assert_eq!((stbn_warm_calls, dfsb_warm_calls), (128, 128));
        assert_eq!(DEFAULT_MMFF_STBN_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(DEFAULT_MMFF_DFSB_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
    }

    #[derive(Clone, Copy, Debug)]
    enum ExpectedAngleLookup {
        Hit { row: usize, value_bits: [u64; 2] },
        Miss,
        MissingDefinition { atom_type: u32 },
    }

    struct MmffAngleLookupSnapshot {
        angle_source_bytes: Vec<u8>,
        angle_source_address: usize,
        def_source_bytes: Vec<u8>,
        def_source_address: usize,
        def_eq_levels: Vec<[u8; 4]>,
        def_rows_address: usize,
        angle_key_rows: [Vec<u8>; 4],
        angle_value_bits: Vec<[u64; 2]>,
        angle_storage_addresses: [usize; 5],
    }

    impl MmffAngleLookupSnapshot {
        fn capture(
            angle_source: &str,
            def_source: &str,
            definitions: &MmffDefCollection,
            collection: &MmffAngleCollection,
        ) -> Self {
            Self {
                angle_source_bytes: angle_source.as_bytes().to_vec(),
                angle_source_address: angle_source.as_ptr() as usize,
                def_source_bytes: def_source.as_bytes().to_vec(),
                def_source_address: def_source.as_ptr() as usize,
                def_eq_levels: definitions
                    .d_params
                    .iter()
                    .map(|row| row.eq_level)
                    .collect(),
                def_rows_address: definitions.d_params.as_ptr() as usize,
                angle_key_rows: [
                    collection.d_i_atom_type.clone(),
                    collection.d_j_atom_type.clone(),
                    collection.d_k_atom_type.clone(),
                    collection.d_angle_type.clone(),
                ],
                angle_value_bits: collection
                    .d_params
                    .iter()
                    .map(|row| [row.ka.to_bits(), row.theta0.to_bits()])
                    .collect(),
                angle_storage_addresses: [
                    collection.d_i_atom_type.as_ptr() as usize,
                    collection.d_j_atom_type.as_ptr() as usize,
                    collection.d_k_atom_type.as_ptr() as usize,
                    collection.d_angle_type.as_ptr() as usize,
                    collection.d_params.as_ptr() as usize,
                ],
            }
        }

        fn assert_unchanged(
            &self,
            angle_source: &str,
            def_source: &str,
            definitions: &MmffDefCollection,
            collection: &MmffAngleCollection,
            context: &str,
        ) {
            assert_eq!(
                angle_source.as_bytes(),
                self.angle_source_bytes,
                "{context}"
            );
            assert_eq!(
                angle_source.as_ptr() as usize,
                self.angle_source_address,
                "{context}"
            );
            assert_eq!(def_source.as_bytes(), self.def_source_bytes, "{context}");
            assert_eq!(
                def_source.as_ptr() as usize,
                self.def_source_address,
                "{context}"
            );
            assert_eq!(
                definitions
                    .d_params
                    .iter()
                    .map(|row| row.eq_level)
                    .collect::<Vec<_>>(),
                self.def_eq_levels,
                "{context}"
            );
            assert_eq!(
                definitions.d_params.as_ptr() as usize,
                self.def_rows_address,
                "{context}"
            );
            assert_eq!(
                collection.d_i_atom_type.as_slice(),
                self.angle_key_rows[0].as_slice(),
                "{context}"
            );
            assert_eq!(
                collection.d_j_atom_type.as_slice(),
                self.angle_key_rows[1].as_slice(),
                "{context}"
            );
            assert_eq!(
                collection.d_k_atom_type.as_slice(),
                self.angle_key_rows[2].as_slice(),
                "{context}"
            );
            assert_eq!(
                collection.d_angle_type.as_slice(),
                self.angle_key_rows[3].as_slice(),
                "{context}"
            );
            assert_eq!(
                collection
                    .d_params
                    .iter()
                    .map(|row| [row.ka.to_bits(), row.theta0.to_bits()])
                    .collect::<Vec<_>>(),
                self.angle_value_bits,
                "{context}"
            );
            assert_eq!(
                [
                    collection.d_i_atom_type.as_ptr() as usize,
                    collection.d_j_atom_type.as_ptr() as usize,
                    collection.d_k_atom_type.as_ptr() as usize,
                    collection.d_angle_type.as_ptr() as usize,
                    collection.d_params.as_ptr() as usize,
                ],
                self.angle_storage_addresses,
                "{context}"
            );
        }
    }

    fn assert_mmff_angle_lookup_result(
        result: Result<Option<&MmffAngle>, MmffAngleLookupError>,
        expected: ExpectedAngleLookup,
        collection: &MmffAngleCollection,
        context: &str,
    ) {
        match (result, expected) {
            (Ok(Some(actual)), ExpectedAngleLookup::Hit { row, value_bits }) => {
                assert!(std::ptr::eq(actual, &collection.d_params[row]), "{context}");
                assert_eq!(
                    [actual.ka.to_bits(), actual.theta0.to_bits()],
                    value_bits,
                    "{context}"
                );
            }
            (Ok(None), ExpectedAngleLookup::Miss) => {}
            (Err(error), ExpectedAngleLookup::MissingDefinition { atom_type }) => {
                assert_eq!(
                    error,
                    MmffAngleLookupError::MissingDefinition { atom_type },
                    "{context}"
                );
                assert_eq!(
                    error.to_string(),
                    format!("missing MMFF definition for atom type {atom_type}"),
                    "{context}"
                );
                assert!(std::error::Error::source(&error).is_none(), "{context}");
            }
            (actual, expected) => panic!("{context}: actual {actual:?}, expected {expected:?}"),
        }
    }

    #[test]
    fn mmff_angle_lookup_matrix_48_queries_preserves_source_order_and_state() {
        const DEF_SOURCE: &str = concat!(
            "*fixed\n",
            "X\t1\t10\t11\t12\t13\n",
            "X\t2\t20\t21\t22\t23\n",
            "X\t3\t30\t31\t32\t33\n",
        );
        const SOURCE: &str = concat!(
            "*fixed\n",
            "0\t10\t7\t20\t1.25\t90\n",
            "1\t10\t7\t20\t2.5\t100\n",
            "1\t10\t7\t20\t9\t110\n",
            "2\t11\t7\t21\t0.5\t120\n",
            "3\t12\t7\t22\t0.75\t130\n",
            "4\t13\t7\t23\t0.125\t140\n",
            "0\t10\t8\t20\t3.75\t150\n",
        );
        const VALUE_BITS: [[u64; 2]; 7] = [
            [0x3ff4_0000_0000_0000, 0x4056_8000_0000_0000],
            [0x4004_0000_0000_0000, 0x4059_0000_0000_0000],
            [0x4022_0000_0000_0000, 0x405b_8000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x405e_0000_0000_0000],
            [0x3fe8_0000_0000_0000, 0x4060_4000_0000_0000],
            [0x3fc0_0000_0000_0000, 0x4061_8000_0000_0000],
            [0x400e_0000_0000_0000, 0x4062_c000_0000_0000],
        ];
        const QUERIES: [(u32, u32, u32, u32, Option<usize>); 12] = [
            (0, 1, 7, 2, Some(0)),
            (0, 2, 7, 1, Some(0)),
            (1, 1, 7, 2, Some(1)),
            (2, 1, 7, 2, Some(3)),
            (3, 1, 7, 2, Some(4)),
            (4, 1, 7, 2, Some(5)),
            (5, 1, 7, 2, None),
            (0, 1, 8, 2, Some(6)),
            (0, 0, 9, 99, None),
            (0, 1, 263, 2, None),
            (256, 1, 7, 2, None),
            (0, 3, 7, 2, None),
        ];

        let source_cases = [
            ("LF", SOURCE.to_owned()),
            ("CRLF", SOURCE.replace('\n', "\r\n")),
            ("doubled TAB", SOURCE.replace('\t', "\t\t")),
            (
                "unterminated final row",
                SOURCE
                    .strip_suffix('\n')
                    .expect("the fixture has a final LF")
                    .to_owned(),
            ),
        ];
        let def_source = DEF_SOURCE.to_owned();
        let definitions =
            MmffDefCollection::from_text(&def_source).expect("the frozen Def fixture is valid");
        assert_eq!(definitions.d_params[0].eq_level, [10, 11, 12, 13]);
        assert_eq!(definitions.d_params[1].eq_level, [20, 21, 22, 23]);
        assert_eq!(definitions.d_params[2].eq_level, [30, 31, 32, 33]);

        let mut get_calls = 0;
        for (format_index, (label, source)) in source_cases.into_iter().enumerate() {
            let collection = MmffAngleCollection::from_text(&source)
                .unwrap_or_else(|error| panic!("{label}: {error}"));
            for &(angle_type, i_atom_type, j_atom_type, k_atom_type, expected_row) in &QUERIES {
                assert!(get_calls < 48, "unexpected extra Angle lookup");
                let expected = match expected_row {
                    Some(6) if format_index == 3 => ExpectedAngleLookup::Miss,
                    Some(row) => ExpectedAngleLookup::Hit {
                        row,
                        value_bits: VALUE_BITS[row],
                    },
                    None => ExpectedAngleLookup::Miss,
                };
                let before = MmffAngleLookupSnapshot::capture(
                    &source,
                    &def_source,
                    &definitions,
                    &collection,
                );
                let result = collection.get(
                    &definitions,
                    angle_type,
                    i_atom_type,
                    j_atom_type,
                    k_atom_type,
                );
                before.assert_unchanged(&source, &def_source, &definitions, &collection, label);
                get_calls += 1;
                assert_mmff_angle_lookup_result(result, expected, &collection, label);
            }
        }
        assert_eq!(get_calls, 48);
    }

    #[test]
    fn mmff_angle_lookup_priority_8_calls_preserve_first_available_stage() {
        const DEF_SOURCE: &str = concat!(
            "*fixed\n",
            "X\t1\t10\t11\t12\t13\n",
            "X\t2\t20\t21\t22\t23\n",
            "X\t3\t30\t31\t32\t33\n",
        );
        const TABLES: [(&str, usize, [u64; 2], usize); 4] = [
            (
                "*fixed\n0\t10\t7\t20\t1.25\t90\n0\t11\t7\t21\t0.5\t120\n0\t12\t7\t22\t0.75\t130\n0\t13\t7\t23\t0.125\t140\n",
                0,
                [0x3ff4_0000_0000_0000, 0x4056_8000_0000_0000],
                4,
            ),
            (
                "*fixed\n0\t11\t7\t21\t0.5\t120\n0\t12\t7\t22\t0.75\t130\n0\t13\t7\t23\t0.125\t140\n",
                3,
                [0x3fe0_0000_0000_0000, 0x405e_0000_0000_0000],
                3,
            ),
            (
                "*fixed\n0\t12\t7\t22\t0.75\t130\n0\t13\t7\t23\t0.125\t140\n",
                4,
                [0x3fe8_0000_0000_0000, 0x4060_4000_0000_0000],
                2,
            ),
            (
                "*fixed\n0\t13\t7\t23\t0.125\t140\n",
                5,
                [0x3fc0_0000_0000_0000, 0x4061_8000_0000_0000],
                1,
            ),
        ];
        const EXPECTED_SOURCE_ROW_INDICES: [usize; 8] = [0, 0, 3, 3, 4, 4, 5, 5];

        let def_source = DEF_SOURCE.to_owned();
        let definitions =
            MmffDefCollection::from_text(&def_source).expect("the frozen Def fixture is valid");
        let mut source_row_indices = Vec::with_capacity(8);
        let mut get_calls = 0;
        for &(table_source, expected_source_row, expected_bits, expected_rows) in &TABLES {
            let source = table_source.to_owned();
            let collection =
                MmffAngleCollection::from_text(&source).expect("the fixed priority table is valid");
            assert_eq!(collection.d_params.len(), expected_rows);
            for &(i_atom_type, k_atom_type) in &[(1_u32, 2_u32), (2, 1)] {
                assert!(get_calls < 8, "unexpected extra stage-priority lookup");
                let before = MmffAngleLookupSnapshot::capture(
                    &source,
                    &def_source,
                    &definitions,
                    &collection,
                );
                let result = collection.get(&definitions, 0, i_atom_type, 7, k_atom_type);
                before.assert_unchanged(
                    &source,
                    &def_source,
                    &definitions,
                    &collection,
                    "stage priority",
                );
                get_calls += 1;
                assert_mmff_angle_lookup_result(
                    result,
                    ExpectedAngleLookup::Hit {
                        row: 0,
                        value_bits: expected_bits,
                    },
                    &collection,
                    "stage priority",
                );
                source_row_indices.push(expected_source_row);
            }
        }
        assert_eq!(get_calls, 8);
        assert_eq!(source_row_indices, EXPECTED_SOURCE_ROW_INDICES);
    }

    #[test]
    fn mmff_angle_lookup_cast_width_5_calls_preserve_full_width_queries() {
        const DEF_SOURCE: &str = concat!(
            "*fixed\n",
            "X\t1\t10\t11\t12\t13\n",
            "X\t2\t20\t21\t22\t23\n",
            "X\t3\t30\t31\t32\t33\n",
        );
        const SOURCE: &str = "256\t266\t263\t276\t-0\t180\n";
        let source = SOURCE.to_owned();
        let def_source = DEF_SOURCE.to_owned();
        let definitions =
            MmffDefCollection::from_text(&def_source).expect("the frozen Def fixture is valid");
        let collection =
            MmffAngleCollection::from_text(&source).expect("the frozen narrowed-key row is valid");
        assert_eq!(collection.d_angle_type.as_slice(), &[0]);
        assert_eq!(collection.d_i_atom_type.as_slice(), &[10]);
        assert_eq!(collection.d_j_atom_type.as_slice(), &[7]);
        assert_eq!(collection.d_k_atom_type.as_slice(), &[20]);
        assert_eq!(
            [
                collection.d_params[0].ka.to_bits(),
                collection.d_params[0].theta0.to_bits(),
            ],
            [0x8000_0000_0000_0000, 0x4066_8000_0000_0000]
        );

        const QUERIES: [(u32, u32, u32, u32, ExpectedAngleLookup); 5] = [
            (
                0,
                1,
                7,
                2,
                ExpectedAngleLookup::Hit {
                    row: 0,
                    value_bits: [0x8000_0000_0000_0000, 0x4066_8000_0000_0000],
                },
            ),
            (256, 1, 7, 2, ExpectedAngleLookup::Miss),
            (0, 1, 263, 2, ExpectedAngleLookup::Miss),
            (
                0,
                257,
                7,
                2,
                ExpectedAngleLookup::MissingDefinition { atom_type: 257 },
            ),
            (
                0,
                1,
                7,
                258,
                ExpectedAngleLookup::MissingDefinition { atom_type: 258 },
            ),
        ];
        let mut get_calls = 0;
        for &(angle_type, i_atom_type, j_atom_type, k_atom_type, expected) in &QUERIES {
            assert!(get_calls < 5, "unexpected extra cast-width lookup");
            let before =
                MmffAngleLookupSnapshot::capture(&source, &def_source, &definitions, &collection);
            let result = collection.get(
                &definitions,
                angle_type,
                i_atom_type,
                j_atom_type,
                k_atom_type,
            );
            before.assert_unchanged(
                &source,
                &def_source,
                &definitions,
                &collection,
                "cast and full-width controls",
            );
            get_calls += 1;
            assert_mmff_angle_lookup_result(
                result,
                expected,
                &collection,
                "cast and full-width controls",
            );
        }
        assert_eq!(get_calls, 5);
    }

    #[test]
    fn mmff_angle_lookup_missing_def_5_calls_keep_typed_error_precedence() {
        const DEF_SOURCE: &str = concat!(
            "*fixed\n",
            "X\t1\t10\t11\t12\t13\n",
            "X\t2\t20\t21\t22\t23\n",
            "X\t3\t30\t31\t32\t33\n",
        );
        const SOURCE: &str = concat!(
            "*fixed\n",
            "0\t10\t7\t20\t1.25\t90\n",
            "1\t10\t7\t20\t2.5\t100\n",
            "1\t10\t7\t20\t9\t110\n",
            "2\t11\t7\t21\t0.5\t120\n",
            "3\t12\t7\t22\t0.75\t130\n",
            "4\t13\t7\t23\t0.125\t140\n",
            "0\t10\t8\t20\t3.75\t150\n",
        );
        const QUERIES: [(u32, u32, u32, u32, u32); 5] = [
            (0, 0, 7, 2, 0),
            (0, 99, 7, 2, 99),
            (0, 1, 7, 0, 0),
            (0, 1, 7, 99, 99),
            (0, 0, 7, 0, 0),
        ];

        let source = SOURCE.to_owned();
        let def_source = DEF_SOURCE.to_owned();
        let definitions =
            MmffDefCollection::from_text(&def_source).expect("the frozen Def fixture is valid");
        let collection =
            MmffAngleCollection::from_text(&source).expect("the frozen Angle fixture is valid");
        let mut get_calls = 0;
        for &(angle_type, i_atom_type, j_atom_type, k_atom_type, expected_atom_type) in &QUERIES {
            assert!(get_calls < 5, "unexpected extra missing-Def lookup");
            let before =
                MmffAngleLookupSnapshot::capture(&source, &def_source, &definitions, &collection);
            let result = collection.get(
                &definitions,
                angle_type,
                i_atom_type,
                j_atom_type,
                k_atom_type,
            );
            before.assert_unchanged(
                &source,
                &def_source,
                &definitions,
                &collection,
                "missing Def safety and precedence",
            );
            get_calls += 1;
            assert_mmff_angle_lookup_result(
                result,
                ExpectedAngleLookup::MissingDefinition {
                    atom_type: expected_atom_type,
                },
                &collection,
                "missing Def safety and precedence",
            );
        }
        assert_eq!(get_calls, 5);
    }

    #[test]
    fn mmff_angle_constructor_formats_keep_literal_key_and_value_arrays() {
        const SOURCE: &str = "*fixed\n0\t10\t7\t20\t1.25\t90\n1\t10\t7\t20\t2.5\t100\n1\t10\t7\t20\t9\t110\n2\t11\t7\t21\t0.5\t120\n3\t12\t7\t22\t0.75\t130\n4\t13\t7\t23\t0.125\t140\n0\t10\t8\t20\t3.75\t150\n";
        let cases = [
            ("LF", SOURCE.to_owned(), 7),
            ("CRLF", SOURCE.replace('\n', "\r\n"), 7),
            ("doubled TAB", SOURCE.replace('\t', "\t\t"), 7),
            (
                "unterminated final row",
                SOURCE
                    .strip_suffix('\n')
                    .expect("fixture has final LF")
                    .to_owned(),
                6,
            ),
        ];
        let expected_angle = [0_u8, 1, 1, 2, 3, 4, 0];
        let expected_i = [10_u8, 10, 10, 11, 12, 13, 10];
        let expected_j = [7_u8, 7, 7, 7, 7, 7, 8];
        let expected_k = [20_u8, 20, 20, 21, 22, 23, 20];
        let expected_values = [
            [0x3ff4_0000_0000_0000, 0x4056_8000_0000_0000],
            [0x4004_0000_0000_0000, 0x4059_0000_0000_0000],
            [0x4022_0000_0000_0000, 0x405b_8000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x405e_0000_0000_0000],
            [0x3fe8_0000_0000_0000, 0x4060_4000_0000_0000],
            [0x3fc0_0000_0000_0000, 0x4061_8000_0000_0000],
            [0x400e_0000_0000_0000, 0x4062_c000_0000_0000],
        ];

        let mut constructor_calls = 0;
        for (label, source, expected_rows) in cases {
            let input_before = source.as_bytes().to_vec();
            let collection = MmffAngleCollection::from_text(&source)
                .unwrap_or_else(|error| panic!("{label}: {error}"));
            constructor_calls += 1;

            assert_eq!(source.as_bytes(), input_before.as_slice(), "{label}");
            assert_eq!(collection.d_params.len(), expected_rows, "{label}");
            assert_eq!(collection.d_i_atom_type.len(), expected_rows, "{label}");
            assert_eq!(collection.d_j_atom_type.len(), expected_rows, "{label}");
            assert_eq!(collection.d_k_atom_type.len(), expected_rows, "{label}");
            assert_eq!(collection.d_angle_type.len(), expected_rows, "{label}");
            assert_eq!(
                collection.d_i_atom_type.as_slice(),
                &expected_i[..expected_rows],
                "{label}"
            );
            assert_eq!(
                collection.d_j_atom_type.as_slice(),
                &expected_j[..expected_rows],
                "{label}"
            );
            assert_eq!(
                collection.d_k_atom_type.as_slice(),
                &expected_k[..expected_rows],
                "{label}"
            );
            assert_eq!(
                collection.d_angle_type.as_slice(),
                &expected_angle[..expected_rows],
                "{label}"
            );
            let actual_values: Vec<_> = collection
                .d_params
                .iter()
                .map(|value| [value.ka.to_bits(), value.theta0.to_bits()])
                .collect();
            assert_eq!(
                actual_values.as_slice(),
                &expected_values[..expected_rows],
                "{label}"
            );
        }
        assert_eq!(constructor_calls, 4);
    }

    #[test]
    fn mmff_angle_constructor_reports_eight_frozen_input_safety_errors() {
        let invalid_rows = [
            (
                "*fixed\nX\t10\t7\t20\t1.25\t90\n",
                0,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                },
            ),
            (
                "*fixed\n0\tX\t7\t20\t1.25\t90\n",
                1,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                },
            ),
            (
                "*fixed\n0\t10\tX\t20\t1.25\t90\n",
                2,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                },
            ),
            (
                "*fixed\n0\t10\t7\tX\t1.25\t90\n",
                3,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                },
            ),
            (
                "*fixed\n0\t10\t7\t20\tX\t90\n",
                4,
                MmffParamParseCause::InvalidFloat {
                    cell: "X".to_owned(),
                },
            ),
            (
                "*fixed\n0\t10\t7\t20\t1.25\tX\n",
                5,
                MmffParamParseCause::InvalidFloat {
                    cell: "X".to_owned(),
                },
            ),
        ];
        let mut error_calls = 0;
        for (source, column, cause) in invalid_rows {
            let input_before = source.as_bytes().to_vec();
            let result = MmffAngleCollection::from_text(source);
            error_calls += 1;
            assert_eq!(source.as_bytes(), input_before.as_slice());
            assert_eq!(
                result,
                Err(MmffParamParseError {
                    table: MmffParamTable::Angle,
                    line: 2,
                    column,
                    cause,
                })
            );
        }

        let missing_theta = "*fixed\n0\t1\t2\t3\t1.25\n";
        let missing_theta_before = missing_theta.as_bytes().to_vec();
        let missing_theta_result = MmffAngleCollection::from_text(missing_theta);
        error_calls += 1;
        assert_eq!(missing_theta.as_bytes(), missing_theta_before.as_slice());
        assert_eq!(
            missing_theta_result,
            Err(MmffParamParseError {
                table: MmffParamTable::Angle,
                line: 2,
                column: 5,
                cause: MmffParamParseCause::MissingToken,
            })
        );

        let blank = "*fixed\n\n";
        let blank_before = blank.as_bytes().to_vec();
        let blank_result = MmffAngleCollection::from_text(blank);
        error_calls += 1;
        assert_eq!(blank.as_bytes(), blank_before.as_slice());
        assert_eq!(
            blank_result,
            Err(MmffParamParseError {
                table: MmffParamTable::Angle,
                line: 2,
                column: 0,
                cause: MmffParamParseCause::EmptyProcessedLine,
            })
        );
        assert_eq!(error_calls, 8);
    }

    #[test]
    fn mmff_angle_constructor_ignores_suffix_and_empty_input_selects_asset() {
        let suffix_source = "*fixed\n0\t10\t7\t20\t1.25\t90\tX\n";
        let suffix_before = suffix_source.as_bytes().to_vec();
        let suffix_collection = MmffAngleCollection::from_text(suffix_source)
            .expect("the source consumes only six cells");
        assert_eq!(suffix_source.as_bytes(), suffix_before.as_slice());
        assert_eq!(suffix_collection.d_i_atom_type.as_slice(), &[10_u8]);
        assert_eq!(suffix_collection.d_j_atom_type.as_slice(), &[7_u8]);
        assert_eq!(suffix_collection.d_k_atom_type.as_slice(), &[20_u8]);
        assert_eq!(suffix_collection.d_angle_type.as_slice(), &[0_u8]);
        assert_eq!(
            [
                suffix_collection.d_params[0].ka.to_bits(),
                suffix_collection.d_params[0].theta0.to_bits(),
            ],
            [0x3ff4_0000_0000_0000, 0x4056_8000_0000_0000]
        );

        const ASSET: &str = include_str!("default_angle.tsv");
        assert_eq!(ASSET.len(), 66833);
        assert_eq!(
            sha256_for_fixed_asset_test(ASSET.as_bytes()),
            [
                0x53, 0xe6, 0xe3, 0x85, 0xeb, 0xa6, 0x7d, 0x02, 0x0e, 0x7f, 0xe0, 0x0f, 0x9c, 0x95,
                0x20, 0x18, 0x1f, 0x6c, 0x96, 0x4d, 0x1f, 0x72, 0x82, 0x7a, 0x63, 0x8d, 0x66, 0x82,
                0x5f, 0x48, 0x3b, 0x70,
            ]
        );
        assert_eq!(
            ASSET.lines().filter(|line| !line.starts_with('*')).count(),
            2342
        );

        let empty_source = String::new();
        let empty_before = empty_source.as_bytes().to_vec();
        let defaults = MmffAngleCollection::from_text(&empty_source)
            .expect("empty constructor input selects the frozen default asset");
        assert_eq!(empty_source.as_bytes(), empty_before.as_slice());
        assert_eq!(defaults.d_i_atom_type.len(), 2342);
        assert_eq!(defaults.d_j_atom_type.len(), 2342);
        assert_eq!(defaults.d_k_atom_type.len(), 2342);
        assert_eq!(defaults.d_angle_type.len(), 2342);
        assert_eq!(defaults.d_params.len(), 2342);

        let assert_sentinel = |index: usize, keys: [u8; 4], values: [u64; 2]| {
            assert_eq!(
                [
                    defaults.d_angle_type[index],
                    defaults.d_i_atom_type[index],
                    defaults.d_j_atom_type[index],
                    defaults.d_k_atom_type[index],
                ],
                keys
            );
            assert_eq!(
                [
                    defaults.d_params[index].ka.to_bits(),
                    defaults.d_params[index].theta0.to_bits(),
                ],
                values
            );
        };
        assert_sentinel(
            0,
            [0, 0, 1, 0],
            [0x0000_0000_0000_0000, 0x405b_3999_9999_999a],
        );
        assert_sentinel(
            1171,
            [0, 1, 20, 26],
            [0x3fe7_126e_978d_4fdf, 0x405d_671a_9fbe_76c9],
        );
        assert_sentinel(
            2341,
            [0, 64, 82, 65],
            [0x3ff4_7ef9_db22_d0e5, 0x405c_3d1e_b851_eb85],
        );
    }

    #[derive(Debug, PartialEq, Eq)]
    struct MmffAngleDefaultSnapshot {
        collection_address: usize,
        key_rows: [Vec<u8>; 4],
        value_bits: Vec<[u64; 2]>,
        storage_addresses: [usize; 5],
        sentinel_rows: Vec<(usize, [u8; 4], [u64; 2])>,
    }

    fn mmff_angle_default_snapshot(collection: &MmffAngleCollection) -> MmffAngleDefaultSnapshot {
        const SENTINEL_INDICES: [usize; 3] = [0, 1171, 2341];
        let sentinel_rows = SENTINEL_INDICES
            .into_iter()
            .filter(|&index| index < collection.d_params.len())
            .map(|index| {
                let row = &collection.d_params[index];
                (
                    row as *const MmffAngle as usize,
                    [
                        collection.d_angle_type[index],
                        collection.d_i_atom_type[index],
                        collection.d_j_atom_type[index],
                        collection.d_k_atom_type[index],
                    ],
                    [row.ka.to_bits(), row.theta0.to_bits()],
                )
            })
            .collect();

        MmffAngleDefaultSnapshot {
            collection_address: collection as *const MmffAngleCollection as usize,
            key_rows: [
                collection.d_i_atom_type.clone(),
                collection.d_j_atom_type.clone(),
                collection.d_k_atom_type.clone(),
                collection.d_angle_type.clone(),
            ],
            value_bits: collection
                .d_params
                .iter()
                .map(|row| [row.ka.to_bits(), row.theta0.to_bits()])
                .collect(),
            storage_addresses: [
                collection.d_i_atom_type.as_ptr() as usize,
                collection.d_j_atom_type.as_ptr() as usize,
                collection.d_k_atom_type.as_ptr() as usize,
                collection.d_angle_type.as_ptr() as usize,
                collection.d_params.as_ptr() as usize,
            ],
            sentinel_rows,
        }
    }

    #[test]
    fn mmff_angle_defaults_asset_sentinels_and_128_warm_borrows() {
        const ASSET: &str = include_str!("default_angle.tsv");
        assert_eq!(ASSET.len(), 66833);
        assert!(ASSET.ends_with('\n'));
        assert_eq!(
            sha256_for_fixed_asset_test(ASSET.as_bytes()),
            [
                0x53, 0xe6, 0xe3, 0x85, 0xeb, 0xa6, 0x7d, 0x02, 0x0e, 0x7f, 0xe0, 0x0f, 0x9c, 0x95,
                0x20, 0x18, 0x1f, 0x6c, 0x96, 0x4d, 0x1f, 0x72, 0x82, 0x7a, 0x63, 0x8d, 0x66, 0x82,
                0x5f, 0x48, 0x3b, 0x70,
            ]
        );
        assert_eq!(fnv1a64(ASSET.as_bytes()), 0xe8f2_e205_1e05_0f47);
        assert_eq!(
            ASSET.lines().filter(|line| !line.starts_with('*')).count(),
            2342
        );

        assert_eq!(DEFAULT_MMFF_ANGLE_CONSTRUCTIONS.load(Ordering::Relaxed), 0);
        let custom_source = String::from("*custom\n0\t10\t7\t20\t7.25\t180\n");
        let custom_source_before = custom_source.as_bytes().to_vec();
        let custom = MmffAngleCollection::from_text(&custom_source)
            .expect("custom Angle construction remains independent");
        let custom_before = mmff_angle_default_snapshot(&custom);
        assert_eq!(custom_source.as_bytes(), custom_source_before.as_slice());
        assert_eq!(custom.d_i_atom_type.as_slice(), &[10]);
        assert_eq!(custom.d_j_atom_type.as_slice(), &[7]);
        assert_eq!(custom.d_k_atom_type.as_slice(), &[20]);
        assert_eq!(custom.d_angle_type.as_slice(), &[0]);
        assert_eq!(
            [
                custom.d_params[0].ka.to_bits(),
                custom.d_params[0].theta0.to_bits()
            ],
            [0x401d_0000_0000_0000, 0x4066_8000_0000_0000]
        );
        assert_eq!(DEFAULT_MMFF_ANGLE_CONSTRUCTIONS.load(Ordering::Relaxed), 0);

        let defaults = default_mmff_angle().expect("fixed Angle default asset must parse");
        assert_eq!(DEFAULT_MMFF_ANGLE_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(defaults.d_i_atom_type.len(), 2342);
        assert_eq!(defaults.d_j_atom_type.len(), 2342);
        assert_eq!(defaults.d_k_atom_type.len(), 2342);
        assert_eq!(defaults.d_angle_type.len(), 2342);
        assert_eq!(defaults.d_params.len(), 2342);

        let assert_sentinel = |index: usize, keys: [u8; 4], values: [u64; 2]| {
            assert_eq!(
                [
                    defaults.d_angle_type[index],
                    defaults.d_i_atom_type[index],
                    defaults.d_j_atom_type[index],
                    defaults.d_k_atom_type[index],
                ],
                keys
            );
            assert_eq!(
                [
                    defaults.d_params[index].ka.to_bits(),
                    defaults.d_params[index].theta0.to_bits(),
                ],
                values
            );
        };
        assert_sentinel(
            0,
            [0, 0, 1, 0],
            [0x0000_0000_0000_0000, 0x405b_3999_9999_999a],
        );
        assert_sentinel(
            1171,
            [0, 1, 20, 26],
            [0x3fe7_126e_978d_4fdf, 0x405d_671a_9fbe_76c9],
        );
        assert_sentinel(
            2341,
            [0, 64, 82, 65],
            [0x3ff4_7ef9_db22_d0e5, 0x405c_3d1e_b851_eb85],
        );

        let defaults_before = mmff_angle_default_snapshot(defaults);
        assert_eq!(defaults_before.sentinel_rows.len(), 3);
        assert_eq!(defaults_before.sentinel_rows[0].1, [0, 0, 1, 0]);
        assert_eq!(
            defaults_before.sentinel_rows[0].2,
            [0x0000_0000_0000_0000, 0x405b_3999_9999_999a]
        );
        assert_eq!(defaults_before.sentinel_rows[1].1, [0, 1, 20, 26]);
        assert_eq!(
            defaults_before.sentinel_rows[1].2,
            [0x3fe7_126e_978d_4fdf, 0x405d_671a_9fbe_76c9]
        );
        assert_eq!(defaults_before.sentinel_rows[2].1, [0, 64, 82, 65]);
        assert_eq!(
            defaults_before.sentinel_rows[2].2,
            [0x3ff4_7ef9_db22_d0e5, 0x405c_3d1e_b851_eb85]
        );
        assert_ne!(
            custom_before.storage_addresses,
            defaults_before.storage_addresses
        );
        assert_ne!(
            custom_before.collection_address,
            defaults_before.collection_address
        );
        assert_eq!(DEFAULT_MMFF_ANGLE_CONSTRUCTIONS.load(Ordering::Relaxed), 1);

        let post_default_custom_source =
            String::from("*custom after default\n1\t11\t7\t21\t8\t170\n");
        let post_default_custom_before = post_default_custom_source.as_bytes().to_vec();
        let post_default_custom = MmffAngleCollection::from_text(&post_default_custom_source)
            .expect("custom construction after default leaves the cache independent");
        let post_default_custom_snapshot = mmff_angle_default_snapshot(&post_default_custom);
        assert_eq!(
            post_default_custom_source.as_bytes(),
            post_default_custom_before.as_slice()
        );
        assert_eq!(DEFAULT_MMFF_ANGLE_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(mmff_angle_default_snapshot(defaults), defaults_before);
        assert_ne!(
            post_default_custom_snapshot.collection_address,
            defaults_before.collection_address
        );

        let warm_calls = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(8);
            for _ in 0..8 {
                let expected = &defaults_before;
                handles.push(scope.spawn(move || {
                    let mut calls = 0;
                    for _ in 0..16 {
                        let before = mmff_angle_default_snapshot(defaults);
                        assert_eq!(&before, expected);
                        let cached =
                            default_mmff_angle().expect("cached Angle default remains valid");
                        calls += 1;
                        let after = mmff_angle_default_snapshot(cached);
                        assert_eq!(
                            cached as *const MmffAngleCollection as usize,
                            expected.collection_address
                        );
                        assert_eq!(&after, &before);
                        assert_eq!(&after, expected);
                        assert_eq!(DEFAULT_MMFF_ANGLE_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
                    }
                    calls
                }));
            }

            handles.into_iter().fold(0, |total, handle| {
                total + handle.join().expect("scoped Angle default worker")
            })
        });
        assert_eq!(warm_calls, 128);
        assert_eq!(mmff_angle_default_snapshot(defaults), defaults_before);
        assert_eq!(mmff_angle_default_snapshot(&custom), custom_before);
        assert_eq!(
            mmff_angle_default_snapshot(&post_default_custom),
            post_default_custom_snapshot
        );
        assert_eq!(DEFAULT_MMFF_ANGLE_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
    }

    const MMFF_OOP_LOOKUP_DEF_SOURCE: &str = concat!(
        "*fixed\n",
        "X\t1\t10\t11\t12\t13\n",
        "X\t2\t20\t21\t22\t23\n",
        "X\t3\t30\t31\t32\t33\n",
    );

    const MMFF_OOP_LOOKUP_SOURCE: &str = concat!(
        "*fixed\n",
        "10\t7\t20\t20\t4.5\n",
        "10\t7\t20\t30\t1.25\n",
        "10\t7\t20\t30\t9\n",
        "11\t7\t21\t31\t0.5\n",
        "12\t7\t22\t32\t0.75\n",
        "13\t7\t23\t33\t0.125\n",
        "10\t8\t20\t30\t3.75\n",
    );

    #[derive(Clone, Copy, Debug)]
    enum ExpectedMmffOopLookup {
        Hit { row: usize, koop_bits: u64 },
        Miss,
        MissingDefinition { atom_type: u32 },
    }

    struct MmffOopLookupSnapshot {
        table_source_bytes: Vec<u8>,
        table_source_address: usize,
        def_source_bytes: Vec<u8>,
        def_source_address: usize,
        def_eq_rows: Vec<[u8; 4]>,
        def_rows_address: usize,
        key_rows: [Vec<u8>; 4],
        koop_bits: Vec<u64>,
        storage_addresses: [usize; 5],
    }

    impl MmffOopLookupSnapshot {
        fn capture(
            table_source: &str,
            def_source: &str,
            definitions: &MmffDefCollection,
            collection: &MmffOopCollection,
        ) -> Self {
            Self {
                table_source_bytes: table_source.as_bytes().to_vec(),
                table_source_address: table_source.as_ptr() as usize,
                def_source_bytes: def_source.as_bytes().to_vec(),
                def_source_address: def_source.as_ptr() as usize,
                def_eq_rows: definitions
                    .d_params
                    .iter()
                    .map(|row| row.eq_level)
                    .collect(),
                def_rows_address: definitions.d_params.as_ptr() as usize,
                key_rows: [
                    collection.d_i_atom_type.clone(),
                    collection.d_j_atom_type.clone(),
                    collection.d_k_atom_type.clone(),
                    collection.d_l_atom_type.clone(),
                ],
                koop_bits: collection
                    .d_params
                    .iter()
                    .map(|row| row.koop.to_bits())
                    .collect(),
                storage_addresses: [
                    collection.d_i_atom_type.as_ptr() as usize,
                    collection.d_j_atom_type.as_ptr() as usize,
                    collection.d_k_atom_type.as_ptr() as usize,
                    collection.d_l_atom_type.as_ptr() as usize,
                    collection.d_params.as_ptr() as usize,
                ],
            }
        }

        fn assert_unchanged(
            &self,
            table_source: &str,
            def_source: &str,
            definitions: &MmffDefCollection,
            collection: &MmffOopCollection,
            context: &str,
        ) {
            assert_eq!(
                table_source.as_bytes(),
                self.table_source_bytes,
                "{context}"
            );
            assert_eq!(
                table_source.as_ptr() as usize,
                self.table_source_address,
                "{context}"
            );
            assert_eq!(def_source.as_bytes(), self.def_source_bytes, "{context}");
            assert_eq!(
                def_source.as_ptr() as usize,
                self.def_source_address,
                "{context}"
            );
            assert_eq!(
                definitions
                    .d_params
                    .iter()
                    .map(|row| row.eq_level)
                    .collect::<Vec<_>>(),
                self.def_eq_rows,
                "{context}"
            );
            assert_eq!(
                definitions.d_params.as_ptr() as usize,
                self.def_rows_address,
                "{context}"
            );
            assert_eq!(
                collection.d_i_atom_type.as_slice(),
                self.key_rows[0].as_slice(),
                "{context}"
            );
            assert_eq!(
                collection.d_j_atom_type.as_slice(),
                self.key_rows[1].as_slice(),
                "{context}"
            );
            assert_eq!(
                collection.d_k_atom_type.as_slice(),
                self.key_rows[2].as_slice(),
                "{context}"
            );
            assert_eq!(
                collection.d_l_atom_type.as_slice(),
                self.key_rows[3].as_slice(),
                "{context}"
            );
            assert_eq!(
                collection
                    .d_params
                    .iter()
                    .map(|row| row.koop.to_bits())
                    .collect::<Vec<_>>(),
                self.koop_bits,
                "{context}"
            );
            assert_eq!(
                [
                    collection.d_i_atom_type.as_ptr() as usize,
                    collection.d_j_atom_type.as_ptr() as usize,
                    collection.d_k_atom_type.as_ptr() as usize,
                    collection.d_l_atom_type.as_ptr() as usize,
                    collection.d_params.as_ptr() as usize,
                ],
                self.storage_addresses,
                "{context}"
            );
        }
    }

    fn assert_mmff_oop_lookup_result(
        result: Result<Option<&MmffOop>, MmffOopLookupError>,
        expected: ExpectedMmffOopLookup,
        collection: &MmffOopCollection,
        context: &str,
    ) {
        match (result, expected) {
            (Ok(Some(actual)), ExpectedMmffOopLookup::Hit { row, koop_bits }) => {
                let expected_row = collection
                    .d_params
                    .get(row)
                    .unwrap_or_else(|| panic!("{context}: expected row {row} is absent"));
                assert!(std::ptr::eq(actual, expected_row), "{context}");
                assert_eq!(actual.koop.to_bits(), koop_bits, "{context}");
            }
            (Ok(None), ExpectedMmffOopLookup::Miss) => {}
            (Err(error), ExpectedMmffOopLookup::MissingDefinition { atom_type }) => {
                assert_eq!(
                    error,
                    MmffOopLookupError::MissingDefinition { atom_type },
                    "{context}"
                );
                assert_eq!(
                    error.to_string(),
                    format!("missing MMFF definition for atom type {atom_type}"),
                    "{context}"
                );
                assert!(std::error::Error::source(&error).is_none(), "{context}");
            }
            (actual, expected) => {
                panic!("{context}: actual {actual:?}, expected {expected:?}")
            }
        }
    }

    fn run_mmff_oop_lookup_with_checkpoint(
        source: &str,
        query: [u32; 4],
        expected: ExpectedMmffOopLookup,
        context: &str,
    ) {
        let table_source = source.to_owned();
        let def_source = MMFF_OOP_LOOKUP_DEF_SOURCE.to_owned();
        let collection = MmffOopCollection::from_text(false, &table_source)
            .unwrap_or_else(|error| panic!("{context}: fixed OOP source failed: {error}"));
        let definitions = MmffDefCollection::from_text(&def_source)
            .unwrap_or_else(|error| panic!("{context}: fixed Def source failed: {error}"));
        let snapshot =
            MmffOopLookupSnapshot::capture(&table_source, &def_source, &definitions, &collection);

        let result = collection.get(&definitions, query[0], query[1], query[2], query[3]);
        snapshot.assert_unchanged(
            &table_source,
            &def_source,
            &definitions,
            &collection,
            context,
        );
        assert_mmff_oop_lookup_result(result, expected, &collection, context);
    }

    fn mmff_oop_constructor_with_input_checkpoint(
        is_mmff_s: bool,
        source: &str,
    ) -> Result<MmffOopCollection, MmffParamParseError> {
        let input_address = source.as_ptr() as usize;
        let input_bytes = source.as_bytes().to_vec();
        let result = MmffOopCollection::from_text(is_mmff_s, source);
        assert_eq!(source.as_ptr() as usize, input_address);
        assert_eq!(source.as_bytes(), input_bytes.as_slice());
        result
    }

    #[test]
    fn mmff_oop_constructor_formats_preserve_literal_arrays_and_bits() {
        const SOURCE: &str = concat!(
            "*fixed\n",
            "10\t7\t20\t20\t4.5\n",
            "10\t7\t20\t30\t1.25\n",
            "10\t7\t20\t30\t9\n",
            "11\t7\t21\t31\t0.5\n",
            "12\t7\t22\t32\t0.75\n",
            "13\t7\t23\t33\t0.125\n",
            "10\t8\t20\t30\t3.75\n",
        );
        let crlf = SOURCE.replace('\n', "\r\n");
        let doubled_tabs = SOURCE.replace('\t', "\t\t");
        let unterminated = SOURCE
            .strip_suffix('\n')
            .expect("the fixed source fixture ends with LF")
            .to_owned();
        let formats = [
            ("LF", SOURCE.to_owned(), 7),
            ("CRLF", crlf, 7),
            ("doubled TAB", doubled_tabs, 7),
            ("unterminated final row", unterminated, 6),
        ];
        const EXPECTED_I: [u8; 7] = [10, 10, 10, 11, 12, 13, 10];
        const EXPECTED_J: [u8; 7] = [7, 7, 7, 7, 7, 7, 8];
        const EXPECTED_K: [u8; 7] = [20, 20, 20, 21, 22, 23, 20];
        const EXPECTED_L: [u8; 7] = [20, 30, 30, 31, 32, 33, 30];
        const EXPECTED_KOOP_BITS: [u64; 7] = [
            0x4012_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x4022_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x3fe8_0000_0000_0000,
            0x3fc0_0000_0000_0000,
            0x400e_0000_0000_0000,
        ];

        let mut constructor_calls = 0;
        for is_mmff_s in [false, true] {
            for (label, source, expected_rows) in &formats {
                assert!(constructor_calls < 8, "unexpected extra format constructor");
                let collection = mmff_oop_constructor_with_input_checkpoint(is_mmff_s, source)
                    .unwrap_or_else(|error| panic!("{label}, is_mmff_s={is_mmff_s}: {error}"));
                constructor_calls += 1;

                assert_eq!(collection.d_i_atom_type.len(), *expected_rows, "{label}");
                assert_eq!(collection.d_j_atom_type.len(), *expected_rows, "{label}");
                assert_eq!(collection.d_k_atom_type.len(), *expected_rows, "{label}");
                assert_eq!(collection.d_l_atom_type.len(), *expected_rows, "{label}");
                assert_eq!(collection.d_params.len(), *expected_rows, "{label}");
                assert_eq!(
                    collection.d_i_atom_type.as_slice(),
                    &EXPECTED_I[..*expected_rows],
                    "{label}"
                );
                assert_eq!(
                    collection.d_j_atom_type.as_slice(),
                    &EXPECTED_J[..*expected_rows],
                    "{label}"
                );
                assert_eq!(
                    collection.d_k_atom_type.as_slice(),
                    &EXPECTED_K[..*expected_rows],
                    "{label}"
                );
                assert_eq!(
                    collection.d_l_atom_type.as_slice(),
                    &EXPECTED_L[..*expected_rows],
                    "{label}"
                );
                let actual_koop_bits: Vec<_> = collection
                    .d_params
                    .iter()
                    .map(|params| params.koop.to_bits())
                    .collect();
                assert_eq!(
                    actual_koop_bits.as_slice(),
                    &EXPECTED_KOOP_BITS[..*expected_rows],
                    "{label}, is_mmff_s={is_mmff_s}"
                );
            }
        }
        assert_eq!(constructor_calls, 8);
    }

    #[test]
    fn mmff_oop_constructor_returns_fourteen_exact_input_safety_errors() {
        let invalid_rows = [
            (
                "*fixed\nX\t7\t20\t20\t4.5\n",
                0,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                },
            ),
            (
                "*fixed\n10\tX\t20\t20\t4.5\n",
                1,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                },
            ),
            (
                "*fixed\n10\t7\tX\t20\t4.5\n",
                2,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                },
            ),
            (
                "*fixed\n10\t7\t20\tX\t4.5\n",
                3,
                MmffParamParseCause::InvalidUnsigned {
                    cell: "X".to_owned(),
                },
            ),
            (
                "*fixed\n10\t7\t20\t20\tX\n",
                4,
                MmffParamParseCause::InvalidFloat {
                    cell: "X".to_owned(),
                },
            ),
            (
                "*fixed\n10\t7\t20\t20\n",
                4,
                MmffParamParseCause::MissingToken,
            ),
            ("*fixed\n\n", 0, MmffParamParseCause::EmptyProcessedLine),
        ];

        let mut error_calls = 0;
        for is_mmff_s in [false, true] {
            for (source, column, cause) in &invalid_rows {
                assert!(error_calls < 14, "unexpected extra invalid constructor");
                let error = mmff_oop_constructor_with_input_checkpoint(is_mmff_s, source)
                    .expect_err("the fixed invalid row must retain its typed error");
                error_calls += 1;
                assert_eq!(
                    error,
                    MmffParamParseError {
                        table: MmffParamTable::Oop,
                        line: 2,
                        column: *column,
                        cause: cause.clone(),
                    },
                    "is_mmff_s={is_mmff_s}, input={source:?}"
                );
            }
        }
        assert_eq!(error_calls, 14);
    }

    #[test]
    fn mmff_oop_constructor_ignores_suffix_and_selects_both_assets() {
        const SUFFIX_SOURCE: &str = "*fixed\n10\t7\t20\t20\t4.5\tX\n";
        let mut suffix_calls = 0;
        for is_mmff_s in [false, true] {
            let collection = mmff_oop_constructor_with_input_checkpoint(is_mmff_s, SUFFIX_SOURCE)
                .expect("the source consumes only five OOP fields");
            suffix_calls += 1;
            assert_eq!(collection.d_i_atom_type.as_slice(), &[10]);
            assert_eq!(collection.d_j_atom_type.as_slice(), &[7]);
            assert_eq!(collection.d_k_atom_type.as_slice(), &[20]);
            assert_eq!(collection.d_l_atom_type.as_slice(), &[20]);
            assert_eq!(collection.d_params[0].koop.to_bits(), 0x4012_0000_0000_0000);
        }
        assert_eq!(suffix_calls, 2);

        const SENTINELS: [(usize, [u8; 4], u64, u64); 5] = [
            (
                0,
                [0, 2, 0, 0],
                0x3f94_7ae1_47ae_147b,
                0x3f94_7ae1_47ae_147b,
            ),
            (
                38,
                [0, 10, 0, 0],
                0xbf94_7ae1_47ae_147b,
                0x3f8e_b851_eb85_1eb8,
            ),
            (
                58,
                [6, 37, 37, 37],
                0x3fa8_9374_bc6a_7efa,
                0x3fa8_9374_bc6a_7efa,
            ),
            (
                66,
                [0, 40, 0, 0],
                0xbf74_7ae1_47ae_147b,
                0x3f9e_b851_eb85_1eb8,
            ),
            (116, [0, 82, 0, 0], 0, 0),
        ];
        let mut empty_asset_calls = 0;
        for is_mmff_s in [false, true] {
            let empty_source = String::new();
            let collection = mmff_oop_constructor_with_input_checkpoint(is_mmff_s, &empty_source)
                .expect("empty input selects its frozen source asset");
            empty_asset_calls += 1;
            assert_eq!(collection.d_i_atom_type.len(), 117);
            assert_eq!(collection.d_j_atom_type.len(), 117);
            assert_eq!(collection.d_k_atom_type.len(), 117);
            assert_eq!(collection.d_l_atom_type.len(), 117);
            assert_eq!(collection.d_params.len(), 117);

            for (row, keys, regular_bits, mmff_s_bits) in SENTINELS {
                assert_eq!(
                    [
                        collection.d_i_atom_type[row],
                        collection.d_j_atom_type[row],
                        collection.d_k_atom_type[row],
                        collection.d_l_atom_type[row],
                    ],
                    keys,
                    "row {row}, is_mmff_s={is_mmff_s}"
                );
                assert_eq!(
                    collection.d_params[row].koop.to_bits(),
                    if is_mmff_s { mmff_s_bits } else { regular_bits },
                    "row {row}, is_mmff_s={is_mmff_s}"
                );
            }
        }
        assert_eq!(empty_asset_calls, 2);
    }

    fn mmff_tor_constructor_with_input_checkpoint(
        is_mmff_s: bool,
        source: &str,
    ) -> Result<MmffTorCollection, MmffParamParseError> {
        let input_address = source.as_ptr() as usize;
        let input_bytes = source.as_bytes().to_vec();
        let result = MmffTorCollection::from_text(is_mmff_s, source);
        assert_eq!(source.as_ptr() as usize, input_address);
        assert_eq!(source.as_bytes(), input_bytes.as_slice());
        result
    }

    #[test]
    fn mmff_tor_constructor_32_frozen_calls_preserve_inputs_and_outputs() {
        const SOURCE: &str = concat!(
            "*fixed\n",
            "2\t10\t7\t8\t20\t3.75\t-0.0\t2.5\n",
            "5\t10\t7\t8\t20\t1.25\t-0.0\t2.5\n",
            "5\t10\t7\t8\t20\t9\t-0.0\t2.5\n",
            "2\t11\t7\t8\t23\t4.5\t-0.0\t2.5\n",
            "5\t11\t7\t8\t23\t0.5\t-0.0\t2.5\n",
            "2\t13\t7\t8\t21\t5.5\t-0.0\t2.5\n",
            "5\t13\t7\t8\t21\t0.75\t-0.0\t2.5\n",
            "2\t13\t7\t8\t23\t6.5\t-0.0\t2.5\n",
            "5\t13\t7\t8\t23\t0.125\t-0.0\t2.5\n",
        );
        const EXPECTED_TOR: [u8; 9] = [2, 5, 5, 2, 5, 2, 5, 2, 5];
        const EXPECTED_I: [u8; 9] = [10, 10, 10, 11, 11, 13, 13, 13, 13];
        const EXPECTED_J: [u8; 9] = [7, 7, 7, 7, 7, 7, 7, 7, 7];
        const EXPECTED_K: [u8; 9] = [8, 8, 8, 8, 8, 8, 8, 8, 8];
        const EXPECTED_L: [u8; 9] = [20, 20, 20, 23, 23, 21, 21, 23, 23];
        const EXPECTED_V1_BITS: [u64; 9] = [
            0x400e_0000_0000_0000,
            0x3ff4_0000_0000_0000,
            0x4022_0000_0000_0000,
            0x4012_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4016_0000_0000_0000,
            0x3fe8_0000_0000_0000,
            0x401a_0000_0000_0000,
            0x3fc0_0000_0000_0000,
        ];
        const NEGATIVE_ZERO_BITS: u64 = 0x8000_0000_0000_0000;
        const TWO_POINT_FIVE_BITS: u64 = 0x4004_0000_0000_0000;

        let crlf = SOURCE.replace('\n', "\r\n");
        let doubled_tabs = SOURCE.replace('\t', "\t\t");
        let unterminated = SOURCE
            .strip_suffix('\n')
            .expect("the fixed table ends with LF")
            .to_owned();
        let formats = [
            ("LF", SOURCE.to_owned(), 9),
            ("CRLF", crlf, 9),
            ("doubled TAB", doubled_tabs, 9),
            ("unterminated final row", unterminated, 8),
        ];

        let mut constructor_calls = 0;
        for is_mmff_s in [false, true] {
            for (label, source, expected_rows) in &formats {
                assert!(constructor_calls < 32, "unexpected extra Tor constructor");
                let result = mmff_tor_constructor_with_input_checkpoint(is_mmff_s, source);
                constructor_calls += 1;
                let collection = result
                    .unwrap_or_else(|error| panic!("{label}, is_mmff_s={is_mmff_s}: {error}"));

                assert_eq!(collection.d_tor_type.len(), *expected_rows, "{label}");
                assert_eq!(collection.d_i_atom_type.len(), *expected_rows, "{label}");
                assert_eq!(collection.d_j_atom_type.len(), *expected_rows, "{label}");
                assert_eq!(collection.d_k_atom_type.len(), *expected_rows, "{label}");
                assert_eq!(collection.d_l_atom_type.len(), *expected_rows, "{label}");
                assert_eq!(collection.d_params.len(), *expected_rows, "{label}");
                assert_eq!(
                    collection.d_tor_type.as_slice(),
                    &EXPECTED_TOR[..*expected_rows],
                    "{label}"
                );
                assert_eq!(
                    collection.d_i_atom_type.as_slice(),
                    &EXPECTED_I[..*expected_rows],
                    "{label}"
                );
                assert_eq!(
                    collection.d_j_atom_type.as_slice(),
                    &EXPECTED_J[..*expected_rows],
                    "{label}"
                );
                assert_eq!(
                    collection.d_k_atom_type.as_slice(),
                    &EXPECTED_K[..*expected_rows],
                    "{label}"
                );
                assert_eq!(
                    collection.d_l_atom_type.as_slice(),
                    &EXPECTED_L[..*expected_rows],
                    "{label}"
                );
                let actual_v1_bits: Vec<_> = collection
                    .d_params
                    .iter()
                    .map(|params| params.v1.to_bits())
                    .collect();
                let actual_v2_bits: Vec<_> = collection
                    .d_params
                    .iter()
                    .map(|params| params.v2.to_bits())
                    .collect();
                let actual_v3_bits: Vec<_> = collection
                    .d_params
                    .iter()
                    .map(|params| params.v3.to_bits())
                    .collect();
                assert_eq!(
                    actual_v1_bits.as_slice(),
                    &EXPECTED_V1_BITS[..*expected_rows],
                    "{label}, is_mmff_s={is_mmff_s}"
                );
                assert_eq!(
                    actual_v2_bits.as_slice(),
                    &[NEGATIVE_ZERO_BITS; 9][..*expected_rows],
                    "{label}, is_mmff_s={is_mmff_s}"
                );
                assert_eq!(
                    actual_v3_bits.as_slice(),
                    &[TWO_POINT_FIVE_BITS; 9][..*expected_rows],
                    "{label}, is_mmff_s={is_mmff_s}"
                );
            }
        }

        const INVALID_CELLS: [(&str, usize); 8] = [
            ("*fixed\nX\t10\t7\t8\t20\t3.75\t-0.0\t2.5\n", 0),
            ("*fixed\n2\tX\t7\t8\t20\t3.75\t-0.0\t2.5\n", 1),
            ("*fixed\n2\t10\tX\t8\t20\t3.75\t-0.0\t2.5\n", 2),
            ("*fixed\n2\t10\t7\tX\t20\t3.75\t-0.0\t2.5\n", 3),
            ("*fixed\n2\t10\t7\t8\tX\t3.75\t-0.0\t2.5\n", 4),
            ("*fixed\n2\t10\t7\t8\t20\tX\t-0.0\t2.5\n", 5),
            ("*fixed\n2\t10\t7\t8\t20\t3.75\tX\t2.5\n", 6),
            ("*fixed\n2\t10\t7\t8\t20\t3.75\t-0.0\tX\n", 7),
        ];
        const MISSING_LAST_FLOAT: &str = "*fixed\n2\t10\t7\t8\t20\t3.75\t-0.0\n";
        const BLANK_PROCESSED_ROW: &str = "*fixed\n\n";
        for is_mmff_s in [false, true] {
            for (source, column) in INVALID_CELLS {
                assert!(constructor_calls < 32, "unexpected extra Tor constructor");
                let result = mmff_tor_constructor_with_input_checkpoint(is_mmff_s, source);
                constructor_calls += 1;
                let cause = if column < 5 {
                    MmffParamParseCause::InvalidUnsigned {
                        cell: "X".to_owned(),
                    }
                } else {
                    MmffParamParseCause::InvalidFloat {
                        cell: "X".to_owned(),
                    }
                };
                assert_eq!(
                    result.expect_err("the fixed invalid cell must return a typed error"),
                    MmffParamParseError {
                        table: MmffParamTable::Tor,
                        line: 2,
                        column,
                        cause,
                    },
                    "is_mmff_s={is_mmff_s}, column={column}"
                );
            }

            for (source, column, cause) in [
                (MISSING_LAST_FLOAT, 7, MmffParamParseCause::MissingToken),
                (
                    BLANK_PROCESSED_ROW,
                    0,
                    MmffParamParseCause::EmptyProcessedLine,
                ),
            ] {
                assert!(constructor_calls < 32, "unexpected extra Tor constructor");
                let result = mmff_tor_constructor_with_input_checkpoint(is_mmff_s, source);
                constructor_calls += 1;
                assert_eq!(
                    result.expect_err("the fixed malformed row must return a typed error"),
                    MmffParamParseError {
                        table: MmffParamTable::Tor,
                        line: 2,
                        column,
                        cause,
                    },
                    "is_mmff_s={is_mmff_s}, input={source:?}"
                );
            }
        }

        const SUFFIX_SOURCE: &str = "*fixed\n2\t10\t7\t8\t20\t3.75\t-0.0\t2.5\tX\n";
        for is_mmff_s in [false, true] {
            assert!(constructor_calls < 32, "unexpected extra Tor constructor");
            let result = mmff_tor_constructor_with_input_checkpoint(is_mmff_s, SUFFIX_SOURCE);
            constructor_calls += 1;
            let collection = result.expect("the source ignores cells after V3");
            assert_eq!(collection.d_tor_type.as_slice(), &[2]);
            assert_eq!(collection.d_i_atom_type.as_slice(), &[10]);
            assert_eq!(collection.d_j_atom_type.as_slice(), &[7]);
            assert_eq!(collection.d_k_atom_type.as_slice(), &[8]);
            assert_eq!(collection.d_l_atom_type.as_slice(), &[20]);
            assert_eq!(collection.d_params[0].v1.to_bits(), 0x400e_0000_0000_0000);
            assert_eq!(collection.d_params[0].v2.to_bits(), NEGATIVE_ZERO_BITS);
            assert_eq!(collection.d_params[0].v3.to_bits(), TWO_POINT_FIVE_BITS);
        }

        const REGULAR_ASSET: &str = include_str!("default_tor.tsv");
        const MMFF_S_ASSET: &str = include_str!("default_tor_s.tsv");
        assert_eq!(REGULAR_ASSET.len(), 39_795);
        assert_eq!(MMFF_S_ASSET.len(), 39_888);
        assert!(REGULAR_ASSET.ends_with('\n'));
        assert!(MMFF_S_ASSET.ends_with('\n'));
        assert_eq!(
            sha256_for_fixed_asset_test(REGULAR_ASSET.as_bytes()),
            [
                0xfe, 0xc5, 0x3b, 0x92, 0x8f, 0xee, 0x6c, 0xf4, 0x89, 0xec, 0x6d, 0xd3, 0xe7, 0x36,
                0x97, 0xf6, 0x26, 0xb9, 0x44, 0xef, 0x4c, 0x9d, 0x76, 0xd8, 0x84, 0x13, 0x97, 0xec,
                0x8c, 0x09, 0x17, 0xa6,
            ]
        );
        assert_eq!(
            sha256_for_fixed_asset_test(MMFF_S_ASSET.as_bytes()),
            [
                0xb5, 0xf2, 0xd7, 0x58, 0x89, 0x30, 0xee, 0x19, 0x37, 0x12, 0xc3, 0x50, 0x12, 0xe5,
                0x4c, 0xfd, 0x83, 0xb5, 0xd3, 0x27, 0xb3, 0xa3, 0xf5, 0x26, 0x04, 0xcc, 0x26, 0x55,
                0x6a, 0x37, 0x57, 0x6a,
            ]
        );
        const SENTINELS: [(usize, [u8; 5], [u64; 3], [u64; 3]); 4] = [
            (
                0,
                [0, 0, 1, 1, 0],
                [0, 0, 0x3fd3_3333_3333_3333],
                [0, 0, 0x3fd3_3333_3333_3333],
            ),
            (
                22,
                [0, 5, 1, 1, 10],
                [0, 0, 0x3fdb_53f7_ced9_1687],
                [0, 0, 0x3fda_c083_126e_978d],
            ),
            (
                463,
                [0, 0, 3, 54, 0],
                [0, 0x4020_0000_0000_0000, 0],
                [0, 0x4020_0000_0000_0000, 0],
            ),
            (
                925,
                [0, 0, 80, 81, 0],
                [0, 0x4010_0000_0000_0000, 0],
                [0, 0x4010_0000_0000_0000, 0],
            ),
        ];
        for is_mmff_s in [false, true] {
            assert!(constructor_calls < 32, "unexpected extra Tor constructor");
            let empty_source = String::new();
            let result = mmff_tor_constructor_with_input_checkpoint(is_mmff_s, &empty_source);
            constructor_calls += 1;
            let collection = result.expect("empty input selects the frozen Tor asset");
            assert_eq!(collection.d_tor_type.len(), 926);
            assert_eq!(collection.d_i_atom_type.len(), 926);
            assert_eq!(collection.d_j_atom_type.len(), 926);
            assert_eq!(collection.d_k_atom_type.len(), 926);
            assert_eq!(collection.d_l_atom_type.len(), 926);
            assert_eq!(collection.d_params.len(), 926);
            for (row, keys, regular_bits, mmff_s_bits) in SENTINELS {
                assert_eq!(
                    [
                        collection.d_tor_type[row],
                        collection.d_i_atom_type[row],
                        collection.d_j_atom_type[row],
                        collection.d_k_atom_type[row],
                        collection.d_l_atom_type[row],
                    ],
                    keys,
                    "row {row}, is_mmff_s={is_mmff_s}"
                );
                let expected_bits = if is_mmff_s { mmff_s_bits } else { regular_bits };
                assert_eq!(
                    [
                        collection.d_params[row].v1.to_bits(),
                        collection.d_params[row].v2.to_bits(),
                        collection.d_params[row].v3.to_bits(),
                    ],
                    expected_bits,
                    "row {row}, is_mmff_s={is_mmff_s}"
                );
            }
        }
        assert_eq!(constructor_calls, 32);
    }

    const MMFF_TOR_LOOKUP_DEF_SOURCE: &str = concat!(
        "*fixed\n",
        "X\t1\t10\t11\t12\t13\n",
        "X\t2\t20\t21\t22\t23\n",
    );

    const MMFF_TOR_LOOKUP_SOURCE: &str = concat!(
        "*fixed\n",
        "2\t10\t7\t8\t20\t3.75\t-0.0\t2.5\n",
        "5\t10\t7\t8\t20\t1.25\t-0.0\t2.5\n",
        "5\t10\t7\t8\t20\t9\t-0.0\t2.5\n",
        "2\t11\t7\t8\t23\t4.5\t-0.0\t2.5\n",
        "5\t11\t7\t8\t23\t0.5\t-0.0\t2.5\n",
        "2\t13\t7\t8\t21\t5.5\t-0.0\t2.5\n",
        "5\t13\t7\t8\t21\t0.75\t-0.0\t2.5\n",
        "2\t13\t7\t8\t23\t6.5\t-0.0\t2.5\n",
        "5\t13\t7\t8\t23\t0.125\t-0.0\t2.5\n",
    );

    type MmffTorLookupQuery = ((u32, u32), u32, u32, u32, u32);

    const MMFF_TOR_BITS_3_75: [u64; 3] = [
        0x400e_0000_0000_0000,
        0x8000_0000_0000_0000,
        0x4004_0000_0000_0000,
    ];
    const MMFF_TOR_BITS_1_25: [u64; 3] = [
        0x3ff4_0000_0000_0000,
        0x8000_0000_0000_0000,
        0x4004_0000_0000_0000,
    ];
    const MMFF_TOR_BITS_0_5: [u64; 3] = [
        0x3fe0_0000_0000_0000,
        0x8000_0000_0000_0000,
        0x4004_0000_0000_0000,
    ];
    const MMFF_TOR_BITS_0_75: [u64; 3] = [
        0x3fe8_0000_0000_0000,
        0x8000_0000_0000_0000,
        0x4004_0000_0000_0000,
    ];
    const MMFF_TOR_BITS_0_125: [u64; 3] = [
        0x3fc0_0000_0000_0000,
        0x8000_0000_0000_0000,
        0x4004_0000_0000_0000,
    ];
    const MMFF_TOR_BITS_4_5: [u64; 3] = [
        0x4012_0000_0000_0000,
        0x8000_0000_0000_0000,
        0x4004_0000_0000_0000,
    ];
    const MMFF_TOR_BITS_5_5: [u64; 3] = [
        0x4016_0000_0000_0000,
        0x8000_0000_0000_0000,
        0x4004_0000_0000_0000,
    ];
    const MMFF_TOR_BITS_6_5: [u64; 3] = [
        0x401a_0000_0000_0000,
        0x8000_0000_0000_0000,
        0x4004_0000_0000_0000,
    ];

    #[derive(Clone, Copy, Debug)]
    enum ExpectedMmffTorLookup {
        Hit {
            returned_type: u32,
            row: usize,
            value_bits: [u64; 3],
        },
        Miss {
            returned_type: u32,
        },
        MissingDefinition {
            atom_type: u32,
        },
        InvalidEquivalentLevel {
            level: usize,
        },
    }

    struct MmffTorLookupSnapshot {
        table_source_bytes: Vec<u8>,
        table_source_address: usize,
        def_source_bytes: Vec<u8>,
        def_source_address: usize,
        def_eq_rows: Vec<[u8; 4]>,
        def_rows_address: usize,
        key_rows: [Vec<u8>; 5],
        value_bits: Vec<[u64; 3]>,
        storage_addresses: [usize; 6],
    }

    impl MmffTorLookupSnapshot {
        fn capture(
            table_source: &str,
            def_source: &str,
            definitions: &MmffDefCollection,
            collection: &MmffTorCollection,
        ) -> Self {
            Self {
                table_source_bytes: table_source.as_bytes().to_vec(),
                table_source_address: table_source.as_ptr() as usize,
                def_source_bytes: def_source.as_bytes().to_vec(),
                def_source_address: def_source.as_ptr() as usize,
                def_eq_rows: definitions
                    .d_params
                    .iter()
                    .map(|row| row.eq_level)
                    .collect(),
                def_rows_address: definitions.d_params.as_ptr() as usize,
                key_rows: [
                    collection.d_tor_type.clone(),
                    collection.d_i_atom_type.clone(),
                    collection.d_j_atom_type.clone(),
                    collection.d_k_atom_type.clone(),
                    collection.d_l_atom_type.clone(),
                ],
                value_bits: collection
                    .d_params
                    .iter()
                    .map(|row| [row.v1.to_bits(), row.v2.to_bits(), row.v3.to_bits()])
                    .collect(),
                storage_addresses: [
                    collection.d_tor_type.as_ptr() as usize,
                    collection.d_i_atom_type.as_ptr() as usize,
                    collection.d_j_atom_type.as_ptr() as usize,
                    collection.d_k_atom_type.as_ptr() as usize,
                    collection.d_l_atom_type.as_ptr() as usize,
                    collection.d_params.as_ptr() as usize,
                ],
            }
        }

        fn record_unchanged(
            &self,
            table_source: &str,
            def_source: &str,
            definitions: &MmffDefCollection,
            collection: &MmffTorCollection,
            context: &str,
            discrepancies: &mut Vec<String>,
        ) {
            if table_source.as_bytes() != self.table_source_bytes.as_slice() {
                discrepancies.push(format!("{context}: table source bytes changed"));
            }
            if table_source.as_ptr() as usize != self.table_source_address {
                discrepancies.push(format!("{context}: table source address changed"));
            }
            if def_source.as_bytes() != self.def_source_bytes.as_slice() {
                discrepancies.push(format!("{context}: Def source bytes changed"));
            }
            if def_source.as_ptr() as usize != self.def_source_address {
                discrepancies.push(format!("{context}: Def source address changed"));
            }

            let actual_def_rows: Vec<_> = definitions
                .d_params
                .iter()
                .map(|row| row.eq_level)
                .collect();
            if actual_def_rows != self.def_eq_rows {
                discrepancies.push(format!("{context}: Def equivalence rows changed"));
            }
            if definitions.d_params.as_ptr() as usize != self.def_rows_address {
                discrepancies.push(format!("{context}: Def buffer address changed"));
            }

            let actual_keys = [
                collection.d_tor_type.as_slice(),
                collection.d_i_atom_type.as_slice(),
                collection.d_j_atom_type.as_slice(),
                collection.d_k_atom_type.as_slice(),
                collection.d_l_atom_type.as_slice(),
            ];
            for (index, (actual, expected)) in actual_keys
                .into_iter()
                .zip(self.key_rows.iter())
                .enumerate()
            {
                if actual != expected.as_slice() {
                    discrepancies.push(format!("{context}: key buffer {index} changed"));
                }
            }

            let actual_value_bits: Vec<_> = collection
                .d_params
                .iter()
                .map(|row| [row.v1.to_bits(), row.v2.to_bits(), row.v3.to_bits()])
                .collect();
            if actual_value_bits != self.value_bits {
                discrepancies.push(format!("{context}: value bits changed"));
            }

            let actual_addresses = [
                collection.d_tor_type.as_ptr() as usize,
                collection.d_i_atom_type.as_ptr() as usize,
                collection.d_j_atom_type.as_ptr() as usize,
                collection.d_k_atom_type.as_ptr() as usize,
                collection.d_l_atom_type.as_ptr() as usize,
                collection.d_params.as_ptr() as usize,
            ];
            if actual_addresses != self.storage_addresses {
                discrepancies.push(format!("{context}: collection buffer address changed"));
            }
        }
    }

    fn record_mmff_tor_lookup_result(
        result: Result<(u32, Option<&MmffTor>), MmffTorLookupError>,
        expected: ExpectedMmffTorLookup,
        collection: &MmffTorCollection,
        context: &str,
        discrepancies: &mut Vec<String>,
    ) {
        match (result, expected) {
            (
                Ok((returned_type, Some(actual))),
                ExpectedMmffTorLookup::Hit {
                    returned_type: expected_type,
                    row,
                    value_bits,
                },
            ) => {
                if returned_type != expected_type {
                    discrepancies.push(format!(
                        "{context}: returned type {returned_type}, expected {expected_type}"
                    ));
                }
                let expected_row = &collection.d_params[row];
                if !std::ptr::eq(actual, expected_row) {
                    discrepancies.push(format!(
                        "{context}: result pointer is not literal row {row}"
                    ));
                }
                let actual_bits = [
                    actual.v1.to_bits(),
                    actual.v2.to_bits(),
                    actual.v3.to_bits(),
                ];
                if actual_bits != value_bits {
                    discrepancies.push(format!(
                        "{context}: row {row} bits {actual_bits:?}, expected {value_bits:?}"
                    ));
                }
            }
            (
                Ok((returned_type, None)),
                ExpectedMmffTorLookup::Miss {
                    returned_type: expected,
                },
            ) if returned_type == expected => {}
            (Err(error), ExpectedMmffTorLookup::MissingDefinition { atom_type }) => {
                let expected_error = MmffTorLookupError::MissingDefinition { atom_type };
                let expected_display = format!("missing MMFF definition for atom type {atom_type}");
                if error != expected_error {
                    discrepancies.push(format!(
                        "{context}: error {error:?}, expected {expected_error:?}"
                    ));
                }
                if error.to_string() != expected_display {
                    discrepancies.push(format!(
                        "{context}: display {:?}, expected {expected_display:?}",
                        error.to_string()
                    ));
                }
                if std::error::Error::source(&error).is_some() {
                    discrepancies.push(format!("{context}: missing-Def error has a source"));
                }
            }
            (Err(error), ExpectedMmffTorLookup::InvalidEquivalentLevel { level }) => {
                let expected_error = MmffTorLookupError::InvalidEquivalentLevel { level };
                let expected_display =
                    format!("MMFF torsion equivalence level {level} is outside 0..4");
                if error != expected_error {
                    discrepancies.push(format!(
                        "{context}: error {error:?}, expected {expected_error:?}"
                    ));
                }
                if error.to_string() != expected_display {
                    discrepancies.push(format!(
                        "{context}: display {:?}, expected {expected_display:?}",
                        error.to_string()
                    ));
                }
                if std::error::Error::source(&error).is_some() {
                    discrepancies.push(format!("{context}: invalid-level error has a source"));
                }
            }
            (actual, expected) => discrepancies.push(format!(
                "{context}: actual {actual:?}, expected {expected:?}"
            )),
        }
    }

    fn run_mmff_tor_lookup_with_checkpoint(
        source: &str,
        query: MmffTorLookupQuery,
        expected: ExpectedMmffTorLookup,
        context: &str,
        lookup_calls: &mut usize,
        discrepancies: &mut Vec<String>,
    ) {
        let table_source = source.to_owned();
        let def_source = MMFF_TOR_LOOKUP_DEF_SOURCE.to_owned();
        let collection = MmffTorCollection::from_text(false, &table_source)
            .unwrap_or_else(|error| panic!("{context}: fixed Tor source failed: {error}"));
        let definitions = MmffDefCollection::from_text(&def_source)
            .unwrap_or_else(|error| panic!("{context}: fixed Def source failed: {error}"));
        let expected_def_rows = [[10, 11, 12, 13], [20, 21, 22, 23]];
        let actual_def_rows: Vec<_> = definitions
            .d_params
            .iter()
            .map(|row| row.eq_level)
            .collect();
        assert_eq!(
            actual_def_rows.as_slice(),
            expected_def_rows.as_slice(),
            "{context}: fixed positional Def fixture"
        );
        if let ExpectedMmffTorLookup::Hit { row, .. } = expected {
            assert!(
                row < collection.d_params.len(),
                "{context}: expected literal row {row} is absent before lookup"
            );
        }
        let snapshot =
            MmffTorLookupSnapshot::capture(&table_source, &def_source, &definitions, &collection);

        let ((first_type, second_type), i_atom_type, j_atom_type, k_atom_type, l_atom_type) = query;
        let result = collection.get(
            &definitions,
            (first_type, second_type),
            i_atom_type,
            j_atom_type,
            k_atom_type,
            l_atom_type,
        );
        snapshot.record_unchanged(
            &table_source,
            &def_source,
            &definitions,
            &collection,
            context,
            discrepancies,
        );
        *lookup_calls += 1;
        record_mmff_tor_lookup_result(result, expected, &collection, context, discrepancies);
    }

    #[test]
    fn mmff_tor_lookup_111_frozen_calls_preserve_source_results_and_state() {
        const FORMAT_QUERIES: [(MmffTorLookupQuery, ExpectedMmffTorLookup); 12] = [
            (
                ((5, 2), 1, 7, 8, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 5,
                    row: 1,
                    value_bits: MMFF_TOR_BITS_1_25,
                },
            ),
            (
                ((5, 2), 2, 8, 7, 1),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 5,
                    row: 1,
                    value_bits: MMFF_TOR_BITS_1_25,
                },
            ),
            (
                ((2, 0), 1, 7, 8, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 2,
                    row: 0,
                    value_bits: MMFF_TOR_BITS_3_75,
                },
            ),
            (
                ((1, 2), 1, 7, 8, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 2,
                    row: 0,
                    value_bits: MMFF_TOR_BITS_3_75,
                },
            ),
            (
                ((4, 0), 1, 7, 8, 2),
                ExpectedMmffTorLookup::Miss { returned_type: 0 },
            ),
            (
                ((5, 0), 1, 7, 8, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 5,
                    row: 1,
                    value_bits: MMFF_TOR_BITS_1_25,
                },
            ),
            (
                ((517, 0), 1, 7, 8, 2),
                ExpectedMmffTorLookup::Miss { returned_type: 0 },
            ),
            (
                ((4, 0), 1, 263, 8, 2),
                ExpectedMmffTorLookup::Miss { returned_type: 0 },
            ),
            (
                ((4, 0), 1, 7, 264, 2),
                ExpectedMmffTorLookup::Miss { returned_type: 0 },
            ),
            (
                ((4, 0), 0, 263, 8, 2),
                ExpectedMmffTorLookup::MissingDefinition { atom_type: 0 },
            ),
            (
                ((4, 0), 1, 7, 8, 0),
                ExpectedMmffTorLookup::MissingDefinition { atom_type: 0 },
            ),
            (
                ((4, 0), 257, 7, 8, 2),
                ExpectedMmffTorLookup::MissingDefinition { atom_type: 257 },
            ),
        ];

        let crlf = MMFF_TOR_LOOKUP_SOURCE.replace('\n', "\r\n");
        let doubled_tabs = MMFF_TOR_LOOKUP_SOURCE.replace('\t', "\t\t");
        let unterminated = MMFF_TOR_LOOKUP_SOURCE
            .strip_suffix('\n')
            .expect("the fixed Tor lookup source ends in LF")
            .to_owned();
        let formats = [
            ("LF", MMFF_TOR_LOOKUP_SOURCE.to_owned()),
            ("CRLF", crlf),
            ("doubled TAB", doubled_tabs),
            ("unterminated final row", unterminated),
        ];
        let mut lookup_calls = 0;
        let mut discrepancies = Vec::new();

        for (format_name, source) in &formats {
            for (query_index, (query, expected)) in FORMAT_QUERIES.iter().copied().enumerate() {
                let context = format!("{format_name} query {}", query_index + 1);
                run_mmff_tor_lookup_with_checkpoint(
                    source,
                    query,
                    expected,
                    &context,
                    &mut lookup_calls,
                    &mut discrepancies,
                );
            }
        }

        const PRIORITY_STAGE_ROWS: [&str; 4] = [
            "4\t10\t7\t8\t20\t1.25\t-0.0\t2.5\n",
            "4\t11\t7\t8\t23\t0.5\t-0.0\t2.5\n",
            "4\t13\t7\t8\t21\t0.75\t-0.0\t2.5\n",
            "4\t13\t7\t8\t23\t0.125\t-0.0\t2.5\n",
        ];
        const PRIORITY_STAGE_SETS: [&[usize]; 4] = [&[0, 1, 2, 3], &[1, 2, 3], &[2, 3], &[3]];
        const PRIORITY_FORWARD: [(usize, [u64; 3]); 4] = [
            (0, MMFF_TOR_BITS_1_25),
            (0, MMFF_TOR_BITS_0_5),
            (0, MMFF_TOR_BITS_0_75),
            (0, MMFF_TOR_BITS_0_125),
        ];
        const PRIORITY_REVERSE: [(usize, [u64; 3]); 4] = [
            (0, MMFF_TOR_BITS_1_25),
            (1, MMFF_TOR_BITS_0_75),
            (0, MMFF_TOR_BITS_0_75),
            (0, MMFF_TOR_BITS_0_125),
        ];

        for (set_index, stages) in PRIORITY_STAGE_SETS.iter().enumerate() {
            let mut source = String::from("*fixed\n");
            for stage in *stages {
                source.push_str(PRIORITY_STAGE_ROWS[*stage]);
            }
            for (direction_name, query_base, expected_rows) in [
                ("forward", (1, 7, 8, 2), PRIORITY_FORWARD[set_index]),
                ("reverse", (2, 8, 7, 1), PRIORITY_REVERSE[set_index]),
            ] {
                for second_type in [0, 2] {
                    let query = (
                        (4, second_type),
                        query_base.0,
                        query_base.1,
                        query_base.2,
                        query_base.3,
                    );
                    let context =
                        format!("priority set {set_index}, {direction_name}, second={second_type}");
                    run_mmff_tor_lookup_with_checkpoint(
                        &source,
                        query,
                        ExpectedMmffTorLookup::Hit {
                            returned_type: 4,
                            row: expected_rows.0,
                            value_bits: expected_rows.1,
                        },
                        &context,
                        &mut lookup_calls,
                        &mut discrepancies,
                    );
                }
            }
        }

        let branch_cases: [(&str, &str, (u32, u32), ExpectedMmffTorLookup); 12] = [
            (
                "B0",
                "*fixed\n5\t10\t7\t8\t20\t1.25\t-0.0\t2.5\n",
                (5, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 5,
                    row: 0,
                    value_bits: MMFF_TOR_BITS_1_25,
                },
            ),
            (
                "B1",
                "*fixed\n5\t11\t7\t8\t23\t0.5\t-0.0\t2.5\n",
                (5, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 5,
                    row: 0,
                    value_bits: MMFF_TOR_BITS_0_5,
                },
            ),
            (
                "B2",
                "*fixed\n5\t13\t7\t8\t21\t0.75\t-0.0\t2.5\n",
                (5, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 5,
                    row: 0,
                    value_bits: MMFF_TOR_BITS_0_75,
                },
            ),
            (
                "B3",
                "*fixed\n5\t13\t7\t8\t23\t0.125\t-0.0\t2.5\n",
                (5, 0),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 5,
                    row: 0,
                    value_bits: MMFF_TOR_BITS_0_125,
                },
            ),
            (
                "B4",
                concat!(
                    "*fixed\n",
                    "2\t10\t7\t8\t20\t3.75\t-0.0\t2.5\n",
                    "5\t13\t7\t8\t23\t0.125\t-0.0\t2.5\n",
                ),
                (5, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 2,
                    row: 0,
                    value_bits: MMFF_TOR_BITS_3_75,
                },
            ),
            (
                "B5",
                "*fixed\n5\t13\t7\t8\t23\t0.125\t-0.0\t2.5\n",
                (5, 2),
                ExpectedMmffTorLookup::InvalidEquivalentLevel { level: 4 },
            ),
            (
                "B6",
                "*fixed\n2\t11\t7\t8\t23\t4.5\t-0.0\t2.5\n",
                (5, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 2,
                    row: 0,
                    value_bits: MMFF_TOR_BITS_4_5,
                },
            ),
            (
                "B7",
                "*fixed\n",
                (4, 0),
                ExpectedMmffTorLookup::Miss { returned_type: 0 },
            ),
            (
                "B8",
                "*fixed\n",
                (5, 0),
                ExpectedMmffTorLookup::Miss { returned_type: 0 },
            ),
            (
                "B9",
                "*fixed\n",
                (5, 2),
                ExpectedMmffTorLookup::InvalidEquivalentLevel { level: 4 },
            ),
            (
                "B10",
                "*fixed\n4\t13\t7\t8\t23\t0.125\t-0.0\t2.5\n",
                (4, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 4,
                    row: 0,
                    value_bits: MMFF_TOR_BITS_0_125,
                },
            ),
            (
                "B11",
                "*fixed\n2\t10\t7\t8\t20\t3.75\t-0.0\t2.5\n",
                (4, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 2,
                    row: 0,
                    value_bits: MMFF_TOR_BITS_3_75,
                },
            ),
        ];

        for (name, source, tor_pair, expected) in branch_cases {
            for repeat in 0..2 {
                let context = format!("branch {name}, repeat {}", repeat + 1);
                run_mmff_tor_lookup_with_checkpoint(
                    source,
                    (tor_pair, 1, 7, 8, 2),
                    expected,
                    &context,
                    &mut lookup_calls,
                    &mut discrepancies,
                );
            }
        }

        const CAST_SOURCE: &str = "*fixed\n517\t266\t263\t264\t276\t0\t1.25\t2.5\n";
        const CAST_QUERIES: [(MmffTorLookupQuery, ExpectedMmffTorLookup); 7] = [
            (
                ((5, 0), 1, 7, 8, 2),
                ExpectedMmffTorLookup::Hit {
                    returned_type: 5,
                    row: 0,
                    value_bits: [
                        0x0000_0000_0000_0000,
                        0x3ff4_0000_0000_0000,
                        0x4004_0000_0000_0000,
                    ],
                },
            ),
            (
                ((517, 0), 1, 7, 8, 2),
                ExpectedMmffTorLookup::Miss { returned_type: 0 },
            ),
            (
                ((4, 0), 1, 263, 8, 2),
                ExpectedMmffTorLookup::Miss { returned_type: 0 },
            ),
            (
                ((4, 0), 1, 7, 264, 2),
                ExpectedMmffTorLookup::Miss { returned_type: 0 },
            ),
            (
                ((4, 0), 1, 263, 264, 2),
                ExpectedMmffTorLookup::Miss { returned_type: 0 },
            ),
            (
                ((5, 0), 257, 7, 8, 2),
                ExpectedMmffTorLookup::MissingDefinition { atom_type: 257 },
            ),
            (
                ((5, 0), 1, 7, 8, 258),
                ExpectedMmffTorLookup::MissingDefinition { atom_type: 258 },
            ),
        ];
        for (index, (query, expected)) in CAST_QUERIES.iter().copied().enumerate() {
            let context = format!("stored-key cast control {}", index + 1);
            run_mmff_tor_lookup_with_checkpoint(
                CAST_SOURCE,
                query,
                expected,
                &context,
                &mut lookup_calls,
                &mut discrepancies,
            );
        }

        const EQUAL_CENTRAL_SOURCE: &str = "*fixed\n5\t10\t7\t7\t20\t1.25\t-0.0\t2.5\n";
        for (direction_name, query) in [
            ("forward", ((5, 0), 1, 7, 7, 2)),
            ("reverse", ((5, 0), 2, 7, 7, 1)),
        ] {
            for repeat in 0..2 {
                let context = format!("equal central {direction_name}, repeat {}", repeat + 1);
                run_mmff_tor_lookup_with_checkpoint(
                    EQUAL_CENTRAL_SOURCE,
                    query,
                    ExpectedMmffTorLookup::Hit {
                        returned_type: 5,
                        row: 0,
                        value_bits: MMFF_TOR_BITS_1_25,
                    },
                    &context,
                    &mut lookup_calls,
                    &mut discrepancies,
                );
            }
        }

        const MISSING_DEFINITION_QUERIES: [(u32, u32); 6] =
            [(0, 2), (99, 2), (1, 0), (1, 99), (0, 0), (99, 99)];
        for (central_j, central_k) in [(263, 8), (264, 8)] {
            for (index, (i_type, l_type)) in MISSING_DEFINITION_QUERIES.iter().copied().enumerate()
            {
                let expected_atom_type = if i_type == 0 || i_type == 99 {
                    i_type
                } else {
                    l_type
                };
                let context = format!(
                    "missing Def {}, central=({central_j},{central_k})",
                    index + 1
                );
                run_mmff_tor_lookup_with_checkpoint(
                    "*fixed\n",
                    ((4, 0), i_type, central_j, central_k, l_type),
                    ExpectedMmffTorLookup::MissingDefinition {
                        atom_type: expected_atom_type,
                    },
                    &context,
                    &mut lookup_calls,
                    &mut discrepancies,
                );
            }
        }

        assert_eq!(lookup_calls, 111, "all frozen actual lookup calls executed");
        assert!(
            discrepancies.is_empty(),
            "source lookup discrepancies: {discrepancies:#?}"
        );
    }

    #[test]
    fn mmff_oop_lookup_matrix_48_queries_preserves_source_order_and_state() {
        const TERMINAL_PERMUTATIONS: [(u32, u32, u32); 6] = [
            (1, 2, 3),
            (1, 3, 2),
            (2, 1, 3),
            (2, 3, 1),
            (3, 1, 2),
            (3, 2, 1),
        ];
        let crlf = MMFF_OOP_LOOKUP_SOURCE.replace('\n', "\r\n");
        let doubled_tabs = MMFF_OOP_LOOKUP_SOURCE.replace('\t', "\t\t");
        let unterminated = MMFF_OOP_LOOKUP_SOURCE
            .strip_suffix('\n')
            .expect("the fixed OOP lookup source ends in LF")
            .to_owned();
        let formats = [
            ("LF", MMFF_OOP_LOOKUP_SOURCE.to_owned(), 7),
            ("CRLF", crlf, 7),
            ("doubled TAB", doubled_tabs, 7),
            ("unterminated final row", unterminated, 6),
        ];
        let mut lookup_calls = 0;

        for (format_label, source, expected_rows) in formats {
            for (permutation_index, (i, k, l)) in TERMINAL_PERMUTATIONS.into_iter().enumerate() {
                assert!(lookup_calls < 48, "unexpected extra format lookup");
                run_mmff_oop_lookup_with_checkpoint(
                    &source,
                    [i, 7, k, l],
                    ExpectedMmffOopLookup::Hit {
                        row: 1,
                        koop_bits: 0x3ff4_0000_0000_0000,
                    },
                    &format!("{format_label} terminal permutation {permutation_index}"),
                );
                lookup_calls += 1;
            }

            let row_6_expectation = if expected_rows == 7 {
                ExpectedMmffOopLookup::Hit {
                    row: 6,
                    koop_bits: 0x400e_0000_0000_0000,
                }
            } else {
                ExpectedMmffOopLookup::Miss
            };
            let remaining_queries = [
                (
                    [1, 7, 2, 2],
                    ExpectedMmffOopLookup::Hit {
                        row: 0,
                        koop_bits: 0x4012_0000_0000_0000,
                    },
                    "first duplicate-equivalence row",
                ),
                (
                    [1, 8, 2, 3],
                    row_6_expectation,
                    "last row with LF-sensitive presence",
                ),
                (
                    [0, 9, 99, 99],
                    ExpectedMmffOopLookup::Miss,
                    "central-j miss precedes missing Def rows",
                ),
                (
                    [1, 263, 2, 3],
                    ExpectedMmffOopLookup::Miss,
                    "full-width central-j miss",
                ),
                (
                    [1, 7, 2, 3],
                    ExpectedMmffOopLookup::Hit {
                        row: 1,
                        koop_bits: 0x3ff4_0000_0000_0000,
                    },
                    "source-first duplicate row",
                ),
                (
                    [3, 7, 3, 3],
                    ExpectedMmffOopLookup::Miss,
                    "valid Def positions with no OOP row",
                ),
            ];
            for (query_index, (query, expected, label)) in remaining_queries.into_iter().enumerate()
            {
                assert!(lookup_calls < 48, "unexpected extra format lookup");
                run_mmff_oop_lookup_with_checkpoint(
                    &source,
                    query,
                    expected,
                    &format!("{format_label} query {query_index}: {label}"),
                );
                lookup_calls += 1;
            }
        }

        assert_eq!(lookup_calls, 48);
    }

    #[test]
    fn mmff_oop_lookup_stage_priority_24_terminal_permutations() {
        const TERMINAL_PERMUTATIONS: [(u32, u32, u32); 6] = [
            (1, 2, 3),
            (1, 3, 2),
            (2, 1, 3),
            (2, 3, 1),
            (3, 1, 2),
            (3, 2, 1),
        ];
        let stage_tables: [(&str, u64); 4] = [
            (
                concat!(
                    "*all stages\n",
                    "10\t7\t20\t30\t1.25\n",
                    "11\t7\t21\t31\t0.5\n",
                    "12\t7\t22\t32\t0.75\n",
                    "13\t7\t23\t33\t0.125\n",
                ),
                0x3ff4_0000_0000_0000,
            ),
            (
                concat!(
                    "*stage one onward\n",
                    "11\t7\t21\t31\t0.5\n",
                    "12\t7\t22\t32\t0.75\n",
                    "13\t7\t23\t33\t0.125\n",
                ),
                0x3fe0_0000_0000_0000,
            ),
            (
                concat!(
                    "*stage two onward\n",
                    "12\t7\t22\t32\t0.75\n",
                    "13\t7\t23\t33\t0.125\n",
                ),
                0x3fe8_0000_0000_0000,
            ),
            (
                "*stage three only\n13\t7\t23\t33\t0.125\n",
                0x3fc0_0000_0000_0000,
            ),
        ];
        let mut lookup_calls = 0;

        for (table_index, (source, expected_koop_bits)) in stage_tables.into_iter().enumerate() {
            for (permutation_index, (i, k, l)) in TERMINAL_PERMUTATIONS.into_iter().enumerate() {
                assert!(lookup_calls < 24, "unexpected extra stage-priority lookup");
                run_mmff_oop_lookup_with_checkpoint(
                    source,
                    [i, 7, k, l],
                    ExpectedMmffOopLookup::Hit {
                        row: 0,
                        koop_bits: expected_koop_bits,
                    },
                    &format!("stage table {table_index}, permutation {permutation_index}"),
                );
                lookup_calls += 1;
            }
        }

        assert_eq!(lookup_calls, 24);
    }

    #[test]
    fn mmff_oop_lookup_cast_and_missing_def_controls_are_exact() {
        const CAST_SOURCE: &str = "266\t263\t276\t286\t-0\n";
        let cast_queries = [
            (
                [1, 7, 2, 3],
                ExpectedMmffOopLookup::Hit {
                    row: 0,
                    koop_bits: 0x8000_0000_0000_0000,
                },
                "wrapped stored keys and negative zero",
            ),
            (
                [1, 263, 2, 3],
                ExpectedMmffOopLookup::Miss,
                "full-width query is not narrowed",
            ),
            (
                [257, 7, 2, 3],
                ExpectedMmffOopLookup::MissingDefinition { atom_type: 257 },
                "full-width missing I definition",
            ),
            (
                [1, 7, 258, 3],
                ExpectedMmffOopLookup::MissingDefinition { atom_type: 258 },
                "full-width missing K definition",
            ),
            (
                [1, 7, 2, 259],
                ExpectedMmffOopLookup::MissingDefinition { atom_type: 259 },
                "full-width missing L definition",
            ),
        ];
        let mut lookup_calls = 0;
        for (query_index, (query, expected, label)) in cast_queries.into_iter().enumerate() {
            assert!(lookup_calls < 12, "unexpected extra cast-control lookup");
            run_mmff_oop_lookup_with_checkpoint(
                CAST_SOURCE,
                query,
                expected,
                &format!("cast query {query_index}: {label}"),
            );
            lookup_calls += 1;
        }

        let missing_definition_queries = [
            (
                [0, 7, 2, 3],
                ExpectedMmffOopLookup::MissingDefinition { atom_type: 0 },
                "I zero is reported before K/L",
            ),
            (
                [99, 7, 2, 3],
                ExpectedMmffOopLookup::MissingDefinition { atom_type: 99 },
                "I out of range",
            ),
            (
                [1, 7, 0, 3],
                ExpectedMmffOopLookup::MissingDefinition { atom_type: 0 },
                "K zero after valid I",
            ),
            (
                [1, 7, 99, 3],
                ExpectedMmffOopLookup::MissingDefinition { atom_type: 99 },
                "K out of range after valid I",
            ),
            (
                [1, 7, 2, 0],
                ExpectedMmffOopLookup::MissingDefinition { atom_type: 0 },
                "L zero after valid I and K",
            ),
            (
                [1, 7, 2, 99],
                ExpectedMmffOopLookup::MissingDefinition { atom_type: 99 },
                "L out of range after valid I and K",
            ),
            (
                [0, 7, 0, 0],
                ExpectedMmffOopLookup::MissingDefinition { atom_type: 0 },
                "I precedence when all Def positions are missing",
            ),
        ];
        for (query_index, (query, expected, label)) in
            missing_definition_queries.into_iter().enumerate()
        {
            assert!(lookup_calls < 12, "unexpected extra missing-Def lookup");
            run_mmff_oop_lookup_with_checkpoint(
                MMFF_OOP_LOOKUP_SOURCE,
                query,
                expected,
                &format!("missing-Def query {query_index}: {label}"),
            );
            lookup_calls += 1;
        }

        assert_eq!(lookup_calls, 12);
    }

    #[derive(Debug, PartialEq, Eq)]
    struct MmffOopDefaultSnapshot {
        collection_address: usize,
        key_rows: [Vec<u8>; 4],
        koop_bits: Vec<u64>,
        storage_addresses: [usize; 5],
        sentinel_rows: Vec<(usize, usize, [u8; 4], u64)>,
    }

    fn mmff_oop_default_snapshot(collection: &MmffOopCollection) -> MmffOopDefaultSnapshot {
        const SENTINEL_INDICES: [usize; 5] = [0, 38, 58, 66, 116];
        let sentinel_rows = SENTINEL_INDICES
            .into_iter()
            .filter(|&index| index < collection.d_params.len())
            .map(|index| {
                let row = &collection.d_params[index];
                (
                    index,
                    row as *const MmffOop as usize,
                    [
                        collection.d_i_atom_type[index],
                        collection.d_j_atom_type[index],
                        collection.d_k_atom_type[index],
                        collection.d_l_atom_type[index],
                    ],
                    row.koop.to_bits(),
                )
            })
            .collect();

        MmffOopDefaultSnapshot {
            collection_address: collection as *const MmffOopCollection as usize,
            key_rows: [
                collection.d_i_atom_type.clone(),
                collection.d_j_atom_type.clone(),
                collection.d_k_atom_type.clone(),
                collection.d_l_atom_type.clone(),
            ],
            koop_bits: collection
                .d_params
                .iter()
                .map(|row| row.koop.to_bits())
                .collect(),
            storage_addresses: [
                collection.d_i_atom_type.as_ptr() as usize,
                collection.d_j_atom_type.as_ptr() as usize,
                collection.d_k_atom_type.as_ptr() as usize,
                collection.d_l_atom_type.as_ptr() as usize,
                collection.d_params.as_ptr() as usize,
            ],
            sentinel_rows,
        }
    }

    #[test]
    fn mmff_oop_defaults_assets_sentinels_and_128_warm_borrows() {
        const REGULAR_ASSET: &str = include_str!("default_oop.tsv");
        const MMFF_S_ASSET: &str = include_str!("default_oop_s.tsv");
        assert_eq!(REGULAR_ASSET.len(), 3045);
        assert_eq!(MMFF_S_ASSET.len(), 3062);
        assert!(REGULAR_ASSET.ends_with('\n'));
        assert!(MMFF_S_ASSET.ends_with('\n'));
        assert_eq!(
            sha256_for_fixed_asset_test(REGULAR_ASSET.as_bytes()),
            [
                0x97, 0x50, 0x17, 0x18, 0x88, 0x43, 0x33, 0xa9, 0xfe, 0x35, 0x46, 0x0d, 0x02, 0x53,
                0x76, 0xf5, 0xeb, 0xd7, 0xf0, 0x64, 0x4c, 0xf2, 0x92, 0x07, 0xa7, 0xdb, 0x6a, 0x07,
                0x8e, 0x37, 0x2b, 0xf9,
            ]
        );
        assert_eq!(
            sha256_for_fixed_asset_test(MMFF_S_ASSET.as_bytes()),
            [
                0x84, 0x88, 0x07, 0xb0, 0x23, 0xcc, 0x31, 0xa2, 0x1a, 0xba, 0xbb, 0x77, 0x9b, 0x93,
                0x4a, 0xbf, 0xf0, 0xaf, 0x85, 0x47, 0xe5, 0xce, 0x8a, 0x6d, 0x2f, 0x85, 0xbc, 0x09,
                0x4b, 0x84, 0x3e, 0x32,
            ]
        );
        assert_eq!(
            REGULAR_ASSET
                .lines()
                .filter(|line| !line.starts_with('*'))
                .count(),
            117
        );
        assert_eq!(
            MMFF_S_ASSET
                .lines()
                .filter(|line| !line.starts_with('*'))
                .count(),
            117
        );

        assert_eq!(DEFAULT_MMFF_OOP_CONSTRUCTIONS.load(Ordering::Relaxed), 0);
        assert_eq!(DEFAULT_MMFF_OOP_S_CONSTRUCTIONS.load(Ordering::Relaxed), 0);

        const CUSTOM_SOURCE: &str = "*custom\n10\t7\t20\t30\t5.5\n";
        let custom_regular = mmff_oop_constructor_with_input_checkpoint(false, CUSTOM_SOURCE)
            .expect("regular custom construction is independent of cached defaults");
        let custom_regular_before = mmff_oop_default_snapshot(&custom_regular);
        let custom_s = mmff_oop_constructor_with_input_checkpoint(true, CUSTOM_SOURCE)
            .expect("MMFFs custom construction is independent of cached defaults");
        let custom_s_before = mmff_oop_default_snapshot(&custom_s);
        for custom in [&custom_regular, &custom_s] {
            assert_eq!(custom.d_i_atom_type.as_slice(), &[10]);
            assert_eq!(custom.d_j_atom_type.as_slice(), &[7]);
            assert_eq!(custom.d_k_atom_type.as_slice(), &[20]);
            assert_eq!(custom.d_l_atom_type.as_slice(), &[30]);
            assert_eq!(custom.d_params[0].koop.to_bits(), 0x4016_0000_0000_0000);
        }

        let empty_custom_regular = mmff_oop_constructor_with_input_checkpoint(false, "")
            .expect("empty custom constructor selects regular source asset independently");
        let empty_custom_regular_before = mmff_oop_default_snapshot(&empty_custom_regular);
        let empty_custom_s = mmff_oop_constructor_with_input_checkpoint(true, "")
            .expect("empty custom constructor selects MMFFs source asset independently");
        let empty_custom_s_before = mmff_oop_default_snapshot(&empty_custom_s);
        assert_eq!(empty_custom_regular_before.key_rows[0].len(), 117);
        assert_eq!(empty_custom_s_before.key_rows[0].len(), 117);
        assert_eq!(DEFAULT_MMFF_OOP_CONSTRUCTIONS.load(Ordering::Relaxed), 0);
        assert_eq!(DEFAULT_MMFF_OOP_S_CONSTRUCTIONS.load(Ordering::Relaxed), 0);

        let regular = default_mmff_oop(false).expect("regular OOP default asset parses");
        assert_eq!(DEFAULT_MMFF_OOP_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(DEFAULT_MMFF_OOP_S_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        let mmff_s = default_mmff_oop(true).expect("MMFFs OOP default asset parses");
        assert_eq!(DEFAULT_MMFF_OOP_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(DEFAULT_MMFF_OOP_S_CONSTRUCTIONS.load(Ordering::Relaxed), 1);

        for defaults in [regular, mmff_s] {
            assert_eq!(defaults.d_i_atom_type.len(), 117);
            assert_eq!(defaults.d_j_atom_type.len(), 117);
            assert_eq!(defaults.d_k_atom_type.len(), 117);
            assert_eq!(defaults.d_l_atom_type.len(), 117);
            assert_eq!(defaults.d_params.len(), 117);
        }
        let regular_before = mmff_oop_default_snapshot(regular);
        let mmff_s_before = mmff_oop_default_snapshot(mmff_s);
        assert_eq!(
            regular_before.key_rows,
            empty_custom_regular_before.key_rows
        );
        assert_eq!(
            regular_before.koop_bits,
            empty_custom_regular_before.koop_bits
        );
        assert_eq!(mmff_s_before.key_rows, empty_custom_s_before.key_rows);
        assert_eq!(mmff_s_before.koop_bits, empty_custom_s_before.koop_bits);

        const SENTINELS: [(usize, [u8; 4], u64, u64); 5] = [
            (
                0,
                [0, 2, 0, 0],
                0x3f94_7ae1_47ae_147b,
                0x3f94_7ae1_47ae_147b,
            ),
            (
                38,
                [0, 10, 0, 0],
                0xbf94_7ae1_47ae_147b,
                0x3f8e_b851_eb85_1eb8,
            ),
            (
                58,
                [6, 37, 37, 37],
                0x3fa8_9374_bc6a_7efa,
                0x3fa8_9374_bc6a_7efa,
            ),
            (
                66,
                [0, 40, 0, 0],
                0xbf74_7ae1_47ae_147b,
                0x3f9e_b851_eb85_1eb8,
            ),
            (116, [0, 82, 0, 0], 0, 0),
        ];
        assert_eq!(regular_before.sentinel_rows.len(), SENTINELS.len());
        assert_eq!(mmff_s_before.sentinel_rows.len(), SENTINELS.len());
        for (position, (index, keys, regular_bits, mmff_s_bits)) in
            SENTINELS.into_iter().enumerate()
        {
            let regular_row = regular_before.sentinel_rows[position];
            let mmff_s_row = mmff_s_before.sentinel_rows[position];
            assert_eq!(regular_row.0, index);
            assert_eq!(regular_row.2, keys);
            assert_eq!(regular_row.3, regular_bits);
            assert_eq!(mmff_s_row.0, index);
            assert_eq!(mmff_s_row.2, keys);
            assert_eq!(mmff_s_row.3, mmff_s_bits);
            assert_eq!(
                regular.d_params[index].koop.to_bits(),
                regular_bits,
                "regular OOP sentinel row {index}"
            );
            assert_eq!(
                mmff_s.d_params[index].koop.to_bits(),
                mmff_s_bits,
                "MMFFs OOP sentinel row {index}"
            );
        }

        assert_ne!(
            regular_before.collection_address,
            mmff_s_before.collection_address
        );
        for (&regular_address, &mmff_s_address) in regular_before
            .storage_addresses
            .iter()
            .zip(&mmff_s_before.storage_addresses)
        {
            assert_ne!(regular_address, mmff_s_address);
        }
        for custom in [&custom_regular_before, &custom_s_before] {
            assert_ne!(custom.collection_address, regular_before.collection_address);
            assert_ne!(custom.collection_address, mmff_s_before.collection_address);
            for &custom_address in &custom.storage_addresses {
                for &default_address in regular_before
                    .storage_addresses
                    .iter()
                    .chain(&mmff_s_before.storage_addresses)
                {
                    assert_ne!(custom_address, default_address);
                }
            }
        }
        assert_eq!(DEFAULT_MMFF_OOP_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(DEFAULT_MMFF_OOP_S_CONSTRUCTIONS.load(Ordering::Relaxed), 1);

        let post_default_custom_regular =
            mmff_oop_constructor_with_input_checkpoint(false, CUSTOM_SOURCE)
                .expect("regular custom construction after defaults remains independent");
        let post_default_custom_regular_before =
            mmff_oop_default_snapshot(&post_default_custom_regular);
        let post_default_custom_s = mmff_oop_constructor_with_input_checkpoint(true, CUSTOM_SOURCE)
            .expect("MMFFs custom construction after defaults remains independent");
        let post_default_custom_s_before = mmff_oop_default_snapshot(&post_default_custom_s);
        for post_custom in [
            &post_default_custom_regular_before,
            &post_default_custom_s_before,
        ] {
            assert_ne!(
                post_custom.collection_address,
                regular_before.collection_address
            );
            assert_ne!(
                post_custom.collection_address,
                mmff_s_before.collection_address
            );
        }
        assert_eq!(
            post_default_custom_regular_before.key_rows,
            custom_regular_before.key_rows
        );
        assert_eq!(
            post_default_custom_regular_before.koop_bits,
            custom_regular_before.koop_bits
        );
        assert_eq!(
            post_default_custom_s_before.key_rows,
            custom_s_before.key_rows
        );
        assert_eq!(
            post_default_custom_s_before.koop_bits,
            custom_s_before.koop_bits
        );
        for post_custom in [
            &post_default_custom_regular_before,
            &post_default_custom_s_before,
        ] {
            for &post_custom_address in &post_custom.storage_addresses {
                for &default_address in regular_before
                    .storage_addresses
                    .iter()
                    .chain(&mmff_s_before.storage_addresses)
                {
                    assert_ne!(post_custom_address, default_address);
                }
            }
        }
        assert_eq!(mmff_oop_default_snapshot(regular), regular_before);
        assert_eq!(mmff_oop_default_snapshot(mmff_s), mmff_s_before);
        assert_eq!(DEFAULT_MMFF_OOP_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(DEFAULT_MMFF_OOP_S_CONSTRUCTIONS.load(Ordering::Relaxed), 1);

        let (regular_warm_calls, mmff_s_warm_calls) = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(8);
            for _ in 0..8 {
                let regular_expected = &regular_before;
                let mmff_s_expected = &mmff_s_before;
                handles.push(scope.spawn(move || {
                    let mut regular_calls = 0;
                    let mut mmff_s_calls = 0;
                    for call_index in 0..16 {
                        let is_mmff_s = call_index % 2 == 1;
                        let regular_before_call = mmff_oop_default_snapshot(regular);
                        let mmff_s_before_call = mmff_oop_default_snapshot(mmff_s);
                        assert_eq!(&regular_before_call, regular_expected);
                        assert_eq!(&mmff_s_before_call, mmff_s_expected);
                        assert_eq!(DEFAULT_MMFF_OOP_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
                        assert_eq!(DEFAULT_MMFF_OOP_S_CONSTRUCTIONS.load(Ordering::Relaxed), 1);

                        let cached = default_mmff_oop(is_mmff_s)
                            .expect("warm OOP default retains its cached parse result");
                        if is_mmff_s {
                            mmff_s_calls += 1;
                            assert!(std::ptr::eq(cached, mmff_s));
                        } else {
                            regular_calls += 1;
                            assert!(std::ptr::eq(cached, regular));
                        }

                        let regular_after_call = mmff_oop_default_snapshot(regular);
                        let mmff_s_after_call = mmff_oop_default_snapshot(mmff_s);
                        assert_eq!(regular_after_call, regular_before_call);
                        assert_eq!(mmff_s_after_call, mmff_s_before_call);
                        assert_eq!(&regular_after_call, regular_expected);
                        assert_eq!(&mmff_s_after_call, mmff_s_expected);
                        assert_eq!(DEFAULT_MMFF_OOP_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
                        assert_eq!(DEFAULT_MMFF_OOP_S_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
                    }
                    (regular_calls, mmff_s_calls)
                }));
            }

            handles
                .into_iter()
                .fold((0, 0), |(regular_total, s_total), handle| {
                    let (regular_calls, s_calls) =
                        handle.join().expect("scoped OOP default worker");
                    (regular_total + regular_calls, s_total + s_calls)
                })
        });
        assert_eq!(regular_warm_calls, 64);
        assert_eq!(mmff_s_warm_calls, 64);
        assert_eq!(regular_warm_calls + mmff_s_warm_calls, 128);
        assert_eq!(mmff_oop_default_snapshot(regular), regular_before);
        assert_eq!(mmff_oop_default_snapshot(mmff_s), mmff_s_before);
        assert_eq!(DEFAULT_MMFF_OOP_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
        assert_eq!(DEFAULT_MMFF_OOP_S_CONSTRUCTIONS.load(Ordering::Relaxed), 1);
    }

    #[derive(Debug, PartialEq, Eq)]
    struct MmffTorDefaultSnapshot {
        collection_address: usize,
        key_rows: [Vec<u8>; 5],
        value_bits: Vec<[u64; 3]>,
        storage_addresses: [usize; 6],
        sentinel_rows: Vec<(usize, usize, [u8; 5], [u64; 3])>,
    }

    fn mmff_tor_default_snapshot(collection: &MmffTorCollection) -> MmffTorDefaultSnapshot {
        const SENTINEL_INDICES: [usize; 4] = [0, 22, 463, 925];
        let sentinel_rows = SENTINEL_INDICES
            .into_iter()
            .filter(|&index| index < collection.d_params.len())
            .map(|index| {
                let row = &collection.d_params[index];
                (
                    index,
                    row as *const MmffTor as usize,
                    [
                        collection.d_tor_type[index],
                        collection.d_i_atom_type[index],
                        collection.d_j_atom_type[index],
                        collection.d_k_atom_type[index],
                        collection.d_l_atom_type[index],
                    ],
                    [row.v1.to_bits(), row.v2.to_bits(), row.v3.to_bits()],
                )
            })
            .collect();

        MmffTorDefaultSnapshot {
            collection_address: collection as *const MmffTorCollection as usize,
            key_rows: [
                collection.d_tor_type.clone(),
                collection.d_i_atom_type.clone(),
                collection.d_j_atom_type.clone(),
                collection.d_k_atom_type.clone(),
                collection.d_l_atom_type.clone(),
            ],
            value_bits: collection
                .d_params
                .iter()
                .map(|row| [row.v1.to_bits(), row.v2.to_bits(), row.v3.to_bits()])
                .collect(),
            storage_addresses: [
                collection.d_tor_type.as_ptr() as usize,
                collection.d_i_atom_type.as_ptr() as usize,
                collection.d_j_atom_type.as_ptr() as usize,
                collection.d_k_atom_type.as_ptr() as usize,
                collection.d_l_atom_type.as_ptr() as usize,
                collection.d_params.as_ptr() as usize,
            ],
            sentinel_rows,
        }
    }

    fn mmff_tor_default_construction_counts() -> (usize, usize) {
        (
            DEFAULT_MMFF_TOR_CONSTRUCTIONS.load(Ordering::Relaxed),
            DEFAULT_MMFF_TOR_S_CONSTRUCTIONS.load(Ordering::Relaxed),
        )
    }

    #[test]
    fn mmff_tor_defaults_assets_sentinels_and_128_warm_borrows() {
        const REGULAR_ASSET: &str = include_str!("default_tor.tsv");
        const MMFF_S_ASSET: &str = include_str!("default_tor_s.tsv");
        assert_eq!(REGULAR_ASSET.len(), 39_795);
        assert_eq!(MMFF_S_ASSET.len(), 39_888);
        assert!(REGULAR_ASSET.ends_with('\n'));
        assert!(MMFF_S_ASSET.ends_with('\n'));
        assert_eq!(
            sha256_for_fixed_asset_test(REGULAR_ASSET.as_bytes()),
            [
                0xfe, 0xc5, 0x3b, 0x92, 0x8f, 0xee, 0x6c, 0xf4, 0x89, 0xec, 0x6d, 0xd3, 0xe7, 0x36,
                0x97, 0xf6, 0x26, 0xb9, 0x44, 0xef, 0x4c, 0x9d, 0x76, 0xd8, 0x84, 0x13, 0x97, 0xec,
                0x8c, 0x09, 0x17, 0xa6,
            ]
        );
        assert_eq!(
            sha256_for_fixed_asset_test(MMFF_S_ASSET.as_bytes()),
            [
                0xb5, 0xf2, 0xd7, 0x58, 0x89, 0x30, 0xee, 0x19, 0x37, 0x12, 0xc3, 0x50, 0x12, 0xe5,
                0x4c, 0xfd, 0x83, 0xb5, 0xd3, 0x27, 0xb3, 0xa3, 0xf5, 0x26, 0x04, 0xcc, 0x26, 0x55,
                0x6a, 0x37, 0x57, 0x6a,
            ]
        );
        assert_eq!(
            REGULAR_ASSET
                .lines()
                .filter(|line| !line.starts_with('*'))
                .count(),
            926
        );
        assert_eq!(
            MMFF_S_ASSET
                .lines()
                .filter(|line| !line.starts_with('*'))
                .count(),
            926
        );

        assert_eq!(mmff_tor_default_construction_counts(), (0, 0));

        const CUSTOM_SOURCE: &str = "*fixed\n5\t10\t7\t8\t20\t1.25\t-0.0\t2.5\n";
        let custom_regular = mmff_tor_constructor_with_input_checkpoint(false, CUSTOM_SOURCE)
            .expect("regular custom Tor construction is independent of cached defaults");
        let custom_regular_before = mmff_tor_default_snapshot(&custom_regular);
        let custom_s = mmff_tor_constructor_with_input_checkpoint(true, CUSTOM_SOURCE)
            .expect("MMFFs custom Tor construction is independent of cached defaults");
        let custom_s_before = mmff_tor_default_snapshot(&custom_s);
        for custom in [&custom_regular, &custom_s] {
            assert_eq!(custom.d_tor_type.as_slice(), &[5]);
            assert_eq!(custom.d_i_atom_type.as_slice(), &[10]);
            assert_eq!(custom.d_j_atom_type.as_slice(), &[7]);
            assert_eq!(custom.d_k_atom_type.as_slice(), &[8]);
            assert_eq!(custom.d_l_atom_type.as_slice(), &[20]);
            assert_eq!(
                [
                    custom.d_params[0].v1.to_bits(),
                    custom.d_params[0].v2.to_bits(),
                    custom.d_params[0].v3.to_bits(),
                ],
                [
                    0x3ff4_0000_0000_0000,
                    0x8000_0000_0000_0000,
                    0x4004_0000_0000_0000,
                ]
            );
        }

        let empty_custom_regular = mmff_tor_constructor_with_input_checkpoint(false, "")
            .expect("empty custom Tor constructor selects the regular asset");
        let empty_custom_regular_before = mmff_tor_default_snapshot(&empty_custom_regular);
        let empty_custom_s = mmff_tor_constructor_with_input_checkpoint(true, "")
            .expect("empty custom Tor constructor selects the MMFFs asset");
        let empty_custom_s_before = mmff_tor_default_snapshot(&empty_custom_s);
        for empty_custom in [&empty_custom_regular, &empty_custom_s] {
            assert_eq!(empty_custom.d_tor_type.len(), 926);
            assert_eq!(empty_custom.d_i_atom_type.len(), 926);
            assert_eq!(empty_custom.d_j_atom_type.len(), 926);
            assert_eq!(empty_custom.d_k_atom_type.len(), 926);
            assert_eq!(empty_custom.d_l_atom_type.len(), 926);
            assert_eq!(empty_custom.d_params.len(), 926);
        }
        assert_eq!(mmff_tor_default_construction_counts(), (0, 0));

        let regular_result = default_mmff_tor(false);
        assert_eq!(mmff_tor_default_construction_counts(), (1, 1));
        let regular = regular_result.expect("regular Tor default asset parses");
        let regular_before = mmff_tor_default_snapshot(regular);

        let regular_before_s_getter = mmff_tor_default_snapshot(regular);
        let s_result = default_mmff_tor(true);
        let regular_after_s_getter = mmff_tor_default_snapshot(regular);
        assert_eq!(regular_after_s_getter, regular_before_s_getter);
        assert_eq!(mmff_tor_default_construction_counts(), (1, 1));
        let mmff_s = s_result.expect("MMFFs Tor default asset parses");
        let mmff_s_before = mmff_tor_default_snapshot(mmff_s);

        for defaults in [regular, mmff_s] {
            assert_eq!(defaults.d_tor_type.len(), 926);
            assert_eq!(defaults.d_i_atom_type.len(), 926);
            assert_eq!(defaults.d_j_atom_type.len(), 926);
            assert_eq!(defaults.d_k_atom_type.len(), 926);
            assert_eq!(defaults.d_l_atom_type.len(), 926);
            assert_eq!(defaults.d_params.len(), 926);
        }
        assert_eq!(
            regular_before.key_rows,
            empty_custom_regular_before.key_rows
        );
        assert_eq!(
            regular_before.value_bits,
            empty_custom_regular_before.value_bits
        );
        assert_eq!(mmff_s_before.key_rows, empty_custom_s_before.key_rows);
        assert_eq!(mmff_s_before.value_bits, empty_custom_s_before.value_bits);

        const SENTINELS: [(usize, [u8; 5], [u64; 3], [u64; 3]); 4] = [
            (
                0,
                [0, 0, 1, 1, 0],
                [0, 0, 0x3fd3_3333_3333_3333],
                [0, 0, 0x3fd3_3333_3333_3333],
            ),
            (
                22,
                [0, 5, 1, 1, 10],
                [0, 0, 0x3fdb_53f7_ced9_1687],
                [0, 0, 0x3fda_c083_126e_978d],
            ),
            (
                463,
                [0, 0, 3, 54, 0],
                [0, 0x4020_0000_0000_0000, 0],
                [0, 0x4020_0000_0000_0000, 0],
            ),
            (
                925,
                [0, 0, 80, 81, 0],
                [0, 0x4010_0000_0000_0000, 0],
                [0, 0x4010_0000_0000_0000, 0],
            ),
        ];
        assert_eq!(regular_before.sentinel_rows.len(), SENTINELS.len());
        assert_eq!(mmff_s_before.sentinel_rows.len(), SENTINELS.len());
        for (sentinel_index, (index, keys, regular_bits, s_bits)) in
            SENTINELS.into_iter().enumerate()
        {
            let regular_row = regular_before.sentinel_rows[sentinel_index];
            let s_row = mmff_s_before.sentinel_rows[sentinel_index];
            assert_eq!(regular_row.0, index);
            assert_eq!(regular_row.2, keys);
            assert_eq!(regular_row.3, regular_bits);
            assert_eq!(s_row.0, index);
            assert_eq!(s_row.2, keys);
            assert_eq!(s_row.3, s_bits);
            assert_eq!(
                [
                    regular.d_params[index].v1.to_bits(),
                    regular.d_params[index].v2.to_bits(),
                    regular.d_params[index].v3.to_bits(),
                ],
                regular_bits,
                "regular Tor sentinel row {index}"
            );
            assert_eq!(
                [
                    mmff_s.d_params[index].v1.to_bits(),
                    mmff_s.d_params[index].v2.to_bits(),
                    mmff_s.d_params[index].v3.to_bits(),
                ],
                s_bits,
                "MMFFs Tor sentinel row {index}"
            );
            assert!(std::ptr::eq(
                &regular.d_params[index],
                &regular.d_params[regular_row.0]
            ));
            assert!(std::ptr::eq(
                &mmff_s.d_params[index],
                &mmff_s.d_params[s_row.0]
            ));
        }

        assert_ne!(
            regular_before.collection_address,
            mmff_s_before.collection_address
        );
        for (&regular_address, &s_address) in regular_before
            .storage_addresses
            .iter()
            .zip(&mmff_s_before.storage_addresses)
        {
            assert_ne!(regular_address, s_address);
        }
        for custom in [
            &custom_regular_before,
            &custom_s_before,
            &empty_custom_regular_before,
            &empty_custom_s_before,
        ] {
            assert_ne!(custom.collection_address, regular_before.collection_address);
            assert_ne!(custom.collection_address, mmff_s_before.collection_address);
            for &custom_address in &custom.storage_addresses {
                for &default_address in regular_before
                    .storage_addresses
                    .iter()
                    .chain(&mmff_s_before.storage_addresses)
                {
                    assert_ne!(custom_address, default_address);
                }
            }
        }

        let regular_after_custom_before = mmff_tor_default_snapshot(regular);
        let s_after_custom_before = mmff_tor_default_snapshot(mmff_s);
        let post_default_custom_regular =
            mmff_tor_constructor_with_input_checkpoint(false, CUSTOM_SOURCE)
                .expect("regular custom Tor construction after defaults remains independent");
        let post_default_custom_regular_after =
            mmff_tor_default_snapshot(&post_default_custom_regular);
        let post_default_custom_s = mmff_tor_constructor_with_input_checkpoint(true, CUSTOM_SOURCE)
            .expect("MMFFs custom Tor construction after defaults remains independent");
        let post_default_custom_s_after = mmff_tor_default_snapshot(&post_default_custom_s);
        assert_eq!(
            post_default_custom_regular_after.key_rows,
            custom_regular_before.key_rows
        );
        assert_eq!(
            post_default_custom_regular_after.value_bits,
            custom_regular_before.value_bits
        );
        assert_eq!(
            post_default_custom_s_after.key_rows,
            custom_s_before.key_rows
        );
        assert_eq!(
            post_default_custom_s_after.value_bits,
            custom_s_before.value_bits
        );
        assert_eq!(
            mmff_tor_default_snapshot(regular),
            regular_after_custom_before
        );
        assert_eq!(mmff_tor_default_snapshot(mmff_s), s_after_custom_before);
        assert_eq!(mmff_tor_default_construction_counts(), (1, 1));
        for custom in [
            &post_default_custom_regular_after,
            &post_default_custom_s_after,
        ] {
            assert_ne!(custom.collection_address, regular_before.collection_address);
            assert_ne!(custom.collection_address, mmff_s_before.collection_address);
            for &custom_address in &custom.storage_addresses {
                for &default_address in regular_before
                    .storage_addresses
                    .iter()
                    .chain(&mmff_s_before.storage_addresses)
                {
                    assert_ne!(custom_address, default_address);
                }
            }
        }

        let (regular_warm_calls, s_warm_calls) = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(8);
            for _ in 0..8 {
                let regular_expected = &regular_before;
                let s_expected = &mmff_s_before;
                handles.push(scope.spawn(move || {
                    let mut regular_calls = 0;
                    let mut s_calls = 0;
                    for call_index in 0..16 {
                        let is_mmff_s = call_index % 2 == 1;
                        let regular_before_call = mmff_tor_default_snapshot(regular);
                        let s_before_call = mmff_tor_default_snapshot(mmff_s);
                        let counts_before_call = mmff_tor_default_construction_counts();
                        assert_eq!(&regular_before_call, regular_expected);
                        assert_eq!(&s_before_call, s_expected);
                        assert_eq!(counts_before_call, (1, 1));

                        let result = default_mmff_tor(is_mmff_s);
                        let regular_after_call = mmff_tor_default_snapshot(regular);
                        let s_after_call = mmff_tor_default_snapshot(mmff_s);
                        let counts_after_call = mmff_tor_default_construction_counts();
                        assert_eq!(regular_after_call, regular_before_call);
                        assert_eq!(s_after_call, s_before_call);
                        assert_eq!(regular_after_call, *regular_expected);
                        assert_eq!(s_after_call, *s_expected);
                        assert_eq!(counts_after_call, counts_before_call);
                        assert_eq!(counts_after_call, (1, 1));

                        let cached = result.expect("warm Tor getter retains cached parse outcome");
                        if is_mmff_s {
                            s_calls += 1;
                            assert!(std::ptr::eq(cached, mmff_s));
                        } else {
                            regular_calls += 1;
                            assert!(std::ptr::eq(cached, regular));
                        }
                    }
                    (regular_calls, s_calls)
                }));
            }

            handles
                .into_iter()
                .fold((0, 0), |(regular_total, s_total), handle| {
                    let (regular_calls, s_calls) =
                        handle.join().expect("scoped Tor default worker");
                    (regular_total + regular_calls, s_total + s_calls)
                })
        });
        assert_eq!(regular_warm_calls, 64);
        assert_eq!(s_warm_calls, 64);
        assert_eq!(regular_warm_calls + s_warm_calls, 128);
        assert_eq!(mmff_tor_default_snapshot(regular), regular_before);
        assert_eq!(mmff_tor_default_snapshot(mmff_s), mmff_s_before);
        assert_eq!(mmff_tor_default_construction_counts(), (1, 1));
    }

    fn mmff_vdw_constructor_with_input_checkpoint(
        source: &str,
        actual_calls: &mut usize,
    ) -> Result<MmffVdwCollection, MmffVdwParseError> {
        let input_address = source.as_ptr();
        let input_bytes = source.as_bytes().to_vec();
        let result = MmffVdwCollection::from_text(source);
        assert_eq!(source.as_ptr(), input_address);
        assert_eq!(source.as_bytes(), input_bytes.as_slice());
        *actual_calls += 1;
        result
    }

    fn assert_mmff_vdw_literal_collection(
        collection: &MmffVdwCollection,
        header_bits: [u64; 5],
        keys: &[u8],
        rows: &[[u64; 5]],
        da_values: &[u8],
    ) {
        assert_eq!(
            [
                collection.power.to_bits(),
                collection.b.to_bits(),
                collection.beta.to_bits(),
                collection.darad.to_bits(),
                collection.daeps.to_bits(),
            ],
            header_bits
        );
        assert_eq!(collection.d_atom_type.as_slice(), keys);
        assert_eq!(collection.d_params.len(), rows.len());
        assert_eq!(da_values.len(), rows.len());
        for (row_index, ((row, expected_bits), &expected_da)) in collection
            .d_params
            .iter()
            .zip(rows)
            .zip(da_values)
            .enumerate()
        {
            assert_eq!(
                [
                    row.alpha_i.to_bits(),
                    row.n_i.to_bits(),
                    row.a_i.to_bits(),
                    row.g_i.to_bits(),
                    row.r_star.to_bits(),
                ],
                *expected_bits,
                "VdW row {row_index}"
            );
            assert_eq!(row.da, expected_da, "VdW row {row_index} DA");
        }
    }

    fn replace_mmff_vdw_field(
        source: &str,
        physical_line_index: usize,
        column: usize,
        replacement: &str,
    ) -> String {
        let mut lines = source.split('\n').map(str::to_owned).collect::<Vec<_>>();
        let mut fields = lines[physical_line_index]
            .split('\t')
            .map(str::to_owned)
            .collect::<Vec<_>>();
        fields[column] = replacement.to_owned();
        lines[physical_line_index] = fields.join("\t");
        lines.join("\n")
    }

    fn assert_mmff_vdw_error_source(error: &MmffVdwParseError) {
        match error {
            MmffVdwParseError::Table(original) => {
                let expected: &(dyn std::error::Error + 'static) = original;
                let source = std::error::Error::source(error)
                    .expect("table error preserves its original parse cause");
                assert!(std::ptr::eq(source, expected));
            }
            MmffVdwParseError::MissingHeader => {
                assert!(std::error::Error::source(error).is_none());
            }
        }
    }

    #[derive(Debug, PartialEq, Eq)]
    struct MmffVdwLookupSnapshot {
        collection_address: usize,
        header_bits: [u64; 5],
        keys: Vec<u8>,
        row_bits: Vec<[u64; 5]>,
        da_values: Vec<u8>,
        key_buffer_address: usize,
        row_buffer_address: usize,
    }

    fn mmff_vdw_lookup_snapshot(collection: &MmffVdwCollection) -> MmffVdwLookupSnapshot {
        MmffVdwLookupSnapshot {
            collection_address: collection as *const MmffVdwCollection as usize,
            header_bits: [
                collection.power.to_bits(),
                collection.b.to_bits(),
                collection.beta.to_bits(),
                collection.darad.to_bits(),
                collection.daeps.to_bits(),
            ],
            keys: collection.d_atom_type.clone(),
            row_bits: collection
                .d_params
                .iter()
                .map(|row| {
                    [
                        row.alpha_i.to_bits(),
                        row.n_i.to_bits(),
                        row.a_i.to_bits(),
                        row.g_i.to_bits(),
                        row.r_star.to_bits(),
                    ]
                })
                .collect(),
            da_values: collection.d_params.iter().map(|row| row.da).collect(),
            key_buffer_address: collection.d_atom_type.as_ptr() as usize,
            row_buffer_address: collection.d_params.as_ptr() as usize,
        }
    }

    #[test]
    fn mmff_vdw_constructor_38_frozen_calls_preserve_inputs_and_outputs() {
        const BASE: &str = concat!(
            "*fixed\n",
            "0.5\t-0.0\t12\t0.8\t0.5\n",
            "257\t4\t2\t3\t-0.0\t-ignored\n",
            "257\t9\t3\t2\t4\tAB\n",
        );
        const CUSTOM_HEADER_BITS: [u64; 5] = [
            0x3fe0_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x4028_0000_0000_0000,
            0x3fe9_9999_9999_999a,
            0x3fe0_0000_0000_0000,
        ];
        const CUSTOM_KEYS: [u8; 2] = [1, 1];
        const CUSTOM_ROWS: [[u64; 5]; 2] = [
            [
                0x4010_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x8000_0000_0000_0000,
                0x4018_0000_0000_0000,
            ],
            [
                0x4022_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x4018_0000_0000_0000,
            ],
        ];
        const CUSTOM_DA: [u8; 2] = [b'-', b'A'];
        const DEFAULT_HEADER_BITS: [u64; 5] = [
            0x3fd0_0000_0000_0000,
            0x3fc9_9999_9999_999a,
            0x4028_0000_0000_0000,
            0x3fe9_9999_9999_999a,
            0x3fe0_0000_0000_0000,
        ];
        const DEFAULT_SENTINELS: [(usize, u8, [u64; 5], u8); 4] = [
            (
                0,
                1,
                [
                    0x3ff0_cccc_cccc_cccd,
                    0x4003_eb85_1eb8_51ec,
                    0x400f_1eb8_51eb_851f,
                    0x3ff4_8312_6e97_8d50,
                    0x400f_807d_4cf4_26b8,
                ],
                b'-',
            ),
            (
                20,
                21,
                [
                    0x3fc3_3333_3333_3333,
                    0x3fe9_9999_9999_999a,
                    0x4010_cccc_cccc_cccd,
                    0x3ff3_5810_624d_d2f2,
                    0x4004_e90f_30bd_1fe1,
                ],
                b'D',
            ),
            (
                81,
                82,
                [
                    0x3fee_6666_6666_6666,
                    0x4006_8f5c_28f5_c28f,
                    0x400f_1eb8_51eb_851f,
                    0x3ff4_8312_6e97_8d50,
                    0x400e_b936_5f82_0f67,
                ],
                b'A',
            ),
            (
                94,
                99,
                [
                    0x3fd6_6666_6666_6666,
                    0x400c_0000_0000_0000,
                    0x4010_0000_0000_0000,
                    0x3ff4_cccc_cccc_cccd,
                    0x4008_9cf6_9f3f_7dd8,
                ],
                b'-',
            ),
        ];

        let mut suffix_lines = BASE.split('\n').map(str::to_owned).collect::<Vec<_>>();
        suffix_lines[1].push_str("\tunused-header-suffix");
        suffix_lines[2].push_str("\tunused-row-suffix");
        suffix_lines[3].push_str("\tunused-row-suffix");
        let formatted_base = suffix_lines.join("\n");
        let crlf = formatted_base.replace('\n', "\r\n");
        let doubled_tabs = formatted_base.replace('\t', "\t\t");
        let unterminated = formatted_base
            .strip_suffix('\n')
            .expect("the fixed table ends with LF")
            .to_owned();
        let format_inputs = [
            ("LF", formatted_base.as_str(), 2),
            ("CRLF", crlf.as_str(), 2),
            ("doubled TAB", doubled_tabs.as_str(), 2),
            ("unterminated final row", unterminated.as_str(), 1),
        ];
        let mut actual_calls = 0;
        for (format, source, expected_rows) in format_inputs {
            for repeat in 0..2 {
                let collection =
                    mmff_vdw_constructor_with_input_checkpoint(source, &mut actual_calls)
                        .expect("the frozen custom VdW format parses");
                assert_mmff_vdw_literal_collection(
                    &collection,
                    CUSTOM_HEADER_BITS,
                    &CUSTOM_KEYS[..expected_rows],
                    &CUSTOM_ROWS[..expected_rows],
                    &CUSTOM_DA[..expected_rows],
                );
                assert_eq!(
                    collection.d_params.len(),
                    expected_rows,
                    "{format} repeat {repeat}"
                );
            }
        }

        for column in 0..5 {
            let source = replace_mmff_vdw_field(BASE, 1, column, "X");
            for repeat in 0..2 {
                let error = mmff_vdw_constructor_with_input_checkpoint(&source, &mut actual_calls)
                    .expect_err("an invalid VdW header field is rejected");
                assert_eq!(
                    error,
                    MmffVdwParseError::Table(MmffParamParseError {
                        table: MmffParamTable::Vdw,
                        line: 2,
                        column,
                        cause: MmffParamParseCause::InvalidFloat {
                            cell: "X".to_owned(),
                        },
                    })
                );
                assert_mmff_vdw_error_source(&error);
            }
        }

        for column in 0..5 {
            let source = replace_mmff_vdw_field(BASE, 2, column, "X");
            for repeat in 0..2 {
                let cause = if column == 0 {
                    MmffParamParseCause::InvalidUnsigned {
                        cell: "X".to_owned(),
                    }
                } else {
                    MmffParamParseCause::InvalidFloat {
                        cell: "X".to_owned(),
                    }
                };
                let error = mmff_vdw_constructor_with_input_checkpoint(&source, &mut actual_calls)
                    .expect_err("an invalid VdW row field is rejected");
                assert_eq!(
                    error,
                    MmffVdwParseError::Table(MmffParamParseError {
                        table: MmffParamTable::Vdw,
                        line: 3,
                        column,
                        cause,
                    })
                );
                assert_mmff_vdw_error_source(&error);
            }
        }

        let mut no_da_lines = BASE.split('\n').map(str::to_owned).collect::<Vec<_>>();
        let mut first_row = no_da_lines[2]
            .split('\t')
            .map(str::to_owned)
            .collect::<Vec<_>>();
        first_row.pop();
        no_da_lines[2] = first_row.join("\t");
        let missing_da = no_da_lines.join("\n");
        for _ in 0..2 {
            let error = mmff_vdw_constructor_with_input_checkpoint(&missing_da, &mut actual_calls)
                .expect_err("a missing VdW DA token is rejected");
            assert_eq!(
                error,
                MmffVdwParseError::Table(MmffParamParseError {
                    table: MmffParamTable::Vdw,
                    line: 3,
                    column: 5,
                    cause: MmffParamParseCause::MissingToken,
                })
            );
            assert_mmff_vdw_error_source(&error);
        }

        for _ in 0..2 {
            let error = mmff_vdw_constructor_with_input_checkpoint("*fixed\n\n", &mut actual_calls)
                .expect_err("an empty processed VdW line is rejected");
            assert_eq!(
                error,
                MmffVdwParseError::Table(MmffParamParseError {
                    table: MmffParamTable::Vdw,
                    line: 2,
                    column: 0,
                    cause: MmffParamParseCause::EmptyProcessedLine,
                })
            );
            assert_mmff_vdw_error_source(&error);
        }

        for _ in 0..2 {
            let error = mmff_vdw_constructor_with_input_checkpoint("*fixed\n", &mut actual_calls)
                .expect_err("a comment-only VdW source has no processed header");
            assert_eq!(error, MmffVdwParseError::MissingHeader);
            assert_eq!(
                error.to_string(),
                "MMFF VdW source has no processed constants record"
            );
            assert_mmff_vdw_error_source(&error);
        }

        for _ in 0..2 {
            let collection = mmff_vdw_constructor_with_input_checkpoint(
                "*fixed\n0.5\t-0.0\t12\t0.8\t0.5\n",
                &mut actual_calls,
            )
            .expect("a terminated header with no rows is a valid empty collection");
            assert_mmff_vdw_literal_collection(&collection, CUSTOM_HEADER_BITS, &[], &[], &[]);
        }

        const DEFAULT_ASSET: &str = include_str!("default_vdw.tsv");
        assert_eq!(DEFAULT_ASSET.len(), 4073);
        assert_eq!(
            sha256_for_fixed_asset_test(DEFAULT_ASSET.as_bytes()),
            [
                0xc9, 0xd6, 0x0e, 0x3d, 0xbb, 0x91, 0x58, 0x7a, 0x02, 0x31, 0x40, 0x55, 0xf6, 0x4b,
                0x03, 0x16, 0x96, 0x67, 0xd0, 0x42, 0xea, 0xbc, 0xe9, 0xd6, 0xe0, 0x33, 0xa4, 0xf0,
                0xb8, 0x7a, 0xa7, 0xca,
            ]
        );
        for _ in 0..2 {
            let collection = mmff_vdw_constructor_with_input_checkpoint("", &mut actual_calls)
                .expect("empty input selects the frozen default VdW asset");
            assert_eq!(collection.d_params.len(), 95);
            assert_eq!(collection.d_atom_type.len(), 95);
            assert_eq!(
                [
                    collection.power.to_bits(),
                    collection.b.to_bits(),
                    collection.beta.to_bits(),
                    collection.darad.to_bits(),
                    collection.daeps.to_bits(),
                ],
                DEFAULT_HEADER_BITS
            );
            for (index, atom_type, expected_bits, expected_da) in DEFAULT_SENTINELS {
                assert_eq!(collection.d_atom_type[index], atom_type);
                let row = &collection.d_params[index];
                assert_eq!(
                    [
                        row.alpha_i.to_bits(),
                        row.n_i.to_bits(),
                        row.a_i.to_bits(),
                        row.g_i.to_bits(),
                        row.r_star.to_bits(),
                    ],
                    expected_bits,
                    "default VdW source row {atom_type}"
                );
                assert_eq!(row.da, expected_da, "default VdW source row {atom_type}");
            }
        }

        assert_eq!(actual_calls, 38);
    }

    #[test]
    fn mmff_vdw_power_24_frozen_calls_match_source_literal_bits() {
        const CELLS: [(&str, &str, u64, u64, u64); 12] = [
            (
                "1.0",
                "0.0",
                0x3ff0_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x4000_0000_0000_0000,
            ),
            (
                "1.0",
                "0.25",
                0x3ff0_0000_0000_0000,
                0x3fd0_0000_0000_0000,
                0x4000_0000_0000_0000,
            ),
            (
                "1.0",
                "0.5",
                0x3ff0_0000_0000_0000,
                0x3fe0_0000_0000_0000,
                0x4000_0000_0000_0000,
            ),
            (
                "1.0",
                "1.0",
                0x3ff0_0000_0000_0000,
                0x3ff0_0000_0000_0000,
                0x4000_0000_0000_0000,
            ),
            (
                "16.0",
                "0.0",
                0x4030_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x4000_0000_0000_0000,
            ),
            (
                "16.0",
                "0.25",
                0x4030_0000_0000_0000,
                0x3fd0_0000_0000_0000,
                0x4010_0000_0000_0000,
            ),
            (
                "16.0",
                "0.5",
                0x4030_0000_0000_0000,
                0x3fe0_0000_0000_0000,
                0x4020_0000_0000_0000,
            ),
            (
                "16.0",
                "1.0",
                0x4030_0000_0000_0000,
                0x3ff0_0000_0000_0000,
                0x4040_0000_0000_0000,
            ),
            (
                "81.0",
                "0.0",
                0x4054_4000_0000_0000,
                0x0000_0000_0000_0000,
                0x4000_0000_0000_0000,
            ),
            (
                "81.0",
                "0.25",
                0x4054_4000_0000_0000,
                0x3fd0_0000_0000_0000,
                0x4018_0000_0000_0000,
            ),
            (
                "81.0",
                "0.5",
                0x4054_4000_0000_0000,
                0x3fe0_0000_0000_0000,
                0x4032_0000_0000_0000,
            ),
            (
                "81.0",
                "1.0",
                0x4054_4000_0000_0000,
                0x3ff0_0000_0000_0000,
                0x4064_4000_0000_0000,
            ),
        ];
        const NEGATIVE_ZERO_BITS: u64 = 0x8000_0000_0000_0000;
        const TWO_BITS: u64 = 0x4000_0000_0000_0000;
        const HEADER_TAIL_BITS: [u64; 4] = [
            NEGATIVE_ZERO_BITS,
            0x4028_0000_0000_0000,
            0x3fe9_9999_9999_999a,
            0x3fe0_0000_0000_0000,
        ];

        let mut actual_calls = 0;
        for (alpha_text, power_text, alpha_bits, power_bits, r_star_bits) in CELLS {
            let source = format!(
                "*fixed\n{power_text}\t-0.0\t12\t0.8\t0.5\n1\t{alpha_text}\t2\t2\t-0.0\tAB\n"
            );
            let header_bits = [
                power_bits,
                HEADER_TAIL_BITS[0],
                HEADER_TAIL_BITS[1],
                HEADER_TAIL_BITS[2],
                HEADER_TAIL_BITS[3],
            ];
            let row_bits = [[
                alpha_bits,
                TWO_BITS,
                TWO_BITS,
                NEGATIVE_ZERO_BITS,
                r_star_bits,
            ]];
            for repeat in 0..2 {
                let collection =
                    mmff_vdw_constructor_with_input_checkpoint(&source, &mut actual_calls)
                        .expect("the source-frozen power cell parses");
                assert_mmff_vdw_literal_collection(
                    &collection,
                    header_bits,
                    &[1],
                    &row_bits,
                    &[b'A'],
                );
                assert_eq!(collection.d_params.len(), 1, "power repeat {repeat}");
            }
        }

        assert_eq!(actual_calls, 24);
    }

    #[test]
    fn mmff_vdw_lookup_72_frozen_calls_preserve_rows_and_inputs() {
        const BASE: &str = concat!(
            "*fixed\n",
            "0.5\t-0.0\t12\t0.8\t0.5\n",
            "257\t4\t2\t3\t-0.0\t-ignored\n",
            "257\t9\t3\t2\t4\tAB\n",
        );
        const HEADER_BITS: [u64; 5] = [
            0x3fe0_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x4028_0000_0000_0000,
            0x3fe9_9999_9999_999a,
            0x3fe0_0000_0000_0000,
        ];
        const KEYS: [u8; 2] = [1, 1];
        const ROW_BITS: [[u64; 5]; 2] = [
            [
                0x4010_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x8000_0000_0000_0000,
                0x4018_0000_0000_0000,
            ],
            [
                0x4022_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x4018_0000_0000_0000,
            ],
        ];
        const DA_VALUES: [u8; 2] = [b'-', b'A'];
        const QUERIES: [u32; 8] = [0, 1, 2, 255, 256, 257, 258, u32::MAX];
        const HEADER_ONLY_QUERIES: [u32; 4] = [0, 1, 257, u32::MAX];
        const FIRST_ROW_BITS: [u64; 5] = [
            0x4010_0000_0000_0000,
            0x4000_0000_0000_0000,
            0x4008_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x4018_0000_0000_0000,
        ];
        const HEADER_ONLY_SOURCE: &str = "*fixed\n0.5\t-0.0\t12\t0.8\t0.5\n";

        let mut suffix_lines = BASE.split('\n').map(str::to_owned).collect::<Vec<_>>();
        suffix_lines[1].push_str("\tunused-header-suffix");
        suffix_lines[2].push_str("\tunused-row-suffix");
        suffix_lines[3].push_str("\tunused-row-suffix");
        let formatted_base = suffix_lines.join("\n");
        let crlf = formatted_base.replace('\n', "\r\n");
        let doubled_tabs = formatted_base.replace('\t', "\t\t");
        let unterminated = formatted_base
            .strip_suffix('\n')
            .expect("the fixed table ends with LF")
            .to_owned();
        let formats = [
            ("LF", formatted_base.as_str(), 2),
            ("CRLF", crlf.as_str(), 2),
            ("doubled TAB", doubled_tabs.as_str(), 2),
            ("unterminated final row", unterminated.as_str(), 1),
        ];
        let mut discrepancies = Vec::new();
        let mut actual_get_calls = 0;

        for (format, source, row_count) in formats {
            let collection = MmffVdwCollection::from_text(source)
                .expect("the frozen custom VdW lookup table parses");
            for query in QUERIES {
                for repeat in 0..2 {
                    let case = format!("{format} query {query} repeat {repeat}");
                    let query_before = query;
                    let input_address_before = source.as_ptr() as usize;
                    let input_bytes_before = source.as_bytes().to_vec();
                    let before = mmff_vdw_lookup_snapshot(&collection);
                    let result = collection.get(query);
                    actual_get_calls += 1;
                    let after = mmff_vdw_lookup_snapshot(&collection);

                    if before != after {
                        discrepancies.push(format!("{case}: collection snapshot changed"));
                    }
                    if source.as_ptr() as usize != input_address_before {
                        discrepancies.push(format!("{case}: input address changed"));
                    }
                    if source.as_bytes() != input_bytes_before.as_slice() {
                        discrepancies.push(format!("{case}: input bytes changed"));
                    }
                    if query != query_before {
                        discrepancies.push(format!("{case}: query value changed"));
                    }
                    if before.header_bits != HEADER_BITS {
                        discrepancies.push(format!("{case}: header bits differ from literals"));
                    }
                    if before.keys.as_slice() != &KEYS[..row_count] {
                        discrepancies.push(format!("{case}: key vector differs from literals"));
                    }
                    if before.row_bits.as_slice() != &ROW_BITS[..row_count] {
                        discrepancies.push(format!("{case}: row bits differ from literals"));
                    }
                    if before.da_values.as_slice() != &DA_VALUES[..row_count] {
                        discrepancies.push(format!("{case}: DA bytes differ from literals"));
                    }

                    if query == 1 {
                        match result {
                            Some(row) => {
                                if !collection
                                    .d_params
                                    .first()
                                    .is_some_and(|first| std::ptr::eq(row, first))
                                {
                                    discrepancies.push(format!(
                                        "{case}: hit did not return first duplicate row"
                                    ));
                                }
                                let row_bits = [
                                    row.alpha_i.to_bits(),
                                    row.n_i.to_bits(),
                                    row.a_i.to_bits(),
                                    row.g_i.to_bits(),
                                    row.r_star.to_bits(),
                                ];
                                if row_bits != FIRST_ROW_BITS || row.da != b'-' {
                                    discrepancies.push(format!(
                                        "{case}: hit differs from frozen row fields"
                                    ));
                                }
                            }
                            None => discrepancies.push(format!("{case}: expected first-row hit")),
                        }
                    } else if result.is_some() {
                        discrepancies.push(format!("{case}: unexpected hit"));
                    }
                }
            }
        }

        let header_only = MmffVdwCollection::from_text(HEADER_ONLY_SOURCE)
            .expect("a processed header with no rows is a valid collection");
        for query in HEADER_ONLY_QUERIES {
            for repeat in 0..2 {
                let case = format!("header-only query {query} repeat {repeat}");
                let query_before = query;
                let input_address_before = HEADER_ONLY_SOURCE.as_ptr() as usize;
                let input_bytes_before = HEADER_ONLY_SOURCE.as_bytes().to_vec();
                let before = mmff_vdw_lookup_snapshot(&header_only);
                let result = header_only.get(query);
                actual_get_calls += 1;
                let after = mmff_vdw_lookup_snapshot(&header_only);

                if before != after {
                    discrepancies.push(format!("{case}: collection snapshot changed"));
                }
                if HEADER_ONLY_SOURCE.as_ptr() as usize != input_address_before {
                    discrepancies.push(format!("{case}: input address changed"));
                }
                if HEADER_ONLY_SOURCE.as_bytes() != input_bytes_before.as_slice() {
                    discrepancies.push(format!("{case}: input bytes changed"));
                }
                if query != query_before {
                    discrepancies.push(format!("{case}: query value changed"));
                }
                if before.header_bits != HEADER_BITS {
                    discrepancies.push(format!("{case}: header bits differ from literals"));
                }
                if !before.keys.is_empty()
                    || !before.row_bits.is_empty()
                    || !before.da_values.is_empty()
                {
                    discrepancies.push(format!("{case}: header-only table contains row state"));
                }
                if result.is_some() {
                    discrepancies.push(format!("{case}: header-only table unexpectedly hit"));
                }
            }
        }

        if actual_get_calls != 72 {
            discrepancies.push(format!(
                "expected 72 actual get calls, got {actual_get_calls}"
            ));
        }
        assert!(
            discrepancies.is_empty(),
            "MMFF VdW lookup regression discrepancies:\n{}",
            discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_vdw_defaults_assets_sentinels_and_128_warm_borrows() {
        const DEFAULT_ASSET: &str = include_str!("default_vdw.tsv");
        const DEFAULT_HEADER_BITS: [u64; 5] = [
            0x3fd0_0000_0000_0000,
            0x3fc9_9999_9999_999a,
            0x4028_0000_0000_0000,
            0x3fe9_9999_9999_999a,
            0x3fe0_0000_0000_0000,
        ];
        const DEFAULT_SENTINELS: [(usize, u32, [u64; 5], u8); 4] = [
            (
                0,
                1,
                [
                    0x3ff0_cccc_cccc_cccd,
                    0x4003_eb85_1eb8_51ec,
                    0x400f_1eb8_51eb_851f,
                    0x3ff4_8312_6e97_8d50,
                    0x400f_807d_4cf4_26b8,
                ],
                b'-',
            ),
            (
                20,
                21,
                [
                    0x3fc3_3333_3333_3333,
                    0x3fe9_9999_9999_999a,
                    0x4010_cccc_cccc_cccd,
                    0x3ff3_5810_624d_d2f2,
                    0x4004_e90f_30bd_1fe1,
                ],
                b'D',
            ),
            (
                81,
                82,
                [
                    0x3fee_6666_6666_6666,
                    0x4006_8f5c_28f5_c28f,
                    0x400f_1eb8_51eb_851f,
                    0x3ff4_8312_6e97_8d50,
                    0x400e_b936_5f82_0f67,
                ],
                b'A',
            ),
            (
                94,
                99,
                [
                    0x3fd6_6666_6666_6666,
                    0x400c_0000_0000_0000,
                    0x4010_0000_0000_0000,
                    0x3ff4_cccc_cccc_cccd,
                    0x4008_9cf6_9f3f_7dd8,
                ],
                b'-',
            ),
        ];
        const MISSING_SOURCE_TYPES: [u32; 4] = [83, 84, 85, 86];
        const WARM_QUERIES: [u32; 8] = [1, 21, 82, 99, 83, 84, 85, 86];
        const CUSTOM_SOURCE: &str = "*fixed\n0.5\t-0.0\t12\t0.8\t0.5\n1\t4\t2\t3\t-0.0\t-\n";

        let mut discrepancies = Vec::new();
        if DEFAULT_ASSET.len() != 4073 {
            discrepancies.push(format!(
                "default asset length is {}, expected 4073",
                DEFAULT_ASSET.len()
            ));
        }
        if sha256_for_fixed_asset_test(DEFAULT_ASSET.as_bytes())
            != [
                0xc9, 0xd6, 0x0e, 0x3d, 0xbb, 0x91, 0x58, 0x7a, 0x02, 0x31, 0x40, 0x55, 0xf6, 0x4b,
                0x03, 0x16, 0x96, 0x67, 0xd0, 0x42, 0xea, 0xbc, 0xe9, 0xd6, 0xe0, 0x33, 0xa4, 0xf0,
                0xb8, 0x7a, 0xa7, 0xca,
            ]
        {
            discrepancies.push("default asset SHA-256 differs from source literal".to_owned());
        }

        if DEFAULT_MMFF_VDW_CONSTRUCTIONS.load(Ordering::Relaxed) != 0 {
            discrepancies.push("VdW default initializer ran before its owner test".to_owned());
        }
        for (index, source) in [CUSTOM_SOURCE, ""].into_iter().enumerate() {
            let collection = MmffVdwCollection::from_text(source)
                .expect("custom VdW construction succeeds before default access");
            drop(collection);
            let count = DEFAULT_MMFF_VDW_CONSTRUCTIONS.load(Ordering::Relaxed);
            if count != 0 {
                discrepancies.push(format!(
                    "pre-initialization constructor {index} changed count to {count}"
                ));
            }
        }

        let default = default_mmff_vdw().expect("the pinned default VdW asset parses");
        if DEFAULT_MMFF_VDW_CONSTRUCTIONS.load(Ordering::Relaxed) != 1 {
            discrepancies.push("first default getter did not initialize exactly once".to_owned());
        }
        let default_address = default as *const MmffVdwCollection as usize;
        let baseline = mmff_vdw_lookup_snapshot(default);
        if default.d_params.len() != 95 || default.d_atom_type.len() != 95 {
            discrepancies.push(format!(
                "default vector lengths are ({}, {}), expected (95, 95)",
                default.d_params.len(),
                default.d_atom_type.len()
            ));
        }
        if baseline.header_bits != DEFAULT_HEADER_BITS {
            discrepancies.push("default header bits differ from frozen source".to_owned());
        }
        for (index, atom_type, expected_bits, expected_da) in DEFAULT_SENTINELS {
            if default.d_atom_type.get(index) != Some(&(atom_type as u8)) {
                discrepancies.push(format!("default source key at row {index} differs"));
            }
            match default.d_params.get(index) {
                Some(row) => {
                    let row_bits = [
                        row.alpha_i.to_bits(),
                        row.n_i.to_bits(),
                        row.a_i.to_bits(),
                        row.g_i.to_bits(),
                        row.r_star.to_bits(),
                    ];
                    if row_bits != expected_bits || row.da != expected_da {
                        discrepancies.push(format!(
                            "default source row {atom_type} differs from frozen bits/DA"
                        ));
                    }
                    if !default
                        .get(atom_type)
                        .is_some_and(|found| std::ptr::eq(found, row))
                    {
                        discrepancies.push(format!(
                            "default source query {atom_type} did not borrow row {index}"
                        ));
                    }
                }
                None => discrepancies.push(format!("default source row {index} is absent")),
            }
        }
        for atom_type in MISSING_SOURCE_TYPES {
            if default.get(atom_type).is_some() {
                discrepancies.push(format!("source type {atom_type} unexpectedly exists"));
            }
        }

        for (index, source) in [CUSTOM_SOURCE, ""].into_iter().enumerate() {
            let collection = MmffVdwCollection::from_text(source)
                .expect("custom VdW construction succeeds after default access");
            drop(collection);
            let count = DEFAULT_MMFF_VDW_CONSTRUCTIONS.load(Ordering::Relaxed);
            if count != 1 {
                discrepancies.push(format!(
                    "post-initialization constructor {index} changed count to {count}"
                ));
            }
        }

        let baseline_ref = &baseline;
        let reports = std::thread::scope(|scope| {
            let handles = (0..8)
                .map(|worker| {
                    scope.spawn(move || {
                        let mut actual_calls = 0;
                        let mut worker_discrepancies = Vec::new();
                        for iteration in 0..16 {
                            let atom_type = WARM_QUERIES[(worker + iteration) % WARM_QUERIES.len()];
                            let before = mmff_vdw_lookup_snapshot(default);
                            let result = default_mmff_vdw();
                            actual_calls += 1;
                            if DEFAULT_MMFF_VDW_CONSTRUCTIONS.load(Ordering::Relaxed) != 1 {
                                worker_discrepancies.push(format!(
                                    "worker {worker} call {iteration}: initializer count changed"
                                ));
                            }
                            if &before != baseline_ref {
                                worker_discrepancies.push(format!(
                                    "worker {worker} call {iteration}: pre-call snapshot changed"
                                ));
                            }
                            match result {
                                Ok(cached) => {
                                    let after = mmff_vdw_lookup_snapshot(cached);
                                    if !std::ptr::eq(cached, default)
                                        || cached as *const MmffVdwCollection as usize
                                            != default_address
                                    {
                                        worker_discrepancies.push(format!(
                                            "worker {worker} call {iteration}: cell address changed"
                                        ));
                                    }
                                    if &after != baseline_ref || &before != &after {
                                        worker_discrepancies.push(format!(
                                            "worker {worker} call {iteration}: full collection snapshot changed"
                                        ));
                                    }
                                    match DEFAULT_SENTINELS
                                        .iter()
                                        .copied()
                                        .find(|(_, source_type, _, _)| *source_type == atom_type)
                                    {
                                        Some((index, _, expected_bits, expected_da)) => {
                                            match cached.get(atom_type) {
                                                Some(row) => {
                                                    if !cached
                                                        .d_params
                                                        .get(index)
                                                        .is_some_and(|expected| {
                                                            std::ptr::eq(row, expected)
                                                        })
                                                    {
                                                        worker_discrepancies.push(format!(
                                                            "worker {worker} call {iteration}: query {atom_type} row pointer changed"
                                                        ));
                                                    }
                                                    let row_bits = [
                                                        row.alpha_i.to_bits(),
                                                        row.n_i.to_bits(),
                                                        row.a_i.to_bits(),
                                                        row.g_i.to_bits(),
                                                        row.r_star.to_bits(),
                                                    ];
                                                    if row_bits != expected_bits
                                                        || row.da != expected_da
                                                    {
                                                        worker_discrepancies.push(format!(
                                                            "worker {worker} call {iteration}: query {atom_type} differs from source sentinel"
                                                        ));
                                                    }
                                                }
                                                None => worker_discrepancies.push(format!(
                                                    "worker {worker} call {iteration}: query {atom_type} missed"
                                                )),
                                            }
                                        }
                                        None => {
                                            if cached.get(atom_type).is_some() {
                                                worker_discrepancies.push(format!(
                                                    "worker {worker} call {iteration}: missing query {atom_type} hit"
                                                ));
                                            }
                                        }
                                    }
                                }
                                Err(error) => worker_discrepancies.push(format!(
                                    "worker {worker} call {iteration}: default getter returned {error}"
                                )),
                            }
                        }
                        (actual_calls, worker_discrepancies)
                    })
                })
                .collect::<Vec<_>>();
            handles
                .into_iter()
                .map(|handle| handle.join().expect("a warm VdW getter thread completes"))
                .collect::<Vec<_>>()
        });

        let mut actual_warm_calls = 0;
        for (worker_calls, worker_discrepancies) in reports {
            actual_warm_calls += worker_calls;
            discrepancies.extend(worker_discrepancies);
        }
        if actual_warm_calls != 128 {
            discrepancies.push(format!(
                "expected 128 actual warm getters, got {actual_warm_calls}"
            ));
        }
        if DEFAULT_MMFF_VDW_CONSTRUCTIONS.load(Ordering::Relaxed) != 1 {
            discrepancies.push("final VdW initializer count is not one".to_owned());
        }
        assert!(
            discrepancies.is_empty(),
            "MMFF VdW default-cache regression discrepancies:\n{}",
            discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_vdw_util_minimum_64_frozen_calls_match_source_bits_and_preserve_inputs() {
        const PROFILE0_SOURCE: &str = concat!(
            "*fixed\n",
            "0.0\t0.2\t12.0\t0.8\t0.5\n",
            "17\t4.0\t2.0\t3.0\t1.25\t-\n",
            "18\t9.0\t3.0\t5.0\t0.75\t-\n",
        );
        const PROFILE1_SOURCE: &str = concat!(
            "*fixed\n",
            "0.0\t-0.0\t-0.0\t-0.0\t-2.0\n",
            "17\t16.0\t4.0\t4.0\t2.0\t-\n",
            "18\t1.0\t1.0\t4.0\t0.5\t-\n",
        );
        const KEYS: [u8; 2] = [17, 18];
        const SOURCE_DA: [u8; 2] = [b'-', b'-'];
        const PROFILE0_HEADER_BITS: [u64; 5] = [
            0x0000_0000_0000_0000,
            0x3fc9_9999_9999_999a,
            0x4028_0000_0000_0000,
            0x3fe9_9999_9999_999a,
            0x3fe0_0000_0000_0000,
        ];
        const PROFILE1_HEADER_BITS: [u64; 5] = [
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0xc000_0000_0000_0000,
        ];
        const PROFILE0_ROW_BITS: [[u64; 5]; 2] = [
            [
                0x4010_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x3ff4_0000_0000_0000,
                0x4008_0000_0000_0000,
            ],
            [
                0x4022_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x4014_0000_0000_0000,
                0x3fe8_0000_0000_0000,
                0x4014_0000_0000_0000,
            ],
        ];
        const PROFILE1_ROW_BITS: [[u64; 5]; 2] = [
            [
                0x4030_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4010_0000_0000_0000,
            ],
            [
                0x3ff0_0000_0000_0000,
                0x3ff0_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x3fe0_0000_0000_0000,
                0x4010_0000_0000_0000,
            ],
        ];
        const DA_VALUES: [u8; 4] = [b'D', b'A', b'-', b'X'];
        const EXPECTED_MINIMUM_BITS: [u64; 32] = [
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4011_b03c_c100_cf53,
            0x4011_b03c_c100_cf53,
            0x4011_b03c_c100_cf53,
            0x4010_0000_0000_0000,
            0x4011_b03c_c100_cf53,
            0x4011_b03c_c100_cf53,
            0x4011_b03c_c100_cf53,
            0x4010_0000_0000_0000,
            0x4011_b03c_c100_cf53,
            0x4011_b03c_c100_cf53,
            0x4011_b03c_c100_cf53,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
            0x4010_0000_0000_0000,
        ];

        let profile_inputs = [
            (PROFILE0_SOURCE, PROFILE0_HEADER_BITS, &PROFILE0_ROW_BITS),
            (PROFILE1_SOURCE, PROFILE1_HEADER_BITS, &PROFILE1_ROW_BITS),
        ];
        let mut discrepancies = Vec::new();
        let mut actual_calls = 0;
        let mut base_cells = 0;

        for (profile_index, (source, expected_header, expected_rows)) in
            profile_inputs.into_iter().enumerate()
        {
            let collection = MmffVdwCollection::from_text(source)
                .expect("the fixed source-shaped VdW profile parses");
            assert_mmff_vdw_literal_collection(
                &collection,
                expected_header,
                &KEYS,
                expected_rows,
                &SOURCE_DA,
            );
            let baseline = mmff_vdw_lookup_snapshot(&collection);
            let source_address = source.as_ptr() as usize;
            let source_bytes = source.as_bytes().to_vec();

            for (i_index, i_da) in DA_VALUES.into_iter().enumerate() {
                for (j_index, j_da) in DA_VALUES.into_iter().enumerate() {
                    let mut i_row = collection.d_params[0];
                    let mut j_row = collection.d_params[1];
                    i_row.da = i_da;
                    j_row.da = j_da;
                    let cell_index = profile_index * 16 + i_index * 4 + j_index;
                    let expected_bits = EXPECTED_MINIMUM_BITS[cell_index];
                    base_cells += 1;

                    for repeat in 0..2 {
                        let case =
                            format!("profile {profile_index} DA {i_da:?}/{j_da:?} repeat {repeat}");
                        let before = mmff_vdw_lookup_snapshot(&collection);
                        let i_address_before = &i_row as *const _ as usize;
                        let j_address_before = &j_row as *const _ as usize;
                        let i_values_before = [
                            i_row.alpha_i.to_bits(),
                            i_row.n_i.to_bits(),
                            i_row.a_i.to_bits(),
                            i_row.g_i.to_bits(),
                            i_row.r_star.to_bits(),
                        ];
                        let j_values_before = [
                            j_row.alpha_i.to_bits(),
                            j_row.n_i.to_bits(),
                            j_row.a_i.to_bits(),
                            j_row.g_i.to_bits(),
                            j_row.r_star.to_bits(),
                        ];
                        let i_da_before = i_row.da;
                        let j_da_before = j_row.da;
                        let source_address_before = source.as_ptr() as usize;
                        let source_bytes_before = source.as_bytes();

                        let actual = super::super::nonbonded::calc_unscaled_vdw_minimum(
                            &collection,
                            &i_row,
                            &j_row,
                        );
                        actual_calls += 1;
                        let after = mmff_vdw_lookup_snapshot(&collection);

                        if before != baseline || after != before {
                            discrepancies.push(format!("{case}: collection snapshot changed"));
                        }
                        if i_values_before != expected_rows[0]
                            || j_values_before != expected_rows[1]
                            || i_da_before != i_da
                            || j_da_before != j_da
                            || [
                                i_row.alpha_i.to_bits(),
                                i_row.n_i.to_bits(),
                                i_row.a_i.to_bits(),
                                i_row.g_i.to_bits(),
                                i_row.r_star.to_bits(),
                            ] != i_values_before
                            || [
                                j_row.alpha_i.to_bits(),
                                j_row.n_i.to_bits(),
                                j_row.a_i.to_bits(),
                                j_row.g_i.to_bits(),
                                j_row.r_star.to_bits(),
                            ] != j_values_before
                            || i_row.da != i_da_before
                            || j_row.da != j_da_before
                        {
                            discrepancies.push(format!("{case}: row values or DA bytes changed"));
                        }
                        if (&i_row as *const _ as usize) != i_address_before
                            || (&j_row as *const _ as usize) != j_address_before
                        {
                            discrepancies.push(format!("{case}: row address changed"));
                        }
                        if source_address_before != source_address
                            || source.as_ptr() as usize != source_address
                            || source_bytes_before != source_bytes.as_slice()
                            || source.as_bytes() != source_bytes.as_slice()
                        {
                            discrepancies.push(format!("{case}: source text changed"));
                        }
                        if actual.to_bits() != expected_bits {
                            discrepancies.push(format!(
                                "{case}: minimum bits {:016x}, expected {:016x}",
                                actual.to_bits(),
                                expected_bits
                            ));
                        }
                    }
                }
            }
        }

        if actual_calls != 64 {
            discrepancies.push(format!(
                "expected 64 actual minimum calls, got {actual_calls}"
            ));
        }
        if base_cells != 32 {
            discrepancies.push(format!("expected 32 base minimum cells, got {base_cells}"));
        }
        assert!(
            discrepancies.is_empty(),
            "MMFF VdW minimum discrepancies:\n{}",
            discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_vdw_util_depth_64_frozen_calls_match_source_bits_and_preserve_inputs() {
        const PROFILE0_SOURCE: &str = concat!(
            "*fixed\n",
            "0.0\t0.2\t12.0\t0.8\t0.5\n",
            "17\t4.0\t2.0\t3.0\t1.25\t-\n",
            "18\t9.0\t3.0\t5.0\t0.75\t-\n",
        );
        const PROFILE1_SOURCE: &str = concat!(
            "*fixed\n",
            "0.0\t-0.0\t-0.0\t-0.0\t-2.0\n",
            "17\t16.0\t4.0\t4.0\t2.0\t-\n",
            "18\t1.0\t1.0\t4.0\t0.5\t-\n",
        );
        const KEYS: [u8; 2] = [17, 18];
        const SOURCE_DA: [u8; 2] = [b'-', b'-'];
        const PROFILE0_HEADER_BITS: [u64; 5] = [
            0x0000_0000_0000_0000,
            0x3fc9_9999_9999_999a,
            0x4028_0000_0000_0000,
            0x3fe9_9999_9999_999a,
            0x3fe0_0000_0000_0000,
        ];
        const PROFILE1_HEADER_BITS: [u64; 5] = [
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0xc000_0000_0000_0000,
        ];
        const PROFILE0_ROW_BITS: [[u64; 5]; 2] = [
            [
                0x4010_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x3ff4_0000_0000_0000,
                0x4008_0000_0000_0000,
            ],
            [
                0x4022_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x4014_0000_0000_0000,
                0x3fe8_0000_0000_0000,
                0x4014_0000_0000_0000,
            ],
        ];
        const PROFILE1_ROW_BITS: [[u64; 5]; 2] = [
            [
                0x4030_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4010_0000_0000_0000,
            ],
            [
                0x3ff0_0000_0000_0000,
                0x3ff0_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x3fe0_0000_0000_0000,
                0x4010_0000_0000_0000,
            ],
        ];
        const DA_VALUES: [u8; 4] = [b'D', b'A', b'-', b'X'];
        const DEPTH_RADIUS_BITS: [u64; 2] = [0x4006_0000_0000_0000, 0x400c_0000_0000_0000];
        const EXPECTED_DEPTH_BITS: [u64; 32] = [
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x4011_f8eb_7d5c_f77a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
            0x3fe0_d1b0_8942_740a,
        ];

        let profile_inputs = [
            (PROFILE0_SOURCE, PROFILE0_HEADER_BITS, &PROFILE0_ROW_BITS),
            (PROFILE1_SOURCE, PROFILE1_HEADER_BITS, &PROFILE1_ROW_BITS),
        ];
        let mut discrepancies = Vec::new();
        let mut actual_calls = 0;
        let mut base_cells = 0;

        for (profile_index, (source, expected_header, expected_rows)) in
            profile_inputs.into_iter().enumerate()
        {
            let collection = MmffVdwCollection::from_text(source)
                .expect("the fixed source-shaped VdW profile parses");
            assert_mmff_vdw_literal_collection(
                &collection,
                expected_header,
                &KEYS,
                expected_rows,
                &SOURCE_DA,
            );
            let baseline = mmff_vdw_lookup_snapshot(&collection);
            let source_address = source.as_ptr() as usize;
            let source_bytes = source.as_bytes().to_vec();
            let depth_radius = if profile_index == 0 {
                2.75_f64
            } else {
                3.5_f64
            };
            let expected_radius_bits = DEPTH_RADIUS_BITS[profile_index];
            assert_eq!(depth_radius.to_bits(), expected_radius_bits);

            for (i_index, i_da) in DA_VALUES.into_iter().enumerate() {
                for (j_index, j_da) in DA_VALUES.into_iter().enumerate() {
                    let mut i_row = collection.d_params[0];
                    let mut j_row = collection.d_params[1];
                    i_row.da = i_da;
                    j_row.da = j_da;
                    let cell_index = profile_index * 16 + i_index * 4 + j_index;
                    let expected_bits = EXPECTED_DEPTH_BITS[cell_index];
                    base_cells += 1;

                    for repeat in 0..2 {
                        let case =
                            format!("profile {profile_index} DA {i_da:?}/{j_da:?} repeat {repeat}");
                        let before = mmff_vdw_lookup_snapshot(&collection);
                        let i_address_before = &i_row as *const _ as usize;
                        let j_address_before = &j_row as *const _ as usize;
                        let i_values_before = [
                            i_row.alpha_i.to_bits(),
                            i_row.n_i.to_bits(),
                            i_row.a_i.to_bits(),
                            i_row.g_i.to_bits(),
                            i_row.r_star.to_bits(),
                        ];
                        let j_values_before = [
                            j_row.alpha_i.to_bits(),
                            j_row.n_i.to_bits(),
                            j_row.a_i.to_bits(),
                            j_row.g_i.to_bits(),
                            j_row.r_star.to_bits(),
                        ];
                        let i_da_before = i_row.da;
                        let j_da_before = j_row.da;
                        let source_address_before = source.as_ptr() as usize;
                        let source_bytes_before = source.as_bytes();
                        let radius_bits_before = depth_radius.to_bits();

                        let actual = super::super::nonbonded::calc_unscaled_vdw_well_depth(
                            depth_radius,
                            &i_row,
                            &j_row,
                        );
                        actual_calls += 1;
                        let after = mmff_vdw_lookup_snapshot(&collection);

                        if before != baseline || after != before {
                            discrepancies.push(format!("{case}: collection snapshot changed"));
                        }
                        if i_values_before != expected_rows[0]
                            || j_values_before != expected_rows[1]
                            || i_da_before != i_da
                            || j_da_before != j_da
                            || [
                                i_row.alpha_i.to_bits(),
                                i_row.n_i.to_bits(),
                                i_row.a_i.to_bits(),
                                i_row.g_i.to_bits(),
                                i_row.r_star.to_bits(),
                            ] != i_values_before
                            || [
                                j_row.alpha_i.to_bits(),
                                j_row.n_i.to_bits(),
                                j_row.a_i.to_bits(),
                                j_row.g_i.to_bits(),
                                j_row.r_star.to_bits(),
                            ] != j_values_before
                            || i_row.da != i_da_before
                            || j_row.da != j_da_before
                        {
                            discrepancies.push(format!("{case}: row values or DA bytes changed"));
                        }
                        if (&i_row as *const _ as usize) != i_address_before
                            || (&j_row as *const _ as usize) != j_address_before
                        {
                            discrepancies.push(format!("{case}: row address changed"));
                        }
                        if source_address_before != source_address
                            || source.as_ptr() as usize != source_address
                            || source_bytes_before != source_bytes.as_slice()
                            || source.as_bytes() != source_bytes.as_slice()
                        {
                            discrepancies.push(format!("{case}: source text changed"));
                        }
                        if radius_bits_before != expected_radius_bits
                            || depth_radius.to_bits() != expected_radius_bits
                        {
                            discrepancies.push(format!("{case}: depth radius input changed"));
                        }
                        if actual.to_bits() != expected_bits {
                            discrepancies.push(format!(
                                "{case}: well-depth bits {:016x}, expected {:016x}",
                                actual.to_bits(),
                                expected_bits
                            ));
                        }
                    }
                }
            }
        }

        if actual_calls != 64 {
            discrepancies.push(format!(
                "expected 64 actual depth calls, got {actual_calls}"
            ));
        }
        if base_cells != 32 {
            discrepancies.push(format!("expected 32 base depth cells, got {base_cells}"));
        }
        assert!(
            discrepancies.is_empty(),
            "MMFF VdW well-depth discrepancies:\n{}",
            discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_vdw_util_scale_72_frozen_calls_match_source_bits_and_preserve_inputs() {
        const PROFILE0_SOURCE: &str = concat!(
            "*fixed\n",
            "0.0\t0.2\t12.0\t0.8\t0.5\n",
            "17\t4.0\t2.0\t3.0\t1.25\t-\n",
            "18\t9.0\t3.0\t5.0\t0.75\t-\n",
        );
        const PROFILE1_SOURCE: &str = concat!(
            "*fixed\n",
            "0.0\t-0.0\t-0.0\t-0.0\t-2.0\n",
            "17\t16.0\t4.0\t4.0\t2.0\t-\n",
            "18\t1.0\t1.0\t4.0\t0.5\t-\n",
        );
        const KEYS: [u8; 2] = [17, 18];
        const SOURCE_DA: [u8; 2] = [b'-', b'-'];
        const PROFILE0_HEADER_BITS: [u64; 5] = [
            0x0000_0000_0000_0000,
            0x3fc9_9999_9999_999a,
            0x4028_0000_0000_0000,
            0x3fe9_9999_9999_999a,
            0x3fe0_0000_0000_0000,
        ];
        const PROFILE1_HEADER_BITS: [u64; 5] = [
            0x0000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x8000_0000_0000_0000,
            0xc000_0000_0000_0000,
        ];
        const PROFILE0_ROW_BITS: [[u64; 5]; 2] = [
            [
                0x4010_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x3ff4_0000_0000_0000,
                0x4008_0000_0000_0000,
            ],
            [
                0x4022_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x4014_0000_0000_0000,
                0x3fe8_0000_0000_0000,
                0x4014_0000_0000_0000,
            ],
        ];
        const PROFILE1_ROW_BITS: [[u64; 5]; 2] = [
            [
                0x4030_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4010_0000_0000_0000,
            ],
            [
                0x3ff0_0000_0000_0000,
                0x3ff0_0000_0000_0000,
                0x4010_0000_0000_0000,
                0x3fe0_0000_0000_0000,
                0x4010_0000_0000_0000,
            ],
        ];
        const DA_VALUES: [u8; 4] = [b'D', b'A', b'-', b'X'];
        const INPUT_BITS: [[u64; 2]; 2] = [
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
        ];
        const EXPECTED_BASE_BITS: [[u64; 2]; 32] = [
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ff6_6666_6666_6667, 0x3fd4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ff6_6666_6666_6667, 0x3fd4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x3ffc_0000_0000_0000, 0x3fe4_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x8000_0000_0000_0000, 0x4008_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x8000_0000_0000_0000, 0x4008_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
            [0x4002_0000_0000_0000, 0xbff8_0000_0000_0000],
        ];
        const CONTROL_CASES: [([u8; 2], [u64; 2], [u64; 2]); 8] = [
            (
                [b'D', b'A'],
                [0x8000_0000_0000_0000, 0x0000_0000_0000_0000],
                [0x8000_0000_0000_0000, 0x0000_0000_0000_0000],
            ),
            (
                [b'D', b'A'],
                [0x7ff0_0000_0000_0000, 0xfff0_0000_0000_0000],
                [0x7ff0_0000_0000_0000, 0xfff0_0000_0000_0000],
            ),
            (
                [b'A', b'D'],
                [0x8000_0000_0000_0000, 0x0000_0000_0000_0000],
                [0x8000_0000_0000_0000, 0x0000_0000_0000_0000],
            ),
            (
                [b'A', b'D'],
                [0x7ff0_0000_0000_0000, 0xfff0_0000_0000_0000],
                [0x7ff0_0000_0000_0000, 0xfff0_0000_0000_0000],
            ),
            (
                [b'D', b'D'],
                [0x8000_0000_0000_0000, 0x0000_0000_0000_0000],
                [0x8000_0000_0000_0000, 0x0000_0000_0000_0000],
            ),
            (
                [b'D', b'D'],
                [0x7ff0_0000_0000_0000, 0xfff0_0000_0000_0000],
                [0x7ff0_0000_0000_0000, 0xfff0_0000_0000_0000],
            ),
            (
                [b'-', b'A'],
                [0x8000_0000_0000_0000, 0x0000_0000_0000_0000],
                [0x8000_0000_0000_0000, 0x0000_0000_0000_0000],
            ),
            (
                [b'-', b'A'],
                [0x7ff0_0000_0000_0000, 0xfff0_0000_0000_0000],
                [0x7ff0_0000_0000_0000, 0xfff0_0000_0000_0000],
            ),
        ];

        fn check_scale_case(
            case: &str,
            collection: &MmffVdwCollection,
            baseline: &MmffVdwLookupSnapshot,
            source: &str,
            i_row: &MmffVdw,
            j_row: &MmffVdw,
            expected_rows: [[u64; 5]; 2],
            expected_da: [u8; 2],
            input_bits: [u64; 2],
            expected_bits: [u64; 2],
            actual_calls: &mut usize,
            discrepancies: &mut Vec<String>,
        ) {
            let before = mmff_vdw_lookup_snapshot(collection);
            let i_address_before = i_row as *const MmffVdw as usize;
            let j_address_before = j_row as *const MmffVdw as usize;
            let i_values_before = [
                i_row.alpha_i.to_bits(),
                i_row.n_i.to_bits(),
                i_row.a_i.to_bits(),
                i_row.g_i.to_bits(),
                i_row.r_star.to_bits(),
            ];
            let j_values_before = [
                j_row.alpha_i.to_bits(),
                j_row.n_i.to_bits(),
                j_row.a_i.to_bits(),
                j_row.g_i.to_bits(),
                j_row.r_star.to_bits(),
            ];
            let i_da_before = i_row.da;
            let j_da_before = j_row.da;
            let source_address_before = source.as_ptr() as usize;
            let source_bytes_before = source.as_bytes().to_vec();
            let mut radius = f64::from_bits(input_bits[0]);
            let mut depth = f64::from_bits(input_bits[1]);
            let radius_address_before = &radius as *const f64 as usize;
            let depth_address_before = &depth as *const f64 as usize;
            let radius_bits_before = radius.to_bits();
            let depth_bits_before = depth.to_bits();

            super::super::nonbonded::scale_vdw_params(
                &mut radius,
                &mut depth,
                collection,
                i_row,
                j_row,
            );
            *actual_calls += 1;
            let after = mmff_vdw_lookup_snapshot(collection);

            if &before != baseline || after != before {
                discrepancies.push(format!("{case}: collection snapshot changed"));
            }
            if i_values_before != expected_rows[0]
                || j_values_before != expected_rows[1]
                || i_da_before != expected_da[0]
                || j_da_before != expected_da[1]
                || [
                    i_row.alpha_i.to_bits(),
                    i_row.n_i.to_bits(),
                    i_row.a_i.to_bits(),
                    i_row.g_i.to_bits(),
                    i_row.r_star.to_bits(),
                ] != i_values_before
                || [
                    j_row.alpha_i.to_bits(),
                    j_row.n_i.to_bits(),
                    j_row.a_i.to_bits(),
                    j_row.g_i.to_bits(),
                    j_row.r_star.to_bits(),
                ] != j_values_before
                || i_row.da != i_da_before
                || j_row.da != j_da_before
            {
                discrepancies.push(format!("{case}: row values or DA bytes changed"));
            }
            if i_row as *const MmffVdw as usize != i_address_before
                || j_row as *const MmffVdw as usize != j_address_before
            {
                discrepancies.push(format!("{case}: row address changed"));
            }
            if source_address_before != source.as_ptr() as usize
                || source.as_bytes() != source_bytes_before.as_slice()
            {
                discrepancies.push(format!("{case}: source text changed"));
            }
            if radius_address_before != &radius as *const f64 as usize
                || depth_address_before != &depth as *const f64 as usize
                || radius_bits_before != input_bits[0]
                || depth_bits_before != input_bits[1]
            {
                discrepancies.push(format!("{case}: scalar input bits or addresses changed"));
            }
            if radius.to_bits() != expected_bits[0] || depth.to_bits() != expected_bits[1] {
                discrepancies.push(format!(
                    "{case}: scaled bits {:016x},{:016x}, expected {:016x},{:016x}",
                    radius.to_bits(),
                    depth.to_bits(),
                    expected_bits[0],
                    expected_bits[1]
                ));
            }
        }

        let profile_inputs = [
            (PROFILE0_SOURCE, PROFILE0_HEADER_BITS, PROFILE0_ROW_BITS),
            (PROFILE1_SOURCE, PROFILE1_HEADER_BITS, PROFILE1_ROW_BITS),
        ];
        let mut discrepancies = Vec::new();
        let mut actual_calls = 0;
        let mut base_cells = 0;
        let mut base_calls = 0;
        let mut control_calls = 0;

        for (profile_index, (source, expected_header, expected_rows)) in
            profile_inputs.into_iter().enumerate()
        {
            let collection = MmffVdwCollection::from_text(source)
                .expect("the fixed source-shaped VdW profile parses");
            assert_mmff_vdw_literal_collection(
                &collection,
                expected_header,
                &KEYS,
                &expected_rows,
                &SOURCE_DA,
            );
            let baseline = mmff_vdw_lookup_snapshot(&collection);

            for (i_index, i_da) in DA_VALUES.into_iter().enumerate() {
                for (j_index, j_da) in DA_VALUES.into_iter().enumerate() {
                    let i_row = MmffVdw {
                        da: i_da,
                        ..collection.d_params[0]
                    };
                    let j_row = MmffVdw {
                        da: j_da,
                        ..collection.d_params[1]
                    };
                    let cell_index = profile_index * 16 + i_index * 4 + j_index;
                    base_cells += 1;

                    for repeat in 0..2 {
                        let case = format!(
                            "base profile {profile_index} DA {i_da:?}/{j_da:?} repeat {repeat}"
                        );
                        check_scale_case(
                            &case,
                            &collection,
                            &baseline,
                            source,
                            &i_row,
                            &j_row,
                            expected_rows,
                            [i_da, j_da],
                            INPUT_BITS[profile_index],
                            EXPECTED_BASE_BITS[cell_index],
                            &mut actual_calls,
                            &mut discrepancies,
                        );
                        base_calls += 1;
                    }
                }
            }

            if profile_index == 0 {
                for (control_index, (da_pair, input_bits, expected_bits)) in
                    CONTROL_CASES.into_iter().enumerate()
                {
                    let i_row = MmffVdw {
                        da: da_pair[0],
                        ..collection.d_params[0]
                    };
                    let j_row = MmffVdw {
                        da: da_pair[1],
                        ..collection.d_params[1]
                    };
                    let case = format!("control {control_index} DA {da_pair:?}");
                    check_scale_case(
                        &case,
                        &collection,
                        &baseline,
                        source,
                        &i_row,
                        &j_row,
                        expected_rows,
                        da_pair,
                        input_bits,
                        expected_bits,
                        &mut actual_calls,
                        &mut discrepancies,
                    );
                    control_calls += 1;
                }
            }
        }

        if actual_calls != 72 {
            discrepancies.push(format!(
                "expected 72 actual scale calls, got {actual_calls}"
            ));
        }
        if base_cells != 32 {
            discrepancies.push(format!("expected 32 base scale cells, got {base_cells}"));
        }
        if base_calls != 64 {
            discrepancies.push(format!("expected 64 base scale calls, got {base_calls}"));
        }
        if control_calls != 8 {
            discrepancies.push(format!(
                "expected 8 supplementary scale calls, got {control_calls}"
            ));
        }
        assert!(
            discrepancies.is_empty(),
            "MMFF VdW scaling discrepancies:\n{}",
            discrepancies.join("\n")
        );
    }
}
