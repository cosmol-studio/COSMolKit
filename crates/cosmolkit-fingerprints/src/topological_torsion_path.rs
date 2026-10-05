//! Pinned RDKit351f8f Torsions.py path-score and explanation owner.
use crate::{AtomCodeError, AtomCodeInput, AtomCodeOptions, AtomPairsParameters, FingerprintError};
use cosmolkit_core::ValenceAssignment;
use cosmolkit_model::{AtomId, MoleculeProperties, TopologyBlock, TopologyValidationError};
use std::fmt;

#[derive(Debug)]
pub enum TopologicalTorsionPathScoreError {
    ZeroSize,
    ShortPath {
        actual: usize,
        required: usize,
    },
    ShortAtomCodes {
        actual: usize,
        required: usize,
    },
    AtomIndexOutOfRange {
        index: usize,
        atom_count: usize,
    },
    AtomCodeUnderflow {
        index: usize,
        code: u32,
        subtract: u32,
    },
    InvalidTopology(TopologyValidationError),
    AtomCode(AtomCodeError),
    PackedCode(FingerprintError),
}
impl fmt::Display for TopologicalTorsionPathScoreError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroSize => f.write_str("size must be greater than zero"),
            Self::ShortPath { .. } => f.write_str("path must contain at least size atom indices"),
            Self::ShortAtomCodes { .. } => {
                f.write_str("atom_codes must contain at least one entry per molecule atom")
            }
            Self::AtomIndexOutOfRange { index, .. } => {
                write!(f, "atom index {index} is out of range")
            }
            Self::AtomCodeUnderflow { .. } => {
                f.write_str("atom code is smaller than the path branch subtraction")
            }
            Self::InvalidTopology(e) => write!(f, "{e}"),
            Self::AtomCode(e) => write!(f, "{e}"),
            Self::PackedCode(e) => write!(f, "{e}"),
        }
    }
}
impl std::error::Error for TopologicalTorsionPathScoreError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidTopology(e) => Some(e),
            Self::AtomCode(e) => Some(e),
            Self::PackedCode(e) => Some(e),
            _ => None,
        }
    }
}

/// Preserve original modern validated u64 transport and source path ordering.
/// No distinctness or connectivity check; unused trailing indices are untouched.
pub fn topological_torsion_path_score(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: Option<&ValenceAssignment>,
    path: &[usize],
    size: usize,
    atom_codes: Option<&[u32]>,
) -> Result<u64, TopologicalTorsionPathScoreError> {
    // RDKit❗❌:   codes = [None] * size
    // RDKit❗❌:   for i in range(size):
    // RDKit❗❌:     if i == 0 or i == size - 1:
    // RDKit❗❌:       sub = 1
    // RDKit❗❌:     else:
    // RDKit❗❌:       sub = 2
    // RDKit❗❌:     if not atomCodes:
    // RDKit❗❌:       codes[i] = Utils.GetAtomCode(mol.GetAtomWithIdx(path[i]), sub)
    // RDKit❗❌:     else:
    // RDKit❗❌:       codes[i] = atomCodes[path[i]] - sub
    // RDKit❗❌:
    // RDKit❗❌:   # "canonicalize" the code vector:
    // RDKit❗❌:   beg = 0
    // RDKit❗❌:   end = len(codes) - 1
    // RDKit❗❌:   while beg < end:
    // RDKit❗❌:     if codes[beg] == codes[end]:
    // RDKit❗❌:       beg += 1
    // RDKit❗❌:       end -= 1
    // RDKit❗❌:     else:
    // RDKit❗❌:       if codes[beg] > codes[end]:
    // RDKit❗❌:         codes.reverse()
    // RDKit❗❌:       break
    // RDKit❗❌:
    // RDKit❗❌:   accum = 0
    // RDKit❗❌:   codeSize = rdMolDescriptors.AtomPairsParameters.codeSize
    // RDKit❗❌:   for i, code in enumerate(codes):
    // RDKit❗❌:     accum |= code << (codeSize * i)
    // RDKit❗❌:   return accum
    // Complexity: one Vec(size), source paired-end packing reused without a
    // reverse clone. The default atom-code input additionally validates topology
    // once in O(A+B); this cost exceeds the source's selected-atom local reads.
    // Reuse the same validated borrowed input through all selected atoms; no
    // whole molecule/block clone, CIP or repeated structural validation.
    // Original d892ec3 binding explicitly validates size/path/custom lengths,
    // selected indexes and unsigned subtraction before calling the C++ u64
    // packed-code helper. Its full-u32 OR truncation and approved width errors
    // are owned by topological_torsion_code; unbounded Python integers differ.
    use TopologicalTorsionPathScoreError as E;
    if size == 0 {
        return Err(E::ZeroSize);
    }
    if path.len() < size {
        return Err(E::ShortPath {
            actual: path.len(),
            required: size,
        });
    }
    let atom_codes = atom_codes.filter(|v| !v.is_empty());
    if let Some(v) = atom_codes {
        if v.len() < topology.atoms.len() {
            return Err(E::ShortAtomCodes {
                actual: v.len(),
                required: topology.atoms.len(),
            });
        }
    }
    let mut codes = Vec::with_capacity(size);
    if let Some(custom) = atom_codes {
        for (position, &index) in path.iter().take(size).enumerate() {
            if index >= topology.atoms.len() {
                return Err(E::AtomIndexOutOfRange {
                    index,
                    atom_count: topology.atoms.len(),
                });
            }
            let subtract = if position == 0 || position + 1 == size {
                1
            } else {
                2
            };
            codes.push(
                custom[index]
                    .checked_sub(subtract)
                    .ok_or(E::AtomCodeUnderflow {
                        index,
                        code: custom[index],
                        subtract,
                    })?,
            );
        }
    } else {
        let mut input =
            AtomCodeInput::from_ref(topology, properties).map_err(E::InvalidTopology)?;
        for (position, &index) in path.iter().take(size).enumerate() {
            if index >= topology.atoms.len() {
                return Err(E::AtomIndexOutOfRange {
                    index,
                    atom_count: topology.atoms.len(),
                });
            }
            let subtract = if position == 0 || position + 1 == size {
                1
            } else {
                2
            };
            let explicit = valence
                .and_then(|v| v.explicit_valence.get(index))
                .map(|v| *v as i8);
            let assignment = crate::atom_code(
                input,
                Some(AtomId::new(index)),
                explicit,
                &AtomCodeOptions {
                    branch_subtract: subtract,
                    include_chirality: false,
                    use_legacy_stereo_perception: true,
                },
            )
            .map_err(E::AtomCode)?;
            codes.push(assignment.code());
            input = assignment.into_input();
        }
    }
    crate::topological_torsion_code(&codes, false).map_err(E::PackedCode)
}

/// Complete source decoding, including size zero and zero chunks after u64 ends.
pub fn explain_path_score(mut score: u64, size: usize) -> Vec<(&'static str, u32, u32)> {
    // RDKit✔️✔️:   codeSize: int = rdMolDescriptors.AtomPairsParameters.codeSize
    // RDKit✔️✔️:   codeMask = (1 << codeSize) - 1
    // RDKit✔️✔️:   res = [None] * size
    // RDKit✔️✔️:   for i in range(size):
    // RDKit✔️✔️:     if i == 0 or i == size - 1:
    // RDKit✔️✔️:       sub = 1
    // RDKit✔️✔️:     else:
    // RDKit✔️✔️:       sub = 2
    // RDKit✔️✔️:     code = score & codeMask
    // RDKit✔️✔️:     score = score >> codeSize
    // RDKit✔️✔️:     symb, nBranch, nPi = Utils.ExplainAtomCode(code)
    // RDKit✔️✔️:     res[i] = (symb, nBranch + sub, nPi)
    // Source O(size) decoding and one result allocation; no topology, sorting,
    // source-code masks or atom-type tables duplicated in this helper.
    let mask = (1u64 << AtomPairsParameters::code_size()) - 1;
    let mut result = Vec::with_capacity(size);
    for position in 0..size {
        let subtract = if position == 0 || position + 1 == size {
            1
        } else {
            2
        };
        let code = score & mask;
        score >>= AtomPairsParameters::code_size();
        // Chirality is false and masked9bit codes cannot produce its error.
        let decoded = crate::AtomCodeExplanation::from_code(code, 0, false)
            .expect("false-chirality source atom-code decoding is infallible");
        result.push((
            decoded.symbol(),
            decoded.branch_count() + subtract,
            decoded.pi_electrons(),
        ));
    }
    result
}
