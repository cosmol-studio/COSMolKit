//! Pinned RDKit Python Utils.ExplainAtomCode over shared source constants.
use cosmolkit_model::Element;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomCodeExplanation {
    symbol: &'static str,
    branch_count: u32,
    pi_electrons: u32,
    chirality: Option<&'static str>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomCodeExplanationError {
    UnknownChirality { code: u8 },
}
impl fmt::Display for AtomCodeExplanationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownChirality { code } => write!(f, "{code}"),
        }
    }
}
impl std::error::Error for AtomCodeExplanationError {}

impl AtomCodeExplanation {
    pub fn from_code(
        mut code: u64,
        _branch_subtract: i64,
        include_chirality: bool,
    ) -> Result<Self, AtomCodeExplanationError> {
        // RDKit❗✔️: def ExplainAtomCode(code, branchSubtract=0, includeChirality=False):
        // RDKit❗✔️:   typeMask = (1 << rdMolDescriptors.AtomPairsParameters.numTypeBits) - 1
        // RDKit❗✔️:   branchMask = (1 << rdMolDescriptors.AtomPairsParameters.numBranchBits) - 1
        // RDKit❗✔️:   piMask = (1 << rdMolDescriptors.AtomPairsParameters.numPiBits) - 1
        // RDKit❗✔️:   chiMask = (1 << rdMolDescriptors.AtomPairsParameters.numChiralBits) - 1
        // RDKit❗✔️:   nBranch = int(code & branchMask)
        // RDKit❗✔️:   code = code >> rdMolDescriptors.AtomPairsParameters.numBranchBits
        // RDKit❗✔️:   nPi = int(code & piMask)
        // RDKit❗✔️:   code = code >> rdMolDescriptors.AtomPairsParameters.numPiBits
        // RDKit❗✔️:   typeIdx = int(code & typeMask)
        // RDKit❗✔️:   if typeIdx < len(rdMolDescriptors.AtomPairsParameters.atomTypes):
        // RDKit❗✔️:     atomNum = rdMolDescriptors.AtomPairsParameters.atomTypes[typeIdx]
        // RDKit❗✔️:     atomSymbol = Chem.GetPeriodicTable().GetElementSymbol(atomNum)
        // RDKit❗✔️:   else:
        // RDKit❗✔️:     atomSymbol = 'X'
        // RDKit❗✔️:   if not includeChirality:
        // RDKit❗✔️:     return (atomSymbol, nBranch, nPi)
        // RDKit❗✔️:   code = code >> rdMolDescriptors.AtomPairsParameters.numTypeBits
        // RDKit❗✔️:   chiDict = {0: '', 1: 'R', 2: 'S'}
        // RDKit❗✔️:   chiCode = int(code & chiMask)
        // RDKit❗✔️:   return (atomSymbol, nBranch, nPi, chiDict[chiCode])
        // Source branchSubtract is intentionally unused; adding it would alter
        // source behavior. Constant-width masks and the sole Element symbol
        // owner preserve O(1) decoding, with no topology or copied tables.
        // The wrapper exposes the entire 16-entry C++ array, including zero:
        // RDKit❗✔️: std::vector<unsigned int> atomPairTypes(
        // RDKit❗✔️:     RDKit::AtomPairs::atomNumberTypes,
        // RDKit❗✔️:     RDKit::AtomPairs::atomNumberTypes +
        // RDKit❗✔️:         sizeof(RDKit::AtomPairs::atomNumberTypes) / sizeof(unsigned int));
        let branch_count =
            (code & ((1 << crate::AtomPairsParameters::num_branch_bits()) - 1)) as u32;
        code >>= crate::AtomPairsParameters::num_branch_bits();
        let pi_electrons = (code & ((1 << crate::AtomPairsParameters::num_pi_bits()) - 1)) as u32;
        code >>= crate::AtomPairsParameters::num_pi_bits();
        let atomic_number = crate::atom_code::ATOM_NUMBER_TYPES
            [(code & ((1 << crate::AtomPairsParameters::num_type_bits()) - 1)) as usize];
        let symbol = Element::from_atomic_number(atomic_number as u8)
            .expect("source atom number table contains valid Elements")
            .symbol();
        let chirality = if include_chirality {
            code >>= crate::AtomPairsParameters::num_type_bits();
            Some(
                match (code & ((1 << crate::AtomPairsParameters::num_chiral_bits()) - 1)) as u8 {
                    0 => "",
                    1 => "R",
                    2 => "S",
                    value => {
                        return Err(AtomCodeExplanationError::UnknownChirality { code: value });
                    }
                },
            )
        } else {
            None
        };
        Ok(Self {
            symbol,
            branch_count,
            pi_electrons,
            chirality,
        })
    }
    pub fn symbol(&self) -> &'static str {
        self.symbol
    }
    pub fn branch_count(&self) -> u32 {
        self.branch_count
    }
    pub fn pi_electrons(&self) -> u32 {
        self.pi_electrons
    }
    pub fn chirality(&self) -> Option<&'static str> {
        self.chirality
    }
}
