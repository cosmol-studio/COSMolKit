//! Typed state written by the modern RDKit-backed CIP labeler.

use std::{fmt, str::FromStr};

use crate::{AtomId, BondId, BondOrder};

/// Reproduce pinned RDKit `ROMol::clearComputedProps()` for stored properties.
///
/// Clearing follows registered computed-property membership rather than
/// property names. This preserves a non-computed property even when a caller
/// deliberately gives it a conventional computed-property name.
pub(crate) fn clear_computed_source_properties(
    topology: &mut crate::molecule::TopologyBlock,
    properties: &mut crate::MoleculeProperties,
) {
    // BEGIN RDKIT CPP FUNCTION ROMol::clearComputedProps (ROMol.cpp:588-603)
    // RDKit✔️✔️: void ROMol::clearComputedProps(bool includeRings) const {
    // RDKit✔️✔️:   RDProps::clearComputedProps();
    // RDKit✔️✔️:   for (auto atom : atoms()) {
    // RDKit✔️✔️:     atom->clearComputedProps();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto bond : bonds()) {
    // RDKit✔️✔️:     bond->clearComputedProps();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ROMol::clearComputedProps
    properties.clear_computed_props();
    for atom in &mut topology.atoms {
        atom.clear_computed_props();
    }
    for bond in &mut topology.bonds {
        bond.clear_computed_props();
    }
}

/// Options for modern molecular-context CIP assignment.
///
/// When both selections are absent or empty, every atom and bond is considered,
/// matching the pinned RDKit wrapper's truth-value dispatch. Once either
/// selection is non-empty, an absent or empty category selects no entries.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct CipLabelOptions {
    atom_indices: Option<Vec<AtomId>>,
    bond_indices: Option<Vec<BondId>>,
    max_recursive_iterations: u32,
}

impl CipLabelOptions {
    #[must_use]
    pub fn with_atoms(mut self, atom_indices: impl IntoIterator<Item = AtomId>) -> Self {
        self.atom_indices = Some(atom_indices.into_iter().collect());
        self
    }

    #[must_use]
    pub fn with_bonds(mut self, bond_indices: impl IntoIterator<Item = BondId>) -> Self {
        self.bond_indices = Some(bond_indices.into_iter().collect());
        self
    }

    #[must_use]
    pub const fn with_max_recursive_iterations(mut self, limit: u32) -> Self {
        self.max_recursive_iterations = limit;
        self
    }

    #[must_use]
    pub fn atom_indices(&self) -> Option<&[AtomId]> {
        self.atom_indices.as_deref()
    }

    #[must_use]
    pub fn bond_indices(&self) -> Option<&[BondId]> {
        self.bond_indices.as_deref()
    }

    #[must_use]
    pub const fn max_recursive_iterations(&self) -> u32 {
        self.max_recursive_iterations
    }
}

/// A descriptor emitted by the supported modern CIP assignment dispatcher.
///
/// Uppercase and lowercase descriptors are distinct. Lowercase variants are
/// pseudoasymmetric descriptors, not aliases for their uppercase forms.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CipDescriptor {
    R,
    S,
    LowerR,
    LowerS,
    E,
    Z,
    M,
    P,
    LowerM,
    LowerP,
}

impl CipDescriptor {
    /// Returns the stable descriptor spelling stored in `_CIPCode`.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        // BEGIN RDKIT CPP FUNCTION emitted descriptor branches (CIPLabeler/Descriptor.h)
        // RDKit✔️✔️:     case Descriptor::R:
        // RDKit✔️✔️:       return "R";
        // RDKit✔️✔️:     case Descriptor::S:
        // RDKit✔️✔️:       return "S";
        // RDKit✔️✔️:     case Descriptor::r:
        // RDKit✔️✔️:       return "r";
        // RDKit✔️✔️:     case Descriptor::s:
        // RDKit✔️✔️:       return "s";
        // RDKit✔️✔️:     case Descriptor::E:
        // RDKit✔️✔️:       return "E";
        // RDKit✔️✔️:     case Descriptor::Z:
        // RDKit✔️✔️:       return "Z";
        // RDKit✔️✔️:     case Descriptor::M:
        // RDKit✔️✔️:       return "M";
        // RDKit✔️✔️:     case Descriptor::P:
        // RDKit✔️✔️:       return "P";
        // RDKit✔️✔️:     case Descriptor::m:
        // RDKit✔️✔️:       return "m";
        // RDKit✔️✔️:     case Descriptor::p:
        // RDKit✔️✔️:       return "p";
        // END RDKIT CPP FUNCTION emitted descriptor branches
        match self {
            Self::R => "R",
            Self::S => "S",
            Self::LowerR => "r",
            Self::LowerS => "s",
            Self::E => "E",
            Self::Z => "Z",
            Self::M => "M",
            Self::P => "P",
            Self::LowerM => "m",
            Self::LowerP => "p",
        }
    }
}

impl fmt::Display for CipDescriptor {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

impl FromStr for CipDescriptor {
    type Err = CipDescriptorError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "R" => Ok(Self::R),
            "S" => Ok(Self::S),
            "r" => Ok(Self::LowerR),
            "s" => Ok(Self::LowerS),
            "E" => Ok(Self::E),
            "Z" => Ok(Self::Z),
            "M" => Ok(Self::M),
            "P" => Ok(Self::P),
            "m" => Ok(Self::LowerM),
            "p" => Ok(Self::LowerP),
            _ => Err(CipDescriptorError::InvalidStoredDescriptor {
                value: value.to_owned(),
            }),
        }
    }
}

/// Error returned when persisted `_CIPCode` state is not a descriptor emitted
/// by the supported modern assignment dispatcher.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CipDescriptorError {
    #[error("invalid stored modern CIP descriptor `{value}`")]
    InvalidStoredDescriptor { value: String },
}

/// Structured failures from the modern source-backed CIP labeler.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CipLabelerError {
    #[error("CIPLabeler atom index {index} is out of range for {atom_count} atoms")]
    AtomIndexOutOfRange { index: usize, atom_count: usize },
    #[error("CIPLabeler bond index {index} is out of range for {bond_count} bonds")]
    BondIndexOutOfRange { index: usize, bond_count: usize },
    #[error("CIPLabeler atom selection mask has length {actual}, expected {expected}")]
    AtomMaskLengthMismatch { actual: usize, expected: usize },
    #[error("CIPLabeler bond selection mask has length {actual}, expected {expected}")]
    BondMaskLengthMismatch { actual: usize, expected: usize },
    #[error("CIPLabeler {kind} index {index} exceeds the source unsigned-int width")]
    SourceIndexWidthExceeded { kind: &'static str, index: usize },
    #[error("CIPLabeler bond {bond} is not incident to atom {atom}")]
    BondNotIncident { bond: usize, atom: usize },
    #[error("CIPLabeler non integer-order bond is not allowed: {order:?}")]
    NonIntegerBondOrder { order: BondOrder },
    #[error("CIPLabeler node {node} is not an endpoint of edge {edge}")]
    EdgeEndpointMismatch { edge: usize, node: usize },
    #[error("Digraph generation failed: more than {limit}nodes found.")]
    TooManyNodes { limit: usize },
    #[error("Max Iterations Exceeded in CIP label calculation")]
    MaxIterationsExceeded,
    #[error("CIPLabeler unexpected up-edge ordering")]
    UnexpectedUpEdgeOrdering,
    #[error("No sequence rule provided")]
    NoSequenceRuleProvided,
    #[error("Descriptor lists should be the same length!")]
    DescriptorListLengthMismatch,
    #[error("Invalid stereo descriptor")]
    InvalidStereoDescriptor,
    #[error("Substituents should be topologically equivalent!")]
    SubstituentsShouldBeTopologicallyEquivalent,
    #[error("Something unexpected!")]
    SomethingUnexpected,
    #[error("Rule4b instance not in rule set")]
    Rule4bInstanceNotInRuleSet,
    #[error("Rule5New instance not in rule set")]
    Rule5NewInstanceNotInRuleSet,
    #[error("Parity vectors must have size 4.")]
    ParityVectorsMustHaveSize4,
    #[error("CIPLabeler Configuration requires at least one focus atom")]
    EmptyConfigurationFoci,
    #[error("CIPLabeler Tetrahedral received a bad atom")]
    BadTetrahedralAtom,
    #[error("CIPLabeler Tetrahedral received a bad config")]
    BadTetrahedralConfig,
    #[error("CIPLabeler Tetrahedral configuration must have 4 carriers")]
    TetrahedralConfigurationMustHave4Carriers,
    #[error("Received a Descriptor that is not supported for atoms")]
    DescriptorNotSupportedForAtoms,
    #[error("Received an invalid Atom Descriptor")]
    InvalidAtomDescriptor,
    #[error("Could not calculate parity! Carrier mismatch")]
    CarrierMismatch,
    #[error("CIPLabeler Sp2Bond received bad foci")]
    BadSp2BondFoci,
    #[error("CIPLabeler Sp2Bond received bad config")]
    BadSp2BondConfig,
    #[error("CIPLabeler Sp2Bond requires a double bond")]
    Sp2BondNotDoubleBond,
    #[error("CIPLabeler Sp2Bond requires defined bond stereo")]
    Sp2BondHasNoDefinedStereo,
    #[error("CIPLabeler Sp2Bond has incorrect number of stereo atoms")]
    IncorrectNumberOfStereoAtoms,
    #[error("Received a Descriptor that is not supported for double bonds")]
    DescriptorNotSupportedForDoubleBonds,
    #[error("Received an invalid Bond Descriptor")]
    InvalidBondDescriptor,
    #[error("CIPLabeler AtropisomerBond received bad foci")]
    BadAtropisomerBondFoci,
    #[error("CIPLabeler AtropisomerBond received bad config")]
    BadAtropisomerBondConfig,
    #[error("Received a Descriptor that is not supported for atropisomer bonds")]
    DescriptorNotSupportedForAtropisomerBonds,
    #[error(
        "RDKit-compatible CIPLabeler cancellation is unsupported because COSMolKit has no source-equivalent ControlCHandler"
    )]
    CancellationUnsupported,
    #[error("invalid CIPLabeler internal state: {detail}")]
    InvalidInternalState { detail: &'static str },
    #[error(transparent)]
    RingFinding(#[from] crate::RingFindingError),
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
    #[error("CIPLabeler could not construct its operation-local molecule view: {0}")]
    Invariant(#[from] crate::InvariantError),
}

pub(crate) fn descriptor_from_property(
    value: Option<&str>,
) -> Result<Option<CipDescriptor>, CipDescriptorError> {
    value.map(str::parse).transpose()
}
