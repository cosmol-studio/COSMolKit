//! Typed CIP descriptors stored on atom and bond properties.

use std::{fmt, str::FromStr};

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

/// Parse a stored `_CIPCode` property without depending on the CIP algorithm.
pub(crate) fn descriptor_from_property(
    value: Option<&str>,
) -> Result<Option<CipDescriptor>, CipDescriptorError> {
    value.map(str::parse).transpose()
}
