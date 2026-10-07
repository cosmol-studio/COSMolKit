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
    LowerE,
    LowerZ,
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
        // RDKit✔️✔️:     case Descriptor::seqTrans:
        // RDKit✔️✔️:       return "e";
        // RDKit✔️✔️:     case Descriptor::seqCis:
        // RDKit✔️✔️:       return "z";
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
            Self::LowerE => "e",
            Self::LowerZ => "z",
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

impl CipDescriptor {
    fn from_bytes(value: &[u8]) -> Result<Self, CipDescriptorError> {
        match value {
            b"R" => Ok(Self::R),
            b"S" => Ok(Self::S),
            b"r" => Ok(Self::LowerR),
            b"s" => Ok(Self::LowerS),
            b"E" => Ok(Self::E),
            b"Z" => Ok(Self::Z),
            b"e" => Ok(Self::LowerE),
            b"z" => Ok(Self::LowerZ),
            b"M" => Ok(Self::M),
            b"P" => Ok(Self::P),
            b"m" => Ok(Self::LowerM),
            b"p" => Ok(Self::LowerP),
            _ => Err(CipDescriptorError::InvalidStoredDescriptor {
                value: value.into(),
            }),
        }
    }
}

impl FromStr for CipDescriptor {
    type Err = CipDescriptorError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        Self::from_bytes(value.as_bytes())
    }
}

/// Error returned when persisted `_CIPCode` state is not a descriptor emitted
/// by the supported modern assignment dispatcher.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CipDescriptorError {
    #[error("invalid stored CIP neighbor order {value:?}: {detail}")]
    InvalidNeighborOrder {
        value: crate::PropertyText,
        detail: String,
    },
    #[error(transparent)]
    Property(#[from] crate::PropertyValueError),
    #[error("invalid stored modern CIP descriptor {value:?}")]
    InvalidStoredDescriptor { value: crate::PropertyText },
}

/// Parse a stored `_CIPCode` property without depending on the CIP algorithm.
pub(crate) fn descriptor_from_property(
    value: Option<&crate::PropertyText>,
) -> Result<Option<CipDescriptor>, CipDescriptorError> {
    value
        .map(|value| CipDescriptor::from_bytes(value.as_bytes()))
        .transpose()
}

/// Decode the existing modern owner's stored neighbor-order representation.
/// Its JSON string serialization predates extraction; keep that exact format,
/// and also read canonical typed integer vectors without string conversion.
pub(crate) fn neighbor_order_from_property(
    value: Option<&crate::PropertyValue>,
) -> Result<Option<Vec<u32>>, CipDescriptorError> {
    value
        .map(|value| {
            if let crate::PropertyValue::String(text) = value {
                serde_json::from_slice(text.as_bytes()).map_err(|error| {
                    CipDescriptorError::InvalidNeighborOrder {
                        value: text.clone(),
                        detail: error.to_string(),
                    }
                })
            } else {
                value
                    .as_int_vector()?
                    .iter()
                    .map(|&index| {
                        u32::try_from(index).map_err(|error| {
                            CipDescriptorError::InvalidNeighborOrder {
                                value: index.to_string().into(),
                                detail: error.to_string(),
                            }
                        })
                    })
                    .collect()
            }
        })
        .transpose()
}
