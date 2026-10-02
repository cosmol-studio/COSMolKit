//! Thin public SMILES construction over the detached parser and chemistry owners.

use std::fmt;

use crate::{Molecule, OperationError};

/// Original-index fragment selection; symbol arrays are indexed by the full molecule.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct FragmentSmilesWriteParams {
    pub smiles: crate::SmilesWriteParams,
    pub atoms: Vec<crate::AtomId>,
    /// None selects all bonds between selected atoms; Some(empty) selects no bonds.
    pub bonds: Option<Vec<crate::BondId>>,
    pub atom_symbols: Option<Vec<String>>,
    pub bond_symbols: Option<Vec<String>>,
}

/// Fragment selection with explicit CX fields and dimension-scoped coordinates.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct FragmentCxSmilesWriteParams {
    pub cx: crate::CxSmilesWriteParams,
    pub atoms: Vec<crate::AtomId>,
    pub bonds: Option<Vec<crate::BondId>>,
    pub atom_symbols: Option<Vec<String>>,
    pub bond_symbols: Option<Vec<String>>,
}

/// Source-preserving serialization failures, without converting causes to strings.
#[derive(Debug)]
pub enum SmilesWriteError {
    Write(cosmolkit_smiles::SmilesParseError),
    Fragment(cosmolkit_smiles::FragmentWriteInputError),
}

impl fmt::Display for SmilesWriteError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Write(error) => write!(formatter, "SMILES serialization failed: {error}"),
            Self::Fragment(error) => {
                write!(formatter, "fragment SMILES serialization failed: {error}")
            }
        }
    }
}

impl std::error::Error for SmilesWriteError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Write(error) => Some(error),
            Self::Fragment(error) => Some(error),
        }
    }
}

impl From<cosmolkit_smiles::SmilesParseError> for SmilesWriteError {
    fn from(error: cosmolkit_smiles::SmilesParseError) -> Self {
        Self::Write(error)
    }
}

impl From<cosmolkit_smiles::FragmentWriteInputError> for SmilesWriteError {
    fn from(error: cosmolkit_smiles::FragmentWriteInputError) -> Self {
        Self::Fragment(error)
    }
}

impl Molecule {
    fn smiles_ring_state(&self) -> Option<&cosmolkit_core::RingInfo> {
        #[cfg(feature = "cap-rings")]
        {
            self.derived_cache_runtime().ring_info()
        }
        #[cfg(not(feature = "cap-rings"))]
        {
            None
        }
    }

    fn smiles_valence_state(&self) -> Option<&cosmolkit_core::ValenceAssignment> {
        self.derived_cache_runtime().valence_assignment()
    }

    fn smiles_view(&self) -> cosmolkit_smiles::SmilesRecordView<'_> {
        cosmolkit_smiles::SmilesRecordView {
            topology: self.topology(),
            coordinates: self.coordinate_block_runtime(),
            properties: self.properties(),
        }
    }

    /// Serialize without changing molecule state or installing writer caches.
    pub fn to_smiles(&self) -> Result<String, SmilesWriteError> {
        self.to_smiles_with_params(&crate::SmilesWriteParams::default())
    }

    pub fn to_smiles_with_params(
        &self,
        params: &crate::SmilesWriteParams,
    ) -> Result<String, SmilesWriteError> {
        Ok(cosmolkit_smiles::write_smiles_with_params(
            self.smiles_view(),
            params,
        )?)
    }

    pub fn to_cx_smiles(&self) -> Result<String, SmilesWriteError> {
        self.to_cx_smiles_with_params(&crate::CxSmilesWriteParams::default())
    }

    /// Auto coordinates require a unique stored set; explicit dimension/ID resolves ambiguity.
    pub fn to_cx_smiles_with_params(
        &self,
        params: &crate::CxSmilesWriteParams,
    ) -> Result<String, SmilesWriteError> {
        Ok(cosmolkit_smiles::write_cx_smiles_with_params(
            self.smiles_view(),
            params,
        )?)
    }

    pub fn to_fragment_smiles(&self, atoms: &[crate::AtomId]) -> Result<String, SmilesWriteError> {
        Ok(cosmolkit_smiles::write_fragment_smiles_output(
            self.smiles_view(),
            &crate::SmilesWriteParams::default(),
            atoms,
            None,
            None,
            None,
            self.smiles_ring_state(),
            self.smiles_valence_state(),
        )?
        .text)
    }

    pub fn to_fragment_smiles_with_params(
        &self,
        params: &FragmentSmilesWriteParams,
    ) -> Result<String, SmilesWriteError> {
        Ok(cosmolkit_smiles::write_fragment_smiles_output(
            self.smiles_view(),
            &params.smiles,
            &params.atoms,
            params.bonds.as_deref(),
            params.atom_symbols.as_deref(),
            params.bond_symbols.as_deref(),
            self.smiles_ring_state(),
            self.smiles_valence_state(),
        )?
        .text)
    }

    pub fn to_fragment_cx_smiles(
        &self,
        atoms: &[crate::AtomId],
    ) -> Result<String, SmilesWriteError> {
        Ok(cosmolkit_smiles::write_fragment_cx_smiles(
            self.smiles_view(),
            &crate::CxSmilesWriteParams::default(),
            atoms,
            None,
            None,
            None,
            self.smiles_ring_state(),
            self.smiles_valence_state(),
        )?)
    }

    pub fn to_fragment_cx_smiles_with_params(
        &self,
        params: &FragmentCxSmilesWriteParams,
    ) -> Result<String, SmilesWriteError> {
        Ok(cosmolkit_smiles::write_fragment_cx_smiles(
            self.smiles_view(),
            &params.cx,
            &params.atoms,
            params.bonds.as_deref(),
            params.atom_symbols.as_deref(),
            params.bond_symbols.as_deref(),
            self.smiles_ring_state(),
            self.smiles_valence_state(),
        )?)
    }

    /// Preserve source ordering and duplicates. Seeds 1..=i32::MAX reseed;
    /// zero and high-bit u32 seeds continue the shared stream, matching the
    /// source's u32-to-i32 cast and positive-only reseeding condition.
    pub fn to_random_smiles(&self, count: u32, seed: u32) -> Result<Vec<String>, SmilesWriteError> {
        self.to_random_smiles_with_params(count, seed, &crate::RandomSmilesWriteParams::default())
    }

    pub fn to_random_smiles_with_params(
        &self,
        count: u32,
        seed: u32,
        params: &crate::RandomSmilesWriteParams,
    ) -> Result<Vec<String>, SmilesWriteError> {
        Ok(cosmolkit_smiles::write_random_smiles_vector(
            self.smiles_view(),
            count,
            seed,
            params,
        )?)
    }
}

/// Structured failure from the public SMILES construction pipeline.
#[derive(Debug)]
pub enum SmilesError {
    /// SMILES grammar, CXSMILES, or detached parser validation failed.
    Parse(cosmolkit_smiles::SmilesParseError),
    /// Source-defined hydrogen removal failed.
    Hydrogen(cosmolkit_core::HydrogenError),
    /// Source-defined sanitization failed.
    Sanitize(cosmolkit_core::SanitizeError),
    /// Source-defined post-chemistry stereo completion failed.
    Stereo(cosmolkit_smiles::SmilesStereoError),
    /// Final live-state validation failed.
    Construction(OperationError),
}

impl fmt::Display for SmilesError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(error) => write!(formatter, "SMILES parsing failed: {error}"),
            Self::Hydrogen(error) => write!(formatter, "SMILES hydrogen removal failed: {error}"),
            Self::Sanitize(error) => write!(formatter, "SMILES sanitization failed: {error}"),
            Self::Stereo(error) => write!(formatter, "SMILES stereo completion failed: {error}"),
            Self::Construction(error) => {
                write!(formatter, "SMILES molecule construction failed: {error}")
            }
        }
    }
}

impl std::error::Error for SmilesError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Parse(error) => Some(error),
            Self::Hydrogen(error) => Some(error),
            Self::Sanitize(error) => Some(error),
            Self::Stereo(error) => Some(error),
            Self::Construction(error) => Some(error),
        }
    }
}

impl From<cosmolkit_smiles::SmilesParseError> for SmilesError {
    fn from(error: cosmolkit_smiles::SmilesParseError) -> Self {
        Self::Parse(error)
    }
}

impl From<cosmolkit_core::HydrogenError> for SmilesError {
    fn from(error: cosmolkit_core::HydrogenError) -> Self {
        Self::Hydrogen(error)
    }
}

impl From<cosmolkit_core::SanitizeError> for SmilesError {
    fn from(error: cosmolkit_core::SanitizeError) -> Self {
        Self::Sanitize(error)
    }
}

impl From<OperationError> for SmilesError {
    fn from(error: OperationError) -> Self {
        Self::Construction(error)
    }
}

impl From<cosmolkit_smiles::SmilesStereoError> for SmilesError {
    fn from(error: cosmolkit_smiles::SmilesStereoError) -> Self {
        Self::Stereo(error)
    }
}

impl Molecule {
    /// Constructs a molecule from SMILES with the pinned source defaults.
    pub fn from_smiles(input: &str) -> Result<Self, SmilesError> {
        Self::from_smiles_with_params(input, &cosmolkit_smiles::SmilesParseParams::default())
    }

    /// Constructs a molecule from SMILES using explicit parser parameters.
    ///
    /// Parsing, sanitization, and hydrogen removal remain in their detached
    /// algorithm owners. A live molecule is installed only after every
    /// requested stage and final structural validation succeeds.
    pub fn from_smiles_with_params(
        input: &str,
        params: &cosmolkit_smiles::SmilesParseParams,
    ) -> Result<Self, SmilesError> {
        let cosmolkit_smiles::SmilesRecord {
            mut topology,
            mut coordinates,
            mut properties,
        } = cosmolkit_smiles::parse_smiles(input, params)?;

        let mut final_valence = None;
        if params.remove_hydrogens {
            let remove_params = cosmolkit_core::RemoveHsParams {
                update_explicit_count: true,
                sanitize: params.sanitize,
                ..cosmolkit_core::RemoveHsParams::default()
            };
            let result = cosmolkit_core::remove_hydrogens_with_params(
                topology,
                coordinates,
                properties,
                &remove_params,
            )?;
            topology = result.topology;
            coordinates = result.coordinates;
            properties = result.properties;
            final_valence = result.final_valence;
        } else if params.sanitize {
            let result = cosmolkit_core::sanitize_topology(
                &topology,
                &cosmolkit_core::SanitizeParams::default(),
            )?;
            topology = result.topology;
            final_valence = result.final_valence;
        }

        let record = cosmolkit_smiles::finalize_smiles_stereo(
            cosmolkit_smiles::SmilesRecord {
                topology,
                coordinates,
                properties,
            },
            params,
            &mut final_valence,
        )?;
        // Stereo may need local non-strict values, but sanitize=false must not
        // revive runtime validity after CK-VALENCE-001 hydrogen removal.
        Self::from_smiles_parts_with_valence(
            record.topology,
            record.coordinates,
            record.properties,
            if params.sanitize { final_valence } else { None },
        )
        .map_err(SmilesError::Construction)
    }
}
