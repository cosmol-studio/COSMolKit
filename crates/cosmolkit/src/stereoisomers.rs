//! Source-owned stereoisomer queries; live construction belongs to runtime.
pub use cosmolkit_stereo::{EnumerationError, StereoisomerOptions, StereoisomerRandomSource};
use std::{error::Error, fmt, sync::Arc};
#[derive(Debug, Clone)]
pub struct EnumerationRunError(Arc<EnumerationError>);
impl PartialEq for EnumerationRunError {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}
impl fmt::Display for EnumerationRunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self.0.as_ref(), f)
    }
}
impl Error for EnumerationRunError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(self.0.as_ref())
    }
}
impl From<EnumerationError> for crate::OperationError {
    fn from(error: EnumerationError) -> Self {
        Self::Enumeration(EnumerationRunError(Arc::new(error)))
    }
}
impl crate::Molecule {
    /// Return the source-defined upper bound using default enumeration options.
    pub fn stereoisomer_count(&self) -> Result<num_bigint::BigUint, EnumerationError> {
        // RDKit❗✔️: def GetStereoisomerCount(m, options=StereoEnumerationOptions()):
        // Default parameter construction and delegation have constant cost.
        self.stereoisomer_count_with_options(&StereoisomerOptions::default())
    }

    pub fn stereoisomer_count_with_options(
        &self,
        options: &StereoisomerOptions,
    ) -> Result<num_bigint::BigUint, EnumerationError> {
        // Original d892 stereoisomer_count delegates _getFlippers on a copy.
        // One explicit detached materialization; no live-state changes or commit.
        cosmolkit_stereo::stereoisomer_count(
            &cosmolkit_smiles::SmilesRecord {
                topology: self.topology().clone(),
                coordinates: self.coordinate_block_runtime().clone(),
                properties: self.properties().clone(),
            },
            options,
        )
    }
}
