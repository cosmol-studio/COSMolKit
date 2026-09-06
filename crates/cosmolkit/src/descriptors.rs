//! Descriptor methods on the canonical runtime molecule.

use crate::{Molecule, OperationError};

fn descriptor_error(
    operation: &'static str,
    error: cosmolkit_descriptors::DescriptorError,
) -> OperationError {
    OperationError::Algorithm {
        operation,
        detail: error.to_string(),
    }
}

impl Molecule {
    /// Returns the RDKit-compatible average molecular weight.
    #[must_use]
    pub fn molecular_weight(&self) -> Result<f64, OperationError> {
        self.molecular_weight_with_options(false)
    }

    /// Returns average molecular weight with an explicit heavy-atom mode.
    #[must_use]
    pub fn molecular_weight_with_options(&self, only_heavy: bool) -> Result<f64, OperationError> {
        cosmolkit_descriptors::molecular_weight_with_options(self.topology(), only_heavy)
            .map_err(|error| descriptor_error("molecular_weight", error))
    }

    /// Returns the RDKit-compatible exact molecular weight.
    #[must_use]
    pub fn exact_molecular_weight(&self) -> Result<f64, OperationError> {
        self.exact_molecular_weight_with_options(false)
    }

    /// Returns exact molecular weight with an explicit heavy-atom mode.
    #[must_use]
    pub fn exact_molecular_weight_with_options(
        &self,
        only_heavy: bool,
    ) -> Result<f64, OperationError> {
        cosmolkit_descriptors::exact_molecular_weight_with_options(self.topology(), only_heavy)
            .map_err(|error| descriptor_error("exact_molecular_weight", error))
    }

    /// Returns the Hill-ordered molecular formula.
    #[must_use]
    pub fn molecular_formula(&self) -> Result<String, OperationError> {
        self.molecular_formula_with_options(false, false)
    }

    /// Returns a molecular formula with isotope formatting controls.
    #[must_use]
    pub fn molecular_formula_with_options(
        &self,
        separate_isotopes: bool,
        abbreviate_h_isotopes: bool,
    ) -> Result<String, OperationError> {
        cosmolkit_descriptors::molecular_formula_with_options(
            self.topology(),
            separate_isotopes,
            abbreviate_h_isotopes,
        )
        .map_err(|error| descriptor_error("molecular_formula", error))
    }
}
