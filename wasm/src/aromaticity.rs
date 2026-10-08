//! Thin projections of the four public aromaticity operations.
use crate::Molecule;
use cosmolkit::{AromaticityParams, OperationError};

impl Molecule {
    pub fn with_assigned_aromaticity(&self) -> Result<Self, OperationError> {
        self.inner
            .borrow()
            .with_assigned_aromaticity()
            .map(|inner| Self {
                inner: inner.into(),
            })
    }

    pub fn with_assigned_aromaticity_with_params(
        &self,
        params: &AromaticityParams,
    ) -> Result<Self, OperationError> {
        self.inner
            .borrow()
            .with_assigned_aromaticity_with_params(params)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }

    pub fn assign_aromaticity_(&self) -> Result<(), OperationError> {
        self.inner.borrow_mut().assign_aromaticity_()
    }

    pub fn assign_aromaticity_with_params_(
        &self,
        params: &AromaticityParams,
    ) -> Result<(), OperationError> {
        self.inner
            .borrow_mut()
            .assign_aromaticity_with_params_(params)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit::{AromaticityModel, Molecule as CoreMolecule, SmilesParseParams};

    #[test]
    fn aromaticity_projection_matches_public_facade_and_preserves_failed_mutations() {
        let raw = CoreMolecule::from_smiles_with_params(
            "C1=CC=CC=C1",
            &SmilesParseParams {
                sanitize: false,
                ..Default::default()
            },
        )
        .unwrap();
        let bound = Molecule {
            inner: raw.clone().into(),
        };
        assert_eq!(
            *bound.with_assigned_aromaticity().unwrap().inner.borrow(),
            raw.with_assigned_aromaticity().unwrap()
        );
        assert_eq!(*bound.inner.borrow(), raw);
        bound.assign_aromaticity_().unwrap();
        assert_eq!(
            *bound.inner.borrow(),
            raw.with_assigned_aromaticity().unwrap()
        );
        for model in [
            AromaticityModel::Rdkit,
            AromaticityModel::Simple,
            AromaticityModel::Mdl,
            AromaticityModel::Mmff94,
            AromaticityModel::Custom,
        ] {
            let params = AromaticityParams { model };
            let bound = Molecule {
                inner: raw.clone().into(),
            };
            let expected = raw.with_assigned_aromaticity_with_params(&params);
            let actual = bound.with_assigned_aromaticity_with_params(&params);
            match (expected, actual) {
                (Ok(expected), Ok(actual)) => {
                    assert_eq!(*actual.inner.borrow(), expected);
                    assert_eq!(*bound.inner.borrow(), raw);
                    bound.assign_aromaticity_with_params_(&params).unwrap();
                    assert_eq!(*bound.inner.borrow(), expected);
                }
                (Err(expected), Err(actual)) => {
                    assert_eq!(actual.to_string(), expected.to_string());
                    assert!(bound.assign_aromaticity_with_params_(&params).is_err());
                    assert_eq!(*bound.inner.borrow(), raw);
                }
                _ => panic!("binding changed public-facade outcome for {model:?}"),
            }
        }
    }
}
