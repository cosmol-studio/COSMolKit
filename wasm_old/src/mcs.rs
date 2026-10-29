//! Borrowed transport to the single public MCS implementation.
use crate::{Molecule, rust as ck};

pub fn maximum_common_substructure(inputs: &[&Molecule]) -> Result<ck::McsResult, ck::McsError> {
    maximum_common_substructure_with_params(inputs, &ck::McsParameters::default())
}

pub fn maximum_common_substructure_with_params(
    inputs: &[&Molecule],
    params: &ck::McsParameters,
) -> Result<ck::McsResult, ck::McsError> {
    let guards: Vec<_> = inputs.iter().map(|mol| mol.inner.borrow()).collect();
    let references: Vec<_> = guards.iter().map(|mol| &**mol).collect();
    ck::maximum_common_substructure_with_params(&references, params)
}
