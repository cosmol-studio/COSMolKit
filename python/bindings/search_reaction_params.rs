use crate::chemical_values::{EnumValue, McsAtomComparator, McsBondComparator};
use cosmolkit as ck;
use pyo3::prelude::*;

#[derive(FromPyObject, IntoPyObject)]
pub enum McsAtomComparatorInput {
    Enum(McsAtomComparator),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum McsBondComparatorInput {
    Enum(McsBondComparator),
    Name(String),
}
impl McsAtomComparatorInput {
    pub(crate) fn into_rust(self) -> PyResult<ck::McsAtomComparator> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(value) => McsAtomComparator::from_name(&value).map(Into::into),
        }
    }
}
impl McsBondComparatorInput {
    pub(crate) fn into_rust(self) -> PyResult<ck::McsBondComparator> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(value) => McsBondComparator::from_name(&value).map(Into::into),
        }
    }
}

ck::project_configuration!(ck::McsParameters);

ck::project_configuration!(ck::ReactionSingleRunParams);

ck::project_configuration!(ck::ReactionWriteParams);
