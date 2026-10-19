//! Typed enum-or-string extraction is native representation adaptation only.

ck::project_configuration!(ck::AromaticityParams);

ck::project_configuration!(ck::Coordinate2DInputParams);

ck::project_configuration!(ck::SdfReadParams);

use crate::chemical_values::{EnumValue, enum_value};
use cosmolkit as ck;
use pyo3::prelude::*;

// One binding declaration supplies the enum, exhaustive Rust conversion and
// string vocabulary. Input unions are ordinary official PyO3 type metadata.
macro_rules! selector {
    ($name:ident, $input:ident, {$($variant:ident),+ $(,)?}) => {
        enum_value!($name, { $($variant),+ });
        impl $input {
            pub(crate) fn into_rust(self) -> PyResult<ck::$name> {
                match self {
                    Self::Enum(value) => Ok(value.into()),
                    Self::Name(name) => $name::from_name(&name).map(Into::into),
                }
            }
        }
    };
}

// Keep official conversion derives on normal Rust items: expanding their
// tuple-variant derives inside macro_rules currently loses local identifier hygiene.
#[derive(FromPyObject, IntoPyObject)]
pub enum ValenceModelInput {
    Enum(ValenceModel),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum AromaticityModelInput {
    Enum(AromaticityModel),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum CoordinateZPolicyInput {
    Enum(CoordinateZPolicy),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum SdfCoordinateModeInput {
    Enum(SdfCoordinateMode),
    Name(String),
}

selector!(ValenceModel, ValenceModelInput, { RdkitLike });
selector!(AromaticityModel, AromaticityModelInput, { Rdkit, Simple, Mdl, Mmff94, Custom });
selector!(CoordinateZPolicy, CoordinateZPolicyInput, { Ignore, RequireZero, Error });
selector!(SdfCoordinateMode, SdfCoordinateModeInput, { Preserve, Require2D, Require3D });

fn default_valence_model() -> ValenceModelInput {
    ValenceModelInput::Enum(ck::ValenceParams::default().model.into())
}

ck::project_configuration!(ck::ValenceParams);
