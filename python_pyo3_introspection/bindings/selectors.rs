//! Typed enum-or-string extraction is native representation adaptation only.

#[pyclass(module = "fixed_projection._native")]
pub struct AromaticityParams {
    pub(crate) inner: ck::AromaticityParams,
}

#[pymethods]
impl AromaticityParams {
    #[new]
    #[pyo3(signature = (*, model=AromaticityModelInput::Enum(ck::AromaticityParams::default().model.into())))]
    fn new(model: AromaticityModelInput) -> PyResult<Self> {
        Ok(Self {
            inner: ck::AromaticityParams {
                model: model.into_rust()?,
            },
        })
    }
    #[getter]
    fn model(&self) -> AromaticityModel {
        self.inner.model.into()
    }
    #[setter]
    fn set_model(&mut self, value: AromaticityModelInput) -> PyResult<()> {
        self.inner.model = value.into_rust()?;
        Ok(())
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct Coordinate2DInputParams {
    pub(crate) inner: ck::Coordinate2DInputParams,
}

#[pymethods]
impl Coordinate2DInputParams {
    #[new]
    #[pyo3(signature = (*, z_policy=CoordinateZPolicyInput::Enum(ck::Coordinate2DInputParams::default().z_policy.into())))]
    fn new(z_policy: CoordinateZPolicyInput) -> PyResult<Self> {
        Ok(Self {
            inner: ck::Coordinate2DInputParams {
                z_policy: z_policy.into_rust()?,
            },
        })
    }
    #[getter]
    fn z_policy(&self) -> CoordinateZPolicy {
        self.inner.z_policy.into()
    }
    #[setter]
    fn set_z_policy(&mut self, value: CoordinateZPolicyInput) -> PyResult<()> {
        self.inner.z_policy = value.into_rust()?;
        Ok(())
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct SdfReadParams {
    pub(crate) inner: ck::SdfReadParams,
}

#[pymethods]
impl SdfReadParams {
    #[new]
    #[pyo3(signature = (*, sanitize=ck::SdfReadParams::default().sanitize, remove_hs=ck::SdfReadParams::default().remove_hs, strict_parsing=ck::SdfReadParams::default().strict_parsing, expand_attachment_points=ck::SdfReadParams::default().expand_attachment_points, process_property_lists=ck::SdfReadParams::default().process_property_lists, coordinate_mode=SdfCoordinateModeInput::Enum(ck::SdfReadParams::default().coordinate_mode.into())))]
    fn new(
        sanitize: bool,
        remove_hs: bool,
        strict_parsing: bool,
        expand_attachment_points: bool,
        process_property_lists: bool,
        coordinate_mode: SdfCoordinateModeInput,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::SdfReadParams {
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode: coordinate_mode.into_rust()?,
            },
        })
    }
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[setter]
    fn set_sanitize(&mut self, value: bool) -> PyResult<()> {
        self.inner.sanitize = value;
        Ok(())
    }
    #[getter]
    fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    #[setter]
    fn set_remove_hs(&mut self, value: bool) -> PyResult<()> {
        self.inner.remove_hs = value;
        Ok(())
    }
    #[getter]
    fn strict_parsing(&self) -> bool {
        self.inner.strict_parsing
    }
    #[setter]
    fn set_strict_parsing(&mut self, value: bool) -> PyResult<()> {
        self.inner.strict_parsing = value;
        Ok(())
    }
    #[getter]
    fn expand_attachment_points(&self) -> bool {
        self.inner.expand_attachment_points
    }
    #[setter]
    fn set_expand_attachment_points(&mut self, value: bool) -> PyResult<()> {
        self.inner.expand_attachment_points = value;
        Ok(())
    }
    #[getter]
    fn process_property_lists(&self) -> bool {
        self.inner.process_property_lists
    }
    #[setter]
    fn set_process_property_lists(&mut self, value: bool) -> PyResult<()> {
        self.inner.process_property_lists = value;
        Ok(())
    }
    #[getter]
    fn coordinate_mode(&self) -> SdfCoordinateMode {
        self.inner.coordinate_mode.into()
    }
    #[setter]
    fn set_coordinate_mode(&mut self, value: SdfCoordinateModeInput) -> PyResult<()> {
        self.inner.coordinate_mode = value.into_rust()?;
        Ok(())
    }
}

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

#[pyclass(module = "fixed_projection._native")]
pub struct ValenceParams {
    pub(crate) inner: ck::ValenceParams,
}

#[pymethods]
impl ValenceParams {
    #[new]
    #[pyo3(signature = (*, model=default_valence_model(), strict=ck::ValenceParams::default().strict))]
    fn new(model: ValenceModelInput, strict: bool) -> PyResult<Self> {
        Ok(Self {
            inner: ck::ValenceParams {
                model: model.into_rust()?,
                strict,
            },
        })
    }
    #[getter]
    fn model(&self) -> ValenceModel {
        self.inner.model.into()
    }
    #[setter]
    fn set_model(&mut self, value: ValenceModelInput) -> PyResult<()> {
        self.inner.model = value.into_rust()?;
        Ok(())
    }
    #[getter]
    fn strict(&self) -> bool {
        self.inner.strict
    }
    #[setter]
    fn set_strict(&mut self, value: bool) {
        self.inner.strict = value;
    }
}
