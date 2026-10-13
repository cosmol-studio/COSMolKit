//! Persistent evaluators own detached native state; bindings do not write to molecules.
use crate::{
    Molecule,
    host_values::{bool_value, point, sequence, type_error, u32_value, usize_value},
};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::error::Error;
use wasm_bindgen::prelude::*;
fn number(value: &JsValue, name: &str) -> Result<f64, JsValue> {
    value.as_f64().ok_or_else(|| type_error(name))
}
fn optional_index(value: &JsValue, name: &str) -> Result<Option<usize>, JsValue> {
    if value.is_null() || value.is_undefined() {
        Ok(None)
    } else {
        usize_value(value, name).map(Some)
    }
}
fn rows(values: &[[f64; 3]]) -> Array {
    values
        .iter()
        .map(|row| row.iter().copied().map(JsValue::from).collect::<Array>())
        .collect()
}
#[wasm_bindgen]
pub struct MmffForceFieldParams {
    inner: ck::MmffForceFieldParams,
}
#[cosmolkit_wasm::javascript_options("MmffForceFieldParamsOptions")]
#[wasm_bindgen]
impl MmffForceFieldParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string")] mmff_variant: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] non_bonded_threshold: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        ignore_interfragment_interactions: JsValue,
    ) -> Result<Self, JsValue> {
        let defaults = ck::MmffForceFieldParams::default();
        let conformer_id = if conformer_id.is_undefined() {
            defaults.conformer_id()
        } else {
            optional_index(&conformer_id, "conformerId")?
        };
        let mmff_variant = if mmff_variant.is_undefined() {
            defaults.mmff_variant().to_owned()
        } else {
            mmff_variant
                .as_string()
                .ok_or_else(|| type_error("mmffVariant"))?
        };
        let non_bonded_threshold = if non_bonded_threshold.is_undefined() {
            defaults.non_bonded_threshold()
        } else {
            number(&non_bonded_threshold, "nonBondedThreshold")?
        };
        let ignore_interfragment_interactions = if ignore_interfragment_interactions.is_undefined()
        {
            defaults.ignore_interfragment_interactions()
        } else {
            bool_value(
                &ignore_interfragment_interactions,
                "ignoreInterfragmentInteractions",
            )?
        };
        Ok(Self {
            inner: ck::MmffForceFieldParams::new(
                conformer_id,
                mmff_variant,
                non_bonded_threshold,
                ignore_interfragment_interactions,
            ),
        })
    }
    #[wasm_bindgen(getter, js_name = conformerId, unchecked_return_type = "number | null")]
    pub fn conformer_id(&self) -> JsValue {
        self.inner
            .conformer_id()
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
    #[wasm_bindgen(getter, js_name = mmffVariant)]
    pub fn mmff_variant(&self) -> String {
        self.inner.mmff_variant().to_owned()
    }
    #[wasm_bindgen(getter, js_name = nonBondedThreshold)]
    pub fn non_bonded_threshold(&self) -> f64 {
        self.inner.non_bonded_threshold()
    }
    #[wasm_bindgen(getter, js_name = ignoreInterfragmentInteractions)]
    pub fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions()
    }
}
#[wasm_bindgen]
pub struct UffForceFieldParams {
    inner: ck::UffForceFieldParams,
}
#[cosmolkit_wasm::javascript_options("UffForceFieldParamsOptions")]
#[wasm_bindgen]
impl UffForceFieldParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] vdw_threshold: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        ignore_interfragment_interactions: JsValue,
    ) -> Result<Self, JsValue> {
        let defaults = ck::UffForceFieldParams::default();
        let conformer_id = if conformer_id.is_undefined() {
            defaults.conformer_id()
        } else {
            optional_index(&conformer_id, "conformerId")?
        };
        let vdw_threshold = if vdw_threshold.is_undefined() {
            defaults.vdw_threshold()
        } else {
            number(&vdw_threshold, "vdwThreshold")?
        };
        let ignore_interfragment_interactions = if ignore_interfragment_interactions.is_undefined()
        {
            defaults.ignore_interfragment_interactions()
        } else {
            bool_value(
                &ignore_interfragment_interactions,
                "ignoreInterfragmentInteractions",
            )?
        };
        Ok(Self {
            inner: ck::UffForceFieldParams::new(
                conformer_id,
                vdw_threshold,
                ignore_interfragment_interactions,
            ),
        })
    }
    #[wasm_bindgen(getter, js_name = conformerId, unchecked_return_type = "number | null")]
    pub fn conformer_id(&self) -> JsValue {
        self.inner
            .conformer_id()
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
    #[wasm_bindgen(getter, js_name = vdwThreshold)]
    pub fn vdw_threshold(&self) -> f64 {
        self.inner.vdw_threshold()
    }
    #[wasm_bindgen(getter, js_name = ignoreInterfragmentInteractions)]
    pub fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions()
    }
}
#[wasm_bindgen]
pub struct ForceFieldMinimizeParams {
    inner: ck::ForceFieldMinimizeParams,
}
#[cosmolkit_wasm::javascript_options("ForceFieldMinimizeParamsOptions")]
#[wasm_bindgen]
impl ForceFieldMinimizeParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] force_tolerance: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] energy_tolerance: JsValue,
    ) -> Result<Self, JsValue> {
        let defaults = ck::ForceFieldMinimizeParams::default();
        let max_iterations = if max_iterations.is_undefined() {
            defaults.max_iterations()
        } else {
            u32_value(&max_iterations, "maxIterations")?
        };
        let force_tolerance = if force_tolerance.is_undefined() {
            defaults.force_tolerance()
        } else {
            number(&force_tolerance, "forceTolerance")?
        };
        let energy_tolerance = if energy_tolerance.is_undefined() {
            defaults.energy_tolerance()
        } else {
            number(&energy_tolerance, "energyTolerance")?
        };
        Ok(Self {
            inner: ck::ForceFieldMinimizeParams::new(
                max_iterations,
                force_tolerance,
                energy_tolerance,
            ),
        })
    }
    #[wasm_bindgen(getter, js_name = maxIterations)]
    pub fn max_iterations(&self) -> u32 {
        self.inner.max_iterations()
    }
    #[wasm_bindgen(getter, js_name = forceTolerance)]
    pub fn force_tolerance(&self) -> f64 {
        self.inner.force_tolerance()
    }
    #[wasm_bindgen(getter, js_name = energyTolerance)]
    pub fn energy_tolerance(&self) -> f64 {
        self.inner.energy_tolerance()
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum MolecularForceFieldErrorKind {
    MissingConformer = 0,
    Parameterization = 1,
    InvalidParameterization = 2,
    InvalidAtomIndex = 3,
    InvalidFixedAtom = 4,
    CoordinateCount = 5,
    CoordinateShape = 6,
    NonFiniteCoordinate = 7,
    InvalidTolerance = 8,
    Preparation = 9,
    Construction = 10,
    Initialization = 11,
    Rings = 12,
    Kernel = 13,
}
impl From<ck::MolecularForceFieldErrorKind> for MolecularForceFieldErrorKind {
    fn from(value: ck::MolecularForceFieldErrorKind) -> Self {
        match value {
            ck::MolecularForceFieldErrorKind::MissingConformer => Self::MissingConformer,
            ck::MolecularForceFieldErrorKind::Parameterization => Self::Parameterization,
            ck::MolecularForceFieldErrorKind::InvalidParameterization => {
                Self::InvalidParameterization
            }
            ck::MolecularForceFieldErrorKind::InvalidAtomIndex => Self::InvalidAtomIndex,
            ck::MolecularForceFieldErrorKind::InvalidFixedAtom => Self::InvalidFixedAtom,
            ck::MolecularForceFieldErrorKind::CoordinateCount => Self::CoordinateCount,
            ck::MolecularForceFieldErrorKind::CoordinateShape => Self::CoordinateShape,
            ck::MolecularForceFieldErrorKind::NonFiniteCoordinate => Self::NonFiniteCoordinate,
            ck::MolecularForceFieldErrorKind::InvalidTolerance => Self::InvalidTolerance,
            ck::MolecularForceFieldErrorKind::Preparation => Self::Preparation,
            ck::MolecularForceFieldErrorKind::Construction => Self::Construction,
            ck::MolecularForceFieldErrorKind::Initialization => Self::Initialization,
            ck::MolecularForceFieldErrorKind::Rings => Self::Rings,
            ck::MolecularForceFieldErrorKind::Kernel => Self::Kernel,
        }
    }
}
#[wasm_bindgen]
pub struct ForceFieldError {
    kind: MolecularForceFieldErrorKind,
    message: String,
    cause: JsValue,
    requested: Option<usize>,
    atom_index: Option<usize>,
    component: Option<usize>,
    actual: Option<usize>,
    expected: Option<usize>,
}
#[wasm_bindgen]
impl ForceFieldError {
    pub fn kind(&self) -> MolecularForceFieldErrorKind {
        self.kind
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(js_name = requested, unchecked_return_type = "number | null")]
    pub fn requested(&self) -> JsValue {
        self.requested
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
    #[wasm_bindgen(js_name = atomIndex, unchecked_return_type = "number | null")]
    pub fn atom_index(&self) -> JsValue {
        self.atom_index
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
    #[wasm_bindgen(js_name = component, unchecked_return_type = "number | null")]
    pub fn component(&self) -> JsValue {
        self.component
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
    #[wasm_bindgen(js_name = actual, unchecked_return_type = "number | null")]
    pub fn actual(&self) -> JsValue {
        self.actual
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
    #[wasm_bindgen(js_name = expected, unchecked_return_type = "number | null")]
    pub fn expected(&self) -> JsValue {
        self.expected
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
}
pub(crate) fn field_error(source: &ck::ForceFieldError) -> Result<JsValue, JsValue> {
    let detail = ForceFieldError {
        kind: source.kind().into(),
        message: source.to_string(),
        cause: match source.source() {
            Some(cause) => {
                if let Some(field) = cause.downcast_ref::<ck::ForceFieldError>() {
                    field_error(field)?
                } else {
                    crate::alignment_values::source_error(cause)?
                }
            }
            None => JsValue::NULL,
        },
        requested: source.requested(),
        atom_index: source.atom_index(),
        component: source.component(),
        actual: source.actual(),
        expected: source.expected(),
    };
    let error: JsValue = js_sys::Error::new(&detail.message).into();
    use crate::alignment_values::set;
    set(&error, "name", "ForceFieldError".into())?;
    set(&error, "domain", "forcefields".into())?;
    set(&error, "kind", (detail.kind as u32).into())?;
    set(&error, "cause", detail.cause.clone())?;
    set(&error, "requested", detail.requested())?;
    set(&error, "atomIndex", detail.atom_index())?;
    set(&error, "component", detail.component())?;
    set(&error, "actual", detail.actual())?;
    set(&error, "expected", detail.expected())?;
    set(&error, "detail", detail.into())?;
    Ok(error)
}
#[wasm_bindgen]
pub struct MmffForceFieldError {
    kind: MolecularForceFieldErrorKind,
    message: String,
    cause: JsValue,
    requested: Option<usize>,
    atom_index: Option<usize>,
    component: Option<usize>,
    actual: Option<usize>,
    expected: Option<usize>,
}
#[wasm_bindgen]
impl MmffForceFieldError {
    pub fn kind(&self) -> MolecularForceFieldErrorKind {
        self.kind
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(js_name = requested, unchecked_return_type = "number | null")]
    pub fn requested(&self) -> JsValue {
        self.requested
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
    #[wasm_bindgen(js_name = atomIndex, unchecked_return_type = "number | null")]
    pub fn atom_index(&self) -> JsValue {
        self.atom_index
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
    #[wasm_bindgen(js_name = component, unchecked_return_type = "number | null")]
    pub fn component(&self) -> JsValue {
        self.component
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
    #[wasm_bindgen(js_name = actual, unchecked_return_type = "number | null")]
    pub fn actual(&self) -> JsValue {
        self.actual
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
    #[wasm_bindgen(js_name = expected, unchecked_return_type = "number | null")]
    pub fn expected(&self) -> JsValue {
        self.expected
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
}
pub(crate) fn mmff_error(source: &ck::MmffForceFieldError) -> Result<JsValue, JsValue> {
    let detail = MmffForceFieldError {
        kind: source.kind().into(),
        message: source.to_string(),
        cause: match source.source() {
            Some(cause) => {
                if let Some(field) = cause.downcast_ref::<ck::ForceFieldError>() {
                    field_error(field)?
                } else {
                    crate::alignment_values::source_error(cause)?
                }
            }
            None => JsValue::NULL,
        },
        requested: source.requested(),
        atom_index: source.atom_index(),
        component: source.component(),
        actual: source.actual(),
        expected: source.expected(),
    };
    let error: JsValue = js_sys::Error::new(&detail.message).into();
    use crate::alignment_values::set;
    set(&error, "name", "MmffForceFieldError".into())?;
    set(&error, "domain", "forcefields".into())?;
    set(&error, "kind", (detail.kind as u32).into())?;
    set(&error, "cause", detail.cause.clone())?;
    set(&error, "requested", detail.requested())?;
    set(&error, "atomIndex", detail.atom_index())?;
    set(&error, "component", detail.component())?;
    set(&error, "actual", detail.actual())?;
    set(&error, "expected", detail.expected())?;
    set(&error, "detail", detail.into())?;
    Ok(error)
}
#[wasm_bindgen]
pub struct UffForceFieldError {
    kind: MolecularForceFieldErrorKind,
    message: String,
    cause: JsValue,
    requested: Option<usize>,
}
#[wasm_bindgen]
impl UffForceFieldError {
    pub fn kind(&self) -> MolecularForceFieldErrorKind {
        self.kind
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(js_name = requested, unchecked_return_type = "number | null")]
    pub fn requested(&self) -> JsValue {
        self.requested
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
}
pub(crate) fn uff_error(source: &ck::UffForceFieldError) -> Result<JsValue, JsValue> {
    let detail = UffForceFieldError {
        kind: source.kind().into(),
        message: source.to_string(),
        cause: match source.source() {
            Some(cause) => {
                if let Some(field) = cause.downcast_ref::<ck::ForceFieldError>() {
                    field_error(field)?
                } else {
                    crate::alignment_values::source_error(cause)?
                }
            }
            None => JsValue::NULL,
        },
        requested: source.requested(),
    };
    let error: JsValue = js_sys::Error::new(&detail.message).into();
    use crate::alignment_values::set;
    set(&error, "name", "UffForceFieldError".into())?;
    set(&error, "domain", "forcefields".into())?;
    set(&error, "kind", (detail.kind as u32).into())?;
    set(&error, "cause", detail.cause.clone())?;
    set(&error, "requested", detail.requested())?;
    set(&error, "detail", detail.into())?;
    Ok(error)
}

#[wasm_bindgen]
pub struct ForceFieldEnergyGradient {
    inner: ck::ForceFieldEnergyGradient,
}
#[wasm_bindgen]
impl ForceFieldEnergyGradient {
    pub fn energy(&self) -> f64 {
        self.inner.energy()
    }
    #[wasm_bindgen(unchecked_return_type = "number[][]")]
    pub fn gradient(&self) -> Array {
        rows(self.inner.gradient())
    }
}
#[wasm_bindgen]
pub struct ForceFieldMinimizeOutcome {
    inner: ck::ForceFieldMinimizeOutcome,
}
#[wasm_bindgen]
impl ForceFieldMinimizeOutcome {
    pub fn converged(&self) -> bool {
        self.inner.converged()
    }
    pub fn iterations(&self) -> u32 {
        self.inner.iterations()
    }
    pub fn energy(&self) -> f64 {
        self.inner.energy()
    }
}
#[wasm_bindgen]
pub struct MolecularForceField {
    inner: ck::MolecularForceField,
}
#[wasm_bindgen]
impl MolecularForceField {
    #[wasm_bindgen(unchecked_return_type = "number[]")]
    pub fn position(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom_id: JsValue,
    ) -> Result<Array, JsValue> {
        self.inner
            .position(ck::AtomId::new(usize_value(&atom_id, "atomId")?))
            .map(|v| v.into_iter().map(JsValue::from).collect())
            .map_err(|e| field_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = setPosition_)]
    pub fn set_position_(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom_id: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array")] position: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_position_(
                ck::AtomId::new(usize_value(&atom_id, "atomId")?),
                point(&position)?,
            )
            .map_err(|e| field_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(unchecked_return_type = "number[][]")]
    pub fn positions(&self) -> Array {
        rows(&self.inner.positions())
    }
    #[wasm_bindgen(js_name = setPositions_)]
    pub fn set_positions_(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "(number[] | Float64Array)[]")] positions: JsValue,
    ) -> Result<(), JsValue> {
        let positions = sequence(&positions, "positions")?
            .iter()
            .map(|v| point(&v))
            .collect::<Result<Vec<_>, _>>()?;
        self.inner
            .set_positions_(&positions)
            .map_err(|e| field_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = fixedAtoms, unchecked_return_type = "number[]")]
    pub fn fixed_atoms(&self) -> Array {
        self.inner
            .fixed_atoms()
            .into_iter()
            .map(|id| JsValue::from_f64(id.index() as f64))
            .collect()
    }
    #[wasm_bindgen(js_name = setFixedAtoms_)]
    pub fn set_fixed_atoms_(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Uint32Array")] atoms: JsValue,
    ) -> Result<(), JsValue> {
        let atoms = sequence(&atoms, "atoms")?
            .iter()
            .map(|v| usize_value(&v, "atom").map(ck::AtomId::new))
            .collect::<Result<Vec<_>, _>>()?;
        self.inner
            .set_fixed_atoms_(&atoms)
            .map_err(|e| field_error(&e).unwrap_or_else(|e| e))
    }
    pub fn energy(&self) -> Result<f64, JsValue> {
        self.inner
            .energy()
            .map_err(|e| field_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(unchecked_return_type = "number[][]")]
    pub fn gradient(&self) -> Result<Array, JsValue> {
        self.inner
            .gradient()
            .map(|v| rows(&v))
            .map_err(|e| field_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = gradientUnconstrained, unchecked_return_type = "number[][]")]
    pub fn gradient_unconstrained(&self) -> Result<Array, JsValue> {
        self.inner
            .gradient_unconstrained()
            .map(|v| rows(&v))
            .map_err(|e| field_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = energyGradient)]
    pub fn energy_gradient(&self) -> Result<ForceFieldEnergyGradient, JsValue> {
        self.inner
            .energy_gradient()
            .map(|inner| ForceFieldEnergyGradient { inner })
            .map_err(|e| field_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = minimize_)]
    pub fn minimize_(
        &mut self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "ForceFieldMinimizeParams | ForceFieldMinimizeParamsOptions"
        )]
        params: JsValue,
    ) -> Result<ForceFieldMinimizeOutcome, JsValue> {
        let params = ForceFieldMinimizeParams::from_js_configuration(&params)?;
        self.minimize_with_params_(&params)
    }
    #[wasm_bindgen(js_name = minimizeWithParams_)]
    pub fn minimize_with_params_(
        &mut self,
        params: &ForceFieldMinimizeParams,
    ) -> Result<ForceFieldMinimizeOutcome, JsValue> {
        self.inner
            .minimize_with_params_(&params.inner)
            .map(|inner| ForceFieldMinimizeOutcome { inner })
            .map_err(|e| field_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = mmffForceField)]
    pub fn mmff_force_field(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "MmffForceFieldParams | MmffForceFieldParamsOptions"
        )]
        params: JsValue,
    ) -> Result<MolecularForceField, JsValue> {
        let params = MmffForceFieldParams::from_js_configuration(&params)?;
        self.mmff_force_field_with_params(&params)
    }
    #[wasm_bindgen(js_name = mmffForceFieldWithParams)]
    pub fn mmff_force_field_with_params(
        &self,
        params: &MmffForceFieldParams,
    ) -> Result<MolecularForceField, JsValue> {
        self.inner
            .mmff_force_field_with_params(&params.inner)
            .map(|inner| MolecularForceField { inner })
            .map_err(|e| mmff_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = uffForceField)]
    pub fn uff_force_field(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "UffForceFieldParams | UffForceFieldParamsOptions"
        )]
        params: JsValue,
    ) -> Result<MolecularForceField, JsValue> {
        let params = UffForceFieldParams::from_js_configuration(&params)?;
        self.uff_force_field_with_params(&params)
    }
    #[wasm_bindgen(js_name = uffForceFieldWithParams)]
    pub fn uff_force_field_with_params(
        &self,
        params: &UffForceFieldParams,
    ) -> Result<MolecularForceField, JsValue> {
        self.inner
            .uff_force_field_with_params(&params.inner)
            .map(|inner| MolecularForceField { inner })
            .map_err(|e| uff_error(&e).unwrap_or_else(|e| e))
    }
}
