//! Aromaticity language projection. All chemistry remains in the public facade.
use crate::Molecule;
use crate::alignment_values::{operation_error, set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use std::sync::Arc;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum AromaticityModel {
    Rdkit,
    Simple,
    Mdl,
    Mmff94,
    Custom,
}

impl AromaticityModel {
    fn core(self) -> ck::AromaticityModel {
        match self {
            Self::Rdkit => ck::AromaticityModel::Rdkit,
            Self::Simple => ck::AromaticityModel::Simple,
            Self::Mdl => ck::AromaticityModel::Mdl,
            Self::Mmff94 => ck::AromaticityModel::Mmff94,
            Self::Custom => ck::AromaticityModel::Custom,
        }
    }
    fn from_core(value: ck::AromaticityModel) -> Self {
        match value {
            ck::AromaticityModel::Rdkit => Self::Rdkit,
            ck::AromaticityModel::Simple => Self::Simple,
            ck::AromaticityModel::Mdl => Self::Mdl,
            ck::AromaticityModel::Mmff94 => Self::Mmff94,
            ck::AromaticityModel::Custom => Self::Custom,
        }
    }
}

#[wasm_bindgen]
pub struct AromaticityParams {
    inner: ck::AromaticityParams,
}

#[wasm_bindgen]
impl AromaticityParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "AromaticityModel | null")] model: Option<
            AromaticityModel,
        >,
    ) -> Self {
        Self {
            inner: model.map_or_else(ck::AromaticityParams::default, |model| {
                ck::AromaticityParams {
                    model: model.core(),
                }
            }),
        }
    }
    #[wasm_bindgen(getter)]
    pub fn model(&self) -> AromaticityModel {
        AromaticityModel::from_core(self.inner.model)
    }
}

#[wasm_bindgen]
pub struct AromaticityError {
    inner: ck::AromaticityError,
}

#[wasm_bindgen]
impl AromaticityError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "aromaticity".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        aromaticity_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter)]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}

fn aromaticity_kind(source: &ck::AromaticityError) -> &'static str {
    use ck::AromaticityError as E;
    match source {
        E::InvalidTopology(..) => "InvalidTopology",
        E::InvalidQueryState(..) => "InvalidQueryState",
        E::RingInfoNotInitialized => "RingInfoNotInitialized",
        E::RingInfoDimensionMismatch { .. } => "RingInfoDimensionMismatch",
        E::RingTableLengthMismatch { .. } => "RingTableLengthMismatch",
        E::RingRowLengthMismatch { .. } => "RingRowLengthMismatch",
        E::RingAtomOutOfRange { .. } => "RingAtomOutOfRange",
        E::RingBondOutOfRange { .. } => "RingBondOutOfRange",
        E::AtomOutOfRange { .. } => "AtomOutOfRange",
        E::BondOutOfRange { .. } => "BondOutOfRange",
        E::ValenceAssignmentLength { .. } => "ValenceAssignmentLength",
        E::InvalidValenceRow { .. } => "InvalidValenceRow",
        E::ExpectedRingBondNotFound { .. } => "ExpectedRingBondNotFound",
        E::UnsupportedModel { .. } => "UnsupportedModel",
        E::UnsupportedState { .. } => "UnsupportedState",
        E::IntegerOverflow { .. } => "IntegerOverflow",
        E::UnexpectedTopologyShape { .. } => "UnexpectedTopologyShape",
        E::AtomIdentityChanged { .. } => "AtomIdentityChanged",
        E::BondIdentityChanged { .. } => "BondIdentityChanged",
        E::AromaticRingCountOutOfRange { .. } => "AromaticRingCountOutOfRange",
        E::RingFinding(..) => "RingFinding",
        E::Valence(..) => "Valence",
        E::Kekulize(..) => "Kekulize",
    }
}

pub(crate) fn aromaticity_error(source: &ck::AromaticityError) -> Result<JsValue, JsValue> {
    use ck::AromaticityError as E;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("AromaticityError");
    let error: JsValue = error.into();
    set(&error, "domain", "aromaticity".into())?;
    set(&error, "kind", aromaticity_kind(source).into())?;
    set(
        &error,
        "detail",
        AromaticityError {
            inner: source.clone(),
        }
        .into(),
    )?;
    let number = |key: &str, value: usize| set(&error, key, JsValue::from_f64(value as f64));
    match source {
        E::RingInfoDimensionMismatch {
            ring_atom_count,
            ring_bond_count,
            topology_atom_count,
            topology_bond_count,
        } => {
            number("ringAtomCount", *ring_atom_count)?;
            number("ringBondCount", *ring_bond_count)?;
            number("topologyAtomCount", *topology_atom_count)?;
            number("topologyBondCount", *topology_bond_count)?;
        }
        E::RingTableLengthMismatch {
            atom_rows,
            bond_rows,
        } => {
            number("atomRows", *atom_rows)?;
            number("bondRows", *bond_rows)?;
        }
        E::RingRowLengthMismatch {
            ring_index,
            atom_count,
            bond_count,
        } => {
            number("ringIndex", *ring_index)?;
            number("atomCount", *atom_count)?;
            number("bondCount", *bond_count)?;
        }
        E::RingAtomOutOfRange {
            ring_index,
            atom,
            atom_count,
        } => {
            number("ringIndex", *ring_index)?;
            number("atom", atom.index())?;
            number("atomCount", *atom_count)?;
        }
        E::RingBondOutOfRange {
            ring_index,
            bond,
            bond_count,
        } => {
            number("ringIndex", *ring_index)?;
            number("bond", bond.index())?;
            number("bondCount", *bond_count)?;
        }
        E::AtomOutOfRange { atom, atom_count } => {
            number("atom", atom.index())?;
            number("atomCount", *atom_count)?;
        }
        E::BondOutOfRange { bond, bond_count } => {
            number("bond", bond.index())?;
            number("bondCount", *bond_count)?;
        }
        E::ValenceAssignmentLength {
            field,
            actual,
            expected,
        } => {
            set(&error, "field", (*field).into())?;
            number("actual", *actual)?;
            number("expected", *expected)?;
        }
        E::InvalidValenceRow { atom, field, value } => {
            number("atom", atom.index())?;
            set(&error, "field", (*field).into())?;
            set(&error, "value", (*value).into())?;
        }
        E::ExpectedRingBondNotFound { begin, end } => {
            number("begin", begin.index())?;
            number("end", end.index())?;
        }
        E::UnsupportedModel { model, detail } => {
            set(
                &error,
                "model",
                JsValue::from_f64(AromaticityModel::from_core(*model) as u32 as f64),
            )?;
            set(&error, "reason", (*detail).into())?;
        }
        E::UnsupportedState { field, detail } => {
            set(&error, "field", (*field).into())?;
            set(&error, "reason", (*detail).into())?;
        }
        E::IntegerOverflow { field } => set(&error, "field", (*field).into())?,
        E::UnexpectedTopologyShape {
            input_atoms,
            input_bonds,
            output_atoms,
            output_bonds,
        } => {
            number("inputAtoms", *input_atoms)?;
            number("inputBonds", *input_bonds)?;
            number("outputAtoms", *output_atoms)?;
            number("outputBonds", *output_bonds)?;
        }
        E::AtomIdentityChanged {
            row,
            expected,
            actual,
        } => {
            number("row", *row)?;
            number("expected", expected.index())?;
            number("actual", actual.index())?;
        }
        E::BondIdentityChanged {
            row,
            expected,
            actual,
        } => {
            number("row", *row)?;
            number("expected", expected.index())?;
            number("actual", actual.index())?;
        }
        E::AromaticRingCountOutOfRange {
            model,
            actual,
            maximum,
        } => {
            set(
                &error,
                "model",
                JsValue::from_f64(AromaticityModel::from_core(*model) as u32 as f64),
            )?;
            number("actual", *actual)?;
            number("maximum", *maximum)?;
        }
        E::InvalidTopology(..)
        | E::InvalidQueryState(..)
        | E::RingInfoNotInitialized
        | E::RingFinding(..)
        | E::Valence(..)
        | E::Kekulize(..) => {}
    }
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = withAssignedAromaticity)]
    pub fn with_assigned_aromaticity(&self) -> Result<Molecule, JsValue> {
        self.inner
            .with_assigned_aromaticity()
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = withAssignedAromaticityWithParams)]
    pub fn with_assigned_aromaticity_with_params(
        &self,
        params: &AromaticityParams,
    ) -> Result<Molecule, JsValue> {
        self.inner
            .with_assigned_aromaticity_with_params(&params.inner)
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = assignAromaticity)]
    pub fn assign_aromaticity_(&self) -> Result<(), JsValue> {
        self.inner
            .assign_aromaticity_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = assignAromaticityWithParams)]
    pub fn assign_aromaticity_with_params_(
        &self,
        params: &AromaticityParams,
    ) -> Result<(), JsValue> {
        self.inner
            .assign_aromaticity_with_params_(&params.inner)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
}
