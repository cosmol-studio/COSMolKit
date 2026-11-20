//! Registered molecule fragment operations retain native ordering and ownership.
use crate::{
    Molecule, fields::Field, fragment_errors::BondCutError, input,
    operation_errors::OperationError, storage::NativeField,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
crate::fields::value_array!(ck::Molecule, Molecules, "ReadonlyArray<Molecule>");

#[wasm_bindgen]
pub struct BondCut(ck::BondCut);
crate::fields::value_field!(ck::BondCut, BondCut, BondCutInput, "BondCut");
crate::fields::value_array!(ck::BondCut, BondCuts, "ReadonlyArray<BondCut>");

#[wasm_bindgen]
impl BondCut {
    /// Supply both attachment-side isotope labels, or omit both to label each dummy
    /// with the index of the atom it replaces. Effective labels must fit in 16 bits.
    #[wasm_bindgen(js_name=atBond)]
    pub fn at_bond(
        bond: Ts<input::IndexInput>,
        begin_label: Option<Ts<input::IndexInput>>,
        end_label: Option<Ts<input::IndexInput>>,
    ) -> Result<Self, JsValue> {
        ck::BondCut::at_bond(
            ck::BondId::new(input::read::<usize>(bond)?),
            begin_label.map(input::read::<usize>).transpose()?,
            end_label.map(input::read::<usize>).transpose()?,
        )
        .map(Self)
        .map_err(BondCutError::from_rust)
    }
}

#[wasm_bindgen]
impl Molecule {
    /// Return the cut molecule; use fragments() to extract connected components.
    #[wasm_bindgen(js_name=withBondsCut)]
    pub fn with_bonds_cut(&self, cuts: BondCuts) -> Result<Self, JsValue> {
        self.0
            .with_bonds_cut(&Vec::<ck::BondCut>::input(cuts)?)
            .map(Self)
            .map_err(OperationError::from_rust)
    }

    /// Apply the same selected cuts to this molecule in place.
    #[wasm_bindgen(js_name=cutBonds_)]
    pub fn cut_bonds_(&mut self, cuts: BondCuts) -> Result<(), JsValue> {
        self.0
            .cut_bonds_(&Vec::<ck::BondCut>::input(cuts)?)
            .map_err(OperationError::from_rust)
    }

    #[wasm_bindgen(js_name=fragments)]
    pub fn fragments(&self) -> Result<Molecules, JsValue> {
        Vec::<ck::Molecule>::get(NativeField::owned(
            self.0.fragments().map_err(OperationError::from_rust)?,
        ))
    }
    #[wasm_bindgen(js_name=largestFragment)]
    pub fn largest_fragment(&self) -> Result<Self, JsValue> {
        self.0
            .largest_fragment()
            .map(Self)
            .map_err(OperationError::from_rust)
    }
}
