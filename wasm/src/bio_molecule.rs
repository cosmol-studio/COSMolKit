//! BIO text-to-molecule composition uses the canonical facade entrypoints.
use crate::{Molecule, bio_hierarchy::BioStructure, fields::Field, input};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
ck::__ck_configuration_BioMoleculeParams!(crate::project_parameters, ck::BioMoleculeParams, {}, {
    BioMoleculeParams
});
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=toMoleculeWithParams)]
    pub fn to_molecule_with_params(
        &self,
        params: BioMoleculeParamsInput,
    ) -> Result<Molecule, JsValue> {
        self.to_molecule(Some(params))
    }
    #[wasm_bindgen(js_name=toMolecule)]
    pub fn to_molecule(&self, params: Option<BioMoleculeParamsInput>) -> Result<Molecule, JsValue> {
        BioMoleculeParams::with_input(params.as_ref(), |p| {
            self.0
                .to_molecule_with_params(p)
                .map(Molecule)
                .map_err(crate::bio_conversion_errors::BioMoleculeError::from_rust)
        })
    }
}
#[wasm_bindgen(js_name=molFromPdb)]
pub fn mol_from_pdb(
    text: Ts<input::TextInput>,
    params: Option<BioMoleculeParamsInput>,
) -> Result<Molecule, JsValue> {
    Molecule::from_pdb(text, params)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fromPdb)]
    pub fn from_pdb(
        text: Ts<input::TextInput>,
        params: Option<BioMoleculeParamsInput>,
    ) -> Result<Self, JsValue> {
        let text = input::read::<String>(text)?;
        BioMoleculeParams::with_input(params.as_ref(), |p| {
            ck::Molecule::from_pdb_with_params(&text, p)
                .map(Self)
                .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        })
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = fromPdbWithParams)]
    pub fn from_pdb_with_params(
        text: Ts<input::TextInput>,
        params: BioMoleculeParamsInput,
    ) -> Result<Self, JsValue> {
        Self::from_pdb(text, Some(params))
    }
    #[wasm_bindgen(js_name = fromMmcifWithParams)]
    pub fn from_mmcif_with_params(
        text: Ts<input::TextInput>,
        params: BioMoleculeParamsInput,
    ) -> Result<Self, JsValue> {
        Self::from_mmcif(text, Some(params))
    }
}
#[wasm_bindgen(js_name=molFromMmcif)]
pub fn mol_from_mmcif(
    text: Ts<input::TextInput>,
    params: Option<BioMoleculeParamsInput>,
) -> Result<Molecule, JsValue> {
    Molecule::from_mmcif(text, params)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fromMmcif)]
    pub fn from_mmcif(
        text: Ts<input::TextInput>,
        params: Option<BioMoleculeParamsInput>,
    ) -> Result<Self, JsValue> {
        let text = input::read::<String>(text)?;
        BioMoleculeParams::with_input(params.as_ref(), |p| {
            ck::Molecule::from_mmcif_with_params(&text, p)
                .map(Self)
                .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        })
    }
}
