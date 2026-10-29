//! Native code computation returns its explicitly updated COW molecule.
use crate::{Molecule, MoleculeInput, input, operation_errors::OperationError};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct AtomPairAtomCodeResult(ck::AtomPairAtomCodeResult);
#[wasm_bindgen]
impl AtomPairAtomCodeResult {
    #[wasm_bindgen(constructor)]
    pub fn new(code: Ts<input::U32Input>, molecule: MoleculeInput) -> Result<Self, JsValue> {
        let code = input::read(code)?;
        Molecule::with_ref(&molecule, |molecule| {
            Ok(Self(ck::AtomPairAtomCodeResult {
                code,
                molecule: molecule.clone(),
            }))
        })
    }
    #[wasm_bindgen(getter)]
    pub fn code(&self) -> u32 {
        self.0.code
    }
    #[wasm_bindgen(getter)]
    pub fn molecule(&self) -> Molecule {
        Molecule(self.0.molecule.clone())
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=withAtomPairAtomCode)]
    pub fn with_atom_pair_atom_code(
        &self,
        atom_id: Ts<input::IndexInput>,
        branch_subtract: Option<Ts<input::U32Input>>,
        include_chirality: Option<Ts<input::BoolInput>>,
        use_legacy_stereo_perception: Option<Ts<input::BoolInput>>,
    ) -> Result<AtomPairAtomCodeResult, JsValue> {
        self.0
            .with_atom_pair_atom_code(
                ck::AtomId::new(input::read(atom_id)?),
                branch_subtract.map(input::read).transpose()?.unwrap_or(0),
                include_chirality
                    .map(input::read)
                    .transpose()?
                    .unwrap_or(false),
                use_legacy_stereo_perception
                    .map(input::read)
                    .transpose()?
                    .unwrap_or(true),
            )
            .map(AtomPairAtomCodeResult)
            .map_err(OperationError::from_rust)
    }
}
