//! Text IO calls use the same native Params/options boundary as molecular APIs.
use crate::{Molecule, input, sdf_errors::SdfError};
use crate::{io_configuration::*, molecular_io_errors::MolecularIoError};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;

ck::__ck_configuration_SdfReadParams!(crate::project_parameters, ck::SdfReadParams, {}, {
    SdfReadParams
});

#[wasm_bindgen(js_name = molFromSdf)]
pub fn mol_from_sdf(
    text: tsify::Ts<input::TextInput>,
    params: Option<SdfReadParamsInput>,
) -> Result<Molecule, JsValue> {
    Molecule::from_sdf(text, params)
}

#[wasm_bindgen(js_name = molFromMol)]
pub fn mol_from_mol(
    text: tsify::Ts<input::TextInput>,
    params: Option<SdfReadParamsInput>,
) -> Result<Molecule, JsValue> {
    Molecule::from_mol(text, params)
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fromXyzBlock)]
    pub fn from_xyz_block(text: tsify::Ts<input::TextInput>) -> Result<Molecule, JsValue> {
        let text = input::read::<String>(text)?;
        ck::Molecule::from_xyz_block(&text)
            .map(Molecule)
            .map_err(MolecularIoError::from_rust)
    }
    #[wasm_bindgen(js_name=fromMol2)]
    pub fn from_mol2(
        text: tsify::Ts<input::TextInput>,
        params: Option<Mol2ReadParamsInput>,
    ) -> Result<Molecule, JsValue> {
        let text = input::read::<String>(text)?;
        Mol2ReadParams::with_input(params.as_ref(), |params| {
            ck::Molecule::from_mol2_with_params(&text, params)
                .map(Molecule)
                .map_err(MolecularIoError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toMol)]
    pub fn to_mol(&self, params: Option<MolBlockWriteParamsInput>) -> Result<String, JsValue> {
        MolBlockWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_mol_with_params(p)
                .map_err(MolecularIoError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toSdf)]
    pub fn to_sdf(&self, params: Option<MolBlockWriteParamsInput>) -> Result<String, JsValue> {
        MolBlockWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_sdf_with_params(p)
                .map_err(MolecularIoError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toSdf2d)]
    pub fn to_sdf_2d(&self, params: Option<MolBlockWriteParamsInput>) -> Result<String, JsValue> {
        MolBlockWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_sdf_2d_with_params(p)
                .map_err(MolecularIoError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toSdf3d)]
    pub fn to_sdf_3d(&self, params: Option<MolBlockWriteParamsInput>) -> Result<String, JsValue> {
        MolBlockWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_sdf_3d_with_params(p)
                .map_err(MolecularIoError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toPdb)]
    pub fn to_pdb(&self, params: Option<PdbWriteParamsInput>) -> Result<String, JsValue> {
        PdbWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_pdb_with_params(p)
                .map_err(MolecularIoError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toXyz)]
    pub fn to_xyz(&self, params: Option<XyzWriteParamsInput>) -> Result<String, JsValue> {
        XyzWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_xyz_with_params(p)
                .map_err(MolecularIoError::from_rust)
        })
    }
    #[wasm_bindgen(js_name = fromSdf)]
    pub fn from_sdf(
        text: tsify::Ts<input::TextInput>,
        params: Option<SdfReadParamsInput>,
    ) -> Result<Molecule, JsValue> {
        let text = input::read::<String>(text)?;
        SdfReadParams::with_input(params.as_ref(), |params| {
            ck::Molecule::from_sdf_with_params(&text, params)
                .map(Molecule)
                .map_err(SdfError::from_rust)
        })
    }

    #[wasm_bindgen(js_name = fromMol)]
    pub fn from_mol(
        text: tsify::Ts<input::TextInput>,
        params: Option<SdfReadParamsInput>,
    ) -> Result<Molecule, JsValue> {
        let text = input::read::<String>(text)?;
        SdfReadParams::with_input(params.as_ref(), |params| {
            ck::Molecule::from_mol_with_params(&text, params)
                .map(Molecule)
                .map_err(SdfError::from_rust)
        })
    }
}

#[wasm_bindgen(js_name=molFromXyzBlock)]
pub fn mol_from_xyz_block(text: tsify::Ts<input::TextInput>) -> Result<Molecule, JsValue> {
    Molecule::from_xyz_block(text)
}
#[wasm_bindgen(js_name=molFromMol2)]
pub fn mol_from_mol2(
    text: tsify::Ts<input::TextInput>,
    params: Option<Mol2ReadParamsInput>,
) -> Result<Molecule, JsValue> {
    Molecule::from_mol2(text, params)
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = fromSdfWithParams)]
    pub fn from_sdf_with_params(
        text: tsify::Ts<input::TextInput>,
        params: SdfReadParamsInput,
    ) -> Result<Molecule, JsValue> {
        Self::from_sdf(text, Some(params))
    }
    #[wasm_bindgen(js_name = toPdbWithParams)]
    pub fn to_pdb_with_params(&self, params: PdbWriteParamsInput) -> Result<String, JsValue> {
        self.to_pdb(Some(params))
    }
    #[wasm_bindgen(js_name = fromMol2WithParams)]
    pub fn from_mol2_with_params(
        text: tsify::Ts<input::TextInput>,
        params: Mol2ReadParamsInput,
    ) -> Result<Molecule, JsValue> {
        Self::from_mol2(text, Some(params))
    }
    #[wasm_bindgen(js_name = fromMolWithParams)]
    pub fn from_mol_with_params(
        text: tsify::Ts<input::TextInput>,
        params: SdfReadParamsInput,
    ) -> Result<Molecule, JsValue> {
        Self::from_mol(text, Some(params))
    }
    #[wasm_bindgen(js_name = toXyzWithParams)]
    pub fn to_xyz_with_params(&self, params: XyzWriteParamsInput) -> Result<String, JsValue> {
        self.to_xyz(Some(params))
    }
    #[wasm_bindgen(js_name = toMolWithParams)]
    pub fn to_mol_with_params(&self, params: MolBlockWriteParamsInput) -> Result<String, JsValue> {
        self.to_mol(Some(params))
    }
    #[wasm_bindgen(js_name = toSdfWithParams)]
    pub fn to_sdf_with_params(&self, params: MolBlockWriteParamsInput) -> Result<String, JsValue> {
        self.to_sdf(Some(params))
    }
    #[wasm_bindgen(js_name = toSdf2dWithParams)]
    pub fn to_sdf_2d_with_params(
        &self,
        params: MolBlockWriteParamsInput,
    ) -> Result<String, JsValue> {
        self.to_sdf_2d(Some(params))
    }
    #[wasm_bindgen(js_name = toSdf3dWithParams)]
    pub fn to_sdf_3d_with_params(
        &self,
        params: MolBlockWriteParamsInput,
    ) -> Result<String, JsValue> {
        self.to_sdf_3d(Some(params))
    }
}
