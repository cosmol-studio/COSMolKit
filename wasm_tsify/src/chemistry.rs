//! Thin molecular value/in-place projections and owned ring query results.
use crate::{Molecule, operation_errors::OperationError};
use cosmolkit as ck;
use serde::Serialize;
use tsify::Tsify;
use wasm_bindgen::prelude::*;
ck::__ck_configuration_AddHsParams!(crate::project_parameters, ck::AddHsParams, {}, {
    AddHsParams
});
ck::__ck_configuration_RemoveHsParams!(crate::project_parameters, ck::RemoveHsParams, {}, {
    RemoveHsParams
});
ck::__ck_configuration_KekulizeParams!(crate::project_parameters, ck::KekulizeParams, {}, {
    KekulizeParams
});
ck::__ck_configuration_StructureTagParams!(
    crate::project_parameters,
    ck::StructureTagParams,
    {},
    { StructureTagParams }
);
#[wasm_bindgen]
impl Molecule {
    #[cfg(feature = "cap-hashing")]
    #[wasm_bindgen(js_name=murckoScaffold)]
    pub fn murcko_scaffold(&self) -> Result<Molecule, JsValue> {
        self.0
            .murcko_scaffold()
            .map(Molecule)
            .map_err(OperationError::from_rust)
    }
    #[cfg(feature = "cap-hashing")]
    #[wasm_bindgen(js_name=netScaffold)]
    pub fn net_scaffold(&self) -> Result<Molecule, JsValue> {
        self.0
            .net_scaffold()
            .map(Molecule)
            .map_err(OperationError::from_rust)
    }
    #[cfg(feature = "cap-hashing")]
    #[wasm_bindgen(js_name=murckoDecompose)]
    pub fn murcko_decompose(&self) -> Result<Molecule, JsValue> {
        self.0
            .murcko_decompose()
            .map(Molecule)
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=withAssignedRings)]
    pub fn with_assigned_rings(&self) -> Result<Molecule, JsValue> {
        self.0
            .with_assigned_rings()
            .map(Molecule)
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=assignRings_)]
    pub fn assign_rings_(&mut self) -> Result<(), JsValue> {
        self.0.assign_rings_().map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=withAssignedRingFamilies)]
    pub fn with_assigned_ring_families(
        &self,
        params: Option<crate::RingSearchParamsInput>,
    ) -> Result<Molecule, JsValue> {
        crate::RingSearchParams::with_input(params.as_ref(), |params| {
            self.0
                .with_assigned_ring_families_with_params(params)
                .map(Molecule)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=assignRingFamilies_)]
    pub fn assign_ring_families_(
        &mut self,
        params: Option<crate::RingSearchParamsInput>,
    ) -> Result<(), JsValue> {
        crate::RingSearchParams::with_input(params.as_ref(), |params| {
            self.0
                .assign_ring_families_with_params_(params)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withHydrogens)]
    pub fn with_hydrogens(&self, params: Option<AddHsParamsInput>) -> Result<Molecule, JsValue> {
        AddHsParams::with_input(params.as_ref(), |params| {
            self.0
                .with_hydrogens_with_params(params)
                .map(Molecule)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=addHydrogens_)]
    pub fn add_hydrogens_(&mut self, params: Option<AddHsParamsInput>) -> Result<(), JsValue> {
        AddHsParams::with_input(params.as_ref(), |params| {
            self.0
                .add_hydrogens_with_params_(params)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withoutHydrogens)]
    pub fn without_hydrogens(
        &self,
        params: Option<RemoveHsParamsInput>,
    ) -> Result<Molecule, JsValue> {
        RemoveHsParams::with_input(params.as_ref(), |params| {
            self.0
                .without_hydrogens_with_params(params)
                .map(Molecule)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=removeHydrogens_)]
    pub fn remove_hydrogens_(
        &mut self,
        params: Option<RemoveHsParamsInput>,
    ) -> Result<(), JsValue> {
        RemoveHsParams::with_input(params.as_ref(), |params| {
            self.0
                .remove_hydrogens_with_params_(params)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withKekulizedBonds)]
    pub fn with_kekulized_bonds(
        &self,
        params: Option<KekulizeParamsInput>,
    ) -> Result<Molecule, JsValue> {
        KekulizeParams::with_input(params.as_ref(), |params| {
            self.0
                .with_kekulized_bonds_with_params(params)
                .map(Molecule)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=kekulizeBonds_)]
    pub fn kekulize_bonds_(&mut self, params: Option<KekulizeParamsInput>) -> Result<(), JsValue> {
        KekulizeParams::with_input(params.as_ref(), |params| {
            self.0
                .kekulize_bonds_with_params_(params)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withChiralTagsFromStructure)]
    pub fn with_chiral_tags_from_structure(
        &self,
        params: Option<StructureTagParamsInput>,
    ) -> Result<Molecule, JsValue> {
        StructureTagParams::with_input(params.as_ref(), |params| {
            self.0
                .with_chiral_tags_from_structure_with_params(params)
                .map(Molecule)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=assignChiralTagsFromStructure_)]
    pub fn assign_chiral_tags_from_structure_(
        &mut self,
        params: Option<StructureTagParamsInput>,
    ) -> Result<(), JsValue> {
        StructureTagParams::with_input(params.as_ref(), |params| {
            self.0
                .assign_chiral_tags_from_structure_with_params_(params)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=symmetrizedSssr)]
    pub fn symmetrized_sssr(
        &self,
        params: Option<crate::RingSearchParamsInput>,
    ) -> Result<RingSearchResult, JsValue> {
        crate::RingSearchParams::with_input(params.as_ref(), |params| {
            self.0
                .symmetrized_sssr_with_params(params)
                .map(RingSearchResult)
                .map_err(OperationError::from_rust)
        })
    }
}
#[derive(Serialize, Tsify)]
#[serde(transparent)]
pub struct RingRows(#[tsify(type = "ReadonlyArray<ReadonlyArray<number>>")] Vec<Vec<usize>>);
#[wasm_bindgen]
pub struct RingSearchResult(ck::RingSearchResult);
#[wasm_bindgen]
impl RingSearchResult {
    #[wasm_bindgen(getter,js_name=numRings)]
    pub fn num_rings(&self) -> usize {
        self.0.num_rings()
    }
    #[wasm_bindgen(getter,js_name=atomRings)]
    pub fn atom_rings(&self) -> Result<tsify::Ts<RingRows>, JsValue> {
        ring_rows(
            self.0
                .atom_rings()
                .iter()
                .map(|row| row.iter().map(|id| id.index()).collect())
                .collect(),
        )
    }
    #[wasm_bindgen(getter,js_name=bondRings)]
    pub fn bond_rings(&self) -> Result<tsify::Ts<RingRows>, JsValue> {
        ring_rows(
            self.0
                .bond_rings()
                .iter()
                .map(|row| row.iter().map(|id| id.index()).collect())
                .collect(),
        )
    }
}
fn ring_rows(rows: Vec<Vec<usize>>) -> Result<tsify::Ts<RingRows>, JsValue> {
    RingRows(rows)
        .into_ts()
        .map_err(|error| JsError::from(error).into())
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = withKekulizedBondsWithParams)]
    pub fn with_kekulized_bonds_with_params(
        &self,
        params: KekulizeParamsInput,
    ) -> Result<Molecule, JsValue> {
        self.with_kekulized_bonds(Some(params))
    }
    #[wasm_bindgen(js_name = kekulizeBondsWithParams_)]
    pub fn kekulize_bonds_with_params_(
        &mut self,
        params: KekulizeParamsInput,
    ) -> Result<(), JsValue> {
        self.kekulize_bonds_(Some(params))
    }
    #[wasm_bindgen(js_name = withChiralTagsFromStructureWithParams)]
    pub fn with_chiral_tags_from_structure_with_params(
        &self,
        params: StructureTagParamsInput,
    ) -> Result<Molecule, JsValue> {
        self.with_chiral_tags_from_structure(Some(params))
    }
    #[wasm_bindgen(js_name = assignChiralTagsFromStructureWithParams_)]
    pub fn assign_chiral_tags_from_structure_with_params_(
        &mut self,
        params: StructureTagParamsInput,
    ) -> Result<(), JsValue> {
        self.assign_chiral_tags_from_structure_(Some(params))
    }
    #[wasm_bindgen(js_name = withHydrogensWithParams)]
    pub fn with_hydrogens_with_params(
        &self,
        params: AddHsParamsInput,
    ) -> Result<Molecule, JsValue> {
        self.with_hydrogens(Some(params))
    }
    #[wasm_bindgen(js_name = addHydrogensWithParams_)]
    pub fn add_hydrogens_with_params_(&mut self, params: AddHsParamsInput) -> Result<(), JsValue> {
        self.add_hydrogens_(Some(params))
    }
    #[wasm_bindgen(js_name = withoutHydrogensWithParams)]
    pub fn without_hydrogens_with_params(
        &self,
        params: RemoveHsParamsInput,
    ) -> Result<Molecule, JsValue> {
        self.without_hydrogens(Some(params))
    }
    #[wasm_bindgen(js_name = removeHydrogensWithParams_)]
    pub fn remove_hydrogens_with_params_(
        &mut self,
        params: RemoveHsParamsInput,
    ) -> Result<(), JsValue> {
        self.remove_hydrogens_(Some(params))
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=symmetrizedSssrWithParams)]
    pub fn symmetrized_sssr_with_params(
        &self,
        params: crate::RingSearchParamsInput,
    ) -> Result<RingSearchResult, JsValue> {
        self.symmetrized_sssr(Some(params))
    }
    #[wasm_bindgen(js_name=withAssignedRingFamiliesWithParams)]
    pub fn with_assigned_ring_families_with_params(
        &self,
        params: crate::RingSearchParamsInput,
    ) -> Result<Molecule, JsValue> {
        self.with_assigned_ring_families(Some(params))
    }
    #[wasm_bindgen(js_name=assignRingFamiliesWithParams_)]
    pub fn assign_ring_families_with_params_(
        &mut self,
        params: crate::RingSearchParamsInput,
    ) -> Result<(), JsValue> {
        self.assign_ring_families_(Some(params))
    }
}
