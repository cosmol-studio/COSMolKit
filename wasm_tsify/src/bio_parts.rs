//! Complete detached structure parts, with native validation on publication.
use crate::{
    bio_assemblies::*,
    bio_entities::*,
    bio_geometry::*,
    bio_hierarchy::{BioCoordinateFormat, BioCoordinateFormatInput, BioStructure},
    bio_metadata::*,
    bio_read_errors::BioStructureError,
    bio_relationships::*,
    bio_rows::*,
    bio_source_state::*,
    fields::Field,
    storage::NativeField,
};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
crate::fields::value_array!(ck::BioModelRow, BioModelRows, "ReadonlyArray<BioModelRow>");
crate::fields::value_array!(ck::BioChainRow, BioChainRows, "ReadonlyArray<BioChainRow>");
crate::fields::value_array!(
    ck::BioResidueRow,
    BioResidueRows,
    "ReadonlyArray<BioResidueRow>"
);
crate::fields::value_array!(ck::BioAtomRow, BioAtomRows, "ReadonlyArray<BioAtomRow>");
#[wasm_bindgen]
pub struct BioStructureParts(ck::BioStructureParts);
crate::fields::value_field!(
    ck::BioStructureParts,
    BioStructureParts,
    BioStructurePartsInput,
    "BioStructureParts"
);
#[wasm_bindgen]
impl BioStructureParts {
    #[wasm_bindgen(constructor)]
    pub fn new(
        input_format: BioCoordinateFormatInput,
        models: BioModelRows,
        chains: BioChainRows,
        residues: BioResidueRows,
        atoms: BioAtomRows,
        entities: BioEntityRows,
        connections: BioConnections,
        cispeps: BioCisPeps,
        mod_residues: BioModResidues,
        helices: BioHelices,
        sheets: BioSheets,
        ncs_operators: BioNcsOperators,
        assemblies: BioAssemblies,
        metadata: BioMetadataInput,
        source_state: BioStructureSourceStateInput,
        coordinates: BioCoordinateBlockInput,
        crystal: Option<BioCrystalInfoInput>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::BioStructureParts {
            input_format: <ck::BioCoordinateFormat as Field>::input(input_format)?,
            models: <Vec<ck::BioModelRow> as Field>::input(models)?,
            chains: <Vec<ck::BioChainRow> as Field>::input(chains)?,
            residues: <Vec<ck::BioResidueRow> as Field>::input(residues)?,
            atoms: <Vec<ck::BioAtomRow> as Field>::input(atoms)?,
            entities: <Vec<ck::BioEntityRow> as Field>::input(entities)?,
            connections: <Vec<ck::BioConnection> as Field>::input(connections)?,
            cispeps: <Vec<ck::BioCisPep> as Field>::input(cispeps)?,
            mod_residues: <Vec<ck::BioModRes> as Field>::input(mod_residues)?,
            helices: <Vec<ck::BioHelix> as Field>::input(helices)?,
            sheets: <Vec<ck::BioSheet> as Field>::input(sheets)?,
            ncs_operators: <Vec<ck::BioNcsOperator> as Field>::input(ncs_operators)?,
            assemblies: <Vec<ck::BioAssembly> as Field>::input(assemblies)?,
            metadata: <ck::BioMetadata as Field>::input(metadata)?,
            source_state: <ck::BioStructureSourceState as Field>::input(source_state)?,
            coordinates: <ck::BioCoordinateBlock as Field>::input(coordinates)?,
            crystal: crystal
                .map(<ck::BioCrystalInfo as Field>::input)
                .transpose()?,
        }))
    }
    #[wasm_bindgen(getter,js_name=inputFormat)]
    pub fn input_format(&self) -> BioCoordinateFormat {
        self.0.input_format.into()
    }
    #[wasm_bindgen(getter,js_name=models)]
    pub fn models(&self) -> Result<BioModelRows, JsValue> {
        <Vec<ck::BioModelRow> as Field>::get(NativeField::owned(self.0.models.clone()))
    }
    #[wasm_bindgen(getter,js_name=chains)]
    pub fn chains(&self) -> Result<BioChainRows, JsValue> {
        <Vec<ck::BioChainRow> as Field>::get(NativeField::owned(self.0.chains.clone()))
    }
    #[wasm_bindgen(getter,js_name=residues)]
    pub fn residues(&self) -> Result<BioResidueRows, JsValue> {
        <Vec<ck::BioResidueRow> as Field>::get(NativeField::owned(self.0.residues.clone()))
    }
    #[wasm_bindgen(getter,js_name=atoms)]
    pub fn atoms(&self) -> Result<BioAtomRows, JsValue> {
        <Vec<ck::BioAtomRow> as Field>::get(NativeField::owned(self.0.atoms.clone()))
    }
    #[wasm_bindgen(getter,js_name=entities)]
    pub fn entities(&self) -> Result<BioEntityRows, JsValue> {
        <Vec<ck::BioEntityRow> as Field>::get(NativeField::owned(self.0.entities.clone()))
    }
    #[wasm_bindgen(getter,js_name=connections)]
    pub fn connections(&self) -> Result<BioConnections, JsValue> {
        <Vec<ck::BioConnection> as Field>::get(NativeField::owned(self.0.connections.clone()))
    }
    #[wasm_bindgen(getter,js_name=cispeps)]
    pub fn cispeps(&self) -> Result<BioCisPeps, JsValue> {
        <Vec<ck::BioCisPep> as Field>::get(NativeField::owned(self.0.cispeps.clone()))
    }
    #[wasm_bindgen(getter,js_name=modResidues)]
    pub fn mod_residues(&self) -> Result<BioModResidues, JsValue> {
        <Vec<ck::BioModRes> as Field>::get(NativeField::owned(self.0.mod_residues.clone()))
    }
    #[wasm_bindgen(getter,js_name=helices)]
    pub fn helices(&self) -> Result<BioHelices, JsValue> {
        <Vec<ck::BioHelix> as Field>::get(NativeField::owned(self.0.helices.clone()))
    }
    #[wasm_bindgen(getter,js_name=sheets)]
    pub fn sheets(&self) -> Result<BioSheets, JsValue> {
        <Vec<ck::BioSheet> as Field>::get(NativeField::owned(self.0.sheets.clone()))
    }
    #[wasm_bindgen(getter,js_name=ncsOperators)]
    pub fn ncs_operators(&self) -> Result<BioNcsOperators, JsValue> {
        <Vec<ck::BioNcsOperator> as Field>::get(NativeField::owned(self.0.ncs_operators.clone()))
    }
    #[wasm_bindgen(getter,js_name=assemblies)]
    pub fn assemblies(&self) -> Result<BioAssemblies, JsValue> {
        <Vec<ck::BioAssembly> as Field>::get(NativeField::owned(self.0.assemblies.clone()))
    }
    #[wasm_bindgen(getter,js_name=metadata)]
    pub fn metadata(&self) -> Result<BioMetadata, JsValue> {
        <ck::BioMetadata as Field>::get(NativeField::owned(self.0.metadata.clone()))
    }
    #[wasm_bindgen(getter,js_name=sourceState)]
    pub fn source_state(&self) -> Result<BioStructureSourceState, JsValue> {
        <ck::BioStructureSourceState as Field>::get(NativeField::owned(self.0.source_state.clone()))
    }
    #[wasm_bindgen(getter,js_name=coordinates)]
    pub fn coordinates(&self) -> Result<BioCoordinateBlock, JsValue> {
        <ck::BioCoordinateBlock as Field>::get(NativeField::owned(self.0.coordinates.clone()))
    }
    #[wasm_bindgen(getter)]
    pub fn crystal(&self) -> Result<Option<BioCrystalInfo>, JsValue> {
        <Option<ck::BioCrystalInfo> as Field>::get(NativeField::owned(self.0.crystal.clone()))
    }
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=fromParts)]
    pub fn from_parts(parts: BioStructurePartsInput) -> Result<Self, JsValue> {
        ck::BioStructure::from_parts(<ck::BioStructureParts as Field>::input(parts)?)
            .map(Self)
            .map_err(BioStructureError::from_rust)
    }
    #[wasm_bindgen(js_name=validateParts)]
    pub fn validate_parts(parts: BioStructurePartsInput) -> Result<(), JsValue> {
        BioStructureParts::with_ref(&parts, |p| {
            ck::BioStructure::validate_parts(p).map_err(BioStructureError::from_rust)
        })
    }
    pub fn validate(&self) -> Result<(), JsValue> {
        self.0.validate().map_err(BioStructureError::from_rust)
    }
    #[wasm_bindgen(js_name=intoParts)]
    pub fn into_parts(&self) -> BioStructureParts {
        BioStructureParts(self.0.clone().into_parts())
    }
}
