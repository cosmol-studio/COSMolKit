//! Assembly values are native detached records, not a second assembly algorithm.
use crate::{
    bio_geometry::BioAssemblyOperators,
    fields::{Field, TextArray},
    input,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
crate::enums::enum_value!(BioAssemblySpecialKind,"BioAssemblySpecialKind"{NotApplicable="not_applicable",CompleteIcosahedral="complete_icosahedral",RepresentativeHelical="representative_helical",CompletePoint="complete_point"});
#[wasm_bindgen]
pub struct BioAssemblyGenerator(pub(crate) ck::BioAssemblyGenerator);
crate::fields::value_field!(
    ck::BioAssemblyGenerator,
    BioAssemblyGenerator,
    BioAssemblyGeneratorInput,
    "BioAssemblyGenerator"
);
crate::fields::value_array!(
    ck::BioAssemblyGenerator,
    BioAssemblyGenerators,
    "ReadonlyArray<BioAssemblyGenerator>"
);
#[wasm_bindgen]
pub struct BioAssembly(pub(crate) ck::BioAssembly);
crate::fields::value_field!(
    ck::BioAssembly,
    BioAssembly,
    BioAssemblyInput,
    "BioAssembly"
);
crate::fields::value_array!(ck::BioAssembly, BioAssemblies, "ReadonlyArray<BioAssembly>");
#[wasm_bindgen]
impl BioAssemblyGenerator {
    #[wasm_bindgen(constructor)]
    pub fn new(
        chains: Option<Ts<TextArray>>,
        subchains: Option<Ts<TextArray>>,
        operators: Option<BioAssemblyOperators>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioAssemblyGenerator::default();
        Ok(Self(ck::BioAssemblyGenerator::new(
            chains
                .map(<Vec<String> as Field>::input)
                .transpose()?
                .unwrap_or(d.chains),
            subchains
                .map(<Vec<String> as Field>::input)
                .transpose()?
                .unwrap_or(d.subchains),
            operators
                .map(<Vec<ck::BioAssemblyOperator> as Field>::input)
                .transpose()?
                .unwrap_or(d.operators),
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn chains(&self) -> Result<Ts<TextArray>, JsValue> {
        <Vec<String> as Field>::get(NativeField::owned(self.0.chains.clone()))
    }
    #[wasm_bindgen(getter)]
    pub fn subchains(&self) -> Result<Ts<TextArray>, JsValue> {
        <Vec<String> as Field>::get(NativeField::owned(self.0.subchains.clone()))
    }
    #[wasm_bindgen(getter)]
    pub fn operators(&self) -> Result<BioAssemblyOperators, JsValue> {
        <Vec<ck::BioAssemblyOperator> as Field>::get(NativeField::owned(self.0.operators.clone()))
    }
}
#[wasm_bindgen]
impl BioAssembly {
    #[wasm_bindgen(constructor)]
    pub fn new(
        name: Ts<input::TextInput>,
        author_determined: Ts<input::BoolInput>,
        software_determined: Ts<input::BoolInput>,
        special_kind: BioAssemblySpecialKindInput,
        oligomeric_count: Ts<input::I32Input>,
        oligomeric_details: Ts<input::TextInput>,
        software_name: Ts<input::TextInput>,
        buried_surface_area: Ts<input::F64Input>,
        surface_area: Ts<input::F64Input>,
        solvent_free_energy_change: Ts<input::F64Input>,
        generators: BioAssemblyGenerators,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::BioAssembly::new(
            input::read::<String>(name)?,
            input::read::<bool>(author_determined)?,
            input::read::<bool>(software_determined)?,
            <ck::BioAssemblySpecialKind as Field>::input(special_kind)?,
            input::read::<i32>(oligomeric_count)?,
            input::read::<String>(oligomeric_details)?,
            input::read::<String>(software_name)?,
            input::read::<f64>(buried_surface_area)?,
            input::read::<f64>(surface_area)?,
            input::read::<f64>(solvent_free_energy_change)?,
            <Vec<ck::BioAssemblyGenerator> as Field>::input(generators)?,
        )))
    }
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        self.0.name.clone()
    }
    #[wasm_bindgen(getter,js_name=authorDetermined)]
    pub fn author_determined(&self) -> bool {
        self.0.author_determined
    }
    #[wasm_bindgen(getter,js_name=softwareDetermined)]
    pub fn software_determined(&self) -> bool {
        self.0.software_determined
    }
    #[wasm_bindgen(getter,js_name=specialKind)]
    pub fn special_kind(&self) -> BioAssemblySpecialKind {
        self.0.special_kind.into()
    }
    #[wasm_bindgen(getter,js_name=oligomericCount)]
    pub fn oligomeric_count(&self) -> i32 {
        self.0.oligomeric_count
    }
    #[wasm_bindgen(getter,js_name=oligomericDetails)]
    pub fn oligomeric_details(&self) -> String {
        self.0.oligomeric_details.clone()
    }
    #[wasm_bindgen(getter,js_name=softwareName)]
    pub fn software_name(&self) -> String {
        self.0.software_name.clone()
    }
    #[wasm_bindgen(getter,js_name=buriedSurfaceArea)]
    pub fn buried_surface_area(&self) -> f64 {
        self.0.buried_surface_area
    }
    #[wasm_bindgen(getter,js_name=surfaceArea)]
    pub fn surface_area(&self) -> f64 {
        self.0.surface_area
    }
    #[wasm_bindgen(getter,js_name=solventFreeEnergyChange)]
    pub fn solvent_free_energy_change(&self) -> f64 {
        self.0.solvent_free_energy_change
    }
    #[wasm_bindgen(getter)]
    pub fn generators(&self) -> Result<BioAssemblyGenerators, JsValue> {
        <Vec<ck::BioAssemblyGenerator> as Field>::get(NativeField::owned(self.0.generators.clone()))
    }
}
