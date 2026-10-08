//! Complete query-preserving record and concrete molecule read boundaries.
use crate::group_values::SubstanceGroup;
use crate::host_values::{bool_value, type_error, u32_value};
use crate::io_errors::{io_error, sdf_error};
use crate::property_values::{MoleculeProperties, fields};
use crate::Molecule;
use crate::query_values::QueryGraph;
use cosmolkit_wasm::rust as ck;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum CoordinateDimension {
    TwoD,
    ThreeD,
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum SdfCoordinateMode {
    Preserve,
    Require2D,
    Require3D,
}
#[wasm_bindgen]
pub struct SdfReadParams {
    pub(crate) inner: ck::SdfReadParams,
}
#[wasm_bindgen]
impl SdfReadParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] sanitize: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_hydrogens: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] strict_parsing: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        expand_attachment_points: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] process_property_lists: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "SdfCoordinateMode")]
        coordinate_mode: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::SdfReadParams::default();
        if !sanitize.is_undefined() {
            inner.sanitize = bool_value(&sanitize, "sanitize")?;
        }
        if !remove_hydrogens.is_undefined() {
            inner.remove_hydrogens = bool_value(&remove_hydrogens, "removeHydrogens")?;
        }
        if !strict_parsing.is_undefined() {
            inner.strict_parsing = bool_value(&strict_parsing, "strictParsing")?;
        }
        if !expand_attachment_points.is_undefined() {
            inner.expand_attachment_points =
                bool_value(&expand_attachment_points, "expandAttachmentPoints")?;
        }
        if !process_property_lists.is_undefined() {
            inner.process_property_lists =
                bool_value(&process_property_lists, "processPropertyLists")?;
        }
        if !coordinate_mode.is_undefined() {
            inner.coordinate_mode = match u32_value(&coordinate_mode, "coordinateMode")? {
                0 => ck::SdfCoordinateMode::Preserve,
                1 => ck::SdfCoordinateMode::Require2D,
                2 => ck::SdfCoordinateMode::Require3D,
                _ => return Err(js_sys::RangeError::new("invalid SdfCoordinateMode").into()),
            };
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=sanitize)]
    pub fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[wasm_bindgen(getter,js_name=removeHydrogens)]
    pub fn remove_hydrogens(&self) -> bool {
        self.inner.remove_hydrogens
    }
    #[wasm_bindgen(getter,js_name=strictParsing)]
    pub fn strict_parsing(&self) -> bool {
        self.inner.strict_parsing
    }
    #[wasm_bindgen(getter,js_name=expandAttachmentPoints)]
    pub fn expand_attachment_points(&self) -> bool {
        self.inner.expand_attachment_points
    }
    #[wasm_bindgen(getter,js_name=processPropertyLists)]
    pub fn process_property_lists(&self) -> bool {
        self.inner.process_property_lists
    }
    #[wasm_bindgen(getter,js_name=coordinateMode,unchecked_return_type="SdfCoordinateMode")]
    pub fn coordinate_mode(&self) -> u32 {
        match self.inner.coordinate_mode {
            ck::SdfCoordinateMode::Preserve => 0,
            ck::SdfCoordinateMode::Require2D => 1,
            ck::SdfCoordinateMode::Require3D => 2,
        }
    }
}
#[wasm_bindgen(
    inline_js = "export function visitSdfReadParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid SdfReadParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitSdfReadParams)]
    fn visit_params(value: &JsValue, visit: &mut dyn FnMut(&SdfReadParams)) -> Result<(), JsValue>;
}
#[wasm_bindgen]
pub struct SdfGraph {
    pub(crate) inner: cosmolkit_wasm::SdfGraph,
}
#[wasm_bindgen]
impl SdfGraph {
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.inner.kind().to_owned()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Molecule | null")]
    pub fn molecule(&self) -> JsValue {
        self.inner.molecule().map_or(JsValue::NULL, |inner| {
            Molecule {
                inner: Arc::new(inner),
            }
            .into()
        })
    }
    #[wasm_bindgen(getter,js_name=queryGraph,unchecked_return_type="QueryGraph | null")]
    pub fn query_graph(&self) -> JsValue {
        self.inner
            .query_graph()
            .map_or(JsValue::NULL, |inner| QueryGraph { inner }.into())
    }
}
#[wasm_bindgen]
pub struct SdfRecord {
    pub(crate) inner: cosmolkit_wasm::SdfRecord,
}
#[wasm_bindgen]
impl SdfRecord {
    #[wasm_bindgen(js_name=fromSdf)]
    pub fn from_sdf(text: &str) -> Result<Self, JsValue> {
        cosmolkit_wasm::SdfRecord::from_sdf(text)
            .map(|inner| Self { inner })
            .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromSdfWithParams)]
    pub fn from_sdf_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::SdfRecord::from_sdf_with_params(text,&params.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SdfReadParams| {
            result = Some(
                cosmolkit_wasm::SdfRecord::from_sdf_with_params(text, &p.inner)
                    .map(|inner| Self { inner })
                    .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    pub fn graph(&self) -> SdfGraph {
        SdfGraph {
            inner: self.inner.graph(),
        }
    }
    pub fn molecule(&self) -> Result<Molecule, JsValue> {
        self.inner
            .molecule()
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=queryGraph)]
    pub fn query_graph(&self) -> Result<QueryGraph, JsValue> {
        self.inner
            .query_graph()
            .map(|inner| QueryGraph { inner })
            .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=dataFields,unchecked_return_type="[string, string][]")]
    pub fn data_fields(&self) -> Result<JsValue, JsValue> {
        fields(&self.inner.data_fields())
    }
    pub fn properties(&self) -> MoleculeProperties {
        MoleculeProperties {
            inner: self.inner.properties(),
        }
    }
    #[wasm_bindgen(js_name=substanceGroups,unchecked_return_type="SubstanceGroup[]")]
    pub fn substance_groups(&self) -> JsValue {
        self.inner
            .substance_groups()
            .into_iter()
            .map(|inner| JsValue::from(SubstanceGroup { inner }))
            .collect::<js_sys::Array>()
            .into()
    }
    #[wasm_bindgen(js_name=sourceCoordinateDim,unchecked_return_type="CoordinateDimension | null")]
    pub fn source_coordinate_dim(&self) -> JsValue {
        self.inner
            .source_coordinate_dim()
            .map_or(JsValue::NULL, |v| {
                JsValue::from(match v {
                    ck::CoordinateDimension::TwoD => 0u32,
                    ck::CoordinateDimension::ThreeD => 1u32,
                })
            })
    }
    #[wasm_bindgen(unchecked_return_type = "string | null")]
    pub fn title(&self) -> Result<JsValue, JsValue> {
        crate::host_values::optional_text(self.inner.title().as_ref())
    }
    pub fn index(&self) -> usize {
        self.inner.index()
    }
    #[wasm_bindgen(js_name=dataField,unchecked_return_type="string | null")]
    pub fn data_field(&self, name: &str) -> Result<JsValue, JsValue> {
        crate::host_values::optional_text(self.inner.data_field(name).as_ref())
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fromSdf)]
    pub fn from_sdf(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::from_sdf(text)
        cosmolkit_wasm::Molecule::from_sdf(text)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromSdfWithParams)]
    pub fn from_sdf_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::from_sdf_with_params(text,&params.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SdfReadParams| {
            result = Some(
                cosmolkit_wasm::Molecule::from_sdf_with_params(text, &p.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=fromMol)]
    pub fn from_mol(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::from_mol(text)
        cosmolkit_wasm::Molecule::from_mol(text)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromMolWithParams)]
    pub fn from_mol_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::from_mol_with_params(text,&params.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SdfReadParams| {
            result = Some(
                cosmolkit_wasm::Molecule::from_mol_with_params(text, &p.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=readMol)]
    pub fn read_mol(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::read_mol(text)
        cosmolkit_wasm::Molecule::read_mol(text)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=readMolWithParams)]
    pub fn read_mol_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::read_mol_with_params(text,&params.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SdfReadParams| {
            result = Some(
                cosmolkit_wasm::Molecule::read_mol_with_params(text, &p.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=readSdf)]
    pub fn read_sdf(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::read_sdf(text)
        cosmolkit_wasm::Molecule::read_sdf(text)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=readSdfWithParams)]
    pub fn read_sdf_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::read_sdf_with_params(text,&params.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SdfReadParams| {
            result = Some(
                cosmolkit_wasm::Molecule::read_sdf_with_params(text, &p.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
