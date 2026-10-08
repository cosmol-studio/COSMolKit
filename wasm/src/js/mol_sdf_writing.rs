//! Complete immutable writer policy and canonical concrete/query record projections.
use crate::Molecule;
use crate::host_values::{bool_value, type_error, u32_value, usize_value};
use crate::io_errors::{io_error, sdf_error};
use crate::property_values::MoleculeProperties;
use crate::query_construction::QueryGraph;
use crate::sdf_reading::SdfRecord;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum SdfFormat {
    V2000,
    V3000,
}
#[wasm_bindgen]
pub struct MolCoordinateSelection {
    inner: ck::MolCoordinateSelection,
}
#[wasm_bindgen]
impl MolCoordinateSelection {
    #[wasm_bindgen(getter,js_name=Auto)]
    pub fn auto() -> Self {
        Self {
            inner: ck::MolCoordinateSelection::Auto,
        }
    }
    #[wasm_bindgen(js_name=TwoD)]
    pub fn two_d(
        #[wasm_bindgen(unchecked_param_type = "number")] id: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::MolCoordinateSelection::TwoD {
                id: usize_value(&id, "id")?,
            },
        })
    }
    #[wasm_bindgen(js_name=ThreeD)]
    pub fn three_d(
        #[wasm_bindgen(unchecked_param_type = "number")] id: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::MolCoordinateSelection::ThreeD {
                id: usize_value(&id, "id")?,
            },
        })
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        match self.inner {
            ck::MolCoordinateSelection::Auto => "Auto",
            ck::MolCoordinateSelection::TwoD { .. } => "TwoD",
            ck::MolCoordinateSelection::ThreeD { .. } => "ThreeD",
        }
        .into()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn id(&self) -> JsValue {
        match self.inner {
            ck::MolCoordinateSelection::Auto => JsValue::NULL,
            ck::MolCoordinateSelection::TwoD { id } | ck::MolCoordinateSelection::ThreeD { id } => {
                JsValue::from(id as u32)
            }
        }
    }
}
#[wasm_bindgen(
    inline_js = "export function visitMolSelection(v,f){try{f(v);}catch(cause){throw new TypeError('invalid MolCoordinateSelection',{cause});}} export function visitMolWriteParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid MolBlockWriteParams',{cause});}} export function visitMolQuery(v,f){try{f(v);}catch(cause){throw new TypeError('invalid QueryGraph',{cause});}} export function visitMolProperties(v,f){try{f(v);}catch(cause){throw new TypeError('invalid MoleculeProperties',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitMolSelection)]
    fn visit_selection(
        v: &JsValue,
        f: &mut dyn FnMut(&MolCoordinateSelection),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitMolWriteParams)]
    fn visit_params(v: &JsValue, f: &mut dyn FnMut(&MolBlockWriteParams)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitMolQuery)]
    fn visit_query(v: &JsValue, f: &mut dyn FnMut(&QueryGraph)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitMolProperties)]
    fn visit_properties(v: &JsValue, f: &mut dyn FnMut(&MoleculeProperties))
    -> Result<(), JsValue>;
}
#[wasm_bindgen]
pub struct MolBlockWriteParams {
    pub(crate) inner: ck::MolBlockWriteParams,
}
#[wasm_bindgen]
impl MolBlockWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "SdfFormat")] format: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] force_2d: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_stereo: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] kekulize: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] precision: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "MolCoordinateSelection")]
        coordinate_selection: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_coordinates: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::MolBlockWriteParams::default();
        if !format.is_undefined() {
            inner.format = match u32_value(&format, "format")? {
                0 => ck::SdfFormat::V2000,
                1 => ck::SdfFormat::V3000,
                _ => return Err(js_sys::RangeError::new("invalid SdfFormat").into()),
            };
        }
        if !force_2d.is_undefined() {
            inner.force_2d = bool_value(&force_2d, "force2d")?;
        }
        if !include_stereo.is_undefined() {
            inner.include_stereo = bool_value(&include_stereo, "includeStereo")?;
        }
        if !kekulize.is_undefined() {
            inner.kekulize = bool_value(&kekulize, "kekulize")?;
        }
        if !precision.is_undefined() {
            inner.precision = usize_value(&precision, "precision")?;
        }
        if !include_coordinates.is_undefined() {
            inner.include_coordinates = bool_value(&include_coordinates, "includeCoordinates")?;
        }
        if !coordinate_selection.is_undefined() {
            visit_selection(&coordinate_selection, &mut |s: &MolCoordinateSelection| {
                inner.coordinate_selection = s.inner;
            })?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter, unchecked_return_type = "SdfFormat")]
    pub fn format(&self) -> u32 {
        match self.inner.format {
            ck::SdfFormat::V2000 => 0,
            ck::SdfFormat::V3000 => 1,
        }
    }
    #[wasm_bindgen(getter,js_name=force2d)]
    pub fn force_2d(&self) -> bool {
        self.inner.force_2d
    }
    #[wasm_bindgen(getter,js_name=includeStereo)]
    pub fn include_stereo(&self) -> bool {
        self.inner.include_stereo
    }
    #[wasm_bindgen(getter)]
    pub fn kekulize(&self) -> bool {
        self.inner.kekulize
    }
    #[wasm_bindgen(getter)]
    pub fn precision(&self) -> usize {
        self.inner.precision
    }
    #[wasm_bindgen(getter,js_name=coordinateSelection)]
    pub fn coordinate_selection(&self) -> MolCoordinateSelection {
        MolCoordinateSelection {
            inner: self.inner.coordinate_selection,
        }
    }
    #[wasm_bindgen(getter,js_name=includeCoordinates)]
    pub fn include_coordinates(&self) -> bool {
        self.inner.include_coordinates
    }
}
#[wasm_bindgen]
impl SdfRecord {
    #[wasm_bindgen(js_name=fromQueryGraph)]
    pub fn from_query_graph(
        #[wasm_bindgen(unchecked_param_type = "QueryGraph")] query: JsValue,
        #[wasm_bindgen(unchecked_param_type = "MoleculeProperties")] properties: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::SdfRecord::from_query_graph(query.inner.clone(),properties.inner.clone())
        let mut result = None;
        visit_query(&query, &mut |q: &QueryGraph| {
            let visit = visit_properties(&properties, &mut |p: &MoleculeProperties| {
                result = Some(
                    cosmolkit_wasm::SdfRecord::from_query_graph(q.inner.clone(), p.inner.clone())
                        .map(|inner| Self { inner })
                        .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e)),
                );
            });
            if let Err(e) = visit {
                result = Some(Err(e));
            }
        })?;
        result.ok_or_else(|| type_error("query"))?
    }
    #[wasm_bindgen(js_name=toMol)]
    pub fn to_mol(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_mol()
        self.inner
            .to_mol()
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toMolWithParams)]
    pub fn to_mol_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "MolBlockWriteParams")] params: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_mol_with_params(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &MolBlockWriteParams| {
            result = Some(
                self.inner
                    .to_mol_with_params(&p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=toSdf)]
    pub fn to_sdf(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_sdf()
        self.inner
            .to_sdf()
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toSdfWithParams)]
    pub fn to_sdf_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "MolBlockWriteParams")] params: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_sdf_with_params(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &MolBlockWriteParams| {
            result = Some(
                self.inner
                    .to_sdf_with_params(&p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=toMol)]
    pub fn to_mol(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_mol()
        self.inner
            .to_mol()
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toMolWithParams)]
    pub fn to_mol_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "MolBlockWriteParams")] params: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_mol_with_params(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &MolBlockWriteParams| {
            result = Some(
                self.inner
                    .to_mol_with_params(&p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=toSdf)]
    pub fn to_sdf(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_sdf()
        self.inner
            .to_sdf()
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toSdfWithParams)]
    pub fn to_sdf_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "MolBlockWriteParams")] params: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_sdf_with_params(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &MolBlockWriteParams| {
            result = Some(
                self.inner
                    .to_sdf_with_params(&p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=toSdf2d)]
    pub fn to_sdf_2d(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_sdf_2d()
        self.inner
            .to_sdf_2d()
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toSdf2dWithParams)]
    pub fn to_sdf_2d_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "MolBlockWriteParams")] params: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_sdf_2d_with_params(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &MolBlockWriteParams| {
            result = Some(
                self.inner
                    .to_sdf_2d_with_params(&p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=toSdf3d)]
    pub fn to_sdf_3d(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_sdf_3d()
        self.inner
            .to_sdf_3d()
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toSdf3dWithParams)]
    pub fn to_sdf_3d_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "MolBlockWriteParams")] params: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_sdf_3d_with_params(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &MolBlockWriteParams| {
            result = Some(
                self.inner
                    .to_sdf_3d_with_params(&p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=writeMol)]
    pub fn write_mol(&self, path: &str) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.write_mol(path)
        self.inner
            .write_mol(path)
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writeMolWithParams)]
    pub fn write_mol_with_params(
        &self,
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "MolBlockWriteParams")] params: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.write_mol_with_params(path,&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &MolBlockWriteParams| {
            result = Some(
                self.inner
                    .write_mol_with_params(path, &p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=writeSdf)]
    pub fn write_sdf(&self, path: &str) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.write_sdf(path)
        self.inner
            .write_sdf(path)
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writeSdfWithParams)]
    pub fn write_sdf_with_params(
        &self,
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "MolBlockWriteParams")] params: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.write_sdf_with_params(path,&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &MolBlockWriteParams| {
            result = Some(
                self.inner
                    .write_sdf_with_params(path, &p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=writeSdfFiles)]
    pub fn write_sdf_files(
        &self,
        directory: &str,
        #[wasm_bindgen(unchecked_param_type = "string | null")] file_name: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.write_sdf_files(directory,file_name.as_deref())
        let file_name = optional_file_name(&file_name)?;
        self.inner
            .write_sdf_files(directory, file_name.as_deref())
            .map(|p| p.display().to_string())
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writeSdfFilesWithParams)]
    pub fn write_sdf_files_with_params(
        &self,
        directory: &str,
        #[wasm_bindgen(unchecked_param_type = "string | null")] file_name: JsValue,
        #[wasm_bindgen(unchecked_param_type = "MolBlockWriteParams")] params: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.write_sdf_files_with_params(directory,file_name.as_deref(),&p.inner)
        let file_name = optional_file_name(&file_name)?;
        let mut result = None;
        visit_params(&params, &mut |p: &MolBlockWriteParams| {
            result = Some(
                self.inner
                    .write_sdf_files_with_params(directory, file_name.as_deref(), &p.inner)
                    .map(|p| p.display().to_string())
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
fn optional_file_name(value: &JsValue) -> Result<Option<String>, JsValue> {
    if value.is_null() || value.is_undefined() {
        Ok(None)
    } else {
        value
            .as_string()
            .map(Some)
            .ok_or_else(|| type_error("fileName"))
    }
}
