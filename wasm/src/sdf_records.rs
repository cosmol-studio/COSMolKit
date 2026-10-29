//! Finalized SDF records retain concrete/query alternatives and native owners.
use crate::{
    Molecule, input,
    io_configuration::{MolBlockWriteParams, MolBlockWriteParamsInput},
    molecular_io::{SdfReadParams, SdfReadParamsInput},
    molecular_io_errors::MolecularIoError,
    property_errors::text,
    query_values::QueryGraph,
    sdf_errors::SdfError,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use serde::Serialize;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;

#[derive(Serialize, Tsify)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum SdfGraph {
    Molecule {
        #[serde(with = "tsify::serde_wasm_bindgen::preserve")]
        #[tsify(type = "Molecule")]
        value: JsValue,
    },
    Query {
        #[serde(with = "tsify::serde_wasm_bindgen::preserve")]
        #[tsify(type = "QueryGraph")]
        value: JsValue,
    },
}
#[derive(Serialize, Tsify)]
#[serde(transparent)]
pub struct SdfDataFields(pub(crate) Vec<(String, String)>);
#[wasm_bindgen]
pub struct SdfRecord(ReadonlyValue<ck::SdfRecord>);
#[wasm_bindgen]
impl SdfRecord {
    #[wasm_bindgen(js_name=fromSdfWithParams)]
    pub fn from_sdf_with_params(
        text: Ts<input::TextInput>,
        params: SdfReadParamsInput,
    ) -> Result<Self, JsValue> {
        Self::from_sdf(text, Some(params))
    }
    #[wasm_bindgen(getter,js_name=substanceGroups)]
    pub fn substance_groups(&self) -> Result<crate::group_values::SubstanceGroupArray, JsValue> {
        use crate::fields::Field;
        Vec::<ck::SubstanceGroup>::get(crate::storage::NativeField::owned(
            self.0.substance_groups().to_vec(),
        ))
    }
    #[wasm_bindgen(js_name=fromQueryGraph)]
    pub fn from_query_graph(
        query: crate::query_values::QueryGraphInput,
        properties: crate::metadata_values::MoleculePropertiesInput,
    ) -> Result<Self, JsValue> {
        QueryGraph::with_ref(&query, |q| {
            crate::metadata_values::MoleculeProperties::with_ref(&properties, |p| {
                ck::SdfRecord::from_query_graph(q.clone(), p.clone())
                    .map(|v| Self(ReadonlyValue::owned(v)))
                    .map_err(SdfError::from_rust)
            })
        })
    }
    #[wasm_bindgen(getter)]
    pub fn properties(&self) -> crate::metadata_values::MoleculeProperties {
        crate::metadata_values::MoleculeProperties(self.0.properties().clone())
    }
    #[wasm_bindgen(js_name=fromSdf)]
    pub fn from_sdf(
        source: Ts<input::TextInput>,
        params: Option<SdfReadParamsInput>,
    ) -> Result<Self, JsValue> {
        let source = input::read::<String>(source)?;
        SdfReadParams::with_input(params.as_ref(), |p| {
            ck::SdfRecord::from_sdf_with_params(&source, p)
                .map(|v| Self(ReadonlyValue::owned(v)))
                .map_err(SdfError::from_rust)
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
    #[wasm_bindgen(getter)]
    pub fn index(&self) -> usize {
        self.0.index()
    }
    #[wasm_bindgen(getter)]
    pub fn title(&self) -> Result<Option<String>, JsValue> {
        self.0.title().map(text).transpose()
    }
    #[wasm_bindgen(js_name=dataField)]
    pub fn data_field(&self, name: Ts<input::TextInput>) -> Result<Option<String>, JsValue> {
        self.0
            .data_field(&input::read::<String>(name)?)
            .map(text)
            .transpose()
    }
    #[wasm_bindgen(getter,js_name=dataFields)]
    pub fn data_fields(&self) -> Result<Ts<SdfDataFields>, JsValue> {
        let values = self
            .0
            .data_fields()
            .iter()
            .map(|(k, v)| Ok((text(k)?, text(v)?)))
            .collect::<Result<Vec<_>, JsValue>>()?;
        SdfDataFields(values)
            .into_ts()
            .map_err(|e| js_sys::Error::new(&e.to_string()).into())
    }
    #[wasm_bindgen(getter,js_name=sourceCoordinateDim)]
    pub fn source_coordinate_dim(&self) -> Option<crate::enums::CoordinateDimension> {
        self.0.source_coordinate_dim().map(Into::into)
    }
    #[wasm_bindgen(getter)]
    pub fn molecule(&self) -> Result<Molecule, JsValue> {
        self.0
            .molecule()
            .map(|m| Molecule(m.copy()))
            .map_err(SdfError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=queryGraph)]
    pub fn query_graph(&self) -> Result<QueryGraph, JsValue> {
        self.0.query_graph().map_err(SdfError::from_rust)?;
        Ok(QueryGraph(self.0.child(|r| match r.graph() {
            ck::SdfGraph::Query(q) => q,
            _ => unreachable!("checked immutable graph variant"),
        })))
    }
    #[wasm_bindgen(getter)]
    pub fn graph(&self) -> Result<Ts<SdfGraph>, JsValue> {
        let value = match self.0.graph() {
            ck::SdfGraph::Molecule(m) => SdfGraph::Molecule {
                value: Molecule(m.copy()).into(),
            },
            ck::SdfGraph::Query(_) => SdfGraph::Query {
                value: self.query_graph()?.into(),
            },
        };
        value
            .into_ts()
            .map_err(|e| js_sys::Error::new(&e.to_string()).into())
    }
}
#[wasm_bindgen]
impl SdfRecord {
    #[wasm_bindgen(js_name=toMolWithParams)]
    pub fn to_mol_with_params(&self, params: MolBlockWriteParamsInput) -> Result<String, JsValue> {
        self.to_mol(Some(params))
    }
    #[wasm_bindgen(js_name=toSdfWithParams)]
    pub fn to_sdf_with_params(&self, params: MolBlockWriteParamsInput) -> Result<String, JsValue> {
        self.to_sdf(Some(params))
    }
}
