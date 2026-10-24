//! Query construction and serialization use canonical facade entrypoints.
use crate::{Molecule, input, query_values::QueryGraph, storage::ReadonlyValue};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen(js_name=writeSmarts)]
pub fn write_smarts(
    query: crate::query_values::QueryGraphInput,
    params: Option<SmartsWriteParamsInput>,
) -> Result<String, JsValue> {
    QueryGraph::with_ref(&query, |query| {
        SmartsWriteParams::with_input(params.as_ref(), |p| {
            ck::write_smarts(query, p)
                .map_err(crate::smarts_errors::SmartsWriteError::from_rust)
                .and_then(|v| crate::property_errors::text(&v))
        })
    })
}
#[wasm_bindgen(js_name=writeCxSmarts)]
pub fn write_cx_smarts(
    query: crate::query_values::QueryGraphInput,
    params: Option<SmartsWriteParamsInput>,
) -> Result<String, JsValue> {
    QueryGraph::with_ref(&query, |query| {
        SmartsWriteParams::with_input(params.as_ref(), |p| {
            ck::write_cx_smarts(query, p)
                .map_err(crate::smarts_errors::SmartsWriteError::from_rust)
                .and_then(|v| crate::property_errors::text(&v))
        })
    })
}
#[wasm_bindgen(js_name=parseSmartsWithParams)]
pub fn parse_smarts_with_params(
    text: tsify::Ts<input::TextInput>,
    params: SmartsParseParamsInput,
) -> Result<QueryGraph, JsValue> {
    QueryGraph::from_smarts(text, Some(params))
}
#[wasm_bindgen]
impl QueryGraph {
    #[wasm_bindgen(js_name=fromSmartsWithParams)]
    pub fn from_smarts_with_params(
        text: tsify::Ts<input::TextInput>,
        params: SmartsParseParamsInput,
    ) -> Result<QueryGraph, JsValue> {
        Self::from_smarts(text, Some(params))
    }
}
#[wasm_bindgen]
pub struct CompiledQuery(pub(crate) ReadonlyValue<ck::CompiledQuery>);
impl CompiledQuery {
    pub(crate) fn from_native(value: ck::CompiledQuery) -> Self {
        Self(ReadonlyValue::owned(value))
    }
}
crate::fields::value_field!(
    ck::CompiledQuery,
    CompiledQuery,
    CompiledQueryInput,
    "CompiledQuery",
    CompiledQuery::from_native
);
#[wasm_bindgen]
impl CompiledQuery {
    #[wasm_bindgen(getter,js_name=numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.0.num_atoms()
    }
    #[wasm_bindgen(getter,js_name=numBonds)]
    pub fn num_bonds(&self) -> usize {
        self.0.num_bonds()
    }
    #[wasm_bindgen(getter,js_name=atomOrder)]
    pub fn atom_order(&self) -> Result<<Vec<usize> as crate::fields::Field>::Output, JsValue> {
        <Vec<usize> as crate::fields::Field>::get(crate::storage::NativeField::owned(
            self.0.atom_order().to_vec(),
        ))
    }
    #[wasm_bindgen(getter)]
    pub fn query(&self) -> QueryGraph {
        QueryGraph(self.0.child(ck::CompiledQuery::query))
    }
}
#[wasm_bindgen(js_name=compileQuery)]
pub fn compile_query(
    query: crate::query_values::QueryGraphInput,
) -> Result<CompiledQuery, JsValue> {
    QueryGraph::with_ref(&query, |query| {
        ck::compile_query(query)
            .map(|value| CompiledQuery(ReadonlyValue::owned(value)))
            .map_err(crate::maccs::QueryCompileError::from_rust)
    })
}
ck::__ck_configuration_SmartsParseParams!(crate::project_parameters, ck::SmartsParseParams, {}, {
    SmartsParseParams
});
ck::__ck_configuration_SmartsWriteParams!(crate::project_parameters, ck::SmartsWriteParams, {}, {
    SmartsWriteParams
});
#[wasm_bindgen]
impl QueryGraph {
    #[wasm_bindgen(js_name=fromSmarts)]
    pub fn from_smarts(
        text: tsify::Ts<input::TextInput>,
        params: Option<SmartsParseParamsInput>,
    ) -> Result<QueryGraph, JsValue> {
        let text = input::read::<String>(text)?;
        SmartsParseParams::with_input(params.as_ref(), |p| {
            ck::search::from_smarts_with_params(&text, p)
                .map(|v| Self(ReadonlyValue::owned(v)))
                .map_err(crate::smarts_parse_errors::SmartsParseError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toSmarts)]
    pub fn to_smarts(&self, params: Option<SmartsWriteParamsInput>) -> Result<String, JsValue> {
        SmartsWriteParams::with_input(params.as_ref(), |p| {
            ck::write_smarts(&self.0, p)
                .map_err(crate::smarts_errors::SmartsWriteError::from_rust)
                .and_then(|v| crate::property_errors::text(&v))
        })
    }
    #[wasm_bindgen(js_name=toCxSmarts)]
    pub fn to_cx_smarts(&self, params: Option<SmartsWriteParamsInput>) -> Result<String, JsValue> {
        SmartsWriteParams::with_input(params.as_ref(), |p| {
            ck::write_cx_smarts(&self.0, p)
                .map_err(crate::smarts_errors::SmartsWriteError::from_rust)
                .and_then(|v| crate::property_errors::text(&v))
        })
    }
    #[wasm_bindgen(js_name=fingerprintPattern)]
    pub fn fingerprint_pattern(
        &self,
        params: Option<crate::pattern::PatternFingerprintParamsInput>,
    ) -> Result<crate::fingerprint_values::Fingerprint, JsValue> {
        crate::pattern::PatternFingerprintParams::with_input(params.as_ref(), |p| {
            ck::fingerprint_pattern_query_with_params(&self.0, p)
                .map(crate::fingerprint_values::Fingerprint)
                .map_err(crate::pattern::PatternFingerprintError::from_rust)
        })
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=toSmarts)]
    pub fn to_smarts(&self, params: Option<SmartsWriteParamsInput>) -> Result<String, JsValue> {
        SmartsWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_smarts_with_params(p)
                .map_err(crate::smarts_errors::SmartsWriteError::from_rust)
                .and_then(|v| crate::property_errors::text(&v))
        })
    }
    #[wasm_bindgen(js_name=toCxSmarts)]
    pub fn to_cx_smarts(&self, params: Option<SmartsWriteParamsInput>) -> Result<String, JsValue> {
        SmartsWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_cx_smarts_with_params(p)
                .map_err(crate::smarts_errors::SmartsWriteError::from_rust)
                .and_then(|v| crate::property_errors::text(&v))
        })
    }
}
#[wasm_bindgen(js_name=parseSmarts)]
pub fn parse_smarts(
    text: tsify::Ts<input::TextInput>,
    params: Option<SmartsParseParamsInput>,
) -> Result<QueryGraph, JsValue> {
    QueryGraph::from_smarts(text, params)
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = toSmartsWithParams)]
    pub fn to_smarts_with_params(&self, params: SmartsWriteParamsInput) -> Result<String, JsValue> {
        self.to_smarts(Some(params))
    }
    #[wasm_bindgen(js_name = toCxSmartsWithParams)]
    pub fn to_cx_smarts_with_params(
        &self,
        params: SmartsWriteParamsInput,
    ) -> Result<String, JsValue> {
        self.to_cx_smarts(Some(params))
    }
}
