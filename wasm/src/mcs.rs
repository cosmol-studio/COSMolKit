//! Canonical MCS query; retained query results share their native owner.
use crate::{
    fields::Field, input, mcs_errors::McsError, query_values::QueryGraph, storage::ReadonlyValue,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
crate::enums::enum_value!(McsAtomComparator,"McsAtomComparator"{AtomCompareAny="any",AtomCompareElements="elements",AtomCompareIsotopes="isotopes",AtomCompareAnyHeavyAtom="any_heavy_atom"});
crate::enums::enum_value!(McsBondComparator,"McsBondComparator"{BondCompareAny="any",BondCompareOrder="order",BondCompareOrderExact="order_exact"});
ck::__ck_configuration_McsAtomCompareParameters!(
    crate::project_parameters,
    ck::McsAtomCompareParameters,
    {},
    { McsAtomCompareParameters }
);
ck::__ck_configuration_McsBondCompareParameters!(
    crate::project_parameters,
    ck::McsBondCompareParameters,
    {},
    { McsBondCompareParameters }
);
ck::__ck_configuration_McsParameters!(crate::project_parameters, ck::McsParameters, {}, {
    McsParameters
});
#[wasm_bindgen]
pub struct McsResult(ReadonlyValue<ck::McsResult>);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ReadonlyMap<string, QueryGraph>")]
    pub type DegenerateQueries;
}
#[wasm_bindgen]
impl McsResult {
    #[wasm_bindgen(getter)]
    pub fn query(&self) -> Option<QueryGraph> {
        self.0.query.as_ref().map(|_| {
            QueryGraph(
                self.0
                    .child(|r| r.query.as_ref().expect("immutable present query")),
            )
        })
    }
    #[wasm_bindgen(getter,js_name=atomCount)]
    pub fn atom_count(&self) -> usize {
        self.0.atom_count
    }
    #[wasm_bindgen(getter,js_name=bondCount)]
    pub fn bond_count(&self) -> usize {
        self.0.bond_count
    }
    #[wasm_bindgen(getter)]
    pub fn completed(&self) -> bool {
        self.0.completed
    }
    #[wasm_bindgen(getter)]
    pub fn smarts(&self) -> Result<String, JsValue> {
        crate::property_errors::text(&self.0.smarts)
    }
    #[wasm_bindgen(getter)]
    pub fn degenerate(&self) -> Result<DegenerateQueries, JsValue> {
        let map = js_sys::Map::new();
        for key in self.0.degenerate.keys() {
            let text = crate::property_errors::text(key)?;
            let key = key.clone();
            let query = QueryGraph(
                self.0
                    .child(move |r| r.degenerate.get(&key).expect("immutable retained query")),
            );
            map.set(&text.into(), &query.into());
        }
        Ok(crate::fields::readonly_map(&map)?.unchecked_into())
    }
}
#[wasm_bindgen(js_name=maximumCommonSubstructure)]
pub fn maximum_common_substructure(
    inputs: crate::fragments::Molecules,
    params: Option<McsParametersInput>,
) -> Result<McsResult, JsValue> {
    let values = Vec::<ck::Molecule>::input(inputs)?;
    let refs = values.iter().collect::<Vec<_>>();
    McsParameters::with_input(params.as_ref(), |p| {
        ck::maximum_common_substructure_with_params(&refs, p)
            .map(|r| McsResult(ReadonlyValue::owned(r)))
            .map_err(McsError::from_rust)
    })
}

#[wasm_bindgen(js_name=maximumCommonSubstructureWithParams)]
pub fn maximum_common_substructure_with_params(
    inputs: crate::fragments::Molecules,
    params: McsParametersInput,
) -> Result<McsResult, JsValue> {
    maximum_common_substructure(inputs, Some(params))
}
