//! Native path fingerprints, their configuration and requested provenance.
use crate::{
    Molecule,
    errors::child,
    errors::{enum_error, exception},
    fingerprint_values::Fingerprint,
    input,
    query_values::{QueryGraph, QueryGraphInput},
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | QueryGraphError | RingFindingError | PathError | FingerprintError"
    )]
    #[derive(Clone)]
    pub type LayeredCause;
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | QueryGraphError | LayeredFingerprintError | FingerprintError"
    )]
    #[derive(Clone)]
    pub type TopologicalCause;
    #[wasm_bindgen(typescript_type = "ReadonlyArray<ReadonlyArray<number>>")]
    pub type AtomBits;
    #[wasm_bindgen(typescript_type = "ReadonlyMap<number, ReadonlyArray<ReadonlyArray<number>>>")]
    pub type PathBitInfo;
}
fn layered_cause(owner: &ReadonlyValue<ck::LayeredFingerprintError>) -> Option<LayeredCause> {
    use ck::LayeredFingerprintError as E;
    Some(match owner.as_ref() {
        E::InvalidArguments { .. } => return None,
        E::Topology(_) => child!(
            owner,
            E::Topology,
            crate::model_errors::TopologyValidationError
        ),
        E::Query(_) => child!(owner, E::Query, crate::query_model_errors::QueryGraphError),
        E::Rings(_) => child!(owner, E::Rings, crate::ring_errors::RingFindingError),
        E::Paths(_) => child!(owner, E::Paths, crate::path_errors::PathError),
        E::Value(_) => child!(owner, E::Value, crate::fingerprint_errors::FingerprintError),
    })
}
enum_error!(LayeredFingerprintError,LayeredFingerprintErrorKind,ck::LayeredFingerprintError,{
    InvalidArguments=>ck::LayeredFingerprintError::InvalidArguments{..},Topology=>ck::LayeredFingerprintError::Topology(_),Query=>ck::LayeredFingerprintError::Query(_),Rings=>ck::LayeredFingerprintError::Rings(_),Paths=>ck::LayeredFingerprintError::Paths(_),Value=>ck::LayeredFingerprintError::Value(_)
},{
    #[wasm_bindgen(getter)] pub fn reason(&self)->Option<String>{match self.inner.as_ref(){ck::LayeredFingerprintError::InvalidArguments{reason}=>Some((*reason).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<LayeredCause>{layered_cause(&self.inner)}
},source=LayeredCause,cause=layered_cause);
fn topological_cause(
    owner: &ReadonlyValue<ck::TopologicalFingerprintError>,
) -> Option<TopologicalCause> {
    use ck::TopologicalFingerprintError as E;
    Some(match owner.as_ref() {
        E::InvalidArguments { .. } | E::OutputNotRequested { .. } => return None,
        E::Topology(_) => child!(
            owner,
            E::Topology,
            crate::model_errors::TopologyValidationError
        ),
        E::Query(_) => child!(owner, E::Query, crate::query_model_errors::QueryGraphError),
        E::Paths(_) => child!(owner, E::Paths, LayeredFingerprintError),
        E::Value(_) => child!(owner, E::Value, crate::fingerprint_errors::FingerprintError),
    })
}
enum_error!(TopologicalFingerprintError,TopologicalFingerprintErrorKind,ck::TopologicalFingerprintError,{
    InvalidArguments=>ck::TopologicalFingerprintError::InvalidArguments{..},OutputNotRequested=>ck::TopologicalFingerprintError::OutputNotRequested{..},Topology=>ck::TopologicalFingerprintError::Topology(_),Query=>ck::TopologicalFingerprintError::Query(_),Paths=>ck::TopologicalFingerprintError::Paths(_),Value=>ck::TopologicalFingerprintError::Value(_)
},{
    #[wasm_bindgen(getter)] pub fn reason(&self)->Option<String>{match self.inner.as_ref(){ck::TopologicalFingerprintError::InvalidArguments{reason}=>Some((*reason).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn field(&self)->Option<String>{match self.inner.as_ref(){ck::TopologicalFingerprintError::OutputNotRequested{field}=>Some((*field).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<TopologicalCause>{topological_cause(&self.inner)}
},source=TopologicalCause,cause=topological_cause);

#[wasm_bindgen]
pub struct LayeredFingerprintLayers(pub(crate) ck::LayeredFingerprintLayers);
crate::fields::value_field!(
    ck::LayeredFingerprintLayers,
    LayeredFingerprintLayers,
    LayeredFingerprintLayersInput,
    "LayeredFingerprintLayers"
);
#[wasm_bindgen]
impl LayeredFingerprintLayers {
    #[wasm_bindgen(js_name=fromBitsRetain)]
    pub fn from_bits_retain(bits: tsify::Ts<input::U32Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::LayeredFingerprintLayers::from_bits_retain(
            input::read(bits)?,
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn bits(&self) -> u32 {
        self.0.bits()
    }
    pub fn contains(&self, other: LayeredFingerprintLayersInput) -> Result<bool, JsValue> {
        Self::with_ref(&other, |other| Ok(self.0.contains(*other)))
    }
    pub fn union(&self, other: LayeredFingerprintLayersInput) -> Result<Self, JsValue> {
        Self::with_ref(&other, |other| Ok(Self(self.0 | *other)))
    }
}
macro_rules! layer_constants {($($name:ident),*)=>{#[wasm_bindgen] impl LayeredFingerprintLayers{$(#[wasm_bindgen(getter,js_name=$name)] pub fn $name()->Self{Self(ck::LayeredFingerprintLayers::$name)})*}};}
#[allow(non_snake_case)]
const _: () = {
    layer_constants!(
        TOPOLOGY,
        BOND_ORDER,
        ATOM_TYPE,
        RING_PRESENCE,
        RING_SIZE,
        AROMATICITY,
        ACTIVE,
        SUBSTRUCTURE,
        ALL_SOURCE_BITS
    );
};
ck::__ck_configuration_LayeredFingerprintParams!(
    crate::project_parameters,
    ck::LayeredFingerprintParams,
    {},
    { LayeredFingerprintParams }
);
ck::__ck_configuration_TopologicalFingerprintParams!(
    crate::project_parameters,
    ck::TopologicalFingerprintParams,
    {},
    { TopologicalFingerprintParams }
);
ck::__ck_configuration_TopologicalFingerprintOutputRequest!(
    crate::project_parameters,
    ck::TopologicalFingerprintOutputRequest,
    {},
    { TopologicalFingerprintOutputRequest }
);

fn array(values: impl IntoIterator<Item = JsValue>) -> JsValue {
    let a = values.into_iter().collect::<js_sys::Array>();
    js_sys::Object::freeze(&a);
    a.into()
}
#[wasm_bindgen]
pub struct LayeredFingerprintResult(pub(crate) ck::LayeredFingerprintResult);
crate::fields::value_field!(
    ck::LayeredFingerprintResult,
    LayeredFingerprintResult,
    LayeredFingerprintResultInput,
    "LayeredFingerprintResult"
);
crate::fields::optional_value_array!(
    ck::LayeredFingerprintResult,
    LayeredFingerprintResult,
    LayeredFingerprintResults,
    "ReadonlyArray<LayeredFingerprintResult | undefined>"
);
#[wasm_bindgen]
impl LayeredFingerprintResult {
    #[wasm_bindgen(getter)]
    pub fn fingerprint(&self) -> Fingerprint {
        Fingerprint(self.0.fingerprint().clone())
    }
    #[wasm_bindgen(getter,js_name=atomCounts)]
    pub fn atom_counts(&self) -> Option<crate::fingerprint_output::AtomCounts> {
        self.0
            .atom_counts()
            .map(|v| array(v.iter().copied().map(JsValue::from)).unchecked_into())
    }
}
#[wasm_bindgen]
pub struct TopologicalFingerprintResult(ck::TopologicalFingerprintResult);

#[derive(serde::Deserialize, tsify::Tsify)]
#[serde(transparent)]
pub struct AtomBitRows(#[tsify(type = "ReadonlyArray<ReadonlyArray<number>>")] Vec<Vec<u32>>);
#[derive(serde::Deserialize, tsify::Tsify)]
#[serde(transparent)]
pub struct BitPathRows(
    #[tsify(type = "ReadonlyMap<number, ReadonlyArray<ReadonlyArray<number>>>")]
    std::collections::BTreeMap<u32, Vec<Vec<i32>>>,
);
#[wasm_bindgen]
pub struct TopologicalFingerprintOutput(ck::TopologicalFingerprintOutput);
crate::fields::value_field!(
    ck::TopologicalFingerprintOutput,
    TopologicalFingerprintOutput,
    TopologicalFingerprintOutputInput,
    "TopologicalFingerprintOutput"
);
#[wasm_bindgen]
impl TopologicalFingerprintOutput {
    #[wasm_bindgen(constructor)]
    pub fn new(
        atom_bits: Option<tsify::Ts<AtomBitRows>>,
        bit_info: Option<tsify::Ts<BitPathRows>>,
    ) -> Result<Self, JsValue> {
        let atom_bits = atom_bits
            .map(|v| input::decode(v).map(|v| v.0))
            .transpose()?;
        let bit_info = bit_info
            .map(|v| {
                if !v.js_value().is_instance_of::<js_sys::Map>() {
                    return Err(js_sys::TypeError::new("expected Map").into());
                }
                input::decode(v).map(|v| v.0)
            })
            .transpose()?;
        Ok(Self(ck::TopologicalFingerprintOutput {
            atom_bits,
            bit_info,
        }))
    }
    #[wasm_bindgen(getter,js_name=atomBits)]
    pub fn atom_bits(&self) -> Option<AtomBits> {
        self.0.atom_bits.as_ref().map(|rows| {
            array(
                rows.iter()
                    .map(|row| array(row.iter().copied().map(JsValue::from))),
            )
            .unchecked_into()
        })
    }
    #[wasm_bindgen(getter,js_name=bitInfo)]
    pub fn bit_info(&self) -> Result<Option<PathBitInfo>, JsValue> {
        self.0
            .bit_info
            .as_ref()
            .map(|rows| {
                let map = js_sys::Map::new();
                for (&key, paths) in rows {
                    map.set(
                        &key.into(),
                        &array(
                            paths
                                .iter()
                                .map(|row| array(row.iter().copied().map(JsValue::from))),
                        ),
                    );
                }
                Ok(crate::fields::readonly_map(&map)?.unchecked_into())
            })
            .transpose()
    }
}
#[wasm_bindgen]
impl TopologicalFingerprintResult {
    #[wasm_bindgen(constructor)]
    pub fn new(
        fingerprint: crate::fingerprint_values::FingerprintInput,
        output: TopologicalFingerprintOutputInput,
    ) -> Result<Self, JsValue> {
        use crate::fields::Field;
        Ok(Self(ck::TopologicalFingerprintResult {
            fingerprint: ck::Fingerprint::input(fingerprint)?,
            output: ck::TopologicalFingerprintOutput::input(output)?,
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn fingerprint(&self) -> Fingerprint {
        Fingerprint(self.0.fingerprint().clone())
    }
    #[wasm_bindgen(getter,js_name=atomBits)]
    pub fn atom_bits(&self) -> Result<AtomBits, JsValue> {
        self.0
            .atom_bits()
            .map(|rows| {
                array(
                    rows.iter()
                        .map(|row| array(row.iter().copied().map(JsValue::from))),
                )
                .unchecked_into()
            })
            .map_err(TopologicalFingerprintError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=bitInfo)]
    pub fn bit_info(&self) -> Result<PathBitInfo, JsValue> {
        let rows = self
            .0
            .bit_info()
            .map_err(TopologicalFingerprintError::from_rust)?;
        let map = js_sys::Map::new();
        for (&key, paths) in rows {
            map.set(
                &key.into(),
                &array(
                    paths
                        .iter()
                        .map(|row| array(row.iter().copied().map(JsValue::from))),
                ),
            );
        }
        Ok(crate::fields::readonly_map(&map)?.unchecked_into())
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintLayered)]
    pub fn fingerprint_layered(
        &self,
        params: Option<LayeredFingerprintParamsInput>,
    ) -> Result<Fingerprint, JsValue> {
        LayeredFingerprintParams::with_input(params.as_ref(), |p| {
            self.0
                .fingerprint_layered_with_params(p)
                .map(Fingerprint)
                .map_err(LayeredFingerprintError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=fingerprintLayeredWithOutput)]
    pub fn fingerprint_layered_with_output(
        &self,
        params: Option<LayeredFingerprintParamsInput>,
    ) -> Result<LayeredFingerprintResult, JsValue> {
        LayeredFingerprintParams::with_input(params.as_ref(), |p| {
            self.0
                .fingerprint_layered_with_output_with_params(p)
                .map(LayeredFingerprintResult)
                .map_err(LayeredFingerprintError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=fingerprintTopological)]
    pub fn fingerprint_topological(
        &self,
        params: Option<TopologicalFingerprintParamsInput>,
    ) -> Result<Fingerprint, JsValue> {
        TopologicalFingerprintParams::with_input(params.as_ref(), |p| {
            self.0
                .fingerprint_topological_with_params(p)
                .map(Fingerprint)
                .map_err(TopologicalFingerprintError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalWithOutput)]
    pub fn fingerprint_topological_with_output(
        &self,
        params: Option<TopologicalFingerprintParamsInput>,
        additional_output: Option<TopologicalFingerprintOutputRequestInput>,
    ) -> Result<TopologicalFingerprintResult, JsValue> {
        TopologicalFingerprintParams::with_input(params.as_ref(), |p| {
            TopologicalFingerprintOutputRequest::with_input(additional_output.as_ref(), |r| {
                self.0
                    .fingerprint_topological_with_output_with_params(p, *r)
                    .map(TopologicalFingerprintResult)
                    .map_err(TopologicalFingerprintError::from_rust)
            })
        })
    }
}
#[wasm_bindgen(js_name=fingerprintLayeredQuery)]
pub fn fingerprint_layered_query(
    query: QueryGraphInput,
    params: Option<LayeredFingerprintParamsInput>,
) -> Result<Fingerprint, JsValue> {
    QueryGraph::with_ref(&query, |q| {
        LayeredFingerprintParams::with_input(params.as_ref(), |p| {
            ck::fingerprint_layered_query_with_params(q, p)
                .map(Fingerprint)
                .map_err(LayeredFingerprintError::from_rust)
        })
    })
}
#[wasm_bindgen(js_name=fingerprintLayeredQueryWithOutput)]
pub fn fingerprint_layered_query_with_output(
    query: QueryGraphInput,
    params: Option<LayeredFingerprintParamsInput>,
) -> Result<LayeredFingerprintResult, JsValue> {
    QueryGraph::with_ref(&query, |q| {
        LayeredFingerprintParams::with_input(params.as_ref(), |p| {
            ck::fingerprint_layered_query_with_output_with_params(q, p)
                .map(LayeredFingerprintResult)
                .map_err(LayeredFingerprintError::from_rust)
        })
    })
}
#[wasm_bindgen(js_name=fingerprintTopologicalQuery)]
pub fn fingerprint_topological_query(
    query: QueryGraphInput,
    params: Option<TopologicalFingerprintParamsInput>,
) -> Result<Fingerprint, JsValue> {
    QueryGraph::with_ref(&query, |q| {
        TopologicalFingerprintParams::with_input(params.as_ref(), |p| {
            ck::fingerprint_topological_query_with_params(q, p)
                .map(Fingerprint)
                .map_err(TopologicalFingerprintError::from_rust)
        })
    })
}
#[wasm_bindgen(js_name=fingerprintTopologicalQueryWithOutput)]
pub fn fingerprint_topological_query_with_output(
    query: QueryGraphInput,
    params: Option<TopologicalFingerprintParamsInput>,
    additional_output: Option<TopologicalFingerprintOutputRequestInput>,
) -> Result<TopologicalFingerprintResult, JsValue> {
    QueryGraph::with_ref(&query, |q| {
        TopologicalFingerprintParams::with_input(params.as_ref(), |p| {
            TopologicalFingerprintOutputRequest::with_input(additional_output.as_ref(), |r| {
                ck::fingerprint_topological_query_with_output_with_params(q, p, *r)
                    .map(TopologicalFingerprintResult)
                    .map_err(TopologicalFingerprintError::from_rust)
            })
        })
    })
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = fingerprintLayeredWithParams)]
    pub fn fingerprint_layered_with_params(
        &self,
        params: LayeredFingerprintParamsInput,
    ) -> Result<Fingerprint, JsValue> {
        self.fingerprint_layered(Some(params))
    }
    #[wasm_bindgen(js_name = fingerprintTopologicalWithParams)]
    pub fn fingerprint_topological_with_params(
        &self,
        params: TopologicalFingerprintParamsInput,
    ) -> Result<Fingerprint, JsValue> {
        self.fingerprint_topological(Some(params))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintLayeredWithOutputWithParams)]
    pub fn fingerprint_layered_with_output_with_params(
        &self,
        params: LayeredFingerprintParamsInput,
    ) -> Result<LayeredFingerprintResult, JsValue> {
        self.fingerprint_layered_with_output(Some(params))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintTopologicalWithOutputWithParams)]
    pub fn fingerprint_topological_with_output_with_params(
        &self,
        params: TopologicalFingerprintParamsInput,
        additional_output: TopologicalFingerprintOutputRequestInput,
    ) -> Result<TopologicalFingerprintResult, JsValue> {
        self.fingerprint_topological_with_output(Some(params), Some(additional_output))
    }
}
#[wasm_bindgen(js_name=fingerprintLayeredQueryWithParams)]
pub fn fingerprint_layered_query_with_params(
    query: QueryGraphInput,
    params: LayeredFingerprintParamsInput,
) -> Result<Fingerprint, JsValue> {
    fingerprint_layered_query(query, Some(params))
}

#[wasm_bindgen(js_name=fingerprintLayeredQueryWithOutputWithParams)]
pub fn fingerprint_layered_query_with_output_with_params(
    query: QueryGraphInput,
    params: LayeredFingerprintParamsInput,
) -> Result<LayeredFingerprintResult, JsValue> {
    fingerprint_layered_query_with_output(query, Some(params))
}

#[wasm_bindgen(js_name=fingerprintTopologicalQueryWithParams)]
pub fn fingerprint_topological_query_with_params(
    query: QueryGraphInput,
    params: TopologicalFingerprintParamsInput,
) -> Result<Fingerprint, JsValue> {
    fingerprint_topological_query(query, Some(params))
}

#[wasm_bindgen(js_name=fingerprintTopologicalQueryWithOutputWithParams)]
pub fn fingerprint_topological_query_with_output_with_params(
    query: QueryGraphInput,
    params: TopologicalFingerprintParamsInput,
    additional_output: TopologicalFingerprintOutputRequestInput,
) -> Result<TopologicalFingerprintResult, JsValue> {
    fingerprint_topological_query_with_output(query, Some(params), Some(additional_output))
}
