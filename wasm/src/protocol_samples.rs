//! Development-only witnesses for shared binding protocols.
use super::*;
// These are protocol samples, not new chemistry API or duplicated native Params.
#[derive(Serialize, Deserialize, Tsify)]
#[serde(rename_all = "snake_case")]
pub enum Mode {
    Preserve,
    Replace,
}

#[derive(Serialize, Deserialize, Tsify)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Outcome {
    Completed { count: u32 },
    Failed { message: String },
}

#[derive(Serialize, Deserialize, Tsify)]
#[serde(deny_unknown_fields)]
pub struct ArrayInput {
    #[serde(deserialize_with = "input::exact_array")]
    pub point: [f64; 3],
    pub indices: Vec<u32>,
    pub mode: Mode,
}

#[wasm_bindgen]
pub fn roundtrip(input: Ts<ArrayInput>) -> Result<Ts<ArrayInput>, JsValue> {
    Ok(input::decode(input)?.into_ts().map_err(JsError::from)?)
}

#[wasm_bindgen]
pub fn outcome(input: Ts<Outcome>) -> Result<Ts<Outcome>, JsValue> {
    Ok(input::decode(input)?.into_ts().map_err(JsError::from)?)
}

/// Protocol fixture exercises the same native cause projection as real queries.
#[wasm_bindgen]
pub fn raise() -> Result<(), JsValue> {
    Err(operation_errors::OperationError::from_rust(
        ck::OperationError::Rings(ck::RingFindingError::InvalidTopology(
            ck::TopologyValidationError::TemplateAttachmentOrder {
                atom: ck::AtomId::new(1),
                source: ck::TemplateAttachmentOrderError::TargetRemoved {
                    position: 2,
                    target: ck::AtomId::new(3),
                },
            },
        )),
    ))
}

/// Owned typed-array output; mutation in JS cannot modify native storage.
#[wasm_bindgen]
pub fn bytes(input: Ts<input::BytesInput>) -> Result<Vec<u8>, JsValue> {
    input::read(input)
}

/// Explicit 64-bit integer boundary: wasm-bindgen generates bigint, not number.
#[wasm_bindgen]
pub fn integer(input: Ts<input::U64Input>) -> Result<u64, JsValue> {
    input::read(input)
}

// Protocol-only native records: the same Configuration-derived projection used
// by business bindings, not bespoke implementations that evade the framework.
#[derive(Clone, Default, ck::Configuration)]
pub struct ProtocolSelection {
    #[configuration(required)]
    pub indices: Vec<usize>,
}
#[derive(Clone, Default, ck::Configuration)]
pub struct ProtocolNested {
    pub child: ck::RingSearchParams,
    pub optional_child: Option<ck::RingSearchParams>,
}
__ck_configuration_ProtocolSelection!(crate::project_parameters, ProtocolSelection, {}, {
    Selection
});
__ck_configuration_ProtocolNested!(crate::project_parameters, ProtocolNested, {}, {
    NestedParams
});
