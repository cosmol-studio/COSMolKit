//! Detached canonical group values; no molecular runtime or group algorithms.
use crate::host_values::{sequence, type_error, usize_value};
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
fn point(value: &JsValue) -> Result<[f64; 3], JsValue> {
    let a = sequence(value, "point")?;
    if a.length() != 3 {
        return Err(type_error("three coordinates"));
    }
    Ok([
        a.get(0).as_f64().ok_or_else(|| type_error("point[0]"))?,
        a.get(1).as_f64().ok_or_else(|| type_error("point[1]"))?,
        a.get(2).as_f64().ok_or_else(|| type_error("point[2]"))?,
    ])
}
fn vector(value: &[f64; 3]) -> JsValue {
    value
        .iter()
        .map(|v| JsValue::from(*v))
        .collect::<js_sys::Array>()
        .into()
}
#[wasm_bindgen]
pub struct SubstanceGroupId {
    pub(crate) inner: ck::SubstanceGroupId,
}
#[wasm_bindgen]
impl SubstanceGroupId {
    pub fn new(
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::SubstanceGroupId::new(usize_value(&index, "index")?),
        })
    }
    pub fn index(&self) -> usize {
        self.inner.index()
    }
}
#[wasm_bindgen]
pub struct SubstanceGroupKind {
    pub(crate) inner: ck::SubstanceGroupKind,
}
#[wasm_bindgen]
impl SubstanceGroupKind {
    #[wasm_bindgen(getter,js_name=Data)]
    pub fn variant_0() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Data,
        }
    }
    #[wasm_bindgen(getter,js_name=Superatom)]
    pub fn variant_1() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Superatom,
        }
    }
    #[wasm_bindgen(getter,js_name=MultipleGroup)]
    pub fn variant_2() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::MultipleGroup,
        }
    }
    #[wasm_bindgen(getter,js_name=StructuralRepeatUnit)]
    pub fn variant_3() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::StructuralRepeatUnit,
        }
    }
    #[wasm_bindgen(getter,js_name=Monomer)]
    pub fn variant_4() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Monomer,
        }
    }
    #[wasm_bindgen(getter,js_name=Copolymer)]
    pub fn variant_5() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Copolymer,
        }
    }
    #[wasm_bindgen(getter,js_name=Crosslink)]
    pub fn variant_6() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Crosslink,
        }
    }
    #[wasm_bindgen(getter,js_name=Graft)]
    pub fn variant_7() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Graft,
        }
    }
    #[wasm_bindgen(getter,js_name=Modification)]
    pub fn variant_8() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Modification,
        }
    }
    #[wasm_bindgen(getter,js_name=Mer)]
    pub fn variant_9() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Mer,
        }
    }
    #[wasm_bindgen(getter,js_name=AnyPolymer)]
    pub fn variant_10() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::AnyPolymer,
        }
    }
    #[wasm_bindgen(getter,js_name=MixtureComponent)]
    pub fn variant_11() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::MixtureComponent,
        }
    }
    #[wasm_bindgen(getter,js_name=Mixture)]
    pub fn variant_12() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Mixture,
        }
    }
    #[wasm_bindgen(getter,js_name=Formulation)]
    pub fn variant_13() -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Formulation,
        }
    }
    #[wasm_bindgen(js_name=Generic)]
    pub fn generic(value: String) -> Self {
        Self {
            inner: ck::SubstanceGroupKind::Generic(value.into()),
        }
    }
    #[wasm_bindgen(getter,js_name=genericValue,unchecked_return_type="string | null")]
    pub fn generic_value(&self) -> Result<JsValue, JsValue> {
        match &self.inner {
            ck::SubstanceGroupKind::Generic(v) => crate::host_values::text(v).map(JsValue::from),
            _ => Ok(JsValue::NULL),
        }
    }
    pub fn equals(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}
#[wasm_bindgen]
pub struct SGroupBracket {
    pub(crate) inner: ck::SGroupBracket,
}
#[wasm_bindgen]
impl SGroupBracket {
    pub fn new(
        #[wasm_bindgen(unchecked_param_type = "[number, number, number][]")] points: JsValue,
    ) -> Result<Self, JsValue> {
        let a = sequence(&points, "points")?;
        if a.length() != 3 {
            return Err(type_error("three points"));
        }
        Ok(Self {
            inner: ck::SGroupBracket::new([
                point(&a.get(0))?,
                point(&a.get(1))?,
                point(&a.get(2))?,
            ]),
        })
    }
    #[wasm_bindgen(unchecked_return_type = "[number, number, number][]")]
    pub fn points(&self) -> JsValue {
        self.inner
            .points()
            .iter()
            .map(vector)
            .collect::<js_sys::Array>()
            .into()
    }
}
#[wasm_bindgen]
pub struct SGroupCState {
    pub(crate) inner: ck::SGroupCState,
}
#[wasm_bindgen]
impl SGroupCState {
    pub fn new(
        #[wasm_bindgen(unchecked_param_type = "number")] bond: JsValue,
        #[wasm_bindgen(unchecked_param_type = "[number, number, number]")] value: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::SGroupCState::new(
                ck::BondId::new(usize_value(&bond, "bond")?),
                point(&value)?,
            ),
        })
    }
    pub fn bond(&self) -> usize {
        self.inner.bond().index()
    }
    #[wasm_bindgen(unchecked_return_type = "[number, number, number]")]
    pub fn vector(&self) -> JsValue {
        vector(self.inner.vector())
    }
}
#[wasm_bindgen]
pub struct SGroupDisplay {
    pub(crate) inner: ck::SGroupDisplay,
}
#[wasm_bindgen]
impl SGroupDisplay {
    #[wasm_bindgen(unchecked_return_type = "SGroupBracket[]")]
    pub fn brackets(&self) -> JsValue {
        self.inner
            .brackets()
            .iter()
            .copied()
            .map(|inner| JsValue::from(SGroupBracket { inner }))
            .collect::<js_sys::Array>()
            .into()
    }
}
#[wasm_bindgen]
pub struct SubstanceGroup {
    pub(crate) inner: ck::SubstanceGroup,
}
#[wasm_bindgen(
    inline_js = "export function visitGroupId(v,f){try{f(v);}catch(cause){throw new TypeError('invalid SubstanceGroupId',{cause});}} export function visitGroupKind(v,f){try{f(v);}catch(cause){throw new TypeError('invalid SubstanceGroupKind',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitGroupId)]
    fn visit_id(value: &JsValue, visit: &mut dyn FnMut(&SubstanceGroupId)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitGroupKind)]
    fn visit_kind(
        value: &JsValue,
        visit: &mut dyn FnMut(&SubstanceGroupKind),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen]
impl SubstanceGroup {
    pub fn new(
        #[wasm_bindgen(unchecked_param_type = "SubstanceGroupId")] id: JsValue,
        #[wasm_bindgen(unchecked_param_type = "SubstanceGroupKind")] kind: JsValue,
    ) -> Result<Self, JsValue> {
        let mut result = None;
        visit_id(&id, &mut |id: &SubstanceGroupId| {
            let mut value = None;
            let visited = visit_kind(&kind, &mut |kind: &SubstanceGroupKind| {
                value = Some(Self {
                    inner: ck::SubstanceGroup::new(id.inner, kind.inner.clone()),
                });
            });
            result = Some(visited.and_then(|()| value.ok_or_else(|| type_error("kind"))));
        })?;
        result.ok_or_else(|| type_error("id"))?
    }

    pub fn id(&self) -> SubstanceGroupId {
        SubstanceGroupId {
            inner: self.inner.id(),
        }
    }
    pub fn kind(&self) -> SubstanceGroupKind {
        SubstanceGroupKind {
            inner: self.inner.kind().clone(),
        }
    }
    #[wasm_bindgen(unchecked_return_type = "SGroupDisplay | null")]
    pub fn display(&self) -> JsValue {
        self.inner.display().map_or(JsValue::NULL, |inner| {
            SGroupDisplay {
                inner: inner.clone(),
            }
            .into()
        })
    }
    #[wasm_bindgen(unchecked_return_type = "SGroupCState[]")]
    pub fn cstates(&self) -> JsValue {
        self.inner
            .cstates()
            .iter()
            .copied()
            .map(|inner| JsValue::from(SGroupCState { inner }))
            .collect::<js_sys::Array>()
            .into()
    }
    #[wasm_bindgen(js_name=headCrossingBonds,unchecked_return_type="number[]")]
    pub fn head_crossing_bonds(&self) -> JsValue {
        self.inner
            .head_crossing_bonds()
            .iter()
            .map(|id| JsValue::from(id.index() as u32))
            .collect::<js_sys::Array>()
            .into()
    }
    #[wasm_bindgen(js_name=crossingBondCorrespondence,unchecked_return_type="number[]")]
    pub fn crossing_bond_correspondence(&self) -> JsValue {
        self.inner
            .crossing_bond_correspondence()
            .iter()
            .map(|id| JsValue::from(id.index() as u32))
            .collect::<js_sys::Array>()
            .into()
    }
}
