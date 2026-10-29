//! S-group values retain their native identities and explicit coordinate shapes.
use crate::{fields::Field, input, property_errors::text};
use cosmolkit as ck;
use serde::{Deserialize, Serialize};
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
#[derive(Serialize, Deserialize, Tsify)]
#[serde(transparent)]
pub struct Point3(#[serde(deserialize_with = "input::exact_array")] pub(crate) [f64; 3]);
#[derive(Serialize, Deserialize, Tsify)]
#[serde(transparent)]
pub struct Point2(#[serde(deserialize_with = "input::exact_array")] pub(crate) [f64; 2]);
#[derive(Serialize, Deserialize, Tsify)]
#[serde(transparent)]
pub struct BracketCoordinates(#[serde(deserialize_with = "input::exact_array")] [Point3; 3]);
#[wasm_bindgen]
pub struct SubstanceGroupId(pub(crate) ck::SubstanceGroupId);
crate::fields::value_field!(
    ck::SubstanceGroupId,
    SubstanceGroupId,
    SubstanceGroupIdInput,
    "SubstanceGroupId"
);
#[wasm_bindgen]
impl SubstanceGroupId {
    #[wasm_bindgen(constructor)]
    pub fn new(index: Ts<input::IndexInput>) -> Result<Self, JsValue> {
        Ok(Self(ck::SubstanceGroupId::new(input::read(index)?)))
    }
    #[wasm_bindgen(getter)]
    pub fn index(&self) -> usize {
        self.0.index()
    }
}
#[wasm_bindgen]
pub struct SubstanceGroupKind(pub(crate) ck::SubstanceGroupKind);
crate::fields::value_field!(
    ck::SubstanceGroupKind,
    SubstanceGroupKind,
    SubstanceGroupKindInput,
    "SubstanceGroupKind"
);
macro_rules! group_kinds{($($name:ident),+)=>{#[wasm_bindgen]impl SubstanceGroupKind{
    $(#[wasm_bindgen(getter,js_name=$name)]#[allow(non_snake_case)]pub fn $name()->Self{Self(ck::SubstanceGroupKind::$name)})+
    #[wasm_bindgen(js_name=Generic)]pub fn generic(value:Ts<input::TextInput>)->Result<Self,JsValue>{Ok(Self(ck::SubstanceGroupKind::Generic(input::read::<String>(value)?.into())))}
    #[wasm_bindgen(getter,js_name=genericValue)]pub fn generic_value(&self)->Result<Option<String>,JsValue>{match &self.0{ck::SubstanceGroupKind::Generic(v)=>text(v).map(Some),$(ck::SubstanceGroupKind::$name=>Ok(None)),+}}
}};}
group_kinds!(
    Data,
    Superatom,
    MultipleGroup,
    StructuralRepeatUnit,
    Monomer,
    Copolymer,
    Crosslink,
    Graft,
    Modification,
    Mer,
    AnyPolymer,
    MixtureComponent,
    Mixture,
    Formulation
);
#[wasm_bindgen]
pub struct SGroupBracket(pub(crate) ck::SGroupBracket);
crate::fields::value_field!(
    ck::SGroupBracket,
    SGroupBracket,
    SGroupBracketInput,
    "SGroupBracket"
);
crate::fields::value_array!(
    ck::SGroupBracket,
    SGroupBracketArray,
    "ReadonlyArray<SGroupBracket>"
);
#[wasm_bindgen]
impl SGroupBracket {
    #[wasm_bindgen(constructor)]
    pub fn new(points: Ts<BracketCoordinates>) -> Result<Self, JsValue> {
        Ok(Self(ck::SGroupBracket::new(
            input::decode(points)?.0.map(|p| p.0),
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn points(&self) -> Result<Ts<BracketCoordinates>, JsValue> {
        BracketCoordinates(self.0.points().map(Point3))
            .into_ts()
            .map_err(wasm_bindgen::JsError::from)
            .map_err(Into::into)
    }
}
#[wasm_bindgen]
pub struct SGroupCState(pub(crate) ck::SGroupCState);
crate::fields::value_field!(
    ck::SGroupCState,
    SGroupCState,
    SGroupCStateInput,
    "SGroupCState"
);
crate::fields::value_array!(
    ck::SGroupCState,
    SGroupCStateArray,
    "ReadonlyArray<SGroupCState>"
);
#[wasm_bindgen]
impl SGroupCState {
    #[wasm_bindgen(constructor)]
    pub fn new(bond: Ts<input::IndexInput>, vector: Ts<Point3>) -> Result<Self, JsValue> {
        Ok(Self(ck::SGroupCState::new(
            ck::BondId::new(input::read(bond)?),
            input::decode(vector)?.0,
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn bond(&self) -> usize {
        self.0.bond().index()
    }
    #[wasm_bindgen(getter)]
    pub fn vector(&self) -> Result<Ts<Point3>, JsValue> {
        Point3(*self.0.vector())
            .into_ts()
            .map_err(wasm_bindgen::JsError::from)
            .map_err(Into::into)
    }
}
#[wasm_bindgen]
pub struct SGroupDisplay(pub(crate) ck::SGroupDisplay);
crate::fields::value_field!(
    ck::SGroupDisplay,
    SGroupDisplay,
    SGroupDisplayInput,
    "SGroupDisplay"
);
#[wasm_bindgen]
impl SGroupDisplay {
    #[wasm_bindgen(constructor)]
    pub fn new(
        brackets: Option<SGroupBracketArray>,
        field_position: Option<Ts<Point2>>,
        display_tag: Option<Ts<input::TextInput>>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::SGroupDisplay {
            brackets: brackets
                .map(Vec::<ck::SGroupBracket>::input)
                .transpose()?
                .unwrap_or_else(|| ck::SGroupDisplay::default().brackets),
            field_position: field_position.map(input::decode).transpose()?.map(|p| p.0),
            display_tag: display_tag
                .map(input::read::<String>)
                .transpose()?
                .map(Into::into),
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn brackets(&self) -> Result<SGroupBracketArray, JsValue> {
        Vec::<ck::SGroupBracket>::get(crate::storage::NativeField::owned(
            self.0.brackets().to_vec(),
        ))
    }
    #[wasm_bindgen(getter,js_name=fieldPosition)]
    pub fn field_position(&self) -> Result<Option<Ts<Point2>>, JsValue> {
        self.0
            .field_position
            .map(|p| {
                Point2(p)
                    .into_ts()
                    .map_err(wasm_bindgen::JsError::from)
                    .map_err(Into::into)
            })
            .transpose()
    }
    #[wasm_bindgen(getter,js_name=displayTag)]
    pub fn display_tag(&self) -> Result<Option<String>, JsValue> {
        self.0.display_tag.as_ref().map(text).transpose()
    }
}
#[wasm_bindgen]
pub struct SubstanceGroup(pub(crate) ck::SubstanceGroup);
crate::fields::value_field!(
    ck::SubstanceGroup,
    SubstanceGroup,
    SubstanceGroupInput,
    "SubstanceGroup"
);
crate::fields::value_array!(
    ck::SubstanceGroup,
    SubstanceGroupArray,
    "ReadonlyArray<SubstanceGroup>"
);
#[wasm_bindgen]
impl SubstanceGroup {
    #[wasm_bindgen(constructor)]
    pub fn new(id: SubstanceGroupIdInput, kind: SubstanceGroupKindInput) -> Result<Self, JsValue> {
        SubstanceGroupId::with_ref(&id, |id| {
            SubstanceGroupKind::with_ref(&kind, |kind| {
                Ok(Self(ck::SubstanceGroup::new(*id, kind.clone())))
            })
        })
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> SubstanceGroupId {
        SubstanceGroupId(self.0.id())
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> SubstanceGroupKind {
        SubstanceGroupKind(self.0.kind().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn display(&self) -> Option<SGroupDisplay> {
        self.0.display().cloned().map(SGroupDisplay)
    }
    #[wasm_bindgen(getter)]
    pub fn cstates(&self) -> Result<SGroupCStateArray, JsValue> {
        Vec::<ck::SGroupCState>::get(crate::storage::NativeField::owned(
            self.0.cstates().to_vec(),
        ))
    }
    #[wasm_bindgen(getter,js_name=headCrossingBonds)]
    pub fn head_crossing_bonds(&self) -> Vec<usize> {
        self.0
            .head_crossing_bonds()
            .iter()
            .map(|v| v.index())
            .collect()
    }
    #[wasm_bindgen(getter,js_name=crossingBondCorrespondence)]
    pub fn crossing_bond_correspondence(&self) -> Vec<usize> {
        self.0
            .crossing_bond_correspondence()
            .iter()
            .map(|v| v.index())
            .collect()
    }
}
