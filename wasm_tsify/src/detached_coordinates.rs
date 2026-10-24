//! Detached coordinate blocks preserve dimensions and original conformer order.
use crate::Molecule;
use crate::{
    builder::{Points2, Points3},
    enums::{CoordinateDimension, CoordinateDimensionInput},
    fields::Field,
    group_values::{Point2, Point3},
    input,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(getter, js_name = conformers3d)]
    pub fn conformers_3d(&self) -> Result<Conformers3D, JsValue> {
        Vec::<ck::Conformer3D>::get(NativeField::owned(self.0.conformers_3d().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn properties(&self) -> crate::metadata_values::MoleculeProperties {
        crate::metadata_values::MoleculeProperties(self.0.properties().clone())
    }
    pub fn property(
        &self,
        key: Ts<input::TextInput>,
    ) -> Result<Option<crate::property_values::PropertyValue>, JsValue> {
        Ok(self
            .0
            .property(&input::read::<String>(key)?)
            .cloned()
            .map(crate::property_values::PropertyValue))
    }
}
#[wasm_bindgen]
pub struct Conformer2D(pub(crate) ck::Conformer2D);
crate::fields::value_field!(
    ck::Conformer2D,
    Conformer2D,
    Conformer2DInput,
    "Conformer2D"
);
crate::fields::value_array!(ck::Conformer2D, Conformers2D, "ReadonlyArray<Conformer2D>");
#[wasm_bindgen]
impl Conformer2D {
    #[wasm_bindgen(constructor)]
    pub fn new(id: Ts<input::IndexInput>, coordinates: Ts<Points2>) -> Result<Self, JsValue> {
        Ok(Self(ck::Conformer2D::new(
            input::read(id)?,
            input::decode(coordinates)?
                .0
                .into_iter()
                .map(|p| p.0)
                .collect(),
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> usize {
        self.0.id()
    }
    #[wasm_bindgen(getter)]
    pub fn coordinates(&self) -> Result<Ts<Points2>, JsValue> {
        Points2(self.0.coordinates().iter().copied().map(Point2).collect())
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
}
#[wasm_bindgen]
pub struct Conformer3D(pub(crate) ck::Conformer3D);
crate::fields::value_field!(
    ck::Conformer3D,
    Conformer3D,
    Conformer3DInput,
    "Conformer3D"
);
crate::fields::value_array!(ck::Conformer3D, Conformers3D, "ReadonlyArray<Conformer3D>");
#[wasm_bindgen]
impl Conformer3D {
    #[wasm_bindgen(constructor)]
    pub fn new(
        id: Ts<input::IndexInput>,
        coordinates: Ts<Points3>,
        is_3d: Ts<input::BoolInput>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::Conformer3D::new(
            input::read(id)?,
            input::decode(coordinates)?
                .0
                .into_iter()
                .map(|p| p.0)
                .collect(),
            input::read(is_3d)?,
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> usize {
        self.0.id()
    }
    #[wasm_bindgen(getter,js_name=is3d)]
    pub fn is_3d(&self) -> bool {
        self.0.is_3d()
    }
    #[wasm_bindgen(getter)]
    pub fn coordinates(&self) -> Result<Ts<Points3>, JsValue> {
        Points3(self.0.coordinates().iter().copied().map(Point3).collect())
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
}
crate::fields::value_array!(
    ck::CoordinateDimension,
    CoordinateDimensions,
    "ReadonlyArray<CoordinateDimension>"
);
#[wasm_bindgen]
pub struct CoordinateBlock(pub(crate) ck::CoordinateBlock);
crate::fields::value_field!(
    ck::CoordinateBlock,
    CoordinateBlock,
    CoordinateBlockInput,
    "CoordinateBlock"
);
#[wasm_bindgen]
impl CoordinateBlock {
    #[wasm_bindgen(constructor)]
    pub fn new(
        conformers_2d: Conformers2D,
        conformers_3d: Conformers3D,
        source_coordinate_dim: Option<CoordinateDimensionInput>,
        source_conformer_order: Option<CoordinateDimensions>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::CoordinateBlock {
            conformers_2d: Vec::<ck::Conformer2D>::input(conformers_2d)?,
            conformers_3d: Vec::<ck::Conformer3D>::input(conformers_3d)?,
            source_coordinate_dim: Option::<ck::CoordinateDimension>::input(source_coordinate_dim)?,
            source_conformer_order: Option::<Vec<ck::CoordinateDimension>>::input(
                source_conformer_order,
            )?,
        }))
    }
    #[wasm_bindgen(getter,js_name=conformers2d)]
    pub fn conformers_2d(&self) -> Result<Conformers2D, JsValue> {
        Vec::<ck::Conformer2D>::get(NativeField::owned(self.0.conformers_2d.clone()))
    }
    #[wasm_bindgen(getter,js_name=conformers3d)]
    pub fn conformers_3d(&self) -> Result<Conformers3D, JsValue> {
        Vec::<ck::Conformer3D>::get(NativeField::owned(self.0.conformers_3d.clone()))
    }
    #[wasm_bindgen(getter,js_name=sourceCoordinateDim)]
    pub fn source_coordinate_dim(&self) -> Option<CoordinateDimension> {
        self.0.source_coordinate_dim.map(Into::into)
    }
    #[wasm_bindgen(getter,js_name=sourceConformerOrder)]
    pub fn source_conformer_order(&self) -> Result<Option<CoordinateDimensions>, JsValue> {
        Option::<Vec<ck::CoordinateDimension>>::get(NativeField::owned(
            self.0.source_conformer_order.clone(),
        ))
    }
}
