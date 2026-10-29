//! Explicit coordinate installation uses native COW operations and exact XYZ rows.
use crate::{
    Molecule, group_values::Point3, input, operation_errors::OperationError,
    transform_errors::Coordinate3DReadError,
};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
#[derive(serde::Serialize, serde::Deserialize, Tsify)]
#[serde(transparent)]
pub struct Coordinates3D(Vec<Point3>);
#[derive(serde::Deserialize, Tsify)]
#[serde(transparent)]
pub struct CoordinateRows(Vec<Vec<f64>>);
ck::__ck_configuration_AtomPositionParams!(
    crate::project_parameters,
    ck::AtomPositionParams,
    {},
    { AtomPositionParams }
);
ck::__ck_configuration_Coordinate3DInputParams!(
    crate::project_parameters,
    ck::Coordinate3DInputParams,
    {},
    { Coordinate3DInputParams }
);
ck::__ck_configuration_Coordinate2DInputParams!(
    crate::project_parameters,
    ck::Coordinate2DInputParams,
    {},
    { Coordinate2DInputParams }
);
ck::__ck_configuration_Replace3DCoordinatesParams!(
    crate::project_parameters,
    ck::Replace3DCoordinatesParams,
    {},
    { Replace3DCoordinatesParams }
);
fn rows(v: Ts<Coordinates3D>) -> Result<Vec<[f64; 3]>, JsValue> {
    Ok(input::decode(v)?.0.into_iter().map(|p| p.0).collect())
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=withAtomPosition)]
    pub fn with_atom_position(
        &self,
        atom: Ts<input::IndexInput>,
        position: Ts<Point3>,
        params: Option<AtomPositionParamsInput>,
    ) -> Result<Self, JsValue> {
        let atom = ck::AtomId::new(input::read(atom)?);
        let point = input::decode(position)?.0;
        AtomPositionParams::with_input(params.as_ref(), |p| {
            self.0
                .with_atom_position_with_params(atom, point, p)
                .map(Self)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=setAtomPosition_)]
    pub fn set_atom_position(
        &mut self,
        atom: Ts<input::IndexInput>,
        position: Ts<Point3>,
        params: Option<AtomPositionParamsInput>,
    ) -> Result<(), JsValue> {
        let atom = ck::AtomId::new(input::read(atom)?);
        let point = input::decode(position)?.0;
        AtomPositionParams::with_input(params.as_ref(), |p| {
            self.0
                .set_atom_position_with_params_(atom, point, p)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=coordinates3d)]
    pub fn coordinates_3d(
        &self,
        id: Option<Ts<input::IndexInput>>,
    ) -> Result<Ts<Coordinates3D>, JsValue> {
        let id = id.map(input::read::<usize>).transpose()?.unwrap_or(0);
        Coordinates3D(
            self.0
                .coordinates_3d(id)
                .map_err(Coordinate3DReadError::from_rust)?
                .iter()
                .copied()
                .map(Point3)
                .collect(),
        )
        .into_ts()
        .map_err(wasm_bindgen::JsError::from)
        .map_err(Into::into)
    }
    #[wasm_bindgen(js_name=withOnly3dConformer)]
    pub fn with_only_3d_conformer(
        &self,
        coordinates: Ts<Coordinates3D>,
        params: Option<Coordinate3DInputParamsInput>,
    ) -> Result<Self, JsValue> {
        let v = rows(coordinates)?;
        Coordinate3DInputParams::with_input(params.as_ref(), |p| {
            self.0
                .with_only_3d_conformer_with_params(v, p)
                .map(Self)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=setOnly3dConformer_)]
    pub fn set_only_3d_conformer(
        &mut self,
        coordinates: Ts<Coordinates3D>,
        params: Option<Coordinate3DInputParamsInput>,
    ) -> Result<usize, JsValue> {
        let v = rows(coordinates)?;
        Coordinate3DInputParams::with_input(params.as_ref(), |p| {
            self.0
                .set_only_3d_conformer_with_params_(v, p)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withAdded3dConformer)]
    pub fn with_added_3d_conformer(
        &self,
        coordinates: Ts<Coordinates3D>,
        params: Option<Coordinate3DInputParamsInput>,
    ) -> Result<Self, JsValue> {
        let v = rows(coordinates)?;
        Coordinate3DInputParams::with_input(params.as_ref(), |p| {
            self.0
                .with_added_3d_conformer_with_params(v, p)
                .map(Self)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=add3dConformer_)]
    pub fn add_3d_conformer(
        &mut self,
        coordinates: Ts<Coordinates3D>,
        params: Option<Coordinate3DInputParamsInput>,
    ) -> Result<usize, JsValue> {
        let v = rows(coordinates)?;
        Coordinate3DInputParams::with_input(params.as_ref(), |p| {
            self.0
                .add_3d_conformer_with_params_(v, p)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=with3dCoordinates)]
    pub fn with_3d_coordinates(
        &self,
        coordinates: Ts<Coordinates3D>,
        params: Option<Replace3DCoordinatesParamsInput>,
    ) -> Result<Self, JsValue> {
        let v = rows(coordinates)?;
        Replace3DCoordinatesParams::with_input(params.as_ref(), |p| {
            self.0
                .with_3d_coordinates_with_params(v, p)
                .map(Self)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=set3dCoordinates_)]
    pub fn set_3d_coordinates(
        &mut self,
        coordinates: Ts<Coordinates3D>,
        params: Option<Replace3DCoordinatesParamsInput>,
    ) -> Result<(), JsValue> {
        let v = rows(coordinates)?;
        Replace3DCoordinatesParams::with_input(params.as_ref(), |p| {
            self.0
                .set_3d_coordinates_with_params_(v, p)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=with2dCoordinateBlock)]
    pub fn with_2d_coordinate_block(
        &self,
        coordinates: Ts<CoordinateRows>,
        params: Option<Coordinate2DInputParamsInput>,
    ) -> Result<Self, JsValue> {
        let v = input::decode(coordinates)?.0;
        Coordinate2DInputParams::with_input(params.as_ref(), |p| {
            self.0
                .with_2d_coordinate_block_with_params(v, p)
                .map(Self)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=set2dCoordinates_)]
    pub fn set_2d_coordinates(
        &mut self,
        coordinates: Ts<CoordinateRows>,
        params: Option<Coordinate2DInputParamsInput>,
    ) -> Result<(), JsValue> {
        let v = input::decode(coordinates)?.0;
        Coordinate2DInputParams::with_input(params.as_ref(), |p| {
            self.0
                .set_2d_coordinates_with_params_(v, p)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withCleared3dConformers)]
    pub fn with_cleared_3d_conformers(&self) -> Result<Self, JsValue> {
        self.0
            .with_cleared_3d_conformers()
            .map(Self)
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=clear3dConformers_)]
    pub fn clear_3d_conformers(&mut self) -> Result<(), JsValue> {
        self.0
            .clear_3d_conformers_()
            .map_err(OperationError::from_rust)
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = withAtomPositionWithParams)]
    pub fn with_atom_position_with_params(
        &self,
        atom: Ts<input::IndexInput>,
        position: Ts<Point3>,
        params: AtomPositionParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_atom_position(atom, position, Some(params))
    }
    #[wasm_bindgen(js_name = with2dCoordinateBlockWithParams)]
    pub fn with_2d_coordinate_block_with_params(
        &self,
        coordinates: Ts<CoordinateRows>,
        params: Coordinate2DInputParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_2d_coordinate_block(coordinates, Some(params))
    }
    #[wasm_bindgen(js_name = with3dCoordinatesWithParams)]
    pub fn with_3d_coordinates_with_params(
        &self,
        coordinates: Ts<Coordinates3D>,
        params: Replace3DCoordinatesParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_3d_coordinates(coordinates, Some(params))
    }
    #[wasm_bindgen(js_name = withAdded3dConformerWithParams)]
    pub fn with_added_3d_conformer_with_params(
        &self,
        coordinates: Ts<Coordinates3D>,
        params: Coordinate3DInputParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_added_3d_conformer(coordinates, Some(params))
    }
    #[wasm_bindgen(js_name = withOnly3dConformerWithParams)]
    pub fn with_only_3d_conformer_with_params(
        &self,
        coordinates: Ts<Coordinates3D>,
        params: Coordinate3DInputParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_only_3d_conformer(coordinates, Some(params))
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=setAtomPositionWithParams_)]
    pub fn set_atom_position_with_params_(
        &mut self,
        atom: Ts<input::IndexInput>,
        position: Ts<Point3>,
        params: AtomPositionParamsInput,
    ) -> Result<(), JsValue> {
        self.set_atom_position(atom, position, Some(params))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=setOnly3dConformerWithParams_)]
    pub fn set_only_3d_conformer_with_params_(
        &mut self,
        coordinates: Ts<Coordinates3D>,
        params: Coordinate3DInputParamsInput,
    ) -> Result<usize, JsValue> {
        self.set_only_3d_conformer(coordinates, Some(params))
    }
    #[wasm_bindgen(js_name=add3dConformerWithParams_)]
    pub fn add_3d_conformer_with_params_(
        &mut self,
        coordinates: Ts<Coordinates3D>,
        params: Coordinate3DInputParamsInput,
    ) -> Result<usize, JsValue> {
        self.add_3d_conformer(coordinates, Some(params))
    }
    #[wasm_bindgen(js_name=set3dCoordinatesWithParams_)]
    pub fn set_3d_coordinates_with_params_(
        &mut self,
        coordinates: Ts<Coordinates3D>,
        params: Replace3DCoordinatesParamsInput,
    ) -> Result<(), JsValue> {
        self.set_3d_coordinates(coordinates, Some(params))
    }
    #[wasm_bindgen(js_name=set2dCoordinatesWithParams_)]
    pub fn set_2d_coordinates_with_params_(
        &mut self,
        coordinates: Ts<CoordinateRows>,
        params: Coordinate2DInputParamsInput,
    ) -> Result<(), JsValue> {
        self.set_2d_coordinates(coordinates, Some(params))
    }
}
