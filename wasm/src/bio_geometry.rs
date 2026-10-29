//! Structural coordinate values; calculation belongs exclusively to BIO Rust.
use crate::{builder::Points3, fields::Field, group_values::Point3, input, storage::NativeField};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
ck::__ck_configuration_BioCrystalCell!(crate::project_parameters, ck::BioCrystalCell, {}, {
    BioCrystalCell
});
#[derive(serde::Serialize, serde::Deserialize, Tsify)]
#[serde(transparent)]
pub struct Matrix3(#[serde(deserialize_with = "input::exact_array")] pub(crate) [Point3; 3]);
#[wasm_bindgen]
pub struct BioTransform(pub(crate) ck::BioTransform);
crate::fields::value_field!(
    ck::BioTransform,
    BioTransform,
    BioTransformInput,
    "BioTransform"
);
crate::fields::value_array!(
    ck::BioTransform,
    BioTransforms,
    "ReadonlyArray<BioTransform>"
);
#[wasm_bindgen]
impl BioTransform {
    #[wasm_bindgen(constructor)]
    pub fn new(
        matrix: Option<Ts<Matrix3>>,
        translation: Option<Ts<Point3>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioTransform::default();
        Ok(Self(ck::BioTransform::new(
            matrix
                .map(input::decode)
                .transpose()?
                .map(|m| m.0.map(|p| p.0))
                .unwrap_or(*d.matrix()),
            translation
                .map(input::decode)
                .transpose()?
                .map(|p| p.0)
                .unwrap_or(*d.translation()),
        )))
    }
    pub fn approx(
        &self,
        other: BioTransformInput,
        epsilon: Ts<input::F64Input>,
    ) -> Result<bool, JsValue> {
        let epsilon = input::read(epsilon)?;
        Self::with_ref(&other, |other| Ok(self.0.approx(other, epsilon)))
    }
    #[wasm_bindgen(getter)]
    pub fn matrix(&self) -> Result<Ts<Matrix3>, JsValue> {
        Matrix3(self.0.matrix().map(Point3))
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter)]
    pub fn translation(&self) -> Result<Ts<Point3>, JsValue> {
        Point3(*self.0.translation())
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
}
#[wasm_bindgen]
pub struct BioCoordinateBlock(pub(crate) ck::BioCoordinateBlock);
crate::fields::value_field!(
    ck::BioCoordinateBlock,
    BioCoordinateBlock,
    BioCoordinateBlockInput,
    "BioCoordinateBlock"
);
#[wasm_bindgen]
impl BioCoordinateBlock {
    #[wasm_bindgen(constructor)]
    pub fn new(positions: Ts<Points3>) -> Result<Self, JsValue> {
        Ok(Self(ck::BioCoordinateBlock::new(
            input::decode(positions)?
                .0
                .into_iter()
                .map(|p| p.0)
                .collect(),
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn positions(&self) -> Result<Ts<Points3>, JsValue> {
        Points3(self.0.positions().iter().copied().map(Point3).collect())
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter)]
    pub fn len(&self) -> usize {
        self.0.len()
    }
    #[wasm_bindgen(getter,js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}
#[wasm_bindgen]
pub struct BioCrystalInfo(pub(crate) ck::BioCrystalInfo);
crate::fields::value_field!(
    ck::BioCrystalInfo,
    BioCrystalInfo,
    BioCrystalInfoInput,
    "BioCrystalInfo"
);
#[wasm_bindgen]
impl BioCrystalInfo {
    #[wasm_bindgen(constructor)]
    pub fn new(
        cell: BioCrystalCellInput,
        space_group_hm: Option<Ts<input::TextInput>>,
        z_pdb: Option<Ts<input::TextInput>>,
        orthogonal: BioTransformInput,
        fractional: BioTransformInput,
        explicit_matrices: Ts<input::BoolInput>,
        cs_count: Ts<input::I16Input>,
        symmetry_images: BioTransforms,
    ) -> Result<Self, JsValue> {
        BioCrystalCell::with_input(Some(&cell), |cell| {
            Ok(Self(ck::BioCrystalInfo::new(
                *cell,
                space_group_hm.map(input::read).transpose()?,
                z_pdb.map(input::read).transpose()?,
                ck::BioTransform::input(orthogonal)?,
                ck::BioTransform::input(fractional)?,
                input::read(explicit_matrices)?,
                input::read(cs_count)?,
                Vec::<ck::BioTransform>::input(symmetry_images)?,
            )))
        })
    }
    #[wasm_bindgen(getter,js_name=spaceGroupNumber)]
    pub fn space_group_number(&self) -> Option<i32> {
        self.0.space_group_number()
    }
    #[wasm_bindgen(getter,js_name=isCrystal)]
    pub fn is_crystal(&self) -> bool {
        self.0.is_crystal()
    }
    #[wasm_bindgen(getter)]
    pub fn cell(&self) -> Result<BioCrystalCell, JsValue> {
        ck::BioCrystalCell::get(NativeField::owned(self.0.cell()))
    }
    #[wasm_bindgen(getter)]
    pub fn orthogonal(&self) -> BioTransform {
        BioTransform(*self.0.orthogonal())
    }
    #[wasm_bindgen(getter)]
    pub fn fractional(&self) -> BioTransform {
        BioTransform(*self.0.fractional())
    }
    #[wasm_bindgen(getter)]
    pub fn volume(&self) -> f64 {
        self.0.volume()
    }
    #[wasm_bindgen(getter,js_name=reciprocalLengths)]
    pub fn reciprocal_lengths(&self) -> Result<Ts<Point3>, JsValue> {
        Point3(*self.0.reciprocal_lengths())
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=reciprocalCosines)]
    pub fn reciprocal_cosines(&self) -> Result<Ts<Point3>, JsValue> {
        Point3(*self.0.reciprocal_cosines())
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=explicitMatrices)]
    pub fn explicit_matrices(&self) -> bool {
        self.0.explicit_matrices()
    }
    #[wasm_bindgen(getter,js_name=csCount)]
    pub fn cs_count(&self) -> i16 {
        self.0.cs_count()
    }
    #[wasm_bindgen(getter,js_name=symmetryImages)]
    pub fn symmetry_images(&self) -> Result<BioTransforms, JsValue> {
        Vec::<ck::BioTransform>::get(NativeField::owned(self.0.symmetry_images().to_vec()))
    }
    #[wasm_bindgen(getter,js_name=spaceGroupHm)]
    pub fn space_group_hm(&self) -> Option<String> {
        self.0.space_group_hm().map(str::to_owned)
    }
    #[wasm_bindgen(getter,js_name=zPdb)]
    pub fn z_pdb(&self) -> Option<String> {
        self.0.z_pdb().map(str::to_owned)
    }
}
#[wasm_bindgen]
pub struct BioNcsOperator(pub(crate) ck::BioNcsOperator);
crate::fields::value_field!(
    ck::BioNcsOperator,
    BioNcsOperator,
    BioNcsOperatorInput,
    "BioNcsOperator"
);
crate::fields::value_array!(
    ck::BioNcsOperator,
    BioNcsOperators,
    "ReadonlyArray<BioNcsOperator>"
);
#[wasm_bindgen]
impl BioNcsOperator {
    #[wasm_bindgen(constructor)]
    pub fn new(
        id: Ts<input::TextInput>,
        given: Ts<input::BoolInput>,
        transform: BioTransformInput,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::BioNcsOperator::new(
            input::read(id)?,
            input::read(given)?,
            ck::BioTransform::input(transform)?,
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> String {
        self.0.id.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn given(&self) -> bool {
        self.0.given
    }
    #[wasm_bindgen(getter)]
    pub fn transform(&self) -> BioTransform {
        BioTransform(self.0.transform)
    }
}
#[wasm_bindgen]
pub struct BioAssemblyOperator(pub(crate) ck::BioAssemblyOperator);
crate::fields::value_field!(
    ck::BioAssemblyOperator,
    BioAssemblyOperator,
    BioAssemblyOperatorInput,
    "BioAssemblyOperator"
);
crate::fields::value_array!(
    ck::BioAssemblyOperator,
    BioAssemblyOperators,
    "ReadonlyArray<BioAssemblyOperator>"
);
#[wasm_bindgen]
impl BioAssemblyOperator {
    #[wasm_bindgen(constructor)]
    pub fn new(
        name: Option<Ts<input::TextInput>>,
        operator_type: Option<Ts<input::TextInput>>,
        transform: BioTransformInput,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::BioAssemblyOperator::new(
            name.map(input::read).transpose()?,
            operator_type.map(input::read).transpose()?,
            ck::BioTransform::input(transform)?,
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> Option<String> {
        self.0.name.clone()
    }
    #[wasm_bindgen(getter,js_name=operatorType)]
    pub fn operator_type(&self) -> Option<String> {
        self.0.operator_type.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn transform(&self) -> BioTransform {
        BioTransform(self.0.transform)
    }
}
