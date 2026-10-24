//! Exact native coordinate failures; no coordinate policy is implemented here.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "TopologyValidationError | CoordinateValidationError")]
    #[derive(Clone)]
    pub type TransformCause;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type CoordinatePayload;
}
fn transform_cause(e: &ck::TransformError) -> Option<TransformCause> {
    match e {
        ck::TransformError::InvalidTopology(v) => Some(
            crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::TransformError::InvalidCoordinates(v) => Some(
            crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(TransformError,TransformErrorKind,ck::TransformError,{
    InvalidTopology=>ck::TransformError::InvalidTopology(_),InvalidCoordinates=>ck::TransformError::InvalidCoordinates(_),AtomOutOfRange=>ck::TransformError::AtomOutOfRange{..},WeightCount=>ck::TransformError::WeightCount{..},NonFiniteWeight=>ck::TransformError::NonFiniteWeight{..},InvalidWeightSum=>ck::TransformError::InvalidWeightSum{..},NoSelectedAtoms=>ck::TransformError::NoSelectedAtoms,NonFinitePoint=>ck::TransformError::NonFinitePoint{..},NonFiniteTransform=>ck::TransformError::NonFiniteTransform{..},EigenDecompositionDidNotConverge=>ck::TransformError::EigenDecompositionDidNotConverge,No3dConformer=>ck::TransformError::No3dConformer,ConformerNotFound=>ck::TransformError::ConformerNotFound{..},AtomsNotBonded=>ck::TransformError::AtomsNotBonded{..},RingBondNotMovable=>ck::TransformError::RingBondNotMovable{..},BothAngleBondsInRing=>ck::TransformError::BothAngleBondsInRing{..},CoincidentCoordinates=>ck::TransformError::CoincidentCoordinates{..},UndefinedRotationAxis=>ck::TransformError::UndefinedRotationAxis,InvalidTargetValue=>ck::TransformError::InvalidTargetValue{..}
},{
    #[wasm_bindgen(getter)]pub fn topology(&self)->Option<TransformCause>{match self.inner.as_ref(){ck::TransformError::InvalidTopology(v)=>Some(crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn coordinates(&self)->Option<CoordinatePayload>{match self.inner.as_ref(){ck::TransformError::InvalidCoordinates(v)=>Some(crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn atom(&self)->Option<usize>{match *self.inner{ck::TransformError::AtomOutOfRange{atom,..}|ck::TransformError::NonFiniteWeight{atom,..}=>Some(atom.index()),_=>None}}
    #[wasm_bindgen(getter)]pub fn role(&self)->Option<String>{match *self.inner{ck::TransformError::AtomOutOfRange{role,..}|ck::TransformError::NonFinitePoint{role,..}=>Some(role.into()),_=>None}}
    #[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match *self.inner{ck::TransformError::AtomOutOfRange{atom_count,..}=>Some(atom_count),_=>None}}
    #[wasm_bindgen(getter)]pub fn actual(&self)->Option<usize>{match *self.inner{ck::TransformError::WeightCount{actual,..}=>Some(actual),_=>None}}
    #[wasm_bindgen(getter)]pub fn required(&self)->Option<usize>{match *self.inner{ck::TransformError::WeightCount{required,..}=>Some(required),_=>None}}
    #[wasm_bindgen(getter)]pub fn value(&self)->Option<f64>{match *self.inner{ck::TransformError::NonFiniteWeight{value,..}|ck::TransformError::InvalidTargetValue{value,..}=>Some(value),_=>None}}
    #[wasm_bindgen(getter)]pub fn sum(&self)->Option<f64>{match *self.inner{ck::TransformError::InvalidWeightSum{sum}=>Some(sum),_=>None}}
    #[wasm_bindgen(getter)]pub fn axis(&self)->Option<String>{match *self.inner{ck::TransformError::NonFinitePoint{axis,..}=>Some(axis.into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn row(&self)->Option<usize>{match *self.inner{ck::TransformError::NonFiniteTransform{row,..}=>Some(row),_=>None}}
    #[wasm_bindgen(getter)]pub fn column(&self)->Option<usize>{match *self.inner{ck::TransformError::NonFiniteTransform{column,..}=>Some(column),_=>None}}
    #[wasm_bindgen(getter,js_name=conformerId)]pub fn conformer_id(&self)->Option<usize>{match *self.inner{ck::TransformError::ConformerNotFound{conformer_id}=>Some(conformer_id),_=>None}}
    #[wasm_bindgen(getter)]pub fn first(&self)->Option<usize>{match *self.inner{ck::TransformError::AtomsNotBonded{first,..}|ck::TransformError::CoincidentCoordinates{first,..}=>Some(first.index()),_=>None}}
    #[wasm_bindgen(getter)]pub fn second(&self)->Option<usize>{match *self.inner{ck::TransformError::AtomsNotBonded{second,..}|ck::TransformError::CoincidentCoordinates{second,..}=>Some(second.index()),_=>None}}
    #[wasm_bindgen(getter,js_name=firstRole)]pub fn first_role(&self)->Option<String>{match *self.inner{ck::TransformError::AtomsNotBonded{first_role,..}|ck::TransformError::CoincidentCoordinates{first_role,..}=>Some(first_role.into()),_=>None}}
    #[wasm_bindgen(getter,js_name=secondRole)]pub fn second_role(&self)->Option<String>{match *self.inner{ck::TransformError::AtomsNotBonded{second_role,..}|ck::TransformError::CoincidentCoordinates{second_role,..}=>Some(second_role.into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn bond(&self)->Option<usize>{match *self.inner{ck::TransformError::RingBondNotMovable{bond}=>Some(bond.index()),_=>None}}
    #[wasm_bindgen(getter,js_name=firstBond)]pub fn first_bond(&self)->Option<usize>{match *self.inner{ck::TransformError::BothAngleBondsInRing{first_bond,..}=>Some(first_bond.index()),_=>None}}
    #[wasm_bindgen(getter,js_name=secondBond)]pub fn second_bond(&self)->Option<usize>{match *self.inner{ck::TransformError::BothAngleBondsInRing{second_bond,..}=>Some(second_bond.index()),_=>None}}
    #[wasm_bindgen(getter)]pub fn quantity(&self)->Option<String>{match *self.inner{ck::TransformError::InvalidTargetValue{quantity,..}=>Some(quantity.into()),_=>None}}
},source=TransformCause,cause=transform_cause);
enum_error!(Coordinate3DReadError,Coordinate3DReadErrorKind,ck::Coordinate3DReadError,{ConformerNotFound=>ck::Coordinate3DReadError::ConformerNotFound{..}},{
    #[wasm_bindgen(getter,js_name=conformerId)]pub fn conformer_id(&self)->usize{let ck::Coordinate3DReadError::ConformerNotFound{conformer_id,..}=*self.inner;conformer_id}
    #[wasm_bindgen(getter)]pub fn count(&self)->usize{let ck::Coordinate3DReadError::ConformerNotFound{count,..}=*self.inner;count}
});
enum_error!(CoordinateInputError,CoordinateInputErrorKind,ck::CoordinateInputError,{
    UnknownZPolicy=>ck::CoordinateInputError::UnknownZPolicy{..},RowCount=>ck::CoordinateInputError::RowCount{..},Shape=>ck::CoordinateInputError::Shape{..},NonFinite=>ck::CoordinateInputError::NonFinite{..},NonZeroZ=>ck::CoordinateInputError::NonZeroZ{..},ThreeColumnsForbidden=>ck::CoordinateInputError::ThreeColumnsForbidden,ConformerNotFound=>ck::CoordinateInputError::ConformerNotFound{..},ConformerIdOverflow=>ck::CoordinateInputError::ConformerIdOverflow{..},InvalidCoordinates=>ck::CoordinateInputError::InvalidCoordinates(_)
},{
    #[wasm_bindgen(getter)]pub fn value(&self)->Option<CoordinateInputValue>{match self.inner.as_ref(){ck::CoordinateInputError::UnknownZPolicy{value}=>Some(JsValue::from(value).unchecked_into()),ck::CoordinateInputError::NonFinite{value,..}=>Some(JsValue::from(*value).unchecked_into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn dimension(&self)->Option<String>{match *self.inner{ck::CoordinateInputError::RowCount{dimension,..}|ck::CoordinateInputError::Shape{dimension,..}|ck::CoordinateInputError::NonFinite{dimension,..}=>Some(dimension.into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn expected(&self)->Option<CoordinateInputValue>{match *self.inner{ck::CoordinateInputError::RowCount{expected,..}=>Some(JsValue::from(expected).unchecked_into()),ck::CoordinateInputError::Shape{expected,..}=>Some(JsValue::from(expected).unchecked_into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn actual(&self)->Option<usize>{match *self.inner{ck::CoordinateInputError::RowCount{actual,..}=>Some(actual),_=>None}}
    #[wasm_bindgen(getter)]pub fn row(&self)->Option<usize>{match *self.inner{ck::CoordinateInputError::Shape{row,..}|ck::CoordinateInputError::NonFinite{row,..}|ck::CoordinateInputError::NonZeroZ{row,..}=>Some(row),_=>None}}
    #[wasm_bindgen(getter)]pub fn columns(&self)->Option<usize>{match *self.inner{ck::CoordinateInputError::Shape{columns,..}=>Some(columns),_=>None}}
    #[wasm_bindgen(getter)]pub fn column(&self)->Option<usize>{match *self.inner{ck::CoordinateInputError::NonFinite{column,..}=>Some(column),_=>None}}
    #[wasm_bindgen(getter)]pub fn z(&self)->Option<f64>{match *self.inner{ck::CoordinateInputError::NonZeroZ{z,..}=>Some(z),_=>None}}
    #[wasm_bindgen(getter,js_name=conformerId)]pub fn conformer_id(&self)->Option<usize>{match *self.inner{ck::CoordinateInputError::ConformerNotFound{conformer_id,..}=>Some(conformer_id),_=>None}}
    #[wasm_bindgen(getter)]pub fn count(&self)->Option<usize>{match *self.inner{ck::CoordinateInputError::ConformerNotFound{count,..}=>Some(count),_=>None}}
    #[wasm_bindgen(getter,js_name=maxId)]pub fn max_id(&self)->Option<usize>{match *self.inner{ck::CoordinateInputError::ConformerIdOverflow{max_id}=>Some(max_id),_=>None}}
    #[wasm_bindgen(getter)]pub fn coordinates(&self)->Option<CoordinatePayload>{match self.inner.as_ref(){ck::CoordinateInputError::InvalidCoordinates(v)=>Some(crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "number | string")]
    pub type CoordinateInputValue;
}
