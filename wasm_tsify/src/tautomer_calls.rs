//! Native enumeration policy, explicit callbacks and finalized molecule values.
use crate::{
    Molecule, callback, fields::Field, input, storage::NativeField,
    tautomer_callbacks::TautomerEnumerationStatus,
};
use cosmolkit as ck;

#[wasm_bindgen]
impl TautomerParams {
    #[wasm_bindgen(js_name = setMaxTautomers)]
    pub fn set_max_tautomers_method(&self, value: <u32 as Field>::Input) -> Result<(), JsValue> {
        self.set_max_tautomers(value)
    }
    #[wasm_bindgen(js_name = setMaxTransforms)]
    pub fn set_max_transforms_method(&self, value: <u32 as Field>::Input) -> Result<(), JsValue> {
        self.set_max_transforms(value)
    }
    #[wasm_bindgen(js_name = setRemoveSp3Stereo)]
    pub fn set_remove_sp3_stereo_method(
        &self,
        value: <bool as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_remove_sp3_stereo(value)
    }
    #[wasm_bindgen(js_name = setRemoveBondStereo)]
    pub fn set_remove_bond_stereo_method(
        &self,
        value: <bool as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_remove_bond_stereo(value)
    }
    #[wasm_bindgen(js_name = setRemoveIsotopicHydrogens)]
    pub fn set_remove_isotopic_hydrogens_method(
        &self,
        value: <bool as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_remove_isotopic_hydrogens(value)
    }
    #[wasm_bindgen(js_name = setReassignStereo)]
    pub fn set_reassign_stereo_method(&self, value: <bool as Field>::Input) -> Result<(), JsValue> {
        self.set_reassign_stereo(value)
    }
}

use std::rc::Rc;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};

macro_rules! project_tautomer {
    (@resume {} @adapter $native:path {
        $(#[$policy_attr:meta])* policy:$policy:ty,
        $(#[$callback_attr:meta])* callback:$callback:ty,
        $(#[$scorer_attr:meta])* scorer:$scorer:ty,
        $(#[$score_attr:meta])* score_params:$score:ty
    } methods {}) => {
        cosmolkit_tsify_macros::parameters!(
            TautomerParams,
            $native {
                #[configuration(flatten)]
                #[projection_accessors]
                policy: $policy,
                #[projection_read(callback().cloned())]
                #[projection_write(set_callback)]
                callback: $callback,
                #[projection_read(scorer().cloned())]
                #[projection_write(set_scorer)]
                scorer: $scorer,
                score_params: $score
            }
        );
        #[wasm_bindgen]
        impl TautomerParams {
            #[wasm_bindgen(js_name = setCallback)]
            pub fn set_callback_method(
                &self,
                value: <$callback as Field>::Input,
            ) -> Result<(), JsValue> {
                self.set_callback(value)
            }
            #[wasm_bindgen(js_name = setScorer)]
            pub fn set_scorer_method(
                &self,
                value: <$scorer as Field>::Input,
            ) -> Result<(), JsValue> {
                self.set_scorer(value)
            }
        }
    };
}
ck::__ck_configuration_TautomerParams!(project_tautomer, ck::TautomerParams, {}, {});

#[wasm_bindgen(js_name = canonicalTautomerFromMolecules)]
pub fn canonical_tautomer_from_molecules(
    molecules: crate::fragments::Molecules,
) -> Result<Molecule, JsValue> {
    let molecules = Vec::<ck::Molecule>::input(molecules)?;
    callback::run(|| {
        ck::canonical_tautomer_from_molecules(&molecules)
            .map(Molecule)
            .map_err(crate::operation_errors::OperationError::from_rust)
    })
}

#[wasm_bindgen(js_name = canonicalTautomerFromMoleculesWithParams)]
pub fn canonical_tautomer_from_molecules_with_params(
    molecules: crate::fragments::Molecules,
    params: TautomerParamsInput,
) -> Result<Molecule, JsValue> {
    let molecules = Vec::<ck::Molecule>::input(molecules)?;
    TautomerParams::with_input(Some(&params), |params| {
        callback::run(|| {
            ck::canonical_tautomer_from_molecules_with_params(&molecules, params)
                .map(Molecule)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    })
}
macro_rules! policy_copy {
    ($($method:ident,$js:ident,$ty:ty);+)=>{
        #[wasm_bindgen]
        impl TautomerParams {
            $(#[wasm_bindgen(js_name=$js)]
            pub fn $method(&self,value:Ts<<$ty as input::CheckedValue>::Wire>)->Result<Self,JsValue>{
                let value=input::read::<$ty>(value)?;
                self.inner.read(|p|Self{inner:NativeField::owned(p.clone().$method(value))})
            })+
        }
    };
}
policy_copy!(
    with_max_tautomers,withMaxTautomers,u32;
    with_max_transforms,withMaxTransforms,u32;
    with_remove_sp3_stereo,withRemoveSp3Stereo,bool;
    with_remove_bond_stereo,withRemoveBondStereo,bool;
    with_remove_isotopic_hydrogens,withRemoveIsotopicHydrogens,bool;
    with_reassign_stereo,withReassignStereo,bool
);
#[wasm_bindgen]
impl TautomerParams {
    #[wasm_bindgen(getter,js_name=transformCount)]
    pub fn transform_count(&self) -> Result<usize, JsValue> {
        self.inner.read(ck::TautomerParams::transform_count)
    }
}
#[wasm_bindgen]
pub struct TautomerEnumeration(pub(crate) Rc<ck::TautomerEnumeration>);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ReadonlyArray<string>")]
    pub type TautomerSmiles;
    #[wasm_bindgen(typescript_type = "IteratorResult<Molecule, undefined>")]
    pub type TautomerIteration;
    #[wasm_bindgen(typescript_type = "IteratorResult<readonly [string, Molecule], undefined>")]
    pub type TautomerEntryIteration;
}
#[wasm_bindgen]
impl TautomerEnumeration {
    pub fn iter(&self) -> TautomerIter {
        self.iterator()
    }
    pub fn entries(&self) -> TautomerEntries {
        TautomerEntries(NativeField::owned((self.0.clone(), 0)))
    }
    #[wasm_bindgen(js_name = "[Symbol.iterator]")]
    pub fn iterator(&self) -> TautomerIter {
        TautomerIter(NativeField::owned((self.0.clone(), 0)))
    }
    #[wasm_bindgen(getter)]
    pub fn len(&self) -> usize {
        self.0.len()
    }
    #[wasm_bindgen(js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    #[wasm_bindgen(getter)]
    pub fn status(&self) -> TautomerEnumerationStatus {
        self.0.status().into()
    }
    pub fn get(&self, index: Ts<input::IndexInput>) -> Result<Option<Molecule>, JsValue> {
        Ok(self
            .0
            .get(input::read::<usize>(index)?)
            .cloned()
            .map(Molecule))
    }
    #[wasm_bindgen(getter,js_name=modifiedAtoms)]
    pub fn modified_atoms(&self) -> Result<Ts<crate::fields::IndexArray>, JsValue> {
        Vec::<usize>::get(NativeField::owned(
            self.0
                .modified_atoms()
                .iter()
                .map(|id| id.index())
                .collect(),
        ))
    }
    #[wasm_bindgen(getter,js_name=modifiedBonds)]
    pub fn modified_bonds(&self) -> Result<Ts<crate::fields::IndexArray>, JsValue> {
        Vec::<usize>::get(NativeField::owned(
            self.0
                .modified_bonds()
                .iter()
                .map(|id| id.index())
                .collect(),
        ))
    }
    #[wasm_bindgen(js_name=canonicalSmiles)]
    pub fn canonical_smiles(&self) -> Result<TautomerSmiles, JsValue> {
        let values = self
            .0
            .canonical_smiles()
            .into_iter()
            .map(|s| crate::property_errors::text(s).map(JsValue::from))
            .collect::<Result<js_sys::Array, JsValue>>()?;
        Ok(js_sys::Object::freeze(&values).unchecked_into())
    }
    #[wasm_bindgen(js_name=canonicalTautomer)]
    pub fn canonical_tautomer(
        &self,
        params: Option<TautomerParamsInput>,
    ) -> Result<Molecule, JsValue> {
        TautomerParams::with_input(params.as_ref(), |p| {
            callback::run(|| {
                self.0
                    .canonical_tautomer_with_params(p)
                    .map(Molecule)
                    .map_err(crate::operation_errors::OperationError::from_rust)
            })
        })
    }
    #[wasm_bindgen(js_name=canonicalTautomerWithParams)]
    pub fn canonical_tautomer_with_params(
        &self,
        params: TautomerParamsInput,
    ) -> Result<Molecule, JsValue> {
        self.canonical_tautomer(Some(params))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=enumerateTautomers)]
    pub fn enumerate_tautomers(
        &self,
        params: Option<TautomerParamsInput>,
    ) -> Result<TautomerEnumeration, JsValue> {
        TautomerParams::with_input(params.as_ref(), |p| {
            callback::run(|| {
                self.0
                    .enumerate_tautomers_with_params(p)
                    .map(|value| TautomerEnumeration(Rc::new(value)))
                    .map_err(crate::operation_errors::OperationError::from_rust)
            })
        })
    }
    #[wasm_bindgen(js_name=enumerateTautomersWithParams)]
    pub fn enumerate_tautomers_with_params(
        &self,
        params: TautomerParamsInput,
    ) -> Result<TautomerEnumeration, JsValue> {
        self.enumerate_tautomers(Some(params))
    }
    #[wasm_bindgen(js_name=canonicalTautomer)]
    pub fn canonical_tautomer(&self, params: Option<TautomerParamsInput>) -> Result<Self, JsValue> {
        TautomerParams::with_input(params.as_ref(), |p| {
            callback::run(|| {
                self.0
                    .canonical_tautomer_with_params(p)
                    .map(Self)
                    .map_err(crate::operation_errors::OperationError::from_rust)
            })
        })
    }
    #[wasm_bindgen(js_name=canonicalTautomerWithParams)]
    pub fn canonical_tautomer_with_params(
        &self,
        params: TautomerParamsInput,
    ) -> Result<Self, JsValue> {
        self.canonical_tautomer(Some(params))
    }
}
#[wasm_bindgen]
pub struct TautomerIter(NativeField<(Rc<ck::TautomerEnumeration>, usize)>);
#[wasm_bindgen]
pub struct TautomerEntries(NativeField<(Rc<ck::TautomerEnumeration>, usize)>);
#[wasm_bindgen]
impl TautomerEntries {
    #[wasm_bindgen(js_name = "[Symbol.iterator]")]
    pub fn iterator(&self) -> Self {
        Self(self.0.clone())
    }
    pub fn next(&self) -> Result<TautomerEntryIteration, JsValue> {
        let value = self.0.write(|(owner, index)| {
            let Some((key, molecule)) = owner.entries().nth(*index) else {
                return Ok(None);
            };
            let key = crate::property_errors::text(key)?;
            let tuple = js_sys::Array::of2(
                &JsValue::from(key),
                &JsValue::from(Molecule(molecule.clone())),
            );
            *index += 1;
            Ok::<_, JsValue>(Some(js_sys::Object::freeze(&tuple).into()))
        })??;
        crate::iterator::result(value).map(JsCast::unchecked_into)
    }
}
#[wasm_bindgen]
impl TautomerIter {
    #[wasm_bindgen(js_name = "[Symbol.iterator]")]
    pub fn iterator(&self) -> Self {
        Self(self.0.clone())
    }
    pub fn next(&self) -> Result<TautomerIteration, JsValue> {
        let value = self.0.write(|(owner, index)| {
            let value = owner
                .get(*index)
                .cloned()
                .map(|m| JsValue::from(Molecule(m)));
            if value.is_some() {
                *index += 1;
            }
            value
        })?;
        Ok(crate::iterator::result(value)?.unchecked_into())
    }
}
