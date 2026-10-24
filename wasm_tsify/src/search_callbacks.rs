//! Native matcher policy and typed, synchronous callback projections.
use crate::reaction::Reaction;
use crate::{
    Molecule, callback,
    fields::Field,
    query_matches::{MatchResult, MatchResults},
    query_values::{QueryGraph, QueryGraphInput},
    storage::NativeField,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use std::sync::Arc;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "((query: QueryAtom, target: Atom) => boolean)")]
    pub type AtomMatchCallback;
    #[wasm_bindgen(typescript_type = "((query: Bond, target: Bond) => boolean)")]
    pub type BondMatchCallback;
    #[wasm_bindgen(
        typescript_type = "((target: Molecule, atomIndices: ReadonlyArray<number>) => boolean)"
    )]
    pub type FinalMatchCallback;
    #[wasm_bindgen(typescript_type = "Readonly<SubstructMatchParams>")]
    pub type ReadonlySubstructMatchParams;
}
macro_rules! callback_field {
    ($native:ty,$input:ident,$create:expr) => {
        impl Field for $native {
            type Input = $input;
            type Output = $input;
            fn input(value: Self::Input) -> Result<Self, JsValue> {
                Self::decode(value.as_ref())
            }
            fn decode(value: &JsValue) -> Result<Self, JsValue> {
                let handle = callback::Handle::new(value)?;
                let callback_handle = Arc::clone(&handle);
                let native: Self = ($create)(callback_handle);
                handle.bind(&native);
                Ok(native)
            }
            fn get(field: NativeField<Self>) -> Result<Self::Output, JsValue> {
                field.read(callback::original)?.map(JsCast::unchecked_into)
            }
        }
    };
}
callback_field!(ck::ExtraAtomCheck, AtomMatchCallback, |handle: Arc<
    callback::Handle,
>| Arc::new(
    move |_: &ck::QueryGraph, query: &ck::QueryAtom, _: &ck::SearchTarget<'_>, atom: &ck::Atom| {
        handle.boolean(|| {
            let args = js_sys::Array::new();
            args.push(&crate::query_values::QueryAtom(ReadonlyValue::owned(query.clone())).into());
            args.push(&crate::atoms::Atom(atom.clone()).into());
            Ok(args)
        })
    }
));
callback_field!(ck::ExtraBondCheck, BondMatchCallback, |handle: Arc<
    callback::Handle,
>| Arc::new(
    move |query: &ck::Bond, target: &ck::Bond| {
        handle.boolean(|| {
            let args = js_sys::Array::new();
            args.push(&crate::atoms::Bond(query.clone()).into());
            args.push(&crate::atoms::Bond(target.clone()).into());
            Ok(args)
        })
    }
));
callback_field!(ck::ExtraFinalCheck, FinalMatchCallback, |handle: Arc<
    callback::Handle,
>| Arc::new(
    move |target: &ck::SearchTarget<'_>, indices: &[usize]| {
        handle.boolean(|| {
            let args = js_sys::Array::new();
            args.push(&Molecule(callback::target(target.input_index())?).into());
            let indices = Vec::<usize>::get(NativeField::owned(indices.to_vec()))?;
            args.push(&indices.js_value());
            Ok(args)
        })
    }
));
ck::__ck_configuration_SubstructMatchParams!(
    crate::project_parameters,
    ck::SubstructMatchParams,
    {},
    { SubstructMatchParams }
);

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=substructMatch)]
    pub fn substruct_match(
        &self,
        query: QueryGraphInput,
        params: Option<SubstructMatchParamsInput>,
    ) -> Result<Option<MatchResult>, JsValue> {
        QueryGraph::with_ref(&query, |q| {
            if params.is_none() {
                return self
                    .0
                    .substruct_match(q)
                    .map(|value| value.map(MatchResult))
                    .map_err(crate::match_errors::SubstructMatchError::from_rust);
            }
            SubstructMatchParams::with_input(params.as_ref(), |p| {
                callback::run_with_targets(
                    p.extra_final_check.as_ref().map(|_| vec![self.0.clone()]),
                    || {
                        self.0
                            .substruct_match_with_params(q, p)
                            .map(|r| r.map(MatchResult))
                            .map_err(crate::match_errors::SubstructMatchError::from_rust)
                    },
                )
            })
        })
    }
    #[wasm_bindgen(js_name=substructMatches)]
    pub fn substruct_matches(
        &self,
        query: QueryGraphInput,
        params: Option<SubstructMatchParamsInput>,
    ) -> Result<MatchResults, JsValue> {
        QueryGraph::with_ref(&query, |q| {
            SubstructMatchParams::with_input(params.as_ref(), |p| {
                callback::run_with_targets(
                    p.extra_final_check.as_ref().map(|_| vec![self.0.clone()]),
                    || {
                        let result = self
                            .0
                            .substruct_matches_with_params(q, p)
                            .map_err(crate::match_errors::SubstructMatchError::from_rust)?;
                        Vec::<ck::MatchResult>::get(NativeField::owned(result))
                    },
                )
            })
        })
    }
    #[wasm_bindgen(js_name=hasSubstructMatch)]
    pub fn has_substruct_match(
        &self,
        query: QueryGraphInput,
        params: Option<SubstructMatchParamsInput>,
    ) -> Result<bool, JsValue> {
        QueryGraph::with_ref(&query, |q| {
            if params.is_none() {
                return self
                    .0
                    .has_substruct_match(q)
                    .map_err(crate::match_errors::SubstructMatchError::from_rust);
            }
            SubstructMatchParams::with_input(params.as_ref(), |p| {
                callback::run_with_targets(
                    p.extra_final_check.as_ref().map(|_| vec![self.0.clone()]),
                    || {
                        self.0
                            .has_substruct_match_with_params(q, p)
                            .map_err(crate::match_errors::SubstructMatchError::from_rust)
                    },
                )
            })
        })
    }
    #[wasm_bindgen(js_name=substructMatchWithParams)]
    pub fn substruct_match_with_params(
        &self,
        query: QueryGraphInput,
        params: SubstructMatchParamsInput,
    ) -> Result<Option<MatchResult>, JsValue> {
        self.substruct_match(query, Some(params))
    }
    #[wasm_bindgen(js_name=substructMatchesWithParams)]
    pub fn substruct_matches_with_params(
        &self,
        query: QueryGraphInput,
        params: SubstructMatchParamsInput,
    ) -> Result<MatchResults, JsValue> {
        self.substruct_matches(query, Some(params))
    }
    #[wasm_bindgen(js_name=hasSubstructMatchWithParams)]
    pub fn has_substruct_match_with_params(
        &self,
        query: QueryGraphInput,
        params: SubstructMatchParamsInput,
    ) -> Result<bool, JsValue> {
        self.has_substruct_match(query, Some(params))
    }
}
#[wasm_bindgen]
impl Reaction {
    #[wasm_bindgen(js_name=withMatchParams)]
    pub fn with_match_params(&self, params: SubstructMatchParamsInput) -> Result<Self, JsValue> {
        SubstructMatchParams::with_input(Some(&params), |p| {
            self.0.read(|r| {
                r.with_match_params(p)
                    .map(Self::from_native)
                    .map_err(crate::reaction_errors::ReactionModelError::from_rust)
            })?
        })
    }
    #[wasm_bindgen(js_name=matchParams)]
    pub fn match_params(&self) -> ReadonlySubstructMatchParams {
        let view = SubstructMatchParams {
            inner: self
                .0
                .readonly_child(|reaction| Some(reaction.match_params())),
        };
        JsValue::from(view).unchecked_into()
    }

    #[wasm_bindgen(js_name=setMatchParams)]
    pub fn set_match_params(&self, params: SubstructMatchParamsInput) -> Result<(), JsValue> {
        SubstructMatchParams::with_input(Some(&params), |p| self.0.write(|r| r.set_match_params(p)))
    }
}
