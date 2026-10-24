//! MOL2 parsing and the closed native finalization source graph.
use crate::errors::{enum_error, error_projection, exception};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use std::error::Error;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | CoordinateValidationError | AtomPropertyError | MoleculePropertyError"
    )]
    #[derive(Clone)]
    pub type Mol2ReadCause;
    #[wasm_bindgen(typescript_type = "PropertyStringError")]
    pub type Mol2Property;
    #[wasm_bindgen(
        typescript_type = "ValenceError | StereoError | SanitizeError | MolPostError | HydrogenError | RingFindingError | CoordinateValidationError | LegacyStereoError | MoleculePropertyError"
    )]
    #[derive(Clone)]
    pub type Mol2PostCause;
}
fn read_source(e: &ck::Mol2ReadError) -> Option<Mol2ReadCause> {
    use ck::Mol2ReadError as E;
    Some(
        match e {
            E::Topology(v) => crate::model_errors::TopologyValidationError::from_rust(v.clone()),
            E::Coordinates(v) => {
                crate::model_errors::CoordinateValidationError::from_rust(v.clone())
            }
            E::AtomProperty(v) => crate::model_errors::AtomPropertyError::from_rust(v.clone()),
            E::MoleculeProperty(v) => {
                crate::model_errors::MoleculePropertyError::from_rust(v.clone())
            }
            E::Parse(_) | E::Unsupported { .. } | E::PropertyString(_) => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(Mol2ReadError,Mol2ReadErrorKind,ck::Mol2ReadError,{
    Parse=>ck::Mol2ReadError::Parse(_),Unsupported=>ck::Mol2ReadError::Unsupported{..},Topology=>ck::Mol2ReadError::Topology(_),Coordinates=>ck::Mol2ReadError::Coordinates(_),PropertyString=>ck::Mol2ReadError::PropertyString(_),AtomProperty=>ck::Mol2ReadError::AtomProperty(_),MoleculeProperty=>ck::Mol2ReadError::MoleculeProperty(_),
},{
    #[wasm_bindgen(getter)] pub fn parse(&self)->Option<String>{match self.inner.as_ref(){ck::Mol2ReadError::Parse(v)=>Some(v.clone()),_=>None}}
    #[wasm_bindgen(getter)] pub fn feature(&self)->Option<String>{match self.inner.as_ref(){ck::Mol2ReadError::Unsupported{feature}=>Some((*feature).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=propertyString)] pub fn property_string(&self)->Option<Mol2Property>{match self.inner.as_ref(){ck::Mol2ReadError::PropertyString(v)=>Some(crate::property_errors::PropertyStringError::from_rust(*v).unchecked_into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<Mol2ReadCause>{read_source(&self.inner)}
},source=Mol2ReadCause,cause=read_source);

fn post_source(owner: &ReadonlyValue<ck::Mol2PostError>) -> Option<Mol2PostCause> {
    let native = owner.source()?;
    // The private native stage constructor has nine source types. Preserve the
    // actual type, independently of stage strings or formatted messages.
    macro_rules! project {
        ($native:ty,$module:ident::$wrapper:ident) => {
            if native.is::<$native>() {
                return Some(
                    crate::$module::$wrapper::from_owner(owner.child(|e| {
                        e.source()
                            .expect("native stage source")
                            .downcast_ref::<$native>()
                            .expect("checked immutable source type")
                    }))
                    .unchecked_into(),
                );
            }
        };
    }
    project!(ck::ValenceError, valence_errors::ValenceError);
    project!(ck::StereoError, stereo_errors::StereoError);
    project!(ck::SanitizeError, sanitize_errors::SanitizeError);
    project!(ck::MolPostError, mol_post_errors::MolPostError);
    project!(ck::HydrogenError, hydrogen_errors::HydrogenError);
    project!(ck::RingFindingError, ring_errors::RingFindingError);
    project!(
        ck::CoordinateValidationError,
        model_errors::CoordinateValidationError
    );
    project!(
        ck::LegacyStereoError,
        legacy_stereo_errors::LegacyStereoError
    );
    Some(
        crate::model_errors::MoleculePropertyError::from_owner(owner.child(|e| {
            e.source()
                .expect("native stage source")
                .downcast_ref::<ck::MoleculePropertyError>()
                .expect("native StereoMetadata source")
        }))
        .unchecked_into(),
    )
}
error_projection!(
    Mol2PostError,
    ck::Mol2PostError,
    {
        #[wasm_bindgen(getter)]
        pub fn stage(&self) -> String {
            self.inner.stage.into()
        }
        #[wasm_bindgen(getter)]
        pub fn source(&self) -> Option<Mol2PostCause> {
            post_source(&self.inner)
        }
    },
    source = Mol2PostCause,
    cause = post_source
);
