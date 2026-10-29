//! Lazy native enumeration and stateful random providers.
use crate::errors::{enum_error, error_projection};
use crate::storage::ReadonlyValue;
use crate::{Molecule, callback, fields::Field, input, storage::NativeField};
use cosmolkit as ck;

#[wasm_bindgen]
impl StereoisomerOptions {
    #[wasm_bindgen(js_name = setTryEmbedding)]
    pub fn set_try_embedding_method(&self, value: <bool as Field>::Input) -> Result<(), JsValue> {
        self.set_try_embedding(value)
    }
    #[wasm_bindgen(js_name = setOnlyUnassigned)]
    pub fn set_only_unassigned_method(&self, value: <bool as Field>::Input) -> Result<(), JsValue> {
        self.set_only_unassigned(value)
    }
    #[wasm_bindgen(js_name = setMaxIsomers)]
    pub fn set_max_isomers_method(&self, value: <usize as Field>::Input) -> Result<(), JsValue> {
        self.set_max_isomers(value)
    }
    #[wasm_bindgen(js_name = setRandomSource)]
    pub fn set_random_source_method(
        &self,
        value: <Option<ck::StereoisomerRandomSource> as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_random_source(value)
    }
    #[wasm_bindgen(js_name = setUnique)]
    pub fn set_unique_method(&self, value: <bool as Field>::Input) -> Result<(), JsValue> {
        self.set_unique(value)
    }
    #[wasm_bindgen(js_name = setOnlyStereoGroups)]
    pub fn set_only_stereo_groups_method(
        &self,
        value: <bool as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_only_stereo_groups(value)
    }
}

use std::error::Error;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "bigint")]
    pub type IntegerSeed;
    #[wasm_bindgen(typescript_type = "((bitCount: number) => bigint)")]
    pub type RandomBitsCallback;
    #[wasm_bindgen(typescript_type = "IteratorResult<Molecule, undefined>")]
    pub type StereoisomerIteration;
    #[wasm_bindgen(
        typescript_type = "EnumerationError | NonNullable<GenerationError['cause']> | SmilesParseError | PropertyValueError | PoisonError"
    )]
    #[derive(Clone)]
    pub type EnumerationCause;
}
fn integer_text(value: &JsValue) -> Result<String, JsValue> {
    let value = value
        .dyn_ref::<js_sys::BigInt>()
        .ok_or_else(|| js_sys::TypeError::new("expected bigint"))?;
    Ok(value.to_string(10)?.into())
}
fn random_provider<T: std::str::FromStr + 'static>(
    value: &JsValue,
) -> Result<Box<dyn FnMut(usize) -> Result<T, String> + Send + Sync>, JsValue> {
    let handle = callback::Handle::new(value)?;
    Ok(Box::new(move |count| {
        let result = handle
            .call(|| Ok(js_sys::Array::of1(&JsValue::from(count as u32))))
            .and_then(|v| integer_text(&v))
            .and_then(|text| {
                text.parse().map_err(|_| {
                    js_sys::TypeError::new("random bits must be a nonnegative bigint").into()
                })
            });
        result.map_err(|error| {
            callback::record(error);
            "JavaScript random provider failed".into()
        })
    }))
}
#[wasm_bindgen]
pub struct StereoisomerRandomSource(ck::StereoisomerRandomSource);
crate::fields::value_field!(
    ck::StereoisomerRandomSource,
    StereoisomerRandomSource,
    StereoisomerRandomSourceInput,
    "StereoisomerRandomSource"
);
#[wasm_bindgen]
impl StereoisomerRandomSource {
    #[wasm_bindgen(js_name=fromIntegerSeed)]
    pub fn from_integer_seed(value: IntegerSeed) -> Result<Self, JsValue> {
        let text = integer_text(value.as_ref())?;
        let seed = text
            .parse()
            .map_err(|_| js_sys::TypeError::new("invalid integer seed"))?;
        Ok(Self(ck::StereoisomerRandomSource::from_integer_seed(seed)))
    }
    #[wasm_bindgen(js_name=fromRandomBits)]
    pub fn from_random_bits(value: RandomBitsCallback) -> Result<Self, JsValue> {
        Ok(Self(ck::StereoisomerRandomSource::from_random_bits(
            random_provider(value.as_ref())?,
        )))
    }
    pub fn copy(&self) -> Self {
        Self(self.0.clone())
    }
}
ck::__ck_configuration_StereoisomerOptions!(
    crate::project_parameters,
    ck::StereoisomerOptions,
    {},
    { StereoisomerOptions }
);

#[wasm_bindgen]
pub struct StereoisomerIterator(NativeField<ck::StereoisomerIterator>);
#[wasm_bindgen]
impl StereoisomerOptions {
    pub fn default() -> Self {
        Self {
            inner: NativeField::owned(ck::StereoisomerOptions::default()),
        }
    }
}
#[wasm_bindgen]
impl StereoisomerIterator {
    #[wasm_bindgen(js_name = "[Symbol.iterator]")]
    pub fn iterator(&self) -> Self {
        Self(self.0.clone())
    }
    pub fn next(&self) -> Result<StereoisomerIteration, JsValue> {
        callback::run(|| {
            let next = self
                .0
                .write(Iterator::next)?
                .transpose()
                .map_err(crate::operation_errors::OperationError::from_rust)?;
            crate::iterator::result(next.map(|m| Molecule(m).into())).map(JsCast::unchecked_into)
        })
    }
    #[wasm_bindgen(getter,js_name=yieldedCount)]
    pub fn yielded_count(&self) -> Result<usize, JsValue> {
        self.0.read(ck::StereoisomerIterator::yielded_count)
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=stereoisomerCount)]
    pub fn stereoisomer_count(
        &self,
        options: Option<StereoisomerOptionsInput>,
    ) -> Result<js_sys::BigInt, JsValue> {
        callback::run(|| {
            StereoisomerOptions::with_input(options.as_ref(), |p| {
                let count = self
                    .0
                    .stereoisomer_count_with_options(p)
                    .map_err(EnumerationError::from_rust)?;
                js_sys::BigInt::new(&JsValue::from(count.to_string())).map_err(Into::into)
            })
        })
    }
    #[wasm_bindgen(js_name=stereoisomerCountWithOptions)]
    pub fn stereoisomer_count_with_options(
        &self,
        options: StereoisomerOptionsInput,
    ) -> Result<js_sys::BigInt, JsValue> {
        self.stereoisomer_count(Some(options))
    }
    #[wasm_bindgen(js_name=enumerateStereoisomers)]
    pub fn enumerate_stereoisomers(
        &self,
        options: Option<StereoisomerOptionsInput>,
    ) -> Result<StereoisomerIterator, JsValue> {
        callback::run(|| {
            StereoisomerOptions::with_input(options.as_ref(), |p| {
                self.0
                    .enumerate_stereoisomers_with_options(p)
                    .map(|v| StereoisomerIterator(NativeField::owned(v)))
                    .map_err(crate::operation_errors::OperationError::from_rust)
            })
        })
    }
    #[wasm_bindgen(js_name=enumerateStereoisomersWithOptions)]
    pub fn enumerate_stereoisomers_with_options(
        &self,
        options: StereoisomerOptionsInput,
    ) -> Result<StereoisomerIterator, JsValue> {
        self.enumerate_stereoisomers(Some(options))
    }
    #[wasm_bindgen(js_name=enumerateStereoisomersWithRandomBits)]
    pub fn enumerate_stereoisomers_with_random_bits(
        &self,
        options: StereoisomerOptionsInput,
        provider: RandomBitsCallback,
    ) -> Result<StereoisomerIterator, JsValue> {
        let provider = random_provider(provider.as_ref())?;
        callback::run(|| {
            StereoisomerOptions::with_input(Some(&options), |p| {
                self.0
                    .enumerate_stereoisomers_with_random_bits(p, provider)
                    .map(|v| StereoisomerIterator(NativeField::owned(v)))
                    .map_err(crate::operation_errors::OperationError::from_rust)
            })
        })
    }
}

error_projection!(PoisonError, std::sync::PoisonError<()>, {});
fn project(owner: ReadonlyValue<dyn Error>) -> EnumerationCause {
    if owner.is::<ck::EnumerationError>() {
        return EnumerationError::from_owner(
            owner.child(|e| e.downcast_ref().expect("checked enumeration error")),
        )
        .unchecked_into();
    }
    if owner.is::<ck::SmilesParseError>() {
        return crate::smiles_parse_errors::SmilesParseError::from_owner(
            owner.child(|e| e.downcast_ref().expect("checked SMILES error")),
        )
        .unchecked_into();
    }
    if owner.is::<std::sync::PoisonError<()>>() {
        return PoisonError::from_owner(
            owner.child(|e| e.downcast_ref().expect("checked poisoned provider")),
        )
        .unchecked_into();
    }
    crate::conformer_errors::project(owner).unchecked_into()
}
fn cause<T: Error + 'static>(owner: &ReadonlyValue<T>) -> Option<EnumerationCause> {
    owner
        .source()
        .map(|_| project(owner.child(|e| e.source().expect("retained source"))))
}
error_projection!(
    EnumerationRunError,
    ck::EnumerationRunError,
    {},
    source = EnumerationCause,
    cause = cause::<ck::EnumerationRunError>
);
enum_error!(EnumerationError,EnumerationErrorKind,ck::EnumerationError,{
    Valence=>ck::EnumerationError::Valence(_),Rings=>ck::EnumerationError::Rings(_),Stereo=>ck::EnumerationError::Stereo(_),DoubleBond=>ck::EnumerationError::DoubleBond(_),SmilesWrite=>ck::EnumerationError::SmilesWrite(_),PotentialStereo=>ck::EnumerationError::PotentialStereo(_),Hydrogens=>ck::EnumerationError::Hydrogens(_),Embedding=>ck::EnumerationError::Embedding(_),AtomProperty=>ck::EnumerationError::AtomProperty(_),BondProperty=>ck::EnumerationError::BondProperty(_),MoleculeProperty=>ck::EnumerationError::MoleculeProperty(_),Topology=>ck::EnumerationError::Topology(_),Coordinates=>ck::EnumerationError::Coordinates(_),InvalidFlipperAtom=>ck::EnumerationError::InvalidFlipperAtom{..},InvalidFlipperBond=>ck::EnumerationError::InvalidFlipperBond{..},RandomBitsSource=>ck::EnumerationError::RandomBitsSource{..},RandomSourcePoisoned=>ck::EnumerationError::RandomSourcePoisoned(_),EmbeddingCoordinates=>ck::EnumerationError::EmbeddingCoordinates
},{
    #[wasm_bindgen(getter)]pub fn atom(&self)->Option<usize>{match *self.inner{ck::EnumerationError::InvalidFlipperAtom{atom}=>Some(atom.index()),_=>None}}
    #[wasm_bindgen(getter)]pub fn bond(&self)->Option<usize>{match *self.inner{ck::EnumerationError::InvalidFlipperBond{bond}=>Some(bond.index()),_=>None}}
    #[wasm_bindgen(getter,js_name=bitCount)]pub fn bit_count(&self)->Option<usize>{match *self.inner{ck::EnumerationError::RandomBitsSource{bit_count,..}=>Some(bit_count),_=>None}}
    #[wasm_bindgen(getter)]pub fn message(&self)->Option<String>{match self.inner.as_ref(){ck::EnumerationError::RandomBitsSource{message,..}=>Some(message.clone()),_=>None}}
    #[wasm_bindgen(getter)]pub fn error(&self)->Option<EnumerationCause>{
        macro_rules! payload {($($variant:ident),+)=>{match self.inner.as_ref(){$(ck::EnumerationError::$variant(_)=>Some(project(self.inner.child(|e|match e{ck::EnumerationError::$variant(v)=>v as &dyn Error,_=>unreachable!("immutable variant")}))),)+_=>None}};}
        payload!(Valence,Rings,Stereo,DoubleBond,SmilesWrite,PotentialStereo,Hydrogens,Embedding,AtomProperty,BondProperty,MoleculeProperty,Topology,Coordinates,RandomSourcePoisoned)
    }
},source=EnumerationCause,cause=cause::<ck::EnumerationError>);
