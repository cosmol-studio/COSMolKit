//! Scoped native loans. Only snapshot() acquires independent native ownership.
use crate::{callback, fields::Field, input, storage::NativeField, storage::ReadonlyValue};
use cosmolkit as ck;
use std::{
    cell::{Cell, RefCell},
    rc::Rc,
    sync::Arc,
};
use tsify::Ts;
use wasm_bindgen::{JsCast, closure::ScopedClosure, prelude::*};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "((source: TautomerMoleculeView, progress: TautomerProgress) => boolean)"
    )]
    pub type EnumerationCallback;
    #[wasm_bindgen(typescript_type = "((source: TautomerMoleculeView) => number)")]
    pub type ScorerCallback;
    #[wasm_bindgen(typescript_type = "TautomerMoleculeView")]
    pub type ViewValue;
    #[wasm_bindgen(typescript_type = "TautomerProgress")]
    pub type ProgressValue;
    #[wasm_bindgen(typescript_type = "readonly [string, TautomerMoleculeView]")]
    pub type TautomerEntry;
}
crate::enums::enum_value!(TautomerEnumerationStatus,"TautomerEnumerationStatus"{Completed="completed",MaxTautomersReached="max_tautomers_reached",MaxTransformsReached="max_transforms_reached",Canceled="canceled"});

#[derive(Clone)]
struct Lease {
    access: js_sys::Function,
    active: Rc<Cell<bool>>,
}
impl Lease {
    fn call(
        &self,
        kind: u32,
        operation: u32,
        index: usize,
        argument: usize,
    ) -> Result<JsValue, JsValue> {
        if !self.active.get() {
            return Err(js_sys::ReferenceError::new(
                "callback view expired; retain snapshot() instead",
            )
            .into());
        }
        self.access.call4(
            &JsValue::UNDEFINED,
            &kind.into(),
            &operation.into(),
            &(index as u32).into(),
            &(argument as u32).into(),
        )
    }
}
struct Expire(Rc<Cell<bool>>);
impl Drop for Expire {
    fn drop(&mut self) {
        self.0.set(false)
    }
}
#[derive(Clone)]
enum ProgressStorage {
    Live(Lease),
    Owned(Rc<RefCell<ck::TautomerProgress<'static>>>),
}
enum ViewStorage {
    Source(Lease),
    Entry(ProgressStorage, usize),
    Owned(RefCell<ck::TautomerMoleculeView<'static>>),
}
#[wasm_bindgen]
pub struct TautomerMoleculeView {
    storage: ViewStorage,
}
#[wasm_bindgen]
pub struct TautomerProgress {
    storage: ProgressStorage,
}

// Private protocol tags do not infer behavior: every arm calls its facade method.
const SNAPSHOT: u32 = 0;
const COUNT: u32 = 1;
const BONDS: u32 = 2;
const SMILES: u32 = 3;
const SCORE: u32 = 4;
const DEGREE: u32 = 5;
const ATOMS: u32 = 6;
const BOND_VALUES: u32 = 7;
const ATOM: u32 = 8;
const BOND: u32 = 9;
const PROPERTIES: u32 = 10;
const METADATA: u32 = 11;
fn view_command(
    view: &mut ck::TautomerMoleculeView<'_>,
    operation: u32,
    index: usize,
) -> Result<JsValue, JsValue> {
    Ok(match operation {
        SNAPSHOT => TautomerMoleculeView {
            storage: ViewStorage::Owned(RefCell::new(view.snapshot())),
        }
        .into(),
        COUNT => (view.num_atoms() as u32).into(),
        BONDS => (view.num_bonds() as u32).into(),
        SMILES => crate::property_errors::text(
            &view
                .to_smiles()
                .map_err(crate::tautomer_errors::TautomerRunError::from_rust)?,
        )?
        .into(),
        SCORE => crate::tautomer_score::TautomerScore(
            view.tautomer_score()
                .map_err(crate::tautomer_errors::TautomerRunError::from_rust)?,
        )
        .into(),
        DEGREE => view
            .atom_degree(ck::AtomId::new(index))
            .map(|v| JsValue::from(v as u32))
            .unwrap_or(JsValue::UNDEFINED),
        ATOMS => Vec::<ck::Atom>::get(NativeField::owned(view.atoms().to_vec()))?.into(),
        BOND_VALUES => Vec::<ck::Bond>::get(NativeField::owned(view.bonds().to_vec()))?.into(),
        ATOM => view
            .atom(ck::AtomId::new(index))
            .map(|v| crate::atoms::Atom(v.clone()).into())
            .unwrap_or(JsValue::UNDEFINED),
        BOND => view
            .bond(ck::BondId::new(index))
            .map(|v| crate::atoms::Bond(v.clone()).into())
            .unwrap_or(JsValue::UNDEFINED),
        PROPERTIES => crate::metadata_values::MoleculeProperties(view.properties().clone()).into(),
        METADATA => js_sys::Object::freeze(
            &view
                .atom_metadata()
                .map_err(crate::valence_errors::ValenceError::from_rust)?
                .into_iter()
                .map(|v| JsValue::from(crate::metadata_queries::AtomMetadata(v)))
                .collect::<js_sys::Array>(),
        )
        .into(),
        _ => unreachable!("private view command"),
    })
}
fn progress_command(
    progress: &mut ck::TautomerProgress<'_>,
    operation: u32,
    index: usize,
) -> Result<JsValue, JsValue> {
    Ok(match operation {
        SNAPSHOT => TautomerProgress {
            storage: ProgressStorage::Owned(Rc::new(RefCell::new(progress.snapshot()))),
        }
        .into(),
        COUNT => (progress.len() as u32).into(),
        2 => TautomerEnumerationStatus::from(progress.status()).into(),
        3 => progress.num_transforms().into(),
        4 => Vec::<usize>::get(NativeField::owned(
            progress
                .modified_atoms()
                .iter()
                .map(|id| id.index())
                .collect(),
        ))?
        .js_value(),
        5 => Vec::<usize>::get(NativeField::owned(
            progress
                .modified_bonds()
                .iter()
                .map(|id| id.index())
                .collect(),
        ))?
        .js_value(),
        6 => progress
            .entry(index)
            .map(|(key, _)| crate::property_errors::text(key).map(JsValue::from))
            .transpose()?
            .unwrap_or(JsValue::UNDEFINED),
        _ => unreachable!("private progress command"),
    })
}
impl ProgressStorage {
    fn command(&self, operation: u32, index: usize) -> Result<JsValue, JsValue> {
        match self {
            Self::Live(lease) => lease.call(1, operation, index, 0),
            Self::Owned(owner) => progress_command(
                &mut *owner
                    .try_borrow_mut()
                    .map_err(|_| js_sys::TypeError::new("progress is currently borrowed"))?,
                operation,
                index,
            ),
        }
    }
    fn view(&self, index: usize, operation: u32, argument: usize) -> Result<JsValue, JsValue> {
        match self {
            Self::Live(lease) => lease.call(2, operation, index, argument),
            Self::Owned(owner) => {
                let mut owner = owner
                    .try_borrow_mut()
                    .map_err(|_| js_sys::TypeError::new("progress is currently borrowed"))?;
                let (_, mut view) = owner
                    .entry(index)
                    .ok_or_else(|| js_sys::ReferenceError::new("candidate is absent"))?;
                view_command(&mut view, operation, argument)
            }
        }
    }
}
impl TautomerMoleculeView {
    fn command(&self, operation: u32, index: usize) -> Result<JsValue, JsValue> {
        match &self.storage {
            ViewStorage::Source(lease) => lease.call(0, operation, index, 0),
            ViewStorage::Entry(progress, entry) => progress.view(*entry, operation, index),
            ViewStorage::Owned(owner) => view_command(
                &mut *owner
                    .try_borrow_mut()
                    .map_err(|_| js_sys::TypeError::new("snapshot is currently borrowed"))?,
                operation,
                index,
            ),
        }
    }
}
#[wasm_bindgen]
impl TautomerMoleculeView {
    pub fn snapshot(&self) -> Result<ViewValue, JsValue> {
        self.command(SNAPSHOT, 0).map(JsCast::unchecked_into)
    }
    #[wasm_bindgen(js_name=toOwned)]
    pub fn to_owned(&self) -> Result<ViewValue, JsValue> {
        self.snapshot()
    }
    #[wasm_bindgen(getter,js_name=numAtoms)]
    pub fn num_atoms(&self) -> Result<usize, JsValue> {
        Ok(self.command(COUNT, 0)?.as_f64().expect("native count") as usize)
    }
    #[wasm_bindgen(getter,js_name=numBonds)]
    pub fn num_bonds(&self) -> Result<usize, JsValue> {
        Ok(self.command(BONDS, 0)?.as_f64().expect("native count") as usize)
    }
    #[wasm_bindgen(js_name=toSmiles)]
    pub fn to_smiles(&self) -> Result<String, JsValue> {
        Ok(self.command(SMILES, 0)?.as_string().expect("native text"))
    }
    #[wasm_bindgen(js_name=tautomerScore)]
    pub fn tautomer_score(&self) -> Result<ScoreValue, JsValue> {
        self.command(SCORE, 0).map(JsCast::unchecked_into)
    }
    #[wasm_bindgen(js_name=atomDegree)]
    pub fn atom_degree(&self, index: Ts<input::IndexInput>) -> Result<Option<usize>, JsValue> {
        Ok(self
            .command(DEGREE, input::read::<usize>(index)?)?
            .as_f64()
            .map(|v| v as usize))
    }
    pub fn atoms(&self) -> Result<crate::atoms::Atoms, JsValue> {
        self.command(ATOMS, 0).map(JsCast::unchecked_into)
    }
    pub fn bonds(&self) -> Result<crate::atoms::Bonds, JsValue> {
        self.command(BOND_VALUES, 0).map(JsCast::unchecked_into)
    }
    pub fn atom(
        &self,
        index: Ts<input::IndexInput>,
    ) -> Result<Option<crate::atoms::AtomInput>, JsValue> {
        let v = self.command(ATOM, input::read::<usize>(index)?)?;
        Ok((!v.is_undefined()).then(|| v.unchecked_into()))
    }
    pub fn bond(
        &self,
        index: Ts<input::IndexInput>,
    ) -> Result<Option<crate::atoms::BondInput>, JsValue> {
        let v = self.command(BOND, input::read::<usize>(index)?)?;
        Ok((!v.is_undefined()).then(|| v.unchecked_into()))
    }
    #[wasm_bindgen(getter)]
    pub fn properties(&self) -> Result<PropertiesValue, JsValue> {
        self.command(PROPERTIES, 0).map(JsCast::unchecked_into)
    }
    #[wasm_bindgen(js_name=atomMetadata)]
    pub fn atom_metadata(&self) -> Result<MetadataValues, JsValue> {
        self.command(METADATA, 0).map(JsCast::unchecked_into)
    }
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "TautomerScore")]
    pub type ScoreValue;
    #[wasm_bindgen(typescript_type = "MoleculeProperties")]
    pub type PropertiesValue;
    #[wasm_bindgen(typescript_type = "ReadonlyArray<AtomMetadata>")]
    pub type MetadataValues;
    #[wasm_bindgen(typescript_type = "TautomerEnumerationStatus")]
    pub type StatusValue;
}
#[wasm_bindgen]
impl TautomerProgress {
    pub fn snapshot(&self) -> Result<ProgressValue, JsValue> {
        self.storage
            .command(SNAPSHOT, 0)
            .map(JsCast::unchecked_into)
    }
    #[wasm_bindgen(js_name=toOwned)]
    pub fn to_owned(&self) -> Result<ProgressValue, JsValue> {
        self.snapshot()
    }
    #[wasm_bindgen(getter)]
    pub fn len(&self) -> Result<usize, JsValue> {
        Ok(self
            .storage
            .command(COUNT, 0)?
            .as_f64()
            .expect("native count") as usize)
    }
    #[wasm_bindgen(js_name=isEmpty)]
    pub fn is_empty(&self) -> Result<bool, JsValue> {
        Ok(self.len()? == 0)
    }
    #[wasm_bindgen(getter)]
    pub fn status(&self) -> Result<StatusValue, JsValue> {
        self.storage.command(2, 0).map(JsCast::unchecked_into)
    }
    #[wasm_bindgen(getter,js_name=numTransforms)]
    pub fn num_transforms(&self) -> Result<u32, JsValue> {
        Ok(self.storage.command(3, 0)?.as_f64().expect("native count") as u32)
    }
    #[wasm_bindgen(getter,js_name=modifiedAtoms)]
    pub fn modified_atoms(&self) -> Result<Ts<crate::fields::IndexArray>, JsValue> {
        self.storage.command(4, 0).map(Ts::new_unchecked)
    }
    #[wasm_bindgen(getter,js_name=modifiedBonds)]
    pub fn modified_bonds(&self) -> Result<Ts<crate::fields::IndexArray>, JsValue> {
        self.storage.command(5, 0).map(Ts::new_unchecked)
    }
    pub fn entry(&self, index: Ts<input::IndexInput>) -> Result<Option<TautomerEntry>, JsValue> {
        self.entry_at(input::read::<usize>(index)?)
    }
    pub fn entries(&self) -> Result<TautomerProgressEntries, JsValue> {
        self.storage.command(COUNT, 0)?;
        Ok(TautomerProgressEntries(NativeField::owned((
            self.storage.clone(),
            0,
        ))))
    }
}
impl TautomerProgress {
    fn entry_at(&self, index: usize) -> Result<Option<TautomerEntry>, JsValue> {
        let key = self.storage.command(6, index)?;
        if key.is_undefined() {
            return Ok(None);
        }
        let view = TautomerMoleculeView {
            storage: ViewStorage::Entry(self.storage.clone(), index),
        };
        Ok(Some(
            js_sys::Object::freeze(&js_sys::Array::of2(&key, &view.into())).unchecked_into(),
        ))
    }
}
#[wasm_bindgen]
pub struct TautomerProgressEntries(NativeField<(ProgressStorage, usize)>);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "IteratorResult<readonly [string, TautomerMoleculeView], undefined>"
    )]
    pub type ProgressIteration;
}
#[wasm_bindgen]
impl TautomerProgressEntries {
    #[wasm_bindgen(js_name = "[Symbol.iterator]")]
    pub fn iterator(&self) -> Self {
        Self(self.0.clone())
    }
    pub fn next(&self) -> Result<ProgressIteration, JsValue> {
        let value = self.0.write(|(storage, index)| {
            let entry = TautomerProgress {
                storage: storage.clone(),
            }
            .entry_at(*index)?;
            if entry.is_some() {
                *index += 1;
            }
            Ok::<_, JsValue>(entry.map(JsValue::from))
        })??;
        Ok(crate::iterator::result(value)?.unchecked_into())
    }
}

fn scope<R>(
    mut source: ck::TautomerMoleculeView<'_>,
    mut progress: Option<ck::TautomerProgress<'_>>,
    run: impl FnOnce(Lease) -> Result<R, JsValue>,
) -> Result<R, JsValue> {
    let active = Rc::new(Cell::new(true));
    let _expire = Expire(active.clone());
    let mut visit =
        move |kind: u32, operation: u32, index: u32, argument: u32| -> Result<JsValue, JsValue> {
            match kind {
                0 => view_command(&mut source, operation, index as usize),
                1 => progress_command(
                    progress.as_mut().expect("progress callback"),
                    operation,
                    index as usize,
                ),
                2 => {
                    let (_, mut view) = progress
                        .as_mut()
                        .expect("progress callback")
                        .entry(index as usize)
                        .ok_or_else(|| js_sys::ReferenceError::new("candidate is absent"))?;
                    view_command(&mut view, operation, argument as usize)
                }
                _ => unreachable!("private lease target"),
            }
        };
    // Borrowed closure cannot outlive either native loan. The live flag gives
    // the public ReferenceError before calling wasm-bindgen's expired closure.
    let access = ScopedClosure::borrow_mut_aborting(&mut visit);
    run(Lease {
        access: access.as_ref().unchecked_ref::<js_sys::Function>().clone(),
        active,
    })
}
struct EnumerationAdapter(Arc<callback::Handle>);
impl ck::TautomerEnumerationCallback for EnumerationAdapter {
    fn should_continue(
        &self,
        source: ck::TautomerMoleculeView<'_>,
        progress: ck::TautomerProgress<'_>,
    ) -> Result<bool, ck::TautomerRunError> {
        self.0.ensure_thread().map_err(callback_failure)?;
        scope(source, Some(progress), |lease| {
            self.0
                .call(|| {
                    let args = js_sys::Array::new();
                    args.push(
                        &TautomerMoleculeView {
                            storage: ViewStorage::Source(lease.clone()),
                        }
                        .into(),
                    );
                    args.push(
                        &TautomerProgress {
                            storage: ProgressStorage::Live(lease),
                        }
                        .into(),
                    );
                    Ok(args)
                })
                .and_then(|v| input::read::<bool>(Ts::new_unchecked(v)))
        })
        .map_err(callback_failure)
    }
}
struct ScorerAdapter(Arc<callback::Handle>);
impl ck::TautomerScorer for ScorerAdapter {
    fn score(&self, source: ck::TautomerMoleculeView<'_>) -> Result<i32, ck::TautomerRunError> {
        self.0.ensure_thread().map_err(callback_failure)?;
        scope(source, None, |lease| {
            self.0
                .call(|| {
                    Ok(js_sys::Array::of1(
                        &TautomerMoleculeView {
                            storage: ViewStorage::Source(lease),
                        }
                        .into(),
                    ))
                })
                .and_then(|v| input::read::<i32>(Ts::new_unchecked(v)))
        })
        .map_err(callback_failure)
    }
}
fn callback_failure(error: JsValue) -> ck::TautomerRunError {
    callback::record(error);
    ck::TautomerRunError::Callback("JavaScript callback failed".into())
}
macro_rules! callback_field {
    ($native:ty,$input:ident,$adapter:ident) => {
        impl Field for $native {
            type Input = $input;
            type Output = $input;
            fn input(value: Self::Input) -> Result<Self, JsValue> {
                Self::decode(value.as_ref())
            }
            fn decode(value: &JsValue) -> Result<Self, JsValue> {
                let h = callback::Handle::new(value)?;
                let native: Self = Arc::new($adapter(h.clone()));
                h.bind(&native);
                Ok(native)
            }
            fn get(value: NativeField<Self>) -> Result<Self::Output, JsValue> {
                value.read(callback::original)?.map(JsCast::unchecked_into)
            }
        }
    };
}
callback_field!(
    Arc<dyn ck::TautomerEnumerationCallback>,
    EnumerationCallback,
    EnumerationAdapter
);
callback_field!(Arc<dyn ck::TautomerScorer>, ScorerCallback, ScorerAdapter);
