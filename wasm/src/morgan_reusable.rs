//! Reusable generators preserve native shared settings and provider capture.
use crate::{
    Molecule,
    fields::Field,
    fingerprint_output::{FingerprintAdditionalOutput, FingerprintAdditionalOutputInput},
    input,
    morgan_configuration::*,
    morgan_errors::MorganReadError,
    storage::NativeField,
};
use cosmolkit as ck;

#[wasm_bindgen]
impl MorganSettings {
    #[wasm_bindgen(js_name = setRadius)]
    pub fn set_radius_method(&mut self, value: <u32 as Field>::Input) -> Result<(), JsValue> {
        self.set_radius(value)
    }
    #[wasm_bindgen(js_name = setOnlyNonzeroInvariants)]
    pub fn set_only_nonzero_invariants_method(
        &mut self,
        value: <bool as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_only_nonzero_invariants(value)
    }
    #[wasm_bindgen(js_name = setIncludeRedundantEnvironments)]
    pub fn set_include_redundant_environments_method(
        &mut self,
        value: <bool as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_include_redundant_environments(value)
    }
    #[wasm_bindgen(js_name = setIncludeChirality)]
    pub fn set_include_chirality_method(
        &mut self,
        value: <bool as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_include_chirality(value)
    }
    #[wasm_bindgen(js_name = setCountSimulation)]
    pub fn set_count_simulation_method(
        &mut self,
        value: <bool as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_count_simulation(value)
    }
    #[wasm_bindgen(js_name = setFpSize)]
    pub fn set_fp_size_method(&mut self, value: <u32 as Field>::Input) -> Result<(), JsValue> {
        self.set_fp_size(value)
    }
    #[wasm_bindgen(js_name = setBitsPerFeature)]
    pub fn set_bits_per_feature_method(
        &mut self,
        value: <u32 as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_bits_per_feature(value)
    }
    #[wasm_bindgen(js_name = setCountBounds)]
    pub fn set_count_bounds_method(
        &mut self,
        value: <Vec<u32> as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_count_bounds(value)
    }
}

use wasm_bindgen::prelude::*;
ck::__ck_configuration_AtomPairAtomInvariantsGenerator!(
    crate::project_parameters,
    ck::AtomPairAtomInvariantsGenerator,
    {},
    { AtomPairAtomInvariantsGenerator }
);
ck::__ck_configuration_MorganBondInvariantsGenerator!(
    crate::project_parameters,
    ck::MorganBondInvariantsGenerator,
    {},
    { MorganBondInvariantsGenerator }
);
impl MorganBondInvariantsGenerator {
    pub(crate) fn with_optional_ref<R>(
        value: Option<&MorganBondInvariantsGeneratorInput>,
        call: impl FnOnce(Option<&ck::MorganBondInvariantsGenerator>) -> Result<R, JsValue>,
    ) -> Result<R, JsValue> {
        match value {
            Some(value) => Self::with_input(Some(value), |v| call(Some(v))),
            None => call(None),
        }
    }
}
#[wasm_bindgen]
impl AtomPairAtomInvariantsGenerator {
    #[wasm_bindgen(js_name=infoString)]
    pub fn info_string(&self) -> Result<String, JsValue> {
        self.inner
            .read(ck::AtomPairAtomInvariantsGenerator::info_string)
    }
    #[wasm_bindgen(js_name=toJson)]
    pub fn to_json(&self) -> Result<String, JsValue> {
        self.inner
            .read(ck::AtomPairAtomInvariantsGenerator::to_json)
    }
}
#[wasm_bindgen]
pub struct MorganAtomInvariantsGenerator(pub(crate) ck::MorganAtomInvariantsGenerator);
crate::fields::value_field!(
    ck::MorganAtomInvariantsGenerator,
    MorganAtomInvariantsGenerator,
    MorganAtomInvariantsGeneratorInput,
    "MorganAtomInvariantsGenerator"
);
#[wasm_bindgen]
impl MorganAtomInvariantsGenerator {
    pub fn connectivity(
        include_ring_membership: Option<tsify::Ts<input::BoolInput>>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::MorganAtomInvariantsGenerator::connectivity(
            include_ring_membership
                .map(input::read::<bool>)
                .transpose()?
                .unwrap_or(true),
        )))
    }
    pub fn features(
        patterns: Option<crate::query_values::QueryGraphArray>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::MorganAtomInvariantsGenerator::features(
            patterns
                .map(<Vec<ck::QueryGraph> as Field>::input)
                .transpose()?,
        )))
    }
    #[wasm_bindgen(js_name=atomPair)]
    pub fn atom_pair(generator: AtomPairAtomInvariantsGeneratorInput) -> Result<Self, JsValue> {
        AtomPairAtomInvariantsGenerator::with_input(Some(&generator), |native| {
            Ok(Self(ck::MorganAtomInvariantsGenerator::atom_pair(*native)))
        })
    }
}
#[wasm_bindgen]
pub struct MorganFingerprintGenerator(pub(crate) ck::MorganFingerprintGenerator);
crate::fields::value_field!(
    ck::MorganFingerprintGenerator,
    MorganFingerprintGenerator,
    MorganFingerprintGeneratorInput,
    "MorganFingerprintGenerator"
);
#[wasm_bindgen]
impl MorganFingerprintGenerator {
    #[wasm_bindgen(constructor)]
    pub fn new(
        params: Option<MorganGeneratorParamsInput>,
        atom_invariants: Option<MorganAtomInvariantsGeneratorInput>,
        bond_invariants: Option<MorganBondInvariantsGeneratorInput>,
    ) -> Result<Self, JsValue> {
        MorganGeneratorParams::with_input(params.as_ref(), |params| {
            MorganAtomInvariantsGenerator::with_optional_ref(atom_invariants.as_ref(), |atom| {
                let make = |bond: Option<&ck::MorganBondInvariantsGenerator>| {
                    ck::MorganFingerprintGenerator::new(Some(params), atom, bond)
                        .map(Self)
                        .map_err(MorganReadError::from_rust)
                };
                match bond_invariants.as_ref() {
                    Some(bond) => MorganBondInvariantsGenerator::with_input(Some(bond), |bond| {
                        make(Some(bond))
                    }),
                    None => make(None),
                }
            })
        })
    }
    #[wasm_bindgen(js_name=fromJson)]
    pub fn from_json(json: tsify::Ts<input::TextInput>) -> Result<Self, JsValue> {
        ck::MorganFingerprintGenerator::from_json(&input::read::<String>(json)?)
            .map(Self)
            .map_err(MorganReadError::from_rust)
    }
    #[wasm_bindgen(getter)]
    pub fn settings(&self) -> MorganSettings {
        MorganSettings(self.0.settings())
    }
    #[wasm_bindgen(js_name=infoString)]
    pub fn info_string(&self) -> Result<String, JsValue> {
        self.0.info_string().map_err(MorganReadError::from_rust)
    }
    #[wasm_bindgen(js_name=toJson)]
    pub fn to_json(&self) -> Result<String, JsValue> {
        self.0
            .to_json()
            .map_err(MorganReadError::from_rust)
            .and_then(|v| crate::property_errors::text(&v))
    }
}
#[wasm_bindgen]
pub struct MorganSettings(ck::MorganSettings);
macro_rules! settings {
    ($($get:ident,$set:ident,$js:ident,$ty:ty);+$(;)?)=>{
        #[wasm_bindgen]
        impl MorganSettings {
            $(#[wasm_bindgen(getter,js_name=$js)] pub fn $get(&self)->Result<<$ty as Field>::Output,JsValue>{
                <$ty as Field>::get(NativeField::owned(self.0.$get().map_err(MorganReadError::from_rust)?))
            }
            #[wasm_bindgen(setter,js_name=$js)] pub fn $set(&mut self,value:<$ty as Field>::Input)->Result<(),JsValue>{
                self.0.$set(<$ty as Field>::input(value)?).map_err(MorganReadError::from_rust)
            })+
        }
    }
}
settings! {
    radius,set_radius,radius,u32;
    only_nonzero_invariants,set_only_nonzero_invariants,onlyNonzeroInvariants,bool;
    include_redundant_environments,set_include_redundant_environments,includeRedundantEnvironments,bool;
    include_chirality,set_include_chirality,includeChirality,bool;
    count_simulation,set_count_simulation,countSimulation,bool;
    fp_size,set_fp_size,fpSize,u32;
    bits_per_feature,set_bits_per_feature,bitsPerFeature,u32;
    count_bounds,set_count_bounds,countBounds,Vec<u32>;
}
#[wasm_bindgen]
impl MorganSettings {
    #[wasm_bindgen(getter)]
    pub fn params(&self) -> Result<MorganGeneratorParams, JsValue> {
        <ck::MorganGeneratorParams as Field>::get(NativeField::owned(
            self.0.params().map_err(MorganReadError::from_rust)?,
        ))
    }
}
macro_rules! call {
    ($method:ident,$js:ident,$output:path) => {
        #[wasm_bindgen]
        impl Molecule {
            #[wasm_bindgen(js_name=$js)]
            pub fn $method(
                &self,
                generator: MorganFingerprintGeneratorInput,
                params: Option<MorganCallParamsInput>,
                additional_output: Option<FingerprintAdditionalOutputInput>,
            ) -> Result<$output, JsValue> {
                MorganFingerprintGenerator::with_ref(&generator, |generator| {
                    MorganCallParams::with_input(params.as_ref(), |params| {
                        FingerprintAdditionalOutput::with_optional(
                            additional_output.as_ref(),
                            |output| {
                                self.0
                                    .$method(generator, Some(params), output)
                                    .map($output)
                                    .map_err(MorganReadError::from_rust)
                            },
                        )
                    })
                })
            }
        }
    };
}
call!(
    fingerprint_morgan_with_generator,
    fingerprintMorganWithGenerator,
    crate::fingerprint_values::Fingerprint
);
call!(
    fingerprint_morgan_count_with_generator,
    fingerprintMorganCountWithGenerator,
    crate::sparse_fingerprints::SparseCountFingerprint32
);
call!(
    fingerprint_morgan_sparse_with_generator,
    fingerprintMorganSparseWithGenerator,
    crate::sparse_fingerprints::SparseBitFingerprint
);
call!(
    fingerprint_morgan_sparse_count_with_generator,
    fingerprintMorganSparseCountWithGenerator,
    crate::sparse_fingerprints::SparseCountFingerprint
);

#[wasm_bindgen]
impl MorganFingerprintGenerator {
    pub fn fingerprints(
        &self,
        molecules: crate::MoleculeArray,
        num_threads: Option<tsify::Ts<input::I32Input>>,
    ) -> Result<crate::fingerprint_batch::FingerprintBatch, JsValue> {
        let molecules = molecules.native_values()?;
        let refs = molecules.iter().map(Option::as_ref).collect::<Vec<_>>();
        self.0
            .fingerprints(
                &refs,
                num_threads
                    .map(input::read::<i32>)
                    .transpose()?
                    .unwrap_or(1),
            )
            .map(crate::fingerprint_batch::FingerprintBatch)
            .map_err(MorganReadError::from_rust)
    }
}
macro_rules! bulk {
    ($method:ident,$js:ident,$array:path) => {
        #[wasm_bindgen]
        impl MorganFingerprintGenerator {
            #[wasm_bindgen(js_name=$js)]
            pub fn $method(
                &self,
                molecules: crate::MoleculeArray,
                num_threads: Option<tsify::Ts<input::I32Input>>,
            ) -> Result<$array, JsValue> {
                let molecules = molecules.native_values()?;
                let refs = molecules.iter().map(Option::as_ref).collect::<Vec<_>>();
                self.0
                    .$method(
                        &refs,
                        num_threads
                            .map(input::read::<i32>)
                            .transpose()?
                            .unwrap_or(1),
                    )
                    .map(<$array>::from_native)
                    .map_err(MorganReadError::from_rust)
            }
        }
    };
}
bulk!(
    counts,
    counts,
    crate::sparse_fingerprints::SparseCountFingerprint32Array
);
bulk!(
    sparse_fingerprints,
    sparseFingerprints,
    crate::sparse_fingerprints::SparseBitFingerprintArray
);
bulk!(
    sparse_counts,
    sparseCounts,
    crate::sparse_fingerprints::SparseCountFingerprintArray
);
