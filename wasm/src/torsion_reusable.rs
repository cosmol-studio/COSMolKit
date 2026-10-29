//! Persistent torsion executors retain their native settings and providers.
use crate::{
    Molecule,
    fields::Field,
    fingerprint_output::{FingerprintAdditionalOutput, FingerprintAdditionalOutputInput},
    input,
    morgan_reusable::{AtomPairAtomInvariantsGenerator, AtomPairAtomInvariantsGeneratorInput},
    pair_torsion::{
        TopologicalTorsionGeneratorParams, TopologicalTorsionGeneratorParamsInput,
        TopologicalTorsionReadError,
    },
    storage::NativeField,
};
use cosmolkit as ck;

#[wasm_bindgen]
impl TopologicalTorsionSettings {
    #[wasm_bindgen(js_name = setTorsionAtomCount)]
    pub fn set_torsion_atom_count_method(
        &mut self,
        value: <u32 as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_torsion_atom_count(value)
    }
    #[wasm_bindgen(js_name = setOnlyShortestPaths)]
    pub fn set_only_shortest_paths_method(
        &mut self,
        value: <bool as Field>::Input,
    ) -> Result<(), JsValue> {
        self.set_only_shortest_paths(value)
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
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=topologicalTorsionIdsWithParams)]
    pub fn topological_torsion_ids_with_params(
        &self,
        torsion_atom_count: tsify::Ts<input::U32Input>,
    ) -> Result<Vec<u64>, JsValue> {
        self.topological_torsion_ids(Some(torsion_atom_count))
    }
    #[wasm_bindgen(js_name=topologicalTorsionIds)]
    pub fn topological_torsion_ids(
        &self,
        torsion_atom_count: Option<tsify::Ts<input::U32Input>>,
    ) -> Result<Vec<u64>, JsValue> {
        match torsion_atom_count {
            Some(v) => self.0.topological_torsion_ids_with_params(input::read(v)?),
            None => self.0.topological_torsion_ids(),
        }
        .map_err(TopologicalTorsionReadError::from_rust)
    }
}
ck::__ck_configuration_LegacyTopologicalTorsionParams!(
    crate::project_parameters,
    ck::LegacyTopologicalTorsionParams,
    {},
    { LegacyTopologicalTorsionParams }
);
macro_rules! legacy_call {
    ($method:ident,$native:ident,$js:ident,$out:path) => {
        #[wasm_bindgen]
        impl Molecule {
            #[wasm_bindgen(js_name=$js)]
            pub fn $method(
                &self,
                params: Option<LegacyTopologicalTorsionParamsInput>,
            ) -> Result<$out, JsValue> {
                LegacyTopologicalTorsionParams::with_input(params.as_ref(), |p| {
                    self.0
                        .$native(p)
                        .map($out)
                        .map_err(TopologicalTorsionReadError::from_rust)
                })
            }
        }
        paste::paste! {
            #[wasm_bindgen]
            impl Molecule {
                #[wasm_bindgen(js_name=[<$js WithParams>])]
                pub fn $native(&self, params: LegacyTopologicalTorsionParamsInput) -> Result<$out, JsValue> {
                    self.$method(Some(params))
                }
            }
        }
    };
}
legacy_call!(
    fingerprint_topological_torsion_legacy,
    fingerprint_topological_torsion_legacy_with_params,
    fingerprintTopologicalTorsionLegacy,
    crate::fingerprint_values::Fingerprint
);
legacy_call!(
    fingerprint_topological_torsion_count_legacy,
    fingerprint_topological_torsion_count_legacy_with_params,
    fingerprintTopologicalTorsionCountLegacy,
    crate::sparse_fingerprints::SparseCountFingerprint
);
legacy_call!(
    fingerprint_topological_torsion_sparse_count_legacy,
    fingerprint_topological_torsion_sparse_count_legacy_with_params,
    fingerprintTopologicalTorsionSparseCountLegacy,
    crate::sparse_fingerprints::SparseCountFingerprint
);
ck::__ck_configuration_TopologicalTorsionCallParams!(
    crate::project_parameters,
    ck::TopologicalTorsionCallParams,
    {},
    { TopologicalTorsionCallParams }
);
#[wasm_bindgen]
pub struct TopologicalTorsionFingerprintGenerator(
    pub(crate) ck::TopologicalTorsionFingerprintGenerator,
);
crate::fields::value_field!(
    ck::TopologicalTorsionFingerprintGenerator,
    TopologicalTorsionFingerprintGenerator,
    TopologicalTorsionFingerprintGeneratorInput,
    "TopologicalTorsionFingerprintGenerator"
);
#[wasm_bindgen]
impl TopologicalTorsionFingerprintGenerator {
    #[wasm_bindgen(constructor)]
    pub fn new(
        params: Option<TopologicalTorsionGeneratorParamsInput>,
        atom_invariants: Option<AtomPairAtomInvariantsGeneratorInput>,
    ) -> Result<Self, JsValue> {
        TopologicalTorsionGeneratorParams::with_input(params.as_ref(), |p| {
            let make = |a| {
                ck::TopologicalTorsionFingerprintGenerator::new(Some(p), a)
                    .map(Self)
                    .map_err(TopologicalTorsionReadError::from_rust)
            };
            match atom_invariants.as_ref() {
                Some(a) => AtomPairAtomInvariantsGenerator::with_input(Some(a), |a| make(Some(*a))),
                None => make(None),
            }
        })
    }
    #[wasm_bindgen(js_name=fromJson)]
    pub fn from_json(json: tsify::Ts<input::TextInput>) -> Result<Self, JsValue> {
        ck::TopologicalTorsionFingerprintGenerator::from_json(&input::read::<String>(json)?)
            .map(Self)
            .map_err(TopologicalTorsionReadError::from_rust)
    }
    #[wasm_bindgen(getter)]
    pub fn settings(&self) -> TopologicalTorsionSettings {
        TopologicalTorsionSettings(self.0.settings())
    }
    #[wasm_bindgen(js_name=infoString)]
    pub fn info_string(&self) -> Result<String, JsValue> {
        self.0
            .info_string()
            .map_err(TopologicalTorsionReadError::from_rust)
    }
    #[wasm_bindgen(js_name=toJson)]
    pub fn to_json(&self) -> Result<String, JsValue> {
        self.0
            .to_json()
            .map_err(TopologicalTorsionReadError::from_rust)
    }
    pub fn fingerprints(
        &self,
        molecules: crate::MoleculeArray,
        num_threads: Option<tsify::Ts<input::I32Input>>,
    ) -> Result<crate::fingerprint_batch::FingerprintBatch, JsValue> {
        let values = molecules.native_values()?;
        let refs = values.iter().map(Option::as_ref).collect::<Vec<_>>();
        self.0
            .fingerprints(
                &refs,
                num_threads
                    .map(input::read::<i32>)
                    .transpose()?
                    .unwrap_or(1),
            )
            .map(crate::fingerprint_batch::FingerprintBatch)
            .map_err(TopologicalTorsionReadError::from_rust)
    }
}
#[wasm_bindgen]
pub struct TopologicalTorsionSettings(ck::TopologicalTorsionSettings);
macro_rules! settings {($($get:ident,$set:ident,$js:ident,$ty:ty);+$(;)?)=>{#[wasm_bindgen]impl TopologicalTorsionSettings{$(
    #[wasm_bindgen(getter,js_name=$js)]pub fn $get(&self)->Result<<$ty as Field>::Output,JsValue>{<$ty as Field>::get(NativeField::owned(self.0.$get().map_err(TopologicalTorsionReadError::from_rust)?))}
    #[wasm_bindgen(setter,js_name=$js)]pub fn $set(&mut self,value:<$ty as Field>::Input)->Result<(),JsValue>{self.0.$set(<$ty as Field>::input(value)?).map_err(TopologicalTorsionReadError::from_rust)}
)+}};}
settings! {
    torsion_atom_count,set_torsion_atom_count,torsionAtomCount,u32;
    only_shortest_paths,set_only_shortest_paths,onlyShortestPaths,bool;
    include_chirality,set_include_chirality,includeChirality,bool;
    count_simulation,set_count_simulation,countSimulation,bool;
    fp_size,set_fp_size,fpSize,u32;
    bits_per_feature,set_bits_per_feature,bitsPerFeature,u32;
    count_bounds,set_count_bounds,countBounds,Vec<u32>;
}
#[wasm_bindgen]
impl TopologicalTorsionSettings {
    #[wasm_bindgen(getter)]
    pub fn params(&self) -> Result<TopologicalTorsionGeneratorParams, JsValue> {
        <ck::TopologicalTorsionGeneratorParams as Field>::get(NativeField::owned(
            self.0
                .params()
                .map_err(TopologicalTorsionReadError::from_rust)?,
        ))
    }
}
macro_rules! call {
    ($method:ident,$js:ident,$out:path) => {
        #[wasm_bindgen]
        impl Molecule {
            #[wasm_bindgen(js_name=$js)]
            pub fn $method(
                &self,
                generator: TopologicalTorsionFingerprintGeneratorInput,
                params: Option<TopologicalTorsionCallParamsInput>,
                additional_output: Option<FingerprintAdditionalOutputInput>,
            ) -> Result<$out, JsValue> {
                TopologicalTorsionFingerprintGenerator::with_ref(&generator, |g| {
                    TopologicalTorsionCallParams::with_input(params.as_ref(), |p| {
                        FingerprintAdditionalOutput::with_optional(
                            additional_output.as_ref(),
                            |o| {
                                self.0
                                    .$method(g, Some(p), o)
                                    .map($out)
                                    .map_err(TopologicalTorsionReadError::from_rust)
                            },
                        )
                    })
                })
            }
        }
    };
}
call!(
    fingerprint_topological_torsion_with_generator,
    fingerprintTopologicalTorsionWithGenerator,
    crate::fingerprint_values::Fingerprint
);
call!(
    fingerprint_topological_torsion_sparse_with_generator,
    fingerprintTopologicalTorsionSparseWithGenerator,
    crate::sparse_fingerprints::SparseBitFingerprint
);
call!(
    fingerprint_topological_torsion_count_with_generator,
    fingerprintTopologicalTorsionCountWithGenerator,
    crate::sparse_fingerprints::SparseCountFingerprint32
);
call!(
    fingerprint_topological_torsion_sparse_count_with_generator,
    fingerprintTopologicalTorsionSparseCountWithGenerator,
    crate::sparse_fingerprints::SparseCountFingerprint
);
macro_rules! bulk {
    ($method:ident,$js:ident,$array:path) => {
        #[wasm_bindgen]
        impl TopologicalTorsionFingerprintGenerator {
            #[wasm_bindgen(js_name=$js)]
            pub fn $method(
                &self,
                molecules: crate::MoleculeArray,
                num_threads: Option<tsify::Ts<input::I32Input>>,
            ) -> Result<$array, JsValue> {
                let values = molecules.native_values()?;
                let refs = values.iter().map(Option::as_ref).collect::<Vec<_>>();
                self.0
                    .$method(
                        &refs,
                        num_threads
                            .map(input::read::<i32>)
                            .transpose()?
                            .unwrap_or(1),
                    )
                    .map(<$array>::from_native)
                    .map_err(TopologicalTorsionReadError::from_rust)
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
