use crate::fingerprint_providers::AtomPairAtomInvariantsGenerator;
use cosmolkit as ck;
use pyo3::prelude::*;

// These generator records share one representation rule: optional count bounds
// select the actual Rust default; a present empty vector remains present/empty.
// No field defaults or fingerprint algorithms are copied into this declaration.
macro_rules! generator_parameters {
    ($name:ident, {$($field:ident: $ty:ty),+ $(,)?}) => {
        #[pyclass(module="fixed_projection._native")]
        pub struct $name { pub(crate) inner: ck::$name }
        paste::paste! {
            #[pymethods]
            impl $name {
                #[new]
                #[pyo3(signature=(*, $($field=ck::$name::default().$field),+, count_bounds=None))]
                fn new($($field: $ty),+, count_bounds: Option<Vec<u32>>) -> Self {
                    Self { inner: ck::$name {
                        $($field),+, count_bounds: count_bounds.unwrap_or_else(|| ck::$name::default().count_bounds),
                    } }
                }
                $(#[getter]
                fn $field(&self) -> $ty { self.inner.$field.clone() }
                #[setter]
                fn [<set_ $field>](&mut self, value: $ty) { self.inner.$field = value; })+
                #[getter]
                fn count_bounds(&self) -> Vec<u32> { self.inner.count_bounds.clone() }
                #[setter]
                fn set_count_bounds(&mut self, value: Option<Vec<u32>>) {
                    self.inner.count_bounds = value.unwrap_or_else(|| ck::$name::default().count_bounds);
                }
                fn info_string(&self) -> String { self.inner.info_string() }
                fn to_json(&self) -> String { self.inner.to_json() }
                fn with_json(&self, json: &str) -> PyResult<Self> {
                    self.inner.with_json(json).map(|inner| Self { inner })
                        .map_err(crate::fingerprint_json::FingerprintJsonError::from_rust)
                }
            }
        }
    };
}

generator_parameters!(MorganParams, {
    radius: u32, include_chirality: bool, use_bond_types: bool,
    include_ring_membership: bool, only_nonzero_invariants: bool,
    include_redundant_environments: bool, fp_size: u32,
    count_simulation: bool, bits_per_feature: u32,
});
generator_parameters!(AtomPairParams, {
    min_distance: u32, max_distance: u32, include_chirality: bool,
    use_2d: bool, count_simulation: bool, fp_size: u32, bits_per_feature: u32,
});
generator_parameters!(TopologicalTorsionParams, {
    torsion_atom_count: u32, only_shortest_paths: bool, include_chirality: bool,
    count_simulation: bool, fp_size: u32, bits_per_feature: u32,
});

// Shared source facade call transport: owned nested generator values and
// genuine optional selections. Supplied mutable Params never become aliases.
macro_rules! call_parameters {
    ($name:ident, $generator:ident) => {
        #[pyclass(module="fixed_projection._native")]
        pub struct $name { pub(crate) inner: ck::$name }
        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature=(*, generator=None, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None, custom_bond_invariants=None, conformer_id=ck::$name::default().conformer_id, atom_invariants_generator=None, use_legacy_stereo_perception=ck::$name::default().use_legacy_stereo_perception))]
            fn new(generator: Option<&$generator>, from_atoms: Option<Vec<u32>>,
                ignore_atoms: Option<Vec<u32>>, custom_atom_invariants: Option<Vec<u32>>,
                custom_bond_invariants: Option<Vec<u32>>, conformer_id: i32,
                atom_invariants_generator: Option<&AtomPairAtomInvariantsGenerator>,
                use_legacy_stereo_perception: bool) -> Self {
                Self { inner: ck::$name {
                    generator: generator.map(|value| value.inner.clone()).unwrap_or_else(|| ck::$name::default().generator),
                    from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants,
                    conformer_id,
                    atom_invariants_generator: atom_invariants_generator.map(|value| value.inner),
                    use_legacy_stereo_perception,
                } }
            }
            #[getter]
            fn generator(&self) -> $generator { $generator { inner: self.inner.generator.clone() } }
            #[setter]
            fn set_generator(&mut self, value: Option<&$generator>) {
                self.inner.generator = value.map(|value| value.inner.clone()).unwrap_or_else(|| ck::$name::default().generator);
            }
            #[getter]
            fn from_atoms(&self) -> Option<Vec<u32>> { self.inner.from_atoms.clone() }
            #[setter]
            fn set_from_atoms(&mut self, value: Option<Vec<u32>>) { self.inner.from_atoms = value; }
            #[getter]
            fn ignore_atoms(&self) -> Option<Vec<u32>> { self.inner.ignore_atoms.clone() }
            #[setter]
            fn set_ignore_atoms(&mut self, value: Option<Vec<u32>>) { self.inner.ignore_atoms = value; }
            #[getter]
            fn custom_atom_invariants(&self) -> Option<Vec<u32>> { self.inner.custom_atom_invariants.clone() }
            #[setter]
            fn set_custom_atom_invariants(&mut self, value: Option<Vec<u32>>) { self.inner.custom_atom_invariants = value; }
            #[getter]
            fn custom_bond_invariants(&self) -> Option<Vec<u32>> { self.inner.custom_bond_invariants.clone() }
            #[setter]
            fn set_custom_bond_invariants(&mut self, value: Option<Vec<u32>>) { self.inner.custom_bond_invariants = value; }
            #[getter]
            fn conformer_id(&self) -> i32 { self.inner.conformer_id }
            #[setter]
            fn set_conformer_id(&mut self, value: i32) { self.inner.conformer_id = value; }
            #[getter]
            fn atom_invariants_generator(&self) -> Option<AtomPairAtomInvariantsGenerator> {
                self.inner.atom_invariants_generator.map(|inner| AtomPairAtomInvariantsGenerator { inner })
            }
            #[setter]
            fn set_atom_invariants_generator(&mut self, value: Option<&AtomPairAtomInvariantsGenerator>) {
                self.inner.atom_invariants_generator = value.map(|value| value.inner);
            }
            #[getter]
            fn use_legacy_stereo_perception(&self) -> bool { self.inner.use_legacy_stereo_perception }
            #[setter]
            fn set_use_legacy_stereo_perception(&mut self, value: bool) { self.inner.use_legacy_stereo_perception = value; }
        }
    }
}
call_parameters!(AtomPairFingerprintParams, AtomPairParams);
call_parameters!(
    TopologicalTorsionFingerprintParams,
    TopologicalTorsionParams
);
