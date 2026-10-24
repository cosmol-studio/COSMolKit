//! Ordinary fingerprint calls: construct native options; never regenerate algorithms.
use crate::fingerprint_configuration::{
    AtomPairFingerprintParams, AtomPairGeneratorParams, MorganGeneratorParams,
    TopologicalTorsionFingerprintParams, TopologicalTorsionGeneratorParams,
};
use crate::fingerprint_providers::AtomPairAtomInvariantsGenerator;
use crate::fingerprint_providers::FingerprintAdditionalOutput;
use crate::fingerprint_values::{
    SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32,
};
use crate::morgan_configuration::{MorganFingerprintParams, MorganInvariants};
use crate::native::{Fingerprint, Molecule};
use crate::parameters::{ConfigurationField, FieldInput, Scalar};
use cosmolkit as ck;
use pyo3::prelude::*;

macro_rules! run_fingerprint {
    ($name:ident, $method:ident, $params:ident, $output:ident, $error:ident) => {
        pub(crate) fn $name(
            molecule: &Molecule,
            py: Python<'_>,
            params: &ck::$params,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<$output> {
            let additional_output = additional_output.map(|value| &mut value.inner);
            py.detach(|| molecule.inner.$method(params, additional_output))
                .map(|inner| $output { inner })
                .map_err(crate::morgan_errors::$error::from_rust)
        }
    };
}

pub(crate) fn morgan_keywords(
    radius: FieldInput<Scalar<u32>>,
    include_chirality: FieldInput<Scalar<bool>>,
    use_bond_types: FieldInput<Scalar<bool>>,
    include_ring_membership: FieldInput<Scalar<bool>>,
    only_nonzero_invariants: FieldInput<Scalar<bool>>,
    include_redundant_environments: FieldInput<Scalar<bool>>,
    fp_size: FieldInput<Scalar<u32>>,
    count_simulation: FieldInput<Scalar<bool>>,
    bits_per_feature: FieldInput<Scalar<u32>>,
    count_bounds: FieldInput<Vec<u32>>,
    from_atoms: FieldInput<Option<Vec<u32>>>,
    ignore_atoms: FieldInput<Option<Vec<u32>>>,
    custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
    custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
    invariants: FieldInput<ConfigurationField<MorganInvariants>>,
) -> ck::MorganFingerprintParams {
    let mut params = ck::MorganFingerprintParams::default();
    if let FieldInput::Value(value) = radius {
        params.set_radius(value);
    }
    if let FieldInput::Value(value) = include_chirality {
        params.set_include_chirality(value);
    }
    if let FieldInput::Value(value) = use_bond_types {
        params.set_use_bond_types(value);
    }
    if let FieldInput::Value(value) = include_ring_membership {
        params.set_include_ring_membership(value);
    }
    if let FieldInput::Value(value) = only_nonzero_invariants {
        params.set_only_nonzero_invariants(value);
    }
    if let FieldInput::Value(value) = include_redundant_environments {
        params.set_include_redundant_environments(value);
    }
    if let FieldInput::Value(value) = fp_size {
        params.set_fp_size(value);
    }
    if let FieldInput::Value(value) = count_simulation {
        params.set_count_simulation(value);
    }
    if let FieldInput::Value(value) = bits_per_feature {
        params.set_bits_per_feature(value);
    }
    if let FieldInput::Value(value) = count_bounds {
        params.set_count_bounds(value);
    }
    if let FieldInput::Value(value) = from_atoms {
        params.from_atoms = value;
    }
    if let FieldInput::Value(value) = ignore_atoms {
        params.ignore_atoms = value;
    }
    if let FieldInput::Value(value) = custom_atom_invariants {
        params.custom_atom_invariants = value;
    }
    if let FieldInput::Value(value) = custom_bond_invariants {
        params.custom_bond_invariants = value;
    }
    if let FieldInput::Value(value) = invariants {
        params.invariants = value;
    }
    params
}

pub(crate) fn atom_pair_keywords(
    min_distance: FieldInput<Scalar<u32>>,
    max_distance: FieldInput<Scalar<u32>>,
    include_chirality: FieldInput<Scalar<bool>>,
    use_2d: FieldInput<Scalar<bool>>,
    count_simulation: FieldInput<Scalar<bool>>,
    fp_size: FieldInput<Scalar<u32>>,
    bits_per_feature: FieldInput<Scalar<u32>>,
    count_bounds: FieldInput<Vec<u32>>,
    from_atoms: FieldInput<Option<Vec<u32>>>,
    ignore_atoms: FieldInput<Option<Vec<u32>>>,
    custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
    custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
    conformer_id: FieldInput<Option<Scalar<usize>>>,
    atom_invariants_generator: FieldInput<
        Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>,
    >,
    use_legacy_stereo_perception: FieldInput<Scalar<bool>>,
) -> ck::AtomPairFingerprintParams {
    let mut params = ck::AtomPairFingerprintParams::default();
    if let FieldInput::Value(value) = min_distance {
        params.set_min_distance(value);
    }
    if let FieldInput::Value(value) = max_distance {
        params.set_max_distance(value);
    }
    if let FieldInput::Value(value) = include_chirality {
        params.set_include_chirality(value);
    }
    if let FieldInput::Value(value) = use_2d {
        params.set_use_2d(value);
    }
    if let FieldInput::Value(value) = count_simulation {
        params.set_count_simulation(value);
    }
    if let FieldInput::Value(value) = fp_size {
        params.set_fp_size(value);
    }
    if let FieldInput::Value(value) = bits_per_feature {
        params.set_bits_per_feature(value);
    }
    if let FieldInput::Value(value) = count_bounds {
        params.set_count_bounds(value);
    }
    if let FieldInput::Value(value) = from_atoms {
        params.from_atoms = value;
    }
    if let FieldInput::Value(value) = ignore_atoms {
        params.ignore_atoms = value;
    }
    if let FieldInput::Value(value) = custom_atom_invariants {
        params.custom_atom_invariants = value;
    }
    if let FieldInput::Value(value) = custom_bond_invariants {
        params.custom_bond_invariants = value;
    }
    if let FieldInput::Value(value) = conformer_id {
        params.conformer_id = value;
    }
    if let FieldInput::Value(value) = atom_invariants_generator {
        params.atom_invariants_generator = value;
    }
    if let FieldInput::Value(value) = use_legacy_stereo_perception {
        params.use_legacy_stereo_perception = value;
    }
    params
}

pub(crate) fn topological_torsion_keywords(
    torsion_atom_count: FieldInput<Scalar<u32>>,
    only_shortest_paths: FieldInput<Scalar<bool>>,
    include_chirality: FieldInput<Scalar<bool>>,
    count_simulation: FieldInput<Scalar<bool>>,
    fp_size: FieldInput<Scalar<u32>>,
    bits_per_feature: FieldInput<Scalar<u32>>,
    count_bounds: FieldInput<Vec<u32>>,
    from_atoms: FieldInput<Option<Vec<u32>>>,
    ignore_atoms: FieldInput<Option<Vec<u32>>>,
    custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
    custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
    atom_invariants_generator: FieldInput<
        Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>,
    >,
    use_legacy_stereo_perception: FieldInput<Scalar<bool>>,
) -> ck::TopologicalTorsionFingerprintParams {
    let mut params = ck::TopologicalTorsionFingerprintParams::default();
    if let FieldInput::Value(value) = torsion_atom_count {
        params.set_torsion_atom_count(value);
    }
    if let FieldInput::Value(value) = only_shortest_paths {
        params.set_only_shortest_paths(value);
    }
    if let FieldInput::Value(value) = include_chirality {
        params.set_include_chirality(value);
    }
    if let FieldInput::Value(value) = count_simulation {
        params.set_count_simulation(value);
    }
    if let FieldInput::Value(value) = fp_size {
        params.set_fp_size(value);
    }
    if let FieldInput::Value(value) = bits_per_feature {
        params.set_bits_per_feature(value);
    }
    if let FieldInput::Value(value) = count_bounds {
        params.set_count_bounds(value);
    }
    if let FieldInput::Value(value) = from_atoms {
        params.from_atoms = value;
    }
    if let FieldInput::Value(value) = ignore_atoms {
        params.ignore_atoms = value;
    }
    if let FieldInput::Value(value) = custom_atom_invariants {
        params.custom_atom_invariants = value;
    }
    if let FieldInput::Value(value) = custom_bond_invariants {
        params.custom_bond_invariants = value;
    }
    if let FieldInput::Value(value) = atom_invariants_generator {
        params.atom_invariants_generator = value;
    }
    if let FieldInput::Value(value) = use_legacy_stereo_perception {
        params.use_legacy_stereo_perception = value;
    }
    params
}

run_fingerprint!(
    run_morgan,
    fingerprint_morgan_with_params,
    MorganFingerprintParams,
    Fingerprint,
    MorganReadError
);
run_fingerprint!(
    run_morgan_count,
    fingerprint_morgan_count_with_params,
    MorganFingerprintParams,
    SparseCountFingerprint32,
    MorganReadError
);
run_fingerprint!(
    run_morgan_sparse,
    fingerprint_morgan_sparse_with_params,
    MorganFingerprintParams,
    SparseBitFingerprint,
    MorganReadError
);
run_fingerprint!(
    run_morgan_sparse_count,
    fingerprint_morgan_sparse_count_with_params,
    MorganFingerprintParams,
    SparseCountFingerprint,
    MorganReadError
);
run_fingerprint!(
    run_atom_pair,
    fingerprint_atom_pair_with_params,
    AtomPairFingerprintParams,
    Fingerprint,
    AtomPairReadError
);
run_fingerprint!(
    run_atom_pair_count,
    fingerprint_atom_pair_count_with_params,
    AtomPairFingerprintParams,
    SparseCountFingerprint32,
    AtomPairReadError
);
run_fingerprint!(
    run_atom_pair_sparse,
    fingerprint_atom_pair_sparse_with_params,
    AtomPairFingerprintParams,
    SparseBitFingerprint,
    AtomPairReadError
);
run_fingerprint!(
    run_atom_pair_sparse_count,
    fingerprint_atom_pair_sparse_count_with_params,
    AtomPairFingerprintParams,
    SparseCountFingerprint,
    AtomPairReadError
);
run_fingerprint!(
    run_topological_torsion,
    fingerprint_topological_torsion_with_params,
    TopologicalTorsionFingerprintParams,
    Fingerprint,
    TopologicalTorsionReadError
);
run_fingerprint!(
    run_topological_torsion_count,
    fingerprint_topological_torsion_count_with_params,
    TopologicalTorsionFingerprintParams,
    SparseCountFingerprint32,
    TopologicalTorsionReadError
);
run_fingerprint!(
    run_topological_torsion_sparse,
    fingerprint_topological_torsion_sparse_with_params,
    TopologicalTorsionFingerprintParams,
    SparseBitFingerprint,
    TopologicalTorsionReadError
);
run_fingerprint!(
    run_topological_torsion_sparse_count,
    fingerprint_topological_torsion_sparse_count_with_params,
    TopologicalTorsionFingerprintParams,
    SparseCountFingerprint,
    TopologicalTorsionReadError
);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fingerprint_providers::FingerprintAdditionalOutput;
    use crate::fingerprint_values::{
        SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32,
    };
    use crate::native::{Fingerprint, Molecule};
    use pyo3::exceptions::PyTypeError;
    use pyo3::prelude::*;
    use pyo3::types::PyDict;

    #[test]
    fn ordinary_fingerprint_families_forward_all_outputs_options_errors_and_metadata() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Molecule::from_smiles("CCCO").unwrap();
            let before = source.to_binary().unwrap();
            let molecule = Py::new(py, Molecule { inner: source.clone() })?;
            macro_rules! family {
                ($name:ident, $params:ident, $generator:ident, $failure:ident) => {{
                    let mut options = ck::$params::default();
                    options.generator.fp_size = 256;
                    options.from_atoms = Some(vec![0, 2]);
                    options.ignore_atoms = Some(vec![]);
                    options.custom_atom_invariants = Some(vec![11, 12, 13, 14]);
                    options.custom_bond_invariants = Some(vec![21, 22, 23]);
                    let params = Py::new(py, $params::from_rust(py, options.clone())?)?;
                    let kwargs = PyDict::new(py);
                    kwargs.set_item("fp_size", 256_u32)?;
                    kwargs.set_item("from_atoms", vec![0_u32, 2])?;
                    kwargs.set_item("ignore_atoms", Vec::<u32>::new())?;
                    kwargs.set_item("custom_atom_invariants", vec![11_u32, 12, 13, 14])?;
                    kwargs.set_item("custom_bond_invariants", vec![21_u32, 22, 23])?;
                    macro_rules! output {
                        ($suffix:tt, $value:ident) => { paste::paste! {
                            let method = stringify!([<fingerprint_ $name $suffix>]);
                            let alias = stringify!([<fingerprint_ $name $suffix _with_params>]);
                            let actual = molecule.bind(py).call_method0(method)?
                                .extract::<Py<$value>>()?;
                            assert_eq!(actual.borrow(py).inner, source.[<fingerprint_ $name $suffix>]().unwrap());
                            let expected = source.[<fingerprint_ $name $suffix _with_params>](&options, None).unwrap();
                            let configured = molecule.bind(py).call_method1(method, (&params,))?
                                .extract::<Py<$value>>()?;
                            assert_eq!(configured.borrow(py).inner, expected);
                            let direct = molecule.bind(py).call_method1(alias, (&params,))?
                                .extract::<Py<$value>>()?;
                            assert_eq!(direct.borrow(py).inner, expected);
                            let keywords = molecule.bind(py).call_method(method, (), Some(&kwargs))?
                                .extract::<Py<$value>>()?;
                            assert_eq!(keywords.borrow(py).inner, expected);
                            assert!(molecule.bind(py).call_method(method, (&params,), Some(&kwargs))
                                .unwrap_err().is_instance_of::<PyTypeError>(py));

                            let mut native_output = ck::FingerprintAdditionalOutput::new();
                            native_output.allocate_atom_counts();
                            native_output.allocate_bit_info_map();
                            native_output.allocate_bit_paths();
                            source.[<fingerprint_ $name $suffix _with_params>](&options, Some(&mut native_output)).unwrap();
                            let mut output = ck::FingerprintAdditionalOutput::new();
                            output.allocate_atom_counts();
                            output.allocate_bit_info_map();
                            output.allocate_bit_paths();
                            let output = Py::new(py, FingerprintAdditionalOutput { inner: output })?;
                            molecule.bind(py).call_method1(method, (&params, &output))?;
                            let output = output.borrow(py);
                            assert_eq!(output.inner.atom_counts(), native_output.atom_counts());
                            assert_eq!(output.inner.bit_info_map(), native_output.bit_info_map());
                            assert_eq!(output.inner.bit_paths(), native_output.bit_paths());

                            let mut invalid = options.clone();
                            invalid.custom_atom_invariants = Some(vec![]);
                            let native_error = source.[<fingerprint_ $name $suffix _with_params>](&invalid, None).err().unwrap();
                            let invalid = Py::new(py, $params::from_rust(py, invalid)?)?;
                            let error = molecule.bind(py).call_method1(method, (&invalid,)).unwrap_err();
                            assert!(error.is_instance(py, <crate::morgan_errors::$failure as crate::error_framework::ExceptionPayload>::exception_type(py)?.as_any()));
                            assert_eq!(error.value(py).str()?.to_str()?, native_error.to_string());
                            let cause = error.value(py).getattr("source")?;
                            assert!(cause.is(&error.value(py).getattr("__cause__")?));
                        } };
                    }
                    output!("", Fingerprint);
                    output!(_count, SparseCountFingerprint32);
                    output!(_sparse, SparseBitFingerprint);
                    output!(_sparse_count, SparseCountFingerprint);
                    assert_eq!(params.borrow(py).to_rust(py)?, options);
                }};
            }
            family!(morgan, MorganFingerprintParams, MorganGeneratorParams, MorganReadError);
            family!(atom_pair, AtomPairFingerprintParams, AtomPairGeneratorParams, AtomPairReadError);
            family!(topological_torsion, TopologicalTorsionFingerprintParams, TopologicalTorsionGeneratorParams, TopologicalTorsionReadError);
            assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
            Ok(())
        }).unwrap();
    }
}
