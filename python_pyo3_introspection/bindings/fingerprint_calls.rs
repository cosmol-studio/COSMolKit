//! Ordinary fingerprint calls: construct native options; never regenerate algorithms.
use crate::fingerprint_configuration::{
    AtomPairFingerprintParams, AtomPairParams, MorganParams, TopologicalTorsionFingerprintParams,
    TopologicalTorsionParams,
};
use crate::fingerprint_providers::AtomPairAtomInvariantsGenerator;
use crate::morgan_configuration::{MorganFingerprintParams, MorganInvariants};
use cosmolkit as ck;

pub(crate) fn morgan_keywords(
    generator: Option<&MorganParams>,
    from_atoms: Option<Vec<u32>>,
    ignore_atoms: Option<Vec<u32>>,
    custom_atom_invariants: Option<Vec<u32>>,
    custom_bond_invariants: Option<Vec<u32>>,
    conformer_id: Option<i32>,
    invariants: Option<&MorganInvariants>,
) -> MorganFingerprintParams {
    let defaults = ck::MorganFingerprintParams::default();
    MorganFingerprintParams {
        inner: ck::MorganFingerprintParams {
            generator: generator
                .map(|value| value.inner.clone())
                .unwrap_or(defaults.generator),
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            conformer_id: conformer_id.unwrap_or(defaults.conformer_id),
            invariants: invariants
                .map(|value| value.inner.clone())
                .unwrap_or(defaults.invariants),
        },
    }
}

pub(crate) fn atom_pair_keywords(
    generator: Option<&AtomPairParams>,
    from_atoms: Option<Vec<u32>>,
    ignore_atoms: Option<Vec<u32>>,
    custom_atom_invariants: Option<Vec<u32>>,
    custom_bond_invariants: Option<Vec<u32>>,
    conformer_id: Option<i32>,
    atom_invariants_generator: Option<&AtomPairAtomInvariantsGenerator>,
    use_legacy_stereo_perception: Option<bool>,
) -> AtomPairFingerprintParams {
    let defaults = ck::AtomPairFingerprintParams::default();
    AtomPairFingerprintParams {
        inner: ck::AtomPairFingerprintParams {
            generator: generator
                .map(|value| value.inner.clone())
                .unwrap_or(defaults.generator),
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            conformer_id: conformer_id.unwrap_or(defaults.conformer_id),
            atom_invariants_generator: atom_invariants_generator.map(|value| value.inner),
            use_legacy_stereo_perception: use_legacy_stereo_perception
                .unwrap_or(defaults.use_legacy_stereo_perception),
        },
    }
}

pub(crate) fn topological_torsion_keywords(
    generator: Option<&TopologicalTorsionParams>,
    from_atoms: Option<Vec<u32>>,
    ignore_atoms: Option<Vec<u32>>,
    custom_atom_invariants: Option<Vec<u32>>,
    custom_bond_invariants: Option<Vec<u32>>,
    conformer_id: Option<i32>,
    atom_invariants_generator: Option<&AtomPairAtomInvariantsGenerator>,
    use_legacy_stereo_perception: Option<bool>,
) -> TopologicalTorsionFingerprintParams {
    let defaults = ck::TopologicalTorsionFingerprintParams::default();
    TopologicalTorsionFingerprintParams {
        inner: ck::TopologicalTorsionFingerprintParams {
            generator: generator
                .map(|value| value.inner.clone())
                .unwrap_or(defaults.generator),
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            conformer_id: conformer_id.unwrap_or(defaults.conformer_id),
            atom_invariants_generator: atom_invariants_generator.map(|value| value.inner),
            use_legacy_stereo_perception: use_legacy_stereo_perception
                .unwrap_or(defaults.use_legacy_stereo_perception),
        },
    }
}

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
                    let generator = Py::new(py, $generator { inner: options.generator.clone() })?;
                    let params = Py::new(py, $params { inner: options.clone() })?;
                    let kwargs = PyDict::new(py);
                    kwargs.set_item("generator", &generator)?;
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
                            let invalid = Py::new(py, $params { inner: invalid })?;
                            let error = molecule.bind(py).call_method1(method, (&invalid,)).unwrap_err();
                            assert!(error.is_instance_of::<crate::morgan_errors::$failure>(py));
                            assert_eq!(error.value(py).str()?.to_str()?, native_error.to_string());
                            let cause = error.value(py).getattr("source")?;
                            assert!(cause.is(&error.value(py).getattr("__cause__")?));
                        } };
                    }
                    output!("", Fingerprint);
                    output!(_count, SparseCountFingerprint32);
                    output!(_sparse, SparseBitFingerprint);
                    output!(_sparse_count, SparseCountFingerprint);
                    assert_eq!(params.borrow(py).inner, options);
                }};
            }
            family!(morgan, MorganFingerprintParams, MorganParams, MorganReadError);
            family!(atom_pair, AtomPairFingerprintParams, AtomPairParams, AtomPairReadError);
            family!(topological_torsion, TopologicalTorsionFingerprintParams, TopologicalTorsionParams, TopologicalTorsionReadError);
            assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
            Ok(())
        }).unwrap();
    }
}
