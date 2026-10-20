//! Reusable Morgan state: each method forwards to the canonical live facade.
use crate::fingerprint_configuration::MorganGeneratorParams;
use crate::fingerprint_providers::{FingerprintAdditionalOutput, MorganBondInvariantsGenerator};
use crate::fingerprint_values::{
    SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32,
};
use crate::morgan_configuration::MorganAtomInvariantsGenerator;
use crate::morgan_errors::MorganReadError;
use crate::native::{Fingerprint, Molecule};
use crate::property_values::decode_source_text;
use cosmolkit as ck;
use pyo3::prelude::*;

macro_rules! molecule_call {
    ($name:ident, $result:ident) => {
        pub(crate) fn $name(
            molecule: &Molecule,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            params: Option<&crate::params::MorganCallParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<$result> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            let additional_output = additional_output.map(|value| &mut value.inner);
            py.detach(|| {
                molecule
                    .inner
                    .$name(&generator.inner, params.as_ref(), additional_output)
            })
            .map(|inner| $result { inner })
            .map_err(MorganReadError::from_rust)
        }
    };
}
molecule_call!(fingerprint_morgan_with_generator, Fingerprint);
molecule_call!(
    fingerprint_morgan_count_with_generator,
    SparseCountFingerprint32
);
molecule_call!(
    fingerprint_morgan_sparse_with_generator,
    SparseBitFingerprint
);
molecule_call!(
    fingerprint_morgan_sparse_count_with_generator,
    SparseCountFingerprint
);

pub(crate) fn call_keywords(
    py: Python<'_>,
    from_atoms: Option<Vec<u32>>,
    ignore_atoms: Option<Vec<u32>>,
    custom_atom_invariants: Option<Vec<u32>>,
    custom_bond_invariants: Option<Vec<u32>>,
    conformer_id: Option<i32>,
) -> PyResult<crate::params::MorganCallParams> {
    crate::params::MorganCallParams::from_rust(
        py,
        ck::MorganCallParams::new(
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            conformer_id.unwrap_or_else(|| ck::MorganCallParams::default().conformer_id),
        ),
    )
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct MorganFingerprintGenerator {
    pub(crate) inner: ck::MorganFingerprintGenerator,
}

#[pymethods]
impl MorganFingerprintGenerator {
    #[new]
    #[pyo3(signature = (params=None, atom_invariants=None, bond_invariants=None))]
    fn construct(
        py: Python<'_>,
        params: Option<&MorganGeneratorParams>,
        atom_invariants: Option<&MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&MorganBondInvariantsGenerator>,
    ) -> PyResult<Self> {
        let bond_invariants =
            pyo3::Python::attach(|py| bond_invariants.map(|value| value.to_rust(py)).transpose())?;

        let params = params.map(|value| value.to_rust(py)).transpose()?;
        ck::MorganFingerprintGenerator::new(
            params.as_ref(),
            atom_invariants.map(|value| &value.inner),
            bond_invariants.as_ref(),
        )
        .map(|inner| Self { inner })
        .map_err(MorganReadError::from_rust)
    }
    #[staticmethod]
    #[pyo3(signature = (params=None, atom_invariants=None, bond_invariants=None))]
    fn new(
        py: Python<'_>,
        params: Option<&MorganGeneratorParams>,
        atom_invariants: Option<&MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&MorganBondInvariantsGenerator>,
    ) -> PyResult<Self> {
        Self::construct(py, params, atom_invariants, bond_invariants)
    }
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        ck::MorganFingerprintGenerator::from_json(json)
            .map(|inner| Self { inner })
            .map_err(MorganReadError::from_rust)
    }
    #[getter]
    fn settings(&self) -> MorganSettings {
        MorganSettings {
            inner: self.inner.settings(),
        }
    }
    fn info_string(&self) -> PyResult<String> {
        self.inner.info_string().map_err(MorganReadError::from_rust)
    }
    fn to_json(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_json()
            .map_err(MorganReadError::from_rust)
            .and_then(|text| decode_source_text(py, &text))
    }
    #[pyo3(signature = (molecules, num_threads=1))]
    fn fingerprints(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<Molecule>>>,
        num_threads: i32,
    ) -> PyResult<crate::sequence::Output<ck::FingerprintBatch>> {
        let guards = molecules
            .iter()
            .map(|value| value.as_ref().map(|value| value.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let inputs = guards
            .iter()
            .map(|value| value.as_ref().map(|value| &value.inner))
            .collect::<Vec<_>>();
        py.detach(|| self.inner.fingerprints(&inputs, num_threads))
            .map(crate::sequence::Output)
            .map_err(MorganReadError::from_rust)
    }
    #[pyo3(signature = (molecules, num_threads=1))]
    fn counts(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<SparseCountFingerprint32>>> {
        let guards = molecules
            .iter()
            .map(|value| value.as_ref().map(|value| value.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let inputs = guards
            .iter()
            .map(|value| value.as_ref().map(|value| &value.inner))
            .collect::<Vec<_>>();
        py.detach(|| self.inner.counts(&inputs, num_threads))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| SparseCountFingerprint32 { inner }))
                    .collect()
            })
            .map_err(MorganReadError::from_rust)
    }
    #[pyo3(signature = (molecules, num_threads=1))]
    fn sparse_fingerprints(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<SparseBitFingerprint>>> {
        let guards = molecules
            .iter()
            .map(|value| value.as_ref().map(|value| value.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let inputs = guards
            .iter()
            .map(|value| value.as_ref().map(|value| &value.inner))
            .collect::<Vec<_>>();
        py.detach(|| self.inner.sparse_fingerprints(&inputs, num_threads))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| SparseBitFingerprint { inner }))
                    .collect()
            })
            .map_err(MorganReadError::from_rust)
    }
    #[pyo3(signature = (molecules, num_threads=1))]
    fn sparse_counts(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<SparseCountFingerprint>>> {
        let guards = molecules
            .iter()
            .map(|value| value.as_ref().map(|value| value.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let inputs = guards
            .iter()
            .map(|value| value.as_ref().map(|value| &value.inner))
            .collect::<Vec<_>>();
        py.detach(|| self.inner.sparse_counts(&inputs, num_threads))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| SparseCountFingerprint { inner }))
                    .collect()
            })
            .map_err(MorganReadError::from_rust)
    }
}

#[pyclass(module = "cosmolkit")]
pub struct MorganSettings {
    inner: ck::MorganSettings,
}

macro_rules! live_settings {
    ($($field:ident: $ty:ty),+ $(,)?) => {
        paste::paste! {
            #[pymethods]
            impl MorganSettings {
                $(#[getter]
                fn $field(&self) -> PyResult<$ty> {
                    self.inner.$field().map_err(MorganReadError::from_rust)
                }
                #[setter($field)]
                fn [<assign_ $field>](&mut self, value: $ty) -> PyResult<()> {
                    self.[<set_ $field>](value)
                }
                fn [<set_ $field>](&mut self, value: $ty) -> PyResult<()> {
                    self.inner.[<set_ $field>](value).map_err(MorganReadError::from_rust)
                })+
                #[getter]
                fn params(&self, py: Python<'_>) -> PyResult<MorganGeneratorParams> {
                    let value = self.inner.params().map_err(MorganReadError::from_rust)?;
                    MorganGeneratorParams::from_rust(py, value)
                }
            }
        }
    };
}
live_settings! {
    radius: u32, only_nonzero_invariants: bool, include_redundant_environments: bool,
    include_chirality: bool, count_simulation: bool, fp_size: u32,
    bits_per_feature: u32, count_bounds: Vec<u32>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fingerprint_json::SerdeJsonError;
    use crate::morgan_errors::{FingerprintPreparationError, MorganError};
    use crate::query_compiler::QueryGraph;
    use crate::smarts_writer::SmartsWriteError;
    use pyo3::exceptions::PyTypeError;
    use pyo3::types::PyDict;

    #[test]
    fn reusable_morgan_forwards_live_settings_scalar_bulk_outputs_and_retained_error_graphs() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Molecule::from_smiles("C[C@H](O)CC").unwrap();
            let before = source.to_binary().unwrap();
            let molecule = Py::new(
                py,
                Molecule {
                    inner: source.clone(),
                },
            )?;
            let mut native_params = ck::MorganGeneratorParams::default();
            native_params.fp_size = 128;
            native_params.include_chirality = true;
            let params = MorganGeneratorParams::from_rust(py, native_params.clone())?;
            let generator = Py::new(
                py,
                MorganFingerprintGenerator::new(py, Some(&params), None, None)?,
            )?;
            let settings = generator.bind(py).call_method0("settings")?;
            settings.setattr("radius", 3)?;
            assert_eq!(generator.borrow(py).inner.settings().radius().unwrap(), 3);
            assert_eq!(params.to_rust(py)?, native_params);
            let json = generator
                .bind(py)
                .call_method0("to_json")?
                .extract::<String>()?;
            let recreated = MorganFingerprintGenerator::from_json(&json)?;
            assert_eq!(
                recreated.inner.settings().params().unwrap(),
                generator.borrow(py).inner.settings().params().unwrap()
            );
            let call = call_keywords(py, Some(vec![0, 2]), None, None, None, None)?;
            let call = Py::new(py, call)?;
            let output = Py::new(
                py,
                FingerprintAdditionalOutput {
                    inner: ck::FingerprintAdditionalOutput::default(),
                },
            )?;
            output.bind(py).call_method0("allocate_atom_counts")?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("from_atoms", vec![0_u32, 2])?;
            for (method, bulk) in [
                ("fingerprint_morgan_with_generator", "fingerprints"),
                ("fingerprint_morgan_count_with_generator", "counts"),
                (
                    "fingerprint_morgan_sparse_with_generator",
                    "sparse_fingerprints",
                ),
                (
                    "fingerprint_morgan_sparse_count_with_generator",
                    "sparse_counts",
                ),
            ] {
                let actual = molecule
                    .bind(py)
                    .call_method1(method, (&generator, &call, &output))?;
                let keywords =
                    molecule
                        .bind(py)
                        .call_method(method, (&generator,), Some(&kwargs))?;
                let defaults = molecule.bind(py).call_method1(method, (&generator,))?;
                match method {
                    "fingerprint_morgan_with_generator" => {
                        let expected = source
                            .fingerprint_morgan_with_generator(
                                &generator.borrow(py).inner,
                                Some(&call.borrow(py).to_rust(py)?),
                                None,
                            )
                            .unwrap();
                        assert_eq!(
                            actual.extract::<Py<Fingerprint>>()?.borrow(py).inner,
                            expected
                        );
                        assert_eq!(
                            keywords.extract::<Py<Fingerprint>>()?.borrow(py).inner,
                            expected
                        );
                    }
                    "fingerprint_morgan_count_with_generator" => {
                        let expected = source
                            .fingerprint_morgan_count_with_generator(
                                &generator.borrow(py).inner,
                                Some(&call.borrow(py).to_rust(py)?),
                                None,
                            )
                            .unwrap();
                        assert_eq!(
                            actual
                                .extract::<Py<SparseCountFingerprint32>>()?
                                .borrow(py)
                                .inner,
                            expected
                        );
                        assert_eq!(
                            keywords
                                .extract::<Py<SparseCountFingerprint32>>()?
                                .borrow(py)
                                .inner,
                            expected
                        );
                    }
                    "fingerprint_morgan_sparse_with_generator" => {
                        let expected = source
                            .fingerprint_morgan_sparse_with_generator(
                                &generator.borrow(py).inner,
                                Some(&call.borrow(py).to_rust(py)?),
                                None,
                            )
                            .unwrap();
                        assert_eq!(
                            actual
                                .extract::<Py<SparseBitFingerprint>>()?
                                .borrow(py)
                                .inner,
                            expected
                        );
                        assert_eq!(
                            keywords
                                .extract::<Py<SparseBitFingerprint>>()?
                                .borrow(py)
                                .inner,
                            expected
                        );
                    }
                    _ => {
                        let expected = source
                            .fingerprint_morgan_sparse_count_with_generator(
                                &generator.borrow(py).inner,
                                Some(&call.borrow(py).to_rust(py)?),
                                None,
                            )
                            .unwrap();
                        assert_eq!(
                            actual
                                .extract::<Py<SparseCountFingerprint>>()?
                                .borrow(py)
                                .inner,
                            expected
                        );
                        assert_eq!(
                            keywords
                                .extract::<Py<SparseCountFingerprint>>()?
                                .borrow(py)
                                .inner,
                            expected
                        );
                    }
                }
                for threads in [1, 2] {
                    let batch = generator.bind(py).call_method1(
                        bulk,
                        (
                            vec![
                                Some(molecule.clone_ref(py)),
                                None,
                                Some(molecule.clone_ref(py)),
                            ],
                            threads,
                        ),
                    )?;
                    assert_eq!(batch.len()?, 3);
                    assert!(batch.get_item(1)?.is_none());
                    // Compare the exposed native values, not their repr or a hash.
                    let first = batch.get_item(0)?;
                    match bulk {
                        "fingerprints" => assert_eq!(
                            first.extract::<Py<Fingerprint>>()?.borrow(py).inner,
                            defaults.extract::<Py<Fingerprint>>()?.borrow(py).inner
                        ),
                        "counts" => assert_eq!(
                            first
                                .extract::<Py<SparseCountFingerprint32>>()?
                                .borrow(py)
                                .inner,
                            defaults
                                .extract::<Py<SparseCountFingerprint32>>()?
                                .borrow(py)
                                .inner
                        ),
                        "sparse_fingerprints" => assert_eq!(
                            first
                                .extract::<Py<SparseBitFingerprint>>()?
                                .borrow(py)
                                .inner,
                            defaults
                                .extract::<Py<SparseBitFingerprint>>()?
                                .borrow(py)
                                .inner
                        ),
                        _ => assert_eq!(
                            first
                                .extract::<Py<SparseCountFingerprint>>()?
                                .borrow(py)
                                .inner,
                            defaults
                                .extract::<Py<SparseCountFingerprint>>()?
                                .borrow(py)
                                .inner
                        ),
                    }
                }
                assert!(
                    molecule
                        .bind(py)
                        .call_method(method, (&generator, &call), Some(&kwargs))
                        .unwrap_err()
                        .is_instance_of::<PyTypeError>(py)
                );
            }
            assert_eq!(
                output.bind(py).call_method0("atom_counts")?.len()?,
                source.num_atoms()
            );
            assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
            drop(generator);
            settings.call_method1("set_radius", (4,))?;
            assert_eq!(settings.getattr("radius")?.extract::<u32>()?, 4);

            let error = MorganFingerprintGenerator::from_json("{bad").err().unwrap();
            assert!(error.is_instance_of::<MorganReadError>(py));
            let generator_error = error.value(py).getattr("source")?;
            assert!(generator_error.is_instance_of::<MorganError>());
            let json = generator_error.getattr("source")?;
            let parse = json.getattr("source")?;
            assert!(parse.is_instance_of::<SerdeJsonError>());
            assert!(parse.getattr("column")?.extract::<usize>()? > 0);
            assert!(parse.is(&json.getattr("__cause__")?));
            // sanitize=false does not imply a missing valence cache: source
            // parsing/cleanup can still prepare it. Exercise this exact native
            // error value without inventing a different SMILES parse policy.
            let native_error = ck::MorganReadError::Preparation(
                ck::FingerprintPreparationError::MissingPreparedValence,
            );
            let display = native_error.to_string();
            let error = MorganReadError::from_rust(native_error);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is_instance_of::<FingerprintPreparationError>()
            );

            let io = std::io::Error::from_raw_os_error(13);
            let raw = io.raw_os_error();
            let error = MorganReadError::from_rust(ck::MorganReadError::Generator(
                ck::MorganError::Worker(ck::FingerprintWorkerError::ThreadSpawn(io)),
            ));
            let io = error
                .value(py)
                .getattr("source")?
                .getattr("source")?
                .getattr("source")?;
            assert!(io.is_instance_of::<crate::io_source_errors::IoError>());
            assert_eq!(io.getattr("raw_os_error")?.extract::<Option<i32>>()?, raw);

            let mut recursive =
                ck::RecursiveStructureQuery::from_query_graph(ck::parse_smarts("CC").unwrap(), 41);
            recursive.set_initialized(true);
            let error = MorganReadError::from_rust(ck::MorganReadError::Generator(
                ck::MorganError::SmartsWrite(ck::SmartsWriteError::UnsupportedAtomQuery {
                    predicate: ck::AtomQueryPredicate::RecursiveSmarts(recursive),
                }),
            ));
            let writer = error.value(py).getattr("source")?.getattr("source")?;
            assert!(writer.is_instance_of::<SmartsWriteError>());
            let recursive = writer.getattr("predicate")?.getattr("value")?;
            assert!(recursive.getattr("initialized")?.extract::<bool>()?);
            let query = recursive
                .getattr("query_graph")?
                .extract::<Py<QueryGraph>>()?;
            drop(error);
            drop(writer);
            drop(recursive);
            assert_eq!(
                query
                    .bind(py)
                    .call_method0("num_atoms")?
                    .extract::<usize>()?,
                2
            );
            Ok(())
        })
        .unwrap();
    }
}
