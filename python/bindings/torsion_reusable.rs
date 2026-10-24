//! Reusable torsion operators forward the sole native generator and live settings.
use crate::fingerprint_configuration::TopologicalTorsionGeneratorParams;
use crate::fingerprint_providers::{AtomPairAtomInvariantsGenerator, FingerprintAdditionalOutput};
use crate::fingerprint_values::{
    SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32,
};
use crate::morgan_errors::TopologicalTorsionReadError;
use crate::native::{Fingerprint, Molecule};
use crate::params::TopologicalTorsionCallParams;
use cosmolkit as ck;
use pyo3::prelude::*;

pub(crate) fn call_keywords(
    py: Python<'_>,
    from_atoms: Option<Vec<u32>>,
    ignore_atoms: Option<Vec<u32>>,
    custom_atom_invariants: Option<Vec<u32>>,
    custom_bond_invariants: Option<Vec<u32>>,
    use_legacy_stereo_perception: Option<bool>,
) -> PyResult<TopologicalTorsionCallParams> {
    let defaults = ck::TopologicalTorsionCallParams::default();
    TopologicalTorsionCallParams::from_rust(
        py,
        ck::TopologicalTorsionCallParams {
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            use_legacy_stereo_perception: use_legacy_stereo_perception
                .unwrap_or(defaults.use_legacy_stereo_perception),
        },
    )
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct TopologicalTorsionFingerprintGenerator {
    pub(crate) inner: ck::TopologicalTorsionFingerprintGenerator,
}
#[pymethods]
impl TopologicalTorsionFingerprintGenerator {
    #[new]
    #[pyo3(signature=(params=None,atom_invariants=None))]
    fn construct(
        py: Python<'_>,
        params: Option<&TopologicalTorsionGeneratorParams>,
        atom_invariants: Option<&AtomPairAtomInvariantsGenerator>,
    ) -> PyResult<Self> {
        let atom_invariants =
            pyo3::Python::attach(|py| atom_invariants.map(|value| value.to_rust(py)).transpose())?;

        let params = params.map(|value| value.to_rust(py)).transpose()?;
        ck::TopologicalTorsionFingerprintGenerator::new(params.as_ref(), atom_invariants)
            .map(|inner| Self { inner })
            .map_err(TopologicalTorsionReadError::from_rust)
    }
    #[staticmethod]
    #[pyo3(signature=(params=None,atom_invariants=None))]
    fn new(
        py: Python<'_>,
        params: Option<&TopologicalTorsionGeneratorParams>,
        atom_invariants: Option<&AtomPairAtomInvariantsGenerator>,
    ) -> PyResult<Self> {
        Self::construct(py, params, atom_invariants)
    }
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        ck::TopologicalTorsionFingerprintGenerator::from_json(json)
            .map(|inner| Self { inner })
            .map_err(TopologicalTorsionReadError::from_rust)
    }
    #[getter]
    fn settings(&self) -> TopologicalTorsionSettings {
        TopologicalTorsionSettings {
            inner: self.inner.settings(),
        }
    }
    fn info_string(&self) -> PyResult<String> {
        self.inner
            .info_string()
            .map_err(TopologicalTorsionReadError::from_rust)
    }
    fn to_json(&self) -> PyResult<String> {
        self.inner
            .to_json()
            .map_err(TopologicalTorsionReadError::from_rust)
    }

    #[pyo3(signature=(molecules,num_threads=1))]
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
            .map_err(TopologicalTorsionReadError::from_rust)
    }

    #[pyo3(signature=(molecules,num_threads=1))]
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
            .map_err(TopologicalTorsionReadError::from_rust)
    }

    #[pyo3(signature=(molecules,num_threads=1))]
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
            .map_err(TopologicalTorsionReadError::from_rust)
    }

    #[pyo3(signature=(molecules,num_threads=1))]
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
            .map_err(TopologicalTorsionReadError::from_rust)
    }
}
#[pyclass(module = "cosmolkit")]
pub struct TopologicalTorsionSettings {
    inner: ck::TopologicalTorsionSettings,
}
macro_rules! live_settings {
    ($($field:ident:$ty:ty),+ $(,)?)=>{
        paste::paste!{
            #[pymethods]
            impl TopologicalTorsionSettings {
                $(#[getter]
                fn $field(&self)->PyResult<$ty>{self.inner.$field().map_err(TopologicalTorsionReadError::from_rust)}
                #[setter($field)]
                fn [<assign_ $field>](&mut self,value:$ty)->PyResult<()>{self.[<set_ $field>](value)}
                fn [<set_ $field>](&mut self,value:$ty)->PyResult<()>{
                    self.inner.[<set_ $field>](value).map_err(TopologicalTorsionReadError::from_rust)
                })+
                #[getter]
                fn params(&self, py: Python<'_>)->PyResult<TopologicalTorsionGeneratorParams>{
                    let value = self.inner.params().map_err(TopologicalTorsionReadError::from_rust)?;
                    TopologicalTorsionGeneratorParams::from_rust(py, value)
                }
            }
        }
    };
}
live_settings! {
    torsion_atom_count:u32,only_shortest_paths:bool,include_chirality:bool,
    count_simulation:bool,fp_size:u32,bits_per_feature:u32,count_bounds:Vec<u32>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fingerprint_json::SerdeJsonError;
    use crate::morgan_errors::TopologicalTorsionError;
    use pyo3::exceptions::PyTypeError;
    use pyo3::types::PyDict;

    #[test]
    fn reusable_torsion_forwards_live_settings_all_scalar_bulk_forms_and_real_errors() {
        Python::initialize();
        Python::attach(|py|->PyResult<()>{
            let source=ck::Molecule::from_smiles("C[C@H](O)CC").unwrap();
            let before=source.to_binary().unwrap();
            let molecule=Py::new(py,Molecule{inner:source.clone()})?;
            let native_params=ck::TopologicalTorsionGeneratorParams{fp_size:128,include_chirality:true,..Default::default()};
            let native_provider=ck::AtomPairAtomInvariantsGenerator{include_chirality:true,topological_torsion_correction:true};
            let params=Py::new(py,TopologicalTorsionGeneratorParams::from_rust(py,native_params.clone())?)?;
            let provider=Py::new(py,AtomPairAtomInvariantsGenerator{inner:crate::parameters::NativeField::Owned(native_provider)})?;
            let generator=py.get_type::<TopologicalTorsionFingerprintGenerator>()
                .call1((&params,&provider))?.extract::<Py<TopologicalTorsionFingerprintGenerator>>()?;
            let native=ck::TopologicalTorsionFingerprintGenerator::new(Some(&native_params),Some(native_provider)).unwrap();
            let settings=generator.bind(py).call_method0("settings")?;
            macro_rules! setting {
                ($field:ident,$value:expr)=>{paste::paste!{
                    settings.setattr(stringify!($field),$value)?;
                    native.settings().[<set_ $field>]($value).unwrap();
                    assert!(settings.getattr(stringify!($field))?.eq(native.settings().$field().unwrap())?);
                    settings.call_method1(stringify!([<set_ $field>]),($value,))?;
                }};
            }
            setting!(torsion_atom_count,3_u32);
            setting!(only_shortest_paths,true);
            setting!(include_chirality,false);
            setting!(count_simulation,false);
            setting!(fp_size,256_u32);
            setting!(bits_per_feature,2_u32);
            setting!(count_bounds,vec![1_u32,3,5]);
            let snapshot=settings.call_method0("params")?.extract::<Py<TopologicalTorsionGeneratorParams>>()?;
            assert_eq!(snapshot.borrow(py).to_rust(py)?,native.settings().params().unwrap());
            assert_eq!(params.borrow(py).to_rust(py)?,native_params);
            assert!(settings.setattr("fp_size","bad").unwrap_err().is_instance_of::<PyTypeError>(py));
            assert_eq!(settings.getattr("fp_size")?.extract::<u32>()?,256);
            assert_eq!(generator.bind(py).call_method0("info_string")?.extract::<String>()?,native.info_string().unwrap());
            let json=generator.bind(py).call_method0("to_json")?.extract::<String>()?;
            assert_eq!(json,native.to_json().unwrap());
            let restored=py.get_type::<TopologicalTorsionFingerprintGenerator>()
                .call_method1("from_json",(json,))?.extract::<Py<TopologicalTorsionFingerprintGenerator>>()?;
            assert_eq!(restored.borrow(py).inner.to_json().unwrap(),native.to_json().unwrap());

            // The explicit provider still includes chirality. Turning only the
            // generator flag off narrows its sparse code space; native Rust
            // rejects the resulting code instead of silently replacing the
            // provider. Preserve that real failure, then select a coherent
            // configuration for the four successful output forms below.
            let expected=source.fingerprint_topological_torsion_sparse_count_with_generator(&native,None,None).unwrap_err();
            assert!(matches!(&expected,ck::TopologicalTorsionReadError::Generator(ck::TopologicalTorsionError::Fingerprint(ck::FingerprintError::SparseIndexOutOfRange{..}))));
            let error=molecule.bind(py).call_method1("fingerprint_topological_torsion_sparse_count_with_generator",(&generator,)).unwrap_err();
            assert!(error.is_instance(py, <TopologicalTorsionReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(error.value(py).str()?.to_str()?,expected.to_string());
            settings.setattr("include_chirality",true)?;
            native.settings().set_include_chirality(true).unwrap();

            let options=ck::TopologicalTorsionCallParams{
                from_atoms:Some(vec![0,4]),ignore_atoms:Some(vec![]),
                custom_atom_invariants:Some(vec![11,12,13,14,15]),
                custom_bond_invariants:Some(vec![21,22,23,24]),..Default::default()
            };
            let call=Py::new(py,TopologicalTorsionCallParams::from_rust(py,options.clone())?)?;
            let kwargs=PyDict::new(py);
            kwargs.set_item("from_atoms",vec![0_u32,4])?;
            kwargs.set_item("ignore_atoms",Vec::<u32>::new())?;
            kwargs.set_item("custom_atom_invariants",vec![11_u32,12,13,14,15])?;
            kwargs.set_item("custom_bond_invariants",vec![21_u32,22,23,24])?;
            macro_rules! output {
                ($suffix:tt,$bulk:ident,$value:ident)=>{paste::paste!{
                    let method=stringify!([<fingerprint_topological_torsion $suffix _with_generator>]);
                    let expected_default=source.[<fingerprint_topological_torsion $suffix _with_generator>](&native,None,None).unwrap();
                    let default=molecule.bind(py).call_method1(method,(&generator,))?.extract::<Py<$value>>()?;
                    assert_eq!(default.borrow(py).inner,expected_default);
                    let expected=source.[<fingerprint_topological_torsion $suffix _with_generator>](&native,Some(&options),None).unwrap();
                    let configured=molecule.bind(py).call_method1(method,(&generator,&call))?.extract::<Py<$value>>()?;
                    assert_eq!(configured.borrow(py).inner,expected);
                    let keywords=molecule.bind(py).call_method(method,(&generator,),Some(&kwargs))?.extract::<Py<$value>>()?;
                    assert_eq!(keywords.borrow(py).inner,expected);
                    assert!(molecule.bind(py).call_method(method,(&generator,&call),Some(&kwargs)).unwrap_err().is_instance_of::<PyTypeError>(py));

                    let mut native_output=ck::FingerprintAdditionalOutput::new();
                    native_output.allocate_atom_counts();native_output.allocate_bit_paths();
                    source.[<fingerprint_topological_torsion $suffix _with_generator>](&native,Some(&options),Some(&mut native_output)).unwrap();
                    let mut output=ck::FingerprintAdditionalOutput::new();
                    output.allocate_atom_counts();output.allocate_bit_paths();
                    let output=Py::new(py,FingerprintAdditionalOutput{inner:output})?;
                    molecule.bind(py).call_method1(method,(&generator,&call,&output))?;
                    assert_eq!(output.borrow(py).inner.atom_counts(),native_output.atom_counts());
                    assert_eq!(output.borrow(py).inner.bit_paths(),native_output.bit_paths());
                    for workers in [1,2]{
                        let actual=generator.bind(py).call_method1(stringify!($bulk),
                            (vec![Some(molecule.clone_ref(py)),None,Some(molecule.clone_ref(py))],workers))?
                            .extract::<Vec<Option<Py<$value>>>>()?;
                        let expected=native.$bulk(&[Some(&source),None,Some(&source)],workers).unwrap();
                        assert_eq!(actual.len(),expected.len());
                        for (actual,expected)in actual.iter().zip(expected.iter()){
                            assert_eq!(actual.as_ref().map(|v|v.borrow(py).inner.clone()),*expected);
                        }
                    }
                    let invalid=ck::TopologicalTorsionCallParams{custom_atom_invariants:Some(vec![]),..Default::default()};
                    let expected=source.[<fingerprint_topological_torsion $suffix _with_generator>](&native,Some(&invalid),None).err().unwrap().to_string();
                    let invalid=Py::new(py,TopologicalTorsionCallParams::from_rust(py,invalid)?)?;
                    let error=molecule.bind(py).call_method1(method,(&generator,&invalid)).unwrap_err();
                    assert!(error.is_instance(py, <TopologicalTorsionReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                    assert_eq!(error.value(py).str()?.to_str()?,expected);
                    assert!(error.value(py).getattr("source")?.is_instance(<TopologicalTorsionError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
                }};
            }
            output!("",fingerprints,Fingerprint);
            output!(_count,counts,SparseCountFingerprint32);
            output!(_sparse,sparse_fingerprints,SparseBitFingerprint);
            output!(_sparse_count,sparse_counts,SparseCountFingerprint);
            assert_eq!(call.borrow(py).to_rust(py)?,options);
            assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(),before);
            drop(generator);
            settings.setattr("torsion_atom_count",4)?;
            assert_eq!(settings.getattr("torsion_atom_count")?.extract::<u32>()?,4);
            let error=TopologicalTorsionFingerprintGenerator::from_json("{bad").err().unwrap();
            let generator_error=error.value(py).getattr("source")?;
            let json_error=generator_error.getattr("source")?;
            let parse=json_error.getattr("source")?;
            assert!(parse.is_instance(<SerdeJsonError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(parse.is(&json_error.getattr("__cause__")?));
            Ok(())
        }).unwrap();
    }
}

pub(crate) fn fingerprint_topological_torsion_with_generator(
    molecule: &Molecule,
    py: Python<'_>,
    generator: &TopologicalTorsionFingerprintGenerator,
    params: Option<&TopologicalTorsionCallParams>,
    additional_output: Option<&mut FingerprintAdditionalOutput>,
) -> PyResult<Fingerprint> {
    let params = params.map(|value| value.to_rust(py)).transpose()?;
    py.detach(|| {
        molecule
            .inner
            .fingerprint_topological_torsion_with_generator(
                &generator.inner,
                params.as_ref(),
                additional_output.map(|v| &mut v.inner),
            )
    })
    .map(|inner| Fingerprint { inner })
    .map_err(TopologicalTorsionReadError::from_rust)
}

pub(crate) fn fingerprint_topological_torsion_count_with_generator(
    molecule: &Molecule,
    py: Python<'_>,
    generator: &TopologicalTorsionFingerprintGenerator,
    params: Option<&TopologicalTorsionCallParams>,
    additional_output: Option<&mut FingerprintAdditionalOutput>,
) -> PyResult<SparseCountFingerprint32> {
    let params = params.map(|value| value.to_rust(py)).transpose()?;
    py.detach(|| {
        molecule
            .inner
            .fingerprint_topological_torsion_count_with_generator(
                &generator.inner,
                params.as_ref(),
                additional_output.map(|v| &mut v.inner),
            )
    })
    .map(|inner| SparseCountFingerprint32 { inner })
    .map_err(TopologicalTorsionReadError::from_rust)
}

pub(crate) fn fingerprint_topological_torsion_sparse_with_generator(
    molecule: &Molecule,
    py: Python<'_>,
    generator: &TopologicalTorsionFingerprintGenerator,
    params: Option<&TopologicalTorsionCallParams>,
    additional_output: Option<&mut FingerprintAdditionalOutput>,
) -> PyResult<SparseBitFingerprint> {
    let params = params.map(|value| value.to_rust(py)).transpose()?;
    py.detach(|| {
        molecule
            .inner
            .fingerprint_topological_torsion_sparse_with_generator(
                &generator.inner,
                params.as_ref(),
                additional_output.map(|v| &mut v.inner),
            )
    })
    .map(|inner| SparseBitFingerprint { inner })
    .map_err(TopologicalTorsionReadError::from_rust)
}

pub(crate) fn fingerprint_topological_torsion_sparse_count_with_generator(
    molecule: &Molecule,
    py: Python<'_>,
    generator: &TopologicalTorsionFingerprintGenerator,
    params: Option<&TopologicalTorsionCallParams>,
    additional_output: Option<&mut FingerprintAdditionalOutput>,
) -> PyResult<SparseCountFingerprint> {
    let params = params.map(|value| value.to_rust(py)).transpose()?;
    py.detach(|| {
        molecule
            .inner
            .fingerprint_topological_torsion_sparse_count_with_generator(
                &generator.inner,
                params.as_ref(),
                additional_output.map(|v| &mut v.inner),
            )
    })
    .map(|inner| SparseCountFingerprint { inner })
    .map_err(TopologicalTorsionReadError::from_rust)
}
