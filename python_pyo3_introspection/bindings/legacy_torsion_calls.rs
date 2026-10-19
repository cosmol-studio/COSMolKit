//! The distinct canonical legacy algorithms; no alternative chemistry in bindings.
use crate::numerical_params::LegacyTopologicalTorsionParams;
use cosmolkit as ck;

pub(crate) fn keywords(
    torsion_atom_count: Option<u32>,
    include_chirality: Option<bool>,
    fp_size: Option<u32>,
    bits_per_entry: Option<u32>,
    from_atoms: Option<Vec<u32>>,
    ignore_atoms: Option<Vec<u32>>,
    custom_atom_invariants: Option<Vec<u32>>,
) -> ck::LegacyTopologicalTorsionParams {
    let defaults = ck::LegacyTopologicalTorsionParams::default();
    ck::LegacyTopologicalTorsionParams::new(
        torsion_atom_count.unwrap_or(defaults.torsion_atom_count),
        include_chirality.unwrap_or(defaults.include_chirality),
        fp_size.unwrap_or(defaults.fp_size),
        bits_per_entry.unwrap_or(defaults.bits_per_entry),
        from_atoms,
        ignore_atoms,
        custom_atom_invariants,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fingerprint_values::SparseCountFingerprint;
    use crate::morgan_errors::TopologicalTorsionReadError;
    use crate::native::{Fingerprint, Molecule};
    use pyo3::exceptions::PyTypeError;
    use pyo3::prelude::*;
    use pyo3::types::PyDict;

    #[test]
    fn legacy_torsion_forms_and_ids_preserve_native_defaults_options_errors_and_source() {
        Python::initialize();
        Python::attach(|py|->PyResult<()>{
            let source=ck::Molecule::from_smiles("CCCO").unwrap();
            let before=source.to_binary().unwrap();
            let molecule=Py::new(py,Molecule{inner:source.clone()})?;
            let options=ck::LegacyTopologicalTorsionParams::new(
                3,true,256,2,Some(vec![0,3]),Some(vec![]),Some(vec![11,12,13,14]));
            let params=Py::new(py,LegacyTopologicalTorsionParams::from_rust(py, options.clone())?)?;
            let kwargs=PyDict::new(py);
            kwargs.set_item("torsion_atom_count",3)?;
            kwargs.set_item("include_chirality",true)?;
            kwargs.set_item("fp_size",256)?;
            kwargs.set_item("bits_per_entry",2)?;
            kwargs.set_item("from_atoms",vec![0,3])?;
            kwargs.set_item("ignore_atoms",Vec::<u32>::new())?;
            kwargs.set_item("custom_atom_invariants",vec![11,12,13,14])?;
            macro_rules! output{
                ($suffix:tt,$value:ident)=>{paste::paste!{
                    let method=stringify!([<fingerprint_topological_torsion $suffix _legacy>]);
                    let alias=stringify!([<fingerprint_topological_torsion $suffix _legacy_with_params>]);
                    let expected_default=source.[<fingerprint_topological_torsion $suffix _legacy>]().unwrap();
                    assert_eq!(molecule.bind(py).call_method0(method)?.extract::<Py<$value>>()?.borrow(py).inner,expected_default);
                    let expected=source.[<fingerprint_topological_torsion $suffix _legacy_with_params>](&options).unwrap();
                    for value in [
                        molecule.bind(py).call_method1(method,(&params,))?,
                        molecule.bind(py).call_method1(alias,(&params,))?,
                        molecule.bind(py).call_method(method,(),Some(&kwargs))?,
                    ]{assert_eq!(value.extract::<Py<$value>>()?.borrow(py).inner,expected);}
                    assert!(molecule.bind(py).call_method(method,(&params,),Some(&kwargs)).unwrap_err().is_instance_of::<PyTypeError>(py));
                    let mut invalid=options.clone();
                    invalid.custom_atom_invariants=Some(vec![]);
                    let expected=source.[<fingerprint_topological_torsion $suffix _legacy_with_params>](&invalid).err().unwrap().to_string();
                    let invalid=Py::new(py,LegacyTopologicalTorsionParams::from_rust(py, invalid)?)?;
                    let error=molecule.bind(py).call_method1(method,(&invalid,)).unwrap_err();
                    assert!(error.is_instance_of::<TopologicalTorsionReadError>(py));
                    assert_eq!(error.value(py).str()?.to_str()?,expected);
                    assert!(error.value(py).getattr("source")?.is(&error.value(py).getattr("__cause__")?));
                }};
            }
            output!("",Fingerprint);
            output!(_count,SparseCountFingerprint);
            output!(_sparse_count,SparseCountFingerprint);
            assert_eq!(molecule.bind(py).call_method0("topological_torsion_ids")?.extract::<Vec<u64>>()?,source.topological_torsion_ids().unwrap());
            for size in [3_u32,4,5]{
                let expected=source.topological_torsion_ids_with_params(size).unwrap();
                assert_eq!(molecule.bind(py).call_method1("topological_torsion_ids",(size,))?.extract::<Vec<u64>>()?,expected);
                assert_eq!(molecule.bind(py).call_method1("topological_torsion_ids_with_params",(size,))?.extract::<Vec<u64>>()?,expected);
            }
            assert_eq!(params.borrow(py).to_rust(py)?,options);
            assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(),before);
            Ok(())
        }).unwrap();
    }
}
