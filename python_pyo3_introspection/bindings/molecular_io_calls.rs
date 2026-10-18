//! Whole molecular text/file adapters; all format and chemistry behavior stays native.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    use pyo3::types::PyDict;
    #[test]
    fn complete_molecular_io_files_records_and_native_suppliers_preserve_options_and_errors() {
        Python::initialize();
        Python::attach(|py|->PyResult<()>{
let source=ck::Molecule::from_smiles("CO").unwrap().with_added_3d_conformer(vec![vec![0.,0.,0.],vec![1.,0.2,0.4]]).unwrap();
let before=source.to_binary().unwrap();
let mol=Py::new(py,Molecule{inner:source.clone()})?;
let write=ck::MolBlockWriteParams{precision:3,include_coordinates:false,..Default::default()};
let p=Py::new(py,crate::io_batch_params::MolBlockWriteParams{inner:write})?;
let kw=PyDict::new(py);kw.set_item("precision",3)?;kw.set_item("include_coordinates",false)?;
macro_rules! configured{($name:ident,$alias:ident)=>{{
let expected=source.$alias(&write);
for actual in [mol.bind(py).call_method1(stringify!($name),(&p,)),mol.bind(py).call_method(stringify!($name),(),Some(&kw)),mol.bind(py).call_method1(stringify!($alias),(&p,))]{
match(&expected,actual){(Ok(a),Ok(b))=>assert_eq!(a,&b.extract::<String>()?),(Err(a),Err(b))=>{assert!(b.is_instance_of::<crate::molecular_io_errors::MolecularIoError>(py));assert_eq!(a.to_string(),b.value(py).str()?.to_str()?);},_=>panic!("configured writer branch differs: {}",stringify!($name))}
}
assert!(mol.bind(py).call_method(stringify!($name),(&p,),Some(&kw)).unwrap_err().is_instance_of::<pyo3::exceptions::PyTypeError>(py));
}}}
configured!(to_mol,to_mol_with_params);configured!(to_sdf,to_sdf_with_params);configured!(to_sdf_2d,to_sdf_2d_with_params);configured!(to_sdf_3d,to_sdf_3d_with_params);
assert_eq!(mol.bind(py).call_method0("to_pdb")?.extract::<String>()?,source.to_pdb().unwrap());
assert_eq!(mol.bind(py).call_method0("to_xyz")?.extract::<String>()?,source.to_xyz().unwrap());
let pdb=source.to_pdb().unwrap();let xyz=source.to_xyz().unwrap();
let expected=ck::Molecule::from_pdb(&pdb).unwrap();
let actual=py.get_type::<Molecule>().call_method1("from_pdb",(&pdb,))?.extract::<Py<Molecule>>()?;
assert_eq!(actual.borrow(py).inner,expected);
assert_eq!(mol_from_pdb(py,&pdb,None,None,None,None,None)?.inner,expected);
assert_eq!(py.get_type::<Molecule>().call_method1("from_xyz_block",(&xyz,))?.extract::<Py<Molecule>>()?.borrow(py).inner,ck::Molecule::from_xyz_block(&xyz).unwrap());
let sd=source.to_sdf().unwrap();let record=Py::new(py,crate::sdf_records::SdfRecord::from_rust(ck::SdfRecord::from_sdf(&sd).unwrap()))?;
for name in ["to_mol","to_sdf"]{let v=record.bind(py).call_method0(name)?.extract::<String>()?;assert!(!v.is_empty());assert!(record.bind(py).call_method(name,(&p,),Some(&kw)).unwrap_err().is_instance_of::<pyo3::exceptions::PyTypeError>(py));}
let directory=std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("artifacts").join(format!("io-complete-{}",std::process::id()));std::fs::create_dir_all(&directory).unwrap();
let path=directory.join("source.sdf");std::fs::write(&path,format!("{sd}{sd}")).unwrap();
let native=ck::SdfDataset::open(path_text(&path)?).unwrap();
let dataset=py.get_type::<crate::sdf_sources::SdfDataset>().call_method1("open",(&path,))?;
assert_eq!(dataset.call_method0("len")?.extract::<usize>()?,native.len());
assert_eq!(dataset.call_method1("record_text",(0,))?.extract::<String>()?,native.record_text(0).unwrap());
let records=dataset.call_method0("iter")?.try_iter()?.collect::<PyResult<Vec<_>>>()?;assert_eq!(records.len(),2);
drop(dataset);assert_eq!(records[0].call_method0("molecule")?.extract::<Py<Molecule>>()?.borrow(py).inner,native.record(0).unwrap().molecule().unwrap().clone());
let stream=py.get_type::<crate::sdf_sources::SdfRecordStream>().call_method1("open",(&path,))?;
assert!(stream.call_method0("next_record")?.extract::<Option<Py<crate::sdf_records::SdfRecord>>>()?.is_some());assert_eq!(stream.call_method0("records_consumed")?.extract::<usize>()?,1);
assert!(stream.call_method0("next_record")?.extract::<Option<Py<crate::sdf_records::SdfRecord>>>()?.is_some());assert!(stream.call_method0("next_record")?.is_none());
let reader=py.get_type::<crate::sdf_sources::SdfReader>().call_method1("open",(&path,))?;
let iterator=reader.call_method1("batches",(1,crate::io_batch_params::BatchErrorMode::Strict,Some(1)))?;
drop(reader);let batches=iterator.try_iter()?.collect::<PyResult<Vec<_>>>()?;assert_eq!(batches.len(),2);
for (name,native_method)in [("write_mol","mol"),("write_sdf","sdf"),("write_xyz","xyz")]{let p=directory.join(format!("out.{native_method}"));mol.bind(py).call_method1(name,(&p,))?;assert!(!std::fs::read(p).unwrap().is_empty());}
for name in ["read_sdf","read_mol"]{assert_eq!(py.get_type::<Molecule>().call_method1(name,(&path,))?.extract::<Py<Molecule>>()?.borrow(py).inner,ck::Molecule::read_sdf(path_text(&path)?).unwrap());}
let missing=directory.join("not-here").join("missing.sdf");let e=py.get_type::<Molecule>().call_method1("read_sdf",(&missing,)).unwrap_err();assert!(e.is_instance_of::<crate::molecular_io_errors::MolecularIoError>(py));let cause=e.value(py).getattr("source")?;assert!(cause.is_instance_of::<crate::io_source_errors::IoError>());assert!(cause.is(&e.value(py).getattr("__cause__")?));drop(e);assert!(!cause.str()?.is_empty()?);
assert_eq!(mol.borrow(py).inner.to_binary().unwrap(),before);
Ok(())
}).unwrap();
    }
}

use crate::io_batch_params::*;
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;
pub(crate) fn path_text(path: &std::path::Path) -> PyResult<&str> {
    path.to_str().ok_or_else(|| {
        pyo3::exceptions::PyUnicodeError::new_err("native path interface requires a UTF-8 path")
    })
}
pub(crate) fn mol_write_keywords(
    format: Option<SdfFormatInput>,
    force_2d: Option<bool>,
    include_stereo: Option<bool>,
    kekulize: Option<bool>,
    precision: Option<usize>,
    coordinate_selection: Option<MolCoordinateSelectionInput>,
    include_coordinates: Option<bool>,
) -> PyResult<crate::io_batch_params::MolBlockWriteParams> {
    let defaults = ck::MolBlockWriteParams::default();
    Ok(crate::io_batch_params::MolBlockWriteParams {
        inner: ck::MolBlockWriteParams {
            format: format
                .map(SdfFormatInput::into_rust)
                .transpose()?
                .unwrap_or(defaults.format),
            force_2d: force_2d.unwrap_or(defaults.force_2d),
            include_stereo: include_stereo.unwrap_or(defaults.include_stereo),
            kekulize: kekulize.unwrap_or(defaults.kekulize),
            precision: precision.unwrap_or(defaults.precision),
            coordinate_selection: coordinate_selection
                .map(MolCoordinateSelectionInput::into_rust)
                .transpose()?
                .unwrap_or(defaults.coordinate_selection),
            include_coordinates: include_coordinates.unwrap_or(defaults.include_coordinates),
        },
    })
}
pub(crate) fn pdb_write_keywords(
    coordinate_selection: Option<MolCoordinateSelectionInput>,
    all_conformers: Option<bool>,
    flavor: Option<u32>,
) -> PyResult<crate::io_batch_params::PdbWriteParams> {
    let defaults = ck::PdbWriteParams::default();
    Ok(crate::io_batch_params::PdbWriteParams {
        inner: ck::PdbWriteParams {
            coordinate_selection: coordinate_selection
                .map(MolCoordinateSelectionInput::into_rust)
                .transpose()?
                .unwrap_or(defaults.coordinate_selection),
            all_conformers: all_conformers.unwrap_or(defaults.all_conformers),
            flavor: flavor.unwrap_or(defaults.flavor),
        },
    })
}
pub(crate) fn mol2_read_keywords(
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    variant: Option<Mol2TypeInput>,
    cleanup_substructures: Option<bool>,
) -> PyResult<crate::io_batch_params::Mol2ReadParams> {
    let defaults = ck::Mol2ReadParams::default();
    Ok(crate::io_batch_params::Mol2ReadParams {
        inner: ck::Mol2ReadParams {
            sanitize: sanitize.unwrap_or(defaults.sanitize),
            remove_hs: remove_hs.unwrap_or(defaults.remove_hs),
            variant: variant
                .map(Mol2TypeInput::into_rust)
                .transpose()?
                .unwrap_or(defaults.variant),
            cleanup_substructures: cleanup_substructures.unwrap_or(defaults.cleanup_substructures),
        },
    })
}
pub(crate) fn xyz_write_keywords(
    conformer_id: Option<usize>,
    precision: Option<u32>,
) -> PyResult<crate::params::XyzWriteParams> {
    let defaults = ck::XyzWriteParams::default();
    Ok(crate::params::XyzWriteParams {
        inner: ck::XyzWriteParams {
            conformer_id: conformer_id.or(defaults.conformer_id),
            precision: precision.unwrap_or(defaults.precision),
        },
    })
}
pub(crate) fn avalon_keywords(
    n_bits: Option<u32>,
    is_query: Option<bool>,
    bit_flags: Option<u32>,
) -> PyResult<crate::numerical_params::AvalonFingerprintParams> {
    let defaults = ck::AvalonFingerprintParams::default();
    Ok(crate::numerical_params::AvalonFingerprintParams {
        inner: ck::AvalonFingerprintParams {
            n_bits: n_bits.unwrap_or(defaults.n_bits),
            is_query: is_query.unwrap_or(defaults.is_query),
            bit_flags: ck::AvalonFingerprintFlags::from_bits_retain(
                bit_flags.unwrap_or(defaults.bit_flags.bits()),
            ),
        },
    })
}
pub(crate) fn bio_keywords(
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    flavor: Option<u32>,
    proximity_bonding: Option<bool>,
) -> PyResult<crate::params::BioMoleculeParams> {
    Ok(crate::bio_molecule_conversion::keywords(
        sanitize,
        remove_hs,
        flavor,
        proximity_bonding,
    ))
}
pub(crate) fn sdf_keywords(
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    strict_parsing: Option<bool>,
    expand_attachment_points: Option<bool>,
    process_property_lists: Option<bool>,
    coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
) -> PyResult<crate::selectors::SdfReadParams> {
    crate::sdf_records::configuration(
        sanitize,
        remove_hs,
        strict_parsing,
        expand_attachment_points,
        process_property_lists,
        coordinate_mode,
    )
}
#[pyfunction]
#[pyo3(signature=(text,params=None,*,sanitize=None,remove_hs=None,flavor=None,proximity_bonding=None))]
pub fn mol_from_pdb(
    py: Python<'_>,
    text: &str,
    params: Option<&crate::params::BioMoleculeParams>,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    flavor: Option<u32>,
    proximity_bonding: Option<bool>,
) -> PyResult<Molecule> {
    if params.is_some()
        && (sanitize.is_some()
            || remove_hs.is_some()
            || flavor.is_some()
            || proximity_bonding.is_some())
    {
        return Err(pyo3::exceptions::PyTypeError::new_err(
            "Params and keyword configuration are mutually exclusive",
        ));
    }
    if sanitize.is_some() || remove_hs.is_some() || flavor.is_some() || proximity_bonding.is_some()
    {
        _mol_from_pdb_keywords(py, text, sanitize, remove_hs, flavor, proximity_bonding)
    } else {
        _mol_from_pdb_params(py, text, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_pdb"}"#]
#[pyo3(signature=(text,params=None))]
pub fn _mol_from_pdb_params(
    py: Python<'_>,
    text: &str,
    params: Option<&crate::params::BioMoleculeParams>,
) -> PyResult<Molecule> {
    Molecule::from_pdb(py, text, params, None, None, None, None)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_pdb"}"#]
#[pyo3(signature=(text,*,sanitize=None,remove_hs=None,flavor=None,proximity_bonding=None))]
pub fn _mol_from_pdb_keywords(
    py: Python<'_>,
    text: &str,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    flavor: Option<u32>,
    proximity_bonding: Option<bool>,
) -> PyResult<Molecule> {
    Molecule::from_pdb(
        py,
        text,
        None,
        sanitize,
        remove_hs,
        flavor,
        proximity_bonding,
    )
}
#[pyfunction]
#[pyo3(signature=(text,params=None,*,sanitize=None,remove_hs=None,flavor=None,proximity_bonding=None))]
pub fn mol_from_mmcif(
    py: Python<'_>,
    text: &str,
    params: Option<&crate::params::BioMoleculeParams>,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    flavor: Option<u32>,
    proximity_bonding: Option<bool>,
) -> PyResult<Molecule> {
    if params.is_some()
        && (sanitize.is_some()
            || remove_hs.is_some()
            || flavor.is_some()
            || proximity_bonding.is_some())
    {
        return Err(pyo3::exceptions::PyTypeError::new_err(
            "Params and keyword configuration are mutually exclusive",
        ));
    }
    if sanitize.is_some() || remove_hs.is_some() || flavor.is_some() || proximity_bonding.is_some()
    {
        _mol_from_mmcif_keywords(py, text, sanitize, remove_hs, flavor, proximity_bonding)
    } else {
        _mol_from_mmcif_params(py, text, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_mmcif"}"#]
#[pyo3(signature=(text,params=None))]
pub fn _mol_from_mmcif_params(
    py: Python<'_>,
    text: &str,
    params: Option<&crate::params::BioMoleculeParams>,
) -> PyResult<Molecule> {
    Molecule::from_mmcif(py, text, params, None, None, None, None)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_mmcif"}"#]
#[pyo3(signature=(text,*,sanitize=None,remove_hs=None,flavor=None,proximity_bonding=None))]
pub fn _mol_from_mmcif_keywords(
    py: Python<'_>,
    text: &str,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    flavor: Option<u32>,
    proximity_bonding: Option<bool>,
) -> PyResult<Molecule> {
    Molecule::from_mmcif(
        py,
        text,
        None,
        sanitize,
        remove_hs,
        flavor,
        proximity_bonding,
    )
}
#[pyfunction]
#[pyo3(signature=(text,params=None,*,sanitize=None,remove_hs=None,variant=None,cleanup_substructures=None))]
pub fn mol_from_mol2(
    py: Python<'_>,
    text: &str,
    params: Option<&crate::io_batch_params::Mol2ReadParams>,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    variant: Option<Mol2TypeInput>,
    cleanup_substructures: Option<bool>,
) -> PyResult<Molecule> {
    if params.is_some()
        && (sanitize.is_some()
            || remove_hs.is_some()
            || variant.is_some()
            || cleanup_substructures.is_some())
    {
        return Err(pyo3::exceptions::PyTypeError::new_err(
            "Params and keyword configuration are mutually exclusive",
        ));
    }
    if sanitize.is_some()
        || remove_hs.is_some()
        || variant.is_some()
        || cleanup_substructures.is_some()
    {
        _mol_from_mol2_keywords(
            py,
            text,
            sanitize,
            remove_hs,
            variant,
            cleanup_substructures,
        )
    } else {
        _mol_from_mol2_params(py, text, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_mol2"}"#]
#[pyo3(signature=(text,params=None))]
pub fn _mol_from_mol2_params(
    py: Python<'_>,
    text: &str,
    params: Option<&crate::io_batch_params::Mol2ReadParams>,
) -> PyResult<Molecule> {
    Molecule::from_mol2(py, text, params, None, None, None, None)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_mol2"}"#]
#[pyo3(signature=(text,*,sanitize=None,remove_hs=None,variant=None,cleanup_substructures=None))]
pub fn _mol_from_mol2_keywords(
    py: Python<'_>,
    text: &str,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    variant: Option<Mol2TypeInput>,
    cleanup_substructures: Option<bool>,
) -> PyResult<Molecule> {
    Molecule::from_mol2(
        py,
        text,
        None,
        sanitize,
        remove_hs,
        variant,
        cleanup_substructures,
    )
}
#[pyfunction]
pub fn mol_from_xyz_block(py: Python<'_>, text: &str) -> PyResult<Molecule> {
    Molecule::from_xyz_block(py, text)
}
