//! Structural reader conveniences forward the same canonical Rust factories.
use crate::bio_read_errors::ProteinReadError;
use crate::params::BioPdbReadParams;
use crate::protein_values::Protein;
use cosmolkit as ck;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use std::sync::Arc;

#[pyfunction]
pub fn protein_from_mmcif(py: Python<'_>, text: &str) -> PyResult<Protein> {
    py.detach(|| ck::Protein::from_mmcif(text))
        .map(|inner| Protein {
            inner: Arc::new(inner),
        })
        .map_err(ProteinReadError::from_rust)
}

#[pyfunction]
#[pyo3(signature=(text,params=crate::arguments::Argument::Missing,*,max_line_length=crate::arguments::Argument::Missing,check_non_ascii=crate::arguments::Argument::Missing,ignore_ter=crate::arguments::Argument::Missing,split_chain_on_ter=crate::arguments::Argument::Missing,skip_remarks=crate::arguments::Argument::Missing))]
pub fn protein_from_pdb(
    py: Python<'_>,
    text: &str,
    params: crate::arguments::Argument<Option<PyRef<'_, BioPdbReadParams>>>,
    max_line_length: crate::arguments::Argument<Option<i32>>,
    check_non_ascii: crate::arguments::Argument<Option<bool>>,
    ignore_ter: crate::arguments::Argument<Option<bool>>,
    split_chain_on_ter: crate::arguments::Argument<Option<bool>>,
    skip_remarks: crate::arguments::Argument<Option<bool>>,
) -> PyResult<Protein> {
    let keyword_form = max_line_length.is_supplied()
        || check_non_ascii.is_supplied()
        || ignore_ter.is_supplied()
        || split_chain_on_ter.is_supplied()
        || skip_remarks.is_supplied();
    crate::arguments::configuration_arguments!(params; max_line_length, check_non_ascii, ignore_ter, split_chain_on_ter, skip_remarks);
    let params = params.as_deref();
    if keyword_form {
        _protein_from_pdb_keywords(
            py,
            text,
            max_line_length,
            check_non_ascii,
            ignore_ter,
            split_chain_on_ter,
            skip_remarks,
        )
    } else {
        _protein_from_pdb_params(py, text, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"protein_from_pdb"}"#]
#[pyo3(signature=(text,params=None))]
pub fn _protein_from_pdb_params(
    py: Python<'_>,
    text: &str,
    params: Option<&BioPdbReadParams>,
) -> PyResult<Protein> {
    let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

    py.detach(|| match params {
        Some(params) => ck::Protein::from_pdb_with_params(text, &params),
        None => ck::Protein::from_pdb(text),
    })
    .map(|inner| Protein {
        inner: Arc::new(inner),
    })
    .map_err(ProteinReadError::from_rust)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"protein_from_pdb"}"#]
#[pyo3(signature=(text,*,max_line_length=None,check_non_ascii=None,ignore_ter=None,split_chain_on_ter=None,skip_remarks=None))]
pub fn _protein_from_pdb_keywords(
    py: Python<'_>,
    text: &str,
    max_line_length: Option<i32>,
    check_non_ascii: Option<bool>,
    ignore_ter: Option<bool>,
    split_chain_on_ter: Option<bool>,
    skip_remarks: Option<bool>,
) -> PyResult<Protein> {
    let params = crate::bio_pdb_reader::keywords(
        max_line_length,
        check_non_ascii,
        ignore_ter,
        split_chain_on_ter,
        skip_remarks,
    );
    _protein_from_pdb_params(py, text, Some(&params))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bio_read_errors::{BioReadError, JsonError, MmjsonReadError, MmjsonStructureError};
    use crate::bio_structure::BioStructure;
    use crate::bio_values::BioReadParams;
    use crate::io_source_errors::IoError;
    use pyo3::types::{PyDict, PyModule};
    use std::error::Error;
    use std::path::PathBuf;

    #[test]
    fn generic_bio_readers_match_native_factories_formats_files_causes_and_owned_results() {
        Python::initialize();
        Python::attach(|py|->PyResult<()> {
            let pdb=concat!(
                "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n",
                "HETATM    2  O   HOH A   2       8.000   9.000  10.000  1.00 20.00           O  \nEND\n");
            let expected_bio=ck::BioStructure::from_pdb(pdb).unwrap();
            let expected_protein=ck::Protein::from_pdb(pdb).unwrap();
            let cif=expected_bio.to_mmcif().unwrap();
            let bio=py.get_type::<BioStructure>();
            let protein=py.get_type::<Protein>();
            let read=Py::new(py,BioReadParams{inner:crate::parameters::NativeField::Owned(ck::BioReadParams{format:ck::BioCoordinateFormat::Pdb,source_name:"explicit\0source".into()})})?;
            let options=PyDict::new(py);options.set_item("format","pdb")?;options.set_item("source_name","explicit\0source")?;
            assert_eq!(bio.call_method1("from_text",(pdb,))?.extract::<Py<BioStructure>>()?.borrow(py).inner,ck::BioStructure::from_text(pdb).unwrap());
            let expected_named_bio=ck::BioStructure::from_text_with_params(pdb,&read.borrow(py).to_rust(py)?).unwrap();
            for actual in [bio.call_method1("from_text",(pdb,&read))?,bio.call_method("from_text",(pdb,),Some(&options))?,bio.call_method1("from_text_with_params",(pdb,&read))?] {
                assert_eq!(actual.extract::<Py<BioStructure>>()?.borrow(py).inner,expected_named_bio);
            }
            assert_eq!(protein.call_method1("from_text",(pdb,))?.extract::<Py<Protein>>()?.borrow(py).inner.as_ref(),&ck::Protein::from_text(pdb).unwrap());
            let expected_named_protein=ck::Protein::from_text_with_params(pdb,&read.borrow(py).to_rust(py)?).unwrap();
            for actual in [protein.call_method1("from_text",(pdb,&read))?,protein.call_method("from_text",(pdb,),Some(&options))?,protein.call_method1("from_text_with_params",(pdb,&read))?] {
                assert_eq!(actual.extract::<Py<Protein>>()?.borrow(py).inner.as_ref(),&expected_named_protein);
            }
            assert_eq!(protein.call_method1("from_pdb",(pdb,))?.extract::<Py<Protein>>()?.borrow(py).inner.as_ref(),&expected_protein);
            let configured=PyDict::new(py);configured.set_item("skip_remarks",true)?;
            let pdb_params=Py::new(py,BioPdbReadParams{inner:crate::parameters::NativeField::Owned(ck::BioPdbReadParams{skip_remarks:true,..Default::default()})})?;
            for actual in [protein.call_method1("from_pdb",(pdb,&pdb_params))?,protein.call_method("from_pdb",(pdb,),Some(&configured))?,protein.call_method1("from_pdb_with_params",(pdb,&pdb_params))?] {
                assert_eq!(actual.extract::<Py<Protein>>()?.borrow(py).inner.as_ref(),&ck::Protein::from_pdb_with_params(pdb,&pdb_params.borrow(py).to_rust(py)?).unwrap());
            }
            assert!(protein.call_method("from_pdb",(pdb,&pdb_params),Some(&configured)).unwrap_err().is_instance_of::<PyTypeError>(py));
            let api=PyModule::new(py,"bio_reader_test")?;
            api.add_function(wrap_pyfunction!(protein_from_pdb,&api)?)?;
            api.add_function(wrap_pyfunction!(protein_from_mmcif,&api)?)?;
            assert_eq!(api.call_method1("protein_from_pdb",(pdb,))?.extract::<Py<Protein>>()?.borrow(py).inner.as_ref(),&expected_protein);
            for actual in [protein.call_method1("from_mmcif",(&cif,))?,api.call_method1("protein_from_mmcif",(&cif,))?] {
                assert_eq!(actual.extract::<Py<Protein>>()?.borrow(py).inner.as_ref(),&ck::Protein::from_mmcif(&cif).unwrap());
            }
            assert!(bio.call_method("from_text",(pdb,&read),Some(&options)).unwrap_err().is_instance_of::<PyTypeError>(py));
            assert!(protein.call_method("from_text",(pdb,&read),Some(&options)).unwrap_err().is_instance_of::<PyTypeError>(py));

            let folder=py.import("tempfile")?.getattr("TemporaryDirectory")?.call0()?;
            let directory=PathBuf::from(folder.getattr("name")?.extract::<String>()?);
            let path=directory.join("source.pdb");std::fs::write(&path,pdb).unwrap();
            let py_path=py.import("pathlib")?.getattr("Path")?.call1((&path,))?;
            let file_options=PyDict::new(py);file_options.set_item("format","pdb")?;
            for actual in [bio.call_method1("read",(&path,))?,bio.call_method1("read",(&py_path,))?,bio.call_method("read",(&py_path,),Some(&file_options))?,bio.call_method1("read_with_format",(&path,"pdb"))?] {
                assert_eq!(actual.extract::<Py<BioStructure>>()?.borrow(py).inner,ck::BioStructure::read(&path).unwrap());
            }
            for actual in [protein.call_method1("read",(&path,))?,protein.call_method1("read",(&py_path,))?,protein.call_method("read",(&py_path,),Some(&file_options))?,protein.call_method1("read_with_format",(&path,"pdb"))?] {
                assert_eq!(actual.extract::<Py<Protein>>()?.borrow(py).inner.as_ref(),&ck::Protein::read(&path).unwrap());
            }
            let missing=directory.join("missing.pdb");
            let native=ck::BioStructure::read(&missing).unwrap_err();
            let error=bio.call_method1("read",(&missing,)).unwrap_err();
            assert!(error.is_instance_of::<BioReadError>(py));
            assert_eq!(error.value(py).str()?.to_str()?,native.to_string());
            assert_eq!(error.value(py).getattr("path")?.extract::<PathBuf>()?,missing);
            let cause=error.value(py).getattr("source")?.unbind();
            assert!(cause.bind(py).is_instance_of::<IoError>());
            assert!(cause.bind(py).is(error.value(py).getattr("__cause__")?.as_any()));
            drop(error);drop(native);
            assert!(cause.bind(py).getattr("kind")?.is_instance_of::<crate::io_source_errors::IoErrorKind>());

            let bad_utf8=directory.join("bytes.pdb");
            std::fs::write(&bad_utf8,[0xff,0]).unwrap();
            let native=ck::BioStructure::read(&bad_utf8).unwrap_err();
            let error=bio.call_method1("read",(&bad_utf8,)).unwrap_err();
            assert_eq!(error.value(py).str()?.to_str()?,native.to_string());
            let source=error.value(py).getattr("source")?.unbind();
            assert!(source.bind(py).is_instance_of::<crate::utf8_errors::FromUtf8Error>());
            drop(error);drop(native);
            assert_eq!(source.bind(py).call_method0("as_bytes")?.extract::<Vec<u8>>()?,vec![0xff,0]);

            let native=ck::BioReadError::ChemComp(Box::new(ck::ChemCompModelError::Atom(ck::ChemCompAtomError::MissingValue{row:7,column:8})));
            let display=native.to_string();
            let error=BioReadError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?,display);
            let source=error.value(py).getattr("source")?;
            assert!(source.is_instance_of::<crate::bio_read_errors::ChemCompModelError>());
            let atom=source.getattr("source")?.unbind();
            assert!(atom.bind(py).is_instance_of::<crate::bio_read_errors::ChemCompAtomError>());
            drop(error);drop(source);
            assert_eq!(atom.bind(py).getattr("row")?.extract::<usize>()?,7);
            assert_eq!(atom.bind(py).getattr("column")?.extract::<usize>()?,8);
            assert!(atom.bind(py).setattr("row",9).is_err());

            for (kind,escaped) in [(ck::MmjsonUnexpectedType::False,crate::bio_read_errors::MmjsonUnexpectedType::False_),
                (ck::MmjsonUnexpectedType::True,crate::bio_read_errors::MmjsonUnexpectedType::True_)] {
                let native=ck::BioReadError::Mmjson(Box::new(ck::MmjsonReadError::Structure{
                    source_name:"json".into(),error:ck::MmjsonStructureError::Text(ck::MmjsonValueError::UnexpectedType(kind))}));
                let display=native.to_string();
                let error=BioReadError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?,display);
                let payload=error.value(py).getattr("source")?.getattr("source")?.getattr("source")?;
                assert!(payload.getattr("value")?.extract::<crate::bio_read_errors::MmjsonUnexpectedType>()?==escaped);
            }

            // mmJSON is an actual native parser with concrete errors, not serde_json.
            for text in ["{\n\"data_x\": ","[1]"] {
                let params=ck::BioReadParams{format:ck::BioCoordinateFormat::Mmjson,source_name:"original-json".into()};
                let native=ck::BioStructure::from_text_with_params(text,&params).unwrap_err();
                let display=native.to_string();
                let native_json=native.source().unwrap().downcast_ref::<ck::MmjsonReadError>().unwrap();
                let native_leaf=native_json.source().unwrap();
                let options=PyDict::new(py);options.set_item("format","mmjson")?;options.set_item("source_name","original-json")?;
                let error=bio.call_method("from_text",(text,),Some(&options)).unwrap_err();
                assert_eq!(error.value(py).str()?.to_str()?,display);
                let source=error.value(py).getattr("source")?;
                assert!(source.is_instance_of::<MmjsonReadError>());
                assert_eq!(source.getattr("source_name")?.extract::<String>()?,"original-json");
                let leaf=source.getattr("source")?;
                assert!(leaf.is(&source.getattr("error")?));
                if let Some(expected)=native_leaf.downcast_ref::<ck::JsonError>() {
                    assert!(leaf.is_instance_of::<JsonError>());
                    assert_eq!(leaf.getattr("line")?.extract::<usize>()?,expected.line());
                    assert_eq!(leaf.getattr("column")?.extract::<usize>()?,expected.column());
                    assert_eq!(leaf.getattr("byte_offset")?.extract::<usize>()?,expected.byte_offset());
                    assert_eq!(leaf.getattr("argument")?.extract::<Option<i32>>()?,expected.argument());
                    assert!(leaf.getattr("kind")?.extract::<crate::bio_read_errors::JsonErrorKind>()?==expected.kind().into());
                    // A native record has kind/position fields, no invented
                    // enum discriminant to freeze. Check its real descriptors.
                    assert!(!leaf.hasattr("variant")?);
                    assert!(leaf.setattr("kind",py.None()).is_err());
                    assert!(leaf.setattr("line",0).is_err());
                } else {
                    assert!(native_leaf.is::<ck::MmjsonStructureError>());assert!(leaf.is_instance_of::<MmjsonStructureError>());
                    assert!(leaf.setattr("variant",py.None()).is_err());
                }
                assert_eq!(leaf.str()?.to_str()?,native_leaf.to_string());
            }
            folder.call_method0("cleanup")?;
            Ok(())
        }).unwrap();
    }
}
