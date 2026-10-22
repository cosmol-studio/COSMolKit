use crate::bio_errors::BioStructureError;
use crate::bio_geometry_values::{BioAssembly, BioCoordinateBlock, BioCrystalInfo, BioNcsOperator};
use crate::bio_metadata_values::BioMetadata;
use crate::bio_relationships::{BioCisPep, BioConnection, BioHelix, BioModRes, BioSheet};
use crate::bio_rows::{BioAtomRow, BioEntityRow, BioResidueRow};
use crate::bio_scalar_metadata::BioStructureSourceState;
use crate::bio_source_values::{AltLocLabel, AltLocRequest};
use crate::bio_values::{
    AtomName, BioChainRow, BioCoordinateFormat, BioCoordinateFormatInput, BioModelRow, BioTransform,
};
use crate::elements::Element;
use cosmolkit as ck;
use pyo3::prelude::*;

// Complete transport of detached Rust parts. No implicit empty tables are
// invented: the Rust type has no Default, so every field is required here.
macro_rules! structure_values {
    (mmcif {$($mmcif:ident: $mmcif_ty:ty),+ $(,)?}, $($list:ident: $value:ident),+ $(,)?) => {
        crate::configuration::type_projection! {
            #[pyclass(module="cosmolkit", frozen, from_py_object)]
            #[derive(Clone)]
            pub struct BioStructureParts {
                pub(crate) inner: ck::BioStructureParts
            }
        }
        #[pymethods]
        impl BioStructureParts {
            #[new]
            #[pyo3(signature=(*, input_format, $($list,)+ metadata, source_state, coordinates, crystal))]
            fn new(input_format: BioCoordinateFormatInput, $($list: Vec<$value>,)+ metadata: BioMetadata,
                source_state: BioStructureSourceState, coordinates: BioCoordinateBlock,
                crystal: Option<BioCrystalInfo>) -> PyResult<Self> {
                Ok(Self { inner: ck::BioStructureParts { input_format: input_format.into_rust()?,
                    $($list: $list.into_iter().map(|value| value.inner).collect(),)+
                    metadata: metadata.inner, source_state: source_state.inner, coordinates: coordinates.inner,
                    crystal: crystal.map(|value| value.inner),
                } })
            }
            #[getter] fn input_format(&self) -> BioCoordinateFormat { self.inner.input_format.into() }
            $(#[getter] fn $list(&self) -> Vec<$value> {
                self.inner.$list.iter().cloned().map(|inner| $value { inner }).collect()
            })+
            #[getter] fn metadata(&self) -> BioMetadata { BioMetadata { inner: self.inner.metadata.clone() } }
            #[getter] fn source_state(&self) -> BioStructureSourceState { BioStructureSourceState { inner: self.inner.source_state.clone() } }
            #[getter] fn coordinates(&self) -> BioCoordinateBlock { BioCoordinateBlock { inner: self.inner.coordinates.clone() } }
            #[getter] fn crystal(&self) -> Option<BioCrystalInfo> { self.inner.crystal.clone().map(|inner| BioCrystalInfo { inner }) }
        }

        #[pyclass(module="cosmolkit", from_py_object)]
        #[derive(Clone)]
        pub struct BioStructure { pub(crate) inner: ck::BioStructure }
        #[pymethods]
        impl BioStructure {

#[staticmethod]
#[pyo3(signature=(text,params=crate::arguments::Argument::Missing,*,max_line_length=crate::arguments::Argument::Missing,check_non_ascii=crate::arguments::Argument::Missing,ignore_ter=crate::arguments::Argument::Missing,split_chain_on_ter=crate::arguments::Argument::Missing,skip_remarks=crate::arguments::Argument::Missing))]
fn from_pdb(py:Python<'_>,text:&str,params:crate::arguments::Argument<Option<PyRef<'_, crate::params::BioPdbReadParams>>>,max_line_length:crate::arguments::Argument<Option<i32>>,check_non_ascii:crate::arguments::Argument<Option<bool>>,ignore_ter:crate::arguments::Argument<Option<bool>>,split_chain_on_ter:crate::arguments::Argument<Option<bool>>,skip_remarks:crate::arguments::Argument<Option<bool>>)->PyResult<Self>{
        let keyword_form = max_line_length.is_supplied() || check_non_ascii.is_supplied() || ignore_ter.is_supplied() || split_chain_on_ter.is_supplied() || skip_remarks.is_supplied();
        crate::arguments::configuration_arguments!(params; max_line_length, check_non_ascii, ignore_ter, split_chain_on_ter, skip_remarks);
        let params = params.as_deref();
        if keyword_form { Self::_from_pdb_keywords(py, text, max_line_length, check_non_ascii, ignore_ter, split_chain_on_ter, skip_remarks) } else { Self::_from_pdb_params(py, text, params) }
    }
#[staticmethod]
fn from_pdb_with_params(py:Python<'_>,text:&str,params:&crate::params::BioPdbReadParams)->PyResult<Self>{
let params = pyo3::Python::attach(|py| params.to_rust(py))?;

py.detach(||ck::BioStructure::from_pdb_with_params(text,&params)).map(|inner|Self{inner}).map_err(crate::bio_pdb_reader::BioPdbReadError::from_rust)
}
#[staticmethod]
#[doc=r#"{"role":"overload","target":"BioStructure.from_pdb"}"#]
#[pyo3(signature=(text,params=None))]
fn _from_pdb_params(py:Python<'_>,text:&str,params:Option<&crate::params::BioPdbReadParams>)->PyResult<Self>{crate::bio_pdb_reader::_bio_from_pdb_params(py,text,params)}
#[staticmethod]
#[doc=r#"{"role":"overload","target":"BioStructure.from_pdb"}"#]
#[pyo3(signature=(text,*,max_line_length=None,check_non_ascii=None,ignore_ter=None,split_chain_on_ter=None,skip_remarks=None))]
fn _from_pdb_keywords(py:Python<'_>,text:&str,max_line_length:Option<i32>,check_non_ascii:Option<bool>,ignore_ter:Option<bool>,split_chain_on_ter:Option<bool>,skip_remarks:Option<bool>)->PyResult<Self>{crate::bio_pdb_reader::_bio_from_pdb_keywords(py,text,max_line_length,check_non_ascii,ignore_ter,split_chain_on_ter,skip_remarks)}

            #[staticmethod]
#[pyo3(signature=(text, params=crate::arguments::Argument::Missing, *, format=crate::arguments::Argument::Missing, source_name=crate::arguments::Argument::Missing))]
            fn from_text(py: Python<'_>, text: &str, params: crate::arguments::Argument<Option<PyRef<'_, crate::bio_values::BioReadParams>>>, format: crate::arguments::Argument<Option<BioCoordinateFormatInput>>, source_name: crate::arguments::Argument<Option<String>>) -> PyResult<Self> {
        let keyword_form = format.is_supplied() || source_name.is_supplied();
        crate::arguments::configuration_arguments!(params; format, source_name);
        let params = params.as_deref();
        if keyword_form { Self::_from_text_keywords(py, text, format, source_name) } else { Self::_from_text_params(py, text, params) }
    }
            #[staticmethod]
            #[doc = r#"{"role":"overload","target":"BioStructure.from_text"}"#]
            #[pyo3(signature=(text, params=None))]
            fn _from_text_params(py: Python<'_>, text: &str, params: Option<&crate::bio_values::BioReadParams>) -> PyResult<Self> {
let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

                py.detach(|| match params {
                    Some(params) => ck::BioStructure::from_text_with_params(text, &params),
                    None => ck::BioStructure::from_text(text),
                }).map(|inner| Self { inner }).map_err(crate::bio_read_errors::BioReadError::from_rust)
            }
            #[staticmethod]
            #[doc = r#"{"role":"overload","target":"BioStructure.from_text"}"#]
            #[pyo3(signature=(text, *, format=None, source_name=None))]
            fn _from_text_keywords(py: Python<'_>, text: &str, format: Option<BioCoordinateFormatInput>, source_name: Option<String>) -> PyResult<Self> {
                let params = crate::bio_read_errors::keywords(format, source_name)?;
                Self::_from_text_params(py, text, Some(&params))
            }
            #[staticmethod]
            fn from_text_with_params(py: Python<'_>, text: &str, params: &crate::bio_values::BioReadParams) -> PyResult<Self> {
                Self::_from_text_params(py, text, Some(params))
            }
            #[staticmethod]
            #[pyo3(signature=(path, *, format=None))]
            fn read(py: Python<'_>, path: std::path::PathBuf, format: Option<BioCoordinateFormatInput>) -> PyResult<Self> {
                let format = format.map(BioCoordinateFormatInput::into_rust).transpose()?;
                py.detach(|| match format {
                    Some(format) => ck::BioStructure::read_with_format(&path, format),
                    None => ck::BioStructure::read(&path),
                }).map(|inner| Self { inner }).map_err(crate::bio_read_errors::BioReadError::from_rust)
            }
            #[staticmethod]
            fn read_with_format(py: Python<'_>, path: std::path::PathBuf, format: BioCoordinateFormatInput) -> PyResult<Self> {
                Self::read(py, path, Some(format))
            }
            #[staticmethod]
            fn from_mmcif(py: Python<'_>, text: &str) -> PyResult<Self> {
                crate::bio_mmcif_reader::bio_from_mmcif(py, text)
            }
            #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, $($mmcif=crate::arguments::Argument::Missing),+))]
            fn to_mmcif(&self, params: crate::arguments::Argument<Option<PyRef<'_, crate::bio_write_configuration::BioMmcifWriteParams>>>, $($mmcif: crate::arguments::Argument<Option<$mmcif_ty>>),+) -> PyResult<String> {
                let configured = false $(|| $mmcif.is_supplied())+;
                crate::arguments::configuration_arguments!(params; $($mmcif),+);
                let params = params.as_deref();
                if configured { self._to_mmcif_keywords($($mmcif),+) }
                else { self._to_mmcif_params(params) }
            }
            #[doc = r#"{"role":"overload","target":"BioStructure.to_mmcif"}"#]
            #[pyo3(signature=(params=None))]
            fn _to_mmcif_params(&self, params: Option<&crate::bio_write_configuration::BioMmcifWriteParams>) -> PyResult<String> {
let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

                match params {
                    Some(params) => self.inner.to_mmcif_with_params(&params),
                    None => self.inner.to_mmcif(),
                }.map_err(crate::bio_mmcif_writer::BioMmcifWriteError::from_rust)
            }
            #[doc = r#"{"role":"overload","target":"BioStructure.to_mmcif"}"#]
            #[pyo3(signature=(*, $($mmcif=None),+))]
            fn _to_mmcif_keywords(&self, $($mmcif: Option<$mmcif_ty>),+) -> PyResult<String> {
                let defaults = ck::BioMmcifWriteParams::default();
                let params = crate::bio_write_configuration::BioMmcifWriteParams { inner: crate::parameters::NativeField::Owned(ck::BioMmcifWriteParams {
                    $($mmcif: $mmcif.unwrap_or(defaults.$mmcif)),+
                } )};
                self._to_mmcif_params(Some(&params))
            }
            fn to_mmcif_with_params(&self, params: &crate::bio_write_configuration::BioMmcifWriteParams) -> PyResult<String> {
let params = pyo3::Python::attach(|py| params.to_rust(py))?;

                self.inner.to_mmcif_with_params(&params).map_err(crate::bio_mmcif_writer::BioMmcifWriteError::from_rust)
            }
            #[pyo3(signature=(path, params=crate::arguments::Argument::Missing, *, $($mmcif=crate::arguments::Argument::Missing),+))]
            fn write_mmcif(&self, path: std::path::PathBuf, params: crate::arguments::Argument<Option<PyRef<'_, crate::bio_write_configuration::BioMmcifWriteParams>>>, $($mmcif: crate::arguments::Argument<Option<$mmcif_ty>>),+) -> PyResult<()> {
                let configured = false $(|| $mmcif.is_supplied())+;
                crate::arguments::configuration_arguments!(params; $($mmcif),+);
                let params = params.as_deref();
                if configured { self._write_mmcif_keywords(path, $($mmcif),+) }
                else { self._write_mmcif_params(path, params) }
            }
            #[doc = r#"{"role":"overload","target":"BioStructure.write_mmcif"}"#]
            #[pyo3(signature=(path, params=None))]
            fn _write_mmcif_params(&self, path: std::path::PathBuf, params: Option<&crate::bio_write_configuration::BioMmcifWriteParams>) -> PyResult<()> {
let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

                match params {
                    Some(params) => self.inner.write_mmcif_with_params(&path, &params),
                    None => self.inner.write_mmcif(&path),
                }.map_err(crate::bio_mmcif_writer::BioMmcifWriteError::from_rust)
            }
            #[doc = r#"{"role":"overload","target":"BioStructure.write_mmcif"}"#]
            #[pyo3(signature=(path, *, $($mmcif=None),+))]
            fn _write_mmcif_keywords(&self, path: std::path::PathBuf, $($mmcif: Option<$mmcif_ty>),+) -> PyResult<()> {
                let defaults = ck::BioMmcifWriteParams::default();
                let params = crate::bio_write_configuration::BioMmcifWriteParams { inner: crate::parameters::NativeField::Owned(ck::BioMmcifWriteParams {
                    $($mmcif: $mmcif.unwrap_or(defaults.$mmcif)),+
                } )};
                self._write_mmcif_params(path, Some(&params))
            }
            fn write_mmcif_with_params(&self, path: std::path::PathBuf, params: &crate::bio_write_configuration::BioMmcifWriteParams) -> PyResult<()> {
let params = pyo3::Python::attach(|py| params.to_rust(py))?;

                self.inner.write_mmcif_with_params(&path, &params).map_err(crate::bio_mmcif_writer::BioMmcifWriteError::from_rust)
            }
            #[staticmethod]
            fn from_parts(parts: BioStructureParts) -> PyResult<Self> {
                ck::BioStructure::from_parts(parts.inner).map(|inner| Self { inner }).map_err(BioStructureError::from_rust)
            }
            #[staticmethod]
            fn validate_parts(parts: &BioStructureParts) -> PyResult<()> {
                ck::BioStructure::validate_parts(&parts.inner).map_err(BioStructureError::from_rust)
            }
            fn validate(&self) -> PyResult<()> { self.inner.validate().map_err(BioStructureError::from_rust) }
            fn selected_atom_ids(&self, selection: &crate::bio_selection_values::BioSelection) -> PyResult<Vec<u32>> {
                self.inner.selected_atom_ids(&selection.inner)
                    .map(|ids| ids.into_iter().map(ck::BioAtomId::value).collect())
                    .map_err(crate::bio_selection_errors::BioSelectionMatchError::from_rust)
            }
            fn with_selection(&self, selection: &crate::bio_selection_values::BioSelection) -> PyResult<Self> {
                self.inner.with_selection(&selection.inner).map(|inner| Self { inner })
                    .map_err(crate::bio_selection_errors::BioOperationError::from_rust)
            }
            fn retain_selection_(&mut self, selection: &crate::bio_selection_values::BioSelection) -> PyResult<()> {
                self.inner.retain_selection_(&selection.inner)
                    .map_err(crate::bio_selection_errors::BioOperationError::from_rust)
            }
#[pyo3(signature=(params=crate::arguments::Argument::Missing, *, sanitize=crate::arguments::Argument::Missing, remove_hs=crate::arguments::Argument::Missing, flavor=crate::arguments::Argument::Missing, proximity_bonding=crate::arguments::Argument::Missing))]
            fn to_molecule(&self, py: Python<'_>, params: crate::arguments::Argument<Option<PyRef<'_, crate::params::BioMoleculeParams>>>, sanitize: crate::arguments::Argument<Option<bool>>, remove_hs: crate::arguments::Argument<Option<bool>>, flavor: crate::arguments::Argument<Option<crate::pdb_flags::PdbReadFlags>>, proximity_bonding: crate::arguments::Argument<Option<bool>>) -> PyResult<crate::native::Molecule> {
        let keyword_form = sanitize.is_supplied() || remove_hs.is_supplied() || flavor.is_supplied() || proximity_bonding.is_supplied();
        crate::arguments::configuration_arguments!(params; sanitize, remove_hs, flavor, proximity_bonding);
        let params = params.as_deref();
        if keyword_form { self._to_molecule_keywords(py, sanitize, remove_hs, flavor, proximity_bonding) } else { self._to_molecule_params(py, params) }
    }
            #[doc = r#"{"role":"overload","target":"BioStructure.to_molecule"}"#]
            #[pyo3(signature=(params=None))]
            fn _to_molecule_params(&self, py: Python<'_>, params: Option<&crate::params::BioMoleculeParams>) -> PyResult<crate::native::Molecule> {
let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

                py.detach(|| match params { Some(params) => self.inner.to_molecule_with_params(&params), None => self.inner.to_molecule() })
                    .map(|inner| crate::native::Molecule { inner }).map_err(crate::bio_molecule_conversion::BioMoleculeError::from_rust)
            }
            #[doc = r#"{"role":"overload","target":"BioStructure.to_molecule"}"#]
            #[pyo3(signature=(*, sanitize=None, remove_hs=None, flavor=None, proximity_bonding=None))]
            fn _to_molecule_keywords(&self, py: Python<'_>, sanitize: Option<bool>, remove_hs: Option<bool>, flavor: Option<crate::pdb_flags::PdbReadFlags>, proximity_bonding: Option<bool>) -> PyResult<crate::native::Molecule> {
                let params = crate::bio_molecule_conversion::keywords(sanitize, remove_hs, flavor, proximity_bonding);
                self._to_molecule_params(py, Some(&params))
            }
            fn to_molecule_with_params(&self, py: Python<'_>, params: &crate::params::BioMoleculeParams) -> PyResult<crate::native::Molecule> {
                self._to_molecule_params(py, Some(params))
            }
            fn with_translated_coordinates(&self, offset: [f64; 3]) -> PyResult<Self> {
                self.inner.with_translated_coordinates(offset).map(|inner| Self { inner })
                    .map_err(crate::bio_selection_errors::BioOperationError::from_rust)
            }
            fn translate_(&mut self, offset: [f64; 3]) -> PyResult<()> {
                self.inner.translate_(offset)
                    .map_err(crate::bio_selection_errors::BioOperationError::from_rust)
            }
#[pyo3(signature=(params=crate::arguments::Argument::Missing, *, ter_records=crate::arguments::Argument::Missing, numbered_ter=crate::arguments::Argument::Missing, ter_ignores_type=crate::arguments::Argument::Missing, preserve_serial=crate::arguments::Argument::Missing, end_record=crate::arguments::Argument::Missing))]
            fn to_pdb(&self, params: crate::arguments::Argument<Option<PyRef<'_, crate::params::BioPdbWriteParams>>>, ter_records: crate::arguments::Argument<Option<bool>>, numbered_ter: crate::arguments::Argument<Option<bool>>, ter_ignores_type: crate::arguments::Argument<Option<bool>>, preserve_serial: crate::arguments::Argument<Option<bool>>, end_record: crate::arguments::Argument<Option<bool>>) -> PyResult<String> {
        let keyword_form = ter_records.is_supplied() || numbered_ter.is_supplied() || ter_ignores_type.is_supplied() || preserve_serial.is_supplied() || end_record.is_supplied();
        crate::arguments::configuration_arguments!(params; ter_records, numbered_ter, ter_ignores_type, preserve_serial, end_record);
        let params = params.as_deref();
        if keyword_form { self._to_pdb_keywords(ter_records, numbered_ter, ter_ignores_type, preserve_serial, end_record) } else { self._to_pdb_params(params) }
    }
            #[doc = r#"{"role":"overload","target":"BioStructure.to_pdb"}"#]
            #[pyo3(signature=(params=None))]
            fn _to_pdb_params(&self, params: Option<&crate::params::BioPdbWriteParams>) -> PyResult<String> {
let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

                match params {
                    Some(params) => self.inner.to_pdb_with_params(&params),
                    None => self.inner.to_pdb(),
                }.map_err(crate::bio_pdb_writer::BioPdbWriteError::from_rust)
            }
            #[doc = r#"{"role":"overload","target":"BioStructure.to_pdb"}"#]
            #[pyo3(signature=(*, ter_records=None, numbered_ter=None, ter_ignores_type=None, preserve_serial=None, end_record=None))]
            fn _to_pdb_keywords(&self, ter_records: Option<bool>, numbered_ter: Option<bool>, ter_ignores_type: Option<bool>, preserve_serial: Option<bool>, end_record: Option<bool>) -> PyResult<String> {
                let params = crate::bio_pdb_writer::keywords(ter_records, numbered_ter, ter_ignores_type, preserve_serial, end_record);
                self._to_pdb_params(Some(&params))
            }
            fn to_pdb_with_params(&self, params: &crate::params::BioPdbWriteParams) -> PyResult<String> {
let params = pyo3::Python::attach(|py| params.to_rust(py))?;

                self.inner.to_pdb_with_params(&params).map_err(crate::bio_pdb_writer::BioPdbWriteError::from_rust)
            }
#[pyo3(signature=(path, params=crate::arguments::Argument::Missing, *, ter_records=crate::arguments::Argument::Missing, numbered_ter=crate::arguments::Argument::Missing, ter_ignores_type=crate::arguments::Argument::Missing, preserve_serial=crate::arguments::Argument::Missing, end_record=crate::arguments::Argument::Missing))]
            fn write_pdb(&self, path: std::path::PathBuf, params: crate::arguments::Argument<Option<PyRef<'_, crate::params::BioPdbWriteParams>>>, ter_records: crate::arguments::Argument<Option<bool>>, numbered_ter: crate::arguments::Argument<Option<bool>>, ter_ignores_type: crate::arguments::Argument<Option<bool>>, preserve_serial: crate::arguments::Argument<Option<bool>>, end_record: crate::arguments::Argument<Option<bool>>) -> PyResult<()> {
        let keyword_form = ter_records.is_supplied() || numbered_ter.is_supplied() || ter_ignores_type.is_supplied() || preserve_serial.is_supplied() || end_record.is_supplied();
        crate::arguments::configuration_arguments!(params; ter_records, numbered_ter, ter_ignores_type, preserve_serial, end_record);
        let params = params.as_deref();
        if keyword_form { self._write_pdb_keywords(path, ter_records, numbered_ter, ter_ignores_type, preserve_serial, end_record) } else { self._write_pdb_params(path, params) }
    }
            #[doc = r#"{"role":"overload","target":"BioStructure.write_pdb"}"#]
            #[pyo3(signature=(path, params=None))]
            fn _write_pdb_params(&self, path: std::path::PathBuf, params: Option<&crate::params::BioPdbWriteParams>) -> PyResult<()> {
let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

                match params {
                    Some(params) => self.inner.write_pdb_with_params(&path, &params),
                    None => self.inner.write_pdb(&path),
                }.map_err(crate::bio_pdb_writer::BioPdbWriteError::from_rust)
            }
            #[doc = r#"{"role":"overload","target":"BioStructure.write_pdb"}"#]
            #[pyo3(signature=(path, *, ter_records=None, numbered_ter=None, ter_ignores_type=None, preserve_serial=None, end_record=None))]
            fn _write_pdb_keywords(&self, path: std::path::PathBuf, ter_records: Option<bool>, numbered_ter: Option<bool>, ter_ignores_type: Option<bool>, preserve_serial: Option<bool>, end_record: Option<bool>) -> PyResult<()> {
                let params = crate::bio_pdb_writer::keywords(ter_records, numbered_ter, ter_ignores_type, preserve_serial, end_record);
                self._write_pdb_params(path, Some(&params))
            }
            fn write_pdb_with_params(&self, path: std::path::PathBuf, params: &crate::params::BioPdbWriteParams) -> PyResult<()> {
let params = pyo3::Python::attach(|py| params.to_rust(py))?;

                self.inner.write_pdb_with_params(&path, &params).map_err(crate::bio_pdb_writer::BioPdbWriteError::from_rust)
            }
            fn protein(&self) -> PyResult<crate::protein_values::Protein> {
                self.inner.protein().map(|inner| crate::protein_values::Protein { inner: std::sync::Arc::new(inner) }).map_err(crate::protein_values::ProteinProjectionError::from_rust)
            }
            /// Return detached parts without invalidating Python aliases of this value.
            /// The canonical Rust conversion operates on a COW snapshot.
            fn into_parts(&self) -> BioStructureParts { BioStructureParts { inner: self.inner.clone().into_parts() } }
            #[getter]
            fn input_format(&self) -> BioCoordinateFormat { self.inner.input_format().into() }
            #[getter] fn models(&self) -> Vec<crate::bio_hierarchy::BioModel> {
                self.inner.models().into_iter().map(|inner| crate::bio_hierarchy::BioModel { inner }).collect()
            }
            #[getter] fn chains(&self) -> Vec<crate::bio_hierarchy::BioChain> {
                self.inner.chains().into_iter().map(|inner| crate::bio_hierarchy::BioChain { inner }).collect()
            }
            #[getter] fn residues(&self) -> Vec<crate::bio_hierarchy::BioResidue> {
                self.inner.residues().into_iter().map(|inner| crate::bio_hierarchy::BioResidue { inner }).collect()
            }
            #[getter] fn atoms(&self) -> Vec<crate::bio_hierarchy::BioAtom> {
                self.inner.atoms().into_iter().map(|inner| crate::bio_hierarchy::BioAtom { inner }).collect()
            }
            #[getter] fn entities(&self) -> Vec<BioEntityRow> { self.inner.entities().iter().cloned().map(|inner| BioEntityRow { inner }).collect() }
            #[getter] fn connections(&self) -> Vec<BioConnection> { self.inner.connections().iter().cloned().map(|inner| BioConnection { inner }).collect() }
            #[getter] fn cispeps(&self) -> Vec<BioCisPep> { self.inner.cispeps().iter().cloned().map(|inner| BioCisPep { inner }).collect() }
            #[getter] fn mod_residues(&self) -> Vec<BioModRes> { self.inner.mod_residues().iter().cloned().map(|inner| BioModRes { inner }).collect() }
            #[getter] fn helices(&self) -> Vec<BioHelix> { self.inner.helices().iter().cloned().map(|inner| BioHelix { inner }).collect() }
            #[getter] fn sheets(&self) -> Vec<BioSheet> { self.inner.sheets().iter().cloned().map(|inner| BioSheet { inner }).collect() }
            #[getter] fn ncs_operators(&self) -> Vec<BioNcsOperator> { self.inner.ncs_operators().iter().cloned().map(|inner| BioNcsOperator { inner }).collect() }
            #[getter] fn assemblies(&self) -> Vec<BioAssembly> { self.inner.assemblies().iter().cloned().map(|inner| BioAssembly { inner }).collect() }
            #[getter]
            fn metadata(&self) -> BioMetadata { BioMetadata { inner: self.inner.metadata().clone() } }
            #[getter]
            fn source_state(&self) -> BioStructureSourceState { BioStructureSourceState { inner: self.inner.source_state().clone() } }
            #[getter]
            fn coordinates(&self) -> BioCoordinateBlock { BioCoordinateBlock { inner: self.inner.coordinates().clone() } }
            #[getter]
            fn crystal(&self) -> Option<BioCrystalInfo> { self.inner.crystal().cloned().map(|inner| BioCrystalInfo { inner }) }
            #[getter]
            fn num_models(&self) -> usize { self.inner.num_models() }
            fn __len__(&self) -> usize { self.inner.num_models() }
            fn model(&self, id: u32) -> Option<crate::bio_hierarchy::BioModel> {
                self.inner.model(ck::BioModelId::new(id)).map(|inner| crate::bio_hierarchy::BioModel { inner })
            }
            fn chain(&self, id: u32) -> Option<crate::bio_hierarchy::BioChain> {
                self.inner.chain(ck::BioChainId::new(id)).map(|inner| crate::bio_hierarchy::BioChain { inner })
            }
            fn residue(&self, id: u32) -> Option<crate::bio_hierarchy::BioResidue> {
                self.inner.residue(ck::BioResidueId::new(id)).map(|inner| crate::bio_hierarchy::BioResidue { inner })
            }
            fn atom(&self, id: u32) -> Option<crate::bio_hierarchy::BioAtom> {
                self.inner.atom(ck::BioAtomId::new(id)).map(|inner| crate::bio_hierarchy::BioAtom { inner })
            }
            fn __getitem__(&self, index: crate::indexing::Integer) -> PyResult<crate::bio_hierarchy::BioModel> {
                let index = crate::indexing::normalize(index, self.inner.num_models())?;
                Ok(crate::bio_hierarchy::BioModel { inner: self.inner.model(ck::BioModelId::new(index as u32)).expect("normalized model index") })
            }
            #[getter]
            fn num_chains(&self) -> usize { self.inner.num_chains() }
            #[getter]
            fn num_residues(&self) -> usize { self.inner.num_residues() }
            #[getter]
            fn num_atoms(&self) -> usize { self.inner.num_atoms() }
            #[getter]
            fn num_entities(&self) -> usize { self.inner.num_entities() }
            #[getter]
            fn name(&self) -> &str { self.inner.name() }
            #[getter]
            fn has_origx(&self) -> bool { self.inner.has_origx() }
            #[getter]
            fn origx(&self) -> BioTransform { BioTransform { inner: *self.inner.origx() } }
            #[getter]
            fn ncs_oper_identity_id(&self) -> Option<&str> { self.inner.ncs_oper_identity_id() }
            #[getter]
            fn resolution(&self) -> f64 { self.inner.resolution() }
            #[getter]
            fn ter_status(&self) -> u8 { self.inner.ter_status() }
            fn atom_position(&self, atom: u32) -> Option<[f64;3]> { self.inner.atom_position(ck::BioAtomId::new(atom)) }
            fn residue_atoms(&self, residue: u32) -> Option<Vec<crate::bio_hierarchy::BioAtom>> {
                self.inner.residue_atoms(ck::BioResidueId::new(residue)).map(|rows| rows.into_iter().map(|inner| crate::bio_hierarchy::BioAtom { inner }).collect())
            }
            fn find_entity(&self, source_id: &str) -> Option<(u32, BioEntityRow)> {
                self.inner.find_entity(source_id).map(|(id, value)| (id.value(), BioEntityRow { inner: value.clone() }))
            }
            fn find_entity_of_subchain(&self, subchain: &str) -> Option<(u32, BioEntityRow)> {
                self.inner.find_entity_of_subchain(subchain).map(|(id, value)| (id.value(), BioEntityRow { inner: value.clone() }))
            }
            #[pyo3(signature=(residue_id, name, request, element=None))]
            fn find_atom(&self, residue_id: u32, name: &AtomName, request: &AltLocRequest, element: Option<&Element>) -> Option<(u32, crate::bio_hierarchy::BioAtom)> {
                self.inner.find_atom(ck::BioResidueId::new(residue_id), name.inner, request.inner, element.map(|value| value.inner)).map(|(id, value)| (id.value(), crate::bio_hierarchy::BioAtom { inner: value }))
            }
            #[pyo3(signature=(residue_id, name, altloc=None))]
            fn atom_by_altloc(&self, residue_id: u32, name: &AtomName, altloc: Option<&AltLocLabel>) -> PyResult<(u32, crate::bio_hierarchy::BioAtom)> {
                self.inner.atom_by_altloc(ck::BioResidueId::new(residue_id), name.inner, altloc.map(|value| value.inner)).map(|(id, value)| (id.value(), crate::bio_hierarchy::BioAtom { inner: value })).map_err(BioStructureError::from_rust)
            }
        }
    };
}
structure_values! {
    mmcif {
        atoms: bool, block_name: bool, entry: bool, database_status: bool, author: bool,
        cell: bool, symmetry: bool, entity: bool, entity_poly: bool, struct_ref: bool,
        chem_comp: bool, exptl: bool, diffrn: bool, reflns: bool, refine: bool,
        title_keywords: bool, ncs: bool, struct_asym: bool, origx: bool, struct_conf: bool,
        struct_sheet: bool, struct_biol: bool, assembly: bool, conn: bool, cis: bool,
        modres: bool, scale: bool, atom_type: bool, entity_poly_seq: bool, tls: bool,
        software: bool, group_pdb: bool, auth_all: bool, prefer_pairs: bool, compact: bool,
        misuse_hash: bool, align_pairs: u16, align_loops: u16
    },
    models: BioModelRow,
    chains: BioChainRow,
    residues: BioResidueRow,
    atoms: BioAtomRow,
    entities: BioEntityRow,
    connections: BioConnection,
    cispeps: BioCisPep,
    mod_residues: BioModRes,
    helices: BioHelix,
    sheets: BioSheet,
    ncs_operators: BioNcsOperator,
    assemblies: BioAssembly,
}
