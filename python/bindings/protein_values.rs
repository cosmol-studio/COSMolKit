use crate::bio_errors::BioStructureError;
use crate::bio_rows::{BioAtomRow, BioResidueRow};
use crate::bio_source_values::{AltLocLabel, ResidueKind};
use crate::bio_structure::BioStructure;
use crate::bio_values::{AtomName, BioChainRow, BioCoordinateFormat, ChainKind, ChainSourceIds};
use crate::elements::Element;
use crate::error_framework::error_projection;
use crate::residue_values::{ResidueCode, ResidueInfo};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;
use std::sync::Arc;

error_projection!(ProteinProjectionError, ProteinProjectionErrorKind, ck::ProteinProjectionError, "bio", {
    Structure(_),
    MissingEntityMapping { source: u32 => |id: &ck::BioEntityId| id.value() },
    NonAminoAcidResidue { index: usize },
}, cause = |inner: &ck::ProteinProjectionError| match inner {
    ck::ProteinProjectionError::Structure(source) => Some(BioStructureError::from_rust(source.clone())),
    ck::ProteinProjectionError::MissingEntityMapping { .. } | ck::ProteinProjectionError::NonAminoAcidResidue { .. } => None,
}, methods = {
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    fn cause(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Option<crate::error_framework::ErrorObject<BioStructureError>>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
});

#[cfg(test)]
mod error_tests {
    use super::*;

    #[test]
    fn source_payload_and_real_cause_remain_distinct() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let error = ProteinProjectionError::from_rust(
                ck::ProteinProjectionError::MissingEntityMapping {
                    source: ck::BioEntityId::new(17),
                },
            );
            assert_eq!(error.value(py).getattr("source")?.extract::<u32>()?, 17);
            assert!(error.value(py).getattr("cause")?.is_none());
            assert!(error.cause(py).is_none());
            let error = ProteinProjectionError::from_rust(ck::ProteinProjectionError::Structure(
                ck::BioStructureError::RowIndexTooLarge { value: 23 },
            ));
            let cause = error.cause(py).unwrap();
            assert!(cause.is_instance(py, <BioStructureError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert!(error.value(py).getattr("cause")?.is(cause.value(py)));
            assert!(error.value(py).getattr("source")?.is_none());
            assert_eq!(cause.value(py).getattr("value")?.extract::<usize>()?, 23);
            assert!(error.value(py).setattr("cause", py.None()).is_err());
            Ok(())
        })
        .unwrap();
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct ProteinSelectionSummary {
        inner: ck::ProteinSelectionSummary,
    }
}
#[pymethods]
impl ProteinSelectionSummary {
    #[getter]
    fn chains(&self) -> usize {
        self.inner.chains
    }
    #[getter]
    fn residues(&self) -> usize {
        self.inner.residues
    }
    #[getter]
    fn atoms(&self) -> usize {
        self.inner.atoms
    }
}

#[pyclass(module = "cosmolkit", from_py_object)]
#[derive(Clone)]
pub struct Protein {
    pub(crate) inner: Arc<ck::Protein>,
}
#[pymethods]
impl Protein {
    #[staticmethod]
    #[pyo3(signature=(text, params=crate::arguments::Argument::Missing, *, max_line_length=crate::arguments::Argument::Missing, check_non_ascii=crate::arguments::Argument::Missing, ignore_ter=crate::arguments::Argument::Missing, split_chain_on_ter=crate::arguments::Argument::Missing, skip_remarks=crate::arguments::Argument::Missing))]
    fn from_pdb(
        py: Python<'_>,
        text: &str,
        params: crate::arguments::Argument<Option<PyRef<'_, crate::params::BioPdbReadParams>>>,
        max_line_length: crate::arguments::Argument<Option<i32>>,
        check_non_ascii: crate::arguments::Argument<Option<bool>>,
        ignore_ter: crate::arguments::Argument<Option<bool>>,
        split_chain_on_ter: crate::arguments::Argument<Option<bool>>,
        skip_remarks: crate::arguments::Argument<Option<bool>>,
    ) -> PyResult<Self> {
        let keyword_form = max_line_length.is_supplied()
            || check_non_ascii.is_supplied()
            || ignore_ter.is_supplied()
            || split_chain_on_ter.is_supplied()
            || skip_remarks.is_supplied();
        crate::arguments::configuration_arguments!(params; max_line_length, check_non_ascii, ignore_ter, split_chain_on_ter, skip_remarks);
        let params = params.as_deref();
        if keyword_form {
            Self::_from_pdb_keywords(
                py,
                text,
                max_line_length,
                check_non_ascii,
                ignore_ter,
                split_chain_on_ter,
                skip_remarks,
            )
        } else {
            Self::_from_pdb_params(py, text, params)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"Protein.from_pdb"}"#]
    #[pyo3(signature=(text, params=None))]
    fn _from_pdb_params(
        py: Python<'_>,
        text: &str,
        params: Option<&crate::params::BioPdbReadParams>,
    ) -> PyResult<Self> {
        crate::bio_reader_calls::_protein_from_pdb_params(py, text, params)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"Protein.from_pdb"}"#]
    #[pyo3(signature=(text, *, max_line_length=None, check_non_ascii=None, ignore_ter=None, split_chain_on_ter=None, skip_remarks=None))]
    fn _from_pdb_keywords(
        py: Python<'_>,
        text: &str,
        max_line_length: Option<i32>,
        check_non_ascii: Option<bool>,
        ignore_ter: Option<bool>,
        split_chain_on_ter: Option<bool>,
        skip_remarks: Option<bool>,
    ) -> PyResult<Self> {
        crate::bio_reader_calls::_protein_from_pdb_keywords(
            py,
            text,
            max_line_length,
            check_non_ascii,
            ignore_ter,
            split_chain_on_ter,
            skip_remarks,
        )
    }
    #[staticmethod]
    fn from_pdb_with_params(
        py: Python<'_>,
        text: &str,
        params: &crate::params::BioPdbReadParams,
    ) -> PyResult<Self> {
        Self::_from_pdb_params(py, text, Some(params))
    }
    #[staticmethod]
    #[pyo3(signature=(text, params=crate::arguments::Argument::Missing, *, format=crate::arguments::Argument::Missing, source_name=crate::arguments::Argument::Missing))]
    fn from_text(
        py: Python<'_>,
        text: &str,
        params: crate::arguments::Argument<Option<PyRef<'_, crate::bio_values::BioReadParams>>>,
        format: crate::arguments::Argument<Option<crate::bio_values::BioCoordinateFormatInput>>,
        source_name: crate::arguments::Argument<Option<String>>,
    ) -> PyResult<Self> {
        let keyword_form = format.is_supplied() || source_name.is_supplied();
        crate::arguments::configuration_arguments!(params; format, source_name);
        let params = params.as_deref();
        if keyword_form {
            Self::_from_text_keywords(py, text, format, source_name)
        } else {
            Self::_from_text_params(py, text, params)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"Protein.from_text"}"#]
    #[pyo3(signature=(text, params=None))]
    fn _from_text_params(
        py: Python<'_>,
        text: &str,
        params: Option<&crate::bio_values::BioReadParams>,
    ) -> PyResult<Self> {
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

        py.detach(|| match params {
            Some(params) => ck::Protein::from_text_with_params(text, &params),
            None => ck::Protein::from_text(text),
        })
        .map(|inner| Self {
            inner: Arc::new(inner),
        })
        .map_err(crate::bio_read_errors::ProteinReadError::from_rust)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"Protein.from_text"}"#]
    #[pyo3(signature=(text, *, format=None, source_name=None))]
    fn _from_text_keywords(
        py: Python<'_>,
        text: &str,
        format: Option<crate::bio_values::BioCoordinateFormatInput>,
        source_name: Option<String>,
    ) -> PyResult<Self> {
        let params = crate::bio_read_errors::keywords(format, source_name)?;
        Self::_from_text_params(py, text, Some(&params))
    }
    #[staticmethod]
    fn from_text_with_params(
        py: Python<'_>,
        text: &str,
        params: &crate::bio_values::BioReadParams,
    ) -> PyResult<Self> {
        Self::_from_text_params(py, text, Some(params))
    }
    #[staticmethod]
    #[pyo3(signature=(path, *, format=None))]
    fn read(
        py: Python<'_>,
        path: std::path::PathBuf,
        format: Option<crate::bio_values::BioCoordinateFormatInput>,
    ) -> PyResult<Self> {
        let format = format
            .map(crate::bio_values::BioCoordinateFormatInput::into_rust)
            .transpose()?;
        py.detach(|| match format {
            Some(format) => ck::Protein::read_with_format(&path, format),
            None => ck::Protein::read(&path),
        })
        .map(|inner| Self {
            inner: Arc::new(inner),
        })
        .map_err(crate::bio_read_errors::ProteinReadError::from_rust)
    }
    #[staticmethod]
    fn read_with_format(
        py: Python<'_>,
        path: std::path::PathBuf,
        format: crate::bio_values::BioCoordinateFormatInput,
    ) -> PyResult<Self> {
        Self::read(py, path, Some(format))
    }
    #[staticmethod]
    fn from_mmcif(py: Python<'_>, text: &str) -> PyResult<Self> {
        crate::bio_reader_calls::protein_from_mmcif(py, text)
    }
    #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, sanitize=crate::arguments::Argument::Missing, remove_hs=crate::arguments::Argument::Missing, flavor=crate::arguments::Argument::Missing, proximity_bonding=crate::arguments::Argument::Missing))]
    fn to_molecule(
        &self,
        py: Python<'_>,
        params: crate::arguments::Argument<Option<PyRef<'_, crate::params::BioMoleculeParams>>>,
        sanitize: crate::arguments::Argument<Option<bool>>,
        remove_hs: crate::arguments::Argument<Option<bool>>,
        flavor: crate::arguments::Argument<Option<crate::pdb_flags::PdbReadFlags>>,
        proximity_bonding: crate::arguments::Argument<Option<bool>>,
    ) -> PyResult<crate::native::Molecule> {
        let keyword_form = sanitize.is_supplied()
            || remove_hs.is_supplied()
            || flavor.is_supplied()
            || proximity_bonding.is_supplied();
        crate::arguments::configuration_arguments!(params; sanitize, remove_hs, flavor, proximity_bonding);
        let params = params.as_deref();
        if keyword_form {
            self._to_molecule_keywords(py, sanitize, remove_hs, flavor, proximity_bonding)
        } else {
            self._to_molecule_params(py, params)
        }
    }
    #[doc = r#"{"role":"overload","target":"Protein.to_molecule"}"#]
    #[pyo3(signature=(params=None))]
    fn _to_molecule_params(
        &self,
        py: Python<'_>,
        params: Option<&crate::params::BioMoleculeParams>,
    ) -> PyResult<crate::native::Molecule> {
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

        py.detach(|| match params {
            Some(params) => self.inner.to_molecule_with_params(&params),
            None => self.inner.to_molecule(),
        })
        .map(|inner| crate::native::Molecule { inner })
        .map_err(crate::bio_molecule_conversion::BioMoleculeError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"Protein.to_molecule"}"#]
    #[pyo3(signature=(*, sanitize=None, remove_hs=None, flavor=None, proximity_bonding=None))]
    fn _to_molecule_keywords(
        &self,
        py: Python<'_>,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        flavor: Option<crate::pdb_flags::PdbReadFlags>,
        proximity_bonding: Option<bool>,
    ) -> PyResult<crate::native::Molecule> {
        let params = crate::bio_molecule_conversion::keywords(
            sanitize,
            remove_hs,
            flavor,
            proximity_bonding,
        );
        self._to_molecule_params(py, Some(&params))
    }
    fn to_molecule_with_params(
        &self,
        py: Python<'_>,
        params: &crate::params::BioMoleculeParams,
    ) -> PyResult<crate::native::Molecule> {
        self._to_molecule_params(py, Some(params))
    }
    fn with_translated_coordinates(&self, offset: [f64; 3]) -> PyResult<Self> {
        self.inner
            .with_translated_coordinates(offset)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(crate::bio_selection_errors::BioOperationError::from_rust)
    }
    fn translate_(&mut self, offset: [f64; 3]) -> PyResult<()> {
        Arc::make_mut(&mut self.inner)
            .translate_(offset)
            .map_err(crate::bio_selection_errors::BioOperationError::from_rust)
    }
    fn selected_atom_ids(
        &self,
        selection: &crate::bio_selection_values::BioSelection,
    ) -> PyResult<Vec<u32>> {
        self.inner
            .selected_atom_ids(&selection.inner)
            .map(|ids| ids.into_iter().map(ck::BioAtomId::value).collect())
            .map_err(crate::bio_selection_errors::BioSelectionMatchError::from_rust)
    }
    fn with_selection(
        &self,
        selection: &crate::bio_selection_values::BioSelection,
    ) -> PyResult<Self> {
        self.inner
            .with_selection(&selection.inner)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(crate::bio_selection_errors::BioOperationError::from_rust)
    }
    fn retain_selection_(
        &mut self,
        selection: &crate::bio_selection_values::BioSelection,
    ) -> PyResult<()> {
        Arc::make_mut(&mut self.inner)
            .retain_selection_(&selection.inner)
            .map_err(crate::bio_selection_errors::BioOperationError::from_rust)
    }
    #[getter]
    fn input_format(&self) -> BioCoordinateFormat {
        self.inner.input_format().into()
    }
    fn as_bio_structure(&self) -> BioStructure {
        BioStructure {
            inner: self.inner.as_bio_structure().clone(),
        }
    }
    /// Return a COW snapshot without invalidating aliases to the Python value.
    fn into_bio_structure(&self) -> BioStructure {
        BioStructure {
            inner: self.inner.as_ref().clone().into_bio_structure(),
        }
    }
    #[getter]
    fn num_models(&self) -> usize {
        self.inner.num_models()
    }
    #[getter]
    fn num_chains(&self) -> usize {
        self.inner.num_chains()
    }
    #[getter]
    fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    fn selection_summary(&self) -> ProteinSelectionSummary {
        ProteinSelectionSummary {
            inner: self.inner.selection_summary(),
        }
    }
    #[getter]
    fn chains(&self) -> Vec<ProteinChainRef> {
        self.inner
            .chains()
            .map(|value| ProteinChainRef {
                owner: self.inner.clone(),
                index: value.id().value() as usize,
            })
            .collect()
    }
    #[getter]
    fn residues(&self) -> Vec<ProteinResidueRef> {
        self.inner
            .residues()
            .map(|value| ProteinResidueRef {
                owner: self.inner.clone(),
                index: value.id().value() as usize,
            })
            .collect()
    }
    #[getter]
    fn atoms(&self) -> Vec<ProteinAtomRef> {
        self.inner
            .atoms()
            .map(|value| ProteinAtomRef {
                owner: self.inner.clone(),
                index: value.id().value() as usize,
            })
            .collect()
    }
    fn chain(&self, index: usize) -> Option<ProteinChainRef> {
        self.inner.chain(index).map(|value| ProteinChainRef {
            owner: self.inner.clone(),
            index: value.id().value() as usize,
        })
    }
    fn __len__(&self) -> usize {
        self.inner.num_chains()
    }
    fn __getitem__(&self, index: crate::indexing::Integer) -> PyResult<ProteinChainRef> {
        let index = crate::indexing::normalize(index, self.inner.num_chains())?;
        self.chain(index)
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err(index))
    }
}

// Each native ref retains the immutable Rust owner and a private index emitted
// by that owner's canonical iterator. Reborrowing is O(1): all three Rust
// iterators override nth to forward directly to their Range cursor. No unsafe
// lifetime extension, duplicated hierarchy storage or chemistry is needed.
macro_rules! owned_view {
    ($name:ident, $rust:ident, $iter:ident) => {
        #[pyclass(module = "cosmolkit", frozen, from_py_object)]
        #[derive(Clone)]
        pub struct $name {
            owner: Arc<ck::Protein>,
            index: usize,
        }
        impl $name {
            fn view(&self) -> ck::$rust<'_> {
                self.owner
                    .$iter()
                    .nth(self.index)
                    .expect("private index belongs to retained immutable protein owner")
            }
        }
    };
}
owned_view!(ProteinChainRef, ProteinChainRef, chains);
owned_view!(ProteinResidueRef, ProteinResidueRef, residues);
owned_view!(ProteinAtomRef, ProteinAtomRef, atoms);

#[pymethods]
impl ProteinChainRef {
    #[getter]
    fn id(&self) -> u32 {
        self.view().id().value()
    }
    #[getter]
    fn row(&self) -> BioChainRow {
        BioChainRow {
            inner: self.view().row().clone(),
        }
    }
    #[getter]
    fn kind(&self) -> ChainKind {
        self.view().kind().into()
    }
    #[getter]
    fn source(&self) -> ChainSourceIds {
        ChainSourceIds {
            inner: self.view().source().clone(),
        }
    }
    #[getter]
    fn residues(&self) -> Vec<ProteinResidueRef> {
        self.view()
            .residues()
            .map(|value| ProteinResidueRef {
                owner: self.owner.clone(),
                index: value.id().value() as usize,
            })
            .collect()
    }
    #[getter]
    fn atoms(&self) -> Vec<ProteinAtomRef> {
        self.view()
            .atoms()
            .map(|value| ProteinAtomRef {
                owner: self.owner.clone(),
                index: value.id().value() as usize,
            })
            .collect()
    }
}
#[pymethods]
impl ProteinResidueRef {
    #[getter]
    fn id(&self) -> u32 {
        self.view().id().value()
    }
    #[getter]
    fn row(&self) -> BioResidueRow {
        BioResidueRow {
            inner: self.view().row().clone(),
        }
    }
    #[getter]
    fn name(&self) -> String {
        self.view().name().as_str().to_owned()
    }
    #[getter]
    fn kind(&self) -> ResidueKind {
        self.view().kind().into()
    }
    #[getter]
    fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: self.view().info(),
        }
    }
    #[getter]
    fn code(&self) -> ResidueCode {
        self.view().code().into()
    }
    #[getter]
    fn one_letter_code(&self) -> char {
        self.view().one_letter_code()
    }
    #[getter]
    fn fasta_code(&self) -> char {
        self.view().fasta_code()
    }
    #[getter]
    fn is_standard(&self) -> bool {
        self.view().is_standard()
    }
    #[getter]
    fn chain(&self) -> ProteinChainRef {
        ProteinChainRef {
            owner: self.owner.clone(),
            index: self.view().chain().id().value() as usize,
        }
    }
    #[getter]
    fn atoms(&self) -> Vec<ProteinAtomRef> {
        self.view()
            .atoms()
            .map(|value| ProteinAtomRef {
                owner: self.owner.clone(),
                index: value.id().value() as usize,
            })
            .collect()
    }
}
#[pymethods]
impl ProteinAtomRef {
    #[getter]
    fn id(&self) -> u32 {
        self.view().id().value()
    }
    #[getter]
    fn row(&self) -> BioAtomRow {
        BioAtomRow {
            inner: self.view().row().clone(),
        }
    }
    #[getter]
    fn name(&self) -> AtomName {
        AtomName {
            inner: self.view().name(),
        }
    }
    #[getter]
    fn element(&self) -> Element {
        Element {
            inner: self.view().element(),
        }
    }
    #[getter]
    fn altloc(&self) -> Option<AltLocLabel> {
        self.view().altloc().map(|inner| AltLocLabel { inner })
    }
    #[getter]
    fn residue(&self) -> ProteinResidueRef {
        ProteinResidueRef {
            owner: self.owner.clone(),
            index: self.view().residue().id().value() as usize,
        }
    }
    #[getter]
    fn position(&self) -> [f64; 3] {
        self.view().position()
    }
}
