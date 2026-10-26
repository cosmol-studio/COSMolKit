use crate::chemical_values::enum_value;
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(EntityKind, { Unknown, Polymer, NonPolymer, Branched, Water });
enum_value!(PolymerKind, { Unknown, PeptideL, PeptideD, Dna, Rna, DnaRnaHybrid,
    SaccharideD, SaccharideL, Pna, CyclicPseudoPeptide, Other });
enum_value!(ResidueKind, { AminoAcid, Dna, Rna, Saccharide, Water, Buffer, Ligand, Unknown });

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Author residue sequence number paired with an optional insertion code.
    pub struct PdbSeqId {
        pub(crate) inner: ck::PdbSeqId,
    }
}
#[pymethods]
impl PdbSeqId {
    #[new]
    #[pyo3(signature=(seq_num, ins_code=None))]
    fn new(seq_num: i32, ins_code: Option<u8>) -> Self {
        Self {
            inner: ck::PdbSeqId::new(seq_num, ins_code),
        }
    }
    #[getter]
    /// Author residue sequence number.
    fn seq_num(&self) -> i32 {
        self.inner.seq_num()
    }
    #[getter]
    /// Insertion-code character, or None if absent.
    fn ins_code(&self) -> Option<char> {
        self.inner.ins_code().map(char::from)
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// A single native alternate-location label byte.
    pub struct AltLocLabel {
        pub(crate) inner: ck::AltLocLabel,
    }
}
#[pymethods]
impl AltLocLabel {
    #[new]
    fn new(value: u8) -> Self {
        Self {
            inner: ck::AltLocLabel::new(value),
        }
    }
    #[getter]
    /// Raw alternate-location label byte.
    fn value(&self) -> u8 {
        self.inner.value()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Original atom serial identifier, distinct from a local atom row index.
    pub struct AtomSourceIds {
        pub(crate) inner: ck::AtomSourceIds,
    }
}
#[pymethods]
impl AtomSourceIds {
    #[new]
    #[pyo3(signature=(serial=None))]
    fn new(serial: Option<i32>) -> Self {
        Self {
            inner: ck::AtomSourceIds::new(serial.map(ck::PdbAtomSerial::new)),
        }
    }
    #[getter]
    /// Source atom serial number, or None if absent.
    fn serial(&self) -> Option<i32> {
        self.inner.serial().map(|value| value.value())
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Source residue numbering, subchain, segment and entity identifiers.
    pub struct ResidueSourceIds {
        pub(crate) inner: ck::ResidueSourceIds,
    }
}
#[pymethods]
impl ResidueSourceIds {
    #[staticmethod]
    #[pyo3(signature=(seq_id=None, label_seq_id=None, segment_id=None, subchain_id=None, label_entity_id=None))]
    /// Create residue source identifiers, returning None when the native identifier validation rejects them.
    fn new(
        seq_id: Option<PdbSeqId>,
        label_seq_id: Option<i32>,
        segment_id: Option<[u8; 4]>,
        subchain_id: Option<String>,
        label_entity_id: Option<String>,
    ) -> Option<Self> {
        ck::ResidueSourceIds::new(
            seq_id.map(|value| value.inner),
            label_seq_id,
            segment_id,
            subchain_id,
            label_entity_id,
        )
        .map(|inner| Self { inner })
    }
    #[getter]
    /// Author residue sequence identifier, or None if absent.
    fn seq_id(&self) -> Option<PdbSeqId> {
        self.inner.seq_id().map(|inner| PdbSeqId { inner })
    }
    #[getter]
    /// mmCIF label sequence identifier, or None if absent.
    fn label_seq_id(&self) -> Option<i32> {
        self.inner.label_seq_id()
    }
    #[getter]
    /// Raw source segment identifier bytes, or None if absent.
    fn segment_id<'py>(&self, py: Python<'py>) -> Option<Bound<'py, pyo3::types::PyBytes>> {
        self.inner
            .segment_id()
            .map(|bytes| pyo3::types::PyBytes::new(py, bytes))
    }
    #[getter]
    /// Source subchain identifier, or None if absent.
    fn subchain_id(&self) -> Option<&str> {
        self.inner.subchain_id()
    }
    #[getter]
    /// mmCIF label entity identifier, or None if absent.
    fn label_entity_id(&self) -> Option<&str> {
        self.inner.label_entity_id()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Entity identifier retained from the source structure file.
    pub struct EntitySourceIds {
        pub(crate) inner: ck::EntitySourceIds,
    }
}
#[pymethods]
impl EntitySourceIds {
    #[new]
    fn new(source_entity_id: String) -> Self {
        Self {
            inner: ck::EntitySourceIds::new(source_entity_id),
        }
    }
    #[getter]
    /// Entity identifier used by the original structure file.
    fn source_entity_id(&self) -> &str {
        self.inner.source_entity_id()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Explicit alternate-location selection: any label or an exact optional label.
    pub struct AltLocRequest {
        pub(crate) inner: ck::AltLocRequest,
    }
}
#[pymethods]
impl AltLocRequest {
    #[staticmethod]
    /// Select any alternate location.
    fn any() -> Self {
        Self {
            inner: ck::AltLocRequest::Any,
        }
    }
    #[staticmethod]
    /// Select exactly the supplied alternate location; None selects the unlabelled location.
    fn exact(label: Option<AltLocLabel>) -> Self {
        Self {
            inner: ck::AltLocRequest::Exact(label.map(|value| value.inner)),
        }
    }
}
