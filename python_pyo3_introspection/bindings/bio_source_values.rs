use crate::chemical_values::enum_value;
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(EntityKind, { Unknown, Polymer, NonPolymer, Branched, Water });
enum_value!(PolymerKind, { Unknown, PeptideL, PeptideD, Dna, Rna, DnaRnaHybrid,
    SaccharideD, SaccharideL, Pna, CyclicPseudoPeptide, Other });
enum_value!(ResidueKind, { AminoAcid, Dna, Rna, Saccharide, Water, Buffer, Ligand, Unknown });

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone, Copy)]
pub struct PdbSeqId {
    pub(crate) inner: ck::PdbSeqId,
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
    fn seq_num(&self) -> i32 {
        self.inner.seq_num()
    }
    fn ins_code(&self) -> Option<char> {
        self.inner.ins_code().map(char::from)
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone, Copy)]
pub struct AltLocLabel {
    pub(crate) inner: ck::AltLocLabel,
}
#[pymethods]
impl AltLocLabel {
    #[new]
    fn new(value: u8) -> Self {
        Self {
            inner: ck::AltLocLabel::new(value),
        }
    }
    fn value(&self) -> u8 {
        self.inner.value()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone, Copy)]
pub struct AtomSourceIds {
    pub(crate) inner: ck::AtomSourceIds,
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
    fn serial(&self) -> Option<i32> {
        self.inner.serial().map(|value| value.value())
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct ResidueSourceIds {
    pub(crate) inner: ck::ResidueSourceIds,
}
#[pymethods]
impl ResidueSourceIds {
    #[staticmethod]
    #[pyo3(signature=(seq_id=None, label_seq_id=None, segment_id=None, subchain_id=None, label_entity_id=None))]
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
    fn seq_id(&self) -> Option<PdbSeqId> {
        self.inner.seq_id().map(|inner| PdbSeqId { inner })
    }
    fn label_seq_id(&self) -> Option<i32> {
        self.inner.label_seq_id()
    }
    fn segment_id<'py>(&self, py: Python<'py>) -> Option<Bound<'py, pyo3::types::PyBytes>> {
        self.inner
            .segment_id()
            .map(|bytes| pyo3::types::PyBytes::new(py, bytes))
    }
    fn subchain_id(&self) -> Option<&str> {
        self.inner.subchain_id()
    }
    fn label_entity_id(&self) -> Option<&str> {
        self.inner.label_entity_id()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct EntitySourceIds {
    pub(crate) inner: ck::EntitySourceIds,
}
#[pymethods]
impl EntitySourceIds {
    #[new]
    fn new(source_entity_id: String) -> Self {
        Self {
            inner: ck::EntitySourceIds::new(source_entity_id),
        }
    }
    fn source_entity_id(&self) -> &str {
        self.inner.source_entity_id()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone, Copy)]
pub struct AltLocRequest {
    pub(crate) inner: ck::AltLocRequest,
}
#[pymethods]
impl AltLocRequest {
    #[staticmethod]
    fn any() -> Self {
        Self {
            inner: ck::AltLocRequest::Any,
        }
    }
    #[staticmethod]
    fn exact(label: Option<AltLocLabel>) -> Self {
        Self {
            inner: ck::AltLocRequest::Exact(label.map(|value| value.inner)),
        }
    }
}
