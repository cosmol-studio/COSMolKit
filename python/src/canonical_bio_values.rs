//! Detached BIO values projected without parsing or reconstruction.
use crate::canonical_bio_binding::PdbSeqId;
use cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Contiguous half-open span of structural row IDs, represented by start and length.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioRowSpan {
    inner: ck::BioRowSpan<ck::BioAtomId>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioRowSpan {
    /// Construct a BioRowSpan value from the supplied inputs.
    #[new]
    fn new(py: Python<'_>, start: u32, len: u32) -> PyResult<Self> {
        ck::BioRowSpan::new(start, len)
            .map(|inner| Self { inner })
            .map_err(|error| crate::canonical_bio_binding::structure_error(py, &error))
    }
    /// Return the first local row ID in the span.
    fn start(&self) -> u32 {
        self.inner.start()
    }
    /// Return the number of rows in the span.
    fn len(&self) -> u32 {
        self.inner.len()
    }
    /// Return the exclusive ending local row ID of the span.
    fn end(&self) -> u32 {
        self.inner.end()
    }
    /// Return whether there are no stored entries.
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

/// Crystallographic cell with edge lengths a/b/c in angstroms and alpha/beta/gamma angles in degrees.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioCrystalCell {
    pub(crate) inner: ck::BioCrystalCell,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioCrystalCell {
    /// Construct a BioCrystalCell value from the supplied inputs.
    #[new]
    #[pyo3(signature=(a=1.0,b=1.0,c=1.0,alpha=90.0,beta=90.0,gamma=90.0))]
    fn new(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        Self {
            inner: ck::BioCrystalCell {
                a,
                b,
                c,
                alpha,
                beta,
                gamma,
            },
        }
    }
    /// Unit-cell edge length a in angstroms.
    #[getter]
    fn a(&self) -> f64 {
        self.inner.a
    }
    /// Unit-cell edge length b in angstroms.
    #[getter]
    fn b(&self) -> f64 {
        self.inner.b
    }
    /// Unit-cell edge length c in angstroms.
    #[getter]
    fn c(&self) -> f64 {
        self.inner.c
    }
    /// Unit-cell angle alpha in degrees.
    #[getter]
    fn alpha(&self) -> f64 {
        self.inner.alpha
    }
    /// Unit-cell angle beta in degrees.
    #[getter]
    fn beta(&self) -> f64 {
        self.inner.beta
    }
    /// Unit-cell angle gamma in degrees.
    #[getter]
    fn gamma(&self) -> f64 {
        self.inner.gamma
    }
}
/// SIFTS mapping from a structural residue to a UniProt accession index and residue number.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioSiftsUnpResidue {
    pub(crate) inner: ck::BioSiftsUnpResidue,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioSiftsUnpResidue {
    /// Construct a BioSiftsUnpResidue value from the supplied inputs.
    #[new]
    #[pyo3(signature=(residue=None, accession_index=0, number=0))]
    fn new(residue: Option<u8>, accession_index: u8, number: u16) -> Self {
        Self {
            inner: ck::BioSiftsUnpResidue::new(residue, accession_index, number),
        }
    }
    /// Return the parent residue reference.
    #[getter]
    fn residue(&self) -> Option<u8> {
        self.inner.residue()
    }
    /// Index of the UniProt accession referenced by this SIFTS mapping.
    #[getter]
    fn accession_index(&self) -> u8 {
        self.inner.accession_index()
    }
    /// UniProt residue number retained by this SIFTS mapping.
    #[getter]
    fn number(&self) -> u16 {
        self.inner.number()
    }
}

/// External sequence/database reference and its source/database residue ranges.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioEntityDbRef {
    pub(crate) inner: ck::BioEntityDbRef,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioEntityDbRef {
    /// Construct a BioEntityDbRef value from the supplied inputs.
    #[new]
    #[pyo3(signature=(*,db_name="",accession_code="",id_code="",isoform="",seq_begin=None,seq_end=None,db_begin=None,db_end=None,label_seq_begin=None,label_seq_end=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        db_name: &str,
        accession_code: &str,
        id_code: &str,
        isoform: &str,
        seq_begin: Option<&PdbSeqId>,
        seq_end: Option<&PdbSeqId>,
        db_begin: Option<&PdbSeqId>,
        db_end: Option<&PdbSeqId>,
        label_seq_begin: Option<i32>,
        label_seq_end: Option<i32>,
    ) -> Self {
        Self {
            inner: ck::BioEntityDbRef {
                db_name: db_name.to_owned(),
                accession_code: accession_code.to_owned(),
                id_code: id_code.to_owned(),
                isoform: isoform.to_owned(),
                seq_begin: seq_begin.map(|v| v.inner).unwrap_or_default(),
                seq_end: seq_end.map(|v| v.inner).unwrap_or_default(),
                db_begin: db_begin.map(|v| v.inner).unwrap_or_default(),
                db_end: db_end.map(|v| v.inner).unwrap_or_default(),
                label_seq_begin,
                label_seq_end,
            },
        }
    }
    /// External database name.
    #[getter]
    fn db_name(&self) -> &str {
        &self.inner.db_name
    }
    /// External database accession code.
    #[getter]
    fn accession_code(&self) -> &str {
        &self.inner.accession_code
    }
    /// External database entry identifier.
    #[getter]
    fn id_code(&self) -> &str {
        &self.inner.id_code
    }
    /// External sequence isoform identifier.
    #[getter]
    fn isoform(&self) -> &str {
        &self.inner.isoform
    }
    /// Beginning source sequence position for this database mapping.
    #[getter]
    fn seq_begin(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.seq_begin,
        }
    }
    /// Ending source sequence position for this database mapping.
    #[getter]
    fn seq_end(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.seq_end,
        }
    }
    /// Beginning external database sequence position.
    #[getter]
    fn db_begin(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.db_begin,
        }
    }
    /// Ending external database sequence position.
    #[getter]
    fn db_end(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.db_end,
        }
    }
    /// Beginning mmCIF label sequence position for this mapping.
    #[getter]
    fn label_seq_begin(&self) -> Option<i32> {
        self.inner.label_seq_begin
    }
    /// Ending mmCIF label sequence position for this mapping.
    #[getter]
    fn label_seq_end(&self) -> Option<i32> {
        self.inner.label_seq_end
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<BioRowSpan>()?;
    module.add_class::<BioCrystalCell>()?;
    module.add_class::<BioSiftsUnpResidue>()?;
    module.add_class::<BioEntityDbRef>()?;
    Ok(())
}
