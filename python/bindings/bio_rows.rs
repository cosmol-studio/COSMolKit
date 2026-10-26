use crate::bio_source_values::{
    AltLocLabel, AtomSourceIds, EntityKind, EntitySourceIds, PdbSeqId, PolymerKind, ResidueKind,
    ResidueSourceIds,
};
use crate::bio_values::{AtomName, BioRowSpan};
use crate::chemical_values::{EnumValue, enum_value};
use crate::elements::Element;
use crate::residue_values::{ResidueInfoKind, ResidueInfoKindInput};
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(BioCalcFlag, { NotSet, NoHydrogen, Determined, Calculated, Dummy });
// PyO3's tuple-enum IntoPyObject derive must expand on a normal Rust item;
// placing these derive inputs inside macro_rules currently breaks hygiene.
#[derive(FromPyObject, IntoPyObject)]
pub enum BioCalcFlagInput {
    Enum(BioCalcFlag),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum EntityKindInput {
    Enum(EntityKind),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum PolymerKindInput {
    Enum(PolymerKind),
    Name(String),
}
macro_rules! enum_input {
    ($input:ident, $value:ident) => {
        impl $input {
            fn into_rust(self) -> PyResult<ck::$value> {
                match self {
                    Self::Enum(value) => Ok(value.into()),
                    Self::Name(name) => $value::from_name(&name).map(Into::into),
                }
            }
        }
    };
}
enum_input!(BioCalcFlagInput, BioCalcFlag);
enum_input!(EntityKindInput, EntityKind);
enum_input!(PolymerKindInput, PolymerKind);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// SIFTS mapping of a structure residue to a UniProt accession index and residue number.
    pub struct BioSiftsUnpResidue {
        pub(crate) inner: ck::BioSiftsUnpResidue,
    }
}
#[pymethods]
impl BioSiftsUnpResidue {
    #[new]
    #[pyo3(signature=(*, residue=ck::BioSiftsUnpResidue::default().residue(), accession_index=ck::BioSiftsUnpResidue::default().accession_index(), number=ck::BioSiftsUnpResidue::default().number()))]
    fn new(residue: Option<u8>, accession_index: u8, number: u16) -> Self {
        Self {
            inner: ck::BioSiftsUnpResidue::new(residue, accession_index, number),
        }
    }
    #[getter]
    /// Optional residue character encoded by the SIFTS mapping.
    fn residue(&self) -> Option<u8> {
        self.inner.residue()
    }
    #[getter]
    /// Index into the entity's SIFTS UniProt accession list.
    fn accession_index(&self) -> u8 {
        self.inner.accession_index()
    }
    #[getter]
    /// Mapped UniProt residue number.
    fn number(&self) -> u16 {
        self.inner.number()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Source sequence range and corresponding external database accession and range.
    pub struct BioEntityDbRef {
        pub(crate) inner: ck::BioEntityDbRef,
    }
}
#[pymethods]
impl BioEntityDbRef {
    #[new]
    #[pyo3(signature=(*, db_name=ck::BioEntityDbRef::default().db_name, accession_code=ck::BioEntityDbRef::default().accession_code, id_code=ck::BioEntityDbRef::default().id_code, isoform=ck::BioEntityDbRef::default().isoform, seq_begin=None, seq_end=None, db_begin=None, db_end=None, label_seq_begin=None, label_seq_end=None))]
    fn new(
        db_name: String,
        accession_code: String,
        id_code: String,
        isoform: String,
        seq_begin: Option<PdbSeqId>,
        seq_end: Option<PdbSeqId>,
        db_begin: Option<PdbSeqId>,
        db_end: Option<PdbSeqId>,
        label_seq_begin: Option<i32>,
        label_seq_end: Option<i32>,
    ) -> Self {
        Self {
            inner: ck::BioEntityDbRef {
                db_name,
                accession_code,
                id_code,
                isoform,
                seq_begin: seq_begin
                    .map(|value| value.inner)
                    .unwrap_or_else(|| ck::BioEntityDbRef::default().seq_begin),
                seq_end: seq_end
                    .map(|value| value.inner)
                    .unwrap_or_else(|| ck::BioEntityDbRef::default().seq_end),
                db_begin: db_begin
                    .map(|value| value.inner)
                    .unwrap_or_else(|| ck::BioEntityDbRef::default().db_begin),
                db_end: db_end
                    .map(|value| value.inner)
                    .unwrap_or_else(|| ck::BioEntityDbRef::default().db_end),
                label_seq_begin,
                label_seq_end,
            },
        }
    }
    #[getter]
    /// Name of the external sequence database.
    fn db_name(&self) -> &str {
        &self.inner.db_name
    }
    #[getter]
    /// Accession code in the external database.
    fn accession_code(&self) -> &str {
        &self.inner.accession_code
    }
    #[getter]
    /// External database entry identifier.
    fn id_code(&self) -> &str {
        &self.inner.id_code
    }
    #[getter]
    /// External database isoform identifier.
    fn isoform(&self) -> &str {
        &self.inner.isoform
    }
    #[getter]
    /// First source residue identifier in the database mapping.
    fn seq_begin(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.seq_begin,
        }
    }
    #[getter]
    /// Last source residue identifier in the database mapping.
    fn seq_end(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.seq_end,
        }
    }
    #[getter]
    /// First database residue identifier in the mapping.
    fn db_begin(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.db_begin,
        }
    }
    #[getter]
    /// Last database residue identifier in the mapping.
    fn db_end(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.db_end,
        }
    }
    #[getter]
    /// First mmCIF label sequence number in the mapping, or None.
    fn label_seq_begin(&self) -> Option<i32> {
        self.inner.label_seq_begin
    }
    #[getter]
    /// Last mmCIF label sequence number in the mapping, or None.
    fn label_seq_end(&self) -> Option<i32> {
        self.inner.label_seq_end
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Detached biological entity metadata, sequence and database references; no owning structure.
    pub struct BioEntityRow {
        pub(crate) inner: ck::BioEntityRow,
    }
}
#[pymethods]
impl BioEntityRow {
    #[new]
    fn new(
        kind: EntityKindInput,
        polymer_kind: PolymerKindInput,
        reflects_microhetero: bool,
        full_sequence: Vec<String>,
        dbrefs: Vec<BioEntityDbRef>,
        sifts_unp_accessions: Vec<String>,
        subchains: Vec<String>,
        source: &EntitySourceIds,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BioEntityRow::new(
                kind.into_rust()?,
                polymer_kind.into_rust()?,
                reflects_microhetero,
                full_sequence,
                dbrefs.into_iter().map(|value| value.inner).collect(),
                sifts_unp_accessions,
                subchains,
                source.inner.clone(),
            ),
        })
    }
    #[staticmethod]
    /// Return the first monomer name from a comma-separated microheterogeneity list.
    fn first_mon(mon_list: &str) -> &str {
        ck::BioEntityRow::first_mon(mon_list)
    }
    #[getter]
    /// Stored native classification of this entity or residue.
    fn kind(&self) -> EntityKind {
        self.inner.kind().into()
    }
    #[getter]
    /// Polymer classification stored for the entity.
    fn polymer_kind(&self) -> PolymerKind {
        self.inner.polymer_kind().into()
    }
    #[getter]
    /// Whether the sequence records microheterogeneity.
    fn reflects_microhetero(&self) -> bool {
        self.inner.reflects_microhetero()
    }
    #[getter]
    /// Independent list of source monomer names for the complete sequence.
    fn full_sequence(&self) -> Vec<String> {
        self.inner.full_sequence().to_vec()
    }
    #[getter]
    /// Independent list of external sequence database references.
    fn dbrefs(&self) -> Vec<BioEntityDbRef> {
        self.inner
            .dbrefs()
            .iter()
            .cloned()
            .map(|inner| BioEntityDbRef { inner })
            .collect()
    }
    #[getter]
    /// UniProt accessions referenced by residue SIFTS mappings.
    fn sifts_unp_accessions(&self) -> Vec<String> {
        self.inner.sifts_unp_accessions().to_vec()
    }
    #[getter]
    /// Source subchain identifiers associated with this entity.
    fn subchains(&self) -> Vec<String> {
        self.inner.subchains().to_vec()
    }
    #[getter]
    /// Source-format identifiers, distinct from local structure row indices.
    fn source(&self) -> EntitySourceIds {
        EntitySourceIds {
            inner: self.inner.source().clone(),
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Detached residue row with a parent chain identifier and atom span; no owning structure.
    pub struct BioResidueRow {
        pub(crate) inner: ck::BioResidueRow,
    }
}
#[pymethods]
impl BioResidueRow {
    #[new]
    #[pyo3(signature=(chain_id, atom_span, name, residue_info_kind, entity_kind, entity_id, het_flag, source, sifts_unp))]
    fn new(
        chain_id: u32,
        atom_span: &BioRowSpan,
        name: &str,
        residue_info_kind: ResidueInfoKindInput,
        entity_kind: EntityKindInput,
        entity_id: Option<u32>,
        het_flag: Option<u8>,
        source: &ResidueSourceIds,
        sifts_unp: BioSiftsUnpResidue,
    ) -> PyResult<Self> {
        let name = ck::ResidueName::from_ascii(name.as_bytes()).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("name is not a valid ResidueName")
        })?;
        let span = atom_span
            .as_rust()
            .map_err(crate::bio_errors::BioStructureError::from_rust)?;
        Ok(Self {
            inner: ck::BioResidueRow::new(
                ck::BioChainId::new(chain_id),
                span,
                name,
                residue_info_kind.into_rust()?,
                entity_kind.into_rust()?,
                entity_id.map(ck::BioEntityId::new),
                het_flag,
                source.inner.clone(),
                sifts_unp.inner,
            ),
        })
    }
    #[getter]
    /// Local identifier of the parent chain.
    fn chain_id(&self) -> u32 {
        self.inner.chain_id().value()
    }
    #[getter]
    /// Half-open span of atom rows belonging to this residue.
    fn atom_span(&self) -> PyResult<BioRowSpan> {
        BioRowSpan::from_rust(self.inner.atom_span())
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    /// Source residue or atom name, preserving its native representation.
    fn name(&self) -> String {
        self.inner.name().as_str().to_owned()
    }
    #[getter]
    /// Residue-information classification stored in the row.
    fn residue_info_kind(&self) -> ResidueInfoKind {
        self.inner.residue_info_kind().into()
    }
    #[getter]
    /// Stored native classification of this entity or residue.
    fn kind(&self) -> ResidueKind {
        self.inner.kind().into()
    }
    #[getter]
    /// Entity classification associated with this residue.
    fn entity_kind(&self) -> EntityKind {
        self.inner.entity_kind().into()
    }
    #[getter]
    /// Parent entity identifier, or None if unassigned.
    fn entity_id(&self) -> Option<u32> {
        self.inner.entity_id().map(|value| value.value())
    }
    #[getter]
    /// Source heteroatom-record flag, or None if absent.
    fn het_flag(&self) -> Option<u8> {
        self.inner.het_flag()
    }
    #[getter]
    /// Source-format identifiers, distinct from local structure row indices.
    fn source(&self) -> ResidueSourceIds {
        ResidueSourceIds {
            inner: self.inner.source().clone(),
        }
    }
    #[getter]
    /// Retained UniProt residue mapping.
    fn sifts_unp(&self) -> BioSiftsUnpResidue {
        BioSiftsUnpResidue {
            inner: self.inner.sifts_unp(),
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Detached atom metadata row. Cartesian positions belong to the structure's separate coordinate block.
    pub struct BioAtomRow {
        pub(crate) inner: ck::BioAtomRow,
    }
}
#[pymethods]
impl BioAtomRow {
    #[new]
    #[pyo3(signature=(residue_id, name, element, isotope_mass_number, altloc, formal_charge, calc_flag, occupancy, b_iso, anisou, tls_group_id, fraction, source))]
    fn new(
        residue_id: u32,
        name: &AtomName,
        element: &Element,
        isotope_mass_number: Option<u16>,
        altloc: Option<AltLocLabel>,
        formal_charge: i8,
        calc_flag: BioCalcFlagInput,
        occupancy: f64,
        b_iso: f64,
        anisou: [f64; 6],
        tls_group_id: i16,
        fraction: f64,
        source: &AtomSourceIds,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BioAtomRow::new(
                ck::BioResidueId::new(residue_id),
                name.inner,
                element.inner,
                isotope_mass_number,
                altloc.map(|value| value.inner),
                formal_charge,
                calc_flag.into_rust()?,
                occupancy,
                b_iso,
                anisou,
                tls_group_id,
                fraction,
                source.inner,
            ),
        })
    }
    #[getter]
    /// Local identifier of the parent residue.
    fn residue_id(&self) -> u32 {
        self.inner.residue_id().value()
    }
    #[getter]
    /// Source residue or atom name, preserving its native representation.
    fn name(&self) -> AtomName {
        AtomName {
            inner: self.inner.name(),
        }
    }
    #[getter]
    /// Chemical element assigned to this atom.
    fn element(&self) -> Element {
        Element {
            inner: self.inner.element(),
        }
    }
    #[getter]
    /// Explicit isotope mass number, or None for natural isotopic composition.
    fn isotope_mass_number(&self) -> Option<u16> {
        self.inner.isotope_mass_number()
    }
    #[getter]
    /// Alternate-location identifier, or None if unlabelled.
    fn altloc(&self) -> Option<AltLocLabel> {
        self.inner.altloc().map(|inner| AltLocLabel { inner })
    }
    #[getter]
    /// Formal charge stored on the atom.
    fn formal_charge(&self) -> i8 {
        self.inner.formal_charge()
    }
    #[getter]
    /// Source classification of whether the atom position was calculated or determined.
    fn calc_flag(&self) -> BioCalcFlag {
        self.inner.calc_flag().into()
    }
    #[getter]
    /// Occupancy reported for this atom.
    fn occupancy(&self) -> f64 {
        self.inner.occupancy()
    }
    #[getter]
    /// Isotropic atomic displacement parameter.
    fn b_iso(&self) -> f64 {
        self.inner.b_iso()
    }
    #[getter]
    /// Six independent components of the anisotropic displacement tensor.
    fn anisou(&self) -> Vec<f64> {
        self.inner.anisou().to_vec()
    }
    #[getter]
    /// TLS refinement group identifier stored for this atom.
    fn tls_group_id(&self) -> i16 {
        self.inner.tls_group_id()
    }
    #[getter]
    /// Isotopic fraction retained by the native atom row.
    fn fraction(&self) -> f64 {
        self.inner.fraction()
    }
    #[getter]
    /// Source-format identifiers, distinct from local structure row indices.
    fn source(&self) -> AtomSourceIds {
        AtomSourceIds {
            inner: *self.inner.source(),
        }
    }
}
