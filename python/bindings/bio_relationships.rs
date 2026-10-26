use crate::chemical_values::{EnumValue, enum_value};
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(BioConnectionKind, { Covale, Disulf, Hydrog, MetalC, Unknown });
enum_value!(BioAsu, { Same, Different, Any });
enum_value!(BioHelixClass, { UnknownHelix, RAlpha, ROmega, RPi, RGamma, R310,
    LAlpha, LOmega, LGamma, Helix27, HelixPolyProlineNone });

#[derive(FromPyObject, IntoPyObject)]
pub enum BioConnectionKindInput {
    Enum(BioConnectionKind),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum BioAsuInput {
    Enum(BioAsu),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum BioHelixClassInput {
    Enum(BioHelixClass),
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
enum_input!(BioConnectionKindInput, BioConnectionKind);
enum_input!(BioAsuInput, BioAsu);
enum_input!(BioHelixClassInput, BioHelixClass);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Source residue address comprising numbering, insertion code, segment and residue name.
    pub struct ResidueAddress {
        pub(crate) inner: ck::ResidueAddress,
    }
}
#[pymethods]
impl ResidueAddress {
    #[staticmethod]
    /// Create a source residue address; return None if its residue name or address fields are invalid.
    fn new(
        sequence_number: Option<i32>,
        insertion_code: Option<u8>,
        segment: &str,
        name: &str,
    ) -> Option<Self> {
        let name = ck::ResidueName::from_ascii(name.as_bytes())?;
        ck::ResidueAddress::new(sequence_number, insertion_code, segment.as_bytes(), name)
            .map(|inner| Self { inner })
    }
    #[getter]
    /// Author residue sequence number, or None when unspecified.
    fn sequence_number(&self) -> Option<i32> {
        self.inner.sequence_number()
    }
    #[getter]
    /// Optional source insertion-code byte.
    fn insertion_code(&self) -> Option<u8> {
        self.inner.insertion_code()
    }
    #[getter]
    /// Source segment identifier.
    fn segment(&self) -> &str {
        self.inner.segment()
    }
    #[getter]
    /// Name or identifier retained in this source record.
    fn name(&self) -> String {
        self.inner.name().as_str().to_owned()
    }
    /// Compare residue addresses using the native matching rules, including the segment.
    fn matches(&self, other: &Self) -> bool {
        self.inner.matches(&other.inner)
    }
    /// Compare residue addresses using native matching rules while ignoring the segment.
    fn matches_without_segment(&self, other: &Self) -> bool {
        self.inner.matches_without_segment(&other.inner)
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Source atom address comprising chain, residue, logical atom name and alternate location.
    pub struct AtomAddress {
        pub(crate) inner: ck::AtomAddress,
    }
}
#[pymethods]
impl AtomAddress {
    #[new]
    #[pyo3(signature=(chain_name, residue, logical_atom_name, altloc=None))]
    fn new(
        chain_name: &str,
        residue: ResidueAddress,
        logical_atom_name: String,
        altloc: Option<u8>,
    ) -> PyResult<Self> {
        let chain = ck::PdbChainId::from_ascii(chain_name.as_bytes()).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("chain_name is not a valid PdbChainId")
        })?;
        Ok(Self {
            inner: ck::AtomAddress::new(chain, residue.inner, logical_atom_name, altloc),
        })
    }
    #[getter]
    /// Source chain identifier.
    fn chain_name(&self) -> String {
        self.inner.chain_name().as_str().to_owned()
    }
    #[getter]
    /// Residue address containing the referenced atom.
    fn residue(&self) -> ResidueAddress {
        ResidueAddress {
            inner: self.inner.residue(),
        }
    }
    #[getter]
    /// Logical atom name used to resolve this address.
    fn logical_atom_name(&self) -> &str {
        self.inner.logical_atom_name()
    }
    #[getter]
    /// Native alternate-location label byte.
    fn altloc(&self) -> u8 {
        self.inner.altloc()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Reported inter-atom connection with partner addresses and symmetry metadata.
    pub struct BioConnection {
        pub(crate) inner: ck::BioConnection,
    }
}
#[pymethods]
impl BioConnection {
    #[new]
    #[pyo3(signature=(*, name=None, link_id=None, kind=None, asu=None, partner1=None, partner2=None, reported_distance=None, reported_sym=None))]
    fn new(
        name: Option<String>,
        link_id: Option<String>,
        kind: Option<BioConnectionKindInput>,
        asu: Option<BioAsuInput>,
        partner1: Option<AtomAddress>,
        partner2: Option<AtomAddress>,
        reported_distance: Option<f64>,
        reported_sym: Option<[i16; 4]>,
    ) -> PyResult<Self> {
        let defaults = ck::BioConnection::default();
        Ok(Self {
            inner: ck::BioConnection {
                name: name.unwrap_or(defaults.name),
                link_id: link_id.unwrap_or(defaults.link_id),
                kind: kind
                    .map(BioConnectionKindInput::into_rust)
                    .transpose()?
                    .unwrap_or(defaults.kind),
                asu: asu
                    .map(BioAsuInput::into_rust)
                    .transpose()?
                    .unwrap_or(defaults.asu),
                partner1: partner1
                    .map(|value| value.inner)
                    .unwrap_or(defaults.partner1),
                partner2: partner2
                    .map(|value| value.inner)
                    .unwrap_or(defaults.partner2),
                reported_distance: reported_distance.unwrap_or(defaults.reported_distance),
                reported_sym: reported_sym.unwrap_or(defaults.reported_sym),
            },
        })
    }
    #[getter]
    /// Name or identifier retained in this source record.
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    /// Chemical link identifier reported by the source.
    fn link_id(&self) -> &str {
        &self.inner.link_id
    }
    #[getter]
    /// Native connection classification.
    fn kind(&self) -> BioConnectionKind {
        self.inner.kind.into()
    }
    #[getter]
    /// Relationship of the connection partners to the asymmetric unit.
    fn asu(&self) -> BioAsu {
        self.inner.asu.into()
    }
    #[getter]
    /// Address of the first connection partner.
    fn partner1(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner1.clone(),
        }
    }
    #[getter]
    /// Address of the second connection partner.
    fn partner2(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner2.clone(),
        }
    }
    #[getter]
    /// Inter-partner distance reported in the source.
    fn reported_distance(&self) -> f64 {
        self.inner.reported_distance
    }
    #[getter]
    /// Four-component symmetry identifier retained for the reported connection.
    fn reported_sym(&self) -> [i16; 4] {
        self.inner.reported_sym
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Cis-peptide record with partner addresses, model and reported angle.
    pub struct BioCisPep {
        pub(crate) inner: ck::BioCisPep,
    }
}
#[pymethods]
impl BioCisPep {
    #[new]
    #[pyo3(signature=(*, partner_c=None, partner_n=None, model_num=None, only_altloc=None, reported_angle=None))]
    fn new(
        partner_c: Option<AtomAddress>,
        partner_n: Option<AtomAddress>,
        model_num: Option<i32>,
        only_altloc: Option<u8>,
        reported_angle: Option<f64>,
    ) -> Self {
        let defaults = ck::BioCisPep::default();
        Self {
            inner: ck::BioCisPep {
                partner_c: partner_c
                    .map(|value| value.inner)
                    .unwrap_or(defaults.partner_c),
                partner_n: partner_n
                    .map(|value| value.inner)
                    .unwrap_or(defaults.partner_n),
                model_num: model_num.unwrap_or(defaults.model_num),
                only_altloc: only_altloc.unwrap_or(defaults.only_altloc),
                reported_angle: reported_angle.unwrap_or(defaults.reported_angle),
            },
        }
    }
    #[getter]
    /// Address of the carbon-side cis-peptide partner.
    fn partner_c(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner_c.clone(),
        }
    }
    #[getter]
    /// Address of the nitrogen-side cis-peptide partner.
    fn partner_n(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner_n.clone(),
        }
    }
    #[getter]
    /// Source model number associated with the cis-peptide record.
    fn model_num(&self) -> i32 {
        self.inner.model_num
    }
    #[getter]
    /// Alternate-location restriction retained in the cis-peptide record.
    fn only_altloc(&self) -> u8 {
        self.inner.only_altloc
    }
    #[getter]
    /// Peptide angle reported by the source, in degrees.
    fn reported_angle(&self) -> f64 {
        self.inner.reported_angle
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Modified-residue record linking a source residue to its parent chemical component.
    pub struct BioModRes {
        pub(crate) inner: ck::BioModRes,
    }
}
#[pymethods]
impl BioModRes {
    #[new]
    #[pyo3(signature=(*, chain_name=None, res_id=None, parent_comp_id=None, mod_id=None, details=None))]
    fn new(
        chain_name: Option<&str>,
        res_id: Option<ResidueAddress>,
        parent_comp_id: Option<String>,
        mod_id: Option<String>,
        details: Option<String>,
    ) -> PyResult<Self> {
        let defaults = ck::BioModRes::default();
        let chain_name = chain_name
            .map(|value| {
                ck::PdbChainId::from_ascii(value.as_bytes()).ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err("chain_name is not a valid PdbChainId")
                })
            })
            .transpose()?
            .unwrap_or(defaults.chain_name);
        Ok(Self {
            inner: ck::BioModRes {
                chain_name,
                res_id: res_id.map(|value| value.inner).unwrap_or(defaults.res_id),
                parent_comp_id: parent_comp_id.unwrap_or(defaults.parent_comp_id),
                mod_id: mod_id.unwrap_or(defaults.mod_id),
                details: details.unwrap_or(defaults.details),
            },
        })
    }
    #[getter]
    /// Source chain identifier.
    fn chain_name(&self) -> &str {
        self.inner.chain_name.as_str()
    }
    #[getter]
    /// Address of the modified residue.
    fn res_id(&self) -> ResidueAddress {
        ResidueAddress {
            inner: self.inner.res_id,
        }
    }
    #[getter]
    /// Parent chemical-component identifier for the modified residue.
    fn parent_comp_id(&self) -> &str {
        &self.inner.parent_comp_id
    }
    #[getter]
    /// Modification identifier.
    fn mod_id(&self) -> &str {
        &self.inner.mod_id
    }
    #[getter]
    /// Source explanation of the residue modification.
    fn details(&self) -> &str {
        &self.inner.details
    }
}

#[pyclass(module = "cosmolkit", from_py_object)]
#[derive(Clone)]
/// Reported helical secondary-structure segment and PDB helix classification.
pub struct BioHelix {
    pub(crate) inner: ck::BioHelix,
}
#[pymethods]
impl BioHelix {
    #[new]
    #[pyo3(signature=(*, start=None, end=None, pdb_helix_class=None, length=None))]
    fn new(
        start: Option<AtomAddress>,
        end: Option<AtomAddress>,
        pdb_helix_class: Option<BioHelixClassInput>,
        length: Option<i32>,
    ) -> PyResult<Self> {
        let defaults = ck::BioHelix::default();
        Ok(Self {
            inner: ck::BioHelix {
                start: start.map(|value| value.inner).unwrap_or(defaults.start),
                end: end.map(|value| value.inner).unwrap_or(defaults.end),
                pdb_helix_class: pdb_helix_class
                    .map(BioHelixClassInput::into_rust)
                    .transpose()?
                    .unwrap_or(defaults.pdb_helix_class),
                length: length.unwrap_or(defaults.length),
            },
        })
    }
    #[getter]
    /// Address marking the beginning of this secondary-structure segment.
    fn start(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.start.clone(),
        }
    }
    #[getter]
    /// Address marking the end of this secondary-structure segment.
    fn end(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.end.clone(),
        }
    }
    #[getter]
    /// PDB helix classification.
    fn pdb_helix_class(&self) -> BioHelixClass {
        self.inner.pdb_helix_class.into()
    }
    #[getter]
    /// Reported helix length in residues.
    fn length(&self) -> i32 {
        self.inner.length
    }
    /// Assign the helix class using the native PDB integer-code conversion.
    fn set_helix_class_as_int(&mut self, n: i32) {
        self.inner.set_helix_class_as_int(n);
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Reported beta strand with its boundaries, orientation and hydrogen-bond anchors.
    pub struct BioStrand {
        pub(crate) inner: ck::BioStrand,
    }
}
#[pymethods]
impl BioStrand {
    #[new]
    fn new(
        start: AtomAddress,
        end: AtomAddress,
        hbond_atom2: AtomAddress,
        hbond_atom1: AtomAddress,
        sense: i32,
        name: String,
    ) -> Self {
        Self {
            inner: ck::BioStrand::new(
                start.inner,
                end.inner,
                hbond_atom2.inner,
                hbond_atom1.inner,
                sense,
                name,
            ),
        }
    }
    #[getter]
    /// Address marking the beginning of this secondary-structure segment.
    fn start(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.start.clone(),
        }
    }
    #[getter]
    /// Address marking the end of this secondary-structure segment.
    fn end(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.end.clone(),
        }
    }
    #[getter]
    /// Second atom address in the reported sheet hydrogen bond.
    fn hbond_atom2(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.hbond_atom2.clone(),
        }
    }
    #[getter]
    /// First atom address in the reported sheet hydrogen bond.
    fn hbond_atom1(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.hbond_atom1.clone(),
        }
    }
    #[getter]
    /// Strand orientation relative to its predecessor: native parallel/antiparallel sense code.
    fn sense(&self) -> i32 {
        self.inner.sense
    }
    #[getter]
    /// Name or identifier retained in this source record.
    fn name(&self) -> &str {
        &self.inner.name
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Named beta sheet containing ordered strand records.
    pub struct BioSheet {
        pub(crate) inner: ck::BioSheet,
    }
}
#[pymethods]
impl BioSheet {
    #[new]
    #[pyo3(signature=(sheet_id, strands=None))]
    fn new(sheet_id: &str, strands: Option<Vec<BioStrand>>) -> Self {
        let mut inner = ck::BioSheet::new(sheet_id);
        if let Some(strands) = strands {
            inner.strands = strands.into_iter().map(|value| value.inner).collect();
        }
        Self { inner }
    }
    #[getter]
    /// Name or identifier retained in this source record.
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    /// Independent strand records belonging to this sheet.
    fn strands(&self) -> Vec<BioStrand> {
        self.inner
            .strands
            .iter()
            .cloned()
            .map(|inner| BioStrand { inner })
            .collect()
    }
}
