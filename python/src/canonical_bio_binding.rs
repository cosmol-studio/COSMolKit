//! Shared hierarchy snapshots and thin Python projections of canonical BIO APIs.

use crate::canonical_bio_residue::{ResidueInfo, enum_member};
use crate::canonical_element_metadata::Element;
use ::cosmolkit as ck;
use pyo3::exceptions::{PyIndexError, PyOSError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::path::PathBuf;
use std::sync::Arc;

pyo3::create_exception!(cosmolkit, BioReadError, PyValueError);
pyo3::create_exception!(cosmolkit, BioMoleculeError, PyValueError);
pyo3::create_exception!(cosmolkit, BioMoleculeConversionError, PyValueError);
pyo3::create_exception!(cosmolkit, BioPdbReadError, BioReadError);
pyo3::create_exception!(cosmolkit, BioMmcifReadError, BioReadError);
pyo3::create_exception!(cosmolkit, ProteinReadError, BioReadError);
pyo3::create_exception!(cosmolkit, BioOperationError, PyValueError);
pyo3::create_exception!(cosmolkit, BioMmcifWriteError, PyValueError);
pyo3::create_exception!(cosmolkit, BioSelectionParseError, PyValueError);
pyo3::create_exception!(cosmolkit, BioPdbWriteError, PyValueError);

fn annotate(
    py: Python<'_>,
    error: PyErr,
    kind: &str,
    source: &(dyn std::error::Error + 'static),
) -> PyErr {
    let fields = || -> PyResult<()> {
        error.value(py).setattr("domain", "bio")?;
        error.value(py).setattr("kind", kind)
    };
    if let Err(e) = fields() {
        return e;
    }
    error.set_cause(
        py,
        source
            .source()
            .map(|e| crate::canonical_values::source_pyerr(py, e)),
    );
    error
}
fn read_error(py: Python<'_>, source: ck::BioReadError) -> PyErr {
    let kind = match &source {
        ck::BioReadError::Io { .. } => "Io",
        ck::BioReadError::Utf8 { .. } => "Utf8",
        ck::BioReadError::Pdb(_) => "Pdb",
        ck::BioReadError::Cif(_) => "Cif",
        ck::BioReadError::Mmcif(_) => "Mmcif",
        ck::BioReadError::Mmjson(_) => "Mmjson",
        ck::BioReadError::ChemComp(_) => "ChemComp",
        ck::BioReadError::WrongFormat { .. } => "WrongFormat",
        ck::BioReadError::UnknownFileFormat(_) => "UnknownFileFormat",
    };
    let error = if let ck::BioReadError::Io { path, source: io } = &source {
        PyOSError::new_err((
            io.raw_os_error(),
            source.to_string(),
            path.to_string_lossy().into_owned(),
        ))
    } else {
        BioReadError::new_err(source.to_string())
    };
    annotate(py, error, kind, &source)
}
fn pdb_error(py: Python<'_>, source: ck::BioPdbReadError) -> PyErr {
    let error = annotate(
        py,
        BioPdbReadError::new_err(source.to_string()),
        "Pdb",
        &source,
    );
    let fields = || -> PyResult<()> {
        error
            .value(py)
            .setattr("line_number", source.line_number())?;
        error
            .value(py)
            .setattr("stage", format!("{:?}", source.stage()))?;
        error.value(py).setattr(
            "record_tag",
            source
                .record_tag()
                .map(|x| String::from_utf8_lossy(&x).into_owned()),
        )
    };
    match fields() {
        Ok(()) => error,
        Err(e) => e,
    }
}
fn mmcif_error(py: Python<'_>, source: ck::BioMmcifReadError) -> PyErr {
    let error = annotate(
        py,
        BioMmcifReadError::new_err(source.to_string()),
        "Mmcif",
        &source,
    );
    match error
        .value(py)
        .setattr("stage", format!("{:?}", source.stage()))
    {
        Ok(()) => error,
        Err(e) => e,
    }
}

fn protein_error(py: Python<'_>, source: ck::ProteinReadError) -> PyErr {
    let kind = match &source {
        ck::ProteinReadError::Structure(_) => "Structure",
        ck::ProteinReadError::Pdb(_) => "Pdb",
        ck::ProteinReadError::Mmcif(_) => "Mmcif",
        ck::ProteinReadError::Projection(_) => "Projection",
    };
    let error = annotate(
        py,
        ProteinReadError::new_err(source.to_string()),
        kind,
        &source,
    );
    let cause = match source {
        ck::ProteinReadError::Structure(e) => read_error(py, e),
        ck::ProteinReadError::Pdb(e) => pdb_error(py, e),
        ck::ProteinReadError::Mmcif(e) => mmcif_error(py, e),
        ck::ProteinReadError::Projection(e) => crate::canonical_values::source_pyerr(py, &e),
    };
    error.set_cause(py, Some(cause));
    error
}

fn operation_error(py: Python<'_>, source: ck::BioOperationError) -> PyErr {
    let kind = match &source {
        ck::BioOperationError::Structure(_) => "Structure",
        ck::BioOperationError::Protein(_) => "Protein",
        ck::BioOperationError::Selection(_) => "Selection",
    };
    annotate(
        py,
        BioOperationError::new_err(source.to_string()),
        kind,
        &source,
    )
}
fn write_error(py: Python<'_>, source: ck::BioMmcifWriteError) -> PyErr {
    let (kind, error) = match &source {
        ck::BioMmcifWriteError::FileWrite { path, source: io } => (
            "FileWrite",
            PyOSError::new_err((
                io.raw_os_error(),
                source.to_string(),
                path.to_string_lossy().into_owned(),
            )),
        ),
        ck::BioMmcifWriteError::Io(io) => (
            "Io",
            PyOSError::new_err((io.raw_os_error(), source.to_string())),
        ),
        ck::BioMmcifWriteError::Structure(_) => {
            ("Structure", BioMmcifWriteError::new_err(source.to_string()))
        }
        ck::BioMmcifWriteError::Cif(_) => ("Cif", BioMmcifWriteError::new_err(source.to_string())),
        ck::BioMmcifWriteError::InvalidText { .. } => (
            "InvalidText",
            BioMmcifWriteError::new_err(source.to_string()),
        ),
    };
    annotate(py, error, kind, &source)
}
fn expand_path(py: Python<'_>, path: PathBuf) -> PyResult<PathBuf> {
    py.import("os.path")?
        .getattr("expanduser")?
        .call1((path,))?
        .extract()
}
fn format_from_code(code: u8) -> PyResult<ck::BioCoordinateFormat> {
    use ck::BioCoordinateFormat as F;
    match code {
        0 => Ok(F::Unknown),
        1 => Ok(F::Detect),
        2 => Ok(F::Pdb),
        3 => Ok(F::Mmcif),
        4 => Ok(F::Mmjson),
        5 => Ok(F::ChemComp),
        _ => Err(PyValueError::new_err(format!(
            "invalid BioCoordinateFormat code {code}"
        ))),
    }
}
fn index(index: isize, len: usize, kind: &str) -> PyResult<usize> {
    let n = if index < 0 {
        index.checked_add(len as isize)
    } else {
        Some(index)
    };
    n.filter(|n| *n >= 0 && (*n as usize) < len)
        .map(|n| n as usize)
        .ok_or_else(|| PyIndexError::new_err(format!("{kind} index out of range")))
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioReadParams {
    inner: ck::BioReadParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioReadParams {
    #[new]
    #[pyo3(signature=(*, format=0, source_name="<string>"))]
    fn new(format: u8, source_name: &str) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BioReadParams {
                format: format_from_code(format)?,
                source_name: source_name.to_owned(),
            },
        })
    }
    #[getter]
    #[gen_stub(override_return_type(type_repr = "BioCoordinateFormat"))]
    fn format<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "BioCoordinateFormat", self.inner.format as i64)
    }
    #[getter]
    fn source_name(&self) -> &str {
        &self.inner.source_name
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioPdbReadParams {
    inner: ck::BioPdbReadParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioPdbReadParams {
    #[new]
    #[pyo3(signature=(*, max_line_length=0, check_non_ascii=false, ignore_ter=false, split_chain_on_ter=false, skip_remarks=false))]
    fn new(
        max_line_length: i32,
        check_non_ascii: bool,
        ignore_ter: bool,
        split_chain_on_ter: bool,
        skip_remarks: bool,
    ) -> Self {
        Self {
            inner: ck::BioPdbReadParams {
                max_line_length,
                check_non_ascii,
                ignore_ter,
                split_chain_on_ter,
                skip_remarks,
            },
        }
    }
    #[getter]
    fn max_line_length(&self) -> i32 {
        self.inner.max_line_length
    }
    #[getter]
    fn check_non_ascii(&self) -> bool {
        self.inner.check_non_ascii
    }
    #[getter]
    fn ignore_ter(&self) -> bool {
        self.inner.ignore_ter
    }
    #[getter]
    fn split_chain_on_ter(&self) -> bool {
        self.inner.split_chain_on_ter
    }
    #[getter]
    fn skip_remarks(&self) -> bool {
        self.inner.skip_remarks
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioSelection {
    inner: ck::BioSelection,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioSelection {
    #[staticmethod]
    fn from_cid(py: Python<'_>, cid: &str) -> PyResult<Self> {
        ck::BioSelection::from_cid(cid)
            .map(|inner| Self { inner })
            .map_err(|e| {
                annotate(
                    py,
                    BioSelectionParseError::new_err(e.to_string()),
                    "Parse",
                    &e,
                )
            })
    }
    fn to_cid(&self) -> String {
        self.inner.to_cid()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct BioMmcifWriteParams {
    inner: ck::BioMmcifWriteParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioMmcifWriteParams {
    #[new]
    #[pyo3(signature=(all_groups=true, atoms=None, block_name=None, entry=None, database_status=None, author=None, cell=None, symmetry=None, entity=None, entity_poly=None, struct_ref=None, chem_comp=None, exptl=None, diffrn=None, reflns=None, refine=None, title_keywords=None, ncs=None, struct_asym=None, origx=None, struct_conf=None, struct_sheet=None, struct_biol=None, assembly=None, conn=None, cis=None, modres=None, scale=None, atom_type=None, entity_poly_seq=None, tls=None, software=None, group_pdb=None, auth_all=None, prefer_pairs=false, compact=false, misuse_hash=false, align_pairs=0, align_loops=0))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        all_groups: bool,
        atoms: Option<bool>,
        block_name: Option<bool>,
        entry: Option<bool>,
        database_status: Option<bool>,
        author: Option<bool>,
        cell: Option<bool>,
        symmetry: Option<bool>,
        entity: Option<bool>,
        entity_poly: Option<bool>,
        struct_ref: Option<bool>,
        chem_comp: Option<bool>,
        exptl: Option<bool>,
        diffrn: Option<bool>,
        reflns: Option<bool>,
        refine: Option<bool>,
        title_keywords: Option<bool>,
        ncs: Option<bool>,
        struct_asym: Option<bool>,
        origx: Option<bool>,
        struct_conf: Option<bool>,
        struct_sheet: Option<bool>,
        struct_biol: Option<bool>,
        assembly: Option<bool>,
        conn: Option<bool>,
        cis: Option<bool>,
        modres: Option<bool>,
        scale: Option<bool>,
        atom_type: Option<bool>,
        entity_poly_seq: Option<bool>,
        tls: Option<bool>,
        software: Option<bool>,
        group_pdb: Option<bool>,
        auth_all: Option<bool>,
        prefer_pairs: bool,
        compact: bool,
        misuse_hash: bool,
        align_pairs: u16,
        align_loops: u16,
    ) -> Self {
        Self {
            inner: ck::BioMmcifWriteParams {
                atoms: atoms.unwrap_or(all_groups),
                block_name: block_name.unwrap_or(all_groups),
                entry: entry.unwrap_or(all_groups),
                database_status: database_status.unwrap_or(all_groups),
                author: author.unwrap_or(all_groups),
                cell: cell.unwrap_or(all_groups),
                symmetry: symmetry.unwrap_or(all_groups),
                entity: entity.unwrap_or(all_groups),
                entity_poly: entity_poly.unwrap_or(all_groups),
                struct_ref: struct_ref.unwrap_or(all_groups),
                chem_comp: chem_comp.unwrap_or(all_groups),
                exptl: exptl.unwrap_or(all_groups),
                diffrn: diffrn.unwrap_or(all_groups),
                reflns: reflns.unwrap_or(all_groups),
                refine: refine.unwrap_or(all_groups),
                title_keywords: title_keywords.unwrap_or(all_groups),
                ncs: ncs.unwrap_or(all_groups),
                struct_asym: struct_asym.unwrap_or(all_groups),
                origx: origx.unwrap_or(all_groups),
                struct_conf: struct_conf.unwrap_or(all_groups),
                struct_sheet: struct_sheet.unwrap_or(all_groups),
                struct_biol: struct_biol.unwrap_or(all_groups),
                assembly: assembly.unwrap_or(all_groups),
                conn: conn.unwrap_or(all_groups),
                cis: cis.unwrap_or(all_groups),
                modres: modres.unwrap_or(all_groups),
                scale: scale.unwrap_or(all_groups),
                atom_type: atom_type.unwrap_or(all_groups),
                entity_poly_seq: entity_poly_seq.unwrap_or(all_groups),
                tls: tls.unwrap_or(all_groups),
                software: software.unwrap_or(all_groups),
                group_pdb: group_pdb.unwrap_or(all_groups),
                auth_all: auth_all.unwrap_or(false),
                prefer_pairs,
                compact,
                misuse_hash,
                align_pairs,
                align_loops,
            },
        }
    }
    #[getter]
    fn atoms(&self) -> bool {
        self.inner.atoms
    }
    #[setter]
    fn set_atoms(&mut self, value: bool) {
        self.inner.atoms = value;
    }
    #[getter]
    fn block_name(&self) -> bool {
        self.inner.block_name
    }
    #[setter]
    fn set_block_name(&mut self, value: bool) {
        self.inner.block_name = value;
    }
    #[getter]
    fn entry(&self) -> bool {
        self.inner.entry
    }
    #[setter]
    fn set_entry(&mut self, value: bool) {
        self.inner.entry = value;
    }
    #[getter]
    fn database_status(&self) -> bool {
        self.inner.database_status
    }
    #[setter]
    fn set_database_status(&mut self, value: bool) {
        self.inner.database_status = value;
    }
    #[getter]
    fn author(&self) -> bool {
        self.inner.author
    }
    #[setter]
    fn set_author(&mut self, value: bool) {
        self.inner.author = value;
    }
    #[getter]
    fn cell(&self) -> bool {
        self.inner.cell
    }
    #[setter]
    fn set_cell(&mut self, value: bool) {
        self.inner.cell = value;
    }
    #[getter]
    fn symmetry(&self) -> bool {
        self.inner.symmetry
    }
    #[setter]
    fn set_symmetry(&mut self, value: bool) {
        self.inner.symmetry = value;
    }
    #[getter]
    fn entity(&self) -> bool {
        self.inner.entity
    }
    #[setter]
    fn set_entity(&mut self, value: bool) {
        self.inner.entity = value;
    }
    #[getter]
    fn entity_poly(&self) -> bool {
        self.inner.entity_poly
    }
    #[setter]
    fn set_entity_poly(&mut self, value: bool) {
        self.inner.entity_poly = value;
    }
    #[getter]
    fn struct_ref(&self) -> bool {
        self.inner.struct_ref
    }
    #[setter]
    fn set_struct_ref(&mut self, value: bool) {
        self.inner.struct_ref = value;
    }
    #[getter]
    fn chem_comp(&self) -> bool {
        self.inner.chem_comp
    }
    #[setter]
    fn set_chem_comp(&mut self, value: bool) {
        self.inner.chem_comp = value;
    }
    #[getter]
    fn exptl(&self) -> bool {
        self.inner.exptl
    }
    #[setter]
    fn set_exptl(&mut self, value: bool) {
        self.inner.exptl = value;
    }
    #[getter]
    fn diffrn(&self) -> bool {
        self.inner.diffrn
    }
    #[setter]
    fn set_diffrn(&mut self, value: bool) {
        self.inner.diffrn = value;
    }
    #[getter]
    fn reflns(&self) -> bool {
        self.inner.reflns
    }
    #[setter]
    fn set_reflns(&mut self, value: bool) {
        self.inner.reflns = value;
    }
    #[getter]
    fn refine(&self) -> bool {
        self.inner.refine
    }
    #[setter]
    fn set_refine(&mut self, value: bool) {
        self.inner.refine = value;
    }
    #[getter]
    fn title_keywords(&self) -> bool {
        self.inner.title_keywords
    }
    #[setter]
    fn set_title_keywords(&mut self, value: bool) {
        self.inner.title_keywords = value;
    }
    #[getter]
    fn ncs(&self) -> bool {
        self.inner.ncs
    }
    #[setter]
    fn set_ncs(&mut self, value: bool) {
        self.inner.ncs = value;
    }
    #[getter]
    fn struct_asym(&self) -> bool {
        self.inner.struct_asym
    }
    #[setter]
    fn set_struct_asym(&mut self, value: bool) {
        self.inner.struct_asym = value;
    }
    #[getter]
    fn origx(&self) -> bool {
        self.inner.origx
    }
    #[setter]
    fn set_origx(&mut self, value: bool) {
        self.inner.origx = value;
    }
    #[getter]
    fn struct_conf(&self) -> bool {
        self.inner.struct_conf
    }
    #[setter]
    fn set_struct_conf(&mut self, value: bool) {
        self.inner.struct_conf = value;
    }
    #[getter]
    fn struct_sheet(&self) -> bool {
        self.inner.struct_sheet
    }
    #[setter]
    fn set_struct_sheet(&mut self, value: bool) {
        self.inner.struct_sheet = value;
    }
    #[getter]
    fn struct_biol(&self) -> bool {
        self.inner.struct_biol
    }
    #[setter]
    fn set_struct_biol(&mut self, value: bool) {
        self.inner.struct_biol = value;
    }
    #[getter]
    fn assembly(&self) -> bool {
        self.inner.assembly
    }
    #[setter]
    fn set_assembly(&mut self, value: bool) {
        self.inner.assembly = value;
    }
    #[getter]
    fn conn(&self) -> bool {
        self.inner.conn
    }
    #[setter]
    fn set_conn(&mut self, value: bool) {
        self.inner.conn = value;
    }
    #[getter]
    fn cis(&self) -> bool {
        self.inner.cis
    }
    #[setter]
    fn set_cis(&mut self, value: bool) {
        self.inner.cis = value;
    }
    #[getter]
    fn modres(&self) -> bool {
        self.inner.modres
    }
    #[setter]
    fn set_modres(&mut self, value: bool) {
        self.inner.modres = value;
    }
    #[getter]
    fn scale(&self) -> bool {
        self.inner.scale
    }
    #[setter]
    fn set_scale(&mut self, value: bool) {
        self.inner.scale = value;
    }
    #[getter]
    fn atom_type(&self) -> bool {
        self.inner.atom_type
    }
    #[setter]
    fn set_atom_type(&mut self, value: bool) {
        self.inner.atom_type = value;
    }
    #[getter]
    fn entity_poly_seq(&self) -> bool {
        self.inner.entity_poly_seq
    }
    #[setter]
    fn set_entity_poly_seq(&mut self, value: bool) {
        self.inner.entity_poly_seq = value;
    }
    #[getter]
    fn tls(&self) -> bool {
        self.inner.tls
    }
    #[setter]
    fn set_tls(&mut self, value: bool) {
        self.inner.tls = value;
    }
    #[getter]
    fn software(&self) -> bool {
        self.inner.software
    }
    #[setter]
    fn set_software(&mut self, value: bool) {
        self.inner.software = value;
    }
    #[getter]
    fn group_pdb(&self) -> bool {
        self.inner.group_pdb
    }
    #[setter]
    fn set_group_pdb(&mut self, value: bool) {
        self.inner.group_pdb = value;
    }
    #[getter]
    fn auth_all(&self) -> bool {
        self.inner.auth_all
    }
    #[setter]
    fn set_auth_all(&mut self, value: bool) {
        self.inner.auth_all = value;
    }
    #[getter]
    fn prefer_pairs(&self) -> bool {
        self.inner.prefer_pairs
    }
    #[setter]
    fn set_prefer_pairs(&mut self, value: bool) {
        self.inner.prefer_pairs = value;
    }
    #[getter]
    fn compact(&self) -> bool {
        self.inner.compact
    }
    #[setter]
    fn set_compact(&mut self, value: bool) {
        self.inner.compact = value;
    }
    #[getter]
    fn misuse_hash(&self) -> bool {
        self.inner.misuse_hash
    }
    #[setter]
    fn set_misuse_hash(&mut self, value: bool) {
        self.inner.misuse_hash = value;
    }
    #[getter]
    fn align_pairs(&self) -> u16 {
        self.inner.align_pairs
    }
    #[setter]
    fn set_align_pairs(&mut self, value: u16) {
        self.inner.align_pairs = value;
    }
    #[getter]
    fn align_loops(&self) -> u16 {
        self.inner.align_loops
    }
    #[setter]
    fn set_align_loops(&mut self, value: u16) {
        self.inner.align_loops = value;
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct BioPdbWriteParams {
    inner: ck::BioPdbWriteParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioPdbWriteParams {
    #[new]
    #[pyo3(signature=(ter_records=true, numbered_ter=true, ter_ignores_type=false, preserve_serial=false, end_record=true))]
    fn new(
        ter_records: bool,
        numbered_ter: bool,
        ter_ignores_type: bool,
        preserve_serial: bool,
        end_record: bool,
    ) -> Self {
        Self {
            inner: ck::BioPdbWriteParams {
                ter_records,
                numbered_ter,
                ter_ignores_type,
                preserve_serial,
                end_record,
            },
        }
    }
    #[getter]
    fn ter_records(&self) -> bool {
        self.inner.ter_records
    }
    #[setter]
    fn set_ter_records(&mut self, value: bool) {
        self.inner.ter_records = value;
    }
    #[getter]
    fn numbered_ter(&self) -> bool {
        self.inner.numbered_ter
    }
    #[setter]
    fn set_numbered_ter(&mut self, value: bool) {
        self.inner.numbered_ter = value;
    }
    #[getter]
    fn ter_ignores_type(&self) -> bool {
        self.inner.ter_ignores_type
    }
    #[setter]
    fn set_ter_ignores_type(&mut self, value: bool) {
        self.inner.ter_ignores_type = value;
    }
    #[getter]
    fn preserve_serial(&self) -> bool {
        self.inner.preserve_serial
    }
    #[setter]
    fn set_preserve_serial(&mut self, value: bool) {
        self.inner.preserve_serial = value;
    }
    #[getter]
    fn end_record(&self) -> bool {
        self.inner.end_record
    }
    #[setter]
    fn set_end_record(&mut self, value: bool) {
        self.inner.end_record = value;
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct BioStructure {
    inner: Arc<ck::BioStructure>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioStructure {
    fn to_molecule_with_params(
        &self,
        py: Python<'_>,
        params: &BioMoleculeParams,
    ) -> PyResult<crate::drawing_binding::Molecule> {
        self.inner
            .to_molecule_with_params(&params.inner)
            .map(crate::drawing_binding::Molecule::from_inner)
            .map_err(|e| molecule_error(py, e))
    }
    #[pyo3(signature=(*,sanitize=true,remove_hs=true,flavor=0,proximity_bonding=true))]
    fn to_molecule(
        &self,
        py: Python<'_>,
        sanitize: bool,
        remove_hs: bool,
        flavor: u32,
        proximity_bonding: bool,
    ) -> PyResult<crate::drawing_binding::Molecule> {
        self.inner
            .to_molecule_with_params(&ck::BioMoleculeParams {
                sanitize,
                remove_hs,
                flavor,
                proximity_bonding,
            })
            .map(crate::drawing_binding::Molecule::from_inner)
            .map_err(|e| molecule_error(py, e))
    }

    #[staticmethod]
    fn from_pdb(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::BioStructure::from_pdb(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| pdb_error(py, e))
    }
    #[staticmethod]
    fn from_pdb_with_params(
        py: Python<'_>,
        text: &str,
        params: &BioPdbReadParams,
    ) -> PyResult<Self> {
        ck::BioStructure::from_pdb_with_params(text, &params.inner)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| pdb_error(py, e))
    }
    #[staticmethod]
    fn from_mmcif(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::BioStructure::from_mmcif(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| mmcif_error(py, e))
    }
    #[staticmethod]
    fn from_text(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::BioStructure::from_text(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| read_error(py, e))
    }
    #[staticmethod]
    fn from_text_with_params(py: Python<'_>, text: &str, params: &BioReadParams) -> PyResult<Self> {
        ck::BioStructure::from_text_with_params(text, &params.inner)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| read_error(py, e))
    }
    #[staticmethod]
    fn read(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        ck::BioStructure::read(&expand_path(py, path)?)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| read_error(py, e))
    }
    #[staticmethod]
    fn read_with_format(py: Python<'_>, path: PathBuf, format: u8) -> PyResult<Self> {
        ck::BioStructure::read_with_format(&expand_path(py, path)?, format_from_code(format)?)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| read_error(py, e))
    }
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn input_format(&self) -> String {
        format!("{:?}", self.inner.input_format())
    }
    fn num_models(&self) -> usize {
        self.inner.num_models()
    }
    fn num_chains(&self) -> usize {
        self.inner.num_chains()
    }
    fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    fn num_entities(&self) -> usize {
        self.inner.num_entities()
    }
    fn models(&self) -> Vec<BioModelRow> {
        (0..self.inner.num_models())
            .map(|index| BioModelRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    fn chains(&self) -> Vec<BioChainRow> {
        (0..self.inner.num_chains())
            .map(|index| BioChainRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    fn residues(&self) -> Vec<BioResidueRow> {
        (0..self.inner.num_residues())
            .map(|index| BioResidueRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    fn atoms(&self) -> Vec<BioAtomRow> {
        (0..self.inner.num_atoms())
            .map(|index| BioAtomRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    fn entities(&self) -> Vec<BioEntityRow> {
        (0..self.inner.num_entities())
            .map(|index| BioEntityRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    fn protein(&self, py: Python<'_>) -> PyResult<Protein> {
        self.inner
            .protein()
            .map(|inner| Protein {
                inner: Arc::new(inner),
            })
            .map_err(|e| {
                annotate(
                    py,
                    ProteinReadError::new_err(e.to_string()),
                    "Projection",
                    &e,
                )
            })
    }
    fn to_mmcif_with_params(
        &self,
        py: Python<'_>,
        params: &BioMmcifWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_mmcif_with_params(&params.inner)
            .map_err(|e| write_error(py, e))
    }
    fn write_mmcif_with_params(
        &self,
        py: Python<'_>,
        path: PathBuf,
        params: &BioMmcifWriteParams,
    ) -> PyResult<()> {
        self.inner
            .write_mmcif_with_params(&expand_path(py, path)?, &params.inner)
            .map_err(|e| write_error(py, e))
    }
    fn to_pdb(&self, py: Python<'_>) -> PyResult<String> {
        self.inner.to_pdb().map_err(|e| pdb_write_error(py, e))
    }
    fn to_pdb_with_params(&self, py: Python<'_>, params: &BioPdbWriteParams) -> PyResult<String> {
        self.inner
            .to_pdb_with_params(&params.inner)
            .map_err(|e| pdb_write_error(py, e))
    }
    fn write_pdb(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        self.inner
            .write_pdb(&expand_path(py, path)?)
            .map_err(|e| pdb_write_error(py, e))
    }
    fn write_pdb_with_params(
        &self,
        py: Python<'_>,
        path: PathBuf,
        params: &BioPdbWriteParams,
    ) -> PyResult<()> {
        self.inner
            .write_pdb_with_params(&expand_path(py, path)?, &params.inner)
            .map_err(|e| pdb_write_error(py, e))
    }
    fn to_mmcif(&self, py: Python<'_>) -> PyResult<String> {
        self.inner.to_mmcif().map_err(|e| write_error(py, e))
    }
    fn write_mmcif(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        self.inner
            .write_mmcif(&expand_path(py, path)?)
            .map_err(|e| write_error(py, e))
    }
    fn with_selection(&self, py: Python<'_>, selection: &BioSelection) -> PyResult<Self> {
        self.inner
            .with_selection(&selection.inner)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| operation_error(py, e))
    }
    fn retain_selection_(&mut self, py: Python<'_>, selection: &BioSelection) -> PyResult<()> {
        Arc::make_mut(&mut self.inner)
            .retain_selection_(&selection.inner)
            .map_err(|e| operation_error(py, e))
    }
    fn with_translated_coordinates(&self, py: Python<'_>, offset: [f64; 3]) -> PyResult<Self> {
        self.inner
            .with_translated_coordinates(offset)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| operation_error(py, e))
    }
    fn translate_(&mut self, py: Python<'_>, offset: [f64; 3]) -> PyResult<()> {
        Arc::make_mut(&mut self.inner)
            .translate_(offset)
            .map_err(|e| operation_error(py, e))
    }
    fn selected_atom_ids(&self, py: Python<'_>, selection: &BioSelection) -> PyResult<Vec<usize>> {
        self.inner
            .selected_atom_ids(&selection.inner)
            .map(|x| x.into_iter().map(|id| id.index()).collect())
            .map_err(|e| {
                annotate(
                    py,
                    BioOperationError::new_err(e.to_string()),
                    "SelectionMatch",
                    &e,
                )
            })
    }
    fn __getitem__(&self, value: isize) -> PyResult<BioModelRow> {
        Ok(BioModelRow {
            inner: self.inner.clone(),
            index: index(value, self.inner.num_models(), "BioStructure model")?,
        })
    }
    fn __len__(&self) -> usize {
        self.inner.num_models()
    }
    fn __repr__(&self) -> String {
        format!(
            "BioStructure(models={}, chains={}, residues={}, atoms={}, entities={})",
            self.inner.num_models(),
            self.inner.num_chains(),
            self.inner.num_residues(),
            self.inner.num_atoms(),
            self.inner.num_entities()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioModelRow {
    inner: Arc<ck::BioStructure>,
    index: usize,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioModelRow {
    fn id(&self) -> usize {
        self.index
    }
    fn source_model_number(&self) -> Option<i32> {
        self.inner.models()[self.index].source_model_number()
    }
    fn chains(&self) -> Vec<BioChainRow> {
        let span = self.inner.models()[self.index].chain_span();
        (span.start() as usize..span.end() as usize)
            .map(|index| BioChainRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    fn __len__(&self) -> usize {
        self.inner.models()[self.index].chain_span().len() as usize
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioChainRow {
    inner: Arc<ck::BioStructure>,
    index: usize,
}
impl BioChainRow {
    fn row(&self) -> &ck::BioChainRow {
        &self.inner.chains()[self.index]
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioChainRow {
    fn id(&self) -> usize {
        self.index
    }
    fn model_id(&self) -> usize {
        self.row().model_id().index()
    }
    fn entity_id(&self) -> Option<usize> {
        self.row().entity_id().map(|id| id.index())
    }
    fn kind(&self) -> String {
        format!("{:?}", self.row().kind())
    }
    fn source(&self) -> ChainSourceIds {
        ChainSourceIds {
            inner: self.row().source().clone(),
        }
    }
    fn residues(&self) -> Vec<BioResidueRow> {
        let span = self.row().residue_span();
        (span.start() as usize..span.end() as usize)
            .map(|index| BioResidueRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    fn atoms(&self) -> Vec<BioAtomRow> {
        let span = self.row().residue_span();
        self.inner.residues()[span.start() as usize..span.end() as usize]
            .iter()
            .flat_map(|r| r.atom_span().start() as usize..r.atom_span().end() as usize)
            .map(|index| BioAtomRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    fn __len__(&self) -> usize {
        self.row().residue_span().len() as usize
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioResidueRow {
    inner: Arc<ck::BioStructure>,
    index: usize,
}
impl BioResidueRow {
    fn row(&self) -> &ck::BioResidueRow {
        &self.inner.residues()[self.index]
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioResidueRow {
    fn id(&self) -> usize {
        self.index
    }
    fn chain_id(&self) -> usize {
        self.row().chain_id().index()
    }
    fn name(&self) -> String {
        self.row().name().as_str().to_owned()
    }
    fn kind(&self) -> String {
        format!("{:?}", self.row().kind())
    }
    fn entity_kind(&self) -> String {
        format!("{:?}", self.row().entity_kind())
    }
    fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: ck::find_residue_info(self.row().name().as_str()),
        }
    }
    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(
            py,
            "ResidueCode",
            i64::from(ck::residue_code(self.row().name().as_str()).as_u16()),
        )
    }
    fn source(&self) -> ResidueSourceIds {
        ResidueSourceIds {
            inner: self.row().source().clone(),
        }
    }
    fn atoms(&self) -> Vec<BioAtomRow> {
        let span = self.row().atom_span();
        (span.start() as usize..span.end() as usize)
            .map(|index| BioAtomRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    fn __len__(&self) -> usize {
        self.row().atom_span().len() as usize
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioAtomRow {
    inner: Arc<ck::BioStructure>,
    index: usize,
}
impl BioAtomRow {
    fn row(&self) -> &ck::BioAtomRow {
        &self.inner.atoms()[self.index]
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioAtomRow {
    fn id(&self) -> usize {
        self.index
    }
    fn residue_id(&self) -> usize {
        self.row().residue_id().index()
    }
    fn name(&self) -> String {
        self.row().name().as_str().trim().to_owned()
    }
    fn element(&self) -> Element {
        Element::from_inner(self.row().element())
    }
    fn element_symbol(&self) -> &'static str {
        self.row().element().symbol()
    }
    fn position(&self) -> Option<(f64, f64, f64)> {
        self.inner
            .atom_position(ck::BioAtomId::new(self.index as u32))
            .map(|[x, y, z]| (x, y, z))
    }
    fn altloc(&self) -> Option<String> {
        self.row()
            .altloc()
            .map(|x| char::from(x.value()).to_string())
    }
    fn occupancy(&self) -> f64 {
        self.row().occupancy()
    }
    fn b_iso(&self) -> f64 {
        self.row().b_iso()
    }
    fn formal_charge(&self) -> i8 {
        self.row().formal_charge()
    }
    fn source(&self) -> AtomSourceIds {
        AtomSourceIds {
            inner: *self.row().source(),
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioEntityRow {
    inner: Arc<ck::BioStructure>,
    index: usize,
}
impl BioEntityRow {
    fn row(&self) -> &ck::BioEntityRow {
        &self.inner.entities()[self.index]
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioEntityRow {
    fn id(&self) -> usize {
        self.index
    }
    fn source(&self) -> EntitySourceIds {
        EntitySourceIds {
            inner: self.row().source().clone(),
        }
    }
    fn kind(&self) -> String {
        format!("{:?}", self.row().kind())
    }
    fn polymer_kind(&self) -> String {
        format!("{:?}", self.row().polymer_kind())
    }
    fn full_sequence(&self) -> Vec<String> {
        self.row().full_sequence().to_vec()
    }
    fn subchains(&self) -> Vec<String> {
        self.row().subchains().to_vec()
    }
    fn __len__(&self) -> usize {
        self.row().full_sequence().len()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ChainSourceIds {
    inner: ck::ChainSourceIds,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ChainSourceIds {
    fn auth_chain_id(&self) -> Option<String> {
        self.inner.auth_chain_id().map(|x| x.as_str().to_owned())
    }
    fn label_asym_id(&self) -> Option<&str> {
        self.inner.label_asym_id()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomSourceIds {
    inner: ck::AtomSourceIds,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl AtomSourceIds {
    fn serial(&self) -> Option<i32> {
        self.inner.serial().map(|x| x.value())
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct EntitySourceIds {
    inner: ck::EntitySourceIds,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl EntitySourceIds {
    fn source_entity_id(&self) -> &str {
        self.inner.source_entity_id()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct PdbSeqId {
    inner: ck::PdbSeqId,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PdbSeqId {
    fn seq_num(&self) -> i32 {
        self.inner.seq_num()
    }
    fn ins_code(&self) -> Option<String> {
        self.inner.ins_code().map(|x| char::from(x).to_string())
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ResidueSourceIds {
    inner: ck::ResidueSourceIds,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ResidueSourceIds {
    fn seq_id(&self) -> Option<PdbSeqId> {
        self.inner.seq_id().map(|inner| PdbSeqId { inner })
    }
    fn label_seq_id(&self) -> Option<i32> {
        self.inner.label_seq_id()
    }
    fn subchain_id(&self) -> Option<&str> {
        self.inner.subchain_id()
    }
    fn label_entity_id(&self) -> Option<&str> {
        self.inner.label_entity_id()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct Protein {
    inner: Arc<ck::Protein>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl Protein {
    fn to_molecule_with_params(
        &self,
        py: Python<'_>,
        params: &BioMoleculeParams,
    ) -> PyResult<crate::drawing_binding::Molecule> {
        self.inner
            .to_molecule_with_params(&params.inner)
            .map(crate::drawing_binding::Molecule::from_inner)
            .map_err(|e| molecule_error(py, e))
    }
    #[pyo3(signature=(*,sanitize=true,remove_hs=true,flavor=0,proximity_bonding=true))]
    fn to_molecule(
        &self,
        py: Python<'_>,
        sanitize: bool,
        remove_hs: bool,
        flavor: u32,
        proximity_bonding: bool,
    ) -> PyResult<crate::drawing_binding::Molecule> {
        self.inner
            .to_molecule_with_params(&ck::BioMoleculeParams {
                sanitize,
                remove_hs,
                flavor,
                proximity_bonding,
            })
            .map(crate::drawing_binding::Molecule::from_inner)
            .map_err(|e| molecule_error(py, e))
    }

    #[staticmethod]
    fn from_pdb(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Protein::from_pdb(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    #[staticmethod]
    fn from_pdb_with_params(
        py: Python<'_>,
        text: &str,
        params: &BioPdbReadParams,
    ) -> PyResult<Self> {
        ck::Protein::from_pdb_with_params(text, &params.inner)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    #[staticmethod]
    fn from_mmcif(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Protein::from_mmcif(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    #[staticmethod]
    fn from_text(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Protein::from_text(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    #[staticmethod]
    fn from_text_with_params(py: Python<'_>, text: &str, params: &BioReadParams) -> PyResult<Self> {
        ck::Protein::from_text_with_params(text, &params.inner)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    #[staticmethod]
    fn read(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        ck::Protein::read(&expand_path(py, path)?)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    #[staticmethod]
    fn read_with_format(py: Python<'_>, path: PathBuf, format: u8) -> PyResult<Self> {
        ck::Protein::read_with_format(&expand_path(py, path)?, format_from_code(format)?)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    fn input_format(&self) -> String {
        format!("{:?}", self.inner.input_format())
    }
    fn num_models(&self) -> usize {
        self.inner.num_models()
    }
    fn num_chains(&self) -> usize {
        self.inner.num_chains()
    }
    fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    fn chains(&self) -> Vec<ProteinChainRef> {
        self.inner
            .chains()
            .map(|x| ProteinChainRef {
                inner: self.inner.clone(),
                index: x.id().index(),
            })
            .collect()
    }
    fn residues(&self) -> Vec<ProteinResidueRef> {
        self.inner
            .residues()
            .map(|x| ProteinResidueRef {
                inner: self.inner.clone(),
                index: x.id().index(),
            })
            .collect()
    }
    fn atoms(&self) -> Vec<ProteinAtomRef> {
        self.inner
            .atoms()
            .map(|x| ProteinAtomRef {
                inner: self.inner.clone(),
                index: x.id().index(),
            })
            .collect()
    }
    fn with_selection(&self, py: Python<'_>, selection: &BioSelection) -> PyResult<Self> {
        self.inner
            .with_selection(&selection.inner)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| operation_error(py, e))
    }
    fn retain_selection_(&mut self, py: Python<'_>, selection: &BioSelection) -> PyResult<()> {
        Arc::make_mut(&mut self.inner)
            .retain_selection_(&selection.inner)
            .map_err(|e| operation_error(py, e))
    }
    fn with_translated_coordinates(&self, py: Python<'_>, offset: [f64; 3]) -> PyResult<Self> {
        self.inner
            .with_translated_coordinates(offset)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| operation_error(py, e))
    }
    fn translate_(&mut self, py: Python<'_>, offset: [f64; 3]) -> PyResult<()> {
        Arc::make_mut(&mut self.inner)
            .translate_(offset)
            .map_err(|e| operation_error(py, e))
    }
    fn selected_atom_ids(&self, py: Python<'_>, selection: &BioSelection) -> PyResult<Vec<usize>> {
        self.inner
            .selected_atom_ids(&selection.inner)
            .map(|x| x.into_iter().map(|id| id.index()).collect())
            .map_err(|e| {
                annotate(
                    py,
                    BioOperationError::new_err(e.to_string()),
                    "SelectionMatch",
                    &e,
                )
            })
    }
    fn __getitem__(&self, value: isize) -> PyResult<ProteinChainRef> {
        Ok(ProteinChainRef {
            inner: self.inner.clone(),
            index: index(value, self.inner.num_chains(), "Protein chain")?,
        })
    }
    fn __len__(&self) -> usize {
        self.inner.num_chains()
    }
    fn __repr__(&self) -> String {
        format!(
            "Protein(chains={}, residues={}, atoms={})",
            self.inner.num_chains(),
            self.inner.num_residues(),
            self.inner.num_atoms()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ProteinChainRef {
    inner: Arc<ck::Protein>,
    index: usize,
}
impl ProteinChainRef {
    fn view(&self) -> ck::ProteinChainRef<'_> {
        self.inner
            .chain(self.index)
            .expect("validated protein chain view")
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ProteinChainRef {
    fn id(&self) -> usize {
        self.index
    }
    fn kind(&self) -> String {
        format!("{:?}", self.view().kind())
    }
    fn residues(&self) -> Vec<ProteinResidueRef> {
        self.view()
            .residues()
            .map(|x| ProteinResidueRef {
                inner: self.inner.clone(),
                index: x.id().index(),
            })
            .collect()
    }
    fn atoms(&self) -> Vec<ProteinAtomRef> {
        self.view()
            .atoms()
            .map(|x| ProteinAtomRef {
                inner: self.inner.clone(),
                index: x.id().index(),
            })
            .collect()
    }
    fn __len__(&self) -> usize {
        self.view().residues().len()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ProteinResidueRef {
    inner: Arc<ck::Protein>,
    index: usize,
}
impl ProteinResidueRef {
    fn view(&self) -> ck::ProteinResidueRef<'_> {
        self.inner
            .residues()
            .nth(self.index)
            .expect("validated protein residue view")
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ProteinResidueRef {
    fn id(&self) -> usize {
        self.index
    }
    fn name(&self) -> String {
        self.view().name().as_str().to_owned()
    }
    fn kind(&self) -> String {
        format!("{:?}", self.view().kind())
    }
    fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: self.view().info(),
        }
    }
    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "ResidueCode", i64::from(self.view().code().as_u16()))
    }
    fn one_letter_code(&self) -> String {
        self.view().one_letter_code().to_string()
    }
    fn fasta_code(&self) -> String {
        self.view().fasta_code().to_string()
    }
    fn canonical_one_letter_code(&self) -> Option<String> {
        self.view()
            .info()
            .canonical_one_letter_code()
            .map(|x| x.to_string())
    }
    #[gen_stub(override_return_type(type_repr = "typing.Optional[ResidueCode]"))]
    fn parent_standard_code<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.view()
            .info()
            .parent_standard_code()
            .map(|x| enum_member(py, "ResidueCode", i64::from(x.as_u16())))
            .transpose()
    }
    fn is_modified_amino_acid(&self) -> bool {
        self.view().info().is_modified_amino_acid()
    }
    fn is_standard(&self) -> bool {
        self.view().is_standard()
    }
    fn atoms(&self) -> Vec<ProteinAtomRef> {
        self.view()
            .atoms()
            .map(|x| ProteinAtomRef {
                inner: self.inner.clone(),
                index: x.id().index(),
            })
            .collect()
    }
    fn __len__(&self) -> usize {
        self.view().atoms().len()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ProteinAtomRef {
    inner: Arc<ck::Protein>,
    index: usize,
}
impl ProteinAtomRef {
    fn view(&self) -> ck::ProteinAtomRef<'_> {
        self.inner
            .atoms()
            .nth(self.index)
            .expect("validated protein atom view")
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ProteinAtomRef {
    fn id(&self) -> usize {
        self.index
    }
    fn name(&self) -> String {
        self.view().name().as_str().trim().to_owned()
    }
    fn element(&self) -> Element {
        Element::from_inner(self.view().element())
    }
    fn element_symbol(&self) -> &'static str {
        self.view().element().symbol()
    }
    fn atomic_number(&self) -> u8 {
        self.view().element().atomic_number()
    }
    fn position(&self) -> (f64, f64, f64) {
        let [x, y, z] = self.view().position();
        (x, y, z)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioMoleculeParams {
    inner: ck::BioMoleculeParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioMoleculeParams {
    #[new]
    #[pyo3(signature=(*,sanitize=true,remove_hs=true,flavor=0,proximity_bonding=true))]
    fn new(sanitize: bool, remove_hs: bool, flavor: u32, proximity_bonding: bool) -> Self {
        Self {
            inner: ck::BioMoleculeParams {
                sanitize,
                remove_hs,
                flavor,
                proximity_bonding,
            },
        }
    }
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[getter]
    fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    #[getter]
    fn flavor(&self) -> u32 {
        self.inner.flavor
    }
    #[getter]
    fn proximity_bonding(&self) -> bool {
        self.inner.proximity_bonding
    }
}
fn molecule_error(py: Python<'_>, source: ck::BioMoleculeError) -> PyErr {
    let kind = match &source {
        ck::BioMoleculeError::Conversion(_) => "Conversion",
        ck::BioMoleculeError::Construction(_) => "Construction",
    };
    let error = annotate(
        py,
        BioMoleculeError::new_err(source.to_string()),
        kind,
        &source,
    );
    if let ck::BioMoleculeError::Conversion(source) = source {
        let kind = match &source {
            ck::BioMoleculeConversionError::Structure(_) => "Structure",
            ck::BioMoleculeConversionError::Topology(_) => "Topology",
            ck::BioMoleculeConversionError::Chemistry(_) => "Chemistry",
            ck::BioMoleculeConversionError::Sanitize(_) => "Sanitize",
            ck::BioMoleculeConversionError::Hydrogens(_) => "Hydrogens",
            ck::BioMoleculeConversionError::Valence(_) => "Valence",
            ck::BioMoleculeConversionError::Stereo(_) => "Stereo",
        };
        error.set_cause(
            py,
            Some(annotate(
                py,
                BioMoleculeConversionError::new_err(source.to_string()),
                kind,
                &source,
            )),
        );
    }
    error
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    let members = PyDict::new(py);
    for (name, code) in [
        ("Unknown", 0),
        ("Detect", 1),
        ("Pdb", 2),
        ("Mmcif", 3),
        ("Mmjson", 4),
        ("ChemComp", 5),
    ] {
        members.set_item(name, code)?;
    }
    let formats = py
        .import("enum")?
        .getattr("IntEnum")?
        .call1(("BioCoordinateFormat", members))?;
    formats.setattr("__module__", "cosmolkit")?;
    module.add("BioCoordinateFormat", formats)?;
    module.add_class::<BioMoleculeParams>()?;
    module.add("BioMoleculeError", py.get_type::<BioMoleculeError>())?;
    module.add(
        "BioMoleculeConversionError",
        py.get_type::<BioMoleculeConversionError>(),
    )?;
    module.add("BioReadError", py.get_type::<BioReadError>())?;
    module.add("BioPdbReadError", py.get_type::<BioPdbReadError>())?;
    module.add("BioMmcifReadError", py.get_type::<BioMmcifReadError>())?;
    module.add("ProteinReadError", py.get_type::<ProteinReadError>())?;
    module.add("BioOperationError", py.get_type::<BioOperationError>())?;
    module.add("BioPdbWriteError", py.get_type::<BioPdbWriteError>())?;
    module.add_class::<BioMmcifWriteParams>()?;
    module.add_class::<BioPdbWriteParams>()?;
    module.add("BioMmcifWriteError", py.get_type::<BioMmcifWriteError>())?;
    module.add(
        "BioSelectionParseError",
        py.get_type::<BioSelectionParseError>(),
    )?;
    module.add_class::<BioReadParams>()?;
    module.add_class::<BioPdbReadParams>()?;
    module.add_class::<BioSelection>()?;
    module.add_class::<BioStructure>()?;
    module.add_class::<Protein>()?;
    module.add_class::<BioModelRow>()?;
    module.add_class::<BioChainRow>()?;
    module.add_class::<BioResidueRow>()?;
    module.add_class::<BioAtomRow>()?;
    module.add_class::<BioEntityRow>()?;
    module.add_class::<ProteinChainRef>()?;
    module.add_class::<ProteinResidueRef>()?;
    module.add_class::<ProteinAtomRef>()?;
    module.add_class::<AtomSourceIds>()?;
    module.add_class::<ChainSourceIds>()?;
    module.add_class::<ResidueSourceIds>()?;
    module.add_class::<EntitySourceIds>()?;
    module.add_class::<PdbSeqId>()?;
    Ok(())
}

fn pdb_write_error(py: Python<'_>, source: ck::BioPdbWriteError) -> PyErr {
    let kind = match &source {
        ck::BioPdbWriteError::InvalidStructure(_) => "InvalidStructure",
        ck::BioPdbWriteError::ChainNameTooLong { .. } => "ChainNameTooLong",
        ck::BioPdbWriteError::NegativeSerial { .. } => "NegativeSerial",
        ck::BioPdbWriteError::SerialIncrementOverflow { .. } => "SerialIncrementOverflow",
        ck::BioPdbWriteError::SerialOffsetOverflow { .. } => "SerialOffsetOverflow",
        ck::BioPdbWriteError::NegativeBase36Value { .. } => "NegativeBase36Value",
        ck::BioPdbWriteError::UnrepresentableField { .. } => "UnrepresentableField",
        ck::BioPdbWriteError::Io { .. } => "Io",
    };
    let error = if let ck::BioPdbWriteError::Io { path, source: io } = &source {
        PyOSError::new_err((
            io.raw_os_error(),
            source.to_string(),
            path.to_string_lossy().into_owned(),
        ))
    } else {
        BioPdbWriteError::new_err(source.to_string())
    };
    annotate(py, error, kind, &source)
}
