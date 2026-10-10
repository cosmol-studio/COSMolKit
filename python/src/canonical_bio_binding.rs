//! Shared hierarchy snapshots and thin Python projections of canonical BIO APIs.

use crate::canonical_bio_residue::{ResidueInfo, enum_member};
use crate::canonical_element_metadata::Element;
use ::cosmolkit as ck;
use pyo3::exceptions::{PyIndexError, PyOSError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::path::PathBuf;
use std::sync::Arc;

fn format_input(value: &Bound<'_, PyAny>) -> PyResult<u8> {
    let code = crate::canonical_atom_bond::enum_code(value, "BioCoordinateFormat")?;
    u8::try_from(code)
        .map_err(|_| pyo3::exceptions::PyOverflowError::new_err("format code is out of range"))
}

fn kind_member<'py>(
    py: Python<'py>,
    name: &str,
    value: impl std::fmt::Debug,
) -> PyResult<Bound<'py, PyAny>> {
    py.import("cosmolkit")?
        .getattr(name)?
        .call1((format!("{value:?}"),))
}

pyo3::create_exception!(
    cosmolkit,
    BioReadError,
    PyValueError,
    "Biological structure input could not be read or parsed."
);
pyo3::create_exception!(
    cosmolkit,
    BioMoleculeError,
    PyValueError,
    "A biological structure could not be projected into a molecular graph."
);
pyo3::create_exception!(
    cosmolkit,
    BioMoleculeConversionError,
    PyValueError,
    "Molecule-to-biological-structure conversion failed."
);
pyo3::create_exception!(
    cosmolkit,
    BioPdbReadError,
    BioReadError,
    "PDB biological structure input could not be read or parsed."
);
pyo3::create_exception!(
    cosmolkit,
    BioMmcifReadError,
    BioReadError,
    "mmCIF biological structure input could not be read or parsed."
);
pyo3::create_exception!(
    cosmolkit,
    ProteinReadError,
    BioReadError,
    "Protein input could not be read into the structural model."
);
pyo3::create_exception!(
    cosmolkit,
    BioOperationError,
    PyValueError,
    "A biological structure operation failed without committing partial changes."
);
pyo3::create_exception!(
    cosmolkit,
    BioMmcifWriteError,
    PyValueError,
    "The biological structure could not be serialized as mmCIF."
);
pyo3::create_exception!(
    cosmolkit,
    BioSelectionParseError,
    PyValueError,
    "The structural selection expression is invalid."
);
pyo3::create_exception!(
    cosmolkit,
    BioPdbWriteError,
    PyValueError,
    "The biological structure could not be serialized as PDB."
);
pyo3::create_exception!(
    cosmolkit,
    BioStructureError,
    PyValueError,
    "Biological hierarchy, coordinates or row references are structurally inconsistent."
);

pub(crate) fn structure_error(py: Python<'_>, source: &ck::BioStructureError) -> PyErr {
    let error = BioStructureError::new_err(source.to_string());
    let fields = || -> PyResult<&str> {
        let value = error.value(py);
        Ok(match source {
            ck::BioStructureError::RowIndexTooLarge { value: index } => {
                value.setattr("value", *index)?;
                "RowIndexTooLarge"
            }
            ck::BioStructureError::RowSpanOverflow { start, len } => {
                value.setattr("start", *start)?;
                value.setattr("len", *len)?;
                "RowSpanOverflow"
            }
            ck::BioStructureError::RowSpanOutOfBounds {
                start,
                len,
                table_len,
            } => {
                value.setattr("start", *start)?;
                value.setattr("len", *len)?;
                value.setattr("table_len", *table_len)?;
                "RowSpanOutOfBounds"
            }
            ck::BioStructureError::TableTooLarge { table, len } => {
                value.setattr("table", *table)?;
                value.setattr("len", *len)?;
                "TableTooLarge"
            }
            ck::BioStructureError::NonContiguousSpan {
                table,
                expected_start,
                actual_start,
            } => {
                value.setattr("table", *table)?;
                value.setattr("expected_start", *expected_start)?;
                value.setattr("actual_start", *actual_start)?;
                "NonContiguousSpan"
            }
            ck::BioStructureError::IncompleteCoverage {
                table,
                covered,
                table_len,
            } => {
                value.setattr("table", *table)?;
                value.setattr("covered", *covered)?;
                value.setattr("table_len", *table_len)?;
                "IncompleteCoverage"
            }
            ck::BioStructureError::ParentMismatch { table, index } => {
                value.setattr("table", *table)?;
                value.setattr("index", *index)?;
                "ParentMismatch"
            }
            ck::BioStructureError::RowReferenceOutOfBounds {
                table,
                index,
                table_len,
            } => {
                value.setattr("table", *table)?;
                value.setattr("index", *index)?;
                value.setattr("table_len", *table_len)?;
                "RowReferenceOutOfBounds"
            }
            ck::BioStructureError::CoordinateCountMismatch {
                atom_count,
                coordinate_count,
            } => {
                value.setattr("atom_count", *atom_count)?;
                value.setattr("coordinate_count", *coordinate_count)?;
                "CoordinateCountMismatch"
            }
            ck::BioStructureError::EntitySubchainMismatch {
                entity_id,
                subchain,
            } => {
                value.setattr("entity_id", entity_id.index())?;
                value.setattr("subchain", subchain)?;
                "EntitySubchainMismatch"
            }
            ck::BioStructureError::EmptyResidueSpan { operation } => {
                value.setattr("operation", *operation)?;
                "EmptyResidueSpan"
            }
            ck::BioStructureError::ImpossibleCrystalAngle => "ImpossibleCrystalAngle",
            ck::BioStructureError::AtomNotFound => "AtomNotFound",
        })
    };
    match fields() {
        Ok(kind) => annotate(py, error, kind, source),
        Err(error) => error,
    }
}

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
            .setattr("_line_number", source.line_number())?;
        error.value(py).setattr(
            "_stage",
            Py::new(
                py,
                crate::canonical_error_values::BioPdbReadStage::from(source.stage()),
            )?,
        )?;
        error.value(py).setattr(
            "_record_tag",
            source.record_tag().map(|x| PyBytes::new(py, &x)),
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
    let fields = || -> PyResult<()> {
        error.value(py).setattr(
            "_stage",
            Py::new(
                py,
                crate::canonical_error_values::BioMmcifReadStage::from(source.stage()),
            )?,
        )
    };
    match fields() {
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

/// Writable configuration for PDB/mmCIF structural text or file reading.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct BioReadParams {
    inner: ck::BioReadParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioReadParams {
    /// Configure PDB/mmCIF structural text or file reading; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(*, format=0, source_name="<string>"))]
    fn new(
        #[pyo3(from_py_with = format_input)]
        #[gen_stub(override_type(
            type_repr = "BioCoordinateFormat | builtins.str | builtins.int"
        ))]
        format: u8,
        source_name: &str,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BioReadParams {
                format: format_from_code(format)?,
                source_name: source_name.to_owned(),
            },
        })
    }
    /// Input/output format selector accepted by the corresponding operation.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "BioCoordinateFormat"))]
    fn format<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "BioCoordinateFormat", self.inner.format as i64)
    }
    /// Source name retained as input provenance.
    #[getter]
    fn source_name(&self) -> &str {
        &self.inner.source_name
    }
}

/// Writable configuration for PDB record handling.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct BioPdbReadParams {
    inner: ck::BioPdbReadParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioPdbReadParams {
    /// Configure PDB record handling; omitted fields use the defaults shown in the signature.
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
    /// Maximum accepted PDB input line length.
    #[getter]
    fn max_line_length(&self) -> i32 {
        self.inner.max_line_length
    }
    /// Whether non-ASCII bytes in PDB input are checked.
    #[getter]
    fn check_non_ascii(&self) -> bool {
        self.inner.check_non_ascii
    }
    /// Whether PDB TER records are ignored.
    #[getter]
    fn ignore_ter(&self) -> bool {
        self.inner.ignore_ter
    }
    /// Whether a TER record starts a separate chain.
    #[getter]
    fn split_chain_on_ter(&self) -> bool {
        self.inner.split_chain_on_ter
    }
    /// Whether PDB REMARK records are omitted from retained source metadata.
    #[getter]
    fn skip_remarks(&self) -> bool {
        self.inner.skip_remarks
    }
}

/// Structural selection expressed using CID syntax; applies to hierarchy rows rather than chemical SMARTS predicates.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioSelection {
    inner: ck::BioSelection,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioSelection {
    /// Parse a structural CID selection expression; invalid syntax raises BioSelectionError.
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
    /// Return the canonical CID text for this structural selection.
    fn to_cid(&self) -> String {
        self.inner.to_cid()
    }
}

/// Writable configuration for mmCIF output category and formatting selection.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct BioMmcifWriteParams {
    inner: ck::BioMmcifWriteParams,
}
#[cosmolkit_macros::python_configuration(existing_setters)]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioMmcifWriteParams {
    /// Configure mmCIF output category and formatting selection; omitted fields use the defaults shown in the signature.
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
    /// Read/write the collective group switch without changing formatting.
    /// A mixed selection reads false; assignment sets all 32 group flags.
    #[getter]
    fn all_groups(&self) -> bool {
        [
            self.inner.atoms,
            self.inner.block_name,
            self.inner.entry,
            self.inner.database_status,
            self.inner.author,
            self.inner.cell,
            self.inner.symmetry,
            self.inner.entity,
            self.inner.entity_poly,
            self.inner.struct_ref,
            self.inner.chem_comp,
            self.inner.exptl,
            self.inner.diffrn,
            self.inner.reflns,
            self.inner.refine,
            self.inner.title_keywords,
            self.inner.ncs,
            self.inner.struct_asym,
            self.inner.origx,
            self.inner.struct_conf,
            self.inner.struct_sheet,
            self.inner.struct_biol,
            self.inner.assembly,
            self.inner.conn,
            self.inner.cis,
            self.inner.modres,
            self.inner.scale,
            self.inner.atom_type,
            self.inner.entity_poly_seq,
            self.inner.tls,
            self.inner.software,
            self.inner.group_pdb,
        ]
        .into_iter()
        .all(|enabled| enabled)
    }
    #[setter]
    fn set_all_groups(&mut self, value: bool) {
        self.inner.atoms = value;
        self.inner.block_name = value;
        self.inner.entry = value;
        self.inner.database_status = value;
        self.inner.author = value;
        self.inner.cell = value;
        self.inner.symmetry = value;
        self.inner.entity = value;
        self.inner.entity_poly = value;
        self.inner.struct_ref = value;
        self.inner.chem_comp = value;
        self.inner.exptl = value;
        self.inner.diffrn = value;
        self.inner.reflns = value;
        self.inner.refine = value;
        self.inner.title_keywords = value;
        self.inner.ncs = value;
        self.inner.struct_asym = value;
        self.inner.origx = value;
        self.inner.struct_conf = value;
        self.inner.struct_sheet = value;
        self.inner.struct_biol = value;
        self.inner.assembly = value;
        self.inner.conn = value;
        self.inner.cis = value;
        self.inner.modres = value;
        self.inner.scale = value;
        self.inner.atom_type = value;
        self.inner.entity_poly_seq = value;
        self.inner.tls = value;
        self.inner.software = value;
        self.inner.group_pdb = value;
    }
    /// Whether the atom_site coordinate rows are included in mmCIF output.
    #[getter]
    fn atoms(&self) -> bool {
        self.inner.atoms
    }
    #[setter]
    fn set_atoms(&mut self, value: Option<bool>) {
        self.inner.atoms = value.unwrap_or(true);
    }
    /// Name used for the mmCIF data block.
    #[getter]
    fn block_name(&self) -> bool {
        self.inner.block_name
    }
    #[setter]
    fn set_block_name(&mut self, value: Option<bool>) {
        self.inner.block_name = value.unwrap_or(true);
    }
    /// Whether mmCIF entry identification is written.
    #[getter]
    fn entry(&self) -> bool {
        self.inner.entry
    }
    #[setter]
    fn set_entry(&mut self, value: Option<bool>) {
        self.inner.entry = value.unwrap_or(true);
    }
    /// Whether mmCIF database-status metadata is written.
    #[getter]
    fn database_status(&self) -> bool {
        self.inner.database_status
    }
    #[setter]
    fn set_database_status(&mut self, value: Option<bool>) {
        self.inner.database_status = value.unwrap_or(true);
    }
    /// Whether mmCIF audit-author metadata is written.
    #[getter]
    fn author(&self) -> bool {
        self.inner.author
    }
    #[setter]
    fn set_author(&mut self, value: Option<bool>) {
        self.inner.author = value.unwrap_or(true);
    }
    /// Whether the crystallographic cell category is included in mmCIF output.
    #[getter]
    fn cell(&self) -> bool {
        self.inner.cell
    }
    #[setter]
    fn set_cell(&mut self, value: Option<bool>) {
        self.inner.cell = value.unwrap_or(true);
    }
    /// Whether mmCIF crystallographic symmetry information is written.
    #[getter]
    fn symmetry(&self) -> bool {
        self.inner.symmetry
    }
    #[setter]
    fn set_symmetry(&mut self, value: Option<bool>) {
        self.inner.symmetry = value.unwrap_or(true);
    }
    /// Whether mmCIF entity definitions are written.
    #[getter]
    fn entity(&self) -> bool {
        self.inner.entity
    }
    #[setter]
    fn set_entity(&mut self, value: Option<bool>) {
        self.inner.entity = value.unwrap_or(true);
    }
    /// Whether mmCIF polymer entity descriptions are written.
    #[getter]
    fn entity_poly(&self) -> bool {
        self.inner.entity_poly
    }
    #[setter]
    fn set_entity_poly(&mut self, value: Option<bool>) {
        self.inner.entity_poly = value.unwrap_or(true);
    }
    /// Whether mmCIF external sequence/database references are written.
    #[getter]
    fn struct_ref(&self) -> bool {
        self.inner.struct_ref
    }
    #[setter]
    fn set_struct_ref(&mut self, value: Option<bool>) {
        self.inner.struct_ref = value.unwrap_or(true);
    }
    /// Whether mmCIF chemical component descriptions are written.
    #[getter]
    fn chem_comp(&self) -> bool {
        self.inner.chem_comp
    }
    #[setter]
    fn set_chem_comp(&mut self, value: Option<bool>) {
        self.inner.chem_comp = value.unwrap_or(true);
    }
    /// Whether mmCIF experimental-method metadata is written.
    #[getter]
    fn exptl(&self) -> bool {
        self.inner.exptl
    }
    #[setter]
    fn set_exptl(&mut self, value: Option<bool>) {
        self.inner.exptl = value.unwrap_or(true);
    }
    /// Whether mmCIF diffraction experiment metadata is written.
    #[getter]
    fn diffrn(&self) -> bool {
        self.inner.diffrn
    }
    #[setter]
    fn set_diffrn(&mut self, value: Option<bool>) {
        self.inner.diffrn = value.unwrap_or(true);
    }
    /// Whether mmCIF reflection statistics are written.
    #[getter]
    fn reflns(&self) -> bool {
        self.inner.reflns
    }
    #[setter]
    fn set_reflns(&mut self, value: Option<bool>) {
        self.inner.reflns = value.unwrap_or(true);
    }
    /// Whether mmCIF refinement statistics are written.
    #[getter]
    fn refine(&self) -> bool {
        self.inner.refine
    }
    #[setter]
    fn set_refine(&mut self, value: Option<bool>) {
        self.inner.refine = value.unwrap_or(true);
    }
    /// Whether mmCIF structure title and keywords are written.
    #[getter]
    fn title_keywords(&self) -> bool {
        self.inner.title_keywords
    }
    #[setter]
    fn set_title_keywords(&mut self, value: Option<bool>) {
        self.inner.title_keywords = value.unwrap_or(true);
    }
    /// Whether mmCIF noncrystallographic symmetry operators are written.
    #[getter]
    fn ncs(&self) -> bool {
        self.inner.ncs
    }
    #[setter]
    fn set_ncs(&mut self, value: Option<bool>) {
        self.inner.ncs = value.unwrap_or(true);
    }
    /// Whether mmCIF asymmetric-unit chain/entity associations are written.
    #[getter]
    fn struct_asym(&self) -> bool {
        self.inner.struct_asym
    }
    #[setter]
    fn set_struct_asym(&mut self, value: Option<bool>) {
        self.inner.struct_asym = value.unwrap_or(true);
    }
    /// Whether original-coordinate transformation metadata is written.
    #[getter]
    fn origx(&self) -> bool {
        self.inner.origx
    }
    #[setter]
    fn set_origx(&mut self, value: Option<bool>) {
        self.inner.origx = value.unwrap_or(true);
    }
    /// Whether mmCIF secondary-structure helix annotations are written.
    #[getter]
    fn struct_conf(&self) -> bool {
        self.inner.struct_conf
    }
    #[setter]
    fn set_struct_conf(&mut self, value: Option<bool>) {
        self.inner.struct_conf = value.unwrap_or(true);
    }
    /// Whether mmCIF secondary-structure sheet annotations are written.
    #[getter]
    fn struct_sheet(&self) -> bool {
        self.inner.struct_sheet
    }
    #[setter]
    fn set_struct_sheet(&mut self, value: Option<bool>) {
        self.inner.struct_sheet = value.unwrap_or(true);
    }
    /// Whether mmCIF biological-structure descriptions are written.
    #[getter]
    fn struct_biol(&self) -> bool {
        self.inner.struct_biol
    }
    #[setter]
    fn set_struct_biol(&mut self, value: Option<bool>) {
        self.inner.struct_biol = value.unwrap_or(true);
    }
    /// Whether mmCIF biological assembly definitions are written.
    #[getter]
    fn assembly(&self) -> bool {
        self.inner.assembly
    }
    #[setter]
    fn set_assembly(&mut self, value: Option<bool>) {
        self.inner.assembly = value.unwrap_or(true);
    }
    /// Whether mmCIF inter-atom connection annotations are written.
    #[getter]
    fn conn(&self) -> bool {
        self.inner.conn
    }
    #[setter]
    fn set_conn(&mut self, value: Option<bool>) {
        self.inner.conn = value.unwrap_or(true);
    }
    /// Whether mmCIF cis-peptide annotations are written.
    #[getter]
    fn cis(&self) -> bool {
        self.inner.cis
    }
    #[setter]
    fn set_cis(&mut self, value: Option<bool>) {
        self.inner.cis = value.unwrap_or(true);
    }
    /// Whether mmCIF modified-residue annotations are written.
    #[getter]
    fn modres(&self) -> bool {
        self.inner.modres
    }
    #[setter]
    fn set_modres(&mut self, value: Option<bool>) {
        self.inner.modres = value.unwrap_or(true);
    }
    /// Whether mmCIF fractional-coordinate transformations are written.
    #[getter]
    fn scale(&self) -> bool {
        self.inner.scale
    }
    #[setter]
    fn set_scale(&mut self, value: Option<bool>) {
        self.inner.scale = value.unwrap_or(true);
    }
    /// Whether the atom_type category is included in mmCIF output.
    #[getter]
    fn atom_type(&self) -> bool {
        self.inner.atom_type
    }
    #[setter]
    fn set_atom_type(&mut self, value: Option<bool>) {
        self.inner.atom_type = value.unwrap_or(true);
    }
    /// Whether mmCIF complete polymer entity sequences are written.
    #[getter]
    fn entity_poly_seq(&self) -> bool {
        self.inner.entity_poly_seq
    }
    #[setter]
    fn set_entity_poly_seq(&mut self, value: Option<bool>) {
        self.inner.entity_poly_seq = value.unwrap_or(true);
    }
    /// Whether mmCIF translation/libration/screw refinement metadata is written.
    #[getter]
    fn tls(&self) -> bool {
        self.inner.tls
    }
    #[setter]
    fn set_tls(&mut self, value: Option<bool>) {
        self.inner.tls = value.unwrap_or(true);
    }
    /// Whether structural software metadata is written.
    #[getter]
    fn software(&self) -> bool {
        self.inner.software
    }
    #[setter]
    fn set_software(&mut self, value: Option<bool>) {
        self.inner.software = value.unwrap_or(true);
    }
    /// Whether atom-site mmCIF output includes the PDB ATOM/HETATM group field.
    #[getter]
    fn group_pdb(&self) -> bool {
        self.inner.group_pdb
    }
    #[setter]
    fn set_group_pdb(&mut self, value: Option<bool>) {
        self.inner.group_pdb = value.unwrap_or(true);
    }
    /// Whether all author identifiers are included in mmCIF output.
    #[getter]
    fn auth_all(&self) -> bool {
        self.inner.auth_all
    }
    #[setter]
    fn set_auth_all(&mut self, value: Option<bool>) {
        self.inner.auth_all = value.unwrap_or(false);
    }
    /// Whether single-row mmCIF categories use key/value pairs instead of loops.
    #[getter]
    fn prefer_pairs(&self) -> bool {
        self.inner.prefer_pairs
    }
    #[setter]
    fn set_prefer_pairs(&mut self, value: bool) {
        self.inner.prefer_pairs = value;
    }
    /// Whether compact mmCIF formatting is enabled.
    #[getter]
    fn compact(&self) -> bool {
        self.inner.compact
    }
    #[setter]
    fn set_compact(&mut self, value: bool) {
        self.inner.compact = value;
    }
    /// Whether compact hash-line formatting is used for mmCIF output.
    #[getter]
    fn misuse_hash(&self) -> bool {
        self.inner.misuse_hash
    }
    #[setter]
    fn set_misuse_hash(&mut self, value: bool) {
        self.inner.misuse_hash = value;
    }
    /// Whether mmCIF key/value pairs are column-aligned.
    #[getter]
    fn align_pairs(&self) -> u16 {
        self.inner.align_pairs
    }
    #[setter]
    fn set_align_pairs(&mut self, value: u16) {
        self.inner.align_pairs = value;
    }
    /// Whether mmCIF loop columns are aligned.
    #[getter]
    fn align_loops(&self) -> u16 {
        self.inner.align_loops
    }
    #[setter]
    fn set_align_loops(&mut self, value: u16) {
        self.inner.align_loops = value;
    }
}

/// Writable configuration for PDB text formatting.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct BioPdbWriteParams {
    inner: ck::BioPdbWriteParams,
}
#[cosmolkit_macros::python_configuration(existing_setters)]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioPdbWriteParams {
    /// Configure PDB text formatting; omitted fields use the defaults shown in the signature.
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
    /// Whether PDB output includes TER records.
    #[getter]
    fn ter_records(&self) -> bool {
        self.inner.ter_records
    }
    #[setter]
    fn set_ter_records(&mut self, value: bool) {
        self.inner.ter_records = value;
    }
    /// Whether generated PDB TER records receive serial numbers.
    #[getter]
    fn numbered_ter(&self) -> bool {
        self.inner.numbered_ter
    }
    #[setter]
    fn set_numbered_ter(&mut self, value: bool) {
        self.inner.numbered_ter = value;
    }
    /// Whether TER output ignores the polymer/non-polymer residue classification.
    #[getter]
    fn ter_ignores_type(&self) -> bool {
        self.inner.ter_ignores_type
    }
    #[setter]
    fn set_ter_ignores_type(&mut self, value: bool) {
        self.inner.ter_ignores_type = value;
    }
    /// Whether atom serial numbers from the input are preserved in PDB output.
    #[getter]
    fn preserve_serial(&self) -> bool {
        self.inner.preserve_serial
    }
    #[setter]
    fn set_preserve_serial(&mut self, value: bool) {
        self.inner.preserve_serial = value;
    }
    /// Whether PDB output includes the final END record.
    #[getter]
    fn end_record(&self) -> bool {
        self.inner.end_record
    }
    #[setter]
    fn set_end_record(&mut self, value: bool) {
        self.inner.end_record = value;
    }
}

/// Fixed-width structural atom name with explicit byte and text access.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomName {
    inner: ck::AtomName,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomName {
    /// Construct an atom name from ASCII bytes; return None for invalid width or non-ASCII input.
    #[staticmethod]
    fn from_ascii(bytes: &Bound<'_, PyBytes>) -> Option<Self> {
        ck::AtomName::from_ascii(bytes.as_bytes()).map(|inner| Self { inner })
    }
    /// Return the stored source bytes without text normalization.
    fn as_bytes<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, self.inner.as_bytes())
    }
    /// Return the text representation of the stored name.
    fn as_str(&self) -> &str {
        self.inner.as_str()
    }
}

/// Stored alternate-location label for a structural atom site.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AltLocLabel {
    inner: ck::AltLocLabel,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AltLocLabel {
    /// Construct a AltLocLabel value from the supplied inputs.
    #[new]
    fn new(value: u8) -> Self {
        Self {
            inner: ck::AltLocLabel::new(value),
        }
    }
    /// Return the stored alternate-location character code.
    fn value(&self) -> u8 {
        self.inner.value()
    }
}

/// Alternate-location selection request: any alternate or an exact label.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AltLocRequest {
    inner: ck::AltLocRequest,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AltLocRequest {
    /// AltLocRequest value selecting any.
    #[classattr]
    #[pyo3(name = "Any")]
    fn any() -> AltLocRequest {
        Self {
            inner: ck::AltLocRequest::Any,
        }
    }
    /// Construct a request for the supplied alternate-location label.
    #[staticmethod]
    #[pyo3(name = "Exact", signature = (altloc))]
    fn exact(altloc: Option<&AltLocLabel>) -> Self {
        Self {
            inner: ck::AltLocRequest::Exact(altloc.map(|label| label.inner)),
        }
    }
}

/// Three-dimensional affine transformation used by structural symmetry and assembly metadata.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioTransform {
    pub(crate) inner: ck::BioTransform,
}

/// Detached atom-ordered Cartesian coordinates for a BIO hierarchy.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioCoordinateBlock {
    inner: Arc<ck::BioStructure>,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioCoordinateBlock {
    /// Return Cartesian (x, y, z) positions in structural atom-row order.
    #[getter]
    fn positions(&self) -> Vec<[f64; 3]> {
        self.inner.coordinates().positions().to_vec()
    }
    fn __len__(&self) -> usize {
        self.inner.coordinates().len()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioTransform {
    /// Return the 3 by 3 linear transformation matrix in row order.
    #[getter]
    fn matrix(&self) -> [[f64; 3]; 3] {
        *self.inner.matrix()
    }
    /// Return the three-component translation vector.
    #[getter]
    fn translation(&self) -> [f64; 3] {
        *self.inner.translation()
    }
    /// Source affine-transform approximate equality, including its asymmetric NaN treatment.
    fn approx(&self, other: &Self, epsilon: f64) -> bool {
        self.inner.approx(&other.inner, epsilon)
    }
}

/// Crystallographic unit-cell and space-group information.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioCrystalInfo {
    // Retain the immutable structure snapshot without copying symmetry images.
    inner: Arc<ck::BioStructure>,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioCrystalInfo {
    /// Crystallographic unit-cell dimensions and angles.
    fn cell(&self) -> crate::canonical_bio_values::BioCrystalCell {
        crate::canonical_bio_values::BioCrystalCell {
            inner: self.inner.crystal().expect("existing crystal view").cell(),
        }
    }
    /// Crystallographic space-group number, when available.
    fn space_group_number(&self) -> Option<i32> {
        self.inner
            .crystal()
            .expect("a crystal view is constructed only for an existing crystal")
            .space_group_number()
    }
}

/// Owned PDB/mmCIF structural hierarchy with models, chains, residues, atoms, coordinates and metadata.
///
/// from_pdb()/from_mmcif() take text; read() takes a path. Local row IDs are distinct
/// from retained source serials and author/label sequence identifiers.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct BioStructure {
    inner: Arc<ck::BioStructure>,
}

/// Detached structural hierarchy, coordinates and metadata used for checked BioStructure construction.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BioStructureParts {
    inner: ck::BioStructureParts,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioStructureParts {
    /// Detected/declared structural input format.
    #[getter]
    fn input_format<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "BioCoordinateFormat", self.inner.input_format as i64)
    }
    fn __repr__(&self) -> String {
        format!("BioStructureParts({:?})", self.inner)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioStructure {
    /// Stored inter-atom connection annotations in source order.
    fn connections(&self) -> Vec<crate::canonical_bio_metadata::BioConnection> {
        self.inner
            .connections()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_bio_metadata::BioConnection { inner })
            .collect()
    }
    /// Stored cis-peptide annotations in source order.
    fn cispeps(&self) -> Vec<crate::canonical_bio_metadata::BioCisPep> {
        self.inner
            .cispeps()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_bio_metadata::BioCisPep { inner })
            .collect()
    }
    /// Modified-residue annotations in source order.
    fn mod_residues(&self) -> Vec<crate::canonical_bio_metadata::BioModRes> {
        self.inner
            .mod_residues()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_bio_metadata::BioModRes { inner })
            .collect()
    }
    /// Secondary-structure helix annotations.
    fn helices(&self) -> Vec<crate::canonical_bio_metadata::BioHelix> {
        self.inner
            .helices()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_bio_metadata::BioHelix { inner })
            .collect()
    }
    /// Secondary-structure sheet annotations.
    fn sheets(&self) -> Vec<crate::canonical_bio_metadata::BioSheet> {
        self.inner
            .sheets()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_bio_metadata::BioSheet { inner })
            .collect()
    }
    /// Noncrystallographic-symmetry transformations.
    fn ncs_operators(&self) -> Vec<crate::canonical_bio_metadata::BioNcsOperator> {
        self.inner
            .ncs_operators()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_bio_metadata::BioNcsOperator { inner })
            .collect()
    }
    /// Biological assembly definitions and their generators.
    fn assemblies(&self) -> Vec<crate::canonical_bio_metadata::BioAssembly> {
        self.inner
            .assemblies()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_bio_metadata::BioAssembly { inner })
            .collect()
    }
    /// Experimental, crystallographic and bibliographic metadata.
    fn metadata(&self) -> crate::canonical_bio_metadata::BioMetadata {
        crate::canonical_bio_metadata::BioMetadata {
            inner: self.inner.metadata().clone(),
        }
    }
    /// Source-format state retained for structural roundtrips.
    fn source_state(&self) -> crate::canonical_bio_metadata::BioStructureSourceState {
        crate::canonical_bio_metadata::BioStructureSourceState {
            inner: self.inner.source_state().clone(),
        }
    }
    /// Construct a value from explicit detached parts; required structural consistency is checked at the public boundary.
    #[staticmethod]
    fn from_parts(py: Python<'_>, parts: &BioStructureParts) -> PyResult<Self> {
        ck::BioStructure::from_parts(parts.inner.clone())
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|error| structure_error(py, &error))
    }
    /// Validate detached structure parts before constructing a live structure.
    #[staticmethod]
    fn validate_parts(py: Python<'_>, parts: &BioStructureParts) -> PyResult<()> {
        ck::BioStructure::validate_parts(&parts.inner).map_err(|error| structure_error(py, &error))
    }
    /// Return detached structure parts, including coordinates, hierarchy and metadata.
    fn into_parts(&self) -> BioStructureParts {
        // Rust consumes an owned value; Python retains its receiver and passes
        // a public clone to that same facade method, preserving every block.
        BioStructureParts {
            inner: self.inner.as_ref().clone().into_parts(),
        }
    }
    /// Find an atom by its structural address and alternate-location request.
    #[pyo3(signature = (residue_id, name, request, element))]
    fn find_atom(
        &self,
        residue_id: u32,
        name: &AtomName,
        request: &AltLocRequest,
        element: Option<&Element>,
    ) -> Option<(usize, BioAtomRow)> {
        self.inner
            .find_atom(
                ck::BioResidueId::new(residue_id),
                name.inner,
                request.inner,
                element.map(|value| value.inner),
            )
            .map(|(id, _)| {
                let index = id.index();
                (
                    index,
                    BioAtomRow {
                        inner: Arc::clone(&self.inner),
                        index,
                    },
                )
            })
    }
    /// Find a residue atom using its name and alternate-location request.
    #[pyo3(signature = (residue_id, name, altloc))]
    fn atom_by_altloc(
        &self,
        py: Python<'_>,
        residue_id: u32,
        name: &AtomName,
        altloc: Option<&AltLocLabel>,
    ) -> PyResult<(usize, BioAtomRow)> {
        self.inner
            .atom_by_altloc(
                ck::BioResidueId::new(residue_id),
                name.inner,
                altloc.map(|label| label.inner),
            )
            .map(|(id, _)| {
                let index = id.index();
                (
                    index,
                    BioAtomRow {
                        inner: Arc::clone(&self.inner),
                        index,
                    },
                )
            })
            .map_err(|error| structure_error(py, &error))
    }
    /// Validate the stored structure/template and return its validation result; invalid data is not silently repaired.
    fn validate(&self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .validate()
            .map_err(|error| structure_error(py, &error))
    }
    /// Return the BIO coordinate block in structural atom-row order.
    fn coordinates(&self) -> BioCoordinateBlock {
        BioCoordinateBlock {
            inner: Arc::clone(&self.inner),
        }
    }
    /// Convert the selected BIO atoms and positions into a new Molecule using the conversion options. Uses the supplied configuration object.
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
    /// Convert the selected BIO atoms and positions into a new Molecule using the conversion options.
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

    /// Parse PDB text into a structural object, preserving hierarchy, coordinates and supported source metadata; this argument is text, not a file path.
    #[staticmethod]
    fn from_pdb(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::BioStructure::from_pdb(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| pdb_error(py, e))
    }
    /// Parse PDB text into a structural object, preserving hierarchy, coordinates and supported source metadata; this argument is text, not a file path. Uses the supplied configuration object.
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
    /// Parse mmCIF text into a structural object, preserving hierarchy, coordinates and supported metadata; this argument is text, not a file path.
    #[staticmethod]
    fn from_mmcif(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::BioStructure::from_mmcif(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| mmcif_error(py, e))
    }
    /// Parse structural text using the requested BioReadParams input format.
    #[staticmethod]
    fn from_text(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::BioStructure::from_text(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| read_error(py, e))
    }
    /// Parse structural text using the requested BioReadParams input format. Uses the supplied configuration object.
    #[staticmethod]
    fn from_text_with_params(py: Python<'_>, text: &str, params: &BioReadParams) -> PyResult<Self> {
        ck::BioStructure::from_text_with_params(text, &params.inner)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| read_error(py, e))
    }
    /// Read a PDB/mmCIF file from a filesystem path, selecting its format through the read options.
    #[staticmethod]
    fn read(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        ck::BioStructure::read(&expand_path(py, path)?)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| read_error(py, e))
    }
    /// Read a structural file from a filesystem path using an explicit input format.
    #[staticmethod]
    fn read_with_format(
        py: Python<'_>,
        path: PathBuf,
        #[pyo3(from_py_with = format_input)]
        #[gen_stub(override_type(
            type_repr = "BioCoordinateFormat | builtins.str | builtins.int"
        ))]
        format: u8,
    ) -> PyResult<Self> {
        ck::BioStructure::read_with_format(&expand_path(py, path)?, format_from_code(format)?)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| read_error(py, e))
    }
    /// Stored name of this value.
    fn name(&self) -> &str {
        self.inner.name()
    }
    /// Detected/declared structural input format.
    fn input_format(&self) -> String {
        format!("{:?}", self.inner.input_format())
    }
    /// Number of structural models.
    fn num_models(&self) -> usize {
        self.inner.num_models()
    }
    /// Number of chains.
    fn num_chains(&self) -> usize {
        self.inner.num_chains()
    }
    /// Number of residues.
    fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }
    /// Number of atoms in the graph or selected structure.
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    /// Number of molecular entities.
    fn num_entities(&self) -> usize {
        self.inner.num_entities()
    }
    /// Whether the structure contains an original-coordinate transformation.
    fn has_origx(&self) -> bool {
        self.inner.has_origx()
    }
    /// Original-coordinate transformation stored in the structure.
    fn origx(&self) -> BioTransform {
        BioTransform {
            inner: *self.inner.origx(),
        }
    }
    /// Identifier of the identity noncrystallographic-symmetry operator, when available.
    fn ncs_oper_identity_id(&self) -> Option<&str> {
        self.inner.ncs_oper_identity_id()
    }
    /// Experimental resolution in angstroms as recorded in the structure.
    fn resolution(&self) -> f64 {
        self.inner.resolution()
    }
    /// Source TER-record handling state.
    fn ter_status(&self) -> u8 {
        self.inner.ter_status()
    }
    /// Return the (x, y, z) position in angstroms for a local atom ID.
    fn atom_position(&self, atom: u32) -> Option<[f64; 3]> {
        self.inner.atom_position(ck::BioAtomId::new(atom))
    }
    /// Return the atom rows belonging to the specified residue.
    fn residue_atoms(&self, residue: u32) -> Option<Vec<BioAtomRow>> {
        let rows = self.inner.residue_atoms(ck::BioResidueId::new(residue))?;
        let start = self.inner.residues()[residue as usize].atom_span().start() as usize;
        Some(
            (start..start + rows.len())
                .map(|index| BioAtomRow {
                    inner: self.inner.clone(),
                    index,
                })
                .collect(),
        )
    }
    /// Crystallographic unit cell and space-group information.
    fn crystal(&self) -> Option<BioCrystalInfo> {
        self.inner.crystal().map(|_| BioCrystalInfo {
            inner: self.inner.clone(),
        })
    }
    /// Find the entity with the requested source entity identifier, or None if absent.
    fn find_entity(&self, source_id: &str) -> Option<(usize, BioEntityRow)> {
        self.inner.find_entity(source_id).map(|(id, _)| {
            (
                id.index(),
                BioEntityRow {
                    inner: self.inner.clone(),
                    index: id.index(),
                },
            )
        })
    }
    /// Find the entity associated with the requested subchain, or None if absent.
    fn find_entity_of_subchain(&self, subchain: &str) -> Option<(usize, BioEntityRow)> {
        self.inner.find_entity_of_subchain(subchain).map(|(id, _)| {
            (
                id.index(),
                BioEntityRow {
                    inner: self.inner.clone(),
                    index: id.index(),
                },
            )
        })
    }
    /// Return structural model rows in stored order.
    fn models(&self) -> Vec<BioModelRow> {
        (0..self.inner.num_models())
            .map(|index| BioModelRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    /// Return chain rows/references in stored hierarchy order.
    fn chains(&self) -> Vec<BioChainRow> {
        (0..self.inner.num_chains())
            .map(|index| BioChainRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    /// Return residue rows/references in stored hierarchy order.
    fn residues(&self) -> Vec<BioResidueRow> {
        (0..self.inner.num_residues())
            .map(|index| BioResidueRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    /// Return atom rows in stored graph/hierarchy order.
    fn atoms(&self) -> Vec<BioAtomRow> {
        (0..self.inner.num_atoms())
            .map(|index| BioAtomRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    /// Return molecular entity rows in stored order.
    fn entities(&self) -> Vec<BioEntityRow> {
        (0..self.inner.num_entities())
            .map(|index| BioEntityRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    /// Return a Protein view over the structure protein hierarchy.
    fn protein(&self, py: Python<'_>) -> PyResult<Protein> {
        self.inner
            .protein()
            .map(|inner| Protein {
                inner: Arc::new(inner),
            })
            .map_err(|e| {
                annotate(
                    py,
                    crate::canonical_registered_errors::ProteinProjectionError::new_err(
                        e.to_string(),
                    ),
                    "Projection",
                    &e,
                )
            })
    }
    /// Return mmCIF text using the selected category/formatting options; does not write a file. Uses the supplied configuration object.
    fn to_mmcif_with_params(
        &self,
        py: Python<'_>,
        params: &BioMmcifWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_mmcif_with_params(&params.inner)
            .map_err(|e| write_error(py, e))
    }
    /// Write MMCIF output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
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
    /// Return PDB text using the selected formatting options; does not write a file.
    fn to_pdb(&self, py: Python<'_>) -> PyResult<String> {
        self.inner.to_pdb().map_err(|e| pdb_write_error(py, e))
    }
    /// Return PDB text using the selected formatting options; does not write a file. Uses the supplied configuration object.
    fn to_pdb_with_params(&self, py: Python<'_>, params: &BioPdbWriteParams) -> PyResult<String> {
        self.inner
            .to_pdb_with_params(&params.inner)
            .map_err(|e| pdb_write_error(py, e))
    }
    /// Write PDB output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    fn write_pdb(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        self.inner
            .write_pdb(&expand_path(py, path)?)
            .map_err(|e| pdb_write_error(py, e))
    }
    /// Write PDB output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
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
    /// Return mmCIF text using the selected category/formatting options; does not write a file.
    fn to_mmcif(&self, py: Python<'_>) -> PyResult<String> {
        self.inner.to_mmcif().map_err(|e| write_error(py, e))
    }
    /// Write MMCIF output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    fn write_mmcif(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        self.inner
            .write_mmcif(&expand_path(py, path)?)
            .map_err(|e| write_error(py, e))
    }
    /// Return a new structure restricted to the selection; leave the source unchanged.
    fn with_selection(&self, py: Python<'_>, selection: &BioSelection) -> PyResult<Self> {
        self.inner
            .with_selection(&selection.inner)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| operation_error(py, e))
    }
    /// Retain only the selected structure rows in place, updating hierarchy, coordinates and references together.
    fn retain_selection_(&mut self, py: Python<'_>, selection: &BioSelection) -> PyResult<()> {
        Arc::make_mut(&mut self.inner)
            .retain_selection_(&selection.inner)
            .map_err(|e| operation_error(py, e))
    }
    /// Return a new structure with every atom position translated by the supplied vector; leave source coordinates unchanged.
    fn with_translated_coordinates(&self, py: Python<'_>, offset: [f64; 3]) -> PyResult<Self> {
        self.inner
            .with_translated_coordinates(offset)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| operation_error(py, e))
    }
    /// Translate all structural atom positions in place by the supplied vector.
    fn translate_(&mut self, py: Python<'_>, offset: [f64; 3]) -> PyResult<()> {
        Arc::make_mut(&mut self.inner)
            .translate_(offset)
            .map_err(|e| operation_error(py, e))
    }
    /// Return the local atom identifiers matching the selection.
    fn selected_atom_ids(&self, py: Python<'_>, selection: &BioSelection) -> PyResult<Vec<usize>> {
        self.inner
            .selected_atom_ids(&selection.inner)
            .map(|x| x.into_iter().map(|id| id.index()).collect())
            .map_err(|e| {
                annotate(
                    py,
                    crate::canonical_registered_errors::BioSelectionMatchError::new_err(
                        e.to_string(),
                    ),
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

/// Read-only structural model row with a local model ID, source model number and child chain span.
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
    /// Zero-based identifier in the owning object; not a PDB serial or residue number.
    fn id(&self) -> usize {
        self.index
    }
    /// Model number recorded in the source structure.
    fn source_model_number(&self) -> Option<i32> {
        self.inner.models()[self.index].source_model_number()
    }
    /// Return chain rows/references in stored hierarchy order.
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

/// Read-only chain row with local parent/entity identifiers, source identifiers and residue/atom spans.
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
    /// Zero-based identifier in the owning object; not a PDB serial or residue number.
    fn id(&self) -> usize {
        self.index
    }
    /// Local identifier of the parent model.
    fn model_id(&self) -> usize {
        self.row().model_id().index()
    }
    /// Local identifier of the associated molecular entity, when present.
    fn entity_id(&self) -> Option<usize> {
        self.row().entity_id().map(|id| id.index())
    }
    /// Classification/discriminant of this value as defined by its owning type.
    #[gen_stub(override_return_type(type_repr = "ChainKind"))]
    fn kind<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        kind_member(py, "ChainKind", self.row().kind())
    }
    /// Source-format identifiers retained separately from local row identifiers.
    fn source(&self) -> ChainSourceIds {
        ChainSourceIds {
            inner: self.row().source().clone(),
        }
    }
    /// Return residue rows/references in stored hierarchy order.
    fn residues(&self) -> Vec<BioResidueRow> {
        let span = self.row().residue_span();
        (span.start() as usize..span.end() as usize)
            .map(|index| BioResidueRow {
                inner: self.inner.clone(),
                index,
            })
            .collect()
    }
    /// Return atom rows in stored graph/hierarchy order.
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

/// Read-only residue row with chain ownership, residue classification and retained source sequence identifiers.
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
    /// Zero-based identifier in the owning object; not a PDB serial or residue number.
    fn id(&self) -> usize {
        self.index
    }
    /// Local identifier of the parent chain.
    fn chain_id(&self) -> usize {
        self.row().chain_id().index()
    }
    /// Stored name of this value.
    fn name(&self) -> String {
        self.row().name().as_str().to_owned()
    }
    /// Classification/discriminant of this value as defined by its owning type.
    #[gen_stub(override_return_type(type_repr = "ResidueKind"))]
    fn kind<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        kind_member(py, "ResidueKind", self.row().kind())
    }
    /// Structural entity category of the residue.
    #[gen_stub(override_return_type(type_repr = "EntityKind"))]
    fn entity_kind<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        kind_member(py, "EntityKind", self.row().entity_kind())
    }
    /// SIFTS/UniProt residue mapping metadata, when present.
    fn sifts_unp(&self) -> crate::canonical_bio_values::BioSiftsUnpResidue {
        crate::canonical_bio_values::BioSiftsUnpResidue {
            inner: self.row().sifts_unp(),
        }
    }
    /// Return the associated residue dictionary information, including its classification and sequence codes.
    fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: ck::find_residue_info(self.row().name().as_str()),
        }
    }
    /// Stored residue/atom code; its interpretation is defined by the owning value type.
    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(
            py,
            "ResidueCode",
            i64::from(ck::residue_code(self.row().name().as_str()).as_u16()),
        )
    }
    /// Source-format identifiers retained separately from local row identifiers.
    fn source(&self) -> ResidueSourceIds {
        ResidueSourceIds {
            inner: self.row().source().clone(),
        }
    }
    /// Return atom rows in stored graph/hierarchy order.
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

/// Read-only structural atom row with residue ownership, chemical identity, Cartesian position, occupancy and displacement factor.
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
    /// Atom-site calculation/refinement flag.
    #[gen_stub(override_return_type(type_repr = "BioCalcFlag"))]
    fn calc_flag<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        kind_member(py, "BioCalcFlag", self.row().calc_flag())
    }
    /// Zero-based identifier in the owning object; not a PDB serial or residue number.
    fn id(&self) -> usize {
        self.index
    }
    /// Local identifier of the parent residue.
    fn residue_id(&self) -> usize {
        self.row().residue_id().index()
    }
    /// Stored name of this value.
    fn name(&self) -> String {
        self.row().name().as_str().trim().to_owned()
    }
    /// Chemical element as an Element value.
    fn element(&self) -> Element {
        Element::from_inner(self.row().element())
    }
    /// Chemical element symbol, such as "C" or "Cl".
    fn element_symbol(&self) -> &'static str {
        self.row().element().symbol()
    }
    /// Cartesian position as an (x, y, z) tuple in angstroms.
    fn position(&self) -> Option<(f64, f64, f64)> {
        self.inner
            .atom_position(ck::BioAtomId::new(self.index as u32))
            .map(|[x, y, z]| (x, y, z))
    }
    /// Alternate-location label, or None when the site is not alternate.
    fn altloc(&self) -> Option<String> {
        self.row()
            .altloc()
            .map(|x| char::from(x.value()).to_string())
    }
    /// Crystallographic atom-site occupancy.
    fn occupancy(&self) -> f64 {
        self.row().occupancy()
    }
    /// Isotropic atomic displacement B factor in square angstroms.
    fn b_iso(&self) -> f64 {
        self.row().b_iso()
    }
    /// Formal charge in units of the elementary charge.
    fn formal_charge(&self) -> i8 {
        self.row().formal_charge()
    }
    /// Source-format identifiers retained separately from local row identifiers.
    fn source(&self) -> AtomSourceIds {
        AtomSourceIds {
            inner: *self.row().source(),
        }
    }
}

/// Read-only molecular entity metadata, including polymer kind, sequence, subchains and database references.
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
    /// Zero-based identifier in the owning object; not a PDB serial or residue number.
    fn id(&self) -> usize {
        self.index
    }
    /// Source-format identifiers retained separately from local row identifiers.
    fn source(&self) -> EntitySourceIds {
        EntitySourceIds {
            inner: self.row().source().clone(),
        }
    }
    /// Classification/discriminant of this value as defined by its owning type.
    #[gen_stub(override_return_type(type_repr = "EntityKind"))]
    fn kind<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        kind_member(py, "EntityKind", self.row().kind())
    }
    /// Polymer category, such as peptide or nucleic acid.
    #[gen_stub(override_return_type(type_repr = "PolymerKind"))]
    fn polymer_kind<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        kind_member(py, "PolymerKind", self.row().polymer_kind())
    }
    /// External database references for the entity in stored order.
    fn dbrefs(&self) -> Vec<crate::canonical_bio_values::BioEntityDbRef> {
        self.row()
            .dbrefs()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_bio_values::BioEntityDbRef { inner })
            .collect()
    }
    /// Complete entity sequence, including residues not resolved in coordinates.
    fn full_sequence(&self) -> Vec<String> {
        self.row().full_sequence().to_vec()
    }
    /// Source subchain identifiers associated with the entity.
    fn subchains(&self) -> Vec<String> {
        self.row().subchains().to_vec()
    }
    fn __len__(&self) -> usize {
        self.row().full_sequence().len()
    }
}

/// Source chain identifiers: author chain name and mmCIF label asymmetric-unit ID.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ChainSourceIds {
    inner: ck::ChainSourceIds,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ChainSourceIds {
    /// Author-provided chain identifier from the input structure.
    fn auth_chain_id(&self) -> Option<String> {
        self.inner.auth_chain_id().map(|x| x.as_str().to_owned())
    }
    /// mmCIF label_asym_id identifying the asymmetric-unit chain.
    fn label_asym_id(&self) -> Option<&str> {
        self.inner.label_asym_id()
    }
}
/// Source atom-site identifiers, including the PDB serial number.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomSourceIds {
    inner: ck::AtomSourceIds,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl AtomSourceIds {
    /// Atom serial number from the source file.
    fn serial(&self) -> Option<i32> {
        self.inner.serial().map(|x| x.value())
    }
}
/// Source entity identifier retained separately from local entity IDs.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct EntitySourceIds {
    inner: ck::EntitySourceIds,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl EntitySourceIds {
    /// Entity identifier used in the source structure.
    fn source_entity_id(&self) -> &str {
        self.inner.source_entity_id()
    }
}
/// Source residue sequence number plus insertion code; not a local residue row ID.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct PdbSeqId {
    pub(crate) inner: ck::PdbSeqId,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PdbSeqId {
    /// Construct a PdbSeqId value from the supplied inputs.
    #[new]
    #[pyo3(signature=(seq_num, ins_code=None))]
    fn new(seq_num: i32, ins_code: Option<u8>) -> Self {
        Self {
            inner: ck::PdbSeqId::new(seq_num, ins_code),
        }
    }
    /// Residue sequence number, when present.
    fn seq_num(&self) -> i32 {
        self.inner.seq_num()
    }
    /// Residue insertion code retained with the sequence number.
    fn ins_code(&self) -> Option<String> {
        self.inner.ins_code().map(|x| char::from(x).to_string())
    }
}
/// Source residue identifiers: author number/insertion code, label number, subchain and entity.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ResidueSourceIds {
    inner: ck::ResidueSourceIds,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ResidueSourceIds {
    /// Author residue sequence identifier, including insertion code.
    fn seq_id(&self) -> Option<PdbSeqId> {
        self.inner.seq_id().map(|inner| PdbSeqId { inner })
    }
    /// mmCIF label residue sequence number, when present.
    fn label_seq_id(&self) -> Option<i32> {
        self.inner.label_seq_id()
    }
    /// Source subchain identifier.
    fn subchain_id(&self) -> Option<&str> {
        self.inner.subchain_id()
    }
    /// mmCIF label entity identifier.
    fn label_entity_id(&self) -> Option<&str> {
        self.inner.label_entity_id()
    }
}

/// Protein-oriented structural hierarchy with chain, residue and atom references. Parsing accepts text through from_* and files through read(); selection preserves the source unless the method ends in an underscore.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct Protein {
    inner: Arc<ck::Protein>,
}
/// Counts of chains, residues and atoms retained by a protein selection.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ProteinSelectionSummary {
    inner: ck::ProteinSelectionSummary,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ProteinSelectionSummary {
    /// Return chain rows/references in stored hierarchy order.
    #[getter]
    fn chains(&self) -> usize {
        self.inner.chains
    }
    /// Return residue rows/references in stored hierarchy order.
    #[getter]
    fn residues(&self) -> usize {
        self.inner.residues
    }
    /// Return the number of atoms in the structural selection.
    #[getter]
    fn atoms(&self) -> usize {
        self.inner.atoms
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl Protein {
    /// Return the underlying structural representation without parsing or generating coordinates.
    fn as_bio_structure(&self) -> BioStructure {
        BioStructure {
            inner: Arc::new(self.inner.as_bio_structure().clone()),
        }
    }
    /// Return the structural representation as a BioStructure value.
    fn into_bio_structure(&self) -> BioStructure {
        BioStructure {
            inner: Arc::new(self.inner.as_ref().clone().into_bio_structure()),
        }
    }
    /// Return the numbers of chains, residues and atoms in the protein selection.
    fn selection_summary(&self) -> ProteinSelectionSummary {
        ProteinSelectionSummary {
            inner: self.inner.selection_summary(),
        }
    }
    /// Return the chain reference for a local chain identifier.
    fn chain(&self, index: usize) -> Option<ProteinChainRef> {
        self.inner.chain(index).map(|row| ProteinChainRef {
            inner: self.inner.clone(),
            index: row.id().index(),
        })
    }
    /// Convert the selected BIO atoms and positions into a new Molecule using the conversion options. Uses the supplied configuration object.
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
    /// Convert the selected BIO atoms and positions into a new Molecule using the conversion options.
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

    /// Parse PDB text into a structural object, preserving hierarchy, coordinates and supported source metadata; this argument is text, not a file path.
    #[staticmethod]
    fn from_pdb(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Protein::from_pdb(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    /// Parse PDB text into a structural object, preserving hierarchy, coordinates and supported source metadata; this argument is text, not a file path. Uses the supplied configuration object.
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
    /// Parse mmCIF text into a structural object, preserving hierarchy, coordinates and supported metadata; this argument is text, not a file path.
    #[staticmethod]
    fn from_mmcif(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Protein::from_mmcif(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    /// Parse structural text using the requested BioReadParams input format.
    #[staticmethod]
    fn from_text(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Protein::from_text(text)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    /// Parse structural text using the requested BioReadParams input format. Uses the supplied configuration object.
    #[staticmethod]
    fn from_text_with_params(py: Python<'_>, text: &str, params: &BioReadParams) -> PyResult<Self> {
        ck::Protein::from_text_with_params(text, &params.inner)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    /// Read a PDB/mmCIF file from a filesystem path, selecting its format through the read options.
    #[staticmethod]
    fn read(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        ck::Protein::read(&expand_path(py, path)?)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    /// Read a structural file from a filesystem path using an explicit input format.
    #[staticmethod]
    fn read_with_format(
        py: Python<'_>,
        path: PathBuf,
        #[pyo3(from_py_with = format_input)]
        #[gen_stub(override_type(
            type_repr = "BioCoordinateFormat | builtins.str | builtins.int"
        ))]
        format: u8,
    ) -> PyResult<Self> {
        ck::Protein::read_with_format(&expand_path(py, path)?, format_from_code(format)?)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| protein_error(py, e))
    }
    /// Detected/declared structural input format.
    fn input_format(&self) -> String {
        format!("{:?}", self.inner.input_format())
    }
    /// Number of structural models.
    fn num_models(&self) -> usize {
        self.inner.num_models()
    }
    /// Number of chains.
    fn num_chains(&self) -> usize {
        self.inner.num_chains()
    }
    /// Number of residues.
    fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }
    /// Number of atoms in the graph or selected structure.
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    /// Return chain rows/references in stored hierarchy order.
    fn chains(&self) -> Vec<ProteinChainRef> {
        self.inner
            .chains()
            .map(|x| ProteinChainRef {
                inner: self.inner.clone(),
                index: x.id().index(),
            })
            .collect()
    }
    /// Return residue rows/references in stored hierarchy order.
    fn residues(&self) -> Vec<ProteinResidueRef> {
        self.inner
            .residues()
            .map(|x| ProteinResidueRef {
                inner: self.inner.clone(),
                index: x.id().index(),
            })
            .collect()
    }
    /// Return atom rows in stored graph/hierarchy order.
    fn atoms(&self) -> Vec<ProteinAtomRef> {
        self.inner
            .atoms()
            .map(|x| ProteinAtomRef {
                inner: self.inner.clone(),
                index: x.id().index(),
            })
            .collect()
    }
    /// Return a new structure restricted to the selection; leave the source unchanged.
    fn with_selection(&self, py: Python<'_>, selection: &BioSelection) -> PyResult<Self> {
        self.inner
            .with_selection(&selection.inner)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| operation_error(py, e))
    }
    /// Retain only the selected protein rows in place, updating hierarchy and coordinates together.
    fn retain_selection_(&mut self, py: Python<'_>, selection: &BioSelection) -> PyResult<()> {
        Arc::make_mut(&mut self.inner)
            .retain_selection_(&selection.inner)
            .map_err(|e| operation_error(py, e))
    }
    /// Return a new protein with translated atom coordinates, preserving hierarchy and metadata.
    fn with_translated_coordinates(&self, py: Python<'_>, offset: [f64; 3]) -> PyResult<Self> {
        self.inner
            .with_translated_coordinates(offset)
            .map(|x| Self { inner: Arc::new(x) })
            .map_err(|e| operation_error(py, e))
    }
    /// Translate protein atom coordinates in place while retaining hierarchy and metadata.
    fn translate_(&mut self, py: Python<'_>, offset: [f64; 3]) -> PyResult<()> {
        Arc::make_mut(&mut self.inner)
            .translate_(offset)
            .map_err(|e| operation_error(py, e))
    }
    /// Return the local atom identifiers matching the selection.
    fn selected_atom_ids(&self, py: Python<'_>, selection: &BioSelection) -> PyResult<Vec<usize>> {
        self.inner
            .selected_atom_ids(&selection.inner)
            .map(|x| x.into_iter().map(|id| id.index()).collect())
            .map_err(|e| {
                annotate(
                    py,
                    crate::canonical_registered_errors::BioSelectionMatchError::new_err(
                        e.to_string(),
                    ),
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

/// Read-only reference to a protein chain and its residue/atom hierarchy.
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
    /// Return the detached row represented by this hierarchy reference.
    fn row(&self) -> BioChainRow {
        BioChainRow {
            inner: Arc::new(self.inner.as_bio_structure().clone()),
            index: self.view().id().index(),
        }
    }
    /// Source-format identifiers retained separately from local row identifiers.
    fn source(&self) -> ChainSourceIds {
        ChainSourceIds {
            inner: self.view().source().clone(),
        }
    }
    /// Zero-based identifier in the owning object; not a PDB serial or residue number.
    fn id(&self) -> usize {
        self.index
    }
    /// Classification/discriminant of this value as defined by its owning type.
    #[gen_stub(override_return_type(type_repr = "ChainKind"))]
    fn kind<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        kind_member(py, "ChainKind", self.view().kind())
    }
    /// Return residue rows/references in stored hierarchy order.
    fn residues(&self) -> Vec<ProteinResidueRef> {
        self.view()
            .residues()
            .map(|x| ProteinResidueRef {
                inner: self.inner.clone(),
                index: x.id().index(),
            })
            .collect()
    }
    /// Return atom rows in stored graph/hierarchy order.
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

/// Read-only protein residue reference with dictionary classification, sequence codes and atom access.
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
    /// Return the detached row represented by this hierarchy reference.
    fn row(&self) -> BioResidueRow {
        BioResidueRow {
            inner: Arc::new(self.inner.as_bio_structure().clone()),
            index: self.view().id().index(),
        }
    }
    /// Return the chain reference for a local chain identifier.
    fn chain(&self) -> ProteinChainRef {
        ProteinChainRef {
            inner: self.inner.clone(),
            index: self.view().chain().id().index(),
        }
    }
    /// Zero-based identifier in the owning object; not a PDB serial or residue number.
    fn id(&self) -> usize {
        self.index
    }
    /// Stored name of this value.
    fn name(&self) -> String {
        self.view().name().as_str().to_owned()
    }
    /// Classification/discriminant of this value as defined by its owning type.
    #[gen_stub(override_return_type(type_repr = "ResidueKind"))]
    fn kind<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        kind_member(py, "ResidueKind", self.view().kind())
    }
    /// Return the associated residue dictionary information, including its classification and sequence codes.
    fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: self.view().info(),
        }
    }
    /// Stored residue/atom code; its interpretation is defined by the owning value type.
    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "ResidueCode", i64::from(self.view().code().as_u16()))
    }
    /// One-letter residue code recorded by the residue dictionary.
    fn one_letter_code(&self) -> String {
        self.view().one_letter_code().to_string()
    }
    /// Residue code suitable for FASTA sequence output.
    fn fasta_code(&self) -> String {
        self.view().fasta_code().to_string()
    }
    /// One-letter code of the corresponding standard parent residue.
    fn canonical_one_letter_code(&self) -> Option<String> {
        self.view()
            .info()
            .canonical_one_letter_code()
            .map(|x| x.to_string())
    }
    /// Standard parent residue code for a modified residue, when known.
    #[gen_stub(override_return_type(type_repr = "typing.Optional[ResidueCode]"))]
    fn parent_standard_code<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.view()
            .info()
            .parent_standard_code()
            .map(|x| enum_member(py, "ResidueCode", i64::from(x.as_u16())))
            .transpose()
    }
    /// Whether the residue dictionary identifies a modified amino acid.
    fn is_modified_amino_acid(&self) -> bool {
        self.view().info().is_modified_amino_acid()
    }
    /// Whether the residue belongs to the standard residue set.
    fn is_standard(&self) -> bool {
        self.view().is_standard()
    }
    /// Return atom rows in stored graph/hierarchy order.
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

/// Read-only protein atom reference with parent residue and Cartesian position.
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
    /// Return the detached row represented by this hierarchy reference.
    fn row(&self) -> BioAtomRow {
        BioAtomRow {
            inner: Arc::new(self.inner.as_bio_structure().clone()),
            index: self.view().id().index(),
        }
    }
    /// Alternate-location label, or None when the site is not alternate.
    fn altloc(&self) -> Option<String> {
        self.view()
            .altloc()
            .map(|label| char::from(label.value()).to_string())
    }
    /// Return the parent residue reference.
    fn residue(&self) -> ProteinResidueRef {
        ProteinResidueRef {
            inner: self.inner.clone(),
            index: self.view().residue().id().index(),
        }
    }
    /// Zero-based identifier in the owning object; not a PDB serial or residue number.
    fn id(&self) -> usize {
        self.index
    }
    /// Stored name of this value.
    fn name(&self) -> String {
        self.view().name().as_str().trim().to_owned()
    }
    /// Chemical element as an Element value.
    fn element(&self) -> Element {
        Element::from_inner(self.view().element())
    }
    /// Chemical element symbol, such as "C" or "Cl".
    fn element_symbol(&self) -> &'static str {
        self.view().element().symbol()
    }
    /// Atomic number (proton count); zero denotes a dummy atom.
    fn atomic_number(&self) -> u8 {
        self.view().element().atomic_number()
    }
    /// Cartesian position as an (x, y, z) tuple in angstroms.
    fn position(&self) -> (f64, f64, f64) {
        let [x, y, z] = self.view().position();
        (x, y, z)
    }
}

/// Writable configuration for structural BIO-to-Molecule conversion.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct BioMoleculeParams {
    inner: ck::BioMoleculeParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BioMoleculeParams {
    /// Configure structural BIO-to-Molecule conversion; omitted fields use the defaults shown in the signature.
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
    /// Apply to a new molecule and return the result: perform the selected chemical sanitization stages. The source molecule is unchanged.
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    /// Whether removable explicit hydrogens are removed during input conversion.
    #[getter]
    fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    /// Bit mask controlling BIO-to-molecule conversion behavior.
    #[getter]
    fn flavor(&self) -> u32 {
        self.inner.flavor
    }
    /// Whether spatially close atoms are considered when building molecular bonds.
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
    crate::canonical_bio_metadata::register(module)?;
    module.add_class::<crate::canonical_error_values::BioPdbReadStage>()?;
    module.add_class::<crate::canonical_error_values::BioMmcifReadStage>()?;
    crate::canonical_error_accessors::attach(
        module.py().get_type::<BioPdbReadError>().as_any(),
        &[
            ("stage", "_stage"),
            ("line_number", "_line_number"),
            ("record_tag", "_record_tag"),
        ],
    )?;
    crate::canonical_error_accessors::attach(
        module.py().get_type::<BioMmcifReadError>().as_any(),
        &[("stage", "_stage")],
    )?;
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
    formats.setattr("__doc__", "Biological coordinate input format. Detect selects format recognition; explicit values select PDB, mmCIF, mmJSON or chemical-component input.")?;
    module.add("BioCoordinateFormat", formats)?;
    module.add_class::<BioMoleculeParams>()?;
    module.add("BioMoleculeError", py.get_type::<BioMoleculeError>())?;
    module.add(
        "BioMoleculeConversionError",
        py.get_type::<BioMoleculeConversionError>(),
    )?;
    module.add("BioReadError", py.get_type::<BioReadError>())?;
    module.add("BioStructureError", py.get_type::<BioStructureError>())?;
    module.add_class::<AtomName>()?;
    module.add_class::<AltLocLabel>()?;
    module.add_class::<AltLocRequest>()?;
    module.add_class::<BioCoordinateBlock>()?;
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
    module.add_class::<BioTransform>()?;
    module.add_class::<BioCrystalInfo>()?;
    module.add_class::<BioStructure>()?;
    module.add_class::<BioStructureParts>()?;
    module.add_class::<Protein>()?;
    module.add_class::<ProteinSelectionSummary>()?;
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
