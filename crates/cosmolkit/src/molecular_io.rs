//! Canonical molecular-format construction and thin IO error projection.
use crate::{Molecule, OperationError, SdfError, SdfReadParams, SdfRecord};
use std::fmt;
use std::fs::File;
use std::path::PathBuf;

pub use cosmolkit_io::{
    Mol2PostError, Mol2ReadError, Mol2Type, XyzReadError, XyzWriteError, XyzWriteParams,
};

pub use cosmolkit_io::{MolBlockWriteParams, MolCoordinateSelection, MolWriteError, SdfFormat};

#[derive(Debug)]
pub enum MolecularIoError {
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    XyzRead(XyzReadError),
    XyzWrite(XyzWriteError),
    MolWrite(MolWriteError),
    Mol2Read(Mol2ReadError),
    Mol2Post(Mol2PostError),
    Sdf(SdfError),
    Construction(OperationError),
    NoRecord {
        format: &'static str,
    },
    Parameter {
        name: &'static str,
        detail: &'static str,
    },
}
impl fmt::Display for MolecularIoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "{}: {}", path.display(), source),
            Self::XyzRead(e) => write!(f, "XYZ read: {e}"),
            Self::XyzWrite(e) => write!(f, "XYZ write: {e}"),
            Self::MolWrite(e) => e.fmt(f),
            Self::Mol2Read(e) => write!(f, "MOL2 read: {e}"),
            Self::Mol2Post(e) => write!(f, "MOL2 finalization: {e}"),
            Self::Sdf(e) => e.fmt(f),
            Self::Construction(e) => e.fmt(f),
            Self::NoRecord { format } => write!(f, "{format} found no molecule record"),
            Self::Parameter { name, detail } => write!(f, "{name}: {detail}"),
        }
    }
}
impl std::error::Error for MolecularIoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            Self::XyzRead(e) => Some(e),
            Self::XyzWrite(e) => Some(e),
            Self::MolWrite(e) => Some(e),
            Self::Mol2Read(e) => Some(e),
            Self::Mol2Post(e) => Some(e),
            Self::Sdf(e) => Some(e),
            Self::Construction(e) => Some(e),
            Self::NoRecord { .. } | Self::Parameter { .. } => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mol2ReadParams {
    pub sanitize: bool,
    pub remove_hydrogens: bool,
    pub variant: Mol2Type,
    pub cleanup_substructures: bool,
}
impl Default for Mol2ReadParams {
    fn default() -> Self {
        Self {
            sanitize: true,
            remove_hydrogens: true,
            variant: Mol2Type::Corina,
            cleanup_substructures: true,
        }
    }
}

pub(crate) fn expanded_path(path: &str) -> Result<PathBuf, MolecularIoError> {
    if path == "~" || path.starts_with("~/") {
        let home = std::env::var_os("HOME").ok_or(MolecularIoError::Parameter {
            name: "path",
            detail: "cannot expand '~': HOME is not set",
        })?;
        let mut result = PathBuf::from(home);
        if let Some(rest) = path.strip_prefix("~/") {
            result.push(rest);
        }
        Ok(result)
    } else {
        Ok(PathBuf::from(path))
    }
}
pub(crate) fn open_file(path: &str) -> Result<File, MolecularIoError> {
    let path = expanded_path(path)?;
    File::open(&path).map_err(|source| MolecularIoError::Io { path, source })
}
fn read_text(path: &str) -> Result<String, MolecularIoError> {
    use std::io::Read;
    let mut file = open_file(path)?;
    let mut text = String::new();
    file.read_to_string(&mut text)
        .map_err(|source| MolecularIoError::Io {
            path: PathBuf::from(path),
            source,
        })?;
    Ok(text)
}
pub(crate) fn write_text(path: &str, text: &str) -> Result<(), MolecularIoError> {
    let path = expanded_path(path)?;
    std::fs::write(&path, text).map_err(|source| MolecularIoError::Io { path, source })
}

impl Molecule {
    /// Source string projection of an atom property; absent atom or key returns None.
    pub fn atom_property_string(
        &self,
        id: crate::AtomId,
        key: &str,
    ) -> Result<Option<String>, crate::PropertyStringError> {
        self.atom(id)
            .and_then(|atom| atom.prop(key))
            .map(cosmolkit_core::property_value_to_string)
            .transpose()
    }
    /// Source string projection of a bond property; absent bond or key returns None.
    pub fn bond_property_string(
        &self,
        id: crate::BondId,
        key: &str,
    ) -> Result<Option<String>, crate::PropertyStringError> {
        self.bond(id)
            .and_then(|bond| bond.prop(key))
            .map(cosmolkit_core::property_value_to_string)
            .transpose()
    }
    pub fn from_xyz_block(text: &str) -> Result<Self, MolecularIoError> {
        let (topology, coordinates, properties) =
            cosmolkit_io::read_xyz_detached(text).map_err(MolecularIoError::XyzRead)?;
        Self::from_validated_parts(topology, coordinates, properties)
            .map_err(MolecularIoError::Construction)
    }
    pub fn read_xyz(path: &str) -> Result<Self, MolecularIoError> {
        Self::from_xyz_block(&read_text(path)?)
    }
    pub fn to_xyz(&self) -> Result<String, MolecularIoError> {
        self.to_xyz_with_params(&XyzWriteParams::default())
    }
    pub fn to_xyz_with_params(&self, params: &XyzWriteParams) -> Result<String, MolecularIoError> {
        cosmolkit_io::write_xyz_detached_with_params(
            self.topology(),
            self.coordinate_block_runtime(),
            self.properties(),
            *params,
        )
        .map_err(MolecularIoError::XyzWrite)
    }
    pub fn write_xyz(&self, path: &str) -> Result<(), MolecularIoError> {
        self.write_xyz_with_params(path, &XyzWriteParams::default())
    }
    pub fn write_xyz_with_params(
        &self,
        path: &str,
        params: &XyzWriteParams,
    ) -> Result<(), MolecularIoError> {
        write_text(path, &self.to_xyz_with_params(params)?)
    }
    pub fn from_mol(text: &str) -> Result<Self, SdfError> {
        Self::from_mol_with_params(text, &SdfReadParams::default())
    }
    pub fn from_mol_with_params(text: &str, params: &SdfReadParams) -> Result<Self, SdfError> {
        let parsed = cosmolkit_io::read_mol_graph_record_detached_with_params(
            text,
            crate::sdf_supplier::data_params(params),
        )?;
        SdfRecord::from_parsed(parsed, params, 0)?
            .molecule()
            .cloned()
            .map_err(|_| SdfError::QueryRecord)
    }
    pub fn read_mol(path: &str) -> Result<Self, MolecularIoError> {
        Self::read_mol_with_params(path, &SdfReadParams::default())
    }
    pub fn read_mol_with_params(
        path: &str,
        params: &SdfReadParams,
    ) -> Result<Self, MolecularIoError> {
        Self::from_mol_with_params(&read_text(path)?, params).map_err(MolecularIoError::Sdf)
    }
    pub fn read_sdf(path: &str) -> Result<Self, MolecularIoError> {
        Self::read_sdf_with_params(path, &SdfReadParams::default())
    }
    pub fn read_sdf_with_params(
        path: &str,
        params: &SdfReadParams,
    ) -> Result<Self, MolecularIoError> {
        let mut reader = crate::SdfRecordStream::open_with_params(path, params)?;
        let record = reader
            .next_record()
            .map_err(MolecularIoError::Sdf)?
            .ok_or(MolecularIoError::NoRecord { format: "SDF" })?;
        record.molecule().cloned().map_err(MolecularIoError::Sdf)
    }
    pub fn from_mol2(text: &str) -> Result<Self, MolecularIoError> {
        Self::from_mol2_with_params(text, &Mol2ReadParams::default())
    }
    pub fn from_mol2_with_params(
        text: &str,
        params: &Mol2ReadParams,
    ) -> Result<Self, MolecularIoError> {
        let record = cosmolkit_io::read_mol2_detached_with_params(
            text,
            cosmolkit_io::Mol2ReadParams {
                variant: params.variant,
                cleanup_substructures: params.cleanup_substructures,
            },
        )
        .map_err(MolecularIoError::Mol2Read)?
        .ok_or(MolecularIoError::NoRecord { format: "MOL2" })?;
        let mut record =
            cosmolkit_io::finish_mol2_record(record, params.sanitize, params.remove_hydrogens)
                .map_err(MolecularIoError::Mol2Post)?;
        let state = record.take_post_state();
        Self::from_parsed_parts_with_derived_state(
            record.topology,
            record.coordinates,
            record.properties,
            state.valence,
            state.rings,
        )
        .map_err(MolecularIoError::Construction)
    }
    pub fn read_mol2(path: &str) -> Result<Self, MolecularIoError> {
        Self::read_mol2_with_params(path, &Mol2ReadParams::default())
    }
    pub fn read_mol2_with_params(
        path: &str,
        params: &Mol2ReadParams,
    ) -> Result<Self, MolecularIoError> {
        Self::from_mol2_with_params(&read_text(path)?, params)
    }
}

impl Molecule {
    pub(crate) fn mol_write_input(&self) -> cosmolkit_io::MolWriteInput<'_> {
        cosmolkit_io::MolWriteInput {
            topology: self.topology(),
            coordinates: self.coordinate_block_runtime(),
            properties: self.properties(),
            rings: self.derived_cache_runtime().valid_ring_info(),
        }
    }
    pub fn to_mol(&self) -> Result<String, MolecularIoError> {
        self.to_mol_with_params(&MolBlockWriteParams::default())
    }
    pub fn to_mol_with_params(
        &self,
        params: &MolBlockWriteParams,
    ) -> Result<String, MolecularIoError> {
        cosmolkit_io::write_mol_block_with_params(self.mol_write_input(), params)
            .map_err(MolecularIoError::MolWrite)
    }
    pub fn to_sdf(&self) -> Result<String, MolecularIoError> {
        self.to_sdf_with_params(&MolBlockWriteParams::default())
    }
    pub fn to_sdf_with_params(
        &self,
        params: &MolBlockWriteParams,
    ) -> Result<String, MolecularIoError> {
        cosmolkit_io::write_sdf_with_params(self.mol_write_input(), params)
            .map_err(MolecularIoError::MolWrite)
    }
    pub fn to_sdf_2d(&self) -> Result<String, MolecularIoError> {
        self.to_sdf_2d_with_params(&MolBlockWriteParams::default())
    }
    pub fn to_sdf_2d_with_params(
        &self,
        params: &MolBlockWriteParams,
    ) -> Result<String, MolecularIoError> {
        cosmolkit_io::write_sdf_2d_with_params(self.mol_write_input(), params)
            .map_err(MolecularIoError::MolWrite)
    }
    pub fn to_sdf_3d(&self) -> Result<String, MolecularIoError> {
        self.to_sdf_3d_with_params(&MolBlockWriteParams::default())
    }
    pub fn to_sdf_3d_with_params(
        &self,
        params: &MolBlockWriteParams,
    ) -> Result<String, MolecularIoError> {
        cosmolkit_io::write_sdf_3d_with_params(self.mol_write_input(), params)
            .map_err(MolecularIoError::MolWrite)
    }
    pub fn write_mol(&self, path: &str) -> Result<(), MolecularIoError> {
        self.write_mol_with_params(path, &MolBlockWriteParams::default())
    }
    pub fn write_mol_with_params(
        &self,
        path: &str,
        params: &MolBlockWriteParams,
    ) -> Result<(), MolecularIoError> {
        write_text(path, &self.to_mol_with_params(params)?)
    }
    pub fn write_sdf(&self, path: &str) -> Result<(), MolecularIoError> {
        self.write_sdf_with_params(path, &MolBlockWriteParams::default())
    }
    pub fn write_sdf_with_params(
        &self,
        path: &str,
        params: &MolBlockWriteParams,
    ) -> Result<(), MolecularIoError> {
        write_text(path, &self.to_sdf_with_params(params)?)
    }
    pub fn write_sdf_files(
        &self,
        directory: &str,
        file_name: Option<&str>,
    ) -> Result<PathBuf, MolecularIoError> {
        self.write_sdf_files_with_params(directory, file_name, &MolBlockWriteParams::default())
    }
    pub fn write_sdf_files_with_params(
        &self,
        directory: &str,
        file_name: Option<&str>,
        params: &MolBlockWriteParams,
    ) -> Result<PathBuf, MolecularIoError> {
        let dir = expanded_path(directory)?;
        if !dir.exists() {
            return Err(MolecularIoError::Io {
                path: dir,
                source: std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    "directory does not exist",
                ),
            });
        }
        if !dir.is_dir() {
            return Err(MolecularIoError::Io {
                path: dir,
                source: std::io::Error::new(
                    std::io::ErrorKind::NotADirectory,
                    "path is not a directory",
                ),
            });
        }
        let name = file_name.unwrap_or("molecule.sdf");
        if name.trim().is_empty() {
            return Err(MolecularIoError::Parameter {
                name: "file_name",
                detail: "file_name cannot be empty",
            });
        }
        let path = dir.join(name);
        let text = self.to_sdf_with_params(params)?;
        std::fs::write(&path, text).map_err(|source| MolecularIoError::Io {
            path: path.clone(),
            source,
        })?;
        Ok(path)
    }
}
