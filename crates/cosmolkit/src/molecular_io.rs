//! Thin molecular IO projections. Chemistry and counted output belong to IO.
use crate::{Molecule, OperationError, SdfError, SdfReadParams, SdfRecord};
pub use cosmolkit_io::{
    Mol2PostError, Mol2ReadError, Mol2Type, MolBlockWriteParams, MolCoordinateSelection,
    MolWriteError, SdfFormat, XyzReadError, XyzWriteError, XyzWriteParams,
};
use std::{
    fmt,
    fs::File,
    io::{Read, Write},
    path::PathBuf,
};

#[derive(Debug)]
pub enum MolecularIoError {
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    XyzRead(XyzReadError),
    XyzWrite(XyzWriteError),
    Construction(OperationError),
    Mol2Read(Mol2ReadError),
    Mol2Post(Mol2PostError),
    MolWrite(MolWriteError),
    Sdf(SdfError),
    NoRecord {
        format: &'static str,
    },
    OutputUtf8 {
        format: &'static str,
        source: std::string::FromUtf8Error,
    },
    Parameter {
        name: &'static str,
        detail: &'static str,
    },
}
impl fmt::Display for MolecularIoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::XyzRead(e) => write!(f, "XYZ read: {e}"),
            Self::XyzWrite(e) => write!(f, "XYZ write: {e}"),
            Self::Construction(e) => e.fmt(f),
            Self::Mol2Read(e) => write!(f, "MOL2 read: {e}"),
            Self::Mol2Post(e) => write!(f, "MOL2 finalization: {e}"),
            Self::MolWrite(e) => e.fmt(f),
            Self::Sdf(e) => e.fmt(f),
            Self::NoRecord { format } => write!(f, "{format} found no molecule record"),
            Self::OutputUtf8 { format, source } => {
                write!(f, "{format} output is not UTF-8: {source}")
            }
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
            Self::Construction(e) => Some(e),
            Self::Mol2Read(e) => Some(e),
            Self::Mol2Post(e) => Some(e),
            Self::MolWrite(e) => Some(e),
            Self::Sdf(e) => Some(e),
            Self::OutputUtf8 { source, .. } => Some(source),
            Self::Parameter { .. } | Self::NoRecord { .. } => None,
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
impl Molecule {
    /// Render an atom property using the source lexical conversion.
    pub fn atom_property_string(
        &self,
        id: crate::AtomId,
        key: &str,
    ) -> Result<Option<crate::PropertyText>, crate::PropertyStringError> {
        self.atom(id)
            .and_then(|atom| atom.prop(key))
            .map(crate::property_value_to_text)
            .transpose()
    }
    /// Render a bond property using the source lexical conversion.
    pub fn bond_property_string(
        &self,
        id: crate::BondId,
        key: &str,
    ) -> Result<Option<crate::PropertyText>, crate::PropertyStringError> {
        self.bond(id)
            .and_then(|bond| bond.prop(key))
            .map(crate::property_value_to_text)
            .transpose()
    }
    /// Parse atom identities and Cartesian coordinates; XYZ does not infer bonds.
    pub fn from_xyz_block(text: &str) -> Result<Self, MolecularIoError> {
        let (topology, coordinates, properties) =
            cosmolkit_io::read_xyz_detached(text).map_err(MolecularIoError::XyzRead)?;
        Self::from_parts(topology, coordinates, properties).map_err(MolecularIoError::Construction)
    }
    pub fn read_xyz(path: &str) -> Result<Self, MolecularIoError> {
        let path = expanded_path(path)?;
        let mut file = File::open(&path).map_err(|source| MolecularIoError::Io {
            path: path.clone(),
            source,
        })?;
        let mut text = String::new();
        file.read_to_string(&mut text)
            .map_err(|source| MolecularIoError::Io { path, source })?;
        Self::from_xyz_block(&text)
    }
    pub fn to_xyz(&self) -> Result<String, MolecularIoError> {
        self.to_xyz_with_params(&XyzWriteParams::default())
    }
    pub fn to_xyz_with_params(&self, params: &XyzWriteParams) -> Result<String, MolecularIoError> {
        let text = cosmolkit_io::write_xyz_detached_with_params(
            self.topology(),
            self.coordinate_block_runtime(),
            self.properties(),
            *params,
        )
        .map_err(MolecularIoError::XyzWrite)?;
        String::from_utf8(text.into_bytes()).map_err(|source| MolecularIoError::OutputUtf8 {
            format: "XYZ",
            source,
        })
    }
    pub fn write_xyz(&self, path: &str) -> Result<(), MolecularIoError> {
        self.write_xyz_with_params(path, &XyzWriteParams::default())
    }
    pub fn write_xyz_with_params(
        &self,
        path: &str,
        params: &XyzWriteParams,
    ) -> Result<(), MolecularIoError> {
        // RDKit❗✔️: void MolToXYZFile(const ROMol &mol, const std::string &fName, int confId,
        // RDKit❗✔️:                   unsigned int precision) {
        // RDKit❗✔️:   std::ofstream outStream(fName);
        // RDKit❗✔️:   if (!outStream) {
        // RDKit❗✔️:     std::ostringstream errout;
        // RDKit❗✔️:     errout << "Bad output file " << fName;
        // RDKit❗✔️:     throw BadFileException(errout.str());
        // RDKit❗✔️:   }
        // RDKit❗✔️:   outStream << MolToXYZBlock(mol, confId, precision);
        // RDKit❗✔️: }
        // Open/truncate before rendering, preserving the source failure order.
        // Raw counted bytes are written directly; the String projection above
        // is the only UTF-8 conversion. O(N + L), no extra output clone.
        let path = expanded_path(path)?;
        let mut file = File::create(&path).map_err(|source| MolecularIoError::Io {
            path: path.clone(),
            source,
        })?;
        let text = cosmolkit_io::write_xyz_detached_with_params(
            self.topology(),
            self.coordinate_block_runtime(),
            self.properties(),
            *params,
        )
        .map_err(MolecularIoError::XyzWrite)?;
        file.write_all(text.as_bytes())
            .map_err(|source| MolecularIoError::Io { path, source })
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

pub(crate) fn open_file(path: &str) -> Result<File, MolecularIoError> {
    let path = expanded_path(path)?;
    File::open(&path).map_err(|source| MolecularIoError::Io { path, source })
}
fn read_text(path: &str) -> Result<String, MolecularIoError> {
    let path = expanded_path(path)?;
    let mut file = File::open(&path).map_err(|source| MolecularIoError::Io {
        path: path.clone(),
        source,
    })?;
    let mut text = String::new();
    file.read_to_string(&mut text)
        .map_err(|source| MolecularIoError::Io { path, source })?;
    Ok(text)
}
pub(crate) fn output_string(
    text: cosmolkit_model::PropertyText,
    format: &'static str,
) -> Result<String, MolecularIoError> {
    String::from_utf8(text.into_bytes())
        .map_err(|source| MolecularIoError::OutputUtf8 { format, source })
}
impl Molecule {
    pub fn from_mol(text: &str) -> Result<Self, SdfError> {
        Self::from_mol_with_params(text, &SdfReadParams::default())
    }
    pub fn from_mol_with_params(text: &str, params: &SdfReadParams) -> Result<Self, SdfError> {
        let parsed = cosmolkit_io::read_mol_graph_record_detached_with_params(
            text,
            crate::sdf_supplier::data_params(params),
        )?;
        SdfRecord::from_parsed(parsed, 0)?
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
            allow_coordinate_generation: cfg!(feature = "cap-depict"),
        }
    }
    pub fn to_mol(&self) -> Result<String, MolecularIoError> {
        self.to_mol_with_params(&MolBlockWriteParams::default())
    }
    pub fn to_mol_with_params(
        &self,
        params: &MolBlockWriteParams,
    ) -> Result<String, MolecularIoError> {
        output_string(
            cosmolkit_io::write_mol_block_with_params(self.mol_write_input(), params)
                .map_err(MolecularIoError::MolWrite)?,
            "MOL",
        )
    }
    pub fn to_sdf(&self) -> Result<String, MolecularIoError> {
        self.to_sdf_with_params(&MolBlockWriteParams::default())
    }
    pub fn to_sdf_with_params(
        &self,
        params: &MolBlockWriteParams,
    ) -> Result<String, MolecularIoError> {
        output_string(
            cosmolkit_io::write_sdf_with_params(self.mol_write_input(), params)
                .map_err(MolecularIoError::MolWrite)?,
            "SDF",
        )
    }
    pub fn to_sdf_2d(&self) -> Result<String, MolecularIoError> {
        self.to_sdf_2d_with_params(&MolBlockWriteParams::default())
    }
    pub fn to_sdf_2d_with_params(
        &self,
        params: &MolBlockWriteParams,
    ) -> Result<String, MolecularIoError> {
        output_string(
            cosmolkit_io::write_sdf_2d_with_params(self.mol_write_input(), params)
                .map_err(MolecularIoError::MolWrite)?,
            "SDF",
        )
    }
    pub fn to_sdf_3d(&self) -> Result<String, MolecularIoError> {
        self.to_sdf_3d_with_params(&MolBlockWriteParams::default())
    }
    pub fn to_sdf_3d_with_params(
        &self,
        params: &MolBlockWriteParams,
    ) -> Result<String, MolecularIoError> {
        output_string(
            cosmolkit_io::write_sdf_3d_with_params(self.mol_write_input(), params)
                .map_err(MolecularIoError::MolWrite)?,
            "SDF",
        )
    }
    pub fn write_mol(&self, path: &str) -> Result<(), MolecularIoError> {
        self.write_mol_with_params(path, &MolBlockWriteParams::default())
    }
    pub fn write_mol_with_params(
        &self,
        path: &str,
        params: &MolBlockWriteParams,
    ) -> Result<(), MolecularIoError> {
        // RDKit❗✔️: void MolToMolFile(const ROMol &mol, const std::string &fName,
        // RDKit❗✔️:                   const MolWriterParams &params, int confId) {
        // RDKit❗✔️:   auto *outStream = new std::ofstream(fName.c_str());
        // RDKit❗✔️:   if (!(*outStream) || outStream->bad()) {
        // RDKit❗✔️:     delete outStream;
        // RDKit❗✔️:     std::ostringstream errout;
        // RDKit❗✔️:     errout << "Bad output file " << fName;
        // RDKit❗✔️:     throw BadFileException(errout.str());
        // RDKit❗✔️:   }
        // RDKit❗✔️:   std::string outString = MolToMolBlock(mol, params, confId);
        // RDKit❗✔️:   *outStream << outString;
        // RDKit❗✔️:   delete outStream;
        // RDKit❗✔️: }
        // File creation precedes source formatting; write counted bytes without decoding.
        // O(N+L) owner formatting plus one output write; no output clone.
        let path = expanded_path(path)?;
        let mut file = File::create(&path).map_err(|source| MolecularIoError::Io {
            path: path.clone(),
            source,
        })?;
        let text = cosmolkit_io::write_mol_block_with_params(self.mol_write_input(), params)
            .map_err(MolecularIoError::MolWrite)?;
        file.write_all(text.as_bytes())
            .map_err(|source| MolecularIoError::Io { path, source })
    }

    pub fn write_sdf(&self, path: &str) -> Result<(), MolecularIoError> {
        self.write_sdf_with_params(path, &MolBlockWriteParams::default())
    }
    pub fn write_sdf_with_params(
        &self,
        path: &str,
        params: &MolBlockWriteParams,
    ) -> Result<(), MolecularIoError> {
        // COSMolKit❗✔️:     fn write_sdf(
        // COSMolKit❗✔️:         &self,
        // COSMolKit❗✔️:         path: &str,
        // COSMolKit❗✔️:         format: Option<&str>,
        // COSMolKit❗✔️:         include_stereo: bool,
        // COSMolKit❗✔️:         kekulize: bool,
        // COSMolKit❗✔️:     ) -> PyResult<()> {
        // COSMolKit❗✔️:         let expanded_path = expand_user_path(path)?;
        // COSMolKit❗✔️:         let fmt = parse_sdf_format(format)?;
        // COSMolKit❗✔️:         let block = molecule_to_sdf_record_string(&self.inner, fmt, include_stereo, kekulize)
        // COSMolKit❗✔️:             .map_err(|err| PyValueError::new_err(format!("write_sdf failed: {err}")))?;
        // COSMolKit❗✔️:         let mut f = File::create(&expanded_path)
        // COSMolKit❗✔️:             .map_err(|e| PyValueError::new_err(format!("write_sdf create failed: {e}")))?;
        // COSMolKit❗✔️:         f.write_all(block.as_bytes())
        // COSMolKit❗✔️:             .map_err(|e| PyValueError::new_err(format!("write_sdf write failed: {e}")))?;
        // COSMolKit❗✔️:         Ok(())
        // COSMolKit❗✔️:     }
        // Preserve the frozen public SDF adapter's render-before-create order.
        // The detached owner returns counted bytes; the file path does not decode them.
        self.write_sdf_path(expanded_path(path)?, params)
    }
    fn write_sdf_path(
        &self,
        path: PathBuf,
        params: &MolBlockWriteParams,
    ) -> Result<(), MolecularIoError> {
        // COSMolKit❗✔️: let block = molecule_to_sdf_record_string(&self.inner, fmt, include_stereo, kekulize)
        // COSMolKit❗✔️:     .map_err(|err| PyValueError::new_err(format!("write_sdf failed: {err}")))?;
        // COSMolKit❗✔️: let mut f = File::create(&expanded_path)
        // COSMolKit❗✔️:     .map_err(|e| PyValueError::new_err(format!("write_sdf create failed: {e}")))?;
        // COSMolKit❗✔️: f.write_all(block.as_bytes())
        // COSMolKit❗✔️:     .map_err(|e| PyValueError::new_err(format!("write_sdf write failed: {e}")))?;
        // Complexity: O(N+L), one owner output buffer, one file write, no output copy.
        let text = cosmolkit_io::write_sdf_with_params(self.mol_write_input(), params)
            .map_err(MolecularIoError::MolWrite)?;
        let mut file = File::create(&path).map_err(|source| MolecularIoError::Io {
            path: path.clone(),
            source,
        })?;
        file.write_all(text.as_bytes())
            .map_err(|source| MolecularIoError::Io { path, source })
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
        // COSMolKit❗✔️:     fn write_sdf_to_directory(
        // COSMolKit❗✔️:         &self,
        // COSMolKit❗✔️:         directory: &str,
        // COSMolKit❗✔️:         file_name: Option<&str>,
        // COSMolKit❗✔️:         format: Option<&str>,
        // COSMolKit❗✔️:         include_stereo: bool,
        // COSMolKit❗✔️:         kekulize: bool,
        // COSMolKit❗✔️:     ) -> PyResult<String> {
        // COSMolKit❗✔️:         let expanded_directory = expand_user_path(directory)?;
        // COSMolKit❗✔️:         let dir = expanded_directory.as_path();
        // COSMolKit❗✔️:         if !dir.exists() {
        // COSMolKit❗✔️:             return Err(PyValueError::new_err(format!(
        // COSMolKit❗✔️:                 "directory does not exist: {directory}"
        // COSMolKit❗✔️:             )));
        // COSMolKit❗✔️:         }
        // COSMolKit❗✔️:         if !dir.is_dir() {
        // COSMolKit❗✔️:             return Err(PyValueError::new_err(format!(
        // COSMolKit❗✔️:                 "path is not a directory: {directory}"
        // COSMolKit❗✔️:             )));
        // COSMolKit❗✔️:         }
        // COSMolKit❗✔️:         let name = file_name.unwrap_or("molecule.sdf");
        // COSMolKit❗✔️:         if name.trim().is_empty() {
        // COSMolKit❗✔️:             return Err(PyValueError::new_err("file_name cannot be empty"));
        // COSMolKit❗✔️:         }
        // COSMolKit❗✔️:         let output = dir.join(name);
        // COSMolKit❗✔️:         let output_str = output
        // COSMolKit❗✔️:             .to_str()
        // COSMolKit❗✔️:             .ok_or_else(|| PyValueError::new_err("output path is not valid UTF-8"))?;
        // COSMolKit❗✔️:         self.write_sdf(output_str, format, include_stereo, kekulize)?;
        // COSMolKit❗✔️:         Ok(output_str.to_string())
        // COSMolKit❗✔️:     }
        // The Rust projection retains the native PathBuf. Python checks the
        // returned path at its String boundary; raw SDF bytes never pass through it.
        // O(path length + N + L); constant directory metadata checks and one join.
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
        self.write_sdf_path(path.clone(), params)?;
        Ok(path)
    }
}
