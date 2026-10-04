//! Public structural objects and their thin algorithm adapters.
//! See dev/bio_architecture.md for ownership and lightweight operation semantics.
mod runtime;
use crate::{
    BioMmcifReadError, BioPdbReadError, BioReadError, BioStructureError, ProteinProjectionError,
};
pub use runtime::{BioStructure, Protein};

/// Thin root CID-selection wrapper (BIO-CID C26) around the one
/// canonical detached selection data. The field stays private: this is a
/// projection type, not a second selection engine; construction goes
/// through the frozen CID API (from_cid, later units) and nothing mutable
/// or parser-shaped is re-exported here.
///
/// Use `structure.with_selection(&selection)` for a new selected value,
/// or `structure.retain_selection_(&selection)` to replace it in place.
/// Both also apply to `Protein`. Selected empty models/chains/residues and
/// source metadata (including assembly references) are retained as in Gemmi;
/// fields omitted by Gemmi's `Structure::empty_copy` reset to their defaults.
pub struct BioSelection {
    pub(crate) data: cosmolkit_bio::BioSelectionData,
}

impl BioSelection {
    /// Registered Experimental constructor (BIO-CID C27): parse real CID
    /// text through the one canonical IO reader (Gemmi `Selection::
    /// Selection(const std::string&)` -> parse_cid, select.cpp:261-263);
    /// failures surface the typed parse-error vocabulary unchanged. A
    /// thin delegate — no parser logic at the root.
    pub fn from_cid(cid: &str) -> Result<Self, crate::BioSelectionParseError> {
        Ok(Self {
            data: cosmolkit_io::read_bio_selection(cid)?,
        })
    }

    /// Registered Experimental serializer (BIO-CID C28): re-serialize
    /// through the canonical IO owner (Gemmi `Selection::str`,
    /// select.cpp:288-327). A thin delegate — no formatter logic at the
    /// root.
    pub fn to_cid(&self) -> String {
        cosmolkit_io::write_bio_selection(&self.data)
    }
}

/// Reading and explicit amino-acid projection failures retain their source.
#[derive(Debug)]
pub enum ProteinReadError {
    Structure(BioReadError),
    Pdb(BioPdbReadError),
    Mmcif(BioMmcifReadError),
    Projection(ProteinProjectionError),
}
impl std::fmt::Display for ProteinReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Structure(e) => write!(f, "{e}"),
            Self::Pdb(e) => write!(f, "{e}"),
            Self::Mmcif(e) => write!(f, "{e}"),
            Self::Projection(e) => write!(f, "{e}"),
        }
    }
}
impl std::error::Error for ProteinReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Structure(e) => e,
            Self::Pdb(e) => e,
            Self::Mmcif(e) => e,
            Self::Projection(e) => e,
        })
    }
}

/// Final data validation failure; a failed operation never replaces its input.
#[derive(Debug)]
pub enum BioOperationError {
    Structure(BioStructureError),
    Protein(ProteinProjectionError),
    Selection(crate::BioSelectionCopyError),
}
impl std::fmt::Display for BioOperationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Structure(e) => write!(f, "{e}"),
            Self::Protein(e) => write!(f, "{e}"),
            Self::Selection(e) => write!(f, "{e}"),
        }
    }
}
impl std::error::Error for BioOperationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Structure(e) => e,
            Self::Protein(e) => e,
            Self::Selection(e) => e,
        })
    }
}

fn selection_impl(
    access: runtime::SelectionAccess<'_>,
    selection: &BioSelection,
) -> Result<runtime::SelectionAccessReplacement, BioOperationError> {
    // Thin domain delegation: only declared borrowed fields cross the boundary.
    // The generated runtime alone installs and validates the six returned blocks.
    let source = cosmolkit_bio::BioStructureCopySource {
        input_format: *access.input_format,
        models: access.models,
        chains: access.chains,
        residues: access.residues,
        atoms: access.atoms,
        entities: access.entities,
        connections: access.connections,
        cispeps: access.cispeps,
        mod_residues: access.mod_residues,
        helices: access.helices,
        sheets: access.sheets,
        metadata: access.metadata,
        source_state: access.source_state,
        coordinates: access.coordinates,
        crystal: access.crystal,
        ncs_operators: access.ncs_operators,
        assemblies: access.assemblies,
    };
    let blocks = cosmolkit_bio::copy_selection_blocks(source, &selection.data)
        .map_err(BioOperationError::Selection)?;
    Ok(runtime::SelectionAccessReplacement {
        models: blocks.models,
        chains: blocks.chains,
        residues: blocks.residues,
        atoms: blocks.atoms,
        source_state: blocks.source_state,
        coordinates: blocks.coordinates,
    })
}

fn translate_coordinates_impl(
    access: runtime::TranslateCoordinatesAccess<'_>,
    offset: [f64; 3],
) -> Result<(), BioOperationError> {
    cosmolkit_bio::translate_coordinates(access.coordinates, offset);
    Ok(())
}

// Compile probes live beside real operation bodies, outside private storage.
#[cfg(cosmolkit_bio_privacy_probe)]
#[allow(dead_code, unused_variables)]
fn bio_privacy_probe(
    access: runtime::TranslateCoordinatesAccess<'_>,
    structure: BioStructure,
    protein: Protein,
) {
    let _ = access.coordinates.positions_mut();
    #[cfg(cosmolkit_bio_privacy_case = "undeclared")]
    let _ = access.atoms;
    #[cfg(cosmolkit_bio_privacy_case = "storage")]
    {
        let _ = structure.data;
        let _ = protein.structure;
    }
}

#[cfg(cosmolkit_bio_privacy_probe)]
#[allow(dead_code, unused_variables)]
fn bio_replacement_privacy_probe(access: runtime::SelectionAccess<'_>) {
    let replacement = runtime::SelectionAccessReplacement {
        models: std::sync::Arc::clone(access.models),
        chains: std::sync::Arc::clone(access.chains),
        residues: std::sync::Arc::clone(access.residues),
        atoms: std::sync::Arc::clone(access.atoms),
        coordinates: std::sync::Arc::clone(access.coordinates),
        source_state: std::sync::Arc::clone(access.source_state),
    };
    #[cfg(cosmolkit_bio_privacy_case = "replacement_readonly")]
    {
        *access.atoms = std::sync::Arc::new(Vec::new());
    }
    #[cfg(cosmolkit_bio_privacy_case = "metadata_readonly")]
    {
        let _ = std::sync::Arc::make_mut(access.metadata);
    }
    #[cfg(cosmolkit_bio_privacy_case = "replacement_undeclared")]
    let _ = replacement.assemblies;
    #[cfg(cosmolkit_bio_privacy_case = "replacement_required")]
    let _ = runtime::SelectionAccessReplacement {
        models: replacement.models,
        chains: replacement.chains,
        residues: replacement.residues,
        atoms: replacement.atoms,
        coordinates: replacement.coordinates,
    };
}
