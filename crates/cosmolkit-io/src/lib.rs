//! File-format IO over detached `cosmolkit-model` values.
//!
//! Internal `molecule` enables molecular formats and their chemistry/query
//! dependencies. `bio` independently enables structural-biology formats and
//! CID selection; it does not activate molecular IO or the search owner.

#[cfg(all(feature = "bio", feature = "molecule"))]
mod bio_molecule;
#[cfg(all(feature = "bio", feature = "molecule"))]
pub use bio_molecule::{
    BioMoleculeConversionError, BioMoleculeParams, bio_structure_to_molecule_parts,
};
#[cfg(feature = "bio")]
mod bio_chemcomp;
#[cfg(feature = "bio")]
mod bio_cid;
#[cfg(feature = "bio")]
mod bio_mmcif;
#[cfg(feature = "bio")]
mod bio_numeric;
#[cfg(feature = "bio")]
mod bio_pdb;
#[cfg(feature = "bio")]
mod bio_read;
#[cfg(feature = "bio")]
mod bio_write;
#[cfg(feature = "bio")]
#[doc(hidden)]
pub mod cif;
#[cfg(feature = "molecule")]
pub mod mol2;
#[cfg(feature = "molecule")]
mod mol_post;
#[cfg(feature = "molecule")]
mod numeric;
#[cfg(feature = "molecule")]
pub mod pdb;
#[cfg(feature = "molecule")]
mod pdb_chemistry;
#[cfg(feature = "molecule")]
pub mod sdf;
#[cfg(feature = "molecule")]
mod sdf_sgroups;
#[cfg(feature = "molecule")]
pub mod xyz;

#[cfg(feature = "bio")]
pub use bio_cid::{
    BioSelectionParseError, SelectionSeqidRangeError, SelectionSyntaxError, read_bio_selection,
    write_bio_selection,
};
#[cfg(feature = "bio")]
pub use bio_mmcif::{BioMmcifReadError, BioMmcifReadStage, read_mmcif_bio_structure};
#[cfg(feature = "bio")]
pub use bio_pdb::{BioPdbReadError, BioPdbReadParams, BioPdbReadStage, read_pdb_bio_structure};
#[cfg(feature = "bio")]
pub use bio_read::{BioReadError, BioReadParams, read_bio_structure, read_bio_structure_file};
#[cfg(feature = "bio")]
pub use bio_write::{
    BioMmcifWriteError, BioMmcifWriteParams, bio_structure_to_mmcif_text,
    write_bio_structure_mmcif_file,
};
#[cfg(feature = "molecule")]
pub use mol_post::{MolPostError, MolPostParams, finish_mol_block_record};
#[cfg(feature = "molecule")]
pub use mol2::{
    Mol2ReadError, Mol2ReadParams, Mol2Record, Mol2Type, read_mol2_detached,
    read_mol2_detached_with_params,
};
#[cfg(feature = "molecule")]
pub use pdb::{
    PdbReadError, PdbReadParams, PdbWriteError, PdbWriteParams, read_pdb_detached,
    read_pdb_detached_with_params, write_pdb_detached, write_pdb_detached_with_params,
};
#[cfg(feature = "molecule")]
pub use pdb_chemistry::{
    PdbPostprocessError, PdbPostprocessParams, apply_standard_pdb_residue_chirality_detached,
    postprocess_pdb_detached,
};

#[cfg(feature = "molecule")]
pub use sdf::{
    MolBlockReadParams, MolBlockRecord, QueryMolBlockRecord, SdfCoordinateMode, SdfDataReadParams,
    SdfGraphDataset, SdfGraphReader, SdfGraphRecord, SdfReadError, SdfRecord, SdfRecordMetadata,
    SdfRecordText, SdfWriteError, index_sdf_records, read_mol_block_detached,
    read_mol_block_detached_with_params, read_sdf_graph_record_detached,
    read_sdf_graph_record_detached_with_params, read_sdf_record_detached,
    read_sdf_record_detached_with_params, read_sdf_record_text, read_sdf_records_detached,
    read_sdf_records_detached_with_params, read_v2000_detached, read_v2000_detached_with_params,
    read_v3000_detached, read_v3000_detached_with_params, write_sdf_record_detached,
    write_v2000_detached, write_v3000_detached,
};
#[cfg(feature = "molecule")]
pub use xyz::{
    XyzReadError, XyzWriteError, XyzWriteParams, read_xyz_detached, write_xyz_detached,
    write_xyz_detached_with_params,
};
#[cfg(feature = "bio")]
mod bio_pdb_write;
#[cfg(feature = "bio")]
pub use bio_pdb_write::{
    BioPdbWriteError, BioPdbWriteParams, bio_structure_to_pdb_text, write_bio_structure_pdb_file,
};
