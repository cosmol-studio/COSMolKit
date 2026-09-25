//! File-format IO over detached `cosmolkit-model` values.

mod bio_pdb;
#[doc(hidden)]
pub mod cif;
pub mod mol2;
mod mol_post;
mod numeric;
pub mod pdb;
mod pdb_chemistry;
pub mod sdf;
mod sdf_sgroups;
pub mod xyz;

pub use mol_post::{MolPostError, MolPostParams, finish_mol_block_record};
pub use mol2::{
    Mol2ReadError, Mol2ReadParams, Mol2Record, Mol2Type, read_mol2_detached,
    read_mol2_detached_with_params,
};
pub use pdb::{
    PdbReadError, PdbReadParams, PdbWriteError, PdbWriteParams, read_pdb_detached,
    read_pdb_detached_with_params, write_pdb_detached, write_pdb_detached_with_params,
};
pub use pdb_chemistry::{
    PdbPostprocessError, PdbPostprocessParams, apply_standard_pdb_residue_chirality_detached,
    postprocess_pdb_detached,
};

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
pub use xyz::{
    XyzReadError, XyzWriteError, XyzWriteParams, read_xyz_detached, write_xyz_detached,
    write_xyz_detached_with_params,
};
