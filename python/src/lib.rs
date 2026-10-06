//! Canonical Python projections of the public COSMolKit facade.
//!
//! Default and explicit drawing selections share one live Molecule class.
//! Historical adapters remain as source evidence, outside module compilation.

mod canonical_atom_bond;
mod canonical_batch;
mod canonical_batch_fingerprint_values;
mod canonical_batch_params;
mod canonical_binary;
mod canonical_bio_binding;
mod canonical_bio_metadata;
mod canonical_bio_residue;
mod canonical_builder;
mod canonical_chemistry_values;
mod canonical_coordinate_input;
mod canonical_descriptor_binding;
mod canonical_detached_blocks;
mod canonical_element_metadata;
mod canonical_error_accessors;
mod canonical_error_values;
mod canonical_fingerprint_values;
mod canonical_group_values;
mod canonical_molecular_hash;
mod canonical_operation_metadata;
mod canonical_potential_stereo;
mod canonical_sdf;
mod canonical_search;
mod canonical_smiles_writer;
mod canonical_stereo_queries;
mod canonical_valence;
mod canonical_values;
mod drawing_binding;
mod user_path;

#[cfg(feature = "stubgen")]
pyo3_stub_gen::define_stub_info_gatherer!(stub_info);

mod mmff_binding;

mod tautomer_binding;

mod canonical_property_values;
mod uff_binding;

mod alignment_binding;
mod canonical_layered;
mod canonical_maccs;
mod canonical_path_score;

mod canonical_pattern;
mod conformer_binding;

mod canonical_topological;
