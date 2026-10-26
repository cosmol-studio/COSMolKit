//! Binding protocol probe; chemistry remains in the public Rust facade.
use cosmolkit as ck;
use serde::{Deserialize, Serialize};
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;

#[cfg(feature = "registry-check")]
mod registry_snapshot;

#[wasm_bindgen]
pub fn version() -> String {
    ck::version().into()
}

#[cfg(feature = "cap-alignment")]
mod alignment_calls;
#[cfg(feature = "cap-alignment")]
mod alignment_errors;
#[cfg(feature = "cap-alignment")]
mod alignment_values;
mod aromaticity_errors;
mod assignment;
#[cfg(feature = "cap-fingerprints")]
mod atom_code;
mod atom_properties;
mod atoms;
mod atropisomer_errors;
mod attachment_values;
#[cfg(feature = "cap-fingerprints")]
mod avalon;
#[cfg(feature = "cap-batch")]
mod batch;
#[cfg(all(feature = "cap-batch", feature = "cap-fingerprints"))]
mod batch_output;
#[cfg(feature = "cap-batch")]
mod batch_queries;
#[cfg(feature = "cap-bio")]
mod bio_assemblies;
#[cfg(feature = "cap-bio")]
mod bio_conversion_errors;
#[cfg(feature = "cap-bio")]
mod bio_entities;
#[cfg(feature = "cap-bio")]
mod bio_geometry;
#[cfg(feature = "cap-bio")]
mod bio_hierarchy;
#[cfg(feature = "cap-bio")]
mod bio_lookup;
#[cfg(feature = "cap-bio")]
mod bio_metadata;
#[cfg(feature = "cap-bio")]
mod bio_molecule;
#[cfg(feature = "cap-bio")]
mod bio_operations;
#[cfg(feature = "cap-bio")]
mod bio_parts;
#[cfg(feature = "cap-bio")]
mod bio_read_errors;
#[cfg(feature = "cap-bio")]
mod bio_relationships;
#[cfg(feature = "cap-bio")]
mod bio_rows;
#[cfg(feature = "cap-bio")]
mod bio_scalar_metadata;
#[cfg(feature = "cap-bio")]
mod bio_source_state;
#[cfg(feature = "cap-bio")]
mod bio_source_values;
#[cfg(feature = "cap-bio")]
mod bio_text;
#[cfg(feature = "cap-bio")]
mod bio_writers;
mod builder;
mod callback;
mod canonical_errors;
mod chemistry;
mod chemistry_problems;
mod cip_labeler_errors;
mod cip_rank_errors;
#[cfg(feature = "cap-conformer")]
mod conformer_calls;
#[cfg(feature = "cap-conformer")]
mod conformer_errors;
mod construction_values;
mod coordinate_input;
mod coordinates;
#[cfg(feature = "cap-depict")]
mod depict_errors;
#[cfg(feature = "cap-descriptors")]
mod descriptor_calls;
#[cfg(feature = "cap-descriptors")]
mod descriptor_errors;
mod detached_coordinates;
mod double_stereo_errors;
#[cfg(feature = "cap-depict")]
mod drawing;
mod enums;
mod errors;
mod external_errors;
mod fields;
#[cfg(feature = "cap-fingerprints")]
mod fingerprint_batch;
#[cfg(feature = "cap-fingerprints")]
mod fingerprint_errors;
#[cfg(feature = "cap-fingerprints")]
mod fingerprint_output;
#[cfg(feature = "cap-fingerprints")]
mod fingerprint_values;
#[cfg(feature = "cap-forcefields")]
mod forcefield_calls;
#[cfg(feature = "cap-forcefields")]
mod forcefield_execution_errors;
#[cfg(feature = "cap-forcefields")]
mod forcefield_kernel_errors;
#[cfg(feature = "cap-forcefields")]
mod forcefield_params;
mod fragment_errors;
mod fragments;
mod group_values;
mod hydrogen_errors;
#[cfg(feature = "cap-inchi")]
mod inchi;
mod input;
mod io_configuration;
mod io_write_errors;
mod iterator;
mod kekulize_errors;
mod legacy_stereo_errors;
#[cfg(feature = "cap-fingerprints")]
mod maccs;
#[cfg(feature = "cap-search")]
mod match_errors;
mod matrix_calls;
mod matrix_errors;
#[cfg(feature = "cap-search")]
mod mcs;
#[cfg(feature = "cap-search")]
mod mcs_errors;
mod metadata_queries;
mod metadata_values;
#[cfg(feature = "cap-forcefields")]
mod mmff_errors;
mod model_errors;
mod mol2_errors;
mod mol_post_errors;
mod mol_post_support_errors;
mod mol_write_errors;
#[cfg(feature = "cap-hashing")]
mod molecular_hash;
mod molecular_io;
mod molecular_io_errors;
#[cfg(feature = "cap-fingerprints")]
mod morgan_calls;
#[cfg(feature = "cap-fingerprints")]
mod morgan_configuration;
#[cfg(feature = "cap-fingerprints")]
mod morgan_errors;
#[cfg(feature = "cap-fingerprints")]
mod morgan_reusable;
mod notation;
mod operation_errors;
mod operation_metadata;
#[cfg(feature = "cap-fingerprints")]
mod pair_torsion;
mod parameters;
mod parser_errors;
#[cfg(feature = "cap-smiles")]
mod path_errors;
#[cfg(feature = "cap-fingerprints")]
mod path_fingerprints;
#[cfg(feature = "cap-fingerprints")]
mod pattern;
#[cfg(feature = "cap-bio")]
mod pdb_read_errors;
mod periodic_errors;
#[cfg(feature = "cap-forcefields")]
mod persistent_forcefields;
mod potential_stereo_errors;
mod property_errors;
mod property_values;
#[cfg(feature = "cap-bio")]
mod protein;
#[cfg(feature = "cap-search")]
mod query_calls;
mod query_errors;
#[cfg(feature = "cap-search")]
mod query_matches;
mod query_model_errors;
mod query_values;
#[cfg(feature = "cap-reaction")]
mod reaction;
#[cfg(feature = "cap-reaction")]
mod reaction_calls;
#[cfg(feature = "cap-reaction")]
mod reaction_errors;
#[cfg(feature = "cap-reaction")]
mod reaction_execution_errors;
#[cfg(feature = "cap-reaction")]
mod reaction_values;
#[cfg(feature = "cap-bio")]
mod residue_values;
mod ring_errors;
mod sanitize;
mod sanitize_errors;
mod sanitize_leaf_errors;
mod sanitize_support_errors;
#[cfg(feature = "cap-hashing")]
mod scaffold_errors;
mod sdf_errors;
mod sdf_records;
#[cfg(feature = "cap-search")]
mod search_callbacks;
#[cfg(feature = "cap-search")]
mod smarts_errors;
#[cfg(feature = "cap-search")]
mod smarts_parse_errors;
mod smiles_errors;
mod smiles_parse_errors;
mod smiles_stereo_errors;
#[cfg(feature = "cap-fingerprints")]
mod sparse_fingerprints;
mod stereo_errors;
mod stereo_groups;
mod stereo_order_errors;
mod stereo_reads;
#[cfg(feature = "cap-stereoisomers")]
mod stereoisomers;
mod storage;
#[cfg(feature = "cap-tautomer")]
mod tautomer_callbacks;
#[cfg(feature = "cap-tautomer")]
mod tautomer_calls;
#[cfg(feature = "cap-tautomer")]
mod tautomer_catalog;
#[cfg(feature = "cap-tautomer")]
mod tautomer_errors;
#[cfg(feature = "cap-tautomer")]
mod tautomer_score;
mod topology_errors;
mod topology_values;
#[cfg(feature = "cap-fingerprints")]
mod torsion_path;
#[cfg(feature = "cap-fingerprints")]
mod torsion_reusable;
mod transform_errors;
#[cfg(feature = "cap-forcefields")]
mod uff_errors;
mod valence_errors;
#[cfg(any(feature = "cap-search", feature = "cap-io"))]
mod wedge_errors;
#[cfg(any(feature = "cap-alignment", feature = "cap-fingerprints"))]
mod worker_errors;
mod xyz_errors;

// The field list and field types come from the native declaration, not a JS list.
ck::__ck_configuration_RingSearchParams!(crate::project_parameters, ck::RingSearchParams, {}, {
    RingSearchParams
});
ck::__ck_configuration_SmilesParseParams!(crate::project_parameters, ck::SmilesParseParams, {}, {
    SmilesParseParams
});
ck::__ck_configuration_SmilesWriteParams!(crate::project_parameters, ck::SmilesWriteParams, {}, {
    SmilesWriteParams
});
ck::__ck_configuration_RandomSmilesWriteParams!(
    crate::project_parameters,
    ck::RandomSmilesWriteParams,
    {},
    { RandomSmilesWriteParams }
);
ck::__ck_configuration_FragmentSmilesWriteParams!(
    crate::project_parameters,
    ck::FragmentSmilesWriteParams,
    {},
    { FragmentSmilesWriteParams }
);

#[wasm_bindgen]
pub struct Molecule(ck::Molecule);
crate::fields::value_field!(ck::Molecule, Molecule, MoleculeInput, "Molecule");
crate::fields::optional_value_array!(
    ck::Molecule,
    Molecule,
    MoleculeArray,
    "ReadonlyArray<Molecule | null | undefined>"
);

#[wasm_bindgen(js_name = molFromSmiles)]
pub fn mol_from_smiles(
    smiles: Ts<input::TextInput>,
    params: Option<SmilesParseParamsInput>,
) -> Result<Molecule, JsValue> {
    Molecule::from_smiles(smiles, params)
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(constructor)]
    pub fn new(smiles: Ts<input::TextInput>) -> Result<Molecule, JsValue> {
        let smiles = input::read::<String>(smiles)?;
        ck::Molecule::from_smiles(&smiles)
            .map(Self)
            .map_err(smiles_errors::SmilesError::from_rust)
    }

    pub fn copy(&self) -> Molecule {
        Self(self.0.copy())
    }

    #[wasm_bindgen(js_name = fromSmiles)]
    pub fn from_smiles(
        smiles: Ts<input::TextInput>,
        params: Option<SmilesParseParamsInput>,
    ) -> Result<Molecule, JsValue> {
        let smiles = input::read::<String>(smiles)?;
        SmilesParseParams::with_input(params.as_ref(), |params| {
            ck::Molecule::from_smiles_with_params(&smiles, params)
                .map(Self)
                .map_err(smiles_errors::SmilesError::from_rust)
        })
    }

    #[wasm_bindgen(js_name = toSmiles)]
    pub fn to_smiles(&self, params: Option<SmilesWriteParamsInput>) -> Result<String, JsValue> {
        SmilesWriteParams::with_input(params.as_ref(), |params| {
            self.0
                .to_smiles_with_params(params)
                .map_err(smiles_errors::SmilesWriteError::from_rust)
                .and_then(|text| property_errors::text(&text))
        })
    }

    #[wasm_bindgen(js_name = toFragmentSmiles)]
    pub fn to_fragment_smiles(
        &self,
        params: FragmentSmilesWriteParamsInput,
    ) -> Result<String, JsValue> {
        FragmentSmilesWriteParams::with_input(Some(&params), |params| {
            self.0
                .to_fragment_smiles_with_params(params)
                .map_err(smiles_errors::SmilesWriteError::from_rust)
                .and_then(|text| property_errors::text(&text))
        })
    }

    #[wasm_bindgen(getter, js_name = numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.0.num_atoms()
    }

    #[wasm_bindgen(getter, js_name = numBonds)]
    pub fn num_bonds(&self) -> usize {
        self.0.num_bonds()
    }

    pub fn atom(
        &self,
        index: Ts<input::IndexInput>,
    ) -> Result<Option<atoms::MoleculeAtom>, JsValue> {
        let index = input::read::<usize>(index)?;
        Ok(self.0.atom(ck::AtomId::new(index)).map(atoms::MoleculeAtom))
    }

    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> Vec<atoms::MoleculeAtom> {
        self.0
            .atoms()
            .into_iter()
            .map(atoms::MoleculeAtom)
            .collect()
    }

    pub fn bond(&self, index: Ts<input::IndexInput>) -> Result<Option<atoms::Bond>, JsValue> {
        let index = input::read::<usize>(index)?;
        Ok(self
            .0
            .bond(ck::BondId::new(index))
            .cloned()
            .map(atoms::Bond))
    }

    #[wasm_bindgen(getter)]
    pub fn bonds(&self) -> Vec<atoms::Bond> {
        self.0.bonds().iter().cloned().map(atoms::Bond).collect()
    }

    #[wasm_bindgen(js_name = toRandomSmiles)]
    pub fn to_random_smiles(
        &self,
        count: Ts<input::U32Input>,
        seed: Ts<input::U32Input>,
        params: Option<RandomSmilesWriteParamsInput>,
    ) -> Result<Vec<String>, JsValue> {
        let count = input::read::<u32>(count)?;
        let seed = input::read::<u32>(seed)?;
        RandomSmilesWriteParams::with_input(params.as_ref(), |params| {
            self.0
                .to_random_smiles_with_params(count, seed, params)
                .map_err(smiles_errors::SmilesWriteError::from_rust)?
                .iter()
                .map(property_errors::text)
                .collect()
        })
    }

    #[wasm_bindgen(js_name = toCxSmiles)]
    pub fn to_cx_smiles(
        &self,
        params: Option<notation::CxSmilesWriteParamsInput>,
    ) -> Result<String, JsValue> {
        notation::CxSmilesWriteParams::with_input(params.as_ref(), |params| {
            self.0
                .to_cx_smiles_with_params(params)
                .map_err(smiles_errors::SmilesWriteError::from_rust)
                .and_then(|text| property_errors::text(&text))
        })
    }

    #[wasm_bindgen(js_name = toFragmentCxSmiles)]
    pub fn to_fragment_cx_smiles(
        &self,
        params: notation::FragmentCxSmilesWriteParamsInput,
    ) -> Result<String, JsValue> {
        notation::FragmentCxSmilesWriteParams::with_input(Some(&params), |params| {
            self.0
                .to_fragment_cx_smiles_with_params(params)
                .map_err(smiles_errors::SmilesWriteError::from_rust)
                .and_then(|text| property_errors::text(&text))
        })
    }

    #[cfg(feature = "cap-descriptors")]
    #[wasm_bindgen(js_name = numRings)]
    pub fn num_rings(&self) -> Result<u32, JsValue> {
        self.0
            .num_rings()
            .map_err(descriptor_errors::DescriptorReadError::from_rust)
    }

    /// Protocol probe: deliver an independently owned COW snapshot, not a borrow.
    pub fn inspect(&self, callback: &InspectCallback) -> Result<bool, JsValue> {
        let value = callback.call1(&JsValue::UNDEFINED, &JsValue::from(self.copy()))?;
        input::read::<bool>(Ts::new_unchecked(value))
    }
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(extends = js_sys::Function, typescript_type = "((molecule: Molecule) => boolean)")]
    pub type InspectCallback;

}

// These are protocol samples, not new chemistry API or duplicated native Params.
#[derive(Serialize, Deserialize, Tsify)]
#[serde(rename_all = "snake_case")]
pub enum Mode {
    Preserve,
    Replace,
}

#[derive(Serialize, Deserialize, Tsify)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Outcome {
    Completed { count: u32 },
    Failed { message: String },
}

#[derive(Serialize, Deserialize, Tsify)]
#[serde(deny_unknown_fields)]
pub struct ArrayInput {
    #[serde(deserialize_with = "input::exact_array")]
    pub point: [f64; 3],
    pub indices: Vec<u32>,
    pub mode: Mode,
}

#[wasm_bindgen]
pub fn roundtrip(input: Ts<ArrayInput>) -> Result<Ts<ArrayInput>, JsValue> {
    Ok(input::decode(input)?.into_ts().map_err(JsError::from)?)
}

#[wasm_bindgen]
pub fn outcome(input: Ts<Outcome>) -> Result<Ts<Outcome>, JsValue> {
    Ok(input::decode(input)?.into_ts().map_err(JsError::from)?)
}

/// Protocol fixture exercises the same native cause projection as real queries.
#[wasm_bindgen]
pub fn raise() -> Result<(), JsValue> {
    Err(operation_errors::OperationError::from_rust(
        ck::OperationError::Rings(ck::RingFindingError::InvalidTopology(
            ck::TopologyValidationError::TemplateAttachmentOrder {
                atom: ck::AtomId::new(1),
                source: ck::TemplateAttachmentOrderError::TargetRemoved {
                    position: 2,
                    target: ck::AtomId::new(3),
                },
            },
        )),
    ))
}

/// Owned typed-array output; mutation in JS cannot modify native storage.
#[wasm_bindgen]
pub fn bytes(input: Ts<input::BytesInput>) -> Result<Vec<u8>, JsValue> {
    input::read(input)
}

/// Explicit 64-bit integer boundary: wasm-bindgen generates bigint, not number.
#[wasm_bindgen]
pub fn integer(input: Ts<input::U64Input>) -> Result<u64, JsValue> {
    input::read(input)
}

// Protocol-only native records: the same Configuration-derived projection used
// by business bindings, not bespoke implementations that evade the framework.
#[derive(Clone, Default, ck::Configuration)]
pub struct ProtocolSelection {
    #[configuration(required)]
    pub indices: Vec<usize>,
}
#[derive(Clone, Default, ck::Configuration)]
pub struct ProtocolNested {
    pub child: ck::RingSearchParams,
    pub optional_child: Option<ck::RingSearchParams>,
}
__ck_configuration_ProtocolSelection!(crate::project_parameters, ProtocolSelection, {}, {
    Selection
});
__ck_configuration_ProtocolNested!(crate::project_parameters, ProtocolNested, {}, {
    NestedParams
});

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = toSmilesWithParams)]
    pub fn to_smiles_with_params(&self, params: SmilesWriteParamsInput) -> Result<String, JsValue> {
        self.to_smiles(Some(params))
    }
    #[wasm_bindgen(js_name=toCxSmilesWithParams)]
    pub fn to_cx_smiles_with_params(
        &self,
        params: notation::CxSmilesWriteParamsInput,
    ) -> Result<String, JsValue> {
        self.to_cx_smiles(Some(params))
    }
    #[wasm_bindgen(js_name=toFragmentSmilesWithParams)]
    pub fn to_fragment_smiles_with_params(
        &self,
        params: FragmentSmilesWriteParamsInput,
    ) -> Result<String, JsValue> {
        self.to_fragment_smiles(params)
    }
    #[wasm_bindgen(js_name=toFragmentCxSmilesWithParams)]
    pub fn to_fragment_cx_smiles_with_params(
        &self,
        params: notation::FragmentCxSmilesWriteParamsInput,
    ) -> Result<String, JsValue> {
        self.to_fragment_cx_smiles(params)
    }
    #[wasm_bindgen(js_name = toRandomSmilesWithParams)]
    pub fn to_random_smiles_with_params(
        &self,
        count: Ts<input::U32Input>,
        seed: Ts<input::U32Input>,
        params: RandomSmilesWriteParamsInput,
    ) -> Result<Vec<String>, JsValue> {
        self.to_random_smiles(count, seed, Some(params))
    }
    #[wasm_bindgen(js_name = fromSmilesWithParams)]
    pub fn from_smiles_with_params(
        smiles: Ts<input::TextInput>,
        params: SmilesParseParamsInput,
    ) -> Result<Molecule, JsValue> {
        Self::from_smiles(smiles, Some(params))
    }
}
