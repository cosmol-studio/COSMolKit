//! Reaction-ready intermediates shared by recursive fragmentation queries.
use crate::{ReactionInput, ReactionProduct};
use cosmolkit_core::{RingInfo, SanitizeError, ValenceAssignment, sanitize_owned_topology};
use cosmolkit_smiles::{SmilesParseError, SmilesRecord, SmilesRecordView, write_smiles};

pub(super) struct Product {
    record: SmilesRecord,
    rings: Option<RingInfo>,
    valence: Option<ValenceAssignment>,
}

impl Product {
    pub(super) fn sanitize(product: ReactionProduct) -> Result<Self, SanitizeError> {
        // RDKit source: Code/GraphMol/MolOps.cpp, sanitizeMol
        // RDKit❗✔️:   mol.clearComputedProps();
        let mut properties = product.properties;
        properties.clear_computed_props()?;
        let mut sanitized = sanitize_owned_topology(product.topology, &Default::default())?;

        // Drop the previous reaction's atom provenance before the SMILES writer
        // copies this topology. Each subsequent reaction assigns its own indices.
        for atom in &mut sanitized.topology.atoms {
            if atom.prop("react_atom_idx").is_some() {
                atom.clear_prop("react_atom_idx")?;
            }
            if atom.prop("react_idx").is_some() {
                atom.clear_prop("react_idx")?;
            }
        }

        Ok(Self {
            record: SmilesRecord {
                topology: sanitized.topology,
                coordinates: product.coordinates,
                properties,
            },
            rings: sanitized.final_rings,
            valence: sanitized.final_valence,
        })
    }

    pub(super) fn as_input(&self) -> ReactionInput<'_> {
        ReactionInput {
            topology: &self.record.topology,
            properties: &self.record.properties,
            coordinates: &self.record.coordinates,
            rings: self.rings.as_ref(),
            valence: self.valence.as_ref(),
        }
    }
}

pub(super) enum Reactant<'a> {
    Source(ReactionInput<'a>),
    Product(Product),
}

impl<'a> Reactant<'a> {
    pub(super) fn as_input(&self) -> ReactionInput<'_> {
        match self {
            Self::Source(input) => *input,
            Self::Product(product) => product.as_input(),
        }
    }
}

pub(super) fn canonical_key<E>(input: ReactionInput<'_>) -> Result<String, E>
where
    E: From<SmilesParseError> + From<std::string::FromUtf8Error>,
{
    let bytes = write_smiles(SmilesRecordView {
        topology: input.topology,
        coordinates: input.coordinates,
        properties: input.properties,
        rings: input.rings,
    })?;
    Ok(String::from_utf8(bytes.into_bytes())?)
}
