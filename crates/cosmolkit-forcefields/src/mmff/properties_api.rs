use cosmolkit_model::TopologyBlock;

use super::mol_properties::{
    MmffAtomProperties, MmffMolProperties, MmffMolPropertiesError, MmffVariant,
};

/// Options for MMFF atom typing and charge assignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmffPropertiesParams {
    pub mmff_variant: String,
}

impl Default for MmffPropertiesParams {
    fn default() -> Self {
        Self {
            mmff_variant: "MMFF94".to_owned(),
        }
    }
}

/// Detached MMFF atom types and charges, in the input atom order.
#[derive(Debug, Clone, PartialEq)]
pub struct MmffProperties {
    valid: bool,
    variant: MmffVariant,
    atoms: Vec<MmffAtomProperties>,
}

impl MmffProperties {
    pub fn is_valid(&self) -> bool {
        self.valid
    }
    pub fn variant(&self) -> MmffVariant {
        self.variant
    }
    pub fn atoms(&self) -> &[MmffAtomProperties] {
        &self.atoms
    }
    pub fn atom_type(&self, index: usize) -> Result<u8, MmffMolPropertiesError> {
        Ok(self.atom(index)?.atom_type)
    }
    pub fn formal_charge(&self, index: usize) -> Result<f64, MmffMolPropertiesError> {
        Ok(self.atom(index)?.formal_charge)
    }
    pub fn partial_charge(&self, index: usize) -> Result<f64, MmffMolPropertiesError> {
        Ok(self.atom(index)?.partial_charge)
    }
    fn atom(&self, index: usize) -> Result<&MmffAtomProperties, MmffMolPropertiesError> {
        self.atoms
            .get(index)
            .ok_or(MmffMolPropertiesError::AtomIndexOutOfRange {
                atom_index: index,
                atoms: self.atoms.len(),
            })
    }
}

pub fn mmff_properties(
    topology: &TopologyBlock,
    mmff_sanitized: bool,
    params: &MmffPropertiesParams,
    supplied_rings: Option<&cosmolkit_core::RingInfo>,
) -> Result<MmffProperties, MmffMolPropertiesError> {
    // Behavior: one properties constructor owns atom typing and charge
    // computation. Return its source validity/variant and move ordered rows.
    // Complexity: no second typing traversal or charge-vector projection.
    let properties = MmffMolProperties::new_prepared(
        topology,
        mmff_sanitized,
        &params.mmff_variant,
        0,
        supplied_rings,
    )?;
    Ok(MmffProperties {
        valid: properties.valid,
        variant: properties.variant,
        atoms: properties.atom_properties,
    })
}

pub fn mmff_has_all_molecule_params(
    topology: &TopologyBlock,
    mmff_sanitized: bool,
    supplied_rings: Option<&cosmolkit_core::RingInfo>,
) -> Result<bool, MmffMolPropertiesError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8,
    // Code/GraphMol/ForceFieldHelpers/Wrap/rdForceFields.cpp:173-177.
    // RDKit❗❌: bool MMFFHasAllMoleculeParams(const ROMol &mol) {
    // RDKit❗❌:   ROMol molCopy(mol);
    // RDKit❗❌:   MMFF::MMFFMolProperties mmffMolProperties(molCopy);
    // RDKit❗❌:   return mmffMolProperties.isValid();
    // RDKit❗❌: }
    // Behavior: source construction owns a detached topology copy, preserving
    // the live input on both success and error. Untyped atoms yield false;
    // preparation/table errors remain typed and are never swallowed.
    // Complexity: the existing detached core preparation may clone topology
    // at individual stages; this is more allocation than the source copy.
    Ok(
        MmffMolProperties::new_prepared(topology, mmff_sanitized, "MMFF94", 0, supplied_rings)?
            .is_valid(),
    )
}
