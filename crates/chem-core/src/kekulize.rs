use crate::Molecule;

/// Kekulization errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KekulizeError {
    ImpossibleAromaticAssignment,
    NotImplemented,
}

/// Convert aromatic bond representation to one concrete alternating form.
pub fn kekulize_in_place(molecule: &mut Molecule) -> Result<(), KekulizeError> {
    let _ = molecule;
    unimplemented!("kekulization is not implemented yet")
}
