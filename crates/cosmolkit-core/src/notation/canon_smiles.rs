use crate::{Molecule, SmilesWriteError, SmilesWriteParams};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomColor {
    Initial,
    Unique,
    Duplicate,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MolStackType {
    Atom,
    Bond,
    Branch,
    RingClosure,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MolStackElem {
    Atom(usize),
    Bond(usize),
    Branch,
    RingClosure(usize),
}

impl MolStackElem {
    #[must_use]
    pub const fn kind(&self) -> MolStackType {
        match self {
            Self::Atom(_) => MolStackType::Atom,
            Self::Bond(_) => MolStackType::Bond,
            Self::Branch => MolStackType::Branch,
            Self::RingClosure(_) => MolStackType::RingClosure,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct FragmentTraversal {
    pub atoms: Vec<usize>,
    pub bonds: Vec<usize>,
}

pub fn rank_mol_atoms(molecule: &Molecule) -> Result<Vec<usize>, SmilesWriteError> {
    crate::canon_rank::rank_mol_atoms(molecule).map_err(Into::into)
}

pub fn build_noncanonical_fragment(
    molecule: &Molecule,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    let mut writer_params = params.clone();
    writer_params.canonical = false;
    super::smiles_write::mol_to_smiles(molecule, &writer_params)
}

pub fn canonicalize_fragment(
    molecule: &Molecule,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    let mut writer_params = params.clone();
    writer_params.canonical = true;
    writer_params.do_random = false;
    super::smiles_write::mol_to_smiles(molecule, &writer_params)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_noncanonical_fragment_forces_noncanonical_writer_path() {
        let molecule = Molecule::from_smiles_with_sanitize("OC.C", false).unwrap();
        let params = SmilesWriteParams {
            canonical: true,
            do_random: false,
            ..Default::default()
        };

        let mut expected_params = params.clone();
        expected_params.canonical = false;
        let expected =
            crate::notation::smiles_write::mol_to_smiles(&molecule, &expected_params).unwrap();

        let actual = build_noncanonical_fragment(&molecule, &params).unwrap();
        assert_eq!(actual, expected);
    }

    #[test]
    fn build_noncanonical_fragment_does_not_mutate_caller_params() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let params = SmilesWriteParams {
            canonical: true,
            ..Default::default()
        };

        let _ = build_noncanonical_fragment(&molecule, &params).unwrap();
        assert!(params.canonical);
    }

    #[test]
    fn canonicalize_fragment_forces_canonical_writer_path() {
        let molecule = Molecule::from_smiles_with_sanitize("OC.C", false).unwrap();
        let params = SmilesWriteParams {
            canonical: false,
            do_random: true,
            ..Default::default()
        };

        let mut expected_params = params.clone();
        expected_params.canonical = true;
        expected_params.do_random = false;
        let expected =
            crate::notation::smiles_write::mol_to_smiles(&molecule, &expected_params).unwrap();

        let actual = canonicalize_fragment(&molecule, &params).unwrap();
        assert_eq!(actual, expected);
    }

    #[test]
    fn canonicalize_fragment_does_not_mutate_caller_params() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let params = SmilesWriteParams {
            canonical: false,
            do_random: true,
            ..Default::default()
        };

        let _ = canonicalize_fragment(&molecule, &params).unwrap();
        assert!(!params.canonical);
        assert!(params.do_random);
    }
}
