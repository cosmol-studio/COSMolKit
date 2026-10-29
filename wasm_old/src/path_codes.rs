//! Canonical path scoring and value-returning atom-code projection.
use crate::Molecule;
use cosmolkit as ck;
pub struct AtomPairAtomCodeResult {
    inner: ck::AtomPairAtomCodeResult,
}
impl AtomPairAtomCodeResult {
    pub fn code(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.code
        self.inner.code
    }
    pub fn molecule(&self) -> Molecule {
        // COSMolKit❗✔️: inner: self.inner.molecule.clone(),
        Molecule {
            inner: self.inner.molecule.clone().into(),
        }
    }
}
impl Molecule {
    pub fn topological_torsion_path_score(
        &self,
        path: &[usize],
        size: usize,
        atom_codes: Option<&[u32]>,
    ) -> Result<u64, ck::TopologicalTorsionPathScoreError> {
        // COSMolKit❗✔️: .topological_torsion_path_score(&path, size, atom_codes.as_deref())
        self.inner
            .borrow()
            .topological_torsion_path_score(path, size, atom_codes)
    }
    pub fn with_atom_pair_atom_code(
        &self,
        atom_id: usize,
        branch_subtract: u32,
        include_chirality: bool,
        use_legacy_stereo_perception: bool,
    ) -> Result<AtomPairAtomCodeResult, ck::OperationError> {
        // COSMolKit❗✔️: ck::AtomId::new(atom_id),
        // COSMolKit❗✔️: branch_subtract,
        // COSMolKit❗✔️: include_chirality,
        // COSMolKit❗✔️: use_legacy_stereo_perception,
        self.inner
            .borrow()
            .with_atom_pair_atom_code(
                ck::AtomId::new(atom_id),
                branch_subtract,
                include_chirality,
                use_legacy_stereo_perception,
            )
            .map(|inner| AtomPairAtomCodeResult { inner })
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn path_codes_forward_full_defaults_values_errors_and_atomicity() {
        let m = Molecule::from_smiles("CCCC").unwrap();
        let core = ck::Molecule::from_smiles("CCCC").unwrap();
        for (path, size, codes) in [
            (vec![0, 1, 2, 3], 4, None),
            (vec![3, 2, 1, 0], 4, Some(vec![33, 34, 34, 33])),
            (vec![0, 1, 2, 3, 999], 4, Some(vec![])),
            (vec![0, 1, 2, 3, 0, 1, 2, 3], 8, None),
        ] {
            assert_eq!(
                m.topological_torsion_path_score(&path, size, codes.as_deref())
                    .unwrap(),
                core.topological_torsion_path_score(&path, size, codes.as_deref())
                    .unwrap()
            );
        }
        for (path, size, codes) in [
            (vec![], 0, None),
            (vec![0], 4, None),
            (vec![0, 1, 2, 3], 4, Some(vec![33])),
            (vec![99], 1, None),
            (vec![0], 1, Some(vec![0, 0, 0, 0])),
            (vec![0, 1, 2, 3, 0, 1, 2, 3, 0], 9, None),
        ] {
            let a = m
                .topological_torsion_path_score(&path, size, codes.as_deref())
                .unwrap_err();
            let b = core
                .topological_torsion_path_score(&path, size, codes.as_deref())
                .unwrap_err();
            assert_eq!(format!("{a:?}"), format!("{b:?}"));
        }
        for text in ["CCO", "F[C@](Cl)(Br)I"] {
            let m = Molecule::from_smiles(text).unwrap();
            let core = ck::Molecule::from_smiles(text).unwrap();
            let before = m.to_smiles().unwrap();
            for atom in 0..usize::try_from(m.num_atoms()).unwrap() {
                for (subtract, chiral, legacy) in
                    [(0, false, true), (1, false, true), (0, true, false)]
                {
                    let got = m.with_atom_pair_atom_code(atom, subtract, chiral, legacy);
                    let expected = core.with_atom_pair_atom_code(
                        ck::AtomId::new(atom),
                        subtract,
                        chiral,
                        legacy,
                    );
                    match (got, expected) {
                        (Ok(a), Ok(b)) => {
                            assert_eq!(a.code(), b.code);
                            assert_eq!(
                                a.molecule().to_smiles().unwrap(),
                                b.molecule.to_smiles().unwrap()
                            );
                        }
                        (Err(a), Err(b)) => assert_eq!(format!("{a:?}"), format!("{b:?}")),
                        _ => panic!("facade discrepancy"),
                    };
                    assert_eq!(m.to_smiles().unwrap(), before);
                }
            }
            assert!(m.with_atom_pair_atom_code(99, 0, false, true).is_err());
            assert_eq!(m.to_smiles().unwrap(), before);
        }
    }
}
