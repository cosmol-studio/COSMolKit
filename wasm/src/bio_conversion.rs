//! Typed BIO conversion bridges preserve canonical errors and move the resulting value.
use crate::Molecule;
use cosmolkit as ck;
pub fn bio_structure_to_molecule(
    value: &ck::BioStructure,
) -> Result<Molecule, ck::BioMoleculeError> {
    // COSMolKit❗✔️: value.to_molecule()
    value.to_molecule().map(|inner| Molecule {
        inner: inner.into(),
    })
}
pub fn bio_structure_to_molecule_with_params(
    value: &ck::BioStructure,
    params: &ck::BioMoleculeParams,
) -> Result<Molecule, ck::BioMoleculeError> {
    // COSMolKit❗✔️: value.to_molecule_with_params(params)
    value.to_molecule_with_params(params).map(|inner| Molecule {
        inner: inner.into(),
    })
}
pub fn protein_to_molecule(value: &ck::Protein) -> Result<Molecule, ck::BioMoleculeError> {
    // COSMolKit❗✔️: value.to_molecule()
    value.to_molecule().map(|inner| Molecule {
        inner: inner.into(),
    })
}
pub fn protein_to_molecule_with_params(
    value: &ck::Protein,
    params: &ck::BioMoleculeParams,
) -> Result<Molecule, ck::BioMoleculeError> {
    // COSMolKit❗✔️: value.to_molecule_with_params(params)
    value.to_molecule_with_params(params).map(|inner| Molecule {
        inner: inner.into(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;
    const PDB: &str = "ATOM      1  CA AALA A   1       1.000   2.000   3.000  0.60 20.00           C  \nATOM      2  CA BALA A   1       4.000   5.000   6.000  0.40 30.00           C  \nATOM      3  N   ALA A   1       7.000   8.000   9.000  1.00 10.00           N  \nEND\n";
    const HYDROGEN: &str = "ATOM      1  CA AALA A   1       1.000   2.000   3.000  0.60 20.00           C  \nATOM      2  CA BALA A   1       4.000   5.000   6.000  0.40 30.00           C  \nATOM      3  N   ALA A   1       7.000   8.000   9.000  1.00 10.00           N  \nATOM      4  H   ALA A   1       7.800   8.000   9.000  1.00 20.00           H  \nEND\n";
    const INVALID: &str = "ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \nATOM      2  F   ALA A   1       1.300   0.000   0.000  1.00 20.00           F  \nATOM      3  F   ALA A   1      -1.300   0.000   0.000  1.00 20.00           F  \nATOM      4  F   ALA A   1       0.000   1.300   0.000  1.00 20.00           F  \nATOM      5  F   ALA A   1       0.000  -1.300   0.000  1.00 20.00           F  \nATOM      6  F   ALA A   1       0.000   0.000   1.300  1.00 20.00           F  \nCONECT    1    2    3    4    5    6\nEND\n";
    #[test]
    fn bio_conversion_both_types_all_parameters_match_canonical_and_leave_sources_unchanged() {
        for text in [PDB, HYDROGEN, "END\n"] {
            let structure = ck::BioStructure::from_pdb(text).unwrap();
            let protein = ck::Protein::from_pdb(text).unwrap();
            let before = structure.to_pdb().unwrap();
            assert_eq!(
                *bio_structure_to_molecule(&structure)
                    .unwrap()
                    .inner
                    .borrow(),
                structure.to_molecule().unwrap()
            );
            assert_eq!(
                *protein_to_molecule(&protein).unwrap().inner.borrow(),
                protein.to_molecule().unwrap()
            );
            for sanitize in [false, true] {
                for remove_hs in [false, true] {
                    for flavor in [0, 1, u32::MAX] {
                        for proximity_bonding in [false, true] {
                            let params = ck::BioMoleculeParams {
                                sanitize,
                                remove_hs,
                                flavor,
                                proximity_bonding,
                            };
                            assert_eq!(
                                *bio_structure_to_molecule_with_params(&structure, &params)
                                    .unwrap()
                                    .inner
                                    .borrow(),
                                structure.to_molecule_with_params(&params).unwrap()
                            );
                            assert_eq!(
                                *protein_to_molecule_with_params(&protein, &params)
                                    .unwrap()
                                    .inner
                                    .borrow(),
                                protein.to_molecule_with_params(&params).unwrap()
                            );
                        }
                    }
                }
            }
            assert_eq!(structure.to_pdb().unwrap(), before);
        }
    }
    #[test]
    fn bio_conversion_failures_preserve_the_actual_owner_error() {
        let structure = ck::BioStructure::from_pdb(INVALID).unwrap();
        let before = structure.to_pdb().unwrap();
        let canonical = structure.to_molecule().unwrap_err();
        let projected = bio_structure_to_molecule(&structure).err().unwrap();
        assert_eq!(format!("{projected:?}"), format!("{canonical:?}"));
        assert!(matches!(projected, ck::BioMoleculeError::Conversion(_)));
        if let ck::BioMoleculeError::Conversion(ref cause) = projected {
            assert!(matches!(
                cause,
                ck::BioMoleculeConversionError::Hydrogens(_)
            ));
            assert!(cause.source().is_none());
        }
        assert_eq!(structure.to_pdb().unwrap(), before);
        let relaxed = ck::BioMoleculeParams {
            sanitize: false,
            ..Default::default()
        };
        assert_eq!(
            *bio_structure_to_molecule_with_params(&structure, &relaxed)
                .unwrap()
                .inner
                .borrow(),
            structure.to_molecule_with_params(&relaxed).unwrap()
        );
    }
}
