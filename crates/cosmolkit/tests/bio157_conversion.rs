//! P5 public conversion proposals; independent p1 review/ROOT ruling pending.
#![cfg(all(feature = "cap-bio", feature = "cap-io"))]
use cosmolkit::{AtomId, BioMoleculeError, BioMoleculeParams, BioStructure, Protein};

const PDB: &str = "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C  \nHETATM    2  C2  LIG A   1       1.400   0.000   0.000  1.00 20.00           C  \nHETATM    3  H1  LIG A   1      -1.000   0.000   0.000  1.00 20.00           H  \nCONECT    1    2    3\nEND\n";

#[test]
fn pipeline_parameter_product_uses_one_checked_molecule_boundary() {
    let source = BioStructure::from_pdb(PDB).unwrap();
    let before = source.clone();
    for sanitize in [false, true] {
        for remove_hs in [false, true] {
            for proximity_bonding in [false, true] {
                let molecule = source
                    .to_molecule_with_params(&BioMoleculeParams {
                        sanitize,
                        remove_hs,
                        proximity_bonding,
                        flavor: 0,
                    })
                    .unwrap();
                assert_eq!(
                    molecule.num_atoms(),
                    if sanitize && remove_hs { 2 } else { 3 }
                );
                assert_eq!(
                    molecule.num_bonds(),
                    if sanitize && remove_hs { 1 } else { 2 }
                );
                let first = molecule.atom(AtomId::new(0)).unwrap();
                assert_eq!(first.pdb_residue_info().unwrap().serial_number(), 1);
                assert_eq!(source, before);
            }
        }
    }
    assert_eq!(source.to_molecule().unwrap().num_atoms(), 2);
}

#[test]
fn empty_hierarchy_and_protein_projection_convert_without_mutating_sources() {
    let empty = BioStructure::from_pdb("END\n").unwrap();
    assert_eq!(empty.to_molecule().unwrap().num_atoms(), 0);
    let protein = Protein::from_pdb(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nEND\n",
    )
    .unwrap();
    assert_eq!(protein.to_molecule().unwrap().num_atoms(), 1);
    assert_eq!(protein.num_atoms(), 1);
}

#[test]
fn all_bio_models_and_typed_failure_stay_at_the_value_boundary() {
    let source=BioStructure::from_pdb("MODEL        1\nHETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C  \nENDMDL\nMODEL        2\nHETATM    2  C2  LIG A   1      10.000   0.000   0.000  1.00 20.00           C  \nENDMDL\nEND\n").unwrap();
    let params = BioMoleculeParams {
        sanitize: false,
        remove_hs: false,
        proximity_bonding: false,
        flavor: 0,
    };
    assert_eq!(source.num_models(), 2);
    assert_eq!(
        source.to_molecule_with_params(&params).unwrap().num_atoms(),
        2
    );
    let _: fn(&BioStructure, &BioMoleculeParams) -> Result<cosmolkit::Molecule, BioMoleculeError> =
        BioStructure::to_molecule_with_params;
}
