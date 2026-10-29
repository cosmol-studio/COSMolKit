//! Binding crate re-export checks for hierarchy traversal and exact source lookups.
use crate::{
    AltLocLabel, AltLocRequest, AtomName, BioAtomId, BioResidueId, BioStructure, BioStructureError,
};
const PDB: &str = "ATOM      1  CA AALA A   1       1.000   2.000   3.000  0.60 20.00           C  \nATOM      2  CA BALA A   1       4.000   5.000   6.000  0.40 30.00           C  \nATOM      3  N   ALA A   1       7.000   8.000   9.000  1.00 10.00           N  \nEND\n";
#[test]
fn hierarchy_reexports_preserve_positions_exact_altloc_and_error_context() {
    let s = BioStructure::from_pdb(PDB).unwrap();
    assert_eq!(s.num_atoms(), 3);
    assert_eq!(
        s.coordinates().positions(),
        &[[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]]
    );
    assert_eq!(s.atom_position(BioAtomId::new(u32::MAX)), None);
    assert!(s.residue_atoms(BioResidueId::new(9)).is_none());
    let name = AtomName::from_ascii(b" CA ").unwrap();
    assert_eq!(
        s.find_atom(BioResidueId::new(0), name, AltLocRequest::Any, None)
            .unwrap()
            .0
            .value(),
        0
    );
    assert_eq!(
        s.atom_by_altloc(BioResidueId::new(0), name, Some(AltLocLabel::new(b'B')))
            .unwrap()
            .0
            .value(),
        1
    );
    assert!(matches!(
        s.atom_by_altloc(BioResidueId::new(0), name, None),
        Err(BioStructureError::AtomNotFound)
    ));
    assert!(matches!(
        s.atom_by_altloc(BioResidueId::new(99), name, None),
        Err(BioStructureError::RowReferenceOutOfBounds {
            table: "residues",
            index: 99,
            table_len: 1
        })
    ));
}
