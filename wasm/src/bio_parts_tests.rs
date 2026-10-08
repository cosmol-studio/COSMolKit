//! Native binding re-export checks for complete parts retention and validation failures.
use crate::{BioCoordinateBlock, BioStructure, BioStructureError};
#[test]
fn bio_parts_reexports_preserve_every_block_and_original_failure() {
    let s = BioStructure::from_pdb("REMARK   1 SOURCE RETAINED\nATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nEND\n").unwrap();
    let parts = s.clone().into_parts();
    BioStructure::validate_parts(&parts).unwrap();
    let restored = BioStructure::from_parts(parts.clone()).unwrap();
    restored.validate().unwrap();
    assert_eq!(restored.clone().into_parts(), parts);
    assert_eq!(s, restored);
    let mut invalid = parts.clone();
    invalid.coordinates = BioCoordinateBlock::new(Vec::new());
    let expected = BioStructureError::CoordinateCountMismatch {
        atom_count: 1,
        coordinate_count: 0,
    };
    assert_eq!(
        BioStructure::validate_parts(&invalid),
        Err(expected.clone())
    );
    assert_eq!(BioStructure::from_parts(invalid).unwrap_err(), expected);
    BioStructure::validate_parts(&parts).unwrap();
    s.validate().unwrap();
}
