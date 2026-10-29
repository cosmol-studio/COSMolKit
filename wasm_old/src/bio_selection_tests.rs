//! Native boundary regressions for BIO COW and structured selection failures.
use crate::{
    BioModelRow, BioOperationError, BioRowModelError, BioRowSpan, BioRowTraverseError,
    BioSelection, BioSelectionCopyCause, BioStructure,
};
use std::error::Error;
const PDB: &str = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nHETATM    2  O   HOH B   2       4.000   5.000   6.000  1.00 20.00           O  \nEND\n";
#[test]
fn bio_selection_reexports_preserve_cow_and_failure_atomicity() {
    let mut s = BioStructure::from_pdb(PDB).unwrap();
    let peer = s.clone();
    let selection = BioSelection::from_cid("/1/A/1").unwrap();
    assert_eq!(
        s.selected_atom_ids(&selection)
            .unwrap()
            .iter()
            .map(|id| id.value())
            .collect::<Vec<_>>(),
        [0]
    );
    let value = s.with_selection(&selection).unwrap();
    assert_eq!(value.num_atoms(), 1);
    s.translate_([1., 2., 3.]).unwrap();
    assert_eq!(peer.coordinates().positions()[0], [1., 2., 3.]);
    assert_eq!(s.coordinates().positions()[0], [2., 4., 6.]);
    s.retain_selection_(&selection).unwrap();
    assert_eq!(s.num_atoms(), 1);
    let mut parts = peer.clone().into_parts();
    parts.models[0] = BioModelRow::new(BioRowSpan::new(0, 2).unwrap(), None);
    let mut malformed = BioStructure::from_parts(parts).unwrap();
    let original = malformed.clone();
    let all = BioSelection::from_cid("/").unwrap();
    let query = malformed.selected_atom_ids(&all).unwrap_err();
    assert!(
        query
            .source()
            .unwrap()
            .downcast_ref::<BioRowTraverseError>()
            .is_some()
    );
    let error = malformed.retain_selection_(&all).unwrap_err();
    let BioOperationError::Selection(copy) = error else {
        panic!("wrong operation error")
    };
    assert!(matches!(
        copy.cause(),
        BioSelectionCopyCause::Traverse(BioRowTraverseError::Model(
            BioRowModelError::MissingModelNumber
        ))
    ));
    assert_eq!(malformed, original);
    assert_eq!(
        malformed.with_selection(&all).unwrap_err().to_string(),
        format!(
            "selection copy failed: {}",
            BioRowModelError::MissingModelNumber
        )
    );
}
