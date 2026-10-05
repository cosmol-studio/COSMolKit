#![cfg(feature = "cap-fingerprints")]
use cosmolkit::{Molecule, TopologicalTorsionPathScoreError as E, explain_path_score};
#[test]
fn pinned_documented_three_atom_scores_reverse_and_explanations() {
    for (smiles, expected, decoded) in [
        (
            "C=CC",
            10506272,
            vec![("C", 1, 0), ("C", 2, 1), ("C", 1, 1)],
        ),
        (
            "C=CO",
            25186344,
            vec![("C", 1, 1), ("C", 2, 1), ("O", 1, 0)],
        ),
    ] {
        let mol = Molecule::from_smiles(smiles).unwrap();
        assert_eq!(
            mol.topological_torsion_path_score(&[0, 1, 2], 3, None)
                .unwrap(),
            expected
        );
        assert_eq!(
            mol.topological_torsion_path_score(&[2, 1, 0], 3, None)
                .unwrap(),
            expected
        );
        assert_eq!(
            mol.topological_torsion_path_score(&[0, 1, 2, 999], 3, Some(&[]))
                .unwrap(),
            expected
        );
        assert_eq!(explain_path_score(expected, 3), decoded);
    }
}
#[test]
fn nonconnected_repeated_indices_and_custom_atoms_retain_source_semantics() {
    let mol = Molecule::from_smiles("C.C.O").unwrap();
    let expected = 32 | (96 << 9) | (32 << 18) | (96u64 << 27);
    assert_eq!(
        mol.topological_torsion_path_score(&[0, 2, 0, 2], 4, None)
            .unwrap(),
        expected
    );
    assert_eq!(
        mol.topological_torsion_path_score(&[2, 0, 2, 0], 4, None)
            .unwrap(),
        expected
    );
    let two = Molecule::from_smiles("CC").unwrap();
    assert_eq!(
        two.topological_torsion_path_score(&[0, 1], 2, Some(&[1, 512]))
            .unwrap(),
        261632
    );
    assert_eq!(
        explain_path_score(261632, 2),
        vec![("B", 1, 0), ("*", 8, 3)]
    );
}
#[test]
fn original_validation_order_and_context_are_structured() {
    let mol = Molecule::from_smiles("CCC").unwrap();
    assert!(matches!(
        mol.topological_torsion_path_score(&[], 0, None),
        Err(E::ZeroSize)
    ));
    assert!(matches!(
        mol.topological_torsion_path_score(&[0], 2, None),
        Err(E::ShortPath {
            actual: 1,
            required: 2
        })
    ));
    assert!(matches!(
        mol.topological_torsion_path_score(&[9], 1, Some(&[4])),
        Err(E::ShortAtomCodes {
            actual: 1,
            required: 3
        })
    ));
    assert!(matches!(
        mol.topological_torsion_path_score(&[3], 1, Some(&[4; 3])),
        Err(E::AtomIndexOutOfRange {
            index: 3,
            atom_count: 3
        })
    ));
    assert!(matches!(
        mol.topological_torsion_path_score(&[0, 1, 2], 3, Some(&[10, 1, 10])),
        Err(E::AtomCodeUnderflow {
            index: 1,
            code: 1,
            subtract: 2
        })
    ));
    assert!(matches!(
        mol.topological_torsion_path_score(&[0; 9], 9, Some(&[10; 3])),
        Err(E::PackedCode(_))
    ));
}
#[test]
fn full_unsigned_transport_and_zero_chunks_past_eight_are_preserved() {
    assert!(explain_path_score(u64::MAX, 0).is_empty());
    assert_eq!(
        explain_path_score(0, 4),
        vec![("B", 1, 0), ("B", 2, 0), ("B", 2, 0), ("B", 1, 0)]
    );
    let decoded = explain_path_score(u64::MAX, 10);
    assert_eq!(decoded[0], ("*", 8, 3));
    assert_eq!(decoded[1], ("*", 9, 3));
    assert_eq!(decoded[8], ("B", 2, 0));
    assert_eq!(decoded[9], ("B", 1, 0));
}
