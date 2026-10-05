#![cfg(feature = "cap-fingerprints")]
use cosmolkit::{AtomCodeExplanation, AtomPairsParameters as P, BINDING_CONTRACT};

#[test]
fn pinned_atom_pair_encoding_constants_and_width_relations() {
    assert_eq!(P::version(), "1.1.0");
    assert_eq!(
        [
            P::num_type_bits(),
            P::num_pi_bits(),
            P::num_branch_bits(),
            P::num_chiral_bits(),
            P::code_size(),
            P::num_path_bits(),
            P::max_path_length(),
            P::num_atom_pair_fingerprint_bits()
        ],
        [4, 2, 3, 2, 9, 5, 31, 23]
    );
}
#[test]
fn complete_source_atom_types_are_owned_and_include_implicit_zero() {
    let expected = vec![5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53, 0];
    let mut first = P::atom_types();
    assert_eq!(first, expected);
    first[0] = 99;
    first.pop();
    assert_eq!(P::atom_types(), expected);
}
#[test]
fn source_type_order_drives_decoding_and_chirality_masks() {
    let symbols = [
        "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "As", "Se", "Br", "Sb", "Te", "I", "*",
    ];
    for (index, symbol) in symbols.into_iter().enumerate() {
        let code = ((index as u64) << (P::num_branch_bits() + P::num_pi_bits()))
            | 7
            | (3 << P::num_branch_bits());
        let decoded = AtomCodeExplanation::from_code(code, 0, false).unwrap();
        assert_eq!(
            (
                decoded.symbol(),
                decoded.branch_count(),
                decoded.pi_electrons()
            ),
            (symbol, 7, 3)
        );
        for (bits, label) in [(0, ""), (1, "R"), (2, "S")] {
            let decoded =
                AtomCodeExplanation::from_code(code | (bits << P::code_size()), 0, true).unwrap();
            assert_eq!(decoded.chirality(), Some(label));
        }
    }
}
#[test]
fn canonical_vocabulary_has_one_type_and_ten_registered_reads() {
    assert_eq!(
        BINDING_CONTRACT
            .iter()
            .filter(|r| r.semantic_id == "types.AtomPairsParameters")
            .count(),
        1
    );
    let rows: Vec<_> = BINDING_CONTRACT
        .iter()
        .filter(|r| r.semantic_id.starts_with("AtomPairsParameters."))
        .collect();
    assert_eq!(rows.len(), 10);
    for name in [
        "version",
        "num_type_bits",
        "num_pi_bits",
        "num_branch_bits",
        "num_chiral_bits",
        "code_size",
        "num_path_bits",
        "max_path_length",
        "num_atom_pair_fingerprint_bits",
        "atom_types",
    ] {
        let row = rows
            .iter()
            .find(|r| r.semantic_id == format!("AtomPairsParameters.{name}"))
            .unwrap();
        assert_eq!(row.python_name, name);
        assert_eq!(row.feature, "cap-fingerprints");
    }
}
