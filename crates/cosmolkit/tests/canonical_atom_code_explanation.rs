//! Proposed fixed public-boundary regressions from pinned Utils.py.
use cosmolkit::{AtomCodeExplanation, AtomCodeExplanationError};

#[test]
fn source_documented_carbon_oxygen_and_chiral_examples() {
    for (code, symbol, branches, pi) in [
        (41, "C", 1, 1),
        (42, "C", 2, 1),
        (43, "C", 3, 1),
        (105, "O", 1, 1),
        (97, "O", 1, 0),
    ] {
        let result = AtomCodeExplanation::from_code(code, 0, false).unwrap();
        assert_eq!(
            (
                result.symbol(),
                result.branch_count(),
                result.pi_electrons()
            ),
            (symbol, branches, pi)
        );
        assert_eq!(result.chirality(), None);
    }
    for (code, chirality) in [(35, ""), (547, "R"), (1059, "S")] {
        let result = AtomCodeExplanation::from_code(code, 0, true).unwrap();
        assert_eq!(result.symbol(), "C");
        assert_eq!(result.branch_count(), 3);
        assert_eq!(result.pi_electrons(), 0);
        assert_eq!(result.chirality(), Some(chirality));
    }
}

#[test]
fn source_unused_branch_subtract_full_array_and_high_bits() {
    for branch_subtract in [i64::MIN, -2, 0, 1, 2, i64::MAX] {
        let result = AtomCodeExplanation::from_code(481, branch_subtract, true).unwrap();
        assert_eq!(
            (
                result.symbol(),
                result.branch_count(),
                result.pi_electrons(),
                result.chirality()
            ),
            ("*", 1, 0, Some(""))
        );
    }
    assert_eq!(
        AtomCodeExplanation::from_code(547 | (1 << 63), 0, true).unwrap(),
        AtomCodeExplanation::from_code(547, 0, true).unwrap()
    );
    let maximal = AtomCodeExplanation::from_code(u64::MAX, 0, false).unwrap();
    assert_eq!(
        (
            maximal.symbol(),
            maximal.branch_count(),
            maximal.pi_electrons(),
            maximal.chirality()
        ),
        ("*", 7, 3, None)
    );
}

#[test]
fn source_unknown_chirality_is_an_error_only_when_requested() {
    assert_eq!(
        AtomCodeExplanation::from_code(1536, 0, true),
        Err(AtomCodeExplanationError::UnknownChirality { code: 3 })
    );
    assert_eq!(
        AtomCodeExplanation::from_code(1536, 0, false)
            .unwrap()
            .symbol(),
        "B"
    );
    assert_eq!(
        AtomCodeExplanation::from_code(u64::MAX, 0, true)
            .unwrap_err()
            .to_string(),
        "3"
    );
}
