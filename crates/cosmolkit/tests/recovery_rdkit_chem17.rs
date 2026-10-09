use cosmolkit::{Molecule, SanitizeOperations, SanitizeParams, SmilesParseParams};
use cosmolkit_core::{KekulizeError, KekulizeParams, RingSearchParams, ValenceModel};
// Exact finite official-source strings, previously probed in RDKit .6 only.
// Recovery-main Rust assertions below remain UNRUN until formal application.
const ROOT7: &str = "c12ccc3cccc4ccc(c1c43)c1c3cccc4c5cccc6c7cccc8c9cccc%10c%11c%12ccc%13cccc%14ccc(c%12c%13%14)c%12c2c1c1c(c43)c(c65)c(c78)c(c9%10)c1c%11%12";
const ROOT30: &str = "c1ccc2c3cccc4c5cccc6c7cccc8c9c%10ccc%11cccc%12ccc(c%10c%11%12)c%10c%11c%12ccc%13cccc%14ccc(c%12c%14%13)c%12c1c2c1c(c43)c(c56)c(c78)c(c9%10)c1c%12%11";
fn raw(s: &str) -> Molecule {
    Molecule::from_smiles_with_params(
        s,
        &SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..SmilesParseParams::default()
        },
    )
    .unwrap()
}
#[test]
fn official_8403_fast_failure_is_retried_canonically_by_registered_sanitize() {
    for input in [ROOT7, ROOT30] {
        let source = raw(input);
        let observer = source.clone();
        let mut topology = source.topology().clone();
        let mut valence = cosmolkit_core::assign_valence_with_options_for_topology(
            &topology,
            ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let mut rings =
            cosmolkit_core::symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
        assert!(
            matches!(
                cosmolkit_core::source_kekulize_attempt(
                    &mut topology,
                    &mut valence,
                    &mut rings,
                    &KekulizeParams {
                        canonical: false,
                        ..KekulizeParams::default()
                    }
                ),
                Err(KekulizeError::NotKekulizable { .. })
            ),
            "source-derived fast branch: {input}"
        );
        let fixed = source
            .sanitize_with_params(&SanitizeParams {
                operations: SanitizeOperations::PROPERTIES
                    | SanitizeOperations::SYMM_RINGS
                    | SanitizeOperations::KEKULIZE,
            })
            .unwrap();
        assert_eq!(fixed.atoms().len(), 62);
        assert!(fixed.bonds().iter().all(|b| !b.is_aromatic()));
        assert_eq!(source, observer);
        assert!(Molecule::from_smiles(input).is_ok());
    }
}
#[test]
fn official_8403_problem_detector_retries_without_reporting_first_failure() {
    for input in [ROOT7, ROOT30] {
        let source = raw(input);
        let observer = source.clone();
        let report = source.detect_chemistry_problems().unwrap();
        assert!(report.problems.is_empty(), "{input}: {report:?}");
        assert_eq!(source, observer);
    }
}
