//! Original fixed TAU conditions, preserved at the canonical runtime boundary.
//! Source lineage: historical chemistry/tautomer.rs 5972-6085, pinned RDKit
//! testTautomer.cpp CHEMBL23979/CHEMBL12724; only API projection is adapted.
#![cfg(all(feature = "cap-tautomer", feature = "cap-smiles"))]
use cosmolkit::{AtomId, BondId, Molecule, TautomerEnumerationStatus, TautomerParams};
use std::collections::BTreeSet;
#[test]
fn enumeration_clears_computed_ring_stereo_before_transforming_candidates_like_rdkit() {
    let molecule = Molecule::from_smiles(
        "N[C@H](C(=O)N1CCCC1)[C@H]1CC[C@H](NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1",
    )
    .expect("parse CHEMBL23979 ring-stereo regression");
    let params = TautomerParams::default().with_reassign_stereo(false);

    let result = molecule
        .enumerate_tautomers_with_params(&params)
        .expect("enumerate CHEMBL23979 ring-stereo regression");

    assert_eq!(result.status(), TautomerEnumerationStatus::Completed);
    assert_eq!(
        result.canonical_smiles(),
        [
            "N=C(C1CC[C@H](NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1)C(O)N1CCCC1",
            "NC(=C(O)N1CCCC1)C1CC[C@H](NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1",
            "NC(=C1CC[C@H](NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1)C(O)N1CCCC1",
            "NC(C(=O)N1CCCC1)C1CC[C@H](NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1",
        ]
    );
    assert_eq!(
        result.modified_atoms(),
        &BTreeSet::from([
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3),
            AtomId::new(9),
        ])
    );
    assert_eq!(
        result.modified_bonds(),
        &BTreeSet::from([
            BondId::new(0),
            BondId::new(1),
            BondId::new(2),
            BondId::new(8),
        ])
    );
    assert_eq!(
        result
            .iter()
            .map(|candidate| candidate.tautomer_score().unwrap().total())
            .collect::<Vec<_>>(),
        [251, 250, 250, 255]
    );
    assert_eq!(
        result
            .canonical_tautomer_with_params(&params)
            .expect("select CHEMBL23979 canonical tautomer")
            .to_smiles()
            .expect("write CHEMBL23979 canonical tautomer"),
        "NC(C(=O)N1CCCC1)C1CCC(NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1"
    );
}

#[test]
fn enumeration_rekeys_stale_computed_double_bond_stereo_like_rdkit_chembl12724() {
    let molecule =
        Molecule::from_smiles("COc1ccc(OC)c(/C=N/N=C(\\N)NO)c1.Cc1ccc(S(=O)(=O)O)cc1").unwrap();
    let params = TautomerParams::default().with_reassign_stereo(false);
    let result = molecule.enumerate_tautomers_with_params(&params).unwrap();
    assert_eq!(
        result.canonical_smiles(),
        [
            "COc1ccc(OC)c(/C=N/N=C(N)NO)c1.Cc1ccc(S(=O)(=O)O)cc1",
            "COc1ccc(OC)c(/C=N/NC(=N)NO)c1.Cc1ccc(S(=O)(=O)O)cc1",
            "COc1ccc(OC)c(/C=N/NC(N)=NO)c1.Cc1ccc(S(=O)(=O)O)cc1",
            "COc1ccc(OC)c(/C=N/NC(N)N=O)c1.Cc1ccc(S(=O)(=O)O)cc1",
            "COc1ccc(OC)c(C=NN=C(N)NO)c1.Cc1ccc(S(=O)(=O)O)cc1",
        ]
    );
    assert_eq!(
        result.modified_atoms(),
        &BTreeSet::from([
            AtomId::new(11),
            AtomId::new(12),
            AtomId::new(13),
            AtomId::new(14),
            AtomId::new(15),
        ])
    );
    assert_eq!(
        result.modified_bonds(),
        &BTreeSet::from([
            BondId::new(11),
            BondId::new(12),
            BondId::new(13),
            BondId::new(14),
        ])
    );
    assert_eq!(result.status(), TautomerEnumerationStatus::Completed);
    assert_eq!(
        result
            .iter()
            .map(|candidate| candidate.tautomer_score().unwrap().total())
            .collect::<Vec<_>>(),
        [509, 509, 513, 506, 509]
    );
    assert_eq!(
        result
            .canonical_tautomer_with_params(&params)
            .expect("select CHEMBL12724 canonical tautomer")
            .to_smiles()
            .expect("write CHEMBL12724 canonical tautomer"),
        "COc1ccc(OC)c(C=NNC(N)=NO)c1.Cc1ccc(S(=O)(=O)O)cc1"
    );
}

struct ScoreFn<F>(F);
impl<F> cosmolkit::TautomerScorer for ScoreFn<F>
where
    F: Fn(cosmolkit::TautomerMoleculeView<'_>) -> Result<i32, cosmolkit::TautomerRunError>
        + Send
        + Sync,
{
    fn score(
        &self,
        value: cosmolkit::TautomerMoleculeView<'_>,
    ) -> Result<i32, cosmolkit::TautomerRunError> {
        (self.0)(value)
    }
}
fn score_params(
    f: impl Fn(cosmolkit::TautomerMoleculeView<'_>) -> Result<i32, cosmolkit::TautomerRunError>
    + Send
    + Sync
    + 'static,
) -> TautomerParams {
    let mut value = TautomerParams::default();
    value.set_scorer(Some(std::sync::Arc::new(ScoreFn(f))));
    value
}
#[test]
fn canonical_selection_iterable_computes_lexical_ties_and_matches_result_path() {
    let inputs = ["CCC", "CC", "C"].map(|s| Molecule::from_smiles(s).unwrap());
    let before = inputs.clone();
    let result =
        cosmolkit::canonical_tautomer_from_molecules_with_params(&inputs, &score_params(|_| Ok(9)))
            .unwrap();
    assert_eq!(result.to_smiles().unwrap(), "C");
    assert_eq!(inputs, before);
    let source = Molecule::from_smiles("CC(C)=O").unwrap();
    let result = source.enumerate_tautomers().unwrap();
    let iterable = result.iter().cloned().collect::<Vec<_>>();
    assert_eq!(
        cosmolkit::canonical_tautomer_from_molecules(&iterable).unwrap(),
        result.canonical_tautomer().unwrap()
    );
}
#[test]
fn iterable_empty_minimum_scores_and_single_score_skip_preserve_source_conditions() {
    use cosmolkit::{OperationError, TautomerRunError};
    assert!(matches!(
        cosmolkit::canonical_tautomer_from_molecules(&[]),
        Err(OperationError::Tautomer(
            TautomerRunError::NoCanonicalTautomer
        ))
    ));
    let inputs = ["C", "CC"].map(|s| Molecule::from_smiles(s).unwrap());
    assert!(matches!(
        cosmolkit::canonical_tautomer_from_molecules_with_params(
            &inputs,
            &score_params(|_| Ok(i32::MIN))
        ),
        Err(OperationError::Tautomer(
            TautomerRunError::NoCanonicalTautomer
        ))
    ));
    let result = cosmolkit::canonical_tautomer_from_molecules_with_params(
        &inputs[..1],
        &score_params(|_| panic!("source size-one branch skips scoring")),
    )
    .unwrap();
    assert_eq!(result.to_smiles().unwrap(), "C");
    assert_eq!(result.property("_StereochemDone"), Some("1"));
}
#[test]
fn iterable_retains_duplicate_values_and_first_tie_identity_without_sorting() {
    let tagged = |name: &str| {
        Molecule::from_smiles("C")
            .unwrap()
            .to_builder()
            .with_property("source".into(), name.into())
            .unwrap()
            .build()
            .unwrap()
    };
    let inputs = [
        tagged("first"),
        tagged("second"),
        Molecule::from_smiles("CC").unwrap(),
    ];
    let calls = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let observed = calls.clone();
    let result = cosmolkit::canonical_tautomer_from_molecules_with_params(
        &inputs,
        &score_params(move |_| {
            observed.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            Ok(9)
        }),
    )
    .unwrap();
    assert_eq!(calls.load(std::sync::atomic::Ordering::SeqCst), 3);
    assert_eq!(result.property("source"), Some("first"));
    assert_eq!(inputs[1].property("source"), Some("second"));
}
