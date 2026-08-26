use cosmolkit_core::{
    AtomPairFingerprintParams, AvalonFingerprintParams, MaccsFingerprintParams, Molecule,
    MorganFingerprintParams, PatternFingerprintParams, TopologicalFingerprintParams,
    TopologicalTorsionFingerprintParams, topological_torsion_fingerprint,
};
use rayon::prelude::*;

#[derive(Clone, Debug, PartialEq, Eq)]
struct Snapshot {
    smiles: String,
    avalon: Vec<usize>,
    topological: Vec<usize>,
    morgan: Vec<usize>,
    atom_pair: Vec<usize>,
    pattern: Vec<usize>,
    tautomeric_pattern: Vec<usize>,
    maccs: Vec<usize>,
    topological_torsion: Vec<usize>,
}

fn snapshot(molecule: &Molecule) -> Snapshot {
    Snapshot {
        smiles: molecule.to_smiles(true).expect("snapshot SMILES"),
        avalon: molecule
            .avalon_fingerprint(&AvalonFingerprintParams::default())
            .expect("Avalon fingerprint")
            .on_bits(),
        topological: molecule
            .topological_fingerprint(&TopologicalFingerprintParams::default())
            .expect("topological fingerprint")
            .on_bits(),
        morgan: molecule
            .morgan_fingerprint(&MorganFingerprintParams::default())
            .expect("Morgan fingerprint")
            .on_bits(),
        atom_pair: molecule
            .atom_pair_fingerprint(&AtomPairFingerprintParams::default())
            .expect("AtomPair fingerprint")
            .on_bits(),
        pattern: molecule
            .pattern_fingerprint(&PatternFingerprintParams::default())
            .expect("Pattern fingerprint")
            .on_bits(),
        tautomeric_pattern: molecule
            .pattern_fingerprint(&PatternFingerprintParams {
                n_bits: 257,
                tautomeric: true,
            })
            .expect("tautomeric Pattern fingerprint")
            .on_bits(),
        maccs: molecule
            .maccs_fingerprint(&MaccsFingerprintParams::default())
            .expect("MACCS fingerprint")
            .on_bits(),
        topological_torsion: topological_torsion_fingerprint(
            molecule,
            &TopologicalTorsionFingerprintParams::default(),
        )
        .expect("Topological Torsion fingerprint")
        .on_bits(),
    }
}

#[test]
fn fingerprint_families_are_order_independent_and_non_mutating() {
    let smiles = [
        "CCO",
        "c1ccccc1O",
        "C[C@H](N)C(=O)O",
        "[2H]O[2H]",
        "C1CCCCC1",
    ];
    let molecules = smiles
        .iter()
        .map(|value| Molecule::from_smiles(value).expect("fixture molecule"))
        .collect::<Vec<_>>();
    let before = molecules.iter().map(snapshot).collect::<Vec<_>>();

    for molecule in &molecules {
        let _ = molecule.avalon_fingerprint(&AvalonFingerprintParams::default());
        let _ = molecule.pattern_fingerprint(&PatternFingerprintParams {
            n_bits: 257,
            tautomeric: true,
        });
        let _ = molecule.morgan_fingerprint(&MorganFingerprintParams::default());
        let _ = molecule.atom_pair_fingerprint(&AtomPairFingerprintParams::default());
        let _ = molecule.topological_fingerprint(&TopologicalFingerprintParams::default());
        let _ = topological_torsion_fingerprint(
            molecule,
            &TopologicalTorsionFingerprintParams::default(),
        );
        let _ = molecule.maccs_fingerprint(&MaccsFingerprintParams::default());
        let _ = molecule.pattern_fingerprint(&PatternFingerprintParams::default());
        let _ = molecule.maccs_fingerprint(&MaccsFingerprintParams::default());
        let _ = molecule.avalon_fingerprint(&AvalonFingerprintParams::default());
    }
    let after = molecules.iter().map(snapshot).collect::<Vec<_>>();
    assert_eq!(
        before, after,
        "fingerprint calls changed molecule state or caches"
    );
}

#[test]
fn fingerprint_families_are_deterministic_under_parallel_repetition() {
    let molecules = ["CCO", "c1ccccc1O", "C[C@H](N)C(=O)O", "C1CCCCC1"]
        .into_par_iter()
        .map(|smiles| Molecule::from_smiles(smiles).expect("fixture molecule"))
        .collect::<Vec<_>>();
    let expected = molecules.iter().map(snapshot).collect::<Vec<_>>();
    for _ in 0..4 {
        let actual = molecules.par_iter().map(snapshot).collect::<Vec<_>>();
        assert_eq!(actual, expected, "parallel fingerprint output drifted");
    }
}
