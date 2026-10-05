#![cfg(all(feature = "cap-smiles", feature = "cap-sanitize"))]
use cosmolkit::{Molecule, SanitizeOperations, SanitizeParams, SmilesParseParams};
pub const CASES: [(usize, &str); 3] = [
    (
        2508,
        "CC(C)(c1ccc(OCC[C@H](O)Cl)cc1)c1ccc(OCC[C@@H](O)Cl)cc1",
    ),
    (
        227,
        "CCC/C=C/C=C/C(O)=N[C@@H](Cc1c[nH]c2cc(Cl)ccc12)C(O)=N[C@H](CC(=N)O)C(O)=N[C@@H](CC(=O)O)C(O)=N[C@@H]1C(O)=NCC(O)=N[C@@H](CCCN)C(O)=N[C@@H](CC(=O)O)C(O)=N[C@H](C)C(O)=N[C@@H](CC(=O)O)C(O)=NCC(O)=N[C@H](C)C(O)=N[C@@H]([C@H](C)CC(=O)O)C(O)=N[C@@H](CC(=O)c2ccc(Cl)cc2N)C(=O)O[C@@H]1C",
    ),
    (
        155,
        "CCC[C@H]1C[C@@H]2CC[C@@H](O2)[C@H](C)C(=O)O[C@@H]([C@H](CC)[C@H]2CC[C@@H](C[C@@H](CCC)N(C)C)O2)[C@H](C)[C@@H]2CC[C@@H](O2)[C@@H](CC)C(=O)O1",
    ),
];
#[test]
fn sanitize_num_arom_constructor_and_operation_transport() {
    for (line, input) in CASES {
        let expected = match line {
            2508 => "2",
            227 => "3",
            155 => "0",
            _ => panic!("unexpected source row"),
        };
        for remove_hydrogens in [false, true] {
            let output = Molecule::from_smiles_with_params(
                input,
                &SmilesParseParams {
                    remove_hydrogens,
                    ..SmilesParseParams::default()
                },
            )
            .unwrap();
            assert_eq!(
                output.properties().prop("numArom"),
                Some(expected),
                "line {line}, remove_hydrogens={remove_hydrogens}"
            );
            assert!(output.properties().is_prop_computed("numArom"));
        }
        let source = Molecule::from_smiles_with_params(
            input,
            &SmilesParseParams {
                sanitize: false,
                remove_hydrogens: false,
                ..SmilesParseParams::default()
            },
        )
        .unwrap();
        let peer = source.clone();
        let output = source
            .sanitize_with_params(&SanitizeParams::default())
            .unwrap();
        assert_eq!(
            output.properties().prop("numArom"),
            Some(expected),
            "line {line}, sanitize operation"
        );
        assert!(output.properties().is_prop_computed("numArom"));
        assert_eq!(source, peer);
        assert_eq!(source.properties().prop("numArom"), None);
        let cleared = output
            .sanitize_with_params(&SanitizeParams {
                operations: SanitizeOperations::NONE,
            })
            .unwrap();
        assert_eq!(cleared.properties().prop("numArom"), None);
        assert_eq!(output.properties().prop("numArom"), Some(expected));
    }
}
