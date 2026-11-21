#![cfg(all(feature = "cap-fragmentation", feature = "cap-smiles"))]
use cosmolkit::{Molecule, RecapLookupError};

fn sorted(mut products: Vec<String>) -> Vec<String> {
    products.sort();
    products
}

#[test]
fn shared_descendant_keeps_both_direct_parents() {
    let source = Molecule::from_smiles("CC(=O)Nc1ccc(NC(=O)C)cc1").unwrap();
    let before = source.to_smiles().unwrap();
    let graph = source.decompose_recap().unwrap();
    assert_eq!(source.to_smiles().unwrap(), before);
    drop(source);
    let root = "CC(=O)Nc1ccc(NC(C)=O)cc1";
    let intermediate = "*Nc1ccc(NC(C)=O)cc1";
    let acetyl = "*C(C)=O";
    let terminal = "*Nc1ccc(N*)cc1";
    assert_eq!(graph.root(), root);
    assert_eq!(
        sorted(graph.children(root).unwrap()),
        [acetyl, intermediate]
    );
    assert_eq!(
        sorted(graph.children(intermediate).unwrap()),
        [acetyl, terminal]
    );
    assert_eq!(sorted(graph.parents(acetyl).unwrap()), [intermediate, root]);
    assert_eq!(sorted(graph.leaves()), [acetyl, terminal]);
    assert!(graph.children(acetyl).unwrap().is_empty());
    assert_eq!(graph.parents(terminal).unwrap(), [intermediate]);
}

#[test]
fn forbidden_sibling_rejects_the_whole_reaction() {
    for (smi, expected) in [
        ("C", vec![]),
        ("CCCCOCCC", vec![]),
        ("CCCCOCCCC", vec!["*CCCC"]),
        ("CCCCOC(C)C", vec!["*C(C)C", "*CCCC"]),
    ] {
        let graph = Molecule::from_smiles(smi)
            .unwrap()
            .decompose_recap()
            .unwrap();
        assert_eq!(
            sorted(graph.children(graph.root()).unwrap()),
            expected,
            "{smi}"
        );
        assert_eq!(sorted(graph.leaves()), expected, "{smi}");
        assert!(graph.parents(graph.root()).unwrap().is_empty());
        for result in [graph.parents("missing"), graph.children("missing")] {
            assert_eq!(
                result.unwrap_err(),
                RecapLookupError::UnknownProduct {
                    product: "missing".into()
                }
            );
        }
    }
}
