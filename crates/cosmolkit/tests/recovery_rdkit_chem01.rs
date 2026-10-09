use cosmolkit::Molecule;
#[test]
fn source_9398_natural_macrocycle_donors() {
    for (s, indices, expected) in [
        (
            "O=C(O)c1cccc2Oc3cncc(n3)Oc3c(C(=O)O)cccc3Oc3cncc(n3)Oc12",
            &[8, 15, 25, 32][..],
            false,
        ),
        ("O=c1ccccc(=O)c(=O)o1", &[1][..], true),
        ("O=c1ccccc(=O)ooo1", &[1][..], false),
    ] {
        let m = Molecule::from_smiles(s).unwrap();
        for &i in indices {
            assert_eq!(m.atoms()[i].is_aromatic(), expected, "{s} atom{i}");
        }
    }
}
