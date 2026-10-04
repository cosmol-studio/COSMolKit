use cosmolkit::Molecule;

#[test]
fn uff_param_p07_rdkit_2026_03_1_fixed_twenty_query_rows() {
    // Fixed external RDKit 2026.03.1 expectations from the packet. Each input
    // is exercised both as parsed and after the public AddHs operation.
    let cases: [(&str, [bool; 2]); 10] = [
        ("", [true, true]),
        ("CC", [true, true]),
        ("CCO", [true, true]),
        ("c1ccccc1", [true, true]),
        ("CC(=O)O", [true, true]),
        ("N", [true, true]),
        ("[Na+].[Cl-]", [true, true]),
        ("*", [false, false]),
        ("[Cu+]", [false, false]),
        ("[SiH4]", [true, true]),
    ];

    let mut query_calls = 0;
    for (smiles, expected_by_hydrogen_option) in cases {
        for (add_hydrogens, expected) in [
            (false, expected_by_hydrogen_option[0]),
            (true, expected_by_hydrogen_option[1]),
        ] {
            let parsed = Molecule::from_smiles(smiles)
                .unwrap_or_else(|error| panic!("fixed SMILES {smiles:?}: {error}"));
            let with_requested_hydrogens = if add_hydrogens {
                parsed
                    .with_hydrogens()
                    .unwrap_or_else(|error| panic!("AddHs for {smiles:?}: {error}"))
            } else {
                parsed
            };
            // Sanitized SMILES already carries final validated valence.
            // AddHs expands topology and explicitly invalidates VALENCE in its
            // operation contract; prepare that new topology once, not again
            // for the unchanged parser-only workflow.
            let prepared = if add_hydrogens {
                with_requested_hydrogens
                    .with_assigned_valence()
                    .unwrap_or_else(|error| panic!("post-AddHs valence for {smiles:?}: {error}"))
            } else {
                with_requested_hydrogens
            };

            let atom_state_before = prepared
                .atoms()
                .iter()
                .map(|atom| {
                    (
                        atom.id(),
                        atom.hybridization(),
                        atom.formal_charge(),
                        atom.is_aromatic(),
                    )
                })
                .collect::<Vec<_>>();
            let bond_state_before = prepared
                .bonds()
                .iter()
                .map(|bond| (bond.id(), bond.is_aromatic(), bond.is_conjugated()))
                .collect::<Vec<_>>();

            let actual = prepared
                .uff_has_all_molecule_params()
                .unwrap_or_else(|error| panic!("UFF query for {smiles:?}: {error}"));
            query_calls += 1;
            assert_eq!(actual, expected, "input={smiles:?}, add_h={add_hydrogens}");

            if matches!(smiles, "*" | "[Cu+]") {
                assert!(
                    !actual,
                    "missing UFF parameter is a successful false result"
                );
            }
            if smiles.is_empty() {
                assert!(actual, "empty source atom traversal returns true");
            }

            assert_eq!(
                prepared
                    .atoms()
                    .iter()
                    .map(|atom| (
                        atom.id(),
                        atom.hybridization(),
                        atom.formal_charge(),
                        atom.is_aromatic(),
                    ))
                    .collect::<Vec<_>>(),
                atom_state_before,
                "query changed source-derived atom typing state for {smiles:?}"
            );
            assert_eq!(
                prepared
                    .bonds()
                    .iter()
                    .map(|bond| (bond.id(), bond.is_aromatic(), bond.is_conjugated()))
                    .collect::<Vec<_>>(),
                bond_state_before,
                "query changed source-derived bond typing state for {smiles:?}"
            );
        }
    }

    assert_eq!(
        query_calls, 20,
        "one actual query call per fixed reference row"
    );
}
