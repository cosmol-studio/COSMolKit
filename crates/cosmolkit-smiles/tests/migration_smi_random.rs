use cosmolkit_smiles::{
    RandomSmilesWriteParams, SmilesParseParams, SmilesWriteParams, parse_smiles,
    write_random_smiles_vector, write_smiles_with_random,
};

// Original source-text fixtures decode only at this observation boundary.
// Invalid UTF-8 fails; the complete byte payload is never substituted.
fn fixture_writer_text(text: cosmolkit_model::PropertyText) -> String {
    String::from_utf8(text.into_bytes()).expect("original writer fixture UTF-8 bytes")
}

#[test]
fn random_writer_and_vector_match_pinned_sequences_and_shared_stream_semantics() {
    let record = parse_smiles("CC(C)(O)F", &SmilesParseParams::default()).expect("parse fixture");
    let before = record.clone();
    let params = SmilesWriteParams {
        canonical: false,
        ..SmilesWriteParams::default()
    };
    let vector_params = RandomSmilesWriteParams::default();
    let random_vector =
        |record: &cosmolkit_smiles::SmilesRecord, count, seed, params: &RandomSmilesWriteParams| {
            write_random_smiles_vector(record, count, seed, params).expect("random vector writer")
        };
    let expected = |values: &[&str]| -> Vec<cosmolkit_model::PropertyText> {
        values
            .iter()
            .map(|value| cosmolkit_model::PropertyText::from(*value))
            .collect()
    };

    // RDKit 2026.03.1, seed 42, canonical=false and doRandom=true emits this
    // branch order. The five-atom tree consumes x1 for its root, x2..x5 in
    // dfsFindCycles and x6..x9 in dfsBuildStack; x10 remains in the shared core
    // stream after the complete writer returns.
    cosmolkit_core::with_rdkit_random_generator(42, |_| ());
    assert_eq!(
        write_smiles_with_random(&record, &params, true)
            .map(fixture_writer_text)
            .expect("random writer"),
        "CC(F)(C)O"
    );
    let continuation = cosmolkit_core::with_rdkit_random_generator(0, |rng| rng.next_u32());
    assert_eq!(continuation, 488_601_845);

    // Pinned SmilesWrite.cpp::MolToRandomSmilesVect, RDKit 2026.03.1.
    // A positive seed is applied before the loop; all rows draw from the same
    // stream and preserve source traversal order.
    let seed42 = expected(&[
        "CC(F)(C)O",
        "CC(C)(F)O",
        "C(O)(C)(C)F",
        "OC(F)(C)C",
        "C(C)(C)(O)F",
    ]);
    assert_eq!(random_vector(&record, 5, 42, &vector_params), seed42);
    assert_eq!(random_vector(&record, 5, 42, &vector_params), seed42);

    // A zero seed preserves state across separate outer writer calls.
    cosmolkit_core::with_rdkit_random_generator(42, |_| ());
    assert_eq!(
        write_smiles_with_random(&record, &params, true)
            .map(fixture_writer_text)
            .expect("random single writer"),
        "CC(F)(C)O"
    );
    assert_eq!(
        random_vector(&record, 4, 0, &vector_params),
        expected(&["CC(C)(F)O", "C(O)(C)(C)F", "OC(F)(C)C", "C(C)(C)(O)F"])
    );

    // A positive seed reseeds even when the requested vector is empty; the
    // following zero-seed call starts at the first seed-1 traversal.
    assert!(random_vector(&record, 0, 1, &vector_params).is_empty());
    assert_eq!(
        random_vector(&record, 3, 0, &vector_params),
        expected(&["C(C)(F)(O)C", "C(C)(C)(O)F", "OC(F)(C)C"])
    );

    // The loaded release oracle uses rdcast<int> as static_cast<int>. High-bit
    // seeds therefore do not reseed; INT_MAX is positive and the minstd_rand
    // modulus normalizes its engine seed to one.
    for high_bit_seed in [0x8000_0000, u32::MAX] {
        assert!(random_vector(&record, 0, 42, &vector_params).is_empty());
        assert_eq!(
            random_vector(&record, 4, high_bit_seed, &vector_params),
            expected(&["CC(F)(C)O", "CC(C)(F)O", "C(O)(C)(C)F", "OC(F)(C)C"])
        );
    }
    assert_eq!(
        random_vector(&record, 4, 0x7fff_ffff, &vector_params),
        expected(&["C(C)(F)(O)C", "C(C)(C)(O)F", "OC(F)(C)C", "CC(O)(C)F"])
    );

    // The vector API retains duplicate strings rather than deduplicating rows.
    let duplicate_record =
        parse_smiles("CC", &SmilesParseParams::default()).expect("duplicate fixture");
    assert_eq!(
        random_vector(&duplicate_record, 5, 42, &vector_params),
        expected(&["CC", "CC", "CC", "CC", "CC"])
    );

    // Each source vector option is forwarded to every noncanonical random
    // writer invocation, including isomeric suppression, Kekule output and
    // explicit bond/hydrogen formatting.
    let stereo_record =
        parse_smiles("[13CH3][C@H](F)Cl", &SmilesParseParams::default()).expect("stereo fixture");
    assert_eq!(
        random_vector(&stereo_record, 3, 42, &vector_params),
        expected(&[
            "F[C@H](Cl)[13CH3]",
            "F[C@H](Cl)[13CH3]",
            "F[C@H](Cl)[13CH3]",
        ])
    );
    let nonisomeric = RandomSmilesWriteParams {
        isomeric_smiles: false,
        ..RandomSmilesWriteParams::default()
    };
    assert_eq!(
        random_vector(&stereo_record, 3, 42, &nonisomeric),
        expected(&["FC(Cl)C", "FC(Cl)C", "FC(Cl)C"])
    );
    let aromatic_record =
        parse_smiles("c1ccccc1", &SmilesParseParams::default()).expect("aromatic fixture");
    let kekule = RandomSmilesWriteParams {
        kekule: true,
        ..RandomSmilesWriteParams::default()
    };
    assert_eq!(
        random_vector(&aromatic_record, 3, 42, &kekule),
        expected(&["C1=CC=CC=C1", "C1C=CC=CC=1", "C1C=CC=CC=1"])
    );
    let chain_record = parse_smiles("CCO", &SmilesParseParams::default()).expect("chain fixture");
    let explicit_bonds = RandomSmilesWriteParams {
        all_bonds_explicit: true,
        ..RandomSmilesWriteParams::default()
    };
    assert_eq!(
        random_vector(&chain_record, 3, 42, &explicit_bonds),
        expected(&["C-C-O", "O-C-C", "C-C-O"])
    );
    let explicit_hydrogens = RandomSmilesWriteParams {
        all_hydrogens_explicit: true,
        ..RandomSmilesWriteParams::default()
    };
    assert_eq!(
        random_vector(&chain_record, 3, 42, &explicit_hydrogens),
        expected(&["[CH3][CH2][OH]", "[OH][CH2][CH3]", "[CH3][CH2][OH]"])
    );

    // RDKit 2026.03.1 MolToRandomSmilesVect, seed 42, count 4, after
    // MolFromSmiles with sanitize=true/removeHs=true. Tuple order is
    // doIsomericSmiles, doKekule, allBondsExplicit, allHsExplicit.
    const OPTION_TUPLES: [(bool, bool, bool, bool); 16] = [
        (false, false, false, false),
        (true, false, false, false),
        (false, true, false, false),
        (true, true, false, false),
        (false, false, true, false),
        (true, false, true, false),
        (false, true, true, false),
        (true, true, true, false),
        (false, false, false, true),
        (true, false, false, true),
        (false, true, false, true),
        (true, true, false, true),
        (false, false, true, true),
        (true, false, true, true),
        (false, true, true, true),
        (true, true, true, true),
    ];
    const CONNECTED_EXPECTED: [[&str; 4]; 16] = [
        [
            "C(F)(C)c1ccccc1",
            "c1(ccccc1)C(F)C",
            "C(C)(F)c1ccccc1",
            "c1ccccc1C(F)C",
        ],
        [
            "[C@H](F)([13CH3])c1ccccc1",
            "c1(ccccc1)[C@@H](F)[13CH3]",
            "[C@@H]([13CH3])(F)c1ccccc1",
            "c1ccccc1[C@@H](F)[13CH3]",
        ],
        [
            "C(F)(C)C1C=CC=CC=1",
            "C1(C=CC=CC=1)C(F)C",
            "C(C)(F)C1C=CC=CC=1",
            "C1=CC=CC=C1C(F)C",
        ],
        [
            "[C@H](F)([13CH3])C1C=CC=CC=1",
            "C1(C=CC=CC=1)[C@@H](F)[13CH3]",
            "[C@@H]([13CH3])(F)C1C=CC=CC=1",
            "C1=CC=CC=C1[C@@H](F)[13CH3]",
        ],
        [
            "C(-F)(-C)-c1:c:c:c:c:c:1",
            "c1(:c:c:c:c:c:1)-C(-F)-C",
            "C(-C)(-F)-c1:c:c:c:c:c:1",
            "c1:c:c:c:c:c:1-C(-F)-C",
        ],
        [
            "[C@H](-F)(-[13CH3])-c1:c:c:c:c:c:1",
            "c1(:c:c:c:c:c:1)-[C@@H](-F)-[13CH3]",
            "[C@@H](-[13CH3])(-F)-c1:c:c:c:c:c:1",
            "c1:c:c:c:c:c:1-[C@@H](-F)-[13CH3]",
        ],
        [
            "C(-F)(-C)-C1-C=C-C=C-C=1",
            "C1(-C=C-C=C-C=1)-C(-F)-C",
            "C(-C)(-F)-C1-C=C-C=C-C=1",
            "C1=C-C=C-C=C-1-C(-F)-C",
        ],
        [
            "[C@H](-F)(-[13CH3])-C1-C=C-C=C-C=1",
            "C1(-C=C-C=C-C=1)-[C@@H](-F)-[13CH3]",
            "[C@@H](-[13CH3])(-F)-C1-C=C-C=C-C=1",
            "C1=C-C=C-C=C-1-[C@@H](-F)-[13CH3]",
        ],
        [
            "[CH]([F])([CH3])[c]1[cH][cH][cH][cH][cH]1",
            "[c]1([cH][cH][cH][cH][cH]1)[CH]([F])[CH3]",
            "[CH]([CH3])([F])[c]1[cH][cH][cH][cH][cH]1",
            "[cH]1[cH][cH][cH][cH][c]1[CH]([F])[CH3]",
        ],
        [
            "[C@H]([F])([13CH3])[c]1[cH][cH][cH][cH][cH]1",
            "[c]1([cH][cH][cH][cH][cH]1)[C@@H]([F])[13CH3]",
            "[C@@H]([13CH3])([F])[c]1[cH][cH][cH][cH][cH]1",
            "[cH]1[cH][cH][cH][cH][c]1[C@@H]([F])[13CH3]",
        ],
        [
            "[CH]([F])([CH3])[C]1[CH]=[CH][CH]=[CH][CH]=1",
            "[C]1([CH]=[CH][CH]=[CH][CH]=1)[CH]([F])[CH3]",
            "[CH]([CH3])([F])[C]1[CH]=[CH][CH]=[CH][CH]=1",
            "[CH]1=[CH][CH]=[CH][CH]=[C]1[CH]([F])[CH3]",
        ],
        [
            "[C@H]([F])([13CH3])[C]1[CH]=[CH][CH]=[CH][CH]=1",
            "[C]1([CH]=[CH][CH]=[CH][CH]=1)[C@@H]([F])[13CH3]",
            "[C@@H]([13CH3])([F])[C]1[CH]=[CH][CH]=[CH][CH]=1",
            "[CH]1=[CH][CH]=[CH][CH]=[C]1[C@@H]([F])[13CH3]",
        ],
        [
            "[CH](-[F])(-[CH3])-[c]1:[cH]:[cH]:[cH]:[cH]:[cH]:1",
            "[c]1(:[cH]:[cH]:[cH]:[cH]:[cH]:1)-[CH](-[F])-[CH3]",
            "[CH](-[CH3])(-[F])-[c]1:[cH]:[cH]:[cH]:[cH]:[cH]:1",
            "[cH]1:[cH]:[cH]:[cH]:[cH]:[c]:1-[CH](-[F])-[CH3]",
        ],
        [
            "[C@H](-[F])(-[13CH3])-[c]1:[cH]:[cH]:[cH]:[cH]:[cH]:1",
            "[c]1(:[cH]:[cH]:[cH]:[cH]:[cH]:1)-[C@@H](-[F])-[13CH3]",
            "[C@@H](-[13CH3])(-[F])-[c]1:[cH]:[cH]:[cH]:[cH]:[cH]:1",
            "[cH]1:[cH]:[cH]:[cH]:[cH]:[c]:1-[C@@H](-[F])-[13CH3]",
        ],
        [
            "[CH](-[F])(-[CH3])-[C]1-[CH]=[CH]-[CH]=[CH]-[CH]=1",
            "[C]1(-[CH]=[CH]-[CH]=[CH]-[CH]=1)-[CH](-[F])-[CH3]",
            "[CH](-[CH3])(-[F])-[C]1-[CH]=[CH]-[CH]=[CH]-[CH]=1",
            "[CH]1=[CH]-[CH]=[CH]-[CH]=[C]-1-[CH](-[F])-[CH3]",
        ],
        [
            "[C@H](-[F])(-[13CH3])-[C]1-[CH]=[CH]-[CH]=[CH]-[CH]=1",
            "[C]1(-[CH]=[CH]-[CH]=[CH]-[CH]=1)-[C@@H](-[F])-[13CH3]",
            "[C@@H](-[13CH3])(-[F])-[C]1-[CH]=[CH]-[CH]=[CH]-[CH]=1",
            "[CH]1=[CH]-[CH]=[CH]-[CH]=[C]-1-[C@@H](-[F])-[13CH3]",
        ],
    ];
    const DISCONNECTED_EXPECTED: [[&str; 4]; 16] = [
        [
            "c1ccccc1.C1CCCCC1.FC(C)Cl",
            "c1ccccc1.C1CCCCC1.ClC(C)F",
            "c1ccccc1.C1CCCCC1.C(F)(Cl)C",
            "c1ccccc1.C1CCCCC1.FC(Cl)C",
        ],
        [
            "c1ccccc1.C1CCCCC1.F[C@@H]([13CH3])Cl",
            "c1ccccc1.C1CCCCC1.Cl[C@H]([13CH3])F",
            "c1ccccc1.C1CCCCC1.[C@@H](F)(Cl)[13CH3]",
            "c1ccccc1.C1CCCCC1.F[C@H](Cl)[13CH3]",
        ],
        [
            "C1=CC=CC=C1.C1CCCCC1.FC(C)Cl",
            "C1=CC=CC=C1.C1CCCCC1.ClC(C)F",
            "C1C=CC=CC=1.C1CCCCC1.C(F)(Cl)C",
            "C1=CC=CC=C1.C1CCCCC1.FC(Cl)C",
        ],
        [
            "C1=CC=CC=C1.C1CCCCC1.F[C@@H]([13CH3])Cl",
            "C1=CC=CC=C1.C1CCCCC1.Cl[C@H]([13CH3])F",
            "C1C=CC=CC=1.C1CCCCC1.[C@@H](F)(Cl)[13CH3]",
            "C1=CC=CC=C1.C1CCCCC1.F[C@H](Cl)[13CH3]",
        ],
        [
            "c1:c:c:c:c:c:1.C1-C-C-C-C-C-1.F-C(-C)-Cl",
            "c1:c:c:c:c:c:1.C1-C-C-C-C-C-1.Cl-C(-C)-F",
            "c1:c:c:c:c:c:1.C1-C-C-C-C-C-1.C(-F)(-Cl)-C",
            "c1:c:c:c:c:c:1.C1-C-C-C-C-C-1.F-C(-Cl)-C",
        ],
        [
            "c1:c:c:c:c:c:1.C1-C-C-C-C-C-1.F-[C@@H](-[13CH3])-Cl",
            "c1:c:c:c:c:c:1.C1-C-C-C-C-C-1.Cl-[C@H](-[13CH3])-F",
            "c1:c:c:c:c:c:1.C1-C-C-C-C-C-1.[C@@H](-F)(-Cl)-[13CH3]",
            "c1:c:c:c:c:c:1.C1-C-C-C-C-C-1.F-[C@H](-Cl)-[13CH3]",
        ],
        [
            "C1=C-C=C-C=C-1.C1-C-C-C-C-C-1.F-C(-C)-Cl",
            "C1=C-C=C-C=C-1.C1-C-C-C-C-C-1.Cl-C(-C)-F",
            "C1-C=C-C=C-C=1.C1-C-C-C-C-C-1.C(-F)(-Cl)-C",
            "C1=C-C=C-C=C-1.C1-C-C-C-C-C-1.F-C(-Cl)-C",
        ],
        [
            "C1=C-C=C-C=C-1.C1-C-C-C-C-C-1.F-[C@@H](-[13CH3])-Cl",
            "C1=C-C=C-C=C-1.C1-C-C-C-C-C-1.Cl-[C@H](-[13CH3])-F",
            "C1-C=C-C=C-C=1.C1-C-C-C-C-C-1.[C@@H](-F)(-Cl)-[13CH3]",
            "C1=C-C=C-C=C-1.C1-C-C-C-C-C-1.F-[C@H](-Cl)-[13CH3]",
        ],
        [
            "[cH]1[cH][cH][cH][cH][cH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[F][CH]([CH3])[Cl]",
            "[cH]1[cH][cH][cH][cH][cH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[Cl][CH]([CH3])[F]",
            "[cH]1[cH][cH][cH][cH][cH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[CH]([F])([Cl])[CH3]",
            "[cH]1[cH][cH][cH][cH][cH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[F][CH]([Cl])[CH3]",
        ],
        [
            "[cH]1[cH][cH][cH][cH][cH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[F][C@@H]([13CH3])[Cl]",
            "[cH]1[cH][cH][cH][cH][cH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[Cl][C@H]([13CH3])[F]",
            "[cH]1[cH][cH][cH][cH][cH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[C@@H]([F])([Cl])[13CH3]",
            "[cH]1[cH][cH][cH][cH][cH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[F][C@H]([Cl])[13CH3]",
        ],
        [
            "[CH]1=[CH][CH]=[CH][CH]=[CH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[F][CH]([CH3])[Cl]",
            "[CH]1=[CH][CH]=[CH][CH]=[CH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[Cl][CH]([CH3])[F]",
            "[CH]1[CH]=[CH][CH]=[CH][CH]=1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[CH]([F])([Cl])[CH3]",
            "[CH]1=[CH][CH]=[CH][CH]=[CH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[F][CH]([Cl])[CH3]",
        ],
        [
            "[CH]1=[CH][CH]=[CH][CH]=[CH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[F][C@@H]([13CH3])[Cl]",
            "[CH]1=[CH][CH]=[CH][CH]=[CH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[Cl][C@H]([13CH3])[F]",
            "[CH]1[CH]=[CH][CH]=[CH][CH]=1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[C@@H]([F])([Cl])[13CH3]",
            "[CH]1=[CH][CH]=[CH][CH]=[CH]1.[CH2]1[CH2][CH2][CH2][CH2][CH2]1.[F][C@H]([Cl])[13CH3]",
        ],
        [
            "[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[F]-[CH](-[CH3])-[Cl]",
            "[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[Cl]-[CH](-[CH3])-[F]",
            "[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[CH](-[F])(-[Cl])-[CH3]",
            "[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[F]-[CH](-[Cl])-[CH3]",
        ],
        [
            "[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[F]-[C@@H](-[13CH3])-[Cl]",
            "[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[Cl]-[C@H](-[13CH3])-[F]",
            "[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[C@@H](-[F])(-[Cl])-[13CH3]",
            "[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[F]-[C@H](-[Cl])-[13CH3]",
        ],
        [
            "[CH]1=[CH]-[CH]=[CH]-[CH]=[CH]-1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[F]-[CH](-[CH3])-[Cl]",
            "[CH]1=[CH]-[CH]=[CH]-[CH]=[CH]-1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[Cl]-[CH](-[CH3])-[F]",
            "[CH]1-[CH]=[CH]-[CH]=[CH]-[CH]=1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[CH](-[F])(-[Cl])-[CH3]",
            "[CH]1=[CH]-[CH]=[CH]-[CH]=[CH]-1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[F]-[CH](-[Cl])-[CH3]",
        ],
        [
            "[CH]1=[CH]-[CH]=[CH]-[CH]=[CH]-1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[F]-[C@@H](-[13CH3])-[Cl]",
            "[CH]1=[CH]-[CH]=[CH]-[CH]=[CH]-1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[Cl]-[C@H](-[13CH3])-[F]",
            "[CH]1-[CH]=[CH]-[CH]=[CH]-[CH]=1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[C@@H](-[F])(-[Cl])-[13CH3]",
            "[CH]1=[CH]-[CH]=[CH]-[CH]=[CH]-1.[CH2]1-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-1.[F]-[C@H](-[Cl])-[13CH3]",
        ],
    ];

    for (fixture, expected_matrix) in [
        ("c1ccccc1[C@H]([13CH3])F", &CONNECTED_EXPECTED),
        (
            "c1ccccc1.C1CCCCC1.[13CH3][C@H](F)Cl",
            &DISCONNECTED_EXPECTED,
        ),
    ] {
        let fixture_record = parse_smiles(fixture, &SmilesParseParams::default())
            .expect("parse random matrix fixture");
        let fixture_before = fixture_record.clone();
        let distinct_tuples = OPTION_TUPLES
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(distinct_tuples.len(), 16, "fixture={fixture}");

        for (options, expected_rows) in OPTION_TUPLES.iter().copied().zip(expected_matrix) {
            let (isomeric_smiles, kekule, all_bonds_explicit, all_hydrogens_explicit) = options;
            let option_params = RandomSmilesWriteParams {
                isomeric_smiles,
                kekule,
                all_bonds_explicit,
                all_hydrogens_explicit,
            };
            assert_eq!(
                random_vector(&fixture_record, 4, 42, &option_params),
                expected(expected_rows),
                "fixture={fixture}, options={options:?}"
            );
        }
        assert_eq!(fixture_record, fixture_before, "fixture={fixture}");
    }

    assert_eq!(record, before);
}
