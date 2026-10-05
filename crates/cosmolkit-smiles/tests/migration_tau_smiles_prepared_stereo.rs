use cosmolkit_smiles::{parse_smiles, write_smiles};
use cosmolkit_types::{BondDirection, BondStereo};
#[test]
fn completed_stereo_marker_does_not_reperceive_a_new_double_bond_from_stale_slashes() {
    // Pinned source TAU case smiles_5000:1399, transform 1,3 (thio)keto/enol f,
    // mapping 12,13,14: native post_stereo -> ordinary MolToSmiles checkpoint.
    let mut record = parse_smiles(
        "O=C(O)CCCCCCC/C=C\\C=C(O)/C=C(O)/C=C/CCO",
        &Default::default(),
    )
    .unwrap();
    record =
        cosmolkit_smiles::finalize_smiles_stereo(record, &Default::default(), &mut None, &mut None)
            .unwrap();
    let bond = &mut record.topology.bonds[12];
    assert_eq!(bond.stereo(), BondStereo::E);
    bond.set_stereo(BondStereo::None).unwrap();
    bond.set_stereo_atoms(None);
    record
        .properties
        .set_computed_prop("_StereochemDone", "1")
        .unwrap();
    assert!(record.topology.bonds.iter().any(|b| matches!(
        b.direction(),
        BondDirection::EndDownRight | BondDirection::EndUpRight
    )));
    let before = record.clone();
    assert_eq!(
        write_smiles(&record).unwrap(),
        "O=C(O)CCCCCCC/C=C\\C=C(O)/C=C(O)/C=C/CCO"
    );
    assert_eq!(record, before);
}
