#![cfg(all(feature = "cap-bio", feature = "cap-serialization"))]
use cosmolkit::*;

const PDB: &str = "MODEL        1\nATOM      1  N  AALA A   7      11.104  13.207   9.900  0.50 20.00           N  \nATOM      2  CA AALA A   7      12.210  13.912  10.555  0.50 21.00           C  \nENDMDL\nMODEL        2\nATOM      1  CA  GLY B   9      15.790  12.920  10.910  1.00 22.00           C  \nENDMDL\nEND\n";

#[test]
fn protein_transfer_preserves_viewer_hierarchy_annotations_and_float_bits() {
    let mut parts = BioStructure::from_pdb(PDB).unwrap().into_parts();
    let point = [
        f64::from_bits(0x3ff0_0000_0000_0001),
        -0.0,
        f64::from_bits(0x7ff8_0000_0000_0123),
    ];
    parts.coordinates.positions_mut()[0] = point;
    parts.source_state.name = "protein viewer".into();
    parts.source_state.conect_map.insert(1, vec![2, 2]);
    parts
        .source_state
        .raw_remarks
        .push("not a viewer field".into());
    parts.metadata.authors.push("not a viewer field".into());
    let address = AtomAddress::new(
        PdbChainId::from_ascii(b"A").unwrap(),
        ResidueAddress::new(Some(7), None, b"", ResidueName::from_ascii(b"ALA").unwrap()).unwrap(),
        "CA",
        Some(b'A'),
    );
    parts.connections.push(BioConnection {
        partner1: address.clone(),
        partner2: address.clone(),
        kind: BioConnectionKind::Disulf,
        reported_distance: 2.125,
        ..Default::default()
    });
    parts.cispeps.push(BioCisPep {
        partner_c: address.clone(),
        partner_n: address.clone(),
        reported_angle: -0.0,
        ..Default::default()
    });
    parts.mod_residues.push(BioModRes {
        chain_name: address.chain_name(),
        res_id: address.residue(),
        parent_comp_id: "ALA".into(),
        mod_id: "display".into(),
        details: "kept".into(),
    });
    parts.helices.push(BioHelix {
        start: address.clone(),
        end: address.clone(),
        pdb_helix_class: BioHelixClass::RAlpha,
        length: 1,
    });
    parts.sheets.push(BioSheet {
        name: "S1".into(),
        strands: vec![BioStrand::new(
            address.clone(),
            address.clone(),
            address.clone(),
            address,
            -1,
            "strand".into(),
        )],
    });
    let source = BioStructure::from_parts(parts).unwrap().protein().unwrap();
    let before = source.as_bio_structure().clone().into_parts();
    let bytes = source.to_transfer_binary().unwrap();
    let restored = Protein::from_transfer_binary(&bytes).unwrap();
    let after = restored.as_bio_structure().clone().into_parts();
    assert_eq!(restored.num_models(), 2);
    assert_eq!(after.models, before.models);
    assert_eq!(after.chains, before.chains);
    assert_eq!(after.residues, before.residues);
    for (actual, expected) in after.atoms.iter().zip(&before.atoms) {
        assert_eq!(actual.residue_id(), expected.residue_id());
        assert_eq!(actual.name(), expected.name());
        assert_eq!(actual.element(), expected.element());
        assert_eq!(actual.altloc(), expected.altloc());
        assert_eq!(actual.occupancy().to_bits(), expected.occupancy().to_bits());
        assert_eq!(actual.b_iso().to_bits(), expected.b_iso().to_bits());
        assert_eq!(actual.source(), expected.source());
    }
    for (actual, expected) in after
        .coordinates
        .positions()
        .iter()
        .zip(before.coordinates.positions())
    {
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
    }
    assert_eq!(after.connections, before.connections);
    assert_eq!(after.cispeps, before.cispeps);
    assert_eq!(after.mod_residues, before.mod_residues);
    assert_eq!(after.helices, before.helices);
    assert_eq!(after.sheets, before.sheets);
    assert_eq!(after.source_state.name, before.source_state.name);
    assert_eq!(
        after.source_state.conect_map,
        before.source_state.conect_map
    );
    assert!(after.source_state.raw_remarks.is_empty());
    assert!(after.metadata.authors.is_empty());
    assert_eq!(source.to_transfer_binary().unwrap(), bytes);
    assert_eq!(restored.to_transfer_binary().unwrap(), bytes);
    drop(source);
    assert_eq!(restored.num_atoms(), 3);
    assert_eq!(
        Protein::from_pdb("END\n")
            .unwrap()
            .to_transfer_binary()
            .and_then(|data| Protein::from_transfer_binary(&data))
            .unwrap()
            .num_atoms(),
        0
    );
}

#[test]
fn protein_transfer_rejects_version_before_payload_and_wrong_or_damaged_formats() {
    let bytes = Protein::from_pdb(PDB)
        .unwrap()
        .to_transfer_binary()
        .unwrap();
    let mut mismatch = bytes.clone();
    mismatch[10] = b'9';
    mismatch.truncate(10 + env!("CARGO_PKG_VERSION").len());
    match Protein::from_transfer_binary(&mismatch).unwrap_err() {
        ProteinTransferError::VersionMismatch { expected, actual } => {
            assert_eq!(expected, env!("CARGO_PKG_VERSION"));
            assert_ne!(actual, expected);
        }
        other => panic!("unexpected error: {other:?}"),
    }
    for end in [0, 8, 10, bytes.len() - 1] {
        assert!(matches!(
            Protein::from_transfer_binary(&bytes[..end]),
            Err(ProteinTransferError::InvalidData { .. })
        ));
    }
    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(Protein::from_transfer_binary(&trailing).is_err());
    assert!(Molecule::from_transfer_binary(&bytes).is_err());
    let molecule = Molecule::new();
    assert!(Protein::from_transfer_binary(&molecule.to_transfer_binary().unwrap()).is_err());
}
