use cosmolkit::{
    Atom, AtomId, AtomSpec, Conformer3D, CoordinateBlock, Element, Molecule, MoleculeProperties,
    PropertyValue, TopologyBlock,
};

#[test]
fn public_archive20_roundtrip_preserves_state_and_does_not_detach_input() {
    let topology = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("__computedProps", "opaque")
                .unwrap()
                .with_computed_prop("rank", PropertyValue::Int(12))
                .unwrap(),
        )],
        vec![],
        vec![],
        vec![],
    )
    .unwrap();
    let coordinates = CoordinateBlock {
        conformers_3d: vec![
            Conformer3D::new(17, vec![[-0.0, 1.25, -2.5]], false).with_prop("note", "kept\0exact"),
        ],
        ..Default::default()
    };
    let properties = MoleculeProperties::default()
        .with_prop("__computedProps", "opaque\0value")
        .unwrap()
        .with_computed_prop("mass", "12")
        .unwrap()
        .with_sdf_data_field("repeat", "first")
        .with_sdf_data_field("repeat", "second");
    let molecule = Molecule::from_parts(topology, coordinates, properties).unwrap();
    let peer = molecule.clone();
    let topology_before = molecule.topology() as *const _;
    let properties_before = molecule.properties() as *const _;
    let bytes = molecule.to_binary().unwrap();
    assert_eq!(&bytes[..12], b"COSMOL\0\0\x02\0\0\0");
    let restored = Molecule::from_binary(&bytes).unwrap();
    assert_eq!(restored.topology(), molecule.topology());
    assert_eq!(restored.properties(), molecule.properties());
    assert_eq!(restored.conformers_3d(), molecule.conformers_3d());
    assert_eq!(
        restored.conformers_3d()[0].coordinates()[0][0].to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(restored.to_binary().unwrap(), bytes);
    assert_eq!(peer.to_binary().unwrap(), bytes);
    assert_eq!(topology_before, molecule.topology() as *const _);
    assert_eq!(properties_before, molecule.properties() as *const _);
    assert!(std::ptr::eq(molecule.topology(), peer.topology()));
    assert!(std::ptr::eq(molecule.properties(), peer.properties()));
}

#[test]
fn public_archive20_rejects_corruption_without_changing_an_existing_molecule() {
    let molecule = Molecule::new();
    let baseline = molecule.to_binary().unwrap();
    for end in 0..baseline.len() {
        assert!(
            Molecule::from_binary(&baseline[..end]).is_err(),
            "cut {end}"
        );
    }
    let mut future = baseline.clone();
    future[8] = 99;
    assert!(Molecule::from_binary(&future).is_err());
    // Structurally valid archive payload, but claims a valid ring cache without
    // its payload. IO decoding is not a substitute for runtime cache validation.
    let mut inconsistent_cache = baseline.clone();
    let mut offset = 14;
    let mut checked = false;
    for _ in 0..u16::from_le_bytes([baseline[12], baseline[13]]) {
        let id = u16::from_le_bytes([baseline[offset], baseline[offset + 1]]);
        let length =
            u32::from_le_bytes(baseline[offset + 6..offset + 10].try_into().unwrap()) as usize;
        let payload = offset + 10;
        if id == 3 {
            // Derived schema 1: four fields, field 0, validity bits 0.
            assert_eq!(&baseline[payload..payload + 3], &[4, 0, 0]);
            inconsistent_cache[payload + 2] = 1;
            checked = true;
            break;
        }
        offset = payload + length;
    }
    assert!(checked);
    assert!(Molecule::from_binary(&inconsistent_cache).is_err());
    assert_eq!(molecule.to_binary().unwrap(), baseline);
}
