#![cfg(all(feature = "cap-io", feature = "cap-smiles"))]
//! Exact original CW/CCW builders, header and code1/6 assertions; explicit force_2d transport.
use cosmolkit::*;

#[test]
fn mol_to_v2000_block_infers_wedge_from_chiral_tag_without_coordinates() {
    // Build a chiral center with 4 explicit non-H neighbors, NO coordinates.
    // pick_bonds_to_wedge should select the bond to the lowest-scored neighbor
    // (lowest atomic number wins when degree and chirality are equal).
    let mut builder = MoleculeBuilder::new().with_name("no-coords-wedge".into());
    let c_chiral =
        builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
    let c_me = builder.add_atom(AtomSpec::new(Element::C));
    let n = builder.add_atom(AtomSpec::new(Element::N));
    let o = builder.add_atom(AtomSpec::new(Element::O));
    let f = builder.add_atom(AtomSpec::new(Element::F));
    builder
        .add_bond(BondSpec::new(c_chiral, c_me, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, n, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, o, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, f, BondOrder::Single))
        .unwrap();
    // Deliberately no coordinates
    let molecule = builder.build().unwrap();
    assert!(molecule.coordinates_2d().is_none());

    let block = canonical_mol(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            include_stereo: true,
            kekulize: false,
            force_2d: true,
            ..Default::default()
        },
    )
    .unwrap();

    println!("generated_original_complete_block={block:?}\n{block}");
    assert!(block.starts_with("no-coords-wedge\n  COSMolKit          2D\n\n"));
    // The C-methyl bond (atom 1 → 2, added first) has the lowest atomic
    // number (6 vs 7/8/9) so pick_bond_to_wedge scores it lowest.
    // TetrahedralCw → BeginWedge → stereo code 1.
    // V2000 bond format: "{:>3}{:>3}{:>3} {:>2}" → stereo at [10..12].
    let has_wedge = block
        .lines()
        .any(|l| l.starts_with("  1  2") && l.len() >= 12 && l[10..12].trim() == "1");
    assert!(
        has_wedge,
        "expected BeginWedge (code 1) on bond 1-2 in output:\n{block}"
    );
}

#[test]
fn mol_to_v2000_block_infers_dash_from_chiral_tag_ccw_without_coordinates() {
    // Same molecule as above but with TetrahedralCcw → BeginDash → code 6.
    let mut builder = MoleculeBuilder::new().with_name("no-coords-dash".into());
    let c_chiral =
        builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw));
    let c_me = builder.add_atom(AtomSpec::new(Element::C));
    let n = builder.add_atom(AtomSpec::new(Element::N));
    let o = builder.add_atom(AtomSpec::new(Element::O));
    let f = builder.add_atom(AtomSpec::new(Element::F));
    builder
        .add_bond(BondSpec::new(c_chiral, c_me, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, n, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, o, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, f, BondOrder::Single))
        .unwrap();
    let molecule = builder.build().unwrap();
    assert!(molecule.coordinates_2d().is_none());

    let block = canonical_mol(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            include_stereo: true,
            kekulize: false,
            force_2d: true,
            ..Default::default()
        },
    )
    .unwrap();

    println!("generated_original_complete_block={block:?}\n{block}");
    assert!(block.starts_with("no-coords-dash\n  COSMolKit          2D\n\n"));
    let has_dash = block
        .lines()
        .any(|l| l.starts_with("  1  2") && l.len() >= 12 && l[10..12].trim() == "6");
    assert!(
        has_dash,
        "expected BeginDash (code 6) on bond 1-2 in output:\n{block}"
    );
}

fn canonical_mol(m: &Molecule, p: &MolBlockWriteParams) -> Result<String, MolecularIoError> {
    let generated = m
        .with_2d_coordinates_with_params(&Coordinate2DParams {
            canonical_orientation: true,
            ..Default::default()
        })
        .map_err(MolecularIoError::Construction)?;
    assert!(generated.coordinates_2d().is_some());
    assert!(m.coordinates_2d().is_none());
    generated.to_mol_with_params(&MolBlockWriteParams {
        coordinate_selection: MolCoordinateSelection::TwoD { id: 0 },
        ..*p
    })
}
