#![cfg(all(
    feature = "cap-depict",
    feature = "cap-smiles",
    feature = "cap-kekulize",
    feature = "cap-hydrogens"
))]

use cosmolkit::{BondOrder, Molecule};

#[test]
fn drawing_missing_2d_and_3d_only_queries_preserve_source_coordinates() {
    use cosmolkit::{Conformer3D, CoordinateBlock};
    let parsed = Molecule::from_smiles("CCO").unwrap();
    for three_d_only in [false, true] {
        let coordinates = CoordinateBlock {
            conformers_3d: if three_d_only {
                vec![
                    Conformer3D::new(
                        29,
                        vec![[0.0, -0.0, 0.0], [1.0, 0.5, 0.25], [2.0, 1.0, 1.25]],
                        true,
                    )
                    .with_prop("source", "independent-3d"),
                ]
            } else {
                Vec::new()
            },
            ..CoordinateBlock::default()
        };
        let source = Molecule::from_parts(
            parsed.topology().clone(),
            coordinates.clone(),
            parsed.properties().clone(),
        )
        .unwrap();
        // Freeze the actual validated live input before the first query.
        // Construction infers source_coordinate_dim from the stored dimensions.
        let coordinates = source.to_builder().coordinates().clone();
        let peer = source.clone();
        let topology = source.topology().clone();
        let properties = source.properties().clone();
        let bits = source
            .conformers_3d()
            .iter()
            .flat_map(|c| c.coordinates().iter().flatten().map(|v| v.to_bits()))
            .collect::<Vec<_>>();
        for molecule in [&source, &peer] {
            assert!(molecule.to_svg(300, 300).unwrap().contains("</svg>"));
            for original in [&source, &peer] {
                assert_eq!(original.topology(), &topology);
                assert_eq!(original.properties(), &properties);
                assert_eq!(original.to_builder().coordinates(), &coordinates);
                assert!(original.coordinates_2d().is_none());
                assert_eq!(
                    original
                        .conformers_3d()
                        .iter()
                        .flat_map(|c| c.coordinates().iter().flatten().map(|v| v.to_bits()))
                        .collect::<Vec<_>>(),
                    bits
                );
            }
            assert_eq!(
                &molecule.to_png(300, 300).unwrap()[..8],
                b"\x89PNG\r\n\x1a\n"
            );
            for original in [&source, &peer] {
                assert_eq!(original.topology(), &topology);
                assert_eq!(original.properties(), &properties);
                assert_eq!(original.to_builder().coordinates(), &coordinates);
                assert!(original.coordinates_2d().is_none());
                assert_eq!(
                    original
                        .conformers_3d()
                        .iter()
                        .flat_map(|c| c.coordinates().iter().flatten().map(|v| v.to_bits()))
                        .collect::<Vec<_>>(),
                    bits
                );
            }
        }
    }
}

// Preparation observations below exercise existing public owner routes.
// Direct detached drawing-preparation assertions remain in drawing_prepare_tests.rs.
fn prepare_molecule_for_drawing(
    molecule: &Molecule,
    force_coords: bool,
) -> Result<Molecule, cosmolkit::OperationError> {
    let mut prepared = molecule.with_kekulized_bonds_with_params(&cosmolkit::KekulizeParams {
        mark_atoms_bonds: false,
        ..cosmolkit::KekulizeParams::default()
    })?;
    let rings = cosmolkit_core::symmetrized_sssr(
        prepared.topology(),
        &cosmolkit_core::RingSearchParams::default(),
    )
    .unwrap();
    let candidates = prepared
        .atoms()
        .iter()
        .filter(|atom| {
            rings.num_atom_rings(atom.id()) > 1
                && matches!(
                    atom.chiral_tag(),
                    cosmolkit::ChiralTag::TetrahedralCw | cosmolkit::ChiralTag::TetrahedralCcw
                )
        })
        .map(|atom| atom.id())
        .collect::<Vec<_>>();
    if !candidates.is_empty() {
        prepared = prepared.with_hydrogens_with_params(&cosmolkit::AddHsParams {
            only_on_atoms: Some(candidates),
            ..cosmolkit::AddHsParams::default()
        })?;
    }
    if force_coords || prepared.coordinates_2d().is_none() {
        prepared = prepared.with_2d_coordinates_with_params(&cosmolkit::Coordinate2DParams {
            canonical_orientation: true,
            ..cosmolkit::Coordinate2DParams::default()
        })?;
    }
    Ok(prepared)
}

#[test]
fn aromatic_outer_line_draw_coordinate_matches_rdkit_rounding_regression() {
    let molecule = Molecule::from_smiles(r"O=C(O)/C(=C\c1ccccc1)Cc1ccccc1").unwrap();
    let svg = molecule.to_svg(300, 300).unwrap();

    assert!(svg.contains("<path class='bond-12 atom-12 atom-13' d='M 76.7,130.0 L 48.7,113.9'"));
}

#[test]
fn overlapping_depiction_atoms_skip_bond_regression() {
    let molecule = Molecule::from_smiles("CC12CC34CN(C(=N)N)CC3(C1)CC2(C)C4.Cl").unwrap();
    let svg = molecule.to_svg(300, 300).unwrap();

    assert!(!svg.contains("<path class='bond-16 atom-13 atom-1'"));
    assert!(svg.contains("<path class='bond-4 atom-4 atom-5' d='M 121.8,184.1 L 111.8,170.9'"));
}

#[test]
fn minimum_svg_font_size_clamps_at_binary64_boundary_regression() {
    let molecule = Molecule::from_smiles(
        "C=C(c1ccc(Oc2ccc(Oc3ccc(C(=C)C4COC5(CCCCC5)OO4)cc3)cc2)cc1)C1COC2(CCCCC2)OO1",
    )
    .unwrap();
    let svg = molecule.to_svg(300, 300).unwrap();

    assert!(svg.contains("class='atom-6' style='font-size:6px"));
    assert!(!svg.contains("style='font-size:5px"));
}

#[test]
fn prepared_drawing_uses_rdkit_canonical_kekule_assignment_regression() {
    let molecule = Molecule::from_smiles(
        "CC1=C(CCNCCN)c2cc3[nH]c(cc4[nH]c(cc5nc(cc1n2)C(C)=C5CCNCCN)c(C)c4CCNCCN)c(CCNCCN)c3C",
    )
    .unwrap();
    let prepared = prepare_molecule_for_drawing(&molecule, true).unwrap();
    let bond_orders: String = prepared
        .bonds()
        .iter()
        .map(|bond| match bond.order() {
            BondOrder::Single => 'S',
            BondOrder::Double => 'D',
            BondOrder::Triple => 'T',
            BondOrder::Aromatic => 'A',
            _ => '?',
        })
        .collect();

    assert_eq!(
        bond_orders,
        "SDSSSSSSSSDSSDSSSSDSDSDSSSDSSSSSSDSSSSSSSSSSSSSSSDSSDSDS"
    );
    let coordinates = prepared.coordinates_2d().unwrap();
    assert_eq!(
        coordinates[0][0].to_bits(),
        (-0.7271872994113224_f64).to_bits()
    );
    assert_eq!(
        coordinates[0][1].to_bits(),
        (-7.508043887184295_f64).to_bits()
    );
    assert_eq!(
        coordinates[49][0].to_bits(),
        10.57352991454218_f64.to_bits()
    );
    assert_eq!(
        coordinates[49][1].to_bits(),
        4.749281830440951_f64.to_bits()
    );

    let svg = molecule.to_svg(300, 300).unwrap();
    assert!(svg.contains("<path class='bond-0 atom-0 atom-1' d='M 140.4,250.0 L 145.7,232.9'"));
}

#[test]
fn polycyclic_chiral_hydrogen_preparation_matches_rdkit_regression() {
    let molecule = Molecule::from_smiles(
        "O=C1CCC(=O)O[C@@H]2[C@@H](O)[C@H](O)[C@@H](COC(=O)CCC(=O)O[C@@H]3[C@@H](O)[C@H](O)[C@@H](CO1)O[C@@H]3O)O[C@@H]2O",
    )
    .unwrap();
    let prepared = prepare_molecule_for_drawing(&molecule, false).unwrap();

    assert_eq!(prepared.num_atoms(), 46);
    assert_eq!(
        prepared
            .atoms()
            .iter()
            .filter(|atom| atom.atomic_number() == 1)
            .count(),
        10
    );

    let svg = molecule.to_svg(300, 300).unwrap();
    assert!(svg.contains("<path class='bond-0 atom-0 atom-1' d='M 189.8,250.4 L 186.4,239.5'"));
}

#[test]
fn test_mol_to_svg_contains_expected_elements() {
    // Build a simple methane molecule
    use cosmolkit::BondSpec;
    use cosmolkit::MoleculeBuilder;
    use cosmolkit::{AtomSpec, Element};

    let mut builder = MoleculeBuilder::new();
    let c = builder.add_atom(AtomSpec::new(Element::C));
    let h1 = builder.add_atom(AtomSpec::new(Element::H));
    let h2 = builder.add_atom(AtomSpec::new(Element::H));
    let h3 = builder.add_atom(AtomSpec::new(Element::H));
    let h4 = builder.add_atom(AtomSpec::new(Element::H));
    builder
        .add_bond(BondSpec::new(c, h1, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c, h2, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c, h3, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c, h4, BondOrder::Single))
        .unwrap();
    let mol = builder.build().expect("build methane");

    let svg = mol.to_svg(300, 300).expect("svg rendering");
    // SVG doc should have expected tags
    assert!(
        svg.starts_with("<?xml"),
        "SVG should start with XML declaration"
    );
    assert!(svg.contains("viewBox"), "SVG should have viewBox");
    assert!(svg.contains("</svg>"), "SVG should have closing tag");
    // Atom text emission is depiction-policy dependent; for this smoke
    // test we only require a valid SVG document with molecular geometry.
    assert!(
        svg.contains("<line") || svg.contains("<path"),
        "SVG should contain bond geometry"
    );
}

#[test]
fn draw_svg_generates_missing_2d_coords_via_registered_operation() {
    use cosmolkit::BondSpec;
    use cosmolkit::MoleculeBuilder;
    use cosmolkit::{AtomSpec, Element};

    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
        .unwrap();
    let mol = builder.build().unwrap();

    let svg = mol.to_svg(300, 300).unwrap();

    assert!(svg.contains("</svg>"));
    assert!(svg.contains("<line") || svg.contains("<path"));
}

#[test]
fn svg_drawer_handles_dense_mapped_polycyclic_regression() {
    let smiles = "[C:12]12([CH:62]([CH3:65])[c:61]3[cH:64][cH:67][cH:68][cH:66][cH:63]3)[CH:20]4[c:30]5[c:40]6[c:49]7[c:57]8[c:60]([c:59]9[c:55]([c:47]([c:44]([c:52]9[c:51]([c:43]%10[c:35]%11[c:25]%12[c:19]%13%14)[c:53]8[c:45]%11[c:39]6[c:29]4%13)[c:34]([c:24]%15[c:15]%16[c:7]%17[c:3]%18%19)[c:33]%10[c:23]%16[c:16]%12[c:8]%18[c:11]%14[c:5]1%20)[c:37]([c:36]%21[c:26]%22[c:18]%23[c:10]%24[c:13]%25[c:6]%26%27)[c:27]%15[c:17]%22[c:9]%17[c:4]%24[c:1]%19[c:2]%20%26)[c:54]([c:46]%21[c:38]%28[c:28]%23[c:21]%25%29)[c:56]%30[c:48]%28[c:41]%31[c:31]%29[c:22]%32[c:14]2%27)[c:58]%30[c:50]7[c:42]%31[c:32]5%32";
    let molecule = Molecule::from_smiles(smiles).expect("parse dense mapped polycyclic SMILES");
    let svg = molecule
        .to_svg(300, 300)
        .expect("draw dense mapped polycyclic SVG");

    assert!(svg.contains("<svg"));
    assert!(svg.contains("width='300px'"));
    assert!(svg.contains("height='300px'"));
    assert!(!svg.contains("NaN"));
    assert!(!svg.contains("inf"));
}
