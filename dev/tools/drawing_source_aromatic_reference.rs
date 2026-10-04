//! One-shot independent published oracle with two pinned source repairs.
//! rustc links only the fresh isolated release rlibs; tests never run this.
use cosmolkit_core::{AtomId, AtomSpec, BondOrder, BondSpec, Element, Molecule,
    MoleculeBuilder, ValenceModel};
use serde_json::json;
use std::{error::Error, fs::{self, OpenOptions}, io::Write, path::Path,
    process::Command};

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), Box<dyn Error>> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    file.write_all(bytes)?;
    Ok(())
}

fn ring_atoms(nitrogen: bool) -> Vec<AtomSpec> {
    (0..6).map(|i| AtomSpec::new(if nitrogen && i == 0 { Element::N } else { Element::C })
        .with_aromatic(true)).collect()
}

fn check_state(molecule: &Molecule, nitrogen: bool, coords: &[[f64; 2]]) {
    assert_eq!((molecule.num_atoms(), molecule.num_bonds()), (6, 6));
    for (i, atom) in molecule.atoms().iter().enumerate() {
        assert_eq!(atom.id().index(), i);
        assert_eq!(atom.atomic_number(), if nitrogen && i == 0 { 7 } else { 6 });
        assert!(atom.is_aromatic());
    }
    for (i, bond) in molecule.bonds().iter().enumerate() {
        assert_eq!(bond.id().index(), i);
        assert_eq!((bond.begin().index(), bond.end().index()), (i, (i + 1) % 6));
        assert_eq!(bond.order(), BondOrder::Aromatic);
        assert!(!bond.is_aromatic());
    }
    assert_eq!(molecule.conformers_2d().len(), 1);
    assert_eq!(molecule.conformers_2d()[0].id(), 0);
    assert!(molecule.conformers_3d().is_empty());
    let bits = |points: &[[f64; 2]]| points.iter()
        .map(|p| [p[0].to_bits(), p[1].to_bits()]).collect::<Vec<_>>();
    assert_eq!(bits(molecule.coordinates_2d().unwrap()), bits(coords));
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 { return Err("usage: driver FRESH_OUTPUT_DIR OLD_FIXTURE_DIR FROZEN_IDENTITY".into()); }
    let out = Path::new(&args[1]);
    if out.exists() { return Err("output directory exists; refusing regeneration".into()); }
    let old = Path::new(&args[2]);
    let identity = fs::read(&args[3])?;
    let _: serde_json::Value = serde_json::from_slice(&identity)?;
    let ring_edges: Vec<_> = (0..6).map(|i| (i, (i + 1) % 6, BondOrder::Aromatic)).collect();
    let ring_coords = vec![[1.5, 0.0], [0.75, 1.299038105676658],
        [-0.75, 1.299038105676658], [-1.5, 0.0],
        [-0.75, -1.299038105676658], [0.75, -1.299038105676658]];
    let cases = vec![
        ("benzene", ring_atoms(false), ring_edges.clone(), ring_coords.clone()),
        ("pyridine", ring_atoms(true), ring_edges, ring_coords),
    ];
    let mut artifacts: Vec<(String, Vec<u8>)> = Vec::new();
    for (label, atoms, bonds, coords) in cases {
        let atom_specs: Vec<_> = atoms.iter().map(|a| format!("{a:#?}")).collect();
        let mut builder = MoleculeBuilder::new();
        for atom in atoms { builder.add_atom(atom); }
        for &(begin,end,order) in &bonds {
            builder.add_bond(BondSpec::new(AtomId::new(begin), AtomId::new(end), order))?;
        }
        builder.set_2d_coordinates(coords.clone())?;
        let molecule: Molecule = builder.build()?;
        let nitrogen = label == "pyridine";
        check_state(&molecule, nitrogen, &coords);
        let before = format!("{molecule:#?}");
        let valence = cosmolkit_core::valence::assign_valence(&molecule, ValenceModel::RdkitLike)?;
        let rings = cosmolkit_core::rings::find_sssr(&molecule)?;
        let input = json!({"label":label,"width":300,"height":300,
            "atom_specs":atom_specs,"atoms":molecule.atoms().iter().map(|a| format!("{a:#?}")).collect::<Vec<_>>(),
            "bonds":molecule.bonds().iter().map(|b| format!("{b:#?}")).collect::<Vec<_>>(),
            "coordinates":coords,"coordinate_bits":coords.iter().map(|p|[p[0].to_bits(),p[1].to_bits()]).collect::<Vec<_>>(),
            "properties":format!("{:#?}",molecule.properties()),
            "source_snapshot":before,"computed_valence":format!("{valence:#?}"),
            "computed_rings":format!("{rings:#?}")});
        let input_bytes = serde_json::to_vec_pretty(&input)?;
        assert_eq!(input_bytes, fs::read(old.join(format!("{label}.input.json")))?, "raw input identity: {label}");
        // Exposed preparation counterpart on a separate copy, not private drawer state.
        let prepared = molecule.clone().with_kekulized_bonds(false)?;
        check_state(&prepared, nitrogen, &coords);
        assert_eq!(format!("{:#?}", prepared.atoms()), format!("{:#?}", molecule.atoms()));
        assert_eq!(format!("{:#?}", prepared.bonds()), format!("{:#?}", molecule.bonds()));
        assert_eq!(format!("{:#?}", prepared.properties()), format!("{:#?}", molecule.properties()));
        assert_eq!(format!("{:#?}", prepared.conformers_2d()), format!("{:#?}", molecule.conformers_2d()));
        assert_eq!(before, format!("{molecule:#?}"), "B1 changed raw source: {label}");
        println!("B0/B1 {label}: A6 flagFalse6 atomTrue6 IDs/endpoints/elements/XYbits/ID0 unchanged; exposed preparation counterpart");
        let svg_result = molecule.to_svg(300,300);
        assert_eq!(before, format!("{molecule:#?}"), "to_svg input changed: {label}");
        check_state(&molecule, nitrogen, &coords);
        let svg = svg_result?;
        let png_result = molecule.to_png(300,300);
        assert_eq!(before, format!("{molecule:#?}"), "to_png input changed: {label}");
        check_state(&molecule, nitrogen, &coords);
        let png = png_result?;
        let pixmap = tiny_skia::Pixmap::decode_png(&png)?;
        assert_eq!((pixmap.width(),pixmap.height()),(300,300));
        artifacts.push((format!("{label}.input.json"),input_bytes));
        artifacts.push((format!("{label}.svg"),svg.into_bytes()));
        artifacts.push((format!("{label}.png"),png));
        artifacts.push((format!("{label}.rgba"),pixmap.data().to_vec()));
        println!("prepared {label}: one public SVG and one public PNG, each source preserved, RGBA300x300");
    }
    artifacts.push(("identity.json".into(), identity));
    artifacts.sort_by(|a,b| a.0.cmp(&b.0));
    fs::create_dir(out)?;
    for (name, bytes) in &artifacts { write_new(&out.join(name), bytes)?; }
    let result = Command::new("sha256sum").current_dir(out)
        .args(artifacts.iter().map(|(name,_)| name)).output()?;
    if !result.status.success() { return Err("reference checksum creation failed".into()); }
    write_new(&out.join("SHA256SUMS"), &result.stdout)?;
    println!("frozen two independent scenes, identity and checksums; no current workspace linked or examined");
    Ok(())
}
