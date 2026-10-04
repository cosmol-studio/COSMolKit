//! Development-only, one-shot preparation against the frozen published library.
//! Compile with rustc and the isolated baseline rlibs; never called by tests.
//! All outputs use create_new: a mismatch cannot refresh a reference.
use std::{error::Error, fs::{self, OpenOptions}, io::Write, path::Path};
use cosmolkit_core::{AtomId, AtomSpec, BondDirection, BondOrder, BondSpec,
    ChiralTag, Element, Molecule, MoleculeBuilder, ValenceModel};
use serde_json::json;

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), Box<dyn Error>> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    file.write_all(bytes)?;
    Ok(())
}

fn ring_atoms(nitrogen: bool) -> Vec<AtomSpec> {
    (0..6).map(|i| AtomSpec::new(if nitrogen && i == 0 { Element::N } else { Element::C })
        .with_aromatic(true)).collect()
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 { return Err("usage: driver FRESH_OUTPUT_DIR PINNED_FONT".into()); }
    let out = Path::new(&args[1]);
    if out.exists() { return Err("output directory already exists; refusing regeneration".into()); }
    let c = || AtomSpec::new(Element::C);
    let single = BondOrder::Single;
    let double = BondOrder::Double;
    let ring_edges: Vec<_> = (0..6).map(|i| (i, (i + 1) % 6, BondOrder::Aromatic)).collect();
    let ring_coords = vec![[1.5, 0.0], [0.75, 1.299038105676658],
        [-0.75, 1.299038105676658], [-1.5, 0.0],
        [-0.75, -1.299038105676658], [0.75, -1.299038105676658]];
    let cases: Vec<(&str, Vec<AtomSpec>, Vec<(usize, usize, BondOrder)>, Vec<[f64; 2]>)> = vec![
        ("empty", vec![], vec![], vec![]),
        ("single_carbon", vec![c()], vec![], vec![[0.0, 0.0]]),
        ("explicit_h_methane", vec![c().with_no_implicit(true),
            AtomSpec::new(Element::H), AtomSpec::new(Element::H),
            AtomSpec::new(Element::H), AtomSpec::new(Element::H)],
            vec![(0,1,single),(0,2,single),(0,3,single),(0,4,single)],
            vec![[0.0,0.0],[1.5,0.0],[0.0,1.5],[-1.5,0.0],[0.0,-1.5]]),
        ("ethane", vec![c(),c()], vec![(0,1,single)], vec![[0.0,0.0],[1.5,0.0]]),
        ("ethanol", vec![c(),c(),AtomSpec::new(Element::O)],
            vec![(0,1,single),(1,2,single)], vec![[0.0,0.0],[1.5,0.0],[2.25,1.299038105676658]]),
        ("benzene", ring_atoms(false), ring_edges.clone(), ring_coords.clone()),
        ("pyridine", ring_atoms(true), ring_edges, ring_coords),
        ("carbonyl", vec![c(),AtomSpec::new(Element::O)], vec![(0,1,double)], vec![[0.0,0.0],[1.5,0.0]]),
        ("isotope", vec![c().with_isotope(13),AtomSpec::new(Element::O)],
            vec![(0,1,single)], vec![[0.0,0.0],[1.5,0.0]]),
        ("formal_charge", vec![AtomSpec::new(Element::N).with_formal_charge(1).with_explicit_hydrogens(4)
            .with_no_implicit(true)], vec![], vec![[0.0,0.0]]),
        ("tetrahedral_wedge", vec![c().with_chiral_tag(ChiralTag::TetrahedralCcw),
            AtomSpec::new(Element::F),AtomSpec::new(Element::CL),
            AtomSpec::new(Element::BR),AtomSpec::new(Element::I)],
            vec![(0,1,single),(0,2,single),(0,3,single),(0,4,single)],
            vec![[0.0,0.0],[1.5,0.0],[0.0,1.5],[-1.5,0.0],[0.0,-1.5]]),
        ("mapped_atoms", vec![c().with_atom_map(7),AtomSpec::new(Element::O).with_atom_map(23)],
            vec![(0,1,single)], vec![[0.0,0.0],[1.5,0.0]]),
    ];
    // Evaluate every baseline output before publishing any file. A baseline
    // error aborts the preparation and never becomes an empty fixture.
    let mut artifacts: Vec<(String, Vec<u8>)> = Vec::new();
    for (label, atoms, bonds, coords) in cases {
        let atom_specs: Vec<_> = atoms.iter().map(|a| format!("{a:#?}")).collect();
        let mut builder = MoleculeBuilder::new();
        for atom in atoms { builder.add_atom(atom); }
        for &(begin,end,order) in &bonds {
            let direction = if label == "tetrahedral_wedge" && end == 1 {
                BondDirection::BeginWedge
            } else { BondDirection::None };
            builder.add_bond(BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
                .with_direction(direction))?;
        }
        builder.set_2d_coordinates(coords.clone())?;
        let molecule: Molecule = builder.build()?;
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
        let svg = molecule.to_svg(300,300)?;
        let png = molecule.to_png(300,300)?;
        assert_eq!(before,format!("{molecule:#?}"),"baseline source changed: {label}");
        let pixmap = tiny_skia::Pixmap::decode_png(&png)?;
        assert_eq!((pixmap.width(),pixmap.height()),(300,300));
        artifacts.push((format!("{label}.input.json"),serde_json::to_vec_pretty(&input)?));
        artifacts.push((format!("{label}.svg"),svg.into_bytes()));
        artifacts.push((format!("{label}.png"),png));
        artifacts.push((format!("{label}.rgba"),pixmap.data().to_vec()));
        println!("prepared {label}: 300x300");
    }
    // Exact original inline font-regression input, independent of chemistry.
    let svg = "<?xml version='1.0' encoding='iso-8859-1'?>\
                   <svg version='1.1' baseProfile='full' xmlns='http://www.w3.org/2000/svg' \
                   width='120px' height='80px' viewBox='0 0 120 80'>\
                   <rect width='120' height='80' fill='#FFFFFF'/>\
                   <text x='12' y='54' style='font-size:48px;font-style:normal;\
                   font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;\
                   text-anchor:start;fill:#000000'>O</text>\
                   </svg>";
    let mut options = usvg::Options::default();
    options.font_family = "Noto Sans".to_owned();
    options.fontdb_mut().load_font_source(usvg::fontdb::Source::Binary(
        std::sync::Arc::new(fs::read(&args[2])?)));
    options.fontdb_mut().set_sans_serif_family("Noto Sans");
    let tree = usvg::Tree::from_str(svg,&options)?;
    let size = tree.size().to_int_size();
    let mut pixmap = tiny_skia::Pixmap::new(size.width(),size.height()).ok_or("font pixmap allocation")?;
    resvg::render(&tree,tiny_skia::Transform::default(),&mut pixmap.as_mut());
    assert_eq!((pixmap.width(),pixmap.height()),(120,80));
    assert!(pixmap.pixels().iter().any(|p|p.alpha()==255 && p.red()<245 && p.green()<245 && p.blue()<245));
    artifacts.push(("font_text.svg".into(),svg.as_bytes().to_vec()));
    artifacts.push(("font_text.png".into(),pixmap.encode_png()?));
    artifacts.push(("font_text.rgba".into(),pixmap.data().to_vec()));
    fs::create_dir_all(out)?;
    for (name,bytes) in artifacts { write_new(&out.join(name),&bytes)?; }
    println!("published 12 supplied-layout cases + exact legacy font case; no migrated library linked");
    Ok(())
}
