//! One-shot fixed highlight references from the immutable published0.3.0 rlibs.
//! Standalone rustc driver; ordinary tests never execute this file.
//! Destination must be absent and every diagnostic file uses create_new.
use cosmolkit_core::{
    AtomId, AtomSpec, BondId, BondOrder, BondSpec, Element, MoleculeBuilder,
    RingFindType, RingInfo, ValenceModel,
};
use serde_json::json;
use std::{
    error::Error,
    fs::{self, OpenOptions},
    io::Write,
    path::Path,
    process::Command,
};

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), Box<dyn Error>> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    file.write_all(bytes)?;
    Ok(())
}

fn sha256(path: &Path) -> Result<String, Box<dyn Error>> {
    let result = Command::new("sha256sum").arg(path).output()?;
    if !result.status.success() {
        return Err(format!("sha256sum failed for {}: {:?}", path.display(), result).into());
    }
    let text = String::from_utf8(result.stdout)?;
    Ok(text.split_whitespace().next().ok_or("missing SHA256")?.to_owned())
}

fn main() -> Result<(), Box<dyn Error>> {
    let output = std::env::args().nth(1).ok_or("usage: driver NEW_OUTPUT_DIRECTORY")?;
    let output = Path::new(&output);
    if output.exists() {
        return Err("destination exists; reference refresh forbidden".into());
    }
    let cases = [
        ("absent", None, None),
        ("atom_one", Some("1"), None),
        ("atom_true", Some("true"), None),
        ("bond_empty", None, Some("")),
        ("bond_zero", None, Some("0")),
        ("bond_one", None, Some("1")),
        ("bond_true", None, Some("true")),
        ("bond_case", None, Some("True")),
    ];
    let mut artifacts = Vec::new();
    let mut inputs = Vec::new();
    let mut calls = 0;
    for (label, atom_note, bond_note) in cases {
        let mut atom0 = AtomSpec::new(Element::C);
        if let Some(note) = atom_note {
            atom0 = atom0.with_prop("_highlight", note);
        }
        let atom1 = AtomSpec::new(Element::C);
        let mut bond = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
        if let Some(note) = bond_note {
            bond = bond.with_prop("_highlight", note);
        }
        let specs = json!({
            "atoms": [format!("{atom0:#?}"), format!("{atom1:#?}")],
            "bond": format!("{bond:#?}"),
        });
        let mut builder = MoleculeBuilder::new();
        assert_eq!(builder.add_atom(atom0), AtomId::new(0));
        assert_eq!(builder.add_atom(atom1), AtomId::new(1));
        assert_eq!(builder.add_bond(bond)?, BondId::new(0));
        assert_eq!(builder.degree(AtomId::new(0)), 1);
        assert_eq!(builder.degree(AtomId::new(1)), 1);
        assert_eq!(builder.neighbor_bonds(AtomId::new(0)), &[BondId::new(0)]);
        assert_eq!(builder.neighbor_bonds(AtomId::new(1)), &[BondId::new(0)]);
        builder.set_2d_coordinates(vec![[0.0, 0.0], [1.5, 0.0]])?;
        let molecule = builder.build()?;
        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.atoms()[0].prop("_highlight"), atom_note);
        assert_eq!(molecule.atoms()[1].prop("_highlight"), None);
        assert_eq!(molecule.bonds()[0].prop("_highlight"), bond_note);
        assert!(molecule.properties().props().is_empty());
        assert_eq!(molecule.atoms()[0].props().len(), usize::from(atom_note.is_some()));
        assert!(molecule.atoms()[1].props().is_empty());
        assert_eq!(molecule.bonds()[0].props().len(), usize::from(bond_note.is_some()));
        assert!(molecule.prop("_highlight").is_none());
        assert!(molecule.atoms().iter().all(|a| a.computed_prop_names().is_empty()));
        assert!(molecule.bonds().iter().all(|b| b.computed_prop_names().is_empty()));
        assert!(molecule.properties().computed_prop_names().is_empty());
        assert_eq!(molecule.coordinates_2d().unwrap(), &[[0.0, 0.0], [1.5, 0.0]]);
        assert_eq!(molecule.conformers_2d().len(), 1);
        assert_eq!(molecule.conformers_2d()[0].id(), 0);
        assert!(molecule.conformers_3d().is_empty());
        assert_eq!(molecule.atoms()[0].id(), AtomId::new(0));
        assert_eq!(molecule.atoms()[1].id(), AtomId::new(1));
        assert_eq!(molecule.bonds()[0].id(), BondId::new(0));
        assert_eq!(molecule.bonds()[0].begin(), AtomId::new(0));
        assert_eq!(molecule.bonds()[0].end(), AtomId::new(1));
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
        assert!(molecule.substance_groups().is_empty());
        assert!(molecule.stereo_groups().is_empty());
        let before = format!("{molecule:#?}");
        let valence = cosmolkit_core::valence::assign_valence(&molecule, ValenceModel::RdkitLike)?;
        assert_eq!(valence.explicit_valence, vec![1, 1]);
        assert_eq!(valence.implicit_hydrogens, vec![3, 3]);
        let rings = cosmolkit_core::rings::symmetrize_sssr(&molecule)?;
        assert_eq!(rings, RingInfo::new(RingFindType::SymmSssr, 2, 1));
        let bit_rows: Vec<_> = molecule.coordinates_2d().unwrap().iter().map(|p| p.map(f64::to_bits)).collect();
        assert_eq!(bit_rows, vec![[0, 0], [4609434218613702656, 0]]);
        assert_eq!(before, format!("{molecule:#?}"), "canonical computations preserve source");
        // Actual legacy public call; count after invocation, before inspecting it.
        let result = molecule.to_svg(300, 300);
        calls += 1;
        let after = format!("{molecule:#?}");
        assert_eq!(before, after, "public source preservation: {label}");
        assert_eq!(bit_rows, molecule.coordinates_2d().unwrap().iter().map(|p| p.map(f64::to_bits)).collect::<Vec<_>>());
        let svg = result?;
        inputs.push(json!({
            "label": label, "width": 300, "height": 300,
            "highlights": {"atom0": atom_note, "bond0": bond_note, "molecule": null},
            "specs": specs,
            "coordinates": [[0.0,0.0],[1.5,0.0]], "coordinate_bits": bit_rows,
            "conformer_2d_ids": [0], "conformer_3d_ids": [],
            "actual_atoms": format!("{:#?}", molecule.atoms()),
            "actual_bonds": format!("{:#?}", molecule.bonds()),
            "actual_properties": format!("{:#?}", molecule.properties()),
            "actual_valence": format!("{valence:#?}"), "actual_rings": format!("{rings:#?}"),
            "source_before": before, "source_after": after, "source_preserved": true,
        }));
        println!("source call{calls}: {label} SVG bytes={} preserved=true", svg.len());
        artifacts.push((label, svg));
    }
    assert_eq!(calls, 8);
    let absent = &artifacts[0].1;
    for i in [1, 2, 3, 4, 7] {
        assert_eq!(&artifacts[i].1, absent, "source default/case control {}", artifacts[i].0);
    }
    assert_eq!(artifacts[5].1, artifacts[6].1, "one and true controls");
    assert_ne!(&artifacts[5].1, absent, "highlight differs");
    // All eight source outputs and controls precede publication.
    fs::create_dir(output)?;
    let mut outputs = Vec::new();
    for (label, svg) in artifacts {
        let name = format!("{label}.svg");
        let path = output.join(&name);
        write_new(&path, svg.as_bytes())?;
        outputs.push(json!({"label": label, "file": name, "bytes": svg.len(), "sha256": sha256(&path)?}));
    }
    let mut links = Vec::new();
    for path in [
        "tools/testdata/drawing_highlight_reference.rs",
        "target/legacy-drawing-reference-build/release/deps/libcosmolkit_core-e2e8e1910a9fd42c.rlib",
        "target/legacy-drawing-reference-build/release/deps/libserde_json-21be7c7129c38d8b.rlib",
        "/tmp/ck-draw-legacy-reference-O8Thas/cosmolkit-core-0.3.0/src/properties/draw.rs",
        "/tmp/ck-draw-legacy-reference-O8Thas/cosmolkit-core-0.3.0/Cargo.lock",
    ] {
        links.push(json!({"path": path, "sha256": sha256(Path::new(path))?}));
    }
    let manifest = json!({
        "schema": "drawing_highlight_legacy_reference_v1",
        "reference": {"package": "cosmolkit-core", "version": "0.3.0", "license": "MIT",
            "isolated_package": "/tmp/ck-draw-legacy-reference-O8Thas/cosmolkit-core-0.3.0",
            "vcs_metadata": "d892ec3507c5b568c5ed5d86ae44e466f7d03855"},
        "source_calls": calls, "source_default_and_case_controls": true, "source_one_equals_true_and_differs_absent": true,
        "migrated_outputs_examined": false, "inputs": inputs, "outputs": outputs, "links": links,
        "profile": {"width":300, "height":300, "coordinates":"supplied 2D", "options":"legacy default"},
        "platform": std::env::consts::ARCH,
        "executable_sha256": sha256(&std::env::current_exe()?)?,
    });
    write_new(&output.join("manifest.json"), &serde_json::to_vec_pretty(&manifest)?)?;
    println!("source8 complete; five default/case and one/true controls verified; no migrated library linked");
    Ok(())
}
