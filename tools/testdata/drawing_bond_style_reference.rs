//! One-shot source-first fixed bond styles; links immutable published0.3.0 only.
//! Destination must be absent. Every reference artifact uses create_new.
use cosmolkit_core::{AtomId, AtomSpec, BondId, BondOrder, BondSpec, Element, MoleculeBuilder, RingFindType, RingInfo, ValenceModel};
use serde_json::json;
use std::{error::Error, fs::{self, OpenOptions}, io::Write, path::Path, process::Command};

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), Box<dyn Error>> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    file.write_all(bytes)?;
    Ok(())
}
fn sha256(path: &Path) -> Result<String, Box<dyn Error>> {
    let result = Command::new("sha256sum").arg(path).output()?;
    if !result.status.success() { return Err(format!("sha256sum failed: {result:?}").into()); }
    Ok(String::from_utf8(result.stdout)?.split_whitespace().next().ok_or("missing SHA256")?.to_owned())
}
fn coordinate_identity(molecule: &cosmolkit_core::Molecule) -> Vec<(usize, Vec<[u64; 2]>)> {
    molecule.conformers_2d().iter().map(|c| (c.id(), c.coordinates().iter().map(|p| p.map(f64::to_bits)).collect())).collect()
}
fn main() -> Result<(), Box<dyn Error>> {
    let destination = std::env::args().nth(1).ok_or("usage: driver NEW_OUTPUT_DIRECTORY")?;
    let output = Path::new(&destination);
    if output.exists() { return Err("destination exists; reference refresh forbidden".into()); }
    let cases = [
        ("triple_fwd", BondOrder::Triple, 0, 1, [3,3], [1,1]),
        ("triple_rev", BondOrder::Triple, 1, 0, [3,3], [1,1]),
        ("quadruple_fwd", BondOrder::Quadruple, 0, 1, [4,4], [0,0]),
        ("quadruple_rev", BondOrder::Quadruple, 1, 0, [4,4], [0,0]),
        ("hydrogen_fwd", BondOrder::Hydrogen, 0, 1, [0,0], [4,4]),
        ("hydrogen_rev", BondOrder::Hydrogen, 1, 0, [0,0], [4,4]),
        ("unspecified_fwd", BondOrder::Unspecified, 0, 1, [0,0], [4,4]),
        ("unspecified_rev", BondOrder::Unspecified, 1, 0, [0,0], [4,4]),
        ("dative_fwd", BondOrder::Dative, 0, 1, [0,1], [4,3]),
        ("dative_rev", BondOrder::Dative, 1, 0, [1,0], [3,4]),
        ("dative_one_fwd", BondOrder::DativeOne, 0, 1, [0,1], [4,3]),
        ("dative_one_rev", BondOrder::DativeOne, 1, 0, [1,0], [3,4]),
    ];
    let mut prepared = Vec::new();
    // Verify ALL canonical old rows against frozen literals before ANY SVG call.
    for (label, order, begin, end, explicit, implicit) in cases {
        let atom0 = AtomSpec::new(Element::C);
        let atom1 = AtomSpec::new(Element::C);
        let bond_spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), order);
        let specs = json!({"atoms":[format!("{atom0:#?}"),format!("{atom1:#?}")], "bond":format!("{bond_spec:#?}")});
        let mut builder = MoleculeBuilder::new();
        assert_eq!(builder.add_atom(atom0), AtomId::new(0));
        assert_eq!(builder.add_atom(atom1), AtomId::new(1));
        assert_eq!(builder.add_bond(bond_spec)?, BondId::new(0));
        for id in [AtomId::new(0),AtomId::new(1)] {
            assert_eq!(builder.degree(id), 1);
            assert_eq!(builder.neighbor_bonds(id), &[BondId::new(0)]);
        }
        builder.set_2d_coordinates(vec![[0.0,0.0],[1.5,0.0]])?;
        let molecule = builder.build()?;
        assert_eq!(molecule.num_atoms(),2);
        assert_eq!(molecule.num_bonds(),1);
        assert_eq!(molecule.atoms()[0].id(),AtomId::new(0));
        assert_eq!(molecule.atoms()[1].id(),AtomId::new(1));
        assert_eq!(molecule.bonds()[0].id(),BondId::new(0));
        assert_eq!(molecule.bonds()[0].begin(),AtomId::new(begin));
        assert_eq!(molecule.bonds()[0].end(),AtomId::new(end));
        assert_eq!(molecule.bonds()[0].order(),order);
        assert!(molecule.atoms().iter().all(|a| a.props().is_empty() && a.computed_prop_names().is_empty() && a.query().is_none()));
        assert!(molecule.bonds().iter().all(|b| b.props().is_empty() && b.computed_prop_names().is_empty() && b.query().is_none()));
        assert!(molecule.properties().props().is_empty());
        assert!(molecule.properties().computed_prop_names().is_empty());
        assert!(molecule.substance_groups().is_empty());
        assert!(molecule.stereo_groups().is_empty());
        assert_eq!(molecule.source_coordinate_dim(),Some(cosmolkit_core::CoordinateDimension::TwoD));
        assert_eq!(coordinate_identity(&molecule),vec![(0,vec![[0,0],[4609434218613702656,0]])]);
        assert!(molecule.conformers_3d().is_empty());
        assert!(molecule.conformers_2d()[0].props().is_empty());
        let original = format!("{molecule:#?}");
        let valence = cosmolkit_core::valence::assign_valence(&molecule,ValenceModel::RdkitLike)?;
        assert_eq!(valence.explicit_valence,explicit,"frozen explicit {label}");
        assert_eq!(valence.implicit_hydrogens,implicit,"frozen implicit {label}");
        let rings = cosmolkit_core::rings::symmetrize_sssr(&molecule)?;
        assert_eq!(rings,RingInfo::new(RingFindType::SymmSssr,2,1));
        assert!(rings.is_initialized());
        assert_eq!(rings.find_type(),RingFindType::SymmSssr);
        assert!(rings.atom_rings().is_empty());
        assert!(rings.bond_rings().is_empty());
        for id in [AtomId::new(0),AtomId::new(1)] {
            assert_eq!(rings.num_atom_rings(id),0);
            assert!(rings.atom_members(id).is_empty());
        }
        assert_eq!(rings.num_bond_rings(BondId::new(0)),0);
        assert!(rings.bond_members(BondId::new(0)).is_empty());
        assert_eq!(original,format!("{molecule:#?}"),"canonical computations preserve {label}");
        println!("source prerequisite {label} begin={begin} end={end} order={order:?} explicit={explicit:?} implicit={implicit:?} rings=emptySymmSSSR");
        prepared.push((label,order,begin,end,molecule,specs,valence,rings));
    }
    assert_eq!(prepared.len(),12);
    let mut inputs = Vec::new();
    let mut artifacts = Vec::new();
    let mut calls = 0;
    for (label,order,begin,end,molecule,specs,valence,rings) in prepared {
        let before = format!("{molecule:#?}");
        let identity = coordinate_identity(&molecule);
        let result = molecule.to_svg(300,300);
        calls += 1;
        let after = format!("{molecule:#?}");
        assert_eq!(before,after,"whole original source preservation {label}");
        assert_eq!(identity,coordinate_identity(&molecule),"coordinate IDs/bits {label}");
        let svg = result?;
        inputs.push(json!({
            "label":label,"semantic_order":format!("{order:?}"),"begin":begin,"end":end,"bond_id":0,"atom_ids":[0,1],
            "width":300,"height":300,"specs":specs,"coordinates":[[0.0,0.0],[1.5,0.0]],
            "coordinate_bits":[[0u64,0u64],[4609434218613702656u64,0u64]],"conformer_2d_ids":[0],"conformer_3d_ids":[],
            "actual_atoms":format!("{:#?}",molecule.atoms()),"actual_bonds":format!("{:#?}",molecule.bonds()),
            "actual_properties":format!("{:#?}",molecule.properties()),
            "explicit_valence":valence.explicit_valence,"implicit_hydrogens":valence.implicit_hydrogens,
            "actual_valence":format!("{valence:#?}"),"actual_rings":format!("{rings:#?}"),
            "source_before":before,"source_after":after,"source_preserved":true,"coordinate_identity_before":identity,
            "coordinate_identity_after":coordinate_identity(&molecule)
        }));
        println!("source call{calls}: {label} SVG bytes={} preserved=true",svg.len());
        artifacts.push((label,svg));
    }
    assert_eq!(calls,12);
    for (a,b) in [(4,6),(5,7),(8,10),(9,11)] {
        assert_eq!(artifacts[a].1,artifacts[b].1,"same-direction source alias {} {}",artifacts[a].0,artifacts[b].0);
    }
    fs::create_dir(output)?;
    let mut outputs = Vec::new();
    for (label,svg) in artifacts {
        let name = format!("{label}.svg");
        let path = output.join(&name);
        write_new(&path,svg.as_bytes())?;
        outputs.push(json!({"label":label,"file":name,"bytes":svg.len(),"sha256":sha256(&path)?}));
    }
    let mut links = Vec::new();
    for path in [
        "tools/testdata/drawing_bond_style_reference.rs",
        "target/legacy-drawing-reference-build/release/deps/libcosmolkit_core-e2e8e1910a9fd42c.rlib",
        "target/legacy-drawing-reference-build/release/deps/libserde_json-21be7c7129c38d8b.rlib",
        "/tmp/ck-draw-legacy-reference-O8Thas/cosmolkit-core-0.3.0/src/properties/draw.rs",
        "/tmp/ck-draw-legacy-reference-O8Thas/cosmolkit-core-0.3.0/src/chemistry/valence.rs",
        "/tmp/ck-draw-legacy-reference-O8Thas/cosmolkit-core-0.3.0/Cargo.lock",
    ] { links.push(json!({"path":path,"sha256":sha256(Path::new(path))?})); }
    let manifest = json!({
        "schema":"drawing_bond_style_legacy_reference_v1","source_calls":calls,"input_count":12,"output_count":12,
        "reference":{"package":"cosmolkit-core","version":"0.3.0","license":"MIT",
        "isolated_package":"/tmp/ck-draw-legacy-reference-O8Thas/cosmolkit-core-0.3.0","vcs_metadata":"d892ec3507c5b568c5ed5d86ae44e466f7d03855"},
        "same_direction_alias_controls":[["hydrogen_fwd","unspecified_fwd"],["hydrogen_rev","unspecified_rev"],["dative_fwd","dative_one_fwd"],["dative_rev","dative_one_rev"]],
        "source_alias_controls_passed":true,"all_canonical_prerequisites_before_any_svg":true,"migrated_outputs_examined":false,
        "inputs":inputs,"outputs":outputs,"links":links,
        "profile":{"width":300,"height":300,"coordinates":"supplied2D","options":"legacy default"},
        "platform":{"arch":std::env::consts::ARCH,"os":std::env::consts::OS},
        "executable_sha256":sha256(&std::env::current_exe()?)?
    });
    write_new(&output.join("manifest.json"),&serde_json::to_vec_pretty(&manifest)?)?;
    println!("source12 complete; all four same-direction aliases exact; no migrated library linked");
    Ok(())
}

