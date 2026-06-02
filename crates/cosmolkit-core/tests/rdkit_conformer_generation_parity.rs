use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::Command;

use cosmolkit_core::{
    EmbedParameters, ForceFieldVec3, Molecule,
    chemistry::distgeom::{embed_molecule, embed_multiple_confs},
    io::{molfile::read_mol_record_from_str_with_params, sdf::SdfReadParams},
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct ConformerGenerationFixtureRecord {
    kind: String,
    fixture: Option<String>,
    source: String,
    bytes: u64,
    sha256: String,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn inventory_path() -> PathBuf {
    repo_root().join("tests/fixtures/conformer_generation/rdkit_inventory.jsonl")
}

fn fixture_root() -> PathBuf {
    repo_root().join("tests/fixtures/conformer_generation")
}

fn load_inventory() -> Vec<ConformerGenerationFixtureRecord> {
    let path = inventory_path();
    let file =
        File::open(&path).unwrap_or_else(|err| panic!("failed to open {}: {err}", path.display()));
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(idx, line)| {
            let line = line.unwrap_or_else(|err| {
                panic!("failed to read {} line {}: {err}", path.display(), idx + 1)
            });
            serde_json::from_str(&line).unwrap_or_else(|err| {
                panic!("failed to parse {} line {}: {err}", path.display(), idx + 1)
            })
        })
        .collect()
}

fn sha256sum(path: &std::path::Path) -> String {
    let output = Command::new("sha256sum")
        .arg(path)
        .output()
        .unwrap_or_else(|err| panic!("failed to run sha256sum for {}: {err}", path.display()));
    assert!(
        output.status.success(),
        "sha256sum failed for {} with status {:?}: {}",
        path.display(),
        output.status.code(),
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap_or_else(|err| {
        panic!(
            "sha256sum output for {} is not UTF-8: {err}",
            path.display()
        )
    });
    stdout
        .split_whitespace()
        .next()
        .unwrap_or_else(|| panic!("sha256sum output for {} was empty", path.display()))
        .to_owned()
}

fn read_rdkit_mol_fixture(source: &str) -> cosmolkit_core::Molecule {
    let path = vendored_fixture_path(source);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read RDKit mol fixture {}: {err}", path.display()));
    read_mol_record_from_str_with_params(
        &text,
        SdfReadParams {
            sanitize: true,
            remove_hs: false,
            process_property_lists: false,
            ..Default::default()
        },
    )
    .unwrap_or_else(|err| {
        panic!(
            "failed to parse RDKit mol fixture {}: {err}",
            path.display()
        )
    })
    .molecule
}

fn vendored_fixture_path(source: &str) -> PathBuf {
    let record = load_inventory()
        .into_iter()
        .find(|record| record.source == source)
        .unwrap_or_else(|| panic!("fixture inventory missing source {source}"));
    let fixture = record.fixture.unwrap_or_else(|| {
        panic!("fixture inventory source {source} does not have a vendored fixture path")
    });
    fixture_root().join(fixture)
}

fn conformer_coords(mol: &Molecule) -> Vec<Vec<[f64; 3]>> {
    mol.conformers_3d()
        .iter()
        .map(|conformer| conformer.coordinates().to_vec())
        .collect()
}

fn assert_all_conformer_coords_finite(label: &str, mol: &Molecule) {
    assert!(
        mol.conformers_3d()
            .iter()
            .flat_map(|conf| conf.coordinates().iter())
            .flat_map(|coord| coord.iter())
            .all(|value| value.is_finite()),
        "{label} produced non-finite conformer coordinates"
    );
}

#[test]
fn conformer_generation_fixture_inventory_has_required_rdkit_sources() {
    let records = load_inventory();
    let sources = records
        .iter()
        .map(|record| record.source.as_str())
        .collect::<Vec<_>>();

    for required in [
        "third_party/rdkit/Code/DistGeom/testDistGeom.cpp",
        "third_party/rdkit/Code/GraphMol/DistGeomHelpers/testDgeomHelpers.cpp",
        "third_party/rdkit/Code/GraphMol/DistGeomHelpers/catch_tests.cpp",
        "third_party/rdkit/Code/GraphMol/DistGeomHelpers/Wrap/testDistGeom.py",
    ] {
        assert!(
            sources.contains(&required),
            "fixture inventory missing required RDKit source {required}"
        );
    }

    assert_eq!(
        records
            .iter()
            .filter(|record| record.kind == "test_data")
            .count(),
        22,
        "fixture inventory must include every DistGeomHelpers/test_data file"
    );
}

#[test]
fn conformer_generation_fixture_inventory_has_vendored_test_data_paths() {
    for record in load_inventory() {
        if record.kind != "test_data" {
            assert!(
                record.fixture.is_none(),
                "source-only record {} should not declare a vendored fixture path",
                record.source
            );
            continue;
        }
        let fixture = record.fixture.as_ref().unwrap_or_else(|| {
            panic!(
                "test_data record missing vendored fixture path: {}",
                record.source
            )
        });
        assert!(
            fixture.starts_with("rdkit/test_data/"),
            "unexpected vendored fixture layout for {}: {}",
            record.source,
            fixture
        );
    }
}

#[test]
fn conformer_generation_single_parity_embeds_rdkit_simple_torsion_fixtures() {
    let cases: [(&str, fn() -> EmbedParameters, &str); 8] = [
        (
            "DG",
            EmbedParameters::default,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.dg.mol",
        ),
        (
            "KDG",
            EmbedParameters::kdg,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.kdg.mol",
        ),
        (
            "ETDG",
            EmbedParameters::etdg,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etdg.mol",
        ),
        (
            "ETDGv2",
            EmbedParameters::etdg_v2,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/torsion.etkdg.v2.mol",
        ),
        (
            "ETKDG",
            EmbedParameters::etkdg,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol",
        ),
        (
            "ETKDGv2",
            EmbedParameters::etkdg_v2,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/torsion.etkdg.v2.mol",
        ),
        (
            "ETKDGv3",
            EmbedParameters::etkdg_v3,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle.etkdgv3.mol",
        ),
        (
            "srETKDGv3",
            EmbedParameters::sr_etkdg_v3,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.smallring.etkdgv3.mol",
        ),
    ];

    for (label, factory, source) in cases {
        let mol = read_rdkit_mol_fixture(source);
        let mut params = factory();
        params.random_seed = 42;
        params.num_threads = 1;

        let atom_count = mol.num_atoms();
        let original_conformer_count = mol.conformers_3d().len();
        let (embedded, conf_id) = embed_molecule(&mol, &mut params)
            .unwrap_or_else(|err| panic!("{label} failed to embed {source}: {err}"));
        assert_eq!(
            conf_id, 0,
            "{label} should match RDKit source tests that expect EmbedMolecule(...) == 0"
        );
        assert_eq!(
            embedded.conformers_3d().len(),
            1,
            "{label} should add exactly one conformer for {source}"
        );
        let conf = &embedded.conformers_3d()[0];
        assert!(conf.is_3d(), "{label} conformer should be 3D for {source}");
        assert_eq!(
            conf.coordinates().len(),
            atom_count,
            "{label} conformer coordinate count should match atom count for {source}"
        );
        assert!(
            conf.coordinates()
                .iter()
                .flat_map(|coord| coord.iter())
                .all(|value| value.is_finite()),
            "{label} conformer contains non-finite coordinates for {source}"
        );
        assert_eq!(
            mol.conformers_3d().len(),
            original_conformer_count,
            "{label} embedding should preserve value semantics for {source}"
        );
    }
}

#[test]
fn conformer_generation_multi_parity_covers_rdkit_multi_conformer_controls() {
    let source =
        "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol";
    let mol = read_rdkit_mol_fixture(source);

    let mut params = EmbedParameters::etkdg();
    params.random_seed = 0xf00d;
    params.num_threads = 1;
    params.track_failures = true;
    params.timeout = 1;
    let original_conformer_count = mol.conformers_3d().len();
    let (embedded, ids) = embed_multiple_confs(&mol, 10, &mut params)
        .unwrap_or_else(|err| panic!("ETKDG multi-conformer embed failed for {source}: {err}"));
    assert_eq!(ids, (0..10).collect::<Vec<i32>>());
    assert_eq!(embedded.conformers_3d().len(), 10);
    assert_eq!(mol.conformers_3d().len(), original_conformer_count);
    assert_eq!(params.failures.iter().sum::<u32>(), 0);

    let mut repeated_params = EmbedParameters::etkdg();
    repeated_params.random_seed = 0xf00d;
    repeated_params.num_threads = 1;
    let (repeated, repeated_ids) = embed_multiple_confs(&mol, 10, &mut repeated_params)
        .unwrap_or_else(|err| panic!("repeat ETKDG embed failed for {source}: {err}"));
    assert_eq!(repeated_ids, ids);
    assert_eq!(conformer_coords(&repeated), conformer_coords(&embedded));

    let mut pruned_params = EmbedParameters::etkdg();
    pruned_params.random_seed = 0xf00d;
    pruned_params.num_threads = 1;
    pruned_params.prune_rms_thresh = 0.5;
    pruned_params.use_symmetry_for_pruning = true;
    let (pruned, pruned_ids) = embed_multiple_confs(&mol, 10, &mut pruned_params)
        .unwrap_or_else(|err| panic!("pruned ETKDG embed failed for {source}: {err}"));
    assert!(!pruned_ids.is_empty());
    assert!(pruned_ids.len() <= 10);
    assert_eq!(pruned.conformers_3d().len(), pruned_ids.len());

    let mut sequential_params = EmbedParameters::etkdg();
    sequential_params.random_seed = 0xf00d;
    sequential_params.num_threads = 1;
    sequential_params.enable_sequential_random_seeds = true;
    let (sequential, sequential_ids) = embed_multiple_confs(&mol, 4, &mut sequential_params)
        .unwrap_or_else(|err| panic!("sequential-seed ETKDG embed failed for {source}: {err}"));
    assert_eq!(sequential_ids, vec![0, 1, 2, 3]);
    let mut sequential_repeat_params = sequential_params.clone();
    let (sequential_repeat, sequential_repeat_ids) =
        embed_multiple_confs(&mol, 4, &mut sequential_repeat_params).unwrap_or_else(|err| {
            panic!("repeat sequential-seed embed failed for {source}: {err}")
        });
    assert_eq!(sequential_repeat_ids, sequential_ids);
    assert_eq!(
        conformer_coords(&sequential_repeat),
        conformer_coords(&sequential)
    );

    let fragment_mol = Molecule::from_smiles("CC.CC").expect("fragment molecule");
    let mut fragment_params = EmbedParameters::etkdg();
    fragment_params.random_seed = 0xf00d;
    fragment_params.num_threads = 1;
    fragment_params.embed_fragments_separately = true;
    let (fragment_embedded, fragment_ids) =
        embed_multiple_confs(&fragment_mol, 3, &mut fragment_params)
            .expect("fragment conformer embedding");
    assert_eq!(fragment_ids, vec![0, 1, 2]);
    assert_eq!(fragment_embedded.conformers_3d().len(), 3);
}

#[test]
fn conformer_generation_parameter_parity_covers_public_parameter_controls() {
    let source =
        "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol";
    let mol = read_rdkit_mol_fixture(source);

    let mut coord_params = EmbedParameters::etkdg_v3();
    coord_params.random_seed = 0xC0FFEE;
    coord_params.num_threads = 1;
    coord_params.use_random_coords = true;
    coord_params.coord_map = Some(BTreeMap::from([
        (0, ForceFieldVec3::new(0.0, 0.0, 0.0)),
        (1, ForceFieldVec3::new(0.0, 0.0, 1.5)),
        (2, ForceFieldVec3::new(0.0, 1.5, 1.5)),
    ]));
    let (coord_embedded, coord_id) = embed_molecule(&mol, &mut coord_params)
        .unwrap_or_else(|err| panic!("coordMap/random-coordinate embed failed: {err}"));
    assert_eq!(coord_id, 0);
    let coord_conf = &coord_embedded.conformers_3d()[0];
    assert_eq!(coord_conf.coordinates()[0], [0.0, 0.0, 0.0]);
    assert_eq!(coord_conf.coordinates()[1], [0.0, 0.0, 1.5]);
    assert_eq!(coord_conf.coordinates()[2], [0.0, 1.5, 1.5]);

    let mut cpci_params = EmbedParameters::etkdg_v3();
    cpci_params.random_seed = 0xC0FFEE;
    cpci_params.num_threads = 1;
    cpci_params.cpci = Some(BTreeMap::from([((0, 3), 0.5), ((1, 4), -0.25)]));
    let (cpci_embedded, cpci_id) = embed_molecule(&mol, &mut cpci_params)
        .unwrap_or_else(|err| panic!("CPCI embed failed: {err}"));
    assert_eq!(cpci_id, 0);
    assert_eq!(cpci_embedded.conformers_3d().len(), 1);

    let mut first_params = EmbedParameters::etkdg();
    first_params.random_seed = 42;
    first_params.num_threads = 1;
    let (with_existing, first_ids) = embed_multiple_confs(&mol, 2, &mut first_params)
        .unwrap_or_else(|err| panic!("initial conformer embed failed: {err}"));
    assert_eq!(first_ids, vec![0, 1]);

    let mut append_params = EmbedParameters::etkdg();
    append_params.random_seed = 43;
    append_params.num_threads = 1;
    append_params.clear_confs = false;
    let (appended, appended_ids) = embed_multiple_confs(&with_existing, 2, &mut append_params)
        .unwrap_or_else(|err| panic!("clearConfs=false embed failed: {err}"));
    assert_eq!(appended_ids, vec![2, 3]);
    assert_eq!(appended.conformers_3d().len(), 4);

    let mut clear_params = EmbedParameters::etkdg();
    clear_params.random_seed = 44;
    clear_params.num_threads = 1;
    clear_params.clear_confs = true;
    let (cleared, cleared_ids) = embed_multiple_confs(&with_existing, 1, &mut clear_params)
        .unwrap_or_else(|err| panic!("clearConfs=true embed failed: {err}"));
    assert_eq!(cleared_ids, vec![0]);
    assert_eq!(cleared.conformers_3d().len(), 1);
}

#[test]
fn conformer_generation_stereo_parity_covers_source_backed_stereo_and_torsion_paths() {
    let chirality = read_rdkit_mol_fixture(
        "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/chirality_failure_test.mol",
    );
    let mut chirality_params = EmbedParameters::etkdg_v3();
    chirality_params.random_seed = 0xf00d;
    chirality_params.num_threads = 1;
    chirality_params.track_failures = true;
    chirality_params.max_iterations = 50;
    let (chirality_embedded, chirality_id) = embed_molecule(&chirality, &mut chirality_params)
        .unwrap_or_else(|err| panic!("chirality-failure fixture returned error: {err}"));
    assert_eq!(
        chirality_id, -1,
        "RDKit catch_tests.cpp expects chirality_failure_test.mol embedding to fail"
    );
    assert!(chirality_embedded.conformers_3d().is_empty());
    assert!(
        chirality_params.failures.iter().sum::<u32>() > 0,
        "RDKit tracks conformer-generation failures for the chirality fixture"
    );

    let double_bond =
        Molecule::from_smiles("O=C1OCC/C=C/CC/C=C/C(=N/OCC(=O)N2CCCCC2)Cc2cc(O)cc(O)c21")
            .expect("RDKit platinum-set double-bond stereo SMILES");
    let mut double_bond_params = EmbedParameters::etkdg_v3();
    double_bond_params.random_seed = 0xf00d + 81;
    double_bond_params.num_threads = 1;
    let (double_bond_embedded, double_bond_id) =
        embed_molecule(&double_bond, &mut double_bond_params)
            .unwrap_or_else(|err| panic!("double-bond stereo embed failed: {err}"));
    assert!(
        double_bond_id >= 0,
        "RDKit catch_tests.cpp expects the platinum-set double-bond stereo molecule to embed"
    );
    assert_all_conformer_coords_finite("double-bond stereo", &double_bond_embedded);

    for (label, factory, source) in [
        (
            "macrocycle ETKDG",
            EmbedParameters::etkdg as fn() -> EmbedParameters,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle.etkdg.mol",
        ),
        (
            "macrocycle ETKDGv3",
            EmbedParameters::etkdg_v3 as fn() -> EmbedParameters,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle.etkdgv3.mol",
        ),
        (
            "small-ring srETKDGv3",
            EmbedParameters::sr_etkdg_v3 as fn() -> EmbedParameters,
            "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.smallring.etkdgv3.mol",
        ),
    ] {
        let mol = read_rdkit_mol_fixture(source);
        let mut params = factory();
        params.random_seed = 42;
        params.num_threads = 1;
        let (embedded, id) = embed_molecule(&mol, &mut params)
            .unwrap_or_else(|err| panic!("{label} fixture embed failed for {source}: {err}"));
        assert_eq!(
            id, 0,
            "RDKit testDgeomHelpers.cpp expects {label} fixture embedding to return 0"
        );
        assert_all_conformer_coords_finite(label, &embedded);
    }

    for (label, smiles) in [("allene", "C=C=C"), ("sulfone", "CS(=O)(=O)C")] {
        let mol = Molecule::from_smiles(smiles).unwrap_or_else(|err| panic!("{label}: {err}"));
        let mut params = EmbedParameters::etkdg();
        params.random_seed = 0xf00d;
        params.num_threads = 1;
        let (embedded, id) = embed_molecule(&mol, &mut params)
            .unwrap_or_else(|err| panic!("{label} linear double-bond path failed: {err}"));
        assert!(
            id >= 0,
            "RDKit FindDoubleBonds source fixture {label} should remain embeddable"
        );
        assert_all_conformer_coords_finite(label, &embedded);
    }

    let amide = Molecule::from_smiles("CC(=O)NC").expect("force-trans-amide fixture");
    let mut trans_params = EmbedParameters::default();
    trans_params.force_trans_amides = true;
    trans_params.random_seed = 0xf00d;
    trans_params.num_threads = 1;
    trans_params.use_exp_torsion_angle_prefs = false;
    trans_params.use_basic_knowledge = true;
    let (trans_embedded, trans_ids) = embed_multiple_confs(&amide, 10, &mut trans_params)
        .unwrap_or_else(|err| panic!("forceTransAmides=true embed failed: {err}"));
    assert_eq!(trans_ids, (0..10).collect::<Vec<i32>>());
    assert_all_conformer_coords_finite("forceTransAmides=true", &trans_embedded);

    let mut relaxed_params = trans_params.clone();
    relaxed_params.force_trans_amides = false;
    let (relaxed_embedded, relaxed_ids) = embed_multiple_confs(&amide, 10, &mut relaxed_params)
        .unwrap_or_else(|err| panic!("forceTransAmides=false embed failed: {err}"));
    assert_eq!(relaxed_ids, (0..10).collect::<Vec<i32>>());
    assert_all_conformer_coords_finite("forceTransAmides=false", &relaxed_embedded);
}

#[test]
fn conformer_generation_fixture_inventory_matches_vendored_fixtures() {
    let records = load_inventory();
    assert_eq!(
        records.len(),
        26,
        "unexpected conformer fixture inventory size"
    );
    for record in records {
        let Some(fixture) = record.fixture.as_ref() else {
            continue;
        };
        let path = fixture_root().join(fixture);
        let metadata = std::fs::metadata(&path).unwrap_or_else(|err| {
            panic!(
                "missing vendored RDKit conformer fixture {}: {err}",
                path.display()
            )
        });
        assert_eq!(
            metadata.len(),
            record.bytes,
            "byte length changed for {}",
            fixture
        );
        assert_eq!(
            sha256sum(&path),
            record.sha256,
            "SHA-256 changed for {}",
            fixture
        );
    }
}
