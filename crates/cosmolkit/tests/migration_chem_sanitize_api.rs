use std::error::Error as _;
use std::path::PathBuf;
use std::process::{Command, Output};

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec,
    ChemistryProblemError, Conformer2D, CoordinateBlock, Element, Molecule, MoleculeOpKind,
    MoleculeOpOutput, MoleculeProperties, OperationDomain, OperationError, ParityPolicy,
    SanitizeError, SanitizeOperations, SanitizeParams, SanitizeStage, StateModel, SupportStatus,
    TopologyBlock, TopologyEditKind, feature_spec, operation_invariant, operation_parity,
    operation_spec, support_matrix,
};

fn topology_from_specs(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        atom_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect(),
        bond_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
            .collect(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap()
}

fn atom_spec(atomic_number: u8) -> AtomSpec {
    AtomSpec::new(Element::from_atomic_number(atomic_number).unwrap())
}

fn bond_spec(begin: usize, end: usize, order: BondOrder) -> BondSpec {
    BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
}

fn alternating_benzene() -> TopologyBlock {
    topology_from_specs(
        (0..6)
            .map(|index| {
                AtomSpec::new(Element::C)
                    .with_prop("atom-note", format!("atom-{index}"))
                    .unwrap()
                    .with_computed_prop("_CIPCode", "R")
                    .unwrap()
            })
            .collect(),
        (0..6)
            .map(|index| {
                BondSpec::new(
                    AtomId::new(index),
                    AtomId::new((index + 1) % 6),
                    if index % 2 == 0 {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    },
                )
                .with_prop("bond-note", format!("bond-{index}"))
                .unwrap()
                .with_computed_prop("_CIPCode", "E")
                .unwrap()
            })
            .collect(),
    )
}

fn molecule(topology: TopologyBlock) -> Molecule {
    let atom_count = topology.atoms.len();
    Molecule::from_parts(
        topology,
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                37,
                (0..atom_count)
                    .map(|index| [index as f64, -(index as f64)])
                    .collect(),
            )],
            ..CoordinateBlock::default()
        },
        MoleculeProperties::default()
            .with_name("sanitize-public")
            .with_prop("source", "preserved")
            .unwrap()
            .with_computed_prop("_CIPComputed", "true")
            .unwrap(),
    )
    .unwrap()
}

fn problem_topology() -> TopologyBlock {
    topology_from_specs(
        vec![
            atom_spec(6),
            atom_spec(8),
            atom_spec(6),
            atom_spec(6),
            atom_spec(9),
            atom_spec(6),
            AtomSpec::new(Element::C).with_aromatic(true),
            AtomSpec::new(Element::C).with_aromatic(true),
            AtomSpec::new(Element::C).with_aromatic(true),
        ],
        vec![
            bond_spec(0, 1, BondOrder::Single),
            bond_spec(1, 2, BondOrder::Single),
            bond_spec(1, 3, BondOrder::Single),
            bond_spec(3, 4, BondOrder::Single),
            bond_spec(4, 5, BondOrder::Single),
            bond_spec(5, 6, BondOrder::Single),
            BondSpec::new(AtomId::new(6), AtomId::new(7), BondOrder::Aromatic).with_aromatic(true),
            BondSpec::new(AtomId::new(7), AtomId::new(8), BondOrder::Aromatic).with_aromatic(true),
            BondSpec::new(AtomId::new(8), AtomId::new(6), BondOrder::Aromatic).with_aromatic(true),
        ],
    )
}

#[test]
fn canonical_signatures_defaults_and_flag_vocabulary_are_exact() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> = Molecule::sanitize;
    let _: for<'a, 'b> fn(&'a Molecule, &'b SanitizeParams) -> Result<Molecule, OperationError> =
        Molecule::sanitize_with_params;
    let _: fn(&Molecule) -> Result<cosmolkit::ChemistryProblemReport, SanitizeError> =
        Molecule::detect_chemistry_problems;
    let _: for<'a, 'b> fn(
        &'a Molecule,
        &'b SanitizeParams,
    ) -> Result<cosmolkit::ChemistryProblemReport, SanitizeError> =
        Molecule::detect_chemistry_problems_with_params;

    assert_eq!(
        SanitizeParams::default().operations,
        SanitizeOperations::ALL
    );
    assert_eq!(SanitizeOperations::NONE.bits(), 0x000);
    assert_eq!(SanitizeOperations::CLEANUP.bits(), 0x001);
    assert_eq!(SanitizeOperations::PROPERTIES.bits(), 0x002);
    assert_eq!(SanitizeOperations::SYMM_RINGS.bits(), 0x004);
    assert_eq!(SanitizeOperations::KEKULIZE.bits(), 0x008);
    assert_eq!(SanitizeOperations::FIND_RADICALS.bits(), 0x010);
    assert_eq!(SanitizeOperations::SET_AROMATICITY.bits(), 0x020);
    assert_eq!(SanitizeOperations::SET_CONJUGATION.bits(), 0x040);
    assert_eq!(SanitizeOperations::SET_HYBRIDIZATION.bits(), 0x080);
    assert_eq!(SanitizeOperations::CLEANUP_CHIRALITY.bits(), 0x100);
    assert_eq!(SanitizeOperations::ADJUST_HS.bits(), 0x200);
    assert_eq!(SanitizeOperations::CLEANUP_ORGANOMETALLICS.bits(), 0x400);
    assert_eq!(SanitizeOperations::CLEANUP_ATROPISOMERS.bits(), 0x800);
    assert_eq!(SanitizeOperations::ALL.bits(), 0x0fff_ffff);
    assert!(matches!(
        SanitizeOperations::from_bits(0x1000),
        Err(SanitizeError::InvalidOperations {
            bits: 0x1000,
            unknown_bits: 0x1000,
        })
    ));
}

#[test]
fn binding_contract_exposes_exactly_the_frozen_eleven_entries() {
    let expected = [
        "types.SanitizeOperations",
        "types.SanitizeStage",
        "types.SanitizeParams",
        "types.SanitizeError",
        "types.ChemistryProblemError",
        "types.ChemistryProblem",
        "types.ChemistryProblemReport",
        "Molecule.sanitize",
        "Molecule.sanitize_with_params",
        "Molecule.detect_chemistry_problems",
        "Molecule.detect_chemistry_problems_with_params",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "sanitize")
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    for row in &rows {
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    for row in &rows[..7] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
    }
    for row in &rows[7..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
    }
    assert_eq!(
        rows[7].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(
        rows[8].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(rows[9].callable.unwrap().state_model, StateModel::ReadOnly);
    assert_eq!(rows[10].callable.unwrap().state_model, StateModel::ReadOnly);
    assert!(rows[9].callable.unwrap().operation_semantic_id.is_none());
    assert!(rows[10].callable.unwrap().operation_semantic_id.is_none());
}

#[test]
fn generated_registry_and_all_four_matrices_share_one_exact_operation() {
    let feature = feature_spec("sanitize").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("sanitize_with_params").unwrap();
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::Local);
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.access.read(), BlockSet::NONE);
    assert_eq!(
        spec.access.write(),
        BlockSet::TOPOLOGY
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE)
    );
    assert_eq!(spec.may_mutate, spec.access.write());
    assert_eq!(spec.auto_remap, BlockSet::NONE);
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.derived_effects.recompute.bits(), 0);
    assert_eq!(spec.derived_effects.preserve.bits(), 1 << 5);
    assert_eq!(spec.derived_effects.invalidate.bits(), 0xdf);
    assert_eq!(format!("{:?}", spec.cip_state), "ClearComputed");
    assert!(spec.io_roundtrip);
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "weak_sanitize_topology_state"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "sanitize_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| row.feature.name == "sanitize")
        .unwrap();
    assert!(std::ptr::eq(support.operation.unwrap(), spec));
}

#[test]
fn value_sanitize_preserves_identity_coordinates_and_ordinary_properties() {
    let source = molecule(alternating_benzene());
    let observer = source.clone();
    let output = source.sanitize().unwrap();
    let explicit = source
        .sanitize_with_params(&SanitizeParams::default())
        .unwrap();

    assert_eq!(output, explicit);
    assert_eq!(source, observer);
    assert_eq!(source.num_atoms(), output.num_atoms());
    assert_eq!(source.num_bonds(), output.num_bonds());
    assert!(
        source
            .atoms()
            .iter()
            .zip(output.atoms())
            .all(|(before, after)| before.id() == after.id())
    );
    assert!(
        source
            .bonds()
            .iter()
            .zip(output.bonds())
            .all(|(before, after)| {
                before.id() == after.id()
                    && before.begin() == after.begin()
                    && before.end() == after.end()
            })
    );
    assert_eq!(source.topology().adjacency, output.topology().adjacency);
    assert!(std::ptr::eq(source.coordinates(), output.coordinates()));
    assert_eq!(source.conformers(), output.conformers());
    assert_eq!(output.property("source"), Some("preserved"));
    assert_eq!(output.atoms()[0].prop("atom-note"), Some("atom-0"));
    assert_eq!(output.bonds()[0].prop("bond-note"), Some("bond-0"));
    assert_eq!(output.property("_CIPComputed"), None);
    assert_eq!(output.atoms()[0].prop("_CIPCode"), None);
    assert_eq!(output.bonds()[0].prop("_CIPCode"), None);
    assert!(output.atoms().iter().all(Atom::is_aromatic));
    assert!(output.bonds().iter().all(Bond::is_aromatic));
    assert!(format!("{output:?}").contains("derived_cache_is_empty: true"));
    assert!(!std::ptr::eq(source.topology(), output.topology()));
    assert!(!std::ptr::eq(source.properties(), output.properties()));
    assert_eq!(source.property("_CIPComputed"), Some("true"));
}

#[test]
fn none_selection_still_clears_only_computed_topology_and_cip_properties() {
    let source = molecule(alternating_benzene());
    let output = source
        .sanitize_with_params(&SanitizeParams {
            operations: SanitizeOperations::NONE,
        })
        .unwrap();
    assert_eq!(output.property("source"), Some("preserved"));
    assert_eq!(output.property("_CIPComputed"), None);
    assert_eq!(output.atoms()[0].prop("atom-note"), Some("atom-0"));
    assert_eq!(output.bonds()[0].prop("bond-note"), Some("bond-0"));
    assert_eq!(output.atoms()[0].prop("_CIPCode"), None);
    assert_eq!(output.bonds()[0].prop("_CIPCode"), None);
    assert_eq!(output.coordinates(), source.coordinates());
    assert!(std::ptr::eq(output.coordinates(), source.coordinates()));
    assert!(output.atoms().iter().all(|atom| !atom.is_aromatic()));
    assert!(output.bonds().iter().all(|bond| !bond.is_aromatic()));
}

#[test]
fn read_only_problem_detection_preserves_source_and_exact_problem_order() {
    let source = molecule(problem_topology());
    let observer = source.clone();
    let report = source.detect_chemistry_problems().unwrap();
    assert_eq!(report.problems.len(), 3);
    assert_eq!(report.problems[0].operation, SanitizeStage::Properties);
    assert!(matches!(
        report.problems[0].error,
        ChemistryProblemError::Valence(_)
    ));
    assert_eq!(report.problems[1].operation, SanitizeStage::Properties);
    assert!(matches!(
        report.problems[1].error,
        ChemistryProblemError::Valence(_)
    ));
    assert_eq!(report.problems[2].operation, SanitizeStage::Kekulize);
    assert!(matches!(
        report.problems[2].error,
        ChemistryProblemError::Kekulize(_)
    ));
    let explicit = source
        .detect_chemistry_problems_with_params(&SanitizeParams::default())
        .unwrap();
    assert_eq!(report, explicit);
    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    assert!(std::ptr::eq(source.coordinates(), observer.coordinates()));
    assert!(std::ptr::eq(source.properties(), observer.properties()));
}

#[test]
fn sanitize_failure_is_structured_and_atomic() {
    let source = molecule(problem_topology());
    let observer = source.clone();
    let error = source.sanitize().unwrap_err();
    assert!(matches!(
        &error,
        OperationError::Sanitize(SanitizeError::Properties {
            stage: SanitizeStage::Properties,
            ..
        })
    ));
    assert!(error.source().is_some());
    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    assert!(std::ptr::eq(source.coordinates(), observer.coordinates()));
    assert!(std::ptr::eq(source.properties(), observer.properties()));
}

fn privacy_cargo_check(case: &str, strict: bool) -> Output {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let workspace = manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .expect("cosmolkit crate must be nested under the workspace crates directory");
    let target = workspace.join("target/sanitize-privacy-compile");
    let inherited = std::env::var("RUSTFLAGS").unwrap_or_default();
    let rustflags = format!(
        "{inherited} --cfg cosmolkit_runtime_privacy_probe \
         --cfg=cosmolkit_runtime_privacy_case=\"{case}\""
    );
    let features = if strict {
        "sanitize,op-contracts-strict,cosmolkit-core/op-contracts-strict"
    } else {
        "sanitize"
    };
    Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()))
        .current_dir(workspace)
        .env("CARGO_TARGET_DIR", target)
        .env("CARGO_INCREMENTAL", "0")
        .env("RUSTFLAGS", rustflags)
        .args([
            "check",
            "--quiet",
            "-p",
            "cosmolkit",
            "--lib",
            "--features",
            features,
        ])
        .output()
        .expect("run the real sanitize operation-body compile-privacy probe")
}

#[test]
fn sanitize_operation_body_sees_only_its_generated_capabilities() {
    for strict in [false, true] {
        let mode = if strict { "strict" } else { "default" };
        let allowed = privacy_cargo_check("sanitize_allowed", strict);
        assert!(
            allowed.status.success(),
            "authorized sanitize capabilities did not compile in {mode} mode:\n{}",
            String::from_utf8_lossy(&allowed.stderr)
        );

        let forbidden = privacy_cargo_check("sanitize_forbidden", strict);
        assert!(
            !forbidden.status.success(),
            "unauthorized sanitize capabilities unexpectedly compiled in {mode} mode"
        );
        let stderr = String::from_utf8_lossy(&forbidden.stderr);
        for rejected_surface in [
            "coordinates",
            "read_topology_runtime",
            "read_properties_runtime",
            "checkout_topology_runtime",
            "install_topology_runtime",
            "spec",
            "source",
            "derived_cache",
            "in_place_target",
            "new_in_place",
            "finish",
            "abort_in_place",
            "finish_in_place",
        ] {
            assert!(
                stderr.contains(rejected_surface),
                "{mode} compile failure did not prove `{rejected_surface}` is hidden:\n{stderr}"
            );
        }
    }
}

#[test]
fn sanitize_has_no_mapping_multiple_output_or_inplace_alias() {
    let spec = operation_spec("sanitize_with_params").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "Molecule");
    assert_eq!(spec.topology_edit, TopologyEditKind::Local);
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.auto_remap, BlockSet::NONE);
    assert!(
        BINDING_CONTRACT
            .iter()
            .filter(|row| row.feature == "sanitize")
            .all(|row| !row.semantic_id.ends_with('_'))
    );
}
