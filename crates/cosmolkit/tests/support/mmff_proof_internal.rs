//! Local MMFF preservation-proof controls and concrete failure diagnostics.
use super::*;
use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};
use std::sync::Arc;

fn source() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let first = builder.add_atom(AtomSpec::new(Element::C));
    let second = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(first, second, BondOrder::Single))
        .unwrap();
    builder
        .add_3d_conformer(vec![[0., 0., 0.], [2., 0., 0.]])
        .unwrap();
    builder.build().unwrap().with_assigned_rings().unwrap()
}
fn parts(source: &Molecule) -> OpParts<'_, crate::ops::WithMmffOptimizedAccess> {
    let spec = crate::MOLECULE_OPS
        .iter()
        .find(|spec| spec.method == "with_mmff_optimized_with_params")
        .unwrap();
    let mut parts = OpParts::new(source, spec).unwrap();
    let mut properties = parts.checkout_properties_runtime().unwrap();
    properties
        .set_computed_prop("_MMFFSanitized", 1_i32)
        .unwrap();
    parts.install_properties_runtime(properties).unwrap();
    parts
}
fn preserved() -> DerivedState {
    DerivedState::RING_FAMILIES.union(DerivedState::COORDINATES)
}
fn prove(
    parts: &mut OpParts<'_, crate::ops::WithMmffOptimizedAccess>,
) -> Result<(), OperationError> {
    parts.prove_preserved_runtime(preserved(), PreservationProof::MmffPreparedOptimization)
}
fn assert_mismatch(
    parts: &mut OpParts<'_, crate::ops::WithMmffOptimizedAccess>,
    field: &str,
    expected: &str,
    actual: &str,
) {
    let error = prove(parts).unwrap_err();
    assert_eq!(
        error,
        OperationError::DerivedPreservationMismatch {
            operation: "with_mmff_optimized_with_params",
            states: preserved(),
            field: field.into(),
            expected: expected.into(),
            actual: actual.into(),
        }
    );
    assert_eq!(
        error.to_string(),
        format!(
            "operation `with_mmff_optimized_with_params` has invalid derived-effect preserve for mask 0x22: {field}: expected {expected}, got {actual}"
        )
    );
}
fn assert_source_unchanged(source: &Molecule, peer: &Molecule) {
    assert_eq!(source, peer);
    assert!(Arc::ptr_eq(
        &source.topology_arc_runtime(),
        &peer.topology_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.coordinates_arc_runtime(),
        &peer.coordinates_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.properties_arc_runtime(),
        &peer.properties_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.derived_cache_arc_runtime(),
        &peer.derived_cache_arc_runtime()
    ));
}
#[test]
fn mmff_proof_allows_preparation_guard_and_finite_3d_position_changes() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut coordinates = parts.checkout_coordinates_runtime().unwrap();
    coordinates.conformers_3d[0].coordinates_mut()[1][0] = 1.5;
    parts.install_coordinates_runtime(coordinates).unwrap();
    prove(&mut parts).unwrap();
    assert_source_unchanged(&source, &peer);
}
#[test]
fn mmff_proof_rejects_unrelated_atom_formal_charge() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut topology = parts.checkout_topology_runtime().unwrap();
    topology.atoms[0].set_formal_charge(1);
    parts.install_topology_runtime(topology).unwrap();
    assert_mismatch(&mut parts, "topology.atoms[0].formal_charge", "0", "1");
    assert_source_unchanged(&source, &peer);
}
#[test]
fn mmff_proof_rejects_conformer_identity_change() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut coordinates = parts.checkout_coordinates_runtime().unwrap();
    let old = coordinates.conformers_3d[0].clone();
    coordinates.conformers_3d[0] =
        crate::Conformer3D::new(17, old.coordinates().to_vec(), old.is_3d());
    parts.install_coordinates_runtime(coordinates).unwrap();
    assert_mismatch(&mut parts, "coordinates.conformers_3d[0].id", "0", "17");
    assert_source_unchanged(&source, &peer);
}
#[test]
fn mmff_proof_reports_bond_endpoint_change() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut topology = parts.checkout_topology_runtime().unwrap();
    let begin = topology.bonds[0].begin();
    let end = topology.bonds[0].end();
    topology.bonds[0].set_endpoints(end, begin);
    parts.install_topology_runtime(topology).unwrap();
    assert_mismatch(
        &mut parts,
        "topology.bonds[0].begin",
        "AtomId(0)",
        "AtomId(1)",
    );
    assert_source_unchanged(&source, &peer);
}
#[test]
fn mmff_proof_reports_atom_count_change() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let _topology = parts.checkout_topology_runtime().unwrap();
    parts
        .install_topology_runtime(cosmolkit_model::TopologyBlock::default())
        .unwrap();
    assert_mismatch(&mut parts, "topology.atoms.len", "2", "0");
    assert_source_unchanged(&source, &peer);
}
#[test]
fn mmff_proof_reports_bond_direction_change() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut topology = parts.checkout_topology_runtime().unwrap();
    topology.bonds[0].set_direction(crate::BondDirection::EndUpRight);
    parts.install_topology_runtime(topology).unwrap();
    assert_mismatch(
        &mut parts,
        "topology.bonds[0].direction",
        "None",
        "EndUpRight",
    );
    assert_source_unchanged(&source, &peer);
}
fn directed_source(smiles: &str, bond: usize, sanitized: bool) -> Molecule {
    let source = Molecule::from_smiles(smiles).unwrap();
    let mut topology = source.topology().clone();
    topology.bonds[bond].set_direction(crate::BondDirection::EndUpRight);
    let mut properties = source.properties().clone();
    if sanitized {
        properties
            .set_computed_prop("_MMFFSanitized", 1_i32)
            .unwrap();
    }
    Molecule::from_parts(
        topology,
        source.coordinate_block_runtime().clone(),
        properties,
    )
    .unwrap()
    .with_assigned_rings()
    .unwrap()
}
#[test]
fn mmff_proof_accepts_source_kekulize_direction_cleanup() {
    for order in [BondOrder::Single, BondOrder::Double, BondOrder::Aromatic] {
        let source = directed_source("c1ccccc1", 0, false);
        let peer = source.clone();
        let mut parts = parts(&source);
        let mut topology = parts.checkout_topology_runtime().unwrap();
        topology.bonds[0].set_direction(crate::BondDirection::None);
        topology.bonds[0].set_order(order);
        parts.install_topology_runtime(topology).unwrap();
        prove(&mut parts).unwrap();
        assert_source_unchanged(&source, &peer);
    }
}
#[test]
fn mmff_proof_rejects_direction_cleanup_without_source_preparation() {
    for (smiles, sanitized) in [("CC", false), ("c1ccccc1", true), ("Cc1ccccc1", false)] {
        let source = directed_source(smiles, 0, sanitized);
        let peer = source.clone();
        let mut parts = parts(&source);
        let mut topology = parts.checkout_topology_runtime().unwrap();
        topology.bonds[0].set_direction(crate::BondDirection::None);
        parts.install_topology_runtime(topology).unwrap();
        assert_mismatch(
            &mut parts,
            "topology.bonds[0].direction",
            "EndUpRight",
            "None",
        );
        assert_source_unchanged(&source, &peer);
    }
}
#[test]
fn mmff_proof_rejects_non_source_direction_replacement() {
    let source = directed_source("c1ccccc1", 0, false);
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut topology = parts.checkout_topology_runtime().unwrap();
    topology.bonds[0].set_direction(crate::BondDirection::EndDownRight);
    parts.install_topology_runtime(topology).unwrap();
    assert_mismatch(
        &mut parts,
        "topology.bonds[0].direction",
        "EndUpRight",
        "EndDownRight",
    );
    assert_source_unchanged(&source, &peer);
}
#[test]
fn mmff_proof_keeps_source_permitted_non_finite_positions() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut coordinates = parts.checkout_coordinates_runtime().unwrap();
    coordinates.conformers_3d[0].coordinates_mut()[1][2] = f64::NAN;
    // Stored positions need not be finite; numerical validation belongs to
    // the algorithm, not this metadata-only preservation proof.
    coordinates.validate_for_atom_count(2).unwrap();
    parts.install_coordinates_runtime(coordinates).unwrap();
    prove(&mut parts).unwrap();
    assert_source_unchanged(&source, &peer);
}
#[test]
fn mmff_proof_rejects_unrelated_property_write() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut properties = parts.checkout_properties_runtime().unwrap();
    properties.set_prop("unrelated", "written").unwrap();
    parts.install_properties_runtime(properties).unwrap();
    assert_mismatch(
        &mut parts,
        "properties.props[\"unrelated\"]",
        "None",
        &format!(
            "Some({:?})",
            cosmolkit_model::PropertyValue::String("written".into())
        ),
    );
    assert_source_unchanged(&source, &peer);
}
#[test]
fn mmff_proof_rejects_missing_preparation_guard() {
    let source = source();
    let peer = source.clone();
    let spec = crate::MOLECULE_OPS
        .iter()
        .find(|spec| spec.method == "with_mmff_optimized_with_params")
        .unwrap();
    let mut parts = OpParts::<crate::ops::WithMmffOptimizedAccess>::new(&source, spec).unwrap();
    let error = prove(&mut parts).unwrap_err();
    assert!(
        matches!(&error, OperationError::DerivedPreservationMismatch { field, expected, actual, .. }
        if field.contains("_MMFFSanitized") && expected == "Some(Int(1))" && actual == "None"),
        "{error}"
    );
    assert_source_unchanged(&source, &peer);
}
#[test]
fn mmff_proof_rejects_preserved_cache_validity_change() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut cache = parts.checkout_derived_cache_runtime().unwrap();
    cache.mark_valid(DerivedState::RING_FAMILIES);
    parts.install_derived_cache_runtime(cache).unwrap();
    let error = prove(&mut parts).unwrap_err();
    assert!(
        matches!(&error, OperationError::DerivedPreservationMismatch { field, .. }
        if field == "derived_cache.preserved_valid_states"),
        "{error}"
    );
    assert_source_unchanged(&source, &peer);
}
#[test]
fn mmff_proof_preserves_actual_signed_zero_coordinate_result() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut coordinates = parts.checkout_coordinates_runtime().unwrap();
    coordinates.conformers_3d[0].coordinates_mut()[0][0] = -0.0;
    parts.install_coordinates_runtime(coordinates).unwrap();
    prove(&mut parts).unwrap();
    assert_eq!(
        parts.current_coordinates_candidate().unwrap().conformers_3d[0].coordinates()[0][0]
            .to_bits(),
        (-0.0f64).to_bits()
    );
    assert!(matches!(parts.coordinates, WorkingBlock::Installed(_)));
    assert_source_unchanged(&source, &peer);
}

#[test]
fn mmff_proof_rejects_equal_text_with_non_source_scalar_tags() {
    for marker in [
        cosmolkit_model::PropertyValue::String("1".into()),
        cosmolkit_model::PropertyValue::UInt(1),
        cosmolkit_model::PropertyValue::Bool(true),
    ] {
        let source = source();
        let peer = source.clone();
        let mut parts = parts(&source);
        let mut properties = parts.checkout_properties_runtime().unwrap();
        properties
            .set_computed_prop("_MMFFSanitized", marker)
            .unwrap();
        parts.install_properties_runtime(properties).unwrap();
        let error = prove(&mut parts).unwrap_err();
        assert!(
            matches!(&error, OperationError::DerivedPreservationMismatch { field, expected, .. }
            if field.contains("_MMFFSanitized") && expected == "Some(Int(1))"),
            "{error}"
        );
        assert_source_unchanged(&source, &peer);
    }
}
#[test]
fn mmff_proof_rejects_coordinate_occurrence_metadata_changes() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut coordinates = parts.checkout_coordinates_runtime().unwrap();
    assert_eq!(
        coordinates.source_conformer_order,
        Some(vec![cosmolkit_model::CoordinateDimension::ThreeD])
    );
    // The canonical builder already records this append. Removing the known
    // occurrence order changes source metadata even for a single dimension.
    coordinates.source_conformer_order = None;
    parts.install_coordinates_runtime(coordinates).unwrap();
    assert_mismatch(
        &mut parts,
        "coordinates.source_conformer_order",
        "Some([ThreeD])",
        "None",
    );
    assert_source_unchanged(&source, &peer);
}

#[test]
fn mmff_proof_allows_identical_canonical_coordinate_occurrence_metadata() {
    let source = source();
    let peer = source.clone();
    let mut parts = parts(&source);
    let mut coordinates = parts.checkout_coordinates_runtime().unwrap();
    assert_eq!(
        coordinates.source_conformer_order,
        Some(vec![cosmolkit_model::CoordinateDimension::ThreeD])
    );
    coordinates.source_conformer_order = Some(vec![cosmolkit_model::CoordinateDimension::ThreeD]);
    parts.install_coordinates_runtime(coordinates).unwrap();
    prove(&mut parts).unwrap();
    assert!(matches!(parts.coordinates, WorkingBlock::Shared));
    assert_source_unchanged(&source, &peer);
}
