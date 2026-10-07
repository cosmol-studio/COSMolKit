//! Proposed controls for the new MMFF proof; existing normative controls unchanged.
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
    assert!(prove(&mut parts).is_err());
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
    assert!(prove(&mut parts).is_err());
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
    assert!(prove(&mut parts).is_err());
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
    assert!(prove(&mut parts).is_err());
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
    assert!(prove(&mut parts).is_err());
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
        assert!(prove(&mut parts).is_err());
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
    assert!(prove(&mut parts).is_err());
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
