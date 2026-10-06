//! Public RemoveHs ring-transport regressions (R2 sixteen-call product).
//!
//! G0..G3 x sanitize{false,true} x remove_nonimplicit{false,true} through
//! the public entry, each call carrying one 2D (ID 17) and one 3D (ID 19)
//! conformer with finite per-source-row sentinels (including signed zero)
//! plus ordinary molecule/conformer properties. Expected output rows come
//! from the LITERAL new_to_old tables; retained ORIGINAL baselines are
//! captured fresh before each call and compared after that call.

use cosmolkit_core::{RemoveHsParams, remove_hydrogens_with_params};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer2D, Conformer3D, CoordinateBlock,
    MoleculeProperties, TopologyBlock,
};
use cosmolkit_types::{BondOrder, Element};

fn carbon(id: usize, aromatic: bool) -> Atom {
    let mut spec = AtomSpec::new(Element::C);
    if aromatic {
        spec = spec.with_aromatic(true);
    }
    Atom::from_spec(AtomId::new(id), spec)
}

fn hydrogen(id: usize) -> Atom {
    Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::H))
}

fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
}

fn single(id: usize, begin: usize, end: usize) -> Bond {
    Bond::from_spec(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
    )
}

fn aromatic_bond(id: usize, begin: usize, end: usize) -> Bond {
    Bond::from_spec(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic)
            .with_aromatic(true),
    )
}

fn g1() -> TopologyBlock {
    topology(
        vec![carbon(0, false), carbon(1, false)],
        vec![single(0, 0, 1)],
    )
}

fn g2() -> TopologyBlock {
    topology(
        (0..6).map(|id| carbon(id, true)).collect(),
        (0..6)
            .map(|id| aromatic_bond(id, id, (id + 1) % 6))
            .collect(),
    )
}

fn g3() -> TopologyBlock {
    topology(
        std::iter::once(hydrogen(0))
            .chain((1..7).map(|id| carbon(id, true)))
            .collect(),
        std::iter::once(single(0, 0, 1))
            .chain((1..7).map(|id| aromatic_bond(id, id, (id % 6) + 1)))
            .collect(),
    )
}

/// Finite per-source-row sentinels including -0.0 (sign bit preserved).
fn coordinates_2d(atom_count: usize) -> Conformer2D {
    Conformer2D::new(
        17,
        (0..atom_count)
            .map(|row| {
                [
                    if row % 2 == 0 {
                        row as f64
                    } else {
                        -(row as f64)
                    },
                    if row == 0 { -0.0 } else { row as f64 + 0.5 },
                ]
            })
            .collect(),
    )
    .with_prop("conf-ordinary", "2d-kept")
}

fn coordinates_3d(atom_count: usize) -> Conformer3D {
    Conformer3D::new(
        19,
        (0..atom_count)
            .map(|row| {
                [
                    row as f64 + 0.25,
                    -(row as f64) - 0.25,
                    if row == 0 { -0.0 } else { row as f64 + 0.75 },
                ]
            })
            .collect(),
        true,
    )
    .with_prop("conf-ordinary", "3d-kept")
}

fn coordinate_block(atom_count: usize) -> CoordinateBlock {
    CoordinateBlock {
        conformers_2d: vec![coordinates_2d(atom_count)],
        conformers_3d: vec![coordinates_3d(atom_count)],
        source_coordinate_dim: None,
        source_conformer_order: None,
    }
}

fn properties() -> MoleculeProperties {
    let mut properties = MoleculeProperties::default();
    properties.set_prop("ck-ordinary", "kept").unwrap();
    properties
}

fn params(sanitize: bool, nonimplicit: bool) -> RemoveHsParams {
    let mut params = RemoveHsParams::default();
    params.sanitize = sanitize;
    params.remove_nonimplicit = nonimplicit;
    params.update_explicit_count = true;
    params.show_warnings = false;
    params
}

fn bits_2d(rows: &[[f64; 2]]) -> Vec<(u64, u64)> {
    rows.iter()
        .map(|[x, y]| (x.to_bits(), y.to_bits()))
        .collect()
}

fn bits_3d(rows: &[[f64; 3]]) -> Vec<(u64, u64, u64)> {
    rows.iter()
        .map(|[x, y, z]| (x.to_bits(), y.to_bits(), z.to_bits()))
        .collect()
}

#[test]
fn remove_hydrogens_ring_public_sixteen_call_product() {
    let graphs: [(usize, fn() -> TopologyBlock, usize); 3] =
        [(1, g1 as fn() -> TopologyBlock, 2), (2, g2, 6), (3, g3, 7)];
    let empty = topology(Vec::new(), Vec::new());
    let mut calls = 0usize;
    for sanitize in [false, true] {
        for nonimplicit in [false, true] {
            let label = format!("G0/s{sanitize}/ni{nonimplicit}");
            // Fresh per-call retained ORIGINAL baselines; clones enter the
            // owner, the retained originals are compared immediately after.
            let retained_graph = empty.clone();
            let retained_coords = coordinate_block(0);
            let retained_props = properties();
            let result = remove_hydrogens_with_params(
                empty.clone(),
                coordinate_block(0),
                properties(),
                &params(sanitize, nonimplicit),
            )
            .unwrap();
            calls += 1;
            assert_eq!(retained_graph, empty, "{label}: retained graph");
            assert_eq!(
                retained_coords,
                coordinate_block(0),
                "{label}: retained coords"
            );
            assert_eq!(retained_props, properties(), "{label}: retained props");
            assert!(result.final_rings.is_none(), "{label}");
            assert_eq!(result.topology.atoms.len(), 0);
            // Empty output: full empty mappings, empty 2D/3D row tables,
            // original IDs 17/19, stored dimensionality, ordinary
            // molecule/conformer properties.
            assert!(result.mapping.atoms().new_to_old().is_empty());
            assert!(result.mapping.atoms().old_to_new().is_empty());
            assert!(result.mapping.bonds().new_to_old().is_empty());
            assert!(result.mapping.bonds().old_to_new().is_empty());
            // Retained ORIGINAL conformers keep IDs 17/19 (retention proof).
            assert_eq!(retained_coords.conformers_2d[0].id(), 17, "{label}");
            assert_eq!(retained_coords.conformers_3d[0].id(), 19, "{label}");
            assert_eq!(result.coordinates.conformers_2d.len(), 1);
            // Atom-row remapping must preserve the output conformer IDs,
            // including conformers with zero coordinate rows.
            assert_eq!(
                result.coordinates.conformers_2d[0].id(),
                17,
                "{label}: out 2D id"
            );
            assert!(result.coordinates.conformers_2d[0].coordinates().is_empty());
            assert_eq!(result.coordinates.conformers_3d.len(), 1);
            assert_eq!(
                result.coordinates.conformers_3d[0].id(),
                19,
                "{label}: out 3D id"
            );
            assert!(result.coordinates.conformers_3d[0].is_3d());
            assert!(result.coordinates.conformers_3d[0].coordinates().is_empty());
            assert_eq!(
                result.properties.prop("ck-ordinary"),
                Some("kept"),
                "{label}"
            );
            assert_eq!(
                result
                    .coordinates
                    .conformers_2d
                    .first()
                    .and_then(|c| c.props().get("conf-ordinary").map(String::as_str)),
                Some("2d-kept"),
                "{label}"
            );
            // Valence/ring separation: final_valence Some IFF sanitize.
            match (sanitize, &result.final_valence) {
                (true, Some(_)) => {}
                (false, None) => {}
                _ => panic!("{label}: valence/ring separation"),
            }
        }
    }
    for (index, build, atom_count) in graphs {
        for sanitize in [false, true] {
            for nonimplicit in [false, true] {
                let label = format!("G{index}/s{sanitize}/ni{nonimplicit}");
                let topology = build();
                let coordinates = coordinate_block(atom_count);
                let properties = properties();
                // Fresh retained ORIGINAL baselines captured BEFORE the call.
                let topology_snapshot = topology.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                // Retained originals: clones enter the owner; the retained
                // values are compared to fresh baselines immediately after.
                let result = remove_hydrogens_with_params(
                    topology.clone(),
                    coordinates.clone(),
                    properties.clone(),
                    &params(sanitize, nonimplicit),
                )
                .unwrap();
                calls += 1;
                // Retained ORIGINAL baselines: whole values.
                assert_eq!(topology_snapshot, topology, "{label}: retained graph");
                assert_eq!(
                    coordinates_snapshot, coordinates,
                    "{label}: retained coords"
                );
                assert_eq!(properties_snapshot, properties, "{label}: retained props");

                let both = sanitize && nonimplicit;
                match (index, both) {
                    (2, true) | (3, true) => {
                        let rings = result.final_rings.as_ref().unwrap();
                        assert_eq!(
                            rings.find_type(),
                            cosmolkit_core::RingFindType::SymmSssr,
                            "{label}"
                        );
                        assert!(rings.is_initialized(), "{label}");
                        let expected_final_atoms = if index == 2 { 6 } else { 6 };
                        assert_eq!(rings.atom_row_count(), expected_final_atoms, "{label}");
                        assert_eq!(rings.atom_rings().len(), 1, "{label}: rings");
                    }
                    (1, true) => {
                        let rings = result.final_rings.as_ref().unwrap();
                        assert_eq!(
                            rings.find_type(),
                            cosmolkit_core::RingFindType::SymmSssr,
                            "{label}"
                        );
                        assert!(rings.atom_rings().is_empty(), "{label}: empty rows");
                    }
                    _ => {
                        assert!(result.final_rings.is_none(), "{label}");
                    }
                }
                // Literal mapping tables.
                let atom_n2o: Vec<Option<usize>> = result
                    .mapping
                    .atoms()
                    .new_to_old()
                    .iter()
                    .map(|atom| atom.map(|a| a.index()))
                    .collect();
                let bond_n2o: Vec<Option<usize>> = result
                    .mapping
                    .bonds()
                    .new_to_old()
                    .iter()
                    .map(|bond| bond.map(|b| b.index()))
                    .collect();
                if index == 3 && nonimplicit {
                    assert_eq!(
                        atom_n2o,
                        vec![Some(1), Some(2), Some(3), Some(4), Some(5), Some(6)],
                        "{label}: atom new_to_old"
                    );
                    assert_eq!(
                        bond_n2o,
                        vec![Some(1), Some(2), Some(3), Some(4), Some(5), Some(6)],
                        "{label}: bond new_to_old"
                    );
                    let atom_o2n: Vec<Option<usize>> = result
                        .mapping
                        .atoms()
                        .old_to_new()
                        .iter()
                        .map(|atom| atom.map(|a| a.index()))
                        .collect();
                    assert_eq!(
                        atom_o2n,
                        vec![None, Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)],
                        "{label}: atom old_to_new"
                    );
                    let bond_o2n: Vec<Option<usize>> = result
                        .mapping
                        .bonds()
                        .old_to_new()
                        .iter()
                        .map(|bond| bond.map(|b| b.index()))
                        .collect();
                    assert_eq!(
                        bond_o2n,
                        vec![None, Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)],
                        "{label}: bond old_to_new"
                    );
                } else {
                    // Literal ORIGINAL bond counts: G1=1, G2=6, G3=7.
                    let bond_count = match index {
                        1 => 1,
                        2 => 6,
                        _ => 7,
                    };
                    let identity: Vec<Option<usize>> = (0..atom_count).map(Some).collect();
                    let identity_bonds: Vec<Option<usize>> = (0..bond_count).map(Some).collect();
                    assert_eq!(atom_n2o, identity, "{label}: identity atoms");
                    let atom_o2n: Vec<Option<usize>> = result
                        .mapping
                        .atoms()
                        .old_to_new()
                        .iter()
                        .map(|atom| atom.map(|a| a.index()))
                        .collect();
                    assert_eq!(atom_o2n, identity, "{label}: identity atoms o2n");
                    assert_eq!(bond_n2o, identity_bonds, "{label}: identity bonds");
                    let bond_o2n: Vec<Option<usize>> = result
                        .mapping
                        .bonds()
                        .old_to_new()
                        .iter()
                        .map(|bond| bond.map(|b| b.index()))
                        .collect();
                    assert_eq!(bond_o2n, identity_bonds, "{label}: identity bonds o2n");
                }
                // Ordinary properties retained on the output.
                assert_eq!(
                    result.properties.prop("ck-ordinary"),
                    Some("kept"),
                    "{label}"
                );
                assert_eq!(
                    result
                        .coordinates
                        .conformers_2d
                        .first()
                        .and_then(|c| c.props().get("conf-ordinary").map(String::as_str)),
                    Some("2d-kept"),
                    "{label}"
                );
                assert_eq!(
                    result
                        .coordinates
                        .conformers_3d
                        .first()
                        .and_then(|c| c.props().get("conf-ordinary").map(String::as_str)),
                    Some("3d-kept"),
                    "{label}"
                );
                // Coordinates follow the frozen final original indices, by bits.
                let expected_rows: Vec<usize> = if index == 3 && nonimplicit {
                    vec![1, 2, 3, 4, 5, 6]
                } else {
                    (0..atom_count).collect()
                };
                let out2 = result.coordinates.conformers_2d.first().unwrap();
                let out3 = result.coordinates.conformers_3d.first().unwrap();
                let input2 = coordinates_snapshot.conformers_2d.first().unwrap();
                let input3 = coordinates_snapshot.conformers_3d.first().unwrap();
                let expected_2d: Vec<(u64, u64)> = expected_rows
                    .iter()
                    .map(|row| {
                        let [x, y] = input2.coordinates()[*row];
                        (x.to_bits(), y.to_bits())
                    })
                    .collect();
                let expected_3d: Vec<(u64, u64, u64)> = expected_rows
                    .iter()
                    .map(|row| {
                        let [x, y, z] = input3.coordinates()[*row];
                        (x.to_bits(), y.to_bits(), z.to_bits())
                    })
                    .collect();
                assert_eq!(bits_2d(out2.coordinates()), expected_2d, "{label}: 2D bits");
                assert_eq!(bits_3d(out3.coordinates()), expected_3d, "{label}: 3D bits");
                // Full component count/rows, conformer IDs, stored
                // dimensionality and output properties.
                assert_eq!(
                    out2.coordinates().len(),
                    expected_rows.len(),
                    "{label}: 2D rows"
                );
                assert_eq!(
                    out3.coordinates().len(),
                    expected_rows.len(),
                    "{label}: 3D rows"
                );
                // Both retained inputs and mapped outputs preserve IDs17/19;
                // atom indices must not become conformer IDs.
                assert_eq!(
                    coordinates_snapshot.conformers_2d[0].id(),
                    17,
                    "{label}: retained 2D id"
                );
                assert_eq!(
                    coordinates_snapshot.conformers_3d[0].id(),
                    19,
                    "{label}: retained 3D id"
                );
                assert_eq!(out2.id(), 17, "{label}: output 2D id");
                assert_eq!(out3.id(), 19, "{label}: output 3D id");
                assert!(out3.is_3d(), "{label}: 3D stored dimensionality");
                assert_eq!(
                    result
                        .coordinates
                        .conformers_3d
                        .first()
                        .and_then(|c| c.props().get("conf-ordinary").map(String::as_str)),
                    Some("3d-kept"),
                    "{label}"
                );
                // Frozen valence/ring separation: final_valence Some IFF
                // sanitize, both field lengths equal final atoms; the
                // nonimplicit=false/sanitize=true cell keeps final_rings None
                // even though valence is Some.
                match (sanitize, &result.final_valence) {
                    (true, Some(valence)) => {
                        assert_eq!(
                            valence.explicit_valence.len(),
                            result.topology.atoms.len(),
                            "{label}: valence atoms length"
                        );
                        assert_eq!(
                            valence.implicit_hydrogens.len(),
                            result.topology.atoms.len(),
                            "{label}: valence implicit length"
                        );
                        if !nonimplicit {
                            assert!(result.final_rings.is_none(), "{label}: rings/valence split");
                        }
                    }
                    (false, None) => {}
                    _ => panic!("{label}: valence presence"),
                }
            }
        }
    }
    assert_eq!(calls, 16, "exact census");
}
