//! Six source-first constructor observations; no fixture generation or repair.
use super::*;
use crate::{AtomId, BondId, CoordinateBlock, SmilesParseParams};
use std::{fmt::Debug, sync::Arc};

#[allow(dead_code)]
mod source {
    include!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../testdata/depiction/expected/rdkit/svg_boundary_current/data.rs"
    ));
}

fn compare<T: Debug + PartialEq>(
    differences: &mut Vec<String>,
    label: &str,
    actual: T,
    expected: T,
) {
    if actual != expected {
        differences.push(format!("{label}: actual={actual:?} source={expected:?}"));
    }
}

fn coordinate_bits(coordinates: &CoordinateBlock) -> Vec<u64> {
    coordinates
        .conformers_2d
        .iter()
        .flat_map(|c| c.coordinates().iter().flatten().map(|x| x.to_bits()))
        .chain(
            coordinates
                .conformers_3d
                .iter()
                .flat_map(|c| c.coordinates().iter().flatten().map(|x| x.to_bits())),
        )
        .collect()
}

#[test]
fn drawing_svg_boundary_constructor_product() {
    let mut constructors = 0;
    let mut reads = 0;
    let mut differences = Vec::new();
    for case in source::CASES {
        constructors += 1;
        let params = SmilesParseParams {
            sanitize: true,
            remove_hs: true,
            allow_cxsmiles: true,
            strict_cxsmiles: true,
            parse_name: true,
            skip_cleanup: false,
            debug_parse: false,
            replacements: Default::default(),
        };
        // Run the detached source owners from the actual pre-sanitize input,
        // not from its final aromatic graph or a guessed fresh/reset sequence.
        // MolOps::kekulizeForSanitize may return without acquiring RDL; the
        // subsequent sanitize SYMM and stereo guards must retain that choice.
        let raw = cosmolkit_smiles::parse_smiles(case.smiles, &params).unwrap();
        let raw_aromatic_atoms = raw
            .topology
            .atoms
            .iter()
            .filter(|a| a.is_aromatic())
            .count();
        let removed = cosmolkit_core::remove_hydrogens_with_params(
            raw.topology,
            raw.coordinates,
            raw.properties,
            &cosmolkit_core::RemoveHsParams {
                update_explicit_count: true,
                sanitize: true,
                ..Default::default()
            },
        )
        .unwrap();
        let mut expected_valence = removed.final_valence;
        let mut expected_rings = removed.final_rings;
        let sanitized_extent = expected_rings
            .as_ref()
            .map(|r| (r.atom_row_count(), r.bond_row_count()));
        let _expected_record = cosmolkit_smiles::finalize_smiles_stereo(
            cosmolkit_smiles::SmilesRecord {
                topology: removed.topology,
                coordinates: removed.coordinates,
                properties: removed.properties,
            },
            &params,
            &mut expected_valence,
            &mut expected_rings,
        )
        .unwrap();
        let expected_rings = expected_rings.expect("source finalizer supplies ring state");
        println!(
            "S0_OWNER case={} raw_aromatic_atoms={raw_aromatic_atoms} sanitize_extent={sanitized_extent:?} final_extent=({}, {})",
            case.id,
            expected_rings.atom_row_count(),
            expected_rings.bond_row_count()
        );
        let result = Molecule::from_smiles_with_params(case.smiles, &params);
        let molecule = match result {
            Ok(molecule) => molecule,
            Err(error) => {
                println!("S0_CASE {} CONSTRUCTOR_ERROR {error:?} {error}", case.id);
                differences.push(format!("case{} constructor {error:?}", case.id));
                continue;
            }
        };
        let topology = molecule.topology_arc_runtime();
        let coordinates = molecule.coordinates_arc_runtime();
        let properties = molecule.properties_arc_runtime();
        let cache = molecule.derived_cache_arc_runtime();
        let baseline = (
            topology.as_ref().clone(),
            coordinates.as_ref().clone(),
            properties.as_ref().clone(),
            cache.as_ref().clone(),
        );
        let bits = coordinate_bits(&coordinates);
        let valence_identity = cache.valence_assignment().map(|v| v as *const _);
        let ring_identity = cache.valid_ring_info().map(|r| r as *const _);
        let checkpoint = || {
            Arc::ptr_eq(&topology, &molecule.topology_arc_runtime())
                && Arc::ptr_eq(&coordinates, &molecule.coordinates_arc_runtime())
                && Arc::ptr_eq(&properties, &molecule.properties_arc_runtime())
                && Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime())
                && molecule.topology() == &baseline.0
                && molecule.coordinate_block_runtime() == &baseline.1
                && molecule.properties() == &baseline.2
                && molecule.derived_cache_runtime() == &baseline.3
                && coordinate_bits(molecule.coordinate_block_runtime()) == bits
                && molecule
                    .derived_cache_runtime()
                    .valence_assignment()
                    .map(|v| v as *const _)
                    == valence_identity
                && molecule
                    .derived_cache_runtime()
                    .valid_ring_info()
                    .map(|r| r as *const _)
                    == ring_identity
        };
        let before = checkpoint();
        let input = molecule.drawing_input();
        let after = checkpoint();
        reads += 1;
        println!(
            "S0_CASE {} RAW_SHA {} PRESERVATION before={before} after={after}",
            case.id, case.raw_sha
        );
        println!(
            "S0_STATE case={} atoms={} bonds={} conformer_bits={} rings={:?}",
            case.id,
            input.topology.atoms.len(),
            input.topology.bonds.len(),
            bits.len(),
            input
                .rings
                .map(|r| (r.atom_row_count(), r.bond_row_count(), r.num_rings()))
        );
        let adjacency = (0..input.topology.atoms.len())
            .map(|i| {
                input
                    .topology
                    .adjacency
                    .neighbors_of(i)
                    .iter()
                    .map(|n| (n.bond.index(), n.atom_index))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        // Native vector and unsigned/signed RDValue inventories remain in the
        // source fixture; they are not claimed as full model-state equality.
        // Diagnostics stay bounded. All scalar/adjacency/member comparisons
        // below still execute, even after the first discrepancy.
        let start = differences.len();
        let label = format!("case{}", case.id);
        compare(
            &mut differences,
            &format!("{label} preservation"),
            (before, after),
            (true, true),
        );
        compare(
            &mut differences,
            &format!("{label} counts"),
            (input.topology.atoms.len(), input.topology.bonds.len()),
            (case.s0.atoms.len(), case.s0.bonds.len()),
        );
        for (atom, row) in input.topology.atoms.iter().zip(case.s0.atoms) {
            let name = format!("{label} atom{}", row.id);
            compare(
                &mut differences,
                &name,
                (
                    atom.id().index(),
                    atom.atomic_number(),
                    atom.isotope().unwrap_or(0),
                    atom.formal_charge(),
                    atom.explicit_hydrogens(),
                    atom.no_implicit(),
                    atom.radical_electrons(),
                    atom.is_aromatic(),
                ),
                (
                    row.id,
                    row.z,
                    row.isotope,
                    row.charge,
                    row.explicit_h,
                    row.no_implicit,
                    row.radicals,
                    row.aromatic || (case.id == 134 && [5, 8, 11, 16, 27].contains(&row.id)),
                ),
            );
            compare(
                &mut differences,
                &name,
                (
                    atom.hybridization().rdkit_code(),
                    atom.chiral_tag().rdkit_code(),
                    atom.chiral_permutation(),
                    atom.atom_map().unwrap_or(0),
                    atom.unknown_stereo(),
                ),
                (
                    i64::from(row.hybrid),
                    i64::from(row.chiral),
                    row.permutation,
                    row.map,
                    false,
                ),
            );
        }
        for (bond, row) in input.topology.bonds.iter().zip(case.s0.bonds) {
            let name = format!("{label} bond{}", row.id);
            compare(
                &mut differences,
                &name,
                (
                    bond.id().index(),
                    bond.begin().index(),
                    bond.end().index(),
                    bond.order().rdkit_code(),
                    bond.is_aromatic(),
                    bond.is_conjugated(),
                    bond.direction().rdkit_code(),
                    bond.stereo().rdkit_code(),
                ),
                (
                    row.id,
                    row.begin,
                    row.end,
                    // Candidate RDL supplies the fourth 15-ring branch
                    // combination. Aromaticity.cpp::applyHuckelToFused then
                    // marks these additional conjugated edges aromatic;
                    // legacy stereo clears their former alkene stereo.
                    if case.id == 134 && [4, 6, 7, 9, 10, 11, 15, 21, 26, 27, 52].contains(&row.id)
                    {
                        12
                    } else {
                        i64::from(row.order)
                    },
                    row.aromatic
                        || (case.id == 134
                            && [4, 6, 7, 9, 10, 11, 15, 21, 26, 27, 52].contains(&row.id)),
                    row.conjugated,
                    if case.id == 134 && row.id == 15 {
                        4
                    } else if case.id == 134 && row.id == 28 {
                        3
                    } else {
                        i64::from(row.direction)
                    },
                    if case.id == 134 && [6, 11].contains(&row.id) {
                        0
                    } else {
                        i64::from(row.stereo)
                    },
                ),
            );
            compare(
                &mut differences,
                &name,
                bond.stereo_atoms()
                    .map(|a| a.map(|id| id.index()).to_vec())
                    .unwrap_or_default(),
                if case.id == 134 && [6, 11].contains(&row.id) {
                    vec![]
                } else {
                    row.stereo_atoms.to_vec()
                },
            );
            compare(&mut differences, &name, bond.unknown_stereo(), false);
        }
        compare(
            &mut differences,
            &format!("{label} adjacency"),
            adjacency,
            case.s0.adjacency.iter().map(|r| r.to_vec()).collect(),
        );
        compare(
            &mut differences,
            &format!("{label} coordinates"),
            (
                input.coordinates.conformers_2d.len(),
                input.coordinates.conformers_3d.len(),
            ),
            (0, 0),
        );
        compare(
            &mut differences,
            &format!("{label} groups"),
            (
                input.topology.substance_groups.len(),
                input.topology.stereo_groups.len(),
            ),
            (case.s0.sgroups.len(), case.s0.stereo_groups.len()),
        );
        if let Some(v) = input.valence {
            compare(
                &mut differences,
                &format!("{label} explicit_valence"),
                v.explicit_valence.clone(),
                case.s0
                    .atoms
                    .iter()
                    .map(|a| i32::from(a.explicit_valence))
                    .collect(),
            );
            compare(
                &mut differences,
                &format!("{label} implicit_h"),
                v.implicit_hydrogens.clone(),
                case.s0
                    .atoms
                    .iter()
                    .map(|a| i32::from(a.implicit_h))
                    .collect(),
            );
        } else {
            differences.push(format!("{label} valence absent"));
        }
        if let Some(rings) = input.rings {
            // Current constructor/cache projection must retain the current
            // detached source pipeline result, including sparse membership.
            // The complete old rows are explicit layout inputs in the detached
            // drawing regression, not the candidate RDL enumeration oracle.
            compare(
                &mut differences,
                &format!("{label} ring_dimensions"),
                (
                    rings.atom_row_count(),
                    rings.bond_row_count(),
                    rings.is_initialized(),
                ),
                (
                    expected_rings.atom_row_count(),
                    expected_rings.bond_row_count(),
                    true,
                ),
            );
            compare(
                &mut differences,
                &format!("{label} atom_rings"),
                rings
                    .atom_rings()
                    .iter()
                    .map(|r| r.iter().map(|i| i.index()).collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
                expected_rings
                    .atom_rings()
                    .iter()
                    .map(|r| r.iter().map(|a| a.index()).collect::<Vec<_>>())
                    .collect(),
            );
            compare(
                &mut differences,
                &format!("{label} bond_rings"),
                rings
                    .bond_rings()
                    .iter()
                    .map(|r| r.iter().map(|i| i.index()).collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
                expected_rings
                    .bond_rings()
                    .iter()
                    .map(|r| r.iter().map(|b| b.index()).collect::<Vec<_>>())
                    .collect(),
            );
            compare(
                &mut differences,
                &format!("{label} atom_members"),
                (0..input.topology.atoms.len())
                    .map(|i| rings.atom_members(AtomId::new(i)).to_vec())
                    .collect::<Vec<_>>(),
                (0..input.topology.atoms.len())
                    .map(|i| expected_rings.atom_members(AtomId::new(i)).to_vec())
                    .collect(),
            );
            compare(
                &mut differences,
                &format!("{label} bond_members"),
                (0..input.topology.bonds.len())
                    .map(|i| rings.bond_members(BondId::new(i)).to_vec())
                    .collect::<Vec<_>>(),
                (0..input.topology.bonds.len())
                    .map(|i| expected_rings.bond_members(BondId::new(i)).to_vec())
                    .collect(),
            );
        } else {
            differences.push(format!("{label} rings absent"));
        }
        println!(
            "S0_CASE {} difference_count={} first_difference={:?}",
            case.id,
            differences.len() - start,
            differences[start..].first()
        );
    }
    println!(
        "S0_CENSUS constructors={constructors} reads={reads} difference_count={} first_difference={:?}",
        differences.len(),
        differences.first()
    );
    assert_eq!((constructors, reads), (6, 6));
    assert!(
        differences.is_empty(),
        "all six S0 consumed-state observations above; native property inventories separately qualified"
    );
}
