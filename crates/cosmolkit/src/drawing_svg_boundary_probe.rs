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
        let result = Molecule::from_smiles_with_params(
            case.smiles,
            &SmilesParseParams {
                sanitize: true,
                remove_hydrogens: true,
                allow_cxsmiles: true,
                strict_cxsmiles: true,
                parse_name: true,
                skip_cleanup: false,
                debug_parse: false,
                replacements: Default::default(),
            },
        );
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
            "S0_TOPOLOGY {:#?}\nS0_COORDINATES {:#?}\nS0_PROPERTIES {:#?}\nS0_CACHE {:#?}\nS0_VALENCE {:#?}\nS0_RINGS {:#?}\nS0_BITS {bits:?}",
            input.topology, input.coordinates, input.properties, cache, input.valence, input.rings
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
        println!("S0_ADJACENCY {adjacency:?}");
        // Complete source-only inventories remain observable. Native vectors and
        // unsigned/signed RDValue distinctions are not full model-state equality.
        println!(
            "S0_NATIVE_PROPERTIES atoms={:?} bonds={:?} molecule={:?}",
            case.s0.atoms.iter().map(|a| a.props).collect::<Vec<_>>(),
            case.s0.bonds.iter().map(|b| b.props).collect::<Vec<_>>(),
            case.s0.props
        );
        for atom in &input.topology.atoms {
            println!(
                "S0_ATOM {} hybrid={} chiral={} props={:?} computed={:?}",
                atom.id().index(),
                atom.hybridization().rdkit_code(),
                atom.chiral_tag().rdkit_code(),
                atom.props(),
                atom.computed_prop_names()
            );
        }
        for bond in &input.topology.bonds {
            println!(
                "S0_BOND {} order={} direction={} stereo={} props={:?} computed={:?}",
                bond.id().index(),
                bond.order().rdkit_code(),
                bond.direction().rdkit_code(),
                bond.stereo().rdkit_code(),
                bond.props(),
                bond.computed_prop_names()
            );
        }
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
                    row.aromatic,
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
                    i64::from(row.order),
                    row.aromatic,
                    row.conjugated,
                    i64::from(row.direction),
                    i64::from(row.stereo),
                ),
            );
            compare(
                &mut differences,
                &name,
                bond.stereo_atoms()
                    .map(|a| a.map(|id| id.index()).to_vec())
                    .unwrap_or_default(),
                row.stereo_atoms.to_vec(),
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
            println!(
                "S0_RING_MEMBERS atoms={:?} bonds={:?}",
                (0..input.topology.atoms.len())
                    .map(|i| rings.atom_members(AtomId::new(i)))
                    .collect::<Vec<_>>(),
                (0..input.topology.bonds.len())
                    .map(|i| rings.bond_members(BondId::new(i)))
                    .collect::<Vec<_>>()
            );
            compare(
                &mut differences,
                &format!("{label} ring_dimensions"),
                (
                    rings.atom_row_count(),
                    rings.bond_row_count(),
                    rings.is_initialized(),
                ),
                (case.s0.atoms.len(), case.s0.bonds.len(), true),
            );
            compare(
                &mut differences,
                &format!("{label} atom_rings"),
                rings
                    .atom_rings()
                    .iter()
                    .map(|r| r.iter().map(|i| i.index()).collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
                case.s0.rings.atom_rows.iter().map(|r| r.to_vec()).collect(),
            );
            compare(
                &mut differences,
                &format!("{label} bond_rings"),
                rings
                    .bond_rings()
                    .iter()
                    .map(|r| r.iter().map(|i| i.index()).collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
                case.s0.rings.bond_rows.iter().map(|r| r.to_vec()).collect(),
            );
            compare(
                &mut differences,
                &format!("{label} atom_members"),
                (0..input.topology.atoms.len())
                    .map(|i| rings.atom_members(AtomId::new(i)).to_vec())
                    .collect::<Vec<_>>(),
                case.s0
                    .rings
                    .atom_members
                    .iter()
                    .map(|r| r.to_vec())
                    .collect(),
            );
            compare(
                &mut differences,
                &format!("{label} bond_members"),
                (0..input.topology.bonds.len())
                    .map(|i| rings.bond_members(BondId::new(i)).to_vec())
                    .collect::<Vec<_>>(),
                case.s0
                    .rings
                    .bond_members
                    .iter()
                    .map(|r| r.to_vec())
                    .collect(),
            );
        } else {
            differences.push(format!("{label} rings absent"));
        }
        println!(
            "S0_CASE {} CONSUMED_DIFFERENCES {:?}",
            case.id,
            &differences[start..]
        );
    }
    println!("S0_CENSUS constructors={constructors} reads={reads} differences={differences:#?}");
    assert_eq!((constructors, reads), (6, 6));
    assert!(
        differences.is_empty(),
        "all six S0 consumed-state observations above; native property inventories separately qualified"
    );
}
