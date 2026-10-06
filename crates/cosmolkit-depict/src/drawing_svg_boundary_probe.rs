//! Fixed source-literal preparation and independent renderer observations.
use super::*;
use crate::draw::render_prepared_svg;
use cosmolkit_core::RingFindType;
use cosmolkit_model::{
    Atom, AtomSpec, Bond, BondSpec, BondStereo, Conformer2D, Element, Hybridization, PropertyValue,
};
use std::{
    fmt::Debug,
    io::Write,
    process::{Command, Stdio},
};

#[allow(dead_code)]
mod source {
    include!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../testdata/depiction/expected/rdkit/svg_boundary_current/data.rs"
    ));
}

fn value(property: &source::NativeProperty) -> Result<Option<PropertyValue>, String> {
    Ok(match property.value {
        source::NativeValue::Int(n) => Some(PropertyValue::Int(
            i32::try_from(n).map_err(|e| e.to_string())?,
        )),
        source::NativeValue::Str(s) => Some(PropertyValue::String(s.to_owned())),
        // Native bookkeeping remains in the fixed source inventory, not a
        // model property. This is an explicit representation qualification.
        source::NativeValue::Strings(_) => None,
    })
}

fn literal(snapshot: &source::Snapshot) -> Result<PreparedDrawing, String> {
    if !snapshot.sgroups.is_empty() || !snapshot.stereo_groups.is_empty() {
        return Err("nonempty groups require existing typed projection".into());
    }
    let mut atoms = Vec::new();
    for row in snapshot.atoms {
        if row.query || row.permutation.is_some() || row.map != 0 {
            return Err(format!("unexpected fixed atom representation {}", row.id));
        }
        let mut spec = AtomSpec::new(Element::from_atomic_number(row.z).ok_or("element")?)
            .with_isotope(row.isotope)
            .with_formal_charge(row.charge)
            .with_explicit_hydrogens(row.explicit_h)
            .with_no_implicit(row.no_implicit)
            .with_radical_electrons(row.radicals)
            .with_aromatic(row.aromatic)
            .with_hybridization(
                Hybridization::from_rdkit_code(i64::from(row.hybrid)).ok_or("hybrid")?,
            )
            .with_chiral_tag(ChiralTag::from_rdkit_code(i64::from(row.chiral)).ok_or("chiral")?);
        for property in row.props {
            if let Some(value) = value(property)? {
                spec = if property.computed {
                    spec.with_computed_prop(property.key, value)
                } else {
                    spec.with_prop(property.key, value)
                }
                .map_err(|e| e.to_string())?;
            } else if property.key != "__computedProps" {
                return Err("unmodeled atom vector property".into());
            }
        }
        atoms.push(Atom::from_spec(AtomId::new(row.id), spec));
    }
    let mut bonds = Vec::new();
    for row in snapshot.bonds {
        if row.query {
            return Err("unexpected fixed bond query".into());
        }
        let mut spec = BondSpec::new(
            AtomId::new(row.begin),
            AtomId::new(row.end),
            BondOrder::from_rdkit_code(i64::from(row.order)).ok_or("order")?,
        )
        .with_aromatic(row.aromatic)
        .with_conjugated(row.conjugated)
        .with_direction(
            BondDirection::from_rdkit_code(i64::from(row.direction)).ok_or("direction")?,
        )
        .with_stereo(BondStereo::from_rdkit_code(i64::from(row.stereo)).ok_or("stereo")?);
        match row.stereo_atoms {
            [] => (),
            [a, b] => spec = spec.with_stereo_atoms(AtomId::new(*a), AtomId::new(*b)),
            _ => return Err("stereo atom count".into()),
        }
        for property in row.props {
            if let Some(value) = value(property)? {
                spec = if property.computed {
                    spec.with_computed_prop(property.key, value)
                } else {
                    spec.with_prop(property.key, value)
                }
                .map_err(|e| e.to_string())?;
            } else if property.key != "__computedProps" {
                return Err("unmodeled bond vector property".into());
            }
        }
        bonds.push(Bond::from_spec(BondId::new(row.id), spec));
    }
    let topology =
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).map_err(|e| e.to_string())?;
    if adjacency(&topology)
        != snapshot
            .adjacency
            .iter()
            .map(|r| r.to_vec())
            .collect::<Vec<_>>()
    {
        return Err("literal adjacency projection".into());
    }
    let mut properties = MoleculeProperties::default();
    for property in snapshot.props {
        let text = match property.value {
            source::NativeValue::Int(n) => n.to_string(),
            source::NativeValue::Str(s) => s.to_owned(),
            source::NativeValue::Strings(_) if property.key == "__computedProps" => continue,
            _ => return Err("unmodeled molecule vector property".into()),
        };
        properties = if property.computed {
            properties.with_computed_prop(property.key, text)
        } else {
            properties.with_prop(property.key, text)
        }
        .map_err(|e| e.to_string())?;
    }
    let mut coordinates = CoordinateBlock::default();
    for conformer in snapshot.conformers {
        if conformer.is3d
            || !conformer.props.is_empty()
            || conformer.xyz_bits.iter().any(|p| p[2] != 0)
        {
            return Err("fixed conformer projection prerequisites".into());
        }
        // The fixed source snapshot enumerates conformers in actual ROMol
        // insertion order. Preserve that fact in the detached projection;
        // coordinates, IDs, flags and all original comparisons stay intact.
        coordinates
            .record_source_conformer_append(cosmolkit_model::CoordinateDimension::TwoD)
            .map_err(|e| e.to_string())?;
        coordinates.conformers_2d.push(Conformer2D::new(
            conformer.id as usize,
            conformer
                .xyz_bits
                .iter()
                .map(|p| [f64::from_bits(p[0]), f64::from_bits(p[1])])
                .collect(),
        ));
    }
    coordinates
        .validate_for_atom_count(topology.atoms.len())
        .map_err(|e| e.to_string())?;
    let valence = ValenceAssignment {
        explicit_valence: snapshot
            .atoms
            .iter()
            .map(|a| i32::from(a.explicit_valence))
            .collect(),
        implicit_hydrogens: snapshot
            .atoms
            .iter()
            .map(|a| i32::from(a.implicit_h))
            .collect(),
    };
    // Modeled input quality is explicit; native find_type/cache is unexposed.
    let mut rings = RingInfo::new(
        RingFindType::SymmSssr,
        topology.atoms.len(),
        topology.bonds.len(),
    );
    if snapshot.rings.atom_rows.len() != snapshot.rings.bond_rows.len()
        || snapshot.rings.num_rings != snapshot.rings.atom_rows.len()
    {
        return Err("source ring counts".into());
    }
    for (i, (a, b)) in snapshot
        .rings
        .atom_rows
        .iter()
        .zip(snapshot.rings.bond_rows)
        .enumerate()
    {
        let count = rings.add_ring(a, b).map_err(|e| e.to_string())?;
        if count != i + 1 {
            return Err("append return count".into());
        }
    }
    let atom_members = (0..topology.atoms.len())
        .map(|i| rings.atom_members(AtomId::new(i)).to_vec())
        .collect::<Vec<_>>();
    let bond_members = (0..topology.bonds.len())
        .map(|i| rings.bond_members(BondId::new(i)).to_vec())
        .collect::<Vec<_>>();
    if rings
        .atom_rings()
        .iter()
        .map(|r| r.iter().map(|i| i.index()).collect::<Vec<_>>())
        .collect::<Vec<_>>()
        != snapshot
            .rings
            .atom_rows
            .iter()
            .map(|r| r.to_vec())
            .collect::<Vec<_>>()
        || rings
            .bond_rings()
            .iter()
            .map(|r| r.iter().map(|i| i.index()).collect::<Vec<_>>())
            .collect::<Vec<_>>()
            != snapshot
                .rings
                .bond_rows
                .iter()
                .map(|r| r.to_vec())
                .collect::<Vec<_>>()
        || atom_members
            != snapshot
                .rings
                .atom_members
                .iter()
                .map(|r| r.to_vec())
                .collect::<Vec<_>>()
        || bond_members
            != snapshot
                .rings
                .bond_members
                .iter()
                .map(|r| r.to_vec())
                .collect::<Vec<_>>()
        || !rings.is_initialized()
        || !rings.is_symm_sssr()
        || (rings.atom_row_count(), rings.bond_row_count())
            != (topology.atoms.len(), topology.bonds.len())
    {
        return Err("literal carrier rows/members/dimensions/quality".into());
    }
    Ok(PreparedDrawing {
        topology,
        coordinates,
        properties,
        valence,
        rings,
    })
}

fn adjacency(topology: &TopologyBlock) -> Vec<Vec<(usize, usize)>> {
    (0..topology.atoms.len())
        .map(|i| {
            topology
                .adjacency
                .neighbors_of(i)
                .iter()
                .map(|n| (n.bond.index(), n.atom_index))
                .collect()
        })
        .collect()
}

fn bits(coordinates: &CoordinateBlock) -> Vec<u64> {
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

fn observe(stage: &str, id: usize, input: &PreparedDrawing) {
    println!(
        "{stage}_CASE {id} TOPOLOGY {:#?}\nCOORDINATES {:#?}\nPROPERTIES {:#?}\nVALENCE {:#?}\nRINGS {:#?}",
        input.topology, input.coordinates, input.properties, input.valence, input.rings
    );
    for a in &input.topology.atoms {
        println!(
            "{stage}_ATOM {} hybrid={} chiral={} props={:?} computed={:?}",
            a.id().index(),
            a.hybridization().rdkit_code(),
            a.chiral_tag().rdkit_code(),
            a.props(),
            a.computed_prop_names()
        );
    }
    for b in &input.topology.bonds {
        println!(
            "{stage}_BOND {} order={} direction={} stereo={} props={:?} computed={:?}",
            b.id().index(),
            b.order().rdkit_code(),
            b.direction().rdkit_code(),
            b.stereo().rdkit_code(),
            b.props(),
            b.computed_prop_names()
        );
    }
    println!(
        "{stage}_ADJACENCY {:?}\n{stage}_BITS {:?}\n{stage}_RING_MEMBERS atoms={:?} bonds={:?}",
        adjacency(&input.topology),
        bits(&input.coordinates),
        (0..input.topology.atoms.len())
            .map(|i| input.rings.atom_members(AtomId::new(i)))
            .collect::<Vec<_>>(),
        (0..input.topology.bonds.len())
            .map(|i| input.rings.bond_members(BondId::new(i)))
            .collect::<Vec<_>>()
    );
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

#[test]
fn drawing_svg_boundary_preparation_product() {
    let mut calls = 0;
    let mut differences = Vec::new();
    for case in source::CASES {
        let start = differences.len();
        let input = match literal(&case.s0) {
            Ok(v) => v,
            Err(e) => {
                differences.push(format!("case{} S0 literal {e}", case.id));
                continue;
            }
        };
        let expected = match literal(&case.s1) {
            Ok(v) => v,
            Err(e) => {
                differences.push(format!("case{} S1 literal {e}", case.id));
                continue;
            }
        };
        observe("S1_INPUT", case.id, &input);
        println!("S1_NATIVE_EXPECTED_CASE {} {:#?}", case.id, case.s1);
        let baseline = (
            input.topology.clone(),
            input.coordinates.clone(),
            input.properties.clone(),
            input.valence.clone(),
            input.rings.clone(),
        );
        let baseline_bits = bits(&input.coordinates);
        let checkpoint = || {
            (
                &input.topology,
                &input.coordinates,
                &input.properties,
                &input.valence,
                &input.rings,
            ) == (
                &baseline.0,
                &baseline.1,
                &baseline.2,
                &baseline.3,
                &baseline.4,
            ) && bits(&input.coordinates) == baseline_bits
        };
        let before = checkpoint();
        let result = prepare(DrawingInput {
            topology: &input.topology,
            coordinates: &input.coordinates,
            properties: &input.properties,
            valence: Some(&input.valence),
            rings: Some(&input.rings),
        });
        let after = checkpoint();
        calls += 1;
        println!(
            "S1_CASE {} PRESERVATION before={before} after={after}",
            case.id
        );
        let label = format!("case{}", case.id);
        match result {
            Ok(output) => {
                observe("S1_OUTPUT", case.id, &output);
                // Every output row/bit is emitted before collecting equality.
                compare(
                    &mut differences,
                    &format!("{label} topology"),
                    &output.topology,
                    &expected.topology,
                );
                compare(
                    &mut differences,
                    &format!("{label} coordinates"),
                    &output.coordinates,
                    &expected.coordinates,
                );
                let actual_bits = bits(&output.coordinates);
                let expected_bits = bits(&expected.coordinates);
                println!(
                    "S1_CASE {} XY_MATCH_COUNT {}/{} FIRST_BIT_DIFF {:?}",
                    case.id,
                    actual_bits
                        .iter()
                        .zip(&expected_bits)
                        .filter(|(a, b)| a == b)
                        .count(),
                    expected_bits.len(),
                    actual_bits
                        .iter()
                        .zip(&expected_bits)
                        .position(|(a, b)| a != b)
                );
                compare(
                    &mut differences,
                    &format!("{label} exact_bits"),
                    actual_bits,
                    expected_bits,
                );
                compare(
                    &mut differences,
                    &format!("{label} properties"),
                    &output.properties,
                    &expected.properties,
                );
                compare(
                    &mut differences,
                    &format!("{label} valence"),
                    &output.valence,
                    &expected.valence,
                );
                compare(
                    &mut differences,
                    &format!("{label} rings"),
                    &output.rings,
                    &expected.rings,
                );
            }
            Err(error) => {
                println!("S1_CASE {} ERROR {error:?} {error}", case.id);
                differences.push(format!("{label} prepare {error:?}"));
            }
        }
        compare(
            &mut differences,
            &format!("{label} preservation"),
            (before, after),
            (true, true),
        );
        println!(
            "S1_CASE {} DIFFERENCES {:?}",
            case.id,
            &differences[start..]
        );
    }
    println!("S1_CENSUS calls={calls} differences={differences:#?}");
    assert_eq!(calls, 6);
    assert!(
        differences.is_empty(),
        "all six actual preparation observations collected"
    );
}

fn namespace_projection(svg: &str) -> String {
    svg.replace(
        "xmlns:rdkit='http://www.rdkit.org/xml'",
        "xmlns:tool='__tool_namespace__'",
    )
    .replace(
        "xmlns:cosmolkit='https://www.cosmol.org'",
        "xmlns:tool='__tool_namespace__'",
    )
    .replace("rdkit:", "tool:")
    .replace("cosmolkit:", "tool:")
}

fn sha256(bytes: &[u8]) -> Result<String, String> {
    // Existing developer checksum executable observes bytes only: no source
    // invocation, file generation, fixture refresh or dependency addition.
    let mut child = Command::new("sha256sum")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| e.to_string())?;
    child
        .stdin
        .take()
        .ok_or("checksum stdin")?
        .write_all(bytes)
        .map_err(|e| e.to_string())?;
    let output = child.wait_with_output().map_err(|e| e.to_string())?;
    if !output.status.success() {
        return Err(format!("checksum status {}", output.status));
    }
    String::from_utf8(output.stdout)
        .map_err(|e| e.to_string())?
        .split_whitespace()
        .next()
        .map(str::to_owned)
        .ok_or("checksum output".into())
}

#[test]
fn drawing_svg_boundary_renderer_product() {
    let mut calls = 0;
    let mut differences = Vec::new();
    for case in source::CASES {
        let start = differences.len();
        let input = match literal(&case.s1) {
            Ok(v) => v,
            Err(e) => {
                differences.push(format!("case{} S1 literal {e}", case.id));
                continue;
            }
        };
        observe("S2_INPUT", case.id, &input);
        let baseline = (
            input.topology.clone(),
            input.coordinates.clone(),
            input.properties.clone(),
            input.valence.clone(),
            input.rings.clone(),
        );
        let baseline_bits = bits(&input.coordinates);
        let checkpoint = || {
            (
                &input.topology,
                &input.coordinates,
                &input.properties,
                &input.valence,
                &input.rings,
            ) == (
                &baseline.0,
                &baseline.1,
                &baseline.2,
                &baseline.3,
                &baseline.4,
            ) && bits(&input.coordinates) == baseline_bits
        };
        let before = checkpoint();
        let result = render_prepared_svg(&input.borrow(), 300, 300);
        let after = checkpoint();
        calls += 1;
        println!(
            "S2_CASE {} PRESERVATION before={before} after={after}",
            case.id
        );
        let label = format!("case{}", case.id);
        match result {
            Ok(svg) => {
                println!(
                    "S2_CASE {} SVG_BEGIN\n{svg}S2_CASE {} SVG_END",
                    case.id, case.id
                );
                let actual = namespace_projection(&svg);
                let expected = namespace_projection(case.svg);
                let first = actual
                    .bytes()
                    .zip(expected.bytes())
                    .position(|(a, b)| a != b)
                    .or_else(|| {
                        (actual.len() != expected.len()).then_some(actual.len().min(expected.len()))
                    });
                println!(
                    "S2_CASE {} RAW_BYTES {} SOURCE_BYTES {} RAW_SHA {:?} SOURCE_SHA {} PROJECTED_BYTES {} SOURCE_PROJECTED_BYTES {} FIRST_BYTE_DIFF {first:?}",
                    case.id,
                    svg.len(),
                    case.svg.len(),
                    sha256(svg.as_bytes()),
                    case.svg_sha,
                    actual.len(),
                    expected.len()
                );
                match sha256(case.svg.as_bytes()) {
                    Ok(hash) => compare(
                        &mut differences,
                        &format!("{label} frozen source checksum"),
                        hash.as_str(),
                        case.svg_sha,
                    ),
                    Err(error) => differences.push(format!("{label} checksum {error}")),
                }
                if actual != expected {
                    differences.push(format!(
                        "{label} SVG first={first:?} actual_bytes={} source_bytes={}",
                        actual.len(),
                        expected.len()
                    ));
                }
            }
            Err(error) => {
                println!("S2_CASE {} ERROR {error:?} {error}", case.id);
                differences.push(format!("{label} render {error:?}"));
            }
        }
        compare(
            &mut differences,
            &format!("{label} preservation"),
            (before, after),
            (true, true),
        );
        println!(
            "S2_CASE {} DIFFERENCES {:?}",
            case.id,
            &differences[start..]
        );
    }
    println!("S2_CENSUS calls={calls} differences={differences:#?}");
    assert_eq!(calls, 6);
    assert!(
        differences.is_empty(),
        "all six independent same-prepared SVG observations collected"
    );
}

#[allow(dead_code)]
mod stage_source {
    include!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../testdata/depiction/expected/rdkit/svg_boundary_current/stages.rs"
    ));
}

#[test]
fn drawing_prepare_stage_boundaries() {
    use super::draw_prepare_stage_probe::{Guard, Observation, coordinate_identity};

    fn check_stage(
        case: usize,
        observation: &Observation,
        expected: &PreparedDrawing,
        differences: &mut Vec<String>,
    ) {
        let label = format!("case{case} {}", observation.stage);
        let actual_bits = bits(&observation.coordinates);
        let expected_bits = bits(&expected.coordinates);
        // Whole actual observation and native source inventories are emitted
        // before equality. No intermediate chemistry/cache update is invoked.
        println!("PREP_STAGE_OBSERVATION case={case} {observation:#?}");
        for atom in &observation.topology.atoms {
            println!(
                "PREP_STAGE_ATOM case={case} stage={} atom={} props={:?} computed={:?}",
                observation.stage,
                atom.id().index(),
                atom.props(),
                atom.computed_prop_names()
            );
        }
        let first_atom = observation
            .topology
            .atoms
            .iter()
            .zip(&expected.topology.atoms)
            .position(|(a, b)| a != b);
        let first_bond = observation
            .topology
            .bonds
            .iter()
            .zip(&expected.topology.bonds)
            .position(|(a, b)| a != b);
        let first_bit = actual_bits
            .iter()
            .zip(&expected_bits)
            .position(|(a, b)| a != b)
            .or_else(|| {
                (actual_bits.len() != expected_bits.len())
                    .then_some(actual_bits.len().min(expected_bits.len()))
            });
        if let Some(i) = first_atom {
            let a = &observation.topology.atoms[i];
            let e = &expected.topology.atoms[i];
            println!(
                "PREP_FIRST_ATOM case={case} stage={} row={i} props_equal={} computed_equal={} ACTUAL {a:#?} SOURCE {e:#?}",
                observation.stage,
                a.props() == e.props(),
                a.computed_prop_names() == e.computed_prop_names()
            );
        }
        if let Some(i) = first_bond {
            let a = &observation.topology.bonds[i];
            let e = &expected.topology.bonds[i];
            let field = if a.begin() != e.begin() || a.end() != e.end() {
                "endpoints"
            } else if a.order() != e.order() {
                "order"
            } else if a.direction() != e.direction() {
                "direction"
            } else {
                "other_bond_field"
            };
            println!(
                "PREP_FIRST_BOND case={case} stage={} row={i} field={field} ACTUAL {a:#?} SOURCE {e:#?}",
                observation.stage
            );
        }
        if let Some(i) = first_bit {
            println!(
                "PREP_FIRST_BIT case={case} stage={} index={i} atom={} axis={} actual={:?} source={:?}",
                observation.stage,
                i / 2,
                i % 2,
                actual_bits.get(i),
                expected_bits.get(i)
            );
        }
        println!(
            "PREP_STAGE_SUMMARY case={case} stage={} topology={} coordinates={} properties={} rings={} bits={} identities={} first_atom={first_atom:?} first_bond={first_bond:?} first_bit={first_bit:?}",
            observation.stage,
            observation.topology == expected.topology,
            observation.coordinates == expected.coordinates,
            observation.properties == expected.properties,
            observation.rings.as_ref() == Some(&expected.rings),
            actual_bits == expected_bits,
            observation.coordinate_identity == coordinate_identity(&expected.coordinates)
        );
        compare(
            differences,
            &format!("{label} topology"),
            &observation.topology,
            &expected.topology,
        );
        compare(
            differences,
            &format!("{label} coordinates"),
            &observation.coordinates,
            &expected.coordinates,
        );
        compare(
            differences,
            &format!("{label} properties"),
            &observation.properties,
            &expected.properties,
        );
        compare(
            differences,
            &format!("{label} rings"),
            observation.rings.as_ref(),
            Some(&expected.rings),
        );
        compare(
            differences,
            &format!("{label} bits"),
            actual_bits,
            expected_bits,
        );
        compare(
            differences,
            &format!("{label} identities"),
            &observation.coordinate_identity,
            &coordinate_identity(&expected.coordinates),
        );
    }

    let mut calls = 0;
    let mut checkpoints = 0;
    let mut differences = Vec::new();
    for (case, stages) in source::CASES.iter().zip(stage_source::CASES) {
        let start = differences.len();
        compare(&mut differences, "source case order", case.id, stages.id);
        println!(
            "PREP_SOURCE_CASE {} candidates={:?} CK_candidates=omitted_existing_list_moves_into_AddHs; INPUT_VALENCE=unchanged_input_not_updated_stage_cache",
            case.id, stages.candidates
        );
        let input = match literal(&case.s0) {
            Ok(v) => v,
            Err(e) => {
                differences.push(format!("case{} S0 literal {e}", case.id));
                continue;
            }
        };
        let baseline = (
            input.topology.clone(),
            input.coordinates.clone(),
            input.properties.clone(),
            input.valence.clone(),
            input.rings.clone(),
        );
        let baseline_bits = bits(&input.coordinates);
        let baseline_identity = coordinate_identity(&input.coordinates);
        let checkpoint = || {
            (
                &input.topology,
                &input.coordinates,
                &input.properties,
                &input.valence,
                &input.rings,
            ) == (
                &baseline.0,
                &baseline.1,
                &baseline.2,
                &baseline.3,
                &baseline.4,
            ) && bits(&input.coordinates) == baseline_bits
                && coordinate_identity(&input.coordinates) == baseline_identity
        };
        let guard = match Guard::arm() {
            Ok(g) => g,
            Err(e) => {
                differences.push(format!("case{} collector {e}", case.id));
                continue;
            }
        };
        let before = checkpoint();
        let result = prepare(DrawingInput {
            topology: &input.topology,
            coordinates: &input.coordinates,
            properties: &input.properties,
            valence: Some(&input.valence),
            rings: Some(&input.rings),
        });
        let after = checkpoint();
        calls += 1;
        let observations = guard.finish();
        checkpoints += observations.len();
        println!(
            "PREP_PRESERVATION case={} before={before} after={after} actual_checkpoints={}",
            case.id,
            observations.len()
        );
        compare(
            &mut differences,
            &format!("case{} input preservation", case.id),
            (before, after),
            (true, true),
        );
        compare(
            &mut differences,
            &format!("case{} stage count", case.id),
            observations.len(),
            5,
        );
        for (index, actual) in observations.iter().enumerate() {
            let Some(native) = stages.stages.get(index) else {
                println!("PREP_EXTRA_STAGE case={} {actual:#?}", case.id);
                differences.push(format!("case{} extra stage{}", case.id, index));
                continue;
            };
            println!(
                "PREP_NATIVE_STAGE case={} name={} {native:#?}",
                case.id,
                stage_source::NAMES[index]
            );
            compare(
                &mut differences,
                &format!("case{} stage order{index}", case.id),
                actual.stage,
                stage_source::NAMES[index],
            );
            match literal(native) {
                Ok(expected) => check_stage(case.id, actual, &expected, &mut differences),
                Err(e) => differences.push(format!("case{} stage{index} literal {e}", case.id)),
            }
        }
        match result {
            Ok(output) => {
                observe("PREP_FINAL", case.id, &output);
                println!("PREP_FINAL_NATIVE case={} {:#?}", case.id, case.s1);
                match literal(&case.s1) {
                    Ok(expected) => {
                        let final_observation = Observation {
                            stage: "returned_final",
                            coordinate_identity: coordinate_identity(&output.coordinates),
                            topology: output.topology,
                            coordinates: output.coordinates,
                            properties: output.properties,
                            rings: Some(output.rings),
                        };
                        check_stage(case.id, &final_observation, &expected, &mut differences);
                        compare(
                            &mut differences,
                            &format!("case{} returned final valence", case.id),
                            &output.valence,
                            &expected.valence,
                        );
                    }
                    Err(e) => differences.push(format!("case{} final source literal {e}", case.id)),
                }
            }
            Err(error) => {
                println!("PREP_ERROR case={} {error:?} {error}", case.id);
                differences.push(format!("case{} prepare {error:?}", case.id));
            }
        }
        println!(
            "PREP_CASE_DIFFERENCES case={} {:#?}",
            case.id,
            &differences[start..]
        );
    }
    println!("PREP_CENSUS calls={calls} checkpoints={checkpoints} differences={differences:#?}");
    assert_eq!(calls, 6, "six actual prepare Result returns");
    assert_eq!(checkpoints, 30, "thirty actual-site snapshots");
    assert!(
        differences.is_empty(),
        "all six cases observed before diagnostic equality"
    );
}
