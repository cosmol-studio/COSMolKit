//! Two independent line22 observations with pre-CK source literals.
use super::*;
use crate::draw::render_prepared_svg;
use cosmolkit_core::RingFindType;
use cosmolkit_model::{
    Atom, AtomSpec, Bond, BondSpec, BondStereo, Conformer2D, Element, Hybridization, PropertyValue,
};

const XY_BITS: [[u64; 2]; 8] = [
    [13828302655841107970, 13832806255468478466],
    [13828302655841107969, 13596367275031527424],
    [13835621005235585024, 0],
    [13828302655841107968, 4609434218613702656],
    [4604930618986332159, 13594115475217842176],
    [4604930618986332164, 13832806255468478466],
    [4604930618986332160, 4609434218613702657],
    [4612248968380809216, 0],
];

fn literal_topology(prepared: bool) -> TopologyBlock {
    let elements = [
        Element::O,
        Element::C,
        Element::C,
        Element::CL,
        Element::C,
        Element::O,
        Element::CL,
        Element::C,
    ];
    let tags = [
        ChiralTag::Unspecified,
        ChiralTag::TetrahedralCcw,
        ChiralTag::Unspecified,
        ChiralTag::Unspecified,
        ChiralTag::TetrahedralCw,
        ChiralTag::Unspecified,
        ChiralTag::Unspecified,
        ChiralTag::Unspecified,
    ];
    let atoms = (0..8)
        .map(|i| {
            let mut spec = AtomSpec::new(elements[i])
                .with_no_implicit([false, true, false, false, true, false, false, false][i])
                .with_chiral_tag(tags[i])
                .with_hybridization(Hybridization::Sp3)
                .with_computed_prop("_CIPRank", PropertyValue::Int([2, 1, 0, 3, 1, 2, 3, 0][i]))
                .unwrap();
            if i == 1 || i == 4 {
                spec = spec
                    .with_computed_prop("_ChiralityPossible", PropertyValue::Int(1))
                    .unwrap()
                    .with_computed_prop("_CIPCode", if i == 1 { "R" } else { "S" })
                    .unwrap();
            }
            Atom::from_spec(AtomId::new(i), spec)
        })
        .collect();
    let bonds = [(0, 1), (1, 2), (1, 3), (1, 4), (4, 5), (4, 6), (4, 7)]
        .into_iter()
        .enumerate()
        .map(|(i, (a, b))| {
            Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single).with_direction(
                    if prepared && (i == 1 || i == 6) {
                        BondDirection::BeginWedge
                    } else {
                        BondDirection::None
                    },
                ),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn literal_valence() -> ValenceAssignment {
    ValenceAssignment {
        explicit_valence: vec![1, 4, 1, 1, 4, 1, 1, 1],
        implicit_hydrogens: vec![1, 0, 3, 0, 0, 1, 0, 3],
    }
}

fn bits(coordinates: &CoordinateBlock) -> Vec<Vec<[u64; 2]>> {
    coordinates
        .conformers_2d
        .iter()
        .map(|c| {
            c.coordinates()
                .iter()
                .map(|p| p.map(f64::to_bits))
                .collect()
        })
        .collect()
}

fn observe(
    stage: &str,
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
) {
    println!(
        "{stage}_TOPOLOGY {topology:#?}\n{stage}_COORDINATES {coordinates:#?}\n{stage}_PROPERTIES {properties:#?}\n{stage}_VALENCE {valence:#?}\n{stage}_RINGS {rings:#?}"
    );
    for a in &topology.atoms {
        println!(
            "{stage}_ATOM {} z={} isotope={:?} charge={} explicit_h={} no_implicit={} radicals={} aromatic={} hybrid={} chiral={} permutation={:?} map={:?} unknown={} props={:?} computed={:?} query=absent_in_concrete_Atom",
            a.id().index(),
            a.atomic_number(),
            a.isotope(),
            a.formal_charge(),
            a.explicit_hydrogens(),
            a.no_implicit(),
            a.radical_electrons(),
            a.is_aromatic(),
            a.hybridization().rdkit_code(),
            a.chiral_tag().rdkit_code(),
            a.chiral_permutation(),
            a.atom_map(),
            a.unknown_stereo(),
            a.props(),
            a.computed_prop_names()
        );
    }
    for b in &topology.bonds {
        println!(
            "{stage}_BOND {} endpoints=({}, {}) order={:?} aromatic={} conjugated={} direction={} stereo={:?} stereo_atoms={:?} unknown={} props={:?} computed={:?} query=absent_in_concrete_Bond",
            b.id().index(),
            b.begin().index(),
            b.end().index(),
            b.order(),
            b.is_aromatic(),
            b.is_conjugated(),
            b.direction().rdkit_code(),
            b.stereo(),
            b.stereo_atoms(),
            b.unknown_stereo(),
            b.props(),
            b.computed_prop_names()
        );
    }
    println!(
        "{stage}_ADJACENCY {:?}",
        (0..topology.atoms.len())
            .map(|i| topology
                .adjacency
                .neighbors_of(i)
                .iter()
                .map(|n| (n.bond.index(), n.atom_index))
                .collect::<Vec<_>>())
            .collect::<Vec<_>>()
    );
    println!(
        "{stage}_RING_MEMBERS atoms={:?} bonds={:?}",
        (0..8)
            .map(|i| rings.atom_members(AtomId::new(i)))
            .collect::<Vec<_>>(),
        (0..7)
            .map(|i| rings.bond_members(BondId::new(i)))
            .collect::<Vec<_>>()
    );
    println!("{stage}_XY_BITS {:?}", bits(coordinates));
}

fn assert_projection(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    prepared: bool,
) {
    // The entire literal topology checks all atom/bond fields and property types,
    // ordered adjacency and empty groups, without expected CK reconstruction.
    assert_eq!(topology, &literal_topology(prepared));
    assert_eq!(valence, &literal_valence());
    assert!(rings.is_initialized());
    assert_eq!((rings.atom_row_count(), rings.bond_row_count()), (8, 7));
    assert!(rings.atom_rings().is_empty() && rings.bond_rings().is_empty());
    for i in 0..8 {
        assert_eq!(rings.num_atom_rings(AtomId::new(i)), 0);
        assert!(rings.atom_members(AtomId::new(i)).is_empty());
    }
    for i in 0..7 {
        assert_eq!(rings.num_bond_rings(BondId::new(i)), 0);
        assert!(rings.bond_members(BondId::new(i)).is_empty());
    }
    for b in &topology.bonds {
        assert_eq!(b.stereo(), BondStereo::None);
        assert_eq!(b.stereo_atoms(), None);
    }
}

#[test]
fn drawing_line22_probe_preparation() {
    let topology = literal_topology(false);
    let coordinates = CoordinateBlock::default();
    let properties = MoleculeProperties::default();
    let valence = literal_valence();
    // Fixed modeled input quality; native source quality is unexposed.
    let rings = RingInfo::new(RingFindType::SymmSssr, 8, 7);
    topology.validate().unwrap();
    coordinates.validate_for_atom_count(8).unwrap();
    check_properties(&properties, 8, 7).unwrap();
    assert_projection(&topology, &valence, &rings, false);
    let baseline = (
        topology.clone(),
        coordinates.clone(),
        properties.clone(),
        valence.clone(),
        rings.clone(),
    );
    let baseline_bits = bits(&coordinates);
    let checkpoint = || {
        assert_eq!(
            (&topology, &coordinates, &properties, &valence, &rings),
            (
                &baseline.0,
                &baseline.1,
                &baseline.2,
                &baseline.3,
                &baseline.4
            )
        );
        assert_eq!(bits(&coordinates), baseline_bits);
    };
    checkpoint();
    let result = prepare(DrawingInput {
        topology: &topology,
        coordinates: &coordinates,
        properties: &properties,
        valence: Some(&valence),
        rings: Some(&rings),
    });
    checkpoint();
    if let Err(error) = &result {
        println!("S1_PREPARATION_ERROR {error:?} {error}");
    }
    let output = result.expect("S1 actual prepare result");
    observe(
        "S1",
        &output.topology,
        &output.coordinates,
        &output.properties,
        &output.valence,
        &output.rings,
    );
    // Complete observations, including all sixteen bits, precede equality.
    assert_projection(&output.topology, &output.valence, &output.rings, true);
    assert_eq!(output.properties, properties);
    assert_eq!(output.coordinates.conformers_2d.len(), 1);
    assert!(output.coordinates.conformers_3d.is_empty());
    assert_eq!(output.coordinates.conformers_2d[0].id(), 0);
    assert_eq!(bits(&output.coordinates), vec![XY_BITS.to_vec()]);
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

#[test]
fn drawing_line22_probe_same_prepared_renderer() {
    let topology = literal_topology(true);
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(
            0,
            XY_BITS.map(|p| p.map(f64::from_bits)).to_vec(),
        )],
        conformers_3d: vec![],
        source_coordinate_dim: None,
        source_conformer_order: Some(vec![cosmolkit_model::CoordinateDimension::TwoD]),
    };
    let properties = MoleculeProperties::default();
    let valence = literal_valence();
    let rings = RingInfo::new(RingFindType::SymmSssr, 8, 7);
    topology.validate().unwrap();
    coordinates.validate_for_atom_count(8).unwrap();
    check_properties(&properties, 8, 7).unwrap();
    assert_projection(&topology, &valence, &rings, true);
    assert_eq!(bits(&coordinates), vec![XY_BITS.to_vec()]);
    observe(
        "S2_INPUT",
        &topology,
        &coordinates,
        &properties,
        &valence,
        &rings,
    );
    let baseline = (
        topology.clone(),
        coordinates.clone(),
        properties.clone(),
        valence.clone(),
        rings.clone(),
    );
    let baseline_bits = bits(&coordinates);
    let checkpoint = || {
        assert_eq!(
            (&topology, &coordinates, &properties, &valence, &rings),
            (
                &baseline.0,
                &baseline.1,
                &baseline.2,
                &baseline.3,
                &baseline.4
            )
        );
        assert_eq!(bits(&coordinates), baseline_bits);
    };
    checkpoint();
    let result = render_prepared_svg(
        &PreparedDrawingInput {
            topology: &topology,
            layout: &coordinates.conformers_2d[0],
            properties: &properties,
            valence: &valence,
            rings: &rings,
        },
        300,
        300,
    );
    checkpoint();
    if let Err(error) = &result {
        println!("S2_RENDER_ERROR {error:?} {error}");
    }
    let svg = result.expect("S2 actual independently same-prepared render result");
    let svg = std::str::from_utf8(&svg).expect("fixed source SVG UTF8");
    println!("S2_SAME_PREPARED_SVG_BEGIN\n{svg}S2_SAME_PREPARED_SVG_END");
    let expected = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../testdata/depiction/expected/rdkit/drawing_state/line22.svg"
    ));
    assert_eq!(namespace_projection(svg), namespace_projection(expected));
}
