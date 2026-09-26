use cosmolkit_core::{RingInfo, find_sssr_from_parts};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CoordinateBlock, TopologyBlock,
};
use cosmolkit_search::{McsCandidateMatchError, McsError, McsParameters, SearchTarget, find_mcs};
use cosmolkit_types::{BondOrder, BondStereo, ChiralTag, Element};

fn topology(elements: &[Element], edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
    let atoms = elements
        .iter()
        .copied()
        .enumerate()
        .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
        .collect();
    let bonds = edges
        .iter()
        .copied()
        .enumerate()
        .map(|(index, (begin, end, order))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .expect("fixed MCS input is structurally valid")
}

fn path(atoms: usize) -> TopologyBlock {
    topology(
        &vec![Element::C; atoms],
        &(0..atoms.saturating_sub(1))
            .map(|index| (index, index + 1, BondOrder::Single))
            .collect::<Vec<_>>(),
    )
}

fn target<'a>(
    topology: &'a TopologyBlock,
    coordinates: &'a CoordinateBlock,
    rings: Option<&'a RingInfo>,
) -> SearchTarget<'a> {
    SearchTarget::new(topology, coordinates, &topology.stereo_groups, rings, None)
}

#[test]
fn q134_entry_rejects_input_threshold_and_required_valence_with_typed_errors() {
    let coordinates = CoordinateBlock::default();
    let single = path(1);
    let view = target(&single, &coordinates, None);
    assert!(matches!(
        find_mcs(&[], &McsParameters::default()),
        Err(McsCandidateMatchError::State(McsError::TooFewInputs {
            count: 0
        }))
    ));
    assert!(matches!(
        find_mcs(&[view], &McsParameters::default()),
        Err(McsCandidateMatchError::State(McsError::TooFewInputs {
            count: 1
        }))
    ));
    let mut too_high = McsParameters::default();
    too_high.threshold = 1.01;
    assert!(matches!(
        find_mcs(&[view, view], &too_high),
        Err(McsCandidateMatchError::State(McsError::ThresholdAboveOne))
    ));
    let mut valence = McsParameters::default();
    valence.atom_compare_parameters.match_valences = true;
    assert!(matches!(
        find_mcs(&[view, view], &valence),
        Err(McsCandidateMatchError::State(
            McsError::MissingValence { .. }
        ))
    ));
}

#[test]
fn q134_entry_empty_single_atom_and_bond_results_follow_source_branches() {
    let coordinates = CoordinateBlock::default();
    let empty = path(0);
    let empty_view = target(&empty, &coordinates, None);
    let result = find_mcs(&[empty_view, empty_view], &McsParameters::default()).unwrap();
    assert_eq!(
        (result.atom_count, result.bond_count, result.completed),
        (0, 0, true)
    );
    assert!(result.query.is_none());
    assert!(result.smarts.is_empty());

    let atom = path(1);
    let atom_view = target(&atom, &coordinates, None);
    let result = find_mcs(&[atom_view, atom_view], &McsParameters::default()).unwrap();
    assert_eq!(
        (result.atom_count, result.bond_count, result.completed),
        (1, 0, true)
    );
    assert_eq!(result.query.as_ref().unwrap().num_atoms(), 1);
    assert!(!result.smarts.is_empty());

    let bond = path(2);
    let bond_view = target(&bond, &coordinates, None);
    let result = find_mcs(&[bond_view, bond_view], &McsParameters::default()).unwrap();
    assert_eq!(
        (result.atom_count, result.bond_count, result.completed),
        (2, 1, true)
    );
    assert_eq!(result.query.as_ref().unwrap().num_bonds(), 1);
    assert!(!result.smarts.is_empty());
}

#[test]
fn q134_entry_threshold_multi_candidate_and_store_all_keep_result_context() {
    let coordinates = CoordinateBlock::default();
    let small = path(2);
    let medium = topology(
        &[Element::C, Element::C, Element::O],
        &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)],
    );
    let large = topology(
        &[Element::C, Element::C, Element::O, Element::C],
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
        ],
    );
    let views = [
        target(&small, &coordinates, None),
        target(&medium, &coordinates, None),
        target(&large, &coordinates, None),
    ];
    let mut params = McsParameters::default();
    params.threshold = 0.5;
    let result = find_mcs(&views, &params).unwrap();
    assert!(result.completed);
    assert_eq!((result.atom_count, result.bond_count), (3, 2));
    assert_eq!(result.query.as_ref().unwrap().num_bonds(), 2);
    // Result carriers are source dummy query atoms; the atom-number predicates
    // and emitted SMARTS retain the oxygen identity from the second query.
    assert!(result.smarts.contains("[#8]"), "{}", result.smarts);

    params.store_all = true;
    let multi_all = find_mcs(&views, &params).unwrap();
    assert_eq!((multi_all.atom_count, multi_all.bond_count), (3, 2));
    assert!(multi_all.query.is_none());
    assert!(
        multi_all
            .degenerate
            .keys()
            .any(|smarts| smarts.contains("[#8]"))
    );

    let equal = path(3);
    let equal_views = [
        target(&equal, &coordinates, None),
        target(&equal, &coordinates, None),
    ];
    let stored = find_mcs(&equal_views, &params).unwrap();
    assert!(stored.completed);
    assert!(stored.bond_count > 0);
    assert!(stored.query.is_none());
    assert!(!stored.degenerate.is_empty());
    assert!(
        stored
            .degenerate
            .values()
            .all(|graph| graph.num_bonds() == stored.bond_count)
    );
}

#[test]
fn q134_entry_fused_ring_and_tetrahedral_options_reach_final_mapping() {
    let coordinates = CoordinateBlock::default();
    let fused = topology(
        &[Element::C; 4],
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (3, 1, BondOrder::Single),
        ],
    );
    let rings = find_sssr_from_parts(fused.atoms.len(), &fused.bonds, &fused.adjacency).unwrap();
    let fused_view = target(&fused, &coordinates, Some(&rings));
    let mut params = McsParameters::default();
    params.bond_compare_parameters.match_fused_rings_strict = true;
    let result = find_mcs(&[fused_view, fused_view], &params).unwrap();
    assert_eq!((result.atom_count, result.bond_count), (4, 5));

    let mut cw = topology(
        &[Element::C, Element::F, Element::CL, Element::BR],
        &[
            (0, 1, BondOrder::Single),
            (0, 2, BondOrder::Single),
            (0, 3, BondOrder::Single),
        ],
    );
    cw.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);
    let mut ccw = cw.clone();
    ccw.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCcw);
    let cw_view = target(&cw, &coordinates, None);
    let ccw_view = target(&ccw, &coordinates, None);
    params.bond_compare_parameters.match_fused_rings_strict = false;
    params.atom_compare_parameters.match_chiral_tag = true;
    let same = find_mcs(&[cw_view, cw_view], &params).unwrap();
    let opposite = find_mcs(&[cw_view, ccw_view], &params).unwrap();
    assert_eq!(same.bond_count, 3);
    assert!(opposite.bond_count < 3);
}

#[test]
fn q134_entry_double_bond_stereo_and_timeout_preserve_partial_result() {
    let coordinates = CoordinateBlock::default();
    let mut cis = topology(
        &[Element::F, Element::C, Element::C, Element::CL],
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
        ],
    );
    cis.bonds[1].set_stereo_atoms(Some([AtomId::new(0), AtomId::new(3)]));
    cis.bonds[1].set_stereo(BondStereo::Cis).unwrap();
    let mut trans = cis.clone();
    trans.bonds[1].set_stereo(BondStereo::Trans).unwrap();
    let cis_view = target(&cis, &coordinates, None);
    let trans_view = target(&trans, &coordinates, None);
    let mut plain = cis.clone();
    plain.bonds[1].set_stereo(BondStereo::None).unwrap();
    let plain_view = target(&plain, &coordinates, None);
    assert_eq!(
        find_mcs(&[plain_view, plain_view], &McsParameters::default())
            .unwrap()
            .bond_count,
        3
    );
    let mut params = McsParameters::default();
    params.atom_compare_parameters.match_chiral_tag = true;
    let mut bond_only = params.clone();
    bond_only.atom_compare_parameters.match_chiral_tag = false;
    bond_only.bond_compare_parameters.match_stereo = true;
    let source_bond_only = find_mcs(&[cis_view, cis_view], &bond_only).unwrap();
    assert_eq!(source_bond_only.bond_count, 3);
    let same = find_mcs(&[cis_view, cis_view], &params).unwrap();
    let opposite = find_mcs(&[cis_view, trans_view], &params).unwrap();
    // The source checks each partial seed. With MatchChiralTag enabled, the
    // two-bond intermediate has exactly one mapped stereo neighbor and is
    // rejected even for equal labels, before the three-bond seed can grow.
    assert_eq!(same.bond_count, 1);
    assert!(opposite.bond_count <= same.bond_count);

    let long = path(1000);
    let long_view = target(&long, &coordinates, None);
    params.atom_compare_parameters.match_chiral_tag = false;
    params.timeout = 1;
    let timed = find_mcs(&[long_view, long_view], &params).unwrap();
    assert!(!timed.completed);
    assert!(timed.atom_count >= 2);
    assert!(timed.bond_count >= 1);
    assert!(timed.query.is_some());

    params.store_all = true;
    let timed_all = find_mcs(&[long_view, long_view], &params).unwrap();
    assert!(!timed_all.completed);
    assert!(timed_all.bond_count >= 1);
    assert!(timed_all.query.is_none());
    assert!(!timed_all.degenerate.is_empty());
    assert!(
        timed_all
            .degenerate
            .values()
            .all(|graph| graph.num_bonds() == timed_all.bond_count)
    );
}
