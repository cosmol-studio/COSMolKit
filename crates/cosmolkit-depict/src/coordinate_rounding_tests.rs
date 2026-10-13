use super::*;
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, ChiralTag, Element, Hybridization,
};

fn molecule() -> TopologyBlock {
    // Exact sanitized graph from pinned RDKit for CC(C)[C@@]12C=C3CC[C@H]4[C@](C)(C(=O)O)CCC[C@]4(C)[C@]3(CC1)OO2.
    let rows = [
        (6, 0, false, 0, 4),
        (6, 0, false, 0, 4),
        (6, 0, false, 0, 4),
        (6, 0, true, 1, 4),
        (6, 0, false, 0, 3),
        (6, 0, false, 0, 3),
        (6, 0, false, 0, 4),
        (6, 0, false, 0, 4),
        (6, 1, true, 1, 4),
        (6, 0, true, 2, 4),
        (6, 0, false, 0, 4),
        (6, 0, false, 0, 3),
        (8, 0, false, 0, 3),
        (8, 0, false, 0, 3),
        (6, 0, false, 0, 4),
        (6, 0, false, 0, 4),
        (6, 0, false, 0, 4),
        (6, 0, true, 2, 4),
        (6, 0, false, 0, 4),
        (6, 0, true, 2, 4),
        (6, 0, false, 0, 4),
        (6, 0, false, 0, 4),
        (8, 0, false, 0, 4),
        (8, 0, false, 0, 4),
    ];
    let ranks = [
        0, 11, 0, 16, 12, 14, 6, 4, 10, 15, 2, 18, 20, 19, 7, 3, 5, 13, 1, 17, 9, 8, 22, 21,
    ];
    let atoms = rows
        .into_iter()
        .enumerate()
        .map(|(i, (z, h, no_implicit, chiral, hybrid))| {
            Atom::from_spec(
                AtomId::new(i),
                AtomSpec::new(Element::from_atomic_number(z).unwrap())
                    .with_explicit_hydrogens(h)
                    .with_no_implicit(no_implicit)
                    .with_chiral_tag(ChiralTag::from_rdkit_code(chiral).unwrap())
                    .with_hybridization(Hybridization::from_rdkit_code(hybrid).unwrap())
                    .with_computed_prop("_CIPRank", cosmolkit_model::PropertyValue::UInt(ranks[i]))
                    .unwrap(),
            )
        })
        .collect();
    let edges = [
        (0, 1, 1),
        (1, 2, 1),
        (1, 3, 1),
        (3, 4, 1),
        (4, 5, 2),
        (5, 6, 1),
        (6, 7, 1),
        (7, 8, 1),
        (8, 9, 1),
        (9, 10, 1),
        (9, 11, 1),
        (11, 12, 2),
        (11, 13, 1),
        (9, 14, 1),
        (14, 15, 1),
        (15, 16, 1),
        (16, 17, 1),
        (17, 18, 1),
        (17, 19, 1),
        (19, 20, 1),
        (20, 21, 1),
        (19, 22, 1),
        (22, 23, 1),
        (21, 3, 1),
        (23, 3, 1),
        (19, 5, 1),
        (17, 8, 1),
    ];
    let bonds = edges
        .into_iter()
        .enumerate()
        .map(|(i, (a, b, order))| {
            Bond::from_spec(
                BondId::new(i),
                BondSpec::new(
                    AtomId::new(a),
                    AtomId::new(b),
                    BondOrder::from_rdkit_code(order).unwrap(),
                )
                .with_conjugated(i == 11 || i == 12),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

#[test]
fn bridge_collision_shortening_matches_source_bits() {
    let topology = molecule();
    let original = topology.clone();
    // Freeze the source regression's ordered ring seeds at the embedding
    // boundary. Candidate RDL enumeration is tested independently; changing
    // its order must not discard this collision-shortening numerical coverage.
    let mut rings = cosmolkit_core::RingInfo::new(
        cosmolkit_core::RingFindType::SymmSssr,
        24,
        topology.bonds.len(),
    );
    for atoms in [
        vec![4, 5, 19, 20, 21, 3],
        vec![4, 5, 19, 22, 23, 3],
        vec![6, 7, 8, 17, 19, 5],
        vec![9, 14, 15, 16, 17, 8],
        vec![20, 21, 3, 23, 22, 19],
    ] {
        let bonds = atoms
            .iter()
            .zip(atoms.iter().cycle().skip(1))
            .take(atoms.len())
            .map(|(&a, &b)| {
                topology
                    .bonds
                    .iter()
                    .find(|bond| {
                        let ends = [bond.begin().index(), bond.end().index()];
                        ends == [a, b] || ends == [b, a]
                    })
                    .unwrap()
                    .id()
                    .index()
            })
            .collect::<Vec<_>>();
        rings.add_ring(&atoms, &bonds).unwrap();
    }
    // Pinned RDKit 2026.03.6 Compute2DCoords, before and after canonical orientation.
    // The previous hypot norm first diverged inside bridge-ring shortening,
    // affecting atoms 22 and 23 before orientation propagated the error.
    let expected = [
        [
            [13828302655841107970, 4615957940093015616],
            [13832806255468478464, 4613032766328682880],
            [13837309855095848960, 4613032766328682878],
            [13828302655841107966, 4608529166701312385],
            [13832806255345191609, 4497254011592217493],
            [13828302655841107974, 13831901203413728794],
            [13832806255642832412, 13836404803032463744],
            [13828302656538523756, 13839329976847128126],
            [4604930618288916367, 13839329976947791423],
            [4609434218090640818, 13840908402835995008],
            [4610607280694724360, 13842571595280796514],
            [4591182818165802240, 13841486023459961737],
            [13830818392256967452, 13840400451656606957],
            [13818665362398948680, 13843149215887283246],
            [4613937817979542231, 13840908402886326656],
            [4615626667926983138, 13839329977149118016],
            [4613937818153896174, 13836404803334453632],
            [4609434218439348704, 13836404803233790337],
            [4612248968380809212, 13831901203757414790],
            [4604930618986332152, 13831901203556088195],
            [4609434218613702657, 0],
            [4604930618986332163, 4608529166701312384],
            [4598992703716521539, 13827064389467794724],
            [13823939069658822261, 4603237882951781131],
        ],
        [
            [4618330813668136569, 13827468725801305994],
            [4617183095829657892, 4601651626131725948],
            [4617682185283041166, 4611127500071139263],
            [4614882959020833481, 4592617356063157832],
            [4612587523427660084, 4608131613708342792],
            [4606621179657817933, 4606084922712388234],
            [13813922715015828520, 4611589407129677340],
            [13833005976333130690, 4610091464235036620],
            [13835002334338998424, 4596839746505323172],
            [13838323809416626063, 13816660091402616800],
            [13840618306494466389, 4589566533895056512],
            [13838489472130455032, 4608887111507692856],
            [13835650650111201210, 4612116850449385690],
            [13840527805886536016, 4611829914212114171],
            [13839321988419559929, 13833043512430756213],
            [13837026552816448800, 13836528652051910222],
            [13832407821138643901, 13835779680604589861],
            [13830268470610759128, 13829540797628156328],
            [4587529272263596032, 13835003363032375292],
            [4602578108344852054, 13826544911838874894],
            [4609496312633298152, 13833505420000915954],
            [4613884780114066931, 13832007476909960504],
            [4608464189178695241, 13822869770738273858],
            [4613053197992290535, 13805562683127943680],
        ],
    ];
    for (canonical_orientation, expected) in [false, true].into_iter().zip(expected) {
        let mut fragments = compute_initial_coordinates(&topology, &rings, None, false).unwrap();
        for fragment in &mut fragments {
            fragment.remove_collisions_bond_and_spiro_flip().unwrap();
            fragment.remove_collisions_open_angles().unwrap();
            fragment.remove_collisions_shorten_bonds().unwrap();
            fragment.remove_collisions_path_angle_expansion().unwrap();
        }
        orient_and_shift_fragments(&mut fragments, canonical_orientation, Some(0));
        let actual = (0..24)
            .map(|atom| {
                let loc = fragments
                    .iter()
                    .find_map(|fragment| fragment.atoms.get(&atom))
                    .unwrap()
                    .loc;
                loc.map(f64::to_bits)
            })
            .collect::<Vec<_>>();
        assert_eq!(
            actual, expected,
            "canonical_orientation={canonical_orientation}"
        );
        assert_eq!(topology, original);
    }
}
