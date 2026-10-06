// Fixed regressions for the private atropisomer cleanup phase.
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
use cosmolkit_types::{BondOrder, Element};

fn topology_from(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn cc_graph() -> TopologyBlock {
    topology_from(
        vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ],
        vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
        )],
    )
}

mod l04 {
    use super::{cc_graph, topology_from};
    use crate::{AtropisomerError, RingFindType, RingInfo, RingSearchParams, find_sssr};
    use crate::{
        atropisomer::cleanup_invalid_atropisomers, hybridization::HybridizationAssignment,
    };
    use cosmolkit_model::{AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec};
    use cosmolkit_types::{BondOrder, BondStereo, Element, Hybridization};

    fn hybs(values: &[Hybridization]) -> HybridizationAssignment {
        HybridizationAssignment {
            values: values.to_vec(),
        }
    }

    fn tagged_cc(stereo: BondStereo) -> cosmolkit_model::TopologyBlock {
        topology_from(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                    .with_stereo(stereo),
            )],
        )
    }

    #[test]
    fn sanitize_ring_l04_valid_sssr_and_symm_retain_acyclic_sp2_tag() {
        for stereo in [BondStereo::AtropCw, BondStereo::AtropCcw] {
            let input = tagged_cc(stereo);
            let hybridization = hybs(&[Hybridization::Sp2, Hybridization::Sp2]);
            let sssr = find_sssr(&input, &RingSearchParams::default()).unwrap();
            let symm = RingInfo::new(RingFindType::SymmSssr, 2, 1);
            for (name, rings) in [("sssr", sssr), ("symm", symm)] {
                let topology_snapshot = input.clone();
                let hybs_snapshot = hybridization.clone();
                let rings_snapshot = rings.clone();
                let output = cleanup_invalid_atropisomers(&input, &hybridization, &rings);
                assert!(output.is_ok(), "{name}/{stereo:?}: expected success");
                let output = output.unwrap();
                // Acyclic tag in no ring with Sp2 endpoints is retained.
                assert_eq!(output.bonds[0].stereo(), stereo, "{name}/{stereo:?}");
                assert_eq!(output, input, "{name}/{stereo:?}: topology changed");
                assert_eq!(input, topology_snapshot, "{name}: input mutated");
                assert_eq!(
                    hybridization, hybs_snapshot,
                    "{name}: hybridization mutated"
                );
                assert_eq!(rings, rings_snapshot, "{name}: rings mutated");
            }
        }
    }

    #[test]
    fn sanitize_ring_l04_no_tag_fast_rejected_before_hybridization_length() {
        // No tagged bonds: the strict rings check still runs eagerly and
        // must reject initialized Fast BEFORE the wrong hybridization
        // length is even examined.
        let input = cc_graph();
        let hybridization = hybs(&[Hybridization::Sp2, Hybridization::Sp2, Hybridization::Sp2]);
        let fast = RingInfo::new(RingFindType::Fast, 2, 1);
        let topology_snapshot = input.clone();
        let hybs_snapshot = hybridization.clone();
        let rings_snapshot = fast.clone();
        assert_eq!(
            cleanup_invalid_atropisomers(&input, &hybridization, &fast),
            Err(AtropisomerError::RingInfoNotSssr)
        );
        assert_eq!(input, topology_snapshot, "input mutated");
        assert_eq!(hybridization, hybs_snapshot, "hybridization mutated");
        assert_eq!(fast, rings_snapshot, "rings mutated");
    }

    #[test]
    fn sanitize_ring_l04_malformed_dimensions_rejected() {
        let input = cc_graph();
        let hybridization = hybs(&[Hybridization::Sp2, Hybridization::Sp2]);
        let atom_rows_wrong = RingInfo::new(RingFindType::Sssr, 3, 1);
        let bond_rows_wrong = RingInfo::new(RingFindType::Sssr, 2, 2);
        let atom_snapshot = atom_rows_wrong.clone();
        assert_eq!(
            cleanup_invalid_atropisomers(&input, &hybridization, &atom_rows_wrong),
            Err(AtropisomerError::RingAtomRowCount {
                actual: 3,
                expected: 2,
            })
        );
        assert_eq!(atom_rows_wrong, atom_snapshot, "rings mutated");
        assert_eq!(
            cleanup_invalid_atropisomers(&input, &hybridization, &bond_rows_wrong),
            Err(AtropisomerError::RingBondRowCount {
                actual: 2,
                expected: 1,
            })
        );
    }

    #[test]
    fn sanitize_ring_l04_invalid_topology_first_in_existing_order() {
        // Stale adjacency for a two-atom topology is invalid, so topology
        // validation fires before the equally-bad rings/hybridization.
        let input = cc_graph();
        let invalid = cosmolkit_model::TopologyBlock {
            adjacency: AdjacencyList::from_topology(0, &[]),
            ..input.clone()
        };
        let hybridization = hybs(&[Hybridization::Sp2, Hybridization::Sp2, Hybridization::Sp2]);
        let fast = RingInfo::new(RingFindType::Fast, 2, 1);
        assert!(matches!(
            cleanup_invalid_atropisomers(&invalid, &hybridization, &fast),
            Err(AtropisomerError::InvalidTopology { .. })
        ));
    }
}
