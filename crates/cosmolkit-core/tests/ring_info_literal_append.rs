use cosmolkit_core::{RingFindType, RingFindingError, RingInfo};
use cosmolkit_model::{AtomId, BondId};

#[test]
fn ring_literal_append_preserves_quality_rows() {
    let mut discrepancies = Vec::new();
    let mut successful_calls = 0;
    let mut rejected_calls = 0;
    let qualities = [
        (RingFindType::OtherOrUnknown, (false, false, false)),
        (RingFindType::Fast, (true, false, false)),
        (RingFindType::Sssr, (true, true, false)),
        (RingFindType::SymmSssr, (true, true, true)),
    ];
    for (quality, predicates) in qualities {
        let mut carrier = RingInfo::new(quality, 5, 5);
        for (step, row, expected_members) in [
            (
                1,
                [0, 1, 2],
                vec![vec![0], vec![0], vec![0], vec![], vec![]],
            ),
            (
                2,
                [2, 1, 0],
                vec![vec![0, 1], vec![0, 1], vec![0, 1], vec![], vec![]],
            ),
        ] {
            let atom_indices = row;
            let bond_indices = row;
            let atom_before = atom_indices;
            let bond_before = bond_indices;
            let result = carrier.add_ring(&atom_indices, &bond_indices);
            successful_calls += 1;
            if atom_indices != atom_before || bond_indices != bond_before {
                discrepancies.push(format!("{quality:?} step{step} input preservation"));
            }
            if result != Ok(step) {
                discrepancies.push(format!("{quality:?} step{step} result {result:?}"));
            }
            let expected_atoms = if step == 1 {
                vec![vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]]
            } else {
                vec![
                    vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)],
                    vec![AtomId::new(2), AtomId::new(1), AtomId::new(0)],
                ]
            };
            let expected_bonds = if step == 1 {
                vec![vec![BondId::new(0), BondId::new(1), BondId::new(2)]]
            } else {
                vec![
                    vec![BondId::new(0), BondId::new(1), BondId::new(2)],
                    vec![BondId::new(2), BondId::new(1), BondId::new(0)],
                ]
            };
            if carrier.atom_rings() != expected_atoms
                || carrier.bond_rings() != expected_bonds
                || carrier.num_rings() != step
            {
                discrepancies.push(format!("{quality:?} step{step} ordered rows/count"));
            }
            for (index, expected) in expected_members.iter().enumerate() {
                if carrier.atom_members(AtomId::new(index)) != expected
                    || carrier.bond_members(BondId::new(index)) != expected
                {
                    discrepancies.push(format!("{quality:?} step{step} membership{index}"));
                }
            }
            if !carrier.is_initialized()
                || carrier.find_type() != quality
                || carrier.atom_row_count() != 5
                || carrier.bond_row_count() != 5
                || (
                    carrier.is_find_fast_or_better(),
                    carrier.is_sssr_or_better(),
                    carrier.is_symm_sssr(),
                ) != predicates
            {
                discrepancies.push(format!("{quality:?} step{step} quality/dimensions"));
            }
        }

        let mut rejected = RingInfo::new(quality, 5, 5);
        let carrier_before = rejected.clone();
        let atom_indices = [0, 1, 2];
        let bond_indices = [0, 1];
        let atom_before = atom_indices;
        let bond_before = bond_indices;
        let result = rejected.add_ring(&atom_indices, &bond_indices);
        rejected_calls += 1;
        if atom_indices != atom_before || bond_indices != bond_before {
            discrepancies.push(format!("{quality:?} rejected input preservation"));
        }
        if result
            != Err(RingFindingError::Value {
                message: "length mismatch",
            })
        {
            discrepancies.push(format!("{quality:?} exact error {result:?}"));
        }
        if rejected != carrier_before
            || !rejected.atom_rings().is_empty()
            || !rejected.bond_rings().is_empty()
            || rejected.num_rings() != 0
            || !rejected.is_initialized()
            || rejected.find_type() != quality
            || rejected.atom_row_count() != 5
            || rejected.bond_row_count() != 5
            || (
                rejected.is_find_fast_or_better(),
                rejected.is_sssr_or_better(),
                rejected.is_symm_sssr(),
            ) != predicates
        {
            discrepancies.push(format!(
                "{quality:?} whole error carrier/empty rows/quality/dimensions"
            ));
        }
        for index in 0..5 {
            if !rejected.atom_members(AtomId::new(index)).is_empty()
                || !rejected.bond_members(BondId::new(index)).is_empty()
            {
                discrepancies.push(format!("{quality:?} rejected membership{index}"));
            }
        }
    }
    assert_eq!((successful_calls, rejected_calls), (8, 4));
    assert!(discrepancies.is_empty(), "{}", discrepancies.join("\n"));
}
