//! SMILES-specific post-chemistry stereo dispatch over detached values.

use cosmolkit_core::{
    DoubleBondStereoError, LegacyStereoError, RingFindType, RingFindingError, RingInfo,
    RingSearchParams, ValenceAssignment, ValenceError, ValenceModel, assign_legacy_stereochemistry,
    assign_valence_state_for_atom_from_parts, assign_valence_with_options_for_topology,
    clear_single_bond_directions, fast_find_rings, set_double_bond_neighbor_directions,
    symmetrized_sssr,
};
use cosmolkit_model::Conformer3D;

use crate::{SmilesParseParams, SmilesRecord};

/// Structured failures from the source SMILES post-parse stereo stage.
#[derive(Debug, thiserror::Error)]
pub enum SmilesStereoError {
    #[error(transparent)]
    Properties(#[from] cosmolkit_model::MoleculePropertyError),
    #[error("prepared valence field {field} has {actual} rows; expected {expected}")]
    ValenceRows {
        field: &'static str,
        actual: usize,
        expected: usize,
    },
    #[error("prepared ring {dimension} membership has {actual} rows; expected {expected}")]
    RingRows {
        dimension: &'static str,
        actual: usize,
        expected: usize,
    },
    #[error(transparent)]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error(transparent)]
    Directions(#[from] DoubleBondStereoError),
    #[error(transparent)]
    Assignment(#[from] LegacyStereoError),
    #[error(transparent)]
    Rings(#[from] RingFindingError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

/// Test-only acquisition counters at the ACTUAL production guard sites.
/// Thread-local, never reset; observations are relative per-call deltas.
#[cfg(test)]
pub(crate) mod ring_probe {
    use std::cell::RefCell;

    use cosmolkit_core::RingFindType;

    use cosmolkit_core::RingInfo;

    use std::cell::Cell;

    thread_local! {
        static SYMM_ACQUISITIONS: Cell<u64> = const { Cell::new(0) };
        static FAST_ACQUISITIONS: Cell<u64> = const { Cell::new(0) };
        static LEGACY_ENTRIES: Cell<u64> = const { Cell::new(0) };
        static ACQUISITION_HISTORY: RefCell<Vec<Acquisition>> =
            const { RefCell::new(Vec::new()) };
    }

    /// One ACTUAL finder-return observation: find type, both nonempty-row
    /// outer-buffer addresses and the full integral ordered state, captured
    /// at the production guard site BEFORE the result is moved into the
    /// carrier. Append-only, never reset; clones exist only under cfg(test).
    pub(crate) struct Acquisition {
        pub(crate) find_type: RingFindType,
        pub(crate) atom_rows_ptr: usize,
        pub(crate) bond_rows_ptr: usize,
        pub(crate) state: RingInfo,
    }

    pub(crate) fn record_symm() {
        SYMM_ACQUISITIONS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn record_fast() {
        FAST_ACQUISITIONS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn observe(result: &RingInfo) {
        ACQUISITION_HISTORY.with(|history| {
            history.borrow_mut().push(Acquisition {
                find_type: result.find_type(),
                atom_rows_ptr: result.atom_rings().as_ptr() as usize,
                bond_rows_ptr: result.bond_rings().as_ptr() as usize,
                state: result.clone(),
            });
        });
    }

    pub(crate) fn history_len() -> usize {
        ACQUISITION_HISTORY.with(|history| history.borrow().len())
    }

    /// The observations appended after `baseline`, in order; each is (find
    /// type, both addresses, state clone).
    pub(crate) fn history_after(baseline: usize) -> Vec<(RingFindType, usize, usize, RingInfo)> {
        ACQUISITION_HISTORY.with(|history| {
            let history = history.borrow();
            history[baseline.min(history.len())..]
                .iter()
                .map(|observation| {
                    (
                        observation.find_type,
                        observation.atom_rows_ptr,
                        observation.bond_rows_ptr,
                        observation.state.clone(),
                    )
                })
                .collect()
        })
    }

    pub(crate) fn record_legacy_entry() {
        LEGACY_ENTRIES.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn symm_acquisitions() -> u64 {
        SYMM_ACQUISITIONS.with(Cell::get)
    }

    pub(crate) fn fast_acquisitions() -> u64 {
        FAST_ACQUISITIONS.with(Cell::get)
    }

    pub(crate) fn legacy_entries() -> u64 {
        LEGACY_ENTRIES.with(Cell::get)
    }
}

// F01: the 144-call finalizer matrix. 3 graphs x 4 flag cells x 2 marker
// states x 6 incoming ring states. Counters are thread-local, never reset;
// all observations are relative per-call deltas at the ACTUAL production
// guard sites.
#[cfg(test)]
mod smiles_stereo_ring_matrix_tests {
    use super::finalize_smiles_stereo;
    use super::ring_probe as probe;
    use crate::{SmilesParseParams, SmilesRecord};
    use cosmolkit_core::{
        RingFindType, RingInfo, ValenceAssignment, ValenceModel,
        assign_valence_with_options_for_topology, fast_find_rings, find_sssr, symmetrized_sssr,
    };
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CoordinateBlock, MoleculeProperties,
        TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    fn single(id: usize, begin: usize, end: usize) -> Bond {
        Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
        )
    }

    fn g0() -> TopologyBlock {
        topology(Vec::new(), Vec::new())
    }

    fn g1() -> TopologyBlock {
        topology(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![single(0, 0, 1)],
        )
    }

    fn g2() -> TopologyBlock {
        topology(
            (0..6)
                .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
                .collect(),
            (0..6).map(|id| single(id, id, (id + 1) % 6)).collect(),
        )
    }

    fn record(topology: TopologyBlock) -> SmilesRecord {
        let mut properties = MoleculeProperties::default();
        properties.set_prop("ck-ordinary", "kept").unwrap();
        SmilesRecord {
            topology,
            coordinates: CoordinateBlock::default(),
            properties,
        }
    }

    fn params(sanitize: bool, remove_hydrogens: bool) -> SmilesParseParams {
        SmilesParseParams {
            sanitize,
            remove_hydrogens,
            ..SmilesParseParams::default()
        }
    }

    fn ring_states(topology: &TopologyBlock) -> Vec<(&'static str, Option<RingInfo>)> {
        vec![
            ("none", None),
            (
                "other-empty",
                Some(RingInfo::new(
                    RingFindType::OtherOrUnknown,
                    topology.atoms.len(),
                    topology.bonds.len(),
                )),
            ),
            (
                "fast-empty",
                Some(RingInfo::new(
                    RingFindType::Fast,
                    topology.atoms.len(),
                    topology.bonds.len(),
                )),
            ),
            ("fast", Some(fast_find_rings(topology).unwrap())),
            (
                "sssr",
                Some(find_sssr(topology, &Default::default()).unwrap()),
            ),
            (
                "symm",
                Some(symmetrized_sssr(topology, &Default::default()).unwrap()),
            ),
        ]
    }

    #[test]
    fn smiles_stereo_ring_matrix_one_hundred_forty_four_calls() {
        let graphs: [(usize, fn() -> TopologyBlock); 3] =
            [(0, g0 as fn() -> TopologyBlock), (1, g1), (2, g2)];
        let mut calls = 0usize;
        for (graph_index, build) in graphs {
            for sanitize in [false, true] {
                for remove_hydrogens in [false, true] {
                    let active = sanitize || remove_hydrogens;
                    for marker in [false, true] {
                        for state_index in 0usize..6 {
                            let label = format!(
                                "G{graph_index}/s{sanitize}/rh{remove_hydrogens}/m{marker}/st{state_index}"
                            );
                            let input = build();
                            // Exact fixture prerequisites.
                            let expected_atoms = [0usize, 2, 6][graph_index];
                            let expected_bonds = [0usize, 1, 6][graph_index];
                            assert_eq!(input.atoms.len(), expected_atoms, "{label}: atoms");
                            assert_eq!(input.bonds.len(), expected_bonds, "{label}: bonds");
                            for atom in &input.atoms {
                                assert_eq!(atom.element(), Element::C, "{label}: element");
                                assert_eq!(
                                    atom.chiral_tag(),
                                    cosmolkit_types::ChiralTag::Unspecified
                                );
                            }
                            for bond in &input.bonds {
                                assert_eq!(bond.order(), BondOrder::Single);
                                assert_eq!(bond.stereo(), cosmolkit_types::BondStereo::None);
                            }
                            let mut base = record(input.clone());
                            if marker {
                                base.properties
                                    .set_prop("_needsDetectBondStereo", "0")
                                    .unwrap();
                            }
                            base.properties
                                .set_computed_prop("_StereochemDone", "0")
                                .unwrap();
                            // Fixture valence cloned per call.
                            let fixture_valence = assign_valence_with_options_for_topology(
                                &input,
                                ValenceModel::RdkitLike,
                                false,
                            )
                            .unwrap();
                            let baseline_record = base.clone();
                            let baseline_valence = fixture_valence.clone();
                            let states = ring_states(&input);
                            let (_, ref state) = states[state_index];
                            let baseline_state = state.clone();
                            let mut prepared_rings = state.clone();
                            // Pointers of the value ACTUALLY PASSED (the
                            // moved-in carrier), captured before the call.
                            let carrier_ptrs = prepared_rings.as_ref().map(|rings| {
                                (
                                    rings.atom_rings().as_ptr() as usize,
                                    rings.bond_rings().as_ptr() as usize,
                                )
                            });
                            let mut prepared_valence = Some(fixture_valence.clone());
                            // ACTUAL pre-call baselines: the retained source
                            // record, the valence actually supplied, and the
                            // carrier actually passed (whole value + both
                            // row-buffer addresses).
                            let retained_source = base.clone();
                            let source_baseline = retained_source.clone();
                            let valence_baseline = prepared_valence.clone();
                            let symm_before = probe::symm_acquisitions();
                            let fast_before = probe::fast_acquisitions();
                            let legacy_before = probe::legacy_entries();
                            let history_before = probe::history_len();
                            let output = finalize_smiles_stereo(
                                base,
                                &params(sanitize, remove_hydrogens),
                                &mut prepared_valence,
                                &mut prepared_rings,
                            )
                            .unwrap_or_else(|error| panic!("{label}: {error:?}"));
                            calls += 1;
                            let symm_delta = probe::symm_acquisitions() - symm_before;
                            let fast_delta = probe::fast_acquisitions() - fast_before;
                            let legacy_delta = probe::legacy_entries() - legacy_before;
                            // Literal acquisition dispositions.
                            let (want_fast, want_symm) = if !active {
                                (0u64, 0u64)
                            } else if marker {
                                if state_index <= 4 { (0, 1) } else { (0, 0) }
                            } else if state_index <= 1 {
                                (1, 0)
                            } else {
                                (0, 0)
                            };
                            assert_eq!(fast_delta, want_fast, "{label}: fast");
                            assert_eq!(symm_delta, want_symm, "{label}: symm");
                            // P0 repair: the RETAINED source record compared
                            // to its own pre-call baseline after every call
                            // (the actual argument record was moved in).
                            assert_eq!(
                                retained_source, source_baseline,
                                "{label}: retained source unchanged"
                            );
                            // P0 repair: the ACTUAL prepared valence (the
                            // argument passed by &mut) compared to its own
                            // pre-call baseline; these nonchiral fixtures
                            // take no legacy H-cleanup refresh.
                            assert_eq!(
                                prepared_valence, valence_baseline,
                                "{label}: actual valence unchanged"
                            );
                            // P0 repair: actual finder-return history —
                            // exact 0/1 count, find type, and the full
                            // ordered observed state equal to the returned
                            // carrier; on nonempty acquisitions BOTH
                            // observed addresses equal the final carrier
                            // addresses (move proof).
                            let observations = probe::history_after(history_before);
                            assert_eq!(
                                observations.len() as u64,
                                want_fast + want_symm,
                                "{label}: observation count"
                            );
                            if observations.len() == 1 {
                                let (
                                    observed_type,
                                    observed_atom_ptr,
                                    observed_bond_ptr,
                                    observed_state,
                                ) = &observations[0];
                                let expected_observed = if want_fast == 1 {
                                    RingFindType::Fast
                                } else {
                                    RingFindType::SymmSssr
                                };
                                assert_eq!(
                                    *observed_type, expected_observed,
                                    "{label}: observed find type"
                                );
                                let observed_carrier = prepared_rings.as_ref().unwrap();
                                assert_eq!(
                                    observed_state, observed_carrier,
                                    "{label}: observed state equals returned carrier"
                                );
                                if !observed_carrier.atom_rings().is_empty() {
                                    assert_eq!(
                                        observed_carrier.atom_rings().as_ptr() as usize,
                                        *observed_atom_ptr,
                                        "{label}: observed atom address equals carrier"
                                    );
                                    assert_eq!(
                                        observed_carrier.bond_rings().as_ptr() as usize,
                                        *observed_bond_ptr,
                                        "{label}: observed bond address equals carrier"
                                    );
                                }
                            }
                            // Retained input baseline: supplied record value
                            // compared via the untouched clone.
                            assert_eq!(
                                baseline_record.topology.atoms.len(),
                                expected_atoms,
                                "{label}: baseline atoms"
                            );
                            assert_eq!(baseline_valence, fixture_valence, "{label}: valence");
                            assert!(output.coordinates.conformers_2d.is_empty());
                            assert!(output.coordinates.conformers_3d.is_empty());
                            assert_eq!(
                                output.properties.prop("ck-ordinary"),
                                Some(&cosmolkit_model::PropertyValue::String("kept".into())),
                                "{label}: sentinel"
                            );
                            if active {
                                assert_eq!(
                                    output.properties.prop("_StereochemDone"),
                                    Some(&cosmolkit_model::PropertyValue::Int(1)),
                                    "{label}: done computed"
                                );
                                assert!(
                                    output.properties.prop("_needsDetectBondStereo").is_none(),
                                    "{label}: marker cleared"
                                );
                                assert!(legacy_delta >= 1, "{label}: legacy entry");
                                assert_eq!(legacy_delta, 1, "{label}: exactly one legacy");
                            } else {
                                assert_eq!(legacy_delta, 0, "{label}: no legacy");
                                assert_eq!(
                                    output.properties.prop("_StereochemDone"),
                                    Some(&cosmolkit_model::PropertyValue::String("0".into())),
                                    "{label}: done retained"
                                );
                                if marker {
                                    assert_eq!(
                                        output.properties.prop("_needsDetectBondStereo"),
                                        Some(&cosmolkit_model::PropertyValue::String("0".into())),
                                        "{label}: marker retained"
                                    );
                                }
                            }
                            if !active {
                                // BOTH-false: exact supplied state preserved.
                                assert_eq!(prepared_rings, baseline_state, "{label}: state");
                                if let Some(rings) = prepared_rings.as_ref() {
                                    if !rings.atom_rings().is_empty() {
                                        assert_eq!(
                                            rings.atom_rings().as_ptr() as usize,
                                            carrier_ptrs.unwrap().0,
                                            "{label}: atom pointer"
                                        );
                                        assert_eq!(
                                            rings.bond_rings().as_ptr() as usize,
                                            carrier_ptrs.unwrap().1,
                                            "{label}: bond pointer"
                                        );
                                    }
                                }
                            }
                            // Every ZERO-ACQUISITION arm: the sufficient
                            // carrier is EXACTLY reused — whole value,
                            // ordered paired rows, memberships, type and
                            // dimensions equal the actual pre-call carrier,
                            // and BOTH nonempty outer-buffer pointers are
                            // preserved (move identity, not absence of
                            // finder-internal allocations).
                            if want_fast + want_symm == 0 {
                                assert_eq!(
                                    prepared_rings, baseline_state,
                                    "{label}: zero-acquisition carrier reused"
                                );
                                if let Some(rings) = prepared_rings.as_ref() {
                                    if !rings.atom_rings().is_empty() {
                                        assert_eq!(
                                            rings.atom_rings().as_ptr() as usize,
                                            carrier_ptrs.unwrap().0,
                                            "{label}: reuse atom pointer"
                                        );
                                        assert_eq!(
                                            rings.bond_rings().as_ptr() as usize,
                                            carrier_ptrs.unwrap().1,
                                            "{label}: reuse bond pointer"
                                        );
                                    }
                                }
                            }
                            // Quality/dimension checks for the final carrier.
                            if active {
                                let rings = prepared_rings.as_ref().unwrap();
                                // Frozen final RingFindType table.
                                let expected_type = if marker {
                                    RingFindType::SymmSssr
                                } else if state_index <= 3 {
                                    RingFindType::Fast
                                } else if state_index == 4 {
                                    RingFindType::Sssr
                                } else {
                                    RingFindType::SymmSssr
                                };
                                assert_eq!(
                                    rings.find_type(),
                                    expected_type,
                                    "{label}: final find type"
                                );
                                assert!(rings.is_initialized(), "{label}: initialized");
                                assert_eq!(
                                    rings.atom_row_count(),
                                    expected_atoms,
                                    "{label}: atom dims"
                                );
                                assert_eq!(
                                    rings.bond_row_count(),
                                    expected_bonds,
                                    "{label}: bond dims"
                                );
                                if graph_index == 2 {
                                    let acquired = fast_delta + symm_delta > 0;
                                    let supplied_nonempty = baseline_state
                                        .as_ref()
                                        .is_some_and(|state| !state.atom_rings().is_empty());
                                    if acquired || supplied_nonempty {
                                        // Real rows: one ring, normalized
                                        // [0..5], all members[0].
                                        assert_eq!(
                                            rings.atom_rings().len(),
                                            1,
                                            "{label}: ring count"
                                        );
                                        let mut atoms_row: Vec<usize> = rings.atom_rings()[0]
                                            .iter()
                                            .map(|atom| atom.index())
                                            .collect();
                                        atoms_row.sort_unstable();
                                        assert_eq!(
                                            atoms_row,
                                            vec![0, 1, 2, 3, 4, 5],
                                            "{label}: normalized atoms"
                                        );
                                        // P0 repair: normalized paired BOND
                                        // row for the real six-ring.
                                        assert_eq!(
                                            rings.bond_rings().len(),
                                            1,
                                            "{label}: bond ring count"
                                        );
                                        let mut bonds_row: Vec<usize> = rings.bond_rings()[0]
                                            .iter()
                                            .map(|bond| bond.index())
                                            .collect();
                                        bonds_row.sort_unstable();
                                        assert_eq!(
                                            bonds_row,
                                            vec![0, 1, 2, 3, 4, 5],
                                            "{label}: normalized bonds"
                                        );
                                    } else {
                                        // The source quality guard TRUSTS the
                                        // supplied initialized-empty rows.
                                        assert!(
                                            rings.atom_rings().is_empty(),
                                            "{label}: supplied empty trusted"
                                        );
                                        // P0 repair: empty paired bond rows
                                        // and empty member entries on the
                                        // trusted empty-supplied branch.
                                        assert!(
                                            rings.bond_rings().is_empty(),
                                            "{label}: empty bond rows"
                                        );
                                        for index in 0..6usize {
                                            assert_eq!(
                                                rings.atom_members(AtomId::new(index)),
                                                &[] as &[usize],
                                                "{label}: empty atom member {index}"
                                            );
                                            assert_eq!(
                                                rings.bond_members(BondId::new(index)),
                                                &[] as &[usize],
                                                "{label}: empty bond member {index}"
                                            );
                                        }
                                    }
                                    if acquired || supplied_nonempty {
                                        for index in 0..6usize {
                                            assert_eq!(
                                                rings.atom_members(AtomId::new(index)),
                                                &[0],
                                                "{label}: member {index}"
                                            );
                                            // P0 repair: every bond member.
                                            assert_eq!(
                                                rings.bond_members(BondId::new(index)),
                                                &[0],
                                                "{label}: bond member {index}"
                                            );
                                        }
                                    }
                                } else {
                                    assert!(rings.atom_rings().is_empty(), "{label}: empty rows");
                                    // P0 repair: empty paired bond rows;
                                    // G1 also proves empty member entries.
                                    assert!(
                                        rings.bond_rings().is_empty(),
                                        "{label}: empty bond rows"
                                    );
                                    if graph_index == 1 {
                                        assert_eq!(
                                            rings.atom_members(AtomId::new(0)),
                                            &[] as &[usize],
                                            "{label}: empty atom member"
                                        );
                                        assert_eq!(
                                            rings.bond_members(BondId::new(0)),
                                            &[] as &[usize],
                                            "{label}: empty bond member"
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 144, "exact census");
    }
}

// F02: eight exact prepared-state failures on G1 with active flags
// sanitize=true/remove=false. Four ring-dimension cases (sufficient-quality
// carriers actually reused) and four valence cases (explicit/implicit rows3,
// rings None). Source error order preserved; no legacy dispatch on error.
#[cfg(test)]
mod smiles_stereo_ring_errors_tests {
    use super::finalize_smiles_stereo;
    use super::ring_probe as probe;
    use crate::{SmilesParseParams, SmilesRecord};
    use cosmolkit_core::{
        RingFindType, RingInfo, ValenceAssignment, ValenceModel,
        assign_valence_with_options_for_topology,
    };
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CoordinateBlock, MoleculeProperties,
        TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    fn g1() -> TopologyBlock {
        TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            )],
            Vec::new(),
            Vec::new(),
        )
        .unwrap()
    }

    fn record(marker: bool) -> SmilesRecord {
        let mut properties = MoleculeProperties::default();
        properties.set_prop("ck-ordinary", "kept").unwrap();
        if marker {
            properties.set_prop("_needsDetectBondStereo", "0").unwrap();
        }
        SmilesRecord {
            topology: g1(),
            coordinates: CoordinateBlock::default(),
            properties,
        }
    }

    fn params() -> SmilesParseParams {
        SmilesParseParams {
            sanitize: true,
            remove_hydrogens: false,
            ..SmilesParseParams::default()
        }
    }

    fn valence() -> ValenceAssignment {
        assign_valence_with_options_for_topology(&g1(), ValenceModel::RdkitLike, false).unwrap()
    }

    #[test]
    fn smiles_stereo_ring_errors_eight_call_product() {
        let mut calls = 0usize;
        // --- Four ring-dimension cases ---------------------------------
        for (case, atom_rows, bond_rows, marker) in [
            ("atom3/bond1", 3usize, 1usize, false),
            ("atom2/bond2", 2, 2, false),
            ("atom3/bond1/marker", 3, 1, true),
            ("atom2/bond2/marker", 2, 2, true),
        ] {
            let label = case.to_string();
            let input = record(marker);
            // Retained source record: whole value, compared to its own
            // fresh pre-call baseline (not a fixture substitute).
            let retained_source = input.clone();
            let source_baseline = retained_source.clone();
            let fixture_valence = valence();
            // Construct malformed-dimension carriers through the public
            // constructor with mismatched counts.
            let malformed = if marker {
                RingInfo::new(RingFindType::SymmSssr, atom_rows, bond_rows)
            } else {
                RingInfo::new(RingFindType::Fast, atom_rows, bond_rows)
            };
            // ACTUAL pre-call carrier and valence baselines.
            let carrier_baseline = malformed.clone();
            let mut prepared_valence = Some(fixture_valence.clone());
            let valence_baseline = prepared_valence.clone();
            let mut prepared_rings = Some(malformed);
            let symm_before = probe::symm_acquisitions();
            let fast_before = probe::fast_acquisitions();
            let legacy_before = probe::legacy_entries();
            let history_before = probe::history_len();
            let error = finalize_smiles_stereo(
                input,
                &params(),
                &mut prepared_valence,
                &mut prepared_rings,
            )
            .unwrap_err();
            calls += 1;
            let expected_dimension = if atom_rows != 2 { "atom" } else { "bond" };
            let (expected_actual, expected_expected) = if atom_rows != 2 {
                (3usize, 2usize)
            } else {
                (2, 1)
            };
            assert!(
                matches!(
                    &error,
                    super::SmilesStereoError::RingRows {
                        dimension,
                        actual,
                        expected
                    } if *dimension == expected_dimension
                        && *actual == expected_actual
                        && *expected == expected_expected
                ),
                "{label}: got {error:?}"
            );
            // Zero acquisitions; no legacy dispatch.
            assert_eq!(probe::symm_acquisitions() - symm_before, 0, "{label}: symm");
            assert_eq!(probe::fast_acquisitions() - fast_before, 0, "{label}: fast");
            assert_eq!(
                probe::legacy_entries() - legacy_before,
                0,
                "{label}: legacy"
            );
            // No finder-return observations were appended.
            assert_eq!(
                probe::history_after(history_before).len(),
                0,
                "{label}: no observations"
            );
            // Retained record: whole value unchanged after the error.
            assert_eq!(
                retained_source, source_baseline,
                "{label}: retained record unchanged"
            );
            // ACTUAL prepared_valence unchanged after the error.
            assert_eq!(
                prepared_valence, valence_baseline,
                "{label}: actual valence unchanged"
            );
            // ACTUAL carrier remains EXACTLY the malformed input.
            assert_eq!(
                prepared_rings,
                Some(carrier_baseline),
                "{label}: carrier unchanged"
            );
        }
        // --- Four valence cases ---------------------------------------
        for (case, field, marker) in [
            ("explicit3", "explicit_valence", false),
            ("implicit3", "implicit_hydrogens", false),
            ("explicit3/marker", "explicit_valence", true),
            ("implicit3/marker", "implicit_hydrogens", true),
        ] {
            let label = case.to_string();
            let input = record(marker);
            let retained_source = input.clone();
            let source_baseline = retained_source.clone();
            let mut malformed = valence();
            if field == "explicit_valence" {
                malformed.explicit_valence.resize(3, 0);
            } else {
                malformed.implicit_hydrogens.resize(3, 0);
            }
            let mut prepared_valence = Some(malformed);
            let valence_baseline = prepared_valence.clone();
            let mut prepared_rings: Option<RingInfo> = None;
            let symm_before = probe::symm_acquisitions();
            let fast_before = probe::fast_acquisitions();
            let legacy_before = probe::legacy_entries();
            let history_before = probe::history_len();
            let error = finalize_smiles_stereo(
                input,
                &params(),
                &mut prepared_valence,
                &mut prepared_rings,
            )
            .unwrap_err();
            calls += 1;
            assert!(
                matches!(
                    &error,
                    super::SmilesStereoError::ValenceRows {
                        field: observed,
                        actual,
                        expected
                    } if *observed == field && *actual == 3 && *expected == 2
                ),
                "{label}: got {error:?}"
            );
            assert_eq!(
                probe::legacy_entries() - legacy_before,
                0,
                "{label}: legacy"
            );
            if marker {
                // Marker-direction Symm acquisition happens BEFORE the
                // valence checks; exactly one.
                assert_eq!(probe::symm_acquisitions() - symm_before, 1, "{label}: symm");
                // Exactly ONE actual finder-return observation; its ordered
                // state equals the returned carrier (no atomicity claim for
                // this branch — rings were replaced before the valence
                // error).
                let observations = probe::history_after(history_before);
                assert_eq!(observations.len(), 1, "{label}: one observation");
                let rings = prepared_rings.as_ref().unwrap();
                assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{label}");
                assert!(rings.atom_rings().is_empty(), "{label}: empty");
                assert!(rings.is_initialized(), "{label}: initialized");
                assert_eq!(rings.atom_row_count(), 2, "{label}: carrier atom dims 2");
                assert_eq!(rings.bond_row_count(), 1, "{label}: carrier bond dims 1");
                assert!(rings.bond_rings().is_empty(), "{label}: bond rows empty");
                assert_eq!(
                    rings.atom_members(cosmolkit_model::AtomId::new(0)),
                    &[] as &[usize],
                    "{label}: atom0 membership empty"
                );
                assert_eq!(
                    rings.bond_members(cosmolkit_model::BondId::new(0)),
                    &[] as &[usize],
                    "{label}: bond0 membership empty"
                );
                let (_, observed_atom_ptr, observed_bond_ptr, observed) = &observations[0];
                assert_eq!(rings, observed, "{label}: observed state equals carrier");
                assert_eq!(
                    rings.atom_rings().as_ptr() as usize,
                    *observed_atom_ptr,
                    "{label}: observed atom address moved"
                );
                assert_eq!(
                    rings.bond_rings().as_ptr() as usize,
                    *observed_bond_ptr,
                    "{label}: observed bond address moved"
                );
            } else {
                assert_eq!(probe::symm_acquisitions() - symm_before, 0, "{label}: symm");
                assert_eq!(
                    probe::history_after(history_before).len(),
                    0,
                    "{label}: no observations"
                );
                assert!(prepared_rings.is_none(), "{label}: rings none");
            }
            assert_eq!(probe::fast_acquisitions() - fast_before, 0, "{label}: fast");
            // Retained record and actual malformed valence unchanged.
            assert_eq!(
                retained_source, source_baseline,
                "{label}: retained record unchanged"
            );
            assert_eq!(
                prepared_valence, valence_baseline,
                "{label}: actual valence unchanged"
            );
        }
        assert_eq!(calls, 8, "exact census");
    }
}

// F03: sixteen source-ordered preparation calls. Real parse -> (RemoveHs if
// remove | else SAN if sanitize) -> finalizer, moving the actual result
// fields. Final quality literal table incl. bothtrue EMPTY Fast (RH's
// original-count guard skips SAN on the empty graph).
#[cfg(test)]
mod smiles_stereo_ring_routes_tests {
    use super::finalize_smiles_stereo;
    use super::ring_probe as probe;
    use crate::{SmilesParseParams, SmilesRecord, parse_smiles};
    use cosmolkit_core::{RingFindType, SanitizeParams, sanitize_topology};
    use cosmolkit_model::AtomId;
    use cosmolkit_model::CoordinateBlock;
    use cosmolkit_types::Element;

    fn profile(sanitize: bool, remove_hydrogens: bool) -> SmilesParseParams {
        SmilesParseParams {
            sanitize,
            remove_hydrogens,
            ..SmilesParseParams::default()
        }
    }

    #[test]
    fn smiles_stereo_ring_routes_sixteen_call_product() {
        let inputs = ["", "CC", "c1ccccc1", "[H]C1CCCCC1"];
        let mut calls = 0usize;
        for input in inputs {
            for (name, sanitize, remove_hydrogens) in [
                ("bothfalse", false, false),
                ("remove-only", false, true),
                ("sanitize-only", true, false),
                ("bothtrue", true, true),
            ] {
                let label = format!("{name}/{input:?}");
                let parse_params = profile(sanitize, remove_hydrogens);
                let parsed = parse_smiles(input, &parse_params).unwrap();
                // Source-ordered preparation, exactly as the constructor.
                let (topology, coordinates, properties, mut prepared_valence, mut prepared_rings) =
                    if remove_hydrogens {
                        let remove_params = cosmolkit_core::RemoveHsParams {
                            update_explicit_count: true,
                            sanitize: parse_params.sanitize,
                            ..cosmolkit_core::RemoveHsParams::default()
                        };
                        let result = cosmolkit_core::remove_hydrogens_with_params(
                            parsed.topology,
                            parsed.coordinates,
                            parsed.properties,
                            &remove_params,
                        )
                        .unwrap();
                        (
                            result.topology,
                            result.coordinates,
                            result.properties,
                            result.final_valence,
                            result.final_rings,
                        )
                    } else if sanitize {
                        let result =
                            sanitize_topology(&parsed.topology, &SanitizeParams::default())
                                .unwrap();
                        let cosmolkit_core::SanitizeAssignment {
                            topology,
                            final_valence,
                            final_rings,
                            ..
                        } = result;
                        let properties = {
                            let mut props = parsed.properties.clone();
                            props.clear_computed_props();
                            props
                        };
                        (
                            topology,
                            parsed.coordinates,
                            properties,
                            final_valence,
                            final_rings,
                        )
                    } else {
                        (
                            parsed.topology,
                            parsed.coordinates,
                            parsed.properties,
                            None,
                            None,
                        )
                    };
                let record = SmilesRecord {
                    topology,
                    coordinates,
                    properties,
                };
                // Retained prepared source-record baseline and actual
                // pre-call carrier baselines (whole value + both nonempty
                // row-buffer addresses).
                let retained_source = record.clone();
                let source_baseline = retained_source.clone();
                let valence_baseline = prepared_valence.clone();
                let carrier_baseline = prepared_rings.clone();
                let carrier_ptrs = prepared_rings.as_ref().map(|rings| {
                    (
                        rings.atom_rings().as_ptr() as usize,
                        rings.bond_rings().as_ptr() as usize,
                    )
                });
                let symm_before = probe::symm_acquisitions();
                let fast_before = probe::fast_acquisitions();
                let history_before = probe::history_len();
                let output = finalize_smiles_stereo(
                    record,
                    &parse_params,
                    &mut prepared_valence,
                    &mut prepared_rings,
                )
                .unwrap_or_else(|error| panic!("{label}: {error:?}"));
                calls += 1;
                let symm_delta = probe::symm_acquisitions() - symm_before;
                let fast_delta = probe::fast_acquisitions() - fast_before;
                let observations = probe::history_after(history_before);
                // Retained prepared record and coordinates unchanged.
                assert_eq!(
                    retained_source, source_baseline,
                    "{label}: retained record unchanged"
                );
                assert_eq!(
                    output.coordinates, retained_source.coordinates,
                    "{label}: coordinates unchanged"
                );
                // Actual supplied Some(valence) unchanged; an incoming None
                // is NOT falsely required to remain None on active arms
                // (the finalizer prepares it nonstrictly).
                if valence_baseline.is_some() {
                    assert_eq!(
                        prepared_valence, valence_baseline,
                        "{label}: actual valence unchanged"
                    );
                }
                // Frozen final-quality literal table.
                match name {
                    "bothfalse" => {
                        assert!(prepared_rings.is_none(), "{label}");
                        assert_eq!(symm_delta, 0, "{label}: no symm");
                        assert_eq!(fast_delta, 0, "{label}: no fast");
                        assert_eq!(observations.len(), 0, "{label}: no obs");
                    }
                    "remove-only" => {
                        let rings = prepared_rings.as_ref().unwrap();
                        assert_eq!(rings.find_type(), RingFindType::Fast, "{label}");
                        // remove-only has sanitize=false, so RH's nested
                        // sanitize never fires and final_rings arrives None
                        // for ALL inputs: the finalizer's Fast guard
                        // acquires exactly ONE canonical Fast.
                        assert_eq!(symm_delta, 0, "{label}: no symm");
                        assert_eq!(fast_delta, 1, "{label}: one fast");
                        assert_eq!(observations.len(), 1, "{label}: one obs");
                    }
                    "sanitize-only" => {
                        let rings = prepared_rings.as_ref().unwrap();
                        assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{label}");
                        assert_eq!(symm_delta, 0, "{label}: no symm");
                        assert_eq!(fast_delta, 0, "{label}: no fast");
                        assert_eq!(observations.len(), 0, "{label}: no obs");
                        assert_eq!(
                            prepared_rings, carrier_baseline,
                            "{label}: sufficient carrier reused"
                        );
                    }
                    _ => {
                        let rings = prepared_rings.as_ref().unwrap();
                        if input.is_empty() {
                            // RH original-count guard skips SAN on EMPTY.
                            assert_eq!(rings.find_type(), RingFindType::Fast, "{label}");
                            // Empty RH run leaves Fast quality; exactly one
                            // Symm acquisition at the marker-less direction
                            // position does not occur (no marker), so the
                            // Fast guard may borrow or acquire once.
                            if carrier_baseline
                                .as_ref()
                                .is_some_and(|base| base.is_find_fast_or_better())
                            {
                                assert_eq!(symm_delta, 0, "{label}: no symm");
                                assert_eq!(fast_delta, 0, "{label}: no fast");
                                assert_eq!(observations.len(), 0, "{label}: no obs");
                                assert_eq!(
                                    prepared_rings, carrier_baseline,
                                    "{label}: sufficient carrier reused"
                                );
                            } else {
                                assert_eq!(symm_delta, 0, "{label}: no symm");
                                assert_eq!(fast_delta, 1, "{label}: one fast");
                                assert_eq!(observations.len(), 1, "{label}: one obs");
                            }
                        } else {
                            assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{label}");
                            // SAN already supplied Symm: reused unchanged.
                            assert_eq!(symm_delta, 0, "{label}: no symm");
                            assert_eq!(fast_delta, 0, "{label}: no fast");
                            assert_eq!(observations.len(), 0, "{label}: no obs");
                            assert_eq!(
                                prepared_rings, carrier_baseline,
                                "{label}: sufficient carrier reused"
                            );
                        }
                    }
                }
                // Final dimensions match the final topology.
                if let Some(rings) = prepared_rings.as_ref() {
                    assert!(rings.is_initialized(), "{label}: initialized");
                    assert_eq!(
                        rings.atom_row_count(),
                        output.topology.atoms.len(),
                        "{label}: atom dims"
                    );
                    assert_eq!(
                        rings.bond_row_count(),
                        output.topology.bonds.len(),
                        "{label}: bond dims"
                    );
                }
                // _StereochemDone disposition.
                if sanitize || remove_hydrogens {
                    assert_eq!(
                        output.properties.prop("_StereochemDone"),
                        Some(&cosmolkit_model::PropertyValue::Int(1)),
                        "{label}"
                    );
                }
                // Literal final topology counts.
                let (want_atoms, want_bonds) = match input {
                    "" => (0usize, 0usize),
                    "CC" => (2, 1),
                    "c1ccccc1" => (6, 6),
                    _ => {
                        if remove_hydrogens {
                            (6, 6)
                        } else {
                            (7, 7)
                        }
                    }
                };
                assert_eq!(output.topology.atoms.len(), want_atoms, "{label}: atoms");
                assert_eq!(output.topology.bonds.len(), want_bonds, "{label}: bonds");
                // Literal final atom identity against the retained
                // original-to-final mapping: original H0 removed only with
                // remove=true; retained C1..6 become C0..5 then, else stay
                // C1..6; elements and (absent) isotopes preserved.
                match input {
                    "" | "CC" => {}
                    "c1ccccc1" => {
                        for index in 0..6usize {
                            assert_eq!(
                                output.topology.atoms[index].element(),
                                Element::C,
                                "{label}: element {index}"
                            );
                            assert!(
                                output.topology.atoms[index].isotope().is_none(),
                                "{label}: isotope {index}"
                            );
                        }
                    }
                    _ => {
                        if remove_hydrogens {
                            for index in 0..6usize {
                                assert_eq!(
                                    output.topology.atoms[index].element(),
                                    Element::C,
                                    "{label}: element {index}"
                                );
                                assert!(
                                    output.topology.atoms[index].isotope().is_none(),
                                    "{label}: isotope {index}"
                                );
                            }
                        } else {
                            assert_eq!(output.topology.atoms[0].element(), Element::H);
                            for index in 1..7usize {
                                assert_eq!(
                                    output.topology.atoms[index].element(),
                                    Element::C,
                                    "{label}: element {index}"
                                );
                                assert!(
                                    output.topology.atoms[index].isotope().is_none(),
                                    "{label}: isotope {index}"
                                );
                            }
                        }
                    }
                }
                // Cycle normalized paired sets and every paired membership
                // entry, including H0 empty; acyclic paired rows empty.
                if let Some(rings) = prepared_rings.as_ref() {
                    match input {
                        "" | "CC" => {
                            assert!(rings.atom_rings().is_empty(), "{label}: rows");
                            assert!(rings.bond_rings().is_empty(), "{label}: bond rows");
                        }
                        "c1ccccc1" => {
                            let mut atoms_row: Vec<usize> = rings.atom_rings()[0]
                                .iter()
                                .map(|atom| atom.index())
                                .collect();
                            atoms_row.sort_unstable();
                            let mut bonds_row: Vec<usize> = rings.bond_rings()[0]
                                .iter()
                                .map(|bond| bond.index())
                                .collect();
                            bonds_row.sort_unstable();
                            assert_eq!(atoms_row, vec![0, 1, 2, 3, 4, 5], "{label}: cycle atoms");
                            assert_eq!(bonds_row, vec![0, 1, 2, 3, 4, 5], "{label}: cycle bonds");
                            for index in 0..6usize {
                                assert_eq!(
                                    rings.atom_members(AtomId::new(index)),
                                    &[0],
                                    "{label}: atom member {index}"
                                );
                                assert_eq!(
                                    rings.bond_members(cosmolkit_model::BondId::new(index)),
                                    &[0],
                                    "{label}: bond member {index}"
                                );
                            }
                        }
                        _ => {
                            let expected: Vec<usize> = if remove_hydrogens {
                                vec![0, 1, 2, 3, 4, 5]
                            } else {
                                vec![1, 2, 3, 4, 5, 6]
                            };
                            let mut atoms_row: Vec<usize> = rings.atom_rings()[0]
                                .iter()
                                .map(|atom| atom.index())
                                .collect();
                            atoms_row.sort_unstable();
                            let mut bonds_row: Vec<usize> = rings.bond_rings()[0]
                                .iter()
                                .map(|bond| bond.index())
                                .collect();
                            bonds_row.sort_unstable();
                            assert_eq!(atoms_row, expected, "{label}: cycle atoms");
                            let expected_bonds: Vec<usize> = if remove_hydrogens {
                                vec![0, 1, 2, 3, 4, 5]
                            } else {
                                vec![1, 2, 3, 4, 5, 6]
                            };
                            assert_eq!(bonds_row, expected_bonds, "{label}: cycle bonds");
                            let upper = if remove_hydrogens { 6 } else { 7 };
                            for index in 0..upper {
                                let expected_members: &[usize] = if !remove_hydrogens && index == 0
                                {
                                    &[]
                                } else {
                                    &[0]
                                };
                                assert_eq!(
                                    rings.atom_members(AtomId::new(index)),
                                    expected_members,
                                    "{label}: atom member {index}"
                                );
                                assert_eq!(
                                    rings.bond_members(cosmolkit_model::BondId::new(index)),
                                    expected_members,
                                    "{label}: bond member {index}"
                                );
                            }
                        }
                    }
                }
                // Nonempty SUFFICIENT incoming carriers keep both addresses
                // (reuse); freshly acquired carriers have no pre-call
                // baseline to preserve.
                if let Some((baseline_atom_ptr, baseline_bond_ptr)) = carrier_ptrs {
                    assert_eq!(
                        prepared_rings.as_ref().unwrap().atom_rings().as_ptr() as usize,
                        baseline_atom_ptr,
                        "{label}: atom address preserved"
                    );
                    assert_eq!(
                        prepared_rings.as_ref().unwrap().bond_rings().as_ptr() as usize,
                        baseline_bond_ptr,
                        "{label}: bond address preserved"
                    );
                }
                let _ = (symm_delta, fast_delta);
                let _ = CoordinateBlock::default();
            }
        }
        assert_eq!(calls, 16, "exact census");
    }
}

/// Complete stereo after the caller has performed the requested source
/// sanitize/RemoveHs stage. Never accepts or constructs a live molecule.
pub fn finalize_smiles_stereo(
    mut record: SmilesRecord,
    params: &SmilesParseParams,
    prepared_valence: &mut Option<ValenceAssignment>,
    prepared_rings: &mut Option<RingInfo>,
) -> Result<SmilesRecord, SmilesStereoError> {
    // RDKit SmilesParse.cpp, MolFromSmiles (2026.03.1):
    // RDKit✔️✔️:   if (res && (params.sanitize || params.removeHs)) {
    // The canonical core owners have already run the preceding removeHs or
    // sanitizeMol branch. Continue with their final topology and coordinates.
    // RDKit✔️✔️:     if (res->hasProp(SmilesParseOps::detail::_needsDetectBondStereo)) {
    // RDKit✔️✔️:       // we encountered either wiggly bond in the CXSMILES,
    // RDKit✔️✔️:       // these need to be handled the same way they were in mol files
    // RDKit✔️✔️:       if (conf || conf3d) {
    // RDKit✔️✔️:         MolOps::clearSingleBondDirFlags(*res);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       MolOps::setDoubleBondNeighborDirections(*res, conf ? conf : conf3d);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res->clearProp(SmilesParseOps::detail::_needsDetectBondStereo);
    // RDKit✔️✔️:     // figure out stereochemistry:
    // RDKit✔️✔️:     bool cleanIt = true, force = true, flagPossible = true;
    // RDKit✔️✔️:     MolOps::assignStereochemistry(*res, cleanIt, force, flagPossible);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     //  we still need to do something about double bond stereochemistry
    // RDKit✔️✔️:     //  (was github issue 337)
    // RDKit✔️✔️:     //  now that atom stereochem has been perceived, the wedging
    // RDKit✔️✔️:     //  information is no longer needed, so we clear
    // RDKit✔️✔️:     //  single bond dir flags:
    // RDKit✔️✔️:     MolOps::clearSingleBondDirFlags(*res, true);
    // RDKit✔️✔️:   }
    // Behavior: retain the marker when BOTH flags are false. Otherwise consume
    // it only after successful direction reconstruction, then run the existing
    // fixed-profile legacy assignment (not merely Cis/Trans -> Z/E relabeling).
    // Complexity: dispatch reuses the unique core algorithms. Borrow an already
    // computed assignment; only a cold call computes all rows. Detached outputs
    // carry no runtime cache authority.
    // The optional XY lift below costs O(V); no extra topology clone is needed
    // at this boundary. Errors discard the owned record without live mutation.
    if !params.sanitize && !params.remove_hydrogens {
        record.topology = clear_single_bond_directions(record.topology, true)?;
        return Ok(record);
    }

    // Ring-state quality is now guard-driven at the exact source positions
    // below; the previous unconditional eager symmetrized_sssr is removed.
    if record.properties.prop("_needsDetectBondStereo").is_some() {
        // Complete pinned source: Chirality.cpp setDoubleBondNeighborDirections
        // ring-quality guard, applied at the actual dispatch position.
        // RDKit✔️✔️:   if (!mol.getRingInfo()->isSymmSssr()) {
        // RDKit✔️✔️:     RDKit::MolOps::symmetrizeSSSR(mol);
        // RDKit✔️✔️:   }
        // Behavior review: the direction owner requires Symm even when no
        // candidate double bond exists and on acyclic/empty graphs. None,
        // reset, Other, Fast and Sssr are REPLACED by one canonical
        // symmetrized_sssr MOVED into the carrier; SymmSssr (including
        // initialized-empty) is borrowed unchanged. Cost review: exactly one
        // acquisition per insufficient state (the source symmetrizeSSSR);
        // sufficient states add no find, clone or row copy.
        let needs_symm = !prepared_rings
            .as_ref()
            .is_some_and(|carrier| carrier.is_symm_sssr());
        if needs_symm {
            let acquired = symmetrized_sssr(&record.topology, &RingSearchParams::default())?;
            #[cfg(test)]
            ring_probe::record_symm();
            #[cfg(test)]
            ring_probe::observe(&acquired);
            *prepared_rings = Some(acquired);
        }
        let rings = prepared_rings
            .as_ref()
            .expect("symmetrized carrier present at the direction owner");
        // Narrow prepared-state safety, atom-then-bond, immediately before
        // the core ring read. Core RingInfo's checked constructors own
        // index/member invariants; only membership DIMENSIONS are checked
        // here, never same-sized stale rows or num_rings()>0 inference.
        let (atom_rows, bond_rows) = (rings.atom_row_count(), rings.bond_row_count());
        if atom_rows != record.topology.atoms.len() {
            return Err(SmilesStereoError::RingRows {
                dimension: "atom",
                actual: atom_rows,
                expected: record.topology.atoms.len(),
            });
        }
        if bond_rows != record.topology.bonds.len() {
            return Err(SmilesStereoError::RingRows {
                dimension: "bond",
                actual: bond_rows,
                expected: record.topology.bonds.len(),
            });
        }
        // RDKit✔️✔️:       if (!testConf->is3D()) {
        // RDKit✔️✔️:         if (conf == nullptr) {  // only take the first 2d conf
        // RDKit✔️✔️:           conf = testConf;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         if (conf3d == nullptr) {  // only take the first 3d conf
        // RDKit✔️✔️:           conf3d = testConf;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit's geometry kernel accepts XYZ also for a 2D conformer. Lift
        // only a borrowed 2D input, prefer it over 3D, and keep stored rows intact.
        let (two_d, three_d) = source_stereo_conformers(&record.coordinates)?;
        let conformer = two_d.as_deref().or(three_d);
        if conformer.is_some() {
            record.topology = clear_single_bond_directions(record.topology, false)?;
        }
        let update = set_double_bond_neighbor_directions(record.topology, rings, conformer)?;
        record.topology = update.topology;
        if update.needs_detect_bond_stereo {
            // RDKit❗✔️:     mol.setProp("_needsDetectBondStereo", 1);
            record
                .properties
                .set_prop("_needsDetectBondStereo", 1_i32)?;
        }
    }
    record.properties.clear_prop("_needsDetectBondStereo")?;
    // assignStereochemistry updates a missing property cache non-strictly.
    // This local assignment is not a claim that unsanitized chemistry passed
    // strict sanitization, nor does it undo CK-VALENCE-001 runtime invalidation.
    if prepared_valence.is_none() {
        *prepared_valence = Some(assign_valence_with_options_for_topology(
            &record.topology,
            ValenceModel::RdkitLike,
            false,
        )?);
    }
    let valence = prepared_valence.as_mut().expect("assigned above");
    // Validate both fields before any core stereo path indexes a prepared row.
    for (field, actual) in [
        ("explicit_valence", valence.explicit_valence.len()),
        ("implicit_hydrogens", valence.implicit_hydrogens.len()),
    ] {
        if actual != record.topology.atoms.len() {
            return Err(SmilesStereoError::ValenceRows {
                field,
                actual,
                expected: record.topology.atoms.len(),
            });
        }
    }
    // Complete pinned source: Chirality.cpp legacyStereoPerception ring guard,
    // applied immediately BEFORE the legacy dispatch (after the valence checks
    // above, preserving the source error order).
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // Behavior review: None/reset/Other are REPLACED by one canonical
    // fast_find_rings MOVED into the carrier; Fast/Sssr/Symm (including
    // initialized-empty) are borrowed unchanged, so the marker path above
    // finds nothing a second time. assign_legacy_stereochemistry's own Fast
    // guard then reuses these rows. Cost review: exactly one acquisition per
    // insufficient state (the source fastFindRings); sufficient states add
    // no find, clone or unconditional final symmetrization.
    let needs_fast = !prepared_rings
        .as_ref()
        .is_some_and(|carrier| carrier.is_find_fast_or_better());
    if needs_fast {
        let acquired = fast_find_rings(&record.topology)?;
        #[cfg(test)]
        ring_probe::record_fast();
        #[cfg(test)]
        ring_probe::observe(&acquired);
        *prepared_rings = Some(acquired);
    }
    let rings = prepared_rings
        .as_ref()
        .expect("fast-or-better carrier present at the legacy owner");
    {
        let (atom_rows, bond_rows) = (rings.atom_row_count(), rings.bond_row_count());
        if atom_rows != record.topology.atoms.len() {
            return Err(SmilesStereoError::RingRows {
                dimension: "atom",
                actual: atom_rows,
                expected: record.topology.atoms.len(),
            });
        }
        if bond_rows != record.topology.bonds.len() {
            return Err(SmilesStereoError::RingRows {
                dimension: "bond",
                actual: bond_rows,
                expected: record.topology.bonds.len(),
            });
        }
    }
    // Remember only possible H-cleanup rows, not a cloned topology or an
    // additional complete valence assignment.
    let cleanup_candidates = record
        .topology
        .atoms
        .iter()
        .filter(|atom| {
            matches!(
                atom.chiral_tag(),
                cosmolkit_model::ChiralTag::TetrahedralCw
                    | cosmolkit_model::ChiralTag::TetrahedralCcw
            ) && atom.explicit_hydrogens() == 1
                && atom.formal_charge() == 0
                && !atom.is_aromatic()
        })
        .map(|atom| atom.id())
        .collect::<Vec<_>>();
    #[cfg(test)]
    ring_probe::record_legacy_entry();
    record.topology = assign_legacy_stereochemistry(record.topology, valence, rings)?;
    // RDKit✔️✔️:       atom->setNumExplicitHs(0);
    // RDKit✔️✔️:       atom->setNoImplicit(false);
    // RDKit✔️✔️:       atom->calcExplicitValence(false);
    // RDKit✔️✔️:       atom->calcImplicitValence(false);
    // Behavior: legacy cleanup can change H state after sanitization. Refresh
    // exactly those changed rows through the existing atom-valence owner.
    // Complexity: O(V) candidate scan, O(tagged atoms) temporary IDs and O(degree)
    // per changed row; no second whole-graph assignment or topology clone.
    for id in cleanup_candidates {
        if record.topology.atoms[id.index()].explicit_hydrogens() == 0 {
            let (explicit, implicit) = assign_valence_state_for_atom_from_parts(
                &record.topology.atoms,
                &record.topology.bonds,
                &record.topology.adjacency,
                id,
                false,
            )?;
            valence.explicit_valence[id.index()] = explicit;
            valence.implicit_hydrogens[id.index()] = implicit;
        }
    }
    // RDKit Chirality.cpp::assignStereochemistry (pinned legacy profile):
    // RDKit✔️✔️:   mol.setProp(common_properties::_StereochemDone, 1, true);
    // Behavior: topology-only core dispatch cannot write molecule properties.
    // Transport the successful computed marker here, so downstream consumers
    // preserve source hasProp() semantics without repeating stereo assignment.
    // Complexity: one property insertion; no topology copy or perception pass.
    record
        .properties
        .set_computed_prop("_StereochemDone", cosmolkit_model::PropertyValue::Int(1))?;
    Ok(record)
}

pub(crate) fn source_stereo_conformers(
    coordinates: &cosmolkit_model::CoordinateBlock,
) -> Result<
    (
        Option<std::borrow::Cow<'_, Conformer3D>>,
        Option<&Conformer3D>,
    ),
    cosmolkit_model::CoordinateValidationError,
> {
    // RDKit❗✔️:   const Conformer *conf = nullptr, *conf3d = nullptr;
    // RDKit❗✔️:   if (res && res->getNumConformers() > 0) {
    // RDKit❗✔️:     for (unsigned int confId = 0; confId < res->getNumConformers(); ++confId) {
    // RDKit❗✔️:       auto *testConf = &res->getConformer(confId);
    // RDKit❗✔️:       if (!testConf->is3D()) {
    // RDKit❗✔️:         if (conf == nullptr) {  // only take the first 2d conf
    // RDKit❗✔️:           conf = testConf;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         if (conf3d == nullptr) {  // only take the first 3d conf
    // RDKit❗✔️:           conf3d = testConf;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (conf != nullptr && conf3d != nullptr) {
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // Behavior: source insertion order is the existing MODEL occurrence order;
    // source is3D is independent of XYZ storage. A false-flag XYZ row remains
    // borrowed with all Z bits; genuine XY alone receives a transient lift.
    // Missing mixed order is a structural error, never an inferred preference.
    // Complexity: source-order scan stops when both flags are found, O(C).
    // XYZ borrows allocate nothing; the existing XY lift costs O(V) once.
    use cosmolkit_model::{CoordinateDimension, CoordinateValidationError};
    let (mut xy, mut xyz) = (0usize, 0usize);
    let (mut two_d, mut three_d) = (None, None);
    let order = coordinates.source_conformer_order.as_deref();
    if order.is_none()
        && !coordinates.conformers_2d.is_empty()
        && !coordinates.conformers_3d.is_empty()
    {
        return Err(CoordinateValidationError::MissingSourceConformerOrder);
    }
    let count = order.map_or(
        coordinates.conformers_2d.len() + coordinates.conformers_3d.len(),
        |order| order.len(),
    );
    for index in 0..count {
        let dimension = order.map_or_else(
            || {
                if coordinates.conformers_2d.is_empty() {
                    CoordinateDimension::ThreeD
                } else {
                    CoordinateDimension::TwoD
                }
            },
            |order| order[index],
        );
        match dimension {
            CoordinateDimension::TwoD => {
                let conformer = coordinates.conformers_2d.get(xy).ok_or(
                    CoordinateValidationError::SourceConformerOrder {
                        two_d: xy + 1,
                        three_d: xyz,
                        expected_two_d: coordinates.conformers_2d.len(),
                        expected_three_d: coordinates.conformers_3d.len(),
                    },
                )?;
                xy += 1;
                if two_d.is_none() {
                    two_d = Some(std::borrow::Cow::Owned(Conformer3D::new(
                        conformer.id(),
                        conformer
                            .coordinates()
                            .iter()
                            .map(|point| [point[0], point[1], 0.0])
                            .collect(),
                        false,
                    )));
                }
            }
            CoordinateDimension::ThreeD => {
                let conformer = coordinates.conformers_3d.get(xyz).ok_or(
                    CoordinateValidationError::SourceConformerOrder {
                        two_d: xy,
                        three_d: xyz + 1,
                        expected_two_d: coordinates.conformers_2d.len(),
                        expected_three_d: coordinates.conformers_3d.len(),
                    },
                )?;
                xyz += 1;
                if conformer.is_3d() {
                    if three_d.is_none() {
                        three_d = Some(conformer);
                    }
                } else if two_d.is_none() {
                    two_d = Some(std::borrow::Cow::Borrowed(conformer));
                }
            }
        }
        if two_d.is_some() && three_d.is_some() {
            break;
        }
    }
    Ok((two_d, three_d))
}

#[cfg(test)]
mod source_conformer_flag_regressions {
    use super::source_stereo_conformers;
    use cosmolkit_model::{
        Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, CoordinateValidationError,
    };

    #[test]
    fn false_flag_xyz_retains_z_bits_and_is_borrowed_in_source_order() {
        let coordinates = CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(9, vec![[1.0, 2.0, 3.0]], true),
                Conformer3D::new(2, vec![[4.0, 5.0, -0.0]], false),
                Conformer3D::new(1, vec![[6.0, 7.0, 8.0]], false),
            ],
            ..Default::default()
        };
        let before = coordinates.clone();
        let (two_d, three_d) = source_stereo_conformers(&coordinates).unwrap();
        let std::borrow::Cow::Borrowed(two_d) = two_d.unwrap() else {
            panic!("XYZ must be borrowed")
        };
        assert!(std::ptr::eq(two_d, &coordinates.conformers_3d[1]));
        assert_eq!(two_d.coordinates()[0][2].to_bits(), (-0.0_f64).to_bits());
        assert!(std::ptr::eq(
            three_d.unwrap(),
            &coordinates.conformers_3d[0]
        ));
        assert_eq!(coordinates, before);
    }

    #[test]
    fn mixed_storage_occurrence_order_selects_the_actual_first_false_flag() {
        for xyz_first in [false, true] {
            let order = if xyz_first {
                vec![
                    CoordinateDimension::ThreeD,
                    CoordinateDimension::TwoD,
                    CoordinateDimension::ThreeD,
                ]
            } else {
                vec![
                    CoordinateDimension::TwoD,
                    CoordinateDimension::ThreeD,
                    CoordinateDimension::ThreeD,
                ]
            };
            let coordinates = CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(7, vec![[1.0, 2.0]])],
                conformers_3d: vec![
                    Conformer3D::new(4, vec![[3.0, 4.0, 0.0001]], false),
                    Conformer3D::new(6, vec![[5.0, 6.0, 7.0]], true),
                ],
                source_conformer_order: Some(order),
                ..Default::default()
            };
            let before = coordinates.clone();
            let (two_d, three_d) = source_stereo_conformers(&coordinates).unwrap();
            let two_d = two_d.unwrap();
            assert_eq!(two_d.id(), if xyz_first { 4 } else { 7 });
            assert_eq!(
                two_d.coordinates()[0][2].to_bits(),
                if xyz_first {
                    0.0001_f64.to_bits()
                } else {
                    0.0_f64.to_bits()
                }
            );
            assert_eq!(three_d.unwrap().id(), 6);
            assert_eq!(coordinates, before);
        }
    }

    #[test]
    fn missing_mixed_order_is_a_typed_error_without_inference() {
        let coordinates = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(0, vec![])],
            conformers_3d: vec![Conformer3D::new(0, vec![], false)],
            ..Default::default()
        };
        let before = coordinates.clone();
        assert!(matches!(
            source_stereo_conformers(&coordinates),
            Err(CoordinateValidationError::MissingSourceConformerOrder)
        ));
        assert_eq!(coordinates, before);
    }
}
