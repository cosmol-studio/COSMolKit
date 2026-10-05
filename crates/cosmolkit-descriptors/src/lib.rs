//! Detached molecular mass, formula, and topology descriptor primitives.
//!
//! Descriptor boundaries borrow detached model values and explicit prepared
//! assignments. They cannot observe or mutate a live runtime molecule or its
//! derived-cache state.

mod chi;
mod connectivity;
pub mod counts;
mod crippen;
mod labute;
pub mod lipinski;
mod mqn;
mod qed;
pub use qed::qed;
pub(crate) mod patterns;
pub mod rings;
pub mod rotatable;
pub mod stereo;
mod tpsa;
mod vsa;

pub use chi::{
    CHI_0_N_VERSION, CHI_0_V_VERSION, CHI_1_N_VERSION, CHI_1_V_VERSION, CHI_2_N_VERSION,
    CHI_2_V_VERSION, CHI_3_N_VERSION, CHI_3_V_VERSION, CHI_4_N_VERSION, CHI_4_V_VERSION,
    CHI_N_N_VERSION, CHI_N_V_VERSION, chi_0_n, chi_0_n_with_state, chi_0_v, chi_0_v_with_state,
    chi_1_n, chi_1_n_with_state, chi_1_v, chi_1_v_with_state, chi_2_n, chi_2_n_with_state, chi_2_v,
    chi_2_v_with_state, chi_3_n, chi_3_n_with_state, chi_3_v, chi_3_v_with_state, chi_4_n,
    chi_4_n_with_state, chi_4_v, chi_4_v_with_state, chi_n_n, chi_n_n_with_state, chi_n_v,
    chi_n_v_with_state,
};
pub use connectivity::{
    HALL_KIER_ALPHA_VERSION, KAPPA_1_VERSION, KAPPA_2_VERSION, KAPPA_3_VERSION, PHI_VERSION, chi_0,
    chi_1, hall_kier_alpha, kappa_1, kappa_2, kappa_3, phi,
};
pub use crippen::{CrippenContributions, crippen_contributions};
pub use crippen::{CrippenParamRow, default_crippen_params};
pub use crippen::{CrippenTotals, crippen_totals};
pub use crippen::{crippen_clogp, crippen_mr};
pub use labute::{
    LabuteAsaContributions, LabuteContributions, labute_asa, labute_asa_contributions,
    labute_contributions,
};
pub use mqn::{MQN_VERSION, mqns};
pub use tpsa::{DescriptorComputedState, tpsa_contributions};
pub use vsa::{assign_contribs_to_bins, slogp_vsa, smr_vsa};

pub use counts::{num_atoms_prepared, num_heavy_atoms_prepared};
pub use lipinski::{
    NUM_AMIDE_BONDS_VERSION, NUM_HBA_VERSION, NUM_HBD_VERSION, NUM_HETEROATOMS_PATTERN,
    NUM_HETEROATOMS_VERSION, fraction_csp3_prepared, lipinski_hba_prepared, lipinski_hbd_prepared,
    num_amide_bonds_prepared, num_hba_prepared, num_hbd_prepared, num_heteroatoms_prepared,
};
pub use rings::{
    NUM_ALIPHATIC_CARBOCYCLES_VERSION, NUM_ALIPHATIC_HETEROCYCLES_VERSION,
    NUM_ALIPHATIC_RINGS_VERSION, NUM_AROMATIC_CARBOCYCLES_VERSION,
    NUM_AROMATIC_HETEROCYCLES_VERSION, NUM_AROMATIC_RINGS_VERSION, NUM_BRIDGEHEAD_ATOMS_VERSION,
    NUM_HETEROCYCLES_VERSION, NUM_RINGS_VERSION, NUM_SATURATED_CARBOCYCLES_VERSION,
    NUM_SATURATED_HETEROCYCLES_VERSION, NUM_SATURATED_RINGS_VERSION, NUM_SPIRO_ATOMS_VERSION,
    num_aliphatic_carbocycles_prepared, num_aliphatic_heterocycles_prepared,
    num_aliphatic_rings_prepared, num_aromatic_carbocycles_prepared,
    num_aromatic_heterocycles_prepared, num_aromatic_rings_prepared, num_bridgehead_atoms_prepared,
    num_bridgehead_atoms_with_ring_info, num_heterocycles_prepared, num_rings_prepared,
    num_saturated_carbocycles_prepared, num_saturated_heterocycles_prepared,
    num_saturated_rings_prepared, num_spiro_atoms_prepared, num_spiro_atoms_with_ring_info,
};
pub use rotatable::{
    NON_RING_AMIDES_PATTERN, NON_STRICT_ROTATABLE_PATTERN, NUM_ROTATABLE_BONDS_VERSION,
    RotatableBondsOptions, STRICT_LINKAGES_BASE_PATTERN, STRICT_ROTATABLE_PATTERN,
    SYMMETRIC_RINGS_PATTERN, TERMINAL_TRIPLE_BONDS_PATTERN, num_rotatable_bonds_prepared,
};
pub use stereo::{
    NUM_ATOM_STEREO_CENTERS_VERSION, NUM_UNSPECIFIED_ATOM_STEREO_CENTERS_VERSION,
    num_atom_stereo_centers_prepared, num_unspecified_atom_stereo_centers_prepared,
};

use std::{borrow::Cow, cmp::Ordering, collections::BTreeMap};

use cosmolkit_core::{
    LegacyStereoError, RingFindingError, RingInfo, ValenceAssignment, ValenceModel,
    assign_valence_with_options_for_topology, atomic_mass, most_common_isotope_mass,
    rdkit_element_symbol, total_hydrogen_count_from_validated,
};
use cosmolkit_model::{AtomId, CoordinateBlock, Element, MoleculeProperties, TopologyBlock};
use cosmolkit_search::{QueryMatchContextError, SmartsParseError, SubstructMatchError};

const RDKIT_ELECTRON_MASS: f64 = 0.00054857991;

/// Detached read input for the source Chi recomputation algorithms.
///
/// The pinned `hkDeltas`/`nVals`/Chi source reads topology and prepared valence
/// only. Coordinates, properties and ring state are not carrier inputs.
/// This value borrows existing rows without validation, preparation, cache
/// authority or identity inference. Each entrypoint retains its original
/// topology/valence structural validation and typed errors.
#[derive(Debug, Clone, Copy)]
pub struct ChiInput<'a> {
    topology: &'a TopologyBlock,
    valence: &'a ValenceAssignment,
}

impl<'a> ChiInput<'a> {
    /// Borrows the exact topology and prepared valence rows read by Chi.
    pub const fn new(topology: &'a TopologyBlock, valence: &'a ValenceAssignment) -> Self {
        Self { topology, valence }
    }

    /// Borrowed detached topology rows.
    pub const fn topology(&self) -> &'a TopologyBlock {
        self.topology
    }

    /// Borrowed prepared valence rows; never recomputed here.
    pub const fn valence(&self) -> &'a ValenceAssignment {
        self.valence
    }
}

/// Borrowed detached FINAL molecule state evaluated by descriptor functions.
///
/// This carrier is the packet-frozen detached input: references to the final
/// `TopologyBlock`, `CoordinateBlock`, `MoleculeProperties`, and the prepared
/// `ValenceAssignment`/`RingInfo` rows. No live `Molecule`, runtime cache
/// authority, or interior mutability enters it, and a carrier instance never
/// infers input identity from pointers. Callers that change chemistry build a
/// new input and pass a fresh [`DescriptorComputedState`] (added with the
/// cache-dependent owners) or explicitly clear it.
#[derive(Debug, Clone, Copy)]
pub struct DescriptorInput<'a> {
    topology: &'a TopologyBlock,
    coordinates: &'a CoordinateBlock,
    properties: &'a MoleculeProperties,
    valence: &'a ValenceAssignment,
    ring_info: &'a RingInfo,
}

/// RINGS-NARROW1-112 frozen detached ring-read boundary: shared test-only
/// fixture preparation and the frozen 18-case literal table. Fixture setup
/// happens OUTSIDE every counted descriptor invocation, through the existing
/// real owners only (detached parser with `remove_hydrogens=false`, the
/// real RemoveHs owner with `update_explicit_count=true, sanitize=true`
/// for the remove policy, `sanitize_topology` with ALL defaults for the
/// keep policy, the existing non-strict RdkitLike assignment owner for the
/// prepared-control input, and one explicit `symmetrized_sssr` with
/// default `RingSearchParams` for the frozen table). This setup makes NO
/// constructor-cache transport or source-stage acceptance claim.
#[cfg(test)]
mod ring_read_boundary_tests {
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, CoordinateBlock,
        Element, MoleculeProperties, TopologyBlock,
    };

    /// Frozen 18-case x 11-literal table (columns R01..R11 in unit order),
    /// fixed BEFORE implementation from RDKit 2026.03.1 and cross-checked
    /// with the pinned classifier predicates; never derived from the tested
    /// kernels. Both H policies use the SAME exact SMILES; keep/remove atom
    /// row counts are the two `usize` fields.
    const RING_READ_CASES: [(&str, usize, usize, [u32; 11]); 18] = [
        ("", 0, 0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        ("C", 1, 1, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        ("C1CCCCC1", 6, 6, [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1]),
        ("C1=CCCCC1", 6, 6, [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]),
        ("c1ccccc1", 6, 6, [1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0]),
        ("n1ccccc1", 6, 6, [1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0]),
        ("N1CCCCC1", 6, 6, [1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0]),
        (
            "C1CCC2(CC1)CCCC2",
            10,
            10,
            [2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 2],
        ),
        ("c1ccc2ccccc2c1", 10, 10, [2, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0]),
        (
            "C12C3C4C1C5C2C3C45",
            8,
            8,
            [6, 0, 0, 6, 6, 0, 0, 0, 6, 0, 6],
        ),
        ("[H]C1CCCCC1", 7, 6, [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1]),
        ("[2H]C1CCCCC1", 7, 7, [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1]),
        ("*1CCCC1", 5, 5, [1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0]),
        ("C1CCC2CCCCC2C1", 10, 10, [2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 2]),
        ("C1CC2CCC1C2", 7, 7, [2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 2]),
        ("c1ccc2CCCCc2c1", 10, 10, [2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0]),
        (
            "c1ccncc1C1CCCCC1",
            12,
            12,
            [2, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1],
        ),
        (
            "N12N3N4N1N5N2N3N45",
            8,
            8,
            [6, 6, 0, 6, 6, 0, 0, 6, 0, 6, 0],
        ),
    ];

    struct RingReadFixture {
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
        assignment: cosmolkit_core::ValenceAssignment,
        ring_info: cosmolkit_core::RingInfo,
    }

    /// Setup-only fixture through the EXISTING owners. The explicit
    /// non-strict RdkitLike assignment below also covers the empty-source
    /// RemoveHs branch that carries no final assignment — through the same
    /// real owner, never a fake default row.
    fn ring_read_fixture(smiles: &str, remove_hydrogens: bool) -> RingReadFixture {
        let params = cosmolkit_smiles::SmilesParseParams {
            remove_hydrogens: false,
            ..Default::default()
        };
        let record =
            cosmolkit_smiles::parse_smiles(smiles, &params).expect("parse frozen ring SMILES");
        let (topology, coordinates, properties) = if remove_hydrogens {
            let result = cosmolkit_core::remove_hydrogens_with_params(
                record.topology,
                record.coordinates,
                record.properties,
                &cosmolkit_core::RemoveHsParams {
                    update_explicit_count: true,
                    sanitize: true,
                    ..cosmolkit_core::RemoveHsParams::default()
                },
            )
            .expect("real RemoveHs owner on the original SMILES");
            (result.topology, result.coordinates, result.properties)
        } else {
            let result = cosmolkit_core::sanitize_topology(
                &record.topology,
                &cosmolkit_core::SanitizeParams::default(),
            )
            .expect("sanitize owner with ALL defaults");
            (result.topology, record.coordinates, record.properties)
        };
        let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
            &topology,
            cosmolkit_core::ValenceModel::RdkitLike,
            false,
        )
        .expect("non-strict RdkitLike assignment for the prepared-control input");
        let ring_info = cosmolkit_core::symmetrized_sssr(
            &topology,
            &cosmolkit_core::RingSearchParams::default(),
        )
        .expect("explicit symmetrized SSSR rows for the frozen table");
        RingReadFixture {
            topology,
            coordinates,
            properties,
            assignment,
            ring_info,
        }
    }

    /// Cube topology (the `C12C3C4C1C5C2C3C45` connectivity) with ALL
    /// eight atoms set to `element`: 8 atoms, 12 single bonds, every
    /// vertex degree three — `find_sssr` yields five rings and
    /// `symmetrized_sssr` yields six on this SAME topology.
    fn cube_topology(element: Element) -> TopologyBlock {
        let atoms: Vec<Atom> = (0..8)
            .map(|row| Atom::from_spec(AtomId::new(row), AtomSpec::new(element)))
            .collect();
        let edges = [
            (0usize, 1usize),
            (1, 2),
            (2, 3),
            (3, 0),
            (3, 4),
            (4, 5),
            (5, 0),
            (5, 6),
            (6, 1),
            (6, 7),
            (7, 2),
            (7, 4),
        ];
        let bonds: Vec<Bond> = edges
            .iter()
            .enumerate()
            .map(|(row, (left, right))| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(*left), AtomId::new(*right), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    /// Per-call prerequisite dimensions, checked immediately before EVERY
    /// counted call: paired ring-row counts, the two EXISTING public
    /// membership-size accessors against the actual table sizes, per-row
    /// atom/bond widths, and BOTH membership sums (atoms and bonds)
    /// against the stored ring rows. The fixture-level
    /// `is_symm_sssr`/`is_sssr_or_better` initialization predicates are
    /// asserted separately at their unchanged sites; these size checks
    /// are copied-source dimension correspondence, NOT a new
    /// initialization enforcement.
    fn assert_membership_dimensions(
        topology: &TopologyBlock,
        ring_info: &cosmolkit_core::RingInfo,
        context: &str,
    ) {
        assert_eq!(
            ring_info.atom_rings().len(),
            ring_info.bond_rings().len(),
            "paired ring-row counts ({context})"
        );
        assert_eq!(
            ring_info.atom_row_count(),
            topology.atoms.len(),
            "atom membership rows == atom rows ({context})"
        );
        assert_eq!(
            ring_info.bond_row_count(),
            topology.bonds.len(),
            "bond membership rows == bond rows ({context})"
        );
        let mut atom_membership = 0usize;
        for row in 0..topology.atoms.len() {
            atom_membership += ring_info.atom_members(AtomId::new(row)).len();
        }
        let atom_row_entries: usize = ring_info.atom_rings().iter().map(|ring| ring.len()).sum();
        assert_eq!(
            atom_membership, atom_row_entries,
            "atom-membership sum ({context})"
        );
        let mut bond_membership = 0usize;
        for row in 0..topology.bonds.len() {
            bond_membership += ring_info.bond_members(BondId::new(row)).len();
        }
        let bond_row_entries: usize = ring_info.bond_rings().iter().map(|ring| ring.len()).sum();
        assert_eq!(
            bond_membership, bond_row_entries,
            "bond-membership sum ({context})"
        );
        for (atom_row, bond_row) in ring_info
            .atom_rings()
            .iter()
            .zip(ring_info.bond_rings().iter())
        {
            assert_eq!(atom_row.len(), bond_row.len(), "row widths ({context})");
        }
    }

    /// Concrete per-call snapshot: the five full input values plus the
    /// CURRENT helper-counter value captured immediately before ONE
    /// counted call. Invoked at each real call site; this is a private
    /// concrete snapshot helper, not a generic runner or dispatch
    /// abstraction.
    struct CallCheckpoint {
        topology: TopologyBlock,
        ring_info: cosmolkit_core::RingInfo,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
        assignment: cosmolkit_core::ValenceAssignment,
        helper_baseline: usize,
    }

    fn capture_call_checkpoint(
        topology: &TopologyBlock,
        ring_info: &cosmolkit_core::RingInfo,
        coordinates: &CoordinateBlock,
        properties: &MoleculeProperties,
        assignment: &cosmolkit_core::ValenceAssignment,
    ) -> CallCheckpoint {
        CallCheckpoint {
            topology: topology.clone(),
            ring_info: ring_info.clone(),
            coordinates: coordinates.clone(),
            properties: properties.clone(),
            assignment: assignment.clone(),
            helper_baseline: crate::VALENCE_HELPER_ENTRIES.with(std::cell::Cell::get),
        }
    }

    /// Immediately after the SAME call: full-value equality of all five
    /// inputs with the captured snapshot (PartialEq on these concrete
    /// fixture types — NOT a generic NaN/signed-zero proof) and the
    /// counter equal to the captured, never-reset baseline. The counter
    /// observation proves only this helper's non-entry across the call;
    /// no-recomputation is established by the one-call delegate shape and
    /// the supplied 5/6 row-set artifacts.
    fn verify_call_checkpoint(
        checkpoint: &CallCheckpoint,
        topology: &TopologyBlock,
        ring_info: &cosmolkit_core::RingInfo,
        coordinates: &CoordinateBlock,
        properties: &MoleculeProperties,
        assignment: &cosmolkit_core::ValenceAssignment,
        context: &str,
    ) {
        assert_eq!(
            topology, &checkpoint.topology,
            "topology full-value equal ({context})"
        );
        assert_eq!(
            ring_info, &checkpoint.ring_info,
            "full rows full-value equal ({context})"
        );
        assert_eq!(
            coordinates, &checkpoint.coordinates,
            "coordinates full-value equal ({context})"
        );
        assert_eq!(
            properties, &checkpoint.properties,
            "properties full-value equal ({context})"
        );
        assert_eq!(
            assignment, &checkpoint.assignment,
            "assignment full-value equal ({context})"
        );
        assert_eq!(
            crate::VALENCE_HELPER_ENTRIES.with(std::cell::Cell::get),
            checkpoint.helper_baseline,
            "helper counter equals captured baseline ({context})"
        );
    }

    #[test]
    fn descriptor_ring_read_r01_product() {
        // R01 frozen product: 18 cases x 2 H policies x 2 routes
        // (narrow delegate + existing prepared) x 2 repeats = 144 calls.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                // Primary initialization prerequisite (unchanged site).
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                // Real H/deuterium row identity — never a bare count proof.
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide original snapshots AFTER
                // setup, retained as SUPPLEMENTARY persistent baselines;
                // the primary proof is the fresh per-call checkpoint.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_rings_with_ring_info(&fixture.ring_info)
                        .expect("narrow delegate");
                    assert_eq!(got, literals[0], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got = crate::num_rings_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[0], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                // Supplementary persistent baselines still hold.
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R01 product census");
    }

    #[test]
    fn descriptor_ring_read_r01_supplied_rows() {
        // R01 supplied-row control: 2 cube topologies x 2 supplied row sets
        // (find_sssr=5, symmetrized_sssr=6) x 2 routes = 8 calls. The SAME
        // topology is kept for both row sets, proving the borrowed supply
        // is honored rather than one ring finder being rerun.
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                // R01 literal = the SUPPLIED row count (5/6 for both
                // elements; the six-row literals are the frozen cube rows).
                let expected = expected_count as u32;
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_rings_with_ring_info(&ring_info).expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got = crate::num_rings_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                // Supplementary persistent baselines still hold.
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R01 supplied-rows census");
    }

    #[test]
    fn descriptor_ring_read_r02_product() {
        // R02 frozen product: 18 cases x 2 H policies x 2 routes
        // (narrow delegate + existing prepared) x 2 repeats = 144 calls.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                // Primary initialization prerequisite (unchanged site).
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide original snapshots AFTER
                // setup, retained as SUPPLEMENTARY persistent baselines;
                // the primary proof is the fresh per-call checkpoint.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_heterocycles_with_ring_info(
                        &fixture.topology,
                        &fixture.ring_info,
                    )
                    .expect("narrow delegate");
                    assert_eq!(got, literals[1], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got = crate::num_heterocycles_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[1], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                // Supplementary persistent baselines still hold.
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R02 product census");
    }

    #[test]
    fn descriptor_ring_read_r02_supplied_rows() {
        // R02 supplied-row control: 2 cubes x 2 row sets x 2 routes = 8
        // calls. R02 literals: carbon cube 0/0 (five/six rows); nitrogen
        // cube 5/6 — the element and row-set distinctions both bite.
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                let expected = match (element, five_rows) {
                    (Element::C, true) => 0,
                    (Element::C, false) => 0,
                    (Element::N, true) => 5,
                    (Element::N, false) => 6,
                    _ => unreachable!("two elements x two row sets"),
                };
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_heterocycles_with_ring_info(&topology, &ring_info)
                    .expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got = crate::num_heterocycles_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                // Supplementary persistent baselines still hold.
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R02 supplied-rows census");
    }

    #[test]
    fn descriptor_ring_read_r03_product() {
        // R03 frozen product: 18 cases x 2 H policies x 2 routes
        // (narrow delegate + existing prepared) x 2 repeats = 144 calls.
        // Corrected pattern: is_symm_sssr is the fixture prerequisite at
        // its unchanged site; per call — dimension correspondence, fresh
        // five-value checkpoint with current counter baseline, full-value
        // verify against the SAME call.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide snapshots, supplementary.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_aromatic_rings_with_ring_info(
                        &fixture.topology,
                        &fixture.ring_info,
                    )
                    .expect("narrow delegate");
                    assert_eq!(got, literals[2], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got = crate::num_aromatic_rings_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[2], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R03 product census");
    }

    #[test]
    fn descriptor_ring_read_r03_supplied_rows() {
        // R03 supplied-row control: 2 cubes x 2 row sets x 2 routes = 8
        // calls. R03 literals: 0 for ALL four cube combos — the synthetic
        // single-bond cubes never set aromatic flags, so every row
        // retracts at its first member; the control still proves the
        // supplied rows are honored through both routes.
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                let expected = 0u32;
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_aromatic_rings_with_ring_info(&topology, &ring_info)
                    .expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got = crate::num_aromatic_rings_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R03 supplied-rows census");
    }

    #[test]
    fn descriptor_ring_read_r04_product() {
        // R04 frozen product: 18 cases x 2 H policies x 2 routes
        // (narrow delegate + existing prepared) x 2 repeats = 144 calls.
        // Corrected pattern: is_symm_sssr fixture prerequisite at its
        // unchanged site; per call — dimension correspondence, fresh
        // five-value checkpoint with current counter baseline, full-value
        // verify against the SAME call.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide snapshots, supplementary.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_saturated_rings_with_ring_info(
                        &fixture.topology,
                        &fixture.ring_info,
                    )
                    .expect("narrow delegate");
                    assert_eq!(got, literals[3], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got = crate::num_saturated_rings_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[3], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R04 product census");
    }

    #[test]
    fn descriptor_ring_read_r04_supplied_rows() {
        // R04 supplied-row control: 2 cubes x 2 row sets x 2 routes = 8
        // calls. R04 literals: (C,5)=5, (C,6)=6, (N,5)=5, (N,6)=6 — every
        // cube bond is Single-ordered and non-aromatic, so the count
        // equals the SUPPLIED row count and both the row-set and element
        // distinctions are exercised.
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                let expected = expected_count as u32;
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_saturated_rings_with_ring_info(&topology, &ring_info)
                    .expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got = crate::num_saturated_rings_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R04 supplied-rows census");
    }

    #[test]
    fn descriptor_ring_read_r05_product() {
        // R05 frozen product: 18 cases x 2 H policies x 2 routes
        // (narrow delegate + existing prepared) x 2 repeats = 144 calls.
        // Corrected pattern: is_symm_sssr fixture prerequisite at its
        // unchanged site; per call — dimension correspondence, fresh
        // five-value checkpoint with current counter baseline, full-value
        // verify against the SAME call.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide snapshots, supplementary.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_aliphatic_rings_with_ring_info(
                        &fixture.topology,
                        &fixture.ring_info,
                    )
                    .expect("narrow delegate");
                    assert_eq!(got, literals[4], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got = crate::num_aliphatic_rings_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[4], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R05 product census");
    }

    #[test]
    fn descriptor_ring_read_r05_supplied_rows() {
        // R05 supplied-row control: 2 cubes x 2 row sets x 2 routes = 8
        // calls. R05 literals: (C,5)=5, (C,6)=6, (N,5)=5, (N,6)=6 — every
        // cube row has at least one non-aromatic member bond, so the
        // count equals the SUPPLIED row count.
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                let expected = expected_count as u32;
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_aliphatic_rings_with_ring_info(&topology, &ring_info)
                    .expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got = crate::num_aliphatic_rings_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R05 supplied-rows census");
    }

    #[test]
    fn descriptor_ring_read_r06_product() {
        // R06 frozen product: 18 cases x 2 H policies x 2 routes x 2
        // repeats = 144 calls. Corrected pattern: is_symm_sssr fixture
        // prerequisite at its unchanged site; per call — dimension
        // correspondence, fresh five-value checkpoint with current counter
        // baseline, full-value verify against the SAME call.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide snapshots, supplementary.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_aromatic_heterocycles_with_ring_info(
                        &fixture.topology,
                        &fixture.ring_info,
                    )
                    .expect("narrow delegate");
                    assert_eq!(got, literals[5], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got =
                        crate::num_aromatic_heterocycles_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[5], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R06 product census");
    }

    #[test]
    fn descriptor_ring_read_r06_supplied_rows() {
        // R06 supplied-row control: 2 cubes x 2 row sets x 2 routes = 8
        // calls. R06 literals: 0 for ALL four cube combos — the
        // non-aromatic single-bond cubes fail the all-aromatic conjunct.
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                let expected = 0u32;
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_aromatic_heterocycles_with_ring_info(&topology, &ring_info)
                    .expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got =
                    crate::num_aromatic_heterocycles_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R06 supplied-rows census");
    }

    #[test]
    fn descriptor_ring_read_r07_product() {
        // R07 frozen product: 18 cases x 2 H policies x 2 routes x 2
        // repeats = 144 calls. Corrected pattern throughout.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide snapshots, supplementary.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_aromatic_carbocycles_with_ring_info(
                        &fixture.topology,
                        &fixture.ring_info,
                    )
                    .expect("narrow delegate");
                    assert_eq!(got, literals[6], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got =
                        crate::num_aromatic_carbocycles_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[6], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R07 product census");
    }

    #[test]
    fn descriptor_ring_read_r07_supplied_rows() {
        // R07 supplied-row control: 2 cubes x 2 row sets x 2 routes = 8
        // calls. R07 literals: 0 for ALL four cube combos — the
        // non-aromatic single-bond cubes fail the aromatic arm.
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                let expected = 0u32;
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_aromatic_carbocycles_with_ring_info(&topology, &ring_info)
                    .expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got = crate::num_aromatic_carbocycles_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R07 supplied-rows census");
    }

    #[test]
    fn descriptor_ring_read_r08_product() {
        // R08 frozen product: 18 cases x 2 H policies x 2 routes x 2
        // repeats = 144 calls. Corrected pattern throughout.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide snapshots, supplementary.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_aliphatic_heterocycles_with_ring_info(
                        &fixture.topology,
                        &fixture.ring_info,
                    )
                    .expect("narrow delegate");
                    assert_eq!(got, literals[7], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got =
                        crate::num_aliphatic_heterocycles_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[7], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R08 product census");
    }

    #[test]
    fn descriptor_ring_read_r08_supplied_rows() {
        // R08 supplied-row control: 2 cubes x 2 row sets x 2 routes = 8
        // calls. R08 literals: (C,5)=0, (C,6)=0, (N,5)=5, (N,6)=6 — every
        // cube row is fully non-aromatic (hasAliph) and the nitrogen cube
        // rows additionally set hasHetero at every endpoint; the carbon
        // cube never sets hasHetero.
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                let expected = match (element, five_rows) {
                    (Element::C, true) => 0,
                    (Element::C, false) => 0,
                    (Element::N, true) => 5,
                    (Element::N, false) => 6,
                    _ => unreachable!("two elements x two row sets"),
                };
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_aliphatic_heterocycles_with_ring_info(&topology, &ring_info)
                    .expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got =
                    crate::num_aliphatic_heterocycles_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R08 supplied-rows census");
    }

    #[test]
    fn descriptor_ring_read_r09_product() {
        // R09 frozen product: 18 cases x 2 H policies x 2 routes x 2
        // repeats = 144 calls. Corrected pattern throughout.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide snapshots, supplementary.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_aliphatic_carbocycles_with_ring_info(
                        &fixture.topology,
                        &fixture.ring_info,
                    )
                    .expect("narrow delegate");
                    assert_eq!(got, literals[8], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got =
                        crate::num_aliphatic_carbocycles_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[8], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R09 product census");
    }

    #[test]
    fn descriptor_ring_read_r09_supplied_rows() {
        // R09 supplied-row control: 2 cubes x 2 row sets x 2 routes = 8
        // calls. R09 literals: (C,5)=5, (C,6)=6, (N,5)=0, (N,6)=0 — every
        // cube row is fully non-aromatic (hasAliph); the carbon cube has
        // no hetero endpoint so rows count; the nitrogen cube sets
        // hasHetero at every member and breaks out (disqualified).
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                let expected = match (element, five_rows) {
                    (Element::C, true) => 5,
                    (Element::C, false) => 6,
                    (Element::N, true) => 0,
                    (Element::N, false) => 0,
                    _ => unreachable!("two elements x two row sets"),
                };
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_aliphatic_carbocycles_with_ring_info(&topology, &ring_info)
                    .expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got =
                    crate::num_aliphatic_carbocycles_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R09 supplied-rows census");
    }

    #[test]
    fn descriptor_ring_read_r10_product() {
        // R10 frozen product: 18 cases x 2 H policies x 2 routes x 2
        // repeats = 144 calls. Corrected pattern throughout.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide snapshots, supplementary.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_saturated_heterocycles_with_ring_info(
                        &fixture.topology,
                        &fixture.ring_info,
                    )
                    .expect("narrow delegate");
                    assert_eq!(got, literals[9], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got =
                        crate::num_saturated_heterocycles_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[9], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R10 product census");
    }

    #[test]
    fn descriptor_ring_read_r10_supplied_rows() {
        // R10 supplied-row control: 2 cubes x 2 row sets x 2 routes = 8
        // calls. R10 literals: (C,5)=0, (C,6)=0, (N,5)=5, (N,6)=6 — every
        // cube row is all-Single non-aromatic (passes the conjunct); the
        // carbon cube never sets the hetero arm; the nitrogen cube sets
        // it at every member.
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                let expected = match (element, five_rows) {
                    (Element::C, true) => 0,
                    (Element::C, false) => 0,
                    (Element::N, true) => 5,
                    (Element::N, false) => 6,
                    _ => unreachable!("two elements x two row sets"),
                };
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_saturated_heterocycles_with_ring_info(&topology, &ring_info)
                    .expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got =
                    crate::num_saturated_heterocycles_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R10 supplied-rows census");
    }

    #[test]
    fn descriptor_ring_read_r11_product() {
        // R11 frozen product: 18 cases x 2 H policies x 2 routes x 2
        // repeats = 144 calls. Corrected pattern throughout.
        let mut calls = 0usize;
        for (smiles, rows_keep, rows_remove, literals) in RING_READ_CASES.iter() {
            for remove_hydrogens in [false, true] {
                let fixture = ring_read_fixture(smiles, remove_hydrogens);
                let expected_rows = if remove_hydrogens {
                    *rows_remove
                } else {
                    *rows_keep
                };
                assert_eq!(
                    fixture.topology.atoms.len(),
                    expected_rows,
                    "atom-row census ({smiles}, remove={remove_hydrogens})"
                );
                assert!(
                    fixture.ring_info.is_symm_sssr(),
                    "initialized symmetrized rows ({smiles}, remove={remove_hydrogens})"
                );
                let hydrogen_rows: Vec<_> = fixture
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.element() == Element::H)
                    .collect();
                match (*smiles, remove_hydrogens) {
                    ("[H]C1CCCCC1", false) => {
                        assert_eq!(hydrogen_rows.len(), 1, "kept neighbor H row");
                        assert_eq!(hydrogen_rows[0].isotope(), None, "protium identity");
                    }
                    ("[H]C1CCCCC1", true) => {
                        assert!(hydrogen_rows.is_empty(), "neighbor H removed by real owner");
                    }
                    ("[2H]C1CCCCC1", _) => {
                        assert_eq!(hydrogen_rows.len(), 1, "deuterium row retained");
                        assert_eq!(hydrogen_rows[0].isotope(), Some(2), "isotope Some(2)");
                    }
                    _ => {}
                }
                // Never-refreshed fixture-wide snapshots, supplementary.
                let topology_snapshot = fixture.topology.clone();
                let ring_snapshot = fixture.ring_info.clone();
                let coordinates_snapshot = fixture.coordinates.clone();
                let properties_snapshot = fixture.properties.clone();
                let assignment_snapshot = fixture.assignment.clone();
                for repeat in 0..2 {
                    let context = format!("{smiles}, remove={remove_hydrogens}, narrow, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let got = crate::num_saturated_carbocycles_with_ring_info(
                        &fixture.topology,
                        &fixture.ring_info,
                    )
                    .expect("narrow delegate");
                    assert_eq!(got, literals[10], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;

                    let context =
                        format!("{smiles}, remove={remove_hydrogens}, prepared, {repeat}");
                    assert_membership_dimensions(&fixture.topology, &fixture.ring_info, &context);
                    let checkpoint = capture_call_checkpoint(
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                    );
                    let input = crate::DescriptorInput::new(
                        &fixture.topology,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &fixture.ring_info,
                    );
                    let got =
                        crate::num_saturated_carbocycles_prepared(&input).expect("prepared route");
                    assert_eq!(got, literals[10], "literal ({context})");
                    verify_call_checkpoint(
                        &checkpoint,
                        &fixture.topology,
                        &fixture.ring_info,
                        &fixture.coordinates,
                        &fixture.properties,
                        &fixture.assignment,
                        &context,
                    );
                    calls += 1;
                }
                assert!(
                    fixture.topology == topology_snapshot
                        && fixture.ring_info == ring_snapshot
                        && fixture.coordinates == coordinates_snapshot
                        && fixture.properties == properties_snapshot
                        && fixture.assignment == assignment_snapshot,
                    "supplementary fixture-wide baselines ({smiles}, remove={remove_hydrogens})"
                );
            }
        }
        assert_eq!(calls, 144, "exact R11 product census");
    }

    #[test]
    fn descriptor_ring_read_r11_supplied_rows() {
        // R11 supplied-row control: 2 cubes x 2 row sets x 2 routes = 8
        // calls. R11 literals: (C,5)=5, (C,6)=6, (N,5)=0, (N,6)=0 — every
        // cube row is all-Single non-aromatic; the carbon cube has
        // all-carbon endpoints so rows count; the nitrogen cube fails
        // the endpoint arm at every member.
        let mut calls = 0usize;
        for element in [Element::C, Element::N] {
            let topology = cube_topology(element);
            assert!(
                topology.atoms.iter().all(|atom| atom.element() == element),
                "cube element prerequisite"
            );
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let assignment = cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )
            .expect("cube assignment");
            for five_rows in [true, false] {
                let ring_info = if five_rows {
                    cosmolkit_core::find_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("find_sssr five rows")
                } else {
                    cosmolkit_core::symmetrized_sssr(
                        &topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .expect("symmetrized six rows")
                };
                let expected_count = if five_rows { 5 } else { 6 };
                assert_eq!(
                    ring_info.atom_rings().len(),
                    expected_count,
                    "paired atom ring count ({element:?}, five_rows={five_rows})"
                );
                assert_eq!(
                    ring_info.bond_rings().len(),
                    expected_count,
                    "paired bond ring count ({element:?}, five_rows={five_rows})"
                );
                assert!(
                    ring_info.is_sssr_or_better(),
                    "initialized SSSR-or-better rows ({element:?}, five_rows={five_rows})"
                );
                let expected = match (element, five_rows) {
                    (Element::C, true) => 5,
                    (Element::C, false) => 6,
                    (Element::N, true) => 0,
                    (Element::N, false) => 0,
                    _ => unreachable!("two elements x two row sets"),
                };
                // Never-refreshed row-set snapshots, supplementary.
                let topology_snapshot = topology.clone();
                let ring_snapshot = ring_info.clone();
                let coordinates_snapshot = coordinates.clone();
                let properties_snapshot = properties.clone();
                let assignment_snapshot = assignment.clone();
                let context = format!("{element:?}, five_rows={five_rows}, narrow");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let got = crate::num_saturated_carbocycles_with_ring_info(&topology, &ring_info)
                    .expect("narrow delegate");
                assert_eq!(got, expected, "narrow ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                calls += 1;

                let context = format!("{element:?}, five_rows={five_rows}, prepared");
                assert_membership_dimensions(&topology, &ring_info, &context);
                let checkpoint = capture_call_checkpoint(
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                );
                let input = crate::DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                let got =
                    crate::num_saturated_carbocycles_prepared(&input).expect("prepared route");
                assert_eq!(got, expected, "prepared ({context})");
                verify_call_checkpoint(
                    &checkpoint,
                    &topology,
                    &ring_info,
                    &coordinates,
                    &properties,
                    &assignment,
                    &context,
                );
                assert!(
                    topology == topology_snapshot
                        && ring_info == ring_snapshot
                        && coordinates == coordinates_snapshot
                        && properties == properties_snapshot
                        && assignment == assignment_snapshot,
                    "supplementary row-set baselines ({element:?}, five_rows={five_rows})"
                );
                calls += 1;
            }
        }
        assert_eq!(calls, 8, "exact R11 supplied-rows census");
    }
}

impl<'a> DescriptorInput<'a> {
    /// Borrows the five prepared FINAL inputs.
    ///
    /// The caller owns correspondence between the blocks and the prepared
    /// rows (this constructor performs no shape validation; each descriptor
    /// function validates exactly the rows it reads, mirroring the source's
    /// per-function TEST_ASSERT/PRECONDITION sites recorded in the D01
    /// audit rather than a shared hidden precheck).
    pub const fn new(
        topology: &'a TopologyBlock,
        coordinates: &'a CoordinateBlock,
        properties: &'a MoleculeProperties,
        valence: &'a ValenceAssignment,
        ring_info: &'a RingInfo,
    ) -> Self {
        Self {
            topology,
            coordinates,
            properties,
            valence,
            ring_info,
        }
    }

    /// Final topology rows.
    pub const fn topology(&self) -> &'a TopologyBlock {
        self.topology
    }

    /// Final detached coordinate block.
    pub const fn coordinates(&self) -> &'a CoordinateBlock {
        self.coordinates
    }

    /// Final detached molecule-level properties.
    pub const fn properties(&self) -> &'a MoleculeProperties {
        self.properties
    }

    /// Prepared final valence rows (borrowed, never recomputed here).
    pub const fn valence(&self) -> &'a ValenceAssignment {
        self.valence
    }

    /// Prepared final ring rows.
    pub const fn ring_info(&self) -> &'a RingInfo {
        self.ring_info
    }
}

/// Source-cause category for fixed-pattern search work inside descriptor
/// functions.
///
/// RDKit compiles each fixed SMARTS once (Lipinski.cpp `ss_matcher` +
/// `pattern_flyweight`, lines 29-74) and aborts a failed compile with
/// `POSTCONDITION(m_matcher, "no matcher")`; matching itself has no error
/// path. COSMolKit models those abort sites as typed causes so a malformed
/// pattern or failing match propagates structurally instead of panicking.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum DescriptorSearchCause {
    /// Fixed SMARTS pattern failed to compile (source POSTCONDITION site).
    #[error("SMARTS compile failed: {0}")]
    Compile(#[from] SmartsParseError),
    /// Substructure match evaluation failed.
    #[error("substructure match failed: {0}")]
    Match(#[from] SubstructMatchError),
    /// Prepared shared query-context validation failed.
    #[error("prepared query context invalid: {0}")]
    Context(#[from] QueryMatchContextError),
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum DescriptorError {
    #[error("descriptor `{function}`: topology edit failed: {source}")]
    TopologyEdit {
        function: &'static str,
        #[source]
        source: cosmolkit_model::TopologyEditError,
    },
    #[error("descriptor `qed`: default RemoveHs did not supply final {field}")]
    MissingFinalHydrogenState { field: &'static str },

    /// A sole path-owner result violates the source TEST_ASSERT invariant.
    #[error(
        "descriptor `{function}`: connectivity path has {actual_rows:?} atom rows; expected {expected_rows}"
    )]
    InvalidConnectivityPath {
        function: &'static str,
        expected_rows: usize,
        actual_rows: Option<usize>,
    },
    /// A detached descriptor input violates local topology invariants.
    #[error("descriptor `{function}`: invalid topology: {source}")]
    InvalidTopology {
        function: &'static str,
        #[source]
        source: cosmolkit_model::TopologyValidationError,
    },
    /// The existing core path owner failed; its cause remains structural.
    #[error("descriptor `{function}`: path source error: {source}")]
    Path {
        function: &'static str,
        #[source]
        source: cosmolkit_core::PathError,
    },
    /// The optional Hall–Kier output sink violates the source PRECONDITION.
    #[error(
        "descriptor `hall_kier_alpha`: contribution sink has {actual} rows; requires at least {minimum}"
    )]
    InvalidHallKierContributionRows { actual: usize, minimum: usize },
    #[error(
        "descriptor `{function}`: cached TPSA contributions exist but the scalar total is missing (include_sulfur_phosphorus={include_sulfur_phosphorus})"
    )]
    MissingComputedScalar {
        function: &'static str,
        include_sulfur_phosphorus: bool,
    },
    /// Cached Labute rows exist but the `_labuteAtomHContrib` read fails.
    #[error(
        "descriptor `{function}`: cached Labute contributions exist but the hydrogen term is missing"
    )]
    MissingLabuteHydrogens { function: &'static str },
    /// Cached Labute rows/H exist but the `_labuteASA` read fails.
    #[error(
        "descriptor `{function}`: cached Labute contributions exist but the ASA total is missing"
    )]
    MissingLabuteAsa { function: &'static str },
    /// Bin-assignment array sizes violate the source PRECONDITIONs.
    #[error(
        "descriptor binning: contribs has {contribs_len} rows, bin_prop has {bin_prop_len}, res needs at least {bins_len}+1"
    )]
    MismatchedBinArrays {
        contribs_len: usize,
        bin_prop_len: usize,
        bins_len: usize,
    },
    /// A Crippen parameter numeric cell failed lexical parsing.
    #[error("descriptor `crippen`: parameter {field} cell {cell:?} is not a valid number")]
    CrippenParamNumeric { field: &'static str, cell: String },
    /// Cached Crippen LogP rows exist but the `_crippenMRContribs` read
    /// fails (source reads BOTH row vectors before length checks).
    #[error(
        "descriptor `{function}`: cached Crippen LogP contributions exist but the MR contributions are missing"
    )]
    MissingCrippenMrContributions { function: &'static str },
    /// Cached Crippen LogP scalar exists but the `_crippenMR` read fails
    /// (source reads BOTH scalars on the warm arm).
    #[error(
        "descriptor `{function}`: cached Crippen LogP total exists but the MR total is missing"
    )]
    MissingCrippenMr { function: &'static str },
    /// An optional Crippen output sink violates the source PRECONDITION
    /// size (checked BEFORE any cache access).
    #[error("descriptor `crippen`: optional {field} sink has {actual} rows; expected {expected}")]
    InvalidCrippenOptionalRows {
        field: &'static str,
        actual: usize,
        expected: usize,
    },
    /// The stored default Crippen pattern for a table row is unexpectedly
    /// None (all 110 are source-proven Some; the Option store retains no
    /// parser cause to preserve).
    #[error("descriptor `crippen`: default pattern for table row {row} is unexpectedly missing")]
    MissingCrippenDefaultPattern { row: u32 },
    /// The hydrogenated working-copy construction failed inside the
    /// modeled branch; the real core HydrogenError is preserved
    /// structurally (borrowed `Error::source`), never flattened.
    #[error("descriptor `{function}`: hydrogen transform failed: {source}")]
    Hydrogens {
        function: &'static str,
        #[source]
        source: cosmolkit_core::HydrogenError,
    },
    #[error("descriptor `{function}`: {field} accumulation exceeds its source integer domain")]
    CountOverflow {
        function: &'static str,
        field: &'static str,
    },
    #[error("descriptor `{function}`: {source}")]
    Valence {
        function: &'static str,
        source: cosmolkit_core::ValenceError,
    },
    #[error("descriptor `{function}`: {field} has {actual} rows; expected {expected}")]
    InvalidValenceRows {
        function: &'static str,
        field: &'static str,
        actual: usize,
        expected: usize,
    },
    #[error("descriptor `{function}` cannot evaluate this topology: {detail}")]
    Unsupported {
        function: &'static str,
        detail: String,
    },
    #[error("descriptor `{function}`: ring source error: {source}")]
    Ring {
        function: &'static str,
        source: RingFindingError,
    },
    #[error("descriptor `{function}`: search source error: {source}")]
    Search {
        function: &'static str,
        source: DescriptorSearchCause,
    },
    #[error("descriptor `{function}`: stereo source error: {source}")]
    Stereo {
        function: &'static str,
        source: LegacyStereoError,
    },
}

pub type DescriptorResult<T> = Result<T, DescriptorError>;

pub(crate) fn prepared_valence<'a>(
    topology: &TopologyBlock,
    assignment: Option<&'a cosmolkit_core::ValenceAssignment>,
    function: &'static str,
) -> DescriptorResult<Cow<'a, cosmolkit_core::ValenceAssignment>> {
    let value = match assignment {
        Some(value) => Cow::Borrowed(value),
        None => Cow::Owned(valence(topology, function)?),
    };
    for (field, actual) in [
        ("explicit_valence", value.explicit_valence.len()),
        ("implicit_hydrogens", value.implicit_hydrogens.len()),
    ] {
        if actual != topology.atoms.len() {
            return Err(DescriptorError::InvalidValenceRows {
                function,
                field,
                actual,
                expected: topology.atoms.len(),
            });
        }
    }
    Ok(value)
}

pub(crate) fn validate_topology(
    topology: &TopologyBlock,
    function: &'static str,
) -> DescriptorResult<()> {
    topology
        .validate()
        .map_err(|error| DescriptorError::Unsupported {
            function,
            detail: error.to_string(),
        })
}

pub(crate) fn valence(
    topology: &TopologyBlock,
    function: &'static str,
) -> DescriptorResult<cosmolkit_core::ValenceAssignment> {
    #[cfg(test)]
    VALENCE_HELPER_ENTRIES.with(|count| count.set(count.get() + 1));
    assign_valence_with_options_for_topology(topology, ValenceModel::RdkitLike, false).map_err(
        |error| DescriptorError::Unsupported {
            function,
            detail: error.to_string(),
        },
    )
}

#[cfg(test)]
thread_local! {
    /// Per-thread count of entries into the cold `valence` helper.
    /// Test-only evidence that ring-only cold descriptors such as
    /// `num_rings` perform ZERO valence-helper entries; not production
    /// instrumentation and not a claim about any other preparation cost.
    pub(crate) static VALENCE_HELPER_ENTRIES: std::cell::Cell<usize> =
        const { std::cell::Cell::new(0) };
}

fn descriptor_atomic_mass(
    element: Element,
    isotope: Option<u16>,
    function: &'static str,
) -> DescriptorResult<f64> {
    atomic_mass(element, isotope).map_err(|error| DescriptorError::Unsupported {
        function,
        detail: error.to_string(),
    })
}

/// Counts nitrogen and oxygen atoms using RDKit's direct Lipinski definition.
///
/// Cold convenience form; delegates to the one kernel in [`lipinski`].
pub fn lipinski_hba(topology: &TopologyBlock) -> DescriptorResult<u32> {
    lipinski::lipinski_hba_kernel(topology)
}

/// Sums hydrogens on nitrogen and oxygen using RDKit's direct Lipinski definition.
///
/// Cold convenience form: computes ONE cold canonical assignment, then the
/// one kernel in [`lipinski`].
pub fn lipinski_hbd(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let assignment = valence(topology, "lipinski_hbd")?;
    lipinski::lipinski_hbd_kernel(topology, &assignment)
}

/// Counts explicit atoms whose atomic number is greater than one.
///
/// Cold convenience form: delegates to the single kernel in [`counts`].
pub fn num_heavy_atoms(topology: &TopologyBlock) -> DescriptorResult<u32> {
    counts::num_heavy_atoms_kernel(topology)
}
/// Counts explicit atoms plus attached implicit/explicit atom-state hydrogens.
///
/// Cold convenience form: computes ONE cold canonical assignment, then the
/// one kernel in [`counts`].
pub fn num_atoms(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let assignment = valence(topology, "num_atoms")?;
    counts::num_atoms_kernel(topology, &assignment)
}

/// Mandatory-borrow total atom count (DQ01): explicit rows plus attached
/// implicit/explicit-property hydrogens from a CALLER-SUPPLIED validated
/// assignment.
///
/// This is the single delegation surface for the public root query
/// `Molecule::total_atom_count`. It narrows the cross-crate input to the
/// two borrowed values the source closure reads; it performs no loop of
/// its own, builds no `DescriptorInput`, never recomputes valence (the
/// kernel's `prepared_valence` takes the `Cow::Borrowed` arm), and adds
/// no clone, ring lookup or output buffering. All behavior, verbatim
/// source anchors (MolDescriptors.cpp:33-37, ROMol.cpp:176-186,
/// `includeNeighbors=false` explicit-H non-double-counting) and typed
/// error transport remain those of [`counts::num_atoms_kernel`].
pub fn total_atom_count_with_valence(
    topology: &TopologyBlock,
    assignment: &cosmolkit_core::ValenceAssignment,
) -> DescriptorResult<u32> {
    counts::num_atoms_kernel(topology, assignment)
}

/// Fraction of carbon atoms whose total degree is four.
///
/// Cold convenience form: computes ONE cold canonical assignment, then the
/// one kernel in [`lipinski`].
pub fn fraction_csp3(topology: &TopologyBlock) -> DescriptorResult<f64> {
    let assignment = valence(topology, "fraction_csp3")?;
    lipinski::fraction_csp3_kernel(topology, &assignment)
}

/// Mandatory-borrow Lipinski donor-hydrogen sum (DQ02): the N/O
/// `getTotalNumHs(true)` hydrogen sum from a CALLER-SUPPLIED validated
/// assignment.
///
/// Single delegation surface for the public root query
/// `Molecule::lipinski_hbd`. No loop of its own, no `DescriptorInput`, no
/// valence recompute (the kernel borrows the supplied assignment), no
/// clone or output buffering. All behavior, verbatim Lipinski.cpp anchors
/// and typed errors remain those of [`lipinski::lipinski_hbd_kernel`];
/// this is NOT the general SMARTS-based NumHBD.
pub fn lipinski_hbd_with_valence(
    topology: &TopologyBlock,
    assignment: &cosmolkit_core::ValenceAssignment,
) -> DescriptorResult<u32> {
    lipinski::lipinski_hbd_kernel(topology, assignment)
}

/// Mandatory-borrow fraction-of-CSP3 carbons (DQ03): carbon rows whose
/// SOURCE total degree (explicit neighbors + implicit/explicit-property
/// hydrogens, `includeNeighbors=false`) is exactly four, from a
/// CALLER-SUPPLIED validated assignment.
///
/// Single delegation surface for the public root query
/// `Molecule::fraction_csp3`. No loop of its own, no `DescriptorInput`, no
/// valence recompute, no clone or output buffering; the criterion is
/// source total degree, never a hybridization flag, and the zero-carbon
/// fallback is the kernel's literal 0.0. All behavior, verbatim
/// Lipinski.cpp / Atom.cpp getTotalDegree anchors and typed errors remain
/// those of [`lipinski::fraction_csp3_kernel`].
pub fn fraction_csp3_with_valence(
    topology: &TopologyBlock,
    assignment: &cosmolkit_core::ValenceAssignment,
) -> DescriptorResult<f64> {
    lipinski::fraction_csp3_kernel(topology, assignment)
}

/// Cold SSSR ring rows for SMARTS-count descriptor inputs.
///
/// One `find_sssr` per cold call with the packet-frozen search parameters;
/// typed [`DescriptorError::Ring`] on failure, never a default empty ring
/// set (the prepared query context requires initialized rings).
pub(crate) fn ring_info(
    topology: &TopologyBlock,
    function: &'static str,
) -> DescriptorResult<RingInfo> {
    cosmolkit_core::find_sssr(
        topology,
        &cosmolkit_core::RingSearchParams {
            include_dative_bonds: false,
            include_hydrogen_bonds: false,
        },
    )
    .map_err(|source| DescriptorError::Ring { function, source })
}

/// General SMARTS-based hydrogen-bond donor count.
///
/// Cold convenience form: computes ONE cold canonical assignment and ONE
/// cold SSSR ring set, then delegates to [`lipinski::num_hbd_prepared`]
/// (retained fixed pattern, default match parameters).
pub fn num_hbd(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let assignment = valence(topology, "num_hbd")?;
    // HBD-PUBLIC: the cold form still performs BOTH of its historical cold
    // preparations (canonical assignment AND the cold SSSR ring set — the
    // latter retained explicitly to preserve the old error surface, since
    // a failing cold ring find failed the old call too), then routes the
    // count through the ONE narrow HBD owner with the supplied
    // topology/assignment; the pattern itself reads no ring predicate.
    let _rings = ring_info(topology, "num_hbd")?;
    lipinski::num_hbd_with_valence(topology, &assignment)
}

/// General SMARTS-based hydrogen-bond donor count over an existing
/// prepared valence assignment (HBD-PUBLIC narrow domain form).
///
/// ONE qualified delegation to the narrow lipinski HBD owner
/// ([`lipinski::num_hbd_with_valence`]): the fixed pattern reads
/// hydrogen-count and valence rows but NO ring predicate, so this form
/// borrows the supplied FINAL topology and prepared valence rows through
/// the ONE narrow borrowed-valence context and matcher path. No ring
/// state is read, gated on, found or fabricated; the assignment is never
/// recomputed, installed or cloned.
pub fn num_hbd_with_valence(
    topology: &TopologyBlock,
    valence: &cosmolkit_core::ValenceAssignment,
) -> DescriptorResult<u32> {
    lipinski::num_hbd_with_valence(topology, valence)
}

/// General SMARTS-based hydrogen-bond acceptor count.
///
/// Cold convenience form: computes ONE cold canonical assignment and ONE
/// cold SSSR ring set, then delegates to [`lipinski::num_hba_prepared`]
/// (retained fixed recursive pattern, default match parameters).
pub fn num_hba(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let assignment = valence(topology, "num_hba")?;
    let rings = ring_info(topology, "num_hba")?;
    let coordinates = CoordinateBlock::default();
    let properties = MoleculeProperties::default();
    let input = DescriptorInput::new(topology, &coordinates, &properties, &assignment, &rings);
    lipinski::num_hba_prepared(&input)
}

/// General SMARTS-based heteroatom count.
///
/// Topology-only SMARTS heteroatom count (RDKit `CalcNumHeteroatoms`).
///
/// D-A resolution: the `[!#6;!#1]` pattern reads ONLY atomic-number rows,
/// so this entry delegates to [`lipinski::num_heteroatoms_topology`] via
/// the narrow topology query context — NO valence assignment, ring find,
/// fabricated chemistry rows or direct element counting; default match
/// parameters (uniquify, maxMatches=1000) retained.
pub fn num_heteroatoms(topology: &TopologyBlock) -> DescriptorResult<u32> {
    lipinski::num_heteroatoms_topology(topology)
}

/// HETERO-PUBLIC domain regressions: the frozen 160-call input-state
/// product, the 12-call maxMatches boundary product, and the raw
/// pentavalent-carbon discriminator over the topology-only owner.
#[cfg(test)]
mod descriptor_heteroatoms_domain_tests {
    use super::*;

    /// Frozen 20-case literal table (SMILES, expected count); identical
    /// under BOTH remove-H policies and raw/sanitized construction.
    const CASES: [(&str, u32); 20] = [
        ("", 0),
        ("C", 0),
        ("CCO", 1),
        ("[NH4+]", 1),
        ("[O-]", 1),
        ("[H][H]", 0),
        ("[2H]O[2H]", 1),
        ("[13CH4]", 0),
        ("N", 1),
        ("O", 1),
        ("C=O", 1),
        ("O=C(N)N", 3),
        ("n1ccccc1", 1),
        ("[nH]1cccc1", 1),
        ("C1CCCCC1", 0),
        ("CC(N)C(=O)O", 3),
        ("[H]N([H])[H]", 1),
        ("C[C+](C)C", 0),
        ("*", 1),
        ("CC#N", 1),
    ];

    /// Frozen literal CONSTRUCTOR prerequisites (input identities from the
    /// frozen SMILES strings and remove-H policy, from the DQ literal
    /// table — never derived from num_heteroatoms output). Ordered atomic
    /// numbers per case, isotopes, and explicit-H row specifications.
    /// (smiles, ordered atomic numbers keep, ordered isotopes keep,
    /// explicit-H rows keep)
    const PREREQS: [(&str, &[u8], &[Option<u16>], &[(usize, u8)]); 20] = [
        ("", &[], &[], &[]),
        ("C", &[6], &[None], &[]),
        ("CCO", &[6, 6, 8], &[None, None, None], &[]),
        ("[NH4+]", &[7], &[None], &[(0, 4)]),
        ("[O-]", &[8], &[None], &[]),
        ("[H][H]", &[1, 1], &[None, None], &[]),
        ("[2H]O[2H]", &[1, 8, 1], &[Some(2), None, Some(2)], &[]),
        ("[13CH4]", &[6], &[Some(13)], &[(0, 4)]),
        ("N", &[7], &[None], &[]),
        ("O", &[8], &[None], &[]),
        ("C=O", &[6, 8], &[None, None], &[]),
        ("O=C(N)N", &[8, 6, 7, 7], &[None, None, None, None], &[]),
        ("n1ccccc1", &[7, 6, 6, 6, 6, 6], &[None; 6], &[]),
        ("[nH]1cccc1", &[7, 6, 6, 6, 6], &[None; 5], &[(0, 1)]),
        ("C1CCCCC1", &[6, 6, 6, 6, 6, 6], &[None; 6], &[]),
        ("CC(N)C(=O)O", &[6, 6, 7, 6, 8, 8], &[None; 6], &[]),
        // Ammonia: H rows survive ONLY remove-H=false (1,7,1,1); removed
        // policy yields just [7] with N explicit-H 3 AFTER removal.
        ("[H]N([H])[H]", &[1, 7, 1, 1], &[None; 4], &[(1, 0)]),
        ("C[C+](C)C", &[6, 6, 6, 6], &[None; 4], &[]),
        ("*", &[0], &[None], &[]),
        ("CC#N", &[6, 6, 7], &[None, None, None], &[]),
    ];

    /// Assert the frozen input identities on a CONSTRUCTED topology row
    /// list (atomic numbers in order, isotope options, explicit-H rows).
    /// Input prerequisite check only — no descriptor call, no expected
    /// count synthesis.
    fn assert_input_prerequisites(
        label: &str,
        smiles: &str,
        topology: &TopologyBlock,
        remove_hydrogens: bool,
    ) {
        let index = PREREQS
            .iter()
            .position(|(candidate, _, _, _)| *candidate == smiles)
            .unwrap_or_else(|| panic!("{label}: unknown prerequisite case {smiles:?}"));
        let (_, atomic_numbers, isotopes, explicit_h) = PREREQS[index];
        // The remove-H=true policy removes plain hydrogen ATOM rows
        // attached to heavy atoms (ammonia H rows) while deuterium rows
        // survive BOTH policies; an all-hydrogen molecule ([H][H]) keeps
        // its rows under both policies (DQ frozen rows table: [H][H] 2/2,
        // [2H]O[2H] 3/3, ammonia 4/1).
        let all_hydrogen = atomic_numbers.iter().all(|&z| z == 1);
        let expected_numbers: Vec<u8> = if remove_hydrogens {
            atomic_numbers
                .iter()
                .zip(isotopes.iter())
                .filter(|(z, isotope)| **z != 1 || isotope.is_some() || all_hydrogen)
                .map(|(&z, _)| z)
                .collect()
        } else {
            atomic_numbers.to_vec()
        };
        assert_eq!(
            topology.atoms.len(),
            expected_numbers.len(),
            "{label}: atom row count"
        );
        // Direct row-by-row verification against the kept expectation.
        let kept: Vec<(u8, Option<u16>)> = atomic_numbers
            .iter()
            .zip(isotopes.iter())
            .filter(|(z, isotope)| {
                !remove_hydrogens || **z != 1 || isotope.is_some() || all_hydrogen
            })
            .map(|(&z, &isotope)| (z, isotope))
            .collect();
        for (row, (expected_z, expected_isotope)) in topology.atoms.iter().zip(kept.iter()) {
            assert_eq!(row.atomic_number(), *expected_z, "{label}: atomic number");
            assert_eq!(row.isotope(), *expected_isotope, "{label}: isotope");
        }
        if remove_hydrogens && smiles == "[H]N([H])[H]" {
            assert_eq!(topology.atoms[0].atomic_number(), 7, "{label}: N kept");
            assert_eq!(topology.atoms[0].explicit_hydrogens(), 3, "{label}: N H=3");
        } else {
            // The other atom-spec H literals are unchanged under BOTH
            // policies; only ammonia changes row index/count on removal.
            for &(row_index, expected_h) in explicit_h {
                assert_eq!(
                    topology.atoms[row_index].explicit_hydrogens(),
                    expected_h,
                    "{label}: explicit H row {row_index}"
                );
            }
        }
    }

    /// Real parser preparation per policy (raw parse / sanitized),
    /// reusing the existing detached owners exactly like the ring
    /// fixtures; no invented sanitize calls.
    fn prepared_topology(smiles: &str, sanitized: bool) -> TopologyBlock {
        let parsed = cosmolkit_smiles::parse_smiles(
            smiles,
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hydrogens: false,
                ..cosmolkit_smiles::SmilesParseParams::default()
            },
        )
        .expect("frozen case parses");
        if sanitized {
            cosmolkit_core::sanitize_topology(
                &parsed.topology,
                &cosmolkit_core::SanitizeParams::default(),
            )
            .expect("frozen case sanitizes")
            .topology
        } else {
            parsed.topology
        }
    }

    #[test]
    fn descriptor_heteroatoms_domain_input_state_product() {
        // 20 x 2 remove-H policies x 2 raw/sanitized x 2 repeats = 160
        // actual calls to the topology-only owner with frozen literals.
        let mut calls = 0usize;
        for (smiles, expected) in CASES {
            for remove_hydrogens in [false, true] {
                // remove-H policy changes the CONSTRUCTOR preparation
                // (real RemoveHs) exactly like the ring fixtures; the
                // literal count is IDENTICAL under both policies.
                let parsed = cosmolkit_smiles::parse_smiles(
                    smiles,
                    &cosmolkit_smiles::SmilesParseParams {
                        sanitize: false,
                        remove_hydrogens: false,
                        ..cosmolkit_smiles::SmilesParseParams::default()
                    },
                )
                .expect("frozen case parses");
                let topology = if remove_hydrogens {
                    let result = cosmolkit_core::remove_hydrogens_with_params(
                        parsed.topology,
                        parsed.coordinates,
                        parsed.properties,
                        &cosmolkit_core::RemoveHsParams {
                            update_explicit_count: true,
                            sanitize: true,
                            ..cosmolkit_core::RemoveHsParams::default()
                        },
                    )
                    .expect("frozen case RemoveHs");
                    result.topology
                } else {
                    parsed.topology
                };
                for sanitized in [false, true] {
                    let final_topology = if sanitized {
                        // For the sanitized arm run the real ALL-default
                        // sanitize on the already-prepared topology.
                        cosmolkit_core::sanitize_topology(
                            &topology,
                            &cosmolkit_core::SanitizeParams::default(),
                        )
                        .expect("frozen case sanitizes")
                        .topology
                    } else {
                        topology.clone()
                    };
                    for _repeat in 0..2 {
                        let label = format!("{smiles}/rh={remove_hydrogens}/s={sanitized}");
                        // Constructor/element/H/isotope prerequisite BEFORE
                        // every invocation (frozen input identities, never
                        // SUT-derived).
                        assert_input_prerequisites(
                            &label,
                            smiles,
                            &final_topology,
                            remove_hydrogens,
                        );
                        // Per-call fresh topology baseline.
                        let baseline = final_topology.clone();
                        let helper_before = VALENCE_HELPER_ENTRIES.with(std::cell::Cell::get);
                        let count = num_heteroatoms(&final_topology)
                            .unwrap_or_else(|error| panic!("{label}: {error:?}"));
                        calls += 1;
                        assert_eq!(count, expected, "{label}");
                        assert_eq!(final_topology, baseline, "{label}: topology unchanged");
                        // Zero helper reentry: the topology-only owner
                        // never enters the valence helper.
                        assert_eq!(
                            VALENCE_HELPER_ENTRIES.with(std::cell::Cell::get),
                            helper_before,
                            "{label}: zero valence helper entries"
                        );
                    }
                }
            }
        }
        assert_eq!(calls, 160, "exact census");
    }

    #[test]
    fn descriptor_heteroatoms_domain_maxmatches_boundary() {
        // Disconnected real N-row topologies: literal counts
        // [0,1,999,1000,1001,1005] clamp at maxMatches=1000.
        let row_counts = [0usize, 1, 999, 1000, 1001, 1005];
        let expected = [0u32, 1, 999, 1000, 1000, 1000];
        let mut calls = 0usize;
        for (row_index, rows) in row_counts.iter().enumerate() {
            for _repeat in 0..2 {
                let atoms = (0..*rows)
                    .map(|id| {
                        cosmolkit_model::Atom::from_spec(
                            cosmolkit_model::AtomId::new(id),
                            cosmolkit_model::AtomSpec::new(cosmolkit_model::Element::N),
                        )
                    })
                    .collect::<Vec<_>>();
                let topology =
                    TopologyBlock::try_from_parts(atoms, Vec::new(), Vec::new(), Vec::new())
                        .expect("legitimate model constructor rows");
                let count = num_heteroatoms(&topology)
                    .unwrap_or_else(|error| panic!("rows={rows}: {error:?}"));
                calls += 1;
                assert_eq!(count, expected[row_index], "rows={rows}");
            }
        }
        assert_eq!(calls, 12, "exact census");
    }

    #[test]
    fn descriptor_heteroatoms_domain_raw_pentavalent_carbon() {
        // Raw C(C)(C)(C)(C)O: the matcher counts only the O row (1)
        // despite the pentavalent carbon; NO sanitize/valence lookup runs
        // before the tested path (zero valence helper entries).
        let parsed = cosmolkit_smiles::parse_smiles(
            "C(C)(C)(C)(C)O",
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hydrogens: false,
                ..cosmolkit_smiles::SmilesParseParams::default()
            },
        )
        .expect("raw parse");
        let helper_before = VALENCE_HELPER_ENTRIES.with(std::cell::Cell::get);
        let count = num_heteroatoms(&parsed.topology).expect("topology-only path");
        assert_eq!(count, 1, "raw pentavalent C: only O counts");
        assert_eq!(
            VALENCE_HELPER_ENTRIES.with(std::cell::Cell::get),
            helper_before,
            "no valence lookup before the tested path"
        );
    }
}

/// HBD-PUBLIC frozen 104-call domain narrow-owner product plus malformed
/// row regressions.
#[cfg(test)]
mod descriptor_hbd_narrow_tests {
    use super::*;

    /// Frozen 13-case literal table, independently recorded by ROOT from
    /// RDKit 2026.03.1 under BOTH removeHs policies BEFORE any CK run.
    /// Isolated S has two H and is 0; ammonium is 1; never SUT-derived.
    const CASES: [(&str, u32); 13] = [
        ("", 0),
        ("CCO", 1),
        ("NCC(=O)O", 2),
        ("[NH4+]", 1),
        ("c1cc[nH]c1", 1),
        ("c1ccccc1", 0),
        ("CC#CC", 0),
        ("[H]N([H])[H]", 1),
        ("CS", 1),
        ("N", 1),
        ("[OH2]", 0),
        ("S", 0),
        ("COC", 0),
    ];

    #[test]
    fn descriptor_hbd_narrow_input_state_product() {
        // 13 literals x 2 remove-H real-owner routes x 2 forms (narrow +
        // existing prepared) x 2 repeats = 104 real HBD calls.
        let mut calls = 0usize;
        for (smiles, expected) in CASES {
            for remove_hydrogens in [false, true] {
                let parsed = cosmolkit_smiles::parse_smiles(
                    smiles,
                    &cosmolkit_smiles::SmilesParseParams {
                        sanitize: false,
                        remove_hydrogens: false,
                        ..cosmolkit_smiles::SmilesParseParams::default()
                    },
                )
                .expect("frozen case parses");
                let topology = if remove_hydrogens {
                    let result = cosmolkit_core::remove_hydrogens_with_params(
                        parsed.topology,
                        parsed.coordinates,
                        parsed.properties,
                        &cosmolkit_core::RemoveHsParams {
                            update_explicit_count: true,
                            sanitize: true,
                            ..cosmolkit_core::RemoveHsParams::default()
                        },
                    )
                    .expect("frozen case RemoveHs");
                    result.topology
                } else {
                    parsed.topology
                };
                // The public/ring-fixture constructor policy: real
                // ALL-default sanitize on the prepared topology, then the
                // real non-strict assignment owner for the FINAL rows.
                let final_topology = cosmolkit_core::sanitize_topology(
                    &topology,
                    &cosmolkit_core::SanitizeParams::default(),
                )
                .expect("frozen case sanitizes")
                .topology;
                let assignment = cosmolkit_core::assign_valence_for_topology(
                    &final_topology,
                    cosmolkit_core::ValenceModel::RdkitLike,
                )
                .expect("frozen case assignment");
                for form in ["narrow", "prepared"] {
                    for repeat in 0..2 {
                        let label = format!("{smiles}/rh={remove_hydrogens}/{form}#{repeat}");
                        // Fresh per-call input snapshots.
                        let topology_before = final_topology.clone();
                        let assignment_before = assignment.clone();
                        let helper_before = VALENCE_HELPER_ENTRIES.with(std::cell::Cell::get);
                        let count = match form {
                            "narrow" => num_hbd_with_valence(&final_topology, &assignment),
                            _ => {
                                // Existing prepared form: real ring rows
                                // through the real owner, full DescriptorInput.
                                let rings = cosmolkit_core::find_sssr(
                                    &final_topology,
                                    &cosmolkit_core::RingSearchParams::default(),
                                )
                                .expect("frozen case rings");
                                let coordinates = CoordinateBlock::default();
                                let properties = MoleculeProperties::default();
                                let input = DescriptorInput::new(
                                    &final_topology,
                                    &coordinates,
                                    &properties,
                                    &assignment,
                                    &rings,
                                );
                                lipinski::num_hbd_prepared(&input)
                            }
                        }
                        .unwrap_or_else(|error| panic!("{label}: {error:?}"));
                        calls += 1;
                        assert_eq!(count, expected, "{label}: literal output");
                        assert_eq!(
                            final_topology, topology_before,
                            "{label}: topology unchanged"
                        );
                        assert_eq!(
                            assignment, assignment_before,
                            "{label}: assignment unchanged"
                        );
                        // Zero helper reentry: neither form recomputes
                        // valence inside the owner.
                        assert_eq!(
                            VALENCE_HELPER_ENTRIES.with(std::cell::Cell::get),
                            helper_before,
                            "{label}: zero valence helper entries"
                        );
                    }
                }
            }
        }
        assert_eq!(calls, 104, "exact 104-call census");
    }

    #[test]
    fn descriptor_hbd_narrow_rejects_malformed_valence_rows() {
        // Two malformed valence lengths reach the exact Context source
        // error through the ACTUAL narrow call; whole supplied values are
        // preserved. The source pattern never consults ring state, so
        // empty/reset ring rows are not a gate on this path.
        let parsed =
            cosmolkit_smiles::parse_smiles("CCO", &cosmolkit_smiles::SmilesParseParams::default())
                .expect("CCO parses");
        let topology = cosmolkit_core::sanitize_topology(
            &parsed.topology,
            &cosmolkit_core::SanitizeParams::default(),
        )
        .expect("CCO sanitizes")
        .topology;
        // The validator checks explicit_valence first, so the non-target
        // field must carry the VALID length in each case.
        for (field, explicit, implicit) in [
            ("explicit_valence", 2usize, 3usize),
            ("implicit_hydrogens", 3, 5),
        ] {
            let malformed = cosmolkit_core::ValenceAssignment {
                explicit_valence: vec![1; explicit],
                implicit_hydrogens: vec![1; implicit],
            };
            let malformed_before = malformed.clone();
            let topology_before = topology.clone();
            let error = num_hbd_with_valence(&topology, &malformed)
                .err()
                .unwrap_or_else(|| panic!("{field} length must be rejected"));
            let DescriptorError::Search {
                function,
                source: DescriptorSearchCause::Context(context),
            } = &error
            else {
                panic!("expected Search/Context, got {error:?}")
            };
            assert_eq!(*function, "num_hbd", "exact function tag");
            let cosmolkit_search::QueryMatchContextError::ValenceRows {
                field: actual_field,
                expected,
                actual,
            } = context
            else {
                panic!("expected ValenceRows, got {context:?}")
            };
            assert_eq!(*actual_field, field, "exact field");
            assert_eq!(*expected, 3, "expected = CCO atom count");
            assert_eq!(
                *actual,
                if field == "explicit_valence" {
                    explicit
                } else {
                    implicit
                },
                "exact actual length"
            );
            assert_eq!(malformed, malformed_before, "whole valence preserved");
            assert_eq!(topology, topology_before, "whole topology preserved");
        }
    }
}

/// General SMARTS-based amide-bond count.
///
/// Cold convenience form: computes ONE cold canonical assignment and ONE
/// cold SSSR ring set, then delegates to [`lipinski::num_amide_bonds_prepared`]
/// (retained fixed pattern `"C(=[O;!R])N"`, default match parameters;
/// `!R` consumes the ring set).
pub fn num_amide_bonds(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let assignment = valence(topology, "num_amide_bonds")?;
    let rings = ring_info(topology, "num_amide_bonds")?;
    let coordinates = CoordinateBlock::default();
    let properties = MoleculeProperties::default();
    let input = DescriptorInput::new(topology, &coordinates, &properties, &assignment, &rings);
    lipinski::num_amide_bonds_prepared(&input)
}

/// Rotatable-bond count under an explicit definition option.
///
/// Cold convenience form: computes ONE cold canonical assignment and ONE
/// cold SSSR ring set, then delegates to
/// [`rotatable::num_rotatable_bonds_prepared`] (retained fixed patterns,
/// default match parameters). `Default` resolves to the pinned build
/// default (`Strict`); all four options are live.
pub fn num_rotatable_bonds(
    topology: &TopologyBlock,
    options: RotatableBondsOptions,
) -> DescriptorResult<u32> {
    let assignment = valence(topology, "num_rotatable_bonds")?;
    let rings = ring_info(topology, "num_rotatable_bonds")?;
    let coordinates = CoordinateBlock::default();
    let properties = MoleculeProperties::default();
    let input = DescriptorInput::new(topology, &coordinates, &properties, &assignment, &rings);
    rotatable::num_rotatable_bonds_prepared(&input, options)
}

/// Number of SSSR rings.
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one [`rings::num_rings_kernel`] — no
/// valence assignment and no dummy coordinate/property blocks are prepared
/// for this ring-only descriptor.
pub fn num_rings(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_rings")?;
    rings::num_rings_kernel(&rings)
}

/// Number of rings read from CALLER-SUPPLIED initialized `RingInfo` rows
/// (narrow detached boundary; no `DescriptorInput`, no ring recomputation).
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_rings_with_ring_info(ring_info: &cosmolkit_core::RingInfo) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:206-209):
    //   unsigned int calcNumRings(const ROMol &mol) {
    //     return mol.getRingInfo()->numRings();
    //   }
    // RDKit✔️✔️: unsigned int calcNumRings(const ROMol &mol) {
    // RDKit✔️✔️:   return mol.getRingInfo()->numRings();
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrow.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` parameter — the value `mol.getRingInfo()` would return;
    // no `ROMol`, topology, coordinates, properties or valence is taken.
    // Behavior review: the source logic (read the final RingInfo's
    // `numRings()`) is implemented by the REUSED
    // [`rings::num_rings_kernel`], which owns its own in-function anchors
    // and the core `RingInfo::numRings` projection; this delegate neither
    // duplicates that algorithm nor promotes the kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one row-count read plus one typed `u32` conversion — the source
    // pointer-fetch + size-read cost class; no allocation beyond the
    // `Result`, no helper entry.
    rings::num_rings_kernel(ring_info)
}

/// Number of heterocycles (ring rows with ANY non-carbon member atom)
/// read from CALLER-SUPPLIED initialized `RingInfo` rows and the paired
/// atom table (narrow detached boundary; no `DescriptorInput`, no ring
/// recomputation).
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_heterocycles_with_ring_info(
    topology: &TopologyBlock,
    ring_info: &cosmolkit_core::RingInfo,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:233-243):
    //   unsigned int calcNumHeterocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->atomRings()) {
    //       for (auto i : iv) {
    //         if (mol.getAtomWithIdx(i)->getAtomicNum() != 6) {
    //           ++res;
    //           break;
    //         }
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->atomRings()) {
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (mol.getAtomWithIdx(i)->getAtomicNum() != 6) {
    // RDKit✔️✔️:         ++res;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrows.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` (the value `mol.getRingInfo()` would return) plus the
    // `topology.atoms` table its member indices address; no coordinates,
    // properties or valence is taken.
    // Behavior review: the source logic (the `atomRings()` single pass,
    // the `getAtomicNum() != 6` ANY-member predicate and the first-match
    // `break`) is implemented by the REUSED
    // [`rings::num_heterocycles_kernel`], which owns its own in-function
    // anchors; this delegate neither duplicates that algorithm nor
    // promotes the kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one early-exit pass over the supplied rows — O(rows) best case,
    // O(total supplied membership) worst case, the source loop shape; no
    // allocation beyond the `Result`, no helper entry.
    rings::num_heterocycles_kernel(ring_info, &topology.atoms)
}

/// Number of aromatic rings (ring rows whose EVERY member bond carries
/// the aromatic flag) read from CALLER-SUPPLIED initialized `RingInfo`
/// rows and the paired bond table (narrow detached boundary; no
/// `DescriptorInput`, no ring recomputation).
///
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_aromatic_rings_with_ring_info(
    topology: &TopologyBlock,
    ring_info: &cosmolkit_core::RingInfo,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:247-257):
    //   unsigned int calcNumAromaticRings(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       ++res;
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           --res;
    //           break;
    //         }
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAromaticRings(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     ++res;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         --res;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrows.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` (the value `mol.getRingInfo()` would return) plus the
    // `topology.bonds` table its member indices address; no atoms,
    // coordinates, properties or valence is taken.
    // Behavior review: the source logic (the `bondRings()` pass, the
    // pre-increment, the first non-aromatic `--res; break` retraction —
    // a row counts only when EVERY member bond is aromatic) is implemented
    // by the REUSED [`rings::num_aromatic_rings_kernel`], which owns its
    // own in-function anchors; this delegate neither duplicates that
    // algorithm nor promotes the kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one early-exit pass over the supplied rows — O(rows) best case,
    // O(total supplied bond-ring membership) worst case, the source loop
    // shape; no allocation beyond the `Result`, no helper entry.
    rings::num_aromatic_rings_kernel(ring_info, &topology.bonds)
}

/// Number of saturated rings (ring rows whose EVERY member bond is
/// Single-ordered AND non-aromatic) read from CALLER-SUPPLIED initialized
/// `RingInfo` rows and the paired bond table (narrow detached boundary;
/// no `DescriptorInput`, no ring recomputation).
///
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_saturated_rings_with_ring_info(
    topology: &TopologyBlock,
    ring_info: &cosmolkit_core::RingInfo,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:261-272):
    //   unsigned int calcNumSaturatedRings(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       ++res;
    //       for (int i : iv) {
    //         if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    //             mol.getBondWithIdx(i)->getIsAromatic()) {
    //           --res;
    //           break;
    //         }
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumSaturatedRings(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     ++res;
    // RDKit✔️✔️:     for (int i : iv) {
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         --res;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrows.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` (the value `mol.getRingInfo()` would return) plus the
    // `topology.bonds` table its member indices address (both the bond
    // kind and the aromatic flag are read there); no atoms, coordinates,
    // properties or valence is taken.
    // Behavior review: the source logic (the `bondRings()` pass, the
    // optimistic increment, the first-failure `--res; break` retraction on
    // the conjunct `getBondType() != Bond::SINGLE || getIsAromatic()` — a
    // row counts only when EVERY member bond is Single-ordered AND
    // non-aromatic) is implemented by the REUSED
    // [`rings::num_saturated_rings_kernel`], which owns its own
    // in-function anchors; this delegate neither duplicates that
    // algorithm nor promotes the kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one early-exit pass over the supplied rows — O(rows) best case,
    // O(total supplied bond-ring membership) worst case, the source loop
    // shape; no allocation beyond the `Result`, no helper entry.
    rings::num_saturated_rings_kernel(ring_info, &topology.bonds)
}

/// Number of aliphatic rings (ring rows with AT LEAST ONE non-aromatic
/// member bond) read from CALLER-SUPPLIED initialized `RingInfo` rows and
/// the paired bond table (narrow detached boundary; no `DescriptorInput`,
/// no ring recomputation).
///
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_aliphatic_rings_with_ring_info(
    topology: &TopologyBlock,
    ring_info: &cosmolkit_core::RingInfo,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:276-286):
    //   unsigned int calcNumAliphaticRings(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           ++res;
    //           break;
    //         }
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAliphaticRings(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         ++res;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrows.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` (the value `mol.getRingInfo()` would return) plus the
    // `topology.bonds` table whose aromatic flags its member indices
    // address; no atoms, bond orders, coordinates, properties or valence
    // is taken.
    // Behavior review: the source logic (the `bondRings()` pass with NO
    // optimistic increment; `++res` INSIDE the inner loop at the FIRST
    // member bond whose flag is false, then `break` — a row counts once
    // iff it has at least one non-aromatic member, the set-complement of
    // the R03 all-aromatic predicate) is implemented by the REUSED
    // [`rings::num_aliphatic_rings_kernel`], which owns its own
    // in-function anchors; this delegate neither duplicates that
    // algorithm nor promotes the kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one early-exit pass over the supplied rows — O(rows) best case,
    // O(total supplied bond-ring membership) worst case (fully-aromatic
    // rows scan to their end), the source loop shape; no allocation
    // beyond the `Result`, no helper entry.
    rings::num_aliphatic_rings_kernel(ring_info, &topology.bonds)
}

/// Number of aromatic heterocycles (ring rows whose EVERY member bond is
/// aromatic AND at least one member bond has a non-carbon endpoint) read
/// from CALLER-SUPPLIED initialized `RingInfo` rows, the paired bond table
/// and the paired atom table (narrow detached boundary; no
/// `DescriptorInput`, no ring recomputation).
///
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_aromatic_heterocycles_with_ring_info(
    topology: &TopologyBlock,
    ring_info: &cosmolkit_core::RingInfo,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:289-311):
    //   unsigned int calcNumAromaticHeterocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool countIt = false;
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           countIt = false;
    //           break;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sink.
    //         if (!countIt &&
    //             (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //              mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    //           countIt = true;
    //         }
    //       }
    //       if (countIt) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAromaticHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool countIt = false;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sink.
    // RDKit✔️✔️:       if (!countIt &&
    // RDKit✔️✔️:           (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:            mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    // RDKit✔️✔️:         countIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (countIt) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrows.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` (the value `mol.getRingInfo()` would return) plus the
    // `topology.bonds` table (aromatic flags; typed begin()/end() AtomIds)
    // and the `topology.atoms` table those endpoints index; no
    // coordinates, properties or valence is taken.
    // Behavior review: the source logic (per-row `countIt` flag forced
    // false with a break at the FIRST non-aromatic member bond; sticky
    // set at the first member bond with a non-carbon endpoint under the
    // `!countIt` perf guard; increment AFTER the loop — a row counts iff
    // every member bond is aromatic AND at least one has a non-carbon
    // endpoint) is implemented by the REUSED
    // [`rings::num_aromatic_heterocycles_kernel`], which owns its own
    // in-function anchors including the verbatim "doofy" comment; this
    // delegate neither duplicates that algorithm nor promotes the
    // kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one pass with an early break at the first non-aromatic member
    // and a sticky flag skipping endpoint checks after the first hetero
    // bond — the source loop shape; no allocation beyond the `Result`,
    // no helper entry.
    rings::num_aromatic_heterocycles_kernel(ring_info, &topology.bonds, &topology.atoms)
}

/// Number of aromatic carbocycles (ring rows whose EVERY member bond is
/// aromatic AND every endpoint of every member bond is carbon) read from
/// CALLER-SUPPLIED initialized `RingInfo` rows, the paired bond table and
/// the paired atom table (narrow detached boundary; no `DescriptorInput`,
/// no ring recomputation).
///
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_aromatic_carbocycles_with_ring_info(
    topology: &TopologyBlock,
    ring_info: &cosmolkit_core::RingInfo,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:313-333):
    //   unsigned int calcNumAromaticCarbocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool countIt = true;
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           countIt = false;
    //           break;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sync.
    //         if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //             mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    //           countIt = false;
    //           break;
    //         }
    //       }
    //       if (countIt) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAromaticCarbocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool countIt = true;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (countIt) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrows.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` (the value `mol.getRingInfo()` would return) plus the
    // `topology.bonds` table (aromatic flags; typed begin()/end()
    // AtomIds) and the `topology.atoms` table those endpoints index; no
    // coordinates, properties or valence is taken.
    // Behavior review: the source logic (per-row `countIt` starts TRUE;
    // the inner scan breaks with `countIt = false` at EITHER the first
    // non-aromatic member bond OR the first member bond with a non-carbon
    // endpoint; increment after the loop — counts fully-aromatic
    // all-carbon rows, the carbocycle complement of the N06 predicate) is
    // implemented by the REUSED
    // [`rings::num_aromatic_carbocycles_kernel`], which owns its own
    // in-function anchors including the verbatim "time sync" comment
    // variant; this delegate neither duplicates that algorithm nor
    // promotes the kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one pass with an early break at the first failing member (either
    // arm) — the source loop shape; no allocation beyond the `Result`,
    // no helper entry.
    rings::num_aromatic_carbocycles_kernel(ring_info, &topology.bonds, &topology.atoms)
}

/// Number of aliphatic heterocycles (ring rows with AT LEAST ONE
/// non-aromatic member bond AND at least one member bond with a
/// non-carbon endpoint, at any positions) read from CALLER-SUPPLIED
/// initialized `RingInfo` rows, the paired bond table and the paired
/// atom table (narrow detached boundary; no `DescriptorInput`, no ring
/// recomputation).
///
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_aliphatic_heterocycles_with_ring_info(
    topology: &TopologyBlock,
    ring_info: &cosmolkit_core::RingInfo,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:335-358):
    //   unsigned int calcNumAliphaticHeterocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool hasAliph = false;
    //       bool hasHetero = false;
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           hasAliph = true;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sink.
    //         if (!hasHetero &&
    //             (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //              mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    //           hasHetero = true;
    //         }
    //       }
    //       if (hasHetero && hasAliph) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAliphaticHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool hasAliph = false;
    // RDKit✔️✔️:     bool hasHetero = false;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         hasAliph = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sink.
    // RDKit✔️✔️:       if (!hasHetero &&
    // RDKit✔️✔️:           (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:            mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    // RDKit✔️✔️:         hasHetero = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (hasHetero && hasAliph) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrows.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` (the value `mol.getRingInfo()` would return) plus the
    // `topology.bonds` table (aromatic flags; typed begin()/end()
    // AtomIds) and the `topology.atoms` table those endpoints index; no
    // coordinates, properties or valence is taken.
    // Behavior review: the source logic (two independent flags —
    // `hasAliph` set by ANY non-aromatic member, sticky guarded
    // `hasHetero` set by any non-carbon-endpoint member; NO break, full
    // member scan; BOTH flags after the loop => ++) is implemented by the
    // REUSED [`rings::num_aliphatic_heterocycles_kernel`], which owns its
    // own in-function anchors including the verbatim doofy comment; this
    // delegate neither duplicates that algorithm nor promotes the
    // kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one full-row two-flag scan — O(total supplied bond-ring
    // membership), the source loop shape (no early exit in the source);
    // no allocation beyond the `Result`, no helper entry.
    rings::num_aliphatic_heterocycles_kernel(ring_info, &topology.bonds, &topology.atoms)
}

/// Number of aliphatic carbocycles (ring rows with AT LEAST ONE
/// non-aromatic member bond AND NO member bond with a non-carbon
/// endpoint) read from CALLER-SUPPLIED initialized `RingInfo` rows, the
/// paired bond table and the paired atom table (narrow detached
/// boundary; no `DescriptorInput`, no ring recomputation).
///
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_aliphatic_carbocycles_with_ring_info(
    topology: &TopologyBlock,
    ring_info: &cosmolkit_core::RingInfo,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:361-382):
    //   unsigned int calcNumAliphaticCarbocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool hasAliph = false;
    //       bool hasHetero = false;
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           hasAliph = true;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sync.
    //         if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //             mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    //           hasHetero = true;
    //           break;
    //         }
    //       }
    //       if (hasAliph && !hasHetero) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAliphaticCarbocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool hasAliph = false;
    // RDKit✔️✔️:     bool hasHetero = false;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         hasAliph = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    // RDKit✔️✔️:         hasHetero = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (hasAliph && !hasHetero) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrows.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` (the value `mol.getRingInfo()` would return) plus the
    // `topology.bonds` table (aromatic flags; typed begin()/end()
    // AtomIds) and the `topology.atoms` table those endpoints index; no
    // coordinates, properties or valence is taken.
    // Behavior review: the source logic (a non-aromatic member sets
    // `hasAliph` WITHOUT break; a non-carbon-endpoint member sets
    // `hasHetero` AND breaks immediately — no sticky guard; after the
    // loop `hasAliph && !hasHetero` increments — a row counts iff it has
    // at least one non-aromatic member and NO hetero-endpoint member at
    // all; the break is outcome-preserving) is implemented by the REUSED
    // [`rings::num_aliphatic_carbocycles_kernel`], which owns its own
    // in-function anchors including the verbatim "time sync" comment;
    // this delegate neither duplicates that algorithm nor promotes the
    // kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one pass with an early break at the first hetero-endpoint member
    // — the source loop shape; no allocation beyond the `Result`, no
    // helper entry.
    rings::num_aliphatic_carbocycles_kernel(ring_info, &topology.bonds, &topology.atoms)
}

/// Number of saturated heterocycles (ring rows whose EVERY member bond is
/// Single-ordered AND non-aromatic AND at least one member bond has a
/// non-carbon endpoint) read from CALLER-SUPPLIED initialized `RingInfo`
/// rows, the paired bond table and the paired atom table (narrow detached
/// boundary; no `DescriptorInput`, no ring recomputation).
///
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_saturated_heterocycles_with_ring_info(
    topology: &TopologyBlock,
    ring_info: &cosmolkit_core::RingInfo,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:385-407):
    //   unsigned int calcNumSaturatedHeterocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool countIt = false;
    //       for (auto i : iv) {
    //         if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    //             mol.getBondWithIdx(i)->getIsAromatic()) {
    //           countIt = false;
    //           break;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sync.
    //         if (!countIt &&
    //             (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //              mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    //           countIt = true;
    //         }
    //       }
    //       if (countIt) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumSaturatedHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool countIt = false;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
    // RDKit✔️✔️:       if (!countIt &&
    // RDKit✔️✔️:           (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:            mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    // RDKit✔️✔️:         countIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (countIt) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrows.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` (the value `mol.getRingInfo()` would return) plus the
    // `topology.bonds` table (kinds AND flags; typed begin()/end()
    // AtomIds) and the `topology.atoms` table those endpoints index; no
    // coordinates, properties or valence is taken.
    // Behavior review: the source logic (retraction at the first member
    // failing the saturated conjunct `getBondType() != Bond::SINGLE ||
    // getIsAromatic()`; the N06-style sticky guarded hetero arm setting
    // countIt at the first non-carbon-endpoint member; increment after
    // the loop) is implemented by the REUSED
    // [`rings::num_saturated_heterocycles_kernel`], which owns its own
    // in-function anchors including the verbatim "time sync" comment;
    // this delegate neither duplicates that algorithm nor promotes the
    // kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one pass with an early break at the first non-saturated member
    // and a sticky flag skipping endpoint checks after the first hetero
    // bond — the source loop shape; no allocation beyond the `Result`,
    // no helper entry.
    rings::num_saturated_heterocycles_kernel(ring_info, &topology.bonds, &topology.atoms)
}

/// Number of saturated carbocycles (ring rows whose EVERY member bond is
/// Single-ordered AND non-aromatic AND every endpoint of every member
/// bond is carbon) read from CALLER-SUPPLIED initialized `RingInfo`
/// rows, the paired bond table and the paired atom table (narrow
/// detached boundary; no `DescriptorInput`, no ring recomputation).
///
/// The complete verbatim pinned source body and the separate
/// borrow/behavior/complexity reviews live IN the function body beside
/// the one reused call.
pub fn num_saturated_carbocycles_with_ring_info(
    topology: &TopologyBlock,
    ring_info: &cosmolkit_core::RingInfo,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:410-432):
    //   unsigned int calcNumSaturatedCarbocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool countIt = true;
    //       for (auto i : iv) {
    //         if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    //             mol.getBondWithIdx(i)->getIsAromatic()) {
    //           countIt = false;
    //           break;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sync.
    //         if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //             mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    //           countIt = false;
    //           break;
    //         }
    //       }
    //       if (countIt) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumSaturatedCarbocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool countIt = true;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (countIt) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // Borrow review: this delegate adds ONLY the caller-supplied borrows.
    // The source's `const ROMol &mol` reference role maps to the borrowed
    // `ring_info` (the value `mol.getRingInfo()` would return) plus the
    // `topology.bonds` table (kinds AND flags; typed begin()/end()
    // AtomIds) and the `topology.atoms` table those endpoints index; no
    // coordinates, properties or valence is taken.
    // Behavior review: the source logic (`countIt` starts TRUE; breaks
    // false at EITHER the first member failing the saturated conjunct
    // `getBondType() != Bond::SINGLE || getIsAromatic()` OR the first
    // member with a non-carbon endpoint; increment after the loop —
    // counts all-Single non-aromatic all-carbon rows, the saturated
    // carbocycle complement of N10) is implemented by the REUSED
    // [`rings::num_saturated_carbocycles_kernel`], which owns its own
    // in-function anchors including the verbatim "time sync" comment;
    // this delegate neither duplicates that algorithm nor promotes the
    // kernel's markers.
    // Complexity review: borrow-only transport; the single qualified call
    // is one pass with an early break at the first failing member
    // (either arm) — the source loop shape; no allocation beyond the
    // `Result`, no helper entry.
    rings::num_saturated_carbocycles_kernel(ring_info, &topology.bonds, &topology.atoms)
}

/// Number of heterocycles (rings with ANY non-carbon member atom).
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one [`rings::num_heterocycles_kernel`] —
/// no valence assignment and no dummy coordinate/property blocks.
pub fn num_heterocycles(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_heterocycles")?;
    rings::num_heterocycles_kernel(&rings, &topology.atoms)
}

/// Number of aromatic rings (ring rows whose EVERY member bond carries the
/// aromatic flag).
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one [`rings::num_aromatic_rings_kernel`] —
/// no valence assignment and no dummy coordinate/property blocks. The
/// predicate reads the supplied bond flags only, never bond orders.
pub fn num_aromatic_rings(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_aromatic_rings")?;
    rings::num_aromatic_rings_kernel(&rings, &topology.bonds)
}

/// Number of saturated rings (ring rows whose EVERY member bond is Single
/// and non-aromatic).
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one [`rings::num_saturated_rings_kernel`]
/// — no valence assignment and no dummy coordinate/property blocks. The
/// predicate reads the supplied bond orders and flags only.
pub fn num_saturated_rings(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_saturated_rings")?;
    rings::num_saturated_rings_kernel(&rings, &topology.bonds)
}

/// Number of aliphatic rings (ring rows with AT LEAST ONE non-aromatic
/// member bond; the complement of the all-aromatic count).
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one [`rings::num_aliphatic_rings_kernel`]
/// — no valence assignment and no dummy coordinate/property blocks. The
/// predicate reads the supplied bond flags only, never bond orders.
pub fn num_aliphatic_rings(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_aliphatic_rings")?;
    rings::num_aliphatic_rings_kernel(&rings, &topology.bonds)
}

/// Number of aromatic heterocycles (fully aromatic rows with at least one
/// non-carbon member bond endpoint).
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one
/// [`rings::num_aromatic_heterocycles_kernel`] — no valence assignment
/// and no dummy coordinate/property blocks.
pub fn num_aromatic_heterocycles(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_aromatic_heterocycles")?;
    rings::num_aromatic_heterocycles_kernel(&rings, &topology.bonds, &topology.atoms)
}

/// Number of aromatic carbocycles (fully aromatic all-carbon rows).
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one
/// [`rings::num_aromatic_carbocycles_kernel`] — no valence assignment
/// and no dummy coordinate/property blocks.
pub fn num_aromatic_carbocycles(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_aromatic_carbocycles")?;
    rings::num_aromatic_carbocycles_kernel(&rings, &topology.bonds, &topology.atoms)
}

/// Number of aliphatic heterocycles (rows with at least one non-aromatic
/// member bond AND at least one hetero-endpoint member bond).
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one
/// [`rings::num_aliphatic_heterocycles_kernel`] — no valence assignment
/// and no dummy coordinate/property blocks.
pub fn num_aliphatic_heterocycles(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_aliphatic_heterocycles")?;
    rings::num_aliphatic_heterocycles_kernel(&rings, &topology.bonds, &topology.atoms)
}

/// Number of aliphatic carbocycles (rows with at least one non-aromatic
/// member bond and NO hetero-endpoint member bond).
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one
/// [`rings::num_aliphatic_carbocycles_kernel`] — no valence assignment
/// and no dummy coordinate/property blocks.
pub fn num_aliphatic_carbocycles(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_aliphatic_carbocycles")?;
    rings::num_aliphatic_carbocycles_kernel(&rings, &topology.bonds, &topology.atoms)
}

/// Number of saturated heterocycles (all-Single non-aromatic rows with at
/// least one hetero-endpoint member bond).
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one
/// [`rings::num_saturated_heterocycles_kernel`] — no valence assignment
/// and no dummy coordinate/property blocks.
pub fn num_saturated_heterocycles(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_saturated_heterocycles")?;
    rings::num_saturated_heterocycles_kernel(&rings, &topology.bonds, &topology.atoms)
}

/// Number of saturated carbocycles (all-Single non-aromatic all-carbon
/// rows).
///
/// Cold convenience form: computes ONLY the ONE cold canonical SSSR ring
/// assignment and delegates to the one
/// [`rings::num_saturated_carbocycles_kernel`] — no valence assignment
/// and no dummy coordinate/property blocks.
pub fn num_saturated_carbocycles(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_saturated_carbocycles")?;
    rings::num_saturated_carbocycles_kernel(&rings, &topology.bonds, &topology.atoms)
}

/// Number of spiro atoms (atoms that are the SINGLE shared member of a
/// ring-row pair).
///
/// Cold convenience form: the ONE canonical SSSR assignment IS the
/// source's ensure-SSSR acquisition step (the absent-ringInfo arm of
/// `calcNumSpiroAtoms`); the pair scan then runs on those rows. No
/// valence assignment and no dummy coordinate/property blocks.
pub fn num_spiro_atoms(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_spiro_atoms")?;
    let mut atoms = Vec::new();
    rings::spiro_atom_ids_kernel(&rings, &mut atoms)
}

/// The ordered unique spiro atom IDs (push order = row-pair scan order),
/// mirroring the source's optional out-parameter with a fresh vector.
///
/// Cold convenience form: same ONE canonical SSSR acquisition as
/// [`num_spiro_atoms`]; no valence assignment and no dummy blocks.
pub fn spiro_atom_ids(topology: &TopologyBlock) -> DescriptorResult<Vec<AtomId>> {
    let rings = ring_info(topology, "num_spiro_atoms")?;
    let mut atoms = Vec::new();
    rings::spiro_atom_ids_kernel(&rings, &mut atoms)?;
    Ok(atoms)
}

/// Number of bridgehead atoms (single-incidence endpoint atoms of
/// multi-bond-shared ring-row pairs).
///
/// Cold convenience form: the ONE canonical SSSR assignment IS the
/// source's ensure-SSSR acquisition step; no valence assignment and no
/// dummy coordinate/property blocks.
pub fn num_bridgehead_atoms(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let rings = ring_info(topology, "num_bridgehead_atoms")?;
    let mut atoms = Vec::new();
    rings::bridgehead_atom_ids_kernel(&rings, &topology.bonds, topology.atoms.len(), &mut atoms)
}

/// The ordered unique bridgehead atom IDs (pair-scan order, ascending
/// atom index within a pair), mirroring the source's optional
/// out-parameter with a fresh vector.
pub fn bridgehead_atom_ids(topology: &TopologyBlock) -> DescriptorResult<Vec<AtomId>> {
    let rings = ring_info(topology, "num_bridgehead_atoms")?;
    let mut atoms = Vec::new();
    rings::bridgehead_atom_ids_kernel(&rings, &topology.bonds, topology.atoms.len(), &mut atoms)?;
    Ok(atoms)
}

/// Number of potential atom stereo centers (possible-chirality atoms).
///
/// Cold convenience form: a bare topology carries NO `_StereochemDone`
/// property, so the source's absent arm always applies — ONE legacy
/// stereo assignment (cleanIt/force/flagPossible all true) on the
/// canonical prepared input, then the count. NOTE: unlike the ring-only
/// descriptors, this cold path USES the valence helper BY SOURCE
/// NECESSITY — `assignStereochemistry` consumes the valence/rank rows;
/// there is no source path that assigns stereochemistry without them.
pub fn num_atom_stereo_centers(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let assignment = valence(topology, "num_atom_stereo_centers")?;
    let rings = ring_info(topology, "num_atom_stereo_centers")?;
    let coordinates = CoordinateBlock::default();
    let properties = MoleculeProperties::default();
    let input = DescriptorInput::new(topology, &coordinates, &properties, &assignment, &rings);
    stereo::num_atom_stereo_centers_prepared(&input)
}

/// Number of potential-but-unspecified atom stereo centers.
///
/// Cold convenience form: identical absent-arm dispatch as
/// [`num_atom_stereo_centers`]; counts possible-chirality atoms whose
/// chiral tag is still `Unspecified`.
pub fn num_unspecified_atom_stereo_centers(topology: &TopologyBlock) -> DescriptorResult<u32> {
    let assignment = valence(topology, "num_unspecified_atom_stereo_centers")?;
    let rings = ring_info(topology, "num_unspecified_atom_stereo_centers")?;
    let coordinates = CoordinateBlock::default();
    let properties = MoleculeProperties::default();
    let input = DescriptorInput::new(topology, &coordinates, &properties, &assignment, &rings);
    stereo::num_unspecified_atom_stereo_centers_prepared(&input)
}

/// Cold convenience TPSA from a bare topology (uncached arm).
///
/// Builds the canonical prepared input (ONE valence + ONE ring
/// assignment — the fact builder consumes both, the same
/// source-necessity as the stereo unit) and delegates to the one
/// [`tpsa::tpsa_atom_contribs_kernel`]. The cached/scalar-dispatch
/// entrypoint is the canonical [`tpsa`] on [`DescriptorInput`] +
/// [`DescriptorComputedState`].
pub fn tpsa_from_topology(topology: &TopologyBlock, include_sand_p: bool) -> DescriptorResult<f64> {
    let assignment = valence(topology, "tpsa")?;
    let rings = ring_info(topology, "tpsa")?;
    let mut contribs = vec![0.0f64; topology.atoms.len()];
    tpsa::tpsa_atom_contribs_kernel(topology, &assignment, &rings, include_sand_p, &mut contribs)
}

/// Topological polar surface area with the detached typed cache (the
/// source `calcTPSA` scalar dispatch, MolSurf.cpp:347-360).
///
/// Borrows an ALREADY-prepared [`DescriptorInput`], dispatches through
/// the one scalar owner, and reads/mutates only the caller's detached
/// [`DescriptorComputedState`]. This function constructs neither the
/// input nor a default state and does not force recompute: the `force`
/// flag is forwarded unchanged.
///
/// # Exact canonical access
///
/// The canonical domain boundary is exactly the two re-exported
/// entrypoints plus this scalar wrapper:
///
/// ```
/// use cosmolkit_descriptors::{
///     DescriptorComputedState, DescriptorInput, DescriptorResult, tpsa,
///     tpsa_contributions,
/// };
///
/// let scalar: fn(
///     &DescriptorInput<'_>,
///     bool,
///     bool,
///     &mut DescriptorComputedState,
/// ) -> DescriptorResult<f64> = tpsa;
/// let contributions: fn(
///     &DescriptorInput<'_>,
///     bool,
///     bool,
///     &mut DescriptorComputedState,
/// ) -> DescriptorResult<Vec<f64>> = tpsa_contributions;
/// ```
///
/// The implementation module stays private to the crate:
///
/// ```compile_fail,E0603
/// // Private implementation-module access: `tpsa` is not a public
/// // module and `tpsa_scalar` is crate-private. This proof is module
/// // privacy only — not full operation/runtime isolation.
/// let f = cosmolkit_descriptors::tpsa::tpsa_scalar;
/// ```
pub fn tpsa(
    input: &DescriptorInput<'_>,
    include_sulfur_phosphorus: bool,
    force: bool,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<f64> {
    tpsa::tpsa_scalar(input, include_sulfur_phosphorus, force, state)
}

/// Calculates average molecular weight using all atoms and implicit Hs.
#[must_use]
pub fn molecular_weight(topology: &TopologyBlock) -> DescriptorResult<f64> {
    molecular_weight_with_options(topology, false)
}

/// RDKit-compatible average molecular weight with an explicit heavy-atom mode.
#[must_use]
pub fn molecular_weight_with_options(
    topology: &TopologyBlock,
    only_heavy: bool,
) -> DescriptorResult<f64> {
    molecular_weight_with_valence(topology, only_heavy, None)
}

/// Detached mass kernel. A supplied final assignment is borrowed, never cloned
/// or recomputed. Cold detached calls explicitly use the canonical core owner.
pub fn molecular_weight_with_valence(
    topology: &TopologyBlock,
    only_heavy: bool,
    assignment: Option<&cosmolkit_core::ValenceAssignment>,
) -> DescriptorResult<f64> {
    // RDKit✔️✔️: double calcAMW(const ROMol &mol, bool onlyHeavy) {
    // RDKit✔️✔️:   return MolOps::getAvgMolWt(mol, onlyHeavy);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double getAvgMolWt(const ROMol &mol, bool onlyHeavy) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   for (const auto &atom : mol.atoms()) {
    // RDKit✔️✔️:     if (!onlyHeavy || atom->getAtomicNum() != 1) res += atom->getMass();
    // RDKit✔️✔️:     if (!onlyHeavy) res += atom->getTotalNumHs() * table->getAtomicWeight(1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    validate_topology(topology, "molecular_weight")?;
    let assignment = if only_heavy {
        None
    } else {
        Some(prepared_valence(topology, assignment, "molecular_weight")?)
    };
    let hydrogen_mass = descriptor_atomic_mass(Element::H, None, "molecular_weight")?;
    let mut result = 0.0;
    for atom in &topology.atoms {
        if !only_heavy || atom.atomic_number() != 1 {
            result += descriptor_atomic_mass(atom.element(), atom.isotope(), "molecular_weight")?;
        }
        if !only_heavy {
            let hydrogens = total_hydrogen_count_from_validated(
                topology,
                assignment.as_ref().expect("H-sensitive branch"),
                atom.id(),
                false,
            )
            .map_err(|source| DescriptorError::Valence {
                function: "molecular_weight",
                source,
            })?;
            result += f64::from(hydrogens) * hydrogen_mass;
        }
    }
    Ok(result)
}

/// Calculates exact molecular weight using the most common isotope mass.
#[must_use]
pub fn exact_molecular_weight(topology: &TopologyBlock) -> DescriptorResult<f64> {
    exact_molecular_weight_with_options(topology, false)
}

/// RDKit-compatible exact molecular weight with an explicit heavy-atom mode.
#[must_use]
pub fn exact_molecular_weight_with_options(
    topology: &TopologyBlock,
    only_heavy: bool,
) -> DescriptorResult<f64> {
    exact_molecular_weight_with_valence(topology, only_heavy, None)
}

/// Detached exact-mass kernel with borrowed prepared valence.
pub fn exact_molecular_weight_with_valence(
    topology: &TopologyBlock,
    only_heavy: bool,
    assignment: Option<&cosmolkit_core::ValenceAssignment>,
) -> DescriptorResult<f64> {
    // RDKit✔️✔️: double calcExactMW(const ROMol &mol, bool onlyHeavy) {
    // RDKit✔️✔️:   return MolOps::getExactMolWt(mol, onlyHeavy);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: const double electronMass = 0.00054857991;
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (atNum != 1 || !onlyHeavy) res += table->getMostCommonIsotopeMass(atNum);
    // RDKit✔️✔️:     res -= constants::electronMass * atom->getFormalCharge();
    // RDKit✔️✔️:     if (!onlyHeavy) nHsToCount += atom->getTotalNumHs(false);
    // RDKit✔️✔️:   }
    validate_topology(topology, "exact_molecular_weight")?;
    let assignment = if only_heavy {
        None
    } else {
        Some(prepared_valence(
            topology,
            assignment,
            "exact_molecular_weight",
        )?)
    };
    let mut result = 0.0;
    let mut hydrogens_to_count = 0_i32;
    for atom in &topology.atoms {
        if atom.atomic_number() != 1 || !only_heavy {
            result += if atom.isotope().is_none() {
                most_common_isotope_mass(atom.element())
            } else {
                descriptor_atomic_mass(atom.element(), atom.isotope(), "exact_molecular_weight")?
            };
            result -= RDKIT_ELECTRON_MASS * f64::from(atom.formal_charge());
        }
        if !only_heavy {
            let hydrogens = total_hydrogen_count_from_validated(
                topology,
                assignment.as_ref().expect("H-sensitive branch"),
                atom.id(),
                false,
            )
            .map_err(|source| DescriptorError::Valence {
                function: "exact_molecular_weight",
                source,
            })? as i32;
            // Signed overflow is source-undefined, not an alternative H count.
            hydrogens_to_count = hydrogens_to_count.checked_add(hydrogens).ok_or(
                DescriptorError::CountOverflow {
                    function: "exact_molecular_weight",
                    field: "hydrogens",
                },
            )?;
        }
    }
    if !only_heavy {
        result += f64::from(hydrogens_to_count) * most_common_isotope_mass(Element::H);
    }
    Ok(result)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct FormulaKey {
    isotope: u32,
    symbol: &'static str,
}

fn hill_compare(left: &FormulaKey, right: &FormulaKey) -> Ordering {
    if left.symbol == "C" {
        return if right.symbol == "C" {
            left.isotope.cmp(&right.isotope)
        } else {
            Ordering::Less
        };
    }
    if right.symbol == "C" {
        return Ordering::Greater;
    }
    if left.symbol == "H" {
        return if right.symbol == "H" {
            left.isotope.cmp(&right.isotope)
        } else {
            Ordering::Less
        };
    }
    if right.symbol == "H" {
        return Ordering::Greater;
    }
    if left.symbol == "D" {
        return Ordering::Less;
    }
    if right.symbol == "D" {
        return Ordering::Greater;
    }
    if left.symbol == "T" {
        return Ordering::Less;
    }
    if right.symbol == "T" {
        return Ordering::Greater;
    }
    left.cmp(right)
}

/// Calculates the Hill-ordered molecular formula.
#[must_use]
pub fn molecular_formula(topology: &TopologyBlock) -> DescriptorResult<String> {
    molecular_formula_with_options(topology, false, false)
}

/// RDKit-compatible formula generation with isotope formatting controls.
#[must_use]
pub fn molecular_formula_with_options(
    topology: &TopologyBlock,
    separate_isotopes: bool,
    abbreviate_h_isotopes: bool,
) -> DescriptorResult<String> {
    molecular_formula_with_valence(topology, separate_isotopes, abbreviate_h_isotopes, None)
}

/// Detached formula kernel with borrowed prepared valence.
pub fn molecular_formula_with_valence(
    topology: &TopologyBlock,
    separate_isotopes: bool,
    abbreviate_h_isotopes: bool,
    assignment: Option<&cosmolkit_core::ValenceAssignment>,
) -> DescriptorResult<String> {
    // RDKit✔️✔️: std::string getMolFormula(const ROMol &mol, bool separateIsotopes,
    // RDKit✔️✔️:                           bool abbreviateHIsotopes) {
    // RDKit✔️✔️:   std::map<std::pair<unsigned int, std::string>, unsigned int> counts;
    // RDKit✔️✔️:   unsigned int nHs = 0;
    validate_topology(topology, "molecular_formula")?;
    let assignment = prepared_valence(topology, assignment, "molecular_formula")?;
    let mut counts = BTreeMap::<FormulaKey, u32>::new();
    let mut charge = 0_i32;
    let mut hydrogens = 0_u32;
    for atom in &topology.atoms {
        let atomic_number = atom.atomic_number();
        let mut key = FormulaKey {
            isotope: 0,
            symbol: rdkit_element_symbol(atomic_number).map_err(|error| {
                DescriptorError::Unsupported {
                    function: "molecular_formula",
                    detail: error.to_string(),
                }
            })?,
        };
        if separate_isotopes {
            let isotope = atom.isotope().map(u32::from).unwrap_or(0);
            if abbreviate_h_isotopes && atomic_number == 1 && (isotope == 2 || isotope == 3) {
                key.symbol = if isotope == 2 { "D" } else { "T" };
            } else {
                key.isotope = isotope;
            }
        }
        *counts.entry(key).or_insert(0) += 1;
        let count = total_hydrogen_count_from_validated(topology, &assignment, atom.id(), false)
            .map_err(|source| DescriptorError::Valence {
                function: "molecular_formula",
                source,
            })?;
        hydrogens = hydrogens
            .checked_add(count)
            .ok_or(DescriptorError::CountOverflow {
                function: "molecular_formula",
                field: "hydrogens",
            })?;
        charge += i32::from(atom.formal_charge());
    }
    if hydrogens != 0 {
        *counts
            .entry(FormulaKey {
                isotope: 0,
                symbol: "H",
            })
            .or_insert(0) += hydrogens;
    }
    let mut keys = counts.keys().copied().collect::<Vec<_>>();
    keys.sort_by(hill_compare);
    let mut result = String::new();
    for key in keys {
        if key.isotope > 0 {
            result.push('[');
            result.push_str(&key.isotope.to_string());
            result.push_str(key.symbol);
            result.push(']');
        } else {
            result.push_str(key.symbol);
        }
        if let Some(count) = counts.get(&key).copied()
            && count > 1
        {
            result.push_str(&count.to_string());
        }
    }
    if charge > 0 {
        result.push('+');
        if charge > 1 {
            result.push_str(&charge.to_string());
        }
    } else if charge < 0 {
        result.push('-');
        if charge < -1 {
            result.push_str(&(-charge).to_string());
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element,
    };

    fn detached_ethanol() -> TopologyBlock {
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            ),
        ];
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    fn smiles_topology(smiles: &str) -> TopologyBlock {
        cosmolkit_smiles::parse_smiles(smiles, &Default::default())
            .expect("parse detached descriptor fixture")
            .topology
    }

    #[test]
    fn detached_descriptors_match_rdkit_ethanol_fixture() {
        let topology = detached_ethanol();
        assert_eq!(molecular_formula(&topology).unwrap(), "C2H6O");
        assert!((molecular_weight(&topology).unwrap() - 46.069).abs() < 1e-9);
        assert!((exact_molecular_weight(&topology).unwrap() - 46.041864812).abs() < 1e-9);
    }

    #[test]
    fn valence_transport_prepared_kernels_borrow_and_reject_missing_rows() {
        let topology = detached_ethanol();
        let assignment = valence(&topology, "fixture").unwrap();
        let borrowed = prepared_valence(&topology, Some(&assignment), "fixture").unwrap();
        assert!(matches!(borrowed, Cow::Borrowed(_)));
        assert!(std::ptr::eq(borrowed.as_ref(), &assignment));
        for heavy in [false, true] {
            assert_eq!(
                molecular_weight_with_valence(&topology, heavy, Some(&assignment)),
                molecular_weight_with_options(&topology, heavy)
            );
            assert_eq!(
                exact_molecular_weight_with_valence(&topology, heavy, Some(&assignment)),
                exact_molecular_weight_with_options(&topology, heavy)
            );
        }
        for separate in [false, true] {
            for abbreviate in [false, true] {
                assert_eq!(
                    molecular_formula_with_valence(
                        &topology,
                        separate,
                        abbreviate,
                        Some(&assignment)
                    ),
                    molecular_formula_with_options(&topology, separate, abbreviate)
                );
            }
        }
        for explicit in [false, true] {
            for length in [0, 2, 4] {
                let mut malformed = assignment.clone();
                if explicit {
                    malformed.explicit_valence.resize(length, 0);
                } else {
                    malformed.implicit_hydrogens.resize(length, 0);
                }
                assert!(molecular_weight_with_valence(&topology, false, Some(&malformed)).is_err());
                assert!(
                    exact_molecular_weight_with_valence(&topology, false, Some(&malformed))
                        .is_err()
                );
                // Heavy-only kernels do not consume H state in the reference.
                assert!(molecular_weight_with_valence(&topology, true, Some(&malformed)).is_ok());
                assert!(
                    exact_molecular_weight_with_valence(&topology, true, Some(&malformed)).is_ok()
                );
                for separate in [false, true] {
                    for abbreviate in [false, true] {
                        assert!(
                            molecular_formula_with_valence(
                                &topology,
                                separate,
                                abbreviate,
                                Some(&malformed)
                            )
                            .is_err()
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn valence_transport_hydrogen_getter_preserves_no_implicit_and_errors() {
        use std::error::Error;
        for no_implicit in [false, true] {
            let topology = TopologyBlock::try_from_parts(
                vec![Atom::from_spec(
                    AtomId::new(0),
                    AtomSpec::new(Element::C).with_no_implicit(no_implicit),
                )],
                vec![],
                vec![],
                vec![],
            )
            .unwrap();
            for implicit in [-1, 0, 1] {
                let assignment = cosmolkit_core::ValenceAssignment {
                    explicit_valence: vec![0],
                    implicit_hydrogens: vec![implicit],
                };
                let average = molecular_weight_with_valence(&topology, false, Some(&assignment));
                let exact =
                    exact_molecular_weight_with_valence(&topology, false, Some(&assignment));
                let formula =
                    molecular_formula_with_valence(&topology, false, false, Some(&assignment));
                if !no_implicit && implicit < 0 {
                    for error in [
                        average.unwrap_err(),
                        exact.unwrap_err(),
                        formula.unwrap_err(),
                    ] {
                        assert!(matches!(&error,DescriptorError::Valence{source:
                            cosmolkit_core::ValenceError::ImplicitValenceCacheNotInitialized{atom},..}
                            if atom.index()==0));
                        assert!(
                            error
                                .source()
                                .unwrap()
                                .downcast_ref::<cosmolkit_core::ValenceError>()
                                .is_some()
                        );
                    }
                } else {
                    let h = if no_implicit { 0 } else { implicit };
                    assert_eq!(average.unwrap(), 12.011 + f64::from(h) * 1.008);
                    assert_eq!(
                        exact.unwrap(),
                        12.0 + f64::from(h) * most_common_isotope_mass(Element::H)
                    );
                    assert_eq!(formula.unwrap(), if h == 0 { "C" } else { "CH" });
                }
            }
        }
    }

    #[test]
    fn detached_lipinski_and_atom_counts_match_pinned_rdkit_cases() {
        const CASES: [(&str, [u32; 4], f64); 6] = [
            ("CCO", [1, 1, 3, 9], 1.0),
            ("NC(=O)C", [2, 2, 4, 9], 0.5),
            ("NC(=O)N", [3, 4, 4, 8], 0.0),
            ("NCC(=O)O", [3, 3, 5, 10], 0.5),
            ("c1ncc[nH]1", [2, 1, 5, 9], 0.0),
            ("[Na+].[Cl-]", [0, 0, 2, 2], 0.0),
        ];

        for (smiles, expected, expected_fraction) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                [
                    lipinski_hba(&topology).unwrap(),
                    lipinski_hbd(&topology).unwrap(),
                    num_heavy_atoms(&topology).unwrap(),
                    num_atoms(&topology).unwrap(),
                ],
                expected,
                "{smiles}"
            );
            assert_eq!(
                fraction_csp3(&topology).unwrap(),
                expected_fraction,
                "{smiles}"
            );
        }
    }

    #[test]
    fn detached_lipinski_donor_count_includes_explicit_hydrogen_neighbors() {
        let mut atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::N)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        atoms.extend(
            (2..4).map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::H))),
        );
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(2),
                BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
            ),
        ];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };

        assert_eq!(lipinski_hba(&topology), Ok(1));
        assert_eq!(lipinski_hbd(&topology), Ok(2));
        assert_eq!(num_heavy_atoms(&topology), Ok(2));
        assert_eq!(num_atoms(&topology), Ok(7));
        assert_eq!(fraction_csp3(&topology), Ok(1.0));
    }

    fn d01_prepared_input() -> (
        TopologyBlock,
        cosmolkit_core::ValenceAssignment,
        cosmolkit_core::RingInfo,
    ) {
        let topology = detached_ethanol();
        let assignment = valence(&topology, "d01_fixture").unwrap();
        let ring_info = cosmolkit_core::find_sssr(
            &topology,
            &cosmolkit_core::RingSearchParams {
                include_dative_bonds: false,
                include_hydrogen_bonds: false,
            },
        )
        .unwrap();
        (topology, assignment, ring_info)
    }

    #[test]
    fn descriptor_d01_input_borrows_all_five_final_inputs_and_is_copy() {
        let (topology, assignment, ring_info) = d01_prepared_input();
        let coordinates = cosmolkit_model::CoordinateBlock::default();
        let properties = cosmolkit_model::MoleculeProperties::default();
        let input = DescriptorInput::new(
            &topology,
            &coordinates,
            &properties,
            &assignment,
            &ring_info,
        );
        assert!(std::ptr::eq(input.topology(), &topology));
        assert!(std::ptr::eq(input.coordinates(), &coordinates));
        assert!(std::ptr::eq(input.properties(), &properties));
        assert!(std::ptr::eq(input.valence(), &assignment));
        assert!(std::ptr::eq(input.ring_info(), &ring_info));
        let copied = input;
        assert!(std::ptr::eq(copied.valence(), &assignment));
    }

    #[test]
    fn descriptor_d01_error_causes_preserve_categories_and_error_source() {
        use std::error::Error;

        // Topology cause: non-sequential atom rows fail validation and reach
        // the descriptor as the existing Unsupported category with detail.
        let bad_topology = TopologyBlock {
            atoms: vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(5), AtomSpec::new(Element::C)),
            ],
            ..TopologyBlock::default()
        };
        match lipinski_hba(&bad_topology) {
            Err(DescriptorError::Unsupported { function, detail }) => {
                assert_eq!(function, "lipinski_hba");
                assert!(detail.contains("has id"));
            }
            other => panic!("expected topology cause, got {other:?}"),
        }

        // Valence cause: a negative implicit row is a typed Valence error
        // with a working Error::source chain, never a silent zero.
        let single = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        let bad_assignment = cosmolkit_core::ValenceAssignment {
            explicit_valence: vec![0],
            implicit_hydrogens: vec![-1],
        };
        let error =
            molecular_weight_with_valence(&single, false, Some(&bad_assignment)).unwrap_err();
        assert!(matches!(
            &error,
            DescriptorError::Valence {
                source: cosmolkit_core::ValenceError::ImplicitValenceCacheNotInitialized { atom },
                ..
            } if atom.index() == 0
        ));
        assert!(
            error
                .source()
                .unwrap()
                .downcast_ref::<cosmolkit_core::ValenceError>()
                .is_some()
        );

        // Periodic-table cause: an isotope absent from the pinned table is
        // the existing Unsupported category carrying the table's own error.
        let exotic = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::H).with_isotope(99),
            )],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        match molecular_weight(&exotic) {
            Err(DescriptorError::Unsupported { function, detail }) => {
                assert_eq!(function, "molecular_weight");
                assert!(detail.contains("unknown isotope 99"));
            }
            other => panic!("expected periodic-table cause, got {other:?}"),
        }

        // Ring cause: a real RingFindingError from the canonical core ring
        // owner (RingInfo::isRingFused PRECONDITION on an out-of-bounds ring
        // index over a zero-ring prepared input) keeps its category and
        // Error::source chain.
        let (_, _, ring_info) = d01_prepared_input();
        let mut ring_rows = ring_info;
        let ring_error = ring_rows.is_ring_fused(0).unwrap_err();
        assert!(matches!(
            &ring_error,
            cosmolkit_core::RingFindingError::Value { message }
                if *message == "ringIdx out of bounds"
        ));
        let error = DescriptorError::Ring {
            function: "descriptor_d01_fixture",
            source: ring_error,
        };
        assert!(error.to_string().contains("ring source error"));
        assert!(
            error
                .source()
                .unwrap()
                .downcast_ref::<cosmolkit_core::RingFindingError>()
                .is_some()
        );

        // Search cause: a real SMARTS parse failure flows through the wired
        // search dependency into the typed cause with a two-level chain.
        let parse_error =
            cosmolkit_search::parse_smarts("((", &cosmolkit_search::SmartsParseParams::default())
                .unwrap_err();
        let error = DescriptorError::Search {
            function: "descriptor_d01_fixture",
            source: DescriptorSearchCause::Compile(parse_error),
        };
        assert!(error.to_string().contains("search source error"));
        let level_one = error.source().unwrap();
        assert!(level_one.downcast_ref::<DescriptorSearchCause>().is_some());
        assert!(
            level_one
                .source()
                .unwrap()
                .downcast_ref::<cosmolkit_search::SmartsParseError>()
                .is_some()
        );

        // Stereo cause: a real validation failure value raised through the
        // legacy stereo owner's transparent cause keeps its chain.
        let stereo_source = bad_topology.validate().unwrap_err();
        let error = DescriptorError::Stereo {
            function: "descriptor_d01_fixture",
            source: LegacyStereoError::InvalidTopology(stereo_source),
        };
        assert!(error.to_string().contains("stereo source error"));
        let level_one = error.source().unwrap();
        assert!(level_one.downcast_ref::<LegacyStereoError>().is_some());
        // The source variant is `#[error(transparent)]`, so Display and the
        // remaining source chain forward to the wrapped topology error; the
        // AtomIdMismatch leaf itself carries no deeper source.
        assert!(level_one.to_string().contains("has id"));
        assert!(level_one.source().is_none());
    }

    #[test]
    fn descriptor_d01_malformed_valence_shape_is_typed_not_zero_not_unsupported() {
        let (topology, assignment, _) = d01_prepared_input();
        for field_explicit in [false, true] {
            for length in [0_usize, 2, 4] {
                let mut malformed = assignment.clone();
                if field_explicit {
                    malformed.explicit_valence.resize(length, 0);
                } else {
                    malformed.implicit_hydrogens.resize(length, 0);
                }
                for heavy in [false, true] {
                    let outcome = molecular_weight_with_valence(&topology, heavy, Some(&malformed));
                    if heavy {
                        // Source getAvgMolWt(onlyHeavy=true) never reads H
                        // state, so a malformed H-bearing shape is not even
                        // observed by the heavy-only kernel.
                        assert!(outcome.is_ok(), "heavy={heavy} len={length}");
                    } else {
                        let error = outcome.unwrap_err();
                        assert_eq!(
                            error,
                            DescriptorError::InvalidValenceRows {
                                function: "molecular_weight",
                                field: if field_explicit {
                                    "explicit_valence"
                                } else {
                                    "implicit_hydrogens"
                                },
                                actual: length,
                                expected: topology.atoms.len(),
                            },
                            "field_explicit={field_explicit} len={length}"
                        );
                        assert!(!matches!(error, DescriptorError::Unsupported { .. }));
                    }
                }
                // The H-sensitive exact-mass kernel rejects the same shapes.
                assert!(
                    exact_molecular_weight_with_valence(&topology, false, Some(&malformed))
                        .is_err()
                );
            }
        }
    }

    #[test]
    fn descriptor_d01_borrowed_kernels_reused_through_descriptor_input() {
        let (topology, assignment, ring_info) = d01_prepared_input();
        let coordinates = cosmolkit_model::CoordinateBlock::default();
        let properties = cosmolkit_model::MoleculeProperties::default();
        let input = DescriptorInput::new(
            &topology,
            &coordinates,
            &properties,
            &assignment,
            &ring_info,
        );
        for heavy in [false, true] {
            assert_eq!(
                molecular_weight_with_valence(&input.topology(), heavy, Some(input.valence())),
                molecular_weight_with_options(&topology, heavy)
            );
            assert_eq!(
                exact_molecular_weight_with_valence(
                    &input.topology(),
                    heavy,
                    Some(input.valence())
                ),
                exact_molecular_weight_with_options(&topology, heavy)
            );
        }
        for separate in [false, true] {
            for abbreviate in [false, true] {
                assert_eq!(
                    molecular_formula_with_valence(
                        &input.topology(),
                        separate,
                        abbreviate,
                        Some(input.valence())
                    ),
                    molecular_formula_with_options(&topology, separate, abbreviate)
                );
            }
        }
    }

    fn prepared_input_for(
        topology: &TopologyBlock,
    ) -> (
        cosmolkit_core::ValenceAssignment,
        cosmolkit_core::RingInfo,
        cosmolkit_model::CoordinateBlock,
        cosmolkit_model::MoleculeProperties,
    ) {
        let assignment = valence(topology, "prepared_fixture").unwrap();
        let ring_info = cosmolkit_core::find_sssr(
            topology,
            &cosmolkit_core::RingSearchParams {
                include_dative_bonds: false,
                include_hydrogen_bonds: false,
            },
        )
        .unwrap();
        (
            assignment,
            ring_info,
            cosmolkit_model::CoordinateBlock::default(),
            cosmolkit_model::MoleculeProperties::default(),
        )
    }

    #[test]
    fn descriptor_d03_direct_n_o_count_neutral_charged_aromatic_empty() {
        // Source Lipinski.cpp:77-86: the ONLY predicate is atomic number 7
        // or 8 — charge, aromaticity, hydrogens, and degree are invisible.
        // This is deliberately NOT the recursive general-HBA SMARTS: a
        // quaternary [N+] counts 1 here while the general pattern (Q03)
        // excludes it.
        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("CCO", 1),
            ("N#N", 2),
            ("[NH4+]", 1),
            ("[N+](C)(C)(C)C", 1),
            ("c1ccccc1", 0),
            ("c1cc[nH]c1", 1),
            ("O=C=O", 2),
            ("CC(=O)N", 2),
            ("[Na+].[Cl-]", 0),
            ("CS", 0),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(lipinski_hba(&topology).unwrap(), *expected, "{smiles}");
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                lipinski_hba_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_d04_attached_h_sum_explicit_implicit_noimplicit_isotope() {
        // Source Lipinski.cpp:88-97: N/O rows sum getTotalNumHs(true) —
        // explicit H + implicit H + hydrogen ATOM neighbors (the neighbor
        // predicate is atomicNum==1, so isotopic [2H] neighbors count).
        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("CCO", 1),
            ("O", 2),
            ("NCC(=O)O", 3),
            ("[NH4+]", 4),
            ("c1cc[nH]c1", 1),
            ("[2H]O", 2),
            ("[N+](C)(C)(C)C", 0),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(lipinski_hbd(&topology).unwrap(), *expected, "{smiles}");
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                lipinski_hbd_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // noImplicit N with two explicit H ATOM neighbors: neighbors count,
        // implicit rows contribute zero by policy, never by defaulting.
        let mut atoms = vec![
            Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::N).with_no_implicit(true),
            ),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        atoms.extend(
            (2..4).map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::H))),
        );
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(2),
                BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
            ),
        ];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
        let input = DescriptorInput::new(
            &topology,
            &coordinates,
            &properties,
            &assignment,
            &ring_info,
        );
        assert_eq!(lipinski_hbd(&topology).unwrap(), 2);
        assert_eq!(lipinski_hbd_prepared(&input).unwrap(), 2);
    }

    #[test]
    fn descriptor_d04_supplied_valence_shapes_reject_missing_rows() {
        let topology = smiles_topology("NCC(=O)O");
        let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
        for field_explicit in [false, true] {
            for length in [0_usize, 2, 4] {
                let mut malformed = assignment.clone();
                if field_explicit {
                    malformed.explicit_valence.resize(length, 0);
                } else {
                    malformed.implicit_hydrogens.resize(length, 0);
                }
                let input = DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &malformed,
                    &ring_info,
                );
                let error = lipinski_hbd_prepared(&input).unwrap_err();
                assert_eq!(
                    error,
                    DescriptorError::InvalidValenceRows {
                        function: "lipinski_hbd",
                        field: if field_explicit {
                            "explicit_valence"
                        } else {
                            "implicit_hydrogens"
                        },
                        actual: length,
                        expected: topology.atoms.len(),
                    },
                    "field_explicit={field_explicit} len={length}"
                );
                assert!(!matches!(error, DescriptorError::Unsupported { .. }));
            }
        }
    }

    #[test]
    fn descriptor_d05_heavy_atom_predicate_empty_c_h_dummy_isotopic_h_mixed() {
        // Source ROMol.cpp:187-196: one pass over atom rows, predicate
        // `atomicNum() > 1` — isotope/charge/aromatic state never consulted.
        const CASES: &[(&str, u32)] =
            &[("", 0), ("C", 1), ("CCO", 3), ("[H][H]", 0), ("[2H]C#N", 2)];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(num_heavy_atoms(&topology).unwrap(), *expected, "{smiles}");
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                num_heavy_atoms_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Manual mixed rows: C counts; dummy (atomicNum 0) and isotopic [2H]
        // (atomicNum 1) do not — the direct source predicate on atomic number.
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::DUMMY)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::H).with_isotope(2)),
        ];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        assert_eq!(num_heavy_atoms(&topology).unwrap(), 1);
        let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
        let input = DescriptorInput::new(
            &topology,
            &coordinates,
            &properties,
            &assignment,
            &ring_info,
        );
        assert_eq!(num_heavy_atoms_prepared(&input).unwrap(), 1);
    }

    #[test]
    fn descriptor_d06_rows_plus_attached_h_explicit_implicit_rows_only() {
        // Source ROMol.cpp:176-186 (!onlyExplicit): rows + sum of
        // getTotalNumHs() with DEFAULT includeNeighbors=false — explicit
        // H-atom rows count once as rows, never again as neighbor state.
        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("C", 5),
            ("CCO", 9),
            ("[NH4+]", 5),
            ("[H][H]", 2),
            ("c1cc[nH]c1", 10),
            ("[2H]O", 3),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(num_atoms(&topology).unwrap(), *expected, "{smiles}");
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                num_atoms_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_d06_supplied_valence_shapes_reject_missing_rows() {
        let topology = smiles_topology("CCO");
        let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
        for field_explicit in [false, true] {
            for length in [0_usize, 1, 2] {
                let mut malformed = assignment.clone();
                if field_explicit {
                    malformed.explicit_valence.resize(length, 0);
                } else {
                    malformed.implicit_hydrogens.resize(length, 0);
                }
                let input = DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &malformed,
                    &ring_info,
                );
                let error = num_atoms_prepared(&input).unwrap_err();
                assert_eq!(
                    error,
                    DescriptorError::InvalidValenceRows {
                        function: "num_atoms",
                        field: if field_explicit {
                            "explicit_valence"
                        } else {
                            "implicit_hydrogens"
                        },
                        actual: length,
                        expected: topology.atoms.len(),
                    },
                    "field_explicit={field_explicit} len={length}"
                );
            }
        }
    }

    #[test]
    fn descriptor_d07_total_degree_criterion_carbon_classes() {
        // Source Lipinski.cpp:212-232: carbon rows only; sp3 criterion is
        // getTotalDegree()==4 (bond neighbors + implicit/explicit-property H),
        // never a hybridization flag; source-defined zero-carbon fallback 0.
        const CASES: &[(&str, f64)] = &[
            ("", 0.0),
            ("CCO", 1.0),
            ("c1ccccc1", 0.0),
            ("[Na+].[Cl-]", 0.0),
            ("CC#CC", 0.5),
            ("C[C+](C)C", 0.75),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(fraction_csp3(&topology).unwrap(), *expected, "{smiles}");
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                fraction_csp3_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_d07_supplied_valence_shapes_reject_missing_rows() {
        let topology = smiles_topology("CCO");
        let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
        for length in [0_usize, 1, 2] {
            let mut malformed = assignment.clone();
            malformed.implicit_hydrogens.resize(length, 0);
            let input =
                DescriptorInput::new(&topology, &coordinates, &properties, &malformed, &ring_info);
            assert_eq!(
                fraction_csp3_prepared(&input).unwrap_err(),
                DescriptorError::InvalidValenceRows {
                    function: "fraction_csp3",
                    field: "implicit_hydrogens",
                    actual: length,
                    expected: topology.atoms.len(),
                }
            );
        }
    }

    #[test]
    fn descriptor_q01_pattern_compiles_once_and_is_retained() {
        // Source Lipinski.cpp:56-58: pattern_flyweight with no_tracking
        // parses each distinct pattern once per process and shares the
        // matcher forever after. Two acquisitions of one pattern must share
        // one compiled query; a distinct pattern compiles separately.
        // A/B use patterns UNIQUE to this test: the compile counter is
        // thread-local (patterns.rs) while the retained map is
        // process-global, so a pattern shared with a concurrent sibling
        // test could pre-populate the map and hide the compile branch on
        // this thread (observed once as left:0/right:1); a unique pattern
        // keeps the exactly-once assertion deterministic.
        const A: &str = "[#15]";
        const B: &str = "[#16]";
        let before = crate::patterns::PATTERN_COMPILES.with(|count| count.get());
        let first = crate::patterns::retained_pattern("q01_retention", A).unwrap();
        let second = crate::patterns::retained_pattern("q01_retention", A).unwrap();
        let after = crate::patterns::PATTERN_COMPILES.with(|count| count.get());
        assert_eq!(after - before, 1, "same pattern must compile exactly once");
        assert!(std::sync::Arc::ptr_eq(&first, &second));
        let other = crate::patterns::retained_pattern("q01_retention", B).unwrap();
        assert!(!std::sync::Arc::ptr_eq(&first, &other));
    }

    #[test]
    fn descriptor_q01_shared_prepared_context_counts_plain_and_recursive_patterns() {
        // Source Lipinski.cpp:39-53 countMatches with default parameters:
        // plain and '$'-recursive patterns are counted against one target
        // with the same default SubstructMatchParameters (uniquify=true,
        // maxMatches=1000, recursionPossible=true). One prepared
        // QueryMatchContext is shared by every pattern.
        use cosmolkit_search::build_prepared_query_match_context;
        const PLAIN: &str = "[#6]";
        const RECURSIVE: &str = "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]";
        for (smiles, plain, recursive) in [
            ("", 0, 0),
            ("CC#CC", 4, 0),
            // biphenyl: c1ccccc1-c1ccccc1
            ("c1ccccc1-c1ccccc1", 12, 1),
        ] {
            let topology = smiles_topology(smiles);
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            let context = build_prepared_query_match_context(
                input.topology(),
                input.ring_info(),
                input.valence(),
            )
            .unwrap();
            assert_eq!(
                crate::patterns::count_pattern_matches_with_context(
                    &input,
                    "q01_shared_context",
                    PLAIN,
                    &context
                )
                .unwrap(),
                plain,
                "plain pattern on {smiles}"
            );
            assert_eq!(
                crate::patterns::count_pattern_matches_with_context(
                    &input,
                    "q01_shared_context",
                    RECURSIVE,
                    &context
                )
                .unwrap(),
                recursive,
                "recursive pattern on {smiles}"
            );
            assert_eq!(
                crate::patterns::count_pattern_matches(&input, "q01_shared_context", PLAIN)
                    .unwrap(),
                plain,
                "single-pattern convenience form on {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_q01_compile_failure_is_typed_error_not_default() {
        // Source Lipinski.cpp:35-37: a failed SmartsToMol hits
        // POSTCONDITION(m_matcher, "no matcher"), i.e. a failed fixed
        // pattern is a fatal programming error, never a silent zero count.
        // COSMolKit surfaces the typed compile cause on every attempt and
        // retains no query for the failure.
        for _ in 0..2 {
            match crate::patterns::retained_pattern("q01_compile_failure", "((") {
                Err(DescriptorError::Search {
                    source: DescriptorSearchCause::Compile(_),
                    ..
                }) => {}
                other => panic!("expected typed compile error, got {other:?}"),
            }
        }
        let topology = smiles_topology("CCO");
        let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
        let input = DescriptorInput::new(
            &topology,
            &coordinates,
            &properties,
            &assignment,
            &ring_info,
        );
        match crate::patterns::count_pattern_matches(&input, "q01_compile_failure", "((") {
            Err(DescriptorError::Search {
                function: "q01_compile_failure",
                source: DescriptorSearchCause::Compile(_),
            }) => {}
            other => panic!("expected typed compile error, got {other:?}"),
        }
    }

    #[test]
    fn descriptor_q02_num_hbd_source_pattern_branches() {
        // Source Lipinski.cpp:197-198: SMARTSCOUNTFUNC(NumHBD,
        // "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]", "2.0.1") —
        // one single-atom query; count = matching atoms. Branches below
        // cover positive (N-v3, O-H1, S-H1), negative (benzene, alkyne),
        // charged ([NH4+] via N&!H0&+1&v4), aromatic ([nH]) and explicit-H
        // folding (H primitive = total H: implicit + explicit neighbors).
        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("CCO", 1),
            ("NCC(=O)O", 2),
            ("[NH4+]", 1),
            ("c1cc[nH]c1", 1),
            ("c1ccccc1", 0),
            ("CC#CC", 0),
            ("[H]N([H])[H]", 1),
            ("CS", 1),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(num_hbd(&topology).unwrap(), *expected, "cold {smiles}");
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                lipinski::num_hbd_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
        assert_eq!(NUM_HBD_VERSION, "2.0.1");
    }

    #[test]
    fn descriptor_q02_num_hbd_distinct_from_lipinski_donor_sum() {
        // The general NumHBD pattern count and the DIRECT Lipinski donor
        // hydrogen sum (Lipinski.cpp:88-97) are different source functions:
        // sulfur matches only the pattern; [NH4+] counts one pattern atom
        // but four attached hydrogens in the direct sum.
        let methanethiol = smiles_topology("CS");
        assert_eq!(num_hbd(&methanethiol).unwrap(), 1);
        assert_eq!(lipinski_hbd(&methanethiol).unwrap(), 0);
        let ammonium = smiles_topology("[NH4+]");
        assert_eq!(num_hbd(&ammonium).unwrap(), 1);
        assert_eq!(lipinski_hbd(&ammonium).unwrap(), 4);
    }

    #[test]
    fn descriptor_q03_num_hba_recursive_pattern_branches() {
        // Source Lipinski.cpp:199-202: SMARTSCOUNTFUNC(NumHBA, "[$([O,S;H1;
        // v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=
        // !@[O,N,P,S])]),$([nH0X2,o,s;+0])]", "2.0.2") — one single-atom
        // query, comma-OR of five recursive subqueries. Fixed expectations
        // derived from the pattern alone (v = total valence, so a carbonyl
        // O is H0;v2 and matches branch 2; the last two cases mirror RDKit
        // test.cpp:2227-2234 exactly).
        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("CCO", 1),           // b1: OH on plain C
            ("COC", 1),           // b2: ether O
            ("[OH2]", 0),         // H2v2 matches neither H1v2 nor H0v2
            ("CC(=O)O", 1),       // acid OH excluded by b1 neighbor; carbonyl O via b2
            ("CC(=O)[O-]", 2),    // b2 carbonyl + b3 anionic O
            ("CC(=O)N", 1),       // b4 excludes amide N; carbonyl O via b2
            ("CC#N", 1),          // b4: nitrile N is v3 with no single-bond path
            ("CN", 1),            // b4: plain amine N
            ("Nc1ccccc1", 1),     // b4: aromatic-attached amine N
            ("c1ccncc1", 1),      // b5: pyridine nH0X2
            ("c1ccoc1", 1),       // b5: furan o
            ("c1cccc(=O)n1C", 1), // RDKit test.cpp:2232 mirror (lactam carbonyl O)
            ("c1cccn1C", 0),      // RDKit test.cpp:2227 mirror (substituted n fails X2)
            ("CC(=O)OCC", 2),     // b2 carbonyl + b2 ester alkoxy O
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(num_hba(&topology).unwrap(), *expected, "cold {smiles}");
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                lipinski::num_hba_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
        assert_eq!(NUM_HBA_VERSION, "2.0.2");
    }

    #[test]
    fn descriptor_q03_num_hba_distinct_from_lipinski_acceptor_sum() {
        // The general NumHBA recursive pattern and the DIRECT Lipinski N+O
        // sum (Lipinski.cpp:77-86) are different source functions: an acid
        // OH counts in neither but its carbonyl O counts once here and both
        // oxygens count in the direct sum; thiophene sulfur counts only in
        // the pattern.
        let acid = smiles_topology("CC(=O)O");
        assert_eq!(num_hba(&acid).unwrap(), 1);
        assert_eq!(lipinski_hba(&acid).unwrap(), 2);
        let thiophene = smiles_topology("c1ccsc1");
        assert_eq!(num_hba(&thiophene).unwrap(), 1);
        assert_eq!(lipinski_hba(&thiophene).unwrap(), 0);
    }

    #[test]
    fn descriptor_q04_num_heteroatoms_literal_pattern_and_rows() {
        // Source Lipinski.cpp:203: SMARTSCOUNTFUNC(NumHeteroatoms,
        // "[!#6;!#1]", "1.0.1") — one single-atom query, `;`-AND of two
        // negated atomic-number primitives. Expectations derive from the
        // pattern alone: graph rows with atomicNum not in {1,6}; implicit H
        // never counts; [2H] is #1 (excluded); dummy (0) satisfies both
        // negations and counts.
        assert_eq!(NUM_HETEROATOMS_PATTERN, "[!#6;!#1]");
        assert_eq!(NUM_HETEROATOMS_VERSION, "1.0.1");
        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("CCO", 1),
            ("c1ccccc1", 0),
            ("[H][H]", 0),
            ("[2H]O", 1),
            ("CCl", 1),
            ("[Na+].[Cl-]", 2),
            ("c1cc[nH]c1", 1),
            ("CC#CC", 0),
            ("[NH4+]", 1),
            ("[OH2]", 1),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_heteroatoms(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                lipinski::num_heteroatoms_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Manual rows: C (#6) and [2H] (#1) excluded; DUMMY (atomicNum 0)
        // satisfies `!#6;!#1` and counts — pure pattern consequence.
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::DUMMY)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::H).with_isotope(2)),
        ];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        assert_eq!(num_heteroatoms(&topology).unwrap(), 1);
        let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
        let input = DescriptorInput::new(
            &topology,
            &coordinates,
            &properties,
            &assignment,
            &ring_info,
        );
        assert_eq!(lipinski::num_heteroatoms_prepared(&input).unwrap(), 1);
    }

    #[test]
    fn descriptor_q05_num_amide_bonds_query_count_and_contrasts() {
        // Source Lipinski.cpp:204: SMARTSCOUNTFUNC(NumAmideBonds,
        // "C(=[O;!R])N", "1.0.0") — three-atom query, non-recursive,
        // default SubstructMatchParameters with uniquify: the count is the
        // number of UNIQUE full match vectors. Urea O=C(N)N yields TWO
        // (the two N mappings are distinct vectors); a lactam counts
        // because its carbonyl O is acyclic (!R) while ring C/N are
        // unconstrained; sulfonamides do not match (carbonyl atom must
        // be carbon).
        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("CC(=O)N", 1),
            ("CC(=O)NC(C)C", 1),
            ("NC(=O)N", 2),
            ("CC(=O)O", 0),
            ("CC(=O)OC", 0),
            ("CC(=O)C", 0),
            ("O=C1NCCC1", 1),
            ("c1ccccc1C(=O)N", 1),
            ("CC(=O)NCC(=O)N", 2),
            ("CC#N", 0),
            ("CS(=O)(=O)N", 0),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_amide_bonds(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                lipinski::num_amide_bonds_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
        assert_eq!(NUM_AMIDE_BONDS_VERSION, "1.0.0");
    }

    #[test]
    fn descriptor_r01_option_enum_four_values_and_pinned_default_resolution() {
        // Source Lipinski.h:39-43: enum NumRotatableBondsOptions { Default,
        // NonStrict, Strict, StrictLinkages }. The four variants are the
        // complete parameter product of this unit; the deprecated bool
        // overload (Lipinski.cpp:185-187, true->Strict / false->NonStrict)
        // is documented source-only and intentionally has no Rust alias.
        let default = RotatableBondsOptions::Default;
        let non_strict = RotatableBondsOptions::NonStrict;
        let strict = RotatableBondsOptions::Strict;
        let strict_linkages = RotatableBondsOptions::StrictLinkages;

        // Four distinct values; Copy + equality round-trip.
        let all = [default, non_strict, strict, strict_linkages];
        for (i, a) in all.iter().enumerate() {
            for (j, b) in all.iter().enumerate() {
                assert_eq!(a == b, i == j);
            }
            let copied = *a;
            assert_eq!(copied, *a);
        }

        // Source default argument is the Default option (Lipinski.h:57-58).
        assert_eq!(RotatableBondsOptions::default(), default);

        // Pinned build default: RDK_USE_STRICT_ROTOR_DEFINITION is ON
        // (third_party/rdkit/CMakeLists.txt:55,586-588), so the anonymous
        // namespace DefaultStrictDefinition is Strict (Lipinski.cpp:99-106)
        // and the dispatch resolves Default -> Strict (Lipinski.cpp:113-115);
        // the other three options pass through unchanged.
        assert_eq!(
            rotatable::resolve_options(default),
            strict,
            "pinned Default must resolve to Strict"
        );
        assert_eq!(rotatable::resolve_options(non_strict), non_strict);
        assert_eq!(rotatable::resolve_options(strict), strict);
        assert_eq!(rotatable::resolve_options(strict_linkages), strict_linkages);
        // Resolution never invents a fifth state: Default is the only
        // option that changes.
        assert_ne!(rotatable::resolve_options(default), non_strict);
        assert_ne!(rotatable::resolve_options(default), strict_linkages);

        // Version literal (Lipinski.cpp:111).
        assert_eq!(NUM_ROTATABLE_BONDS_VERSION, "3.2.0");
    }

    #[test]
    fn descriptor_r02_non_strict_rotatable_count_exclusions() {
        // Source Lipinski.cpp:116: NonStrict branch counts matches of
        // "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]" — each endpoint atom must not
        // be attached to any triple bond and must not be degree-1; the
        // bond must be (single or aromatic) and not in a ring.
        assert_eq!(
            NON_STRICT_ROTATABLE_PATTERN,
            "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]"
        );

        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("CC", 0),          // both atoms degree-1 (terminal)
            ("CCC", 0),         // every single bond touches a terminal CH3
            ("CCCC", 1),        // central bond: both endpoints D2
            ("CCOCC", 2),       // both mid bonds
            ("CC#CC", 0),       // triple-adjacent endpoints
            ("CCC#CCC", 0),     // both candidate bonds flank the triple
            ("C1CCCCC1", 0),    // every bond is a ring bond
            ("C1CCCCC1CCC", 2), // exocyclic rC-CH2 AND CH2-CH2 both count;
            // only the terminal CH3 bond is excluded
            ("c1ccccc1-c1ccccc1", 1), // chain single bond between two aromatic
            // ring atoms, both D3 non-triple
            ("CC(=O)NC", 1), // NonStrict COUNTS the amide C-N bond
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_rotatable_bonds(&topology, RotatableBondsOptions::NonStrict).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rotatable::num_rotatable_bonds_prepared(&input, RotatableBondsOptions::NonStrict)
                    .unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Staged surface: after R06 every option is live. StrictLinkages
        // on CCC: base=0 (every single bond touches a terminal CH3) ->
        // early return 0. (The loop originally covered Default/Strict
        // while R02 was newest, then only StrictLinkages while R03-R05
        // were newest; R06 replaced the staged-Unsupported assertion with
        // the live value — see descriptor_r06_ test.)
        let topology = smiles_topology("CCC");
        assert_eq!(
            num_rotatable_bonds(&topology, RotatableBondsOptions::StrictLinkages).unwrap(),
            0
        );
    }

    #[test]
    fn descriptor_r03_strict_rotatable_count_exclusions() {
        // Source Lipinski.cpp:119-128: Strict branch counts matches of the
        // assembled strict SMARTS (five C++ adjacent literals -> ONE
        // string). Both endpoints exclude triple-adjacent, degree-1,
        // CF3/CCl3/CBr3, tert-butyl, and methyl atoms; the FIRST endpoint
        // additionally excludes amide/ester/amidinium carbonyl carbons and
        // their N/O/S partners (both match directions are blocked for an
        // amide bond: the carbonyl-C direction and the N-partner direction).
        assert_eq!(
            STRICT_ROTATABLE_PATTERN,
            "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CH3])]"
        );

        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("CC", 0),                // terminal (degree-1) endpoints
            ("CCCC", 1),              // plain butane central bond still counts
            ("C1CCCCC1", 0),          // ring bonds
            ("c1ccccc1-c1ccccc1", 1), // Strict still counts the ring-ring
            // chain link (only StrictLinkages
            // removes it)
            ("CC(=O)NC", 0),   // amide C-N excluded (NonStrict: 1)
            ("CC(=O)OC", 0),   // ester C(=O)-O excluded (NonStrict: 1)
            ("CCC(F)(F)F", 0), // terminal CF3 carbon has heavy
            // degree 4; only the !$(C(F)(F)F)
            // arm kills the C1-C2 bond
            // (NonStrict: 1)
            ("CC(C)(C)C(C)C", 0), // quaternary C with three CH3 neighbors
                                  // hits the tBu arm (NonStrict: 1)
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_rotatable_bonds(&topology, RotatableBondsOptions::Strict).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rotatable::num_rotatable_bonds_prepared(&input, RotatableBondsOptions::Strict)
                    .unwrap(),
                *expected,
                "prepared {smiles}"
            );
            // Pinned build default: Default resolves to Strict
            // (RDK_USE_STRICT_ROTOR_DEFINITION ON), so the default-argument
            // behavior equals the explicit Strict count on every probe.
            assert_eq!(
                num_rotatable_bonds(&topology, RotatableBondsOptions::Default).unwrap(),
                *expected,
                "default-arg {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_r04_strict_linkages_base_and_symmetric_rings() {
        // Source Lipinski.cpp:144,146-149,154-163: StrictLinkages staged
        // arithmetic (R04 scope): base count of "[!$([D1&!#1])]-,:;!@[!$([D1&!#1])]"
        // — endpoint excluded only when degree-1 AND not hydrogen (an
        // explicit-H endpoint does NOT block) — early return 0 without the
        // symmetric-ring evaluation when the base is empty, then ordered
        // subtraction of the symmetric aromatic-6-ring link count clamped
        // at zero. This staged function is the real caller stage the final
        // StrictLinkages branch will invoke; probes avoid R05/R06-only
        // behavior (no terminal triple bonds adjacent to candidate bonds,
        // no non-ring amides).
        assert_eq!(
            STRICT_LINKAGES_BASE_PATTERN,
            "[!$([D1&!#1])]-,:;!@[!$([D1&!#1])]"
        );
        assert_eq!(
            SYMMETRIC_RINGS_PATTERN,
            "[a;r6;$(a(-,:;!@[a;r6])(a[!#1])a[!#1])]-,:;!@[a;r6;$(a(-,:;!@[a;r6])(a[!#1])a)]"
        );

        const CASES: &[(&str, i32)] = &[
            ("", 0),                          // base=0 -> early return
            ("CC", 0),                        // both endpoints D1 heavy -> base=0
            ("CCCC", 1),                      // base=1, symRings=0
            ("c1ccccc1-c1ccccc1", 0),         // base=1, symRings=1 (biaryl link)
            ("CC1CCCC(C)C1-C1C(C)CCCC1C", 1), // aliphatic ring-C link stays
            // rotatable (source comment
            // example); methyl bonds D1,
            // ring bonds !@-excluded
            ("Cc1cccc(C)c1-c1c(C)cccc1", 0), // base=1, symRings=1 (source
            // comment: non-rotatable)
            // Explicit-H nuance of the base pattern's `!#1`: hydrogen
            // endpoints are D1 but #1, so the negated [D1&!#1] does NOT
            // exclude them — all six C-H bonds plus the C-C bond match
            // (7 distinct vectors; each H is a distinct atom).
            ("[H]C([H])([H])C([H])([H])[H]", 7),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rotatable::strict_linkages_stage(&input).unwrap(),
                *expected,
                "staged {smiles}"
            );
        }

        // Staging discipline: after R06 all options are live.
        // StrictLinkages on CCCC: base=1, no symmetric rings, triples, or
        // amides -> 1 (source-evidenced by the R06 liveness; the earlier
        // staged-Unsupported assertion was replaced).
        let topology = smiles_topology("CCCC");
        assert_eq!(
            num_rotatable_bonds(&topology, RotatableBondsOptions::StrictLinkages).unwrap(),
            1
        );
    }

    #[test]
    fn descriptor_r05_strict_linkages_triple_bond_correction() {
        // Source Lipinski.cpp:150,162-167: after the symmetric-ring stage,
        // subtract the count of "C#[#6,#7]" — any C#C or C#N triple bond
        // (NO degree constraint despite the variable name; the StrictLinkages
        // base pattern does not exclude triple-adjacent atoms, so single
        // chain bonds flanking a triple must be removed here) — then clamp
        // at zero. Composes the R04 stage output.
        assert_eq!(TERMINAL_TRIPLE_BONDS_PATTERN, "C#[#6,#7]");

        const CASES: &[(&str, i32)] = &[
            ("CCCC", 1),   // no triple: unchanged from R04 stage
            ("CCC#CC", 0), // base=1 (propargylic C1-C2),
            // symRings=0, triple=1 -> 0
            // (R04 stage value was 1: composition
            // visible)
            ("CCCC#CC", 1), // base=2 (both propargylic bonds),
            // triple=1 -> 1
            ("CCC#N", 0), // base=1 (C1-C2), triple=1 (C#N) -> 0
            ("CC#CC", 0), // base=0 -> early return; the triple
            // stage is never reached
            ("c1ccccc1-c1ccccc1", 0), // base=1, symRings=1, triple=0 -> 0
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rotatable::strict_linkages_stage(&input).unwrap(),
                *expected,
                "staged {smiles}"
            );
        }

        // Staging discipline: after R06 all options are live.
        // StrictLinkages on CCC#CC: base=1, triple=1 -> 0 (R06 liveness
        // replaced the staged-Unsupported assertion).
        let topology = smiles_topology("CCC#CC");
        assert_eq!(
            num_rotatable_bonds(&topology, RotatableBondsOptions::StrictLinkages).unwrap(),
            0
        );
    }

    #[test]
    fn descriptor_r06_strict_linkages_shared_atom_amide_correction() {
        // Source Lipinski.cpp:145,168-189: ordered uniquified matches of
        // "[C&!R](=O)NC" walked with an atoms-seen bitset — a match
        // decrements only when NONE of its atoms was marked by an earlier
        // match AND res > 0; every atom of every match is marked in
        // matcher order regardless (no break on overlap or zero, no
        // reorder/dedup, no fragment counting, no match-count subtraction).
        assert_eq!(NON_RING_AMIDES_PATTERN, "[C&!R](=O)NC");

        const CASES: &[(&str, u32)] = &[
            ("CC(=O)NC", 0), // zero saturation: base=1, one distinct
            // amide decrement floors the result at 0
            ("CC(=O)NCC(=O)NC", 2), // disjoint amides: base=4, both
            // matches atom-disjoint -> 2 decrements
            ("CC(=O)N(C)C(=O)C", 1), // shared-N amides: FOUR matches
            // ((C1,O2,N3,C4),(C1,O2,N3,C5),
            // (C5,O6,N3,C1),(C5,O6,N3,C4)) but
            // overlap propagation lets ONLY the
            // first decrement: base=2 -> 1 (naive
            // match-count subtraction would clamp
            // to 0; distinct-fragment counting
            // would also give 0; source = 1)
            ("CCCC", 1),              // no amides: unchanged
            ("CCC#CC", 0),            // triple correction unchanged
            ("c1ccccc1-c1ccccc1", 0), // symmetric rings unchanged
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_rotatable_bonds(&topology, RotatableBondsOptions::StrictLinkages).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rotatable::num_rotatable_bonds_prepared(
                    &input,
                    RotatableBondsOptions::StrictLinkages
                )
                .unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Packed-word boundary: a 65-atom chain (61 chain carbons +
        // carbonyl C + O + N + methyl C) whose amide atoms sit at indices
        // 61..64, spanning the 64-bit word-0/word-1 boundary of the
        // atoms-seen bitset. base = 60 chain bonds (C1-C2 .. C60-C61) +
        // the C61-N63 bond = 61; one distinct amide decrement -> 60.
        let long_amide = format!("{}C(=O)NC", "C".repeat(61));
        let topology = smiles_topology(&long_amide);
        assert_eq!(topology.atoms.len(), 65);
        assert_eq!(
            num_rotatable_bonds(&topology, RotatableBondsOptions::StrictLinkages).unwrap(),
            60,
            "packed-word boundary amide"
        );

        // Test-only, descriptor-owned evidence: ONE prepared context build
        // per whole StrictLinkages evaluation (all four pattern passes
        // share it; the counter sees every descriptor-side build).
        let before = patterns::CONTEXT_BUILDS.with(|count| count.get());
        assert_eq!(
            num_rotatable_bonds(&topology, RotatableBondsOptions::StrictLinkages).unwrap(),
            60
        );
        let after = patterns::CONTEXT_BUILDS.with(|count| count.get());
        assert_eq!(
            after - before,
            1,
            "a full StrictLinkages evaluation must build exactly ONE context"
        );
    }

    #[test]
    fn descriptor_r07_complete_rotor_dispatch_matrix() {
        // Complete calcNumRotatableBonds dispatch (Lipinski.cpp:108-194):
        // all four modes against pinned RDKit test cases. Expectations are
        // the pattern-derived literals from the Step 148 audit:
        //  - CC(C)(C)c1cc(O)c(cc1O)C(C)(C)C is RDKit's own default-build
        //    sentinel (test.cpp:370-382): the C++ treats a count of 2 as
        //    PROOF of a NonStrict build, so the pinned Strict default must
        //    NOT return 2 there (derived: 0).
        //  - c1ccccc1c1ccc(CCC)cc1 and c1cc[nH]c1c1[nH]c(CCC)cc1 and CCCC
        //    are the catch_tests.cpp Github#5104 literals (NonStrict=3 /
        //    Strict=3 / StrictLinkages=2; StrictLinkages=3; Strict=1).
        const CASES: &[(&str, u32, u32, u32, u32)] = &[
            // (smiles, NonStrict, Strict, Default, StrictLinkages)
            ("CC(C)(C)c1cc(O)c(cc1O)C(C)(C)C", 2, 0, 0, 2),
            ("c1ccccc1c1ccc(CCC)cc1", 3, 3, 3, 2),
            ("c1cc[nH]c1c1[nH]c(CCC)cc1", 3, 3, 3, 3),
            ("CCCC", 1, 1, 1, 1),
        ];
        for (smiles, non_strict, strict, default, strict_linkages) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_rotatable_bonds(&topology, RotatableBondsOptions::NonStrict).unwrap(),
                *non_strict,
                "NonStrict {smiles}"
            );
            assert_eq!(
                num_rotatable_bonds(&topology, RotatableBondsOptions::Strict).unwrap(),
                *strict,
                "Strict {smiles}"
            );
            // Pinned build default: Default == Strict on every case
            // (RDK_USE_STRICT_ROTOR_DEFINITION ON), including the sentinel
            // probe where 2 would prove a NonStrict default.
            assert_eq!(
                num_rotatable_bonds(&topology, RotatableBondsOptions::Default).unwrap(),
                *default,
                "Default {smiles}"
            );
            assert_eq!(
                num_rotatable_bonds(&topology, RotatableBondsOptions::StrictLinkages).unwrap(),
                *strict_linkages,
                "StrictLinkages {smiles}"
            );
        }
        assert_eq!(NUM_ROTATABLE_BONDS_VERSION, "3.2.0");
    }

    #[test]
    fn descriptor_r07_detached_state_fixtures_partial_sanitization_and_explicit_h() {
        // Supervisor corrections R07-CLOSE + R07-STAGE: the pinned
        // catch_tests.cpp:194-239 states as FIXED DETACHED FIXTURES from
        // existing read-only owners. STAGE FACT (R07-STAGE Action 2):
        // cosmolkit-smiles::parse_smiles builds RAW PARSER ROWS and never
        // executes sanitize_topology — the smiles_topology helper returns
        // that raw topology. The blocks labeled RAW-PARSER below are
        // supplementary raw-row coverage only (no sanitizer-execution
        // claim); the PINNED A/B pristine/mutated fixtures run on the FULL
        // core sanitize owner output. Expectations frozen from source
        // asserts + pattern facts BEFORE execution (receipt R07-STAGE
        // Actions 1-2); every case runs all four modes on BOTH the cold
        // public path and the prepared route from the SAME fixture topology.
        let check = |topology: &TopologyBlock, label: &str, expected: [u32; 4]| {
            let modes = [
                RotatableBondsOptions::NonStrict,
                RotatableBondsOptions::Strict,
                RotatableBondsOptions::Default,
                RotatableBondsOptions::StrictLinkages,
            ];
            for (mode, expected) in modes.iter().zip(expected) {
                assert_eq!(
                    num_rotatable_bonds(topology, *mode).unwrap(),
                    expected,
                    "cold {label} {:?}",
                    mode
                );
            }
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(topology);
            let input =
                DescriptorInput::new(topology, &coordinates, &properties, &assignment, &ring_info);
            for (mode, expected) in modes.iter().zip(expected) {
                assert_eq!(
                    rotatable::num_rotatable_bonds_prepared(&input, *mode).unwrap(),
                    expected,
                    "prepared {label} {:?}",
                    mode
                );
            }
        };

        // (A-raw, supplementary) RAW-PARSER rows of c1ccccc1c1ccc(CCC)cc1:
        // the raw parser emits an Aromatic-order bridge(5,6); reversing the
        // order to Single is a raw-row kind flip. Counts stay 3/3/3/2.
        let biphenyl = smiles_topology("c1ccccc1c1ccc(CCC)cc1");
        let (a_id, b_id) = (AtomId::new(5), AtomId::new(6));
        let find_link = |topology: &TopologyBlock| {
            topology
                .bonds
                .iter()
                .find(|bond| {
                    (bond.begin(), bond.end()) == (a_id, b_id)
                        || (bond.begin(), bond.end()) == (b_id, a_id)
                })
                .expect("inter-ring bond (5,6) exists")
                .clone()
        };
        let pristine_link = find_link(&biphenyl);
        // Getter-free order assertions via set_order identity: setting
        // Aromatic is identity <=> the bond already carries Aromatic order.
        let mut aromatic_probe = pristine_link.clone();
        aromatic_probe.set_order(BondOrder::Aromatic);
        assert_eq!(
            aromatic_probe, pristine_link,
            "owner fact: bond 5-6 is Aromatic-order"
        );
        let mut expected_single = pristine_link.clone();
        expected_single.set_order(BondOrder::Single);
        let mut single_kind_biphenyl = biphenyl.clone();
        let mut mutated = false;
        for bond in &mut single_kind_biphenyl.bonds {
            if (bond.begin(), bond.end()) == (a_id, b_id)
                || (bond.begin(), bond.end()) == (b_id, a_id)
            {
                bond.set_order(BondOrder::Single);
                mutated = true;
            }
        }
        assert!(mutated);
        // State-existence: exactly the type-only mutation happened.
        assert_eq!(find_link(&single_kind_biphenyl), expected_single);
        assert_ne!(find_link(&single_kind_biphenyl), pristine_link);
        check(
            &single_kind_biphenyl,
            "biphenyl bond5-6 SINGLE kind",
            [3, 3, 3, 2],
        );
        // The aromatic-kind side (RDKit's mutated state) equals our base.
        check(
            &biphenyl,
            "biphenyl bond5-6 AROMATIC kind (base)",
            [3, 3, 3, 2],
        );

        // (B-raw, supplementary) RAW-PARSER rows of the r5-linked molecule,
        // same raw-row kind flip on bridge(4,5); counts stay 3/3/3/3.
        let r5_linked = smiles_topology("c1cc[nH]c1c1[nH]c(CCC)cc1");
        let (a_id, b_id) = (AtomId::new(4), AtomId::new(5));
        let find_r5_link = |topology: &TopologyBlock| {
            topology
                .bonds
                .iter()
                .find(|bond| {
                    (bond.begin(), bond.end()) == (a_id, b_id)
                        || (bond.begin(), bond.end()) == (b_id, a_id)
                })
                .expect("r5-link bond (4,5) exists")
                .clone()
        };
        let pristine_r5 = find_r5_link(&r5_linked);
        let mut r5_aromatic_probe = pristine_r5.clone();
        r5_aromatic_probe.set_order(BondOrder::Aromatic);
        assert_eq!(
            r5_aromatic_probe, pristine_r5,
            "owner fact: r5-link bond 4-5 is Aromatic-order"
        );
        let mut single_kind_r5 = r5_linked.clone();
        for bond in &mut single_kind_r5.bonds {
            if (bond.begin(), bond.end()) == (a_id, b_id)
                || (bond.begin(), bond.end()) == (b_id, a_id)
            {
                bond.set_order(BondOrder::Single);
            }
        }
        let mut expected_r5_single = pristine_r5.clone();
        expected_r5_single.set_order(BondOrder::Single);
        assert_eq!(find_r5_link(&single_kind_r5), expected_r5_single);
        assert_ne!(find_r5_link(&single_kind_r5), pristine_r5);
        check(&single_kind_r5, "r5-link bond4-5 SINGLE kind", [3, 3, 3, 3]);
        check(
            &r5_linked,
            "r5-link bond4-5 AROMATIC kind (base)",
            [3, 3, 3, 3],
        );

        // (C, supplementary) CCCC with all ten hydrogens as EXPLICIT RAW
        // PARSER ROWS (parse_smiles with remove_hydrogens:false; NO
        // sanitize_topology execution — params.sanitize does not run the
        // core sanitizer). Source analogue: MolOps::addHs output rows.
        // State existence: 14 topology atoms (4 C + 10 explicit H).
        // Expectations unchanged (frozen): Strict = 1 (SOURCE-asserted:
        // `!$([CH3])` fires on the terminal carbons because the H primitive
        // folds their three explicit-H neighbors), NonStrict = 3,
        // Default = 1, StrictLinkages = 13 (base `!$([D1&!#1])` does NOT
        // exclude degree-1 hydrogens: 10 C-H + 3 C-C).
        let mut keep_h = cosmolkit_smiles::SmilesParseParams::default();
        keep_h.remove_hydrogens = false;
        let explicit_h_butane = cosmolkit_smiles::parse_smiles(
            "[H]C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H]",
            &keep_h,
        )
        .expect("parse explicit-H butane fixture")
        .topology;
        assert_eq!(explicit_h_butane.atoms.len(), 14);
        check(&explicit_h_butane, "CCCC addHs", [3, 1, 1, 13]);

        // (D) c1ccccc1c1ccc(CCC)cc1 parsed with sanitize=false, then
        // sanitized with every implemented stage except KEKULIZE. MASK
        // FACT (R07-STAGE Action 2): RDKit SANITIZE_ALL (MolOps.h:534,
        // 0xFFFFFFF) INCLUDES reserved bits, exactly like CK ALL; the
        // named-bit mask below is a SUPPORTED CK CONSTRUCTION of the same
        // currently implemented selected stages (everything except
        // kekulize), NOT the identical raw RDKit mask value (CK from_bits
        // rejects arbitrary unknown-bit patterns such as ALL^KEKULIZE).
        let mut no_sanitize = cosmolkit_smiles::SmilesParseParams::default();
        no_sanitize.sanitize = false;
        no_sanitize.remove_hydrogens = false;
        let raw = cosmolkit_smiles::parse_smiles("c1ccccc1c1ccc(CCC)cc1", &no_sanitize)
            .expect("parse unsanitized biphenyl fixture")
            .topology;
        use cosmolkit_core::{SanitizeOperations, SanitizeParams, sanitize_topology};
        let named_all = SanitizeOperations::CLEANUP
            | SanitizeOperations::PROPERTIES
            | SanitizeOperations::SYMM_RINGS
            | SanitizeOperations::KEKULIZE
            | SanitizeOperations::FIND_RADICALS
            | SanitizeOperations::SET_AROMATICITY
            | SanitizeOperations::SET_CONJUGATION
            | SanitizeOperations::SET_HYBRIDIZATION
            | SanitizeOperations::CLEANUP_CHIRALITY
            | SanitizeOperations::ADJUST_HS
            | SanitizeOperations::CLEANUP_ORGANOMETALLICS
            | SanitizeOperations::CLEANUP_ATROPISOMERS;
        let all_minus_kekulize =
            SanitizeOperations::from_bits(named_all.bits() ^ SanitizeOperations::KEKULIZE.bits())
                .expect("valid reduced operation set");
        let partially_sanitized = sanitize_topology(
            &raw,
            &SanitizeParams {
                operations: all_minus_kekulize,
            },
        )
        .expect("partial sanitization (all implemented stages minus kekulize) succeeds")
        .topology;
        // Pre/post state assertions justified by source: kekulize is the
        // only implemented order-rewriting stage, so its EXCLUSION means
        // the fully-sanitized output (with kekulize) must DIFFER from the
        // partially-sanitized one — the two sanitize paths genuinely
        // produce distinct states on this input.
        let fully_sanitized_biphenyl = sanitize_topology(
            &raw,
            &SanitizeParams {
                operations: SanitizeOperations::ALL,
            },
        )
        .expect("full sanitization succeeds")
        .topology;
        assert_ne!(
            fully_sanitized_biphenyl, partially_sanitized,
            "kekulize exclusion must produce a distinct state"
        );
        check(
            &fully_sanitized_biphenyl,
            "biphenyl FULL sanitize",
            [3, 3, 3, 2],
        );
        check(
            &partially_sanitized,
            "biphenyl all-stages-minus-KEKULIZE",
            [3, 3, 3, 2],
        );

        // (A-pinned) The actual catch_tests "basics" mutation fixture:
        // FULL core sanitize output of the raw parse, then the bridge
        // bond-type-only mutation. PREREQUISITE (RDKit-justified, asserted
        // getter-free via set_* identity): the fully finalized/sanitized
        // inter-ring bridge(5,6) is Single AND non-aromatic (the pinned
        // pristine observable; per source, GetUnspecifiedBondType returns
        // AROMATIC when both endpoints are aromatic, so the SINGLE state
        // is a finalize/sanitize-stage outcome, not a parser-stage
        // attribution). If this prerequisite fails, the source-shaped
        // state is not producible through the existing owners and must be
        // REPORTED (no manual chemistry repair, no reversed mutation, no
        // weakened assertion).
        let bridge_is_single_nonaromatic = |topology: &TopologyBlock, a: usize, b: usize| {
            let (a, b) = (AtomId::new(a), AtomId::new(b));
            let bridge = topology
                .bonds
                .iter()
                .find(|bond| {
                    (bond.begin(), bond.end()) == (a, b) || (bond.begin(), bond.end()) == (b, a)
                })
                .expect("bridge bond exists")
                .clone();
            let mut single_probe = bridge.clone();
            single_probe.set_order(BondOrder::Single);
            let mut nonaromatic_probe = bridge.clone();
            nonaromatic_probe.set_aromatic(false);
            single_probe == bridge && nonaromatic_probe == bridge
        };
        assert!(
            bridge_is_single_nonaromatic(&fully_sanitized_biphenyl, 5, 6),
            "PREREQUISITE: fully-sanitized bridge(5,6) must be Single and \
             non-aromatic (RDKit pristine state)"
        );
        // Pristine: 3/3/3/2 (catch "basics" first three CHECKs).
        check(
            &fully_sanitized_biphenyl,
            "biphenyl FULL sanitize (pristine)",
            [3, 3, 3, 2],
        );
        // Mutated: clone the pristine, set ONLY the bridge order to
        // Aromatic (setBondType analogue), assert is_aromatic stays false
        // and the whole snapshot equals the pristine with exactly that one
        // bond-order change.
        let mut mutated_biphenyl = fully_sanitized_biphenyl.clone();
        let mut expected_biphenyl = fully_sanitized_biphenyl.clone();
        for (topology, is_expected) in [
            (&mut mutated_biphenyl, false),
            (&mut expected_biphenyl, true),
        ] {
            let _ = is_expected;
            for bond in &mut topology.bonds {
                if (bond.begin(), bond.end()) == (AtomId::new(5), AtomId::new(6))
                    || (bond.begin(), bond.end()) == (AtomId::new(6), AtomId::new(5))
                {
                    bond.set_order(BondOrder::Aromatic);
                }
            }
        }
        assert_eq!(
            mutated_biphenyl, expected_biphenyl,
            "snapshot: only the bridge bond order changed"
        );
        assert_ne!(
            mutated_biphenyl, fully_sanitized_biphenyl,
            "mutation is a real state change"
        );
        assert!(
            bridge_is_single_nonaromatic(
                // is_aromatic still false on the MUTATED bridge: only order
                // changed, verified by the nonaromatic identity probe below.
                &{
                    let mut probe_topology = mutated_biphenyl.clone();
                    for bond in &mut probe_topology.bonds {
                        if (bond.begin(), bond.end()) == (AtomId::new(5), AtomId::new(6))
                            || (bond.begin(), bond.end()) == (AtomId::new(6), AtomId::new(5))
                        {
                            bond.set_aromatic(false);
                        }
                    }
                    probe_topology
                },
                5,
                6
            ) || {
                // Direct getter-free check: setting is_aromatic(false) on the
                // mutated bridge is an identity <=> is_aromatic stayed false.
                let (a, b) = (AtomId::new(5), AtomId::new(6));
                let mutated_bridge = mutated_biphenyl
                    .bonds
                    .iter()
                    .find(|bond| {
                        (bond.begin(), bond.end()) == (a, b) || (bond.begin(), bond.end()) == (b, a)
                    })
                    .unwrap()
                    .clone();
                let mut probe = mutated_bridge.clone();
                probe.set_aromatic(false);
                probe == mutated_bridge
            }
        );
        check(
            &mutated_biphenyl,
            "biphenyl FULL sanitize + bridge AROMATIC (mutated)",
            [3, 3, 3, 2],
        );

        // (B-pinned) r5-linked molecule: full sanitize, prerequisite,
        // pristine + bridge(4,5)-order-only mutation, both 3/3/3/3
        // (StrictLinkages=3 source-asserted under the mutation).
        let r5_raw = cosmolkit_smiles::parse_smiles("c1cc[nH]c1c1[nH]c(CCC)cc1", &no_sanitize)
            .expect("parse unsanitized r5-link fixture")
            .topology;
        let fully_sanitized_r5 = sanitize_topology(
            &r5_raw,
            &SanitizeParams {
                operations: SanitizeOperations::ALL,
            },
        )
        .expect("full sanitize of r5-link succeeds")
        .topology;
        assert!(
            bridge_is_single_nonaromatic(&fully_sanitized_r5, 4, 5),
            "PREREQUISITE: fully-sanitized r5-link bridge(4,5) must be \
             Single and non-aromatic"
        );
        check(
            &fully_sanitized_r5,
            "r5-link FULL sanitize (pristine)",
            [3, 3, 3, 3],
        );
        let mut mutated_r5 = fully_sanitized_r5.clone();
        let mut expected_r5 = fully_sanitized_r5.clone();
        for topology in [&mut mutated_r5, &mut expected_r5] {
            for bond in &mut topology.bonds {
                if (bond.begin(), bond.end()) == (AtomId::new(4), AtomId::new(5))
                    || (bond.begin(), bond.end()) == (AtomId::new(5), AtomId::new(4))
                {
                    bond.set_order(BondOrder::Aromatic);
                }
            }
        }
        assert_eq!(
            mutated_r5, expected_r5,
            "snapshot: only the bridge order changed"
        );
        assert_ne!(
            mutated_r5, fully_sanitized_r5,
            "mutation is a real state change"
        );
        let (a, b) = (AtomId::new(4), AtomId::new(5));
        let mutated_r5_bridge = mutated_r5
            .bonds
            .iter()
            .find(|bond| {
                (bond.begin(), bond.end()) == (a, b) || (bond.begin(), bond.end()) == (b, a)
            })
            .unwrap()
            .clone();
        let mut r5_probe = mutated_r5_bridge.clone();
        r5_probe.set_aromatic(false);
        assert_eq!(
            r5_probe, mutated_r5_bridge,
            "is_aromatic stayed false after mutation"
        );
        check(
            &mutated_r5,
            "r5-link FULL sanitize + bridge AROMATIC (mutated)",
            [3, 3, 3, 3],
        );
    }

    #[test]
    fn descriptor_n01_num_rings_supplied_rows_and_ring_shapes() {
        // Source Lipinski.cpp:205-209: NumRingsVersion "1.0.1";
        // calcNumRings = the count of the molecule's SUPPLIED RingInfo
        // rows (RingInfo::numRings, reused from the core owner). The
        // descriptor reads the prepared input's ring rows — it NEVER
        // recalculates SSSR.
        assert_eq!(NUM_RINGS_VERSION, "1.0.1");

        // Ring-shape expectations derived from SSSR structure (source
        // semantics): empty = 0; single 6-ring = 1; fused (naphthalene,
        // 2 fused 6-rings: SSSR size = bonds - atoms + 1 = 11 - 10 + 1)
        // = 2; spiro (spiro[4.5]decane: two rings sharing exactly one
        // atom) = 2.
        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("C1CCCCC1", 1),
            ("c1ccc2ccccc2c1", 2),
            ("C1CCC2(CC1)CCCCC2", 2),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(num_rings(&topology).unwrap(), *expected, "cold {smiles}");
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_rings_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
            // Supplied-rows evidence: the prepared result is exactly the
            // SUPPLIED RingInfo row count (rows read, not recomputed — the
            // prepared kernel has no topology-derived recompute path).
            assert_eq!(
                rings::num_rings_prepared(&input).unwrap(),
                u32::try_from(input.ring_info().num_rings()).unwrap(),
                "supplied rows {smiles}"
            );
        }

        // N01-CLOSE distinguishing fixture: cubane C12C3C4C1C5C2C3C45 —
        // the SAME topology given two DIFFERENT legitimate ring row sets.
        // Source-backed derivation (frozen BEFORE any execution):
        // E - V + 1 = 12 - 8 + 1 = 5 independent cycles (SSSR), while the
        // symmetrized set retains all six square faces. The prepared
        // kernel must report the SUPPLIED rows (5 with find_sssr, 6 with
        // symmetrized_sssr) — this distinguishes row-reading from any
        // SSSR recomputation, which could only ever produce 5. The cold
        // wrapper computes its one canonical SSSR assignment -> 5.
        // migration_core_sssr.rs's cubane regression (5/6) is
        // supplementary reuse evidence, not the oracle here.
        let cubane = smiles_topology("C12C3C4C1C5C2C3C45");
        let params = cosmolkit_core::RingSearchParams {
            include_dative_bonds: false,
            include_hydrogen_bonds: false,
        };
        let sssr_rows = cosmolkit_core::find_sssr(&cubane, &params).unwrap();
        assert_eq!(sssr_rows.num_rings(), 5, "frozen derivation: E-V+1=5");
        let symmetric_rows = cosmolkit_core::symmetrized_sssr(&cubane, &params).unwrap();
        assert_eq!(symmetric_rows.num_rings(), 6, "frozen derivation: 6 faces");
        let assignment = valence(&cubane, "n01_cubane").unwrap();
        let (coordinates, properties) = (
            cosmolkit_model::CoordinateBlock::default(),
            cosmolkit_model::MoleculeProperties::default(),
        );
        let sssr_input =
            DescriptorInput::new(&cubane, &coordinates, &properties, &assignment, &sssr_rows);
        assert_eq!(rings::num_rings_prepared(&sssr_input).unwrap(), 5);
        let symmetric_input = DescriptorInput::new(
            &cubane,
            &coordinates,
            &properties,
            &assignment,
            &symmetric_rows,
        );
        assert_eq!(rings::num_rings_prepared(&symmetric_input).unwrap(), 6);
        // Cold wrapper: one canonical SSSR assignment -> 5, with ZERO
        // entries into the cold valence helper (test-only thread-local
        // counter; before/after difference around the cold invocation).
        let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(num_rings(&cubane).unwrap(), 5);
        let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(
            after - before,
            0,
            "cold num_rings must not enter the valence helper"
        );
        // Same zero-entry evidence around a second cold call on another
        // fixture (single ring) to pin the cost property, not one molecule.
        let cyclohexane = smiles_topology("C1CCCCC1");
        let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(num_rings(&cyclohexane).unwrap(), 1);
        let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(
            after - before,
            0,
            "cold num_rings must not enter the valence helper"
        );
    }

    #[test]
    fn descriptor_n02_num_heterocycles_any_member_non_carbon() {
        // Source Lipinski.cpp:232-242: NumHeterocyclesVersion "1.0.0";
        // calcNumHeterocycles counts a SUPPLIED atom-ring row ONCE when ANY
        // member atom has getAtomicNum() != 6 (first non-carbon member
        // breaks the inner loop). getAtomicNum() != 6 projects to the typed
        // element() != Element::C, so dummy atoms (atomicNum 0) count as
        // hetero members exactly as in the source.
        assert_eq!(NUM_HETEROCYCLES_VERSION, "1.0.0");

        const CASES: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("C1CCCCC1", 0),          // carbocycle: all members C
            ("c1ccc2ccccc2c1", 0),    // fused naphthalene: both rows all-C
            ("C1CCNCC1", 1),          // piperidine: N member -> 1
            ("C1COCC1", 1),           // tetrahydropyran: O member -> 1
            ("c1ccoc1", 1),           // furan: aromatic O member -> 1
            ("C1CCNCC1CC1CCCCC1", 1), // bicyclic piperidine + carbocycle:
            // two distinct rows judged
            // independently; only the N-ring has a
            // non-carbon member -> 1
            ("C1CCNCC1CC1CCOCC1", 2), // BOTH rings hetero (N ring + O
                                      // ring): two independent hetero rows -> 2
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_heterocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_heterocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Dummy-member ring (real dummy): `*` parses to Element::DUMMY
        // (atomic number 0); atomic-number0 satisfies getAtomicNum() != 6,
        // so the row counts. Raw parser rows with remove_hydrogens:false.
        let mut keep_h = cosmolkit_smiles::SmilesParseParams::default();
        keep_h.remove_hydrogens = false;
        let dummy_ring = cosmolkit_smiles::parse_smiles("*1CCCC1", &keep_h)
            .expect("parse real dummy-member ring fixture")
            .topology;
        assert_eq!(
            dummy_ring.atoms[0].element().atomic_number(),
            0,
            "real dummy atomic number 0 asserted"
        );
        let before_dummy = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(
            num_heterocycles(&dummy_ring).unwrap(),
            1,
            "dummy atomic-number-0 member counts as hetero"
        );
        let after_dummy = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(after_dummy - before_dummy, 0);
        let (a, r, c, p) = prepared_input_for(&dummy_ring);
        let dummy_input = DescriptorInput::new(&dummy_ring, &c, &p, &a, &r);
        assert_eq!(rings::num_heterocycles_prepared(&dummy_input).unwrap(), 1);

        // Metal regression ([Cu] = atomic number 29, NOT dummy — the
        // earlier "dummy (atomicNum 0)" description of this fixture was
        // WRONG and is withdrawn; Cu covers the ordinary non-carbon
        // metal branch only).
        let copper_ring = cosmolkit_smiles::parse_smiles("[Cu]1CCCC1", &keep_h)
            .expect("parse copper-member ring fixture")
            .topology;
        assert_eq!(
            copper_ring.atoms[0].element().atomic_number(),
            29,
            "Cu atomic number 29 asserted"
        );
        assert_eq!(
            num_heterocycles(&copper_ring).unwrap(),
            1,
            "non-carbon metal member counts as hetero"
        );

        // N02-CLOSE frozen position-by-identity product: base
        // five-member all-carbon cycle C1CCCC1; for each ring-row
        // position 0..4 crossed with atomic numbers {6,7,8,0}, a detached
        // clone with exactly that one atom's element set via the existing
        // Atom::set_element. Literal expectations {0,1,1,1} independent of
        // position. Selected atomic number asserted BEFORE evaluation;
        // exactly 40 descriptor calls (20 fixtures x cold+prepared).
        let base_cycle = smiles_topology("C1CCCC1");
        let params = cosmolkit_core::RingSearchParams {
            include_dative_bonds: false,
            include_hydrogen_bonds: false,
        };
        let base_rows = cosmolkit_core::find_sssr(&base_cycle, &params)
            .expect("base cycle SSSR")
            .atom_rings()
            .to_vec();
        assert_eq!(base_rows.len(), 1, "one ring row");
        let row: Vec<usize> = base_rows[0].iter().map(|atom_id| atom_id.index()).collect();
        assert_eq!(row.len(), 5, "five-member row");
        let mut descriptor_calls = 0u32;
        for &position in &row {
            // Frozen literal (Z, expected) pairs — no derived expectation.
            for (atomic_number, expected) in [(6u8, 0u32), (7, 1), (8, 1), (0, 1)] {
                let mut fixture = base_cycle.clone();
                fixture.atoms[position].set_element(
                    Element::from_atomic_number(atomic_number)
                        .expect("fixed valid vocabulary entry"),
                );
                assert_eq!(
                    fixture.atoms[position].element().atomic_number(),
                    atomic_number,
                    "selected atom identity asserted before evaluation"
                );
                // Per-product-cold counter checkpoint (before/after delta
                // must be 0 for EACH of the 20 cold calls).
                let before_cold = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    num_heterocycles(&fixture).unwrap(),
                    expected,
                    "cold position {position} Z={atomic_number}"
                );
                let after_cold = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    after_cold - before_cold,
                    0,
                    "product cold position {position} Z={atomic_number} must not enter the valence helper"
                );
                descriptor_calls += 1;
                let (assignment, ring_info, coordinates, properties) = prepared_input_for(&fixture);
                let input = DescriptorInput::new(
                    &fixture,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                assert_eq!(
                    rings::num_heterocycles_prepared(&input).unwrap(),
                    expected,
                    "prepared position {position} Z={atomic_number}"
                );
                descriptor_calls += 1;
            }
        }
        assert_eq!(descriptor_calls, 40, "exactly 40 descriptor calls");

        // N02-CLOSE frozen shared/member cases (cold + prepared, literal
        // expectations; the shared atom's membership in BOTH ring rows is
        // asserted for the first four, not merely the presence of two
        // rings):
        //   C12(CCC1)CCC2    -> 0 (shared C atom in both rows, all-carbon)
        //   [N+]12(CCC1)CCC2 -> 2 (shared N+ member: BOTH rows hetero)
        //   C12(CCN1)CCC2    -> 1 (N member only in one row)
        //   *12(CCC1)CCC2    -> 2 (shared dummy-0 member: BOTH rows hetero)
        //   N1NCCC1          -> 1 (two non-carbon members in ONE row count once)
        //   NC1CCCC1         -> 0 (hetero atom outside the ring rows)
        const SHARED: &[(&str, u32)] = &[
            ("C12(CCC1)CCC2", 0),
            ("[N+]12(CCC1)CCC2", 2),
            ("C12(CCN1)CCC2", 1),
            ("*12(CCC1)CCC2", 2),
            ("N1NCCC1", 1),
            ("NC1CCCC1", 0),
        ];
        for (smiles, expected) in SHARED {
            let topology = smiles_topology(smiles);
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            // Shared-row membership evidence for the first four cases: atom
            // index 0 belongs to BOTH supplied ring rows.
            if smiles.starts_with("C12") || smiles.starts_with("[N+]") || smiles.starts_with('*') {
                assert_eq!(ring_info.atom_rings().len(), 2, "two rows {smiles}");
                let shared = AtomId::new(0);
                assert!(
                    ring_info
                        .atom_rings()
                        .iter()
                        .all(|row| row.contains(&shared)),
                    "shared atom 0 belongs to both ring rows {smiles}"
                );
            }
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
            assert_eq!(
                num_heterocycles(&topology).unwrap(),
                *expected,
                "cold shared {smiles}"
            );
            let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
            assert_eq!(
                after - before,
                0,
                "cold num_heterocycles must not enter the valence helper"
            );
            assert_eq!(
                rings::num_heterocycles_prepared(&input).unwrap(),
                *expected,
                "prepared shared {smiles}"
            );
        }

        // N02-CLOSE distinguishing supplied-ring fixture: the all-N cube
        // N12N3N4N1N5N2N3N45 (same cube graph as the N01-CLOSE cubane, all
        // eight atom identities N — every atomic number asserted 7). For
        // this SAME topology, find_sssr gives 5 independent cycles and
        // symmetrized_sssr retains 6 square faces; every row is hetero, so
        // the prepared kernel must report the SUPPLIED rows (5 then 6) —
        // no recomputation can pass both assertions. Cold computes its one
        // canonical SSSR assignment -> 5. Zero valence-helper entries
        // around the cold call.
        let all_n_cube = smiles_topology("N12N3N4N1N5N2N3N45");
        assert!(
            all_n_cube
                .atoms
                .iter()
                .all(|atom| atom.element().atomic_number() == 7)
        );
        let cube_sssr = cosmolkit_core::find_sssr(&all_n_cube, &params).unwrap();
        assert_eq!(cube_sssr.num_rings(), 5, "frozen: E-V+1=5");
        let cube_symmetric = cosmolkit_core::symmetrized_sssr(&all_n_cube, &params).unwrap();
        assert_eq!(cube_symmetric.num_rings(), 6, "frozen: 6 faces");
        let cube_assignment = valence(&all_n_cube, "n02_cube").unwrap();
        let (cube_coordinates, cube_properties) = (
            cosmolkit_model::CoordinateBlock::default(),
            cosmolkit_model::MoleculeProperties::default(),
        );
        let cube_ss_input = DescriptorInput::new(
            &all_n_cube,
            &cube_coordinates,
            &cube_properties,
            &cube_assignment,
            &cube_sssr,
        );
        assert_eq!(rings::num_heterocycles_prepared(&cube_ss_input).unwrap(), 5);
        let cube_sym_input = DescriptorInput::new(
            &all_n_cube,
            &cube_coordinates,
            &cube_properties,
            &cube_assignment,
            &cube_symmetric,
        );
        assert_eq!(
            rings::num_heterocycles_prepared(&cube_sym_input).unwrap(),
            6
        );
        let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(num_heterocycles(&all_n_cube).unwrap(), 5);
        let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(
            after - before,
            0,
            "cold num_heterocycles must not enter the valence helper (cube)"
        );

        // Supplied-rows evidence: prepared result equals the independent
        // any-non-carbon row count over the SAME supplied rows and atoms.
        let piperidine = smiles_topology("C1CCNCC1");
        let (assignment, ring_info, coordinates, properties) = prepared_input_for(&piperidine);
        let input = DescriptorInput::new(
            &piperidine,
            &coordinates,
            &properties,
            &assignment,
            &ring_info,
        );
        let independent = ring_info
            .atom_rings()
            .iter()
            .filter(|row| {
                row.iter()
                    .any(|atom_id| piperidine.atoms[atom_id.index()].element() != Element::C)
            })
            .count();
        assert_eq!(
            rings::num_heterocycles_prepared(&input).unwrap(),
            u32::try_from(independent).unwrap()
        );
    }

    #[test]
    fn descriptor_n03_aromatic_rings_all_bond_flag_product() {
        // Source Lipinski.cpp:246-257: NumAromaticRingsVersion "1.0.0";
        // calcNumAromaticRings iterates the SUPPLIED bondRings() rows,
        // optimistically ++res per row, then --res at the FIRST member
        // bond with !getIsAromatic() and breaks: a row counts iff EVERY
        // member bond carries the aromatic FLAG — the boolean flag only,
        // never BondOrder, ring size or any chemical classification.
        // Frozen product: one benzene ring row (six member bonds), every
        // mask 0..63 over the six is_aromatic flags (member position p
        // <-> mask bit p) via the existing setter, cold AND prepared per
        // mask = 128 descriptor calls. Expectation table is literally 63
        // zeros followed by one 1: only mask 63 (all six flags set)
        // counts. Mixed masks cover every first-false retraction position
        // without case selection.
        assert_eq!(NUM_AROMATIC_RINGS_VERSION, "1.0.0");

        let pristine = smiles_topology("c1ccccc1");
        // ONE supplied six-bond ring row, read from the same cold owner
        // the wrappers delegate through.
        let base_rings = ring_info(&pristine, "descriptor_n03").unwrap();
        assert_eq!(base_rings.num_rings(), 1, "benzene: exactly one ring row");
        let row = base_rings.bond_rings()[0].clone();
        assert_eq!(row.len(), 6, "benzene row has six member bonds");

        // Literal source-derived expectation table: masks 0..62 -> 0,
        // mask 63 -> 1.
        #[rustfmt::skip]
        const EXPECTED: [u32; 64] = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        ];

        let mut descriptor_calls = 0usize;
        for mask in 0u32..64 {
            // Clone detached topology per mask; change ONLY the six member
            // bonds' is_aromatic flags (member position p <-> mask bit p).
            let mut masked = pristine.clone();
            for (bit, bond_id) in row.iter().enumerate() {
                masked.bonds[bond_id.index()].set_aromatic(((mask >> bit) & 1) == 1);
            }
            // Snapshot proof that ONLY those six flags changed: forcing
            // the SAME six flags to the same value on both blocks must
            // restore whole-block equality — BondOrder, endpoints,
            // adjacency, sgroups and every atom row included. This stays
            // independent of the pristine parse's own flag state.
            let mut restored_masked = masked.clone();
            let mut restored_pristine = pristine.clone();
            for bond_id in &row {
                restored_masked.bonds[bond_id.index()].set_aromatic(true);
                restored_pristine.bonds[bond_id.index()].set_aromatic(true);
            }
            assert_eq!(
                restored_masked, restored_pristine,
                "mask {mask}: only the six is_aromatic flags changed"
            );

            let expected = EXPECTED[mask as usize];
            // COLD: bracket each of the 64 cold calls with the test-only
            // valence-helper counter; delta must be literally 0.
            let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
            assert_eq!(
                num_aromatic_rings(&masked).unwrap(),
                expected,
                "cold mask {mask}"
            );
            let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
            assert_eq!(
                after - before,
                0,
                "cold mask {mask}: no valence-helper entry"
            );
            descriptor_calls += 1;
            // PREPARED: canonical matching input from the existing owners
            // (valence + SSSR + empty coordinate/property blocks).
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&masked);
            let input =
                DescriptorInput::new(&masked, &coordinates, &properties, &assignment, &ring_info);
            assert_eq!(
                rings::num_aromatic_rings_prepared(&input).unwrap(),
                expected,
                "prepared mask {mask}"
            );
            descriptor_calls += 1;
        }
        assert_eq!(descriptor_calls, 128, "64 masks x (cold + prepared)");

        // Source-shaped fixed cases (outside the 128-call product counter):
        // empty 0, acyclic 0, benzene 1 (all six parser-provided flags
        // aromatic), fused naphthalene 2 (both supplied rows fully
        // aromatic).
        const FIXED: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("c1ccccc1", 1),
            ("c1ccc2ccccc2c1", 2),
        ];
        for (smiles, expected) in FIXED {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_aromatic_rings(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_aromatic_rings_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_n04_saturated_rings_order_and_flag_arms() {
        // Source Lipinski.cpp:260-272: NumSaturatedRingsVersion "1.0.0";
        // calcNumSaturatedRings iterates the SUPPLIED bondRings() rows,
        // optimistically ++res per row, then --res at the FIRST member bond
        // with getBondType() != Bond::SINGLE OR getIsAromatic() and breaks.
        // The two arms are independent: a Single-ORDERED bond carrying the
        // aromatic flag retracts via the flag arm; a Double/Triple/Aromatic
        // ORDER retracts via the order arm. Frozen product: cyclohexane row
        // (six all-single members), each fault in {Double order, Triple
        // order, Single+aromatic-flag} at EACH of the six member positions
        // -> literally 0 for every fault; all-single control -> 1; dative
        // contrast: a Dative-ordered member removes the row from ring
        // finding entirely (include_dative_bonds:false) -> 0 by absence.
        assert_eq!(NUM_SATURATED_RINGS_VERSION, "1.0.0");

        let pristine = smiles_topology("C1CCCCC1");
        let base_rings = ring_info(&pristine, "descriptor_n04").unwrap();
        assert_eq!(base_rings.num_rings(), 1, "cyclohexane: exactly one row");
        let row = base_rings.bond_rings()[0].clone();
        assert_eq!(row.len(), 6, "cyclohexane row has six member bonds");

        // All-single control: 1 on BOTH paths (source all-single rule).
        let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(num_saturated_rings(&pristine).unwrap(), 1, "control cold");
        let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(after - before, 0, "control cold: no valence helper");
        let (assignment, control_ring_info, coordinates, properties) =
            prepared_input_for(&pristine);
        let input = DescriptorInput::new(
            &pristine,
            &coordinates,
            &properties,
            &assignment,
            &control_ring_info,
        );
        assert_eq!(
            rings::num_saturated_rings_prepared(&input).unwrap(),
            1,
            "control prepared"
        );

        // Complete small product: 6 member positions x 3 faults = 18
        // clones, each expected LITERALLY 0 (any failing member retracts
        // the row on its own arm). Position p <-> row position p.
        const FAULTS: &[&str] = &["double_order", "triple_order", "flag_on_single"];
        let mut descriptor_calls = 0usize;
        for (position, bond_id) in row.iter().enumerate() {
            for fault in FAULTS {
                let mut masked = pristine.clone();
                match *fault {
                    "double_order" => {
                        masked.bonds[bond_id.index()].set_order(BondOrder::Double);
                    }
                    "triple_order" => {
                        masked.bonds[bond_id.index()].set_order(BondOrder::Triple);
                    }
                    "flag_on_single" => {
                        masked.bonds[bond_id.index()].set_aromatic(true);
                    }
                    _ => unreachable!("frozen fault set"),
                }
                // Snapshot proof ONLY the one chosen field changed: revert
                // the single mutation and restore whole-block equality.
                let mut restored = masked.clone();
                match *fault {
                    "double_order" | "triple_order" => {
                        restored.bonds[bond_id.index()].set_order(BondOrder::Single);
                    }
                    "flag_on_single" => {
                        restored.bonds[bond_id.index()].set_aromatic(false);
                    }
                    _ => unreachable!("frozen fault set"),
                }
                assert_eq!(
                    restored, pristine,
                    "position {position} fault {fault}: only the chosen field changed"
                );

                let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    num_saturated_rings(&masked).unwrap(),
                    0,
                    "cold position {position} fault {fault}"
                );
                let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    after - before,
                    0,
                    "cold position {position} fault {fault}: no valence helper"
                );
                descriptor_calls += 1;
                let (assignment, ring_info, coordinates, properties) = prepared_input_for(&masked);
                let input = DescriptorInput::new(
                    &masked,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                assert_eq!(
                    rings::num_saturated_rings_prepared(&input).unwrap(),
                    0,
                    "prepared position {position} fault {fault}"
                );
                descriptor_calls += 1;
            }
        }
        assert_eq!(descriptor_calls, 36, "18 fault clones x (cold + prepared)");

        // Dative contrast: a Dative-ordered member is excluded from ring
        // finding by the canonical params (include_dative_bonds:false), so
        // the six-cycle degenerates to a path: NO bond-ring row exists and
        // the count is 0 by row ABSENCE (a different mechanism than the
        // order/flag retraction arms above).
        let mut dative = pristine.clone();
        dative.bonds[row[0].index()].set_order(BondOrder::Dative);
        let dative_rings = ring_info(&dative, "descriptor_n04_dative").unwrap();
        assert_eq!(
            dative_rings.num_rings(),
            0,
            "dative member removes the ring row entirely"
        );
        let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(num_saturated_rings(&dative).unwrap(), 0, "dative cold");
        let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(after - before, 0, "dative cold: no valence helper");
        let (assignment, dative_ring_info, coordinates, properties) = prepared_input_for(&dative);
        let input = DescriptorInput::new(
            &dative,
            &coordinates,
            &properties,
            &assignment,
            &dative_ring_info,
        );
        assert_eq!(
            rings::num_saturated_rings_prepared(&input).unwrap(),
            0,
            "dative prepared"
        );

        // Fixed source-shaped cases (outside the 36-call product counter):
        // empty/acyclic 0; all-single carbocycle and piperidine 1; one
        // Double retracts (cyclohexene); raw parsed benzene 0 (Aromatic
        // order and aromatic flags fail the source predicate; this helper
        // does not sanitize or kekulize); norbornane 2 (two SSSR all-single
        // rows); two disjoint all-single rings 2.
        const FIXED: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("C1CCCCC1", 1),
            ("C1CCNCC1", 1),
            ("C1=CCCCC1", 0),
            ("c1ccccc1", 0),
            ("C1CC2CCC1C2", 2),
            ("C1CCNCC1CC1CCCCC1", 2),
        ];
        for (smiles, expected) in FIXED {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_saturated_rings(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_saturated_rings_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_n05_aliphatic_rings_any_non_aromatic_member() {
        // Source Lipinski.cpp:275-286: NumAliphaticRingsVersion "1.0.0";
        // calcNumAliphaticRings iterates the SUPPLIED bondRings() rows with
        // NO optimistic increment — the ++res sits INSIDE the inner loop at
        // the FIRST member bond with !getIsAromatic(), then break: a row
        // counts ONCE iff it has AT LEAST ONE non-aromatic member bond.
        // This is the negation of the ALL-aromatic predicate (N03), NOT a
        // synonym of "not saturated": cyclohexene is unsaturated yet
        // counts here (its Double bond is non-aromatic); a fully aromatic
        // row counts 0. Frozen product: benzene row masks 0..63 over the
        // six is_aromatic flags — masks 0..62 -> 1, only mask 63 (all six
        // aromatic) -> 0, cold AND prepared per mask = 128 calls.
        assert_eq!(NUM_ALIPHATIC_RINGS_VERSION, "1.0.0");

        let pristine = smiles_topology("c1ccccc1");
        let base_rings = ring_info(&pristine, "descriptor_n05").unwrap();
        assert_eq!(base_rings.num_rings(), 1, "benzene: exactly one row");
        let row = base_rings.bond_rings()[0].clone();
        assert_eq!(row.len(), 6, "benzene row has six member bonds");

        // Literal source-derived expectation table: masks 0..62 -> 1,
        // mask 63 -> 0 (inverse of the N03 all-aromatic table).
        #[rustfmt::skip]
        const EXPECTED: [u32; 64] = [
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
        ];

        let mut descriptor_calls = 0usize;
        for mask in 0u32..64 {
            let mut masked = pristine.clone();
            for (bit, bond_id) in row.iter().enumerate() {
                masked.bonds[bond_id.index()].set_aromatic(((mask >> bit) & 1) == 1);
            }
            // Snapshot proof ONLY the six flags changed (force the same
            // six flags on both blocks, then whole-block equality).
            let mut restored_masked = masked.clone();
            let mut restored_pristine = pristine.clone();
            for bond_id in &row {
                restored_masked.bonds[bond_id.index()].set_aromatic(true);
                restored_pristine.bonds[bond_id.index()].set_aromatic(true);
            }
            assert_eq!(
                restored_masked, restored_pristine,
                "mask {mask}: only the six is_aromatic flags changed"
            );

            let expected = EXPECTED[mask as usize];
            let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
            assert_eq!(
                num_aliphatic_rings(&masked).unwrap(),
                expected,
                "cold mask {mask}"
            );
            let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
            assert_eq!(
                after - before,
                0,
                "cold mask {mask}: no valence-helper entry"
            );
            descriptor_calls += 1;
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&masked);
            let input =
                DescriptorInput::new(&masked, &coordinates, &properties, &assignment, &ring_info);
            assert_eq!(
                rings::num_aliphatic_rings_prepared(&input).unwrap(),
                expected,
                "prepared mask {mask}"
            );
            descriptor_calls += 1;
        }
        assert_eq!(descriptor_calls, 128, "64 masks x (cold + prepared)");

        // "Not saturated synonym" discriminators (source aromatic
        // negation, NOT the complement of the saturated count):
        // cyclohexene is UNSATURATED (num_saturated_rings == 0) yet still
        // counts as aliphatic here (its Double bond is non-aromatic);
        // benzene is aromatic (num_aromatic_rings == 1) and NOT aliphatic;
        // cyclohexane is BOTH saturated and aliphatic.
        let cyclohexane = smiles_topology("C1CCCCC1");
        let cyclohexene = smiles_topology("C1=CCCCC1");
        let benzene = smiles_topology("c1ccccc1");
        assert_eq!(num_saturated_rings(&cyclohexene).unwrap(), 0);
        assert_eq!(num_aliphatic_rings(&cyclohexene).unwrap(), 1);
        assert_eq!(num_aromatic_rings(&benzene).unwrap(), 1);
        assert_eq!(num_aliphatic_rings(&benzene).unwrap(), 0);
        assert_eq!(num_saturated_rings(&cyclohexane).unwrap(), 1);
        assert_eq!(num_aliphatic_rings(&cyclohexane).unwrap(), 1);

        // Fixed source-shaped cases (outside the 128-call product counter).
        const FIXED: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("C1CCCCC1", 1),
            ("C1=CCCCC1", 1),
            ("c1ccccc1", 0),
            ("c1ccc2ccccc2c1", 0),
            ("C1CCNCC1", 1),
            ("C12C3C4C1C5C2C3C45", 5),
        ];
        for (smiles, expected) in FIXED {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_aliphatic_rings(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_aliphatic_rings_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_n06_aromatic_heterocycles_all_four_combinations() {
        // Source Lipinski.cpp:288-311: NumAromaticHeterocyclesVersion
        // "1.0.0"; calcNumAromaticHeterocycles counts a SUPPLIED row iff
        // EVERY member bond is aromatic (first non-aromatic bond forces
        // countIt=false and breaks) AND at least one member bond has a
        // non-carbon endpoint (atomicNum != 6, sticky countIt=true). The
        // source's "checking each atom twice, kind of doofy" comment is
        // preserved verbatim in the kernel. All-4 combinations below; the
        // complete small product = two bases (all-carbon benzene row,
        // one-nitrogen pyridine row) x 64 six-flag masks each.
        assert_eq!(NUM_AROMATIC_HETEROCYCLES_VERSION, "1.0.0");

        // All-4 aromatic x hetero membership combinations (cold+prepared).
        const CASES: &[(&str, u32)] = &[
            ("c1ccncc1", 1), // aromatic + hetero (N) -> counts
            ("c1ccccc1", 0), // aromatic + all-carbon -> no hetero bond
            ("C1CCNCC1", 0), // non-aromatic + hetero -> first flag breaks
            ("C1CCCCC1", 0), // non-aromatic + all-carbon -> breaks first
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_aromatic_heterocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_aromatic_heterocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Complete small parameter product: for EACH base (benzene =
        // all-carbon row; pyridine = one N member among six atoms), apply
        // every 6-flag mask to the row's member bonds. Benzene: ALL 64
        // masks -> 0 (even all-aromatic lacks a hetero endpoint).
        // Pyridine: only mask 63 (all six aromatic) -> 1; every other
        // mask has a non-aromatic bond that breaks the scan.
        #[rustfmt::skip]
        const EXPECTED_BENZENE: [u32; 64] = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        #[rustfmt::skip]
        const EXPECTED_PYRIDINE: [u32; 64] = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        ];

        let benzene = smiles_topology("c1ccccc1");
        let pyridine = smiles_topology("c1ccncc1");
        let benzene_base = ring_info(&benzene, "descriptor_n06_benzene").unwrap();
        let pyridine_base = ring_info(&pyridine, "descriptor_n06_pyridine").unwrap();
        assert_eq!(benzene_base.num_rings(), 1);
        assert_eq!(pyridine_base.num_rings(), 1);
        let benzene_row = benzene_base.bond_rings()[0].clone();
        let pyridine_row = pyridine_base.bond_rings()[0].clone();
        assert_eq!(benzene_row.len(), 6);
        assert_eq!(pyridine_row.len(), 6);

        let mut descriptor_calls = 0usize;
        let bases: [(&TopologyBlock, &Vec<_>, &[u32; 64], &str); 2] = [
            (&benzene, &benzene_row, &EXPECTED_BENZENE, "benzene"),
            (&pyridine, &pyridine_row, &EXPECTED_PYRIDINE, "pyridine"),
        ];
        for (pristine, row, expected_table, base_name) in bases {
            for mask in 0u32..64 {
                let mut masked = pristine.clone();
                for (bit, bond_id) in row.iter().enumerate() {
                    masked.bonds[bond_id.index()].set_aromatic(((mask >> bit) & 1) == 1);
                }
                // Snapshot proof ONLY the six flags changed.
                let mut restored_masked = masked.clone();
                let mut restored_pristine = pristine.clone();
                for bond_id in row.iter() {
                    restored_masked.bonds[bond_id.index()].set_aromatic(true);
                    restored_pristine.bonds[bond_id.index()].set_aromatic(true);
                }
                assert_eq!(
                    restored_masked, restored_pristine,
                    "{base_name} mask {mask}: only the six flags changed"
                );

                let expected = expected_table[mask as usize];
                let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    num_aromatic_heterocycles(&masked).unwrap(),
                    expected,
                    "cold {base_name} mask {mask}"
                );
                let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    after - before,
                    0,
                    "cold {base_name} mask {mask}: no valence helper"
                );
                descriptor_calls += 1;
                let (assignment, ring_info, coordinates, properties) = prepared_input_for(&masked);
                let input = DescriptorInput::new(
                    &masked,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                assert_eq!(
                    rings::num_aromatic_heterocycles_prepared(&input).unwrap(),
                    expected,
                    "prepared {base_name} mask {mask}"
                );
                descriptor_calls += 1;
            }
        }
        assert_eq!(
            descriptor_calls, 256,
            "2 bases x 64 masks x (cold + prepared)"
        );

        // Fixed source-shaped extras (outside the 256-call product counter).
        const FIXED: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("c1ccoc1", 1),           // furan: aromatic O member
            ("c1ccsc1", 1),           // thiophene: aromatic S member
            ("c1nccnc1", 1),          // pyrimidine: two N members, still ONE row
            ("c1ccc2ccccc2c1", 0),    // naphthalene: aromatic but all-carbon
            ("C1CCNCC1CC1CCCCC1", 0), // both rows non-aromatic
        ];
        for (smiles, expected) in FIXED {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_aromatic_heterocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_aromatic_heterocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_n07_aromatic_carbocycles_all_four_combinations() {
        // Source Lipinski.cpp:312-333: NumAromaticCarbocyclesVersion
        // "1.0.0"; calcNumAromaticCarbocycles counts a SUPPLIED row iff
        // EVERY member bond is aromatic AND every member bond endpoint is
        // carbon (countIt starts TRUE; either failure arm breaks with
        // false). The inline comment here reads "big time sync" (vs N06's
        // "time sink") — preserved verbatim in the kernel. All-4 aromatic
        // x carbon-only combinations below; complete small product = two
        // bases (benzene all-C row, pyridine one-N row) x 64 six-flag
        // masks each.
        assert_eq!(NUM_AROMATIC_CARBOCYCLES_VERSION, "1.0.0");

        // All-4 aromatic x carbon-only combinations (cold+prepared).
        const CASES: &[(&str, u32)] = &[
            ("c1ccccc1", 1), // aromatic + all-carbon -> counts
            ("c1ccncc1", 0), // aromatic + hetero (N) -> endpoint arm breaks
            ("C1CCCCC1", 0), // non-aromatic + all-carbon -> flag arm breaks
            ("C1CCNCC1", 0), // non-aromatic + hetero -> flag arm first
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_aromatic_carbocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_aromatic_carbocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Complete small parameter product: benzene row — only mask 63
        // (all six aromatic, all-C) -> 1; pyridine row — ALL 64 masks -> 0
        // (the N endpoint breaks the second arm even when all aromatic).
        #[rustfmt::skip]
        const EXPECTED_BENZENE: [u32; 64] = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        ];
        #[rustfmt::skip]
        const EXPECTED_PYRIDINE: [u32; 64] = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];

        let benzene = smiles_topology("c1ccccc1");
        let pyridine = smiles_topology("c1ccncc1");
        let benzene_base = ring_info(&benzene, "descriptor_n07_benzene").unwrap();
        let pyridine_base = ring_info(&pyridine, "descriptor_n07_pyridine").unwrap();
        assert_eq!(benzene_base.num_rings(), 1);
        assert_eq!(pyridine_base.num_rings(), 1);
        let benzene_row = benzene_base.bond_rings()[0].clone();
        let pyridine_row = pyridine_base.bond_rings()[0].clone();
        assert_eq!(benzene_row.len(), 6);
        assert_eq!(pyridine_row.len(), 6);

        let mut descriptor_calls = 0usize;
        let bases: [(&TopologyBlock, &Vec<_>, &[u32; 64], &str); 2] = [
            (&benzene, &benzene_row, &EXPECTED_BENZENE, "benzene"),
            (&pyridine, &pyridine_row, &EXPECTED_PYRIDINE, "pyridine"),
        ];
        for (pristine, row, expected_table, base_name) in bases {
            for mask in 0u32..64 {
                let mut masked = pristine.clone();
                for (bit, bond_id) in row.iter().enumerate() {
                    masked.bonds[bond_id.index()].set_aromatic(((mask >> bit) & 1) == 1);
                }
                let mut restored_masked = masked.clone();
                let mut restored_pristine = pristine.clone();
                for bond_id in row.iter() {
                    restored_masked.bonds[bond_id.index()].set_aromatic(true);
                    restored_pristine.bonds[bond_id.index()].set_aromatic(true);
                }
                assert_eq!(
                    restored_masked, restored_pristine,
                    "{base_name} mask {mask}: only the six flags changed"
                );

                let expected = expected_table[mask as usize];
                let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    num_aromatic_carbocycles(&masked).unwrap(),
                    expected,
                    "cold {base_name} mask {mask}"
                );
                let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    after - before,
                    0,
                    "cold {base_name} mask {mask}: no valence helper"
                );
                descriptor_calls += 1;
                let (assignment, ring_info, coordinates, properties) = prepared_input_for(&masked);
                let input = DescriptorInput::new(
                    &masked,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                assert_eq!(
                    rings::num_aromatic_carbocycles_prepared(&input).unwrap(),
                    expected,
                    "prepared {base_name} mask {mask}"
                );
                descriptor_calls += 1;
            }
        }
        assert_eq!(
            descriptor_calls, 256,
            "2 bases x 64 masks x (cold + prepared)"
        );

        // Fixed source-shaped extras (outside the 256-call product counter).
        const FIXED: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("c1ccc2ccccc2c1", 2),    // naphthalene: both rows all-C aromatic
            ("c1ccoc1", 0),           // furan: hetero endpoint
            ("c1ccsc1", 0),           // thiophene: hetero endpoint
            ("C1CCNCC1CC1CCCCC1", 0), // non-aromatic rows
        ];
        for (smiles, expected) in FIXED {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_aromatic_carbocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_aromatic_carbocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_n08_aliphatic_heterocycles_all_four_combinations() {
        // Source Lipinski.cpp:334-358: NumAliphaticHeterocyclesVersion
        // "1.0.0"; calcNumAliphaticHeterocycles counts a SUPPLIED row iff
        // it has BOTH at least one NON-aromatic member bond (hasAliph, any
        // position) AND at least one hetero-endpoint member bond (hasHetero,
        // sticky, any position). The inner scan visits ALL members — no
        // break (the two flags can be set by DIFFERENT bonds). Complete
        // small product = two bases (benzene all-C row, pyridine one-N row)
        // x 64 six-flag masks each.
        assert_eq!(NUM_ALIPHATIC_HETEROCYCLES_VERSION, "1.0.0");

        // All-4 aliphatic x hetero combinations (cold+prepared).
        const CASES: &[(&str, u32)] = &[
            ("C1CCNCC1", 1), // non-aromatic + hetero (N) -> both flags
            ("c1ccncc1", 0), // fully aromatic + hetero -> no aliph
            ("C1CCCCC1", 0), // non-aromatic + all-carbon -> no hetero
            ("c1ccccc1", 0), // fully aromatic + all-carbon -> neither
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_aliphatic_heterocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_aliphatic_heterocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Complete small parameter product: pyridine (hetero present) —
        // every mask EXCEPT 63 (all aromatic => no aliph) counts 1;
        // benzene (no hetero) — ALL 64 masks -> 0.
        #[rustfmt::skip]
        const EXPECTED_PYRIDINE: [u32; 64] = [
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
        ];
        #[rustfmt::skip]
        const EXPECTED_BENZENE: [u32; 64] = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];

        let benzene = smiles_topology("c1ccccc1");
        let pyridine = smiles_topology("c1ccncc1");
        let benzene_base = ring_info(&benzene, "descriptor_n08_benzene").unwrap();
        let pyridine_base = ring_info(&pyridine, "descriptor_n08_pyridine").unwrap();
        assert_eq!(benzene_base.num_rings(), 1);
        assert_eq!(pyridine_base.num_rings(), 1);
        let benzene_row = benzene_base.bond_rings()[0].clone();
        let pyridine_row = pyridine_base.bond_rings()[0].clone();
        assert_eq!(benzene_row.len(), 6);
        assert_eq!(pyridine_row.len(), 6);

        let mut descriptor_calls = 0usize;
        let bases: [(&TopologyBlock, &Vec<_>, &[u32; 64], &str); 2] = [
            (&pyridine, &pyridine_row, &EXPECTED_PYRIDINE, "pyridine"),
            (&benzene, &benzene_row, &EXPECTED_BENZENE, "benzene"),
        ];
        for (pristine, row, expected_table, base_name) in bases {
            for mask in 0u32..64 {
                let mut masked = pristine.clone();
                for (bit, bond_id) in row.iter().enumerate() {
                    masked.bonds[bond_id.index()].set_aromatic(((mask >> bit) & 1) == 1);
                }
                let mut restored_masked = masked.clone();
                let mut restored_pristine = pristine.clone();
                for bond_id in row.iter() {
                    restored_masked.bonds[bond_id.index()].set_aromatic(true);
                    restored_pristine.bonds[bond_id.index()].set_aromatic(true);
                }
                assert_eq!(
                    restored_masked, restored_pristine,
                    "{base_name} mask {mask}: only the six flags changed"
                );

                let expected = expected_table[mask as usize];
                let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    num_aliphatic_heterocycles(&masked).unwrap(),
                    expected,
                    "cold {base_name} mask {mask}"
                );
                let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    after - before,
                    0,
                    "cold {base_name} mask {mask}: no valence helper"
                );
                descriptor_calls += 1;
                let (assignment, ring_info, coordinates, properties) = prepared_input_for(&masked);
                let input = DescriptorInput::new(
                    &masked,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                assert_eq!(
                    rings::num_aliphatic_heterocycles_prepared(&input).unwrap(),
                    expected,
                    "prepared {base_name} mask {mask}"
                );
                descriptor_calls += 1;
            }
        }
        assert_eq!(
            descriptor_calls, 256,
            "2 bases x 64 masks x (cold + prepared)"
        );

        // Fixed source-shaped extras (outside the 256-call product counter).
        const FIXED: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("c1ccoc1", 0),  // furan: fully aromatic -> no aliph
            ("c1ccsc1", 0),  // thiophene: fully aromatic -> no aliph
            ("C1COCCN1", 1), // morpholine: all non-aromatic, O+N
            ("C1CNCCN1", 1), // piperazine: all non-aromatic, two N
        ];
        for (smiles, expected) in FIXED {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_aliphatic_heterocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_aliphatic_heterocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_n09_aliphatic_carbocycles_all_four_combinations() {
        // Source Lipinski.cpp:360-382: NumAliphaticCarbocyclesVersion
        // "1.0.0"; calcNumAliphaticCarbocycles counts a SUPPLIED row iff
        // it has at least one NON-aromatic member bond (hasAliph, set
        // without break) AND NO hetero-endpoint member bond (a hetero
        // endpoint sets hasHetero and BREAKS immediately — the row is
        // already disqualified, so the early exit cannot change the
        // outcome). Fixtures here are parse-stage rows (aromatic flags
        // and Aromatic bond orders come from parsing); no sanitize or
        // kekulize stage is claimed. Complete small product = two bases
        // (benzene all-C row, pyridine one-N row) x 64 six-flag masks.
        assert_eq!(NUM_ALIPHATIC_CARBOCYCLES_VERSION, "1.0.0");

        // All-4 aliphatic x carbon-only combinations (cold+prepared).
        const CASES: &[(&str, u32)] = &[
            ("C1CCCCC1", 1), // non-aromatic + all-carbon -> both conditions
            ("c1ccncc1", 0), // aromatic + hetero (N) -> hetero breaks
            ("c1ccccc1", 0), // aromatic + all-carbon -> no aliph member
            ("C1CCNCC1", 0), // non-aromatic + hetero -> hetero breaks
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_aliphatic_carbocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_aliphatic_carbocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Complete small parameter product: benzene (all-C) — every mask
        // EXCEPT 63 (all aromatic => no aliph) counts 1; pyridine (one N)
        // — ALL 64 masks -> 0 (the N-endpoint bond breaks immediately).
        #[rustfmt::skip]
        const EXPECTED_BENZENE: [u32; 64] = [
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
        ];
        #[rustfmt::skip]
        const EXPECTED_PYRIDINE: [u32; 64] = [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];

        let benzene = smiles_topology("c1ccccc1");
        let pyridine = smiles_topology("c1ccncc1");
        let benzene_base = ring_info(&benzene, "descriptor_n09_benzene").unwrap();
        let pyridine_base = ring_info(&pyridine, "descriptor_n09_pyridine").unwrap();
        assert_eq!(benzene_base.num_rings(), 1);
        assert_eq!(pyridine_base.num_rings(), 1);
        let benzene_row = benzene_base.bond_rings()[0].clone();
        let pyridine_row = pyridine_base.bond_rings()[0].clone();
        assert_eq!(benzene_row.len(), 6);
        assert_eq!(pyridine_row.len(), 6);

        let mut descriptor_calls = 0usize;
        let bases: [(&TopologyBlock, &Vec<_>, &[u32; 64], &str); 2] = [
            (&benzene, &benzene_row, &EXPECTED_BENZENE, "benzene"),
            (&pyridine, &pyridine_row, &EXPECTED_PYRIDINE, "pyridine"),
        ];
        for (pristine, row, expected_table, base_name) in bases {
            for mask in 0u32..64 {
                let mut masked = pristine.clone();
                for (bit, bond_id) in row.iter().enumerate() {
                    masked.bonds[bond_id.index()].set_aromatic(((mask >> bit) & 1) == 1);
                }
                let mut restored_masked = masked.clone();
                let mut restored_pristine = pristine.clone();
                for bond_id in row.iter() {
                    restored_masked.bonds[bond_id.index()].set_aromatic(true);
                    restored_pristine.bonds[bond_id.index()].set_aromatic(true);
                }
                assert_eq!(
                    restored_masked, restored_pristine,
                    "{base_name} mask {mask}: only the six flags changed"
                );

                let expected = expected_table[mask as usize];
                let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    num_aliphatic_carbocycles(&masked).unwrap(),
                    expected,
                    "cold {base_name} mask {mask}"
                );
                let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    after - before,
                    0,
                    "cold {base_name} mask {mask}: no valence helper"
                );
                descriptor_calls += 1;
                let (assignment, ring_info, coordinates, properties) = prepared_input_for(&masked);
                let input = DescriptorInput::new(
                    &masked,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                assert_eq!(
                    rings::num_aliphatic_carbocycles_prepared(&input).unwrap(),
                    expected,
                    "prepared {base_name} mask {mask}"
                );
                descriptor_calls += 1;
            }
        }
        assert_eq!(
            descriptor_calls, 256,
            "2 bases x 64 masks x (cold + prepared)"
        );

        // Fixed source-shaped extras (outside the 256-call product counter).
        const FIXED: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("C1=CCCCC1", 1),          // cyclohexene: non-aromatic all-C
            ("C12C3C4C1C5C2C3C45", 5), // cubane: five all-single all-C rows
            ("c1ccc2ccccc2c1", 0),     // naphthalene: fully aromatic
            ("c1ccoc1", 0),            // furan: hetero endpoint
            ("c1ccsc1", 0),            // thiophene: hetero endpoint
        ];
        for (smiles, expected) in FIXED {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_aliphatic_carbocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_aliphatic_carbocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_n10_saturated_heterocycles_all_four_combinations() {
        // Source Lipinski.cpp:384-407: NumSaturatedHeterocyclesVersion
        // "1.0.0"; calcNumSaturatedHeterocycles counts a SUPPLIED row iff
        // EVERY member bond is Single-ordered AND non-aromatic (the N04
        // saturated conjunct — both arms independent, first failure breaks
        // with countIt=false) AND at least one member bond has a
        // non-carbon endpoint (sticky countIt=true arm). Fixtures are
        // parse-stage rows; no sanitize/kekulize stage is claimed.
        assert_eq!(NUM_SATURATED_HETEROCYCLES_VERSION, "1.0.0");

        // All-4 saturated x hetero combinations (cold+prepared).
        const CASES: &[(&str, u32)] = &[
            ("C1CCNCC1", 1), // saturated + hetero (N) -> counts
            ("C1CCCCC1", 0), // saturated + all-carbon -> no hetero bond
            ("c1ccncc1", 0), // aromatic + hetero -> order/flag arm breaks
            ("c1ccccc1", 0), // aromatic + all-carbon -> order arm breaks
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_saturated_heterocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_saturated_heterocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Complete small product: piperidine row (hetero present), each
        // saturated-fault in {Double order, Triple order, flag-on-Single}
        // at EACH of the six member positions -> literally 0 (the
        // saturated conjunct breaks first regardless of the hetero arm);
        // all-Single control -> 1.
        let pristine = smiles_topology("C1CCNCC1");
        let base_rings = ring_info(&pristine, "descriptor_n10").unwrap();
        assert_eq!(base_rings.num_rings(), 1, "piperidine: exactly one row");
        let row = base_rings.bond_rings()[0].clone();
        assert_eq!(row.len(), 6, "piperidine row has six member bonds");

        let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(
            num_saturated_heterocycles(&pristine).unwrap(),
            1,
            "control cold"
        );
        let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(after - before, 0, "control cold: no valence helper");
        let (assignment, control_ring_info, coordinates, properties) =
            prepared_input_for(&pristine);
        let input = DescriptorInput::new(
            &pristine,
            &coordinates,
            &properties,
            &assignment,
            &control_ring_info,
        );
        assert_eq!(
            rings::num_saturated_heterocycles_prepared(&input).unwrap(),
            1,
            "control prepared"
        );

        const FAULTS: &[&str] = &["double_order", "triple_order", "flag_on_single"];
        let mut descriptor_calls = 0usize;
        for (position, bond_id) in row.iter().enumerate() {
            for fault in FAULTS {
                let mut masked = pristine.clone();
                match *fault {
                    "double_order" => {
                        masked.bonds[bond_id.index()].set_order(BondOrder::Double);
                    }
                    "triple_order" => {
                        masked.bonds[bond_id.index()].set_order(BondOrder::Triple);
                    }
                    "flag_on_single" => {
                        masked.bonds[bond_id.index()].set_aromatic(true);
                    }
                    _ => unreachable!("frozen fault set"),
                }
                let mut restored = masked.clone();
                match *fault {
                    "double_order" | "triple_order" => {
                        restored.bonds[bond_id.index()].set_order(BondOrder::Single);
                    }
                    "flag_on_single" => {
                        restored.bonds[bond_id.index()].set_aromatic(false);
                    }
                    _ => unreachable!("frozen fault set"),
                }
                assert_eq!(
                    restored, pristine,
                    "position {position} fault {fault}: only the chosen field changed"
                );

                let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    num_saturated_heterocycles(&masked).unwrap(),
                    0,
                    "cold position {position} fault {fault}"
                );
                let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    after - before,
                    0,
                    "cold position {position} fault {fault}: no valence helper"
                );
                descriptor_calls += 1;
                let (assignment, ring_info, coordinates, properties) = prepared_input_for(&masked);
                let input = DescriptorInput::new(
                    &masked,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                assert_eq!(
                    rings::num_saturated_heterocycles_prepared(&input).unwrap(),
                    0,
                    "prepared position {position} fault {fault}"
                );
                descriptor_calls += 1;
            }
        }
        assert_eq!(descriptor_calls, 36, "18 fault clones x (cold + prepared)");

        // Fixed source-shaped extras (outside the 36-call product counter).
        const FIXED: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("C1COCCN1", 1),          // morpholine: saturated + O/N
            ("C1CNCCN1", 1),          // piperazine: saturated + two N
            ("C1=CCCCC1", 0),         // cyclohexene: Double retracts
            ("c1ccoc1", 0),           // furan: aromatic retracts
            ("c1ccsc1", 0),           // thiophene: aromatic retracts
            ("C1CCNCC1CC1CCCCC1", 1), // N-ring counts; carbocycle does not
        ];
        for (smiles, expected) in FIXED {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_saturated_heterocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_saturated_heterocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_n11_saturated_carbocycles_all_four_combinations() {
        // Source Lipinski.cpp:409-432: NumSaturatedCarbocyclesVersion
        // "1.0.0"; calcNumSaturatedCarbocycles counts a SUPPLIED row iff
        // EVERY member bond is Single-ordered, non-aromatic AND has two
        // carbon endpoints (countIt starts true; break-with-false at the
        // FIRST failure on either arm). DISTINCT from the N04 saturated
        // count, which ignores endpoints: piperidine counts 1 there but 0
        // here. Fixtures are parse-stage rows.
        assert_eq!(NUM_SATURATED_CARBOCYCLES_VERSION, "1.0.0");

        // All-4 saturated x carbon-only combinations (cold+prepared).
        const CASES: &[(&str, u32)] = &[
            ("C1CCCCC1", 1), // saturated + all-carbon -> counts
            ("C1CCNCC1", 0), // saturated + hetero (N) -> endpoint arm
            ("c1ccccc1", 0), // aromatic + all-carbon -> saturated arm
            ("c1ccncc1", 0), // aromatic + hetero -> saturated arm first
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_saturated_carbocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_saturated_carbocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // Complete small product: cyclohexane row (all-carbon), each
        // saturated-fault in {Double order, Triple order, flag-on-Single}
        // at EACH of the six member positions -> literally 0; all-Single
        // control -> 1.
        let pristine = smiles_topology("C1CCCCC1");
        let base_rings = ring_info(&pristine, "descriptor_n11").unwrap();
        assert_eq!(base_rings.num_rings(), 1, "cyclohexane: exactly one row");
        let row = base_rings.bond_rings()[0].clone();
        assert_eq!(row.len(), 6, "cyclohexane row has six member bonds");

        let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(
            num_saturated_carbocycles(&pristine).unwrap(),
            1,
            "control cold"
        );
        let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(after - before, 0, "control cold: no valence helper");
        let (assignment, control_ring_info, coordinates, properties) =
            prepared_input_for(&pristine);
        let input = DescriptorInput::new(
            &pristine,
            &coordinates,
            &properties,
            &assignment,
            &control_ring_info,
        );
        assert_eq!(
            rings::num_saturated_carbocycles_prepared(&input).unwrap(),
            1,
            "control prepared"
        );

        const FAULTS: &[&str] = &["double_order", "triple_order", "flag_on_single"];
        let mut descriptor_calls = 0usize;
        for (position, bond_id) in row.iter().enumerate() {
            for fault in FAULTS {
                let mut masked = pristine.clone();
                match *fault {
                    "double_order" => {
                        masked.bonds[bond_id.index()].set_order(BondOrder::Double);
                    }
                    "triple_order" => {
                        masked.bonds[bond_id.index()].set_order(BondOrder::Triple);
                    }
                    "flag_on_single" => {
                        masked.bonds[bond_id.index()].set_aromatic(true);
                    }
                    _ => unreachable!("frozen fault set"),
                }
                let mut restored = masked.clone();
                match *fault {
                    "double_order" | "triple_order" => {
                        restored.bonds[bond_id.index()].set_order(BondOrder::Single);
                    }
                    "flag_on_single" => {
                        restored.bonds[bond_id.index()].set_aromatic(false);
                    }
                    _ => unreachable!("frozen fault set"),
                }
                assert_eq!(
                    restored, pristine,
                    "position {position} fault {fault}: only the chosen field changed"
                );

                let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    num_saturated_carbocycles(&masked).unwrap(),
                    0,
                    "cold position {position} fault {fault}"
                );
                let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
                assert_eq!(
                    after - before,
                    0,
                    "cold position {position} fault {fault}: no valence helper"
                );
                descriptor_calls += 1;
                let (assignment, ring_info, coordinates, properties) = prepared_input_for(&masked);
                let input = DescriptorInput::new(
                    &masked,
                    &coordinates,
                    &properties,
                    &assignment,
                    &ring_info,
                );
                assert_eq!(
                    rings::num_saturated_carbocycles_prepared(&input).unwrap(),
                    0,
                    "prepared position {position} fault {fault}"
                );
                descriptor_calls += 1;
            }
        }
        assert_eq!(descriptor_calls, 36, "18 fault clones x (cold + prepared)");

        // Fixed source-shaped extras (outside the 36-call product counter).
        const FIXED: &[(&str, u32)] = &[
            ("", 0),
            ("CCCCCC", 0),
            ("C1COCCN1", 0),           // morpholine: hetero endpoints
            ("C1CNCCN1", 0),           // piperazine: hetero endpoints
            ("C1=CCCCC1", 0),          // cyclohexene: Double retracts
            ("C12C3C4C1C5C2C3C45", 5), // cubane: five all-single all-C rows
            ("c1ccoc1", 0),            // furan: saturated arm
            ("C1CCNCC1CC1CCCCC1", 1),  // carbocycle counts; N-ring does not
        ];
        for (smiles, expected) in FIXED {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_saturated_carbocycles(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_saturated_carbocycles_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }
    }

    #[test]
    fn descriptor_n12_spiro_atoms_acquisition_ids_and_out_state() {
        // Source Lipinski.cpp:435-464: NumSpiroAtomsVersion "1.0.0";
        // calcNumSpiroAtoms ensures SSSR-or-better ring info (absent or
        // weaker-than-SSSR => recompute), then scans ALL unordered atom-ring
        // row pairs: a pair sharing EXACTLY ONE atom contributes that atom,
        // deduped (std::find) into the optional out vector whose length is
        // returned — a prepopulated vector keeps its entries, dedups against
        // them and inflates the count. Fused pairs (|i^j|=2) and larger
        // shares never contribute.
        assert_eq!(NUM_SPIRO_ATOMS_VERSION, "1.0.0");

        // Fixed contrasts, cold + prepared (the prepared path here supplies
        // SSSR-or-better rows, read AS SUPPLIED — the family convention).
        const CASES: &[(&str, u32)] = &[
            ("C1CCC2(C1)CCCCC2", 1),   // spiro[4.5]decane: one shared atom
            ("c1ccc2ccccc2c1", 0),     // fused naphthalene: |i^j| = 2
            ("C1CC2CCC1C2", 0),        // norbornane: bridgehead pair shared
            ("C12C3C4C1C5C2C3C45", 0), // cubane: face pairs share 0 or 2
            ("C1CCNCC1CC1CCCCC1", 0),  // disjoint rings
            ("C1CCCCC1", 0),           // single ring: no pairs
            ("", 0),
            ("CCCCCC", 0),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_spiro_atoms(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_spiro_atoms_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // SSSR ACQUISITION, weaker-than-SSSR arm: a prepared input whose
        // supplied RingInfo is OtherOrUnknown (empty rows, NOT
        // sssr-or-better) must trigger the ONE canonical recompute and
        // still return the correct count. The absent arm is the cold path
        // above. The sssr-or-better arm is the normal prepared path above.
        let spiro = smiles_topology("C1CCC2(C1)CCCCC2");
        let weak_rings = cosmolkit_core::RingInfo::new(
            cosmolkit_core::RingFindType::OtherOrUnknown,
            spiro.atoms.len(),
            spiro.bonds.len(),
        );
        assert!(
            !weak_rings.is_sssr_or_better(),
            "fixture trigger: weaker-than-SSSR ring info"
        );
        let assignment = valence(&spiro, "descriptor_n12_weak").unwrap();
        let weak_coordinates = cosmolkit_model::CoordinateBlock::default();
        let weak_properties = cosmolkit_model::MoleculeProperties::default();
        let weak_input = DescriptorInput::new(
            &spiro,
            &weak_coordinates,
            &weak_properties,
            &assignment,
            &weak_rings,
        );
        assert_eq!(
            rings::num_spiro_atoms_prepared(&weak_input).unwrap(),
            1,
            "weaker-than-SSSR acquisition recomputes and counts the spiro atom"
        );

        // EXACT ordered unique atom IDs: the single returned ID is the one
        // atom present in BOTH atom-ring rows (identity derived from the
        // row data, not from the function under test).
        let base_rings = ring_info(&spiro, "descriptor_n12_ids").unwrap();
        let rows = base_rings.atom_rings();
        assert_eq!(rows.len(), 2, "spiro[4.5]decane: exactly two rows");
        let shared: Vec<_> = rows[0]
            .iter()
            .filter(|atom_id| rows[1].contains(atom_id))
            .collect();
        assert_eq!(shared.len(), 1, "the rows share exactly one atom");
        let expected_spiro = *shared[0];
        let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        let ids = spiro_atom_ids(&spiro).unwrap();
        let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(after - before, 0, "cold spiro_atom_ids: no valence helper");
        assert_eq!(ids.len(), 1);
        assert_eq!(ids[0], expected_spiro, "exact spiro atom identity");
        assert_eq!(num_spiro_atoms(&spiro).unwrap(), 1);
        assert!(spiro_atom_ids(&smiles_topology("")).unwrap().is_empty());

        // PREPOPULATED output state (the source's out-parameter semantics,
        // exercised at the kernel boundary): a preexisting entry is kept,
        // deduped against, and INCLUDED in the returned length.
        let mut dedup_against = vec![expected_spiro];
        assert_eq!(
            rings::spiro_atom_ids_kernel(&base_rings, &mut dedup_against).unwrap(),
            1,
            "prepopulated with the true spiro atom: no duplicate push"
        );
        assert_eq!(dedup_against.len(), 1);
        let mut foreign = vec![AtomId::new(0)];
        assert_eq!(
            rings::spiro_atom_ids_kernel(&base_rings, &mut foreign).unwrap(),
            2,
            "prepopulated with a non-spiro atom: kept + appended, count 2"
        );
        assert_eq!(foreign.len(), 2);
        assert_eq!(foreign[0], AtomId::new(0), "preexisting entry preserved");
        assert_eq!(foreign[1], expected_spiro, "append order after preexisting");

        // Complete small parameter product: every pair-intersection
        // cardinality class over the disjoint-pair fixture rows — single
        // ring (no pairs) and disjoint rings (|i^j| = 0) both yield 0; the
        // |i^j| = 1 class is the spiro fixture above; |i^j| >= 2 classes
        // are naphthalene/norbornane/cubane above.
    }

    #[test]
    fn descriptor_n13_bridgehead_atoms_overlap_order_and_dedup() {
        // Source Lipinski.cpp:465-503: NumBridgeheadAtomsVersion "2.0.0";
        // calcNumBridgeheadAtoms ensures SSSR-or-better ring info, then for
        // every unordered bondRings row pair sharing MORE THAN ONE bond,
        // counts endpoint incidences of the SHARED bonds in a FRESH
        // per-pair array; atoms with incidence EXACTLY ONE (the chain-end
        // atoms of the shared-bond chain) are bridgeheads, deduped into the
        // optional out vector whose length is returned. Fused pairs sharing
        // a single bond NEVER qualify.
        assert_eq!(NUM_BRIDGEHEAD_ATOMS_VERSION, "2.0.0");

        // Fixed contrasts, cold + prepared (supplied SSSR-or-better rows).
        const CASES: &[(&str, u32)] = &[
            ("C1CC2CCC1C2", 2),        // norbornane: two rings share 2 bonds
            ("c1ccc2ccccc2c1", 0),     // naphthalene: single shared bond
            ("C1CCC2(C1)CCCCC2", 0),   // spiro: 0 shared bonds (1 atom)
            ("C12C3C4C1C5C2C3C45", 0), // cubane: faces share 0 or 1 bond
            ("C1CCNCC1CC1CCCCC1", 0),  // disjoint rings
            ("C1CCCCC1", 0),           // single ring: no pairs
            ("", 0),
            ("CCCCCC", 0),
        ];
        for (smiles, expected) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_bridgehead_atoms(&topology).unwrap(),
                *expected,
                "cold {smiles}"
            );
            let (assignment, ring_info, coordinates, properties) = prepared_input_for(&topology);
            let input = DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                &assignment,
                &ring_info,
            );
            assert_eq!(
                rings::num_bridgehead_atoms_prepared(&input).unwrap(),
                *expected,
                "prepared {smiles}"
            );
        }

        // SSSR ACQUISITION, weaker-than-SSSR arm: a prepared input with an
        // OtherOrUnknown (empty, not sssr-or-better) RingInfo must trigger
        // the ONE canonical recompute and still find the two norbornane
        // bridgeheads. The absent arm is the cold path; the sssr-or-better
        // arm is the normal prepared path above.
        let norbornane = smiles_topology("C1CC2CCC1C2");
        let weak_rings = cosmolkit_core::RingInfo::new(
            cosmolkit_core::RingFindType::OtherOrUnknown,
            norbornane.atoms.len(),
            norbornane.bonds.len(),
        );
        assert!(
            !weak_rings.is_sssr_or_better(),
            "fixture trigger: weaker-than-SSSR ring info"
        );
        let assignment = valence(&norbornane, "descriptor_n13_weak").unwrap();
        let weak_coordinates = cosmolkit_model::CoordinateBlock::default();
        let weak_properties = cosmolkit_model::MoleculeProperties::default();
        let weak_input = DescriptorInput::new(
            &norbornane,
            &weak_coordinates,
            &weak_properties,
            &assignment,
            &weak_rings,
        );
        assert_eq!(
            rings::num_bridgehead_atoms_prepared(&weak_input).unwrap(),
            2,
            "weaker-than-SSSR acquisition recomputes and finds both bridgeheads"
        );

        // EXACT ordered unique atom IDs: derive the expected pair from the
        // row data — the qualifying pair's shared bonds form a chain; the
        // two atoms touching exactly ONE shared bond are the bridgeheads,
        // emitted in ASCENDING atom index within the pair (the source's
        // atomCounts scan order).
        let base_rings = ring_info(&norbornane, "descriptor_n13_ids").unwrap();
        let rows = base_rings.bond_rings();
        assert_eq!(rows.len(), 2, "norbornane: exactly two rows");
        let shared: Vec<_> = rows[0]
            .iter()
            .filter(|bond_id| rows[1].contains(bond_id))
            .collect();
        assert!(shared.len() > 1, "norbornane rows share more than one bond");
        let mut incidence = vec![0usize; norbornane.atoms.len()];
        for bond_id in &shared {
            let bond = &norbornane.bonds[bond_id.index()];
            incidence[bond.begin().index()] += 1;
            incidence[bond.end().index()] += 1;
        }
        let mut expected: Vec<AtomId> = incidence
            .iter()
            .enumerate()
            .filter(|(_, count)| **count == 1)
            .map(|(index, _)| AtomId::new(index))
            .collect();
        expected.sort_by_key(|atom_id| atom_id.index());
        assert_eq!(expected.len(), 2, "exactly two chain-end atoms");
        let before = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        let ids = bridgehead_atom_ids(&norbornane).unwrap();
        let after = VALENCE_HELPER_ENTRIES.with(|count| count.get());
        assert_eq!(
            after - before,
            0,
            "cold bridgehead_atom_ids: no valence helper"
        );
        assert_eq!(ids, expected, "exact bridgehead identities and order");

        // PREPOPULATED output state (source out-parameter semantics at the
        // kernel boundary): preexisting entries are kept, deduped against,
        // and INCLUDED in the returned length.
        let mut preseeded_one = vec![expected[0]];
        assert_eq!(
            rings::bridgehead_atom_ids_kernel(
                &base_rings,
                &norbornane.bonds,
                norbornane.atoms.len(),
                &mut preseeded_one
            )
            .unwrap(),
            2,
            "one preexisting bridgehead: kept, other appended, count 2"
        );
        assert_eq!(preseeded_one.len(), 2);
        let mut preseeded_both = expected.clone();
        assert_eq!(
            rings::bridgehead_atom_ids_kernel(
                &base_rings,
                &norbornane.bonds,
                norbornane.atoms.len(),
                &mut preseeded_both
            )
            .unwrap(),
            2,
            "both preexisting bridgeheads: no duplicates pushed"
        );
        assert_eq!(preseeded_both.len(), 2);
    }

    #[test]
    fn descriptor_s01_stereo_count_assigned_state_dispatch() {
        // Source Lipinski.cpp:503-553: hasStereoAssigned is a PRESENCE
        // check on the molecule-level _StereochemDone property (any value,
        // including "0", counts as assigned); ABSENT triggers ONE
        // copy-and-assign with cleanIt=true, force=true,
        // flagPossible=true via the legacy core owner. Counts read the
        // _ChiralityPossible atom computed prop; the unspecified variant
        // adds chiral_tag() == Unspecified. NOTE: the cold path uses the
        // valence helper BY SOURCE NECESSITY (assignStereochemistry
        // consumes the valence/rank rows), so the ring-family
        // zero-valence-counter brackets do NOT apply to this unit.
        assert_eq!(NUM_ATOM_STEREO_CENTERS_VERSION, "1.0.1");
        assert_eq!(NUM_UNSPECIFIED_ATOM_STEREO_CENTERS_VERSION, "1.0.1");

        // ABSENT arm, cold (bare topology carries no _StereochemDone):
        // the ONE legacy assignment runs with clean/force/flagPossible
        // all true. Frozen from source + legacy-owner semantics BEFORE
        // running: the specified-tag alpha carbon of C[C@H](N)C(=O)O is
        // a possible center whose tag is kept (1 possible, 0
        // unspecified); the untagged alpha carbon is possible AND
        // unspecified; ethanol and empty have none.
        const COLD: &[(&str, u32, u32)] = &[
            ("C[C@H](N)C(=O)O", 1, 0),
            ("C[CH](N)C(=O)O", 1, 1),
            ("CCO", 0, 0),
            ("", 0, 0),
        ];
        for (smiles, centers, unspecified) in COLD {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_atom_stereo_centers(&topology).unwrap(),
                *centers,
                "cold centers {smiles}"
            );
            assert_eq!(
                num_unspecified_atom_stereo_centers(&topology).unwrap(),
                *unspecified,
                "cold unspecified {smiles}"
            );
        }

        // PRESENCE arms (prepared): propane whose MIDDLE atom carries a
        // hand-set stale `_ChiralityPossible` computed prop. With
        // `_StereochemDone` present — value "0" AND value "1" — the
        // assigned arm reads the SUPPLIED atoms unchanged (1 possible,
        // 1 unspecified). With `_StereochemDone` ABSENT the assignment
        // runs with cleanIt=true: the middle carbon of propane is NOT a
        // possible center and the stale prop is cleared (0/0). The
        // stale-prop discriminator proves which arm executed.
        let mut stale = smiles_topology("CCC");
        stale.atoms[1]
            .set_computed_prop("_ChiralityPossible", "1")
            .expect("set stale chirality prop");
        let assignment = valence(&stale, "descriptor_s01").unwrap();
        let rings = ring_info(&stale, "descriptor_s01").unwrap();
        let coordinates = CoordinateBlock::default();
        for done_value in ["0", "1"] {
            let properties = MoleculeProperties::default()
                .with_prop("_StereochemDone", done_value)
                .expect("set _StereochemDone");
            let input =
                DescriptorInput::new(&stale, &coordinates, &properties, &assignment, &rings);
            assert_eq!(
                stereo::num_atom_stereo_centers_prepared(&input).unwrap(),
                1,
                "present-{done_value}: supplied state read unchanged"
            );
            assert_eq!(
                stereo::num_unspecified_atom_stereo_centers_prepared(&input).unwrap(),
                1,
                "present-{done_value}: supplied tag Unspecified"
            );
        }
        let absent_properties = MoleculeProperties::default();
        let absent_input = DescriptorInput::new(
            &stale,
            &coordinates,
            &absent_properties,
            &assignment,
            &rings,
        );
        assert_eq!(
            stereo::num_atom_stereo_centers_prepared(&absent_input).unwrap(),
            0,
            "absent: cleanIt=true clears the stale prop on a non-center"
        );
        assert_eq!(
            stereo::num_unspecified_atom_stereo_centers_prepared(&absent_input).unwrap(),
            0,
            "absent: no possible centers remain"
        );

        // Complete small parameter product over the dispatch axes:
        // {_StereochemDone absent, present-0, present-1} x {specified tag,
        // unspecified tag, no tag} — the cold COLD table covers the
        // no-marker rows for all three tag classes; the stale-prop block
        // covers all three marker states on a fixed tag class; together
        // the 3x3 product is enumerated with literal expectations.
    }

    #[test]
    fn descriptor_s02_chirality_property_count_tag_matrix() {
        // Source Lipinski.cpp:509-530: numAtomStereoCenters counts the
        // `_ChiralityPossible` PROPERTY ONLY — it is never a chiral-tag or
        // potential-stereo re-enumeration. Frozen matrix BEFORE running,
        // derived from the source dispatch + the legacy owner semantics:
        // LEGAL TAGGED center (tag kept by cleanIt=true, prop set) => 1;
        // UNTAGGED possible center => 1; ILLEGAL tag on a non-center =>
        // stripped by cleanIt on the absent arm, and on the PRESENT arm
        // the tag contributes NOTHING because the property is absent (the
        // tag-vs-property discriminator); stale properties on non-centers
        // count as-is on the PRESENT arm (property read, no enumeration).
        assert_eq!(NUM_ATOM_STEREO_CENTERS_VERSION, "1.0.1");

        // ABSENT-marker arm (cold; the ONE legacy assignment runs):
        // tagged-legal 1, untagged-possible 1, no-center 0, empty 0, and
        // two adjacent LEGAL centers (threonine) 2.
        const COLD: &[(&str, u32)] = &[
            ("C[C@H](N)C(=O)O", 1),         // tagged-legal alanine
            ("C[CH](N)C(=O)O", 1),          // untagged possible
            ("C[C@H](O)[C@H](N)C(=O)O", 2), // two tagged-legal centers
            ("CCO", 0),
            ("", 0),
        ];
        for (smiles, centers) in COLD {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_atom_stereo_centers(&topology).unwrap(),
                *centers,
                "cold {smiles}"
            );
        }

        // PRESENT-marker arm on an ILLEGALLY tagged propane: marker
        // present => supplied state read unchanged. The middle carbon
        // carries TetrahedralCw but NO _ChiralityPossible => 0 — a
        // tag-enumerating implementation would return 1 here.
        let mut illegal = smiles_topology("CCC");
        illegal.atoms[1].set_chiral_tag(cosmolkit_model::ChiralTag::TetrahedralCw);
        let assignment = valence(&illegal, "descriptor_s02").unwrap();
        let rings = ring_info(&illegal, "descriptor_s02").unwrap();
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default()
            .with_prop("_StereochemDone", "1")
            .expect("set _StereochemDone");
        let input = DescriptorInput::new(&illegal, &coordinates, &properties, &assignment, &rings);
        assert_eq!(
            stereo::num_atom_stereo_centers_prepared(&input).unwrap(),
            0,
            "present: tag without the property contributes nothing"
        );
        // ABSENT-marker arm on the SAME fixture: cleanIt=true strips the
        // illegal tag and force=true assigns no property => 0/0.
        let absent_properties = MoleculeProperties::default();
        let absent_input = DescriptorInput::new(
            &illegal,
            &coordinates,
            &absent_properties,
            &assignment,
            &rings,
        );
        assert_eq!(
            stereo::num_atom_stereo_centers_prepared(&absent_input).unwrap(),
            0,
            "absent: illegal tag stripped, no property set"
        );
        assert_eq!(
            stereo::num_unspecified_atom_stereo_centers_prepared(&absent_input).unwrap(),
            0,
            "absent: no possible centers remain"
        );

        // PROPERTY-not-enumeration discriminator, present arm: stale
        // `_ChiralityPossible` props on TWO non-center ethanol atoms count
        // as-is (2) with the marker present — the count reads properties,
        // it never re-derives potential stereo.
        let mut stale = smiles_topology("CCO");
        stale.atoms[0]
            .set_computed_prop("_ChiralityPossible", "1")
            .expect("stale prop 0");
        stale.atoms[2]
            .set_computed_prop("_ChiralityPossible", "1")
            .expect("stale prop 2");
        let stale_assignment = valence(&stale, "descriptor_s02_stale").unwrap();
        let stale_rings = ring_info(&stale, "descriptor_s02_stale").unwrap();
        let stale_input = DescriptorInput::new(
            &stale,
            &coordinates,
            &properties,
            &stale_assignment,
            &stale_rings,
        );
        assert_eq!(
            stereo::num_atom_stereo_centers_prepared(&stale_input).unwrap(),
            2,
            "present: stale properties on non-centers counted as-is"
        );
    }

    #[test]
    fn descriptor_n12_multiple_ids_order_product() {
        // N12-VERIFY frozen contract, fixture 1: the REAL disconnected
        // two-spiro topology C1CCC2(C1)CCCCC2.C1CCC2(C1)CCCCC2 (20
        // atoms, 4 perceived rings, spiro atoms 3 and 13 — independently
        // checked with RDKit 2026.3.1). All 24 row-index permutations x 6
        // literal preseeds = 144 REAL spiro_atom_ids_kernel calls. Pair
        // iteration is i<j; the intersection follows the first row's
        // member order; append only a single intersection member not
        // already in the CURRENT output; preexisting duplicate entries
        // are NOT cleaned up. Expected vectors are the frozen literals
        // below, never derived from the implementation.
        let topology = smiles_topology("C1CCC2(C1)CCCCC2.C1CCC2(C1)CCCCC2");
        assert_eq!(topology.atoms.len(), 20, "fixture: 20 atoms");
        let canonical = ring_info(&topology, "descriptor_n12_verify1").unwrap();
        let atom_rows = canonical.atom_rings();
        let bond_rows = canonical.bond_rings();
        assert_eq!(atom_rows.len(), 4, "fixture: 4 perceived rings");
        assert_eq!(bond_rows.len(), 4);

        // Membership prerequisites, proven from the canonical rows BEFORE
        // any kernel call: exactly two rows contain atom 3 (label them
        // A0/A1), exactly two contain atom 13 (B0/B1), each spiro atom is
        // shared by EXACTLY its two component rows, and no row contains
        // both spiro atoms.
        let spiro = |atom_index: usize| {
            atom_rows
                .iter()
                .enumerate()
                .filter(|(_, row)| row.contains(&AtomId::new(atom_index)))
                .map(|(row_index, _)| row_index)
                .collect::<Vec<_>>()
        };
        let a_rows = spiro(3);
        let b_rows = spiro(13);
        assert_eq!(a_rows.len(), 2, "atom 3 in exactly two rows (A0/A1)");
        assert_eq!(b_rows.len(), 2, "atom 13 in exactly two rows (B0/B1)");
        // Frozen N12-CLOSE canonical component-order prerequisite: the
        // A component rows come FIRST in the canonical owner's row list
        // (source input order), then the B component rows.
        assert_eq!(a_rows, vec![0, 1], "A0/A1 are canonical rows 0 and 1");
        assert_eq!(b_rows, vec![2, 3], "B0/B1 are canonical rows 2 and 3");
        for row in a_rows.iter().chain(b_rows.iter()) {
            let contains_both = atom_rows[*row].contains(&AtomId::new(3))
                && atom_rows[*row].contains(&AtomId::new(13));
            assert!(!contains_both, "no row contains both spiro atoms");
        }
        // A0/A1/B0/B1 are the canonical indices; keep the ACTUAL rows.
        let labels = [a_rows[0], a_rows[1], b_rows[0], b_rows[1]];

        // All 24 permutations of the four labeled rows, frozen literally
        // in this order: the first 12 enumerate A before B in pair order
        // (output [3,13]); the last 12 enumerate B first ([13,3]).
        #[rustfmt::skip]
        const PERMUTATIONS: [[usize; 4]; 24] = [
            [0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1],
            [0, 3, 1, 2], [0, 3, 2, 1], [1, 0, 2, 3], [1, 0, 3, 2],
            [1, 2, 0, 3], [1, 2, 3, 0], [1, 3, 0, 2], [1, 3, 2, 0],
            [2, 0, 1, 3], [2, 0, 3, 1], [2, 1, 0, 3], [2, 1, 3, 0],
            [2, 3, 0, 1], [2, 3, 1, 0], [3, 0, 1, 2], [3, 0, 2, 1],
            [3, 1, 0, 2], [3, 1, 2, 0], [3, 2, 0, 1], [3, 2, 1, 0],
        ];

        // Six literal preseed rows. For each preseed the A-first and
        // B-first expected vectors/counts are SEPARATE literals — never
        // computed by scanning intersections or calling a descriptor.
        struct Seed {
            preseed: Vec<usize>,
            a_first: Vec<usize>,
            b_first: Vec<usize>,
        }
        #[rustfmt::skip]
        let seeds: [Seed; 6] = [
            Seed { preseed: vec![], a_first: vec![3, 13], b_first: vec![13, 3] },
            Seed { preseed: vec![3], a_first: vec![3, 13], b_first: vec![3, 13] },
            Seed { preseed: vec![13], a_first: vec![13, 3], b_first: vec![13, 3] },
            Seed { preseed: vec![0], a_first: vec![0, 3, 13], b_first: vec![0, 13, 3] },
            Seed { preseed: vec![13, 3], a_first: vec![13, 3], b_first: vec![13, 3] },
            Seed { preseed: vec![3, 3], a_first: vec![3, 3, 13], b_first: vec![3, 3, 13] },
        ];

        let mut real_calls = 0usize;
        for (permutation_index, permutation) in PERMUTATIONS.iter().enumerate() {
            let a_first_output = permutation_index < 12;
            // Transport PAIRED perceived atom/bond rows under this
            // permutation via the existing core ring owner; the rows are
            // the real perceived rows, only their order changes.
            let permuted_atom_rows: Vec<Vec<AtomId>> = permutation
                .iter()
                .map(|label| atom_rows[labels[*label]].clone())
                .collect();
            let permuted_bond_rows: Vec<Vec<BondId>> = permutation
                .iter()
                .map(|label| bond_rows[labels[*label]].clone())
                .collect();
            let transported = cosmolkit_core::ring_info_from_selected_rows(
                topology.atoms.len(),
                topology.bonds.len(),
                &permuted_atom_rows,
                &permuted_bond_rows,
            )
            .expect("transport permuted rows");
            let snapshot: Vec<Vec<AtomId>> = transported.atom_rings().to_vec();
            let bond_snapshot: Vec<Vec<BondId>> = transported.bond_rings().to_vec();

            for seed in &seeds {
                let mut out: Vec<AtomId> = seed.preseed.iter().map(|i| AtomId::new(*i)).collect();
                let expected_indexes = if a_first_output {
                    &seed.a_first
                } else {
                    &seed.b_first
                };
                let expected: Vec<AtomId> =
                    expected_indexes.iter().map(|i| AtomId::new(*i)).collect();
                let expected_count = u32::try_from(expected.len()).expect("literal count");

                let returned =
                    rings::spiro_atom_ids_kernel(&transported, &mut out).expect("kernel call");
                assert_eq!(
                    returned, expected_count,
                    "perm {permutation_index} seed {:?}: return",
                    seed.preseed
                );
                assert_eq!(
                    out, expected,
                    "perm {permutation_index} seed {:?}: full vector",
                    seed.preseed
                );
                assert_eq!(
                    transported.atom_rings(),
                    &snapshot[..],
                    "perm {permutation_index} seed {:?}: rows unchanged",
                    seed.preseed
                );
                assert_eq!(
                    transported.bond_rings(),
                    &bond_snapshot[..],
                    "perm {permutation_index} seed {:?}: bond rows unchanged",
                    seed.preseed
                );
                real_calls += 1;
            }
        }
        assert_eq!(real_calls, 144, "24 permutations x 6 preseeds");

        // Supplementary real cold checks on the actual supplied-row order
        // (row-permutation coverage is the matrix above; this asserts the
        // canonical count/ID set separately).
        assert_eq!(num_spiro_atoms(&topology).unwrap(), 2);
        // RAW cold order asserted BEFORE any sorting: the actual
        // supplied-row order must yield [3, 13] directly (A component
        // enumerated first). A failure here is a source-order
        // discrepancy to report, not a reason to derive the expectation
        // from the output.
        let cold_ids = spiro_atom_ids(&topology).unwrap();
        assert_eq!(
            cold_ids,
            vec![AtomId::new(3), AtomId::new(13)],
            "raw cold IDs in actual supplied-row order"
        );
        // Supplementary sorted-set check (ID membership, order-free).
        let mut cold_ids = spiro_atom_ids(&topology).unwrap();
        cold_ids.sort_by_key(|atom_id| atom_id.index());
        assert_eq!(cold_ids, vec![AtomId::new(3), AtomId::new(13)]);
    }

    #[test]
    fn descriptor_n12_repeated_pair_hits_product() {
        // N12-VERIFY frozen contract, fixture 2: the REAL topology
        // C1CC2CCC1C23CCCC3 (11 atoms, 3 perceived rings; canonical
        // member sets {0,1,2,5,6}, {2,3,4,5,6}, {6,7,8,9,10} —
        // independently checked with RDKit 2026.3.1). The two norbornane
        // rows share {2,5,6} (THREE atoms -> not a spiro pair); the third
        // row intersects EACH norbornane row in EXACTLY atom 6 — TWO
        // distinct qualifying pair hits for the SAME atom 6, deduped to a
        // single output entry [6]. All 6 permutations x 5 literal
        // preseeds = 30 REAL spiro_atom_ids_kernel calls with full
        // vector/count assertions and unchanged-row checks in EVERY case.
        let topology = smiles_topology("C1CC2CCC1C23CCCC3");
        assert_eq!(topology.atoms.len(), 11, "fixture: 11 atoms");
        let canonical = ring_info(&topology, "descriptor_n12_verify2").unwrap();
        let atom_rows = canonical.atom_rings();
        let bond_rows = canonical.bond_rings();
        assert_eq!(atom_rows.len(), 3, "fixture: 3 perceived rings");
        assert_eq!(bond_rows.len(), 3);

        // Literal member-set prerequisites: match the perceived rows to
        // the frozen canonical sets (retaining the rows' real member
        // order), then prove the TWO distinct pair hits for 6 explicitly.
        let to_set = |row: &[AtomId]| -> Vec<usize> {
            let mut set: Vec<usize> = row.iter().map(|atom_id| atom_id.index()).collect();
            set.sort_unstable();
            set
        };
        let find_row = |set: &[usize]| -> usize {
            atom_rows
                .iter()
                .position(|row| to_set(row) == set)
                .expect("match canonical member set to a perceived row")
        };
        let row_a = find_row(&[0, 1, 2, 5, 6]);
        let row_b = find_row(&[2, 3, 4, 5, 6]);
        let row_c = find_row(&[6, 7, 8, 9, 10]);
        let a: &[AtomId] = &atom_rows[row_a];
        let b: &[AtomId] = &atom_rows[row_b];
        let c: &[AtomId] = &atom_rows[row_c];
        let shared_ab: Vec<usize> = a
            .iter()
            .filter(|atom_id| b.contains(atom_id))
            .map(|atom_id| atom_id.index())
            .collect();
        let mut shared_ab_set = shared_ab.clone();
        shared_ab_set.sort_unstable();
        assert_eq!(
            shared_ab.len(),
            3,
            "norbornane pair shares THREE atoms: not a spiro pair"
        );
        assert_eq!(
            shared_ab_set,
            vec![2, 5, 6],
            "norbornane pair shares the SET 2,5,6 (row-order preserved in \
             shared_ab: {shared_ab:?})"
        );
        let shared_ac: Vec<usize> = a
            .iter()
            .filter(|atom_id| c.contains(atom_id))
            .map(|atom_id| atom_id.index())
            .collect();
        let shared_bc: Vec<usize> = b
            .iter()
            .filter(|atom_id| c.contains(atom_id))
            .map(|atom_id| atom_id.index())
            .collect();
        assert_eq!(
            shared_ac,
            vec![6],
            "hit 1 of 2: row A x row C intersect in exactly atom 6"
        );
        assert_eq!(
            shared_bc,
            vec![6],
            "hit 2 of 2: row B x row C intersect in exactly atom 6"
        );

        // All 6 permutations of the three rows, frozen literally.
        #[rustfmt::skip]
        const PERMUTATIONS: [[usize; 3]; 6] = [
            [0, 1, 2], [0, 2, 1], [1, 0, 2],
            [1, 2, 0], [2, 0, 1], [2, 1, 0],
        ];
        let labels = [row_a, row_b, row_c];

        // Five literal preseed rows with frozen expected vectors/counts
        // (identical across permutations: the single spiro atom 6 with
        // preexisting duplicates preserved).
        let seeds: [(Vec<usize>, Vec<usize>); 5] = [
            (vec![], vec![6]),
            (vec![6], vec![6]),
            (vec![0], vec![0, 6]),
            (vec![6, 6], vec![6, 6]),
            (vec![0, 6], vec![0, 6]),
        ];

        let mut real_calls = 0usize;
        for permutation in &PERMUTATIONS {
            let permuted_atom_rows: Vec<Vec<AtomId>> = permutation
                .iter()
                .map(|label| atom_rows[labels[*label]].clone())
                .collect();
            let permuted_bond_rows: Vec<Vec<BondId>> = permutation
                .iter()
                .map(|label| bond_rows[labels[*label]].clone())
                .collect();
            let transported = cosmolkit_core::ring_info_from_selected_rows(
                topology.atoms.len(),
                topology.bonds.len(),
                &permuted_atom_rows,
                &permuted_bond_rows,
            )
            .expect("transport permuted rows");
            let snapshot: Vec<Vec<AtomId>> = transported.atom_rings().to_vec();
            let bond_snapshot: Vec<Vec<BondId>> = transported.bond_rings().to_vec();

            for (preseed, expected_indexes) in &seeds {
                let mut out: Vec<AtomId> = preseed.iter().map(|i| AtomId::new(*i)).collect();
                let expected: Vec<AtomId> =
                    expected_indexes.iter().map(|i| AtomId::new(*i)).collect();
                let expected_count = u32::try_from(expected.len()).expect("literal count");
                let returned =
                    rings::spiro_atom_ids_kernel(&transported, &mut out).expect("kernel call");
                assert_eq!(
                    returned, expected_count,
                    "perm {permutation:?} seed {preseed:?}: return"
                );
                assert_eq!(
                    out, expected,
                    "perm {permutation:?} seed {preseed:?}: full vector"
                );
                assert_eq!(
                    transported.atom_rings(),
                    &snapshot[..],
                    "perm {permutation:?} seed {preseed:?}: rows unchanged"
                );
                assert_eq!(
                    transported.bond_rings(),
                    &bond_snapshot[..],
                    "perm {permutation:?} seed {preseed:?}: bond rows unchanged"
                );
                real_calls += 1;
            }
        }
        assert_eq!(real_calls, 30, "6 permutations x 5 preseeds");

        // Supplementary real cold checks on this exact topology.
        assert_eq!(num_spiro_atoms(&topology).unwrap(), 1);
        assert_eq!(spiro_atom_ids(&topology).unwrap(), vec![AtomId::new(6)]);
    }

    #[test]
    fn descriptor_s03_unspecified_property_tag_product() {
        // Source Lipinski.cpp:530-553: numUnspecifiedAtomStereoCenters
        // shares the SAME assigned-state dispatch as numAtomStereoCenters
        // (presence-mode _StereochemDone guard; absent => one copy +
        // assignStereochemistry cleanIt=true force=true flagPossible=true)
        // and counts an atom iff it has _ChiralityPossible AND
        // chiral_tag() == Unspecified. Frozen all4 property x tag product
        // on prepared fixtures BEFORE running (literal expectations):
        // (possible, unspecified) => 1; (possible, specified) => 0;
        // (no-property, unspecified) => 0; (no-property, specified) => 0.
        assert_eq!(NUM_UNSPECIFIED_ATOM_STEREO_CENTERS_VERSION, "1.0.1");

        // ABSENT arm, cold: the ONE shared assignment runs. Tagged alanine
        // center => possible AND specified => 0; untagged possible center
        // => 1; threonine's two tagged centers => 0; none/empty => 0.
        const COLD: &[(&str, u32)] = &[
            ("C[C@H](N)C(=O)O", 0),
            ("C[CH](N)C(=O)O", 1),
            ("C[C@H](O)[C@H](N)C(=O)O", 0),
            ("CCO", 0),
            ("", 0),
        ];
        for (smiles, unspecified) in COLD {
            let topology = smiles_topology(smiles);
            assert_eq!(
                num_unspecified_atom_stereo_centers(&topology).unwrap(),
                *unspecified,
                "cold {smiles}"
            );
        }

        // PRESENT arm on real propane fixtures with _StereochemDone="1":
        // the supplied state is read unchanged (shared dispatch — no
        // assignment runs), so the all4 product is pinned by hand-set
        // state. Middle atom 1 carries each combination in turn.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default()
            .with_prop("_StereochemDone", "1")
            .expect("set _StereochemDone");
        // (possible, unspecified): property set, tag Unspecified => 1.
        let mut both = smiles_topology("CCC");
        both.atoms[1]
            .set_computed_prop("_ChiralityPossible", "1")
            .expect("set possible");
        // (possible, specified): property set, TetrahedralCw => 0.
        let mut possible_specified = both.clone();
        possible_specified.atoms[1].set_chiral_tag(cosmolkit_model::ChiralTag::TetrahedralCw);
        // (no-property, unspecified): no property, default tag => 0.
        let neither = smiles_topology("CCC");
        // (no-property, specified): no property, TetrahedralCw => 0 — the
        // tag alone never counts (the S02 discriminator on the shared
        // dispatch).
        let mut only_tag = smiles_topology("CCC");
        only_tag.atoms[1].set_chiral_tag(cosmolkit_model::ChiralTag::TetrahedralCw);

        let cases: [(&str, &TopologyBlock, u32); 4] = [
            ("possible+unspecified", &both, 1),
            ("possible+specified", &possible_specified, 0),
            ("no-property+unspecified", &neither, 0),
            ("no-property+specified", &only_tag, 0),
        ];
        for (label, topology, expected) in cases {
            let assignment = valence(topology, "descriptor_s03").unwrap();
            let rings = ring_info(topology, "descriptor_s03").unwrap();
            let input =
                DescriptorInput::new(topology, &coordinates, &properties, &assignment, &rings);
            assert_eq!(
                stereo::num_unspecified_atom_stereo_centers_prepared(&input).unwrap(),
                expected,
                "present arm {label}"
            );
        }
    }

    #[test]
    fn descriptor_d02_average_and_exact_weight_heavy_product() {
        // Pinned atomic_data.cpp constants (avg weight; most-common isotope
        // mass) and MolProps.cpp accumulation order, hand-derived per atom.
        const W_H: f64 = 1.008;
        const M_H: f64 = 1.007825032;
        const E: f64 = 0.00054857991;
        struct Case {
            smiles: &'static str,
            avg: f64,
            avg_heavy: f64,
            exact: f64,
            exact_heavy: f64,
        }
        let cases = [
            Case {
                smiles: "",
                avg: 0.0,
                avg_heavy: 0.0,
                exact: 0.0,
                exact_heavy: 0.0,
            },
            Case {
                smiles: "CCO",
                avg: 12.011 + 3.0 * W_H + 12.011 + 2.0 * W_H + 15.999 + 1.0 * W_H,
                avg_heavy: 12.011 + 12.011 + 15.999,
                exact: 12.0 + 12.0 + 15.99491462 + 6.0 * M_H,
                exact_heavy: 12.0 + 12.0 + 15.99491462,
            },
            Case {
                smiles: "[Na+].[Cl-]",
                avg: 22.99 + 35.453,
                avg_heavy: 22.99 + 35.453,
                exact: 22.98976928 - 1.0 * E + 34.96885268 + 1.0 * E,
                exact_heavy: 22.98976928 - 1.0 * E + 34.96885268 + 1.0 * E,
            },
            Case {
                smiles: "[NH4+]",
                avg: 14.007 + 4.0 * W_H,
                avg_heavy: 14.007,
                exact: 14.003074 - 1.0 * E + 4.0 * M_H,
                exact_heavy: 14.003074 - 1.0 * E,
            },
            Case {
                smiles: "c1cc[nH]c1",
                avg: 12.011
                    + 1.0 * W_H
                    + 12.011
                    + 1.0 * W_H
                    + 12.011
                    + 1.0 * W_H
                    + 14.007
                    + 1.0 * W_H
                    + 12.011
                    + 1.0 * W_H,
                avg_heavy: 12.011 + 12.011 + 12.011 + 14.007 + 12.011,
                exact: 12.0 + 12.0 + 12.0 + 14.003074 + 12.0 + 5.0 * M_H,
                exact_heavy: 12.0 + 12.0 + 12.0 + 14.003074 + 12.0,
            },
        ];
        for case in cases {
            let topology = smiles_topology(case.smiles);
            assert_eq!(
                molecular_weight(&topology).unwrap(),
                case.avg,
                "avg {}",
                case.smiles
            );
            assert_eq!(
                molecular_weight_with_options(&topology, true).unwrap(),
                case.avg_heavy,
                "avg-heavy {}",
                case.smiles
            );
            assert_eq!(
                exact_molecular_weight(&topology).unwrap(),
                case.exact,
                "exact {}",
                case.smiles
            );
            assert_eq!(
                exact_molecular_weight_with_options(&topology, true).unwrap(),
                case.exact_heavy,
                "exact-heavy {}",
                case.smiles
            );
        }
    }

    #[test]
    fn descriptor_d02_formula_option_product_and_special_orders() {
        let empty = smiles_topology("");
        for separate in [false, true] {
            for abbreviate in [false, true] {
                assert_eq!(
                    molecular_formula_with_options(&empty, separate, abbreviate).unwrap(),
                    "",
                    "empty {separate}/{abbreviate}"
                );
            }
        }

        // [13CH4]: the isotope key appears only under separateIsotopes;
        // abbreviate affects only D/T hydrogen atoms, so both values agree.
        let methane13 = smiles_topology("[13CH4]");
        assert_eq!(molecular_formula(&methane13).unwrap(), "CH4");
        assert_eq!(
            molecular_formula_with_options(&methane13, false, true).unwrap(),
            "CH4"
        );
        assert_eq!(
            molecular_formula_with_options(&methane13, true, false).unwrap(),
            "[13C]H4"
        );
        assert_eq!(
            molecular_formula_with_options(&methane13, true, true).unwrap(),
            "[13C]H4"
        );

        // [2H]C: D abbreviation only under separate x abbreviate; the
        // isotope-numbered H key sorts after plain H.
        let methane_d = smiles_topology("[2H]C");
        assert_eq!(molecular_formula(&methane_d).unwrap(), "CH4");
        assert_eq!(
            molecular_formula_with_options(&methane_d, true, false).unwrap(),
            "CH3[2H]"
        );
        assert_eq!(
            molecular_formula_with_options(&methane_d, true, true).unwrap(),
            "CH3D"
        );

        // Charges: +1/-1 print no magnitude, |q|>1 prints the number; Hill
        // order puts C, then H, then the remaining symbols alphabetically.
        assert_eq!(
            molecular_formula(&smiles_topology("[NH4+]")).unwrap(),
            "H4N+"
        );
        assert_eq!(
            molecular_formula(&smiles_topology("[Cl-]C")).unwrap(),
            "CH3Cl-"
        );
        assert_eq!(
            molecular_formula(&smiles_topology("[Mg+2]")).unwrap(),
            "Mg+2"
        );
        assert_eq!(molecular_formula(&smiles_topology("[O-2]")).unwrap(), "O-2");
        assert_eq!(molecular_formula(&smiles_topology("CCO")).unwrap(), "C2H6O");
        assert_eq!(
            molecular_formula(&smiles_topology("[Na+].[Cl-]")).unwrap(),
            "ClNa"
        );
    }

    #[test]
    fn descriptor_d02_explicit_hydrogen_topology_and_arithmetic_order() {
        const W_H: f64 = 1.008;
        const M_H: f64 = 1.007825032;
        let mut atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::N)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        atoms.extend(
            (2..4).map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::H))),
        );
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(2),
                BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
            ),
        ];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        // Atom order N, C, H, H with implicit H rows 0, 3, 0, 0: explicit H
        // atoms carry their own mass and never add further H counts.
        assert_eq!(
            molecular_weight(&topology).unwrap(),
            14.007 + 0.0 * W_H + 12.011 + 3.0 * W_H + 1.008 + 0.0 * W_H + 1.008
        );
        assert_eq!(
            exact_molecular_weight(&topology).unwrap(),
            14.003074 + 12.0 + 1.007825032 + 1.007825032 + 3.0 * M_H
        );
        assert_eq!(molecular_formula(&topology).unwrap(), "CH5N");
        assert_eq!(
            molecular_formula_with_options(&topology, true, true).unwrap(),
            "CH5N"
        );
    }

    #[test]
    fn descriptor_d02_heavy_only_reads_no_h_state_and_nonstrict_overvalence() {
        // Heavy-only kernels never prepare H state: an exotic hydrogen
        // isotope absent from the pinned table fails the H-sensitive kernels
        // at the periodic-table cause while onlyHeavy succeeds.
        const W_CL: f64 = 35.453;
        const M_CL: f64 = 34.96885268;
        let exotic = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::H).with_isotope(99)),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            )],
            vec![],
            vec![],
        )
        .unwrap();
        assert!(molecular_weight(&exotic).is_err());
        assert!(exact_molecular_weight(&exotic).is_err());
        assert_eq!(
            molecular_weight_with_options(&exotic, true).unwrap(),
            12.011
        );
        assert_eq!(
            exact_molecular_weight_with_options(&exotic, true).unwrap(),
            12.0
        );

        // Pinned non-strict implicit-valence semantics (calculateImplicitValence
        // `res = 0` branch): an over-valent C yields implicit H zero WITHOUT
        // an error, so the H-sensitive kernels equal the heavy-only sums.
        let mut atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))];
        atoms.extend(
            (1..6).map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::CL))),
        );
        let bonds = (0..5)
            .map(|index| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(0), AtomId::new(index + 1), BondOrder::Single),
                )
            })
            .collect::<Vec<_>>();
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        assert_eq!(molecular_weight(&topology).unwrap(), 12.011 + 5.0 * W_CL);
        assert_eq!(
            exact_molecular_weight(&topology).unwrap(),
            12.0 + 5.0 * M_CL
        );
    }

    // ---- T01 TPSA bond/H fact builder (MolSurf.cpp:103-333) ----

    fn t01_pair(
        first: Element,
        second: Element,
        order: BondOrder,
        aromatic: bool,
    ) -> TopologyBlock {
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(first)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(second)),
        ];
        let mut spec = BondSpec::new(AtomId::new(0), AtomId::new(1), order);
        if aromatic {
            spec = spec.with_aromatic(true);
        }
        let bonds = vec![Bond::from_spec(BondId::new(0), spec)];
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    fn t01_star(atom: Element, order: BondOrder, neighbors: usize) -> TopologyBlock {
        let mut atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(atom))];
        let mut bonds = Vec::new();
        for i in 0..neighbors {
            atoms.push(Atom::from_spec(
                AtomId::new(i + 1),
                AtomSpec::new(Element::C),
            ));
            bonds.push(Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(0), AtomId::new(i + 1), order),
            ));
        }
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    fn t01_zero_assignment(len: usize) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: vec![0; len],
            implicit_hydrogens: vec![0; len],
        }
    }

    fn t01_run(
        topology: &TopologyBlock,
        assignment: &ValenceAssignment,
        include_sand_p: bool,
    ) -> (f64, Vec<f64>) {
        let rings = cosmolkit_core::RingInfo::new(
            cosmolkit_core::RingFindType::OtherOrUnknown,
            topology.atoms.len(),
            topology.bonds.len(),
        );
        let mut contribs = vec![0.0f64; topology.atoms.len()];
        let res = tpsa::tpsa_atom_contribs_kernel(
            topology,
            assignment,
            &rings,
            include_sand_p,
            &mut contribs,
        )
        .expect("t01 fact kernel");
        (res, contribs)
    }

    #[test]
    fn descriptor_t01_explicit_h_precedence_product() {
        // Source MolSurf.cpp:120-135: a bond whose begin OR end atom is
        // hydrogen NEVER reaches the aromatic/order arms (begin==H =>
        // nNbrs[end] -= 1, nHs[end] += 1; else end==H symmetric). Frozen
        // literal expectations (zero implicit Hs, no 3-ring, sandp=false):
        // every H-edge cell has net nNbrs=0 and nHs=1 on the heavy atom,
        // so N falls to 30.5 - 0*8.2 + 1*1.5 = 32.0 and O to
        // 28.5 - 0*8.6 + 1*1.5 = 30.0. If precedence were broken, the
        // aromatic arm would give N 30.5-8.2 = 22.3 / O 19.9 and the
        // DOUBLE arm would hit O 17.07 (nHs==0, chg==0, nDoub==1).
        // Each row's last field is the HEAVY atom index (H-first rows
        // carry the heavy atom at index 1).
        for (first, second, order, aromatic, expected, heavy) in [
            (Element::H, Element::N, BondOrder::Single, false, 32.0, 1),
            (Element::H, Element::N, BondOrder::Aromatic, true, 32.0, 1),
            (Element::N, Element::H, BondOrder::Aromatic, true, 32.0, 0),
            (Element::H, Element::O, BondOrder::Single, false, 30.0, 1),
            (Element::H, Element::O, BondOrder::Double, false, 30.0, 1),
            (Element::O, Element::H, BondOrder::Double, false, 30.0, 0),
            (Element::H, Element::O, BondOrder::Aromatic, true, 30.0, 1),
        ] {
            let topology = t01_pair(first, second, order, aromatic);
            let assignment = t01_zero_assignment(2);
            let (res, contribs) = t01_run(&topology, &assignment, false);
            assert!(
                (res - expected).abs() < 1e-9,
                "{first:?}-{second:?} {order:?} arom={aromatic}: {res} != {expected}"
            );
            let mut want = [0.0f64, 0.0f64];
            want[heavy] = expected;
            for (got, expected_contrib) in contribs.iter().zip(want) {
                assert!(
                    (got - expected_contrib).abs() < 1e-9,
                    "contribs {contribs:?} != {want:?}"
                );
            }
        }
    }

    #[test]
    fn descriptor_t01_bond_class_product_and_dative_counts_nothing() {
        // Source MolSurf.cpp:135-149: non-H bonds classify by isAromatic
        // then SINGLE/DOUBLE/TRIPLE; the switch default counts NOTHING.
        // Frozen [O, C] cells with one implicit H on O: Single hits the
        // nNbrs==1 arm (nHs==1 && chg==0 && nSing==1) => 20.23; Double,
        // Triple, Aromatic and Dative all fall to the fallback
        // 28.5 - 1*8.6 + 1*1.5 = 21.4 (Dative reaches the default arm;
        // the adjacency still carries the edge so degree stays 1, the
        // typed atom->getDegree() projection).
        for (order, aromatic, expected, exact) in [
            (BondOrder::Single, false, 20.23, true),
            (BondOrder::Double, false, 21.4, false),
            (BondOrder::Triple, false, 21.4, false),
            (BondOrder::Aromatic, true, 21.4, false),
            (BondOrder::Dative, false, 21.4, false),
        ] {
            let topology = t01_pair(Element::O, Element::C, order, aromatic);
            let mut assignment = t01_zero_assignment(2);
            assignment.implicit_hydrogens[0] = 1;
            let (res, contribs) = t01_run(&topology, &assignment, false);
            if exact {
                assert_eq!(res, expected, "{order:?}");
                assert_eq!(contribs[0], expected);
            } else {
                assert!(
                    (res - expected).abs() < 1e-9,
                    "{order:?}: {res} != {expected}"
                );
                assert!((contribs[0] - expected).abs() < 1e-9);
            }
            assert_eq!(contribs[1], 0.0, "carbon filtered");
        }
    }

    #[test]
    fn descriptor_t01_explicit_h_no_double_count_product() {
        // Source MolSurf.cpp:123-128 + 160: the bond pass tallies H-atom
        // neighbors into nHs, then the preamble adds
        // atom->getTotalNumHs() with includeNeighbors=FALSE — implicit +
        // atom-spec explicit Hs only, never the explicit neighbors again.
        // Frozen methylammonium-shaped facts: N(+1) with three explicit
        // H atom bonds and one C single, zero implicit => nHs=3, net
        // nNbrs=1 (degree 4 minus 3), nSing=1, chg=1 => literal arm
        // 27.64. Wrong include_neighbors=TRUE would give nHs=6 =>
        // fallback 30.5-8.2+9.0 = 31.3; skipping the bond pass would
        // give nNbrs=4/nSing=4 => the chg==1 case-4 arm 0.0.
        let mut atoms = vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::N).with_formal_charge(1),
        )];
        let mut bonds = Vec::new();
        for i in 0..3 {
            atoms.push(Atom::from_spec(
                AtomId::new(i + 1),
                AtomSpec::new(Element::H),
            ));
            bonds.push(Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(0), AtomId::new(i + 1), BondOrder::Single),
            ));
        }
        atoms.push(Atom::from_spec(AtomId::new(4), AtomSpec::new(Element::C)));
        bonds.push(Bond::from_spec(
            BondId::new(3),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ));
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(5);
        let (res, contribs) = t01_run(&topology, &assignment, false);
        assert_eq!(res, 27.64);
        assert_eq!(contribs, vec![27.64, 0.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn descriptor_t01_atom_spec_explicit_hydrogens_counted() {
        // Source MolSurf.cpp:160 (nHs[i] += atom->getTotalNumHs()):
        // getTotalNumHs includes the atom-spec explicit H count
        // (d_numExplicitHs) with no H atom nodes. Frozen: N with
        // with_explicit_hydrogens(2), one C single, zero implicit =>
        // nHs=2, nNbrs=1, nSing=1, chg=0 => literal arm 26.02; an
        // implicit-only reading would give nHs=0 => fallback 22.3.
        let atoms = vec![
            Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::N).with_explicit_hydrogens(2),
            ),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(2);
        let (res, contribs) = t01_run(&topology, &assignment, false);
        assert_eq!(res, 26.02);
        assert_eq!(contribs, vec![26.02, 0.0]);
    }

    #[test]
    fn descriptor_t01_three_ring_facts_real_molecules() {
        // Source MolSurf.cpp:162: in3Ring = isAtomInRingOfSize(i, 3)
        // feeds the in3Ring/!in3Ring arm pairs. Frozen cold-wrapper
        // literals (valence owner for implicit Hs, SSSR owner for the
        // 3-ring rows): aziridine N (nHs=1, nSing=2, in3Ring) 21.94 vs
        // dimethylamine 12.03; oxirane O (nSing=2, in3Ring) 12.53 vs
        // dimethyl ether 9.23; 1-methylaziridine N (nSing=3, nHs=0,
        // in3Ring) 3.01 vs trimethylamine 3.24; ethanol O arm 20.23;
        // acetonitrile N arm (nHs=0, chg=0, nTrip=1) 23.79.
        for (smiles, expected) in [
            ("C1CN1", 21.94),
            ("CNC", 12.03),
            ("C1CO1", 12.53),
            ("COC", 9.23),
            ("C1CN1C", 3.01),
            ("CN(C)C", 3.24),
            ("CCO", 20.23),
            ("CC#N", 23.79),
        ] {
            let topology = smiles_topology(smiles);
            assert_eq!(
                tpsa_from_topology(&topology, false).unwrap(),
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn descriptor_t01_include_sand_p_product() {
        // Source MolSurf.cpp:154-156: the element filter admits N/O
        // always and P/S only when includeSandP; P/S tables have NO
        // fallback (tmp starts 0.0, MolSurf.cpp:278/298). Frozen
        // hand-built cells (zero implicit): P with three C singles =>
        // filtered 0.0 (false) vs literal arm nHs==0 && chg==0 &&
        // nSing==3 => 13.59 (true); S with one H bond + one C single =>
        // 0.0 (false) vs nHs==1 && chg==0 && nSing==1 => 38.80 (true).
        let phosphine = t01_star(Element::P, BondOrder::Single, 3);
        let assignment = t01_zero_assignment(4);
        let (res, contribs) = t01_run(&phosphine, &assignment, false);
        assert_eq!((res, contribs.clone()), (0.0, vec![0.0; 4]));
        let (res, contribs) = t01_run(&phosphine, &assignment, true);
        assert_eq!(res, 13.59);
        assert_eq!(contribs, vec![13.59, 0.0, 0.0, 0.0]);

        let thiol = t01_pair(Element::S, Element::H, BondOrder::Single, false);
        let mut atoms = thiol.atoms.clone();
        atoms.push(Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::C)));
        let mut bonds = thiol.bonds.clone();
        bonds.push(Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
        ));
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(3);
        let (res, contribs) = t01_run(&topology, &assignment, false);
        assert_eq!((res, contribs.clone()), (0.0, vec![0.0; 3]));
        let (res, contribs) = t01_run(&topology, &assignment, true);
        assert_eq!(res, 38.80);
        assert_eq!(contribs, vec![38.80, 0.0, 0.0]);
    }

    #[test]
    fn descriptor_t01_fallback_clamp_and_prepared_consistency() {
        // Source MolSurf.cpp:234-241 / 279-286: the N/O linear fallbacks
        // clamp below zero. Frozen cells (zero implicit, no 3-ring): N
        // with four C singles and chg=0 misses the case-4 arm (needs
        // chg==1) => 30.5-32.8 = -2.3 clamped to 0.0; O with four C
        // singles => 28.5-34.4 clamped to 0.0; the chg==1 N cell hits
        // the literal case-4 arm 0.0. Prepared/cold consistency on
        // dimethyl ether: kernel with one valence + one SSSR set equals
        // the cold wrapper and res equals the contribs sum.
        let neutral = t01_star(Element::N, BondOrder::Single, 4);
        let assignment = t01_zero_assignment(5);
        let (res, contribs) = t01_run(&neutral, &assignment, false);
        assert_eq!((res, contribs.clone()), (0.0, vec![0.0; 5]));

        let oxygen = t01_star(Element::O, BondOrder::Single, 4);
        let assignment = t01_zero_assignment(5);
        let (res, _) = t01_run(&oxygen, &assignment, false);
        assert_eq!(res, 0.0);

        let mut charged = neutral.clone();
        charged.atoms[0] = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::N).with_formal_charge(1),
        );
        let assignment = t01_zero_assignment(5);
        let (res, contribs) = t01_run(&charged, &assignment, false);
        assert_eq!((res, contribs.clone()), (0.0, vec![0.0; 5]));

        let topology = smiles_topology("COC");
        let prepared_valence = valence(&topology, "t01").unwrap();
        let prepared_rings = ring_info(&topology, "t01").unwrap();
        let mut contribs = vec![0.0f64; topology.atoms.len()];
        let kernel = tpsa::tpsa_atom_contribs_kernel(
            &topology,
            &prepared_valence,
            &prepared_rings,
            false,
            &mut contribs,
        )
        .unwrap();
        assert_eq!(kernel, 9.23);
        assert_eq!(tpsa_from_topology(&topology, false).unwrap(), kernel);
        assert!((contribs.iter().sum::<f64>() - kernel).abs() < 1e-9);
    }

    // ---- T02 TPSA nitrogen table (MolSurf.cpp:168-251) ----

    /// N(0) + `h_bonds` H-atom Single bonds + one C neighbor per
    /// (order, aromatic) pair; zero-implicit assignment (net nNbrs =
    /// heavy bond count).
    fn t02_n(
        chg: i8,
        spec_h: u8,
        h_bonds: usize,
        bonds: &[(BondOrder, bool)],
    ) -> (TopologyBlock, ValenceAssignment) {
        let mut spec = AtomSpec::new(Element::N).with_explicit_hydrogens(spec_h);
        if chg != 0 {
            spec = spec.with_formal_charge(chg);
        }
        let mut atoms = vec![Atom::from_spec(AtomId::new(0), spec)];
        let mut block_bonds = Vec::new();
        for _ in 0..h_bonds {
            atoms.push(Atom::from_spec(
                AtomId::new(atoms.len()),
                AtomSpec::new(Element::H),
            ));
            block_bonds.push(Bond::from_spec(
                BondId::new(block_bonds.len()),
                BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(atoms.len() - 1),
                    BondOrder::Single,
                ),
            ));
        }
        for (order, aromatic) in bonds {
            atoms.push(Atom::from_spec(
                AtomId::new(atoms.len()),
                AtomSpec::new(Element::C),
            ));
            let mut bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(atoms.len() - 1), *order);
            if *aromatic {
                bond_spec = bond_spec.with_aromatic(true);
            }
            block_bonds.push(Bond::from_spec(BondId::new(block_bonds.len()), bond_spec));
        }
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &block_bonds),
            atoms,
            bonds: block_bonds,
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(topology.atoms.len());
        (topology, assignment)
    }

    #[test]
    fn descriptor_t02_real_molecule_arms() {
        // Cold tpsa() literals where the valence/SSSR owners produce
        // the arm facts deterministically (12 of the 26 arms).
        for (smiles, expected) in [
            ("CC#N", 23.79),         // nHs0 chg0 trip1
            ("C=N", 23.85),          // nHs1 chg0 doub1
            ("CN", 26.02),           // nHs2 chg0 sing1
            ("CN=C", 12.36),         // nHs0 chg0 sing1 doub1
            ("C1CN1", 21.94),        // nHs1 chg0 sing2 in3Ring
            ("CNC", 12.03),          // nHs1 chg0 sing2 !in3Ring
            ("C1CN1C", 3.01),        // nHs0 chg0 sing3 in3Ring
            ("CN(C)C", 3.24),        // nHs0 chg0 sing3 !in3Ring
            ("c1ccncc1", 12.89),     // nHs0 chg0 arom2
            ("c1cc[nH]c1", 15.79),   // nHs1 chg0 arom2
            ("Cn1cccc1", 4.93),      // nHs0 chg0 sing1 arom2
            ("[N+](C)(C)(C)C", 0.0), // case4: nHs0 sing4 chg1
        ] {
            assert_eq!(
                tpsa_from_topology(&smiles_topology(smiles), false).unwrap(),
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn descriptor_t02_hand_built_arm_and_fallback_product() {
        // Remaining 14 arms via hand-built fact topologies (zero-implicit
        // assignment + weak RingInfo isolates the fact builder), then
        // fallback/clamp cells. (chg, spec_h, h_bonds, bonds, expected,
        // exact) — literal source values.
        for (chg, spec_h, h_bonds, bonds, expected, exact) in [
            (
                1i8,
                2u8,
                0usize,
                &[(BondOrder::Double, false)][..],
                25.59f64,
                true,
            ),
            (1, 0, 3, &[(BondOrder::Single, false)], 27.64, true),
            (
                0,
                0,
                0,
                &[(BondOrder::Triple, false), (BondOrder::Double, false)],
                13.60,
                true,
            ),
            (
                1,
                0,
                0,
                &[(BondOrder::Triple, false), (BondOrder::Single, false)],
                4.36,
                true,
            ),
            (
                1,
                1,
                0,
                &[(BondOrder::Double, false), (BondOrder::Single, false)],
                13.97,
                true,
            ),
            (
                1,
                2,
                0,
                &[(BondOrder::Single, false), (BondOrder::Single, false)],
                16.61,
                true,
            ),
            (
                1,
                1,
                0,
                &[(BondOrder::Aromatic, true), (BondOrder::Aromatic, true)],
                14.14,
                true,
            ),
            (
                0,
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Double, false),
                    (BondOrder::Double, false),
                ],
                11.68,
                true,
            ),
            (
                1,
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Double, false),
                ],
                3.01,
                true,
            ),
            (
                1,
                1,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                4.44,
                true,
            ),
            (
                0,
                0,
                0,
                &[
                    (BondOrder::Aromatic, true),
                    (BondOrder::Aromatic, true),
                    (BondOrder::Aromatic, true),
                ],
                4.41,
                true,
            ),
            (
                0,
                0,
                0,
                &[
                    (BondOrder::Double, false),
                    (BondOrder::Aromatic, true),
                    (BondOrder::Aromatic, true),
                ],
                8.39,
                true,
            ),
            (
                1,
                0,
                0,
                &[
                    (BondOrder::Aromatic, true),
                    (BondOrder::Aromatic, true),
                    (BondOrder::Aromatic, true),
                ],
                4.10,
                true,
            ),
            (
                1,
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Aromatic, true),
                    (BondOrder::Aromatic, true),
                ],
                3.88,
                true,
            ),
            // fallback 30.5 - 8.2 + 0 (nNbrs=1), 30.5 - 16.4 (nNbrs=2),
            // clamped case nNbrs=4:
            (0, 0, 0, &[(BondOrder::Single, false)], 22.3, false),
            (
                0,
                0,
                0,
                &[(BondOrder::Single, false), (BondOrder::Single, false)],
                14.1,
                false,
            ),
            (
                0,
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                0.0,
                false,
            ),
        ] {
            let (topology, assignment) = t02_n(chg, spec_h, h_bonds, bonds);
            let (res, _) = t01_run(&topology, &assignment, false);
            if exact {
                assert_eq!(res, expected, "chg={chg} spec_h={spec_h} h={h_bonds}");
            } else {
                assert!(
                    (res - expected).abs() < 1e-9,
                    "chg={chg} spec_h={spec_h} h={h_bonds}: {res} != {expected}"
                );
            }
        }
    }

    // ---- T03 TPSA oxygen table (MolSurf.cpp:252-280) ----

    #[test]
    fn descriptor_t03_real_molecule_arms() {
        // All six literal oxygen branches via cold tpsa() (source
        // literals): case1 17.07 (carbonyl), 20.23 (alcohol OH), 23.06
        // (single-bond O-); case2 12.53 (oxirane sing2 in3Ring), 9.23
        // (ether sing2 !in3Ring), 13.14 (furan arom2). Acetate sums the
        // carbonyl + alkoxide arms: 17.07 + 23.06 = 40.13 (tolerance —
        // f64 sum of two table literals).
        for (smiles, expected, exact) in [
            ("C=O", 17.07, true),
            ("CO", 20.23, true),
            ("[O-]C(=O)C", 40.13, false),
            ("C1CO1", 12.53, true),
            ("COC", 9.23, true),
            ("c1ccoc1", 13.14, true),
        ] {
            let res = tpsa_from_topology(&smiles_topology(smiles), false).unwrap();
            if exact {
                assert_eq!(res, expected, "{smiles}");
            } else {
                assert!(
                    (res - expected).abs() < 1e-9,
                    "{smiles}: {res} != {expected}"
                );
            }
        }
    }

    #[test]
    fn descriptor_t03_fallback_and_clamp_cells() {
        // Source MolSurf.cpp:275-283: fallback 28.5 - nNbrs*8.6 + nHs*1.5
        // with the < 0 => 0.0 clamp. Hand-built fact cells (zero implicit,
        // weak RingInfo): O with one TRIPLE C (nNbrs=1, no matching arm)
        // => 19.9; bare O with one spec H (nNbrs=0) => 30.0; O with four
        // single C bonds => 28.5 - 34.4 clamped => 0.0.
        let triple = t01_pair(Element::O, Element::C, BondOrder::Triple, false);
        let assignment = t01_zero_assignment(2);
        let (res, _) = t01_run(&triple, &assignment, false);
        assert!((res - 19.9).abs() < 1e-9, "triple fallback: {res}");

        let atoms = vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::O).with_explicit_hydrogens(1),
        )];
        let bare = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(1);
        let (res, _) = t01_run(&bare, &assignment, false);
        assert!((res - 30.0).abs() < 1e-9, "nNbrs=0 fallback: {res}");

        let four = t01_star(Element::O, BondOrder::Single, 4);
        let assignment = t01_zero_assignment(5);
        let (res, _) = t01_run(&four, &assignment, false);
        assert_eq!(res, 0.0, "clamped fallback");
    }

    #[test]
    fn descriptor_t01_complete_explicit_h_fact_product() {
        // Frozen TPSA-FACT-CLOSE product: {N,O} x {H first,H second} x
        // {Single,Double,Triple,Aromatic,Dative} x is_aromatic {false,true}
        // x include_sand_p {false,true} = 80 REAL kernel calls. Zero
        // implicit and zero atom-spec explicit Hs, charge 0, empty
        // supplied ring rows, exactly one H-edge. Source H precedence
        // gives net degree 0 and nHs 1 on the heavy atom, so every cell
        // lands on the fallback: N 32.0 (0x4040000000000000), O 30.0
        // (0x403e000000000000); the H atom keeps its zero-filled +0.0
        // (0x0) entry. Bit constants frozen from source arithmetic
        // before execution, never from implementation output.
        const N_32_BITS: u64 = 0x4040_0000_0000_0000;
        const O_30_BITS: u64 = 0x403e_0000_0000_0000;
        const H_PLUS_ZERO_BITS: u64 = 0x0;
        let mut calls = 0u32;
        for heavy_el in [Element::N, Element::O] {
            for h_first in [false, true] {
                for order in [
                    BondOrder::Single,
                    BondOrder::Double,
                    BondOrder::Triple,
                    BondOrder::Aromatic,
                    BondOrder::Dative,
                ] {
                    for aromatic in [false, true] {
                        for include_sand_p in [false, true] {
                            let (first, second) = if h_first {
                                (Element::H, heavy_el)
                            } else {
                                (heavy_el, Element::H)
                            };
                            let topology = t01_pair(first, second, order, aromatic);
                            let heavy_idx = usize::from(h_first);
                            let h_idx = 1 - heavy_idx;
                            // Prerequisites BEFORE the call.
                            assert_eq!(topology.atoms[h_idx].element(), Element::H);
                            assert_eq!(topology.atoms[heavy_idx].element(), heavy_el);
                            assert_eq!(topology.bonds.len(), 1);
                            assert_eq!(topology.bonds[0].begin().index(), 0);
                            assert_eq!(topology.bonds[0].end().index(), 1);
                            assert_eq!(topology.bonds[0].order(), order);
                            assert_eq!(topology.bonds[0].is_aromatic(), aromatic);
                            assert_eq!(topology.atoms[h_idx].explicit_hydrogens(), 0);
                            assert_eq!(topology.atoms[heavy_idx].explicit_hydrogens(), 0);
                            assert_eq!(topology.atoms[heavy_idx].formal_charge(), 0);
                            let assignment = t01_zero_assignment(2);
                            let rings = cosmolkit_core::RingInfo::new(
                                cosmolkit_core::RingFindType::OtherOrUnknown,
                                2,
                                1,
                            );
                            let topology_snapshot = topology.clone();
                            let assignment_snapshot = assignment.clone();
                            let mut contribs = vec![0.0f64; 2];
                            let res = tpsa::tpsa_atom_contribs_kernel(
                                &topology,
                                &assignment,
                                &rings,
                                include_sand_p,
                                &mut contribs,
                            )
                            .expect("explicit-H kernel call");
                            let heavy_bits = if heavy_el == Element::N {
                                N_32_BITS
                            } else {
                                O_30_BITS
                            };
                            assert_eq!(contribs[heavy_idx].to_bits(), heavy_bits);
                            assert_eq!(contribs[h_idx].to_bits(), H_PLUS_ZERO_BITS);
                            assert_eq!(res.to_bits(), heavy_bits);
                            // Inputs unchanged by the call.
                            assert_eq!(topology, topology_snapshot);
                            assert_eq!(assignment, assignment_snapshot);
                            calls += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 80);
    }

    #[test]
    fn descriptor_t01_complete_non_h_fact_product() {
        // Frozen TPSA-FACT-CLOSE product: {O first,O second} x the SAME
        // five orders x independent aromatic flag {false,true} x
        // include_sand_p {false,true} = 40 REAL kernel calls. Other
        // endpoint C; exactly ONE implicit H on O; zero atom-spec H;
        // charge 0; empty supplied ring rows. ONLY Single+flag=false
        // hits the literal arm 20.23 (bits 0x40343ae147ae147b); every
        // other cell reaches the source ordered fallback
        // (28.5 - 8.6) + 1.5 (bits 0x4035666666666666). The carbon row
        // keeps its zero-filled +0.0 (bits 0x0). Bit constants frozen
        // from source arithmetic before execution, never from
        // implementation output.
        const O_20_23_BITS: u64 = 0x4034_3ae1_47ae_147b;
        const FALLBACK_21_4_BITS: u64 = 0x4035_6666_6666_6666;
        const C_PLUS_ZERO_BITS: u64 = 0x0;
        let mut calls = 0u32;
        let mut single_plain_cells = 0u32;
        for o_first in [false, true] {
            for order in [
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Triple,
                BondOrder::Aromatic,
                BondOrder::Dative,
            ] {
                for aromatic in [false, true] {
                    for include_sand_p in [false, true] {
                        let (first, second) = if o_first {
                            (Element::O, Element::C)
                        } else {
                            (Element::C, Element::O)
                        };
                        let topology = t01_pair(first, second, order, aromatic);
                        // O sits at index 0 when O-first, index 1 otherwise.
                        let o_idx = 1 - usize::from(o_first);
                        let c_idx = usize::from(o_first);
                        // Prerequisites BEFORE the call.
                        assert_eq!(topology.atoms[o_idx].element(), Element::O);
                        assert_eq!(topology.atoms[c_idx].element(), Element::C);
                        assert_eq!(topology.bonds.len(), 1);
                        assert_eq!(topology.bonds[0].begin().index(), 0);
                        assert_eq!(topology.bonds[0].end().index(), 1);
                        assert_eq!(topology.bonds[0].order(), order);
                        assert_eq!(topology.bonds[0].is_aromatic(), aromatic);
                        assert_eq!(topology.atoms[o_idx].explicit_hydrogens(), 0);
                        assert_eq!(topology.atoms[o_idx].formal_charge(), 0);
                        assert_eq!(topology.atoms[c_idx].formal_charge(), 0);
                        let mut assignment = t01_zero_assignment(2);
                        assignment.implicit_hydrogens[o_idx] = 1;
                        assignment.implicit_hydrogens[c_idx] = 0;
                        let rings = cosmolkit_core::RingInfo::new(
                            cosmolkit_core::RingFindType::OtherOrUnknown,
                            2,
                            1,
                        );
                        let topology_snapshot = topology.clone();
                        let assignment_snapshot = assignment.clone();
                        let mut contribs = vec![0.0f64; 2];
                        let res = tpsa::tpsa_atom_contribs_kernel(
                            &topology,
                            &assignment,
                            &rings,
                            include_sand_p,
                            &mut contribs,
                        )
                        .expect("non-H kernel call");
                        let expected_o_bits = if order == BondOrder::Single && !aromatic {
                            single_plain_cells += 1;
                            O_20_23_BITS
                        } else {
                            FALLBACK_21_4_BITS
                        };
                        assert_eq!(contribs[o_idx].to_bits(), expected_o_bits);
                        assert_eq!(contribs[c_idx].to_bits(), C_PLUS_ZERO_BITS);
                        assert_eq!(res.to_bits(), expected_o_bits);
                        // Inputs unchanged by the call.
                        assert_eq!(topology, topology_snapshot);
                        assert_eq!(assignment, assignment_snapshot);
                        calls += 1;
                    }
                }
            }
        }
        assert_eq!(calls, 40);
        assert_eq!(single_plain_cells, 4);
    }

    // ---- T04 TPSA phosphorus table (MolSurf.cpp:281-304) ----

    /// P(0) with optional charge/spec-H plus one C neighbor per
    /// (order, aromatic) pair; zero-implicit assignment.
    fn t04_p(
        chg: i8,
        spec_h: u8,
        bonds: &[(BondOrder, bool)],
    ) -> (TopologyBlock, ValenceAssignment) {
        let mut spec = AtomSpec::new(Element::P).with_explicit_hydrogens(spec_h);
        if chg != 0 {
            spec = spec.with_formal_charge(chg);
        }
        let mut atoms = vec![Atom::from_spec(AtomId::new(0), spec)];
        let mut block_bonds = Vec::new();
        for (order, aromatic) in bonds {
            atoms.push(Atom::from_spec(
                AtomId::new(atoms.len()),
                AtomSpec::new(Element::C),
            ));
            let mut bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(atoms.len() - 1), *order);
            if *aromatic {
                bond_spec = bond_spec.with_aromatic(true);
            }
            block_bonds.push(Bond::from_spec(BondId::new(block_bonds.len()), bond_spec));
        }
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &block_bonds),
            atoms,
            bonds: block_bonds,
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(topology.atoms.len());
        (topology, assignment)
    }

    #[test]
    fn descriptor_t04_phosphorus_arms_zero_fallback_and_sandp_product() {
        // Source MolSurf.cpp:281-304: FOUR literal arms (34.14 case2
        // sing1+doub1; 13.59 case3 sing3; 23.47 case3 nHs1 sing2+doub1;
        // 9.81 case4 sing3+doub1), tmp starts 0.0 with NO fallback —
        // every non-matching cell keeps the source ZERO. The whole arm
        // is guarded by includeSandP: P is filtered (contribution +0.0)
        // when false. Hand-built zero-implicit cells; expectations are
        // the source literals, never implementation output.
        for (chg, spec_h, bonds, sandp, expected) in [
            (
                0i8,
                0u8,
                &[(BondOrder::Single, false), (BondOrder::Double, false)][..],
                true,
                34.14f64,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                true,
                13.59,
            ),
            (
                0,
                1,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Double, false),
                ],
                true,
                23.47,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Double, false),
                ],
                true,
                9.81,
            ),
            // includeSandP = false filters P entirely: the 34.14 facts
            // contribute +0.0.
            (
                0,
                0,
                &[(BondOrder::Single, false), (BondOrder::Double, false)],
                false,
                0.0,
            ),
            // Source-zero fallback: non-matching facts keep tmp = 0.0
            // (P has NO linear fallback). sing2 nHs0 (case2 mismatch),
            // nNbrs=1 (no case 1), nNbrs=5 (default), nHs1+chg1 sing3
            // (case3 arms require chg == 0).
            (
                0,
                0,
                &[(BondOrder::Single, false), (BondOrder::Single, false)],
                true,
                0.0,
            ),
            (0, 0, &[(BondOrder::Single, false)], true, 0.0),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                true,
                0.0,
            ),
            (
                1,
                1,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                true,
                0.0,
            ),
        ] {
            let (topology, assignment) = t04_p(chg, spec_h, bonds);
            let (res, contribs) = t01_run(&topology, &assignment, sandp);
            assert_eq!(res, expected, "chg={chg} spec_h={spec_h} sandp={sandp}");
            assert_eq!(contribs[0], expected, "P row carries the contribution");
            for (idx, contrib) in contribs.iter().enumerate().skip(1) {
                assert_eq!(*contrib, 0.0, "C row {idx} filtered");
            }
        }
    }

    // ---- T05 TPSA sulfur table (MolSurf.cpp:305-337) ----

    #[test]
    fn descriptor_t05_sulfur_arms_zero_fallback_and_sandp_product() {
        // Source MolSurf.cpp:305-337: SEVEN literal arms (32.09 doub1;
        // 38.80 sing1+H; 25.30 sing2; 28.24 arom2; 21.70 arom2+doub1;
        // 19.21 sing2+doub1; 8.38 sing2+doub2), tmp starts 0.0 with NO
        // fallback — non-matching cells keep the source ZERO; the whole
        // arm is guarded by includeSandP (filtered +0.0 when false).
        // Hand-built zero-implicit S cells; expectations are the source
        // literals, never implementation output.
        for (chg, spec_h, bonds, sandp, expected) in [
            // Seven literal arms, include_sand_p = true.
            (0i8, 0u8, &[(BondOrder::Double, false)][..], true, 32.09f64),
            (0, 1, &[(BondOrder::Single, false)], true, 38.80),
            (
                0,
                0,
                &[(BondOrder::Single, false), (BondOrder::Single, false)],
                true,
                25.30,
            ),
            (
                0,
                0,
                &[(BondOrder::Aromatic, true), (BondOrder::Aromatic, true)],
                true,
                28.24,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Aromatic, true),
                    (BondOrder::Aromatic, true),
                    (BondOrder::Double, false),
                ],
                true,
                21.70,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Double, false),
                ],
                true,
                19.21,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Double, false),
                    (BondOrder::Double, false),
                ],
                true,
                8.38,
            ),
            // includeSandP = false filters S entirely (32.09 facts -> +0.0).
            (0, 0, &[(BondOrder::Double, false)], false, 0.0),
            // Source-zero fallback: doub1 with chg=1 (case1 mismatch) and
            // nNbrs=5 (default) keep tmp = 0.0 — S has NO linear fallback.
            (1, 0, &[(BondOrder::Double, false)], true, 0.0),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                true,
                0.0,
            ),
        ] {
            let (topology, assignment) = t04_p(chg, spec_h, bonds);
            // The T04 builder parametrizes only the element via the spec;
            // rebuild the S row by swapping atom 0's element.
            let mut sulfur = topology.clone();
            sulfur.atoms[0] = Atom::from_spec(AtomId::new(0), {
                let mut spec = AtomSpec::new(Element::S).with_explicit_hydrogens(spec_h);
                if chg != 0 {
                    spec = spec.with_formal_charge(chg);
                }
                spec
            });
            let (res, contribs) = t01_run(&sulfur, &assignment, sandp);
            assert_eq!(res, expected, "chg={chg} spec_h={spec_h} sandp={sandp}");
            assert_eq!(contribs[0], expected, "S row carries the contribution");
            for (idx, contrib) in contribs.iter().enumerate().skip(1) {
                assert_eq!(*contrib, 0.0, "C row {idx} filtered");
            }
        }
    }

    #[test]
    fn descriptor_t05_complete_flag_product() {
        // Frozen T05 clarification product: ELEVEN source fact shapes x
        // sandp{false,true} = 22 REAL kernel calls. Seven literal arms
        // (32.09 doub1; 38.80 sing1+H; 25.30 sing2; 28.24 arom2; 21.70
        // arom2+doub1; 19.21 sing2+doub1; 8.38 sing2+doub2) and four
        // true-zero cells (bare S; [S]; [SSS]; [SSSSS]). Central S,
        // otherwise C, implicit H 0, empty supplied ring rows, all
        // non-Aromatic flags false. sandp=false => every row +0.0;
        // sandp=true => the source literals. Bit assertions use the
        // source literal's own bits (no kernel-derived value); +0.0 is
        // 0x0. Fixture facts asserted before every call; topology and
        // assignment asserted unchanged after every call; calls == 22.
        const PLUS_ZERO_BITS: u64 = 0x0;
        let mut calls = 0u32;
        for (chg, spec_h, bonds, true_total) in [
            (0i8, 0u8, &[][..], 0.0f64),
            (0, 0, &[(BondOrder::Double, false)][..], 32.09),
            (0, 1, &[(BondOrder::Single, false)], 38.80),
            (
                0,
                0,
                &[(BondOrder::Single, false), (BondOrder::Single, false)],
                25.30,
            ),
            (
                0,
                0,
                &[(BondOrder::Aromatic, true), (BondOrder::Aromatic, true)],
                28.24,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Aromatic, true),
                    (BondOrder::Aromatic, true),
                    (BondOrder::Double, false),
                ],
                21.70,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Double, false),
                ],
                19.21,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Double, false),
                    (BondOrder::Double, false),
                ],
                8.38,
            ),
            (0, 0, &[(BondOrder::Single, false)], 0.0),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                0.0,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                0.0,
            ),
        ] {
            for include_sand_p in [false, true] {
                let (topology, assignment) = t04_p(chg, spec_h, bonds);
                let mut sulfur = topology;
                let mut spec = AtomSpec::new(Element::S).with_explicit_hydrogens(spec_h);
                if chg != 0 {
                    spec = spec.with_formal_charge(chg);
                }
                sulfur.atoms[0] = Atom::from_spec(AtomId::new(0), spec);
                // Fixture facts BEFORE the call.
                assert_eq!(sulfur.atoms[0].element(), Element::S);
                assert_eq!(sulfur.atoms[0].formal_charge(), chg);
                assert_eq!(sulfur.atoms[0].explicit_hydrogens(), spec_h);
                assert_eq!(sulfur.bonds.len(), bonds.len());
                for (bond, (order, aromatic)) in sulfur.bonds.iter().zip(bonds) {
                    assert_eq!(bond.order(), *order);
                    assert_eq!(bond.is_aromatic(), *aromatic);
                }
                for (idx, atom) in sulfur.atoms.iter().enumerate().skip(1) {
                    assert_eq!(atom.element(), Element::C, "row {idx}");
                    assert_eq!(atom.formal_charge(), 0);
                }
                assert_eq!(assignment.implicit_hydrogens, vec![0; sulfur.atoms.len()]);
                let rings = cosmolkit_core::RingInfo::new(
                    cosmolkit_core::RingFindType::OtherOrUnknown,
                    sulfur.atoms.len(),
                    sulfur.bonds.len(),
                );
                let topology_snapshot = sulfur.clone();
                let assignment_snapshot = assignment.clone();
                let mut contribs = vec![0.0f64; sulfur.atoms.len()];
                let res = tpsa::tpsa_atom_contribs_kernel(
                    &sulfur,
                    &assignment,
                    &rings,
                    include_sand_p,
                    &mut contribs,
                )
                .expect("t05 flag-product kernel call");
                let expected_bits = if include_sand_p {
                    true_total.to_bits()
                } else {
                    PLUS_ZERO_BITS
                };
                assert_eq!(res.to_bits(), expected_bits, "sandp={include_sand_p}");
                assert_eq!(contribs[0].to_bits(), expected_bits, "S row");
                for (idx, contrib) in contribs.iter().enumerate().skip(1) {
                    assert_eq!(contrib.to_bits(), PLUS_ZERO_BITS, "C row {idx}");
                }
                assert_eq!(sulfur, topology_snapshot);
                assert_eq!(assignment, assignment_snapshot);
                calls += 1;
            }
        }
        assert_eq!(calls, 22);
    }

    #[test]
    fn descriptor_t04_complete_guard_product() {
        // Frozen T04-GUARD product: the SAME eight source fact shapes x
        // sandp{false,true} = 16 REAL kernel calls; the old nine-call
        // test remains untouched as supplementary coverage. Rows
        // (charge, specH, heavy bonds; true total): [SD]=34.14,
        // [SSS]=13.59, [SSD]+H=23.47, [SSSD]=9.81, [SS]=0.0, [S]=0.0,
        // [SSSSS]=0.0, chg1+H+[SSS]=0.0. All aromatic flags false,
        // implicit H 0, central P otherwise C, empty supplied ring
        // rows. sandp=false => EVERY row +0.0; sandp=true => the
        // source literals. Bit assertions use the source literal's own
        // bits (no kernel-derived value); +0.0 = 0x0.
        const PLUS_ZERO_BITS: u64 = 0x0;
        let mut calls = 0u32;
        for (chg, spec_h, bonds, true_total) in [
            (
                0i8,
                0u8,
                &[(BondOrder::Single, false), (BondOrder::Double, false)][..],
                34.14f64,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                13.59,
            ),
            (
                0,
                1,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Double, false),
                ],
                23.47,
            ),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Double, false),
                ],
                9.81,
            ),
            (
                0,
                0,
                &[(BondOrder::Single, false), (BondOrder::Single, false)],
                0.0,
            ),
            (0, 0, &[(BondOrder::Single, false)], 0.0),
            (
                0,
                0,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                0.0,
            ),
            (
                1,
                1,
                &[
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                    (BondOrder::Single, false),
                ],
                0.0,
            ),
        ] {
            for include_sand_p in [false, true] {
                let (topology, assignment) = t04_p(chg, spec_h, bonds);
                // Fixture facts BEFORE the call.
                assert_eq!(topology.atoms[0].element(), Element::P);
                assert_eq!(topology.atoms[0].formal_charge(), chg);
                assert_eq!(topology.atoms[0].explicit_hydrogens(), spec_h);
                assert_eq!(topology.bonds.len(), bonds.len());
                for (bond, (order, aromatic)) in topology.bonds.iter().zip(bonds) {
                    assert_eq!(bond.order(), *order);
                    assert_eq!(bond.is_aromatic(), *aromatic);
                }
                for (idx, atom) in topology.atoms.iter().enumerate().skip(1) {
                    assert_eq!(atom.element(), Element::C, "row {idx}");
                    assert_eq!(atom.formal_charge(), 0);
                }
                assert_eq!(assignment.implicit_hydrogens, vec![0; topology.atoms.len()]);
                let rings = cosmolkit_core::RingInfo::new(
                    cosmolkit_core::RingFindType::OtherOrUnknown,
                    topology.atoms.len(),
                    topology.bonds.len(),
                );
                let topology_snapshot = topology.clone();
                let assignment_snapshot = assignment.clone();
                let mut contribs = vec![0.0f64; topology.atoms.len()];
                let res = tpsa::tpsa_atom_contribs_kernel(
                    &topology,
                    &assignment,
                    &rings,
                    include_sand_p,
                    &mut contribs,
                )
                .expect("t04 guard-product kernel call");
                let expected_bits = if include_sand_p {
                    true_total.to_bits()
                } else {
                    PLUS_ZERO_BITS
                };
                assert_eq!(res.to_bits(), expected_bits, "sandp={include_sand_p}");
                assert_eq!(contribs[0].to_bits(), expected_bits, "P row");
                for (idx, contrib) in contribs.iter().enumerate().skip(1) {
                    assert_eq!(contrib.to_bits(), PLUS_ZERO_BITS, "C row {idx}");
                }
                assert_eq!(topology, topology_snapshot);
                assert_eq!(assignment, assignment_snapshot);
                calls += 1;
            }
        }
        assert_eq!(calls, 16);
    }

    // ---- T06 TPSA ordered contributions composition (MolSurf.cpp:153-345) ----

    #[test]
    fn descriptor_t06_ordered_composition_graphs_and_sandp_product() {
        // Frozen T06 semantics: contributions accumulate STRICTLY in
        // atom-index order; skipped rows keep +0.0; the sandp guard
        // changes ONLY the S/P rows' admission. Bit assertions use the
        // source literal's own bits; +0.0 = 0x0.
        const PLUS_ZERO_BITS: u64 = 0x0;

        // (a) EMPTY graph: 0 atoms, 0 bonds -> res +0.0, no rows.
        let empty = TopologyBlock::default();
        let assignment = t01_zero_assignment(0);
        let (res, contribs) = t01_run(&empty, &assignment, false);
        assert_eq!(res.to_bits(), PLUS_ZERO_BITS);
        assert!(contribs.is_empty());
        let (res, contribs) = t01_run(&empty, &assignment, true);
        assert_eq!(res.to_bits(), PLUS_ZERO_BITS);
        assert!(contribs.is_empty());

        // (b) ORDINARY ethanol-shaped graph: rows [C +0.0, C +0.0,
        // O 20.23] in INDEX order; total = ordered sum (20.23 bits),
        // both sandp values (guard has no effect on N/O rows).
        let mut atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
        ];
        atoms.truncate(3);
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            ),
        ];
        let ordinary = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let mut assignment = t01_zero_assignment(3);
        assignment.implicit_hydrogens[2] = 1;
        let snapshot = (ordinary.clone(), assignment.clone());
        for sandp in [false, true] {
            let (res, contribs) = t01_run(&ordinary, &assignment, sandp);
            assert_eq!(contribs[0].to_bits(), PLUS_ZERO_BITS, "C row 0");
            assert_eq!(contribs[1].to_bits(), PLUS_ZERO_BITS, "C row 1");
            assert_eq!(contribs[2].to_bits(), 20.23f64.to_bits(), "O row 2");
            assert_eq!(res.to_bits(), 20.23f64.to_bits(), "sandp={sandp}");
        }
        assert_eq!(ordinary, snapshot.0);
        assert_eq!(assignment, snapshot.1);

        // (c)+(d) S and P graphs x sandp: guard flips ONLY the S/P row.
        // S: specH1 + one C single -> 38.80 vs +0.0. P: three C singles
        // -> 13.59 vs +0.0.
        let (s_topology, s_assignment) = t04_p(0, 1, &[(BondOrder::Single, false)]);
        let mut s_graph = s_topology;
        s_graph.atoms[0] = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::S).with_explicit_hydrogens(1),
        );
        let s_snapshot = (s_graph.clone(), s_assignment.clone());
        for sandp in [false, true] {
            let (res, contribs) = t01_run(&s_graph, &s_assignment, sandp);
            let expected = if sandp {
                38.80f64.to_bits()
            } else {
                PLUS_ZERO_BITS
            };
            assert_eq!(contribs[0].to_bits(), expected, "S row sandp={sandp}");
            assert_eq!(contribs[1].to_bits(), PLUS_ZERO_BITS, "C row");
            assert_eq!(res.to_bits(), expected);
        }
        assert_eq!(s_graph, s_snapshot.0);
        assert_eq!(s_assignment, s_snapshot.1);

        let (p_graph, p_assignment) = t04_p(
            0,
            0,
            &[
                (BondOrder::Single, false),
                (BondOrder::Single, false),
                (BondOrder::Single, false),
            ],
        );
        let p_snapshot = (p_graph.clone(), p_assignment.clone());
        for sandp in [false, true] {
            let (res, contribs) = t01_run(&p_graph, &p_assignment, sandp);
            let expected = if sandp {
                13.59f64.to_bits()
            } else {
                PLUS_ZERO_BITS
            };
            assert_eq!(contribs[0].to_bits(), expected, "P row sandp={sandp}");
            for (idx, contrib) in contribs.iter().enumerate().skip(1) {
                assert_eq!(contrib.to_bits(), PLUS_ZERO_BITS, "C row {idx}");
            }
            assert_eq!(res.to_bits(), expected);
        }
        assert_eq!(p_graph, p_snapshot.0);
        assert_eq!(p_assignment, p_snapshot.1);

        // T01-T05 caller closure: cold tpsa() == one prepared kernel call.
        let ethanol = smiles_topology("CCO");
        let prepared_valence = valence(&ethanol, "t06").unwrap();
        let prepared_rings = ring_info(&ethanol, "t06").unwrap();
        let mut contribs = vec![0.0f64; ethanol.atoms.len()];
        let kernel = tpsa::tpsa_atom_contribs_kernel(
            &ethanol,
            &prepared_valence,
            &prepared_rings,
            false,
            &mut contribs,
        )
        .unwrap();
        assert_eq!(kernel.to_bits(), 20.23f64.to_bits());
        assert_eq!(contribs[2].to_bits(), 20.23f64.to_bits());
        assert_eq!(
            tpsa_from_topology(&ethanol, false).unwrap().to_bits(),
            kernel.to_bits()
        );
    }

    #[test]
    fn descriptor_t06_three_active_components_order_product() {
        // Frozen T06-ORDER product: three ACTIVE components in six
        // permutations discriminate floating summation ORDER (the prior
        // graphs had at most one nonzero row). A = N(chg0,specH0)-Triple-C
        // -> 23.79; B = N(chg0,specH2)-Single-C -> 26.02; C =
        // O(chg0,specH1)-Single-C -> 20.23; every carbon +0.0. Eligible
        // atom of component k lands at row 2*k, its carbon at 2*k+1, with
        // fresh consistent IDs/edges/adjacency per permutation. Six
        // orders x two sandp x two routes (kernel + existing cold
        // wrapper) = 24 REAL chemistry calls. Frozen total bits from the
        // supervisor's independent C++ source-arithmetic diagnostic
        // (/tmp/ck-t06-sum-reference, exit0): 012=0x4051828f5c28f5c3,
        // 021=0x4051828f5c28f5c2, 102=0x4051828f5c28f5c3,
        // 120=0x4051828f5c28f5c2, 201=0x4051828f5c28f5c2,
        // 210=0x4051828f5c28f5c2. No sum(), no oracle.
        const ORDER_TOTAL_BITS: [[u64; 2]; 3] = [
            [0x4051_828f_5c28_f5c3, 0x4051_828f_5c28_f5c2], // 012, 021
            [0x4051_828f_5c28_f5c3, 0x4051_828f_5c28_f5c2], // 102, 120
            [0x4051_828f_5c28_f5c2, 0x4051_828f_5c28_f5c2], // 201, 210
        ];
        // (element, spec_h, bond order, literal contribution)
        let components: [(Element, u8, BondOrder, f64); 3] = [
            (Element::N, 0, BondOrder::Triple, 23.79),
            (Element::N, 2, BondOrder::Single, 26.02),
            (Element::O, 1, BondOrder::Single, 20.23),
        ];
        let orders: [[usize; 3]; 6] = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        let mut calls = 0u32;
        for (order_idx, order) in orders.iter().enumerate() {
            // Fresh consistent IDs/edges/adjacency for THIS order.
            let mut atoms = Vec::new();
            let mut bonds = Vec::new();
            for (k, &comp) in order.iter().enumerate() {
                let (element, spec_h, bond_order, _) = components[comp];
                atoms.push(Atom::from_spec(
                    AtomId::new(2 * k),
                    AtomSpec::new(element).with_explicit_hydrogens(spec_h),
                ));
                atoms.push(Atom::from_spec(
                    AtomId::new(2 * k + 1),
                    AtomSpec::new(Element::C),
                ));
                bonds.push(Bond::from_spec(
                    BondId::new(k),
                    BondSpec::new(AtomId::new(2 * k), AtomId::new(2 * k + 1), bond_order),
                ));
            }
            let topology = TopologyBlock {
                adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
                atoms,
                bonds,
                ..TopologyBlock::default()
            };
            // Prerequisites BEFORE any call: elements, specH, charge,
            // bond orders/flags, degree, supplied zero-implicit assignment.
            for (k, &comp) in order.iter().enumerate() {
                let (element, spec_h, bond_order, _) = components[comp];
                assert_eq!(topology.atoms[2 * k].element(), element, "row {k}");
                assert_eq!(topology.atoms[2 * k].explicit_hydrogens(), spec_h);
                assert_eq!(topology.atoms[2 * k].formal_charge(), 0);
                assert_eq!(topology.bonds[k].order(), bond_order);
                assert!(!topology.bonds[k].is_aromatic());
                assert_eq!(topology.bonds[k].begin().index(), 2 * k);
                assert_eq!(topology.bonds[k].end().index(), 2 * k + 1);
                assert_eq!(topology.adjacency.neighbors_of(2 * k).len(), 1);
                assert_eq!(topology.atoms[2 * k + 1].element(), Element::C);
            }
            let assignment = t01_zero_assignment(6);
            let rings =
                cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 6, 3);
            for sandp_idx in 0..2usize {
                let include_sand_p = sandp_idx == 1;
                let expected_bits = ORDER_TOTAL_BITS[order_idx / 2][order_idx % 2];
                let topology_snapshot = topology.clone();
                let assignment_snapshot = assignment.clone();
                // Route 1: kernel with zero-filled rows; assert the
                // ENTIRE output vector by literal contribution bits.
                let mut contribs = vec![0.0f64; 6];
                let res = tpsa::tpsa_atom_contribs_kernel(
                    &topology,
                    &assignment,
                    &rings,
                    include_sand_p,
                    &mut contribs,
                )
                .expect("order-product kernel call");
                for (k, &comp) in order.iter().enumerate() {
                    assert_eq!(
                        contribs[2 * k].to_bits(),
                        components[comp].3.to_bits(),
                        "order {order:?} row {}",
                        2 * k
                    );
                    assert_eq!(
                        contribs[2 * k + 1].to_bits(),
                        0x0,
                        "carbon row {}",
                        2 * k + 1
                    );
                }
                assert_eq!(res.to_bits(), expected_bits, "kernel order {order:?}");
                assert_eq!(topology, topology_snapshot);
                assert_eq!(assignment, assignment_snapshot);
                calls += 1;
                // Route 2: existing cold wrapper (scalar only).
                let scalar = tpsa_from_topology(&topology, include_sand_p)
                    .expect("order-product cold wrapper call");
                assert_eq!(scalar.to_bits(), expected_bits, "wrapper order {order:?}");
                assert_eq!(topology, topology_snapshot);
                calls += 1;
            }
        }
        assert_eq!(calls, 24);
    }

    // ---- T07-CACHE contribution guard regressions (MolSurf.cpp:103-118) ----

    #[test]
    fn descriptor_t07_contributions_presence_force_and_flag_product() {
        // Frozen table: tpsa_contributions x sandp{false,true} x
        // force{false,true} x all FOUR independent presence states
        // {neither, scalar-only, contributions-only, both}. Unforced
        // BOTH => stored rows (literal sentinel [7,8,9]) unchanged;
        // unforced contributions-only => typed MissingComputedScalar;
        // every other cell recomputes ethanol rows [0,0,20.23] (bits)
        // and publishes both values. Zero VALENCE_HELPER_ENTRIES
        // reentry (supplied FINAL valence is borrowed, never
        // recomputed). Other-slot preservation asserted throughout.
        let topology = smiles_topology("CCO");
        let assignment = valence(&topology, "t07c").unwrap();
        let rings = ring_info(&topology, "t07c").unwrap();
        let properties = MoleculeProperties::default();
        let coordinates = CoordinateBlock::default();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        const O_20_23_BITS: u64 = 0x4034_3ae1_47ae_147b;
        const SENTINEL_ROWS: [f64; 3] = [7.0, 8.0, 9.0];
        let seed_both = |state: &mut tpsa::DescriptorComputedState, sandp: bool| {
            let slot = state.slot_mut_for_tests(sandp);
            slot.scalar = Some(123.0);
            slot.contributions = Some(SENTINEL_ROWS.to_vec());
        };
        let seed_scalar_only = |state: &mut tpsa::DescriptorComputedState, sandp: bool| {
            state.slot_mut_for_tests(sandp).scalar = Some(456.0);
        };
        let seed_contribs_only = |state: &mut tpsa::DescriptorComputedState, sandp: bool| {
            state.slot_mut_for_tests(sandp).contributions = Some(SENTINEL_ROWS.to_vec());
        };
        let mut calls = 0u32;
        for sandp in [false, true] {
            for force in [false, true] {
                for presence in 0..4 {
                    let mut state = tpsa::DescriptorComputedState::default();
                    match presence {
                        1 => seed_scalar_only(&mut state, sandp),
                        2 => seed_contribs_only(&mut state, sandp),
                        3 => seed_both(&mut state, sandp),
                        _ => {}
                    }
                    // Seed the OTHER slot to prove preservation.
                    seed_both(&mut state, !sandp);
                    let state_snapshot = state.clone();
                    let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                    let result = tpsa::tpsa_contributions(&input, sandp, force, &mut state);
                    calls += 1;
                    assert_eq!(
                        VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                        helper_before,
                        "zero valence-helper reentry"
                    );
                    match (force, presence) {
                        (false, 3) => {
                            let rows = result.expect("warm both serves stored rows");
                            assert_eq!(rows, SENTINEL_ROWS.to_vec());
                            // State unchanged on the warm hit.
                            assert_eq!(state, state_snapshot);
                        }
                        (false, 2) => {
                            let err = result.expect_err("contributions-only errors");
                            assert!(matches!(
                                err,
                                DescriptorError::MissingComputedScalar {
                                    function: "tpsa_contributions",
                                    include_sulfur_phosphorus: s
                                } if s == sandp
                            ));
                            assert_eq!(state, state_snapshot);
                        }
                        _ => {
                            let rows = result.expect("recompute");
                            assert_eq!(rows[0].to_bits(), 0x0);
                            assert_eq!(rows[1].to_bits(), 0x0);
                            assert_eq!(rows[2].to_bits(), O_20_23_BITS);
                            let slot = state.slot_for_tests(sandp);
                            assert_eq!(slot.scalar.map(f64::to_bits), Some(O_20_23_BITS));
                            let cached = slot.contributions.as_ref().unwrap();
                            assert_eq!(cached[2].to_bits(), O_20_23_BITS);
                        }
                    }
                    // Other slot untouched in EVERY arm.
                    let other = state.slot_for_tests(!sandp);
                    assert_eq!(other.scalar, Some(123.0));
                    assert_eq!(other.contributions, Some(SENTINEL_ROWS.to_vec()));
                }
            }
        }
        assert_eq!(calls, 16);
    }

    #[test]
    fn descriptor_t07_contributions_invalid_valence_atomicity() {
        // Typed cold failure: malformed supplied FINAL valence rows keep
        // the known DescriptorError::InvalidValenceRows cause and leave
        // the computed state UNCHANGED (publish-after-success).
        let topology = smiles_topology("CCO");
        let good = valence(&topology, "t07v").unwrap();
        let mut malformed = good.clone();
        malformed.implicit_hydrogens.resize(2, 0);
        let rings = ring_info(&topology, "t07v").unwrap();
        let properties = MoleculeProperties::default();
        let coordinates = CoordinateBlock::default();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &malformed, &rings);
        let mut state = tpsa::DescriptorComputedState::default();
        let snapshot = state.clone();
        let err = tpsa::tpsa_contributions(&input, false, false, &mut state)
            .expect_err("malformed rows must fail typed");
        assert!(matches!(
            err,
            DescriptorError::InvalidValenceRows {
                function: "tpsa_contributions",
                ..
            }
        ));
        assert_eq!(state, snapshot, "failure leaves state unchanged");
    }

    // ---- T07-CACHE scalar guard regressions (MolSurf.cpp:347-360) ----

    #[test]
    fn descriptor_t07_scalar_presence_force_and_flag_product() {
        // tpsa (canonical scalar dispatch) x sandp{false,true} x
        // force{false,true} x all FOUR independent presence states.
        // Unforced scalar-only / both => the STORED sentinel scalar
        // (123.0 bits when sandp=false, 456.0 bits when sandp=true)
        // with NO recompute; unforced contributions-only => typed
        // MissingComputedScalar via the contribution owner; every
        // other cell recomputes the ethanol scalar 20.23 (bits).
        // Other-slot preservation and zero valence-helper reentry
        // asserted in every arm.
        let topology = smiles_topology("CCO");
        let assignment = valence(&topology, "t07s").unwrap();
        let rings = ring_info(&topology, "t07s").unwrap();
        let properties = MoleculeProperties::default();
        let coordinates = CoordinateBlock::default();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        const O_20_23_BITS: u64 = 0x4034_3ae1_47ae_147b;
        let sentinel = |sandp: bool| if sandp { 456.0f64 } else { 123.0f64 };
        let mut calls = 0u32;
        for sandp in [false, true] {
            for force in [false, true] {
                for presence in 0..4 {
                    // presence: 0 neither, 1 scalar-only,
                    // 2 contributions-only, 3 both.
                    let mut state = tpsa::DescriptorComputedState::default();
                    if presence == 1 || presence == 3 {
                        state.slot_mut_for_tests(sandp).scalar = Some(sentinel(sandp));
                    }
                    if presence == 2 || presence == 3 {
                        state.slot_mut_for_tests(sandp).contributions =
                            Some(vec![0.0, 0.0, sentinel(sandp)]);
                    }
                    // Other slot fully seeded for preservation.
                    let other = state.slot_mut_for_tests(!sandp);
                    other.scalar = Some(sentinel(!sandp));
                    other.contributions = Some(vec![1.0, 2.0, 3.0]);
                    let snapshot = state.clone();
                    let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                    let result = tpsa(&input, sandp, force, &mut state);
                    calls += 1;
                    assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), helper_before);
                    match (force, presence) {
                        (false, 1) | (false, 3) => {
                            let res = result.expect("warm scalar hit");
                            assert_eq!(res.to_bits(), sentinel(sandp).to_bits());
                            // Warm scalar leaves the state UNCHANGED.
                            assert_eq!(state, snapshot);
                        }
                        (false, 2) => {
                            let err = result.expect_err("contributions-only errors");
                            assert!(matches!(
                                err,
                                DescriptorError::MissingComputedScalar {
                                    include_sulfur_phosphorus: s,
                                    ..
                                } if s == sandp
                            ));
                            assert_eq!(state, snapshot);
                        }
                        _ => {
                            let res = result.expect("recompute");
                            assert_eq!(res.to_bits(), O_20_23_BITS);
                            let slot = state.slot_for_tests(sandp);
                            assert_eq!(slot.scalar.map(f64::to_bits), Some(O_20_23_BITS));
                        }
                    }
                    let untouched = state.slot_for_tests(!sandp);
                    assert_eq!(untouched.scalar, Some(sentinel(!sandp)));
                    assert_eq!(untouched.contributions, Some(vec![1.0, 2.0, 3.0]));
                }
            }
        }
        assert_eq!(calls, 16);
    }

    // ---- T07-CACHE composed product + lifecycle sequences ----

    #[test]
    fn descriptor_t07_cache_product_cold_warm_cloned() {
        // Frozen 24-call product: 2 API kinds x 2 sandp x 2 force x 3
        // initial states {empty, warm, cloned-warm}. Warm slots use the
        // literal sentinels 123.0(false)/456.0(true) and rows [7,8,9];
        // force=false returns stored results unchanged; force=true and
        // cold recompute the ethanol literals [0,0,20.23]/20.23.
        let topology = smiles_topology("CCO");
        let assignment = valence(&topology, "t07p").unwrap();
        let rings = ring_info(&topology, "t07p").unwrap();
        let properties = MoleculeProperties::default();
        let coordinates = CoordinateBlock::default();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        const O_20_23_BITS: u64 = 0x4034_3ae1_47ae_147b;
        let sentinel = |sandp: bool| if sandp { 456.0f64 } else { 123.0f64 };
        let make_warm = |sandp: bool| {
            let mut state = tpsa::DescriptorComputedState::default();
            let slot = state.slot_mut_for_tests(sandp);
            slot.scalar = Some(sentinel(sandp));
            slot.contributions = Some(vec![7.0, 8.0, 9.0]);
            state
        };
        // Frozen T07-MATRIX other-slot sentinels (789.0 / [4,5,6]) and
        // full-state checkpoints for EVERY call.
        let seed_other = |state: &mut tpsa::DescriptorComputedState, sandp: bool| {
            let other = state.slot_mut_for_tests(!sandp);
            other.scalar = Some(789.0);
            other.contributions = Some(vec![4.0, 5.0, 6.0]);
        };
        let mut calls = 0u32;
        for api in 0..2 {
            for sandp in [false, true] {
                for force in [false, true] {
                    for initial in 0..3 {
                        let mut state = match initial {
                            0 => tpsa::DescriptorComputedState::default(),
                            1 => make_warm(sandp),
                            // cloned-warm: Clone preserves the warm slot.
                            _ => make_warm(sandp).clone(),
                        };
                        seed_other(&mut state, sandp);
                        let topology_snapshot = topology.clone();
                        let assignment_snapshot = assignment.clone();
                        let properties_snapshot = properties.clone();
                        let rings_snapshot = rings.clone();
                        let state_snapshot = state.clone();
                        let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                        let warm_hit = !force && initial != 0;
                        match api {
                            0 => {
                                let rows =
                                    tpsa::tpsa_contributions(&input, sandp, force, &mut state)
                                        .expect("contributions");
                                if warm_hit {
                                    assert_eq!(rows, vec![7.0, 8.0, 9.0]);
                                } else {
                                    // Full literal row bits, not only [2].
                                    assert_eq!(rows.len(), 3);
                                    assert_eq!(rows[0].to_bits(), 0x0);
                                    assert_eq!(rows[1].to_bits(), 0x0);
                                    assert_eq!(rows[2].to_bits(), O_20_23_BITS);
                                }
                            }
                            _ => {
                                let res = tpsa(&input, sandp, force, &mut state).expect("scalar");
                                let expected = if warm_hit {
                                    sentinel(sandp).to_bits()
                                } else {
                                    O_20_23_BITS
                                };
                                assert_eq!(res.to_bits(), expected);
                            }
                        }
                        calls += 1;
                        assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), helper_before);
                        assert_eq!(topology, topology_snapshot);
                        assert_eq!(assignment, assignment_snapshot);
                        assert_eq!(properties, properties_snapshot);
                        assert_eq!(rings, rings_snapshot);
                        let other = state.slot_for_tests(!sandp);
                        assert_eq!(other.scalar.map(f64::to_bits), Some(789.0f64.to_bits()));
                        assert_eq!(
                            other
                                .contributions
                                .as_ref()
                                .map(|rows| rows.iter().map(|v| v.to_bits()).collect::<Vec<_>>()),
                            Some(
                                [4.0f64, 5.0, 6.0]
                                    .iter()
                                    .map(|v| v.to_bits())
                                    .collect::<Vec<_>>()
                            )
                        );
                        if warm_hit {
                            // Warm success (either API) preserves the
                            // ENTIRE state.
                            assert_eq!(state, state_snapshot);
                        } else {
                            // Cold/forced success: selected scalar bits
                            // and FULL cached row vector by bits.
                            let slot = state.slot_for_tests(sandp);
                            assert_eq!(slot.scalar.map(f64::to_bits), Some(O_20_23_BITS));
                            let cached = slot.contributions.as_ref().unwrap();
                            assert_eq!(cached.len(), 3);
                            assert_eq!(cached[0].to_bits(), 0x0);
                            assert_eq!(cached[1].to_bits(), 0x0);
                            assert_eq!(cached[2].to_bits(), O_20_23_BITS);
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 24);
    }

    #[test]
    fn descriptor_t07_cache_product_independent_presence_and_sequences() {
        // Frozen 32-call product: 4 presence states x 2 API kinds x 2
        // sandp x 2 force. Unforced scalar-only: scalar API hits,
        // contribution API recomputes. Unforced contributions-only:
        // BOTH APIs error MissingComputedScalar. Forced: all states
        // recompute. Followed by the ordered lifecycle sequences:
        // cross-key cold fill, warm clone, explicit clear, forced
        // overwrite, typed cold invalid-valence leaving state
        // unchanged.
        let topology = smiles_topology("CCO");
        let assignment = valence(&topology, "t07q").unwrap();
        let rings = ring_info(&topology, "t07q").unwrap();
        let properties = MoleculeProperties::default();
        let coordinates = CoordinateBlock::default();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        const O_20_23_BITS: u64 = 0x4034_3ae1_47ae_147b;
        // Frozen T07-MATRIX other-slot sentinels and full-state
        // checkpoints for EVERY call.
        let seed_other = |state: &mut tpsa::DescriptorComputedState, sandp: bool| {
            let other = state.slot_mut_for_tests(!sandp);
            other.scalar = Some(789.0);
            other.contributions = Some(vec![4.0, 5.0, 6.0]);
        };
        let mut calls = 0u32;
        for presence in 0..4 {
            for api in 0..2 {
                for sandp in [false, true] {
                    for force in [false, true] {
                        let mut state = tpsa::DescriptorComputedState::default();
                        if presence == 1 || presence == 3 {
                            state.slot_mut_for_tests(sandp).scalar = Some(123.0);
                        }
                        if presence == 2 || presence == 3 {
                            state.slot_mut_for_tests(sandp).contributions =
                                Some(vec![7.0, 8.0, 9.0]);
                        }
                        seed_other(&mut state, sandp);
                        let topology_snapshot = topology.clone();
                        let assignment_snapshot = assignment.clone();
                        let properties_snapshot = properties.clone();
                        let rings_snapshot = rings.clone();
                        let state_snapshot = state.clone();
                        let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                        match api {
                            0 => {
                                let result =
                                    tpsa::tpsa_contributions(&input, sandp, force, &mut state);
                                if !force && presence == 2 {
                                    assert!(matches!(
                                        result.unwrap_err(),
                                        DescriptorError::MissingComputedScalar { .. }
                                    ));
                                } else if !force && presence == 3 {
                                    assert_eq!(result.unwrap(), vec![7.0, 8.0, 9.0]);
                                } else {
                                    let rows = result.unwrap();
                                    assert_eq!(rows.len(), 3);
                                    assert_eq!(rows[0].to_bits(), 0x0);
                                    assert_eq!(rows[1].to_bits(), 0x0);
                                    assert_eq!(rows[2].to_bits(), O_20_23_BITS);
                                }
                            }
                            _ => {
                                let result = tpsa(&input, sandp, force, &mut state);
                                if !force && presence == 2 {
                                    assert!(matches!(
                                        result.unwrap_err(),
                                        DescriptorError::MissingComputedScalar { .. }
                                    ));
                                } else if !force && presence == 1 {
                                    assert_eq!(result.unwrap().to_bits(), 123.0f64.to_bits());
                                } else if !force && presence == 3 {
                                    // Source calcTPSA warm arm (MolSurf.cpp:
                                    // 351-355): !force && scalar-present
                                    // returns the STORED scalar regardless
                                    // of contributions presence.
                                    assert_eq!(result.unwrap().to_bits(), 123.0f64.to_bits());
                                } else {
                                    assert_eq!(result.unwrap().to_bits(), O_20_23_BITS);
                                }
                            }
                        }
                        calls += 1;
                        assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), helper_before);
                        assert_eq!(topology, topology_snapshot);
                        assert_eq!(assignment, assignment_snapshot);
                        assert_eq!(properties, properties_snapshot);
                        assert_eq!(rings, rings_snapshot);
                        let other = state.slot_for_tests(!sandp);
                        assert_eq!(other.scalar.map(f64::to_bits), Some(789.0f64.to_bits()));
                        assert_eq!(
                            other
                                .contributions
                                .as_ref()
                                .map(|rows| rows.iter().map(|v| v.to_bits()).collect::<Vec<_>>()),
                            Some(
                                [4.0f64, 5.0, 6.0]
                                    .iter()
                                    .map(|v| v.to_bits())
                                    .collect::<Vec<_>>()
                            )
                        );
                        // Full-state/result obligations: warm success and
                        // contributions-only error preserve the ENTIRE
                        // state; cold/forced success leaves selected scalar
                        // bits and FULL cached rows by bits; contribution
                        // return rows use full literal bits.
                        match api {
                            0 => {
                                if !force && (presence == 2 || presence == 3) {
                                    if presence == 2 {
                                        // Error arm preserves the ENTIRE state.
                                        assert_eq!(state, state_snapshot);
                                    } else {
                                        assert_eq!(state, state_snapshot);
                                    }
                                } else {
                                    let slot = state.slot_for_tests(sandp);
                                    assert_eq!(slot.scalar.map(f64::to_bits), Some(O_20_23_BITS));
                                    let cached = slot.contributions.as_ref().unwrap();
                                    assert_eq!(cached.len(), 3);
                                    assert_eq!(cached[0].to_bits(), 0x0);
                                    assert_eq!(cached[1].to_bits(), 0x0);
                                    assert_eq!(cached[2].to_bits(), O_20_23_BITS);
                                }
                            }
                            _ => {
                                if !force && (presence == 1 || presence == 3) {
                                    assert_eq!(state, state_snapshot);
                                } else if !force && presence == 2 {
                                    // Error arm preserves the ENTIRE state.
                                    assert_eq!(state, state_snapshot);
                                } else {
                                    let slot = state.slot_for_tests(sandp);
                                    assert_eq!(slot.scalar.map(f64::to_bits), Some(O_20_23_BITS));
                                    let cached = slot.contributions.as_ref().unwrap();
                                    assert_eq!(cached.len(), 3);
                                    assert_eq!(cached[0].to_bits(), 0x0);
                                    assert_eq!(cached[1].to_bits(), 0x0);
                                    assert_eq!(cached[2].to_bits(), O_20_23_BITS);
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 32);

        // Sequence 1: cross-key cold fill then warm in BOTH keys.
        let mut state = tpsa::DescriptorComputedState::default();
        assert_eq!(
            tpsa(&input, false, false, &mut state).unwrap().to_bits(),
            O_20_23_BITS
        );
        assert_eq!(
            tpsa(&input, true, false, &mut state).unwrap().to_bits(),
            O_20_23_BITS
        );
        assert_eq!(
            state.slot_for_tests(false).scalar.unwrap().to_bits(),
            O_20_23_BITS
        );
        assert_eq!(
            state.slot_for_tests(true).scalar.unwrap().to_bits(),
            O_20_23_BITS
        );

        // Sequence 2: warm clone serves BOTH APIs without recompute.
        let mut clone_state = state.clone();
        let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
        assert_eq!(
            tpsa(&input, false, false, &mut clone_state)
                .unwrap()
                .to_bits(),
            O_20_23_BITS
        );
        assert_eq!(
            tpsa::tpsa_contributions(&input, true, false, &mut clone_state).unwrap()[2].to_bits(),
            O_20_23_BITS
        );
        assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), helper_before);

        // Sequence 3: explicit clear empties BOTH slots -> recompute.
        state.clear();
        assert_eq!(state, tpsa::DescriptorComputedState::default());
        assert_eq!(
            tpsa(&input, false, false, &mut state).unwrap().to_bits(),
            O_20_23_BITS
        );

        // Sequence 4: forced overwrite replaces the sentinel scalar.
        state.slot_mut_for_tests(false).scalar = Some(999.0);
        assert_eq!(
            tpsa(&input, false, true, &mut state).unwrap().to_bits(),
            O_20_23_BITS
        );
        assert_eq!(
            state.slot_for_tests(false).scalar.unwrap().to_bits(),
            O_20_23_BITS
        );

        // Sequence 5: typed COLD invalid-valence leaves state unchanged.
        // The state is warm here; per the source order (hasProp guard at
        // MolSurf.cpp:351-355 PRECEDES computation) an unforced call
        // would serve the warm scalar, so force=true drives the cold
        // path through the contribution owner and its shape check.
        let good = valence(&topology, "t07r").unwrap();
        let mut malformed = good.clone();
        malformed.explicit_valence.resize(2, 0);
        let bad_input =
            DescriptorInput::new(&topology, &coordinates, &properties, &malformed, &rings);
        let snapshot = state.clone();
        assert!(matches!(
            tpsa(&bad_input, false, true, &mut state).unwrap_err(),
            DescriptorError::InvalidValenceRows { .. }
        ));
        assert_eq!(state, snapshot);
    }

    #[test]
    fn descriptor_t07_cache_product_malformed_warm_force_guard() {
        // Frozen T07-MATRIX guard product: 2 API kinds x 2 sandp x
        // force{false,true} = 8 REAL calls on ethanol whose supplied
        // FINAL assignment has explicit_valence shortened to 2 (required
        // 3) with the full implicit vector unchanged. Selected slot is
        // seeded scalar 777.0 + rows [7,8,9]; OTHER slot 789.0 +
        // [4,5,6]. !force serves the STORED warm values despite the
        // malformed unused assignment (the source hasProp guard at
        // MolSurf.cpp:351-355 precedes any computation); force=true
        // returns exactly InvalidValenceRows{function:
        // "tpsa_contributions", field: "explicit_valence", actual: 2,
        // expected: 3} (prepared_valence's shape check). Every call
        // leaves the ENTIRE cache and all supplied inputs unchanged and
        // has zero valence-helper reentry; count == 8.
        let topology = smiles_topology("CCO");
        let good = valence(&topology, "t07g").unwrap();
        let mut malformed = good.clone();
        malformed.explicit_valence.resize(2, 0);
        let rings = ring_info(&topology, "t07g").unwrap();
        let properties = MoleculeProperties::default();
        let coordinates = CoordinateBlock::default();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &malformed, &rings);
        let mut calls = 0u32;
        for api in 0..2 {
            for sandp in [false, true] {
                for force in [false, true] {
                    let mut state = tpsa::DescriptorComputedState::default();
                    let selected = state.slot_mut_for_tests(sandp);
                    selected.scalar = Some(777.0);
                    selected.contributions = Some(vec![7.0, 8.0, 9.0]);
                    let other = state.slot_mut_for_tests(!sandp);
                    other.scalar = Some(789.0);
                    other.contributions = Some(vec![4.0, 5.0, 6.0]);
                    let topology_snapshot = topology.clone();
                    let assignment_snapshot = malformed.clone();
                    let properties_snapshot = properties.clone();
                    let rings_snapshot = rings.clone();
                    let state_snapshot = state.clone();
                    let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                    match api {
                        0 => {
                            let result = tpsa::tpsa_contributions(&input, sandp, force, &mut state);
                            if !force {
                                // Warm rows served; malformed assignment
                                // never read.
                                assert_eq!(result.unwrap(), vec![7.0, 8.0, 9.0]);
                            } else {
                                assert!(matches!(
                                    result.unwrap_err(),
                                    DescriptorError::InvalidValenceRows {
                                        function: "tpsa_contributions",
                                        field: "explicit_valence",
                                        actual: 2,
                                        expected: 3,
                                    }
                                ));
                            }
                        }
                        _ => {
                            let result = tpsa(&input, sandp, force, &mut state);
                            if !force {
                                assert_eq!(result.unwrap().to_bits(), 777.0f64.to_bits());
                            } else {
                                assert!(matches!(
                                    result.unwrap_err(),
                                    DescriptorError::InvalidValenceRows {
                                        function: "tpsa_contributions",
                                        field: "explicit_valence",
                                        actual: 2,
                                        expected: 3,
                                    }
                                ));
                            }
                        }
                    }
                    calls += 1;
                    assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), helper_before);
                    assert_eq!(topology, topology_snapshot);
                    assert_eq!(malformed, assignment_snapshot);
                    assert_eq!(properties, properties_snapshot);
                    assert_eq!(rings, rings_snapshot);
                    assert_eq!(state, state_snapshot, "entire cache unchanged");
                }
            }
        }
        assert_eq!(calls, 8);
    }

    // ---- L01 Labute ASA (MolSurf.cpp:25-86) ----

    #[test]
    fn descriptor_l01_enum_scale_mapping_and_association() {
        // Frozen L01 semantics from the Step-340 audit. Two-atom
        // hand-built topologies with asserted radius inputs; every
        // expected value is hand-derived from the SOURCE arithmetic
        // (bij = Ri+Rj - facts[enum]; dij = min(max(|Ri-Rj|, bij),
        // Ri+Rj); Vi[begin] += Rj^2 - (Ri-dij)^2/dij; Vi[end] += Ri^2 -
        // (Rj-dij)^2/dij; final Vi[i] = PI*Ri*(4*Ri - Vi[i])), with the
        // radii coming from the existing rdkit_rb0 owner and the
        // per-arm bij/dij intermediates asserted from source formulas,
        // never from the implementation under test.
        let topology = smiles_topology("CC"); // ethane: C-C single
        let assignment = valence(&topology, "l01").unwrap();
        let rings = ring_info(&topology, "l01").unwrap();
        let properties = MoleculeProperties::default();
        let coordinates = CoordinateBlock::default();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let rc = cosmolkit_core::rdkit_rb0(6);
        let rh = cosmolkit_core::rdkit_rb0(1);

        // ONE bond, SINGLE (facts[1] = 0.0): bij = rc+rc - 0.0;
        // |Ri-Rj| = 0 so dij = max(0, bij) = bij (bij <= Ri+Rj holds).
        let bij_single = rc + rc - 0.0;
        let dij_single = (rc - rc).abs().max(bij_single).min(rc + rc);
        let contrib_begin = rc * rc - (rc - dij_single) * (rc - dij_single) / dij_single;
        let mut vi0 = contrib_begin;
        let vi0_final = std::f64::consts::PI * rc * (4.0 * rc - vi0);
        vi0 = vi0_final;
        let mut state = DescriptorComputedState::default();
        let result =
            labute_contributions(&input, false, false, &mut state).expect("single-bond labute");
        assert_eq!(result.atoms.len(), 2);
        // Both rows get the SAME accumulation by symmetry; compare with
        // tolerance because the final formula multiplies PI.
        assert!((result.atoms[0] - vi0).abs() < 1e-9);
        assert!((result.atoms[1] - vi0).abs() < 1e-9);
        assert_eq!(result.hydrogens, 0.0);
        // ASA total equals the row sum (no H term).
        let expected_asa = vi0 + vi0;
        let asa = labute_asa(&input, false, false, &mut state).unwrap();
        assert!((asa - expected_asa).abs() < 1e-9);

        // include_hydrogens pass on the same fixture: every atom gets
        // the implicit-H Vi increment and hContrib accumulates the
        // mirror term, transformed only above the 1e-4 gate.
        let bij_h = rc + rh;
        let dij_h = (rc - rh).abs().max(bij_h).min(rc + rh);
        let vi_inc = rh * rh - (rc - dij_h) * (rc - dij_h) / dij_h;
        let hc_raw = rc * rc - (rh - dij_h) * (rh - dij_h) / dij_h;
        let hc_transformed = std::f64::consts::PI * rh * (4.0 * rh - 2.0 * hc_raw);
        let mut state = DescriptorComputedState::default();
        let result = labute_contributions(&input, true, false, &mut state).expect("labute with H");
        let vi_h_final = std::f64::consts::PI * rc * (4.0 * rc - contrib_begin - vi_inc);
        assert!((result.atoms[0] - vi_h_final).abs() < 1e-9);
        assert!((result.atoms[1] - vi_h_final).abs() < 1e-9);
        // Per-atom mirror terms sum before the gate.
        assert!((result.hydrogens - hc_transformed).abs() < 1e-9);

        // ENUM-SCALE mapping product on hand-built [C, O] pairs: the
        // literal scale values {Single: 0.0, Double: 0.2, Triple: 0.3,
        // Aromatic: 0.1, Dative/other: none} follow facts[enumValue]
        // with bondType < 4 semantics; assert bij/dij from the source
        // formula per arm via the arm-specific scale.
        let ro = cosmolkit_core::rdkit_rb0(8);
        for (order, aromatic, scale) in [
            (BondOrder::Single, false, 0.0),
            (BondOrder::Double, false, 0.2),
            (BondOrder::Triple, false, 0.3),
            (BondOrder::Aromatic, true, 0.1),
            (
                BondOrder::Dative,
                false,
                0.0, /* none: >= 4 no scaling */
            ),
        ] {
            let mut atoms = vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
            ];
            atoms.truncate(2);
            let mut bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), order);
            if aromatic {
                bond_spec = bond_spec.with_aromatic(true);
            }
            let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
            let pair = TopologyBlock {
                adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
                atoms,
                bonds,
                ..TopologyBlock::default()
            };
            let pair_assignment = valence(&pair, "l01p").unwrap();
            let pair_rings = ring_info(&pair, "l01p").unwrap();
            let pair_input = DescriptorInput::new(
                &pair,
                &coordinates,
                &properties,
                &pair_assignment,
                &pair_rings,
            );
            // Source bij per arm: rc+ro-scale; dij clamped.
            let bij = rc + ro - scale;
            let dij = (rc - ro).abs().max(bij).min(rc + ro);
            let c_begin = ro * ro - (rc - dij) * (rc - dij) / dij;
            let c_end = rc * rc - (ro - dij) * (ro - dij) / dij;
            let v0 = std::f64::consts::PI * rc * (4.0 * rc - c_begin);
            let v1 = std::f64::consts::PI * ro * (4.0 * ro - c_end);
            let mut state = DescriptorComputedState::default();
            let result = labute_contributions(&pair_input, false, false, &mut state)
                .expect("enum-arm labute");
            assert!((result.atoms[0] - v0).abs() < 1e-9, "{order:?}");
            assert!((result.atoms[1] - v1).abs() < 1e-9, "{order:?}");
        }

        // ZERO-bond cell of the mapping product: bare C atom, Vi[i]=0
        // then PI*Ri*4Ri association; hContrib stays 0 and (below the
        // 1e-4 gate) is neither transformed nor added.
        let lone = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms: vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let lone_assignment = valence(&lone, "l01z").unwrap();
        let lone_rings = ring_info(&lone, "l01z").unwrap();
        let lone_input = DescriptorInput::new(
            &lone,
            &coordinates,
            &properties,
            &lone_assignment,
            &lone_rings,
        );
        let mut state = DescriptorComputedState::default();
        let result =
            labute_contributions(&lone_input, true, false, &mut state).expect("zero-bond labute");
        let v_lone = std::f64::consts::PI * rc * 4.0 * rc;
        // With include_hydrogens the implicit-H pass still runs on the
        // lone atom; use the include pass formula.
        let bij_l = rc + rh;
        let dij_l = (rc - rh).abs().max(bij_l).min(rc + rh);
        let v_lone_h = std::f64::consts::PI
            * rc
            * (4.0 * rc - (rh * rh - (rc - dij_l) * (rc - dij_l) / dij_l));
        assert!((result.atoms[0] - v_lone_h).abs() < 1e-9);
        let hc_l = rc * rc - (rh - dij_l) * (rh - dij_l) / dij_l;
        let hc_l_t = std::f64::consts::PI * rh * (4.0 * rh - hc_l);
        assert!((result.hydrogens - hc_l_t).abs() < 1e-9);
        let _ = (v_lone, hc_l);
        // State published all three values (rows, H, ASA).
        assert!(state.labute_slot().asa.is_some());
        assert!(state.labute_slot().rows.is_some());
        assert!(state.labute_slot().hydrogens.is_some());
    }

    #[test]
    fn descriptor_l01_complete_enum_aromatic_h_product() {
        // Frozen L01-ENUM product: ALL 22 named enum variants in numeric
        // source order x aromatic{false,true} x include_hydrogens{false,
        // true} = 88 actual contribution calls, plus a two-carbon
        // no-bond contrast x include_hydrogens = 2 calls; each of the 90
        // cases uses a FRESH empty state, calls labute_contributions
        // once, then labute_asa once with the SAME flags/state (!force):
        // 180 actual entry calls. Radii: C 0.77, H 0.33 (pinned
        // atomic_data.cpp). Expectations are the frozen literal bits from
        // the supervisor's standalone source-arithmetic reference, never
        // a formula/helper copying the implementation.
        const NO_BOND: (u64, u64, u64) = (0x401d_cd6a_6270_2d7d, 0x0, 0x402d_cd6a_6270_2d7d);
        const NO_BOND_H: (u64, u64, u64) = (
            0x401d_b4e4_7726_7dcc,
            0x3ff4_1b85_1ccc_a042,
            0x4030_1c2a_8d60_08ea,
        );
        const ZERO_SCALE: (u64, u64, u64) = (0x401b_ca6e_1564_c404, 0x0, 0x402b_ca6e_1564_c404);
        const ZERO_SCALE_H: (u64, u64, u64) = (
            0x401b_b1e8_2a1b_1454,
            0x3ff4_1b85_1ccc_a042,
            0x402e_3558_cdb4_a85c,
        );
        const SCALE_01: (u64, u64, u64) = (0x401b_14f3_0426_147b, 0x0, 0x402b_14f3_0426_147b);
        const SCALE_01_H: (u64, u64, u64) = (
            0x401a_fc6d_18dc_64cb,
            0x3ff4_1b85_1ccc_a042,
            0x402d_7fdd_bc75_f8d3,
        );
        const SCALE_02: (u64, u64, u64) = (0x401a_695a_6f57_dbab, 0x0, 0x402a_695a_6f57_dbab);
        const SCALE_02_H: (u64, u64, u64) = (
            0x401a_50d4_840e_2bfb,
            0x3ff4_1b85_1ccc_a042,
            0x402c_d445_27a7_c003,
        );
        const SCALE_03: (u64, u64, u64) = (0x4019_ca08_8ddb_809d, 0x0, 0x4029_ca08_8ddb_809d);
        const SCALE_03_H: (u64, u64, u64) = (
            0x4019_b182_a291_d0ec,
            0x3ff4_1b85_1ccc_a042,
            0x402c_34f3_462b_64f4,
        );

        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        use cosmolkit_model::BondOrder as BO;
        // Literal Bond.h enum codes and frozen non-aromatic outputs.
        // This table is independent of the implementation's branch logic;
        // aromatic=true has the fixed SCALE_01/SCALE_01_H source outputs.
        let variants = [
            (BO::Unspecified, 0, [SCALE_01, SCALE_01_H]),
            (BO::Single, 1, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::Double, 2, [SCALE_02, SCALE_02_H]),
            (BO::Triple, 3, [SCALE_03, SCALE_03_H]),
            (BO::Quadruple, 4, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::Quintuple, 5, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::Hextuple, 6, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::OneAndHalf, 7, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::TwoAndHalf, 8, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::ThreeAndHalf, 9, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::FourAndHalf, 10, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::FiveAndHalf, 11, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::Aromatic, 12, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::Ionic, 13, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::Hydrogen, 14, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::ThreeCenter, 15, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::DativeOne, 16, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::Dative, 17, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::DativeLeft, 18, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::DativeRight, 19, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::Other, 20, [ZERO_SCALE, ZERO_SCALE_H]),
            (BO::Zero, 21, [ZERO_SCALE, ZERO_SCALE_H]),
        ];
        let mut entry_calls = 0u32;
        let mut contribution_calls = 0u32;
        let rc = cosmolkit_core::rdkit_rb0(6);
        assert!((rc - 0.77).abs() < 1e-12, "pinned C radius");
        assert!(
            (cosmolkit_core::rdkit_rb0(1) - 0.33).abs() < 1e-12,
            "pinned H radius"
        );
        let mut run_case = |topology: &TopologyBlock,
                            order: Option<(BO, i64)>,
                            aromatic: bool,
                            include_h: bool,
                            expected: (u64, u64, u64)| {
            // Hand-built supplied FINAL assignment: the Labute owner
            // reads only topology radii/bonds (no valence H counts),
            // and the canonical valence owner legitimately rejects the
            // hand-built enum fixtures the product must exercise.
            let assignment = t01_zero_assignment(2);
            let rings = cosmolkit_core::RingInfo::new(
                cosmolkit_core::RingFindType::OtherOrUnknown,
                2,
                topology.bonds.len(),
            );
            let input =
                DescriptorInput::new(topology, &coordinates, &properties, &assignment, &rings);
            // Prerequisites: two carbons, enum code, independent flag,
            // exact bond presence/endpoints when bonded.
            assert_eq!(topology.atoms.len(), 2);
            assert_eq!(topology.atoms[0].element(), Element::C);
            assert_eq!(topology.atoms[1].element(), Element::C);
            if let Some((order, code)) = order {
                assert_eq!(topology.bonds.len(), 1);
                assert_eq!(topology.bonds[0].order(), order);
                assert_eq!(topology.bonds[0].order().rdkit_code(), code);
                assert_eq!(topology.bonds[0].is_aromatic(), aromatic);
                assert_eq!(topology.bonds[0].begin().index(), 0);
                assert_eq!(topology.bonds[0].end().index(), 1);
            } else {
                assert!(topology.bonds.is_empty());
            }
            let topology_snapshot = topology.clone();
            let coordinates_snapshot = coordinates.clone();
            let assignment_snapshot = assignment.clone();
            let properties_snapshot = properties.clone();
            let rings_snapshot = rings.clone();
            let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
            let mut state = DescriptorComputedState::default();
            let result = labute_contributions(&input, include_h, false, &mut state)
                .expect("enum-product contributions");
            contribution_calls += 1;
            entry_calls += 1;
            assert_eq!(result.atoms.len(), 2);
            assert_eq!(result.atoms[0].to_bits(), expected.0, "row0");
            assert_eq!(result.atoms[1].to_bits(), expected.0, "row1 equal");
            assert_eq!(result.hydrogens.to_bits(), expected.1, "H");
            let asa = labute_asa(&input, include_h, false, &mut state).expect("enum-product asa");
            entry_calls += 1;
            assert_eq!(asa.to_bits(), expected.2, "asa total");
            assert_eq!(
                VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                helper_before,
                "zero valence-helper reentry"
            );
            assert_eq!(*topology, topology_snapshot);
            assert_eq!(coordinates, coordinates_snapshot);
            assert_eq!(assignment, assignment_snapshot);
            assert_eq!(properties, properties_snapshot);
            assert_eq!(rings, rings_snapshot);
        };
        for (order, code, non_aromatic_outputs) in variants {
            for aromatic in [false, true] {
                // Aromatic flag on the Aromatic enum is the parse-stage
                // default pairing; the independent flag product exercises
                // the flag axis regardless of enum.
                for include_h in [false, true] {
                    let mut bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), order);
                    if aromatic {
                        bond_spec = bond_spec.with_aromatic(true);
                    }
                    let atoms = vec![
                        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                    ];
                    let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
                    let topology = TopologyBlock {
                        adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
                        atoms,
                        bonds,
                        ..TopologyBlock::default()
                    };
                    let outputs = if aromatic {
                        [SCALE_01, SCALE_01_H]
                    } else {
                        non_aromatic_outputs
                    };
                    run_case(
                        &topology,
                        Some((order, code)),
                        aromatic,
                        include_h,
                        outputs[usize::from(include_h)],
                    );
                }
            }
        }
        // No-bond contrast: two carbons, no bond; outputs MUST differ
        // from the actual ZERO bond cell.
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let no_bond = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        for include_h in [false, true] {
            let expected = if include_h { NO_BOND_H } else { NO_BOND };
            run_case(&no_bond, None, false, include_h, expected);
            // Discriminator: the actual ZERO bond has DIFFERENT bits.
            let zero_expected = [ZERO_SCALE, ZERO_SCALE_H][usize::from(include_h)];
            assert_ne!(expected.0, zero_expected.0, "no-bond != ZERO bond");
            assert_ne!(expected.2, zero_expected.2);
        }
        assert_eq!(contribution_calls, 90);
        assert_eq!(entry_calls, 180);
    }

    #[test]
    fn descriptor_l02_include_hs_false_true_product() {
        // Frozen L02 expectations (Step-348 audit): literal bits from
        // the standalone source-arithmetic reference for MolSurf.cpp
        // 58-81 over the pinned rb0 radii (C 0.77, H 0.33), never from
        // the implementation under test.
        const ETHANE_NO_H_ROW: u64 = 0x401b_ca6e_1564_c404;
        const ETHANE_NO_H_H: u64 = 0x0;
        const ETHANE_NO_H_TOTAL: u64 = 0x402b_ca6e_1564_c404;
        const ETHANE_H_ROW: u64 = 0x401b_b1e8_2a1b_1454;
        const ETHANE_H_RAW: u64 = 0x3fbb_98c7_e282_40c0;
        const ETHANE_H_RETURNED: u64 = 0x3ff4_1b85_1ccc_a042;
        const ETHANE_H_TOTAL: u64 = 0x402e_3558_cdb4_a85c;
        assert_eq!(cosmolkit_core::rdkit_rb0(6).to_bits(), 0.77_f64.to_bits());
        assert_eq!(cosmolkit_core::rdkit_rb0(1).to_bits(), 0.33_f64.to_bits());
        // Gate discriminator: the raw per-atom sum differs from the
        // PI*Rj*(4Rj-hContrib) transformed value actually returned.
        assert_ne!(ETHANE_H_RAW, ETHANE_H_RETURNED);

        let bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        for (include_h, row, hydrogen, total) in [
            (false, ETHANE_NO_H_ROW, ETHANE_NO_H_H, ETHANE_NO_H_TOTAL),
            (true, ETHANE_H_ROW, ETHANE_H_RETURNED, ETHANE_H_TOTAL),
        ] {
            let assignment = t01_zero_assignment(2);
            let rings =
                cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 1);
            let input =
                DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
            let topology_snapshot = topology.clone();
            let coordinates_snapshot = coordinates.clone();
            let assignment_snapshot = assignment.clone();
            let properties_snapshot = properties.clone();
            let rings_snapshot = rings.clone();
            let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
            let mut state = DescriptorComputedState::default();
            let result =
                labute_contributions(&input, include_h, false, &mut state).expect("ethane labute");
            assert_eq!(result.atoms.len(), 2, "one row per explicit atom");
            assert_eq!(result.atoms[0].to_bits(), row);
            assert_eq!(result.atoms[1].to_bits(), row);
            assert_eq!(result.hydrogens.to_bits(), hydrogen);
            let asa = labute_asa(&input, include_h, false, &mut state).unwrap();
            assert_eq!(asa.to_bits(), total);
            assert_eq!(topology, topology_snapshot);
            assert_eq!(coordinates, coordinates_snapshot);
            assert_eq!(assignment, assignment_snapshot);
            assert_eq!(properties, properties_snapshot);
            assert_eq!(rings, rings_snapshot);
            assert_eq!(
                VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                helper_before,
                "no valence-helper reentry"
            );
        }
    }

    #[test]
    fn descriptor_l02_once_per_explicit_atom_not_h_multiplication() {
        // The source includeHs loop runs ONCE PER EXPLICIT ATOM over
        // atom indices; it never multiplies by attached (implicit or
        // explicit) H counts. Frozen bits from the standalone reference.
        const ETHANE_H_ROW: u64 = 0x401b_b1e8_2a1b_1454;
        const ETHANE_H_RETURNED: u64 = 0x3ff4_1b85_1ccc_a042;
        const ETHANE_H_TOTAL: u64 = 0x402e_3558_cdb4_a85c;
        const NODE_NO_H_C: u64 = 0x401d_cd6a_6270_2d7d;
        const NODE_NO_H_HROW: u64 = 0x3ff5_e548_ef81_6fca;
        const NODE_NO_H_TOTAL: u64 = 0x4021_a35e_4f28_44b8;
        const NODE_H_C: u64 = 0x401d_b4e4_7726_7dcc;
        const NODE_H_HROW: u64 = 0x3ff6_d382_6f71_d157;
        const NODE_H_RETURNED: u64 = 0x3ff5_eea0_8617_6993;
        const NODE_H_TOTAL: u64 = 0x4024_72b6_9a44_6643;

        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        // (a) plain C-C and (b) C(H3)-C(H3): identical per-atom rows,
        // identical H term, identical total — the loop counted TWO
        // atoms in both, never 2*(1+3) H sites.
        for explicit_h in [0, 3] {
            let bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
            let atoms = vec![
                Atom::from_spec(
                    AtomId::new(0),
                    AtomSpec::new(Element::C).with_explicit_hydrogens(explicit_h),
                ),
                Atom::from_spec(
                    AtomId::new(1),
                    AtomSpec::new(Element::C).with_explicit_hydrogens(explicit_h),
                ),
            ];
            let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
            let topology = TopologyBlock {
                adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
                atoms,
                bonds,
                ..TopologyBlock::default()
            };
            let assignment = t01_zero_assignment(2);
            let rings =
                cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 1);
            let input =
                DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
            let mut state = DescriptorComputedState::default();
            let result = labute_contributions(&input, true, false, &mut state)
                .expect("explicit-H variant labute");
            assert_eq!(result.atoms.len(), 2, "explicit_h = {explicit_h}");
            assert_eq!(result.atoms[0].to_bits(), ETHANE_H_ROW);
            assert_eq!(result.atoms[1].to_bits(), ETHANE_H_ROW);
            assert_eq!(result.hydrogens.to_bits(), ETHANE_H_RETURNED);
            let asa = labute_asa(&input, true, false, &mut state).unwrap();
            assert_eq!(asa.to_bits(), ETHANE_H_TOTAL);
        }

        // (c) An explicit H NODE counts as exactly ONE atom with its
        // own row; the pass adds one iteration for it, nothing more.
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::H)),
        ];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        for (include_h, c_row, h_row, hydrogen, total) in [
            (false, NODE_NO_H_C, NODE_NO_H_HROW, 0x0_u64, NODE_NO_H_TOTAL),
            (true, NODE_H_C, NODE_H_HROW, NODE_H_RETURNED, NODE_H_TOTAL),
        ] {
            let assignment = t01_zero_assignment(2);
            let rings =
                cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 0);
            let input =
                DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
            let mut state = DescriptorComputedState::default();
            let result =
                labute_contributions(&input, include_h, false, &mut state).expect("C-H node");
            assert_eq!(result.atoms.len(), 2, "one row per node, H node included");
            assert_eq!(result.atoms[0].to_bits(), c_row);
            assert_eq!(result.atoms[1].to_bits(), h_row);
            assert_eq!(result.hydrogens.to_bits(), hydrogen);
            let asa = labute_asa(&input, include_h, false, &mut state).unwrap();
            assert_eq!(asa.to_bits(), total);
        }
    }

    #[test]
    fn descriptor_l02_fabs_threshold_gate_both_arms() {
        // Frozen fabs(hContrib) > 1e-4 gate semantics (MolSurf.cpp:78).
        // BELOW threshold: hContrib stays the UNTRANSFORMED raw sum and
        // is NOT added to res (still cached). ABOVE: transformed and
        // added. The below-threshold fixture is the smallest genuine
        // nonzero lattice over the pinned rb0 radii: 6 atoms of radius
        // 0.7 (He3 N3), 10 O, 3 F, 2 Ne -> raw sum -4.609e-05.
        const T_ROW_R07_H: u64 = 0x4018_9a28_d28c_97b0;
        const T_ROW_O_H: u64 = 0x4015_e79e_d526_ee3d;
        const T_ROW_F_H: u64 = 0x4012_d14d_48fb_df16;
        const T_ROW_R07: u64 = 0x4018_a14d_57b3_73df;
        const T_ROW_O: u64 = 0x4015_e548_ef81_6fca;
        const T_ROW_F: u64 = 0x4012_c3e1_898e_3ce3;
        const T_RAW_H: u64 = 0xbf08_2a0d_ae6b_e000;
        const T_TOTAL_WITHOUT_H: u64 = 0x405d_8516_2c2d_da93;
        const T_TOTAL_NO_H: u64 = 0x405d_84ae_8b55_4b39;
        // Gate-open discriminators (what a wrongly-open gate WOULD
        // produce): transformed H and total-with-H.
        const T_H_WOULD_BE: u64 = 0x3ff5_e57b_09fc_3908;
        const T_TOTAL_WOULD_BE: u64 = 0x405d_dcac_1855_cb77;
        const ETHANE_H_RAW: u64 = 0x3fbb_98c7_e282_40c0;
        const ETHANE_H_RETURNED: u64 = 0x3ff4_1b85_1ccc_a042;
        const ETHANE_H_TOTAL: u64 = 0x402e_3558_cdb4_a85c;
        const CM_ROW: u64 = 0x0;
        const CM_GATE_OPEN_H: u64 = 0x3ff5_e548_ef81_6fca;
        // Pinned radii prerequisites (atomic_data.cpp rows).
        assert_eq!(cosmolkit_core::rdkit_rb0(2).to_bits(), 0.7_f64.to_bits());
        assert_eq!(cosmolkit_core::rdkit_rb0(7).to_bits(), 0.7_f64.to_bits());
        assert_eq!(cosmolkit_core::rdkit_rb0(10).to_bits(), 0.7_f64.to_bits());
        assert_eq!(cosmolkit_core::rdkit_rb0(8).to_bits(), 0.66_f64.to_bits());
        assert_eq!(cosmolkit_core::rdkit_rb0(9).to_bits(), 0.611_f64.to_bits());
        assert_eq!(cosmolkit_core::rdkit_rb0(96).to_bits(), 0.0_f64.to_bits());
        // Prerequisite: the raw sum is genuinely inside the threshold.
        assert!(f64::from_bits(T_RAW_H).abs() <= 1e-4);
        assert!(f64::from_bits(ETHANE_H_RAW).abs() > 1e-4);

        let elements = [
            Element::HE,
            Element::HE,
            Element::HE,
            Element::N,
            Element::N,
            Element::N,
            Element::O,
            Element::O,
            Element::O,
            Element::O,
            Element::O,
            Element::O,
            Element::O,
            Element::O,
            Element::O,
            Element::O,
            Element::F,
            Element::F,
            Element::F,
            Element::NE,
            Element::NE,
        ];
        let atoms: Vec<Atom> = elements
            .iter()
            .enumerate()
            .map(|(i, el)| Atom::from_spec(AtomId::new(i), AtomSpec::new(*el)))
            .collect();
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let expected_row = |i: usize, include_h: bool| -> u64 {
            match elements[i] {
                Element::HE | Element::N | Element::NE => {
                    if include_h {
                        T_ROW_R07_H
                    } else {
                        T_ROW_R07
                    }
                }
                Element::O => {
                    if include_h {
                        T_ROW_O_H
                    } else {
                        T_ROW_O
                    }
                }
                _ => {
                    if include_h {
                        T_ROW_F_H
                    } else {
                        T_ROW_F
                    }
                }
            }
        };
        for (include_h, hydrogen, total) in [
            (true, T_RAW_H, T_TOTAL_WITHOUT_H),
            (false, 0x0, T_TOTAL_NO_H),
        ] {
            let assignment = t01_zero_assignment(elements.len());
            let rings = cosmolkit_core::RingInfo::new(
                cosmolkit_core::RingFindType::OtherOrUnknown,
                elements.len(),
                0,
            );
            let input =
                DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
            let mut state = DescriptorComputedState::default();
            let result = labute_contributions(&input, include_h, false, &mut state)
                .expect("threshold fixture");
            assert_eq!(result.atoms.len(), elements.len());
            for (i, row) in result.atoms.iter().enumerate() {
                assert_eq!(row.to_bits(), expected_row(i, include_h), "row {i}");
            }
            // Below-threshold semantics: hydrogens IS the raw sum,
            // untransformed, and the total EXCLUDES it.
            assert_eq!(result.hydrogens.to_bits(), hydrogen);
            assert_ne!(result.hydrogens.to_bits(), T_H_WOULD_BE);
            let asa = labute_asa(&input, include_h, false, &mut state).unwrap();
            assert_eq!(asa.to_bits(), total);
            assert_ne!(asa.to_bits(), T_TOTAL_WOULD_BE);
        }

        // ABOVE-threshold arm: ethane includeHs — raw 0.1078 is past
        // the gate, so hydrogens is the transformed value and the
        // total includes it.
        let bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
        let ethane = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(2);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 1);
        let input = DescriptorInput::new(&ethane, &coordinates, &properties, &assignment, &rings);
        let mut state = DescriptorComputedState::default();
        let result = labute_contributions(&input, true, false, &mut state).expect("above gate");
        assert_eq!(result.hydrogens.to_bits(), ETHANE_H_RETURNED);
        assert_ne!(result.hydrogens.to_bits(), ETHANE_H_RAW);
        let asa = labute_asa(&input, true, false, &mut state).unwrap();
        assert_eq!(asa.to_bits(), ETHANE_H_TOTAL);

        // Exact-zero secondary cell: the pinned table gives Curium
        // Rb0 = 0.0 (atomic_data.cpp Cm row), so raw hContrib is
        // exactly +0.0 — at-or-below threshold on BOTH comparisons:
        // untransformed, not added. A wrongly-open gate would return
        // PI*0.33*1.32 and add it.
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::CM)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::CM)),
        ];
        let curium = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        for include_h in [false, true] {
            let assignment = t01_zero_assignment(2);
            let rings =
                cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 0);
            let input =
                DescriptorInput::new(&curium, &coordinates, &properties, &assignment, &rings);
            let mut state = DescriptorComputedState::default();
            let result =
                labute_contributions(&input, include_h, false, &mut state).expect("Cm cell");
            assert_eq!(result.atoms.len(), 2);
            assert_eq!(result.atoms[0].to_bits(), CM_ROW);
            assert_eq!(result.atoms[1].to_bits(), CM_ROW);
            assert_eq!(result.hydrogens.to_bits(), 0x0);
            assert_ne!(result.hydrogens.to_bits(), CM_GATE_OPEN_H);
            let asa = labute_asa(&input, include_h, false, &mut state).unwrap();
            assert_eq!(asa.to_bits(), 0x0);
            assert_ne!(asa.to_bits(), CM_GATE_OPEN_H);
        }
    }

    #[test]
    fn descriptor_l02_empty_molecule() {
        // nAtoms = 0: empty rows, exactly +0.0 H term and total, and
        // the slot is still published after success.
        let topology = TopologyBlock::default();
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        for include_h in [false, true] {
            let assignment = t01_zero_assignment(0);
            let rings =
                cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 0, 0);
            let input =
                DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
            let mut state = DescriptorComputedState::default();
            let result = labute_contributions(&input, include_h, false, &mut state).expect("empty");
            assert!(result.atoms.is_empty());
            assert_eq!(result.hydrogens.to_bits(), 0x0);
            let asa = labute_asa(&input, include_h, false, &mut state).unwrap();
            assert_eq!(asa.to_bits(), 0x0);
            let slot = state.labute_slot();
            let rows = slot.rows.clone().expect("rows published");
            let hydrogen = slot.hydrogens.expect("H published");
            assert!(rows.is_empty());
            assert_eq!(hydrogen.to_bits(), 0x0);
            assert_eq!(slot.asa, Some(0.0));
        }
    }

    #[test]
    fn descriptor_l03_call_sequences_first_result_reuse() {
        // Frozen L03 semantics (Step-356 audit): the source cache key is
        // hasProp(_labuteAtomContribs) — NOT includeHs-keyed, NOT
        // force-keyed. First-result reuse: every later non-force call is
        // served the FIRST computed values regardless of its includeHs
        // argument; force bypasses and OVERWRITES; Clone carries the cache.
        // Frozen literals from the L02 source-arithmetic tables.
        const FALSE_ROW: u64 = 0x401b_ca6e_1564_c404;
        const FALSE_H: u64 = 0x0;
        const FALSE_TOTAL: u64 = 0x402b_ca6e_1564_c404;
        const TRUE_ROW: u64 = 0x401b_b1e8_2a1b_1454;
        const TRUE_H: u64 = 0x3ff4_1b85_1ccc_a042;
        const TRUE_TOTAL: u64 = 0x402e_3558_cdb4_a85c;

        let bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let assignment = t01_zero_assignment(2);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 1);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let topology_snapshot = topology.clone();
        let coordinates_snapshot = coordinates.clone();
        let assignment_snapshot = assignment.clone();
        let properties_snapshot = properties.clone();
        let rings_snapshot = rings.clone();

        let check = |result: &crate::labute::LabuteContributions,
                     asa: f64,
                     row: u64,
                     hydrogen: u64,
                     total: u64| {
            assert_eq!(result.atoms.len(), 2);
            assert_eq!(result.atoms[0].to_bits(), row);
            assert_eq!(result.atoms[1].to_bits(), row);
            assert_eq!(result.hydrogens.to_bits(), hydrogen);
            assert_eq!(asa.to_bits(), total);
        };

        // Sequence 1: cold(false) then warm(true) — the warm call is
        // served the FALSE-computed values (first-result reuse).
        let mut state = DescriptorComputedState::default();
        let r = labute_contributions(&input, false, false, &mut state).unwrap();
        check(
            &r,
            labute_asa(&input, false, false, &mut state).unwrap(),
            FALSE_ROW,
            FALSE_H,
            FALSE_TOTAL,
        );
        let r = labute_contributions(&input, true, false, &mut state).unwrap();
        check(
            &r,
            labute_asa(&input, true, false, &mut state).unwrap(),
            FALSE_ROW,
            FALSE_H,
            FALSE_TOTAL,
        );

        // Sequence 2: cold(true) then warm(false) — served the
        // TRUE-computed values.
        let mut state = DescriptorComputedState::default();
        let r = labute_contributions(&input, true, false, &mut state).unwrap();
        check(
            &r,
            labute_asa(&input, true, false, &mut state).unwrap(),
            TRUE_ROW,
            TRUE_H,
            TRUE_TOTAL,
        );
        let r = labute_contributions(&input, false, false, &mut state).unwrap();
        check(
            &r,
            labute_asa(&input, false, false, &mut state).unwrap(),
            TRUE_ROW,
            TRUE_H,
            TRUE_TOTAL,
        );

        // Sequence 3: Clone carries the cache — a cloned state serves the
        // carried values, and the original keeps its own.
        let mut state = DescriptorComputedState::default();
        let r = labute_contributions(&input, true, false, &mut state).unwrap();
        check(
            &r,
            labute_asa(&input, true, false, &mut state).unwrap(),
            TRUE_ROW,
            TRUE_H,
            TRUE_TOTAL,
        );
        let mut cloned = state.clone();
        let r = labute_contributions(&input, false, false, &mut cloned).unwrap();
        check(
            &r,
            labute_asa(&input, false, false, &mut cloned).unwrap(),
            TRUE_ROW,
            TRUE_H,
            TRUE_TOTAL,
        );
        assert_eq!(cloned, state);

        // Sequence 4: force bypasses the guard and OVERWRITES the stored
        // values; later non-force calls serve the NEW ones.
        let mut state = DescriptorComputedState::default();
        let r = labute_contributions(&input, false, false, &mut state).unwrap();
        check(
            &r,
            labute_asa(&input, false, false, &mut state).unwrap(),
            FALSE_ROW,
            FALSE_H,
            FALSE_TOTAL,
        );
        let r = labute_contributions(&input, true, true, &mut state).unwrap();
        check(
            &r,
            labute_asa(&input, true, false, &mut state).unwrap(),
            TRUE_ROW,
            TRUE_H,
            TRUE_TOTAL,
        );
        let r = labute_contributions(&input, false, false, &mut state).unwrap();
        check(
            &r,
            labute_asa(&input, false, false, &mut state).unwrap(),
            TRUE_ROW,
            TRUE_H,
            TRUE_TOTAL,
        );

        // Inputs untouched and no valence-helper reentry anywhere above.
        assert_eq!(topology, topology_snapshot);
        assert_eq!(coordinates, coordinates_snapshot);
        assert_eq!(assignment, assignment_snapshot);
        assert_eq!(properties, properties_snapshot);
        assert_eq!(rings, rings_snapshot);
        assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), 0);
    }

    #[test]
    fn descriptor_l03_frozen_presence_routes_64() {
        // Exact frozen MolSurf.cpp:27-34 / 89-101 guard table, not an
        // expected-value calculation copied from the implementation.
        #[derive(Clone, Copy, Debug)]
        enum Action {
            Compute,
            MissingH,
            MissingAsa,
            Serve,
        }
        use Action::*;
        const CONTRIBUTION: [Action; 8] = [
            Compute, MissingH, Compute, MissingAsa, Compute, MissingH, Compute, Serve,
        ];
        const SCALAR: [Action; 8] = [
            Compute, MissingH, Compute, MissingAsa, Serve, Serve, Serve, Serve,
        ];
        // Unspecified=0 -> source scale0.1, independently frozen before L03.
        const COMPUTED: [(u64, u64, u64); 2] = [
            (0x401b14f30426147b, 0, 0x402b14f30426147b),
            (0x401afc6d18dc64cb, 0x3ff41b851ccca042, 0x402d7fddbc75f8d3),
        ];
        #[derive(Debug)]
        enum Observation {
            Contributions(LabuteContributions),
            Scalar(f64),
        }

        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified),
        )];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(2, &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let assignment = t01_zero_assignment(2);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 1);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut calls = [0u32; 2];
        for scalar in [false, true] {
            for mask in 0..8 {
                for force in [false, true] {
                    for include_h in [false, true] {
                        assert_eq!(topology.atoms.len(), 2);
                        assert!(topology.atoms.iter().all(|a| a.element() == Element::C));
                        assert_eq!(topology.bonds[0].order().rdkit_code(), 0);
                        assert!(!topology.bonds[0].is_aromatic());
                        assert_eq!(topology.bonds[0].begin().index(), 0);
                        assert_eq!(topology.bonds[0].end().index(), 1);
                        let mut state = DescriptorComputedState::default();
                        let slot = state.labute_slot_mut();
                        slot.rows = (mask & 1 != 0).then(|| vec![4.0, 5.0]);
                        slot.hydrogens = (mask & 2 != 0).then_some(7.0);
                        slot.asa = (mask & 4 != 0).then_some(789.0);
                        for (sandp, value, rows) in [
                            (false, 123.0, vec![1.0, 2.0]),
                            (true, 456.0, vec![3.0, 4.0]),
                        ] {
                            let slot = state.slot_mut_for_tests(sandp);
                            slot.scalar = Some(value);
                            slot.contributions = Some(rows);
                        }
                        let before_state = state.clone();
                        let before_inputs = (
                            topology.clone(),
                            coordinates.clone(),
                            properties.clone(),
                            assignment.clone(),
                            rings.clone(),
                        );
                        let before_counter = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                        calls[usize::from(scalar)] += 1;
                        let outcome = if scalar {
                            labute_asa(&input, include_h, force, &mut state)
                                .map(Observation::Scalar)
                        } else {
                            labute_contributions(&input, include_h, force, &mut state)
                                .map(Observation::Contributions)
                        };
                        let action = if force {
                            Compute
                        } else if scalar {
                            SCALAR[mask]
                        } else {
                            CONTRIBUTION[mask]
                        };
                        let mut expected_state = before_state.clone();
                        match action {
                            MissingH => assert_eq!(
                                outcome.unwrap_err(),
                                DescriptorError::MissingLabuteHydrogens {
                                    function: "labute_contributions"
                                }
                            ),
                            MissingAsa => assert_eq!(
                                outcome.unwrap_err(),
                                DescriptorError::MissingLabuteAsa {
                                    function: "labute_contributions"
                                }
                            ),
                            Compute | Serve => {
                                let (rows, hydrogen, asa) = match action {
                                    Compute => {
                                        let (r, h, total) = COMPUTED[usize::from(include_h)];
                                        (vec![r, r], h, total)
                                    }
                                    Serve => (
                                        vec![4.0f64.to_bits(), 5.0f64.to_bits()],
                                        7.0f64.to_bits(),
                                        789.0f64.to_bits(),
                                    ),
                                    _ => unreachable!(),
                                };
                                match outcome.expect("literal source success action") {
                                    Observation::Contributions(result) => {
                                        assert_eq!(result.atoms.len(), rows.len());
                                        assert_eq!(
                                            result
                                                .atoms
                                                .iter()
                                                .map(|v| v.to_bits())
                                                .collect::<Vec<_>>(),
                                            rows
                                        );
                                        assert_eq!(result.hydrogens.to_bits(), hydrogen);
                                    }
                                    Observation::Scalar(result) => {
                                        assert_eq!(result.to_bits(), asa)
                                    }
                                }
                                if matches!(action, Compute) {
                                    let slot = expected_state.labute_slot_mut();
                                    slot.rows =
                                        Some(rows.iter().copied().map(f64::from_bits).collect());
                                    slot.hydrogens = Some(f64::from_bits(hydrogen));
                                    slot.asa = Some(f64::from_bits(asa));
                                }
                            }
                        }
                        assert_eq!(
                            state, expected_state,
                            "complete cache, route={scalar} mask={mask}"
                        );
                        let slot = state.labute_slot();
                        let expected_slot = expected_state.labute_slot();
                        assert_eq!(
                            slot.rows
                                .as_ref()
                                .map(|rows| rows.iter().map(|v| v.to_bits()).collect::<Vec<_>>()),
                            expected_slot
                                .rows
                                .as_ref()
                                .map(|rows| rows.iter().map(|v| v.to_bits()).collect::<Vec<_>>())
                        );
                        assert_eq!(
                            slot.hydrogens.map(f64::to_bits),
                            expected_slot.hydrogens.map(f64::to_bits)
                        );
                        assert_eq!(
                            slot.asa.map(f64::to_bits),
                            expected_slot.asa.map(f64::to_bits)
                        );
                        for (sandp, value, rows) in [
                            (false, 123.0f64, [1.0f64, 2.0]),
                            (true, 456.0f64, [3.0f64, 4.0]),
                        ] {
                            let slot = state.slot_for_tests(sandp);
                            assert_eq!(slot.scalar.map(f64::to_bits), Some(value.to_bits()));
                            assert_eq!(
                                slot.contributions
                                    .as_ref()
                                    .map(|v| v.iter().map(|v| v.to_bits()).collect::<Vec<_>>()),
                                Some(rows.map(f64::to_bits).to_vec())
                            );
                        }
                        assert_eq!(
                            (
                                topology.clone(),
                                coordinates.clone(),
                                properties.clone(),
                                assignment.clone(),
                                rings.clone()
                            ),
                            before_inputs
                        );
                        assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), before_counter);
                    }
                }
            }
        }
        assert_eq!(calls, [32, 32]);
    }

    #[test]
    fn descriptor_l03_frozen_h_sequences_24() {
        // Four first-H/second-H pairs; each performs exactly six real calls:
        // cold, warm, cloned warm, forced clone, warm clone, cleared-clone cold.
        const COMPUTED: [(u64, u64, u64); 2] = [
            (0x401b14f30426147b, 0, 0x402b14f30426147b),
            (0x401afc6d18dc64cb, 0x3ff41b851ccca042, 0x402d7fddbc75f8d3),
        ];
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified),
        )];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(2, &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let assignment = t01_zero_assignment(2);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 1);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let seed_tpsa = |state: &mut DescriptorComputedState| {
            for (sandp, value, rows) in [
                (false, 123.0, vec![1.0, 2.0]),
                (true, 456.0, vec![3.0, 4.0]),
            ] {
                let slot = state.slot_mut_for_tests(sandp);
                slot.scalar = Some(value);
                slot.contributions = Some(rows);
            }
        };
        let mut calls = 0u32;
        let mut check_call = |state: &mut DescriptorComputedState,
                              include_h: bool,
                              force: bool,
                              expected_h: bool| {
            let before = (
                topology.clone(),
                coordinates.clone(),
                properties.clone(),
                assignment.clone(),
                rings.clone(),
            );
            let before_counter = VALENCE_HELPER_ENTRIES.with(|c| c.get());
            let before_state = state.clone();
            let expected = COMPUTED[usize::from(expected_h)];
            calls += 1;
            let result = labute_contributions(&input, include_h, force, state).unwrap();
            assert_eq!(result.atoms.len(), 2);
            assert_eq!(
                result.atoms.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                vec![expected.0, expected.0]
            );
            assert_eq!(result.hydrogens.to_bits(), expected.1);
            let mut expected_state = before_state.clone();
            let slot = expected_state.labute_slot_mut();
            slot.rows = Some(vec![f64::from_bits(expected.0); 2]);
            slot.hydrogens = Some(f64::from_bits(expected.1));
            slot.asa = Some(f64::from_bits(expected.2));
            assert_eq!(*state, expected_state);
            assert_eq!(state.labute_slot().asa.map(f64::to_bits), Some(expected.2));
            assert_eq!(
                state.labute_slot().hydrogens.map(f64::to_bits),
                Some(expected.1)
            );
            assert_eq!(
                state
                    .labute_slot()
                    .rows
                    .as_ref()
                    .map(|v| v.iter().map(|v| v.to_bits()).collect::<Vec<_>>()),
                Some(vec![expected.0; 2])
            );
            if !force && before_state.labute_slot().rows.is_some() {
                assert_eq!(*state, before_state);
            }
            for sandp in [false, true] {
                assert_eq!(
                    state.slot_for_tests(sandp),
                    before_state.slot_for_tests(sandp)
                );
            }
            assert_eq!(
                (
                    topology.clone(),
                    coordinates.clone(),
                    properties.clone(),
                    assignment.clone(),
                    rings.clone()
                ),
                before
            );
            assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), before_counter);
        };
        for first_h in [false, true] {
            for second_h in [false, true] {
                let mut state = DescriptorComputedState::default();
                seed_tpsa(&mut state);
                check_call(&mut state, first_h, false, first_h);
                check_call(&mut state, second_h, false, first_h);
                let original_state = state.clone();
                let mut peer = state.clone();
                check_call(&mut peer, second_h, false, first_h);
                assert_eq!(state, original_state);
                check_call(&mut peer, second_h, true, second_h);
                assert_eq!(state, original_state);
                check_call(&mut peer, first_h, false, second_h);
                assert_eq!(state, original_state);
                peer.clear();
                assert_eq!(peer, DescriptorComputedState::default());
                assert_eq!(state, original_state);
                seed_tpsa(&mut peer);
                check_call(&mut peer, first_h, false, first_h);
                assert_eq!(state, original_state);
            }
        }
        assert_eq!(calls, 24);
    }

    #[test]
    fn descriptor_l03_contribution_presence_product_32() {
        // Exact 8x2x2 = 32-call presence product for labute_contributions.
        // Sentinel seeds prove warm arms SERVE STORED values and failure
        // arms are ATOMIC (slot unchanged); untouched TPSA checks prove
        // Labute never mutates the TPSA family's slots.
        const SENT_ROWS: [f64; 2] = [789.0, 789.0];
        const SENT_H: f64 = 123.0;
        const SENT_ASA: f64 = 456.0;
        const TPSA_ROWS: [f64; 2] = [999.0, 999.0];
        const TPSA_SCALAR: f64 = 888.0;
        const FALSE_ROW: u64 = 0x401b_ca6e_1564_c404;
        const FALSE_H: u64 = 0x0;
        const TRUE_ROW: u64 = 0x401b_b1e8_2a1b_1454;
        const TRUE_H: u64 = 0x3ff4_1b85_1ccc_a042;

        let bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let assignment = t01_zero_assignment(2);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 1);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);

        let mut calls = 0u32;
        for (rows_present, h_present, asa_present) in [
            (false, false, false),
            (false, false, true),
            (false, true, false),
            (false, true, true),
            (true, false, false),
            (true, false, true),
            (true, true, false),
            (true, true, true),
        ] {
            for force in [false, true] {
                for include_h in [false, true] {
                    calls += 1;
                    let mut state = DescriptorComputedState::default();
                    {
                        let slot = state.labute_slot_mut();
                        slot.rows = rows_present.then(|| SENT_ROWS.to_vec());
                        slot.hydrogens = h_present.then_some(SENT_H);
                        slot.asa = asa_present.then_some(SENT_ASA);
                    }
                    for sandp in [false, true] {
                        let tpsa = state.slot_mut_for_tests(sandp);
                        tpsa.contributions = Some(TPSA_ROWS.to_vec());
                        tpsa.scalar = Some(TPSA_SCALAR);
                    }
                    let (row, hydrogen) = if include_h {
                        (TRUE_ROW, TRUE_H)
                    } else {
                        (FALSE_ROW, FALSE_H)
                    };
                    let outcome = labute_contributions(&input, include_h, force, &mut state);
                    if force || !rows_present {
                        // Cold arm (guard bypassed or rows key absent):
                        // computed literals returned, ALL three fields
                        // (over)written with the computed values.
                        let result = outcome.expect("cold compute");
                        assert_eq!(result.atoms[0].to_bits(), row);
                        assert_eq!(result.atoms[1].to_bits(), row);
                        assert_eq!(result.hydrogens.to_bits(), hydrogen);
                        let slot = state.labute_slot();
                        assert_eq!(slot.rows.as_deref(), Some(result.atoms.as_slice()));
                        assert_eq!(
                            slot.hydrogens.map(f64::to_bits),
                            Some(result.hydrogens.to_bits())
                        );
                        assert!(slot.asa.is_some());
                    } else if !h_present {
                        // Warm rows hit, ordered read fails at the H
                        // property FIRST: MissingLabuteHydrogens (never
                        // Asa, never a recompute), slot ATOMIC.
                        assert_eq!(
                            outcome,
                            Err(DescriptorError::MissingLabuteHydrogens {
                                function: "labute_contributions"
                            })
                        );
                        let slot = state.labute_slot();
                        assert_eq!(slot.rows, Some(SENT_ROWS.to_vec()));
                        assert_eq!(slot.hydrogens, None);
                        assert_eq!(slot.asa, asa_present.then_some(SENT_ASA));
                    } else if !asa_present {
                        // Rows + H served, the ASA read fails next:
                        // MissingLabuteAsa, slot ATOMIC.
                        assert_eq!(
                            outcome,
                            Err(DescriptorError::MissingLabuteAsa {
                                function: "labute_contributions"
                            })
                        );
                        let slot = state.labute_slot();
                        assert_eq!(slot.rows, Some(SENT_ROWS.to_vec()));
                        assert_eq!(slot.hydrogens, Some(SENT_H));
                        assert_eq!(slot.asa, None);
                    } else {
                        // Fully warm: sentinel rows/H SERVED as stored.
                        let result = outcome.expect("warm serve");
                        assert_eq!(result.atoms, SENT_ROWS.to_vec());
                        assert_eq!(result.hydrogens.to_bits(), SENT_H.to_bits());
                        let slot = state.labute_slot();
                        assert_eq!(slot.rows, Some(SENT_ROWS.to_vec()));
                        assert_eq!(slot.hydrogens, Some(SENT_H));
                        assert_eq!(slot.asa, Some(SENT_ASA));
                    }
                    // Untouched TPSA checks: both TPSA family slots are
                    // bit-identical sentinel seeds after EVERY call.
                    for sandp in [false, true] {
                        let tpsa = state.slot_for_tests(sandp);
                        assert_eq!(tpsa.contributions, Some(TPSA_ROWS.to_vec()));
                        assert_eq!(tpsa.scalar, Some(TPSA_SCALAR));
                    }
                }
            }
        }
        assert_eq!(calls, 32);
    }

    #[test]
    fn descriptor_l03_scalar_presence_product_32() {
        // Exact 8x2x2 = 32-call presence product for labute_asa. The
        // scalar warm guard keys on the ASA presence alone (source
        // hasProp(_labuteASA)); the cold path dispatches through the ONE
        // contribution owner, whose ordered reads fix the error sequence.
        const SENT_ROWS: [f64; 2] = [789.0, 789.0];
        const SENT_H: f64 = 123.0;
        const SENT_ASA: f64 = 456.0;
        const TPSA_ROWS: [f64; 2] = [999.0, 999.0];
        const TPSA_SCALAR: f64 = 888.0;
        const FALSE_TOTAL: u64 = 0x402b_ca6e_1564_c404;
        const TRUE_TOTAL: u64 = 0x402e_3558_cdb4_a85c;

        let bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let assignment = t01_zero_assignment(2);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 1);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut calls = 0u32;
        for (rows_present, h_present, asa_present) in [
            (false, false, false),
            (false, false, true),
            (false, true, false),
            (false, true, true),
            (true, false, false),
            (true, false, true),
            (true, true, false),
            (true, true, true),
        ] {
            for force in [false, true] {
                for include_h in [false, true] {
                    calls += 1;
                    let mut state = DescriptorComputedState::default();
                    {
                        let slot = state.labute_slot_mut();
                        slot.rows = rows_present.then(|| SENT_ROWS.to_vec());
                        slot.hydrogens = h_present.then_some(SENT_H);
                        slot.asa = asa_present.then_some(SENT_ASA);
                    }
                    for sandp in [false, true] {
                        let tpsa = state.slot_mut_for_tests(sandp);
                        tpsa.contributions = Some(TPSA_ROWS.to_vec());
                        tpsa.scalar = Some(TPSA_SCALAR);
                    }
                    let total = if include_h { TRUE_TOTAL } else { FALSE_TOTAL };
                    let outcome = labute_asa(&input, include_h, force, &mut state);
                    if !force && asa_present {
                        assert_eq!(outcome.unwrap().to_bits(), SENT_ASA.to_bits());
                        let slot = state.labute_slot();
                        assert_eq!(slot.rows, rows_present.then(|| SENT_ROWS.to_vec()));
                        assert_eq!(slot.hydrogens, h_present.then_some(SENT_H));
                        assert_eq!(slot.asa, Some(SENT_ASA));
                    } else if force || !rows_present {
                        assert_eq!(outcome.unwrap().to_bits(), total);
                        let slot = state.labute_slot();
                        assert!(slot.rows.is_some());
                        assert!(slot.hydrogens.is_some());
                        assert!(slot.asa.is_some());
                    } else if !h_present {
                        assert_eq!(
                            outcome,
                            Err(DescriptorError::MissingLabuteHydrogens {
                                function: "labute_contributions"
                            })
                        );
                        let slot = state.labute_slot();
                        assert_eq!(slot.rows, Some(SENT_ROWS.to_vec()));
                        assert_eq!(slot.hydrogens, None);
                        assert_eq!(slot.asa, None);
                    } else {
                        assert_eq!(
                            outcome,
                            Err(DescriptorError::MissingLabuteAsa {
                                function: "labute_contributions"
                            })
                        );
                        let slot = state.labute_slot();
                        assert_eq!(slot.rows, Some(SENT_ROWS.to_vec()));
                        assert_eq!(slot.hydrogens, Some(SENT_H));
                        assert_eq!(slot.asa, None);
                    }
                    for sandp in [false, true] {
                        let tpsa = state.slot_for_tests(sandp);
                        assert_eq!(tpsa.contributions, Some(TPSA_ROWS.to_vec()));
                        assert_eq!(tpsa.scalar, Some(TPSA_SCALAR));
                    }
                }
            }
        }
        assert_eq!(calls, 32);
    }
    #[test]
    fn descriptor_l04_scalar_matches_complete_owner_output() {
        // Frozen L04 header (Step-364 audit + supervisor closure): the
        // C-C UNSPECIFIED bond (facts[0] = 0.1), sentinel family
        // rows[4,5]/H 7/ASA 789, DISTINCT TPSA seeds per slot, full
        // per-call input/counter/state bit assertions, and the in-body
        // calcLabuteASA anchors. Literals from the frozen enum table.
        const U_FALSE_ROW: u64 = 0x401b_14f3_0426_147b;
        const U_FALSE_H: u64 = 0x0;
        const U_FALSE_TOTAL: u64 = 0x402b_14f3_0426_147b;
        const U_TRUE_ROW: u64 = 0x401a_fc6d_18dc_64cb;
        const U_TRUE_H: u64 = 0x3ff4_1b85_1ccc_a042;
        const U_TRUE_TOTAL: u64 = 0x402d_7fdd_bc75_f8d3;
        const TPSA_CONV_SCALAR: f64 = 123.0;
        const TPSA_CONV_ROWS: [f64; 2] = [1.0, 2.0];
        const TPSA_SANDP_SCALAR: f64 = 456.0;
        const TPSA_SANDP_ROWS: [f64; 2] = [3.0, 4.0];

        let bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified);
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let assignment = t01_zero_assignment(2);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 1);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let topology_snapshot = topology.clone();
        let coordinates_snapshot = coordinates.clone();
        let assignment_snapshot = assignment.clone();
        let properties_snapshot = properties.clone();
        let rings_snapshot = rings.clone();

        for (include_h, row, hydrogen, total) in [
            (false, U_FALSE_ROW, U_FALSE_H, U_FALSE_TOTAL),
            (true, U_TRUE_ROW, U_TRUE_H, U_TRUE_TOTAL),
        ] {
            for force in [false, true] {
                let mut state = DescriptorComputedState::default();
                state.slot_mut_for_tests(false).scalar = Some(TPSA_CONV_SCALAR);
                state.slot_mut_for_tests(false).contributions = Some(TPSA_CONV_ROWS.to_vec());
                state.slot_mut_for_tests(true).scalar = Some(TPSA_SANDP_SCALAR);
                state.slot_mut_for_tests(true).contributions = Some(TPSA_SANDP_ROWS.to_vec());
                let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                // Cold scalar (fresh or force): frozen total bits, and
                // the owner's COMPLETE output beside it — by bits.
                let asa = labute_asa(&input, include_h, force, &mut state).unwrap();
                assert_eq!(asa.to_bits(), total);
                {
                    let slot = state.labute_slot();
                    let rows = slot.rows.as_ref().expect("rows published");
                    assert_eq!(rows.len(), 2);
                    assert_eq!(rows[0].to_bits(), row);
                    assert_eq!(rows[1].to_bits(), row);
                    assert_eq!(slot.hydrogens.expect("H published").to_bits(), hydrogen);
                    assert_eq!(slot.asa.expect("ASA published").to_bits(), total);
                }
                // Warm contributions on the SAME state: one owner run,
                // scalar + atom/H contributions served together.
                let served = labute_contributions(&input, include_h, false, &mut state).unwrap();
                assert_eq!(served.atoms[0].to_bits(), row);
                assert_eq!(served.atoms[1].to_bits(), row);
                assert_eq!(served.hydrogens.to_bits(), hydrogen);
                let asa2 = labute_asa(&input, include_h, false, &mut state).unwrap();
                assert_eq!(asa2.to_bits(), total);
                // Full per-call input/counter/state checks.
                assert_eq!(topology, topology_snapshot);
                assert_eq!(coordinates, coordinates_snapshot);
                assert_eq!(assignment, assignment_snapshot);
                assert_eq!(properties, properties_snapshot);
                assert_eq!(rings, rings_snapshot);
                assert_eq!(
                    VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                    helper_before,
                    "no valence-helper reentry"
                );
                assert_eq!(state.slot_for_tests(false).scalar, Some(TPSA_CONV_SCALAR));
                assert_eq!(
                    state.slot_for_tests(false).contributions,
                    Some(TPSA_CONV_ROWS.to_vec())
                );
                assert_eq!(state.slot_for_tests(true).scalar, Some(TPSA_SANDP_SCALAR));
                assert_eq!(
                    state.slot_for_tests(true).contributions,
                    Some(TPSA_SANDP_ROWS.to_vec())
                );
            }
        }
    }
    #[test]
    fn descriptor_l04_no_duplicated_area_algorithm() {
        // The scalar projection contains ZERO area arithmetic (verbatim
        // anchors inside labute_asa): with rows[4,5] and H 7 stored but
        // the ASA absent, labute_asa dispatches into the ONE owner whose
        // ordered ASA read fails — it never re-derives the total from
        // the stored rows. Clear() then empties every slot before a
        // fresh cold run.
        const SENT_ROWS: [f64; 2] = [4.0, 5.0];
        const SENT_H: f64 = 7.0;
        const SENT_ASA: f64 = 789.0;
        const TPSA_CONV_SCALAR: f64 = 123.0;
        const TPSA_CONV_ROWS: [f64; 2] = [1.0, 2.0];
        const TPSA_SANDP_SCALAR: f64 = 456.0;
        const TPSA_SANDP_ROWS: [f64; 2] = [3.0, 4.0];
        const U_FALSE_TOTAL: u64 = 0x402b_14f3_0426_147b;
        const U_FALSE_ROW: u64 = 0x401b_14f3_0426_147b;

        let bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified);
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let assignment = t01_zero_assignment(2);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 1);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let topology_snapshot = topology.clone();
        let coordinates_snapshot = coordinates.clone();
        let assignment_snapshot = assignment.clone();
        let properties_snapshot = properties.clone();
        let rings_snapshot = rings.clone();

        let mut state = DescriptorComputedState::default();
        {
            let slot = state.labute_slot_mut();
            slot.rows = Some(SENT_ROWS.to_vec());
            slot.hydrogens = Some(SENT_H);
        }
        state.slot_mut_for_tests(false).scalar = Some(TPSA_CONV_SCALAR);
        state.slot_mut_for_tests(false).contributions = Some(TPSA_CONV_ROWS.to_vec());
        state.slot_mut_for_tests(true).scalar = Some(TPSA_SANDP_SCALAR);
        state.slot_mut_for_tests(true).contributions = Some(TPSA_SANDP_ROWS.to_vec());
        let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
        assert_eq!(
            labute_asa(&input, false, false, &mut state),
            Err(DescriptorError::MissingLabuteAsa {
                function: "labute_contributions"
            })
        );
        {
            let slot = state.labute_slot();
            assert_eq!(slot.rows, Some(SENT_ROWS.to_vec()));
            assert_eq!(slot.hydrogens, Some(SENT_H));
            assert_eq!(slot.asa, None);
        }
        assert_eq!(topology, topology_snapshot);
        assert_eq!(coordinates, coordinates_snapshot);
        assert_eq!(assignment, assignment_snapshot);
        assert_eq!(properties, properties_snapshot);
        assert_eq!(rings, rings_snapshot);
        assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), helper_before);
        assert_eq!(state.slot_for_tests(false).scalar, Some(TPSA_CONV_SCALAR));
        assert_eq!(
            state.slot_for_tests(false).contributions,
            Some(TPSA_CONV_ROWS.to_vec())
        );
        assert_eq!(state.slot_for_tests(true).scalar, Some(TPSA_SANDP_SCALAR));
        assert_eq!(
            state.slot_for_tests(true).contributions,
            Some(TPSA_SANDP_ROWS.to_vec())
        );

        // Fully-seeded warm serve: stored ASA only, nothing recomputed.
        state.labute_slot_mut().asa = Some(SENT_ASA);
        assert_eq!(
            labute_asa(&input, true, false, &mut state)
                .unwrap()
                .to_bits(),
            SENT_ASA.to_bits()
        );
        {
            let slot = state.labute_slot();
            assert_eq!(slot.rows, Some(SENT_ROWS.to_vec()));
            assert_eq!(slot.hydrogens, Some(SENT_H));
            assert_eq!(slot.asa, Some(SENT_ASA));
        }

        // Clear sequence: every family slot empties, then a fresh cold
        // run republishes the complete owner output by bits.
        state.clear();
        {
            let slot = state.labute_slot();
            assert_eq!(slot.rows, None);
            assert_eq!(slot.hydrogens, None);
            assert_eq!(slot.asa, None);
        }
        assert_eq!(state.slot_for_tests(false).scalar, None);
        assert_eq!(state.slot_for_tests(false).contributions, None);
        assert_eq!(state.slot_for_tests(true).scalar, None);
        assert_eq!(state.slot_for_tests(true).contributions, None);
        let asa = labute_asa(&input, false, false, &mut state).unwrap();
        assert_eq!(asa.to_bits(), U_FALSE_TOTAL);
        {
            let slot = state.labute_slot();
            let rows = slot.rows.as_ref().expect("rows republished");
            assert_eq!(rows[0].to_bits(), U_FALSE_ROW);
            assert_eq!(rows[1].to_bits(), U_FALSE_ROW);
            assert_eq!(slot.hydrogens.expect("H republished").to_bits(), 0x0);
            assert_eq!(slot.asa.expect("ASA republished").to_bits(), U_FALSE_TOTAL);
        }
        assert_eq!(topology, topology_snapshot);
        assert_eq!(coordinates, coordinates_snapshot);
        assert_eq!(assignment, assignment_snapshot);
        assert_eq!(properties, properties_snapshot);
        assert_eq!(rings, rings_snapshot);
        assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), helper_before);
    }

    #[test]
    fn descriptor_c01_table_shape_and_label_sequence() {
        // Frozen C01 semantics (Step-372 audit): the ordered default
        // table has EXACTLY 110 rows; the label sequence below is the
        // literal source order (72 labels, repeats = per-label counts),
        // including O12 before O7 and S2 before S1 before S3. idx is
        // the 0-based position. Nothing here derives from element
        // arithmetic.
        const EXPECTED_LABELS: [&str; 110] = [
            "C1", "C1", "C1", "C2", "C2", "C3", "C3", "C4", "C4", "C5", "C6", "C6", "C6", "C6",
            "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19",
            "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C26", "C26", "C26", "C27", "CS",
            "H1", "H2", "H2", "H2", "H3", "H3", "H4", "H4", "HS", "N1", "N2", "N3", "N4", "N5",
            "N6", "N7", "N8", "N8", "N9", "N10", "N11", "N12", "N13", "N13", "N13", "N14", "N14",
            "N14", "NS", "O1", "O2", "O3", "O4", "O5", "O5", "O6", "O6", "O12", "O7", "O8", "O9",
            "O9", "O9", "O9", "O9", "O10", "O10", "O10", "O11", "OS", "F", "Cl", "Br", "I", "Hal",
            "Hal", "Hal", "P", "S2", "S2", "S1", "S3", "Me1", "Me1", "Me1", "Me1", "Me1", "Me1",
            "Me2", "Me2", "Me2",
        ];
        let rows = default_crippen_params();
        assert_eq!(rows.len(), 110);
        for (position, (row, expected)) in rows.iter().zip(EXPECTED_LABELS).enumerate() {
            assert_eq!(row.label, expected, "position {position}");
            assert_eq!(row.idx, position as u32);
            assert!(!row.smarts.is_empty(), "smarts present at {position}");
        }
    }

    #[test]
    fn descriptor_c01_exact_fields_and_order_flips() {
        // Frozen literal source rows at the audited positions: the C4
        // row split across two C++ literals, the N10 blank-MR row, the
        // O12-before-O7 order flip boundary (O12 blank MR, O7 literal
        // zero), the Hal blank-MR rows, the S2/S1/S3 flip, and the final
        // Me2 row. All values are literal source cells.
        let rows = default_crippen_params();
        // Row 0: first C1.
        assert_eq!(rows[0].label, "C1");
        assert_eq!(rows[0].smarts, "[CH4]");
        assert_eq!(rows[0].logp, 0.1441);
        assert_eq!(rows[0].mr, 2.503);
        // Row 8: the C4 row assembled from two C++ string literals.
        assert_eq!(rows[8].label, "C4");
        assert_eq!(
            rows[8].smarts,
            "[CH0X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])([A;!#1])[A;!#1]"
        );
        assert_eq!(rows[8].logp, -0.2051);
        assert_eq!(rows[8].mr, 2.731);
        // Row 58: N10 blank MR.
        assert_eq!(rows[58].label, "N10");
        assert_eq!(rows[58].smarts, "[NH3,NH2,NH;+,+2,+3]");
        assert_eq!(rows[58].logp, -1.95);
        assert_eq!(rows[58].mr, 0.0);
        // Rows 76/77: O12 BEFORE O7 (intentional flip); O12's MR cell is
        // blank (0.0), O7's MR cell is a literal 0.
        assert_eq!(rows[76].label, "O12");
        assert_eq!(rows[76].smarts, "[O-]C(=O)");
        assert_eq!(rows[76].logp, -1.326);
        assert_eq!(rows[76].mr, 0.0);
        assert_eq!(rows[77].label, "O7");
        assert_eq!(rows[77].smarts, "[OX1;-,-2,-3][!#1;!N;!S]");
        assert_eq!(rows[77].logp, -1.189);
        assert_eq!(rows[77].mr, 0.0);
        // Rows 93-95: the three Hal rows, blank MR each.
        for (i, smarts) in [
            "[#9,#17,#35,#53;-]",
            "[#53;+,+2,+3]",
            "[+;#3,#11,#19,#37,#55]",
        ]
        .into_iter()
        .enumerate()
        {
            assert_eq!(rows[93 + i].label, "Hal");
            assert_eq!(rows[93 + i].smarts, smarts);
            assert_eq!(rows[93 + i].logp, -2.996);
            assert_eq!(rows[93 + i].mr, 0.0);
        }
        // Rows 97/99/100: S2, then S1, then S3 (intentional flip).
        assert_eq!(rows[97].label, "S2");
        assert_eq!(rows[97].smarts, "[S;-,-2,-3,-4,+1,+2,+3,+5,+6]");
        assert_eq!(rows[97].logp, -0.0024);
        assert_eq!(rows[97].mr, 7.365);
        assert_eq!(rows[98].label, "S2");
        assert_eq!(rows[99].label, "S1");
        assert_eq!(rows[99].smarts, "[S;A]");
        assert_eq!(rows[99].logp, 0.6482);
        assert_eq!(rows[99].mr, 7.591);
        assert_eq!(rows[100].label, "S3");
        assert_eq!(rows[100].smarts, "[s;a]");
        assert_eq!(rows[100].logp, 0.6237);
        assert_eq!(rows[100].mr, 6.691);
        // Row 109: the final Me2 row, blank MR.
        assert_eq!(rows[109].label, "Me2");
        assert_eq!(rows[109].smarts, "[#72,#73,#74,#75,#76,#77,#78,#79,#80]");
        assert_eq!(rows[109].logp, -0.0025);
        assert_eq!(rows[109].mr, 0.0);
    }

    #[test]
    fn descriptor_c01_blank_mr_cells_and_asset_text() {
        // The nine blank-MR rows parse to 0.0 with their logP intact,
        // and the committed asset carries the exact header line, 111
        // lines, no data line starting '#', and the source notes.
        let rows = default_crippen_params();
        let blank_mr: Vec<(usize, &str, f64)> = rows
            .iter()
            .enumerate()
            .filter(|(_, r)| {
                crate::crippen::DEFAULT_PARAM_DATA
                    .lines()
                    .nth(1 + r.idx as usize)
                    .is_some_and(|line| line.split('\t').nth(3) == Some(""))
            })
            .map(|(i, r)| (i, r.label, r.logp))
            .collect();
        let expected_blank: [(usize, &str, f64); 9] = [
            (58, "N10", -1.95),
            (60, "N12", -1.119),
            (76, "O12", -1.326),
            (93, "Hal", -2.996),
            (94, "Hal", -2.996),
            (95, "Hal", -2.996),
            (107, "Me2", -0.0025),
            (108, "Me2", -0.0025),
            (109, "Me2", -0.0025),
        ];
        assert_eq!(blank_mr, expected_blank);
        for (position, label, _) in expected_blank {
            assert_eq!(rows[position].label, label);
            assert_eq!(rows[position].mr.to_bits(), 0, "blank MR at {position}");
        }
        let asset = crate::crippen::DEFAULT_PARAM_DATA;
        assert!(asset.starts_with("#ID\tSMARTS\tlogP\tMR\t"));
        assert_eq!(
            asset.lines().next(),
            Some("#ID\tSMARTS\tlogP\tMR\tNotes/Questions")
        );
        // Frozen independently from pinned defaultParamData before this edit.
        // FNV-1a is an asset-change regression tripwire, not a claim of
        // cryptographic identity or a substitute for source review.
        assert_eq!(asset.len(), 3912);
        let checksum = asset.bytes().fold(0xcbf29ce484222325u64, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3)
        });
        assert_eq!(checksum, 0xf5a1a2441e4d64b6);
        let lines: Vec<&str> = asset.lines().collect();
        assert_eq!(lines.len(), 111);
        assert!(lines[1..].iter().all(|l| !l.starts_with('#')));
        assert!(asset.contains("order flip here intentional"));
        assert!(asset.contains("Expanded definition of (pseudo-)ionic S"));
        assert!(asset.contains("Footnote h indicates these should be here?"));
    }

    #[test]
    fn descriptor_c02_numeric_cell_product() {
        // Frozen C02 semantics (Step-380 audit): LogP has NO fallback
        // (malformed -> typed error, empty -> 0.0); the fallback exists
        // ONLY on the MR side (empty OR malformed -> 0.0). Complete
        // small cell product over the same cells for both fields.
        for cell in ["", "0", "-1.95", "0.1441", "14.02", "abc", "1.5x"] {
            let logp = crate::crippen::parse_crippen_logp(cell);
            let mr = crate::crippen::parse_crippen_mr(cell);
            match cell {
                "" => {
                    assert_eq!(logp, Ok(0.0), "logp empty cell");
                    assert_eq!(mr.to_bits(), 0.0_f64.to_bits(), "mr empty cell");
                }
                "0" => {
                    assert_eq!(logp, Ok(0.0));
                    assert_eq!(mr.to_bits(), 0.0_f64.to_bits());
                }
                "-1.95" => {
                    assert_eq!(logp, Ok(-1.95));
                    assert_eq!(mr, -1.95);
                }
                "0.1441" => {
                    assert_eq!(logp, Ok(0.1441));
                    assert_eq!(mr, 0.1441);
                }
                "14.02" => {
                    assert_eq!(logp, Ok(14.02));
                    assert_eq!(mr, 14.02);
                }
                malformed => {
                    // LogP: typed failure preserving the offending cell —
                    // NEVER the 0.0 fallback.
                    assert_eq!(
                        logp,
                        Err(DescriptorError::CrippenParamNumeric {
                            field: "logp",
                            cell: malformed.to_string(),
                        }),
                        "malformed logp cell {malformed:?}"
                    );
                    assert_ne!(logp, Ok(0.0));
                    // MR: the ONLY source fallback applies.
                    assert_eq!(mr.to_bits(), 0.0_f64.to_bits());
                }
            }
        }
    }

    #[test]
    fn descriptor_c02_committed_table_avoids_error_path() {
        // Every committed asset LogP cell parses Ok (the typed error
        // path is unreachable for the pinned table); MR cells are all
        // finite. Uses the raw asset cells, not the parsed rows.
        let asset = crate::crippen::DEFAULT_PARAM_DATA;
        for (line_no, line) in asset.lines().enumerate().skip(1) {
            let cells: Vec<&str> = line.split('\t').collect();
            let logp_cell = cells[2];
            let mr_cell = cells[3];
            assert!(
                crate::crippen::parse_crippen_logp(logp_cell).is_ok(),
                "line {line_no}"
            );
            assert!(crate::crippen::parse_crippen_mr(mr_cell).is_finite());
        }
    }

    #[test]
    fn descriptor_c03_default_patterns_compiled_once_and_retained() {
        // Frozen C03 semantics (Step-388 audit): the DEFAULT collection
        // compiles ALL 110 ordered patterns exactly once (flyweight
        // construct-once) and retains them; every later acquisition
        // borrows the SAME graphs and compiles nothing. The pinned
        // table's SMARTS are all valid — every entry is Some.
        let rows = default_crippen_params();
        let first = crate::crippen::crippen_patterns();
        assert_eq!(first.len(), rows.len());
        assert!(first.len() == 110);
        assert!(
            first.iter().all(|p| p.is_some()),
            "every pinned Crippen SMARTS compiles"
        );
        let before = crate::crippen::CRIPPEN_SMARTS_COMPILES.with(|count| count.get());
        let second = crate::crippen::crippen_patterns();
        for (index, (a, b)) in first.iter().zip(second).enumerate() {
            let (a, b) = (a.as_ref().unwrap(), b.as_ref().unwrap());
            assert!(
                std::sync::Arc::ptr_eq(a, b),
                "retained graph at {index} is the same allocation"
            );
        }
        let after = crate::crippen::CRIPPEN_SMARTS_COMPILES.with(|count| count.get());
        assert_eq!(after, before, "no reparse on repeated default acquisition");
    }

    #[test]
    fn descriptor_c03_custom_ownership_and_null_retention() {
        // Frozen C03 semantics: a CUSTOM table builds a caller-owned
        // collection (default OnceLock untouched); a row whose SMARTS
        // fails to compile RETAINS None — the source's unchecked
        // SmartsToMol nullptr — without panicking or erroring at
        // acquisition.
        let custom =
            "#ID\tSMARTS\tlogP\tMR\tNotes/Questions\nX1\t[CH4]\t0.1\t0.2\t\nX2\t[QQ\t0.1\t0.2\t\n";
        let before = crate::crippen::CRIPPEN_SMARTS_COMPILES.with(|count| count.get());
        let custom_patterns = crate::crippen::crippen_patterns_from(custom);
        let after = crate::crippen::CRIPPEN_SMARTS_COMPILES.with(|count| count.get());
        assert_eq!(custom_patterns.len(), 2);
        assert!(custom_patterns[0].is_some(), "valid custom row compiles");
        assert!(
            custom_patterns[1].is_none(),
            "malformed custom SMARTS retains the null (unchecked SmartsToMol)"
        );
        assert_eq!(after - before, 2, "one compile per custom row, no more");
        // The default collection is independent and unchanged. Build the
        // default FIRST so the later borrow measures compile-free
        // regardless of test execution order (the OnceLock is
        // process-global; the counter is per-thread).
        let _ = crate::crippen::crippen_patterns();
        let before_default = crate::crippen::CRIPPEN_SMARTS_COMPILES.with(|count| count.get());
        let default = crate::crippen::crippen_patterns();
        assert_eq!(default.len(), 110);
        assert!(default.iter().all(|p| p.is_some()));
        let after2 = crate::crippen::CRIPPEN_SMARTS_COMPILES.with(|count| count.get());
        assert_eq!(after2, before_default, "default borrow compiles nothing");
    }

    #[test]
    fn descriptor_c03_headerless_and_comment_order_product() {
        // Frozen C03-ROWS product (plan TOP): 3 prefix states x2 row
        // orders x2 inter-row comment states = 12 real
        // crippen_patterns_from calls. Expectations are the LITERAL
        // specified order — never actual output. Every table ends
        // with a newline (EOF behavior stays out of this regression).
        const CARBON_ROW: &str = "C\t[CH4]\t0.1\t0.2\t\n";
        const OXYGEN_ROW: &str = "O\t[#8]\t0.3\t0.4\t\n";
        const PREFIXES: [&str; 3] = ["", "#header\n", "#header\n#before\n"];
        let mut calls = 0u32;
        let mut total_compiles = 0u32;
        for prefix in PREFIXES {
            for (first_row, second_row) in [(CARBON_ROW, OXYGEN_ROW), (OXYGEN_ROW, CARBON_ROW)] {
                for middle in ["", "#middle\n"] {
                    let table = format!("{prefix}{first_row}{middle}{second_row}");
                    let before = crate::crippen::CRIPPEN_SMARTS_COMPILES.with(|count| count.get());
                    let patterns = crate::crippen::crippen_patterns_from(&table);
                    let after = crate::crippen::CRIPPEN_SMARTS_COMPILES.with(|count| count.get());
                    calls += 1;
                    // Exactly TWO Some graphs, in the literal row order
                    // (prefix comments skipped; headerless first row
                    // KEPT; inter-row comment skipped).
                    assert_eq!(patterns.len(), 2, "table {table:?}");
                    let first = patterns[0]
                        .as_ref()
                        .expect("first row graph (headerless first row kept)");
                    let second = patterns[1].as_ref().expect("second row graph");
                    let first_element = if first_row == CARBON_ROW {
                        Element::C
                    } else {
                        Element::O
                    };
                    let second_element = if first_row == CARBON_ROW {
                        Element::O
                    } else {
                        Element::C
                    };
                    assert_eq!(first.atoms()[0].element(), Some(first_element));
                    assert_eq!(second.atoms()[0].element(), Some(second_element));
                    // Per-call compile bracket: exactly 2 compiles.
                    assert_eq!(
                        after - before,
                        2,
                        "one compile per real row, table {table:?}"
                    );
                    total_compiles += 2;
                }
            }
        }
        assert_eq!(calls, 12);
        assert_eq!(total_compiles, 24);
    }

    #[test]
    fn descriptor_c04_first_match_typing_and_untyped() {
        // Frozen C04 semantics (Step-396 audit), all expectations literal
        // source-table values: FIRST table row wins; only the first query
        // atom of a match is typed; untyped atoms keep 0.0/0/""; the loop
        // stops once every atom is typed.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();

        // Methane: the carbon is typed by row 0 `[CH4]` (C1) even though
        // `CS [#6]` (row 38) also matches — FIRST TABLE MATCH.
        let topology = smiles_topology("C");
        let assignment = valence(&topology, "c04").unwrap();
        let rings = ring_info(&topology, "c04").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut st_c = tpsa::DescriptorComputedState::default();
        let mut types = vec![0u32; input.topology().atoms.len()];
        let mut labels = vec![String::new(); input.topology().atoms.len()];
        let c = crippen_contributions(&input, true, &mut st_c, Some(&mut types), Some(&mut labels))
            .unwrap();
        assert_eq!(c.logp.len(), 1);
        assert_eq!(labels[0], "C1");
        assert_eq!(types[0], 0);
        assert_eq!(c.logp[0], 0.1441);
        assert_eq!(c.molar_refractivity[0], 2.503);

        // Ethane: BOTH carbons typed by row 1 `[CH3]C` (C1) — all atoms
        // typed by early rows, so `CS` never runs (early stop).
        let topology = smiles_topology("CC");
        let assignment = valence(&topology, "c04").unwrap();
        let rings = ring_info(&topology, "c04").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut st_c = tpsa::DescriptorComputedState::default();
        let mut types = vec![0u32; input.topology().atoms.len()];
        let mut labels = vec![String::new(); input.topology().atoms.len()];
        let c = crippen_contributions(&input, true, &mut st_c, Some(&mut types), Some(&mut labels))
            .unwrap();
        for i in 0..2 {
            assert_eq!(labels[i], "C1", "atom {i}");
            assert_eq!(c.logp[i], 0.1441);
            assert_eq!(c.molar_refractivity[i], 2.503);
        }

        // Methanol: the carbon is typed by C3 `[CH3][N,O,...]` — typing
        // lands on the FIRST QUERY ATOM ([CH3], the carbon), not the
        // oxygen atom of the match; the oxygen itself is typed O2.
        let topology = smiles_topology("CO");
        let assignment = valence(&topology, "c04").unwrap();
        let rings = ring_info(&topology, "c04").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut st_c = tpsa::DescriptorComputedState::default();
        let mut types = vec![0u32; input.topology().atoms.len()];
        let mut labels = vec![String::new(); input.topology().atoms.len()];
        let c = crippen_contributions(&input, true, &mut st_c, Some(&mut types), Some(&mut labels))
            .unwrap();
        assert_eq!(labels[0], "C3");
        assert_eq!(types[0], 5);
        assert_eq!(c.logp[0], -0.2035);
        assert_eq!(c.molar_refractivity[0], 2.753);
        assert_eq!(labels[1], "O2");
        assert_eq!(c.logp[1], -0.2893);
        assert_eq!(c.molar_refractivity[1], 0.8238);

        // Lone Curium: NO table row matches — untyped entries keep the
        // initial 0.0 / 0 / "" values.
        let atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::CM))];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(1);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 1, 0);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut st_c = tpsa::DescriptorComputedState::default();
        let mut types = vec![0u32; input.topology().atoms.len()];
        let mut labels = vec![String::new(); input.topology().atoms.len()];
        let c = crippen_contributions(&input, true, &mut st_c, Some(&mut types), Some(&mut labels))
            .unwrap();
        assert_eq!(c.logp[0].to_bits(), 0.0_f64.to_bits());
        assert_eq!(c.molar_refractivity[0].to_bits(), 0.0_f64.to_bits());
        assert_eq!(types[0], 0);
        assert_eq!(labels[0], "");
    }

    #[test]
    fn descriptor_c02_numeric_source_product() {
        // Frozen C02-NUMERIC product (plan TOP): 12 literal bases x3
        // signs x3 outer affixes x2 field routes = 216 calls with
        // literal source-bit expectations only.
        const BASES: [(&str, Option<u64>); 12] = [
            ("0", Some(0x0000_0000_0000_0000)),
            ("1.5", Some(0x3ff8_0000_0000_0000)),
            ("1e3", Some(0x408f_4000_0000_0000)),
            ("1e309", None),
            ("1e-999", Some(0x0000_0000_0000_0000)),
            ("inf", Some(0x7ff0_0000_0000_0000)),
            ("InFiNiTy", Some(0x7ff0_0000_0000_0000)),
            ("nan", Some(0x7ff8_0000_0000_0000)),
            ("nan()", Some(0x7ff8_0000_0000_0000)),
            ("nAn(foo)", Some(0x7ff8_0000_0000_0000)),
            ("abc", None),
            ("1.5x", None),
        ];
        let mut logp_calls = 0u32;
        let mut mr_calls = 0u32;
        for (base, expected_bits) in BASES {
            for sign in ["", "+", "-"] {
                for (lead, trail) in [("", ""), (" ", ""), ("", " ")] {
                    let cell = format!("{lead}{sign}{base}{trail}");
                    let logp = crate::crippen::parse_crippen_logp(&cell);
                    let mr = crate::crippen::parse_crippen_mr(&cell);
                    logp_calls += 1;
                    mr_calls += 1;
                    let accepted = if lead.is_empty() && trail.is_empty() {
                        expected_bits.map(|bits| bits | u64::from(sign == "-") << 63)
                    } else {
                        None
                    };
                    match accepted {
                        Some(bits) => {
                            let v = logp.expect("accepted LogP");
                            assert_eq!(v.to_bits(), bits, "cell {cell:?}");
                            assert_eq!(mr.to_bits(), bits, "cell {cell:?}");
                            if bits & 0x7fff_ffff_ffff_ffff == 0x7ff8_0000_0000_0000 {
                                assert!(v.is_nan(), "cell {cell:?}");
                                assert!(mr.is_nan(), "cell {cell:?}");
                            }
                        }
                        None => {
                            assert_eq!(
                                logp,
                                Err(DescriptorError::CrippenParamNumeric {
                                    field: "logp",
                                    cell: cell.clone(),
                                }),
                                "cell {cell:?}"
                            );
                            assert_eq!(mr.to_bits(), 0.0_f64.to_bits(), "cell {cell:?}");
                        }
                    }
                }
            }
        }
        assert_eq!(logp_calls, 108);
        assert_eq!(mr_calls, 108);
        assert_eq!(logp_calls + mr_calls, 216);
    }

    #[test]
    fn descriptor_c02_numeric_source_boundaries() {
        // Frozen: 14 literal cells x2 routes = 28 calls.
        const CELLS: [(&str, Option<u64>); 14] = [
            ("", Some(0x0000_0000_0000_0000)),
            ("1e", None),
            ("1e+", None),
            ("1e-", None),
            (".", None),
            ("+", None),
            ("-", None),
            ("0x1p0", None),
            ("1,5", None),
            ("1_5", None),
            ("1.7976931348623157e308", Some(0x7fef_ffff_ffff_ffff)),
            ("1.7976931348623159e308", None),
            ("2.2250738585072014e-308", Some(0x0010_0000_0000_0000)),
            ("5e-324", Some(0x0000_0000_0000_0001)),
        ];
        let mut calls = 0u32;
        for (cell, expected_bits) in CELLS {
            let logp = crate::crippen::parse_crippen_logp(cell);
            let mr = crate::crippen::parse_crippen_mr(cell);
            calls += 2;
            match expected_bits {
                Some(bits) => {
                    assert_eq!(
                        logp.expect("accepted LogP").to_bits(),
                        bits,
                        "cell {cell:?}"
                    );
                    assert_eq!(mr.to_bits(), bits, "cell {cell:?}");
                }
                None => {
                    assert_eq!(
                        logp,
                        Err(DescriptorError::CrippenParamNumeric {
                            field: "logp",
                            cell: cell.to_string(),
                        }),
                        "cell {cell:?}"
                    );
                    assert_eq!(mr.to_bits(), 0.0_f64.to_bits(), "cell {cell:?}");
                }
            }
        }
        assert_eq!(calls, 28);
    }

    #[test]
    fn descriptor_c02_numeric_source_nan_payload_bytes() {
        // Frozen: four payload literals x3 signs x2 routes = 24 calls,
        // all ACCEPTED as the canonical quiet NaN with exact sign bits;
        // payload bytes are preserved verbatim, never interpreted.
        const PAYLOADS: [&str; 4] = ["nan(\u{3bb})", "nan((x))", "nan( a )", "nan(\0)"];
        let mut calls = 0u32;
        for payload in PAYLOADS {
            for sign in ["", "+", "-"] {
                let cell = format!("{sign}{payload}");
                let logp = crate::crippen::parse_crippen_logp(&cell);
                let mr = crate::crippen::parse_crippen_mr(&cell);
                calls += 2;
                let expected = if sign == "-" {
                    0xfff8_0000_0000_0000
                } else {
                    0x7ff8_0000_0000_0000
                };
                let v = logp.expect("nan(payload) accepted");
                assert!(v.is_nan(), "cell {cell:?}");
                assert_eq!(v.to_bits(), expected, "cell {cell:?}");
                assert!(mr.is_nan(), "cell {cell:?}");
                assert_eq!(mr.to_bits(), expected, "cell {cell:?}");
            }
        }
        assert_eq!(calls, 24);
    }

    #[test]
    fn descriptor_c05_type_index_label_projection() {
        // Frozen C05 semantics (Step-404 audit), all expectations literal
        // SOURCE table values: the type index is the table ROW index
        // (param.idx) and the label the row label, written in the same
        // first-match order as logp/mr; untyped atoms keep the initial
        // 0/"" values; the owned boundary always returns all four
        // length-nAtoms vectors.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();

        // Methane: FIRST table row [CH4] = idx 0, label "C1".
        let topology = smiles_topology("C");
        let assignment = valence(&topology, "c05").unwrap();
        let rings = ring_info(&topology, "c05").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut st_c = tpsa::DescriptorComputedState::default();
        let mut types = vec![0u32; input.topology().atoms.len()];
        let mut labels = vec![String::new(); input.topology().atoms.len()];
        let c = crippen_contributions(&input, true, &mut st_c, Some(&mut types), Some(&mut labels))
            .unwrap();
        assert_eq!(types, vec![0u32]);
        assert_eq!(labels, vec!["C1"]);
        assert_eq!(c.logp.len(), 1);
        assert_eq!(c.molar_refractivity.len(), 1);

        // Ethane: both carbons by row 1 [CH3]C — label "C1" but ROW
        // INDEX 1 (the typed projection follows the same FIRST-match
        // rule: [CH4] row 0 does not match a bonded carbon).
        let topology = smiles_topology("CC");
        let assignment = valence(&topology, "c05").unwrap();
        let rings = ring_info(&topology, "c05").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut st_c = tpsa::DescriptorComputedState::default();
        let mut types = vec![0u32; input.topology().atoms.len()];
        let mut labels = vec![String::new(); input.topology().atoms.len()];
        let c = crippen_contributions(&input, true, &mut st_c, Some(&mut types), Some(&mut labels))
            .unwrap();
        assert_eq!(types, vec![1u32, 1]);
        assert_eq!(labels, vec!["C1", "C1"]);

        // Methanol: C typed by C3 (row 5, label "C3"), O typed by O2
        // (row 69, label "O2") — the ORDERING discriminator: different
        // atoms carry different first-match rows in one call.
        let topology = smiles_topology("CO");
        let assignment = valence(&topology, "c05").unwrap();
        let rings = ring_info(&topology, "c05").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut st_c = tpsa::DescriptorComputedState::default();
        let mut types = vec![0u32; input.topology().atoms.len()];
        let mut labels = vec![String::new(); input.topology().atoms.len()];
        let c = crippen_contributions(&input, true, &mut st_c, Some(&mut types), Some(&mut labels))
            .unwrap();
        assert_eq!(types, vec![5u32, 69]);
        assert_eq!(labels, vec!["C3", "O2"]);

        // Lone Curium: NO row matches — types stay 0, labels stay ""
        // (the caller-initial values the source leaves untouched).
        let atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::CM))];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(1);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 1, 0);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut st_c = tpsa::DescriptorComputedState::default();
        let mut types = vec![0u32; input.topology().atoms.len()];
        let mut labels = vec![String::new(); input.topology().atoms.len()];
        let c = crippen_contributions(&input, true, &mut st_c, Some(&mut types), Some(&mut labels))
            .unwrap();
        assert_eq!(types, vec![0u32]);
        assert_eq!(labels, vec![""]);
    }

    #[test]
    fn descriptor_c06_totals_ordered_contributions() {
        // Frozen C06 semantics (Step-412 audit): totals by EXACT BITS
        // from the frozen C04/C05 row literals summed left-to-right in
        // ATOM ORDER (the source's iterator accumulation). The totals
        // owner operates on the molecule AS GIVEN (includeHs=false).
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();

        // Methane: 0.1441 + 2.503 (single atom).
        let topology = smiles_topology("C");
        let assignment = valence(&topology, "c06").unwrap();
        let rings = ring_info(&topology, "c06").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut fresh_st = tpsa::DescriptorComputedState::default();
        let totals_c = crippen_totals(&input, false, false, &mut fresh_st).unwrap();
        let (logp, mr) = (totals_c.logp, totals_c.molar_refractivity);
        assert_eq!(logp.to_bits(), 0.1441_f64.to_bits());
        assert_eq!(mr.to_bits(), 2.503_f64.to_bits());

        // Ethane: 0.1441 + 0.1441 (row [CH3]C for both carbons) and
        // 2.503 + 2.503 — left-to-right in index order.
        let topology = smiles_topology("CC");
        let assignment = valence(&topology, "c06").unwrap();
        let rings = ring_info(&topology, "c06").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut fresh_st = tpsa::DescriptorComputedState::default();
        let totals_c = crippen_totals(&input, false, false, &mut fresh_st).unwrap();
        let (logp, mr) = (totals_c.logp, totals_c.molar_refractivity);
        let expected_logp = 0.0_f64 + 0.1441 + 0.1441;
        let expected_mr = 0.0_f64 + 2.503 + 2.503;
        assert_eq!(logp.to_bits(), expected_logp.to_bits());
        assert_eq!(mr.to_bits(), expected_mr.to_bits());

        // Methanol: C3 (-0.2035) + O2 (-0.2893); 2.753 + 0.8238.
        let topology = smiles_topology("CO");
        let assignment = valence(&topology, "c06").unwrap();
        let rings = ring_info(&topology, "c06").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut fresh_st = tpsa::DescriptorComputedState::default();
        let totals_c = crippen_totals(&input, false, false, &mut fresh_st).unwrap();
        let (logp, mr) = (totals_c.logp, totals_c.molar_refractivity);
        let expected_logp = 0.0_f64 + -0.2035 + -0.2893;
        let expected_mr = 0.0_f64 + 2.753 + 0.8238;
        assert_eq!(logp.to_bits(), expected_logp.to_bits());
        assert_eq!(mr.to_bits(), expected_mr.to_bits());

        // Lone Curium: NO row matches -> both totals exactly +0.0 bits.
        let atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::CM))];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(1);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 1, 0);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut fresh_st = tpsa::DescriptorComputedState::default();
        let totals_c = crippen_totals(&input, false, false, &mut fresh_st).unwrap();
        let (logp, mr) = (totals_c.logp, totals_c.molar_refractivity);
        assert_eq!(logp.to_bits(), 0.0_f64.to_bits());
        assert_eq!(mr.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn descriptor_c06_force_and_idempotence_product() {
        // Frozen C06 sequences (Step-412 audit): force false/true x
        // repeated/clone calls. Every arm below carries a FRESH state,
        // so each is a cold computation and the frozen observable is:
        // every call returns bit-identical totals regardless of force
        // or repetition (idempotence), matching the source's cold path
        // (the shared-key scalar cache itself is the canonical
        // crippen_totals guard's arm, exercised by the C08 selectors).
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("CCO");
        let assignment = valence(&topology, "c06p").unwrap();
        let rings = ring_info(&topology, "c06p").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());

        let mut fresh_st = tpsa::DescriptorComputedState::default();
        let totals_c = crippen_totals(&input, false, false, &mut fresh_st).unwrap();
        let (first_logp, first_mr) = (totals_c.logp, totals_c.molar_refractivity);
        for force in [false, true] {
            let mut fresh_st = tpsa::DescriptorComputedState::default();
            let totals_c = crippen_totals(&input, false, force, &mut fresh_st).unwrap();
            let (logp, mr) = (totals_c.logp, totals_c.molar_refractivity);
            assert_eq!(logp.to_bits(), first_logp.to_bits(), "force={force}");
            assert_eq!(mr.to_bits(), first_mr.to_bits(), "force={force}");
            // Repeated call: bit-identical (ordered cold idempotence).
            let mut fresh_st2 = tpsa::DescriptorComputedState::default();
            let totals_c2 = crippen_totals(&input, false, force, &mut fresh_st2).unwrap();
            let (logp2, mr2) = (totals_c2.logp, totals_c2.molar_refractivity);
            assert_eq!(logp2.to_bits(), first_logp.to_bits());
            assert_eq!(mr2.to_bits(), first_mr.to_bits());
        }
        // Clone semantics: the stateless-totals identity — a clone of
        // every input block produces identical totals (no hidden state).
        let assignment2 = assignment.clone();
        let rings2 = rings.clone();
        let input2 =
            DescriptorInput::new(&topology, &coordinates, &properties, &assignment2, &rings2);
        let mut fresh_st = tpsa::DescriptorComputedState::default();
        let totals_c = crippen_totals(&input2, false, false, &mut fresh_st).unwrap();
        let (logp, mr) = (totals_c.logp, totals_c.molar_refractivity);
        assert_eq!(logp.to_bits(), first_logp.to_bits());
        assert_eq!(mr.to_bits(), first_mr.to_bits());
        assert_eq!(
            VALENCE_HELPER_ENTRIES.with(|c| c.get()),
            helper_before,
            "no valence-helper reentry across all totals calls"
        );
    }

    #[test]
    fn descriptor_c07_totals_with_hydrogens() {
        // Frozen C07 semantics (Step-420 audit): includeHs=true applies
        // the canonical core AddHs owner with (false,false) on a working
        // copy, types through the H-family table rows, and accumulates
        // ordered — asserted by exact bits from the frozen table rows
        // (never element arithmetic). The ORIGINAL input blocks must be
        // untouched (the source never harms the molecule).
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();

        // Ethane "CC": both carbons are CH3 groups (3 implicit H each).
        // With Hs added, each carbon types by C1-row1 [CH3]C (idx 1,
        // logp 0.1441 / mr 2.503) and each of the 6 H nodes types by
        // H1 [#1][#6,#1] (idx 39, logp 0.123 / mr 1.057).
        let topology = smiles_topology("CC");
        let assignment = valence(&topology, "c07").unwrap();
        let rings = ring_info(&topology, "c07").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let topology_snapshot = topology.clone();
        let assignment_snapshot = assignment.clone();
        let rings_snapshot = rings.clone();
        let mut fresh_st = tpsa::DescriptorComputedState::default();
        let totals_h = crippen_totals(&input, true, false, &mut fresh_st).unwrap();
        let (logp, mr) = (totals_h.logp, totals_h.molar_refractivity);
        // Sequential left-to-right accumulation in ATOM ORDER (two C
        // rows then six H rows) — never a pre-multiplied 6*term.
        let expected_logp =
            0.0_f64 + 0.1441 + 0.1441 + 0.123 + 0.123 + 0.123 + 0.123 + 0.123 + 0.123;
        let expected_mr = 0.0_f64 + 2.503 + 2.503 + 1.057 + 1.057 + 1.057 + 1.057 + 1.057 + 1.057;
        assert_eq!(logp.to_bits(), expected_logp.to_bits());
        assert_eq!(mr.to_bits(), expected_mr.to_bits());
        // The ORIGINAL input blocks are untouched (working copy only).
        assert_eq!(topology, topology_snapshot);
        assert_eq!(assignment, assignment_snapshot);
        assert_eq!(rings, rings_snapshot);

        // Delegation discriminator: a molecule whose implicit H count is
        // ZERO (a lone Curium atom) — the AddHs transform is a no-op,
        // so the totals equal the C06 `crippen_totals` bits exactly.
        let atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::CM))];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(1);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 1, 0);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut fresh_st_h = tpsa::DescriptorComputedState::default();
        let totals_h = crippen_totals(&input, true, false, &mut fresh_st_h).unwrap();
        let (logp_h, mr_h) = (totals_h.logp, totals_h.molar_refractivity);
        let mut fresh_st = tpsa::DescriptorComputedState::default();
        let totals_c = crippen_totals(&input, false, false, &mut fresh_st).unwrap();
        let (logp, mr) = (totals_c.logp, totals_c.molar_refractivity);
        assert_eq!(logp_h.to_bits(), logp.to_bits());
        assert_eq!(mr_h.to_bits(), mr.to_bits());
        assert_eq!(logp_h.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn descriptor_c07_explicit_and_isotopic_hydrogens() {
        // Explicit-H case: a molecule that ALREADY has an explicit H
        // node — the AddHs owner adds only the missing implicit Hs.
        // Methane written with one explicit H ([H]C(H)(H)H) has no
        // implicit Hs left, so the totals equal the fully-added form.
        // Isotopic case: an isotope label on a heavy atom does NOT
        // shift its table row (the matcher sees the same element and
        // environment; H rows differ by H-node presence only).
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();

        let explicit = smiles_topology("[H]C([H])([H])[H]");
        let assignment_e = valence(&explicit, "c07e").unwrap();
        let rings_e = ring_info(&explicit, "c07e").unwrap();
        let input_e = DescriptorInput::new(
            &explicit,
            &coordinates,
            &properties,
            &assignment_e,
            &rings_e,
        );
        let mut fresh_st_e = tpsa::DescriptorComputedState::default();
        let totals_e = crippen_totals(&input_e, true, false, &mut fresh_st_e)
            .expect("explicit-H methane totals");
        let (logp_explicit, mr_explicit) = (totals_e.logp, totals_e.molar_refractivity);
        // No implicit Hs to add: the totals equal the no-transform form.
        let mut fresh_st_p = tpsa::DescriptorComputedState::default();
        let totals_p = crippen_totals(&input_e, false, false, &mut fresh_st_p).unwrap();
        let (logp_plain, mr_plain) = (totals_p.logp, totals_p.molar_refractivity);
        assert_eq!(logp_explicit.to_bits(), logp_plain.to_bits());
        assert_eq!(mr_explicit.to_bits(), mr_plain.to_bits());

        // Isotopic heavy atom: a 13C methane — the C row is UNCHANGED
        // ([CH4] matches the isotope-labeled carbon the same way), and
        // the added Hs type by the same H1 rows; totals equal the
        // plain methane-with-Hs totals bit-for-bit.
        let bond_spec = BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single);
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C).with_isotope(13)),
        ];
        let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
        let isotopic = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        // Real derived assignment so the carbons carry their true
        // implicit-H counts (a zero assignment would add no Hs).
        let assignment_i = valence(&isotopic, "c07i").unwrap();
        let rings_i = ring_info(&isotopic, "c07i").unwrap();
        let input_i = DescriptorInput::new(
            &isotopic,
            &coordinates,
            &properties,
            &assignment_i,
            &rings_i,
        );
        let mut fresh_st_i = tpsa::DescriptorComputedState::default();
        let totals_i = crippen_totals(&input_i, true, false, &mut fresh_st_i)
            .expect("isotopic C-C totals with Hs");
        let (logp_iso, mr_iso) = (totals_i.logp, totals_i.molar_refractivity);
        // Plain C-C ([CH3]C + H1 rows) reference bits for comparison,
        // sequential accumulation in atom order.
        let expected_logp =
            0.0_f64 + 0.1441 + 0.1441 + 0.123 + 0.123 + 0.123 + 0.123 + 0.123 + 0.123;
        let expected_mr = 0.0_f64 + 2.503 + 2.503 + 1.057 + 1.057 + 1.057 + 1.057 + 1.057 + 1.057;
        assert_eq!(logp_iso.to_bits(), expected_logp.to_bits());
        assert_eq!(mr_iso.to_bits(), expected_mr.to_bits());
    }

    #[test]
    fn descriptor_c08_shared_engine_and_call_order_matrix() {
        // Frozen C08 semantics (Step-428 audit): clogp/mr are pure
        // projections of the ONE shared totals engine (header-default
        // includeHs=true/force=false); exact H-mode bits from the frozen
        // source table. The 4-arm matrix uses fresh states on the same
        // input; the separate shared-state pair checks warm projections.
        // The force arm exercises the
        // engine directly (the projections take no options, matching
        // the source) with a fresh state.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("CCO");
        let assignment = valence(&topology, "c08").unwrap();
        let rings = ring_info(&topology, "c08").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());

        // Reference: the shared canonical engine at the HEADER DEFAULTS
        // (includeHs=true, force=false — Crippen.h:72-75).
        let mut ref_state = tpsa::DescriptorComputedState::default();
        let ref_totals = crippen_totals(&input, true, false, &mut ref_state).unwrap();
        let ref_logp = ref_totals.logp.to_bits();
        let ref_mr = ref_totals.molar_refractivity.to_bits();

        // The 4-arm call-order matrix: each projection called twice in
        // interleaved order — every call returns the shared-engine bits.
        let mut calls = 0u32;
        for (first_clogp, first_mr) in [(true, true), (true, false), (false, true), (false, false)]
        {
            if first_clogp {
                let mut st = tpsa::DescriptorComputedState::default();
                assert_eq!(crippen_clogp(&input, &mut st).unwrap().to_bits(), ref_logp);
                calls += 1;
            }
            if first_mr {
                let mut st = tpsa::DescriptorComputedState::default();
                assert_eq!(crippen_mr(&input, &mut st).unwrap().to_bits(), ref_mr);
                calls += 1;
            }
            // The trailing halves in reverse order.
            let mut st = tpsa::DescriptorComputedState::default();
            assert_eq!(crippen_mr(&input, &mut st).unwrap().to_bits(), ref_mr);
            let mut st = tpsa::DescriptorComputedState::default();
            assert_eq!(crippen_clogp(&input, &mut st).unwrap().to_bits(), ref_logp);
            calls += 2;
        }
        assert_eq!(calls, 12);

        // These two fresh-state engine calls exercise each force setting;
        // they are cold computations, not warm-cache evidence. The
        // following projections share the already-populated reference
        // state and should return its scalar pair without preparation.
        for force in [false, true] {
            let mut st = tpsa::DescriptorComputedState::default();
            let totals = crippen_totals(&input, true, force, &mut st).unwrap();
            assert_eq!(totals.logp.to_bits(), ref_logp, "force={force}");
            assert_eq!(totals.molar_refractivity.to_bits(), ref_mr, "force={force}");
        }
        let mut shared = ref_state.clone();
        assert_eq!(
            crippen_clogp(&input, &mut shared).unwrap().to_bits(),
            ref_logp
        );
        assert_eq!(crippen_mr(&input, &mut shared).unwrap().to_bits(), ref_mr);
        // One-shared-engine census: the 12 cold H-mode projections +
        // 2 force arms + 1 reference each legitimately prepare ONE
        // final valence (the frozen hydrogen contract allows exactly
        // one per cold H computation); the SHARED warm projections
        // after them add none.
        assert_eq!(
            VALENCE_HELPER_ENTRIES.with(|c| c.get()),
            helper_before + 15,
            "one final valence per cold H call; none on warm serves"
        );

        // Frozen literal anchors (methane/methanol/Cm from the C06
        // tables) prove the halves independently, not just equality.
        let topology = smiles_topology("C");
        let assignment = valence(&topology, "c08b").unwrap();
        let rings = ring_info(&topology, "c08b").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        // Frozen H-bits reference: methane H-mode defaults.
        let mut st = tpsa::DescriptorComputedState::default();
        assert_eq!(
            crippen_clogp(&input, &mut st).unwrap().to_bits(),
            0x3fe4_5aee_631f_8a09_u64
        );
        let mut st = tpsa::DescriptorComputedState::default();
        assert_eq!(
            crippen_mr(&input, &mut st).unwrap().to_bits(),
            0x401a_ec8b_4395_8106_u64
        );
        let topology = smiles_topology("CO");
        let assignment = valence(&topology, "c08c").unwrap();
        let rings = ring_info(&topology, "c08c").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        // Frozen H-bits reference: methanol H-mode defaults.
        let mut st = tpsa::DescriptorComputedState::default();
        assert_eq!(
            crippen_clogp(&input, &mut st).unwrap().to_bits(),
            0xbfd9_0e56_0418_9375_u64
        );
        let mut st = tpsa::DescriptorComputedState::default();
        assert_eq!(
            crippen_mr(&input, &mut st).unwrap().to_bits(),
            0x4020_491d_14e3_bcd3_u64
        );
        let atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::CM))];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(1);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 1, 0);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut st = tpsa::DescriptorComputedState::default();
        assert_eq!(
            crippen_clogp(&input, &mut st).unwrap().to_bits(),
            0.0_f64.to_bits()
        );
        let mut st = tpsa::DescriptorComputedState::default();
        assert_eq!(
            crippen_mr(&input, &mut st).unwrap().to_bits(),
            0.0_f64.to_bits()
        );
    }

    #[test]
    fn descriptor_v01_upper_bound_and_accumulation() {
        // Frozen V01 semantics (Step-436 audit): upper_bound placement
        // (equality -> the NEXT bin; below-all -> bin 0; above-all ->
        // the final slot; strictly-between), source accumulation order
        // (left-to-right += observable in the final bits), and the
        // empty-bin/empty-contribs behavior. Literal frozen bits.
        let bins = [0.1_f64, 0.2, 0.3];

        // Equality: bVal == bins[1] (0.2) -> upper_bound returns 2
        // (the FIRST index greater than 0.2), NOT 1.
        let res = assign_contribs_to_bins(&[5.0], &[0.2], &bins).unwrap();
        assert_eq!(res.len(), 4);
        assert_eq!(res, vec![0.0, 0.0, 5.0, 0.0]);

        // Below all: bVal = 0.05 -> bin 0.
        let res = assign_contribs_to_bins(&[7.0], &[0.05], &bins).unwrap();
        assert_eq!(res, vec![7.0, 0.0, 0.0, 0.0]);

        // Above all: bVal = 0.5 -> the final slot (index 3).
        let res = assign_contribs_to_bins(&[9.0], &[0.5], &bins).unwrap();
        assert_eq!(res, vec![0.0, 0.0, 0.0, 9.0]);

        // Strictly between edges 0.1 and 0.2 -> bin 1.
        let res = assign_contribs_to_bins(&[3.0], &[0.15], &bins).unwrap();
        assert_eq!(res, vec![0.0, 3.0, 0.0, 0.0]);

        // Source accumulation ORDER: three atoms hit ONE bin — the
        // final bits equal the left-to-right += sum (0.1+0.2+0.3 in
        // that order), never a reordered or pre-summed total.
        let res = assign_contribs_to_bins(&[0.1, 0.2, 0.3], &[0.05, 0.05, 0.05], &bins).unwrap();
        let expected = 0.0_f64 + 0.1 + 0.2 + 0.3;
        assert_eq!(res[0].to_bits(), expected.to_bits());
        assert_eq!(res[1..], vec![0.0, 0.0, 0.0]);

        // Empty bins: upper_bound on empty = 0 -> every atom in res[0],
        // res length 1.
        let res = assign_contribs_to_bins(&[1.0, 2.0], &[0.5, 0.5], &[]).unwrap();
        assert_eq!(res.len(), 1);
        let expected = 0.0_f64 + 1.0 + 2.0;
        assert_eq!(res[0].to_bits(), expected.to_bits());

        // Empty contribs: res of length bins.len()+1, all +0.0 bits.
        let res = assign_contribs_to_bins(&[], &[], &bins).unwrap();
        assert_eq!(res.len(), 4);
        assert!(res.iter().all(|v| v.to_bits() == 0.0_f64.to_bits()));
    }

    #[test]
    fn descriptor_v01_unsorted_bins_and_preconditions() {
        // Frozen V01 semantics (Step-436 audit): the UNSORTED custom
        // bins cell freezes the partition-point result literally — the
        // direct upper_bound analog on a non-partitioned range (both
        // compute the first index where the binary-search walk's
        // predicate flips); NO sorting and NO sortedness validation
        // exist in the source or the owner. The PRECONDITION
        // violation returns the typed MismatchedBinArrays error.
        // Unsorted bins [0.3, 0.1, 0.2], v = 0.15: the walk finds
        // mid=1 (0.1 <= 0.15 true -> lo=2), mid=2 (0.2 <= 0.15
        // false -> hi=2) -> index 2.
        let unsorted = [0.3_f64, 0.1, 0.2];
        let res = assign_contribs_to_bins(&[1.0], &[0.15], &unsorted).unwrap();
        assert_eq!(res.len(), 4);
        assert_eq!(res, vec![0.0, 0.0, 1.0, 0.0]);

        // PRECONDITION: contribs.len() != bin_prop.len() -> the typed
        // error, never a panic.
        let err = assign_contribs_to_bins(&[1.0, 2.0], &[1.0], &[0.5]).unwrap_err();
        assert_eq!(
            err,
            DescriptorError::MismatchedBinArrays {
                contribs_len: 2,
                bin_prop_len: 1,
                bins_len: 1
            }
        );
    }

    #[test]
    fn descriptor_v02_default_bins_and_owner_delegation() {
        // Frozen V02 semantics (Step-444 audit): default-None selects
        // the FROZEN 11-boundary blist -> EXACTLY 12 outputs; the VSA
        // amounts come from the Labute owner with
        // include_hydrogens=true (source `tmp` = the discarded
        // hydrogen lump); the binning property is the CRIPPEN logP row
        // per atom; placement follows V01 upper_bound semantics.
        // Hand-derived bits from the pinned radii (C 0.77, O 0.66,
        // H 0.33) and the frozen Crippen table literals
        // ([CH4] 0.1441, C3 -0.2035, O2 -0.2893).

        // The frozen default boundary list itself (source blist[11]).
        assert_eq!(
            vsa::DEFAULT_SLOGP_BINS,
            [-0.4, -0.2, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6]
        );

        // Methane: one C (r=0.77), no bonds, H pass with rj=0.33.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("C");
        let assignment = valence(&topology, "v02a").unwrap();
        let rings = ring_info(&topology, "v02a").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut state = tpsa::DescriptorComputedState::default();
        let res = slogp_vsa(&input, None, false, &mut state).unwrap();
        assert_eq!(res.len(), 12);
        let dij_c_h = (0.77_f64 - 0.33).abs().max(0.77 + 0.33).min(0.77 + 0.33);
        let vi_c = 0.0_f64 + (0.33 * 0.33 - (0.77 - dij_c_h) * (0.77 - dij_c_h) / dij_c_h);
        let vsa_c = std::f64::consts::PI * 0.77 * (4.0 * 0.77 - vi_c);
        // logp row 0.1441: first boundary greater than it is 0.15
        // (index 4) — everything else stays +0.0.
        assert_eq!(res[4].to_bits(), vsa_c.to_bits());
        for (i, v) in res.iter().enumerate() {
            if i != 4 {
                assert_eq!(v.to_bits(), 0.0_f64.to_bits(), "bin {i}");
            }
        }

        // Same-owner discriminator: a manual composition over the same
        // two owners plus the V01 binning equals the wrapper output.
        let lab = labute_contributions(&input, true, false, &mut state).unwrap();
        let mut st_c = tpsa::DescriptorComputedState::default();
        let crip = crippen_contributions(&input, true, &mut st_c, None, None).unwrap();
        let manual = assign_contribs_to_bins(&lab.atoms, &crip.logp, &vsa::DEFAULT_SLOGP_BINS);
        assert_eq!(manual.unwrap(), res);

        // Methanol: C (0.77) - O (0.66) single bond, both logp rows
        // (-0.2035, -0.2893) land in bin 1 (first boundary greater
        // than each is -0.2); accumulation order is atom order
        // (C then O) — observable in the bits.
        let topology = smiles_topology("CO");
        let assignment = valence(&topology, "v02b").unwrap();
        let rings = ring_info(&topology, "v02b").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut state = tpsa::DescriptorComputedState::default();
        let res = slogp_vsa(&input, None, false, &mut state).unwrap();
        assert_eq!(res.len(), 12);
        let dij_bond = (0.77_f64 - 0.66).abs().max(0.77 + 0.66).min(0.77 + 0.66);
        let dij_o_h = (0.66_f64 - 0.33).abs().max(0.66 + 0.33).min(0.66 + 0.33);
        let mut vi_c = 0.0_f64;
        vi_c += 0.66 * 0.66 - (0.77 - dij_bond) * (0.77 - dij_bond) / dij_bond;
        vi_c += 0.33 * 0.33 - (0.77 - dij_c_h) * (0.77 - dij_c_h) / dij_c_h;
        let mut vi_o = 0.0_f64;
        vi_o += 0.77 * 0.77 - (0.66 - dij_bond) * (0.66 - dij_bond) / dij_bond;
        vi_o += 0.33 * 0.33 - (0.66 - dij_o_h) * (0.66 - dij_o_h) / dij_o_h;
        let vsa_c = std::f64::consts::PI * 0.77 * (4.0 * 0.77 - vi_c);
        let vsa_o = std::f64::consts::PI * 0.66 * (4.0 * 0.66 - vi_o);
        let expected_bin1 = vsa_c + vsa_o;
        assert_eq!(res[1].to_bits(), expected_bin1.to_bits());
        for (i, v) in res.iter().enumerate() {
            if i != 1 {
                assert_eq!(v.to_bits(), 0.0_f64.to_bits(), "bin {i}");
            }
        }
    }

    #[test]
    fn descriptor_v02_custom_bins_and_edges() {
        // Frozen V02 semantics (Step-444 audit): custom Some(bins)
        // copies the caller's boundaries (m -> m+1 outputs, empty ->
        // 1); a logp row exactly ON an edge lands in the bin AFTER the
        // edge (V01 upper_bound equality semantics).
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("C");
        let assignment = valence(&topology, "v02c").unwrap();
        let rings = ring_info(&topology, "v02c").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut state = tpsa::DescriptorComputedState::default();

        let dij_c_h = (0.77_f64 - 0.33).abs().max(0.77 + 0.33).min(0.77 + 0.33);
        let vi_c = 0.0_f64 + (0.33 * 0.33 - (0.77 - dij_c_h) * (0.77 - dij_c_h) / dij_c_h);
        let vsa_c = std::f64::consts::PI * 0.77 * (4.0 * 0.77 - vi_c);

        // Edge equality: logp row 0.1441 equals the LAST boundary ->
        // the final slot (index 2 of 3).
        let res = slogp_vsa(&input, Some(&[-0.4, 0.1441]), false, &mut state).unwrap();
        assert_eq!(res.len(), 3);
        assert_eq!(res[2].to_bits(), vsa_c.to_bits());
        assert_eq!(res[0].to_bits(), 0.0_f64.to_bits());
        assert_eq!(res[1].to_bits(), 0.0_f64.to_bits());

        // Edge equality on the FIRST boundary -> index 1 (after the
        // edge, not 0).
        let res = slogp_vsa(&input, Some(&[0.1441, 0.5]), false, &mut state).unwrap();
        assert_eq!(res.len(), 3);
        assert_eq!(res[1].to_bits(), vsa_c.to_bits());

        // One boundary below the row -> the final slot.
        let res = slogp_vsa(&input, Some(&[0.0]), false, &mut state).unwrap();
        assert_eq!(res.len(), 2);
        assert_eq!(res[1].to_bits(), vsa_c.to_bits());

        // Empty custom bins -> one output accumulating everything.
        let res = slogp_vsa(&input, Some(&[]), false, &mut state).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].to_bits(), vsa_c.to_bits());

        // None vs empty-Some discriminator: 12 vs 1 outputs.
        let res_none = slogp_vsa(&input, None, false, &mut state).unwrap();
        assert_eq!(res_none.len(), 12);
        assert_eq!(res_none[4].to_bits(), vsa_c.to_bits());
    }

    #[test]
    fn descriptor_v03_default_bins_and_owner_delegation() {
        // Frozen V03 semantics (Step-452 audit): default-None selects
        // the FROZEN 9-boundary blist -> EXACTLY 10 outputs; the
        // binning property is the CRIPPEN MR row per atom (the logP
        // rows are the DISCARDED vector — the mirror of V02); same
        // Labute owner with include_hydrogens=true. Hand-derived bits
        // from the pinned radii (C 0.77, O 0.66, H 0.33) and the
        // frozen Crippen table literals ([CH4] 2.503, C3 2.753, O2
        // 0.8238).

        // The frozen default boundary list itself (source blist[9]).
        assert_eq!(
            vsa::DEFAULT_SMR_BINS,
            [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0]
        );

        // Methane: MR row 2.503 -> strictly between 2.45 and 2.75 ->
        // bin 4 (first boundary greater is 2.75); everything else
        // stays +0.0.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("C");
        let assignment = valence(&topology, "v03a").unwrap();
        let rings = ring_info(&topology, "v03a").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut state = tpsa::DescriptorComputedState::default();
        let res = smr_vsa(&input, None, false, &mut state).unwrap();
        assert_eq!(res.len(), 10);
        let dij_c_h = (0.77_f64 - 0.33).abs().max(0.77 + 0.33).min(0.77 + 0.33);
        let vi_c = 0.0_f64 + (0.33 * 0.33 - (0.77 - dij_c_h) * (0.77 - dij_c_h) / dij_c_h);
        let vsa_c = std::f64::consts::PI * 0.77 * (4.0 * 0.77 - vi_c);
        assert_eq!(res[4].to_bits(), vsa_c.to_bits());
        for (i, v) in res.iter().enumerate() {
            if i != 4 {
                assert_eq!(v.to_bits(), 0.0_f64.to_bits(), "bin {i}");
            }
        }

        // Same-owner discriminator keyed on .mr: a manual composition
        // over the same two owners plus the V01 binning equals the
        // wrapper output.
        let lab = labute_contributions(&input, true, false, &mut state).unwrap();
        let mut st_c = tpsa::DescriptorComputedState::default();
        let crip = crippen_contributions(&input, true, &mut st_c, None, None).unwrap();
        let manual =
            assign_contribs_to_bins(&lab.atoms, &crip.molar_refractivity, &vsa::DEFAULT_SMR_BINS);
        assert_eq!(manual.unwrap(), res);

        // Methanol: MR rows 2.753 (C3) and 0.8238 (O2) land in TWO
        // DIFFERENT bins (5 and 0 — unlike V02's shared bin 1):
        // res[5] = VSA_C bits, res[0] = VSA_O bits.
        let topology = smiles_topology("CO");
        let assignment = valence(&topology, "v03b").unwrap();
        let rings = ring_info(&topology, "v03b").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut state = tpsa::DescriptorComputedState::default();
        let res = smr_vsa(&input, None, false, &mut state).unwrap();
        assert_eq!(res.len(), 10);
        let dij_bond = (0.77_f64 - 0.66).abs().max(0.77 + 0.66).min(0.77 + 0.66);
        let dij_o_h = (0.66_f64 - 0.33).abs().max(0.66 + 0.33).min(0.66 + 0.33);
        let mut vi_c = 0.0_f64;
        vi_c += 0.66 * 0.66 - (0.77 - dij_bond) * (0.77 - dij_bond) / dij_bond;
        vi_c += 0.33 * 0.33 - (0.77 - dij_c_h) * (0.77 - dij_c_h) / dij_c_h;
        let mut vi_o = 0.0_f64;
        vi_o += 0.77 * 0.77 - (0.66 - dij_bond) * (0.66 - dij_bond) / dij_bond;
        vi_o += 0.33 * 0.33 - (0.66 - dij_o_h) * (0.66 - dij_o_h) / dij_o_h;
        let vsa_c = std::f64::consts::PI * 0.77 * (4.0 * 0.77 - vi_c);
        let vsa_o = std::f64::consts::PI * 0.66 * (4.0 * 0.66 - vi_o);
        assert_eq!(res[5].to_bits(), vsa_c.to_bits());
        assert_eq!(res[0].to_bits(), vsa_o.to_bits());
        for (i, v) in res.iter().enumerate() {
            if i != 0 && i != 5 {
                assert_eq!(v.to_bits(), 0.0_f64.to_bits(), "bin {i}");
            }
        }
    }

    #[test]
    fn descriptor_v03_custom_bins_and_edges() {
        // Frozen V03 semantics (Step-452 audit): custom Some(bins)
        // copies the caller's boundaries (m -> m+1 outputs, empty ->
        // 1); an MR row exactly ON an edge lands in the bin AFTER the
        // edge; a boundary ABOVE the row lands in bin 0.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("C");
        let assignment = valence(&topology, "v03c").unwrap();
        let rings = ring_info(&topology, "v03c").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut state = tpsa::DescriptorComputedState::default();

        let dij_c_h = (0.77_f64 - 0.33).abs().max(0.77 + 0.33).min(0.77 + 0.33);
        let vi_c = 0.0_f64 + (0.33 * 0.33 - (0.77 - dij_c_h) * (0.77 - dij_c_h) / dij_c_h);
        let vsa_c = std::f64::consts::PI * 0.77 * (4.0 * 0.77 - vi_c);

        // Edge equality: MR row 2.503 equals the LAST boundary -> the
        // final slot (index 2 of 3).
        let res = smr_vsa(&input, Some(&[2.0, 2.503]), false, &mut state).unwrap();
        assert_eq!(res.len(), 3);
        assert_eq!(res[2].to_bits(), vsa_c.to_bits());
        assert_eq!(res[0].to_bits(), 0.0_f64.to_bits());
        assert_eq!(res[1].to_bits(), 0.0_f64.to_bits());

        // Edge equality on the FIRST boundary -> index 1 (after the
        // edge, not 0).
        let res = smr_vsa(&input, Some(&[2.503, 3.0]), false, &mut state).unwrap();
        assert_eq!(res.len(), 3);
        assert_eq!(res[1].to_bits(), vsa_c.to_bits());

        // Boundary ABOVE the row -> bin 0 (the below-all slot).
        let res = smr_vsa(&input, Some(&[4.0]), false, &mut state).unwrap();
        assert_eq!(res.len(), 2);
        assert_eq!(res[0].to_bits(), vsa_c.to_bits());

        // Empty custom bins -> one output accumulating everything.
        let res = smr_vsa(&input, Some(&[]), false, &mut state).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].to_bits(), vsa_c.to_bits());

        // None vs empty-Some discriminator: 10 vs 1 outputs.
        let res_none = smr_vsa(&input, None, false, &mut state).unwrap();
        assert_eq!(res_none.len(), 10);
        assert_eq!(res_none[4].to_bits(), vsa_c.to_bits());
    }

    #[test]
    fn descriptor_crippen_fix_numeric_success_and_error_allocation() {
        // Frozen CRIPPEN-CLOSE6 semantics: pin BOTH outcomes of the
        // numeric wrapper after the lazy-error fix — every success
        // cell returns its exact pinned-Boost bits with NO error
        // constructed (laziness is code-inspected; the test pins that
        // successes never produce the typed error), and every failing
        // cell returns CrippenParamNumeric preserving the ORIGINAL
        // cell. All existing descriptor_c02_numeric_source_
        // assertions remain untouched alongside.

        // Success cells: empty -> +0.0.
        assert_eq!(
            crate::crippen::parse_crippen_logp("").unwrap().to_bits(),
            0.0_f64.to_bits()
        );

        // Ordinary decimal: the table's own [CH4] logP cell.
        assert_eq!(
            crate::crippen::parse_crippen_logp("0.1441")
                .unwrap()
                .to_bits(),
            0.1441_f64.to_bits()
        );

        // Signed zero: "-0.0" keeps the negative zero bits.
        assert_eq!(
            crate::crippen::parse_crippen_logp("-0.0")
                .unwrap()
                .to_bits(),
            (-0.0_f64).to_bits()
        );

        // inf family: exact-case-insensitive tokens, signs preserved.
        assert_eq!(
            crate::crippen::parse_crippen_logp("inf").unwrap().to_bits(),
            f64::INFINITY.to_bits()
        );
        assert_eq!(
            crate::crippen::parse_crippen_logp("-INF")
                .unwrap()
                .to_bits(),
            f64::NEG_INFINITY.to_bits()
        );
        assert_eq!(
            crate::crippen::parse_crippen_logp("Infinity")
                .unwrap()
                .to_bits(),
            f64::INFINITY.to_bits()
        );

        // nan family: quiet-NaN bits; payload ignored; sign flips the
        // sign bit only.
        assert_eq!(
            crate::crippen::parse_crippen_logp("nan").unwrap().to_bits(),
            0x7ff8_0000_0000_0000_u64
        );
        assert_eq!(
            crate::crippen::parse_crippen_logp("nan(payload)")
                .unwrap()
                .to_bits(),
            0x7ff8_0000_0000_0000_u64
        );
        assert_eq!(
            crate::crippen::parse_crippen_logp("-NAN(x9)")
                .unwrap()
                .to_bits(),
            0xfff8_0000_0000_0000_u64
        );

        // Error cells: the ORIGINAL cell is preserved verbatim in the
        // typed error — malformed text and decimal overflow alike.
        let err = crate::crippen::parse_crippen_logp("abc").unwrap_err();
        assert_eq!(
            err,
            DescriptorError::CrippenParamNumeric {
                field: "logp",
                cell: "abc".to_string(),
            }
        );
        let err = crate::crippen::parse_crippen_logp("1e309").unwrap_err();
        assert_eq!(
            err,
            DescriptorError::CrippenParamNumeric {
                field: "logp",
                cell: "1e309".to_string(),
            }
        );

        // The MR-side source fallback: malformed cells store 0.0.
        assert_eq!(
            crate::crippen::parse_crippen_mr("abc").to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            crate::crippen::parse_crippen_mr("1e309").to_bits(),
            0.0_f64.to_bits()
        );
    }

    #[test]
    fn descriptor_crippen_fix_state_sixteen_presence_clone_clear() {
        // Frozen CRIPPEN-CLOSE12 semantics: all 2^4 Crippen presence
        // patterns; every clone preserves bits/order and INDEPENDENT
        // presence per field; clear() empties all four Crippen fields
        // AND the TPSA/Labute slots; the Hydrogens variant exposes the
        // real core error through a BORROWED Error::source downcast;
        // unrelated descriptor slots are preserved by Crippen state
        // operations.
        let seed_rows_lp = [11.5_f64, 12.5];
        let seed_rows_mr = [21.5_f64, 22.5];
        let seed_lp = 789.0_f64;
        let seed_mr = 456.0_f64;
        let bits = |v: &[f64]| -> Vec<u64> { v.iter().map(|x| x.to_bits()).collect() };

        let mut patterns = 0u32;
        for lp_rows in [false, true] {
            for mr_rows in [false, true] {
                for lp in [false, true] {
                    for mr in [false, true] {
                        patterns += 1;
                        let mut state = tpsa::DescriptorComputedState::default();
                        {
                            let slot = state.crippen_slot_mut();
                            slot.logp_rows = lp_rows.then(|| seed_rows_lp.to_vec());
                            slot.mr_rows = mr_rows.then(|| seed_rows_mr.to_vec());
                            slot.logp = lp.then_some(seed_lp);
                            slot.mr = mr.then_some(seed_mr);
                        }
                        // Presence discrimination: each field's
                        // presence follows ONLY its own pattern bit.
                        let slot = state.crippen_slot();
                        assert_eq!(slot.logp_rows.is_some(), lp_rows);
                        assert_eq!(slot.mr_rows.is_some(), mr_rows);
                        assert_eq!(slot.logp.is_some(), lp);
                        assert_eq!(slot.mr.is_some(), mr);
                        // Clone preserves bits/order/presence exactly.
                        let clone = state.clone();
                        let c = clone.crippen_slot();
                        assert_eq!(
                            c.logp_rows.as_ref().map(|v| bits(v)).as_deref(),
                            lp_rows.then(|| bits(&seed_rows_lp)).as_deref()
                        );
                        assert_eq!(
                            c.mr_rows.as_ref().map(|v| bits(v)).as_deref(),
                            mr_rows.then(|| bits(&seed_rows_mr)).as_deref()
                        );
                        assert_eq!(c.logp.map(f64::to_bits), lp.then(|| seed_lp.to_bits()));
                        assert_eq!(c.mr.map(f64::to_bits), mr.then(|| seed_mr.to_bits()));
                        // Mutating the clone leaves the original
                        // untouched (no shared mutation authority).
                        let mut mutated = clone;
                        mutated.crippen_slot_mut().logp = Some(-1.0);
                        assert_eq!(state.crippen_slot().logp, lp.then_some(seed_lp));
                    }
                }
            }
        }
        assert_eq!(patterns, 16);

        // Unrelated slots survive Crippen seeding and cloning by bits.
        let mut state = tpsa::DescriptorComputedState::default();
        state.slot_mut_for_tests(false).scalar = Some(123.0);
        state.slot_mut_for_tests(false).contributions = Some(vec![1.0, 2.0]);
        state.slot_mut_for_tests(true).scalar = Some(456.0);
        state.slot_mut_for_tests(true).contributions = Some(vec![3.0, 4.0]);
        {
            let slot = state.labute_slot_mut();
            slot.rows = Some(vec![4.0, 5.0]);
            slot.hydrogens = Some(7.0);
            slot.asa = Some(789.0);
        }
        state.crippen_slot_mut().logp = Some(seed_lp);
        let clone = state.clone();
        assert_eq!(
            clone.slot_for_tests(false).scalar.map(f64::to_bits),
            Some(123.0_f64.to_bits())
        );
        assert_eq!(
            clone
                .slot_for_tests(true)
                .contributions
                .as_ref()
                .map(|v| bits(v))
                .as_deref(),
            Some(&bits(&[3.0, 4.0])[..])
        );
        assert_eq!(
            clone.labute_slot().asa.map(f64::to_bits),
            Some(789.0_f64.to_bits())
        );
        assert_eq!(
            clone.crippen_slot().logp.map(f64::to_bits),
            Some(seed_lp.to_bits())
        );

        // Explicit clear empties BOTH TPSA slots, all THREE Labute
        // presences AND all FOUR Crippen presences.
        state.clear();
        let empty = tpsa::DescriptorComputedState::default();
        assert_eq!(state, empty);

        // The Hydrogens variant preserves the real core error through
        // a BORROWED Error::source downcast — never text/Unsupported.
        let err = DescriptorError::Hydrogens {
            function: "crippen_totals",
            source: cosmolkit_core::HydrogenError::Unsupported {
                operation: "add_hydrogens_with_params",
                reason: "state fixture",
            },
        };
        use std::error::Error;
        let down = err
            .source()
            .and_then(|s| s.downcast_ref::<cosmolkit_core::HydrogenError>())
            .expect("borrowed HydrogenError source");
        assert_eq!(
            *down,
            cosmolkit_core::HydrogenError::Unsupported {
                operation: "add_hydrogens_with_params",
                reason: "state fixture",
            }
        );
    }

    #[test]
    fn descriptor_crippen_fix_kernel_sinks_packing_counters() {
        // Frozen CRIPPEN-CLOSE18 semantics: the ONE cold kernel with
        // all four optional-sink combinations on CC (numeric rows
        // identical across combos; typed sinks get the table values);
        // untyped Cm keeps CALLER sentinels (no invented table types);
        // the packed 64/65 boundary (word 1 addresses atom 64); and
        // exactly one kernel + one context per cold invocation.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();

        // CC: both atoms typed by the FIRST matching row (C1, idx 1 —
        // the source first-match rule; supervisor reference bits for
        // the no-H totals: logP 0x3fd271de69ad42c4, MR
        // 0x40140624dd2f1aa0).
        let topology = smiles_topology("CC");
        let assignment = valence(&topology, "cc18a").unwrap();
        let rings = ring_info(&topology, "cc18a").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);

        let kernel_before = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
        let context_before = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
        let mut lp = vec![7.0; 2];
        let mut mr = vec![8.0; 2];
        crate::crippen::crippen_cold_kernel(&input, &mut lp, &mut mr, None, None).unwrap();
        assert_eq!(
            crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - kernel_before,
            1
        );
        assert_eq!(
            crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - context_before,
            1
        );
        // Numeric rows were ZERO-INITIALIZED by the kernel, then
        // first-match typed; the left-to-right no-H totals reproduce
        // the supervisor's independent reference bits.
        let x = lp[0];
        let y = mr[0];
        assert_eq!(lp, vec![x, x]);
        assert_eq!(mr, vec![y, y]);
        assert_eq!((0.0_f64 + x + x).to_bits(), 0x3fd2_71de_69ad_42c4_u64);
        assert_eq!((0.0_f64 + y + y).to_bits(), 0x4014_0624_dd2f_1aa0_u64);

        // The four sink combinations: numeric rows identical; typed
        // sinks overwritten; None sinks simply not supplied.
        for (with_types, with_labels) in
            [(false, false), (true, false), (false, true), (true, true)]
        {
            let mut lp2 = vec![0.0; 2];
            let mut mr2 = vec![0.0; 2];
            let mut types = vec![91u32, 92];
            let mut labels = vec!["sentinelA".to_string(), "sentinelB".to_string()];
            crate::crippen::crippen_cold_kernel(
                &input,
                &mut lp2,
                &mut mr2,
                with_types.then_some(&mut types),
                with_labels.then_some(&mut labels),
            )
            .unwrap();
            assert_eq!(lp2, lp, "numeric rows must not depend on sinks");
            assert_eq!(mr2, mr);
            if with_types {
                assert_eq!(types, vec![1u32, 1]);
            }
            if with_labels {
                assert_eq!(labels, vec!["C1".to_string(), "C1".to_string()]);
            }
        }

        // Untyped Cm: caller sentinels survive — no invented table
        // types or labels; numeric rows stay the kernel's zeros.
        let atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::CM))];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let assignment = t01_zero_assignment(1);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 1, 0);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut lp = vec![7.0; 1];
        let mut mr = vec![8.0; 1];
        let mut types = vec![91u32];
        let mut labels = vec!["sentinelB".to_string()];
        crate::crippen::crippen_cold_kernel(
            &input,
            &mut lp,
            &mut mr,
            Some(&mut types),
            Some(&mut labels),
        )
        .unwrap();
        assert_eq!(lp, vec![0.0]);
        assert_eq!(mr, vec![0.0]);
        assert_eq!(types, vec![91u32]);
        assert_eq!(labels, vec!["sentinelB".to_string()]);

        // Packed layout contract: 64 atoms -> ONE word, 65 -> TWO.
        assert_eq!(64_usize.div_ceil(64), 1);
        assert_eq!(65_usize.div_ceil(64), 2);

        // The 65-carbon chain: atom index 64 lives in word 1 bit 0 —
        // it must be typed exactly like atom 0 (both terminal CH3).
        let topology = smiles_topology(&"C".repeat(65));
        let assignment = valence(&topology, "cc18b").unwrap();
        let rings = ring_info(&topology, "cc18b").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut lp = vec![0.0; 65];
        let mut mr = vec![0.0; 65];
        let mut types = vec![777u32; 65];
        crate::crippen::crippen_cold_kernel(&input, &mut lp, &mut mr, Some(&mut types), None)
            .unwrap();
        assert_ne!(types[0], 777);
        assert_eq!(types[64], types[0]);
        assert_eq!(lp[64].to_bits(), lp[0].to_bits());
        assert_eq!(mr[64].to_bits(), mr[0].to_bits());
        // And a 64-carbon chain stays inside the single word.
        let topology = smiles_topology(&"C".repeat(64));
        let assignment = valence(&topology, "cc18c").unwrap();
        let rings = ring_info(&topology, "cc18c").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut lp = vec![0.0; 64];
        let mut mr = vec![0.0; 64];
        let mut types = vec![777u32; 64];
        crate::crippen::crippen_cold_kernel(&input, &mut lp, &mut mr, Some(&mut types), None)
            .unwrap();
        assert_eq!(types[63], types[0]);
        assert_eq!(lp[63].to_bits(), lp[0].to_bits());
    }

    #[test]
    fn descriptor_crippen_fix_kernel_complete_product() {
        // Frozen CRIPPEN-KERNEL-PROOF4: exactly 20 REAL cold calls —
        // five shapes {empty, CC, Cm, 64-C-chain, 65-C-chain} x four
        // optional-sink combinations {none, types, labels, both}. The
        // literal table comes from Crippen.cpp's cold loop and default
        // rows 0 ([CH4]), 1 ([CH3]C), 2 ([CH2](C)C) — all C1 with
        // 0.1441/2.503 — never from a tested invocation.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let bits = |v: &[f64]| -> Vec<u64> { v.iter().map(|x| x.to_bits()).collect() };

        // Build the five shapes with source-shaped prepared inputs.
        let mut shapes: Vec<(&str, TopologyBlock, usize, Option<u64>, usize, Vec<u32>)> =
            Vec::new();
        let t = smiles_topology("");
        shapes.push(("empty", t, 0, None, 1, vec![]));
        let t = smiles_topology("CC");
        shapes.push(("CC", t, 2, Some(3), 2, vec![1, 1]));
        let t = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms: vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::CM))],
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        shapes.push(("Cm", t, 1, Some(1), 110, vec![91]));
        let t = smiles_topology(&"C".repeat(64));
        let mut chain64 = vec![1u32];
        chain64.extend(std::iter::repeat_n(2, 62));
        chain64.push(1);
        shapes.push(("64-C-chain", t, 64, Some(u64::MAX), 3, chain64));
        let t = smiles_topology(&"C".repeat(65));
        let mut chain65 = vec![1u32];
        chain65.extend(std::iter::repeat_n(2, 63));
        chain65.push(1);
        shapes.push(("65-C-chain", t, 65, Some(1), 3, chain65));

        let mut real_calls = 0u32;
        let expected_packed_words = [0usize, 1, 1, 1, 2];
        for (shape_index, (name, topology, atoms, last_word, row_visits, expected_types)) in
            shapes.into_iter().enumerate()
        {
            // Atom count asserted BEFORE any call.
            assert_eq!(topology.atoms.len(), atoms, "{name} atom count");
            let assignment = valence(&topology, "kp4").unwrap();
            let rings = ring_info(&topology, "kp4").unwrap();
            // Baseline AFTER this shape's fixture preparation: any
            // later delta would be the KERNEL re-preparing (it never
            // does — the counter is a cfg(test) instrumentation of
            // the shared prepared-valence helper).
            let valence_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
            let input =
                DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
            // Structural input snapshot (immutably borrowed, but the
            // content is re-checked after every call).
            let snapshot_atoms: Vec<u8> = topology
                .atoms
                .iter()
                .map(|a| a.element().atomic_number())
                .collect();
            let snapshot_bonds = topology.bonds.len();

            for (with_types, with_labels) in
                [(false, false), (true, false), (false, true), (true, true)]
            {
                real_calls += 1;
                let mut lp = vec![7.0_f64; atoms];
                let mut mr = vec![8.0_f64; atoms];
                let mut types = vec![91u32; atoms];
                let mut labels = vec!["sentinel-capacity".to_string(); atoms];
                // Capture EVERY label's pointer and capacity before.
                let ptrs: Vec<*const u8> = labels.iter().map(|l| l.as_ptr()).collect();
                let caps: Vec<usize> = labels.iter().map(|l| l.capacity()).collect();

                let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                let w0 = crate::crippen::CRIPPEN_PACKED_WORDS.with(|c| c.get());
                let v0 = crate::crippen::CRIPPEN_ROW_VISITS.with(|c| c.get());
                let lw0 = crate::crippen::CRIPPEN_LABEL_WRITES.with(|c| c.get());
                let valence0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                let topology0 = topology.clone();
                let coordinates0 = coordinates.clone();
                let properties0 = properties.clone();
                let assignment0 = assignment.clone();
                let rings0 = rings.clone();

                crate::crippen::crippen_cold_kernel(
                    &input,
                    &mut lp,
                    &mut mr,
                    with_types.then_some(&mut types),
                    with_labels.then_some(&mut labels),
                )
                .unwrap();

                assert_eq!(topology, topology0, "{name} full topology");
                assert_eq!(coordinates, coordinates0, "{name} coordinates");
                assert_eq!(properties, properties0, "{name} properties");
                assert_eq!(assignment, assignment0, "{name} supplied valence");
                assert_eq!(rings, rings0, "{name} supplied rings");
                assert_eq!(
                    VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                    valence0,
                    "{name} no descriptor-valence helper entry per call"
                );

                // Per-call counter deltas.
                assert_eq!(
                    crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - k0,
                    1
                );
                assert_eq!(
                    crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - c0,
                    1
                );
                assert_eq!(
                    crate::crippen::CRIPPEN_PACKED_WORDS.with(|c| c.get()) - w0,
                    expected_packed_words[shape_index],
                    "{name} packed words"
                );
                assert_eq!(
                    crate::crippen::CRIPPEN_ROW_VISITS.with(|c| c.get()) - v0,
                    row_visits,
                    "{name} row visits"
                );
                assert_eq!(
                    crate::crippen::CRIPPEN_LAST_WORD.with(|c| c.get()),
                    last_word,
                    "{name} last initial word"
                );
                let typed_count = expected_types.iter().filter(|&&t| t != 91).count();
                let expect_writes = if with_labels { typed_count } else { 0 };
                assert_eq!(
                    crate::crippen::CRIPPEN_LABEL_WRITES.with(|c| c.get()) - lw0,
                    expect_writes,
                    "{name} label writes"
                );

                // Numeric rows by literal bits: every TYPED atom is a
                // C1 row (0.1441/2.503); the untyped Cm keeps the
                // kernel's zeros; empty has no rows.
                let expect_lp: Vec<u64> = expected_types
                    .iter()
                    .map(|&t| {
                        if t == 91 {
                            0.0_f64.to_bits()
                        } else {
                            0.1441_f64.to_bits()
                        }
                    })
                    .collect();
                let expect_mr: Vec<u64> = expected_types
                    .iter()
                    .map(|&t| {
                        if t == 91 {
                            0.0_f64.to_bits()
                        } else {
                            2.503_f64.to_bits()
                        }
                    })
                    .collect();
                assert_eq!(bits(&lp), expect_lp, "{name} logP rows");
                assert_eq!(bits(&mr), expect_mr, "{name} MR rows");

                // Optional outputs: typed or untouched sentinel.
                if with_types {
                    assert_eq!(types, expected_types, "{name} typed indices");
                } else {
                    assert_eq!(types, vec![91u32; atoms], "{name} absent type sink");
                }
                assert_eq!(labels.len(), atoms, "{name} label rows");
                for (i, label) in labels.iter().enumerate() {
                    assert_eq!(label.as_ptr(), ptrs[i], "{name} label {i} pointer");
                    assert_eq!(label.capacity(), caps[i], "{name} label {i} capacity");
                }
                if with_labels {
                    // Vec length, every String pointer/capacity and
                    // content are stable-or-typed per the sink rule.
                    assert_eq!(labels.len(), atoms);
                    for (i, label) in labels.iter().enumerate() {
                        assert_eq!(label.as_ptr(), ptrs[i], "{name} label {i} pointer");
                        assert_eq!(label.capacity(), caps[i], "{name} label {i} capacity");
                        if expected_types[i] == 91 {
                            assert_eq!(label, "sentinel-capacity", "{name} untouched");
                        } else {
                            assert_eq!(label, "C1", "{name} typed label");
                        }
                    }
                } else {
                    for label in &labels {
                        assert_eq!(label, "sentinel-capacity", "{name} absent label sink");
                    }
                }

                // Input snapshot: nothing moved.
                assert_eq!(snapshot_bonds, topology.bonds.len(), "{name} bonds");
                assert_eq!(
                    snapshot_atoms,
                    topology
                        .atoms
                        .iter()
                        .map(|a| a.element().atomic_number())
                        .collect::<Vec<_>>(),
                    "{name} atoms"
                );
            }
            assert_eq!(
                VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                valence_before,
                "{name}"
            );
        }
        assert_eq!(real_calls, 20);
        // Across all twenty calls the kernel never re-prepares the
        // descriptor valence context (checked per shape above).
    }

    #[test]
    fn descriptor_crippen_fix_guard_malformed_preparation_product() {
        // Frozen CRIPPEN-CLOSE24 semantics: the 24-call
        // length{x} presence{x} force product on a malformed-valence
        // input — only the valid-length warm cache hit succeeds;
        // missing-MR takes precedence over length; every other arm
        // reaches the cold kernel and returns its REAL prepared error
        // — plus nine explicit optional-sink precondition rows that
        // fire BEFORE warm service.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("CC");
        // Malformed valence: three rows for a two-atom topology.
        let assignment = t01_zero_assignment(3);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 0);
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);

        let seed_lp = [11.0_f64, 12.0, 13.0];
        let seed_mr = [21.0_f64, 22.0, 23.0];
        let mut calls = 0u32;
        for len in [0usize, 2, 3] {
            for lp_present in [false, true] {
                for mr_present in [false, true] {
                    for force in [false, true] {
                        calls += 1;
                        let mut state = tpsa::DescriptorComputedState::default();
                        if lp_present {
                            state.crippen_slot_mut().logp_rows = Some(seed_lp[..len].to_vec());
                        }
                        if mr_present {
                            state.crippen_slot_mut().mr_rows = Some(seed_mr[..len].to_vec());
                        }
                        let result = crate::crippen::crippen_contribution_guard(
                            &input, force, &mut state, None, None,
                        );
                        if !force && lp_present && !mr_present {
                            // Missing MR takes precedence over length.
                            assert_eq!(
                                result.unwrap_err(),
                                DescriptorError::MissingCrippenMrContributions {
                                    function: "crippen_contributions"
                                }
                            );
                        } else if !force && lp_present && mr_present && len == 2 {
                            // The ONLY success: the valid-length hit.
                            let (lp, mr) = result.unwrap();
                            assert_eq!(lp, vec![11.0, 12.0]);
                            assert_eq!(mr, vec![21.0, 22.0]);
                        } else {
                            // Every other arm goes cold and returns the
                            // REAL prepared error (malformed valence),
                            // preserving prior state.
                            match result.unwrap_err() {
                                DescriptorError::Search { .. } => {}
                                other => panic!("expected prepared error, got {other:?}"),
                            }
                        }
                        // State preserved on every errored arm.
                        if !(!force && lp_present && mr_present && len == 2) {
                            let slot = state.crippen_slot();
                            assert_eq!(slot.logp_rows, lp_present.then(|| seed_lp[..len].to_vec()));
                            assert_eq!(slot.mr_rows, mr_present.then(|| seed_mr[..len].to_vec()));
                            assert_eq!(slot.logp, None);
                            assert_eq!(slot.mr, None);
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 24);

        // Optional invalid sink sizes {0,1,3} x {types,labels,both} are
        // explicit precondition rows that fire BEFORE warm service:
        // with a fully valid warm cache seeded, a bad sink still fails
        // with the typed shape error, never serving the cache.
        let good_assignment = valence(&topology, "cc24").unwrap();
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 0);
        let good_input = DescriptorInput::new(
            &topology,
            &coordinates,
            &properties,
            &good_assignment,
            &rings,
        );
        for bad_len in [0usize, 1, 3] {
            for (with_types, with_labels) in [(true, false), (false, true), (true, true)] {
                let mut state = tpsa::DescriptorComputedState::default();
                state.crippen_slot_mut().logp_rows = Some(vec![11.0, 12.0]);
                state.crippen_slot_mut().mr_rows = Some(vec![21.0, 22.0]);
                let mut types = vec![91u32; bad_len];
                let mut labels = vec!["sentinelA".to_string(); bad_len];
                let err = crate::crippen::crippen_contribution_guard(
                    &good_input,
                    false,
                    &mut state,
                    with_types.then_some(&mut types),
                    with_labels.then_some(&mut labels),
                )
                .unwrap_err();
                assert_eq!(
                    err,
                    DescriptorError::InvalidCrippenOptionalRows {
                        field: if with_types {
                            "atom_types"
                        } else {
                            "atom_labels"
                        },
                        actual: bad_len,
                        expected: 2,
                    }
                );
                // The warm cache is untouched by the precondition.
                let slot = state.crippen_slot();
                assert_eq!(slot.logp_rows, Some(vec![11.0, 12.0]));
                assert_eq!(slot.mr_rows, Some(vec![21.0, 22.0]));
            }
        }
    }

    #[test]
    fn descriptor_crippen_fix_contributions_ninety_six_call_product() {
        // Frozen CRIPPEN-CLOSE30 semantics: the complete 96-call
        // product on WELL-FORMED CC — row-length {0,2,3} x LogP
        // presence x MR presence x force x four optional-sink
        // combinations. Seed rows literal [11,12,13]/[21,22,23];
        // sinks [91,92]/["sentinelA","sentinelB"]. !force +
        // LP-present + MR-absent => missing-MR; !force + both-present
        // + length2 => exact seeded rows with sinks UNTOUCHED (zero
        // owner/preparation counter deltas); every other arm goes
        // cold => two 0.1441/2.503 rows with requested sinks typed
        // row1/C1. Assert calls == 96.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("CC");
        let assignment = valence(&topology, "cc30").unwrap();
        let rings = ring_info(&topology, "cc30").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        // Original, never-refreshed complete input baselines. The
        // independent per-call baselines below additionally bracket
        // each invocation together with its state and counters.
        let topology_snapshot = topology.clone();
        let coordinates_snapshot = coordinates.clone();
        let properties_snapshot = properties.clone();
        let assignment_snapshot = assignment.clone();
        let rings_snapshot = rings.clone();

        let seed_lp = [11.0_f64, 12.0, 13.0];
        let seed_mr = [21.0_f64, 22.0, 23.0];
        let mut calls = 0u32;
        for len in [0usize, 2, 3] {
            for lp_present in [false, true] {
                for mr_present in [false, true] {
                    for force in [false, true] {
                        for (with_types, with_labels) in
                            [(false, false), (true, false), (false, true), (true, true)]
                        {
                            calls += 1;
                            let mut state = tpsa::DescriptorComputedState::default();
                            if lp_present {
                                state.crippen_slot_mut().logp_rows = Some(seed_lp[..len].to_vec());
                            }
                            if mr_present {
                                state.crippen_slot_mut().mr_rows = Some(seed_mr[..len].to_vec());
                            }
                            let mut types = vec![91u32, 92];
                            let mut labels = vec!["sentinelA".to_string(), "sentinelB".to_string()];
                            let state_before = state.clone();

                            let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                            let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                            let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());

                            let topology_before = topology.clone();
                            let coordinates_before = coordinates.clone();
                            let properties_before = properties.clone();
                            let assignment_before = assignment.clone();
                            let rings_before = rings.clone();
                            let result = crippen_contributions(
                                &input,
                                force,
                                &mut state,
                                with_types.then_some(&mut types),
                                with_labels.then_some(&mut labels),
                            );

                            // Frozen CONTRIB-PROOF: the COMPLETE
                            // five-input comparison runs immediately
                            // AFTER every call, BEFORE branching — no
                            // arm can skip it.
                            assert_eq!(topology, topology_before, "per-call topology");
                            assert_eq!(coordinates, coordinates_before, "per-call coordinates");
                            assert_eq!(properties, properties_before, "per-call properties");
                            assert_eq!(assignment, assignment_before, "per-call valence");
                            assert_eq!(rings, rings_before, "per-call rings");
                            assert_eq!(topology, topology_snapshot, "topology");
                            assert_eq!(coordinates, coordinates_snapshot, "coordinates");
                            assert_eq!(properties, properties_snapshot, "properties");
                            assert_eq!(assignment, assignment_snapshot, "valence");
                            assert_eq!(rings, rings_snapshot, "rings");

                            let warm_hit = !force && lp_present && mr_present && len == 2;
                            if !force && lp_present && !mr_present {
                                assert_eq!(
                                    result.unwrap_err(),
                                    DescriptorError::MissingCrippenMrContributions {
                                        function: "crippen_contributions"
                                    }
                                );
                                // Frozen CONTRIB-PROOF: zero all-three
                                // counters, ENTIRE state equal to the
                                // pre-call clone, both sink sentinels
                                // unchanged regardless of requested
                                // sinks.
                                assert_eq!(
                                    crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()),
                                    k0
                                );
                                assert_eq!(
                                    crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()),
                                    c0
                                );
                                assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), v0);
                                assert_eq!(state, state_before);
                                assert_eq!(types, vec![91u32, 92]);
                                assert_eq!(
                                    labels,
                                    vec!["sentinelA".to_string(), "sentinelB".to_string()]
                                );
                            } else if warm_hit {
                                let contribs = result.unwrap();
                                assert_eq!(contribs.logp, vec![11.0, 12.0]);
                                assert_eq!(contribs.molar_refractivity, vec![21.0, 22.0]);
                                assert_eq!(types, vec![91u32, 92]);
                                assert_eq!(
                                    labels,
                                    vec!["sentinelA".to_string(), "sentinelB".to_string()]
                                );
                                assert_eq!(
                                    crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()),
                                    k0
                                );
                                assert_eq!(
                                    crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()),
                                    c0
                                );
                                assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), v0);
                                let slot = state.crippen_slot();
                                assert_eq!(slot.logp_rows, Some(vec![11.0, 12.0]));
                                assert_eq!(slot.mr_rows, Some(vec![21.0, 22.0]));
                                // Frozen CONTRIB-PROOF: the warm hit
                                // keeps ENTIRE state bit-identical
                                // (the pair above IS the seeded clone)
                                // — no continue, the common checks run.
                                assert_eq!(state, state_before);
                            } else {
                                let contribs = result.unwrap();
                                assert_eq!(contribs.logp, vec![0.1441_f64, 0.1441]);
                                assert_eq!(contribs.molar_refractivity, vec![2.503_f64, 2.503]);
                                if with_types {
                                    assert_eq!(types, vec![1u32, 1]);
                                } else {
                                    assert_eq!(types, vec![91u32, 92]);
                                }
                                if with_labels {
                                    assert_eq!(labels, vec!["C1".to_string(), "C1".to_string()]);
                                } else {
                                    assert_eq!(
                                        labels,
                                        vec!["sentinelA".to_string(), "sentinelB".to_string()]
                                    );
                                }
                                assert_eq!(
                                    crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - k0,
                                    1
                                );
                                // Frozen CONTRIB-PROOF: context delta
                                // 1, valence delta 0, ENTIRE state
                                // equals the pre-call clone with ONLY
                                // the two contribution fields replaced
                                // by the frozen cold pair.
                                assert_eq!(
                                    crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - c0,
                                    1
                                );
                                assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), v0);
                                let mut expected_state = state_before.clone();
                                expected_state.crippen_slot_mut().logp_rows =
                                    Some(vec![0.1441, 0.1441]);
                                expected_state.crippen_slot_mut().mr_rows =
                                    Some(vec![2.503, 2.503]);
                                assert_eq!(state, expected_state);
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 96);
    }

    #[test]
    fn descriptor_crippen_fix_hydrogens_working_copy_owner() {
        // Frozen CRIPPEN-CLOSE42 semantics: the ONE hydrogenated
        // working-copy owner — implicit/explicit/isotopic/zero-added-H
        // fixtures with the frozen supervisor H-bits reference table,
        // an ACTUAL invalid-coordinate AddHs failure returning the
        // real HydrogenError structurally, paired retained ring-row
        // order/member-size assertions, and complete input-storage
        // snapshots.
        let properties = MoleculeProperties::default();
        let bits = |x: f64| x.to_bits();

        // Frozen H-bits reference (independent pinned-RDKit source
        // table): implicit methane and ethane.
        for (smiles, h_logp_bits, h_mr_bits) in [
            ("C", 0x3fe4_5aee_631f_8a09_u64, 0x401a_ec8b_4395_8106_u64),
            ("CC", 0x3ff0_6b50_b0f2_7bb3_u64, 0x4026_b22d_0e56_041a_u64),
        ] {
            let coordinates = CoordinateBlock::default();
            let topology = smiles_topology(smiles);
            let assignment = valence(&topology, "cc42").unwrap();
            let rings = ring_info(&topology, "cc42").unwrap();
            let input =
                DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
            let snapshots = (
                topology.clone(),
                coordinates.clone(),
                properties.clone(),
                assignment.clone(),
                rings.clone(),
            );
            let mut parent = tpsa::DescriptorComputedState::default();
            let (logp, mr) =
                crate::crippen::crippen_hydrogenated_totals(&input, true, &mut parent).unwrap();
            assert_eq!(bits(logp), h_logp_bits, "{smiles} H logP bits");
            assert_eq!(bits(mr), h_mr_bits, "{smiles} H MR bits");
            // Parent publishes ONLY the two scalars; rows stay absent.
            let slot = parent.crippen_slot();
            assert_eq!(slot.logp.map(bits), Some(h_logp_bits));
            assert_eq!(slot.mr.map(bits), Some(h_mr_bits));
            assert_eq!(slot.logp_rows, None);
            assert_eq!(slot.mr_rows, None);
            // Input storage snapshots: all five preserved.
            assert_eq!(topology, snapshots.0);
            assert_eq!(coordinates, snapshots.1);
            assert_eq!(properties, snapshots.2);
            assert_eq!(assignment, snapshots.3);
            assert_eq!(rings, snapshots.4);
        }

        // Zero-added-H fixture: explicit-H methane adds NOTHING —
        // every hydrogen is already explicit; the totals equal the
        // implicit methane H bits (same 1 C1 + 4 H1 rows).
        let coordinates = CoordinateBlock::default();
        let topology = smiles_topology("[H]C([H])([H])[H]");
        let assignment = valence(&topology, "cc42z").unwrap();
        let rings = ring_info(&topology, "cc42z").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut parent = tpsa::DescriptorComputedState::default();
        let (logp, mr) =
            crate::crippen::crippen_hydrogenated_totals(&input, true, &mut parent).unwrap();
        assert_eq!(bits(logp), 0x3fe4_5aee_631f_8a09_u64);
        assert_eq!(bits(mr), 0x401a_ec8b_4395_8106_u64);

        // Isotopic explicit hydrogen: the computation runs through the
        // real owner and publishes finite totals with the parent
        // scalar-only contract.
        let topology = smiles_topology("[2H]C");
        let assignment = valence(&topology, "cc42i").unwrap();
        let rings = ring_info(&topology, "cc42i").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        let mut parent = tpsa::DescriptorComputedState::default();
        let (logp, mr) =
            crate::crippen::crippen_hydrogenated_totals(&input, true, &mut parent).unwrap();
        assert!(logp.is_finite());
        assert!(mr.is_finite());
        assert!(parent.crippen_slot().logp.is_some());
        assert_eq!(parent.crippen_slot().logp_rows, None);

        // Paired retained ring transport: cyclopropane — after addHs
        // the ONE original ring keeps its row ORDER and MEMBER SIZE
        // at the extended counts; no fabricated rings.
        let topology = smiles_topology("C1CC1");
        let assignment = valence(&topology, "cc42r").unwrap();
        let rings = ring_info(&topology, "cc42r").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        assert_eq!(rings.num_rings(), 1);
        let atom_ring_before: Vec<Vec<AtomId>> = rings.atom_rings().to_vec();
        let bond_ring_before: Vec<Vec<BondId>> = rings.bond_rings().to_vec();
        let mut parent = tpsa::DescriptorComputedState::default();
        crate::crippen::crippen_hydrogenated_totals(&input, true, &mut parent).unwrap();
        // The input rings are untouched (transport only reads them).
        assert_eq!(rings.atom_rings(), atom_ring_before.as_slice());
        assert_eq!(rings.bond_rings(), bond_ring_before.as_slice());

        // ACTUAL invalid-coordinate AddHs failure: one atom, one 2D
        // conformer with ZERO rows — the REAL owner returns
        // HydrogenError::InvalidCoordinates structurally; the parent
        // state and every input block are preserved.
        let atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let bad_coordinates = CoordinateBlock {
            conformers_2d: vec![cosmolkit_model::Conformer2D::new(0, vec![])],
            ..CoordinateBlock::default()
        };
        let assignment = t01_zero_assignment(1);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 1, 0);
        let input = DescriptorInput::new(
            &topology,
            &bad_coordinates,
            &properties,
            &assignment,
            &rings,
        );
        let snapshots = (
            topology.clone(),
            bad_coordinates.clone(),
            properties.clone(),
            assignment.clone(),
            rings.clone(),
        );
        let mut parent = tpsa::DescriptorComputedState::default();
        parent.crippen_slot_mut().logp = Some(123.0);
        let parent_before = parent.clone();
        use std::error::Error;
        let err =
            crate::crippen::crippen_hydrogenated_totals(&input, true, &mut parent).unwrap_err();
        let DescriptorError::Hydrogens { function, .. } = &err else {
            panic!("expected structural Hydrogens error, got {err:?}")
        };
        assert_eq!(*function, "crippen_totals");
        let down = err
            .source()
            .and_then(|s| s.downcast_ref::<cosmolkit_core::HydrogenError>())
            .expect("borrowed HydrogenError source");
        assert!(matches!(
            down,
            cosmolkit_core::HydrogenError::InvalidCoordinates(_)
        ));
        // Parent state and every input preserved on failure.
        assert_eq!(parent, parent_before);
        assert_eq!(topology, snapshots.0);
        assert_eq!(bad_coordinates, snapshots.1);
        assert_eq!(properties, snapshots.2);
        assert_eq!(assignment, snapshots.3);
        assert_eq!(rings, snapshots.4);
    }

    #[test]
    fn descriptor_crippen_fix_totals_thirty_two_call_product() {
        // Frozen CRIPPEN-CLOSE48 semantics: the 32-call product on CC
        // — scalar presences 2x2 x force 2 x includeHs 2 x seeded
        // parent contribution rows 2. LP-nonforce/MR-present serves
        // 789/456; LP-nonforce/MR-absent is the typed missing-MR;
        // otherwise no-H nonforce with valid seeded rows [11,12]/
        // [21,22] totals 23/43; all other no-H arms use the frozen
        // source no-H bits; every H arm uses the frozen H bits AND
        // leaves the parent rows unchanged (or absent) even when valid
        // parent contribution seeds exist.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("CC");
        let assignment = valence(&topology, "cc48").unwrap();
        let rings = ring_info(&topology, "cc48").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);

        const NO_H_LOGP: u64 = 0x3fd2_71de_69ad_42c4;
        const NO_H_MR: u64 = 0x4014_0624_dd2f_1aa0;
        const H_LOGP: u64 = 0x3ff0_6b50_b0f2_7bb3;
        const H_MR: u64 = 0x4026_b22d_0e56_041a;

        let mut calls = 0u32;
        for lp_present in [false, true] {
            for mr_present in [false, true] {
                for force in [false, true] {
                    for include_h in [false, true] {
                        for seeded_rows in [false, true] {
                            calls += 1;
                            let mut state = tpsa::DescriptorComputedState::default();
                            state.crippen_slot_mut().logp = lp_present.then_some(789.0);
                            state.crippen_slot_mut().mr = mr_present.then_some(456.0);
                            if seeded_rows {
                                state.crippen_slot_mut().logp_rows = Some(vec![11.0, 12.0]);
                                state.crippen_slot_mut().mr_rows = Some(vec![21.0, 22.0]);
                            }
                            state.slot_mut_for_tests(false).scalar = Some(123.0);
                            state.slot_mut_for_tests(false).contributions = Some(vec![1.0, 2.0]);
                            state.slot_mut_for_tests(true).scalar = Some(456.0);
                            state.slot_mut_for_tests(true).contributions = Some(vec![3.0, 4.0]);
                            state.labute_slot_mut().rows = Some(vec![4.0, 5.0]);
                            state.labute_slot_mut().hydrogens = Some(7.0);
                            state.labute_slot_mut().asa = Some(789.0);
                            let before = state.clone();
                            let snapshots = (
                                topology.clone(),
                                coordinates.clone(),
                                properties.clone(),
                                assignment.clone(),
                                rings.clone(),
                            );
                            let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                            let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                            let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                            let result = crippen_totals(&input, include_h, force, &mut state);
                            assert_eq!(topology, snapshots.0);
                            assert_eq!(coordinates, snapshots.1);
                            assert_eq!(properties, snapshots.2);
                            assert_eq!(assignment, snapshots.3);
                            assert_eq!(rings, snapshots.4);
                            let mut expected = before.clone();
                            let scalar_hit = !force && lp_present;
                            let contribution_hit = !include_h && !force && seeded_rows;
                            let cold = !scalar_hit && !contribution_hit;
                            if !scalar_hit {
                                let (lp_bits, mr_bits) = if include_h {
                                    (H_LOGP, H_MR)
                                } else if contribution_hit {
                                    (23.0_f64.to_bits(), 43.0_f64.to_bits())
                                } else {
                                    expected.crippen_slot_mut().logp_rows =
                                        Some(vec![0.1441, 0.1441]);
                                    expected.crippen_slot_mut().mr_rows = Some(vec![2.503, 2.503]);
                                    (NO_H_LOGP, NO_H_MR)
                                };
                                expected.crippen_slot_mut().logp = Some(f64::from_bits(lp_bits));
                                expected.crippen_slot_mut().mr = Some(f64::from_bits(mr_bits));
                            }
                            assert_eq!(state, expected);
                            assert_eq!(
                                crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - k0,
                                usize::from(cold)
                            );
                            assert_eq!(
                                crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - c0,
                                usize::from(cold)
                            );
                            assert_eq!(
                                VALENCE_HELPER_ENTRIES.with(|c| c.get()) - v0,
                                usize::from(!scalar_hit && include_h)
                            );
                            if !force && lp_present {
                                if mr_present {
                                    // Warm scalar serve — state bit-identical.
                                    let totals = result.unwrap();
                                    assert_eq!(totals.logp.to_bits(), 789.0_f64.to_bits());
                                    assert_eq!(
                                        totals.molar_refractivity.to_bits(),
                                        456.0_f64.to_bits()
                                    );
                                    assert_eq!(state, before);
                                } else {
                                    let err = result.unwrap_err();
                                    assert_eq!(
                                        err,
                                        DescriptorError::MissingCrippenMr {
                                            function: "crippen_totals"
                                        }
                                    );
                                    assert_eq!(state, before);
                                }
                                continue;
                            }
                            let totals = result.unwrap();
                            let slot = state.crippen_slot();
                            if include_h {
                                // Frozen H bits; parent rows unchanged
                                // (seeded pair kept or absent), even
                                // though valid seeds exist; scalars
                                // published after success.
                                assert_eq!(totals.logp.to_bits(), H_LOGP);
                                assert_eq!(totals.molar_refractivity.to_bits(), H_MR);
                                assert_eq!(slot.logp_rows, seeded_rows.then(|| vec![11.0, 12.0]));
                                assert_eq!(slot.mr_rows, seeded_rows.then(|| vec![21.0, 22.0]));
                                assert_eq!(slot.logp.map(f64::to_bits), Some(H_LOGP));
                                assert_eq!(slot.mr.map(f64::to_bits), Some(H_MR));
                            } else if !force && seeded_rows {
                                // Contribution warm hit: 11+12 / 21+22.
                                assert_eq!(
                                    totals.logp.to_bits(),
                                    (0.0_f64 + 11.0 + 12.0).to_bits()
                                );
                                assert_eq!(
                                    totals.molar_refractivity.to_bits(),
                                    (0.0_f64 + 21.0 + 22.0).to_bits()
                                );
                                assert_eq!(slot.logp_rows, Some(vec![11.0, 12.0]));
                                assert_eq!(slot.mr_rows, Some(vec![21.0, 22.0]));
                            } else {
                                // Cold no-H: frozen table bits; rows
                                // and scalars published as the cold pair.
                                assert_eq!(totals.logp.to_bits(), NO_H_LOGP);
                                assert_eq!(totals.molar_refractivity.to_bits(), NO_H_MR);
                                assert_eq!(slot.logp_rows, Some(vec![0.1441_f64, 0.1441]));
                                assert_eq!(slot.mr_rows, Some(vec![2.503_f64, 2.503]));
                                assert_eq!(slot.logp.map(f64::to_bits), Some(NO_H_LOGP));
                                assert_eq!(slot.mr.map(f64::to_bits), Some(NO_H_MR));
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 32);
    }

    #[test]
    fn descriptor_crippen_fix_hydrogens_sixteen_fixture_product() {
        // Frozen HYDROGEN-PROOF4 test 1 of 3: 16 calls = four literal
        // fixtures x force x parent-contribution-seed. Per-call:
        // five-input snapshots, observation-driven extended
        // counts/prefix IDs/rows, child == default(), deltas 1/1/1,
        // parent == before-clone with ONLY scalars replaced.
        let properties = MoleculeProperties::default();
        let cases: &[(&str, usize, usize, usize, usize, u64, u64)] = &[
            (
                "C",
                1,
                0,
                5,
                4,
                0x3fe4_5aee_631f_8a09,
                0x401a_ec8b_4395_8106,
            ),
            (
                "CC",
                2,
                1,
                8,
                7,
                0x3ff0_6b50_b0f2_7bb3,
                0x4026_b22d_0e56_041a,
            ),
            (
                "[H]C([H])([H])[H]",
                5,
                4,
                5,
                4,
                0x3fe4_5aee_631f_8a09,
                0x401a_ec8b_4395_8106,
            ),
            (
                "[2H]C",
                2,
                1,
                5,
                4,
                0x3fe4_5aee_631f_8a09,
                0x401a_ec8b_4395_8106,
            ),
        ];
        let mut calls = 0u32;
        for (smiles, orig_a, orig_b, ext_a, ext_b, lp_bits, mr_bits) in cases {
            let coordinates = CoordinateBlock::default();
            let topology = smiles_topology(smiles);
            assert_eq!(topology.atoms.len(), *orig_a, "{smiles} atoms");
            assert_eq!(topology.bonds.len(), *orig_b, "{smiles} bonds");
            if *smiles == "[2H]C" {
                // Deuterium prerequisite and preservation.
                assert_eq!(topology.atoms[0].isotope(), Some(2));
            }
            let assignment = valence(&topology, "hp4a").unwrap();
            let rings = ring_info(&topology, "hp4a").unwrap();
            let input =
                DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
            for force in [false, true] {
                for seed in [false, true] {
                    calls += 1;
                    let mut parent = tpsa::DescriptorComputedState::default();
                    parent.crippen_slot_mut().logp = Some(789.0);
                    parent.crippen_slot_mut().mr = Some(456.0);
                    if seed {
                        parent.crippen_slot_mut().logp_rows = Some(vec![777.0; *orig_a]);
                        parent.crippen_slot_mut().mr_rows = Some(vec![888.0; *orig_a]);
                    }
                    parent.slot_mut_for_tests(false).scalar = Some(123.0);
                    parent.slot_mut_for_tests(false).contributions = Some(vec![1.0, 2.0]);
                    parent.slot_mut_for_tests(true).scalar = Some(456.0);
                    parent.slot_mut_for_tests(true).contributions = Some(vec![3.0, 4.0]);
                    {
                        let slot = parent.labute_slot_mut();
                        slot.rows = Some(vec![4.0, 5.0]);
                        slot.hydrogens = Some(7.0);
                        slot.asa = Some(789.0);
                    }
                    let before = parent.clone();
                    let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                    let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                    let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                    let snapshots = (
                        topology.clone(),
                        coordinates.clone(),
                        properties.clone(),
                        assignment.clone(),
                        rings.clone(),
                    );
                    let (lp, mr) =
                        crate::crippen::crippen_hydrogenated_totals(&input, force, &mut parent)
                            .unwrap();
                    assert_eq!(lp.to_bits(), *lp_bits, "{smiles} LP bits");
                    assert_eq!(mr.to_bits(), *mr_bits, "{smiles} MR bits");
                    // Five inputs unchanged.
                    assert_eq!(topology, snapshots.0);
                    assert_eq!(coordinates, snapshots.1);
                    assert_eq!(properties, snapshots.2);
                    assert_eq!(assignment, snapshots.3);
                    assert_eq!(rings, snapshots.4);
                    // Deuterium preserved.
                    if *smiles == "[2H]C" {
                        assert_eq!(topology.atoms[0].isotope(), Some(2));
                    }
                    // ACTUAL observed consumer input.
                    let obs = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                        .with(|slot| slot.borrow().clone())
                        .expect("observation set");
                    assert_eq!(obs.atom_count, *ext_a, "{smiles} extended atoms");
                    assert_eq!(obs.bond_count, *ext_b, "{smiles} extended bonds");
                    assert_eq!(
                        obs.original_atom_ids,
                        topology.atoms.iter().map(|a| a.id()).collect::<Vec<_>>()
                    );
                    assert_eq!(
                        obs.original_bond_ids,
                        topology
                            .bonds
                            .iter()
                            .map(|bond| bond.id())
                            .collect::<Vec<_>>()
                    );
                    assert_eq!(
                        obs.original_isotopes,
                        topology
                            .atoms
                            .iter()
                            .map(|atom| atom.isotope())
                            .collect::<Vec<_>>()
                    );
                    if *smiles == "[2H]C" {
                        assert_eq!(obs.original_isotopes[0], Some(2));
                    }
                    assert_eq!(obs.atom_membership_count, *ext_a);
                    assert_eq!(obs.bond_membership_count, *ext_b);
                    assert_eq!(obs.atom_rings, rings.atom_rings());
                    assert_eq!(obs.bond_rings, rings.bond_rings());
                    assert_eq!(obs.child_state, tpsa::DescriptorComputedState::default());
                    // Parent: before-clone with ONLY scalars replaced.
                    let mut expected = before.clone();
                    expected.crippen_slot_mut().logp = Some(lp);
                    expected.crippen_slot_mut().mr = Some(mr);
                    assert_eq!(parent, expected);
                    // Deltas 1/1/1.
                    assert_eq!(
                        crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - k0,
                        1
                    );
                    assert_eq!(
                        crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - c0,
                        1
                    );
                    assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()) - v0, 1);
                }
            }
        }
        assert_eq!(calls, 16);
    }

    #[test]
    fn descriptor_crippen_fix_hydrogens_ring_transport_product() {
        // Frozen HYDROGEN-PROOF4 test 2 of 3: 4 calls = C1CC12CC2
        // (5/6 -> 13/14, TWO ring rows) x force x paired-row-order
        // (canonical order vs BOTH tables reversed). The observed
        // TRANSPORTED rows equal the full ordered supplied rows —
        // proof of the actual consumer input, not input immutability.
        let properties = MoleculeProperties::default();
        let coordinates = CoordinateBlock::default();
        let topology = smiles_topology("C1CC12CC2");
        assert_eq!(topology.atoms.len(), 5);
        assert_eq!(topology.bonds.len(), 6);
        let assignment = valence(&topology, "hp4b").unwrap();
        let canonical_rings = ring_info(&topology, "hp4b").unwrap();
        assert_eq!(canonical_rings.num_rings(), 2);
        let canonical_atom_rows: Vec<Vec<AtomId>> = canonical_rings.atom_rings().to_vec();
        let canonical_bond_rows: Vec<Vec<BondId>> = canonical_rings.bond_rings().to_vec();
        // BOTH tables reversed through the EXISTING selected-rows
        // constructor.
        let reversed_rings = cosmolkit_core::ring_info_from_selected_rows(
            5,
            6,
            &canonical_atom_rows
                .iter()
                .rev()
                .cloned()
                .collect::<Vec<_>>(),
            &canonical_bond_rows
                .iter()
                .rev()
                .cloned()
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let frozen: Vec<(Vec<Vec<AtomId>>, Vec<Vec<BondId>>)> = vec![
            (canonical_atom_rows.clone(), canonical_bond_rows.clone()),
            (
                reversed_rings.atom_rings().to_vec(),
                reversed_rings.bond_rings().to_vec(),
            ),
        ];
        let mut calls = 0u32;
        for (atom_rows, bond_rows) in &frozen {
            for force in [false, true] {
                calls += 1;
                let supplied =
                    cosmolkit_core::ring_info_from_selected_rows(5, 6, atom_rows, bond_rows)
                        .unwrap();
                let input = DescriptorInput::new(
                    &topology,
                    &coordinates,
                    &properties,
                    &assignment,
                    &supplied,
                );
                let snapshots = (
                    topology.clone(),
                    coordinates.clone(),
                    properties.clone(),
                    assignment.clone(),
                    supplied.clone(),
                );
                let mut parent = tpsa::DescriptorComputedState::default();
                parent.crippen_slot_mut().logp = Some(789.0);
                parent.crippen_slot_mut().mr = Some(456.0);
                parent.crippen_slot_mut().logp_rows = Some(vec![777.0; 5]);
                parent.crippen_slot_mut().mr_rows = Some(vec![888.0; 5]);
                parent.slot_mut_for_tests(false).scalar = Some(123.0);
                parent.slot_mut_for_tests(false).contributions = Some(vec![1.0, 2.0]);
                parent.slot_mut_for_tests(true).scalar = Some(456.0);
                parent.slot_mut_for_tests(true).contributions = Some(vec![3.0, 4.0]);
                parent.labute_slot_mut().rows = Some(vec![4.0, 5.0]);
                parent.labute_slot_mut().hydrogens = Some(7.0);
                parent.labute_slot_mut().asa = Some(789.0);
                let before = parent.clone();
                let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                let (lp, mr) =
                    crate::crippen::crippen_hydrogenated_totals(&input, force, &mut parent)
                        .unwrap();
                assert_eq!(lp.to_bits(), 0x3ff8_f765_fd8a_daba_u64);
                assert_eq!(mr.to_bits(), 0x4034_e6a7_ef9d_b22c_u64);
                assert_eq!(topology, snapshots.0);
                assert_eq!(coordinates, snapshots.1);
                assert_eq!(properties, snapshots.2);
                assert_eq!(assignment, snapshots.3);
                assert_eq!(supplied, snapshots.4);
                let obs = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                    .with(|slot| slot.borrow().clone())
                    .expect("observation set");
                assert_eq!(obs.atom_count, 13);
                assert_eq!(obs.bond_count, 14);
                // The TRANSPORTED rows equal the full ordered SUPPLIED
                // rows — the actual consumer input.
                assert_eq!(obs.atom_rings, *atom_rows);
                assert_eq!(obs.bond_rings, *bond_rows);
                assert_eq!(
                    obs.original_atom_ids,
                    topology.atoms.iter().map(|a| a.id()).collect::<Vec<_>>()
                );
                assert_eq!(
                    obs.original_bond_ids,
                    topology
                        .bonds
                        .iter()
                        .map(|bond| bond.id())
                        .collect::<Vec<_>>()
                );
                assert_eq!(obs.atom_membership_count, 13);
                assert_eq!(obs.bond_membership_count, 14);
                assert_eq!(obs.child_state, tpsa::DescriptorComputedState::default());
                let mut expected = before.clone();
                expected.crippen_slot_mut().logp = Some(lp);
                expected.crippen_slot_mut().mr = Some(mr);
                assert_eq!(parent, expected);
                assert_eq!(
                    crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - k0,
                    1
                );
                assert_eq!(
                    crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - c0,
                    1
                );
                assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()) - v0, 1);
            }
        }
        assert_eq!(calls, 4);
    }

    #[test]
    fn descriptor_crippen_fix_hydrogens_error_product() {
        // Frozen HYDROGEN-PROOF4 test 3 of 3: 4 calls = the REAL
        // invalid-coordinate fixture x force x parent-seed. The
        // structural error downcasts to HydrogenError and its stored
        // InvalidCoordinates reason; parent/inputs/observation
        // unchanged; all three counter deltas 0.
        let properties = MoleculeProperties::default();
        let atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms,
            bonds: Vec::new(),
            ..TopologyBlock::default()
        };
        let bad_coordinates = CoordinateBlock {
            conformers_2d: vec![cosmolkit_model::Conformer2D::new(0, vec![])],
            ..CoordinateBlock::default()
        };
        let assignment = t01_zero_assignment(1);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 1, 0);
        let input = DescriptorInput::new(
            &topology,
            &bad_coordinates,
            &properties,
            &assignment,
            &rings,
        );
        use std::error::Error;
        let mut calls = 0u32;
        for force in [false, true] {
            for seed in [false, true] {
                calls += 1;
                let mut parent = tpsa::DescriptorComputedState::default();
                parent.crippen_slot_mut().logp = Some(789.0);
                parent.crippen_slot_mut().mr = Some(456.0);
                if seed {
                    parent.crippen_slot_mut().logp_rows = Some(vec![777.0]);
                    parent.crippen_slot_mut().mr_rows = Some(vec![888.0]);
                }
                parent.slot_mut_for_tests(false).scalar = Some(123.0);
                parent.slot_mut_for_tests(false).contributions = Some(vec![1.0, 2.0]);
                parent.slot_mut_for_tests(true).scalar = Some(456.0);
                parent.slot_mut_for_tests(true).contributions = Some(vec![3.0, 4.0]);
                parent.labute_slot_mut().rows = Some(vec![4.0, 5.0]);
                parent.labute_slot_mut().hydrogens = Some(7.0);
                parent.labute_slot_mut().asa = Some(789.0);
                let before = parent.clone();
                let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                let observation_before = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                    .with(|slot| slot.borrow().clone());
                let snapshots = (
                    topology.clone(),
                    bad_coordinates.clone(),
                    properties.clone(),
                    assignment.clone(),
                    rings.clone(),
                );
                let err = crate::crippen::crippen_hydrogenated_totals(&input, force, &mut parent)
                    .unwrap_err();
                let DescriptorError::Hydrogens { function, .. } = &err else {
                    panic!("expected Hydrogens, got {err:?}")
                };
                assert_eq!(*function, "crippen_totals");
                let down = err
                    .source()
                    .and_then(|s| s.downcast_ref::<cosmolkit_core::HydrogenError>())
                    .expect("borrowed HydrogenError");
                let cosmolkit_core::HydrogenError::InvalidCoordinates(reason) = down else {
                    panic!("expected InvalidCoordinates, got {down:?}")
                };
                assert_eq!(
                    reason,
                    &cosmolkit_model::CoordinateValidationError::RowCount {
                        dimension: "2D",
                        conformer: 0,
                        rows: 0,
                        atom_count: 1,
                    }
                );
                assert_eq!(parent, before);
                assert_eq!(topology, snapshots.0);
                assert_eq!(bad_coordinates, snapshots.1);
                assert_eq!(properties, snapshots.2);
                assert_eq!(assignment, snapshots.3);
                assert_eq!(rings, snapshots.4);
                // Core currently makes HydrogenError an Error leaf. The
                // concrete coordinate cause is preserved in its variant
                // and asserted above; do not claim a second source link.
                assert!(down.source().is_none());
                // Early AddHs error preserves the prior observation verbatim.
                let obs = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                    .with(|slot| slot.borrow().clone());
                assert_eq!(obs, observation_before);
                assert_eq!(crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()), k0);
                assert_eq!(crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()), c0);
                assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), v0);
            }
        }
        assert_eq!(calls, 4);
    }

    #[test]
    fn descriptor_crippen_fix_projections_mode_sequences() {
        // Frozen PROJECTION-PROOF contract: 16 rows x 3 REAL calls =
        // 48 product calls + 1 clear + 8 cold reference = 57 TOTAL.
        // Every call brackets complete five inputs, full state and
        // kernel/context/valence counters with LITERAL expectations
        // from the source bit table, never production outputs.
        let properties = MoleculeProperties::default();
        let coordinates = CoordinateBlock::default();
        let topology = smiles_topology("CC");
        let assignment = valence(&topology, "cc54").unwrap();
        let rings = ring_info(&topology, "cc54").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);

        const H_LOGP: u64 = 0x3ff0_6b50_b0f2_7bb3;
        const H_MR: u64 = 0x4026_b22d_0e56_041a;
        const NO_H_LOGP: u64 = 0x3fd2_71de_69ad_42c4;
        const NO_H_MR: u64 = 0x4014_0624_dd2f_1aa0;
        let literal = |h: bool| {
            if h {
                (H_LOGP, H_MR)
            } else {
                (NO_H_LOGP, NO_H_MR)
            }
        };

        let mut product_calls = 0u32;
        for first_h in [false, true] {
            for second_h in [false, true] {
                for second_force in [false, true] {
                    for kind_logp in [true, false] {
                        // Seed ONLY unrelated state; Crippen all None.
                        let mut state = tpsa::DescriptorComputedState::default();
                        state.slot_mut_for_tests(false).scalar = Some(123.0);
                        state.slot_mut_for_tests(false).contributions = Some(vec![1.0, 2.0]);
                        state.slot_mut_for_tests(true).scalar = Some(456.0);
                        state.slot_mut_for_tests(true).contributions = Some(vec![3.0, 4.0]);
                        {
                            let slot = state.labute_slot_mut();
                            slot.rows = Some(vec![4.0, 5.0]);
                            slot.hydrogens = Some(7.0);
                            slot.asa = Some(789.0);
                        }
                        // Capture the PRE-call state and observation.
                        let snapshots = (
                            topology.clone(),
                            coordinates.clone(),
                            properties.clone(),
                            assignment.clone(),
                            rings.clone(),
                        );
                        let pre_first = state.clone();
                        let obs_before_1 = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                            .with(|slot| slot.borrow().clone());
                        // Call 1: first totals (always nonforce).
                        let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                        let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                        let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                        let first = crippen_totals(&input, first_h, false, &mut state).unwrap();
                        product_calls += 1;
                        assert_eq!(first.logp.to_bits(), literal(first_h).0, "1st literal");
                        assert_eq!(
                            first.molar_refractivity.to_bits(),
                            literal(first_h).1,
                            "1st literal"
                        );
                        // Build expected from the PRE-call clone; never
                        // clone post-call state as an expectation.
                        let mut expected = pre_first.clone();
                        {
                            let slot = expected.crippen_slot_mut();
                            slot.logp = Some(f64::from_bits(literal(first_h).0));
                            slot.mr = Some(f64::from_bits(literal(first_h).1));
                            if !first_h {
                                slot.logp_rows = Some(vec![0.1441, 0.1441]);
                                slot.mr_rows = Some(vec![2.503, 2.503]);
                            }
                        }
                        assert_eq!(state, expected, "1st whole state");
                        // Observation: warm/no-H unchanged.
                        let obs_after_1 = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                            .with(|slot| slot.borrow().clone());
                        if first_h {
                            let observed = obs_after_1.as_ref().expect("cold H observation set");
                            assert_eq!((observed.atom_count, observed.bond_count), (8, 7));
                            assert_eq!(observed.child_state, DescriptorComputedState::default());
                        } else {
                            assert_eq!(obs_after_1, obs_before_1, "no-H obs unchanged");
                        }
                        assert_eq!(
                            crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - k0,
                            1
                        );
                        assert_eq!(
                            crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - c0,
                            1
                        );
                        assert_eq!(
                            VALENCE_HELPER_ENTRIES.with(|c| c.get()) - v0,
                            if first_h { 1 } else { 0 }
                        );
                        assert_eq!(topology, snapshots.0);
                        assert_eq!(coordinates, snapshots.1);
                        assert_eq!(properties, snapshots.2);
                        assert_eq!(assignment, snapshots.3);
                        assert_eq!(rings, snapshots.4);
                        // Capture the PRE-second-call state and observation.
                        let snapshots = (
                            topology.clone(),
                            coordinates.clone(),
                            properties.clone(),
                            assignment.clone(),
                            rings.clone(),
                        );
                        let pre_second = state.clone();
                        let obs_before_2 = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                            .with(|slot| slot.borrow().clone());
                        // Call 2: second totals.
                        let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                        let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                        let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                        let second =
                            crippen_totals(&input, second_h, second_force, &mut state).unwrap();
                        product_calls += 1;
                        if !second_force {
                            // Warm: FIRST mode literal pair regardless
                            // of requested mode; state unchanged.
                            assert_eq!(second.logp.to_bits(), literal(first_h).0, "warm 1st pair");
                            assert_eq!(second.molar_refractivity.to_bits(), literal(first_h).1);
                            assert_eq!(second.logp.to_bits(), first.logp.to_bits());
                            assert_eq!(
                                second.molar_refractivity.to_bits(),
                                first.molar_refractivity.to_bits()
                            );
                            assert_eq!(state, pre_second);
                            assert_eq!(crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()), k0);
                            assert_eq!(crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()), c0);
                            assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), v0);
                        } else {
                            // Force: SECOND mode literal pair.
                            assert_eq!(second.logp.to_bits(), literal(second_h).0, "force 2nd");
                            assert_eq!(second.molar_refractivity.to_bits(), literal(second_h).1);
                            let mut expected = pre_second.clone();
                            {
                                let slot = expected.crippen_slot_mut();
                                slot.logp = Some(f64::from_bits(literal(second_h).0));
                                slot.mr = Some(f64::from_bits(literal(second_h).1));
                                if !second_h {
                                    slot.logp_rows = Some(vec![0.1441, 0.1441]);
                                    slot.mr_rows = Some(vec![2.503, 2.503]);
                                }
                            }
                            assert_eq!(state, expected, "force whole state");
                            assert_eq!(
                                crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - k0,
                                1
                            );
                            assert_eq!(
                                crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - c0,
                                1
                            );
                            assert_eq!(
                                VALENCE_HELPER_ENTRIES.with(|c| c.get()) - v0,
                                if second_h { 1 } else { 0 }
                            );
                        }
                        // All five inputs after the second call.
                        assert_eq!(topology, snapshots.0);
                        assert_eq!(coordinates, snapshots.1);
                        assert_eq!(properties, snapshots.2);
                        assert_eq!(assignment, snapshots.3);
                        assert_eq!(rings, snapshots.4);
                        // Observation: nonforce unchanged; cold H replaced.
                        let obs_after_2 = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                            .with(|slot| slot.borrow().clone());
                        if !second_force || !second_h {
                            assert_eq!(obs_after_2, obs_before_2, "2nd obs unchanged");
                        } else {
                            let observed =
                                obs_after_2.as_ref().expect("cold H 2nd observation set");
                            assert_eq!((observed.atom_count, observed.bond_count), (8, 7));
                            assert_eq!(observed.child_state, DescriptorComputedState::default());
                        }
                        // Capture the PRE-projection state and observation.
                        let snapshots = (
                            topology.clone(),
                            coordinates.clone(),
                            properties.clone(),
                            assignment.clone(),
                            rings.clone(),
                        );
                        let pre_proj = state.clone();
                        let obs_before_3 = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                            .with(|slot| slot.borrow().clone());
                        // Call 3: default projection (warm serve).
                        let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                        let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                        let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                        let served_pair = if second_force {
                            literal(second_h)
                        } else {
                            literal(first_h)
                        };
                        if kind_logp {
                            let lp = crippen_clogp(&input, &mut state).unwrap();
                            assert_eq!(lp.to_bits(), served_pair.0, "proj literal");
                            assert_eq!(lp.to_bits(), second.logp.to_bits());
                        } else {
                            let mr = crippen_mr(&input, &mut state).unwrap();
                            assert_eq!(mr.to_bits(), served_pair.1, "proj literal");
                            assert_eq!(mr.to_bits(), second.molar_refractivity.to_bits());
                        }
                        product_calls += 1;
                        assert_eq!(state, pre_proj, "proj state unchanged");
                        let obs_after_3 = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                            .with(|slot| slot.borrow().clone());
                        assert_eq!(obs_after_3, obs_before_3, "proj obs unchanged");
                        assert_eq!(crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()), k0);
                        assert_eq!(crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()), c0);
                        assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()), v0);
                        assert_eq!(topology, snapshots.0);
                        assert_eq!(coordinates, snapshots.1);
                        assert_eq!(properties, snapshots.2);
                        assert_eq!(assignment, snapshots.3);
                        assert_eq!(rings, snapshots.4);
                    }
                }
            }
        }
        assert_eq!(product_calls, 48, "exact product census");

        // Supplementary call 1 of 9: clear + cold default projection.
        let mut state = tpsa::DescriptorComputedState::default();
        state.slot_mut_for_tests(false).scalar = Some(123.0);
        state.slot_mut_for_tests(false).contributions = Some(vec![1.0, 2.0]);
        state.slot_mut_for_tests(true).scalar = Some(456.0);
        state.slot_mut_for_tests(true).contributions = Some(vec![3.0, 4.0]);
        {
            let slot = state.labute_slot_mut();
            slot.rows = Some(vec![4.0, 5.0]);
            slot.hydrogens = Some(7.0);
            slot.asa = Some(789.0);
        }
        state.crippen_slot_mut().logp = Some(789.0);
        state.crippen_slot_mut().mr = Some(456.0);
        state.clear();
        assert_eq!(state, tpsa::DescriptorComputedState::default(), "clear all");
        // Clone: complete state snapshot BEFORE mutating the clone.
        let original_state = state.clone();
        let mut clone = state.clone();
        clone.crippen_slot_mut().logp = Some(1.0);
        assert_eq!(state, original_state);
        assert_eq!(state, tpsa::DescriptorComputedState::default());
        let snapshots = (
            topology.clone(),
            coordinates.clone(),
            properties.clone(),
            assignment.clone(),
            rings.clone(),
        );
        let before_clear_call = state.clone();
        let _observation_before_clear_call =
            crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION.with(|slot| slot.borrow().clone());
        let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
        let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
        let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
        let lp = crippen_clogp(&input, &mut state).unwrap();
        assert_eq!(lp.to_bits(), H_LOGP, "clear restores H default");
        let mut expected = before_clear_call;
        expected.crippen_slot_mut().logp = Some(f64::from_bits(H_LOGP));
        expected.crippen_slot_mut().mr = Some(f64::from_bits(H_MR));
        assert_eq!(
            state, expected,
            "cold H repopulates scalars only, rows None"
        );
        assert_eq!(
            crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - k0,
            1
        );
        assert_eq!(
            crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - c0,
            1
        );
        assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()) - v0, 1);
        assert_eq!(topology, snapshots.0);
        assert_eq!(coordinates, snapshots.1);
        assert_eq!(properties, snapshots.2);
        assert_eq!(assignment, snapshots.3);
        assert_eq!(rings, snapshots.4);
        let clear_observed = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
            .with(|slot| slot.borrow().clone())
            .expect("cold cleared-state H observation");
        assert_eq!(
            (clear_observed.atom_count, clear_observed.bond_count),
            (8, 7)
        );
        assert_eq!(
            clear_observed.child_state,
            DescriptorComputedState::default()
        );

        // Supplementary calls 2-9: cold defaults on the frozen four
        // reference molecules — TWO fresh-state projections each.
        let mut supplementary = 1u32;
        let table: &[(&str, u64, u64, usize, usize)] = &[
            ("C", 0x3fe4_5aee_631f_8a09, 0x401a_ec8b_4395_8106, 5, 4),
            ("CC", H_LOGP, H_MR, 8, 7),
            ("CO", 0xbfd9_0e56_0418_9375, 0x4020_491d_14e3_bcd3, 6, 5),
            ("CCO", 0xbf56_f006_8db8_bb00, 0x4029_8504_816f_006a, 9, 8),
        ];
        for (smiles, lp_bits, mr_bits, h_atoms, h_bonds) in table {
            let mol_topology = smiles_topology(smiles);
            let mol_assignment = valence(&mol_topology, "pp4").unwrap();
            let mol_rings = ring_info(&mol_topology, "pp4").unwrap();
            let mol_input = DescriptorInput::new(
                &mol_topology,
                &coordinates,
                &properties,
                &mol_assignment,
                &mol_rings,
            );
            for use_clogp in [true, false] {
                let mol_snapshots = (
                    mol_topology.clone(),
                    coordinates.clone(),
                    properties.clone(),
                    mol_assignment.clone(),
                    mol_rings.clone(),
                );
                let mut st = tpsa::DescriptorComputedState::default();
                st.slot_mut_for_tests(false).scalar = Some(123.0);
                st.slot_mut_for_tests(false).contributions = Some(vec![1.0, 2.0]);
                st.slot_mut_for_tests(true).scalar = Some(456.0);
                st.slot_mut_for_tests(true).contributions = Some(vec![3.0, 4.0]);
                {
                    let slot = st.labute_slot_mut();
                    slot.rows = Some(vec![4.0, 5.0]);
                    slot.hydrogens = Some(7.0);
                    slot.asa = Some(789.0);
                }
                let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                // Capture observation BEFORE — never reset it.
                let _obs_before = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                    .with(|slot| slot.borrow().clone());
                let reference_before = st.clone();
                let out = if use_clogp {
                    crippen_clogp(&mol_input, &mut st).unwrap()
                } else {
                    crippen_mr(&mol_input, &mut st).unwrap()
                };
                supplementary += 1;
                let mut reference_expected = reference_before;
                reference_expected.crippen_slot_mut().logp = Some(f64::from_bits(*lp_bits));
                reference_expected.crippen_slot_mut().mr = Some(f64::from_bits(*mr_bits));
                assert_eq!(st, reference_expected, "{smiles} complete reference state");
                let expected_half = if use_clogp { *lp_bits } else { *mr_bits };
                assert_eq!(out.to_bits(), expected_half, "{smiles} literal");
                let slot = st.crippen_slot();
                assert_eq!(slot.logp.map(f64::to_bits), Some(*lp_bits));
                assert_eq!(slot.mr.map(f64::to_bits), Some(*mr_bits));
                assert_eq!(slot.logp_rows, None);
                assert_eq!(slot.mr_rows, None);
                assert_eq!(
                    st.slot_for_tests(false).scalar.map(f64::to_bits),
                    Some(123.0_f64.to_bits())
                );
                assert_eq!(
                    st.slot_for_tests(true).scalar.map(f64::to_bits),
                    Some(456.0_f64.to_bits())
                );
                assert_eq!(
                    st.labute_slot().asa.map(f64::to_bits),
                    Some(789.0_f64.to_bits())
                );
                assert_eq!(
                    st.slot_for_tests(false).contributions,
                    Some(vec![1.0_f64, 2.0])
                );
                assert_eq!(
                    st.slot_for_tests(true).contributions,
                    Some(vec![3.0_f64, 4.0])
                );
                assert_eq!(
                    crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - k0,
                    1
                );
                assert_eq!(
                    crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - c0,
                    1
                );
                assert_eq!(VALENCE_HELPER_ENTRIES.with(|c| c.get()) - v0, 1);
                assert_eq!(mol_topology, mol_snapshots.0);
                assert_eq!(coordinates, mol_snapshots.1);
                assert_eq!(properties, mol_snapshots.2);
                assert_eq!(mol_assignment, mol_snapshots.3);
                assert_eq!(mol_rings, mol_snapshots.4);
                let obs = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                    .with(|slot| slot.borrow().clone())
                    .expect("cold H observation (not reset)");
                assert_eq!(obs.atom_count, *h_atoms, "{smiles} H atoms");
                assert_eq!(obs.bond_count, *h_bonds, "{smiles} H bonds");
                assert_eq!(obs.child_state, DescriptorComputedState::default());
            }
        }
        assert_eq!(supplementary, 9, "exact supplementary census");
        assert_eq!(product_calls + supplementary, 57, "exact total census");
    }

    #[test]
    fn descriptor_crippen_fix_vsa_routing_product() {
        // Frozen VSA-PROOF contract: exactly 24 REAL calls (force2 x
        // rows2 x bins3 x kind2), counted AFTER each invocation; each
        // kind starts with its OWN fresh state. The output oracle is
        // the literal area bits 401db4e477267dcc — independently
        // derived from the pinned reference, never the Rust owner.
        const AREA_BITS: u64 = 0x401d_b4e4_7726_7dcc;
        const LABUTE_H_BITS: u64 = 0x3ff5_0067_0627_0806;
        const LABUTE_ASA_BITS: u64 = 0x4021_7a7f_1c58_1fe7;
        let properties = MoleculeProperties::default();
        let coordinates = CoordinateBlock::default();
        let topology = smiles_topology("C");
        let assignment = valence(&topology, "cc60").unwrap();
        let rings = ring_info(&topology, "cc60").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        // Supplementary arithmetic check (retained from the original).
        let dij_c_h = (0.77_f64 - 0.33).abs().max(0.77 + 0.33).min(0.77 + 0.33);
        let vi_c = 0.0_f64 + (0.33 * 0.33 - (0.77 - dij_c_h) * (0.77 - dij_c_h) / dij_c_h);
        let methane_area = std::f64::consts::PI * 0.77 * (4.0 * 0.77 - vi_c);
        assert_eq!(methane_area.to_bits(), AREA_BITS, "arithmetic check");

        let bin_cases: &[(Option<&[f64]>, usize, usize, usize, usize, usize, usize)] = &[
            // (bins, expected_len_lp, expected_len_mr, cold_lp_dest, warm_lp_dest, cold_mr_dest, warm_mr_dest)
            (None, 12, 10, 4, 11, 4, 9),
            (Some(&[1.0, 10.0]), 3, 3, 0, 1, 1, 2),
            (Some(&[]), 1, 1, 0, 0, 0, 0),
        ];

        let mut calls = 0u32;
        for force in [false, true] {
            for rows_present in [false, true] {
                for &(bins, len_lp, len_mr, cold_lp, warm_lp, cold_mr, warm_mr) in bin_cases {
                    for use_slogp in [true, false] {
                        // Each kind starts with its OWN state.
                        let mut state = tpsa::DescriptorComputedState::default();
                        state.slot_mut_for_tests(false).scalar = Some(123.0);
                        state.slot_mut_for_tests(false).contributions = Some(vec![1.0, 2.0]);
                        state.slot_mut_for_tests(true).scalar = Some(456.0);
                        state.slot_mut_for_tests(true).contributions = Some(vec![3.0, 4.0]);
                        state.crippen_slot_mut().logp = Some(789.0);
                        state.crippen_slot_mut().mr = Some(456.0);
                        if rows_present {
                            state.crippen_slot_mut().logp_rows = Some(vec![1.25]);
                            state.crippen_slot_mut().mr_rows = Some(vec![20.0]);
                        }
                        // Per-call five-input baselines (frozen
                        // contract: BEFORE EVERY call clone all five
                        // complete inputs — never outside the loop).
                        let topology_baseline = topology.clone();
                        let coordinates_baseline = coordinates.clone();
                        let properties_baseline = properties.clone();
                        let assignment_baseline = assignment.clone();
                        let rings_baseline = rings.clone();
                        let pre = state.clone();
                        let obs_before = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                            .with(|slot| slot.borrow().clone());
                        let k0 = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get());
                        let c0 = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get());
                        let v0 = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                        // THE call.
                        let res = if use_slogp {
                            slogp_vsa(&input, bins, force, &mut state).unwrap()
                        } else {
                            smr_vsa(&input, bins, force, &mut state).unwrap()
                        };
                        calls += 1;
                        let warm = !force && rows_present;
                        let expected_len = if use_slogp { len_lp } else { len_mr };
                        let dest = if warm {
                            if use_slogp { warm_lp } else { warm_mr }
                        } else if use_slogp {
                            cold_lp
                        } else {
                            cold_mr
                        };
                        assert_eq!(res.len(), expected_len, "length");
                        assert_eq!(res[dest].to_bits(), AREA_BITS, "area bits at dest {dest}");
                        for (i, v) in res.iter().enumerate() {
                            if i != dest {
                                assert_eq!(v.to_bits(), 0.0_f64.to_bits(), "bin {i} zero");
                            }
                        }
                        // ENTIRE expected state from the PRE-call clone.
                        let mut expected = pre.clone();
                        {
                            let slot = expected.labute_slot_mut();
                            slot.rows = Some(vec![f64::from_bits(AREA_BITS)]);
                            slot.hydrogens = Some(f64::from_bits(LABUTE_H_BITS));
                            slot.asa = Some(f64::from_bits(LABUTE_ASA_BITS));
                        }
                        if !warm {
                            // Cold publishes BOTH table rows (the
                            // canonical contribution owner always
                            // publishes logP AND MR rows).
                            expected.crippen_slot_mut().logp_rows = Some(vec![0.1441]);
                            expected.crippen_slot_mut().mr_rows = Some(vec![2.503]);
                        }
                        // Crippen scalars stay as seeded (VSA doesn't touch them).
                        // TPSA slots stay as seeded.
                        assert_eq!(state, expected, "whole state after call {calls}");
                        // Counter deltas: 0/0 warm, 1/1 cold; valence 0 always.
                        let kd = crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()) - k0;
                        let cd = crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()) - c0;
                        if warm {
                            assert_eq!(kd, 0, "warm kernel delta");
                            assert_eq!(cd, 0, "warm context delta");
                        } else {
                            assert_eq!(kd, 1, "cold kernel delta");
                            assert_eq!(cd, 1, "cold context delta");
                        }
                        assert_eq!(
                            VALENCE_HELPER_ENTRIES.with(|c| c.get()) - v0,
                            0,
                            "valence delta"
                        );
                        // Observation unchanged (VSA never triggers H expansion).
                        let obs_after = crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION
                            .with(|slot| slot.borrow().clone());
                        assert_eq!(obs_after, obs_before, "obs unchanged");
                        // All five inputs unchanged.
                        assert_eq!(topology, topology_baseline, "topology");
                        assert_eq!(coordinates, coordinates_baseline, "coordinates");
                        assert_eq!(properties, properties_baseline, "properties");
                        assert_eq!(assignment, assignment_baseline, "valence");
                        assert_eq!(rings, rings_baseline, "rings");
                    }
                }
            }
        }
        assert_eq!(calls, 24, "exact census");
    }

    #[test]
    fn descriptor_crippen_fix_scalar_presence_force_product() {
        // Frozen CRIPPEN-CLOSE36 semantics: the eight-call scalar
        // presence x force discriminator on the no-H route with a
        // MALFORMED prepared input — LP-present + nonforce returns
        // 789/456 when MR is present, missing-MR otherwise; all other
        // arms return the cold discriminator (None); the state stays
        // BIT-IDENTICAL through every call (the guard never mutates
        // and performs no preparation on any arm). Separate real
        // cold-path controls repeat the product on well-formed CC.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("CC");
        let assignment = t01_zero_assignment(3);
        let rings =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 2, 0);
        let _malformed =
            DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);

        let mut calls = 0u32;
        for lp_present in [false, true] {
            for mr_present in [false, true] {
                for force in [false, true] {
                    calls += 1;
                    let mut state = tpsa::DescriptorComputedState::default();
                    state.crippen_slot_mut().logp = lp_present.then_some(789.0);
                    state.crippen_slot_mut().mr = mr_present.then_some(456.0);
                    let before = state.clone();
                    let result = crate::crippen::crippen_scalar_guard(force, &state);
                    if !force && lp_present && mr_present {
                        assert_eq!(result.unwrap(), Some((789.0_f64, 456.0_f64)));
                    } else if !force && lp_present && !mr_present {
                        let err = result.unwrap_err();
                        assert_eq!(
                            err,
                            DescriptorError::MissingCrippenMr {
                                function: "crippen_totals"
                            }
                        );
                    } else {
                        assert_eq!(result.unwrap(), None, "cold discriminator");
                    }
                    assert_eq!(state, before);
                }
            }
        }
        assert_eq!(calls, 8);

        // Separate real cold-path controls on WELL-FORMED CC: the same
        // eight combinations with seeded rows also present — the
        // discriminator is purely presence/force-based and the state
        // (including rows) stays bit-identical.
        let topology = smiles_topology("CC");
        let assignment = valence(&topology, "cc36").unwrap();
        let rings = ring_info(&topology, "cc36").unwrap();
        let _well_formed =
            DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);
        for lp_present in [false, true] {
            for mr_present in [false, true] {
                for force in [false, true] {
                    let mut state = tpsa::DescriptorComputedState::default();
                    state.crippen_slot_mut().logp = lp_present.then_some(1.5);
                    state.crippen_slot_mut().mr = mr_present.then_some(2.5);
                    state.crippen_slot_mut().logp_rows = Some(vec![0.1441, 0.1441]);
                    state.crippen_slot_mut().mr_rows = Some(vec![2.503, 2.503]);
                    let before = state.clone();
                    let result = crate::crippen::crippen_scalar_guard(force, &state);
                    if !force && lp_present && mr_present {
                        assert_eq!(result.unwrap(), Some((1.5_f64, 2.5_f64)));
                    } else if !force && lp_present {
                        assert!(matches!(
                            result.unwrap_err(),
                            DescriptorError::MissingCrippenMr { .. }
                        ));
                    } else {
                        assert_eq!(result.unwrap(), None);
                    }
                    assert_eq!(state, before);
                }
            }
        }
    }

    #[test]
    fn descriptor_a01_owner_linkage_and_retained_patterns() {
        // Frozen A01 semantics (Step-460 audit): every exposed
        // convenience API equals its ONE detached *_prepared owner
        // through the SAME shared context helpers (no parallel
        // arithmetic); the retained pattern constants keep their
        // source identities; the strict/non-strict amide split keeps
        // the r02/r03 frozen counts.
        let topology = smiles_topology("CCO");
        let assignment = valence(&topology, "a01a").unwrap();
        let rings = ring_info(&topology, "a01a").unwrap();
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);

        // Owner-linkage matrix: wrapper == owner on the same shared
        // context, for every family that exposes both levels.
        assert_eq!(
            num_atoms(&topology).unwrap(),
            num_atoms_prepared(&input).unwrap()
        );
        assert_eq!(
            num_heavy_atoms(&topology).unwrap(),
            num_heavy_atoms_prepared(&input).unwrap()
        );
        assert_eq!(
            lipinski_hba(&topology).unwrap(),
            lipinski_hba_prepared(&input).unwrap()
        );
        assert_eq!(
            lipinski_hbd(&topology).unwrap(),
            lipinski_hbd_prepared(&input).unwrap()
        );
        assert_eq!(
            num_hbd(&topology).unwrap(),
            num_hbd_prepared(&input).unwrap()
        );
        assert_eq!(
            num_hba(&topology).unwrap(),
            num_hba_prepared(&input).unwrap()
        );
        assert_eq!(
            num_heteroatoms(&topology).unwrap(),
            num_heteroatoms_prepared(&input).unwrap()
        );
        assert_eq!(
            num_amide_bonds(&topology).unwrap(),
            num_amide_bonds_prepared(&input).unwrap()
        );
        assert_eq!(
            fraction_csp3(&topology).unwrap().to_bits(),
            fraction_csp3_prepared(&input).unwrap().to_bits()
        );
        assert_eq!(
            num_rings(&topology).unwrap(),
            num_rings_prepared(&input).unwrap()
        );
        assert_eq!(
            num_rotatable_bonds(&topology, RotatableBondsOptions::Strict).unwrap(),
            num_rotatable_bonds_prepared(&input, RotatableBondsOptions::Strict).unwrap()
        );
        assert_eq!(
            num_atom_stereo_centers(&topology).unwrap(),
            num_atom_stereo_centers_prepared(&input).unwrap()
        );
        assert_eq!(
            num_unspecified_atom_stereo_centers(&topology).unwrap(),
            num_unspecified_atom_stereo_centers_prepared(&input).unwrap()
        );
        let mut state = tpsa::DescriptorComputedState::default();
        assert_eq!(
            tpsa_from_topology(&topology, false).unwrap().to_bits(),
            tpsa(&input, false, false, &mut state).unwrap().to_bits()
        );

        // Retained patterns keep their frozen source identities.
        assert_eq!(
            NON_STRICT_ROTATABLE_PATTERN,
            "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]"
        );
        assert!(!STRICT_ROTATABLE_PATTERN.is_empty());
        assert!(!STRICT_LINKAGES_BASE_PATTERN.is_empty());
        assert!(!SYMMETRIC_RINGS_PATTERN.is_empty());
        assert_eq!(TERMINAL_TRIPLE_BONDS_PATTERN, "C#[#6,#7]");
        assert_eq!(NON_RING_AMIDES_PATTERN, "[C&!R](=O)NC");
        assert_eq!(NUM_HETEROATOMS_PATTERN, "[!#6;!#1]");

        // The amide strict/non-strict split (r02/r03 frozen counts):
        // NonStrict counts the amide C-N bond, Strict and Default
        // (which resolves to the pinned Strict build) do not.
        let topology = smiles_topology("CC(=O)NC");
        assert_eq!(
            num_rotatable_bonds(&topology, RotatableBondsOptions::NonStrict).unwrap(),
            1
        );
        assert_eq!(
            num_rotatable_bonds(&topology, RotatableBondsOptions::Strict).unwrap(),
            0
        );
        assert_eq!(
            num_rotatable_bonds(&topology, RotatableBondsOptions::Default).unwrap(),
            0
        );
    }

    #[test]
    fn descriptor_a01_cold_warm_and_force_semantics() {
        // Frozen A01 semantics (Step-460 audit, corrected per
        // supervisor): the Crippen force mirror (both arms use FRESH
        // states here, so both are cold computations and the identical
        // bits are cold-repeat evidence, not warm-cache evidence — the
        // four-field cache is the canonical crippen_contributions/
        // crippen_totals guards' arm); the Labute
        // cold/warm/force identity through the
        // typed slot; the tpsa clear() semantics; and the VSA default
        // products equal the same-owner manual composition.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("CO");
        let assignment = valence(&topology, "a01b").unwrap();
        let rings = ring_info(&topology, "a01b").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);

        // The Crippen force mirror: false and true arms produce the
        // identical bits (both arms carry FRESH states — cold repeats,
        // not warm-cache evidence; the warm/cold split is the canonical
        // guards' cache arm, exercised by the C06/C08 selector tests).
        let mut fresh_st_c = tpsa::DescriptorComputedState::default();
        let totals_cold = crippen_totals(&input, false, false, &mut fresh_st_c).unwrap();
        let (logp_cold, mr_cold) = (totals_cold.logp, totals_cold.molar_refractivity);
        let mut fresh_st_f = tpsa::DescriptorComputedState::default();
        let totals_forced = crippen_totals(&input, false, true, &mut fresh_st_f).unwrap();
        let (logp_forced, mr_forced) = (totals_forced.logp, totals_forced.molar_refractivity);
        assert_eq!(logp_cold.to_bits(), logp_forced.to_bits());
        assert_eq!(mr_cold.to_bits(), mr_forced.to_bits());

        // Labute cold/warm/force identity through the typed slot.
        let mut state = tpsa::DescriptorComputedState::default();
        let cold = labute_contributions(&input, true, false, &mut state).unwrap();
        let warm = labute_contributions(&input, true, false, &mut state).unwrap();
        let forced = labute_contributions(&input, true, true, &mut state).unwrap();
        assert_eq!(cold.atoms, warm.atoms);
        assert_eq!(cold.atoms, forced.atoms);
        assert_eq!(cold.hydrogens.to_bits(), warm.hydrogens.to_bits());
        assert_eq!(cold.hydrogens.to_bits(), forced.hydrogens.to_bits());

        // tpsa clear() empties the computed slots; recomputation
        // reproduces the same rows.
        let rows_cold = tpsa_contributions(&input, false, false, &mut state).unwrap();
        state.clear();
        let rows_after = tpsa_contributions(&input, false, false, &mut state).unwrap();
        assert_eq!(rows_cold, rows_after);
        state.clear();
        state.clear();

        // The VSA wrappers link to the same owners: their default
        // products equal the manual same-owner composition.
        let slogp = slogp_vsa(&input, None, false, &mut state).unwrap();
        let smr = smr_vsa(&input, None, false, &mut state).unwrap();
        assert_eq!(slogp.len(), 12);
        assert_eq!(smr.len(), 10);
        let lab = labute_contributions(&input, true, false, &mut state).unwrap();
        let mut st_c = tpsa::DescriptorComputedState::default();
        let crip = crippen_contributions(&input, true, &mut st_c, None, None).unwrap();
        assert_eq!(
            assign_contribs_to_bins(&lab.atoms, &crip.logp, &vsa::DEFAULT_SLOGP_BINS).unwrap(),
            slogp
        );
        assert_eq!(
            assign_contribs_to_bins(&lab.atoms, &crip.molar_refractivity, &vsa::DEFAULT_SMR_BINS)
                .unwrap(),
            smr
        );
    }

    #[test]
    fn descriptor_crippen_fix_surface_signature_controls() {
        // Frozen CLOSE66 surface controls: the canonical Crippen
        // contract is EXACTLY the crate-root re-exports. The four
        // entrypoints keep their frozen signatures (fn-pointer type
        // equality compiles only on an exact match), the result
        // structs expose exactly the two frozen fields, and the
        // parameter table keeps the `&'static [CrippenParamRow]`
        // shape. The E0603 private-module/helper proofs are the
        // named compile_fail doctests on the four entries; no
        // root/registry exposure is added here.
        let _: fn(
            &DescriptorInput<'_>,
            bool,
            &mut tpsa::DescriptorComputedState,
            Option<&mut Vec<u32>>,
            Option<&mut Vec<String>>,
        ) -> DescriptorResult<CrippenContributions> = crippen_contributions;
        let _: fn(
            &DescriptorInput<'_>,
            bool,
            bool,
            &mut tpsa::DescriptorComputedState,
        ) -> DescriptorResult<CrippenTotals> = crippen_totals;
        let _: fn(
            &DescriptorInput<'_>,
            &mut tpsa::DescriptorComputedState,
        ) -> DescriptorResult<f64> = crippen_clogp;
        let _: fn(
            &DescriptorInput<'_>,
            &mut tpsa::DescriptorComputedState,
        ) -> DescriptorResult<f64> = crippen_mr;
        let _: fn() -> &'static [CrippenParamRow] = default_crippen_params;

        // Result-struct field shapes: exactly the two frozen fields
        // per struct (rows for contributions, scalars for totals).
        let contribs = CrippenContributions {
            logp: Vec::new(),
            molar_refractivity: Vec::new(),
        };
        let _: Vec<f64> = contribs.logp;
        let _: Vec<f64> = contribs.molar_refractivity;
        let totals_value = CrippenTotals {
            logp: 0.0,
            molar_refractivity: 0.0,
        };
        let _: f64 = totals_value.logp;
        let _: f64 = totals_value.molar_refractivity;

        // The frozen 110-row default table.
        assert_eq!(default_crippen_params().len(), 110);

        // One REAL call through each entry binds the frozen
        // signatures to actual methane behavior: the no-H totals are
        // the C04 row literals; the H-mode projections match the H
        // arm of the SAME totals owner (header defaults,
        // Crippen.h:72-75).
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("C");
        let assignment = valence(&topology, "c66").unwrap();
        let rings = ring_info(&topology, "c66").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);

        let mut st_rows = tpsa::DescriptorComputedState::default();
        let rows = crippen_contributions(&input, false, &mut st_rows, None, None).unwrap();
        assert_eq!(rows.logp[0], 0.1441);
        assert_eq!(rows.molar_refractivity[0], 2.503);

        let mut st_no_h = tpsa::DescriptorComputedState::default();
        let no_h = crippen_totals(&input, false, false, &mut st_no_h).unwrap();
        assert_eq!(no_h.logp.to_bits(), 0.1441_f64.to_bits());
        assert_eq!(no_h.molar_refractivity.to_bits(), 2.503_f64.to_bits());

        let mut st_h = tpsa::DescriptorComputedState::default();
        let with_h = crippen_totals(&input, true, false, &mut st_h).unwrap();
        let mut st_clogp = tpsa::DescriptorComputedState::default();
        let clogp = crippen_clogp(&input, &mut st_clogp).unwrap();
        let mut st_mr = tpsa::DescriptorComputedState::default();
        let mr = crippen_mr(&input, &mut st_mr).unwrap();
        assert_eq!(clogp.to_bits(), with_h.logp.to_bits());
        assert_eq!(mr.to_bits(), with_h.molar_refractivity.to_bits());
    }

    #[test]
    fn descriptor_crippen_fix_sequences_clone_clear_force_atomicity() {
        // Frozen SEQUENCE-PROOF semantics (CCO; pinned Crippen.cpp
        // 35-133 + independently queried RDKit 2026.03.1): the TEN-call
        // sequence with its frozen (kernel, context, valence) delta
        // table — cold (1,1,0); warm (0,0,0); clone warm (0,0,0);
        // force (1,1,0); row-fed no-H totals (0,0,0); fresh H totals
        // (1,1,1); warm H totals (0,0,0); bad-types error (0,0,0);
        // missing-MR-rows error (0,0,0); post-clear (1,1,0). "Zero"
        // deltas mean zero kernel/context/valence REENTRY, not zero
        // allocation or total work: warm hits copy the cached rows, and
        // the shared context/matcher owners allocate at their own call
        // sites. BEFORE EVERY call the five complete supplied inputs,
        // the full destination state, all three counters and the
        // observation are captured; AFTER EVERY call the five inputs
        // are compared unchanged, the exact counter deltas asserted and
        // the WHOLE state compared against its pre-call clone plus the
        // frozen expected publication.
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let topology = smiles_topology("CCO");
        let assignment = valence(&topology, "c72").unwrap();
        let rings = ring_info(&topology, "c72").unwrap();
        let input = DescriptorInput::new(&topology, &coordinates, &properties, &assignment, &rings);

        // Frozen literals (pinned table rows + no-H/H scalar bits).
        const LP_ROWS: [f64; 3] = [0.1441, -0.2035, -0.2893];
        const MR_ROWS: [f64; 3] = [2.503, 2.753, 0.8238];
        const NO_H_LOGP_BITS: u64 = 0xbfd65119ce075f70;
        const NO_H_MR_BITS: u64 = 0x401851b71758e21a;
        const H_LOGP_BITS: u64 = 0xbf56f0068db8bb00;
        const H_MR_BITS: u64 = 0x40298504816f006a;

        let counters = || {
            (
                crate::crippen::CRIPPEN_KERNEL_CALLS.with(|c| c.get()),
                crate::crippen::CRIPPEN_CONTEXT_CALLS.with(|c| c.get()),
                VALENCE_HELPER_ENTRIES.with(|c| c.get()),
            )
        };
        let observation =
            || crate::crippen::HYDROGEN_WORKING_COPY_OBSERVATION.with(|s| s.borrow().clone());
        let snapshot = || {
            (
                topology.clone(),
                coordinates.clone(),
                properties.clone(),
                assignment.clone(),
                rings.clone(),
            )
        };
        let check_inputs = |s: &(
            TopologyBlock,
            CoordinateBlock,
            MoleculeProperties,
            ValenceAssignment,
            cosmolkit_core::RingInfo,
        )| {
            assert_eq!(topology, s.0, "input topology preserved");
            assert_eq!(coordinates, s.1, "input coordinates preserved");
            assert_eq!(properties, s.2, "input properties preserved");
            assert_eq!(assignment, s.3, "input valence preserved");
            assert_eq!(rings, s.4, "input ring info preserved");
        };
        let mut calls = 0usize;

        // Call 1 — cold contributions (1,1,0) on a state whose OTHER
        // slots hold sentinels (retention proven by the WHOLE-state
        // expected comparisons below, which derive from pre-call
        // clones).
        let mut state = tpsa::DescriptorComputedState::default();
        state.slot_mut_for_tests(false).scalar = Some(1.25);
        state.slot_mut_for_tests(false).contributions = Some(vec![9.0, 9.0, 9.0]);
        state.slot_mut_for_tests(true).scalar = Some(2.5);
        state.labute_slot_mut().rows = Some(vec![7.0, 7.0, 7.0]);
        state.labute_slot_mut().hydrogens = Some(7.5);
        state.labute_slot_mut().asa = Some(7.75);
        let snap = snapshot();
        let state_before = state.clone();
        let obs0 = observation();
        let (k0, c0, v0) = counters();
        let cold = crippen_contributions(&input, false, &mut state, None, None).unwrap();
        calls += 1;
        check_inputs(&snap);
        assert_eq!(counters(), (k0 + 1, c0 + 1, v0), "cold (1,1,0)");
        assert_eq!(observation(), obs0, "non-H call preserves observation");
        let mut expected = state_before.clone();
        expected.crippen_slot_mut().logp_rows = Some(LP_ROWS.to_vec());
        expected.crippen_slot_mut().mr_rows = Some(MR_ROWS.to_vec());
        assert_eq!(state, expected, "cold publishes ONLY the frozen rows");
        assert_eq!(cold.logp, LP_ROWS);
        assert_eq!(cold.molar_refractivity, MR_ROWS);

        // Call 2 — warm hit on the SAME state (0,0,0): zero
        // kernel/context/valence reentry (not zero work — the cached
        // rows are copied).
        let snap = snapshot();
        let state_before = state.clone();
        let obs0 = observation();
        let (k, c, v) = counters();
        let warm = crippen_contributions(&input, false, &mut state, None, None).unwrap();
        calls += 1;
        check_inputs(&snap);
        assert_eq!(counters(), (k, c, v), "warm (0,0,0)");
        assert_eq!(observation(), obs0);
        assert_eq!(state, state_before, "warm preserves the entire state");
        assert_eq!(warm, cold, "warm rows are bit-identical copies");

        // Call 3 — the CLONE preserves the cache and warm-hits (0,0,0).
        let mut clone_state = state.clone();
        let snap = snapshot();
        let clone_before = clone_state.clone();
        let obs0 = observation();
        let (k, c, v) = counters();
        let cloned = crippen_contributions(&input, false, &mut clone_state, None, None).unwrap();
        calls += 1;
        check_inputs(&snap);
        assert_eq!(counters(), (k, c, v), "clone warm (0,0,0)");
        assert_eq!(observation(), obs0);
        assert_eq!(clone_state, clone_before, "clone warm preserves state");
        assert_eq!(cloned, cold);

        // Call 4 — force=true bypasses the cache (1,1,0) and
        // republishes the frozen rows only; the ORIGINAL state (kept
        // independently) is untouched by the clone's forced refresh.
        let original_guard = state.clone();
        let snap = snapshot();
        let clone_before = clone_state.clone();
        let obs0 = observation();
        let (k, c, v) = counters();
        let forced = crippen_contributions(&input, true, &mut clone_state, None, None).unwrap();
        calls += 1;
        check_inputs(&snap);
        assert_eq!(counters(), (k + 1, c + 1, v), "force (1,1,0)");
        assert_eq!(observation(), obs0);
        let mut expected = clone_before.clone();
        expected.crippen_slot_mut().logp_rows = Some(LP_ROWS.to_vec());
        expected.crippen_slot_mut().mr_rows = Some(MR_ROWS.to_vec());
        assert_eq!(
            clone_state, expected,
            "force republishes the frozen rows only"
        );
        assert_eq!(forced, cold, "forced rows identical bits");
        assert_eq!(
            state, original_guard,
            "the clone's forced refresh cannot mutate the original state's fields"
        );

        // Call 5 — the warm ROW cache feeds the no-H totals path
        // (0,0,0): scalar miss -> row hit -> ordered sums -> ONLY the
        // two frozen no-H scalars published (rows and unrelated slots
        // retained).
        let snap = snapshot();
        let state_before = state.clone();
        let obs0 = observation();
        let (k, c, v) = counters();
        let totals = crippen_totals(&input, false, false, &mut state).unwrap();
        calls += 1;
        check_inputs(&snap);
        assert_eq!(counters(), (k, c, v), "row-fed no-H totals (0,0,0)");
        assert_eq!(observation(), obs0);
        let mut expected = state_before.clone();
        expected.crippen_slot_mut().logp = Some(f64::from_bits(NO_H_LOGP_BITS));
        expected.crippen_slot_mut().mr = Some(f64::from_bits(NO_H_MR_BITS));
        assert_eq!(
            state, expected,
            "no-H totals publish ONLY the frozen no-H scalars"
        );
        assert_eq!(totals.logp.to_bits(), NO_H_LOGP_BITS);
        assert_eq!(totals.molar_refractivity.to_bits(), NO_H_MR_BITS);
        // Supplementary arithmetic checks (the frozen bits above are
        // the actual oracles).
        let mut manual_logp = 0.0;
        for value in &cold.logp {
            manual_logp += value;
        }
        let mut manual_mr = 0.0;
        for value in &cold.molar_refractivity {
            manual_mr += value;
        }
        assert_eq!(totals.logp.to_bits(), manual_logp.to_bits());
        assert_eq!(totals.molar_refractivity.to_bits(), manual_mr.to_bits());

        // Call 6 — fresh H totals (1,1,1): the scalar guard must MISS
        // (scalars cleared), the parent's rows and other slots are
        // RETAINED, and ONLY the H scalar bits are published.
        let mut h_state = state.clone();
        h_state.crippen_slot_mut().logp = None;
        h_state.crippen_slot_mut().mr = None;
        let snap = snapshot();
        let h_before = h_state.clone();
        let (k, c, v) = counters();
        let h_cold = crippen_totals(&input, true, false, &mut h_state).unwrap();
        calls += 1;
        check_inputs(&snap);
        assert_eq!(counters(), (k + 1, c + 1, v + 1), "fresh H totals (1,1,1)");
        let mut expected = h_before.clone();
        expected.crippen_slot_mut().logp = Some(f64::from_bits(H_LOGP_BITS));
        expected.crippen_slot_mut().mr = Some(f64::from_bits(H_MR_BITS));
        assert_eq!(
            h_state, expected,
            "H totals publish ONLY the H scalars, retaining parent rows and other slots"
        );
        assert_eq!(h_cold.logp.to_bits(), H_LOGP_BITS);
        assert_eq!(h_cold.molar_refractivity.to_bits(), H_MR_BITS);
        // The ACTUAL child observation (source-defined equality checks
        // only — no value-inequality work inference).
        let obs_h = observation().expect("H working-copy observation is Some");
        assert_eq!(obs_h.atom_count, 9, "ethanol extended atoms");
        assert_eq!(obs_h.bond_count, 8, "ethanol extended bonds");
        assert_eq!(
            obs_h.atom_membership_count, 9,
            "membership dimensions match the atom count"
        );
        assert_eq!(
            obs_h.bond_membership_count, 8,
            "membership dimensions match the bond count"
        );
        assert_eq!(
            obs_h.original_atom_ids,
            vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]
        );
        assert_eq!(
            obs_h.original_bond_ids,
            vec![BondId::new(0), BondId::new(1)]
        );
        assert_eq!(obs_h.original_isotopes, vec![None, None, None]);
        assert!(obs_h.atom_rings.is_empty(), "empty paired atom ring rows");
        assert!(obs_h.bond_rings.is_empty(), "empty paired bond ring rows");
        assert_eq!(
            obs_h.child_state,
            tpsa::DescriptorComputedState::default(),
            "fresh cleared child state at the pre-contribution site"
        );

        // Call 7 — warm H totals (0,0,0): the scalar guard returns
        // with zero reentry; the observation is UNCHANGED by equality
        // (repeated identical inputs may yield equal observations — no
        // inequality is source-defined).
        let snap = snapshot();
        let h_before = h_state.clone();
        let (k, c, v) = counters();
        let obs_h_ref = observation();
        let h_warm = crippen_totals(&input, true, false, &mut h_state).unwrap();
        calls += 1;
        check_inputs(&snap);
        assert_eq!(counters(), (k, c, v), "warm H totals (0,0,0)");
        assert_eq!(h_state, h_before, "warm H preserves the entire state");
        assert_eq!(
            h_warm, h_cold,
            "cold-vs-warm identical bits (supplementary)"
        );
        assert_eq!(
            observation(),
            obs_h_ref,
            "no second working copy (equality)"
        );

        // Call 8 — bad optional-types error (0,0,0): the PRECONDITION
        // fires before any cache access; the WHOLE state, the sink
        // contents and the observation are preserved.
        let snap = snapshot();
        let before = state.clone();
        let obs0 = observation();
        let (k, c, v) = counters();
        let mut bad_types = vec![0u32; 2];
        let bad_types_before = bad_types.clone();
        let err = crippen_contributions(&input, false, &mut state, Some(&mut bad_types), None)
            .unwrap_err();
        calls += 1;
        check_inputs(&snap);
        assert_eq!(
            err,
            DescriptorError::InvalidCrippenOptionalRows {
                field: "atom_types",
                actual: 2,
                expected: 3,
            }
        );
        assert_eq!(counters(), (k, c, v), "precondition error (0,0,0)");
        assert_eq!(state, before, "precondition failure is state-atomic");
        assert_eq!(
            bad_types, bad_types_before,
            "optional sink contents preserved"
        );
        assert_eq!(observation(), obs0, "error preserves observation");

        // Call 9 — missing-MR-rows error (0,0,0): the warm read fails
        // BEFORE any kernel/context/valence work; WHOLE state and
        // observation preserved.
        let mut partial = state.clone();
        partial.crippen_slot_mut().mr_rows = None;
        let snap = snapshot();
        let partial_before = partial.clone();
        let obs0 = observation();
        let (k, c, v) = counters();
        let err = crippen_contributions(&input, false, &mut partial, None, None).unwrap_err();
        calls += 1;
        check_inputs(&snap);
        assert_eq!(
            err,
            DescriptorError::MissingCrippenMrContributions {
                function: "crippen_contributions"
            }
        );
        assert_eq!(counters(), (k, c, v), "warm-read error (0,0,0)");
        assert_eq!(
            partial, partial_before,
            "MR-missing failure is state-atomic"
        );
        assert_eq!(observation(), obs0, "error preserves observation");

        // Call 10 — clear() resets the WHOLE state to default (not
        // just selected fields); the next call is cold (1,1,0) and
        // rebuilds ONLY the Crippen rows.
        state.clear();
        assert_eq!(
            state,
            tpsa::DescriptorComputedState::default(),
            "clear resets the WHOLE state to default"
        );
        let snap = snapshot();
        let state_before = state.clone();
        let obs0 = observation();
        let (k, c, v) = counters();
        let after_clear = crippen_contributions(&input, false, &mut state, None, None).unwrap();
        calls += 1;
        check_inputs(&snap);
        assert_eq!(counters(), (k + 1, c + 1, v), "post-clear (1,1,0)");
        assert_eq!(observation(), obs0, "non-H call preserves observation");
        let mut expected = state_before.clone();
        expected.crippen_slot_mut().logp_rows = Some(LP_ROWS.to_vec());
        expected.crippen_slot_mut().mr_rows = Some(MR_ROWS.to_vec());
        assert_eq!(state, expected, "post-clear rebuilds ONLY the Crippen rows");
        assert_eq!(after_clear, cold, "post-clear rows identical bits");
        assert_eq!(calls, 10, "exact ten-call census");
    }

    /// Frozen DQ-PUBLIC reference table (installed RDKit 2026.03.1,
    /// sanitize=true, both removeHs policies, recorded BEFORE
    /// implementation): (SMILES, explicit rows with removeHs=false/true,
    /// heavy count, total count, Lipinski HBA, Lipinski HBD, CSP3 f64
    /// bits). Except the shown row-length contrast, all five outputs are
    /// identical under both parser policies; the policy-true row collapse
    /// comes from the real RemoveHs owner on the SAME original SMILES —
    /// neighbor H rows are removed, isotope-labeled H and lone H2 rows are
    /// retained.
    const PUBLIC_QUERY_CASES: [(&str, usize, usize, u32, u32, u32, u32, u64); 20] = [
        ("", 0, 0, 0, 0, 0, 0, 0x0000_0000_0000_0000),
        ("C", 1, 1, 1, 5, 0, 0, 0x3ff0_0000_0000_0000),
        ("CCO", 3, 3, 3, 9, 1, 1, 0x3ff0_0000_0000_0000),
        ("[NH4+]", 1, 1, 1, 5, 1, 4, 0x0000_0000_0000_0000),
        ("[O-]", 1, 1, 1, 1, 1, 0, 0x0000_0000_0000_0000),
        ("[H][H]", 2, 2, 0, 2, 0, 0, 0x0000_0000_0000_0000),
        ("[2H]O[2H]", 3, 3, 1, 3, 1, 2, 0x0000_0000_0000_0000),
        ("[13CH4]", 1, 1, 1, 5, 0, 0, 0x3ff0_0000_0000_0000),
        ("N", 1, 1, 1, 4, 1, 3, 0x0000_0000_0000_0000),
        ("O", 1, 1, 1, 3, 1, 2, 0x0000_0000_0000_0000),
        ("C=O", 2, 2, 2, 4, 1, 0, 0x0000_0000_0000_0000),
        ("O=C(N)N", 4, 4, 4, 8, 3, 4, 0x0000_0000_0000_0000),
        ("n1ccccc1", 6, 6, 6, 11, 1, 0, 0x0000_0000_0000_0000),
        ("[nH]1cccc1", 5, 5, 5, 10, 1, 1, 0x0000_0000_0000_0000),
        ("C1CCCCC1", 6, 6, 6, 18, 0, 0, 0x3ff0_0000_0000_0000),
        ("CC(N)C(=O)O", 6, 6, 6, 13, 3, 3, 0x3fe5_5555_5555_5555),
        ("[H]N([H])[H]", 4, 1, 1, 4, 1, 3, 0x0000_0000_0000_0000),
        ("C[C+](C)C", 4, 4, 4, 13, 0, 0, 0x3fe8_0000_0000_0000),
        ("*", 1, 1, 0, 1, 0, 0, 0x0000_0000_0000_0000),
        ("CC#N", 3, 3, 3, 6, 1, 0, 0x3fe0_0000_0000_0000),
    ];

    /// DQ fixture setup through the EXISTING detached owners: parse the
    /// ORIGINAL frozen SMILES; for remove_hydrogens=true run the core
    /// `remove_hydrogens_with_params` with `update_explicit_count=true`,
    /// `sanitize=true` on the parsed blocks; otherwise run the core
    /// `sanitize_topology` with default parameters. The final assignment is
    /// borrowed from the real setup result. Setup only — adapter calls
    /// never recompute.
    fn dq_fixture(
        smiles: &str,
        remove_hydrogens: bool,
    ) -> (
        cosmolkit_model::TopologyBlock,
        cosmolkit_core::ValenceAssignment,
    ) {
        let params = cosmolkit_smiles::SmilesParseParams {
            remove_hydrogens: false,
            ..Default::default()
        };
        let record = cosmolkit_smiles::parse_smiles(smiles, &params).expect("parse frozen SMILES");
        if remove_hydrogens {
            let result = cosmolkit_core::remove_hydrogens_with_params(
                record.topology,
                record.coordinates,
                record.properties,
                &cosmolkit_core::RemoveHsParams {
                    update_explicit_count: true,
                    sanitize: true,
                    ..cosmolkit_core::RemoveHsParams::default()
                },
            )
            .expect("RemoveHs owner on the original SMILES");
            let assignment = result
                .final_valence
                .expect("sanitized RemoveHs carries the final assignment");
            (result.topology, assignment)
        } else {
            let result = cosmolkit_core::sanitize_topology(
                &record.topology,
                &cosmolkit_core::SanitizeParams::default(),
            )
            .expect("sanitize owner on the original SMILES");
            let assignment = result
                .final_valence
                .expect("sanitize carries the final assignment");
            (result.topology, assignment)
        }
    }

    #[test]
    fn descriptor_public_dq01_total_atom_count_with_valence() {
        // DQ01 frozen product: all 40 sanitized case/policy rows through the
        // mandatory-borrow delegate. Valence is prepared ONLY during fixture
        // setup through the existing owner; every actual adapter call
        // asserts a ZERO valence-helper delta and unchanged borrowed inputs.
        let mut calls = 0usize;
        for (smiles, rows_false, rows_true, _, total, _, _, _) in PUBLIC_QUERY_CASES {
            for remove_hydrogens in [false, true] {
                let expected_rows = if remove_hydrogens {
                    rows_true
                } else {
                    rows_false
                };
                let (topology, assignment) = dq_fixture(smiles, remove_hydrogens);
                assert_eq!(
                    topology.atoms.len(),
                    expected_rows,
                    "{smiles} policy={remove_hydrogens} row count"
                );
                let topology_before = topology.clone();
                let assignment_before = assignment.clone();
                let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                let result = total_atom_count_with_valence(&topology, &assignment).unwrap();
                calls += 1;
                assert_eq!(
                    VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                    helper_before,
                    "{smiles} policy={remove_hydrogens}: adapter reassigns no valence"
                );
                assert_eq!(topology, topology_before, "borrowed topology unchanged");
                assert_eq!(
                    assignment, assignment_before,
                    "borrowed assignment unchanged"
                );
                assert_eq!(result, total, "{smiles} policy={remove_hydrogens}");
            }
        }
        assert_eq!(calls, 40, "exact 40-call census");

        // Malformed borrowed rows on CCO (3 atoms): explicit/implicit row
        // lengths 2/4 — four typed errors, inputs preserved, no valence
        // work. No fake healthy assignment is supplied to obtain results.
        let topology = smiles_topology("CCO");
        let healthy = valence(&topology, "dq01").unwrap();
        for (field, actual) in [
            ("explicit_valence", 2_usize),
            ("explicit_valence", 4_usize),
            ("implicit_hydrogens", 2_usize),
            ("implicit_hydrogens", 4_usize),
        ] {
            let mut bad = healthy.clone();
            let rows = if field == "explicit_valence" {
                &healthy.explicit_valence
            } else {
                &healthy.implicit_hydrogens
            };
            let mut value = rows.clone();
            value.resize(actual, *rows.last().unwrap_or(&0));
            if field == "explicit_valence" {
                bad.explicit_valence = value;
            } else {
                bad.implicit_hydrogens = value;
            }
            let topology_before = topology.clone();
            let bad_before = bad.clone();
            let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
            let err = total_atom_count_with_valence(&topology, &bad).unwrap_err();
            assert_eq!(
                VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                helper_before,
                "malformed call performs no valence work"
            );
            assert_eq!(topology, topology_before);
            assert_eq!(bad, bad_before);
            assert_eq!(
                err,
                DescriptorError::InvalidValenceRows {
                    function: "num_atoms",
                    field,
                    actual,
                    expected: 3,
                },
                "{field}={actual}"
            );
        }
    }

    #[test]
    fn descriptor_public_dq02_lipinski_hbd_with_valence() {
        // DQ02 frozen product: all 40 sanitized case/policy rows through the
        // mandatory-borrow donor-hydrogen-sum delegate (N/O rows only; the
        // hydrogen sum, not the donor-atom count). Valence is prepared ONLY
        // during fixture setup; each adapter call asserts a ZERO
        // valence-helper delta and unchanged borrowed inputs. Frozen
        // discriminators retained: [2H]O[2H]=2 (isotopic H NEIGHBORS count),
        // [NH4+]=4, [nH]1cccc1=1, C=O=0, [H]N([H])[H]=3 under BOTH policies.
        let mut calls = 0usize;
        for (smiles, rows_false, rows_true, _, _, _, hbd, _) in PUBLIC_QUERY_CASES {
            for remove_hydrogens in [false, true] {
                let expected_rows = if remove_hydrogens {
                    rows_true
                } else {
                    rows_false
                };
                let (topology, assignment) = dq_fixture(smiles, remove_hydrogens);
                assert_eq!(
                    topology.atoms.len(),
                    expected_rows,
                    "{smiles} policy={remove_hydrogens} row count"
                );
                let topology_before = topology.clone();
                let assignment_before = assignment.clone();
                let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                let result = lipinski_hbd_with_valence(&topology, &assignment).unwrap();
                calls += 1;
                assert_eq!(
                    VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                    helper_before,
                    "{smiles} policy={remove_hydrogens}: adapter reassigns no valence"
                );
                assert_eq!(topology, topology_before, "borrowed topology unchanged");
                assert_eq!(
                    assignment, assignment_before,
                    "borrowed assignment unchanged"
                );
                assert_eq!(result, hbd, "{smiles} policy={remove_hydrogens}");
            }
        }
        assert_eq!(calls, 40, "exact 40-call census");

        // Malformed borrowed rows on CCO: the same four typed errors with
        // the DQ02 function tag; inputs preserved; zero valence work.
        let topology = smiles_topology("CCO");
        let healthy = valence(&topology, "dq02").unwrap();
        for (field, actual) in [
            ("explicit_valence", 2_usize),
            ("explicit_valence", 4_usize),
            ("implicit_hydrogens", 2_usize),
            ("implicit_hydrogens", 4_usize),
        ] {
            let mut bad = healthy.clone();
            let rows = if field == "explicit_valence" {
                &healthy.explicit_valence
            } else {
                &healthy.implicit_hydrogens
            };
            let mut value = rows.clone();
            value.resize(actual, *rows.last().unwrap_or(&0));
            if field == "explicit_valence" {
                bad.explicit_valence = value;
            } else {
                bad.implicit_hydrogens = value;
            }
            let topology_before = topology.clone();
            let bad_before = bad.clone();
            let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
            let err = lipinski_hbd_with_valence(&topology, &bad).unwrap_err();
            assert_eq!(
                VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                helper_before,
                "malformed call performs no valence work"
            );
            assert_eq!(topology, topology_before);
            assert_eq!(bad, bad_before);
            assert_eq!(
                err,
                DescriptorError::InvalidValenceRows {
                    function: "lipinski_hbd",
                    field,
                    actual,
                    expected: 3,
                },
                "{field}={actual}"
            );
        }
    }

    #[test]
    fn descriptor_public_dq03_fraction_csp3_with_valence() {
        // DQ03 frozen product: all 40 sanitized case/policy rows through the
        // mandatory-borrow CSP3-fraction delegate, compared by LITERAL f64
        // BITS. The criterion is source total degree (neighbors + H with
        // includeNeighbors=false), never a hybridization flag. Frozen
        // discriminators retained: zero-carbon rows return +0.0 bits
        // (empty, [O-], [H][H], [2H]O[2H], N, O, C=O, O=C(N)N, aromatics,
        // wildcard), the charged carbon C[C+](C)C yields 0x3fe8000000000000
        // (3/4), and CC#N yields 0x3fe0000000000000 (1/2).
        let mut calls = 0usize;
        for (smiles, rows_false, rows_true, _, _, _, _, csp3_bits) in PUBLIC_QUERY_CASES {
            for remove_hydrogens in [false, true] {
                let expected_rows = if remove_hydrogens {
                    rows_true
                } else {
                    rows_false
                };
                let (topology, assignment) = dq_fixture(smiles, remove_hydrogens);
                assert_eq!(
                    topology.atoms.len(),
                    expected_rows,
                    "{smiles} policy={remove_hydrogens} row count"
                );
                let topology_before = topology.clone();
                let assignment_before = assignment.clone();
                let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
                let result = fraction_csp3_with_valence(&topology, &assignment).unwrap();
                calls += 1;
                assert_eq!(
                    VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                    helper_before,
                    "{smiles} policy={remove_hydrogens}: adapter reassigns no valence"
                );
                assert_eq!(topology, topology_before, "borrowed topology unchanged");
                assert_eq!(
                    assignment, assignment_before,
                    "borrowed assignment unchanged"
                );
                assert_eq!(
                    result.to_bits(),
                    csp3_bits,
                    "{smiles} policy={remove_hydrogens} exact bits"
                );
            }
        }
        assert_eq!(calls, 40, "exact 40-call census");

        // Malformed borrowed rows on CCO: four typed errors with the DQ03
        // function tag; inputs preserved; zero valence work.
        let topology = smiles_topology("CCO");
        let healthy = valence(&topology, "dq03").unwrap();
        for (field, actual) in [
            ("explicit_valence", 2_usize),
            ("explicit_valence", 4_usize),
            ("implicit_hydrogens", 2_usize),
            ("implicit_hydrogens", 4_usize),
        ] {
            let mut bad = healthy.clone();
            let rows = if field == "explicit_valence" {
                &healthy.explicit_valence
            } else {
                &healthy.implicit_hydrogens
            };
            let mut value = rows.clone();
            value.resize(actual, *rows.last().unwrap_or(&0));
            if field == "explicit_valence" {
                bad.explicit_valence = value;
            } else {
                bad.implicit_hydrogens = value;
            }
            let topology_before = topology.clone();
            let bad_before = bad.clone();
            let helper_before = VALENCE_HELPER_ENTRIES.with(|c| c.get());
            let err = fraction_csp3_with_valence(&topology, &bad).unwrap_err();
            assert_eq!(
                VALENCE_HELPER_ENTRIES.with(|c| c.get()),
                helper_before,
                "malformed call performs no valence work"
            );
            assert_eq!(topology, topology_before);
            assert_eq!(bad, bad_before);
            assert_eq!(
                err,
                DescriptorError::InvalidValenceRows {
                    function: "fraction_csp3",
                    field,
                    actual,
                    expected: 3,
                },
                "{field}={actual}"
            );
        }
    }
}
