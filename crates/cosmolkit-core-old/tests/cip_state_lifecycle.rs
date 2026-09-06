use cosmolkit_core::{
    AddHsParams, AtomSpec, BondOrder, BondSpec, ChiralTag, CipStatePolicy, Element, MOLECULE_OPS,
    Molecule, MoleculeBuilder, OperationDomain, RemoveHsParams, SanitizeOps,
    chemistry::tautomer::TautomerEnumerator, mol_from_binary, mol_to_binary,
};

#[derive(Debug, Clone, PartialEq, Eq)]
struct CipState {
    computed: Option<String>,
    computed_is_computed: bool,
    atom_codes: Vec<Option<String>>,
    atom_code_is_computed: Vec<bool>,
    atom_neighbor_orders: Vec<Option<String>>,
    atom_neighbor_order_is_computed: Vec<bool>,
    atom_ranks: Vec<Option<String>>,
    atom_rank_is_computed: Vec<bool>,
    bond_codes: Vec<Option<String>>,
    bond_code_is_computed: Vec<bool>,
    bond_neighbor_orders: Vec<Option<String>>,
    bond_neighbor_order_is_computed: Vec<bool>,
}

impl CipState {
    fn from_molecule(molecule: &Molecule) -> Self {
        Self {
            computed: molecule.prop("_CIPComputed").map(str::to_owned),
            computed_is_computed: molecule.is_prop_computed("_CIPComputed"),
            atom_codes: molecule
                .atoms()
                .iter()
                .map(|atom| atom.prop("_CIPCode").map(str::to_owned))
                .collect(),
            atom_code_is_computed: molecule
                .atoms()
                .iter()
                .map(|atom| atom.is_prop_computed("_CIPCode"))
                .collect(),
            atom_neighbor_orders: molecule
                .atoms()
                .iter()
                .map(|atom| atom.prop("_CIPNeighborOrder").map(str::to_owned))
                .collect(),
            atom_neighbor_order_is_computed: molecule
                .atoms()
                .iter()
                .map(|atom| atom.is_prop_computed("_CIPNeighborOrder"))
                .collect(),
            atom_ranks: molecule
                .atoms()
                .iter()
                .map(|atom| atom.prop("_CIPRank").map(str::to_owned))
                .collect(),
            atom_rank_is_computed: molecule
                .atoms()
                .iter()
                .map(|atom| atom.is_prop_computed("_CIPRank"))
                .collect(),
            bond_codes: molecule
                .bonds()
                .iter()
                .map(|bond| bond.prop("_CIPCode").map(str::to_owned))
                .collect(),
            bond_code_is_computed: molecule
                .bonds()
                .iter()
                .map(|bond| bond.is_prop_computed("_CIPCode"))
                .collect(),
            bond_neighbor_orders: molecule
                .bonds()
                .iter()
                .map(|bond| bond.prop("_CIPNeighborOrder").map(str::to_owned))
                .collect(),
            bond_neighbor_order_is_computed: molecule
                .bonds()
                .iter()
                .map(|bond| bond.is_prop_computed("_CIPNeighborOrder"))
                .collect(),
        }
    }

    fn assert_computed_cleared(&self) {
        assert_eq!(self.computed, None);
        assert!(!self.computed_is_computed);
        assert!(self.atom_neighbor_orders.iter().all(Option::is_none));
        assert!(
            self.atom_neighbor_order_is_computed
                .iter()
                .all(|computed| !computed)
        );
        assert!(self.atom_ranks.iter().all(Option::is_none));
        assert!(self.atom_rank_is_computed.iter().all(|computed| !computed));
        assert!(self.bond_neighbor_orders.iter().all(Option::is_none));
        assert!(
            self.bond_neighbor_order_is_computed
                .iter()
                .all(|computed| !computed)
        );
    }
}

fn lifecycle_seed() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let center = builder.add_atom(
        AtomSpec::new(Element::C)
            .with_chiral_tag(ChiralTag::TetrahedralCcw)
            .with_prop("_CIPCode", "S")
            .with_computed_prop("_CIPNeighborOrder", "[1,2,3,4]")
            .with_computed_prop("_CIPRank", "17"),
    );
    let fluorine = builder.add_atom(AtomSpec::new(Element::F));
    let chlorine = builder.add_atom(AtomSpec::new(Element::CL));
    let bromine = builder.add_atom(AtomSpec::new(Element::BR));
    let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
    for (offset, neighbor) in [fluorine, chlorine, bromine, hydrogen]
        .into_iter()
        .enumerate()
    {
        let mut bond = BondSpec::new(center, neighbor, BondOrder::Single);
        if offset == 0 {
            bond = bond
                .with_prop("_CIPCode", "E")
                .with_computed_prop("_CIPNeighborOrder", "[0,2]");
        }
        builder.add_bond(bond).expect("lifecycle bond must build");
    }
    builder
        .add_3d_conformer(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ])
        .expect("lifecycle coordinates must align");
    builder
        .build()
        .expect("lifecycle molecule must build")
        .with_assigned_valence()
        .expect("lifecycle molecule valence must be assignable")
        .with_cip_labels()
        .expect("lifecycle molecule CIP state must be assignable")
}

fn run_topology_operation(method: &str, molecule: &Molecule) -> Vec<Molecule> {
    let single = match method {
        "with_hydrogens_with_params" => molecule.with_hydrogens_with_params(AddHsParams::default()),
        "without_hydrogens_with_sanitize" => molecule.without_hydrogens_with_sanitize(true),
        "without_hydrogens_with_params" => {
            molecule.without_hydrogens_with_params(RemoveHsParams::default(), true)
        }
        "with_kekulized_bonds" => molecule.with_kekulized_bonds(false),
        "sanitize_with_ops" => molecule.sanitize_with_ops(SanitizeOps::ALL),
        "with_assigned_valence_strict" => molecule.with_assigned_valence_strict(true),
        "with_assigned_rings" => molecule.with_assigned_rings(),
        "with_assigned_ring_families" => molecule.with_assigned_ring_families(),
        "with_assigned_aromaticity" => molecule.with_assigned_aromaticity(),
        "with_assigned_radicals" => molecule.with_assigned_radicals(),
        "with_chiral_tags_from_structure" => molecule.with_chiral_tags_from_structure(-1, true),
        "with_cip_labels_with_options" => molecule.with_cip_labels(),
        "enumerate_tautomers_with_options" => {
            return TautomerEnumerator::new()
                .enumerate(molecule)
                .unwrap_or_else(|error| panic!("topology operation `{method}` failed: {error}"))
                .molecules();
        }
        other => panic!("topology operation `{other}` is missing from the CIP lifecycle runner"),
    };
    vec![single.unwrap_or_else(|error| panic!("topology operation `{method}` failed: {error}"))]
}

#[test]
fn cip_state_lifecycle_matrix_executes_every_registered_topology_operation() {
    let source = lifecycle_seed();
    let source_state = CipState::from_molecule(&source);
    let mut executed = 0usize;

    for operation in MOLECULE_OPS
        .iter()
        .copied()
        .filter(|operation| operation.domain == OperationDomain::Topology)
    {
        executed += 1;
        let results = run_topology_operation(operation.method, &source);
        assert!(
            !results.is_empty(),
            "{} emitted no branches",
            operation.method
        );
        for result in results {
            let state = CipState::from_molecule(&result);
            match operation.cip_state {
                CipStatePolicy::Preserve => {
                    assert_eq!(state, source_state, "{}", operation.method)
                }
                CipStatePolicy::ClearComputed => {
                    state.assert_computed_cleared();
                    assert_eq!(state.atom_codes.first(), source_state.atom_codes.first());
                    assert_eq!(state.bond_codes.first(), source_state.bond_codes.first());
                }
                CipStatePolicy::Assign => assert_eq!(state.computed.as_deref(), Some("1")),
                CipStatePolicy::TautomerSourceTransition => {
                    assert_eq!(operation.method, "enumerate_tautomers_with_options");
                }
            }
        }
        assert_eq!(
            source,
            lifecycle_seed(),
            "{} mutated its source",
            operation.method
        );
    }

    assert_eq!(
        executed,
        MOLECULE_OPS
            .iter()
            .filter(|operation| operation.domain == OperationDomain::Topology)
            .count()
    );
}

#[test]
fn cip_state_lifecycle_clone_and_binary_roundtrip_preserve_complete_state() {
    let labeled = lifecycle_seed()
        .with_cip_labels()
        .expect("modern CIP assignment must succeed");
    let expected = CipState::from_molecule(&labeled);

    assert!(expected.computed_is_computed);
    assert!(
        expected
            .atom_code_is_computed
            .iter()
            .all(|computed| !computed)
    );
    assert!(
        expected
            .atom_neighbor_order_is_computed
            .iter()
            .zip(&expected.atom_neighbor_orders)
            .all(|(computed, value)| value.is_none() || *computed)
    );
    assert!(
        expected
            .atom_rank_is_computed
            .iter()
            .zip(&expected.atom_ranks)
            .all(|(computed, value)| value.is_none() || *computed)
    );
    assert!(
        expected
            .bond_code_is_computed
            .iter()
            .all(|computed| !computed)
    );
    assert!(
        expected
            .bond_neighbor_order_is_computed
            .iter()
            .zip(&expected.bond_neighbor_orders)
            .all(|(computed, value)| value.is_none() || *computed)
    );

    assert_eq!(CipState::from_molecule(&labeled.clone()), expected);
    let encoded = mol_to_binary(&labeled).expect("CIP molecule must serialize");
    let restored = mol_from_binary(&encoded).expect("CIP molecule must deserialize");
    assert_eq!(CipState::from_molecule(&restored), expected);
}

#[test]
fn cip_state_lifecycle_registry_has_one_assignment_owner() {
    let assignment_owners = MOLECULE_OPS
        .iter()
        .filter(|operation| operation.cip_state == CipStatePolicy::Assign)
        .map(|operation| operation.method)
        .collect::<Vec<_>>();
    assert_eq!(assignment_owners, ["with_cip_labels_with_options"]);
}

#[test]
fn computed_clearing_uses_membership_instead_of_cip_property_names() {
    let mut builder = MoleculeBuilder::new();
    let carbon = builder.add_atom(
        AtomSpec::new(Element::C)
            .with_prop("_CIPCode", "user-code")
            .with_prop("_CIPNeighborOrder", "user-neighbors")
            .with_prop("_CIPRank", "user-rank"),
    );
    let oxygen = builder.add_atom(AtomSpec::new(Element::O));
    builder
        .add_bond(
            BondSpec::new(carbon, oxygen, BondOrder::Single)
                .with_prop("_CIPCode", "user-bond-code")
                .with_prop("_CIPNeighborOrder", "user-bond-neighbors"),
        )
        .expect("membership regression bond must build");
    let molecule = builder
        .build()
        .expect("membership regression molecule must build");

    let sanitized = molecule
        .sanitize()
        .expect("non-computed CIP-named properties must not violate the clear-computed contract");

    assert_eq!(sanitized.atoms()[0].prop("_CIPCode"), Some("user-code"));
    assert_eq!(
        sanitized.atoms()[0].prop("_CIPNeighborOrder"),
        Some("user-neighbors")
    );
    assert_eq!(sanitized.atoms()[0].prop("_CIPRank"), Some("user-rank"));
    assert_eq!(
        sanitized.bonds()[0].prop("_CIPCode"),
        Some("user-bond-code")
    );
    assert_eq!(
        sanitized.bonds()[0].prop("_CIPNeighborOrder"),
        Some("user-bond-neighbors")
    );
}
