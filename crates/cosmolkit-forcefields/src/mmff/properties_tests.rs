//! All retained original detached properties/parameter/hydrogen unit cases.
use super::*;
use cosmolkit_model::{
    AtomId, AtomSpec, BondId, BondOrder, BondSpec, Conformer3D, CoordinateBlock, Element,
};
// Test-only detached fixture storage. No live TestInput or domain owner is
// reimplemented. Every old literal atom/bond/conformer row is retained.
#[derive(Clone, Debug, Default)]
struct TestInput {
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
}
impl std::ops::Deref for TestInput {
    type Target = TopologyBlock;
    fn deref(&self) -> &TopologyBlock {
        &self.topology
    }
}
impl TestInput {
    fn new() -> Self {
        Self::default()
    }
    fn atoms(&self) -> &[cosmolkit_model::Atom] {
        &self.topology.atoms
    }
    fn num_atoms(&self) -> usize {
        self.topology.atoms.len()
    }
    fn num_bonds(&self) -> usize {
        self.topology.bonds.len()
    }
}
struct TestFixtureBuilder(TestInput);
impl TestFixtureBuilder {
    fn new() -> Self {
        Self(TestInput::default())
    }
    fn add_atom(&mut self, spec: AtomSpec) -> AtomId {
        {
            let id = AtomId::new(self.0.topology.atoms.len());
            self.0
                .topology
                .atoms
                .push(cosmolkit_model::Atom::from_spec(id, spec));
            self.0.topology.adjacency = cosmolkit_model::AdjacencyList::from_topology(
                self.0.topology.atoms.len(),
                &self.0.topology.bonds,
            );
            id
        }
    }
    fn add_bond(&mut self, spec: BondSpec) -> Result<BondId, cosmolkit_model::TopologyEditError> {
        {
            let id = BondId::new(self.0.topology.bonds.len());
            self.0
                .topology
                .bonds
                .push(cosmolkit_model::Bond::from_spec(id, spec));
            self.0.topology.adjacency = cosmolkit_model::AdjacencyList::from_topology(
                self.0.topology.atoms.len(),
                &self.0.topology.bonds,
            );
            Ok(id)
        }
    }
    fn add_3d_conformer(&mut self, rows: Vec<[f64; 3]>) -> Result<(), String> {
        let id = self.0.coordinates.conformers_3d.len();
        self.add_conformer(Conformer3D::new(id, rows, true))
    }
    fn add_conformer(&mut self, row: Conformer3D) -> Result<(), String> {
        row.validate_for_atom_count(self.0.num_atoms())
            .map_err(|e| e.to_string())?;
        self.0.coordinates.conformers_3d.push(row);
        Ok(())
    }
    fn build(self) -> Result<TestInput, String> {
        self.0.topology.validate().map_err(|e| e.to_string())?;
        Ok(self.0)
    }
}
fn mmff_props_for_atom_types(atom_types: &[u8]) -> MmffMolProperties {
    mmff_props_for_molecule_and_atom_types(TestInput::new(), atom_types)
}

fn mmff_props_for_molecule_and_atom_types(
    molecule: TestInput,
    atom_types: &[u8],
) -> MmffMolProperties {
    mmff_props_for_atom_properties(
        molecule,
        &atom_types
            .iter()
            .copied()
            .map(|atom_type| MmffAtomProperties {
                atom_type,
                formal_charge: 0.0,
                partial_charge: 0.0,
            })
            .collect::<Vec<_>>(),
    )
}

fn mmff_props_for_atom_properties(
    molecule: TestInput,
    atom_properties: &[MmffAtomProperties],
) -> MmffMolProperties {
    let num_bonds = molecule.num_bonds();
    MmffMolProperties {
        topology: molecule.topology,
        valid: true,
        variant: MmffVariant::Mmff94,
        bond_term: true,
        angle_term: true,
        stretch_bend_term: true,
        oop_term: true,
        torsion_term: true,
        vdw_term: true,
        ele_term: true,
        dielectric_constant: 1.0,
        dielectric_model: MMFF_DIELECTRIC_CONSTANT,
        verbosity: 0,
        atom_properties: atom_properties.to_vec(),
        aromatic_ring_count: 0,
        acquired_rings: None,
    }
}

fn two_atom_molecule(first: Element, second: Element, order: BondOrder) -> TestInput {
    two_atom_molecule_with_aromaticity(first, second, order, false)
}

fn two_atom_molecule_with_aromaticity(
    first: Element,
    second: Element,
    order: BondOrder,
    aromatic: bool,
) -> TestInput {
    let mut builder = TestFixtureBuilder::new();
    let first = builder.add_atom(AtomSpec::new(first));
    let second = builder.add_atom(AtomSpec::new(second));
    builder
        .add_bond(BondSpec::new(first, second, order).with_aromatic(aromatic))
        .expect("test molecule bond endpoints are valid");
    builder.build().expect("test molecule is valid")
}

fn hydrogen_neighbor_type(neighbor_element: Element, neighbor_atom_type: u8) -> u8 {
    let molecule = two_atom_molecule(Element::H, neighbor_element, BondOrder::Single);
    let atom_properties = [
        MmffAtomProperties::default(),
        MmffAtomProperties {
            atom_type: neighbor_atom_type,
            formal_charge: 0.0,
            partial_charge: 0.0,
        },
    ];
    set_mmff_hydrogen_type(&molecule, &atom_properties, &molecule.atoms()[0])
        .expect("hydrogen atom typing succeeds")
        .atom_type
}

fn oxygen_type_six_environment(
    ipso_element: Element,
    terminal: Option<(Element, BondOrder)>,
) -> TestInput {
    let mut builder = TestFixtureBuilder::new();
    let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
    let oxygen = builder.add_atom(AtomSpec::new(Element::O));
    let ipso = builder.add_atom(AtomSpec::new(ipso_element));
    builder
        .add_bond(BondSpec::new(hydrogen, oxygen, BondOrder::Single))
        .expect("test H-O bond endpoints are valid");
    builder
        .add_bond(BondSpec::new(oxygen, ipso, BondOrder::Single))
        .expect("test O-ipso bond endpoints are valid");
    if let Some((terminal_element, order)) = terminal {
        let terminal = builder.add_atom(AtomSpec::new(terminal_element));
        builder
            .add_bond(BondSpec::new(ipso, terminal, order))
            .expect("test ipso-terminal bond endpoints are valid");
    }
    builder.build().expect("test molecule is valid")
}

fn hydrogen_type_for_oxygen_type_six_environment(molecule: &TestInput) -> u8 {
    let mut atom_properties = vec![MmffAtomProperties::default(); molecule.num_atoms()];
    atom_properties[1].atom_type = 6;
    set_mmff_hydrogen_type(molecule, &atom_properties, &molecule.atoms()[0])
        .expect("oxygen-bound hydrogen atom typing succeeds")
        .atom_type
}

fn three_atom_angle_molecule(first: Element, second: Element, third: Element) -> TestInput {
    let mut builder = TestFixtureBuilder::new();
    let first = builder.add_atom(AtomSpec::new(first));
    let second = builder.add_atom(AtomSpec::new(second));
    let third = builder.add_atom(AtomSpec::new(third));
    builder
        .add_bond(BondSpec::new(first, second, BondOrder::Single))
        .expect("test molecule bond endpoints are valid");
    builder
        .add_bond(BondSpec::new(second, third, BondOrder::Single))
        .expect("test molecule bond endpoints are valid");
    builder.build().expect("test molecule is valid")
}

fn four_atom_torsion_molecule(
    first: Element,
    second: Element,
    third: Element,
    fourth: Element,
) -> TestInput {
    let mut builder = TestFixtureBuilder::new();
    let first = builder.add_atom(AtomSpec::new(first));
    let second = builder.add_atom(AtomSpec::new(second));
    let third = builder.add_atom(AtomSpec::new(third));
    let fourth = builder.add_atom(AtomSpec::new(fourth));
    for (begin, end) in [(first, second), (second, third), (third, fourth)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    builder.build().expect("test molecule is valid")
}

fn four_atom_oop_molecule(
    first: Element,
    central: Element,
    third: Element,
    fourth: Element,
) -> TestInput {
    let mut builder = TestFixtureBuilder::new();
    let first = builder.add_atom(AtomSpec::new(first));
    let central = builder.add_atom(AtomSpec::new(central));
    let third = builder.add_atom(AtomSpec::new(third));
    let fourth = builder.add_atom(AtomSpec::new(fourth));
    for (begin, end) in [(first, central), (central, third), (central, fourth)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    builder.build().expect("test molecule is valid")
}

fn triangle_molecule() -> TestInput {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    for (begin, end) in [(a0, a1), (a1, a2), (a2, a0)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    builder.build().expect("test molecule is valid")
}

fn pentagon_molecule() -> TestInput {
    let mut builder = TestFixtureBuilder::new();
    let atoms = (0..5)
        .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
        .collect::<Vec<_>>();
    for idx in 0..5 {
        builder
            .add_bond(BondSpec::new(
                atoms[idx],
                atoms[(idx + 1) % 5],
                BondOrder::Single,
            ))
            .expect("test molecule bond endpoints are valid");
    }
    builder.build().expect("test molecule is valid")
}

fn square_with_diagonal_molecule() -> TestInput {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    for (begin, end) in [(a0, a1), (a1, a2), (a2, a3), (a3, a0), (a0, a2)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    builder.build().expect("test molecule is valid")
}

fn square_molecule() -> TestInput {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    for (begin, end) in [(a0, a1), (a1, a2), (a2, a3), (a3, a0)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    builder.build().expect("test molecule is valid")
}

#[test]
fn mmff_mol_properties_atom_properties_default_matches_rdkit_header() {
    let props = MmffAtomProperties::default();

    assert_eq!(props.atom_type, 0);
    assert_eq!(props.formal_charge, 0.0);
    assert_eq!(props.partial_charge, 0.0);
}

#[test]
fn mmff_mol_properties_constructor_variant_only_accepts_exact_lowercase_s_like_rdkit() {
    let molecule = TestInput::new();

    let mmff94s = MmffMolProperties::new(&molecule, false, "MMFF94s", 0).unwrap();
    let mmff94_upper_s = MmffMolProperties::new(&molecule, false, "MMFF94S", 0).unwrap();

    assert_eq!(mmff94s.mmff_variant(), MmffVariant::Mmff94s);
    assert_eq!(mmff94_upper_s.mmff_variant(), MmffVariant::Mmff94);
}

#[test]
fn mmff_hydrogen_type_covers_source_neighbor_element_matrix() {
    for (name, element, neighbor_type, expected) in [
        ("carbon", Element::C, 1, 5),
        ("silicon", Element::SI, 19, 5),
        ("phosphorus", Element::P, 25, 71),
        ("sulfur", Element::S, 15, 71),
        ("unmodeled fluorine neighbor", Element::F, 11, 0),
    ] {
        assert_eq!(
            hydrogen_neighbor_type(element, neighbor_type),
            expected,
            "{name} hydrogen type mismatch"
        );
    }
}

#[test]
fn mmff_hydrogen_type_covers_every_source_nitrogen_type_case() {
    for neighbor_type in [8, 39, 62, 67, 68] {
        assert_eq!(
            hydrogen_neighbor_type(Element::N, neighbor_type),
            23,
            "nitrogen type {neighbor_type} must produce HNR/HPYL/HNOX type 23"
        );
    }
    for neighbor_type in [34, 54, 55, 56, 58, 81] {
        assert_eq!(
            hydrogen_neighbor_type(Element::N, neighbor_type),
            36,
            "nitrogen type {neighbor_type} must produce cationic hydrogen type 36"
        );
    }
    assert_eq!(hydrogen_neighbor_type(Element::N, 9), 27);
    for neighbor_type in [0, 10, 40, 80, u8::MAX] {
        assert_eq!(
            hydrogen_neighbor_type(Element::N, neighbor_type),
            28,
            "default nitrogen type {neighbor_type} must produce type 28"
        );
    }
}

#[test]
fn mmff_hydrogen_type_covers_source_oxygen_type_dispatch() {
    for (neighbor_type, expected) in [(49, 50), (51, 52), (70, 31)] {
        assert_eq!(
            hydrogen_neighbor_type(Element::O, neighbor_type),
            expected,
            "oxygen type {neighbor_type} dispatch mismatch"
        );
    }
    for neighbor_type in [0, 1, 5, 7, 69, 71, u8::MAX] {
        assert_eq!(
            hydrogen_neighbor_type(Element::O, neighbor_type),
            21,
            "default oxygen type {neighbor_type} must produce generic hydroxyl type 21"
        );
    }
}

#[test]
fn mmff_hydrogen_type_covers_oxygen_type_six_environment_matrix() {
    struct Case {
        name: &'static str,
        ipso: Element,
        terminal: Option<(Element, BondOrder)>,
        expected: u8,
    }

    let cases = [
        Case {
            name: "HOCO",
            ipso: Element::C,
            terminal: Some((Element::O, BondOrder::Double)),
            expected: 24,
        },
        Case {
            name: "HOP",
            ipso: Element::P,
            terminal: None,
            expected: 24,
        },
        Case {
            name: "HOCC double",
            ipso: Element::C,
            terminal: Some((Element::C, BondOrder::Double)),
            expected: 29,
        },
        Case {
            name: "HOCN double",
            ipso: Element::C,
            terminal: Some((Element::N, BondOrder::Double)),
            expected: 29,
        },
        Case {
            name: "HOCC aromatic",
            ipso: Element::C,
            terminal: Some((Element::C, BondOrder::Aromatic)),
            expected: 29,
        },
        Case {
            name: "HOCN aromatic",
            ipso: Element::C,
            terminal: Some((Element::N, BondOrder::Aromatic)),
            expected: 29,
        },
        Case {
            name: "HOS",
            ipso: Element::S,
            terminal: None,
            expected: 33,
        },
        Case {
            name: "generic alcohol",
            ipso: Element::C,
            terminal: Some((Element::C, BondOrder::Single)),
            expected: 21,
        },
    ];

    for case in cases {
        let molecule = oxygen_type_six_environment(case.ipso, case.terminal);
        assert_eq!(
            hydrogen_type_for_oxygen_type_six_environment(&molecule),
            case.expected,
            "{} environment mismatch",
            case.name
        );
    }
}

#[test]
fn mmff_hydrogen_type_prefers_hoco_over_hocc_like_source_switch_break() {
    let mut builder = TestFixtureBuilder::new();
    let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
    let oxygen = builder.add_atom(AtomSpec::new(Element::O));
    let carbonyl_carbon = builder.add_atom(AtomSpec::new(Element::C));
    let carbonyl_oxygen = builder.add_atom(AtomSpec::new(Element::O));
    let alkene_carbon = builder.add_atom(AtomSpec::new(Element::C));
    for (begin, end, order) in [
        (hydrogen, oxygen, BondOrder::Single),
        (oxygen, carbonyl_carbon, BondOrder::Single),
        (carbonyl_carbon, carbonyl_oxygen, BondOrder::Double),
        (carbonyl_carbon, alkene_carbon, BondOrder::Double),
    ] {
        builder
            .add_bond(BondSpec::new(begin, end, order))
            .expect("test priority molecule bond endpoints are valid");
    }
    let molecule = builder.build().expect("test priority molecule is valid");

    assert_eq!(hydrogen_type_for_oxygen_type_six_environment(&molecule), 24);
}

#[test]
fn mmff_mol_properties_constructor_marks_unbound_hydrogen_invalid() {
    let mut builder = TestFixtureBuilder::new();
    builder.add_atom(AtomSpec::new(Element::H));
    let molecule = builder
        .build()
        .expect("isolated hydrogen molecule is valid");
    let props = MmffMolProperties::new(&molecule, false, "MMFF94", MMFF_VERBOSITY_NONE)
        .expect("isolated hydrogen returns invalid MMFF properties");

    assert!(!props.is_valid());
    assert_eq!(props.get_mmff_atom_type(0).unwrap(), 0);
}

#[test]
fn mmff_mol_properties_get_atom_type_returns_stored_atom_type() {
    let props = mmff_props_for_atom_types(&[1, 37, 62]);

    assert_eq!(props.get_mmff_atom_type(1).unwrap(), 37);
}

#[test]
fn mmff_mol_properties_get_atom_type_accepts_last_valid_index() {
    let props = mmff_props_for_atom_types(&[4, 8, 12]);

    assert_eq!(props.get_mmff_atom_type(2).unwrap(), 12);
}

#[test]
fn mmff_mol_properties_get_atom_type_reports_out_of_range_index() {
    let props = mmff_props_for_atom_types(&[7, 9]);

    let err = props.get_mmff_atom_type(2).unwrap_err();

    match err {
        MmffMolPropertiesError::AtomIndexOutOfRange { atom_index, atoms } => {
            assert_eq!(atom_index, 2);
            assert_eq!(atoms, 2);
        }
        other => panic!("expected atom-index out-of-range error, got {other:?}"),
    }
}

#[test]
fn mmff_mol_properties_get_formal_charge_returns_stored_formal_charge() {
    let props = mmff_props_for_atom_properties(
        TestInput::new(),
        &[
            MmffAtomProperties {
                atom_type: 1,
                formal_charge: 0.0,
                partial_charge: 0.0,
            },
            MmffAtomProperties {
                atom_type: 37,
                formal_charge: -1.25,
                partial_charge: 0.0,
            },
        ],
    );

    assert_eq!(props.get_mmff_formal_charge(1).unwrap(), -1.25);
}

#[test]
fn mmff_mol_properties_get_formal_charge_accepts_last_valid_index() {
    let props = mmff_props_for_atom_properties(
        TestInput::new(),
        &[
            MmffAtomProperties {
                atom_type: 4,
                formal_charge: 0.5,
                partial_charge: 0.0,
            },
            MmffAtomProperties {
                atom_type: 8,
                formal_charge: 1.75,
                partial_charge: 0.0,
            },
        ],
    );

    assert_eq!(props.get_mmff_formal_charge(1).unwrap(), 1.75);
}

#[test]
fn mmff_mol_properties_get_formal_charge_reports_out_of_range_index() {
    let props = mmff_props_for_atom_properties(
        TestInput::new(),
        &[MmffAtomProperties {
            atom_type: 7,
            formal_charge: -0.5,
            partial_charge: 0.0,
        }],
    );

    let err = props.get_mmff_formal_charge(1).unwrap_err();

    match err {
        MmffMolPropertiesError::AtomIndexOutOfRange { atom_index, atoms } => {
            assert_eq!(atom_index, 1);
            assert_eq!(atoms, 1);
        }
        other => panic!("expected atom-index out-of-range error, got {other:?}"),
    }
}

#[test]
fn mmff_mol_properties_get_partial_charge_returns_stored_partial_charge() {
    let props = mmff_props_for_atom_properties(
        TestInput::new(),
        &[
            MmffAtomProperties {
                atom_type: 1,
                formal_charge: 0.0,
                partial_charge: 0.125,
            },
            MmffAtomProperties {
                atom_type: 37,
                formal_charge: 0.0,
                partial_charge: -0.625,
            },
        ],
    );

    assert_eq!(props.get_mmff_partial_charge(1).unwrap(), -0.625);
}

#[test]
fn mmff_mol_properties_get_partial_charge_accepts_last_valid_index() {
    let props = mmff_props_for_atom_properties(
        TestInput::new(),
        &[
            MmffAtomProperties {
                atom_type: 4,
                formal_charge: 0.0,
                partial_charge: -0.25,
            },
            MmffAtomProperties {
                atom_type: 8,
                formal_charge: 0.0,
                partial_charge: 0.875,
            },
        ],
    );

    assert_eq!(props.get_mmff_partial_charge(1).unwrap(), 0.875);
}

#[test]
fn mmff_mol_properties_get_partial_charge_reports_out_of_range_index() {
    let props = mmff_props_for_atom_properties(
        TestInput::new(),
        &[MmffAtomProperties {
            atom_type: 7,
            formal_charge: 0.0,
            partial_charge: -0.125,
        }],
    );

    let err = props.get_mmff_partial_charge(1).unwrap_err();

    match err {
        MmffMolPropertiesError::AtomIndexOutOfRange { atom_index, atoms } => {
            assert_eq!(atom_index, 1);
            assert_eq!(atoms, 1);
        }
        other => panic!("expected atom-index out-of-range error, got {other:?}"),
    }
}

#[test]
fn mmff_mol_properties_get_bond_stretch_params_returns_tabulated_params() {
    let molecule = two_atom_molecule(Element::C, Element::H, BondOrder::Single);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 5]);

    let params = props
        .get_mmff_bond_stretch_params(0, 1)
        .expect("default MMFF bond table parses")
        .expect("C-H atom types have tabulated bond-stretch params");

    assert_eq!(params.0, 0);
    assert_eq!(params.1.kb, 4.766);
    assert_eq!(params.1.r0, 1.093);
}

#[test]
fn mmff_mol_properties_get_bond_stretch_params_accepts_reversed_atom_order() {
    let molecule = two_atom_molecule(Element::C, Element::H, BondOrder::Single);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 5]);

    let params = props
        .get_mmff_bond_stretch_params(1, 0)
        .expect("default MMFF bond table parses")
        .expect("reverse lookup finds the same undirected bond");

    assert_eq!(params.0, 0);
    assert_eq!(params.1.kb, 4.766);
    assert_eq!(params.1.r0, 1.093);
}

#[test]
fn mmff_mol_properties_get_bond_stretch_params_returns_none_when_invalid() {
    let molecule = two_atom_molecule(Element::C, Element::H, BondOrder::Single);
    let mut props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 5]);
    props.valid = false;

    let params = props
        .get_mmff_bond_stretch_params(0, 1)
        .expect("invalid properties return RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_bond_stretch_params_returns_none_without_bond() {
    let mut builder = TestFixtureBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    builder.add_atom(AtomSpec::new(Element::H));
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 5]);

    let params = props
        .get_mmff_bond_stretch_params(0, 1)
        .expect("unbonded atoms return RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_bond_stretch_params_reports_atom_index_out_of_range() {
    let molecule = two_atom_molecule(Element::C, Element::H, BondOrder::Single);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 5]);

    let err = props.get_mmff_bond_stretch_params(0, 2).unwrap_err();

    match err {
        MmffMolPropertiesError::AtomIndexOutOfRange { atom_index, atoms } => {
            assert_eq!(atom_index, 2);
            assert_eq!(atoms, 2);
        }
        other => panic!("expected atom-index out-of-range error, got {other:?}"),
    }
}

#[test]
fn mmff_mol_properties_get_bond_stretch_params_matches_rdkit_c_o_empirical_rule() {
    let molecule = two_atom_molecule(Element::C, Element::O, BondOrder::Single);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[60, 6]);

    let params = props
        .get_mmff_bond_stretch_params(0, 1)
        .expect("C-O empirical lookup succeeds")
        .expect("bond has empirical parameters");

    assert_eq!(params.0, 0);
    assert!((params.1.r0 - 1.405).abs() < 1.0e-12);
    assert!((params.1.kb - 5.129115902527102).abs() < 1.0e-12);
}

#[test]
fn mmff_mol_properties_get_bond_stretch_params_matches_rdkit_c_c_empirical_rule() {
    let molecule = two_atom_molecule(Element::C, Element::C, BondOrder::Single);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[63, 22]);

    let params = props
        .get_mmff_bond_stretch_params(0, 1)
        .expect("C-C empirical lookup succeeds")
        .expect("bond has empirical parameters");

    assert_eq!(params.0, 0);
    assert!((params.1.r0 - 1.54).abs() < 1.0e-12);
    assert!((params.1.kb - 3.403846905179782).abs() < 1.0e-12);
}

#[test]
fn mmff_mol_properties_empirical_bond_rule_accepts_reversed_atom_order() {
    let molecule = two_atom_molecule(Element::C, Element::O, BondOrder::Single);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[60, 6]);

    let forward = props.get_mmff_bond_stretch_params(0, 1).unwrap().unwrap();
    let reversed = props.get_mmff_bond_stretch_params(1, 0).unwrap().unwrap();

    assert_eq!(forward, reversed);
}

#[test]
fn mmff_mol_properties_empirical_bond_rule_uses_herschbach_laurie_without_bndk() {
    let molecule = two_atom_molecule(Element::LI, Element::LI, BondOrder::Single);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1]);
    let mmff_prop = default_mmff_prop().unwrap();

    let params = props
        .get_mmff_bond_stretch_empirical_rule_params(0, 1, mmff_prop)
        .expect("Li-Li uses the Herschbach-Laurie fallback");

    assert!((params.r0 - 2.68).abs() < 1.0e-12);
    assert!((params.kb - 0.5904545052481254).abs() < 1.0e-12);
}

#[test]
fn mmff_mol_properties_get_angle_bend_params_returns_tabulated_params() {
    let molecule = three_atom_angle_molecule(Element::H, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[5, 1, 1]);

    let params = props
        .get_mmff_angle_bend_params(0, 1, 2)
        .expect("default MMFF angle table parses")
        .expect("H-C-C atom types have tabulated angle-bend params");

    assert_eq!(params.0, 0);
    assert_eq!(params.1.ka, 0.636);
    assert_eq!(params.1.theta0, 110.549);
}

#[test]
fn mmff_mol_properties_get_angle_bend_params_accepts_reversed_endpoint_order() {
    let molecule = three_atom_angle_molecule(Element::H, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[5, 1, 1]);

    let params = props
        .get_mmff_angle_bend_params(2, 1, 0)
        .expect("default MMFF angle table parses")
        .expect("reverse endpoint lookup finds the same angle params");

    assert_eq!(params.0, 0);
    assert_eq!(params.1.ka, 0.636);
    assert_eq!(params.1.theta0, 110.549);
}

#[test]
fn mmff_mol_properties_get_angle_bend_params_handles_three_membered_ring_angle_type() {
    let molecule = triangle_molecule();
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[8, 8, 8]);

    let params = props
        .get_mmff_angle_bend_params(0, 1, 2)
        .expect("default MMFF angle table parses")
        .expect("N-N-N three-membered ring angle has tabulated params");

    assert_eq!(params.0, 3);
    assert_eq!(params.1.ka, 0.230);
    assert_eq!(params.1.theta0, 60.000);
}

#[test]
fn mmff_mol_properties_get_angle_bend_params_handles_four_membered_ring_angle_type() {
    let molecule = square_molecule();
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[3, 3, 3, 3]);

    let params = props
        .get_mmff_angle_bend_params(0, 1, 2)
        .expect("default MMFF angle table parses")
        .expect("C-C-C four-membered ring angle has tabulated params");

    assert_eq!(params.0, 8);
    assert_eq!(params.1.ka, 1.280);
    assert_eq!(params.1.theta0, 89.965);
}

#[test]
fn mmff_mol_properties_get_angle_bend_params_returns_none_when_invalid() {
    let molecule = three_atom_angle_molecule(Element::H, Element::C, Element::C);
    let mut props = mmff_props_for_molecule_and_atom_types(molecule, &[5, 1, 1]);
    props.valid = false;

    let params = props
        .get_mmff_angle_bend_params(0, 1, 2)
        .expect("invalid properties return RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_angle_bend_params_returns_none_without_first_bond() {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::H));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(a1, a2, BondOrder::Single))
        .expect("test molecule bond endpoints are valid");
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[5, 1, 1]);

    let params = props
        .get_mmff_angle_bend_params(a0.index(), a1.index(), a2.index())
        .expect("missing first bond returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_angle_bend_params_returns_none_without_second_bond() {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::H));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
        .expect("test molecule bond endpoints are valid");
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[5, 1, 1]);

    let params = props
        .get_mmff_angle_bend_params(a0.index(), a1.index(), a2.index())
        .expect("missing second bond returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_angle_bend_params_reports_atom_index_out_of_range() {
    let molecule = three_atom_angle_molecule(Element::H, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[5, 1]);

    let err = props.get_mmff_angle_bend_params(0, 1, 2).unwrap_err();

    match err {
        MmffMolPropertiesError::AtomIndexOutOfRange { atom_index, atoms } => {
            assert_eq!(atom_index, 2);
            assert_eq!(atoms, 2);
        }
        other => panic!("expected atom-index out-of-range error, got {other:?}"),
    }
}

#[test]
fn mmff_mol_properties_get_angle_bend_params_matches_rdkit_old_row_empirical_ka() {
    let molecule = three_atom_angle_molecule(Element::F, Element::C, Element::CL);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[11, 1, 12]);

    let forward = props
        .get_mmff_angle_bend_params(0, 1, 2)
        .expect("default MMFF angle table parses")
        .expect("F-C-Cl angle receives empirical force constant");
    let reversed = props
        .get_mmff_angle_bend_params(2, 1, 0)
        .expect("default MMFF angle table parses")
        .expect("reversed F-C-Cl angle receives empirical force constant");

    assert_eq!(forward.0, 0);
    assert!((forward.1.ka - 1.2566039721725888).abs() < 1.0e-12);
    assert_eq!(forward.1.theta0, 108.9);
    assert_eq!(forward, reversed);
}

#[test]
fn mmff_angle_empirical_rule_covers_rest_value_and_small_ring_branches() {
    let bond = |r0| MmffBond { kb: 0.0, r0 };
    let prop = |atno, crd, val, mltb, linh| MmffProp {
        atno,
        crd,
        val,
        pilp: 0,
        mltb,
        arom: 0,
        linh,
        sbmb: 0,
    };
    let cases = [
        (
            three_atom_angle_molecule(Element::F, Element::C, Element::CL),
            prop(6, 4, 4, 0, 0),
            bond(1.36),
            bond(1.773),
            109.45,
            1.244006518144849,
        ),
        (
            three_atom_angle_molecule(Element::H, Element::O, Element::C),
            prop(8, 2, 2, 0, 0),
            bond(0.97),
            bond(1.43),
            105.0,
            0.938397748236572,
        ),
        (
            three_atom_angle_molecule(Element::C, Element::C, Element::C),
            prop(6, 2, 4, 0, 1),
            bond(1.54),
            bond(1.54),
            180.0,
            0.36380963203127237,
        ),
        (
            three_atom_angle_molecule(Element::H, Element::N, Element::C),
            prop(7, 3, 3, 0, 0),
            bond(1.01),
            bond(1.47),
            107.0,
            0.731386150508525,
        ),
        (
            three_atom_angle_molecule(Element::C, Element::P, Element::C),
            prop(15, 3, 3, 0, 0),
            bond(1.84),
            bond(1.84),
            92.0,
            1.2252479685179425,
        ),
        (
            triangle_molecule(),
            prop(6, 4, 4, 0, 0),
            bond(1.54),
            bond(1.54),
            60.0,
            0.16371433441407263,
        ),
        (
            square_molecule(),
            prop(6, 4, 4, 0, 0),
            bond(1.54),
            bond(1.54),
            90.0,
            1.2369527489063263,
        ),
    ];

    for (molecule, central_prop, bond_1, bond_2, expected_theta0, expected_ka) in cases {
        let props = mmff_props_for_molecule_and_atom_types(
            molecule,
            &[central_prop.atno, central_prop.atno, central_prop.atno],
        );
        let params = props.get_mmff_angle_bend_empirical_rule_params(
            None,
            &central_prop,
            &bond_1,
            &bond_2,
            0,
            1,
            2,
        );

        assert_eq!(params.theta0, expected_theta0);
        assert!((params.ka - expected_ka).abs() < 1.0e-12);
    }

    let molecule = three_atom_angle_molecule(Element::F, Element::C, Element::CL);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[11, 1, 12]);
    let params = props.get_mmff_angle_bend_empirical_rule_params(
        Some(&MmffAngle {
            ka: 0.0,
            theta0: 108.9,
        }),
        &prop(6, 4, 4, 0, 0),
        &bond(1.36),
        &bond(1.773),
        0,
        1,
        2,
    );
    assert_eq!(params.theta0, 108.9);
    assert!((params.ka - 1.2566039721725888).abs() < 1.0e-12);
}

#[test]
fn mmff_mol_properties_get_stretch_bend_params_returns_tabulated_params() {
    let molecule = three_atom_angle_molecule(Element::C, Element::C, Element::H);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 5]);

    let params = props
        .get_mmff_stretch_bend_params(0, 1, 2)
        .expect("default MMFF stretch-bend tables parse")
        .expect("C-C-H atom types have tabulated stretch-bend params");

    assert_eq!(params.0, 0);
    assert_eq!(params.1.kba_ijk, 0.227);
    assert_eq!(params.1.kba_kji, 0.070);
    assert_eq!(params.2[0].kb, 4.258);
    assert_eq!(params.2[0].r0, 1.508);
    assert_eq!(params.2[1].kb, 4.766);
    assert_eq!(params.2[1].r0, 1.093);
    assert_eq!(params.3.ka, 0.636);
    assert_eq!(params.3.theta0, 110.549);
}

#[test]
fn mmff_mol_properties_get_stretch_bend_params_swaps_terminal_params_like_rdkit() {
    let molecule = three_atom_angle_molecule(Element::H, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[5, 1, 1]);

    let params = props
        .get_mmff_stretch_bend_params(0, 1, 2)
        .expect("default MMFF stretch-bend tables parse")
        .expect("H-C-C atom types have tabulated stretch-bend params");

    assert_eq!(params.0, 0);
    assert_eq!(params.1.kba_ijk, 0.070);
    assert_eq!(params.1.kba_kji, 0.227);
}

#[test]
fn mmff_mol_properties_get_stretch_bend_params_uses_default_fallback_params() {
    let molecule = three_atom_angle_molecule(Element::C, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 4]);

    let params = props
        .get_mmff_stretch_bend_params(0, 1, 2)
        .expect("default MMFF stretch-bend tables parse")
        .expect("missing explicit stretch-bend row uses default row params");

    assert_eq!(params.0, 0);
    assert_eq!(params.1.kba_ijk, 0.30);
    assert_eq!(params.1.kba_kji, 0.30);
    assert_eq!(params.3.ka, 1.006);
    assert_eq!(params.3.theta0, 110.265);
}

#[test]
fn mmff_mol_properties_get_stretch_bend_params_returns_none_when_invalid() {
    let molecule = three_atom_angle_molecule(Element::C, Element::C, Element::H);
    let mut props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 5]);
    props.valid = false;

    let params = props
        .get_mmff_stretch_bend_params(0, 1, 2)
        .expect("invalid properties return RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_stretch_bend_params_returns_none_without_first_bond() {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::H));
    builder
        .add_bond(BondSpec::new(a1, a2, BondOrder::Single))
        .expect("test molecule bond endpoints are valid");
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 5]);

    let params = props
        .get_mmff_stretch_bend_params(a0.index(), a1.index(), a2.index())
        .expect("missing first bond returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_stretch_bend_params_returns_none_without_second_bond() {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::H));
    builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
        .expect("test molecule bond endpoints are valid");
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 5]);

    let params = props
        .get_mmff_stretch_bend_params(a0.index(), a1.index(), a2.index())
        .expect("missing second bond returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_stretch_bend_params_returns_none_for_linear_central_atom() {
    let molecule = three_atom_angle_molecule(Element::C, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 4, 1]);

    let params = props
        .get_mmff_stretch_bend_params(0, 1, 2)
        .expect("linear central atom returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_stretch_bend_params_reports_atom_index_out_of_range() {
    let molecule = three_atom_angle_molecule(Element::C, Element::C, Element::H);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1]);

    let err = props.get_mmff_stretch_bend_params(0, 1, 2).unwrap_err();

    match err {
        MmffMolPropertiesError::AtomIndexOutOfRange { atom_index, atoms } => {
            assert_eq!(atom_index, 2);
            assert_eq!(atoms, 2);
        }
        other => panic!("expected atom-index out-of-range error, got {other:?}"),
    }
}

#[test]
fn mmff_mol_properties_get_stretch_bend_params_uses_bond_empirical_fallback() {
    let molecule = three_atom_angle_molecule(Element::C, Element::O, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 7, 1]);

    let params = props
        .get_mmff_stretch_bend_params(0, 1, 2)
        .expect("ported bond empirical fallback should succeed")
        .expect("C-O-C stretch-bend parameters should be available");

    assert_eq!(params.0, 0);
    assert_eq!(params.1.kba_ijk, 0.3);
    assert_eq!(params.1.kba_kji, 0.3);
    assert_eq!(params.2[0], params.2[1]);
    assert_eq!(params.2[0].r0, 1.405);
    assert_eq!(params.2[0].kb, 5.129115902527102);
    assert_eq!(params.3.theta0, 120.0);
    assert_eq!(params.3.ka, 1.1806979441809595);
}

#[test]
fn mmff_mol_properties_get_stretch_bend_params_uses_angle_empirical_fallback() {
    let molecule = three_atom_angle_molecule(Element::F, Element::C, Element::CL);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[11, 1, 12]);

    let params = props
        .get_mmff_stretch_bend_params(0, 1, 2)
        .expect("default MMFF stretch-bend tables parse")
        .expect("F-C-Cl stretch-bend receives empirical angle parameters");

    assert!((params.3.ka - 1.2566039721725888).abs() < 1.0e-12);
    assert_eq!(params.3.theta0, 108.9);
    assert_eq!(
        params.2[0],
        MmffBond {
            kb: 6.011,
            r0: 1.36
        }
    );
    assert_eq!(
        params.2[1],
        MmffBond {
            kb: 2.974,
            r0: 1.773
        }
    );
}

#[test]
fn mmff_mol_properties_get_torsion_params_returns_tabulated_params() {
    let molecule = four_atom_torsion_molecule(Element::C, Element::C, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1, 1]);

    let params = props
        .get_mmff_torsion_params(0, 1, 2, 3)
        .expect("default MMFF torsion table parses")
        .expect("C-C-C-C atom types have tabulated torsion params");

    assert_eq!(params.0, 0);
    assert_eq!(params.1.v1, 0.103);
    assert_eq!(params.1.v2, 0.681);
    assert_eq!(params.1.v3, 0.332);
}

#[test]
fn mmff_mol_properties_get_torsion_params_uses_mmff94s_torsion_table() {
    let molecule = four_atom_torsion_molecule(Element::H, Element::C, Element::C, Element::F);
    let mut props = mmff_props_for_molecule_and_atom_types(molecule, &[5, 1, 1, 10]);
    props.variant = MmffVariant::Mmff94s;

    let params = props
        .get_mmff_torsion_params(0, 1, 2, 3)
        .expect("default MMFF94s torsion table parses")
        .expect("H-C-C-F atom types have tabulated MMFF94s torsion params");

    assert_eq!(params.0, 0);
    assert_eq!(params.1.v1, 0.0);
    assert_eq!(params.1.v2, 0.0);
    assert_eq!(params.1.v3, 0.418);
}

#[test]
fn mmff_mol_properties_get_torsion_params_returns_none_for_zero_torsion_params() {
    let molecule = four_atom_torsion_molecule(Element::C, Element::C, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 2, 3]);

    let params = props
        .get_mmff_torsion_params(0, 1, 2, 3)
        .expect("default MMFF torsion table parses");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_torsion_params_returns_none_when_invalid() {
    let molecule = four_atom_torsion_molecule(Element::C, Element::C, Element::C, Element::C);
    let mut props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1, 1]);
    props.valid = false;

    let params = props
        .get_mmff_torsion_params(0, 1, 2, 3)
        .expect("invalid properties return RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_torsion_params_returns_none_without_first_bond() {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    for (begin, end) in [(a1, a2), (a2, a3)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1, 1]);

    let params = props
        .get_mmff_torsion_params(a0.index(), a1.index(), a2.index(), a3.index())
        .expect("missing first bond returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_torsion_params_returns_none_without_middle_bond() {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    for (begin, end) in [(a0, a1), (a2, a3)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1, 1]);

    let params = props
        .get_mmff_torsion_params(a0.index(), a1.index(), a2.index(), a3.index())
        .expect("missing middle bond returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_torsion_params_returns_none_without_last_bond() {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    for (begin, end) in [(a0, a1), (a1, a2)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1, 1]);

    let params = props
        .get_mmff_torsion_params(a0.index(), a1.index(), a2.index(), a3.index())
        .expect("missing last bond returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_torsion_params_reports_atom_index_out_of_range() {
    let molecule = four_atom_torsion_molecule(Element::C, Element::C, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1]);

    let err = props.get_mmff_torsion_params(0, 1, 2, 3).unwrap_err();

    match err {
        MmffMolPropertiesError::AtomIndexOutOfRange { atom_index, atoms } => {
            assert_eq!(atom_index, 3);
            assert_eq!(atoms, 3);
        }
        other => panic!("expected atom-index out-of-range error, got {other:?}"),
    }
}

#[test]
fn mmff_mol_properties_get_torsion_params_dispatches_to_empirical_fallback() {
    let molecule = four_atom_torsion_molecule(Element::C, Element::C, Element::C, Element::MG);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1, 99]);

    // Pinned vector definitions contain95 rows; Params.h dereferences
    // type99 before empirical fallback. Preserve the original synthetic
    // fixture and propagate its structural cause instead of guessing data.
    assert!(matches!(
        props.get_mmff_torsion_params(0, 1, 2, 3),
        Err(MmffMolPropertiesError::Params(MmffParamError::Torsion(
            MmffTorLookupError::MissingDefinition { atom_type: 99 }
        )))
    ));
    // Retain the original empirical formula's numerical coverage at its
    // direct defined helper boundary; central property rows are valid.
    let params = props
        .get_mmff_torsion_empirical_rule_params(1, 2, default_mmff_prop().unwrap())
        .unwrap();
    assert_eq!(params.v1, 0.0);
    assert_eq!(params.v2, 0.0);
    assert!((params.v3 - 2.12 / 9.0).abs() < 1.0e-12);
}

#[test]
fn mmff_torsion_empirical_rules_cover_all_source_branches() {
    struct Case {
        name: &'static str,
        elements: [Element; 2],
        atom_types: [u8; 2],
        order: BondOrder,
        aromatic: bool,
        expected: [f64; 3],
    }

    let cases = [
        Case {
            name: "rule_a_linear",
            elements: [Element::C, Element::C],
            atom_types: [4, 4],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 0.0, 0.0],
        },
        Case {
            name: "rule_b_aromatic",
            elements: [Element::C, Element::N],
            atom_types: [37, 38],
            order: BondOrder::Aromatic,
            aromatic: true,
            expected: [0.0, 3.0, 0.0],
        },
        Case {
            name: "rule_c_double",
            elements: [Element::C, Element::C],
            atom_types: [2, 2],
            order: BondOrder::Double,
            aromatic: false,
            expected: [0.0, 12.0, 0.0],
        },
        Case {
            name: "rule_d_coordination_four",
            elements: [Element::C, Element::C],
            atom_types: [1, 1],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 0.0, 2.12 / 9.0],
        },
        Case {
            name: "rule_e_zero",
            elements: [Element::C, Element::C],
            atom_types: [1, 2],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 0.0, 0.0],
        },
        Case {
            name: "rule_e_nonzero",
            elements: [Element::C, Element::N],
            atom_types: [1, 8],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 0.0, 3.18_f64.sqrt() / 6.0],
        },
        Case {
            name: "rule_f_zero_row_103_boundary",
            elements: [Element::N, Element::C],
            atom_types: [63, 22],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 0.0, 0.0],
        },
        Case {
            name: "rule_f_nonzero",
            elements: [Element::N, Element::C],
            atom_types: [8, 1],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 0.0, 3.18_f64.sqrt() / 6.0],
        },
        Case {
            name: "rule_g_case_1",
            elements: [Element::O, Element::S],
            atom_types: [59, 44],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 0.0, 0.0],
        },
        Case {
            name: "rule_g_case_2_mltb_one",
            elements: [Element::O, Element::C],
            atom_types: [59, 2],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 6.0, 0.0],
        },
        Case {
            name: "rule_g_case_2_same_third_period",
            elements: [Element::S, Element::S],
            atom_types: [15, 17],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 2.25, 0.0],
        },
        Case {
            name: "rule_g_case_2_cross_period",
            elements: [Element::S, Element::C],
            atom_types: [15, 2],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 0.9 * 2.5_f64.sqrt(), 0.0],
        },
        Case {
            name: "rule_g_case_3",
            elements: [Element::C, Element::O],
            atom_types: [2, 59],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 6.0, 0.0],
        },
        Case {
            name: "rule_g_case_4",
            elements: [Element::C, Element::S],
            atom_types: [41, 17],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 2.4 * 2.5_f64.sqrt(), 0.0],
        },
        Case {
            name: "rule_g_case_5",
            elements: [Element::C, Element::C],
            atom_types: [2, 2],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 1.8, 0.0],
        },
        Case {
            name: "rule_h_oxygen_sulfur",
            elements: [Element::O, Element::S],
            atom_types: [6, 15],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, -4.0, 0.0],
        },
        Case {
            name: "rule_h_general",
            elements: [Element::N, Element::N],
            atom_types: [8, 8],
            order: BondOrder::Single,
            aromatic: false,
            expected: [0.0, 0.0, 0.375],
        },
    ];
    let mmff_prop = default_mmff_prop().expect("default MMFFProp parses");

    for case in cases {
        let molecule = two_atom_molecule_with_aromaticity(
            case.elements[0],
            case.elements[1],
            case.order,
            case.aromatic,
        );
        let props = mmff_props_for_molecule_and_atom_types(molecule, &case.atom_types);
        let params = props
            .get_mmff_torsion_empirical_rule_params(0, 1, mmff_prop)
            .unwrap_or_else(|err| panic!("{} empirical lookup failed: {err}", case.name));
        let actual = [params.v1, params.v2, params.v3];

        for (coefficient, (actual, expected)) in ["V1", "V2", "V3"]
            .into_iter()
            .zip(actual.into_iter().zip(case.expected))
        {
            assert!(
                (actual - expected).abs() < 1.0e-12,
                "{} {coefficient} mismatch: actual={actual}, expected={expected}",
                case.name
            );
        }
    }
}

#[test]
fn mmff_mol_properties_get_torsion_type_uses_type_two_single_central_rule() {
    let molecule = four_atom_torsion_molecule(Element::C, Element::C, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[37, 37, 1, 1]);
    let mmff_prop = default_mmff_prop().expect("default MMFFProp parses");

    let torsion_type = props
        .get_mmff_torsion_type(0, 1, 2, 3, mmff_prop)
        .expect("torsion type is computed from tabulated atom properties");

    assert_eq!(torsion_type, (2, 0));
}

#[test]
fn mmff_mol_properties_get_torsion_type_uses_four_membered_ring_type() {
    let molecule = square_molecule();
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1, 1]);
    let mmff_prop = default_mmff_prop().expect("default MMFFProp parses");

    let torsion_type = props
        .get_mmff_torsion_type(0, 1, 2, 3, mmff_prop)
        .expect("four-membered ring torsion type is computed");

    assert_eq!(torsion_type, (4, 0));
}

#[test]
fn mmff_mol_properties_get_torsion_type_keeps_base_type_for_fused_three_ring_guard() {
    let molecule = square_with_diagonal_molecule();
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1, 1]);
    let mmff_prop = default_mmff_prop().expect("default MMFFProp parses");

    let torsion_type = props
        .get_mmff_torsion_type(0, 1, 2, 3, mmff_prop)
        .expect("guarded four-membered torsion type is computed");

    assert_eq!(torsion_type, (0, 0));
}

#[test]
fn mmff_mol_properties_get_torsion_type_uses_five_membered_ring_type_with_carbon_type() {
    let molecule = pentagon_molecule();
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1, 1, 1]);
    let mmff_prop = default_mmff_prop().expect("default MMFFProp parses");

    let torsion_type = props
        .get_mmff_torsion_type(0, 1, 2, 3, mmff_prop)
        .expect("five-membered ring torsion type is computed");

    assert_eq!(torsion_type, (5, 0));
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_returns_tabulated_params() {
    let molecule = four_atom_oop_molecule(Element::C, Element::C, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 2, 1, 2]);

    let params = props
        .get_mmff_oop_bend_params(0, 1, 2, 3)
        .expect("default MMFF OOP table parses")
        .expect("C94 1-2-1-2 atom types have tabulated OOP params");

    assert_eq!(params.koop, 0.030);
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_sorts_outer_atom_types() {
    let molecule = four_atom_oop_molecule(Element::C, Element::C, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[2, 2, 1, 1]);

    let params = props
        .get_mmff_oop_bend_params(0, 1, 2, 3)
        .expect("default MMFF OOP table parses")
        .expect("OOP lookup matches source collection outer-atom sorting");

    assert_eq!(params.koop, 0.030);
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_uses_mmff94_default_row() {
    let molecule = four_atom_oop_molecule(Element::C, Element::N, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 10, 1, 1]);

    let params = props
        .get_mmff_oop_bend_params(0, 1, 2, 3)
        .expect("default MMFF OOP table parses")
        .expect("*-10-*-* default MMFF OOP row is present");

    assert_eq!(params.koop, -0.020);
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_uses_mmff94s_oop_table() {
    let molecule = four_atom_oop_molecule(Element::C, Element::N, Element::C, Element::C);
    let mut props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 10, 1, 1]);
    props.variant = MmffVariant::Mmff94s;

    let params = props
        .get_mmff_oop_bend_params(0, 1, 2, 3)
        .expect("default MMFF94s OOP table parses")
        .expect("*-10-*-* default MMFF94s OOP row is present");

    assert_eq!(params.koop, 0.015);
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_returns_zero_tabulated_params() {
    let molecule = four_atom_oop_molecule(Element::C, Element::O, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 8, 1, 1]);

    let params = props
        .get_mmff_oop_bend_params(0, 1, 2, 3)
        .expect("default MMFF OOP table parses")
        .expect("zero-valued *-8-*-* OOP row is still a table hit");

    assert_eq!(params.koop, 0.0);
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_returns_none_when_invalid() {
    let molecule = four_atom_oop_molecule(Element::C, Element::C, Element::C, Element::C);
    let mut props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 2, 1, 2]);
    props.valid = false;

    let params = props
        .get_mmff_oop_bend_params(0, 1, 2, 3)
        .expect("invalid properties return RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_returns_none_without_first_bond() {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    for (begin, end) in [(a1, a2), (a1, a3)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 2, 1, 2]);

    let params = props
        .get_mmff_oop_bend_params(a0.index(), a1.index(), a2.index(), a3.index())
        .expect("missing first OOP bond returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_returns_none_without_second_bond() {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    for (begin, end) in [(a0, a1), (a1, a3)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 2, 1, 2]);

    let params = props
        .get_mmff_oop_bend_params(a0.index(), a1.index(), a2.index(), a3.index())
        .expect("missing second OOP bond returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_returns_none_without_third_bond() {
    let mut builder = TestFixtureBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    for (begin, end) in [(a0, a1), (a1, a2)] {
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
    }
    let molecule = builder.build().expect("test molecule is valid");
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 2, 1, 2]);

    let params = props
        .get_mmff_oop_bend_params(a0.index(), a1.index(), a2.index(), a3.index())
        .expect("missing third OOP bond returns RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_reports_atom_index_out_of_range() {
    let molecule = four_atom_oop_molecule(Element::C, Element::C, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 2, 1]);

    let err = props.get_mmff_oop_bend_params(0, 1, 2, 3).unwrap_err();

    match err {
        MmffMolPropertiesError::AtomIndexOutOfRange { atom_index, atoms } => {
            assert_eq!(atom_index, 3);
            assert_eq!(atoms, 3);
        }
        other => panic!("expected atom-index out-of-range error, got {other:?}"),
    }
}

#[test]
fn mmff_mol_properties_get_oop_bend_params_returns_none_without_table_row() {
    let molecule = four_atom_oop_molecule(Element::C, Element::C, Element::C, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 1, 1, 1]);

    let params = props
        .get_mmff_oop_bend_params(0, 1, 2, 3)
        .expect("default MMFF OOP table parses");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_vdw_params_returns_unscaled_non_da_pair() {
    let props = mmff_props_for_atom_types(&[1, 1]);

    let params = props
        .get_mmff_vdw_params(0, 1)
        .expect("default MMFF VdW table parses")
        .expect("atom type 1 has tabulated VdW params");

    assert_eq!((params.r_ij_star_unscaled * 1000.0).round() as i32, 3938);
    assert_eq!((params.epsilon_unscaled * 1000.0).round() as i32, 68);
    assert_eq!(params.r_ij_star, params.r_ij_star_unscaled);
    assert_eq!(params.epsilon, params.epsilon_unscaled);
}

#[test]
fn mmff_mol_properties_get_vdw_params_scales_da_pair_like_rdkit() {
    let props = mmff_props_for_atom_types(&[8, 23]);

    let params = props
        .get_mmff_vdw_params(0, 1)
        .expect("default MMFF VdW table parses")
        .expect("atom types 8 and 23 have tabulated VdW params");

    assert_eq!((params.r_ij_star_unscaled * 1000.0).round() as i32, 3321);
    assert_eq!((params.epsilon_unscaled * 1000.0).round() as i32, 34);
    assert_eq!((params.r_ij_star * 1000.0).round() as i32, 2657);
    assert_eq!((params.epsilon * 1000.0).round() as i32, 17);
}

#[test]
fn mmff_mol_properties_get_vdw_params_scales_reversed_da_pair() {
    let props = mmff_props_for_atom_types(&[23, 8]);

    let params = props
        .get_mmff_vdw_params(0, 1)
        .expect("default MMFF VdW table parses")
        .expect("atom types 23 and 8 have tabulated VdW params");

    assert_eq!((params.r_ij_star_unscaled * 1000.0).round() as i32, 3321);
    assert_eq!((params.epsilon_unscaled * 1000.0).round() as i32, 34);
    assert_eq!((params.r_ij_star * 1000.0).round() as i32, 2657);
    assert_eq!((params.epsilon * 1000.0).round() as i32, 17);
}

#[test]
fn mmff_mol_properties_get_vdw_params_returns_none_when_invalid() {
    let mut props = mmff_props_for_atom_types(&[8, 23]);
    props.valid = false;

    let params = props
        .get_mmff_vdw_params(0, 1)
        .expect("invalid properties return RDKit false equivalent");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_vdw_params_returns_none_without_first_table_row() {
    let props = mmff_props_for_atom_types(&[83, 1]);

    let params = props
        .get_mmff_vdw_params(0, 1)
        .expect("default MMFF VdW table parses");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_vdw_params_returns_none_without_second_table_row() {
    let props = mmff_props_for_atom_types(&[1, 83]);

    let params = props
        .get_mmff_vdw_params(0, 1)
        .expect("default MMFF VdW table parses");

    assert_eq!(params, None);
}

#[test]
fn mmff_mol_properties_get_vdw_params_reports_first_atom_index_out_of_range() {
    let props = mmff_props_for_atom_types(&[1]);

    let err = props.get_mmff_vdw_params(1, 0).unwrap_err();

    match err {
        MmffMolPropertiesError::AtomIndexOutOfRange { atom_index, atoms } => {
            assert_eq!(atom_index, 1);
            assert_eq!(atoms, 1);
        }
        other => panic!("expected atom-index out-of-range error, got {other:?}"),
    }
}

#[test]
fn mmff_mol_properties_get_vdw_params_reports_second_atom_index_out_of_range() {
    let props = mmff_props_for_atom_types(&[1]);

    let err = props.get_mmff_vdw_params(0, 1).unwrap_err();

    match err {
        MmffMolPropertiesError::AtomIndexOutOfRange { atom_index, atoms } => {
            assert_eq!(atom_index, 1);
            assert_eq!(atoms, 1);
        }
        other => panic!("expected atom-index out-of-range error, got {other:?}"),
    }
}
#[test]
fn mmff_properties_defined_si_si_missing_table_uses_empirical_fallback() {
    let molecule = four_atom_torsion_molecule(Element::C, Element::SI, Element::SI, Element::C);
    let props = mmff_props_for_molecule_and_atom_types(molecule, &[1, 19, 19, 1]);
    let defs = default_mmff_def().unwrap();
    let table = default_mmff_tor(false).unwrap();
    assert_eq!(table.get(defs, (0, 0), 1, 19, 19, 1).unwrap().1, None);
    let actual = props.get_mmff_torsion_params(0, 1, 2, 3).unwrap().unwrap();
    let empirical = props
        .get_mmff_torsion_empirical_rule_params(1, 2, default_mmff_prop().unwrap())
        .unwrap();
    assert_eq!(actual.0, 0);
    assert_eq!(actual.1, empirical);
}

#[test]
fn mmff_mol_properties_constructor_initializes_empty_molecule_defaults() {
    let molecule = TestInput::new();
    let props = MmffMolProperties::new(&molecule, false, "MMFF94", MMFF_VERBOSITY_HIGH)
        .expect("empty molecule has no atom typing work");
    assert!(props.is_valid());
    assert_eq!(props.mmff_variant(), MmffVariant::Mmff94);
    assert!(props.bond_term);
    assert!(props.angle_term);
    assert!(props.stretch_bend_term);
    assert!(props.oop_term);
    assert!(props.torsion_term);
    assert!(props.vdw_term);
    assert!(props.ele_term);
    assert_eq!(props.dielectric_constant, 1.0);
    assert_eq!(props.dielectric_model, MMFF_DIELECTRIC_CONSTANT);
    assert_eq!(props.verbosity, MMFF_VERBOSITY_HIGH);
    assert!(props.atom_properties.is_empty());
    // _MMFFSanitized lives in the caller's property block. Its guarded
    // computed installation and existing-value preservation are exercised
    // by mmff_live_source_property_guard_ring_rows_and_unchanged_block_sharing.
}
