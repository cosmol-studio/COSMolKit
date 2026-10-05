//! Fixed public metadata vectors from RDKit 2026.03.1, BSD-3-Clause.
//! Source revision: 351f8f378f8ad6bbd517980c38896e66bf907af8.
//! Source: Code/GraphMol/atomic_data.cpp, periodicTableAtomData; select the
//! first numeric row for each atomic number 0..=118, retaining valence order.
//! Source SHA-256: 7f9cee6e430b60d303a0a7fa9e33c45e5c529ee204f86afeab6a20f68b6b0631.
//! Literal expectations are fixed at development time. Tests never invoke
//! an oracle, parse upstream source or derive expected values from core.

use cosmolkit::{
    BINDING_CONTRACT, BindingItem, BindingOwner, BindingTypeRole, Element, ElementInfo,
    FunctionStatus,
};
#[cfg(feature = "cap-valence")]
use cosmolkit::{BindingDefault, BindingKind, StateModel, element_info};

fn entry(id: &str) -> &'static cosmolkit::BindingContractEntry {
    let mut matches = BINDING_CONTRACT.iter().filter(|row| row.semantic_id == id);
    let row = matches.next().unwrap_or_else(|| panic!("missing {id}"));
    assert!(matches.next().is_none(), "duplicate registration for {id}");
    row
}

#[test]
fn canonical_metadata_types_remain_available_without_capabilities() {
    for (id, name, role) in [
        ("types.Element", "Element", BindingTypeRole::Value),
        ("types.ElementInfo", "ElementInfo", BindingTypeRole::Result),
    ] {
        let row = entry(id);
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.rust_path.replace(' ', ""), format!("crate::{name}"));
        assert_eq!(row.python_name, name);
        assert_eq!(row.javascript_name, name);
        assert_eq!(row.feature, "metadata");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.type_role, Some(role));
        assert_eq!(row.callable, None);
    }
    let _: Option<ElementInfo> = None;
    for number in 0..=118 {
        assert_eq!(
            Element::from_atomic_number(number).unwrap().atomic_number(),
            number
        );
    }
    for number in 119..=u8::MAX {
        assert!(Element::from_atomic_number(number).is_none());
    }
}

#[cfg(not(feature = "cap-valence"))]
#[test]
fn element_info_callable_is_absent_without_its_capability() {
    assert!(
        BINDING_CONTRACT
            .iter()
            .all(|row| row.semantic_id != "module.element_info")
    );
}

#[cfg(feature = "cap-valence")]
#[test]
fn element_info_registration_has_the_exact_public_contract() {
    let _: fn(Element) -> ElementInfo = element_info;
    let row = entry("module.element_info");
    assert_eq!(row.item, BindingItem::Callable);
    assert_eq!(row.owner, BindingOwner::Module);
    assert_eq!(row.rust_path.replace(' ', ""), "crate::element_info");
    assert_eq!(row.python_name, "element_info");
    assert_eq!(row.javascript_name, "elementInfo");
    assert_eq!(row.feature, "cap-valence");
    assert_eq!(row.status, FunctionStatus::Experimental);
    assert_eq!(row.type_role, None);
    let callable = row.callable.unwrap();
    assert_eq!(callable.kind, BindingKind::Module);
    assert_eq!(callable.receiver, None);
    assert_eq!(callable.state_model, StateModel::ReadOnly);
    assert_eq!(callable.operation_semantic_id, None);
    assert_eq!(callable.error_type, None);
    assert_eq!(callable.output_type.replace(' ', ""), "crate::ElementInfo");
    assert_eq!(callable.parameters.len(), 1);
    let parameter = callable.parameters[0];
    assert_eq!(parameter.name, "element");
    assert_eq!(parameter.type_name.replace(' ', ""), "crate::Element");
    assert_eq!(parameter.default, BindingDefault::Required);
}

#[cfg(feature = "cap-valence")]
#[test]
fn all_119_elements_return_exact_source_metadata_through_the_facade() {
    assert_eq!(EXPECTED.len(), 119);
    for (index, &(number, symbol, period, outer_electrons, valences, rb0, weight)) in
        EXPECTED.iter().enumerate()
    {
        assert_eq!(
            usize::from(number),
            index,
            "fixed vectors must be contiguous"
        );
        let element = Element::from_atomic_number(number).unwrap();
        let actual = element_info(element);
        assert_eq!(actual.element, element, "element identity at {number}");
        assert_eq!(actual.atomic_number, number, "atomic number at {number}");
        assert_eq!(actual.symbol, symbol, "symbol at {number}");
        assert_eq!(actual.period, period, "period at {number}");
        assert_eq!(
            actual.outer_electrons, outer_electrons,
            "outer electrons at {number}"
        );
        assert_eq!(
            actual.valences, valences,
            "complete ordered valences at {number}"
        );
        assert_eq!(actual.rb0.to_bits(), rb0.to_bits(), "Rb0 at {number}");
        assert_eq!(
            actual.atomic_weight.to_bits(),
            weight.to_bits(),
            "weight at {number}"
        );
    }
}

#[cfg(feature = "cap-valence")]
#[test]
fn returned_metadata_borrows_stable_immutable_data() {
    for number in 0..=118 {
        let element = Element::from_atomic_number(number).unwrap();
        let first = element_info(element);
        let second = element_info(element);
        assert_eq!(first, second);
        assert!(std::ptr::eq(first.symbol, second.symbol));
        assert!(std::ptr::eq(first.valences, second.valences));
        let _: &'static str = first.symbol;
        let _: &'static [i32] = first.valences;
    }
    // Vocabulary aliases use the canonical numeric row, not a second table.
    assert_eq!(
        element_info(Element::from_symbol("Uut").unwrap()),
        element_info(Element::NH)
    );
    assert_eq!(
        element_info(Element::from_symbol("Uup").unwrap()),
        element_info(Element::MC)
    );
}

#[cfg(feature = "cap-valence")]
type ExpectedRow = (u8, &'static str, u8, i32, &'static [i32], f64, f64);

#[cfg(feature = "cap-valence")]
const EXPECTED: [ExpectedRow; 119] = [
    (0, "*", 0, 0, &[-1], 0.0, 0.0),
    (1, "H", 1, 1, &[1], 0.33, 1.008),
    (2, "He", 1, 2, &[0], 0.7, 4.003),
    (3, "Li", 2, 1, &[1, -1], 1.23, 6.941),
    (4, "Be", 2, 2, &[2], 0.9, 9.012),
    (5, "B", 2, 3, &[3], 0.82, 10.812),
    (6, "C", 2, 4, &[4], 0.77, 12.011),
    (7, "N", 2, 5, &[3], 0.7, 14.007),
    (8, "O", 2, 6, &[2], 0.66, 15.999),
    (9, "F", 2, 7, &[1], 0.611, 18.998),
    (10, "Ne", 2, 8, &[0], 0.7, 20.18),
    (11, "Na", 3, 1, &[1, -1], 1.54, 22.99),
    (12, "Mg", 3, 2, &[2, -1], 1.36, 24.305),
    (13, "Al", 3, 3, &[3], 1.18, 26.982),
    (14, "Si", 3, 4, &[4], 0.937, 28.086),
    (15, "P", 3, 5, &[3, 5], 0.89, 30.974),
    (16, "S", 3, 6, &[2, 4, 6], 1.04, 32.067),
    (17, "Cl", 3, 7, &[1], 0.997, 35.453),
    (18, "Ar", 3, 8, &[0], 1.74, 39.948),
    (19, "K", 4, 1, &[1, -1], 2.03, 39.098),
    (20, "Ca", 4, 2, &[2, -1], 1.74, 40.078),
    (21, "Sc", 4, 3, &[-1], 1.44, 44.956),
    (22, "Ti", 4, 4, &[-1], 1.32, 47.867),
    (23, "V", 4, 5, &[-1], 1.22, 50.944),
    (24, "Cr", 4, 6, &[-1], 1.18, 51.996),
    (25, "Mn", 4, 7, &[-1], 1.17, 54.938),
    (26, "Fe", 4, 8, &[-1], 1.17, 55.845),
    (27, "Co", 4, 9, &[-1], 1.16, 58.933),
    (28, "Ni", 4, 10, &[-1], 1.15, 58.693),
    (29, "Cu", 4, 11, &[-1], 1.17, 63.546),
    (30, "Zn", 4, 2, &[-1], 1.25, 65.39),
    (31, "Ga", 4, 3, &[3], 1.26, 69.723),
    (32, "Ge", 4, 4, &[4], 1.188, 72.61),
    (33, "As", 4, 5, &[3, 5], 1.2, 74.922),
    (34, "Se", 4, 6, &[2, 4, 6], 1.17, 78.96),
    (35, "Br", 4, 7, &[1], 1.167, 79.904),
    (36, "Kr", 4, 8, &[0], 1.91, 83.8),
    (37, "Rb", 5, 1, &[1, -1], 2.16, 85.468),
    (38, "Sr", 5, 2, &[2, -1], 1.91, 87.62),
    (39, "Y", 5, 3, &[-1], 1.62, 88.906),
    (40, "Zr", 5, 4, &[-1], 1.45, 91.224),
    (41, "Nb", 5, 5, &[-1], 1.34, 92.906),
    (42, "Mo", 5, 6, &[-1], 1.3, 95.94),
    (43, "Tc", 5, 7, &[-1], 1.27, 98.0),
    (44, "Ru", 5, 8, &[-1], 1.25, 101.07),
    (45, "Rh", 5, 9, &[-1], 1.25, 102.906),
    (46, "Pd", 5, 10, &[-1], 1.28, 106.42),
    (47, "Ag", 5, 11, &[-1], 1.34, 107.868),
    (48, "Cd", 5, 2, &[-1], 1.48, 112.412),
    (49, "In", 5, 3, &[3], 1.44, 114.818),
    (50, "Sn", 5, 4, &[2, 4], 1.385, 118.711),
    (51, "Sb", 5, 5, &[3, 5], 1.4, 121.76),
    (52, "Te", 5, 6, &[2, 4, 6], 1.378, 127.6),
    (53, "I", 5, 7, &[1, 3, 5], 1.387, 126.904),
    (54, "Xe", 5, 8, &[0, 2, 4, 6], 1.98, 131.29),
    (55, "Cs", 6, 1, &[1], 2.35, 132.905),
    (56, "Ba", 6, 2, &[2, -1], 1.98, 137.328),
    (57, "La", 6, 3, &[-1], 1.69, 138.906),
    (58, "Ce", 6, 4, &[-1], 1.83, 140.116),
    (59, "Pr", 6, 3, &[-1], 1.82, 140.908),
    (60, "Nd", 6, 4, &[-1], 1.81, 144.24),
    (61, "Pm", 6, 5, &[-1], 1.8, 145.0),
    (62, "Sm", 6, 6, &[-1], 1.8, 150.36),
    (63, "Eu", 6, 7, &[-1], 1.99, 151.964),
    (64, "Gd", 6, 8, &[-1], 1.79, 157.25),
    (65, "Tb", 6, 9, &[-1], 1.76, 158.925),
    (66, "Dy", 6, 10, &[-1], 1.75, 162.5),
    (67, "Ho", 6, 11, &[-1], 1.74, 164.93),
    (68, "Er", 6, 12, &[-1], 1.73, 167.26),
    (69, "Tm", 6, 13, &[-1], 1.72, 168.934),
    (70, "Yb", 6, 14, &[-1], 1.94, 173.04),
    (71, "Lu", 6, 15, &[-1], 1.72, 174.967),
    (72, "Hf", 6, 4, &[-1], 1.44, 178.49),
    (73, "Ta", 6, 5, &[-1], 1.34, 180.948),
    (74, "W", 6, 6, &[-1], 1.3, 183.84),
    (75, "Re", 6, 7, &[-1], 1.28, 186.207),
    (76, "Os", 6, 8, &[-1], 1.26, 190.23),
    (77, "Ir", 6, 9, &[-1], 1.27, 192.217),
    (78, "Pt", 6, 10, &[-1], 1.3, 195.078),
    (79, "Au", 6, 11, &[-1], 1.34, 196.967),
    (80, "Hg", 6, 2, &[-1], 1.49, 200.59),
    (81, "Tl", 6, 3, &[-1], 1.48, 204.383),
    (82, "Pb", 6, 4, &[2, 4], 1.48, 207.2),
    (83, "Bi", 6, 5, &[3, 5], 1.45, 208.98),
    (84, "Po", 6, 6, &[2, 4, 6], 1.46, 209.0),
    (85, "At", 6, 7, &[1, 3, 5], 1.45, 210.0),
    (86, "Rn", 6, 8, &[0], 2.4, 222.0),
    (87, "Fr", 7, 1, &[1], 2.0, 223.0),
    (88, "Ra", 7, 2, &[2, -1], 1.9, 226.0),
    (89, "Ac", 7, 3, &[-1], 1.88, 227.0),
    (90, "Th", 7, 4, &[-1], 1.79, 232.038),
    (91, "Pa", 7, 3, &[-1], 1.61, 231.036),
    (92, "U", 7, 4, &[-1], 1.58, 238.029),
    (93, "Np", 7, 5, &[-1], 1.55, 237.0),
    (94, "Pu", 7, 6, &[-1], 1.53, 244.0),
    (95, "Am", 7, 7, &[-1], 1.07, 243.0),
    (96, "Cm", 7, 8, &[-1], 0.0, 247.0),
    (97, "Bk", 7, 9, &[-1], 0.0, 247.0),
    (98, "Cf", 7, 10, &[-1], 0.0, 251.0),
    (99, "Es", 7, 11, &[-1], 0.0, 252.0),
    (100, "Fm", 7, 12, &[-1], 0.0, 257.0),
    (101, "Md", 7, 13, &[-1], 0.0, 258.0),
    (102, "No", 7, 14, &[-1], 0.0, 259.0),
    (103, "Lr", 7, 15, &[-1], 0.0, 262.0),
    (104, "Rf", 7, 2, &[-1], 0.0, 267.0),
    (105, "Db", 7, 2, &[-1], 0.0, 268.0),
    (106, "Sg", 7, 2, &[-1], 0.0, 269.0),
    (107, "Bh", 7, 2, &[-1], 0.0, 270.0),
    (108, "Hs", 7, 2, &[-1], 0.0, 269.0),
    (109, "Mt", 7, 2, &[-1], 0.0, 278.0),
    (110, "Ds", 7, 2, &[-1], 0.0, 281.0),
    (111, "Rg", 7, 2, &[-1], 0.0, 281.0),
    (112, "Cn", 7, 2, &[-1], 0.0, 285.0),
    (113, "Nh", 7, 2, &[-1], 0.0, 284.0),
    (114, "Fl", 7, 2, &[-1], 0.0, 289.0),
    (115, "Mc", 7, 2, &[-1], 0.0, 288.0),
    (116, "Lv", 7, 2, &[-1], 0.0, 293.0),
    (117, "Ts", 7, 2, &[-1], 0.0, 292.0),
    (118, "Og", 7, 2, &[-1], 0.0, 294.0),
];
