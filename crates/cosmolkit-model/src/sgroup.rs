use std::collections::BTreeMap;

use crate::{AtomId, BondId};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum SubstanceGroupValidationError {
    IdMismatch {
        position: usize,
        id: SubstanceGroupId,
    },
    AtomOutOfRange {
        sgroup: SubstanceGroupId,
        atom: AtomId,
        atom_count: usize,
    },
    BondOutOfRange {
        sgroup: SubstanceGroupId,
        bond: BondId,
        bond_count: usize,
    },
    ParentOutOfRange {
        sgroup: SubstanceGroupId,
        parent: SubstanceGroupId,
    },
}

pub(crate) fn validate_substance_groups(
    groups: &[SubstanceGroup],
    atom_count: usize,
    bond_count: usize,
) -> Result<(), SubstanceGroupValidationError> {
    let group_count = groups.len();
    for (position, group) in groups.iter().enumerate() {
        if group.id() != SubstanceGroupId::new(position) {
            return Err(SubstanceGroupValidationError::IdMismatch {
                position,
                id: group.id(),
            });
        }
        for atom in group.atoms().iter().chain(group.parent_atoms()) {
            if atom.index() >= atom_count {
                return Err(SubstanceGroupValidationError::AtomOutOfRange {
                    sgroup: group.id(),
                    atom: *atom,
                    atom_count,
                });
            }
        }
        for point in group.attach_points() {
            for atom in std::iter::once(point.atom).chain(point.leaving_atom) {
                if atom.index() >= atom_count {
                    return Err(SubstanceGroupValidationError::AtomOutOfRange {
                        sgroup: group.id(),
                        atom,
                        atom_count,
                    });
                }
            }
        }
        for bond in group
            .bonds()
            .iter()
            .chain(group.cstates().iter().map(|state| &state.bond))
            .chain(group.head_crossing_bonds())
            .chain(group.crossing_bond_correspondence())
        {
            if bond.index() >= bond_count {
                return Err(SubstanceGroupValidationError::BondOutOfRange {
                    sgroup: group.id(),
                    bond: *bond,
                    bond_count,
                });
            }
        }
        if let Some(parent) = group.parent()
            && parent.index() >= group_count
        {
            return Err(SubstanceGroupValidationError::ParentOutOfRange {
                sgroup: group.id(),
                parent,
            });
        }
    }
    Ok(())
}

fn remove_first<T: PartialEq>(container: &mut Vec<T>, element: &T) -> bool {
    // RDKit✔️✔️: auto pos = std::find(container.begin(), container.end(), element);
    // RDKit✔️✔️: if (pos != container.end()) {
    // RDKit✔️✔️:   container.erase(pos);
    // RDKit✔️✔️: }
    let Some(position) = container.iter().position(|candidate| candidate == element) else {
        return false;
    };
    container.remove(position);
    true
}

/// Substance-group identity inside a molecule.
///
/// RDKit SDF/MolBlock parsing can produce SGroups before they are interpreted
/// by higher-level chemistry operations. Keep this as explicit molecule state;
/// do not flatten SGroup data into string properties without human-author
/// approval, because that loses atom/bond membership and processing semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SubstanceGroupId(usize);

impl SubstanceGroupId {
    #[must_use]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    #[must_use]
    pub const fn index(&self) -> usize {
        self.0
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SubstanceGroupKind {
    Data,
    Superatom,
    MultipleGroup,
    StructuralRepeatUnit,
    Monomer,
    Copolymer,
    Crosslink,
    Graft,
    Modification,
    Mer,
    AnyPolymer,
    MixtureComponent,
    Mixture,
    Formulation,
    Generic(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SGroupBondRole {
    Crossing,
    Contained,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SGroupBracket {
    pub points: [[f64; 3]; 3],
}

impl SGroupBracket {
    #[must_use]
    pub const fn new(points: [[f64; 3]; 3]) -> Self {
        Self { points }
    }

    #[must_use]
    pub const fn points(&self) -> &[[f64; 3]; 3] {
        &self.points
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SGroupCState {
    pub bond: BondId,
    pub vector: [f64; 3],
}

impl SGroupCState {
    #[must_use]
    pub const fn new(bond: BondId, vector: [f64; 3]) -> Self {
        Self { bond, vector }
    }

    #[must_use]
    pub const fn bond(&self) -> BondId {
        self.bond
    }

    #[must_use]
    pub const fn vector(&self) -> &[f64; 3] {
        &self.vector
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct SGroupDisplay {
    pub brackets: Vec<SGroupBracket>,
    pub field_position: Option<[f64; 2]>,
    pub display_tag: Option<String>,
}

impl SGroupDisplay {
    #[must_use]
    pub fn brackets(&self) -> &[SGroupBracket] {
        &self.brackets
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SGroupData {
    pub field_name: Option<String>,
    pub field_type: Option<String>,
    pub field_info: Option<String>,
    pub field_display: Option<String>,
    pub units: Option<String>,
    pub query_type: Option<String>,
    pub query_op: Option<String>,
    pub values: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SGroupConnection {
    HeadToHead,
    HeadToTail,
    Either,
    Unknown(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SGroupBracketStyle {
    Bracket,
    Parenthesis,
    None,
    Unknown(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SGroupAttachPoint {
    pub atom: AtomId,
    pub leaving_atom: Option<AtomId>,
    pub label: Option<String>,
    pub order: Option<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SubstanceGroup {
    id: SubstanceGroupId,
    rdkit_sequence_id: Option<u32>,
    external_id: Option<u32>,
    kind: SubstanceGroupKind,
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
    bond_roles: BTreeMap<BondId, SGroupBondRole>,
    head_crossing_bonds: Vec<BondId>,
    crossing_bond_correspondence: Vec<BondId>,
    parent_atoms: Vec<AtomId>,
    parent: Option<SubstanceGroupId>,
    label: Option<String>,
    connection: Option<SGroupConnection>,
    subtype: Option<String>,
    bracket_style: Option<SGroupBracketStyle>,
    expansion_state: Option<String>,
    class: Option<String>,
    component_number: Option<u32>,
    display: Option<SGroupDisplay>,
    data: Option<SGroupData>,
    attach_points: Vec<SGroupAttachPoint>,
    cstates: Vec<SGroupCState>,
    props: BTreeMap<String, String>,
    data_fields: Vec<String>,
}

impl SubstanceGroup {
    #[must_use]
    pub fn new(id: SubstanceGroupId, kind: SubstanceGroupKind) -> Self {
        Self {
            id,
            rdkit_sequence_id: None,
            external_id: None,
            kind,
            atoms: Vec::new(),
            bonds: Vec::new(),
            bond_roles: BTreeMap::new(),
            head_crossing_bonds: Vec::new(),
            crossing_bond_correspondence: Vec::new(),
            parent_atoms: Vec::new(),
            parent: None,
            label: None,
            connection: None,
            subtype: None,
            bracket_style: None,
            expansion_state: None,
            class: None,
            component_number: None,
            display: None,
            data: None,
            attach_points: Vec::new(),
            cstates: Vec::new(),
            props: BTreeMap::new(),
            data_fields: Vec::new(),
        }
    }

    #[must_use]
    pub const fn id(&self) -> SubstanceGroupId {
        self.id
    }

    #[must_use]
    pub const fn external_id(&self) -> Option<u32> {
        self.external_id
    }

    #[must_use]
    pub const fn rdkit_sequence_id(&self) -> Option<u32> {
        self.rdkit_sequence_id
    }

    #[must_use]
    pub const fn kind(&self) -> &SubstanceGroupKind {
        &self.kind
    }

    #[must_use]
    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }

    #[must_use]
    pub fn bonds(&self) -> &[BondId] {
        &self.bonds
    }

    /// Ordered crossing-bond references encoded by V3000 `XBHEAD` and CX
    /// polymer head-crossing state.
    #[must_use]
    pub fn head_crossing_bonds(&self) -> &[BondId] {
        &self.head_crossing_bonds
    }

    /// Ordered crossing-bond correspondence references encoded by V3000
    /// `XBCORR` and CX polymer tail-crossing state.
    #[must_use]
    pub fn crossing_bond_correspondence(&self) -> &[BondId] {
        &self.crossing_bond_correspondence
    }

    #[must_use]
    pub fn bond_role(&self, bond: BondId) -> SGroupBondRole {
        self.bond_roles
            .get(&bond)
            .copied()
            .unwrap_or(SGroupBondRole::Crossing)
    }

    #[must_use]
    pub fn parent_atoms(&self) -> &[AtomId] {
        &self.parent_atoms
    }

    #[must_use]
    pub const fn parent(&self) -> Option<SubstanceGroupId> {
        self.parent
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn data_fields(&self) -> &[String] {
        &self.data_fields
    }

    #[must_use]
    pub fn label(&self) -> Option<&str> {
        self.label.as_deref()
    }

    #[must_use]
    pub const fn connection(&self) -> Option<&SGroupConnection> {
        self.connection.as_ref()
    }

    #[must_use]
    pub fn subtype(&self) -> Option<&str> {
        self.subtype.as_deref()
    }

    #[must_use]
    pub const fn bracket_style(&self) -> Option<&SGroupBracketStyle> {
        self.bracket_style.as_ref()
    }

    #[must_use]
    pub const fn display(&self) -> Option<&SGroupDisplay> {
        self.display.as_ref()
    }

    #[must_use]
    pub fn expansion_state(&self) -> Option<&str> {
        self.expansion_state.as_deref()
    }

    #[must_use]
    pub fn class(&self) -> Option<&str> {
        self.class.as_deref()
    }

    #[must_use]
    pub const fn component_number(&self) -> Option<u32> {
        self.component_number
    }

    #[must_use]
    pub const fn data(&self) -> Option<&SGroupData> {
        self.data.as_ref()
    }

    #[must_use]
    pub fn attach_points(&self) -> &[SGroupAttachPoint] {
        &self.attach_points
    }

    #[must_use]
    pub fn cstates(&self) -> &[SGroupCState] {
        &self.cstates
    }

    #[must_use]
    pub const fn with_rdkit_sequence_id(mut self, rdkit_sequence_id: u32) -> Self {
        self.rdkit_sequence_id = Some(rdkit_sequence_id);
        self
    }

    #[must_use]
    pub const fn with_external_id(mut self, external_id: u32) -> Self {
        self.external_id = Some(external_id);
        self
    }

    #[must_use]
    pub const fn with_parent(mut self, parent: SubstanceGroupId) -> Self {
        self.parent = Some(parent);
        self
    }

    #[must_use]
    pub fn with_atoms(mut self, atoms: Vec<AtomId>) -> Self {
        self.atoms = atoms;
        self
    }

    #[must_use]
    pub fn with_bonds(mut self, bonds: Vec<BondId>) -> Self {
        self.bond_roles
            .retain(|bond, _| bonds.iter().any(|candidate| candidate == bond));
        self.bonds = bonds;
        self
    }

    #[must_use]
    pub fn with_bond_role(mut self, bond: BondId, role: SGroupBondRole) -> Self {
        if self.bonds.iter().any(|candidate| *candidate == bond) {
            self.bond_roles.insert(bond, role);
        }
        self
    }

    #[must_use]
    pub fn with_head_crossing_bonds(mut self, bonds: Vec<BondId>) -> Self {
        self.head_crossing_bonds = bonds;
        self
    }

    #[must_use]
    pub fn with_crossing_bond_correspondence(mut self, bonds: Vec<BondId>) -> Self {
        self.crossing_bond_correspondence = bonds;
        self
    }

    #[must_use]
    pub fn with_parent_atoms(mut self, parent_atoms: Vec<AtomId>) -> Self {
        self.parent_atoms = parent_atoms;
        self
    }

    #[must_use]
    pub fn with_label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    #[must_use]
    pub fn with_connection(mut self, connection: SGroupConnection) -> Self {
        self.connection = Some(connection);
        self
    }

    #[must_use]
    pub fn with_subtype(mut self, subtype: impl Into<String>) -> Self {
        self.subtype = Some(subtype.into());
        self
    }

    #[must_use]
    pub fn with_bracket_style(mut self, bracket_style: SGroupBracketStyle) -> Self {
        self.bracket_style = Some(bracket_style);
        self
    }

    #[must_use]
    pub fn with_display(mut self, display: SGroupDisplay) -> Self {
        self.display = Some(display);
        self
    }

    #[must_use]
    pub fn with_expansion_state(mut self, expansion_state: impl Into<String>) -> Self {
        self.expansion_state = Some(expansion_state.into());
        self
    }

    #[must_use]
    pub fn with_class(mut self, class: impl Into<String>) -> Self {
        self.class = Some(class.into());
        self
    }

    #[must_use]
    pub const fn with_component_number(mut self, component_number: u32) -> Self {
        self.component_number = Some(component_number);
        self
    }

    #[must_use]
    pub fn with_data(mut self, data: SGroupData) -> Self {
        self.data = Some(data);
        self
    }

    #[must_use]
    pub fn with_attach_points(mut self, attach_points: Vec<SGroupAttachPoint>) -> Self {
        self.attach_points = attach_points;
        self
    }

    #[must_use]
    pub fn with_cstates(mut self, cstates: Vec<SGroupCState>) -> Self {
        self.cstates = cstates;
        self
    }

    #[must_use]
    pub fn with_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub fn with_data_field(mut self, value: impl Into<String>) -> Self {
        self.data_fields.push(value.into());
        self
    }

    pub fn push_data_field(&mut self, value: impl Into<String>) {
        self.data_fields.push(value.into());
    }

    pub fn push_atom(&mut self, atom: AtomId) {
        self.atoms.push(atom);
    }

    pub fn push_bond(&mut self, bond: BondId) {
        self.push_bond_with_role(bond, SGroupBondRole::Crossing);
    }

    pub fn push_bond_with_role(&mut self, bond: BondId, role: SGroupBondRole) {
        self.bonds.push(bond);
        if role != SGroupBondRole::Crossing {
            self.bond_roles.insert(bond, role);
        }
    }

    pub fn push_head_crossing_bond(&mut self, bond: BondId) {
        self.head_crossing_bonds.push(bond);
    }

    pub fn push_crossing_bond_correspondence(&mut self, bond: BondId) {
        self.crossing_bond_correspondence.push(bond);
    }

    pub fn remove_atom(&mut self, atom: AtomId) {
        remove_first(&mut self.atoms, &atom);
    }

    pub fn remove_bond(&mut self, bond: BondId) {
        if remove_first(&mut self.bonds, &bond)
            && !self.bonds.iter().any(|candidate| *candidate == bond)
        {
            self.bond_roles.remove(&bond);
        }
    }

    pub fn remove_parent_atom(&mut self, atom: AtomId) {
        remove_first(&mut self.parent_atoms, &atom);
    }

    pub fn clear_attach_point_leaving_atom(&mut self, atom: AtomId) {
        for attach_point in &mut self.attach_points {
            if attach_point.leaving_atom == Some(atom) {
                attach_point.leaving_atom = None;
            }
        }
    }

    #[allow(dead_code)]
    pub fn push_parent_atom(&mut self, atom: AtomId) {
        self.parent_atoms.push(atom);
    }

    pub fn set_external_id(&mut self, external_id: u32) {
        self.external_id = Some(external_id);
    }

    #[allow(dead_code)]
    pub fn set_parent(&mut self, parent: SubstanceGroupId) {
        self.parent = Some(parent);
    }

    #[allow(dead_code)]
    pub fn set_rdkit_sequence_id(&mut self, rdkit_sequence_id: u32) {
        self.rdkit_sequence_id = Some(rdkit_sequence_id);
    }

    pub fn set_id(&mut self, id: SubstanceGroupId) {
        self.id = id;
    }

    #[allow(dead_code)]
    pub fn set_label(&mut self, label: impl Into<String>) {
        self.label = Some(label.into());
    }

    #[allow(dead_code)]
    pub fn set_connection(&mut self, connection: SGroupConnection) {
        self.connection = Some(connection);
    }

    #[allow(dead_code)]
    pub fn set_subtype(&mut self, subtype: impl Into<String>) {
        self.subtype = Some(subtype.into());
    }

    #[allow(dead_code)]
    pub fn set_bracket_style(&mut self, bracket_style: SGroupBracketStyle) {
        self.bracket_style = Some(bracket_style);
    }

    #[allow(dead_code)]
    pub fn set_expansion_state(&mut self, expansion_state: impl Into<String>) {
        self.expansion_state = Some(expansion_state.into());
    }

    #[allow(dead_code)]
    pub fn set_class(&mut self, class: impl Into<String>) {
        self.class = Some(class.into());
    }

    #[allow(dead_code)]
    pub fn set_component_number(&mut self, component_number: u32) {
        self.component_number = Some(component_number);
    }

    #[allow(dead_code)]
    pub fn display_mut(&mut self) -> &mut SGroupDisplay {
        self.display.get_or_insert_with(SGroupDisplay::default)
    }

    #[allow(dead_code)]
    pub fn data_mut(&mut self) -> &mut SGroupData {
        self.data.get_or_insert_with(SGroupData::default)
    }

    #[allow(dead_code)]
    pub fn push_attach_point(&mut self, attach_point: SGroupAttachPoint) {
        self.attach_points.push(attach_point);
    }

    #[allow(dead_code)]
    pub fn push_cstate(&mut self, cstate: SGroupCState) {
        self.cstates.push(cstate);
    }

    pub fn set_prop(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.props.insert(key.into(), value.into());
    }

    pub fn clear_prop(&mut self, key: &str) {
        self.props.remove(key);
    }

    #[must_use]
    pub fn includes_atom(&self, atom: AtomId) -> bool {
        // RDKit✔️✔️: if (std::find(d_atoms.begin(), d_atoms.end(), atomIdx) != d_atoms.end()) {
        // RDKit✔️✔️:   return true;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: if (std::find(d_patoms.begin(), d_patoms.end(), atomIdx) != d_patoms.end()) {
        // RDKit✔️✔️:   return true;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: for (const auto &ap : d_saps) {
        // RDKit✔️✔️:   if (ap.aIdx == atomIdx || ap.lvIdx == rdcast<int>(atomIdx)) {
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: return false;
        self.atoms.contains(&atom)
            || self.parent_atoms.contains(&atom)
            || self.attach_points.iter().any(|attach_point| {
                attach_point.atom == atom || attach_point.leaving_atom == Some(atom)
            })
    }

    #[must_use]
    pub fn includes_bond(&self, bond: BondId) -> bool {
        // RDKit✔️✔️: if (std::find(d_bonds.begin(), d_bonds.end(), bondIdx) != d_bonds.end()) {
        // RDKit✔️✔️:   return true;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: for (const auto &cs : d_cstates) {
        // RDKit✔️✔️:   if (cs.bondIdx == bondIdx) {
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: return false;
        self.bonds.contains(&bond)
            || self.cstates.iter().any(|cstate| cstate.bond == bond)
            || self.head_crossing_bonds.contains(&bond)
            || self.crossing_bond_correspondence.contains(&bond)
    }

    pub fn can_remap_without_parent(
        &self,
        atom_map: &[Option<AtomId>],
        bond_map: &[Option<BondId>],
    ) -> bool {
        self.atoms
            .iter()
            .all(|atom| atom_map.get(atom.index()).is_some_and(Option::is_some))
            && self
                .bonds
                .iter()
                .all(|bond| bond_map.get(bond.index()).is_some_and(Option::is_some))
            && self
                .head_crossing_bonds
                .iter()
                .all(|bond| bond_map.get(bond.index()).is_some_and(Option::is_some))
            && self
                .crossing_bond_correspondence
                .iter()
                .all(|bond| bond_map.get(bond.index()).is_some_and(Option::is_some))
            && self
                .parent_atoms
                .iter()
                .all(|atom| atom_map.get(atom.index()).is_some_and(Option::is_some))
            && self.attach_points.iter().all(|attach_point| {
                atom_map
                    .get(attach_point.atom.index())
                    .is_some_and(Option::is_some)
                    && attach_point.leaving_atom.is_none_or(|leaving_atom| {
                        atom_map
                            .get(leaving_atom.index())
                            .is_some_and(Option::is_some)
                    })
            })
            && self.cstates.iter().all(|cstate| {
                bond_map
                    .get(cstate.bond.index())
                    .is_some_and(Option::is_some)
            })
    }

    pub fn remapped(
        &self,
        id: SubstanceGroupId,
        atom_map: &[Option<AtomId>],
        bond_map: &[Option<BondId>],
        sgroup_map: &[Option<SubstanceGroupId>],
    ) -> Option<Self> {
        let atoms: Option<Vec<_>> = self
            .atoms
            .iter()
            .map(|atom| atom_map.get(atom.index()).and_then(|x| *x))
            .collect();
        let bonds: Option<Vec<_>> = self
            .bonds
            .iter()
            .map(|bond| bond_map.get(bond.index()).and_then(|x| *x))
            .collect();
        let head_crossing_bonds: Option<Vec<_>> = self
            .head_crossing_bonds
            .iter()
            .map(|bond| bond_map.get(bond.index()).and_then(|x| *x))
            .collect();
        let crossing_bond_correspondence: Option<Vec<_>> = self
            .crossing_bond_correspondence
            .iter()
            .map(|bond| bond_map.get(bond.index()).and_then(|x| *x))
            .collect();
        let parent_atoms: Option<Vec<_>> = self
            .parent_atoms
            .iter()
            .map(|atom| atom_map.get(atom.index()).and_then(|x| *x))
            .collect();
        let parent = match self.parent {
            Some(parent) => sgroup_map.get(parent.index()).and_then(|x| *x),
            None => None,
        };
        if self.parent.is_some() && parent.is_none() {
            return None;
        }
        let attach_points: Option<Vec<_>> = self
            .attach_points
            .iter()
            .map(|attach_point| {
                let atom = atom_map.get(attach_point.atom.index()).and_then(|x| *x)?;
                let leaving_atom = match attach_point.leaving_atom {
                    Some(leaving_atom) => {
                        Some(atom_map.get(leaving_atom.index()).and_then(|x| *x)?)
                    }
                    None => None,
                };
                Some(SGroupAttachPoint {
                    atom,
                    leaving_atom,
                    label: attach_point.label.clone(),
                    order: attach_point.order,
                })
            })
            .collect();
        let cstates: Option<Vec<_>> = self
            .cstates
            .iter()
            .map(|cstate| {
                let bond = bond_map.get(cstate.bond.index()).and_then(|x| *x)?;
                Some(SGroupCState {
                    bond,
                    vector: cstate.vector,
                })
            })
            .collect();
        let mut bond_roles = BTreeMap::new();
        for (old_bond, role) in &self.bond_roles {
            let new_bond = bond_map.get(old_bond.index()).and_then(|x| *x)?;
            if *role != SGroupBondRole::Crossing {
                bond_roles.insert(new_bond, *role);
            }
        }
        Some(Self {
            id,
            rdkit_sequence_id: self.rdkit_sequence_id,
            external_id: self.external_id,
            kind: self.kind.clone(),
            atoms: atoms?,
            bonds: bonds?,
            bond_roles,
            head_crossing_bonds: head_crossing_bonds?,
            crossing_bond_correspondence: crossing_bond_correspondence?,
            parent_atoms: parent_atoms?,
            parent,
            label: self.label.clone(),
            connection: self.connection.clone(),
            subtype: self.subtype.clone(),
            bracket_style: self.bracket_style.clone(),
            expansion_state: self.expansion_state.clone(),
            class: self.class.clone(),
            component_number: self.component_number,
            display: self.display.clone(),
            data: self.data.clone(),
            attach_points: attach_points?,
            cstates: cstates?,
            props: self.props.clone(),
            data_fields: self.data_fields.clone(),
        })
    }
}

/// Relationship between the configurations stored in an enhanced stereo group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StereoGroupKind {
    Absolute,
    Or,
    And,
}

/// Detached enhanced-stereo membership and source ID metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    id: Option<u32>,
    kind: StereoGroupKind,
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
}

impl StereoGroup {
    #[must_use]
    pub fn new(kind: StereoGroupKind, atoms: Vec<AtomId>, bonds: Vec<BondId>) -> Self {
        Self {
            id: None,
            kind,
            atoms,
            bonds,
        }
    }

    #[must_use]
    pub const fn with_id(mut self, id: u32) -> Self {
        self.id = Some(id);
        self
    }

    #[must_use]
    pub const fn id(&self) -> Option<u32> {
        self.id
    }

    #[must_use]
    pub const fn kind(&self) -> StereoGroupKind {
        self.kind
    }

    #[must_use]
    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }

    #[must_use]
    pub fn bonds(&self) -> &[BondId] {
        &self.bonds
    }

    pub fn push_atom(&mut self, atom: AtomId) {
        self.atoms.push(atom);
    }

    pub fn push_bond(&mut self, bond: BondId) {
        self.bonds.push(bond);
    }

    pub fn remove_atom(&mut self, atom: AtomId) {
        // RDKit✔️✔️: auto atomPos = findAtom(group);
        // RDKit✔️✔️: if (atomPos != group.d_atoms.end()) {
        // RDKit✔️✔️:   group.d_atoms.erase(atomPos);
        // RDKit✔️✔️: }
        remove_first(&mut self.atoms, &atom);
    }

    pub fn remove_bond(&mut self, bond: BondId) {
        // RDKit✔️✔️: auto bondPos = findBond(group);
        // RDKit✔️✔️: if (bondPos != group.d_bonds.end()) {
        // RDKit✔️✔️:   group.d_bonds.erase(bondPos);
        // RDKit✔️✔️: }
        remove_first(&mut self.bonds, &bond);
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        // RDKit✔️✔️: return gp.getAtoms().empty() &&
        // RDKit✔️✔️:        gp.getBonds().empty();
        self.atoms.is_empty() && self.bonds.is_empty()
    }

    /// Returns a fully remapped group, or `None` if any membership is lost.
    #[must_use]
    pub fn remapped(
        &self,
        atom_map: &[Option<AtomId>],
        bond_map: &[Option<BondId>],
    ) -> Option<Self> {
        let atoms: Option<Vec<_>> = self
            .atoms
            .iter()
            .map(|atom| atom_map.get(atom.index()).and_then(|mapped| *mapped))
            .collect();
        let bonds: Option<Vec<_>> = self
            .bonds
            .iter()
            .map(|bond| bond_map.get(bond.index()).and_then(|mapped| *mapped))
            .collect();
        Some(Self {
            id: self.id,
            kind: self.kind,
            atoms: atoms?,
            bonds: bonds?,
        })
    }
}
