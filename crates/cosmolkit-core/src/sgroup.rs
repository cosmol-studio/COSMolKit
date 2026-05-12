use std::collections::BTreeMap;

use crate::{AtomId, BondId};

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
    pub const fn index(self) -> usize {
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
    pub p1: [f64; 2],
    pub p2: [f64; 2],
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SGroupCState {
    pub bond: BondId,
    pub vector: [f64; 2],
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct SGroupDisplay {
    pub brackets: Vec<SGroupBracket>,
    pub field_position: Option<[f64; 2]>,
    pub display_tag: Option<String>,
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

    pub(crate) fn push_data_field(&mut self, value: impl Into<String>) {
        self.data_fields.push(value.into());
    }

    pub(crate) fn push_atom(&mut self, atom: AtomId) {
        self.atoms.push(atom);
    }

    pub(crate) fn push_bond(&mut self, bond: BondId) {
        self.push_bond_with_role(bond, SGroupBondRole::Crossing);
    }

    pub(crate) fn push_bond_with_role(&mut self, bond: BondId, role: SGroupBondRole) {
        self.bonds.push(bond);
        if role != SGroupBondRole::Crossing {
            self.bond_roles.insert(bond, role);
        }
    }

    pub(crate) fn remove_atom(&mut self, atom: AtomId) {
        self.atoms.retain(|candidate| *candidate != atom);
    }

    pub(crate) fn remove_bond(&mut self, bond: BondId) {
        self.bonds.retain(|candidate| *candidate != bond);
        self.bond_roles.remove(&bond);
    }

    pub(crate) fn remove_parent_atom(&mut self, atom: AtomId) {
        self.parent_atoms.retain(|candidate| *candidate != atom);
    }

    pub(crate) fn clear_attach_point_leaving_atom(&mut self, atom: AtomId) {
        for attach_point in &mut self.attach_points {
            if attach_point.leaving_atom == Some(atom) {
                attach_point.leaving_atom = None;
            }
        }
    }

    #[allow(dead_code)]
    pub(crate) fn push_parent_atom(&mut self, atom: AtomId) {
        self.parent_atoms.push(atom);
    }

    pub(crate) fn set_external_id(&mut self, external_id: u32) {
        self.external_id = Some(external_id);
    }

    #[allow(dead_code)]
    pub(crate) fn set_parent(&mut self, parent: SubstanceGroupId) {
        self.parent = Some(parent);
    }

    #[allow(dead_code)]
    pub(crate) fn set_rdkit_sequence_id(&mut self, rdkit_sequence_id: u32) {
        self.rdkit_sequence_id = Some(rdkit_sequence_id);
    }

    #[allow(dead_code)]
    pub(crate) fn set_label(&mut self, label: impl Into<String>) {
        self.label = Some(label.into());
    }

    #[allow(dead_code)]
    pub(crate) fn set_connection(&mut self, connection: SGroupConnection) {
        self.connection = Some(connection);
    }

    #[allow(dead_code)]
    pub(crate) fn set_subtype(&mut self, subtype: impl Into<String>) {
        self.subtype = Some(subtype.into());
    }

    #[allow(dead_code)]
    pub(crate) fn set_bracket_style(&mut self, bracket_style: SGroupBracketStyle) {
        self.bracket_style = Some(bracket_style);
    }

    #[allow(dead_code)]
    pub(crate) fn set_expansion_state(&mut self, expansion_state: impl Into<String>) {
        self.expansion_state = Some(expansion_state.into());
    }

    #[allow(dead_code)]
    pub(crate) fn set_class(&mut self, class: impl Into<String>) {
        self.class = Some(class.into());
    }

    #[allow(dead_code)]
    pub(crate) fn set_component_number(&mut self, component_number: u32) {
        self.component_number = Some(component_number);
    }

    #[allow(dead_code)]
    pub(crate) fn display_mut(&mut self) -> &mut SGroupDisplay {
        self.display.get_or_insert_with(SGroupDisplay::default)
    }

    #[allow(dead_code)]
    pub(crate) fn data_mut(&mut self) -> &mut SGroupData {
        self.data.get_or_insert_with(SGroupData::default)
    }

    #[allow(dead_code)]
    pub(crate) fn push_attach_point(&mut self, attach_point: SGroupAttachPoint) {
        self.attach_points.push(attach_point);
    }

    #[allow(dead_code)]
    pub(crate) fn push_cstate(&mut self, cstate: SGroupCState) {
        self.cstates.push(cstate);
    }

    pub(crate) fn set_prop(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.props.insert(key.into(), value.into());
    }

    pub(crate) fn can_remap_without_parent(
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

    pub(crate) fn remapped(
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
