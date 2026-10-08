use crate::{MoleculePropertyError, PropertyText, PropertyValue};
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
    Generic(PropertyText),
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
    pub display_tag: Option<PropertyText>,
}

impl SGroupDisplay {
    #[must_use]
    pub fn brackets(&self) -> &[SGroupBracket] {
        &self.brackets
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SGroupData {
    pub field_name: Option<PropertyText>,
    pub field_type: Option<PropertyText>,
    pub field_info: Option<PropertyText>,
    pub field_display: Option<PropertyText>,
    pub units: Option<PropertyText>,
    pub query_type: Option<PropertyText>,
    pub query_op: Option<PropertyText>,
    pub values: Vec<PropertyText>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SGroupConnection {
    HeadToHead,
    HeadToTail,
    Either,
    Unknown(PropertyText),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SGroupBracketStyle {
    Bracket,
    Parenthesis,
    None,
    Unknown(PropertyText),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SGroupAttachPoint {
    pub atom: AtomId,
    pub leaving_atom: Option<AtomId>,
    pub label: Option<PropertyText>,
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
    label: Option<PropertyText>,
    connection: Option<SGroupConnection>,
    subtype: Option<PropertyText>,
    bracket_style: Option<SGroupBracketStyle>,
    expansion_state: Option<PropertyText>,
    class: Option<PropertyText>,
    component_number: Option<u32>,
    display: Option<SGroupDisplay>,
    data: Option<SGroupData>,
    attach_points: Vec<SGroupAttachPoint>,
    cstates: Vec<SGroupCState>,
    props: crate::property_value::PropertyStore,
    data_fields: Vec<PropertyText>,
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
            props: crate::property_value::PropertyStore::new(),
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

    /// Actual sparse role entry; distinguishes absence from explicit Crossing.
    #[doc(hidden)]
    #[must_use]
    pub fn explicit_bond_role(&self, bond: BondId) -> Option<SGroupBondRole> {
        // Detached read only, O(log R); constructors/setters/remapping keep
        // entries restricted to existing member bonds. No runtime authority.
        self.bond_roles.get(&bond).copied()
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

    /// Explicit role entries without materializing implicit Crossing defaults.
    #[doc(hidden)]
    pub fn stored_bond_roles(&self) -> &BTreeMap<BondId, SGroupBondRole> {
        &self.bond_roles
    }

    #[must_use]
    pub const fn parent(&self) -> Option<SubstanceGroupId> {
        self.parent
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<PropertyText, PropertyValue> {
        self.props.values()
    }

    #[must_use]
    pub fn data_fields(&self) -> &[PropertyText] {
        &self.data_fields
    }

    #[must_use]
    pub fn label(&self) -> Option<&PropertyText> {
        self.label.as_ref()
    }

    #[must_use]
    pub const fn connection(&self) -> Option<&SGroupConnection> {
        self.connection.as_ref()
    }

    #[must_use]
    pub fn subtype(&self) -> Option<&PropertyText> {
        self.subtype.as_ref()
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
    pub fn expansion_state(&self) -> Option<&PropertyText> {
        self.expansion_state.as_ref()
    }

    #[must_use]
    pub fn class(&self) -> Option<&PropertyText> {
        self.class.as_ref()
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
    pub fn with_label(mut self, label: impl Into<PropertyText>) -> Self {
        self.label = Some(label.into());
        self
    }

    #[must_use]
    pub fn with_connection(mut self, connection: SGroupConnection) -> Self {
        self.connection = Some(connection);
        self
    }

    #[must_use]
    pub fn with_subtype(mut self, subtype: impl Into<PropertyText>) -> Self {
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
    pub fn with_expansion_state(mut self, expansion_state: impl Into<PropertyText>) -> Self {
        self.expansion_state = Some(expansion_state.into());
        self
    }

    #[must_use]
    pub fn with_class(mut self, class: impl Into<PropertyText>) -> Self {
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

    pub fn with_prop(
        mut self,
        key: impl Into<PropertyText>,
        value: impl Into<PropertyValue>,
    ) -> Result<Self, MoleculePropertyError> {
        self.set_prop(key, value)?;
        Ok(self)
    }

    #[must_use]
    pub fn with_data_field(mut self, value: impl Into<PropertyText>) -> Self {
        self.data_fields.push(value.into());
        self
    }

    pub fn push_data_field(&mut self, value: impl Into<PropertyText>) {
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
    pub fn set_label(&mut self, label: impl Into<PropertyText>) {
        self.label = Some(label.into());
    }

    #[allow(dead_code)]
    pub fn set_connection(&mut self, connection: SGroupConnection) {
        self.connection = Some(connection);
    }

    #[allow(dead_code)]
    pub fn set_subtype(&mut self, subtype: impl Into<PropertyText>) {
        self.subtype = Some(subtype.into());
    }

    #[allow(dead_code)]
    pub fn set_bracket_style(&mut self, bracket_style: SGroupBracketStyle) {
        self.bracket_style = Some(bracket_style);
    }

    #[allow(dead_code)]
    pub fn set_expansion_state(&mut self, expansion_state: impl Into<PropertyText>) {
        self.expansion_state = Some(expansion_state.into());
    }

    #[allow(dead_code)]
    pub fn set_class(&mut self, class: impl Into<PropertyText>) {
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

    pub fn set_prop(
        &mut self,
        key: impl Into<PropertyText>,
        value: impl Into<PropertyValue>,
    ) -> Result<(), MoleculePropertyError> {
        // BEGIN COMPLETE RDProps::setProp
        // RDKit✔️❌: void setProp(const std::string_view key, T val, bool computed = false) const {
        // RDKit✔️❌:     if(key.empty()) {
        // RDKit✔️❌:       throw ValueErrorException("Cannot set property with empty key");
        // RDKit✔️❌:     }
        // RDKit✔️❌:     if (computed) {
        // RDKit✔️❌:       STR_VECT compLst;
        // RDKit✔️❌:       getPropIfPresent(RDKit::detail::computedPropName, compLst);
        // RDKit✔️❌:       if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
        // RDKit✔️❌:         compLst.emplace_back(key);
        // RDKit✔️❌:         d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:     d_props.setVal(key, val);
        // RDKit✔️❌:   }
        // END COMPLETE RDProps::setProp
        // Source default computed=false; delegate the sole typed store, retain
        // tag/bytes/order on replacement and propagate its empty-key failure.
        // Computed assignment is a separate source capability, not inferred
        // from an underscore key or the value kind. Store cost is tree-backed.
        self.props
            .set(key.into(), value.into())
            .map_err(MoleculePropertyError::from)
    }

    pub fn clear_prop(&mut self, key: impl AsRef<[u8]>) -> Result<(), MoleculePropertyError> {
        // BEGIN COMPLETE RDProps::clearProp
        // RDKit✔️❌: void clearProp(const std::string_view key) const {
        // RDKit✔️❌:     STR_VECT compLst;
        // RDKit✔️❌:     if (getPropIfPresent(RDKit::detail::computedPropName, compLst)) {
        // RDKit✔️❌:       auto svi = std::find(compLst.begin(), compLst.end(), key);
        // RDKit✔️❌:       if (svi != compLst.end()) {
        // RDKit✔️❌:         compLst.erase(svi);
        // RDKit✔️❌:         d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:     d_props.clearVal(key);
        // RDKit✔️❌:   }
        // END COMPLETE RDProps::clearProp
        // Delegate actual computed-list read/removal order and typed failure.
        self.props.clear(key).map_err(MoleculePropertyError::from)
    }

    /// Borrow detached typed properties in source insertion order.
    #[doc(hidden)]
    pub fn property_records(&self) -> impl Iterator<Item = (&PropertyText, &PropertyValue)> + '_ {
        self.props.ordered()
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
        self.bonds.contains(&bond) || self.cstates.iter().any(|cstate| cstate.bond == bond)
    }

    fn adjust_source_removed_atom(&mut self, atom: AtomId) {
        // RDKit❗❌: bool SubstanceGroup::adjustToRemovedAtom(unsigned int atomIdx) {
        // RDKit❗❌:   bool res = false;
        // RDKit❗❌:   for (auto &aid : d_atoms) {
        // RDKit❗❌:     if (aid == atomIdx) {
        // RDKit❗❌:       throw SubstanceGroupException(
        // RDKit❗❌:           "adjustToRemovedAtom() called on SubstanceGroup which contains the "
        // RDKit❗❌:           "atom");
        // RDKit❗❌:     }
        // RDKit❗❌:     if (aid > atomIdx) {
        // RDKit❗❌:       res = true;
        // RDKit❗❌:       --aid;
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // RDKit❗❌:   for (auto &aid : d_patoms) {
        // RDKit❗❌:     if (aid == atomIdx) {
        // RDKit❗❌:       throw SubstanceGroupException(
        // RDKit❗❌:           "adjustToRemovedAtom() called on SubstanceGroup which contains the "
        // RDKit❗❌:           "atom");
        // RDKit❗❌:     }
        // RDKit❗❌:     if (aid > atomIdx) {
        // RDKit❗❌:       res = true;
        // RDKit❗❌:       --aid;
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // RDKit❗❌:   for (auto &ap : d_saps) {
        // RDKit❗❌:     if (ap.aIdx == atomIdx || ap.lvIdx == rdcast<int>(atomIdx)) {
        // RDKit❗❌:       throw SubstanceGroupException(
        // RDKit❗❌:           "adjustToRemovedAtom() called on SubstanceGroup which contains the "
        // RDKit❗❌:           "atom");
        // RDKit❗❌:     }
        // RDKit❗❌:     if (ap.aIdx > atomIdx) {
        // RDKit❗❌:       res = true;
        // RDKit❗❌:       --ap.aIdx;
        // RDKit❗❌:     }
        // RDKit❗❌:     if (ap.lvIdx > rdcast<int>(atomIdx)) {
        // RDKit❗❌:       res = true;
        // RDKit❗❌:       --ap.lvIdx;
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // RDKit❗❌:
        // RDKit❗❌:   return res;
        // RDKit❗❌: }
        // RDKit❗❌:
        for id in self.atoms.iter_mut().chain(&mut self.parent_atoms) {
            if id.index() > atom.index() {
                *id = AtomId::new(id.index() - 1);
            }
        }
        for point in &mut self.attach_points {
            if point.atom.index() > atom.index() {
                point.atom = AtomId::new(point.atom.index() - 1);
            }
            if let Some(leaving) = &mut point.leaving_atom {
                if leaving.index() > atom.index() {
                    *leaving = AtomId::new(leaving.index() - 1);
                }
            }
        }
    }

    fn adjust_source_removed_bond(&mut self, bond: BondId) {
        // RDKit❗❌: bool SubstanceGroup::adjustToRemovedBond(unsigned int bondIdx) {
        // RDKit❗❌:   bool res = false;
        // RDKit❗❌:   for (auto &bid : d_bonds) {
        // RDKit❗❌:     if (bid == bondIdx) {
        // RDKit❗❌:       throw SubstanceGroupException(
        // RDKit❗❌:           "adjustToRemovedBond() called on SubstanceGroup which contains the "
        // RDKit❗❌:           "bond");
        // RDKit❗❌:     }
        // RDKit❗❌:     if (bid > bondIdx) {
        // RDKit❗❌:       res = true;
        // RDKit❗❌:       --bid;
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // RDKit❗❌:   for (auto &cs : d_cstates) {
        // RDKit❗❌:     if (cs.bondIdx == bondIdx) {
        // RDKit❗❌:       throw SubstanceGroupException(
        // RDKit❗❌:           "adjustToRemovedBond() called on SubstanceGroup which contains the "
        // RDKit❗❌:           "bond");
        // RDKit❗❌:     }
        // RDKit❗❌:     if (cs.bondIdx > bondIdx) {
        // RDKit❗❌:       res = true;
        // RDKit❗❌:       --cs.bondIdx;
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // RDKit❗❌:   return res;
        // RDKit❗❌: }
        // Caller excludes every matching bond and CState first. Preserve raw
        // XBHEAD/XBCORR and roles: the native helper does not adjust those.
        for bid in &mut self.bonds {
            if bid.index() > bond.index() {
                *bid = BondId::new(bid.index() - 1);
            }
        }
        for cs in &mut self.cstates {
            if cs.bond.index() > bond.index() {
                cs.bond = BondId::new(cs.bond.index() - 1);
            }
        }
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

    /// Source RWMol insertion offsets only explicit graph references. Raw
    /// hierarchy/XBHEAD/XBCORR properties retain their source values.
    #[doc(hidden)]
    pub fn with_inserted_offsets(
        mut self,
        id: SubstanceGroupId,
        atom_offset: usize,
        bond_offset: usize,
    ) -> Self {
        // RDKit❗❌: void insertSubstanceGroups(RWMol &mol, const RWMol &other,
        // RDKit❗❌:                            unsigned int origNumAtoms,
        // RDKit❗❌:                            unsigned int origNumBonds) {
        // RDKit❗❌:   for (auto sgroup : getSubstanceGroups(other)) {
        // RDKit❗❌:     sgroup.setOwningMol(&mol);
        // RDKit❗❌:
        // RDKit❗❌:     // update the sgroup's atom and bond indices
        // RDKit❗❌:     auto atom_indices = sgroup.getAtoms();
        // RDKit❗❌:     std::transform(atom_indices.begin(), atom_indices.end(),
        // RDKit❗❌:                    atom_indices.begin(),
        // RDKit❗❌:                    [&origNumAtoms](unsigned int old_index) {
        // RDKit❗❌:                      return origNumAtoms + old_index;
        // RDKit❗❌:                    });
        // RDKit❗❌:     sgroup.setAtoms(atom_indices);
        // RDKit❗❌:
        // RDKit❗❌:     auto bond_indices = sgroup.getBonds();
        // RDKit❗❌:     std::transform(bond_indices.begin(), bond_indices.end(),
        // RDKit❗❌:                    bond_indices.begin(),
        // RDKit❗❌:                    [&origNumBonds](unsigned int old_index) {
        // RDKit❗❌:                      return origNumBonds + old_index;
        // RDKit❗❌:                    });
        // RDKit❗❌:     sgroup.setBonds(bond_indices);
        // RDKit❗❌:
        // RDKit❗❌:     // patoms
        // RDKit❗❌:     auto patom_indices = sgroup.getParentAtoms();
        // RDKit❗❌:     std::transform(patom_indices.begin(), patom_indices.end(),
        // RDKit❗❌:                    patom_indices.begin(),
        // RDKit❗❌:                    [&origNumAtoms](unsigned int old_index) {
        // RDKit❗❌:                      return origNumAtoms + old_index;
        // RDKit❗❌:                    });
        // RDKit❗❌:     sgroup.setParentAtoms(patom_indices);
        // RDKit❗❌:
        // RDKit❗❌:     // cstates (these are references, can be updated in place)
        // RDKit❗❌:     for (auto &cstate : sgroup.getCStates()) {
        // RDKit❗❌:       cstate.bondIdx = origNumBonds + cstate.bondIdx;
        // RDKit❗❌:     }
        // RDKit❗❌:
        // RDKit❗❌:     // attachment points (can also be updated in place)
        // RDKit❗❌:     for (auto &sap : sgroup.getAttachPoints()) {
        // RDKit❗❌:       sap.aIdx = origNumAtoms + sap.aIdx;
        // RDKit❗❌:       if (sap.lvIdx != -1) {
        // RDKit❗❌:         sap.lvIdx = static_cast<int>(origNumAtoms + sap.lvIdx);
        // RDKit❗❌:       }
        // RDKit❗❌:     }
        // RDKit❗❌:
        // RDKit❗❌:     addSubstanceGroup(mol, sgroup);
        // RDKit❗❌:   }
        // RDKit❗❌: }
        // One consumed group is the loop element of the canonical append
        // caller, which preserves source collection order and owns the clone.
        // Native setters validate shifted references immediately; this model
        // adapter retains later detached graph validation and usize arithmetic.
        // Unsigned32 wrap, signed SAP negatives other than sentinel -1, owning
        // pointer/index errors and setter failure timing remain explicit gaps.
        // Parent hierarchy and generic XBHEAD/XBCORR entries are not shifted:
        // the native function changes only the reference vectors shown above.
        self.id = id;
        for atom in self.atoms.iter_mut().chain(self.parent_atoms.iter_mut()) {
            *atom = AtomId::new(atom.index() + atom_offset);
        }
        for bond in &mut self.bonds {
            *bond = BondId::new(bond.index() + bond_offset);
        }
        // Ordered bond-role index reconstruction costs O(R log R); SOURCE
        // directly offsets linear member vectors. This extra index is retained
        // as canonical detached metadata, without claiming cost equivalence.
        self.bond_roles = self
            .bond_roles
            .into_iter()
            .map(|(bond, role)| (BondId::new(bond.index() + bond_offset), role))
            .collect();
        for state in &mut self.cstates {
            state.bond = BondId::new(state.bond.index() + bond_offset);
        }
        for point in &mut self.attach_points {
            point.atom = AtomId::new(point.atom.index() + atom_offset);
            if let Some(atom) = point.leaving_atom {
                point.leaving_atom = Some(AtomId::new(atom.index() + atom_offset));
            }
        }
        // Explicit-reference loops remain linear. The extra ordered role
        // index cost is stated above; unrelated generic properties stay raw.
        self
    }

    pub fn remapped(
        &self,
        id: SubstanceGroupId,
        atom_map: &[Option<AtomId>],
        bond_map: &[Option<BondId>],
        sgroup_map: &[Option<SubstanceGroupId>],
    ) -> Option<Self> {
        // BEGIN RDKIT CPP FUNCTION SubstanceGroup::adjustToRemovedAtom
        // RDKit❗✔️: bool SubstanceGroup::adjustToRemovedAtom(unsigned int atomIdx) {
        // RDKit❗✔️:   bool res = false;
        // RDKit❗✔️:   for (auto &aid : d_atoms) {
        // RDKit❗✔️:     if (aid == atomIdx) {
        // RDKit❗✔️:       throw SubstanceGroupException(
        // RDKit❗✔️:           "adjustToRemovedAtom() called on SubstanceGroup which contains the "
        // RDKit❗✔️:           "atom");
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (aid > atomIdx) {
        // RDKit❗✔️:       res = true;
        // RDKit❗✔️:       --aid;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   for (auto &aid : d_patoms) {
        // RDKit❗✔️:     if (aid == atomIdx) {
        // RDKit❗✔️:       throw SubstanceGroupException(
        // RDKit❗✔️:           "adjustToRemovedAtom() called on SubstanceGroup which contains the "
        // RDKit❗✔️:           "atom");
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (aid > atomIdx) {
        // RDKit❗✔️:       res = true;
        // RDKit❗✔️:       --aid;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   for (auto &ap : d_saps) {
        // RDKit❗✔️:     if (ap.aIdx == atomIdx || ap.lvIdx == rdcast<int>(atomIdx)) {
        // RDKit❗✔️:       throw SubstanceGroupException(
        // RDKit❗✔️:           "adjustToRemovedAtom() called on SubstanceGroup which contains the "
        // RDKit❗✔️:           "atom");
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (ap.aIdx > atomIdx) {
        // RDKit❗✔️:       res = true;
        // RDKit❗✔️:       --ap.aIdx;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (ap.lvIdx > rdcast<int>(atomIdx)) {
        // RDKit❗✔️:       res = true;
        // RDKit❗✔️:       --ap.lvIdx;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION
        // BEGIN RDKIT CPP FUNCTION SubstanceGroup::adjustToRemovedBond
        // RDKit❗✔️: bool SubstanceGroup::adjustToRemovedBond(unsigned int bondIdx) {
        // RDKit❗✔️:   bool res = false;
        // RDKit❗✔️:   for (auto &bid : d_bonds) {
        // RDKit❗✔️:     if (bid == bondIdx) {
        // RDKit❗✔️:       throw SubstanceGroupException(
        // RDKit❗✔️:           "adjustToRemovedBond() called on SubstanceGroup which contains the "
        // RDKit❗✔️:           "bond");
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (bid > bondIdx) {
        // RDKit❗✔️:       res = true;
        // RDKit❗✔️:       --bid;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   for (auto &cs : d_cstates) {
        // RDKit❗✔️:     if (cs.bondIdx == bondIdx) {
        // RDKit❗✔️:       throw SubstanceGroupException(
        // RDKit❗✔️:           "adjustToRemovedBond() called on SubstanceGroup which contains the "
        // RDKit❗✔️:           "bond");
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (cs.bondIdx > bondIdx) {
        // RDKit❗✔️:       res = true;
        // RDKit❗✔️:       --cs.bondIdx;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION
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

/// Detached enhanced-stereo membership with independent read and write IDs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    /// Source/read identifier, when present in the input representation.
    id: Option<u32>,
    /// Output/write identifier assigned by the source writer rules.
    write_id: u32,
    kind: StereoGroupKind,
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
}

impl StereoGroup {
    #[must_use]
    pub fn new(kind: StereoGroupKind, atoms: Vec<AtomId>, bonds: Vec<BondId>) -> Self {
        Self {
            // RDKit❗✔️: d_readId{readId} {}
            // Keep the model's absent-ID state distinct from explicit zero.
            id: None,
            // RDKit✔️✔️: unsigned d_writeId = 0u;
            write_id: 0,
            kind,
            atoms,
            bonds,
        }
    }

    /// Sets the source/read identifier without changing the output/write ID.
    #[must_use]
    pub const fn with_id(mut self, id: u32) -> Self {
        // RDKit✔️✔️: d_readId{readId} {}
        self.id = Some(id);
        self
    }

    /// Returns the source/read identifier, when one was represented.
    #[must_use]
    pub const fn id(&self) -> Option<u32> {
        self.id
    }

    /// Returns the output/write identifier.
    #[doc(hidden)]
    #[must_use]
    pub const fn write_id(&self) -> u32 {
        // RDKit✔️✔️: unsigned getWriteId() const { return d_writeId; }
        self.write_id
    }

    /// Sets the output/write identifier while retaining the source/read ID.
    #[doc(hidden)]
    #[must_use]
    pub const fn with_write_id(mut self, write_id: u32) -> Self {
        // RDKit❗✔️: void setWriteId(unsigned id) { d_writeId = id; }
        // This detached builder consumes and returns the value instead of
        // mutating a borrowed StereoGroup.
        self.write_id = write_id;
        self
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

    pub(crate) fn shift_source_atom_ids_after_removed(&mut self, removed: AtomId) {
        for atom in &mut self.atoms {
            if atom.index() > removed.index() {
                *atom = AtomId::new(atom.index() - 1);
            }
        }
    }

    pub(crate) fn remap_source_retained_atom_ids(&mut self, aliases: &[Option<AtomId>]) {
        for atom in &mut self.atoms {
            if let Some(Some(mapped)) = aliases.get(atom.index()) {
                *atom = *mapped;
            }
        }
    }

    pub(crate) fn remap_source_retained_bond_ids(&mut self, aliases: &[Option<BondId>]) {
        for bond in &mut self.bonds {
            if let Some(Some(mapped)) = aliases.get(bond.index()) {
                *bond = *mapped;
            }
        }
    }

    pub(crate) fn shift_source_bond_ids_after_removed(&mut self, removed: BondId) {
        // Native StereoGroup retains Bond pointers while RWMol changes indices.
        // Stable detached IDs must track those same retained objects explicitly.
        for bond in &mut self.bonds {
            if bond.index() > removed.index() {
                *bond = BondId::new(bond.index() - 1);
            }
        }
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
        // BEGIN RDKIT CPP FUNCTION Subset.cpp::copySelectedStereoGroups
        // RDKit✔️✔️: extracted_stereo_groups.push_back({stereo_group.getGroupType(),
        // RDKit✔️✔️:                                    std::move(atoms), std::move(bonds),
        // RDKit✔️✔️:                                    stereo_group.getReadId()});
        // RDKit✔️✔️: extracted_stereo_groups.back().setWriteId(stereo_group.getWriteId());
        // END RDKIT CPP FUNCTION Subset.cpp::copySelectedStereoGroups
        Some(Self {
            id: self.id,
            // RDKit✔️✔️: extracted_stereo_groups.back().setWriteId(stereo_group.getWriteId());
            write_id: self.write_id,
            kind: self.kind,
            atoms: atoms?,
            bonds: bonds?,
        })
    }
}

/// Returns the detached enhanced-stereo output ID used by domain writers.
#[doc(hidden)]
#[must_use]
pub const fn stereo_group_write_id(group: &StereoGroup) -> u32 {
    // RDKit✔️✔️: unsigned getWriteId() const { return d_writeId; }
    group.write_id
}

/// Sets the detached enhanced-stereo output ID used by domain writers.
#[doc(hidden)]
pub fn set_stereo_group_write_id(group: &mut StereoGroup, write_id: u32) {
    // RDKit✔️✔️: void setWriteId(unsigned id) { d_writeId = id; }
    group.write_id = write_id;
}

#[cfg(test)]
mod stereo_write_id_tests {
    use super::*;

    #[test]
    fn stereo_write_id_default_and_explicit_values_are_independent_from_read_id() {
        let default_group = StereoGroup::new(StereoGroupKind::Or, vec![AtomId::new(0)], vec![]);
        assert_eq!(stereo_group_write_id(&default_group), 0);
        assert_eq!(default_group.id(), None);

        let mut group = default_group.with_id(7);
        assert_eq!(group.id(), Some(7));
        assert_eq!(stereo_group_write_id(&group), 0);

        set_stereo_group_write_id(&mut group, 0);
        assert_eq!(group.id(), Some(7));
        assert_eq!(stereo_group_write_id(&group), 0);

        set_stereo_group_write_id(&mut group, 19);
        assert_eq!(group.id(), Some(7));
        assert_eq!(stereo_group_write_id(&group), 19);
    }

    #[test]
    fn stereo_write_id_clone_equality_and_remap_preserve_stored_state() {
        let mut group = StereoGroup::new(
            StereoGroupKind::And,
            vec![AtomId::new(0)],
            vec![BondId::new(0)],
        )
        .with_id(7);
        set_stereo_group_write_id(&mut group, 19);

        let cloned = group.clone();
        assert_eq!(cloned, group);
        assert_eq!(stereo_group_write_id(&cloned), 19);

        let mut different_write_state = cloned.clone();
        set_stereo_group_write_id(&mut different_write_state, 20);
        assert_ne!(different_write_state, group);
        assert_eq!(different_write_state.id(), group.id());

        let remapped = group
            .remapped(&[Some(AtomId::new(3))], &[Some(BondId::new(4))])
            .expect("all group members have mappings");
        assert_eq!(remapped.id(), Some(7));
        assert_eq!(stereo_group_write_id(&remapped), 19);
        assert_eq!(remapped.atoms(), &[AtomId::new(3)]);
        assert_eq!(remapped.bonds(), &[BondId::new(4)]);

        assert_eq!(group.remapped(&[None], &[Some(BondId::new(4))]), None);
        assert_eq!(group.remapped(&[Some(AtomId::new(3))], &[None]), None);
    }
}

#[cfg(test)]
mod cf3d_sgids_model_1_tests {
    use super::{StereoGroup, StereoGroupKind};
    use crate::{AtomId, BondId};

    #[test]
    fn cf3d_sgids_model_1_default_and_setters_keep_identity_axes_independent() {
        let default = StereoGroup::new(StereoGroupKind::Or, vec![], vec![]);
        assert_eq!(default.id(), None);
        assert_eq!(default.write_id(), 0);

        let read = default.clone().with_id(7);
        assert_eq!(read.id(), Some(7));
        assert_eq!(read.write_id(), 0);

        let both = read.clone().with_write_id(9);
        assert_eq!(both.id(), Some(7));
        assert_eq!(both.write_id(), 9);
        assert_eq!(both, read.clone().with_write_id(9));
        assert_ne!(both, read);

        let changed_read = both.clone().with_id(0);
        assert_eq!(changed_read.id(), Some(0));
        assert_eq!(changed_read.write_id(), 9);
        assert_ne!(default, default.clone().with_id(0));
    }

    #[test]
    fn cf3d_sgids_model_1_all_read_write_combinations_clone_and_remap() {
        let atom_map = [Some(AtomId::new(4)), Some(AtomId::new(2))];
        let bond_map = [Some(BondId::new(3))];
        for read_id in [None, Some(0), Some(7)] {
            for write_id in [0, 9] {
                let mut source = StereoGroup::new(
                    StereoGroupKind::And,
                    vec![AtomId::new(0), AtomId::new(1)],
                    vec![BondId::new(0)],
                );
                if let Some(read_id) = read_id {
                    source = source.with_id(read_id);
                }
                source = source.with_write_id(write_id);
                let unchanged = source.clone();

                let cloned = source.clone();
                assert_eq!(cloned.id(), read_id);
                assert_eq!(cloned.write_id(), write_id);
                assert_eq!(cloned, source);

                let remapped = source.remapped(&atom_map, &bond_map).unwrap();
                assert_eq!(remapped.id(), read_id);
                assert_eq!(remapped.write_id(), write_id);
                assert_eq!(remapped.atoms(), &[AtomId::new(4), AtomId::new(2)]);
                assert_eq!(remapped.bonds(), &[BondId::new(3)]);
                assert_eq!(remapped.kind(), StereoGroupKind::And);
                assert_eq!(source, unchanged);
            }
        }
    }
}

/// Source assignment merges multiple ABS groups, retaining non-ABS order and
/// source's reverse concatenation of ABS members without sorting/deduplication.
#[doc(hidden)]
pub fn merge_absolute_stereo_groups(groups: Vec<StereoGroup>) -> Vec<StereoGroup> {
    // RDKit❗✔️: void ROMol::setStereoGroups(std::vector<StereoGroup> stereo_groups) {
    // RDKit❗✔️:   auto is_abs = [](const auto &sg) {
    // RDKit❗✔️:     return sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE;
    // RDKit❗✔️:   };
    // RDKit❗✔️:
    // RDKit❗✔️:   // if there's more than one ABS group, merge them
    // RDKit❗✔️:   if (auto num_abs = std::ranges::count_if(stereo_groups, is_abs);
    // RDKit❗✔️:       num_abs <= 1) {
    // RDKit❗✔️:     d_stereo_groups = std::move(stereo_groups);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     std::vector<Atom *> abs_atoms;
    // RDKit❗✔️:     std::vector<Bond *> abs_bonds;
    // RDKit❗✔️:     std::vector<StereoGroup> new_stereo_groups;
    // RDKit❗✔️:     new_stereo_groups.reserve(stereo_groups.size() - num_abs + 1);
    // RDKit❗✔️:     for (auto &&sg : stereo_groups) {
    // RDKit❗✔️:       if (is_abs(sg)) {
    // RDKit❗✔️:         auto &other_atoms = sg.getAtoms();
    // RDKit❗✔️:         auto &other_bonds = sg.getBonds();
    // RDKit❗✔️:         abs_atoms.insert(abs_atoms.begin(), other_atoms.begin(),
    // RDKit❗✔️:                          other_atoms.end());
    // RDKit❗✔️:         abs_bonds.insert(abs_bonds.begin(), other_bonds.begin(),
    // RDKit❗✔️:                          other_bonds.end());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         new_stereo_groups.push_back(std::move(sg));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     new_stereo_groups.emplace_back(StereoGroupType::STEREO_ABSOLUTE,
    // RDKit❗✔️:                                    std::move(abs_atoms), std::move(abs_bonds));
    // RDKit❗✔️:     d_stereo_groups = std::move(new_stereo_groups);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Same two linear group scans and source prepend-vector insertion cost.
    // Owned groups preserve non-ABS IDs/bonds without extra deep clones.
    let count = groups
        .iter()
        .filter(|g| g.kind == StereoGroupKind::Absolute)
        .count();
    if count <= 1 {
        return groups;
    }
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut result = Vec::with_capacity(groups.len() - count + 1);
    for group in groups {
        if group.kind == StereoGroupKind::Absolute {
            atoms.splice(0..0, group.atoms);
            bonds.splice(0..0, group.bonds);
        } else {
            result.push(group);
        }
    }
    result.push(StereoGroup::new(StereoGroupKind::Absolute, atoms, bonds));
    result
}

/// Source graph insertion of enhanced groups, including ordered ABS merging.
#[doc(hidden)]
pub fn insert_stereo_groups(
    existing: &[StereoGroup],
    incoming: &[StereoGroup],
    atom_offset: usize,
    bond_offset: usize,
) -> Vec<StereoGroup> {
    // RDKit❗❌: void insertStereoGroups(RWMol &mol, const ROMol &other,
    // RDKit❗❌:                         unsigned int origNumAtoms, unsigned int origNumBonds) {
    // RDKit❗❌:   if (other.getStereoGroups().empty()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<RDKit::Atom *> abs_atoms;
    // RDKit❗❌:   std::vector<RDKit::Bond *> abs_bonds;
    // RDKit❗❌:   std::vector<RDKit::StereoGroup> new_groups;
    // RDKit❗❌:   new_groups.reserve(mol.getStereoGroups().size());
    // RDKit❗❌:   for (const auto &sg : mol.getStereoGroups()) {
    // RDKit❗❌:     // The sdf specification forbids more than one ABS stereo group, but we
    // RDKit❗❌:     // don't enforce that in our code. But if we see more than one ABS groups
    // RDKit❗❌:     // here, just merge the atoms and bonds in them into one group. Other stereo
    // RDKit❗❌:     // groups are just forwarded.
    // RDKit❗❌:     if (sg.getGroupType() == RDKit::StereoGroupType::STEREO_ABSOLUTE) {
    // RDKit❗❌:       auto &atoms = sg.getAtoms();
    // RDKit❗❌:       auto &bonds = sg.getBonds();
    // RDKit❗❌:       abs_atoms.insert(abs_atoms.end(), atoms.begin(), atoms.end());
    // RDKit❗❌:       abs_bonds.insert(abs_bonds.end(), bonds.begin(), bonds.end());
    // RDKit❗❌:     } else {
    // RDKit❗❌:       new_groups.emplace_back(sg);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto &sg : other.getStereoGroups()) {
    // RDKit❗❌:     // update the stereo group's atom and bond indices
    // RDKit❗❌:     std::vector<RDKit::Atom *> new_atoms;
    // RDKit❗❌:     std::vector<RDKit::Bond *> new_bonds;
    // RDKit❗❌:     for (auto atom : sg.getAtoms()) {
    // RDKit❗❌:       auto idx = atom->getIdx() + origNumAtoms;
    // RDKit❗❌:       new_atoms.push_back(mol.getAtomWithIdx(idx));
    // RDKit❗❌:     }
    // RDKit❗❌:     for (auto bond : sg.getBonds()) {
    // RDKit❗❌:       auto idx = bond->getIdx() + origNumBonds;
    // RDKit❗❌:       new_bonds.push_back(mol.getBondWithIdx(idx));
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // Collect all ABS atoms and bonds so they are added as a single group
    // RDKit❗❌:     if (sg.getGroupType() == RDKit::StereoGroupType::STEREO_ABSOLUTE) {
    // RDKit❗❌:       abs_atoms.insert(abs_atoms.end(), new_atoms.begin(), new_atoms.end());
    // RDKit❗❌:       abs_bonds.insert(abs_bonds.end(), new_bonds.begin(), new_bonds.end());
    // RDKit❗❌:     } else {
    // RDKit❗❌:       RDKit::StereoGroup new_group(sg.getGroupType(), new_atoms, new_bonds,
    // RDKit❗❌:                                    sg.getReadId());
    // RDKit❗❌:       // default write ID to 0 to avoid id clashes. We can use
    // RDKit❗❌:       // assignStereoGroupIds() later on to assign new IDs
    // RDKit❗❌:       new_group.setWriteId(0);
    // RDKit❗❌:       new_groups.push_back(std::move(new_group));
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!abs_atoms.empty() || !abs_bonds.empty()) {
    // RDKit❗❌:     new_groups.emplace_back(RDKit::StereoGroupType::STEREO_ABSOLUTE, abs_atoms,
    // RDKit❗❌:                             abs_bonds);
    // RDKit❗❌:   }
    // RDKit❗❌:   mol.setStereoGroups(new_groups);
    // RDKit❗❌: }
    // Source group order, duplicate member multiplicity, forward ABS append,
    // existing read/write IDs and incoming write-ID reset to zero are retained.
    // Stable detached IDs replace owning-pointer resolution; usize offsets do
    // not reproduce unsigned32 overflow, retained for final source-width review.
    // Cost: nonempty work is linear in all groups/members like native copying.
    // Empty incoming still deep-copies existing groups to return an owned Vec,
    // whereas native returns immediately; overall cost marker records this
    // existing boundary loss. No sorting/deduplication or heuristic ID repair.
    if incoming.is_empty() {
        return existing.to_vec();
    }
    let mut abs_atoms = Vec::new();
    let mut abs_bonds = Vec::new();
    let mut result = Vec::with_capacity(existing.len());
    for group in existing {
        if group.kind == StereoGroupKind::Absolute {
            abs_atoms.extend_from_slice(&group.atoms);
            abs_bonds.extend_from_slice(&group.bonds);
        } else {
            result.push(group.clone());
        }
    }
    for group in incoming {
        let atoms: Vec<_> = group
            .atoms
            .iter()
            .map(|a| AtomId::new(a.index() + atom_offset))
            .collect();
        let bonds: Vec<_> = group
            .bonds
            .iter()
            .map(|b| BondId::new(b.index() + bond_offset))
            .collect();
        if group.kind == StereoGroupKind::Absolute {
            abs_atoms.extend(atoms);
            abs_bonds.extend(bonds);
        } else {
            let mut group_copy = StereoGroup::new(group.kind, atoms, bonds);
            if let Some(id) = group.id {
                group_copy = group_copy.with_id(id);
            }
            result.push(group_copy);
        }
    }
    if !abs_atoms.is_empty() || !abs_bonds.is_empty() {
        result.push(StereoGroup::new(
            StereoGroupKind::Absolute,
            abs_atoms,
            abs_bonds,
        ));
    }
    merge_absolute_stereo_groups(result)
}

#[cfg(test)]
mod source_insert_stereo_groups_complete_tests {
    use super::*;

    fn group(
        kind: StereoGroupKind,
        atoms: &[usize],
        bonds: &[usize],
        read: u32,
        write: u32,
    ) -> StereoGroup {
        StereoGroup::new(
            kind,
            atoms.iter().copied().map(AtomId::new).collect(),
            bonds.iter().copied().map(BondId::new).collect(),
        )
        .with_id(read)
        .with_write_id(write)
    }

    #[test]
    fn empty_incoming_preserves_multiple_abs_groups_without_assignment_merge() {
        let existing = vec![
            group(StereoGroupKind::Absolute, &[2], &[1], 7, 8),
            group(StereoGroupKind::Absolute, &[0, 0], &[], 9, 10),
        ];
        let result = insert_stereo_groups(&existing, &[], 100, 200);
        assert_eq!(result, existing);
        assert_ne!(result.as_ptr(), existing.as_ptr());
        assert_eq!(result[0].write_id(), 8);
        assert_eq!(result[1].write_id(), 10);
    }

    #[test]
    fn nonempty_insertion_keeps_member_order_duplicates_and_independent_id_rules() {
        let existing = vec![
            group(StereoGroupKind::Absolute, &[2, 2], &[1], 3, 4),
            group(StereoGroupKind::And, &[1], &[0], 5, 6),
            group(StereoGroupKind::Absolute, &[0], &[2], 7, 8),
        ];
        let incoming = vec![
            group(StereoGroupKind::Or, &[1, 0], &[2, 0], 9, 10),
            group(StereoGroupKind::Absolute, &[0, 1, 0], &[1], 11, 12),
            group(StereoGroupKind::And, &[], &[], 13, 14),
        ];
        let result = insert_stereo_groups(&existing, &incoming, 10, 20);
        assert_eq!(result.len(), 4);
        assert_eq!(result[0], existing[1]);
        assert_eq!(result[1].kind(), StereoGroupKind::Or);
        assert_eq!(result[1].id(), Some(9));
        assert_eq!(result[1].write_id(), 0);
        assert_eq!(result[1].atoms(), [AtomId::new(11), AtomId::new(10)]);
        assert_eq!(result[1].bonds(), [BondId::new(22), BondId::new(20)]);
        assert_eq!(result[2].kind(), StereoGroupKind::And);
        assert_eq!(result[2].id(), Some(13));
        assert_eq!(result[2].write_id(), 0);
        assert!(result[2].is_empty());
        assert_eq!(result[3].kind(), StereoGroupKind::Absolute);
        assert_eq!(result[3].id(), None);
        assert_eq!(result[3].write_id(), 0);
        assert_eq!(
            result[3].atoms(),
            [
                AtomId::new(2),
                AtomId::new(2),
                AtomId::new(0),
                AtomId::new(10),
                AtomId::new(11),
                AtomId::new(10)
            ]
        );
        assert_eq!(
            result[3].bonds(),
            [BondId::new(1), BondId::new(2), BondId::new(21)]
        );
        assert_eq!(existing[1].write_id(), 6);
        assert_eq!(incoming[0].write_id(), 10);
    }

    #[test]
    fn nonempty_empty_abs_input_drops_empty_abs_groups_and_keeps_bond_only_abs() {
        let empty = group(StereoGroupKind::Absolute, &[], &[], 5, 6);
        assert!(insert_stereo_groups(&[empty.clone()], &[empty.clone()], 3, 4).is_empty());
        let bond_only = group(StereoGroupKind::Absolute, &[], &[0, 0], 7, 8);
        let result = insert_stereo_groups(&[empty], &[bond_only], 3, 4);
        assert_eq!(result.len(), 1);
        assert!(result[0].atoms().is_empty());
        assert_eq!(result[0].bonds(), [BondId::new(4), BondId::new(4)]);
        assert_eq!(result[0].id(), None);
        assert_eq!(result[0].write_id(), 0);
    }
}

#[cfg(test)]
mod source_insert_substance_groups_complete_tests {
    use super::*;

    #[test]
    fn insertion_shifts_explicit_refs_and_preserves_raw_properties_and_sparse_payloads() {
        let source = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(2), AtomId::new(0), AtomId::new(2)])
            .with_parent_atoms(vec![AtomId::new(1), AtomId::new(1)])
            .with_bonds(vec![BondId::new(1), BondId::new(0), BondId::new(1)])
            .with_bond_role(BondId::new(1), SGroupBondRole::Contained)
            .with_parent(SubstanceGroupId::new(0))
            .with_head_crossing_bonds(vec![BondId::new(0)])
            .with_crossing_bond_correspondence(vec![BondId::new(1), BondId::new(0)])
            .with_cstates(vec![SGroupCState::new(BondId::new(1), [-0.0, 2.0, 3.0])])
            .with_attach_points(vec![
                SGroupAttachPoint {
                    atom: AtomId::new(2),
                    leaving_atom: None,
                    label: Some(PropertyText::from_bytes(b"label\0\xff")),
                    order: Some(u32::MAX),
                },
                SGroupAttachPoint {
                    atom: AtomId::new(0),
                    leaving_atom: Some(AtomId::new(1)),
                    label: None,
                    order: None,
                },
            ])
            .with_rdkit_sequence_id(11)
            .with_external_id(19)
            .with_prop("PARENT", PropertyValue::UInt(2))
            .unwrap()
            .with_prop("XBHEAD", PropertyValue::IntVector(vec![0, 1]))
            .unwrap()
            .with_prop("XBCORR", PropertyValue::IntVector(vec![1, 0]))
            .unwrap();
        let shifted = source
            .clone()
            .with_inserted_offsets(SubstanceGroupId::new(8), 10, 20);
        assert_eq!(shifted.id(), SubstanceGroupId::new(8));
        assert_eq!(
            shifted.atoms(),
            [AtomId::new(12), AtomId::new(10), AtomId::new(12)]
        );
        assert_eq!(shifted.parent_atoms(), [AtomId::new(11), AtomId::new(11)]);
        assert_eq!(
            shifted.bonds(),
            [BondId::new(21), BondId::new(20), BondId::new(21)]
        );
        assert_eq!(
            shifted.bond_roles.get(&BondId::new(21)),
            Some(&SGroupBondRole::Contained)
        );
        assert_eq!(shifted.cstates[0].bond, BondId::new(21));
        assert_eq!(shifted.cstates[0].vector[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(shifted.attach_points[0].atom, AtomId::new(12));
        assert_eq!(shifted.attach_points[0].leaving_atom, None);
        assert_eq!(shifted.attach_points[1].atom, AtomId::new(10));
        assert_eq!(shifted.attach_points[1].leaving_atom, Some(AtomId::new(11)));
        assert_eq!(
            shifted.attach_points[0].label,
            source.attach_points[0].label
        );
        assert_eq!(shifted.attach_points[0].order, Some(u32::MAX));
        assert_eq!(shifted.props, source.props);
        assert_eq!(shifted.parent, source.parent);
        assert_eq!(shifted.head_crossing_bonds, source.head_crossing_bonds);
        assert_eq!(
            shifted.crossing_bond_correspondence,
            source.crossing_bond_correspondence
        );
        assert_eq!(shifted.rdkit_sequence_id, source.rdkit_sequence_id);
        assert_eq!(shifted.external_id, source.external_id);
        assert_eq!(
            source.atoms(),
            [AtomId::new(2), AtomId::new(0), AtomId::new(2)]
        );
    }

    #[test]
    fn empty_and_zero_offset_insertions_keep_payload_and_do_not_assign_index_property() {
        let source = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
            .with_prop("index", PropertyValue::UInt(91))
            .unwrap();
        let shifted = source
            .clone()
            .with_inserted_offsets(SubstanceGroupId::new(4), 11, 17);
        assert_eq!(shifted.id(), SubstanceGroupId::new(4));
        assert_eq!(shifted.props, source.props);
        assert!(shifted.atoms.is_empty());
        assert!(shifted.bonds.is_empty());
        assert!(shifted.cstates.is_empty());
        assert!(shifted.attach_points.is_empty());
        assert_eq!(
            source.clone().with_inserted_offsets(source.id(), 0, 0),
            source
        );
    }
}

/// The single source-shaped deletion helper, borrowing canonical detached groups.
fn remove_source_groups_referencing<E: From<crate::TopologyEditError>>(
    groups: &mut Vec<SubstanceGroup>,
    includes: impl Fn(&SubstanceGroup) -> bool,
    adjust: impl Fn(&mut SubstanceGroup),
    uint_reader: &mut (impl FnMut(&PropertyValue) -> Result<u32, E> + ?Sized),
) -> Result<(), E> {
    // RDKit❗❌: bool removedParentInHierarchy(
    // RDKit❗❌:     unsigned int idx, const std::vector<SubstanceGroup> &sgs,
    // RDKit❗❌:     const boost::dynamic_bitset<> &toRemove,
    // RDKit❗❌:     const std::map<unsigned int, unsigned int> &indexLookup) {
    // RDKit❗❌:   PRECONDITION(idx < sgs.size(), "cannot find SubstanceGroup");
    // RDKit❗❌:   if (toRemove[idx]) {
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int parent;
    // RDKit❗❌:   if (sgs[idx].getPropIfPresent("PARENT", parent)) {
    // RDKit❗❌:     auto piter = indexLookup.find(parent);
    // RDKit❗❌:     if (piter != indexLookup.end()) {
    // RDKit❗❌:       return removedParentInHierarchy(piter->second, sgs, toRemove,
    // RDKit❗❌:                                       indexLookup);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return false;
    // RDKit❗❌: }
    // RDKit❗❌:
    // RDKit❗❌: template <bool INCLUDES_METHOD(SubstanceGroup &, unsigned int),
    // RDKit❗❌:           void ADJUST_METHOD(SubstanceGroup &, unsigned int)>
    // RDKit❗❌: void removeSubstanceGroupsReferencing(RWMol &mol, unsigned int idx) {
    // RDKit❗❌:   auto &sgs = getSubstanceGroups(mol);
    // RDKit❗❌:   if (!sgs.empty()) {
    // RDKit❗❌:     // first collect the ones that should be removed
    // RDKit❗❌:     boost::dynamic_bitset<> toRemove(sgs.size());
    // RDKit❗❌:     unsigned int nRemoved = 0;
    // RDKit❗❌:     bool parentsPresent = false;
    // RDKit❗❌:     for (unsigned int i = 0; i < sgs.size(); ++i) {
    // RDKit❗❌:       if (!parentsPresent && sgs[i].hasProp("PARENT")) {
    // RDKit❗❌:         parentsPresent = true;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (INCLUDES_METHOD(sgs[i], idx)) {
    // RDKit❗❌:         toRemove.set(i);
    // RDKit❗❌:         ++nRemoved;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // if we're going to be removing anything and there are PARENTS present,
    // RDKit❗❌:     // we need to build a lookup map between index->position in original array
    // RDKit❗❌:     std::map<unsigned int, unsigned int> indexLookup;
    // RDKit❗❌:     if (parentsPresent && nRemoved) {
    // RDKit❗❌:       for (unsigned int i = 0; i < sgs.size(); ++i) {
    // RDKit❗❌:         unsigned int index;
    // RDKit❗❌:         if (sgs[i].getPropIfPresent("index", index)) {
    // RDKit❗❌:           indexLookup[index] = i;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     // now go through and keep everything that shouldn't be removed
    // RDKit❗❌:     // and who doesn't have a PARENT that should be removed in their hierarchy
    // RDKit❗❌:     std::vector<SubstanceGroup> newsgs;
    // RDKit❗❌:     newsgs.reserve(sgs.size() - nRemoved);
    // RDKit❗❌:     unsigned int i = 0;
    // RDKit❗❌:     for (auto &&sg : sgs) {
    // RDKit❗❌:       if (!toRemove[i]) {
    // RDKit❗❌:         // we might be keeping it. Check the parent
    // RDKit❗❌:         if (!parentsPresent || !sg.hasProp("PARENT")) {
    // RDKit❗❌:           ADJUST_METHOD(sg, idx);
    // RDKit❗❌:           newsgs.push_back(std::move(sg));
    // RDKit❗❌:         } else if (parentsPresent) {
    // RDKit❗❌:           unsigned int parent;
    // RDKit❗❌:           // has our parent been removed?
    // RDKit❗❌:           if (sg.getPropIfPresent("PARENT", parent)) {
    // RDKit❗❌:             auto piter = indexLookup.find(parent);
    // RDKit❗❌:             bool keepIt = false;
    // RDKit❗❌:             if (piter == indexLookup.end()) {
    // RDKit❗❌:               // our parent isn't around, so it isn't being removed
    // RDKit❗❌:               // note: this is an odd case and probably shouldn't happen, but
    // RDKit❗❌:               // this isn't the place to enforce that
    // RDKit❗❌:               keepIt = true;
    // RDKit❗❌:             } else if (!toRemove[piter->second]) {
    // RDKit❗❌:               // our parent isn't being removed, recursively check up through
    // RDKit❗❌:               // parents to see if we find any that are being removed:
    // RDKit❗❌:               if (!removedParentInHierarchy(piter->second, sgs, toRemove,
    // RDKit❗❌:                                             indexLookup)) {
    // RDKit❗❌:                 keepIt = true;
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:             if (keepIt) {
    // RDKit❗❌:               ADJUST_METHOD(sg, idx);
    // RDKit❗❌:               newsgs.push_back(std::move(sg));
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       ++i;
    // RDKit❗❌:     }
    // RDKit❗❌:     sgs = std::move(newsgs);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Native move-from state on later conversion failure is not representable:
    // prior field adjustments remain, but ownership stays with this vector.
    // Dense detached IDs/typed parents must be remapped after source selection;
    // raw index/PARENT and crossing-property values remain unchanged.
    if groups.is_empty() {
        return Ok(());
    }
    let mut removed = Vec::with_capacity(groups.len());
    let mut parents_present = false;
    let mut n_removed = 0;
    for group in groups.iter() {
        if !parents_present && group.props.get(b"PARENT".as_slice()).is_some() {
            parents_present = true;
        }
        let drop = includes(group);
        n_removed += usize::from(drop);
        removed.push(drop);
    }
    let mut index_lookup = BTreeMap::new();
    if parents_present && n_removed != 0 {
        for (i, group) in groups.iter().enumerate() {
            if let Some(value) = group.props.get(b"index".as_slice()) {
                index_lookup.insert(uint_reader(value)?, i);
            }
        }
    }
    let mut keep = vec![false; groups.len()];
    for i in 0..groups.len() {
        if removed[i] {
            continue;
        }
        let retain = if !parents_present || !groups[i].props.get(b"PARENT".as_slice()).is_some() {
            true
        } else {
            let mut parent = uint_reader(groups[i].props.get(b"PARENT".as_slice()).unwrap())?;
            let mut depth = 0;
            loop {
                let Some(&position) = index_lookup.get(&parent) else {
                    break true;
                };
                if removed[position] {
                    break false;
                }
                let Some(value) = groups[position].props.get(b"PARENT".as_slice()) else {
                    break true;
                };
                if depth == groups.len() {
                    return Err(crate::TopologyEditError::SourceGroupParentCycle {
                        group: groups[i].id,
                    }
                    .into());
                }
                depth += 1;
                parent = uint_reader(value)?;
            }
        };
        if retain {
            adjust(&mut groups[i]);
            keep[i] = true;
        }
    }
    let mut ids = vec![None; groups.len()];
    let mut next = 0;
    for (i, retain) in keep.iter().copied().enumerate() {
        if retain {
            ids[i] = Some(SubstanceGroupId::new(next));
            next += 1;
        }
    }
    let mut i = 0;
    groups.retain_mut(|group| {
        let old_position = i;
        i += 1;
        if let Some(id) = ids[old_position] {
            group.id = id;
            group.parent = group
                .parent
                .and_then(|p| ids.get(p.index()).copied().flatten());
            true
        } else {
            false
        }
    });
    Ok(())
}

pub(crate) fn remove_source_groups_referencing_bond<E: From<crate::TopologyEditError>>(
    groups: &mut Vec<SubstanceGroup>,
    bond: BondId,
    uint_reader: &mut (impl FnMut(&PropertyValue) -> Result<u32, E> + ?Sized),
) -> Result<(), E> {
    // RDKit❗❌: removeSubstanceGroupsReferencing<includesBond, removedBond>(mol, idx);
    remove_source_groups_referencing(
        groups,
        |g| g.includes_bond(bond),
        |g| g.adjust_source_removed_bond(bond),
        uint_reader,
    )
}
pub(crate) fn remove_source_groups_referencing_atom<E: From<crate::TopologyEditError>>(
    groups: &mut Vec<SubstanceGroup>,
    atom: AtomId,
    uint_reader: &mut (impl FnMut(&PropertyValue) -> Result<u32, E> + ?Sized),
) -> Result<(), E> {
    // RDKit❗❌: removeSubstanceGroupsReferencing<includesAtom, removedAtom>(mol, idx);
    remove_source_groups_referencing(
        groups,
        |g| g.includes_atom(atom),
        |g| g.adjust_source_removed_atom(atom),
        uint_reader,
    )
}
