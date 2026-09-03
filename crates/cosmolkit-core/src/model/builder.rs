// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use crate::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, BondStereo, Conformer2D, Conformer3D, Molecule,
    MoleculeBuildError, MoleculeProperties, StereoGroup, SubstanceGroup,
    molecule::{
        AtomMapping, BondMapping, CoordinateBlock, MoleculeCapabilities, TopologyBlock,
        TopologyMapping, TopologyTrust,
    },
};

/// One-shot molecule construction API.
///
/// Builders are mutable by design. `Molecule` is not. If future code needs to
/// assemble parser output, it should use this type or a narrow parser-private
/// wrapper around it.
#[derive(Debug, Default, Clone)]
pub struct MoleculeBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    /// Build-time adjacency list, auto-maintained by `add_bond` and
    /// `add_atom`. Eliminates the need for callers to manually track
    /// adjacency during incremental construction (cf. SMILES parser's
    /// removed `atom_bonds`).
    ///
    /// `adjacency[i]` contains the `BondId`s of all bonds incident to
    /// atom `i`.  Invariant: `adjacency.len() == atoms.len()`.
    adjacency: Vec<Vec<BondId>>,
    substance_groups: Vec<SubstanceGroup>,
    stereo_groups: Vec<StereoGroup>,
    conformers_2d: Vec<Conformer2D>,
    conformers_3d: Vec<Conformer3D>,
    properties: MoleculeProperties,
    capabilities: MoleculeCapabilities,
}

impl MoleculeBuilder {
    #[must_use]
    pub fn new() -> Self {
        Self {
            capabilities: MoleculeCapabilities::new(TopologyTrust::TrustedGraph),
            ..Self::default()
        }
    }

    pub fn add_atom(&mut self, spec: AtomSpec) -> AtomId {
        let id = AtomId::new(self.atoms.len());
        self.atoms.push(Atom::from_spec(id, spec));
        for conformer in &mut self.conformers_2d {
            conformer.push_coord([0.0, 0.0]);
        }
        self.adjacency.push(Vec::new());
        id
    }

    pub fn add_bond(&mut self, spec: BondSpec) -> Result<BondId, MoleculeBuildError> {
        validate_bond_spec(self.atoms.len(), &spec)?;
        let id = BondId::new(self.bonds.len());
        let begin = spec.begin();
        let end = spec.end();
        self.bonds.push(Bond::from_spec(id, spec));
        self.adjacency[begin.index()].push(id);
        self.adjacency[end.index()].push(id);
        Ok(id)
    }

    pub fn set_atom_formal_charge(
        &mut self,
        atom: AtomId,
        formal_charge: i8,
    ) -> Result<(), MoleculeBuildError> {
        let Some(atom) = self.atoms.get_mut(atom.index()) else {
            return Err(MoleculeBuildError::InvalidMoleculeState {
                message: "atom index out of range while setting formal charge".to_string(),
            });
        };
        atom.set_formal_charge(formal_charge);
        Ok(())
    }

    /// Returns the number of bonds incident to `atom`.
    pub fn degree(&self, atom: AtomId) -> usize {
        self.adjacency.get(atom.index()).map_or(0, Vec::len)
    }

    /// Returns the bonds incident to `atom`.
    pub fn neighbor_bonds(&self, atom: AtomId) -> &[BondId] {
        self.adjacency
            .get(atom.index())
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub(crate) fn bond(&self, bond: BondId) -> Option<&Bond> {
        self.bonds.get(bond.index())
    }

    pub(crate) fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    pub(crate) fn atom_mut(&mut self, atom: AtomId) -> Option<&mut Atom> {
        self.atoms.get_mut(atom.index())
    }

    pub(crate) fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    pub(crate) fn bond_mut(&mut self, bond: BondId) -> Option<&mut Bond> {
        self.bonds.get_mut(bond.index())
    }

    pub(crate) fn atoms_mut(&mut self) -> &mut [Atom] {
        &mut self.atoms
    }

    pub(crate) fn bonds_mut(&mut self) -> &mut [Bond] {
        &mut self.bonds
    }

    pub(crate) fn bond_between_atoms(&self, begin: AtomId, end: AtomId) -> Option<BondId> {
        let begin_bonds = self.adjacency.get(begin.index())?;
        begin_bonds.iter().copied().find(|bond_id| {
            let bond = &self.bonds[bond_id.index()];
            (bond.begin() == begin && bond.end() == end)
                || (bond.begin() == end && bond.end() == begin)
        })
    }

    pub(crate) fn set_bond_order(
        &mut self,
        bond: BondId,
        order: crate::BondOrder,
    ) -> Result<(), MoleculeBuildError> {
        let Some(bond) = self.bonds.get_mut(bond.index()) else {
            return Err(MoleculeBuildError::InvalidMoleculeState {
                message: "bond index out of range while setting bond order".to_string(),
            });
        };
        bond.set_order(order);
        Ok(())
    }

    pub(crate) fn remove_bond_between_atoms(
        &mut self,
        begin: AtomId,
        end: AtomId,
    ) -> Option<BondId> {
        let bond_id = self.bond_between_atoms(begin, end)?;
        self.remove_bond_for_construction(bond_id);
        Some(bond_id)
    }

    fn remove_bond_for_construction(&mut self, bond_to_remove: BondId) {
        self.bonds.remove(bond_to_remove.index());
        for (index, bond) in self.bonds.iter_mut().enumerate() {
            bond.set_id_for_construction(BondId::new(index));
        }
        self.adjacency = vec![Vec::new(); self.atoms.len()];
        for bond in &self.bonds {
            self.adjacency[bond.begin().index()].push(bond.id());
            self.adjacency[bond.end().index()].push(bond.id());
        }
    }

    pub(crate) fn remove_atoms_for_construction(
        &mut self,
        atoms_to_remove: &[AtomId],
    ) -> Vec<Option<AtomId>> {
        let mut remove = vec![false; self.atoms.len()];
        for atom in atoms_to_remove {
            if let Some(slot) = remove.get_mut(atom.index()) {
                *slot = true;
            }
        }

        let mut atom_mapping = vec![None; self.atoms.len()];
        let mut atoms = Vec::with_capacity(self.atoms.len().saturating_sub(atoms_to_remove.len()));
        for atom in self.atoms.drain(..) {
            if remove[atom.id().index()] {
                continue;
            }
            let new_id = AtomId::new(atoms.len());
            atom_mapping[atom.id().index()] = Some(new_id);
            atoms.push(atom.with_id(new_id));
        }

        let mut bonds = Vec::new();
        for bond in self.bonds.drain(..) {
            let Some(begin) = atom_mapping[bond.begin().index()] else {
                continue;
            };
            let Some(end) = atom_mapping[bond.end().index()] else {
                continue;
            };
            let stereo_atoms = bond.stereo_atoms().and_then(|[begin_ref, end_ref]| {
                Some([
                    atom_mapping[begin_ref.index()]?,
                    atom_mapping[end_ref.index()]?,
                ])
            });
            let new_id = BondId::new(bonds.len());
            bonds.push(bond.remapped(new_id, begin, end, stereo_atoms));
        }

        self.atoms = atoms;
        self.bonds = bonds;
        self.adjacency = vec![Vec::new(); self.atoms.len()];
        for bond in &self.bonds {
            self.adjacency[bond.begin().index()].push(bond.id());
            self.adjacency[bond.end().index()].push(bond.id());
        }
        for conformer in &mut self.conformers_2d {
            let old_coords = conformer.coordinates().to_vec();
            let coords = old_coords
                .into_iter()
                .enumerate()
                .filter_map(|(idx, coord)| (!remove[idx]).then_some(coord))
                .collect();
            *conformer = Conformer2D::new(conformer.id(), coords);
        }
        for conformer in &mut self.conformers_3d {
            let old_coords = conformer.coordinates().to_vec();
            let coords = old_coords
                .into_iter()
                .enumerate()
                .filter_map(|(idx, coord)| (!remove[idx]).then_some(coord))
                .collect();
            *conformer = Conformer3D::new(conformer.id(), coords, conformer.is_3d());
        }
        atom_mapping
    }

    pub(crate) fn renumber_atoms_for_construction(
        &mut self,
        old_atom_order: &[AtomId],
    ) -> Result<TopologyMapping, MoleculeBuildError> {
        if old_atom_order.len() != self.atoms.len() {
            return Err(MoleculeBuildError::InvalidMoleculeState {
                message: "atom renumbering order length does not match atom count".to_string(),
            });
        }
        let mut seen = vec![false; self.atoms.len()];
        for atom in old_atom_order {
            if atom.index() >= self.atoms.len() || seen[atom.index()] {
                return Err(MoleculeBuildError::InvalidMoleculeState {
                    message: "atom renumbering order is not a permutation".to_string(),
                });
            }
            seen[atom.index()] = true;
        }

        let mut atom_old_to_new = vec![None; self.atoms.len()];
        let mut atom_new_to_old = Vec::with_capacity(self.atoms.len());
        let mut atoms = Vec::with_capacity(self.atoms.len());
        for (new_idx, old_atom) in old_atom_order.iter().copied().enumerate() {
            let new_id = AtomId::new(new_idx);
            atom_old_to_new[old_atom.index()] = Some(new_id);
            atom_new_to_old.push(Some(old_atom));
            atoms.push(self.atoms[old_atom.index()].clone().with_id(new_id));
        }

        let mut bond_old_to_new = vec![None; self.bonds.len()];
        let mut bond_new_to_old = Vec::with_capacity(self.bonds.len());
        let mut bonds = Vec::with_capacity(self.bonds.len());
        for bond in &self.bonds {
            let begin = atom_old_to_new[bond.begin().index()]
                .expect("atom permutation maps every old atom");
            let end =
                atom_old_to_new[bond.end().index()].expect("atom permutation maps every old atom");
            let stereo_atoms = bond.stereo_atoms().and_then(|[left, right]| {
                Some([
                    atom_old_to_new[left.index()]?,
                    atom_old_to_new[right.index()]?,
                ])
            });
            let new_id = BondId::new(bonds.len());
            bond_old_to_new[bond.id().index()] = Some(new_id);
            bond_new_to_old.push(Some(bond.id()));
            bonds.push(bond.clone().remapped(new_id, begin, end, stereo_atoms));
        }

        for conformer in &mut self.conformers_2d {
            let coords = old_atom_order
                .iter()
                .map(|atom| conformer.coordinates()[atom.index()])
                .collect();
            *conformer = Conformer2D::new(conformer.id(), coords);
        }
        for conformer in &mut self.conformers_3d {
            let coords = old_atom_order
                .iter()
                .map(|atom| conformer.coordinates()[atom.index()])
                .collect();
            *conformer = Conformer3D::new(conformer.id(), coords, conformer.is_3d());
        }

        self.substance_groups = self
            .substance_groups
            .iter()
            .enumerate()
            .filter_map(|(old_index, sgroup)| {
                let sgroup_map = (0..self.substance_groups.len())
                    .map(crate::SubstanceGroupId::new)
                    .map(Some)
                    .collect::<Vec<_>>();
                sgroup.remapped(
                    crate::SubstanceGroupId::new(old_index),
                    &atom_old_to_new,
                    &bond_old_to_new,
                    &sgroup_map,
                )
            })
            .collect();
        self.stereo_groups = self
            .stereo_groups
            .iter()
            .filter_map(|group| group.remapped(&atom_old_to_new, &bond_old_to_new))
            .collect();

        self.properties
            .remap_topology(&atom_new_to_old, &bond_new_to_old);

        self.atoms = atoms;
        self.bonds = bonds;
        self.adjacency = vec![Vec::new(); self.atoms.len()];
        for bond in &self.bonds {
            self.adjacency[bond.begin().index()].push(bond.id());
            self.adjacency[bond.end().index()].push(bond.id());
        }

        Ok(TopologyMapping {
            atoms: AtomMapping {
                old_to_new: atom_old_to_new,
                new_to_old: atom_new_to_old,
            },
            bonds: BondMapping {
                old_to_new: bond_old_to_new,
                new_to_old: bond_new_to_old,
            },
        })
    }

    pub fn set_2d_coordinates(&mut self, coords: Vec<[f64; 2]>) -> Result<(), MoleculeBuildError> {
        if coords.len() != self.atoms.len() {
            return Err(MoleculeBuildError::CoordinateRowCount {
                rows: coords.len(),
                atom_count: self.atoms.len(),
            });
        }
        self.conformers_2d.clear();
        self.conformers_2d.push(Conformer2D::new(0, coords));
        Ok(())
    }

    pub fn add_2d_conformer(&mut self, coords: Vec<[f64; 2]>) -> Result<usize, MoleculeBuildError> {
        if coords.len() != self.atoms.len() {
            return Err(MoleculeBuildError::CoordinateRowCount {
                rows: coords.len(),
                atom_count: self.atoms.len(),
            });
        }
        let id = self.conformers_2d.len();
        self.conformers_2d.push(Conformer2D::new(id, coords));
        Ok(id)
    }

    pub fn add_3d_conformer(&mut self, coords: Vec<[f64; 3]>) -> Result<usize, MoleculeBuildError> {
        if coords.len() != self.atoms.len() {
            return Err(MoleculeBuildError::ConformerRowCount {
                rows: coords.len(),
                atom_count: self.atoms.len(),
            });
        }
        let id = self.conformers_3d.len();
        self.conformers_3d.push(Conformer3D::new(id, coords, true));
        Ok(id)
    }

    /// Add a conformer while preserving its explicit id and dimensionality.
    pub fn add_conformer(&mut self, conformer: Conformer3D) -> Result<(), MoleculeBuildError> {
        if conformer.coordinates().len() != self.atoms.len() {
            return Err(MoleculeBuildError::ConformerRowCount {
                rows: conformer.coordinates().len(),
                atom_count: self.atoms.len(),
            });
        }
        self.conformers_3d.push(conformer);
        Ok(())
    }

    pub(crate) fn conformers_3d_len(&self) -> usize {
        self.conformers_3d.len()
    }

    pub(crate) fn conformers_3d(&self) -> &[Conformer3D] {
        &self.conformers_3d
    }

    pub fn add_substance_group(
        &mut self,
        substance_group: SubstanceGroup,
    ) -> Result<(), MoleculeBuildError> {
        validate_substance_group(self.atoms.len(), self.bonds.len(), &substance_group)?;
        self.substance_groups.push(substance_group);
        Ok(())
    }

    pub(crate) fn substance_groups_len(&self) -> usize {
        self.substance_groups.len()
    }

    pub(crate) fn substance_groups(&self) -> &[SubstanceGroup] {
        &self.substance_groups
    }

    pub(crate) fn substance_group_mut(&mut self, index: usize) -> Option<&mut SubstanceGroup> {
        self.substance_groups.get_mut(index)
    }

    pub fn add_stereo_group(
        &mut self,
        stereo_group: StereoGroup,
    ) -> Result<(), MoleculeBuildError> {
        validate_stereo_group(self.atoms.len(), self.bonds.len(), &stereo_group)?;
        self.stereo_groups.push(stereo_group);
        Ok(())
    }

    pub(crate) fn stereo_groups_mut(&mut self) -> &mut [StereoGroup] {
        &mut self.stereo_groups
    }

    pub(crate) fn stereo_groups_len(&self) -> usize {
        self.stereo_groups.len()
    }

    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.properties = self.properties.with_name(name);
        self
    }

    #[must_use]
    pub fn with_property(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.properties = self.properties.with_prop(key, value);
        self
    }

    #[must_use]
    pub fn with_sdf_data_field(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.properties = self.properties.with_sdf_data_field(key, value);
        self
    }

    #[must_use]
    #[allow(dead_code)]
    pub(crate) fn with_properties(mut self, properties: MoleculeProperties) -> Self {
        self.properties = properties;
        self
    }

    #[must_use]
    pub fn with_topology_trust(mut self, topology_trust: TopologyTrust) -> Self {
        self.capabilities = self.capabilities.with_topology_trust(topology_trust);
        self
    }

    pub(crate) fn set_topology_trust(&mut self, topology_trust: TopologyTrust) {
        self.capabilities = self.capabilities.with_topology_trust(topology_trust);
    }

    pub fn build(self) -> Result<Molecule, MoleculeBuildError> {
        for bond in &self.bonds {
            validate_bond_spec(self.atoms.len(), &bond_spec_from_bond(bond))?;
        }
        for substance_group in &self.substance_groups {
            validate_substance_group(self.atoms.len(), self.bonds.len(), substance_group)?;
        }
        for stereo_group in &self.stereo_groups {
            validate_stereo_group(self.atoms.len(), self.bonds.len(), stereo_group)?;
        }
        let source_coordinate_dim = if self.conformers_3d.iter().any(Conformer3D::is_3d) {
            Some(crate::CoordinateDimension::ThreeD)
        } else if !self.conformers_2d.is_empty() || !self.conformers_3d.is_empty() {
            Some(crate::CoordinateDimension::TwoD)
        } else {
            None
        };
        Molecule::from_blocks_with_capabilities(
            TopologyBlock {
                atoms: self.atoms,
                bonds: self.bonds,
                adjacency: crate::AdjacencyList::default(),
                substance_groups: self.substance_groups,
                stereo_groups: self.stereo_groups,
            },
            CoordinateBlock {
                conformers_2d: self.conformers_2d,
                conformers_3d: self.conformers_3d,
                source_coordinate_dim,
            },
            self.properties,
            self.capabilities,
        )
        .map_err(|err| match err {
            crate::InvariantError::InvalidBondEndpoint {
                begin,
                end,
                atom_count,
                ..
            } => MoleculeBuildError::BondEndpointOutOfRange {
                begin,
                end,
                atom_count,
            },
            crate::InvariantError::SelfLoopBond { atom, .. } => {
                MoleculeBuildError::SelfLoopBond { atom }
            }
            crate::InvariantError::CoordinateRowCount { rows, atom_count } => {
                MoleculeBuildError::CoordinateRowCount { rows, atom_count }
            }
            crate::InvariantError::ConformerRowCount {
                rows, atom_count, ..
            } => MoleculeBuildError::ConformerRowCount { rows, atom_count },
            crate::InvariantError::InvalidSubstanceGroupAtom {
                atom, atom_count, ..
            } => MoleculeBuildError::SubstanceGroupAtomOutOfRange { atom, atom_count },
            crate::InvariantError::InvalidSubstanceGroupBond {
                bond, bond_count, ..
            } => MoleculeBuildError::SubstanceGroupBondOutOfRange { bond, bond_count },
            crate::InvariantError::InvalidSubstanceGroupParent { parent, .. } => {
                MoleculeBuildError::SubstanceGroupParentOutOfRange { parent }
            }
            crate::InvariantError::SubstanceGroupIndexMismatch { actual, .. } => {
                MoleculeBuildError::SubstanceGroupParentOutOfRange { parent: actual }
            }
            other => MoleculeBuildError::InvalidMoleculeState {
                message: other.to_string(),
            },
        })
    }
}

impl Molecule {
    #[must_use]
    pub fn to_builder(&self) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new().with_topology_trust(self.topology_trust());

        for atom in self.atoms() {
            let mut spec = AtomSpec::new(atom.element())
                .with_formal_charge(atom.formal_charge())
                .with_explicit_hydrogens(atom.explicit_hydrogens())
                .with_chiral_tag(atom.chiral_tag())
                .with_unknown_stereo(atom.unknown_stereo())
                .with_implicit_hydrogen(atom.implicit_hydrogen())
                .with_tracked_isotopic_hydrogens(atom.tracked_isotopic_hydrogens().to_vec())
                .with_aromatic(atom.is_aromatic())
                .with_no_implicit(atom.no_implicit())
                .with_radical_electrons(atom.radical_electrons())
                .with_hybridization(atom.hybridization());
            if let Some(chiral_permutation) = atom.chiral_permutation() {
                spec = spec.with_chiral_permutation(chiral_permutation);
            }
            if let Some(mol_parity) = atom.mol_parity() {
                spec = spec.with_mol_parity(mol_parity);
            }
            if let Some(mol_inversion_flag) = atom.mol_inversion_flag() {
                spec = spec.with_mol_inversion_flag(mol_inversion_flag);
            }
            if let Some(isotope) = atom.isotope() {
                spec = spec.with_isotope(isotope);
            }
            if let Some(atom_map) = atom.atom_map() {
                spec = spec.with_atom_map(atom_map);
            }
            if let Some(query) = atom.query().cloned() {
                spec = spec.with_query(query);
            }
            if let Some(info) = atom.pdb_residue_info().cloned() {
                spec = spec.with_pdb_residue_info(info);
            }
            for (key, value) in atom.props() {
                spec = if atom.is_prop_computed(key) {
                    spec.with_computed_prop(key.clone(), value.clone())
                } else {
                    spec.with_prop(key.clone(), value.clone())
                };
            }
            builder.add_atom(spec);
        }

        for bond in self.bonds() {
            let mut spec = BondSpec::new(bond.begin(), bond.end(), bond.order())
                .with_aromatic(bond.is_aromatic())
                .with_conjugated(bond.is_conjugated())
                .with_direction(bond.direction())
                .with_stereo(bond.stereo())
                .with_unknown_stereo(bond.unknown_stereo());
            if let Some([begin_ref, end_ref]) = bond.stereo_atoms() {
                spec = spec.with_stereo_atoms(begin_ref, end_ref);
            }
            if let Some(query) = bond.query().cloned() {
                spec = spec.with_query(query);
            }
            for (key, value) in bond.props() {
                spec = if bond.is_prop_computed(key) {
                    spec.with_computed_prop(key.clone(), value.clone())
                } else {
                    spec.with_prop(key.clone(), value.clone())
                };
            }
            builder
                .add_bond(spec)
                .expect("existing molecule bond topology must be buildable");
        }

        for conformer in self.conformers_2d() {
            builder
                .add_2d_conformer(conformer.coordinates().to_vec())
                .expect("existing molecule 2D conformer row count must be valid");
        }
        for conformer in self.conformers_3d() {
            builder
                .add_conformer(conformer.clone())
                .expect("existing molecule 3D conformer row count must be valid");
        }
        for substance_group in self.substance_groups() {
            builder
                .add_substance_group(substance_group.clone())
                .expect("existing molecule substance group must be buildable");
        }
        for stereo_group in self.stereo_groups() {
            builder
                .add_stereo_group(stereo_group.clone())
                .expect("existing molecule stereo group must be buildable");
        }

        builder.with_properties(self.properties().clone())
    }
}

pub(crate) fn validate_substance_group(
    atom_count: usize,
    bond_count: usize,
    substance_group: &SubstanceGroup,
) -> Result<(), MoleculeBuildError> {
    for atom in substance_group.atoms() {
        if atom.index() >= atom_count {
            return Err(MoleculeBuildError::SubstanceGroupAtomOutOfRange {
                atom: *atom,
                atom_count,
            });
        }
    }
    for atom in substance_group.parent_atoms() {
        if atom.index() >= atom_count {
            return Err(MoleculeBuildError::SubstanceGroupAtomOutOfRange {
                atom: *atom,
                atom_count,
            });
        }
    }
    for bond in substance_group.bonds() {
        if bond.index() >= bond_count {
            return Err(MoleculeBuildError::SubstanceGroupBondOutOfRange {
                bond: *bond,
                bond_count,
            });
        }
    }
    for attach_point in substance_group.attach_points() {
        if attach_point.atom.index() >= atom_count {
            return Err(MoleculeBuildError::SubstanceGroupAtomOutOfRange {
                atom: attach_point.atom,
                atom_count,
            });
        }
        if let Some(leaving_atom) = attach_point.leaving_atom
            && leaving_atom.index() >= atom_count
        {
            return Err(MoleculeBuildError::SubstanceGroupAtomOutOfRange {
                atom: leaving_atom,
                atom_count,
            });
        }
    }
    for cstate in substance_group.cstates() {
        if cstate.bond.index() >= bond_count {
            return Err(MoleculeBuildError::SubstanceGroupBondOutOfRange {
                bond: cstate.bond,
                bond_count,
            });
        }
    }
    Ok(())
}

pub(crate) fn validate_stereo_group(
    atom_count: usize,
    bond_count: usize,
    stereo_group: &StereoGroup,
) -> Result<(), MoleculeBuildError> {
    for atom in stereo_group.atoms() {
        if atom.index() >= atom_count {
            return Err(MoleculeBuildError::SubstanceGroupAtomOutOfRange {
                atom: *atom,
                atom_count,
            });
        }
    }
    for bond in stereo_group.bonds() {
        if bond.index() >= bond_count {
            return Err(MoleculeBuildError::SubstanceGroupBondOutOfRange {
                bond: *bond,
                bond_count,
            });
        }
    }
    Ok(())
}

fn bond_spec_from_bond(bond: &Bond) -> BondSpec {
    let mut spec = BondSpec::new(bond.begin(), bond.end(), bond.order())
        .with_aromatic(bond.is_aromatic())
        .with_conjugated(bond.is_conjugated())
        .with_direction(bond.direction())
        .with_stereo(bond.stereo());
    if let Some([begin_ref, end_ref]) = bond.stereo_atoms() {
        spec = spec.with_stereo_atoms(begin_ref, end_ref);
    }
    spec
}

pub(crate) fn validate_bond_spec(
    atom_count: usize,
    spec: &BondSpec,
) -> Result<(), MoleculeBuildError> {
    if spec.begin() == spec.end() {
        return Err(MoleculeBuildError::SelfLoopBond { atom: spec.begin() });
    }
    if spec.begin().index() >= atom_count || spec.end().index() >= atom_count {
        return Err(MoleculeBuildError::BondEndpointOutOfRange {
            begin: spec.begin(),
            end: spec.end(),
            atom_count,
        });
    }
    if let Some([begin_ref, end_ref]) = spec.stereo_atoms()
        && (begin_ref.index() >= atom_count || end_ref.index() >= atom_count)
    {
        return Err(MoleculeBuildError::BondEndpointOutOfRange {
            begin: begin_ref,
            end: end_ref,
            atom_count,
        });
    }
    // BEGIN RDKIT CPP FUNCTION Bond::setStereo
    // RDKit✔️✔️: void setStereo(BondStereo what) {
    // RDKit✔️✔️:   PRECONDITION(((what != STEREOCIS && what != STEREOTRANS) ||
    // RDKit✔️✔️:                 getStereoAtoms().size() == 2),
    // RDKit✔️✔️:                "Stereo atoms should be specified before specifying CIS/TRANS "
    // RDKit✔️✔️:                "bond stereochemistry")
    // RDKit✔️✔️:   d_stereo = what;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Bond::setStereo
    if matches!(spec.stereo(), BondStereo::Cis | BondStereo::Trans) && spec.stereo_atoms().is_none()
    {
        return Err(MoleculeBuildError::BondStereoAtomsRequired {
            stereo: spec.stereo(),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::{
        AtomQueryPredicate, BondDirection, BondOrder, BondQueryPredicate, BondSpec, BondStereo,
        ChiralTag, Element, Hybridization, MoleculeBuilder, QueryNode, SubstanceGroup,
        SubstanceGroupId, SubstanceGroupKind,
    };

    #[test]
    fn builder_preserves_atom_and_bond_state() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(
            crate::AtomSpec::new(Element::C)
                .with_formal_charge(-1)
                .with_explicit_hydrogens(1)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_aromatic(true)
                .with_isotope(13)
                .with_atom_map(7)
                .with_no_implicit(true)
                .with_radical_electrons(1)
                .with_hybridization(Hybridization::Sp2)
                .with_query(QueryNode::predicate(AtomQueryPredicate::Any)),
        );
        let oxygen = builder.add_atom(
            crate::AtomSpec::new(Element::O)
                .with_formal_charge(1)
                .with_explicit_hydrogens(2),
        );
        let dummy = builder.add_atom(crate::AtomSpec::new(Element::DUMMY));
        builder
            .add_bond(
                BondSpec::new(carbon, oxygen, BondOrder::Double)
                    .with_aromatic(true)
                    .with_conjugated(true)
                    .with_direction(BondDirection::EndUpRight)
                    .with_stereo(BondStereo::Cis)
                    .with_stereo_atoms(carbon, oxygen)
                    .with_query(QueryNode::predicate(BondQueryPredicate::Any)),
            )
            .unwrap();
        builder
            .set_2d_coordinates(vec![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
            .unwrap();
        builder
            .add_3d_conformer(vec![[1.0, 2.0, 0.0], [3.0, 4.0, 0.5], [5.0, 6.0, 1.0]])
            .unwrap();
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                    .with_atoms(vec![carbon, oxygen])
                    .with_prop("FIELDNAME", "example")
                    .with_data_field("payload"),
            )
            .unwrap();
        builder = builder
            .with_property("_MolFileInfo", "info")
            .with_sdf_data_field("ID", "cmpd-1");

        let molecule = builder.build().unwrap();
        let carbon = molecule.atom(carbon).unwrap();
        assert_eq!(carbon.atomic_number(), 6);
        assert_eq!(carbon.formal_charge(), -1);
        assert_eq!(carbon.explicit_hydrogens(), 1);
        assert_eq!(carbon.chiral_tag(), ChiralTag::TetrahedralCw);
        assert!(carbon.is_aromatic());
        assert_eq!(carbon.isotope(), Some(13));
        assert_eq!(carbon.atom_map(), Some(7));
        assert!(carbon.no_implicit());
        assert_eq!(carbon.radical_electrons(), 1);
        assert_eq!(carbon.hybridization(), Hybridization::Sp2);
        assert_eq!(
            carbon.query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::Any))
        );

        let oxygen = molecule.atom(oxygen).unwrap();
        assert_eq!(oxygen.atomic_number(), 8);
        assert_eq!(oxygen.formal_charge(), 1);
        assert_eq!(oxygen.explicit_hydrogens(), 2);
        assert_eq!(molecule.atom(dummy).unwrap().atomic_number(), 0);

        let bond = &molecule.bonds()[0];
        assert_eq!(bond.begin(), carbon.id());
        assert_eq!(bond.end(), oxygen.id());
        assert_eq!(bond.order(), BondOrder::Double);
        assert!(bond.is_aromatic());
        assert!(bond.is_conjugated());
        assert_eq!(bond.direction(), BondDirection::EndUpRight);
        assert_eq!(bond.stereo(), BondStereo::Cis);
        assert_eq!(bond.stereo_atoms(), Some([carbon.id(), oxygen.id()]));
        assert_eq!(
            bond.query(),
            Some(&QueryNode::predicate(BondQueryPredicate::Any))
        );
        assert_eq!(
            molecule.coordinates_2d().unwrap(),
            &[[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
        );
        assert_eq!(
            molecule.conformers_3d()[0].coordinates(),
            &[[1.0, 2.0, 0.0], [3.0, 4.0, 0.5], [5.0, 6.0, 1.0]]
        );
        assert_eq!(
            molecule.source_coordinate_dim(),
            Some(crate::CoordinateDimension::ThreeD)
        );
        assert_eq!(molecule.substance_groups().len(), 1);
        assert_eq!(
            molecule.substance_groups()[0].atoms(),
            &[carbon.id(), oxygen.id()]
        );
        assert_eq!(molecule.prop("_MolFileInfo"), Some("info"));
        assert_eq!(
            molecule.properties().sdf_data_fields(),
            &[("ID".to_string(), "cmpd-1".to_string())]
        );
    }
}
