// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use crate::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, BondStereo, Conformer3D, Molecule,
    MoleculeBuildError, MoleculeProperties, StereoGroup, SubstanceGroup,
    molecule::{CoordinateBlock, TopologyBlock},
};

/// One-shot molecule construction API.
///
/// Builders are mutable by design. `Molecule` is not. If future code needs to
/// assemble parser output, it should use this type or a narrow parser-private
/// wrapper around it.
#[derive(Debug, Default)]
pub struct MoleculeBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    substance_groups: Vec<SubstanceGroup>,
    stereo_groups: Vec<StereoGroup>,
    coords_2d: Option<Vec<[f64; 2]>>,
    conformers_3d: Vec<Conformer3D>,
    properties: MoleculeProperties,
}

impl MoleculeBuilder {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_atom(&mut self, spec: AtomSpec) -> AtomId {
        let id = AtomId::new(self.atoms.len());
        self.atoms.push(Atom::from_spec(id, spec));
        if let Some(coords) = &mut self.coords_2d {
            coords.push([0.0, 0.0]);
        }
        id
    }

    pub fn add_bond(&mut self, spec: BondSpec) -> Result<BondId, MoleculeBuildError> {
        validate_bond_spec(self.atoms.len(), &spec)?;
        let id = BondId::new(self.bonds.len());
        self.bonds.push(Bond::from_spec(id, spec));
        Ok(id)
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
        if let Some(coords) = &mut self.coords_2d {
            let old_coords = std::mem::take(coords);
            *coords = old_coords
                .into_iter()
                .enumerate()
                .filter_map(|(idx, coord)| (!remove[idx]).then_some(coord))
                .collect();
        }
        for conformer in &mut self.conformers_3d {
            let old_coords = conformer.coords().to_vec();
            let coords = old_coords
                .into_iter()
                .enumerate()
                .filter_map(|(idx, coord)| (!remove[idx]).then_some(coord))
                .collect();
            *conformer = Conformer3D::new(conformer.id(), coords, conformer.is_3d());
        }
        atom_mapping
    }

    pub fn set_2d_coordinates(&mut self, coords: Vec<[f64; 2]>) -> Result<(), MoleculeBuildError> {
        if coords.len() != self.atoms.len() {
            return Err(MoleculeBuildError::CoordinateRowCount {
                rows: coords.len(),
                atom_count: self.atoms.len(),
            });
        }
        self.coords_2d = Some(coords);
        Ok(())
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

    pub(crate) fn add_conformer(
        &mut self,
        conformer: Conformer3D,
    ) -> Result<(), MoleculeBuildError> {
        if conformer.coords().len() != self.atoms.len() {
            return Err(MoleculeBuildError::ConformerRowCount {
                rows: conformer.coords().len(),
                atom_count: self.atoms.len(),
            });
        }
        self.conformers_3d.push(conformer);
        Ok(())
    }

    pub(crate) fn conformers_3d_len(&self) -> usize {
        self.conformers_3d.len()
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
        } else if self.coords_2d.is_some() || !self.conformers_3d.is_empty() {
            Some(crate::CoordinateDimension::TwoD)
        } else {
            None
        };
        Molecule::from_blocks(
            TopologyBlock {
                atoms: self.atoms,
                bonds: self.bonds,
                substance_groups: self.substance_groups,
                stereo_groups: self.stereo_groups,
            },
            CoordinateBlock {
                coords_2d: self.coords_2d,
                conformers_3d: self.conformers_3d,
                source_coordinate_dim,
            },
            self.properties,
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
            molecule.coords_2d().unwrap(),
            &[[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
        );
        assert_eq!(
            molecule.conformers_3d()[0].coords(),
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
