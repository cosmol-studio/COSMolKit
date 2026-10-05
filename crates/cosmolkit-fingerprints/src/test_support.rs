//! Shared detached source-prepared fixed-test fixture.
use crate::AtomPairPreparedInput;
use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer3D, Element,
};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};
pub(crate) struct TestMolecule {
    pub(crate) topology: TopologyBlock,
    pub(crate) properties: MoleculeProperties,
    pub(crate) coordinates: CoordinateBlock,
    pub(crate) valence: ValenceAssignment,
    pub(crate) rings: RingInfo,
}
impl TestMolecule {
    pub(crate) fn prepared(topology: TopologyBlock, coordinates: CoordinateBlock) -> Self {
        let properties = MoleculeProperties::default();
        let assignment = cosmolkit_core::sanitize_topology(
            &topology,
            &cosmolkit_core::SanitizeParams::default(),
        )
        .unwrap();
        Self {
            topology: assignment.topology,
            properties,
            coordinates,
            valence: assignment.final_valence.unwrap(),
            rings: assignment.final_rings.unwrap(),
        }
    }
    pub(crate) fn from_smiles(smiles: &str) -> Result<Self, cosmolkit_smiles::SmilesParseError> {
        let params = cosmolkit_smiles::SmilesParseParams::default();
        let record = cosmolkit_smiles::parse_smiles(smiles, &params)?;
        let assignment = cosmolkit_core::sanitize_topology(
            &record.topology,
            &cosmolkit_core::SanitizeParams::default(),
        )
        .unwrap();
        let mut valence = assignment.final_valence;
        let mut rings = assignment.final_rings;
        let finalized = cosmolkit_smiles::finalize_smiles_stereo(
            cosmolkit_smiles::SmilesRecord {
                topology: assignment.topology,
                ..record
            },
            &params,
            &mut valence,
            &mut rings,
        )
        .unwrap();
        Ok(Self {
            topology: finalized.topology,
            properties: finalized.properties,
            coordinates: finalized.coordinates,
            valence: valence.unwrap(),
            rings: rings.unwrap(),
        })
    }
    pub(crate) fn new() -> Self {
        Self::prepared(TopologyBlock::default(), CoordinateBlock::default())
    }
    pub(crate) fn num_atoms(&self) -> usize {
        self.topology.atoms.len()
    }
    pub(crate) fn input(&self) -> AtomPairPreparedInput<'_> {
        AtomPairPreparedInput {
            topology: &self.topology,
            properties: &self.properties,
            coordinates: &self.coordinates,
            valence: &self.valence,
            rings: &self.rings,
            use_legacy_stereo_perception: true,
        }
    }
}
#[derive(Default)]
pub(crate) struct TestBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    coordinates: CoordinateBlock,
}
impl TestBuilder {
    pub(crate) fn add_atom(&mut self, spec: AtomSpec) -> AtomId {
        let id = AtomId::new(self.atoms.len());
        self.atoms.push(Atom::from_spec(id, spec));
        id
    }
    pub(crate) fn add_bond(&mut self, spec: BondSpec) -> Result<(), ()> {
        self.bonds
            .push(Bond::from_spec(BondId::new(self.bonds.len()), spec));
        Ok(())
    }
    pub(crate) fn add_conformer(&mut self, conformer: Conformer3D) -> Result<(), ()> {
        self.coordinates.conformers_3d.push(conformer);
        Ok(())
    }
    pub(crate) fn build(self) -> Result<TestMolecule, cosmolkit_model::TopologyValidationError> {
        let topology =
            TopologyBlock::try_from_parts(self.atoms, self.bonds, Vec::new(), Vec::new())?;
        Ok(TestMolecule::prepared(topology, self.coordinates))
    }
}
