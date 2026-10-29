//! Native detached builder; publication and COW remain Rust-owned.
use crate::{
    Molecule,
    construction_values::{AtomSpec, AtomSpecInput, BondSpec, BondSpecInput},
    enums::BondOrderInput,
    fields::Field,
    group_values::{Point2, Point3, SubstanceGroup, SubstanceGroupInput},
    input,
    metadata_values::{MoleculeProperties, MoleculePropertiesInput},
    operation_errors::OperationError,
    stereo_groups::{StereoGroup, StereoGroupInput},
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
#[derive(serde::Serialize, serde::Deserialize, Tsify)]
#[serde(transparent)]
pub struct Points2(pub(crate) Vec<Point2>);
#[derive(serde::Serialize, serde::Deserialize, Tsify)]
#[serde(transparent)]
pub struct Points3(pub(crate) Vec<Point3>);
#[wasm_bindgen]
pub struct MoleculeBuilder(pub(crate) ck::MoleculeBuilder);
#[wasm_bindgen]
impl MoleculeBuilder {
    #[wasm_bindgen(js_name=fromParts)]
    pub fn from_parts(
        topology: crate::topology_values::TopologyBlockInput,
        coordinates: crate::detached_coordinates::CoordinateBlockInput,
        properties: MoleculePropertiesInput,
    ) -> Result<Self, JsValue> {
        crate::topology_values::TopologyBlock::with_ref(&topology, |t| {
            crate::detached_coordinates::CoordinateBlock::with_ref(&coordinates, |c| {
                MoleculeProperties::with_ref(&properties, |p| {
                    Ok(Self(ck::MoleculeBuilder::from_parts(
                        t.clone(),
                        c.clone(),
                        p.clone(),
                    )))
                })
            })
        })
    }
    #[wasm_bindgen(getter)]
    pub fn coordinates(&self) -> crate::detached_coordinates::CoordinateBlock {
        crate::detached_coordinates::CoordinateBlock(self.0.coordinates().clone())
    }
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self(ck::MoleculeBuilder::new())
    }
    pub fn build(&self) -> Result<Molecule, JsValue> {
        self.0
            .clone()
            .build()
            .map(Molecule)
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=addAtom)]
    pub fn add_atom(&mut self, atom: AtomSpecInput) -> Result<usize, JsValue> {
        AtomSpec::with_ref(&atom, |v| Ok(self.0.add_atom(v.clone()).index()))
    }
    #[wasm_bindgen(js_name=addBond)]
    pub fn add_bond(&mut self, bond: BondSpecInput) -> Result<usize, JsValue> {
        BondSpec::with_ref(&bond, |v| {
            self.0
                .add_bond(v.clone())
                .map(ck::BondId::index)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=setAtomFormalCharge)]
    pub fn set_atom_formal_charge(
        &mut self,
        atom: Ts<input::IndexInput>,
        charge: Ts<input::I8Input>,
    ) -> Result<(), JsValue> {
        self.0
            .set_atom_formal_charge(ck::AtomId::new(input::read(atom)?), input::read(charge)?)
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=setBondOrder)]
    pub fn set_bond_order(
        &mut self,
        bond: Ts<input::IndexInput>,
        order: BondOrderInput,
    ) -> Result<(), JsValue> {
        self.0
            .set_bond_order(
                ck::BondId::new(input::read(bond)?),
                ck::BondOrder::input(order)?,
            )
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=removeBondBetweenAtoms)]
    pub fn remove_bond_between_atoms(
        &mut self,
        begin: Ts<input::IndexInput>,
        end: Ts<input::IndexInput>,
    ) -> Result<bool, JsValue> {
        self.0
            .remove_bond_between_atoms(
                ck::AtomId::new(input::read(begin)?),
                ck::AtomId::new(input::read(end)?),
            )
            .map_err(OperationError::from_rust)
    }
    pub fn degree(&self, atom: Ts<input::IndexInput>) -> Result<usize, JsValue> {
        self.0
            .degree(ck::AtomId::new(input::read(atom)?))
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=neighborBonds)]
    pub fn neighbor_bonds(
        &self,
        atom: Ts<input::IndexInput>,
    ) -> Result<<Vec<ck::BondId> as Field>::Output, JsValue> {
        let v = self
            .0
            .neighbor_bonds(ck::AtomId::new(input::read(atom)?))
            .map_err(OperationError::from_rust)?;
        Vec::<ck::BondId>::get(NativeField::owned(v))
    }
    #[wasm_bindgen(js_name=bondBetweenAtoms)]
    pub fn bond_between_atoms(
        &self,
        begin: Ts<input::IndexInput>,
        end: Ts<input::IndexInput>,
    ) -> Result<Option<usize>, JsValue> {
        self.0
            .bond_between_atoms(
                ck::AtomId::new(input::read(begin)?),
                ck::AtomId::new(input::read(end)?),
            )
            .map(|v| v.map(ck::BondId::index))
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=set2dCoordinates)]
    pub fn set_2d_coordinates(&mut self, coordinates: Ts<Points2>) -> Result<(), JsValue> {
        self.0
            .set_2d_coordinates(
                input::decode(coordinates)?
                    .0
                    .into_iter()
                    .map(|p| p.0)
                    .collect(),
            )
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=add2dConformer)]
    pub fn add_2d_conformer(&mut self, coordinates: Ts<Points2>) -> Result<usize, JsValue> {
        self.0
            .add_2d_conformer(
                input::decode(coordinates)?
                    .0
                    .into_iter()
                    .map(|p| p.0)
                    .collect(),
            )
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=add3dConformer)]
    pub fn add_3d_conformer(&mut self, coordinates: Ts<Points3>) -> Result<usize, JsValue> {
        self.0
            .add_3d_conformer(
                input::decode(coordinates)?
                    .0
                    .into_iter()
                    .map(|p| p.0)
                    .collect(),
            )
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=addSubstanceGroup)]
    pub fn add_substance_group(&mut self, group: SubstanceGroupInput) -> Result<usize, JsValue> {
        SubstanceGroup::with_ref(&group, |v| {
            self.0
                .add_substance_group(v.clone())
                .map(|v| v.index())
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=addStereoGroup)]
    pub fn add_stereo_group(&mut self, group: StereoGroupInput) -> Result<usize, JsValue> {
        StereoGroup::with_ref(&group, |v| {
            self.0
                .add_stereo_group(v.clone())
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> Result<crate::atoms::Atoms, JsValue> {
        Vec::<ck::Atom>::get(NativeField::owned(self.0.atoms().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn bonds(&self) -> Result<crate::atoms::Bonds, JsValue> {
        Vec::<ck::Bond>::get(NativeField::owned(self.0.bonds().to_vec()))
    }
    #[wasm_bindgen(getter,js_name=stereoGroups)]
    pub fn stereo_groups(&self) -> Result<crate::stereo_groups::StereoGroups, JsValue> {
        Vec::<ck::StereoGroup>::get(NativeField::owned(self.0.stereo_groups().to_vec()))
    }
    #[wasm_bindgen(getter,js_name=substanceGroups)]
    pub fn substance_groups(&self) -> Result<<Vec<ck::SubstanceGroup> as Field>::Output, JsValue> {
        Vec::<ck::SubstanceGroup>::get(NativeField::owned(self.0.substance_groups().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn properties(&self) -> MoleculeProperties {
        MoleculeProperties(self.0.properties().clone())
    }
    #[wasm_bindgen(js_name=withProperties)]
    pub fn with_properties(&self, p: MoleculePropertiesInput) -> Result<Self, JsValue> {
        MoleculeProperties::with_ref(&p, |p| Ok(Self(self.0.clone().with_properties(p.clone()))))
    }
    #[wasm_bindgen(js_name=withName)]
    pub fn with_name(&self, name: Ts<input::TextInput>) -> Result<Self, JsValue> {
        Ok(Self(self.0.clone().with_name(input::read::<String>(name)?)))
    }
    #[wasm_bindgen(js_name=withProperty)]
    pub fn with_property(
        &self,
        key: Ts<input::TextInput>,
        value: Ts<input::TextInput>,
    ) -> Result<Self, JsValue> {
        self.0
            .clone()
            .with_property(input::read::<String>(key)?, input::read::<String>(value)?)
            .map(Self)
            .map_err(OperationError::from_rust)
    }
    #[wasm_bindgen(js_name=withSdfDataField)]
    pub fn with_sdf_data_field(
        &self,
        key: Ts<input::TextInput>,
        value: Ts<input::TextInput>,
    ) -> Result<Self, JsValue> {
        Ok(Self(self.0.clone().with_sdf_data_field(
            input::read::<String>(key)?,
            input::read::<String>(value)?,
        )))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=toBuilder)]
    pub fn to_builder(&self) -> MoleculeBuilder {
        MoleculeBuilder(self.0.to_builder())
    }
}
