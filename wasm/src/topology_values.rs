//! A topology view retains native molecular COW storage until explicitly copied.
use crate::{Molecule, fields::Field, storage::NativeField};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
enum TopologyOwner {
    Molecule(ck::Molecule),
    Detached(ck::TopologyBlock),
}
impl std::ops::Deref for TopologyOwner {
    type Target = ck::TopologyBlock;
    fn deref(&self) -> &ck::TopologyBlock {
        match self {
            Self::Molecule(m) => m.topology(),
            Self::Detached(t) => t,
        }
    }
}
#[wasm_bindgen]
pub struct TopologyBlock(TopologyOwner);
crate::fields::value_field!(
    ck::TopologyBlock,
    TopologyBlock,
    TopologyBlockInput,
    "TopologyBlock",
    TopologyBlock::owned
);
impl TopologyBlock {
    fn owned(value: ck::TopologyBlock) -> Self {
        Self(TopologyOwner::Detached(value))
    }
    fn native(&self) -> &ck::TopologyBlock {
        &self.0
    }
}
#[wasm_bindgen]
impl TopologyBlock {
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> Result<crate::atoms::Atoms, JsValue> {
        Vec::<ck::Atom>::get(NativeField::owned(self.native().atoms.clone()))
    }
    #[wasm_bindgen(getter)]
    pub fn bonds(&self) -> Result<crate::atoms::Bonds, JsValue> {
        Vec::<ck::Bond>::get(NativeField::owned(self.native().bonds.clone()))
    }
    #[wasm_bindgen(getter,js_name=substanceGroups)]
    pub fn substance_groups(&self) -> Result<<Vec<ck::SubstanceGroup> as Field>::Output, JsValue> {
        Vec::<ck::SubstanceGroup>::get(NativeField::owned(self.native().substance_groups.clone()))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(getter)]
    pub fn topology(&self) -> TopologyBlock {
        TopologyBlock(TopologyOwner::Molecule(self.0.clone()))
    }
    #[wasm_bindgen(js_name=fromParts)]
    pub fn from_parts(
        topology: TopologyBlockInput,
        coordinates: crate::detached_coordinates::CoordinateBlockInput,
        properties: crate::metadata_values::MoleculePropertiesInput,
    ) -> Result<Self, JsValue> {
        TopologyBlock::with_ref(&topology, |t| {
            crate::detached_coordinates::CoordinateBlock::with_ref(&coordinates, |c| {
                crate::metadata_values::MoleculeProperties::with_ref(&properties, |p| {
                    ck::Molecule::from_parts(t.clone(), c.clone(), p.clone())
                        .map(Self)
                        .map_err(crate::operation_errors::OperationError::from_rust)
                })
            })
        })
    }
}
