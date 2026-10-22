//! Detached construction snapshots and checked canonical builder delegation.
use crate::atom_bond::{Atom, Bond};
use crate::group_values::SubstanceGroup;
use crate::property_values::{MoleculeProperties, PropertyValue};
use crate::{Molecule, host_values::*};
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Map};
use wasm_bindgen::prelude::*;

fn operation<T>(result: Result<T, ck::OperationError>) -> Result<T, JsValue> {
    result.map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
}
fn ids(value: &JsValue, name: &str) -> Result<Vec<usize>, JsValue> {
    sequence(value, name)?
        .iter()
        .map(|v| usize_value(&v, name))
        .collect()
}
fn array_ids(rows: impl Iterator<Item = usize>) -> JsValue {
    rows.map(|v| JsValue::from(v as u32))
        .collect::<Array>()
        .into()
}
pub(crate) fn matrix(value: &JsValue) -> Result<Vec<Vec<f64>>, JsValue> {
    sequence(value, "coordinates")?
        .iter()
        .map(|row| {
            sequence(&row, "coordinate row")?
                .iter()
                .map(|v| v.as_f64().ok_or_else(|| type_error("coordinate component")))
                .collect()
        })
        .collect()
}
pub(crate) fn fixed_matrix<const D: usize>(value: &JsValue) -> Result<Vec<[f64; D]>, JsValue> {
    sequence(value, "coordinates")?
        .iter()
        .map(|row| {
            let row = sequence(&row, "coordinate row")?;
            if row.length() as usize != D {
                return Err(js_sys::RangeError::new("incorrect coordinate row width").into());
            }
            let mut values = [0.0; D];
            for (index, value) in values.iter_mut().enumerate() {
                *value = row.get(index as u32).as_f64()
                    .ok_or_else(|| type_error("coordinate component"))?;
            }
            Ok(values)
        })
        .collect()
}
fn coordinate_array<const D: usize>(rows: &[[f64; D]]) -> JsValue {
    rows.iter()
        .map(|r| r.iter().map(|v| JsValue::from(*v)).collect::<Array>())
        .collect::<Array>()
        .into()
}
fn props_map(
    props: &std::collections::BTreeMap<ck::PropertyText, ck::PropertyText>,
) -> Result<Map, JsValue> {
    let out = Map::new();
    for (k, v) in props {
        out.set(&text(k)?.into(), &text(v)?.into());
    }
    Ok(out)
}

#[wasm_bindgen]
pub struct AtomSpec {
    inner: ck::AtomSpec,
}
#[wasm_bindgen]
impl AtomSpec {
    #[wasm_bindgen(constructor)]
    pub fn new(element: &crate::element_boundary::Element) -> Self {
        Self {
            inner: ck::AtomSpec::new(element.inner.into()),
        }
    }
    #[wasm_bindgen(js_name=withFormalCharge)]
    pub fn with_formal_charge(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: self.inner.clone().with_formal_charge(integer(
                &value,
                "formalCharge",
                i8::MIN as f64,
                i8::MAX as f64,
            )? as i8),
        })
    }
    #[wasm_bindgen(js_name=withExplicitHydrogens)]
    pub fn with_explicit_hydrogens(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: self.inner.clone().with_explicit_hydrogens(integer(
                &value,
                "explicitHydrogens",
                0.,
                u8::MAX as f64,
            )? as u8),
        })
    }
    #[wasm_bindgen(js_name=withAtomMap)]
    pub fn with_atom_map(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: self
                .inner
                .clone()
                .with_atom_map(u32_value(&value, "atomMap")?),
        })
    }
    #[wasm_bindgen(js_name=withIsotope)]
    pub fn with_isotope(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: self
                .inner
                .clone()
                .with_isotope(integer(&value, "isotope", 0., u16::MAX as f64)? as u16),
        })
    }
    #[wasm_bindgen(js_name=withNoImplicit)]
    pub fn with_no_implicit(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: self
                .inner
                .clone()
                .with_no_implicit(bool_value(&value, "noImplicit")?),
        })
    }
}
pub(crate) fn bond_order(value: &JsValue) -> Result<ck::BondOrder, JsValue> {
    if let Some(name) = value.as_string() {
        return ck::BondOrder::from_rdkit_name(&name).ok_or_else(|| type_error("BondOrder"));
    }
    ck::BondOrder::from_rdkit_code(integer(value, "BondOrder", 0., i64::MAX as f64)? as i64)
        .ok_or_else(|| js_sys::RangeError::new("invalid BondOrder").into())
}
#[wasm_bindgen]
pub struct BondSpec {
    inner: ck::BondSpec,
}
#[wasm_bindgen]
impl BondSpec {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_param_type = "number")] begin: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] end: JsValue,
        #[wasm_bindgen(unchecked_param_type = "BondOrder | string")] order: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::BondSpec::new(
                ck::AtomId::new(usize_value(&begin, "begin")?),
                ck::AtomId::new(usize_value(&end, "end")?),
                bond_order(&order)?,
            ),
        })
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum StereoGroupKind {
    Absolute,
    Or,
    And,
}
#[wasm_bindgen]
pub struct StereoGroup {
    inner: ck::StereoGroup,
}
#[wasm_bindgen]
impl StereoGroup {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_param_type = "StereoGroupKind | string")] kind: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number[]")] atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[]")] bonds: JsValue,
    ) -> Result<Self, JsValue> {
        let kind = if let Some(s) = kind.as_string() {
            match s.as_str() {
                "Absolute" => ck::StereoGroupKind::Absolute,
                "Or" => ck::StereoGroupKind::Or,
                "And" => ck::StereoGroupKind::And,
                _ => return Err(type_error("StereoGroupKind")),
            }
        } else {
            match integer(&kind, "StereoGroupKind", 0., 2.)? as u8 {
                0 => ck::StereoGroupKind::Absolute,
                1 => ck::StereoGroupKind::Or,
                2 => ck::StereoGroupKind::And,
                _ => unreachable!(),
            }
        };
        let atoms = ids(&atoms, "atoms")?
            .into_iter()
            .map(ck::AtomId::new)
            .collect();
        let bonds = if bonds.is_undefined() {
            Vec::new()
        } else {
            ids(&bonds, "bonds")?
                .into_iter()
                .map(ck::BondId::new)
                .collect()
        };
        ck::StereoGroup::new(kind, atoms, bonds)
            .map(|inner| Self { inner })
            .map_err(|e| js_sys::RangeError::new(&e.to_string()).into())
    }
    pub fn kind(&self) -> StereoGroupKind {
        match self.inner.kind() {
            ck::StereoGroupKind::Absolute => StereoGroupKind::Absolute,
            ck::StereoGroupKind::Or => StereoGroupKind::Or,
            ck::StereoGroupKind::And => StereoGroupKind::And,
        }
    }
    #[wasm_bindgen(unchecked_return_type = "number[]")]
    pub fn atoms(&self) -> JsValue {
        array_ids(self.inner.atoms().iter().map(|id| id.index()))
    }
    #[wasm_bindgen(unchecked_return_type = "number[]")]
    pub fn bonds(&self) -> JsValue {
        array_ids(self.inner.bonds().iter().map(|id| id.index()))
    }
}
fn atoms(rows: &[ck::Atom]) -> JsValue {
    rows.iter()
        .cloned()
        .map(|inner| JsValue::from(Atom { inner }))
        .collect::<Array>()
        .into()
}
fn bonds(rows: &[ck::Bond]) -> JsValue {
    rows.iter()
        .cloned()
        .map(|inner| JsValue::from(Bond { inner }))
        .collect::<Array>()
        .into()
}
fn groups(rows: &[ck::SubstanceGroup]) -> JsValue {
    rows.iter()
        .cloned()
        .map(|inner| JsValue::from(SubstanceGroup { inner }))
        .collect::<Array>()
        .into()
}
fn stereo_groups(rows: &[ck::StereoGroup]) -> JsValue {
    rows.iter()
        .cloned()
        .map(|inner| JsValue::from(StereoGroup { inner }))
        .collect::<Array>()
        .into()
}

#[wasm_bindgen]
pub struct TopologyBlock {
    inner: ck::TopologyBlock,
}
#[wasm_bindgen]
impl TopologyBlock {
    #[wasm_bindgen(getter, unchecked_return_type = "Atom[]")]
    pub fn atoms(&self) -> JsValue {
        atoms(&self.inner.atoms)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Bond[]")]
    pub fn bonds(&self) -> JsValue {
        bonds(&self.inner.bonds)
    }
    #[wasm_bindgen(getter,js_name=substanceGroups,unchecked_return_type="SubstanceGroup[]")]
    pub fn substance_groups(&self) -> JsValue {
        groups(&self.inner.substance_groups)
    }
    #[wasm_bindgen(getter,js_name=stereoGroups,unchecked_return_type="StereoGroup[]")]
    pub fn stereo_groups(&self) -> JsValue {
        stereo_groups(&self.inner.stereo_groups)
    }
}
#[wasm_bindgen]
pub struct Conformer2D {
    inner: ck::Conformer2D,
}
#[wasm_bindgen]
impl Conformer2D {
    pub fn id(&self) -> usize {
        self.inner.id()
    }
    #[wasm_bindgen(unchecked_return_type = "[number, number][]")]
    pub fn coordinates(&self) -> JsValue {
        coordinate_array(self.inner.coordinates())
    }
    #[wasm_bindgen(unchecked_return_type = "Map<string, string>")]
    pub fn props(&self) -> Result<Map, JsValue> {
        props_map(self.inner.props())
    }
}
#[wasm_bindgen]
pub struct Conformer3D {
    inner: ck::Conformer3D,
}
#[wasm_bindgen]
impl Conformer3D {
    pub fn id(&self) -> usize {
        self.inner.id()
    }
    #[wasm_bindgen(js_name=is3d)]
    pub fn is_3d(&self) -> bool {
        self.inner.is_3d()
    }
    #[wasm_bindgen(unchecked_return_type = "[number, number, number][]")]
    pub fn coordinates(&self) -> JsValue {
        coordinate_array(self.inner.coordinates())
    }
    #[wasm_bindgen(unchecked_return_type = "Map<string, string>")]
    pub fn props(&self) -> Result<Map, JsValue> {
        props_map(self.inner.props())
    }
}
#[wasm_bindgen]
pub struct CoordinateBlock {
    inner: ck::CoordinateBlock,
}
#[wasm_bindgen]
impl CoordinateBlock {
    #[wasm_bindgen(getter,js_name=conformers2d,unchecked_return_type="Conformer2D[]")]
    pub fn conformers_2d(&self) -> JsValue {
        self.inner
            .conformers_2d
            .iter()
            .cloned()
            .map(|inner| JsValue::from(Conformer2D { inner }))
            .collect::<Array>()
            .into()
    }
    #[wasm_bindgen(getter,js_name=conformers3d,unchecked_return_type="Conformer3D[]")]
    pub fn conformers_3d(&self) -> JsValue {
        conformers_3d(&self.inner.conformers_3d)
    }
    #[wasm_bindgen(getter,js_name=sourceCoordinateDim,unchecked_return_type="CoordinateDimension | null")]
    pub fn source_coordinate_dim(&self) -> JsValue {
        self.inner.source_coordinate_dim.map_or(JsValue::NULL, |v| {
            JsValue::from(match v {
                ck::CoordinateDimension::TwoD => 0u32,
                ck::CoordinateDimension::ThreeD => 1u32,
            })
        })
    }
    #[wasm_bindgen(getter,js_name=sourceConformerOrder,unchecked_return_type="CoordinateDimension[] | null")]
    pub fn source_conformer_order(&self) -> JsValue {
        self.inner
            .source_conformer_order
            .as_ref()
            .map_or(JsValue::NULL, |rows| {
                rows.iter()
                    .map(|v| {
                        JsValue::from(match v {
                            ck::CoordinateDimension::TwoD => 0u32,
                            ck::CoordinateDimension::ThreeD => 1u32,
                        })
                    })
                    .collect::<Array>()
                    .into()
            })
    }
}
fn conformers_3d(rows: &[ck::Conformer3D]) -> JsValue {
    rows.iter()
        .cloned()
        .map(|inner| JsValue::from(Conformer3D { inner }))
        .collect::<Array>()
        .into()
}

#[wasm_bindgen]
pub struct MoleculeBuilder {
    inner: ck::MoleculeBuilder,
}
#[wasm_bindgen]
impl MoleculeBuilder {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            inner: ck::MoleculeBuilder::new(),
        }
    }
    #[wasm_bindgen(js_name=fromParts)]
    pub fn from_parts(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        properties: &MoleculeProperties,
    ) -> Self {
        Self {
            inner: ck::MoleculeBuilder::from_parts(
                topology.inner.clone(),
                coordinates.inner.clone(),
                properties.inner.clone(),
            ),
        }
    }
    pub fn build(&self) -> Result<Molecule, JsValue> {
        operation(self.inner.clone().build()).map(|inner| Molecule {
            inner: cosmolkit_wasm::Molecule::from(inner).into(),
        })
    }
    #[wasm_bindgen(js_name=addAtom)]
    pub fn add_atom(&mut self, atom: &AtomSpec) -> usize {
        self.inner.add_atom(atom.inner.clone()).index()
    }
    #[wasm_bindgen(js_name=addBond)]
    pub fn add_bond(&mut self, bond: &BondSpec) -> Result<usize, JsValue> {
        operation(self.inner.add_bond(bond.inner.clone())).map(|id| id.index())
    }
    #[wasm_bindgen(js_name=setAtomFormalCharge)]
    pub fn set_atom_formal_charge(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] charge: JsValue,
    ) -> Result<(), JsValue> {
        operation(self.inner.set_atom_formal_charge(
            ck::AtomId::new(usize_value(&atom, "atom")?),
            integer(&charge, "formalCharge", i8::MIN as f64, i8::MAX as f64)? as i8,
        ))
    }
    #[wasm_bindgen(js_name=setBondOrder)]
    pub fn set_bond_order(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] bond: JsValue,
        #[wasm_bindgen(unchecked_param_type = "BondOrder | string")] order: JsValue,
    ) -> Result<(), JsValue> {
        operation(self.inner.set_bond_order(
            ck::BondId::new(usize_value(&bond, "bond")?),
            bond_order(&order)?,
        ))
    }
    #[wasm_bindgen(js_name=removeBondBetweenAtoms)]
    pub fn remove_bond_between_atoms(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] begin: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] end: JsValue,
    ) -> Result<bool, JsValue> {
        operation(self.inner.remove_bond_between_atoms(
            ck::AtomId::new(usize_value(&begin, "begin")?),
            ck::AtomId::new(usize_value(&end, "end")?),
        ))
    }
    pub fn degree(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
    ) -> Result<usize, JsValue> {
        operation(
            self.inner
                .degree(ck::AtomId::new(usize_value(&atom, "atom")?)),
        )
    }
    #[wasm_bindgen(js_name=neighborBonds,unchecked_return_type="number[]")]
    pub fn neighbor_bonds(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
    ) -> Result<JsValue, JsValue> {
        operation(
            self.inner
                .neighbor_bonds(ck::AtomId::new(usize_value(&atom, "atom")?)),
        )
        .map(|v| array_ids(v.iter().map(|id| id.index())))
    }
    #[wasm_bindgen(js_name=bondBetweenAtoms,unchecked_return_type="number | null")]
    pub fn bond_between_atoms(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] begin: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] end: JsValue,
    ) -> Result<JsValue, JsValue> {
        operation(self.inner.bond_between_atoms(
            ck::AtomId::new(usize_value(&begin, "begin")?),
            ck::AtomId::new(usize_value(&end, "end")?),
        ))
        .map(|id| id.map_or(JsValue::NULL, |id| JsValue::from(id.index() as u32)))
    }
    #[wasm_bindgen(unchecked_return_type = "Atom[]")]
    pub fn atoms(&self) -> JsValue {
        atoms(self.inner.atoms())
    }
    #[wasm_bindgen(unchecked_return_type = "Bond[]")]
    pub fn bonds(&self) -> JsValue {
        bonds(self.inner.bonds())
    }
    #[wasm_bindgen(js_name=substanceGroups,unchecked_return_type="SubstanceGroup[]")]
    pub fn substance_groups(&self) -> JsValue {
        groups(self.inner.substance_groups())
    }
    #[wasm_bindgen(js_name=stereoGroups,unchecked_return_type="StereoGroup[]")]
    pub fn stereo_groups(&self) -> JsValue {
        stereo_groups(self.inner.stereo_groups())
    }
    pub fn coordinates(&self) -> CoordinateBlock {
        CoordinateBlock {
            inner: self.inner.coordinates().clone(),
        }
    }
    pub fn properties(&self) -> MoleculeProperties {
        MoleculeProperties {
            inner: self.inner.properties().clone(),
        }
    }
    #[wasm_bindgen(js_name=set2dCoordinates)]
    pub fn set_2d_coordinates(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
    ) -> Result<(), JsValue> {
        operation(self.inner.set_2d_coordinates(fixed_matrix(&coordinates)?))
    }
    #[wasm_bindgen(js_name=add2dConformer)]
    pub fn add_2d_conformer(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
    ) -> Result<usize, JsValue> {
        operation(self.inner.add_2d_conformer(fixed_matrix(&coordinates)?))
    }
    #[wasm_bindgen(js_name=add3dConformer)]
    pub fn add_3d_conformer(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
    ) -> Result<usize, JsValue> {
        operation(self.inner.add_3d_conformer(fixed_matrix(&coordinates)?))
    }
    #[wasm_bindgen(js_name=addSubstanceGroup)]
    pub fn add_substance_group(&mut self, group: &SubstanceGroup) -> Result<usize, JsValue> {
        operation(self.inner.add_substance_group(group.inner.clone())).map(|id| id.index())
    }
    #[wasm_bindgen(js_name=addStereoGroup)]
    pub fn add_stereo_group(&mut self, group: &StereoGroup) -> Result<usize, JsValue> {
        operation(self.inner.add_stereo_group(group.inner.clone()))
    }
    #[wasm_bindgen(js_name=withName)]
    pub fn with_name(&self, name: String) -> Self {
        Self {
            inner: self.inner.clone().with_name(name),
        }
    }
    #[wasm_bindgen(js_name=withProperty)]
    pub fn with_property(&self, key: String, value: String) -> Result<Self, JsValue> {
        operation(self.inner.clone().with_property(key, value)).map(|inner| Self { inner })
    }
    #[wasm_bindgen(js_name=withSdfDataField)]
    pub fn with_sdf_data_field(&self, key: String, value: String) -> Self {
        Self {
            inner: self.inner.clone().with_sdf_data_field(key, value),
        }
    }
    #[wasm_bindgen(js_name=withProperties)]
    pub fn with_properties(&self, properties: &MoleculeProperties) -> Self {
        Self {
            inner: self.inner.clone().with_properties(properties.inner.clone()),
        }
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fromParts)]
    pub fn from_parts(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        properties: &MoleculeProperties,
    ) -> Result<Molecule, JsValue> {
        operation(cosmolkit_wasm::Molecule::from_parts(
            topology.inner.clone(),
            coordinates.inner.clone(),
            properties.inner.clone(),
        ))
        .map(|inner| Molecule {
            inner: inner.into(),
        })
    }
    #[wasm_bindgen(js_name=toBuilder)]
    pub fn to_builder(&self) -> MoleculeBuilder {
        MoleculeBuilder {
            inner: self.inner.to_builder(),
        }
    }
    pub fn topology(&self) -> TopologyBlock {
        TopologyBlock {
            inner: self.inner.topology(),
        }
    }
    pub fn properties(&self) -> MoleculeProperties {
        MoleculeProperties {
            inner: self.inner.properties(),
        }
    }
    #[wasm_bindgen(unchecked_return_type = "PropertyValue | null")]
    pub fn property(&self, key: &str) -> JsValue {
        self.inner
            .property(key)
            .map_or(JsValue::NULL, |inner| PropertyValue { inner }.into())
    }
    #[wasm_bindgen(js_name=conformers3d,unchecked_return_type="Conformer3D[]")]
    pub fn conformers_3d(&self) -> JsValue {
        conformers_3d(&self.inner.conformers_3d())
    }
    #[wasm_bindgen(js_name=has2dCoordinates)]
    pub fn has_2d_coordinates(&self) -> bool {
        self.inner.has_2d_coordinates()
    }
    #[wasm_bindgen(js_name=hasValenceViolation)]
    pub fn has_valence_violation(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom_id: JsValue,
    ) -> Result<bool, JsValue> {
        self.inner
            .has_valence_violation(ck::AtomId::new(usize_value(&atom_id, "atomId")?))
            .map_err(|e| crate::valence_errors::valence_error(&e).unwrap_or_else(|e| e))
    }
}
