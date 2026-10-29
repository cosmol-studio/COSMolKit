//! Read-only detached atom and bond projections; no live mutable storage.
use crate::{Molecule, host_values::usize_value};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum ChiralTag {
    Unspecified = 0,
    TetrahedralCw = 1,
    TetrahedralCcw = 2,
    Other = 3,
    Tetrahedral = 4,
    Allene = 5,
    SquarePlanar = 6,
    TrigonalBipyramidal = 7,
    Octahedral = 8,
}
impl From<ck::ChiralTag> for ChiralTag {
    fn from(value: ck::ChiralTag) -> Self {
        match value {
            ck::ChiralTag::Unspecified => Self::Unspecified,
            ck::ChiralTag::TetrahedralCw => Self::TetrahedralCw,
            ck::ChiralTag::TetrahedralCcw => Self::TetrahedralCcw,
            ck::ChiralTag::Other => Self::Other,
            ck::ChiralTag::Tetrahedral => Self::Tetrahedral,
            ck::ChiralTag::Allene => Self::Allene,
            ck::ChiralTag::SquarePlanar => Self::SquarePlanar,
            ck::ChiralTag::TrigonalBipyramidal => Self::TrigonalBipyramidal,
            ck::ChiralTag::Octahedral => Self::Octahedral,
        }
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum Hybridization {
    Unspecified = 0,
    S = 1,
    Sp = 2,
    Sp2 = 3,
    Sp3 = 4,
    Sp2d = 5,
    Sp3d = 6,
    Sp3d2 = 7,
    Other = 8,
}
impl From<ck::Hybridization> for Hybridization {
    fn from(value: ck::Hybridization) -> Self {
        match value {
            ck::Hybridization::Unspecified => Self::Unspecified,
            ck::Hybridization::S => Self::S,
            ck::Hybridization::Sp => Self::Sp,
            ck::Hybridization::Sp2 => Self::Sp2,
            ck::Hybridization::Sp3 => Self::Sp3,
            ck::Hybridization::Sp2d => Self::Sp2d,
            ck::Hybridization::Sp3d => Self::Sp3d,
            ck::Hybridization::Sp3d2 => Self::Sp3d2,
            ck::Hybridization::Other => Self::Other,
        }
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum BondDirection {
    None = 0,
    BeginWedge = 1,
    BeginDash = 2,
    EndDownRight = 3,
    EndUpRight = 4,
    EitherDouble = 5,
    Unknown = 6,
}
impl From<ck::BondDirection> for BondDirection {
    fn from(value: ck::BondDirection) -> Self {
        match value {
            ck::BondDirection::None => Self::None,
            ck::BondDirection::BeginWedge => Self::BeginWedge,
            ck::BondDirection::BeginDash => Self::BeginDash,
            ck::BondDirection::EndDownRight => Self::EndDownRight,
            ck::BondDirection::EndUpRight => Self::EndUpRight,
            ck::BondDirection::EitherDouble => Self::EitherDouble,
            ck::BondDirection::Unknown => Self::Unknown,
        }
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum BondStereo {
    None = 0,
    Any = 1,
    Z = 2,
    E = 3,
    Cis = 4,
    Trans = 5,
    AtropCw = 6,
    AtropCcw = 7,
}
impl From<ck::BondStereo> for BondStereo {
    fn from(value: ck::BondStereo) -> Self {
        match value {
            ck::BondStereo::None => Self::None,
            ck::BondStereo::Any => Self::Any,
            ck::BondStereo::Z => Self::Z,
            ck::BondStereo::E => Self::E,
            ck::BondStereo::Cis => Self::Cis,
            ck::BondStereo::Trans => Self::Trans,
            ck::BondStereo::AtropCw => Self::AtropCw,
            ck::BondStereo::AtropCcw => Self::AtropCcw,
        }
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum CipDescriptor {
    R = 0,
    S = 1,
    LowerR = 2,
    LowerS = 3,
    E = 4,
    Z = 5,
    LowerE = 6,
    LowerZ = 7,
    M = 8,
    P = 9,
    LowerM = 10,
    LowerP = 11,
}
impl From<ck::CipDescriptor> for CipDescriptor {
    fn from(value: ck::CipDescriptor) -> Self {
        match value {
            ck::CipDescriptor::R => Self::R,
            ck::CipDescriptor::S => Self::S,
            ck::CipDescriptor::LowerR => Self::LowerR,
            ck::CipDescriptor::LowerS => Self::LowerS,
            ck::CipDescriptor::E => Self::E,
            ck::CipDescriptor::Z => Self::Z,
            ck::CipDescriptor::LowerE => Self::LowerE,
            ck::CipDescriptor::LowerZ => Self::LowerZ,
            ck::CipDescriptor::M => Self::M,
            ck::CipDescriptor::P => Self::P,
            ck::CipDescriptor::LowerM => Self::LowerM,
            ck::CipDescriptor::LowerP => Self::LowerP,
        }
    }
}

#[wasm_bindgen]
pub struct Atom {
    pub(crate) inner: ck::Atom,
}
#[wasm_bindgen]
impl Atom {
    pub fn id(&self) -> usize {
        self.inner.id().index()
    }
    #[wasm_bindgen(unchecked_return_type = "Element")]
    pub fn element(&self) -> Result<JsValue, JsValue> {
        crate::element_boundary::Element::from_atomic_number(self.inner.atomic_number() as f64)
    }
    #[wasm_bindgen(js_name = atomicNumber)]
    pub fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    #[wasm_bindgen(js_name = formalCharge)]
    pub fn formal_charge(&self) -> i8 {
        self.inner.formal_charge()
    }
    #[wasm_bindgen(js_name = explicitHydrogens)]
    pub fn explicit_hydrogens(&self) -> u8 {
        self.inner.explicit_hydrogens()
    }
    #[wasm_bindgen(js_name = chiralTagCode)]
    pub fn chiral_tag_code(&self) -> f64 {
        self.inner.chiral_tag_code() as f64
    }
    #[wasm_bindgen(js_name = chiralTagName)]
    pub fn chiral_tag_name(&self) -> String {
        self.inner.chiral_tag_name().into()
    }
    #[wasm_bindgen(js_name = isAromatic)]
    pub fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    #[wasm_bindgen(js_name = noImplicit)]
    pub fn no_implicit(&self) -> bool {
        self.inner.no_implicit()
    }
    #[wasm_bindgen(js_name = radicalElectrons)]
    pub fn radical_electrons(&self) -> u8 {
        self.inner.radical_electrons()
    }
    #[wasm_bindgen(js_name = chiralTag)]
    pub fn chiral_tag(&self) -> ChiralTag {
        self.inner.chiral_tag().into()
    }
    #[wasm_bindgen(js_name = hybridization)]
    pub fn hybridization(&self) -> Hybridization {
        self.inner.hybridization().into()
    }
    #[wasm_bindgen(js_name = isotope, unchecked_return_type = "number | null")]
    pub fn isotope(&self) -> JsValue {
        self.inner
            .isotope()
            .map(u32::from)
            .map_or(JsValue::NULL, JsValue::from)
    }
    #[wasm_bindgen(js_name = atomMap, unchecked_return_type = "number | null")]
    pub fn atom_map(&self) -> JsValue {
        self.inner.atom_map().map_or(JsValue::NULL, JsValue::from)
    }
    #[wasm_bindgen(js_name = cipRank, unchecked_return_type = "number | null")]
    pub fn cip_rank(&self) -> Result<JsValue, JsValue> {
        self.inner
            .cip_rank()
            .map(|v| v.map_or(JsValue::NULL, JsValue::from))
            .map_err(|e| crate::property_values::property_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = cipDescriptor, unchecked_return_type = "CipDescriptor | null")]
    pub fn cip_descriptor(&self) -> Result<JsValue, JsValue> {
        self.inner
            .cip_descriptor()
            .map(|v| {
                v.map_or(JsValue::NULL, |v| {
                    JsValue::from(CipDescriptor::from(v) as u32)
                })
            })
            .map_err(|e| cip_error(e))
    }
    #[wasm_bindgen(js_name = cipNeighborOrder, unchecked_return_type = "number[] | null")]
    pub fn cip_neighbor_order(&self) -> Result<JsValue, JsValue> {
        self.inner
            .cip_neighbor_order()
            .map(|v| {
                v.map_or(JsValue::NULL, |v| {
                    v.into_iter().map(JsValue::from).collect::<Array>().into()
                })
            })
            .map_err(|e| cip_error(e))
    }
    #[wasm_bindgen(js_name = templateAttachmentOrder, unchecked_return_type = "TemplateAttachmentOrder | null")]
    pub fn template_attachment_order(&self) -> JsValue {
        self.inner
            .template_attachment_order()
            .cloned()
            .map_or(JsValue::NULL, |inner| {
                TemplateAttachmentOrder { inner }.into()
            })
    }
}
#[wasm_bindgen]
pub struct Bond {
    pub(crate) inner: ck::Bond,
}
#[wasm_bindgen]
impl Bond {
    pub fn id(&self) -> usize {
        self.inner.id().index()
    }
    pub fn begin(&self) -> usize {
        self.inner.begin().index()
    }
    pub fn end(&self) -> usize {
        self.inner.end().index()
    }
    #[wasm_bindgen(js_name = cipDescriptor, unchecked_return_type = "CipDescriptor | null")]
    pub fn cip_descriptor(&self) -> Result<JsValue, JsValue> {
        self.inner
            .cip_descriptor()
            .map(|v| {
                v.map_or(JsValue::NULL, |v| {
                    JsValue::from(CipDescriptor::from(v) as u32)
                })
            })
            .map_err(|e| cip_error(e))
    }
    #[wasm_bindgen(js_name = cipNeighborOrder, unchecked_return_type = "number[] | null")]
    pub fn cip_neighbor_order(&self) -> Result<JsValue, JsValue> {
        self.inner
            .cip_neighbor_order()
            .map(|v| {
                v.map_or(JsValue::NULL, |v| {
                    v.into_iter().map(JsValue::from).collect::<Array>().into()
                })
            })
            .map_err(|e| cip_error(e))
    }
    #[wasm_bindgen(js_name = orderCode)]
    pub fn order_code(&self) -> f64 {
        self.inner.order_code() as f64
    }
    #[wasm_bindgen(js_name = orderName)]
    pub fn order_name(&self) -> String {
        self.inner.order_name().into()
    }
    #[wasm_bindgen(js_name = directionCode)]
    pub fn direction_code(&self) -> f64 {
        self.inner.direction_code() as f64
    }
    #[wasm_bindgen(js_name = directionName)]
    pub fn direction_name(&self) -> String {
        self.inner.direction_name().into()
    }
    #[wasm_bindgen(js_name = stereoCode)]
    pub fn stereo_code(&self) -> f64 {
        self.inner.stereo_code() as f64
    }
    #[wasm_bindgen(js_name = stereoName)]
    pub fn stereo_name(&self) -> String {
        self.inner.stereo_name().into()
    }
    #[wasm_bindgen(js_name = isAromatic)]
    pub fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    #[wasm_bindgen(js_name = direction)]
    pub fn direction(&self) -> BondDirection {
        self.inner.direction().into()
    }
    #[wasm_bindgen(js_name = stereo)]
    pub fn stereo(&self) -> BondStereo {
        self.inner.stereo().into()
    }
    #[wasm_bindgen(unchecked_return_type = "BondOrder")]
    pub fn order(&self) -> u32 {
        self.inner.order().rdkit_code() as u32
    }
    #[wasm_bindgen(js_name = stereoAtoms, unchecked_return_type = "number[] | null")]
    pub fn stereo_atoms(&self) -> JsValue {
        self.inner.stereo_atoms().map_or(JsValue::NULL, |v| {
            v.into_iter()
                .map(|id| JsValue::from_f64(id.index() as f64))
                .collect::<Array>()
                .into()
        })
    }
}
#[wasm_bindgen]
pub struct AtomMetadata {
    pub(crate) inner: ck::AtomMetadata,
}
#[wasm_bindgen]
impl AtomMetadata {
    #[wasm_bindgen(js_name = degree)]
    pub fn degree(&self) -> usize {
        self.inner.degree
    }
    #[wasm_bindgen(js_name = explicitValence)]
    pub fn explicit_valence(&self) -> i32 {
        self.inner.explicit_valence
    }
    #[wasm_bindgen(js_name = implicitHydrogens)]
    pub fn implicit_hydrogens(&self) -> i32 {
        self.inner.implicit_hydrogens
    }
    #[wasm_bindgen(js_name = totalHydrogens)]
    pub fn total_hydrogens(&self) -> i32 {
        self.inner.total_hydrogens
    }
    #[wasm_bindgen(js_name = totalValence)]
    pub fn total_valence(&self) -> i32 {
        self.inner.total_valence
    }
}
#[wasm_bindgen]
impl Molecule {
    pub fn atoms(&self) -> Vec<Atom> {
        self.inner
            .atoms()
            .into_iter()
            .map(|inner| Atom { inner })
            .collect()
    }
    pub fn bonds(&self) -> Vec<Bond> {
        self.inner
            .bonds()
            .into_iter()
            .map(|inner| Bond { inner })
            .collect()
    }
    #[wasm_bindgen(unchecked_return_type = "Atom | null")]
    pub fn atom(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<JsValue, JsValue> {
        Ok(self
            .inner
            .atom(ck::AtomId::new(usize_value(&index, "atom")?))
            .map_or(JsValue::NULL, |inner| Atom { inner }.into()))
    }
    #[wasm_bindgen(unchecked_return_type = "Bond | null")]
    pub fn bond(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<JsValue, JsValue> {
        Ok(self
            .inner
            .bond(ck::BondId::new(usize_value(&index, "bond")?))
            .map_or(JsValue::NULL, |inner| Bond { inner }.into()))
    }
    #[wasm_bindgen(js_name = atomMetadata)]
    pub fn atom_metadata(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] recalculate: JsValue,
    ) -> Result<Vec<AtomMetadata>, JsValue> {
        let recalculate = if recalculate.is_undefined() {
            true
        } else {
            crate::host_values::bool_value(&recalculate, "recalculate")?
        };
        self.inner
            .atom_metadata(recalculate)
            .map(|v| v.into_iter().map(|inner| AtomMetadata { inner }).collect())
            .map_err(|e| crate::valence_errors::valence_error(&e).unwrap_or_else(|e| e))
    }
}

#[wasm_bindgen]
pub struct CipDescriptorError {
    inner: ck::CipDescriptorError,
}
#[wasm_bindgen]
impl CipDescriptorError {
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        match self.inner {
            ck::CipDescriptorError::InvalidNeighborOrder { .. } => "InvalidNeighborOrder",
            ck::CipDescriptorError::InvalidStoredDescriptor { .. } => "InvalidStoredDescriptor",
            ck::CipDescriptorError::Property(_) => "Property",
        }
        .into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn value(&self) -> Result<JsValue, JsValue> {
        match &self.inner {
            ck::CipDescriptorError::InvalidNeighborOrder { value, .. }
            | ck::CipDescriptorError::InvalidStoredDescriptor { value } => {
                crate::host_values::text(value).map(JsValue::from)
            }
            _ => Ok(JsValue::NULL),
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        use std::error::Error;
        self.inner
            .source()
            .map_or(Ok(JsValue::NULL), crate::alignment_values::source_error)
    }
}
fn cip_error(inner: ck::CipDescriptorError) -> JsValue {
    let detail = CipDescriptorError { inner };
    let error: JsValue = js_sys::Error::new(&detail.message()).into();
    let result = (|| {
        crate::alignment_values::set(&error, "name", "CipDescriptorError".into())?;
        crate::alignment_values::set(&error, "domain", "cip".into())?;
        crate::alignment_values::set(&error, "kind", detail.kind().into())?;
        crate::alignment_values::set(&error, "cause", detail.cause()?)?;
        crate::alignment_values::set(&error, "detail", detail.into())?;
        Ok::<_, JsValue>(error)
    })();
    result.unwrap_or_else(|e| e)
}
#[wasm_bindgen]
pub struct TemplateAttachment {
    pub(crate) inner: ck::TemplateAttachment,
}
#[wasm_bindgen]
impl TemplateAttachment {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_param_type = "number")] target: JsValue,
        label: String,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::TemplateAttachment::new(
                ck::AtomId::new(usize_value(&target, "target")?),
                label,
            ),
        })
    }
    pub fn target(&self) -> usize {
        self.inner.target().index()
    }
    pub fn label(&self) -> String {
        self.inner.label().into()
    }
}
#[wasm_bindgen(inline_js = "export function visitAttachment(value, visit) { visit(value); }")]
extern "C" {
    #[wasm_bindgen(catch, js_name = visitAttachment)]
    fn visit_attachment(
        value: &JsValue,
        visit: &mut dyn FnMut(&TemplateAttachment),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen]
pub struct TemplateAttachmentOrder {
    pub(crate) inner: ck::TemplateAttachmentOrder,
}
#[wasm_bindgen]
impl TemplateAttachmentOrder {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_param_type = "TemplateAttachment[]")] entries: JsValue,
    ) -> Result<Self, JsValue> {
        let entries = crate::host_values::sequence(&entries, "entries")?
            .iter()
            .map(|entry| {
                let mut result = None;
                visit_attachment(&entry, &mut |entry: &TemplateAttachment| {
                    result = Some(entry.inner.clone());
                })?;
                result.ok_or_else(|| crate::host_values::type_error("entry"))
            })
            .collect::<Result<Vec<_>, JsValue>>()?;
        ck::TemplateAttachmentOrder::new(entries)
            .map(|inner| Self { inner })
            .map_err(attachment_error)
    }
    pub fn entries(&self) -> Vec<TemplateAttachment> {
        self.inner
            .entries()
            .iter()
            .cloned()
            .map(|inner| TemplateAttachment { inner })
            .collect()
    }
}
#[wasm_bindgen]
pub struct TemplateAttachmentOrderError {
    inner: ck::TemplateAttachmentOrderError,
}
#[wasm_bindgen]
impl TemplateAttachmentOrderError {
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        use ck::TemplateAttachmentOrderError as E;
        match self.inner {
            E::Empty => "Empty",
            E::DuplicateTarget { .. } => "DuplicateTarget",
            E::DuplicateLabel { .. } => "DuplicateLabel",
            E::TargetOutOfRange { .. } => "TargetOutOfRange",
            E::MappingTargetOutOfRange { .. } => "MappingTargetOutOfRange",
            E::TargetRemoved { .. } => "TargetRemoved",
        }
        .into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, js_name = firstPosition, unchecked_return_type = "number | null")]
    pub fn first_position(&self) -> JsValue {
        use ck::TemplateAttachmentOrderError as E;
        match &self.inner {
            E::DuplicateTarget { first_position, .. }
            | E::DuplicateLabel { first_position, .. } => JsValue::from_f64(*first_position as f64),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, js_name = duplicatePosition, unchecked_return_type = "number | null")]
    pub fn duplicate_position(&self) -> JsValue {
        use ck::TemplateAttachmentOrderError as E;
        match &self.inner {
            E::DuplicateTarget {
                duplicate_position, ..
            }
            | E::DuplicateLabel {
                duplicate_position, ..
            } => JsValue::from_f64(*duplicate_position as f64),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, js_name = position, unchecked_return_type = "number | null")]
    pub fn position(&self) -> JsValue {
        use ck::TemplateAttachmentOrderError as E;
        match &self.inner {
            E::TargetOutOfRange { position, .. }
            | E::MappingTargetOutOfRange { position, .. }
            | E::TargetRemoved { position, .. } => JsValue::from_f64(*position as f64),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, js_name = atomCount, unchecked_return_type = "number | null")]
    pub fn atom_count(&self) -> JsValue {
        use ck::TemplateAttachmentOrderError as E;
        match &self.inner {
            E::TargetOutOfRange { atom_count, .. } => JsValue::from_f64(*atom_count as f64),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, js_name = mappingLen, unchecked_return_type = "number | null")]
    pub fn mapping_len(&self) -> JsValue {
        use ck::TemplateAttachmentOrderError as E;
        match &self.inner {
            E::MappingTargetOutOfRange { mapping_len, .. } => {
                JsValue::from_f64(*mapping_len as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn target(&self) -> JsValue {
        use ck::TemplateAttachmentOrderError as E;
        match self.inner {
            E::DuplicateTarget { target, .. }
            | E::TargetOutOfRange { target, .. }
            | E::MappingTargetOutOfRange { target, .. }
            | E::TargetRemoved { target, .. } => JsValue::from_f64(target.index() as f64),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn label(&self) -> JsValue {
        match &self.inner {
            ck::TemplateAttachmentOrderError::DuplicateLabel { label, .. } => label.into(),
            _ => JsValue::NULL,
        }
    }
}
fn attachment_error(inner: ck::TemplateAttachmentOrderError) -> JsValue {
    let detail = TemplateAttachmentOrderError { inner };
    let error: JsValue = js_sys::Error::new(&detail.message()).into();
    let result = (|| {
        use crate::alignment_values::set;
        set(&error, "name", "TemplateAttachmentOrderError".into())?;
        set(&error, "domain", "template_attachment".into())?;
        set(&error, "kind", detail.kind().into())?;
        set(&error, "firstPosition", detail.first_position())?;
        set(&error, "duplicatePosition", detail.duplicate_position())?;
        set(&error, "target", detail.target())?;
        set(&error, "label", detail.label())?;
        set(&error, "detail", detail.into())?;
        Ok::<_, JsValue>(error)
    })();
    result.unwrap_or_else(|e| e)
}
