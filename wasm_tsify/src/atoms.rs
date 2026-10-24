//! Retained molecular atom values and detached atom/bond reads.
use crate::enums::{BondDirection, BondOrder, BondStereo, ChiralTag, CipDescriptor, Hybridization};
use crate::input;
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "readonly [number, number]")]
    pub type StereoAtomPair;
}

macro_rules! index_field {
    ($native:ty) => {
        impl crate::fields::Field for $native {
            type Input = Ts<input::IndexInput>;
            type Output = usize;
            fn input(value: Self::Input) -> Result<Self, JsValue> {
                input::read::<usize>(value).map(Self::new)
            }
            fn decode(value: &JsValue) -> Result<Self, JsValue> {
                Self::input(Ts::new_unchecked(value.clone()))
            }
            fn get(value: crate::storage::NativeField<Self>) -> Result<usize, JsValue> {
                value.read(|value| value.index())
            }
        }
        impl crate::fields::Field for Vec<$native> {
            type Input = <Vec<usize> as crate::fields::Field>::Input;
            type Output = <Vec<usize> as crate::fields::Field>::Output;
            fn input(value: Self::Input) -> Result<Self, JsValue> {
                <Vec<usize> as crate::fields::Field>::input(value)
                    .map(|values| values.into_iter().map(<$native>::new).collect())
            }
            fn decode(value: &JsValue) -> Result<Self, JsValue> {
                Self::input(Ts::new_unchecked(value.clone()))
            }
            fn get(value: crate::storage::NativeField<Self>) -> Result<Self::Output, JsValue> {
                let indices =
                    value.read(|values| values.iter().map(|value| value.index()).collect())?;
                <Vec<usize> as crate::fields::Field>::get(crate::storage::NativeField::owned(
                    indices,
                ))
            }
        }
    };
}
index_field!(ck::AtomId);
index_field!(ck::BondId);

#[wasm_bindgen]
#[derive(Clone)]
pub struct Element(pub(crate) ck::Element);
crate::fields::value_field!(ck::Element, Element, ElementInput, "Element");

#[wasm_bindgen]
impl Element {
    #[wasm_bindgen(js_name = fromSymbol)]
    pub fn from_symbol(symbol: tsify::Ts<crate::input::TextInput>) -> Result<Element, JsValue> {
        crate::input::read::<String>(symbol)?
            .parse::<ck::Element>()
            .map(Self)
            .map_err(crate::errors::ElementParseError::from_rust)
    }
    #[wasm_bindgen(js_name = fromAtomicNumber)]
    pub fn from_atomic_number(
        number: tsify::Ts<crate::input::U8Input>,
    ) -> Result<Option<Element>, JsValue> {
        Ok(ck::Element::from_atomic_number(crate::input::read::<u8>(number)?).map(Self))
    }
    #[wasm_bindgen(getter, js_name = atomicNumber)]
    pub fn atomic_number(&self) -> u8 {
        self.0.atomic_number()
    }
    #[wasm_bindgen(getter)]
    pub fn symbol(&self) -> String {
        self.0.symbol().to_owned()
    }
}

#[wasm_bindgen]
#[derive(Clone)]
pub struct Atom(pub(crate) ck::Atom);
crate::fields::value_field!(ck::Atom, Atom, AtomInput, "Atom");
crate::fields::value_array!(ck::Atom, Atoms, "ReadonlyArray<Atom>");

#[wasm_bindgen]
#[derive(Clone)]
pub struct MoleculeAtom(pub(crate) ck::MoleculeAtom);

macro_rules! atom_accessors {
    ($name:ident) => {
        #[wasm_bindgen]
        impl $name {
            #[wasm_bindgen(getter)]
            pub fn id(&self) -> usize {
                self.0.id().index()
            }
            #[wasm_bindgen(getter)]
            pub fn element(&self) -> Element {
                Element(self.0.element())
            }
            #[wasm_bindgen(getter, js_name = atomicNumber)]
            pub fn atomic_number(&self) -> u8 {
                self.0.atomic_number()
            }
            #[wasm_bindgen(getter, js_name = formalCharge)]
            pub fn formal_charge(&self) -> i8 {
                self.0.formal_charge()
            }
            #[wasm_bindgen(getter, js_name = chiralTag)]
            pub fn chiral_tag(&self) -> ChiralTag {
                self.0.chiral_tag().into()
            }
            #[wasm_bindgen(getter, js_name = chiralTagCode)]
            pub fn chiral_tag_code(&self) -> i64 {
                self.0.chiral_tag_code()
            }
            #[wasm_bindgen(getter, js_name = chiralTagName)]
            pub fn chiral_tag_name(&self) -> String {
                self.0.chiral_tag_name().to_owned()
            }
            #[wasm_bindgen(getter)]
            pub fn isotope(&self) -> Option<u16> {
                self.0.isotope()
            }
            #[wasm_bindgen(getter, js_name = atomMap)]
            pub fn atom_map(&self) -> Option<u32> {
                self.0.atom_map()
            }
            #[wasm_bindgen(getter, js_name = isAromatic)]
            pub fn is_aromatic(&self) -> bool {
                self.0.is_aromatic()
            }
            #[wasm_bindgen(getter, js_name = explicitHydrogens)]
            pub fn explicit_hydrogens(&self) -> u8 {
                self.0.explicit_hydrogens()
            }
            #[wasm_bindgen(getter, js_name = noImplicit)]
            pub fn no_implicit(&self) -> bool {
                self.0.no_implicit()
            }
            #[wasm_bindgen(getter, js_name = radicalElectrons)]
            pub fn radical_electrons(&self) -> u8 {
                self.0.radical_electrons()
            }
            #[wasm_bindgen(getter)]
            pub fn hybridization(&self) -> Hybridization {
                self.0.hybridization().into()
            }
            #[wasm_bindgen(getter, js_name = cipDescriptor)]
            pub fn cip_descriptor(&self) -> Result<Option<CipDescriptor>, JsValue> {
                self.0
                    .cip_descriptor()
                    .map(|value| value.map(Into::into))
                    .map_err(crate::model_errors::CipDescriptorError::from_rust)
            }
            #[wasm_bindgen(getter, js_name = cipNeighborOrder)]
            pub fn cip_neighbor_order(&self) -> Result<Option<Vec<u32>>, JsValue> {
                self.0
                    .cip_neighbor_order()
                    .map_err(crate::model_errors::CipDescriptorError::from_rust)
            }
            #[wasm_bindgen(getter, js_name = cipRank)]
            pub fn cip_rank(&self) -> Result<Option<u32>, JsValue> {
                self.0
                    .cip_rank()
                    .map_err(crate::model_errors::PropertyValueError::from_rust)
            }
            #[wasm_bindgen(getter,js_name=templateAttachmentOrder)]
            pub fn template_attachment_order(
                &self,
            ) -> Result<Option<crate::attachment_values::TemplateAttachmentOrder>, JsValue> {
                self.0
                    .template_attachment_order()
                    .map(|v| {
                        <ck::TemplateAttachmentOrder as crate::fields::Field>::get(
                            crate::storage::NativeField::owned(v.clone()),
                        )
                    })
                    .transpose()
            }
        }
    };
}
atom_accessors!(Atom);
atom_accessors!(MoleculeAtom);

#[wasm_bindgen]
impl MoleculeAtom {
    #[wasm_bindgen(getter,js_name=explicitValence)]
    pub fn explicit_valence(&self) -> Result<i32, JsValue> {
        self.0
            .explicit_valence()
            .map_err(crate::valence_errors::ValenceError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=implicitHydrogens)]
    pub fn implicit_hydrogens(&self) -> Result<i32, JsValue> {
        self.0
            .implicit_hydrogens()
            .map_err(crate::valence_errors::ValenceError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=totalHydrogens)]
    pub fn total_hydrogens(&self) -> Result<i32, JsValue> {
        self.0
            .total_hydrogens()
            .map_err(crate::valence_errors::ValenceError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=totalValence)]
    pub fn total_valence(&self) -> Result<i32, JsValue> {
        self.0
            .total_valence()
            .map_err(crate::valence_errors::ValenceError::from_rust)
    }
    pub fn value(&self) -> Atom {
        Atom(self.0.value().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn degree(&self) -> usize {
        self.0.degree()
    }
}

#[wasm_bindgen]
#[derive(Clone)]
pub struct Bond(pub(crate) ck::Bond);
crate::fields::value_field!(ck::Bond, Bond, BondInput, "Bond");
crate::fields::value_array!(ck::Bond, Bonds, "ReadonlyArray<Bond>");

#[wasm_bindgen]
impl Bond {
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> usize {
        self.0.id().index()
    }
    #[wasm_bindgen(getter)]
    pub fn begin(&self) -> usize {
        self.0.begin().index()
    }
    #[wasm_bindgen(getter)]
    pub fn end(&self) -> usize {
        self.0.end().index()
    }
    #[wasm_bindgen(getter)]
    pub fn order(&self) -> BondOrder {
        self.0.order().into()
    }
    #[wasm_bindgen(getter, js_name = orderCode)]
    pub fn order_code(&self) -> i64 {
        self.0.order_code()
    }
    #[wasm_bindgen(getter, js_name = orderName)]
    pub fn order_name(&self) -> String {
        self.0.order_name().to_owned()
    }
    #[wasm_bindgen(getter)]
    pub fn direction(&self) -> BondDirection {
        self.0.direction().into()
    }
    #[wasm_bindgen(getter, js_name = directionCode)]
    pub fn direction_code(&self) -> i64 {
        self.0.direction_code()
    }
    #[wasm_bindgen(getter, js_name = directionName)]
    pub fn direction_name(&self) -> String {
        self.0.direction_name().to_owned()
    }
    #[wasm_bindgen(getter)]
    pub fn stereo(&self) -> BondStereo {
        self.0.stereo().into()
    }
    #[wasm_bindgen(getter, js_name = stereoCode)]
    pub fn stereo_code(&self) -> i64 {
        self.0.stereo_code()
    }
    #[wasm_bindgen(getter, js_name = stereoName)]
    pub fn stereo_name(&self) -> String {
        self.0.stereo_name().to_owned()
    }
    #[wasm_bindgen(getter,js_name=stereoAtoms)]
    pub fn stereo_atoms(&self) -> Option<StereoAtomPair> {
        self.0.stereo_atoms().map(|[a, b]| {
            js_sys::Object::freeze(&js_sys::Array::of2(
                &JsValue::from(a.index() as u32),
                &JsValue::from(b.index() as u32),
            ))
            .unchecked_into()
        })
    }
    #[wasm_bindgen(getter,js_name=cipDescriptor)]
    pub fn cip_descriptor(&self) -> Result<Option<CipDescriptor>, JsValue> {
        self.0
            .cip_descriptor()
            .map(|v| v.map(Into::into))
            .map_err(crate::model_errors::CipDescriptorError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=cipNeighborOrder)]
    pub fn cip_neighbor_order(&self) -> Result<Option<Vec<u32>>, JsValue> {
        self.0
            .cip_neighbor_order()
            .map_err(crate::model_errors::CipDescriptorError::from_rust)
    }
    #[wasm_bindgen(getter, js_name = isAromatic)]
    pub fn is_aromatic(&self) -> bool {
        self.0.is_aromatic()
    }
}
