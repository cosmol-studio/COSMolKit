//! Native valence diagnostics; payload message is distinct from Error.message.
use crate::{
    enums::{BondOrder, ValencePhase},
    errors::{enum_error, exception},
    model_errors::TopologyValidationError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    #[derive(Clone)]
    pub type ValenceCause;
}
pub(crate) fn cause(value: &ck::ValenceError) -> Option<ValenceCause> {
    match value {
        ck::ValenceError::InvalidTopology { source } => {
            Some(TopologyValidationError::from_rust(source.clone()).unchecked_into())
        }
        ck::ValenceError::PiElectronExplicitValenceCacheNotInitialized { .. }
        | ck::ValenceError::PiElectronInvariant { .. }
        | ck::ValenceError::InvalidValence { .. }
        | ck::ValenceError::AtomOutOfRange { .. }
        | ck::ValenceError::AdjacencyAtomOutOfRange { .. }
        | ck::ValenceError::AdjacencyBondOutOfRange { .. }
        | ck::ValenceError::AdjacencyEndpointMismatch { .. }
        | ck::ValenceError::InvalidExplicitValenceInput { .. }
        | ck::ValenceError::PeriodicTableLookup { .. }
        | ck::ValenceError::ExplicitValenceCacheNotInitialized { .. }
        | ck::ValenceError::ImplicitValenceCacheNotInitialized { .. }
        | ck::ValenceError::HydrogenCountOverflow { .. }
        | ck::ValenceError::BadBondType { .. } => None,
    }
}
enum_error!(ValenceError, ValenceErrorKind, ck::ValenceError, {
    PiElectronExplicitValenceCacheNotInitialized => ck::ValenceError::PiElectronExplicitValenceCacheNotInitialized { .. },
    PiElectronInvariant => ck::ValenceError::PiElectronInvariant { .. }, InvalidValence => ck::ValenceError::InvalidValence { .. },
    InvalidTopology => ck::ValenceError::InvalidTopology { .. }, AtomOutOfRange => ck::ValenceError::AtomOutOfRange { .. },
    AdjacencyAtomOutOfRange => ck::ValenceError::AdjacencyAtomOutOfRange { .. }, AdjacencyBondOutOfRange => ck::ValenceError::AdjacencyBondOutOfRange { .. },
    AdjacencyEndpointMismatch => ck::ValenceError::AdjacencyEndpointMismatch { .. }, InvalidExplicitValenceInput => ck::ValenceError::InvalidExplicitValenceInput { .. },
    PeriodicTableLookup => ck::ValenceError::PeriodicTableLookup { .. }, ExplicitValenceCacheNotInitialized => ck::ValenceError::ExplicitValenceCacheNotInitialized { .. },
    ImplicitValenceCacheNotInitialized => ck::ValenceError::ImplicitValenceCacheNotInitialized { .. }, HydrogenCountOverflow => ck::ValenceError::HydrogenCountOverflow { .. },
    BadBondType => ck::ValenceError::BadBondType { .. },
}, {
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner {
        ck::ValenceError::PiElectronExplicitValenceCacheNotInitialized { atom } | ck::ValenceError::PiElectronInvariant { atom, .. } | ck::ValenceError::InvalidValence { atom, .. }
        | ck::ValenceError::AtomOutOfRange { atom, .. } | ck::ValenceError::AdjacencyAtomOutOfRange { atom, .. } | ck::ValenceError::AdjacencyBondOutOfRange { atom, .. }
        | ck::ValenceError::AdjacencyEndpointMismatch { atom, .. } | ck::ValenceError::InvalidExplicitValenceInput { atom, .. } | ck::ValenceError::ExplicitValenceCacheNotInitialized { atom }
        | ck::ValenceError::ImplicitValenceCacheNotInitialized { atom } | ck::ValenceError::HydrogenCountOverflow { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = explicitValence)] pub fn explicit_valence(&self) -> Option<u32> { match *self.inner { ck::ValenceError::PiElectronInvariant { explicit_valence, .. } => Some(explicit_valence), _ => None } }
    #[wasm_bindgen(getter, js_name = physicalBonds)] pub fn physical_bonds(&self) -> Option<u32> { match *self.inner { ck::ValenceError::PiElectronInvariant { physical_bonds, .. } => Some(physical_bonds), _ => None } }
    #[wasm_bindgen(getter, js_name = atomicNumber)] pub fn atomic_number(&self) -> Option<u8> { match *self.inner { ck::ValenceError::InvalidValence { atomic_number, .. } | ck::ValenceError::PeriodicTableLookup { atomic_number, .. } => Some(atomic_number), _ => None } }
    #[wasm_bindgen(getter, js_name = formalCharge)] pub fn formal_charge(&self) -> Option<i8> { match *self.inner { ck::ValenceError::InvalidValence { formal_charge, .. } => Some(formal_charge), _ => None } }
    #[wasm_bindgen(getter)] pub fn phase(&self) -> Option<ValencePhase> { match *self.inner { ck::ValenceError::InvalidValence { phase, .. } => Some(phase.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn calculated(&self) -> Option<i32> { match *self.inner { ck::ValenceError::InvalidValence { calculated, .. } => calculated, _ => None } }
    #[wasm_bindgen(getter)] pub fn reason(&self) -> Option<String> { match *self.inner { ck::ValenceError::InvalidValence { reason, .. } => Some(reason.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn message(&self) -> Option<String> { match self.inner.as_ref() { ck::ValenceError::InvalidValence { message, .. } => Some(message.clone()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::ValenceError::AtomOutOfRange { atom_count, .. } | ck::ValenceError::AdjacencyAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = neighborAtom)] pub fn neighbor_atom(&self) -> Option<usize> { match *self.inner { ck::ValenceError::AdjacencyAtomOutOfRange { neighbor_atom, .. } | ck::ValenceError::AdjacencyEndpointMismatch { neighbor_atom, .. } => Some(neighbor_atom), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner { ck::ValenceError::AdjacencyBondOutOfRange { bond, .. } | ck::ValenceError::AdjacencyEndpointMismatch { bond, .. } => Some(bond.index()), ck::ValenceError::BadBondType { bond, .. } => bond.map(ck::BondId::index), _ => None } }
    #[wasm_bindgen(getter, js_name = bondCount)] pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::ValenceError::AdjacencyBondOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn begin(&self) -> Option<usize> { match *self.inner { ck::ValenceError::AdjacencyEndpointMismatch { begin, .. } => Some(begin.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn end(&self) -> Option<usize> { match *self.inner { ck::ValenceError::AdjacencyEndpointMismatch { end, .. } => Some(end.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn value(&self) -> Option<i32> { match *self.inner { ck::ValenceError::InvalidExplicitValenceInput { value, .. } => Some(value), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::ValenceError::PeriodicTableLookup { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn explicit(&self) -> Option<u32> { match *self.inner { ck::ValenceError::HydrogenCountOverflow { explicit, .. } => Some(explicit), _ => None } }
    #[wasm_bindgen(getter)] pub fn implicit(&self) -> Option<u32> { match *self.inner { ck::ValenceError::HydrogenCountOverflow { implicit, .. } => Some(implicit), _ => None } }
    #[wasm_bindgen(getter, js_name = neighborHydrogens)] pub fn neighbor_hydrogens(&self) -> Option<usize> { match *self.inner { ck::ValenceError::HydrogenCountOverflow { neighbor_hydrogens, .. } => Some(neighbor_hydrogens), _ => None } }
    #[wasm_bindgen(getter)] pub fn order(&self) -> Option<BondOrder> { match *self.inner { ck::ValenceError::BadBondType { order, .. } => Some(order.into()), _ => None } }
}, source = ValenceCause, cause = cause);
