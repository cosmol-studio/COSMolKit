//! Public, dependency-light chemical vocabulary shared by COSMolKit crates.
//!
//! This crate intentionally has no dependency on a molecule owner, parser,
//! operation runtime, or chemistry algorithm implementation.

mod vocabulary;

pub use vocabulary::{
    BondDirection, BondOrder, BondStereo, ChiralTag, ELEMENTS, ELEMENTS_WITH_DUMMY, Element,
    ElementInfo, ElementParseError, Hybridization,
};
