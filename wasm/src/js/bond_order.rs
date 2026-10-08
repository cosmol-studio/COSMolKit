//! Canonical bond-order vocabulary, shared by all projections.
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum BondOrder {
    Unspecified = 0,
    Single = 1,
    Double = 2,
    Triple = 3,
    Quadruple = 4,
    Quintuple = 5,
    Hextuple = 6,
    OneAndHalf = 7,
    TwoAndHalf = 8,
    ThreeAndHalf = 9,
    FourAndHalf = 10,
    FiveAndHalf = 11,
    Aromatic = 12,
    Ionic = 13,
    Hydrogen = 14,
    ThreeCenter = 15,
    DativeOne = 16,
    Dative = 17,
    DativeLeft = 18,
    DativeRight = 19,
    Other = 20,
    Zero = 21,
}
