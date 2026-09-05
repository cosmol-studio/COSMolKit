use std::{fmt, str::FromStr};

use serde::{Deserialize, Deserializer, Serialize, Serializer};
use thiserror::Error;

/// Stable chemical-element identity.
///
/// The private atomic number is always in the inclusive range `0..=118`.
/// Zero is the source-compatible dummy atom (`*`); `1..=118` are H through
/// Og. Use [`Element::from_atomic_number`] or [`FromStr`] for checked external
/// construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Element {
    atomic_number: u8,
}

/// Error returned when a string is not a source-recognized element symbol.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[error("unknown element symbol '{input}'")]
pub struct ElementParseError {
    input: String,
}

impl ElementParseError {
    #[must_use]
    pub fn input(&self) -> &str {
        &self.input
    }
}

/// Source-aligned periodic-table metadata for an [`Element`].
///
/// This is a read-only projection of the periodic-table data used by the
/// chemistry core. `valences` preserves the source sentinel values, including
/// `-1`; callers must not reinterpret those values as an exhaustive IUPAC
/// oxidation-state table.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
#[non_exhaustive]
pub struct ElementInfo {
    pub element: Element,
    pub symbol: &'static str,
    pub atomic_number: u8,
    pub period: u8,
    pub outer_electrons: i32,
    pub valences: &'static [i32],
    /// RDKit's source-defined `Rb0` bond radius in angstroms.
    pub rb0: f64,
    pub atomic_weight: f64,
}

impl Element {
    pub const DUMMY: Self = Self { atomic_number: 0 };
    // Period 1
    pub const H: Self = Self { atomic_number: 1 };
    pub const HE: Self = Self { atomic_number: 2 };
    // Period 2
    pub const LI: Self = Self { atomic_number: 3 };
    pub const BE: Self = Self { atomic_number: 4 };
    pub const B: Self = Self { atomic_number: 5 };
    pub const C: Self = Self { atomic_number: 6 };
    pub const N: Self = Self { atomic_number: 7 };
    pub const O: Self = Self { atomic_number: 8 };
    pub const F: Self = Self { atomic_number: 9 };
    pub const NE: Self = Self { atomic_number: 10 };
    // Period 3
    pub const NA: Self = Self { atomic_number: 11 };
    pub const MG: Self = Self { atomic_number: 12 };
    pub const AL: Self = Self { atomic_number: 13 };
    pub const SI: Self = Self { atomic_number: 14 };
    pub const P: Self = Self { atomic_number: 15 };
    pub const S: Self = Self { atomic_number: 16 };
    pub const CL: Self = Self { atomic_number: 17 };
    pub const AR: Self = Self { atomic_number: 18 };
    // Period 4
    pub const K: Self = Self { atomic_number: 19 };
    pub const CA: Self = Self { atomic_number: 20 };
    pub const SC: Self = Self { atomic_number: 21 };
    pub const TI: Self = Self { atomic_number: 22 };
    pub const V: Self = Self { atomic_number: 23 };
    pub const CR: Self = Self { atomic_number: 24 };
    pub const MN: Self = Self { atomic_number: 25 };
    pub const FE: Self = Self { atomic_number: 26 };
    pub const CO: Self = Self { atomic_number: 27 };
    pub const NI: Self = Self { atomic_number: 28 };
    pub const CU: Self = Self { atomic_number: 29 };
    pub const ZN: Self = Self { atomic_number: 30 };
    pub const GA: Self = Self { atomic_number: 31 };
    pub const GE: Self = Self { atomic_number: 32 };
    pub const AS: Self = Self { atomic_number: 33 };
    pub const SE: Self = Self { atomic_number: 34 };
    pub const BR: Self = Self { atomic_number: 35 };
    pub const KR: Self = Self { atomic_number: 36 };
    // Period 5
    pub const RB: Self = Self { atomic_number: 37 };
    pub const SR: Self = Self { atomic_number: 38 };
    pub const Y: Self = Self { atomic_number: 39 };
    pub const ZR: Self = Self { atomic_number: 40 };
    pub const NB: Self = Self { atomic_number: 41 };
    pub const MO: Self = Self { atomic_number: 42 };
    pub const TC: Self = Self { atomic_number: 43 };
    pub const RU: Self = Self { atomic_number: 44 };
    pub const RH: Self = Self { atomic_number: 45 };
    pub const PD: Self = Self { atomic_number: 46 };
    pub const AG: Self = Self { atomic_number: 47 };
    pub const CD: Self = Self { atomic_number: 48 };
    pub const IN: Self = Self { atomic_number: 49 };
    pub const SN: Self = Self { atomic_number: 50 };
    pub const SB: Self = Self { atomic_number: 51 };
    pub const TE: Self = Self { atomic_number: 52 };
    pub const I: Self = Self { atomic_number: 53 };
    pub const XE: Self = Self { atomic_number: 54 };
    // Period 6
    pub const CS: Self = Self { atomic_number: 55 };
    pub const BA: Self = Self { atomic_number: 56 };
    // Lanthanides
    pub const LA: Self = Self { atomic_number: 57 };
    pub const CE: Self = Self { atomic_number: 58 };
    pub const PR: Self = Self { atomic_number: 59 };
    pub const ND: Self = Self { atomic_number: 60 };
    pub const PM: Self = Self { atomic_number: 61 };
    pub const SM: Self = Self { atomic_number: 62 };
    pub const EU: Self = Self { atomic_number: 63 };
    pub const GD: Self = Self { atomic_number: 64 };
    pub const TB: Self = Self { atomic_number: 65 };
    pub const DY: Self = Self { atomic_number: 66 };
    pub const HO: Self = Self { atomic_number: 67 };
    pub const ER: Self = Self { atomic_number: 68 };
    pub const TM: Self = Self { atomic_number: 69 };
    pub const YB: Self = Self { atomic_number: 70 };
    pub const LU: Self = Self { atomic_number: 71 };
    // Period 6 continued
    pub const HF: Self = Self { atomic_number: 72 };
    pub const TA: Self = Self { atomic_number: 73 };
    pub const W: Self = Self { atomic_number: 74 };
    pub const RE: Self = Self { atomic_number: 75 };
    pub const OS: Self = Self { atomic_number: 76 };
    pub const IR: Self = Self { atomic_number: 77 };
    pub const PT: Self = Self { atomic_number: 78 };
    pub const AU: Self = Self { atomic_number: 79 };
    pub const HG: Self = Self { atomic_number: 80 };
    pub const TL: Self = Self { atomic_number: 81 };
    pub const PB: Self = Self { atomic_number: 82 };
    pub const BI: Self = Self { atomic_number: 83 };
    pub const PO: Self = Self { atomic_number: 84 };
    pub const AT: Self = Self { atomic_number: 85 };
    pub const RN: Self = Self { atomic_number: 86 };
    // Period 7
    pub const FR: Self = Self { atomic_number: 87 };
    pub const RA: Self = Self { atomic_number: 88 };
    // Actinides
    pub const AC: Self = Self { atomic_number: 89 };
    pub const TH: Self = Self { atomic_number: 90 };
    pub const PA: Self = Self { atomic_number: 91 };
    pub const U: Self = Self { atomic_number: 92 };
    pub const NP: Self = Self { atomic_number: 93 };
    pub const PU: Self = Self { atomic_number: 94 };
    pub const AM: Self = Self { atomic_number: 95 };
    pub const CM: Self = Self { atomic_number: 96 };
    pub const BK: Self = Self { atomic_number: 97 };
    pub const CF: Self = Self { atomic_number: 98 };
    pub const ES: Self = Self { atomic_number: 99 };
    pub const FM: Self = Self { atomic_number: 100 };
    pub const MD: Self = Self { atomic_number: 101 };
    pub const NO: Self = Self { atomic_number: 102 };
    pub const LR: Self = Self { atomic_number: 103 };
    // Period 7 continued
    pub const RF: Self = Self { atomic_number: 104 };
    pub const DB: Self = Self { atomic_number: 105 };
    pub const SG: Self = Self { atomic_number: 106 };
    pub const BH: Self = Self { atomic_number: 107 };
    pub const HS: Self = Self { atomic_number: 108 };
    pub const MT: Self = Self { atomic_number: 109 };
    pub const DS: Self = Self { atomic_number: 110 };
    pub const RG: Self = Self { atomic_number: 111 };
    pub const CN: Self = Self { atomic_number: 112 };
    pub const NH: Self = Self { atomic_number: 113 };
    pub const FL: Self = Self { atomic_number: 114 };
    pub const MC: Self = Self { atomic_number: 115 };
    pub const LV: Self = Self { atomic_number: 116 };
    pub const TS: Self = Self { atomic_number: 117 };
    pub const OG: Self = Self { atomic_number: 118 };

    /// Construct an element from a checked atomic number.
    ///
    /// `0` is the dummy atom and `1..=118` are the real elements. Values above
    /// 118 are rejected rather than creating an invalid element identity.
    #[must_use]
    pub const fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        if atomic_number <= 118 {
            Some(Self { atomic_number })
        } else {
            None
        }
    }

    #[must_use]
    pub const fn atomic_number(self) -> u8 {
        self.atomic_number
    }

    /// Construct an element from a source-recognized symbol.
    ///
    /// Parsing is case-sensitive, matching the source periodic-table API.
    /// Legacy source aliases such as `Uut` and `Uup` are accepted and
    /// canonicalize to `Nh` and `Mc` respectively.
    #[must_use]
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        crate::periodic_table::atomic_number_from_symbol(symbol).and_then(Self::from_atomic_number)
    }

    /// Return the canonical element symbol (`*`, `H` through `Og`).
    #[must_use]
    pub fn symbol(self) -> &'static str {
        crate::periodic_table::symbol(self.atomic_number).expect("Element's private atomic number is always in 0..=118")
    }

    /// Return the source-aligned periodic-table metadata used by COSMolKit.
    #[must_use]
    pub fn info(self) -> ElementInfo {
        let atomic_number = self.atomic_number;
        ElementInfo {
            element: self,
            symbol: self.symbol(),
            atomic_number,
            period: crate::periodic_table::period(atomic_number)
                .expect("Element's private atomic number is always in 0..=118"),
            outer_electrons: crate::periodic_table::outer_electrons(atomic_number)
                .expect("Element's private atomic number is always in 0..=118"),
            valences: crate::periodic_table::valences(atomic_number)
                .expect("Element's private atomic number is always in 0..=118"),
            rb0: crate::periodic_table::rb0(atomic_number),
            atomic_weight: crate::periodic_table::atomic_weight(atomic_number),
        }
    }

    /// Iterate over the 118 real elements, excluding the dummy atom.
    pub fn iter() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator {
        ELEMENTS.iter().copied()
    }

    /// Iterate over the dummy atom followed by all 118 real elements.
    pub fn iter_with_dummy() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator {
        ELEMENTS_WITH_DUMMY.iter().copied()
    }
}

impl fmt::Display for Element {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.symbol())
    }
}

impl FromStr for Element {
    type Err = ElementParseError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        Self::from_symbol(input).ok_or_else(|| ElementParseError {
            input: input.to_string(),
        })
    }
}

impl Serialize for Element {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(self.symbol())
    }
}

impl<'de> Deserialize<'de> for Element {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let symbol = String::deserialize(deserializer)?;
        symbol.parse().map_err(serde::de::Error::custom)
    }
}

const fn element_array<const N: usize>(first_atomic_number: u8) -> [Element; N] {
    let mut elements = [Element::DUMMY; N];
    let mut index = 0;
    while index < N {
        elements[index] = Element {
            atomic_number: first_atomic_number + index as u8,
        };
        index += 1;
    }
    elements
}

/// All 118 real elements in ascending atomic-number order (H through Og).
pub static ELEMENTS: [Element; 118] = element_array::<118>(1);

/// The dummy atom followed by all real elements in atomic-number order.
pub static ELEMENTS_WITH_DUMMY: [Element; 119] = element_array::<119>(0);
