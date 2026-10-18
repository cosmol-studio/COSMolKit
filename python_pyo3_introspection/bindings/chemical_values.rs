use cosmolkit as ck;
use pyo3::prelude::*;

pub(crate) trait EnumValue: Sized {
    fn from_name(name: &str) -> PyResult<Self>;
}

pub(crate) fn spelling(name: &str) -> String {
    let chars: Vec<_> = name.chars().collect();
    let mut result = String::new();
    for (index, ch) in chars.iter().copied().enumerate() {
        if index > 0 {
            let previous = chars[index - 1];
            let boundary = (ch.is_ascii_uppercase()
                && (previous.is_ascii_lowercase()
                    || previous.is_ascii_uppercase()
                        && chars.get(index + 1).is_some_and(char::is_ascii_lowercase)))
                || ch.is_ascii_digit() && previous.is_ascii_lowercase();
            if boundary {
                result.push('_');
            }
        }
        result.push(ch.to_ascii_lowercase());
    }
    result
}

// Shared binding representation only. Every variant conversion remains
// exhaustive against the canonical Rust type; metadata may escape Python keywords.
macro_rules! enum_value {
    ($name:ident, {$($(#[$attribute:meta])* $variant:ident),+ $(,)?}) => {
        enum_value!(@mapped $name, { $($(#[$attribute])* $variant => $variant),+ });
    };
    ($name:ident, {$($variant:ident => $source:ident),+ $(,)?}) => {
        enum_value!(@mapped $name, { $($variant => $source),+ });
    };
    (@mapped $name:ident, {$($(#[$attribute:meta])* $variant:ident => $source:ident),+ $(,)?}) => {
        #[pyo3::pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
        #[derive(Clone, Copy, PartialEq)]
        pub enum $name { $($(#[$attribute])* $variant = const {
            let value = ck::$name::$source as i128;
            assert!(value >= isize::MIN as i128 && value <= isize::MAX as i128);
            value as isize
        }),+ }

        impl From<ck::$name> for $name {
            fn from(value: ck::$name) -> Self {
                match value { $(ck::$name::$source => Self::$variant),+ }
            }
        }
        impl From<$name> for ck::$name {
            fn from(value: $name) -> Self {
                match value { $($name::$variant => Self::$source),+ }
            }
        }
        impl crate::chemical_values::EnumValue for $name {
            fn from_name(name: &str) -> pyo3::PyResult<Self> {
                static VALUES: std::sync::OnceLock<Vec<(String, $name)>> = std::sync::OnceLock::new();
                let values = VALUES.get_or_init(|| vec![$((crate::chemical_values::spelling(stringify!($variant)), Self::$variant)),+]);
                values.iter().find(|(text, _)| text == name).map(|(_, value)| *value)
                    .ok_or_else(|| pyo3::exceptions::PyValueError::new_err(format!("unknown {}: {name}", stringify!($name))))
            }
        }
    };
}
pub(crate) use enum_value;

enum_value!(BondOrder, {
    Unspecified, Single, Double, Triple, Quadruple, Quintuple, Hextuple,
    OneAndHalf, TwoAndHalf, ThreeAndHalf, FourAndHalf, FiveAndHalf, Aromatic,
    Ionic, Hydrogen, ThreeCenter, DativeOne, Dative, DativeLeft, DativeRight,
    Other, Zero,
});
enum_value!(ChiralTag, {
    Unspecified, TetrahedralCw, TetrahedralCcw, Other, Tetrahedral, Allene,
    SquarePlanar, TrigonalBipyramidal, Octahedral,
});
enum_value!(Hybridization, { Unspecified, S, Sp, Sp2, Sp3, Sp2d, Sp3d, Sp3d2, Other });
enum_value!(BondDirection, { #[pyo3(name = "None_")] None, BeginWedge, BeginDash, EndDownRight, EndUpRight, EitherDouble, Unknown });
enum_value!(BondStereo, { #[pyo3(name = "None_")] None, Any, Z, E, Cis, Trans, AtropCw, AtropCcw });
enum_value!(McsAtomComparator, {
    Any => AtomCompareAny,
    Elements => AtomCompareElements,
    Isotopes => AtomCompareIsotopes,
    AnyHeavyAtom => AtomCompareAnyHeavyAtom,
});
enum_value!(McsBondComparator, {
    Any => BondCompareAny,
    Order => BondCompareOrder,
    OrderExact => BondCompareOrderExact,
});

#[derive(FromPyObject, IntoPyObject)]
pub enum BondOrderInput {
    Enum(BondOrder),
    Name(String),
}
impl BondOrderInput {
    pub(crate) fn into_rust(self) -> PyResult<ck::BondOrder> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(name) => BondOrder::from_name(&name).map(Into::into),
        }
    }
}
