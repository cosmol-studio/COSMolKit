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
    ($name:ident, input = $input:ty, {$($variant:ident),+ $(,)?}) => {
        enum_value!(@mapped $name, input = $input, { $($variant => $variant),+ });
    };
    ($name:ident, {$($(#[$($attribute:tt)*])* $variant:ident $([$($alias:ident),*])? $(text = $text:literal)?),+ $(,)?}
        $(, constructor = $constructor:ident)? $(, methods = {$($method:item)*})?) => {
        enum_value!(@mapped $name, { $($(#[$($attribute)*])* $variant => $variant $([$($alias),*])? $(text = $text)?),+ }
            $(, constructor = $constructor)? $(, methods = {$($method)*})?);
    };
    ($name:ident, {$($variant:ident => $source:ident $([$($alias:ident),*])? $(text = $text:literal)?),+ $(,)?}
        $(, constructor = $constructor:ident)? $(, methods = {$($method:item)*})?) => {
        enum_value!(@mapped $name, { $($variant => $source $([$($alias),*])? $(text = $text)?),+ }
            $(, constructor = $constructor)? $(, methods = {$($method)*})?);
    };
    (@text $variant:ident) => { crate::chemical_values::spelling(stringify!($variant)) };
    (@text $variant:ident, $text:literal) => { String::from($text) };
    (@mapped $name:ident, $(input = $input:ty,)? {$($(#[$($attribute:tt)*])* $variant:ident => $source:ident $([$($alias:ident),*])? $(text = $text:literal)?),+ $(,)?}
        $(, constructor = $constructor:ident)? $(, methods = {$($method:item)*})?) => {
        #[pyo3::pyclass(module = "cosmolkit", eq, eq_int, from_py_object)]
        #[derive(Clone, Copy, PartialEq)]
        pub enum $name { $($(#[$($attribute)*])* $variant = const {
            let value = ck::$name::$source as i128;
            assert!(value >= isize::MIN as i128 && value <= isize::MAX as i128);
            value as isize
        }),+ }

        $crate::enums::enum_protocol!(@select [$($constructor)?] $name,
            {$($variant => ($crate::enums::public_name!($variant; $(#[$($attribute)*])*), [$($($alias),*)?])),+},
            methods = {$($($method)*)?});

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
        impl crate::parameters::NativeEnum for $name {
            type Native = ck::$name;
            fn from_native(value: Self::Native) -> Self { value.into() }
            fn into_native(self) -> Self::Native { self.into() }
            $(
                const INPUT: pyo3::inspect::PyStaticExpr = <$input as pyo3::FromPyObject<'_, '_>>::INPUT_TYPE;
                fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
                    value.extract::<$input>()?.into_rust()
                }
            )?
        }
        impl crate::configuration::FieldProjection for ck::$name {
            type Adapter = crate::parameters::EnumField<$name>;
        }
        impl crate::chemical_values::EnumValue for $name {
            fn from_name(name: &str) -> pyo3::PyResult<Self> {
                static VALUES: std::sync::OnceLock<Vec<(String, $name)>> = std::sync::OnceLock::new();
                let values = VALUES.get_or_init(|| vec![$((enum_value!(@text $variant $(, $text)?), Self::$variant)),+]);
                values.iter().find(|(text, _)| text == name).map(|(_, value)| *value)
                    .ok_or_else(|| pyo3::exceptions::PyValueError::new_err(format!("unknown {}: {name}", stringify!($name))))
            }
        }
    };
}
pub(crate) use enum_value;

enum_value!(BondOrder, {
    Unspecified [UNSPECIFIED], Single [SINGLE], Double [DOUBLE], Triple [TRIPLE], Quadruple [QUADRUPLE], Quintuple [QUINTUPLE], Hextuple [HEXTUPLE],
    OneAndHalf [ONEANDAHALF], TwoAndHalf [TWOANDAHALF], ThreeAndHalf [THREEANDAHALF], FourAndHalf [FOURANDAHALF], FiveAndHalf [FIVEANDAHALF], Aromatic [AROMATIC],
    Ionic [IONIC], Hydrogen [HYDROGEN], ThreeCenter [THREECENTER], DativeOne [DATIVEONE], Dative [DATIVE], DativeLeft [DATIVEL], DativeRight [DATIVER],
    Other [OTHER], Zero [ZERO],
});
enum_value!(ChiralTag, {
    Unspecified [CHI_UNSPECIFIED], TetrahedralCw [CHI_TETRAHEDRAL_CW], TetrahedralCcw [CHI_TETRAHEDRAL_CCW], Other [CHI_OTHER], Tetrahedral [CHI_TETRAHEDRAL], Allene [CHI_ALLENE],
    SquarePlanar [CHI_SQUAREPLANAR], TrigonalBipyramidal [CHI_TRIGONALBIPYRAMIDAL], Octahedral [CHI_OCTAHEDRAL],
});
enum_value!(Hybridization, { Unspecified [UNSPECIFIED], S, Sp [SP], Sp2 [SP2], Sp3 [SP3], Sp2d [SP2D], Sp3d [SP3D], Sp3d2 [SP3D2], Other [OTHER] });

enum_value!(BondDirection, { #[pyo3(name = "None_")] None [NONE], BeginWedge [BEGINWEDGE], BeginDash [BEGINDASH], EndDownRight [ENDDOWNRIGHT], EndUpRight [ENDUPRIGHT], EitherDouble [EITHERDOUBLE], Unknown [UNKNOWN] });
enum_value!(BondStereo, { #[pyo3(name = "None_")] None [STEREONONE], Any [STEREOANY], Z [STEREOZ], E [STEREOE], Cis [STEREOCIS], Trans [STEREOTRANS], AtropCw [STEREOATROPCW], AtropCcw [STEREOATROPCCW] });
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
