//! One declaration owns each enum's native conversion and JavaScript vocabulary.
use cosmolkit as ck;
use wasm_bindgen::prelude::*;

macro_rules! enum_value {
    ($name:ident, $ts:literal { $($variant:ident $(=> $number:literal)? = $text:literal),+ $(,)? }) => {
        #[wasm_bindgen]
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub enum $name {
            $($variant $(= $number)?),+
        }

        // wasm-bindgen requires literal discriminants. The declaration order
        // remains compile-time checked against the actual native enum values.
        const _: () = {
            $(assert!($name::$variant as i64 == ck::$name::$variant as i64);)+
        };

        impl From<ck::$name> for $name {
            fn from(value: ck::$name) -> Self {
                match value {
                    $(ck::$name::$variant => Self::$variant),+
                }
            }
        }
        impl TryFrom<$name> for ck::$name {
            type Error = JsValue;
            fn try_from(value: $name) -> Result<Self, JsValue> {
                match value {
                    $($name::$variant => Ok(Self::$variant),)+
                }
            }
        }
        paste::paste! {
            #[wasm_bindgen]
            extern "C" {
                #[wasm_bindgen(typescript_type = $ts)]
                pub type [<$name Input>];
            }
            impl crate::fields::Field for ck::$name {
                type Input = [<$name Input>];
                type Output = $name;
                fn input(value: Self::Input) -> Result<Self, JsValue> {
                    Self::decode(value.as_ref())
                }
                fn decode(raw: &JsValue) -> Result<Self, JsValue> {
                    let value = crate::input::read::<u32>(tsify::Ts::new_unchecked(raw.clone()))?;
                    match value {
                        $(value if value == $name::$variant as u32 => Ok(Self::$variant),)+
                        _ => Err(
                            js_sys::TypeError::new(concat!("unknown ", stringify!($name), " value")).into()
                        ),
                    }
                }
                fn get(value: crate::storage::NativeField<Self>) -> Result<$name, JsValue> {
                    value.read(|value| (*value).into())
                }
            }
        }
    };
}

pub(crate) use enum_value;

#[cfg(feature = "cap-transforms")]
#[wasm_bindgen(js_name = coordinateZPolicyFromName)]
pub fn coordinate_z_policy_from_name(
    text: tsify::Ts<crate::input::TextInput>,
) -> Result<CoordinateZPolicy, JsValue> {
    let text = crate::input::read::<String>(text)?;
    ck::CoordinateZPolicy::from_name(&text)
        .map(Into::into)
        .map_err(crate::transform_errors::CoordinateInputError::from_rust)
}

enum_value!(BondOrder, "BondOrder" {
    Unspecified = "unspecified", Single = "single", Double = "double", Triple = "triple",
    Quadruple = "quadruple", Quintuple = "quintuple", Hextuple = "hextuple",
    OneAndHalf = "one_and_half", TwoAndHalf = "two_and_half", ThreeAndHalf = "three_and_half",
    FourAndHalf = "four_and_half", FiveAndHalf = "five_and_half", Aromatic = "aromatic",
    Ionic = "ionic", Hydrogen = "hydrogen", ThreeCenter = "three_center",
    DativeOne = "dative_one", Dative = "dative", DativeLeft = "dative_left",
    DativeRight = "dative_right", Other = "other", Zero = "zero",
});
enum_value!(ChiralTag, "ChiralTag" {
    Unspecified = "unspecified", TetrahedralCw = "tetrahedral_cw", TetrahedralCcw = "tetrahedral_ccw",
    Other = "other", Tetrahedral = "tetrahedral", Allene = "allene", SquarePlanar = "square_planar",
    TrigonalBipyramidal = "trigonal_bipyramidal", Octahedral = "octahedral",
});
enum_value!(Hybridization, "Hybridization" {
    Unspecified = "unspecified", S = "s", Sp = "sp", Sp2 = "sp_2", Sp3 = "sp_3",
    Sp2d = "sp_2d", Sp3d = "sp_3d", Sp3d2 = "sp_3d_2", Other = "other",
});
enum_value!(BondDirection, "BondDirection" {
    None = "none", BeginWedge = "begin_wedge", BeginDash = "begin_dash",
    EndDownRight = "end_down_right", EndUpRight = "end_up_right",
    EitherDouble = "either_double", Unknown = "unknown",
});
enum_value!(BondStereo, "BondStereo" {
    None = "none", Any = "any", Z = "z", E = "e", Cis = "cis", Trans = "trans",
    AtropCw = "atrop_cw", AtropCcw = "atrop_ccw",
});
enum_value!(CipDescriptor, "CipDescriptor" {
    R = "R", S = "S", LowerR = "r", LowerS = "s", E = "E", Z = "Z",
    LowerE = "e", LowerZ = "z", M = "M", P = "P", LowerM = "m", LowerP = "p",
});
enum_value!(MoleculeOpOutput, "MoleculeOpOutput" { Single = "single", Multiple = "multiple", LazyMultiple = "lazy_multiple" });
enum_value!(ParityPolicy, "ParityPolicy" { NotApplicable = "not_applicable", RequiredWhenSupported = "required_when_supported", RequiredNow = "required_now" });
enum_value!(ValencePhase, "ValencePhase" { EffectiveAtomicNumber = "effective_atomic_number", Explicit = "explicit", Implicit = "implicit" });
enum_value!(SdfPropertyListTarget, "SdfPropertyListTarget" { Atom = "atom", Bond = "bond" });
#[cfg(feature = "cap-io")]
enum_value!(SdfCoordinateMode, "SdfCoordinateMode" { Preserve = "preserve", Require2D = "require_2d", Require3D = "require_3d" });
enum_value!(CoordinateDimension, "CoordinateDimension" { TwoD = "two_d", ThreeD = "three_d" });
#[cfg(feature = "cap-io")]
enum_value!(SdfFormat, "SdfFormat" { V2000 = "v2000", V3000 = "v3000" });
#[cfg(feature = "cap-io")]
enum_value!(Mol2Type, "Mol2Type" { Corina = "corina" });
#[cfg(feature = "cap-descriptors")]
enum_value!(RotatableBondsOptions, "RotatableBondsOptions" { Default = "default", NonStrict = "non_strict", Strict = "strict", StrictLinkages = "strict_linkages" });
#[cfg(feature = "cap-transforms")]
enum_value!(CoordinateZPolicy, "CoordinateZPolicy" { Ignore = "ignore", RequireZero = "require_zero", Error = "error" });

enum_value!(SubstanceGroupAtomList, "SubstanceGroupAtomList" { Atoms = "atoms", ParentAtoms = "parent_atoms" });
#[cfg(feature = "cap-valence")]
enum_value!(ValenceModel, "ValenceModel" { RdkitLike = "rdkit_like" });
#[cfg(feature = "cap-inchi")]
enum_value!(InchiDiagnosticLevel, "InchiDiagnosticLevel" { Warning = "warning", Error = "error" });
#[cfg(feature = "cap-inchi")]
enum_value!(InchiErrorKind, "InchiErrorKind" { AllocationFailed = "allocation_failed", UnsupportedState = "unsupported_state", InvalidInput = "invalid_input", InvalidSourceOutput = "invalid_source_output", SanitizeFailed = "sanitize_failed", Toolkit = "toolkit", SourcePort = "source_port" });
enum_value!(StereoGroupKind, "StereoGroupKind" { Absolute = "absolute", Or = "or", And = "and" });
#[cfg(feature = "cap-sanitize")]
enum_value!(SanitizeStage, "SanitizeStage" {
    None => 0 = "none", Cleanup => 1 = "cleanup", Properties => 2 = "properties",
    SymmRings => 4 = "symm_rings", Kekulize => 8 = "kekulize", FindRadicals => 16 = "find_radicals",
    SetAromaticity => 32 = "set_aromaticity", SetConjugation => 64 = "set_conjugation",
    SetHybridization => 128 = "set_hybridization", CleanupChirality => 256 = "cleanup_chirality",
    AdjustHs => 512 = "adjust_hs", CleanupOrganometallics => 1024 = "cleanup_organometallics",
    CleanupAtropisomers => 2048 = "cleanup_atropisomers"
});
#[cfg(feature = "cap-aromaticity")]
enum_value!(AromaticityModel, "AromaticityModel" { Rdkit = "rdkit", Simple = "simple", Mdl = "mdl", Mmff94 = "mmff94", Custom = "custom" });
#[cfg(feature = "cap-stereo")]
enum_value!(PotentialStereoType, "PotentialStereoType" { AtomTetrahedral = "atom_tetrahedral", AtomSquarePlanar = "atom_square_planar", AtomTrigonalBipyramidal = "atom_trigonal_bipyramidal", AtomOctahedral = "atom_octahedral", BondDouble = "bond_double", BondCumuleneEven = "bond_cumulene_even", BondAtropisomer = "bond_atropisomer" });
#[cfg(feature = "cap-stereo")]
enum_value!(PotentialStereoSpecified, "PotentialStereoSpecified" { Unspecified = "unspecified", Specified = "specified", Unknown = "unknown" });
#[cfg(feature = "cap-stereo")]
enum_value!(PotentialStereoDescriptor, "PotentialStereoDescriptor" { None = "none", TetrahedralClockwise = "tetrahedral_clockwise", TetrahedralCounterclockwise = "tetrahedral_counterclockwise", BondCis = "bond_cis", BondTrans = "bond_trans", BondAtropCw = "bond_atrop_cw", BondAtropCcw = "bond_atrop_ccw" });
