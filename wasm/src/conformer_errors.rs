//! Native embedding errors, with exact source ownership and payloads.
use crate::errors::{enum_error, error_projection};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use std::error::Error;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "BoundsQueryError | GenerationError | GraphBoundsError | BoundsMatrixError | UffBoundsError | MatrixError | ValenceError | TopologyValidationError | CoordinateValidationError | RingFindingError | ConjugationError | HybridizationError | LegacyStereoError | MoleculeFragmentsError | MissingExplicitHydrogensError | PruningError | EmbedSeedError | GenerationFailure | AtropisomerError | ConformerOptimizerError | ConformerError | ThreadCountError | IoError | CrystalffTorsionPreferencesError | HydrogenError | QueryGraphError | AlignmentError | SubstructMatchError | QueryMatchContextError | UffParameterError | UffBondMathError | UffBuilderError | ForceFieldKernelError | UffInversionContributionError | StereoGroupError | TemplateAttachmentOrderError | AtomPropertyError | BondValueError | MoleculePropertyError | MappingValidationError | CanonicalRankError | PropertyStringError | RequiredPropertyStringError | StereoOrderError | DoubleBondStereoError | PotentialStereoError | TopologyEditError | PropertyIntReadError | PropertyUIntReadError | PropertyULongReadError | SmilesStereoError | ParserCarrierError | ParserStereoOrderError | KekulizeError | SanitizeError | QueryCompileError | SmartsParseError | MoleculeFragmentsFailure | PropertyValueError | UIntLexicalReadError"
    )]
    #[derive(Clone)]
    pub type ConformerCause;
}
pub(crate) fn project(owner: ReadonlyValue<dyn Error>) -> ConformerCause {
    macro_rules! selected {
        ($native:ty,$projection:path) => {
            if owner.is::<$native>() {
                return <$projection>::from_owner(
                    owner.child(|e| e.downcast_ref().expect("checked conformer source")),
                )
                .unchecked_into();
            }
        };
    }
    macro_rules! own {($($name:ident),+)=>{$(selected!(ck::$name,$name);)+};}
    own!(
        BoundsQueryError,
        GenerationError,
        GraphBoundsError,
        BoundsMatrixError,
        UffBoundsError,
        MissingExplicitHydrogensError,
        PruningError,
        EmbedSeedError,
        GenerationFailure,
        ConformerOptimizerError,
        ConformerError,
        CrystalffTorsionPreferencesError
    );
    selected!(ck::MatrixError, crate::matrix_errors::MatrixError);
    selected!(ck::ValenceError, crate::valence_errors::ValenceError);
    selected!(
        ck::TopologyValidationError,
        crate::model_errors::TopologyValidationError
    );
    selected!(
        ck::CoordinateValidationError,
        crate::model_errors::CoordinateValidationError
    );
    selected!(ck::RingFindingError, crate::ring_errors::RingFindingError);
    selected!(
        ck::ConjugationError,
        crate::sanitize_support_errors::ConjugationError
    );
    selected!(
        ck::HybridizationError,
        crate::sanitize_leaf_errors::HybridizationError
    );
    selected!(
        ck::LegacyStereoError,
        crate::legacy_stereo_errors::LegacyStereoError
    );
    selected!(
        ck::MoleculeFragmentsError,
        crate::fragment_errors::MoleculeFragmentsError
    );
    selected!(
        ck::AtropisomerError,
        crate::atropisomer_errors::AtropisomerError
    );
    selected!(ck::ThreadCountError, crate::worker_errors::ThreadCountError);
    selected!(ck::HydrogenError, crate::hydrogen_errors::HydrogenError);
    selected!(
        ck::QueryGraphError,
        crate::query_model_errors::QueryGraphError
    );
    selected!(ck::AlignmentError, crate::alignment_errors::AlignmentError);
    selected!(
        ck::SubstructMatchError,
        crate::match_errors::SubstructMatchError
    );
    selected!(
        ck::QueryMatchContextError,
        crate::match_errors::QueryMatchContextError
    );
    selected!(ck::UffParameterError, crate::uff_errors::UffParameterError);
    selected!(ck::UffBondMathError, crate::uff_errors::UffBondMathError);
    selected!(ck::UffBuilderError, crate::uff_errors::UffBuilderError);
    selected!(
        ck::ForceFieldKernelError,
        crate::forcefield_kernel_errors::ForceFieldKernelError
    );
    selected!(
        ck::UffInversionContributionError,
        crate::uff_errors::UffInversionContributionError
    );
    selected!(std::io::Error, crate::external_errors::IoError);
    selected!(
        ck::UIntLexicalReadError,
        crate::property_errors::UIntLexicalReadError
    );
    selected!(
        ck::PropertyValueError,
        crate::model_errors::PropertyValueError
    );
    selected!(ck::StereoGroupError, crate::model_errors::StereoGroupError);
    selected!(
        ck::TemplateAttachmentOrderError,
        crate::model_errors::TemplateAttachmentOrderError
    );
    selected!(
        ck::AtomPropertyError,
        crate::model_errors::AtomPropertyError
    );
    selected!(ck::BondValueError, crate::model_errors::BondValueError);
    selected!(
        ck::MoleculePropertyError,
        crate::model_errors::MoleculePropertyError
    );
    selected!(
        ck::MappingValidationError,
        crate::model_errors::MappingValidationError
    );
    selected!(
        ck::CanonicalRankError,
        crate::canonical_errors::CanonicalRankError
    );
    selected!(
        ck::PropertyStringError,
        crate::property_errors::PropertyStringError
    );
    selected!(
        ck::RequiredPropertyStringError,
        crate::property_errors::RequiredPropertyStringError
    );
    selected!(
        ck::StereoOrderError,
        crate::stereo_order_errors::StereoOrderError
    );
    selected!(
        ck::DoubleBondStereoError,
        crate::double_stereo_errors::DoubleBondStereoError
    );
    selected!(
        ck::PotentialStereoError,
        crate::potential_stereo_errors::PotentialStereoError
    );
    selected!(
        ck::TopologyEditError,
        crate::topology_errors::TopologyEditError
    );
    selected!(
        ck::PropertyIntReadError,
        crate::property_errors::PropertyIntReadError
    );
    selected!(
        ck::PropertyUIntReadError,
        crate::property_errors::PropertyUIntReadError
    );
    selected!(
        ck::PropertyULongReadError,
        crate::property_errors::PropertyULongReadError
    );
    selected!(
        ck::SmilesStereoError,
        crate::smiles_stereo_errors::SmilesStereoError
    );
    selected!(
        ck::ParserCarrierError,
        crate::parser_errors::ParserCarrierError
    );
    selected!(
        ck::ParserStereoOrderError,
        crate::parser_errors::ParserStereoOrderError
    );
    selected!(ck::KekulizeError, crate::kekulize_errors::KekulizeError);
    selected!(ck::SanitizeError, crate::sanitize_errors::SanitizeError);
    selected!(
        ck::QueryCompileError,
        crate::query_compile_error::QueryCompileError
    );
    selected!(
        ck::SmartsParseError,
        crate::smarts_parse_errors::SmartsParseError
    );
    selected!(
        ck::MoleculeFragmentsFailure,
        crate::fragment_errors::MoleculeFragmentsFailure
    );
    unreachable!("unmodeled native embedding source")
}
fn cause<T: Error + 'static>(owner: &ReadonlyValue<T>) -> Option<ConformerCause> {
    owner
        .source()
        .map(|_| project(owner.child(|e| e.source().expect("retained native source"))))
}
// Tuple variants retain their immediate payload independently of transparent
// Error::source(), which can deliberately skip that value.
macro_rules! payload {
    ($owner:expr,$native:ident,[$($variant:ident),+])=>{
        match $owner.as_ref(){
            $(ck::$native::$variant(_)=>Some(project($owner.child(|e|match e{
                ck::$native::$variant(v)=>v as &dyn Error,
                _=>unreachable!("immutable native variant"),
            }))),)+
            _=>None,
        }
    };
}
error_projection!(
    ConformerRunError,
    ck::ConformerRunError,
    {},
    source = ConformerCause,
    cause = cause::<ck::ConformerRunError>
);
error_projection!(
    MissingExplicitHydrogensError,
    ck::MissingExplicitHydrogensError,
    {},
    source = ConformerCause,
    cause = cause::<ck::MissingExplicitHydrogensError>
);
crate::enums::enum_value!(MatrixAxis,"MatrixAxis"{Row="row",Column="column"});
crate::enums::enum_value!(BoundKind,"BoundKind"{Upper="upper",Lower="lower"});
enum_error!(ConformerError,ConformerErrorKind,ck::ConformerError,{
    Unsupported=>ck::ConformerError::Unsupported,CannotNormalizeZeroLengthVector=>ck::ConformerError::CannotNormalizeZeroLengthVector,GenerationFailed=>ck::ConformerError::GenerationFailed(_),InvalidEmbedParametersJson=>ck::ConformerError::InvalidEmbedParametersJson(_)
},{#[wasm_bindgen(getter)]pub fn message(&self)->Option<String>{match self.inner.as_ref(){ck::ConformerError::GenerationFailed(v)|ck::ConformerError::InvalidEmbedParametersJson(v)=>Some(v.clone()),_=>None}}});
enum_error!(BoundsMatrixError,BoundsMatrixErrorKind,ck::BoundsMatrixError,{
    DimensionOverflow=>ck::BoundsMatrixError::DimensionOverflow{..},AllocationFailed=>ck::BoundsMatrixError::AllocationFailed{..},DataLengthTooShort=>ck::BoundsMatrixError::DataLengthTooShort{..},IndexOutOfBounds=>ck::BoundsMatrixError::IndexOutOfBounds{..},InvalidBound=>ck::BoundsMatrixError::InvalidBound{..}
},{
    #[wasm_bindgen(getter)]pub fn dimension(&self)->Option<usize>{match *self.inner{ck::BoundsMatrixError::DimensionOverflow{dimension}|ck::BoundsMatrixError::IndexOutOfBounds{dimension,..}=>Some(dimension),_=>None}}
    #[wasm_bindgen(getter)]pub fn elements(&self)->Option<usize>{match *self.inner{ck::BoundsMatrixError::AllocationFailed{elements}=>Some(elements),_=>None}}
    #[wasm_bindgen(getter)]pub fn expected(&self)->Option<usize>{match *self.inner{ck::BoundsMatrixError::DataLengthTooShort{expected,..}=>Some(expected),_=>None}}
    #[wasm_bindgen(getter)]pub fn actual(&self)->Option<usize>{match *self.inner{ck::BoundsMatrixError::DataLengthTooShort{actual,..}=>Some(actual),_=>None}}
    #[wasm_bindgen(getter)]pub fn index(&self)->Option<usize>{match *self.inner{ck::BoundsMatrixError::IndexOutOfBounds{index,..}=>Some(index),_=>None}}
    #[wasm_bindgen(getter)]pub fn axis(&self)->Option<MatrixAxis>{match *self.inner{ck::BoundsMatrixError::IndexOutOfBounds{axis,..}=>Some(axis.into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn kind(&self)->Option<BoundKind>{match *self.inner{ck::BoundsMatrixError::InvalidBound{kind,..}=>Some(kind.into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn row(&self)->Option<usize>{match *self.inner{ck::BoundsMatrixError::InvalidBound{row,..}=>Some(row),_=>None}}
    #[wasm_bindgen(getter)]pub fn column(&self)->Option<usize>{match *self.inner{ck::BoundsMatrixError::InvalidBound{column,..}=>Some(column),_=>None}}
});
enum_error!(EmbedSeedError,EmbedSeedErrorKind,ck::EmbedSeedError,{InvalidRandomSeed=>ck::EmbedSeedError::InvalidRandomSeed{..},ResultInvariant=>ck::EmbedSeedError::ResultInvariant{..}},{
    #[wasm_bindgen(getter)]pub fn seed(&self)->i32{match *self.inner{ck::EmbedSeedError::InvalidRandomSeed{seed}|ck::EmbedSeedError::ResultInvariant{seed,..}=>seed}}
    #[wasm_bindgen(getter)]pub fn index(&self)->Option<usize>{match *self.inner{ck::EmbedSeedError::ResultInvariant{index,..}=>Some(index),_=>None}}
    #[wasm_bindgen(getter)]pub fn value(&self)->Option<i32>{match *self.inner{ck::EmbedSeedError::ResultInvariant{value,..}=>Some(value),_=>None}}
});
enum_error!(GraphBoundsError,GraphBoundsErrorKind,ck::GraphBoundsError,{
    Bounds=>ck::GraphBoundsError::Bounds(_),Uff=>ck::GraphBoundsError::Uff(_),Matrix=>ck::GraphBoundsError::Matrix(_),Valence=>ck::GraphBoundsError::Valence(_),InvalidBounds=>ck::GraphBoundsError::InvalidBounds(_),Detail=>ck::GraphBoundsError::Detail(_),Input=>ck::GraphBoundsError::Input(_),AtomicNumber=>ck::GraphBoundsError::AtomicNumber(_)
},{
    #[wasm_bindgen(getter)]pub fn error(&self)->Option<ConformerCause>{payload!(self.inner,GraphBoundsError,[Bounds,Uff,Matrix,Valence])}
    #[wasm_bindgen(getter)]pub fn message(&self)->Option<String>{match self.inner.as_ref(){ck::GraphBoundsError::InvalidBounds(v)|ck::GraphBoundsError::Detail(v)=>Some(v.clone()),ck::GraphBoundsError::Input(v)=>Some((*v).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=atomicNumber)]pub fn atomic_number(&self)->Option<u8>{match *self.inner{ck::GraphBoundsError::AtomicNumber(v)=>Some(v),_=>None}}
},source=ConformerCause,cause=cause::<ck::GraphBoundsError>);
enum_error!(BoundsQueryError,BoundsQueryErrorKind,ck::BoundsQueryError,{
    Topology=>ck::BoundsQueryError::Topology(_),Valence=>ck::BoundsQueryError::Valence(_),Rings=>ck::BoundsQueryError::Rings(_),Conjugation=>ck::BoundsQueryError::Conjugation(_),Hybridization=>ck::BoundsQueryError::Hybridization(_),Bounds=>ck::BoundsQueryError::Bounds(_)
},{#[wasm_bindgen(getter)]pub fn error(&self)->Option<ConformerCause>{payload!(self.inner,BoundsQueryError,[Topology,Valence,Rings,Conjugation,Hybridization,Bounds])}},source=ConformerCause,cause=cause::<ck::BoundsQueryError>);
enum_error!(PruningError,PruningErrorKind,ck::PruningError,{
    Hydrogens=>ck::PruningError::Hydrogens(_),Projection=>ck::PruningError::Projection(_),Terminal=>ck::PruningError::Terminal(_),Match=>ck::PruningError::Match(_),Context=>ck::PruningError::Context(_),Input=>ck::PruningError::Input(_)
},{#[wasm_bindgen(getter)]pub fn input(&self)->Option<String>{match self.inner.as_ref(){ck::PruningError::Input(v)=>Some((*v).into()),_=>None}}#[wasm_bindgen(getter)]pub fn error(&self)->Option<ConformerCause>{payload!(self.inner,PruningError,[Hydrogens,Projection,Terminal,Match,Context])}},source=ConformerCause,cause=cause::<ck::PruningError>);
enum_error!(GenerationFailure,GenerationFailureKind,ck::GenerationFailure,{TopologyBounds=>ck::GenerationFailure::TopologyBounds(_),TorsionPreferences=>ck::GenerationFailure::TorsionPreferences(_)},{#[wasm_bindgen(getter)]pub fn error(&self)->Option<ConformerCause>{payload!(self.inner,GenerationFailure,[TopologyBounds,TorsionPreferences])}},source=ConformerCause,cause=cause::<ck::GenerationFailure>);
enum_error!(GenerationError,GenerationErrorKind,ck::GenerationError,{
    Topology=>ck::GenerationError::Topology(_),Coordinates=>ck::GenerationError::Coordinates(_),Valence=>ck::GenerationError::Valence(_),Rings=>ck::GenerationError::Rings(_),Conjugation=>ck::GenerationError::Conjugation(_),Hybridization=>ck::GenerationError::Hybridization(_),Stereo=>ck::GenerationError::Stereo(_),Fragments=>ck::GenerationError::Fragments(_),Hydrogens=>ck::GenerationError::Hydrogens(_),Pruning=>ck::GenerationError::Pruning(_),Seed=>ck::GenerationError::Seed(_),GenerationFailed=>ck::GenerationError::GenerationFailed(_),Atropisomer=>ck::GenerationError::Atropisomer(_),Bounds=>ck::GenerationError::Bounds(_),ForceField=>ck::GenerationError::ForceField(_),Numeric=>ck::GenerationError::Numeric(_),ThreadCount=>ck::GenerationError::ThreadCount(_),WorkerLaunch=>ck::GenerationError::WorkerLaunch(_),Input=>ck::GenerationError::Input(_),CoordinateMapIndex=>ck::GenerationError::CoordinateMapIndex{..},CoordinateMapNonFinite=>ck::GenerationError::CoordinateMapNonFinite{..},UnsupportedImplementation=>ck::GenerationError::UnsupportedImplementation{..}
},{
    #[wasm_bindgen(getter)]pub fn error(&self)->Option<ConformerCause>{payload!(self.inner,GenerationError,[Topology,Coordinates,Valence,Rings,Conjugation,Hybridization,Stereo,Fragments,Hydrogens,Pruning,Seed,GenerationFailed,Atropisomer,Bounds,ForceField,Numeric,ThreadCount,WorkerLaunch])}
    #[wasm_bindgen(getter)]pub fn input(&self)->Option<String>{match self.inner.as_ref(){ck::GenerationError::Input(v)=>Some((*v).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=atomIndex)]pub fn atom_index(&self)->Option<i32>{match *self.inner{ck::GenerationError::CoordinateMapIndex{atom_index,..}|ck::GenerationError::CoordinateMapNonFinite{atom_index}=>Some(atom_index),_=>None}}
    #[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match *self.inner{ck::GenerationError::CoordinateMapIndex{atom_count,..}=>Some(atom_count),_=>None}}
    #[wasm_bindgen(getter)]pub fn implementation(&self)->Option<String>{match self.inner.as_ref(){ck::GenerationError::UnsupportedImplementation{implementation}=>Some((*implementation).into()),_=>None}}
},source=ConformerCause,cause=cause::<ck::GenerationError>);

#[wasm_bindgen]
pub enum UffBoundsCauseKind {
    Parameter,
    BondOrder,
    BondMath,
    AssignmentLength,
}
#[wasm_bindgen]
pub struct UffBoundsCause {
    inner: ReadonlyValue<ck::UffBoundsError>,
}
#[wasm_bindgen]
impl UffBoundsCause {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> UffBoundsCauseKind {
        match self.inner.cause() {
            ck::UffBoundsCause::Parameter(_) => UffBoundsCauseKind::Parameter,
            ck::UffBoundsCause::BondOrder(_) => UffBoundsCauseKind::BondOrder,
            ck::UffBoundsCause::BondMath(_) => UffBoundsCauseKind::BondMath,
            ck::UffBoundsCause::AssignmentLength { .. } => UffBoundsCauseKind::AssignmentLength,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn field(&self) -> Option<String> {
        match self.inner.cause() {
            ck::UffBoundsCause::AssignmentLength { field, .. } => Some((*field).into()),
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn actual(&self) -> Option<usize> {
        match self.inner.cause() {
            ck::UffBoundsCause::AssignmentLength { actual, .. } => Some(*actual),
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn expected(&self) -> Option<usize> {
        match self.inner.cause() {
            ck::UffBoundsCause::AssignmentLength { expected, .. } => Some(*expected),
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn error(&self) -> Option<ConformerCause> {
        cause(&self.inner)
    }
}
error_projection!(
    UffBoundsError,
    ck::UffBoundsError,
    {
        #[wasm_bindgen(getter)]
        pub fn cause(&self) -> UffBoundsCause {
            UffBoundsCause {
                inner: self.inner.clone(),
            }
        }
    },
    source = ConformerCause,
    cause = cause::<ck::UffBoundsError>
);
#[wasm_bindgen]
pub enum ConformerOptimizerCauseKind {
    Input,
    Construction,
}
#[wasm_bindgen]
pub enum DistGeomForceFieldErrorKind {
    Kernel,
    Inversion,
    IndexOverflow,
}
#[wasm_bindgen]
pub struct DistGeomForceFieldError {
    inner: ReadonlyValue<ck::DistGeomForceFieldError>,
}
#[wasm_bindgen]
impl DistGeomForceFieldError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> DistGeomForceFieldErrorKind {
        match *self.inner {
            ck::DistGeomForceFieldError::Kernel(_) => DistGeomForceFieldErrorKind::Kernel,
            ck::DistGeomForceFieldError::Inversion(_) => DistGeomForceFieldErrorKind::Inversion,
            ck::DistGeomForceFieldError::IndexOverflow(_) => {
                DistGeomForceFieldErrorKind::IndexOverflow
            }
        }
    }
    #[wasm_bindgen(getter)]
    pub fn index(&self) -> Option<usize> {
        match *self.inner {
            ck::DistGeomForceFieldError::IndexOverflow(v) => Some(v),
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn error(&self) -> Option<ConformerCause> {
        match *self.inner {
            ck::DistGeomForceFieldError::Kernel(_) => {
                Some(project(self.inner.child(|v| match v {
                    ck::DistGeomForceFieldError::Kernel(e) => e as &dyn Error,
                    _ => unreachable!("immutable variant"),
                })))
            }
            ck::DistGeomForceFieldError::Inversion(_) => {
                Some(project(self.inner.child(|v| match v {
                    ck::DistGeomForceFieldError::Inversion(e) => e as &dyn Error,
                    _ => unreachable!("immutable variant"),
                })))
            }
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct ConformerOptimizerCause {
    inner: ReadonlyValue<ck::ConformerOptimizerCause>,
}
#[wasm_bindgen]
impl ConformerOptimizerCause {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> ConformerOptimizerCauseKind {
        match *self.inner {
            ck::ConformerOptimizerCause::Input(_) => ConformerOptimizerCauseKind::Input,
            ck::ConformerOptimizerCause::Construction(_) => {
                ConformerOptimizerCauseKind::Construction
            }
        }
    }
    #[wasm_bindgen(getter)]
    pub fn input(&self) -> Option<String> {
        match *self.inner {
            ck::ConformerOptimizerCause::Input(v) => Some(v.into()),
            _ => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn construction(&self) -> Option<DistGeomForceFieldError> {
        match *self.inner {
            ck::ConformerOptimizerCause::Construction(_) => Some(DistGeomForceFieldError {
                inner: self.inner.child(|v| match v {
                    ck::ConformerOptimizerCause::Construction(e) => e,
                    _ => unreachable!("immutable variant"),
                }),
            }),
            _ => None,
        }
    }
}
error_projection!(
    ConformerOptimizerError,
    ck::ConformerOptimizerError,
    {
        #[wasm_bindgen(getter)]
        pub fn cause(&self) -> ConformerOptimizerCause {
            ConformerOptimizerCause {
                inner: self.inner.child(ck::ConformerOptimizerError::cause),
            }
        }
    },
    source = ConformerCause,
    cause = cause::<ck::ConformerOptimizerError>
);
enum_error!(CrystalffTorsionPreferencesError,CrystalffTorsionPreferencesErrorKind,ck::CrystalffTorsionPreferencesError,{
    InvalidExperimentalTorsionVersion=>ck::CrystalffTorsionPreferencesError::InvalidExperimentalTorsionVersion{..},MissingSmartsPattern=>ck::CrystalffTorsionPreferencesError::MissingSmartsPattern{..},IncompleteParameterLine=>ck::CrystalffTorsionPreferencesError::IncompleteParameterLine{..},InvalidIntegerToken=>ck::CrystalffTorsionPreferencesError::InvalidIntegerToken{..},InvalidFloatToken=>ck::CrystalffTorsionPreferencesError::InvalidFloatToken{..},InvalidParameterSourceFormat=>ck::CrystalffTorsionPreferencesError::InvalidParameterSourceFormat{..},UnsupportedEscapeSequence=>ck::CrystalffTorsionPreferencesError::UnsupportedEscapeSequence{..},UnterminatedQuotedLiteral=>ck::CrystalffTorsionPreferencesError::UnterminatedQuotedLiteral{..},TrailingEscape=>ck::CrystalffTorsionPreferencesError::TrailingEscape{..},EmptyMolecule=>ck::CrystalffTorsionPreferencesError::EmptyMolecule,PreparedTarget=>ck::CrystalffTorsionPreferencesError::PreparedTarget(_),Match=>ck::CrystalffTorsionPreferencesError::Match(_),QueryBuildFailed=>ck::CrystalffTorsionPreferencesError::QueryBuildFailed{..},RingInfoUnavailable=>ck::CrystalffTorsionPreferencesError::RingInfoUnavailable{..},MissingCentralBond=>ck::CrystalffTorsionPreferencesError::MissingCentralBond{..}
},{
    #[wasm_bindgen(getter)]pub fn error(&self)->Option<ConformerCause>{payload!(self.inner,CrystalffTorsionPreferencesError,[PreparedTarget,Match])}
    #[wasm_bindgen(getter)]pub fn version(&self)->Option<u32>{match *self.inner{ck::CrystalffTorsionPreferencesError::InvalidExperimentalTorsionVersion{version}=>Some(version),_=>None}}
    #[wasm_bindgen(getter)]pub fn line(&self)->Option<String>{match self.inner.as_ref(){ck::CrystalffTorsionPreferencesError::MissingSmartsPattern{line}|ck::CrystalffTorsionPreferencesError::IncompleteParameterLine{line}|ck::CrystalffTorsionPreferencesError::InvalidIntegerToken{line,..}|ck::CrystalffTorsionPreferencesError::InvalidFloatToken{line,..}=>Some(line.clone()),_=>None}}
    #[wasm_bindgen(getter)]pub fn token(&self)->Option<String>{match self.inner.as_ref(){ck::CrystalffTorsionPreferencesError::InvalidIntegerToken{token,..}|ck::CrystalffTorsionPreferencesError::InvalidFloatToken{token,..}=>Some(token.clone()),_=>None}}
    #[wasm_bindgen(getter,js_name=constantName)]pub fn constant_name(&self)->Option<String>{match self.inner.as_ref(){ck::CrystalffTorsionPreferencesError::InvalidParameterSourceFormat{constant_name}|ck::CrystalffTorsionPreferencesError::UnsupportedEscapeSequence{constant_name,..}|ck::CrystalffTorsionPreferencesError::UnterminatedQuotedLiteral{constant_name}|ck::CrystalffTorsionPreferencesError::TrailingEscape{constant_name}=>Some((*constant_name).into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn escape(&self)->Option<String>{match *self.inner{ck::CrystalffTorsionPreferencesError::UnsupportedEscapeSequence{escape,..}=>Some(escape.to_string()),_=>None}}
    #[wasm_bindgen(getter)]pub fn smarts(&self)->Option<String>{match self.inner.as_ref(){ck::CrystalffTorsionPreferencesError::QueryBuildFailed{smarts,..}|ck::CrystalffTorsionPreferencesError::MissingCentralBond{smarts,..}=>Some(smarts.clone()),_=>None}}
    #[wasm_bindgen(getter)]pub fn reason(&self)->Option<String>{match self.inner.as_ref(){ck::CrystalffTorsionPreferencesError::QueryBuildFailed{reason,..}|ck::CrystalffTorsionPreferencesError::RingInfoUnavailable{reason}=>Some(reason.clone()),_=>None}}
    #[wasm_bindgen(getter)]pub fn aid2(&self)->Option<usize>{match *self.inner{ck::CrystalffTorsionPreferencesError::MissingCentralBond{aid2,..}=>Some(aid2),_=>None}}
    #[wasm_bindgen(getter)]pub fn aid3(&self)->Option<usize>{match *self.inner{ck::CrystalffTorsionPreferencesError::MissingCentralBond{aid3,..}=>Some(aid3),_=>None}}
},source=ConformerCause,cause=cause::<ck::CrystalffTorsionPreferencesError>);
