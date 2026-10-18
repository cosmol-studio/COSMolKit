//! Readonly native query values retained by a real SMARTS writer failure.
use crate::chemical_values::{BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization};
use crate::error_framework::Source;
use crate::query_compiler::QueryGraph;
use crate::retained_error::RetainedError;
use cosmolkit as ck;
use pyo3::prelude::*;

// This is binding representation, not query interpretation. Exhaustive matches
// keep the Python discriminants in sync with the actual native variants.
macro_rules! value_kind {
    ($kind:ident, $rust:ident, {$($tag:ident $({$($field:ident),+ $(,)?})? $(($value:tt))?),+ $(,)?}) => {
        #[pyclass(module="fixed_projection._native", eq, eq_int, from_py_object)]
        #[derive(Clone, Copy, PartialEq)]
        pub enum $kind { $($tag),+ }
        impl From<&ck::$rust> for $kind {
            fn from(value: &ck::$rust) -> Self {
                match value { $($crate::error_framework::variant_pattern!(ck::$rust, $tag $({$($field),+})? $(($value))?) => Self::$tag),+ }
            }
        }
    };
}
pub(crate) use value_kind;

value_kind!(AtomQueryPredicateKind, AtomQueryPredicate, {
    Any, AtomicNumber(_), AtomType { atomic_number, aromatic }, AtomicNumberIn(_),
    AtomicNumberNotIn(_), FormalCharge(_), NegativeFormalCharge(_), NumRadicalElectrons(_),
    HasChiralTag, MissingChiralTag, Isotope(_), HydrogenCount(_), HasImplicitHydrogen,
    ImplicitHydrogenCount(_), ImplicitHydrogenCountLessEqual(_), ImplicitValence(_),
    ExplicitValence(_), ExplicitDegree(_), ExplicitDegreeLessEqual(_), NonHydrogenDegree(_),
    NonHydrogenDegreeLessEqual(_), NonHydrogenDegreeGreaterEqual(_), HeavyAtomDegree(_),
    NumHeteroatomNeighbors(_), HasHeteroatomNeighbors, NumAliphaticHeteroatomNeighbors(_),
    HasAliphaticHeteroatomNeighbors, RingBondCount(_), RingBondCountLessEqual(_), HasRingBond,
    IsBridgehead, IsAromatic(_), IsUnsaturated, RecursiveSmarts(_), HasProperty(_),
    PropertyValue { name, value }, RGroupLabel(_), MolFileAlias(_), HybridizationMatch(_),
    TotalDegree(_), TotalDegreeLessEqual(_), TotalDegreeGreaterEqual(_), TotalValence(_),
    TotalValenceLessEqual(_), TotalValenceGreaterEqual(_), InRing, NumAtomRings(_),
    InRingOfSize(_), InRingOfSizeLessEqual(_), InRingOfSizeGreaterEqual(_), SmallestRingSize(_),
    SmallestRingSizeLessEqual(_), SmallestRingSizeGreaterEqual(_), Mass(_), ChiralTagMatch(_),
    ChiralPermutationMatch(_), DegreeLessEqual(_), DegreeGreaterEqual(_), Range(_), UnsupportedFeature(_),
});
value_kind!(BondQueryPredicateKind, BondQueryPredicate, {
    Any, Order(_), OrderIn(_), IsAromatic(_), IsInRing(_), Direction(_), Stereo(_), HasStereo,
    IsConjugated, NumRingBonds(_), InRingOfSize(_), MinRingSize(_), NumRingBondsGreaterEqual(_),
    NumRingBondsLessEqual(_), MolFileQueryCode(_), HasProperty(_), PropertyValue { name, value },
    UnsupportedFeature(_),
});
value_kind!(AtomRangeBoundsKind, AtomRangeBounds, {
    LessEqual(_), GreaterEqual(_), Inclusive { lower, upper, lower_open, upper_open },
});
value_kind!(AtomRangeDataFunctionKind, AtomRangeDataFunction, {
    ExplicitDegree, NonHydrogenDegree, TotalDegree, TotalValence, NumAtomRings,
    NumHeteroatomNeighbors, NumAliphaticHeteroatomNeighbors, MinRingSize, RingBondCount,
    ImplicitHydrogenCount, FormalCharge, NegativeFormalCharge,
    AtomRingSize { lower, upper, lower_open, upper_open },
});

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct AtomQueryPredicate {
    pub(crate) owner: RetainedError<ck::SmartsWriteError>,
}
impl AtomQueryPredicate {
    fn inner(&self) -> &ck::AtomQueryPredicate {
        let ck::SmartsWriteError::UnsupportedAtomQuery { predicate } = self.owner.as_ref() else {
            unreachable!("atom predicate view retains its native variant")
        };
        predicate
    }
}
type AtomValue = Source<
    i64,
    Source<
        bool,
        Source<
            Vec<u8>,
            Source<
                String,
                Source<
                    Hybridization,
                    Source<ChiralTag, Source<RecursiveStructureQuery, AtomRangeQuery>>,
                >,
            >,
        >,
    >,
>;
#[pymethods]
impl AtomQueryPredicate {
    #[getter]
    fn variant(&self) -> AtomQueryPredicateKind {
        self.inner().into()
    }
    #[getter]
    fn atomic_number(&self) -> Option<u8> {
        match self.inner() {
            ck::AtomQueryPredicate::AtomType { atomic_number, .. } => Some(*atomic_number),
            _ => None,
        }
    }
    #[getter]
    fn aromatic(&self) -> Option<bool> {
        match self.inner() {
            ck::AtomQueryPredicate::AtomType { aromatic, .. } => Some(*aromatic),
            _ => None,
        }
    }
    #[getter]
    fn name(&self) -> Option<String> {
        match self.inner() {
            ck::AtomQueryPredicate::PropertyValue { name, .. } => Some(name.clone()),
            _ => None,
        }
    }
    #[getter]
    fn value(&self) -> Option<AtomValue> {
        use ck::AtomQueryPredicate as P;
        match self.inner() {
            P::FormalCharge(v)
            | P::NegativeFormalCharge(v)
            | P::Isotope(v)
            | P::HydrogenCount(v)
            | P::ImplicitHydrogenCount(v)
            | P::ImplicitValence(v)
            | P::ExplicitValence(v)
            | P::ExplicitDegree(v)
            | P::NumHeteroatomNeighbors(v)
            | P::NumAliphaticHeteroatomNeighbors(v)
            | P::RingBondCount(v)
            | P::TotalDegree(v)
            | P::TotalValence(v)
            | P::NumAtomRings(v)
            | P::InRingOfSize(v)
            | P::SmallestRingSize(v) => Some(Source::First(i64::from(*v))),
            P::AtomicNumber(v)
            | P::NumRadicalElectrons(v)
            | P::ImplicitHydrogenCountLessEqual(v)
            | P::ExplicitDegreeLessEqual(v)
            | P::RingBondCountLessEqual(v)
            | P::TotalDegreeLessEqual(v)
            | P::TotalDegreeGreaterEqual(v)
            | P::TotalValenceLessEqual(v)
            | P::TotalValenceGreaterEqual(v)
            | P::InRingOfSizeLessEqual(v)
            | P::InRingOfSizeGreaterEqual(v)
            | P::SmallestRingSizeLessEqual(v)
            | P::SmallestRingSizeGreaterEqual(v)
            | P::DegreeLessEqual(v)
            | P::DegreeGreaterEqual(v) => Some(Source::First(i64::from(*v))),
            P::NonHydrogenDegree(v)
            | P::NonHydrogenDegreeLessEqual(v)
            | P::NonHydrogenDegreeGreaterEqual(v)
            | P::HeavyAtomDegree(v)
            | P::RGroupLabel(v)
            | P::ChiralPermutationMatch(v) => Some(Source::First(i64::from(*v))),
            P::Mass(v) => Some(Source::First(i64::from(*v))),
            P::IsAromatic(v) => Some(Source::Second(Source::First(*v))),
            P::AtomicNumberIn(v) | P::AtomicNumberNotIn(v) => {
                Some(Source::Second(Source::Second(Source::First(v.clone()))))
            }
            P::HasProperty(v) | P::MolFileAlias(v) | P::PropertyValue { value: v, .. } => Some(
                Source::Second(Source::Second(Source::Second(Source::First(v.clone())))),
            ),
            P::UnsupportedFeature(v) => Some(Source::Second(Source::Second(Source::Second(
                Source::First((*v).into()),
            )))),
            P::HybridizationMatch(v) => Some(Source::Second(Source::Second(Source::Second(
                Source::Second(Source::First((*v).into())),
            )))),
            P::ChiralTagMatch(v) => Some(Source::Second(Source::Second(Source::Second(
                Source::Second(Source::Second(Source::First((*v).into()))),
            )))),
            P::RecursiveSmarts(_) => Some(Source::Second(Source::Second(Source::Second(
                Source::Second(Source::Second(Source::Second(Source::First(
                    RecursiveStructureQuery {
                        owner: self.owner.clone(),
                    },
                )))),
            )))),
            P::Range(v) => Some(Source::Second(Source::Second(Source::Second(
                Source::Second(Source::Second(Source::Second(Source::Second(
                    AtomRangeQuery { inner: v.clone() },
                )))),
            )))),
            P::Any
            | P::AtomType { .. }
            | P::HasChiralTag
            | P::MissingChiralTag
            | P::HasImplicitHydrogen
            | P::HasHeteroatomNeighbors
            | P::HasAliphaticHeteroatomNeighbors
            | P::HasRingBond
            | P::IsBridgehead
            | P::IsUnsaturated
            | P::InRing => None,
        }
    }
}

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct BondQueryPredicate {
    pub(crate) owner: RetainedError<ck::SmartsWriteError>,
}
type BondValue = Source<
    BondOrder,
    Source<
        Vec<BondOrder>,
        Source<bool, Source<BondDirection, Source<BondStereo, Source<i64, String>>>>,
    >,
>;
impl BondQueryPredicate {
    fn inner(&self) -> &ck::BondQueryPredicate {
        let ck::SmartsWriteError::UnsupportedBondQuery { predicate } = self.owner.as_ref() else {
            unreachable!("bond predicate view retains its native variant")
        };
        predicate
    }
}
#[pymethods]
impl BondQueryPredicate {
    #[getter]
    fn variant(&self) -> BondQueryPredicateKind {
        self.inner().into()
    }
    #[getter]
    fn name(&self) -> Option<String> {
        match self.inner() {
            ck::BondQueryPredicate::PropertyValue { name, .. } => Some(name.clone()),
            _ => None,
        }
    }
    #[getter]
    fn value(&self) -> Option<BondValue> {
        use ck::BondQueryPredicate as P;
        match self.inner() {
            P::Order(v) => Some(Source::First((*v).into())),
            P::OrderIn(v) => Some(Source::Second(Source::First(
                v.iter().copied().map(Into::into).collect(),
            ))),
            P::IsAromatic(v) | P::IsInRing(v) => {
                Some(Source::Second(Source::Second(Source::First(*v))))
            }
            P::Direction(v) => Some(Source::Second(Source::Second(Source::Second(
                Source::First((*v).into()),
            )))),
            P::Stereo(v) => Some(Source::Second(Source::Second(Source::Second(
                Source::Second(Source::First((*v).into())),
            )))),
            P::NumRingBonds(v) | P::InRingOfSize(v) | P::MinRingSize(v) => {
                Some(Source::Second(Source::Second(Source::Second(
                    Source::Second(Source::Second(Source::First(i64::from(*v)))),
                ))))
            }
            P::NumRingBondsGreaterEqual(v) | P::NumRingBondsLessEqual(v) => {
                Some(Source::Second(Source::Second(Source::Second(
                    Source::Second(Source::Second(Source::First(i64::from(*v)))),
                ))))
            }
            P::MolFileQueryCode(v) => Some(Source::Second(Source::Second(Source::Second(
                Source::Second(Source::Second(Source::First(i64::from(*v)))),
            )))),
            P::HasProperty(v) | P::PropertyValue { value: v, .. } => {
                Some(Source::Second(Source::Second(Source::Second(
                    Source::Second(Source::Second(Source::Second(v.clone()))),
                ))))
            }
            P::UnsupportedFeature(v) => Some(Source::Second(Source::Second(Source::Second(
                Source::Second(Source::Second(Source::Second((*v).into()))),
            )))),
            P::Any | P::HasStereo | P::IsConjugated => None,
        }
    }
}

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct RecursiveStructureQuery {
    owner: RetainedError<ck::SmartsWriteError>,
}
impl RecursiveStructureQuery {
    fn inner(&self) -> &ck::RecursiveStructureQuery {
        let ck::SmartsWriteError::UnsupportedAtomQuery {
            predicate: ck::AtomQueryPredicate::RecursiveSmarts(inner),
        } = self.owner.as_ref()
        else {
            unreachable!("recursive predicate view retains its native variant")
        };
        inner
    }
}
#[pymethods]
impl RecursiveStructureQuery {
    #[getter]
    fn query_graph(&self) -> Option<QueryGraph> {
        self.inner().query_graph()?;
        Some(QueryGraph::from_writer_recursive(self.owner.clone()))
    }
    #[getter]
    fn source_smarts(&self, py: Python<'_>) -> PyResult<Option<String>> {
        self.inner()
            .source_smarts()
            .map(|value| crate::property_values::decode_source_text(py, value))
            .transpose()
    }
    #[getter]
    fn serial_number(&self) -> u32 {
        self.inner().serial_number()
    }
    #[getter]
    fn initialized(&self) -> bool {
        self.inner().is_initialized()
    }
}

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct AtomRangeQuery {
    inner: ck::AtomRangeQuery,
}
#[pymethods]
impl AtomRangeQuery {
    #[getter]
    fn bounds(&self) -> AtomRangeBounds {
        AtomRangeBounds {
            inner: self.inner.bounds(),
        }
    }
    #[getter]
    fn data_function(&self) -> AtomRangeDataFunction {
        AtomRangeDataFunction {
            inner: self.inner.data_function(),
        }
    }
}

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct AtomRangeBounds {
    inner: ck::AtomRangeBounds,
}
#[pymethods]
impl AtomRangeBounds {
    #[getter]
    fn variant(&self) -> AtomRangeBoundsKind {
        (&self.inner).into()
    }
    #[getter]
    fn value(&self) -> Option<i32> {
        match self.inner {
            ck::AtomRangeBounds::LessEqual(v) | ck::AtomRangeBounds::GreaterEqual(v) => Some(v),
            _ => None,
        }
    }
    #[getter]
    fn lower(&self) -> Option<i32> {
        match self.inner {
            ck::AtomRangeBounds::Inclusive { lower, .. } => Some(lower),
            _ => None,
        }
    }
    #[getter]
    fn upper(&self) -> Option<i32> {
        match self.inner {
            ck::AtomRangeBounds::Inclusive { upper, .. } => Some(upper),
            _ => None,
        }
    }
    #[getter]
    fn lower_open(&self) -> Option<bool> {
        match self.inner {
            ck::AtomRangeBounds::Inclusive { lower_open, .. } => Some(lower_open),
            _ => None,
        }
    }
    #[getter]
    fn upper_open(&self) -> Option<bool> {
        match self.inner {
            ck::AtomRangeBounds::Inclusive { upper_open, .. } => Some(upper_open),
            _ => None,
        }
    }
}

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct AtomRangeDataFunction {
    inner: ck::AtomRangeDataFunction,
}
#[pymethods]
impl AtomRangeDataFunction {
    #[getter]
    fn variant(&self) -> AtomRangeDataFunctionKind {
        (&self.inner).into()
    }
    #[getter]
    fn lower(&self) -> Option<i32> {
        match self.inner {
            ck::AtomRangeDataFunction::AtomRingSize { lower, .. } => Some(lower),
            _ => None,
        }
    }
    #[getter]
    fn upper(&self) -> Option<i32> {
        match self.inner {
            ck::AtomRangeDataFunction::AtomRingSize { upper, .. } => Some(upper),
            _ => None,
        }
    }
    #[getter]
    fn lower_open(&self) -> Option<bool> {
        match self.inner {
            ck::AtomRangeDataFunction::AtomRingSize { lower_open, .. } => Some(lower_open),
            _ => None,
        }
    }
    #[getter]
    fn upper_open(&self) -> Option<bool> {
        match self.inner {
            ck::AtomRangeDataFunction::AtomRingSize { upper_open, .. } => Some(upper_open),
            _ => None,
        }
    }
}
