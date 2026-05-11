#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryNode<T> {
    Predicate(T),
    And(Vec<QueryNode<T>>),
    Or(Vec<QueryNode<T>>),
    Not(Box<QueryNode<T>>),
}

impl<T> QueryNode<T> {
    #[must_use]
    pub fn predicate(predicate: T) -> Self {
        Self::Predicate(predicate)
    }

    #[must_use]
    pub fn and(children: Vec<QueryNode<T>>) -> Self {
        Self::And(children)
    }

    #[must_use]
    pub fn or(children: Vec<QueryNode<T>>) -> Self {
        Self::Or(children)
    }

    #[must_use]
    pub fn not(child: QueryNode<T>) -> Self {
        Self::Not(Box::new(child))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AtomQueryPredicate {
    Any,
    AtomicNumber(u8),
    AtomicNumberIn(Vec<u8>),
    AtomicNumberNotIn(Vec<u8>),
    FormalCharge(i8),
    Isotope(u16),
    ImplicitHydrogenCount(u8),
    ImplicitHydrogenCountLessEqual(u8),
    ExplicitDegree(u8),
    ExplicitDegreeLessEqual(u8),
    NonHydrogenDegree(u32),
    RingBondCount(u8),
    RingBondCountLessEqual(u8),
    RingBondCountNeedsScan,
    IsAromatic(bool),
    IsUnsaturated,
    RecursiveSmarts(String),
    RGroupLabel(u32),
    MolFileAlias(String),
    UnsupportedFeature(&'static str),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BondQueryPredicate {
    Any,
    Order(crate::BondOrder),
    OrderIn(Vec<crate::BondOrder>),
    IsAromatic(bool),
    IsInRing(bool),
    MolFileQueryCode(u32),
    UnsupportedFeature(&'static str),
}
