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
    Null,
    AtomicNumEquals(u8),
    AtomicNumIn(Vec<u8>),
    AtomicNumNotIn(Vec<u8>),
    Aromatic,
    Aliphatic,
    FormalChargeEquals(i8),
    IsotopeEquals(u16),
    ImplicitHCountEquals(u8),
    ImplicitHCountLessEqual(u8),
    RingBondCountEquals(u32),
    RingBondCountLessEqual(u32),
    RingBondCountNeedScan,
    Unsaturated,
    TotalValenceEquals(i32),
    SubstitutionCountEquals(i32),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BondQueryPredicate {
    Null,
    SingleOrDouble,
    SingleOrAromatic,
    DoubleOrAromatic,
    IsInRing,
}
