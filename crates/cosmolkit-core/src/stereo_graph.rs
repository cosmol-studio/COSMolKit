//! Borrowed graph access for the single drawing-stereo and ring algorithms.
//! Query predicates stay in their canonical carriers; no graph is materialized.

use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, Bond, BondId, NeighborRef, PropertyText, PropertyValue, QueryAtom,
    QueryBond, QueryGraph, TopologyBlock,
};
use cosmolkit_types::ChiralTag;
use std::collections::BTreeMap;

/// Read-only canonical atom attributes required by drawing stereo.
#[doc(hidden)]
pub trait StereoAtomAccess: sealed::AtomCarrier {
    fn id(&self) -> AtomId;
    fn chiral_tag(&self) -> ChiralTag;
    fn atomic_number(&self) -> u8;
    fn explicit_hydrogens(&self) -> u8;
    fn no_implicit(&self) -> bool;
    fn prop(&self, key: &str) -> Option<&PropertyValue>;
    fn props(&self) -> &BTreeMap<PropertyText, PropertyValue>;
}

macro_rules! stereo_atom_access {
    ($carrier:ty) => {
        impl StereoAtomAccess for $carrier {
            fn id(&self) -> AtomId {
                <$carrier>::id(self)
            }
            fn chiral_tag(&self) -> ChiralTag {
                <$carrier>::chiral_tag(self)
            }
            fn atomic_number(&self) -> u8 {
                <$carrier>::atomic_number(self)
            }
            fn explicit_hydrogens(&self) -> u8 {
                <$carrier>::explicit_hydrogens(self)
            }
            fn no_implicit(&self) -> bool {
                <$carrier>::no_implicit(self)
            }
            fn prop(&self, key: &str) -> Option<&PropertyValue> {
                <$carrier>::prop(self, key)
            }
            fn props(&self) -> &BTreeMap<PropertyText, PropertyValue> {
                <$carrier>::props(self)
            }
        }
    };
}
stereo_atom_access!(Atom);
stereo_atom_access!(QueryAtom);

#[derive(Clone, Copy)]
pub enum BondRows<'a> {
    Concrete(&'a [Bond]),
    Query(&'a [QueryBond]),
}
impl<'a> BondRows<'a> {
    pub fn len(self) -> usize {
        match self {
            Self::Concrete(v) => v.len(),
            Self::Query(v) => v.len(),
        }
    }
    pub fn is_empty(self) -> bool {
        self.len() == 0
    }
    pub fn get(self, index: usize) -> Option<&'a Bond> {
        match self {
            Self::Concrete(v) => v.get(index),
            Self::Query(v) => v.get(index).map(QueryBond::bond),
        }
    }
    pub fn iter(self) -> impl ExactSizeIterator<Item = &'a Bond> {
        (0..self.len()).map(move |i| self.get(i).expect("bounded borrowed bond row"))
    }
}
impl std::ops::Index<usize> for BondRows<'_> {
    type Output = Bond;
    fn index(&self, i: usize) -> &Bond {
        self.get(i).expect("borrowed bond index in range")
    }
}
impl<'a> IntoIterator for BondRows<'a> {
    type Item = &'a Bond;
    type IntoIter = BondRowIterator<'a>;
    fn into_iter(self) -> Self::IntoIter {
        BondRowIterator {
            rows: self,
            index: 0,
        }
    }
}
pub struct BondRowIterator<'a> {
    rows: BondRows<'a>,
    index: usize,
}
impl<'a> Iterator for BondRowIterator<'a> {
    type Item = &'a Bond;
    fn next(&mut self) -> Option<Self::Item> {
        let result = self.rows.get(self.index);
        self.index += usize::from(result.is_some());
        result
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.rows.len() - self.index;
        (remaining, Some(remaining))
    }
}
impl ExactSizeIterator for BondRowIterator<'_> {}
impl std::iter::FusedIterator for BondRowIterator<'_> {}

#[derive(Clone, Copy)]
pub enum NeighborRows<'a> {
    Concrete(&'a [NeighborRef]),
    Query(&'a [(usize, usize)]),
}
impl NeighborRows<'_> {
    pub fn len(self) -> usize {
        match self {
            Self::Concrete(v) => v.len(),
            Self::Query(v) => v.len(),
        }
    }
    pub fn is_empty(self) -> bool {
        self.len() == 0
    }
    pub fn iter(self) -> impl ExactSizeIterator<Item = NeighborRef> {
        (0..self.len()).map(move |i| match self {
            Self::Concrete(v) => v[i],
            Self::Query(v) => NeighborRef {
                atom_index: v[i].0,
                bond: BondId::new(v[i].1),
            },
        })
    }
}
#[derive(Clone, Copy)]
pub enum GraphAdjacency<'a> {
    Concrete(&'a AdjacencyList),
    Query(&'a QueryGraph),
}
impl<'a> GraphAdjacency<'a> {
    pub fn neighbors_of(self, atom: usize) -> NeighborRows<'a> {
        match self {
            Self::Concrete(v) => NeighborRows::Concrete(v.neighbors_of(atom)),
            Self::Query(v) => NeighborRows::Query(&v.adjacency()[atom]),
        }
    }
}

mod sealed {
    pub trait Graph {}
    impl Graph for cosmolkit_model::TopologyBlock {}
    impl Graph for cosmolkit_model::QueryGraph {}
    pub trait AtomCarrier {}
    impl AtomCarrier for cosmolkit_model::Atom {}
    impl AtomCarrier for cosmolkit_model::QueryAtom {}
}
#[derive(Debug, thiserror::Error)]
pub enum StereoGraphError {
    #[error(transparent)]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error(transparent)]
    Query(#[from] cosmolkit_model::QueryGraphError),
}
pub trait StereoGraphAccess: sealed::Graph {
    type Atom: StereoAtomAccess;
    fn atoms(&self) -> &[Self::Atom];
    fn bonds(&self) -> BondRows<'_>;
    fn adjacency(&self) -> GraphAdjacency<'_>;
    fn validate(&self) -> Result<(), StereoGraphError>;
}
impl StereoGraphAccess for TopologyBlock {
    type Atom = Atom;
    fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
    fn bonds(&self) -> BondRows<'_> {
        BondRows::Concrete(&self.bonds)
    }
    fn adjacency(&self) -> GraphAdjacency<'_> {
        GraphAdjacency::Concrete(&self.adjacency)
    }
    fn validate(&self) -> Result<(), StereoGraphError> {
        TopologyBlock::validate(self).map_err(Into::into)
    }
}
impl StereoGraphAccess for QueryGraph {
    type Atom = QueryAtom;
    fn atoms(&self) -> &[QueryAtom] {
        QueryGraph::atoms(self)
    }
    fn bonds(&self) -> BondRows<'_> {
        BondRows::Query(QueryGraph::bonds(self))
    }
    fn adjacency(&self) -> GraphAdjacency<'_> {
        GraphAdjacency::Query(self)
    }
    fn validate(&self) -> Result<(), StereoGraphError> {
        QueryGraph::validate(self).map_err(Into::into)
    }
}
