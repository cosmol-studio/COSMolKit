use crate::query::{AtomQueryPredicate, BondQueryPredicate, QueryNode};

pub fn push_atom_query(
    target: &mut Option<QueryNode<AtomQueryPredicate>>,
    query: QueryNode<AtomQueryPredicate>,
) {
    match target.take() {
        None => *target = Some(query),
        Some(QueryNode::And(mut children)) => {
            children.push(query);
            *target = Some(QueryNode::And(children));
        }
        Some(existing) => *target = Some(QueryNode::and(vec![existing, query])),
    }
}

pub fn push_bond_query(
    target: &mut Option<QueryNode<BondQueryPredicate>>,
    query: QueryNode<BondQueryPredicate>,
) {
    match target.take() {
        None => *target = Some(query),
        Some(QueryNode::And(mut children)) => {
            children.push(query);
            *target = Some(QueryNode::And(children));
        }
        Some(existing) => *target = Some(QueryNode::and(vec![existing, query])),
    }
}
