use cosmolkit_model::{BondQueryPredicate, QueryNode};

pub(crate) fn bond_predicate_has_type_query(predicate: &QueryNode<BondQueryPredicate>) -> bool {
    // BEGIN RDKIT CPP FUNCTION QueryOps::hasBondTypeQuery
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT bool hasBondTypeQuery(
    // RDKit✔️✔️:     const Queries::Query<int, Bond const *, true> &qry) {
    // RDKit✔️✔️:   const auto df = qry.getDescription();
    // RDKit✔️✔️:   const auto dt = qry.getTypeLabel();
    // RDKit✔️✔️:   // is this a bond order query?
    // RDKit✔️✔️:   if (dt == "BondOrder" ||
    // RDKit✔️✔️:       (dt.empty() &&
    // RDKit✔️✔️:        std::find(bondOrderQueryFunctions.begin(), bondOrderQueryFunctions.end(),
    // RDKit✔️✔️:                  df) != bondOrderQueryFunctions.end())) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto &child :
    // RDKit✔️✔️:        boost::make_iterator_range(qry.beginChildren(), qry.endChildren())) {
    // RDKit✔️✔️:     if (hasBondTypeQuery(*child)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION QueryOps::hasBondTypeQuery
    // Behavior review: every modeled `Order`/`OrderIn` leaf carries the source
    // BondOrder type label; all Boolean node kinds recurse without changing
    // leaf identity, while Any and unrelated leaves remain false.
    // Complexity review: one depth-first visit over the query tree, O(nodes)
    // time and O(depth) call stack, matches the source recursive traversal and
    // performs no allocation or cloning.
    match predicate {
        QueryNode::Predicate(BondQueryPredicate::Order(_) | BondQueryPredicate::OrderIn(_)) => true,
        QueryNode::Predicate(_) => false,
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
            children.iter().any(bond_predicate_has_type_query)
        }
        QueryNode::Not(child) => bond_predicate_has_type_query(child),
    }
}

pub(crate) fn bond_predicate_is_order_query(predicate: &QueryNode<BondQueryPredicate>) -> bool {
    // BEGIN RDKIT CPP FUNCTION MolOps::isBondOrderQuery
    // RDKit✔️✔️: bool isBondOrderQuery(const Bond *bond) {
    // RDKit✔️✔️:   if (bond->hasQuery()) {
    // RDKit✔️✔️:     auto q = dynamic_cast<const QueryBond *>(bond)->getQuery();
    // RDKit✔️✔️:     // complex bond type queries are also bond order queries!
    // RDKit✔️✔️:     if (q->getTypeLabel() == "BondOrder" ||
    // RDKit✔️✔️:         QueryOps::hasComplexBondTypeQuery(*q)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::isBondOrderQuery
    // Behavior review: in the modeled AST a root Order/OrderIn leaf is the
    // source type-label case, and any Boolean wrapper containing such a leaf
    // is the modeled complex-bond-type case. Other explicit queries stay false.
    // Complexity review: delegation retains the source O(nodes) recursive
    // classification with no allocation or additional graph traversal.
    bond_predicate_has_type_query(predicate)
}
