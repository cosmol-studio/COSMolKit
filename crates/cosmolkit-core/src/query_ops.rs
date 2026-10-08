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

// QueryNode::Not carries the canonical outer source negation flag, not a
// source child edge. Source description remains the underlying node identity.
fn source_bond_query_root(
    mut predicate: &QueryNode<BondQueryPredicate>,
) -> &QueryNode<BondQueryPredicate> {
    while let QueryNode::Not(child) = predicate {
        predicate = child;
    }
    predicate
}

fn bond_predicate_has_complex_type_query_with_seen(
    predicate: &QueryNode<BondQueryPredicate>,
    mut seen_bond_order: bool,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION QueryOps::hasComplexBondTypeQueryHelper complete source
    // RDKit✔️❌: bool hasComplexBondTypeQueryHelper(
    // RDKit✔️❌:     const Queries::Query<int, Bond const *, true> &qry, bool seenBondOrder) {
    // RDKit✔️❌:   const auto df = qry.getDescription();
    // RDKit✔️❌:   bool isBondOrder = (df == "BondOrder");
    // RDKit✔️❌:   // is this a bond order query?
    // RDKit✔️❌:   if (std::find(bondOrderQueryFunctions.begin(), bondOrderQueryFunctions.end(),
    // RDKit✔️❌:                 df) != bondOrderQueryFunctions.end()) {
    // RDKit✔️❌:     if (seenBondOrder || !isBondOrder || qry.getNegation()) {
    // RDKit✔️❌:       return true;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (const auto &child :
    // RDKit✔️❌:        boost::make_iterator_range(qry.beginChildren(), qry.endChildren())) {
    // RDKit✔️❌:     if (hasComplexBondTypeQueryHelper(*child, seenBondOrder | isBondOrder)) {
    // RDKit✔️❌:       return true;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (child->getDescription() == "BondOrder") {
    // RDKit✔️❌:       seenBondOrder = true;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return false;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION QueryOps::hasComplexBondTypeQueryHelper complete source
    // BEGIN RDKIT CPP FUNCTION QueryOps::bondOrderQueryFunctions complete source
    // RDKit✔️✔️: const std::vector<std::string> bondOrderQueryFunctions{
    // RDKit✔️✔️:     std::string("BondOrder"), std::string("SingleOrAromaticBond"),
    // RDKit✔️✔️:     std::string("DoubleOrAromaticBond"), std::string("SingleOrDoubleBond"),
    // RDKit✔️✔️:     std::string("SingleOrDoubleOrAromaticBond")};
    // END RDKIT CPP FUNCTION QueryOps::bondOrderQueryFunctions complete source
    // Order is exactly the plain BondOrder description. OrderIn is the
    // canonical legacy fixed-order-function family (the four existing
    // source-backed factories); its value list is not inspected by source
    // description classification. Other modeled leaves/Boolean roots have
    // descriptions outside that legacy family. Negation comes from the sole
    // typed flag accessor; composite negation does not propagate to children.
    // Source iteration and seen-state updates are deliberately ordered: only
    // a DIRECT child's BondOrder description updates later sibling state.
    // O(query nodes), O(depth) stack, no AST allocation/copy. Peeling typed
    // negation boxes and re-reading a child's root after recursion is extra
    // metadata traversal compared with native O(1) fields, marked worse.
    let root = source_bond_query_root(predicate);
    let is_bond_order = matches!(root, QueryNode::Predicate(BondQueryPredicate::Order(_)));
    let is_legacy_order = matches!(
        root,
        QueryNode::Predicate(BondQueryPredicate::Order(_) | BondQueryPredicate::OrderIn(_))
    );
    if is_legacy_order && (seen_bond_order || !is_bond_order || predicate.is_negated()) {
        return true;
    }
    if let QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) = root {
        for child in children {
            if bond_predicate_has_complex_type_query_with_seen(
                child,
                seen_bond_order | is_bond_order,
            ) {
                return true;
            }
            if matches!(
                source_bond_query_root(child),
                QueryNode::Predicate(BondQueryPredicate::Order(_))
            ) {
                seen_bond_order = true;
            }
        }
    }
    false
}

fn bond_predicate_has_complex_type_query(predicate: &QueryNode<BondQueryPredicate>) -> bool {
    // BEGIN RDKIT CPP FUNCTION QueryOps::hasComplexBondTypeQuery(query) complete source
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT bool hasComplexBondTypeQuery(
    // RDKit✔️✔️:     const Queries::Query<int, Bond const *, true> &qry) {
    // RDKit✔️✔️:   return hasComplexBondTypeQueryHelper(qry, false);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION QueryOps::hasComplexBondTypeQuery(query) complete source
    bond_predicate_has_complex_type_query_with_seen(predicate, false)
}

pub(crate) fn bond_has_complex_type_query(bond: &cosmolkit_model::Bond) -> bool {
    // BEGIN RDKIT CPP FUNCTION QueryOps::hasComplexBondTypeQuery(bond) complete source
    // RDKit✔️✔️: inline bool hasComplexBondTypeQuery(const Bond &bond) {
    // RDKit✔️✔️:   if (!bond.hasQuery()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return hasComplexBondTypeQuery(*bond.getQuery());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION QueryOps::hasComplexBondTypeQuery(bond) complete source
    // Existing explicit Bond.query is the source identity; ordinary unspecified
    // bond order is not manufactured into a query. No missing-query fallback.
    bond.query()
        .is_some_and(bond_predicate_has_complex_type_query)
}

/// Apply the sole source query classifier to the canonical query row.
#[doc(hidden)]
pub fn query_bond_has_complex_type_query(bond: &cosmolkit_model::QueryBond) -> bool {
    !bond.predicate_is_carrier_derived() && bond_predicate_has_complex_type_query(bond.predicate())
}
