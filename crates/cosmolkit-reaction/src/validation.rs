use crate::{Reaction, ReactionInitializationError, ReactionRole, ReactionValidationError};
use cosmolkit_model::{AtomId, AtomQueryPredicate, QueryAtom, QueryNode};
use std::collections::{BTreeMap, VecDeque};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReactionValidationParams {
    pub silent: bool,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionValidationSeverity {
    Warning,
    Error,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionValidationIssueKind {
    MissingReactants,
    MissingProducts,
    DuplicateReactantMap,
    UnmappedReactant,
    DuplicateProductMap,
    MissingProductMapReactant,
    UnmappedProduct,
    UnmappedReactantMaps,
    MultipleCharge,
    MultipleHydrogenCount,
    MultipleMass,
    MultipleIsotope,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionValidationIssue {
    pub kind: ReactionValidationIssueKind,
    pub severity: ReactionValidationSeverity,
    pub role: Option<ReactionRole>,
    pub template: Option<usize>,
    pub atom: Option<AtomId>,
    pub map: Option<i32>,
    pub maps: Vec<i32>,
    pub detail: String,
}
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReactionValidationReport {
    pub warnings: Vec<ReactionValidationIssue>,
    pub errors: Vec<ReactionValidationIssue>,
}
impl ReactionValidationReport {
    pub fn num_warnings(&self) -> usize {
        self.warnings.len()
    }
    pub fn num_errors(&self) -> usize {
        self.errors.len()
    }
    pub fn is_valid(&self) -> bool {
        self.errors.is_empty()
    }
    fn record(&mut self, params: &ReactionValidationParams, issue: ReactionValidationIssue) {
        if !params.silent {
            eprintln!("{}", issue.detail);
        }
        match issue.severity {
            ReactionValidationSeverity::Warning => self.warnings.push(issue),
            ReactionValidationSeverity::Error => self.errors.push(issue),
        }
    }
}
fn issue(
    kind: ReactionValidationIssueKind,
    severity: ReactionValidationSeverity,
    role: Option<ReactionRole>,
    template: Option<usize>,
    atom: Option<AtomId>,
    map: Option<i32>,
    detail: String,
) -> ReactionValidationIssue {
    ReactionValidationIssue {
        kind,
        severity,
        role,
        template,
        atom,
        map,
        maps: Vec::new(),
        detail,
    }
}

pub(crate) fn atom_map(
    atom: &QueryAtom,
    role: ReactionRole,
    template: usize,
) -> Result<Option<i32>, ReactionValidationError> {
    // Canonical typed slots are authoritative. If the slot is absent, retain
    // the source generic-property read, including bad_any_cast and overflow.
    if let Some(value) = atom.atom_map() {
        return i32::try_from(value)
            .map(Some)
            .map_err(|_| ReactionValidationError::MapOverflow {
                role,
                template,
                atom: atom.id(),
                property: "molAtomMapNumber",
                value,
            });
    }
    atom.prop("molAtomMapNumber")
        .map(|value| {
            cosmolkit_core::property_value_to_int(value).map_err(|source| {
                ReactionValidationError::Property {
                    role,
                    template,
                    atom: atom.id(),
                    property: "molAtomMapNumber",
                    source,
                }
            })
        })
        .transpose()
}

fn unnegated(mut node: &QueryNode<AtomQueryPredicate>) -> &QueryNode<AtomQueryPredicate> {
    // Not encodes the source query's negation flag, not an extra query child.
    while let QueryNode::Not(child) = node {
        node = child;
    }
    node
}

fn validate_prepared(
    reaction: &mut Reaction,
    params: &ReactionValidationParams,
) -> Result<ReactionValidationReport, ReactionValidationError> {
    // RDKit❗✔️: bool ChemicalReaction::validate(unsigned int &numWarnings,
    // RDKit❗✔️:                                 unsigned int &numErrors, bool silent) const {
    // RDKit❗✔️:   bool res = true;
    // RDKit❗✔️:   numWarnings = 0;
    // RDKit❗✔️:   numErrors = 0;
    // RDKit❗✔️:
    // RDKit❗✔️:   if (!this->getNumReactantTemplates()) {
    // RDKit❗✔️:     if (!silent) {
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog) << "reaction has no reactants\n";
    // RDKit❗✔️:     }
    // RDKit❗✔️:     numErrors++;
    // RDKit❗✔️:     res = false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (!this->getNumProductTemplates()) {
    // RDKit❗✔️:     if (!silent) {
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog) << "reaction has no products\n";
    // RDKit❗✔️:     }
    // RDKit❗✔️:     numErrors++;
    // RDKit❗✔️:     res = false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<int> mapNumbersSeen;
    // RDKit❗✔️:   std::map<int, const Atom *> reactingAtoms;
    // RDKit❗✔️:   unsigned int molIdx = 0;
    // RDKit❗✔️:   for (auto molIter = this->beginReactantTemplates();
    // RDKit❗✔️:        molIter != this->endReactantTemplates(); ++molIter) {
    // RDKit❗✔️:     bool thisMolMapped = false;
    // RDKit❗✔️:     for (ROMol::AtomIterator atomIt = (*molIter)->beginAtoms();
    // RDKit❗✔️:          atomIt != (*molIter)->endAtoms(); ++atomIt) {
    // RDKit❗✔️:       int mapNum;
    // RDKit❗✔️:       if ((*atomIt)->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗✔️:                                       mapNum)) {
    // RDKit❗✔️:         thisMolMapped = true;
    // RDKit❗✔️:         if (std::find(mapNumbersSeen.begin(), mapNumbersSeen.end(), mapNum) !=
    // RDKit❗✔️:             mapNumbersSeen.end()) {
    // RDKit❗✔️:           if (!silent) {
    // RDKit❗✔️:             BOOST_LOG(rdErrorLog) << "reactant atom-mapping number " << mapNum
    // RDKit❗✔️:                                   << " found multiple times.\n";
    // RDKit❗✔️:           }
    // RDKit❗✔️:           numErrors++;
    // RDKit❗✔️:           res = false;
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           mapNumbersSeen.push_back(mapNum);
    // RDKit❗✔️:           reactingAtoms[mapNum] = *atomIt;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (!thisMolMapped) {
    // RDKit❗✔️:       if (!silent) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "reactant " << molIdx << " has no mapped atoms.\n";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       numWarnings++;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     molIdx++;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<int> productNumbersSeen;
    // RDKit❗✔️:   molIdx = 0;
    // RDKit❗✔️:   for (auto molIter = this->beginProductTemplates();
    // RDKit❗✔️:        molIter != this->endProductTemplates(); ++molIter) {
    // RDKit❗✔️:     // clear out some possible cached properties to prevent
    // RDKit❗✔️:     // misleading warnings
    // RDKit❗✔️:     for (ROMol::AtomIterator atomIt = (*molIter)->beginAtoms();
    // RDKit❗✔️:          atomIt != (*molIter)->endAtoms(); ++atomIt) {
    // RDKit❗✔️:       if ((*atomIt)->hasProp(common_properties::_QueryFormalCharge)) {
    // RDKit❗✔️:         (*atomIt)->clearProp(common_properties::_QueryFormalCharge);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if ((*atomIt)->hasProp(common_properties::_QueryHCount)) {
    // RDKit❗✔️:         (*atomIt)->clearProp(common_properties::_QueryHCount);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if ((*atomIt)->hasProp(common_properties::_QueryMass)) {
    // RDKit❗✔️:         (*atomIt)->clearProp(common_properties::_QueryMass);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if ((*atomIt)->hasProp(common_properties::_QueryIsotope)) {
    // RDKit❗✔️:         (*atomIt)->clearProp(common_properties::_QueryIsotope);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     bool thisMolMapped = false;
    // RDKit❗✔️:     for (ROMol::AtomIterator atomIt = (*molIter)->beginAtoms();
    // RDKit❗✔️:          atomIt != (*molIter)->endAtoms(); ++atomIt) {
    // RDKit❗✔️:       int mapNum;
    // RDKit❗✔️:       if ((*atomIt)->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗✔️:                                       mapNum)) {
    // RDKit❗✔️:         thisMolMapped = true;
    // RDKit❗✔️:         bool seenAlready =
    // RDKit❗✔️:             std::find(productNumbersSeen.begin(), productNumbersSeen.end(),
    // RDKit❗✔️:                       mapNum) != productNumbersSeen.end();
    // RDKit❗✔️:         if (seenAlready) {
    // RDKit❗✔️:           if (!silent) {
    // RDKit❗✔️:             BOOST_LOG(rdWarningLog) << "product atom-mapping number " << mapNum
    // RDKit❗✔️:                                     << " found multiple times.\n";
    // RDKit❗✔️:           }
    // RDKit❗✔️:           numWarnings++;
    // RDKit❗✔️:           // ------------
    // RDKit❗✔️:           //   Always check to see if the atoms connectivity changes independent
    // RDKit❗✔️:           //   if it is mapped multiple times
    // RDKit❗✔️:           // ------------
    // RDKit❗✔️:           const Atom *rAtom = reactingAtoms[mapNum];
    // RDKit❗✔️:           CHECK_INVARIANT(rAtom, "missing atom");
    // RDKit❗✔️:           if (rAtom->getDegree() != (*atomIt)->getDegree()) {
    // RDKit❗✔️:             (*atomIt)->setProp(common_properties::_ReactionDegreeChanged, 1);
    // RDKit❗✔️:           }
    // RDKit❗✔️:
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           productNumbersSeen.push_back(mapNum);
    // RDKit❗✔️:         }
    // RDKit❗✔️:         auto ivIt =
    // RDKit❗✔️:             std::find(mapNumbersSeen.begin(), mapNumbersSeen.end(), mapNum);
    // RDKit❗✔️:         if (ivIt == mapNumbersSeen.end()) {
    // RDKit❗✔️:           if (!seenAlready) {
    // RDKit❗✔️:             if (!silent) {
    // RDKit❗✔️:               BOOST_LOG(rdWarningLog) << "product atom-mapping number "
    // RDKit❗✔️:                                       << mapNum << " not found in reactants.\n";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             numWarnings++;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           mapNumbersSeen.erase(ivIt);
    // RDKit❗✔️:
    // RDKit❗✔️:           // ------------
    // RDKit❗✔️:           //   The atom is mapped, check to see if its connectivity changes
    // RDKit❗✔️:           // ------------
    // RDKit❗✔️:           const Atom *rAtom = reactingAtoms[mapNum];
    // RDKit❗✔️:           CHECK_INVARIANT(rAtom, "missing atom");
    // RDKit❗✔️:           if (rAtom->getDegree() != (*atomIt)->getDegree()) {
    // RDKit❗✔️:             (*atomIt)->setProp(common_properties::_ReactionDegreeChanged, 1);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       // ------------
    // RDKit❗✔️:       //    Deal with queries
    // RDKit❗✔️:       // ------------
    // RDKit❗✔️:       if ((*atomIt)->hasQuery()) {
    // RDKit❗✔️:         std::list<const Atom::QUERYATOM_QUERY *> queries;
    // RDKit❗✔️:         queries.push_back((*atomIt)->getQuery());
    // RDKit❗✔️:         while (!queries.empty()) {
    // RDKit❗✔️:           const Atom::QUERYATOM_QUERY *query = queries.front();
    // RDKit❗✔️:           queries.pop_front();
    // RDKit❗✔️:           for (auto qIter = query->beginChildren();
    // RDKit❗✔️:                qIter != query->endChildren(); ++qIter) {
    // RDKit❗✔️:             queries.push_back((*qIter).get());
    // RDKit❗✔️:           }
    // RDKit❗✔️:           if (query->getDescription() == "AtomFormalCharge" ||
    // RDKit❗✔️:               query->getDescription() == "AtomNegativeFormalCharge") {
    // RDKit❗✔️:             int qval;
    // RDKit❗✔️:             int neg =
    // RDKit❗✔️:                 query->getDescription() == "AtomNegativeFormalCharge" ? -1 : 1;
    // RDKit❗✔️:             if ((*atomIt)->getPropIfPresent(
    // RDKit❗✔️:                     common_properties::_QueryFormalCharge, qval) &&
    // RDKit❗✔️:                 (neg * qval) !=
    // RDKit❗✔️:                     static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
    // RDKit❗✔️:               if (!silent) {
    // RDKit❗✔️:                 BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:                     << "atom " << (*atomIt)->getIdx() << " in product "
    // RDKit❗✔️:                     << molIdx << " has multiple charge specifications.\n";
    // RDKit❗✔️:               }
    // RDKit❗✔️:               numWarnings++;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               int neg = query->getDescription() == "AtomNegativeFormalCharge"
    // RDKit❗✔️:                             ? -1
    // RDKit❗✔️:                             : 1;
    // RDKit❗✔️:               (*atomIt)->setProp(
    // RDKit❗✔️:                   common_properties::_QueryFormalCharge,
    // RDKit❗✔️:                   neg *
    // RDKit❗✔️:                       static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal());
    // RDKit❗✔️:             }
    // RDKit❗✔️:           } else if (query->getDescription() == "AtomHCount") {
    // RDKit❗✔️:             int qval;
    // RDKit❗✔️:             if ((*atomIt)->getPropIfPresent(common_properties::_QueryHCount,
    // RDKit❗✔️:                                             qval) &&
    // RDKit❗✔️:                 qval !=
    // RDKit❗✔️:                     static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
    // RDKit❗✔️:               if (!silent) {
    // RDKit❗✔️:                 BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:                     << "atom " << (*atomIt)->getIdx() << " in product "
    // RDKit❗✔️:                     << molIdx << " has multiple H count specifications.\n";
    // RDKit❗✔️:               }
    // RDKit❗✔️:               numWarnings++;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               (*atomIt)->setProp(
    // RDKit❗✔️:                   common_properties::_QueryHCount,
    // RDKit❗✔️:                   static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal());
    // RDKit❗✔️:             }
    // RDKit❗✔️:           } else if (query->getDescription() == "AtomMass") {
    // RDKit❗✔️:             int qval;
    // RDKit❗✔️:             if ((*atomIt)->getPropIfPresent(common_properties::_QueryMass,
    // RDKit❗✔️:                                             qval) &&
    // RDKit❗✔️:                 qval !=
    // RDKit❗✔️:                     static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
    // RDKit❗✔️:               if (!silent) {
    // RDKit❗✔️:                 BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:                     << "atom " << (*atomIt)->getIdx() << " in product "
    // RDKit❗✔️:                     << molIdx << " has multiple mass specifications.\n";
    // RDKit❗✔️:               }
    // RDKit❗✔️:               numWarnings++;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               (*atomIt)->setProp(
    // RDKit❗✔️:                   common_properties::_QueryMass,
    // RDKit❗✔️:                   static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal() /
    // RDKit❗✔️:                       massIntegerConversionFactor);
    // RDKit❗✔️:             }
    // RDKit❗✔️:           } else if (query->getDescription() == "AtomIsotope") {
    // RDKit❗✔️:             int qval;
    // RDKit❗✔️:             if ((*atomIt)->getPropIfPresent(common_properties::_QueryIsotope,
    // RDKit❗✔️:                                             qval) &&
    // RDKit❗✔️:                 qval !=
    // RDKit❗✔️:                     static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
    // RDKit❗✔️:               if (!silent) {
    // RDKit❗✔️:                 BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:                     << "atom " << (*atomIt)->getIdx() << " in product "
    // RDKit❗✔️:                     << molIdx << " has multiple isotope specifications.\n";
    // RDKit❗✔️:               }
    // RDKit❗✔️:               numWarnings++;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               (*atomIt)->setProp(
    // RDKit❗✔️:                   common_properties::_QueryIsotope,
    // RDKit❗✔️:                   static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal());
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (!thisMolMapped) {
    // RDKit❗✔️:       if (!silent) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "product " << molIdx << " has no mapped atoms.\n";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       numWarnings++;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     molIdx++;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!mapNumbersSeen.empty()) {
    // RDKit❗✔️:     if (!silent) {
    // RDKit❗✔️:       std::ostringstream ostr;
    // RDKit❗✔️:       ostr
    // RDKit❗✔️:           << "mapped atoms in the reactants were not mapped in the products.\n";
    // RDKit❗✔️:       ostr << "  unmapped numbers are: ";
    // RDKit❗✔️:       for (std::vector<int>::const_iterator ivIt = mapNumbersSeen.begin();
    // RDKit❗✔️:            ivIt != mapNumbersSeen.end(); ++ivIt) {
    // RDKit❗✔️:         ostr << *ivIt << " ";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       ostr << "\n";
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << ostr.str();
    // RDKit❗✔️:     }
    // RDKit❗✔️:     numWarnings++;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    use ReactionRole::{Product, Reactant};
    use ReactionValidationIssueKind as K;
    use ReactionValidationSeverity::{Error, Warning};
    let mut report = ReactionValidationReport::default();
    if reaction.reactants.is_empty() {
        report.record(
            params,
            issue(
                K::MissingReactants,
                Error,
                Some(Reactant),
                None,
                None,
                None,
                "reaction has no reactants".into(),
            ),
        );
    }
    if reaction.products.is_empty() {
        report.record(
            params,
            issue(
                K::MissingProducts,
                Error,
                Some(Product),
                None,
                None,
                None,
                "reaction has no products".into(),
            ),
        );
    }
    let mut maps_seen = Vec::new();
    let mut reacting = BTreeMap::new();
    for (template, graph) in reaction.reactants.iter().enumerate() {
        let mut mapped = false;
        for atom in graph.atoms() {
            if let Some(map) = atom_map(atom, Reactant, template)? {
                mapped = true;
                if maps_seen.contains(&map) {
                    report.record(
                        params,
                        issue(
                            K::DuplicateReactantMap,
                            Error,
                            Some(Reactant),
                            Some(template),
                            Some(atom.id()),
                            Some(map),
                            format!("reactant atom-mapping number {map} found multiple times."),
                        ),
                    );
                } else {
                    maps_seen.push(map);
                    reacting.insert(map, graph.adjacency()[atom.id().index()].len());
                }
            }
        }
        if !mapped {
            report.record(
                params,
                issue(
                    K::UnmappedReactant,
                    Warning,
                    Some(Reactant),
                    Some(template),
                    None,
                    None,
                    format!("reactant {template} has no mapped atoms."),
                ),
            );
        }
    }
    let mut product_seen = Vec::new();
    for (template, graph) in reaction.products.iter_mut().enumerate() {
        for atom in graph.atoms_mut() {
            for key in [
                "_QueryFormalCharge",
                "_QueryHCount",
                "_QueryMass",
                "_QueryIsotope",
            ] {
                atom.clear_prop(key);
            }
        }
        let mut mapped = false;
        for row in 0..graph.num_atoms() {
            let id = AtomId::new(row);
            if let Some(map) = atom_map(&graph.atoms()[row], Product, template)? {
                mapped = true;
                let seen_already = product_seen.contains(&map);
                if seen_already {
                    report.record(
                        params,
                        issue(
                            K::DuplicateProductMap,
                            Warning,
                            Some(Product),
                            Some(template),
                            Some(id),
                            Some(map),
                            format!("product atom-mapping number {map} found multiple times."),
                        ),
                    );
                    let degree =
                        reacting
                            .get(&map)
                            .ok_or(ReactionValidationError::MissingReactingAtom {
                                template,
                                atom: id,
                                map,
                            })?;
                    if *degree != graph.adjacency()[row].len() {
                        graph.atoms_mut()[row]
                            .set_prop("_ReactionDegreeChanged", 1i32)
                            .map_err(|source| ReactionValidationError::Annotation {
                                template,
                                atom: id,
                                property: "_ReactionDegreeChanged",
                                source,
                            })?;
                    }
                } else {
                    product_seen.push(map);
                }
                if let Some(pos) = maps_seen.iter().position(|candidate| *candidate == map) {
                    maps_seen.remove(pos);
                    let degree =
                        reacting
                            .get(&map)
                            .ok_or(ReactionValidationError::MissingReactingAtom {
                                template,
                                atom: id,
                                map,
                            })?;
                    if *degree != graph.adjacency()[row].len() {
                        graph.atoms_mut()[row]
                            .set_prop("_ReactionDegreeChanged", 1i32)
                            .map_err(|source| ReactionValidationError::Annotation {
                                template,
                                atom: id,
                                property: "_ReactionDegreeChanged",
                                source,
                            })?;
                    }
                } else if !seen_already {
                    report.record(
                        params,
                        issue(
                            K::MissingProductMapReactant,
                            Warning,
                            Some(Product),
                            Some(template),
                            Some(id),
                            Some(map),
                            format!("product atom-mapping number {map} not found in reactants."),
                        ),
                    );
                }
            }
            if graph.atoms()[row].predicate_is_carrier_derived() {
                continue;
            }
            // Annotations are evaluated breadth first. Keep duplicate/negated
            // constraints and the source first-value retention on conflict.
            let atom = &mut graph.atoms_mut()[row];
            let mut annotation_values: BTreeMap<&'static str, i32> = BTreeMap::new();
            let mut queue = VecDeque::from([atom.predicate()]);
            let mut writes = Vec::new();
            while let Some(node) = queue.pop_front() {
                let node = unnegated(node);
                if let QueryNode::And(children)
                | QueryNode::Or(children)
                | QueryNode::Xor(children) = node
                {
                    queue.extend(children);
                }
                let QueryNode::Predicate(predicate) = node else {
                    continue;
                };
                let (key, value, sign, store, kind, description) =
                    match predicate {
                        AtomQueryPredicate::FormalCharge(value) => (
                            "_QueryFormalCharge",
                            *value,
                            1i32,
                            *value,
                            K::MultipleCharge,
                            "charge",
                        ),
                        AtomQueryPredicate::NegativeFormalCharge(value) => (
                            "_QueryFormalCharge",
                            *value,
                            -1i32,
                            value.checked_neg().ok_or(
                                ReactionValidationError::QueryArithmetic { template, atom: id },
                            )?,
                            K::MultipleCharge,
                            "charge",
                        ),
                        AtomQueryPredicate::HydrogenCount(value) => (
                            "_QueryHCount",
                            *value,
                            1,
                            *value,
                            K::MultipleHydrogenCount,
                            "H count",
                        ),
                        AtomQueryPredicate::Mass(value) => (
                            "_QueryMass",
                            i32::from(*value) * 1000,
                            1,
                            i32::from(*value),
                            K::MultipleMass,
                            "mass",
                        ),
                        AtomQueryPredicate::Isotope(value) => (
                            "_QueryIsotope",
                            *value,
                            1,
                            *value,
                            K::MultipleIsotope,
                            "isotope",
                        ),
                        _ => continue,
                    };
                let conflict = if let Some(previous) = annotation_values.get(key) {
                    previous
                        .checked_mul(sign)
                        .ok_or(ReactionValidationError::QueryArithmetic { template, atom: id })?
                        != value
                } else {
                    false
                };
                if conflict {
                    report.record(params, issue(kind, Warning, Some(Product), Some(template), Some(id), None, format!("atom {row} in product {template} has multiple {description} specifications.")));
                } else {
                    annotation_values.insert(key, store);
                    writes.push((key, store));
                }
            }
            for (property, value) in writes {
                atom.set_prop(property, value).map_err(|source| {
                    ReactionValidationError::Annotation {
                        template,
                        atom: id,
                        property,
                        source,
                    }
                })?;
            }
        }
        if !mapped {
            report.record(
                params,
                issue(
                    K::UnmappedProduct,
                    Warning,
                    Some(Product),
                    Some(template),
                    None,
                    None,
                    format!("product {template} has no mapped atoms."),
                ),
            );
        }
    }
    if !maps_seen.is_empty() {
        let mut detail = String::from(
            "mapped atoms in the reactants were not mapped in the products.\n  unmapped numbers are: ",
        );
        for map in &maps_seen {
            use std::fmt::Write;
            write!(&mut detail, "{map} ").expect("writing to a String is infallible");
        }
        let mut item = issue(
            K::UnmappedReactantMaps,
            Warning,
            Some(Reactant),
            None,
            None,
            None,
            detail,
        );
        item.maps = maps_seen;
        report.record(params, item);
    }
    Ok(report)
}

/// Query validation uses private detached annotations and leaves the caller
/// unchanged under the approved immutable Reaction contract.
#[doc(hidden)]
pub fn validate_reaction(
    reaction: &Reaction,
    params: &ReactionValidationParams,
) -> Result<ReactionValidationReport, ReactionValidationError> {
    let mut prepared = reaction.clone();
    validate_prepared(&mut prepared, params)
}

#[doc(hidden)]
pub fn initialize_reaction(
    reaction: &Reaction,
    params: &ReactionValidationParams,
) -> Result<Reaction, ReactionInitializationError> {
    // RDKit❗✔️: void ChemicalReaction::initReactantMatchers(bool silent) {
    // RDKit❗✔️:   unsigned int nWarnings, nErrors;
    // RDKit❗✔️:   if (!this->validate(nWarnings, nErrors, silent)) {
    // RDKit❗✔️:     BOOST_LOG(rdErrorLog) << "initialization failed\n";
    // RDKit❗✔️:     this->df_needsInit = true;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     this->df_needsInit = false;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️:
    // D1: source mutations belong to the returned value, never the caller.
    let mut prepared = reaction.clone();
    let report = validate_prepared(&mut prepared, params)
        .map_err(|source| ReactionInitializationError::Validation { source })?;
    if !report.is_valid() {
        prepared.needs_init = true;
        eprintln!("initialization failed");
        return Err(ReactionInitializationError::Invalid { report });
    }
    prepared.needs_init = false;
    Ok(prepared)
}
