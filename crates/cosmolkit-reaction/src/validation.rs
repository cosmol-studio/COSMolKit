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
        (self.warnings.len() as u32) as usize
    }
    pub fn num_errors(&self) -> usize {
        (self.errors.len() as u32) as usize
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
    atom: &impl crate::materialize::ReactionAtomPropertyRead,
    role: ReactionRole,
    template: usize,
) -> Result<Option<i32>, ReactionValidationError> {
    // Canonical typed slots are authoritative. If the slot is absent, retain
    // the source generic-property read, including bad_any_cast and overflow.
    if let Some(value) = atom.typed_map() {
        return i32::try_from(value)
            .map(Some)
            .map_err(|_| ReactionValidationError::MapOverflow {
                role,
                template,
                atom: atom.property_id(),
                property: "molAtomMapNumber",
                value,
            });
    }
    atom.property("molAtomMapNumber")
        .map(|value| {
            cosmolkit_core::property_value_to_int(value).map_err(|source| {
                ReactionValidationError::Property {
                    role,
                    template,
                    atom: atom.property_id(),
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
    report: &mut ReactionValidationReport,
) -> Result<bool, ReactionValidationError> {
    // RDKit❗❌: bool ChemicalReaction::validate(unsigned int &numWarnings,
    // RDKit❗❌:                                 unsigned int &numErrors, bool silent) const {
    // RDKit❗❌:   bool res = true;
    // RDKit❗❌:   numWarnings = 0;
    // RDKit❗❌:   numErrors = 0;
    // RDKit❗❌:
    // RDKit❗❌:   if (!this->getNumReactantTemplates()) {
    // RDKit❗❌:     if (!silent) {
    // RDKit❗❌:       BOOST_LOG(rdErrorLog) << "reaction has no reactants\n";
    // RDKit❗❌:     }
    // RDKit❗❌:     numErrors++;
    // RDKit❗❌:     res = false;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (!this->getNumProductTemplates()) {
    // RDKit❗❌:     if (!silent) {
    // RDKit❗❌:       BOOST_LOG(rdErrorLog) << "reaction has no products\n";
    // RDKit❗❌:     }
    // RDKit❗❌:     numErrors++;
    // RDKit❗❌:     res = false;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<int> mapNumbersSeen;
    // RDKit❗❌:   std::map<int, const Atom *> reactingAtoms;
    // RDKit❗❌:   unsigned int molIdx = 0;
    // RDKit❗❌:   for (auto molIter = this->beginReactantTemplates();
    // RDKit❗❌:        molIter != this->endReactantTemplates(); ++molIter) {
    // RDKit❗❌:     bool thisMolMapped = false;
    // RDKit❗❌:     for (ROMol::AtomIterator atomIt = (*molIter)->beginAtoms();
    // RDKit❗❌:          atomIt != (*molIter)->endAtoms(); ++atomIt) {
    // RDKit❗❌:       int mapNum;
    // RDKit❗❌:       if ((*atomIt)->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗❌:                                       mapNum)) {
    // RDKit❗❌:         thisMolMapped = true;
    // RDKit❗❌:         if (std::find(mapNumbersSeen.begin(), mapNumbersSeen.end(), mapNum) !=
    // RDKit❗❌:             mapNumbersSeen.end()) {
    // RDKit❗❌:           if (!silent) {
    // RDKit❗❌:             BOOST_LOG(rdErrorLog) << "reactant atom-mapping number " << mapNum
    // RDKit❗❌:                                   << " found multiple times.\n";
    // RDKit❗❌:           }
    // RDKit❗❌:           numErrors++;
    // RDKit❗❌:           res = false;
    // RDKit❗❌:         } else {
    // RDKit❗❌:           mapNumbersSeen.push_back(mapNum);
    // RDKit❗❌:           reactingAtoms[mapNum] = *atomIt;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!thisMolMapped) {
    // RDKit❗❌:       if (!silent) {
    // RDKit❗❌:         BOOST_LOG(rdWarningLog)
    // RDKit❗❌:             << "reactant " << molIdx << " has no mapped atoms.\n";
    // RDKit❗❌:       }
    // RDKit❗❌:       numWarnings++;
    // RDKit❗❌:     }
    // RDKit❗❌:     molIdx++;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<int> productNumbersSeen;
    // RDKit❗❌:   molIdx = 0;
    // RDKit❗❌:   for (auto molIter = this->beginProductTemplates();
    // RDKit❗❌:        molIter != this->endProductTemplates(); ++molIter) {
    // RDKit❗❌:     // clear out some possible cached properties to prevent
    // RDKit❗❌:     // misleading warnings
    // RDKit❗❌:     for (ROMol::AtomIterator atomIt = (*molIter)->beginAtoms();
    // RDKit❗❌:          atomIt != (*molIter)->endAtoms(); ++atomIt) {
    // RDKit❗❌:       if ((*atomIt)->hasProp(common_properties::_QueryFormalCharge)) {
    // RDKit❗❌:         (*atomIt)->clearProp(common_properties::_QueryFormalCharge);
    // RDKit❗❌:       }
    // RDKit❗❌:       if ((*atomIt)->hasProp(common_properties::_QueryHCount)) {
    // RDKit❗❌:         (*atomIt)->clearProp(common_properties::_QueryHCount);
    // RDKit❗❌:       }
    // RDKit❗❌:       if ((*atomIt)->hasProp(common_properties::_QueryMass)) {
    // RDKit❗❌:         (*atomIt)->clearProp(common_properties::_QueryMass);
    // RDKit❗❌:       }
    // RDKit❗❌:       if ((*atomIt)->hasProp(common_properties::_QueryIsotope)) {
    // RDKit❗❌:         (*atomIt)->clearProp(common_properties::_QueryIsotope);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     bool thisMolMapped = false;
    // RDKit❗❌:     for (ROMol::AtomIterator atomIt = (*molIter)->beginAtoms();
    // RDKit❗❌:          atomIt != (*molIter)->endAtoms(); ++atomIt) {
    // RDKit❗❌:       int mapNum;
    // RDKit❗❌:       if ((*atomIt)->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗❌:                                       mapNum)) {
    // RDKit❗❌:         thisMolMapped = true;
    // RDKit❗❌:         bool seenAlready =
    // RDKit❗❌:             std::find(productNumbersSeen.begin(), productNumbersSeen.end(),
    // RDKit❗❌:                       mapNum) != productNumbersSeen.end();
    // RDKit❗❌:         if (seenAlready) {
    // RDKit❗❌:           if (!silent) {
    // RDKit❗❌:             BOOST_LOG(rdWarningLog) << "product atom-mapping number " << mapNum
    // RDKit❗❌:                                     << " found multiple times.\n";
    // RDKit❗❌:           }
    // RDKit❗❌:           numWarnings++;
    // RDKit❗❌:           // ------------
    // RDKit❗❌:           //   Always check to see if the atoms connectivity changes independent
    // RDKit❗❌:           //   if it is mapped multiple times
    // RDKit❗❌:           // ------------
    // RDKit❗❌:           const Atom *rAtom = reactingAtoms[mapNum];
    // RDKit❗❌:           CHECK_INVARIANT(rAtom, "missing atom");
    // RDKit❗❌:           if (rAtom->getDegree() != (*atomIt)->getDegree()) {
    // RDKit❗❌:             (*atomIt)->setProp(common_properties::_ReactionDegreeChanged, 1);
    // RDKit❗❌:           }
    // RDKit❗❌:
    // RDKit❗❌:         } else {
    // RDKit❗❌:           productNumbersSeen.push_back(mapNum);
    // RDKit❗❌:         }
    // RDKit❗❌:         auto ivIt =
    // RDKit❗❌:             std::find(mapNumbersSeen.begin(), mapNumbersSeen.end(), mapNum);
    // RDKit❗❌:         if (ivIt == mapNumbersSeen.end()) {
    // RDKit❗❌:           if (!seenAlready) {
    // RDKit❗❌:             if (!silent) {
    // RDKit❗❌:               BOOST_LOG(rdWarningLog) << "product atom-mapping number "
    // RDKit❗❌:                                       << mapNum << " not found in reactants.\n";
    // RDKit❗❌:             }
    // RDKit❗❌:             numWarnings++;
    // RDKit❗❌:           }
    // RDKit❗❌:         } else {
    // RDKit❗❌:           mapNumbersSeen.erase(ivIt);
    // RDKit❗❌:
    // RDKit❗❌:           // ------------
    // RDKit❗❌:           //   The atom is mapped, check to see if its connectivity changes
    // RDKit❗❌:           // ------------
    // RDKit❗❌:           const Atom *rAtom = reactingAtoms[mapNum];
    // RDKit❗❌:           CHECK_INVARIANT(rAtom, "missing atom");
    // RDKit❗❌:           if (rAtom->getDegree() != (*atomIt)->getDegree()) {
    // RDKit❗❌:             (*atomIt)->setProp(common_properties::_ReactionDegreeChanged, 1);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       // ------------
    // RDKit❗❌:       //    Deal with queries
    // RDKit❗❌:       // ------------
    // RDKit❗❌:       if ((*atomIt)->hasQuery()) {
    // RDKit❗❌:         std::list<const Atom::QUERYATOM_QUERY *> queries;
    // RDKit❗❌:         queries.push_back((*atomIt)->getQuery());
    // RDKit❗❌:         while (!queries.empty()) {
    // RDKit❗❌:           const Atom::QUERYATOM_QUERY *query = queries.front();
    // RDKit❗❌:           queries.pop_front();
    // RDKit❗❌:           for (auto qIter = query->beginChildren();
    // RDKit❗❌:                qIter != query->endChildren(); ++qIter) {
    // RDKit❗❌:             queries.push_back((*qIter).get());
    // RDKit❗❌:           }
    // RDKit❗❌:           if (query->getDescription() == "AtomFormalCharge" ||
    // RDKit❗❌:               query->getDescription() == "AtomNegativeFormalCharge") {
    // RDKit❗❌:             int qval;
    // RDKit❗❌:             int neg =
    // RDKit❗❌:                 query->getDescription() == "AtomNegativeFormalCharge" ? -1 : 1;
    // RDKit❗❌:             if ((*atomIt)->getPropIfPresent(
    // RDKit❗❌:                     common_properties::_QueryFormalCharge, qval) &&
    // RDKit❗❌:                 (neg * qval) !=
    // RDKit❗❌:                     static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
    // RDKit❗❌:               if (!silent) {
    // RDKit❗❌:                 BOOST_LOG(rdWarningLog)
    // RDKit❗❌:                     << "atom " << (*atomIt)->getIdx() << " in product "
    // RDKit❗❌:                     << molIdx << " has multiple charge specifications.\n";
    // RDKit❗❌:               }
    // RDKit❗❌:               numWarnings++;
    // RDKit❗❌:             } else {
    // RDKit❗❌:               int neg = query->getDescription() == "AtomNegativeFormalCharge"
    // RDKit❗❌:                             ? -1
    // RDKit❗❌:                             : 1;
    // RDKit❗❌:               (*atomIt)->setProp(
    // RDKit❗❌:                   common_properties::_QueryFormalCharge,
    // RDKit❗❌:                   neg *
    // RDKit❗❌:                       static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal());
    // RDKit❗❌:             }
    // RDKit❗❌:           } else if (query->getDescription() == "AtomHCount") {
    // RDKit❗❌:             int qval;
    // RDKit❗❌:             if ((*atomIt)->getPropIfPresent(common_properties::_QueryHCount,
    // RDKit❗❌:                                             qval) &&
    // RDKit❗❌:                 qval !=
    // RDKit❗❌:                     static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
    // RDKit❗❌:               if (!silent) {
    // RDKit❗❌:                 BOOST_LOG(rdWarningLog)
    // RDKit❗❌:                     << "atom " << (*atomIt)->getIdx() << " in product "
    // RDKit❗❌:                     << molIdx << " has multiple H count specifications.\n";
    // RDKit❗❌:               }
    // RDKit❗❌:               numWarnings++;
    // RDKit❗❌:             } else {
    // RDKit❗❌:               (*atomIt)->setProp(
    // RDKit❗❌:                   common_properties::_QueryHCount,
    // RDKit❗❌:                   static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal());
    // RDKit❗❌:             }
    // RDKit❗❌:           } else if (query->getDescription() == "AtomMass") {
    // RDKit❗❌:             int qval;
    // RDKit❗❌:             if ((*atomIt)->getPropIfPresent(common_properties::_QueryMass,
    // RDKit❗❌:                                             qval) &&
    // RDKit❗❌:                 qval !=
    // RDKit❗❌:                     static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
    // RDKit❗❌:               if (!silent) {
    // RDKit❗❌:                 BOOST_LOG(rdWarningLog)
    // RDKit❗❌:                     << "atom " << (*atomIt)->getIdx() << " in product "
    // RDKit❗❌:                     << molIdx << " has multiple mass specifications.\n";
    // RDKit❗❌:               }
    // RDKit❗❌:               numWarnings++;
    // RDKit❗❌:             } else {
    // RDKit❗❌:               (*atomIt)->setProp(
    // RDKit❗❌:                   common_properties::_QueryMass,
    // RDKit❗❌:                   static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal() /
    // RDKit❗❌:                       massIntegerConversionFactor);
    // RDKit❗❌:             }
    // RDKit❗❌:           } else if (query->getDescription() == "AtomIsotope") {
    // RDKit❗❌:             int qval;
    // RDKit❗❌:             if ((*atomIt)->getPropIfPresent(common_properties::_QueryIsotope,
    // RDKit❗❌:                                             qval) &&
    // RDKit❗❌:                 qval !=
    // RDKit❗❌:                     static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
    // RDKit❗❌:               if (!silent) {
    // RDKit❗❌:                 BOOST_LOG(rdWarningLog)
    // RDKit❗❌:                     << "atom " << (*atomIt)->getIdx() << " in product "
    // RDKit❗❌:                     << molIdx << " has multiple isotope specifications.\n";
    // RDKit❗❌:               }
    // RDKit❗❌:               numWarnings++;
    // RDKit❗❌:             } else {
    // RDKit❗❌:               (*atomIt)->setProp(
    // RDKit❗❌:                   common_properties::_QueryIsotope,
    // RDKit❗❌:                   static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal());
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!thisMolMapped) {
    // RDKit❗❌:       if (!silent) {
    // RDKit❗❌:         BOOST_LOG(rdWarningLog)
    // RDKit❗❌:             << "product " << molIdx << " has no mapped atoms.\n";
    // RDKit❗❌:       }
    // RDKit❗❌:       numWarnings++;
    // RDKit❗❌:     }
    // RDKit❗❌:     molIdx++;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!mapNumbersSeen.empty()) {
    // RDKit❗❌:     if (!silent) {
    // RDKit❗❌:       std::ostringstream ostr;
    // RDKit❗❌:       ostr
    // RDKit❗❌:           << "mapped atoms in the reactants were not mapped in the products.\n";
    // RDKit❗❌:       ostr << "  unmapped numbers are: ";
    // RDKit❗❌:       for (std::vector<int>::const_iterator ivIt = mapNumbersSeen.begin();
    // RDKit❗❌:            ivIt != mapNumbersSeen.end(); ++ivIt) {
    // RDKit❗❌:         ostr << *ivIt << " ";
    // RDKit❗❌:       }
    // RDKit❗❌:       ostr << "\n";
    // RDKit❗❌:       BOOST_LOG(rdWarningLog) << ostr.str();
    // RDKit❗❌:     }
    // RDKit❗❌:     numWarnings++;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    use ReactionRole::{Product, Reactant};
    use ReactionValidationIssueKind as K;
    use ReactionValidationSeverity::{Error, Warning};
    *report = ReactionValidationReport::default();
    if reaction.num_reactant_templates() == 0 {
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
    if reaction.num_product_templates() == 0 {
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
        for (row, atom) in graph.atoms().iter().enumerate() {
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
                    reacting.insert(map, (template, row));
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
                if crate::management::atom_has_property_source(atom, key.as_bytes()) {
                    atom.clear_prop(key)
                        .map_err(|source| ReactionValidationError::Annotation {
                            template,
                            atom: atom.id(),
                            property: key,
                            source,
                        })?;
                }
            }
        }
        let mut mapped = false;
        for row in 0..graph.num_atoms() {
            let id = graph.atoms()[row].id();
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
                    let reactant = &reaction.reactants[degree.0];
                    let react_atom = reactant.atoms()[degree.1].id();
                    let react_degree = reactant.try_atom_degree(react_atom).ok_or(
                        ReactionValidationError::AtomBounds {
                            role: Reactant,
                            template: degree.0,
                            atom: react_atom,
                            atom_count: reactant.num_atoms(),
                        },
                    )?;
                    let prod_degree =
                        graph
                            .try_atom_degree(id)
                            .ok_or(ReactionValidationError::AtomBounds {
                                role: Product,
                                template,
                                atom: id,
                                atom_count: graph.num_atoms(),
                            })?;
                    if react_degree != prod_degree {
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
                    let reactant = &reaction.reactants[degree.0];
                    let react_atom = reactant.atoms()[degree.1].id();
                    let react_degree = reactant.try_atom_degree(react_atom).ok_or(
                        ReactionValidationError::AtomBounds {
                            role: Reactant,
                            template: degree.0,
                            atom: react_atom,
                            atom_count: reactant.num_atoms(),
                        },
                    )?;
                    let prod_degree =
                        graph
                            .try_atom_degree(id)
                            .ok_or(ReactionValidationError::AtomBounds {
                                role: Product,
                                template,
                                atom: id,
                                atom_count: graph.num_atoms(),
                            })?;
                    if react_degree != prod_degree {
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
            let (query, properties) = graph.atoms_mut()[row].predicate_and_properties_mut();
            let mut queue = VecDeque::from([query]);
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
                let (property, value, sign, divisor, kind, description) = match predicate {
                    AtomQueryPredicate::FormalCharge(value) => (
                        "_QueryFormalCharge",
                        *value,
                        1i32,
                        1,
                        K::MultipleCharge,
                        "charge",
                    ),
                    AtomQueryPredicate::NegativeFormalCharge(value) => (
                        "_QueryFormalCharge",
                        *value,
                        -1i32,
                        1,
                        K::MultipleCharge,
                        "charge",
                    ),
                    AtomQueryPredicate::HydrogenCount(value) => (
                        "_QueryHCount",
                        *value,
                        1,
                        1,
                        K::MultipleHydrogenCount,
                        "H count",
                    ),
                    AtomQueryPredicate::Mass(value) => (
                        "_QueryMass",
                        i32::from(*value) * 1000,
                        1,
                        1000,
                        K::MultipleMass,
                        "mass",
                    ),
                    AtomQueryPredicate::Isotope(value) => {
                        ("_QueryIsotope", *value, 1, 1, K::MultipleIsotope, "isotope")
                    }
                    _ => continue,
                };
                let previous = properties
                    .get(property.as_bytes())
                    .map(|previous| {
                        cosmolkit_core::property_value_to_int(previous).map_err(|source| {
                            ReactionValidationError::Property {
                                role: Product,
                                template,
                                atom: id,
                                property,
                                source,
                            }
                        })
                    })
                    .transpose()?;
                let conflict = if let Some(previous) = previous {
                    previous
                        .checked_mul(sign)
                        .ok_or(ReactionValidationError::QueryArithmetic { template, atom: id })?
                        != value
                } else {
                    false
                };
                if conflict {
                    report.record(params, issue(kind, Warning, Some(Product), Some(template), Some(id), None,
                        format!("atom {id} in product {template} has multiple {description} specifications.")));
                } else {
                    // Native writes immediately, and only this branch evaluates
                    // neg*getVal. Keep the mass comparison against the unscaled
                    // query value even though the stored annotation is divided.
                    let store = value
                        .checked_mul(sign)
                        .ok_or(ReactionValidationError::QueryArithmetic { template, atom: id })?
                        / divisor;
                    properties
                        .set(property.into(), cosmolkit_model::PropertyValue::Int(store))
                        .map_err(|source| {
                            let source = match source {
                                cosmolkit_model::PropertyStoreError::EmptyKey => {
                                    cosmolkit_model::AtomPropertyError::EmptyKey
                                }
                                cosmolkit_model::PropertyStoreError::ComputedListKind(source) => {
                                    cosmolkit_model::AtomPropertyError::ComputedListKind(source)
                                }
                            };
                            ReactionValidationError::Annotation {
                                template,
                                atom: id,
                                property,
                                source,
                            }
                        })?;
                }
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
    Ok(report.errors.is_empty())
}

/// Query validation uses private detached annotations and leaves the caller
/// unchanged under the approved immutable Reaction contract.
#[doc(hidden)]
pub fn validate_reaction(
    reaction: &Reaction,
    params: &ReactionValidationParams,
) -> Result<ReactionValidationReport, ReactionValidationError> {
    let mut prepared = reaction.clone();
    let mut report = ReactionValidationReport::default();
    validate_prepared(&mut prepared, params, &mut report)?;
    Ok(report)
}

pub(crate) fn init_reactant_matchers_source(
    reaction: &mut Reaction,
    params: &ReactionValidationParams,
    report: &mut ReactionValidationReport,
) -> Result<(), ReactionValidationError> {
    // RDKit❗✔️: void ChemicalReaction::initReactantMatchers(bool silent) {
    // RDKit❗✔️:   unsigned int nWarnings, nErrors;
    // RDKit❗✔️:   if (!this->validate(nWarnings, nErrors, silent)) {
    // RDKit❗✔️:     BOOST_LOG(rdErrorLog) << "initialization failed\n";
    // RDKit❗✔️:     this->df_needsInit = true;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     this->df_needsInit = false;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // The source returns normally after a false validation result. Validation
    // exceptions leave the prior initialization flag and mutation prefix intact.
    if !validate_prepared(reaction, params, report)? {
        eprintln!("initialization failed");
        reaction.needs_init = true;
    } else {
        reaction.needs_init = false;
    }
    Ok(())
}

#[doc(hidden)]
pub fn initialize_reaction(
    reaction: &Reaction,
    params: &ReactionValidationParams,
) -> Result<Reaction, ReactionInitializationError> {
    // The approved immutable projection adapts the source's normal false
    // validation return to the existing structured Invalid report.
    let mut prepared = reaction.clone();
    let mut report = ReactionValidationReport::default();
    init_reactant_matchers_source(&mut prepared, params, &mut report)
        .map_err(|source| ReactionInitializationError::Validation { source })?;
    if !report.is_valid() {
        return Err(ReactionInitializationError::Invalid { report });
    }
    Ok(prepared)
}

#[cfg(test)]
mod complete_reaction_validate_source_tests {
    use super::*;
    use AtomQueryPredicate as P;
    use QueryNode as Q;
    use cosmolkit_model::{Atom, AtomSpec, BondId, BondSpec, PropertyValue, QueryBond, QueryGraph};
    use cosmolkit_types::{BondOrder, Element};

    fn atom(row: usize, map: Option<i32>, query: Q<P>) -> QueryAtom {
        let mut atom = QueryAtom::new(AtomId::new(row), AtomSpec::new(Element::C));
        atom.set_predicate(query);
        if let Some(map) = map {
            atom.set_prop("molAtomMapNumber", PropertyValue::Int(map))
                .unwrap();
        }
        atom
    }
    fn graph(atoms: Vec<QueryAtom>, edges: &[(usize, usize)]) -> QueryGraph {
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                QueryBond::new(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect();
        QueryGraph::from_parts(atoms, bonds, [], vec![], vec![], vec![]).unwrap()
    }
    fn single(map: Option<i32>, query: Q<P>) -> QueryGraph {
        graph(vec![atom(0, map, query)], &[])
    }
    fn paired(product: QueryGraph) -> Reaction {
        let mut reaction = Reaction::new();
        reaction.reactants = vec![single(Some(1), Q::Predicate(P::Any))];
        reaction.products = vec![product];
        reaction
    }
    fn validate(
        reaction: &mut Reaction,
        report: &mut ReactionValidationReport,
    ) -> Result<bool, ReactionValidationError> {
        validate_prepared(reaction, &ReactionValidationParams { silent: true }, report)
    }
    fn props<'a>(reaction: &'a Reaction, key: &str) -> Option<&'a PropertyValue> {
        reaction.products[0].atoms()[0].prop(key)
    }

    #[test]
    fn resets_existing_outputs_and_preserves_source_missing_role_error_order() {
        let mut reaction = Reaction::new();
        let stale = issue(
            ReactionValidationIssueKind::DuplicateProductMap,
            ReactionValidationSeverity::Warning,
            None,
            None,
            None,
            Some(99),
            "stale".into(),
        );
        let mut report = ReactionValidationReport {
            warnings: vec![stale],
            errors: vec![],
        };
        assert!(!validate(&mut reaction, &mut report).unwrap());
        assert_eq!(report.num_warnings(), 0);
        assert_eq!(report.num_errors(), 2);
        assert_eq!(
            report.errors.iter().map(|i| i.kind).collect::<Vec<_>>(),
            [
                ReactionValidationIssueKind::MissingReactants,
                ReactionValidationIssueKind::MissingProducts
            ]
        );
    }

    #[test]
    fn remaining_map_numbers_keep_first_occurrence_and_vector_erase_order() {
        let mut reaction = Reaction::new();
        reaction.reactants = [9, -2, 0, 9]
            .into_iter()
            .map(|m| single(Some(m), Q::Predicate(P::Any)))
            .collect();
        reaction.products = vec![single(Some(-2), Q::Predicate(P::Any))];
        let mut report = ReactionValidationReport::default();
        assert!(!validate(&mut reaction, &mut report).unwrap());
        assert_eq!(report.errors.len(), 1);
        assert_eq!(report.warnings.len(), 1);
        assert_eq!(report.warnings[0].maps, [9, 0]);
    }

    #[test]
    fn duplicate_product_missing_reactant_retains_output_warning_prefix_on_invariant_error() {
        let mut reaction = paired(single(Some(5), Q::Predicate(P::Any)));
        reaction
            .products
            .push(single(Some(5), Q::Predicate(P::Any)));
        let mut report = ReactionValidationReport::default();
        assert!(matches!(
            validate(&mut reaction, &mut report),
            Err(ReactionValidationError::MissingReactingAtom {
                template: 1,
                map: 5,
                ..
            })
        ));
        assert_eq!(
            report.warnings.iter().map(|i| i.kind).collect::<Vec<_>>(),
            [
                ReactionValidationIssueKind::MissingProductMapReactant,
                ReactionValidationIssueKind::DuplicateProductMap
            ]
        );
        assert_eq!(report.errors.len(), 0);
    }

    #[test]
    fn duplicate_reactant_uses_first_atom_object_for_later_degree_comparison() {
        let mut reaction = paired(single(Some(1), Q::Predicate(P::Any)));
        reaction.reactants.push(graph(
            vec![
                atom(0, Some(1), Q::Predicate(P::Any)),
                atom(1, None, Q::Predicate(P::Any)),
            ],
            &[(0, 1)],
        ));
        let mut report = ReactionValidationReport::default();
        assert!(!validate(&mut reaction, &mut report).unwrap());
        assert_eq!(report.errors.len(), 1);
        assert_eq!(report.warnings.len(), 0);
        assert_eq!(props(&reaction, "_ReactionDegreeChanged"), None);
    }

    #[test]
    fn degree_change_is_written_for_every_duplicate_and_existing_equal_degree_flag_is_kept() {
        let mut reaction = paired(single(Some(1), Q::Predicate(P::Any)));
        reaction.products[0].atoms_mut()[0]
            .set_prop("_ReactionDegreeChanged", 99i32)
            .unwrap();
        let mut report = ReactionValidationReport::default();
        assert!(validate(&mut reaction, &mut report).unwrap());
        assert_eq!(
            props(&reaction, "_ReactionDegreeChanged"),
            Some(&PropertyValue::Int(99))
        );
        reaction.reactants[0] = graph(
            vec![
                atom(0, Some(1), Q::Predicate(P::Any)),
                atom(1, None, Q::Predicate(P::Any)),
            ],
            &[(0, 1)],
        );
        reaction
            .products
            .push(single(Some(1), Q::Predicate(P::Any)));
        assert!(validate(&mut reaction, &mut report).unwrap());
        assert_eq!(
            props(&reaction, "_ReactionDegreeChanged"),
            Some(&PropertyValue::Int(1))
        );
        assert_eq!(
            reaction.products[1].atoms()[0].prop("_ReactionDegreeChanged"),
            Some(&PropertyValue::Int(1))
        );
        assert_eq!(
            report.warnings[0].kind,
            ReactionValidationIssueKind::DuplicateProductMap
        );
    }

    #[test]
    fn breadth_first_annotations_keep_negation_flags_and_native_repeated_mass_warning() {
        let query = Q::And(vec![
            Q::Or(vec![
                Q::Predicate(P::FormalCharge(2)),
                Q::Predicate(P::HydrogenCount(3)),
            ]),
            Q::Predicate(P::FormalCharge(1)),
            Q::Predicate(P::Mass(12)),
            Q::Predicate(P::Mass(12)),
            Q::Predicate(P::Isotope(7)),
            Q::Predicate(P::Isotope(8)),
            Q::Not(Box::new(Q::Predicate(P::FormalCharge(1)))),
        ]);
        let mut reaction = paired(single(Some(1), query));
        let mut report = ReactionValidationReport::default();
        assert!(validate(&mut reaction, &mut report).unwrap());
        assert_eq!(
            report.warnings.iter().map(|i| i.kind).collect::<Vec<_>>(),
            [
                ReactionValidationIssueKind::MultipleMass,
                ReactionValidationIssueKind::MultipleIsotope,
                ReactionValidationIssueKind::MultipleCharge
            ]
        );
        for (key, value) in [
            ("_QueryFormalCharge", 1),
            ("_QueryMass", 12),
            ("_QueryIsotope", 7),
            ("_QueryHCount", 3),
        ] {
            assert_eq!(props(&reaction, key), Some(&PropertyValue::Int(value)));
        }
        let keys = reaction.products[0].atoms()[0]
            .property_records(true, true)
            .unwrap()
            .map(|(key, _)| key.as_bytes().to_vec())
            .filter(|key| key.starts_with(b"_Query"))
            .collect::<Vec<_>>();
        assert_eq!(
            keys,
            [
                b"_QueryFormalCharge".to_vec(),
                b"_QueryMass".to_vec(),
                b"_QueryIsotope".to_vec(),
                b"_QueryHCount".to_vec()
            ]
        );
    }

    #[test]
    fn source_query_write_is_immediate_before_later_undefined_arithmetic_error() {
        let query = Q::And(vec![
            Q::Predicate(P::HydrogenCount(2)),
            Q::Predicate(P::NegativeFormalCharge(i32::MIN)),
            Q::Predicate(P::Isotope(3)),
        ]);
        let mut reaction = paired(single(Some(1), query));
        let mut report = ReactionValidationReport::default();
        assert!(matches!(
            validate(&mut reaction, &mut report),
            Err(ReactionValidationError::QueryArithmetic { .. })
        ));
        assert_eq!(
            props(&reaction, "_QueryHCount"),
            Some(&PropertyValue::Int(2))
        );
        assert_eq!(props(&reaction, "_QueryFormalCharge"), None);
        assert_eq!(props(&reaction, "_QueryIsotope"), None);
    }

    #[test]
    fn conflicting_negative_minimum_does_not_eagerly_evaluate_source_skipped_write() {
        let query = Q::And(vec![
            Q::Predicate(P::FormalCharge(0)),
            Q::Predicate(P::NegativeFormalCharge(i32::MIN)),
            Q::Predicate(P::HydrogenCount(3)),
        ]);
        let mut reaction = paired(single(Some(1), query));
        let mut report = ReactionValidationReport::default();
        assert!(validate(&mut reaction, &mut report).unwrap());
        assert_eq!(report.warnings.len(), 1);
        assert_eq!(
            report.warnings[0].kind,
            ReactionValidationIssueKind::MultipleCharge
        );
        assert_eq!(
            props(&reaction, "_QueryFormalCharge"),
            Some(&PropertyValue::Int(0))
        );
        assert_eq!(
            props(&reaction, "_QueryHCount"),
            Some(&PropertyValue::Int(3))
        );
    }

    #[test]
    fn cache_clear_only_when_present_and_propagates_real_computed_metadata_error() {
        let mut reaction = paired(single(Some(1), Q::Predicate(P::Any)));
        reaction.products[0].atoms_mut()[0]
            .set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        let mut report = ReactionValidationReport::default();
        assert!(validate(&mut reaction, &mut report).unwrap());
        reaction.products[0].atoms_mut()[0]
            .set_prop("_QueryFormalCharge", 7i32)
            .unwrap();
        assert!(matches!(
            validate(&mut reaction, &mut report),
            Err(ReactionValidationError::Annotation {
                property: "_QueryFormalCharge",
                ..
            })
        ));
        assert_eq!(
            props(&reaction, "_QueryFormalCharge"),
            Some(&PropertyValue::Int(7))
        );
        assert_eq!(report.warnings.len(), 0);
    }

    #[test]
    fn entire_current_product_cache_clear_precedes_map_read_but_later_templates_are_untouched() {
        let mut product = graph(
            vec![
                atom(0, None, Q::Predicate(P::Any)),
                atom(1, None, Q::Predicate(P::Any)),
            ],
            &[],
        );
        product.atoms_mut()[0]
            .set_prop("molAtomMapNumber", PropertyValue::Bool(false))
            .unwrap();
        product.atoms_mut()[0]
            .set_prop("_QueryFormalCharge", 7i32)
            .unwrap();
        product.atoms_mut()[1]
            .set_prop("_QueryHCount", 99i32)
            .unwrap();
        let mut reaction = paired(product);
        let mut later = single(Some(1), Q::Predicate(P::Any));
        later.atoms_mut()[0].set_prop("_QueryMass", 88i32).unwrap();
        reaction.products.push(later);
        let mut report = ReactionValidationReport::default();
        assert!(matches!(
            validate(&mut reaction, &mut report),
            Err(ReactionValidationError::Property { .. })
        ));
        assert_eq!(props(&reaction, "_QueryFormalCharge"), None);
        assert_eq!(reaction.products[0].atoms()[1].prop("_QueryHCount"), None);
        assert_eq!(
            reaction.products[1].atoms()[0].prop("_QueryMass"),
            Some(&PropertyValue::Int(88))
        );
    }

    #[test]
    fn derived_carrier_placeholder_is_not_a_source_query() {
        let mut atom = QueryAtom::from_carrier_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Q::Predicate(P::FormalCharge(8)),
        );
        atom.set_atom_map(Some(1));
        let mut reaction = paired(graph(vec![atom], &[]));
        let mut report = ReactionValidationReport::default();
        assert!(validate(&mut reaction, &mut report).unwrap());
        assert_eq!(props(&reaction, "_QueryFormalCharge"), None);
        assert!(reaction.products[0].atoms()[0].predicate_is_carrier_derived());
    }
    mod complete_init_reactant_matchers_source_tests {
        use super::*;

        fn init(
            reaction: &mut Reaction,
            report: &mut ReactionValidationReport,
        ) -> Result<(), ReactionValidationError> {
            init_reactant_matchers_source(
                reaction,
                &ReactionValidationParams { silent: true },
                report,
            )
        }

        #[test]
        fn normal_validation_success_resets_the_flag_from_either_prior_state() {
            for prior in [false, true] {
                let mut reaction = paired(single(Some(1), Q::Predicate(P::HydrogenCount(2))));
                reaction.needs_init = prior;
                reaction.match_params.max_matches = 37;
                reaction.match_params.use_chirality = true;
                reaction.match_params.atom_properties = vec!["source-init-retains-options".into()];
                let property_storage = reaction.match_params.atom_properties.as_ptr();
                let mut report = ReactionValidationReport::default();
                init(&mut reaction, &mut report).unwrap();
                assert!(!reaction.needs_init);
                assert!(report.is_valid());
                assert_eq!(
                    props(&reaction, "_QueryHCount"),
                    Some(&PropertyValue::Int(2))
                );
                assert_eq!(reaction.match_params.max_matches, 37);
                assert!(reaction.match_params.use_chirality);
                assert_eq!(
                    reaction.match_params.atom_properties,
                    ["source-init-retains-options"]
                );
                assert_eq!(
                    reaction.match_params.atom_properties.as_ptr(),
                    property_storage
                );
            }
        }

        #[test]
        fn normal_validation_false_returns_success_and_keeps_annotations_from_either_prior_state() {
            for prior in [false, true] {
                let mut reaction = paired(single(Some(1), Q::Predicate(P::HydrogenCount(2))));
                reaction.reactants.push(reaction.reactants[0].clone());
                reaction.needs_init = prior;
                let mut report = ReactionValidationReport::default();
                init(&mut reaction, &mut report).unwrap();
                assert!(reaction.needs_init);
                assert_eq!(report.num_errors(), 1);
                assert_eq!(
                    report.errors[0].kind,
                    ReactionValidationIssueKind::DuplicateReactantMap
                );
                assert_eq!(
                    props(&reaction, "_QueryHCount"),
                    Some(&PropertyValue::Int(2))
                );
            }
        }

        #[test]
        fn validation_exception_preserves_prior_flag_and_actual_annotation_prefix() {
            for prior in [false, true] {
                let mut reaction = paired(single(
                    Some(1),
                    Q::And(vec![
                        Q::Predicate(P::HydrogenCount(2)),
                        Q::Predicate(P::NegativeFormalCharge(i32::MIN)),
                        Q::Predicate(P::Isotope(13)),
                    ]),
                ));
                reaction.needs_init = prior;
                let mut report = ReactionValidationReport::default();
                assert!(matches!(
                    init(&mut reaction, &mut report),
                    Err(ReactionValidationError::QueryArithmetic { .. })
                ));
                assert_eq!(reaction.needs_init, prior);
                assert_eq!(
                    props(&reaction, "_QueryHCount"),
                    Some(&PropertyValue::Int(2))
                );
                assert_eq!(props(&reaction, "_QueryIsotope"), None);
            }
        }

        #[test]
        fn each_initialization_resets_previous_validation_outputs() {
            let mut report = ReactionValidationReport::default();
            init(&mut Reaction::new(), &mut report).unwrap();
            assert_eq!(report.num_errors(), 2);
            let mut reaction = paired(single(Some(1), Q::Predicate(P::Any)));
            init(&mut reaction, &mut report).unwrap();
            assert!(report.is_valid());
            assert_eq!(report.num_warnings(), 0);
            assert!(!reaction.needs_init);
        }

        #[test]
        fn immutable_projection_delegates_and_returns_structured_invalid_without_mutating_input() {
            let mut reaction = paired(single(Some(1), Q::Predicate(P::HydrogenCount(2))));
            reaction.needs_init = false;
            let output =
                initialize_reaction(&reaction, &ReactionValidationParams { silent: true }).unwrap();
            assert!(!output.needs_init);
            assert_eq!(props(&output, "_QueryHCount"), Some(&PropertyValue::Int(2)));
            assert_eq!(props(&reaction, "_QueryHCount"), None);
            reaction.reactants.push(reaction.reactants[0].clone());
            assert!(
                matches!(initialize_reaction(&reaction, &ReactionValidationParams { silent: true }),
                Err(ReactionInitializationError::Invalid { report }) if report.num_errors() == 1)
            );
            assert!(!reaction.needs_init);
            assert_eq!(props(&reaction, "_QueryHCount"), None);
        }
    }
}
