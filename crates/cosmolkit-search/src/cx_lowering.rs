//! QueryGraph-native lowering for representation-independent CX records.

use cosmolkit_cx::{
    CxAtomConstraint, CxCoordinateBondKind, CxCountConstraint, CxDoubleBondStereoKind, CxRecord,
    CxStereoGroupKind, CxWedgeDirection, ParsedCxExtensions,
};
use cosmolkit_model::{
    AtomId, AtomQueryPredicate, BondDirection, BondOrder, BondStereo, Conformer3D, QueryGraph,
    QueryNode, StereoGroup, StereoGroupKind,
};

const QUERY_SCAN_MAGIC_VALUE: u32 = 0xDEAD_BEEF;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CxQueryLoweringError {
    #[error("CX atom index {index} is outside the query graph")]
    AtomIndex { index: usize },
    #[error("CX bond index {index} is outside the query graph")]
    BondIndex { index: usize },
    #[error("CX coordinate count {actual} does not match query atom count {expected}")]
    CoordinateCount { actual: usize, expected: usize },
    #[error("CX coordinate bond atom {atom} is not an endpoint of bond {bond}")]
    BondAtomMismatch { atom: usize, bond: usize },
    #[error("CX record has no QueryGraph representation: {record}")]
    UnsupportedRecord { record: &'static str },
    #[error("query graph is invalid after CX lowering: {0}")]
    InvalidGraph(String),
}

fn append_atom_predicate(
    graph: &mut QueryGraph,
    atom: usize,
    predicate: AtomQueryPredicate,
) -> Result<(), CxQueryLoweringError> {
    let query_atom = graph
        .atom_mut(atom)
        .ok_or(CxQueryLoweringError::AtomIndex { index: atom })?;
    let current = std::mem::replace(query_atom.predicate_mut(), QueryNode::and(Vec::new()));
    *query_atom.predicate_mut() = QueryNode::and(vec![current, QueryNode::predicate(predicate)]);
    Ok(())
}

/// Apply parsed CX records directly to the canonical query value.
///
/// Parsing remains owned by `cosmolkit-cx`; this function owns only the
/// destination semantics. It never projects query data through a concrete
/// molecule.
pub fn apply_cx_to_query_graph(
    graph: &mut QueryGraph,
    parsed: &ParsedCxExtensions,
) -> Result<(), CxQueryLoweringError> {
    for record in parsed.records() {
        match record {
            CxRecord::Coordinates(coordinates) => {
                let values = coordinates
                    .values
                    .iter()
                    .map(|value| value.unwrap_or([0.0; 3]))
                    .collect::<Vec<_>>();
                if values.len() != graph.num_atoms() {
                    return Err(CxQueryLoweringError::CoordinateCount {
                        actual: values.len(),
                        expected: graph.num_atoms(),
                    });
                }
                graph
                    .add_conformer_3d(Conformer3D::new(
                        coordinates.conformer,
                        values,
                        coordinates.is_3d,
                    ))
                    .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
            }
            CxRecord::AtomLabels(values) => {
                for (index, value) in values.iter().enumerate() {
                    if let Some(value) = value {
                        let atom = graph
                            .atom_mut(index)
                            .ok_or(CxQueryLoweringError::AtomIndex { index })?;
                        atom.set_prop("atomLabel", value);
                    }
                }
            }
            CxRecord::AtomValues(values) => {
                for (index, value) in values.iter().enumerate() {
                    if let Some(value) = value {
                        let atom = graph
                            .atom_mut(index)
                            .ok_or(CxQueryLoweringError::AtomIndex { index })?;
                        atom.set_prop("molFileValue", value);
                    }
                }
            }
            CxRecord::AtomProperties(properties) => {
                for property in properties {
                    let atom =
                        graph
                            .atom_mut(property.atom)
                            .ok_or(CxQueryLoweringError::AtomIndex {
                                index: property.atom,
                            })?;
                    atom.set_prop(property.name.clone(), property.value.clone());
                }
            }
            CxRecord::CoordinateBonds(annotation) => {
                let order = match annotation.kind {
                    CxCoordinateBondKind::Dative => BondOrder::Dative,
                    CxCoordinateBondKind::Hydrogen => BondOrder::Hydrogen,
                };
                for reference in &annotation.bonds {
                    let bond = graph.bonds_mut().get_mut(reference.bond).ok_or(
                        CxQueryLoweringError::BondIndex {
                            index: reference.bond,
                        },
                    )?;
                    if bond.begin().index() != reference.atom
                        && bond.end().index() != reference.atom
                    {
                        return Err(CxQueryLoweringError::BondAtomMismatch {
                            atom: reference.atom,
                            bond: reference.bond,
                        });
                    }
                    bond.bond_mut().set_order(order);
                }
            }
            CxRecord::ZeroBonds(indices) => {
                for &index in indices {
                    let bond = graph
                        .bonds_mut()
                        .get_mut(index)
                        .ok_or(CxQueryLoweringError::BondIndex { index })?;
                    bond.bond_mut().set_order(BondOrder::Zero);
                }
            }
            CxRecord::Unsaturation(indices) => {
                for &index in indices {
                    append_atom_predicate(graph, index, AtomQueryPredicate::IsUnsaturated)?;
                }
            }
            CxRecord::RingBonds(constraints) => {
                for constraint in constraints {
                    let predicate = match constraint.constraint {
                        CxCountConstraint::Exact(value) => AtomQueryPredicate::RingBondCount(
                            i32::try_from(value)
                                .expect("CX ring-bond equality is parser-bounded to 0, 2, or 3"),
                        ),
                        CxCountConstraint::LessEqual(value) => {
                            AtomQueryPredicate::RingBondCountLessEqual(value as u8)
                        }
                        CxCountConstraint::QueryScan => {
                            AtomQueryPredicate::RingBondCount(QUERY_SCAN_MAGIC_VALUE as i32)
                        }
                    };
                    append_atom_predicate(graph, constraint.atom, predicate)?;
                }
            }
            CxRecord::Substitution(constraints) => {
                for CxAtomConstraint { atom, constraint } in constraints {
                    let predicate = match constraint {
                        CxCountConstraint::Exact(value) => {
                            AtomQueryPredicate::NonHydrogenDegree(*value)
                        }
                        CxCountConstraint::LessEqual(value) => {
                            AtomQueryPredicate::NonHydrogenDegreeLessEqual(*value)
                        }
                        CxCountConstraint::QueryScan => {
                            AtomQueryPredicate::NonHydrogenDegree(QUERY_SCAN_MAGIC_VALUE)
                        }
                    };
                    append_atom_predicate(graph, *atom, predicate)?;
                }
            }
            CxRecord::EnhancedStereo(stereo) => {
                let kind = match stereo.kind {
                    CxStereoGroupKind::Absolute => StereoGroupKind::Absolute,
                    CxStereoGroupKind::Or => StereoGroupKind::Or,
                    CxStereoGroupKind::And => StereoGroupKind::And,
                };
                let atoms = stereo
                    .atoms
                    .iter()
                    .map(|&index| {
                        (index < graph.num_atoms())
                            .then_some(AtomId::new(index))
                            .ok_or(CxQueryLoweringError::AtomIndex { index })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                if !atoms.is_empty() {
                    graph.add_stereo_group(
                        StereoGroup::new(kind, atoms, Vec::new()).with_id(stereo.group_id),
                    );
                }
            }
            CxRecord::WedgedBonds(wedges) => {
                for wedge in wedges {
                    let bond = graph
                        .bonds_mut()
                        .get_mut(wedge.bond)
                        .ok_or(CxQueryLoweringError::BondIndex { index: wedge.bond })?;
                    let atom = AtomId::new(wedge.atom);
                    if bond.begin() != atom && bond.end() != atom {
                        return Err(CxQueryLoweringError::BondAtomMismatch {
                            atom: wedge.atom,
                            bond: wedge.bond,
                        });
                    }
                    if bond.begin() != atom {
                        let begin = bond.begin();
                        bond.bond_mut().set_endpoints(atom, begin);
                    }
                    let (configuration, direction) = match wedge.direction {
                        CxWedgeDirection::Unknown => ("2", BondDirection::Unknown),
                        CxWedgeDirection::BeginWedge => ("1", BondDirection::BeginWedge),
                        CxWedgeDirection::BeginDash => ("3", BondDirection::BeginDash),
                    };
                    bond.bond_mut().set_prop("_MolFileBondCfg", configuration);
                    bond.bond_mut().set_direction(direction);
                }
            }
            CxRecord::DoubleBondStereo(stereo) => {
                let value = match stereo.stereo {
                    CxDoubleBondStereoKind::Any => BondStereo::Any,
                    CxDoubleBondStereoKind::Cis => BondStereo::Cis,
                    CxDoubleBondStereoKind::Trans => BondStereo::Trans,
                };
                for &index in &stereo.bonds {
                    let bond = graph
                        .bonds_mut()
                        .get_mut(index)
                        .ok_or(CxQueryLoweringError::BondIndex { index })?;
                    bond.bond_mut().set_stereo(value);
                }
            }
            CxRecord::Radicals(radicals) => {
                for radical in radicals {
                    let atom =
                        graph
                            .atom_mut(radical.atom)
                            .ok_or(CxQueryLoweringError::AtomIndex {
                                index: radical.atom,
                            })?;
                    atom.set_radical_electrons(radical.electrons);
                }
            }
            CxRecord::LinkNodes(_) => {
                return Err(CxQueryLoweringError::UnsupportedRecord {
                    record: "link node",
                });
            }
            CxRecord::DataSGroup(_) => {
                return Err(CxQueryLoweringError::UnsupportedRecord {
                    record: "data SGroup",
                });
            }
            CxRecord::SGroupHierarchy(_) => {
                return Err(CxQueryLoweringError::UnsupportedRecord {
                    record: "SGroup hierarchy",
                });
            }
            CxRecord::PolymerSGroup(_) => {
                return Err(CxQueryLoweringError::UnsupportedRecord {
                    record: "polymer SGroup",
                });
            }
            CxRecord::VariableAttachments(_) => {
                return Err(CxQueryLoweringError::UnsupportedRecord {
                    record: "variable attachment",
                });
            }
            CxRecord::Unknown(_) => {}
        }
    }
    graph
        .validate()
        .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query_behavior::{
        make_atom_in_ring_of_size_query, make_atom_min_ring_size_query,
        make_atom_ring_bond_count_query, make_atom_ring_query,
    };
    use cosmolkit_cx::{CxAtomConstraint, CxCountConstraint, CxRingBond, ParsedCxExtensions};
    use cosmolkit_model::{AtomSpec, BondId, BondSpec};
    use cosmolkit_types::Element;

    fn graph() -> QueryGraph {
        let atoms = vec![
            cosmolkit_model::QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
            cosmolkit_model::QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![cosmolkit_model::QueryBond::new(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        QueryGraph::from_parts(
            atoms,
            bonds,
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap()
    }

    #[test]
    fn lowers_query_constraints_without_concrete_projection() {
        let parsed = ParsedCxExtensions::new(
            vec![
                CxRecord::Unsaturation(vec![0]),
                CxRecord::Substitution(vec![CxAtomConstraint {
                    atom: 1,
                    constraint: CxCountConstraint::Exact(1),
                }]),
                CxRecord::AtomLabels(vec![Some("left".to_owned()), None]),
            ],
            8,
        );
        let mut query = graph();
        apply_cx_to_query_graph(&mut query, &parsed).unwrap();
        assert_eq!(query.atom(0).unwrap().prop("atomLabel"), Some("left"));
        assert!(matches!(
            query.atom(0).unwrap().predicate(),
            QueryNode::And(children) if children.len() == 2
        ));
        assert!(matches!(
            query.atom(1).unwrap().predicate(),
            QueryNode::And(children) if children.len() == 2
        ));
    }

    #[test]
    fn q07e_ring_factories_and_cx_lowering_preserve_i32_targets_and_sentinels() {
        let maximum = 2_147_483_639;
        for target in [0, 255, 256, maximum] {
            assert_eq!(
                make_atom_ring_query(target),
                QueryNode::predicate(AtomQueryPredicate::NumAtomRings(target))
            );
            assert_eq!(
                make_atom_in_ring_of_size_query(target),
                QueryNode::predicate(AtomQueryPredicate::InRingOfSize(target))
            );
            assert_eq!(
                make_atom_min_ring_size_query(target),
                QueryNode::predicate(AtomQueryPredicate::SmallestRingSize(target))
            );
            assert_eq!(
                make_atom_ring_bond_count_query(target),
                QueryNode::predicate(AtomQueryPredicate::RingBondCount(target))
            );
        }
        assert_eq!(
            make_atom_ring_query(-1),
            QueryNode::predicate(AtomQueryPredicate::NumAtomRings(-1))
        );
        assert_eq!(
            make_atom_ring_bond_count_query(i32::MIN),
            QueryNode::predicate(AtomQueryPredicate::RingBondCount(i32::MIN))
        );

        fn contains_predicate(
            node: &QueryNode<AtomQueryPredicate>,
            expected: &AtomQueryPredicate,
        ) -> bool {
            match node {
                QueryNode::Predicate(predicate) => predicate == expected,
                QueryNode::And(children) | QueryNode::Or(children) => children
                    .iter()
                    .any(|child| contains_predicate(child, expected)),
                _ => false,
            }
        }

        let mut query = graph();
        let carrier_before = query.atom(0).unwrap().try_to_atom().unwrap();
        let parsed = ParsedCxExtensions::new(
            vec![CxRecord::RingBonds(vec![
                CxRingBond {
                    atom: 0,
                    constraint: CxCountConstraint::Exact(3),
                },
                CxRingBond {
                    atom: 0,
                    constraint: CxCountConstraint::QueryScan,
                },
                CxRingBond {
                    atom: 1,
                    constraint: CxCountConstraint::LessEqual(4),
                },
            ])],
            8,
        );
        apply_cx_to_query_graph(&mut query, &parsed).unwrap();

        let atom_zero = query.atom(0).unwrap();
        assert!(contains_predicate(
            atom_zero.predicate(),
            &AtomQueryPredicate::RingBondCount(3)
        ));
        assert!(contains_predicate(
            atom_zero.predicate(),
            &AtomQueryPredicate::RingBondCount(0xDEAD_BEEF_u32 as i32)
        ));
        assert_eq!(atom_zero.try_to_atom().unwrap(), carrier_before);
        assert!(!atom_zero.predicate_is_carrier_derived());
        assert!(contains_predicate(
            query.atom(1).unwrap().predicate(),
            &AtomQueryPredicate::RingBondCountLessEqual(4)
        ));
    }

    #[test]
    fn rejects_concrete_only_cx_records() {
        let parsed = ParsedCxExtensions::new(vec![CxRecord::LinkNodes(Vec::new())], 4);
        let mut query = graph();
        assert_eq!(
            apply_cx_to_query_graph(&mut query, &parsed),
            Err(CxQueryLoweringError::UnsupportedRecord {
                record: "link node"
            })
        );
    }

    #[test]
    fn rejects_out_of_range_cx_indices() {
        let parsed = ParsedCxExtensions::new(
            vec![CxRecord::Radicals(vec![cosmolkit_cx::CxRadical {
                atom: 3,
                electrons: 1,
            }])],
            4,
        );
        let mut query = graph();
        assert_eq!(
            apply_cx_to_query_graph(&mut query, &parsed),
            Err(CxQueryLoweringError::AtomIndex { index: 3 })
        );
    }
}
