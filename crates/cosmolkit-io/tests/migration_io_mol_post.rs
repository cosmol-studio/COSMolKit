//! Molfile postprocessing regressions against pinned RDKit 2026.03.1.

use cosmolkit_io::{
    MolBlockRecord, MolPostError, MolPostParams, QueryMolBlockRecord, finish_mol_block_record,
    read_mol_block_detached,
};
use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondQueryPredicate, BondSpec,
    Conformer3D, CoordinateBlock, CoordinateDimension, MoleculeProperties, QueryAtom, QueryBond,
    QueryGraph, QueryNode, RecursiveStructureQuery, SGroupData, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Element, Hybridization};

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>, groups: Vec<SubstanceGroup>) -> TopologyBlock {
    TopologyBlock::try_from_parts(atoms, bonds, groups, Vec::new()).unwrap()
}

fn concrete(topology: TopologyBlock) -> MolBlockRecord {
    MolBlockRecord::Concrete {
        topology,
        coordinates: CoordinateBlock::default(),
        properties: MoleculeProperties::default(),
    }
}

fn concrete_or_explicit_query(topology: TopologyBlock, query: bool) -> MolBlockRecord {
    if !query {
        return concrete(topology);
    }
    let atoms = topology
        .atoms
        .into_iter()
        .enumerate()
        .map(|(index, atom)| {
            let atomic_number = if index == 0 {
                Element::N.atomic_number()
            } else {
                atom.element().atomic_number()
            };
            QueryAtom::from_parts(
                atom,
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atomic_number)),
            )
        })
        .collect();
    let bonds = topology
        .bonds
        .into_iter()
        .map(|bond| QueryBond::from_parts(bond, QueryNode::predicate(BondQueryPredicate::Any)))
        .collect();
    MolBlockRecord::Query(QueryMolBlockRecord {
        query: QueryGraph::from_parts(
            atoms,
            bonds,
            Default::default(),
            vec![],
            vec![],
            topology.stereo_groups,
        )
        .unwrap(),
        substance_groups: topology.substance_groups,
        properties: MoleculeProperties::default(),
        source_coordinate_dim: None,
    })
}

fn dat_group(
    id: usize,
    field_name: Option<&str>,
    query_type: Option<&str>,
    query_op: Option<&str>,
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
    values: &[&str],
) -> SubstanceGroup {
    SubstanceGroup::new(SubstanceGroupId::new(id), SubstanceGroupKind::Data)
        .with_atoms(atoms)
        .with_bonds(bonds)
        .with_data(SGroupData {
            field_name: field_name.map(str::to_owned),
            query_type: query_type.map(str::to_owned),
            query_op: query_op.map(str::to_owned),
            values: values.iter().map(|value| (*value).to_owned()).collect(),
            ..SGroupData::default()
        })
}

fn unsanitized() -> MolPostParams {
    MolPostParams {
        sanitize: false,
        remove_hs: false,
        expand_attachment_points: false,
    }
}

fn v3000_attachment(value: Option<&str>, subst: Option<i32>) -> String {
    let attachment = value.map_or(String::new(), |value| format!(" ATTCHPT={value}"));
    let subst = subst.map_or(String::new(), |value| format!(" SUBST={value}"));
    format!(
        "attachment\n  COSMolKit\n\n  0  0  0     0  0            999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 2 1 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0\n\
M  V30 2 O 1 0 0 0{attachment}{subst}\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n"
    )
}

fn v2000_attachment(value: u8) -> String {
    let carbon = format!(
        "{:>10.4}{:>10.4}{:>10.4} C   0  0  0  0  0  0  0  0  0  0  0  0",
        0.0, 0.0, 0.0
    );
    let oxygen = format!(
        "{:>10.4}{:>10.4}{:>10.4} O   0  0  0  0  0  0  0  0  0  0  0  0",
        1.0, 0.0, 0.0
    );
    format!(
        "attachment\n  COSMolKit\n\n  2  1  0  0  0  0            999 V2000\n{carbon}\n{oxygen}\n  1  2  1  0  0  0  0\nM  APO  1   2 {:>3}\nM  END\n",
        value
    )
}

#[test]
fn mol_post_attachment_v3000_value_and_final_classification() {
    for (value, labels) in [
        (None, &[][..]),
        (Some("0"), &[][..]),
        (Some("1"), &["1"][..]),
        (Some("2"), &["2"][..]),
        (Some("-1"), &["1", "2"][..]),
        (Some("3"), &[][..]),
        (Some("+1"), &[][..]),
    ] {
        let input = v3000_attachment(value, None);
        let parsed = read_mol_block_detached(&input).unwrap();
        let original = parsed.clone();
        let finished = finish_mol_block_record(
            parsed,
            false,
            MolPostParams {
                sanitize: false,
                remove_hs: false,
                expand_attachment_points: true,
            },
        )
        .unwrap();
        if labels.is_empty() {
            let MolBlockRecord::Concrete { topology, .. } = finished else {
                panic!("no append must not promote: {value:?}");
            };
            assert_eq!(topology.atoms.len(), 2);
            assert_eq!(topology.bonds.len(), 1);
        } else {
            let MolBlockRecord::Query(record) = finished else {
                panic!("source null query must promote: {value:?}");
            };
            assert_eq!(record.query.num_atoms(), 2 + labels.len());
            assert_eq!(record.query.num_bonds(), 1 + labels.len());
            for (offset, label) in labels.iter().enumerate() {
                let atom = &record.query.atoms()[2 + offset];
                assert_eq!(
                    atom.predicate(),
                    &QueryNode::predicate(AtomQueryPredicate::Any)
                );
                assert!(!atom.predicate_is_carrier_derived());
                assert_eq!(atom.prop("_fromAttchpt"), Some(*label));
                assert_eq!(
                    record.query.bonds()[1 + offset].bond().order(),
                    BondOrder::Single
                );
                assert!(record.query.bonds()[1 + offset].predicate_is_carrier_derived());
            }
            assert_eq!(record.query.atoms()[1].prop("molAttachPoint"), None);
            let coordinate_rows = record.query.coordinates_2d().map(<[_]>::len).or_else(|| {
                record
                    .query
                    .conformers_3d()
                    .first()
                    .map(|conformer| conformer.coordinates().len())
            });
            assert_eq!(coordinate_rows, Some(2 + labels.len()));
        }
        assert_eq!(original, read_mol_block_detached(&input).unwrap());
    }
}

#[test]
fn mol_post_attachment_v2000_and_existing_query_obey_source_order() {
    for (value, count) in [(1, 3), (2, 3), (3, 4)] {
        let parsed = read_mol_block_detached(&v2000_attachment(value)).unwrap();
        let finished = finish_mol_block_record(
            parsed,
            false,
            MolPostParams {
                sanitize: false,
                remove_hs: false,
                expand_attachment_points: true,
            },
        )
        .unwrap();
        let MolBlockRecord::Query(record) = finished else {
            panic!("V2000 APO must append source null-query atoms");
        };
        assert_eq!(record.query.num_atoms(), count);
    }
    let parsed = read_mol_block_detached(&v3000_attachment(Some("1"), Some(1))).unwrap();
    let finished = finish_mol_block_record(
        parsed,
        false,
        MolPostParams {
            sanitize: false,
            remove_hs: false,
            expand_attachment_points: true,
        },
    )
    .unwrap();
    let MolBlockRecord::Query(record) = finished else {
        panic!("SUBST is a query");
    };
    assert_eq!(record.query.num_atoms(), 3);
    assert_eq!(
        record.query.atoms()[2].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::Any)
    );
    assert_eq!(record.query.atoms()[1].prop("molAttachPoint"), None);

    // ProcessMolProps maps SUBST=-2 to the degree *after* attachment expansion.
    let parsed = read_mol_block_detached(&v3000_attachment(Some("1"), Some(-2))).unwrap();
    let finished = finish_mol_block_record(
        parsed,
        false,
        MolPostParams {
            sanitize: false,
            remove_hs: false,
            expand_attachment_points: true,
        },
    )
    .unwrap();
    let MolBlockRecord::Query(record) = finished else {
        panic!("SUBST is a query");
    };
    assert_eq!(record.query.num_atoms(), 3);
    assert_eq!(
        record.query.atoms()[1].predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
            QueryNode::predicate(AtomQueryPredicate::ExplicitDegree(2)),
        ]),
    );
}

#[test]
fn mol_post_attachment_options_and_invalid_local_value_are_atomic() {
    for sanitize in [false, true] {
        for remove_hs in [false, true] {
            let parsed = read_mol_block_detached(&v3000_attachment(Some("1"), None)).unwrap();
            let source = parsed.clone();
            let expanded = finish_mol_block_record(
                parsed,
                false,
                MolPostParams {
                    sanitize,
                    remove_hs,
                    expand_attachment_points: true,
                },
            )
            .unwrap();
            let MolBlockRecord::Query(record) = expanded else {
                panic!("attachment must promote");
            };
            assert_eq!(
                record.query.num_atoms(),
                3,
                "sanitize={sanitize} remove_hs={remove_hs}"
            );
            assert_eq!(
                source,
                read_mol_block_detached(&v3000_attachment(Some("1"), None)).unwrap()
            );
            let disabled = finish_mol_block_record(
                source,
                false,
                MolPostParams {
                    sanitize,
                    remove_hs,
                    expand_attachment_points: false,
                },
            )
            .unwrap();
            assert!(matches!(disabled, MolBlockRecord::Concrete { .. }));
        }
    }
    let mut invalid = concrete(topology(vec![atom(0, Element::C)], vec![], vec![]));
    if let MolBlockRecord::Concrete { topology, .. } = &mut invalid {
        topology.atoms[0]
            .set_prop("molAttachPoint", "nonsense")
            .unwrap();
    }
    let original = invalid.clone();
    assert!(matches!(
        finish_mol_block_record(invalid.clone(), false, MolPostParams { expand_attachment_points: true, ..MolPostParams::default() }),
        Err(MolPostError::AttachmentValue { atom, value }) if atom == AtomId::new(0) && value == "nonsense"
    ));
    assert_eq!(invalid, original);
}

#[test]
fn mol_post_explicit_valence_prepass_runs_for_concrete_and_query_before_properties() {
    // MolFileParser.cpp::finishMolProcessing calls calcExplicitValence(false)
    // for every atom before ProcessMolProps, regardless of sanitize/removeHs.
    // Bond.cpp::getBondTypeAsDouble rejects OTHER with "Bad bond type".
    for query in [false, true] {
        for substitution in [None, Some("300")] {
            let mut first = atom(0, Element::C);
            if let Some(value) = substitution {
                first.set_prop("molSubstCount", value).unwrap();
            }
            let input = concrete_or_explicit_query(
                topology(
                    vec![first, atom(1, Element::C)],
                    vec![Bond::from_spec(
                        BondId::new(0),
                        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Other),
                    )],
                    vec![],
                ),
                query,
            );
            for sanitize in [false, true] {
                for remove_hs in [false, true] {
                    let result = finish_mol_block_record(
                        input.clone(),
                        false,
                        MolPostParams {
                            sanitize,
                            remove_hs,
                            expand_attachment_points: false,
                        },
                    );
                    assert!(
                        matches!(result, Err(MolPostError::Processing(ref message)) if message == "Bad bond type"),
                        "query={query} substitution={substitution:?} sanitize={sanitize} remove_hs={remove_hs}: {result:?}"
                    );
                    assert_eq!(
                        input,
                        concrete_or_explicit_query(
                            topology(
                                vec![
                                    {
                                        let mut atom = atom(0, Element::C);
                                        if let Some(value) = substitution {
                                            atom.set_prop("molSubstCount", value).unwrap();
                                        }
                                        atom
                                    },
                                    atom(1, Element::C),
                                ],
                                vec![Bond::from_spec(
                                    BondId::new(0),
                                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Other),
                                )],
                                vec![],
                            ),
                            query,
                        ),
                    );
                }
            }
        }
    }

    // An ordinary, source-valid attachment input must still reach expansion.
    let valid = read_mol_block_detached(&v3000_attachment(Some("1"), None)).unwrap();
    let result = finish_mol_block_record(
        valid,
        false,
        MolPostParams {
            sanitize: false,
            remove_hs: false,
            expand_attachment_points: true,
        },
    )
    .unwrap();
    assert!(matches!(result, MolBlockRecord::Query(record) if record.query.num_atoms() == 3));
}

fn v3000_tetrahedral_wedge(z: &str, cfg: u8) -> String {
    format!(
        "test\n  COSMolKit\n\n  0  0  0     0  0            999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 5 4 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 {z} 0\n\
M  V30 2 F 1 0 0 0\n\
M  V30 3 Cl 0 1 0 0\n\
M  V30 4 Br -1 0 0 0\n\
M  V30 5 I 0 -1 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2 CFG={cfg}\n\
M  V30 2 1 1 3\n\
M  V30 3 1 1 4\n\
M  V30 4 1 1 5\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n"
    )
}

fn v3000_tetrahedral_case(
    ligands: [&str; 4],
    query_center: bool,
    z: &str,
    cfg: u8,
    absolute_group: bool,
) -> String {
    let query = if query_center { " RBCNT=2" } else { "" };
    let collection = if absolute_group {
        "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/STEABS ATOMS=(1 1)\nM  V30 END COLLECTION\n"
    } else {
        ""
    };
    format!(
        "test\n  COSMolKit\n\n  0  0  0     0  0            999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 5 4 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 {z} 0{query}\n\
M  V30 2 {} 1 0 0 0\n\
M  V30 3 {} 0 1 0 0\n\
M  V30 4 {} -1 0 0 0\n\
M  V30 5 {} 0 -1 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2 CFG={cfg}\n\
M  V30 2 1 1 3\n\
M  V30 3 1 1 4\n\
M  V30 4 1 1 5\n\
M  V30 END BOND\n\
{collection}M  V30 END CTAB\n\
M  END\n",
        ligands[0], ligands[1], ligands[2], ligands[3]
    )
}

fn finished_center_tag_and_group_count(
    block: &str,
    sanitize: bool,
    remove_hs: bool,
) -> (ChiralTag, usize, CoordinateBlock) {
    let parsed = read_mol_block_detached(block).expect("tetrahedral mol block");
    let finished = finish_mol_block_record(
        parsed,
        true,
        MolPostParams {
            sanitize,
            remove_hs,
            expand_attachment_points: false,
        },
    )
    .expect("tetrahedral mol-post");
    match finished {
        MolBlockRecord::Concrete {
            topology,
            coordinates,
            ..
        } => (
            topology.atoms[0].chiral_tag(),
            topology.stereo_groups.len(),
            coordinates,
        ),
        MolBlockRecord::Query(record) => (
            record.query.atoms()[0].chiral_tag(),
            record.query.stereo_groups().len(),
            record.query.coordinate_block(record.source_coordinate_dim),
        ),
    }
}

fn v3000_alkene_stereo_case(query_atom: bool) -> String {
    let query = if query_atom { " RBCNT=1" } else { "" };
    format!(
        "test\n  COSMolKit\n\n  0  0  0     0  0            999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 4 3 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 F -1 1 0 0\n\
M  V30 2 C 0 0 0 0{query}\n\
M  V30 3 C 1 0 0 0\n\
M  V30 4 F 2 -1 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 2 2 2 3\n\
M  V30 3 1 3 4\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n"
    )
}

fn v3000_benzene_query_case() -> String {
    "benzene query\n  COSMolKit\n\n  0  0  0     0  0            999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 6 6 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 1 0 0 0 RBCNT=2\n\
M  V30 2 C 0.5 0.866 0 0\n\
M  V30 3 C -0.5 0.866 0 0\n\
M  V30 4 C -1 0 0 0\n\
M  V30 5 C -0.5 -0.866 0 0\n\
M  V30 6 C 0.5 -0.866 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 2 1 2\n\
M  V30 2 1 2 3\n\
M  V30 3 2 3 4\n\
M  V30 4 1 4 5\n\
M  V30 5 2 5 6\n\
M  V30 6 1 6 1\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n"
        .to_owned()
}

fn v3000_iterative_legacy_stereo_case(query_record: bool) -> String {
    let query = if query_record { " RBCNT=1" } else { "" };
    format!(
        "iterative legacy stereo\n  RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 9 8 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 F -1.299038 0.750000 0.000000 0\n\
M  V30 2 C 0.000000 0.000000 0.000000 0{query}\n\
M  V30 3 Cl -0.000000 -1.500000 0.000000 0\n\
M  V30 4 C 1.299038 0.750000 0.000000 0\n\
M  V30 5 Br 0.549038 2.049038 0.000000 0\n\
M  V30 6 I 2.049038 -0.549038 0.000000 0\n\
M  V30 7 C 2.598076 1.500000 0.000000 0\n\
M  V30 8 F 2.598076 3.000000 0.000000 0\n\
M  V30 9 Cl 3.897114 0.750000 0.000000 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 2 1 CFG=1\n\
M  V30 2 1 2 3\n\
M  V30 3 1 2 4\n\
M  V30 4 1 4 5 CFG=3\n\
M  V30 5 1 4 6\n\
M  V30 6 1 4 7\n\
M  V30 7 1 7 8 CFG=1\n\
M  V30 8 1 7 9\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n"
    )
}

fn v3000_ring_special_case(query_record: bool) -> String {
    let query = if query_record { " RBCNT=2" } else { "" };
    format!(
        "ring special legacy stereo\n  RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 8 8 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 3.000000 0.000000 0.000000 0\n\
M  V30 2 C 1.500000 0.000000 0.000000 0{query}\n\
M  V30 3 C 0.750000 -1.299038 0.000000 0\n\
M  V30 4 C -0.750000 -1.299038 0.000000 0\n\
M  V30 5 C -1.500000 0.000000 0.000000 0\n\
M  V30 6 C -3.000000 0.000000 0.000000 0\n\
M  V30 7 C -0.750000 1.299038 0.000000 0\n\
M  V30 8 C 0.750000 1.299038 0.000000 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 2 1 CFG=1\n\
M  V30 2 1 2 3\n\
M  V30 3 1 3 4\n\
M  V30 4 1 4 5\n\
M  V30 5 1 5 6 CFG=1\n\
M  V30 6 1 5 7\n\
M  V30 7 1 7 8\n\
M  V30 8 1 8 2\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n"
    )
}

fn v3000_query_hydrogen(query_hydrogen: bool, substitution_count: Option<i32>) -> String {
    let carbon_query = if query_hydrogen { "" } else { " RBCNT=1" };
    let hydrogen_query = if query_hydrogen { " RBCNT=1" } else { "" };
    let substitution = substitution_count
        .map(|value| format!(" SUBST={value}"))
        .unwrap_or_default();
    format!(
        "query hydrogen\n  COSMolKit\n\n  0  0  0     0  0            999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 2 1 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0{carbon_query}{substitution}\n\
M  V30 2 H 1 0 0 0{hydrogen_query}\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n"
    )
}

fn finish_v3000_wedge(z: &str, cfg: u8) -> (TopologyBlock, CoordinateBlock) {
    let parsed = read_mol_block_detached(&v3000_tetrahedral_wedge(z, cfg)).unwrap();
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        ..
    } = finish_mol_block_record(parsed, true, unsanitized()).unwrap()
    else {
        panic!("concrete record expected")
    };
    (topology, coordinates)
}

#[test]
fn mol_post_false_flag_xyz_stereo_preserves_carrier_bits_and_uses_effective_2d_branch() {
    let (xy_topology, xy_coordinates) = finish_v3000_wedge("0", 1);
    assert_eq!(xy_topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
    assert_eq!(xy_coordinates.conformers_2d.len(), 1);
    assert!(xy_coordinates.conformers_3d.is_empty());

    for (z_text, expected_bits) in [
        ("0.0005", 0.0005_f64.to_bits()),
        ("-0", (-0.0_f64).to_bits()),
    ] {
        let (topology, coordinates) = finish_v3000_wedge(z_text, 1);
        assert_eq!(topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert!(coordinates.conformers_2d.is_empty());
        assert_eq!(coordinates.conformers_3d.len(), 1);
        assert!(!coordinates.conformers_3d[0].is_3d());
        assert_eq!(
            coordinates.conformers_3d[0].coordinates()[0][2].to_bits(),
            expected_bits
        );
        assert_eq!(
            coordinates.source_coordinate_dim,
            Some(CoordinateDimension::TwoD)
        );
        assert!(
            topology
                .bonds
                .iter()
                .all(|bond| bond.direction() == BondDirection::None)
        );
    }
}

#[test]
fn mol_post_final_stereochemistry_cleans_duplicate_ligands_only_in_sanitize_branch() {
    for query_center in [false, true] {
        let block = v3000_tetrahedral_case(["F", "F", "Br", "I"], query_center, "0.0005", 1, true);
        for (sanitize, remove_hs, expected_tag, expected_groups) in [
            (false, false, ChiralTag::TetrahedralCcw, 1),
            (false, true, ChiralTag::TetrahedralCcw, 1),
            (true, false, ChiralTag::Unspecified, 0),
            (true, true, ChiralTag::Unspecified, 0),
        ] {
            let (tag, groups, coordinates) =
                finished_center_tag_and_group_count(&block, sanitize, remove_hs);
            assert_eq!(
                tag, expected_tag,
                "query_center={query_center}, sanitize={sanitize}, remove_hs={remove_hs}"
            );
            assert_eq!(
                groups, expected_groups,
                "query_center={query_center}, sanitize={sanitize}, remove_hs={remove_hs}"
            );
            assert!(coordinates.conformers_2d.is_empty());
            assert_eq!(coordinates.conformers_3d.len(), 1);
            assert!(!coordinates.conformers_3d[0].is_3d());
            assert_eq!(
                coordinates.conformers_3d[0].coordinates()[0][2].to_bits(),
                0.0005_f64.to_bits()
            );
        }
    }
}

#[test]
fn mol_post_final_stereochemistry_retains_valid_false_flag_xyz_wedge_and_group() {
    for query_center in [false, true] {
        for (cfg, expected_tag) in [
            (1, ChiralTag::TetrahedralCcw),
            (3, ChiralTag::TetrahedralCw),
        ] {
            let block =
                v3000_tetrahedral_case(["F", "Cl", "Br", "I"], query_center, "-0", cfg, true);
            for remove_hs in [false, true] {
                let (tag, groups, coordinates) =
                    finished_center_tag_and_group_count(&block, true, remove_hs);
                assert_eq!(tag, expected_tag);
                assert_eq!(groups, 1);
                assert_eq!(
                    coordinates.conformers_3d[0].coordinates()[0][2].to_bits(),
                    (-0.0_f64).to_bits()
                );
                assert!(!coordinates.conformers_3d[0].is_3d());
            }
        }
    }
}

#[test]
fn mol_post_final_stereochemistry_assigns_coordinate_double_bond_stereo_after_sanitize() {
    for query_atom in [false, true] {
        let block = v3000_alkene_stereo_case(query_atom);
        for (sanitize, expected) in [
            (false, BondStereo::None),
            // Fixed RDKit 2026.03.1 legacy assignStereochemistry reports
            // STEREOE for this coordinate ordering; Trans is the separate
            // source enum value used by setBondStereoFromDirections.
            (true, BondStereo::E),
        ] {
            let parsed = read_mol_block_detached(&block).expect("alkene mol block");
            let finished = finish_mol_block_record(
                parsed,
                true,
                MolPostParams {
                    sanitize,
                    remove_hs: false,
                    expand_attachment_points: false,
                },
            )
            .expect("alkene mol-post");
            let stereo = match finished {
                MolBlockRecord::Concrete { topology, .. } => topology.bonds[1].stereo(),
                MolBlockRecord::Query(record) => record.query.bonds()[1].bond().stereo(),
            };
            assert_eq!(
                stereo, expected,
                "query_atom={query_atom}, sanitize={sanitize}"
            );
        }
    }
}

#[test]
fn mol_post_legacy_closure_iteratively_reranks_resolved_stereo_for_concrete_and_query() {
    for query_record in [false, true] {
        let parsed = read_mol_block_detached(&v3000_iterative_legacy_stereo_case(query_record))
            .expect("iterative legacy stereo block");
        let finished = finish_mol_block_record(
            parsed,
            true,
            MolPostParams {
                sanitize: true,
                remove_hs: false,
                expand_attachment_points: false,
            },
        )
        .expect("iterative legacy stereo mol-post");
        let (tags, codes) = match finished {
            MolBlockRecord::Concrete { topology, .. } => (
                topology
                    .atoms
                    .iter()
                    .map(Atom::chiral_tag)
                    .collect::<Vec<_>>(),
                topology
                    .atoms
                    .iter()
                    .map(|atom| atom.prop("_CIPCode").map(str::to_owned))
                    .collect::<Vec<_>>(),
            ),
            MolBlockRecord::Query(record) => (
                record
                    .query
                    .atoms()
                    .iter()
                    .map(|atom| atom.chiral_tag())
                    .collect::<Vec<_>>(),
                record
                    .query
                    .atoms()
                    .iter()
                    .map(|atom| atom.prop("_CIPCode").map(str::to_owned))
                    .collect::<Vec<_>>(),
            ),
        };
        assert_eq!(
            [tags[1], tags[3], tags[6]],
            [
                ChiralTag::TetrahedralCcw,
                ChiralTag::TetrahedralCcw,
                ChiralTag::TetrahedralCw,
            ],
            "query_record={query_record}"
        );
        // Pinned RDKit 2026.03.1 legacy mode gives Concrete [R,R,S] but
        // Query (atom 2 promoted by RBCNT=1) [R,S,S] for this exact fixture.
        let expected_codes: [Option<&str>; 3] = if query_record {
            [Some("R"), Some("S"), Some("S")]
        } else {
            [Some("R"), Some("R"), Some("S")]
        };
        assert_eq!(
            [
                codes[1].as_deref(),
                codes[3].as_deref(),
                codes[6].as_deref()
            ],
            expected_codes,
            "query_record={query_record}"
        );
    }
}

#[test]
fn mol_post_legacy_closure_retains_source_ring_special_cases_for_concrete_and_query() {
    for query_record in [false, true] {
        let parsed = read_mol_block_detached(&v3000_ring_special_case(query_record))
            .expect("ring-special legacy stereo block");
        let finished = finish_mol_block_record(
            parsed,
            true,
            MolPostParams {
                sanitize: true,
                remove_hs: false,
                expand_attachment_points: false,
            },
        )
        .expect("ring-special legacy stereo mol-post");
        let atoms = match &finished {
            MolBlockRecord::Concrete { topology, .. } => &topology.atoms,
            MolBlockRecord::Query(record) => {
                assert_eq!(record.query.atoms().len(), 8);
                assert_eq!(
                    record.query.atoms()[1].chiral_tag(),
                    ChiralTag::TetrahedralCw
                );
                assert_eq!(
                    record.query.atoms()[4].chiral_tag(),
                    ChiralTag::TetrahedralCw
                );
                assert_eq!(record.query.atoms()[1].prop("_ringStereoAtoms"), Some("5"));
                assert_eq!(record.query.atoms()[4].prop("_ringStereoAtoms"), Some("2"));
                continue;
            }
        };
        assert_eq!(atoms[1].chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(atoms[4].chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(atoms[1].prop("_ringStereoAtoms"), Some("5"));
        assert_eq!(atoms[4].prop("_ringStereoAtoms"), Some("2"));
    }
}

#[test]
fn mol_post_legacy_closure_sets_stereochem_done_only_in_sanitize_branch() {
    for query_record in [false, true] {
        for sanitize in [false, true] {
            let parsed = read_mol_block_detached(&v3000_iterative_legacy_stereo_case(query_record))
                .expect("stereochem-done block");
            let finished = finish_mol_block_record(
                parsed,
                true,
                MolPostParams {
                    sanitize,
                    remove_hs: false,
                    expand_attachment_points: false,
                },
            )
            .expect("stereochem-done mol-post");
            let properties = match &finished {
                MolBlockRecord::Concrete { properties, .. } => properties,
                MolBlockRecord::Query(record) => &record.properties,
            };
            assert_eq!(
                properties.prop("_StereochemDone"),
                sanitize.then_some("1"),
                "query_record={query_record}, sanitize={sanitize}"
            );
            assert_eq!(properties.is_prop_computed("_StereochemDone"), sanitize);
        }
    }
}

#[test]
fn mol_post_legacy_closure_cleans_directional_state_for_concrete_and_query() {
    for query in [false, true] {
        let input = topology(
            (0..4).map(|index| atom(index, Element::C)).collect(),
            vec![
                Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double),
                ),
                Bond::from_spec(
                    BondId::new(1),
                    BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single)
                        .with_direction(BondDirection::EndUpRight),
                ),
                Bond::from_spec(
                    BondId::new(2),
                    BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
                ),
            ],
            vec![],
        );
        let finished = finish_mol_block_record(
            concrete_or_explicit_query(input, query),
            false,
            MolPostParams {
                sanitize: true,
                remove_hs: false,
                expand_attachment_points: false,
            },
        )
        .unwrap();
        match finished {
            MolBlockRecord::Concrete { topology, .. } => {
                assert_eq!(topology.bonds[0].stereo(), BondStereo::None);
                assert_eq!(topology.bonds[1].direction(), BondDirection::None);
            }
            MolBlockRecord::Query(record) => {
                assert_eq!(
                    record.query.atoms()[0].predicate(),
                    &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
                );
                assert_eq!(record.query.bonds()[0].bond().stereo(), BondStereo::None);
                assert_eq!(
                    record.query.bonds()[1].bond().direction(),
                    BondDirection::None
                );
                assert_eq!(
                    record.query.bonds()[1].predicate(),
                    &QueryNode::predicate(BondQueryPredicate::Any)
                );
            }
        }
    }
}

#[test]
fn mol_post_legacy_closure_cleans_atrop_groups_before_general_groups_for_both_records() {
    for query in [false, true] {
        let input = TopologyBlock::try_from_parts(
            vec![
                atom(0, Element::C),
                atom(1, Element::C),
                atom(2, Element::C),
                atom(3, Element::C),
            ],
            vec![
                Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
                ),
                Bond::from_spec(
                    BondId::new(1),
                    BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single)
                        .with_stereo(BondStereo::AtropCcw),
                ),
                Bond::from_spec(
                    BondId::new(2),
                    BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Double),
                ),
            ],
            vec![],
            vec![
                StereoGroup::new(
                    StereoGroupKind::Absolute,
                    vec![AtomId::new(1), AtomId::new(2)],
                    vec![],
                )
                .with_id(17),
            ],
        )
        .unwrap();
        let finished = finish_mol_block_record(
            concrete_or_explicit_query(input, query),
            false,
            MolPostParams {
                sanitize: true,
                remove_hs: false,
                expand_attachment_points: false,
            },
        )
        .unwrap();
        let groups = match &finished {
            MolBlockRecord::Concrete { topology, .. } => &topology.stereo_groups,
            MolBlockRecord::Query(record) => {
                assert_eq!(
                    record.query.atoms()[0].predicate(),
                    &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
                );
                assert_eq!(
                    record.query.bonds()[0].predicate(),
                    &QueryNode::predicate(BondQueryPredicate::Any)
                );
                record.query.stereo_groups()
            }
        };
        assert_eq!(groups.len(), 1, "query={query}");
        assert_eq!(groups[0].atoms(), &[], "query={query}");
        assert_eq!(groups[0].bonds(), &[BondId::new(1)], "query={query}");
        assert_eq!(groups[0].id(), None, "query={query}");
    }
}

#[test]
fn mol_post_false_flag_xyz_stereo_preserves_wedge_dash_orientation_before_clearing() {
    for (cfg, expected) in [
        (1, ChiralTag::TetrahedralCcw),
        (3, ChiralTag::TetrahedralCw),
    ] {
        let (topology, coordinates) = finish_v3000_wedge("0.0005", cfg);
        assert_eq!(topology.atoms[0].chiral_tag(), expected);
        assert!(!coordinates.conformers_3d[0].is_3d());
        assert!(
            topology
                .bonds
                .iter()
                .all(|bond| bond.direction() == BondDirection::None)
        );
    }
}

#[test]
fn mol_post_false_flag_xyz_stereo_feeds_atropisomer_owner_before_direction_clearing() {
    let atoms = (0..4)
        .map(|index| {
            let hybridization = if index == 1 || index == 2 {
                Hybridization::Sp2
            } else {
                Hybridization::Sp3
            };
            Atom::from_spec(
                AtomId::new(index),
                AtomSpec::new(Element::C).with_hybridization(hybridization),
            )
        })
        .collect();
    let bonds = vec![
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single)
                .with_direction(BondDirection::BeginWedge),
        ),
        Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
        ),
        Bond::from_spec(
            BondId::new(2),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
        ),
    ];
    let coordinates = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(
            7,
            vec![
                [0.0, 1.0, -0.0],
                [0.0, 0.0, 0.0005],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
            ],
            false,
        )],
        source_coordinate_dim: Some(CoordinateDimension::TwoD),
        ..CoordinateBlock::default()
    };
    let record = MolBlockRecord::Concrete {
        topology: topology(atoms, bonds, vec![]),
        coordinates: coordinates.clone(),
        properties: MoleculeProperties::default(),
    };
    let MolBlockRecord::Concrete {
        topology,
        coordinates: output_coordinates,
        ..
    } = finish_mol_block_record(record, false, unsanitized()).unwrap()
    else {
        panic!("concrete record expected")
    };
    assert_eq!(topology.bonds[1].stereo(), BondStereo::AtropCcw);
    assert_eq!(topology.bonds[0].direction(), BondDirection::None);
    assert_eq!(output_coordinates, coordinates);
}

#[test]
fn mol_post_query_closure_false_flag_xyz_runs_stereo_before_direction_clearing() {
    let block = v3000_tetrahedral_wedge("0.0005", 1).replace("0 0.0005 0", "0 0.0005 0 RBCNT=1");
    let parsed = read_mol_block_detached(&block).expect("query wedge record");
    let MolBlockRecord::Query(record) =
        finish_mol_block_record(parsed, true, unsanitized()).expect("query mol-post")
    else {
        panic!("RBCNT must retain a query record")
    };
    assert_eq!(
        record.query.atoms()[0].chiral_tag(),
        ChiralTag::TetrahedralCcw
    );
    assert_eq!(
        record.query.conformers_3d()[0].coordinates()[0][2].to_bits(),
        0.0005_f64.to_bits()
    );
    assert!(!record.query.conformers_3d()[0].is_3d());
    assert!(
        record
            .query
            .bonds()
            .iter()
            .all(|bond| bond.bond().direction() == BondDirection::None)
    );
}

#[test]
fn q05_query_identity_composition_parsed_hydrogen_distinguishes_explicit_and_carrier_rows() {
    for sanitize in [false, true] {
        for remove_hs in [false, true] {
            for (query_hydrogen, expected_atoms) in [
                (false, if sanitize && remove_hs { 1 } else { 2 }),
                (true, 2),
            ] {
                let parsed = read_mol_block_detached(&v3000_query_hydrogen(query_hydrogen, None))
                    .expect("query hydrogen record");
                let MolBlockRecord::Query(record) = finish_mol_block_record(
                    parsed,
                    false,
                    MolPostParams {
                        sanitize,
                        remove_hs,
                        expand_attachment_points: false,
                    },
                )
                .expect("query hydrogen mol-post") else {
                    panic!("query record expected")
                };
                assert_eq!(
                    record.query.num_atoms(),
                    expected_atoms,
                    "sanitize={sanitize}, remove_hs={remove_hs}, query_hydrogen={query_hydrogen}"
                );
                assert_eq!(
                    record
                        .query
                        .coordinates_2d()
                        .expect("source 2D coordinates")
                        .len(),
                    expected_atoms
                );
            }
        }
    }
}

#[test]
fn q05_query_identity_composition_mol_post_parameter_matrix_preserves_state_and_applies_hydrogen_removal_only_after_sanitize()
 {
    for query_record in [false, true] {
        let block = if query_record {
            v3000_query_hydrogen(false, None)
        } else {
            v3000_query_hydrogen(false, None).replace(" RBCNT=1", "")
        };
        let mut source = read_mol_block_detached(&block).expect("parameter-matrix record");
        let retained_group = SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::Generic("SUP".to_owned()),
        )
        .with_atoms(vec![AtomId::new(0), AtomId::new(1)]);
        match &mut source {
            MolBlockRecord::Concrete {
                topology,
                properties,
                ..
            } => {
                topology.substance_groups = vec![retained_group.clone()];
                *properties = properties.clone().with_prop("retained", "matrix").unwrap();
            }
            MolBlockRecord::Query(record) => {
                record.substance_groups = vec![retained_group.clone()];
                record.properties = record
                    .properties
                    .clone()
                    .with_prop("retained", "matrix")
                    .unwrap();
            }
        }
        let snapshot = source.clone();

        for sanitize in [false, true] {
            for remove_hs in [false, true] {
                let expected_atoms = if sanitize && remove_hs { 1 } else { 2 };
                let finished = finish_mol_block_record(
                    source.clone(),
                    true,
                    MolPostParams {
                        sanitize,
                        remove_hs,
                        expand_attachment_points: false,
                    },
                )
                .expect("parameter-matrix mol-post");
                let expected_group_atoms = if expected_atoms == 1 {
                    vec![AtomId::new(0)]
                } else {
                    vec![AtomId::new(0), AtomId::new(1)]
                };

                match finished {
                    MolBlockRecord::Concrete {
                        topology,
                        coordinates,
                        properties,
                    } => {
                        assert!(!query_record);
                        assert_eq!(topology.atoms.len(), expected_atoms);
                        assert_eq!(
                            coordinates.conformers_2d[0].coordinates().len(),
                            expected_atoms
                        );
                        assert_eq!(properties.prop("retained"), Some("matrix"));
                        assert_eq!(topology.substance_groups.len(), 1);
                        assert_eq!(topology.substance_groups[0].atoms(), expected_group_atoms);
                    }
                    MolBlockRecord::Query(record) => {
                        assert!(query_record);
                        assert_eq!(record.query.num_atoms(), expected_atoms);
                        assert_eq!(record.query.coordinates_2d().unwrap().len(), expected_atoms);
                        assert_eq!(record.properties.prop("retained"), Some("matrix"));
                        assert_eq!(record.substance_groups.len(), 1);
                        assert_eq!(record.substance_groups[0].atoms(), expected_group_atoms);
                    }
                }
                assert_eq!(
                    source, snapshot,
                    "source changed for sanitize={sanitize}, remove_hs={remove_hs}, query_record={query_record}"
                );
            }
        }
    }
}

#[test]
fn q05_query_identity_composition_ordinary_hydrogen_removal_remaps_typed_sgroups() {
    let parsed = read_mol_block_detached(&v3000_query_hydrogen(false, None))
        .expect("query molecule with concrete hydrogen");
    let MolBlockRecord::Query(mut record) = parsed else {
        panic!("query record expected")
    };
    record.substance_groups = vec![
        SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::Generic("SUP".to_owned()),
        )
        .with_atoms(vec![AtomId::new(0), AtomId::new(1)]),
    ];
    let MolBlockRecord::Query(record) = finish_mol_block_record(
        MolBlockRecord::Query(record),
        false,
        MolPostParams {
            sanitize: true,
            remove_hs: true,
            expand_attachment_points: false,
        },
    )
    .expect("query hydrogen removal with typed SGroup") else {
        panic!("query record expected")
    };
    assert_eq!(record.query.num_atoms(), 1);
    assert_eq!(record.substance_groups.len(), 1);
    assert_eq!(record.substance_groups[0].atoms(), [AtomId::new(0)]);
}

#[test]
fn mol_post_query_closure_substitution_count_converts_and_composes_in_source_order() {
    for (value, expected) in [
        (-1, AtomQueryPredicate::ExplicitDegree(0)),
        (-2, AtomQueryPredicate::ExplicitDegree(1)),
        (1, AtomQueryPredicate::ExplicitDegree(1)),
        (6, AtomQueryPredicate::ExplicitDegreeLessEqual(6)),
    ] {
        let parsed = read_mol_block_detached(&v3000_query_hydrogen(false, Some(value)))
            .expect("SUBST query record");
        let MolBlockRecord::Query(record) =
            finish_mol_block_record(parsed, false, unsanitized()).expect("SUBST mol-post")
        else {
            panic!("SUBST must produce a query record")
        };
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    QueryNode::predicate(AtomQueryPredicate::RingBondCount(1)),
                ]),
                QueryNode::predicate(expected),
            ]),
            "SUBST={value}"
        );
    }
}

#[test]
fn mol_post_query_closure_failure_is_atomic_for_query_records() {
    let group = dat_group(
        0,
        Some("HYD"),
        None,
        None,
        vec![AtomId::new(0)],
        vec![],
        &["256"],
    );
    let query = QueryGraph::from_parts(
        vec![QueryAtom::from_parts(
            atom(0, Element::C),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        )],
        vec![],
        Default::default(),
        vec![],
        vec![],
        vec![],
    )
    .unwrap();
    let source = MolBlockRecord::Query(QueryMolBlockRecord {
        query,
        substance_groups: vec![group],
        properties: MoleculeProperties::default(),
        source_coordinate_dim: None,
    });
    let snapshot = source.clone();
    assert_eq!(
        finish_mol_block_record(source, false, unsanitized()),
        Err(MolPostError::Representation("HYD count outside u8"))
    );
    let MolBlockRecord::Query(snapshot) = snapshot else {
        unreachable!()
    };
    assert_eq!(snapshot.query.atoms()[0].explicit_hydrogens(), 0);
    assert_eq!(snapshot.substance_groups.len(), 1);
}

#[test]
fn mol_post_query_closure_atom_and_dat_queries_follow_source_order_then_complete_scan() {
    let atom = Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C)
            .with_prop("molSubstCount", "1")
            .unwrap(),
    );
    let smart = dat_group(
        0,
        None,
        Some("SMARTSQ"),
        None,
        vec![AtomId::new(0)],
        vec![],
        &["[#7]"],
    );
    let charge = dat_group(
        1,
        Some("ZCH"),
        None,
        None,
        vec![AtomId::new(0)],
        vec![],
        &["-1"],
    );
    let record = MolBlockRecord::Concrete {
        topology: topology(vec![atom], vec![], vec![smart, charge]),
        coordinates: CoordinateBlock::default(),
        properties: MoleculeProperties::default()
            .with_prop("_NeedsQueryScan", "1")
            .unwrap(),
    };
    let MolBlockRecord::Query(record) =
        finish_mol_block_record(record, false, unsanitized()).expect("ordered query closure")
    else {
        panic!("SMARTSQ must produce a query record")
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
    );
    assert_eq!(record.query.atoms()[0].formal_charge(), -1);
    assert_eq!(record.query.prop("_NeedsQueryScan"), None);
    assert!(record.substance_groups.is_empty());
}

#[test]
fn mol_post_params_match_source_defaults_and_noop_expansion() {
    assert_eq!(
        MolPostParams::default(),
        MolPostParams {
            sanitize: true,
            remove_hs: true,
            expand_attachment_points: false,
        }
    );
    let record = concrete(topology(vec![atom(0, Element::C)], vec![], vec![]));
    let source = record.clone();
    assert!(matches!(
        finish_mol_block_record(
            record,
            false,
            MolPostParams {
                expand_attachment_points: true,
                ..MolPostParams::default()
            }
        ),
        Ok(MolBlockRecord::Concrete { .. })
    ));
    assert_eq!(
        source,
        concrete(topology(vec![atom(0, Element::C)], vec![], vec![]))
    );
}

#[test]
fn mol_post_mol_tot_valence_sentinels_clear_property_and_respect_zbo_h() {
    for value in ["15", "-1"] {
        let spec = AtomSpec::new(Element::C)
            .with_explicit_hydrogens(3)
            .with_prop("molTotValence", value)
            .unwrap();
        let record = concrete(topology(
            vec![Atom::from_spec(AtomId::new(0), spec)],
            vec![],
            vec![],
        ));
        let MolBlockRecord::Concrete { topology, .. } =
            finish_mol_block_record(record, false, unsanitized()).unwrap()
        else {
            panic!("concrete record expected")
        };
        assert_eq!(topology.atoms[0].explicit_hydrogens(), 0);
        assert!(topology.atoms[0].no_implicit());
        assert_eq!(topology.atoms[0].prop("molTotValence"), None);
    }

    let spec = AtomSpec::new(Element::C)
        .with_explicit_hydrogens(2)
        .with_prop("molTotValence", "4")
        .unwrap()
        .with_prop("_ZBO_H", "1")
        .unwrap();
    let record = concrete(topology(
        vec![Atom::from_spec(AtomId::new(0), spec)],
        vec![],
        vec![],
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        finish_mol_block_record(record, false, unsanitized()).unwrap()
    else {
        panic!("concrete record expected")
    };
    assert_eq!(topology.atoms[0].explicit_hydrogens(), 2);
    assert!(!topology.atoms[0].no_implicit());
    assert_eq!(topology.atoms[0].prop("molTotValence"), None);
}

#[test]
fn mol_post_processes_atom_properties_before_hyd_group() {
    // ProcessMolProps applies molTotValence before processSGroups. HYD then
    // overwrites the H count but does not undo noImplicit set by the earlier
    // valence branch.
    let spec = AtomSpec::new(Element::C)
        .with_prop("molTotValence", "4")
        .unwrap();
    let hyd = dat_group(
        0,
        Some("HYD"),
        None,
        None,
        vec![AtomId::new(0)],
        vec![],
        &["2"],
    );
    let record = concrete(topology(
        vec![Atom::from_spec(AtomId::new(0), spec)],
        vec![],
        vec![hyd],
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        finish_mol_block_record(record, false, unsanitized()).unwrap()
    else {
        panic!("concrete record expected")
    };
    assert_eq!(topology.atoms[0].explicit_hydrogens(), 2);
    assert!(topology.atoms[0].no_implicit());
    assert_eq!(topology.atoms[0].prop("_ZBO_H"), Some("1"));
    assert_eq!(topology.atoms[0].prop("molTotValence"), None);
    assert!(topology.substance_groups.is_empty());
}

#[test]
fn mol_post_dat_actions_remove_only_recognized_groups() {
    let atoms = vec![atom(0, Element::C), atom(1, Element::N)];
    let bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified)
            .with_prop("kept", "yes")
            .unwrap(),
    );
    let groups = vec![
        dat_group(
            0,
            Some("MRV_COORDINATE_BOND_TYPE"),
            None,
            None,
            vec![],
            vec![],
            &["1", "ignored-extra"],
        ),
        dat_group(
            1,
            Some("ZCH"),
            None,
            None,
            vec![AtomId::new(0), AtomId::new(1)],
            vec![],
            &["-1;2"],
        ),
        dat_group(
            2,
            Some("UNKNOWN"),
            None,
            None,
            vec![AtomId::new(0)],
            vec![],
            &["kept"],
        ),
    ];
    let record = concrete(topology(atoms, vec![bond], groups));
    let MolBlockRecord::Concrete { topology, .. } =
        finish_mol_block_record(record, false, unsanitized()).unwrap()
    else {
        panic!("concrete record expected")
    };
    assert_eq!(topology.bonds[0].order(), BondOrder::Dative);
    assert_eq!(topology.bonds[0].prop("kept"), Some("yes"));
    assert_eq!(topology.atoms[0].formal_charge(), -1);
    assert_eq!(topology.atoms[1].formal_charge(), 2);
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(topology.substance_groups[0].id(), SubstanceGroupId::new(0));
}

#[test]
fn mol_post_short_zch_and_hyd_fields_are_source_noops_but_groups_are_consumed() {
    for field in ["ZCH", "HYD"] {
        let group = dat_group(
            0,
            Some(field),
            None,
            None,
            vec![AtomId::new(0), AtomId::new(1)],
            vec![],
            &["3"],
        );
        let record = concrete(topology(
            vec![atom(0, Element::C), atom(1, Element::N)],
            vec![],
            vec![group],
        ));
        let MolBlockRecord::Concrete { topology, .. } =
            finish_mol_block_record(record, false, unsanitized()).unwrap()
        else {
            panic!("concrete record expected")
        };
        assert_eq!(topology.atoms[0].formal_charge(), 0);
        assert_eq!(topology.atoms[1].formal_charge(), 0);
        assert_eq!(topology.atoms[0].explicit_hydrogens(), 0);
        assert_eq!(topology.atoms[1].explicit_hydrogens(), 0);
        assert!(topology.substance_groups.is_empty());
    }
}

#[test]
fn mol_post_mrv_implicit_h_requires_an_aromatic_bond() {
    let aromatic = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Aromatic).with_aromatic(true),
    );
    let group = dat_group(
        0,
        Some("MRV_IMPLICIT_H"),
        None,
        None,
        vec![AtomId::new(0), AtomId::new(2)],
        vec![],
        &["IMPL_H2", "OTHER9"],
    );
    let record = concrete(topology(
        vec![
            atom(0, Element::C),
            atom(1, Element::C),
            atom(2, Element::N),
        ],
        vec![aromatic],
        vec![group],
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        finish_mol_block_record(record, false, unsanitized()).unwrap()
    else {
        panic!("concrete record expected")
    };
    assert_eq!(topology.atoms[0].explicit_hydrogens(), 2);
    assert_eq!(topology.atoms[2].explicit_hydrogens(), 0);
    assert!(topology.substance_groups.is_empty());
}

#[test]
fn mol_post_smartsq_builds_typed_query_and_preserves_unconsumed_groups() {
    let smart = dat_group(
        0,
        None,
        Some("SMARTSQ"),
        None,
        vec![AtomId::new(0)],
        vec![],
        &["[#7]"],
    );
    let unknown = dat_group(
        1,
        None,
        Some("OTHER"),
        None,
        vec![AtomId::new(0)],
        vec![],
        &["[#8]"],
    );
    let record = concrete(topology(
        vec![atom(0, Element::C)],
        vec![],
        vec![smart, unknown],
    ));
    let MolBlockRecord::Query(record) =
        finish_mol_block_record(record, false, unsanitized()).unwrap()
    else {
        panic!("query record expected")
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
    );
    assert_eq!(record.query.atoms()[0].prop("MRV SMA"), Some("[#7]"));
    assert_eq!(record.query.atoms()[0].prop("_MolFileAtomQuery"), Some("1"));
    assert_eq!(record.substance_groups.len(), 1);
    assert_eq!(record.substance_groups[0].id(), SubstanceGroupId::new(0));
}

#[test]
fn mol_post_smartsq_non_equals_is_consumed_without_replacing_predicate() {
    let smart = dat_group(
        0,
        None,
        Some("SQ"),
        Some("!="),
        vec![AtomId::new(0)],
        vec![],
        &["[#7]"],
    );
    let record = concrete(topology(vec![atom(0, Element::C)], vec![], vec![smart]));
    let MolBlockRecord::Query(record) =
        finish_mol_block_record(record, false, unsanitized()).unwrap()
    else {
        panic!("query record expected")
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6))
    );
    assert!(record.substance_groups.is_empty());
}

#[test]
fn mol_post_query_scan_is_cleared_and_completed() {
    let base = atom(0, Element::C);
    let query_atom = QueryAtom::from_parts(
        base,
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    let query = QueryGraph::from_parts(
        vec![query_atom],
        Vec::<QueryBond>::new(),
        [("_NeedsQueryScan".to_owned(), "1".to_owned())]
            .into_iter()
            .collect(),
        vec![],
        vec![],
        vec![],
    )
    .unwrap();
    let record = MolBlockRecord::Query(QueryMolBlockRecord {
        query,
        substance_groups: vec![],
        properties: MoleculeProperties::default(),
        source_coordinate_dim: None,
    });
    let MolBlockRecord::Query(record) =
        finish_mol_block_record(record, false, unsanitized()).unwrap()
    else {
        panic!("query record expected")
    };
    assert_eq!(record.query.prop("_NeedsQueryScan"), None);
    assert_eq!(record.query.num_atoms(), 1);
}

fn explicit_query_record(atoms: Vec<QueryAtom>, bonds: Vec<QueryBond>) -> MolBlockRecord {
    MolBlockRecord::Query(QueryMolBlockRecord {
        query: QueryGraph::from_parts(atoms, bonds, Default::default(), vec![], vec![], vec![])
            .unwrap(),
        substance_groups: vec![],
        properties: MoleculeProperties::default(),
        source_coordinate_dim: None,
    })
}

#[test]
fn mol_post_explicit_query_provenance_preserves_unmarked_atom_predicates_and_matching() {
    let predicates = [
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
        QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
            QueryNode::predicate(AtomQueryPredicate::FormalCharge(1)),
        ]),
        QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(
            RecursiveStructureQuery::from_query_graph(
                QueryGraph::from_parts(
                    vec![QueryAtom::from_parts(
                        atom(0, Element::N),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
                    )],
                    vec![],
                    Default::default(),
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
                17,
            )
            .with_source_smarts("[$([#7])]"),
        )),
    ];

    for predicate in predicates {
        for sanitize in [false, true] {
            for remove_hs in [false, true] {
                let record = explicit_query_record(
                    vec![QueryAtom::from_parts(
                        atom(0, Element::C),
                        predicate.clone(),
                    )],
                    vec![],
                );
                let MolBlockRecord::Query(finished) = finish_mol_block_record(
                    record,
                    false,
                    MolPostParams {
                        sanitize,
                        remove_hs,
                        expand_attachment_points: false,
                    },
                )
                .unwrap() else {
                    panic!("explicit query record expected")
                };
                assert_eq!(
                    finished.query.atoms()[0].predicate(),
                    &predicate,
                    "sanitize={sanitize}, remove_hs={remove_hs}"
                );
            }
        }
    }

    let record = explicit_query_record(
        vec![QueryAtom::from_parts(
            atom(0, Element::C),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
        )],
        vec![],
    );
    let MolBlockRecord::Query(finished) =
        finish_mol_block_record(record, false, unsanitized()).unwrap()
    else {
        panic!("explicit query record expected")
    };
    let nitrogen = topology(vec![atom(0, Element::N)], vec![], vec![]);
    let carbon = topology(vec![atom(0, Element::C)], vec![], vec![]);
    assert_eq!(
        cosmolkit_search::QueryGraphOperator::new(&finished.query)
            .matches(&nitrogen)
            .unwrap()
            .len(),
        1
    );
    assert!(
        cosmolkit_search::QueryGraphOperator::new(&finished.query)
            .matches(&carbon)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn mol_post_explicit_query_provenance_preserves_unmarked_bond_predicates() {
    let predicates = [
        QueryNode::predicate(BondQueryPredicate::Any),
        QueryNode::or(vec![
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
        ]),
        QueryNode::and(vec![
            QueryNode::predicate(BondQueryPredicate::Any),
            QueryNode::not(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
        ]),
    ];
    for predicate in predicates {
        for sanitize in [false, true] {
            let raw_bond = Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            );
            let record = explicit_query_record(
                vec![
                    QueryAtom::from_parts(
                        atom(0, Element::C),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    ),
                    QueryAtom::from_parts(
                        atom(1, Element::C),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    ),
                ],
                vec![QueryBond::from_parts(raw_bond, predicate.clone())],
            );
            let MolBlockRecord::Query(finished) = finish_mol_block_record(
                record,
                false,
                MolPostParams {
                    sanitize,
                    remove_hs: false,
                    expand_attachment_points: false,
                },
            )
            .unwrap() else {
                panic!("explicit query record expected")
            };
            assert_eq!(
                finished.query.bonds()[0].predicate(),
                &predicate,
                "sanitize={sanitize}"
            );
        }
    }
}

#[test]
fn q05_query_identity_composition_explicit_query_provenance_survives_remapping_and_errors_are_atomic()
 {
    let atoms = vec![
        QueryAtom::from_parts(
            atom(0, Element::C),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
        ),
        QueryAtom::from_parts(
            atom(1, Element::H),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
        ),
    ];
    let bond = QueryBond::from_parts(
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        ),
        QueryNode::predicate(BondQueryPredicate::Any),
    );
    let source = explicit_query_record(atoms, vec![bond]);
    for sanitize in [false, true] {
        for remove_hs in [false, true] {
            let MolBlockRecord::Query(finished) = finish_mol_block_record(
                source.clone(),
                false,
                MolPostParams {
                    sanitize,
                    remove_hs,
                    expand_attachment_points: false,
                },
            )
            .unwrap() else {
                panic!("explicit query record expected")
            };
            assert_eq!(
                finished.query.num_atoms(),
                2,
                "sanitize={sanitize}, remove_hs={remove_hs}"
            );
            assert_eq!(
                finished.query.atoms()[0].predicate(),
                &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
            );
            assert_eq!(
                finished.query.atoms()[1].predicate(),
                &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1))
            );
            assert_eq!(
                finished.query.bonds()[0].predicate(),
                &QueryNode::predicate(BondQueryPredicate::Any)
            );
        }
    }

    assert!(matches!(
        finish_mol_block_record(
            source.clone(),
            false,
            MolPostParams {
                sanitize: false,
                remove_hs: false,
                expand_attachment_points: true,
            },
        ),
        Ok(MolBlockRecord::Query(_))
    ));
    let MolBlockRecord::Query(source_record) = source else {
        unreachable!()
    };
    assert_eq!(
        source_record.query.atoms()[0].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
    );
    assert_eq!(
        source_record.query.bonds()[0].predicate(),
        &QueryNode::predicate(BondQueryPredicate::Any)
    );
}

#[test]
fn mol_post_preserves_input_value_on_processing_error() {
    let group = dat_group(
        0,
        Some("HYD"),
        None,
        None,
        vec![AtomId::new(0)],
        vec![],
        &["256"],
    );
    let record = concrete(topology(vec![atom(0, Element::C)], vec![], vec![group]));
    let source = record.clone();
    assert_eq!(
        finish_mol_block_record(record, false, unsanitized()),
        Err(MolPostError::Representation("HYD count outside u8"))
    );
    let MolBlockRecord::Concrete { topology, .. } = source else {
        panic!("concrete record expected")
    };
    assert_eq!(topology.atoms[0].explicit_hydrogens(), 0);
    assert_eq!(topology.substance_groups.len(), 1);
}

#[test]
fn mol_post_query_bond_state_survives_dat_processing() {
    let atoms = vec![
        QueryAtom::from_parts(
            atom(0, Element::C),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        ),
        QueryAtom::from_parts(
            atom(1, Element::N),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
        ),
    ];
    let raw_bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified)
            .with_prop("_MolFileBondQuery", "1")
            .unwrap(),
    );
    let bonds = vec![QueryBond::from_parts(
        raw_bond,
        QueryNode::predicate(BondQueryPredicate::Any),
    )];
    let query =
        QueryGraph::from_parts(atoms, bonds, Default::default(), vec![], vec![], vec![]).unwrap();
    let group = dat_group(
        0,
        Some("MRV_COORDINATE_BOND_TYPE"),
        None,
        None,
        vec![],
        vec![],
        &["1"],
    );
    let record = MolBlockRecord::Query(QueryMolBlockRecord {
        query,
        substance_groups: vec![group],
        properties: MoleculeProperties::default(),
        source_coordinate_dim: None,
    });
    let MolBlockRecord::Query(record) =
        finish_mol_block_record(record, false, unsanitized()).unwrap()
    else {
        panic!("query record expected")
    };
    assert_eq!(record.query.bonds()[0].bond().order(), BondOrder::Dative);
    assert_eq!(
        record.query.bonds()[0].predicate(),
        &QueryNode::predicate(BondQueryPredicate::Any)
    );
    assert!(record.substance_groups.is_empty());
}

#[test]
fn mol_post_query_predicate_sync_rebuilds_synthesized_aromatic_carriers() {
    let parsed = read_mol_block_detached(&v3000_benzene_query_case()).expect("benzene query");
    let MolBlockRecord::Query(record) = finish_mol_block_record(
        parsed,
        true,
        MolPostParams {
            sanitize: true,
            remove_hs: false,
            expand_attachment_points: false,
        },
    )
    .expect("benzene query mol-post") else {
        panic!("RBCNT must retain a query record")
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)),
        ])
    );
    for atom in &record.query.atoms()[1..] {
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6))
        );
        assert!(atom.is_aromatic());
    }
    for bond in record.query.bonds() {
        assert_eq!(bond.bond().order(), BondOrder::Aromatic);
        assert!(bond.bond().is_aromatic());
        assert_eq!(
            bond.predicate(),
            &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic))
        );
    }
    let target = TopologyBlock::try_from_parts(
        record
            .query
            .atoms()
            .iter()
            .map(|atom| {
                atom.try_to_atom()
                    .expect("these concrete query carriers retain Element identity")
            })
            .collect(),
        record
            .query
            .bonds()
            .iter()
            .map(|bond| bond.bond().clone())
            .collect(),
        vec![],
        record.query.stereo_groups().to_vec(),
    )
    .expect("aromatic target topology");
    assert!(
        !cosmolkit_search::QueryGraphOperator::new(&record.query)
            .matches(&target)
            .expect("query match")
            .is_empty(),
        "final carrier predicates must match the equivalent sanitized target"
    );
}

#[test]
fn mol_post_query_predicate_sync_preserves_explicit_source_query_bond() {
    let block = "explicit query bond\n  COSMolKit\n\n  0  0  0     0  0            999 V3000\n\
M  V30 BEGIN CTAB\n\
M  V30 COUNTS 2 1 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0\n\
M  V30 2 C 1 0 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 6 1 2\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n";
    for sanitize in [false, true] {
        let parsed = read_mol_block_detached(block).expect("explicit bond query");
        let MolBlockRecord::Query(record) = finish_mol_block_record(
            parsed,
            false,
            MolPostParams {
                sanitize,
                remove_hs: false,
                expand_attachment_points: false,
            },
        )
        .expect("explicit bond query mol-post") else {
            panic!("query bond must retain query record")
        };
        assert_eq!(record.query.bonds().len(), 1);
        assert_eq!(
            record.query.bonds()[0].predicate(),
            &QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ]))
        );
        assert_eq!(
            record.query.bonds()[0].bond().prop("_MolFileBondQuery"),
            Some("1")
        );
    }
}
