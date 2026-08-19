use super::*;

use crate::{
    AtomQueryPredicate, AtomSpec, BondDirection, BondOrder, BondQueryPredicate, BondSpec,
    BondStereo, ChiralTag, Element, QueryNode, SGroupAttachPoint, SGroupBondRole, SGroupBracket,
    SGroupBracketStyle, SGroupCState, SGroupConnection, SGroupData, SGroupDisplay, StereoGroup,
    StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
};
use crate::{MOLBLOCK_IO_FEATURE, Molecule, UnsupportedFeatureError};

fn charged_isotope_molecule() -> Molecule {
    let mut builder = Molecule::builder().with_name("charged");
    builder.add_atom(
        AtomSpec::new(Element::C)
            .with_formal_charge(-1)
            .with_isotope(13)
            .with_prop("molFileValue", "payload")
            .with_prop("molFileAlias", "AliasLabel"),
    );
    builder.set_2d_coordinates(vec![[1.25, -2.5]]).unwrap();
    builder.build().unwrap()
}

#[test]
fn molblock_writer_rejects_out_of_range_effective_atomic_number_like_rdkit() {
    let molecule = Molecule::from_smiles(
        "[C+9]([F])([F])([F])([F])([F])([F])([F])([F])([F])([F])([F])([F])[F]",
    )
    .unwrap();
    let error = mol_to_v2000_block(&molecule).unwrap_err();
    assert!(matches!(
        error,
        MolWriteError::Value(ref reason) if reason == "Atomic number not found"
    ));
}

#[test]
fn non_default_valence_skips_periodic_table_lookup_when_implicit_hs_are_allowed() {
    let mut builder = Molecule::builder();
    let atom_id = builder.add_atom(AtomSpec::new(Element::C).with_formal_charge(9));
    let molecule = builder.build().unwrap();
    let atom = &molecule.atoms()[atom_id.index()];
    assert!(!atom.no_implicit());
    let valence = molblock_valence_assignment(&molecule).unwrap();
    assert!(!has_non_default_valence(&molecule, atom, &valence).unwrap());
}

fn ethene_molecule() -> Molecule {
    let mut builder = Molecule::builder().with_name("ethene");
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Double))
        .unwrap();
    builder
        .set_2d_coordinates(vec![[-0.75, 0.0], [0.75, -0.0]])
        .unwrap();
    builder.build().unwrap()
}

fn sodium_chloride_molecule() -> Molecule {
    let mut builder = Molecule::builder().with_name("salt");
    builder.add_atom(AtomSpec::new(Element::from_atomic_number(11).unwrap()).with_formal_charge(1));
    builder
        .add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()).with_formal_charge(-1));
    builder
        .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
        .unwrap();
    builder.build().unwrap()
}

fn zbo_extension_molecule() -> Molecule {
    let mut builder = Molecule::builder()
        .with_name("zbo")
        .with_property("_MolFileInfoLine", "  COSMolKit          2D");
    let carbon = builder.add_atom(
        AtomSpec::new(Element::C)
            .with_formal_charge(-1)
            .with_explicit_hydrogens(2)
            .with_prop("_ZBO_H", "1"),
    );
    let rgroup = builder.add_atom(
        AtomSpec::new(Element::DUMMY)
            .with_prop("_MolFileRLabel", "7")
            .with_prop("_MolFile_PXA", " payload"),
    );
    builder
        .add_bond(BondSpec::new(carbon, rgroup, BondOrder::Zero))
        .unwrap();
    builder
        .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
        .unwrap();
    builder.build().unwrap()
}

fn atom_list_query_molecule() -> Molecule {
    let mut builder = Molecule::builder().with_name("atom-list");
    builder.add_atom(
        AtomSpec::new(Element::DUMMY)
            .with_query(QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(
                vec![6, 7],
            )))
            .with_prop("molFileValue", "[#6,#7]"),
    );
    builder.add_atom(AtomSpec::new(Element::C).with_query(QueryNode::predicate(
        AtomQueryPredicate::AtomicNumberNotIn(vec![8, 16]),
    )));
    builder
        .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
        .unwrap();
    builder.build().unwrap()
}

fn query_bond_molecule() -> Molecule {
    let mut builder = Molecule::builder().with_name("query-bond");
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    let a4 = builder.add_atom(AtomSpec::new(Element::C));
    let a5 = builder.add_atom(AtomSpec::new(Element::C));
    let a6 = builder.add_atom(AtomSpec::new(Element::C));
    let a7 = builder.add_atom(AtomSpec::new(Element::C));
    let a8 = builder.add_atom(AtomSpec::new(Element::C));
    let a9 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(
            BondSpec::new(a0, a1, BondOrder::Unspecified)
                .with_query(QueryNode::predicate(BondQueryPredicate::Any)),
        )
        .unwrap();
    builder
        .add_bond(
            BondSpec::new(a2, a3, BondOrder::Unspecified).with_query(QueryNode::predicate(
                BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Double]),
            )),
        )
        .unwrap();
    builder
        .add_bond(
            BondSpec::new(a4, a5, BondOrder::Unspecified).with_query(QueryNode::predicate(
                BondQueryPredicate::MolFileQueryCode(42),
            )),
        )
        .unwrap();
    builder
        .add_bond(
            BondSpec::new(a6, a7, BondOrder::Single)
                .with_query(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
        )
        .unwrap();
    builder
        .add_bond(
            BondSpec::new(a8, a9, BondOrder::Unspecified).with_query(QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                    BondOrder::Single,
                    BondOrder::Aromatic,
                ])),
                QueryNode::predicate(BondQueryPredicate::IsInRing(false)),
            ])),
        )
        .unwrap();
    builder
        .set_2d_coordinates(vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [4.0, 0.0],
            [5.0, 0.0],
            [6.0, 0.0],
            [7.0, 0.0],
            [8.0, 0.0],
            [9.0, 0.0],
        ])
        .unwrap();
    builder.build().unwrap()
}

fn chiral_flag_molecule() -> Molecule {
    let mut builder = Molecule::builder()
        .with_name("chiral-flag")
        .with_property("_MolFileChiralFlag", "1");
    builder.add_atom(AtomSpec::new(Element::C));
    builder.set_2d_coordinates(vec![[0.0, 0.0]]).unwrap();
    builder.build().unwrap()
}

fn aromatic_benzene_molecule() -> Molecule {
    let mut builder = Molecule::builder().with_name("benzene");
    let atoms = (0..6)
        .map(|_| builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true)))
        .collect::<Vec<_>>();
    for idx in 0..6 {
        builder
            .add_bond(
                BondSpec::new(atoms[idx], atoms[(idx + 1) % 6], BondOrder::Aromatic)
                    .with_aromatic(true),
            )
            .unwrap();
    }
    builder
        .set_2d_coordinates(vec![
            [1.0, 0.0],
            [0.5, 0.866],
            [-0.5, 0.866],
            [-1.0, 0.0],
            [-0.5, -0.866],
            [0.5, -0.866],
        ])
        .unwrap();
    builder.build().unwrap()
}

fn chiral_state_molecule() -> Molecule {
    let mut builder = Molecule::builder().with_name("chiral-state");
    builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
    builder.set_2d_coordinates(vec![[0.0, 0.0]]).unwrap();
    builder.build().unwrap()
}

fn v3000_bond_cfg_molecule() -> Molecule {
    let mut builder = Molecule::builder().with_name("bond-cfg");
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(
            BondSpec::new(a0, a1, BondOrder::Single).with_direction(BondDirection::BeginWedge),
        )
        .unwrap();
    builder
        .add_bond(
            BondSpec::new(a1, a2, BondOrder::Double)
                .with_direction(BondDirection::EitherDouble)
                .with_stereo(BondStereo::Any),
        )
        .unwrap();
    builder
        .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
        .unwrap();
    builder.build().unwrap()
}

fn sgroup_molecule() -> Molecule {
    let mut builder = Molecule::builder().with_name("sgroup");
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::O));
    let b0 = builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
        .unwrap();
    let sup = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
        .with_external_id(7)
        .with_atoms(vec![a0, a1])
        .with_parent_atoms(vec![a0])
        .with_bonds(vec![b0])
        .with_label("Me")
        .with_subtype("ALT")
        .with_connection(SGroupConnection::HeadToTail)
        .with_expansion_state("E")
        .with_bracket_style(SGroupBracketStyle::Bracket)
        .with_display(SGroupDisplay {
            brackets: vec![SGroupBracket {
                p1: [0.0, 1.0],
                p2: [2.0, 3.0],
            }],
            ..SGroupDisplay::default()
        })
        .with_cstates(vec![SGroupCState {
            bond: b0,
            vector: [0.5, 0.25],
        }])
        .with_attach_points(vec![SGroupAttachPoint {
            atom: a0,
            leaving_atom: Some(a1),
            label: Some("AP".to_string()),
            order: None,
        }]);
    let dat = SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
        .with_parent(SubstanceGroupId::new(0))
        .with_component_number(5)
        .with_class("CLASS")
        .with_bracket_style(SGroupBracketStyle::Parenthesis)
        .with_data(SGroupData {
            field_name: Some("FIELD".to_string()),
            field_type: Some("T".to_string()),
            field_info: Some("INFO".to_string()),
            field_display: Some("display spec".to_string()),
            query_type: Some("Q".to_string()),
            query_op: Some("OP".to_string()),
            values: vec!["first valuesecond value".to_string()],
            ..SGroupData::default()
        });
    builder.add_substance_group(sup).unwrap();
    builder.add_substance_group(dat).unwrap();
    builder
        .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
        .unwrap();
    builder.build().unwrap()
}

#[test]
fn mol_to_v2000_block_writes_basic_topology_coordinates_and_properties() {
    let molecule = charged_isotope_molecule();

    let block = mol_to_v2000_block(&molecule).unwrap();

    assert!(block.starts_with("charged\n  COSMolKit          2D\n\n"));
    assert!(block.contains("  1  0  0  0  0  0  0  0  0  0999 V2000\n"));
    assert!(block.contains("    1.2500   -2.5000    0.0000 C   0"));
    assert!(block.contains("M  CHG  1   1  -1\n"));
    assert!(block.contains("M  ISO  1   1  13\n"));
    assert!(block.contains("V    1 payload\n"));
    assert!(block.contains("A    1\nAliasLabel\n"));
    assert!(block.ends_with("M  END\n"));
}

#[test]
fn mol_to_v2000_block_pads_atom_symbol_like_rdkit_atomgetmolfilesymbol() {
    let molecule = ethene_molecule();

    let block = mol_to_v2000_block(&molecule).unwrap();

    assert!(
        block.contains("   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n")
    );
    assert!(
        block.contains("    0.7500   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n")
    );
}

#[test]
fn mol_to_v2000_block_writes_isolated_ion_valence_like_rdkit_getmolfileatomproperties() {
    let molecule = sodium_chloride_molecule();

    let block = mol_to_v2000_block(&molecule).unwrap();

    assert!(
        block.contains("    0.0000    0.0000    0.0000 Na  0  0  0  0  0 15  0  0  0  0  0  0\n")
    );
    assert!(
        block.contains("    1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n")
    );
    assert!(block.contains("M  CHG  2   1   1   2  -1\n"));
}

#[test]
fn mol_to_v3000_block_writes_isolated_ion_valence_like_rdkit_getv3000molfileatomline() {
    let molecule = sodium_chloride_molecule();

    let block = mol_to_v3000_block(&molecule).unwrap();

    assert!(block.contains("M  V30 1 Na 0.000000 0.000000 0.000000 0 CHG=1 VAL=-1\n"));
    assert!(block.contains("M  V30 2 Cl 1.000000 0.000000 0.000000 0 CHG=-1\n"));
}

#[test]
fn molfile_total_valence_field_tracks_rdkit_r_dummy_after_v3000_read() {
    let input = concat!(
        "\n",
        "     RDKit          2D\n",
        "\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 COUNTS 2 1 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 R -0.750000 0.000000 0.000000 1 VAL=1\n",
        "M  V30 2 C 0.750000 -0.000000 0.000000 0\n",
        "M  V30 END ATOM\n",
        "M  V30 BEGIN BOND\n",
        "M  V30 1 1 1 2\n",
        "M  V30 END BOND\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );

    let record = crate::io::molfile::read_mol_record_from_str(input).unwrap();
    let molecule = &record.molecule;
    let atom = &molecule.atoms()[0];

    assert_eq!(atom.atomic_number(), 0);
    assert_eq!(atom.prop("dummyLabel"), Some("R"));
    assert_eq!(atom.prop("_MolFileRLabel"), None);
    assert_eq!(atom.atom_map(), Some(1));
    assert_eq!(atom.query(), None);
    assert!(atom.no_implicit());
    assert_eq!(atom.explicit_hydrogens(), 0);
    assert_eq!(molecule.topology_block().adjacency.neighbors_of(0).len(), 1);
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
    let valence = molblock_valence_assignment(molecule).unwrap();
    assert!(has_non_default_valence(molecule, atom, &valence).unwrap());
    assert_eq!(molfile_total_valence_field(molecule, atom).unwrap(), 1);
}

#[test]
fn molblock_writer_marks_rdkit_sulfur_valence_six_in_organic_subset_case() {
    let input = concat!(
        "\n",
        "     RDKit          2D\n",
        "\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 COUNTS 27 29 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 O 2.914477 -0.465275 0.000000 0\n",
        "M  V30 2 S 3.048886 1.028691 0.000000 0 VAL=6\n",
        "M  V30 3 O 3.183296 2.522657 0.000000 0\n",
        "M  V30 4 C 1.554921 1.163100 0.000000 0\n",
        "M  V30 5 C 0.691536 -0.063507 0.000000 0\n",
        "M  V30 6 C -0.802430 0.070902 0.000000 0\n",
        "M  V30 7 O -1.433011 1.431919 0.000000 0\n",
        "M  V30 8 N -1.665815 -1.155706 0.000000 0\n",
        "M  V30 9 C -3.159781 -1.021296 0.000000 0\n",
        "M  V30 10 C -3.790362 0.339721 0.000000 0\n",
        "M  V30 11 C -5.284328 0.474130 0.000000 0\n",
        "M  V30 12 C -6.147713 -0.752477 0.000000 0\n",
        "M  V30 13 N -7.641679 -0.618068 0.000000 0\n",
        "M  V30 14 C -8.272260 0.742949 0.000000 0\n",
        "M  V30 15 O -7.408875 1.969557 0.000000 0\n",
        "M  V30 16 C -9.766226 0.877359 0.000000 0\n",
        "M  V30 17 C -5.517132 -2.113494 0.000000 0\n",
        "M  V30 18 C -4.023166 -2.247904 0.000000 0\n",
        "M  V30 19 C 4.542852 0.894281 0.000000 0\n",
        "M  V30 20 C 5.406237 2.120889 0.000000 0\n",
        "M  V30 21 C 6.900203 1.986480 0.000000 0\n",
        "M  V30 22 C 7.530784 0.625463 0.000000 0\n",
        "M  V30 23 C 6.667399 -0.601145 0.000000 0\n",
        "M  V30 24 N 7.001229 -2.063526 0.000000 0\n",
        "M  V30 25 O 5.713581 -2.832918 0.000000 0\n",
        "M  V30 26 N 4.583941 -1.846047 0.000000 0\n",
        "M  V30 27 C 5.173433 -0.466736 0.000000 0\n",
        "M  V30 END ATOM\n",
        "M  V30 BEGIN BOND\n",
        "M  V30 1 2 1 2\n",
        "M  V30 2 2 2 3\n",
        "M  V30 3 1 2 4\n",
        "M  V30 4 1 4 5\n",
        "M  V30 5 1 5 6\n",
        "M  V30 6 2 6 7\n",
        "M  V30 7 1 6 8\n",
        "M  V30 8 1 8 9\n",
        "M  V30 9 1 9 10\n",
        "M  V30 10 2 10 11\n",
        "M  V30 11 1 11 12\n",
        "M  V30 12 1 12 13\n",
        "M  V30 13 1 13 14\n",
        "M  V30 14 2 14 15\n",
        "M  V30 15 1 14 16\n",
        "M  V30 16 2 12 17\n",
        "M  V30 17 1 17 18\n",
        "M  V30 18 1 2 19\n",
        "M  V30 19 2 19 20\n",
        "M  V30 20 1 20 21\n",
        "M  V30 21 2 21 22\n",
        "M  V30 22 1 22 23\n",
        "M  V30 23 2 23 24\n",
        "M  V30 24 1 24 25\n",
        "M  V30 25 1 25 26\n",
        "M  V30 26 2 26 27\n",
        "M  V30 27 2 18 9\n",
        "M  V30 28 1 27 19\n",
        "M  V30 29 1 27 23\n",
        "M  V30 END BOND\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );

    let record = crate::io::molfile::read_mol_record_from_str(input).unwrap();
    let molecule = &record.molecule;

    let v2000 = mol_to_v2000_block(molecule).unwrap();
    let v3000 = mol_to_v3000_block(molecule).unwrap();

    assert!(
        v2000.contains("    3.0489    1.0287    0.0000 S   0  0  0  0  0  6  0  0  0  0  0  0\n")
    );
    assert!(v3000.contains("M  V30 2 S 3.048886 1.028691 0.000000 0 VAL=6\n"));
}

#[test]
fn molblock_writer_emits_cfg_for_rdkit_nontetrahedral_cu_center_in_3d_v3000_case() {
    let input = concat!(
        "\n",
        "     RDKit          3D\n",
        "\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 COUNTS 49 56 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 O 1.773199 -2.434971 -1.968751 0\n",
        "M  V30 2 C 2.960214 -2.358387 -1.523438 0\n",
        "M  V30 3 O 4.024209 -3.031780 -2.108819 0\n",
        "M  V30 4 Na 3.703515 -4.215160 -3.771975 0 VAL=1\n",
        "M  V30 5 C 3.196494 -1.521121 -0.344202 0\n",
        "M  V30 6 C 1.995702 -0.866080 0.176991 0\n",
        "M  V30 7 C 1.020444 -1.534533 0.916017 0\n",
        "M  V30 8 C 1.218857 -2.565929 1.848886 0\n",
        "M  V30 9 C 2.440536 -3.313019 2.031054 0\n",
        "M  V30 10 O 3.396304 -2.876893 2.920297 0\n",
        "M  V30 11 Na 5.116837 -3.996632 3.127426 0 VAL=1\n",
        "M  V30 12 O 2.658233 -4.375695 1.383407 0\n",
        "M  V30 13 C -0.000843 -2.712686 2.552710 0\n",
        "M  V30 14 C -0.280146 -3.681137 3.638974 0\n",
        "M  V30 15 C -0.887759 -1.764946 2.018356 0\n",
        "M  V30 16 C -2.092602 -1.343302 2.531843 0\n",
        "M  V30 17 C -2.932462 -0.429252 1.995405 0\n",
        "M  V30 18 N -2.756374 0.258420 0.898218 0\n",
        "M  V30 19 C -3.790520 0.874731 0.426566 0\n",
        "M  V30 20 C -4.865061 0.591442 1.367697 0\n",
        "M  V30 21 C -6.222541 1.104323 1.252495 0\n",
        "M  V30 22 O -7.074951 0.808646 2.127270 0\n",
        "M  V30 23 C -4.354470 -0.201504 2.326318 0\n",
        "M  V30 24 C -5.062290 -0.752115 3.501864 0\n",
        "M  V30 25 C -4.398848 -0.277176 4.787976 0\n",
        "M  V30 26 C -3.830809 1.399467 -0.806225 0\n",
        "M  V30 27 C -3.025408 0.993357 -1.894058 0\n",
        "M  V30 28 C -3.210553 1.122385 -3.276505 0\n",
        "M  V30 29 C -4.452100 1.449439 -3.922945 0\n",
        "M  V30 30 C -4.477517 1.531397 -5.223528 0\n",
        "M  V30 31 C -1.958483 0.868145 -3.896666 0\n",
        "M  V30 32 C -1.663773 0.890404 -5.348214 0\n",
        "M  V30 33 C -1.070687 0.596398 -2.867414 0\n",
        "M  V30 34 N -1.759441 0.672006 -1.777402 0 CHG=-1\n",
        "M  V30 35 N -0.217195 -1.170504 1.099634 0 CHG=-1\n",
        "M  V30 36 Cu -1.002261 0.310774 0.008258 0 CHG=2 VAL=4\n",
        "M  V30 37 N 0.690956 0.827364 -0.556579 0\n",
        "M  V30 38 C 0.977701 1.315867 -1.799854 0\n",
        "M  V30 39 C 0.326220 0.553642 -2.774122 0\n",
        "M  V30 40 C 1.879333 0.504370 -0.105601 0\n",
        "M  V30 41 C 2.629189 1.759753 -0.214384 0\n",
        "M  V30 42 C 1.918814 2.387942 -1.490749 0\n",
        "M  V30 43 C 2.948260 2.617556 -2.547185 0\n",
        "M  V30 44 C 4.079729 1.585089 -0.454607 0\n",
        "M  V30 45 C 4.714563 2.960733 -0.541839 0\n",
        "M  V30 46 C 4.439139 3.615664 0.780770 0\n",
        "M  V30 47 O 4.900111 4.911501 1.028356 0\n",
        "M  V30 48 Na 4.576947 5.881101 2.787487 0 VAL=1\n",
        "M  V30 49 O 3.801588 3.030910 1.680790 0\n",
        "M  V30 END ATOM\n",
        "M  V30 BEGIN BOND\n",
        "M  V30 1 2 1 2\n",
        "M  V30 2 1 2 3\n",
        "M  V30 3 1 3 4\n",
        "M  V30 4 1 2 5\n",
        "M  V30 5 1 5 6\n",
        "M  V30 6 2 6 7\n",
        "M  V30 7 1 7 8\n",
        "M  V30 8 1 8 9\n",
        "M  V30 9 1 9 10\n",
        "M  V30 10 1 10 11\n",
        "M  V30 11 2 9 12\n",
        "M  V30 12 2 8 13\n",
        "M  V30 13 1 13 14\n",
        "M  V30 14 1 13 15\n",
        "M  V30 15 2 15 16\n",
        "M  V30 16 1 16 17\n",
        "M  V30 17 2 17 18\n",
        "M  V30 18 1 18 19\n",
        "M  V30 19 1 19 20\n",
        "M  V30 20 1 20 21\n",
        "M  V30 21 2 21 22\n",
        "M  V30 22 2 20 23\n",
        "M  V30 23 1 23 24\n",
        "M  V30 24 1 24 25\n",
        "M  V30 25 2 19 26\n",
        "M  V30 26 1 26 27\n",
        "M  V30 27 2 27 28\n",
        "M  V30 28 1 28 29\n",
        "M  V30 29 2 29 30\n",
        "M  V30 30 1 28 31\n",
        "M  V30 31 1 31 32\n",
        "M  V30 32 2 31 33\n",
        "M  V30 33 1 33 34\n",
        "M  V30 34 1 7 35\n",
        "M  V30 35 9 35 36\n",
        "M  V30 36 9 37 36\n",
        "M  V30 37 1 37 38\n",
        "M  V30 38 2 38 39\n",
        "M  V30 39 2 37 40\n",
        "M  V30 40 1 40 41\n",
        "M  V30 41 1 41 42\n",
        "M  V30 42 1 42 43\n",
        "M  V30 43 1 41 44\n",
        "M  V30 44 1 44 45\n",
        "M  V30 45 1 45 46\n",
        "M  V30 46 1 46 47\n",
        "M  V30 47 1 47 48\n",
        "M  V30 48 2 46 49\n",
        "M  V30 49 1 40 6\n",
        "M  V30 50 1 35 15\n",
        "M  V30 51 1 23 17\n",
        "M  V30 52 9 18 36\n",
        "M  V30 53 1 34 27\n",
        "M  V30 54 1 39 33\n",
        "M  V30 55 9 34 36\n",
        "M  V30 56 1 42 38\n",
        "M  V30 END BOND\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );

    let record = crate::io::molfile::read_mol_record_from_str(input).unwrap();
    let block = mol_to_v3000_block(&record.molecule).unwrap();

    assert!(block.contains("M  V30 36 Cu -1.002261 0.310774 0.008258 0 CFG=2 CHG=2 VAL=4\n"));
}

#[test]
fn molblock_writer_accepts_rdkit_collapsed_hydrogen_roundtrip_from_v3000_source() {
    let input = concat!(
        "\n",
        "     RDKit          2D\n",
        "\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 COUNTS 4 3 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 C -1.299038 -0.750000 0.000000 7 MASS=13\n",
        "M  V30 2 C 0.000000 0.000000 0.000000 0\n",
        "M  V30 3 F 1.299038 -0.750000 0.000000 0\n",
        "M  V30 4 Cl 0.000000 1.500000 0.000000 0\n",
        "M  V30 END ATOM\n",
        "M  V30 BEGIN BOND\n",
        "M  V30 1 1 2 1 CFG=1\n",
        "M  V30 2 1 2 3\n",
        "M  V30 3 1 2 4\n",
        "M  V30 END BOND\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );

    let record = crate::io::molfile::read_mol_record_from_str(input).unwrap();
    let molecule = &record.molecule;

    let v2000 = mol_to_v2000_block(molecule).unwrap();
    let v3000 = mol_to_v3000_block(molecule).unwrap();

    assert!(
        v2000.contains("   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  7  0  0\n")
    );
    assert!(v2000.contains("  2  1  1  1\n"));
    assert!(v2000.contains("M  ISO  1   1  13\n"));
    assert!(v3000.contains("M  V30 1 C -1.299038 -0.750000 0.000000 7 MASS=13\n"));
    assert!(v3000.contains("M  V30 1 1 2 1 CFG=1\n"));
}

#[test]
fn mol_to_2d_sdf_record_appends_sdf_data_fields_and_delimiter() {
    let molecule = charged_isotope_molecule().with_sdf_data_field("ID", "cmpd-1");

    let record = mol_to_2d_sdf_record(&molecule, SdfFormat::V2000).unwrap();

    assert!(record.contains("M  END\n>  <ID>  \ncmpd-1\n\n$$$$\n"));
}

#[test]
fn molblock_write_params_route_v2000_and_sdf_record_paths() {
    let molecule = charged_isotope_molecule().with_sdf_data_field("ID", "cmpd-1");
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        force_2d: true,
        ..Default::default()
    };

    let block = mol_to_mol_block_with_params(&molecule, &params).unwrap();
    let record = mol_to_sdf_record_with_params(&molecule, &params).unwrap();

    assert!(block.contains("  COSMolKit          2D\n"));
    assert!(record.ends_with(">  <ID>  \ncmpd-1\n\n$$$$\n"));
}

#[test]
fn molblock_write_params_auto_upgrade_dative_molecule_to_v3000_like_rdkit() {
    let mut builder = Molecule::builder().with_name("dative");
    let n1 = builder.add_atom(AtomSpec::new(Element::N));
    let cu = builder.add_atom(AtomSpec::new(Element::from_atomic_number(29).unwrap()));
    let n2 = builder.add_atom(AtomSpec::new(Element::N));
    builder
        .add_bond(BondSpec::new(n1, cu, BondOrder::Dative))
        .unwrap();
    builder
        .add_bond(BondSpec::new(n2, cu, BondOrder::Dative))
        .unwrap();
    builder
        .set_2d_coordinates(vec![[-1.5, 0.0], [0.0, 0.0], [1.5, 0.0]])
        .unwrap();
    let molecule = builder.build().unwrap();

    let block = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            force_2d: true,
            kekulize: false,
            include_stereo: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(block.contains("999 V3000\n"));
    assert!(block.contains("M  V30 1 9 1 2\n"));
    assert!(block.contains("M  V30 2 9 3 2\n"));
}

#[test]
fn molblock_writer_emits_zero_bond_code_for_rdkit_quadruple_bond_case() {
    let mut builder = Molecule::builder().with_name("quadruple");
    let rh1 = builder.add_atom(AtomSpec::new(Element::from_atomic_number(45).unwrap()));
    let rh2 = builder.add_atom(AtomSpec::new(Element::from_atomic_number(45).unwrap()));
    builder.atom_mut(rh1).unwrap().set_formal_charge(-1);
    builder.atom_mut(rh2).unwrap().set_formal_charge(-1);
    builder
        .add_bond(BondSpec::new(rh1, rh2, BondOrder::Quadruple))
        .unwrap();
    builder
        .set_2d_coordinates(vec![[-0.75, 0.0], [0.75, 0.0]])
        .unwrap();
    let molecule = builder.build().unwrap();

    let v2000 = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            force_2d: true,
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let v3000 = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V3000,
            force_2d: true,
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(v2000.contains("  1  2  0  0\n"));
    assert!(v3000.contains("M  V30 1 0 1 2\n"));
}

#[test]
fn molblock_writer_uses_cached_total_valence_for_rdkit_sih_case() {
    let input = concat!(
        "\n",
        "     RDKit          3D\n",
        "\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 COUNTS 4 3 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 C -0.879252 -1.519843 -0.145444 0\n",
        "M  V30 2 Si 0.008612 0.002178 0.440476 0 VAL=4\n",
        "M  V30 3 C 1.774885 -0.007414 -0.148732 0\n",
        "M  V30 4 C -0.904244 1.525079 -0.146301 0\n",
        "M  V30 END ATOM\n",
        "M  V30 BEGIN BOND\n",
        "M  V30 1 1 1 2\n",
        "M  V30 2 1 2 3\n",
        "M  V30 3 1 2 4\n",
        "M  V30 END BOND\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );

    let record = crate::io::molfile::read_mol_record_from_str(input).unwrap();
    let molecule = &record.molecule;

    let v2000 = mol_to_mol_block_with_params(
        molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let v3000 = mol_to_mol_block_with_params(
        molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V3000,
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(v2000.contains("    0.0086    0.0022    0.4405 Si  0  0  0  0  0  4"));
    assert!(v3000.contains("M  V30 2 Si 0.008612 0.002178 0.440476 0 VAL=4\n"));
}

#[test]
fn molblock_writer_marks_rdkit_unspecified_double_bond_as_crossed() {
    let input = concat!(
        "\n",
        "     RDKit          2D\n",
        "\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 COUNTS 4 3 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 F 1.979613 -0.136500 0.000000 0\n",
        "M  V30 2 C 0.599379 0.450827 0.000000 0\n",
        "M  V30 3 C -0.599379 -0.450827 0.000000 0\n",
        "M  V30 4 Cl -1.979613 0.136500 0.000000 0\n",
        "M  V30 END ATOM\n",
        "M  V30 BEGIN BOND\n",
        "M  V30 1 1 1 2\n",
        "M  V30 2 2 2 3 CFG=2\n",
        "M  V30 3 1 3 4\n",
        "M  V30 END BOND\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );

    let record = crate::io::molfile::read_mol_record_from_str(input).unwrap();
    let molecule = &record.molecule;

    let v2000 = mol_to_mol_block_with_params(
        molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let v3000 = mol_to_mol_block_with_params(
        molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V3000,
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(v2000.contains("  2  3  2  3\n"));
    assert!(v3000.contains("M  V30 2 2 2 3 CFG=2\n"));
}

#[test]
fn molblock_writer_marks_unspecified_exocyclic_imine_as_crossed_like_rdkit() {
    let molecule = Molecule::from_smiles("CN=c1ccn(C)cc1")
        .unwrap()
        .with_2d_coordinates()
        .unwrap();
    let v2000 = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            include_stereo: true,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let v3000 = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V3000,
            include_stereo: true,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(v2000.contains("  2  3  2  3\n"));
    assert!(v3000.contains("M  V30 2 2 2 3 CFG=2\n"));
}

#[test]
fn molblock_writer_crosses_only_rdkit_eligible_unspecified_ring_double_bonds() {
    for (smiles, should_cross) in [("C1C=CCCCC1", false), ("C1C=CCCCCC1", true)] {
        let molecule = Molecule::from_smiles(smiles)
            .unwrap()
            .with_2d_coordinates()
            .unwrap();
        let v3000 = mol_to_mol_block_with_params(
            &molecule,
            &MolBlockWriteParams {
                format: SdfFormat::V3000,
                include_stereo: true,
                kekulize: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(v3000.contains(" CFG=2\n"), should_cross, "{smiles}");
    }
}

#[test]
fn molblock_writer_does_not_cross_implicit_hydrogen_imine_like_rdkit() {
    let molecule = Molecule::from_smiles("NC(O)=N").unwrap();

    let v2000 = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            include_stereo: true,
            ..Default::default()
        },
    )
    .unwrap();
    let v3000 = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V3000,
            include_stereo: true,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(v2000.contains("  2  4  2  0"));
    assert!(!v3000.contains(" CFG=2\n"));
}

#[test]
fn molblock_write_params_route_v3000_precision() {
    let molecule = charged_isotope_molecule();
    let params = MolBlockWriteParams {
        format: SdfFormat::V3000,
        force_2d: true,
        precision: 2,
        ..Default::default()
    };

    let record = mol_to_sdf_record_with_params(&molecule, &params).unwrap();

    assert!(record.contains("M  V30 1 C 1.25 -2.50 0.00 0 CHG=-1 MASS=13\n"));
    assert!(record.ends_with("$$$$\n"));
}

#[test]
fn molblock_write_params_kekulize_prepares_temporary_molecule() {
    let molecule = aromatic_benzene_molecule();
    let kekulized = mol_to_v2000_block(&molecule).unwrap();
    let aromatic = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            force_2d: true,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();

    // RDKit prepareMol() kekulizes the temporary molecule when requested;
    // the original aromatic molecule remains unchanged.
    assert!(kekulized.lines().any(|line| line.starts_with("  1  2  1")));
    assert!(aromatic.lines().any(|line| line.starts_with("  1  2  4")));
    assert!(molecule.bonds().iter().all(|bond| bond.is_aromatic()));
}

#[test]
fn molblock_write_params_include_stereo_false_still_writes_existing_non_parity_output() {
    let molecule = chiral_state_molecule();
    // RDKit's includeStereo flag does not suppress existing atom-line output here.
    // With 2D coordinates and no wedge bond information, both code paths succeed
    // and serialize the same atom record.
    let stereo = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            force_2d: true,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let no_stereo = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            force_2d: true,
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(stereo.contains("    0.0000    0.0000    0.0000 C   0"));
    assert!(no_stereo.contains("    0.0000    0.0000    0.0000 C   0"));
    assert!(stereo.ends_with("M  END\n"));
    assert!(no_stereo.ends_with("M  END\n"));
}

#[test]
fn molblock_writer_generates_missing_2d_coords_via_registered_operation() {
    let mut builder = Molecule::builder().with_name("needs-2d");
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
        .unwrap();
    let mol = builder.build().unwrap();

    let block = mol_to_v2000_block(&mol).unwrap();

    assert!(block.contains("V2000"));
    assert!(block.contains("  2  1 "));
}

#[test]
fn mol_to_v2000_block_writes_rgroup_pxa_and_zero_bond_extensions() {
    let molecule = zbo_extension_molecule();

    let block = mol_to_v2000_block(&molecule).unwrap();

    assert!(block.starts_with("zbo\n  COSMolKit          2D\n\n"));
    assert!(block.contains("M  RGP  1   2   7\n"));
    assert!(block.contains("M  ZBO  1   1   0\n"));
    assert!(block.contains("M  HYD  2   1   2   2   0\n"));
    assert!(block.contains("M  ZCH  1   1  -1\n"));
    assert!(block.contains("M  PXA   2 payload\n"));
}

#[test]
fn mol_to_v2000_block_writes_atom_list_query_lines() {
    let molecule = atom_list_query_molecule();

    let block = mol_to_v2000_block(&molecule).unwrap();

    assert!(block.contains("    0.0000    0.0000    0.0000 L   0"));
    assert!(block.contains("V    1 [#6,#7]\n"));
    assert!(block.contains("M  ALS   1  2 F C   N   \n"));
    assert!(block.contains("M  ALS   2  2 T O   S   \n"));
}

#[test]
fn mol_to_v2000_block_writes_supported_query_bond_type_codes() {
    let molecule = query_bond_molecule();

    let block = mol_to_v2000_block(&molecule).unwrap();

    assert!(block.contains("  1  2  8  0\n"));
    assert!(block.contains("  3  4  5  0\n"));
    assert!(block.contains("  5  6 42  0\n"));
    assert!(block.contains("  7  8  1  0  0  1\n"));
    assert!(block.contains("  9 10  6  0  0  2\n"));
}

#[test]
fn mol_to_v2000_block_writes_supported_bond_stereo_codes() {
    let molecule = v3000_bond_cfg_molecule();

    let block = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            force_2d: true,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();
    let without_stereo = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            force_2d: true,
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(block.contains("  1  2  1  1\n"));
    assert!(block.contains("  2  3  2  3\n"));
    assert!(without_stereo.contains("  1  2  1  1\n"));
    assert!(without_stereo.contains("  2  3  2  3\n"));
}

#[test]
fn mol_to_v2000_block_writes_molfile_chiral_flag_in_counts_line() {
    let molecule = chiral_flag_molecule();

    let block = mol_to_v2000_block(&molecule).unwrap();

    assert!(block.contains("  1  0  0  0  1  0  0  0  0  0999 V2000\n"));
}

#[test]
fn mol_to_v2000_block_writes_sgroup_lines() {
    let molecule = sgroup_molecule();

    let block = mol_to_v2000_block(&molecule).unwrap();

    assert!(block.contains("  2  1  0  2  0  0  0  0  0  0999 V2000\n"));
    assert!(block.contains("M  STY  2   1 SUP   2 DAT\n"));
    assert!(block.contains("M  SLB  1   1   7\n"));
    assert!(block.contains("M  SST  1   1 ALT\n"));
    assert!(block.contains("M  SCN  1   1 HT \n"));
    assert!(block.contains("M  SDS EXP  1   1\n"));
    assert!(block.contains("M  SPL  1   2   1\n"));
    assert!(block.contains("M  SNC  1   2   5\n"));
    assert!(block.contains("M  SBT  2   1   0   2   1\n"));
    assert!(block.contains("M  SAL   1  2   1   2\n"));
    assert!(block.contains("M  SPA   1  1   1\n"));
    assert!(block.contains("M  SBL   1  1   1\n"));
    assert!(block.contains("M  SDI   1  4    0.0000    1.0000    2.0000    3.0000\n"));
    assert!(block.contains("M  SMT   1 Me\n"));
    assert!(block.contains("M  SBV   1   1    0.5000    0.2500\n"));
    assert!(block.contains(
        "M  SDT   2 FIELD                         T INFO                Q OP             \n"
    ));
    assert!(block.contains("M  SDD   2 display spec\n"));
    assert!(block.contains("M  SED   2 first valuesecond value\n"));
    assert!(block.contains("M  SAP   1  1   1   2 AP\n"));
    assert!(block.contains("M  SCL   2 CLASS\n"));
}

#[test]
fn mol_to_v3000_block_writes_basic_ctab_atom_bond_blocks() {
    let mut builder = Molecule::builder()
        .with_name("v3000")
        .with_property("_MolFileChiralFlag", "1");
    let carbon = builder.add_atom(
        AtomSpec::new(Element::C)
            .with_formal_charge(-1)
            .with_isotope(13),
    );
    let oxygen = builder.add_atom(AtomSpec::new(Element::O));
    builder
        .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Double))
        .unwrap();
    builder
        .set_2d_coordinates(vec![[1.25, -2.5], [3.0, 4.0]])
        .unwrap();
    let molecule = builder.build().unwrap();

    let block = mol_to_2d_sdf_record(&molecule, SdfFormat::V3000).unwrap();

    assert!(block.starts_with("v3000\n  COSMolKit          2D\n\n"));
    assert!(block.contains("  0  0  0  0  0  0  0  0  0  0999 V3000\n"));
    assert!(block.contains("M  V30 BEGIN CTAB\n"));
    assert!(block.contains("M  V30 COUNTS 2 1 0 0 1\n"));
    assert!(block.contains("M  V30 BEGIN ATOM\n"));
    assert!(block.contains("M  V30 1 C 1.250000 -2.500000 0.000000 0 CHG=-1 MASS=13\n"));
    assert!(block.contains("M  V30 2 O 3.000000 4.000000 0.000000 0\n"));
    assert!(block.contains("M  V30 BEGIN BOND\n"));
    assert!(block.contains("M  V30 1 2 1 2\n"));
    assert!(block.contains("M  V30 END CTAB\nM  END\n$$$$\n"));
}

#[test]
fn mol_to_v3000_block_writes_supported_bond_cfg_values() {
    let molecule = v3000_bond_cfg_molecule();

    let block = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V3000,
            force_2d: true,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(block.contains("M  V30 1 1 1 2 CFG=1\n"));
    assert!(block.contains("M  V30 2 2 2 3 CFG=2\n"));
}

#[test]
fn mol_to_v3000_block_roundtrips_rdkit_atom_cfg_and_bond_cfg_from_3d_source() {
    let input = concat!(
        "\n",
        "     RDKit          3D\n",
        "\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 COUNTS 4 3 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 C -0.892286 -1.105272 -0.093154 7 MASS=13\n",
        "M  V30 2 C -0.086565 0.078288 0.363403 0 CFG=2\n",
        "M  V30 3 F -0.620779 1.215600 -0.194045 0\n",
        "M  V30 4 Cl 1.599630 -0.188615 -0.076204 0\n",
        "M  V30 END ATOM\n",
        "M  V30 BEGIN BOND\n",
        "M  V30 1 1 2 1 CFG=3\n",
        "M  V30 2 1 2 3\n",
        "M  V30 3 1 2 4\n",
        "M  V30 END BOND\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );

    let record = crate::io::molfile::read_mol_record_from_str(input).unwrap();
    let block = mol_to_v3000_block(&record.molecule).unwrap();

    assert!(block.contains("M  V30 2 C -0.086565 0.078288 0.363403 0 CFG=2\n"));
    assert!(block.contains("M  V30 1 1 2 1 CFG=3\n"));
    assert!(block.contains("M  V30 1 C -0.892286 -1.105272 -0.093154 7 MASS=13\n"));
}

#[test]
fn molblock_writer_ignores_raw_double_bond_neighbor_dirs_in_molfile_output() {
    let input = concat!(
        "\n",
        "     RDKit          2D\n",
        "\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 COUNTS 4 3 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 F -1.979613 -0.136500 0.000000 0\n",
        "M  V30 2 C -0.599379 0.450827 0.000000 0\n",
        "M  V30 3 C 0.599379 -0.450827 0.000000 0\n",
        "M  V30 4 F 1.979613 0.136500 0.000000 0\n",
        "M  V30 END ATOM\n",
        "M  V30 BEGIN BOND\n",
        "M  V30 1 1 1 2\n",
        "M  V30 2 2 2 3\n",
        "M  V30 3 1 3 4\n",
        "M  V30 END BOND\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );

    let record = crate::io::molfile::read_mol_record_from_str(input).unwrap();
    let v2000 = mol_to_v2000_block(&record.molecule).unwrap();
    let v3000 = mol_to_v3000_block(&record.molecule).unwrap();

    assert!(v2000.contains("  1  2  1  0\n"));
    assert!(v2000.contains("  2  3  2  0\n"));
    assert!(v2000.contains("  3  4  1  0\n"));
    assert!(v3000.contains("M  V30 1 1 1 2\n"));
    assert!(v3000.contains("M  V30 2 2 2 3\n"));
    assert!(v3000.contains("M  V30 3 1 3 4\n"));
    assert!(!v3000.contains("CFG="));
}

#[test]
fn mol_to_v3000_block_writes_sgroup_lines() {
    let molecule = sgroup_molecule();

    let block = mol_to_v3000_block(&molecule).unwrap();

    assert!(block.contains("M  V30 COUNTS 2 1 2 0 0\n"));
    assert!(block.contains("M  V30 BEGIN SGROUP\n"));
    assert!(block.contains("M  V30 1 SUP 7"));
    assert!(block.contains("ATOMS=(2 1 2)"));
    assert!(block.contains("XBONDS=(1 1)"));
    assert!(block.contains("PATOMS=(1 1)"));
    assert!(block.contains("SUBTYPE=ALT"));
    assert!(block.contains("CONNECT=HT"));
    assert!(block.contains("LABEL=Me"));
    assert!(block.contains("BRKXYZ=(9 0.0000 1.0000 0 2.0000 3.0000 0 0 0 0)"));
    assert!(block.contains("CSTATE=(4 1 0.5000 0.2500 0)"));
    assert_eq!(block.matches("CSTATE=(").count(), 1);
    assert!(block.contains("SAP=(3 1 2 AP)"));
    assert!(block.contains("M  V30 2 DAT 0"));
    assert!(block.contains("PARENT=1"));
    assert!(block.contains("COMPNO=5"));
    assert!(block.contains("FIELDNAME=FIELD"));
    assert!(block.contains("FIELDINFO=INFO"));
    assert!(block.contains("FIELDDISP=\"display spec\""));
    assert!(block.contains("QUERYTYPE=Q"));
    assert!(block.contains("QUERYOP=OP"));
    assert!(block.contains("FIELDDATA=\"first valuesecond value\""));
    assert!(block.contains("CLASS=CLASS"));
    assert!(block.contains("BRKTYP=PAREN"));
    assert!(block.contains("M  V30 END SGROUP\n"));
}

#[test]
fn mol_to_v3000_block_writes_zero_bond_sgroups() {
    let molecule = zbo_extension_molecule();

    let block = mol_to_v3000_block(&molecule).unwrap();

    assert!(block.contains("M  V30 COUNTS 2 1 3 0 0\n"));
    assert!(block.contains("M  V30 1 1 1 2\n"));
    assert!(block.contains("M  V30 BEGIN SGROUP\n"));
    assert!(block.contains("M  V30 1 DAT 0 ATOMS=(2 1 2) XBONDS=(1 1) FIELDNAME=ZBO\n"));
    assert!(block.contains("M  V30 2 DAT 0 ATOMS=(2 1 2) FIELDNAME=HYD FIELDDATA=\"2;0\"\n"));
    assert!(block.contains("M  V30 3 DAT 0 ATOMS=(2 1 2) FIELDNAME=ZCH FIELDDATA=\"-1;0\"\n"));
    assert!(block.contains("M  V30 END SGROUP\n"));
}

#[test]
fn mol_to_v3000_block_writes_enhanced_stereo_collections() {
    let mut builder = Molecule::builder().with_name("collections");
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::O));
    builder
        .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
        .unwrap();
    builder
        .add_stereo_group(StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![a0],
            Vec::new(),
        ))
        .unwrap();
    builder
        .add_stereo_group(StereoGroup::new(StereoGroupKind::Or, vec![a1], Vec::new()).with_id(2))
        .unwrap();
    let molecule = builder.build().unwrap();

    let block = mol_to_v3000_block(&molecule).unwrap();

    assert!(block.contains("M  V30 BEGIN COLLECTION\n"));
    assert!(block.contains("M  V30 MDLV30/STEABS ATOMS=(1 1)\n"));
    assert!(block.contains("M  V30 MDLV30/STEREL2 ATOMS=(1 2)\n"));
    assert!(block.contains("M  V30 END COLLECTION\n"));
}

#[test]
fn mol_to_v2000_block_infers_wedge_from_chiral_tag_without_coordinates() {
    // Build a chiral center with 4 explicit non-H neighbors, NO coordinates.
    // pick_bonds_to_wedge should select the bond to the lowest-scored neighbor
    // (lowest atomic number wins when degree and chirality are equal).
    let mut builder = Molecule::builder().with_name("no-coords-wedge");
    let c_chiral =
        builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
    let c_me = builder.add_atom(AtomSpec::new(Element::C));
    let n = builder.add_atom(AtomSpec::new(Element::N));
    let o = builder.add_atom(AtomSpec::new(Element::O));
    let f = builder.add_atom(AtomSpec::new(Element::F));
    builder
        .add_bond(BondSpec::new(c_chiral, c_me, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, n, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, o, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, f, BondOrder::Single))
        .unwrap();
    // Deliberately no coordinates
    let molecule = builder.build().unwrap();
    assert!(molecule.coordinates_2d().is_none());

    let block = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            include_stereo: true,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(block.starts_with("no-coords-wedge\n  COSMolKit          2D\n\n"));
    // The C-methyl bond (atom 1 → 2, added first) has the lowest atomic
    // number (6 vs 7/8/9) so pick_bond_to_wedge scores it lowest.
    // TetrahedralCw → BeginWedge → stereo code 1.
    // V2000 bond format: "{:>3}{:>3}{:>3} {:>2}" → stereo at [10..12].
    let has_wedge = block
        .lines()
        .any(|l| l.starts_with("  1  2") && l.len() >= 12 && l[10..12].trim() == "1");
    assert!(
        has_wedge,
        "expected BeginWedge (code 1) on bond 1-2 in output:\n{block}"
    );
}

#[test]
fn mol_to_v2000_block_infers_dash_from_chiral_tag_ccw_without_coordinates() {
    // Same molecule as above but with TetrahedralCcw → BeginDash → code 6.
    let mut builder = Molecule::builder().with_name("no-coords-dash");
    let c_chiral =
        builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw));
    let c_me = builder.add_atom(AtomSpec::new(Element::C));
    let n = builder.add_atom(AtomSpec::new(Element::N));
    let o = builder.add_atom(AtomSpec::new(Element::O));
    let f = builder.add_atom(AtomSpec::new(Element::F));
    builder
        .add_bond(BondSpec::new(c_chiral, c_me, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, n, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, o, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c_chiral, f, BondOrder::Single))
        .unwrap();
    let molecule = builder.build().unwrap();
    assert!(molecule.coordinates_2d().is_none());

    let block = mol_to_mol_block_with_params(
        &molecule,
        &MolBlockWriteParams {
            format: SdfFormat::V2000,
            include_stereo: true,
            kekulize: false,
            ..Default::default()
        },
    )
    .unwrap();

    assert!(block.starts_with("no-coords-dash\n  COSMolKit          2D\n\n"));
    let has_dash = block
        .lines()
        .any(|l| l.starts_with("  1  2") && l.len() >= 12 && l[10..12].trim() == "6");
    assert!(
        has_dash,
        "expected BeginDash (code 6) on bond 1-2 in output:\n{block}"
    );
}
