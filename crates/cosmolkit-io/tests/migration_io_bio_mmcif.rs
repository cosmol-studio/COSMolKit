use cosmolkit_bio::{
    BioConnectionKind, BioCoordinateFormat, BioHelixClass, BioSoftwareClassification,
};
use cosmolkit_io::{BioMmcifReadStage, read_mmcif_bio_structure};

const FULL_FEATURE_MMCIF: &str =
    include_str!("../../../testdata/bio/fixtures/gemmi_full_feature_sample.cif");

fn full_feature_mmcif_with_completion_categories() -> String {
    // The repository-authored base fixture already exercises block identity,
    // cell/transforms, authors, entity sequence and DBRefs, source identities,
    // atom sites, NCS, SIFTS, connections, cis-peptides, MODRES, and assemblies.
    // These additional rows cover the remaining completed category readers in
    // the source dispatcher; the final alias is deliberately resolved only
    // after the hierarchy and all address-bearing categories have been read.
    let mut text = FULL_FEATURE_MMCIF.replace("CYS", "~A");
    text.push_str(concat!(
        "\nloop_\n",
        "_refine.pdbx_refine_id\n",
        "_refine.ls_d_res_high\n",
        "_refine.ls_d_res_low\n",
        "_refine.ls_percent_reflns_obs\n",
        "_refine.ls_number_reflns_obs\n",
        "_refine.ls_number_reflns_R_work\n",
        "_refine.ls_number_reflns_R_free\n",
        "_refine.ls_R_factor_obs\n",
        "_refine.ls_R_factor_R_work\n",
        "_refine.ls_R_factor_R_free\n",
        "R1 2.4 3.6 94.5 1200 1100 100 0.21 0.19 0.24\n",
        "loop_\n",
        "_pdbx_refine_tls.id\n",
        "_pdbx_refine_tls.pdbx_refine_id\n",
        "_pdbx_refine_tls.T[1][1]\n",
        "_pdbx_refine_tls.T[2][2]\n",
        "_pdbx_refine_tls.T[3][3]\n",
        "_pdbx_refine_tls.T[1][2]\n",
        "_pdbx_refine_tls.T[1][3]\n",
        "_pdbx_refine_tls.T[2][3]\n",
        "_pdbx_refine_tls.L[1][1]\n",
        "_pdbx_refine_tls.L[2][2]\n",
        "_pdbx_refine_tls.L[3][3]\n",
        "_pdbx_refine_tls.L[1][2]\n",
        "_pdbx_refine_tls.L[1][3]\n",
        "_pdbx_refine_tls.L[2][3]\n",
        "_pdbx_refine_tls.S[1][1]\n",
        "_pdbx_refine_tls.S[1][2]\n",
        "_pdbx_refine_tls.S[1][3]\n",
        "_pdbx_refine_tls.S[2][1]\n",
        "_pdbx_refine_tls.S[2][2]\n",
        "_pdbx_refine_tls.S[2][3]\n",
        "_pdbx_refine_tls.S[3][1]\n",
        "_pdbx_refine_tls.S[3][2]\n",
        "_pdbx_refine_tls.S[3][3]\n",
        "_pdbx_refine_tls.origin_x\n",
        "_pdbx_refine_tls.origin_y\n",
        "_pdbx_refine_tls.origin_z\n",
        "tls1 R1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24\n",
        "loop_\n",
        "_pdbx_refine_tls_group.refine_tls_id\n",
        "_pdbx_refine_tls_group.beg_auth_asym_id\n",
        "_pdbx_refine_tls_group.beg_auth_seq_id\n",
        "_pdbx_refine_tls_group.beg_PDB_ins_code\n",
        "_pdbx_refine_tls_group.end_auth_seq_id\n",
        "_pdbx_refine_tls_group.end_PDB_ins_code\n",
        "_pdbx_refine_tls_group.selection_details\n",
        "tls1 X 101 . 102 . 'TLS selection'\n",
        "loop_\n",
        "_exptl.method\n",
        "_exptl.crystals_number\n",
        "'x-ray diffraction' 1\n",
        "loop_\n",
        "_exptl_crystal.id\n",
        "_exptl_crystal.description\n",
        "X 'fixture crystal'\n",
        "loop_\n",
        "_diffrn.id\n",
        "_diffrn.crystal_id\n",
        "_diffrn.ambient_temp\n",
        "D1 X 100.5\n",
        "loop_\n",
        "_diffrn_detector.diffrn_id\n",
        "_diffrn_detector.pdbx_collection_date\n",
        "_diffrn_detector.detector\n",
        "_diffrn_detector.type\n",
        "_diffrn_detector.details\n",
        "D1 2026-01-02 detector CCD optics\n",
        "loop_\n",
        "_diffrn_radiation.diffrn_id\n",
        "_diffrn_radiation.pdbx_scattering_type\n",
        "_diffrn_radiation.pdbx_monochromatic_or_laue_m_l\n",
        "_diffrn_radiation.monochromator\n",
        "D1 x-ray M 'Si(111)'\n",
        "loop_\n",
        "_diffrn_source.diffrn_id\n",
        "_diffrn_source.source\n",
        "_diffrn_source.type\n",
        "_diffrn_source.pdbx_synchrotron_site\n",
        "_diffrn_source.pdbx_synchrotron_beamline\n",
        "_diffrn_source.pdbx_wavelength_list\n",
        "D1 'storage ring' synchrotron Diamond I24 '0.98, 1.00'\n",
        "loop_\n",
        "_reflns.pdbx_diffrn_id\n",
        "_reflns.number_obs\n",
        "_reflns.d_resolution_high\n",
        "_reflns.d_resolution_low\n",
        "_reflns.percent_possible_obs\n",
        "_reflns.pdbx_redundancy\n",
        "_reflns.pdbx_Rmerge_I_obs\n",
        "_reflns.pdbx_Rsym_value\n",
        "_reflns.pdbx_netI_over_sigmaI\n",
        "D1 1200 2.25 3.75 97.5 1.4 0.02 0.03 8.5\n",
        "loop_\n",
        "_software.name\n",
        "_software.classification\n",
        "_software.version\n",
        "_software.date\n",
        "_software.description\n",
        "_software.contact_author\n",
        "_software.contact_author_email\n",
        "refiner 'model building' 1.0 2026-01-02 'model refinement' Ada ada@example.test\n",
        "loop_\n",
        "_atom_site_anisotrop.id\n",
        "_atom_site_anisotrop.U[1][1]\n",
        "_atom_site_anisotrop.U[2][2]\n",
        "_atom_site_anisotrop.U[3][3]\n",
        "_atom_site_anisotrop.U[1][2]\n",
        "_atom_site_anisotrop.U[1][3]\n",
        "_atom_site_anisotrop.U[2][3]\n",
        "1 1 2 3 4 5 6\n",
        "loop_\n",
        "_struct_conf.conf_type_id\n",
        "_struct_conf.beg_auth_asym_id\n",
        "_struct_conf.beg_label_comp_id\n",
        "_struct_conf.beg_auth_seq_id\n",
        "_struct_conf.pdbx_beg_PDB_ins_code\n",
        "_struct_conf.end_auth_asym_id\n",
        "_struct_conf.end_label_comp_id\n",
        "_struct_conf.end_auth_seq_id\n",
        "_struct_conf.pdbx_end_PDB_ins_code\n",
        "_struct_conf.pdbx_PDB_helix_class\n",
        "_struct_conf.pdbx_PDB_helix_length\n",
        "HELX_P X ~A 101 ? X ~A 102 ? 1 2\n",
        "loop_\n",
        "_struct_sheet.id\n",
        "S1\n",
        "loop_\n",
        "_struct_sheet_range.sheet_id\n",
        "_struct_sheet_range.id\n",
        "_struct_sheet_range.beg_auth_asym_id\n",
        "_struct_sheet_range.beg_label_comp_id\n",
        "_struct_sheet_range.beg_auth_seq_id\n",
        "_struct_sheet_range.pdbx_beg_PDB_ins_code\n",
        "_struct_sheet_range.end_auth_asym_id\n",
        "_struct_sheet_range.end_label_comp_id\n",
        "_struct_sheet_range.end_auth_seq_id\n",
        "_struct_sheet_range.pdbx_end_PDB_ins_code\n",
        "S1 A X ~A 101 ? X ~A 102 ?\n",
        "loop_\n",
        "_struct_sheet_order.sheet_id\n",
        "_struct_sheet_order.range_id_2\n",
        "_struct_sheet_order.sense\n",
        "S1 A P\n",
        "loop_\n",
        "_pdbx_struct_sheet_hbond.sheet_id\n",
        "_pdbx_struct_sheet_hbond.range_id_2\n",
        "_pdbx_struct_sheet_hbond.range_1_auth_asym_id\n",
        "_pdbx_struct_sheet_hbond.range_1_label_comp_id\n",
        "_pdbx_struct_sheet_hbond.range_1_auth_seq_id\n",
        "_pdbx_struct_sheet_hbond.range_1_PDB_ins_code\n",
        "_pdbx_struct_sheet_hbond.range_1_label_atom_id\n",
        "_pdbx_struct_sheet_hbond.range_2_auth_asym_id\n",
        "_pdbx_struct_sheet_hbond.range_2_label_comp_id\n",
        "_pdbx_struct_sheet_hbond.range_2_auth_seq_id\n",
        "_pdbx_struct_sheet_hbond.range_2_PDB_ins_code\n",
        "_pdbx_struct_sheet_hbond.range_2_label_atom_id\n",
        "S1 A X ~A 101 ? SG X ~A 102 ? SG\n",
        "loop_\n",
        "_chem_comp.id\n",
        "_chem_comp.three_letter_code\n",
        "~A CYS\n",
    ));
    text
}

#[test]
fn c30_reader_materializes_pinned_full_feature_fixture_in_source_order() {
    let text = full_feature_mmcif_with_completion_categories();
    let structure = read_mmcif_bio_structure(&text, "c30-full.cif").unwrap();
    structure.validate().unwrap();

    // Gemmi mmcif.hpp::make_structure selects the first coordinate block;
    // populate_structure_from_block's final CCD pass changes every matching
    // residue/address/sequence value back to its full source code.
    assert_eq!(structure.input_format(), BioCoordinateFormat::Mmcif);
    assert_eq!(structure.source_state().name, "demo");
    assert_eq!(structure.source_state().resolution, 2.4);
    assert_eq!(structure.source_state().info["_entry.id"], "9XYZ");
    assert_eq!(structure.metadata().authors, ["DOE, J.", "SMITH, A."]);

    assert_eq!(structure.models().len(), 1);
    assert_eq!(structure.chains().len(), 1);
    assert_eq!(structure.residues().len(), 2);
    assert_eq!(structure.atoms().len(), 2);
    assert_eq!(
        structure.coordinates().positions(),
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
    );
    assert_eq!(structure.atoms()[0].name().as_bytes(), b"SG");
    assert_eq!(
        structure.atoms()[0].anisou(),
        &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    );
    assert_eq!(structure.residues()[0].name().as_str(), "CYS");
    assert_eq!(structure.residues()[1].name().as_str(), "CYS");
    assert_eq!(structure.entities().len(), 1);
    assert_eq!(
        structure.entities()[0].full_sequence(),
        ["CYS".to_owned(), "CYS".to_owned()]
    );
    assert_eq!(structure.entities()[0].dbrefs().len(), 1);
    assert_eq!(structure.entities()[0].sifts_unp_accessions(), ["P12345"]);
    assert_eq!(structure.residues()[0].sifts_unp().residue(), Some(b'C'));
    assert_eq!(structure.residues()[0].sifts_unp().number(), 15);

    let crystal = structure.crystal().unwrap();
    assert_eq!(crystal.cell().a, 10.0);
    assert_eq!(crystal.space_group_hm(), Some("P 1"));
    assert_eq!(crystal.fractional().translation(), &[1.0, 2.0, 3.0]);
    assert_eq!(
        structure.source_state().origx.translation(),
        &[4.0, 5.0, 6.0]
    );
    // Pinned Gemmi omits identity transforms from the NCS vector and stores
    // only the identity operation's id as source metadata.
    assert_eq!(structure.ncs_operators().len(), 1);
    assert_eq!(structure.ncs_operators()[0].id, "N");
    assert_eq!(
        structure.source_state().info.get("_struct_ncs_oper.id"),
        Some(&"I".to_owned())
    );

    assert_eq!(structure.metadata().experiments.len(), 1);
    let experiment = &structure.metadata().experiments[0];
    assert_eq!(experiment.method, "x-ray diffraction");
    assert_eq!(experiment.number_of_crystals, 1);
    assert_eq!(experiment.diffraction_ids, ["D1"]);
    assert_eq!(experiment.unique_reflections, 1200);
    assert_eq!(experiment.reflections.resolution_high, 2.25);
    assert_eq!(structure.metadata().crystals.len(), 1);
    let diffraction = &structure.metadata().crystals[0].diffractions[0];
    assert_eq!(diffraction.id, "D1");
    assert_eq!(diffraction.temperature, 100.5);
    assert_eq!(diffraction.detector, "detector");
    assert_eq!(diffraction.detector_make, "CCD");
    assert_eq!(diffraction.optics, "optics");
    assert_eq!(diffraction.mono_or_laue, b'M');
    assert_eq!(diffraction.beamline, "I24");

    assert_eq!(structure.metadata().refinement.len(), 1);
    let refinement = &structure.metadata().refinement[0];
    assert_eq!(refinement.id, "R1");
    assert_eq!(refinement.basic.reflection_count, 1200);
    assert_eq!(refinement.tls_groups.len(), 1);
    assert_eq!(refinement.tls_groups[0].t, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    assert_eq!(refinement.tls_groups[0].origin, [22.0, 23.0, 24.0]);
    assert_eq!(refinement.tls_groups[0].selections[0].chain.as_str(), "X");
    assert_eq!(
        refinement.tls_groups[0].selections[0].details,
        "TLS selection"
    );
    assert_eq!(
        structure.metadata().software[0].classification,
        BioSoftwareClassification::ModelBuilding
    );

    assert_eq!(structure.connections().len(), 1);
    let connection = &structure.connections()[0];
    assert_eq!(connection.kind, BioConnectionKind::Disulf);
    assert_eq!(connection.partner1.residue().name().as_str(), "CYS");
    assert_eq!(connection.partner2.residue().name().as_str(), "CYS");
    assert_eq!(structure.cispeps().len(), 1);
    assert_eq!(
        structure.cispeps()[0].partner_c.residue().name().as_str(),
        "CYS"
    );
    assert_eq!(structure.mod_residues().len(), 1);
    assert_eq!(structure.mod_residues()[0].res_id.name().as_str(), "MSE");
    assert_eq!(structure.helices().len(), 1);
    assert_eq!(
        structure.helices()[0].pdb_helix_class,
        BioHelixClass::RAlpha
    );
    assert_eq!(
        structure.helices()[0].start.residue().name().as_str(),
        "CYS"
    );
    assert_eq!(structure.sheets().len(), 1);
    assert_eq!(structure.sheets()[0].strands.len(), 1);
    assert_eq!(structure.sheets()[0].strands[0].sense, 1);
    assert_eq!(
        structure.sheets()[0].strands[0]
            .hbond_atom1
            .logical_atom_name(),
        "SG"
    );
    assert_eq!(
        structure.sheets()[0].strands[0]
            .start
            .residue()
            .name()
            .as_str(),
        "CYS"
    );
    assert_eq!(structure.assemblies().len(), 1);
    assert_eq!(structure.assemblies()[0].generators.len(), 1);
    assert_eq!(structure.assemblies()[0].generators[0].operators.len(), 2);
}

#[test]
fn c30_reader_uses_auth_atom_name_then_source_label_fallback_without_padding() {
    let text = concat!(
        "data_atom_names\n",
        "loop_\n",
        "_atom_site.id\n",
        "_atom_site.group_PDB\n",
        "_atom_site.type_symbol\n",
        "_atom_site.label_atom_id\n",
        "_atom_site.auth_atom_id\n",
        "_atom_site.label_alt_id\n",
        "_atom_site.label_comp_id\n",
        "_atom_site.label_asym_id\n",
        "_atom_site.label_entity_id\n",
        "_atom_site.label_seq_id\n",
        "_atom_site.pdbx_PDB_ins_code\n",
        "_atom_site.Cartn_x\n",
        "_atom_site.Cartn_y\n",
        "_atom_site.Cartn_z\n",
        "_atom_site.occupancy\n",
        "_atom_site.B_iso_or_equiv\n",
        "_atom_site.pdbx_formal_charge\n",
        "_atom_site.auth_seq_id\n",
        "_atom_site.auth_comp_id\n",
        "_atom_site.auth_asym_id\n",
        "_atom_site.pdbx_PDB_model_num\n",
        "_atom_site.calc_flag\n",
        "_atom_site.pdbx_tls_group_id\n",
        "_atom_site.ccp4_deuterium_fraction\n",
        "1 ATOM C LABEL_NAME CA . GLY L 1 1 . 0 0 0 1 20 . 10 GLY AUTH 1 . . .\n",
        "2 ATOM C ' C1 ' . B GLY L 1 1 . 1 0 0 1 20 . 10 GLY AUTH 1 . . .\n",
    );
    let structure = read_mmcif_bio_structure(text, "atom-names.cif").unwrap();
    structure.validate().unwrap();
    assert_eq!(structure.atoms().len(), 2);
    assert_eq!(structure.residues().len(), 1);
    assert_eq!(structure.atoms()[0].name().as_bytes(), b"CA");
    assert_eq!(structure.atoms()[0].name().as_str(), "CA");
    assert_eq!(structure.atoms()[1].name().as_bytes(), b" C1 ");
    assert_eq!(structure.atoms()[1].name().as_str(), " C1 ");
}

#[test]
fn c30_reader_does_not_enter_the_separate_chemcomp_coordinate_reader() {
    let text = concat!(
        "data_component_only\n",
        "_chem_comp.id LIG\n",
        "_chem_comp.three_letter_code LIG\n",
        "loop_\n",
        "_chem_comp_atom.atom_id\n",
        "_chem_comp_atom.type_symbol\n",
        "_chem_comp_atom.charge\n",
        "_chem_comp_atom.x\n",
        "_chem_comp_atom.y\n",
        "_chem_comp_atom.z\n",
        "C1 C 0 1.0 2.0 3.0\n",
    );
    let structure = read_mmcif_bio_structure(text, "chemcomp-only.cif").unwrap();
    structure.validate().unwrap();
    assert_eq!(structure.input_format(), BioCoordinateFormat::Mmcif);
    assert!(structure.models().is_empty());
    assert!(structure.atoms().is_empty());
    assert!(structure.coordinates().is_empty());
}

#[test]
fn c30_reader_rejects_later_coordinate_blocks_before_materialization() {
    let text = "data_first\n_entry.id first\ndata_second\n_atom_site.id 1\n";
    let error = read_mmcif_bio_structure(text, "later-coordinates.cif").unwrap_err();
    assert_eq!(error.stage(), BioMmcifReadStage::CoordinateBlock);
    assert!(error.to_string().contains("block #2"));
    assert!(error.to_string().contains("later-coordinates.cif"));
}
