use cosmolkit_core::{RingInfo, ValenceAssignment, ValenceModel};
use cosmolkit_model::{AtomId, BondId, StereoGroup, StereoGroupKind};
use cosmolkit_smiles::{
    CxCoordinateSelection, CxSmilesFields, CxSmilesWriteParams, SmilesParseParams, SmilesRecord,
    SmilesWriteOutput, SmilesWriteParams, finalize_smiles_stereo, parse_smiles,
    write_fragment_cx_smiles, write_fragment_smiles_output,
};

fn parsed_record(smiles: &str) -> SmilesRecord {
    parse_smiles(smiles, &SmilesParseParams::default())
        .unwrap_or_else(|error| panic!("failed to parse fixed fragment input {smiles:?}: {error}"))
}

fn finalized_record(smiles: &str) -> SmilesRecord {
    let parse_params = SmilesParseParams::default();
    let parsed = parse_smiles(smiles, &parse_params)
        .unwrap_or_else(|error| panic!("failed to parse fixed fragment input {smiles:?}: {error}"));
    let remove_params = cosmolkit_core::RemoveHsParams {
        update_explicit_count: true,
        sanitize: parse_params.sanitize,
        ..cosmolkit_core::RemoveHsParams::default()
    };
    let prepared = cosmolkit_core::remove_hydrogens_with_params(
        parsed.topology,
        parsed.coordinates,
        parsed.properties,
        &remove_params,
    )
    .unwrap_or_else(|error| panic!("failed source-order H preparation: {error}"));
    // Move the REAL RH-prepared ring state into the finalizer carrier.
    let mut ring_carrier = prepared.final_rings;
    finalize_smiles_stereo(
        SmilesRecord {
            topology: prepared.topology,
            coordinates: prepared.coordinates,
            properties: prepared.properties,
        },
        &parse_params,
        &mut None,
        &mut ring_carrier,
    )
    .unwrap_or_else(|error| panic!("failed pinned stereo finalization: {error}"))
}

fn sanitized_finalized_record(smiles: &str) -> SmilesRecord {
    let parse_params = SmilesParseParams::default();
    let parsed = parse_smiles(smiles, &parse_params)
        .unwrap_or_else(|error| panic!("failed to parse fixed fragment input {smiles:?}: {error}"));
    let sanitized = cosmolkit_core::sanitize_topology(
        &parsed.topology,
        &cosmolkit_core::SanitizeParams::default(),
    )
    .unwrap_or_else(|error| panic!("failed source sanitize stage: {error}"));
    let remove_params = cosmolkit_core::RemoveHsParams {
        update_explicit_count: true,
        sanitize: parse_params.sanitize,
        ..cosmolkit_core::RemoveHsParams::default()
    };
    let prepared = cosmolkit_core::remove_hydrogens_with_params(
        sanitized.topology,
        parsed.coordinates,
        parsed.properties,
        &remove_params,
    )
    .unwrap_or_else(|error| panic!("failed source-order H preparation: {error}"));
    let mut ring_carrier = prepared.final_rings;
    finalize_smiles_stereo(
        SmilesRecord {
            topology: prepared.topology,
            coordinates: prepared.coordinates,
            properties: prepared.properties,
        },
        &parse_params,
        &mut None,
        &mut ring_carrier,
    )
    .unwrap_or_else(|error| panic!("failed pinned stereo finalization: {error}"))
}

fn all_indices(count: usize) -> Vec<usize> {
    (0..count).collect()
}

fn write_fragment(
    record: &SmilesRecord,
    params: &SmilesWriteParams,
    atoms: &[usize],
    bonds: Option<&[usize]>,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
    source_rings: Option<&RingInfo>,
    existing_valence: Option<&ValenceAssignment>,
) -> SmilesWriteOutput {
    let atom_ids = atoms.iter().copied().map(AtomId::new).collect::<Vec<_>>();
    let bond_ids = bonds.map(|ids| ids.iter().copied().map(BondId::new).collect::<Vec<_>>());
    write_fragment_smiles_output(
        record,
        params,
        &atom_ids,
        bond_ids.as_deref(),
        atom_symbols,
        bond_symbols,
        source_rings,
        existing_valence,
    )
    .unwrap_or_else(|error| panic!("fixed fragment writer input failed: {error}"))
}

fn write_fragment_cx(record: &SmilesRecord, fields: CxSmilesFields, atoms: &[usize]) -> String {
    let atom_ids = atoms.iter().copied().map(AtomId::new).collect::<Vec<_>>();
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams::default(),
        fields,
        coordinate_selection: CxCoordinateSelection::Auto,
    };
    write_fragment_cx_smiles(record, &params, &atom_ids, None, None, None, None, None)
        .unwrap_or_else(|error| panic!("fixed fragment CX writer input failed: {error}"))
}

fn assert_output(
    output: &SmilesWriteOutput,
    text: &str,
    atom_order: &[usize],
    bond_order: &[usize],
) {
    assert_eq!(output.text, text);
    assert_eq!(
        output.atom_order,
        atom_order
            .iter()
            .copied()
            .map(AtomId::new)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        output.bond_order,
        bond_order
            .iter()
            .copied()
            .map(BondId::new)
            .collect::<Vec<_>>()
    );
}

fn add_group(
    record: &mut SmilesRecord,
    kind: StereoGroupKind,
    atom_indices: &[usize],
    bond_indices: &[usize],
) {
    record.topology.stereo_groups.push(
        StereoGroup::new(
            kind,
            atom_indices.iter().copied().map(AtomId::new).collect(),
            bond_indices.iter().copied().map(BondId::new).collect(),
        )
        .with_id(7),
    );
}

#[test]
fn fragment_entry_matches_pinned_canonical_and_isomeric_profiles_with_original_maps() {
    // Pinned RDKit 2026.03.1 MolFragmentToSmiles, all original rows selected.
    let record = finalized_record("F[C@H](Cl)OC(F)Cl");
    let before = record.clone();
    let atoms = all_indices(record.topology.atoms.len());
    let natural_atoms = [0, 1, 2, 3, 4, 5, 6];
    let natural_bonds = [0, 1, 2, 3, 4, 5];

    for (canonical, isomeric, text, atom_order, bond_order) in [
        (
            false,
            false,
            "FC(Cl)OC(F)Cl",
            &natural_atoms[..],
            &natural_bonds[..],
        ),
        (
            false,
            true,
            "F[C@H](Cl)OC(F)Cl",
            &natural_atoms[..],
            &natural_bonds[..],
        ),
        (
            true,
            false,
            "FC(Cl)OC(F)Cl",
            &natural_atoms[..],
            &natural_bonds[..],
        ),
        (
            true,
            true,
            "FC(Cl)O[C@H](F)Cl",
            &[5, 4, 6, 3, 1, 0, 2],
            &[4, 5, 3, 2, 0, 1],
        ),
    ] {
        let params = SmilesWriteParams {
            canonical,
            do_isomeric_smiles: isomeric,
            ..SmilesWriteParams::default()
        };
        let output = write_fragment(&record, &params, &atoms, None, None, None, None, None);
        assert_output(&output, text, atom_order, bond_order);
        assert_eq!(record, before, "writer changed the shared input record");
    }
}

#[test]
fn fragment_entry_matches_pinned_enhanced_atom_bond_and_mixed_groups() {
    // Each group is constructed against the same parsed graph and source row
    // IDs as the pinned RDKit StereoGroup API; the noncanonical profile checks
    // the independent no-normalization branch.
    let cases: [(
        &str,
        StereoGroupKind,
        &[usize],
        &[usize],
        &str,
        &[usize],
        &[usize],
        &str,
        &[usize],
        &[usize],
    ); 3] = [
        (
            "F[C@@H](Cl)C[C@H](F)Cl",
            StereoGroupKind::Or,
            &[1, 4],
            &[],
            "F[C@H](Cl)C[C@@H](F)Cl",
            &[5, 4, 6, 3, 1, 0, 2],
            &[4, 5, 3, 2, 0, 1],
            "FC(Cl)CC(F)Cl",
            &[0, 1, 2, 3, 4, 5, 6],
            &[0, 1, 2, 3, 4, 5],
        ),
        (
            "F/C=C\\Cl",
            StereoGroupKind::And,
            &[],
            &[1],
            "F/C=C\\Cl",
            &[0, 1, 2, 3],
            &[0, 1, 2],
            "FC=CCl",
            &[0, 1, 2, 3],
            &[0, 1, 2],
        ),
        (
            "F[C@@H](Cl)/C=C\\Cl",
            StereoGroupKind::And,
            &[1],
            &[3],
            "F[C@H](Cl)/C=C\\Cl",
            &[0, 1, 2, 3, 4, 5],
            &[0, 1, 2, 3, 4],
            "FC(Cl)C=CCl",
            &[0, 1, 2, 3, 4, 5],
            &[0, 1, 2, 3, 4],
        ),
    ];

    for (
        smiles,
        kind,
        group_atoms,
        group_bonds,
        canonical_text,
        canonical_atoms,
        canonical_bonds,
        noncanonical_text,
        noncanonical_atoms,
        noncanonical_bonds,
    ) in cases
    {
        let mut record = finalized_record(smiles);
        add_group(&mut record, kind, group_atoms, group_bonds);
        let before = record.clone();
        let atoms = all_indices(record.topology.atoms.len());

        let canonical = write_fragment(
            &record,
            &SmilesWriteParams::default(),
            &atoms,
            None,
            None,
            None,
            None,
            None,
        );
        assert_output(&canonical, canonical_text, canonical_atoms, canonical_bonds);

        let noncanonical_params = SmilesWriteParams {
            canonical: false,
            do_isomeric_smiles: false,
            ..SmilesWriteParams::default()
        };
        let noncanonical = write_fragment(
            &record,
            &noncanonical_params,
            &atoms,
            None,
            None,
            None,
            None,
            None,
        );
        assert_output(
            &noncanonical,
            noncanonical_text,
            noncanonical_atoms,
            noncanonical_bonds,
        );
        assert_eq!(record, before, "writer changed the grouped input {smiles}");
    }
}

#[test]
fn fragment_entry_matches_pinned_selected_root_and_disconnected_components() {
    let rooted = finalized_record("CC(O)N");
    let rooted_before = rooted.clone();
    let rooted_params = SmilesWriteParams {
        rooted_at_atom: Some(AtomId::new(2)),
        ..SmilesWriteParams::default()
    };
    let rooted_output = write_fragment(
        &rooted,
        &rooted_params,
        &[0, 1, 2, 3],
        None,
        None,
        None,
        None,
        None,
    );
    assert_output(&rooted_output, "OC(C)N", &[2, 1, 0, 3], &[1, 0, 2]);
    assert_eq!(rooted, rooted_before);

    let disconnected = finalized_record("CCO.CN");
    let disconnected_before = disconnected.clone();
    let disconnected_output = write_fragment(
        &disconnected,
        &SmilesWriteParams::default(),
        &[0, 2, 4],
        None,
        None,
        None,
        None,
        None,
    );
    assert_output(&disconnected_output, "C.N.O", &[0, 4, 2], &[]);
    assert_eq!(disconnected, disconnected_before);
}

#[test]
fn fragment_entry_uses_full_index_custom_symbols_and_s60_map_policy() {
    let symbols_record = finalized_record("CCO.CN");
    let symbols_before = symbols_record.clone();
    let atom_symbols = ["Q", "R", "S", "T", "U"]
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let bond_symbols = ["x", "y", "z"]
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let symbols_output = write_fragment(
        &symbols_record,
        &SmilesWriteParams::default(),
        &[0, 1, 2],
        None,
        Some(&atom_symbols),
        Some(&bond_symbols),
        None,
        None,
    );
    assert_output(&symbols_output, "QxRyS", &[0, 1, 2], &[0, 1]);
    assert_eq!(symbols_record, symbols_before);

    let mapped = finalized_record("[NH2:1]c1ccccc1");
    let mapped_before = mapped.clone();
    let atoms = all_indices(mapped.topology.atoms.len());
    for (ignore_atom_maps, text, atom_order, bond_order) in [
        (
            false,
            "c1ccc([NH2:1])cc1",
            &[4, 3, 2, 1, 0, 6, 5][..],
            &[3, 2, 1, 0, 6, 5, 4][..],
        ),
        (
            true,
            "[NH2:1]c1ccccc1",
            &[0, 1, 2, 3, 4, 5, 6][..],
            &[0, 1, 2, 3, 4, 5, 6][..],
        ),
    ] {
        let params = SmilesWriteParams {
            ignore_atom_map_numbers: ignore_atom_maps,
            ..SmilesWriteParams::default()
        };
        let output = write_fragment(&mapped, &params, &atoms, None, None, None, None, None);
        assert_output(&output, text, atom_order, bond_order);
        assert_eq!(mapped, mapped_before);
    }
}

#[test]
fn fragment_entry_kekulizes_full_and_partial_rings_at_source_masks() {
    let full_ring = finalized_record("c1ccccc1");
    let full_before = full_ring.clone();
    let full_params = SmilesWriteParams {
        do_kekule: true,
        ..SmilesWriteParams::default()
    };
    let full_output = write_fragment(
        &full_ring,
        &full_params,
        &[0, 1, 2, 3, 4, 5],
        None,
        None,
        None,
        None,
        None,
    );
    assert_output(
        &full_output,
        "C1=CC=CC=C1",
        &[2, 1, 0, 5, 4, 3],
        &[1, 0, 5, 4, 3, 2],
    );
    assert_eq!(full_ring, full_before);

    let partial_ring = finalized_record("c1ccccc1");
    let partial_before = partial_ring.clone();
    let selected_bonds = [0, 1, 2, 3, 4];
    let partial_output = write_fragment(
        &partial_ring,
        &full_params,
        &[0, 1, 2, 3, 4, 5],
        Some(&selected_bonds),
        None,
        None,
        None,
        None,
    );
    assert_output(
        &partial_output,
        "C:C:C:C:C:C",
        &[0, 1, 2, 3, 4, 5],
        &[0, 1, 2, 3, 4],
    );
    assert_eq!(partial_ring, partial_before);
}

#[test]
fn fragment_entry_transfers_s59_cache_ring_and_stereo_preparation_cases() {
    let smiles = "C1[C@H](F)CC[C@H](Cl)C1";
    let expected_isomeric = "F[C@H]1CC[C@@H](Cl)CC1";
    let expected_isomeric_atoms = [2, 1, 3, 4, 5, 6, 7, 0];
    let expected_isomeric_bonds = [1, 2, 3, 4, 5, 6, 7, 0];
    let atoms = [0, 1, 2, 3, 4, 5, 6, 7];

    let mut legacy = parsed_record(smiles);
    legacy.properties.clear_prop("_StereochemDone");
    let legacy_before = legacy.clone();
    let output = write_fragment(
        &legacy,
        &SmilesWriteParams::default(),
        &atoms,
        None,
        None,
        None,
        None,
        None,
    );
    assert_output(
        &output,
        expected_isomeric,
        &expected_isomeric_atoms,
        &expected_isomeric_bonds,
    );
    assert_eq!(legacy, legacy_before);

    let nonisomeric_params = SmilesWriteParams {
        do_isomeric_smiles: false,
        ..SmilesWriteParams::default()
    };
    let nonisomeric = write_fragment(
        &legacy,
        &nonisomeric_params,
        &atoms,
        None,
        None,
        None,
        None,
        None,
    );
    assert_output(
        &nonisomeric,
        "FC1CCC(Cl)CC1",
        &[2, 1, 0, 7, 5, 6, 4, 3],
        &[1, 0, 7, 6, 5, 4, 3, 2],
    );
    assert_eq!(legacy, legacy_before);

    let mut prepared = finalized_record(smiles);
    prepared
        .properties
        .set_prop("_StereochemDone", "0")
        .unwrap();
    prepared.topology.atoms[1]
        .set_prop("_CIPCode", "R")
        .unwrap();
    prepared.topology.atoms[5]
        .set_prop("_CIPCode", "S")
        .unwrap();
    let prepared_before = prepared.clone();
    let source_rings = cosmolkit_core::fast_find_rings(&prepared.topology).unwrap();
    let source_rings_before = source_rings.clone();
    let existing_valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &prepared.topology,
        ValenceModel::RdkitLike,
        false,
    )
    .unwrap();
    let existing_valence_before = existing_valence.clone();
    let output = write_fragment(
        &prepared,
        &SmilesWriteParams::default(),
        &atoms,
        None,
        None,
        None,
        Some(&source_rings),
        Some(&existing_valence),
    );
    assert_output(
        &output,
        expected_isomeric,
        &expected_isomeric_atoms,
        &expected_isomeric_bonds,
    );
    assert_eq!(prepared, prepared_before);
    assert_eq!(source_rings, source_rings_before);
    assert_eq!(existing_valence, existing_valence_before);
}

#[test]
fn fragment_cx_filters_fields_and_preserves_source_zero_sgroup_rows() {
    // Pinned RDKit 2026.03.1 MolFragmentToCXSmiles: get_sgroup_data_block
    // and get_sgroup_polymer_block use zero-initialized reverse atom maps.
    let data = parsed_record("CCO |SgD:0,2:FIELD:VALUE::::|");
    let data_before = data.clone();
    assert_eq!(write_fragment_cx(&data, CxSmilesFields::NONE, &[0]), "C");
    assert_eq!(
        write_fragment_cx(&data, CxSmilesFields::SGROUPS, &[0]),
        "C |SgD:0,0:FIELD:VALUE::::|"
    );
    assert_eq!(
        data, data_before,
        "CX fragment writing changed data SGroup input"
    );

    let polymer = parsed_record("CCO |Sg:n:0,2::eu:::|");
    let polymer_before = polymer.clone();
    assert_eq!(write_fragment_cx(&polymer, CxSmilesFields::NONE, &[1]), "C");
    assert_eq!(
        write_fragment_cx(&polymer, CxSmilesFields::POLYMER, &[1]),
        "C |Sg:n:0,0::eu:::|"
    );
    assert_eq!(
        polymer, polymer_before,
        "CX fragment writing changed fully excluded polymer SGroup input"
    );
}

#[test]
fn fragment_writer_gates_tagged_nonpotential_chirality_by_source_marker_presence() {
    use cosmolkit_types::ChiralTag;

    let mut template = parsed_record("FNC");
    template.topology.atoms[1].set_chiral_tag(ChiralTag::TetrahedralCw);

    for (clean_stereo, done_marker_present, expected) in [
        (false, false, "CNF"),
        (true, false, "CNF"),
        (false, true, "C[N@@H]F"),
        (true, true, "C[N@@H]F"),
    ] {
        let mut input = template.clone();
        if done_marker_present {
            input.properties.set_prop("_StereochemDone", "0").unwrap();
        } else {
            input.properties.clear_prop("_StereochemDone");
        }
        let before = input.clone();
        let params = SmilesWriteParams {
            canonical: false,
            clean_stereo,
            rooted_at_atom: Some(AtomId::new(2)),
            ..SmilesWriteParams::default()
        };
        let output = write_fragment(&input, &params, &[0, 1, 2], None, None, None, None, None);

        assert_eq!(
            output.text, expected,
            "clean_stereo={clean_stereo}, _StereochemDone present={done_marker_present}"
        );
        assert_eq!(
            input, before,
            "writer changed the nonpotential tagged input"
        );
    }
}

#[test]
fn fragment_writer_preserves_valid_tetrahedral_and_nontetrahedral_guard_paths() {
    let valid = parsed_record("F[C@H](Cl)Br");
    for (clean_stereo, done_marker_present) in
        [(false, false), (true, false), (false, true), (true, true)]
    {
        let mut input = valid.clone();
        if done_marker_present {
            input.properties.set_prop("_StereochemDone", "0").unwrap();
        } else {
            input.properties.clear_prop("_StereochemDone");
        }
        let before = input.clone();
        let params = SmilesWriteParams {
            canonical: false,
            clean_stereo,
            rooted_at_atom: Some(AtomId::new(0)),
            ..SmilesWriteParams::default()
        };
        let output = write_fragment(&input, &params, &[0, 1, 2, 3], None, None, None, None, None);

        assert_eq!(
            output.text, "F[C@H](Cl)Br",
            "clean_stereo={clean_stereo}, _StereochemDone present={done_marker_present}"
        );
        assert_eq!(input, before, "writer changed the valid tetrahedral input");
    }

    let nontetrahedral = parsed_record("[Pt@SP1](F)(Cl)(Br)I");
    let before = nontetrahedral.clone();
    let output = write_fragment(
        &nontetrahedral,
        &SmilesWriteParams::default(),
        &all_indices(nontetrahedral.topology.atoms.len()),
        None,
        None,
        None,
        None,
        None,
    );
    assert_eq!(output.text, "[F][Pt@SP1]([Cl])([Br])[I]");
    assert_eq!(
        nontetrahedral, before,
        "writer changed non-tetrahedral input"
    );
}

#[test]
fn fragment_writer_reuses_retained_ring_rows_for_nitrogen_and_skips_broken_centers() {
    let mut ring_nitrogen = sanitized_finalized_record("F[N@]1C(C)C1");
    ring_nitrogen
        .properties
        .set_prop("_StereochemDone", "0")
        .unwrap();
    let ring_before = ring_nitrogen.clone();
    let source_rings = cosmolkit_core::fast_find_rings(&ring_nitrogen.topology).unwrap();
    let source_rings_before = source_rings.clone();
    let existing_valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &ring_nitrogen.topology,
        ValenceModel::RdkitLike,
        false,
    )
    .unwrap();
    let existing_valence_before = existing_valence.clone();
    let ring_params = SmilesWriteParams {
        canonical: false,
        rooted_at_atom: Some(AtomId::new(3)),
        ..SmilesWriteParams::default()
    };
    let ring_output = write_fragment(
        &ring_nitrogen,
        &ring_params,
        &all_indices(ring_nitrogen.topology.atoms.len()),
        None,
        None,
        None,
        Some(&source_rings),
        Some(&existing_valence),
    );
    assert_eq!(ring_output.text, "CC1[N@@](F)C1");
    assert_eq!(ring_nitrogen, ring_before);
    assert_eq!(source_rings, source_rings_before);
    assert_eq!(existing_valence, existing_valence_before);

    let mut cut = parsed_record("F[C@H](Cl)Br");
    cut.properties.set_prop("_StereochemDone", "0").unwrap();
    let cut_before = cut.clone();
    let cut_params = SmilesWriteParams {
        canonical: false,
        rooted_at_atom: Some(AtomId::new(0)),
        ..SmilesWriteParams::default()
    };
    let cut_output = write_fragment(
        &cut,
        &cut_params,
        &[0, 1, 2],
        Some(&[0, 1]),
        None,
        None,
        None,
        None,
    );
    assert_output(&cut_output, "FCCl", &[0, 1, 2], &[0, 1]);
    assert_eq!(cut, cut_before, "writer changed selected-fragment input");
}
