// BIO-CID C26: BioSelection type + parse/match error vocabulary exposed
// and registered as Experimental supporting values/errors. Expectations
// derive independently from dev/public_api_design.md (registry describes
// real public types; Experimental support; declarations are not
// implementations) and the frozen CID contract: a thin private-field root
// wrapper, no parser helpers or mutable storage re-exported.

use cosmolkit::binding_contract::{BINDING_CONTRACT, BindingItem, BindingTypeRole, FunctionStatus};

#[test]
fn bio_cid_c26_registry_type_identities() {
    let mut found = 0;
    for entry in BINDING_CONTRACT {
        match entry.semantic_id {
            "types.BioSelection" => {
                assert_eq!(entry.item, BindingItem::Type);
                assert_eq!(entry.type_role, Some(BindingTypeRole::Value));
                assert_eq!(entry.feature, "cap-bio");
                assert_eq!(entry.status, FunctionStatus::Experimental);
                assert!(entry.rust_path.ends_with("BioSelection"));
                assert_eq!(entry.python_name, "BioSelection");
                assert_eq!(entry.javascript_name, "BioSelection");
                found += 1;
            }
            "types.BioSelectionParseError" => {
                assert_eq!(entry.item, BindingItem::Type);
                assert_eq!(entry.type_role, Some(BindingTypeRole::Error));
                assert_eq!(entry.feature, "cap-bio");
                assert_eq!(entry.status, FunctionStatus::Experimental);
                found += 1;
            }
            "types.BioSelectionMatchError" => {
                assert_eq!(entry.item, BindingItem::Type);
                assert_eq!(entry.type_role, Some(BindingTypeRole::Error));
                assert_eq!(entry.feature, "cap-bio");
                assert_eq!(entry.status, FunctionStatus::Experimental);
                found += 1;
            }
            _ => {}
        }
    }
    assert_eq!(found, 3, "exactly the three C26 type entries");
}

#[test]
fn bio_cid_c26_error_vocabulary_surface() {
    use cosmolkit::{BioSelectionMatchError, BioSelectionParseError, SelectionSyntaxError};
    use std::error::Error;

    // The parse-error vocabulary is exposed with its exact source accessor
    // surface (method-item presence) — no construction surface leaks
    // (fields are private; only the canonical reader constructs them).
    let cid_of: fn(&SelectionSyntaxError) -> &str = SelectionSyntaxError::cid;
    let pos_of: fn(&SelectionSyntaxError) -> usize = SelectionSyntaxError::pos;
    let info_of: fn(&SelectionSyntaxError) -> Option<&str> = SelectionSyntaxError::info;
    let _ = (cid_of, pos_of, info_of);
    fn assert_error<T: Error>() {}
    assert_error::<BioSelectionParseError>();
    assert_error::<BioSelectionMatchError>();
}

#[test]
fn bio_cid_c27_from_cid_valid_malformed_and_registry() {
    use cosmolkit::binding_contract::{BINDING_CONTRACT, BindingKind};
    use cosmolkit::{BioSelection, BioSelectionParseError, SelectionSyntaxError};

    // Valid constructors across every grammar stage; malformed inputs
    // carry the exact source byte observations (C14 evidence rules).
    for valid in [
        "/",
        "//",
        "A",
        "//A",
        "/2/A",
        "A/3-4",
        "A/(ALA)",
        "A//CA[C]",
        "A//:A",
        ";polymer",
        ";solvent",
        "A/;q>0.5",
        "/2/A/14-20/CA[C]:B",
    ] {
        BioSelection::from_cid(valid).unwrap_or_else(|e| panic!("{valid}: {e}"));
    }
    let error = match BioSelection::from_cid("A/(ALA") {
        Err(error) => error,
        Ok(_) => panic!("expected malformed CID to fail"),
    };
    match &error {
        BioSelectionParseError::Syntax(syntax) => {
            assert_eq!(syntax.cid(), "A/(ALA");
            // Exact wrong_syntax observation: the residue-stage terminal
            // gate reports offset 0 with the stage note (select.cpp:199-201
            // passes 0, not the field end).
            assert_eq!(syntax.pos(), 0);
        }
        other => panic!("expected Syntax, got {other:?}"),
    }
    let _ = error.to_string();

    // Registry: the registered constructor has the exact frozen
    // signature/projection/error vocabulary.
    let entry = BINDING_CONTRACT
        .iter()
        .find(|e| e.semantic_id == "BioSelection.from_cid")
        .expect("from_cid registered");
    let callable = entry.callable.expect("from_cid callable contract");
    assert_eq!(callable.kind, BindingKind::Static);
    assert_eq!(entry.feature, "cap-bio");
    assert_eq!(entry.python_name, "from_cid");
    assert_eq!(entry.javascript_name, "fromCid");
}

#[test]
fn bio_cid_c28_to_cid_roundtrips_and_registry() {
    use cosmolkit::BioSelection;
    use cosmolkit::binding_contract::{BINDING_CONTRACT, BindingKind, StateModel};

    // Serialization through the registered thin delegate reproduces the
    // C22 owner outputs (source Selection::str semantics), including the
    // documented non-roundtrips (blank-icode dots, wildcard residues).
    let rows: [(&str, &str); 8] = [
        ("/", "//*//"),
        ("//", "////"),
        ("A", "//A//"),
        ("//A", "//A//"),
        ("/2/A", "/2/A//"),
        ("A/3-4", "//A/3.-4./"),
        ("A/(ALA)", "//A/(ALA)/"),
        ("/2/A/14-20/CA[C]:B", "/2/A/14.-20./CA[C]:B"),
    ];
    for (cid, expected) in rows {
        assert_eq!(
            BioSelection::from_cid(cid).unwrap().to_cid(),
            expected,
            "{cid}"
        );
    }

    let entry = BINDING_CONTRACT
        .iter()
        .find(|e| e.semantic_id == "BioSelection.to_cid")
        .expect("to_cid registered");
    let callable = entry.callable.expect("to_cid callable contract");
    assert_eq!(callable.kind, BindingKind::Instance);
    assert_eq!(callable.state_model, StateModel::ReadOnly);
    assert_eq!(entry.python_name, "to_cid");
    assert_eq!(entry.javascript_name, "toCid");
}

#[test]
fn bio_cid_c29_public_pdb_mmcif_selection_matrices_and_sharing() {
    use cosmolkit::binding_contract::{BINDING_CONTRACT, BindingKind, StateModel};
    use cosmolkit::{BioSelection, BioStructure};

    // Public PDB matrix: real fixture -> real CID text -> original IDs
    // (expectations independently fixed from the fixture text and the
    // Gemmi matching closure; see the C25 matrix in selection_receipt.md).
    let pdb = BioStructure::read_with_format(
        std::path::Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../testdata/bio/fixtures/gemmi_full_feature_sample.pdb"
        )),
        cosmolkit::BioCoordinateFormat::Pdb,
    )
    .unwrap();
    let pdb_rows: [(&str, &[u32]); 10] = [
        ("A/3-4", &[0, 1]),
        ("A/7", &[2]),
        ("A/(ALA)", &[2]),
        ("A/*", &[0, 1, 2, 3, 4]),
        ("/1/A", &[0, 1, 2, 3, 4]),
        ("/2/A", &[]),
        ("A//SG[S]", &[0, 1]),
        ("A//O[O]", &[3]),
        ("A/;q>0.5", &[0, 1, 2, 3, 4]),
        ("A/;b<10", &[]),
    ];
    // Immutable sharing peer: a clone shares the Arc-backed row tables, so
    // pointer identity against the peer proves the queries preserve the
    // PERSISTENT storage (no detach/copy of the tables); this says
    // nothing about transient allocations inside a query call
    // (BIO-CID-SHARE-CLOSE wording).
    let pdb_peer = pdb.clone();
    let pointers = |s: &BioStructure| {
        (
            s.models().as_ptr(),
            s.models().len(),
            s.chains().as_ptr(),
            s.chains().len(),
            s.residues().as_ptr(),
            s.residues().len(),
            s.atoms().as_ptr(),
            s.atoms().len(),
        )
    };
    let positions = |s: &BioStructure| -> Vec<Option<[u64; 3]>> {
        (0..s.atoms().len())
            .map(|index| {
                s.atom_position(cosmolkit::BioAtomId::new(index as u32))
                    .map(|xyz| [xyz[0].to_bits(), xyz[1].to_bits(), xyz[2].to_bits()])
            })
            .collect()
    };
    let pdb_before = format!("{pdb:?}"); // presentation check only
    let pdb_pointers_before = pointers(&pdb);
    let pdb_positions_before = positions(&pdb);
    // Complete checkpoint closure (BIO-CID-SHARE-CLOSE): current four
    // row pointers AND lengths equal the ORIGINAL tuple, equal the PEER
    // tuple, four explicit peer ptr::eq identities, and the original
    // coordinate-bit snapshot. Applied immediately BEFORE and AFTER
    // every query call below; originals are never refreshed.
    let pdb_check = |label: &str| {
        let current = pointers(&pdb);
        let peer_now = pointers(&pdb_peer);
        assert_eq!(current, pdb_pointers_before, "original pointers {label}");
        assert_eq!(current, peer_now, "peer tuple {label}");
        assert!(
            std::ptr::eq(current.0, pdb_peer.models().as_ptr()),
            "models {label}"
        );
        assert!(
            std::ptr::eq(current.2, pdb_peer.chains().as_ptr()),
            "chains {label}"
        );
        assert!(
            std::ptr::eq(current.4, pdb_peer.residues().as_ptr()),
            "residues {label}"
        );
        assert!(
            std::ptr::eq(current.6, pdb_peer.atoms().as_ptr()),
            "atoms {label}"
        );
        assert_eq!(positions(&pdb), pdb_positions_before, "bits {label}");
    };
    for (cid, expected) in pdb_rows {
        pdb_check(&format!("pre {cid}"));
        let ids: Vec<u32> = pdb
            .selected_atom_ids(&BioSelection::from_cid(cid).unwrap())
            .unwrap()
            .iter()
            .map(|id| id.value())
            .collect();
        assert_eq!(ids, expected, "pdb {cid}");
        pdb_check(&format!("post {cid}"));
        // Exact slice pointer + length identity after EVERY query, against
        // both the pre-query state and the sharing peer.
        let after = pointers(&pdb);
        assert_eq!(after, pdb_pointers_before, "pointers {cid}");
        assert!(
            std::ptr::eq(after.0, pdb_peer.models().as_ptr()),
            "models {cid}"
        );
        assert!(
            std::ptr::eq(after.2, pdb_peer.chains().as_ptr()),
            "chains {cid}"
        );
        assert!(
            std::ptr::eq(after.4, pdb_peer.residues().as_ptr()),
            "residues {cid}"
        );
        assert!(
            std::ptr::eq(after.6, pdb_peer.atoms().as_ptr()),
            "atoms {cid}"
        );
        // Exact coordinate bits per existing atom ID after every query.
        assert_eq!(positions(&pdb), pdb_positions_before, "bits {cid}");
    }
    // Debug snapshot retained as a presentation check (it alone proves
    // neither pointer identity nor f64 bits — BIO-CID-SHARE correction).
    assert_eq!(pdb_before, format!("{pdb:?}"));
    pdb_check("pre repeat");
    assert_eq!(
        pdb.selected_atom_ids(&BioSelection::from_cid("A//SG[S]").unwrap())
            .unwrap()
            .len(),
        2
    );
    pdb_check("post repeat");
    assert_eq!(positions(&pdb), pdb_positions_before, "bits repeat");

    // Public mmCIF matrix: auth-precedence chain X, auth seqids 101/102.
    let cif = BioStructure::read_with_format(
        std::path::Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../testdata/bio/fixtures/gemmi_full_feature_sample.cif"
        )),
        cosmolkit::BioCoordinateFormat::Mmcif,
    )
    .unwrap();
    let cif_rows: [(&str, &[u32]); 5] = [
        ("X/101-102", &[0, 1]),
        ("X/1-2", &[]),
        ("X/(CYS)", &[0, 1]),
        (";polymer", &[0, 1]),
        (";solvent", &[]),
    ];
    let cif_peer = cif.clone();
    let cif_pointers_before = pointers(&cif);
    let cif_positions_before = positions(&cif);
    let cif_before = format!("{cif:?}"); // presentation check only
    // cif_check: the same complete checkpoint set as pdb_check, bound to
    // the CIF structure, its immutable peer and its fixed snapshots.
    let cif_check = |label: &str| {
        let current = pointers(&cif);
        let peer_now = pointers(&cif_peer);
        assert_eq!(
            current, cif_pointers_before,
            "cif original pointers {label}"
        );
        assert_eq!(current, peer_now, "cif peer tuple {label}");
        assert!(
            std::ptr::eq(current.0, cif_peer.models().as_ptr()),
            "cif models {label}"
        );
        assert!(
            std::ptr::eq(current.2, cif_peer.chains().as_ptr()),
            "cif chains {label}"
        );
        assert!(
            std::ptr::eq(current.4, cif_peer.residues().as_ptr()),
            "cif residues {label}"
        );
        assert!(
            std::ptr::eq(current.6, cif_peer.atoms().as_ptr()),
            "cif atoms {label}"
        );
        assert_eq!(positions(&cif), cif_positions_before, "cif bits {label}");
    };
    for (cid, expected) in cif_rows {
        cif_check(&format!("pre {cid}"));
        let ids: Vec<u32> = cif
            .selected_atom_ids(&BioSelection::from_cid(cid).unwrap())
            .unwrap()
            .iter()
            .map(|id| id.value())
            .collect();
        assert_eq!(ids, expected, "cif {cid}");
        cif_check(&format!("post {cid}"));
        let after = pointers(&cif);
        assert_eq!(after, cif_pointers_before, "cif pointers {cid}");
        assert!(
            std::ptr::eq(after.0, cif_peer.models().as_ptr()),
            "cif models {cid}"
        );
        assert!(
            std::ptr::eq(after.2, cif_peer.chains().as_ptr()),
            "cif chains {cid}"
        );
        assert!(
            std::ptr::eq(after.4, cif_peer.residues().as_ptr()),
            "cif residues {cid}"
        );
        assert!(
            std::ptr::eq(after.6, cif_peer.atoms().as_ptr()),
            "cif atoms {cid}"
        );
        assert_eq!(positions(&cif), cif_positions_before, "cif bits {cid}");
    }
    assert_eq!(cif_before, format!("{cif:?}"));

    // Registry: read-only instance delegate with the frozen projections.
    let entry = BINDING_CONTRACT
        .iter()
        .find(|e| e.semantic_id == "BioStructure.selected_atom_ids")
        .expect("selected_atom_ids registered");
    let callable = entry.callable.expect("callable contract");
    assert_eq!(callable.kind, BindingKind::Instance);
    assert_eq!(callable.state_model, StateModel::ReadOnly);
    assert_eq!(entry.python_name, "selected_atom_ids");
    assert_eq!(entry.javascript_name, "selectedAtomIds");
}

#[test]
fn bio_cid_c30_protein_projection_and_same_selector_behavior() {
    use cosmolkit::{BioSelection, BioStructure, Protein};

    // Protein projection over the mixed PDB fixture: the projected
    // hierarchy keeps the polymer rows only (CYS3 SG, CYS4 SG, ALA7 CA —
    // projected original IDs 0,1,2), so solvent/HETATM selectors exclude
    // everything while the same selectors keep their C29 semantics on
    // the full structure.
    let structure = BioStructure::read_with_format(
        std::path::Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../testdata/bio/fixtures/gemmi_full_feature_sample.pdb"
        )),
        cosmolkit::BioCoordinateFormat::Pdb,
    )
    .unwrap();
    let protein = structure.protein().unwrap();

    let protein_rows: [(&str, &[u32]); 7] = [
        ("A/*", &[0, 1, 2]),
        ("A//SG[S]", &[0, 1]),
        ("A/(ALA)", &[2]),
        ("A/;q>0.5", &[0, 1, 2]),
        ("A/(HOH)", &[]),  // excluded by the projection, not matched
        (";solvent", &[]), // no solvent rows remain in the projection
        ("A/;b<10", &[]),
    ];
    // Immutable cloned Protein sharing peer: the projection's Arc-backed
    // tables are shared with the clone, and every selector preserves the
    // four actual row-slice pointers/lengths and the coordinate bits
    // through the borrowed BioStructure (BIO-CID-SHARE).
    let peer = protein.clone();
    let pointers = |s: &BioStructure| {
        (
            s.models().as_ptr(),
            s.models().len(),
            s.chains().as_ptr(),
            s.chains().len(),
            s.residues().as_ptr(),
            s.residues().len(),
            s.atoms().as_ptr(),
            s.atoms().len(),
        )
    };
    let positions = |s: &BioStructure| -> Vec<Option<[u64; 3]>> {
        (0..s.atoms().len())
            .map(|index| {
                s.atom_position(cosmolkit::BioAtomId::new(index as u32))
                    .map(|xyz| [xyz[0].to_bits(), xyz[1].to_bits(), xyz[2].to_bits()])
            })
            .collect()
    };
    let before = format!("{protein:?}"); // presentation check only
    let pointers_before = pointers(protein.as_bio_structure());
    let positions_before = positions(protein.as_bio_structure());
    // Complete checkpoint closure (BIO-CID-SHARE-CLOSE): the four row
    // pointers AND lengths equal the ORIGINAL tuple, equal the PEER
    // tuple, four explicit peer ptr::eq identities, and the original
    // coordinate-bit snapshot — checked immediately BEFORE and AFTER
    // every one of the 21 query calls below (including the checkpoint
    // BETWEEN the Protein call and the borrowed-BioStructure call of
    // each parity comparison); originals are never refreshed, and this
    // proves persistent-storage preservation, not absence of transient
    // allocations.
    let check = |label: &str| {
        let current = pointers(protein.as_bio_structure());
        let peer_now = pointers(peer.as_bio_structure());
        assert_eq!(current, pointers_before, "original pointers {label}");
        assert_eq!(current, peer_now, "peer tuple {label}");
        assert!(
            std::ptr::eq(current.0, peer.as_bio_structure().models().as_ptr()),
            "models {label}"
        );
        assert!(
            std::ptr::eq(current.2, peer.as_bio_structure().chains().as_ptr()),
            "chains {label}"
        );
        assert!(
            std::ptr::eq(current.4, peer.as_bio_structure().residues().as_ptr()),
            "residues {label}"
        );
        assert!(
            std::ptr::eq(current.6, peer.as_bio_structure().atoms().as_ptr()),
            "atoms {label}"
        );
        assert_eq!(
            positions(protein.as_bio_structure()),
            positions_before,
            "bits {label}"
        );
    };
    for (cid, expected) in protein_rows {
        check(&format!("pre {cid}"));
        let ids: Vec<u32> = protein
            .selected_atom_ids(&BioSelection::from_cid(cid).unwrap())
            .unwrap()
            .iter()
            .map(|id| id.value())
            .collect();
        assert_eq!(ids, expected, "protein {cid}");
        check(&format!("post {cid}"));
    }
    assert_eq!(before, format!("{protein:?}"));

    // Same-selector behavior: the Protein delegate equals the query over
    // its OWN borrowed BioStructure for every selector above, with the
    // same pointer/bits identity preserved across each comparison.
    for (cid, _) in protein_rows {
        check(&format!("parity pre protein {cid}"));
        let via_protein: Vec<u32> = protein
            .selected_atom_ids(&BioSelection::from_cid(cid).unwrap())
            .unwrap()
            .iter()
            .map(|id| id.value())
            .collect();
        // Intermediate checkpoint BETWEEN the two calls of the comparison.
        check(&format!("parity between {cid}"));
        let via_structure: Vec<u32> = protein
            .as_bio_structure()
            .selected_atom_ids(&BioSelection::from_cid(cid).unwrap())
            .unwrap()
            .iter()
            .map(|id| id.value())
            .collect();
        assert_eq!(via_protein, via_structure, "parity {cid}");
        check(&format!("parity post structure {cid}"));
    }
    // Mixed-structure exclusion contrast: the FULL structure still sees
    // the HOH row under the same selector (C29 semantics preserved).
    assert_eq!(
        structure
            .selected_atom_ids(&BioSelection::from_cid("A/(HOH)").unwrap())
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn bio_cid_c26_wrapper_is_thin_and_private_fielded() {
    // The wrapper is a projection type around the canonical detached
    // data: size identity proves no extra state beyond the inner data.
    use cosmolkit::BioSelection;
    assert_eq!(
        std::mem::size_of::<BioSelection>(),
        std::mem::size_of::<cosmolkit_bio::BioSelectionData>()
    );
}
