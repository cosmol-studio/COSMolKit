use cosmolkit::{
    AltLocLabel, AtomName, AtomSourceIds, BINDING_CONTRACT, BindingItem, BindingKind, BindingOwner,
    BindingTypeRole, ChainSourceIds, EntitySourceIds, FunctionStatus, PdbAtomSerial, PdbChainId,
    PdbSeqId, ResidueCode, ResidueInfo, ResidueInfoKind, ResidueName, ResidueSequenceError,
    ResidueSourceIds, StateModel, UNKNOWN_TABULATED_RESIDUE_INDEX, expand_one_letter,
    expand_one_letter_sequence, find_residue_info, find_residue_info_index, residue_code,
    residue_info, residue_info_checked,
};

fn compact(value: &str) -> String {
    value
        .chars()
        .filter(|character| !character.is_whitespace())
        .collect()
}

#[test]
fn canonical_residue_function_signatures_are_public() {
    let _: fn(usize) -> ResidueInfo = residue_info;
    let _: fn(usize) -> Option<ResidueInfo> = residue_info_checked;
    let _: fn(&str) -> usize = find_residue_info_index;
    let _: fn(&str) -> ResidueInfo = find_residue_info;
    let _: fn(&str) -> ResidueCode = residue_code;
    let _: fn(char, ResidueInfoKind) -> Option<&'static str> = expand_one_letter;
    let _: fn(&str, ResidueInfoKind) -> Result<Vec<String>, ResidueSequenceError> =
        expand_one_letter_sequence;
}

/// BIO-PARENT1-44: the three restored legacy metadata helpers keep their
/// frozen signatures on the public facade value, produce the exact eleven
/// frozen source-row outputs, and their registry metadata matches the
/// declared feature/output/error/receiver/status/state contract.
#[test]
fn legacy_residue_metadata_helpers_are_public_with_frozen_outputs() {
    use cosmolkit::FunctionStatus;

    // Signature controls on the public facade re-export.
    let canonical: fn(ResidueInfo) -> Option<char> = ResidueInfo::canonical_one_letter_code;
    let parent: fn(ResidueInfo) -> Option<ResidueCode> = ResidueInfo::parent_standard_code;
    let modified: fn(ResidueInfo) -> bool = ResidueInfo::is_modified_amino_acid;

    // Frozen eleven source rows (index -> expected triple), annex-derived.
    const ROWS: [(usize, &str, Option<char>, Option<ResidueCode>, bool); 11] = [
        (17, "MSE", Some('M'), Some(ResidueCode::MET), true),
        (45, "HYP", Some('P'), Some(ResidueCode::PRO), true),
        (29, "SEP", Some('S'), Some(ResidueCode::SER), true),
        (20, "PRO", Some('P'), Some(ResidueCode::PRO), false),
        (97, "3FG", None, None, true),
        (154, "HOH", None, None, false),
        (5, "ASX", Some('B'), Some(ResidueCode::ASX), false),
        (25, "UNK", Some('X'), Some(ResidueCode::UNK), false),
        (27, "SEC", Some('U'), Some(ResidueCode::SEC), false),
        (28, "PYL", Some('O'), Some(ResidueCode::PYL), false),
        (39, "DAL", Some('A'), Some(ResidueCode::ALA), true),
    ];
    for (index, name, expected_char, expected_parent, expected_modified) in ROWS {
        let info = residue_info(index);
        assert_eq!(info.name, name, "{name}: row identity");
        assert_eq!(canonical(info), expected_char, "{name}: canonical");
        assert_eq!(parent(info), expected_parent, "{name}: parent");
        assert_eq!(modified(info), expected_modified, "{name}: modified");
    }

    // Registry metadata controls: exactly the three frozen entries with
    // feature/output/error/owned-receiver/status/state declarations.
    let contract = BINDING_CONTRACT;
    for (id, python, javascript) in [
        (
            "ResidueInfo.canonical_one_letter_code",
            "canonical_one_letter_code",
            "canonicalOneLetterCode",
        ),
        (
            "ResidueInfo.parent_standard_code",
            "parent_standard_code",
            "parentStandardCode",
        ),
        (
            "ResidueInfo.is_modified_amino_acid",
            "is_modified_amino_acid",
            "isModifiedAminoAcid",
        ),
    ] {
        let entry = contract
            .iter()
            .find(|entry| entry.semantic_id == id)
            .unwrap_or_else(|| panic!("missing registry entry {id}"));
        assert_eq!(entry.python_name, python, "{id}: python");
        assert_eq!(entry.javascript_name, javascript, "{id}: javascript");
        assert_eq!(entry.feature, "cap-bio", "{id}: feature");
        assert_eq!(entry.status, FunctionStatus::Experimental, "{id}: status");
        let callable = entry
            .callable
            .as_ref()
            .unwrap_or_else(|| panic!("{id}: entry must be a real callable, not a type entry"));
        assert!(
            matches!(callable.kind, cosmolkit::BindingKind::Instance),
            "{id}: kind"
        );
        assert!(
            matches!(callable.receiver, Some(cosmolkit::BindingReceiver::Owned)),
            "{id}: owned receiver"
        );
        assert!(callable.parameters.is_empty(), "{id}: parameters");
        assert!(callable.error_type.is_none(), "{id}: no error");
        assert!(
            matches!(callable.state_model, cosmolkit::StateModel::ValueReturning),
            "{id}: value_returning state"
        );
        assert!(
            callable.operation_semantic_id.is_none(),
            "{id}: no operation"
        );
    }
}

/// BIO-KIND1-28: the restored legacy residue-kind name helper is public on
/// the facade value with const evaluation, the frozen four literal
/// controls, standard Display, and its real binding entry/metadata/ID
/// placement (no owner-suite duplication, no bindings claim).
#[test]
fn legacy_residue_kind_name_is_public_with_frozen_controls() {
    use cosmolkit::FunctionStatus;

    // Signature + const-evaluation controls.
    let name_fn: fn(ResidueInfoKind) -> &'static str = ResidueInfoKind::name;
    const _CONST_PROOF: &str = ResidueInfoKind::Hoh.name();
    assert_eq!(_CONST_PROOF, "HOH");

    // Frozen four literal controls.
    assert_eq!(name_fn(ResidueInfoKind::Unknown), "UNKNOWN");
    assert_eq!(name_fn(ResidueInfoKind::Aa), "AA");
    assert_eq!(name_fn(ResidueInfoKind::Maa), "MAA");
    assert_eq!(name_fn(ResidueInfoKind::Hoh), "HOH");

    // Standard Display through the ONE name owner.
    assert_eq!(format!("{}", ResidueInfoKind::Els), "ELS");

    // Real binding entry + metadata + ID placement.
    let contract = BINDING_CONTRACT;
    let entry = contract
        .iter()
        .find(|entry| entry.semantic_id == "ResidueInfoKind.name")
        .expect("missing registry entry ResidueInfoKind.name");
    assert_eq!(entry.python_name, "name", "python");
    assert_eq!(entry.javascript_name, "name", "javascript");
    assert_eq!(entry.feature, "cap-bio", "feature");
    assert_eq!(entry.status, FunctionStatus::Experimental, "status");
    let callable = entry
        .callable
        .as_ref()
        .expect("entry must be a real callable");
    assert!(
        matches!(callable.kind, cosmolkit::BindingKind::Instance),
        "kind"
    );
    assert!(
        matches!(callable.receiver, Some(cosmolkit::BindingReceiver::Owned)),
        "owned receiver"
    );
    assert!(callable.parameters.is_empty(), "parameters");
    assert_eq!(callable.output_type, "& 'static str", "output");
    assert!(callable.error_type.is_none(), "no error");
    assert!(
        matches!(callable.state_model, cosmolkit::StateModel::ValueReturning),
        "value_returning state"
    );
    assert!(callable.operation_semantic_id.is_none(), "no operation");
}

/// BIO-PARSE1-36: the restored legacy permissive ResidueCode text parser
/// is public on the facade value with the exact FromStr error type, the
/// input accessor fn-pointer, six literal successes, two literal errors,
/// a terminated error source chain, and the real binding entries for the
/// parse-error type and its accessor at their actual registry position.
#[test]
fn residue_identity_and_info_serde_are_public_with_frozen_wire() {
    use crate::compact;
    use cosmolkit::{ResidueCode, ResidueIdentity};

    // Exact signatures through fn pointers, including the const accessors
    // (const fns coerce to the same fn-pointer shape) and the generic
    // `new` selected by its registered monomorphic instantiation.
    let new_fn: fn(String) -> ResidueIdentity = ResidueIdentity::new;
    let name_fn: fn(&ResidueIdentity) -> &str = ResidueIdentity::name;
    let code_fn: fn(&ResidueIdentity) -> ResidueCode = ResidueIdentity::code;
    let info_fn: fn(&ResidueIdentity) -> cosmolkit::ResidueInfo = ResidueIdentity::info;
    let tabulated_fn: fn(&ResidueIdentity) -> bool = ResidueIdentity::is_tabulated;

    let alias = new_fn("wat".to_string());
    assert_eq!(name_fn(&alias), "wat");
    assert_eq!(u16::from(code_fn(&alias).as_u16()), 154);
    assert_eq!(info_fn(&alias).code, code_fn(&alias));
    assert!(tabulated_fn(&alias));

    let unknown = new_fn("CUSTOM".to_string());
    assert_eq!(name_fn(&unknown), "CUSTOM");
    assert_eq!(u16::from(code_fn(&unknown).as_u16()), 367);
    assert!(!tabulated_fn(&unknown));

    // Raw unknown and alias round trips through the facade.
    for (raw, ordinal) in [("CUSTOM", 367u16), ("wat", 154), ("al", 367)] {
        let identity = new_fn(raw.to_string());
        let wire = serde_json::to_string(&identity).unwrap();
        assert_eq!(wire, serde_json::to_string(raw).unwrap());
        let round: ResidueIdentity = serde_json::from_str(&wire).unwrap();
        assert_eq!(round.name(), raw);
        assert_eq!(u16::from(round.code().as_u16()), ordinal);
    }

    // ResidueInfo seven-field example through the facade type.
    let value = serde_json::to_value(info_fn(&new_fn("MSE".to_string()))).unwrap();
    assert_eq!(value["code"], "MSE");
    assert_eq!(value["name"], "MSE");
    assert_eq!(value["kind"], "AA");
    assert_eq!(value["linking_type"], 1);
    assert_eq!(value["one_letter_code"], "m");
    assert_eq!(value["hydrogen_count"], 11);
    assert_eq!(value.as_object().map(|m| m.len()), Some(7));

    // Registry rows: exact feature/status/kind/receiver/output/state.
    let registry = BINDING_CONTRACT;
    let type_row = registry
        .iter()
        .find(|entry| entry.semantic_id == "types.ResidueIdentity")
        .expect("missing registry entry types.ResidueIdentity");
    assert_eq!(type_row.item, BindingItem::Type);
    assert_eq!(type_row.feature, "cap-bio");
    assert_eq!(type_row.status, FunctionStatus::Experimental);
    assert_eq!(type_row.type_role, Some(BindingTypeRole::Value));
    assert!(type_row.callable.is_none());

    let expect_instance = |semantic_id: &str, output: &str| {
        let row = registry
            .iter()
            .find(|entry| entry.semantic_id == semantic_id)
            .unwrap_or_else(|| panic!("missing {semantic_id}"));
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.feature, "cap-bio");
        assert_eq!(row.status, FunctionStatus::Experimental);
        let callable = row.callable.as_ref().unwrap();
        assert_eq!(callable.kind, BindingKind::Instance);
        assert!(matches!(
            callable.receiver,
            Some(cosmolkit::BindingReceiver::Shared)
        ));
        assert_eq!(compact(output), compact(&callable.output_type));
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
    };
    expect_instance("ResidueIdentity.name", "&str");
    expect_instance("ResidueIdentity.code", "crate::ResidueCode");
    expect_instance("ResidueIdentity.info", "crate::ResidueInfo");
    expect_instance("ResidueIdentity.is_tabulated", "bool");

    let new_row = registry
        .iter()
        .find(|entry| entry.semantic_id == "ResidueIdentity.new")
        .expect("missing registry entry ResidueIdentity.new");
    assert_eq!(new_row.item, BindingItem::Callable);
    assert_eq!(new_row.feature, "cap-bio");
    assert_eq!(new_row.status, FunctionStatus::Experimental);
    let callable = new_row.callable.as_ref().unwrap();
    assert_eq!(callable.kind, BindingKind::Static);
    assert_eq!(
        compact("crate::ResidueIdentity"),
        compact(&callable.output_type)
    );
    assert_eq!(callable.state_model, StateModel::ValueReturning);
    assert_eq!(callable.operation_semantic_id, None);

    // Frozen registry order: the six identity rows precede types.ResidueInfo.
    let ids: Vec<&str> = registry.iter().map(|entry| &*entry.semantic_id).collect();
    let base = ids
        .iter()
        .position(|id| *id == "types.ResidueIdentity")
        .expect("types.ResidueIdentity present");
    assert_eq!(
        &ids[base..base + 6],
        [
            "types.ResidueIdentity",
            "ResidueIdentity.new",
            "ResidueIdentity.name",
            "ResidueIdentity.code",
            "ResidueIdentity.info",
            "ResidueIdentity.is_tabulated",
        ]
    );
    assert_eq!(ids[base + 6], "types.ResidueInfo");
}

#[test]
fn legacy_residue_code_from_str_is_public_with_frozen_controls() {
    use cosmolkit::{BindingKind, BindingReceiver, ResidueCode, ResidueCodeParseError, StateModel};
    use std::str::FromStr;

    // Exact FromStr Err type + input fn-pointer signature.
    fn assert_from_str_err_is_parse_error()
    where
        ResidueCode: FromStr<Err = ResidueCodeParseError>,
    {
    }
    assert_from_str_err_is_parse_error();
    let input_fn: fn(&ResidueCodeParseError) -> &str = ResidueCodeParseError::input;

    // Six literal successes ALA/ala/H2O/TRY/UNKNOWN/UNK.
    assert_eq!(ResidueCode::from_str("ALA").map(|c| c as u16), Ok(0));
    assert_eq!(ResidueCode::from_str("ala").map(|c| c as u16), Ok(0));
    assert_eq!(ResidueCode::from_str("H2O").map(|c| c as u16), Ok(154));
    assert_eq!(ResidueCode::from_str("TRY").map(|c| c as u16), Ok(23));
    assert_eq!(ResidueCode::from_str("UNKNOWN").map(|c| c as u16), Ok(367));
    assert_eq!(ResidueCode::from_str("UNK").map(|c| c as u16), Ok(25));

    // Two literal errors " ALA" and "ALA\0" with full error proofs.
    for input in [" ALA", "ALA\0"] {
        let error = ResidueCode::from_str(input)
            .err()
            .unwrap_or_else(|| panic!("{input:?}: expected Err"));
        assert_eq!(input_fn(&error), input, "{input:?}: exact input bytes");
        assert_eq!(
            error.to_string(),
            format!("unknown residue code name '{input}'"),
            "{input:?}: literal Display"
        );
        let dynamic: &dyn std::error::Error = &error;
        assert!(dynamic.source().is_none(), "{input:?}: source None");
    }

    // Real binding contracts and registry order for the two new entries.
    let contract = BINDING_CONTRACT;
    let type_entry = contract
        .iter()
        .find(|entry| entry.semantic_id == "types.ResidueCodeParseError")
        .expect("missing registry entry types.ResidueCodeParseError");
    assert_eq!(type_entry.python_name, "ResidueCodeParseError", "python");
    assert_eq!(
        type_entry.javascript_name, "ResidueCodeParseError",
        "javascript"
    );
    assert_eq!(type_entry.feature, "cap-bio", "feature");
    assert_eq!(
        type_entry.status,
        cosmolkit::FunctionStatus::Experimental,
        "status"
    );
    assert!(
        matches!(type_entry.item, cosmolkit::BindingItem::Type),
        "type item"
    );

    let callable_entry = contract
        .iter()
        .find(|entry| entry.semantic_id == "ResidueCodeParseError.input")
        .expect("missing registry entry ResidueCodeParseError.input");
    assert_eq!(callable_entry.python_name, "input", "python");
    assert_eq!(callable_entry.javascript_name, "input", "javascript");
    assert_eq!(callable_entry.feature, "cap-bio", "feature");
    assert_eq!(
        callable_entry.status,
        cosmolkit::FunctionStatus::Experimental,
        "status"
    );
    let callable = callable_entry
        .callable
        .as_ref()
        .expect("entry must be a real callable");
    assert!(
        matches!(callable.kind, BindingKind::Instance),
        "instance kind"
    );
    assert!(
        matches!(callable.receiver, Some(BindingReceiver::Shared)),
        "shared receiver"
    );
    assert!(callable.parameters.is_empty(), "parameters");
    assert_eq!(callable.output_type, "& str", "output");
    assert!(callable.error_type.is_none(), "no error");
    assert!(
        matches!(callable.state_model, StateModel::ReadOnly),
        "read_only state"
    );
    assert!(callable.operation_semantic_id.is_none(), "no operation");

    // Registry order: both entries sit immediately after types.ResidueCode,
    // before types.ResidueInfo, preserving every prior relative order.
    let ids: Vec<&str> = contract.iter().map(|e| e.semantic_id).collect();
    let code_pos = ids
        .iter()
        .position(|id| *id == "types.ResidueCode")
        .expect("types.ResidueCode present");
    let parse_error_pos = ids
        .iter()
        .position(|id| *id == "types.ResidueCodeParseError")
        .expect("types.ResidueCodeParseError present");
    let input_pos = ids
        .iter()
        .position(|id| *id == "ResidueCodeParseError.input")
        .expect("ResidueCodeParseError.input present");
    let info_pos = ids
        .iter()
        .position(|id| *id == "types.ResidueInfo")
        .expect("types.ResidueInfo present");
    assert_eq!(
        parse_error_pos,
        code_pos + 1,
        "type entry after ResidueCode"
    );
    assert_eq!(input_pos, parse_error_pos + 1, "input after its type");
    assert_eq!(
        info_pos,
        input_pos + 7,
        "ResidueInfo follows the six authorized identity rows"
    );
}

#[test]
fn public_residue_lookup_preserves_known_unknown_and_checked_boundaries() {
    let alanine = residue_info(0);
    assert_eq!(alanine.code, ResidueCode::ALA);
    assert_eq!(alanine.name, "ALA");
    assert_eq!(alanine.kind, ResidueInfoKind::Aa);
    assert_eq!(alanine.one_letter_code, 'A');
    assert!(alanine.is_amino_acid());
    assert!(alanine.is_peptide_linking());

    assert_eq!(find_residue_info_index("ALA"), 0);
    assert_eq!(find_residue_info_index("ala"), 0);
    assert_eq!(find_residue_info("TRY"), residue_info(23));
    assert_eq!(residue_code("H2O"), ResidueCode::HOH);

    let unknown = residue_info(UNKNOWN_TABULATED_RESIDUE_INDEX);
    assert_eq!(unknown.code, ResidueCode::UNKNOWN);
    assert!(!unknown.found());
    assert_eq!(find_residue_info("not-a-residue"), unknown);
    assert_eq!(
        residue_info_checked(UNKNOWN_TABULATED_RESIDUE_INDEX),
        Some(unknown)
    );
    assert_eq!(
        residue_info_checked(UNKNOWN_TABULATED_RESIDUE_INDEX + 1),
        None
    );
    assert_eq!(residue_info_checked(usize::MAX), None);
}

#[test]
fn public_one_letter_expansion_preserves_order_special_tokens_and_errors() {
    assert_eq!(expand_one_letter('a', ResidueInfoKind::Aa), Some("ALA"));
    assert_eq!(expand_one_letter('T', ResidueInfoKind::Dna), Some("DT"));
    assert_eq!(expand_one_letter('T', ResidueInfoKind::Rna), None);
    assert_eq!(expand_one_letter('J', ResidueInfoKind::Aa), None);
    assert_eq!(
        expand_one_letter_sequence(" A\tC\n(MSE) ", ResidueInfoKind::Aa).unwrap(),
        ["ALA", "CYS", "MSE"]
    );
    assert_eq!(
        expand_one_letter_sequence("A(MSE", ResidueInfoKind::Aa),
        Err(ResidueSequenceError::UnmatchedParenthesis)
    );
    assert_eq!(
        expand_one_letter_sequence("AJ", ResidueInfoKind::Aa),
        Err(ResidueSequenceError::UnexpectedLetter {
            kind: "peptide",
            letter: 'J',
            source_code: 74,
        })
    );
}

#[test]
fn source_identifier_values_are_usable_through_the_facade() {
    let atom_name = AtomName::from_ascii(b" CA ").unwrap();
    let cif_atom_name = AtomName::from_ascii(b"CA").unwrap();
    let residue_name = ResidueName::from_ascii(b"MSE ").unwrap();
    let auth_chain = PdbChainId::from_ascii(b"A").unwrap();
    let sequence = PdbSeqId::new(-12, Some(b'B'));
    let serial = PdbAtomSerial::new(17);
    let altloc = AltLocLabel::new(b'A');
    let atom_source = AtomSourceIds::new(Some(serial));
    let residue_source = ResidueSourceIds::new(
        Some(sequence),
        Some(42),
        Some(*b"SEG "),
        Some("subchain".to_owned()),
        Some("entity".to_owned()),
    )
    .unwrap();
    let chain_source = ChainSourceIds::new(Some(auth_chain), Some("label_asym".to_owned()));
    let entity_source = EntitySourceIds::new("entity".to_owned());

    assert_eq!(atom_name.as_str(), " CA ");
    assert_eq!(cif_atom_name.as_bytes(), b"CA");
    assert_ne!(cif_atom_name, atom_name);
    assert_eq!(residue_name.as_str(), "MSE ");
    assert_eq!(sequence.seq_num(), -12);
    assert_eq!(sequence.ins_code(), Some(b'B'));
    assert_eq!(serial.value(), 17);
    assert_eq!(altloc.value(), b'A');
    assert_eq!(atom_source.serial(), Some(serial));
    assert_eq!(residue_source.seq_id(), Some(sequence));
    assert_eq!(residue_source.label_seq_id(), Some(42));
    assert_eq!(chain_source.auth_chain_id(), Some(auth_chain));
    assert_eq!(chain_source.label_asym_id(), Some("label_asym"));
    assert_eq!(entity_source.source_entity_id(), "entity");
}

#[test]
fn binding_contract_matches_the_public_bio_residue_surface() {
    let type_rows = [
        (
            "types.ResidueInfoKind",
            "ResidueInfoKind",
            BindingTypeRole::Value,
        ),
        ("types.ResidueCode", "ResidueCode", BindingTypeRole::Value),
        ("types.ResidueInfo", "ResidueInfo", BindingTypeRole::Value),
        (
            "types.PdbAtomSerial",
            "PdbAtomSerial",
            BindingTypeRole::Value,
        ),
        ("types.PdbChainId", "PdbChainId", BindingTypeRole::Value),
        ("types.PdbSeqId", "PdbSeqId", BindingTypeRole::Value),
        ("types.AtomName", "AtomName", BindingTypeRole::Value),
        ("types.ResidueName", "ResidueName", BindingTypeRole::Value),
        ("types.AltLocLabel", "AltLocLabel", BindingTypeRole::Value),
        (
            "types.AtomSourceIds",
            "AtomSourceIds",
            BindingTypeRole::Value,
        ),
        (
            "types.ResidueSourceIds",
            "ResidueSourceIds",
            BindingTypeRole::Value,
        ),
        (
            "types.ChainSourceIds",
            "ChainSourceIds",
            BindingTypeRole::Value,
        ),
        (
            "types.EntitySourceIds",
            "EntitySourceIds",
            BindingTypeRole::Value,
        ),
        (
            "types.ResidueSequenceError",
            "ResidueSequenceError",
            BindingTypeRole::Error,
        ),
    ];
    for (semantic_id, rust_name, role) in type_rows {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == semantic_id)
            .unwrap();
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.type_role, Some(role));
        assert_eq!(compact(row.rust_path), format!("crate::{rust_name}"));
        assert_eq!(row.python_name, rust_name);
        assert_eq!(row.javascript_name, rust_name);
        assert_eq!(row.feature, "cap-bio");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert!(row.callable.is_none());
    }

    let callable_rows = [
        (
            "module.residue_info",
            "residue_info",
            "residueInfo",
            1,
            "crate::ResidueInfo",
            None,
        ),
        (
            "module.residue_info_checked",
            "residue_info_checked",
            "residueInfoChecked",
            1,
            "Option<crate::ResidueInfo>",
            None,
        ),
        (
            "module.find_residue_info_index",
            "find_residue_info_index",
            "findResidueInfoIndex",
            1,
            "usize",
            None,
        ),
        (
            "module.find_residue_info",
            "find_residue_info",
            "findResidueInfo",
            1,
            "crate::ResidueInfo",
            None,
        ),
        (
            "module.residue_code",
            "residue_code",
            "residueCode",
            1,
            "crate::ResidueCode",
            None,
        ),
        (
            "module.expand_one_letter",
            "expand_one_letter",
            "expandOneLetter",
            2,
            "Option<&'staticstr>",
            None,
        ),
        (
            "module.expand_one_letter_sequence",
            "expand_one_letter_sequence",
            "expandOneLetterSequence",
            2,
            "Vec<String>",
            Some("crate::ResidueSequenceError"),
        ),
    ];
    for (semantic_id, rust_name, javascript, parameter_count, output, error) in callable_rows {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == semantic_id)
            .unwrap();
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Module);
        assert_eq!(compact(row.rust_path), format!("crate::{rust_name}"));
        assert_eq!(row.python_name, rust_name);
        assert_eq!(row.javascript_name, javascript);
        assert_eq!(row.feature, "cap-bio");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert!(row.type_role.is_none());
        let callable = row.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Module);
        assert_eq!(callable.parameters.len(), parameter_count);
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
        assert_eq!(compact(callable.output_type), output);
        assert_eq!(callable.error_type.map(compact), error.map(str::to_owned));
    }
}
