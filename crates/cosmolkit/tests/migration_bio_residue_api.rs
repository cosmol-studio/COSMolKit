use cosmolkit::{
    AltLocLabel, AtomName, AtomSourceIds, BINDING_CONTRACT, BindingExposure, BindingItem,
    BindingKind, BindingOwner, BindingParity, BindingSupport, BindingTypeRole, ChainSourceIds,
    EntitySourceIds, PdbAtomSerial, PdbChainId, PdbSeqId, ResidueCode, ResidueInfo,
    ResidueInfoKind, ResidueName, ResidueSequenceError, ResidueSourceIds, StateModel,
    UNKNOWN_TABULATED_RESIDUE_INDEX, expand_one_letter, expand_one_letter_sequence,
    find_residue_info, find_residue_info_index, residue_code, residue_info, residue_info_checked,
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
    let atom_name = AtomName::from_ascii(*b" CA ").unwrap();
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
        assert_eq!(row.feature, "bio");
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
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
        assert_eq!(row.feature, "bio");
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
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
