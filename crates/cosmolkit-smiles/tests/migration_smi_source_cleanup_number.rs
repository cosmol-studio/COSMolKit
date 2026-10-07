use cosmolkit_smiles::{SmilesParseParams, parse_smiles};
#[test]
fn pinned_atom_map_number_rejects_conservative_signed_overflow_boundary() {
    for value in [
        "2147483640",
        "2147483646",
        "2147483647",
        "2147483648",
        "4294967295",
    ] {
        assert!(
            parse_smiles(&format!("[13CH3:{value}]"), &SmilesParseParams::default()).is_err(),
            "{value}"
        );
    }
    for value in [0_u32, 7, 2147483639] {
        let record =
            parse_smiles(&format!("[13CH3:{value}]"), &SmilesParseParams::default()).unwrap();
        assert_eq!(record.topology.atoms[0].atom_map(), Some(value));
    }
}
#[test]
fn pinned_smiles_cleanup_preserves_original_sgroup_cx_index() {
    for text in [
        "CC |SgD:0:FIELD:DATA:QUERY:INFO:TAG:|",
        "CC |SgD:0:FIELD:DATA::::|",
    ] {
        let record = parse_smiles(text, &SmilesParseParams::default()).unwrap();
        assert_eq!(record.topology.substance_groups.len(), 1);
        assert_eq!(
            record.topology.substance_groups[0]
                .props()
                .get("_cxsmilesindex".as_bytes())
                .unwrap(),
            &cosmolkit_model::PropertyValue::UInt(0),
            "{text}"
        );
    }
}
