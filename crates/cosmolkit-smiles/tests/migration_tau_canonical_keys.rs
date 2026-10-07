use cosmolkit_smiles::{
    SmilesParseParams, SmilesRecord, finalize_smiles_stereo, parse_smiles, write_smiles,
};
// Original source-text fixtures decode only at this observation boundary.
// Invalid UTF-8 fails; the complete byte payload is never substituted.
fn fixture_writer_text(text: cosmolkit_model::PropertyText) -> String {
    String::from_utf8(text.into_bytes()).expect("original writer fixture UTF-8 bytes")
}

#[test]
fn source_tautomer_canonical_keys_retain_the_pinned_fragment_order() {
    // Pinned TAU focused expected default rows: these are already canonical
    // keys, so canonical serialization of the same graph must retain them.
    for text in [
        "O=C(O)C(CO)=NC(=O)C1=C(O)C(O)CC=C1",
        "O=C(O)C(CO)=NC(=O)C1=CCCC(O)=C1O",
    ] {
        let params = SmilesParseParams::default();
        let record = parse_smiles(text, &params).unwrap();
        let result = cosmolkit_core::remove_hydrogens_with_params(
            record.topology,
            record.coordinates,
            record.properties,
            &cosmolkit_core::RemoveHsParams {
                sanitize: true,
                update_explicit_count: true,
                ..Default::default()
            },
        )
        .unwrap();
        let mut valence = result.final_valence;
        let mut rings = result.final_rings;
        let record = finalize_smiles_stereo(
            SmilesRecord {
                topology: result.topology,
                coordinates: result.coordinates,
                properties: result.properties,
            },
            &params,
            &mut valence,
            &mut rings,
        )
        .unwrap();
        assert_eq!(
            write_smiles(&record).map(fixture_writer_text).unwrap(),
            text
        );
    }
}
