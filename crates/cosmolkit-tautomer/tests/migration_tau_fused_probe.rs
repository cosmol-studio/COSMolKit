use cosmolkit_tautomer::*;
use serde_json::json;
use std::{fs::File, io::Write, path::Path, sync::Mutex};
struct Recorder {
    state: Mutex<(usize, File)>,
}
impl TautomerEnumerationCallback for Recorder {
    fn should_continue(
        &self,
        _: TautomerRecordView<'_>,
        result: TautomerProgress<'_>,
    ) -> Result<bool, TautomerRunError> {
        let mut state = self
            .state
            .lock()
            .map_err(|e| TautomerRunError::Callback(e.to_string()))?;
        let row = json!({"step":state.0,"keys":result.entries().map(|(k,_)|k).collect::<Vec<_>>(),"modified_atoms":result.modified_atoms().iter().map(|a|a.index()).collect::<Vec<_>>(),"modified_bonds":result.modified_bonds().iter().map(|b|b.index()).collect::<Vec<_>>()});
        writeln!(&mut state.1, "{}", row).map_err(|e| TautomerRunError::Callback(e.to_string()))?;
        state.0 += 1;
        Ok(true)
    }
}
#[test]
fn original_pcs_fused_ring_boundary_with_preapplication_probe() {
    let input = "Cc1nc2c(nc1C)C(=O)C1=C(C2=O)C2C=CC1CC2";
    let params = cosmolkit_smiles::SmilesParseParams::default();
    let parsed = cosmolkit_smiles::parse_smiles(input, &params).unwrap();
    let removed = cosmolkit_core::remove_hydrogens_with_params(
        parsed.topology,
        parsed.coordinates,
        parsed.properties,
        &cosmolkit_core::RemoveHsParams {
            update_explicit_count: true,
            sanitize: true,
            ..Default::default()
        },
    )
    .unwrap();
    let mut valence = removed.final_valence;
    let mut rings = removed.final_rings;
    let source = cosmolkit_smiles::finalize_smiles_stereo(
        cosmolkit_smiles::SmilesRecord {
            topology: removed.topology,
            coordinates: removed.coordinates,
            properties: removed.properties,
        },
        &params,
        &mut valence,
        &mut rings,
    )
    .unwrap();
    let before = source.clone();
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/TAU-fused-probe1-rust-callbacks.jsonl");
    let recorder = Recorder {
        state: Mutex::new((0, File::create(path).unwrap())),
    };
    let result = enumerate_with_catalog(
        TautomerRecordView {
            topology: &source.topology,
            coordinates: &source.coordinates,
            properties: &source.properties,
            valence: valence.as_ref(),
            rings: rings.as_ref(),
        },
        &TautomerCatalog::current().unwrap(),
        TautomerParams::default(),
        Some(&recorder),
    )
    .unwrap();
    assert_eq!(source, before);
    assert_eq!(
        result.status,
        TautomerEnumerationStatus::MaxTransformsReached
    );
    assert_eq!(result.entries.len(), 272);
}
