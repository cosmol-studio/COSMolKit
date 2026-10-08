use super::registry::{Input, Record, Value};

pub fn run(input: &Input) -> Result<Record, String> {
    match input {
        Input::MolAlign(row) => crate::molalign::run(row),
        Input::Fingerprint(row) => crate::fingerprints::run(row),
        Input::SmilesWrite(row) => crate::smiles_write::run(row),
        Input::Search(row) => crate::search::run(row),
        Input::Uff(_) => crate::uff::run(input),
        Input::Mmff(row) => crate::mmff::run(row),
        Input::PersistentForceField(row) => crate::persistent_forcefields::run(row),
        Input::Molecular { .. } => crate::molecular::run(input),
        Input::BioPdbOutput { case, profile } => {
            let structure = match case.format {
                crate::registry::BioPdbCorpusFormat::Cif => {
                    cosmolkit::BioStructure::from_mmcif(&case.text)
                        .map_err(|error| error.to_string())
                }
                crate::registry::BioPdbCorpusFormat::Pdb => {
                    cosmolkit::BioStructure::from_pdb(&case.text).map_err(|error| error.to_string())
                }
            };
            let structure = match structure {
                Ok(structure) => structure,
                Err(error) => {
                    return Ok(Record {
                        input: input.clone(),
                        output: Value::BioPdbOutput(crate::registry::BioPdbOutputValue {
                            text: String::new(),
                            error: Some(crate::registry::BioPdbOutputError::Parse {
                                format: case.format,
                                message: error.to_string(),
                            }),
                        }),
                    });
                }
            };
            let params = cosmolkit::BioPdbWriteParams {
                ter_records: profile.ter_records,
                numbered_ter: profile.numbered_ter,
                ter_ignores_type: profile.ter_ignores_type,
                preserve_serial: profile.preserve_serial,
                end_record: profile.end_record,
            };
            let (text, error) = match structure.to_pdb_with_params(&params) {
                Ok(t) => (t, None),
                Err(e) => (
                    String::new(),
                    Some(crate::registry::BioPdbOutputError::Write {
                        message: e.to_string(),
                    }),
                ),
            };
            Ok(Record {
                input: input.clone(),
                output: Value::BioPdbOutput(super::registry::BioPdbOutputValue { text, error }),
            })
        }
    }
}
