use cosmolkit_bio::BioCoordinateFormat as F;
use cosmolkit_io::{BioReadError as E, BioReadParams, read_bio_structure, read_bio_structure_file};
use std::error::Error;

const INPUTS: [(&str, &str, F); 5] = [
    ("pdb", "HEADER    TEST\n", F::Pdb),
    ("cif", "data_demo\n_entry.id DEMO\n", F::Mmcif),
    (
        "json",
        r#"{"data_demo":{"entry":{"id":["DEMO"]}}}"#,
        F::Mmjson,
    ),
    (
        "chemcomp",
        "data_comp_CMP\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n",
        F::Mmcif,
    ),
    ("invalid", "", F::Unknown),
];
const FORMATS: [F; 6] = [
    F::Unknown,
    F::Detect,
    F::Pdb,
    F::Mmcif,
    F::Mmjson,
    F::ChemComp,
];

#[derive(Clone, Copy, Debug)]
enum Outcome {
    Success(F, usize),
    Wrong(F),
    UnknownFile,
    PdbError,
    CifError,
    MmjsonError,
    ChemCompError,
}

// Ordered by FORMATS × INPUTS. Derived from mmread.hpp's independent memory
// and file switches and their PDB/CIF/mmJSON/chemcomp callees; see receipt.
const MEMORY: [[Outcome; 5]; 6] = [
    [
        Outcome::Success(F::Pdb, 0),
        Outcome::Success(F::Mmcif, 0),
        Outcome::Success(F::Mmcif, 0),
        Outcome::Success(F::ChemComp, 1),
        Outcome::Wrong(F::Unknown),
    ],
    [
        Outcome::Success(F::Pdb, 0),
        Outcome::Success(F::Mmcif, 0),
        Outcome::Success(F::Mmcif, 0),
        Outcome::Success(F::ChemComp, 1),
        Outcome::Wrong(F::Unknown),
    ],
    [
        Outcome::Success(F::Pdb, 0),
        Outcome::PdbError,
        Outcome::PdbError,
        Outcome::PdbError,
        Outcome::Success(F::Pdb, 0),
    ],
    [
        Outcome::CifError,
        Outcome::Success(F::Mmcif, 0),
        Outcome::CifError,
        Outcome::Success(F::ChemComp, 1),
        Outcome::CifError,
    ],
    [
        Outcome::MmjsonError,
        Outcome::MmjsonError,
        Outcome::Success(F::Mmcif, 0),
        Outcome::MmjsonError,
        Outcome::MmjsonError,
    ],
    [Outcome::Wrong(F::ChemComp); 5],
];
const FILE: [[Outcome; 5]; 6] = [
    [Outcome::UnknownFile; 5],
    MEMORY[1],
    MEMORY[2],
    [
        Outcome::CifError,
        Outcome::Success(F::Mmcif, 0),
        Outcome::CifError,
        Outcome::Success(F::Mmcif, 0),
        Outcome::CifError,
    ],
    [
        Outcome::MmjsonError,
        Outcome::MmjsonError,
        Outcome::Success(F::Mmjson, 0),
        Outcome::MmjsonError,
        Outcome::MmjsonError,
    ],
    [
        Outcome::CifError,
        Outcome::ChemCompError,
        Outcome::CifError,
        Outcome::Success(F::ChemComp, 1),
        Outcome::CifError,
    ],
];

fn assert_route(
    result: Result<cosmolkit_bio::BioStructureData, E>,
    expected: Outcome,
    source_name: &str,
    path: Option<&std::path::Path>,
) {
    match expected {
        Outcome::Success(format, atoms) => {
            let data =
                result.unwrap_or_else(|error| panic!("{expected:?}/{source_name}: {error:?}"));
            data.validate().unwrap();
            assert_eq!(data.input_format(), format, "{source_name}");
            assert_eq!(data.atoms().len(), atoms, "{source_name}");
            assert_eq!(data.coordinates().positions().len(), atoms, "{source_name}");
        }
        Outcome::Wrong(format) => assert!(
            matches!(result, Err(E::WrongFormat { source_name: actual, format: found }) if actual == source_name && found == format),
            "{expected:?}/{source_name}"
        ),
        Outcome::UnknownFile => assert!(
            matches!(result, Err(E::UnknownFileFormat(actual)) if Some(actual.as_path()) == path),
            "{expected:?}/{source_name}"
        ),
        Outcome::PdbError => assert!(
            matches!(result, Err(E::Pdb(ref error)) if error.source().is_some()),
            "{expected:?}/{source_name}"
        ),
        Outcome::CifError => assert!(
            matches!(result, Err(E::Cif(ref error)) if error.source() == source_name),
            "{expected:?}/{source_name}"
        ),
        Outcome::MmjsonError => assert!(
            matches!(result, Err(E::Mmjson(ref error)) if error.source().is_some()),
            "{expected:?}/{source_name}"
        ),
        Outcome::ChemCompError => assert!(
            matches!(result, Err(E::ChemComp(ref error)) if error.to_string() == "Not a chem_comp format."),
            "{expected:?}/{source_name}: {result:?}"
        ),
    }
}

#[test]
fn bio_read_detached_six_formats_cross_five_contents_memory_and_file() {
    // mmread.hpp: memory Unknown/Detect classify contents; file Unknown
    // classifies suffix; file Detect delegates to memory after loading.
    let dir = tempfile::tempdir().unwrap();
    for (content, (name, text, _)) in INPUTS.into_iter().enumerate() {
        let path = dir.path().join(format!("{name}.xyz"));
        std::fs::write(&path, text).unwrap();
        for (request, requested) in FORMATS.into_iter().enumerate() {
            let memory = read_bio_structure(
                text,
                &BioReadParams {
                    format: requested,
                    source_name: format!("{name}.memory"),
                },
            );
            assert_route(
                memory,
                MEMORY[request][content],
                &format!("{name}.memory"),
                None,
            );
            let file = read_bio_structure_file(&path, requested);
            assert_route(
                file,
                FILE[request][content],
                &path.to_string_lossy(),
                Some(&path),
            );
        }
    }
}

#[test]
fn bio_read_detached_extension_precedence_and_file_error_sources() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("disguised.cif");
    std::fs::write(&path, INPUTS[0].1).unwrap();
    assert!(matches!(
        read_bio_structure_file(&path, F::Unknown),
        Err(E::Cif(_))
    ));
    assert_eq!(
        read_bio_structure_file(&path, F::Detect)
            .unwrap()
            .input_format(),
        F::Pdb
    );
    assert_eq!(
        read_bio_structure_file(&path, F::Pdb)
            .unwrap()
            .input_format(),
        F::Pdb
    );

    let unknown = dir.path().join("missing.xyz");
    assert!(
        matches!(read_bio_structure_file(&unknown, F::Unknown), Err(E::UnknownFileFormat(ref path)) if path == &unknown)
    );
    let missing = dir.path().join("missing.cif");
    let error = read_bio_structure_file(&missing, F::Unknown).unwrap_err();
    assert!(matches!(error, E::Io { .. }));
    assert_eq!(
        error
            .source()
            .unwrap()
            .downcast_ref::<std::io::Error>()
            .unwrap()
            .kind(),
        std::io::ErrorKind::NotFound
    );
    let invalid = dir.path().join("invalid.json");
    std::fs::write(&invalid, b"data_demo\n\xff").unwrap();
    let error = read_bio_structure_file(&invalid, F::Unknown).unwrap_err();
    assert!(matches!(error, E::Utf8 { .. }));
    assert!(error.source().is_some());
}
