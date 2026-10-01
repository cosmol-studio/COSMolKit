#![cfg(feature = "cap-bio")]
use std::alloc::GlobalAlloc;

use cosmolkit::{
    BINDING_CONTRACT, BindingKind, BioStructure, FunctionStatus, Protein, ProteinAtomIter,
    ProteinChainIter, ProteinChainRef, ProteinResidueIter, ProteinResidueRef,
};

const CIF: &str = include_str!("../../../testdata/bio/fixtures/gemmi_full_feature_sample.cif");

thread_local! {
    static COUNTING: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static ALLOCATIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

struct CurrentThreadAllocator;

// SAFETY: All allocation/deallocation is delegated unchanged to System;
// the TLS probe only counts allocations on the calling thread and never allocates.
unsafe impl std::alloc::GlobalAlloc for CurrentThreadAllocator {
    unsafe fn alloc(&self, layout: std::alloc::Layout) -> *mut u8 {
        COUNTING.with(|enabled| {
            if enabled.get() {
                ALLOCATIONS.with(|count| count.set(count.get() + 1));
            }
        });
        // SAFETY: same layout and allocation contract as the caller.
        unsafe { std::alloc::System.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: std::alloc::Layout) -> *mut u8 {
        COUNTING.with(|enabled| {
            if enabled.get() {
                ALLOCATIONS.with(|count| count.set(count.get() + 1));
            }
        });
        // SAFETY: same layout and allocation contract as the caller.
        unsafe { std::alloc::System.alloc_zeroed(layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: std::alloc::Layout, size: usize) -> *mut u8 {
        COUNTING.with(|enabled| {
            if enabled.get() {
                ALLOCATIONS.with(|count| count.set(count.get() + 1));
            }
        });
        // SAFETY: same pointer, layout and size contract as the caller.
        unsafe { std::alloc::System.realloc(ptr, layout, size) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: std::alloc::Layout) {
        // SAFETY: same pointer and layout contract as the caller.
        unsafe { std::alloc::System.dealloc(ptr, layout) }
    }
}

#[global_allocator]
static ALLOCATOR: CurrentThreadAllocator = CurrentThreadAllocator;

fn measure_allocations<T>(work: impl FnOnce() -> T) -> (T, usize) {
    // Initialize both TLS slots before entering the measured region.
    COUNTING.with(|_| {});
    ALLOCATIONS.with(|count| count.set(0));
    COUNTING.with(|enabled| enabled.set(true));
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            COUNTING.with(|enabled| enabled.set(false));
        }
    }
    let reset = Reset;
    let result = work();
    drop(reset);
    (result, ALLOCATIONS.with(std::cell::Cell::get))
}

#[test]
fn bio_legacy_i06_six_traversals_create_step_and_exhaust_without_allocations() {
    let structure = BioStructure::from_mmcif(CIF).unwrap();
    let protein = structure.protein().unwrap();
    let chain = protein.chain(0).unwrap();
    let residue = chain.residues().next().unwrap();
    let empty = BioStructure::from_mmcif("data_empty\n")
        .unwrap()
        .protein()
        .unwrap();

    // Fixed source has one chain, two CYS residues and two atoms (IDs 0, 1).
    // Results are held in stack arrays; no Vec is constructed in a measured closure.
    macro_rules! check_ids {
        ($iterator:expr, $expected:expr) => {{
            let ((ids, count), allocations) = measure_allocations(|| {
                let mut iter = $iterator;
                let first = iter.next().map(|item| item.id().value());
                let mut ids = [u32::MAX; 8];
                let mut count = 0;
                if let Some(id) = first {
                    ids[count] = id;
                    count += 1;
                }
                for item in iter.by_ref() {
                    ids[count] = item.id().value();
                    count += 1;
                }
                assert!(iter.next().is_none());
                (ids, count)
            });
            assert_eq!(allocations, 0);
            let expected: &[u32] = $expected;
            assert_eq!(&ids[..count], expected);
        }};
    }
    check_ids!(protein.chains(), &[0]);
    check_ids!(protein.residues(), &[0, 1]);
    check_ids!(protein.atoms(), &[0, 1]);
    check_ids!(chain.residues(), &[0, 1]);
    check_ids!(chain.atoms(), &[0, 1]);
    check_ids!(residue.atoms(), &[0]);
    check_ids!(empty.chains(), &[]);
    check_ids!(empty.residues(), &[]);
    check_ids!(empty.atoms(), &[]);
    let (early, count) = measure_allocations(|| {
        let mut iter = protein.atoms();
        iter.next().unwrap().id().value()
    });
    assert_eq!((early, count), (0, 0));
}

#[test]
fn bio_legacy_i05_canonical_borrowed_iterators_and_registry() {
    let _: for<'a> fn(&'a Protein) -> ProteinChainIter<'a> = Protein::chains;
    let _: for<'a> fn(&'a Protein) -> ProteinResidueIter<'a> = Protein::residues;
    let _: for<'a> fn(&'a Protein) -> ProteinAtomIter<'a> = Protein::atoms;
    fn borrowed<'a>(chain: ProteinChainRef<'a>, residue: ProteinResidueRef<'a>) {
        let _: ProteinResidueIter<'a> = chain.residues();
        let _: ProteinAtomIter<'a> = chain.atoms();
        let _: ProteinAtomIter<'a> = residue.atoms();
    }
    let _: for<'a> fn(ProteinChainRef<'a>, ProteinResidueRef<'a>) = borrowed;
    for (id, rust_type, python, javascript) in [
        (
            "types.ProteinChainIter",
            "ProteinChainIter",
            "ProteinChainIter",
            "ProteinChainIter",
        ),
        (
            "types.ProteinResidueIter",
            "ProteinResidueIter",
            "ProteinResidueIter",
            "ProteinResidueIter",
        ),
        (
            "types.ProteinAtomIter",
            "ProteinAtomIter",
            "ProteinAtomIter",
            "ProteinAtomIter",
        ),
    ] {
        let matches = BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id == id)
            .collect::<Vec<_>>();
        assert_eq!(matches.len(), 1, "{id}");
        assert!(matches[0].rust_path.contains(rust_type));
        assert_eq!(matches[0].python_name, python);
        assert_eq!(matches[0].javascript_name, javascript);
        assert_eq!(matches[0].feature, "cap-bio");
        assert_eq!(matches[0].status, FunctionStatus::Experimental);
    }
    for id in [
        "Protein.chains",
        "Protein.residues",
        "Protein.atoms",
        "ProteinChainRef.residues",
        "ProteinChainRef.atoms",
        "ProteinResidueRef.atoms",
    ] {
        let matches = BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id == id)
            .collect::<Vec<_>>();
        assert_eq!(matches.len(), 1, "{id}");
        assert_eq!(matches[0].feature, "cap-bio");
        assert_eq!(matches[0].status, FunctionStatus::Experimental);
    }
    let structure = BioStructure::from_mmcif(CIF).unwrap();
    let protein = structure.protein().unwrap();
    let mut chains = protein.chains();
    let first = chains.next().unwrap();
    assert_eq!(first.id().value(), 0);
    let first_residue = first.residues().next().unwrap();
    assert_eq!(first_residue.id(), protein.residues().next().unwrap().id());
    assert_eq!(
        first_residue.atoms().next().unwrap().id(),
        protein.atoms().next().unwrap().id()
    );
    assert_eq!(
        first.atoms().next().unwrap().id(),
        protein.atoms().next().unwrap().id()
    );
    assert_eq!(protein.chains().count(), protein.num_chains());
    assert_eq!(protein.residues().count(), protein.num_residues());
    assert_eq!(protein.atoms().count(), protein.num_atoms());
    assert_eq!(
        first_residue.atoms().next().unwrap().position(),
        protein.atoms().next().unwrap().position()
    );
}

#[test]
fn bio_pdbscope_public_text_signatures_defaults_registry_and_bytes() {
    // Compiler-checked method signatures through the public facade.
    let _: for<'a> fn(&'a BioStructure) -> Result<String, cosmolkit::BioMmcifWriteError> =
        BioStructure::to_mmcif;
    let _: for<'a> fn(
        &'a BioStructure,
        &cosmolkit::BioMmcifWriteParams,
    ) -> Result<String, cosmolkit::BioMmcifWriteError> = BioStructure::to_mmcif_with_params;
    // Pinned parameter defaults: the exact seven coordinate-profile controls.
    let defaults = cosmolkit::BioMmcifWriteParams::default();
    assert!(defaults.group_pdb);
    assert!(!defaults.auth_all);
    assert!(!defaults.prefer_pairs);
    assert!(!defaults.compact);
    assert!(!defaults.misuse_hash);
    assert_eq!((defaults.align_pairs, defaults.align_loops), (0, 0));
    // Registry: one row per callable and per type, Experimental and bio-gated.
    for (id, python, javascript) in [
        ("BioStructure.to_mmcif", "to_mmcif", "toMmcif"),
        (
            "BioStructure.to_mmcif_with_params",
            "to_mmcif_with_params",
            "toMmcifWithParams",
        ),
        (
            "types.BioMmcifWriteParams",
            "BioMmcifWriteParams",
            "BioMmcifWriteParams",
        ),
        (
            "types.BioMmcifWriteError",
            "BioMmcifWriteError",
            "BioMmcifWriteError",
        ),
    ] {
        let matches = BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id == id)
            .collect::<Vec<_>>();
        assert_eq!(matches.len(), 1, "{id}");
        assert_eq!(matches[0].feature, "cap-bio");
        assert_eq!(matches[0].status, FunctionStatus::Experimental);
        assert_eq!(matches[0].python_name, python);
        assert_eq!(matches[0].javascript_name, javascript);
    }
    let callable = BINDING_CONTRACT
        .iter()
        .find(|row| row.semantic_id == "BioStructure.to_mmcif")
        .unwrap()
        .callable
        .unwrap();
    assert_eq!(callable.kind, BindingKind::Instance);
    assert!(callable.parameters.is_empty());
    let with_params = BINDING_CONTRACT
        .iter()
        .find(|row| row.semantic_id == "BioStructure.to_mmcif_with_params")
        .unwrap()
        .callable
        .unwrap();
    assert_eq!(with_params.kind, BindingKind::Instance);
    assert_eq!(with_params.parameters.len(), 1);

    // Pinned-profile bytes through the public method: the real fixture's
    // first atom row and preamble, derived from the source layout rules.
    let structure = BioStructure::from_mmcif(CIF).unwrap();
    let text = structure.to_mmcif().unwrap();
    let first_line = text.lines().next().unwrap();
    assert!(first_line.starts_with("data_"), "header line: {first_line}");
    assert!(text.contains("\nloop_\n_atom_site.group_PDB\n_atom_site.id\n"));
    assert!(text.contains("_entry.id "));
    // No non-coordinate category is emitted even for the full-feature input.
    for forbidden in [
        "_cell.",
        "_symmetry.",
        "_struct_ncs_oper.",
        "_struct_conn.",
        "_atom_type.",
        "_pdbx_struct_assembly.",
    ] {
        assert!(!text.contains(forbidden), "{forbidden} must be absent");
    }
    // The borrowed structure is unchanged and both entries agree on bytes.
    assert_eq!(
        structure
            .to_mmcif_with_params(&cosmolkit::BioMmcifWriteParams::default())
            .unwrap(),
        text
    );
    assert_eq!(
        structure.num_atoms(),
        BioStructure::from_mmcif(CIF).unwrap().num_atoms()
    );
    // The zero-model early return still yields the bare header document.
    let empty = BioStructure::from_mmcif("data_demo\n_entry.id DEMO\n").unwrap();
    assert_eq!(empty.to_mmcif().unwrap(), "data_\n");
}

#[test]
fn bio_pdbscope_public_file_output_errors_input_and_protein_projection() {
    // Compiler-checked signatures through the public facade.
    let _: for<'a> fn(
        &'a BioStructure,
        &std::path::Path,
    ) -> Result<(), cosmolkit::BioMmcifWriteError> = BioStructure::write_mmcif;
    let _: for<'a> fn(
        &'a BioStructure,
        &std::path::Path,
        &cosmolkit::BioMmcifWriteParams,
    ) -> Result<(), cosmolkit::BioMmcifWriteError> = BioStructure::write_mmcif_with_params;
    for id in [
        "BioStructure.write_mmcif",
        "BioStructure.write_mmcif_with_params",
    ] {
        let matches = BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id == id)
            .collect::<Vec<_>>();
        assert_eq!(matches.len(), 1, "{id}");
        assert_eq!(matches[0].feature, "cap-bio");
        assert_eq!(matches[0].status, FunctionStatus::Experimental);
        assert_eq!(matches[0].callable.unwrap().kind, BindingKind::Instance);
    }

    let structure = BioStructure::from_mmcif(CIF).unwrap();
    let text_before = structure.to_mmcif().unwrap();
    let directory = std::env::temp_dir();
    let path = directory.join(format!("ck-bio-pdb-public-file-{}.cif", std::process::id()));
    // Default and explicit-default output are byte-identical.
    structure.write_mmcif(&path).unwrap();
    let default_bytes = std::fs::read(&path).unwrap();
    structure
        .write_mmcif_with_params(&path, &cosmolkit::BioMmcifWriteParams::default())
        .unwrap();
    assert_eq!(std::fs::read(&path).unwrap(), default_bytes);
    let text = String::from_utf8(default_bytes).unwrap();
    assert!(text.starts_with("data_"));
    assert!(text.contains("_atom_site."));
    // Coordinate-only: no lossless whole-structure claim; forbidden
    // categories stay absent on the full-feature fixture.
    for forbidden in [
        "_cell.",
        "_struct_ncs_oper.",
        "_struct_conn.",
        "_atom_type.",
    ] {
        assert!(!text.contains(forbidden), "{forbidden} must be absent");
    }
    // The input object is unchanged: the same bytes serialize after every
    // write (whole-object equality is unusable on this fixture because it
    // legitimately carries NaN in its assembly metadata).
    assert_eq!(structure.to_mmcif().unwrap(), text_before);
    // Path/source errors retain the destination path.
    let missing = std::path::PathBuf::from("/nonexistent-ck-bio-public-dir/out.cif");
    let error = structure.write_mmcif(&missing).unwrap_err();
    match &error {
        cosmolkit::BioMmcifWriteError::FileWrite { path, source } => {
            assert_eq!(path, &missing);
            assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
        }
        other => panic!("expected FileWrite, got {other:?}"),
    }
    // Protein projection routes through as_bio_structure; no alias method.
    let protein = structure.protein().unwrap();
    let protein_path = directory.join(format!("ck-bio-pdb-public-prot-{}.cif", std::process::id()));
    protein
        .as_bio_structure()
        .write_mmcif(&protein_path)
        .unwrap();
    let protein_text = String::from_utf8(std::fs::read(&protein_path).unwrap()).unwrap();
    assert!(protein_text.contains("_atom_site."));
    assert!(protein.num_atoms() <= structure.num_atoms());
    let _ = std::fs::remove_file(&path);
    let _ = std::fs::remove_file(&protein_path);
}

#[test]
fn bio_pdbscope_matrix_empty_structure_all_288_combinations() {
    // Zero-model early return: only the misuse_hash fences can appear; the
    // other seven controls cannot change the bytes.
    let empty = BioStructure::from_mmcif("data_demo\n_entry.id DEMO\n").unwrap();
    for params in matrix_params() {
        let expected = if params.misuse_hash {
            "data_\n#\n#\n"
        } else {
            "data_\n"
        };
        assert_eq!(
            empty.to_mmcif_with_params(&params).unwrap(),
            expected,
            "params {params:?}"
        );
    }
}

#[test]
fn bio_pdbscope_matrix_single_model_all_288_combinations() {
    let structure = BioStructure::from_parts(matrix_f2_parts()).unwrap();
    // Hand-derived source-profile values (see receipt step 122 spec).
    let rows = vec![
        vec![
            "1", "C", "CA", "A", "ALA", "A", ".", "3", "B", "1", "2.5", "3", "1", "20", "-1", "7",
            "A", "1",
        ],
        vec![
            "2", "N", "N", ".", "ALA", "A", ".", "3", "B", "-4", "5", "6", "0.5", "10.5", "?", "7",
            "A", "1",
        ],
    ];
    let aniso = vec![vec![
        "2", "N", "0.05", "0.01", "0.02", "0.03", "0.04", "0.06",
    ]];
    for params in matrix_params() {
        let expected = matrix_expected(
            &params,
            "1abc",
            Some(("_entry.id", "1abc")),
            &rows,
            &["ATOM", "ATOM"],
            &aniso,
            false,
        );
        assert_eq!(
            structure.to_mmcif_with_params(&params).unwrap(),
            expected,
            "params {params:?}"
        );
    }
}

#[test]
fn bio_pdbscope_matrix_multimodel_all_288_combinations() {
    let structure = BioStructure::from_parts(matrix_f3_parts()).unwrap();
    let rows = vec![
        vec![
            "1", "P", "P", ".", "DA", "P1", ".", "10", "?", "7", "8", "9", "1", "15", "?", "10",
            "A", "1", "c", "2", "0.25",
        ],
        vec![
            "2", "O", "O", ".", "HOH", "W1", ".", ".", "?", "1.5", "2.5", "3.5", "1", "30", "?",
            "?", "B", "3", ".", "?", "0",
        ],
        vec![
            "3", "C", "C1", ".", "LIG", "L1", "ligent", "5", "?", "2", "1", "4", "1", "12.5", "-2",
            "5", "C", "3", ".", "?", "0",
        ],
    ];
    let aniso = vec![vec![
        "3", "C", "0.1", "0.2", "0.3", "0.04", "0.05", "0.06", "3",
    ]];
    for params in matrix_params() {
        let expected = matrix_expected(
            &params,
            "multi",
            Some(("_entry.id", "M1")),
            &rows,
            &["ATOM", "HETATM", "HETATM"],
            &aniso,
            true,
        );
        assert_eq!(
            structure.to_mmcif_with_params(&params).unwrap(),
            expected,
            "params {params:?}"
        );
    }
}

#[test]
fn bio_pdbscope_matrix_parsed_values_and_order_roundtrip() {
    // Independent parsed-value oracle: the emitted default-profile text is
    // re-read through the delivered mmCIF reader; every value below is
    // pinned from the fixture definition, not from writer output.
    let single = BioStructure::from_parts(matrix_f2_parts()).unwrap();
    let reread = BioStructure::from_mmcif(&single.to_mmcif().unwrap()).unwrap();
    assert_eq!(reread.atoms().len(), 2);
    assert_eq!(reread.residues().len(), 1);
    assert_eq!(reread.models().len(), 1);
    assert_eq!(reread.input_format(), cosmolkit::BioCoordinateFormat::Mmcif);
    let names: Vec<String> = reread
        .atoms()
        .iter()
        .map(|a| a.name().as_str().to_string())
        .collect();
    assert_eq!(names, ["CA", "N"]);
    let serials: Vec<i32> = reread
        .atoms()
        .iter()
        .map(|a| a.source().serial().map(|s| s.value()).unwrap_or(-1))
        .collect();
    assert_eq!(serials, [1, 2]);
    assert_eq!(reread.residues()[0].name().as_str(), "ALA");
    assert_eq!(reread.residues()[0].source().label_seq_id(), Some(3));
    assert_eq!(
        reread.residues()[0]
            .source()
            .seq_id()
            .map(|s| (s.seq_num(), s.ins_code())),
        Some((7, Some(b'B')))
    );
    assert_eq!(reread.atoms()[0].altloc().map(|l| l.value()), Some(b'A'));
    assert_eq!(reread.atoms()[0].formal_charge(), -1);
    assert_eq!(reread.atoms()[1].b_iso(), 10.5);
    // The reader narrows anisotropic components to f32 (source SMat33<float>),
    // so the exact roundtripped values are the f32-narrowed constants.
    assert_eq!(
        reread.atoms()[1].anisou(),
        &[
            0.05f32 as f64,
            0.01f32 as f64,
            0.02f32 as f64,
            0.03f32 as f64,
            0.04f32 as f64,
            0.06f32 as f64,
        ]
    );
    assert_eq!(reread.coordinates().positions()[0], [1.0, 2.5, 3.0]);
    assert_eq!(reread.coordinates().positions()[1], [-4.0, 5.0, 6.0]);

    let multi = BioStructure::from_parts(matrix_f3_parts()).unwrap();
    let reread = BioStructure::from_mmcif(&multi.to_mmcif().unwrap()).unwrap();
    assert_eq!(reread.models().len(), 2);
    let model_nums: Vec<Option<i32>> = reread
        .models()
        .iter()
        .map(|m| m.source_model_number())
        .collect();
    assert_eq!(model_nums, [Some(1), Some(3)]);
    let names: Vec<String> = reread
        .atoms()
        .iter()
        .map(|a| a.name().as_str().to_string())
        .collect();
    assert_eq!(names, ["P", "O", "C1"]);
    let residues: Vec<String> = reread
        .residues()
        .iter()
        .map(|r| r.name().as_str().to_string())
        .collect();
    assert_eq!(residues, ["DA", "HOH", "LIG"]);
    assert_eq!(
        reread.atoms()[0].calc_flag(),
        cosmolkit::BioCalcFlag::Calculated
    );
    assert_eq!(reread.atoms()[0].tls_group_id(), 2);
    assert_eq!(reread.atoms()[0].fraction(), 0.25);
    assert_eq!(reread.atoms()[2].formal_charge(), -2);
    assert_eq!(
        reread.atoms()[2].anisou(),
        &[
            0.1f32 as f64,
            0.2f32 as f64,
            0.3f32 as f64,
            0.04f32 as f64,
            0.05f32 as f64,
            0.06f32 as f64,
        ]
    );
    // Layout variants must not change parsed values.
    for params in [
        cosmolkit::BioMmcifWriteParams {
            prefer_pairs: true,
            ..cosmolkit::BioMmcifWriteParams::default()
        },
        cosmolkit::BioMmcifWriteParams {
            align_pairs: 33,
            align_loops: 33,
            compact: true,
            ..cosmolkit::BioMmcifWriteParams::default()
        },
    ] {
        let again =
            BioStructure::from_mmcif(&single.to_mmcif_with_params(&params).unwrap()).unwrap();
        assert_eq!(
            again.atoms()[1].anisou(),
            &[
                0.05f32 as f64,
                0.01f32 as f64,
                0.02f32 as f64,
                0.03f32 as f64,
                0.04f32 as f64,
                0.06f32 as f64,
            ]
        );
        assert_eq!(again.atoms().len(), 2);
    }
}

fn matrix_params() -> impl Iterator<Item = cosmolkit::BioMmcifWriteParams> {
    (0..32u32).flat_map(|bits| {
        let bools = [
            bits & 1 != 0,
            bits & 2 != 0,
            bits & 4 != 0,
            bits & 8 != 0,
            bits & 16 != 0,
        ];
        [0u16, 1, 33].into_iter().flat_map(move |align_pairs| {
            [0u16, 1, 33]
                .into_iter()
                .map(move |align_loops| cosmolkit::BioMmcifWriteParams {
                    group_pdb: bools[0],
                    auth_all: bools[1],
                    prefer_pairs: bools[2],
                    compact: bools[3],
                    misuse_hash: bools[4],
                    align_pairs,
                    align_loops,
                })
        })
    })
}

/// Assemble the expected document from the pinned source layout rules
/// (to_cif.hpp write_out_pair / write_out_loop / write_cif_block_to_stream
/// / should_be_separated_) over hand-derived fixture values — never over
/// writer output.
#[allow(clippy::too_many_arguments)]
fn matrix_expected(
    params: &cosmolkit::BioMmcifWriteParams,
    block: &str,
    entry: Option<(&str, &str)>,
    rows: &[Vec<&str>],
    groups: &[&str],
    aniso: &[Vec<&str>],
    aniso_model_column: bool,
) -> String {
    // Atom-loop tags per the source column rules.
    let mut tags: Vec<String> = [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "pdbx_formal_charge",
        "auth_atom_id",
        "auth_comp_id",
        "auth_seq_id",
        "auth_asym_id",
        "pdbx_PDB_model_num",
    ]
    .iter()
    .map(|s| format!("_atom_site.{s}"))
    .collect();
    let optional = rows.first().map_or(0, |r| r.len().saturating_sub(18));
    for suffix in ["calc_flag", "pdbx_tls_group_id", "ccp4_deuterium_fraction"]
        .iter()
        .take(optional)
    {
        tags.push(format!("_atom_site.{suffix}"));
    }
    if !params.auth_all {
        let tail = tags.split_off(tags.len() - 5 - optional);
        let keep_auth = &tail[2..];
        tags.extend_from_slice(keep_auth);
    }
    if !params.group_pdb {
        tags.remove(0);
    }
    let mut aniso_tags: Vec<String> = [
        "id",
        "type_symbol",
        "U[1][1]",
        "U[2][2]",
        "U[3][3]",
        "U[1][2]",
        "U[1][3]",
        "U[2][3]",
    ]
    .iter()
    .map(|s| format!("_atom_site_anisotrop.{s}"))
    .collect();
    if aniso_model_column {
        aniso_tags.push("_atom_site_anisotrop.pdbx_PDB_model_num".to_string());
    }

    let pair_line = |tag: &str, value: &str| -> String {
        let mut line = format!("{tag} ");
        if tag.len() < params.align_pairs as usize {
            line.push_str(&" ".repeat(params.align_pairs as usize - tag.len()));
        }
        line.push_str(value);
        line.push('\n');
        line
    };
    let loop_block = |tags: &[String], rows: &[Vec<&str>]| -> String {
        if params.prefer_pairs && rows.len() == 1 {
            let mut out = String::new();
            for (tag, value) in tags.iter().zip(&rows[0]) {
                out.push_str(&pair_line(tag, value));
            }
            return out;
        }
        let mut out = String::from("loop_");
        for tag in tags {
            out.push('\n');
            out.push_str(tag);
        }
        let widths: Vec<usize> = if params.align_loops > 0 {
            (0..tags.len())
                .map(|col| {
                    rows.iter()
                        .map(|r| r[col].len())
                        .max()
                        .unwrap_or(0)
                        .min(params.align_loops as usize)
                })
                .collect()
        } else {
            vec![0; tags.len()]
        };
        for row in rows {
            out.push('\n');
            for col in 0..tags.len() {
                out.push_str(row[col]);
                if col != tags.len() - 1 {
                    if row[col].len() < widths[col] {
                        out.push_str(&" ".repeat(widths[col] - row[col].len()));
                    }
                    out.push(' ');
                }
            }
        }
        out.push('\n');
        out
    };
    fn category(tag: &str) -> &str {
        tag.split('.').next().unwrap_or("")
    }

    let mut items: Vec<(String, String)> = Vec::new();
    if let Some((tag, value)) = entry {
        items.push((tag.to_string(), value.to_string()));
    }
    // Effective value rows follow the source column order: optional
    // group_PDB first, the two auth duplicates after the charge column.
    let effective_rows: Vec<Vec<&str>> = rows
        .iter()
        .zip(groups)
        .map(|(row, group)| {
            let mut values: Vec<&str> = Vec::with_capacity(row.len() + 3);
            if params.group_pdb {
                values.push(group);
            }
            values.extend_from_slice(row);
            if params.auth_all {
                let atom_id = row[2];
                let comp_id = row[4];
                values.insert(values.len() - 3 - optional, atom_id);
                values.insert(values.len() - 3 - optional, comp_id);
            }
            values
        })
        .collect();
    let atom_text = loop_block(&tags, &effective_rows);
    let aniso_text = loop_block(&aniso_tags, aniso);

    let mut out = format!("data_{block}\n");
    if params.misuse_hash {
        out.push_str("#\n");
    }
    let mut prev_pair_category: Option<String> = None;
    for (tag, value) in &items {
        if !params.compact
            && prev_pair_category
                .as_deref()
                .is_some_and(|c| c != category(tag))
        {
            out.push_str(if params.misuse_hash { "#\n" } else { "\n" });
        }
        out.push_str(&pair_line(tag, value));
        prev_pair_category = Some(category(tag).to_string());
    }
    // The atom loop is a Loop item: always separated from any previous item.
    if !atom_text.is_empty() {
        if !params.compact && prev_pair_category.is_some() {
            out.push_str(if params.misuse_hash { "#\n" } else { "\n" });
        }
        out.push_str(&atom_text);
        prev_pair_category = None;
    }
    if !aniso_text.is_empty() {
        if !params.compact {
            out.push_str(if params.misuse_hash { "#\n" } else { "\n" });
        }
        out.push_str(&aniso_text);
    }
    if params.misuse_hash {
        out.push_str("#\n");
    }
    out
}

fn matrix_f2_parts() -> cosmolkit::BioStructureParts {
    use cosmolkit::{
        AtomName, AtomSourceIds, BioAtomRow, BioCalcFlag, BioChainId, BioChainRow,
        BioCoordinateBlock, BioCoordinateFormat, BioModelId, BioModelRow, BioResidueId,
        BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructureParts, BioStructureSourceState,
        ChainKind, ChainSourceIds, EntityKind, PdbChainId, PdbSeqId, ResidueInfoKind, ResidueName,
        ResidueSourceIds,
    };
    let atom = |residue: u32,
                name: &[u8],
                element: cosmolkit::Element,
                altloc: Option<u8>,
                charge: i8,
                occ: f64,
                b_iso: f64,
                anisou: [f64; 6]|
     -> BioAtomRow {
        BioAtomRow::new(
            BioResidueId::new(residue),
            AtomName::from_ascii(name).unwrap(),
            element,
            None,
            altloc.map(cosmolkit::AltLocLabel::new),
            charge,
            BioCalcFlag::NotSet,
            occ,
            b_iso,
            anisou,
            -1,
            0.0,
            AtomSourceIds::default(),
        )
    };
    BioStructureParts {
        input_format: BioCoordinateFormat::Mmcif,
        models: vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))],
        chains: vec![BioChainRow::new(
            BioModelId::new(0),
            None,
            BioRowSpan::new(0, 1).unwrap(),
            ChainKind::Protein,
            ChainSourceIds::new(
                Some(PdbChainId::from_ascii(b"A").unwrap()),
                Some("A".to_string()),
            ),
        )],
        residues: vec![BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(0, 2).unwrap(),
            ResidueName::from_ascii(b"ALA").unwrap(),
            ResidueInfoKind::Aa,
            EntityKind::Polymer,
            None,
            None,
            ResidueSourceIds::new(
                Some(PdbSeqId::new(7, Some(b'B'))),
                Some(3),
                None,
                Some("A".to_string()),
                None,
            )
            .unwrap(),
            BioSiftsUnpResidue::default(),
        )],
        atoms: vec![
            atom(
                0,
                b"CA",
                cosmolkit::Element::C,
                Some(b'A'),
                -1,
                1.0,
                20.0,
                [0.0; 6],
            ),
            atom(
                0,
                b"N",
                cosmolkit::Element::N,
                None,
                0,
                0.5,
                10.5,
                [0.05, 0.01, 0.02, 0.03, 0.04, 0.06],
            ),
        ],
        entities: Vec::new(),
        connections: Vec::new(),
        cispeps: Vec::new(),
        mod_residues: Vec::new(),
        helices: Vec::new(),
        sheets: Vec::new(),
        metadata: cosmolkit::BioMetadata::default(),
        source_state: BioStructureSourceState {
            name: "1abc".to_string(),
            ..BioStructureSourceState::default()
        },
        coordinates: BioCoordinateBlock::new(vec![[1.0, 2.5, 3.0], [-4.0, 5.0, 6.0]]),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    }
}

fn matrix_f3_parts() -> cosmolkit::BioStructureParts {
    use cosmolkit::Element;
    use cosmolkit::{
        AtomName, AtomSourceIds, BioAtomRow, BioCalcFlag, BioChainId, BioChainRow,
        BioCoordinateBlock, BioCoordinateFormat, BioEntityRow, BioModelId, BioModelRow,
        BioResidueId, BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructureParts,
        BioStructureSourceState, ChainKind, ChainSourceIds, EntityKind, EntitySourceIds,
        PdbChainId, PdbSeqId, PolymerKind, ResidueInfoKind, ResidueName, ResidueSourceIds,
    };
    use std::collections::BTreeMap;
    let atom = |residue: u32,
                name: &[u8],
                element: Element,
                charge: i8,
                calc: BioCalcFlag,
                b_iso: f64,
                anisou: [f64; 6],
                tls: i16,
                fraction: f64|
     -> BioAtomRow {
        BioAtomRow::new(
            BioResidueId::new(residue),
            AtomName::from_ascii(name).unwrap(),
            element,
            None,
            None,
            charge,
            calc,
            1.0,
            b_iso,
            anisou,
            tls,
            fraction,
            AtomSourceIds::default(),
        )
    };
    let residue = |chain: u32,
                   span: (u32, u32),
                   name: &[u8],
                   info: ResidueInfoKind,
                   kind: EntityKind,
                   het: Option<u8>,
                   seq: (i32, Option<u8>),
                   label_seq: Option<i32>,
                   subchain: &str|
     -> BioResidueRow {
        BioResidueRow::new(
            BioChainId::new(chain),
            BioRowSpan::new(span.0, span.1).unwrap(),
            ResidueName::from_ascii(name).unwrap(),
            info,
            kind,
            None,
            het,
            ResidueSourceIds::new(
                Some(PdbSeqId::new(seq.0, seq.1)),
                label_seq,
                None,
                Some(subchain.to_string()),
                None,
            )
            .unwrap(),
            BioSiftsUnpResidue::default(),
        )
    };
    let chain = |model: u32, span: (u32, u32), auth: u8, subchain: &str| -> BioChainRow {
        BioChainRow::new(
            BioModelId::new(model),
            None,
            BioRowSpan::new(span.0, span.1).unwrap(),
            ChainKind::Unknown,
            ChainSourceIds::new(
                Some(PdbChainId::from_ascii(&[auth]).unwrap()),
                Some(subchain.to_string()),
            ),
        )
    };
    BioStructureParts {
        input_format: BioCoordinateFormat::Mmcif,
        models: vec![
            BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1)),
            BioModelRow::new(BioRowSpan::new(1, 2).unwrap(), Some(3)),
        ],
        chains: vec![
            chain(0, (0, 1), b'A', "P1"),
            chain(1, (1, 1), b'B', "W1"),
            chain(1, (2, 1), b'C', "L1"),
        ],
        residues: vec![
            residue(
                0,
                (0, 1),
                b"DA",
                ResidueInfoKind::Dna,
                EntityKind::Polymer,
                None,
                (10, None),
                Some(10),
                "P1",
            ),
            residue(
                1,
                (1, 1),
                b"HOH",
                ResidueInfoKind::Hoh,
                EntityKind::Water,
                Some(b'H'),
                (i32::MIN, None),
                None,
                "W1",
            ),
            residue(
                2,
                (2, 1),
                b"LIG",
                ResidueInfoKind::Unknown,
                EntityKind::NonPolymer,
                None,
                (5, None),
                Some(5),
                "L1",
            ),
        ],
        atoms: vec![
            atom(
                0,
                b"P",
                Element::P,
                0,
                BioCalcFlag::Calculated,
                15.0,
                [0.0; 6],
                2,
                0.25,
            ),
            atom(
                1,
                b"O",
                Element::O,
                0,
                BioCalcFlag::NotSet,
                30.0,
                [0.0; 6],
                -1,
                0.0,
            ),
            atom(
                2,
                b"C1",
                Element::C,
                -2,
                BioCalcFlag::NotSet,
                12.5,
                [0.1, 0.2, 0.3, 0.04, 0.05, 0.06],
                -1,
                0.0,
            ),
        ],
        entities: vec![BioEntityRow::new(
            EntityKind::NonPolymer,
            PolymerKind::Unknown,
            false,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            vec!["L1".to_string()],
            EntitySourceIds::new("ligent".to_string()),
        )],
        connections: Vec::new(),
        cispeps: Vec::new(),
        mod_residues: Vec::new(),
        helices: Vec::new(),
        sheets: Vec::new(),
        metadata: cosmolkit::BioMetadata::default(),
        source_state: BioStructureSourceState {
            name: "multi".to_string(),
            info: BTreeMap::from([("_entry.id".to_string(), "M1".to_string())]),
            has_d_fraction: true,
            ..BioStructureSourceState::default()
        },
        coordinates: BioCoordinateBlock::new(vec![
            [7.0, 8.0, 9.0],
            [1.5, 2.5, 3.5],
            [2.0, 1.0, 4.0],
        ]),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    }
}

#[test]
fn bio_legacy_n03_public_borrow_registry_and_protein_projection() {
    let _: for<'a> fn(&'a BioStructure) -> Option<&'a str> = BioStructure::ncs_oper_identity_id;
    let row = BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == "BioStructure.ncs_oper_identity_id")
        .unwrap();
    assert_eq!(row.feature, "cap-bio");
    assert_eq!(row.status, FunctionStatus::Experimental);
    assert_eq!(row.python_name, "ncs_oper_identity_id");
    assert_eq!(row.javascript_name, "ncsOperIdentityId");
    assert_eq!(row.callable.unwrap().kind, BindingKind::Instance);
    assert!(row.callable.unwrap().parameters.is_empty());
    let absent = BioStructure::from_mmcif("data_demo\n_entry.id DEMO\n").unwrap();
    assert_eq!(absent.ncs_oper_identity_id(), None);
    let structure = BioStructure::from_mmcif(CIF).unwrap();
    assert_eq!(structure.ncs_oper_identity_id(), Some("I"));
    assert_eq!(
        structure.ncs_oper_identity_id().unwrap().as_ptr(),
        structure.source_state().info["_struct_ncs_oper.id"].as_ptr()
    );
    let protein = structure.protein().unwrap();
    assert_eq!(protein.as_bio_structure().ncs_oper_identity_id(), Some("I"));
}
