//! Reusable Rust source-port boundary for the IUPAC InChI engine.
//!
//! The crate exposes the four audited toolkit-neutral scalar InChI operations.
//! It is independent of `cosmolkit-core`, RDKit, Python, and external
//! executables so completed functionality can be reused by other Rust projects.
//!
//! See `dev/archive/plans/rdkit_inchi_full_port_plan.md` in the workspace for
//! the required source inventory, marker, and parity contract.

// The generated ABI constants intentionally retain the official C literals.
// Replacing them with Rust standard-library constants would weaken the
// source-reproduction audit without changing their compiled values.
#[allow(clippy::approx_constant)]
mod source_types;

// These control-flow shapes and assignments are direct reproductions of the
// official C source. Newer Clippy releases diagnose them as Rust idioms, but
// rewriting them would obscure the line-level source correspondence.
#[allow(
    clippy::almost_swapped,
    clippy::approx_constant,
    clippy::never_loop,
    clippy::while_immutable_condition
)]
mod source;

pub use source::rdkit::inchi::{
    AdapterAtom as InchiAtom, AdapterBond as InchiBond, AdapterDiagnostic as InchiDiagnostic,
    AdapterDiagnosticLevel as InchiDiagnosticLevel, AdapterGraphError as InchiGraphError, AdapterMol as InchiMolecule,
    AdapterToolkitError as InchiToolkitError, BondDirection as InchiBondDirection, BondStereo as InchiBondStereo,
    BondType as InchiBondType, ChiralTag as InchiChiralTag, InchiToMolToolkit, MolToInchiToolkit,
};

use source::rdkit::inchi::{
    InchiToInchiKeyError, InchiToMolError, MolToInchiError, MolToInchiKeyError, SourceInchiGenerationEngine,
    SourceInchiKeyEngine, SourceInchiStructureEngine,
};
use source_types::{SourceHeap, SourceHeapError};

/// Stable category for failures at the toolkit-neutral InChI boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InchiErrorKind {
    AllocationFailed,
    UnsupportedState,
    InvalidInput,
    InvalidSourceOutput,
    SanitizeFailed,
    Toolkit,
    SourcePort,
}

/// Structured failure from one of the four scalar InChI operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InchiError {
    pub operation: &'static str,
    pub kind: InchiErrorKind,
    pub detail: String,
}

impl std::fmt::Display for InchiError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "{} failed ({:?}): {}",
            self.operation, self.kind, self.detail
        )
    }
}

impl std::error::Error for InchiError {}

/// Return fields written by the official InChI generation or parsing API.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct InchiReturnValues {
    pub return_code: i32,
    pub message: Vec<u8>,
    pub log: Vec<u8>,
    pub aux_info: Vec<u8>,
}

/// Output of [`mol_to_inchi`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MolToInchiOutput {
    pub inchi: Vec<u8>,
    pub return_values: InchiReturnValues,
    pub diagnostics: Vec<InchiDiagnostic>,
}

/// Output of [`mol_to_inchi_key`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MolToInchiKeyOutput {
    pub key: Vec<u8>,
    pub diagnostics: Vec<InchiDiagnostic>,
}

/// Output of [`inchi_to_inchi_key`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InchiToInchiKeyOutput {
    pub key: Vec<u8>,
    pub diagnostics: Vec<InchiDiagnostic>,
}

/// Output of [`mol_from_inchi`].
#[derive(Clone, Debug, PartialEq)]
pub struct MolFromInchiOutput {
    pub molecule: Option<InchiMolecule>,
    pub return_values: InchiReturnValues,
    pub diagnostics: Vec<InchiDiagnostic>,
}

fn source_error(operation: &'static str, error: SourceHeapError) -> InchiError {
    let kind = match error {
        SourceHeapError::AllocationFailed => InchiErrorKind::AllocationFailed,
        SourceHeapError::UnsupportedSourceBehavior => InchiErrorKind::UnsupportedState,
        _ => InchiErrorKind::SourcePort,
    };
    InchiError {
        operation,
        kind,
        detail: format!("{error:?}"),
    }
}

fn toolkit_error(operation: &'static str, error: source::rdkit::inchi::AdapterToolkitError) -> InchiError {
    InchiError {
        operation,
        kind: if error.kind == "MolSanitizeException" {
            InchiErrorKind::SanitizeFailed
        } else {
            InchiErrorKind::Toolkit
        },
        detail: format!("{}: {}", error.kind, error.message),
    }
}

fn generation_error(operation: &'static str, error: MolToInchiError) -> InchiError {
    match error {
        MolToInchiError::Source(error) => source_error(operation, error),
        MolToInchiError::Toolkit(error) => toolkit_error(operation, error),
        MolToInchiError::InvalidOptions => InchiError {
            operation,
            kind: InchiErrorKind::InvalidInput,
            detail: "invalid InChI option bytes".to_owned(),
        },
        MolToInchiError::InvalidConformer => InchiError {
            operation,
            kind: InchiErrorKind::InvalidInput,
            detail: "the first conformer row count does not match the atom count".to_owned(),
        },
        MolToInchiError::ElementSymbolTooLong => InchiError {
            operation,
            kind: InchiErrorKind::InvalidInput,
            detail: "the toolkit element symbol does not fit the official InChI field".to_owned(),
        },
    }
}

fn key_error(operation: &'static str, error: InchiToInchiKeyError) -> InchiError {
    match error {
        InchiToInchiKeyError::Source(error) => source_error(operation, error),
        InchiToInchiKeyError::InvalidSourceOutput(detail) => InchiError {
            operation,
            kind: InchiErrorKind::InvalidSourceOutput,
            detail: detail.to_owned(),
        },
    }
}

fn structure_error(operation: &'static str, error: InchiToMolError) -> InchiError {
    match error {
        InchiToMolError::Source(error) => source_error(operation, error),
        InchiToMolError::Toolkit(error) => toolkit_error(operation, error),
        InchiToMolError::InvalidSourceOutput(detail) => InchiError {
            operation,
            kind: InchiErrorKind::InvalidSourceOutput,
            detail: detail.to_owned(),
        },
        InchiToMolError::Cleanup(error) => InchiError {
            operation,
            kind: InchiErrorKind::InvalidSourceOutput,
            detail: format!("{error:?}"),
        },
        InchiToMolError::BondIndex(error) => InchiError {
            operation,
            kind: InchiErrorKind::InvalidSourceOutput,
            detail: format!("{error:?}"),
        },
    }
}

fn nul_terminated_options(options: Option<&[u8]>) -> Result<Option<Vec<u8>>, InchiError> {
    let Some(options) = options else {
        return Ok(None);
    };
    let mut terminated = options.to_vec();
    terminated.try_reserve(1).map_err(|_| InchiError {
        operation: "mol_to_inchi",
        kind: InchiErrorKind::AllocationFailed,
        detail: "AllocationFailed".to_owned(),
    })?;
    terminated.push(0);
    Ok(Some(terminated))
}

/// Generates an InChI from the toolkit-neutral molecular graph.
pub fn mol_to_inchi(
    toolkit: &mut impl MolToInchiToolkit,
    molecule: &InchiMolecule,
    options: Option<&[u8]>,
) -> Result<MolToInchiOutput, InchiError> {
    let options = nul_terminated_options(options)?;
    let mut heap = SourceHeap::default();
    let mut engine = SourceInchiGenerationEngine::new(&mut heap);
    let mut return_values = source::rdkit::inchi::ExtraInchiReturnValues::default();
    let output =
        source::rdkit::inchi::mol_to_inchi(&mut engine, toolkit, molecule, &mut return_values, options.as_deref())
            .map_err(|error| generation_error("mol_to_inchi", error))?;
    Ok(MolToInchiOutput {
        inchi: output.inchi,
        return_values: InchiReturnValues {
            return_code: return_values.return_code,
            message: return_values.message,
            log: return_values.log,
            aux_info: return_values.aux_info,
        },
        diagnostics: output.diagnostics,
    })
}

/// Generates an InChIKey from the toolkit-neutral molecular graph.
pub fn mol_to_inchi_key(
    toolkit: &mut impl MolToInchiToolkit,
    molecule: &InchiMolecule,
    options: Option<&[u8]>,
) -> Result<MolToInchiKeyOutput, InchiError> {
    let options = nul_terminated_options(options).map_err(|mut error| {
        error.operation = "mol_to_inchi_key";
        error
    })?;
    let mut generation_heap = SourceHeap::default();
    let mut key_heap = SourceHeap::default();
    let mut generation_engine = SourceInchiGenerationEngine::new(&mut generation_heap);
    let mut key_engine = SourceInchiKeyEngine::new(&mut key_heap);
    let output = source::rdkit::inchi::mol_to_inchi_key(
        &mut generation_engine,
        &mut key_engine,
        toolkit,
        molecule,
        options.as_deref(),
    )
    .map_err(|error| match error {
        MolToInchiKeyError::MolToInchi(error) => generation_error("mol_to_inchi_key", error),
        MolToInchiKeyError::InchiToInchiKey(error) => key_error("mol_to_inchi_key", error),
    })?;
    Ok(MolToInchiKeyOutput {
        key: output.key,
        diagnostics: output.diagnostics,
    })
}

/// Generates an InChIKey directly from an InChI byte string.
pub fn inchi_to_inchi_key(inchi: &[u8]) -> Result<InchiToInchiKeyOutput, InchiError> {
    let mut heap = SourceHeap::default();
    let mut engine = SourceInchiKeyEngine::new(&mut heap);
    let output = source::rdkit::inchi::inchi_to_inchi_key(&mut engine, inchi)
        .map_err(|error| key_error("inchi_to_inchi_key", error))?;
    Ok(InchiToInchiKeyOutput {
        key: output.key,
        diagnostics: output.diagnostics,
    })
}

/// Parses an InChI into the toolkit-neutral molecular graph.
pub fn mol_from_inchi(
    toolkit: &mut impl InchiToMolToolkit,
    inchi: &[u8],
    sanitize: bool,
    remove_hs: bool,
) -> Result<MolFromInchiOutput, InchiError> {
    let mut heap = SourceHeap::default();
    let mut engine = SourceInchiStructureEngine::new(&mut heap);
    let mut return_values = source::rdkit::inchi::ExtraInchiReturnValues::default();
    let output =
        source::rdkit::inchi::inchi_to_mol(&mut engine, toolkit, inchi, &mut return_values, sanitize, remove_hs)
            .map_err(|error| structure_error("mol_from_inchi", error))?;
    Ok(MolFromInchiOutput {
        molecule: output.molecule,
        return_values: InchiReturnValues {
            return_code: return_values.return_code,
            message: return_values.message,
            log: return_values.log,
            aux_info: return_values.aux_info,
        },
        diagnostics: output.diagnostics,
    })
}

#[cfg(test)]
mod scalar_api_tests {
    use super::*;

    #[derive(Default)]
    struct ScalarToolkit {
        fail_generation: bool,
    }

    fn toolkit_failure(call: &'static str) -> InchiToolkitError {
        InchiToolkitError {
            kind: "scalar test toolkit error",
            message: call.to_owned(),
        }
    }

    fn atomic_number(element: &[u8]) -> Option<i32> {
        match element {
            b"H" => Some(1),
            b"B" => Some(5),
            b"C" => Some(6),
            b"N" => Some(7),
            b"O" => Some(8),
            b"F" => Some(9),
            b"P" => Some(15),
            b"S" => Some(16),
            b"Cl" => Some(17),
            b"Br" => Some(35),
            b"I" => Some(53),
            _ => None,
        }
    }

    fn element_symbol(atomic_number: i32) -> Option<&'static [u8]> {
        match atomic_number {
            1 => Some(b"H"),
            5 => Some(b"B"),
            6 => Some(b"C"),
            7 => Some(b"N"),
            8 => Some(b"O"),
            9 => Some(b"F"),
            15 => Some(b"P"),
            16 => Some(b"S"),
            17 => Some(b"Cl"),
            35 => Some(b"Br"),
            53 => Some(b"I"),
            _ => None,
        }
    }

    fn average_weight(atomic_number: i32) -> Option<f64> {
        match atomic_number {
            1 => Some(1.008),
            5 => Some(10.81),
            6 => Some(12.011),
            7 => Some(14.007),
            8 => Some(15.999),
            9 => Some(18.998),
            15 => Some(30.974),
            16 => Some(32.06),
            17 => Some(35.45),
            35 => Some(79.904),
            53 => Some(126.904),
            _ => None,
        }
    }

    impl MolToInchiToolkit for ScalarToolkit {
        fn needs_update_property_cache(&mut self, _molecule: &InchiMolecule) -> Result<bool, InchiToolkitError> {
            if self.fail_generation {
                return Err(toolkit_failure("needs_update_property_cache"));
            }
            Ok(false)
        }

        fn update_property_cache(
            &mut self,
            _molecule: &mut InchiMolecule,
            strict: bool,
        ) -> Result<(), InchiToolkitError> {
            assert!(!strict);
            Ok(())
        }

        fn kekulize(&mut self, _molecule: &mut InchiMolecule, mark_atoms_bonds: bool) -> Result<(), InchiToolkitError> {
            assert!(!mark_atoms_bonds);
            Ok(())
        }

        fn element_symbol(&mut self, atomic_number: i32) -> Result<Vec<u8>, InchiToolkitError> {
            element_symbol(atomic_number)
                .map(<[u8]>::to_vec)
                .ok_or_else(|| toolkit_failure("element_symbol"))
        }

        fn atomic_weight(&mut self, atomic_number: i32) -> Result<f64, InchiToolkitError> {
            average_weight(atomic_number).ok_or_else(|| toolkit_failure("atomic_weight"))
        }

        fn total_num_hydrogens(&mut self, molecule: &InchiMolecule, atom_index: u32) -> Result<u32, InchiToolkitError> {
            molecule
                .atoms()
                .get(atom_index as usize)
                .map(|atom| atom.num_explicit_hydrogens)
                .ok_or_else(|| toolkit_failure("total_num_hydrogens"))
        }

        fn calc_implicit_valence(
            &mut self,
            _molecule: &mut InchiMolecule,
            _atom_index: u32,
        ) -> Result<i32, InchiToolkitError> {
            Ok(0)
        }

        fn total_degree(&mut self, molecule: &InchiMolecule, atom_index: u32) -> Result<u32, InchiToolkitError> {
            Ok(molecule
                .bonds()
                .iter()
                .filter(|bond| bond.begin_atom_index() == atom_index || bond.end_atom_index() == atom_index)
                .count() as u32)
        }
    }

    impl InchiToMolToolkit for ScalarToolkit {
        fn atomic_number(&mut self, element: &[u8]) -> Result<i32, InchiToolkitError> {
            atomic_number(element).ok_or_else(|| toolkit_failure("atomic_number"))
        }

        fn average_atomic_weight(&mut self, atomic_number: i32) -> Result<f64, InchiToolkitError> {
            average_weight(atomic_number).ok_or_else(|| toolkit_failure("average_atomic_weight"))
        }

        fn update_property_cache(
            &mut self,
            _molecule: &mut InchiMolecule,
            strict: bool,
        ) -> Result<(), InchiToolkitError> {
            assert!(!strict);
            Ok(())
        }

        fn assign_atom_cip_ranks(&mut self, molecule: &mut InchiMolecule) -> Result<Vec<u32>, InchiToolkitError> {
            Ok((0..molecule.atoms().len() as u32).collect())
        }

        fn remove_hydrogens(&mut self, _molecule: &mut InchiMolecule) -> Result<(), InchiToolkitError> {
            Ok(())
        }

        fn sanitize_molecule(&mut self, _molecule: &mut InchiMolecule) -> Result<(), InchiToolkitError> {
            Ok(())
        }

        fn assign_stereochemistry(
            &mut self,
            _molecule: &mut InchiMolecule,
            clean_it: bool,
            force: bool,
        ) -> Result<(), InchiToolkitError> {
            assert!(clean_it);
            assert!(force);
            Ok(())
        }
    }

    fn carbon(isotope: u32, chiral_tag: InchiChiralTag) -> InchiAtom {
        InchiAtom {
            atomic_number: 6,
            isotope,
            chiral_tag,
            ..InchiAtom::default()
        }
    }

    fn bromochlorofluoromethane(chiral_tag: InchiChiralTag) -> InchiMolecule {
        let mut center = carbon(0, chiral_tag);
        center.num_explicit_hydrogens = 1;
        center.no_implicit = true;
        InchiMolecule::try_from_graph(
            vec![
                center,
                InchiAtom {
                    atomic_number: 9,
                    ..InchiAtom::default()
                },
                InchiAtom {
                    atomic_number: 17,
                    ..InchiAtom::default()
                },
                InchiAtom {
                    atomic_number: 35,
                    ..InchiAtom::default()
                },
            ],
            vec![
                InchiBond::new(0, 1, InchiBondType::Single),
                InchiBond::new(0, 2, InchiBondType::Single),
                InchiBond::new(0, 3, InchiBondType::Single),
            ],
            Vec::new(),
        )
        .unwrap()
    }

    #[test]
    fn inchi_core_scalar_api__inchi_key_success_and_error_diagnostic() {
        let success = inchi_to_inchi_key(b"InChI=1S/CH4/h1H4").unwrap();
        assert_eq!(success.key, b"VNWKTOKETHGBQD-UHFFFAOYSA-N");
        assert!(success.diagnostics.is_empty());

        let invalid = inchi_to_inchi_key(b"").unwrap();
        assert!(invalid.key.is_empty());
        assert_eq!(invalid.diagnostics.len(), 1);
        assert_eq!(invalid.diagnostics[0].level, InchiDiagnosticLevel::Error);
        assert_eq!(
            invalid.diagnostics[0].message,
            "Invalid InChI prefix in generating InChI Key\n"
        );
    }

    #[test]
    fn inchi_core_scalar_api__generation_preserves_input_options_coordinates_and_isotope() {
        let molecule = InchiMolecule::try_from_graph(
            vec![carbon(13, InchiChiralTag::Unspecified)],
            Vec::new(),
            vec![vec![[1.25, -2.5, 3.75]], vec![[9.0, 9.0, 9.0]]],
        )
        .unwrap();
        let before = molecule.clone();
        let output = mol_to_inchi(&mut ScalarToolkit::default(), &molecule, Some(b"-AuxNone")).unwrap();

        assert_eq!(molecule, before);
        assert_eq!(output.return_values.return_code, 0);
        assert_eq!(output.inchi, b"InChI=1S/CH4/h1H4/i1+1");
        assert!(output.return_values.aux_info.is_empty());
        assert!(output.diagnostics.is_empty());
    }

    #[test]
    fn inchi_core_scalar_api__generation_preserves_single_center_relative_and_racemic_stereo() {
        let molecule = bromochlorofluoromethane(InchiChiralTag::TetrahedralCw);

        let absolute = mol_to_inchi(&mut ScalarToolkit::default(), &molecule, Some(b"-AuxNone")).unwrap();
        let relative = mol_to_inchi(&mut ScalarToolkit::default(), &molecule, Some(b"-AuxNone -SRel")).unwrap();
        let racemic = mol_to_inchi(&mut ScalarToolkit::default(), &molecule, Some(b"-AuxNone -SRac")).unwrap();

        assert_eq!(absolute.inchi, b"InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1");
        assert_eq!(relative.inchi, b"InChI=1/CHBrClF/c2-1(3)4/h1H/t1-/s2");
        assert_eq!(racemic.inchi, b"InChI=1/CHBrClF/c2-1(3)4/h1H/t1-/s3");
    }

    #[test]
    fn inchi_core_scalar_api__generation_warning_and_structured_errors() {
        let chiral =
            InchiMolecule::try_from_graph(vec![carbon(0, InchiChiralTag::TetrahedralCw)], Vec::new(), Vec::new())
                .unwrap();
        let warning = mol_to_inchi(&mut ScalarToolkit::default(), &chiral, None).unwrap();
        assert_eq!(warning.diagnostics.len(), 1);
        assert_eq!(warning.diagnostics[0].level, InchiDiagnosticLevel::Warning);
        assert_eq!(
            warning.diagnostics[0].message,
            "tetrahedral chirality on atom with <3 or >4 neighbors will be ignored.\n"
        );

        let invalid_conformer = InchiMolecule::try_from_graph(
            vec![carbon(0, InchiChiralTag::Unspecified)],
            Vec::new(),
            vec![Vec::new()],
        )
        .unwrap();
        let error = mol_to_inchi(&mut ScalarToolkit::default(), &invalid_conformer, None).unwrap_err();
        assert_eq!(error.operation, "mol_to_inchi");
        assert_eq!(error.kind, InchiErrorKind::InvalidInput);

        let toolkit_error = mol_to_inchi(&mut ScalarToolkit { fail_generation: true }, &chiral, None).unwrap_err();
        assert_eq!(toolkit_error.kind, InchiErrorKind::Toolkit);

        assert_eq!(
            source_error("mol_to_inchi", SourceHeapError::AllocationFailed).kind,
            InchiErrorKind::AllocationFailed
        );
        assert_eq!(
            source_error("mol_from_inchi", SourceHeapError::UnsupportedSourceBehavior).kind,
            InchiErrorKind::UnsupportedState
        );
    }

    #[test]
    fn inchi_core_scalar_api__graph_validation_and_topology_replacement_are_atomic() {
        let invalid = InchiMolecule::try_from_graph(
            vec![carbon(0, InchiChiralTag::Unspecified)],
            vec![InchiBond::new(0, 1, InchiBondType::Single)],
            Vec::new(),
        )
        .unwrap_err();
        assert_eq!(invalid.bond_index, 0);
        assert_eq!(invalid.atom_index, 1);
        assert_eq!(invalid.atom_count, 1);

        let mut molecule =
            InchiMolecule::try_from_graph(vec![carbon(0, InchiChiralTag::Unspecified)], Vec::new(), Vec::new())
                .unwrap();
        let before = molecule.clone();
        assert!(
            molecule
                .replace_graph(
                    vec![carbon(0, InchiChiralTag::Unspecified)],
                    vec![InchiBond::new(0, 2, InchiBondType::Single)],
                    Vec::new(),
                )
                .is_err()
        );
        assert_eq!(molecule, before);
    }

    #[test]
    fn inchi_core_scalar_api__mol_from_inchi_preserves_isotope_and_stereo_fields() {
        let mut toolkit = ScalarToolkit::default();
        let isotopic = mol_from_inchi(&mut toolkit, b"InChI=1S/CHBrClF/c2-1(3)4/t1-/m0/s1/i1+1", true, false).unwrap();
        assert_eq!(isotopic.return_values.return_code, 1);
        let molecule = isotopic.molecule.expect("successful parse must return a graph");
        assert_eq!(molecule.atoms().len(), 4);
        assert_eq!(molecule.bonds().len(), 3);
        let center = &molecule.atoms()[0];
        assert_eq!(center.atomic_number, 6);
        assert_eq!(center.isotope, 13);
        assert!(matches!(
            center.chiral_tag,
            InchiChiralTag::TetrahedralCw | InchiChiralTag::TetrahedralCcw
        ));
        assert!(molecule.conformers().is_empty());
    }

    #[test]
    fn inchi_core_scalar_api__mol_from_inchi_parses_phosphoserine_graph() {
        const PHOSPHOSERINE: &[u8] =
            b"InChI=1S/C3H8NO6P/c4-2(3(5)6)1-10-11(7,8)9/h2H,1,4H2,(H,5,6)(H2,7,8,9)/t2-/m0/s1";

        let unsanitized = mol_from_inchi(&mut ScalarToolkit::default(), PHOSPHOSERINE, false, false)
            .expect("source-defined phosphoserine InChI must parse");
        let molecule = unsanitized.molecule.expect("successful parse must return a graph");
        assert_eq!(molecule.atoms().len(), 12);
        assert_eq!(molecule.bonds().len(), 11);
        assert_eq!(
            molecule.atoms().iter().filter(|atom| atom.atomic_number == 1).count(),
            1
        );

        assert_eq!(molecule.atoms().iter().map(|atom| atom.formal_charge).sum::<i32>(), 0);
        assert!(matches!(
            molecule.atoms()[1].chiral_tag,
            InchiChiralTag::TetrahedralCw | InchiChiralTag::TetrahedralCcw
        ));
        assert_eq!(
            molecule
                .bonds()
                .iter()
                .filter(|bond| bond.bond_type == InchiBondType::Double)
                .count(),
            2
        );
    }

    #[test]
    fn inchi_core_scalar_api__mol_from_inchi_parses_guanidinium_alkaloid_hydrochloride() {
        const GUANIDINIUM_ALKALOID_HYDROCHLORIDE: &[u8] = b"InChI=1S/C22H33N3O4.ClH/c1-3-16-8-4-5-11-21(29-16)13-15-9-10-17-18(19(26)27)22(12-6-7-14(2)28-22)24-20(23-21)25(15)17;/h4,8,14-18H,3,5-7,9-13H2,1-2H3,(H2,23,24,26,27);1H/t14-,15+,16+,17-,18+,21+,22-;/m1./s1";

        let parsed = mol_from_inchi(
            &mut ScalarToolkit::default(),
            GUANIDINIUM_ALKALOID_HYDROCHLORIDE,
            false,
            false,
        )
        .expect("source-defined guanidinium alkaloid hydrochloride must parse");
        let molecule = parsed.molecule.expect("successful parse must return a graph");
        assert_eq!(molecule.atoms().len(), 35);
        assert_eq!(molecule.bonds().len(), 38);
        assert_eq!(
            molecule.atoms().iter().filter(|atom| atom.atomic_number == 1).count(),
            5
        );
    }
}

#[cfg(test)]
mod test_support;
