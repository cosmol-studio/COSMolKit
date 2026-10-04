//! Private storage and generated transactions; sibling operation bodies see only access structs.
use crate::*;
use cosmolkit_bio::BioStructureData;

#[derive(Debug, Clone, PartialEq)]
pub struct BioStructure {
    data: BioStructureData,
}
/// An amino-acid-only projection with one owned, already filtered [`BioStructure`].
///
/// `as_bio_structure` borrows that structure; `into_bio_structure` consumes
/// the wrapper. Neither operation restores non-protein rows removed when the
/// projection was created.
///
/// ```
/// use cosmolkit::{BioStructure, Protein};
/// let pdb = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n";
/// let protein = Protein::from_pdb(pdb).unwrap();
/// let view: &BioStructure = protein.as_bio_structure();
/// assert_eq!(view.num_atoms(), protein.num_atoms());
/// assert_eq!(protein.selection_summary().atoms, protein.num_atoms());
/// let owned: BioStructure = protein.into_bio_structure();
/// assert_eq!(owned.num_atoms(), 1);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Protein {
    structure: BioStructure,
}

impl BioStructure {
    /// Validate and construct a structure from detached BIO parts.
    pub fn from_parts(parts: BioStructureParts) -> Result<Self, BioStructureError> {
        Ok(Self {
            data: BioStructureData::from_parts(parts)?,
        })
    }
    pub fn validate_parts(parts: &BioStructureParts) -> Result<(), BioStructureError> {
        BioStructureData::validate_parts(parts)
    }
    pub fn validate(&self) -> Result<(), BioStructureError> {
        self.data.validate()
    }
    pub fn into_parts(self) -> BioStructureParts {
        self.data.into_parts()
    }
    pub fn input_format(&self) -> BioCoordinateFormat {
        self.data.input_format()
    }

    /// Registered Experimental detached query (BIO-CID C29): original
    /// atom IDs accepted by `selection`, in source order, through the
    /// canonical BIO-CID C24 lazy cursor over this structure's borrowed
    /// row tables. Read-only — no mutation operation, no COW; the thin
    /// root delegate never duplicates the matcher.
    pub fn selected_atom_ids(
        &self,
        selection: &crate::BioSelection,
    ) -> Result<Vec<cosmolkit_bio::BioAtomId>, crate::BioSelectionMatchError> {
        selection.data.selected_bio_atom_ids(
            self.data.models(),
            self.data.chains(),
            self.data.residues(),
            self.data.atoms(),
            self.data.input_format(),
        )
    }
    /// Read structural text with the IO owner's explicit dispatch and source context.
    pub fn from_text_with_params(text: &str, params: &BioReadParams) -> Result<Self, BioReadError> {
        Ok(Self {
            data: cosmolkit_io::read_bio_structure(text, params)?,
        })
    }
    /// Read structural text using the default content-detection parameters.
    pub fn from_text(text: &str) -> Result<Self, BioReadError> {
        Self::from_text_with_params(text, &BioReadParams::default())
    }
    /// Read a regular structural file using the selected coordinate format.
    pub fn read_with_format(
        path: &std::path::Path,
        format: BioCoordinateFormat,
    ) -> Result<Self, BioReadError> {
        Ok(Self {
            data: cosmolkit_io::read_bio_structure_file(path, format)?,
        })
    }
    /// Read a regular structural file, selecting its format from its extension.
    pub fn read(path: &std::path::Path) -> Result<Self, BioReadError> {
        Self::read_with_format(path, BioCoordinateFormat::Unknown)
    }
    /// Experimental Gemmi-aligned PDB text reader; source scope and limits are inherited.
    pub fn from_pdb(text: &str) -> Result<Self, BioPdbReadError> {
        Self::from_pdb_with_params(text, &BioPdbReadParams::default())
    }
    pub fn from_pdb_with_params(
        text: &str,
        params: &BioPdbReadParams,
    ) -> Result<Self, BioPdbReadError> {
        Ok(Self {
            data: cosmolkit_io::read_pdb_bio_structure(text, "<string>", params)?,
        })
    }
    /// Experimental coordinate mmCIF reader; no chemical Molecule conversion.
    pub fn from_mmcif(text: &str) -> Result<Self, BioMmcifReadError> {
        Ok(Self {
            data: cosmolkit_io::read_mmcif_bio_structure(text, "<string>")?,
        })
    }
    /// Serialize structure coordinates to mmCIF text with explicit options.
    ///
    /// Experimental: this emits ONLY the coordinate profile — the block
    /// name/`_entry.id`, `_atom_site` and the conditional
    /// `_atom_site_anisotrop`. Crystal/symmetry, NCS, assembly, connection,
    /// cis-peptide, refinement and all other source categories are never
    /// written; no lossless roundtrip is claimed.
    pub fn to_mmcif_with_params(
        &self,
        params: &crate::BioMmcifWriteParams,
    ) -> Result<String, crate::BioMmcifWriteError> {
        cosmolkit_io::bio_structure_to_mmcif_text(&self.data, params)
    }
    /// Serialize structure coordinates to mmCIF text with default options.
    ///
    /// Experimental: same coordinate-only category set and exclusions as
    /// [`BioStructure::to_mmcif_with_params`].
    pub fn to_mmcif(&self) -> Result<String, crate::BioMmcifWriteError> {
        self.to_mmcif_with_params(&crate::BioMmcifWriteParams::default())
    }
    /// Write structure coordinates to an mmCIF file with explicit options.
    ///
    /// Experimental: same coordinate-only category set and exclusions as
    /// [`BioStructure::to_mmcif_with_params`]; the document is fully
    /// serialized before the destination is created or truncated, and write
    /// failures retain the path and underlying IO error.
    pub fn write_mmcif_with_params(
        &self,
        path: &std::path::Path,
        params: &crate::BioMmcifWriteParams,
    ) -> Result<(), crate::BioMmcifWriteError> {
        cosmolkit_io::write_bio_structure_mmcif_file(&self.data, params, path)
    }
    /// Write structure coordinates to an mmCIF file with default options.
    ///
    /// Experimental: same coordinate-only category set and exclusions as
    /// [`BioStructure::to_mmcif_with_params`].
    pub fn write_mmcif(&self, path: &std::path::Path) -> Result<(), crate::BioMmcifWriteError> {
        self.write_mmcif_with_params(path, &crate::BioMmcifWriteParams::default())
    }

    /// Serialize structure coordinates to PDB text with explicit options.
    ///
    /// Experimental: this emits ONLY the seven-record coordinate
    /// projection — ATOM, HETATM, ANISOU, TER, MODEL, ENDMDL and the
    /// optional END record. HEADER, TITLE, SEQRES, SSBOND, LINK, CISPEP,
    /// CRYST1, ORIGX, SCALE, MTRIX, HETNAM, CONECT, assemblies, refinement
    /// and all other source categories are never written; no lossless
    /// roundtrip is claimed.
    pub fn to_pdb_with_params(
        &self,
        params: &crate::BioPdbWriteParams,
    ) -> Result<String, crate::BioPdbWriteError> {
        cosmolkit_io::bio_structure_to_pdb_text(&self.data, params)
    }

    /// Serialize structure coordinates to PDB text with default options.
    ///
    /// Experimental: same seven-record coordinate projection and
    /// exclusions as [`BioStructure::to_pdb_with_params`].
    pub fn to_pdb(&self) -> Result<String, crate::BioPdbWriteError> {
        self.to_pdb_with_params(&crate::BioPdbWriteParams::default())
    }

    /// Write structure coordinates to a PDB file with explicit options.
    ///
    /// Experimental: same seven-record coordinate projection as
    /// [`BioStructure::to_pdb_with_params`]; the text is fully produced
    /// before the destination is created or truncated, and write failures
    /// retain the path and underlying IO error.
    pub fn write_pdb_with_params(
        &self,
        path: &std::path::Path,
        params: &crate::BioPdbWriteParams,
    ) -> Result<(), crate::BioPdbWriteError> {
        cosmolkit_io::write_bio_structure_pdb_file(&self.data, path, params)
    }

    /// Write structure coordinates to a PDB file with default options.
    ///
    /// Experimental: same seven-record coordinate projection as
    /// [`BioStructure::to_pdb_with_params`].
    pub fn write_pdb(&self, path: &std::path::Path) -> Result<(), crate::BioPdbWriteError> {
        self.write_pdb_with_params(path, &crate::BioPdbWriteParams::default())
    }
    pub fn protein(&self) -> Result<Protein, ProteinProjectionError> {
        Ok(Protein {
            structure: BioStructure {
                data: self.data.protein()?.into_structure(),
            },
        })
    }
    fn operation_data_mut(&mut self) -> &mut BioStructureData {
        &mut self.data
    }
    fn validate_operation(&self) -> Result<(), BioOperationError> {
        self.data.validate().map_err(BioOperationError::Structure)
    }
    pub fn models(&self) -> &[cosmolkit_bio::BioModelRow] {
        self.data.models()
    }
    pub fn num_models(&self) -> usize {
        self.data.models().len()
    }
    pub fn chains(&self) -> &[cosmolkit_bio::BioChainRow] {
        self.data.chains()
    }
    pub fn num_chains(&self) -> usize {
        self.data.chains().len()
    }
    pub fn residues(&self) -> &[cosmolkit_bio::BioResidueRow] {
        self.data.residues()
    }
    pub fn num_residues(&self) -> usize {
        self.data.residues().len()
    }
    pub fn atoms(&self) -> &[cosmolkit_bio::BioAtomRow] {
        self.data.atoms()
    }
    pub fn atom_position(&self, atom: BioAtomId) -> Option<[f64; 3]> {
        self.data.atom_position(atom)
    }
    pub fn residue_atoms(&self, residue: BioResidueId) -> Option<&[BioAtomRow]> {
        self.data.residue_atoms(residue)
    }
    pub fn num_atoms(&self) -> usize {
        self.data.atoms().len()
    }
    pub fn entities(&self) -> &[cosmolkit_bio::BioEntityRow] {
        self.data.entities()
    }
    pub fn num_entities(&self) -> usize {
        self.data.entities().len()
    }
    pub fn connections(&self) -> &[cosmolkit_bio::BioConnection] {
        self.data.connections()
    }
    pub fn cispeps(&self) -> &[cosmolkit_bio::BioCisPep] {
        self.data.cispeps()
    }
    pub fn mod_residues(&self) -> &[cosmolkit_bio::BioModRes] {
        self.data.mod_residues()
    }
    pub fn helices(&self) -> &[cosmolkit_bio::BioHelix] {
        self.data.helices()
    }
    pub fn sheets(&self) -> &[cosmolkit_bio::BioSheet] {
        self.data.sheets()
    }
    pub fn metadata(&self) -> &cosmolkit_bio::BioMetadata {
        self.data.metadata()
    }
    pub fn source_state(&self) -> &cosmolkit_bio::BioStructureSourceState {
        self.data.source_state()
    }
    pub fn name(&self) -> &str {
        &self.data.source_state().name
    }
    pub fn has_origx(&self) -> bool {
        self.data.source_state().has_origx
    }
    pub fn origx(&self) -> &BioTransform {
        &self.data.source_state().origx
    }
    /// Borrow the last identity NCS operation ID retained in source metadata.
    pub fn ncs_oper_identity_id(&self) -> Option<&str> {
        self.data.ncs_oper_identity_id()
    }
    pub fn resolution(&self) -> f64 {
        self.data.source_state().resolution
    }
    pub fn ter_status(&self) -> u8 {
        self.data.source_state().ter_status
    }
    pub fn coordinates(&self) -> &cosmolkit_bio::BioCoordinateBlock {
        self.data.coordinates()
    }
    pub fn crystal(&self) -> Option<&cosmolkit_bio::BioCrystalInfo> {
        self.data.crystal()
    }
    pub fn ncs_operators(&self) -> &[cosmolkit_bio::BioNcsOperator] {
        self.data.ncs_operators()
    }
    pub fn assemblies(&self) -> &[cosmolkit_bio::BioAssembly] {
        self.data.assemblies()
    }
    pub fn find_entity(&self, source_id: &str) -> Option<(BioEntityId, &BioEntityRow)> {
        self.data.find_entity(source_id)
    }
    pub fn find_entity_of_subchain(&self, subchain: &str) -> Option<(BioEntityId, &BioEntityRow)> {
        self.data.find_entity_of_subchain(subchain)
    }
    pub fn find_atom(
        &self,
        residue_id: BioResidueId,
        name: AtomName,
        request: AltLocRequest,
        element: Option<Element>,
    ) -> Option<(BioAtomId, &BioAtomRow)> {
        self.data.find_atom(residue_id, name, request, element)
    }
    pub fn atom_by_altloc(
        &self,
        residue_id: BioResidueId,
        name: AtomName,
        altloc: Option<AltLocLabel>,
    ) -> Result<(BioAtomId, &BioAtomRow), BioStructureError> {
        self.data.atom_by_altloc(residue_id, name, altloc)
    }
}

impl Protein {
    /// Return the original input format retained by the underlying structure.
    ///
    /// This reads source metadata; it does not infer a format from the coordinates.
    pub fn input_format(&self) -> BioCoordinateFormat {
        self.structure.input_format()
    }

    /// Borrow the already filtered structural hierarchy without cloning it.
    ///
    /// The view is immutable; it cannot become a mutable structure reference.
    ///
    /// ```compile_fail
    /// use cosmolkit::{BioStructure, Protein};
    /// fn invalid(protein: &mut Protein) {
    ///     let _: &mut BioStructure = protein.as_bio_structure();
    /// }
    /// ```
    pub fn as_bio_structure(&self) -> &BioStructure {
        &self.structure
    }

    /// Registered Experimental detached query (BIO-CID C30): the SAME
    /// canonical query as BioStructure::selected_atom_ids (C29), reached
    /// through the protein's existing borrowed BioStructure — never a
    /// duplicate matcher.
    pub fn selected_atom_ids(
        &self,
        selection: &crate::BioSelection,
    ) -> Result<Vec<cosmolkit_bio::BioAtomId>, crate::BioSelectionMatchError> {
        self.as_bio_structure().selected_atom_ids(selection)
    }

    /// Consume this protein and return its already filtered structural hierarchy.
    pub fn into_bio_structure(self) -> BioStructure {
        self.structure
    }

    /// Read a structure and project its amino-acid hierarchy.
    pub fn from_text_with_params(
        text: &str,
        params: &BioReadParams,
    ) -> Result<Self, ProteinReadError> {
        BioStructure::from_text_with_params(text, params)
            .map_err(ProteinReadError::Structure)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    /// Detect structural content and retain only its amino-acid hierarchy.
    pub fn from_text(text: &str) -> Result<Self, ProteinReadError> {
        BioStructure::from_text(text)
            .map_err(ProteinReadError::Structure)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    /// Read a structural file with explicit format and project its amino acids.
    pub fn read_with_format(
        path: &std::path::Path,
        format: BioCoordinateFormat,
    ) -> Result<Self, ProteinReadError> {
        BioStructure::read_with_format(path, format)
            .map_err(ProteinReadError::Structure)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    /// Read a structural file by extension and project its amino acids.
    pub fn read(path: &std::path::Path) -> Result<Self, ProteinReadError> {
        BioStructure::read(path)
            .map_err(ProteinReadError::Structure)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    pub fn from_pdb(text: &str) -> Result<Self, ProteinReadError> {
        Self::from_pdb_with_params(text, &BioPdbReadParams::default())
    }
    pub fn from_pdb_with_params(
        text: &str,
        params: &BioPdbReadParams,
    ) -> Result<Self, ProteinReadError> {
        BioStructure::from_pdb_with_params(text, params)
            .map_err(ProteinReadError::Pdb)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    pub fn from_mmcif(text: &str) -> Result<Self, ProteinReadError> {
        BioStructure::from_mmcif(text)
            .map_err(ProteinReadError::Mmcif)?
            .protein()
            .map_err(ProteinReadError::Projection)
    }
    fn operation_data_mut(&mut self) -> &mut BioStructureData {
        &mut self.structure.data
    }
    fn validate_operation(&self) -> Result<(), BioOperationError> {
        cosmolkit_bio::validate_protein_structure(&self.structure.data)
            .map_err(BioOperationError::Protein)
    }
    pub fn num_models(&self) -> usize {
        self.structure.models().len()
    }
    pub fn num_chains(&self) -> usize {
        self.structure.chains().len()
    }
    pub fn num_residues(&self) -> usize {
        self.structure.residues().len()
    }
    pub fn num_atoms(&self) -> usize {
        self.structure.atoms().len()
    }
    pub fn selection_summary(&self) -> ProteinSelectionSummary {
        cosmolkit_bio::protein_selection_summary(&self.structure.data)
    }
    pub fn chains(&self) -> ProteinChainIter<'_> {
        cosmolkit_bio::protein_chains(&self.structure.data)
    }
    pub fn residues(&self) -> ProteinResidueIter<'_> {
        cosmolkit_bio::protein_residues(&self.structure.data)
    }
    pub fn atoms(&self) -> ProteinAtomIter<'_> {
        cosmolkit_bio::protein_atoms(&self.structure.data)
    }
    pub fn chain(&self, index: usize) -> Option<ProteinChainRef<'_>> {
        cosmolkit_bio::protein_chain(&self.structure.data, index)
    }
}

#[allow(dead_code)]
#[derive(Debug)]
pub(super) struct BioOperationSpec {
    pub id: &'static str,
    pub target: &'static str,
    pub value_method: &'static str,
    pub inplace_method: &'static str,
    pub read: &'static [&'static str],
    pub write: &'static [&'static str],
    pub replace: &'static [&'static str],
}
cosmolkit_macros::bio_structure_ops! {
    selection(selection: &crate::BioSelection) {
        targets: [BioStructure, Protein],
        value: with_selection,
        inplace: retain_selection_,
        body: super::selection_impl,
        access: SelectionAccess,
        read: [
            input_format: cosmolkit_bio::BioCoordinateFormat,
            entities: std::sync::Arc<Vec<cosmolkit_bio::BioEntityRow>>,
            connections: std::sync::Arc<Vec<cosmolkit_bio::BioConnection>>,
            cispeps: std::sync::Arc<Vec<cosmolkit_bio::BioCisPep>>,
            mod_residues: std::sync::Arc<Vec<cosmolkit_bio::BioModRes>>,
            helices: std::sync::Arc<Vec<cosmolkit_bio::BioHelix>>,
            sheets: std::sync::Arc<Vec<cosmolkit_bio::BioSheet>>,
            metadata: std::sync::Arc<cosmolkit_bio::BioMetadata>,
            crystal: std::sync::Arc<Option<cosmolkit_bio::BioCrystalInfo>>,
            ncs_operators: std::sync::Arc<Vec<cosmolkit_bio::BioNcsOperator>>,
            assemblies: std::sync::Arc<Vec<cosmolkit_bio::BioAssembly>>,
        ],
        write: [],
        replace: [
            models: Vec<cosmolkit_bio::BioModelRow>,
            chains: Vec<cosmolkit_bio::BioChainRow>,
            residues: Vec<cosmolkit_bio::BioResidueRow>,
            atoms: Vec<cosmolkit_bio::BioAtomRow>,
            source_state: cosmolkit_bio::BioStructureSourceState,
            coordinates: cosmolkit_bio::BioCoordinateBlock,
        ],
    }

    translate_coordinates(offset: [f64; 3]) {
        targets: [BioStructure, Protein],
        value: with_translated_coordinates,
        inplace: translate_,
        body: super::translate_coordinates_impl,
        access: TranslateCoordinatesAccess,
        read: [],
        write: [coordinates: cosmolkit_bio::BioCoordinateBlock],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    const ATOM: &str =
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n";

    mod probes {
        use super::*;
        cosmolkit_macros::bio_structure_ops! {
            replacement_failure(mode: u8) {
                targets: [BioStructure, Protein],
                value: probe_replacement,
                inplace: probe_replacement_,
                body: super::replacement_failure_body,
                access: ReplacementFailureAccess,
                read: [],
                write: [],
                replace: [coordinates: BioCoordinateBlock],
            }
            failure(mode: u8) {
                targets: [BioStructure, Protein],
                value: probe_failure,
                inplace: probe_failure_,
                body: super::fail_body,
                access: FailureAccess,
                read: [atoms: Vec<cosmolkit_bio::BioAtomRow>],
                write: [coordinates: cosmolkit_bio::BioCoordinateBlock],
            }
            protein_guard(rows: Vec<cosmolkit_bio::BioResidueRow>) {
                targets: [Protein],
                value: probe_residues,
                inplace: probe_residues_,
                body: super::replace_residues,
                access: ResidueAccess,
                read: [],
                write: [residues: Vec<cosmolkit_bio::BioResidueRow>],
            }
        }
    }
    fn replacement_failure_body(
        access: probes::ReplacementFailureAccess<'_>,
        mode: u8,
    ) -> Result<probes::ReplacementFailureAccessReplacement, BioOperationError> {
        match mode {
            0 => Ok(probes::ReplacementFailureAccessReplacement {
                coordinates: Arc::clone(access.coordinates),
            }),
            1 => Err(BioOperationError::Structure(
                BioStructureError::AtomNotFound,
            )),
            2 => Ok(probes::ReplacementFailureAccessReplacement {
                coordinates: Arc::new(BioCoordinateBlock::default()),
            }),
            _ => panic!("intentional replacement rollback probe"),
        }
    }

    #[test]
    fn bio_selection_public_replacement_error_validation_unwind_are_atomic() {
        let source = BioStructure::from_pdb(ATOM).unwrap();
        let protein = source.protein().unwrap();
        assert_eq!(source.probe_replacement(0).unwrap(), source);
        assert_eq!(protein.probe_replacement(0).unwrap(), protein);
        for mode in [1, 2, 3] {
            let mut value = source.clone();
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                value.probe_replacement_(mode)
            }));
            if mode == 3 {
                assert!(result.is_err());
            } else {
                assert!(result.unwrap().is_err());
            }
            assert_eq!(value, source);
            assert!(Arc::ptr_eq(
                &value.data.coordinates,
                &source.data.coordinates
            ));
            assert!(Arc::ptr_eq(&value.data.atoms, &source.data.atoms));
            let mut value = protein.clone();
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                value.probe_replacement_(mode)
            }));
            if mode == 3 {
                assert!(result.is_err());
            } else {
                assert!(result.unwrap().is_err());
            }
            assert_eq!(value, protein);
            assert!(Arc::ptr_eq(
                &value.structure.data.coordinates,
                &protein.structure.data.coordinates
            ));
        }
    }

    fn fail_body(access: probes::FailureAccess<'_>, mode: u8) -> Result<(), BioOperationError> {
        assert_eq!(access.atoms.len(), access.coordinates.len());
        if mode == 0 {
            return Ok(());
        } // Declared writable does not mean mandatory writing.
        access.coordinates.positions_mut()[0][0] = 99.0;
        match mode {
            1 => Err(BioOperationError::Structure(
                BioStructureError::AtomNotFound,
            )),
            2 => {
                *access.coordinates = BioCoordinateBlock::default();
                Ok(())
            }
            _ => panic!("intentional rollback probe"),
        }
    }
    fn replace_residues(
        access: probes::ResidueAccess<'_>,
        rows: Vec<cosmolkit_bio::BioResidueRow>,
    ) -> Result<(), BioOperationError> {
        *access.residues = rows;
        Ok(())
    }

    #[test]
    fn bio_read_support_defaults_root_types_and_older_error_paths() {
        use std::error::Error;
        let params: crate::BioReadParams = Default::default();
        assert_eq!(params.format, BioCoordinateFormat::Unknown);
        assert_eq!(params.source_name, "<string>");
        let _: Option<crate::BioReadError> = None;
        let structure_error = crate::BioReadError::WrongFormat {
            source_name: "named.cif".into(),
            format: BioCoordinateFormat::ChemComp,
        };
        let projected = ProteinReadError::Structure(structure_error);
        assert!(projected.to_string().contains("named.cif"));
        assert!(
            projected
                .source()
                .unwrap()
                .downcast_ref::<crate::BioReadError>()
                .is_some()
        );
        assert!(matches!(
            Protein::from_pdb("ATOM\n"),
            Err(ProteinReadError::Pdb(_))
        ));
        assert!(matches!(
            Protein::from_mmcif("invalid cif"),
            Err(ProteinReadError::Mmcif(_))
        ));
    }

    #[test]
    fn bio_ops_both_objects_share_impl_and_only_copy_writable_blocks() {
        let source = BioStructure::from_pdb(ATOM).unwrap();
        let snapshot = source.clone();
        assert!(Arc::ptr_eq(
            &source.data.coordinates,
            &snapshot.data.coordinates
        ));
        let transformed = source.with_translated_coordinates([2.0, 3.0, 4.0]).unwrap();
        let mut inplace = source.clone();
        inplace.translate_([2.0, 3.0, 4.0]).unwrap();
        assert_eq!(transformed, inplace);
        assert_eq!(source, snapshot);
        assert_eq!(source.coordinates().positions(), &[[1.0, 2.0, 3.0]]);
        assert_eq!(transformed.coordinates().positions(), &[[3.0, 5.0, 7.0]]);
        assert!(Arc::ptr_eq(&source.data.atoms, &transformed.data.atoms));
        assert!(Arc::ptr_eq(
            &source.data.residues,
            &transformed.data.residues
        ));
        assert!(Arc::ptr_eq(
            &source.data.metadata,
            &transformed.data.metadata
        ));
        assert!(!Arc::ptr_eq(
            &source.data.coordinates,
            &transformed.data.coordinates
        ));
        let source = source.protein().unwrap();
        let snapshot = source.clone();
        let transformed = source.with_translated_coordinates([2.0, 3.0, 4.0]).unwrap();
        let mut inplace = source.clone();
        inplace.translate_([2.0, 3.0, 4.0]).unwrap();
        assert_eq!(inplace, transformed);
        assert_eq!(source, snapshot);
        assert_eq!(
            transformed.atoms().next().unwrap().position(),
            [3.0, 5.0, 7.0]
        );
        assert!(Arc::ptr_eq(
            &source.structure.data.atoms,
            &transformed.structure.data.atoms
        ));
        assert!(!Arc::ptr_eq(
            &source.structure.data.coordinates,
            &transformed.structure.data.coordinates
        ));
        assert_eq!(BIO_STRUCTURE_OPS.len(), 4);
        for spec in BIO_STRUCTURE_OPS
            .iter()
            .filter(|spec| spec.id == "translate_coordinates")
        {
            assert_eq!(spec.write, ["coordinates"]);
            assert!(spec.read.is_empty());
            for method in [spec.value_method, spec.inplace_method] {
                let id = format!("{}.{}", spec.target, method);
                let row = BINDING_CONTRACT
                    .iter()
                    .find(|row| row.semantic_id == id)
                    .unwrap();
                assert_eq!(row.callable.unwrap().operation_semantic_id, Some(method));
                assert_eq!(
                    row.callable.unwrap().state_model,
                    if method == spec.value_method {
                        StateModel::ValueReturning
                    } else {
                        StateModel::InPlace
                    }
                );
            }
        }
    }

    #[test]
    fn bio_selection_public_both_types_value_inplace_metadata_and_empty_parents() {
        let pdb = concat!(
            "MODEL        1\n",
            "ATOM     41  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n",
            "ENDMDL\nMODEL        2\n",
            "ATOM     42  N   GLY B   2       4.000   5.000   6.000  1.00 20.00           N  \n",
            "ENDMDL\nEND\n"
        );
        let mut source = BioStructure::from_pdb(pdb).unwrap();
        Arc::make_mut(&mut source.data.assemblies).push(BioAssembly::new(
            "assembly".into(),
            true,
            false,
            BioAssemblySpecialKind::NotApplicable,
            2,
            "dimer".into(),
            "source".into(),
            f64::NAN,
            -0.0,
            f64::MIN_POSITIVE,
            vec![BioAssemblyGenerator::new(
                vec!["A".into(), "NEVER_SELECTED".into()],
                vec![],
                vec![],
            )],
        ));
        let state = Arc::make_mut(&mut source.data.source_state);
        state.name = "source".into();
        state.conect_map.insert(41, vec![42]);
        state.has_d_fraction = true;
        state.non_ascii_line = 3;
        state.ter_status = b'Z';
        state.raw_remarks = vec!["duplicate".into(), "duplicate".into()];
        state.resolution = -0.0;
        let peer = source.clone();
        let positions = source.data.coordinates.positions().to_vec();
        // Literal source hierarchy: two models, each with one chain/residue/atom.
        // Gemmi retains selected empty parents, even when no atoms survive.
        let profiles = [
            ("/", [2, 2, 2, 2], vec![0, 1]),
            ("/99", [0, 0, 0, 0], vec![]),
            ("//A", [2, 1, 1, 1], vec![0]),
            ("//Z", [2, 0, 0, 0], vec![]),
            ("//*/(ALA)", [2, 2, 1, 1], vec![0]),
            ("//*//CA", [2, 2, 2, 1], vec![0]),
            ("//*//ZZ", [2, 2, 2, 0], vec![]),
        ];
        for protein_target in [false, true] {
            for (cid, counts, indices) in &profiles {
                let selection = crate::BioSelection::from_cid(cid).unwrap();
                let protein = Protein {
                    structure: source.clone(),
                };
                let output = if protein_target {
                    let value = protein.with_selection(&selection).unwrap();
                    let mut inplace = protein.clone();
                    inplace.retain_selection_(&selection).unwrap();
                    assert_eq!(inplace.num_atoms(), value.num_atoms());
                    assert_eq!(
                        inplace.as_bio_structure().coordinates(),
                        value.as_bio_structure().coordinates()
                    );
                    value.into_bio_structure()
                } else {
                    let value = source.with_selection(&selection).unwrap();
                    let mut inplace = source.clone();
                    inplace.retain_selection_(&selection).unwrap();
                    assert_eq!(inplace.models(), value.models());
                    assert_eq!(inplace.chains(), value.chains());
                    assert_eq!(inplace.residues(), value.residues());
                    assert_eq!(inplace.atoms(), value.atoms());
                    assert_eq!(inplace.coordinates(), value.coordinates());
                    value
                };
                assert_eq!(
                    [
                        output.num_models(),
                        output.num_chains(),
                        output.num_residues(),
                        output.num_atoms()
                    ],
                    *counts,
                    "{cid}"
                );
                assert_eq!(output.input_format(), source.input_format());
                for (i, original) in indices.iter().enumerate() {
                    assert_eq!(
                        output.data.coordinates.positions()[i].map(f64::to_bits),
                        positions[*original].map(f64::to_bits)
                    );
                    assert_eq!(
                        output.data.atoms[i].source(),
                        source.data.atoms[*original].source()
                    );
                }
                assert!(Arc::ptr_eq(&output.data.entities, &source.data.entities));
                assert!(Arc::ptr_eq(
                    &output.data.connections,
                    &source.data.connections
                ));
                assert!(Arc::ptr_eq(&output.data.cispeps, &source.data.cispeps));
                assert!(Arc::ptr_eq(
                    &output.data.mod_residues,
                    &source.data.mod_residues
                ));
                assert!(Arc::ptr_eq(&output.data.helices, &source.data.helices));
                assert!(Arc::ptr_eq(&output.data.sheets, &source.data.sheets));
                assert!(Arc::ptr_eq(&output.data.metadata, &source.data.metadata));
                assert!(Arc::ptr_eq(&output.data.crystal, &source.data.crystal));
                assert!(Arc::ptr_eq(
                    &output.data.ncs_operators,
                    &source.data.ncs_operators
                ));
                assert!(Arc::ptr_eq(
                    &output.data.assemblies,
                    &source.data.assemblies
                ));
                assert_eq!(
                    output.data.assemblies[0].generators[0].chains,
                    ["A", "NEVER_SELECTED"]
                );
                assert_eq!(
                    output.data.source_state.raw_remarks,
                    ["duplicate", "duplicate"]
                );
                assert_eq!(
                    output.data.source_state.resolution.to_bits(),
                    (-0.0f64).to_bits()
                );
                assert!(output.data.source_state.conect_map.is_empty());
                assert!(!output.data.source_state.has_d_fraction);
                assert_eq!(output.data.source_state.non_ascii_line, 0);
                assert_eq!(output.data.source_state.ter_status, 0);
                output.validate().unwrap();
                assert!(Arc::ptr_eq(&source.data.atoms, &peer.data.atoms));
                assert!(Arc::ptr_eq(
                    &source.data.coordinates,
                    &peer.data.coordinates
                ));
                assert_eq!(source.data.coordinates.positions(), positions);
                assert_eq!(
                    source.data.source_state.conect_map.get(&41),
                    Some(&vec![42])
                );
            }
        }
    }

    #[test]
    fn bio_selection_public_typed_errors_rollback_both_receivers() {
        use std::error::Error;
        let mut source = BioStructure::from_pdb(ATOM).unwrap();
        Arc::make_mut(&mut source.data.models)[0] =
            BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), None);
        let peer = source.clone();
        let selection = crate::BioSelection::from_cid("/").unwrap();
        for inplace in [false, true] {
            let mut value = source.clone();
            let error = if inplace {
                value.retain_selection_(&selection).unwrap_err()
            } else {
                value.with_selection(&selection).unwrap_err()
            };
            let BioOperationError::Selection(copy) = &error else {
                panic!("{error}")
            };
            assert!(matches!(
                copy.cause(),
                crate::BioSelectionCopyCause::Traverse(crate::BioRowTraverseError::Model(
                    crate::BioRowModelError::MissingModelNumber
                ))
            ));
            assert!(
                error
                    .source()
                    .unwrap()
                    .downcast_ref::<crate::BioSelectionCopyError>()
                    .is_some()
            );
            assert!(
                copy.source()
                    .unwrap()
                    .source()
                    .unwrap()
                    .downcast_ref::<crate::BioRowModelError>()
                    .is_some()
            );
            assert_eq!(value, source);
            assert!(Arc::ptr_eq(&value.data.atoms, &peer.data.atoms));
            let mut protein = Protein {
                structure: source.clone(),
            };
            let error = if inplace {
                protein.retain_selection_(&selection).unwrap_err()
            } else {
                protein.with_selection(&selection).unwrap_err()
            };
            assert!(matches!(error, BioOperationError::Selection(_)));
            assert_eq!(protein.structure, source);
            assert!(Arc::ptr_eq(
                &protein.structure.data.coordinates,
                &peer.data.coordinates
            ));
        }
    }

    #[test]
    fn bio_ops_failure_validation_and_unwind_are_atomic_on_both_targets() {
        let source = BioStructure::from_pdb(ATOM).unwrap();
        let protein = source.protein().unwrap();
        assert_eq!(source.probe_failure(0).unwrap(), source);
        assert_eq!(protein.probe_failure(0).unwrap(), protein);
        for mode in [1, 2] {
            let mut value = source.clone();
            assert!(value.probe_failure_(mode).is_err());
            assert_eq!(value, source);
            assert!(Arc::ptr_eq(
                &value.data.coordinates,
                &source.data.coordinates
            ));
            let mut value = protein.clone();
            assert!(value.probe_failure_(mode).is_err());
            assert_eq!(value, protein);
        }
        let mut value = source.clone();
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| value.probe_failure_(3)))
                .is_err()
        );
        assert_eq!(value, source);
        let mut value = protein.clone();
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| value.probe_failure_(3)))
                .is_err()
        );
        assert_eq!(value, protein);
    }

    #[test]
    fn bio_ops_protein_cannot_commit_non_amino_acid_rows() {
        let source = BioStructure::from_pdb(ATOM).unwrap().protein().unwrap();
        let water = BioStructure::from_pdb(&ATOM.replace("ALA", "HOH")).unwrap();
        let mut value = source.clone();
        let error = value
            .probe_residues_(water.residues().to_vec())
            .unwrap_err();
        assert!(matches!(
            error,
            BioOperationError::Protein(ProteinProjectionError::NonAminoAcidResidue { index: 0 })
        ));
        assert_eq!(value, source);
    }

    #[test]
    fn protein_recovery_storage_single_wrapper_preserves_projection_and_cow() {
        let mixed = format!(
            "{ATOM}HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 20.00           O  \n"
        );
        let source = BioStructure::from_pdb(&mixed).unwrap();
        let protein = source.protein().unwrap();
        assert_eq!(source.atoms().len(), 2);
        assert_eq!(protein.atoms().len(), 1);
        assert_eq!(protein.structure.atoms().len(), 1);
        let snapshot = protein.clone();
        assert!(Arc::ptr_eq(
            &protein.structure.data.atoms,
            &snapshot.structure.data.atoms
        ));
        let translated = protein
            .with_translated_coordinates([1.0, 0.0, 0.0])
            .unwrap();
        let mut inplace = protein.clone();
        inplace.translate_([1.0, 0.0, 0.0]).unwrap();
        assert_eq!(translated, inplace);
        assert_eq!(protein, snapshot);
        assert_eq!(
            translated.atoms().next().unwrap().position(),
            [2.0, 2.0, 3.0]
        );
        assert!(Arc::ptr_eq(
            &protein.structure.data.atoms,
            &translated.structure.data.atoms
        ));
        assert!(!Arc::ptr_eq(
            &protein.structure.data.coordinates,
            &translated.structure.data.coordinates
        ));
    }

    #[test]
    fn protein_recovery_storage_failed_finalization_and_unwind_preserve_receiver() {
        let source = BioStructure::from_pdb(ATOM).unwrap().protein().unwrap();
        let snapshot = source.clone();
        let water = BioStructure::from_pdb(&ATOM.replace("ALA", "HOH")).unwrap();
        let mut candidate = source.clone();
        assert!(matches!(
            candidate.probe_residues_(water.residues().to_vec()),
            Err(BioOperationError::Protein(
                ProteinProjectionError::NonAminoAcidResidue { index: 0 }
            ))
        ));
        assert_eq!(candidate, snapshot);
        assert!(Arc::ptr_eq(
            &candidate.structure.data.atoms,
            &snapshot.structure.data.atoms
        ));
        let mut candidate = source.clone();
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                candidate.probe_failure_(3)
            }))
            .is_err()
        );
        assert_eq!(candidate, snapshot);
    }
}

#[cfg(test)]
mod bio_pdb_output_storage_tests {
    use super::BioStructure;
    use crate::BioPdbWriteParams;
    use cosmolkit_bio::BioMetadata;
    use cosmolkit_bio::BioStructureSourceState;
    use std::sync::Arc;

    const PDB_TEXT: &str = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nTER       2      ALA A   1                                                      \nEND";

    /// Coordinate float-bit snapshot from both receivers.
    struct BitsSnapshot {
        positions: Vec<[u64; 3]>,
        source_state_name_len: usize,
        models_len: usize,
        chains_len: usize,
        residues_len: usize,
        atoms_len: usize,
        entities_len: usize,
        connections_len: usize,
        metadata_is_default: bool,
        source_state_is_default: bool,
        crystal_is_none: bool,
        ncs_len: usize,
        assemblies_len: usize,
    }

    fn capture(structure: &BioStructure) -> BitsSnapshot {
        let parts = structure.clone().into_parts();
        BitsSnapshot {
            positions: parts
                .coordinates
                .positions()
                .iter()
                .map(|p| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()])
                .collect(),
            source_state_name_len: parts.source_state.name.len(),
            models_len: parts.models.len(),
            chains_len: parts.chains.len(),
            residues_len: parts.residues.len(),
            atoms_len: parts.atoms.len(),
            entities_len: parts.entities.len(),
            connections_len: parts.connections.len(),
            metadata_is_default: parts.metadata == BioMetadata::default(),
            source_state_is_default: parts.source_state == BioStructureSourceState::default(),
            crystal_is_none: parts.crystal.is_none(),
            ncs_len: parts.ncs_operators.len(),
            assemblies_len: parts.assemblies.len(),
        }
    }

    fn compare(
        label: &str,
        before: &BitsSnapshot,
        after: &BitsSnapshot,
        discrepancies: &mut Vec<String>,
    ) {
        if before.positions != after.positions {
            discrepancies.push(format!("{label}: coordinate float bits mutated"));
        }
        if before.models_len != after.models_len
            || before.chains_len != after.chains_len
            || before.residues_len != after.residues_len
            || before.atoms_len != after.atoms_len
        {
            discrepancies.push(format!("{label}: hierarchy lengths mutated"));
        }
        if before.entities_len != after.entities_len
            || before.connections_len != after.connections_len
        {
            discrepancies.push(format!("{label}: entity/connection lengths mutated"));
        }
        if before.metadata_is_default != after.metadata_is_default
            || before.source_state_is_default != after.source_state_is_default
        {
            discrepancies.push(format!("{label}: metadata/source_state mutated"));
        }
        if before.crystal_is_none != after.crystal_is_none
            || before.ncs_len != after.ncs_len
            || before.assemblies_len != after.assemblies_len
        {
            discrepancies.push(format!("{label}: crystal/ncs/assemblies mutated"));
        }
        if before.source_state_name_len != after.source_state_name_len {
            discrepancies.push(format!("{label}: source_state name mutated"));
        }
    }

    /// Arc sharing identity: all 16 Arc fields of BioStructureData must
    /// be pointer-equal between the original and its clone (proving the
    /// clone shares, not deep-copies), AND unchanged (same pointers)
    /// after each write call.
    fn arc_identities(structure: &BioStructure) -> [usize; 16] {
        let d = &structure.data;
        [
            Arc::as_ptr(&d.models) as usize,
            Arc::as_ptr(&d.chains) as usize,
            Arc::as_ptr(&d.residues) as usize,
            Arc::as_ptr(&d.atoms) as usize,
            Arc::as_ptr(&d.entities) as usize,
            Arc::as_ptr(&d.connections) as usize,
            Arc::as_ptr(&d.cispeps) as usize,
            Arc::as_ptr(&d.mod_residues) as usize,
            Arc::as_ptr(&d.helices) as usize,
            Arc::as_ptr(&d.sheets) as usize,
            Arc::as_ptr(&d.metadata) as usize,
            Arc::as_ptr(&d.source_state) as usize,
            Arc::as_ptr(&d.coordinates) as usize,
            Arc::as_ptr(&d.crystal) as usize,
            Arc::as_ptr(&d.ncs_operators) as usize,
            Arc::as_ptr(&d.assemblies) as usize,
        ]
    }

    /// Owned COMPLETE payload baseline from into_parts(): every field
    /// including input_format, not just lengths or default flags.
    fn payload_baseline(structure: &BioStructure) -> cosmolkit_bio::BioStructureParts {
        structure.clone().into_parts()
    }

    fn compare_payloads(
        label: &str,
        before: &cosmolkit_bio::BioStructureParts,
        after: &cosmolkit_bio::BioStructureParts,
        discrepancies: &mut Vec<String>,
    ) {
        if before.input_format != after.input_format {
            discrepancies.push(format!("{label}: input_format changed"));
        }
        if before.models != after.models
            || before.chains != after.chains
            || before.residues != after.residues
            || before.atoms != after.atoms
        {
            discrepancies.push(format!("{label}: hierarchy rows changed"));
        }
        if before.coordinates != after.coordinates {
            discrepancies.push(format!("{label}: coordinate block changed"));
        }
        if before.metadata != after.metadata || before.source_state != after.source_state {
            discrepancies.push(format!("{label}: metadata/source_state changed"));
        }
        if before.entities != after.entities
            || before.connections != after.connections
            || before.cispeps != after.cispeps
            || before.mod_residues != after.mod_residues
            || before.helices != after.helices
            || before.sheets != after.sheets
        {
            discrepancies.push(format!("{label}: typed metadata rows changed"));
        }
        if before.crystal != after.crystal
            || before.ncs_operators != after.ncs_operators
            || before.assemblies != after.assemblies
        {
            discrepancies.push(format!("{label}: crystal/ncs/assemblies changed"));
        }
    }

    /// Storage proof: to_pdb()/to_pdb_with_params() MUST NOT mutate any
    /// BioStructure internal state. 12 calls across methods/profiles with
    /// PRE-paired checkpoints on BOTH original and cloned peer, coordinate
    /// float bits, ALL shared-block dimensions, and Err inclusion.
    #[test]
    fn bio_pdb_output_storage_unchanged_across_12_calls() {
        let original = BioStructure::from_pdb(PDB_TEXT).expect("parses");
        let peer = original.clone();

        let params_set = [
            BioPdbWriteParams::default(),
            BioPdbWriteParams {
                end_record: false,
                ..Default::default()
            },
            BioPdbWriteParams {
                ter_records: false,
                ..Default::default()
            },
        ];

        let mut calls = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        let dir = std::env::temp_dir().join("ck_pdb_storage_test");
        std::fs::create_dir_all(&dir).unwrap();

        // 12 calls: to_pdb ×3 profiles + to_pdb_with_params ×3 profiles
        // + write_pdb_with_params ×3 profiles + write_pdb ×1 (default)
        // + write_pdb_with_params to directory (Err) ×1 + peer write ×1.

        // Calls 1-3: to_pdb default ×3 (idempotent repeated reads)
        for rep in 0..3 {
            let label = format!("to_pdb[{rep}]");
            let orig_before = capture(&original);
            let peer_before = capture(&peer);
            let result = original.to_pdb();
            calls += 1;
            let orig_after = capture(&original);
            let peer_after = capture(&peer);
            compare(
                &format!("{label}/orig"),
                &orig_before,
                &orig_after,
                &mut discrepancies,
            );
            compare(
                &format!("{label}/peer"),
                &peer_before,
                &peer_after,
                &mut discrepancies,
            );
            if let Err(e) = &result {
                discrepancies.push(format!("{label}: unexpected Err {e:?}"));
            }
        }

        // Calls 4-6: to_pdb_with_params ×3 profiles
        for (i, params) in params_set.iter().enumerate() {
            let label = format!("to_pdb_with_params[{i}]");
            let orig_before = capture(&original);
            let peer_before = capture(&peer);
            let result = original.to_pdb_with_params(params);
            calls += 1;
            let orig_after = capture(&original);
            let peer_after = capture(&peer);
            compare(
                &format!("{label}/orig"),
                &orig_before,
                &orig_after,
                &mut discrepancies,
            );
            compare(
                &format!("{label}/peer"),
                &peer_before,
                &peer_after,
                &mut discrepancies,
            );
            if let Err(e) = &result {
                discrepancies.push(format!("{label}: unexpected Err {e:?}"));
            }
        }

        // Calls 7-9: write_pdb_with_params ×3 profiles (Ok results)
        for (i, params) in params_set.iter().enumerate() {
            let label = format!("write_pdb_with_params[{i}]");
            let path = dir.join(format!("storage_{i}.pdb"));
            let orig_before = capture(&original);
            let peer_before = capture(&peer);
            let result = original.write_pdb_with_params(&path, params);
            calls += 1;
            let orig_after = capture(&original);
            let peer_after = capture(&peer);
            compare(
                &format!("{label}/orig"),
                &orig_before,
                &orig_after,
                &mut discrepancies,
            );
            compare(
                &format!("{label}/peer"),
                &peer_before,
                &peer_after,
                &mut discrepancies,
            );
            if let Err(e) = &result {
                discrepancies.push(format!("{label}: unexpected Err {e:?}"));
            }
            let _ = std::fs::remove_file(&path);
        }

        // Call 10: write_pdb default (short form)
        {
            let label = "write_pdb";
            let path = dir.join("storage_default.pdb");
            let orig_before = capture(&original);
            let peer_before = capture(&peer);
            let result = original.write_pdb(&path);
            calls += 1;
            let orig_after = capture(&original);
            let peer_after = capture(&peer);
            compare(
                &format!("{label}/orig"),
                &orig_before,
                &orig_after,
                &mut discrepancies,
            );
            compare(
                &format!("{label}/peer"),
                &peer_before,
                &peer_after,
                &mut discrepancies,
            );
            if let Err(e) = &result {
                discrepancies.push(format!("{label}: unexpected Err {e:?}"));
            }
            let _ = std::fs::remove_file(&path);
        }

        // Call 11: write_pdb to a DIRECTORY (Err result — storage must
        // still be unchanged even on error).
        {
            let label = "write_pdb_to_dir";
            let dir_target = dir.join("target_dir");
            std::fs::create_dir_all(&dir_target).unwrap();
            let orig_before = capture(&original);
            let peer_before = capture(&peer);
            let result = original.write_pdb(&dir_target);
            calls += 1;
            let orig_after = capture(&original);
            let peer_after = capture(&peer);
            compare(
                &format!("{label}/orig"),
                &orig_before,
                &orig_after,
                &mut discrepancies,
            );
            compare(
                &format!("{label}/peer"),
                &peer_before,
                &peer_after,
                &mut discrepancies,
            );
            // Err is EXPECTED here (writing to a directory).
            if result.is_ok() {
                discrepancies.push(format!("{label}: unexpected Ok"));
            }
            let _ = std::fs::remove_dir(&dir_target);
        }

        // Call 12: peer.write_pdb (verify the PEER also stays unchanged)
        {
            let label = "peer_write_pdb";
            let path = dir.join("peer_storage.pdb");
            let orig_before = capture(&original);
            let peer_before = capture(&peer);
            let result = peer.write_pdb(&path);
            calls += 1;
            let orig_after = capture(&original);
            let peer_after = capture(&peer);
            compare(
                &format!("{label}/orig"),
                &orig_before,
                &orig_after,
                &mut discrepancies,
            );
            compare(
                &format!("{label}/peer"),
                &peer_before,
                &peer_after,
                &mut discrepancies,
            );
            if let Err(e) = &result {
                discrepancies.push(format!("{label}: unexpected Err {e:?}"));
            }
            let _ = std::fs::remove_file(&path);
        }

        assert_eq!(calls, 12, "exact 12 storage-proof calls");
        assert!(
            discrepancies.is_empty(),
            "storage discrepancies: {discrepancies:?}"
        );
    }

    /// C00-C11 receiver census: 12 additional storage-proof calls using
    /// DIFFERENT structures (fixture-indexed), with Arc identity, complete
    /// payload baselines, and paired sharing checks.
    #[test]
    fn bio_pdb_output_storage_c00_c11_receivers() {
        // Build simple structures for the 12-receiver census. Each uses
        // a slightly different PDB text to ensure independent state.
        let receivers: Vec<BioStructure> = (0..12)
            .map(|i| {
                let text = format!(
                    "ATOM      1  CA  ALA A   {idx:>4}       1.000   2.000   3.000  1.00 20.00           C  \nTER       2      ALA A   {idx:>4}                                                      \nEND",
                    idx = i + 1
                );
                BioStructure::from_pdb(&text).expect("receiver parses")
            })
            .collect();

        let mut calls = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for (index, structure) in receivers.iter().enumerate() {
            let label = format!("receiver[{index}]");
            let peer = structure.clone();

            // Arc sharing: original and clone share ALL 16 Arc pointers.
            let orig_arcs_before = arc_identities(structure);
            let peer_arcs = arc_identities(&peer);
            if orig_arcs_before != peer_arcs {
                discrepancies.push(format!("{label}: clone does NOT share Arc pointers"));
            }

            // Complete payload baseline.
            let orig_payload_before = payload_baseline(structure);
            let peer_payload_before = payload_baseline(&peer);

            // input_format verification.
            if structure.data.input_format != cosmolkit_bio::BioCoordinateFormat::Pdb {
                discrepancies.push(format!("{label}: unexpected input_format"));
            }

            // One call per receiver.
            let result = structure.to_pdb();
            calls += 1;

            // AFTER: Arc pointers unchanged.
            let orig_arcs_after = arc_identities(structure);
            if orig_arcs_before != orig_arcs_after {
                discrepancies.push(format!("{label}: Arc pointers changed"));
            }

            // AFTER: complete payload unchanged.
            let orig_payload_after = payload_baseline(structure);
            compare_payloads(
                &format!("{label}/orig"),
                &orig_payload_before,
                &orig_payload_after,
                &mut discrepancies,
            );
            let peer_payload_after = payload_baseline(&peer);
            compare_payloads(
                &format!("{label}/peer"),
                &peer_payload_before,
                &peer_payload_after,
                &mut discrepancies,
            );

            // AFTER: peer Arc still shared with original.
            let peer_arcs_after = arc_identities(&peer);
            if orig_arcs_after != peer_arcs_after {
                discrepancies.push(format!("{label}: peer Arc sharing broken"));
            }

            if let Err(e) = &result {
                discrepancies.push(format!("{label}: unexpected Err {e:?}"));
            }
        }

        assert_eq!(calls, 12, "exact 12 C00-C11 receiver calls");
        assert!(
            discrepancies.is_empty(),
            "receiver discrepancies: {discrepancies:?}"
        );
    }
}
