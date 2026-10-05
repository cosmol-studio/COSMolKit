//! Amino-acid-only projection and read-only hierarchy views.

use std::fmt;

use cosmolkit_types::Element;

use crate::{
    AltLocLabel, AtomName, BioAtomId, BioAtomRow, BioChainId, BioChainRow, BioCoordinateBlock,
    BioEntityId, BioModelId, BioModelRow, BioResidueId, BioResidueRow, BioRowSpan,
    BioStructureData, BioStructureError, BioStructureParts, ChainKind, ChainSourceIds, ResidueCode,
    ResidueInfo, ResidueInfoKind, ResidueKind, ResidueName, find_residue_info,
};

#[derive(Debug, Clone, PartialEq)]
pub struct ProteinData {
    structure: BioStructureData,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProteinSelectionSummary {
    pub chains: usize,
    pub residues: usize,
    pub atoms: usize,
}

#[must_use]
pub fn protein_selection_summary(structure: &BioStructureData) -> ProteinSelectionSummary {
    ProteinSelectionSummary {
        chains: structure.chains().len(),
        residues: structure.residues().len(),
        atoms: structure.atoms().len(),
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinProjectionError {
    Structure(BioStructureError),
    MissingEntityMapping { source: BioEntityId },
    NonAminoAcidResidue { index: usize },
}

impl fmt::Display for ProteinProjectionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonAminoAcidResidue { index } => write!(
                formatter,
                "protein contains non-amino-acid residue at row {index}"
            ),
            Self::Structure(error) => write!(formatter, "protein projection failed: {error}"),
            Self::MissingEntityMapping { source } => write!(
                formatter,
                "protein projection lost source entity row {}",
                source.value()
            ),
        }
    }
}

impl std::error::Error for ProteinProjectionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structure(error) => Some(error),
            Self::MissingEntityMapping { .. } | Self::NonAminoAcidResidue { .. } => None,
        }
    }
}

impl From<BioStructureError> for ProteinProjectionError {
    fn from(value: BioStructureError) -> Self {
        Self::Structure(value)
    }
}

/// Source-ordered borrowed protein chains without an intermediate collection.
#[derive(Debug, Clone)]
pub struct ProteinChainIter<'a> {
    structure: &'a BioStructureData,
    cursor: std::ops::Range<usize>,
}

impl<'a> ProteinChainIter<'a> {
    fn new(structure: &'a BioStructureData, cursor: std::ops::Range<usize>) -> Self {
        Self { structure, cursor }
    }
}

impl<'a> Iterator for ProteinChainIter<'a> {
    type Item = ProteinChainRef<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        // Legacy: (0..self.structure.chains().len()).map(|index| ProteinChainRef {
        // Legacy:     protein: self,
        // Legacy:     chain_id: ChainId::new(index as u32),
        // Legacy: })
        // Behavior: one borrowed row reference per source-order cursor step.
        // Complexity: constant-time next, no materialized Vec or heap allocation.
        self.cursor.next().map(|index| ProteinChainRef {
            structure: self.structure,
            chain_id: BioChainId::new(index as u32),
        })
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        // The legacy Range::map iterator forwards nth to Range::nth.
        // Preserve that O(1) skip and exhaustion behavior in this borrowed view.
        self.cursor.nth(n).map(|index| Self::Item {
            structure: self.structure,
            chain_id: BioChainId::new(index as u32),
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.cursor.size_hint()
    }
}

impl ExactSizeIterator for ProteinChainIter<'_> {}
impl std::iter::FusedIterator for ProteinChainIter<'_> {}

/// Borrowed source-ordered protein residues over one contiguous hierarchy span.
#[derive(Debug, Clone)]
pub struct ProteinResidueIter<'a> {
    structure: &'a BioStructureData,
    cursor: std::ops::Range<usize>,
}

impl<'a> ProteinResidueIter<'a> {
    fn new(structure: &'a BioStructureData, cursor: std::ops::Range<usize>) -> Self {
        Self { structure, cursor }
    }
}

impl<'a> Iterator for ProteinResidueIter<'a> {
    type Item = ProteinResidueRef<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        // Legacy: (span.start as usize..span.end() as usize).map(move |index| ProteinResidueRef {
        // Legacy:     protein,
        // Legacy:     residue_id: ResidueId::new(index as u32),
        // Legacy: })
        // Behavior: the validated chain span and whole-structure range share
        // the same source-order traversal without collecting skipped rows.
        // Complexity: constant time per row and no per-next allocations.
        self.cursor.next().map(|index| ProteinResidueRef {
            structure: self.structure,
            residue_id: BioResidueId::new(index as u32),
        })
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        // The legacy Range::map iterator forwards nth to Range::nth.
        // Preserve that O(1) skip and exhaustion behavior in this borrowed view.
        self.cursor.nth(n).map(|index| Self::Item {
            structure: self.structure,
            residue_id: BioResidueId::new(index as u32),
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.cursor.size_hint()
    }
}

impl ExactSizeIterator for ProteinResidueIter<'_> {}
impl std::iter::FusedIterator for ProteinResidueIter<'_> {}

/// Borrowed atom traversal over a validated contiguous residue or chain span.
#[derive(Debug, Clone)]
pub struct ProteinAtomIter<'a> {
    structure: &'a BioStructureData,
    cursor: std::ops::Range<usize>,
}

impl<'a> ProteinAtomIter<'a> {
    fn new(structure: &'a BioStructureData, cursor: std::ops::Range<usize>) -> Self {
        Self { structure, cursor }
    }

    fn from_residues(structure: &'a BioStructureData, residues: std::ops::Range<usize>) -> Self {
        let start = structure
            .residues()
            .get(residues.start)
            .map_or(0, |row| row.atom_span().start() as usize);
        let end = if residues.is_empty() {
            start
        } else {
            structure.residues()[residues.end - 1].atom_span().end() as usize
        };
        Self::new(structure, start..end)
    }
}

impl<'a> Iterator for ProteinAtomIter<'a> {
    type Item = ProteinAtomRef<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        // Legacy: (span.start as usize..span.end() as usize).map(move |index| ProteinAtomRef {
        // Legacy:     protein,
        // Legacy:     atom_id: AtomId::new(index as u32),
        // Legacy: })
        // Behavior: validated contiguous atom spans preserve source order,
        // including zero-width residue spans between populated residues.
        // Complexity: O(1) per row, no per-residue Vec or heap allocation.
        self.cursor.next().map(|index| ProteinAtomRef {
            structure: self.structure,
            atom_id: BioAtomId::new(index as u32),
        })
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        // The legacy Range::map iterator forwards nth to Range::nth.
        // Preserve that O(1) skip and exhaustion behavior in this borrowed view.
        self.cursor.nth(n).map(|index| Self::Item {
            structure: self.structure,
            atom_id: BioAtomId::new(index as u32),
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.cursor.size_hint()
    }
}

impl ExactSizeIterator for ProteinAtomIter<'_> {}
impl std::iter::FusedIterator for ProteinAtomIter<'_> {}

#[cfg(test)]
mod bio_legacy_iterator_tests {
    use super::{BioStructureData, ProteinAtomIter, ProteinChainIter, ProteinResidueIter};
    use crate::{
        AtomName, AtomSourceIds, BioAtomRow, BioCalcFlag, BioChainId, BioChainRow,
        BioCoordinateBlock, BioCoordinateFormat, BioModelId, BioModelRow, BioResidueId,
        BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructureParts, ChainKind,
        ChainSourceIds, EntityKind, ResidueInfoKind, ResidueName, ResidueSourceIds,
    };
    use cosmolkit_types::Element;

    fn structure(chain_count: usize) -> BioStructureData {
        let span = |start, len| BioRowSpan::new(start, len).unwrap();
        BioStructureData::from_parts(BioStructureParts {
            input_format: BioCoordinateFormat::Unknown,
            models: if chain_count == 0 {
                vec![]
            } else {
                vec![BioModelRow::new(span(0, chain_count as u32), Some(1))]
            },
            chains: (0..chain_count)
                .map(|_| {
                    BioChainRow::new(
                        BioModelId::new(0),
                        None,
                        BioRowSpan::new(0, 0).unwrap(),
                        ChainKind::Protein,
                        ChainSourceIds::default(),
                    )
                })
                .collect(),
            residues: vec![],
            atoms: vec![],
            entities: vec![],
            connections: vec![],
            cispeps: vec![],
            mod_residues: vec![],
            helices: vec![],
            sheets: vec![],
            metadata: Default::default(),
            source_state: Default::default(),
            coordinates: BioCoordinateBlock::default(),
            crystal: None,
            ncs_operators: vec![],
            assemblies: vec![],
        })
        .unwrap()
    }

    fn with_residues() -> BioStructureData {
        let mut parts = structure(0).into_parts();
        parts.models = vec![
            BioModelRow::new(BioRowSpan::new(0, 2).unwrap(), Some(1)),
            BioModelRow::new(BioRowSpan::new(2, 1).unwrap(), Some(2)),
        ];
        parts.chains = [(0, 0, 2), (0, 2, 0), (1, 2, 1)]
            .into_iter()
            .map(|(model, start, len)| {
                BioChainRow::new(
                    BioModelId::new(model),
                    None,
                    BioRowSpan::new(start, len).unwrap(),
                    ChainKind::Protein,
                    ChainSourceIds::default(),
                )
            })
            .collect();
        parts.residues = [(0, "ALA"), (0, "GLY"), (2, "SER")]
            .into_iter()
            .map(|(chain, name)| {
                BioResidueRow::new(
                    BioChainId::new(chain),
                    BioRowSpan::new(0, 0).unwrap(),
                    ResidueName::from_ascii(name.as_bytes()).unwrap(),
                    ResidueInfoKind::Aa,
                    EntityKind::Polymer,
                    None,
                    None,
                    ResidueSourceIds::default(),
                    BioSiftsUnpResidue::default(),
                )
            })
            .collect();
        BioStructureData::from_parts(parts).unwrap()
    }

    fn with_atoms() -> BioStructureData {
        let mut parts = with_residues().into_parts();
        for (index, (start, len)) in [(0, 1), (1, 0), (1, 1)].into_iter().enumerate() {
            parts.residues[index] = BioResidueRow::new(
                BioChainId::new(if index == 2 { 2 } else { 0 }),
                BioRowSpan::new(start, len).unwrap(),
                ResidueName::from_ascii(if index == 0 {
                    b"ALA"
                } else if index == 1 {
                    b"GLY"
                } else {
                    b"SER"
                })
                .unwrap(),
                ResidueInfoKind::Aa,
                EntityKind::Polymer,
                None,
                None,
                ResidueSourceIds::default(),
                BioSiftsUnpResidue::default(),
            );
        }
        parts.atoms = [(0, b"CA".as_slice()), (2, b"CB".as_slice())]
            .into_iter()
            .map(|(residue, name)| {
                BioAtomRow::new(
                    BioResidueId::new(residue),
                    AtomName::from_ascii(name).unwrap(),
                    Element::C,
                    None,
                    None,
                    0,
                    BioCalcFlag::NotSet,
                    1.0,
                    20.0,
                    [0.0; 6],
                    -1,
                    0.0,
                    AtomSourceIds::default(),
                )
            })
            .collect();
        parts.coordinates = BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        BioStructureData::from_parts(parts).unwrap()
    }

    #[test]
    fn bio_legacy_i03_atom_iterator_empty_interspersed_and_borrowed_positions() {
        let data = with_atoms();
        let mut empty = ProteinAtomIter::from_residues(&data, 1..2);
        assert!(empty.next().is_none());
        assert!(empty.next().is_none());
        let mut first_chain = ProteinAtomIter::from_residues(&data, 0..2);
        assert_eq!(first_chain.next().unwrap().position(), [1.0, 2.0, 3.0]);
        assert!(first_chain.next().is_none());
        let mut all = ProteinAtomIter::new(&data, 0..2);
        assert_eq!(all.next().unwrap().id().value(), 0);
        let last = all.next().unwrap();
        assert_eq!(last.id().value(), 1);
        assert_eq!(last.position(), [4.0, 5.0, 6.0]);
        assert!(all.next().is_none());
        assert!(all.next().is_none());
    }

    #[test]
    fn bio_legacy_i04_six_borrowed_traversals_are_source_ordered_and_interleavable() {
        let data = with_atoms();
        let protein = data.protein().unwrap();
        let mut chains = protein.chains();
        assert_eq!(chains.next().unwrap().id().value(), 0);
        let chain = protein.chain(0).unwrap();
        let mut chain_residues = chain.residues();
        assert_eq!(chain_residues.next().unwrap().id().value(), 0);
        let residue = protein.residues().next().unwrap();
        assert_eq!(residue.id().value(), 0);
        assert_eq!(residue.atoms().next().unwrap().position(), [1.0, 2.0, 3.0]);
        assert_eq!(
            protein
                .atoms()
                .map(|atom| atom.id().value())
                .collect::<Vec<_>>(),
            [0, 1]
        );
        assert_eq!(chain.atoms().next().unwrap().position(), [1.0, 2.0, 3.0]);
        assert_eq!(chain_residues.next().unwrap().id().value(), 1);
        assert_eq!(chains.next().unwrap().id().value(), 1);
        assert!(chains.next().is_none());
        assert_eq!(
            protein.chain(1).unwrap().atoms().next().unwrap().position(),
            [4.0, 5.0, 6.0]
        );
        assert!(std::ptr::eq(
            residue.row(),
            &protein.structure().residues()[0]
        ));
    }

    #[test]
    fn bio_legacy_i02_residue_iterator_whole_empty_chain_and_model_gap() {
        let data = with_residues();
        let mut all = ProteinResidueIter::new(&data, 0..data.residues().len());
        let mut copy = all.clone();
        assert_eq!(all.next().unwrap().id().value(), 0);
        assert_eq!(copy.next().unwrap().id().value(), 0);
        assert_eq!(
            all.map(|residue| residue.id().value()).collect::<Vec<_>>(),
            [1, 2]
        );
        assert_eq!(
            copy.map(|residue| residue.id().value()).collect::<Vec<_>>(),
            [1, 2]
        );
        let mut empty = ProteinResidueIter::new(&data, 2..2);
        assert!(empty.next().is_none());
        assert!(empty.next().is_none());
        let mut last = ProteinResidueIter::new(&data, 2..3);
        assert_eq!(last.next().unwrap().name().as_str(), "SER");
        assert!(last.next().is_none());
        assert!(last.next().is_none());
    }

    #[test]
    fn bio_legacy_i01_chain_iterator_empty_single_many_and_independent_cursors() {
        for count in [0, 1, 3] {
            let data = structure(count);
            let mut first = ProteinChainIter::new(&data, 0..data.chains().len());
            let mut second = first.clone();
            for index in 0..count {
                assert_eq!(first.next().unwrap().id().value(), index as u32);
                assert_eq!(second.next().unwrap().id().value(), index as u32);
            }
            assert!(first.next().is_none());
            assert!(first.next().is_none());
            assert!(second.next().is_none());
        }
        let data = structure(3);
        let mut partial = ProteinChainIter::new(&data, 0..3);
        assert_eq!(partial.next().unwrap().id().value(), 0);
        assert_eq!(partial.next().unwrap().id().value(), 1);
        assert_eq!(partial.size_hint(), (1, Some(1)));
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ProteinChainRef<'a> {
    structure: &'a BioStructureData,
    chain_id: BioChainId,
}

#[derive(Debug, Clone, Copy)]
pub struct ProteinResidueRef<'a> {
    structure: &'a BioStructureData,
    residue_id: BioResidueId,
}

#[derive(Debug, Clone, Copy)]
pub struct ProteinAtomRef<'a> {
    structure: &'a BioStructureData,
    atom_id: BioAtomId,
}

impl BioStructureData {
    pub fn protein(&self) -> Result<ProteinData, ProteinProjectionError> {
        ProteinData::project(self)
    }
}

impl ProteinData {
    /// Detached editing only; callers must validate before publishing a result.
    pub fn structure_mut(&mut self) -> &mut BioStructureData {
        &mut self.structure
    }

    pub fn structure(&self) -> &BioStructureData {
        &self.structure
    }

    #[must_use]
    pub fn into_structure(self) -> BioStructureData {
        self.structure
    }

    pub fn validate(&self) -> Result<(), ProteinProjectionError> {
        validate_protein_structure(&self.structure)
    }

    fn project(source: &BioStructureData) -> Result<Self, ProteinProjectionError> {
        let retained_residues: Vec<bool> = source
            .residues()
            .iter()
            .map(|residue| is_amino_acid_kind(residue.residue_info_kind()))
            .collect();

        let mut retained_chains = vec![false; source.chains().len()];
        for (chain_index, chain) in source.chains().iter().enumerate() {
            let range = span_range(chain.residue_span(), source.residues())?;
            retained_chains[chain_index] = retained_residues[range].iter().any(|retain| *retain);
        }

        let mut retained_models = vec![false; source.models().len()];
        for (model_index, model) in source.models().iter().enumerate() {
            let range = span_range(model.chain_span(), source.chains())?;
            retained_models[model_index] = retained_chains[range].iter().any(|retain| *retain);
        }

        let mut retained_entities = vec![false; source.entities().len()];
        for (chain_index, chain) in source.chains().iter().enumerate() {
            if !retained_chains[chain_index] {
                continue;
            }
            if let Some(entity_id) = chain.entity_id() {
                retained_entities[entity_id.index()] = true;
            }
            for residue_index in span_range(chain.residue_span(), source.residues())? {
                if !retained_residues[residue_index] {
                    continue;
                }
                if let Some(entity_id) = source.residues()[residue_index].entity_id() {
                    retained_entities[entity_id.index()] = true;
                }
            }
        }

        let mut entity_map = vec![None; source.entities().len()];
        let mut entities = Vec::new();
        for (source_index, (entity, retain)) in source
            .entities()
            .iter()
            .zip(retained_entities.iter())
            .enumerate()
        {
            if *retain {
                let projected_id = BioEntityId::new(checked_row_index(entities.len())?);
                entity_map[source_index] = Some(projected_id);
                entities.push(entity.clone());
            }
        }

        let mut models = Vec::new();
        let mut chains = Vec::new();
        let mut residues = Vec::new();
        let mut atoms = Vec::new();
        let mut positions = Vec::new();

        for (model_index, model) in source.models().iter().enumerate() {
            if !retained_models[model_index] {
                continue;
            }
            let projected_model_id = BioModelId::new(checked_row_index(models.len())?);
            let chain_start = chains.len();

            for chain_index in span_range(model.chain_span(), source.chains())? {
                if !retained_chains[chain_index] {
                    continue;
                }
                let source_chain = &source.chains()[chain_index];
                let projected_chain_id = BioChainId::new(checked_row_index(chains.len())?);
                let residue_start = residues.len();

                for residue_index in span_range(source_chain.residue_span(), source.residues())? {
                    if !retained_residues[residue_index] {
                        continue;
                    }
                    let source_residue = &source.residues()[residue_index];
                    let projected_residue_id =
                        BioResidueId::new(checked_row_index(residues.len())?);
                    let atom_start = atoms.len();

                    for atom_index in span_range(source_residue.atom_span(), source.atoms())? {
                        let source_atom = &source.atoms()[atom_index];
                        atoms.push(BioAtomRow::new(
                            projected_residue_id,
                            source_atom.name(),
                            source_atom.element(),
                            source_atom.isotope_mass_number(),
                            source_atom.altloc(),
                            source_atom.formal_charge(),
                            source_atom.calc_flag(),
                            source_atom.occupancy(),
                            source_atom.b_iso(),
                            *source_atom.anisou(),
                            source_atom.tls_group_id(),
                            source_atom.fraction(),
                            *source_atom.source(),
                        ));
                        positions.push(source.coordinates().positions()[atom_index]);
                    }

                    residues.push(BioResidueRow::new(
                        projected_chain_id,
                        BioRowSpan::from_usize(atom_start, atoms.len() - atom_start)?,
                        source_residue.name(),
                        source_residue.residue_info_kind(),
                        source_residue.entity_kind(),
                        remap_entity(source_residue.entity_id(), &entity_map)?,
                        source_residue.het_flag(),
                        source_residue.source().clone(),
                        source_residue.sifts_unp(),
                    ));
                }

                chains.push(BioChainRow::new(
                    projected_model_id,
                    remap_entity(source_chain.entity_id(), &entity_map)?,
                    BioRowSpan::from_usize(residue_start, residues.len() - residue_start)?,
                    ChainKind::Protein,
                    source_chain.source().clone(),
                ));
            }

            models.push(BioModelRow::new(
                BioRowSpan::from_usize(chain_start, chains.len() - chain_start)?,
                model.source_model_number(),
            ));
        }

        let structure = BioStructureData::from_parts(BioStructureParts {
            input_format: source.input_format(),
            models,
            chains,
            residues,
            atoms,
            entities,
            // Gemmi✔️✔️:     st.name = name;
            // Gemmi✔️✔️:     st.connections = connections;
            // Gemmi✔️✔️:     st.cispeps = cispeps;
            // Gemmi✔️✔️:     st.mod_residues = mod_residues;
            // Gemmi✔️✔️:     st.helices = helices;
            // Gemmi✔️✔️:     st.sheets = sheets;
            // Gemmi✔️✔️:     st.meta = meta;
            // Gemmi✔️✔️:     st.info = info;
            // Gemmi✔️✔️:     st.raw_remarks = raw_remarks;
            // Gemmi✔️✔️:     st.has_origx = has_origx;
            // Gemmi✔️✔️:     st.origx = origx;
            // Gemmi✔️✔️:     st.resolution = resolution;
            // Behavior review: M23 preserves source-address relationships and
            // source metadata verbatim through the projection, even when an
            // address names a row excluded by the protein-only hierarchy.
            // Assemblies remain omitted above under the existing projection
            // rule; source serial connectivity remains metadata, not bonds.
            // Complexity review: these clones are linear in their retained
            // owned values, as Gemmi's `empty_copy` deep-copies the same
            // aggregates; no full BioStructureData clone or remapping scan occurs.
            connections: source.connections().to_vec(),
            cispeps: source.cispeps().to_vec(),
            mod_residues: source.mod_residues().to_vec(),
            helices: source.helices().to_vec(),
            sheets: source.sheets().to_vec(),
            metadata: source.metadata().clone(),
            // Gemmi❗✔️:   std::map<int, std::vector<int>> conect_map;
            // Gemmi❗✔️:   bool has_d_fraction = false;  // uses Refmac's ccp4_deuterium_fraction
            // Gemmi❗✔️:   int non_ascii_line = 0;  // first PDB line with non-ASCII bytes, or 0
            // Gemmi❗✔️:   char ter_status = '\0';
            // M23 explicitly freezes preservation of these source-state fields
            // although `Structure::empty_copy` does not copy all of them. This
            // is the approved BIO projection contract, not a claim of exact
            // Gemmi empty_copy behavior. In particular, CONECT serial ids are
            // never reinterpreted as projected row ids or chemical bonds.
            source_state: source.source_state().clone(),
            coordinates: BioCoordinateBlock::new(positions),
            crystal: source.crystal().cloned(),
            ncs_operators: source.ncs_operators().to_vec(),
            // A generator referencing an excluded chain is not a valid
            // biological assembly of the projection. ProteinData is deliberately
            // not a lossless structure/writer value, so assemblies are absent
            // instead of being guessed or partially rewritten.
            assemblies: Vec::new(),
        })?;

        Ok(Self { structure })
    }

    #[must_use]
    pub fn num_models(&self) -> usize {
        self.structure.models().len()
    }

    #[must_use]
    pub fn num_chains(&self) -> usize {
        self.structure.chains().len()
    }

    #[must_use]
    pub fn num_residues(&self) -> usize {
        self.structure.residues().len()
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.structure.atoms().len()
    }

    #[must_use]
    pub fn selection_summary(&self) -> ProteinSelectionSummary {
        protein_selection_summary(&self.structure)
    }

    #[must_use]
    pub fn chains(&self) -> ProteinChainIter<'_> {
        protein_chains(&self.structure)
    }

    #[must_use]
    pub fn chain(&self, index: usize) -> Option<ProteinChainRef<'_>> {
        protein_chain(&self.structure, index)
    }

    #[must_use]
    pub fn residues(&self) -> ProteinResidueIter<'_> {
        protein_residues(&self.structure)
    }

    #[must_use]
    pub fn atoms(&self) -> ProteinAtomIter<'_> {
        protein_atoms(&self.structure)
    }
}

#[must_use]
pub fn protein_chains(structure: &BioStructureData) -> ProteinChainIter<'_> {
    ProteinChainIter::new(structure, 0..structure.chains().len())
}

#[must_use]
pub fn protein_chain(structure: &BioStructureData, index: usize) -> Option<ProteinChainRef<'_>> {
    (index < structure.chains().len()).then(|| ProteinChainRef {
        structure,
        chain_id: BioChainId::new(index as u32),
    })
}

#[must_use]
pub fn protein_residues(structure: &BioStructureData) -> ProteinResidueIter<'_> {
    ProteinResidueIter::new(structure, 0..structure.residues().len())
}

#[must_use]
pub fn protein_atoms(structure: &BioStructureData) -> ProteinAtomIter<'_> {
    ProteinAtomIter::new(structure, 0..structure.atoms().len())
}

/// Validate the already-projected hierarchy without constructing another owner.
pub fn validate_protein_structure(
    structure: &BioStructureData,
) -> Result<(), ProteinProjectionError> {
    structure.validate()?;
    for (index, residue) in structure.residues().iter().enumerate() {
        if !is_amino_acid_kind(residue.residue_info_kind()) {
            return Err(ProteinProjectionError::NonAminoAcidResidue { index });
        }
    }
    Ok(())
}

impl<'a> ProteinChainRef<'a> {
    #[must_use]
    pub const fn id(self) -> BioChainId {
        self.chain_id
    }

    #[must_use]
    pub fn row(self) -> &'a BioChainRow {
        &self.structure.chains()[self.chain_id.index()]
    }

    #[must_use]
    pub fn kind(self) -> ChainKind {
        self.row().kind()
    }

    #[must_use]
    pub fn source(self) -> &'a ChainSourceIds {
        self.row().source()
    }

    #[must_use]
    pub fn residues(self) -> ProteinResidueIter<'a> {
        let span = self.row().residue_span();
        ProteinResidueIter::new(self.structure, span.start() as usize..span.end() as usize)
    }

    #[must_use]
    pub fn atoms(self) -> ProteinAtomIter<'a> {
        let span = self.row().residue_span();
        ProteinAtomIter::from_residues(self.structure, span.start() as usize..span.end() as usize)
    }
}

impl<'a> ProteinResidueRef<'a> {
    #[must_use]
    pub const fn id(self) -> BioResidueId {
        self.residue_id
    }

    #[must_use]
    pub fn row(self) -> &'a BioResidueRow {
        &self.structure.residues()[self.residue_id.index()]
    }

    #[must_use]
    pub fn name(self) -> ResidueName {
        self.row().name()
    }

    #[must_use]
    pub fn kind(self) -> ResidueKind {
        self.row().kind()
    }

    #[must_use]
    pub fn info(self) -> ResidueInfo {
        find_residue_info(self.name().as_str())
    }

    #[must_use]
    pub fn code(self) -> ResidueCode {
        self.info().code
    }

    #[must_use]
    pub fn one_letter_code(self) -> char {
        self.info().one_letter_code
    }

    #[must_use]
    pub fn fasta_code(self) -> char {
        self.info().fasta_code()
    }

    #[must_use]
    pub fn is_standard(self) -> bool {
        self.info().is_standard()
    }

    #[must_use]
    pub fn chain(self) -> ProteinChainRef<'a> {
        ProteinChainRef {
            structure: self.structure,
            chain_id: self.row().chain_id(),
        }
    }

    #[must_use]
    pub fn atoms(self) -> ProteinAtomIter<'a> {
        let span = self.row().atom_span();
        ProteinAtomIter::new(self.structure, span.start() as usize..span.end() as usize)
    }
}

impl<'a> ProteinAtomRef<'a> {
    #[must_use]
    pub const fn id(self) -> BioAtomId {
        self.atom_id
    }

    #[must_use]
    pub fn row(self) -> &'a BioAtomRow {
        &self.structure.atoms()[self.atom_id.index()]
    }

    #[must_use]
    pub fn name(self) -> AtomName {
        self.row().name()
    }

    #[must_use]
    pub fn element(self) -> Element {
        self.row().element()
    }

    #[must_use]
    pub fn altloc(self) -> Option<AltLocLabel> {
        self.row().altloc()
    }

    #[must_use]
    pub fn residue(self) -> ProteinResidueRef<'a> {
        ProteinResidueRef {
            structure: self.structure,
            residue_id: self.row().residue_id(),
        }
    }

    #[must_use]
    pub fn position(self) -> [f64; 3] {
        self.structure.coordinates().positions()[self.atom_id.index()]
    }
}

fn is_amino_acid_kind(kind: ResidueInfoKind) -> bool {
    // Gemmi✔️✔️: bool is_amino_acid() const {
    // Gemmi✔️✔️:   return kind == ResidueKind::AA || kind == ResidueKind::AAD ||
    // Gemmi✔️✔️:          kind == ResidueKind::PAA || kind == ResidueKind::MAA;
    // Gemmi✔️✔️: }
    // Behavior review: the four source kinds map one-for-one to the Rust enum.
    // Complexity review: both implementations perform a constant number of
    // enum comparisons without allocation.
    matches!(
        kind,
        ResidueInfoKind::Aa | ResidueInfoKind::Aad | ResidueInfoKind::Paa | ResidueInfoKind::Maa
    )
}

fn checked_row_index(value: usize) -> Result<u32, BioStructureError> {
    u32::try_from(value).map_err(|_| BioStructureError::RowIndexTooLarge { value })
}

fn span_range<I, T>(
    span: BioRowSpan<I>,
    rows: &[T],
) -> Result<std::ops::Range<usize>, BioStructureError> {
    span.slice(rows)?;
    Ok(span.start() as usize..span.end() as usize)
}

fn remap_entity(
    source: Option<BioEntityId>,
    entity_map: &[Option<BioEntityId>],
) -> Result<Option<BioEntityId>, ProteinProjectionError> {
    source
        .map(|source_id| {
            entity_map
                .get(source_id.index())
                .copied()
                .flatten()
                .ok_or(ProteinProjectionError::MissingEntityMapping { source: source_id })
        })
        .transpose()
}
