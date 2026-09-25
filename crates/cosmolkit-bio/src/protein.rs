//! Amino-acid-only projection and read-only hierarchy views.

use std::fmt;

use cosmolkit_types::Element;

use crate::{
    AltLocLabel, AtomName, BioAtomId, BioAtomRow, BioChainId, BioChainRow, BioCoordinateBlock,
    BioEntityId, BioModelId, BioModelRow, BioResidueId, BioResidueRow, BioRowSpan, BioStructure,
    BioStructureError, BioStructureParts, ChainKind, ChainSourceIds, ResidueCode, ResidueInfo,
    ResidueInfoKind, ResidueKind, ResidueName, find_residue_info,
};

#[derive(Debug, Clone, PartialEq)]
pub struct Protein {
    structure: BioStructure,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinProjectionError {
    Structure(BioStructureError),
    MissingEntityMapping { source: BioEntityId },
}

impl fmt::Display for ProteinProjectionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
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
            Self::MissingEntityMapping { .. } => None,
        }
    }
}

impl From<BioStructureError> for ProteinProjectionError {
    fn from(value: BioStructureError) -> Self {
        Self::Structure(value)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ProteinChainRef<'a> {
    protein: &'a Protein,
    chain_id: BioChainId,
}

#[derive(Debug, Clone, Copy)]
pub struct ProteinResidueRef<'a> {
    protein: &'a Protein,
    residue_id: BioResidueId,
}

#[derive(Debug, Clone, Copy)]
pub struct ProteinAtomRef<'a> {
    protein: &'a Protein,
    atom_id: BioAtomId,
}

impl BioStructure {
    pub fn protein(&self) -> Result<Protein, ProteinProjectionError> {
        Protein::project(self)
    }
}

impl Protein {
    fn project(source: &BioStructure) -> Result<Self, ProteinProjectionError> {
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

        let structure = BioStructure::from_parts(BioStructureParts {
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
            // aggregates; no full BioStructure clone or remapping scan occurs.
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
            // biological assembly of the projection. Protein is deliberately
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
    pub fn chains(&self) -> Vec<ProteinChainRef<'_>> {
        (0..self.num_chains())
            .map(|index| ProteinChainRef {
                protein: self,
                chain_id: BioChainId::new(index as u32),
            })
            .collect()
    }

    #[must_use]
    pub fn chain(&self, index: usize) -> Option<ProteinChainRef<'_>> {
        (index < self.num_chains()).then(|| ProteinChainRef {
            protein: self,
            chain_id: BioChainId::new(index as u32),
        })
    }

    #[must_use]
    pub fn residues(&self) -> Vec<ProteinResidueRef<'_>> {
        (0..self.num_residues())
            .map(|index| ProteinResidueRef {
                protein: self,
                residue_id: BioResidueId::new(index as u32),
            })
            .collect()
    }

    #[must_use]
    pub fn atoms(&self) -> Vec<ProteinAtomRef<'_>> {
        (0..self.num_atoms())
            .map(|index| ProteinAtomRef {
                protein: self,
                atom_id: BioAtomId::new(index as u32),
            })
            .collect()
    }
}

impl<'a> ProteinChainRef<'a> {
    #[must_use]
    pub const fn id(self) -> BioChainId {
        self.chain_id
    }

    #[must_use]
    pub fn row(self) -> &'a BioChainRow {
        &self.protein.structure.chains()[self.chain_id.index()]
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
    pub fn residues(self) -> Vec<ProteinResidueRef<'a>> {
        let span = self.row().residue_span();
        (span.start() as usize..span.end() as usize)
            .map(|index| ProteinResidueRef {
                protein: self.protein,
                residue_id: BioResidueId::new(index as u32),
            })
            .collect()
    }

    #[must_use]
    pub fn atoms(self) -> Vec<ProteinAtomRef<'a>> {
        self.residues()
            .into_iter()
            .flat_map(ProteinResidueRef::atoms)
            .collect()
    }
}

impl<'a> ProteinResidueRef<'a> {
    #[must_use]
    pub const fn id(self) -> BioResidueId {
        self.residue_id
    }

    #[must_use]
    pub fn row(self) -> &'a BioResidueRow {
        &self.protein.structure.residues()[self.residue_id.index()]
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
            protein: self.protein,
            chain_id: self.row().chain_id(),
        }
    }

    #[must_use]
    pub fn atoms(self) -> Vec<ProteinAtomRef<'a>> {
        let span = self.row().atom_span();
        (span.start() as usize..span.end() as usize)
            .map(|index| ProteinAtomRef {
                protein: self.protein,
                atom_id: BioAtomId::new(index as u32),
            })
            .collect()
    }
}

impl<'a> ProteinAtomRef<'a> {
    #[must_use]
    pub const fn id(self) -> BioAtomId {
        self.atom_id
    }

    #[must_use]
    pub fn row(self) -> &'a BioAtomRow {
        &self.protein.structure.atoms()[self.atom_id.index()]
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
            protein: self.protein,
            residue_id: self.row().residue_id(),
        }
    }

    #[must_use]
    pub fn position(self) -> [f64; 3] {
        self.protein.structure.coordinates().positions()[self.atom_id.index()]
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
