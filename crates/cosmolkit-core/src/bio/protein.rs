use crate::bio::resinfo::{ResidueCode, ResidueInfo};
use crate::bio::{
    AtomId, AtomRow, BioStructure, ChainId, ChainKind, ChainRow, CoordinateBlock, ModelId,
    ModelRow, ResidueId, ResidueKind, ResidueRow, RowSpan,
};
#[cfg(feature = "io")]
use crate::{BioCoorFormat, BioReadError};

#[derive(Debug, Clone, PartialEq)]
pub struct Protein {
    structure: BioStructure,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProteinSelectionSummary {
    pub chains: usize,
    pub residues: usize,
    pub atoms: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct ProteinChainRef<'a> {
    protein: &'a Protein,
    chain_id: ChainId,
}

#[derive(Debug, Clone, Copy)]
pub struct ProteinResidueRef<'a> {
    protein: &'a Protein,
    residue_id: ResidueId,
}

#[derive(Debug, Clone, Copy)]
pub struct ProteinAtomRef<'a> {
    protein: &'a Protein,
    atom_id: AtomId,
}

impl Protein {
    #[must_use]
    pub fn project_from_bio_structure(structure: &BioStructure) -> Self {
        Self {
            structure: project_protein_bio_structure(structure),
        }
    }

    #[must_use]
    pub fn as_bio_structure(&self) -> &BioStructure {
        &self.structure
    }

    /// Returns the already protein-filtered structural value owned by this
    /// projection. This does not recover ligands, waters, nucleic acids, or
    /// other rows removed when the `Protein` was created.
    #[must_use]
    pub fn into_bio_structure(self) -> BioStructure {
        self.structure
    }

    #[cfg(feature = "io")]
    pub fn from_pdb_str(text: &str) -> Result<Self, BioReadError> {
        let structure = BioStructure::from_pdb_str(text)?;
        Ok(Self::project_from_bio_structure(&structure))
    }

    #[cfg(feature = "io")]
    pub fn from_pdb(path: &str) -> Result<Self, BioReadError> {
        let structure = BioStructure::from_pdb(path)?;
        Ok(Self::project_from_bio_structure(&structure))
    }

    #[cfg(feature = "io")]
    pub fn from_mmcif_str(text: &str, path: &str) -> Result<Self, BioReadError> {
        Self::from_str_with_format(text, path, BioCoorFormat::Mmcif)
    }

    #[cfg(feature = "io")]
    pub fn from_mmcif(path: &str) -> Result<Self, BioReadError> {
        let structure = BioStructure::from_mmcif(path)?;
        Ok(Self::project_from_bio_structure(&structure))
    }

    #[cfg(feature = "io")]
    pub fn from_str_with_format(
        text: &str,
        path: &str,
        format: BioCoorFormat,
    ) -> Result<Self, BioReadError> {
        let structure = BioStructure::from_str_with_format(text, path, format)?;
        Ok(Self::project_from_bio_structure(&structure))
    }

    #[cfg(feature = "io")]
    pub fn from_structure_str(text: &str, path: &str) -> Result<Self, BioReadError> {
        let structure = BioStructure::from_structure_str(text, path)?;
        Ok(Self::project_from_bio_structure(&structure))
    }

    #[must_use]
    pub fn selection_summary(&self) -> ProteinSelectionSummary {
        ProteinSelectionSummary {
            chains: self.structure.num_chains(),
            residues: self.structure.num_residues(),
            atoms: self.structure.num_atoms(),
        }
    }

    #[must_use]
    pub fn num_models(&self) -> usize {
        self.structure.num_models()
    }

    #[must_use]
    pub fn num_chains(&self) -> usize {
        self.structure.num_chains()
    }

    #[must_use]
    pub fn num_residues(&self) -> usize {
        self.structure.num_residues()
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.structure.num_atoms()
    }

    #[must_use]
    pub fn chains(&self) -> impl Iterator<Item = ProteinChainRef<'_>> {
        (0..self.structure.chains().len()).map(|index| ProteinChainRef {
            protein: self,
            chain_id: ChainId::new(index as u32),
        })
    }

    #[must_use]
    pub fn chain(&self, index: usize) -> Option<ProteinChainRef<'_>> {
        (index < self.structure.chains().len()).then_some(ProteinChainRef {
            protein: self,
            chain_id: ChainId::new(index as u32),
        })
    }

    #[must_use]
    pub fn residues(&self) -> impl Iterator<Item = ProteinResidueRef<'_>> {
        (0..self.structure.residues().len()).map(|index| ProteinResidueRef {
            protein: self,
            residue_id: ResidueId::new(index as u32),
        })
    }

    #[must_use]
    pub fn atoms(&self) -> impl Iterator<Item = ProteinAtomRef<'_>> {
        (0..self.structure.atoms().len()).map(|index| ProteinAtomRef {
            protein: self,
            atom_id: AtomId::new(index as u32),
        })
    }
}

impl<'a> ProteinChainRef<'a> {
    #[must_use]
    pub fn id(self) -> ChainId {
        self.chain_id
    }

    #[must_use]
    pub fn row(self) -> &'a ChainRow {
        &self.protein.structure.chains()[self.chain_id.index() as usize]
    }

    #[must_use]
    pub fn kind(self) -> ChainKind {
        self.row().kind
    }

    #[must_use]
    pub fn residues(self) -> impl Iterator<Item = ProteinResidueRef<'a>> {
        let span = self.row().residue_span;
        let protein = self.protein;
        (span.start as usize..span.end() as usize).map(move |index| ProteinResidueRef {
            protein,
            residue_id: ResidueId::new(index as u32),
        })
    }

    #[must_use]
    pub fn atoms(self) -> impl Iterator<Item = ProteinAtomRef<'a>> {
        self.residues().flat_map(ProteinResidueRef::atoms)
    }
}

impl<'a> ProteinResidueRef<'a> {
    #[must_use]
    pub fn id(self) -> ResidueId {
        self.residue_id
    }

    #[must_use]
    pub fn row(self) -> &'a ResidueRow {
        &self.protein.structure.residues()[self.residue_id.index() as usize]
    }

    #[must_use]
    pub fn kind(self) -> ResidueKind {
        self.row().kind
    }

    #[must_use]
    pub fn name(self) -> &'a str {
        self.row().name.as_str()
    }

    #[must_use]
    pub fn code(self) -> ResidueCode {
        self.info().code
    }

    #[must_use]
    pub fn info(self) -> ResidueInfo {
        crate::bio::resinfo::find_tabulated_residue(self.name())
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
            chain_id: self.row().chain_id,
        }
    }

    #[must_use]
    pub fn atoms(self) -> impl Iterator<Item = ProteinAtomRef<'a>> {
        let span = self.row().atom_span;
        let protein = self.protein;
        (span.start as usize..span.end() as usize).map(move |index| ProteinAtomRef {
            protein,
            atom_id: AtomId::new(index as u32),
        })
    }
}

impl<'a> ProteinAtomRef<'a> {
    #[must_use]
    pub fn id(self) -> AtomId {
        self.atom_id
    }

    #[must_use]
    pub fn row(self) -> &'a AtomRow {
        &self.protein.structure.atoms()[self.atom_id.index() as usize]
    }

    #[must_use]
    pub fn residue(self) -> ProteinResidueRef<'a> {
        ProteinResidueRef {
            protein: self.protein,
            residue_id: self.row().residue_id,
        }
    }

    #[must_use]
    pub fn position(self) -> Option<[f64; 3]> {
        self.protein.structure.atom_position(self.atom_id)
    }
}

fn project_protein_bio_structure(source: &BioStructure) -> BioStructure {
    let mut structure = BioStructure::new();
    structure.name = source.name().to_string();
    structure.input_format = source.input_format();
    structure.metadata = source.metadata().clone();
    structure.crystal = source.crystal().cloned();
    structure.resolution = source.resolution();
    structure.ter_status = source.ter_status();
    structure.mod_residues = source.mod_residues().to_vec();
    structure.assemblies = source.assemblies().to_vec();
    structure.origx = *source.origx();
    structure.has_origx = source.has_origx();
    structure.ncs_operators = source.ncs_operators().to_vec();
    structure.ncs_oper_identity_id = source.ncs_oper_identity_id().map(str::to_string);

    let mut positions = Vec::new();

    for model in source.models() {
        let chain_start = structure.chains.len() as u32;
        let source_chain_start = model.chain_span.start as usize;
        let source_chain_end = model.chain_span.end() as usize;

        for chain_index in source_chain_start..source_chain_end {
            let chain = &source.chains()[chain_index];
            let residue_start = structure.residues.len() as u32;
            let source_residue_start = chain.residue_span.start as usize;
            let source_residue_end = chain.residue_span.end() as usize;

            for residue_index in source_residue_start..source_residue_end {
                let residue = &source.residues()[residue_index];
                if residue.kind != ResidueKind::AminoAcid {
                    continue;
                }

                let atom_start = structure.atoms.len() as u32;
                let source_atom_start = residue.atom_span.start as usize;
                let source_atom_end = residue.atom_span.end() as usize;

                for atom_index in source_atom_start..source_atom_end {
                    let atom = &source.atoms()[atom_index];
                    let mut copied = atom.clone();
                    copied.residue_id = ResidueId::new(structure.residues.len() as u32);
                    structure.atoms.push(copied);
                    positions.push(source.coordinates().positions()[atom_index]);
                }

                let atom_len = structure.atoms.len() as u32 - atom_start;
                if atom_len == 0 {
                    continue;
                }

                let mut copied_residue = residue.clone();
                copied_residue.chain_id = ChainId::new(structure.chains.len() as u32);
                copied_residue.atom_span = RowSpan::new(atom_start, atom_len);
                structure.residues.push(copied_residue);
            }

            let residue_len = structure.residues.len() as u32 - residue_start;
            if residue_len == 0 {
                continue;
            }

            let mut copied_chain = chain.clone();
            copied_chain.model_id = ModelId::new(structure.models.len() as u32);
            copied_chain.residue_span = RowSpan::new(residue_start, residue_len);
            copied_chain.kind = ChainKind::Protein;
            structure.chains.push(copied_chain);
        }

        let chain_len = structure.chains.len() as u32 - chain_start;
        if chain_len == 0 {
            continue;
        }

        let copied_model = ModelRow {
            chain_span: RowSpan::new(chain_start, chain_len),
            source_model_number: model.source_model_number,
        };
        structure.models.push(copied_model);
    }

    structure.coordinates = CoordinateBlock { positions };
    structure
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bio::ResidueKind;

    const MIXED_PDB: &str = "\
ATOM      1  N   MET A   1      11.104  13.207   9.900  1.00 20.00           N
ATOM      2  CA  MET A   1      12.210  13.912  10.555  1.00 20.00           C
ATOM      3  C   MET A   1      13.470  13.079  10.413  1.00 20.00           C
ATOM      4  N   GLY A   2      14.530  13.650  10.980  1.00 20.00           N
ATOM      5  CA  GLY A   2      15.790  12.920  10.910  1.00 20.00           C
HETATM    6  O   HOH A   3      18.000  10.000   8.000  1.00 10.00           O
HETATM    7  C1  LIG B   1      18.500  11.000   8.500  1.00 10.00           C
";

    #[test]
    fn protein_projection_keeps_only_amino_acid_residues() {
        let structure = BioStructure::from_pdb_str(MIXED_PDB).unwrap();
        let protein = Protein::project_from_bio_structure(&structure);

        assert_eq!(protein.num_models(), 1);
        assert_eq!(protein.num_chains(), 1);
        assert_eq!(protein.num_residues(), 2);
        assert_eq!(protein.num_atoms(), 5);
        assert_eq!(
            protein.chain(0).expect("projected protein chain").kind(),
            ChainKind::Protein
        );
        assert!(
            protein
                .as_bio_structure()
                .residues()
                .iter()
                .all(|row| row.kind == ResidueKind::AminoAcid)
        );
    }

    #[test]
    fn protein_chain_residue_atom_navigation_is_hierarchical() {
        let structure = BioStructure::from_pdb_str(MIXED_PDB).unwrap();
        let protein = Protein::project_from_bio_structure(&structure);
        let chain = protein.chain(0).expect("projected protein chain");

        let residues: Vec<_> = chain.residues().collect();
        assert_eq!(residues.len(), 2);
        assert_eq!(residues[0].name(), "MET");
        assert_eq!(residues[1].name(), "GLY");

        let residue_atoms: Vec<_> = residues[0].atoms().collect();
        assert_eq!(residue_atoms.len(), 3);
        assert_eq!(residue_atoms[0].row().name.0, [b' ', b'N', b' ', b' ']);
        assert!(residue_atoms[0].position().is_some());
    }
}
