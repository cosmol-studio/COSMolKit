//! RDKit-compatible input preparation, separate from the official InChI engine.
//!
//! Implementors supply already Kekulized working data. Atom and bond indices and
//! neighbor iteration order must remain stable throughout generation.
use crate::{InchiAtom, InchiBondDirection, InchiBondStereo, InchiBondType, InchiToolkitError};

/// Generate using caller-prepared data and the canonical RDKit input encoder.
/// This does not update the property cache or Kekulize the input.
pub fn generate(
    graph: &mut impl PreparedGraph,
    options: Option<&[u8]>,
) -> Result<crate::MolToInchiOutput, crate::InchiError> {
    use crate::source::rdkit::inchi;
    let mut heap = crate::source_types::SourceHeap::default();
    let mut engine = inchi::SourceInchiGenerationEngine::new(&mut heap);
    let mut return_values = inchi::ExtraInchiReturnValues::default();
    let result = inchi::generate_prepared(
        &mut engine,
        graph,
        &mut return_values,
        options.map(inchi::GenerationOptions::Bytes),
    )
    .map_err(|error| crate::generation_error("mol_to_inchi", error))?;
    Ok(crate::MolToInchiOutput {
        inchi: result.inchi,
        return_values: crate::InchiReturnValues {
            return_code: return_values.return_code,
            message: return_values.message,
            log: return_values.log,
            aux_info: return_values.aux_info,
        },
        diagnostics: result.diagnostics,
    })
}

/// Scalar bond data consumed by the RDKit adapter; no owned graph or stereo list.
pub struct BondInput {
    pub begin_atom_index: u32,
    pub end_atom_index: u32,
    pub bond_type: InchiBondType,
    pub direction: InchiBondDirection,
    pub stereo: InchiBondStereo,
    pub stereo_atoms: Option<[u32; 2]>,
}

/// Narrow graph access used by RDKit reverse cleanup and input encoding.
pub trait Graph {
    fn atom_count(&self) -> usize;
    fn bond_count(&self) -> usize;
    fn atom(&self, index: usize) -> InchiAtom;
    fn bond(&self, index: usize) -> BondInput;
    fn degree(&self, index: usize) -> usize;
    fn neighbor(&self, index: usize, offset: usize) -> (u32, u32);
    fn neighbors(&self, index: usize) -> impl Iterator<Item = (u32, u32)> {
        (0..self.degree(index)).map(move |offset| self.neighbor(index, offset))
    }
}

/// Mutable, detached input after the caller's cache update and Kekulé step.
/// The adapter modifies only charges, bond orders and implicit-valence cache.
pub trait CleanupGraph: Graph {
    fn set_charge(&mut self, index: usize, charge: i8);
    fn set_bond_type(&mut self, index: usize, order: InchiBondType);
}

/// Generation access, including the caller-owned valence calculations.
pub trait PreparedGraph: CleanupGraph {
    fn conformer_len(&self) -> Option<usize>;
    fn position(&self, index: usize) -> [f64; 3];
    fn element_symbol(&mut self, number: i32) -> Result<Vec<u8>, InchiToolkitError>;
    fn atomic_weight(&mut self, number: i32) -> Result<f64, InchiToolkitError>;
    fn total_num_hydrogens(&mut self, index: u32) -> Result<u32, InchiToolkitError>;
    fn calc_implicit_valence(&mut self, index: u32) -> Result<i32, InchiToolkitError>;
    fn total_degree(&mut self, index: u32) -> Result<u32, InchiToolkitError>;
}
