use std::collections::{BTreeSet, HashSet, VecDeque};

use crate::source::api::inchi_dll::{FreeINCHI, FreeStructFromINCHI, GetINCHI, GetStructFromINCHI};
use crate::source::base::ichiparm::InchiBuildMetadata;
use crate::source::base::ikey_dll::GetINCHIKeyFromINCHI;
use crate::source_types::{
    FILE, INCHIKEY_EMPTY_INPUT, INCHIKEY_INVALID_INCHI, INCHIKEY_INVALID_INCHI_PREFIX, INCHIKEY_INVALID_STD_INCHI,
    INCHIKEY_NOT_ENOUGH_MEMORY, INCHIKEY_OK, INCHIKEY_UNKNOWN_ERROR, ISOTOPIC_SHIFT_FLAG, MAXVAL, NO_ATOM, SourceHeap,
    SourceHeapError, SourceMutPointer, clock_t, inchi_Atom, inchi_Input, inchi_InputINCHI, inchi_Output,
    inchi_OutputStruct, inchi_Stereo0D, tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE, tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1DOWN,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1EITHER, tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1UP,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2DOWN, tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2UP,
    tagINCHIBondType_INCHI_BOND_TYPE_ALTERN, tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE,
    tagINCHIStereoParity0D_INCHI_PARITY_EVEN, tagINCHIStereoParity0D_INCHI_PARITY_NONE,
    tagINCHIStereoParity0D_INCHI_PARITY_ODD, tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED,
    tagINCHIStereoType0D_INCHI_StereoType_Allene, tagINCHIStereoType0D_INCHI_StereoType_DoubleBond,
    tagINCHIStereoType0D_INCHI_StereoType_None, tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral,
    tagRetValGetINCHI_inchi_Ret_OKAY, tagRetValGetINCHI_inchi_Ret_WARNING,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum BondDirection {
    None = 0,
    BeginWedge = 1,
    BeginDash = 2,
    EndDownRight = 3,
    EndUpRight = 4,
    EitherDouble = 5,
    Unknown = 6,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum BondType {
    Unspecified = 0,
    Single = 1,
    Double = 2,
    Triple = 3,
    Quadruple = 4,
    Quintuple = 5,
    Hextuple = 6,
    OneAndAHalf = 7,
    TwoAndAHalf = 8,
    ThreeAndAHalf = 9,
    FourAndAHalf = 10,
    FiveAndAHalf = 11,
    Aromatic = 12,
    Ionic = 13,
    Hydrogen = 14,
    ThreeCenter = 15,
    DativeOne = 16,
    Dative = 17,
    DativeL = 18,
    DativeR = 19,
    Other = 20,
    Zero = 21,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(u8)]
pub enum BondStereo {
    #[default]
    None = 0,
    Any = 1,
    Z = 2,
    E = 3,
    Cis = 4,
    Trans = 5,
    AtropCw = 6,
    AtropCcw = 7,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(u8)]
pub enum ChiralTag {
    #[default]
    Unspecified = 0,
    TetrahedralCw = 1,
    TetrahedralCcw = 2,
    Other = 3,
    Tetrahedral = 4,
    Allene = 5,
    SquarePlanar = 6,
    TrigonalBipyramidal = 7,
    Octahedral = 8,
}

#[derive(Clone, Debug)]
pub struct AdapterAtom {
    pub atomic_number: i32,
    pub formal_charge: i32,
    pub num_explicit_hydrogens: u32,
    pub is_aromatic: bool,
    pub isotope: u32,
    pub num_radical_electrons: u32,
    pub no_implicit: bool,
    pub chiral_tag: ChiralTag,
    pub cip_rank: Option<u32>,
    /// RDKit's per-atom property-cache value. This is deliberately separate
    /// from the current graph: `InchiToMol(..., sanitize=false)` can return a
    /// cache computed before `cleanUp()` rewrites a bond.
    pub cached_explicit_valence: Option<i32>,
    pub cached_implicit_valence: Option<i32>,
}

// Property-cache entries are derived observations, not part of the adapter's
// chemical graph identity. Keeping them out of equality preserves the
// pre-cache `AdapterAtom` equality contract; parity tests inspect the cached
// values explicitly where RDKit exposes them.
impl PartialEq for AdapterAtom {
    fn eq(&self, other: &Self) -> bool {
        self.atomic_number == other.atomic_number
            && self.formal_charge == other.formal_charge
            && self.num_explicit_hydrogens == other.num_explicit_hydrogens
            && self.is_aromatic == other.is_aromatic
            && self.isotope == other.isotope
            && self.num_radical_electrons == other.num_radical_electrons
            && self.no_implicit == other.no_implicit
            && self.chiral_tag == other.chiral_tag
            && self.cip_rank == other.cip_rank
    }
}

impl Eq for AdapterAtom {}

impl Default for AdapterAtom {
    fn default() -> Self {
        Self {
            atomic_number: 0,
            formal_charge: 0,
            num_explicit_hydrogens: 0,
            is_aromatic: false,
            isotope: 0,
            num_radical_electrons: 0,
            no_implicit: false,
            chiral_tag: ChiralTag::Unspecified,
            cip_rank: None,
            cached_explicit_valence: None,
            cached_implicit_valence: None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AdapterBond {
    pub(crate) begin_atom_index: u32,
    pub(crate) end_atom_index: u32,
    pub bond_type: BondType,
    pub direction: BondDirection,
    pub is_aromatic: bool,
    pub stereo: BondStereo,
    pub stereo_atoms: Vec<u32>,
}

impl Default for AdapterBond {
    fn default() -> Self {
        Self {
            begin_atom_index: 0,
            end_atom_index: 0,
            bond_type: BondType::Unspecified,
            direction: BondDirection::None,
            is_aromatic: false,
            stereo: BondStereo::None,
            stereo_atoms: Vec::new(),
        }
    }
}

impl AdapterBond {
    pub fn new(begin_atom_index: u32, end_atom_index: u32, bond_type: BondType) -> Self {
        Self {
            begin_atom_index,
            end_atom_index,
            bond_type,
            ..Self::default()
        }
    }

    pub fn begin_atom_index(&self) -> u32 {
        self.begin_atom_index
    }

    pub fn end_atom_index(&self) -> u32 {
        self.end_atom_index
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdapterMol {
    pub(crate) atoms: Vec<AdapterAtom>,
    pub(crate) bonds: Vec<AdapterBond>,
    pub(crate) conformers: Vec<Vec<[f64; 3]>>,
    adjacency: Vec<Vec<(u32, u32)>>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AdapterDiagnosticLevel {
    Warning,
    Error,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AdapterDiagnostic {
    pub level: AdapterDiagnosticLevel,
    pub message: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct ExtraInchiReturnValues {
    pub(crate) return_code: i32,
    pub(crate) message: Vec<u8>,
    pub(crate) log: Vec<u8>,
    pub(crate) aux_info: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AdapterInchiStructureOutput {
    pub(crate) return_code: i32,
    pub(crate) atoms: Vec<inchi_Atom>,
    pub(crate) stereo0d: Vec<inchi_Stereo0D>,
    pub(crate) message: Option<Vec<u8>>,
    pub(crate) log: Option<Vec<u8>>,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct InchiToMolResult {
    pub(crate) molecule: Option<AdapterMol>,
    pub(crate) diagnostics: Vec<AdapterDiagnostic>,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AdapterInchiGenerationInput {
    pub(crate) atoms: Vec<inchi_Atom>,
    pub(crate) stereo0d: Vec<inchi_Stereo0D>,
    pub(crate) options_with_nul: Option<Vec<u8>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AdapterInchiGenerationOutput {
    pub(crate) return_code: i32,
    pub(crate) inchi: Option<Vec<u8>>,
    pub(crate) message: Option<Vec<u8>>,
    pub(crate) log: Option<Vec<u8>>,
    pub(crate) aux_info: Option<Vec<u8>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MolToInchiResult {
    pub(crate) inchi: Vec<u8>,
    pub(crate) diagnostics: Vec<AdapterDiagnostic>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MolToInchiKeyResult {
    pub(crate) key: Vec<u8>,
    pub(crate) diagnostics: Vec<AdapterDiagnostic>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AdapterInchiKeyOutput {
    pub(crate) status: i32,
    pub(crate) key_buffer: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct InchiToInchiKeyResult {
    pub(crate) key: Vec<u8>,
    pub(crate) diagnostics: Vec<AdapterDiagnostic>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum InchiToInchiKeyError {
    Source(SourceHeapError),
    InvalidSourceOutput(&'static str),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum MolToInchiKeyError {
    MolToInchi(MolToInchiError),
    InchiToInchiKey(InchiToInchiKeyError),
}

impl From<MolToInchiError> for MolToInchiKeyError {
    fn from(value: MolToInchiError) -> Self {
        Self::MolToInchi(value)
    }
}

impl From<InchiToInchiKeyError> for MolToInchiKeyError {
    fn from(value: InchiToInchiKeyError) -> Self {
        Self::InchiToInchiKey(value)
    }
}

impl From<SourceHeapError> for InchiToInchiKeyError {
    fn from(value: SourceHeapError) -> Self {
        Self::Source(value)
    }
}

pub(crate) trait InchiKeyEngine {
    fn get_inchi_key(&mut self, inchi_with_nul: &[u8]) -> Result<AdapterInchiKeyOutput, SourceHeapError>;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum MolToInchiError {
    Source(SourceHeapError),
    Toolkit(AdapterToolkitError),
    InvalidOptions,
    InvalidConformer,
    ElementSymbolTooLong,
}

impl From<SourceHeapError> for MolToInchiError {
    fn from(value: SourceHeapError) -> Self {
        Self::Source(value)
    }
}

impl From<AdapterToolkitError> for MolToInchiError {
    fn from(value: AdapterToolkitError) -> Self {
        Self::Toolkit(value)
    }
}

pub(crate) trait InchiGenerationEngine {
    fn get_inchi(
        &mut self,
        input: &AdapterInchiGenerationInput,
    ) -> Result<AdapterInchiGenerationOutput, SourceHeapError>;

    fn free_inchi(&mut self) -> Result<(), SourceHeapError>;
}

pub trait MolToInchiToolkit {
    fn needs_update_property_cache(&mut self, molecule: &AdapterMol) -> Result<bool, AdapterToolkitError>;

    fn update_property_cache(&mut self, molecule: &mut AdapterMol, strict: bool) -> Result<(), AdapterToolkitError>;

    fn kekulize(&mut self, molecule: &mut AdapterMol, mark_atoms_bonds: bool) -> Result<(), AdapterToolkitError>;

    fn element_symbol(&mut self, atomic_number: i32) -> Result<Vec<u8>, AdapterToolkitError>;

    fn atomic_weight(&mut self, atomic_number: i32) -> Result<f64, AdapterToolkitError>;

    fn total_num_hydrogens(&mut self, molecule: &AdapterMol, atom_index: u32) -> Result<u32, AdapterToolkitError>;

    fn calc_implicit_valence(&mut self, molecule: &mut AdapterMol, atom_index: u32)
    -> Result<i32, AdapterToolkitError>;

    fn total_degree(&mut self, molecule: &AdapterMol, atom_index: u32) -> Result<u32, AdapterToolkitError>;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AdapterToolkitError {
    pub kind: &'static str,
    pub message: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum InchiToMolError {
    Source(SourceHeapError),
    Toolkit(AdapterToolkitError),
    Cleanup(AdapterCleanup5Error),
    InvalidSourceOutput(&'static str),
    BondIndex(BondIndexRangeError),
}

impl From<SourceHeapError> for InchiToMolError {
    fn from(value: SourceHeapError) -> Self {
        Self::Source(value)
    }
}

impl From<AdapterToolkitError> for InchiToMolError {
    fn from(value: AdapterToolkitError) -> Self {
        Self::Toolkit(value)
    }
}

impl From<AdapterCleanup5Error> for InchiToMolError {
    fn from(value: AdapterCleanup5Error) -> Self {
        Self::Cleanup(value)
    }
}

impl From<BondIndexRangeError> for InchiToMolError {
    fn from(value: BondIndexRangeError) -> Self {
        Self::BondIndex(value)
    }
}

pub(crate) trait InchiStructureEngine {
    fn get_struct_from_inchi(&mut self, inchi_with_nul: &[u8]) -> Result<AdapterInchiStructureOutput, SourceHeapError>;

    fn free_struct_from_inchi(&mut self) -> Result<(), SourceHeapError>;
}

pub trait InchiToMolToolkit {
    fn atomic_number(&mut self, element: &[u8]) -> Result<i32, AdapterToolkitError>;

    fn average_atomic_weight(&mut self, atomic_number: i32) -> Result<f64, AdapterToolkitError>;

    fn update_property_cache(&mut self, molecule: &mut AdapterMol, strict: bool) -> Result<(), AdapterToolkitError>;

    fn assign_atom_cip_ranks(&mut self, molecule: &mut AdapterMol) -> Result<Vec<u32>, AdapterToolkitError>;

    fn remove_hydrogens(&mut self, molecule: &mut AdapterMol) -> Result<(), AdapterToolkitError>;

    fn sanitize_molecule(&mut self, molecule: &mut AdapterMol) -> Result<(), AdapterToolkitError>;

    /// Synchronizes toolkit-native state after RDKit's in-place `cleanUp` pass.
    ///
    /// The default is intentionally a no-op for adapters that operate directly
    /// on `AdapterMol`. Toolkits with a separate native molecule representation
    /// use this exact source boundary to retain the single-molecule ownership
    /// model of RDKit's `InchiToMol` implementation.
    fn synchronize_after_cleanup(&mut self, _molecule: &AdapterMol) -> Result<(), AdapterToolkitError> {
        Ok(())
    }

    fn assign_stereochemistry(
        &mut self,
        molecule: &mut AdapterMol,
        clean_it: bool,
        force: bool,
    ) -> Result<(), AdapterToolkitError>;
}

pub(crate) struct SourceInchiStructureEngine<'a> {
    pub(crate) heap: &'a mut SourceHeap,
    pub(crate) stdout: SourceMutPointer<FILE>,
    pub(crate) build: InchiBuildMetadata<'a>,
    pub(crate) clock_result: clock_t,
    pending_output: Option<inchi_OutputStruct>,
}

impl<'a> SourceInchiStructureEngine<'a> {
    pub(crate) fn new(heap: &'a mut SourceHeap) -> Self {
        Self {
            heap,
            stdout: SourceMutPointer::null(),
            build: InchiBuildMetadata {
                compiler: "gcc",
                date: "Jan 01 1970",
                time: "00:00:00",
            },
            clock_result: 1_000,
            pending_output: None,
        }
    }
}

fn source_c_text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Result<Vec<u8>, SourceHeapError> {
    let bytes = heap.slice(pointer.as_const())?;
    let end = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    Ok(bytes[..end].iter().map(|byte| *byte as u8).collect())
}

impl InchiStructureEngine for SourceInchiStructureEngine<'_> {
    fn get_struct_from_inchi(&mut self, inchi_with_nul: &[u8]) -> Result<AdapterInchiStructureOutput, SourceHeapError> {
        if self.pending_output.is_some() {
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        }
        let inchi = self
            .heap
            .allocate_model_storage(inchi_with_nul.iter().map(|byte| *byte as i8).collect())?;
        let options = self.heap.allocate_model_storage(vec![0_i8])?;
        let input = inchi_InputINCHI {
            szInChI: inchi,
            szOptions: options,
        };
        let mut output = inchi_OutputStruct::default();
        let return_code = GetStructFromINCHI(
            self.heap,
            Some(&input),
            &mut output,
            self.stdout,
            self.build,
            self.clock_result,
        );
        self.heap.free(inchi)?;
        self.heap.free(options)?;
        let return_code = return_code?;

        let message = if output.szMessage.is_null() {
            None
        } else {
            Some(source_c_text(self.heap, output.szMessage)?)
        };
        let log = if output.szLog.is_null() {
            None
        } else {
            Some(source_c_text(self.heap, output.szLog)?)
        };
        let atoms =
            if return_code == tagRetValGetINCHI_inchi_Ret_OKAY || return_code == tagRetValGetINCHI_inchi_Ret_WARNING {
                let count = usize::try_from(output.num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                if count == 0 {
                    Vec::new()
                } else {
                    self.heap
                        .slice(output.atom.as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec()
                }
            } else {
                Vec::new()
            };
        let stereo0d =
            if return_code == tagRetValGetINCHI_inchi_Ret_OKAY || return_code == tagRetValGetINCHI_inchi_Ret_WARNING {
                let count = usize::try_from(output.num_stereo0D).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                if count == 0 {
                    Vec::new()
                } else {
                    self.heap
                        .slice(output.stereo0D.as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec()
                }
            } else {
                Vec::new()
            };
        self.pending_output = Some(output);
        Ok(AdapterInchiStructureOutput {
            return_code,
            atoms,
            stereo0d,
            message,
            log,
        })
    }

    fn free_struct_from_inchi(&mut self) -> Result<(), SourceHeapError> {
        let mut output = self
            .pending_output
            .take()
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
        FreeStructFromINCHI(self.heap, Some(&mut output))
    }
}

pub(crate) struct SourceInchiGenerationEngine<'a> {
    pub(crate) heap: &'a mut SourceHeap,
    pub(crate) stdout: SourceMutPointer<FILE>,
    pub(crate) build: InchiBuildMetadata<'a>,
    pub(crate) clock_result: clock_t,
    pending_output: Option<PendingInchiGenerationOutput>,
}

impl<'a> SourceInchiGenerationEngine<'a> {
    pub(crate) fn new(heap: &'a mut SourceHeap) -> Self {
        Self {
            heap,
            stdout: SourceMutPointer::null(),
            build: InchiBuildMetadata {
                compiler: "gcc",
                date: "Jan 01 1970",
                time: "00:00:00",
            },
            clock_result: 1_000,
            pending_output: None,
        }
    }
}

struct PendingInchiGenerationOutput {
    output: inchi_Output,
    atoms: SourceMutPointer<inchi_Atom>,
    stereo0d: SourceMutPointer<inchi_Stereo0D>,
    options: SourceMutPointer<i8>,
}

fn free_generation_inputs(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inchi_Atom>,
    stereo0d: SourceMutPointer<inchi_Stereo0D>,
    options: SourceMutPointer<i8>,
) -> Result<(), SourceHeapError> {
    let mut first_error = None;
    for result in [heap.free(options), heap.free(stereo0d), heap.free(atoms)] {
        if first_error.is_none() {
            first_error = result.err();
        }
    }
    first_error.map_or(Ok(()), Err)
}

impl InchiGenerationEngine for SourceInchiGenerationEngine<'_> {
    fn get_inchi(
        &mut self,
        input: &AdapterInchiGenerationInput,
    ) -> Result<AdapterInchiGenerationOutput, SourceHeapError> {
        if self.pending_output.is_some() {
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        }
        let atoms = if input.atoms.is_empty() {
            SourceMutPointer::null()
        } else {
            self.heap.allocate(input.atoms.clone())?
        };
        let stereo0d = if input.stereo0d.is_empty() {
            SourceMutPointer::null()
        } else {
            match self.heap.allocate(input.stereo0d.clone()) {
                Ok(stereo0d) => stereo0d,
                Err(error) => {
                    self.heap.free(atoms)?;
                    return Err(error);
                }
            }
        };
        let options = if let Some(options) = &input.options_with_nul {
            match self.heap.allocate(options.iter().map(|byte| *byte as i8).collect()) {
                Ok(options) => options,
                Err(error) => {
                    self.heap.free(stereo0d)?;
                    self.heap.free(atoms)?;
                    return Err(error);
                }
            }
        } else {
            SourceMutPointer::null()
        };
        let source_input = inchi_Input {
            atom: atoms,
            stereo0D: stereo0d,
            szOptions: options,
            num_atoms: input.atoms.len() as i16,
            num_stereo0D: input.stereo0d.len() as i16,
        };
        let mut output = inchi_Output::default();
        let return_code = GetINCHI(
            self.heap,
            Some(&source_input),
            Some(&mut output),
            0,
            self.stdout,
            self.build,
            self.clock_result,
        );
        let return_code = match return_code {
            Ok(return_code) => return_code,
            Err(error) => {
                let output_cleanup = FreeINCHI(self.heap, Some(&mut output));
                let input_cleanup = free_generation_inputs(self.heap, atoms, stereo0d, options);
                output_cleanup?;
                input_cleanup?;
                return Err(error);
            }
        };
        let copied = (|| {
            Ok(AdapterInchiGenerationOutput {
                return_code,
                inchi: (!output.szInChI.is_null())
                    .then(|| source_c_text(self.heap, output.szInChI))
                    .transpose()?,
                message: (!output.szMessage.is_null())
                    .then(|| source_c_text(self.heap, output.szMessage))
                    .transpose()?,
                log: (!output.szLog.is_null())
                    .then(|| source_c_text(self.heap, output.szLog))
                    .transpose()?,
                aux_info: (!output.szAuxInfo.is_null())
                    .then(|| source_c_text(self.heap, output.szAuxInfo))
                    .transpose()?,
            })
        })();
        match copied {
            Ok(copied) => {
                self.pending_output = Some(PendingInchiGenerationOutput {
                    output,
                    atoms,
                    stereo0d,
                    options,
                });
                Ok(copied)
            }
            Err(error) => {
                let cleanup_output = FreeINCHI(self.heap, Some(&mut output));
                let cleanup_inputs = free_generation_inputs(self.heap, atoms, stereo0d, options);
                cleanup_output?;
                cleanup_inputs?;
                Err(error)
            }
        }
    }

    fn free_inchi(&mut self) -> Result<(), SourceHeapError> {
        let mut pending = self
            .pending_output
            .take()
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
        let output_result = FreeINCHI(self.heap, Some(&mut pending.output));
        let input_result = free_generation_inputs(self.heap, pending.atoms, pending.stereo0d, pending.options);
        output_result?;
        input_result
    }
}

pub(crate) struct SourceInchiKeyEngine<'a> {
    pub(crate) heap: &'a mut SourceHeap,
}

impl<'a> SourceInchiKeyEngine<'a> {
    pub(crate) fn new(heap: &'a mut SourceHeap) -> Self {
        Self { heap }
    }
}

impl InchiKeyEngine for SourceInchiKeyEngine<'_> {
    fn get_inchi_key(&mut self, inchi_with_nul: &[u8]) -> Result<AdapterInchiKeyOutput, SourceHeapError> {
        let input = self
            .heap
            .allocate_model_storage(inchi_with_nul.iter().map(|byte| *byte as i8).collect())?;
        let key = match self.heap.allocate_model_storage(vec![0_i8; 29]) {
            Ok(pointer) => pointer,
            Err(error) => {
                self.heap.free(input)?;
                return Err(error);
            }
        };
        let xtra1 = match self.heap.allocate_model_storage(vec![0_i8; 65]) {
            Ok(pointer) => pointer,
            Err(error) => {
                self.heap.free(key)?;
                self.heap.free(input)?;
                return Err(error);
            }
        };
        let xtra2 = match self.heap.allocate_model_storage(vec![0_i8; 65]) {
            Ok(pointer) => pointer,
            Err(error) => {
                self.heap.free(xtra1)?;
                self.heap.free(key)?;
                self.heap.free(input)?;
                return Err(error);
            }
        };
        let status = GetINCHIKeyFromINCHI(self.heap, input.as_const(), 0, 0, key, xtra1, xtra2);
        let copied = status.and_then(|status| {
            let key_buffer = if status == INCHIKEY_OK as i32 {
                self.heap
                    .slice(key.as_const())?
                    .get(..29)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect()
            } else {
                Vec::new()
            };
            Ok(AdapterInchiKeyOutput { status, key_buffer })
        });
        let mut cleanup_error = None;
        for result in [
            self.heap.free(xtra2),
            self.heap.free(xtra1),
            self.heap.free(key),
            self.heap.free(input),
        ] {
            if cleanup_error.is_none() {
                cleanup_error = result.err();
            }
        }
        if let Some(error) = cleanup_error {
            return Err(error);
        }
        copied
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct BondIndexRangeError {
    pub(crate) kind: &'static str,
    pub(crate) expression: &'static str,
    pub(crate) index: u32,
    pub(crate) upper_bound: u32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AdapterGraphError {
    pub bond_index: usize,
    pub atom_index: u32,
    pub atom_count: usize,
}

impl AdapterMol {
    pub fn try_from_graph(
        atoms: Vec<AdapterAtom>,
        bonds: Vec<AdapterBond>,
        conformers: Vec<Vec<[f64; 3]>>,
    ) -> Result<Self, AdapterGraphError> {
        let mut adjacency = vec![Vec::new(); atoms.len()];
        for (bond_index, bond) in bonds.iter().enumerate() {
            for atom_index in [bond.begin_atom_index, bond.end_atom_index] {
                let Ok(index) = usize::try_from(atom_index) else {
                    return Err(AdapterGraphError {
                        bond_index,
                        atom_index,
                        atom_count: atoms.len(),
                    });
                };
                if index >= atoms.len() {
                    return Err(AdapterGraphError {
                        bond_index,
                        atom_index,
                        atom_count: atoms.len(),
                    });
                }
            }
            let begin = bond.begin_atom_index as usize;
            let end = bond.end_atom_index as usize;
            let graph_bond_index = u32::try_from(bond_index).map_err(|_| AdapterGraphError {
                bond_index,
                atom_index: bond.begin_atom_index,
                atom_count: atoms.len(),
            })?;
            adjacency[begin].push((bond.end_atom_index, graph_bond_index));
            adjacency[end].push((bond.begin_atom_index, graph_bond_index));
        }
        Ok(Self {
            atoms,
            bonds,
            conformers,
            adjacency,
        })
    }

    pub fn atoms(&self) -> &[AdapterAtom] {
        &self.atoms
    }

    pub fn atom_properties_mut(&mut self) -> &mut [AdapterAtom] {
        &mut self.atoms
    }

    pub fn bonds(&self) -> &[AdapterBond] {
        &self.bonds
    }

    pub fn bond_properties_mut(&mut self) -> &mut [AdapterBond] {
        &mut self.bonds
    }

    pub fn conformers(&self) -> &[Vec<[f64; 3]>] {
        &self.conformers
    }

    pub fn conformers_mut(&mut self) -> impl Iterator<Item = &mut [[f64; 3]]> {
        self.conformers.iter_mut().map(std::vec::Vec::as_mut_slice)
    }

    pub fn replace_graph(
        &mut self,
        atoms: Vec<AdapterAtom>,
        bonds: Vec<AdapterBond>,
        conformers: Vec<Vec<[f64; 3]>>,
    ) -> Result<(), AdapterGraphError> {
        *self = Self::try_from_graph(atoms, bonds, conformers)?;
        Ok(())
    }

    pub(crate) fn from_graph(atoms: Vec<AdapterAtom>, bonds: Vec<AdapterBond>) -> Self {
        Self::try_from_graph(atoms, bonds, Vec::new())
            .expect("internal adapter fixture graph must have valid bond endpoints")
    }

    fn checked_bond_index(&self, index: i32) -> Result<usize, BondIndexRangeError> {
        // RDKit's getBondWithIdx takes unsigned int, so GCC/Linux converts a
        // negative pair member modulo 2^32 before URANGE_CHECK observes it.
        let index = index as u32;
        let upper_bound = self.bonds.len() as u32;
        if index >= upper_bound {
            return Err(BondIndexRangeError {
                kind: "Range Error",
                expression: "idx",
                index,
                upper_bound,
            });
        }
        Ok(index as usize)
    }
}

pub(crate) fn assign_bond_dirs(
    mol: &mut AdapterMol,
    z_bond_pairs: &[(i32, i32)],
    e_bond_pairs: &[(i32, i32)],
) -> Result<bool, BondIndexRangeError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:91 assignBondDirs
    // RDKit✔️✔️: bool assignBondDirs(RWMol &mol, INT_PAIR_VECT &zBondPairs,
    // RDKit✔️✔️:                     INT_PAIR_VECT &eBondPairs) {
    // RDKit✔️✔️:   // bonds to assign
    // RDKit✔️✔️:   std::set<int> pending;
    // RDKit✔️✔️:   for (const auto &pair : zBondPairs) {
    // RDKit✔️✔️:     pending.insert(pair.first);
    // RDKit✔️✔️:     pending.insert(pair.second);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto &pair : eBondPairs) {
    // RDKit✔️✔️:     pending.insert(pair.first);
    // RDKit✔️✔️:     pending.insert(pair.second);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // a queue for pending assignments
    // RDKit✔️✔️:   typedef std::queue<std::pair<int, Bond::BondDir>> ASSIGNMENTQTYPE;
    // RDKit✔️✔️:   ASSIGNMENTQTYPE queue;
    // RDKit✔️✔️:   // in a loop, modify one bond at a time, until all bonds are assigned
    // RDKit✔️✔️:   while (!pending.empty() || !queue.empty()) {
    // RDKit✔️✔️:     if (queue.empty()) {
    // RDKit✔️✔️:       // pumping one bond from pending to queue
    // RDKit✔️✔️:       queue.push(std::make_pair(*(pending.begin()), Bond::ENDUPRIGHT));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // pop one entry from queue and do the actual assignment
    // RDKit✔️✔️:       int curBondIdx;
    // RDKit✔️✔️:       Bond::BondDir dir;
    // RDKit✔️✔️:       boost::tie(curBondIdx, dir) = queue.front();
    // RDKit✔️✔️:       queue.pop();
    // RDKit✔️✔️:       Bond *bond = mol.getBondWithIdx(curBondIdx);
    // RDKit✔️✔️:       // is it assigned already?
    // RDKit✔️✔️:       if (bond->getBondDir() != Bond::NONE) {
    // RDKit✔️✔️:         // assigned. then check conflict
    // RDKit✔️✔️:         if (bond->getBondDir() != dir) {
    // RDKit✔️✔️:           // not doable
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         // assign since it's not assigned yet
    // RDKit✔️✔️:         bond->setBondDir(dir);
    // RDKit✔️✔️:         auto searchItr = pending.find(curBondIdx);
    // RDKit✔️✔️:         if (searchItr != pending.end()) {
    // RDKit✔️✔️:           pending.erase(searchItr);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // find all affecting bonds and add to queue by going thru all rules
    // RDKit✔️✔️:         Bond::BondDir otherDir =
    // RDKit✔️✔️:             dir == Bond::ENDUPRIGHT ? Bond::ENDDOWNRIGHT : Bond::ENDUPRIGHT;
    // RDKit✔️✔️:         // same routine for zBondPairs and eBondPairs
    // RDKit✔️✔️:         // use a switch _ to go through both by setting _ to 0 and then 1
    // RDKit✔️✔️:         for (int _ = 0; _ < 2; _++) {
    // RDKit✔️✔️:           INT_PAIR_VECT *_rules = _ == 0 ? &zBondPairs : &eBondPairs;
    // RDKit✔️✔️:           Bond::BondDir _dir = _ == 0 ? dir : otherDir;
    // RDKit✔️✔️:           for (const auto &pair : *_rules) {
    // RDKit✔️✔️:             int other = -1;
    // RDKit✔️✔️:             if (pair.first == curBondIdx) {
    // RDKit✔️✔️:               other = pair.second;
    // RDKit✔️✔️:             } else if (pair.second == curBondIdx) {
    // RDKit✔️✔️:               other = pair.first;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             // a match?
    // RDKit✔️✔️:             if (other != curBondIdx && other != -1) {
    // RDKit✔️✔️:               Bond *otherBond = mol.getBondWithIdx(other);
    // RDKit✔️✔️:               // check if it is assigned
    // RDKit✔️✔️:               if (otherBond->getBondDir() != Bond::NONE) {
    // RDKit✔️✔️:                 // assigned. check conflict
    // RDKit✔️✔️:                 if (otherBond->getBondDir() != _dir) {
    // RDKit✔️✔️:                   // not doable
    // RDKit✔️✔️:                   return false;
    // RDKit✔️✔️:                 }
    // RDKit✔️✔️:               } else {
    // RDKit✔️✔️:                 // not assigned, then add to queue
    // RDKit✔️✔️:                 queue.push(std::make_pair(otherBond->getIdx(), _dir));
    // RDKit✔️✔️:               }  // end if otherBond's bond direction check
    // RDKit✔️✔️:             }  // end if there is a match
    // RDKit✔️✔️:           }  // end loop over pairs in _rules
    // RDKit✔️✔️:         }  // end for _ to go thru rule sets
    // RDKit✔️✔️:       }  // end if this bond is assigned
    // RDKit✔️✔️:     }  // end if queue is empty
    // RDKit✔️✔️:   }  // end while on pending set and queue
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: assignBondDirs
    // BEGIN RDKIT ACTIVE CONFIGURATION: assignBondDirs
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: BondDir values come from Code/GraphMol/Bond.h:81-91.
    // RDKit✔️✔️: getBondWithIdx uses unsigned-int URANGE_CHECK from ROMol.cpp:317-329.
    // END RDKIT ACTIVE CONFIGURATION: assignBondDirs

    let mut pending = BTreeSet::new();
    for &(first, second) in z_bond_pairs {
        pending.insert(first);
        pending.insert(second);
    }
    for &(first, second) in e_bond_pairs {
        pending.insert(first);
        pending.insert(second);
    }

    let mut queue = VecDeque::new();
    while !pending.is_empty() || !queue.is_empty() {
        if queue.is_empty() {
            queue.push_back((
                *pending.first().expect("pending is nonempty"),
                BondDirection::EndUpRight,
            ));
        } else {
            let (current_bond_index, direction) = queue.pop_front().expect("queue is nonempty");
            let current_position = mol.checked_bond_index(current_bond_index)?;
            let current_direction = mol.bonds[current_position].direction;

            if current_direction != BondDirection::None {
                if current_direction != direction {
                    return Ok(false);
                }
            } else {
                mol.bonds[current_position].direction = direction;
                pending.remove(&current_bond_index);

                let other_direction = if direction == BondDirection::EndUpRight {
                    BondDirection::EndDownRight
                } else {
                    BondDirection::EndUpRight
                };

                for (rules, rule_direction) in [(z_bond_pairs, direction), (e_bond_pairs, other_direction)] {
                    for &(first, second) in rules {
                        let other = if first == current_bond_index {
                            second
                        } else if second == current_bond_index {
                            first
                        } else {
                            -1
                        };

                        if other != current_bond_index && other != -1 {
                            let other_position = mol.checked_bond_index(other)?;
                            let other_bond_direction = mol.bonds[other_position].direction;
                            if other_bond_direction != BondDirection::None {
                                if other_bond_direction != rule_direction {
                                    return Ok(false);
                                }
                            } else {
                                queue.push_back((other_position as i32, rule_direction));
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(true)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn find_alternating_bonds(
    mol: &AdapterMol,
    current_atom_index: u32,
    desired_atomic_number: i32,
    desired_atom_charge: i32,
    desired_next_bond_type: BondType,
    desired_ending_bond_type: BondType,
    current_path_length: u32,
    max_path_length: u32,
    last_bond_index: Option<u32>,
    path: &mut Vec<u32>,
    visited: &mut BTreeSet<i32>,
) -> Option<u32> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:207 findAlternatingBonds
    // RDKit✔️✔️: Atom *findAlternatingBonds(
    // RDKit✔️✔️:     ROMol &mol, Atom *current, int desiredAtomicNumber, int desiredAtomCharge,
    // RDKit✔️✔️:     Bond::BondType desiredNextBondType, Bond::BondType desiredEndingBondType,
    // RDKit✔️✔️:     unsigned int currentPathLength, unsigned int maxPathLength, Bond *lastBond,
    // RDKit✔️✔️:     /*OUT*/ std::stack<Bond *> &path, std::set<int> &_visited) {
    // RDKit✔️✔️:   // memory for what has been visited
    // RDKit✔️✔️:   if (lastBond == nullptr) {
    // RDKit✔️✔️:     _visited.clear();
    // RDKit✔️✔️:     while (!path.empty()) {
    // RDKit✔️✔️:       path.pop();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   _visited.insert(current->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // for (int i = 0; i < currentPathLength; i ++)
    // RDKit✔️✔️:   //  std::cerr << ".";
    // RDKit✔️✔️:   // std::cerr << (int) current->getIdx() << "("
    // RDKit✔️✔️:   //  << (int) current->getAtomicNum()
    // RDKit✔️✔️:   //  << ")" << std::endl;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // is this atom the desired one?
    // RDKit✔️✔️:   if (lastBond && current->getAtomicNum() == desiredAtomicNumber &&
    // RDKit✔️✔️:       lastBond->getBondType() == desiredEndingBondType &&
    // RDKit✔️✔️:       current->getFormalCharge() == desiredAtomCharge) {
    // RDKit✔️✔️:     // Yes! But am I better than the existing one - if one exists?
    // RDKit✔️✔️:     if (path.size() == 0 || path.size() > currentPathLength) {
    // RDKit✔️✔️:       // Yes! clear the path and repopulate it
    // RDKit✔️✔️:       while (!path.empty()) {
    // RDKit✔️✔️:         path.pop();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // add myself to the path
    // RDKit✔️✔️:       path.push(lastBond);
    // RDKit✔️✔️:       return current;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // I am no better than the existing one. This will also cause the
    // RDKit✔️✔️:       // path search to not continue down
    // RDKit✔️✔️:       return nullptr;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // searching too far?
    // RDKit✔️✔️:   if (maxPathLength <= currentPathLength) {
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // continue searching down
    // RDKit✔️✔️:   RWMol::ADJ_ITER nid, end;
    // RDKit✔️✔️:   Atom *target = nullptr, *temp;
    // RDKit✔️✔️:   for (boost::tie(nid, end) = mol.getAtomNeighbors(current); nid != end;
    // RDKit✔️✔️:        nid++) {
    // RDKit✔️✔️:     if (_visited.find(*nid) != _visited.end()) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // check whether bond is valid for search to go down through it
    // RDKit✔️✔️:     Bond *bond = mol.getBondBetweenAtoms(current->getIdx(), *nid);
    // RDKit✔️✔️:     if (bond->getBondType() == desiredNextBondType) {
    // RDKit✔️✔️:       // recursive call: for all ways to extend the path, ask each to try
    // RDKit✔️✔️:       // enhancing the current best path (stored in <path>)
    // RDKit✔️✔️:       // by setting SINGLE as the default, we allow a very special case to
    // RDKit✔️✔️:       // be supported: a TRIPLE bond followed by a SINGLE bond
    // RDKit✔️✔️:       // This is used in _Valence5NCleanUp2
    // RDKit✔️✔️:       Bond::BondType nextBondType = Bond::SINGLE;
    // RDKit✔️✔️:       if (desiredNextBondType == Bond::SINGLE) {
    // RDKit✔️✔️:         nextBondType = Bond::DOUBLE;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if ((temp = findAlternatingBonds(
    // RDKit✔️✔️:                mol, mol.getAtomWithIdx(*nid), desiredAtomicNumber,
    // RDKit✔️✔️:                desiredAtomCharge, nextBondType, desiredEndingBondType,
    // RDKit✔️✔️:                currentPathLength + 1, maxPathLength, bond, path, _visited)) !=
    // RDKit✔️✔️:           nullptr) {
    // RDKit✔️✔️:         target = temp;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (desiredEndingBondType != Bond::SINGLE &&
    // RDKit✔️✔️:                desiredEndingBondType != Bond::DOUBLE &&
    // RDKit✔️✔️:                bond->getBondType() == desiredEndingBondType) {
    // RDKit✔️✔️:       // try recursive call limited to one level down to see whether
    // RDKit✔️✔️:       // this can serve as the last leg of the path. This is done only if
    // RDKit✔️✔️:       // the desiredEndingBondType is not part of the alternating bonds
    // RDKit✔️✔️:       if ((temp = findAlternatingBonds(
    // RDKit✔️✔️:                mol, mol.getAtomWithIdx(*nid), desiredAtomicNumber,
    // RDKit✔️✔️:                desiredAtomCharge, Bond::UNSPECIFIED, /* no next */
    // RDKit✔️✔️:                desiredEndingBondType, currentPathLength + 1,
    // RDKit✔️✔️:                0, /* this limits the recursion */
    // RDKit✔️✔️:                bond, path, _visited))) {
    // RDKit✔️✔️:         target = temp;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // about the return
    // RDKit✔️✔️:   if (target != nullptr) {
    // RDKit✔️✔️:     if (lastBond) {
    // RDKit✔️✔️:       path.push(lastBond);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return target;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: findAlternatingBonds
    // BEGIN RDKIT ACTIVE CONFIGURATION: findAlternatingBonds
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: MolGraph is boost::adjacency_list<vecS, vecS, undirectedS>, so neighbor
    // RDKit✔️✔️: iteration follows each vertex's stored out-edge order.
    // END RDKIT ACTIVE CONFIGURATION: findAlternatingBonds

    if last_bond_index.is_none() {
        visited.clear();
        path.clear();
    }
    let current = &mol.atoms[current_atom_index as usize];
    visited.insert(current_atom_index as i32);

    if let Some(last_bond_index) = last_bond_index {
        let last_bond = &mol.bonds[last_bond_index as usize];
        if current.atomic_number == desired_atomic_number
            && last_bond.bond_type == desired_ending_bond_type
            && current.formal_charge == desired_atom_charge
        {
            if path.is_empty() || path.len() > current_path_length as usize {
                path.clear();
                path.push(last_bond_index);
                return Some(current_atom_index);
            }
            return None;
        }
    }

    if max_path_length <= current_path_length {
        return None;
    }

    let mut target = None;
    for &(neighbor_index, bond_index) in &mol.adjacency[current_atom_index as usize] {
        if visited.contains(&(neighbor_index as i32)) {
            continue;
        }
        let bond = &mol.bonds[bond_index as usize];
        if bond.bond_type == desired_next_bond_type {
            let next_bond_type = if desired_next_bond_type == BondType::Single {
                BondType::Double
            } else {
                BondType::Single
            };
            if let Some(found) = find_alternating_bonds(
                mol,
                neighbor_index,
                desired_atomic_number,
                desired_atom_charge,
                next_bond_type,
                desired_ending_bond_type,
                current_path_length.wrapping_add(1),
                max_path_length,
                Some(bond_index),
                path,
                visited,
            ) {
                target = Some(found);
            }
        } else if desired_ending_bond_type != BondType::Single
            && desired_ending_bond_type != BondType::Double
            && bond.bond_type == desired_ending_bond_type
            && let Some(found) = find_alternating_bonds(
                mol,
                neighbor_index,
                desired_atomic_number,
                desired_atom_charge,
                BondType::Unspecified,
                desired_ending_bond_type,
                current_path_length.wrapping_add(1),
                0,
                Some(bond_index),
                path,
                visited,
            )
        {
            target = Some(found);
        }
    }

    if let Some(target) = target {
        if let Some(last_bond_index) = last_bond_index {
            path.push(last_bond_index);
        }
        return Some(target);
    }
    None
}

pub(crate) fn get_num_double_bonded_negatively_charged_neighboring_si(mol: &AdapterMol, atom_index: u32) -> i32 {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:306 getNumDoubleBondedNegativelyChargedNeighboringSi
    // RDKit✔️✔️: int getNumDoubleBondedNegativelyChargedNeighboringSi(ROMol &mol, Atom *a) {
    // RDKit✔️✔️:   RWMol::ADJ_ITER nid1, end1;
    // RDKit✔️✔️:   boost::tie(nid1, end1) = mol.getAtomNeighbors(a);
    // RDKit✔️✔️:   int nSi = 0;
    // RDKit✔️✔️:   int thisId = a->getIdx();
    // RDKit✔️✔️:   while (nid1 != end1) {
    // RDKit✔️✔️:     Atom *nbr = mol.getAtomWithIdx(*nid1);
    // RDKit✔️✔️:     Bond *bond = mol.getBondBetweenAtoms(*nid1, thisId);
    // RDKit✔️✔️:     if (nbr->getAtomicNum() == 14 && nbr->getFormalCharge() == -1 &&
    // RDKit✔️✔️:         bond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       nSi++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     nid1++;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nSi;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: getNumDoubleBondedNegativelyChargedNeighboringSi
    // BEGIN RDKIT ACTIVE CONFIGURATION: getNumDoubleBondedNegativelyChargedNeighboringSi
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: Valid ROMol adjacency contains one entry per incident bond and the source
    // RDKit✔️✔️: count remains representable as int for the modeled in-memory graph.
    // END RDKIT ACTIVE CONFIGURATION: getNumDoubleBondedNegativelyChargedNeighboringSi

    let mut silicon_count = 0_i32;
    for &(neighbor_index, bond_index) in &mol.adjacency[atom_index as usize] {
        let neighbor = &mol.atoms[neighbor_index as usize];
        let bond = &mol.bonds[bond_index as usize];
        if neighbor.atomic_number == 14 && neighbor.formal_charge == -1 && bond.bond_type == BondType::Double {
            silicon_count += 1;
        }
    }
    silicon_count
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AdapterValenceError {
    pub(crate) kind: &'static str,
    pub(crate) message: &'static str,
    pub(crate) bond_index: u32,
    pub(crate) bond_type: BondType,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum AdapterCleanup5Error {
    Precondition {
        kind: &'static str,
        expression: &'static str,
        message: &'static str,
    },
    Valence(AdapterValenceError),
}

impl From<AdapterValenceError> for AdapterCleanup5Error {
    fn from(error: AdapterValenceError) -> Self {
        Self::Valence(error)
    }
}

fn cleanup_explicit_valence(mol: &mut AdapterMol, atom_index: u32) -> Result<i32, AdapterValenceError> {
    let mut valence = f64::from(mol.atoms[atom_index as usize].num_explicit_hydrogens);
    for &(_, bond_index) in &mol.adjacency[atom_index as usize] {
        let bond = &mol.bonds[bond_index as usize];
        let contribution = match bond.bond_type {
            BondType::Unspecified | BondType::Ionic | BondType::Hydrogen | BondType::Zero => 0.0,
            BondType::Single => 1.0,
            BondType::Double => 2.0,
            BondType::Triple => 3.0,
            BondType::Quadruple => 4.0,
            BondType::Quintuple => 5.0,
            BondType::Hextuple => 6.0,
            BondType::OneAndAHalf | BondType::Aromatic => 1.5,
            BondType::TwoAndAHalf => 2.5,
            BondType::ThreeAndAHalf => 3.5,
            BondType::FourAndAHalf => 4.5,
            BondType::FiveAndAHalf => 5.5,
            BondType::DativeOne | BondType::Dative => {
                if bond.end_atom_index == atom_index {
                    1.0
                } else {
                    0.0
                }
            }
            BondType::ThreeCenter | BondType::DativeL | BondType::DativeR | BondType::Other => {
                return Err(AdapterValenceError {
                    kind: "ValueErrorException",
                    message: "Bad bond type",
                    bond_index,
                    bond_type: bond.bond_type,
                });
            }
        };
        valence += contribution;
    }
    let valence = (valence + 0.1).round() as i32;
    // RDKit✔️✔️: atom->calcExplicitValence(false);
    // `calcExplicitValence()` writes `d_explicitValence`. Some clean-up
    // branches intentionally do not call it again after changing a bond, so
    // retaining this value is observable when sanitization is disabled.
    mol.atoms[atom_index as usize].cached_explicit_valence = Some(valence);
    Ok(valence)
}

fn bond_index_between(mol: &AdapterMol, first: u32, second: u32) -> Option<u32> {
    mol.adjacency[first as usize]
        .iter()
        .find_map(|&(neighbor, bond_index)| (neighbor == second).then_some(bond_index))
}

fn valence4n_cleanup1_matches(mol: &AdapterMol) -> Vec<[u32; 5]> {
    let mut seen_atom_sets = HashSet::<[u32; 5]>::new();
    let mut matches = Vec::<[u32; 5]>::new();

    for query_2 in 0..mol.atoms.len() as u32 {
        if mol.atoms[query_2 as usize].atomic_number != 50 {
            continue;
        }
        for &(query_1, bond_12) in &mol.adjacency[query_2 as usize] {
            if mol.atoms[query_1 as usize].atomic_number != 7
                || mol.bonds[bond_12 as usize].bond_type != BondType::Double
            {
                continue;
            }
            for &(query_3, bond_23) in &mol.adjacency[query_2 as usize] {
                if query_3 == query_1
                    || mol.atoms[query_3 as usize].atomic_number != 7
                    || mol.bonds[bond_23 as usize].bond_type != BondType::Double
                {
                    continue;
                }
                for &(query_0, bond_01) in &mol.adjacency[query_1 as usize] {
                    if query_0 == query_2
                        || query_0 == query_3
                        || mol.atoms[query_0 as usize].atomic_number != 6
                        || mol.bonds[bond_01 as usize].bond_type != BondType::Single
                    {
                        continue;
                    }
                    for &(query_4, bond_34) in &mol.adjacency[query_3 as usize] {
                        if query_4 == query_0
                            || query_4 == query_1
                            || query_4 == query_2
                            || mol.atoms[query_4 as usize].atomic_number != 7
                            || mol.bonds[bond_34 as usize].bond_type != BondType::Single
                        {
                            continue;
                        }
                        let Some(bond_40) = bond_index_between(mol, query_4, query_0) else {
                            continue;
                        };
                        if mol.bonds[bond_40 as usize].bond_type != BondType::Double {
                            continue;
                        }

                        let mapping = [query_0, query_1, query_2, query_3, query_4];
                        let mut atom_set = mapping;
                        atom_set.sort_unstable();
                        if seen_atom_sets.insert(atom_set) {
                            matches.push(mapping);
                            if matches.len() == 1000 {
                                return matches;
                            }
                        }
                    }
                }
            }
        }
    }

    matches
}

fn valence5n_cleanup6_matches(mol: &AdapterMol) -> Vec<[u32; 7]> {
    let mut seen_atom_sets = HashSet::<[u32; 7]>::new();
    let mut matches = Vec::<[u32; 7]>::new();

    for query_2 in 0..mol.atoms.len() as u32 {
        if mol.atoms[query_2 as usize].atomic_number != 50 {
            continue;
        }
        for &(query_1, bond_12) in &mol.adjacency[query_2 as usize] {
            if mol.atoms[query_1 as usize].atomic_number != 6
                || mol.bonds[bond_12 as usize].bond_type != BondType::Double
            {
                continue;
            }
            for &(query_3, bond_23) in &mol.adjacency[query_2 as usize] {
                if query_3 == query_1
                    || mol.atoms[query_3 as usize].atomic_number != 6
                    || mol.bonds[bond_23 as usize].bond_type != BondType::Double
                {
                    continue;
                }
                for &(query_6, bond_26) in &mol.adjacency[query_2 as usize] {
                    if query_6 == query_1
                        || query_6 == query_3
                        || mol.atoms[query_6 as usize].atomic_number != 6
                        || mol.bonds[bond_26 as usize].bond_type != BondType::Single
                    {
                        continue;
                    }
                    for &(query_0, bond_01) in &mol.adjacency[query_1 as usize] {
                        if [query_2, query_3, query_6].contains(&query_0)
                            || mol.atoms[query_0 as usize].atomic_number != 6
                            || mol.bonds[bond_01 as usize].bond_type != BondType::Single
                        {
                            continue;
                        }
                        for &(query_4, _) in &mol.adjacency[query_3 as usize] {
                            if [query_0, query_1, query_2, query_6].contains(&query_4)
                                || mol.atoms[query_4 as usize].atomic_number != 6
                            {
                                continue;
                            }
                            for &(query_5, bond_45) in &mol.adjacency[query_4 as usize] {
                                if [query_0, query_1, query_2, query_3, query_6].contains(&query_5)
                                    || mol.atoms[query_5 as usize].atomic_number != 7
                                    || mol.bonds[bond_45 as usize].bond_type != BondType::Single
                                {
                                    continue;
                                }
                                let Some(bond_50) = bond_index_between(mol, query_5, query_0) else {
                                    continue;
                                };
                                if mol.bonds[bond_50 as usize].bond_type != BondType::Double {
                                    continue;
                                }

                                let mapping = [query_0, query_1, query_2, query_3, query_4, query_5, query_6];
                                let mut atom_set = mapping;
                                atom_set.sort_unstable();
                                if seen_atom_sets.insert(atom_set) {
                                    matches.push(mapping);
                                    if matches.len() == 1000 {
                                        return matches;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    matches
}

fn valence5n_cleanup7_matches(mol: &AdapterMol) -> Vec<[u32; 7]> {
    let mut seen_atom_sets = HashSet::<[u32; 7]>::new();
    let mut matches = Vec::<[u32; 7]>::new();

    for query_2 in 0..mol.atoms.len() as u32 {
        if mol.atoms[query_2 as usize].atomic_number != 50 {
            continue;
        }
        for &(query_1, bond_12) in &mol.adjacency[query_2 as usize] {
            if mol.atoms[query_1 as usize].atomic_number != 6
                || mol.bonds[bond_12 as usize].bond_type != BondType::Double
            {
                continue;
            }
            for &(query_3, bond_23) in &mol.adjacency[query_2 as usize] {
                if query_3 == query_1
                    || mol.atoms[query_3 as usize].atomic_number != 7
                    || mol.bonds[bond_23 as usize].bond_type != BondType::Double
                {
                    continue;
                }
                for &(query_6, bond_26) in &mol.adjacency[query_2 as usize] {
                    if [query_1, query_3].contains(&query_6)
                        || mol.atoms[query_6 as usize].atomic_number != 6
                        || mol.bonds[bond_26 as usize].bond_type != BondType::Single
                    {
                        continue;
                    }
                    for &(query_0, _) in &mol.adjacency[query_1 as usize] {
                        if [query_2, query_3, query_6].contains(&query_0)
                            || mol.atoms[query_0 as usize].atomic_number != 6
                        {
                            continue;
                        }
                        for &(query_4, bond_34) in &mol.adjacency[query_3 as usize] {
                            if [query_0, query_1, query_2, query_6].contains(&query_4)
                                || mol.atoms[query_4 as usize].atomic_number != 6
                                || mol.bonds[bond_34 as usize].bond_type != BondType::Single
                            {
                                continue;
                            }
                            for &(query_5, bond_45) in &mol.adjacency[query_4 as usize] {
                                if [query_0, query_1, query_2, query_3, query_6].contains(&query_5)
                                    || mol.atoms[query_5 as usize].atomic_number != 8
                                    || mol.bonds[bond_45 as usize].bond_type != BondType::Single
                                {
                                    continue;
                                }
                                let Some(bond_50) = bond_index_between(mol, query_5, query_0) else {
                                    continue;
                                };
                                if mol.bonds[bond_50 as usize].bond_type != BondType::Single {
                                    continue;
                                }

                                let mapping = [query_0, query_1, query_2, query_3, query_4, query_5, query_6];
                                let mut atom_set = mapping;
                                atom_set.sort_unstable();
                                if seen_atom_sets.insert(atom_set) {
                                    matches.push(mapping);
                                    if matches.len() == 1000 {
                                        return matches;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    matches
}

fn valence5n_cleanup8_matches(mol: &AdapterMol) -> Vec<[u32; 6]> {
    let mut seen_atom_sets = HashSet::<[u32; 6]>::new();
    let mut matches = Vec::<[u32; 6]>::new();
    for query_5 in 0..mol.atoms.len() as u32 {
        if mol.atoms[query_5 as usize].atomic_number != 50 {
            continue;
        }
        for &(query_0, bond_50) in &mol.adjacency[query_5 as usize] {
            if mol.atoms[query_0 as usize].atomic_number != 6
                || mol.bonds[bond_50 as usize].bond_type != BondType::Double
            {
                continue;
            }
            for &(query_1, bond_01) in &mol.adjacency[query_0 as usize] {
                if query_1 == query_5
                    || mol.atoms[query_1 as usize].atomic_number != 7
                    || mol.bonds[bond_01 as usize].bond_type != BondType::Single
                {
                    continue;
                }
                for &(query_4, bond_40) in &mol.adjacency[query_0 as usize] {
                    if [query_1, query_5].contains(&query_4)
                        || mol.atoms[query_4 as usize].atomic_number != 7
                        || mol.bonds[bond_40 as usize].bond_type != BondType::Single
                    {
                        continue;
                    }
                    for &(query_2, bond_12) in &mol.adjacency[query_1 as usize] {
                        if [query_0, query_4, query_5].contains(&query_2)
                            || mol.atoms[query_2 as usize].atomic_number != 6
                            || mol.bonds[bond_12 as usize].bond_type != BondType::Double
                        {
                            continue;
                        }
                        for &(query_3, bond_23) in &mol.adjacency[query_2 as usize] {
                            if [query_0, query_1, query_4, query_5].contains(&query_3)
                                || mol.atoms[query_3 as usize].atomic_number != 7
                                || mol.bonds[bond_23 as usize].bond_type != BondType::Single
                            {
                                continue;
                            }
                            let Some(bond_34) = bond_index_between(mol, query_3, query_4) else {
                                continue;
                            };
                            if mol.bonds[bond_34 as usize].bond_type != BondType::Double {
                                continue;
                            }
                            let mapping = [query_0, query_1, query_2, query_3, query_4, query_5];
                            let mut atom_set = mapping;
                            atom_set.sort_unstable();
                            if seen_atom_sets.insert(atom_set) {
                                matches.push(mapping);
                                if matches.len() == 1000 {
                                    return matches;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    matches
}

fn valence5n_cleanup9_matches(mol: &AdapterMol) -> Vec<[u32; 6]> {
    let mut seen_atom_sets = HashSet::<[u32; 6]>::new();
    let mut matches = Vec::<[u32; 6]>::new();
    for query_5 in 0..mol.atoms.len() as u32 {
        if mol.atoms[query_5 as usize].atomic_number != 50 {
            continue;
        }
        for &(query_0, bond_50) in &mol.adjacency[query_5 as usize] {
            if mol.atoms[query_0 as usize].atomic_number != 6
                || mol.bonds[bond_50 as usize].bond_type != BondType::Double
            {
                continue;
            }
            for &(query_1, bond_01) in &mol.adjacency[query_0 as usize] {
                if query_1 == query_5
                    || mol.atoms[query_1 as usize].atomic_number != 7
                    || mol.bonds[bond_01 as usize].bond_type != BondType::Single
                {
                    continue;
                }
                for &(query_4, bond_40) in &mol.adjacency[query_0 as usize] {
                    if [query_1, query_5].contains(&query_4)
                        || mol.atoms[query_4 as usize].atomic_number != 6
                        || mol.bonds[bond_40 as usize].bond_type != BondType::Single
                    {
                        continue;
                    }
                    for &(query_2, bond_12) in &mol.adjacency[query_1 as usize] {
                        if [query_0, query_4, query_5].contains(&query_2)
                            || mol.atoms[query_2 as usize].atomic_number != 7
                            || mol.bonds[bond_12 as usize].bond_type != BondType::Double
                        {
                            continue;
                        }
                        for &(query_3, bond_23) in &mol.adjacency[query_2 as usize] {
                            if [query_0, query_1, query_4, query_5].contains(&query_3)
                                || mol.atoms[query_3 as usize].atomic_number != 6
                                || mol.bonds[bond_23 as usize].bond_type != BondType::Single
                            {
                                continue;
                            }
                            let Some(bond_34) = bond_index_between(mol, query_3, query_4) else {
                                continue;
                            };
                            if mol.bonds[bond_34 as usize].bond_type != BondType::Double {
                                continue;
                            }
                            let mapping = [query_0, query_1, query_2, query_3, query_4, query_5];
                            let mut atom_set = mapping;
                            atom_set.sort_unstable();
                            if seen_atom_sets.insert(atom_set) {
                                matches.push(mapping);
                                if matches.len() == 1000 {
                                    return matches;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    matches
}

fn valence5n_cleanupa_matches(mol: &AdapterMol) -> Vec<[u32; 2]> {
    let mut seen_atom_sets = HashSet::<[u32; 2]>::new();
    let mut matches = Vec::<[u32; 2]>::new();
    for query_0 in 0..mol.atoms.len() as u32 {
        if mol.atoms[query_0 as usize].atomic_number != 7 {
            continue;
        }
        for &(query_1, bond_01) in &mol.adjacency[query_0 as usize] {
            if mol.atoms[query_1 as usize].atomic_number != 7
                || mol.bonds[bond_01 as usize].bond_type != BondType::Double
            {
                continue;
            }
            let mapping = [query_0, query_1];
            let mut atom_set = mapping;
            atom_set.sort_unstable();
            if seen_atom_sets.insert(atom_set) {
                matches.push(mapping);
                if matches.len() == 1000 {
                    return matches;
                }
            }
        }
    }
    matches
}

pub(crate) fn valence4n_cleanup1(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:324 _Valence4NCleanUp1
    // RDKit✔️✔️: bool _Valence4NCleanUp1(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   // replace the N- with Sn
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != -1 ||
    // RDKit✔️✔️:       atom->calcExplicitValence(false) != 4) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   atom->setAtomicNum(50);
    // RDKit✔️✔️:   atom->setFormalCharge(0);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // substructure matching
    // RDKit✔️✔️:   auto *query = new RWMol();
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 0
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);   // 1
    // RDKit✔️✔️:   query->addAtom(new Atom(50), false, true);  // 2
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);   // 3
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);   // 4
    // RDKit✔️✔️:   query->addBond(0, 1, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(1, 2, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(2, 3, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(3, 4, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(4, 0, Bond::DOUBLE);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<MatchVectType> fgpMatches;
    // RDKit✔️✔️:   SubstructMatch(mol, *query, fgpMatches);
    // RDKit✔️✔️:   delete query;
    // RDKit✔️✔️:   // no action if none or more than one match was found
    // RDKit✔️✔️:   if (fgpMatches.size() != 1) {
    // RDKit✔️✔️:     atom->setAtomicNum(7);
    // RDKit✔️✔️:     atom->setFormalCharge(-1);
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // collect matching atoms
    // RDKit✔️✔️:   int map[5];
    // RDKit✔️✔️:   MatchVectType match = fgpMatches[0];
    // RDKit✔️✔️:   for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
    // RDKit✔️✔️:        mi++) {
    // RDKit✔️✔️:     map[mi->first] = mi->second;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // flip bonds
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[0], map[1])->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[1], map[2])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[2], map[3])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[3], map[4])->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[4], map[0])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   // change the problematic N-
    // RDKit✔️✔️:   atom->setAtomicNum(7);
    // RDKit✔️✔️:   atom->setFormalCharge(-1);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence4NCleanUp1
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence4NCleanUp1
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: The selected InchiToMol caller sets explicit H and noImplicit, leaves atom
    // RDKit✔️✔️: aromatic flags false before cleanUp, and creates single, double, triple, or
    // RDKit✔️✔️: aromatic bonds. SubstructMatch uses default uniquify-by-atom-set behavior.
    // END RDKIT ACTIVE CONFIGURATION: _Valence4NCleanUp1

    let atom_position = atom_index as usize;
    if mol.atoms[atom_position].atomic_number != 7
        || mol.atoms[atom_position].formal_charge != -1
        || cleanup_explicit_valence(mol, atom_index)? != 4
    {
        return Ok(false);
    }

    mol.atoms[atom_position].atomic_number = 50;
    mol.atoms[atom_position].formal_charge = 0;
    let matches = valence4n_cleanup1_matches(mol);
    if matches.len() != 1 {
        mol.atoms[atom_position].atomic_number = 7;
        mol.atoms[atom_position].formal_charge = -1;
        return Ok(false);
    }

    let mapping = matches[0];
    for (first, second, bond_type) in [
        (mapping[0], mapping[1], BondType::Double),
        (mapping[1], mapping[2], BondType::Single),
        (mapping[2], mapping[3], BondType::Single),
        (mapping[3], mapping[4], BondType::Double),
        (mapping[4], mapping[0], BondType::Single),
    ] {
        let bond_index =
            bond_index_between(mol, first, second).expect("the selected substructure match contains every query bond");
        mol.bonds[bond_index as usize].bond_type = bond_type;
    }
    mol.atoms[atom_position].atomic_number = 7;
    mol.atoms[atom_position].formal_charge = -1;
    Ok(true)
}

pub(crate) fn valence4n_cleanup2(mol: &mut AdapterMol, atom_index: u32) -> bool {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:376 _Valence4NCleanUp2
    // RDKit✔️✔️: bool _Valence4NCleanUp2(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 7, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 1,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:   if (target == nullptr) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   stack.top()->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   atom->setFormalCharge(0);
    // RDKit✔️✔️:   target->setFormalCharge(-1);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence4NCleanUp2
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence4NCleanUp2
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: The selected cleanUp caller supplies a valid atom pointer in a valid RWMol.
    // RDKit✔️✔️: findAlternatingBonds is the source-exact adapter helper closed above; with
    // RDKit✔️✔️: maxPathLength 1 this call can select only a directly double-bonded neutral N.
    // END RDKIT ACTIVE CONFIGURATION: _Valence4NCleanUp2

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    let Some(target_index) = find_alternating_bonds(
        mol,
        atom_index,
        7,
        0,
        BondType::Double,
        BondType::Double,
        0,
        1,
        None,
        &mut path,
        &mut visited,
    ) else {
        return false;
    };

    let selected_bond_index = *path.last().expect("a found target has a path bond");
    mol.bonds[selected_bond_index as usize].bond_type = BondType::Single;
    mol.atoms[atom_index as usize].formal_charge = 0;
    mol.atoms[target_index as usize].formal_charge = -1;
    true
}

pub(crate) fn valence5n_cleanup1(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:393 _Valence5NCleanUp1
    // RDKit✔️✔️: bool _Valence5NCleanUp1(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 7, 1, Bond::DOUBLE, Bond::DOUBLE, 0, 5,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:   if (target == nullptr) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   target->setFormalCharge(0);
    // RDKit✔️✔️:   target->calcExplicitValence(false);
    // RDKit✔️✔️:   while (!stack.empty()) {
    // RDKit✔️✔️:     if (stack.top()->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       stack.top()->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       stack.top()->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     stack.pop();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   atom->setFormalCharge(1);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUp1
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp1
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: The selected cleanUp caller supplies a valid atom pointer in a valid RWMol;
    // RDKit✔️✔️: findAlternatingBonds and calcExplicitValence(false) use the source-backed
    // RDKit✔️✔️: adapter behaviors above, including adjacency order and bond-type exceptions.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp1

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    let Some(target_index) = find_alternating_bonds(
        mol,
        atom_index,
        7,
        1,
        BondType::Double,
        BondType::Double,
        0,
        5,
        None,
        &mut path,
        &mut visited,
    ) else {
        return Ok(false);
    };

    mol.atoms[target_index as usize].formal_charge = 0;
    cleanup_explicit_valence(mol, target_index)?;
    while let Some(bond_index) = path.pop() {
        let bond = &mut mol.bonds[bond_index as usize];
        bond.bond_type = if bond.bond_type == BondType::Double {
            BondType::Single
        } else {
            BondType::Double
        };
    }
    mol.atoms[atom_index as usize].formal_charge = 1;
    Ok(true)
}

pub(crate) fn valence5n_cleanup2(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:417 _Valence5NCleanUp2
    // RDKit✔️✔️: bool _Valence5NCleanUp2(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 7, -1, Bond::TRIPLE, Bond::SINGLE, 0, 2,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:   if (target == nullptr) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Bond *bond = stack.top();
    // RDKit✔️✔️:   bond->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   if (bond->getBeginAtomIdx() == atom->getIdx()) {
    // RDKit✔️✔️:     mol.getAtomWithIdx(bond->getEndAtomIdx())->setFormalCharge(-1);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     mol.getAtomWithIdx(bond->getBeginAtomIdx())->setFormalCharge(-1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   stack.pop();
    // RDKit✔️✔️:   stack.top()->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:   target->setFormalCharge(0);
    // RDKit✔️✔️:   target->calcExplicitValence(false);
    // RDKit✔️✔️:   atom->calcExplicitValence(false);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUp2
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp2
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: The selected cleanUp caller supplies a valid atom pointer in a valid RWMol;
    // RDKit✔️✔️: findAlternatingBonds returns the two-bond stack in source stack order for the
    // RDKit✔️✔️: required TRIPLE-then-SINGLE path, and calcExplicitValence(false) can throw.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp2

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    let Some(target_index) = find_alternating_bonds(
        mol,
        atom_index,
        7,
        -1,
        BondType::Triple,
        BondType::Single,
        0,
        2,
        None,
        &mut path,
        &mut visited,
    ) else {
        return Ok(false);
    };

    let root_bond_index = *path.last().expect("a found target has a root-side bond");
    mol.bonds[root_bond_index as usize].bond_type = BondType::Single;
    let root_bond = &mol.bonds[root_bond_index as usize];
    let intermediate_atom_index = if root_bond.begin_atom_index == atom_index {
        root_bond.end_atom_index
    } else {
        root_bond.begin_atom_index
    };
    mol.atoms[intermediate_atom_index as usize].formal_charge = -1;

    path.pop();
    let target_bond_index = *path.last().expect("a found target has a target-side bond");
    mol.bonds[target_bond_index as usize].bond_type = BondType::Double;
    mol.atoms[target_index as usize].formal_charge = 0;
    cleanup_explicit_valence(mol, target_index)?;
    cleanup_explicit_valence(mol, atom_index)?;
    Ok(true)
}

pub(crate) fn valence5n_cleanup3(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:443 _Valence5NCleanUp3
    // RDKit✔️✔️: bool _Valence5NCleanUp3(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 7, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 1,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:   if (target == nullptr) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // we are double bonded to a neighboring N. Check to see if we are also
    // RDKit✔️✔️:   // double bonded to an O. If so, we don't want to mess with the other N
    // RDKit✔️✔️:   // this occurs because the InChI code produces this structure:
    // RDKit✔️✔️:   //  CN(=O)=N(=O)C
    // RDKit✔️✔️:   // and we don't want to mess with that.
    // RDKit✔️✔️:   // this was github #1572
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::stack<Bond *> stack2;
    // RDKit✔️✔️:   std::set<int> _visited2;
    // RDKit✔️✔️:   Atom *target2 =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 8, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 1,
    // RDKit✔️✔️:                            nullptr, stack2, _visited2);
    // RDKit✔️✔️:   if (target2 == nullptr) {
    // RDKit✔️✔️:     target->setFormalCharge(-1);
    // RDKit✔️✔️:     target->calcExplicitValence(false);
    // RDKit✔️✔️:     stack.top()->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:     atom->setFormalCharge(1);
    // RDKit✔️✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUp3
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp3
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: The selected cleanUp caller supplies a valid atom pointer in a valid RWMol;
    // RDKit✔️✔️: both maxPathLength-1 searches use the closed direct-neighbor behavior, and
    // RDKit✔️✔️: calcExplicitValence(false) preserves source mutation-before-exception order.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp3

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    let Some(target_index) = find_alternating_bonds(
        mol,
        atom_index,
        7,
        0,
        BondType::Double,
        BondType::Double,
        0,
        1,
        None,
        &mut path,
        &mut visited,
    ) else {
        return Ok(false);
    };

    let mut oxygen_path = Vec::new();
    let mut oxygen_visited = BTreeSet::new();
    let oxygen_target = find_alternating_bonds(
        mol,
        atom_index,
        8,
        0,
        BondType::Double,
        BondType::Double,
        0,
        1,
        None,
        &mut oxygen_path,
        &mut oxygen_visited,
    );
    if oxygen_target.is_none() {
        mol.atoms[target_index as usize].formal_charge = -1;
        cleanup_explicit_valence(mol, target_index)?;
        let target_bond_index = *path.last().expect("a found target has a path bond");
        mol.bonds[target_bond_index as usize].bond_type = BondType::Single;
        mol.atoms[atom_index as usize].formal_charge = 1;
        cleanup_explicit_valence(mol, atom_index)?;
    }
    Ok(true)
}

pub(crate) fn valence5n_cleanup4(mol: &mut AdapterMol, atom_index: u32) -> bool {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:477 _Valence5NCleanUp4
    // RDKit✔️✔️: bool _Valence5NCleanUp4(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   RWMol::ADJ_ITER nid1, end1;
    // RDKit✔️✔️:   int nSi = 0;
    // RDKit✔️✔️:   int thisId = atom->getIdx();
    // RDKit✔️✔️:   Atom *nbrs[2];
    // RDKit✔️✔️:   Bond *bonds[2];
    // RDKit✔️✔️:   boost::tie(nid1, end1) = mol.getAtomNeighbors(atom);
    // RDKit✔️✔️:   while (nid1 != end1) {
    // RDKit✔️✔️:     Atom *nbr = mol.getAtomWithIdx(*nid1);
    // RDKit✔️✔️:     Bond *bond = mol.getBondBetweenAtoms(*nid1, thisId);
    // RDKit✔️✔️:     if (nbr->getAtomicNum() == 14 && nbr->getFormalCharge() == -1 &&
    // RDKit✔️✔️:         bond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       if (nSi >= 2) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       nbrs[nSi] = nbr;
    // RDKit✔️✔️:       bonds[nSi] = bond;
    // RDKit✔️✔️:       nSi++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++nid1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (nSi != 2) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   nbrs[0]->setFormalCharge(0);
    // RDKit✔️✔️:   nbrs[1]->setFormalCharge(0);
    // RDKit✔️✔️:   bonds[0]->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   bonds[1]->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:
    // RDKit✔️✔️: #if 0
    // RDKit✔️✔️:       // FIX
    // RDKit✔️✔️:       // not clear why this is here, but it almost definitely shouldn't be
    // RDKit✔️✔️:       Atom* s = NULL;
    // RDKit✔️✔️:       Atom* c = NULL;
    // RDKit✔️✔️:       Bond* sc_bond;
    // RDKit✔️✔️:       ROMol::VERTEX_ITER atBegin,atEnd;
    // RDKit✔️✔️:       boost::tie(atBegin,atEnd) = mol.getVertices();
    // RDKit✔️✔️:       while (atBegin != atEnd) {
    // RDKit✔️✔️:           ATOM_SPTR at2 = mol[*atBegin];
    // RDKit✔️✔️:           if (at2->getAtomicNum() == 16 && at2->getFormalCharge() == 1) {
    // RDKit✔️✔️:             boost::tie(nid1, end1) = mol.getAtomNeighbors(at2);
    // RDKit✔️✔️:            while (nid1 != end1) {
    // RDKit✔️✔️:              Atom* nbr = mol.getAtomWithIdx(*nid1);
    // RDKit✔️✔️:              Bond* bond = mol.getBondBetweenAtoms(*nid1, at2->getIdx());
    // RDKit✔️✔️:              if (nbr->getAtomicNum() == 6 && nbr->getFormalCharge() == 0 &&
    // RDKit✔️✔️:                  bond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:                s = &(*at2);
    // RDKit✔️✔️:                c = nbr;
    // RDKit✔️✔️:                sc_bond = bond;
    // RDKit✔️✔️:                break;
    // RDKit✔️✔️:              }
    // RDKit✔️✔️:              ++nid1
    // RDKit✔️✔️:            }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           ++atBegin;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (s == NULL) return false;
    // RDKit✔️✔️:       s->setFormalCharge(0);
    // RDKit✔️✔️:       c->setFormalCharge(-1);
    // RDKit✔️✔️:       sc_bond->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:       atom->setFormalCharge(0);
    // RDKit✔️✔️: #endif
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUp4
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp4
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no target-dependent active branch.
    // RDKit✔️✔️: The entire source #if 0 sulfur/carbon block is preprocessor-inactive;
    // RDKit✔️✔️: production behavior ends after the two silicon charge and bond mutations.
    // RDKit✔️✔️: A valid RWMol has at most one bond for an atom pair, adjacency order follows
    // RDKit✔️✔️: stored out-edge order, and nSi cannot exceed two before the early return.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp4

    let mut silicon_neighbor_indices = [0_u32; 2];
    let mut silicon_bond_indices = [0_u32; 2];
    let mut silicon_count = 0_usize;

    for &(neighbor_index, bond_index) in &mol.adjacency[atom_index as usize] {
        let neighbor = &mol.atoms[neighbor_index as usize];
        let bond = &mol.bonds[bond_index as usize];
        if neighbor.atomic_number == 14 && neighbor.formal_charge == -1 && bond.bond_type == BondType::Double {
            if silicon_count >= 2 {
                return false;
            }
            silicon_neighbor_indices[silicon_count] = neighbor_index;
            silicon_bond_indices[silicon_count] = bond_index;
            silicon_count += 1;
        }
    }

    if silicon_count != 2 {
        return false;
    }
    mol.atoms[silicon_neighbor_indices[0] as usize].formal_charge = 0;
    mol.atoms[silicon_neighbor_indices[1] as usize].formal_charge = 0;
    mol.bonds[silicon_bond_indices[0] as usize].bond_type = BondType::Single;
    mol.bonds[silicon_bond_indices[1] as usize].bond_type = BondType::Single;
    true
}

pub(crate) fn valence5n_cleanup5(
    mol: &mut AdapterMol,
    atom_index: u32,
    atomic_number: i32,
) -> Result<bool, AdapterCleanup5Error> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:544 _Valence5NCleanUp5
    // RDKit✔️✔️: bool _Valence5NCleanUp5(RWMol &mol, Atom *atom, int atomicNum) {
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       atomicNum == 8 || atomicNum == 16 || atomicNum == 9 || atomicNum == 17,
    // RDKit✔️✔️:       "this cleanup looks for O or S or Cl or F");
    // RDKit✔️✔️:   std::stack<Bond *> stackCharged, stackUncharged, *stack;
    // RDKit✔️✔️:   // try search for valence-5 N connected to O or S, determined by the
    // RDKit✔️✔️:   // <atomicNum> parameter with alternating
    // RDKit✔️✔️:   // bonds if there is a charged Oxygen and an uncharged one both
    // RDKit✔️✔️:   // connected to our N through alternating bonds, strip the charge
    // RDKit✔️✔️:   // and hydrogen from the charged one, and the use the uncharged
    // RDKit✔️✔️:   // one in our procedure
    // RDKit✔️✔️:   // see InChI for PubChem compound 10775236:
    // RDKit✔️✔️:   //   CC(C1=CC=CC=N1=C2C(OC)=O)CC2=[OH+]
    // RDKit✔️✔️:   // is converted into
    // RDKit✔️✔️:   //   COC(O)=C1[n+]2ccccc2C(C)CC1=O
    // RDKit✔️✔️:   Atom *unchargedOxygen, *chargedOxygen;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   unchargedOxygen =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, atomicNum, 0, Bond::DOUBLE, Bond::DOUBLE,
    // RDKit✔️✔️:                            0, 7, nullptr, stackUncharged, _visited);
    // RDKit✔️✔️:   chargedOxygen =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, atomicNum, 1, Bond::DOUBLE, Bond::DOUBLE,
    // RDKit✔️✔️:                            0, 7, nullptr, stackCharged, _visited);
    // RDKit✔️✔️:   if (unchargedOxygen == nullptr && chargedOxygen == nullptr) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   stack = &stackUncharged;
    // RDKit✔️✔️:   if (unchargedOxygen == nullptr) {
    // RDKit✔️✔️:     stack = &stackCharged;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (unchargedOxygen && chargedOxygen) {
    // RDKit✔️✔️:     // both exists. fix the charged oxygen now by set it to neutral
    // RDKit✔️✔️:     // with its hydrogen taken and moved later to the uncharged one
    // RDKit✔️✔️:     CHECK_INVARIANT(chargedOxygen->getFormalCharge() == 1,
    // RDKit✔️✔️:                     "expecting +1 charge");
    // RDKit✔️✔️:     chargedOxygen->setFormalCharge(0);
    // RDKit✔️✔️:     chargedOxygen->setNumExplicitHs(0);  // this hydrogen will be
    // RDKit✔️✔️:                                          // added to the uncharged
    // RDKit✔️✔️:                                          // oxygen later
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (unchargedOxygen || chargedOxygen) {
    // RDKit✔️✔️:     // set charge on N
    // RDKit✔️✔️:     atom->setFormalCharge(1);
    // RDKit✔️✔️:     // switch all bonds
    // RDKit✔️✔️:     Bond *b;
    // RDKit✔️✔️:     while (!stack->empty()) {
    // RDKit✔️✔️:       b = stack->top();
    // RDKit✔️✔️:       if (b->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:         b->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         b->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       stack->pop();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (unchargedOxygen && chargedOxygen) {
    // RDKit✔️✔️:       // both charged and uncharged oxygen are found, the uncharged
    // RDKit✔️✔️:       // remains uncharged and take the hydrogen from the charged
    // RDKit✔️✔️:       // one
    // RDKit✔️✔️:       unchargedOxygen->setNumExplicitHs(1);
    // RDKit✔️✔️:     } else if (unchargedOxygen) {
    // RDKit✔️✔️:       // if only uncharged oxygen is found, not the oxygen has -1
    // RDKit✔️✔️:       // charge
    // RDKit✔️✔️:       unchargedOxygen->setFormalCharge(-1);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // if only charged oxygen is found, it's neutral now (and
    // RDKit✔️✔️:       // keeps its hydrogen)
    // RDKit✔️✔️:       chargedOxygen->setFormalCharge(0);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (chargedOxygen) {
    // RDKit✔️✔️:       chargedOxygen->calcExplicitValence(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (unchargedOxygen) {
    // RDKit✔️✔️:       unchargedOxygen->calcExplicitValence(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUp5
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp5
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: PRECONDITION is active and rejects values other than O, S, F, and Cl;
    // RDKit✔️✔️: each root search resets the shared visited set and has a maximum depth of 7.
    // RDKit✔️✔️: CHECK_INVARIANT is source-unreachable after a successful +1-charge search;
    // RDKit✔️✔️: calcExplicitValence(false) preserves charged-before-uncharged exception order.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp5

    if !matches!(atomic_number, 8 | 16 | 9 | 17) {
        return Err(AdapterCleanup5Error::Precondition {
            kind: "Pre-condition Violation",
            expression: "atomicNum == 8 || atomicNum == 16 || atomicNum == 9 || atomicNum == 17",
            message: "this cleanup looks for O or S or Cl or F",
        });
    }

    let mut uncharged_path = Vec::new();
    let mut visited = BTreeSet::new();
    let uncharged_target = find_alternating_bonds(
        mol,
        atom_index,
        atomic_number,
        0,
        BondType::Double,
        BondType::Double,
        0,
        7,
        None,
        &mut uncharged_path,
        &mut visited,
    );
    let mut charged_path = Vec::new();
    let charged_target = find_alternating_bonds(
        mol,
        atom_index,
        atomic_number,
        1,
        BondType::Double,
        BondType::Double,
        0,
        7,
        None,
        &mut charged_path,
        &mut visited,
    );
    if uncharged_target.is_none() && charged_target.is_none() {
        return Ok(false);
    }

    if let (Some(_), Some(charged_index)) = (uncharged_target, charged_target) {
        mol.atoms[charged_index as usize].formal_charge = 0;
        mol.atoms[charged_index as usize].num_explicit_hydrogens = 0;
    }

    mol.atoms[atom_index as usize].formal_charge = 1;
    let selected_path = if uncharged_target.is_some() {
        &mut uncharged_path
    } else {
        &mut charged_path
    };
    while let Some(bond_index) = selected_path.pop() {
        let bond = &mut mol.bonds[bond_index as usize];
        bond.bond_type = if bond.bond_type == BondType::Double {
            BondType::Single
        } else {
            BondType::Double
        };
    }

    match (uncharged_target, charged_target) {
        (Some(uncharged_index), Some(_)) => {
            mol.atoms[uncharged_index as usize].num_explicit_hydrogens = 1;
        }
        (Some(uncharged_index), None) => {
            mol.atoms[uncharged_index as usize].formal_charge = -1;
        }
        (None, Some(charged_index)) => {
            mol.atoms[charged_index as usize].formal_charge = 0;
        }
        (None, None) => unreachable!("no-target case returned above"),
    }
    if let Some(charged_index) = charged_target {
        cleanup_explicit_valence(mol, charged_index)?;
    }
    if let Some(uncharged_index) = uncharged_target {
        cleanup_explicit_valence(mol, uncharged_index)?;
    }
    Ok(true)
}

pub(crate) fn valence5n_cleanup6(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:625 _Valence5NCleanUp6
    // RDKit✔️✔️: bool _Valence5NCleanUp6(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   // replace the N with Sn
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != 0 ||
    // RDKit✔️✔️:       atom->calcExplicitValence(false) != 5) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   atom->setAtomicNum(50);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // substructure matching
    // RDKit✔️✔️:   auto *query = new RWMol();
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 0
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 1
    // RDKit✔️✔️:   query->addAtom(new Atom(50), false, true);  // 2
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 3
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 4
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);   // 5
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 6
    // RDKit✔️✔️:   query->addBond(0, 1, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(1, 2, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(2, 3, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(3, 4, Bond::UNSPECIFIED);
    // RDKit✔️✔️:   query->addBond(4, 5, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(5, 0, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(2, 6, Bond::SINGLE);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<MatchVectType> fgpMatches;
    // RDKit✔️✔️:   SubstructMatch(mol, *query, fgpMatches);
    // RDKit✔️✔️:   delete query;
    // RDKit✔️✔️:   // no action if none or more than one match was found
    // RDKit✔️✔️:   if (fgpMatches.size() != 1) {
    // RDKit✔️✔️:     atom->setAtomicNum(7);
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // collect matching atoms
    // RDKit✔️✔️:   int map[7];
    // RDKit✔️✔️:   MatchVectType match = fgpMatches[0];
    // RDKit✔️✔️:   for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
    // RDKit✔️✔️:        mi++) {
    // RDKit✔️✔️:     map[mi->first] = mi->second;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // flip bonds
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[0], map[1])->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[1], map[2])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[4], map[5])->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[5], map[0])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   // change the problematic N
    // RDKit✔️✔️:   atom->setAtomicNum(7);
    // RDKit✔️✔️:   atom->setFormalCharge(1);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUp6
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp6
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: Plain Atom matching compares atomic numbers; plain UNSPECIFIED Bond matching
    // RDKit✔️✔️: accepts every bond type; default SubstructMatch uniquifies by matched atom set.
    // RDKit✔️✔️: The whole molecule is searched, so the unique match need not contain atom.
    // RDKit✔️✔️: calcExplicitValence(false) is evaluated only after atomic number and charge.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp6

    let atom_position = atom_index as usize;
    if mol.atoms[atom_position].atomic_number != 7
        || mol.atoms[atom_position].formal_charge != 0
        || cleanup_explicit_valence(mol, atom_index)? != 5
    {
        return Ok(false);
    }

    mol.atoms[atom_position].atomic_number = 50;
    let matches = valence5n_cleanup6_matches(mol);
    if matches.len() != 1 {
        mol.atoms[atom_position].atomic_number = 7;
        return Ok(false);
    }

    let mapping = matches[0];
    for (first, second, bond_type) in [
        (mapping[0], mapping[1], BondType::Double),
        (mapping[1], mapping[2], BondType::Single),
        (mapping[4], mapping[5], BondType::Double),
        (mapping[5], mapping[0], BondType::Single),
    ] {
        let bond_index =
            bond_index_between(mol, first, second).expect("the selected substructure match contains every query bond");
        mol.bonds[bond_index as usize].bond_type = bond_type;
    }
    mol.atoms[atom_position].atomic_number = 7;
    mol.atoms[atom_position].formal_charge = 1;
    Ok(true)
}

pub(crate) fn valence5n_cleanup7(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:679 _Valence5NCleanUp7
    // RDKit✔️✔️: bool _Valence5NCleanUp7(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   // is it connected to O via alternating bonds?
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 8, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 5,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:   if (target == nullptr) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // replace the N with Sn
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != 0 ||
    // RDKit✔️✔️:       atom->calcExplicitValence(false) != 5) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   atom->setAtomicNum(50);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // substructure matching
    // RDKit✔️✔️:   auto *query = new RWMol();
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 0
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 1
    // RDKit✔️✔️:   query->addAtom(new Atom(50), false, true);  // 2
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);   // 3
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 4
    // RDKit✔️✔️:   query->addAtom(new Atom(8), false, true);   // 5
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 6
    // RDKit✔️✔️:   query->addBond(0, 1, Bond::UNSPECIFIED);
    // RDKit✔️✔️:   query->addBond(1, 2, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(2, 3, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(3, 4, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(4, 5, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(5, 0, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(2, 6, Bond::SINGLE);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<MatchVectType> fgpMatches;
    // RDKit✔️✔️:   SubstructMatch(mol, *query, fgpMatches);
    // RDKit✔️✔️:   delete query;
    // RDKit✔️✔️:   // no action if none or more than one match was found
    // RDKit✔️✔️:   if (fgpMatches.size() != 1) {
    // RDKit✔️✔️:     atom->setAtomicNum(7);
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // collect matching atoms
    // RDKit✔️✔️:   int map[7];
    // RDKit✔️✔️:   MatchVectType match = fgpMatches[0];
    // RDKit✔️✔️:   for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
    // RDKit✔️✔️:        mi++) {
    // RDKit✔️✔️:     map[mi->first] = mi->second;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // flip bonds
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[1], map[2])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   Bond *b;
    // RDKit✔️✔️:   while (!stack.empty()) {
    // RDKit✔️✔️:     b = stack.top();
    // RDKit✔️✔️:     if (b->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       b->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       b->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     stack.pop();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // set charge on oxygen
    // RDKit✔️✔️:   target->setFormalCharge(-1);
    // RDKit✔️✔️:   // change the problematic N
    // RDKit✔️✔️:   atom->setAtomicNum(7);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUp7
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp7
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: findAlternatingBonds runs before every atom precondition and keeps the
    // RDKit✔️✔️: shortest discovered stack under source adjacency iteration and shared visitation.
    // RDKit✔️✔️: Plain query atoms compare atomic numbers, UNSPECIFIED accepts every modeled
    // RDKit✔️✔️: bond type, and default SubstructMatch uniquifies by matched atom set.
    // RDKit✔️✔️: The query searches the whole molecule and need not contain atom or target.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp7

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    let Some(target_index) = find_alternating_bonds(
        mol,
        atom_index,
        8,
        0,
        BondType::Double,
        BondType::Double,
        0,
        5,
        None,
        &mut path,
        &mut visited,
    ) else {
        return Ok(false);
    };

    let atom_position = atom_index as usize;
    if mol.atoms[atom_position].atomic_number != 7
        || mol.atoms[atom_position].formal_charge != 0
        || cleanup_explicit_valence(mol, atom_index)? != 5
    {
        return Ok(false);
    }

    mol.atoms[atom_position].atomic_number = 50;
    let matches = valence5n_cleanup7_matches(mol);
    if matches.len() != 1 {
        mol.atoms[atom_position].atomic_number = 7;
        return Ok(false);
    }

    let mapping = matches[0];
    let query_bond = bond_index_between(mol, mapping[1], mapping[2])
        .expect("the selected substructure match contains query bond 1-2");
    mol.bonds[query_bond as usize].bond_type = BondType::Single;
    while let Some(bond_index) = path.pop() {
        let bond = &mut mol.bonds[bond_index as usize];
        bond.bond_type = if bond.bond_type == BondType::Double {
            BondType::Single
        } else {
            BondType::Double
        };
    }
    mol.atoms[target_index as usize].formal_charge = -1;
    mol.atoms[atom_position].atomic_number = 7;
    Ok(true)
}

pub(crate) fn valence5n_cleanup8(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:750 _Valence5NCleanUp8
    // RDKit✔️✔️: bool _Valence5NCleanUp8(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   // replace the N with Sn
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != 0 ||
    // RDKit✔️✔️:       atom->calcExplicitValence(false) != 5) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   atom->setAtomicNum(50);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // substructure matching
    // RDKit✔️✔️:   auto *query = new RWMol();
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 0
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);   // 1
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 2
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);   // 3
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);   // 4
    // RDKit✔️✔️:   query->addAtom(new Atom(50), false, true);  // 5
    // RDKit✔️✔️:   query->addBond(0, 1, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(1, 2, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(2, 3, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(3, 4, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(4, 0, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(5, 0, Bond::DOUBLE);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<MatchVectType> fgpMatches;
    // RDKit✔️✔️:   SubstructMatch(mol, *query, fgpMatches);
    // RDKit✔️✔️:   delete query;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (fgpMatches.size() != 1) {
    // RDKit✔️✔️:     atom->setAtomicNum(7);
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // collect matching atoms
    // RDKit✔️✔️:   int map[6];
    // RDKit✔️✔️:   MatchVectType match = fgpMatches[0];
    // RDKit✔️✔️:   for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
    // RDKit✔️✔️:        mi++) {
    // RDKit✔️✔️:     map[mi->first] = mi->second;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // flip bonds
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[1], map[2])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[2], map[3])->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[3], map[4])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[4], map[0])->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[5], map[0])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   mol.getAtomWithIdx(map[1])->setFormalCharge(-1);
    // RDKit✔️✔️:   // change the problematic N
    // RDKit✔️✔️:   atom->setAtomicNum(7);
    // RDKit✔️✔️:   atom->setFormalCharge(1);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUp8
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp8
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: Plain query atoms compare atomic numbers and default SubstructMatch
    // RDKit✔️✔️: searches the whole molecule and uniquifies by matched atom set.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp8

    let atom_position = atom_index as usize;
    if mol.atoms[atom_position].atomic_number != 7
        || mol.atoms[atom_position].formal_charge != 0
        || cleanup_explicit_valence(mol, atom_index)? != 5
    {
        return Ok(false);
    }
    mol.atoms[atom_position].atomic_number = 50;
    let matches = valence5n_cleanup8_matches(mol);
    if matches.len() != 1 {
        mol.atoms[atom_position].atomic_number = 7;
        return Ok(false);
    }
    let mapping = matches[0];
    for (first, second, bond_type) in [
        (mapping[1], mapping[2], BondType::Single),
        (mapping[2], mapping[3], BondType::Double),
        (mapping[3], mapping[4], BondType::Single),
        (mapping[4], mapping[0], BondType::Double),
        (mapping[5], mapping[0], BondType::Single),
    ] {
        let bond_index =
            bond_index_between(mol, first, second).expect("the selected substructure match contains every query bond");
        mol.bonds[bond_index as usize].bond_type = bond_type;
    }
    mol.atoms[mapping[1] as usize].formal_charge = -1;
    mol.atoms[atom_position].atomic_number = 7;
    mol.atoms[atom_position].formal_charge = 1;
    Ok(true)
}

pub(crate) fn valence5n_cleanup9(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:804 _Valence5NCleanUp9
    // RDKit✔️✔️: bool _Valence5NCleanUp9(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   // replace the N with Sn
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != 0 ||
    // RDKit✔️✔️:       atom->calcExplicitValence(false) != 5) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   atom->setAtomicNum(50);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // substructure matching
    // RDKit✔️✔️:   auto *query = new RWMol();
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 0
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);   // 1
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);   // 2
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 3
    // RDKit✔️✔️:   query->addAtom(new Atom(6), false, true);   // 4
    // RDKit✔️✔️:   query->addAtom(new Atom(50), false, true);  // 5
    // RDKit✔️✔️:   query->addBond(0, 1, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(1, 2, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(2, 3, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(3, 4, Bond::DOUBLE);
    // RDKit✔️✔️:   query->addBond(4, 0, Bond::SINGLE);
    // RDKit✔️✔️:   query->addBond(5, 0, Bond::DOUBLE);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<MatchVectType> fgpMatches;
    // RDKit✔️✔️:   SubstructMatch(mol, *query, fgpMatches);
    // RDKit✔️✔️:   delete query;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (fgpMatches.size() != 1) {
    // RDKit✔️✔️:     atom->setAtomicNum(7);
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // collect matching atoms
    // RDKit✔️✔️:   int map[6];
    // RDKit✔️✔️:   MatchVectType match = fgpMatches[0];
    // RDKit✔️✔️:   for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
    // RDKit✔️✔️:        mi++) {
    // RDKit✔️✔️:     map[mi->first] = mi->second;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // flip bonds
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[0], map[1])->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[1], map[2])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   mol.getBondBetweenAtoms(map[5], map[0])->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   mol.getAtomWithIdx(map[2])->setFormalCharge(-1);
    // RDKit✔️✔️:   // change the problematic N
    // RDKit✔️✔️:   atom->setAtomicNum(7);
    // RDKit✔️✔️:   atom->setFormalCharge(1);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUp9
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp9
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: Plain query atoms compare atomic numbers and default SubstructMatch
    // RDKit✔️✔️: searches the whole molecule and uniquifies by matched atom set.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUp9

    let atom_position = atom_index as usize;
    if mol.atoms[atom_position].atomic_number != 7
        || mol.atoms[atom_position].formal_charge != 0
        || cleanup_explicit_valence(mol, atom_index)? != 5
    {
        return Ok(false);
    }
    mol.atoms[atom_position].atomic_number = 50;
    let matches = valence5n_cleanup9_matches(mol);
    if matches.len() != 1 {
        mol.atoms[atom_position].atomic_number = 7;
        return Ok(false);
    }
    let mapping = matches[0];
    for (first, second, bond_type) in [
        (mapping[0], mapping[1], BondType::Double),
        (mapping[1], mapping[2], BondType::Single),
        (mapping[5], mapping[0], BondType::Single),
    ] {
        let bond_index =
            bond_index_between(mol, first, second).expect("the selected substructure match contains every query bond");
        mol.bonds[bond_index as usize].bond_type = bond_type;
    }
    mol.atoms[mapping[2] as usize].formal_charge = -1;
    mol.atoms[atom_position].atomic_number = 7;
    mol.atoms[atom_position].formal_charge = 1;
    Ok(true)
}

pub(crate) fn valence5n_cleanupa(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:855 _Valence5NCleanUpA
    // RDKit✔️✔️: bool _Valence5NCleanUpA(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   // replace the N with Sn
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != 0 ||
    // RDKit✔️✔️:       atom->calcExplicitValence(false) != 5) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // first find the N=N
    // RDKit✔️✔️:   auto *query = new RWMol();
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);  // 0
    // RDKit✔️✔️:   query->addAtom(new Atom(7), false, true);  // 1
    // RDKit✔️✔️:   query->addBond(0, 1, Bond::DOUBLE);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<MatchVectType> fgpMatches;
    // RDKit✔️✔️:   SubstructMatch(mol, *query, fgpMatches);
    // RDKit✔️✔️:   delete query;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (fgpMatches.size() == 0) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::stack<Bond *> bestPath;
    // RDKit✔️✔️:   for (const auto &match : fgpMatches) {
    // RDKit✔️✔️:     // does the match contains the current atom?
    // RDKit✔️✔️:     if (match[0].second == static_cast<int>(atom->getIdx()) ||
    // RDKit✔️✔️:         match[1].second == static_cast<int>(atom->getIdx())) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // set both matched N to Sn
    // RDKit✔️✔️:     mol.getAtomWithIdx(match[0].second)->setAtomicNum(50);
    // RDKit✔️✔️:     mol.getAtomWithIdx(match[1].second)->setAtomicNum(50);
    // RDKit✔️✔️:     // now search the path from current atom to these atoms
    // RDKit✔️✔️:     std::stack<Bond *> stack;
    // RDKit✔️✔️:     std::set<int> _visited;
    // RDKit✔️✔️:     Atom *target =
    // RDKit✔️✔️:         findAlternatingBonds(mol, atom, 50, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 9,
    // RDKit✔️✔️:                              nullptr, stack, _visited);
    // RDKit✔️✔️:     if (target && (bestPath.empty() || stack.size() < bestPath.size())) {
    // RDKit✔️✔️:       bestPath = stack;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     mol.getAtomWithIdx(match[0].second)->setAtomicNum(7);
    // RDKit✔️✔️:     mol.getAtomWithIdx(match[1].second)->setAtomicNum(7);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!bestPath.empty()) {
    // RDKit✔️✔️:     while (!bestPath.empty()) {
    // RDKit✔️✔️:       Bond *bond = bestPath.top();
    // RDKit✔️✔️:       if (bond->getBondType() == Bond::SINGLE) {
    // RDKit✔️✔️:         bond->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         bond->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       bestPath.pop();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     atom->setFormalCharge(1);
    // RDKit✔️✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUpA
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUpA
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: Default SubstructMatch uniquifies N=N mappings by atom set;
    // RDKit✔️✔️: findAlternatingBonds uses source adjacency order and shared visited state per match.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUpA

    let atom_position = atom_index as usize;
    if mol.atoms[atom_position].atomic_number != 7
        || mol.atoms[atom_position].formal_charge != 0
        || cleanup_explicit_valence(mol, atom_index)? != 5
    {
        return Ok(false);
    }
    let matches = valence5n_cleanupa_matches(mol);
    if matches.is_empty() {
        return Ok(false);
    }

    let mut best_path = Vec::<u32>::new();
    for mapping in matches {
        if mapping.contains(&atom_index) {
            continue;
        }
        mol.atoms[mapping[0] as usize].atomic_number = 50;
        mol.atoms[mapping[1] as usize].atomic_number = 50;
        let mut path = Vec::new();
        let mut visited = BTreeSet::new();
        let target = find_alternating_bonds(
            mol,
            atom_index,
            50,
            0,
            BondType::Double,
            BondType::Double,
            0,
            9,
            None,
            &mut path,
            &mut visited,
        );
        if target.is_some() && (best_path.is_empty() || path.len() < best_path.len()) {
            best_path = path;
        }
        mol.atoms[mapping[0] as usize].atomic_number = 7;
        mol.atoms[mapping[1] as usize].atomic_number = 7;
    }

    if best_path.is_empty() {
        return Ok(false);
    }
    while let Some(bond_index) = best_path.pop() {
        let bond = &mut mol.bonds[bond_index as usize];
        bond.bond_type = if bond.bond_type == BondType::Single {
            BondType::Double
        } else {
            BondType::Single
        };
    }
    mol.atoms[atom_position].formal_charge = 1;
    cleanup_explicit_valence(mol, atom_index)?;
    Ok(true)
}

pub(crate) fn valence5n_cleanupb(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:916 _Valence5NCleanUpB
    // RDKit✔️✔️: bool _Valence5NCleanUpB(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 6, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 1,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:   if (target == nullptr) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   target->setFormalCharge(-1);
    // RDKit✔️✔️:   target->calcExplicitValence(false);
    // RDKit✔️✔️:   stack.top()->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   atom->setFormalCharge(1);
    // RDKit✔️✔️:   atom->calcExplicitValence(false);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5NCleanUpB
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUpB
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: findAlternatingBonds uses source adjacency order, exact zero charge,
    // RDKit✔️✔️: a one-bond maximum path, and shared visited/path state.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5NCleanUpB

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    let Some(target_index) = find_alternating_bonds(
        mol,
        atom_index,
        6,
        0,
        BondType::Double,
        BondType::Double,
        0,
        1,
        None,
        &mut path,
        &mut visited,
    ) else {
        return Ok(false);
    };

    mol.atoms[target_index as usize].formal_charge = -1;
    cleanup_explicit_valence(mol, target_index)?;
    let bond_index = *path
        .last()
        .expect("a successful one-bond alternating search records its ending bond");
    mol.bonds[bond_index as usize].bond_type = BondType::Single;
    mol.atoms[atom_index as usize].formal_charge = 1;
    cleanup_explicit_valence(mol, atom_index)?;
    Ok(true)
}

pub(crate) fn valence7s_cleanup1(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:941 _Valence7SCleanUp1
    // RDKit✔️✔️: bool _Valence7SCleanUp1(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 16 || atom->getFormalCharge() != -1 ||
    // RDKit✔️✔️:       atom->calcExplicitValence(false) != 7) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   int aid = atom->getIdx();
    // RDKit✔️✔️:   int neighborsC = 0;
    // RDKit✔️✔️:   int neighborsO = 0;
    // RDKit✔️✔️:   RWMol::ADJ_ITER nid, nid1, end1;
    // RDKit✔️✔️:   boost::tie(nid1, end1) = mol.getAtomNeighbors(atom);
    // RDKit✔️✔️:   nid = end1;
    // RDKit✔️✔️:   while (nid1 != end1) {
    // RDKit✔️✔️:     Atom *otherAtom = mol.getAtomWithIdx(*nid1);
    // RDKit✔️✔️:     if (otherAtom->getAtomicNum() == 8) {
    // RDKit✔️✔️:       if (mol.getBondBetweenAtoms(*nid1, aid)->getBondType() != Bond::DOUBLE) {
    // RDKit✔️✔️:         neighborsO = 100;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         nid = nid1;
    // RDKit✔️✔️:         neighborsO++;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (otherAtom->getAtomicNum() == 6) {
    // RDKit✔️✔️:       if (mol.getBondBetweenAtoms(*nid1, aid)->getBondType() != Bond::SINGLE) {
    // RDKit✔️✔️:         neighborsC = 100;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         neighborsC++;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       neighborsC = 100;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     nid1++;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (nid != end1 && (neighborsC == 1 || neighborsO == 3)) {
    // RDKit✔️✔️:     mol.getBondBetweenAtoms(*nid, aid)->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:     Atom *otherAtom = mol.getAtomWithIdx(*nid);
    // RDKit✔️✔️:     otherAtom->setFormalCharge(-1);
    // RDKit✔️✔️:     atom->setFormalCharge(0);
    // RDKit✔️✔️:     otherAtom->calcExplicitValence(false);
    // RDKit✔️✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence7SCleanUp1
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence7SCleanUp1
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: Neighbor iteration follows stored adjacency order and the selected
    // RDKit✔️✔️: oxygen is the last valid double-bonded oxygen encountered.
    // END RDKIT ACTIVE CONFIGURATION: _Valence7SCleanUp1

    let atom_position = atom_index as usize;
    if mol.atoms[atom_position].atomic_number != 16
        || mol.atoms[atom_position].formal_charge != -1
        || cleanup_explicit_valence(mol, atom_index)? != 7
    {
        return Ok(false);
    }

    let mut neighbors_c = 0_i32;
    let mut neighbors_o = 0_i32;
    let mut selected_oxygen = None;
    for &(neighbor_index, bond_index) in &mol.adjacency[atom_position] {
        let other_atomic_number = mol.atoms[neighbor_index as usize].atomic_number;
        if other_atomic_number == 8 {
            if mol.bonds[bond_index as usize].bond_type != BondType::Double {
                neighbors_o = 100;
                break;
            }
            selected_oxygen = Some((neighbor_index, bond_index));
            neighbors_o = neighbors_o.wrapping_add(1);
        } else if other_atomic_number == 6 {
            if mol.bonds[bond_index as usize].bond_type != BondType::Single {
                neighbors_c = 100;
                break;
            }
            neighbors_c = neighbors_c.wrapping_add(1);
        } else {
            neighbors_c = 100;
            break;
        }
    }

    let Some((oxygen_index, oxygen_bond_index)) = selected_oxygen else {
        return Ok(false);
    };
    if neighbors_c != 1 && neighbors_o != 3 {
        return Ok(false);
    }
    mol.bonds[oxygen_bond_index as usize].bond_type = BondType::Single;
    mol.atoms[oxygen_index as usize].formal_charge = -1;
    mol.atoms[atom_position].formal_charge = 0;
    cleanup_explicit_valence(mol, oxygen_index)?;
    cleanup_explicit_valence(mol, atom_index)?;
    Ok(true)
}

pub(crate) fn valence7s_cleanup2(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:989 _Valence7SCleanUp2
    // RDKit✔️✔️: bool _Valence7SCleanUp2(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 16 || atom->getFormalCharge() != -1 ||
    // RDKit✔️✔️:       atom->calcExplicitValence(false) != 7) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 7, 0, Bond::DOUBLE, Bond::TRIPLE, 0, 3,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:   if (target) {
    // RDKit✔️✔️:     while (!stack.empty()) {
    // RDKit✔️✔️:       Bond *bond = stack.top();
    // RDKit✔️✔️:       if (bond->getBondType() == Bond::SINGLE) {
    // RDKit✔️✔️:         bond->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:       } else if (bond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:         bond->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:       } else if (bond->getBondType() == Bond::TRIPLE) {
    // RDKit✔️✔️:         bond->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       stack.pop();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     atom->setFormalCharge(0);
    // RDKit✔️✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence7SCleanUp2
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence7SCleanUp2
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: findAlternatingBonds clears root search state, uses global visited
    // RDKit✔️✔️: adjacency order, retains the first equal-length path, and records
    // RDKit✔️✔️: the selected path in reverse traversal order on the source stack.
    // END RDKIT ACTIVE CONFIGURATION: _Valence7SCleanUp2

    let atom_position = atom_index as usize;
    if mol.atoms[atom_position].atomic_number != 16
        || mol.atoms[atom_position].formal_charge != -1
        || cleanup_explicit_valence(mol, atom_index)? != 7
    {
        return Ok(false);
    }

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    if find_alternating_bonds(
        mol,
        atom_index,
        7,
        0,
        BondType::Double,
        BondType::Triple,
        0,
        3,
        None,
        &mut path,
        &mut visited,
    )
    .is_none()
    {
        return Ok(false);
    }

    while let Some(bond_index) = path.pop() {
        let bond = &mut mol.bonds[bond_index as usize];
        bond.bond_type = match bond.bond_type {
            BondType::Single => BondType::Double,
            BondType::Double => BondType::Single,
            BondType::Triple => BondType::Double,
            bond_type => bond_type,
        };
    }
    mol.atoms[atom_position].formal_charge = 0;
    cleanup_explicit_valence(mol, atom_index)?;
    Ok(true)
}

pub(crate) fn valence7s_cleanup3(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:1021 _Valence7SCleanUp3
    // RDKit✔️✔️: bool _Valence7SCleanUp3(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 16 || atom->getFormalCharge() != -1 ||
    // RDKit✔️✔️:       atom->calcExplicitValence(false) != 7) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 7, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 1,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:   if (target) {
    // RDKit✔️✔️:     stack.top()->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:     target->setFormalCharge(-1);
    // RDKit✔️✔️:     atom->setFormalCharge(0);
    // RDKit✔️✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence7SCleanUp3
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence7SCleanUp3
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: The one-bond search retains the first matching neutral nitrogen
    // RDKit✔️✔️: in source adjacency order and the source never recalculates target valence.
    // END RDKIT ACTIVE CONFIGURATION: _Valence7SCleanUp3

    let atom_position = atom_index as usize;
    if mol.atoms[atom_position].atomic_number != 16
        || mol.atoms[atom_position].formal_charge != -1
        || cleanup_explicit_valence(mol, atom_index)? != 7
    {
        return Ok(false);
    }

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    let Some(target_index) = find_alternating_bonds(
        mol,
        atom_index,
        7,
        0,
        BondType::Double,
        BondType::Double,
        0,
        1,
        None,
        &mut path,
        &mut visited,
    ) else {
        return Ok(false);
    };
    let bond_index = *path
        .last()
        .expect("a successful one-bond alternating search records its ending bond");
    mol.bonds[bond_index as usize].bond_type = BondType::Single;
    mol.atoms[target_index as usize].formal_charge = -1;
    mol.atoms[atom_position].formal_charge = 0;
    cleanup_explicit_valence(mol, atom_index)?;
    Ok(true)
}

pub(crate) fn valence8s_cleanup1(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:1044 _Valence8SCleanUp1
    // RDKit✔️✔️: bool _Valence8SCleanUp1(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 16 || atom->getFormalCharge() != -1 ||
    // RDKit✔️✔️:       atom->calcExplicitValence(false) != 7) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 7, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 9,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!target) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (!stack.empty()) {
    // RDKit✔️✔️:     if (stack.top()->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       stack.top()->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       stack.top()->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     stack.pop();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   target->setFormalCharge(-1);
    // RDKit✔️✔️:   target->calcExplicitValence(false);
    // RDKit✔️✔️:   target->setNumExplicitHs(0);
    // RDKit✔️✔️:   atom->setFormalCharge(0);
    // RDKit✔️✔️:   atom->calcExplicitValence(false);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence8SCleanUp1
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence8SCleanUp1
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: findAlternatingBonds accepts paths through length nine, retains
    // RDKit✔️✔️: the first equal-length target, and stores the path in reverse order.
    // END RDKIT ACTIVE CONFIGURATION: _Valence8SCleanUp1

    let atom_position = atom_index as usize;
    if mol.atoms[atom_position].atomic_number != 16
        || mol.atoms[atom_position].formal_charge != -1
        || cleanup_explicit_valence(mol, atom_index)? != 7
    {
        return Ok(false);
    }

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    let Some(target_index) = find_alternating_bonds(
        mol,
        atom_index,
        7,
        0,
        BondType::Double,
        BondType::Double,
        0,
        9,
        None,
        &mut path,
        &mut visited,
    ) else {
        return Ok(false);
    };

    while let Some(bond_index) = path.pop() {
        let bond = &mut mol.bonds[bond_index as usize];
        bond.bond_type = if bond.bond_type == BondType::Double {
            BondType::Single
        } else {
            BondType::Double
        };
    }
    mol.atoms[target_index as usize].formal_charge = -1;
    cleanup_explicit_valence(mol, target_index)?;
    mol.atoms[target_index as usize].num_explicit_hydrogens = 0;
    mol.atoms[atom_position].formal_charge = 0;
    cleanup_explicit_valence(mol, atom_index)?;
    Ok(true)
}

pub(crate) fn valence8cl_cleanup1(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:1079 _Valence8ClCleanUp1
    // RDKit✔️✔️: bool _Valence8ClCleanUp1(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   if (atom->calcExplicitValence(false) != 8 || atom->getFormalCharge() != -1) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   int aid = atom->getIdx();
    // RDKit✔️✔️:   bool neighborsAllO = true;
    // RDKit✔️✔️:   RWMol::ADJ_ITER nid1, end1;
    // RDKit✔️✔️:   boost::tie(nid1, end1) = mol.getAtomNeighbors(atom);
    // RDKit✔️✔️:   while (nid1 != end1) {
    // RDKit✔️✔️:     if (mol.getAtomWithIdx(*nid1)->getAtomicNum() != 8) {
    // RDKit✔️✔️:       neighborsAllO = false;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     nid1++;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (neighborsAllO) {
    // RDKit✔️✔️:     atom->setFormalCharge(3);
    // RDKit✔️✔️:     boost::tie(nid1, end1) = mol.getAtomNeighbors(atom);
    // RDKit✔️✔️:     while (nid1 != end1) {
    // RDKit✔️✔️:       Bond *b = mol.getBondBetweenAtoms(aid, *nid1);
    // RDKit✔️✔️:       if (b->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:         b->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:         Atom *otherAtom = mol.getAtomWithIdx(*nid1);
    // RDKit✔️✔️:         otherAtom->setFormalCharge(-1);
    // RDKit✔️✔️:         otherAtom->calcExplicitValence(false);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       nid1++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence8ClCleanUp1
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence8ClCleanUp1
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: The source function has no atomic-number guard, evaluates valence
    // RDKit✔️✔️: before charge, and treats an empty neighbor range as all oxygen.
    // END RDKIT ACTIVE CONFIGURATION: _Valence8ClCleanUp1

    if cleanup_explicit_valence(mol, atom_index)? != 8 || mol.atoms[atom_index as usize].formal_charge != -1 {
        return Ok(false);
    }

    if mol.adjacency[atom_index as usize]
        .iter()
        .any(|&(neighbor_index, _)| mol.atoms[neighbor_index as usize].atomic_number != 8)
    {
        return Ok(false);
    }

    mol.atoms[atom_index as usize].formal_charge = 3;
    let neighbor_count = mol.adjacency[atom_index as usize].len();
    for neighbor_position in 0..neighbor_count {
        let (neighbor_index, bond_index) = mol.adjacency[atom_index as usize][neighbor_position];
        if mol.bonds[bond_index as usize].bond_type == BondType::Double {
            mol.bonds[bond_index as usize].bond_type = BondType::Single;
            mol.atoms[neighbor_index as usize].formal_charge = -1;
            cleanup_explicit_valence(mol, neighbor_index)?;
        }
    }
    cleanup_explicit_valence(mol, atom_index)?;
    Ok(true)
}

pub(crate) fn valence5cl_cleanup1(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:1114 _Valence5ClCleanUp1
    // RDKit✔️✔️: bool _Valence5ClCleanUp1(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   if (atom->calcExplicitValence(false) != 6 || atom->getFormalCharge() != 1) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 8, -1, Bond::SINGLE, Bond::SINGLE, 0, 1,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:   if (!target) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   stack.top()->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:   atom->setFormalCharge(0);
    // RDKit✔️✔️:   target->setFormalCharge(0);
    // RDKit✔️✔️:   atom->calcExplicitValence(false);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence5ClCleanUp1
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence5ClCleanUp1
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: The function has no atomic-number guard, keeps the first equal-length
    // RDKit✔️✔️: adjacency target, and deliberately does not recalculate target valence.
    // END RDKIT ACTIVE CONFIGURATION: _Valence5ClCleanUp1

    if cleanup_explicit_valence(mol, atom_index)? != 6 || mol.atoms[atom_index as usize].formal_charge != 1 {
        return Ok(false);
    }

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    let Some(target_index) = find_alternating_bonds(
        mol,
        atom_index,
        8,
        -1,
        BondType::Single,
        BondType::Single,
        0,
        1,
        None,
        &mut path,
        &mut visited,
    ) else {
        return Ok(false);
    };

    let selected_bond_index = *path.last().expect("a found target has a path bond");
    mol.bonds[selected_bond_index as usize].bond_type = BondType::Double;
    mol.atoms[atom_index as usize].formal_charge = 0;
    mol.atoms[target_index as usize].formal_charge = 0;
    cleanup_explicit_valence(mol, atom_index)?;
    Ok(true)
}

pub(crate) fn valence3cl_cleanup1(mol: &mut AdapterMol, atom_index: u32) -> Result<bool, AdapterValenceError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:1134 _Valence3ClCleanUp1
    // RDKit✔️✔️: bool _Valence3ClCleanUp1(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   if (atom->calcExplicitValence(false) != 3 || atom->getFormalCharge() != 0) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::stack<Bond *> stack;
    // RDKit✔️✔️:   std::set<int> _visited;
    // RDKit✔️✔️:   Atom *target =
    // RDKit✔️✔️:       findAlternatingBonds(mol, atom, 16, 0, Bond::TRIPLE, Bond::TRIPLE, 0, 1,
    // RDKit✔️✔️:                            nullptr, stack, _visited);
    // RDKit✔️✔️:   if (!target) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   stack.top()->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:   atom->calcExplicitValence(false);
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: _Valence3ClCleanUp1
    // BEGIN RDKIT ACTIVE CONFIGURATION: _Valence3ClCleanUp1
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: The function has no atomic-number guard and deliberately does not
    // RDKit✔️✔️: recalculate target valence after changing the triple bond to single.
    // END RDKIT ACTIVE CONFIGURATION: _Valence3ClCleanUp1

    if cleanup_explicit_valence(mol, atom_index)? != 3 || mol.atoms[atom_index as usize].formal_charge != 0 {
        return Ok(false);
    }

    let mut path = Vec::new();
    let mut visited = BTreeSet::new();
    if find_alternating_bonds(
        mol,
        atom_index,
        16,
        0,
        BondType::Triple,
        BondType::Triple,
        0,
        1,
        None,
        &mut path,
        &mut visited,
    )
    .is_none()
    {
        return Ok(false);
    }

    let selected_bond_index = *path.last().expect("a found target has a path bond");
    mol.bonds[selected_bond_index as usize].bond_type = BondType::Single;
    cleanup_explicit_valence(mol, atom_index)?;
    Ok(true)
}

pub(crate) fn clean_up(mol: &mut AdapterMol) -> Result<(), AdapterCleanup5Error> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:1151 cleanUp
    // RDKit✔️✔️: void cleanUp(RWMol &mol) {
    // RDKit✔️✔️:   ROMol::AtomIterator ai;
    // RDKit✔️✔️:   bool aromHolder;
    // RDKit✔️✔️:   for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ++ai) {
    // RDKit✔️✔️:     switch ((*ai)->getAtomicNum()) {
    // RDKit✔️✔️:       case 7:
    // RDKit✔️✔️:         if ((*ai)->calcExplicitValence(false) == 4) {
    // RDKit✔️✔️:           if (_Valence4NCleanUp1(mol, *ai)) {
    // RDKit✔️✔️:             continue;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           if ((*ai)->getFormalCharge() == -1) {
    // RDKit✔️✔️:             if (_Valence4NCleanUp2(mol, *ai)) {
    // RDKit✔️✔️:               continue;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         if ((*ai)->getFormalCharge()) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         aromHolder = (*ai)->getIsAromatic();
    // RDKit✔️✔️:         (*ai)->setIsAromatic(0);
    // RDKit✔️✔️:
    // RDKit✔️✔️:         if ((*ai)->calcExplicitValence(false) == 5) {
    // RDKit✔️✔️:           // rings CN1=CCN=CC=1, CN1=NCOCC=1, [N]=C1N=CN=N1, [N]=C1C=CN=N1
    // RDKit✔️✔️:           (_Valence5NCleanUp6(mol, *ai)) || (_Valence5NCleanUp7(mol, *ai)) ||
    // RDKit✔️✔️:               (_Valence5NCleanUp8(mol, *ai)) ||
    // RDKit✔️✔️:               (_Valence5NCleanUp9(mol, *ai)) ||
    // RDKit✔️✔️:               (_Valence5NCleanUpA(mol, *ai)) ||
    // RDKit✔️✔️:               // try search for valence-5 N connected to a N+
    // RDKit✔️✔️:               (_Valence5NCleanUp1(mol, *ai)) ||
    // RDKit✔️✔️:               // connected to N- through a tiple then single bond
    // RDKit✔️✔️:               (_Valence5NCleanUp2(mol, *ai)) ||
    // RDKit✔️✔️:               // directly to a N
    // RDKit✔️✔️:               (_Valence5NCleanUp3(mol, *ai)) ||
    // RDKit✔️✔️:               // to two Si- via double bonds
    // RDKit✔️✔️:               (_Valence5NCleanUp4(mol, *ai)) ||
    // RDKit✔️✔️:               // alternating bonds to O
    // RDKit✔️✔️:               (_Valence5NCleanUp5(mol, *ai, 8)) ||
    // RDKit✔️✔️:               // alternating bonds to S
    // RDKit✔️✔️:               (_Valence5NCleanUp5(mol, *ai, 16)) ||
    // RDKit✔️✔️:               // alternating bonds to S
    // RDKit✔️✔️:               (_Valence5NCleanUp5(mol, *ai, 9)) ||
    // RDKit✔️✔️:               // alternating bonds to S
    // RDKit✔️✔️:               (_Valence5NCleanUp5(mol, *ai, 17)) ||
    // RDKit✔️✔️:               // last resort
    // RDKit✔️✔️:               (_Valence5NCleanUpB(mol, *ai));
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (aromHolder) {
    // RDKit✔️✔️:           (*ai)->setIsAromatic(1);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case 17:
    // RDKit✔️✔️:         if ((*ai)->calcExplicitValence(false) == 8 &&
    // RDKit✔️✔️:             _Valence8ClCleanUp1(mol, *ai)) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if ((*ai)->calcExplicitValence(false) == 5 &&
    // RDKit✔️✔️:             _Valence5ClCleanUp1(mol, *ai)) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if ((*ai)->calcExplicitValence(false) == 3 &&
    // RDKit✔️✔️:             _Valence3ClCleanUp1(mol, *ai)) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case 16:
    // RDKit✔️✔️:         if ((*ai)->calcExplicitValence(false) == 7) {
    // RDKit✔️✔️:           if (_Valence7SCleanUp1(mol, *ai)) {
    // RDKit✔️✔️:             continue;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           if (_Valence7SCleanUp2(mol, *ai)) {
    // RDKit✔️✔️:             continue;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           if (_Valence7SCleanUp3(mol, *ai)) {
    // RDKit✔️✔️:             continue;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           _Valence8SCleanUp1(mol, *ai);
    // RDKit✔️✔️:         } else if ((*ai)->calcExplicitValence(false) == 8) {
    // RDKit✔️✔️:           _Valence8SCleanUp1(mol, *ai);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case 35:
    // RDKit✔️✔️:         if ((*ai)->calcExplicitValence(false) == 3 &&
    // RDKit✔️✔️:             (*ai)->getFormalCharge() == 0) {
    // RDKit✔️✔️:           // connected to Se. Example: PubChem 10787526
    // RDKit✔️✔️:           if ((*ai)->getDegree() == 1) {
    // RDKit✔️✔️:             RWMol::ADJ_ITER nid, end;
    // RDKit✔️✔️:             boost::tie(nid, end) = mol.getAtomNeighbors(*ai);
    // RDKit✔️✔️:             if (mol.getAtomWithIdx(*nid)->getAtomicNum() == 34) {
    // RDKit✔️✔️:               mol.getBondBetweenAtoms((*ai)->getIdx(), *nid)
    // RDKit✔️✔️:                   ->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     }  // end the switch block
    // RDKit✔️✔️:   }  // end the for loop that iterates over atoms
    // RDKit✔️✔️: }  // end cleanUp
    // END RDKIT C++ FUNCTION: cleanUp
    // BEGIN RDKIT ACTIVE CONFIGURATION: cleanUp
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️✔️: Atom iteration and adjacency follow stored RDKit order; all cleanup callees
    // RDKit✔️✔️: are the source-exact adapter helpers closed above.
    // END RDKIT ACTIVE CONFIGURATION: cleanUp

    for atom_index in 0..mol.atoms.len() as u32 {
        match mol.atoms[atom_index as usize].atomic_number {
            7 => {
                if cleanup_explicit_valence(mol, atom_index)? == 4 {
                    if valence4n_cleanup1(mol, atom_index)? {
                        continue;
                    }
                    if mol.atoms[atom_index as usize].formal_charge == -1 && valence4n_cleanup2(mol, atom_index) {
                        continue;
                    }
                    continue;
                }

                if mol.atoms[atom_index as usize].formal_charge != 0 {
                    continue;
                }
                let aromatic = mol.atoms[atom_index as usize].is_aromatic;
                mol.atoms[atom_index as usize].is_aromatic = false;

                if cleanup_explicit_valence(mol, atom_index)? == 5 {
                    let mut cleaned = valence5n_cleanup6(mol, atom_index)?;
                    if !cleaned {
                        cleaned = valence5n_cleanup7(mol, atom_index)?;
                    }
                    if !cleaned {
                        cleaned = valence5n_cleanup8(mol, atom_index)?;
                    }
                    if !cleaned {
                        cleaned = valence5n_cleanup9(mol, atom_index)?;
                    }
                    if !cleaned {
                        cleaned = valence5n_cleanupa(mol, atom_index)?;
                    }
                    if !cleaned {
                        cleaned = valence5n_cleanup1(mol, atom_index)?;
                    }
                    if !cleaned {
                        cleaned = valence5n_cleanup2(mol, atom_index)?;
                    }
                    if !cleaned {
                        cleaned = valence5n_cleanup3(mol, atom_index)?;
                    }
                    if !cleaned {
                        cleaned = valence5n_cleanup4(mol, atom_index);
                    }
                    for atomic_number in [8, 16, 9, 17] {
                        if !cleaned {
                            cleaned = valence5n_cleanup5(mol, atom_index, atomic_number)?;
                        }
                    }
                    if !cleaned {
                        let _ = valence5n_cleanupb(mol, atom_index)?;
                    }
                }
                if aromatic {
                    mol.atoms[atom_index as usize].is_aromatic = true;
                }
            }
            17 => {
                if cleanup_explicit_valence(mol, atom_index)? == 8 && valence8cl_cleanup1(mol, atom_index)? {
                    continue;
                }
                if cleanup_explicit_valence(mol, atom_index)? == 5 && valence5cl_cleanup1(mol, atom_index)? {
                    continue;
                }
                if cleanup_explicit_valence(mol, atom_index)? == 3 && valence3cl_cleanup1(mol, atom_index)? {
                    continue;
                }
            }
            16 => {
                if cleanup_explicit_valence(mol, atom_index)? == 7 {
                    if valence7s_cleanup1(mol, atom_index)? {
                        continue;
                    }
                    if valence7s_cleanup2(mol, atom_index)? {
                        continue;
                    }
                    if valence7s_cleanup3(mol, atom_index)? {
                        continue;
                    }
                    let _ = valence8s_cleanup1(mol, atom_index)?;
                } else if cleanup_explicit_valence(mol, atom_index)? == 8 {
                    let _ = valence8s_cleanup1(mol, atom_index)?;
                }
            }
            35 => {
                if cleanup_explicit_valence(mol, atom_index)? == 3
                    && mol.atoms[atom_index as usize].formal_charge == 0
                    && mol.adjacency[atom_index as usize].len() == 1
                {
                    let (neighbor_index, bond_index) = mol.adjacency[atom_index as usize][0];
                    if mol.atoms[neighbor_index as usize].atomic_number == 34 {
                        mol.bonds[bond_index as usize].bond_type = BondType::Single;
                    }
                }
            }
            _ => {}
        }
    }
    Ok(())
}

fn adapter_element_name(atom: &inchi_Atom) -> Result<Vec<u8>, InchiToMolError> {
    let end = atom
        .elname
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(InchiToMolError::InvalidSourceOutput(
            "inchi atom element name is not NUL-terminated",
        ))?;
    Ok(atom.elname[..end].iter().map(|byte| *byte as u8).collect())
}

fn adapter_add_atom(molecule: &mut AdapterMol, atom: AdapterAtom) -> u32 {
    let index = molecule.atoms.len() as u32;
    molecule.atoms.push(atom);
    molecule.adjacency.push(Vec::new());
    index
}

fn adapter_add_bond(molecule: &mut AdapterMol, bond: AdapterBond) -> u32 {
    let index = molecule.bonds.len() as u32;
    let begin = bond.begin_atom_index as usize;
    let end = bond.end_atom_index as usize;
    molecule.adjacency[begin].push((bond.end_atom_index, index));
    molecule.adjacency[end].push((bond.begin_atom_index, index));
    molecule.bonds.push(bond);
    index
}

fn adapter_count_swaps_to_interconvert(reference: &[u32], mut probe: Vec<u32>) -> Result<u32, InchiToMolError> {
    if reference.len() != probe.len() {
        return Err(InchiToMolError::InvalidSourceOutput(
            "tetrahedral bond ordering size mismatch",
        ));
    }
    let mut swaps = 0_u32;
    for (position, expected) in reference.iter().copied().enumerate() {
        if probe[position] != expected {
            let Some(found) = probe[position..]
                .iter()
                .position(|value| *value == expected)
                .map(|offset| position + offset)
            else {
                return Err(InchiToMolError::InvalidSourceOutput(
                    "tetrahedral bond ordering does not contain source bond",
                ));
            };
            probe.swap(position, found);
            swaps += 1;
        }
    }
    Ok(swaps)
}

fn adapter_perturbation_order(molecule: &AdapterMol, atom_index: u32, probe: &[u32]) -> Result<u32, InchiToMolError> {
    let reference: Vec<u32> = molecule
        .adjacency
        .get(atom_index as usize)
        .ok_or(InchiToMolError::InvalidSourceOutput(
            "tetrahedral central atom index is outside molecule",
        ))?
        .iter()
        .map(|entry| entry.1)
        .collect();
    adapter_count_swaps_to_interconvert(&reference, probe.to_vec())
}

fn adapter_invert_chirality(atom: &mut AdapterAtom) {
    atom.chiral_tag = match atom.chiral_tag {
        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
        ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
        value => value,
    };
}

pub(crate) fn inchi_to_mol(
    engine: &mut impl InchiStructureEngine,
    toolkit: &mut impl InchiToMolToolkit,
    inchi: &[u8],
    return_values: &mut ExtraInchiReturnValues,
    sanitize: bool,
    remove_hs: bool,
) -> Result<InchiToMolResult, InchiToMolError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:1254 InchiToMol
    // RDKit✔️❌: RWMol *InchiToMol(const std::string &inchi, ExtraInchiReturnValues &rv,
    // RDKit✔️❌:                   bool sanitize, bool removeHs) {
    // RDKit✔️✔️:   // input
    // RDKit✔️✔️:   std::vector<char> _inchi;
    // RDKit✔️✔️:   _inchi.reserve(inchi.size() + 1);
    // RDKit✔️✔️:   std::copy(inchi.begin(), inchi.end(), std::back_inserter(_inchi));
    // RDKit✔️✔️:   _inchi.push_back('\0');

    // RDKit✔️❌:   char options[1] = "";
    // RDKit✔️❌:   inchi_InputINCHI inchiInput;
    // RDKit✔️❌:   inchiInput.szInChI = _inchi.data();
    // RDKit✔️❌:   inchiInput.szOptions = options;

    // RDKit✔️❌:   // creating RWMol for return
    // RDKit✔️❌:   RWMol *m = nullptr;
    // RDKit✔️❌:   {
    // RDKit✔️❌:     // output structure
    // RDKit✔️❌:     inchi_OutputStruct inchiOutput;
    // RDKit✔️❌:     // DLL call
    // RDKit✔️❌:     int retcode = GetStructFromINCHI(&inchiInput, &inchiOutput);

    // RDKit✔️❌:     // prepare output
    // RDKit✔️❌:     rv.returnCode = retcode;
    // RDKit✔️❌:     if (inchiOutput.szMessage) {
    // RDKit✔️❌:       rv.messagePtr = std::string(inchiOutput.szMessage);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (inchiOutput.szLog) {
    // RDKit✔️❌:       rv.logPtr = std::string(inchiOutput.szLog);
    // RDKit✔️❌:     }

    // RDKit✔️❌:     // for isotopes of H
    // RDKit✔️❌:     typedef std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
    // RDKit✔️❌:         ISOTOPES_t;
    // RDKit✔️❌:     ISOTOPES_t isotopes;
    // RDKit✔️❌:     if (retcode == inchi_Ret_OKAY || retcode == inchi_Ret_WARNING) {
    // RDKit✔️❌:       m = new RWMol;
    // RDKit✔️❌:       std::vector<unsigned int> indexToAtomIndexMapping;
    // RDKit✔️❌:       PeriodicTable *periodicTable = PeriodicTable::getTable();
    // RDKit✔️❌:       unsigned int nAtoms = inchiOutput.num_atoms;
    // RDKit✔️❌:       for (unsigned int i = 0; i < nAtoms; i++) {
    // RDKit✔️❌:         inchi_Atom *inchiAtom = &(inchiOutput.atom[i]);
    // RDKit✔️❌:         // use element name to set atomic number
    // RDKit✔️❌:         int atomicNumber = periodicTable->getAtomicNumber(inchiAtom->elname);
    // RDKit✔️❌:         Atom *atom = new Atom(atomicNumber);
    // RDKit✔️❌:         double averageWeight = atom->getMass();
    // RDKit✔️❌:         int refWeight = static_cast<int>(averageWeight + 0.5);
    // RDKit✔️❌:         int isotope = 0;
    // RDKit✔️❌:         if (inchiAtom->isotopic_mass) {
    // RDKit✔️❌:           isotope = inchiAtom->isotopic_mass - ISOTOPIC_SHIFT_FLAG;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (isotope) {
    // RDKit✔️❌:           atom->setIsotope(isotope + refWeight);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         // set charge
    // RDKit✔️❌:         atom->setFormalCharge(inchiAtom->charge);
    // RDKit✔️❌:         // set radical
    // RDKit✔️❌:         if (inchiAtom->radical) {
    // RDKit✔️❌:           if (inchiAtom->radical != 3 && inchiAtom->radical != 2) {
    // RDKit✔️❌:             BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:                 << "expect radical to be either 2 or 3 while getting "
    // RDKit✔️❌:                 << inchiAtom->radical << ". Ignore radical." << std::endl;
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             atom->setNumRadicalElectrons(inchiAtom->radical - 1);
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:         // number of hydrogens
    // RDKit✔️❌:         atom->setNumExplicitHs(inchiAtom->num_iso_H[0]);
    // RDKit✔️❌:         if (inchiAtom->num_iso_H[1]) {
    // RDKit✔️❌:           isotopes.push_back(std::make_tuple(1, i, inchiAtom->num_iso_H[1]));
    // RDKit✔️❌:         } else if (inchiAtom->num_iso_H[2]) {
    // RDKit✔️❌:           isotopes.push_back(std::make_tuple(2, i, inchiAtom->num_iso_H[2]));
    // RDKit✔️❌:         } else if (inchiAtom->num_iso_H[3]) {
    // RDKit✔️❌:           isotopes.push_back(std::make_tuple(3, i, inchiAtom->num_iso_H[3]));
    // RDKit✔️❌:         }
    // RDKit✔️❌:         // at this point the molecule has all Hs it should have. Set the
    // RDKit✔️❌:         // noImplicit flag so
    // RDKit✔️❌:         // we don't end up with extras later (this was github #562):
    // RDKit✔️❌:         atom->setNoImplicit(true);
    // RDKit✔️❌:         // add atom to molecule
    // RDKit✔️❌:         unsigned int aid = m->addAtom(atom, false, true);
    // RDKit✔️❌:         indexToAtomIndexMapping.push_back(aid);
    // RDKit❌❌: #ifdef DEBUG
    // RDKit❌❌:         BOOST_LOG(rdWarningLog)
    // RDKit❌❌:             << "adding " << aid << ":" << atom->getAtomicNum() << ":"
    // RDKit❌❌:             << (int)inchiAtom->num_iso_H[0]
    // RDKit❌❌:             << " charge: " << (int)inchiAtom->charge << std::endl;
    // RDKit❌❌: #endif
    // RDKit✔️❌:       }

    // RDKit✔️❌:       // adding bonds
    // RDKit✔️❌:       std::set<std::pair<unsigned int, unsigned int>> bondRegister;
    // RDKit✔️❌:       for (unsigned int i = 0; i < nAtoms; i++) {
    // RDKit✔️❌:         inchi_Atom *inchiAtom = &(inchiOutput.atom[i]);
    // RDKit✔️❌:         unsigned int nBonds = inchiAtom->num_bonds;
    // RDKit✔️❌:         for (unsigned int b = 0; b < nBonds; b++) {
    // RDKit✔️❌:           unsigned int nbr = inchiAtom->neighbor[b];
    // RDKit✔️❌:           // check register to avoid duplication
    // RDKit✔️❌:           if (bondRegister.find(std::make_pair(i, nbr)) != bondRegister.end() ||
    // RDKit✔️❌:               bondRegister.find(std::make_pair(nbr, i)) != bondRegister.end()) {
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           bondRegister.insert(std::make_pair(i, nbr));
    // RDKit✔️❌:           Bond *bond = nullptr;
    // RDKit✔️❌:           // bond type
    // RDKit✔️❌:           if ((unsigned int)inchiAtom->bond_type[b] <= INCHI_BOND_TYPE_TRIPLE) {
    // RDKit✔️❌:             bond = new Bond((Bond::BondType)inchiAtom->bond_type[b]);
    // RDKit✔️❌:           } else if ((unsigned int)inchiAtom->bond_type[b] ==
    // RDKit✔️❌:                      INCHI_BOND_TYPE_ALTERN) {
    // RDKit✔️❌:             BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:                 << "receive ALTERN bond type which should be avoided. "
    // RDKit✔️❌:                 << "This is treated as aromatic." << std::endl;
    // RDKit✔️❌:             bond = new Bond(Bond::AROMATIC);
    // RDKit✔️❌:             bond->setIsAromatic(true);
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             BOOST_LOG(rdErrorLog) << "illegal bond type ("
    // RDKit✔️❌:                                   << (unsigned int)inchiAtom->bond_type[b]
    // RDKit✔️❌:                                   << ") in InChI" << std::endl;
    // RDKit✔️❌:             FreeStructFromINCHI(&inchiOutput);
    // RDKit✔️❌:             delete m;
    // RDKit✔️❌:             return nullptr;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           // bond ends
    // RDKit✔️❌:           bond->setBeginAtomIdx(indexToAtomIndexMapping[i]);
    // RDKit✔️❌:           bond->setEndAtomIdx(indexToAtomIndexMapping[nbr]);
    // RDKit✔️❌:           // bond stereo
    // RDKit✔️❌:           switch (inchiAtom->bond_stereo[b]) {
    // RDKit✔️❌:             case INCHI_BOND_STEREO_NONE:
    // RDKit✔️❌:               break;
    // RDKit✔️❌:             case INCHI_BOND_STEREO_SINGLE_1UP:
    // RDKit✔️❌:             case INCHI_BOND_STEREO_SINGLE_2DOWN:
    // RDKit✔️❌:               bond->setBondDir(Bond::BEGINWEDGE);
    // RDKit✔️❌:               break;
    // RDKit✔️❌:             case INCHI_BOND_STEREO_SINGLE_1DOWN:
    // RDKit✔️❌:             case INCHI_BOND_STEREO_SINGLE_2UP:
    // RDKit✔️❌:               bond->setBondDir(Bond::BEGINDASH);
    // RDKit✔️❌:               break;
    // RDKit✔️❌:             case INCHI_BOND_STEREO_SINGLE_1EITHER:
    // RDKit✔️❌:               bond->setBondDir(Bond::UNKNOWN);
    // RDKit✔️❌:               break;
    // RDKit✔️❌:             case INCHI_BOND_STEREO_DOUBLE_EITHER:
    // RDKit✔️❌:               bond->setBondDir(Bond::EITHERDOUBLE);
    // RDKit✔️❌:               break;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           // add bond
    // RDKit✔️❌:           m->addBond(bond, true);
    // RDKit❌❌: #ifdef DEBUG
    // RDKit❌❌:           BOOST_LOG(rdWarningLog)
    // RDKit❌❌:               << "adding " << (int)bond->getBeginAtomIdx() << "("
    // RDKit❌❌:               << m->getAtomWithIdx(bond->getBeginAtomIdx())->getAtomicNum()
    // RDKit❌❌:               << ")"
    // RDKit❌❌:               << "-" << (int)bond->getEndAtomIdx() << "("
    // RDKit❌❌:               << m->getAtomWithIdx(bond->getEndAtomIdx())->getAtomicNum() << ")"
    // RDKit❌❌:               << "[" << (int)bond->getBondType() << "]" << std::endl;
    // RDKit❌❌: #endif
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }

    // RDKit✔️❌:       // adding isotopes at the end
    // RDKit✔️❌:       for (auto &ii : isotopes) {
    // RDKit✔️❌:         auto [isotope, aid, repeat] = ii;
    // RDKit✔️❌:         aid = indexToAtomIndexMapping[aid];
    // RDKit✔️❌:         for (unsigned int i = 0; i < repeat; i++) {
    // RDKit✔️❌:           // create atom
    // RDKit✔️❌:           Atom *atom = new Atom;
    // RDKit✔️❌:           atom->setAtomicNum(1);
    // RDKit✔️❌:           // set mass
    // RDKit✔️❌:           atom->setIsotope(isotope);
    // RDKit✔️❌:           int j = m->addAtom(atom, false, true);
    // RDKit✔️❌:           // add bond
    // RDKit✔️❌:           Bond *bond = new Bond(Bond::SINGLE);
    // RDKit✔️❌:           bond->setEndAtomIdx(aid);
    // RDKit✔️❌:           bond->setBeginAtomIdx(j);
    // RDKit✔️❌:           m->addBond(bond, true);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }

    // RDKit✔️❌:       // basic topological structure is ready. calculate valence
    // RDKit✔️❌:       m->updatePropertyCache(false);

    // RDKit✔️❌:       // 0Dstereo
    // RDKit✔️❌:       INT_PAIR_VECT eBondPairs;
    // RDKit✔️❌:       INT_PAIR_VECT zBondPairs;
    // RDKit✔️❌:       unsigned int numStereo0D = inchiOutput.num_stereo0D;
    // RDKit✔️❌:       if (numStereo0D) {
    // RDKit✔️❌:         // calculate CIPCode as they might be used
    // RDKit✔️❌:         UINT_VECT ranks;
    // RDKit✔️❌:         Chirality::assignAtomCIPRanks(*m, ranks);
    // RDKit✔️❌:         for (unsigned int i = 0; i < numStereo0D; i++) {
    // RDKit✔️❌:           inchi_Stereo0D *stereo0DPtr = inchiOutput.stereo0D + i;
    // RDKit✔️❌:           if (stereo0DPtr->parity == INCHI_PARITY_NONE ||
    // RDKit✔️❌:               stereo0DPtr->parity == INCHI_PARITY_UNDEFINED) {
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           switch (stereo0DPtr->type) {
    // RDKit✔️❌:             case INCHI_StereoType_None:
    // RDKit✔️❌:               break;
    // RDKit✔️❌:             case INCHI_StereoType_DoubleBond: {
    // RDKit✔️❌:               // find the bond
    // RDKit✔️❌:               unsigned left = indexToAtomIndexMapping[stereo0DPtr->neighbor[1]];
    // RDKit✔️❌:               unsigned right =
    // RDKit✔️❌:                   indexToAtomIndexMapping[stereo0DPtr->neighbor[2]];
    // RDKit✔️❌:               int originalLeftNbr =
    // RDKit✔️❌:                   indexToAtomIndexMapping[stereo0DPtr->neighbor[0]];
    // RDKit✔️❌:               int originalRightNbr =
    // RDKit✔️❌:                   indexToAtomIndexMapping[stereo0DPtr->neighbor[3]];

    // RDKit✔️❌:               Bond *bond = m->getBondBetweenAtoms(left, right);
    // RDKit✔️❌:               if (!bond) {
    // RDKit✔️❌:                 // Likely to be allene stereochemistry, which we don't handle.
    // RDKit✔️❌:                 BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:                     << "Extended double-bond stereochemistry (e.g. C=C=C=C) "
    // RDKit✔️❌:                        "ignored"
    // RDKit✔️❌:                     << std::endl;
    // RDKit✔️❌:                 continue;
    // RDKit✔️❌:               }
    // RDKit✔️❌:               // also find neighboring atoms. Note we cannot use what InChI
    // RDKit✔️❌:               // returned in stereo0DPtr->neighbor as there can be hydrogen in
    // RDKit✔️❌:               // it, which is later removed and is therefore not reliable. Plus,
    // RDKit✔️❌:               // InChI seems to use lower CIPRank-neighbors rather than
    // RDKit✔️❌:               // higher-CIPRank ones (hence the use of hydrogen neighbor).
    // RDKit✔️❌:               // However, if the neighbors we selected differ from what are in
    // RDKit✔️❌:               // stereo0DPtr->neighbor, we might also need to switch E and Z

    // RDKit✔️❌:               auto findNbrAtoms = [&m, &ranks](unsigned ref) {
    // RDKit✔️❌:                 int nbr = -1;
    // RDKit✔️❌:                 int extraNbr = -1;
    // RDKit✔️❌:                 int cip = -1;
    // RDKit✔️❌:                 int _cip = -1;
    // RDKit✔️❌:                 for (auto bond : m->atomBonds(m->getAtomWithIdx(ref))) {
    // RDKit✔️❌:                   if (bond->getBondType() != Bond::SINGLE &&
    // RDKit✔️❌:                       bond->getBondType() != Bond::AROMATIC) {
    // RDKit✔️❌:                     continue;
    // RDKit✔️❌:                   }
    // RDKit✔️❌:                   auto atom = bond->getOtherAtomIdx(ref);
    // RDKit✔️❌:                   if ((_cip = ranks[atom]) > cip) {
    // RDKit✔️❌:                     if (nbr >= 0) {
    // RDKit✔️❌:                       extraNbr = nbr;
    // RDKit✔️❌:                     }
    // RDKit✔️❌:                     nbr = atom;
    // RDKit✔️❌:                     cip = _cip;
    // RDKit✔️❌:                   } else {
    // RDKit✔️❌:                     extraNbr = atom;
    // RDKit✔️❌:                   }
    // RDKit✔️❌:                 }
    // RDKit✔️❌:                 return std::make_pair(nbr, extraNbr);
    // RDKit✔️❌:               };
    // RDKit✔️❌:               auto [leftNbr, extraLeftNbr] = findNbrAtoms(left);
    // RDKit✔️❌:               auto [rightNbr, extraRightNbr] = findNbrAtoms(right);

    // RDKit✔️❌:               if (leftNbr < 0 || rightNbr < 0) {
    // RDKit✔️❌:                 BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:                     << "Ignoring stereochemistry on double-bond without appropriate neighbors"
    // RDKit✔️❌:                     << std::endl;
    // RDKit✔️❌:                 continue;
    // RDKit✔️❌:               }

    // RDKit✔️❌:               bool switchEZ = false;
    // RDKit✔️❌:               if ((originalLeftNbr == leftNbr &&
    // RDKit✔️❌:                    originalRightNbr != rightNbr) ||
    // RDKit✔️❌:                   (originalLeftNbr != leftNbr &&
    // RDKit✔️❌:                    originalRightNbr == rightNbr)) {
    // RDKit✔️❌:                 switchEZ = true;
    // RDKit✔️❌:               }

    // RDKit✔️❌:               char parity = stereo0DPtr->parity;
    // RDKit✔️❌:               if (parity == INCHI_PARITY_ODD && switchEZ) {
    // RDKit✔️❌:                 parity = INCHI_PARITY_EVEN;
    // RDKit✔️❌:               } else if (parity == INCHI_PARITY_EVEN && switchEZ) {
    // RDKit✔️❌:                 parity = INCHI_PARITY_ODD;
    // RDKit✔️❌:               }

    // RDKit✔️❌:               auto findBondPairs = [&m, &zBondPairs, &eBondPairs](
    // RDKit✔️❌:                                        unsigned ref, int nbr, int extraNbr) {
    // RDKit✔️❌:                 auto bond = m->getBondBetweenAtoms(ref, nbr);
    // RDKit✔️❌:                 if (extraNbr >= 0) {
    // RDKit✔️❌:                   // modifier to track whether bond is reversed
    // RDKit✔️❌:                   int modifier = -1;
    // RDKit✔️❌:                   if (bond->getBeginAtomIdx() != ref) {
    // RDKit✔️❌:                     modifier *= -1;
    // RDKit✔️❌:                   }
    // RDKit✔️❌:                   auto extraBond = m->getBondBetweenAtoms(ref, extraNbr);
    // RDKit✔️❌:                   if (extraBond->getBeginAtomIdx() != ref) {
    // RDKit✔️❌:                     modifier *= -1;
    // RDKit✔️❌:                   }
    // RDKit✔️❌:                   if (modifier == 1) {
    // RDKit✔️❌:                     zBondPairs.push_back(
    // RDKit✔️❌:                         std::make_pair(bond->getIdx(), extraBond->getIdx()));
    // RDKit✔️❌:                   } else {
    // RDKit✔️❌:                     eBondPairs.push_back(
    // RDKit✔️❌:                         std::make_pair(bond->getIdx(), extraBond->getIdx()));
    // RDKit✔️❌:                   }
    // RDKit✔️❌:                 }
    // RDKit✔️❌:                 return bond;
    // RDKit✔️❌:               };

    // RDKit✔️❌:               auto leftBond = findBondPairs(left, leftNbr, extraLeftNbr);
    // RDKit✔️❌:               auto rightBond = findBondPairs(right, rightNbr, extraRightNbr);

    // RDKit✔️❌:               int modifier = -1;  // modifier to track whether bond is reversed
    // RDKit✔️❌:               if (leftBond->getBeginAtomIdx() != left) {
    // RDKit✔️❌:                 modifier *= -1;
    // RDKit✔️❌:               }
    // RDKit✔️❌:               if (rightBond->getBeginAtomIdx() != right) {
    // RDKit✔️❌:                 modifier *= -1;
    // RDKit✔️❌:               }

    // RDKit✔️❌:               if (parity == INCHI_PARITY_ODD) {
    // RDKit✔️❌:                 bond->setStereo(Bond::STEREOZ);
    // RDKit✔️❌:                 if (modifier == 1) {
    // RDKit✔️❌:                   eBondPairs.push_back(
    // RDKit✔️❌:                       std::make_pair(leftBond->getIdx(), rightBond->getIdx()));
    // RDKit✔️❌:                 } else {
    // RDKit✔️❌:                   zBondPairs.push_back(
    // RDKit✔️❌:                       std::make_pair(leftBond->getIdx(), rightBond->getIdx()));
    // RDKit✔️❌:                 }
    // RDKit✔️❌:               } else if (parity == INCHI_PARITY_EVEN) {
    // RDKit✔️❌:                 bond->setStereo(Bond::STEREOE);
    // RDKit✔️❌:                 if (modifier == 1) {
    // RDKit✔️❌:                   zBondPairs.push_back(
    // RDKit✔️❌:                       std::make_pair(leftBond->getIdx(), rightBond->getIdx()));
    // RDKit✔️❌:                 } else {
    // RDKit✔️❌:                   eBondPairs.push_back(
    // RDKit✔️❌:                       std::make_pair(leftBond->getIdx(), rightBond->getIdx()));
    // RDKit✔️❌:                 }
    // RDKit✔️❌:               } else if (parity == INCHI_PARITY_NONE) {
    // RDKit✔️❌:                 bond->setStereo(Bond::STEREONONE);
    // RDKit✔️❌:               } else {
    // RDKit✔️❌:                 bond->setStereo(Bond::STEREOANY);
    // RDKit✔️❌:               }
    // RDKit✔️❌:               // set the stereo atoms for the double bond
    // RDKit✔️❌:               bond->getStereoAtoms().push_back(leftNbr);
    // RDKit✔️❌:               bond->getStereoAtoms().push_back(rightNbr);
    // RDKit✔️❌:               break;
    // RDKit✔️❌:             }
    // RDKit✔️❌:             case INCHI_StereoType_Tetrahedral: {
    // RDKit✔️❌:               unsigned int c =
    // RDKit✔️❌:                   indexToAtomIndexMapping[stereo0DPtr->central_atom];
    // RDKit✔️❌:               Atom *atom = m->getAtomWithIdx(c);
    // RDKit✔️❌:               // find number of swaps for the members
    // RDKit✔️❌:               int nSwaps = 0;
    // RDKit✔️❌:               unsigned int nid = 0;
    // RDKit✔️❌:               if (stereo0DPtr->neighbor[0] == stereo0DPtr->central_atom) {
    // RDKit✔️❌:                 // 3-neighbor case
    // RDKit✔️❌:                 nid = 1;
    // RDKit✔️❌:                 if (atom->getDegree() == 3) {
    // RDKit✔️❌:                   // this happens with chiral three-coordinate S
    // RDKit✔️❌:                   nSwaps = 1;
    // RDKit✔️❌:                 }
    // RDKit✔️❌:               }
    // RDKit✔️❌:               // if (atom->getTotalNumHs(true) == 1)
    // RDKit✔️❌:               //  nSwaps = 1;
    // RDKit✔️❌:               // std::cerr<<"build atom: "<<c<<" "<<atom->getTotalNumHs(true);
    // RDKit✔️❌:               std::list<int> neighbors;
    // RDKit✔️❌:               for (; nid < 4; nid++) {
    // RDKit✔️❌:                 unsigned end =
    // RDKit✔️❌:                     indexToAtomIndexMapping[stereo0DPtr->neighbor[nid]];
    // RDKit✔️❌:                 Bond *bond = m->getBondBetweenAtoms(c, end);
    // RDKit✔️❌:                 neighbors.push_back(bond->getIdx());
    // RDKit✔️❌:                 // std::cerr<<" "<<end<<"("<<bond->getIdx()<<")";
    // RDKit✔️❌:               }
    // RDKit✔️❌:               nSwaps += atom->getPerturbationOrder(neighbors);
    // RDKit✔️❌:               // std::cerr<<" swaps: "<<nSwaps<<" parity: "<<
    // RDKit✔️❌:               //  (stereo0DPtr->parity==INCHI_PARITY_EVEN?"even":"odd")<<std::endl;
    // RDKit✔️❌:               if (stereo0DPtr->parity == INCHI_PARITY_ODD) {
    // RDKit✔️❌:                 atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️❌:               } else {
    // RDKit✔️❌:                 atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    // RDKit✔️❌:               }
    // RDKit✔️❌:               if (nSwaps % 2) {
    // RDKit✔️❌:                 atom->invertChirality();
    // RDKit✔️❌:               }
    // RDKit✔️❌:               break;
    // RDKit✔️❌:             }
    // RDKit✔️❌:             case INCHI_StereoType_Allene:
    // RDKit✔️❌:               BOOST_LOG(rdWarningLog) << "Allene-style stereochemistry is not "
    // RDKit✔️❌:                                          "supported yet and will be ignored."
    // RDKit✔️❌:                                       << std::endl;
    // RDKit✔️❌:               break;
    // RDKit✔️❌:             default:
    // RDKit✔️❌:               BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:                   << "Unrecognized stereo0D type (" << (int)stereo0DPtr->type
    // RDKit✔️❌:                   << ") is ignored!" << std::endl;
    // RDKit✔️❌:           }  // end switch stereotype
    // RDKit✔️❌:         }  // end for loop over all stereo0D entries
    // RDKit✔️❌:         // set the bond directions
    // RDKit✔️❌:         if (!assignBondDirs(*m, zBondPairs, eBondPairs)) {
    // RDKit✔️❌:           BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:               << "Cannot assign bond directions!" << std::endl;
    // RDKit✔️❌:           ;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }  // end if (if stereo0D presents)
    // RDKit✔️❌:     }  // end if (if return code is success)

    // RDKit✔️❌:     // clean up
    // RDKit✔️❌:     FreeStructFromINCHI(&inchiOutput);
    // RDKit✔️❌:   }

    // RDKit✔️❌:   // clean up the molecule to be acceptable to RDKit
    // RDKit✔️❌:   if (m) {
    // RDKit✔️❌:     cleanUp(*m);
    // RDKit✔️❌:     try {
    // RDKit✔️❌:       if (sanitize) {
    // RDKit✔️❌:         if (removeHs) {
    // RDKit✔️❌:           MolOps::removeHs(*m);
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           MolOps::sanitizeMol(*m);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } catch (const MolSanitizeException &) {
    // RDKit✔️❌:       delete m;
    // RDKit✔️❌:       throw;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // call assignStereochemistry just to be safe; otherwise, MolToSmiles may
    // RDKit✔️❌:     // overwrite E/Z and/or bond direction on double bonds.
    // RDKit✔️❌:     MolOps::assignStereochemistry(*m, true, true);
    // RDKit✔️❌:   }

    // RDKit✔️❌:   return m;
    // RDKit✔️❌: }

    // END RDKIT C++ FUNCTION: InchiToMol
    // BEGIN RDKIT ACTIVE CONFIGURATION: InchiToMol
    // RDKit✔️❌: Pinned RDKit 2026.03.1; GCC/Linux; DEBUG is not defined.
    // RDKit✔️❌: The neutral engine and toolkit traits preserve the source call boundary without
    // RDKit✔️❌: linking RDKit or official C in production. SourceHeap output cloning and checked
    // RDKit✔️❌: graph access are known to cost more than the native pointer implementation.
    // END RDKIT ACTIVE CONFIGURATION: InchiToMol

    let mut inchi_with_nul = Vec::with_capacity(inchi.len() + 1);
    inchi_with_nul.extend_from_slice(inchi);
    inchi_with_nul.push(0);

    let output = engine.get_struct_from_inchi(&inchi_with_nul)?;
    return_values.return_code = output.return_code;
    if let Some(message) = &output.message {
        return_values.message.clone_from(message);
    }
    if let Some(log) = &output.log {
        return_values.log.clone_from(log);
    }

    let mut diagnostics = Vec::new();
    let mut molecule = None;
    if output.return_code == tagRetValGetINCHI_inchi_Ret_OKAY
        || output.return_code == tagRetValGetINCHI_inchi_Ret_WARNING
    {
        let mut built = AdapterMol::from_graph(Vec::new(), Vec::new());
        let mut isotopes = Vec::<(u32, u32, u32)>::new();
        let mut index_to_atom_index_mapping = Vec::with_capacity(output.atoms.len());

        for (index, inchi_atom) in output.atoms.iter().enumerate() {
            let element = adapter_element_name(inchi_atom)?;
            let atomic_number = toolkit.atomic_number(&element)?;
            let average_weight = toolkit.average_atomic_weight(atomic_number)?;
            let reference_weight = (average_weight + 0.5) as i32;
            let mut isotope = 0_i32;
            if inchi_atom.isotopic_mass != 0 {
                isotope = i32::from(inchi_atom.isotopic_mass) - ISOTOPIC_SHIFT_FLAG as i32;
            }

            let mut atom = AdapterAtom {
                atomic_number,
                formal_charge: i32::from(inchi_atom.charge),
                num_explicit_hydrogens: inchi_atom.num_iso_H[0] as u32,
                isotope: if isotope != 0 {
                    isotope.wrapping_add(reference_weight) as u32
                } else {
                    0
                },
                no_implicit: true,
                ..AdapterAtom::default()
            };
            if inchi_atom.radical != 0 {
                if inchi_atom.radical != 3 && inchi_atom.radical != 2 {
                    diagnostics.push(AdapterDiagnostic {
                        level: AdapterDiagnosticLevel::Warning,
                        message: format!(
                            "expect radical to be either 2 or 3 while getting {}. Ignore radical.",
                            char::from(inchi_atom.radical as u8)
                        ),
                    });
                } else {
                    atom.num_radical_electrons = (inchi_atom.radical - 1) as u32;
                }
            }

            if inchi_atom.num_iso_H[1] != 0 {
                isotopes.push((1, index as u32, inchi_atom.num_iso_H[1] as u32));
            } else if inchi_atom.num_iso_H[2] != 0 {
                isotopes.push((2, index as u32, inchi_atom.num_iso_H[2] as u32));
            } else if inchi_atom.num_iso_H[3] != 0 {
                isotopes.push((3, index as u32, inchi_atom.num_iso_H[3] as u32));
            }
            let atom_index = adapter_add_atom(&mut built, atom);
            index_to_atom_index_mapping.push(atom_index);
        }

        let mut bond_register = BTreeSet::<(u32, u32)>::new();
        for (index, inchi_atom) in output.atoms.iter().enumerate() {
            let bond_count = inchi_atom.num_bonds as u32;
            for bond_offset in 0..bond_count {
                let offset = usize::try_from(bond_offset)
                    .map_err(|_| InchiToMolError::InvalidSourceOutput("bond offset exceeds usize"))?;
                if offset >= inchi_atom.neighbor.len() {
                    return Err(InchiToMolError::InvalidSourceOutput(
                        "InChI atom bond count exceeds fixed source array",
                    ));
                }
                let neighbor = inchi_atom.neighbor[offset] as u32;
                let index = index as u32;
                if bond_register.contains(&(index, neighbor)) || bond_register.contains(&(neighbor, index)) {
                    continue;
                }
                bond_register.insert((index, neighbor));

                let raw_bond_type = inchi_atom.bond_type[offset] as u32;
                let (bond_type, is_aromatic) = if raw_bond_type <= tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE as u32 {
                    (
                        match raw_bond_type {
                            0 => BondType::Unspecified,
                            1 => BondType::Single,
                            2 => BondType::Double,
                            3 => BondType::Triple,
                            _ => unreachable!("bounded by source comparison"),
                        },
                        false,
                    )
                } else if raw_bond_type == tagINCHIBondType_INCHI_BOND_TYPE_ALTERN as u32 {
                    diagnostics.push(AdapterDiagnostic {
                        level: AdapterDiagnosticLevel::Warning,
                        message: "receive ALTERN bond type which should be avoided. This is treated as aromatic."
                            .to_owned(),
                    });
                    (BondType::Aromatic, true)
                } else {
                    diagnostics.push(AdapterDiagnostic {
                        level: AdapterDiagnosticLevel::Error,
                        message: format!("illegal bond type ({raw_bond_type}) in InChI"),
                    });
                    engine.free_struct_from_inchi()?;
                    return Ok(InchiToMolResult {
                        molecule: None,
                        diagnostics,
                    });
                };

                let begin_atom_index =
                    *index_to_atom_index_mapping
                        .get(index as usize)
                        .ok_or(InchiToMolError::InvalidSourceOutput(
                            "source atom index is outside atom mapping",
                        ))?;
                let end_atom_index =
                    *index_to_atom_index_mapping
                        .get(neighbor as usize)
                        .ok_or(InchiToMolError::InvalidSourceOutput(
                            "source neighbor index is outside atom mapping",
                        ))?;
                let raw_stereo = i32::from(inchi_atom.bond_stereo[offset]);
                let direction = match raw_stereo {
                    value if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE => BondDirection::None,
                    value
                        if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1UP
                            || value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2DOWN =>
                    {
                        BondDirection::BeginWedge
                    }
                    value
                        if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1DOWN
                            || value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2UP =>
                    {
                        BondDirection::BeginDash
                    }
                    value if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1EITHER => BondDirection::Unknown,
                    value if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER => {
                        BondDirection::EitherDouble
                    }
                    _ => BondDirection::None,
                };
                adapter_add_bond(
                    &mut built,
                    AdapterBond {
                        begin_atom_index,
                        end_atom_index,
                        bond_type,
                        direction,
                        is_aromatic,
                        ..AdapterBond::default()
                    },
                );
            }
        }

        for (isotope, source_atom_index, repeat) in isotopes {
            let attached_atom_index = *index_to_atom_index_mapping.get(source_atom_index as usize).ok_or(
                InchiToMolError::InvalidSourceOutput("isotope-H parent index is outside atom mapping"),
            )?;
            for _ in 0..repeat {
                let hydrogen_index = adapter_add_atom(
                    &mut built,
                    AdapterAtom {
                        atomic_number: 1,
                        isotope,
                        ..AdapterAtom::default()
                    },
                );
                adapter_add_bond(
                    &mut built,
                    AdapterBond {
                        begin_atom_index: hydrogen_index,
                        end_atom_index: attached_atom_index,
                        bond_type: BondType::Single,
                        ..AdapterBond::default()
                    },
                );
            }
        }

        toolkit.update_property_cache(&mut built, false)?;

        let mut e_bond_pairs = Vec::<(i32, i32)>::new();
        let mut z_bond_pairs = Vec::<(i32, i32)>::new();
        if !output.stereo0d.is_empty() {
            let ranks = toolkit.assign_atom_cip_ranks(&mut built)?;
            if ranks.len() < built.atoms.len() {
                return Err(InchiToMolError::InvalidSourceOutput(
                    "toolkit CIP rank vector is shorter than the molecule",
                ));
            }
            for (atom, rank) in built.atoms.iter_mut().zip(ranks.iter().copied()) {
                atom.cip_rank = Some(rank);
            }

            for stereo in &output.stereo0d {
                let parity = i32::from(stereo.parity);
                if parity == tagINCHIStereoParity0D_INCHI_PARITY_NONE as i32
                    || parity == tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED as i32
                {
                    continue;
                }
                match i32::from(stereo.type_) {
                    value if value == tagINCHIStereoType0D_INCHI_StereoType_None as i32 => {}
                    value if value == tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i32 => {
                        let mapped = |source_index: i16| {
                            index_to_atom_index_mapping
                                .get(source_index as u16 as usize)
                                .copied()
                                .ok_or(InchiToMolError::InvalidSourceOutput(
                                    "double-bond stereo atom index is outside atom mapping",
                                ))
                        };
                        let left = mapped(stereo.neighbor[1])?;
                        let right = mapped(stereo.neighbor[2])?;
                        let original_left_neighbor = mapped(stereo.neighbor[0])? as i32;
                        let original_right_neighbor = mapped(stereo.neighbor[3])? as i32;

                        let Some(double_bond_index) = bond_index_between(&built, left, right) else {
                            diagnostics.push(AdapterDiagnostic {
                                level: AdapterDiagnosticLevel::Warning,
                                message: "Extended double-bond stereochemistry (e.g. C=C=C=C) ignored".to_owned(),
                            });
                            continue;
                        };

                        let find_neighbors = |molecule: &AdapterMol, reference: u32| -> (i32, i32) {
                            let mut neighbor = -1_i32;
                            let mut extra_neighbor = -1_i32;
                            let mut cip = -1_i32;
                            for &(other_atom, bond_index) in &molecule.adjacency[reference as usize] {
                                let bond = &molecule.bonds[bond_index as usize];
                                if bond.bond_type != BondType::Single && bond.bond_type != BondType::Aromatic {
                                    continue;
                                }
                                // GCC converts the unsigned rank to the source `int _cip`
                                // before comparison.
                                let candidate_cip = ranks[other_atom as usize] as i32;
                                if candidate_cip > cip {
                                    if neighbor >= 0 {
                                        extra_neighbor = neighbor;
                                    }
                                    neighbor = other_atom as i32;
                                    cip = candidate_cip;
                                } else {
                                    extra_neighbor = other_atom as i32;
                                }
                            }
                            (neighbor, extra_neighbor)
                        };
                        let (left_neighbor, extra_left_neighbor) = find_neighbors(&built, left);
                        let (right_neighbor, extra_right_neighbor) = find_neighbors(&built, right);
                        if left_neighbor < 0 || right_neighbor < 0 {
                            diagnostics.push(AdapterDiagnostic {
                                level: AdapterDiagnosticLevel::Warning,
                                message: "Ignoring stereochemistry on double-bond without appropriate neighbors"
                                    .to_owned(),
                            });
                            continue;
                        }

                        let switch_ez = (original_left_neighbor == left_neighbor
                            && original_right_neighbor != right_neighbor)
                            || (original_left_neighbor != left_neighbor && original_right_neighbor == right_neighbor);
                        let mut effective_parity = parity;
                        if effective_parity == tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32 && switch_ez {
                            effective_parity = tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i32;
                        } else if effective_parity == tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i32 && switch_ez {
                            effective_parity = tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32;
                        }

                        let mut find_bond_pairs =
                            |reference: u32, neighbor: i32, extra_neighbor: i32| -> Result<u32, InchiToMolError> {
                                let selected = bond_index_between(&built, reference, neighbor as u32).ok_or(
                                    InchiToMolError::InvalidSourceOutput("selected stereo neighbor has no bond"),
                                )?;
                                if extra_neighbor >= 0 {
                                    let mut modifier = -1_i32;
                                    if built.bonds[selected as usize].begin_atom_index != reference {
                                        modifier *= -1;
                                    }
                                    let extra = bond_index_between(&built, reference, extra_neighbor as u32).ok_or(
                                        InchiToMolError::InvalidSourceOutput("extra stereo neighbor has no bond"),
                                    )?;
                                    if built.bonds[extra as usize].begin_atom_index != reference {
                                        modifier *= -1;
                                    }
                                    if modifier == 1 {
                                        z_bond_pairs.push((selected as i32, extra as i32));
                                    } else {
                                        e_bond_pairs.push((selected as i32, extra as i32));
                                    }
                                }
                                Ok(selected)
                            };
                        let left_bond = find_bond_pairs(left, left_neighbor, extra_left_neighbor)?;
                        let right_bond = find_bond_pairs(right, right_neighbor, extra_right_neighbor)?;

                        let mut modifier = -1_i32;
                        if built.bonds[left_bond as usize].begin_atom_index != left {
                            modifier *= -1;
                        }
                        if built.bonds[right_bond as usize].begin_atom_index != right {
                            modifier *= -1;
                        }

                        let bond = &mut built.bonds[double_bond_index as usize];
                        if effective_parity == tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32 {
                            bond.stereo = BondStereo::Z;
                            if modifier == 1 {
                                e_bond_pairs.push((left_bond as i32, right_bond as i32));
                            } else {
                                z_bond_pairs.push((left_bond as i32, right_bond as i32));
                            }
                        } else if effective_parity == tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i32 {
                            bond.stereo = BondStereo::E;
                            if modifier == 1 {
                                z_bond_pairs.push((left_bond as i32, right_bond as i32));
                            } else {
                                e_bond_pairs.push((left_bond as i32, right_bond as i32));
                            }
                        } else if effective_parity == tagINCHIStereoParity0D_INCHI_PARITY_NONE as i32 {
                            bond.stereo = BondStereo::None;
                        } else {
                            bond.stereo = BondStereo::Any;
                        }
                        bond.stereo_atoms.push(left_neighbor as u32);
                        bond.stereo_atoms.push(right_neighbor as u32);
                    }
                    value if value == tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i32 => {
                        let central = *index_to_atom_index_mapping
                            .get(stereo.central_atom as u16 as usize)
                            .ok_or(InchiToMolError::InvalidSourceOutput(
                                "tetrahedral central atom is outside atom mapping",
                            ))?;
                        let mut swaps = 0_u32;
                        let mut neighbor_position = 0_usize;
                        if stereo.neighbor[0] == stereo.central_atom {
                            neighbor_position = 1;
                            if built.adjacency[central as usize].len() == 3 {
                                swaps = 1;
                            }
                        }
                        let mut neighbors = Vec::new();
                        while neighbor_position < 4 {
                            let end = *index_to_atom_index_mapping
                                .get(stereo.neighbor[neighbor_position] as u16 as usize)
                                .ok_or(InchiToMolError::InvalidSourceOutput(
                                    "tetrahedral neighbor is outside atom mapping",
                                ))?;
                            let bond = bond_index_between(&built, central, end).ok_or(
                                InchiToMolError::InvalidSourceOutput("tetrahedral neighbor has no central bond"),
                            )?;
                            neighbors.push(bond);
                            neighbor_position += 1;
                        }
                        swaps += adapter_perturbation_order(&built, central, &neighbors)?;
                        built.atoms[central as usize].chiral_tag =
                            if parity == tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32 {
                                ChiralTag::TetrahedralCcw
                            } else {
                                ChiralTag::TetrahedralCw
                            };
                        if swaps % 2 != 0 {
                            adapter_invert_chirality(&mut built.atoms[central as usize]);
                        }
                    }
                    value if value == tagINCHIStereoType0D_INCHI_StereoType_Allene as i32 => {
                        diagnostics.push(AdapterDiagnostic {
                            level: AdapterDiagnosticLevel::Warning,
                            message: "Allene-style stereochemistry is not supported yet and will be ignored."
                                .to_owned(),
                        });
                    }
                    value => {
                        diagnostics.push(AdapterDiagnostic {
                            level: AdapterDiagnosticLevel::Warning,
                            message: format!("Unrecognized stereo0D type ({value}) is ignored!"),
                        });
                    }
                }
            }
            if !assign_bond_dirs(&mut built, &z_bond_pairs, &e_bond_pairs)? {
                diagnostics.push(AdapterDiagnostic {
                    level: AdapterDiagnosticLevel::Warning,
                    message: "Cannot assign bond directions!".to_owned(),
                });
            }
        }
        molecule = Some(built);
    }

    engine.free_struct_from_inchi()?;

    if let Some(mut built) = molecule {
        clean_up(&mut built)?;
        toolkit.synchronize_after_cleanup(&built)?;
        if sanitize {
            if remove_hs {
                toolkit.remove_hydrogens(&mut built)?;
            } else {
                toolkit.sanitize_molecule(&mut built)?;
            }
        }
        toolkit.assign_stereochemistry(&mut built, true, true)?;
        molecule = Some(built);
    }

    Ok(InchiToMolResult { molecule, diagnostics })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FixOptionSymbolError {
    InputIsNotNulTerminated,
    OutputIsTooSmall,
}

pub(crate) fn fix_option_symbol(input: &[u8], output: &mut [u8]) -> Result<(), FixOptionSymbolError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:1674 fixOptionSymbol
    // RDKit✔️✔️: void fixOptionSymbol(const char *in, char *out) {
    // RDKit✔️✔️:   unsigned int i;
    // RDKit✔️✔️:   for (i = 0; i < strlen(in); i++) {
    // RDKit✔️✔️: #ifdef _WIN32
    // RDKit✔️✔️:     if (in[i] == '-') {
    // RDKit✔️✔️:       out[i] = '/';
    // RDKit✔️✔️:
    // RDKit✔️✔️: #else
    // RDKit✔️✔️:     if (in[i] == '/') {
    // RDKit✔️✔️:       out[i] = '-';
    // RDKit✔️✔️:
    // RDKit✔️✔️: #endif
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       out[i] = in[i];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   out[i] = '\0';
    // RDKit✔️✔️: }
    // END RDKIT C++ FUNCTION: fixOptionSymbol
    // BEGIN RDKIT ACTIVE CONFIGURATION: fixOptionSymbol
    // RDKit✔️✔️: Pinned RDKit 2026.03.1; GCC/Linux; `_WIN32` is not defined.
    // RDKit✔️✔️: The production callers provide distinct input/output buffers, a terminating NUL,
    // RDKit✔️✔️: and `strlen(in) + 1` writable output bytes. The inactive Windows branch is retained
    // RDKit✔️✔️: only in the verbatim source frame. Rust validates those caller preconditions.
    // END RDKIT ACTIVE CONFIGURATION: fixOptionSymbol

    let input_len = input
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(FixOptionSymbolError::InputIsNotNulTerminated)?;
    if output.len() <= input_len {
        return Err(FixOptionSymbolError::OutputIsTooSmall);
    }
    for (index, byte) in input[..input_len].iter().copied().enumerate() {
        output[index] = if byte == b'/' { b'-' } else { byte };
    }
    output[input_len] = 0;
    Ok(())
}

fn r_cleanup_matches(molecule: &AdapterMol) -> Vec<[u32; 5]> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/Code/GraphMol/Atom.cpp:706 Atom::Match
    // RDKit✔️❌: bool Atom::Match(Atom const *what) const {
    // RDKit✔️❌:   PRECONDITION(what, "bad query atom");
    // RDKit✔️❌:   bool res = getAtomicNum() == what->getAtomicNum();
    // RDKit✔️❌:
    // RDKit✔️❌:   // special dummy--dummy match case:
    // RDKit✔️❌:   //   [*] matches [*],[1*],[2*],etc.
    // RDKit✔️❌:   //   [1*] only matches [*] and [1*]
    // RDKit✔️❌:   if (res) {
    // RDKit✔️❌:     if (!this->getAtomicNum()) {
    // RDKit✔️❌:       // this is the new behavior, based on the isotopes:
    // RDKit✔️❌:       int tgt = this->getIsotope();
    // RDKit✔️❌:       int test = what->getIsotope();
    // RDKit✔️❌:       if (tgt && test && tgt != test) {
    // RDKit✔️❌:         res = false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       // standard atom-atom match: The general rule here is that if this atom
    // RDKit✔️❌:       // has a property that
    // RDKit✔️❌:       // deviates from the default, then the other atom should match that value.
    // RDKit✔️❌:       if ((this->getFormalCharge() &&
    // RDKit✔️❌:            this->getFormalCharge() != what->getFormalCharge()) ||
    // RDKit✔️❌:           (this->getIsotope() && this->getIsotope() != what->getIsotope()) ||
    // RDKit✔️❌:           (this->getNumRadicalElectrons() &&
    // RDKit✔️❌:            this->getNumRadicalElectrons() != what->getNumRadicalElectrons())) {
    // RDKit✔️❌:         res = false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // END RDKIT C++ FUNCTION: Atom::Match

    let mut seen_atom_sets = HashSet::<[u32; 5]>::new();
    let mut matches = Vec::new();

    for query_0 in 0..molecule.atoms.len() as u32 {
        if molecule.atoms[query_0 as usize].atomic_number != 8 || molecule.atoms[query_0 as usize].formal_charge != -1 {
            continue;
        }
        let mut central_atoms = molecule.adjacency[query_0 as usize]
            .iter()
            .filter_map(|&(neighbor, bond)| {
                (molecule.atoms[neighbor as usize].atomic_number == 17
                    && molecule.atoms[neighbor as usize].formal_charge == 3
                    && molecule.bonds[bond as usize].bond_type == BondType::Single)
                    .then_some(neighbor)
            })
            .collect::<Vec<_>>();
        central_atoms.sort_unstable();
        for query_1 in central_atoms {
            let mut oxygen_atoms = molecule.adjacency[query_1 as usize]
                .iter()
                .filter_map(|&(neighbor, bond)| {
                    (neighbor != query_0
                        && molecule.atoms[neighbor as usize].atomic_number == 8
                        && molecule.bonds[bond as usize].bond_type == BondType::Single)
                        .then_some(neighbor)
                })
                .collect::<Vec<_>>();
            oxygen_atoms.sort_unstable();
            for &query_2 in &oxygen_atoms {
                if molecule.atoms[query_2 as usize].formal_charge != -1 {
                    continue;
                }
                for &query_3 in &oxygen_atoms {
                    if query_3 == query_2 || molecule.atoms[query_3 as usize].formal_charge != -1 {
                        continue;
                    }
                    for &query_4 in &oxygen_atoms {
                        if query_4 == query_2 || query_4 == query_3 {
                            continue;
                        }
                        let mapping = [query_0, query_1, query_2, query_3, query_4];
                        let mut atom_set = mapping;
                        atom_set.sort_unstable();
                        if seen_atom_sets.insert(atom_set) {
                            matches.push(mapping);
                            if matches.len() == 1000 {
                                return matches;
                            }
                        }
                    }
                }
            }
        }
    }
    matches
}

pub(crate) fn r_clean_up(molecule: &mut AdapterMol) {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:1694 rCleanUp
    // RDKit✔️❌: void rCleanUp(RWMol &mol) {
    // RDKit✔️❌:   RWMol *q = SmilesToMol("[O-][Cl+3]([O-])([O-])O");
    // RDKit✔️❌:   std::vector<MatchVectType> fgpMatches;
    // RDKit✔️❌:   SubstructMatch(mol, *q, fgpMatches);
    // RDKit✔️❌:   delete q;
    // RDKit✔️❌:   // replace all matches
    // RDKit✔️❌:   for (auto match : fgpMatches) {
    // RDKit✔️❌:     // collect matching atoms
    // RDKit✔️❌:     int map[5];
    // RDKit✔️❌:     for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
    // RDKit✔️❌:          mi++) {
    // RDKit✔️❌:       map[mi->first] = mi->second;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // check charges
    // RDKit✔️❌:     if (mol.getAtomWithIdx(map[1])->getFormalCharge() != 3) {
    // RDKit✔️❌:       return;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     int unchargedFound = -1;
    // RDKit✔️❌:     for (int i = 0; i < 5; i++) {
    // RDKit✔️❌:       if (i == 1) {
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       Atom *o = mol.getAtomWithIdx(map[i]);
    // RDKit✔️❌:       if (o->getFormalCharge() == 0) {
    // RDKit✔️❌:         if (unchargedFound != -1) {
    // RDKit✔️❌:           return;  // too many uncharged oxygen
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           unchargedFound = i;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     // flip bonds and remove charges
    // RDKit✔️❌:     for (int i = 0; i < 5; i++) {
    // RDKit✔️❌:       if (i == 1) {
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (i == unchargedFound) {
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (unchargedFound == -1 && i == 0) {
    // RDKit✔️❌:         mol.getBondBetweenAtoms(map[1], map[i])->setBondType(Bond::SINGLE);
    // RDKit✔️❌:         mol.getAtomWithIdx(map[i])->setFormalCharge(-1);
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       mol.getBondBetweenAtoms(map[1], map[i])->setBondType(Bond::DOUBLE);
    // RDKit✔️❌:       mol.getAtomWithIdx(map[i])->setFormalCharge(0);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     mol.getAtomWithIdx(map[1])->setFormalCharge(0);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return;
    // RDKit✔️❌: }
    // END RDKIT C++ FUNCTION: rCleanUp
    // BEGIN RDKIT ACTIVE CONFIGURATION: rCleanUp
    // RDKit✔️❌: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️❌: `SmilesToMol` produces plain atoms. `Atom::Match` requires atomic number and each
    // RDKit✔️❌: non-default query charge, isotope, and radical value; therefore the three `O-`
    // RDKit✔️❌: atoms and `Cl+3` match charges exactly, while the final neutral query oxygen's
    // RDKit✔️❌: default zero charge is a wildcard. Bond matching requires exact single bonds.
    // RDKit✔️❌: Default unique-by-target-atom-set and 1000-result `SubstructMatch` completes before
    // RDKit✔️❌: this function performs its retained post-match charge checks and mutations.
    // RDKit✔️❌: The dedicated safe matcher preserves source ordering but has known extra temporary
    // RDKit✔️❌: allocation and repeated sorting compared with RDKit's substructure matcher.
    // END RDKIT ACTIVE CONFIGURATION: rCleanUp

    let matches = r_cleanup_matches(molecule);
    for mapping in matches {
        if molecule.atoms[mapping[1] as usize].formal_charge != 3 {
            return;
        }
        let mut uncharged_found = None;
        for index in [0_usize, 2, 3, 4] {
            if molecule.atoms[mapping[index] as usize].formal_charge == 0 {
                if uncharged_found.is_some() {
                    return;
                }
                uncharged_found = Some(index);
            }
        }

        for index in [0_usize, 2, 3, 4] {
            if uncharged_found == Some(index) {
                continue;
            }
            let bond_index = bond_index_between(molecule, mapping[1], mapping[index])
                .expect("rCleanUp match contains every query bond");
            if uncharged_found.is_none() && index == 0 {
                molecule.bonds[bond_index as usize].bond_type = BondType::Single;
                molecule.atoms[mapping[index] as usize].formal_charge = -1;
                continue;
            }
            molecule.bonds[bond_index as usize].bond_type = BondType::Double;
            molecule.atoms[mapping[index] as usize].formal_charge = 0;
        }
        molecule.atoms[mapping[1] as usize].formal_charge = 0;
    }
}

#[rustfmt::skip]
pub(crate) fn mol_to_inchi(
    engine: &mut impl InchiGenerationEngine,
    toolkit: &mut impl MolToInchiToolkit,
    molecule: &AdapterMol,
    return_values: &mut ExtraInchiReturnValues,
    options: Option<&[u8]>,
) -> Result<MolToInchiResult, MolToInchiError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:1747 MolToInchi
    // RDKit✔️❌: complete source frame follows verbatim.
    /*
    std::string MolToInchi(const ROMol &mol, ExtraInchiReturnValues &rv,
                           const char *options) {
      std::unique_ptr<RWMol> m{new RWMol(mol)};
    
      // assign stereochem:
      if (mol.needsUpdatePropertyCache()) {
        m->updatePropertyCache(false);
      }
      // kekulize
      MolOps::Kekulize(*m, false);
    
      // "reverse" cleanup: undo some clean up done by RDKit
      rCleanUp(*m);
    
      unsigned int nAtoms = m->getNumAtoms();
      unsigned int nBonds = m->getNumBonds();
    
      // Make array of inchi_atom (storage space)
      std::unique_ptr<inchi_Atom[]> inchiAtoms(new inchi_Atom[nAtoms]);
      // and a vector for stereo0D
      std::vector<inchi_Stereo0D> stereo0DEntries;
    
      PeriodicTable *periodicTable = PeriodicTable::getTable();
      // Fill inchi_Atom's by atoms in RWMol
      for (unsigned int i = 0; i < nAtoms; i++) {
        Atom *atom = m->getAtomWithIdx(i);
        inchiAtoms[i].num_bonds = 0;
    
        // coordinates
        if (!m->getNumConformers()) {
          inchiAtoms[i].x = 0;
          inchiAtoms[i].y = 0;
          inchiAtoms[i].z = 0;
        } else {
          auto conformerIter = m->beginConformers();
          RDGeom::Point3D coord = (*conformerIter)->getAtomPos(i);
          inchiAtoms[i].x = coord[0];
          inchiAtoms[i].y = coord[1];
          inchiAtoms[i].z = coord[2];
        }
    
        // element name
        unsigned int atomicNumber = atom->getAtomicNum();
        std::string elementName = periodicTable->getElementSymbol(atomicNumber);
        strcpy(inchiAtoms[i].elname, elementName.c_str());
    
        // isotopes
        int isotope = atom->getIsotope();
        if (isotope) {
          inchiAtoms[i].isotopic_mass =
              ISOTOPIC_SHIFT_FLAG + isotope -
              static_cast<int>(periodicTable->getAtomicWeight(atomicNumber) + 0.5);
        } else {
          // check explicit iso property. If this is set, we have a 0 offset
          // Example: CHEMBL220875
          // if (atom->getIsotope()){
          //  inchiAtoms[i].isotopic_mass = ISOTOPIC_SHIFT_FLAG + 0;
          //} else {
          inchiAtoms[i].isotopic_mass = 0;
          //}
        }
    
        // charge
        inchiAtoms[i].charge = atom->getFormalCharge();
    
        // number of iso H
        int nHs = -1;
        switch (atom->getAtomicNum()) {
          case 6:
          case 7:
          case 8:
          case 9:
          case 17:
          case 35:
          case 53:
            nHs = -1;
            break;
          default:
            nHs = atom->getTotalNumHs();
        }
        inchiAtoms[i].num_iso_H[0] = nHs;
        inchiAtoms[i].num_iso_H[1] = 0;
        inchiAtoms[i].num_iso_H[2] = 0;
        inchiAtoms[i].num_iso_H[3] = 0;
    
        // radical
        inchiAtoms[i].radical = 0;
        if (atom->getNumRadicalElectrons()) {
          // the direct specification of radicals in InChI is tricky since they use
          // the MDL representation (singlet, double, triplet) and we just have the
          // number of unpaired electrons. Instead we set the number of implicit Hs
          // here, that together with the atom identity and charge should be
          // sufficient
          inchiAtoms[i].num_iso_H[0] = atom->getTotalNumHs();
        } else {
        }
    
        // convert tetrahedral chirality info to Stereo0D
        if (atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CCW ||
            atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CW) {
          atom->calcImplicitValence();
          if (auto tval = atom->getTotalDegree(); tval < 3 || tval > 4) {
            BOOST_LOG(rdWarningLog)
                << "tetrahedral chirality on atom with <3 or >4 neighbors will be ignored."
                << std::endl;
    
            continue;
          }
          inchi_Stereo0D stereo0D;
          stereo0D.central_atom = i;
          stereo0D.type = INCHI_StereoType_Tetrahedral;
          ROMol::ADJ_ITER nbrIter, endNbrIter;
          boost::tie(nbrIter, endNbrIter) = m->getAtomNeighbors(atom);
          std::vector<std::pair<unsigned int, unsigned int>> neighbors;
          while (nbrIter != endNbrIter) {
            int cip = 0;
            // if (m->getAtomWithIdx(*nbrIter)->hasProp("_CIPRank"))
            //   m->getAtomWithIdx(*nbrIter)->getProp("_CIPRank", cip);
            neighbors.emplace_back(cip, *nbrIter);
            ++nbrIter;
          }
          // std::sort(neighbors.begin(), neighbors.end());
          unsigned char nid = 0;
          // std::cerr<<" at: "<<atom->getIdx();
          for (const auto &p : neighbors) {
            stereo0D.neighbor[nid++] = p.second;
          }
          if (nid == 3) {
            // std::cerr<<" nid==3, reorder";
            // std::cerr<<" "<<i;
            for (; nid > 0; nid--) {
              stereo0D.neighbor[nid] = stereo0D.neighbor[nid - 1];
              // std::cerr<<" "<<stereo0D.neighbor[nid];
            }
            stereo0D.neighbor[0] = i;
          }
          // std::cerr<<std::endl;
          Atom::ChiralType chiralTag;
          if ((chiralTag = atom->getChiralTag()) != Atom::CHI_UNSPECIFIED) {
            bool pushIt = false;
            if (atom->getDegree() == 4) {
              if (chiralTag == Atom::CHI_TETRAHEDRAL_CW) {
                stereo0D.parity = INCHI_PARITY_EVEN;
                pushIt = true;
              } else {
                stereo0D.parity = INCHI_PARITY_ODD;
                pushIt = true;
              }
            } else {
              // std::cerr<<"tag: "<<chiralTag<<std::endl;
              if (chiralTag == Atom::CHI_TETRAHEDRAL_CCW) {
                stereo0D.parity = INCHI_PARITY_EVEN;
                pushIt = true;
              } else if (chiralTag == Atom::CHI_TETRAHEDRAL_CW) {
                stereo0D.parity = INCHI_PARITY_ODD;
                pushIt = true;
              } else {
                BOOST_LOG(rdWarningLog)
                    << "unrecognized chirality tag (" << chiralTag << ") on atom "
                    << i << " is ignored." << std::endl;
              }
            }
            if (pushIt) {
              // this was github #296
              // with molecules like C[S@@](=O)C(C)(C)C the stereochem of the sulfur
              // from
              // the inchi comes back reversed if we don't have wedged bonds. There
              // must
              // be something with the way S stereochem is being handled that I'm
              // not
              // getting.
              // There's something of an explanation at around line 258 of
              // inchi_api.h
              // but that didn't help that much.
              // For want of a better idea, detect this pattern
              // and flip the stereochem:
              // if(atom->getAtomicNum()==16 &&
              //    atom->getDegree()==3 &&
              //    atom->getValence(Atom::ValenceType::EXPLICIT)==4){
              //   if(stereo0D.parity==INCHI_PARITY_EVEN){
              //     stereo0D.parity=INCHI_PARITY_ODD;
              //   } else if(stereo0D.parity==INCHI_PARITY_ODD){
              //     stereo0D.parity=INCHI_PARITY_EVEN;
              //   }
              // }
              stereo0DEntries.push_back(stereo0D);
            }
    
          } else {
            // std::string molParity;
            // atom->getProp("molParity", molParity);
            // if (molParity == "2") {
            //  stereo0D.parity = INCHI_PARITY_EVEN;
            //  stereo0DEntries.push_back(stereo0D);
            //} else if (molParity == "1") {
            //  stereo0D.parity = INCHI_PARITY_ODD;
            //  stereo0DEntries.push_back(stereo0D);
            //} else if (molParity == "0") {
            //  stereo0D.parity = INCHI_PARITY_NONE;
            //  stereo0DEntries.push_back(stereo0D);
            //} else if (molParity == "3") {
            //  stereo0D.parity = INCHI_PARITY_UNKNOWN;
            //  stereo0DEntries.push_back(stereo0D);
            //} else {
            //  BOOST_LOG(rdWarningLog) << "unrecognized parity on atom "
            //    << molParity << " is ignored." << std::endl;
            //}
          }
        }
      }
    
      // read bond info
      for (unsigned int i = 0; i < nBonds; i++) {
        Bond *bond = m->getBondWithIdx(i);
        unsigned int atomIndex1 = bond->getBeginAtomIdx();
        unsigned int atomIndex2 = bond->getEndAtomIdx();
        int bondDirectionModifier = 1;
        // update only for the atom having smaller index
        if (atomIndex1 > atomIndex2) {
          std::swap(atomIndex1, atomIndex2);
          bondDirectionModifier = -1;
        }
    
        // neighbor
        unsigned int idx = inchiAtoms[atomIndex1].num_bonds;
        // The InChI code has a max number of neighbors allowed:
        if (idx >= MAXVAL) {
          BOOST_LOG(rdErrorLog)
              << " atom " << atomIndex1 << " has too many bonds: " << idx
              << ". The InChI library supports at most " << MAXVAL << std::endl;
          return "";
        }
        inchiAtoms[atomIndex1].neighbor[idx] = atomIndex2;
    
        // bond type
        Bond::BondType bondType = bond->getBondType();
        if (bondType > Bond::TRIPLE) {
          BOOST_LOG(rdWarningLog) << "bond type above 3 (" << bondType
                                  << ") is treated as unspecified!" << std::endl;
          bondType = Bond::UNSPECIFIED;
        }
        inchiAtoms[atomIndex1].bond_type[idx] = bondType;
    
        // stereo
        Bond::BondDir bondDirection = bond->getBondDir();
        switch (bondDirection) {
          case Bond::BEGINWEDGE:
            inchiAtoms[atomIndex1].bond_stereo[idx] =
                bondDirectionModifier * INCHI_BOND_STEREO_SINGLE_1UP;
            break;
          case Bond::BEGINDASH:
            inchiAtoms[atomIndex1].bond_stereo[idx] =
                bondDirectionModifier * INCHI_BOND_STEREO_SINGLE_1DOWN;
            break;
          case Bond::EITHERDOUBLE:
            inchiAtoms[atomIndex1].bond_stereo[idx] =
                INCHI_BOND_STEREO_DOUBLE_EITHER;
            break;
          case Bond::UNKNOWN:
            inchiAtoms[atomIndex1].bond_stereo[idx] =
                bondDirectionModifier * INCHI_BOND_STEREO_SINGLE_1EITHER;
            break;
          case Bond::NONE:
          default:
            inchiAtoms[atomIndex1].bond_stereo[idx] = INCHI_BOND_STEREO_NONE;
        }
    
        // double bond stereochemistry
        // single bond in the big ring will get E/Z assigned as well. Though rdkit
        // will eventually remove it, I added it any way
        if (  // bondType == Bond::DOUBLE and
            bond->getStereo() > Bond::STEREOANY &&
            bond->getStereoAtoms().size() >= 2) {
          inchi_Stereo0D stereo0D;
          if (bond->getStereo() == Bond::STEREOZ ||
              bond->getStereo() == Bond::STEREOCIS) {
            stereo0D.parity = INCHI_PARITY_ODD;
          } else {
            stereo0D.parity = INCHI_PARITY_EVEN;
          }
          stereo0D.neighbor[0] = bond->getStereoAtoms()[0];
          stereo0D.neighbor[3] = bond->getStereoAtoms()[1];
          stereo0D.neighbor[1] = atomIndex1;
          stereo0D.neighbor[2] = atomIndex2;
          if (!m->getBondBetweenAtoms(stereo0D.neighbor[0], stereo0D.neighbor[1])) {
            std::swap(stereo0D.neighbor[0], stereo0D.neighbor[3]);
          }
          stereo0D.central_atom = NO_ATOM;
          stereo0D.type = INCHI_StereoType_DoubleBond;
          stereo0DEntries.push_back(stereo0D);
        } else if (bond->getStereo() == Bond::STEREOANY) {
          // have to treat STEREOANY separately because RDKit will clear out
          // StereoAtoms information.
          // Here we just change the coordinates of the two end atoms - to bring
          // them really close - so that InChI will not try to infer stereobond
          // info from coordinates.
          inchiAtoms[atomIndex1].x = inchiAtoms[atomIndex2].x;
          inchiAtoms[atomIndex1].y = inchiAtoms[atomIndex2].y;
          inchiAtoms[atomIndex1].z = inchiAtoms[atomIndex2].z;
        }
    
        // number of bonds
        inchiAtoms[atomIndex1].num_bonds++;
      }
    
      // create stereo0D
      std::unique_ptr<inchi_Stereo0D[]> stereo0Ds;
      if (stereo0DEntries.size()) {
        stereo0Ds.reset(new inchi_Stereo0D[stereo0DEntries.size()]);
        for (unsigned int i = 0; i < stereo0DEntries.size(); i++) {
          stereo0Ds[i] = stereo0DEntries[i];
        }
      }
    
      // create input
      inchi_Input input;
      input.atom = inchiAtoms.get();
      input.stereo0D = stereo0Ds.get();
      std::unique_ptr<char[]> _options;
      if (options) {
        _options.reset(new char[strlen(options) + 1]);
        fixOptionSymbol(options, _options.get());
        input.szOptions = _options.get();
      } else {
        input.szOptions = nullptr;
      }
      input.num_atoms = nAtoms;
      input.num_stereo0D = stereo0DEntries.size();
    
      // create output
      inchi_Output output;
    
      // call DLL
      std::string inchi;
      {
        int retcode = GetINCHI(&input, &output);
    
        // generate output
        rv.returnCode = retcode;
        if (output.szInChI) {
          inchi = std::string(output.szInChI);
        }
        if (output.szMessage) {
          rv.messagePtr = std::string(output.szMessage);
        }
        if (output.szLog) {
          rv.logPtr = std::string(output.szLog);
        }
        if (output.szAuxInfo) {
          rv.auxInfoPtr = std::string(output.szAuxInfo);
        }
    
        // clean up
        FreeINCHI(&output);
      }
    
      return inchi;
    }
    */
    // END RDKIT C++ FUNCTION: MolToInchi
    // BEGIN RDKIT ACTIVE CONFIGURATION: MolToInchi
    // RDKit✔️❌: Pinned RDKit 2026.03.1; GCC/Linux; no function-local preprocessor branch.
    // RDKit✔️❌: The engine contract invokes the in-tree Rust GetINCHI/FreeINCHI port; C++ is test-only.
    // RDKit✔️❌: Adapter conformers preserve source order and only the first conformer is consumed.
    // RDKit✔️❌: Cloning and owned vectors add known allocation overhead relative to native RDKit storage.
    // END RDKIT ACTIVE CONFIGURATION: MolToInchi

    let mut diagnostics = Vec::new();
    let mut working = molecule.clone();
    if toolkit.needs_update_property_cache(molecule)? {
        toolkit.update_property_cache(&mut working, false)?;
    }
    toolkit.kekulize(&mut working, false)?;
    r_clean_up(&mut working);

    let atom_count = working.atoms.len();
    if let Some(conformer) = working.conformers.first()
        && conformer.len() != atom_count
    {
        return Err(MolToInchiError::InvalidConformer);
    }
    let mut inchi_atoms = vec![inchi_Atom::default(); atom_count];
    let mut stereo0d_entries = Vec::new();

    for index in 0..atom_count {
        let atomic_number = working.atoms[index].atomic_number;
        let isotope = working.atoms[index].isotope as i32;
        let formal_charge = working.atoms[index].formal_charge;
        let num_radical_electrons = working.atoms[index].num_radical_electrons;
        let chiral_tag = working.atoms[index].chiral_tag;
        let output = &mut inchi_atoms[index];
        output.num_bonds = 0;
        if let Some(conformer) = working.conformers.first() {
            [output.x, output.y, output.z] = conformer[index];
        }

        let element = toolkit.element_symbol(atomic_number)?;
        if element.len() >= output.elname.len() || element.contains(&0) {
            return Err(MolToInchiError::ElementSymbolTooLong);
        }
        for (target, byte) in output.elname.iter_mut().zip(element.iter().copied()) {
            *target = byte as i8;
        }
        output.elname[element.len()] = 0;

        output.isotopic_mass = if isotope != 0 {
            let reference_weight = (toolkit.atomic_weight(atomic_number)? + 0.5) as i32;
            (ISOTOPIC_SHIFT_FLAG as i32 + isotope - reference_weight) as i16
        } else {
            0
        };
        output.charge = formal_charge as i8;

        let hydrogen_count = match atomic_number {
            6 | 7 | 8 | 9 | 17 | 35 | 53 => -1,
            _ => toolkit.total_num_hydrogens(&working, index as u32)? as i32,
        };
        output.num_iso_H = [hydrogen_count as i8, 0, 0, 0];
        output.radical = 0;
        if num_radical_electrons != 0 {
            output.num_iso_H[0] =
                toolkit.total_num_hydrogens(&working, index as u32)? as i8;
        }

        if matches!(
            chiral_tag,
            ChiralTag::TetrahedralCcw | ChiralTag::TetrahedralCw
        ) {
            toolkit.calc_implicit_valence(&mut working, index as u32)?;
            let total_degree = toolkit.total_degree(&working, index as u32)?;
            if !(3..=4).contains(&total_degree) {
                diagnostics.push(AdapterDiagnostic {
                    level: AdapterDiagnosticLevel::Warning,
                    message: "tetrahedral chirality on atom with <3 or >4 neighbors will be ignored.\n"
                        .to_owned(),
                });
                continue;
            }
            let mut stereo = inchi_Stereo0D {
                central_atom: index as i16,
                type_: tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8,
                ..inchi_Stereo0D::default()
            };
            let mut neighbor_count = 0_usize;
            for &(neighbor, _) in &working.adjacency[index] {
                stereo.neighbor[neighbor_count] = neighbor as i16;
                neighbor_count += 1;
            }
            if neighbor_count == 3 {
                for position in (1..=3).rev() {
                    stereo.neighbor[position] = stereo.neighbor[position - 1];
                }
                stereo.neighbor[0] = index as i16;
            }
            let degree = working.adjacency[index].len();
            stereo.parity = if degree == 4 {
                match chiral_tag {
                    ChiralTag::TetrahedralCw => {
                        tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i8
                    }
                    ChiralTag::TetrahedralCcw => {
                        tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8
                    }
                    _ => continue,
                }
            } else {
                match chiral_tag {
                    ChiralTag::TetrahedralCcw => {
                        tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i8
                    }
                    ChiralTag::TetrahedralCw => {
                        tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8
                    }
                    other => {
                        diagnostics.push(AdapterDiagnostic {
                            level: AdapterDiagnosticLevel::Warning,
                            message: format!(
                                "unrecognized chirality tag ({}) on atom {} is ignored.\n",
                                other as u8, index
                            ),
                        });
                        continue;
                    }
                }
            };
            stereo0d_entries.push(stereo);
        }
    }

    for bond in &working.bonds {
        let mut atom_index_1 = bond.begin_atom_index;
        let mut atom_index_2 = bond.end_atom_index;
        let mut direction_modifier = 1_i32;
        if atom_index_1 > atom_index_2 {
            std::mem::swap(&mut atom_index_1, &mut atom_index_2);
            direction_modifier = -1;
        }
        let atom_output = &mut inchi_atoms[atom_index_1 as usize];
        let slot = atom_output.num_bonds as u32;
        if slot >= MAXVAL {
            diagnostics.push(AdapterDiagnostic {
                level: AdapterDiagnosticLevel::Error,
                message: format!(
                    " atom {} has too many bonds: {}. The InChI library supports at most {}\n",
                    atom_index_1, slot, MAXVAL
                ),
            });
            return Ok(MolToInchiResult {
                inchi: Vec::new(),
                diagnostics,
            });
        }
        let slot = slot as usize;
        atom_output.neighbor[slot] = atom_index_2 as i16;
        let mut bond_type = bond.bond_type;
        if bond_type as u8 > BondType::Triple as u8 {
            diagnostics.push(AdapterDiagnostic {
                level: AdapterDiagnosticLevel::Warning,
                message: format!(
                    "bond type above 3 ({}) is treated as unspecified!\n",
                    bond_type as u8
                ),
            });
            bond_type = BondType::Unspecified;
        }
        atom_output.bond_type[slot] = bond_type as i8;
        atom_output.bond_stereo[slot] = match bond.direction {
            BondDirection::BeginWedge => (direction_modifier
                * tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1UP)
                as i8,
            BondDirection::BeginDash => (direction_modifier
                * tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1DOWN)
                as i8,
            BondDirection::EitherDouble => {
                tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER as i8
            }
            BondDirection::Unknown => (direction_modifier
                * tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1EITHER)
                as i8,
            _ => tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE as i8,
        };

        if bond.stereo as u8 > BondStereo::Any as u8 && bond.stereo_atoms.len() >= 2 {
            let mut stereo = inchi_Stereo0D {
                parity: if matches!(bond.stereo, BondStereo::Z | BondStereo::Cis) {
                    tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8
                } else {
                    tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i8
                },
                neighbor: [
                    bond.stereo_atoms[0] as i16,
                    atom_index_1 as i16,
                    atom_index_2 as i16,
                    bond.stereo_atoms[1] as i16,
                ],
                central_atom: NO_ATOM as i16,
                type_: tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
            };
            if bond_index_between(
                &working,
                stereo.neighbor[0] as u32,
                stereo.neighbor[1] as u32,
            )
            .is_none()
            {
                stereo.neighbor.swap(0, 3);
            }
            stereo0d_entries.push(stereo);
        } else if bond.stereo == BondStereo::Any {
            let source = [
                inchi_atoms[atom_index_2 as usize].x,
                inchi_atoms[atom_index_2 as usize].y,
                inchi_atoms[atom_index_2 as usize].z,
            ];
            inchi_atoms[atom_index_1 as usize].x = source[0];
            inchi_atoms[atom_index_1 as usize].y = source[1];
            inchi_atoms[atom_index_1 as usize].z = source[2];
        }
        inchi_atoms[atom_index_1 as usize].num_bonds += 1;
    }

    let options_with_nul = if let Some(options) = options {
        let length = options
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(MolToInchiError::InvalidOptions)?;
        let mut converted = vec![0_u8; length + 1];
        fix_option_symbol(options, &mut converted)
            .map_err(|_| MolToInchiError::InvalidOptions)?;
        Some(converted)
    } else {
        None
    };
    let input = AdapterInchiGenerationInput {
        atoms: inchi_atoms,
        stereo0d: stereo0d_entries,
        options_with_nul,
    };
    let output = engine.get_inchi(&input)?;
    return_values.return_code = output.return_code;
    let inchi = output.inchi.unwrap_or_default();
    if let Some(message) = output.message {
        return_values.message = message;
    }
    if let Some(log) = output.log {
        return_values.log = log;
    }
    if let Some(aux_info) = output.aux_info {
        return_values.aux_info = aux_info;
    }
    engine.free_inchi()?;
    Ok(MolToInchiResult { inchi, diagnostics })
}

#[rustfmt::skip]
pub(crate) fn inchi_to_inchi_key(
    engine: &mut impl InchiKeyEngine,
    inchi: &[u8],
) -> Result<InchiToInchiKeyResult, InchiToInchiKeyError> {
    // BEGIN RDKIT C++ FUNCTION: third_party/rdkit/External/INCHI-API/inchi.cpp:2145 InchiToInchiKey
    // RDKit✔️❌: complete source frame follows verbatim.
    /*
    std::string InchiToInchiKey(const std::string &inchi) {
      char inchiKey[29];
      char xtra1[65], xtra2[65];
      int ret = 0;
      ret = GetINCHIKeyFromINCHI(inchi.c_str(), 0, 0, inchiKey, xtra1, xtra2);
      std::string error;
      switch (ret) {
        case INCHIKEY_OK:
          return std::string(inchiKey);
        case INCHIKEY_UNKNOWN_ERROR:
          error = "Unknown error";
          break;
        case INCHIKEY_EMPTY_INPUT:
          error = "Empty input";
          break;
        case INCHIKEY_INVALID_INCHI_PREFIX:
          error = "Invalid InChI prefix";
          break;
        case INCHIKEY_NOT_ENOUGH_MEMORY:
          error = "Not enough memory";
          break;
        case INCHIKEY_INVALID_INCHI:
          error = "Invalid input InChI string";
          break;
        case INCHIKEY_INVALID_STD_INCHI:
          error = "Invalid standard InChI string";
          break;
      }
      BOOST_LOG(rdErrorLog) << error << " in generating InChI Key" << std::endl;
      return std::string();
    }
    */
    // END RDKIT C++ FUNCTION: InchiToInchiKey
    // BEGIN RDKIT ACTIVE CONFIGURATION: InchiToInchiKey
    // RDKit✔️❌: Pinned RDKit 2026.03.1; GCC/Linux; no conditional source branch.
    // RDKit✔️❌: `std::string::c_str()` supplies the complete byte sequence plus one terminal
    // RDKit✔️❌: NUL; the C callee observes only the prefix before the first embedded NUL.
    // RDKit✔️❌: `GetINCHIKeyFromINCHI` is the in-tree Rust source port; production does not
    // RDKit✔️❌: call FFI or an external executable. `BOOST_LOG` is captured as a diagnostic.
    // RDKit✔️❌: The fixed 29-byte key buffer must contain the source-required terminal NUL
    // RDKit✔️❌: on success; an invalid scripted source result is a structured error.
    // END RDKIT ACTIVE CONFIGURATION: InchiToInchiKey

    let mut inchi_with_nul = Vec::with_capacity(inchi.len() + 1);
    inchi_with_nul.extend_from_slice(inchi);
    inchi_with_nul.push(0);
    let output = engine.get_inchi_key(&inchi_with_nul)?;
    if output.status == INCHIKEY_OK as i32 {
        let nul = output
            .key_buffer
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(InchiToInchiKeyError::InvalidSourceOutput(
                "successful InChIKey output is not NUL-terminated",
            ))?;
        return Ok(InchiToInchiKeyResult {
            key: output.key_buffer[..nul].to_vec(),
            diagnostics: Vec::new(),
        });
    }

    let error = match output.status {
        status if status == INCHIKEY_UNKNOWN_ERROR as i32 => b"Unknown error".as_slice(),
        status if status == INCHIKEY_EMPTY_INPUT as i32 => b"Empty input".as_slice(),
        status if status == INCHIKEY_INVALID_INCHI_PREFIX as i32 => {
            b"Invalid InChI prefix".as_slice()
        }
        status if status == INCHIKEY_NOT_ENOUGH_MEMORY as i32 => b"Not enough memory".as_slice(),
        status if status == INCHIKEY_INVALID_INCHI as i32 => {
            b"Invalid input InChI string".as_slice()
        }
        status if status == INCHIKEY_INVALID_STD_INCHI as i32 => {
            b"Invalid standard InChI string".as_slice()
        }
        _ => b"".as_slice(),
    };
    let mut message = Vec::with_capacity(error.len() + b" in generating InChI Key\n".len());
    message.extend_from_slice(error);
    message.extend_from_slice(b" in generating InChI Key\n");
    Ok(InchiToInchiKeyResult {
        key: Vec::new(),
        diagnostics: vec![AdapterDiagnostic {
            level: AdapterDiagnosticLevel::Error,
            message: String::from_utf8(message)
                .expect("RDKit InChIKey diagnostic is source-defined ASCII"),
        }],
    })
}

#[rustfmt::skip]
pub(crate) fn mol_to_inchi_key(
    generation_engine: &mut impl InchiGenerationEngine,
    key_engine: &mut impl InchiKeyEngine,
    toolkit: &mut impl MolToInchiToolkit,
    molecule: &AdapterMol,
    options: Option<&[u8]>,
) -> Result<MolToInchiKeyResult, MolToInchiKeyError> {
    // BEGIN RDKIT C++ HEADER FUNCTION: third_party/rdkit/External/INCHI-API/inchi.h:107 MolToInchiKey
    // RDKit✔️❌: complete source frame follows verbatim.
    /*
    inline std::string MolToInchiKey(const ROMol &mol, const char *options = NULL) {
      ExtraInchiReturnValues rv;
      return InchiToInchiKey(MolToInchi(mol, rv, options));
    };
    */
    // END RDKIT C++ HEADER FUNCTION: MolToInchiKey
    // BEGIN RDKIT ACTIVE CONFIGURATION: MolToInchiKey
    // RDKit✔️❌: Pinned RDKit 2026.03.1; GCC/Linux; header-defined inline behavior with no
    // RDKit✔️❌: conditional source branch. C++ fully evaluates `MolToInchi` before entering
    // RDKit✔️❌: `InchiToInchiKey`; an exception from either call immediately leaves the wrapper.
    // RDKit✔️❌: The local `ExtraInchiReturnValues` is not observable by the caller. Production
    // RDKit✔️❌: invokes only the two in-tree Rust source ports and performs no format transit.
    // END RDKIT ACTIVE CONFIGURATION: MolToInchiKey

    let mut return_values = ExtraInchiReturnValues::default();
    let generated = mol_to_inchi(
        generation_engine,
        toolkit,
        molecule,
        &mut return_values,
        options,
    )?;
    let keyed = inchi_to_inchi_key(key_engine, &generated.inchi)?;
    let mut diagnostics = generated.diagnostics;
    diagnostics.extend(keyed.diagnostics);
    Ok(MolToInchiKeyResult {
        key: keyed.key,
        diagnostics,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone)]
    struct ScriptedInchiKeyEngine {
        output: Result<AdapterInchiKeyOutput, SourceHeapError>,
        seen_inputs: Vec<Vec<u8>>,
    }

    impl ScriptedInchiKeyEngine {
        fn with_status(status: i32) -> Self {
            Self {
                output: Ok(AdapterInchiKeyOutput {
                    status,
                    key_buffer: Vec::new(),
                }),
                seen_inputs: Vec::new(),
            }
        }
    }

    impl InchiKeyEngine for ScriptedInchiKeyEngine {
        fn get_inchi_key(&mut self, inchi_with_nul: &[u8]) -> Result<AdapterInchiKeyOutput, SourceHeapError> {
            self.seen_inputs.push(inchi_with_nul.to_vec());
            self.output.clone()
        }
    }

    #[test]
    fn source_port__inchi__inchitoinchikey__line_2145() {
        let mut success_buffer = vec![0xA5; 29];
        let expected_key = b"VNWKTOKETHGBQD-UHFFFAOYSA-N";
        success_buffer[..expected_key.len()].copy_from_slice(expected_key);
        success_buffer[expected_key.len()] = 0;
        let mut success = ScriptedInchiKeyEngine {
            output: Ok(AdapterInchiKeyOutput {
                status: INCHIKEY_OK as i32,
                key_buffer: success_buffer,
            }),
            seen_inputs: Vec::new(),
        };
        assert_eq!(
            inchi_to_inchi_key(&mut success, b"InChI=1S/CH4/h1H4").unwrap(),
            InchiToInchiKeyResult {
                key: expected_key.to_vec(),
                diagnostics: Vec::new(),
            }
        );
        assert_eq!(success.seen_inputs, [b"InChI=1S/CH4/h1H4\0".to_vec()]);

        let failures: &[(i32, &str)] = &[
            (INCHIKEY_UNKNOWN_ERROR as i32, "Unknown error"),
            (INCHIKEY_EMPTY_INPUT as i32, "Empty input"),
            (INCHIKEY_INVALID_INCHI_PREFIX as i32, "Invalid InChI prefix"),
            (INCHIKEY_NOT_ENOUGH_MEMORY as i32, "Not enough memory"),
            (INCHIKEY_INVALID_INCHI as i32, "Invalid input InChI string"),
            (INCHIKEY_INVALID_STD_INCHI as i32, "Invalid standard InChI string"),
            (i32::MAX, ""),
        ];
        for &(status, prefix) in failures {
            let mut engine = ScriptedInchiKeyEngine::with_status(status);
            assert_eq!(
                inchi_to_inchi_key(&mut engine, b"").unwrap(),
                InchiToInchiKeyResult {
                    key: Vec::new(),
                    diagnostics: vec![AdapterDiagnostic {
                        level: AdapterDiagnosticLevel::Error,
                        message: format!("{prefix} in generating InChI Key\n"),
                    }],
                },
                "status {status}"
            );
            assert_eq!(engine.seen_inputs, [vec![0]], "status {status}");
        }

        let embedded = b"InChI=1S/CH4/h1H4\0ignored";
        let mut embedded_engine = ScriptedInchiKeyEngine::with_status(INCHIKEY_EMPTY_INPUT as i32);
        let _ = inchi_to_inchi_key(&mut embedded_engine, embedded).unwrap();
        let mut expected_input = embedded.to_vec();
        expected_input.push(0);
        assert_eq!(embedded_engine.seen_inputs, [expected_input]);

        let mut engine_error = ScriptedInchiKeyEngine {
            output: Err(SourceHeapError::AllocationFailed),
            seen_inputs: Vec::new(),
        };
        assert_eq!(
            inchi_to_inchi_key(&mut engine_error, b"InChI=1S/CH4/h1H4"),
            Err(InchiToInchiKeyError::Source(SourceHeapError::AllocationFailed))
        );

        let mut invalid_success = ScriptedInchiKeyEngine {
            output: Ok(AdapterInchiKeyOutput {
                status: INCHIKEY_OK as i32,
                key_buffer: vec![b'X'; 29],
            }),
            seen_inputs: Vec::new(),
        };
        assert_eq!(
            inchi_to_inchi_key(&mut invalid_success, b"InChI=1S/CH4/h1H4"),
            Err(InchiToInchiKeyError::InvalidSourceOutput(
                "successful InChIKey output is not NUL-terminated"
            ))
        );

        let mut heap = SourceHeap::default();
        let result = {
            let mut source_engine = SourceInchiKeyEngine { heap: &mut heap };
            inchi_to_inchi_key(&mut source_engine, b"InChI=1S/CH4/h1H4").unwrap()
        };
        assert_eq!(result.key, expected_key);
        assert!(result.diagnostics.is_empty());
        assert_eq!(heap.live_allocation_count(), 0);
        assert_eq!(heap.live_source_allocation_count(), 0);
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__inchitoinchikey__exact() {
        use serde_json::{Value, json};
        use std::process::Command;

        fn bytes(value: &Value, field: &str, case_id: &str) -> Vec<u8> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: {field} must be an array"))
                .iter()
                .map(|byte| {
                    u8::try_from(
                        byte.as_u64()
                            .unwrap_or_else(|| panic!("{case_id}: {field} byte must be unsigned")),
                    )
                    .unwrap_or_else(|_| panic!("{case_id}: {field} byte exceeds u8"))
                })
                .collect()
        }

        let repository_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(std::path::Path::parent)
            .expect("crate must be below repository root");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(runner)
            .arg("--inchi-to-inchi-key-records")
            .current_dir(repository_root)
            .output()
            .expect("pinned RDKit C++ oracle must run");
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let records = output.lines().collect::<Vec<_>>();
        assert_eq!(records.len(), 9, "InchiToInchiKey oracle case count changed");
        let mut case_ids = BTreeSet::new();

        for line in records {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            let case_id = official["case_id"].as_str().expect("case_id must be a string");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");
            let input = bytes(&official["input"]["inchi"], "input.inchi", case_id);
            let status = i32::try_from(
                official["input"]["scripted_status"]
                    .as_i64()
                    .unwrap_or_else(|| panic!("{case_id}: scripted_status must be signed")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: scripted_status exceeds i32"));
            let key_buffer = bytes(
                &official["input"]["scripted_key_buffer"],
                "input.scripted_key_buffer",
                case_id,
            );
            let mut engine = ScriptedInchiKeyEngine {
                output: Ok(AdapterInchiKeyOutput {
                    status,
                    key_buffer: key_buffer.clone(),
                }),
                seen_inputs: Vec::new(),
            };
            let result = inchi_to_inchi_key(&mut engine, &input)
                .unwrap_or_else(|error| panic!("{case_id}: Rust wrapper failed: {error:?}"));
            assert_eq!(engine.seen_inputs.len(), 1, "{case_id}: engine calls");
            let seen = &engine.seen_inputs[0];
            let captured_end = seen
                .iter()
                .position(|byte| *byte == 0)
                .expect("wrapper must append a terminal NUL");
            let captured_inchi = seen[..captured_end].to_vec();
            let error_text = result
                .diagnostics
                .iter()
                .map(|diagnostic| {
                    assert_eq!(
                        diagnostic.level,
                        AdapterDiagnosticLevel::Error,
                        "{case_id}: diagnostic stream"
                    );
                    diagnostic.message.as_str()
                })
                .collect::<String>();
            let actual = json!({
                "schema_version": "cosmolkit-inchi-rdkit-cpp-v1",
                "rdkit_version": "2026.03.1",
                "source_sha256": "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f",
                "source_fragment_sha256": "b2d19205873e8444cef0ce660bbcf36bbf35cc528e77802b40b9b28ec15e67ae",
                "operation": "InchiToInchiKey",
                "case_id": case_id,
                "input": {
                    "inchi": input,
                    "scripted_status": status,
                    "scripted_key_buffer": key_buffer,
                },
                "output": {
                    "key": result.key,
                    "error_text": error_text,
                    "captured_inchi": captured_inchi,
                    "captured_xtra1": 0,
                    "captured_xtra2": 0,
                    "call_count": engine.seen_inputs.len(),
                },
            });
            assert_eq!(actual, official, "{case_id}: exact record mismatch");
        }

        assert_eq!(
            case_ids,
            BTreeSet::from([
                "empty-input".to_owned(),
                "invalid-inchi-embedded-nul".to_owned(),
                "invalid-prefix".to_owned(),
                "invalid-standard-inchi".to_owned(),
                "not-enough-memory".to_owned(),
                "success-first-nul".to_owned(),
                "success-full-key".to_owned(),
                "unknown-error".to_owned(),
                "unknown-status".to_owned(),
            ])
        );
    }

    #[test]
    fn source_port__inchi__moltoinchikey__line_107() {
        let molecule = AdapterMol::from_graph(
            vec![
                AdapterAtom {
                    atomic_number: 6,
                    ..AdapterAtom::default()
                },
                AdapterAtom {
                    atomic_number: 8,
                    ..AdapterAtom::default()
                },
            ],
            vec![AdapterBond {
                begin_atom_index: 0,
                end_atom_index: 1,
                bond_type: BondType::Quadruple,
                ..AdapterBond::default()
            }],
        );
        let molecule_before = molecule.clone();
        let mut generation = ScriptedGenerationEngine::default();
        generation.output.inchi = Some(b"InChI=1S/generated".to_vec());
        let mut key = ScriptedInchiKeyEngine::with_status(INCHIKEY_INVALID_INCHI as i32);
        let mut toolkit = RecordingGenerationToolkit::default();
        let result = mol_to_inchi_key(
            &mut generation,
            &mut key,
            &mut toolkit,
            &molecule,
            Some(b"/FixedH\0ignored"),
        )
        .unwrap();
        assert_eq!(
            result,
            MolToInchiKeyResult {
                key: Vec::new(),
                diagnostics: vec![
                    AdapterDiagnostic {
                        level: AdapterDiagnosticLevel::Warning,
                        message: "bond type above 3 (4) is treated as unspecified!\n".to_owned(),
                    },
                    AdapterDiagnostic {
                        level: AdapterDiagnosticLevel::Error,
                        message: "Invalid input InChI string in generating InChI Key\n".to_owned(),
                    },
                ],
            }
        );
        assert_eq!(generation.calls, ["GetINCHI", "FreeINCHI"]);
        assert_eq!(generation.seen_inputs.len(), 1);
        assert_eq!(generation.seen_inputs[0].options_with_nul, Some(b"-FixedH\0".to_vec()));
        assert_eq!(key.seen_inputs, [b"InChI=1S/generated\0".to_vec()]);
        assert_eq!(molecule, molecule_before);

        let mut success_generation = ScriptedGenerationEngine::default();
        let mut success_key = ScriptedInchiKeyEngine {
            output: Ok(AdapterInchiKeyOutput {
                status: INCHIKEY_OK as i32,
                key_buffer: b"KEY\0ignored".to_vec(),
            }),
            seen_inputs: Vec::new(),
        };
        let mut success_toolkit = RecordingGenerationToolkit::default();
        let success = mol_to_inchi_key(
            &mut success_generation,
            &mut success_key,
            &mut success_toolkit,
            &AdapterMol::from_graph(Vec::new(), Vec::new()),
            None,
        )
        .unwrap();
        assert_eq!(success.key, b"KEY");
        assert!(success.diagnostics.is_empty());

        let mut failed_generation = ScriptedGenerationEngine {
            get_error: Some(SourceHeapError::AllocationFailed),
            ..ScriptedGenerationEngine::default()
        };
        let mut uncalled_key = ScriptedInchiKeyEngine::with_status(INCHIKEY_OK as i32);
        let mut failed_toolkit = RecordingGenerationToolkit::default();
        assert_eq!(
            mol_to_inchi_key(
                &mut failed_generation,
                &mut uncalled_key,
                &mut failed_toolkit,
                &AdapterMol::from_graph(Vec::new(), Vec::new()),
                None,
            ),
            Err(MolToInchiKeyError::MolToInchi(MolToInchiError::Source(
                SourceHeapError::AllocationFailed
            )))
        );
        assert!(uncalled_key.seen_inputs.is_empty());

        let mut generated_before_key_error = ScriptedGenerationEngine::default();
        let mut failed_key = ScriptedInchiKeyEngine {
            output: Err(SourceHeapError::PointerOutOfBounds),
            seen_inputs: Vec::new(),
        };
        let mut key_error_toolkit = RecordingGenerationToolkit::default();
        assert_eq!(
            mol_to_inchi_key(
                &mut generated_before_key_error,
                &mut failed_key,
                &mut key_error_toolkit,
                &AdapterMol::from_graph(Vec::new(), Vec::new()),
                None,
            ),
            Err(MolToInchiKeyError::InchiToInchiKey(InchiToInchiKeyError::Source(
                SourceHeapError::PointerOutOfBounds
            )))
        );
        assert_eq!(generated_before_key_error.calls, ["GetINCHI", "FreeINCHI"]);
        assert_eq!(failed_key.seen_inputs.len(), 1);
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__moltoinchikey__exact() {
        use serde_json::{Value, json};
        use std::cell::RefCell;
        use std::process::Command;
        use std::rc::Rc;

        struct OrderedGenerationEngine {
            inner: ScriptedGenerationEngine,
            calls: Rc<RefCell<Vec<&'static str>>>,
        }

        impl InchiGenerationEngine for OrderedGenerationEngine {
            fn get_inchi(
                &mut self,
                input: &AdapterInchiGenerationInput,
            ) -> Result<AdapterInchiGenerationOutput, SourceHeapError> {
                self.calls.borrow_mut().push("GetINCHI");
                self.inner.get_inchi(input)
            }

            fn free_inchi(&mut self) -> Result<(), SourceHeapError> {
                self.calls.borrow_mut().push("FreeINCHI");
                self.inner.free_inchi()
            }
        }

        struct OrderedKeyEngine {
            inner: ScriptedInchiKeyEngine,
            calls: Rc<RefCell<Vec<&'static str>>>,
        }

        impl InchiKeyEngine for OrderedKeyEngine {
            fn get_inchi_key(&mut self, inchi_with_nul: &[u8]) -> Result<AdapterInchiKeyOutput, SourceHeapError> {
                self.calls.borrow_mut().push("GetINCHIKeyFromINCHI");
                self.inner.get_inchi_key(inchi_with_nul)
            }
        }

        fn optional_bytes(value: &Value, case_id: &str) -> Option<Vec<u8>> {
            (!value.is_null()).then(|| mol_to_inchi_oracle_bytes(value, case_id))
        }

        fn original_atoms(molecule: &AdapterMol) -> Vec<Value> {
            molecule
                .atoms
                .iter()
                .enumerate()
                .map(|(index, atom)| {
                    json!({
                        "index": index,
                        "atomic_number": atom.atomic_number,
                        "formal_charge": atom.formal_charge,
                        "num_explicit_hydrogens": atom.num_explicit_hydrogens,
                        "is_aromatic": atom.is_aromatic,
                        "isotope": atom.isotope,
                        "num_radical_electrons": atom.num_radical_electrons,
                        "no_implicit": atom.no_implicit,
                        "chiral_tag": atom.chiral_tag as u8,
                        "cip_rank": atom.cip_rank,
                    })
                })
                .collect()
        }

        fn original_bonds(molecule: &AdapterMol) -> Vec<Value> {
            molecule
                .bonds
                .iter()
                .enumerate()
                .map(|(index, bond)| {
                    json!({
                        "index": index,
                        "begin_atom_index": bond.begin_atom_index,
                        "end_atom_index": bond.end_atom_index,
                        "bond_type": bond.bond_type as u8,
                        "direction": bond.direction as u8,
                        "is_aromatic": bond.is_aromatic,
                        "stereo": bond.stereo as u8,
                        "stereo_atoms": bond.stereo_atoms,
                    })
                })
                .collect()
        }

        let repository_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(std::path::Path::parent)
            .expect("crate must be below repository root");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(runner)
            .arg("--mol-to-inchi-key-records")
            .current_dir(repository_root)
            .output()
            .expect("pinned RDKit C++ oracle must run");
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let records = output.lines().collect::<Vec<_>>();
        assert_eq!(records.len(), 9, "MolToInchiKey oracle case count changed");
        let mut case_ids = BTreeSet::new();

        for line in records {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            let case_id = official["case_id"].as_str().expect("case_id must be a string");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "27b4eaef1714869c42dbf2998807018a03389b4f9ce40438e843248ebfc3614e"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "f6cc64a50170bb62e6c7a97325815f1399a6d3c74a9bc1b4a8bf69da90ccae51"
            );
            assert_eq!(official["operation"], "MolToInchiKey");
            assert_eq!(
                official["dependency_fragment_sha256"],
                json!({
                    "MolToInchi": "000d4eff4353e06aaa802057674fa24f91726872b0c0f9bb593d7a1937721ac2",
                    "InchiToInchiKey": "b2d19205873e8444cef0ce660bbcf36bbf35cc528e77802b40b9b28ec15e67ae",
                })
            );

            let input = &official["input"];
            assert_eq!(input.as_object().unwrap().len(), 10, "{case_id}: input fields");
            let atom_fields = input["atoms"].as_array().expect("atoms must be an array");
            let mut hydrogens = Vec::with_capacity(atom_fields.len());
            let mut total_degrees = Vec::with_capacity(atom_fields.len());
            let atoms = atom_fields
                .iter()
                .map(|field| {
                    assert_eq!(field.as_object().unwrap().len(), 10, "{case_id}: atom fields");
                    hydrogens.push(mol_to_inchi_oracle_u32(
                        &field["total_hydrogens"],
                        "total_hydrogens",
                        case_id,
                    ));
                    total_degrees.push(if field["total_degree"].is_null() {
                        None
                    } else {
                        Some(mol_to_inchi_oracle_u32(&field["total_degree"], "total_degree", case_id))
                    });
                    AdapterAtom {
                        atomic_number: mol_to_inchi_oracle_i32(&field["atomic_number"], "atomic_number", case_id),
                        formal_charge: mol_to_inchi_oracle_i32(&field["formal_charge"], "formal_charge", case_id),
                        num_explicit_hydrogens: mol_to_inchi_oracle_u32(
                            &field["explicit_hydrogens"],
                            "explicit_hydrogens",
                            case_id,
                        ),
                        is_aromatic: field["aromatic"].as_bool().expect("aromatic must be bool"),
                        isotope: mol_to_inchi_oracle_u32(&field["isotope"], "isotope", case_id),
                        num_radical_electrons: mol_to_inchi_oracle_u32(
                            &field["radical_electrons"],
                            "radical_electrons",
                            case_id,
                        ),
                        no_implicit: field["no_implicit"].as_bool().expect("no_implicit must be bool"),
                        chiral_tag: mol_to_inchi_oracle_chiral_tag(&field["chiral_tag"], case_id),
                        cip_rank: None,
                        cached_explicit_valence: None,
                        cached_implicit_valence: None,
                    }
                })
                .collect::<Vec<_>>();
            let bonds = input["bonds"]
                .as_array()
                .expect("bonds must be an array")
                .iter()
                .map(|field| AdapterBond {
                    begin_atom_index: mol_to_inchi_oracle_u32(&field["begin_atom_index"], "begin_atom_index", case_id),
                    end_atom_index: mol_to_inchi_oracle_u32(&field["end_atom_index"], "end_atom_index", case_id),
                    bond_type: mol_to_inchi_oracle_bond_type(&field["bond_type"], case_id),
                    direction: mol_to_inchi_oracle_direction(&field["direction"], case_id),
                    is_aromatic: field["aromatic"].as_bool().expect("aromatic must be bool"),
                    stereo: mol_to_inchi_oracle_stereo(&field["stereo"], case_id),
                    stereo_atoms: field["stereo_atoms"]
                        .as_array()
                        .expect("stereo_atoms must be an array")
                        .iter()
                        .map(|value| mol_to_inchi_oracle_u32(value, "stereo_atom", case_id))
                        .collect(),
                })
                .collect::<Vec<_>>();
            let mut molecule = AdapterMol::from_graph(atoms, bonds);
            molecule.conformers = input["conformers"]
                .as_array()
                .expect("conformers must be an array")
                .iter()
                .map(|conformer| {
                    conformer
                        .as_array()
                        .expect("conformer must be an array")
                        .iter()
                        .map(|point| {
                            let point = point.as_array().expect("coordinate must be an array");
                            assert_eq!(point.len(), 3, "{case_id}: coordinate width");
                            [
                                point[0].as_f64().expect("x must be numeric"),
                                point[1].as_f64().expect("y must be numeric"),
                                point[2].as_f64().expect("z must be numeric"),
                            ]
                        })
                        .collect()
                })
                .collect();
            let molecule_before = molecule.clone();

            let scripted = &input["scripted_generation_output"];
            let throw_on_get = scripted["throw_on_get"].as_bool().expect("throw_on_get must be bool");
            let throw_on_free = scripted["throw_on_free"].as_bool().expect("throw_on_free must be bool");
            let shared_calls = Rc::new(RefCell::new(Vec::new()));
            let mut generation = OrderedGenerationEngine {
                inner: ScriptedGenerationEngine {
                    output: AdapterInchiGenerationOutput {
                        return_code: mol_to_inchi_oracle_i32(&scripted["return_code"], "return_code", case_id),
                        inchi: optional_bytes(&scripted["inchi"], case_id),
                        message: optional_bytes(&scripted["message"], case_id),
                        log: optional_bytes(&scripted["log"], case_id),
                        aux_info: optional_bytes(&scripted["aux_info"], case_id),
                    },
                    get_error: throw_on_get.then_some(SourceHeapError::AllocationFailed),
                    free_error: throw_on_free.then_some(SourceHeapError::MissingAllocation),
                    ..ScriptedGenerationEngine::default()
                },
                calls: Rc::clone(&shared_calls),
            };
            let throw_on_key = input["throw_on_key"].as_bool().expect("throw_on_key must be bool");
            let mut key = OrderedKeyEngine {
                inner: ScriptedInchiKeyEngine {
                    output: if throw_on_key {
                        Err(SourceHeapError::PointerOutOfBounds)
                    } else {
                        Ok(AdapterInchiKeyOutput {
                            status: mol_to_inchi_oracle_i32(
                                &input["scripted_key_status"],
                                "scripted_key_status",
                                case_id,
                            ),
                            key_buffer: mol_to_inchi_oracle_bytes(&input["scripted_key_buffer"], case_id),
                        })
                    },
                    seen_inputs: Vec::new(),
                },
                calls: Rc::clone(&shared_calls),
            };
            let mut toolkit = RecordingGenerationToolkit {
                needs_update: input["needs_update_property_cache"]
                    .as_bool()
                    .expect("needs_update_property_cache must be bool"),
                fail_on: input["generation_fail_on"].as_str().map(str::to_owned),
                hydrogens,
                total_degrees,
                ..RecordingGenerationToolkit::default()
            };
            let options = if input["options"].is_null() {
                None
            } else {
                Some(mol_to_inchi_oracle_bytes(&input["options"], case_id))
            };
            let result = mol_to_inchi_key(&mut generation, &mut key, &mut toolkit, &molecule, options.as_deref());

            let (status, exception_kind, exception_message, key_bytes, diagnostics) = match result {
                Ok(result) => ("returned", Value::Null, Value::Null, result.key, result.diagnostics),
                Err(MolToInchiKeyError::MolToInchi(MolToInchiError::Toolkit(error))) => (
                    "exception",
                    json!("MolSanitizeException"),
                    json!(error.message),
                    Vec::new(),
                    Vec::new(),
                ),
                Err(MolToInchiKeyError::MolToInchi(MolToInchiError::Source(_))) if throw_on_get => (
                    "exception",
                    json!("runtime_error"),
                    json!("GetINCHI"),
                    Vec::new(),
                    Vec::new(),
                ),
                Err(MolToInchiKeyError::MolToInchi(MolToInchiError::Source(_))) if throw_on_free => (
                    "exception",
                    json!("runtime_error"),
                    json!("FreeINCHI"),
                    Vec::new(),
                    Vec::new(),
                ),
                Err(MolToInchiKeyError::InchiToInchiKey(InchiToInchiKeyError::Source(_))) if throw_on_key => (
                    "exception",
                    json!("runtime_error"),
                    json!("GetINCHIKeyFromINCHI"),
                    Vec::new(),
                    Vec::new(),
                ),
                Err(error) => panic!("{case_id}: unexpected Rust error {error:?}"),
            };
            let warning_text = diagnostics
                .iter()
                .filter(|diagnostic| diagnostic.level == AdapterDiagnosticLevel::Warning)
                .map(|diagnostic| diagnostic.message.as_str())
                .collect::<String>();
            let error_text = diagnostics
                .iter()
                .filter(|diagnostic| diagnostic.level == AdapterDiagnosticLevel::Error)
                .map(|diagnostic| diagnostic.message.as_str())
                .collect::<String>();
            let captured = generation.inner.seen_inputs.first();
            let captured_atoms = captured
                .map(|input| {
                    input
                        .atoms
                        .iter()
                        .map(|atom| {
                            let element = atom
                                .elname
                                .iter()
                                .take_while(|byte| **byte != 0)
                                .map(|byte| *byte as u8)
                                .collect::<Vec<_>>();
                            let count = usize::try_from(atom.num_bonds).expect("num_bonds must be nonnegative");
                            json!({
                                "x": atom.x,
                                "y": atom.y,
                                "z": atom.z,
                                "element": String::from_utf8(element).expect("element must be UTF-8"),
                                "isotopic_mass": atom.isotopic_mass,
                                "charge": atom.charge,
                                "radical": atom.radical,
                                "num_iso_h": atom.num_iso_H,
                                "num_bonds": atom.num_bonds,
                                "neighbors": atom.neighbor[..count],
                                "bond_types": atom.bond_type[..count],
                                "bond_stereo": atom.bond_stereo[..count],
                            })
                        })
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            let captured_stereo0d = captured
                .map(|input| {
                    input
                        .stereo0d
                        .iter()
                        .map(|stereo| {
                            json!({
                                "neighbors": stereo.neighbor,
                                "central_atom": stereo.central_atom,
                                "type": stereo.type_,
                                "parity": stereo.parity,
                            })
                        })
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            let captured_options = captured
                .and_then(|input| input.options_with_nul.as_ref())
                .map(|options| {
                    let nul = options
                        .iter()
                        .position(|byte| *byte == 0)
                        .expect("captured options must be NUL terminated");
                    options[..nul].to_vec()
                });
            let captured_key_inchi = key
                .inner
                .seen_inputs
                .first()
                .map(|input| {
                    let nul = input
                        .iter()
                        .position(|byte| *byte == 0)
                        .expect("captured InChI must be NUL terminated");
                    input[..nul].to_vec()
                })
                .unwrap_or_default();
            let key_call_count = key.inner.seen_inputs.len();
            let captured_key_xtra = if key_call_count == 0 { -1 } else { 0 };
            let actual_output = json!({
                "status": status,
                "exception_kind": exception_kind,
                "exception_message": exception_message,
                "key": key_bytes,
                "warning_text": warning_text,
                "error_text": error_text,
                "toolkit_calls": toolkit.calls,
                "adapter_calls": shared_calls.borrow().clone(),
                "captured_atoms": captured_atoms,
                "captured_stereo0d": captured_stereo0d,
                "captured_options": captured_options,
                "captured_key_inchi": captured_key_inchi,
                "captured_key_xtra1": captured_key_xtra,
                "captured_key_xtra2": captured_key_xtra,
                "generation_get_count": generation.inner.calls.iter().filter(|call| **call == "GetINCHI").count(),
                "generation_free_count": generation.inner.calls.iter().filter(|call| **call == "FreeINCHI").count(),
                "generation_outstanding_outputs": 0,
                "key_call_count": key_call_count,
                "original_atoms": original_atoms(&molecule),
                "original_bonds": original_bonds(&molecule),
                "original_conformers": molecule.conformers,
            });
            assert_eq!(molecule, molecule_before, "{case_id}: original molecule mutated");
            assert_eq!(
                actual_output, official["output"],
                "{case_id}: exact MolToInchiKey record mismatch"
            );
        }

        assert_eq!(
            case_ids,
            BTreeSet::from([
                "empty-generation-output".to_owned(),
                "free-exception-short-circuit".to_owned(),
                "generation-return-code-does-not-short-circuit".to_owned(),
                "generation-warning-before-key-error".to_owned(),
                "get-exception-short-circuit".to_owned(),
                "key-exception-after-cleanup".to_owned(),
                "maxval-empty-return-still-keys".to_owned(),
                "success-complete-forwarding".to_owned(),
                "toolkit-exception-short-circuit".to_owned(),
            ])
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_inchi_core_scalar_exact() {
        use std::sync::{Arc, Barrier};

        let tests: [fn(); 4] = [
            rdkit_cpp_oracle__moltoinchi__exact,
            rdkit_cpp_oracle__moltoinchikey__exact,
            rdkit_cpp_oracle__inchitoinchikey__exact,
            rdkit_cpp_oracle__inchitomol__exact,
        ];
        let barrier = Arc::new(Barrier::new(tests.len()));
        std::thread::scope(|scope| {
            let handles = tests
                .into_iter()
                .map(|test| {
                    let barrier = Arc::clone(&barrier);
                    scope.spawn(move || {
                        barrier.wait();
                        test();
                    })
                })
                .collect::<Vec<_>>();
            for handle in handles {
                handle.join().expect("pinned RDKit scalar oracle worker must not panic");
            }
        });
    }

    fn molecule(directions: &[BondDirection]) -> AdapterMol {
        AdapterMol {
            atoms: Vec::new(),
            bonds: directions
                .iter()
                .copied()
                .map(|direction| AdapterBond {
                    begin_atom_index: 0,
                    end_atom_index: 0,
                    bond_type: BondType::Unspecified,
                    direction,
                    ..AdapterBond::default()
                })
                .collect(),
            conformers: Vec::new(),
            adjacency: Vec::new(),
        }
    }

    fn graph(atoms: &[(i32, i32)], bonds: &[(u32, u32, BondType)]) -> AdapterMol {
        AdapterMol::from_graph(
            atoms
                .iter()
                .map(|&(atomic_number, formal_charge)| AdapterAtom {
                    atomic_number,
                    formal_charge,
                    num_explicit_hydrogens: 0,
                    is_aromatic: false,
                    ..AdapterAtom::default()
                })
                .collect(),
            bonds
                .iter()
                .map(|&(begin_atom_index, end_atom_index, bond_type)| AdapterBond {
                    begin_atom_index,
                    end_atom_index,
                    bond_type,
                    direction: BondDirection::None,
                    ..AdapterBond::default()
                })
                .collect(),
        )
    }

    #[test]
    fn cleanup_explicit_valence_records_source_property_cache_write() {
        let mut mol = graph(&[(7, 0), (8, 0)], &[(0, 1, BondType::Double)]);
        mol.atoms[0].num_explicit_hydrogens = 3;

        assert_eq!(mol.atoms[0].cached_explicit_valence, None);
        assert_eq!(cleanup_explicit_valence(&mut mol, 0), Ok(5));
        assert_eq!(mol.atoms[0].cached_explicit_valence, Some(5));
    }

    #[derive(Clone)]
    struct ScriptedStructureEngine {
        output: AdapterInchiStructureOutput,
        seen_inputs: Vec<Vec<u8>>,
        free_count: usize,
        get_error: Option<SourceHeapError>,
        free_error: Option<SourceHeapError>,
    }

    impl ScriptedStructureEngine {
        fn new(output: AdapterInchiStructureOutput) -> Self {
            Self {
                output,
                seen_inputs: Vec::new(),
                free_count: 0,
                get_error: None,
                free_error: None,
            }
        }
    }

    impl InchiStructureEngine for ScriptedStructureEngine {
        fn get_struct_from_inchi(
            &mut self,
            inchi_with_nul: &[u8],
        ) -> Result<AdapterInchiStructureOutput, SourceHeapError> {
            self.seen_inputs.push(inchi_with_nul.to_vec());
            if let Some(error) = self.get_error {
                return Err(error);
            }
            Ok(self.output.clone())
        }

        fn free_struct_from_inchi(&mut self) -> Result<(), SourceHeapError> {
            self.free_count += 1;
            if let Some(error) = self.free_error {
                return Err(error);
            }
            Ok(())
        }
    }

    #[derive(Default)]
    struct RecordingToolkit {
        calls: Vec<&'static str>,
        ranks: Option<Vec<u32>>,
        fail_on: Option<&'static str>,
    }

    impl RecordingToolkit {
        fn record(&mut self, call: &'static str) -> Result<(), AdapterToolkitError> {
            self.calls.push(call);
            if self.fail_on == Some(call) {
                return Err(AdapterToolkitError {
                    kind: "scripted toolkit error",
                    message: call.to_owned(),
                });
            }
            Ok(())
        }
    }

    impl InchiToMolToolkit for RecordingToolkit {
        fn atomic_number(&mut self, element: &[u8]) -> Result<i32, AdapterToolkitError> {
            self.record("atomic_number")?;
            match element {
                b"H" => Ok(1),
                b"C" => Ok(6),
                b"N" => Ok(7),
                b"O" => Ok(8),
                b"S" => Ok(16),
                b"Cl" => Ok(17),
                b"Br" => Ok(35),
                _ => Err(AdapterToolkitError {
                    kind: "Post-condition Violation",
                    message: format!("Element '{}' not found", String::from_utf8_lossy(element)),
                }),
            }
        }

        fn average_atomic_weight(&mut self, atomic_number: i32) -> Result<f64, AdapterToolkitError> {
            self.record("average_atomic_weight")?;
            match atomic_number {
                1 => Ok(1.008),
                6 => Ok(12.011),
                7 => Ok(14.007),
                8 => Ok(15.999),
                16 => Ok(32.06),
                17 => Ok(35.45),
                35 => Ok(79.904),
                _ => Err(AdapterToolkitError {
                    kind: "Pre-condition Violation",
                    message: "Atomic number not found".to_owned(),
                }),
            }
        }

        fn update_property_cache(
            &mut self,
            _molecule: &mut AdapterMol,
            strict: bool,
        ) -> Result<(), AdapterToolkitError> {
            assert!(!strict);
            self.record("update_property_cache")
        }

        fn assign_atom_cip_ranks(&mut self, molecule: &mut AdapterMol) -> Result<Vec<u32>, AdapterToolkitError> {
            self.record("assign_atom_cip_ranks")?;
            Ok(self
                .ranks
                .clone()
                .unwrap_or_else(|| (0..molecule.atoms.len() as u32).collect()))
        }

        fn remove_hydrogens(&mut self, _molecule: &mut AdapterMol) -> Result<(), AdapterToolkitError> {
            self.record("remove_hydrogens")
        }

        fn sanitize_molecule(&mut self, _molecule: &mut AdapterMol) -> Result<(), AdapterToolkitError> {
            self.record("sanitize_molecule")
        }

        fn assign_stereochemistry(
            &mut self,
            _molecule: &mut AdapterMol,
            clean_it: bool,
            force: bool,
        ) -> Result<(), AdapterToolkitError> {
            assert!(clean_it);
            assert!(force);
            self.record("assign_stereochemistry")
        }
    }

    fn source_atom(element: &[u8]) -> inchi_Atom {
        let mut atom = inchi_Atom::default();
        for (destination, source) in atom.elname.iter_mut().zip(element.iter().copied()) {
            *destination = source as i8;
        }
        atom.elname[element.len()] = 0;
        atom
    }

    fn connect_source_atoms(
        atoms: &mut [inchi_Atom],
        first: usize,
        second: usize,
        bond_type: i8,
        first_stereo: i8,
        second_stereo: i8,
    ) {
        let first_offset = atoms[first].num_bonds as usize;
        atoms[first].neighbor[first_offset] = second as i16;
        atoms[first].bond_type[first_offset] = bond_type;
        atoms[first].bond_stereo[first_offset] = first_stereo;
        atoms[first].num_bonds += 1;

        let second_offset = atoms[second].num_bonds as usize;
        atoms[second].neighbor[second_offset] = first as i16;
        atoms[second].bond_type[second_offset] = bond_type;
        atoms[second].bond_stereo[second_offset] = second_stereo;
        atoms[second].num_bonds += 1;
    }

    fn scripted_output(atoms: Vec<inchi_Atom>, stereo0d: Vec<inchi_Stereo0D>) -> AdapterInchiStructureOutput {
        AdapterInchiStructureOutput {
            return_code: tagRetValGetINCHI_inchi_Ret_OKAY,
            atoms,
            stereo0d,
            message: Some(b"message".to_vec()),
            log: Some(b"log".to_vec()),
        }
    }

    fn run_scripted(
        output: AdapterInchiStructureOutput,
        toolkit: &mut RecordingToolkit,
        sanitize: bool,
        remove_hs: bool,
    ) -> (
        Result<InchiToMolResult, InchiToMolError>,
        ScriptedStructureEngine,
        ExtraInchiReturnValues,
    ) {
        let mut engine = ScriptedStructureEngine::new(output);
        let mut return_values = ExtraInchiReturnValues {
            return_code: -77,
            message: b"old-message".to_vec(),
            log: b"old-log".to_vec(),
            aux_info: b"old-aux".to_vec(),
        };
        let result = inchi_to_mol(
            &mut engine,
            toolkit,
            b"InChI=1S/test\0tail",
            &mut return_values,
            sanitize,
            remove_hs,
        );
        (result, engine, return_values)
    }

    #[derive(Clone)]
    struct ScriptedGenerationEngine {
        output: AdapterInchiGenerationOutput,
        seen_inputs: Vec<AdapterInchiGenerationInput>,
        calls: Vec<&'static str>,
        get_error: Option<SourceHeapError>,
        free_error: Option<SourceHeapError>,
    }

    impl Default for ScriptedGenerationEngine {
        fn default() -> Self {
            Self {
                output: AdapterInchiGenerationOutput {
                    return_code: 17,
                    inchi: Some(b"InChI=1S/scripted".to_vec()),
                    message: Some(b"new-message".to_vec()),
                    log: Some(b"new-log".to_vec()),
                    aux_info: Some(b"new-aux".to_vec()),
                },
                seen_inputs: Vec::new(),
                calls: Vec::new(),
                get_error: None,
                free_error: None,
            }
        }
    }

    impl InchiGenerationEngine for ScriptedGenerationEngine {
        fn get_inchi(
            &mut self,
            input: &AdapterInchiGenerationInput,
        ) -> Result<AdapterInchiGenerationOutput, SourceHeapError> {
            self.calls.push("GetINCHI");
            self.seen_inputs.push(input.clone());
            if let Some(error) = self.get_error {
                return Err(error);
            }
            Ok(self.output.clone())
        }

        fn free_inchi(&mut self) -> Result<(), SourceHeapError> {
            self.calls.push("FreeINCHI");
            if let Some(error) = self.free_error {
                return Err(error);
            }
            Ok(())
        }
    }

    #[derive(Default)]
    struct RecordingGenerationToolkit {
        calls: Vec<String>,
        needs_update: bool,
        fail_on: Option<String>,
        hydrogens: Vec<u32>,
        total_degrees: Vec<Option<u32>>,
        symbol_override: Option<Vec<u8>>,
        weight_override: Option<f64>,
    }

    impl RecordingGenerationToolkit {
        fn record(&mut self, call: &'static str) -> Result<(), AdapterToolkitError> {
            self.calls.push(call.to_owned());
            if self.fail_on.as_deref() == Some(call) {
                return Err(AdapterToolkitError {
                    kind: "scripted toolkit error",
                    message: call.to_owned(),
                });
            }
            Ok(())
        }

        fn hydrogen_count(&self, atom_index: u32) -> u32 {
            self.hydrogens.get(atom_index as usize).copied().unwrap_or(0)
        }
    }

    impl MolToInchiToolkit for RecordingGenerationToolkit {
        fn needs_update_property_cache(&mut self, _molecule: &AdapterMol) -> Result<bool, AdapterToolkitError> {
            self.record("needs_update_property_cache")?;
            Ok(self.needs_update)
        }

        fn update_property_cache(
            &mut self,
            _molecule: &mut AdapterMol,
            strict: bool,
        ) -> Result<(), AdapterToolkitError> {
            assert!(!strict);
            self.record("update_property_cache")
        }

        fn kekulize(&mut self, _molecule: &mut AdapterMol, mark_atoms_bonds: bool) -> Result<(), AdapterToolkitError> {
            assert!(!mark_atoms_bonds);
            self.record("kekulize")
        }

        fn element_symbol(&mut self, atomic_number: i32) -> Result<Vec<u8>, AdapterToolkitError> {
            self.record("element_symbol")?;
            if let Some(symbol) = &self.symbol_override {
                return Ok(symbol.clone());
            }
            match atomic_number {
                1 => Ok(b"H".to_vec()),
                6 => Ok(b"C".to_vec()),
                7 => Ok(b"N".to_vec()),
                8 => Ok(b"O".to_vec()),
                9 => Ok(b"F".to_vec()),
                15 => Ok(b"P".to_vec()),
                17 => Ok(b"Cl".to_vec()),
                35 => Ok(b"Br".to_vec()),
                53 => Ok(b"I".to_vec()),
                _ => Err(AdapterToolkitError {
                    kind: "Pre-condition Violation",
                    message: "Atomic number not found".to_owned(),
                }),
            }
        }

        fn atomic_weight(&mut self, atomic_number: i32) -> Result<f64, AdapterToolkitError> {
            self.record("atomic_weight")?;
            if let Some(weight) = self.weight_override {
                return Ok(weight);
            }
            match atomic_number {
                1 => Ok(1.008),
                6 => Ok(12.011),
                7 => Ok(14.007),
                8 => Ok(15.999),
                9 => Ok(18.998),
                15 => Ok(30.5),
                17 => Ok(35.45),
                35 => Ok(79.904),
                53 => Ok(126.904),
                _ => Err(AdapterToolkitError {
                    kind: "Pre-condition Violation",
                    message: "Atomic number not found".to_owned(),
                }),
            }
        }

        fn total_num_hydrogens(&mut self, _molecule: &AdapterMol, atom_index: u32) -> Result<u32, AdapterToolkitError> {
            self.record("total_num_hydrogens")?;
            Ok(self.hydrogen_count(atom_index))
        }

        fn calc_implicit_valence(
            &mut self,
            _molecule: &mut AdapterMol,
            _atom_index: u32,
        ) -> Result<i32, AdapterToolkitError> {
            self.record("calc_implicit_valence")?;
            Ok(0)
        }

        fn total_degree(&mut self, molecule: &AdapterMol, atom_index: u32) -> Result<u32, AdapterToolkitError> {
            self.record("total_degree")?;
            Ok(self
                .total_degrees
                .get(atom_index as usize)
                .and_then(|degree| *degree)
                .unwrap_or_else(|| {
                    molecule.adjacency[atom_index as usize].len() as u32 + self.hydrogen_count(atom_index)
                }))
        }
    }

    fn generation_values() -> ExtraInchiReturnValues {
        ExtraInchiReturnValues {
            return_code: -77,
            message: b"old-message".to_vec(),
            log: b"old-log".to_vec(),
            aux_info: b"old-aux".to_vec(),
        }
    }

    fn run_generation(
        molecule: &AdapterMol,
        toolkit: &mut RecordingGenerationToolkit,
        engine: &mut ScriptedGenerationEngine,
        options: Option<&[u8]>,
    ) -> (Result<MolToInchiResult, MolToInchiError>, ExtraInchiReturnValues) {
        let mut values = generation_values();
        let result = mol_to_inchi(engine, toolkit, molecule, &mut values, options);
        (result, values)
    }

    #[test]
    fn source_port__inchi__moltoinchi__line_1747() {
        let mut atoms = [6, 7, 8, 9, 17, 35, 53, 15]
            .into_iter()
            .map(|atomic_number| AdapterAtom {
                atomic_number,
                ..AdapterAtom::default()
            })
            .collect::<Vec<_>>();
        atoms[0].isotope = 13;
        atoms[0].formal_charge = 130;
        atoms[0].num_radical_electrons = 1;
        atoms[7].isotope = 31;
        atoms[7].formal_charge = -130;
        let atom_molecule = AdapterMol::from_graph(atoms, Vec::new());
        let atom_before = atom_molecule.clone();
        let mut atom_toolkit = RecordingGenerationToolkit {
            needs_update: true,
            hydrogens: vec![2, 3, 4, 5, 6, 7, 8, 9],
            ..RecordingGenerationToolkit::default()
        };
        let mut atom_engine = ScriptedGenerationEngine::default();
        atom_engine.output.message = None;
        atom_engine.output.aux_info = None;
        let (atom_result, atom_values) = run_generation(
            &atom_molecule,
            &mut atom_toolkit,
            &mut atom_engine,
            Some(b"/AuxNone -FixedH\0ignored"),
        );
        assert_eq!(
            atom_result,
            Ok(MolToInchiResult {
                inchi: b"InChI=1S/scripted".to_vec(),
                diagnostics: Vec::new(),
            })
        );
        assert_eq!(atom_molecule, atom_before);
        assert_eq!(atom_engine.calls, ["GetINCHI", "FreeINCHI"]);
        assert_eq!(atom_engine.seen_inputs.len(), 1);
        let atom_input = &atom_engine.seen_inputs[0];
        assert_eq!(atom_input.options_with_nul, Some(b"-AuxNone -FixedH\0".to_vec()));
        assert!(atom_input.stereo0d.is_empty());
        assert_eq!(atom_input.atoms.len(), 8);
        assert!(
            atom_input
                .atoms
                .iter()
                .all(|atom| (atom.x, atom.y, atom.z) == (0.0, 0.0, 0.0))
        );
        assert_eq!(atom_input.atoms[0].isotopic_mass, 10001);
        assert_eq!(atom_input.atoms[7].isotopic_mass, 10000);
        assert_eq!(atom_input.atoms[0].charge, -126);
        assert_eq!(atom_input.atoms[7].charge, 126);
        assert_eq!(atom_input.atoms[0].num_iso_H, [2, 0, 0, 0]);
        for index in 1..7 {
            assert_eq!(atom_input.atoms[index].num_iso_H, [-1, 0, 0, 0]);
        }
        assert_eq!(atom_input.atoms[7].num_iso_H, [9, 0, 0, 0]);
        assert!(atom_input.atoms.iter().all(|atom| atom.radical == 0));
        assert_eq!(
            atom_toolkit.calls[..3],
            ["needs_update_property_cache", "update_property_cache", "kekulize"]
        );
        assert_eq!(
            atom_toolkit
                .calls
                .iter()
                .filter(|call| call.as_str() == "total_num_hydrogens")
                .count(),
            2
        );
        assert_eq!(atom_values.return_code, 17);
        assert_eq!(atom_values.message, b"old-message");
        assert_eq!(atom_values.log, b"new-log");
        assert_eq!(atom_values.aux_info, b"old-aux");

        let cleanup_molecule = chlorine_star(3, [-1, 0, -1, -1]);
        let cleanup_before = cleanup_molecule.clone();
        let mut cleanup_toolkit = RecordingGenerationToolkit::default();
        let mut cleanup_engine = ScriptedGenerationEngine::default();
        let (cleanup_result, _) = run_generation(&cleanup_molecule, &mut cleanup_toolkit, &mut cleanup_engine, None);
        assert!(cleanup_result.is_ok());
        assert_eq!(cleanup_molecule, cleanup_before);
        assert_eq!(
            cleanup_engine.seen_inputs[0]
                .atoms
                .iter()
                .map(|atom| atom.charge)
                .collect::<Vec<_>>(),
            vec![0, 0, 0, 0, 0]
        );
        assert_eq!(
            cleanup_engine.seen_inputs[0].atoms[0].bond_type[0],
            BondType::Double as i8
        );
        assert_eq!(
            cleanup_engine.seen_inputs[0].atoms[1].bond_type[..3],
            [BondType::Single as i8, BondType::Double as i8, BondType::Double as i8]
        );

        let mut conformer_molecule = AdapterMol::from_graph(
            vec![AdapterAtom {
                atomic_number: 1,
                ..AdapterAtom::default()
            }],
            Vec::new(),
        );
        conformer_molecule.conformers = vec![vec![[1.25, -2.5, 3.75]], vec![[9.0; 3]]];
        let mut conformer_toolkit = RecordingGenerationToolkit::default();
        let mut conformer_engine = ScriptedGenerationEngine::default();
        let (conformer_result, _) =
            run_generation(&conformer_molecule, &mut conformer_toolkit, &mut conformer_engine, None);
        assert!(conformer_result.is_ok());
        assert_eq!(
            (
                conformer_engine.seen_inputs[0].atoms[0].x,
                conformer_engine.seen_inputs[0].atoms[0].y,
                conformer_engine.seen_inputs[0].atoms[0].z,
            ),
            (1.25, -2.5, 3.75)
        );
        assert_eq!(conformer_engine.seen_inputs[0].options_with_nul, None);
        assert!(
            !conformer_toolkit
                .calls
                .iter()
                .any(|call| call == "update_property_cache")
        );

        let mut null_output_toolkit = RecordingGenerationToolkit::default();
        let mut null_output_engine = ScriptedGenerationEngine {
            output: AdapterInchiGenerationOutput {
                return_code: 0,
                inchi: None,
                message: None,
                log: None,
                aux_info: None,
            },
            ..ScriptedGenerationEngine::default()
        };
        let (null_output_result, null_output_values) = run_generation(
            &AdapterMol::from_graph(Vec::new(), Vec::new()),
            &mut null_output_toolkit,
            &mut null_output_engine,
            None,
        );
        assert_eq!(
            null_output_result,
            Ok(MolToInchiResult {
                inchi: Vec::new(),
                diagnostics: Vec::new(),
            })
        );
        assert_eq!(null_output_values.return_code, 0);
        assert_eq!(null_output_values.message, b"old-message");
        assert_eq!(null_output_values.log, b"old-log");
        assert_eq!(null_output_values.aux_info, b"old-aux");

        for (degree, tag, total_degree, expected_parity, expected_neighbors) in [
            (
                4,
                ChiralTag::TetrahedralCw,
                4,
                Some(tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i8),
                [1, 2, 3, 4],
            ),
            (
                4,
                ChiralTag::TetrahedralCcw,
                4,
                Some(tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8),
                [1, 2, 3, 4],
            ),
            (
                3,
                ChiralTag::TetrahedralCw,
                3,
                Some(tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8),
                [0, 1, 2, 3],
            ),
            (
                3,
                ChiralTag::TetrahedralCcw,
                3,
                Some(tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i8),
                [0, 1, 2, 3],
            ),
            (2, ChiralTag::TetrahedralCw, 2, None, [0; 4]),
            (4, ChiralTag::TetrahedralCcw, 5, None, [0; 4]),
        ] {
            let mut tetra_atoms = vec![
                AdapterAtom {
                    atomic_number: 6,
                    chiral_tag: tag,
                    ..AdapterAtom::default()
                };
                degree + 1
            ];
            for atom in &mut tetra_atoms[1..] {
                atom.chiral_tag = ChiralTag::Unspecified;
            }
            let tetra_bonds = (1..=degree)
                .map(|neighbor| AdapterBond {
                    begin_atom_index: 0,
                    end_atom_index: neighbor as u32,
                    bond_type: BondType::Single,
                    ..AdapterBond::default()
                })
                .collect();
            let tetra_molecule = AdapterMol::from_graph(tetra_atoms, tetra_bonds);
            let mut tetra_toolkit = RecordingGenerationToolkit {
                total_degrees: std::iter::once(Some(total_degree))
                    .chain(std::iter::repeat_n(None, degree))
                    .collect(),
                ..RecordingGenerationToolkit::default()
            };
            let mut tetra_engine = ScriptedGenerationEngine::default();
            let (tetra_result, _) = run_generation(&tetra_molecule, &mut tetra_toolkit, &mut tetra_engine, None);
            let tetra_result = tetra_result.unwrap();
            let stereo = &tetra_engine.seen_inputs[0].stereo0d;
            match expected_parity {
                Some(parity) => {
                    assert_eq!(stereo.len(), 1);
                    assert_eq!(stereo[0].central_atom, 0);
                    assert_eq!(stereo[0].type_, tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8);
                    assert_eq!(stereo[0].parity, parity);
                    assert_eq!(stereo[0].neighbor, expected_neighbors);
                    assert!(tetra_result.diagnostics.is_empty());
                }
                None => {
                    assert!(stereo.is_empty());
                    assert_eq!(
                        tetra_result.diagnostics,
                        [AdapterDiagnostic {
                            level: AdapterDiagnosticLevel::Warning,
                            message: "tetrahedral chirality on atom with <3 or >4 neighbors will be ignored.\n"
                                .to_owned(),
                        }]
                    );
                }
            }
            let calc_position = tetra_toolkit
                .calls
                .iter()
                .position(|call| call == "calc_implicit_valence")
                .unwrap();
            let degree_position = tetra_toolkit
                .calls
                .iter()
                .position(|call| call == "total_degree")
                .unwrap();
            assert!(calc_position < degree_position);
        }

        let all_bond_types = [
            BondType::Unspecified,
            BondType::Single,
            BondType::Double,
            BondType::Triple,
            BondType::Quadruple,
            BondType::Quintuple,
            BondType::Hextuple,
            BondType::OneAndAHalf,
            BondType::TwoAndAHalf,
            BondType::ThreeAndAHalf,
            BondType::FourAndAHalf,
            BondType::FiveAndAHalf,
            BondType::Aromatic,
            BondType::Ionic,
            BondType::Hydrogen,
            BondType::ThreeCenter,
            BondType::DativeOne,
            BondType::Dative,
            BondType::DativeL,
            BondType::DativeR,
            BondType::Other,
            BondType::Zero,
        ];
        for bond_type in all_bond_types {
            let bond_molecule = AdapterMol::from_graph(
                vec![
                    AdapterAtom {
                        atomic_number: 6,
                        ..AdapterAtom::default()
                    };
                    2
                ],
                vec![AdapterBond {
                    begin_atom_index: 0,
                    end_atom_index: 1,
                    bond_type,
                    ..AdapterBond::default()
                }],
            );
            let mut bond_toolkit = RecordingGenerationToolkit::default();
            let mut bond_engine = ScriptedGenerationEngine::default();
            let (bond_result, _) = run_generation(&bond_molecule, &mut bond_toolkit, &mut bond_engine, None);
            let bond_result = bond_result.unwrap();
            assert_eq!(bond_engine.seen_inputs[0].atoms[0].num_bonds, 1);
            assert_eq!(bond_engine.seen_inputs[0].atoms[0].neighbor[0], 1);
            assert_eq!(
                bond_engine.seen_inputs[0].atoms[0].bond_type[0],
                if bond_type as u8 <= BondType::Triple as u8 {
                    bond_type as i8
                } else {
                    BondType::Unspecified as i8
                }
            );
            assert_eq!(
                bond_result.diagnostics.len(),
                usize::from(bond_type as u8 > BondType::Triple as u8)
            );
        }

        for direction in [
            BondDirection::None,
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
            BondDirection::EndDownRight,
            BondDirection::EndUpRight,
            BondDirection::EitherDouble,
            BondDirection::Unknown,
        ] {
            for reversed in [false, true] {
                let direction_molecule = AdapterMol::from_graph(
                    vec![
                        AdapterAtom {
                            atomic_number: 6,
                            ..AdapterAtom::default()
                        };
                        2
                    ],
                    vec![AdapterBond {
                        begin_atom_index: u32::from(reversed),
                        end_atom_index: u32::from(!reversed),
                        bond_type: BondType::Single,
                        direction,
                        ..AdapterBond::default()
                    }],
                );
                let mut direction_toolkit = RecordingGenerationToolkit::default();
                let mut direction_engine = ScriptedGenerationEngine::default();
                let (direction_result, _) =
                    run_generation(&direction_molecule, &mut direction_toolkit, &mut direction_engine, None);
                assert!(direction_result.is_ok());
                let modifier = if reversed { -1 } else { 1 };
                let expected = match direction {
                    BondDirection::BeginWedge => modifier * tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1UP,
                    BondDirection::BeginDash => modifier * tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1DOWN,
                    BondDirection::EitherDouble => tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER,
                    BondDirection::Unknown => modifier * tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1EITHER,
                    _ => tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE,
                };
                assert_eq!(direction_engine.seen_inputs[0].atoms[0].bond_stereo[0], expected as i8);
            }
        }

        for (stereo_tag, expected_parity) in [
            (BondStereo::Z, tagINCHIStereoParity0D_INCHI_PARITY_ODD),
            (BondStereo::Cis, tagINCHIStereoParity0D_INCHI_PARITY_ODD),
            (BondStereo::E, tagINCHIStereoParity0D_INCHI_PARITY_EVEN),
            (BondStereo::Trans, tagINCHIStereoParity0D_INCHI_PARITY_EVEN),
            (BondStereo::AtropCw, tagINCHIStereoParity0D_INCHI_PARITY_EVEN),
            (BondStereo::AtropCcw, tagINCHIStereoParity0D_INCHI_PARITY_EVEN),
        ] {
            let stereo_molecule = AdapterMol::from_graph(
                vec![
                    AdapterAtom {
                        atomic_number: 6,
                        ..AdapterAtom::default()
                    };
                    4
                ],
                vec![
                    AdapterBond {
                        begin_atom_index: 0,
                        end_atom_index: 1,
                        bond_type: BondType::Single,
                        ..AdapterBond::default()
                    },
                    AdapterBond {
                        begin_atom_index: 1,
                        end_atom_index: 2,
                        bond_type: BondType::Double,
                        stereo: stereo_tag,
                        stereo_atoms: vec![0, 3],
                        ..AdapterBond::default()
                    },
                    AdapterBond {
                        begin_atom_index: 2,
                        end_atom_index: 3,
                        bond_type: BondType::Single,
                        ..AdapterBond::default()
                    },
                ],
            );
            let mut stereo_toolkit = RecordingGenerationToolkit::default();
            let mut stereo_engine = ScriptedGenerationEngine::default();
            let (stereo_result, _) = run_generation(&stereo_molecule, &mut stereo_toolkit, &mut stereo_engine, None);
            assert!(stereo_result.is_ok());
            assert_eq!(stereo_engine.seen_inputs[0].stereo0d.len(), 1);
            assert_eq!(
                stereo_engine.seen_inputs[0].stereo0d[0],
                inchi_Stereo0D {
                    neighbor: [0, 1, 2, 3],
                    central_atom: NO_ATOM as i16,
                    type_: tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
                    parity: expected_parity as i8,
                }
            );
        }

        let short_stereo_molecule = AdapterMol::from_graph(
            vec![
                AdapterAtom {
                    atomic_number: 6,
                    ..AdapterAtom::default()
                };
                2
            ],
            vec![AdapterBond {
                begin_atom_index: 0,
                end_atom_index: 1,
                bond_type: BondType::Double,
                stereo: BondStereo::E,
                stereo_atoms: vec![0],
                ..AdapterBond::default()
            }],
        );
        let mut short_stereo_toolkit = RecordingGenerationToolkit::default();
        let mut short_stereo_engine = ScriptedGenerationEngine::default();
        let (short_stereo_result, _) = run_generation(
            &short_stereo_molecule,
            &mut short_stereo_toolkit,
            &mut short_stereo_engine,
            None,
        );
        assert!(short_stereo_result.is_ok());
        assert!(short_stereo_engine.seen_inputs[0].stereo0d.is_empty());

        let mut any_molecule = AdapterMol::from_graph(
            vec![
                AdapterAtom {
                    atomic_number: 6,
                    ..AdapterAtom::default()
                };
                2
            ],
            vec![AdapterBond {
                begin_atom_index: 1,
                end_atom_index: 0,
                bond_type: BondType::Double,
                stereo: BondStereo::Any,
                stereo_atoms: vec![0, 1],
                ..AdapterBond::default()
            }],
        );
        any_molecule.conformers = vec![vec![[1.0, 2.0, 3.0], [7.0, 8.0, 9.0]]];
        let mut any_toolkit = RecordingGenerationToolkit::default();
        let mut any_engine = ScriptedGenerationEngine::default();
        let (any_result, _) = run_generation(&any_molecule, &mut any_toolkit, &mut any_engine, None);
        assert!(any_result.is_ok());
        assert!(any_engine.seen_inputs[0].stereo0d.is_empty());
        assert_eq!(
            (
                any_engine.seen_inputs[0].atoms[0].x,
                any_engine.seen_inputs[0].atoms[0].y,
                any_engine.seen_inputs[0].atoms[0].z,
            ),
            (7.0, 8.0, 9.0)
        );

        let mut swap_molecule = AdapterMol::from_graph(
            vec![
                AdapterAtom {
                    atomic_number: 6,
                    ..AdapterAtom::default()
                };
                4
            ],
            vec![AdapterBond {
                begin_atom_index: 1,
                end_atom_index: 2,
                bond_type: BondType::Double,
                stereo: BondStereo::E,
                stereo_atoms: vec![3, 0],
                ..AdapterBond::default()
            }],
        );
        swap_molecule.bonds.push(AdapterBond {
            begin_atom_index: 0,
            end_atom_index: 1,
            bond_type: BondType::Single,
            ..AdapterBond::default()
        });
        swap_molecule = AdapterMol::from_graph(swap_molecule.atoms, swap_molecule.bonds);
        let mut swap_toolkit = RecordingGenerationToolkit::default();
        let mut swap_engine = ScriptedGenerationEngine::default();
        let (swap_result, _) = run_generation(&swap_molecule, &mut swap_toolkit, &mut swap_engine, None);
        assert!(swap_result.is_ok());
        assert_eq!(swap_engine.seen_inputs[0].stereo0d[0].neighbor, [0, 1, 2, 3]);

        let max_atoms = vec![
            AdapterAtom {
                atomic_number: 6,
                ..AdapterAtom::default()
            };
            MAXVAL as usize + 2
        ];
        let max_bonds = (1..=MAXVAL + 1)
            .map(|end_atom_index| AdapterBond {
                begin_atom_index: 0,
                end_atom_index,
                bond_type: BondType::Single,
                ..AdapterBond::default()
            })
            .collect();
        let max_molecule = AdapterMol::from_graph(max_atoms, max_bonds);
        let mut max_toolkit = RecordingGenerationToolkit::default();
        let mut max_engine = ScriptedGenerationEngine::default();
        let (max_result, max_values) = run_generation(&max_molecule, &mut max_toolkit, &mut max_engine, None);
        assert_eq!(
            max_result,
            Ok(MolToInchiResult {
                inchi: Vec::new(),
                diagnostics: vec![AdapterDiagnostic {
                    level: AdapterDiagnosticLevel::Error,
                    message: format!(
                        " atom 0 has too many bonds: {}. The InChI library supports at most {}\n",
                        MAXVAL, MAXVAL
                    ),
                }],
            })
        );
        assert!(max_engine.calls.is_empty());
        assert_eq!(max_values, generation_values());

        let mut invalid_conformer = conformer_molecule.clone();
        invalid_conformer.conformers[0].clear();
        let mut invalid_toolkit = RecordingGenerationToolkit::default();
        let mut invalid_engine = ScriptedGenerationEngine::default();
        let (invalid_result, invalid_values) =
            run_generation(&invalid_conformer, &mut invalid_toolkit, &mut invalid_engine, None);
        assert_eq!(invalid_result, Err(MolToInchiError::InvalidConformer));
        assert!(invalid_engine.calls.is_empty());
        assert_eq!(invalid_values, generation_values());

        for symbol in [b"Carbon".to_vec(), b"C\0".to_vec()] {
            let mut symbol_toolkit = RecordingGenerationToolkit {
                symbol_override: Some(symbol),
                ..RecordingGenerationToolkit::default()
            };
            let mut symbol_engine = ScriptedGenerationEngine::default();
            let (symbol_result, _) = run_generation(&conformer_molecule, &mut symbol_toolkit, &mut symbol_engine, None);
            assert_eq!(symbol_result, Err(MolToInchiError::ElementSymbolTooLong));
            assert!(symbol_engine.calls.is_empty());
        }

        let mut option_toolkit = RecordingGenerationToolkit::default();
        let mut option_engine = ScriptedGenerationEngine::default();
        let (option_result, _) = run_generation(
            &conformer_molecule,
            &mut option_toolkit,
            &mut option_engine,
            Some(b"-NoTerminator"),
        );
        assert_eq!(option_result, Err(MolToInchiError::InvalidOptions));
        assert!(option_engine.calls.is_empty());

        let error_cases = [
            (
                "update_property_cache",
                true,
                AdapterAtom {
                    atomic_number: 6,
                    ..AdapterAtom::default()
                },
            ),
            (
                "kekulize",
                false,
                AdapterAtom {
                    atomic_number: 6,
                    ..AdapterAtom::default()
                },
            ),
            (
                "element_symbol",
                false,
                AdapterAtom {
                    atomic_number: 6,
                    ..AdapterAtom::default()
                },
            ),
            (
                "atomic_weight",
                false,
                AdapterAtom {
                    atomic_number: 6,
                    isotope: 13,
                    ..AdapterAtom::default()
                },
            ),
            (
                "total_num_hydrogens",
                false,
                AdapterAtom {
                    atomic_number: 15,
                    ..AdapterAtom::default()
                },
            ),
            (
                "calc_implicit_valence",
                false,
                AdapterAtom {
                    atomic_number: 6,
                    chiral_tag: ChiralTag::TetrahedralCw,
                    ..AdapterAtom::default()
                },
            ),
            (
                "total_degree",
                false,
                AdapterAtom {
                    atomic_number: 6,
                    chiral_tag: ChiralTag::TetrahedralCw,
                    ..AdapterAtom::default()
                },
            ),
        ];
        for (failure, needs_update, atom) in error_cases {
            let error_molecule = AdapterMol::from_graph(vec![atom], Vec::new());
            let mut error_toolkit = RecordingGenerationToolkit {
                needs_update,
                fail_on: Some(failure.to_owned()),
                total_degrees: vec![Some(3)],
                ..RecordingGenerationToolkit::default()
            };
            let mut error_engine = ScriptedGenerationEngine::default();
            let (error_result, error_values) =
                run_generation(&error_molecule, &mut error_toolkit, &mut error_engine, None);
            assert!(matches!(error_result, Err(MolToInchiError::Toolkit(_))));
            assert!(error_engine.calls.is_empty());
            assert_eq!(error_values, generation_values());
        }

        let empty_molecule = AdapterMol::from_graph(Vec::new(), Vec::new());
        let mut get_error_toolkit = RecordingGenerationToolkit::default();
        let mut get_error_engine = ScriptedGenerationEngine {
            get_error: Some(SourceHeapError::AllocationFailed),
            ..ScriptedGenerationEngine::default()
        };
        let (get_error_result, get_error_values) =
            run_generation(&empty_molecule, &mut get_error_toolkit, &mut get_error_engine, None);
        assert_eq!(
            get_error_result,
            Err(MolToInchiError::Source(SourceHeapError::AllocationFailed))
        );
        assert_eq!(get_error_engine.calls, ["GetINCHI"]);
        assert_eq!(get_error_values, generation_values());

        let mut free_error_toolkit = RecordingGenerationToolkit::default();
        let mut free_error_engine = ScriptedGenerationEngine {
            free_error: Some(SourceHeapError::MissingAllocation),
            ..ScriptedGenerationEngine::default()
        };
        let (free_error_result, free_error_values) =
            run_generation(&empty_molecule, &mut free_error_toolkit, &mut free_error_engine, None);
        assert_eq!(
            free_error_result,
            Err(MolToInchiError::Source(SourceHeapError::MissingAllocation))
        );
        assert_eq!(free_error_engine.calls, ["GetINCHI", "FreeINCHI"]);
        assert_eq!(free_error_values.return_code, 17);
        assert_eq!(free_error_values.message, b"new-message");
        assert_eq!(free_error_values.log, b"new-log");
        assert_eq!(free_error_values.aux_info, b"new-aux");

        let allocation_input = AdapterInchiGenerationInput {
            atoms: vec![inchi_Atom::default()],
            stereo0d: vec![inchi_Stereo0D::default()],
            options_with_nul: Some(b"-AuxNone\0".to_vec()),
        };
        for successful_allocations in 0..=2 {
            let mut heap = SourceHeap::default();
            heap.fail_after_allocations(successful_allocations);
            let mut source_engine = SourceInchiGenerationEngine {
                heap: &mut heap,
                stdout: SourceMutPointer::null(),
                build: InchiBuildMetadata {
                    compiler: "gcc",
                    date: "Jan 01 1970",
                    time: "00:00:00",
                },
                clock_result: 1_000,
                pending_output: None,
            };
            assert_eq!(
                source_engine.get_inchi(&allocation_input),
                Err(SourceHeapError::AllocationFailed)
            );
            assert!(source_engine.pending_output.is_none());
            drop(source_engine);
            assert_eq!(heap.live_source_allocation_count(), 0);
            assert_eq!(heap.live_allocation_count(), 0);
        }
    }

    fn mol_to_inchi_oracle_u32(value: &serde_json::Value, field: &str, case_id: &str) -> u32 {
        u32::try_from(
            value
                .as_u64()
                .unwrap_or_else(|| panic!("{case_id}: {field} must be an unsigned integer")),
        )
        .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds u32"))
    }

    fn mol_to_inchi_oracle_i32(value: &serde_json::Value, field: &str, case_id: &str) -> i32 {
        i32::try_from(
            value
                .as_i64()
                .unwrap_or_else(|| panic!("{case_id}: {field} must be a signed integer")),
        )
        .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds i32"))
    }

    fn mol_to_inchi_oracle_bytes(value: &serde_json::Value, case_id: &str) -> Vec<u8> {
        value
            .as_array()
            .unwrap_or_else(|| panic!("{case_id}: byte field must be an array"))
            .iter()
            .map(|byte| {
                u8::try_from(
                    byte.as_u64()
                        .unwrap_or_else(|| panic!("{case_id}: byte must be unsigned")),
                )
                .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
            })
            .collect()
    }

    fn mol_to_inchi_oracle_bond_type(value: &serde_json::Value, case_id: &str) -> BondType {
        match mol_to_inchi_oracle_u32(value, "bond_type", case_id) {
            0 => BondType::Unspecified,
            1 => BondType::Single,
            2 => BondType::Double,
            3 => BondType::Triple,
            4 => BondType::Quadruple,
            5 => BondType::Quintuple,
            6 => BondType::Hextuple,
            7 => BondType::OneAndAHalf,
            8 => BondType::TwoAndAHalf,
            9 => BondType::ThreeAndAHalf,
            10 => BondType::FourAndAHalf,
            11 => BondType::FiveAndAHalf,
            12 => BondType::Aromatic,
            13 => BondType::Ionic,
            14 => BondType::Hydrogen,
            15 => BondType::ThreeCenter,
            16 => BondType::DativeOne,
            17 => BondType::Dative,
            18 => BondType::DativeL,
            19 => BondType::DativeR,
            20 => BondType::Other,
            21 => BondType::Zero,
            other => panic!("{case_id}: unknown bond type {other}"),
        }
    }

    fn mol_to_inchi_oracle_direction(value: &serde_json::Value, case_id: &str) -> BondDirection {
        match mol_to_inchi_oracle_u32(value, "direction", case_id) {
            0 => BondDirection::None,
            1 => BondDirection::BeginWedge,
            2 => BondDirection::BeginDash,
            3 => BondDirection::EndDownRight,
            4 => BondDirection::EndUpRight,
            5 => BondDirection::EitherDouble,
            6 => BondDirection::Unknown,
            other => panic!("{case_id}: unknown bond direction {other}"),
        }
    }

    fn mol_to_inchi_oracle_stereo(value: &serde_json::Value, case_id: &str) -> BondStereo {
        match mol_to_inchi_oracle_u32(value, "stereo", case_id) {
            0 => BondStereo::None,
            1 => BondStereo::Any,
            2 => BondStereo::Z,
            3 => BondStereo::E,
            4 => BondStereo::Cis,
            5 => BondStereo::Trans,
            6 => BondStereo::AtropCw,
            7 => BondStereo::AtropCcw,
            other => panic!("{case_id}: unknown bond stereo {other}"),
        }
    }

    fn mol_to_inchi_oracle_chiral_tag(value: &serde_json::Value, case_id: &str) -> ChiralTag {
        match mol_to_inchi_oracle_u32(value, "chiral_tag", case_id) {
            0 => ChiralTag::Unspecified,
            1 => ChiralTag::TetrahedralCw,
            2 => ChiralTag::TetrahedralCcw,
            other => panic!("{case_id}: unknown chiral tag {other}"),
        }
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__moltoinchi__exact() {
        use serde_json::{Value, json};
        use std::process::Command;

        let repository_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(std::path::Path::parent)
            .expect("crate must be below repository root");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(runner)
            .arg("--mol-to-inchi-records")
            .current_dir(repository_root)
            .output()
            .expect("pinned RDKit C++ oracle must run");
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let records = output.lines().collect::<Vec<_>>();
        assert_eq!(records.len(), 65, "MolToInchi oracle case count changed");

        for line in records {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            let case_id = official["case_id"].as_str().expect("case_id must be a string");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "000d4eff4353e06aaa802057674fa24f91726872b0c0f9bb593d7a1937721ac2"
            );
            assert_eq!(official["operation"], "MolToInchi");
            assert_eq!(official.as_object().unwrap().len(), 8, "{case_id}: record fields");

            let input = &official["input"];
            assert_eq!(input.as_object().unwrap().len(), 7, "{case_id}: input fields");
            let atom_fields = input["atoms"].as_array().expect("atoms must be an array");
            let mut hydrogens = Vec::with_capacity(atom_fields.len());
            let mut total_degrees = Vec::with_capacity(atom_fields.len());
            let atoms = atom_fields
                .iter()
                .map(|field| {
                    assert_eq!(field.as_object().unwrap().len(), 10, "{case_id}: atom fields");
                    hydrogens.push(mol_to_inchi_oracle_u32(
                        &field["total_hydrogens"],
                        "total_hydrogens",
                        case_id,
                    ));
                    total_degrees.push(if field["total_degree"].is_null() {
                        None
                    } else {
                        Some(mol_to_inchi_oracle_u32(&field["total_degree"], "total_degree", case_id))
                    });
                    AdapterAtom {
                        atomic_number: mol_to_inchi_oracle_i32(&field["atomic_number"], "atomic_number", case_id),
                        formal_charge: mol_to_inchi_oracle_i32(&field["formal_charge"], "formal_charge", case_id),
                        num_explicit_hydrogens: mol_to_inchi_oracle_u32(
                            &field["explicit_hydrogens"],
                            "explicit_hydrogens",
                            case_id,
                        ),
                        is_aromatic: field["aromatic"].as_bool().expect("aromatic must be bool"),
                        isotope: mol_to_inchi_oracle_u32(&field["isotope"], "isotope", case_id),
                        num_radical_electrons: mol_to_inchi_oracle_u32(
                            &field["radical_electrons"],
                            "radical_electrons",
                            case_id,
                        ),
                        no_implicit: field["no_implicit"].as_bool().expect("no_implicit must be bool"),
                        chiral_tag: mol_to_inchi_oracle_chiral_tag(&field["chiral_tag"], case_id),
                        cip_rank: None,
                        cached_explicit_valence: None,
                        cached_implicit_valence: None,
                    }
                })
                .collect::<Vec<_>>();
            let bonds = input["bonds"]
                .as_array()
                .expect("bonds must be an array")
                .iter()
                .map(|field| {
                    assert_eq!(field.as_object().unwrap().len(), 7, "{case_id}: bond fields");
                    AdapterBond {
                        begin_atom_index: mol_to_inchi_oracle_u32(
                            &field["begin_atom_index"],
                            "begin_atom_index",
                            case_id,
                        ),
                        end_atom_index: mol_to_inchi_oracle_u32(&field["end_atom_index"], "end_atom_index", case_id),
                        bond_type: mol_to_inchi_oracle_bond_type(&field["bond_type"], case_id),
                        direction: mol_to_inchi_oracle_direction(&field["direction"], case_id),
                        is_aromatic: field["aromatic"].as_bool().expect("aromatic must be bool"),
                        stereo: mol_to_inchi_oracle_stereo(&field["stereo"], case_id),
                        stereo_atoms: field["stereo_atoms"]
                            .as_array()
                            .expect("stereo_atoms must be an array")
                            .iter()
                            .map(|value| mol_to_inchi_oracle_u32(value, "stereo_atom", case_id))
                            .collect(),
                    }
                })
                .collect::<Vec<_>>();
            let mut molecule = AdapterMol::from_graph(atoms, bonds);
            molecule.conformers = input["conformers"]
                .as_array()
                .expect("conformers must be an array")
                .iter()
                .map(|conformer| {
                    conformer
                        .as_array()
                        .expect("conformer must be an array")
                        .iter()
                        .map(|point| {
                            let point = point.as_array().expect("coordinate must be an array");
                            assert_eq!(point.len(), 3, "{case_id}: coordinate width");
                            [
                                point[0].as_f64().expect("x must be numeric"),
                                point[1].as_f64().expect("y must be numeric"),
                                point[2].as_f64().expect("z must be numeric"),
                            ]
                        })
                        .collect()
                })
                .collect();
            let molecule_before = molecule.clone();

            let scripted = &input["scripted_output"];
            assert_eq!(
                scripted.as_object().unwrap().len(),
                7,
                "{case_id}: scripted output fields"
            );
            let optional_bytes = |value: &Value| (!value.is_null()).then(|| mol_to_inchi_oracle_bytes(value, case_id));
            let throw_on_get = scripted["throw_on_get"].as_bool().expect("throw_on_get must be bool");
            let throw_on_free = scripted["throw_on_free"].as_bool().expect("throw_on_free must be bool");
            let mut engine = ScriptedGenerationEngine {
                output: AdapterInchiGenerationOutput {
                    return_code: mol_to_inchi_oracle_i32(&scripted["return_code"], "return_code", case_id),
                    inchi: optional_bytes(&scripted["inchi"]),
                    message: optional_bytes(&scripted["message"]),
                    log: optional_bytes(&scripted["log"]),
                    aux_info: optional_bytes(&scripted["aux_info"]),
                },
                get_error: throw_on_get.then_some(SourceHeapError::AllocationFailed),
                free_error: throw_on_free.then_some(SourceHeapError::MissingAllocation),
                ..ScriptedGenerationEngine::default()
            };
            let fail_on = input["fail_on"].as_str().map(str::to_owned);
            let mut toolkit = RecordingGenerationToolkit {
                needs_update: input["needs_update_property_cache"]
                    .as_bool()
                    .expect("needs_update_property_cache must be bool"),
                fail_on,
                hydrogens,
                total_degrees,
                ..RecordingGenerationToolkit::default()
            };
            let options = if input["options"].is_null() {
                None
            } else {
                Some(mol_to_inchi_oracle_bytes(&input["options"], case_id))
            };
            let mut return_values = generation_values();
            let result = mol_to_inchi(
                &mut engine,
                &mut toolkit,
                &molecule,
                &mut return_values,
                options.as_deref(),
            );

            let (status, exception_kind, exception_message, inchi, diagnostics) = match result {
                Ok(result) => ("returned", Value::Null, Value::Null, result.inchi, result.diagnostics),
                Err(MolToInchiError::Toolkit(error)) => (
                    "exception",
                    json!("MolSanitizeException"),
                    json!(error.message),
                    Vec::new(),
                    Vec::new(),
                ),
                Err(MolToInchiError::Source(_)) if throw_on_get => (
                    "exception",
                    json!("runtime_error"),
                    json!("GetINCHI"),
                    Vec::new(),
                    Vec::new(),
                ),
                Err(MolToInchiError::Source(_)) if throw_on_free => (
                    "exception",
                    json!("runtime_error"),
                    json!("FreeINCHI"),
                    Vec::new(),
                    Vec::new(),
                ),
                Err(error) => panic!("{case_id}: unexpected Rust error {error:?}"),
            };
            let warning_text = diagnostics
                .iter()
                .filter(|diagnostic| diagnostic.level == AdapterDiagnosticLevel::Warning)
                .map(|diagnostic| diagnostic.message.as_str())
                .collect::<String>();
            let error_text = diagnostics
                .iter()
                .filter(|diagnostic| diagnostic.level == AdapterDiagnosticLevel::Error)
                .map(|diagnostic| diagnostic.message.as_str())
                .collect::<String>();
            let captured = engine.seen_inputs.first();
            let captured_atoms = captured
                .map(|input| {
                    input
                        .atoms
                        .iter()
                        .map(|atom| {
                            let element = atom
                                .elname
                                .iter()
                                .take_while(|byte| **byte != 0)
                                .map(|byte| *byte as u8)
                                .collect::<Vec<_>>();
                            let count = usize::try_from(atom.num_bonds).expect("num_bonds must be nonnegative");
                            json!({
                                "x": atom.x,
                                "y": atom.y,
                                "z": atom.z,
                                "element": String::from_utf8(element).expect("element must be UTF-8"),
                                "isotopic_mass": atom.isotopic_mass,
                                "charge": atom.charge,
                                "radical": atom.radical,
                                "num_iso_h": atom.num_iso_H,
                                "num_bonds": atom.num_bonds,
                                "neighbors": atom.neighbor[..count],
                                "bond_types": atom.bond_type[..count],
                                "bond_stereo": atom.bond_stereo[..count],
                            })
                        })
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            let captured_stereo = captured
                .map(|input| {
                    input
                        .stereo0d
                        .iter()
                        .map(|stereo| {
                            json!({
                                "neighbors": stereo.neighbor,
                                "central_atom": stereo.central_atom,
                                "type": stereo.type_,
                                "parity": stereo.parity,
                            })
                        })
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            let captured_options = captured
                .and_then(|input| input.options_with_nul.as_ref())
                .map(|options| {
                    let nul = options
                        .iter()
                        .position(|byte| *byte == 0)
                        .expect("captured options must be NUL terminated");
                    options[..nul].to_vec()
                });
            let original_atoms = molecule
                .atoms
                .iter()
                .enumerate()
                .map(|(index, atom)| {
                    json!({
                        "index": index,
                        "atomic_number": atom.atomic_number,
                        "formal_charge": atom.formal_charge,
                        "num_explicit_hydrogens": atom.num_explicit_hydrogens,
                        "is_aromatic": atom.is_aromatic,
                        "isotope": atom.isotope,
                        "num_radical_electrons": atom.num_radical_electrons,
                        "no_implicit": atom.no_implicit,
                        "chiral_tag": atom.chiral_tag as u8,
                        "cip_rank": atom.cip_rank,
                    })
                })
                .collect::<Vec<_>>();
            let original_bonds = molecule
                .bonds
                .iter()
                .enumerate()
                .map(|(index, bond)| {
                    json!({
                        "index": index,
                        "begin_atom_index": bond.begin_atom_index,
                        "end_atom_index": bond.end_atom_index,
                        "bond_type": bond.bond_type as u8,
                        "direction": bond.direction as u8,
                        "is_aromatic": bond.is_aromatic,
                        "stereo": bond.stereo as u8,
                        "stereo_atoms": bond.stereo_atoms,
                    })
                })
                .collect::<Vec<_>>();
            let actual = json!({
                "status": status,
                "exception_kind": exception_kind,
                "exception_message": exception_message,
                "inchi": inchi,
                "return_values": {
                    "return_code": return_values.return_code,
                    "message": return_values.message,
                    "log": return_values.log,
                    "aux_info": return_values.aux_info,
                },
                "warning_text": warning_text,
                "error_text": error_text,
                "toolkit_calls": toolkit.calls,
                "generation_calls": engine.calls,
                "captured_atoms": captured_atoms,
                "captured_stereo0d": captured_stereo,
                "captured_options": captured_options,
                "get_count": engine.calls.iter().filter(|call| **call == "GetINCHI").count(),
                "free_count": engine.calls.iter().filter(|call| **call == "FreeINCHI").count(),
                "outstanding_outputs": 0,
                "original_atoms": original_atoms,
                "original_bonds": original_bonds,
                "original_conformers": molecule.conformers,
            });
            assert_eq!(molecule, molecule_before, "{case_id}: original molecule mutated");
            assert_eq!(
                actual, official["output"],
                "{case_id}: Rust output differs from pinned C++ MolToInchi"
            );
        }
    }

    #[test]
    fn source_port__inchi__fixoptionsymbol__line_1674() {
        let mut empty_output = [0xa5_u8; 3];
        assert_eq!(fix_option_symbol(b"\0tail", &mut empty_output), Ok(()));
        assert_eq!(empty_output, [0, 0xa5, 0xa5]);

        let input = b"/AuxNone -FixedH //keep\xff\0ignored";
        let mut output = [0xa5_u8; 32];
        assert_eq!(fix_option_symbol(input, &mut output), Ok(()));
        let expected = b"-AuxNone -FixedH --keep\xff\0";
        assert_eq!(&output[..expected.len()], expected);
        assert!(output[expected.len()..].iter().all(|byte| *byte == 0xa5));

        let mut exact_output = [0xa5_u8; 4];
        assert_eq!(fix_option_symbol(b"/-x\0", &mut exact_output), Ok(()));
        assert_eq!(exact_output, *b"--x\0");

        let mut missing_nul_output = [0xa5_u8; 8];
        assert_eq!(
            fix_option_symbol(b"/option", &mut missing_nul_output),
            Err(FixOptionSymbolError::InputIsNotNulTerminated)
        );
        assert_eq!(missing_nul_output, [0xa5; 8]);

        let mut short_output = [0xa5_u8; 3];
        assert_eq!(
            fix_option_symbol(b"/ab\0", &mut short_output),
            Err(FixOptionSymbolError::OutputIsTooSmall)
        );
        assert_eq!(short_output, [0xa5; 3]);
    }

    fn chlorine_star(center_charge: i32, oxygen_charges: [i32; 4]) -> AdapterMol {
        AdapterMol::from_graph(
            vec![
                AdapterAtom {
                    atomic_number: 8,
                    formal_charge: oxygen_charges[0],
                    isotope: 18,
                    ..AdapterAtom::default()
                },
                AdapterAtom {
                    atomic_number: 17,
                    formal_charge: center_charge,
                    num_explicit_hydrogens: 1,
                    ..AdapterAtom::default()
                },
                AdapterAtom {
                    atomic_number: 8,
                    formal_charge: oxygen_charges[1],
                    ..AdapterAtom::default()
                },
                AdapterAtom {
                    atomic_number: 8,
                    formal_charge: oxygen_charges[2],
                    ..AdapterAtom::default()
                },
                AdapterAtom {
                    atomic_number: 8,
                    formal_charge: oxygen_charges[3],
                    ..AdapterAtom::default()
                },
            ],
            (0_u32..4)
                .map(|offset| AdapterBond {
                    begin_atom_index: 1,
                    end_atom_index: if offset == 0 { 0 } else { offset + 1 },
                    bond_type: BondType::Single,
                    direction: BondDirection::BeginDash,
                    stereo: BondStereo::Any,
                    stereo_atoms: vec![4, 3],
                    ..AdapterBond::default()
                })
                .collect(),
        )
    }

    #[test]
    fn source_port__inchi__rcleanup__line_1694() {
        let mut no_match = graph(&[(17, 3), (8, -1)], &[(0, 1, BondType::Double)]);
        let no_match_before = no_match.clone();
        r_clean_up(&mut no_match);
        assert_eq!(no_match, no_match_before);

        let mut wrong_center = chlorine_star(2, [-1, -1, -1, 0]);
        assert!(r_cleanup_matches(&wrong_center).is_empty());
        let wrong_center_before = wrong_center.clone();
        r_clean_up(&mut wrong_center);
        assert_eq!(wrong_center, wrong_center_before);

        let mut two_neutral = chlorine_star(3, [0, -1, 0, -1]);
        assert!(r_cleanup_matches(&two_neutral).is_empty());
        let two_neutral_before = two_neutral.clone();
        r_clean_up(&mut two_neutral);
        assert_eq!(two_neutral, two_neutral_before);

        for charges in [[-2, -3, 1, 2], [-2, -3, 1, 0]] {
            let mut charge_miss = chlorine_star(3, charges);
            assert!(r_cleanup_matches(&charge_miss).is_empty());
            let before = charge_miss.clone();
            r_clean_up(&mut charge_miss);
            assert_eq!(charge_miss, before);
        }

        let mut one_neutral = chlorine_star(3, [-1, 0, -1, -1]);
        assert_eq!(r_cleanup_matches(&one_neutral).len(), 1);
        let one_neutral_before = one_neutral.clone();
        r_clean_up(&mut one_neutral);
        assert_eq!(
            one_neutral
                .atoms
                .iter()
                .map(|atom| atom.formal_charge)
                .collect::<Vec<_>>(),
            vec![0, 0, 0, 0, 0]
        );
        assert_eq!(
            one_neutral.bonds.iter().map(|bond| bond.bond_type).collect::<Vec<_>>(),
            vec![BondType::Double, BondType::Single, BondType::Double, BondType::Double]
        );
        for (after, before) in one_neutral.bonds.iter().zip(&one_neutral_before.bonds) {
            assert_eq!(after.direction, before.direction);
            assert_eq!(after.stereo, before.stereo);
            assert_eq!(after.stereo_atoms, before.stereo_atoms);
        }
        assert_eq!(one_neutral.atoms[0].isotope, 18);
        assert_eq!(one_neutral.atoms[1].num_explicit_hydrogens, 1);

        for charges in [[-1, -1, -1, -1], [-1, -1, -1, 2]] {
            let mut no_neutral = chlorine_star(3, charges);
            assert_eq!(r_cleanup_matches(&no_neutral).len(), 1);
            r_clean_up(&mut no_neutral);
            assert_eq!(
                no_neutral
                    .atoms
                    .iter()
                    .map(|atom| atom.formal_charge)
                    .collect::<Vec<_>>(),
                vec![-1, 0, 0, 0, 0]
            );
            assert_eq!(
                no_neutral.bonds.iter().map(|bond| bond.bond_type).collect::<Vec<_>>(),
                vec![BondType::Single, BondType::Double, BondType::Double, BondType::Double]
            );
        }

        let first = chlorine_star(3, [-1, -1, -1, 0]);
        let mut atoms = first.atoms.clone();
        atoms.push(AdapterAtom {
            atomic_number: 8,
            ..AdapterAtom::default()
        });
        let mut bonds = first.bonds.clone();
        bonds.push(AdapterBond {
            begin_atom_index: 1,
            end_atom_index: 5,
            bond_type: BondType::Single,
            direction: BondDirection::EndUpRight,
            is_aromatic: true,
            stereo: BondStereo::Trans,
            stereo_atoms: vec![2, 4],
        });
        let mut partial = AdapterMol::from_graph(atoms, bonds);
        let partial_before = partial.clone();
        assert_eq!(r_cleanup_matches(&partial).len(), 2);
        r_clean_up(&mut partial);
        assert_eq!(partial.atoms[1].formal_charge, 0);
        assert_eq!(partial.atoms[5], partial_before.atoms[5]);
        assert_eq!(partial.bonds[4], partial_before.bonds[4]);

        let second = chlorine_star(3, [-1, -1, -1, 0]);
        let mut atoms = first.atoms;
        atoms.extend(second.atoms.iter().cloned());
        let mut bonds = first.bonds;
        bonds.extend(second.bonds.iter().cloned().map(|mut bond| {
            bond.begin_atom_index += 5;
            bond.end_atom_index += 5;
            bond
        }));
        let mut disconnected = AdapterMol::from_graph(atoms, bonds);
        assert_eq!(r_cleanup_matches(&disconnected).len(), 2);
        r_clean_up(&mut disconnected);
        assert_eq!(disconnected.atoms[1].formal_charge, 0);
        assert_eq!(disconnected.atoms[6].formal_charge, 0);
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__rcleanup__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::{Value, json};

        fn unsigned(value: &Value, field: &str, case_id: &str) -> u32 {
            u32::try_from(
                value
                    .as_u64()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be unsigned")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: {field} must fit u32"))
        }

        fn signed(value: &Value, field: &str, case_id: &str) -> i32 {
            i32::try_from(
                value
                    .as_i64()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be signed")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: {field} must fit i32"))
        }

        fn bond_type(value: &Value, case_id: &str) -> BondType {
            match unsigned(value, "bond_type", case_id) {
                0 => BondType::Unspecified,
                1 => BondType::Single,
                2 => BondType::Double,
                3 => BondType::Triple,
                4 => BondType::Quadruple,
                5 => BondType::Quintuple,
                6 => BondType::Hextuple,
                7 => BondType::OneAndAHalf,
                8 => BondType::TwoAndAHalf,
                9 => BondType::ThreeAndAHalf,
                10 => BondType::FourAndAHalf,
                11 => BondType::FiveAndAHalf,
                12 => BondType::Aromatic,
                13 => BondType::Ionic,
                14 => BondType::Hydrogen,
                15 => BondType::ThreeCenter,
                16 => BondType::DativeOne,
                17 => BondType::Dative,
                18 => BondType::DativeL,
                19 => BondType::DativeR,
                20 => BondType::Other,
                21 => BondType::Zero,
                other => panic!("{case_id}: unknown bond type {other}"),
            }
        }

        fn direction(value: &Value, case_id: &str) -> BondDirection {
            match unsigned(value, "direction", case_id) {
                0 => BondDirection::None,
                1 => BondDirection::BeginWedge,
                2 => BondDirection::BeginDash,
                3 => BondDirection::EndDownRight,
                4 => BondDirection::EndUpRight,
                5 => BondDirection::EitherDouble,
                6 => BondDirection::Unknown,
                other => panic!("{case_id}: unknown bond direction {other}"),
            }
        }

        fn stereo(value: &Value, case_id: &str) -> BondStereo {
            match unsigned(value, "stereo", case_id) {
                0 => BondStereo::None,
                1 => BondStereo::Any,
                2 => BondStereo::Z,
                3 => BondStereo::E,
                4 => BondStereo::Cis,
                5 => BondStereo::Trans,
                6 => BondStereo::AtropCw,
                7 => BondStereo::AtropCcw,
                other => panic!("{case_id}: unknown bond stereo {other}"),
            }
        }

        fn chiral_tag(value: &Value, case_id: &str) -> ChiralTag {
            match unsigned(value, "chiral_tag", case_id) {
                0 => ChiralTag::Unspecified,
                1 => ChiralTag::TetrahedralCw,
                2 => ChiralTag::TetrahedralCcw,
                3 => ChiralTag::Other,
                other => panic!("{case_id}: unknown chiral tag {other}"),
            }
        }

        fn molecule_from_value(value: &Value, case_id: &str) -> AdapterMol {
            let atom_values = value["atoms"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: atoms must be an array"));
            let atoms = atom_values
                .iter()
                .enumerate()
                .map(|(index, atom)| {
                    assert_eq!(
                        unsigned(&atom["index"], "atom index", case_id),
                        index as u32,
                        "{case_id}: contiguous atom indices"
                    );
                    AdapterAtom {
                        atomic_number: signed(&atom["atomic_number"], "atomic_number", case_id),
                        formal_charge: signed(&atom["formal_charge"], "formal_charge", case_id),
                        num_explicit_hydrogens: unsigned(
                            &atom["num_explicit_hydrogens"],
                            "num_explicit_hydrogens",
                            case_id,
                        ),
                        is_aromatic: atom["is_aromatic"]
                            .as_bool()
                            .unwrap_or_else(|| panic!("{case_id}: atom aromatic must be bool")),
                        isotope: unsigned(&atom["isotope"], "isotope", case_id),
                        num_radical_electrons: unsigned(
                            &atom["num_radical_electrons"],
                            "num_radical_electrons",
                            case_id,
                        ),
                        no_implicit: atom["no_implicit"]
                            .as_bool()
                            .unwrap_or_else(|| panic!("{case_id}: no_implicit must be bool")),
                        chiral_tag: chiral_tag(&atom["chiral_tag"], case_id),
                        cip_rank: if atom["cip_rank"].is_null() {
                            None
                        } else {
                            Some(unsigned(&atom["cip_rank"], "cip_rank", case_id))
                        },
                        cached_explicit_valence: None,
                        cached_implicit_valence: None,
                    }
                })
                .collect();
            let bond_values = value["bonds"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: bonds must be an array"));
            let bonds = bond_values
                .iter()
                .enumerate()
                .map(|(index, bond)| {
                    assert_eq!(
                        unsigned(&bond["index"], "bond index", case_id),
                        index as u32,
                        "{case_id}: contiguous bond indices"
                    );
                    AdapterBond {
                        begin_atom_index: unsigned(&bond["begin_atom_index"], "begin_atom_index", case_id),
                        end_atom_index: unsigned(&bond["end_atom_index"], "end_atom_index", case_id),
                        bond_type: bond_type(&bond["bond_type"], case_id),
                        direction: direction(&bond["direction"], case_id),
                        is_aromatic: bond["is_aromatic"]
                            .as_bool()
                            .unwrap_or_else(|| panic!("{case_id}: bond aromatic must be bool")),
                        stereo: stereo(&bond["stereo"], case_id),
                        stereo_atoms: bond["stereo_atoms"]
                            .as_array()
                            .unwrap_or_else(|| panic!("{case_id}: stereo_atoms must be an array"))
                            .iter()
                            .map(|atom| unsigned(atom, "stereo atom", case_id))
                            .collect(),
                    }
                })
                .collect();
            AdapterMol::from_graph(atoms, bonds)
        }

        fn molecule_value(molecule: &AdapterMol) -> Value {
            json!({
                "atoms": molecule.atoms.iter().enumerate().map(|(index, atom)| json!({
                    "index": index,
                    "atomic_number": atom.atomic_number,
                    "formal_charge": atom.formal_charge,
                    "num_explicit_hydrogens": atom.num_explicit_hydrogens,
                    "is_aromatic": atom.is_aromatic,
                    "isotope": atom.isotope,
                    "num_radical_electrons": atom.num_radical_electrons,
                    "no_implicit": atom.no_implicit,
                    "chiral_tag": atom.chiral_tag as u8,
                    "cip_rank": atom.cip_rank,
                })).collect::<Vec<_>>(),
                "bonds": molecule.bonds.iter().enumerate().map(|(index, bond)| json!({
                    "index": index,
                    "begin_atom_index": bond.begin_atom_index,
                    "end_atom_index": bond.end_atom_index,
                    "bond_type": bond.bond_type as u8,
                    "direction": bond.direction as u8,
                    "is_aromatic": bond.is_aromatic,
                    "stereo": bond.stereo as u8,
                    "stereo_atoms": bond.stereo_atoms,
                })).collect::<Vec<_>>(),
            })
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--r-clean-up-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let mut case_ids = BTreeSet::new();

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "81165fe80fd36e56b94f44fc7f062b0971bdd86d96fec1ebbdd947c954c64858"
            );
            assert_eq!(official["operation"], "rCleanUp");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let mut molecule = molecule_from_value(&official["input"], case_id);
            let input = molecule_value(&molecule);
            r_clean_up(&mut molecule);
            let output_molecule = molecule_value(&molecule);
            let actual = json!({
                "schema_version": "cosmolkit-inchi-rdkit-cpp-v1",
                "rdkit_version": "2026.03.1",
                "source_sha256": "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f",
                "source_fragment_sha256": "81165fe80fd36e56b94f44fc7f062b0971bdd86d96fec1ebbdd947c954c64858",
                "dependency_sha256": {
                    "atom_cpp": "bfa51d29a6b18c4bfa373b2803ee79499b7ffbfa40d5838efb273bba2078961e",
                    "smiles_parse_cpp": "270b0457fbee294a9a1a1d7ada856538d4faec96fd782ffff9d2c5bb9ec0ca59",
                    "substruct_match_cpp": "a9c2829c484c98df285fa503870351d1ac8d947cfe46ad38f3cb6664f1c9ede2",
                    "substruct_match_h": "9c0a6d4c56f8fcbc8b4c3144ffd2893547780ebba5602858b2c5415beaaca1cc",
                    "substruct_utils_cpp": "76dbf84efc3b628f206de39b213373d55031e2762b116026c2fd05581af95c5f",
                },
                "operation": "rCleanUp",
                "case_id": case_id,
                "input": input,
                "output": {
                    "status": "returned",
                    "atoms": output_molecule["atoms"].clone(),
                    "bonds": output_molecule["bonds"].clone(),
                },
            });
            assert_eq!(actual, official, "{case_id}: complete observable object");
        }

        assert_eq!(
            case_ids,
            BTreeSet::from([
                "all-negative".to_owned(),
                "arbitrary-nonmatching-nonzero".to_owned(),
                "arbitrary-nonmatching-one-neutral".to_owned(),
                "insufficient-negative-oxygens".to_owned(),
                "neutral-target-0".to_owned(),
                "neutral-target-1".to_owned(),
                "neutral-target-2".to_owned(),
                "neutral-target-3".to_owned(),
                "nonzero-wildcard".to_owned(),
                "overlapping-precollected-matches".to_owned(),
                "preserved-complete-fields".to_owned(),
                "reversed-bond-endpoints".to_owned(),
                "topology-bond-type-miss".to_owned(),
                "two-disconnected-matches".to_owned(),
                "wrong-center-charge".to_owned(),
            ])
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__fixoptionsymbol__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::{Value, json};

        fn bytes(value: &Value, field: &str, case_id: &str) -> Vec<u8> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: {field} must be an array"))
                .iter()
                .map(|byte| {
                    u8::try_from(
                        byte.as_u64()
                            .unwrap_or_else(|| panic!("{case_id}: {field} byte must be unsigned")),
                    )
                    .unwrap_or_else(|_| panic!("{case_id}: {field} byte must fit u8"))
                })
                .collect()
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--fix-option-symbol-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let mut case_ids = BTreeSet::new();

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "5770d6ab655210bd67e5351ef3ebe39aadc0c8b8e96c1b6e2ec0935acb8844d9"
            );
            assert_eq!(official["operation"], "fixOptionSymbol");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let input = bytes(&official["input"]["input_bytes"], "input_bytes", case_id);
            let mut output_bytes = bytes(&official["input"]["output_bytes"], "initial output_bytes", case_id);
            fix_option_symbol(&input, &mut output_bytes)
                .unwrap_or_else(|error| panic!("{case_id}: Rust failed: {error:?}"));

            let actual = json!({
                "schema_version": "cosmolkit-inchi-rdkit-cpp-v1",
                "rdkit_version": "2026.03.1",
                "source_sha256": "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f",
                "source_fragment_sha256": "5770d6ab655210bd67e5351ef3ebe39aadc0c8b8e96c1b6e2ec0935acb8844d9",
                "operation": "fixOptionSymbol",
                "case_id": case_id,
                "input": {
                    "input_bytes": input,
                    "output_bytes": bytes(
                        &official["input"]["output_bytes"],
                        "initial output_bytes",
                        case_id,
                    ),
                },
                "output": {
                    "input_bytes": bytes(
                        &official["output"]["input_bytes"],
                        "output input_bytes",
                        case_id,
                    ),
                    "output_bytes": output_bytes,
                },
            });
            assert_eq!(actual, official, "{case_id}: complete observable object");
        }

        assert_eq!(
            case_ids,
            BTreeSet::from([
                "all-nonzero-bytes".to_owned(),
                "empty".to_owned(),
                "exact-capacity".to_owned(),
                "first-nul-stops".to_owned(),
                "slash-and-ordinary".to_owned(),
            ])
        );
    }

    #[test]
    fn source_port__inchi__inchitomol__line_1254() {
        let mut failure_toolkit = RecordingToolkit::default();
        let failure_output = AdapterInchiStructureOutput {
            return_code: 2,
            atoms: Vec::new(),
            stereo0d: Vec::new(),
            message: None,
            log: Some(b"new-log".to_vec()),
        };
        let (failure, failure_engine, failure_values) = run_scripted(failure_output, &mut failure_toolkit, true, true);
        assert_eq!(failure.unwrap().molecule, None);
        assert_eq!(failure_engine.seen_inputs, vec![b"InChI=1S/test\0tail\0".to_vec()]);
        assert_eq!(failure_engine.free_count, 1);
        assert_eq!(failure_values.return_code, 2);
        assert_eq!(failure_values.message, b"old-message");
        assert_eq!(failure_values.log, b"new-log");
        assert_eq!(failure_values.aux_info, b"old-aux");
        assert!(failure_toolkit.calls.is_empty());

        let mut atoms: Vec<_> = (0..10).map(|_| source_atom(b"C")).collect();
        atoms[0].isotopic_mass = (ISOTOPIC_SHIFT_FLAG + 1) as i16;
        atoms[0].charge = -2;
        atoms[0].radical = 2;
        atoms[0].num_iso_H = [2, 2, 3, 4];
        atoms[1].radical = 3;
        atoms[2].radical = 4;
        atoms[2].num_iso_H[2] = 1;
        atoms[3].num_iso_H[3] = 1;
        for (bond_type, stereo) in [
            (0_i8, 0_i8),
            (1, 1),
            (2, -6),
            (3, 6),
            (4, -1),
            (1, 4),
            (2, 3),
            (1, -4),
            (1, 99),
        ] {
            let second = atoms[0].num_bonds as usize + 1;
            connect_source_atoms(&mut atoms, 0, second, bond_type, stereo, 0);
        }
        let mut mapping_toolkit = RecordingToolkit::default();
        let (mapping_result, mapping_engine, mapping_values) =
            run_scripted(scripted_output(atoms, Vec::new()), &mut mapping_toolkit, false, true);
        let mapping_result = mapping_result.unwrap();
        let mapping = mapping_result.molecule.unwrap();
        assert_eq!(mapping_engine.free_count, 1);
        assert_eq!(mapping_values.return_code, 0);
        assert_eq!(mapping_values.message, b"message");
        assert_eq!(mapping_values.log, b"log");
        assert_eq!(mapping.atoms[0].isotope, 13);
        assert_eq!(mapping.atoms[0].formal_charge, -2);
        assert_eq!(mapping.atoms[0].num_radical_electrons, 1);
        assert_eq!(mapping.atoms[0].num_explicit_hydrogens, 2);
        assert!(mapping.atoms[0].no_implicit);
        assert_eq!(mapping.atoms[1].num_radical_electrons, 2);
        assert_eq!(mapping.atoms[2].num_radical_electrons, 0);
        assert_eq!(mapping.atoms.len(), 14);
        assert_eq!(mapping.bonds.len(), 13);
        assert_eq!(
            mapping.bonds[..9].iter().map(|bond| bond.bond_type).collect::<Vec<_>>(),
            vec![
                BondType::Unspecified,
                BondType::Single,
                BondType::Double,
                BondType::Triple,
                BondType::Aromatic,
                BondType::Single,
                BondType::Double,
                BondType::Single,
                BondType::Single,
            ]
        );
        assert!(mapping.bonds[4].is_aromatic);
        assert_eq!(
            mapping.bonds[..9].iter().map(|bond| bond.direction).collect::<Vec<_>>(),
            vec![
                BondDirection::None,
                BondDirection::BeginWedge,
                BondDirection::BeginWedge,
                BondDirection::BeginDash,
                BondDirection::BeginDash,
                BondDirection::Unknown,
                BondDirection::EitherDouble,
                BondDirection::None,
                BondDirection::None,
            ]
        );
        assert!(
            mapping_result
                .diagnostics
                .iter()
                .any(|diagnostic| { diagnostic.message.contains("Ignore radical") })
        );
        assert!(
            mapping_result
                .diagnostics
                .iter()
                .any(|diagnostic| { diagnostic.message.contains("treated as aromatic") })
        );
        assert_eq!(mapping_toolkit.calls.last(), Some(&"assign_stereochemistry"));
        assert!(!mapping_toolkit.calls.contains(&"remove_hydrogens"));

        for (remove_hs, expected_call) in [(true, "remove_hydrogens"), (false, "sanitize_molecule")] {
            let mut toolkit = RecordingToolkit::default();
            let (result, engine, _) = run_scripted(
                scripted_output(vec![source_atom(b"C")], Vec::new()),
                &mut toolkit,
                true,
                remove_hs,
            );
            assert!(result.unwrap().molecule.is_some());
            assert_eq!(engine.free_count, 1);
            assert!(toolkit.calls.contains(&expected_call));
            assert_eq!(toolkit.calls.last(), Some(&"assign_stereochemistry"));
        }

        let mut illegal_atoms = vec![source_atom(b"C"), source_atom(b"C")];
        connect_source_atoms(&mut illegal_atoms, 0, 1, 5, 0, 0);
        let mut illegal_toolkit = RecordingToolkit::default();
        let (illegal, illegal_engine, _) = run_scripted(
            scripted_output(illegal_atoms, Vec::new()),
            &mut illegal_toolkit,
            true,
            true,
        );
        let illegal = illegal.unwrap();
        assert!(illegal.molecule.is_none());
        assert_eq!(illegal_engine.free_count, 1);
        assert_eq!(illegal.diagnostics.last().unwrap().level, AdapterDiagnosticLevel::Error);
        assert!(
            illegal_toolkit
                .calls
                .iter()
                .all(|call| { *call != "update_property_cache" && *call != "assign_stereochemistry" })
        );

        let stereo_graph = || {
            let mut atoms: Vec<_> = (0..6).map(|_| source_atom(b"C")).collect();
            connect_source_atoms(&mut atoms, 0, 1, 2, 0, 0);
            connect_source_atoms(&mut atoms, 0, 2, 1, 0, 0);
            connect_source_atoms(&mut atoms, 0, 3, 1, 0, 0);
            connect_source_atoms(&mut atoms, 1, 4, 1, 0, 0);
            connect_source_atoms(&mut atoms, 1, 5, 1, 0, 0);
            atoms
        };
        for (parity, original_left, expected_stereo) in [
            (1_i8, 2_i16, BondStereo::Z),
            (2, 2, BondStereo::E),
            (1, 3, BondStereo::E),
            (3, 2, BondStereo::Any),
        ] {
            let stereo = inchi_Stereo0D {
                neighbor: [original_left, 0, 1, 4],
                central_atom: 0,
                type_: tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
                parity,
            };
            let mut toolkit = RecordingToolkit {
                ranks: Some(vec![0, 0, 9, 3, 8, 2]),
                ..RecordingToolkit::default()
            };
            let (result, engine, _) = run_scripted(
                scripted_output(stereo_graph(), vec![stereo]),
                &mut toolkit,
                false,
                false,
            );
            let molecule = result.unwrap().molecule.unwrap();
            assert_eq!(engine.free_count, 1);
            assert_eq!(molecule.bonds[0].stereo, expected_stereo);
            assert_eq!(molecule.bonds[0].stereo_atoms, vec![2, 4]);
            assert!(
                molecule.bonds[1..]
                    .iter()
                    .any(|bond| { matches!(bond.direction, BondDirection::EndDownRight | BondDirection::EndUpRight) })
            );
        }

        for parity in [
            tagINCHIStereoParity0D_INCHI_PARITY_NONE as i8,
            tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED as i8,
        ] {
            let stereo = inchi_Stereo0D {
                neighbor: [2, 0, 1, 4],
                central_atom: 0,
                type_: tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
                parity,
            };
            let mut toolkit = RecordingToolkit::default();
            let (result, _, _) = run_scripted(
                scripted_output(stereo_graph(), vec![stereo]),
                &mut toolkit,
                false,
                false,
            );
            assert_eq!(result.unwrap().molecule.unwrap().bonds[0].stereo, BondStereo::None);
        }

        let absent_double = inchi_Stereo0D {
            neighbor: [0, 0, 1, 1],
            central_atom: 0,
            type_: tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
            parity: 1,
        };
        let mut absent_toolkit = RecordingToolkit::default();
        let (absent_result, _, _) = run_scripted(
            scripted_output(vec![source_atom(b"C"), source_atom(b"C")], vec![absent_double]),
            &mut absent_toolkit,
            false,
            false,
        );
        assert!(
            absent_result
                .unwrap()
                .diagnostics
                .iter()
                .any(|diagnostic| { diagnostic.message.contains("Extended double-bond") })
        );

        let mut insufficient_atoms = vec![source_atom(b"C"), source_atom(b"C"), source_atom(b"C")];
        connect_source_atoms(&mut insufficient_atoms, 0, 1, 2, 0, 0);
        connect_source_atoms(&mut insufficient_atoms, 0, 2, 1, 0, 0);
        let insufficient_stereo = inchi_Stereo0D {
            neighbor: [2, 0, 1, 2],
            central_atom: 0,
            type_: tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
            parity: 1,
        };
        let mut insufficient_toolkit = RecordingToolkit::default();
        let (insufficient, _, _) = run_scripted(
            scripted_output(insufficient_atoms, vec![insufficient_stereo]),
            &mut insufficient_toolkit,
            false,
            false,
        );
        assert!(
            insufficient
                .unwrap()
                .diagnostics
                .iter()
                .any(|diagnostic| { diagnostic.message.contains("without appropriate neighbors") })
        );

        let mut tetra_atoms: Vec<_> = (0..5).map(|_| source_atom(b"C")).collect();
        for neighbor in 1..5 {
            connect_source_atoms(&mut tetra_atoms, 0, neighbor, 1, 0, 0);
        }
        let tetra = inchi_Stereo0D {
            neighbor: [1, 2, 3, 4],
            central_atom: 0,
            type_: tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8,
            parity: tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8,
        };
        let mut tetra_toolkit = RecordingToolkit::default();
        let (tetra_result, _, _) = run_scripted(
            scripted_output(tetra_atoms, vec![tetra]),
            &mut tetra_toolkit,
            false,
            false,
        );
        assert_eq!(
            tetra_result.unwrap().molecule.unwrap().atoms[0].chiral_tag,
            ChiralTag::TetrahedralCcw
        );

        let mut tetra3_atoms: Vec<_> = (0..4).map(|_| source_atom(b"S")).collect();
        for neighbor in 1..4 {
            connect_source_atoms(&mut tetra3_atoms, 0, neighbor, 1, 0, 0);
        }
        let tetra3 = inchi_Stereo0D {
            neighbor: [0, 1, 2, 3],
            central_atom: 0,
            type_: tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8,
            parity: tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8,
        };
        let mut tetra3_toolkit = RecordingToolkit::default();
        let (tetra3_result, _, _) = run_scripted(
            scripted_output(tetra3_atoms, vec![tetra3]),
            &mut tetra3_toolkit,
            false,
            false,
        );
        assert_eq!(
            tetra3_result.unwrap().molecule.unwrap().atoms[0].chiral_tag,
            ChiralTag::TetrahedralCw
        );

        for (stereo_type, expected_text) in [
            (tagINCHIStereoType0D_INCHI_StereoType_None as i8, None),
            (tagINCHIStereoType0D_INCHI_StereoType_Allene as i8, Some("Allene-style")),
            (99, Some("Unrecognized stereo0D type")),
        ] {
            let stereo = inchi_Stereo0D {
                neighbor: [0; 4],
                central_atom: 0,
                type_: stereo_type,
                parity: 1,
            };
            let mut toolkit = RecordingToolkit::default();
            let (result, _, _) = run_scripted(
                scripted_output(vec![source_atom(b"C")], vec![stereo]),
                &mut toolkit,
                false,
                false,
            );
            let diagnostics = result.unwrap().diagnostics;
            match expected_text {
                Some(text) => assert!(diagnostics.iter().any(|entry| entry.message.contains(text))),
                None => assert!(diagnostics.is_empty()),
            }
        }

        let mut toolkit_error = RecordingToolkit {
            fail_on: Some("update_property_cache"),
            ..RecordingToolkit::default()
        };
        let (error, error_engine, _) = run_scripted(
            scripted_output(vec![source_atom(b"C")], Vec::new()),
            &mut toolkit_error,
            false,
            false,
        );
        assert!(matches!(error, Err(InchiToMolError::Toolkit(_))));
        assert_eq!(
            error_engine.free_count, 0,
            "source exception occurs before FreeStructFromINCHI"
        );

        let mut sanitize_error = RecordingToolkit {
            fail_on: Some("sanitize_molecule"),
            ..RecordingToolkit::default()
        };
        let (error, error_engine, _) = run_scripted(
            scripted_output(vec![source_atom(b"C")], Vec::new()),
            &mut sanitize_error,
            true,
            false,
        );
        assert!(matches!(error, Err(InchiToMolError::Toolkit(_))));
        assert_eq!(
            error_engine.free_count, 1,
            "sanitize exception occurs after FreeStructFromINCHI"
        );

        let mut bad_element = source_atom(b"C");
        bad_element.elname = [b'X' as i8, 0, 0, 0, 0, 0];
        let mut bad_element_toolkit = RecordingToolkit::default();
        let (error, error_engine, _) = run_scripted(
            scripted_output(vec![bad_element], Vec::new()),
            &mut bad_element_toolkit,
            false,
            false,
        );
        assert!(matches!(error, Err(InchiToMolError::Toolkit(_))));
        assert_eq!(error_engine.free_count, 0);

        assert_eq!(
            adapter_count_swaps_to_interconvert(&[0, 1, 2, 3], vec![1, 0, 2, 3]),
            Ok(1)
        );
        assert_eq!(
            adapter_count_swaps_to_interconvert(&[0, 1, 2, 3], vec![1, 2, 3, 0]),
            Ok(3)
        );
        assert_eq!(
            adapter_count_swaps_to_interconvert(&[0, 1, 2, 3], vec![1, 2, 0, 3]),
            Ok(2)
        );

        let mut source_heap = SourceHeap::default();
        source_heap.trace_source_allocations();
        {
            let mut source_engine = SourceInchiStructureEngine {
                heap: &mut source_heap,
                stdout: SourceMutPointer::null(),
                build: InchiBuildMetadata {
                    compiler: "gcc",
                    date: "Jan 01 1970",
                    time: "00:00:00",
                },
                clock_result: 1_000,
                pending_output: None,
            };
            let mut source_toolkit = RecordingToolkit::default();
            let mut source_values = ExtraInchiReturnValues::default();
            let source_result = inchi_to_mol(
                &mut source_engine,
                &mut source_toolkit,
                b"InChI=1S/CH4/h1H4",
                &mut source_values,
                false,
                false,
            )
            .unwrap();
            let source_molecule = source_result.molecule.unwrap();
            assert_eq!(source_values.return_code, tagRetValGetINCHI_inchi_Ret_OKAY);
            assert_eq!(source_molecule.atoms.len(), 1);
            assert_eq!(source_molecule.atoms[0].atomic_number, 6);
            assert_eq!(source_molecule.atoms[0].num_explicit_hydrogens, 4);
            assert!(source_engine.pending_output.is_none());
        }
        // The closed GetStructFromINCHI call graph retains one T_GROUP allocation on
        // this source path; its own focused tests require this exact allocation state.
        assert_eq!(source_heap.live_source_allocation_count(), 1);
        assert_eq!(
            source_heap.live_source_allocations_of::<crate::source_types::T_GROUP>(),
            1
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__inchitomol__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::{Value, json};

        fn integer(value: &Value, field: &str, case_id: &str) -> i64 {
            value
                .as_i64()
                .unwrap_or_else(|| panic!("{case_id}: {field} must be a signed integer"))
        }

        fn unsigned(value: &Value, field: &str, case_id: &str) -> u64 {
            value
                .as_u64()
                .unwrap_or_else(|| panic!("{case_id}: {field} must be an unsigned integer"))
        }

        fn parse_source_atom(value: &Value, case_id: &str) -> inchi_Atom {
            let mut atom = inchi_Atom::default();
            let element = value["element"]
                .as_str()
                .unwrap_or_else(|| panic!("{case_id}: atom element must be text"))
                .as_bytes();
            assert!(element.len() < atom.elname.len(), "{case_id}: element too long");
            for (destination, source) in atom.elname.iter_mut().zip(element.iter().copied()) {
                *destination = source as i8;
            }
            atom.elname[element.len()] = 0;

            let neighbors = value["neighbors"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: neighbors must be an array"));
            let bond_types = value["bond_types"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: bond_types must be an array"));
            let bond_stereo = value["bond_stereo"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: bond_stereo must be an array"));
            let count = usize::try_from(integer(&value["num_bonds"], "num_bonds", case_id))
                .expect("source bond count must fit usize");
            assert_eq!(neighbors.len(), count, "{case_id}: neighbor count");
            assert_eq!(bond_types.len(), count, "{case_id}: bond type count");
            assert_eq!(bond_stereo.len(), count, "{case_id}: bond stereo count");
            atom.num_bonds = i16::try_from(count).expect("source bond count must fit i16");
            for offset in 0..count {
                atom.neighbor[offset] =
                    i16::try_from(integer(&neighbors[offset], "neighbor", case_id)).expect("neighbor must fit i16");
                atom.bond_type[offset] =
                    i8::try_from(integer(&bond_types[offset], "bond_type", case_id)).expect("bond type must fit i8");
                atom.bond_stereo[offset] = i8::try_from(integer(&bond_stereo[offset], "bond_stereo", case_id))
                    .expect("bond stereo must fit i8");
            }
            let isotope_h = value["num_iso_h"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: num_iso_h must be an array"));
            assert_eq!(isotope_h.len(), 4, "{case_id}: isotope-H slot count");
            for (destination, source) in atom.num_iso_H.iter_mut().zip(isotope_h) {
                *destination =
                    i8::try_from(integer(source, "num_iso_h", case_id)).expect("isotope-H count must fit i8");
            }
            atom.isotopic_mass = i16::try_from(integer(&value["isotopic_mass"], "isotopic_mass", case_id))
                .expect("isotopic mass must fit i16");
            atom.radical = i8::try_from(integer(&value["radical"], "radical", case_id)).expect("radical must fit i8");
            atom.charge = i8::try_from(integer(&value["charge"], "charge", case_id)).expect("charge must fit i8");
            atom
        }

        fn parse_source_stereo(value: &Value, case_id: &str) -> inchi_Stereo0D {
            let neighbors = value["neighbors"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: stereo neighbors must be an array"));
            assert_eq!(neighbors.len(), 4, "{case_id}: stereo neighbor count");
            let mut result = inchi_Stereo0D::default();
            for (destination, source) in result.neighbor.iter_mut().zip(neighbors) {
                *destination =
                    i16::try_from(integer(source, "stereo neighbor", case_id)).expect("stereo neighbor must fit i16");
            }
            result.central_atom = i16::try_from(integer(&value["central_atom"], "central_atom", case_id))
                .expect("central atom must fit i16");
            result.type_ =
                i8::try_from(integer(&value["type"], "stereo type", case_id)).expect("stereo type must fit i8");
            result.parity =
                i8::try_from(integer(&value["parity"], "parity", case_id)).expect("stereo parity must fit i8");
            result
        }

        fn molecule_json(molecule: &AdapterMol) -> Value {
            let atoms = molecule
                .atoms
                .iter()
                .enumerate()
                .map(|(index, atom)| {
                    json!({
                        "index": index,
                        "atomic_number": atom.atomic_number,
                        "formal_charge": atom.formal_charge,
                        "num_explicit_hydrogens": atom.num_explicit_hydrogens,
                        "is_aromatic": atom.is_aromatic,
                        "isotope": atom.isotope,
                        "num_radical_electrons": atom.num_radical_electrons,
                        "no_implicit": atom.no_implicit,
                        "chiral_tag": atom.chiral_tag as u8,
                        "cip_rank": atom.cip_rank,
                    })
                })
                .collect::<Vec<_>>();
            let bonds = molecule
                .bonds
                .iter()
                .enumerate()
                .map(|(index, bond)| {
                    json!({
                        "index": index,
                        "begin_atom_index": bond.begin_atom_index,
                        "end_atom_index": bond.end_atom_index,
                        "bond_type": bond.bond_type as u8,
                        "direction": bond.direction as u8,
                        "is_aromatic": bond.is_aromatic,
                        "stereo": bond.stereo as u8,
                        "stereo_atoms": bond.stereo_atoms,
                    })
                })
                .collect::<Vec<_>>();
            json!({
                "atom_count": molecule.atoms.len(),
                "bond_count": molecule.bonds.len(),
                "atom_fields": atoms,
                "bond_fields": bonds,
            })
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--inchi-to-mol-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let mut case_ids = BTreeSet::new();

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "3765efd5f4f8a855a2c2231c1d62a81e7d5880da76804f47f400ed8f5a6464a1"
            );
            assert_eq!(official["operation"], "InchiToMol");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let input = &official["input"];
            let atoms = input["atoms"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: atoms must be an array"))
                .iter()
                .map(|atom| parse_source_atom(atom, case_id))
                .collect();
            let stereo0d = input["stereo0d"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: stereo0d must be an array"))
                .iter()
                .map(|stereo| parse_source_stereo(stereo, case_id))
                .collect();
            let source_output = AdapterInchiStructureOutput {
                return_code: i32::try_from(integer(&input["return_code"], "return_code", case_id))
                    .expect("return code must fit i32"),
                atoms,
                stereo0d,
                message: input["message"].as_str().map(|value| value.as_bytes().to_vec()),
                log: input["log"].as_str().map(|value| value.as_bytes().to_vec()),
            };
            let mut engine = ScriptedStructureEngine::new(source_output);
            if input["throw_on_get"] == true {
                engine.get_error = Some(SourceHeapError::UnsupportedSourceBehavior);
            }
            if input["throw_on_free"] == true {
                engine.free_error = Some(SourceHeapError::UnsupportedSourceBehavior);
            }
            let ranks = input["ranks"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: ranks must be an array"))
                .iter()
                .map(|rank| u32::try_from(unsigned(rank, "rank", case_id)).expect("rank must fit u32"))
                .collect::<Vec<_>>();
            let fail_on = match input["fail_on"].as_str() {
                None => None,
                Some("atomic_number") => Some("atomic_number"),
                Some("average_atomic_weight") => Some("average_atomic_weight"),
                Some("update_property_cache") => Some("update_property_cache"),
                Some("assign_atom_cip_ranks") => Some("assign_atom_cip_ranks"),
                Some("sanitize_molecule") => Some("sanitize_molecule"),
                Some("remove_hydrogens") => Some("remove_hydrogens"),
                Some("assign_stereochemistry") => Some("assign_stereochemistry"),
                Some(value) => panic!("{case_id}: unknown fail_on value {value}"),
            };
            let mut toolkit = RecordingToolkit {
                ranks: (!ranks.is_empty()).then_some(ranks),
                fail_on,
                ..RecordingToolkit::default()
            };
            let mut return_values = ExtraInchiReturnValues {
                return_code: -77,
                message: b"old-message".to_vec(),
                log: b"old-log".to_vec(),
                aux_info: b"old-aux".to_vec(),
            };
            let inchi = input["inchi"]
                .as_str()
                .unwrap_or_else(|| panic!("{case_id}: inchi must be text"));
            let rust = inchi_to_mol(
                &mut engine,
                &mut toolkit,
                inchi.as_bytes(),
                &mut return_values,
                input["sanitize"] == true,
                input["remove_hs"] == true,
            );

            let (status, exception, molecule, diagnostics) = match &rust {
                Ok(result) => (
                    "return",
                    Value::Null,
                    result.molecule.as_ref().map(molecule_json).unwrap_or(Value::Null),
                    result.diagnostics.as_slice(),
                ),
                Err(InchiToMolError::Toolkit(error)) => {
                    let (kind, message) = if fail_on.is_some() {
                        ("MolSanitizeException", error.message.clone())
                    } else {
                        (error.kind, error.message.clone())
                    };
                    (
                        "exception",
                        json!({"kind": kind, "message": message}),
                        Value::Null,
                        &[][..],
                    )
                }
                Err(InchiToMolError::Source(_)) => {
                    let message = if input["throw_on_get"] == true {
                        "GetStructFromINCHI"
                    } else {
                        "FreeStructFromINCHI"
                    };
                    (
                        "exception",
                        json!({"kind": "std::exception", "message": message}),
                        Value::Null,
                        &[][..],
                    )
                }
                Err(error) => panic!("{case_id}: unexpected Rust error {error:?}"),
            };
            let warning_log = diagnostics
                .iter()
                .filter(|diagnostic| diagnostic.level == AdapterDiagnosticLevel::Warning)
                .map(|diagnostic| format!("{}\n", diagnostic.message))
                .collect::<String>();
            let error_log = diagnostics
                .iter()
                .filter(|diagnostic| diagnostic.level == AdapterDiagnosticLevel::Error)
                .map(|diagnostic| format!("{}\n", diagnostic.message))
                .collect::<String>();
            let seen_inchi = engine
                .seen_inputs
                .first()
                .map(|bytes| {
                    let bytes = bytes.strip_suffix(&[0]).unwrap_or(bytes);
                    String::from_utf8(bytes.to_vec()).expect("seen InChI must be UTF-8")
                })
                .unwrap_or_default();
            let output = json!({
                "status": status,
                "exception": exception,
                "return_values": {
                    "return_code": return_values.return_code,
                    "message": String::from_utf8(return_values.message.clone()).expect("message must be UTF-8"),
                    "log": String::from_utf8(return_values.log.clone()).expect("log must be UTF-8"),
                    "aux_info": String::from_utf8(return_values.aux_info.clone()).expect("aux info must be UTF-8"),
                },
                "molecule": molecule,
                "warning_log": warning_log,
                "error_log": error_log,
                "calls": toolkit.calls,
                "get_count": engine.seen_inputs.len(),
                "free_count": engine.free_count,
                "outstanding_outputs": if input["throw_on_get"] == true || engine.free_count != 0 { 0 } else { 1 },
                "seen_inchi": seen_inchi,
                "seen_options": "",
            });
            assert_eq!(official["output"], output, "{case_id}");
        }

        let expected_case_ids = [
            "assign-bond-dirs-conflict",
            "assign-stereo-exception",
            "atom-bond-isotope-radical-mapping",
            "atomic-number-exception",
            "average-weight-exception",
            "cip-rank-exception",
            "double-e",
            "double-gcc-rank-boundary",
            "double-missing-neighbor",
            "double-switch-ez",
            "double-unknown",
            "double-z",
            "empty-success",
            "extended-double-missing-bond",
            "failure-return",
            "free-struct-exception",
            "get-struct-exception",
            "illegal-bond-type",
            "invalid-element-exception",
            "parity-none",
            "parity-undefined",
            "property-cache-exception",
            "remove-hs-exception",
            "sanitize-exception",
            "simple-no-sanitize",
            "simple-remove-hs",
            "simple-sanitize",
            "stereo-allene",
            "stereo-none",
            "stereo-unrecognized",
            "tetrahedral-four-neighbor",
            "tetrahedral-three-neighbor",
            "warning-return",
        ];
        assert_eq!(
            case_ids,
            expected_case_ids.iter().map(|case_id| (*case_id).to_owned()).collect()
        );
    }

    fn cleanup6_graph(root_atomic_number: i32, bridge_type: BondType) -> AdapterMol {
        graph(
            &[
                (6, -2),
                (6, 3),
                (root_atomic_number, 0),
                (6, -4),
                (6, 5),
                (7, -6),
                (6, 7),
            ],
            &[
                (1, 0, BondType::Single),
                (2, 1, BondType::Double),
                (2, 3, BondType::Double),
                (4, 3, bridge_type),
                (4, 5, BondType::Single),
                (0, 5, BondType::Double),
                (6, 2, BondType::Single),
            ],
        )
    }

    fn cleanup7_graph(root_atomic_number: i32, bridge_type: BondType) -> AdapterMol {
        graph(
            &[
                (6, -2),
                (6, 3),
                (root_atomic_number, 0),
                (7, -4),
                (6, 5),
                (8, -6),
                (6, 7),
                (8, 0),
            ],
            &[
                (1, 0, bridge_type),
                (2, 1, BondType::Double),
                (2, 3, BondType::Double),
                (3, 4, BondType::Single),
                (4, 5, BondType::Single),
                (5, 0, BondType::Single),
                (6, 2, BondType::Single),
                (4, 7, BondType::Double),
            ],
        )
    }

    fn cleanup8_graph(root_atomic_number: i32) -> AdapterMol {
        let mut molecule = graph(
            &[(6, -2), (7, 3), (6, -4), (7, 5), (7, -6), (root_atomic_number, 0)],
            &[
                (0, 1, BondType::Single),
                (1, 2, BondType::Double),
                (2, 3, BondType::Single),
                (3, 4, BondType::Double),
                (4, 0, BondType::Single),
                (5, 0, BondType::Double),
            ],
        );
        molecule.atoms[5].num_explicit_hydrogens = 3;
        molecule
    }

    fn cleanup9_graph(root_atomic_number: i32) -> AdapterMol {
        let mut molecule = graph(
            &[(6, -2), (7, 3), (7, -4), (6, 5), (6, -6), (root_atomic_number, 0)],
            &[
                (1, 0, BondType::Single),
                (1, 2, BondType::Double),
                (3, 2, BondType::Single),
                (3, 4, BondType::Double),
                (0, 4, BondType::Single),
                (5, 0, BondType::Double),
            ],
        );
        molecule.atoms[5].num_explicit_hydrogens = 3;
        molecule
    }

    fn cleanupa_path_graph(path_length: usize) -> AdapterMol {
        assert!(path_length > 0);
        let target_index = path_length as u32;
        let partner_index = target_index + 1;
        let mut atoms = vec![(7, 0)];
        atoms.extend((1..path_length).map(|_| (6, 0)));
        atoms.extend([(7, 0), (7, 0), (6, 0), (6, 0), (6, 0)]);
        let mut bonds = Vec::new();
        for index in 0..path_length {
            bonds.push((
                index as u32,
                index as u32 + 1,
                if index % 2 == 0 {
                    BondType::Double
                } else {
                    BondType::Single
                },
            ));
        }
        bonds.push((target_index, partner_index, BondType::Double));
        bonds.extend([
            (0, partner_index + 1, BondType::Single),
            (0, partner_index + 2, BondType::Single),
            (0, partner_index + 3, BondType::Single),
        ]);
        graph(&atoms, &bonds)
    }

    fn disjoint_union(first: &AdapterMol, second: &AdapterMol) -> AdapterMol {
        let offset = first.atoms.len() as u32;
        let mut atoms = first.atoms.clone();
        atoms.extend(second.atoms.iter().cloned());
        let mut bonds = first.bonds.clone();
        bonds.extend(second.bonds.iter().cloned().map(|mut bond| {
            bond.begin_atom_index += offset;
            bond.end_atom_index += offset;
            bond
        }));
        AdapterMol::from_graph(atoms, bonds)
    }

    fn directions(molecule: &AdapterMol) -> Vec<BondDirection> {
        molecule.bonds.iter().map(|bond| bond.direction).collect()
    }

    fn oracle_i32(value: &serde_json::Value, field: &str, case_id: &str) -> i32 {
        i32::try_from(
            value
                .as_i64()
                .unwrap_or_else(|| panic!("{case_id}: {field} must be i32")),
        )
        .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds i32"))
    }

    fn oracle_u32(value: &serde_json::Value, field: &str, case_id: &str) -> u32 {
        u32::try_from(
            value
                .as_u64()
                .unwrap_or_else(|| panic!("{case_id}: {field} must be u32")),
        )
        .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds u32"))
    }

    fn oracle_bond_type(value: &serde_json::Value, case_id: &str) -> BondType {
        match oracle_u32(value, "bond_type", case_id) {
            0 => BondType::Unspecified,
            1 => BondType::Single,
            2 => BondType::Double,
            3 => BondType::Triple,
            4 => BondType::Quadruple,
            5 => BondType::Quintuple,
            6 => BondType::Hextuple,
            7 => BondType::OneAndAHalf,
            8 => BondType::TwoAndAHalf,
            9 => BondType::ThreeAndAHalf,
            10 => BondType::FourAndAHalf,
            11 => BondType::FiveAndAHalf,
            12 => BondType::Aromatic,
            13 => BondType::Ionic,
            14 => BondType::Hydrogen,
            15 => BondType::ThreeCenter,
            16 => BondType::DativeOne,
            17 => BondType::Dative,
            18 => BondType::DativeL,
            19 => BondType::DativeR,
            20 => BondType::Other,
            21 => BondType::Zero,
            value => panic!("{case_id}: unknown BondType {value}"),
        }
    }

    fn oracle_direction(value: &serde_json::Value, case_id: &str) -> BondDirection {
        match oracle_u32(value, "direction", case_id) {
            0 => BondDirection::None,
            1 => BondDirection::BeginWedge,
            2 => BondDirection::BeginDash,
            3 => BondDirection::EndDownRight,
            4 => BondDirection::EndUpRight,
            5 => BondDirection::EitherDouble,
            6 => BondDirection::Unknown,
            value => panic!("{case_id}: unknown BondDirection {value}"),
        }
    }

    fn oracle_atoms(value: &serde_json::Value, case_id: &str) -> Vec<AdapterAtom> {
        value
            .as_array()
            .unwrap_or_else(|| panic!("{case_id}: atom_fields must be an array"))
            .iter()
            .enumerate()
            .map(|(index, field)| {
                assert_eq!(
                    oracle_u32(&field["index"], "atom index", case_id),
                    index as u32,
                    "{case_id}"
                );
                AdapterAtom {
                    atomic_number: oracle_i32(&field["atomic_number"], "atomic_number", case_id),
                    formal_charge: oracle_i32(&field["formal_charge"], "formal_charge", case_id),
                    num_explicit_hydrogens: oracle_u32(
                        &field["num_explicit_hydrogens"],
                        "num_explicit_hydrogens",
                        case_id,
                    ),
                    is_aromatic: field["is_aromatic"]
                        .as_bool()
                        .unwrap_or_else(|| panic!("{case_id}: is_aromatic must be bool")),
                    ..AdapterAtom::default()
                }
            })
            .collect()
    }

    fn oracle_bonds(value: &serde_json::Value, case_id: &str) -> Vec<AdapterBond> {
        value
            .as_array()
            .unwrap_or_else(|| panic!("{case_id}: bond_fields must be an array"))
            .iter()
            .enumerate()
            .map(|(index, field)| {
                assert_eq!(
                    oracle_u32(&field["index"], "bond index", case_id),
                    index as u32,
                    "{case_id}"
                );
                AdapterBond {
                    begin_atom_index: oracle_u32(&field["begin_atom_index"], "begin_atom_index", case_id),
                    end_atom_index: oracle_u32(&field["end_atom_index"], "end_atom_index", case_id),
                    bond_type: oracle_bond_type(&field["bond_type"], case_id),
                    direction: oracle_direction(&field["direction"], case_id),
                    ..AdapterBond::default()
                }
            })
            .collect()
    }

    fn assert_valence_cleanup_oracle(
        selector: &str,
        source_fragment_sha256: &str,
        operation: &str,
        expected_case_ids: &[&str],
        rust_operation: fn(&mut AdapterMol, u32) -> Result<bool, AdapterValenceError>,
    ) {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg(selector)
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let mut case_ids = BTreeSet::new();

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(official["source_fragment_sha256"], source_fragment_sha256);
            assert_eq!(official["operation"], operation);
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let mut molecule = AdapterMol::from_graph(
                oracle_atoms(&official["input"]["atom_fields"], case_id),
                oracle_bonds(&official["input"]["bond_fields"], case_id),
            );
            let molecule_before = molecule.clone();
            let atom_index = oracle_u32(&official["input"]["atom_index"], "atom_index", case_id);
            let rust = rust_operation(&mut molecule, atom_index);

            match official["output"]["status"].as_str() {
                Some("return") => {
                    assert_eq!(
                        rust,
                        Ok(official["output"]["result"]
                            .as_bool()
                            .unwrap_or_else(|| panic!("{case_id}: result must be bool"))),
                        "{case_id}"
                    );
                    assert!(official["output"]["exception"].is_null(), "{case_id}");
                    assert_eq!(official["output"]["diagnostics"], Value::Array(vec![]));
                }
                Some("exception") => {
                    let exception = &official["output"]["exception"];
                    let expected = AdapterValenceError {
                        kind: "ValueErrorException",
                        message: "Bad bond type",
                        bond_index: oracle_u32(&exception["bond_index"], "bond_index", case_id),
                        bond_type: oracle_bond_type(&exception["bond_type"], case_id),
                    };
                    assert_eq!(exception["kind"], expected.kind, "{case_id}");
                    assert_eq!(exception["message"], expected.message, "{case_id}");
                    assert_eq!(rust, Err(expected.clone()), "{case_id}");
                    assert_eq!(official["output"]["result"], Value::Null, "{case_id}");
                    let diagnostics = official["output"]["diagnostics"]
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: diagnostics must be array"));
                    assert_eq!(diagnostics.len(), 1, "{case_id}");
                    assert_eq!(diagnostics[0]["kind"], expected.kind, "{case_id}");
                    assert_eq!(diagnostics[0]["message"], expected.message, "{case_id}");
                    assert_eq!(
                        oracle_u32(&diagnostics[0]["bond_index"], "bond_index", case_id),
                        expected.bond_index,
                        "{case_id}"
                    );
                    assert_eq!(
                        oracle_bond_type(&diagnostics[0]["bond_type"], case_id),
                        expected.bond_type,
                        "{case_id}"
                    );
                }
                status => panic!("{case_id}: unexpected status {status:?}"),
            }

            assert_eq!(
                official["output"]["graph_unchanged"].as_bool(),
                Some(molecule == molecule_before),
                "{case_id}"
            );
            assert_eq!(
                oracle_u32(&official["output"]["atom_count"], "atom_count", case_id),
                molecule.atoms.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                oracle_u32(&official["output"]["bond_count"], "bond_count", case_id),
                molecule.bonds.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                oracle_atoms(&official["output"]["atom_fields"], case_id),
                molecule.atoms,
                "{case_id}"
            );
            assert_eq!(
                oracle_bonds(&official["output"]["bond_fields"], case_id),
                molecule.bonds,
                "{case_id}"
            );
            assert_eq!(
                official["output"]["stereo_fields"]["bond_directions"],
                Value::Array(
                    molecule
                        .bonds
                        .iter()
                        .map(|bond| Value::from(bond.direction as u8))
                        .collect()
                ),
                "{case_id}"
            );
            assert_eq!(official["output"]["properties"], Value::Array(vec![]));
        }

        assert_eq!(
            case_ids,
            expected_case_ids.iter().map(|case_id| (*case_id).to_owned()).collect()
        );
    }

    fn assert_clean_up_oracle() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--clean-up-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let mut case_ids = BTreeSet::new();

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "25458d16f3e2888aed0a60ba3b1c8f0a4e435a912215af7da93fbc491b34b371"
            );
            assert_eq!(official["operation"], "cleanUp");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let mut molecule = AdapterMol::from_graph(
                oracle_atoms(&official["input"]["atom_fields"], case_id),
                oracle_bonds(&official["input"]["bond_fields"], case_id),
            );
            let molecule_before = molecule.clone();
            let rust = clean_up(&mut molecule);

            assert_eq!(official["output"]["result"], Value::Null, "{case_id}");
            match official["output"]["status"].as_str() {
                Some("return") => {
                    assert_eq!(rust, Ok(()), "{case_id}");
                    assert!(official["output"]["exception"].is_null(), "{case_id}");
                    assert_eq!(official["output"]["diagnostics"], Value::Array(vec![]), "{case_id}");
                }
                Some("exception") => {
                    let exception = &official["output"]["exception"];
                    let expected = AdapterValenceError {
                        kind: "ValueErrorException",
                        message: "Bad bond type",
                        bond_index: oracle_u32(&exception["bond_index"], "bond_index", case_id),
                        bond_type: oracle_bond_type(&exception["bond_type"], case_id),
                    };
                    assert_eq!(exception["kind"], expected.kind, "{case_id}");
                    assert_eq!(exception["message"], expected.message, "{case_id}");
                    assert_eq!(rust, Err(AdapterCleanup5Error::Valence(expected.clone())), "{case_id}");
                    let diagnostics = official["output"]["diagnostics"]
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: diagnostics must be array"));
                    assert_eq!(diagnostics.len(), 1, "{case_id}");
                    assert_eq!(diagnostics[0]["kind"], expected.kind, "{case_id}");
                    assert_eq!(diagnostics[0]["message"], expected.message, "{case_id}");
                    assert_eq!(
                        oracle_u32(&diagnostics[0]["bond_index"], "bond_index", case_id),
                        expected.bond_index,
                        "{case_id}"
                    );
                    assert_eq!(
                        oracle_bond_type(&diagnostics[0]["bond_type"], case_id),
                        expected.bond_type,
                        "{case_id}"
                    );
                }
                status => panic!("{case_id}: unexpected status {status:?}"),
            }

            assert_eq!(
                official["output"]["graph_unchanged"].as_bool(),
                Some(molecule == molecule_before),
                "{case_id}"
            );
            assert_eq!(
                oracle_u32(&official["output"]["atom_count"], "atom_count", case_id),
                molecule.atoms.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                oracle_u32(&official["output"]["bond_count"], "bond_count", case_id),
                molecule.bonds.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                oracle_atoms(&official["output"]["atom_fields"], case_id),
                molecule.atoms,
                "{case_id}"
            );
            assert_eq!(
                oracle_bonds(&official["output"]["bond_fields"], case_id),
                molecule.bonds,
                "{case_id}"
            );
            assert_eq!(
                official["output"]["stereo_fields"]["bond_directions"],
                Value::Array(
                    molecule
                        .bonds
                        .iter()
                        .map(|bond| Value::from(bond.direction as u8))
                        .collect()
                ),
                "{case_id}"
            );
            assert_eq!(official["output"]["properties"], Value::Array(vec![]));
        }

        assert_eq!(
            case_ids,
            [
                "empty",
                "ignored-aromatic-carbon",
                "entry-valence-exception",
                "valence4-cleanup1",
                "valence4-cleanup2",
                "valence4-no-match",
                "charged-nitrogen",
                "valence5-all-miss",
                "valence5-cleanup1",
                "valence5-cleanup2",
                "valence5-cleanup3",
                "valence5-cleanup4",
                "valence5-cleanup5-oxygen",
                "valence5-cleanup5-sulfur",
                "valence5-cleanup5-fluorine",
                "valence5-cleanup5-chlorine",
                "valence5-cleanupb",
                "ordered-partial-error",
                "chlorine8-empty-neighbors",
                "chlorine8-all-oxygen",
                "chlorine5-dispatch-helper-mismatch",
                "chlorine3-cleanup",
                "sulfur7-cleanup1",
                "sulfur7-cleanup2",
                "sulfur7-cleanup3",
                "sulfur7-all-miss",
                "sulfur7-long-path-cleanup",
                "sulfur8-source-defined-no-op",
                "bromine-selenium",
                "bromine-charged-miss",
                "bromine-degree-miss",
                "bromine-element-miss",
                "valence5-cleanup6",
                "valence5-cleanup7",
                "valence5-cleanup8",
                "valence5-cleanup9",
                "valence5-cleanupa",
            ]
            .into_iter()
            .map(str::to_owned)
            .collect()
        );
    }

    fn assert_valence_cleanup5_oracle() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--valence5n-cleanup5-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );

        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let expected_case_ids = BTreeSet::from([
            "precondition-invalid".to_owned(),
            "no-target-oxygen".to_owned(),
            "no-target-sulfur".to_owned(),
            "no-target-fluorine".to_owned(),
            "no-target-chlorine".to_owned(),
            "only-uncharged-alternating".to_owned(),
            "only-charged-reversed".to_owned(),
            "both-targets-neutral-path-only".to_owned(),
            "depth-seven-inclusive".to_owned(),
            "over-depth-limit".to_owned(),
            "only-uncharged-valence-exception".to_owned(),
            "only-charged-valence-exception".to_owned(),
            "both-charged-first-exception".to_owned(),
            "both-uncharged-second-exception".to_owned(),
        ]);
        let mut case_ids = BTreeSet::new();

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "572b129f3625cbffe00811214403e660abdcf283c582ff90287011755698b73e"
            );
            assert_eq!(official["operation"], "_Valence5NCleanUp5");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let mut molecule = AdapterMol::from_graph(
                oracle_atoms(&official["input"]["atom_fields"], case_id),
                oracle_bonds(&official["input"]["bond_fields"], case_id),
            );
            let molecule_before = molecule.clone();
            let atom_index = oracle_u32(&official["input"]["atom_index"], "atom_index", case_id);
            let atomic_number = oracle_i32(&official["input"]["atomic_number"], "atomic_number", case_id);
            let rust = valence5n_cleanup5(&mut molecule, atom_index, atomic_number);

            match official["output"]["status"].as_str() {
                Some("return") => {
                    assert_eq!(
                        rust,
                        Ok(official["output"]["result"]
                            .as_bool()
                            .unwrap_or_else(|| panic!("{case_id}: result must be bool"))),
                        "{case_id}"
                    );
                    assert!(official["output"]["exception"].is_null(), "{case_id}");
                    assert_eq!(official["output"]["diagnostics"], Value::Array(vec![]));
                }
                Some("exception") => {
                    let exception = &official["output"]["exception"];
                    let kind = exception["kind"]
                        .as_str()
                        .unwrap_or_else(|| panic!("{case_id}: exception kind must be text"));
                    let expected = if kind == "Pre-condition Violation" {
                        assert_eq!(
                            exception["expression"],
                            "atomicNum == 8 || atomicNum == 16 || atomicNum == 9 || atomicNum == 17",
                            "{case_id}"
                        );
                        assert_eq!(
                            exception["message"], "this cleanup looks for O or S or Cl or F",
                            "{case_id}"
                        );
                        AdapterCleanup5Error::Precondition {
                            kind: "Pre-condition Violation",
                            expression: "atomicNum == 8 || atomicNum == 16 || atomicNum == 9 || atomicNum == 17",
                            message: "this cleanup looks for O or S or Cl or F",
                        }
                    } else {
                        assert_eq!(kind, "ValueErrorException", "{case_id}");
                        assert_eq!(exception["message"], "Bad bond type", "{case_id}");
                        AdapterCleanup5Error::Valence(AdapterValenceError {
                            kind: "ValueErrorException",
                            message: "Bad bond type",
                            bond_index: oracle_u32(&exception["bond_index"], "bond_index", case_id),
                            bond_type: oracle_bond_type(&exception["bond_type"], case_id),
                        })
                    };
                    assert_eq!(rust, Err(expected), "{case_id}");
                    assert_eq!(official["output"]["result"], Value::Null, "{case_id}");
                    let diagnostics = official["output"]["diagnostics"]
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: diagnostics must be array"));
                    assert_eq!(diagnostics.len(), 1, "{case_id}");
                    assert_eq!(diagnostics[0], *exception, "{case_id}");
                }
                status => panic!("{case_id}: unexpected status {status:?}"),
            }

            assert_eq!(
                official["output"]["graph_unchanged"].as_bool(),
                Some(molecule == molecule_before),
                "{case_id}"
            );
            assert_eq!(
                oracle_u32(&official["output"]["atom_count"], "atom_count", case_id),
                molecule.atoms.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                oracle_u32(&official["output"]["bond_count"], "bond_count", case_id),
                molecule.bonds.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                oracle_atoms(&official["output"]["atom_fields"], case_id),
                molecule.atoms,
                "{case_id}"
            );
            assert_eq!(
                oracle_bonds(&official["output"]["bond_fields"], case_id),
                molecule.bonds,
                "{case_id}"
            );
            assert_eq!(
                official["output"]["stereo_fields"]["bond_directions"],
                Value::Array(
                    molecule
                        .bonds
                        .iter()
                        .map(|bond| Value::from(bond.direction as u8))
                        .collect()
                ),
                "{case_id}"
            );
            assert_eq!(official["output"]["properties"], Value::Array(vec![]));
        }

        assert_eq!(case_ids, expected_case_ids);
    }

    #[test]
    fn source_port__inchi__assignbonddirs__line_91() {
        use BondDirection::{BeginDash, BeginWedge, EitherDouble, EndDownRight, EndUpRight, None, Unknown};

        let preserved = [
            None,
            BeginWedge,
            BeginDash,
            EndDownRight,
            EndUpRight,
            EitherDouble,
            Unknown,
        ];
        let mut mol = molecule(&preserved);
        assert_eq!(assign_bond_dirs(&mut mol, &[], &[]), Ok(true));
        assert_eq!(directions(&mol), preserved);

        let mut mol = molecule(&[None; 3]);
        assert_eq!(assign_bond_dirs(&mut mol, &[(2, 0), (0, 1)], &[]), Ok(true));
        assert_eq!(directions(&mol), [EndUpRight; 3]);

        let mut mol = molecule(&[None; 3]);
        assert_eq!(assign_bond_dirs(&mut mol, &[], &[(2, 0), (0, 1)]), Ok(true));
        assert_eq!(directions(&mol), [EndUpRight, EndDownRight, EndDownRight]);

        let mut mol = molecule(&[None; 6]);
        assert_eq!(assign_bond_dirs(&mut mol, &[(5, 4), (3, 2)], &[(1, 0)]), Ok(true));
        assert_eq!(
            directions(&mol),
            [EndUpRight, EndDownRight, EndUpRight, EndUpRight, EndUpRight, EndUpRight,]
        );

        let mut mol = molecule(&[None; 2]);
        assert_eq!(assign_bond_dirs(&mut mol, &[(0, 1), (1, 0), (0, 1)], &[]), Ok(true));
        assert_eq!(directions(&mol), [EndUpRight, EndUpRight]);

        let mut mol = molecule(&[None; 3]);
        assert_eq!(assign_bond_dirs(&mut mol, &[(1, 1)], &[]), Ok(true));
        assert_eq!(directions(&mol), [None, EndUpRight, None]);

        let mut mol = molecule(&[None; 2]);
        assert_eq!(assign_bond_dirs(&mut mol, &[(0, 1)], &[(0, 1)]), Ok(false));
        assert_eq!(directions(&mol), [EndUpRight, EndUpRight]);

        let mut mol = molecule(&[None; 3]);
        assert_eq!(assign_bond_dirs(&mut mol, &[(0, 1), (1, 2)], &[(0, 2)]), Ok(false));
        assert_eq!(directions(&mol), [EndUpRight, EndUpRight, EndDownRight]);

        let mut mol = molecule(&[BeginWedge]);
        assert_eq!(assign_bond_dirs(&mut mol, &[(0, 0)], &[]), Ok(false));
        assert_eq!(directions(&mol), [BeginWedge]);

        let mut mol = molecule(&[None, EndDownRight]);
        assert_eq!(assign_bond_dirs(&mut mol, &[(0, 1)], &[]), Ok(false));
        assert_eq!(directions(&mol), [EndUpRight, EndDownRight]);

        let mut mol = molecule(&[None, EndDownRight]);
        assert_eq!(assign_bond_dirs(&mut mol, &[], &[(0, 1)]), Ok(false));
        assert_eq!(directions(&mol), [EndUpRight, EndDownRight]);

        let mut mol = molecule(&[None]);
        assert_eq!(
            assign_bond_dirs(&mut mol, &[(-1, -1)], &[]),
            Err(BondIndexRangeError {
                kind: "Range Error",
                expression: "idx",
                index: u32::MAX,
                upper_bound: 1,
            })
        );
        assert_eq!(directions(&mol), [None]);

        let mut mol = molecule(&[None; 2]);
        assert_eq!(
            assign_bond_dirs(&mut mol, &[(0, 3)], &[]),
            Err(BondIndexRangeError {
                kind: "Range Error",
                expression: "idx",
                index: 3,
                upper_bound: 2,
            })
        );
        assert_eq!(directions(&mol), [EndUpRight, None]);
    }

    #[test]
    fn source_port__inchi__findalternatingbonds__line_195() {
        use BondType::{Double, Single, Triple};

        let isolated = graph(&[(6, 0)], &[]);
        let mut path = vec![7, 8];
        let mut visited = BTreeSet::from([7, 8]);
        assert_eq!(
            find_alternating_bonds(&isolated, 0, 7, 0, Single, Single, 0, 0, None, &mut path, &mut visited,),
            None
        );
        assert!(path.is_empty());
        assert_eq!(visited, BTreeSet::from([0]));

        let direct = graph(&[(6, 0), (7, 0)], &[(0, 1, Single)]);
        path = vec![0, 0];
        visited = BTreeSet::from([0]);
        assert_eq!(
            find_alternating_bonds(&direct, 1, 7, 0, Double, Single, 1, 0, Some(0), &mut path, &mut visited,),
            Some(1)
        );
        assert_eq!(path, [0]);
        assert_eq!(visited, BTreeSet::from([0, 1]));

        path = vec![0];
        visited = BTreeSet::from([0]);
        assert_eq!(
            find_alternating_bonds(&direct, 1, 7, 0, Double, Single, 1, 9, Some(0), &mut path, &mut visited,),
            None
        );
        assert_eq!(path, [0]);

        let alternating = graph(
            &[(6, 0), (6, 0), (6, 0), (7, 0)],
            &[(0, 1, Single), (1, 2, Double), (2, 3, Single)],
        );
        path.clear();
        visited.clear();
        assert_eq!(
            find_alternating_bonds(
                &alternating,
                0,
                7,
                0,
                Single,
                Single,
                0,
                4,
                None,
                &mut path,
                &mut visited,
            ),
            Some(3)
        );
        assert_eq!(path, [2, 1, 0]);
        assert_eq!(visited, BTreeSet::from([0, 1, 2, 3]));

        let double_start = graph(&[(6, 0), (6, 0), (7, 0)], &[(0, 1, Double), (1, 2, Single)]);
        path.clear();
        visited.clear();
        assert_eq!(
            find_alternating_bonds(
                &double_start,
                0,
                7,
                0,
                Double,
                Single,
                0,
                2,
                None,
                &mut path,
                &mut visited,
            ),
            Some(2)
        );
        assert_eq!(path, [1, 0]);

        let triple_then_single = graph(&[(7, 0), (6, -1), (7, -1)], &[(0, 1, Triple), (1, 2, Single)]);
        path.clear();
        visited.clear();
        assert_eq!(
            find_alternating_bonds(
                &triple_then_single,
                0,
                7,
                -1,
                Triple,
                Single,
                0,
                2,
                None,
                &mut path,
                &mut visited,
            ),
            Some(2)
        );
        assert_eq!(path, [1, 0]);

        let special_ending = graph(&[(6, 0), (7, 0)], &[(0, 1, Triple)]);
        path.clear();
        visited.clear();
        assert_eq!(
            find_alternating_bonds(
                &special_ending,
                0,
                7,
                0,
                Double,
                Triple,
                0,
                3,
                None,
                &mut path,
                &mut visited,
            ),
            Some(1)
        );
        assert_eq!(path, [0]);

        let special_ending_miss = graph(&[(6, 0), (6, 0), (7, 0)], &[(0, 1, Triple), (1, 2, Single)]);
        path.clear();
        visited.clear();
        assert_eq!(
            find_alternating_bonds(
                &special_ending_miss,
                0,
                7,
                0,
                Double,
                Triple,
                0,
                3,
                None,
                &mut path,
                &mut visited,
            ),
            None
        );
        assert!(path.is_empty());
        assert_eq!(visited, BTreeSet::from([0, 1]));

        let disabled_special_ending = graph(&[(6, 0), (7, 0)], &[(0, 1, Double)]);
        path.clear();
        visited.clear();
        assert_eq!(
            find_alternating_bonds(
                &disabled_special_ending,
                0,
                7,
                0,
                Single,
                Double,
                0,
                2,
                None,
                &mut path,
                &mut visited,
            ),
            None
        );
        assert_eq!(visited, BTreeSet::from([0]));

        let cycle = graph(
            &[(6, 0), (6, 0), (6, 0)],
            &[(0, 1, Single), (1, 2, Double), (2, 0, Single)],
        );
        path.clear();
        visited.clear();
        assert_eq!(
            find_alternating_bonds(&cycle, 0, 7, 0, Single, Single, 0, 9, None, &mut path, &mut visited,),
            None
        );
        assert!(path.is_empty());
        assert_eq!(visited, BTreeSet::from([0, 1, 2]));

        let shorter_later = graph(
            &[(6, 0), (6, 0), (6, 0), (7, 0), (7, 0)],
            &[(0, 1, Single), (1, 2, Double), (2, 3, Single), (0, 4, Single)],
        );
        path.clear();
        visited.clear();
        assert_eq!(
            find_alternating_bonds(
                &shorter_later,
                0,
                7,
                0,
                Single,
                Single,
                0,
                4,
                None,
                &mut path,
                &mut visited,
            ),
            Some(4)
        );
        assert_eq!(path, [3]);
        assert_eq!(visited, BTreeSet::from([0, 1, 2, 3, 4]));

        let equal_paths = graph(&[(6, 0), (7, 0), (7, 0)], &[(0, 1, Single), (0, 2, Single)]);
        path.clear();
        visited.clear();
        assert_eq!(
            find_alternating_bonds(
                &equal_paths,
                0,
                7,
                0,
                Single,
                Single,
                0,
                1,
                None,
                &mut path,
                &mut visited,
            ),
            Some(1)
        );
        assert_eq!(path, [0]);
        assert_eq!(visited, BTreeSet::from([0, 1, 2]));

        let charge_miss = graph(&[(6, 0), (7, 1)], &[(0, 1, Single)]);
        path.clear();
        visited.clear();
        assert_eq!(
            find_alternating_bonds(
                &charge_miss,
                0,
                7,
                0,
                Single,
                Single,
                0,
                1,
                None,
                &mut path,
                &mut visited,
            ),
            None
        );
        assert!(path.is_empty());
    }

    #[test]
    fn source_port__inchi__getnumdoublebondednegativelychargedneighboringsi__line_306() {
        use BondType::{Double, Single, Triple};

        let isolated = graph(&[(6, 0)], &[]);
        let isolated_before = isolated.clone();
        assert_eq!(get_num_double_bonded_negatively_charged_neighboring_si(&isolated, 0), 0);
        assert_eq!(isolated, isolated_before);

        let molecule = graph(
            &[(6, 0), (14, -1), (14, 0), (13, -1), (14, -1), (14, -1), (14, -1)],
            &[
                (0, 1, Double),
                (0, 2, Double),
                (0, 3, Double),
                (0, 4, Single),
                (5, 0, Double),
                (0, 6, Triple),
            ],
        );
        let molecule_before = molecule.clone();
        assert_eq!(get_num_double_bonded_negatively_charged_neighboring_si(&molecule, 0), 2);
        assert_eq!(get_num_double_bonded_negatively_charged_neighboring_si(&molecule, 1), 0);
        assert_eq!(molecule, molecule_before);
    }

    #[test]
    fn source_port__inchi__valence4ncleanup1__line_324() {
        use BondDirection::BeginDash;
        use BondType::{Double, Single, ThreeCenter, Triple};

        let wrong_element = graph(&[(6, -1), (6, 0)], &[(0, 1, ThreeCenter)]);
        let mut wrong_element_actual = wrong_element.clone();
        assert_eq!(valence4n_cleanup1(&mut wrong_element_actual, 0), Ok(false));
        assert_eq!(wrong_element_actual, wrong_element);

        let wrong_charge = graph(&[(7, 0), (6, 0)], &[(0, 1, ThreeCenter)]);
        let mut wrong_charge_actual = wrong_charge.clone();
        assert_eq!(valence4n_cleanup1(&mut wrong_charge_actual, 0), Ok(false));
        assert_eq!(wrong_charge_actual, wrong_charge);

        let wrong_valence = graph(&[(7, -1), (6, 0)], &[(0, 1, Triple)]);
        let mut wrong_valence_actual = wrong_valence.clone();
        assert_eq!(valence4n_cleanup1(&mut wrong_valence_actual, 0), Ok(false));
        assert_eq!(wrong_valence_actual, wrong_valence);

        let unsupported = graph(&[(7, -1), (6, 0)], &[(0, 1, ThreeCenter)]);
        let mut unsupported_actual = unsupported.clone();
        assert_eq!(
            valence4n_cleanup1(&mut unsupported_actual, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(unsupported_actual, unsupported);

        let mut explicit_h = graph(&[(7, -1), (6, 0)], &[(0, 1, Triple)]);
        explicit_h.atoms[0].num_explicit_hydrogens = 1;
        let explicit_h_before = explicit_h.clone();
        assert_eq!(valence4n_cleanup1(&mut explicit_h, 0), Ok(false));
        assert_eq!(explicit_h, explicit_h_before);

        let mut no_match = graph(&[(7, -1), (7, 0), (7, 0)], &[(0, 1, Double), (2, 0, Double)]);
        let no_match_before = no_match.clone();
        assert_eq!(valence4n_cleanup1(&mut no_match, 0), Ok(false));
        assert_eq!(no_match, no_match_before);

        let mut unique = graph(
            &[(6, 4), (7, 2), (7, -1), (7, -3), (7, 1), (8, 0)],
            &[
                (1, 0, Single),
                (2, 1, Double),
                (3, 2, Double),
                (4, 3, Single),
                (0, 4, Double),
                (0, 5, Single),
            ],
        );
        unique.bonds[0].direction = BeginDash;
        let unique_atoms_before = unique.atoms.clone();
        assert_eq!(valence4n_cleanup1(&mut unique, 2), Ok(true));
        assert_eq!(unique.atoms, unique_atoms_before);
        assert_eq!(
            unique.bonds.iter().map(|bond| bond.bond_type).collect::<Vec<_>>(),
            vec![Double, Single, Single, Double, Single, Single]
        );
        assert_eq!(unique.bonds[0].direction, BeginDash);

        let mut multiple = graph(
            &[(6, 0), (7, 0), (7, -1), (7, 0), (7, 0), (6, 0), (7, 0)],
            &[
                (0, 1, Single),
                (1, 2, Double),
                (2, 3, Double),
                (3, 4, Single),
                (4, 0, Double),
                (1, 5, Single),
                (5, 6, Double),
                (6, 3, Single),
            ],
        );
        let multiple_before = multiple.clone();
        assert_eq!(valence4n_cleanup1(&mut multiple, 2), Ok(false));
        assert_eq!(multiple, multiple_before);

        let mut remote = graph(
            &[(7, -1), (8, 0), (8, 0), (6, 0), (7, 0), (50, 0), (7, 0), (7, 0)],
            &[
                (0, 1, Double),
                (2, 0, Double),
                (3, 4, Single),
                (4, 5, Double),
                (5, 6, Double),
                (6, 7, Single),
                (7, 3, Double),
            ],
        );
        let remote_atoms_before = remote.atoms.clone();
        assert_eq!(valence4n_cleanup1(&mut remote, 0), Ok(true));
        assert_eq!(remote.atoms, remote_atoms_before);
        assert_eq!(
            remote.bonds.iter().map(|bond| bond.bond_type).collect::<Vec<_>>(),
            vec![Double, Double, Double, Single, Single, Double, Single]
        );
    }

    #[test]
    fn source_port__inchi__valence4ncleanup2__line_376() {
        use BondDirection::{BeginWedge, EndDownRight};
        use BondType::{Double, Single, Triple};

        let no_target = graph(
            &[(7, -1), (6, 0), (7, 1), (7, 0), (7, 0)],
            &[(0, 1, Double), (0, 2, Double), (0, 3, Single), (0, 4, Triple)],
        );
        let mut no_target_actual = no_target.clone();
        assert!(!valence4n_cleanup2(&mut no_target_actual, 0));
        assert_eq!(no_target_actual, no_target);

        let mut unique = graph(&[(7, -1), (7, 0), (8, -2)], &[(1, 0, Double), (0, 2, Single)]);
        unique.atoms[0].num_explicit_hydrogens = 2;
        unique.atoms[1].num_explicit_hydrogens = 1;
        unique.bonds[0].direction = BeginWedge;
        unique.bonds[1].direction = EndDownRight;
        assert!(valence4n_cleanup2(&mut unique, 0));
        assert_eq!(unique.atoms[0].formal_charge, 0);
        assert_eq!(unique.atoms[1].formal_charge, -1);
        assert_eq!(unique.atoms[2].formal_charge, -2);
        assert_eq!(unique.atoms[0].num_explicit_hydrogens, 2);
        assert_eq!(unique.atoms[1].num_explicit_hydrogens, 1);
        assert_eq!(unique.bonds[0].bond_type, Single);
        assert_eq!(unique.bonds[0].direction, BeginWedge);
        assert_eq!(unique.bonds[1].bond_type, Single);
        assert_eq!(unique.bonds[1].direction, EndDownRight);

        let mut multiple = graph(&[(7, -1), (7, 0), (7, 0)], &[(0, 2, Double), (0, 1, Double)]);
        assert!(valence4n_cleanup2(&mut multiple, 0));
        assert_eq!(multiple.atoms[0].formal_charge, 0);
        assert_eq!(multiple.atoms[1].formal_charge, 0);
        assert_eq!(multiple.atoms[2].formal_charge, -1);
        assert_eq!(multiple.bonds[0].bond_type, Single);
        assert_eq!(multiple.bonds[1].bond_type, Double);
    }

    #[test]
    fn source_port__inchi__valence5ncleanup1__line_393() {
        use BondDirection::{BeginDash, EndUpRight};
        use BondType::{Double, Single, ThreeCenter};

        let no_target = graph(&[(7, 0), (7, 0), (7, 1)], &[(0, 1, Double), (0, 2, Single)]);
        let mut no_target_actual = no_target.clone();
        assert_eq!(valence5n_cleanup1(&mut no_target_actual, 0), Ok(false));
        assert_eq!(no_target_actual, no_target);

        let mut direct = graph(&[(6, -3), (7, 1), (8, 2)], &[(1, 0, Double), (1, 2, Single)]);
        direct.atoms[0].num_explicit_hydrogens = 2;
        direct.atoms[1].num_explicit_hydrogens = 1;
        direct.bonds[0].direction = BeginDash;
        direct.bonds[1].direction = EndUpRight;
        assert_eq!(valence5n_cleanup1(&mut direct, 0), Ok(true));
        assert_eq!(direct.atoms[0].formal_charge, 1);
        assert_eq!(direct.atoms[1].formal_charge, 0);
        assert_eq!(direct.atoms[2].formal_charge, 2);
        assert_eq!(direct.atoms[0].num_explicit_hydrogens, 2);
        assert_eq!(direct.atoms[1].num_explicit_hydrogens, 1);
        assert_eq!(direct.bonds[0].bond_type, Single);
        assert_eq!(direct.bonds[0].direction, BeginDash);
        assert_eq!(direct.bonds[1].bond_type, Single);
        assert_eq!(direct.bonds[1].direction, EndUpRight);

        let mut alternating = graph(
            &[(7, -1), (6, 0), (6, 0), (7, 1)],
            &[(0, 1, Double), (1, 2, Single), (2, 3, Double)],
        );
        assert_eq!(valence5n_cleanup1(&mut alternating, 0), Ok(true));
        assert_eq!(alternating.atoms[0].formal_charge, 1);
        assert_eq!(alternating.atoms[3].formal_charge, 0);
        assert_eq!(
            alternating.bonds.iter().map(|bond| bond.bond_type).collect::<Vec<_>>(),
            [Single, Double, Single]
        );

        let over_limit = graph(
            &[(7, -1), (6, 0), (6, 0), (6, 0), (6, 0), (6, 0), (7, 1)],
            &[
                (0, 1, Double),
                (1, 2, Single),
                (2, 3, Double),
                (3, 4, Single),
                (4, 5, Double),
                (5, 6, Single),
            ],
        );
        let mut over_limit_actual = over_limit.clone();
        assert_eq!(valence5n_cleanup1(&mut over_limit_actual, 0), Ok(false));
        assert_eq!(over_limit_actual, over_limit);

        let mut valence_error = graph(&[(7, -1), (7, 1), (6, 0)], &[(0, 1, Double), (1, 2, ThreeCenter)]);
        let valence_error_before = valence_error.clone();
        assert_eq!(
            valence5n_cleanup1(&mut valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 1,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(valence_error.atoms[0], valence_error_before.atoms[0]);
        assert_eq!(valence_error.atoms[1].formal_charge, 0);
        assert_eq!(valence_error.atoms[2], valence_error_before.atoms[2]);
        assert_eq!(valence_error.bonds, valence_error_before.bonds);
    }

    #[test]
    fn source_port__inchi__valence5ncleanup2__line_417() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, ThreeCenter, Triple};

        let no_target = graph(
            &[(7, 0), (6, 4), (7, -1), (7, 0)],
            &[(0, 1, Double), (1, 2, Single), (0, 3, Triple)],
        );
        let mut no_target_actual = no_target.clone();
        assert_eq!(valence5n_cleanup2(&mut no_target_actual, 0), Ok(false));
        assert_eq!(no_target_actual, no_target);

        let too_long = graph(
            &[(7, 0), (6, 0), (6, 0), (7, -1)],
            &[(0, 1, Triple), (1, 2, Single), (2, 3, Double)],
        );
        let mut too_long_actual = too_long.clone();
        assert_eq!(valence5n_cleanup2(&mut too_long_actual, 0), Ok(false));
        assert_eq!(too_long_actual, too_long);

        let mut begin_is_root = graph(
            &[(7, 3), (6, 4), (7, -1), (8, 2)],
            &[(0, 1, Triple), (1, 2, Single), (1, 3, Single)],
        );
        begin_is_root.atoms[0].num_explicit_hydrogens = 1;
        begin_is_root.atoms[1].num_explicit_hydrogens = 2;
        begin_is_root.atoms[2].num_explicit_hydrogens = 3;
        begin_is_root.bonds[0].direction = BeginWedge;
        begin_is_root.bonds[1].direction = EndDownRight;
        begin_is_root.bonds[2].direction = BeginDash;
        assert_eq!(valence5n_cleanup2(&mut begin_is_root, 0), Ok(true));
        assert_eq!(begin_is_root.atoms[0].formal_charge, 3);
        assert_eq!(begin_is_root.atoms[1].formal_charge, -1);
        assert_eq!(begin_is_root.atoms[2].formal_charge, 0);
        assert_eq!(begin_is_root.atoms[3].formal_charge, 2);
        assert_eq!(begin_is_root.atoms[0].num_explicit_hydrogens, 1);
        assert_eq!(begin_is_root.atoms[1].num_explicit_hydrogens, 2);
        assert_eq!(begin_is_root.atoms[2].num_explicit_hydrogens, 3);
        assert_eq!(begin_is_root.bonds[0].bond_type, Single);
        assert_eq!(begin_is_root.bonds[0].direction, BeginWedge);
        assert_eq!(begin_is_root.bonds[1].bond_type, Double);
        assert_eq!(begin_is_root.bonds[1].direction, EndDownRight);
        assert_eq!(begin_is_root.bonds[2].bond_type, Single);
        assert_eq!(begin_is_root.bonds[2].direction, BeginDash);

        let mut end_is_root = graph(
            &[(7, 0), (6, 0), (7, -1), (7, -1)],
            &[(1, 0, Triple), (2, 1, Single), (1, 3, Single)],
        );
        end_is_root.bonds[0].direction = EndUpRight;
        assert_eq!(valence5n_cleanup2(&mut end_is_root, 0), Ok(true));
        assert_eq!(end_is_root.atoms[0].formal_charge, 0);
        assert_eq!(end_is_root.atoms[1].formal_charge, -1);
        assert_eq!(end_is_root.atoms[2].formal_charge, 0);
        assert_eq!(end_is_root.atoms[3].formal_charge, -1);
        assert_eq!(end_is_root.bonds[0].bond_type, Single);
        assert_eq!(end_is_root.bonds[0].direction, EndUpRight);
        assert_eq!(end_is_root.bonds[1].bond_type, Double);
        assert_eq!(end_is_root.bonds[2].bond_type, Single);

        let mut target_valence_error = graph(
            &[(7, 0), (6, 0), (7, -1), (6, 0)],
            &[(0, 1, Triple), (1, 2, Single), (2, 3, ThreeCenter)],
        );
        let target_error_before = target_valence_error.clone();
        assert_eq!(
            valence5n_cleanup2(&mut target_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 2,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(target_valence_error.atoms[0], target_error_before.atoms[0]);
        assert_eq!(target_valence_error.atoms[1].formal_charge, -1);
        assert_eq!(target_valence_error.atoms[2].formal_charge, 0);
        assert_eq!(target_valence_error.atoms[3], target_error_before.atoms[3]);
        assert_eq!(target_valence_error.bonds[0].bond_type, Single);
        assert_eq!(target_valence_error.bonds[1].bond_type, Double);
        assert_eq!(target_valence_error.bonds[2], target_error_before.bonds[2]);

        let mut root_valence_error = graph(
            &[(7, 0), (6, 0), (7, -1), (6, 0)],
            &[(0, 1, Triple), (1, 2, Single), (0, 3, ThreeCenter)],
        );
        assert_eq!(
            valence5n_cleanup2(&mut root_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 2,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(root_valence_error.atoms[0].formal_charge, 0);
        assert_eq!(root_valence_error.atoms[1].formal_charge, -1);
        assert_eq!(root_valence_error.atoms[2].formal_charge, 0);
        assert_eq!(root_valence_error.bonds[0].bond_type, Single);
        assert_eq!(root_valence_error.bonds[1].bond_type, Double);
        assert_eq!(root_valence_error.bonds[2].bond_type, ThreeCenter);
    }

    #[test]
    fn source_port__inchi__valence5ncleanup3__line_443() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, ThreeCenter, Triple};

        let no_nitrogen = graph(
            &[(7, 0), (7, 1), (8, 0), (7, 0)],
            &[(0, 1, Double), (0, 2, Double), (0, 3, Single)],
        );
        let mut no_nitrogen_actual = no_nitrogen.clone();
        assert_eq!(valence5n_cleanup3(&mut no_nitrogen_actual, 0), Ok(false));
        assert_eq!(no_nitrogen_actual, no_nitrogen);

        let mut oxygen_guard = graph(
            &[(7, -3), (7, 0), (8, 0), (6, 2)],
            &[(1, 0, Double), (0, 2, Double), (0, 3, Triple)],
        );
        oxygen_guard.atoms[0].num_explicit_hydrogens = 1;
        oxygen_guard.atoms[1].num_explicit_hydrogens = 2;
        oxygen_guard.atoms[2].num_explicit_hydrogens = 3;
        oxygen_guard.bonds[0].direction = BeginWedge;
        oxygen_guard.bonds[1].direction = EndDownRight;
        oxygen_guard.bonds[2].direction = BeginDash;
        let oxygen_guard_before = oxygen_guard.clone();
        assert_eq!(valence5n_cleanup3(&mut oxygen_guard, 0), Ok(true));
        assert_eq!(oxygen_guard, oxygen_guard_before);

        let mut mutate = graph(
            &[(7, -3), (7, 0), (8, 1), (6, 2)],
            &[(1, 0, Double), (0, 2, Double), (0, 3, Single)],
        );
        mutate.atoms[0].num_explicit_hydrogens = 1;
        mutate.atoms[1].num_explicit_hydrogens = 2;
        mutate.bonds[0].direction = EndUpRight;
        mutate.bonds[1].direction = EndDownRight;
        assert_eq!(valence5n_cleanup3(&mut mutate, 0), Ok(true));
        assert_eq!(mutate.atoms[0].formal_charge, 1);
        assert_eq!(mutate.atoms[1].formal_charge, -1);
        assert_eq!(mutate.atoms[2].formal_charge, 1);
        assert_eq!(mutate.atoms[3].formal_charge, 2);
        assert_eq!(mutate.atoms[0].num_explicit_hydrogens, 1);
        assert_eq!(mutate.atoms[1].num_explicit_hydrogens, 2);
        assert_eq!(mutate.bonds[0].bond_type, Single);
        assert_eq!(mutate.bonds[0].direction, EndUpRight);
        assert_eq!(mutate.bonds[1].bond_type, Double);
        assert_eq!(mutate.bonds[1].direction, EndDownRight);
        assert_eq!(mutate.bonds[2].bond_type, Single);

        let mut target_valence_error = graph(&[(7, 0), (7, 0), (6, 0)], &[(0, 1, Double), (1, 2, ThreeCenter)]);
        let target_error_before = target_valence_error.clone();
        assert_eq!(
            valence5n_cleanup3(&mut target_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 1,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(target_valence_error.atoms[0], target_error_before.atoms[0]);
        assert_eq!(target_valence_error.atoms[1].formal_charge, -1);
        assert_eq!(target_valence_error.atoms[2], target_error_before.atoms[2]);
        assert_eq!(target_valence_error.bonds, target_error_before.bonds);

        let mut root_valence_error = graph(&[(7, 0), (7, 0), (6, 0)], &[(0, 1, Double), (0, 2, ThreeCenter)]);
        root_valence_error.bonds[0].direction = BeginDash;
        assert_eq!(
            valence5n_cleanup3(&mut root_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 1,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(root_valence_error.atoms[0].formal_charge, 1);
        assert_eq!(root_valence_error.atoms[1].formal_charge, -1);
        assert_eq!(root_valence_error.bonds[0].bond_type, Single);
        assert_eq!(root_valence_error.bonds[0].direction, BeginDash);
        assert_eq!(root_valence_error.bonds[1].bond_type, ThreeCenter);
    }

    #[test]
    fn source_port__inchi__valence5ncleanup4__line_477() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, Triple};

        let no_matching_silicon = graph(
            &[(7, 3), (14, -1), (14, 0), (6, -1)],
            &[(0, 1, Single), (0, 2, Double), (0, 3, Double)],
        );
        let mut no_matching_silicon_actual = no_matching_silicon.clone();
        assert!(!valence5n_cleanup4(&mut no_matching_silicon_actual, 0));
        assert_eq!(no_matching_silicon_actual, no_matching_silicon);

        let one_matching_silicon = graph(
            &[(7, 3), (14, -1), (14, -1), (14, 0)],
            &[(0, 1, Double), (0, 2, Single), (0, 3, Double)],
        );
        let mut one_matching_silicon_actual = one_matching_silicon.clone();
        assert!(!valence5n_cleanup4(&mut one_matching_silicon_actual, 0));
        assert_eq!(one_matching_silicon_actual, one_matching_silicon);

        let mut two_matching_silicons = graph(
            &[(7, 3), (6, 0), (14, -1), (14, 0), (14, -1), (16, 1), (6, 0)],
            &[
                (0, 1, Double),
                (2, 0, Double),
                (0, 3, Double),
                (4, 0, Double),
                (5, 6, Double),
            ],
        );
        two_matching_silicons.atoms[2].num_explicit_hydrogens = 1;
        two_matching_silicons.atoms[4].num_explicit_hydrogens = 2;
        two_matching_silicons.bonds[0].direction = BeginWedge;
        two_matching_silicons.bonds[1].direction = EndUpRight;
        two_matching_silicons.bonds[2].direction = BeginDash;
        two_matching_silicons.bonds[3].direction = EndDownRight;
        two_matching_silicons.bonds[4].direction = BeginWedge;
        let two_matching_silicons_before = two_matching_silicons.clone();
        assert!(valence5n_cleanup4(&mut two_matching_silicons, 0));
        assert_eq!(two_matching_silicons.atoms[0], two_matching_silicons_before.atoms[0]);
        assert_eq!(two_matching_silicons.atoms[1], two_matching_silicons_before.atoms[1]);
        assert_eq!(two_matching_silicons.atoms[2].formal_charge, 0);
        assert_eq!(two_matching_silicons.atoms[2].num_explicit_hydrogens, 1);
        assert_eq!(two_matching_silicons.atoms[3], two_matching_silicons_before.atoms[3]);
        assert_eq!(two_matching_silicons.atoms[4].formal_charge, 0);
        assert_eq!(two_matching_silicons.atoms[4].num_explicit_hydrogens, 2);
        assert_eq!(two_matching_silicons.atoms[5], two_matching_silicons_before.atoms[5]);
        assert_eq!(two_matching_silicons.atoms[6], two_matching_silicons_before.atoms[6]);
        assert_eq!(two_matching_silicons.bonds[0], two_matching_silicons_before.bonds[0]);
        assert_eq!(two_matching_silicons.bonds[1].bond_type, Single);
        assert_eq!(two_matching_silicons.bonds[1].direction, EndUpRight);
        assert_eq!(two_matching_silicons.bonds[2], two_matching_silicons_before.bonds[2]);
        assert_eq!(two_matching_silicons.bonds[3].bond_type, Single);
        assert_eq!(two_matching_silicons.bonds[3].direction, EndDownRight);
        assert_eq!(two_matching_silicons.bonds[4], two_matching_silicons_before.bonds[4]);

        let three_matching_silicons = graph(
            &[(7, 3), (14, -1), (14, -1), (14, -1)],
            &[(1, 0, Double), (0, 2, Double), (3, 0, Double)],
        );
        let mut three_matching_silicons_actual = three_matching_silicons.clone();
        assert!(!valence5n_cleanup4(&mut three_matching_silicons_actual, 0));
        assert_eq!(three_matching_silicons_actual, three_matching_silicons);

        let mut exactly_two_before_late_nonmatch = graph(
            &[(7, -4), (14, -1), (14, -1), (14, -1)],
            &[(0, 1, Double), (0, 2, Double), (0, 3, Triple)],
        );
        exactly_two_before_late_nonmatch.bonds[0].direction = BeginDash;
        exactly_two_before_late_nonmatch.bonds[1].direction = EndUpRight;
        let late_nonmatch_before = exactly_two_before_late_nonmatch.clone();
        assert!(valence5n_cleanup4(&mut exactly_two_before_late_nonmatch, 0));
        assert_eq!(exactly_two_before_late_nonmatch.atoms[0].formal_charge, -4);
        assert_eq!(exactly_two_before_late_nonmatch.atoms[1].formal_charge, 0);
        assert_eq!(exactly_two_before_late_nonmatch.atoms[2].formal_charge, 0);
        assert_eq!(exactly_two_before_late_nonmatch.atoms[3].formal_charge, -1);
        assert_eq!(exactly_two_before_late_nonmatch.bonds[0].bond_type, Single);
        assert_eq!(exactly_two_before_late_nonmatch.bonds[0].direction, BeginDash);
        assert_eq!(exactly_two_before_late_nonmatch.bonds[1].bond_type, Single);
        assert_eq!(exactly_two_before_late_nonmatch.bonds[1].direction, EndUpRight);
        assert_eq!(exactly_two_before_late_nonmatch.bonds[2], late_nonmatch_before.bonds[2]);
    }

    #[test]
    fn source_port__inchi__valence5ncleanup5__line_544() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, ThreeCenter};

        let precondition_input = graph(&[(7, 0), (8, 0)], &[(0, 1, Double)]);
        let mut precondition_actual = precondition_input.clone();
        assert_eq!(
            valence5n_cleanup5(&mut precondition_actual, 0, 7),
            Err(AdapterCleanup5Error::Precondition {
                kind: "Pre-condition Violation",
                expression: "atomicNum == 8 || atomicNum == 16 || atomicNum == 9 || atomicNum == 17",
                message: "this cleanup looks for O or S or Cl or F",
            })
        );
        assert_eq!(precondition_actual, precondition_input);

        for atomic_number in [8, 16, 9, 17] {
            let no_target = graph(&[(7, -2), (6, 1)], &[(0, 1, Double)]);
            let mut no_target_actual = no_target.clone();
            assert_eq!(valence5n_cleanup5(&mut no_target_actual, 0, atomic_number), Ok(false));
            assert_eq!(no_target_actual, no_target);
        }

        let mut uncharged = graph(
            &[(7, -3), (6, 2), (6, -2), (8, 0)],
            &[(0, 1, Double), (1, 2, Single), (2, 3, Double)],
        );
        uncharged.atoms[0].num_explicit_hydrogens = 1;
        uncharged.atoms[3].num_explicit_hydrogens = 4;
        uncharged.bonds[0].direction = BeginDash;
        uncharged.bonds[1].direction = EndUpRight;
        uncharged.bonds[2].direction = BeginWedge;
        assert_eq!(valence5n_cleanup5(&mut uncharged, 0, 8), Ok(true));
        assert_eq!(uncharged.atoms[0].formal_charge, 1);
        assert_eq!(uncharged.atoms[0].num_explicit_hydrogens, 1);
        assert_eq!(uncharged.atoms[1].formal_charge, 2);
        assert_eq!(uncharged.atoms[2].formal_charge, -2);
        assert_eq!(uncharged.atoms[3].formal_charge, -1);
        assert_eq!(uncharged.atoms[3].num_explicit_hydrogens, 4);
        assert_eq!(uncharged.bonds[0].bond_type, Single);
        assert_eq!(uncharged.bonds[0].direction, BeginDash);
        assert_eq!(uncharged.bonds[1].bond_type, Double);
        assert_eq!(uncharged.bonds[1].direction, EndUpRight);
        assert_eq!(uncharged.bonds[2].bond_type, Single);
        assert_eq!(uncharged.bonds[2].direction, BeginWedge);

        let mut charged = graph(&[(7, -3), (16, 1)], &[(1, 0, Double)]);
        charged.atoms[1].num_explicit_hydrogens = 3;
        charged.bonds[0].direction = EndDownRight;
        assert_eq!(valence5n_cleanup5(&mut charged, 0, 16), Ok(true));
        assert_eq!(charged.atoms[0].formal_charge, 1);
        assert_eq!(charged.atoms[1].formal_charge, 0);
        assert_eq!(charged.atoms[1].num_explicit_hydrogens, 3);
        assert_eq!(charged.bonds[0].bond_type, Single);
        assert_eq!(charged.bonds[0].direction, EndDownRight);

        let mut both = graph(
            &[(7, -3), (9, 0), (9, 1), (6, 4)],
            &[(0, 1, Double), (2, 0, Double), (0, 3, Single)],
        );
        both.atoms[1].num_explicit_hydrogens = 5;
        both.atoms[2].num_explicit_hydrogens = 2;
        both.bonds[0].direction = BeginWedge;
        both.bonds[1].direction = EndUpRight;
        let both_before = both.clone();
        assert_eq!(valence5n_cleanup5(&mut both, 0, 9), Ok(true));
        assert_eq!(both.atoms[0].formal_charge, 1);
        assert_eq!(both.atoms[1].formal_charge, 0);
        assert_eq!(both.atoms[1].num_explicit_hydrogens, 1);
        assert_eq!(both.atoms[2].formal_charge, 0);
        assert_eq!(both.atoms[2].num_explicit_hydrogens, 0);
        assert_eq!(both.atoms[3], both_before.atoms[3]);
        assert_eq!(both.bonds[0].bond_type, Single);
        assert_eq!(both.bonds[0].direction, BeginWedge);
        assert_eq!(both.bonds[1], both_before.bonds[1]);
        assert_eq!(both.bonds[2], both_before.bonds[2]);

        let mut depth_seven = graph(
            &[(7, 0), (6, 0), (6, 0), (6, 0), (6, 0), (6, 0), (6, 0), (17, 0)],
            &[
                (0, 1, Double),
                (1, 2, Single),
                (2, 3, Double),
                (3, 4, Single),
                (4, 5, Double),
                (5, 6, Single),
                (6, 7, Double),
            ],
        );
        assert_eq!(valence5n_cleanup5(&mut depth_seven, 0, 17), Ok(true));
        assert_eq!(depth_seven.atoms[0].formal_charge, 1);
        assert_eq!(depth_seven.atoms[7].formal_charge, -1);
        assert_eq!(
            depth_seven.bonds.iter().map(|bond| bond.bond_type).collect::<Vec<_>>(),
            [Single, Double, Single, Double, Single, Double, Single]
        );

        let over_depth = graph(
            &[
                (7, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (8, 0),
            ],
            &[
                (0, 1, Double),
                (1, 2, Single),
                (2, 3, Double),
                (3, 4, Single),
                (4, 5, Double),
                (5, 6, Single),
                (6, 7, Double),
                (7, 8, Single),
                (8, 9, Double),
            ],
        );
        let mut over_depth_actual = over_depth.clone();
        assert_eq!(valence5n_cleanup5(&mut over_depth_actual, 0, 8), Ok(false));
        assert_eq!(over_depth_actual, over_depth);

        let mut only_uncharged_error = graph(&[(7, 0), (8, 0), (6, 0)], &[(0, 1, Double), (1, 2, ThreeCenter)]);
        assert_eq!(
            valence5n_cleanup5(&mut only_uncharged_error, 0, 8),
            Err(AdapterCleanup5Error::Valence(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 1,
                bond_type: ThreeCenter,
            }))
        );
        assert_eq!(only_uncharged_error.atoms[0].formal_charge, 1);
        assert_eq!(only_uncharged_error.atoms[1].formal_charge, -1);
        assert_eq!(only_uncharged_error.bonds[0].bond_type, Single);

        let mut only_charged_error = graph(&[(7, 0), (8, 1), (6, 0)], &[(0, 1, Double), (1, 2, ThreeCenter)]);
        only_charged_error.atoms[1].num_explicit_hydrogens = 2;
        assert_eq!(
            valence5n_cleanup5(&mut only_charged_error, 0, 8),
            Err(AdapterCleanup5Error::Valence(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 1,
                bond_type: ThreeCenter,
            }))
        );
        assert_eq!(only_charged_error.atoms[0].formal_charge, 1);
        assert_eq!(only_charged_error.atoms[1].formal_charge, 0);
        assert_eq!(only_charged_error.atoms[1].num_explicit_hydrogens, 2);
        assert_eq!(only_charged_error.bonds[0].bond_type, Single);

        let mut charged_first_error = graph(
            &[(7, 0), (8, 0), (8, 1), (6, 0), (6, 0)],
            &[(0, 1, Double), (0, 2, Double), (1, 3, ThreeCenter), (2, 4, ThreeCenter)],
        );
        charged_first_error.atoms[1].num_explicit_hydrogens = 3;
        charged_first_error.atoms[2].num_explicit_hydrogens = 4;
        assert_eq!(
            valence5n_cleanup5(&mut charged_first_error, 0, 8),
            Err(AdapterCleanup5Error::Valence(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 3,
                bond_type: ThreeCenter,
            }))
        );
        assert_eq!(charged_first_error.atoms[0].formal_charge, 1);
        assert_eq!(charged_first_error.atoms[1].num_explicit_hydrogens, 1);
        assert_eq!(charged_first_error.atoms[2].formal_charge, 0);
        assert_eq!(charged_first_error.atoms[2].num_explicit_hydrogens, 0);
        assert_eq!(charged_first_error.bonds[0].bond_type, Single);
        assert_eq!(charged_first_error.bonds[1].bond_type, Double);

        let mut uncharged_second_error = graph(
            &[(7, 0), (8, 0), (8, 1), (6, 0)],
            &[(0, 1, Double), (0, 2, Double), (1, 3, ThreeCenter)],
        );
        assert_eq!(
            valence5n_cleanup5(&mut uncharged_second_error, 0, 8),
            Err(AdapterCleanup5Error::Valence(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 2,
                bond_type: ThreeCenter,
            }))
        );
        assert_eq!(uncharged_second_error.atoms[0].formal_charge, 1);
        assert_eq!(uncharged_second_error.atoms[1].num_explicit_hydrogens, 1);
        assert_eq!(uncharged_second_error.atoms[2].formal_charge, 0);
        assert_eq!(uncharged_second_error.atoms[2].num_explicit_hydrogens, 0);
        assert_eq!(uncharged_second_error.bonds[0].bond_type, Single);
        assert_eq!(uncharged_second_error.bonds[1].bond_type, Double);
    }

    #[test]
    fn source_port__inchi__valence5ncleanup6__line_625() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{
            Aromatic, Dative, DativeL, DativeOne, DativeR, Double, FiveAndAHalf, FourAndAHalf, Hextuple, Hydrogen,
            Ionic, OneAndAHalf, Other, Quadruple, Quintuple, Single, ThreeAndAHalf, ThreeCenter, Triple, TwoAndAHalf,
            Unspecified, Zero,
        };

        let wrong_element = graph(&[(6, 0), (6, 0)], &[(0, 1, ThreeCenter)]);
        let mut wrong_element_actual = wrong_element.clone();
        assert_eq!(valence5n_cleanup6(&mut wrong_element_actual, 0), Ok(false));
        assert_eq!(wrong_element_actual, wrong_element);

        let wrong_charge = graph(&[(7, 1), (6, 0)], &[(0, 1, ThreeCenter)]);
        let mut wrong_charge_actual = wrong_charge.clone();
        assert_eq!(valence5n_cleanup6(&mut wrong_charge_actual, 0), Ok(false));
        assert_eq!(wrong_charge_actual, wrong_charge);

        let wrong_valence = graph(
            &[(7, 0), (6, 0), (6, 0), (6, 0), (6, 0)],
            &[(0, 1, Single), (0, 2, Single), (0, 3, Single), (0, 4, Single)],
        );
        let mut wrong_valence_actual = wrong_valence.clone();
        assert_eq!(valence5n_cleanup6(&mut wrong_valence_actual, 0), Ok(false));
        assert_eq!(wrong_valence_actual, wrong_valence);

        let valence_error = graph(&[(7, 0), (6, 0)], &[(0, 1, ThreeCenter)]);
        let mut valence_error_actual = valence_error.clone();
        assert_eq!(
            valence5n_cleanup6(&mut valence_error_actual, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(valence_error_actual, valence_error);

        for atom_index in [0_usize, 1, 3, 4, 5, 6] {
            let mut no_match = cleanup6_graph(7, Triple);
            no_match.atoms[atom_index].atomic_number = 8;
            let before = no_match.clone();
            assert_eq!(valence5n_cleanup6(&mut no_match, 2), Ok(false));
            assert_eq!(no_match, before);
        }
        for bond_index in [0_usize, 4, 5] {
            let mut no_match = cleanup6_graph(7, Triple);
            no_match.bonds[bond_index].bond_type = if bond_index == 5 { Single } else { Double };
            let before = no_match.clone();
            assert_eq!(valence5n_cleanup6(&mut no_match, 2), Ok(false));
            assert_eq!(no_match, before);
        }
        for root_bond_types in [[Triple, Single, Single], [Double, Single, Double]] {
            let mut no_match = cleanup6_graph(7, Triple);
            no_match.bonds[1].bond_type = root_bond_types[0];
            no_match.bonds[2].bond_type = root_bond_types[1];
            no_match.bonds[6].bond_type = root_bond_types[2];
            let before = no_match.clone();
            assert_eq!(valence5n_cleanup6(&mut no_match, 2), Ok(false));
            assert_eq!(no_match, before);
        }

        let bridge_types = [
            Unspecified,
            Single,
            Double,
            Triple,
            Quadruple,
            Quintuple,
            Hextuple,
            OneAndAHalf,
            TwoAndAHalf,
            ThreeAndAHalf,
            FourAndAHalf,
            FiveAndAHalf,
            Aromatic,
            Ionic,
            Hydrogen,
            ThreeCenter,
            DativeOne,
            Dative,
            DativeL,
            DativeR,
            Other,
            Zero,
        ];
        for (case_index, bridge_type) in bridge_types.into_iter().enumerate() {
            let mut unique = cleanup6_graph(7, bridge_type);
            unique.atoms[0].num_explicit_hydrogens = 1;
            unique.atoms[4].num_explicit_hydrogens = 2;
            unique.atoms[5].num_explicit_hydrogens = 3;
            unique.bonds[0].direction = BeginWedge;
            unique.bonds[1].direction = BeginDash;
            unique.bonds[2].direction = EndDownRight;
            unique.bonds[3].direction = EndUpRight;
            let mut expected = unique.clone();
            expected.atoms[2].formal_charge = 1;
            expected.bonds[0].bond_type = Double;
            expected.bonds[1].bond_type = Single;
            expected.bonds[4].bond_type = Double;
            expected.bonds[5].bond_type = Single;

            assert_eq!(valence5n_cleanup6(&mut unique, 2), Ok(true), "bridge case {case_index}");
            assert_eq!(unique, expected, "bridge case {case_index}");
            assert_eq!(unique.bonds[3].bond_type, bridge_type);
        }

        let first = cleanup6_graph(7, Triple);
        let second = cleanup6_graph(50, Aromatic);
        let multiple = disjoint_union(&first, &second);
        let mut multiple_actual = multiple.clone();
        assert_eq!(valence5n_cleanup6(&mut multiple_actual, 2), Ok(false));
        assert_eq!(multiple_actual, multiple);

        let target = graph(
            &[(7, 0), (6, 0), (6, 0), (6, 0), (6, 0), (6, 0)],
            &[
                (0, 1, Single),
                (2, 0, Single),
                (0, 3, Single),
                (4, 0, Single),
                (0, 5, Single),
            ],
        );
        let existing_tin_match = cleanup6_graph(50, Triple);
        let mut unrelated_match = disjoint_union(&target, &existing_tin_match);
        let mut expected = unrelated_match.clone();
        expected.atoms[0].formal_charge = 1;
        expected.bonds[5].bond_type = Double;
        expected.bonds[6].bond_type = Single;
        expected.bonds[9].bond_type = Double;
        expected.bonds[10].bond_type = Single;
        assert_eq!(valence5n_cleanup6(&mut unrelated_match, 0), Ok(true));
        assert_eq!(unrelated_match, expected);
        assert_eq!(unrelated_match.atoms[8].atomic_number, 50);
    }

    #[test]
    fn source_port__inchi__valence5ncleanup7__line_679() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{
            Aromatic, Dative, DativeL, DativeOne, DativeR, Double, FiveAndAHalf, FourAndAHalf, Hextuple, Hydrogen,
            Ionic, OneAndAHalf, Other, Quadruple, Quintuple, Single, ThreeAndAHalf, ThreeCenter, Triple, TwoAndAHalf,
            Unspecified, Zero,
        };

        let no_target = graph(&[(7, 0), (6, 0)], &[(0, 1, ThreeCenter)]);
        let mut no_target_actual = no_target.clone();
        assert_eq!(valence5n_cleanup7(&mut no_target_actual, 0), Ok(false));
        assert_eq!(no_target_actual, no_target);

        let mut wrong_element = cleanup7_graph(6, Triple);
        let wrong_element_before = wrong_element.clone();
        assert_eq!(valence5n_cleanup7(&mut wrong_element, 2), Ok(false));
        assert_eq!(wrong_element, wrong_element_before);

        let mut wrong_charge = cleanup7_graph(7, Triple);
        wrong_charge.atoms[2].formal_charge = 1;
        let wrong_charge_before = wrong_charge.clone();
        assert_eq!(valence5n_cleanup7(&mut wrong_charge, 2), Ok(false));
        assert_eq!(wrong_charge, wrong_charge_before);

        let mut wrong_valence = cleanup7_graph(7, Triple);
        wrong_valence.atoms[2].num_explicit_hydrogens = 1;
        let wrong_valence_before = wrong_valence.clone();
        assert_eq!(valence5n_cleanup7(&mut wrong_valence, 2), Ok(false));
        assert_eq!(wrong_valence, wrong_valence_before);

        let mut valence_error = graph(
            &[(7, 0), (6, 0), (6, 0), (8, 0), (6, 0)],
            &[(0, 1, Double), (1, 2, Single), (2, 3, Double), (0, 4, ThreeCenter)],
        );
        let valence_error_before = valence_error.clone();
        assert_eq!(
            valence5n_cleanup7(&mut valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 3,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(valence_error, valence_error_before);

        for atom_index in [0_usize, 1, 3, 4, 5, 6] {
            let mut no_match = cleanup7_graph(7, Triple);
            no_match.atoms[atom_index].atomic_number = 9;
            let before = no_match.clone();
            assert_eq!(valence5n_cleanup7(&mut no_match, 2), Ok(false));
            assert_eq!(no_match, before, "atom predicate {atom_index}");
        }
        for bond_index in [4_usize, 5] {
            let mut no_match = cleanup7_graph(7, Triple);
            no_match.bonds[bond_index].bond_type = Double;
            let before = no_match.clone();
            assert_eq!(valence5n_cleanup7(&mut no_match, 2), Ok(false));
            assert_eq!(no_match, before, "bond predicate {bond_index}");
        }

        let bridge_types = [
            Unspecified,
            Single,
            Double,
            Triple,
            Quadruple,
            Quintuple,
            Hextuple,
            OneAndAHalf,
            TwoAndAHalf,
            ThreeAndAHalf,
            FourAndAHalf,
            FiveAndAHalf,
            Aromatic,
            Ionic,
            Hydrogen,
            ThreeCenter,
            DativeOne,
            Dative,
            DativeL,
            DativeR,
            Other,
            Zero,
        ];
        for (case_index, bridge_type) in bridge_types.into_iter().enumerate() {
            let mut unique = cleanup7_graph(7, bridge_type);
            unique.atoms[0].num_explicit_hydrogens = 1;
            unique.atoms[5].num_explicit_hydrogens = 2;
            unique.bonds[0].direction = BeginWedge;
            unique.bonds[1].direction = BeginDash;
            unique.bonds[2].direction = EndDownRight;
            unique.bonds[7].direction = EndUpRight;
            let mut expected = unique.clone();
            expected.bonds[1].bond_type = Single;
            expected.bonds[2].bond_type = Single;
            expected.bonds[3].bond_type = Double;
            expected.bonds[7].bond_type = Single;
            expected.atoms[7].formal_charge = -1;

            assert_eq!(valence5n_cleanup7(&mut unique, 2), Ok(true), "bridge case {case_index}");
            assert_eq!(unique, expected, "bridge case {case_index}");
            assert_eq!(unique.bonds[0].bond_type, bridge_type);
        }

        let mut over_depth = cleanup7_graph(7, Triple);
        over_depth.bonds.pop();
        over_depth = AdapterMol::from_graph(over_depth.atoms, over_depth.bonds);
        over_depth.atoms.extend([
            AdapterAtom {
                atomic_number: 6,
                formal_charge: 0,
                num_explicit_hydrogens: 0,
                is_aromatic: false,
                ..AdapterAtom::default()
            },
            AdapterAtom {
                atomic_number: 6,
                formal_charge: 0,
                num_explicit_hydrogens: 0,
                is_aromatic: false,
                ..AdapterAtom::default()
            },
            AdapterAtom {
                atomic_number: 6,
                formal_charge: 0,
                num_explicit_hydrogens: 0,
                is_aromatic: false,
                ..AdapterAtom::default()
            },
        ]);
        let mut over_depth_bonds = over_depth.bonds.clone();
        over_depth_bonds.extend([
            AdapterBond {
                begin_atom_index: 4,
                end_atom_index: 8,
                bond_type: Double,
                direction: BondDirection::None,
                ..AdapterBond::default()
            },
            AdapterBond {
                begin_atom_index: 8,
                end_atom_index: 9,
                bond_type: Single,
                direction: BondDirection::None,
                ..AdapterBond::default()
            },
            AdapterBond {
                begin_atom_index: 9,
                end_atom_index: 10,
                bond_type: Double,
                direction: BondDirection::None,
                ..AdapterBond::default()
            },
        ]);
        over_depth = AdapterMol::from_graph(over_depth.atoms, over_depth_bonds);
        let over_depth_before = over_depth.clone();
        assert_eq!(valence5n_cleanup7(&mut over_depth, 2), Ok(false));
        assert_eq!(over_depth, over_depth_before);

        let first = cleanup7_graph(7, Triple);
        let second = cleanup7_graph(50, Aromatic);
        let multiple = disjoint_union(&first, &second);
        let mut multiple_actual = multiple.clone();
        assert_eq!(valence5n_cleanup7(&mut multiple_actual, 2), Ok(false));
        assert_eq!(multiple_actual, multiple);

        let target_component = graph(
            &[(7, 0), (6, 0), (6, 0), (8, 0), (6, 0), (6, 0)],
            &[
                (0, 1, Double),
                (1, 2, Single),
                (2, 3, Double),
                (0, 4, Double),
                (0, 5, Single),
            ],
        );
        let match_component = cleanup7_graph(50, Aromatic);
        let mut unrelated_match = disjoint_union(&target_component, &match_component);
        let mut expected = unrelated_match.clone();
        expected.bonds[0].bond_type = Single;
        expected.bonds[1].bond_type = Double;
        expected.bonds[2].bond_type = Single;
        expected.bonds[6].bond_type = Single;
        expected.atoms[3].formal_charge = -1;
        assert_eq!(valence5n_cleanup7(&mut unrelated_match, 0), Ok(true));
        assert_eq!(unrelated_match, expected);
        assert_eq!(unrelated_match.atoms[8].atomic_number, 50);
    }

    #[test]
    fn source_port__inchi__valence5ncleanup8__line_750() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, ThreeCenter};

        let mut wrong_element = cleanup8_graph(6);
        let before = wrong_element.clone();
        assert_eq!(valence5n_cleanup8(&mut wrong_element, 5), Ok(false));
        assert_eq!(wrong_element, before);
        let mut wrong_charge = cleanup8_graph(7);
        wrong_charge.atoms[5].formal_charge = 1;
        let before = wrong_charge.clone();
        assert_eq!(valence5n_cleanup8(&mut wrong_charge, 5), Ok(false));
        assert_eq!(wrong_charge, before);
        let mut wrong_valence = cleanup8_graph(7);
        wrong_valence.atoms[5].num_explicit_hydrogens = 2;
        let before = wrong_valence.clone();
        assert_eq!(valence5n_cleanup8(&mut wrong_valence, 5), Ok(false));
        assert_eq!(wrong_valence, before);

        let mut valence_error = graph(&[(7, 0), (6, 0)], &[(0, 1, ThreeCenter)]);
        let before = valence_error.clone();
        assert_eq!(
            valence5n_cleanup8(&mut valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(valence_error, before);

        for atom_index in 0_usize..5 {
            let mut no_match = cleanup8_graph(7);
            no_match.atoms[atom_index].atomic_number = 8;
            let before = no_match.clone();
            assert_eq!(valence5n_cleanup8(&mut no_match, 5), Ok(false));
            assert_eq!(no_match, before, "atom predicate {atom_index}");
        }
        for bond_index in 0_usize..5 {
            let mut no_match = cleanup8_graph(7);
            no_match.bonds[bond_index].bond_type = if no_match.bonds[bond_index].bond_type == Single {
                Double
            } else {
                Single
            };
            let before = no_match.clone();
            assert_eq!(valence5n_cleanup8(&mut no_match, 5), Ok(false));
            assert_eq!(no_match, before, "bond predicate {bond_index}");
        }
        let mut root_bond_miss = cleanup8_graph(7);
        root_bond_miss.bonds[5].bond_type = Single;
        root_bond_miss.atoms[5].num_explicit_hydrogens = 4;
        let before = root_bond_miss.clone();
        assert_eq!(valence5n_cleanup8(&mut root_bond_miss, 5), Ok(false));
        assert_eq!(root_bond_miss, before);

        let mut unique = cleanup8_graph(7);
        unique.bonds[0].direction = BeginWedge;
        unique.bonds[1].direction = BeginDash;
        unique.bonds[2].direction = EndDownRight;
        unique.bonds[5].direction = EndUpRight;
        let mut expected = unique.clone();
        expected.bonds[1].bond_type = Single;
        expected.bonds[2].bond_type = Double;
        expected.bonds[3].bond_type = Single;
        expected.bonds[4].bond_type = Double;
        expected.bonds[5].bond_type = Single;
        expected.atoms[1].formal_charge = -1;
        expected.atoms[5].formal_charge = 1;
        assert_eq!(valence5n_cleanup8(&mut unique, 5), Ok(true));
        assert_eq!(unique, expected);

        let first = cleanup8_graph(7);
        let second = cleanup8_graph(50);
        let multiple = disjoint_union(&first, &second);
        let mut multiple_actual = multiple.clone();
        assert_eq!(valence5n_cleanup8(&mut multiple_actual, 5), Ok(false));
        assert_eq!(multiple_actual, multiple);

        let target_component = graph(
            &[(7, 0), (6, 0), (6, 0), (6, 0), (6, 0)],
            &[(0, 1, Double), (0, 2, Single), (0, 3, Single), (0, 4, Single)],
        );
        let remote = cleanup8_graph(50);
        let mut remote_match = disjoint_union(&target_component, &remote);
        let mut expected = remote_match.clone();
        expected.bonds[5].bond_type = Single;
        expected.bonds[6].bond_type = Double;
        expected.bonds[7].bond_type = Single;
        expected.bonds[8].bond_type = Double;
        expected.bonds[9].bond_type = Single;
        expected.atoms[6].formal_charge = -1;
        expected.atoms[0].formal_charge = 1;
        assert_eq!(valence5n_cleanup8(&mut remote_match, 0), Ok(true));
        assert_eq!(remote_match, expected);
        assert_eq!(remote_match.atoms[10].atomic_number, 50);
    }

    #[test]
    fn source_port__inchi__valence5ncleanup9__line_804() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, ThreeCenter};

        let mut wrong_element = cleanup9_graph(6);
        let before = wrong_element.clone();
        assert_eq!(valence5n_cleanup9(&mut wrong_element, 5), Ok(false));
        assert_eq!(wrong_element, before);
        let mut wrong_charge = cleanup9_graph(7);
        wrong_charge.atoms[5].formal_charge = 1;
        let before = wrong_charge.clone();
        assert_eq!(valence5n_cleanup9(&mut wrong_charge, 5), Ok(false));
        assert_eq!(wrong_charge, before);
        let mut wrong_valence = cleanup9_graph(7);
        wrong_valence.atoms[5].num_explicit_hydrogens = 2;
        let before = wrong_valence.clone();
        assert_eq!(valence5n_cleanup9(&mut wrong_valence, 5), Ok(false));
        assert_eq!(wrong_valence, before);

        let mut valence_error = graph(&[(7, 0), (6, 0)], &[(0, 1, ThreeCenter)]);
        let before = valence_error.clone();
        assert_eq!(
            valence5n_cleanup9(&mut valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(valence_error, before);

        for atom_index in 0_usize..5 {
            let mut no_match = cleanup9_graph(7);
            no_match.atoms[atom_index].atomic_number = 8;
            let before = no_match.clone();
            assert_eq!(valence5n_cleanup9(&mut no_match, 5), Ok(false));
            assert_eq!(no_match, before, "atom predicate {atom_index}");
        }
        for bond_index in 0_usize..5 {
            let mut no_match = cleanup9_graph(7);
            no_match.bonds[bond_index].bond_type = if no_match.bonds[bond_index].bond_type == Single {
                Double
            } else {
                Single
            };
            let before = no_match.clone();
            assert_eq!(valence5n_cleanup9(&mut no_match, 5), Ok(false));
            assert_eq!(no_match, before, "bond predicate {bond_index}");
        }
        let mut root_bond_miss = cleanup9_graph(7);
        root_bond_miss.bonds[5].bond_type = Single;
        root_bond_miss.atoms[5].num_explicit_hydrogens = 4;
        let before = root_bond_miss.clone();
        assert_eq!(valence5n_cleanup9(&mut root_bond_miss, 5), Ok(false));
        assert_eq!(root_bond_miss, before);

        let mut unique = cleanup9_graph(7);
        unique.bonds[0].direction = BeginWedge;
        unique.bonds[1].direction = BeginDash;
        unique.bonds[2].direction = EndDownRight;
        unique.bonds[5].direction = EndUpRight;
        let mut expected = unique.clone();
        expected.bonds[0].bond_type = Double;
        expected.bonds[1].bond_type = Single;
        expected.bonds[5].bond_type = Single;
        expected.atoms[2].formal_charge = -1;
        expected.atoms[5].formal_charge = 1;
        assert_eq!(valence5n_cleanup9(&mut unique, 5), Ok(true));
        assert_eq!(unique, expected);

        let first = cleanup9_graph(7);
        let second = cleanup9_graph(50);
        let multiple = disjoint_union(&first, &second);
        let mut multiple_actual = multiple.clone();
        assert_eq!(valence5n_cleanup9(&mut multiple_actual, 5), Ok(false));
        assert_eq!(multiple_actual, multiple);

        let target_component = graph(
            &[(7, 0), (6, 0), (6, 0), (6, 0), (6, 0)],
            &[(0, 1, Double), (0, 2, Single), (0, 3, Single), (0, 4, Single)],
        );
        let remote = cleanup9_graph(50);
        let mut remote_match = disjoint_union(&target_component, &remote);
        let mut expected = remote_match.clone();
        expected.bonds[4].bond_type = Double;
        expected.bonds[5].bond_type = Single;
        expected.bonds[9].bond_type = Single;
        expected.atoms[7].formal_charge = -1;
        expected.atoms[0].formal_charge = 1;
        assert_eq!(valence5n_cleanup9(&mut remote_match, 0), Ok(true));
        assert_eq!(remote_match, expected);
        assert_eq!(remote_match.atoms[10].atomic_number, 50);
    }

    #[test]
    fn source_port__inchi__valence5ncleanupa__line_855() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, ThreeCenter};

        let mut wrong_element = graph(&[(6, 0), (6, 0)], &[(0, 1, ThreeCenter)]);
        let before = wrong_element.clone();
        assert_eq!(valence5n_cleanupa(&mut wrong_element, 0), Ok(false));
        assert_eq!(wrong_element, before);

        let mut wrong_charge = cleanupa_path_graph(1);
        wrong_charge.atoms[0].formal_charge = 1;
        let before = wrong_charge.clone();
        assert_eq!(valence5n_cleanupa(&mut wrong_charge, 0), Ok(false));
        assert_eq!(wrong_charge, before);

        let mut wrong_valence = cleanupa_path_graph(1);
        wrong_valence.atoms[0].num_explicit_hydrogens = 1;
        let before = wrong_valence.clone();
        assert_eq!(valence5n_cleanupa(&mut wrong_valence, 0), Ok(false));
        assert_eq!(wrong_valence, before);

        let mut valence_error = graph(&[(7, 0), (6, 0)], &[(0, 1, ThreeCenter)]);
        let before = valence_error.clone();
        assert_eq!(
            valence5n_cleanupa(&mut valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(valence_error, before);

        let mut no_match = graph(
            &[(7, 0), (6, 0), (6, 0), (6, 0), (6, 0)],
            &[(0, 1, Double), (0, 2, Single), (0, 3, Single), (0, 4, Single)],
        );
        let before = no_match.clone();
        assert_eq!(valence5n_cleanupa(&mut no_match, 0), Ok(false));
        assert_eq!(no_match, before);

        let mut only_match_contains_root = graph(
            &[(7, 0), (7, 0), (6, 0), (6, 0), (6, 0)],
            &[(1, 0, Double), (0, 2, Single), (0, 3, Single), (0, 4, Single)],
        );
        let before = only_match_contains_root.clone();
        assert_eq!(valence5n_cleanupa(&mut only_match_contains_root, 0), Ok(false));
        assert_eq!(only_match_contains_root, before);

        let root = graph(
            &[(7, 0), (6, 0), (6, 0), (6, 0), (6, 0)],
            &[(0, 1, Double), (0, 2, Single), (0, 3, Single), (0, 4, Single)],
        );
        let pair = graph(&[(7, 0), (7, 0)], &[(1, 0, Double)]);
        let mut no_path = disjoint_union(&root, &pair);
        let before = no_path.clone();
        assert_eq!(valence5n_cleanupa(&mut no_path, 0), Ok(false));
        assert_eq!(no_path, before);

        let mut direct = cleanupa_path_graph(1);
        direct.bonds[0].direction = BeginWedge;
        direct.bonds[1].direction = BeginDash;
        direct.bonds[2].direction = EndDownRight;
        direct.bonds[3].direction = EndUpRight;
        let mut expected = direct.clone();
        expected.bonds[0].bond_type = Single;
        expected.atoms[0].formal_charge = 1;
        assert_eq!(valence5n_cleanupa(&mut direct, 0), Ok(true));
        assert_eq!(direct, expected);
        assert!(direct.atoms.iter().all(|atom| atom.atomic_number != 50));

        let mut at_limit = cleanupa_path_graph(9);
        let mut expected = at_limit.clone();
        for bond in &mut expected.bonds[..9] {
            bond.bond_type = if bond.bond_type == Single { Double } else { Single };
        }
        expected.atoms[0].formal_charge = 1;
        assert_eq!(valence5n_cleanupa(&mut at_limit, 0), Ok(true));
        assert_eq!(at_limit, expected);

        let mut over_limit = cleanupa_path_graph(11);
        let before = over_limit.clone();
        assert_eq!(valence5n_cleanupa(&mut over_limit, 0), Ok(false));
        assert_eq!(over_limit, before);

        let mut shortest = graph(
            &[(7, 0), (6, 0), (6, 0), (7, 0), (7, 0), (7, 0), (7, 0), (6, 0)],
            &[
                (0, 1, Double),
                (1, 2, Single),
                (2, 3, Double),
                (3, 4, Double),
                (0, 5, Double),
                (5, 6, Double),
                (0, 7, Single),
            ],
        );
        let mut expected = shortest.clone();
        expected.bonds[4].bond_type = Single;
        expected.atoms[0].formal_charge = 1;
        assert_eq!(valence5n_cleanupa(&mut shortest, 0), Ok(true));
        assert_eq!(shortest, expected);

        let mut equal_paths = graph(
            &[
                (7, 0),
                (6, 0),
                (6, 0),
                (7, 0),
                (7, 0),
                (6, 0),
                (6, 0),
                (7, 0),
                (7, 0),
                (6, 0),
            ],
            &[
                (0, 1, Double),
                (1, 2, Single),
                (2, 3, Double),
                (3, 4, Double),
                (0, 5, Double),
                (5, 6, Single),
                (6, 7, Double),
                (7, 8, Double),
                (0, 9, Single),
            ],
        );
        let mut expected = equal_paths.clone();
        for bond in &mut expected.bonds[..3] {
            bond.bond_type = if bond.bond_type == Single { Double } else { Single };
        }
        expected.atoms[0].formal_charge = 1;
        assert_eq!(valence5n_cleanupa(&mut equal_paths, 0), Ok(true));
        assert_eq!(equal_paths, expected);
    }

    #[test]
    fn source_port__inchi__valence5ncleanupb__line_916() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight};
        use BondType::{Double, Single, ThreeCenter};

        let mut no_target = graph(&[(7, -4), (6, 0), (8, 0)], &[(0, 1, Single), (0, 2, ThreeCenter)]);
        let before = no_target.clone();
        assert_eq!(valence5n_cleanupb(&mut no_target, 0), Ok(false));
        assert_eq!(no_target, before);

        let mut charged_carbon = graph(&[(16, 8), (6, 1)], &[(0, 1, Double)]);
        let before = charged_carbon.clone();
        assert_eq!(valence5n_cleanupb(&mut charged_carbon, 0), Ok(false));
        assert_eq!(charged_carbon, before);

        let mut over_depth = graph(&[(7, 0), (7, 0), (6, 0)], &[(0, 1, Double), (1, 2, Double)]);
        let before = over_depth.clone();
        assert_eq!(valence5n_cleanupb(&mut over_depth, 0), Ok(false));
        assert_eq!(over_depth, before);

        let mut success = graph(&[(8, -5), (6, 0)], &[(0, 1, Double)]);
        success.atoms[0].num_explicit_hydrogens = 2;
        success.atoms[1].num_explicit_hydrogens = 3;
        success.bonds[0].direction = BeginWedge;
        let mut expected = success.clone();
        expected.atoms[0].formal_charge = 1;
        expected.atoms[1].formal_charge = -1;
        expected.bonds[0].bond_type = Single;
        assert_eq!(valence5n_cleanupb(&mut success, 0), Ok(true));
        assert_eq!(success, expected);

        let mut first_in_adjacency = graph(&[(15, 3), (6, 0), (6, 0)], &[(0, 2, Double), (0, 1, Double)]);
        first_in_adjacency.bonds[0].direction = BeginDash;
        first_in_adjacency.bonds[1].direction = EndDownRight;
        let mut expected = first_in_adjacency.clone();
        expected.atoms[0].formal_charge = 1;
        expected.atoms[2].formal_charge = -1;
        expected.bonds[0].bond_type = Single;
        assert_eq!(valence5n_cleanupb(&mut first_in_adjacency, 0), Ok(true));
        assert_eq!(first_in_adjacency, expected);

        let mut target_valence_error = graph(&[(7, 0), (6, 0), (8, 0)], &[(0, 1, Double), (1, 2, ThreeCenter)]);
        let mut expected = target_valence_error.clone();
        expected.atoms[1].formal_charge = -1;
        assert_eq!(
            valence5n_cleanupb(&mut target_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 1,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(target_valence_error, expected);

        let mut root_valence_error = graph(&[(7, -3), (6, 0), (8, 0)], &[(0, 1, Double), (0, 2, ThreeCenter)]);
        let mut expected = root_valence_error.clone();
        expected.atoms[0].formal_charge = 1;
        expected.atoms[1].formal_charge = -1;
        expected.bonds[0].bond_type = Single;
        assert_eq!(
            valence5n_cleanupb(&mut root_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 1,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(root_valence_error, expected);
    }

    #[test]
    fn source_port__inchi__valence7scleanup1__line_941() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, ThreeCenter};

        let mut wrong_element = graph(&[(15, -1), (8, 0)], &[(0, 1, ThreeCenter)]);
        let before = wrong_element.clone();
        assert_eq!(valence7s_cleanup1(&mut wrong_element, 0), Ok(false));
        assert_eq!(wrong_element, before);

        let mut wrong_charge = graph(&[(16, 0), (8, 0)], &[(0, 1, ThreeCenter)]);
        let before = wrong_charge.clone();
        assert_eq!(valence7s_cleanup1(&mut wrong_charge, 0), Ok(false));
        assert_eq!(wrong_charge, before);

        let mut wrong_valence = graph(&[(16, -1), (8, 0)], &[(0, 1, Double)]);
        let before = wrong_valence.clone();
        assert_eq!(valence7s_cleanup1(&mut wrong_valence, 0), Ok(false));
        assert_eq!(wrong_valence, before);

        let mut entry_valence_error = graph(&[(16, -1), (8, 0)], &[(0, 1, ThreeCenter)]);
        let before = entry_valence_error.clone();
        assert_eq!(
            valence7s_cleanup1(&mut entry_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(entry_valence_error, before);

        for (case_id, neighbor, bond_type, explicit_hydrogens) in [
            ("oxygen-nondouble", 8, Single, 6),
            ("carbon-nonsingle", 6, Double, 5),
            ("other-element", 7, Single, 6),
        ] {
            let mut molecule = graph(&[(16, -1), (neighbor, 0)], &[(0, 1, bond_type)]);
            molecule.atoms[0].num_explicit_hydrogens = explicit_hydrogens;
            let before = molecule.clone();
            assert_eq!(valence7s_cleanup1(&mut molecule, 0), Ok(false), "{case_id}");
            assert_eq!(molecule, before, "{case_id}");
        }

        let mut no_oxygen = graph(&[(16, -1), (6, 0)], &[(0, 1, Single)]);
        no_oxygen.atoms[0].num_explicit_hydrogens = 6;
        let before = no_oxygen.clone();
        assert_eq!(valence7s_cleanup1(&mut no_oxygen, 0), Ok(false));
        assert_eq!(no_oxygen, before);

        let mut criteria_false = graph(&[(16, -1), (8, 0), (8, 0)], &[(0, 1, Double), (0, 2, Double)]);
        criteria_false.atoms[0].num_explicit_hydrogens = 3;
        let before = criteria_false.clone();
        assert_eq!(valence7s_cleanup1(&mut criteria_false, 0), Ok(false));
        assert_eq!(criteria_false, before);

        let mut one_carbon = graph(&[(16, -1), (8, 4), (6, -5)], &[(0, 1, Double), (0, 2, Single)]);
        one_carbon.atoms[0].num_explicit_hydrogens = 4;
        one_carbon.atoms[1].num_explicit_hydrogens = 2;
        one_carbon.bonds[0].direction = BeginWedge;
        one_carbon.bonds[1].direction = BeginDash;
        let mut expected = one_carbon.clone();
        expected.bonds[0].bond_type = Single;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[1].formal_charge = -1;
        assert_eq!(valence7s_cleanup1(&mut one_carbon, 0), Ok(true));
        assert_eq!(one_carbon, expected);

        let mut three_oxygen = graph(
            &[(16, -1), (8, 2), (8, 3), (8, 4)],
            &[(0, 2, Double), (0, 1, Double), (0, 3, Double)],
        );
        three_oxygen.atoms[0].num_explicit_hydrogens = 1;
        three_oxygen.bonds[0].direction = EndDownRight;
        three_oxygen.bonds[1].direction = EndUpRight;
        let mut expected = three_oxygen.clone();
        expected.bonds[2].bond_type = Single;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[3].formal_charge = -1;
        assert_eq!(valence7s_cleanup1(&mut three_oxygen, 0), Ok(true));
        assert_eq!(three_oxygen, expected);

        let mut oxygen_valence_error = graph(
            &[(16, -1), (6, 0), (8, 0), (7, 0)],
            &[(0, 1, Single), (0, 2, Double), (2, 3, ThreeCenter)],
        );
        oxygen_valence_error.atoms[0].num_explicit_hydrogens = 4;
        let mut expected = oxygen_valence_error.clone();
        expected.bonds[1].bond_type = Single;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[2].formal_charge = -1;
        assert_eq!(
            valence7s_cleanup1(&mut oxygen_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 2,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(oxygen_valence_error, expected);
    }

    #[test]
    fn source_port__inchi__valence7scleanup2__line_989() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, ThreeCenter, Triple};

        let mut wrong_element = graph(&[(15, -1), (7, 0)], &[(0, 1, ThreeCenter)]);
        let before = wrong_element.clone();
        assert_eq!(valence7s_cleanup2(&mut wrong_element, 0), Ok(false));
        assert_eq!(wrong_element, before);

        let mut wrong_charge = graph(&[(16, 0), (7, 0)], &[(0, 1, ThreeCenter)]);
        let before = wrong_charge.clone();
        assert_eq!(valence7s_cleanup2(&mut wrong_charge, 0), Ok(false));
        assert_eq!(wrong_charge, before);

        let mut wrong_valence = graph(&[(16, -1), (7, 0)], &[(0, 1, Double)]);
        let before = wrong_valence.clone();
        assert_eq!(valence7s_cleanup2(&mut wrong_valence, 0), Ok(false));
        assert_eq!(wrong_valence, before);

        let mut entry_valence_error = graph(&[(16, -1), (7, 0)], &[(0, 1, ThreeCenter)]);
        let before = entry_valence_error.clone();
        assert_eq!(
            valence7s_cleanup2(&mut entry_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(entry_valence_error, before);

        let mut no_target = graph(&[(16, -1), (6, 0)], &[(0, 1, Double)]);
        no_target.atoms[0].num_explicit_hydrogens = 5;
        let before = no_target.clone();
        assert_eq!(valence7s_cleanup2(&mut no_target, 0), Ok(false));
        assert_eq!(no_target, before);

        let mut over_depth = graph(
            &[(16, -1), (6, 0), (6, 0), (6, 0), (7, 0)],
            &[(0, 1, Double), (1, 2, Single), (2, 3, Double), (3, 4, Triple)],
        );
        over_depth.atoms[0].num_explicit_hydrogens = 5;
        let before = over_depth.clone();
        assert_eq!(valence7s_cleanup2(&mut over_depth, 0), Ok(false));
        assert_eq!(over_depth, before);

        let mut direct = graph(&[(16, -1), (7, 0)], &[(0, 1, Triple)]);
        direct.atoms[0].num_explicit_hydrogens = 4;
        direct.atoms[1].num_explicit_hydrogens = 3;
        direct.bonds[0].direction = BeginWedge;
        let mut expected = direct.clone();
        expected.atoms[0].formal_charge = 0;
        expected.bonds[0].bond_type = Double;
        assert_eq!(valence7s_cleanup2(&mut direct, 0), Ok(true));
        assert_eq!(direct, expected);

        let mut two_bonds = graph(&[(16, -1), (6, -3), (7, 0)], &[(0, 1, Double), (1, 2, Triple)]);
        two_bonds.atoms[0].num_explicit_hydrogens = 5;
        two_bonds.bonds[0].direction = BeginDash;
        two_bonds.bonds[1].direction = EndDownRight;
        let mut expected = two_bonds.clone();
        expected.atoms[0].formal_charge = 0;
        expected.bonds[0].bond_type = Single;
        expected.bonds[1].bond_type = Double;
        assert_eq!(valence7s_cleanup2(&mut two_bonds, 0), Ok(true));
        assert_eq!(two_bonds, expected);

        let mut at_limit = graph(
            &[(16, -1), (6, 2), (6, -4), (7, 0)],
            &[(0, 1, Double), (1, 2, Single), (2, 3, Triple)],
        );
        at_limit.atoms[0].num_explicit_hydrogens = 5;
        at_limit.bonds[0].direction = BeginWedge;
        at_limit.bonds[1].direction = EndUpRight;
        at_limit.bonds[2].direction = BeginDash;
        let mut expected = at_limit.clone();
        expected.atoms[0].formal_charge = 0;
        expected.bonds[0].bond_type = Single;
        expected.bonds[1].bond_type = Double;
        expected.bonds[2].bond_type = Double;
        assert_eq!(valence7s_cleanup2(&mut at_limit, 0), Ok(true));
        assert_eq!(at_limit, expected);

        let mut first_equal_path = graph(&[(16, -1), (7, 0), (7, 0)], &[(0, 2, Triple), (0, 1, Triple)]);
        first_equal_path.atoms[0].num_explicit_hydrogens = 1;
        first_equal_path.bonds[0].direction = EndDownRight;
        first_equal_path.bonds[1].direction = EndUpRight;
        let mut expected = first_equal_path.clone();
        expected.atoms[0].formal_charge = 0;
        expected.bonds[0].bond_type = Double;
        assert_eq!(valence7s_cleanup2(&mut first_equal_path, 0), Ok(true));
        assert_eq!(first_equal_path, expected);
    }

    #[test]
    fn source_port__inchi__valence7scleanup3__line_1021() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, ThreeCenter};

        let mut wrong_element = graph(&[(15, -1), (7, 0)], &[(0, 1, ThreeCenter)]);
        let before = wrong_element.clone();
        assert_eq!(valence7s_cleanup3(&mut wrong_element, 0), Ok(false));
        assert_eq!(wrong_element, before);

        let mut wrong_charge = graph(&[(16, 0), (7, 0)], &[(0, 1, ThreeCenter)]);
        let before = wrong_charge.clone();
        assert_eq!(valence7s_cleanup3(&mut wrong_charge, 0), Ok(false));
        assert_eq!(wrong_charge, before);

        let mut wrong_valence = graph(&[(16, -1), (7, 0)], &[(0, 1, Double)]);
        let before = wrong_valence.clone();
        assert_eq!(valence7s_cleanup3(&mut wrong_valence, 0), Ok(false));
        assert_eq!(wrong_valence, before);

        let mut entry_valence_error = graph(&[(16, -1), (7, 0)], &[(0, 1, ThreeCenter)]);
        let before = entry_valence_error.clone();
        assert_eq!(
            valence7s_cleanup3(&mut entry_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(entry_valence_error, before);

        let mut no_target = graph(&[(16, -1), (6, 0)], &[(0, 1, Double)]);
        no_target.atoms[0].num_explicit_hydrogens = 5;
        let before = no_target.clone();
        assert_eq!(valence7s_cleanup3(&mut no_target, 0), Ok(false));
        assert_eq!(no_target, before);

        let mut charged_target = graph(&[(16, -1), (7, 1)], &[(0, 1, Double)]);
        charged_target.atoms[0].num_explicit_hydrogens = 5;
        let before = charged_target.clone();
        assert_eq!(valence7s_cleanup3(&mut charged_target, 0), Ok(false));
        assert_eq!(charged_target, before);

        let mut over_depth = graph(&[(16, -1), (6, 0), (7, 0)], &[(0, 1, Double), (1, 2, Double)]);
        over_depth.atoms[0].num_explicit_hydrogens = 5;
        let before = over_depth.clone();
        assert_eq!(valence7s_cleanup3(&mut over_depth, 0), Ok(false));
        assert_eq!(over_depth, before);

        let mut success = graph(&[(16, -1), (7, 0)], &[(0, 1, Double)]);
        success.atoms[0].num_explicit_hydrogens = 5;
        success.atoms[1].num_explicit_hydrogens = 3;
        success.bonds[0].direction = BeginWedge;
        let mut expected = success.clone();
        expected.bonds[0].bond_type = Single;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[1].formal_charge = -1;
        assert_eq!(valence7s_cleanup3(&mut success, 0), Ok(true));
        assert_eq!(success, expected);

        let mut first_equal_path = graph(&[(16, -1), (7, 0), (7, 0)], &[(0, 2, Double), (0, 1, Double)]);
        first_equal_path.atoms[0].num_explicit_hydrogens = 3;
        first_equal_path.bonds[0].direction = BeginDash;
        first_equal_path.bonds[1].direction = EndDownRight;
        let mut expected = first_equal_path.clone();
        expected.bonds[0].bond_type = Single;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[2].formal_charge = -1;
        assert_eq!(valence7s_cleanup3(&mut first_equal_path, 0), Ok(true));
        assert_eq!(first_equal_path, expected);

        let mut target_valence_not_recomputed =
            graph(&[(16, -1), (7, 0), (6, -4)], &[(0, 1, Double), (1, 2, ThreeCenter)]);
        target_valence_not_recomputed.atoms[0].num_explicit_hydrogens = 5;
        target_valence_not_recomputed.bonds[0].direction = EndUpRight;
        let mut expected = target_valence_not_recomputed.clone();
        expected.bonds[0].bond_type = Single;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[1].formal_charge = -1;
        assert_eq!(valence7s_cleanup3(&mut target_valence_not_recomputed, 0), Ok(true));
        assert_eq!(target_valence_not_recomputed, expected);
    }

    #[test]
    fn source_port__inchi__valence8scleanup1__line_1044() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight, EndUpRight};
        use BondType::{Double, Single, ThreeCenter};

        let mut wrong_element = graph(&[(15, -1), (7, 0)], &[(0, 1, ThreeCenter)]);
        let before = wrong_element.clone();
        assert_eq!(valence8s_cleanup1(&mut wrong_element, 0), Ok(false));
        assert_eq!(wrong_element, before);

        let mut wrong_charge = graph(&[(16, 0), (7, 0)], &[(0, 1, ThreeCenter)]);
        let before = wrong_charge.clone();
        assert_eq!(valence8s_cleanup1(&mut wrong_charge, 0), Ok(false));
        assert_eq!(wrong_charge, before);

        let mut wrong_valence = graph(&[(16, -1), (7, 0)], &[(0, 1, Double)]);
        let before = wrong_valence.clone();
        assert_eq!(valence8s_cleanup1(&mut wrong_valence, 0), Ok(false));
        assert_eq!(wrong_valence, before);

        let mut entry_valence_error = graph(&[(16, -1), (7, 0)], &[(0, 1, ThreeCenter)]);
        let before = entry_valence_error.clone();
        assert_eq!(
            valence8s_cleanup1(&mut entry_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(entry_valence_error, before);

        let mut no_target = graph(&[(16, -1), (6, 0)], &[(0, 1, Double)]);
        no_target.atoms[0].num_explicit_hydrogens = 5;
        let before = no_target.clone();
        assert_eq!(valence8s_cleanup1(&mut no_target, 0), Ok(false));
        assert_eq!(no_target, before);

        let mut charged_target = graph(&[(16, -1), (7, 1)], &[(0, 1, Double)]);
        charged_target.atoms[0].num_explicit_hydrogens = 5;
        let before = charged_target.clone();
        assert_eq!(valence8s_cleanup1(&mut charged_target, 0), Ok(false));
        assert_eq!(charged_target, before);

        let mut direct = graph(&[(16, -1), (7, 0)], &[(0, 1, Double)]);
        direct.atoms[0].num_explicit_hydrogens = 5;
        direct.atoms[1].num_explicit_hydrogens = 3;
        direct.bonds[0].direction = BeginWedge;
        let mut expected = direct.clone();
        expected.bonds[0].bond_type = Single;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[1].formal_charge = -1;
        expected.atoms[1].num_explicit_hydrogens = 0;
        assert_eq!(valence8s_cleanup1(&mut direct, 0), Ok(true));
        assert_eq!(direct, expected);

        let mut alternating = graph(
            &[(16, -1), (6, 4), (6, -3), (7, 0)],
            &[(0, 1, Double), (1, 2, Single), (2, 3, Double)],
        );
        alternating.atoms[0].num_explicit_hydrogens = 5;
        alternating.atoms[3].num_explicit_hydrogens = 2;
        alternating.bonds[0].direction = BeginDash;
        alternating.bonds[1].direction = EndDownRight;
        alternating.bonds[2].direction = EndUpRight;
        let mut expected = alternating.clone();
        expected.bonds[0].bond_type = Single;
        expected.bonds[1].bond_type = Double;
        expected.bonds[2].bond_type = Single;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[3].formal_charge = -1;
        expected.atoms[3].num_explicit_hydrogens = 0;
        assert_eq!(valence8s_cleanup1(&mut alternating, 0), Ok(true));
        assert_eq!(alternating, expected);

        let mut at_limit = graph(
            &[
                (16, -1),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (7, 0),
            ],
            &[
                (0, 1, Double),
                (1, 2, Single),
                (2, 3, Double),
                (3, 4, Single),
                (4, 5, Double),
                (5, 6, Single),
                (6, 7, Double),
                (7, 8, Single),
                (8, 9, Double),
            ],
        );
        at_limit.atoms[0].num_explicit_hydrogens = 5;
        at_limit.atoms[9].num_explicit_hydrogens = 4;
        for (index, bond) in at_limit.bonds.iter_mut().enumerate() {
            bond.direction = if index % 2 == 0 { EndUpRight } else { EndDownRight };
        }
        let mut expected = at_limit.clone();
        for bond in &mut expected.bonds {
            bond.bond_type = if bond.bond_type == Double { Single } else { Double };
        }
        expected.atoms[0].formal_charge = 0;
        expected.atoms[9].formal_charge = -1;
        expected.atoms[9].num_explicit_hydrogens = 0;
        assert_eq!(valence8s_cleanup1(&mut at_limit, 0), Ok(true));
        assert_eq!(at_limit, expected);

        let mut over_limit = graph(
            &[
                (16, -1),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (6, 0),
                (7, 0),
            ],
            &[
                (0, 1, Double),
                (1, 2, Single),
                (2, 3, Double),
                (3, 4, Single),
                (4, 5, Double),
                (5, 6, Single),
                (6, 7, Double),
                (7, 8, Single),
                (8, 9, Double),
                (9, 10, Single),
                (10, 11, Double),
            ],
        );
        over_limit.atoms[0].num_explicit_hydrogens = 5;
        let before = over_limit.clone();
        assert_eq!(valence8s_cleanup1(&mut over_limit, 0), Ok(false));
        assert_eq!(over_limit, before);

        let mut first_equal_path = graph(&[(16, -1), (7, 0), (7, 0)], &[(0, 2, Double), (0, 1, Double)]);
        first_equal_path.atoms[0].num_explicit_hydrogens = 3;
        first_equal_path.atoms[2].num_explicit_hydrogens = 2;
        first_equal_path.bonds[0].direction = EndUpRight;
        first_equal_path.bonds[1].direction = EndDownRight;
        let mut expected = first_equal_path.clone();
        expected.bonds[0].bond_type = Single;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[2].formal_charge = -1;
        expected.atoms[2].num_explicit_hydrogens = 0;
        assert_eq!(valence8s_cleanup1(&mut first_equal_path, 0), Ok(true));
        assert_eq!(first_equal_path, expected);

        let mut target_valence_error = graph(&[(16, -1), (7, 0), (6, -4)], &[(0, 1, Double), (1, 2, ThreeCenter)]);
        target_valence_error.atoms[0].num_explicit_hydrogens = 5;
        target_valence_error.atoms[1].num_explicit_hydrogens = 4;
        target_valence_error.bonds[0].direction = BeginDash;
        let mut expected = target_valence_error.clone();
        expected.bonds[0].bond_type = Single;
        expected.atoms[1].formal_charge = -1;
        assert_eq!(
            valence8s_cleanup1(&mut target_valence_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 1,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(target_valence_error, expected);
    }

    #[test]
    fn source_port__inchi__valence8clcleanup1__line_1079() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight};
        use BondType::{Double, Single, ThreeCenter};

        let mut wrong_valence = graph(&[(17, -1)], &[]);
        wrong_valence.atoms[0].num_explicit_hydrogens = 7;
        let before = wrong_valence.clone();
        assert_eq!(valence8cl_cleanup1(&mut wrong_valence, 0), Ok(false));
        assert_eq!(wrong_valence, before);

        let mut wrong_charge = graph(&[(17, 0)], &[]);
        wrong_charge.atoms[0].num_explicit_hydrogens = 8;
        let before = wrong_charge.clone();
        assert_eq!(valence8cl_cleanup1(&mut wrong_charge, 0), Ok(false));
        assert_eq!(wrong_charge, before);

        let mut valence_before_charge_error = graph(&[(17, 0), (8, 0)], &[(0, 1, ThreeCenter)]);
        let before = valence_before_charge_error.clone();
        assert_eq!(
            valence8cl_cleanup1(&mut valence_before_charge_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(valence_before_charge_error, before);

        let mut non_oxygen = graph(&[(17, -1), (7, 0)], &[(0, 1, Double)]);
        non_oxygen.atoms[0].num_explicit_hydrogens = 6;
        let before = non_oxygen.clone();
        assert_eq!(valence8cl_cleanup1(&mut non_oxygen, 0), Ok(false));
        assert_eq!(non_oxygen, before);

        let mut empty_neighbors = graph(&[(15, -1)], &[]);
        empty_neighbors.atoms[0].num_explicit_hydrogens = 8;
        assert_eq!(valence8cl_cleanup1(&mut empty_neighbors, 0), Ok(true));
        assert_eq!(empty_neighbors.atoms[0].atomic_number, 15);
        assert_eq!(empty_neighbors.atoms[0].formal_charge, 3);
        assert_eq!(empty_neighbors.atoms[0].num_explicit_hydrogens, 8);

        let mut mixed_bonds = graph(&[(17, -1), (8, 4), (8, 5)], &[(0, 1, Double), (0, 2, Single)]);
        mixed_bonds.atoms[0].num_explicit_hydrogens = 5;
        mixed_bonds.atoms[1].num_explicit_hydrogens = 2;
        mixed_bonds.atoms[2].num_explicit_hydrogens = 3;
        mixed_bonds.bonds[0].direction = BeginWedge;
        mixed_bonds.bonds[1].direction = BeginDash;
        let mut expected = mixed_bonds.clone();
        expected.atoms[0].formal_charge = 3;
        expected.atoms[1].formal_charge = -1;
        expected.bonds[0].bond_type = Single;
        assert_eq!(valence8cl_cleanup1(&mut mixed_bonds, 0), Ok(true));
        assert_eq!(mixed_bonds, expected);

        let mut ordered_partial_error = graph(
            &[(17, -1), (8, 2), (8, 3), (6, 0)],
            &[(0, 1, Double), (0, 2, Double), (2, 3, ThreeCenter)],
        );
        ordered_partial_error.atoms[0].num_explicit_hydrogens = 4;
        ordered_partial_error.bonds[0].direction = EndDownRight;
        ordered_partial_error.bonds[1].direction = BeginDash;
        let mut expected = ordered_partial_error.clone();
        expected.atoms[0].formal_charge = 3;
        expected.atoms[1].formal_charge = -1;
        expected.atoms[2].formal_charge = -1;
        expected.bonds[0].bond_type = Single;
        expected.bonds[1].bond_type = Single;
        assert_eq!(
            valence8cl_cleanup1(&mut ordered_partial_error, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 2,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(ordered_partial_error, expected);
    }

    #[test]
    fn source_port__inchi__valence5clcleanup1__line_1114() {
        use BondDirection::{BeginDash, BeginWedge, EndDownRight};
        use BondType::{Double, Single, ThreeCenter};

        let mut entry_valence_exception = graph(&[(17, 0), (8, -1)], &[(0, 1, ThreeCenter)]);
        let before = entry_valence_exception.clone();
        assert_eq!(
            valence5cl_cleanup1(&mut entry_valence_exception, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(entry_valence_exception, before);

        let mut wrong_valence = graph(&[(17, 1), (8, -1)], &[(0, 1, Single)]);
        wrong_valence.atoms[0].num_explicit_hydrogens = 4;
        let before = wrong_valence.clone();
        assert_eq!(valence5cl_cleanup1(&mut wrong_valence, 0), Ok(false));
        assert_eq!(wrong_valence, before);

        let mut wrong_charge = graph(&[(17, 0), (8, -1)], &[(0, 1, Single)]);
        wrong_charge.atoms[0].num_explicit_hydrogens = 5;
        let before = wrong_charge.clone();
        assert_eq!(valence5cl_cleanup1(&mut wrong_charge, 0), Ok(false));
        assert_eq!(wrong_charge, before);

        let mut wrong_element = graph(&[(17, 1), (7, -1)], &[(0, 1, Single)]);
        wrong_element.atoms[0].num_explicit_hydrogens = 5;
        let before = wrong_element.clone();
        assert_eq!(valence5cl_cleanup1(&mut wrong_element, 0), Ok(false));
        assert_eq!(wrong_element, before);

        let mut wrong_target_charge = graph(&[(17, 1), (8, 0)], &[(0, 1, Single)]);
        wrong_target_charge.atoms[0].num_explicit_hydrogens = 5;
        let before = wrong_target_charge.clone();
        assert_eq!(valence5cl_cleanup1(&mut wrong_target_charge, 0), Ok(false));
        assert_eq!(wrong_target_charge, before);

        let mut wrong_bond = graph(&[(17, 1), (8, -1)], &[(0, 1, Double)]);
        wrong_bond.atoms[0].num_explicit_hydrogens = 4;
        let before = wrong_bond.clone();
        assert_eq!(valence5cl_cleanup1(&mut wrong_bond, 0), Ok(false));
        assert_eq!(wrong_bond, before);

        let mut success_without_element_guard =
            graph(&[(15, 1), (8, -1), (6, 4)], &[(0, 1, Single), (1, 2, ThreeCenter)]);
        success_without_element_guard.atoms[0].num_explicit_hydrogens = 5;
        success_without_element_guard.atoms[1].num_explicit_hydrogens = 3;
        success_without_element_guard.bonds[0].direction = BeginWedge;
        success_without_element_guard.bonds[1].direction = BeginDash;
        let mut expected = success_without_element_guard.clone();
        expected.bonds[0].bond_type = Double;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[1].formal_charge = 0;
        assert_eq!(valence5cl_cleanup1(&mut success_without_element_guard, 0), Ok(true));
        assert_eq!(success_without_element_guard, expected);

        let mut first_equal_target = graph(&[(17, 1), (8, -1), (8, -1)], &[(0, 1, Single), (0, 2, Single)]);
        first_equal_target.atoms[0].num_explicit_hydrogens = 4;
        first_equal_target.atoms[1].num_explicit_hydrogens = 2;
        first_equal_target.atoms[2].num_explicit_hydrogens = 3;
        first_equal_target.bonds[0].direction = EndDownRight;
        first_equal_target.bonds[1].direction = BeginDash;
        let mut expected = first_equal_target.clone();
        expected.bonds[0].bond_type = Double;
        expected.atoms[0].formal_charge = 0;
        expected.atoms[1].formal_charge = 0;
        assert_eq!(valence5cl_cleanup1(&mut first_equal_target, 0), Ok(true));
        assert_eq!(first_equal_target, expected);
    }

    #[test]
    fn source_port__inchi__valence3clcleanup1__line_1134() {
        use BondDirection::{BeginDash, BeginWedge};
        use BondType::{Double, Single, ThreeCenter, Triple};

        let mut entry_valence_exception = graph(&[(17, 1), (16, 0)], &[(0, 1, ThreeCenter)]);
        let before = entry_valence_exception.clone();
        assert_eq!(
            valence3cl_cleanup1(&mut entry_valence_exception, 0),
            Err(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            })
        );
        assert_eq!(entry_valence_exception, before);

        let mut wrong_valence = graph(&[(17, 0), (16, 0)], &[(0, 1, Double)]);
        let before = wrong_valence.clone();
        assert_eq!(valence3cl_cleanup1(&mut wrong_valence, 0), Ok(false));
        assert_eq!(wrong_valence, before);

        let mut wrong_charge = graph(&[(17, 1), (16, 0)], &[(0, 1, Triple)]);
        let before = wrong_charge.clone();
        assert_eq!(valence3cl_cleanup1(&mut wrong_charge, 0), Ok(false));
        assert_eq!(wrong_charge, before);

        let mut wrong_element = graph(&[(17, 0), (8, 0)], &[(0, 1, Triple)]);
        let before = wrong_element.clone();
        assert_eq!(valence3cl_cleanup1(&mut wrong_element, 0), Ok(false));
        assert_eq!(wrong_element, before);

        let mut wrong_target_charge = graph(&[(17, 0), (16, 1)], &[(0, 1, Triple)]);
        let before = wrong_target_charge.clone();
        assert_eq!(valence3cl_cleanup1(&mut wrong_target_charge, 0), Ok(false));
        assert_eq!(wrong_target_charge, before);

        let mut wrong_bond = graph(&[(17, 0), (16, 0)], &[(0, 1, Double)]);
        wrong_bond.atoms[0].num_explicit_hydrogens = 1;
        let before = wrong_bond.clone();
        assert_eq!(valence3cl_cleanup1(&mut wrong_bond, 0), Ok(false));
        assert_eq!(wrong_bond, before);

        let mut success_without_element_guard =
            graph(&[(15, 0), (16, 0), (6, 4)], &[(0, 1, Triple), (1, 2, ThreeCenter)]);
        success_without_element_guard.atoms[1].num_explicit_hydrogens = 3;
        success_without_element_guard.bonds[0].direction = BeginWedge;
        success_without_element_guard.bonds[1].direction = BeginDash;
        let mut expected = success_without_element_guard.clone();
        expected.bonds[0].bond_type = Single;
        assert_eq!(valence3cl_cleanup1(&mut success_without_element_guard, 0), Ok(true));
        assert_eq!(success_without_element_guard, expected);
    }

    #[test]
    fn source_port__inchi__cleanup__line_1151() {
        use BondDirection::{BeginDash, BeginWedge};
        use BondType::{Double, Single, ThreeCenter, Triple};

        let mut empty = graph(&[], &[]);
        assert_eq!(clean_up(&mut empty), Ok(()));

        let mut ignored = graph(&[(6, 2)], &[]);
        ignored.atoms[0].is_aromatic = true;
        let before = ignored.clone();
        assert_eq!(clean_up(&mut ignored), Ok(()));
        assert_eq!(ignored, before);

        let mut valence_error = graph(&[(7, 0), (6, 0)], &[(0, 1, ThreeCenter)]);
        valence_error.atoms[0].is_aromatic = true;
        let before = valence_error.clone();
        assert_eq!(
            clean_up(&mut valence_error),
            Err(AdapterCleanup5Error::Valence(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 0,
                bond_type: ThreeCenter,
            }))
        );
        assert_eq!(valence_error, before);

        let mut valence4_cleanup1 = graph(
            &[(6, 4), (7, 2), (7, -1), (7, -3), (7, 1), (8, 0)],
            &[
                (1, 0, Single),
                (2, 1, Double),
                (3, 2, Double),
                (4, 3, Single),
                (0, 4, Double),
                (0, 5, Single),
            ],
        );
        valence4_cleanup1.atoms[2].is_aromatic = true;
        valence4_cleanup1.bonds[0].direction = BeginDash;
        let mut expected = valence4_cleanup1.clone();
        for (bond_index, bond_type) in [Double, Single, Single, Double, Single].into_iter().enumerate() {
            expected.bonds[bond_index].bond_type = bond_type;
        }
        assert_eq!(clean_up(&mut valence4_cleanup1), Ok(()));
        assert_eq!(valence4_cleanup1, expected);

        let mut valence4_cleanup2 = graph(&[(7, -1), (7, 0), (8, -2)], &[(1, 0, Double), (0, 2, Single)]);
        valence4_cleanup2.atoms[0].num_explicit_hydrogens = 1;
        valence4_cleanup2.atoms[0].is_aromatic = true;
        let mut expected = valence4_cleanup2.clone();
        expected.atoms[0].formal_charge = 0;
        expected.atoms[1].formal_charge = -1;
        expected.bonds[0].bond_type = Single;
        assert_eq!(clean_up(&mut valence4_cleanup2), Ok(()));
        assert_eq!(valence4_cleanup2, expected);

        let mut valence4_no_match = graph(&[(7, 0)], &[]);
        valence4_no_match.atoms[0].num_explicit_hydrogens = 4;
        valence4_no_match.atoms[0].is_aromatic = true;
        let before = valence4_no_match.clone();
        assert_eq!(clean_up(&mut valence4_no_match), Ok(()));
        assert_eq!(valence4_no_match, before);

        let mut charged_nitrogen = graph(&[(7, 1)], &[]);
        charged_nitrogen.atoms[0].is_aromatic = true;
        let before = charged_nitrogen.clone();
        assert_eq!(clean_up(&mut charged_nitrogen), Ok(()));
        assert_eq!(charged_nitrogen, before);

        let mut all_valence5_miss = graph(&[(7, 0)], &[]);
        all_valence5_miss.atoms[0].num_explicit_hydrogens = 5;
        all_valence5_miss.atoms[0].is_aromatic = true;
        let before = all_valence5_miss.clone();
        assert_eq!(clean_up(&mut all_valence5_miss), Ok(()));
        assert_eq!(all_valence5_miss, before);

        let mut cleanup6 = cleanup6_graph(7, Triple);
        cleanup6.atoms[2].is_aromatic = true;
        assert_eq!(clean_up(&mut cleanup6), Ok(()));
        assert_eq!(cleanup6.atoms[2].formal_charge, 1);
        assert!(cleanup6.atoms[2].is_aromatic);
        assert_eq!(cleanup6.bonds[0].bond_type, Double);
        assert_eq!(cleanup6.bonds[1].bond_type, Single);
        assert_eq!(cleanup6.bonds[4].bond_type, Double);
        assert_eq!(cleanup6.bonds[5].bond_type, Single);

        let mut cleanup7 = cleanup7_graph(7, Triple);
        cleanup7.atoms[2].is_aromatic = true;
        assert_eq!(clean_up(&mut cleanup7), Ok(()));
        assert_eq!(cleanup7.atoms[2].formal_charge, 0);
        assert_eq!(cleanup7.atoms[7].formal_charge, -1);
        assert_eq!(cleanup7.bonds[1].bond_type, Single);
        assert_eq!(cleanup7.bonds[2].bond_type, Single);
        assert_eq!(cleanup7.bonds[3].bond_type, Double);
        assert_eq!(cleanup7.bonds[7].bond_type, Single);
        assert!(cleanup7.atoms[2].is_aromatic);

        let mut cleanup8 = cleanup8_graph(7);
        cleanup8.atoms[5].is_aromatic = true;
        assert_eq!(clean_up(&mut cleanup8), Ok(()));
        assert_eq!(cleanup8.atoms[5].formal_charge, 1);
        assert!(cleanup8.atoms[5].is_aromatic);

        let mut cleanup9 = cleanup9_graph(7);
        cleanup9.atoms[5].is_aromatic = true;
        assert_eq!(clean_up(&mut cleanup9), Ok(()));
        assert_eq!(cleanup9.atoms[5].formal_charge, 1);
        assert!(cleanup9.atoms[5].is_aromatic);

        let mut cleanupa = cleanupa_path_graph(1);
        cleanupa.atoms[0].is_aromatic = true;
        assert_eq!(clean_up(&mut cleanupa), Ok(()));
        assert_eq!(cleanupa.atoms[0].formal_charge, 1);
        assert!(cleanupa.atoms[0].is_aromatic);
        assert_eq!(cleanupa.bonds[0].bond_type, Single);

        let mut cleanup1 = graph(&[(7, 0), (7, 1)], &[(0, 1, Double)]);
        cleanup1.atoms[0].num_explicit_hydrogens = 3;
        cleanup1.atoms[0].is_aromatic = true;
        assert_eq!(clean_up(&mut cleanup1), Ok(()));
        assert_eq!(cleanup1.atoms[0].formal_charge, 1);
        assert_eq!(cleanup1.atoms[1].formal_charge, 0);
        assert_eq!(cleanup1.bonds[0].bond_type, Single);
        assert!(cleanup1.atoms[0].is_aromatic);

        let mut cleanup2 = graph(&[(7, 0), (6, 0), (7, -1)], &[(0, 1, Triple), (1, 2, Single)]);
        cleanup2.atoms[0].num_explicit_hydrogens = 2;
        cleanup2.atoms[0].is_aromatic = true;
        assert_eq!(clean_up(&mut cleanup2), Ok(()));
        assert_eq!(cleanup2.atoms[1].formal_charge, -1);
        assert_eq!(cleanup2.atoms[2].formal_charge, 0);
        assert_eq!(cleanup2.bonds[0].bond_type, Single);
        assert_eq!(cleanup2.bonds[1].bond_type, Double);
        assert!(cleanup2.atoms[0].is_aromatic);

        let mut cleanup3 = graph(&[(7, 0), (7, 0)], &[(0, 1, Double)]);
        cleanup3.atoms[0].num_explicit_hydrogens = 3;
        cleanup3.atoms[0].is_aromatic = true;
        assert_eq!(clean_up(&mut cleanup3), Ok(()));
        assert_eq!(cleanup3.atoms[0].formal_charge, 1);
        assert_eq!(cleanup3.atoms[1].formal_charge, -1);
        assert_eq!(cleanup3.bonds[0].bond_type, Single);
        assert!(cleanup3.atoms[0].is_aromatic);

        let mut cleanup4 = graph(&[(7, 0), (14, -1), (14, -1)], &[(0, 1, Double), (0, 2, Double)]);
        cleanup4.atoms[0].num_explicit_hydrogens = 1;
        cleanup4.atoms[0].is_aromatic = true;
        assert_eq!(clean_up(&mut cleanup4), Ok(()));
        assert_eq!(cleanup4.atoms[1].formal_charge, 0);
        assert_eq!(cleanup4.atoms[2].formal_charge, 0);
        assert_eq!(cleanup4.bonds[0].bond_type, Single);
        assert_eq!(cleanup4.bonds[1].bond_type, Single);
        assert!(cleanup4.atoms[0].is_aromatic);

        for target_atomic_number in [8, 16, 9, 17] {
            let mut cleanup5 = graph(&[(7, 0), (target_atomic_number, 0)], &[(0, 1, Double)]);
            cleanup5.atoms[0].num_explicit_hydrogens = 3;
            cleanup5.atoms[0].is_aromatic = true;
            assert_eq!(clean_up(&mut cleanup5), Ok(()));
            assert_eq!(cleanup5.atoms[0].formal_charge, 1);
            assert_eq!(cleanup5.atoms[1].formal_charge, -1);
            assert_eq!(cleanup5.bonds[0].bond_type, Single);
            assert!(cleanup5.atoms[0].is_aromatic);
        }

        let mut cleanupb = graph(&[(7, 0), (6, 0)], &[(0, 1, Double)]);
        cleanupb.atoms[0].num_explicit_hydrogens = 3;
        cleanupb.atoms[0].is_aromatic = true;
        assert_eq!(clean_up(&mut cleanupb), Ok(()));
        assert_eq!(cleanupb.atoms[0].formal_charge, 1);
        assert_eq!(cleanupb.atoms[1].formal_charge, -1);
        assert_eq!(cleanupb.bonds[0].bond_type, Single);
        assert!(cleanupb.atoms[0].is_aromatic);

        let mut partial_error = graph(&[(7, 0), (8, 0), (6, 0)], &[(0, 1, Double), (1, 2, ThreeCenter)]);
        partial_error.atoms[0].num_explicit_hydrogens = 3;
        partial_error.atoms[0].is_aromatic = true;
        assert_eq!(
            clean_up(&mut partial_error),
            Err(AdapterCleanup5Error::Valence(AdapterValenceError {
                kind: "ValueErrorException",
                message: "Bad bond type",
                bond_index: 1,
                bond_type: ThreeCenter,
            }))
        );
        assert!(!partial_error.atoms[0].is_aromatic);
        assert_eq!(partial_error.atoms[0].formal_charge, 1);
        assert_eq!(partial_error.atoms[1].formal_charge, -1);
        assert_eq!(partial_error.bonds[0].bond_type, Single);

        let mut chlorine8 = graph(&[(17, -1)], &[]);
        chlorine8.atoms[0].num_explicit_hydrogens = 8;
        assert_eq!(clean_up(&mut chlorine8), Ok(()));
        assert_eq!(chlorine8.atoms[0].formal_charge, 3);

        let mut chlorine5_dispatch_mismatch = graph(&[(17, 1)], &[]);
        chlorine5_dispatch_mismatch.atoms[0].num_explicit_hydrogens = 5;
        let before = chlorine5_dispatch_mismatch.clone();
        assert_eq!(clean_up(&mut chlorine5_dispatch_mismatch), Ok(()));
        assert_eq!(chlorine5_dispatch_mismatch, before);

        let mut chlorine3 = graph(&[(17, 0), (16, 0)], &[(0, 1, Triple)]);
        chlorine3.bonds[0].direction = BeginWedge;
        assert_eq!(clean_up(&mut chlorine3), Ok(()));
        assert_eq!(chlorine3.bonds[0].bond_type, Single);
        assert_eq!(chlorine3.bonds[0].direction, BeginWedge);

        let mut sulfur7_cleanup1 = graph(&[(16, -1), (8, 4), (6, -5)], &[(0, 1, Double), (0, 2, Single)]);
        sulfur7_cleanup1.atoms[0].num_explicit_hydrogens = 4;
        sulfur7_cleanup1.atoms[1].num_explicit_hydrogens = 2;
        sulfur7_cleanup1.bonds[0].direction = BeginWedge;
        assert_eq!(clean_up(&mut sulfur7_cleanup1), Ok(()));
        assert_eq!(sulfur7_cleanup1.atoms[0].formal_charge, 0);
        assert_eq!(sulfur7_cleanup1.atoms[1].formal_charge, -1);
        assert_eq!(sulfur7_cleanup1.bonds[0].bond_type, Single);

        let mut sulfur7_cleanup2 = graph(&[(16, -1), (7, 0)], &[(0, 1, Triple)]);
        sulfur7_cleanup2.atoms[0].num_explicit_hydrogens = 4;
        assert_eq!(clean_up(&mut sulfur7_cleanup2), Ok(()));
        assert_eq!(sulfur7_cleanup2.atoms[0].formal_charge, 0);
        assert_eq!(sulfur7_cleanup2.bonds[0].bond_type, Double);

        let mut sulfur7_cleanup3 = graph(&[(16, -1), (7, 0)], &[(0, 1, Double)]);
        sulfur7_cleanup3.atoms[0].num_explicit_hydrogens = 5;
        sulfur7_cleanup3.atoms[1].num_explicit_hydrogens = 3;
        assert_eq!(clean_up(&mut sulfur7_cleanup3), Ok(()));
        assert_eq!(sulfur7_cleanup3.atoms[0].formal_charge, 0);
        assert_eq!(sulfur7_cleanup3.atoms[1].formal_charge, -1);
        assert_eq!(sulfur7_cleanup3.bonds[0].bond_type, Single);

        let mut sulfur7_all_miss = graph(&[(16, -1)], &[]);
        sulfur7_all_miss.atoms[0].num_explicit_hydrogens = 7;
        let before = sulfur7_all_miss.clone();
        assert_eq!(clean_up(&mut sulfur7_all_miss), Ok(()));
        assert_eq!(sulfur7_all_miss, before);

        let mut sulfur8 = graph(&[(16, -1), (7, 0)], &[(0, 1, Double)]);
        sulfur8.atoms[0].num_explicit_hydrogens = 6;
        sulfur8.atoms[1].num_explicit_hydrogens = 3;
        let before = sulfur8.clone();
        assert_eq!(clean_up(&mut sulfur8), Ok(()));
        assert_eq!(sulfur8, before);

        let mut bromine_selenium = graph(&[(35, 0), (34, 0)], &[(0, 1, Triple)]);
        bromine_selenium.bonds[0].direction = BeginDash;
        assert_eq!(clean_up(&mut bromine_selenium), Ok(()));
        assert_eq!(bromine_selenium.bonds[0].bond_type, Single);
        assert_eq!(bromine_selenium.bonds[0].direction, BeginDash);

        for mut bromine_miss in [
            graph(&[(35, 1), (34, 0)], &[(0, 1, Triple)]),
            graph(&[(35, 0), (34, 0), (6, 0)], &[(0, 1, Double), (0, 2, Single)]),
            graph(&[(35, 0), (6, 0)], &[(0, 1, Triple)]),
        ] {
            let before = bromine_miss.clone();
            assert_eq!(clean_up(&mut bromine_miss), Ok(()));
            assert_eq!(bromine_miss, before);
        }
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanup3__exact() {
        assert_valence_cleanup_oracle(
            "--valence5n-cleanup3-records",
            "bf11af3caae7ff7d1dd7b381defd449e6cb9b193e9dab328a35a6ddd5d73e25d",
            "_Valence5NCleanUp3",
            &[
                "charged-oxygen-does-not-guard",
                "no-neutral-nitrogen",
                "oxygen-guard-true-no-mutation",
                "root-valence-exception-after-bond-and-charges",
                "target-valence-exception-after-charge",
            ],
            valence5n_cleanup3,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanup4__exact() {
        assert_valence_cleanup_oracle(
            "--valence5n-cleanup4-records",
            "701347673195ee11dec11ba3013f6eb10746d62e0a879ce63d0de98ec0642746",
            "_Valence5NCleanUp4",
            &[
                "exactly-two-reversed-preserves-fields",
                "exactly-two-with-late-nonmatch",
                "no-match-all-predicates",
                "one-match",
                "three-matches-early-return",
            ],
            |molecule, atom_index| Ok(valence5n_cleanup4(molecule, atom_index)),
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanup5__exact() {
        assert_valence_cleanup5_oracle();
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanup6__exact() {
        let mut case_ids = vec![
            "wrong-element-short-circuit".to_owned(),
            "wrong-charge-short-circuit".to_owned(),
            "wrong-valence".to_owned(),
            "unsupported-valence-bond-exception".to_owned(),
            "zero-substructure-matches".to_owned(),
            "atom-predicate-miss-0".to_owned(),
            "atom-predicate-miss-1".to_owned(),
            "atom-predicate-miss-3".to_owned(),
            "atom-predicate-miss-4".to_owned(),
            "atom-predicate-miss-5".to_owned(),
            "atom-predicate-miss-6".to_owned(),
            "nonroot-bond-predicate-miss-0".to_owned(),
            "nonroot-bond-predicate-miss-4".to_owned(),
            "nonroot-bond-predicate-miss-5".to_owned(),
            "root-bonds-triple-single-single".to_owned(),
            "root-bonds-double-single-double".to_owned(),
        ];
        case_ids.extend((0..=21).map(|bond_type| format!("unspecified-bridge-type-{bond_type}")));
        case_ids.extend([
            "multiple-unique-atom-set-matches".to_owned(),
            "unique-match-does-not-contain-argument-atom".to_owned(),
        ]);
        let case_id_refs = case_ids.iter().map(String::as_str).collect::<Vec<_>>();
        assert_valence_cleanup_oracle(
            "--valence5n-cleanup6-records",
            "33f731b27422f52bbe7b9304df66f3126a325a3458a577bcf9b29c01d4b9f47d",
            "_Valence5NCleanUp6",
            &case_id_refs,
            valence5n_cleanup6,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanup7__exact() {
        let mut case_ids = vec![
            "no-target-before-preconditions".to_owned(),
            "wrong-element-after-search".to_owned(),
            "wrong-charge-after-search".to_owned(),
            "wrong-valence-after-search".to_owned(),
            "unsupported-valence-bond-exception-after-search".to_owned(),
            "zero-substructure-matches".to_owned(),
        ];
        case_ids.extend([0_u32, 1, 3, 4, 5, 6].map(|atom_index| format!("atom-predicate-miss-{atom_index}")));
        case_ids.extend([4_u32, 5].map(|bond_index| format!("nonpath-bond-predicate-miss-{bond_index}")));
        case_ids.extend((0..=21).map(|bond_type| format!("unspecified-bridge-type-{bond_type}")));
        case_ids.extend([
            "depth-five-nontarget-stops-search".to_owned(),
            "multiple-unique-atom-set-matches".to_owned(),
            "unique-match-does-not-contain-argument-or-target".to_owned(),
        ]);
        let case_id_refs = case_ids.iter().map(String::as_str).collect::<Vec<_>>();
        assert_valence_cleanup_oracle(
            "--valence5n-cleanup7-records",
            "52035cacfe7ced70b44df5ac007a41f101a968fbee2a5338abf5ff81a4efa562",
            "_Valence5NCleanUp7",
            &case_id_refs,
            valence5n_cleanup7,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanup8__exact() {
        let mut case_ids = vec![
            "wrong-element".to_owned(),
            "wrong-charge".to_owned(),
            "wrong-valence".to_owned(),
            "unsupported-valence-bond-exception".to_owned(),
        ];
        case_ids.extend((0..5).map(|atom_index| format!("atom-predicate-miss-{atom_index}")));
        case_ids.extend((0..5).map(|bond_index| format!("bond-predicate-miss-{bond_index}")));
        case_ids.extend([
            "root-bond-predicate-miss".to_owned(),
            "unique-match".to_owned(),
            "multiple-unique-atom-set-matches".to_owned(),
            "unique-match-does-not-contain-argument-atom".to_owned(),
        ]);
        let case_id_refs = case_ids.iter().map(String::as_str).collect::<Vec<_>>();
        assert_valence_cleanup_oracle(
            "--valence5n-cleanup8-records",
            "9874b6974fcf056f0a7313671fcb32d818b00e9d66c7e11934a5919569440c65",
            "_Valence5NCleanUp8",
            &case_id_refs,
            valence5n_cleanup8,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanup9__exact() {
        let mut case_ids = vec![
            "wrong-element".to_owned(),
            "wrong-charge".to_owned(),
            "wrong-valence".to_owned(),
            "unsupported-valence-bond-exception".to_owned(),
        ];
        case_ids.extend((0..5).map(|atom_index| format!("atom-predicate-miss-{atom_index}")));
        case_ids.extend((0..5).map(|bond_index| format!("bond-predicate-miss-{bond_index}")));
        case_ids.extend([
            "root-bond-predicate-miss".to_owned(),
            "unique-match".to_owned(),
            "multiple-unique-atom-set-matches".to_owned(),
            "unique-match-does-not-contain-argument-atom".to_owned(),
        ]);
        let case_id_refs = case_ids.iter().map(String::as_str).collect::<Vec<_>>();
        assert_valence_cleanup_oracle(
            "--valence5n-cleanup9-records",
            "c2b7f87dfdd510b485982edfc6a4f0327ca042e86b5e35eee7484e2193096ef0",
            "_Valence5NCleanUp9",
            &case_id_refs,
            valence5n_cleanup9,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanupa__exact() {
        assert_valence_cleanup_oracle(
            "--valence5n-cleanupa-records",
            "2d4518bea1a8abb8d5736f1c784fa7e30cb033950dcd6ba5da13846bb46d6392",
            "_Valence5NCleanUpA",
            &[
                "wrong-element",
                "wrong-charge",
                "wrong-valence",
                "unsupported-valence-bond-exception",
                "zero-substructure-matches",
                "only-match-contains-root",
                "remote-match-without-path",
                "direct-path",
                "path-length-9",
                "path-length-11",
                "later-shorter-path",
                "equal-path-retains-first",
            ],
            valence5n_cleanupa,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanupb__exact() {
        assert_valence_cleanup_oracle(
            "--valence5n-cleanupb-records",
            "9d097bf179047ec55e31d8f45ffe573f449e3d5828833e5430d2faf355a5defc",
            "_Valence5NCleanUpB",
            &[
                "no-target",
                "charged-carbon",
                "over-depth",
                "success",
                "first-in-adjacency",
                "target-valence-exception",
                "root-valence-exception",
            ],
            valence5n_cleanupb,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence7scleanup1__exact() {
        assert_valence_cleanup_oracle(
            "--valence7s-cleanup1-records",
            "a03efbad1fb09fa47513380a625fefccda34d23e3d4f5fd0b2f2b2cf3f8d076e",
            "_Valence7SCleanUp1",
            &[
                "wrong-element",
                "wrong-charge",
                "wrong-valence",
                "entry-valence-exception",
                "oxygen-nondouble",
                "carbon-nonsingle",
                "other-element",
                "no-oxygen",
                "criteria-false",
                "one-carbon",
                "three-oxygen",
                "oxygen-valence-exception",
            ],
            valence7s_cleanup1,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence7scleanup2__exact() {
        assert_valence_cleanup_oracle(
            "--valence7s-cleanup2-records",
            "d8e8ca7577874896efc2f1df36738d0bee34a8dcabb35493be0095017a2e9618",
            "_Valence7SCleanUp2",
            &[
                "wrong-element",
                "wrong-charge",
                "wrong-valence",
                "entry-valence-exception",
                "no-target",
                "over-depth",
                "direct",
                "two-bonds",
                "at-limit",
                "first-equal-path",
            ],
            valence7s_cleanup2,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence7scleanup3__exact() {
        assert_valence_cleanup_oracle(
            "--valence7s-cleanup3-records",
            "b61785dfecb61f3b567183f96ec1290ef56f576436bfcd07a808907ebbf4b4ed",
            "_Valence7SCleanUp3",
            &[
                "wrong-element",
                "wrong-charge",
                "wrong-valence",
                "entry-valence-exception",
                "no-target",
                "charged-target",
                "over-depth",
                "success",
                "first-equal-path",
                "target-valence-not-recomputed",
            ],
            valence7s_cleanup3,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence8scleanup1__exact() {
        assert_valence_cleanup_oracle(
            "--valence8s-cleanup1-records",
            "e8bda5f9c53e83b3f5a0cce7f26861384606e709e5ef66e7adedac5c0ebcfd78",
            "_Valence8SCleanUp1",
            &[
                "wrong-element",
                "wrong-charge",
                "wrong-valence",
                "entry-valence-exception",
                "no-target",
                "charged-target",
                "direct",
                "alternating",
                "at-limit",
                "over-limit",
                "first-equal-path",
                "target-valence-exception",
            ],
            valence8s_cleanup1,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence8clcleanup1__exact() {
        assert_valence_cleanup_oracle(
            "--valence8cl-cleanup1-records",
            "8e97a538eab77c78c2fb02be47d47cf6be989761cc280cf4857b9500885b1f3a",
            "_Valence8ClCleanUp1",
            &[
                "wrong-valence",
                "wrong-charge",
                "valence-before-charge-exception",
                "non-oxygen",
                "empty-neighbors-no-element-guard",
                "mixed-bonds",
                "ordered-partial-error",
            ],
            valence8cl_cleanup1,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5clcleanup1__exact() {
        assert_valence_cleanup_oracle(
            "--valence5cl-cleanup1-records",
            "102c9d4a05c61706c9b0de30d61fc725fc394f80d532810d0b94b0103d6ff250",
            "_Valence5ClCleanUp1",
            &[
                "entry-valence-exception",
                "wrong-valence",
                "wrong-charge",
                "wrong-element",
                "wrong-target-charge",
                "wrong-bond",
                "success-no-element-guard-target-not-recomputed",
                "first-equal-target",
            ],
            valence5cl_cleanup1,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence3clcleanup1__exact() {
        assert_valence_cleanup_oracle(
            "--valence3cl-cleanup1-records",
            "078efac5f5acece4f6323b171f6aab2b1795873b55a505e56fb44aa54fa7d982",
            "_Valence3ClCleanUp1",
            &[
                "entry-valence-exception",
                "wrong-valence",
                "wrong-charge",
                "wrong-element",
                "wrong-target-charge",
                "wrong-bond",
                "success-no-element-guard-target-not-recomputed",
            ],
            valence3cl_cleanup1,
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__cleanup__exact() {
        assert_clean_up_oracle();
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanup1__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--valence5n-cleanup1-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let mut case_ids = BTreeSet::new();

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "58b4273fc255b5ed398a0d81b30069e58d68478d083b499f41a77ec508250f06"
            );
            assert_eq!(official["operation"], "_Valence5NCleanUp1");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let mut molecule = AdapterMol::from_graph(
                oracle_atoms(&official["input"]["atom_fields"], case_id),
                oracle_bonds(&official["input"]["bond_fields"], case_id),
            );
            let molecule_before = molecule.clone();
            let atom_index = oracle_u32(&official["input"]["atom_index"], "atom_index", case_id);
            let rust = valence5n_cleanup1(&mut molecule, atom_index);

            match official["output"]["status"].as_str() {
                Some("return") => {
                    assert_eq!(
                        rust,
                        Ok(official["output"]["result"]
                            .as_bool()
                            .unwrap_or_else(|| panic!("{case_id}: result must be bool"))),
                        "{case_id}"
                    );
                    assert!(official["output"]["exception"].is_null(), "{case_id}");
                    assert_eq!(official["output"]["diagnostics"], Value::Array(vec![]));
                }
                Some("exception") => {
                    let exception = &official["output"]["exception"];
                    let expected = AdapterValenceError {
                        kind: "ValueErrorException",
                        message: "Bad bond type",
                        bond_index: oracle_u32(&exception["bond_index"], "bond_index", case_id),
                        bond_type: oracle_bond_type(&exception["bond_type"], case_id),
                    };
                    assert_eq!(exception["kind"], expected.kind, "{case_id}");
                    assert_eq!(exception["message"], expected.message, "{case_id}");
                    assert_eq!(rust, Err(expected.clone()), "{case_id}");
                    assert_eq!(official["output"]["result"], Value::Null, "{case_id}");
                    let diagnostics = official["output"]["diagnostics"]
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: diagnostics must be array"));
                    assert_eq!(diagnostics.len(), 1, "{case_id}");
                    assert_eq!(diagnostics[0]["kind"], expected.kind, "{case_id}");
                    assert_eq!(diagnostics[0]["message"], expected.message, "{case_id}");
                    assert_eq!(
                        oracle_u32(&diagnostics[0]["bond_index"], "bond_index", case_id),
                        expected.bond_index,
                        "{case_id}"
                    );
                    assert_eq!(
                        oracle_bond_type(&diagnostics[0]["bond_type"], case_id),
                        expected.bond_type,
                        "{case_id}"
                    );
                }
                status => panic!("{case_id}: unexpected status {status:?}"),
            }

            assert_eq!(
                official["output"]["graph_unchanged"].as_bool(),
                Some(molecule == molecule_before),
                "{case_id}"
            );
            assert_eq!(
                oracle_u32(&official["output"]["atom_count"], "atom_count", case_id),
                molecule.atoms.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                oracle_u32(&official["output"]["bond_count"], "bond_count", case_id),
                molecule.bonds.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                oracle_atoms(&official["output"]["atom_fields"], case_id),
                molecule.atoms,
                "{case_id}"
            );
            assert_eq!(
                oracle_bonds(&official["output"]["bond_fields"], case_id),
                molecule.bonds,
                "{case_id}"
            );
            assert_eq!(
                official["output"]["stereo_fields"]["bond_directions"],
                Value::Array(
                    molecule
                        .bonds
                        .iter()
                        .map(|bond| Value::from(bond.direction as u8))
                        .collect()
                ),
                "{case_id}"
            );
            assert_eq!(official["output"]["properties"], Value::Array(vec![]));
        }

        assert_eq!(
            case_ids,
            BTreeSet::from([
                "alternating-path".to_owned(),
                "direct-reversed-preserves-fields".to_owned(),
                "no-target".to_owned(),
                "over-depth-limit".to_owned(),
                "valence-exception-partial-mutation".to_owned(),
            ])
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence5ncleanup2__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--valence5n-cleanup2-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let mut case_ids = BTreeSet::new();

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "07654fddc42017f1aa76c39ec11013301f92e460b199218c6a2f3bf6bcfb4484"
            );
            assert_eq!(official["operation"], "_Valence5NCleanUp2");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let mut molecule = AdapterMol::from_graph(
                oracle_atoms(&official["input"]["atom_fields"], case_id),
                oracle_bonds(&official["input"]["bond_fields"], case_id),
            );
            let molecule_before = molecule.clone();
            let atom_index = oracle_u32(&official["input"]["atom_index"], "atom_index", case_id);
            let rust = valence5n_cleanup2(&mut molecule, atom_index);

            match official["output"]["status"].as_str() {
                Some("return") => {
                    assert_eq!(
                        rust,
                        Ok(official["output"]["result"]
                            .as_bool()
                            .unwrap_or_else(|| panic!("{case_id}: result must be bool"))),
                        "{case_id}"
                    );
                    assert!(official["output"]["exception"].is_null(), "{case_id}");
                    assert_eq!(official["output"]["diagnostics"], Value::Array(vec![]));
                }
                Some("exception") => {
                    let exception = &official["output"]["exception"];
                    let expected = AdapterValenceError {
                        kind: "ValueErrorException",
                        message: "Bad bond type",
                        bond_index: oracle_u32(&exception["bond_index"], "bond_index", case_id),
                        bond_type: oracle_bond_type(&exception["bond_type"], case_id),
                    };
                    assert_eq!(exception["kind"], expected.kind, "{case_id}");
                    assert_eq!(exception["message"], expected.message, "{case_id}");
                    assert_eq!(rust, Err(expected.clone()), "{case_id}");
                    assert_eq!(official["output"]["result"], Value::Null, "{case_id}");
                    let diagnostics = official["output"]["diagnostics"]
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: diagnostics must be array"));
                    assert_eq!(diagnostics.len(), 1, "{case_id}");
                    assert_eq!(diagnostics[0]["kind"], expected.kind, "{case_id}");
                    assert_eq!(diagnostics[0]["message"], expected.message, "{case_id}");
                    assert_eq!(
                        oracle_u32(&diagnostics[0]["bond_index"], "bond_index", case_id),
                        expected.bond_index,
                        "{case_id}"
                    );
                    assert_eq!(
                        oracle_bond_type(&diagnostics[0]["bond_type"], case_id),
                        expected.bond_type,
                        "{case_id}"
                    );
                }
                status => panic!("{case_id}: unexpected status {status:?}"),
            }

            assert_eq!(
                official["output"]["graph_unchanged"].as_bool(),
                Some(molecule == molecule_before),
                "{case_id}"
            );
            assert_eq!(
                oracle_u32(&official["output"]["atom_count"], "atom_count", case_id),
                molecule.atoms.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                oracle_u32(&official["output"]["bond_count"], "bond_count", case_id),
                molecule.bonds.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                oracle_atoms(&official["output"]["atom_fields"], case_id),
                molecule.atoms,
                "{case_id}"
            );
            assert_eq!(
                oracle_bonds(&official["output"]["bond_fields"], case_id),
                molecule.bonds,
                "{case_id}"
            );
            assert_eq!(
                official["output"]["stereo_fields"]["bond_directions"],
                Value::Array(
                    molecule
                        .bonds
                        .iter()
                        .map(|bond| Value::from(bond.direction as u8))
                        .collect()
                ),
                "{case_id}"
            );
            assert_eq!(official["output"]["properties"], Value::Array(vec![]));
        }

        assert_eq!(
            case_ids,
            BTreeSet::from([
                "begin-is-root-preserves-fields".to_owned(),
                "multiple-targets-first-adjacency-reversed-root".to_owned(),
                "no-target".to_owned(),
                "over-depth-limit".to_owned(),
                "root-valence-exception-after-target-check".to_owned(),
                "target-valence-exception-partial-mutation".to_owned(),
            ])
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__getnumdoublebondednegativelychargedneighboringsi__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        fn i32_field(value: &Value, field: &str, case_id: &str) -> i32 {
            i32::try_from(
                value
                    .as_i64()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be i32")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds i32"))
        }

        fn u32_field(value: &Value, field: &str, case_id: &str) -> u32 {
            u32::try_from(
                value
                    .as_u64()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be u32")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds u32"))
        }

        fn bond_type(value: &Value, case_id: &str) -> BondType {
            match u32_field(value, "bond_type", case_id) {
                0 => BondType::Unspecified,
                1 => BondType::Single,
                2 => BondType::Double,
                3 => BondType::Triple,
                4 => BondType::Quadruple,
                5 => BondType::Quintuple,
                6 => BondType::Hextuple,
                7 => BondType::OneAndAHalf,
                8 => BondType::TwoAndAHalf,
                9 => BondType::ThreeAndAHalf,
                10 => BondType::FourAndAHalf,
                11 => BondType::FiveAndAHalf,
                12 => BondType::Aromatic,
                13 => BondType::Ionic,
                14 => BondType::Hydrogen,
                15 => BondType::ThreeCenter,
                16 => BondType::DativeOne,
                17 => BondType::Dative,
                18 => BondType::DativeL,
                19 => BondType::DativeR,
                20 => BondType::Other,
                21 => BondType::Zero,
                value => panic!("{case_id}: unknown BondType {value}"),
            }
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--neighboring-si-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );

        let output = String::from_utf8(oracle.stdout).expect("RDKit C++ oracle output must be UTF-8");
        let mut case_ids = BTreeSet::new();
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "a5d70ac89adb467c691ba0a8a16878e951833d598c7f7d15a1a9d241f0c17374"
            );
            assert_eq!(
                official["operation"],
                "getNumDoubleBondedNegativelyChargedNeighboringSi"
            );
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let atom_fields = official["input"]["atom_fields"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: atom_fields must be an array"));
            let atoms = atom_fields
                .iter()
                .enumerate()
                .map(|(index, field)| {
                    assert_eq!(
                        u32_field(&field["index"], "atom index", case_id),
                        index as u32,
                        "{case_id}"
                    );
                    AdapterAtom {
                        atomic_number: i32_field(&field["atomic_number"], "atomic_number", case_id),
                        formal_charge: i32_field(&field["formal_charge"], "formal_charge", case_id),
                        num_explicit_hydrogens: 0,
                        is_aromatic: false,
                        ..AdapterAtom::default()
                    }
                })
                .collect();
            let bond_fields = official["input"]["bond_fields"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: bond_fields must be an array"));
            let bonds = bond_fields
                .iter()
                .enumerate()
                .map(|(index, field)| {
                    assert_eq!(
                        u32_field(&field["index"], "bond index", case_id),
                        index as u32,
                        "{case_id}"
                    );
                    assert_eq!(field["direction"], 0, "{case_id}");
                    AdapterBond {
                        begin_atom_index: u32_field(&field["begin_atom_index"], "begin_atom_index", case_id),
                        end_atom_index: u32_field(&field["end_atom_index"], "end_atom_index", case_id),
                        bond_type: bond_type(&field["bond_type"], case_id),
                        direction: BondDirection::None,
                        ..AdapterBond::default()
                    }
                })
                .collect();
            let molecule = AdapterMol::from_graph(atoms, bonds);
            let molecule_before = molecule.clone();
            let atom_index = u32_field(&official["input"]["atom_index"], "atom_index", case_id);
            let rust = get_num_double_bonded_negatively_charged_neighboring_si(&molecule, atom_index);

            assert_eq!(
                rust,
                i32_field(&official["output"]["count"], "count", case_id),
                "{case_id}"
            );
            assert_eq!(molecule, molecule_before, "{case_id}: Rust graph mutation");
            assert_eq!(official["output"]["status"], "return", "{case_id}");
            assert!(official["output"]["exception"].is_null(), "{case_id}");
            assert_eq!(official["output"]["graph_unchanged"], true, "{case_id}");
            assert_eq!(
                official["output"]["atom_fields"], official["input"]["atom_fields"],
                "{case_id}"
            );
            assert_eq!(
                official["output"]["bond_fields"], official["input"]["bond_fields"],
                "{case_id}"
            );
            assert_eq!(
                official["output"]["atom_count"].as_u64(),
                Some(molecule.atoms.len() as u64),
                "{case_id}"
            );
            assert_eq!(
                official["output"]["bond_count"].as_u64(),
                Some(molecule.bonds.len() as u64),
                "{case_id}"
            );
            assert_eq!(
                official["output"]["stereo_fields"]["bond_directions"],
                Value::Array(vec![Value::from(0); molecule.bonds.len()]),
                "{case_id}"
            );
            assert_eq!(official["output"]["properties"], Value::Array(vec![]));
            assert_eq!(official["output"]["diagnostics"], Value::Array(vec![]));
        }
        assert_eq!(
            case_ids,
            BTreeSet::from([
                "all-miss".to_owned(),
                "each-condition-and-two-matches".to_owned(),
                "four-matches".to_owned(),
                "isolated".to_owned(),
                "single-match-reversed-endpoints".to_owned(),
            ])
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence4ncleanup1__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        fn i32_field(value: &Value, field: &str, case_id: &str) -> i32 {
            i32::try_from(
                value
                    .as_i64()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be i32")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds i32"))
        }

        fn u32_field(value: &Value, field: &str, case_id: &str) -> u32 {
            u32::try_from(
                value
                    .as_u64()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be u32")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds u32"))
        }

        fn bond_type(value: &Value, case_id: &str) -> BondType {
            match u32_field(value, "bond_type", case_id) {
                0 => BondType::Unspecified,
                1 => BondType::Single,
                2 => BondType::Double,
                3 => BondType::Triple,
                4 => BondType::Quadruple,
                5 => BondType::Quintuple,
                6 => BondType::Hextuple,
                7 => BondType::OneAndAHalf,
                8 => BondType::TwoAndAHalf,
                9 => BondType::ThreeAndAHalf,
                10 => BondType::FourAndAHalf,
                11 => BondType::FiveAndAHalf,
                12 => BondType::Aromatic,
                13 => BondType::Ionic,
                14 => BondType::Hydrogen,
                15 => BondType::ThreeCenter,
                16 => BondType::DativeOne,
                17 => BondType::Dative,
                18 => BondType::DativeL,
                19 => BondType::DativeR,
                20 => BondType::Other,
                21 => BondType::Zero,
                value => panic!("{case_id}: unknown BondType {value}"),
            }
        }

        fn direction(value: &Value, case_id: &str) -> BondDirection {
            match u32_field(value, "direction", case_id) {
                0 => BondDirection::None,
                1 => BondDirection::BeginWedge,
                2 => BondDirection::BeginDash,
                3 => BondDirection::EndDownRight,
                4 => BondDirection::EndUpRight,
                5 => BondDirection::EitherDouble,
                6 => BondDirection::Unknown,
                value => panic!("{case_id}: unknown BondDirection {value}"),
            }
        }

        fn parse_atoms(value: &Value, case_id: &str) -> Vec<AdapterAtom> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: atom_fields must be an array"))
                .iter()
                .enumerate()
                .map(|(index, field)| {
                    assert_eq!(
                        u32_field(&field["index"], "atom index", case_id),
                        index as u32,
                        "{case_id}"
                    );
                    AdapterAtom {
                        atomic_number: i32_field(&field["atomic_number"], "atomic_number", case_id),
                        formal_charge: i32_field(&field["formal_charge"], "formal_charge", case_id),
                        num_explicit_hydrogens: u32_field(
                            &field["num_explicit_hydrogens"],
                            "num_explicit_hydrogens",
                            case_id,
                        ),
                        is_aromatic: false,
                        ..AdapterAtom::default()
                    }
                })
                .collect()
        }

        fn parse_bonds(value: &Value, case_id: &str) -> Vec<AdapterBond> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: bond_fields must be an array"))
                .iter()
                .enumerate()
                .map(|(index, field)| {
                    assert_eq!(
                        u32_field(&field["index"], "bond index", case_id),
                        index as u32,
                        "{case_id}"
                    );
                    AdapterBond {
                        begin_atom_index: u32_field(&field["begin_atom_index"], "begin_atom_index", case_id),
                        end_atom_index: u32_field(&field["end_atom_index"], "end_atom_index", case_id),
                        bond_type: bond_type(&field["bond_type"], case_id),
                        direction: direction(&field["direction"], case_id),
                        ..AdapterBond::default()
                    }
                })
                .collect()
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--valence4n-cleanup1-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let mut case_ids = BTreeSet::new();

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "11132300ed0a447a05c38e8eed0e710762514b9180aad984f12acc697255d87d"
            );
            assert_eq!(official["operation"], "_Valence4NCleanUp1");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let atoms = parse_atoms(&official["input"]["atom_fields"], case_id);
            let bonds = parse_bonds(&official["input"]["bond_fields"], case_id);
            let mut molecule = AdapterMol::from_graph(atoms, bonds);
            let molecule_before = molecule.clone();
            let atom_index = u32_field(&official["input"]["atom_index"], "atom_index", case_id);
            let rust = valence4n_cleanup1(&mut molecule, atom_index);

            match official["output"]["status"].as_str() {
                Some("return") => {
                    assert_eq!(
                        rust,
                        Ok(official["output"]["result"]
                            .as_bool()
                            .unwrap_or_else(|| panic!("{case_id}: result must be bool"))),
                        "{case_id}"
                    );
                    assert!(official["output"]["exception"].is_null(), "{case_id}");
                    assert_eq!(official["output"]["diagnostics"], Value::Array(vec![]));
                }
                Some("exception") => {
                    let exception = &official["output"]["exception"];
                    assert_eq!(exception["kind"], "ValueErrorException", "{case_id}");
                    assert_eq!(exception["message"], "Bad bond type", "{case_id}");
                    let expected = AdapterValenceError {
                        kind: "ValueErrorException",
                        message: "Bad bond type",
                        bond_index: u32_field(&exception["bond_index"], "bond_index", case_id),
                        bond_type: bond_type(&exception["bond_type"], case_id),
                    };
                    assert_eq!(rust, Err(expected.clone()), "{case_id}");
                    let diagnostics = official["output"]["diagnostics"]
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: diagnostics must be array"));
                    assert_eq!(diagnostics.len(), 1, "{case_id}");
                    assert_eq!(diagnostics[0]["kind"], expected.kind, "{case_id}");
                    assert_eq!(diagnostics[0]["message"], expected.message, "{case_id}");
                    assert_eq!(
                        u32_field(&diagnostics[0]["bond_index"], "bond_index", case_id),
                        expected.bond_index,
                        "{case_id}"
                    );
                    assert_eq!(
                        bond_type(&diagnostics[0]["bond_type"], case_id),
                        expected.bond_type,
                        "{case_id}"
                    );
                    assert!(official["output"]["result"].is_null(), "{case_id}");
                }
                status => panic!("{case_id}: unexpected status {status:?}"),
            }

            assert_eq!(
                official["output"]["graph_unchanged"].as_bool(),
                Some(molecule == molecule_before),
                "{case_id}"
            );
            assert_eq!(
                u32_field(&official["output"]["atom_count"], "atom_count", case_id),
                molecule.atoms.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                u32_field(&official["output"]["bond_count"], "bond_count", case_id),
                molecule.bonds.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                parse_atoms(&official["output"]["atom_fields"], case_id),
                molecule.atoms,
                "{case_id}"
            );
            assert_eq!(
                parse_bonds(&official["output"]["bond_fields"], case_id),
                molecule.bonds,
                "{case_id}"
            );
            assert_eq!(
                official["output"]["stereo_fields"]["bond_directions"],
                Value::Array(
                    molecule
                        .bonds
                        .iter()
                        .map(|bond| Value::from(bond.direction as u8))
                        .collect()
                ),
                "{case_id}"
            );
            assert_eq!(official["output"]["properties"], Value::Array(vec![]));
        }

        assert_eq!(
            case_ids,
            BTreeSet::from([
                "explicit-h-valence-four-no-match".to_owned(),
                "multiple-unique-atom-set-matches".to_owned(),
                "unique-match-does-not-contain-argument-atom".to_owned(),
                "unique-match-reversed-endpoints-extra-edge-and-charges".to_owned(),
                "unsupported-valence-bond-exception".to_owned(),
                "wrong-charge-short-circuit".to_owned(),
                "wrong-element-short-circuit".to_owned(),
                "wrong-valence".to_owned(),
                "zero-substructure-matches".to_owned(),
            ])
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__valence4ncleanup2__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        fn i32_field(value: &Value, field: &str, case_id: &str) -> i32 {
            i32::try_from(
                value
                    .as_i64()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be i32")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds i32"))
        }

        fn u32_field(value: &Value, field: &str, case_id: &str) -> u32 {
            u32::try_from(
                value
                    .as_u64()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be u32")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds u32"))
        }

        fn bond_type(value: &Value, case_id: &str) -> BondType {
            match u32_field(value, "bond_type", case_id) {
                0 => BondType::Unspecified,
                1 => BondType::Single,
                2 => BondType::Double,
                3 => BondType::Triple,
                4 => BondType::Quadruple,
                5 => BondType::Quintuple,
                6 => BondType::Hextuple,
                7 => BondType::OneAndAHalf,
                8 => BondType::TwoAndAHalf,
                9 => BondType::ThreeAndAHalf,
                10 => BondType::FourAndAHalf,
                11 => BondType::FiveAndAHalf,
                12 => BondType::Aromatic,
                13 => BondType::Ionic,
                14 => BondType::Hydrogen,
                15 => BondType::ThreeCenter,
                16 => BondType::DativeOne,
                17 => BondType::Dative,
                18 => BondType::DativeL,
                19 => BondType::DativeR,
                20 => BondType::Other,
                21 => BondType::Zero,
                value => panic!("{case_id}: unknown BondType {value}"),
            }
        }

        fn direction(value: &Value, case_id: &str) -> BondDirection {
            match u32_field(value, "direction", case_id) {
                0 => BondDirection::None,
                1 => BondDirection::BeginWedge,
                2 => BondDirection::BeginDash,
                3 => BondDirection::EndDownRight,
                4 => BondDirection::EndUpRight,
                5 => BondDirection::EitherDouble,
                6 => BondDirection::Unknown,
                value => panic!("{case_id}: unknown BondDirection {value}"),
            }
        }

        fn parse_atoms(value: &Value, case_id: &str) -> Vec<AdapterAtom> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: atom_fields must be an array"))
                .iter()
                .enumerate()
                .map(|(index, field)| {
                    assert_eq!(
                        u32_field(&field["index"], "atom index", case_id),
                        index as u32,
                        "{case_id}"
                    );
                    AdapterAtom {
                        atomic_number: i32_field(&field["atomic_number"], "atomic_number", case_id),
                        formal_charge: i32_field(&field["formal_charge"], "formal_charge", case_id),
                        num_explicit_hydrogens: u32_field(
                            &field["num_explicit_hydrogens"],
                            "num_explicit_hydrogens",
                            case_id,
                        ),
                        is_aromatic: false,
                        ..AdapterAtom::default()
                    }
                })
                .collect()
        }

        fn parse_bonds(value: &Value, case_id: &str) -> Vec<AdapterBond> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: bond_fields must be an array"))
                .iter()
                .enumerate()
                .map(|(index, field)| {
                    assert_eq!(
                        u32_field(&field["index"], "bond index", case_id),
                        index as u32,
                        "{case_id}"
                    );
                    AdapterBond {
                        begin_atom_index: u32_field(&field["begin_atom_index"], "begin_atom_index", case_id),
                        end_atom_index: u32_field(&field["end_atom_index"], "end_atom_index", case_id),
                        bond_type: bond_type(&field["bond_type"], case_id),
                        direction: direction(&field["direction"], case_id),
                        ..AdapterBond::default()
                    }
                })
                .collect()
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--valence4n-cleanup2-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output = String::from_utf8(oracle.stdout).expect("oracle output must be UTF-8 JSONL");
        let mut case_ids = BTreeSet::new();

        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "c2fe8b4125473236e4fa1c638d93c23365c1c7df416f17537259dcc24c4d7587"
            );
            assert_eq!(official["operation"], "_Valence4NCleanUp2");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let atoms = parse_atoms(&official["input"]["atom_fields"], case_id);
            let bonds = parse_bonds(&official["input"]["bond_fields"], case_id);
            let mut molecule = AdapterMol::from_graph(atoms, bonds);
            let molecule_before = molecule.clone();
            let atom_index = u32_field(&official["input"]["atom_index"], "atom_index", case_id);
            let rust = valence4n_cleanup2(&mut molecule, atom_index);

            assert_eq!(official["output"]["status"], "return", "{case_id}");
            assert_eq!(official["output"]["result"].as_bool(), Some(rust), "{case_id}");
            assert!(official["output"]["exception"].is_null(), "{case_id}");
            assert_eq!(official["output"]["diagnostics"], Value::Array(vec![]), "{case_id}");
            assert_eq!(
                official["output"]["graph_unchanged"].as_bool(),
                Some(molecule == molecule_before),
                "{case_id}"
            );
            assert_eq!(
                u32_field(&official["output"]["atom_count"], "atom_count", case_id),
                molecule.atoms.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                u32_field(&official["output"]["bond_count"], "bond_count", case_id),
                molecule.bonds.len() as u32,
                "{case_id}"
            );
            assert_eq!(
                parse_atoms(&official["output"]["atom_fields"], case_id),
                molecule.atoms,
                "{case_id}"
            );
            assert_eq!(
                parse_bonds(&official["output"]["bond_fields"], case_id),
                molecule.bonds,
                "{case_id}"
            );
            assert_eq!(
                official["output"]["stereo_fields"]["bond_directions"],
                Value::Array(
                    molecule
                        .bonds
                        .iter()
                        .map(|bond| Value::from(bond.direction as u8))
                        .collect()
                ),
                "{case_id}"
            );
            assert_eq!(official["output"]["properties"], Value::Array(vec![]), "{case_id}");
        }

        assert_eq!(
            case_ids,
            BTreeSet::from([
                "all-target-predicates-miss".to_owned(),
                "isolated-no-target".to_owned(),
                "miss-before-later-target".to_owned(),
                "multiple-targets-first-insertion".to_owned(),
                "multiple-targets-reversed-insertion".to_owned(),
                "unique-target-reversed-endpoints-preserves-fields".to_owned(),
            ])
        );
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__assignbonddirs__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        fn parse_direction(value: &Value, case_id: &str) -> BondDirection {
            match value
                .as_u64()
                .unwrap_or_else(|| panic!("{case_id}: direction must be unsigned"))
            {
                0 => BondDirection::None,
                1 => BondDirection::BeginWedge,
                2 => BondDirection::BeginDash,
                3 => BondDirection::EndDownRight,
                4 => BondDirection::EndUpRight,
                5 => BondDirection::EitherDouble,
                6 => BondDirection::Unknown,
                direction => panic!("{case_id}: unknown BondDir {direction}"),
            }
        }

        fn parse_directions(value: &Value, case_id: &str) -> Vec<BondDirection> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: directions must be an array"))
                .iter()
                .map(|value| parse_direction(value, case_id))
                .collect()
        }

        fn parse_pairs(value: &Value, case_id: &str) -> Vec<(i32, i32)> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: pairs must be an array"))
                .iter()
                .map(|pair| {
                    let pair = pair
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: pair must be an array"));
                    assert_eq!(pair.len(), 2, "{case_id}: pair length");
                    let first = i32::try_from(
                        pair[0]
                            .as_i64()
                            .unwrap_or_else(|| panic!("{case_id}: first must be i32")),
                    )
                    .unwrap_or_else(|_| panic!("{case_id}: first exceeds i32"));
                    let second = i32::try_from(
                        pair[1]
                            .as_i64()
                            .unwrap_or_else(|| panic!("{case_id}: second must be i32")),
                    )
                    .unwrap_or_else(|_| panic!("{case_id}: second exceeds i32"));
                    (first, second)
                })
                .collect()
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--assign-bond-dirs-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );

        let output = String::from_utf8(oracle.stdout).expect("RDKit C++ oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(official["operation"], "assignBondDirs");
            let case_id = official["case_id"].as_str().expect("case_id must be text");

            let initial = parse_directions(&official["input"]["initial_directions"], case_id);
            let z_pairs = parse_pairs(&official["input"]["z_pairs"], case_id);
            let e_pairs = parse_pairs(&official["input"]["e_pairs"], case_id);
            let z_before = z_pairs.clone();
            let e_before = e_pairs.clone();
            let mut molecule = molecule(&initial);
            let rust = assign_bond_dirs(&mut molecule, &z_pairs, &e_pairs);

            assert_eq!(official["output"]["z_pairs_unchanged"], true, "{case_id}");
            assert_eq!(official["output"]["e_pairs_unchanged"], true, "{case_id}");
            assert_eq!(z_pairs, z_before, "{case_id}");
            assert_eq!(e_pairs, e_before, "{case_id}");
            assert_eq!(
                official["output"]["bond_count"].as_u64(),
                Some(molecule.bonds.len() as u64),
                "{case_id}"
            );

            let expected_directions =
                parse_directions(&official["output"]["stereo_fields"]["bond_directions"], case_id);
            assert_eq!(directions(&molecule), expected_directions, "{case_id}");
            let bond_fields = official["output"]["bond_fields"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: bond_fields must be an array"));
            assert_eq!(bond_fields.len(), molecule.bonds.len(), "{case_id}");
            for (index, (field, bond)) in bond_fields.iter().zip(&molecule.bonds).enumerate() {
                assert_eq!(field["index"].as_u64(), Some(index as u64), "{case_id}");
                assert_eq!(
                    parse_direction(&field["direction"], case_id),
                    bond.direction,
                    "{case_id}"
                );
            }
            assert_eq!(official["output"]["atom_fields"], Value::Array(vec![]));
            assert_eq!(official["output"]["properties"], Value::Array(vec![]));

            match official["output"]["status"].as_str().expect("status must be text") {
                "return" => {
                    let expected = official["output"]["result"]
                        .as_bool()
                        .unwrap_or_else(|| panic!("{case_id}: result must be bool"));
                    assert_eq!(rust, Ok(expected), "{case_id}");
                    assert!(official["output"]["exception"].is_null(), "{case_id}");
                    assert_eq!(official["output"]["diagnostics"], Value::Array(vec![]), "{case_id}");
                }
                "exception" => {
                    assert!(official["output"]["result"].is_null(), "{case_id}");
                    let exception = &official["output"]["exception"];
                    let index = u32::try_from(
                        exception["index"]
                            .as_u64()
                            .unwrap_or_else(|| panic!("{case_id}: index must be u32")),
                    )
                    .unwrap_or_else(|_| panic!("{case_id}: index exceeds u32"));
                    let upper_bound = u32::try_from(
                        exception["upper_bound"]
                            .as_u64()
                            .unwrap_or_else(|| panic!("{case_id}: upper_bound must be u32")),
                    )
                    .unwrap_or_else(|_| panic!("{case_id}: upper_bound exceeds u32"));
                    assert_eq!(exception["kind"], "Range Error", "{case_id}");
                    assert_eq!(exception["expression"], "idx", "{case_id}");
                    assert_eq!(exception["detail"], format!("{index} < {upper_bound}"), "{case_id}");
                    assert_eq!(
                        rust,
                        Err(BondIndexRangeError {
                            kind: "Range Error",
                            expression: "idx",
                            index,
                            upper_bound,
                        }),
                        "{case_id}"
                    );
                    let diagnostics = official["output"]["diagnostics"]
                        .as_array()
                        .unwrap_or_else(|| panic!("{case_id}: diagnostics must be an array"));
                    assert_eq!(diagnostics.len(), 1, "{case_id}");
                    assert_eq!(diagnostics[0]["kind"], "Range Error", "{case_id}");
                    assert_eq!(diagnostics[0]["expression"], "idx", "{case_id}");
                    assert_eq!(
                        diagnostics[0]["detail"],
                        format!("{index} < {upper_bound}"),
                        "{case_id}"
                    );
                }
                status => panic!("{case_id}: unexpected status {status}"),
            }
            record_count += 1;
        }
        assert_eq!(record_count, 13);
    }

    #[test]
    #[ignore = "requires the pinned vendored RDKit source; run explicitly with --ignored"]
    fn rdkit_cpp_oracle__findalternatingbonds__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        fn parse_i32(value: &Value, field: &str, case_id: &str) -> i32 {
            i32::try_from(
                value
                    .as_i64()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be an integer")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds i32"))
        }

        fn parse_u32(value: &Value, field: &str, case_id: &str) -> u32 {
            u32::try_from(
                value
                    .as_u64()
                    .unwrap_or_else(|| panic!("{case_id}: {field} must be unsigned")),
            )
            .unwrap_or_else(|_| panic!("{case_id}: {field} exceeds u32"))
        }

        fn parse_bond_type(value: &Value, case_id: &str) -> BondType {
            match parse_u32(value, "bond_type", case_id) {
                0 => BondType::Unspecified,
                1 => BondType::Single,
                2 => BondType::Double,
                3 => BondType::Triple,
                4 => BondType::Quadruple,
                5 => BondType::Quintuple,
                6 => BondType::Hextuple,
                7 => BondType::OneAndAHalf,
                8 => BondType::TwoAndAHalf,
                9 => BondType::ThreeAndAHalf,
                10 => BondType::FourAndAHalf,
                11 => BondType::FiveAndAHalf,
                12 => BondType::Aromatic,
                13 => BondType::Ionic,
                14 => BondType::Hydrogen,
                15 => BondType::ThreeCenter,
                16 => BondType::DativeOne,
                17 => BondType::Dative,
                18 => BondType::DativeL,
                19 => BondType::DativeR,
                20 => BondType::Other,
                21 => BondType::Zero,
                bond_type => panic!("{case_id}: unknown BondType {bond_type}"),
            }
        }

        fn parse_direction(value: &Value, case_id: &str) -> BondDirection {
            match parse_u32(value, "direction", case_id) {
                0 => BondDirection::None,
                1 => BondDirection::BeginWedge,
                2 => BondDirection::BeginDash,
                3 => BondDirection::EndDownRight,
                4 => BondDirection::EndUpRight,
                5 => BondDirection::EitherDouble,
                6 => BondDirection::Unknown,
                direction => panic!("{case_id}: unknown BondDir {direction}"),
            }
        }

        fn parse_u32_array(value: &Value, field: &str, case_id: &str) -> Vec<u32> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: {field} must be an array"))
                .iter()
                .map(|value| parse_u32(value, field, case_id))
                .collect()
        }

        fn parse_visited(value: &Value, case_id: &str) -> BTreeSet<i32> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: visited must be an array"))
                .iter()
                .map(|value| parse_i32(value, "visited", case_id))
                .collect()
        }

        fn parse_atoms(value: &Value, case_id: &str) -> Vec<AdapterAtom> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: atom_fields must be an array"))
                .iter()
                .enumerate()
                .map(|(index, field)| {
                    assert_eq!(
                        parse_u32(&field["index"], "atom index", case_id),
                        index as u32,
                        "{case_id}: atom index"
                    );
                    AdapterAtom {
                        atomic_number: parse_i32(&field["atomic_number"], "atomic_number", case_id),
                        formal_charge: parse_i32(&field["formal_charge"], "formal_charge", case_id),
                        num_explicit_hydrogens: 0,
                        is_aromatic: false,
                        ..AdapterAtom::default()
                    }
                })
                .collect()
        }

        fn parse_bonds(value: &Value, case_id: &str) -> Vec<AdapterBond> {
            value
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: bond_fields must be an array"))
                .iter()
                .enumerate()
                .map(|(index, field)| {
                    assert_eq!(
                        parse_u32(&field["index"], "bond index", case_id),
                        index as u32,
                        "{case_id}: bond index"
                    );
                    AdapterBond {
                        begin_atom_index: parse_u32(&field["begin_atom_index"], "begin_atom_index", case_id),
                        end_atom_index: parse_u32(&field["end_atom_index"], "end_atom_index", case_id),
                        bond_type: parse_bond_type(&field["bond_type"], case_id),
                        direction: parse_direction(&field["direction"], case_id),
                        ..AdapterBond::default()
                    }
                })
                .collect()
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/rdkit_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--find-alternating-bonds-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "pinned RDKit C++ oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );

        let output = String::from_utf8(oracle.stdout).expect("RDKit C++ oracle output must be UTF-8");
        let mut case_ids = BTreeSet::new();
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-rdkit-cpp-v1");
            assert_eq!(official["rdkit_version"], "2026.03.1");
            assert_eq!(
                official["source_sha256"],
                "104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f"
            );
            assert_eq!(
                official["source_fragment_sha256"],
                "ef0a05d3a018d6928d66af0d5238dea5bb22742e233f1fe8229f62578fe7d003"
            );
            assert_eq!(official["operation"], "findAlternatingBonds");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            assert!(case_ids.insert(case_id.to_owned()), "duplicate {case_id}");

            let input = &official["input"];
            let atoms = parse_atoms(&input["atom_fields"], case_id);
            let bonds = parse_bonds(&input["bond_fields"], case_id);
            let mut molecule = AdapterMol::from_graph(atoms, bonds);
            let molecule_before = molecule.clone();
            let current_atom_index = parse_u32(&input["current_atom_index"], "current_atom_index", case_id);
            let desired_atomic_number = parse_i32(&input["desired_atomic_number"], "desired_atomic_number", case_id);
            let desired_atom_charge = parse_i32(&input["desired_atom_charge"], "desired_atom_charge", case_id);
            let desired_next_bond_type = parse_bond_type(&input["desired_next_bond_type"], case_id);
            let desired_ending_bond_type = parse_bond_type(&input["desired_ending_bond_type"], case_id);
            let current_path_length = parse_u32(&input["current_path_length"], "current_path_length", case_id);
            let max_path_length = parse_u32(&input["max_path_length"], "max_path_length", case_id);
            let last_bond_index = if input["last_bond_index"].is_null() {
                None
            } else {
                Some(parse_u32(&input["last_bond_index"], "last_bond_index", case_id))
            };
            let mut path = parse_u32_array(&input["initial_path"], "initial_path", case_id);
            let mut visited = parse_visited(&input["initial_visited"], case_id);

            let rust_target = find_alternating_bonds(
                &molecule,
                current_atom_index,
                desired_atomic_number,
                desired_atom_charge,
                desired_next_bond_type,
                desired_ending_bond_type,
                current_path_length,
                max_path_length,
                last_bond_index,
                &mut path,
                &mut visited,
            );

            let expected_target = if official["output"]["target_atom_index"].is_null() {
                None
            } else {
                Some(parse_u32(
                    &official["output"]["target_atom_index"],
                    "target_atom_index",
                    case_id,
                ))
            };
            assert_eq!(rust_target, expected_target, "{case_id}: target");
            assert_eq!(
                path,
                parse_u32_array(&official["output"]["path"], "path", case_id),
                "{case_id}: path"
            );
            assert_eq!(
                visited,
                parse_visited(&official["output"]["visited"], case_id),
                "{case_id}: visited"
            );
            assert_eq!(molecule, molecule_before, "{case_id}: Rust graph mutation");
            assert_eq!(official["output"]["graph_unchanged"], true, "{case_id}");
            assert_eq!(official["output"]["status"], "return", "{case_id}");
            assert!(official["output"]["exception"].is_null(), "{case_id}");
            assert_eq!(official["output"]["diagnostics"], Value::Array(vec![]), "{case_id}");
            assert_eq!(official["output"]["properties"], Value::Array(vec![]), "{case_id}");
            assert_eq!(
                official["output"]["atom_count"].as_u64(),
                Some(molecule.atoms.len() as u64),
                "{case_id}"
            );
            assert_eq!(
                official["output"]["bond_count"].as_u64(),
                Some(molecule.bonds.len() as u64),
                "{case_id}"
            );
            assert_eq!(
                parse_atoms(&official["output"]["atom_fields"], case_id),
                molecule.atoms,
                "{case_id}: atom fields"
            );
            assert_eq!(
                parse_bonds(&official["output"]["bond_fields"], case_id),
                molecule.bonds,
                "{case_id}: bond fields"
            );
            let expected_directions = official["output"]["stereo_fields"]["bond_directions"]
                .as_array()
                .unwrap_or_else(|| panic!("{case_id}: bond_directions must be an array"))
                .iter()
                .map(|value| parse_direction(value, case_id))
                .collect::<Vec<_>>();
            assert_eq!(directions(&molecule), expected_directions, "{case_id}");

            // The source accepts a mutable ROMol reference but this function is read-only.
            // Keep the local binding mutable to mirror that API and prove no state changed.
            let _ = &mut molecule;
        }

        let expected_case_ids = BTreeSet::from([
            "atomic-number-miss".to_owned(),
            "cycle-global-visited".to_owned(),
            "direct-target-rejects-equal-path".to_owned(),
            "direct-target-replaces-longer-path".to_owned(),
            "double-start-defaults-next-to-single".to_owned(),
            "ending-bond-type-miss".to_owned(),
            "equal-path-keeps-first-insertion".to_owned(),
            "equal-path-reversed-insertion".to_owned(),
            "formal-charge-miss".to_owned(),
            "near-unsigned-path-limit".to_owned(),
            "negative-atomic-number".to_owned(),
            "nonalternating-ending-success".to_owned(),
            "nonalternating-ending-target-miss".to_owned(),
            "nonroot-preserves-visited-and-skips-neighbor".to_owned(),
            "root-clears-state-at-cutoff".to_owned(),
            "shorter-path-found-later".to_owned(),
            "single-double-ending-disables-special-branch".to_owned(),
            "single-double-single-alternation".to_owned(),
            "triple-then-single-special-case".to_owned(),
            "unmatched-edge-type".to_owned(),
        ]);
        assert_eq!(case_ids, expected_case_ids);
    }
}
