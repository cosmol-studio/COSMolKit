#![allow(
    dead_code,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_imports
)]

use std::any::Any;
use std::collections::BTreeMap;
use std::fmt;
use std::marker::PhantomData;

mod generated;

pub(crate) use generated::*;

/// Source-level replacement for C `void` in typed allocation handles.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct SourceVoid;

/// Owned byte stream replacing the `FILE *` branch of `INCHI_IOSTREAM`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct SourceFile {
    pub(crate) bytes: Vec<u8>,
    pub(crate) position: u64,
    pub(crate) error: bool,
    pub(crate) eof: bool,
    pub(crate) is_standard_stream: bool,
}

/// Explicit argument state replacing C `va_list` in the formatting port.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct SourceVaList {
    pub(crate) arguments: Vec<SourceFormatArgument>,
    pub(crate) position: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum SourceFormatArgument {
    Signed(i64),
    Unsigned(u64),
    Float(f64),
    Byte(i8),
    Bytes(SourceConstPointer<i8>),
    Pointer(SourceConstPointer<SourceVoid>),
    MutSigned(SourceMutPointer<i32>),
}

/// Pointer states produced by the C options tokenizer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SourceArgvPointer {
    Null,
    EmptyLiteral,
    Command(SourceMutPointer<i8>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct AllocationId(u64);

pub(crate) struct SourceMutPointer<T> {
    allocation: Option<AllocationId>,
    element_offset: u64,
    marker: PhantomData<fn() -> T>,
}

pub(crate) struct SourceConstPointer<T> {
    allocation: Option<AllocationId>,
    element_offset: u64,
    marker: PhantomData<fn() -> T>,
}

macro_rules! impl_source_pointer_value_traits {
    ($pointer:ident) => {
        impl<T> Copy for $pointer<T> {}

        impl<T> Clone for $pointer<T> {
            fn clone(&self) -> Self {
                *self
            }
        }

        impl<T> fmt::Debug for $pointer<T> {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter
                    .debug_struct(stringify!($pointer))
                    .field("allocation", &self.allocation)
                    .field("element_offset", &self.element_offset)
                    .finish()
            }
        }

        impl<T> PartialEq for $pointer<T> {
            fn eq(&self, other: &Self) -> bool {
                self.allocation == other.allocation && self.element_offset == other.element_offset
            }
        }

        impl<T> Eq for $pointer<T> {}

        impl<T> Default for $pointer<T> {
            fn default() -> Self {
                Self::null()
            }
        }

        impl<T> $pointer<T> {
            pub(crate) const fn null() -> Self {
                Self {
                    allocation: None,
                    element_offset: 0,
                    marker: PhantomData,
                }
            }

            pub(crate) const fn is_null(self) -> bool {
                self.allocation.is_none()
            }

            pub(crate) fn offset(self, elements: i64) -> Result<Self, SourceHeapError> {
                let magnitude = elements.unsigned_abs();
                let element_offset = if elements.is_negative() {
                    self.element_offset.checked_sub(magnitude)
                } else {
                    self.element_offset.checked_add(magnitude)
                }
                .ok_or(SourceHeapError::PointerOffsetOverflow)?;
                Ok(Self {
                    element_offset,
                    ..self
                })
            }
        }
    };
}

impl_source_pointer_value_traits!(SourceMutPointer);
impl_source_pointer_value_traits!(SourceConstPointer);

impl<T> SourceMutPointer<T> {
    pub(crate) const fn as_const(self) -> SourceConstPointer<T> {
        SourceConstPointer {
            allocation: self.allocation,
            element_offset: self.element_offset,
            marker: PhantomData,
        }
    }

    pub(crate) fn difference(self, origin: Self) -> Result<i64, SourceHeapError> {
        let allocation = self.allocation.ok_or(SourceHeapError::NullPointer)?;
        if origin.allocation != Some(allocation) {
            return Err(SourceHeapError::PointerAllocationMismatch);
        }
        let difference = i128::from(self.element_offset) - i128::from(origin.element_offset);
        i64::try_from(difference).map_err(|_| SourceHeapError::PointerDifferenceOverflow)
    }
}

impl<T> SourceConstPointer<T> {
    /// Reconstitutes the mutable allocation handle where the C source casts
    /// away a const-qualified owner before `free`.
    pub(crate) const fn as_mut(self) -> SourceMutPointer<T> {
        SourceMutPointer {
            allocation: self.allocation,
            element_offset: self.element_offset,
            marker: PhantomData,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SourceHeapError {
    AllocationIdExhausted,
    AllocationSizeOverflow,
    AllocationElementCountOutOfRange,
    AllocationFailed,
    SourceIntegerOverflow,
    MissingNulTerminator,
    InvalidSourceTextEncoding,
    PointerAllocationMismatch,
    PointerDifferenceOverflow,
    NullPointer,
    MissingAllocation,
    AllocationTypeMismatch,
    PointerOffsetOverflow,
    PointerOutOfBounds,
    FreeOfInteriorPointer,
    UnsupportedSourceBehavior,
}

/// Owns allocations referenced by source pointer values.
///
/// Handles preserve C nullability, allocation identity, aliasing, interior
/// pointers, and one-past pointers without exposing native addresses.
#[derive(Default)]
pub(crate) struct SourceHeap {
    next_id: u64,
    allocations: BTreeMap<AllocationId, Box<dyn Any>>,
    source_errno: i32,
    #[cfg(test)]
    allocations_before_failure: Option<u64>,
    #[cfg(test)]
    source_allocation_calls: u64,
}

const INP_ATOM_GCC_LP64_SIZE: usize = 176;

fn inp_atom_gcc_lp64_bytes(atom: &inp_ATOM) -> [u8; INP_ATOM_GCC_LP64_SIZE] {
    let mut bytes = [0_u8; INP_ATOM_GCC_LP64_SIZE];
    for (target, &value) in bytes[0..6].iter_mut().zip(&atom.elname) {
        *target = value as u8;
    }
    bytes[6] = atom.el_number;
    for (index, value) in atom.neighbor.iter().enumerate() {
        bytes[8 + index * 2..10 + index * 2].copy_from_slice(&value.to_ne_bytes());
    }
    bytes[48..50].copy_from_slice(&atom.orig_at_number.to_ne_bytes());
    bytes[50..52].copy_from_slice(&atom.orig_compt_at_numb.to_ne_bytes());
    for (target, &value) in bytes[52..72].iter_mut().zip(&atom.bond_stereo) {
        *target = value as u8;
    }
    bytes[72..92].copy_from_slice(&atom.bond_type);
    bytes[92] = atom.valence as u8;
    bytes[93] = atom.chem_bonds_valence as u8;
    bytes[94] = atom.num_H as u8;
    for (target, &value) in bytes[95..98].iter_mut().zip(&atom.num_iso_H) {
        *target = value as u8;
    }
    bytes[98] = atom.iso_atw_diff as u8;
    bytes[99] = atom.charge as u8;
    bytes[100] = atom.radical as u8;
    bytes[101] = atom.bAmbiguousStereo as u8;
    bytes[102] = atom.cFlags as u8;
    bytes[104..106].copy_from_slice(&atom.at_type.to_ne_bytes());
    bytes[106..108].copy_from_slice(&atom.component.to_ne_bytes());
    bytes[108..110].copy_from_slice(&atom.endpoint.to_ne_bytes());
    bytes[110..112].copy_from_slice(&atom.c_point.to_ne_bytes());
    bytes[112..120].copy_from_slice(&atom.x.to_ne_bytes());
    bytes[120..128].copy_from_slice(&atom.y.to_ne_bytes());
    bytes[128..136].copy_from_slice(&atom.z.to_ne_bytes());
    bytes[136] = atom.bUsed0DParity as u8;
    bytes[137] = atom.p_parity as u8;
    for (index, value) in atom.p_orig_at_num.iter().enumerate() {
        bytes[138 + index * 2..140 + index * 2].copy_from_slice(&value.to_ne_bytes());
    }
    for (target, &value) in bytes[146..149].iter_mut().zip(&atom.sb_ord) {
        *target = value as u8;
    }
    for (target, &value) in bytes[149..152].iter_mut().zip(&atom.sn_ord) {
        *target = value as u8;
    }
    for (target, &value) in bytes[152..155].iter_mut().zip(&atom.sb_parity) {
        *target = value as u8;
    }
    for (index, value) in atom.sn_orig_at_num.iter().enumerate() {
        bytes[156 + index * 2..158 + index * 2].copy_from_slice(&value.to_ne_bytes());
    }
    bytes[162] = atom.bCutVertex as u8;
    bytes[164..166].copy_from_slice(&atom.nRingSystem.to_ne_bytes());
    bytes[166..168].copy_from_slice(&atom.nNumAtInRingSystem.to_ne_bytes());
    bytes[168..170].copy_from_slice(&atom.nBlockSystem.to_ne_bytes());
    bytes
}

fn overwrite_inp_atom_from_gcc_lp64_bytes(
    atom: &mut inp_ATOM,
    bytes: &[u8; INP_ATOM_GCC_LP64_SIZE],
) {
    for (target, &value) in atom.elname.iter_mut().zip(&bytes[0..6]) {
        *target = value as i8;
    }
    atom.el_number = bytes[6];
    for (index, target) in atom.neighbor.iter_mut().enumerate() {
        *target = u16::from_ne_bytes([bytes[8 + index * 2], bytes[9 + index * 2]]);
    }
    atom.orig_at_number = u16::from_ne_bytes([bytes[48], bytes[49]]);
    atom.orig_compt_at_numb = u16::from_ne_bytes([bytes[50], bytes[51]]);
    for (target, &value) in atom.bond_stereo.iter_mut().zip(&bytes[52..72]) {
        *target = value as i8;
    }
    atom.bond_type.copy_from_slice(&bytes[72..92]);
    atom.valence = bytes[92] as i8;
    atom.chem_bonds_valence = bytes[93] as i8;
    atom.num_H = bytes[94] as i8;
    for (target, &value) in atom.num_iso_H.iter_mut().zip(&bytes[95..98]) {
        *target = value as i8;
    }
    atom.iso_atw_diff = bytes[98] as i8;
    atom.charge = bytes[99] as i8;
    atom.radical = bytes[100] as i8;
    atom.bAmbiguousStereo = bytes[101] as i8;
    atom.cFlags = bytes[102] as i8;
    atom.at_type = u16::from_ne_bytes([bytes[104], bytes[105]]);
    atom.component = u16::from_ne_bytes([bytes[106], bytes[107]]);
    atom.endpoint = u16::from_ne_bytes([bytes[108], bytes[109]]);
    atom.c_point = u16::from_ne_bytes([bytes[110], bytes[111]]);
    atom.x = f64::from_ne_bytes(bytes[112..120].try_into().expect("fixed field width"));
    atom.y = f64::from_ne_bytes(bytes[120..128].try_into().expect("fixed field width"));
    atom.z = f64::from_ne_bytes(bytes[128..136].try_into().expect("fixed field width"));
    atom.bUsed0DParity = bytes[136] as i8;
    atom.p_parity = bytes[137] as i8;
    for (index, target) in atom.p_orig_at_num.iter_mut().enumerate() {
        *target = u16::from_ne_bytes([bytes[138 + index * 2], bytes[139 + index * 2]]);
    }
    for (target, &value) in atom.sb_ord.iter_mut().zip(&bytes[146..149]) {
        *target = value as i8;
    }
    for (target, &value) in atom.sn_ord.iter_mut().zip(&bytes[149..152]) {
        *target = value as i8;
    }
    for (target, &value) in atom.sb_parity.iter_mut().zip(&bytes[152..155]) {
        *target = value as i8;
    }
    for (index, target) in atom.sn_orig_at_num.iter_mut().enumerate() {
        *target = u16::from_ne_bytes([bytes[156 + index * 2], bytes[157 + index * 2]]);
    }
    atom.bCutVertex = bytes[162] as i8;
    atom.nRingSystem = u16::from_ne_bytes([bytes[164], bytes[165]]);
    atom.nNumAtInRingSystem = u16::from_ne_bytes([bytes[166], bytes[167]]);
    atom.nBlockSystem = u16::from_ne_bytes([bytes[168], bytes[169]]);
}

pub(crate) fn copy_inp_atom_gcc_lp64_byte_prefix(
    destination: &mut [inp_ATOM],
    source: &[inp_ATOM],
    byte_count: usize,
) -> Result<(), SourceHeapError> {
    let touched_atoms = byte_count.div_ceil(INP_ATOM_GCC_LP64_SIZE);
    if touched_atoms > destination.len() || touched_atoms > source.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let mut remaining = byte_count;
    for index in 0..touched_atoms {
        let source_bytes = inp_atom_gcc_lp64_bytes(&source[index]);
        let mut destination_bytes = inp_atom_gcc_lp64_bytes(&destination[index]);
        let copied = remaining.min(INP_ATOM_GCC_LP64_SIZE);
        destination_bytes[..copied].copy_from_slice(&source_bytes[..copied]);
        overwrite_inp_atom_from_gcc_lp64_bytes(&mut destination[index], &destination_bytes);
        remaining -= copied;
    }
    Ok(())
}

impl SourceHeap {
    pub(crate) const fn source_errno(&self) -> i32 {
        self.source_errno
    }

    pub(crate) fn set_source_errno(&mut self, value: i32) {
        self.source_errno = value;
    }

    #[cfg(test)]
    pub(crate) fn fail_after_allocations(&mut self, successful_allocations: u64) {
        self.allocations_before_failure = Some(successful_allocations);
        self.source_allocation_calls = 0;
    }

    #[cfg(test)]
    pub(crate) fn trace_source_allocations(&mut self) {
        self.allocations_before_failure = None;
        self.source_allocation_calls = 0;
    }

    #[cfg(test)]
    pub(crate) fn source_allocation_calls(&self) -> u64 {
        self.source_allocation_calls
    }

    #[cfg(test)]
    pub(crate) fn live_allocation_count(&self) -> usize {
        self.allocations.len()
    }

    #[cfg(test)]
    pub(crate) fn live_allocations_of<T: 'static>(&self) -> usize {
        self.allocations
            .values()
            .filter(|allocation| allocation.is::<Vec<T>>())
            .count()
    }

    pub(crate) fn allocate<T: 'static>(
        &mut self,
        values: Vec<T>,
    ) -> Result<SourceMutPointer<T>, SourceHeapError> {
        #[cfg(test)]
        {
            self.source_allocation_calls += 1;
        }
        #[cfg(test)]
        if self.allocations_before_failure == Some(0) {
            self.allocations_before_failure = None;
            return Err(SourceHeapError::AllocationFailed);
        }
        #[cfg(test)]
        if let Some(remaining) = &mut self.allocations_before_failure {
            *remaining -= 1;
        }
        self.allocate_model_storage(values)
    }

    pub(crate) fn allocate_model_storage<T: 'static>(
        &mut self,
        values: Vec<T>,
    ) -> Result<SourceMutPointer<T>, SourceHeapError> {
        let id = AllocationId(self.next_id);
        self.next_id = self
            .next_id
            .checked_add(1)
            .ok_or(SourceHeapError::AllocationIdExhausted)?;
        self.allocations.insert(id, Box::new(values));
        Ok(SourceMutPointer {
            allocation: Some(id),
            element_offset: 0,
            marker: PhantomData,
        })
    }

    pub(crate) fn slice<T: 'static>(
        &self,
        pointer: SourceConstPointer<T>,
    ) -> Result<&[T], SourceHeapError> {
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let values = self
            .allocations
            .get(&id)
            .ok_or(SourceHeapError::MissingAllocation)?
            .downcast_ref::<Vec<T>>()
            .ok_or(SourceHeapError::AllocationTypeMismatch)?;
        let offset = usize::try_from(pointer.element_offset)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        values
            .get(offset..)
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }

    pub(crate) fn slice_mut<T: 'static>(
        &mut self,
        pointer: SourceMutPointer<T>,
    ) -> Result<&mut [T], SourceHeapError> {
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let values = self
            .allocations
            .get_mut(&id)
            .ok_or(SourceHeapError::MissingAllocation)?
            .downcast_mut::<Vec<T>>()
            .ok_or(SourceHeapError::AllocationTypeMismatch)?;
        let offset = usize::try_from(pointer.element_offset)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        values
            .get_mut(offset..)
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }

    pub(crate) fn with_slice_mut_and_heap<T: 'static, R>(
        &mut self,
        pointer: SourceMutPointer<T>,
        operation: impl FnOnce(&mut [T], &SourceHeap) -> Result<R, SourceHeapError>,
    ) -> Result<R, SourceHeapError> {
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let mut allocation = self
            .allocations
            .remove(&id)
            .ok_or(SourceHeapError::MissingAllocation)?;
        let offset = match usize::try_from(pointer.element_offset) {
            Ok(offset) => offset,
            Err(_) => {
                self.allocations.insert(id, allocation);
                return Err(SourceHeapError::PointerOutOfBounds);
            }
        };
        let result = match allocation.downcast_mut::<Vec<T>>() {
            Some(values) => match values.get_mut(offset..) {
                Some(values) => operation(values, self),
                None => Err(SourceHeapError::PointerOutOfBounds),
            },
            None => Err(SourceHeapError::AllocationTypeMismatch),
        };
        self.allocations.insert(id, allocation);
        result
    }

    pub(crate) fn free<T: 'static>(
        &mut self,
        pointer: SourceMutPointer<T>,
    ) -> Result<(), SourceHeapError> {
        let Some(id) = pointer.allocation else {
            return Ok(());
        };
        if pointer.element_offset != 0 {
            return Err(SourceHeapError::FreeOfInteriorPointer);
        }
        let allocation = self
            .allocations
            .get(&id)
            .ok_or(SourceHeapError::MissingAllocation)?;
        if !allocation.is::<Vec<T>>() {
            return Err(SourceHeapError::AllocationTypeMismatch);
        }
        self.allocations.remove(&id);
        Ok(())
    }
}

/// Owned bit representation of `union tagSplitLong` for the pinned LP64,
/// little-endian source profile.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct tagSplitLong {
    pub(crate) bits: u64,
}

impl tagSplitLong {
    pub(crate) const fn from_ulong(value: u64) -> Self {
        Self { bits: value }
    }

    pub(crate) const fn ulong(self) -> u64 {
        self.bits
    }

    pub(crate) const fn ushort(self, index: usize) -> u16 {
        debug_assert!(index < 2);
        ((self.bits >> (index * 16)) & u16::MAX as u64) as u16
    }
}

pub(crate) type SU_LONG = tagSplitLong;

/// Owned bit representation of `union BnsAltPath`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct BnsAltPath {
    words: [u32; 2],
}

impl BnsAltPath {
    pub(crate) const fn flow(self, index: usize) -> i32 {
        self.words[index] as i32
    }

    pub(crate) const fn set_flow(&mut self, index: usize, value: i32) {
        self.words[index] = value as u32;
    }

    pub(crate) const fn number(self) -> i32 {
        self.words[0] as i32
    }

    pub(crate) const fn set_number(&mut self, value: i32) {
        self.words[0] = value as u32;
    }

    pub(crate) const fn ineigh(self, index: usize) -> u16 {
        debug_assert!(index < 2);
        ((self.words[0] >> (index * 16)) & u16::MAX as u32) as u16
    }

    pub(crate) const fn set_ineigh(&mut self, index: usize, value: u16) {
        let shift = index * 16;
        let mask = (u16::MAX as u32) << shift;
        self.words[0] = (self.words[0] & !mask) | ((value as u32) << shift);
    }
}

pub(crate) type BNS_ALT_PATH = BnsAltPath;

/// Owned state for the local `stbsp__context` definition in `stb_sprintf.h`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct stbsp__context {
    pub(crate) buf: SourceMutPointer<i8>,
    pub(crate) count: i32,
    pub(crate) length: i32,
    pub(crate) tmp: [i8; STB_SPRINTF_MIN as usize],
}

impl Default for stbsp__context {
    fn default() -> Self {
        Self {
            buf: SourceMutPointer::null(),
            count: 0,
            length: 0,
            tmp: [0; STB_SPRINTF_MIN as usize],
        }
    }
}

pub(crate) type UserActionCallback = Option<fn() -> i32>;
pub(crate) type ConsoleQuitCallback = Option<fn() -> i32>;
pub(crate) type SourceCompareCallback =
    Option<fn(SourceConstPointer<SourceVoid>, SourceConstPointer<SourceVoid>) -> i32>;

pub(crate) type CHECK_CENTERPOINT = generated::local_ichiqueu::CHECK_CENTERPOINT;
pub(crate) type CHECK_DFS_RING = generated::local_ichiqueu::CHECK_DFS_RING;
pub(crate) type CHECK_DFS_PATH = generated::local_ichiqueu::CHECK_DFS_PATH;
pub(crate) type CHECK_DFS_CENTERPOINT = generated::local_ichiqueu::CHECK_DFS_CENTERPOINT;

#[cfg(test)]
mod tests {
    use super::*;

    fn user_action() -> i32 {
        17
    }

    fn compare_source_pointers(
        left: SourceConstPointer<SourceVoid>,
        right: SourceConstPointer<SourceVoid>,
    ) -> i32 {
        i32::from(left != right)
    }

    fn check_centerpoint(_atom: SourceMutPointer<inp_ATOM>, atom_index: i32) -> i32 {
        atom_index + 1
    }

    #[test]
    fn inchi_source_types_and_constants() {
        assert_eq!(std::mem::size_of::<S_CHAR>(), 1);
        assert_eq!(std::mem::size_of::<U_CHAR>(), 1);
        assert_eq!(std::mem::size_of::<S_SHORT>(), 2);
        assert_eq!(std::mem::size_of::<U_SHORT>(), 2);
        assert_eq!(std::mem::size_of::<AT_NUM>(), 2);
        assert_eq!(std::mem::size_of::<AT_RANK>(), 2);
        assert_eq!(std::mem::size_of::<INCHI_MODE>(), 8);
        assert_eq!(S_CHAR::MIN, -128);
        assert_eq!(U_CHAR::MAX, 255);
        assert_eq!(S_SHORT::MIN, -32_768);
        assert_eq!(U_SHORT::MAX, 65_535);

        assert_eq!(CURRENT_VER, b"1.07.5\0");
        assert_eq!(MAX_ATOMS, 32_766);
        assert_eq!(MAXVAL, 20);
        assert_eq!(NO_ATOM, -1);
        assert_eq!(INCHIKEY_OK, 0);
        assert_eq!(INCHIKEY_UNKNOWN_ERROR, 1);
        assert_eq!(INCHIKEY_EMPTY_INPUT, 2);
        assert_eq!(INCHIKEY_INVALID_INCHI_PREFIX, 3);
        assert_eq!(INCHIKEY_NOT_ENOUGH_MEMORY, 4);
        assert_eq!(INCHIKEY_INVALID_INCHI, 20);
        assert_eq!(INCHIKEY_INVALID_STD_INCHI, 21);

        let molfile: INPUT_TYPE = tagInputType_INPUT_MOLFILE;
        let unknown_input_type: INPUT_TYPE = 0xfeed_beef;
        assert_eq!(molfile, 1);
        assert_eq!(unknown_input_type, 0xfeed_beef);
        let reverse_wedge: inchi_BondStereo2D = tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2DOWN;
        let unknown_bond_stereo: inchi_BondStereo2D = i32::MIN;
        assert_eq!(reverse_wedge, -6);
        assert_eq!(unknown_bond_stereo, i32::MIN);
        assert_eq!(tagRetValGetINCHI_inchi_Ret_BREAK, -100);
        assert_eq!(tagRetValGetINCHI_inchi_Ret_BUSY, 5);

        let atom = inchi_Atom::default();
        assert_eq!(atom.neighbor, [0; MAXVAL as usize]);
        assert_eq!(atom.bond_type, [0; MAXVAL as usize]);
        assert_eq!(atom.elname, [0; 6]);
        assert_eq!(atom.num_iso_H, [0; 4]);
        assert_eq!((atom.x, atom.y, atom.z), (0.0, 0.0, 0.0));

        let input = inchi_Input::default();
        assert!(input.atom.is_null());
        assert!(input.stereo0D.is_null());
        assert!(input.szOptions.is_null());
        assert_eq!((input.num_atoms, input.num_stereo0D), (0, 0));

        let norm_atom = NORM_ATOM::default();
        assert_eq!(norm_atom.neighbor, [0; MAXVAL as usize]);
        assert_eq!(norm_atom.bond_stereo, [0; MAXVAL as usize]);
        assert_eq!(norm_atom.p_orig_at_num, [0; 4]);
        assert_eq!(norm_atom.sn_orig_at_num, [0; 3]);

        let input_parameters = INPUT_PARMS::default();
        assert_eq!(input_parameters.szSdfDataHeader, [0; 65]);
        assert!(
            input_parameters
                .path
                .iter()
                .all(|pointer| pointer.is_null())
        );
        assert!(input_parameters.pSdfLabel.is_null());
        assert!(input_parameters.pSdfValue.is_null());

        let local_alt_path = local_ichi_bns::ALT_PATH_CHANGES::default();
        assert_eq!(local_alt_path.nOldCapsVert, [[0; MAXVAL as usize + 1]; 2]);
        assert_eq!(local_alt_path.bSetOldCapsVert, [0; 2]);
        let formatting_context = stbsp__context::default();
        assert!(formatting_context.buf.is_null());
        assert_eq!(formatting_context.tmp, [0; STB_SPRINTF_MIN as usize]);

        let null_mut = SourceMutPointer::<i32>::default();
        let null_const = SourceConstPointer::<i32>::default();
        assert!(null_mut.is_null());
        assert!(null_const.is_null());

        let mut heap = SourceHeap::default();
        let allocation = heap.allocate(vec![10_i32, 20, 30]).unwrap();
        let alias = allocation;
        let interior = alias.offset(1).unwrap();
        heap.slice_mut(interior).unwrap()[0] = 25;
        assert_eq!(heap.slice(allocation.as_const()).unwrap(), &[10, 25, 30]);
        assert_eq!(heap.slice(interior.as_const()).unwrap(), &[25, 30]);
        let one_past = allocation.as_const().offset(3).unwrap();
        assert!(heap.slice(one_past).unwrap().is_empty());
        assert_eq!(
            heap.slice(allocation.as_const().offset(4).unwrap()),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            allocation.offset(-1),
            Err(SourceHeapError::PointerOffsetOverflow)
        );
        assert_eq!(
            heap.free(interior),
            Err(SourceHeapError::FreeOfInteriorPointer)
        );
        assert_eq!(heap.free(SourceMutPointer::<i32>::null()), Ok(()));
        assert_eq!(
            heap.slice(SourceConstPointer::<i32>::null()),
            Err(SourceHeapError::NullPointer)
        );
        heap.free(allocation).unwrap();
        assert_eq!(
            heap.slice(alias.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let bytes = heap.allocate(vec![1_u8, 2]).unwrap();
        let wrong_type = SourceConstPointer::<u16> {
            allocation: bytes.allocation,
            element_offset: 0,
            marker: PhantomData,
        };
        assert_eq!(
            heap.slice(wrong_type),
            Err(SourceHeapError::AllocationTypeMismatch)
        );
        heap.free(bytes).unwrap();

        let split = tagSplitLong::from_ulong(0x1122_3344_5566_7788);
        assert_eq!(split.ulong(), 0x1122_3344_5566_7788);
        assert_eq!(split.ushort(0), 0x7788);
        assert_eq!(split.ushort(1), 0x5566);

        let mut alt_path = BnsAltPath::default();
        alt_path.set_flow(1, -7);
        assert_eq!(alt_path.flow(1), -7);
        alt_path.set_number(0x1234_5678);
        assert_eq!(alt_path.number(), 0x1234_5678);
        assert_eq!(alt_path.ineigh(0), 0x5678);
        assert_eq!(alt_path.ineigh(1), 0x1234);
        alt_path.set_ineigh(1, 0xabcd);
        assert_eq!(alt_path.number() as u32, 0xabcd_5678);

        let action: UserActionCallback = Some(user_action);
        let quit: ConsoleQuitCallback = Some(user_action);
        let compare: SourceCompareCallback = Some(compare_source_pointers);
        let centerpoint: CHECK_CENTERPOINT = Some(check_centerpoint);
        assert_eq!(action.unwrap()(), 17);
        assert_eq!(quit.unwrap()(), 17);
        assert_eq!(
            compare.unwrap()(SourceConstPointer::null(), SourceConstPointer::null()),
            0
        );
        assert_eq!(centerpoint.unwrap()(SourceMutPointer::null(), 8), 9);

        assert_eq!(local_readinch::ISOLATED_ATOM, 15_u32);
        assert_eq!(local_inchi_dll_b::ISOLATED_ATOM, -15_i32);
        assert_eq!(local_readinch::MAX_CHAIN_LEN, 20_u32);
        assert_eq!(local_inchi_dll_b::MAX_CHAIN_LEN, 20_u32);
    }
}
