#![allow(
    dead_code,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_imports
)]

use std::any::TypeId;
#[cfg(test)]
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::marker::PhantomData;
use std::ptr::NonNull;

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

/// Identity of a `T_GROUP_INFO *` output selected during reverse-InChI work.
///
/// Active production paths assign either `NULL` or `&pStruct->One_ti`.
/// `External` retains an existing caller pointer identity when a source error
/// returns before assigning the output.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum SourceTGroupInfoPointer {
    #[default]
    Null,
    StructureOne,
    External(u64),
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

            #[inline(always)]
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

            /// Advances a source pointer after the owning port has proved the
            /// complete C allocation range once.
            ///
            /// # Safety
            ///
            /// `elements` must not overflow the element offset, and the
            /// resulting pointer must remain within the same live allocation.
            #[inline(always)]
            pub(crate) unsafe fn add_unchecked(self, elements: u64) -> Self {
                debug_assert!(self.element_offset.checked_add(elements).is_some());
                Self {
                    element_offset: self.element_offset.wrapping_add(elements),
                    ..self
                }
            }

            /// Returns the element offset from the allocation base after a
            /// source constructor has proved that the pointer is live.
            ///
            /// # Safety
            ///
            /// The pointer must address a live allocation and its element
            /// offset must fit in `usize`.
            #[inline(always)]
            pub(crate) unsafe fn allocation_offset_unchecked(self) -> usize {
                debug_assert!(self.allocation.is_some());
                debug_assert!(usize::try_from(self.element_offset).is_ok());
                self.element_offset as usize
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

    pub(crate) const fn cast<U>(self) -> SourceMutPointer<U> {
        SourceMutPointer {
            allocation: self.allocation,
            element_offset: self.element_offset,
            marker: PhantomData,
        }
    }

    #[inline(always)]
    pub(crate) fn difference(self, origin: Self) -> Result<i64, SourceHeapError> {
        let allocation = self.allocation.ok_or(SourceHeapError::NullPointer)?;
        if origin.allocation != Some(allocation) {
            return Err(SourceHeapError::PointerAllocationMismatch);
        }
        if self.element_offset >= origin.element_offset {
            let difference = self.element_offset - origin.element_offset;
            i64::try_from(difference).map_err(|_| SourceHeapError::PointerDifferenceOverflow)
        } else {
            let difference = origin.element_offset - self.element_offset;
            let difference = i64::try_from(difference).map_err(|_| SourceHeapError::PointerDifferenceOverflow)?;
            Ok(-difference)
        }
    }

    /// Returns the non-negative C element offset after a caller has proved
    /// that both pointers address the same live allocation and that `self`
    /// does not precede `origin`.
    ///
    /// # Safety
    ///
    /// The pointers must have the same non-null allocation identity,
    /// `self.element_offset >= origin.element_offset`, and the difference must
    /// fit in `usize`.
    #[inline(always)]
    pub(crate) unsafe fn forward_difference_unchecked(self, origin: Self) -> usize {
        debug_assert!(self.allocation.is_some());
        debug_assert_eq!(self.allocation, origin.allocation);
        debug_assert!(self.element_offset >= origin.element_offset);
        let difference = self.element_offset - origin.element_offset;
        debug_assert!(usize::try_from(difference).is_ok());
        difference as usize
    }

    pub(crate) const fn allocation_identity(self) -> Option<u64> {
        match self.allocation {
            Some(id) => Some(id.0),
            None => None,
        }
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
struct AllocationArena {
    slots: Vec<Option<AllocationSlot>>,
    live_count: usize,
}

struct AllocationSlot {
    pointer: NonNull<()>,
    len: usize,
    capacity: usize,
    type_id: TypeId,
    drop_vec: unsafe fn(NonNull<()>, usize, usize),
    provenance: AllocationProvenance,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum AllocationProvenance {
    #[default]
    None,
    ContiguousNeighborPointerTable {
        storage: AllocationId,
        row_count: usize,
        neighbor_index_bound: usize,
    },
    IndexBound {
        prefix_len: usize,
        upper_bound: usize,
    },
}

impl AllocationSlot {
    fn new<T: 'static>(mut value: Vec<T>) -> Self {
        let pointer = NonNull::new(value.as_mut_ptr())
            .expect("Vec pointers, including dangling empty pointers, are non-null")
            .cast();
        let len = value.len();
        let capacity = value.capacity();
        std::mem::forget(value);
        Self {
            pointer,
            len,
            capacity,
            type_id: TypeId::of::<Vec<T>>(),
            drop_vec: drop_erased_vec::<T>,
            provenance: AllocationProvenance::None,
        }
    }

    #[inline(always)]
    fn downcast_ref<T: 'static>(&self) -> Option<&[T]> {
        if self.type_id != TypeId::of::<Vec<T>>() {
            return None;
        }
        // SAFETY: `type_id` and the raw parts are recorded from the same Vec,
        // and the allocation remains owned by this slot.
        Some(unsafe { std::slice::from_raw_parts(self.pointer.cast::<T>().as_ptr(), self.len) })
    }

    #[inline(always)]
    fn downcast_mut<T: 'static>(&mut self) -> Option<&mut [T]> {
        if self.type_id != TypeId::of::<Vec<T>>() {
            return None;
        }
        // A mutable view may overwrite source pointers whose construction
        // established allocation relationships recorded below.
        self.provenance = AllocationProvenance::None;
        // SAFETY: as above; `&mut self` also guarantees unique access.
        Some(unsafe { std::slice::from_raw_parts_mut(self.pointer.cast::<T>().as_ptr(), self.len) })
    }

    fn is<T: 'static>(&self) -> bool {
        self.type_id == TypeId::of::<Vec<T>>()
    }
}

impl Drop for AllocationSlot {
    fn drop(&mut self) {
        // SAFETY: the matching function was installed while erasing this Vec.
        unsafe { (self.drop_vec)(self.pointer, self.len, self.capacity) };
    }
}

unsafe fn drop_erased_vec<T>(pointer: NonNull<()>, len: usize, capacity: usize) {
    // SAFETY: the caller supplies the unchanged raw parts of a Vec<T>.
    drop(unsafe { Vec::from_raw_parts(pointer.cast::<T>().as_ptr(), len, capacity) });
}

/// A validated view into an allocation whose `Vec` buffer is stable.
///
/// This is intentionally not `Clone` or `Copy`: callers that create writable
/// views must prove that their allocation identities are distinct. The heap
/// never resizes allocation buffers, so the view remains valid until its
/// allocation is freed.
pub(crate) struct StableSourceSlice<T> {
    pointer: NonNull<T>,
    len: usize,
}

/// A validated read-only view into an allocation whose `Vec` buffer is stable.
///
/// Unlike an ordinary heap slice, this view does not keep the heap borrowed.
/// It is intended for source loops that read one allocation while mutating
/// separate allocations through the heap.
pub(crate) struct StableSourceConstSlice<T> {
    pointer: NonNull<T>,
    len: usize,
}

impl<T> StableSourceConstSlice<T> {
    #[inline(always)]
    pub(crate) const fn len(&self) -> usize {
        self.len
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn get(&self, index: usize) -> Result<&T, SourceHeapError> {
        if index >= self.len {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // SAFETY: construction validates the allocation type and bounds. The
        // caller keeps the allocation live while this view is in use.
        Ok(unsafe { &*self.pointer.as_ptr().add(index) })
    }

    /// Returns an element after the owning source port has proved its C array
    /// contract once at workspace construction.
    ///
    /// # Safety
    ///
    /// `index` must be smaller than `self.len`, and the allocation must remain
    /// live without overlapping mutation for the returned borrow.
    #[inline(always)]
    pub(crate) unsafe fn get_unchecked(&self, index: usize) -> &T {
        // SAFETY: upheld by the caller's source-array contract.
        unsafe { &*self.pointer.as_ptr().add(index) }
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn prefix(&self, len: usize) -> Result<&[T], SourceHeapError> {
        if len > self.len {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // SAFETY: the requested prefix is within the validated allocation.
        Ok(unsafe { std::slice::from_raw_parts(self.pointer.as_ptr(), len) })
    }
}

impl<T> StableSourceSlice<T> {
    #[inline(always)]
    pub(crate) const fn len(&self) -> usize {
        self.len
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn get(&self, index: usize) -> Result<&T, SourceHeapError> {
        if index >= self.len {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // SAFETY: construction validates the allocation type and bounds. The
        // caller of `stable_slice_mut` keeps the allocation live and unique.
        Ok(unsafe { &*self.pointer.as_ptr().add(index) })
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn get_mut(&mut self, index: usize) -> Result<&mut T, SourceHeapError> {
        if index >= self.len {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // SAFETY: as above, with uniqueness enforced by the workspace that
        // owns this non-cloneable handle.
        Ok(unsafe { &mut *self.pointer.as_ptr().add(index) })
    }

    /// Returns an element after the owning source port has proved its C array
    /// contract once at workspace construction.
    ///
    /// # Safety
    ///
    /// `index` must be smaller than `self.len`, and mutable handles to the
    /// allocation must remain confined to the workspace that owns this view.
    #[inline(always)]
    pub(crate) unsafe fn get_unchecked(&self, index: usize) -> &T {
        // SAFETY: upheld by the caller's source-array contract.
        unsafe { &*self.pointer.as_ptr().add(index) }
    }

    /// Returns a writable element after the owning source port has proved its
    /// C array contract once at workspace construction.
    ///
    /// # Safety
    ///
    /// `index` must be smaller than `self.len`, and this handle must be the
    /// unique writable view of the element for the returned borrow.
    #[inline(always)]
    pub(crate) unsafe fn get_unchecked_mut(&mut self, index: usize) -> &mut T {
        // SAFETY: upheld by the caller's source-array and alias contracts.
        unsafe { &mut *self.pointer.as_ptr().add(index) }
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn prefix_mut(&mut self, len: usize) -> Result<&mut [T], SourceHeapError> {
        if len > self.len {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        // SAFETY: the requested prefix is within the validated allocation and
        // the mutable handle is unique for the duration of this operation.
        Ok(unsafe { std::slice::from_raw_parts_mut(self.pointer.as_ptr(), len) })
    }

    #[inline(always)]
    pub(crate) fn fill(&mut self, value: T)
    where
        T: Clone,
    {
        // SAFETY: construction validated the complete allocation suffix and
        // the mutable handle is unique for the duration of this operation.
        unsafe { std::slice::from_raw_parts_mut(self.pointer.as_ptr(), self.len) }.fill(value);
    }
}

impl AllocationArena {
    fn index(id: AllocationId) -> Result<usize, SourceHeapError> {
        usize::try_from(id.0).map_err(|_| SourceHeapError::MissingAllocation)
    }

    fn insert<T: 'static>(&mut self, id: AllocationId, allocation: Vec<T>) -> Result<(), SourceHeapError> {
        let index = usize::try_from(id.0).map_err(|_| SourceHeapError::AllocationIdExhausted)?;
        if index != self.slots.len() {
            return Err(SourceHeapError::AllocationIdExhausted);
        }
        self.slots
            .try_reserve(1)
            .map_err(|_| SourceHeapError::AllocationFailed)?;
        self.slots.push(Some(AllocationSlot::new(allocation)));
        self.live_count += 1;
        Ok(())
    }

    fn get(&self, id: AllocationId) -> Option<&AllocationSlot> {
        self.slots.get(Self::index(id).ok()?)?.as_ref()
    }

    fn get_mut(&mut self, id: AllocationId) -> Option<&mut AllocationSlot> {
        self.slots.get_mut(Self::index(id).ok()?)?.as_mut()
    }

    fn take(&mut self, id: AllocationId) -> Option<AllocationSlot> {
        let allocation = self.slots.get_mut(Self::index(id).ok()?)?.take()?;
        self.live_count -= 1;
        Some(allocation)
    }

    fn restore(&mut self, id: AllocationId, allocation: AllocationSlot) {
        let index = Self::index(id).expect("an allocated ID fits usize");
        let slot = self.slots.get_mut(index).expect("allocated slot exists");
        debug_assert!(slot.is_none());
        *slot = Some(allocation);
        self.live_count += 1;
    }

    fn remove(&mut self, id: AllocationId) -> Option<AllocationSlot> {
        self.take(id)
    }

    fn len(&self) -> usize {
        self.live_count
    }

    #[cfg(test)]
    fn values(&self) -> impl Iterator<Item = &AllocationSlot> {
        self.slots.iter().filter_map(Option::as_ref)
    }
}

#[derive(Default)]
pub(crate) struct SourceHeap {
    next_id: u64,
    allocations: AllocationArena,
    source_errno: i32,
    #[cfg(test)]
    allocations_before_failure: Option<u64>,
    #[cfg(test)]
    source_allocation_calls: u64,
    #[cfg(test)]
    live_source_allocations: BTreeSet<AllocationId>,
    #[cfg(test)]
    live_source_allocation_types: BTreeMap<AllocationId, &'static str>,
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

fn overwrite_inp_atom_from_gcc_lp64_bytes(atom: &mut inp_ATOM, bytes: &[u8; INP_ATOM_GCC_LP64_SIZE]) {
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
    pub(crate) fn live_source_allocation_count(&self) -> usize {
        self.live_source_allocations.len()
    }

    #[cfg(test)]
    pub(crate) fn live_source_allocations_of<T: 'static>(&self) -> usize {
        self.live_source_allocations
            .iter()
            .filter(|id| {
                self.allocations
                    .get(**id)
                    .is_some_and(|allocation| allocation.is::<T>())
            })
            .count()
    }

    #[cfg(test)]
    pub(crate) fn live_source_allocation_types(&self) -> Vec<&'static str> {
        self.live_source_allocation_types.values().copied().collect()
    }

    #[cfg(test)]
    pub(crate) fn live_allocations_of<T: 'static>(&self) -> usize {
        self.allocations
            .values()
            .filter(|allocation| allocation.is::<T>())
            .count()
    }

    pub(crate) fn allocate<T: 'static>(&mut self, values: Vec<T>) -> Result<SourceMutPointer<T>, SourceHeapError> {
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
        let pointer = self.allocate_model_storage(values)?;
        #[cfg(test)]
        self.live_source_allocations
            .insert(pointer.allocation.expect("fresh allocation is non-null"));
        #[cfg(test)]
        self.live_source_allocation_types.insert(
            pointer.allocation.expect("fresh allocation is non-null"),
            std::any::type_name::<T>(),
        );
        Ok(pointer)
    }

    pub(crate) fn allocate_model_storage<T: 'static>(
        &mut self,
        values: Vec<T>,
    ) -> Result<SourceMutPointer<T>, SourceHeapError> {
        let id = AllocationId(self.next_id);
        let next_id = self
            .next_id
            .checked_add(1)
            .ok_or(SourceHeapError::AllocationIdExhausted)?;
        self.allocations.insert(id, values)?;
        self.next_id = next_id;
        Ok(SourceMutPointer {
            allocation: Some(id),
            element_offset: 0,
            marker: PhantomData,
        })
    }

    /// Records the contiguous pointer-table layout produced by the two
    /// Official InChI `CreateNeighList*` constructors.
    pub(crate) fn record_contiguous_neighbor_layout<P: 'static, S: 'static>(
        &mut self,
        pointer_table: SourceMutPointer<P>,
        storage: SourceMutPointer<S>,
        row_count: usize,
        neighbor_index_bound: usize,
    ) -> Result<(), SourceHeapError> {
        if pointer_table.element_offset != 0 || storage.element_offset != 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let table_id = pointer_table.allocation.ok_or(SourceHeapError::NullPointer)?;
        let storage_id = storage.allocation.ok_or(SourceHeapError::NullPointer)?;
        if table_id == storage_id {
            return Err(SourceHeapError::PointerAllocationMismatch);
        }
        let storage_slot = self
            .allocations
            .get(storage_id)
            .ok_or(SourceHeapError::MissingAllocation)?;
        if !storage_slot.is::<S>() {
            return Err(SourceHeapError::AllocationTypeMismatch);
        }
        let table_slot = self
            .allocations
            .get_mut(table_id)
            .ok_or(SourceHeapError::MissingAllocation)?;
        if !table_slot.is::<P>() {
            return Err(SourceHeapError::AllocationTypeMismatch);
        }
        if table_slot.len < row_count {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        table_slot.provenance = AllocationProvenance::ContiguousNeighborPointerTable {
            storage: storage_id,
            row_count,
            neighbor_index_bound,
        };
        Ok(())
    }

    /// Returns the base storage pointer when a live source constructor proved
    /// the complete contiguous neighbor layout. Manually assembled pointer
    /// tables deliberately return `None` and retain their checked path.
    pub(crate) fn proven_contiguous_neighbor_storage<P: 'static, S: 'static>(
        &self,
        pointer_table: SourceConstPointer<P>,
        requested_rows: usize,
        available_neighbor_values: usize,
    ) -> Option<SourceMutPointer<S>> {
        if pointer_table.element_offset != 0 {
            return None;
        }
        let table_id = pointer_table.allocation?;
        let table_slot = self.allocations.get(table_id)?;
        if !table_slot.is::<P>() || table_slot.len < requested_rows {
            return None;
        }
        let AllocationProvenance::ContiguousNeighborPointerTable {
            storage,
            row_count,
            neighbor_index_bound,
        } = table_slot.provenance
        else {
            return None;
        };
        if row_count < requested_rows || neighbor_index_bound > available_neighbor_values {
            return None;
        }
        let storage_slot = self.allocations.get(storage)?;
        if !storage_slot.is::<S>() {
            return None;
        }
        Some(SourceMutPointer {
            allocation: Some(storage),
            element_offset: 0,
            marker: PhantomData,
        })
    }

    /// Records a source-level proof that every value in the prefix is smaller
    /// than `upper_bound`. Callers establish this either by direct source
    /// initialization or by a complete checked scan.
    pub(crate) fn record_index_bound<T: 'static>(
        &mut self,
        pointer: SourceMutPointer<T>,
        prefix_len: usize,
        upper_bound: usize,
    ) -> Result<(), SourceHeapError> {
        if pointer.element_offset != 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let slot = self.allocations.get_mut(id).ok_or(SourceHeapError::MissingAllocation)?;
        if !slot.is::<T>() {
            return Err(SourceHeapError::AllocationTypeMismatch);
        }
        if slot.len < prefix_len {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        slot.provenance = AllocationProvenance::IndexBound {
            prefix_len,
            upper_bound,
        };
        Ok(())
    }

    pub(crate) fn has_proven_index_bound<T: 'static>(
        &self,
        pointer: SourceConstPointer<T>,
        requested_prefix_len: usize,
        available_values: usize,
    ) -> bool {
        if pointer.element_offset != 0 {
            return false;
        }
        let Some(id) = pointer.allocation else {
            return false;
        };
        let Some(slot) = self.allocations.get(id) else {
            return false;
        };
        if !slot.is::<T>() || slot.len < requested_prefix_len {
            return false;
        }
        matches!(
            slot.provenance,
            AllocationProvenance::IndexBound {
                prefix_len,
                upper_bound,
            } if prefix_len >= requested_prefix_len && upper_bound <= available_values
        )
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn slice<T: 'static>(&self, pointer: SourceConstPointer<T>) -> Result<&[T], SourceHeapError> {
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let values = self
            .allocations
            .get(id)
            .ok_or(SourceHeapError::MissingAllocation)?
            .downcast_ref::<T>()
            .ok_or(SourceHeapError::AllocationTypeMismatch)?;
        let offset = usize::try_from(pointer.element_offset).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        values.get(offset..).ok_or(SourceHeapError::PointerOutOfBounds)
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn slice_mut<T: 'static>(&mut self, pointer: SourceMutPointer<T>) -> Result<&mut [T], SourceHeapError> {
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let values = self
            .allocations
            .get_mut(id)
            .ok_or(SourceHeapError::MissingAllocation)?
            .downcast_mut::<T>()
            .ok_or(SourceHeapError::AllocationTypeMismatch)?;
        let offset = usize::try_from(pointer.element_offset).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        values.get_mut(offset..).ok_or(SourceHeapError::PointerOutOfBounds)
    }

    /// Validates a source pointer once for a tightly scoped read-heavy loop.
    ///
    /// # Safety
    ///
    /// The allocation must remain live and its buffer must not be resized
    /// while the returned view is used. Any mutable access to the allocation
    /// must not overlap a reference returned by `get`.
    pub(crate) unsafe fn stable_slice<T: 'static>(
        &self,
        pointer: SourceConstPointer<T>,
    ) -> Result<StableSourceConstSlice<T>, SourceHeapError> {
        let values = self.slice(pointer)?;
        Ok(StableSourceConstSlice {
            pointer: NonNull::new(values.as_ptr().cast_mut()).expect("slice pointers are non-null"),
            len: values.len(),
        })
    }

    /// Validates a mutable source pointer once for a tightly scoped hot loop.
    ///
    /// # Safety
    ///
    /// The allocation must remain live while the returned view is used. Any
    /// other writable stable view used at the same time must refer to a
    /// different allocation, and ordinary heap borrows of this allocation may
    /// not overlap a reference returned by the view.
    pub(crate) unsafe fn stable_slice_mut<T: 'static>(
        &mut self,
        pointer: SourceMutPointer<T>,
    ) -> Result<StableSourceSlice<T>, SourceHeapError> {
        let values = self.slice_mut(pointer)?;
        Ok(StableSourceSlice {
            pointer: NonNull::new(values.as_mut_ptr()).expect("slice pointers are non-null"),
            len: values.len(),
        })
    }

    /// Creates a mutable view for a source operation that only permutes a
    /// prefix whose index bound has already been proved.
    ///
    /// # Safety
    ///
    /// In addition to the requirements of `stable_slice_mut`, every write
    /// through the returned view must preserve the recorded prefix and bound.
    pub(crate) unsafe fn stable_index_bounded_slice_mut<T: 'static>(
        &mut self,
        pointer: SourceMutPointer<T>,
        prefix_len: usize,
        upper_bound: usize,
    ) -> Result<StableSourceSlice<T>, SourceHeapError> {
        if !self.has_proven_index_bound(pointer.as_const(), prefix_len, upper_bound) {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let slot = self.allocations.get_mut(id).ok_or(SourceHeapError::MissingAllocation)?;
        let values = slot
            .downcast_ref::<T>()
            .ok_or(SourceHeapError::AllocationTypeMismatch)?;
        let offset = usize::try_from(pointer.element_offset).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let values = values.get(offset..).ok_or(SourceHeapError::PointerOutOfBounds)?;
        Ok(StableSourceSlice {
            pointer: NonNull::new(values.as_ptr().cast_mut()).expect("slice pointers are non-null"),
            len: values.len(),
        })
    }

    #[track_caller]
    pub(crate) fn with_slice_mut_and_heap<T: 'static, R>(
        &mut self,
        pointer: SourceMutPointer<T>,
        operation: impl FnOnce(&mut [T], &SourceHeap) -> Result<R, SourceHeapError>,
    ) -> Result<R, SourceHeapError> {
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let mut allocation = self.allocations.remove(id).ok_or(SourceHeapError::MissingAllocation)?;
        let offset = match usize::try_from(pointer.element_offset) {
            Ok(offset) => offset,
            Err(_) => {
                self.allocations.restore(id, allocation);
                return Err(SourceHeapError::PointerOutOfBounds);
            }
        };
        let result = match allocation.downcast_mut::<T>() {
            Some(values) => match values.get_mut(offset..) {
                Some(values) => operation(values, self),
                None => Err(SourceHeapError::PointerOutOfBounds),
            },
            None => Err(SourceHeapError::AllocationTypeMismatch),
        };
        self.allocations.restore(id, allocation);
        result
    }

    #[track_caller]
    pub(crate) fn with_slice_mut_and_heap_mut<T: 'static, R>(
        &mut self,
        pointer: SourceMutPointer<T>,
        operation: impl FnOnce(&mut [T], &mut SourceHeap) -> Result<R, SourceHeapError>,
    ) -> Result<R, SourceHeapError> {
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let mut allocation = self.allocations.remove(id).ok_or(SourceHeapError::MissingAllocation)?;
        let offset = match usize::try_from(pointer.element_offset) {
            Ok(offset) => offset,
            Err(_) => {
                self.allocations.restore(id, allocation);
                return Err(SourceHeapError::PointerOutOfBounds);
            }
        };
        let result = match allocation.downcast_mut::<T>() {
            Some(values) => match values.get_mut(offset..) {
                Some(values) => operation(values, self),
                None => Err(SourceHeapError::PointerOutOfBounds),
            },
            None => Err(SourceHeapError::AllocationTypeMismatch),
        };
        self.allocations.restore(id, allocation);
        result
    }

    pub(crate) fn with_two_slices_mut_and_optional_third<A: 'static, B: 'static, C: 'static, R>(
        &mut self,
        first: SourceMutPointer<A>,
        second: SourceMutPointer<B>,
        third: Option<SourceMutPointer<C>>,
        operation: impl FnOnce(&mut [A], &mut [B], Option<&mut [C]>) -> Result<R, SourceHeapError>,
    ) -> Result<R, SourceHeapError> {
        let first_id = first.allocation.ok_or(SourceHeapError::NullPointer)?;
        let second_id = second.allocation.ok_or(SourceHeapError::NullPointer)?;
        let third_id = third.and_then(|pointer| pointer.allocation);
        if first_id == second_id || third_id.is_some_and(|id| id == first_id || id == second_id) {
            return Err(SourceHeapError::PointerAllocationMismatch);
        }

        let mut first_allocation = self
            .allocations
            .remove(first_id)
            .ok_or(SourceHeapError::MissingAllocation)?;
        let Some(mut second_allocation) = self.allocations.remove(second_id) else {
            self.allocations.restore(first_id, first_allocation);
            return Err(SourceHeapError::MissingAllocation);
        };
        let mut third_allocation = if let Some(id) = third_id {
            match self.allocations.remove(id) {
                Some(allocation) => Some((id, allocation)),
                None => {
                    self.allocations.restore(first_id, first_allocation);
                    self.allocations.restore(second_id, second_allocation);
                    return Err(SourceHeapError::MissingAllocation);
                }
            }
        } else {
            None
        };

        let result = (|| {
            let first_offset =
                usize::try_from(first.element_offset).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let second_offset =
                usize::try_from(second.element_offset).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let first_values = first_allocation
                .downcast_mut::<A>()
                .ok_or(SourceHeapError::AllocationTypeMismatch)?
                .get_mut(first_offset..)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let second_values = second_allocation
                .downcast_mut::<B>()
                .ok_or(SourceHeapError::AllocationTypeMismatch)?
                .get_mut(second_offset..)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let third_values = match (&mut third_allocation, third) {
                (Some((_, allocation)), Some(pointer)) => {
                    let offset =
                        usize::try_from(pointer.element_offset).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    Some(
                        allocation
                            .downcast_mut::<C>()
                            .ok_or(SourceHeapError::AllocationTypeMismatch)?
                            .get_mut(offset..)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    )
                }
                _ => None,
            };
            operation(first_values, second_values, third_values)
        })();

        self.allocations.restore(first_id, first_allocation);
        self.allocations.restore(second_id, second_allocation);
        if let Some((id, allocation)) = third_allocation {
            self.allocations.restore(id, allocation);
        }
        result
    }

    pub(crate) fn free<T: 'static>(&mut self, pointer: SourceMutPointer<T>) -> Result<(), SourceHeapError> {
        let Some(id) = pointer.allocation else {
            return Ok(());
        };
        if pointer.element_offset != 0 {
            return Err(SourceHeapError::FreeOfInteriorPointer);
        }
        let allocation = self.allocations.get(id).ok_or(SourceHeapError::MissingAllocation)?;
        if !allocation.is::<T>() {
            return Err(SourceHeapError::AllocationTypeMismatch);
        }
        self.allocations.remove(id);
        #[cfg(test)]
        self.live_source_allocations.remove(&id);
        #[cfg(test)]
        self.live_source_allocation_types.remove(&id);
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

    fn compare_source_pointers(left: SourceConstPointer<SourceVoid>, right: SourceConstPointer<SourceVoid>) -> i32 {
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
        assert!(input_parameters.path.iter().all(|pointer| pointer.is_null()));
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
        assert_eq!(allocation.offset(-1), Err(SourceHeapError::PointerOffsetOverflow));
        assert_eq!(heap.free(interior), Err(SourceHeapError::FreeOfInteriorPointer));
        assert_eq!(heap.free(SourceMutPointer::<i32>::null()), Ok(()));
        assert_eq!(
            heap.slice(SourceConstPointer::<i32>::null()),
            Err(SourceHeapError::NullPointer)
        );
        heap.free(allocation).unwrap();
        assert_eq!(heap.slice(alias.as_const()), Err(SourceHeapError::MissingAllocation));

        let bytes = heap.allocate(vec![1_u8, 2]).unwrap();
        let wrong_type = SourceConstPointer::<u16> {
            allocation: bytes.allocation,
            element_offset: 0,
            marker: PhantomData,
        };
        assert_eq!(heap.slice(wrong_type), Err(SourceHeapError::AllocationTypeMismatch));
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

    #[test]
    fn index_bound_proof_survives_only_bound_preserving_mutation() {
        let mut heap = SourceHeap::default();
        let values = heap.allocate(vec![2_u16, 0, 1]).unwrap();
        heap.record_index_bound(values, 3, 3).unwrap();

        assert!(heap.has_proven_index_bound(values.as_const(), 3, 3));
        assert!(heap.has_proven_index_bound(values.as_const(), 2, 4));
        assert!(!heap.has_proven_index_bound(values.as_const(), 3, 2));

        // SAFETY: swapping two entries preserves the recorded prefix bound.
        let mut bounded = unsafe { heap.stable_index_bounded_slice_mut(values, 3, 3).unwrap() };
        let first = *bounded.get(0).unwrap();
        let second = *bounded.get(1).unwrap();
        *bounded.get_mut(0).unwrap() = second;
        *bounded.get_mut(1).unwrap() = first;
        drop(bounded);
        assert!(heap.has_proven_index_bound(values.as_const(), 3, 3));

        heap.slice_mut(values).unwrap()[0] = 9;
        assert!(!heap.has_proven_index_bound(values.as_const(), 3, 3));
    }

    #[test]
    fn stable_const_view_observes_source_equivalent_in_place_refresh() {
        let mut heap = SourceHeap::default();
        let values = heap.allocate(vec![10_i32, 20]).unwrap();

        // SAFETY: the allocation remains live and fixed-capacity. Each read
        // borrow ends before the source-equivalent in-place refresh below.
        let view = unsafe { heap.stable_slice(values.as_const()).unwrap() };
        assert_eq!(*view.get(0).unwrap(), 10);

        heap.slice_mut(values).unwrap().copy_from_slice(&[30, 40]);

        assert_eq!(*view.get(0).unwrap(), 30);
        assert_eq!(*view.get(1).unwrap(), 40);
    }
}
