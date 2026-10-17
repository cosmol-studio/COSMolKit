#![allow(
    dead_code,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_imports
)]

use std::any::{Any, TypeId};
use std::collections::HashMap;
#[cfg(test)]
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::marker::PhantomData;
use std::num::NonZeroU64;
use std::ptr::NonNull;

mod generated;
mod clock;

pub(crate) use generated::*;
pub(crate) use clock::{LiveSourceClock, SourceClock};

/// Source-level replacement for C `void` in typed allocation handles.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct SourceVoid;

/// Owned byte stream replacing the `FILE *` branch of `INCHI_IOSTREAM`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct SourceFile {
    pub(crate) bytes: Vec<u8>,
    pub(crate) position: u64,
    pub(crate) error: bool,
    /// Sticky stdio F_ERR; `error` retains the existing forced IO failure.
    pub(crate) error_indicator: bool,
    pub(crate) eof: bool,
    pub(crate) is_standard_stream: bool,
}

// The extra private flag occupies existing padding, not another allocation.
const _: () = assert!(std::mem::size_of::<SourceFile>() <= 40);

impl SourceFile {
    pub(crate) fn read_byte(&mut self) -> Result<i32, SourceHeapError> {
        // SUPPLEMENTAL SOURCE: musl/src/stdio/__uflow.c9..10
        // 	if (!__toread(f) && f->read(f, &c, 1)==1) return c;
        // 	return EOF;
        // SUPPLEMENTAL SOURCE: musl/src/stdio/__stdio_read.c14..17
        // 	if (cnt <= 0) {
        // 		f->flags |= cnt ? F_ERR : F_EOF;
        // 		return 0;
        // 	}
        // SUPPLEMENTAL SOURCE: musl/src/stdio/__toread.c13
        // 	return (f->flags & F_EOF) ? EOF : 0;
        // INCHI✔️❗: Existing unbuffered owned-byte backend: a prior sticky
        // indicator never blocks a read; a new forced IO failure sets it and
        // consumes no byte. The native caller independently applies ferror.
        // O(1), no allocation/clone; FILE-ERROR-STATE whole costs stay OPEN.
        let index = usize::try_from(self.position)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if self.eof {
            return Ok(-1);
        }
        if self.error {
            self.error_indicator = true;
            return Ok(-1);
        }
        if let Some(value) = self.bytes.get(index).copied() {
            self.position = self.position.wrapping_add(1);
            Ok(i32::from(value))
        } else {
            self.eof = true;
            Ok(-1)
        }
    }
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
    /// Actual Rust address for diagnostic %p only; never dereferenced.
    Pointer(*const SourceVoid),
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

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct AllocationId {
    owner: NonZeroU64,
    field: u64,
}

impl AllocationId {
    const fn logical(self) -> u64 {
        self.owner.get() - 1
    }

    fn subobject_field(self, row: u64, kind: u64) -> Result<Self, SourceHeapError> {
        // VERBATIM SOURCE: ichicano.c1147 / ichiprt1.c2635
        //         pCG->m_pNeighborsForSort = at[i].neighbor;
        //             pCG->m_pn_RankForSort = at[j].neighbor;
        // INCHI✔️❗: A subobject descriptor, not an allocation or an address.
        // Both owner structs exceed four bytes, so every realizable Vec row
        // fits below 2^62. Preserve a separate relative element offset.
        // O(1), no allocation/copy; INCHI-PERF-REAL-ATOM-FIELD-VIEWS OPEN.
        if self.field != 0 || row >= (1_u64 << 62) {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        Ok(Self { field: (kind << 62) | row, ..self })
    }
}

impl fmt::Debug for AllocationId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.field == 0 {
            formatter.debug_tuple("AllocationId").field(&self.logical()).finish()
        } else {
            formatter.debug_struct("AllocationId")
                .field("owner", &self.logical()).field("field", &self.field).finish()
        }
    }
}

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

// Field metadata occupies the old Option discriminant space. This is a
// private footprint invariant, not a new test or a claim about the C ABI.
const _: () = assert!(std::mem::size_of::<SourceConstPointer<u16>>()
    <= std::mem::size_of::<(Option<u64>, u64)>());

macro_rules! impl_neighbor_field_pointer {
    ($owner:ty, $kind:expr) => {
        impl SourceConstPointer<$owner> {
            pub(crate) fn neighbor_field(self, index: i32)
                -> Result<SourceConstPointer<AT_NUMB>, SourceHeapError>
            {
                let pointer = self.offset(i64::from(index))?;
                let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
                Ok(SourceConstPointer {
                    allocation: Some(id.subobject_field(pointer.element_offset, $kind)?),
                    element_offset: 0,
                    marker: PhantomData,
                })
            }
        }
    };
}

impl_neighbor_field_pointer!(sp_ATOM, 1);
impl_neighbor_field_pointer!(inp_ATOM, 2);

impl SourceConstPointer<INCHIGEN_CONTROL> {
    pub(crate) fn input_parameters_field(self) -> Result<SourceMutPointer<INPUT_PARMS>, SourceHeapError> {
        // SOURCE: inchi_dll_a.c INCHIGEN_Setup
        //     INPUT_PARMS *ip = &( HGen->InpParms );
        // INCHI✔️❗: Actual parent field, not a copied/model allocation.
        // O(1), no extra pointer footprint; REAL-SDF-HEADER-VIEW costs OPEN.
        let id = self.allocation.ok_or(SourceHeapError::NullPointer)?;
        Ok(SourceMutPointer {
            allocation: Some(id.subobject_field(self.element_offset, 3)?),
            element_offset: 0,
            marker: PhantomData,
        })
    }
}

impl SourceConstPointer<INPUT_PARMS> {
    pub(crate) fn sdf_header_field(self) -> Result<SourceMutPointer<i8>, SourceHeapError> {
        // SOURCE: ichiparm.c1266
        //                     ip->pSdfLabel = ip->szSdfDataHeader;
        // INCHI✔️❗: char65 is the first member of the actual parameter
        // record, including a parameter subobject in generator storage.
        // O(1), no allocation/copy; REAL-SDF-HEADER-VIEW costs OPEN.
        let id = self.allocation.ok_or(SourceHeapError::NullPointer)?;
        let allocation = if id.field >> 62 == 3 {
            if self.element_offset != 0 { return Err(SourceHeapError::PointerOutOfBounds); }
            id
        } else {
            id.subobject_field(self.element_offset, 3)?
        };
        Ok(SourceMutPointer { allocation: Some(allocation), element_offset: 0, marker: PhantomData })
    }
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
            let difference = i64::try_from(difference)
                .map_err(|_| SourceHeapError::PointerDifferenceOverflow)?;
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
            Some(id) => Some(id.logical()),
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
    slots: HashMap<usize, Option<AllocationSlot>>,
    issued_ids: u64,
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
    fn corresponding_integer_storage<T: 'static>(&self) -> bool {
        // VERBATIM SOURCE FRAGMENT: ichisize.h45..49 corresponding short types
        // typedef unsigned short AT_NUMB;
        // typedef unsigned short AT_RANK;
        // #define AT_RANK_MASK   ((AT_RANK)~0)
        // 
        // typedef signed short NUM_H;
        // INCHI✔️❗: C468 corresponding signed/unsigned primitive storage.
        // Same size/alignment/stride, all bits valid; no scalar conversion.
        // Exact-type callers bypass this O(1) check; whole costs remain OPEN.
        let requested = TypeId::of::<Vec<T>>();
        macro_rules! pair {
            ($signed:ty, $unsigned:ty) => {
                (self.type_id == TypeId::of::<Vec<$signed>>()
                    && requested == TypeId::of::<Vec<$unsigned>>())
                    || (self.type_id == TypeId::of::<Vec<$unsigned>>()
                        && requested == TypeId::of::<Vec<$signed>>())
            };
        }
        pair!(i8, u8) || pair!(i16, u16) || pair!(i32, u32)
            || pair!(i64, u64) || pair!(isize, usize)
    }

    #[inline(always)]
    fn downcast_ref<T: 'static>(&self) -> Option<&[T]> {
        if self.type_id != TypeId::of::<Vec<T>>()
            && !self.corresponding_integer_storage::<T>()
        {
            return None;
        }
        // SAFETY: metadata/raw parts retain the original owned Vec. The view
        // is either that exact element type or its recognized corresponding
        // primitive integer: identical size/alignment, every bit pattern valid.
        // Length/stride stay unchanged; the stored ORIGINAL drop type stays.
        Some(unsafe { std::slice::from_raw_parts(self.pointer.cast::<T>().as_ptr(), self.len) })
    }

    #[inline(always)]
    fn downcast_mut<T: 'static>(&mut self) -> Option<&mut [T]> {
        if self.type_id != TypeId::of::<Vec<T>>()
            && !self.corresponding_integer_storage::<T>()
        {
            return None;
        }
        // A mutable view may overwrite source pointers whose construction
        // established allocation relationships recorded below.
        self.provenance = AllocationProvenance::None;
        // SAFETY: the exact/corresponding primitive proof above applies;
        // `&mut self` guarantees unique ordinary access and writable stable
        // callers retain their existing distinct/live allocation obligations.
        Some(unsafe { std::slice::from_raw_parts_mut(self.pointer.cast::<T>().as_ptr(), self.len) })
    }

    fn input_parameters_field_ref(&self, row: usize) -> Result<&INPUT_PARMS, SourceHeapError> {
        // SOURCE: inchi_dll_a.c: INPUT_PARMS *ip = &( HGen->InpParms );
        // INCHI✔️❗: Resolve the live owner on every reached view; no clone.
        // O(1), REAL-SDF-HEADER-VIEW dispatch costs OPEN.
        if self.is::<INPUT_PARMS>() {
            self.downcast_ref::<INPUT_PARMS>()
                .ok_or(SourceHeapError::AllocationTypeMismatch)?
                .get(row).ok_or(SourceHeapError::PointerOutOfBounds)
        } else {
            Ok(&self.downcast_ref::<INCHIGEN_CONTROL>()
                .ok_or(SourceHeapError::AllocationTypeMismatch)?
                .get(row).ok_or(SourceHeapError::PointerOutOfBounds)?.InpParms)
        }
    }

    fn input_parameters_field_mut(&mut self, row: usize) -> Result<&mut INPUT_PARMS, SourceHeapError> {
        // SOURCE: inchi_dll_a.c: INPUT_PARMS *ip = &( HGen->InpParms );
        // INCHI✔️❗: Same actual field as input_parameters_field_ref, with exclusive
        // parent access. No registration/snapshot/unsafe; dispatch costs OPEN.
        if self.is::<INPUT_PARMS>() {
            self.downcast_mut::<INPUT_PARMS>()
                .ok_or(SourceHeapError::AllocationTypeMismatch)?
                .get_mut(row).ok_or(SourceHeapError::PointerOutOfBounds)
        } else {
            Ok(&mut self.downcast_mut::<INCHIGEN_CONTROL>()
                .ok_or(SourceHeapError::AllocationTypeMismatch)?
                .get_mut(row).ok_or(SourceHeapError::PointerOutOfBounds)?.InpParms)
        }
    }

    fn copy_disjoint_field_rows<T: Copy + 'static>(
        &mut self,
        destination_field: u64,
        source_field: u64,
        destination_offset: usize,
        source_offset: usize,
        count: usize,
    ) -> Result<(), SourceHeapError> {
        // SOURCE: util.c mystrncpy
        //     if (len)
        //     {
        //         memmove(target, source, len);
        //     }
        // INCHI✔️❗: Distinct admitted fields of these supported owner types
        // occupy distinct Vec rows. Split the actual rows once, then copy the
        // borrowed arrays; same-field overlap keeps copy_within below.
        // O(count), O(1) scratch, no allocation/clone/per-byte arena lookup.
        // Typed dispatch/proof constants stay REAL-SDF-HEADER-VIEW OPEN.
        let kind = destination_field >> 62;
        if kind == 0 || kind != source_field >> 62 {
            return Err(SourceHeapError::PointerAllocationMismatch);
        }
        let mask = (1_u64 << 62) - 1;
        let destination_row = usize::try_from(destination_field & mask)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let source_row = usize::try_from(source_field & mask)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        macro_rules! copy_rows {
            ($parent:ty, $($field:ident).+, $length:literal, $original:ty, $opposite:ty) => {{
                let rows = self.downcast_mut::<$parent>()
                    .ok_or(SourceHeapError::AllocationTypeMismatch)?;
                if destination_row == source_row {
                    return Err(SourceHeapError::PointerAllocationMismatch);
                }
                if destination_row.max(source_row) >= rows.len() {
                    return Err(SourceHeapError::PointerOutOfBounds);
                }
                let (destination, source) = if destination_row < source_row {
                    let (left, right) = rows.split_at_mut(source_row);
                    (&mut left[destination_row], &right[0])
                } else {
                    let (left, right) = rows.split_at_mut(destination_row);
                    (&mut right[0], &left[source_row])
                };
                let source: &dyn Any = if TypeId::of::<T>() == TypeId::of::<$opposite>() {
                    let array: &[$opposite; $length] = bytemuck::cast_slice::<$original, $opposite>(
                        &source.$($field).+,
                    ).try_into().map_err(|_| SourceHeapError::AllocationTypeMismatch)?;
                    array
                } else {
                    &source.$($field).+
                };
                let destination: &mut dyn Any = if TypeId::of::<T>() == TypeId::of::<$opposite>() {
                    let array: &mut [$opposite; $length] = bytemuck::cast_slice_mut::<$original, $opposite>(
                        &mut destination.$($field).+,
                    ).try_into().map_err(|_| SourceHeapError::AllocationTypeMismatch)?;
                    array
                } else {
                    &mut destination.$($field).+
                };
                let source = source.downcast_ref::<[T; $length]>()
                    .ok_or(SourceHeapError::AllocationTypeMismatch)?
                    .get(source_offset..).and_then(|values| values.get(..count))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let destination = destination.downcast_mut::<[T; $length]>()
                    .ok_or(SourceHeapError::AllocationTypeMismatch)?
                    .get_mut(destination_offset..).and_then(|values| values.get_mut(..count))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                destination.copy_from_slice(source);
                Ok(())
            }};
        }
        match kind {
            1 => copy_rows!(sp_ATOM, neighbor, 20, u16, i16),
            2 => copy_rows!(inp_ATOM, neighbor, 20, u16, i16),
            3 if self.is::<INPUT_PARMS>() => copy_rows!(INPUT_PARMS, szSdfDataHeader, 65, i8, u8),
            3 => copy_rows!(INCHIGEN_CONTROL, InpParms.szSdfDataHeader, 65, i8, u8),
            _ => Err(SourceHeapError::PointerAllocationMismatch),
        }
    }

    fn view_ref<T: 'static>(&self, field: u64) -> Result<&[T], SourceHeapError> {
        // VERBATIM SOURCE: ichicans.c564/2370, ichiprt1.c2635
        //             pCG->m_pNeighborsForSort = at[i].neighbor;
        //         pCG->m_pNeighborsForSort = at[i].neighbor;
        //             pCG->m_pn_RankForSort = at[j].neighbor;
        // INCHI✔️❗: Borrow the real field on each reached read, preserving
        // later owner writes and parent lifetime. No snapshot/registration.
        // O(1), no allocation; field/type dispatch costs remain C572 OPEN.
        if field == 0 {
            return self.downcast_ref::<T>().ok_or(SourceHeapError::AllocationTypeMismatch);
        }
        let row = usize::try_from(field & ((1_u64 << 62) - 1))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if field >> 62 == 3 {
            // SOURCE ichiparm.c1266: ip->pSdfLabel = ip->szSdfDataHeader;
            // INCHI✔️❗: The live record or its actual char65 first field,
            // selected by typed access; no copied storage (C580 costs OPEN).
            let parameters = self.input_parameters_field_ref(row)?;
            if TypeId::of::<T>() == TypeId::of::<INPUT_PARMS>() {
                return (parameters as &dyn Any).downcast_ref::<T>()
                    .map(std::slice::from_ref).ok_or(SourceHeapError::AllocationTypeMismatch);
            }
            let header: &dyn Any = if TypeId::of::<T>() == TypeId::of::<u8>() {
                let bytes: &[u8; 65] = bytemuck::cast_slice::<i8, u8>(&parameters.szSdfDataHeader)
                    .try_into().map_err(|_| SourceHeapError::AllocationTypeMismatch)?;
                bytes
            } else {
                &parameters.szSdfDataHeader
            };
            return header.downcast_ref::<[T; 65]>().map(|array| array.as_slice())
                .ok_or(SourceHeapError::AllocationTypeMismatch);
        }
        let neighbors = match field >> 62 {
            1 => &self.downcast_ref::<sp_ATOM>()
                .ok_or(SourceHeapError::AllocationTypeMismatch)?
                .get(row).ok_or(SourceHeapError::PointerOutOfBounds)?.neighbor,
            2 => &self.downcast_ref::<inp_ATOM>()
                .ok_or(SourceHeapError::AllocationTypeMismatch)?
                .get(row).ok_or(SourceHeapError::PointerOutOfBounds)?.neighbor,
            _ => return Err(SourceHeapError::AllocationTypeMismatch),
        };
        let array: &dyn Any = if TypeId::of::<T>() == TypeId::of::<i16>() {
            bytemuck::cast_ref::<[u16; 20], [i16; 20]>(neighbors)
        } else {
            neighbors
        };
        array.downcast_ref::<[T; 20]>().map(|array| array.as_slice())
            .ok_or(SourceHeapError::AllocationTypeMismatch)
    }

    fn view_mut<T: 'static>(&mut self, field: u64) -> Result<&mut [T], SourceHeapError> {
        // Same actual subobject as view_ref, with exclusive parent access.
        // Existing typed mutation invalidates parent provenance. No unsafe
        // conversion or independently owned field storage; C572 costs OPEN.
        if field == 0 {
            return self.downcast_mut::<T>().ok_or(SourceHeapError::AllocationTypeMismatch);
        }
        let row = usize::try_from(field & ((1_u64 << 62) - 1))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if field >> 62 == 3 {
            // SOURCE ichiparm.c1266: ip->pSdfLabel = ip->szSdfDataHeader;
            // INCHI✔️❗: Same live subobject and exclusive owner borrow;
            // no copied storage/unsafe conversion (C580 costs OPEN).
            let parameters = self.input_parameters_field_mut(row)?;
            if TypeId::of::<T>() == TypeId::of::<INPUT_PARMS>() {
                return (parameters as &mut dyn Any).downcast_mut::<T>()
                    .map(std::slice::from_mut).ok_or(SourceHeapError::AllocationTypeMismatch);
            }
            let header: &mut dyn Any = if TypeId::of::<T>() == TypeId::of::<u8>() {
                let bytes: &mut [u8; 65] = bytemuck::cast_slice_mut::<i8, u8>(&mut parameters.szSdfDataHeader)
                    .try_into().map_err(|_| SourceHeapError::AllocationTypeMismatch)?;
                bytes
            } else {
                &mut parameters.szSdfDataHeader
            };
            return header.downcast_mut::<[T; 65]>().map(|array| array.as_mut_slice())
                .ok_or(SourceHeapError::AllocationTypeMismatch);
        }
        let neighbors = match field >> 62 {
            1 => &mut self.downcast_mut::<sp_ATOM>()
                .ok_or(SourceHeapError::AllocationTypeMismatch)?
                .get_mut(row).ok_or(SourceHeapError::PointerOutOfBounds)?.neighbor,
            2 => &mut self.downcast_mut::<inp_ATOM>()
                .ok_or(SourceHeapError::AllocationTypeMismatch)?
                .get_mut(row).ok_or(SourceHeapError::PointerOutOfBounds)?.neighbor,
            _ => return Err(SourceHeapError::AllocationTypeMismatch),
        };
        let array: &mut dyn Any = if TypeId::of::<T>() == TypeId::of::<i16>() {
            bytemuck::cast_mut::<[u16; 20], [i16; 20]>(neighbors)
        } else {
            neighbors
        };
        array.downcast_mut::<[T; 20]>().map(|array| array.as_mut_slice())
            .ok_or(SourceHeapError::AllocationTypeMismatch)
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
        usize::try_from(id.logical()).map_err(|_| SourceHeapError::MissingAllocation)
    }

    fn insert<T: 'static>(
        &mut self,
        id: AllocationId,
        allocation: Vec<T>,
    ) -> Result<(), SourceHeapError> {
        let index = usize::try_from(id.logical()).map_err(|_| SourceHeapError::AllocationIdExhausted)?;
        if id.logical() != self.issued_ids || id.field != 0 {
            return Err(SourceHeapError::AllocationIdExhausted);
        }
        let issued_ids = self
            .issued_ids
            .checked_add(1)
            .ok_or(SourceHeapError::AllocationIdExhausted)?;
        // INCHI✔️❗: #define inchi_malloc   malloc
        // INCHI✔️❗: #define inchi_calloc   calloc
        // Native allocation macros: mode.h:1171..1175. Keep allocation IDs
        // monotonic, but retain descriptors only for live or detached buffers.
        // Expected/amortized O(1) registration; hashing/rehash costs remain
        // INCHI-PERF-ARENA OPEN, not native allocator performance equivalence.
        self.slots
            .try_reserve(1)
            .map_err(|_| SourceHeapError::AllocationFailed)?;
        self.slots.insert(index, Some(AllocationSlot::new(allocation)));
        self.issued_ids = issued_ids;
        self.live_count += 1;
        Ok(())
    }

    fn get(&self, id: AllocationId) -> Option<&AllocationSlot> {
        self.slots.get(&Self::index(id).ok()?)?.as_ref()
    }

    fn get_mut(&mut self, id: AllocationId) -> Option<&mut AllocationSlot> {
        self.slots.get_mut(&Self::index(id).ok()?)?.as_mut()
    }

    fn take(&mut self, id: AllocationId) -> Option<AllocationSlot> {
        let allocation = self.slots.get_mut(&Self::index(id).ok()?)?.take()?;
        self.live_count -= 1;
        Some(allocation)
    }

    fn restore(&mut self, id: AllocationId, allocation: AllocationSlot) {
        let index = Self::index(id).expect("an allocated ID fits usize");
        let slot = self.slots.get_mut(&index).expect("allocated slot exists");
        debug_assert!(slot.is_none());
        *slot = Some(allocation);
        self.live_count += 1;
    }

    fn remove(&mut self, id: AllocationId) -> Option<AllocationSlot> {
        self.take(id)
    }

    fn erase(&mut self, id: AllocationId) -> Option<AllocationSlot> {
        // INCHI✔️❗: #define inchi_free(X)  do{ if(X) free(X); }while(0)
        // mode.h:1180..1182: final release must not retain a historical slot.
        // Temporary remove/restore keeps its key, so restoration never allocates
        // and compaction cannot invalidate a detached allocation's reservation.
        let index = Self::index(id).ok()?;
        self.slots.get(&index)?.as_ref()?;
        let allocation = self.slots.remove(&index)?.take()?;
        self.live_count -= 1;
        let retained = self.slots.len();
        if retained == 0 {
            self.slots = HashMap::new();
        } else if self.slots.capacity() > 32 && retained < self.slots.capacity() / 4 {
            // Geometric compaction is amortized O(1) per release, O(retained)
            // on a compaction call. Native free cannot fail: reserve failure
            // leaves the existing capacity intact, not a chemistry error.
            // This failure-retained peak and hashing/temporary-table costs
            // remain explicit limitations in INCHI-PERF-ARENA.
            let mut compacted = HashMap::new();
            if compacted.try_reserve(retained).is_ok() {
                compacted.extend(std::mem::take(&mut self.slots));
                self.slots = compacted;
            }
        }
        Some(allocation)
    }

    fn len(&self) -> usize {
        self.live_count
    }

    #[cfg(test)]
    fn values(&self) -> impl Iterator<Item = &AllocationSlot> {
        self.slots.values().filter_map(Option::as_ref)
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

// A temporary loan retains its reserved arena key until the actual owner is
// returned. The reservation survives nested allocations/free compaction.
struct SourceAllocationLoan<'a> {
    heap: &'a mut SourceHeap,
    id: AllocationId,
    allocation: Option<AllocationSlot>,
}

impl Drop for SourceAllocationLoan<'_> {
    fn drop(&mut self) {
        // INCHI✔️❗: Restore the same actual owner on Result and unwinding.
        // O(1) expected arena return, no allocation/copy/catch or panic change;
        // C584 lifetime/check costs remain OPEN.
        if let Some(allocation) = self.allocation.take() {
            self.heap.allocations.restore(self.id, allocation);
        }
    }
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
    // Native ichinorm.c5920 (GCC/Linux LP64 sizeof(at)==8):
    //             memcpy(at + num_atoms, at_H, num_added_explicit_H * sizeof(at));
    // INCHI✔️❌: Preserve the source byte count, including partial fields.
    // Full rows copy their canonical values directly; only one partial final
    // row needs the existing LP64 field-byte projection. O(byte_count)/O(1),
    // no allocation; partial-row packing cost stays EXPLICIT-H-BYTE-COPY OPEN.
    let touched_atoms = byte_count.div_ceil(INP_ATOM_GCC_LP64_SIZE);
    if touched_atoms > destination.len() || touched_atoms > source.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let complete_atoms = byte_count / INP_ATOM_GCC_LP64_SIZE;
    destination[..complete_atoms].copy_from_slice(&source[..complete_atoms]);
    let partial_bytes = byte_count % INP_ATOM_GCC_LP64_SIZE;
    if partial_bytes != 0 {
        let source_bytes = inp_atom_gcc_lp64_bytes(&source[complete_atoms]);
        let mut destination_bytes = inp_atom_gcc_lp64_bytes(&destination[complete_atoms]);
        destination_bytes[..partial_bytes].copy_from_slice(&source_bytes[..partial_bytes]);
        overwrite_inp_atom_from_gcc_lp64_bytes(
            &mut destination[complete_atoms],
            &destination_bytes,
        );
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
        self.live_source_allocation_types
            .values()
            .copied()
            .collect()
    }

    #[cfg(test)]
    pub(crate) fn live_allocations_of<T: 'static>(&self) -> usize {
        self.allocations
            .values()
            .filter(|allocation| allocation.is::<T>())
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
        let next_id = self
            .next_id
            .checked_add(1)
            .ok_or(SourceHeapError::AllocationIdExhausted)?;
        let id = AllocationId {
            owner: NonZeroU64::new(next_id).ok_or(SourceHeapError::AllocationIdExhausted)?,
            field: 0,
        };
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
        let table_id = pointer_table
            .allocation
            .ok_or(SourceHeapError::NullPointer)?;
        let storage_id = storage.allocation.ok_or(SourceHeapError::NullPointer)?;
        if table_id.field != 0 || storage_id.field != 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        if table_id.owner == storage_id.owner {
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
        if table_id.field != 0 {
            return None;
        }
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
        if pointer.element_offset != 0 || pointer.allocation.is_some_and(|id| id.field != 0) {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let slot = self
            .allocations
            .get_mut(id)
            .ok_or(SourceHeapError::MissingAllocation)?;
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
        if pointer.element_offset != 0 || pointer.allocation.is_some_and(|id| id.field != 0) {
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
    pub(crate) fn slice<T: 'static>(
        &self,
        pointer: SourceConstPointer<T>,
    ) -> Result<&[T], SourceHeapError> {
        let (values, offset) = self.slice_with_origin(pointer)?;
        Ok(&values[offset..])
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn slice_with_origin<T: 'static>(
        &self,
        pointer: SourceConstPointer<T>,
    ) -> Result<(&[T], usize), SourceHeapError> {
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let values = self
            .allocations
            .get(id)
            .ok_or(SourceHeapError::MissingAllocation)?
            .view_ref::<T>(id.field)?;
        let offset = usize::try_from(pointer.element_offset)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️❗: Retain the actual allocation prefix for signed source
        // subscripts, e.g. strutil.c2409. No copy/allocation; C537 cost OPEN.
        if offset > values.len() { return Err(SourceHeapError::PointerOutOfBounds); }
        Ok((values, offset))
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn slice_mut<T: 'static>(
        &mut self,
        pointer: SourceMutPointer<T>,
    ) -> Result<&mut [T], SourceHeapError> {
        let (values, offset) = self.slice_mut_with_origin(pointer)?;
        Ok(&mut values[offset..])
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn slice_mut_with_origin<T: 'static>(
        &mut self,
        pointer: SourceMutPointer<T>,
    ) -> Result<(&mut [T], usize), SourceHeapError> {
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let values = self
            .allocations
            .get_mut(id)
            .ok_or(SourceHeapError::MissingAllocation)?
            .view_mut::<T>(id.field)?;
        let offset = usize::try_from(pointer.element_offset)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        // INCHI✔️❗: Retain the actual allocation prefix for signed source
        // subscripts, e.g. strutil.c2409. No copy/allocation; C537 cost OPEN.
        if offset > values.len() { return Err(SourceHeapError::PointerOutOfBounds); }
        Ok((values, offset))
    }

    /// Resolves a source signed subscript, including a negative offset from
    /// an interior pointer. Invalid negative accesses retain the prior error.
    #[track_caller]
    #[inline(always)]
    pub(crate) fn element_signed<T: 'static>(
        &self,
        pointer: SourceConstPointer<T>,
        index: i32,
    ) -> Result<&T, SourceHeapError> {
        // INCHI✔️❗: representative native signed subscript, ichicans.c:166:
        //     if (at[cur_atom].parity)
        // Positive admission matches slice/get; negative interior offsets
        // resolve within the same allocation. O(1), no allocation or clone;
        // checked heap/type/offset cost remains independently unresolved.
        if index >= 0 {
            return self
                .slice(pointer)?
                .get(index as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds);
        }
        let pointer = pointer
            .offset(i64::from(index))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        self.slice(pointer)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }

    #[track_caller]
    #[inline(always)]
    pub(crate) fn element_signed_mut<T: 'static>(
        &mut self,
        pointer: SourceMutPointer<T>,
        index: i32,
    ) -> Result<&mut T, SourceHeapError> {
        // INCHI✔️❗: representative native signed lvalue, ichicans.c:171:
        //     visited[cur_atom] = 1;
        // Preserve positive mutable admission/provenance and all prior error
        // categories; invalid negative accesses were PointerOutOfBounds.
        // O(1), no allocation/clone; whole checked-access cost unresolved.
        if index >= 0 {
            return self
                .slice_mut(pointer)?
                .get_mut(index as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds);
        }
        let pointer = pointer
            .offset(i64::from(index))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        self.slice_mut(pointer)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)
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
        let slot = self
            .allocations
            .get_mut(id)
            .ok_or(SourceHeapError::MissingAllocation)?;
        let values = slot
            .downcast_ref::<T>()
            .ok_or(SourceHeapError::AllocationTypeMismatch)?;
        let offset = usize::try_from(pointer.element_offset)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let values = values
            .get(offset..)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        Ok(StableSourceSlice {
            pointer: NonNull::new(values.as_ptr().cast_mut()).expect("slice pointers are non-null"),
            len: values.len(),
        })
    }

    /// Copies source-defined nonoverlapping typed prefixes without hiding a
    /// disjoint source range in the same modeled allocation.
    #[track_caller]
    pub(crate) fn copy_nonoverlapping_prefix<T: Clone + 'static>(
        &mut self,
        destination: SourceMutPointer<T>,
        source: SourceConstPointer<T>,
        count: usize,
        destination_bounds_first: bool,
    ) -> Result<(), SourceHeapError> {
        // INCHI✔️❗: representative verbatim source copy anchors:
        // ichicano.c:1343
        // memcpy(nCurrRank, nTempRank, num_at_tg * sizeof(nCurrRank[0]));
        // ichisort.c:989..990
        // memcpy(nNewAtomNumber, nPrevAtomNumber, num_atoms * sizeof(nNewAtomNumber[0]));
        // memcpy(nNewRank, nPrevRank, num_atoms * sizeof(nNewRank[0]));
        // O(count) typed copy, O(1) scratch, no allocation. Distinct allocations
        // retain extraction; same-allocation disjoint ranges use safe splitting.
        // Checked admission/dispatch cost remains INCHI-PERF-COPY-ALIAS OPEN.
        if destination.allocation_identity() != source.as_mut().allocation_identity() {
            return self.with_slice_mut_and_heap(destination, |target, heap| {
                if destination_bounds_first && count > target.len() {
                    return Err(SourceHeapError::PointerOutOfBounds);
                }
                let values = heap
                    .slice(source)?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                target
                    .get_mut(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone_from_slice(values);
                Ok(())
            });
        }

        // Admit the destination first, including mutable provenance invalidation,
        // exactly as extraction does. Callers retain their prefix-check order.
        let destination_len = self.slice_mut(destination)?.len();
        if destination_bounds_first && count > destination_len {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        self.slice(source)?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if count > destination_len {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        if count == 0 {
            return Ok(());
        }
        let difference = destination.difference(source.as_mut())?;
        let separation = usize::try_from(difference.unsigned_abs())
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if separation < count {
            // memcpy overlap is outside the native defined domain. Preserve
            // the previous modeled overlap error; existing self-copy gates stay.
            return Err(SourceHeapError::MissingAllocation);
        }
        if difference > 0 {
            let values = self.slice_mut(source.as_mut())?;
            let (source_values, destination_values) = values.split_at_mut(separation);
            destination_values[..count].clone_from_slice(&source_values[..count]);
        } else {
            let values = self.slice_mut(destination)?;
            let (destination_values, source_values) = values.split_at_mut(separation);
            destination_values[..count].clone_from_slice(&source_values[..count]);
        }
        Ok(())
    }

    /// Copies the two admitted prefixes of one allocation without a snapshot.
    #[track_caller]
    pub(crate) fn copy_prefix_within_allocation<T: Copy + 'static>(
        &mut self,
        destination: SourceMutPointer<T>,
        source: SourceConstPointer<T>,
        count: usize,
    ) -> Result<(), SourceHeapError> {
        // INCHI✔️❗: ichimak2.c:1269:
        // memcpy(pCanonRank, pCanonRankInv, num_at_tg * sizeof(pCanonRank[0]));
        // Source-defined disjoint ranges and prior safe modeled overlap both
        // copy in O(count)/O(1) scratch with no allocation. Source admission
        // precedes destination admission, exactly as the replaced snapshots.
        // Whole checked-arena/expanded-handle cost remains independently OPEN.
        self.slice(source)?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if destination.allocation != source.allocation {
            let (Some(destination_id), Some(source_id)) = (destination.allocation, source.allocation) else {
                return Err(SourceHeapError::PointerAllocationMismatch);
            };
            if destination_id.owner != source_id.owner
                || destination_id.field >> 62 != source_id.field >> 62
                || destination_id.field == 0
            {
                return Err(SourceHeapError::PointerAllocationMismatch);
            }
            // Different real subobjects of one parent can still be valid
            // memmove operands. Admit source first, destination before writes.
            self.slice_mut(destination)?.get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let source_offset = usize::try_from(source.element_offset)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let destination_offset = usize::try_from(destination.element_offset)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            return self.allocations.get_mut(destination_id)
                .ok_or(SourceHeapError::MissingAllocation)?
                .copy_disjoint_field_rows::<T>(destination_id.field, source_id.field,
                    destination_offset, source_offset, count);
        }
        let source_offset = usize::try_from(source.element_offset)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let base = SourceMutPointer {
            element_offset: 0,
            ..destination
        };
        let values = self.slice_mut(base)?;
        let destination_offset = usize::try_from(destination.element_offset)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        values
            .get(destination_offset..)
            .and_then(|values| values.get(..count))
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let source_end = source_offset
            .checked_add(count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        values.copy_within(source_offset..source_end, destination_offset);
        Ok(())
    }

    #[track_caller]
    pub(crate) fn with_slice_mut_and_heap<T: 'static, R>(
        &mut self,
        pointer: SourceMutPointer<T>,
        operation: impl FnOnce(&mut [T], &SourceHeap) -> Result<R, SourceHeapError>,
    ) -> Result<R, SourceHeapError> {
        let id = pointer.allocation.ok_or(SourceHeapError::NullPointer)?;
        let mut allocation = self
            .allocations
            .remove(id)
            .ok_or(SourceHeapError::MissingAllocation)?;
        let offset = match usize::try_from(pointer.element_offset) {
            Ok(offset) => offset,
            Err(_) => {
                self.allocations.restore(id, allocation);
                return Err(SourceHeapError::PointerOutOfBounds);
            }
        };
        let result = match allocation.view_mut::<T>(id.field) {
            Ok(values) => match values.get_mut(offset..) {
                Some(values) => operation(values, self),
                None => Err(SourceHeapError::PointerOutOfBounds),
            },
            Err(error) => Err(error),
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
        let allocation = self.allocations.remove(id)
            .ok_or(SourceHeapError::MissingAllocation)?;
        let mut loan = SourceAllocationLoan { heap: self, id, allocation: Some(allocation) };
        let offset = usize::try_from(pointer.element_offset)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let allocation = loan.allocation.as_mut().expect("loan owns its detached allocation");
        match allocation.view_mut::<T>(id.field) {
            Ok(values) => match values.get_mut(offset..) {
                Some(values) => operation(values, &mut *loan.heap),
                None => Err(SourceHeapError::PointerOutOfBounds),
            },
            Err(error) => Err(error),
        }
    }

    pub(crate) fn with_two_slices_mut_and_optional_third<A: 'static, B: 'static, C: 'static, R>(
        &mut self,
        first: SourceMutPointer<A>,
        second: SourceMutPointer<B>,
        third: Option<SourceMutPointer<C>>,
        operation: impl FnOnce(&mut [A], &mut [B], Option<&mut [C]>) -> Result<R, SourceHeapError>,
    ) -> Result<R, SourceHeapError> {
        self.with_slice_mut_and_optional_slices(first, Some(second), third, |first, second, third| {
            operation(first, second.ok_or(SourceHeapError::NullPointer)?, third)
        })
    }

    pub(crate) fn with_slice_mut_and_optional_slices<A: 'static, B: 'static, C: 'static, R>(
        &mut self,
        first: SourceMutPointer<A>,
        second: Option<SourceMutPointer<B>>,
        third: Option<SourceMutPointer<C>>,
        operation: impl FnOnce(&mut [A], Option<&mut [B]>, Option<&mut [C]>) -> Result<R, SourceHeapError>,
    ) -> Result<R, SourceHeapError> {
        // INCHI✔️❗: Same safe extraction/removal/restoration algorithm as the
        // original three-slot entry; an actual external fixed-array second
        // target needs no arena identity. Preserve optional-third NULL, all
        // presence-before-offset/type errors and pairwise alias checks.
        // No new backend/model/loan/unsafe/allocation; C542 check costs OPEN.
        let first_id = first.allocation.ok_or(SourceHeapError::NullPointer)?;
        let second_id = second.map(|p| p.allocation.ok_or(SourceHeapError::NullPointer)).transpose()?;
        let third_id = third.and_then(|p| p.allocation);
        if second_id.is_some_and(|id| id.owner == first_id.owner)
            || third_id.is_some_and(|id| id.owner == first_id.owner
                || second_id.is_some_and(|second| id.owner == second.owner))
        {
            return Err(SourceHeapError::PointerAllocationMismatch);
        }
        let mut first_allocation = self.allocations.remove(first_id)
            .ok_or(SourceHeapError::MissingAllocation)?;
        let mut second_allocation = if let Some(id) = second_id {
            match self.allocations.remove(id) {
                Some(allocation) => Some((id, allocation)),
                None => {
                    self.allocations.restore(first_id, first_allocation);
                    return Err(SourceHeapError::MissingAllocation);
                }
            }
        } else { None };
        let mut third_allocation = if let Some(id) = third_id {
            match self.allocations.remove(id) {
                Some(allocation) => Some((id, allocation)),
                None => {
                    self.allocations.restore(first_id, first_allocation);
                    if let Some((id, allocation)) = second_allocation {
                        self.allocations.restore(id, allocation);
                    }
                    return Err(SourceHeapError::MissingAllocation);
                }
            }
        } else { None };
        let result = (|| {
            let first_offset = usize::try_from(first.element_offset)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let second_offset = second.map(|p| usize::try_from(p.element_offset)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)).transpose()?;
            let first_values = first_allocation.view_mut::<A>(first_id.field)?
                .get_mut(first_offset..).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let second_values = match (&mut second_allocation, second_offset) {
                (Some((id, allocation)), Some(offset)) => Some(allocation.view_mut::<B>(id.field)?
                    .get_mut(offset..).ok_or(SourceHeapError::PointerOutOfBounds)?),
                _ => None,
            };
            let third_values = match (&mut third_allocation, third) {
                (Some((id, allocation)), Some(pointer)) => {
                    let offset = usize::try_from(pointer.element_offset)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    Some(allocation.view_mut::<C>(id.field)?
                        .get_mut(offset..).ok_or(SourceHeapError::PointerOutOfBounds)?)
                }
                _ => None,
            };
            operation(first_values, second_values, third_values)
        })();
        self.allocations.restore(first_id, first_allocation);
        if let Some((id, allocation)) = second_allocation { self.allocations.restore(id, allocation); }
        if let Some((id, allocation)) = third_allocation { self.allocations.restore(id, allocation); }
        result
    }


    pub(crate) fn free<T: 'static>(
        &mut self,
        pointer: SourceMutPointer<T>,
    ) -> Result<(), SourceHeapError> {
        let Some(id) = pointer.allocation else {
            return Ok(());
        };
        if pointer.element_offset != 0 || id.field != 0 {
            return Err(SourceHeapError::FreeOfInteriorPointer);
        }
        // VERBATIM SOURCE FRAGMENT: mode.h1179..1182 native untyped base free
        // #endif
        // #ifndef inchi_free
        // #define inchi_free(X)  do{ if(X) free(X); }while(0)
        // #endif
        // INCHI✔️❗: C467 access type does not constrain release of the base.
        // Erase drops via the ORIGINAL stored Vec destructor, never hint T.
        // Preserve missing/reserved errors without an extra get/type check;
        // no unsafe change, scan or clone. Arena/whole-call cost stays OPEN.
        self.allocations.erase(id)
            .ok_or(SourceHeapError::MissingAllocation)?;
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
