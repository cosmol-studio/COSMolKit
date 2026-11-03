//! Official basic InChI generation and reading, without toolkit processing.
//!
//! Inputs use the official InChI atom and stereo conventions. This boundary
//! neither kekulizes nor sanitizes a molecule. Existing graph-based APIs retain
//! their separate RDKit-compatible preprocessing and toolkit callbacks.

use std::ffi::CStr;

use crate::source::api::inchi_dll::{
    FreeINCHI, FreeStructFromINCHI, SourceInchiInputView, SourceStructTextInput,
    get_inchi_with_input, get_struct_from_inchi_input,
};
use crate::source::base::ichiparm::InchiBuildMetadata;
use crate::source_types::{
    FILE, SourceClock, SourceHeap, SourceHeapError, SourceMutPointer, inchi_Output,
    inchi_OutputStruct,
};
use crate::{InchiError, InchiErrorKind};

/// Official atom input record, including its fixed-capacity neighbor arrays.
pub use crate::source_types::inchi_Atom as Atom;
/// Official zero-dimensional stereo record (indices, type and parity).
pub use crate::source_types::inchi_Stereo0D as Stereo0D;

/// Borrowed basic GetINCHI input. Options use the official spelling unchanged.
#[derive(Clone, Copy, Debug)]
pub struct Input<'a> {
    pub atoms: &'a [Atom],
    pub stereo0d: &'a [Stereo0D],
    pub options: Option<&'a CStr>,
}

/// Official return status and independently optional output byte strings.
/// Missing fields remain absent; warnings and native errors retain their status.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Output {
    pub return_code: i32,
    pub inchi: Option<Vec<u8>>,
    pub message: Option<Vec<u8>>,
    pub log: Option<Vec<u8>>,
    pub aux_info: Option<Vec<u8>>,
}

/// Borrowed basic GetStructFromINCHI input, using official option spellings.
#[derive(Clone, Copy, Debug)]
pub struct ReadInput<'a> {
    pub inchi: &'a CStr,
    pub options: Option<&'a CStr>,
}

/// Owned basic GetStructFromINCHI output, without molecular-toolkit conversion.
///
/// Atom indices, hydrogens, isotopes and all zero-dimensional stereo records
/// retain the official conventions. Counts are the lengths of the two arrays.
/// Native errors remain observable in `return_code`, message, log and warnings;
/// they are not replaced by a toolkit error or a fallback molecule.
#[derive(Clone, Debug, PartialEq)]
pub struct ReadOutput {
    pub return_code: i32,
    pub atoms: Vec<Atom>,
    pub stereo0d: Vec<Stereo0D>,
    pub message: Option<Vec<u8>>,
    pub log: Option<Vec<u8>>,
    /// Official INCHIDIFF flags: first index selects reconnected/normal (0)
    /// or disconnected (1); second selects fixed-H (0) or main/mobile-H (1).
    pub warning_flags: [[u64; 2]; 2],
}

/// Read an InChI through the official basic engine, without a toolkit.
///
/// No sanitization, hydrogen removal or RDKit stereo conversion is performed.
/// Output owns its arrays and text; native allocations are released before
/// returning. `Err` denotes a Rust transport/engine failure, while native
/// InChI return codes (including invalid input) are retained in [`ReadOutput`].
/// Polymer/V3000 extension records belong to the separate upstream Ex API and
/// are not represented by this basic entry point.
pub fn read(input: ReadInput<'_>) -> Result<ReadOutput, InchiError> {
    // Official inchi_api.h:
    // INCHI✔️✔️: EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetStructFromINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct );
    // INCHI✔️✔️: EXPIMP_TEMPLATE INCHI_API void INCHI_DECL FreeStructFromINCHI( inchi_OutputStruct *out );
    // Ownership wrapper only: invoke the canonical source port once, copy each
    // returned array/text once, and free native buffers even if copying fails.
    // O(A + S + diagnostic bytes); no graph conversion or chemistry callback.
    let mut heap = SourceHeap::default();
    let options = input.options.map_or(c"", |options| options);
    let view = SourceStructTextInput::borrowed(
        bytemuck::cast_slice(input.inchi.to_bytes_with_nul()),
        bytemuck::cast_slice(options.to_bytes_with_nul()),
    );
    let mut output = inchi_OutputStruct::default();
    let result = (|| {
        let return_code = get_struct_from_inchi_input(
            &mut heap,
            Some(view),
            &mut output,
            SourceMutPointer::null(),
            InchiBuildMetadata {
                compiler: "gcc",
                date: "Jan 01 1970",
                time: "00:00:00",
            },
            crate::source_types::LiveSourceClock,
        )?;
        Ok(ReadOutput {
            return_code,
            atoms: copy_output_rows(&heap, output.atom, output.num_atoms)?,
            stereo0d: copy_output_rows(&heap, output.stereo0D, output.num_stereo0D)?,
            message: (!output.szMessage.is_null())
                .then(|| source_c_text(&heap, output.szMessage))
                .transpose()?,
            log: (!output.szLog.is_null())
                .then(|| source_c_text(&heap, output.szLog))
                .transpose()?,
            warning_flags: output.WarningFlags,
        })
    })();
    let cleanup = FreeStructFromINCHI(&mut heap, Some(&mut output));
    result
        .and_then(|result| cleanup.map(|()| result))
        .map_err(|error| crate::source_error("read", error))
}

fn copy_output_rows<T: Clone + 'static>(
    heap: &SourceHeap,
    pointer: SourceMutPointer<T>,
    count: i16,
) -> Result<Vec<T>, SourceHeapError> {
    let count = usize::try_from(count).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if count == 0 {
        return Ok(Vec::new());
    }
    let rows = heap
        .slice(pointer.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| SourceHeapError::AllocationFailed)?;
    result.extend_from_slice(rows);
    Ok(result)
}

/// Execute the official basic generator without calling any external toolkit.
/// Arrays are borrowed for the call; returned buffers have independent ownership.
pub fn generate(input: Input<'_>) -> Result<Output, InchiError> {
    if input.atoms.len() > i16::MAX as usize || input.stereo0d.len() > i16::MAX as usize {
        return Err(InchiError {
            operation: "generate",
            kind: InchiErrorKind::InvalidInput,
            detail: "atom and stereo counts must fit the official signed 16-bit input counts"
                .into(),
        });
    }
    let mut heap = SourceHeap::default();
    let view = SourceInchiInputView::borrowed_basic(
        input.atoms,
        input.stereo0d,
        input.options.map(CStr::to_bytes_with_nul),
    );
    let result: Result<Output, SourceHeapError> = (|| {
        let (copied, mut allocated) = generate_output(
            &mut heap,
            &view,
            SourceMutPointer::null(),
            InchiBuildMetadata {
                compiler: "gcc",
                date: "Jan 01 1970",
                time: "00:00:00",
            },
            crate::source_types::LiveSourceClock,
        )?;
        FreeINCHI(&mut heap, Some(&mut allocated))?;
        Ok(copied)
    })();
    result.map_err(|error| crate::generation_error("generate", error.into()))
}

pub(crate) fn source_c_text(
    heap: &SourceHeap,
    pointer: SourceMutPointer<i8>,
) -> Result<Vec<u8>, SourceHeapError> {
    let bytes = heap.slice(pointer.as_const())?;
    let end = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(end)
        .map_err(|_| SourceHeapError::AllocationFailed)?;
    result.extend_from_slice(bytemuck::cast_slice(&bytes[..end]));
    Ok(result)
}

pub(crate) fn generate_output(
    heap: &mut SourceHeap,
    input: &SourceInchiInputView<'_>,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
    clock: impl SourceClock,
) -> Result<(Output, inchi_Output), SourceHeapError> {
    // Official inchi_api.h:
    // INCHI✔️❗: EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetINCHI( inchi_Input *inp, inchi_Output *out );
    // INCHI✔️❗: EXPIMP_TEMPLATE INCHI_API void INCHI_DECL FreeINCHI( inchi_Output *out );
    // One canonical GetINCHI invocation; borrow input arrays and copy each
    // present output once before FreeINCHI. SourceHeap costs remain unchanged.
    let mut output = inchi_Output::default();
    let copied = (|| {
        let return_code = get_inchi_with_input(
            heap,
            Some(input),
            Some(&mut output),
            0,
            stdout,
            build,
            clock,
        )?;
        Ok(Output {
            return_code,
            inchi: (!output.szInChI.is_null())
                .then(|| source_c_text(heap, output.szInChI))
                .transpose()?,
            message: (!output.szMessage.is_null())
                .then(|| source_c_text(heap, output.szMessage))
                .transpose()?,
            log: (!output.szLog.is_null())
                .then(|| source_c_text(heap, output.szLog))
                .transpose()?,
            aux_info: (!output.szAuxInfo.is_null())
                .then(|| source_c_text(heap, output.szAuxInfo))
                .transpose()?,
        })
    })();
    match copied {
        Ok(copied) => Ok((copied, output)),
        Err(error) => {
            FreeINCHI(heap, Some(&mut output))?;
            Err(error)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn official_input_generates_methane_without_toolkit() {
        let mut atom = Atom::default();
        atom.elname[0] = b'C' as i8;
        atom.num_iso_H[0] = 4;
        let output = generate(Input {
            atoms: &[atom],
            stereo0d: &[],
            options: None,
        })
        .unwrap();
        assert_eq!(output.return_code, 0);
        assert_eq!(
            output.inchi.as_deref(),
            Some(b"InChI=1S/CH4/h1H4".as_slice())
        );
        assert_eq!(
            output.aux_info.as_deref(),
            Some(b"AuxInfo=1/0/N:1/rA:1C/rB:/rC:;".as_slice())
        );
    }

    #[test]
    fn official_input_rejects_truncated_counts() {
        let atoms = vec![Atom::default(); i16::MAX as usize + 1];
        let error = generate(Input {
            atoms: &atoms,
            stereo0d: &[],
            options: None,
        })
        .unwrap_err();
        assert_eq!(error.kind, InchiErrorKind::InvalidInput);
    }
}
