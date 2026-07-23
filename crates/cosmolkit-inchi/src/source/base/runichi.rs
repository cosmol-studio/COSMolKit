use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichiprt1::OrigStruct_FillOut;
use crate::source::base::mol_fmt4::OrigAtData_WriteToSDfile;
use crate::source_types::{
    _IS_ERROR, _IS_OKAY, CANON_GLOBALS, FLAG_INP_AT_CHIRAL, INCHI_IOS_STRING,
    INCHI_IOS_TYPE_STRING, INCHI_IOSTREAM, INCHI_OUT_SDFILE_ATOMS_DT, INCHI_OUT_SDFILE_ONLY,
    INPUT_PARMS, ORIG_ATOM_DATA, ORIG_STRUCT, STRUCT_DATA, SourceConstPointer, SourceHeap,
    SourceHeapError, SourceMutPointer,
};

fn runichi_c_text(
    heap: &SourceHeap,
    pointer: SourceConstPointer<i8>,
) -> Result<Vec<u8>, SourceHeapError> {
    if pointer.is_null() {
        return Ok(Vec::new());
    }
    let bytes = heap.slice(pointer)?;
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    Ok(bytes[..length].iter().map(|byte| *byte as u8).collect())
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_SaveMolfile(
    heap: &mut SourceHeap,
    original_input: &ORIG_ATOM_DATA,
    structure_data: &STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    input_number: i64,
    output: &mut INCHI_IOSTREAM,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:564 OrigAtData_SaveMolfile
    // INCHI✔️❌: int OrigAtData_SaveMolfile( ORIG_ATOM_DATA  *orig_inp_data,
    // INCHI✔️❌:                             STRUCT_DATA     *sd,
    // INCHI✔️❌:                             INPUT_PARMS     *ip,
    // INCHI✔️❌:                             long            num_inp,
    // INCHI✔️❌:                             INCHI_IOSTREAM  *out_file )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return _IS_OKAY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         char szNumber[256];
    // INCHI✔️❌:         sprintf(szNumber, "Structure #%ld. %s%s%s%s", num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue));
    // INCHI✔️❌:         ret = OrigAtData_WriteToSDfile( orig_inp_data, out_file, szNumber, NULL,
    // INCHI✔️❌:                                         ( sd->bChiralFlag&FLAG_INP_AT_CHIRAL ) ? 1 : 0,
    // INCHI✔️❌:                                         ( ip->bINChIOutputOptions&INCHI_OUT_SDFILE_ATOMS_DT ) ? 1 : 0,
    // INCHI✔️❌:                                         ip->pSdfLabel, ip->pSdfValue );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OrigAtData_SaveMolfile
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_SaveMolfile
    // INCHI✔️❌: #define SDF_LBL_VAL(L, V) ((L) && (L)[0]) ? gsSpace : gsEmpty, ((L) && (L)[0]) ? L : gsEmpty, ((L) && (L)[0]) ? (((V) && (V)[0]) ? gsEqual : gsSpace) : gsEmpty, ((V) && (V)[0]) ? V : ((L) && (L)[0]) ? gsMissing : gsEmpty
    // INCHI✔️❌: const char gsMissing[] = "is missing"; const char gsEmpty[] = "";
    // INCHI✔️❌: const char gsSpace[] = " "; const char gsEqual[] = "=";
    // INCHI✔️❌: #define INCHI_OUT_SDFILE_ONLY 0x0010
    // INCHI✔️❌: #define INCHI_OUT_SDFILE_ATOMS_DT 0x0800
    // INCHI✔️❌: #define FLAG_INP_AT_CHIRAL 1
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; sizeof(long) == 8
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_SaveMolfile

    if input_parameters.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 == 0 {
        return Ok(_IS_OKAY as i32);
    }

    let label = runichi_c_text(heap, input_parameters.pSdfLabel.as_const())?;
    let value = runichi_c_text(heap, input_parameters.pSdfValue.as_const())?;
    let mut title = format!("Structure #{input_number}. ").into_bytes();
    if !label.is_empty() {
        title.push(b' ');
        title.extend_from_slice(&label);
        if value.is_empty() {
            title.push(b' ');
            title.extend_from_slice(b"is missing");
        } else {
            title.push(b'=');
            title.extend_from_slice(&value);
        }
    } else if !value.is_empty() {
        title.extend_from_slice(&value);
    }
    if title.len() >= 256 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let title = heap.allocate_model_storage(
        title
            .iter()
            .copied()
            .chain(std::iter::once(0))
            .map(|byte| byte as i8)
            .collect(),
    )?;
    let result = OrigAtData_WriteToSDfile(
        heap,
        original_input,
        Some(output),
        SourceMutPointer::null(),
        title.as_const(),
        SourceConstPointer::null(),
        i32::from(structure_data.bChiralFlag & FLAG_INP_AT_CHIRAL as i32 != 0),
        i32::from(input_parameters.bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT as i32 != 0),
        input_parameters.pSdfLabel.as_const(),
        input_parameters.pSdfValue.as_const(),
    );
    heap.free(title)?;
    result
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OrigAtData_StoreNativeInput<'a>(
    heap: &mut SourceHeap,
    canon_globals: &mut CANON_GLOBALS,
    return_code: &mut i32,
    structure_data: &mut STRUCT_DATA,
    _input_parameters: &INPUT_PARMS,
    original_atom_data: &mut ORIG_ATOM_DATA,
    original_structure: &'a mut ORIG_STRUCT,
) -> Result<&'a mut ORIG_STRUCT, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:593 OrigAtData_StoreNativeInput
    // INCHI✔️✔️: ORIG_STRUCT * OrigAtData_StoreNativeInput( CANON_GLOBALS    *pCG,
    // INCHI✔️✔️:                                            int              *nRet,
    // INCHI✔️✔️:                                            STRUCT_DATA      *sd,
    // INCHI✔️✔️:                                            INPUT_PARMS      *ip,
    // INCHI✔️✔️:                                            ORIG_ATOM_DATA   *orig_inp_data,
    // INCHI✔️✔️:                                            ORIG_STRUCT      *pOrigStruct )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*    v. 1.05 always create and fill OrigStruc as it may be used to store e.g. polymer info    */
    // INCHI✔️✔️:     /*    If normal AuxInfo is requested, create full reversibility information from native inp data
    // INCHI✔️✔️:     if ( ip->bINChIOutputOptions & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO))
    // INCHI✔️✔️:         return NULL; */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (OrigStruct_FillOut( pCG, orig_inp_data, pOrigStruct, sd ))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         AddErrorMessage( sd->pStrErrStruct, "Cannot interpret reversibility information" );
    // INCHI✔️✔️:         sd->nStructReadError = 99;
    // INCHI✔️✔️:         sd->nErrorType = _IS_ERROR;
    // INCHI✔️✔️:         *nRet = _IS_ERROR;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return pOrigStruct;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OrigAtData_StoreNativeInput
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_StoreNativeInput
    // INCHI✔️✔️: #define _IS_ERROR 2
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_StoreNativeInput

    if OrigStruct_FillOut(
        heap,
        canon_globals,
        original_atom_data,
        original_structure,
        structure_data,
    )? != 0
    {
        let message = b"Cannot interpret reversibility information\0".map(|byte| byte as i8);
        AddErrorMessage(
            Some(&mut structure_data.pStrErrStruct),
            Some(&message),
        )?;
        structure_data.nStructReadError = 99;
        structure_data.nErrorType = _IS_ERROR as i32;
        *return_code = _IS_ERROR as i32;
    }

    Ok(original_structure)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::inp_ATOM;

    fn string_stream() -> INCHI_IOSTREAM {
        INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING::default(),
            ..INCHI_IOSTREAM::default()
        }
    }

    fn stream_bytes(heap: &SourceHeap, stream: &INCHI_IOSTREAM) -> Vec<u8> {
        if stream.s.pStr.is_null() {
            Vec::new()
        } else {
            heap.slice(stream.s.pStr.as_const()).unwrap()[..stream.s.nUsedLength as usize]
                .iter()
                .map(|byte| *byte as u8)
                .collect()
        }
    }

    fn c_text(heap: &mut SourceHeap, bytes: &[u8]) -> SourceMutPointer<i8> {
        heap.allocate_model_storage(
            bytes
                .iter()
                .copied()
                .chain(std::iter::once(0))
                .map(|byte| byte as i8)
                .collect(),
        )
        .unwrap()
    }

    #[test]
    fn source_port__runichi__origatdata_savemolfile__line_564() {
        let mut disabled_heap = SourceHeap::default();
        let disabled = INPUT_PARMS {
            pSdfLabel: SourceMutPointer::null(),
            pSdfValue: SourceMutPointer::null(),
            ..INPUT_PARMS::default()
        };
        let mut disabled_stream = string_stream();
        assert_eq!(
            OrigAtData_SaveMolfile(
                &mut disabled_heap,
                &ORIG_ATOM_DATA::default(),
                &STRUCT_DATA::default(),
                &disabled,
                i64::MIN,
                &mut disabled_stream,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(stream_bytes(&disabled_heap, &disabled_stream), b"");

        for (label_text, value_text, expected_title) in [
            (b"".as_slice(), b"".as_slice(), "Structure #-7. "),
            (b"".as_slice(), b"V".as_slice(), "Structure #-7. V"),
            (
                b"L".as_slice(),
                b"".as_slice(),
                "Structure #-7.  L is missing",
            ),
            (b"L".as_slice(), b"V".as_slice(), "Structure #-7.  L=V"),
        ] {
            let mut heap = SourceHeap::default();
            let label = c_text(&mut heap, label_text);
            let value = c_text(&mut heap, value_text);
            let parameters = INPUT_PARMS {
                pSdfLabel: label,
                pSdfValue: value,
                bINChIOutputOptions: INCHI_OUT_SDFILE_ONLY as i32,
                ..INPUT_PARMS::default()
            };
            let mut stream = string_stream();
            assert_eq!(
                OrigAtData_SaveMolfile(
                    &mut heap,
                    &ORIG_ATOM_DATA::default(),
                    &STRUCT_DATA::default(),
                    &parameters,
                    -7,
                    &mut stream,
                ),
                Ok(0)
            );
            let output = String::from_utf8(stream_bytes(&heap, &stream)).unwrap();
            assert_eq!(output.lines().next(), Some(expected_title));
            assert!(output.ends_with("$$$$\n"));
            if value_text.is_empty() {
                assert!(!output.contains("> <"));
            } else {
                let field = if label_text.is_empty() { "ID" } else { "L" };
                assert!(output.contains(&format!("> <{field}>\n V\n\n")));
            }
        }

        let mut heap = SourceHeap::default();
        let label = c_text(&mut heap, b"ISO");
        let value = c_text(&mut heap, b"D");
        let mut atom = inp_ATOM::default();
        atom.elname[..2].copy_from_slice(&[b'H' as i8, 0]);
        atom.el_number = 1;
        atom.iso_atw_diff = 2;
        atom.valence = 1;
        atom.chem_bonds_valence = 1;
        atom.bond_type[0] = 1;
        let atoms = heap.allocate_model_storage(vec![atom]).unwrap();
        let input = ORIG_ATOM_DATA {
            at: atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let parameters = INPUT_PARMS {
            pSdfLabel: label,
            pSdfValue: value,
            bINChIOutputOptions: (INCHI_OUT_SDFILE_ONLY | INCHI_OUT_SDFILE_ATOMS_DT) as i32,
            ..INPUT_PARMS::default()
        };
        let structure = STRUCT_DATA {
            bChiralFlag: FLAG_INP_AT_CHIRAL as i32,
            ..STRUCT_DATA::default()
        };
        let mut stream = string_stream();
        assert_eq!(
            OrigAtData_SaveMolfile(&mut heap, &input, &structure, &parameters, 9, &mut stream),
            Ok(0)
        );
        let output = String::from_utf8(stream_bytes(&heap, &stream)).unwrap();
        assert!(output.starts_with("Structure #9.  ISO=D\n"));
        assert!(output.contains("  1  0  0  0  1"));
        assert!(output.contains(" D   0  0"));
    }
}
