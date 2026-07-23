use crate::source::api::inchi_dll_b::add_source_error;
use crate::source::api::inchi_dll_ext::SetExtOrigAtDataByInChIExtInput;
use crate::source::base::ichimak2::WriteCoord;
use crate::source::base::mol2atom::FreeOrigAtData;
use crate::source::base::readinch::Extract0DParities;
use crate::source::base::runichi2::TreatErrorsInReadTheStructure;
use crate::source::base::util::{
    detect_unusual_el_valence, extract_charges_and_radicals, extract_h_atoms,
    get_atomic_mass_from_elnum, get_num_H, get_periodic_table_number, inchi_calloc, inchi_free,
    is_el_a_metal, is_in_the_list, mystrncpy_slice, n_bonds_val_to_metal,
};
use crate::source_types::{
    AB_PARITY_UNDF, AB_PARITY_UNKN, BOND_TYPE_ALTERN, BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE,
    BOND_TYPE_TRIPLE, FLAG_INP_AT_CHIRAL, FLAG_INP_AT_NONCHIRAL, FLAG_SET_INP_AT_CHIRAL,
    INCHI_IOSTREAM, INPUT_PARMS, ISOTOPIC_SHIFT_FLAG, ISOTOPIC_SHIFT_MAX, LOG_MASK_NO_WARN,
    MAX_ATOMS, MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BONDS, MAXVAL, MOL_COORD, NO_ATOM,
    NORMALLY_ALLOWED_INP_MAX_ATOMS, NUM_H_ISOTOPES, ORIG_ATOM_DATA, RADICAL_DOUBLET,
    RADICAL_TRIPLET, REQ_MODE_CHIR_FLG_STEREO, REQ_MODE_DIFF_UU_STEREO, REQ_MODE_RACEMIC_STEREO,
    REQ_MODE_RELATIVE_STEREO, REQ_MODE_STEREO, STEREO_DBLE_EITHER, STEREO_SNGL_DOWN,
    STEREO_SNGL_EITHER, STEREO_SNGL_UP, STRUCT_DATA, SourceArgvPointer, SourceHeap,
    SourceHeapError, SourceMutPointer, inchi_Atom, inchi_InputEx, inchi_Output, inchi_OutputStruct,
    inchi_Stereo0D, inp_ATOM, local_util::ERR_ELEM,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1DOWN,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1EITHER,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1UP,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2DOWN,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2EITHER,
    tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2UP, tagINCHIBondType_INCHI_BOND_TYPE_ALTERN,
    tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE, tagINCHIBondType_INCHI_BOND_TYPE_SINGLE,
    tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE, tagINCHIRadical_INCHI_RADICAL_DOUBLET,
    tagINCHIRadical_INCHI_RADICAL_NONE, tagINCHIRadical_INCHI_RADICAL_SINGLET,
    tagINCHIRadical_INCHI_RADICAL_TRIPLET, tagINCHIStereoParity0D_INCHI_PARITY_EVEN,
    tagINCHIStereoParity0D_INCHI_PARITY_ODD, tagINCHIStereoType0D_INCHI_StereoType_Allene,
    tagINCHIStereoType0D_INCHI_StereoType_DoubleBond,
    tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral,
};

const SOURCE_SIZEOF_INCHI_ATOM: u64 = 3 * 8 + 20 * 2 + 20 + 20 + 6 + 2 + 4 + 2 + 1 + 1;
const SOURCE_SIZEOF_INCHI_STEREO0D: u64 = 5 * 2 + 2;

#[allow(non_snake_case)]
pub(crate) fn SetNumImplicitH(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    number_of_atoms: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:994 SetNumImplicitH
    // INCHI✔❌: void SetNumImplicitH( inp_ATOM* at, int num_atoms )
    // INCHI✔❌: {
    // INCHI✔❌:     int bNonMetal;
    // INCHI✔❌:     int a1/*, n1*/;
    // INCHI✔❌:
    // INCHI✔❌:     /* special valences */
    // INCHI✔❌:     for (bNonMetal = 0; bNonMetal < 2; bNonMetal++)
    // INCHI✔❌:     {
    // INCHI✔❌:         for (a1 = 0; a1 < num_atoms; a1++)
    // INCHI✔❌:         {
    // INCHI✔❌:             int bHasMetalNeighbor /*, j*/;
    // INCHI✔❌:             if (bNonMetal != is_el_a_metal( at[a1].el_number ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 continue; /* first process all metals, after that all non-metals */
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             bHasMetalNeighbor = 0;
    // INCHI✔❌:             /***********************************************************************
    // INCHI✔❌:              *  Set number of hydrogen atoms
    // INCHI✔❌:              */
    // INCHI✔❌:             at[a1].num_H = get_num_H( at[a1].elname,
    // INCHI✔❌:                                       at[a1].num_H,
    // INCHI✔❌:                                       at[a1].num_iso_H,
    // INCHI✔❌:                                       at[a1].charge,
    // INCHI✔❌:                                       at[a1].radical,
    // INCHI✔❌:                                       at[a1].chem_bonds_valence,
    // INCHI✔❌:                                       0, /* instead of valence entered by the user: it does not exist here*/
    // INCHI✔❌:                                       ( at[a1].at_type & 1 )  /* bAliased */,
    // INCHI✔❌:                                       !( at[a1].at_type & 2 ) /* bDoNotAddH */,
    // INCHI✔❌:                                       bHasMetalNeighbor );
    // INCHI✔❌:             at[a1].at_type = 0;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌: }
    // END INCHI C FUNCTION: SetNumImplicitH

    for nonmetal in 0..2 {
        for atom_index in 0..number_of_atoms {
            let atom = heap.slice(atoms.as_const())?[atom_index as usize].clone();
            if nonmetal != is_el_a_metal(atom.el_number as i32)? {
                continue;
            }
            let hydrogen_count = get_num_H(
                Some(&atom.elname),
                atom.num_H as i32,
                Some(&atom.num_iso_H),
                atom.charge as i32,
                atom.radical as i32,
                atom.chem_bonds_valence as i32,
                0,
                (atom.at_type & 1) as i32,
                i32::from(atom.at_type & 2 == 0),
                0,
            )?;
            let output = &mut heap.slice_mut(atoms)?[atom_index as usize];
            output.num_H = hydrogen_count as i8;
            output.at_type = 0;
        }
    }
    Ok(())
}

pub(crate) fn FreeINCHI(
    heap: &mut SourceHeap,
    output: Option<&mut inchi_Output>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:153 FreeINCHI
    // INCHI✔❌: void INCHI_DECL FreeINCHI( inchi_Output *out )
    // INCHI✔❌: {
    // INCHI✔❌:     if (!out)
    // INCHI✔❌:     {
    // INCHI✔❌:         return;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (out->szInChI)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_free( out->szInChI );
    // INCHI✔❌:     }
    // INCHI✔❌:     if (out->szLog)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_free( out->szLog );
    // INCHI✔❌:     }
    // INCHI✔❌:     if (out->szMessage)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_free( out->szMessage );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌: }
    // END INCHI C FUNCTION: FreeINCHI

    let Some(output) = output else {
        return Ok(());
    };
    if !output.szInChI.is_null() {
        inchi_free(heap, output.szInChI)?;
    }
    if !output.szLog.is_null() {
        inchi_free(heap, output.szLog)?;
    }
    if !output.szMessage.is_null() {
        inchi_free(heap, output.szMessage)?;
    }
    *output = inchi_Output::default();
    Ok(())
}

pub(crate) fn FreeStdINCHI(
    heap: &mut SourceHeap,
    output: Option<&mut inchi_Output>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:183 FreeStdINCHI
    // INCHI✔❌: void INCHI_DECL FreeStdINCHI( inchi_Output *out )
    // INCHI✔❌: {
    // INCHI✔❌:     FreeINCHI( out );
    // INCHI✔❌: }
    // END INCHI C FUNCTION: FreeStdINCHI

    FreeINCHI(heap, output)
}

#[allow(non_snake_case)]
pub(crate) fn FreeStructFromINCHI(
    heap: &mut SourceHeap,
    output: Option<&mut inchi_OutputStruct>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:208 FreeStructFromINCHI
    // INCHI✔️❌: void INCHI_DECL FreeStructFromINCHI( inchi_OutputStruct *out )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (!out)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (out->atom)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( out->atom );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (out->stereo0D)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( out->stereo0D );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (out->szLog)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( out->szLog );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (out->szMessage)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( out->szMessage );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FreeStructFromINCHI

    let Some(output) = output else {
        return Ok(());
    };
    if !output.atom.is_null() {
        inchi_free(heap, output.atom)?;
    }
    if !output.stereo0D.is_null() {
        inchi_free(heap, output.stereo0D)?;
    }
    if !output.szLog.is_null() {
        inchi_free(heap, output.szLog)?;
    }
    if !output.szMessage.is_null() {
        inchi_free(heap, output.szMessage)?;
    }
    *output = inchi_OutputStruct::default();
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn FreeStructFromStdINCHI(
    heap: &mut SourceHeap,
    output: Option<&mut inchi_OutputStruct>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:195 FreeStructFromStdINCHI
    // INCHI✔️❌: void INCHI_DECL FreeStructFromStdINCHI( inchi_OutputStruct *out )
    // INCHI✔️❌: {
    // INCHI✔️❌:     FreeStructFromINCHI( out );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FreeStructFromStdINCHI

    FreeStructFromINCHI(heap, output)
}

#[allow(non_snake_case)]
pub(crate) fn GetStringLength(
    heap: &SourceHeap,
    pointer: SourceMutPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2086 GetStringLength
    // INCHI✔️❌: int INCHI_DECL GetStringLength( char *p )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return (int) strlen( p );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetStringLength

    if pointer.is_null() {
        return Ok(0);
    }
    let length = heap
        .slice(pointer.as_const())?
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    i32::try_from(length).map_err(|_| SourceHeapError::SourceIntegerOverflow)
}

pub(crate) fn parse_options_string(
    heap: &mut SourceHeap,
    command_pointer: SourceMutPointer<i8>,
    arguments: &mut [SourceArgvPointer],
    max_arguments: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1037 parse_options_string
    // INCHI✔❌: int parse_options_string( char *cmd, const char *argv[], int maxargs )
    // INCHI✔❌: {
    // INCHI✔❌:     char    *p;
    // INCHI✔❌:     char    *pArgCurChar;
    // INCHI✔❌:     int      bInsideQuotes;
    // INCHI✔❌:     int      bCopyCharToArg;
    // INCHI✔❌:     int      nNumBackSlashes;
    // INCHI✔❌:     int      i;
    // INCHI✔❌:
    // INCHI✔❌:     i = 0;
    // INCHI✔❌:     argv[i++] = ""; /* zeroth argument is not used */
    // INCHI✔❌:     p = cmd;
    // INCHI✔❌:     bInsideQuotes = 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* arguments, one by one */
    // INCHI✔❌:     while (i < maxargs - 1)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* bypass spaces */
    // INCHI✔❌:         while (*p == ' ' || *p == '\t')
    // INCHI✔❌:         {
    // INCHI✔❌:             p++;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (!*p)
    // INCHI✔❌:         {
    // INCHI✔❌:             break;
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         /* scan an argument */
    // INCHI✔❌:         argv[i++] = pArgCurChar = p;     /* store preliminary ptr to arg */
    // INCHI✔❌:
    // INCHI✔❌:         while (1)
    // INCHI✔❌:         {
    // INCHI✔❌:             bCopyCharToArg = 1;
    // INCHI✔❌:             nNumBackSlashes = 0;
    // INCHI✔❌:             while (*p == '\\')
    // INCHI✔❌:             {
    // INCHI✔❌:                 ++p;
    // INCHI✔❌:                 ++nNumBackSlashes;
    // INCHI✔❌:             }
    // INCHI✔❌:
    // INCHI✔❌:             /* each pair of backslashes => one backslash; one more backslash => literal quote */
    // INCHI✔❌:             if (*p == '\"')
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* one " found */
    // INCHI✔❌:                 if (nNumBackSlashes % 2 == 0)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (bInsideQuotes)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (*( p + 1 ) == '\"')
    // INCHI✔❌:                         {
    // INCHI✔❌:                             p++;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             bCopyCharToArg = 0;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else
    // INCHI✔❌:                     {
    // INCHI✔❌:                         bCopyCharToArg = 0;
    // INCHI✔❌:                     }
    // INCHI✔❌:                     bInsideQuotes = !bInsideQuotes;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 nNumBackSlashes /= 2;          /* divide nNumBackSlashes by two */
    // INCHI✔❌:             }
    // INCHI✔❌:             while (nNumBackSlashes--)
    // INCHI✔❌:             {
    // INCHI✔❌:                 *pArgCurChar++ = '\\';
    // INCHI✔❌:             }
    // INCHI✔❌:             if (!*p)
    // INCHI✔❌:             {
    // INCHI✔❌:                 break;
    // INCHI✔❌:             }
    // INCHI✔❌:             if (!bInsideQuotes && ( *p == ' ' || *p == '\t' ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 p++;
    // INCHI✔❌:                 /* move to the next char because this char may become
    // INCHI✔❌:                  * zero due to  *pArgCurChar++ = '\0'; line below */
    // INCHI✔❌:                 break;
    // INCHI✔❌:             }
    // INCHI✔❌:             if (bCopyCharToArg)
    // INCHI✔❌:             {
    // INCHI✔❌:                 *pArgCurChar++ = *p;
    // INCHI✔❌:             }
    // INCHI✔❌:             ++p;
    // INCHI✔❌:         }
    // INCHI✔❌:         *pArgCurChar++ = '\0';  /* argument zero termination */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* The last argument is NULL */
    // INCHI✔❌:     argv[i] = NULL;
    // INCHI✔❌:
    // INCHI✔❌:     return i;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: parse_options_string

    let max_arguments =
        usize::try_from(max_arguments).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if max_arguments < 2 || arguments.len() < max_arguments {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let command = heap.slice_mut(command_pointer)?;
    if !command.contains(&0) {
        return Err(SourceHeapError::MissingNulTerminator);
    }

    let mut argument_index = 0_usize;
    arguments[argument_index] = SourceArgvPointer::EmptyLiteral;
    argument_index += 1;
    let mut input_index = 0_usize;
    let mut inside_quotes = false;

    while argument_index < max_arguments - 1 {
        while matches!(command.get(input_index), Some(byte) if *byte == b' ' as i8 || *byte == b'\t' as i8)
        {
            input_index += 1;
        }
        let input = *command
            .get(input_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if input == 0 {
            break;
        }

        let argument_start = input_index;
        arguments[argument_index] = SourceArgvPointer::Command(command_pointer.offset(
            i64::try_from(argument_start).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
        )?);
        argument_index += 1;
        let mut output_index = input_index;

        loop {
            let mut copy_character = true;
            let mut backslash_count = 0_usize;
            while command.get(input_index) == Some(&(b'\\' as i8)) {
                input_index += 1;
                backslash_count += 1;
            }

            if command.get(input_index) == Some(&(b'"' as i8)) {
                if backslash_count % 2 == 0 {
                    if inside_quotes {
                        if command.get(input_index + 1) == Some(&(b'"' as i8)) {
                            input_index += 1;
                        } else {
                            copy_character = false;
                        }
                    } else {
                        copy_character = false;
                    }
                    inside_quotes = !inside_quotes;
                }
                backslash_count /= 2;
            }
            for _ in 0..backslash_count {
                command[output_index] = b'\\' as i8;
                output_index += 1;
            }
            let input = *command
                .get(input_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if input == 0 {
                break;
            }
            if !inside_quotes && (input == b' ' as i8 || input == b'\t' as i8) {
                input_index += 1;
                break;
            }
            if copy_character {
                command[output_index] = input;
                output_index += 1;
            }
            input_index += 1;
        }
        command[output_index] = 0;
    }

    arguments[argument_index] = SourceArgvPointer::Null;
    i32::try_from(argument_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn SetAtomProperties(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    coordinates: SourceMutPointer<MOL_COORD>,
    input_atoms: SourceMutPointer<inchi_Atom>,
    atom_index: i32,
    dimensions: Option<&mut i32>,
    error_text: Option<SourceMutPointer<i8>>,
    error: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1139 SetAtomProperties
    // INCHI✔❌: int SetAtomProperties( inp_ATOM *at,
    // INCHI✔❌:                        MOL_COORD *szCoord,
    // INCHI✔❌:                        inchi_Atom *ati,
    // INCHI✔❌:                        int a1,
    // INCHI✔❌:                        int *nDim,
    // INCHI✔❌:                        char *pStrErr,
    // INCHI✔❌:                        int *err )
    // INCHI✔❌: {
    // INCHI✔❌:     S_CHAR      cRadical;
    // INCHI✔❌:
    // INCHI✔❌:     /* element, check later */
    // INCHI✔❌:     strcpy( at[a1].elname, ati[a1].elname );
    // INCHI✔❌:
    // INCHI✔❌:     /* charge */
    // INCHI✔❌:     at[a1].charge = ati[a1].charge;
    // INCHI✔❌:
    // INCHI✔❌:     /* radical */
    // INCHI✔❌:     switch (ati[a1].radical)
    // INCHI✔❌:     {
    // INCHI✔❌:         case   INCHI_RADICAL_NONE:
    // INCHI✔❌:             cRadical = 0;
    // INCHI✔❌:             break;
    // INCHI✔❌:         case   INCHI_RADICAL_SINGLET:
    // INCHI✔❌: #if( SINGLET_IS_TRIPLET == 1) /* 'singlet' means two electrons make a lone pair instead of 2 bonds*/
    // INCHI✔❌:                               /* its effect on valence is same as the effect of a triplet */
    // INCHI✔❌:             cRadical = RADICAL_TRIPLET;
    // INCHI✔❌: #else
    // INCHI✔❌:             cRadical = RADICAL_SINGLET;
    // INCHI✔❌: #endif
    // INCHI✔❌:             break;
    // INCHI✔❌:         case   INCHI_RADICAL_DOUBLET:
    // INCHI✔❌:             cRadical = RADICAL_DOUBLET;
    // INCHI✔❌:             break;
    // INCHI✔❌:         case   INCHI_RADICAL_TRIPLET:
    // INCHI✔❌:             cRadical = RADICAL_TRIPLET;
    // INCHI✔❌:             break;
    // INCHI✔❌:         default:
    // INCHI✔❌:         {
    // INCHI✔❌:             char szRadicalType[16];
    // INCHI✔❌:             int nRad = ati[a1].radical;
    // INCHI✔❌:             while (nRad > RADICAL_TRIPLET)
    // INCHI✔❌:             {
    // INCHI✔❌:                 nRad -= 2;
    // INCHI✔❌:             }
    // INCHI✔❌:             sprintf( szRadicalType, "%d->%d", ati[a1].radical, nRad );
    // INCHI✔❌:             TREAT_ERR( *err, 0, "Radical center type replaced:" );
    // INCHI✔❌:             TREAT_ERR( *err, 0, szRadicalType );
    // INCHI✔❌:             cRadical = nRad;
    // INCHI✔❌:             if (nRad < 0)
    // INCHI✔❌:             {
    // INCHI✔❌:                 *err |= 8; /*  Unrecognized Radical replaced with non-radical */
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         break;
    // INCHI✔❌:     }
    // INCHI✔❌:     at[a1].radical = cRadical;
    // INCHI✔❌:
    // INCHI✔❌:     /* coordinates */
    // INCHI✔❌:     at[a1].x = ati[a1].x;
    // INCHI✔❌:     at[a1].y = ati[a1].y;
    // INCHI✔❌:     at[a1].z = ati[a1].z;
    // INCHI✔❌:
    // INCHI✔❌:     if (szCoord)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* store text coordinates */
    // INCHI✔❌:         char str[32];
    // INCHI✔❌:         MOL_COORD * coord_p = szCoord + a1;
    // INCHI✔❌:         WriteCoord( str, ati[a1].x );
    // INCHI✔❌:         memcpy( *coord_p, str, 10 );
    // INCHI✔❌:         WriteCoord( str, ati[a1].y );
    // INCHI✔❌:         memcpy( *coord_p + 10, str, 10 );
    // INCHI✔❌:         WriteCoord( str, ati[a1].z );
    // INCHI✔❌:         memcpy( *coord_p + 20, str, 10 );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (MIN_BOND_LENGTH < fabs( ati[a1].x ) || MIN_BOND_LENGTH < fabs( ati[a1].y ) || MIN_BOND_LENGTH < fabs( ati[a1].z ))
    // INCHI✔❌:     {
    // INCHI✔❌:         if (MIN_BOND_LENGTH < fabs( ati[a1].z ))
    // INCHI✔❌:         {
    // INCHI✔❌:             *nDim |= 3;
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             *nDim |= 2;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* orig. at. number */
    // INCHI✔❌:     at[a1].orig_at_number = a1 + 1;
    // INCHI✔❌:     return 0;
    // INCHI✔❌:
    // INCHI✔❌: #undef MIN_BOND_LENGTH
    // INCHI✔❌: }
    // END INCHI C FUNCTION: SetAtomProperties

    let atom_index_usize =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let input = heap
        .slice(input_atoms.as_const())?
        .get(atom_index_usize)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let element_length = input
        .elname
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;

    let radical = match i32::from(input.radical) {
        value if value == tagINCHIRadical_INCHI_RADICAL_NONE as i32 => 0,
        value if value == tagINCHIRadical_INCHI_RADICAL_SINGLET as i32 => RADICAL_TRIPLET as i32,
        value if value == tagINCHIRadical_INCHI_RADICAL_DOUBLET as i32 => RADICAL_DOUBLET as i32,
        value if value == tagINCHIRadical_INCHI_RADICAL_TRIPLET as i32 => RADICAL_TRIPLET as i32,
        value => {
            let mut replacement = value;
            while replacement > RADICAL_TRIPLET as i32 {
                replacement -= 2;
            }
            let error = error.ok_or(SourceHeapError::NullPointer)?;
            add_source_error(heap, error_text, "Radical center type replaced:")?;
            add_source_error(heap, error_text, &format!("{value}->{replacement}"))?;
            if replacement < 0 {
                *error |= 8;
            }
            replacement
        }
    };

    {
        let target = heap
            .slice_mut(atoms)?
            .get_mut(atom_index_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        target.elname[..=element_length].copy_from_slice(&input.elname[..=element_length]);
        target.charge = input.charge;
        target.radical = radical as i8;
        target.x = input.x;
        target.y = input.y;
        target.z = input.z;
        target.orig_at_number = u16::try_from(
            atom_index
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?,
        )
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    }

    if !coordinates.is_null() {
        let coordinate = heap
            .slice_mut(coordinates)?
            .get_mut(atom_index_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for (offset, value) in [(0, input.x), (10, input.y), (20, input.z)] {
            let mut text = [0_i8; 32];
            WriteCoord(&mut text, value)?;
            coordinate[offset..offset + 10].copy_from_slice(&text[..10]);
        }
    }

    const MIN_BOND_LENGTH: f64 = 1.0e-6;
    if MIN_BOND_LENGTH < input.x.abs()
        || MIN_BOND_LENGTH < input.y.abs()
        || MIN_BOND_LENGTH < input.z.abs()
    {
        let dimensions = dimensions.ok_or(SourceHeapError::NullPointer)?;
        if MIN_BOND_LENGTH < input.z.abs() {
            *dimensions |= 3;
        } else {
            *dimensions |= 2;
        }
    }
    Ok(0)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn SetBondProperties(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    input_atoms: SourceMutPointer<inchi_Atom>,
    atom_index: i32,
    bond_index: i32,
    num_atoms: i32,
    num_bonds: Option<&mut i32>,
    error_text: Option<SourceMutPointer<i8>>,
    mut error: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1235 SetBondProperties
    // INCHI✔❌: int SetBondProperties( inp_ATOM *at,
    // INCHI✔❌:                        inchi_Atom *ati,
    // INCHI✔❌:                        int a1,
    // INCHI✔❌:                        int j,
    // INCHI✔❌:                        int nNumAtoms,
    // INCHI✔❌:                        int *nNumBonds,
    // INCHI✔❌:                        char *pStrErr,
    // INCHI✔❌:                        int *err )
    // INCHI✔❌: {
    // INCHI✔❌:     int a2;
    // INCHI✔❌:     S_CHAR     cBondType, cStereoType1, cStereoType2;
    // INCHI✔❌:     AT_NUMB   *p1, *p2;
    // INCHI✔❌:     int        n1, n2;
    // INCHI✔❌:
    // INCHI✔❌:     /* bond type */
    // INCHI✔❌:     switch (ati[a1].bond_type[j])
    // INCHI✔❌:     {
    // INCHI✔❌:         case INCHI_BOND_TYPE_SINGLE:
    // INCHI✔❌:             cBondType = BOND_TYPE_SINGLE;
    // INCHI✔❌:             break;
    // INCHI✔❌:         case INCHI_BOND_TYPE_DOUBLE:
    // INCHI✔❌:             cBondType = BOND_TYPE_DOUBLE;
    // INCHI✔❌:             break;
    // INCHI✔❌:         case INCHI_BOND_TYPE_TRIPLE:
    // INCHI✔❌:             cBondType = BOND_TYPE_TRIPLE;
    // INCHI✔❌:             break;
    // INCHI✔❌:         case INCHI_BOND_TYPE_ALTERN:
    // INCHI✔❌:             cBondType = BOND_TYPE_ALTERN;
    // INCHI✔❌:             break;
    // INCHI✔❌:         default:
    // INCHI✔❌:         {
    // INCHI✔❌:             char szBondType[16];
    // INCHI✔❌:             sprintf( szBondType, "%d", ati[a1].bond_type[j] );
    // INCHI✔❌:             TREAT_ERR( *err, 0, "Unrecognized bond type:" );
    // INCHI✔❌:             TREAT_ERR( *err, 0, szBondType );
    // INCHI✔❌:             *err |= 8; /*  Unrecognized Bond type replaced with single bond */
    // INCHI✔❌:             cBondType = BOND_TYPE_SINGLE;
    // INCHI✔❌:         }
    // INCHI✔❌:         break;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* 2D stereo */
    // INCHI✔❌:
    // INCHI✔❌:     switch (ati[a1].bond_stereo[j])
    // INCHI✔❌:     {
    // INCHI✔❌:     /* stereocenter-related; positive: the sharp end points to this atom  */
    // INCHI✔❌:         case   INCHI_BOND_STEREO_NONE:
    // INCHI✔❌:             cStereoType1 = 0;
    // INCHI✔❌:             cStereoType2 = 0;
    // INCHI✔❌:             break;
    // INCHI✔❌:         case   INCHI_BOND_STEREO_SINGLE_1UP:
    // INCHI✔❌:             cStereoType1 = STEREO_SNGL_UP;
    // INCHI✔❌:             cStereoType2 = -STEREO_SNGL_UP;
    // INCHI✔❌:             break;
    // INCHI✔❌:         case   INCHI_BOND_STEREO_SINGLE_1EITHER:
    // INCHI✔❌:             cStereoType1 = STEREO_SNGL_EITHER;
    // INCHI✔❌:             cStereoType2 = -STEREO_SNGL_EITHER;
    // INCHI✔❌:             break;
    // INCHI✔❌:         case   INCHI_BOND_STEREO_SINGLE_1DOWN:
    // INCHI✔❌:             cStereoType1 = STEREO_SNGL_DOWN;
    // INCHI✔❌:             cStereoType2 = -STEREO_SNGL_DOWN;
    // INCHI✔❌:             break;
    // INCHI✔❌:         /* stereocenter-related; negative: the sharp end points to the opposite atom  */
    // INCHI✔❌:         case   INCHI_BOND_STEREO_SINGLE_2UP:
    // INCHI✔❌:             cStereoType1 = -STEREO_SNGL_UP;
    // INCHI✔❌:             cStereoType2 = STEREO_SNGL_UP;
    // INCHI✔❌:             break;
    // INCHI✔❌:         case   INCHI_BOND_STEREO_SINGLE_2EITHER:
    // INCHI✔❌:             cStereoType1 = -STEREO_SNGL_EITHER;
    // INCHI✔❌:             cStereoType2 = STEREO_SNGL_EITHER;
    // INCHI✔❌:             break;
    // INCHI✔❌:         case   INCHI_BOND_STEREO_SINGLE_2DOWN:
    // INCHI✔❌:             cStereoType1 = -STEREO_SNGL_DOWN;
    // INCHI✔❌:             cStereoType2 = STEREO_SNGL_DOWN;
    // INCHI✔❌:             break;
    // INCHI✔❌:         /* stereobond-related */
    // INCHI✔❌:         case   INCHI_BOND_STEREO_DOUBLE_EITHER:
    // INCHI✔❌:         case  -INCHI_BOND_STEREO_DOUBLE_EITHER:
    // INCHI✔❌:             cStereoType1 = STEREO_DBLE_EITHER;
    // INCHI✔❌:             cStereoType2 = STEREO_DBLE_EITHER;
    // INCHI✔❌:             break;
    // INCHI✔❌:         default:
    // INCHI✔❌:         {
    // INCHI✔❌:             char szBondType[16];
    // INCHI✔❌:             sprintf( szBondType, "%d", ati[a1].bond_stereo[j] );
    // INCHI✔❌:             TREAT_ERR( *err, 0, "Unrecognized bond stereo:" );
    // INCHI✔❌:             TREAT_ERR( *err, 0, szBondType );
    // INCHI✔❌:             *err |= 8; /*  Unrecognized Bond stereo replaced with non-stereo bond */
    // INCHI✔❌:             cStereoType1 = 0;
    // INCHI✔❌:             cStereoType2 = 0;
    // INCHI✔❌:         }
    // INCHI✔❌:         break;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* neighbor */
    // INCHI✔❌:     if (ati[a1].neighbor[j] < 0 || ati[a1].neighbor[j] >= nNumAtoms)
    // INCHI✔❌:     {
    // INCHI✔❌:         *err |= 1; /*  bond for impossible atom number(s); ignored */
    // INCHI✔❌:         TREAT_ERR( *err, 0, "Bond to nonexistent atom" );
    // INCHI✔❌:         goto err_exit;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     a2 = (AT_NUMB) ati[a1].neighbor[j];
    // INCHI✔❌:     if (a2 == a1)
    // INCHI✔❌:     {
    // INCHI✔❌:         *err |= 1; /*  bond for impossible atom number(s); ignored */
    // INCHI✔❌:         TREAT_ERR( *err, 0, "Atom has a bond to itself" );
    // INCHI✔❌:         goto err_exit;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* consistency check; locate the bond in the opposite atom */
    // INCHI✔❌:     p1 = is_in_the_list( at[a1].neighbor, (AT_NUMB) a2, at[a1].valence );
    // INCHI✔❌:     p2 = is_in_the_list( at[a2].neighbor, (AT_NUMB) a1, at[a2].valence );
    // INCHI✔❌:
    // INCHI✔❌:     if (p1 && p2)
    // INCHI✔❌:     {
    // INCHI✔❌:         n1 = (int) ( p1 - at[a1].neighbor );
    // INCHI✔❌:         n2 = (int) ( p2 - at[a2].neighbor );
    // INCHI✔❌:         if ((n1 + 1 < at[a1].valence &&
    // INCHI✔❌:              is_in_the_list( at[a1].neighbor + n1 + 1, (AT_NUMB) a2, at[a1].valence - n1 - 1 ))
    // INCHI✔❌:              ||
    // INCHI✔❌:              (n2 + 1 < at[a2].valence &&
    // INCHI✔❌:              is_in_the_list( at[a2].neighbor + n2 + 1, (AT_NUMB) a1, at[a2].valence - n2 - 1 ))) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔❌:         {
    // INCHI✔❌:             TREAT_ERR( *err, 0, "Multiple bonds between two atoms" );
    // INCHI✔❌:             *err |= 2; /*  multiple bonds between atoms */
    // INCHI✔❌:         }
    // INCHI✔❌:         else if (n1 < at[a1].valence && n2 < at[a2].valence &&
    // INCHI✔❌:              cBondType == at[a2].bond_type[n2] &&
    // INCHI✔❌:              cBondType == at[a1].bond_type[n1] &&
    // INCHI✔❌:              cStereoType1 == at[a1].bond_stereo[n1] &&
    // INCHI✔❌:              cStereoType2 == at[a2].bond_stereo[n2])
    // INCHI✔❌:         {
    // INCHI✔❌:             /*TREAT_ERR (*err, 0, "Duplicated bond(s) between two atoms");*/
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             TREAT_ERR( *err, 0, "Multiple bonds between two atoms" );
    // INCHI✔❌:             *err |= 2; /*  multiple bonds between atoms */
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (( p1 || p2 ) &&
    // INCHI✔❌:         ( p1 || at[a1].valence < MAXVAL ) &&
    // INCHI✔❌:         ( p2 || at[a2].valence < MAXVAL ))
    // INCHI✔❌:     {
    // INCHI✔❌:         n1 = p1 ? (int) ( p1 - at[a1].neighbor ) : at[a1].valence++;
    // INCHI✔❌:         n2 = p2 ? (int) ( p2 - at[a2].neighbor ) : at[a2].valence++;
    // INCHI✔❌:         /* the bond is present in one atom only: possibly program error */
    // INCHI✔❌:         if ((p1 && ( cBondType != at[a1].bond_type[n1] || at[a1].bond_stereo[n1] != cStereoType1 )) ||
    // INCHI✔❌:              (p2 && ( cBondType != at[a2].bond_type[n2] || at[a2].bond_stereo[n2] != cStereoType2 ))) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔❌:         {
    // INCHI✔❌:             TREAT_ERR( *err, 0, "Multiple bonds between two atoms" );
    // INCHI✔❌:             *err |= 2; /*  multiple bonds between atoms */
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             TREAT_ERR( *err, 0, "Duplicated bond(s) between two atoms" );
    // INCHI✔❌:             /* warning */
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (!p1 && !p2 && at[a1].valence < MAXVAL && at[a2].valence < MAXVAL)
    // INCHI✔❌:     {
    // INCHI✔❌:         n1 = at[a1].valence++;
    // INCHI✔❌:         n2 = at[a2].valence++;
    // INCHI✔❌:         ( *nNumBonds )++;
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         char szMsg[64];
    // INCHI✔❌:         *err |= 4; /*  too large number of bonds. Some bonds ignored. */
    // INCHI✔❌:         sprintf(szMsg, "Atom '%s' has more than %d bonds",
    // INCHI✔❌:             at[a1].valence >= MAXVAL ? at[a1].elname : at[a2].elname, MAXVAL);
    // INCHI✔❌:         TREAT_ERR( *err, 0, szMsg );
    // INCHI✔❌:         goto err_exit;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* store the connection */
    // INCHI✔❌:
    // INCHI✔❌:     /* bond type */ /* djb-rwth: fixing buffer overruns */
    // INCHI✔❌:     if ((n1 < MAXVAL) && (n2 < MAXVAL))
    // INCHI✔❌:     {
    // INCHI✔❌:         at[a1].bond_type[n1] =
    // INCHI✔❌:             at[a2].bond_type[n2] = cBondType;
    // INCHI✔❌:         /* connection */
    // INCHI✔❌:         at[a1].neighbor[n1] = (AT_NUMB)a2;
    // INCHI✔❌:         at[a2].neighbor[n2] = (AT_NUMB)a1;
    // INCHI✔❌:         /* stereo */
    // INCHI✔❌:         at[a1].bond_stereo[n1] = cStereoType1; /*  >0: the wedge (pointed) end is at this atom */
    // INCHI✔❌:         at[a2].bond_stereo[n2] = cStereoType2; /*  <0: the wedge (pointed) end is at the opposite atom */
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         goto err_exit;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return 0;
    // INCHI✔❌:
    // INCHI✔❌: err_exit:
    // INCHI✔❌:
    // INCHI✔❌:     return 1;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: SetBondProperties

    let atom_index_usize =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let bond_index_usize =
        usize::try_from(bond_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let input = heap
        .slice(input_atoms.as_const())?
        .get(atom_index_usize)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let source_bond_type = *input
        .bond_type
        .get(bond_index_usize)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let source_stereo = *input
        .bond_stereo
        .get(bond_index_usize)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let source_neighbor = i32::from(
        *input
            .neighbor
            .get(bond_index_usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );

    let bond_type = match i32::from(source_bond_type) {
        value if value == tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i32 => BOND_TYPE_SINGLE as i8,
        value if value == tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE as i32 => BOND_TYPE_DOUBLE as i8,
        value if value == tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE as i32 => BOND_TYPE_TRIPLE as i8,
        value if value == tagINCHIBondType_INCHI_BOND_TYPE_ALTERN as i32 => BOND_TYPE_ALTERN as i8,
        value => {
            add_source_error(heap, error_text, "Unrecognized bond type:")?;
            add_source_error(heap, error_text, &value.to_string())?;
            *error.as_deref_mut().ok_or(SourceHeapError::NullPointer)? |= 8;
            BOND_TYPE_SINGLE as i8
        }
    };

    let (stereo1, stereo2) = match i32::from(source_stereo) {
        value if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE => (0, 0),
        value if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1UP => {
            (STEREO_SNGL_UP as i8, -(STEREO_SNGL_UP as i8))
        }
        value if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1EITHER => {
            (STEREO_SNGL_EITHER as i8, -(STEREO_SNGL_EITHER as i8))
        }
        value if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1DOWN => {
            (STEREO_SNGL_DOWN as i8, -(STEREO_SNGL_DOWN as i8))
        }
        value if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2UP => {
            (-(STEREO_SNGL_UP as i8), STEREO_SNGL_UP as i8)
        }
        value if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2EITHER => {
            (-(STEREO_SNGL_EITHER as i8), STEREO_SNGL_EITHER as i8)
        }
        value if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2DOWN => {
            (-(STEREO_SNGL_DOWN as i8), STEREO_SNGL_DOWN as i8)
        }
        value
            if value == tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER
                || value == -tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER =>
        {
            (STEREO_DBLE_EITHER as i8, STEREO_DBLE_EITHER as i8)
        }
        value => {
            add_source_error(heap, error_text, "Unrecognized bond stereo:")?;
            add_source_error(heap, error_text, &value.to_string())?;
            *error.as_deref_mut().ok_or(SourceHeapError::NullPointer)? |= 8;
            (0, 0)
        }
    };

    if source_neighbor < 0 || source_neighbor >= num_atoms {
        *error.as_deref_mut().ok_or(SourceHeapError::NullPointer)? |= 1;
        add_source_error(heap, error_text, "Bond to nonexistent atom")?;
        return Ok(1);
    }
    let neighbor_atom = source_neighbor as u16;
    if i32::from(neighbor_atom) == atom_index {
        *error.as_deref_mut().ok_or(SourceHeapError::NullPointer)? |= 1;
        add_source_error(heap, error_text, "Atom has a bond to itself")?;
        return Ok(1);
    }
    let neighbor_index = usize::from(neighbor_atom);
    let (atom1, atom2) = {
        let target = heap.slice(atoms.as_const())?;
        (
            target
                .get(atom_index_usize)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
            target
                .get(neighbor_index)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let valence1 = i32::from(atom1.valence);
    let valence2 = i32::from(atom2.valence);
    let position1 = is_in_the_list(Some(&atom1.neighbor), neighbor_atom, valence1)?;
    let position2 = is_in_the_list(Some(&atom2.neighbor), atom_index as u16, valence2)?;
    let max_valence = MAXVAL as i32;

    let (position1, position2) = if let (Some(position1), Some(position2)) = (position1, position2)
    {
        let index1 =
            i32::try_from(position1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let index2 =
            i32::try_from(position2).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let duplicate1 = index1 + 1 < valence1
            && is_in_the_list(
                Some(&atom1.neighbor[position1 + 1..]),
                neighbor_atom,
                valence1 - index1 - 1,
            )?
            .is_some();
        let duplicate2 = index2 + 1 < valence2
            && is_in_the_list(
                Some(&atom2.neighbor[position2 + 1..]),
                atom_index as u16,
                valence2 - index2 - 1,
            )?
            .is_some();
        if duplicate1 || duplicate2 {
            add_source_error(heap, error_text, "Multiple bonds between two atoms")?;
            *error.as_deref_mut().ok_or(SourceHeapError::NullPointer)? |= 2;
        } else if !(index1 < valence1
            && index2 < valence2
            && bond_type as u8 == atom2.bond_type[position2]
            && bond_type as u8 == atom1.bond_type[position1]
            && stereo1 == atom1.bond_stereo[position1]
            && stereo2 == atom2.bond_stereo[position2])
        {
            add_source_error(heap, error_text, "Multiple bonds between two atoms")?;
            *error.as_deref_mut().ok_or(SourceHeapError::NullPointer)? |= 2;
        }
        (position1, position2)
    } else if (position1.is_some() || position2.is_some())
        && (position1.is_some() || valence1 < max_valence)
        && (position2.is_some() || valence2 < max_valence)
    {
        let index1 = position1.unwrap_or(
            usize::try_from(valence1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
        );
        let index2 = position2.unwrap_or(
            usize::try_from(valence2).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
        );
        {
            let target = heap.slice_mut(atoms)?;
            if position1.is_none() {
                target[atom_index_usize].valence += 1;
            }
            if position2.is_none() {
                target[neighbor_index].valence += 1;
            }
        }
        let mismatch = position1.is_some()
            && (bond_type as u8 != atom1.bond_type[index1] || atom1.bond_stereo[index1] != stereo1)
            || position2.is_some()
                && (bond_type as u8 != atom2.bond_type[index2]
                    || atom2.bond_stereo[index2] != stereo2);
        if mismatch {
            add_source_error(heap, error_text, "Multiple bonds between two atoms")?;
            *error.as_deref_mut().ok_or(SourceHeapError::NullPointer)? |= 2;
        } else {
            add_source_error(heap, error_text, "Duplicated bond(s) between two atoms")?;
        }
        (index1, index2)
    } else if position1.is_none()
        && position2.is_none()
        && valence1 < max_valence
        && valence2 < max_valence
    {
        let index1 =
            usize::try_from(valence1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let index2 =
            usize::try_from(valence2).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        {
            let target = heap.slice_mut(atoms)?;
            target[atom_index_usize].valence += 1;
            target[neighbor_index].valence += 1;
        }
        let num_bonds = num_bonds.ok_or(SourceHeapError::NullPointer)?;
        *num_bonds = num_bonds
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        (index1, index2)
    } else {
        *error.as_deref_mut().ok_or(SourceHeapError::NullPointer)? |= 4;
        let element = if valence1 >= max_valence {
            &atom1.elname
        } else {
            &atom2.elname
        };
        let element_length = element
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let element = String::from_utf8(
            element[..element_length]
                .iter()
                .map(|byte| *byte as u8)
                .collect(),
        )
        .map_err(|_| SourceHeapError::InvalidSourceTextEncoding)?;
        add_source_error(
            heap,
            error_text,
            &format!("Atom '{element}' has more than {MAXVAL} bonds"),
        )?;
        return Ok(1);
    };

    if position1 >= MAXVAL as usize || position2 >= MAXVAL as usize {
        return Ok(1);
    }
    let target = heap.slice_mut(atoms)?;
    target[atom_index_usize].bond_type[position1] = bond_type as u8;
    target[neighbor_index].bond_type[position2] = bond_type as u8;
    target[atom_index_usize].neighbor[position1] = neighbor_atom;
    target[neighbor_index].neighbor[position2] = atom_index as u16;
    target[atom_index_usize].bond_stereo[position1] = stereo1;
    target[neighbor_index].bond_stereo[position2] = stereo2;
    Ok(0)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn SetAtomAndBondProperties(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    input_atoms: SourceMutPointer<inchi_Atom>,
    atom_index: i32,
    do_not_add_hydrogen: i32,
    error_text: Option<SourceMutPointer<i8>>,
    mut error: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1439 SetAtomAndBondProperties
    // INCHI✔❌: int SetAtomAndBondProperties( inp_ATOM *at,
    // INCHI✔❌:                               inchi_Atom *ati,
    // INCHI✔❌:                               int a1,
    // INCHI✔❌:                               int bDoNotAddH,
    // INCHI✔❌:                               char *pStrErr,
    // INCHI✔❌:                               int *err )
    // INCHI✔❌: {
    // INCHI✔❌:     int valence, chem_valence, num_alt_bonds, j, n1;
    // INCHI✔❌:     int nRadical, nCharge;
    // INCHI✔❌:     static int el_number_H = 0;
    // INCHI✔❌:
    // INCHI✔❌:     if (!el_number_H)
    // INCHI✔❌:     {
    // INCHI✔❌:         el_number_H = EL_NUMBER_H;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     nRadical = nCharge = 0;
    // INCHI✔❌:     valence = at[a1].valence;
    // INCHI✔❌:     chem_valence = num_alt_bonds = 0;
    // INCHI✔❌:     for (j = 0; j < valence; j++)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (at[a1].bond_type[j] <= BOND_TYPE_TRIPLE)
    // INCHI✔❌:         {
    // INCHI✔❌:             chem_valence += at[a1].bond_type[j];
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             num_alt_bonds++;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     switch (num_alt_bonds)
    // INCHI✔❌:     {
    // INCHI✔❌:         case 0:
    // INCHI✔❌:             break;
    // INCHI✔❌:         case 2:
    // INCHI✔❌:             chem_valence += 3; /* -C= */
    // INCHI✔❌:             break;
    // INCHI✔❌:         case 3:
    // INCHI✔❌:             chem_valence += 4;  /* >C= */
    // INCHI✔❌:             break;
    // INCHI✔❌:         default:
    // INCHI✔❌:         {
    // INCHI✔❌:             char szMsg[64];
    // INCHI✔❌:             *err |= 8; /*  wrong number of alt. bonds */
    // INCHI✔❌:             sprintf(szMsg, "Atom '%s' has %d alternating bonds",
    // INCHI✔❌:                 at[a1].elname, num_alt_bonds);
    // INCHI✔❌:             TREAT_ERR( *err, 0, szMsg );
    // INCHI✔❌:         }
    // INCHI✔❌:         break;
    // INCHI✔❌:     }
    // INCHI✔❌:     at[a1].chem_bonds_valence = chem_valence;
    // INCHI✔❌:
    // INCHI✔❌:     /* aliased hydrogen atoms */
    // INCHI✔❌:     if (ERR_ELEM == ( n1 = get_periodic_table_number( at[a1].elname ) ))
    // INCHI✔❌:     {
    // INCHI✔❌:         /*  Case when elname contains more than 1 element: extract number of H if possible */
    // INCHI✔❌:         if (extract_charges_and_radicals( at[a1].elname, &nRadical, &nCharge ))
    // INCHI✔❌:         {
    // INCHI✔❌:             if ((nRadical && at[a1].radical && nRadical != at[a1].radical) ||
    // INCHI✔❌:                  (nCharge && at[a1].charge && nCharge != at[a1].charge)) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔❌:             {
    // INCHI✔❌:                 TREAT_ERR( *err, 0, "Ignored charge/radical redefinition:" );
    // INCHI✔❌:                 TREAT_ERR( *err, 0, ati[a1].elname );
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 if (nRadical)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     at[a1].radical = nRadical;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (nCharge)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     at[a1].charge = nCharge;
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         at[a1].num_H = extract_H_atoms( at[a1].elname, at[a1].num_iso_H );
    // INCHI✔❌:         if (!at[a1].elname[0] && NUMH( at, a1 ))
    // INCHI✔❌:         {
    // INCHI✔❌:             /* alias contains only H. Added 2004-07-21, fixed 2004-07-22
    // INCHI✔❌:              * move the heaviest isotope to the "central atom"
    // INCHI✔❌:              * Note: this must be consistent with H-H treatment in remove_terminal_HDT()
    // INCHI✔❌:              */
    // INCHI✔❌:             strcpy( at[a1].elname, "H" );
    // INCHI✔❌:             if (NUM_ISO_H( at, a1 ))
    // INCHI✔❌:             {
    // INCHI✔❌:                 for (j = NUM_H_ISOTOPES - 1; 0 <= j; j--)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     if (at[a1].num_iso_H[j])
    // INCHI✔❌:                     {
    // INCHI✔❌:                         at[a1].num_iso_H[j] --;
    // INCHI✔❌:                         at[a1].iso_atw_diff = 1 + j;
    // INCHI✔❌:                         break;
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 at[a1].num_H--;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:
    // INCHI✔❌:         if (ERR_ELEM == ( n1 = get_periodic_table_number( at[a1].elname ) ))
    // INCHI✔❌:         {
    // INCHI✔❌:             n1 = 0;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (n1)
    // INCHI✔❌:         {
    // INCHI✔❌:             at[a1].at_type |= 1; /* "Aliased" atom: data in the element name */
    // INCHI✔❌:             TREAT_ERR( *err, 0, "Parsed compound atom(s):" );
    // INCHI✔❌:             TREAT_ERR( *err, 0, ati[a1].elname );
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     at[a1].el_number = (U_CHAR) n1;
    // INCHI✔❌:     if (!n1)
    // INCHI✔❌:     {
    // INCHI✔❌:         *err |= 64; /*  Unrecognized aromatic bond(s) replaced with single */
    // INCHI✔❌:         TREAT_ERR( *err, 0, "Unknown element(s):" );
    // INCHI✔❌:         TREAT_ERR( *err, 0, at[a1].elname );
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         /* replace explicit D or T with isotopic H (added 2003-06-02) */
    // INCHI✔❌:         if (el_number_H == n1 && !at[a1].iso_atw_diff)
    // INCHI✔❌:         {
    // INCHI✔❌:             switch (at[a1].elname[0])
    // INCHI✔❌:             {
    // INCHI✔❌:                 case 'D':
    // INCHI✔❌:                     at[a1].iso_atw_diff = 2;
    // INCHI✔❌:                     mystrncpy( at[a1].elname, "H", sizeof( at->elname ) );
    // INCHI✔❌:                     break;
    // INCHI✔❌:                 case 'T':
    // INCHI✔❌:                     at[a1].iso_atw_diff = 3;
    // INCHI✔❌:                     mystrncpy( at[a1].elname, "H", sizeof( at->elname ) );
    // INCHI✔❌:                     break;
    // INCHI✔❌:                 case 'H':
    // INCHI✔❌:                     if (1 <= ati[a1].isotopic_mass)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         AT_NUM iso_atw_diff;
    // INCHI✔❌:                         if (ISOTOPIC_SHIFT_FLAG - ISOTOPIC_SHIFT_MAX <= ati[a1].isotopic_mass &&
    // INCHI✔❌:                              ISOTOPIC_SHIFT_FLAG + ISOTOPIC_SHIFT_MAX >= ati[a1].isotopic_mass)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             /* ati[a1].isotopic_mass is isotopic iso_atw_diff + ISOTOPIC_SHIFT_FLAG */
    // INCHI✔❌:                             iso_atw_diff = ati[a1].isotopic_mass - ISOTOPIC_SHIFT_FLAG;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                         {
    // INCHI✔❌:                             /* ati[a1].isotopic_mass is isotopic mass */
    // INCHI✔❌:                             int iso_atw = get_atomic_mass_from_elnum( (int) at[a1].el_number );
    // INCHI✔❌:                             iso_atw_diff = ati[a1].isotopic_mass - iso_atw;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         if (iso_atw_diff >= 0)
    // INCHI✔❌:                             iso_atw_diff++;
    // INCHI✔❌:                         /* reproduce Bug04: allowed non-terminal H heavier than T */
    // INCHI✔❌:                         if (1 <= iso_atw_diff &&
    // INCHI✔❌:                             ( at[a1].valence != 1 || iso_atw_diff <= NUM_H_ISOTOPES ))
    // INCHI✔❌:                         {
    // INCHI✔❌:                             at[a1].iso_atw_diff = (S_CHAR) iso_atw_diff;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {/* isotopic shift */
    // INCHI✔❌:             if (ati[a1].isotopic_mass)
    // INCHI✔❌:             {
    // INCHI✔❌:                 AT_NUM iso_atw_diff;
    // INCHI✔❌:                 if (ISOTOPIC_SHIFT_FLAG - ISOTOPIC_SHIFT_MAX <= ati[a1].isotopic_mass &&
    // INCHI✔❌:                      ISOTOPIC_SHIFT_FLAG + ISOTOPIC_SHIFT_MAX >= ati[a1].isotopic_mass)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* ati[a1].isotopic_mass is isotopic iso_atw_diff + ISOTOPIC_SHIFT_FLAG */
    // INCHI✔❌:                     iso_atw_diff = ati[a1].isotopic_mass - ISOTOPIC_SHIFT_FLAG;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* ati[a1].isotopic_mass is isotopic mass */
    // INCHI✔❌:                     iso_atw_diff = get_atomic_mass_from_elnum( (int) at[a1].el_number );
    // INCHI✔❌:                     iso_atw_diff = ati[a1].isotopic_mass - iso_atw_diff;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 if (iso_atw_diff >= 0)
    // INCHI✔❌:                     iso_atw_diff++;
    // INCHI✔❌:                 at[a1].iso_atw_diff = (S_CHAR) iso_atw_diff;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* add implicit hydrogen atoms flag */
    // INCHI✔❌:     if (ati[a1].num_iso_H[0] == -1)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (!bDoNotAddH)
    // INCHI✔❌:         {
    // INCHI✔❌:             at[a1].at_type |= 2; /* user requested to add H */
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         at[a1].num_H = ati[a1].num_iso_H[0];
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     for (j = 0; j < NUM_H_ISOTOPES; j++)
    // INCHI✔❌:     {
    // INCHI✔❌:         at[a1].num_iso_H[j] = ati[a1].num_iso_H[j + 1];
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (num_alt_bonds)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* atom has aromatic bonds AND the chemical valence is not known */
    // INCHI✔❌:         int num_H = NUMH( at, a1 );
    // INCHI✔❌:         int chem_valence_alt = at[a1].chem_bonds_valence + num_H;
    // INCHI✔❌:         int bUnusualValenceArom =
    // INCHI✔❌:             detect_unusual_el_valence( (int) at[a1].el_number, at[a1].charge,
    // INCHI✔❌:                                         at[a1].radical, chem_valence_alt,
    // INCHI✔❌:                                         num_H, at[a1].valence );
    // INCHI✔❌:         int bUnusualValenceNoArom =
    // INCHI✔❌:             detect_unusual_el_valence( (int) at[a1].el_number, at[a1].charge,
    // INCHI✔❌:                                         at[a1].radical, chem_valence_alt - 1,
    // INCHI✔❌:                                         num_H, at[a1].valence );
    // INCHI✔❌:         if (bUnusualValenceArom && !bUnusualValenceNoArom && 0 == nBondsValToMetal( at, a1 ))
    // INCHI✔❌:         {
    // INCHI✔❌:             /* typically NH in 5-member aromatic ring */
    // INCHI✔❌:             at[a1].chem_bonds_valence--;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return 0;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: SetAtomAndBondProperties

    const EL_NUMBER_H: i32 = 1;
    let atom_index_usize =
        usize::try_from(atom_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let input = heap
        .slice(input_atoms.as_const())?
        .get(atom_index_usize)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut atom = heap
        .slice(atoms.as_const())?
        .get(atom_index_usize)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    let valence = i32::from(atom.valence);
    let mut chemical_valence = 0_i32;
    let mut alternating_bonds = 0_i32;
    for bond_index in 0..valence {
        let bond_index =
            usize::try_from(bond_index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let bond_type = i32::from(
            *atom
                .bond_type
                .get(bond_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        if bond_type <= BOND_TYPE_TRIPLE as i32 {
            chemical_valence = chemical_valence
                .checked_add(bond_type)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        } else {
            alternating_bonds = alternating_bonds
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        }
    }
    match alternating_bonds {
        0 => {}
        2 => {
            chemical_valence = chemical_valence
                .checked_add(3)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        }
        3 => {
            chemical_valence = chemical_valence
                .checked_add(4)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        }
        _ => {
            *error.as_deref_mut().ok_or(SourceHeapError::NullPointer)? |= 8;
            let element = fixed_source_text(&atom.elname)?;
            add_source_error(
                heap,
                error_text,
                &format!("Atom '{element}' has {alternating_bonds} alternating bonds"),
            )?;
        }
    }
    atom.chem_bonds_valence = chemical_valence as i8;

    let mut element_number = get_periodic_table_number(Some(&atom.elname))?;
    if element_number == ERR_ELEM {
        let mut parsed_radical = 0_i32;
        let mut parsed_charge = 0_i32;
        if extract_charges_and_radicals(
            Some(&mut atom.elname),
            Some(&mut parsed_radical),
            Some(&mut parsed_charge),
        )? != 0
        {
            if parsed_radical != 0 && atom.radical != 0 && parsed_radical != i32::from(atom.radical)
                || parsed_charge != 0 && atom.charge != 0 && parsed_charge != i32::from(atom.charge)
            {
                add_source_error(heap, error_text, "Ignored charge/radical redefinition:")?;
                add_source_error(heap, error_text, &fixed_source_text(&input.elname)?)?;
            } else {
                if parsed_radical != 0 {
                    atom.radical = parsed_radical as i8;
                }
                if parsed_charge != 0 {
                    atom.charge = parsed_charge as i8;
                }
            }
        }

        atom.num_H = extract_h_atoms(Some(&mut atom.elname), Some(&mut atom.num_iso_H))? as i8;
        let isotope_hydrogens = atom
            .num_iso_H
            .iter()
            .fold(0_i32, |sum, count| sum + i32::from(*count));
        if atom.elname[0] == 0 && i32::from(atom.num_H) + isotope_hydrogens != 0 {
            atom.elname[0] = b'H' as i8;
            atom.elname[1] = 0;
            if isotope_hydrogens != 0 {
                for isotope in (0..NUM_H_ISOTOPES as usize).rev() {
                    if atom.num_iso_H[isotope] != 0 {
                        atom.num_iso_H[isotope] = atom.num_iso_H[isotope].wrapping_sub(1);
                        atom.iso_atw_diff = (1 + isotope) as i8;
                        break;
                    }
                }
            } else {
                atom.num_H = atom.num_H.wrapping_sub(1);
            }
        }

        element_number = get_periodic_table_number(Some(&atom.elname))?;
        if element_number == ERR_ELEM {
            element_number = 0;
        }
        if element_number != 0 {
            atom.at_type |= 1;
            add_source_error(heap, error_text, "Parsed compound atom(s):")?;
            add_source_error(heap, error_text, &fixed_source_text(&input.elname)?)?;
        }
    }

    atom.el_number = element_number as u8;
    if element_number == 0 {
        *error.as_deref_mut().ok_or(SourceHeapError::NullPointer)? |= 64;
        add_source_error(heap, error_text, "Unknown element(s):")?;
        add_source_error(heap, error_text, &fixed_source_text(&atom.elname)?)?;
    } else if element_number == EL_NUMBER_H && atom.iso_atw_diff == 0 {
        match atom.elname[0] as u8 {
            b'D' => {
                atom.iso_atw_diff = 2;
                mystrncpy_slice(Some(&mut atom.elname), Some(&[b'H' as i8, 0]), 6)?;
            }
            b'T' => {
                atom.iso_atw_diff = 3;
                mystrncpy_slice(Some(&mut atom.elname), Some(&[b'H' as i8, 0]), 6)?;
            }
            b'H' if input.isotopic_mass >= 1 => {
                let mut isotope_difference =
                    isotope_mass_difference(atom.el_number, input.isotopic_mass)?;
                if isotope_difference >= 0 {
                    isotope_difference = isotope_difference.wrapping_add(1);
                }
                if isotope_difference >= 1
                    && (atom.valence != 1 || isotope_difference <= NUM_H_ISOTOPES as i16)
                {
                    atom.iso_atw_diff = isotope_difference as i8;
                }
            }
            _ => {}
        }
    } else if input.isotopic_mass != 0 {
        let mut isotope_difference = isotope_mass_difference(atom.el_number, input.isotopic_mass)?;
        if isotope_difference >= 0 {
            isotope_difference = isotope_difference.wrapping_add(1);
        }
        atom.iso_atw_diff = isotope_difference as i8;
    }

    if input.num_iso_H[0] == -1 {
        if do_not_add_hydrogen == 0 {
            atom.at_type |= 2;
        }
    } else {
        atom.num_H = input.num_iso_H[0];
    }
    atom.num_iso_H.copy_from_slice(&input.num_iso_H[1..]);

    if alternating_bonds != 0 {
        let hydrogen_count = i32::from(atom.num_H)
            + atom
                .num_iso_H
                .iter()
                .map(|count| i32::from(*count))
                .sum::<i32>();
        let aromatic_valence = i32::from(atom.chem_bonds_valence)
            .checked_add(hydrogen_count)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        let unusual_aromatic = detect_unusual_el_valence(
            i32::from(atom.el_number),
            i32::from(atom.charge),
            i32::from(atom.radical),
            aromatic_valence,
            hydrogen_count,
            i32::from(atom.valence),
        )?;
        let unusual_nonaromatic = detect_unusual_el_valence(
            i32::from(atom.el_number),
            i32::from(atom.charge),
            i32::from(atom.radical),
            aromatic_valence
                .checked_sub(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?,
            hydrogen_count,
            i32::from(atom.valence),
        )?;
        if unusual_aromatic != 0 && unusual_nonaromatic == 0 {
            heap.slice_mut(atoms)?[atom_index_usize] = atom.clone();
            if n_bonds_val_to_metal(Some(heap.slice(atoms.as_const())?), atom_index)? == 0 {
                atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_sub(1);
            }
        }
    }

    heap.slice_mut(atoms)?[atom_index_usize] = atom;
    Ok(0)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn InpAtom0DToInchiAtom(
    heap: &mut SourceHeap,
    input_atoms: SourceMutPointer<inp_ATOM>,
    input_atom_count: i32,
    output_atom_count: Option<&mut i16>,
    output_atoms: Option<&mut SourceMutPointer<inchi_Atom>>,
    output_stereo_count: Option<&mut i16>,
    output_stereo: Option<&mut SourceMutPointer<inchi_Stereo0D>>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1670 InpAtom0DToInchiAtom
    // INCHI✔❌: int InpAtom0DToInchiAtom( inp_ATOM *at,
    // INCHI✔❌:                           int num_inp_atoms,
    // INCHI✔❌:                           AT_NUM *num_atoms,
    // INCHI✔❌:                           inchi_Atom **atom,
    // INCHI✔❌:                           AT_NUM *num_stereo0D,
    // INCHI✔❌:                           inchi_Stereo0D **stereo0D )
    // INCHI✔❌: {
    // INCHI✔❌:     int num_stereo_centers, num_stereo_bonds, num_inp_stereo0D, i, m, m1, m2, n, ret = 0;
    // INCHI✔❌:
    // INCHI✔❌:     /* count stereobonds, allenes. cumulenes. and stereoatoms */
    // INCHI✔❌:     num_stereo_centers = num_stereo_bonds = ret = 0;
    // INCHI✔❌:
    // INCHI✔❌:     *atom = NULL;
    // INCHI✔❌:     *num_atoms = 0;
    // INCHI✔❌:     *stereo0D = NULL;
    // INCHI✔❌:     *num_stereo0D = 0;
    // INCHI✔❌:
    // INCHI✔❌:     for (i = 0; i < num_inp_atoms; i++)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (at[i].p_parity)
    // INCHI✔❌:         {
    // INCHI✔❌:             /* stereocenter */
    // INCHI✔❌:             num_stereo_centers++;
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             for (m = 0; m < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m]; m++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ;
    // INCHI✔❌:             }
    // INCHI✔❌:             num_stereo_bonds += m;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     num_stereo_bonds /= 2;
    // INCHI✔❌:     num_inp_stereo0D = num_stereo_bonds + num_stereo_centers;
    // INCHI✔❌:
    // INCHI✔❌:     if (num_inp_atoms > 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         *atom = (inchi_Atom *) inchi_calloc( num_inp_atoms, sizeof( ( *atom )[0] ) );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     *num_atoms = num_inp_atoms;
    // INCHI✔❌:
    // INCHI✔❌:     if (num_inp_stereo0D > 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         *stereo0D = (inchi_Stereo0D *) inchi_calloc( num_inp_stereo0D, sizeof( ( *stereo0D )[0] ) );
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if ((num_inp_atoms && !( *atom )) || (num_inp_stereo0D > 0 && !( *stereo0D ))) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔❌:     {
    // INCHI✔❌:         /* allocation failed */
    // INCHI✔❌:         ret = -1;
    // INCHI✔❌:         goto exit_function;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* copy atom properties */
    // INCHI✔❌:     for (i = 0; i < num_inp_atoms; i++)
    // INCHI✔❌:     {
    // INCHI✔❌:         ( *atom )[i].num_bonds = at[i].valence;
    // INCHI✔❌:         for (m = 0; m < at[i].valence; m++)
    // INCHI✔❌:         {
    // INCHI✔❌:             ( *atom )[i].bond_type[m] = at[i].bond_type[m];
    // INCHI✔❌:             ( *atom )[i].neighbor[m] = at[i].neighbor[m];
    // INCHI✔❌:         }
    // INCHI✔❌:         ( *atom )[i].charge = at[i].charge;
    // INCHI✔❌: #if USE_BCF
    // INCHI✔❌:         memcpy_s( ( *atom )[i].elname, ATOM_EL_LEN, at[i].elname, ATOM_EL_LEN ); /* djb-rwth: function replaced with its safe C11 variant */
    // INCHI✔❌: #else
    // INCHI✔❌:         memcpy( ( *atom )[i].elname, at[i].elname, ATOM_EL_LEN );
    // INCHI✔❌: #endif
    // INCHI✔❌:         if (at[i].iso_atw_diff)
    // INCHI✔❌:         {
    // INCHI✔❌:             ( *atom )[i].isotopic_mass = ISOTOPIC_SHIFT_FLAG + ( at[i].iso_atw_diff > 0 ? at[i].iso_atw_diff - 1 : at[i].iso_atw_diff );
    // INCHI✔❌:         }
    // INCHI✔❌:         ( *atom )[i].num_iso_H[0] = at[i].num_H;
    // INCHI✔❌:         for (m = 0; m < NUM_H_ISOTOPES; m++)
    // INCHI✔❌:         {
    // INCHI✔❌:             ( *atom )[i].num_iso_H[m + 1] = at[i].num_iso_H[m];
    // INCHI✔❌:         }
    // INCHI✔❌:         ( *atom )[i].radical = at[i].radical;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /* stereo */
    // INCHI✔❌:     for (i = n = 0; i < num_inp_atoms; i++)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (at[i].p_parity)
    // INCHI✔❌:         {
    // INCHI✔❌:             if (n < num_inp_stereo0D)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ( *stereo0D )[n].central_atom = i;
    // INCHI✔❌:                 ( *stereo0D )[n].parity = at[i].p_parity;
    // INCHI✔❌:                 ( *stereo0D )[n].type = INCHI_StereoType_Tetrahedral;
    // INCHI✔❌:                 for (m = 0; m < MAX_NUM_STEREO_ATOM_NEIGH; m++)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     ( *stereo0D )[n].neighbor[m] = at[i].p_orig_at_num[m] - 1;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 n++;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:             {
    // INCHI✔❌:                 ret |= 1;
    // INCHI✔❌:                 break;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:         {
    // INCHI✔❌:             for (m1 = 0; m1 < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m1]; m1++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 /* find the opposite atom at the other end of double bond, allene, or cumulene */
    // INCHI✔❌:                 int chain[12], len = 0, nxt_neigh, nxt, cur;
    // INCHI✔❌:                 cur = chain[len++] = i;
    // INCHI✔❌:                 nxt_neigh = at[cur].sb_ord[m1];
    // INCHI✔❌:
    // INCHI✔❌:                 do
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* add next atom */
    // INCHI✔❌:                     chain[len++] = nxt = at[cur].neighbor[nxt_neigh];
    // INCHI✔❌:                     nxt_neigh = ( at[nxt].neighbor[0] == cur );
    // INCHI✔❌:                     cur = nxt;
    // INCHI✔❌:                     /* find nxt_neigh */
    // INCHI✔❌:                 }
    // INCHI✔❌:                 while (!at[cur].sb_parity[0] &&
    // INCHI✔❌:                        len < 12 &&
    // INCHI✔❌:                        at[cur].valence == 2);
    // INCHI✔❌:
    // INCHI✔❌:                 if (at[cur].sb_parity[0] && len <= 4 && i < cur /* count bonds only one time */)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     /* double bond, cumulene, or allene has been found */
    // INCHI✔❌:                     for (m2 = 0; m2 < MAX_NUM_STEREO_BONDS && at[cur].sb_parity[m2]; m2++)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         if (chain[len - 2] == at[cur].neighbor[(int) at[cur].sb_ord[m2]])
    // INCHI✔❌:                         {
    // INCHI✔❌:                             if (n < num_inp_stereo0D)
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 int parity1 = at[i].sb_parity[m1];
    // INCHI✔❌:                                 int parity2 = at[cur].sb_parity[m2];
    // INCHI✔❌:                                 int parity;
    // INCHI✔❌:                                 if (( INCHI_PARITY_ODD == parity1 || INCHI_PARITY_EVEN == parity1 ) &&
    // INCHI✔❌:                                     ( INCHI_PARITY_ODD == parity2 || INCHI_PARITY_EVEN == parity2 ))
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     /* well-defined parity */
    // INCHI✔❌:                                     parity = ( parity1 == parity2 ) ? INCHI_PARITY_EVEN : INCHI_PARITY_ODD;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 else
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     parity = inchi_max( parity1, parity2 );
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 ( *stereo0D )[n].central_atom = ( len == 3 ) ? chain[1] : NO_ATOM;
    // INCHI✔❌:                                 ( *stereo0D )[n].parity = parity;
    // INCHI✔❌:                                 ( *stereo0D )[n].type = len == 3 ? INCHI_StereoType_Allene : INCHI_StereoType_DoubleBond;
    // INCHI✔❌:                                 ( *stereo0D )[n].neighbor[0] = at[i].sn_orig_at_num[m1] - 1;
    // INCHI✔❌:                                 ( *stereo0D )[n].neighbor[1] = i;
    // INCHI✔❌:                                 ( *stereo0D )[n].neighbor[2] = cur;
    // INCHI✔❌:                                 ( *stereo0D )[n].neighbor[3] = at[cur].sn_orig_at_num[m2] - 1;
    // INCHI✔❌:                                 n++;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             else
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 ret |= 1;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             break;
    // INCHI✔❌:                         }
    // INCHI✔❌:                     }
    // INCHI✔❌:                 }
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     *num_stereo0D = n;
    // INCHI✔❌:
    // INCHI✔❌: exit_function:
    // INCHI✔❌:     if (ret < 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (*atom)
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_free( *atom );
    // INCHI✔❌:         }
    // INCHI✔❌:         if (*stereo0D)
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_free( *stereo0D );
    // INCHI✔❌:         }
    // INCHI✔❌:         *atom = NULL;
    // INCHI✔❌:         *stereo0D = NULL;
    // INCHI✔❌:         *num_atoms = 0;
    // INCHI✔❌:         *num_stereo0D = 0;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return ret;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: InpAtom0DToInchiAtom

    let output_atoms = output_atoms.ok_or(SourceHeapError::NullPointer)?;
    let output_atom_count = output_atom_count.ok_or(SourceHeapError::NullPointer)?;
    let output_stereo = output_stereo.ok_or(SourceHeapError::NullPointer)?;
    let output_stereo_count = output_stereo_count.ok_or(SourceHeapError::NullPointer)?;
    *output_atoms = SourceMutPointer::null();
    *output_atom_count = 0;
    *output_stereo = SourceMutPointer::null();
    *output_stereo_count = 0;

    let inputs = if input_atom_count > 0 {
        let count = usize::try_from(input_atom_count)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        heap.slice(input_atoms.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    } else {
        Vec::new()
    };

    let mut stereo_centers = 0_i32;
    let mut stereo_bond_ends = 0_i32;
    for input in &inputs {
        if input.p_parity != 0 {
            stereo_centers = stereo_centers
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        } else {
            let mut count = 0_i32;
            while count < MAX_NUM_STEREO_BONDS as i32 && input.sb_parity[count as usize] != 0 {
                count += 1;
            }
            stereo_bond_ends = stereo_bond_ends
                .checked_add(count)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
        }
    }
    let input_stereo_count = stereo_bond_ends / 2 + stereo_centers;

    if input_atom_count > 0 {
        *output_atoms = match inchi_calloc(heap, input_atom_count as u64, SOURCE_SIZEOF_INCHI_ATOM)
        {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    }
    *output_atom_count = input_atom_count as i16;
    if input_stereo_count > 0 {
        *output_stereo = match inchi_calloc(
            heap,
            input_stereo_count as u64,
            SOURCE_SIZEOF_INCHI_STEREO0D,
        ) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    }

    if input_atom_count != 0 && output_atoms.is_null()
        || input_stereo_count > 0 && output_stereo.is_null()
    {
        if !output_atoms.is_null() {
            inchi_free(heap, *output_atoms)?;
        }
        if !output_stereo.is_null() {
            inchi_free(heap, *output_stereo)?;
        }
        *output_atoms = SourceMutPointer::null();
        *output_stereo = SourceMutPointer::null();
        *output_atom_count = 0;
        *output_stereo_count = 0;
        return Ok(-1);
    }

    if input_atom_count > 0 {
        let atoms = heap.slice_mut(*output_atoms)?;
        for (index, input) in inputs.iter().enumerate() {
            let output = &mut atoms[index];
            output.num_bonds = i16::from(input.valence);
            for bond in 0..i32::from(input.valence) {
                let bond =
                    usize::try_from(bond).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                output.bond_type[bond] = input.bond_type[bond] as i8;
                output.neighbor[bond] = input.neighbor[bond] as i16;
            }
            output.charge = input.charge;
            output.elname.copy_from_slice(&input.elname);
            if input.iso_atw_diff != 0 {
                let adjustment = if input.iso_atw_diff > 0 {
                    i32::from(input.iso_atw_diff) - 1
                } else {
                    i32::from(input.iso_atw_diff)
                };
                output.isotopic_mass = (ISOTOPIC_SHIFT_FLAG as i32 + adjustment) as i16;
            }
            output.num_iso_H[0] = input.num_H;
            output.num_iso_H[1..].copy_from_slice(&input.num_iso_H);
            output.radical = input.radical;
        }
    }

    let mut return_value = 0_i32;
    let mut written_stereo = 0_i32;
    for (index, input) in inputs.iter().enumerate() {
        if input.p_parity != 0 {
            if written_stereo < input_stereo_count {
                let stereo = &mut heap.slice_mut(*output_stereo)?[written_stereo as usize];
                stereo.central_atom = index as i16;
                stereo.parity = input.p_parity;
                stereo.type_ = tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8;
                for neighbor in 0..MAX_NUM_STEREO_ATOM_NEIGH as usize {
                    stereo.neighbor[neighbor] =
                        (i32::from(input.p_orig_at_num[neighbor]) - 1) as i16;
                }
                written_stereo += 1;
            } else {
                return_value |= 1;
                break;
            }
            continue;
        }

        for first_order in 0..MAX_NUM_STEREO_BONDS as usize {
            if input.sb_parity[first_order] == 0 {
                break;
            }
            let mut chain = [0_i32; 12];
            let mut chain_length = 1_usize;
            let mut current = index;
            chain[0] = index as i32;
            let mut next_neighbor = i32::from(input.sb_ord[first_order]);
            loop {
                let next_neighbor_index = usize::try_from(next_neighbor)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let next = usize::from(
                    *inputs[current]
                        .neighbor
                        .get(next_neighbor_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let chain_slot = chain
                    .get_mut(chain_length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                *chain_slot = next as i32;
                chain_length += 1;
                let next_atom = inputs
                    .get(next)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                next_neighbor = i32::from(next_atom.neighbor[0] == current as u16);
                current = next;
                if next_atom.sb_parity[0] != 0 || chain_length >= 12 || next_atom.valence != 2 {
                    break;
                }
            }

            if inputs[current].sb_parity[0] == 0 || chain_length > 4 || index >= current {
                continue;
            }
            for second_order in 0..MAX_NUM_STEREO_BONDS as usize {
                if inputs[current].sb_parity[second_order] == 0 {
                    break;
                }
                let opposite_order = usize::try_from(inputs[current].sb_ord[second_order])
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let opposite_neighbor = *inputs[current]
                    .neighbor
                    .get(opposite_order)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if chain[chain_length - 2] != i32::from(opposite_neighbor) {
                    continue;
                }
                if written_stereo < input_stereo_count {
                    let first_parity = i32::from(input.sb_parity[first_order]);
                    let second_parity = i32::from(inputs[current].sb_parity[second_order]);
                    let parity = if matches!(
                        first_parity,
                        value if value == tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32
                            || value == tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i32
                    ) && matches!(
                        second_parity,
                        value if value == tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32
                            || value == tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i32
                    ) {
                        if first_parity == second_parity {
                            tagINCHIStereoParity0D_INCHI_PARITY_EVEN as i32
                        } else {
                            tagINCHIStereoParity0D_INCHI_PARITY_ODD as i32
                        }
                    } else {
                        first_parity.max(second_parity)
                    };
                    let stereo = &mut heap.slice_mut(*output_stereo)?[written_stereo as usize];
                    stereo.central_atom = if chain_length == 3 {
                        chain[1] as i16
                    } else {
                        NO_ATOM as i16
                    };
                    stereo.parity = parity as i8;
                    stereo.type_ = if chain_length == 3 {
                        tagINCHIStereoType0D_INCHI_StereoType_Allene as i8
                    } else {
                        tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8
                    };
                    stereo.neighbor[0] = (i32::from(input.sn_orig_at_num[first_order]) - 1) as i16;
                    stereo.neighbor[1] = index as i16;
                    stereo.neighbor[2] = current as i16;
                    stereo.neighbor[3] =
                        (i32::from(inputs[current].sn_orig_at_num[second_order]) - 1) as i16;
                    written_stereo += 1;
                } else {
                    return_value |= 1;
                }
                break;
            }
        }
    }
    *output_stereo_count = written_stereo as i16;
    Ok(return_value)
}

fn fixed_source_text(bytes: &[i8]) -> Result<String, SourceHeapError> {
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    String::from_utf8(bytes[..length].iter().map(|byte| *byte as u8).collect())
        .map_err(|_| SourceHeapError::InvalidSourceTextEncoding)
}

fn isotope_mass_difference(element_number: u8, isotopic_mass: i16) -> Result<i16, SourceHeapError> {
    let mass = i32::from(isotopic_mass);
    let lower = i32::try_from(ISOTOPIC_SHIFT_FLAG - ISOTOPIC_SHIFT_MAX)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let upper = i32::try_from(ISOTOPIC_SHIFT_FLAG + ISOTOPIC_SHIFT_MAX)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let difference = if lower <= mass && mass <= upper {
        mass - ISOTOPIC_SHIFT_FLAG as i32
    } else {
        mass - get_atomic_mass_from_elnum(i32::from(element_number))?
    };
    Ok(difference as i16)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn ExtractOneStructure(
    heap: &mut SourceHeap,
    structure: &mut STRUCT_DATA,
    input_parameters: &mut INPUT_PARMS,
    _title: Option<SourceMutPointer<i8>>,
    input: Option<&inchi_InputEx>,
    log_file: Option<&mut INCHI_IOSTREAM>,
    output_file: Option<&mut INCHI_IOSTREAM>,
    problem_file: Option<&mut INCHI_IOSTREAM>,
    original_input: &mut ORIG_ATOM_DATA,
    input_number: &mut i64,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1863 ExtractOneStructure
    // INCHI✔️❌: int ExtractOneStructure( STRUCT_DATA *sd,
    // INCHI✔️❌:                          INPUT_PARMS *ip,
    // INCHI✔️❌:                          char *szTitle,
    // INCHI✔️❌:                          inchi_InputEx *inp,
    // INCHI✔️❌:                          INCHI_IOSTREAM *log_file,
    // INCHI✔️❌:                          INCHI_IOSTREAM *out_file,
    // INCHI✔️❌:                          INCHI_IOSTREAM *prb_file,
    // INCHI✔️❌:                          ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️❌:                          long *num_inp )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int         *err = &sd->nStructReadError;
    // INCHI✔️❌:     char        *pStrErr = sd->pStrErrStruct;
    // INCHI✔️❌:     inp_ATOM    *at = NULL;
    // INCHI✔️❌:     MOL_COORD   *szCoord = NULL;
    // INCHI✔️❌:     inchi_Atom  *ati = NULL;
    // INCHI✔️❌:     int       nNumAtoms = 0;
    // INCHI✔️❌:     int       a1, j, valence, nDim, nNumBonds, nRet = 0, max_num_at;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* vABParityUnknown holds actual value of an internal constant signifying       */
    // INCHI✔️❌:     /* unknown parity: either the same as for undefined parity (default==standard)  */
    // INCHI✔️❌:     /*  or a specific one (non-std; requested by SLUUD switch).                     */
    // INCHI✔️❌:     int vABParityUnknown = AB_PARITY_UNDF;
    // INCHI✔️❌:     if (0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Make labels for unknown and undefined stereo different */
    // INCHI✔️❌:         vABParityUnknown = AB_PARITY_UNKN;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /********************************************************
    // INCHI✔️❌:      *
    // INCHI✔️❌:      *   Extract the structure
    // INCHI✔️❌:      *
    // INCHI✔️❌:      ********************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     FreeOrigAtData( orig_inp_data );
    // INCHI✔️❌:     nDim = 0;
    // INCHI✔️❌:     nNumBonds = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!inp || ( nNumAtoms = inp->num_atoms ) <= 0 || !( ati = inp->atom ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         TREAT_ERR( *err, 0, "Empty structure" );
    // INCHI✔️❌:         *err = 98;
    // INCHI✔️❌:         goto err_exit;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     max_num_at = ip->bLargeMolecules ? MAX_ATOMS : NORMALLY_ALLOWED_INP_MAX_ATOMS;
    // INCHI✔️❌:     if (nNumAtoms >= max_num_at)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         TREAT_ERR( *err, 0, "Too many atoms [did you forget 'LargeMolecules' switch?]" );
    // INCHI✔️❌:         *err = 70;
    // INCHI✔️❌:         orig_inp_data->num_inp_atoms = -1;
    // INCHI✔️❌:         goto err_exit;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     at = (inp_ATOM  *) inchi_calloc( nNumAtoms, sizeof( at[0] ) );
    // INCHI✔️❌:     szCoord = (MOL_COORD *) inchi_calloc( inchi_max( nNumAtoms, 1 ), sizeof( MOL_COORD ) );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!at || !szCoord)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔️❌:         *err = -1;
    // INCHI✔️❌:         goto err_exit;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /********************************************************
    // INCHI✔️❌:      *
    // INCHI✔️❌:      *   Extract typical for Molfile structural data
    // INCHI✔️❌:      *
    // INCHI✔️❌:      ********************************************************/
    // INCHI✔️❌:     /* extract atoms and bonds */
    // INCHI✔️❌:     for (a1 = 0; a1 < nNumAtoms; a1++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* extract atoms */
    // INCHI✔️❌:         SetAtomProperties( at, szCoord, ati, a1, &nDim, pStrErr, err );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (*err)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto err_exit;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* extract connections */
    // INCHI✔️❌:         valence = ati[a1].num_bonds;
    // INCHI✔️❌:         for (j = 0; j < valence; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             SetBondProperties( at, ati, a1, j, nNumAtoms, &nNumBonds, pStrErr, err );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (*err)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto err_exit;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     orig_inp_data->num_inp_atoms = nNumAtoms;
    // INCHI✔️❌:     orig_inp_data->num_inp_bonds = nNumBonds;
    // INCHI✔️❌:     orig_inp_data->num_dimensions = nDim;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* extract elements, chemical valences, implicit H, isotopic shifts */
    // INCHI✔️❌:     for (a1 = 0; a1 < nNumAtoms; a1++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* set temp flags in at[a1].at_type */
    // INCHI✔️❌:         /* (1: data in atom name; 2: request to add H) */
    // INCHI✔️❌:         SetAtomAndBondProperties( at,
    // INCHI✔️❌:                                   ati,
    // INCHI✔️❌:                                   a1,
    // INCHI✔️❌:                                   ip->bDoNotAddH,
    // INCHI✔️❌:                                   pStrErr,
    // INCHI✔️❌:                                   err );
    // INCHI✔️❌:         if (*err)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto err_exit;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* clear temp flags in at[].at_type; add implicit H */
    // INCHI✔️❌:     SetNumImplicitH( at, nNumAtoms );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (*err)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto err_exit;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /********************************************************
    // INCHI✔️❌:      *
    // INCHI✔️❌:      *   Extract the 0D parities (typical for CML)
    // INCHI✔️❌:      *
    // INCHI✔️❌:      ********************************************************/
    // INCHI✔️❌:     Extract0DParities( at,
    // INCHI✔️❌:                        nNumAtoms,
    // INCHI✔️❌:                        inp->stereo0D,
    // INCHI✔️❌:                        inp->num_stereo0D,
    // INCHI✔️❌:                        pStrErr,
    // INCHI✔️❌:                        err,
    // INCHI✔️❌:                        vABParityUnknown );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (*err)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto err_exit;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     orig_inp_data->at = at;
    // INCHI✔️❌:     at = NULL;
    // INCHI✔️❌:     orig_inp_data->num_dimensions = nDim;
    // INCHI✔️❌:     orig_inp_data->num_inp_atoms = nNumAtoms;
    // INCHI✔️❌:     orig_inp_data->num_inp_bonds = nNumBonds;
    // INCHI✔️❌:     orig_inp_data->szCoord = szCoord;
    // INCHI✔️❌:     szCoord = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* chiral flag */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* *****************************************************************************
    // INCHI✔️❌:      * Chiral flags are set in:
    // INCHI✔️❌:      * - ReadTheStructure() inchi-1, wInChI
    // INCHI✔️❌:      * - e_IchiMain.c -- main()               -- C example of calling InChI dll
    // INCHI✔️❌:      * - inchi_dll.c  ExtractOneStructure -- InChI dll code (here)
    // INCHI✔️❌:      *******************************************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:     if (( ip->nMode & REQ_MODE_CHIR_FLG_STEREO ) && ( ip->nMode & REQ_MODE_STEREO ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (ip->bChiralFlag & FLAG_SET_INP_AT_CHIRAL)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* absolute stereo */
    // INCHI✔️❌:             ip->nMode &= ~( REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO );
    // INCHI✔️❌:             sd->bChiralFlag &= ~FLAG_INP_AT_NONCHIRAL;
    // INCHI✔️❌:             sd->bChiralFlag |= FLAG_INP_AT_CHIRAL; /* write AuxInfo as chiral */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         /*if ( ip->bChiralFlag & FLAG_SET_INP_AT_NONCHIRAL )*/
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* relative stereo */
    // INCHI✔️❌:             ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO );
    // INCHI✔️❌:             ip->nMode |= REQ_MODE_RELATIVE_STEREO;
    // INCHI✔️❌:             sd->bChiralFlag &= ~FLAG_INP_AT_CHIRAL;
    // INCHI✔️❌:             sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL; /* write AuxInfo as non-chiral */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (ip->bChiralFlag & FLAG_SET_INP_AT_CHIRAL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sd->bChiralFlag &= ~FLAG_INP_AT_NONCHIRAL;
    // INCHI✔️❌:         sd->bChiralFlag |= FLAG_INP_AT_CHIRAL; /* write AuxInfo as chiral */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (ip->bChiralFlag & FLAG_SET_INP_AT_NONCHIRAL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sd->bChiralFlag &= ~FLAG_INP_AT_CHIRAL;
    // INCHI✔️❌:         sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL; /* write AuxInfo as non-chiral */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* v. 1.05 extensions  */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int res = SetExtOrigAtDataByInChIExtInput( &orig_inp_data->polymer,
    // INCHI✔️❌:                                                    &orig_inp_data->v3000,
    // INCHI✔️❌:                                                    inp->polymer,
    // INCHI✔️❌:                                                    inp->v3000,
    // INCHI✔️❌:                                                    orig_inp_data->num_inp_atoms );
    // INCHI✔️❌:         if (res)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             TREAT_ERR( res, 0, "General error on treating polymers" );
    // INCHI✔️❌:             *err = -1;
    // INCHI✔️❌:             goto err_exit;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     *num_inp += 1;
    // INCHI✔️❌:
    // INCHI✔️❌: err_exit:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at)
    // INCHI✔️❌:     {   /* if not moved to orig_inp_data/then nullified */
    // INCHI✔️❌:         inchi_free( at );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (szCoord)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( szCoord );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nRet = TreatErrorsInReadTheStructure( sd, ip, LOG_MASK_NO_WARN, NULL,
    // INCHI✔️❌:                                           log_file, out_file, prb_file,
    // INCHI✔️❌:                                           orig_inp_data, num_inp );
    // INCHI✔️❌:
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ExtractOneStructure

    let unknown_parity = if input_parameters.nMode & REQ_MODE_DIFF_UU_STEREO as u64 != 0 {
        AB_PARITY_UNKN as i32
    } else {
        AB_PARITY_UNDF as i32
    };
    FreeOrigAtData(heap, Some(original_input))?;
    let mut dimensions = 0;
    let mut num_bonds = 0;
    let mut atoms = SourceMutPointer::null();
    let mut coordinates = SourceMutPointer::null();
    let error_buffer = heap.allocate_model_storage(structure.pStrErrStruct.to_vec())?;

    'extract: {
        let Some(input) = input else {
            add_source_error(heap, Some(error_buffer), "Empty structure")?;
            structure.nStructReadError = 98;
            break 'extract;
        };
        let num_atoms = i32::from(input.num_atoms);
        if num_atoms <= 0 || input.atom.is_null() {
            add_source_error(heap, Some(error_buffer), "Empty structure")?;
            structure.nStructReadError = 98;
            break 'extract;
        }
        let maximum = if input_parameters.bLargeMolecules != 0 {
            MAX_ATOMS as i32
        } else {
            NORMALLY_ALLOWED_INP_MAX_ATOMS as i32
        };
        if num_atoms >= maximum {
            add_source_error(
                heap,
                Some(error_buffer),
                "Too many atoms [did you forget 'LargeMolecules' switch?]",
            )?;
            structure.nStructReadError = 70;
            original_input.num_inp_atoms = -1;
            break 'extract;
        }

        atoms = inchi_calloc::<inp_ATOM>(
            heap,
            num_atoms as u64,
            std::mem::size_of::<inp_ATOM>() as u64,
        )
        .unwrap_or_else(|_| SourceMutPointer::null());
        coordinates = inchi_calloc::<MOL_COORD>(heap, num_atoms.max(1) as u64, 32)
            .unwrap_or_else(|_| SourceMutPointer::null());
        if atoms.is_null() || coordinates.is_null() {
            add_source_error(heap, Some(error_buffer), "Out of RAM")?;
            structure.nStructReadError = -1;
            break 'extract;
        }

        for atom_index in 0..num_atoms {
            SetAtomProperties(
                heap,
                atoms,
                coordinates,
                input.atom,
                atom_index,
                Some(&mut dimensions),
                Some(error_buffer),
                Some(&mut structure.nStructReadError),
            )?;
            if structure.nStructReadError != 0 {
                break 'extract;
            }
            let valence =
                i32::from(heap.slice(input.atom.as_const())?[atom_index as usize].num_bonds);
            for bond_index in 0..valence {
                SetBondProperties(
                    heap,
                    atoms,
                    input.atom,
                    atom_index,
                    bond_index,
                    num_atoms,
                    Some(&mut num_bonds),
                    Some(error_buffer),
                    Some(&mut structure.nStructReadError),
                )?;
            }
            if structure.nStructReadError != 0 {
                break 'extract;
            }
        }

        original_input.num_inp_atoms = num_atoms;
        original_input.num_inp_bonds = num_bonds;
        original_input.num_dimensions = dimensions;
        for atom_index in 0..num_atoms {
            SetAtomAndBondProperties(
                heap,
                atoms,
                input.atom,
                atom_index,
                input_parameters.bDoNotAddH,
                Some(error_buffer),
                Some(&mut structure.nStructReadError),
            )?;
            if structure.nStructReadError != 0 {
                break 'extract;
            }
        }
        SetNumImplicitH(heap, atoms, num_atoms)?;
        if structure.nStructReadError != 0 {
            break 'extract;
        }
        Extract0DParities(
            heap,
            atoms,
            num_atoms,
            input.stereo0D,
            i32::from(input.num_stereo0D),
            Some(error_buffer),
            Some(&mut structure.nStructReadError),
            unknown_parity,
        )?;
        if structure.nStructReadError != 0 {
            break 'extract;
        }

        original_input.at = atoms;
        atoms = SourceMutPointer::null();
        original_input.num_dimensions = dimensions;
        original_input.num_inp_atoms = num_atoms;
        original_input.num_inp_bonds = num_bonds;
        original_input.szCoord = coordinates;
        coordinates = SourceMutPointer::null();

        if input_parameters.nMode & REQ_MODE_CHIR_FLG_STEREO as u64 != 0
            && input_parameters.nMode & REQ_MODE_STEREO as u64 != 0
        {
            if input_parameters.bChiralFlag & FLAG_SET_INP_AT_CHIRAL as i32 != 0 {
                input_parameters.nMode &=
                    !((REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO) as u64);
                structure.bChiralFlag &= !(FLAG_INP_AT_NONCHIRAL as i32);
                structure.bChiralFlag |= FLAG_INP_AT_CHIRAL as i32;
            } else {
                input_parameters.nMode &= !(REQ_MODE_RACEMIC_STEREO as u64);
                input_parameters.nMode |= REQ_MODE_RELATIVE_STEREO as u64;
                structure.bChiralFlag &= !(FLAG_INP_AT_CHIRAL as i32);
                structure.bChiralFlag |= FLAG_INP_AT_NONCHIRAL as i32;
            }
        } else if input_parameters.bChiralFlag & FLAG_SET_INP_AT_CHIRAL as i32 != 0 {
            structure.bChiralFlag &= !(FLAG_INP_AT_NONCHIRAL as i32);
            structure.bChiralFlag |= FLAG_INP_AT_CHIRAL as i32;
        } else if input_parameters.bChiralFlag
            & crate::source_types::FLAG_SET_INP_AT_NONCHIRAL as i32
            != 0
        {
            structure.bChiralFlag &= !(FLAG_INP_AT_CHIRAL as i32);
            structure.bChiralFlag |= FLAG_INP_AT_NONCHIRAL as i32;
        }

        let extension_result = SetExtOrigAtDataByInChIExtInput(
            heap,
            &mut original_input.polymer,
            &mut original_input.v3000,
            input.polymer,
            input.v3000,
            original_input.num_inp_atoms,
        )?;
        if extension_result != 0 {
            add_source_error(
                heap,
                Some(error_buffer),
                "General error on treating polymers",
            )?;
            structure.nStructReadError = -1;
            break 'extract;
        }
        *input_number = input_number
            .checked_add(1)
            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    }

    if !atoms.is_null() {
        inchi_free(heap, atoms)?;
    }
    if !coordinates.is_null() {
        inchi_free(heap, coordinates)?;
    }
    let error_bytes = heap.slice(error_buffer.as_const())?;
    if error_bytes.len() < structure.pStrErrStruct.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    structure.pStrErrStruct.copy_from_slice(&error_bytes[..256]);
    heap.free(error_buffer)?;

    TreatErrorsInReadTheStructure(
        heap,
        structure,
        input_parameters,
        LOG_MASK_NO_WARN as i32,
        None,
        log_file,
        output_file,
        problem_file,
        original_input,
        input_number,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::STR_ERR_LEN;
    use crate::test_support::allocate_source_fixture;

    #[test]
    fn source_port__inchi_dll__setnumimplicith__line_994() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            SetNumImplicitH(&mut heap, SourceMutPointer::null(), i32::MIN),
            Ok(())
        );
        assert_eq!(
            SetNumImplicitH(&mut heap, SourceMutPointer::null(), 0),
            Ok(())
        );

        let mut recalculate = inp_ATOM::default();
        recalculate.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
        recalculate.el_number = 6;
        recalculate.chem_bonds_valence = 1;
        recalculate.at_type = 2;

        let mut preserve = inp_ATOM::default();
        preserve.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
        preserve.el_number = 6;
        preserve.num_H = 5;
        preserve.at_type = 0;

        let mut aliased = inp_ATOM::default();
        aliased.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
        aliased.el_number = 6;
        aliased.num_H = 4;
        aliased.at_type = 1;

        let mut metal = inp_ATOM::default();
        metal.elname[..3].copy_from_slice(&[b'S' as i8, b'c' as i8, 0]);
        metal.el_number = 21;
        metal.num_H = 2;
        metal.at_type = 2;

        let atoms = heap
            .allocate(vec![recalculate, preserve, aliased, metal])
            .unwrap();
        SetNumImplicitH(&mut heap, atoms, 4).unwrap();
        let output = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(
            output.iter().map(|atom| atom.num_H).collect::<Vec<_>>(),
            vec![3, 5, 4, 2]
        );
        assert!(output.iter().all(|atom| atom.at_type == 0));
    }

    fn parse_case(
        text: &[u8],
        max_arguments: i32,
    ) -> (i32, Vec<String>, Vec<SourceArgvPointer>, Vec<i64>, Vec<u8>) {
        let mut heap = SourceHeap::default();
        let command = allocate_source_fixture(
            &mut heap,
            text.iter().map(|byte| *byte as i8).collect::<Vec<_>>(),
        );
        let mut arguments = vec![SourceArgvPointer::Null; max_arguments as usize];
        let count =
            parse_options_string(&mut heap, command, &mut arguments, max_arguments).unwrap();
        let values = arguments[..count as usize]
            .iter()
            .map(|argument| match argument {
                SourceArgvPointer::EmptyLiteral => String::new(),
                SourceArgvPointer::Command(pointer) => {
                    let bytes = heap.slice(pointer.as_const()).unwrap();
                    let length = bytes.iter().position(|byte| *byte == 0).unwrap();
                    String::from_utf8(bytes[..length].iter().map(|byte| *byte as u8).collect())
                        .unwrap()
                }
                SourceArgvPointer::Null => panic!("NULL inside returned argv"),
            })
            .collect();
        let offsets = arguments[..=count as usize]
            .iter()
            .map(|argument| match argument {
                SourceArgvPointer::Null => -1,
                SourceArgvPointer::EmptyLiteral => -2,
                SourceArgvPointer::Command(pointer) => pointer.difference(command).unwrap(),
            })
            .collect();
        let mutated = heap
            .slice(command.as_const())
            .unwrap()
            .iter()
            .map(|byte| *byte as u8)
            .collect();
        inchi_free(&mut heap, command).unwrap();
        (count, values, arguments, offsets, mutated)
    }

    fn source_error_text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> String {
        let bytes = heap.slice(pointer.as_const()).unwrap();
        let length = bytes.iter().position(|byte| *byte == 0).unwrap();
        String::from_utf8(bytes[..length].iter().map(|byte| *byte as u8).collect()).unwrap()
    }

    fn named_atom(name: u8) -> inp_ATOM {
        inp_ATOM {
            elname: [name as i8, 0, 0, 0, 0, 0],
            ..inp_ATOM::default()
        }
    }

    fn source_name(text: &str) -> [i8; 6] {
        assert!(text.len() < 6);
        let mut name = [0_i8; 6];
        for (target, source) in name.iter_mut().zip(text.bytes()) {
            *target = source as i8;
        }
        name
    }

    #[test]
    fn source_port__inchi_dll__extractonestructure__line_1863() {
        let mut heap = SourceHeap::default();
        let mut structure = STRUCT_DATA::default();
        let mut parameters = INPUT_PARMS::default();
        let mut original = ORIG_ATOM_DATA::default();
        let mut number = 0;
        assert_eq!(
            ExtractOneStructure(
                &mut heap,
                &mut structure,
                &mut parameters,
                None,
                None,
                None,
                None,
                None,
                &mut original,
                &mut number,
            ),
            Ok(crate::source_types::_IS_ERROR as i32)
        );
        assert_eq!(structure.nStructReadError, 98);
        assert_eq!(structure.nErrorType, crate::source_types::_IS_ERROR as i32);
        assert_eq!(
            fixed_source_text(&structure.pStrErrStruct).unwrap(),
            "Empty structure"
        );
        assert_eq!(number, 0);

        let atom_pointer = allocate_source_fixture(
            &mut heap,
            vec![inchi_Atom {
                elname: source_name("C"),
                num_iso_H: [-1, 0, 0, 0],
                ..inchi_Atom::default()
            }],
        );
        let too_large = inchi_InputEx {
            atom: atom_pointer,
            num_atoms: NORMALLY_ALLOWED_INP_MAX_ATOMS as i16,
            ..inchi_InputEx::default()
        };
        structure = STRUCT_DATA::default();
        assert_eq!(
            ExtractOneStructure(
                &mut heap,
                &mut structure,
                &mut parameters,
                None,
                Some(&too_large),
                None,
                None,
                None,
                &mut original,
                &mut number,
            ),
            Ok(crate::source_types::_IS_ERROR as i32)
        );
        assert_eq!(structure.nStructReadError, 70);
        assert_eq!(original.num_inp_atoms, -1);

        let valid = inchi_InputEx {
            atom: atom_pointer,
            num_atoms: 1,
            ..inchi_InputEx::default()
        };
        structure = STRUCT_DATA::default();
        parameters.nMode = (REQ_MODE_CHIR_FLG_STEREO
            | REQ_MODE_STEREO
            | REQ_MODE_RELATIVE_STEREO
            | REQ_MODE_RACEMIC_STEREO) as u64;
        parameters.bChiralFlag = FLAG_SET_INP_AT_CHIRAL as i32;
        assert_eq!(
            ExtractOneStructure(
                &mut heap,
                &mut structure,
                &mut parameters,
                None,
                Some(&valid),
                None,
                None,
                None,
                &mut original,
                &mut number,
            ),
            Ok(crate::source_types::_IS_OKAY as i32)
        );
        assert_eq!(number, 1);
        assert_eq!(original.num_inp_atoms, 1);
        assert_eq!(original.num_inp_bonds, 0);
        assert!(!original.at.is_null());
        assert!(!original.szCoord.is_null());
        assert_eq!(heap.slice(original.at.as_const()).unwrap()[0].num_H, 4);
        assert_eq!(structure.bChiralFlag, FLAG_INP_AT_CHIRAL as i32);
        assert_eq!(
            parameters.nMode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO) as u64,
            0
        );

        let mut failed_heap = SourceHeap::default();
        let failed_atom = allocate_source_fixture(
            &mut failed_heap,
            vec![inchi_Atom {
                elname: source_name("C"),
                ..inchi_Atom::default()
            }],
        );
        let failed_input = inchi_InputEx {
            atom: failed_atom,
            num_atoms: 1,
            ..inchi_InputEx::default()
        };
        failed_heap.fail_after_allocations(0);
        let mut failed_structure = STRUCT_DATA::default();
        let mut failed_original = ORIG_ATOM_DATA::default();
        let mut failed_number = 0;
        assert_eq!(
            ExtractOneStructure(
                &mut failed_heap,
                &mut failed_structure,
                &mut INPUT_PARMS::default(),
                None,
                Some(&failed_input),
                None,
                None,
                None,
                &mut failed_original,
                &mut failed_number,
            ),
            Ok(crate::source_types::_IS_FATAL as i32)
        );
        assert_eq!(failed_structure.nStructReadError, -1);
        assert_eq!(
            fixed_source_text(&failed_structure.pStrErrStruct).unwrap(),
            "Out of RAM"
        );
        assert!(failed_original.at.is_null());
        assert!(failed_original.szCoord.is_null());
    }

    #[test]
    fn source_port__inchi_dll__freestructfrominchi__line_208() {
        let mut heap = SourceHeap::default();
        assert_eq!(FreeStructFromINCHI(&mut heap, None), Ok(()));
        let atoms = allocate_source_fixture(&mut heap, vec![inchi_Atom::default()]);
        let stereo = allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default()]);
        let message = allocate_source_fixture(&mut heap, vec![b'M' as i8, 0]);
        let log = allocate_source_fixture(&mut heap, vec![b'L' as i8, 0]);
        let mut output = inchi_OutputStruct {
            atom: atoms,
            stereo0D: stereo,
            num_atoms: 1,
            num_stereo0D: 1,
            szMessage: message,
            szLog: log,
            WarningFlags: [[1, 2], [3, 4]],
        };
        assert_eq!(FreeStructFromINCHI(&mut heap, Some(&mut output)), Ok(()));
        assert_eq!(output, inchi_OutputStruct::default());
        for missing in [
            heap.slice(atoms.as_const()).map(|_| ()),
            heap.slice(stereo.as_const()).map(|_| ()),
            heap.slice(message.as_const()).map(|_| ()),
            heap.slice(log.as_const()).map(|_| ()),
        ] {
            assert_eq!(missing, Err(SourceHeapError::MissingAllocation));
        }

        let mut already_empty = inchi_OutputStruct::default();
        assert_eq!(
            FreeStructFromINCHI(&mut heap, Some(&mut already_empty)),
            Ok(())
        );
    }

    #[test]
    fn source_port__inchi_dll__freestructfromstdinchi__line_195() {
        let mut heap = SourceHeap::default();
        assert_eq!(FreeStructFromStdINCHI(&mut heap, None), Ok(()));
        let atoms = allocate_source_fixture(&mut heap, vec![inchi_Atom::default()]);
        let message = allocate_source_fixture(&mut heap, vec![b'E' as i8, 0]);
        let mut output = inchi_OutputStruct {
            atom: atoms,
            num_atoms: 1,
            szMessage: message,
            WarningFlags: [[u64::MAX; 2]; 2],
            ..inchi_OutputStruct::default()
        };
        assert_eq!(FreeStructFromStdINCHI(&mut heap, Some(&mut output)), Ok(()));
        assert_eq!(output, inchi_OutputStruct::default());
        assert_eq!(
            heap.slice(atoms.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(message.as_const()).map(|_| ()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__inchi_dll__getstringlength__line_2086() {
        let mut heap = SourceHeap::default();
        assert_eq!(GetStringLength(&heap, SourceMutPointer::null()), Ok(0));
        let empty = allocate_source_fixture(&mut heap, vec![0_i8]);
        let embedded =
            allocate_source_fixture(&mut heap, vec![b'a' as i8, b'b' as i8, 0, b'c' as i8, 0]);
        assert_eq!(GetStringLength(&heap, empty), Ok(0));
        assert_eq!(GetStringLength(&heap, embedded), Ok(2));
        let unterminated = allocate_source_fixture(&mut heap, vec![b'x' as i8]);
        assert_eq!(
            GetStringLength(&heap, unterminated),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__inchi_dll__freeinchi__line_153() {
        let mut heap = SourceHeap::default();
        assert_eq!(FreeINCHI(&mut heap, None), Ok(()));

        let inchi = allocate_source_fixture(&mut heap, vec![b'I' as i8, 0]);
        let aux = allocate_source_fixture(&mut heap, vec![b'A' as i8, 0]);
        let message = allocate_source_fixture(&mut heap, vec![b'M' as i8, 0]);
        let log = allocate_source_fixture(&mut heap, vec![b'L' as i8, 0]);
        let mut output = inchi_Output {
            szInChI: inchi,
            szAuxInfo: aux,
            szMessage: message,
            szLog: log,
        };
        assert_eq!(FreeINCHI(&mut heap, Some(&mut output)), Ok(()));
        assert_eq!(output, inchi_Output::default());
        assert!(matches!(
            heap.slice(inchi.as_const()),
            Err(SourceHeapError::MissingAllocation)
        ));
        assert!(matches!(
            heap.slice(log.as_const()),
            Err(SourceHeapError::MissingAllocation)
        ));
        assert!(matches!(
            heap.slice(message.as_const()),
            Err(SourceHeapError::MissingAllocation)
        ));
        assert_eq!(heap.slice(aux.as_const()).unwrap(), &[b'A' as i8, 0]);

        let mut null_output = inchi_Output::default();
        assert_eq!(FreeINCHI(&mut heap, Some(&mut null_output)), Ok(()));
        assert_eq!(null_output, inchi_Output::default());
        inchi_free(&mut heap, aux).unwrap();
    }

    #[test]
    fn source_port__inchi_dll__freestdinchi__line_183() {
        let mut heap = SourceHeap::default();
        assert_eq!(FreeStdINCHI(&mut heap, None), Ok(()));

        let inchi = allocate_source_fixture(&mut heap, vec![b'I' as i8, 0]);
        let aux = allocate_source_fixture(&mut heap, vec![b'A' as i8, 0]);
        let message = allocate_source_fixture(&mut heap, vec![b'M' as i8, 0]);
        let log = allocate_source_fixture(&mut heap, vec![b'L' as i8, 0]);
        let mut output = inchi_Output {
            szInChI: inchi,
            szAuxInfo: aux,
            szMessage: message,
            szLog: log,
        };
        assert_eq!(FreeStdINCHI(&mut heap, Some(&mut output)), Ok(()));
        assert_eq!(output, inchi_Output::default());
        for released in [inchi, message, log] {
            assert!(matches!(
                heap.slice(released.as_const()),
                Err(SourceHeapError::MissingAllocation)
            ));
        }
        assert_eq!(heap.slice(aux.as_const()).unwrap(), &[b'A' as i8, 0]);
        inchi_free(&mut heap, aux).unwrap();
    }

    fn source_stereo_chain(length: usize, first_parity: i8, second_parity: i8) -> Vec<inp_ATOM> {
        assert!((2..=12).contains(&length));
        let mut atoms = vec![inp_ATOM::default(); length];
        for (index, atom) in atoms.iter_mut().enumerate() {
            atom.elname = source_name("C");
            if index == 0 {
                atom.valence = 1;
                atom.neighbor[0] = 1;
            } else if index + 1 == length {
                atom.valence = 1;
                atom.neighbor[0] = (index - 1) as u16;
            } else {
                atom.valence = 2;
                atom.neighbor[0] = (index - 1) as u16;
                atom.neighbor[1] = (index + 1) as u16;
            }
        }
        atoms[0].sb_parity[0] = first_parity;
        atoms[0].sb_ord[0] = 0;
        atoms[0].sn_orig_at_num[0] = 10;
        atoms[length - 1].sb_parity[0] = second_parity;
        atoms[length - 1].sb_ord[0] = 0;
        atoms[length - 1].sn_orig_at_num[0] = 20;
        atoms
    }

    #[test]
    fn source_port__inchi_dll__inpatom0dtoinchiatom__line_1670() {
        {
            let mut heap = SourceHeap::default();
            let old_atoms = allocate_source_fixture(&mut heap, vec![inchi_Atom::default()]);
            let old_stereo = allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default()]);
            let mut atom_count = 7;
            let mut atoms = old_atoms;
            let mut stereo_count = 9;
            let mut stereo = old_stereo;
            assert_eq!(
                InpAtom0DToInchiAtom(
                    &mut heap,
                    SourceMutPointer::null(),
                    0,
                    Some(&mut atom_count),
                    Some(&mut atoms),
                    Some(&mut stereo_count),
                    Some(&mut stereo),
                ),
                Ok(0)
            );
            assert_eq!((atom_count, atoms.is_null()), (0, true));
            assert_eq!((stereo_count, stereo.is_null()), (0, true));
            assert!(heap.slice(old_atoms.as_const()).is_ok());
            assert!(heap.slice(old_stereo.as_const()).is_ok());
            inchi_free(&mut heap, old_atoms).unwrap();
            inchi_free(&mut heap, old_stereo).unwrap();
        }

        {
            let mut heap = SourceHeap::default();
            let mut atom_count = 7;
            let mut atoms = SourceMutPointer::null();
            let mut stereo_count = 9;
            let mut stereo = SourceMutPointer::null();
            assert_eq!(
                InpAtom0DToInchiAtom(
                    &mut heap,
                    SourceMutPointer::null(),
                    -1,
                    Some(&mut atom_count),
                    Some(&mut atoms),
                    Some(&mut stereo_count),
                    Some(&mut stereo),
                ),
                Ok(-1)
            );
            assert_eq!((atom_count, atoms.is_null()), (0, true));
            assert_eq!((stereo_count, stereo.is_null()), (0, true));
        }

        {
            let mut heap = SourceHeap::default();
            let mut input = inp_ATOM {
                elname: [b'C' as i8, 0, b'X' as i8, b'Y' as i8, b'Z' as i8, 0],
                valence: 2,
                charge: -2,
                iso_atw_diff: 2,
                num_H: 3,
                num_iso_H: [1, 2, 3],
                radical: 2,
                p_parity: 2,
                p_orig_at_num: [1, 2, 3, 4],
                ..inp_ATOM::default()
            };
            input.bond_type[..2].copy_from_slice(&[1, 4]);
            input.neighbor[..2].copy_from_slice(&[1, u16::MAX]);
            let input_pointer = allocate_source_fixture(&mut heap, vec![input.clone()]);
            let mut atom_count = 0;
            let mut atoms = SourceMutPointer::null();
            let mut stereo_count = 0;
            let mut stereo = SourceMutPointer::null();
            assert_eq!(
                InpAtom0DToInchiAtom(
                    &mut heap,
                    input_pointer,
                    1,
                    Some(&mut atom_count),
                    Some(&mut atoms),
                    Some(&mut stereo_count),
                    Some(&mut stereo),
                ),
                Ok(0)
            );
            assert_eq!((atom_count, stereo_count), (1, 1));
            let output = &heap.slice(atoms.as_const()).unwrap()[0];
            assert_eq!(output.num_bonds, 2);
            assert_eq!(&output.bond_type[..2], &[1, 4]);
            assert_eq!(&output.neighbor[..2], &[1, -1]);
            assert_eq!(output.charge, -2);
            assert_eq!(output.elname, input.elname);
            assert_eq!(output.isotopic_mass, 10_001);
            assert_eq!(output.num_iso_H, [3, 1, 2, 3]);
            assert_eq!(output.radical, 2);
            assert_eq!((output.x, output.y, output.z), (0.0, 0.0, 0.0));
            let stereo_output = &heap.slice(stereo.as_const()).unwrap()[0];
            assert_eq!(stereo_output.central_atom, 0);
            assert_eq!(stereo_output.parity, 2);
            assert_eq!(
                stereo_output.type_,
                tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8
            );
            assert_eq!(stereo_output.neighbor, [0, 1, 2, 3]);
            inchi_free(&mut heap, atoms).unwrap();
            inchi_free(&mut heap, stereo).unwrap();
        }

        for (length, first, second, expected_parity, expected_type, expected_center) in [
            (
                2_usize,
                1_i8,
                1_i8,
                2_i8,
                tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
                -1_i16,
            ),
            (
                2,
                1,
                2,
                1,
                tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
                -1,
            ),
            (
                3,
                3,
                4,
                4,
                tagINCHIStereoType0D_INCHI_StereoType_Allene as i8,
                1,
            ),
            (
                4,
                2,
                2,
                2,
                tagINCHIStereoType0D_INCHI_StereoType_DoubleBond as i8,
                -1,
            ),
        ] {
            let mut heap = SourceHeap::default();
            let input = source_stereo_chain(length, first, second);
            let input_pointer = allocate_source_fixture(&mut heap, input);
            let mut atom_count = 0;
            let mut atoms = SourceMutPointer::null();
            let mut stereo_count = 0;
            let mut stereo = SourceMutPointer::null();
            assert_eq!(
                InpAtom0DToInchiAtom(
                    &mut heap,
                    input_pointer,
                    length as i32,
                    Some(&mut atom_count),
                    Some(&mut atoms),
                    Some(&mut stereo_count),
                    Some(&mut stereo),
                ),
                Ok(0)
            );
            assert_eq!((atom_count, stereo_count), (length as i16, 1));
            let output = &heap.slice(stereo.as_const()).unwrap()[0];
            assert_eq!(output.parity, expected_parity);
            assert_eq!(output.type_, expected_type);
            assert_eq!(output.central_atom, expected_center);
            assert_eq!(output.neighbor, [9, 0, length as i16 - 1, 19]);
            inchi_free(&mut heap, atoms).unwrap();
            inchi_free(&mut heap, stereo).unwrap();
        }

        {
            let mut heap = SourceHeap::default();
            let input = source_stereo_chain(5, 1, 2);
            let input_pointer = allocate_source_fixture(&mut heap, input);
            let mut atom_count = 0;
            let mut atoms = SourceMutPointer::null();
            let mut stereo_count = -1;
            let mut stereo = SourceMutPointer::null();
            assert_eq!(
                InpAtom0DToInchiAtom(
                    &mut heap,
                    input_pointer,
                    5,
                    Some(&mut atom_count),
                    Some(&mut atoms),
                    Some(&mut stereo_count),
                    Some(&mut stereo),
                ),
                Ok(0)
            );
            assert_eq!((atom_count, stereo_count), (5, 0));
            assert!(!stereo.is_null());
            assert_eq!(
                heap.slice(stereo.as_const()).unwrap()[0],
                inchi_Stereo0D::default()
            );
            inchi_free(&mut heap, atoms).unwrap();
            inchi_free(&mut heap, stereo).unwrap();
        }

        {
            let mut heap = SourceHeap::default();
            let mut input = source_stereo_chain(2, 1, 1);
            input[0].sb_parity[1] = 2;
            input[0].sb_ord[1] = 0;
            input[0].sn_orig_at_num[1] = 11;
            let input_pointer = allocate_source_fixture(&mut heap, input);
            let mut atom_count = 0;
            let mut atoms = SourceMutPointer::null();
            let mut stereo_count = 0;
            let mut stereo = SourceMutPointer::null();
            assert_eq!(
                InpAtom0DToInchiAtom(
                    &mut heap,
                    input_pointer,
                    2,
                    Some(&mut atom_count),
                    Some(&mut atoms),
                    Some(&mut stereo_count),
                    Some(&mut stereo),
                ),
                Ok(1)
            );
            assert_eq!((atom_count, stereo_count), (2, 1));
            inchi_free(&mut heap, atoms).unwrap();
            inchi_free(&mut heap, stereo).unwrap();
        }

        for successful_allocations in [0_u64, 1] {
            let mut heap = SourceHeap::default();
            let input_pointer = allocate_source_fixture(
                &mut heap,
                vec![inp_ATOM {
                    p_parity: 1,
                    ..inp_ATOM::default()
                }],
            );
            heap.fail_after_allocations(successful_allocations);
            let mut atom_count = 5;
            let mut atoms = SourceMutPointer::null();
            let mut stereo_count = 6;
            let mut stereo = SourceMutPointer::null();
            assert_eq!(
                InpAtom0DToInchiAtom(
                    &mut heap,
                    input_pointer,
                    1,
                    Some(&mut atom_count),
                    Some(&mut atoms),
                    Some(&mut stereo_count),
                    Some(&mut stereo),
                ),
                Ok(-1)
            );
            assert_eq!(heap.source_allocation_calls(), 2);
            assert_eq!((atom_count, atoms.is_null()), (0, true));
            assert_eq!((stereo_count, stereo.is_null()), (0, true));
        }

        let mut heap = SourceHeap::default();
        assert_eq!(
            InpAtom0DToInchiAtom(
                &mut heap,
                SourceMutPointer::null(),
                0,
                None,
                None,
                None,
                None,
            ),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__inchi_dll__setatomandbondproperties__line_1439() {
        for (bond_types, expected_valence, expected_error, expected_text) in [
            (vec![1_u8, 2, 3], 6_i8, 0, ""),
            (vec![4_u8, 4], 3, 0, ""),
            (vec![4_u8, 4, 4], 4, 0, ""),
            (vec![4_u8], 0, 8, "Atom 'C' has 1 alternating bonds"),
            (
                vec![4_u8, 4, 4, 4],
                0,
                8,
                "Atom 'C' has 4 alternating bonds",
            ),
        ] {
            let mut heap = SourceHeap::default();
            let mut atom = named_atom(b'C');
            atom.valence = bond_types.len() as i8;
            atom.bond_type[..bond_types.len()].copy_from_slice(&bond_types);
            atom.neighbor[..bond_types.len()].fill(1);
            let mut metal = named_atom(b'F');
            metal.el_number = 26;
            let atoms = allocate_source_fixture(&mut heap, vec![atom, metal]);
            let input = inchi_Atom {
                elname: source_name("C"),
                ..inchi_Atom::default()
            };
            let inputs = allocate_source_fixture(&mut heap, vec![input, inchi_Atom::default()]);
            let errors = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
            let mut error = 0;
            assert_eq!(
                SetAtomAndBondProperties(
                    &mut heap,
                    atoms,
                    inputs,
                    0,
                    0,
                    Some(errors),
                    Some(&mut error),
                ),
                Ok(0)
            );
            let atom = &heap.slice(atoms.as_const()).unwrap()[0];
            assert_eq!(atom.chem_bonds_valence, expected_valence);
            assert_eq!(error, expected_error);
            assert_eq!(source_error_text(&heap, errors), expected_text);
        }

        for (
            alias,
            initial_charge,
            initial_radical,
            expected_name,
            expected_number,
            expected_h,
            expected_isotope,
            expected_charge,
            expected_radical,
            expected_type,
            expected_error,
            expected_text,
        ) in [
            (
                "CH3+",
                0,
                0,
                "C",
                6,
                3,
                0,
                1,
                0,
                3,
                0,
                "Parsed compound atom(s): CH3+",
            ),
            (
                "CH3+",
                2,
                0,
                "C",
                6,
                3,
                0,
                2,
                0,
                3,
                0,
                "Ignored charge/radical redefinition: CH3+; Parsed compound atom(s):",
            ),
            (
                "N^H",
                0,
                3,
                "N",
                7,
                1,
                0,
                0,
                3,
                3,
                0,
                "Ignored charge/radical redefinition: N^H; Parsed compound atom(s):",
            ),
            (
                "H2",
                0,
                0,
                "H",
                1,
                1,
                0,
                0,
                0,
                3,
                0,
                "Parsed compound atom(s): H2",
            ),
            (
                "D2",
                0,
                0,
                "H",
                1,
                0,
                2,
                0,
                0,
                3,
                0,
                "Parsed compound atom(s): D2",
            ),
            (
                "QH2",
                0,
                0,
                "Q",
                0,
                2,
                0,
                0,
                0,
                2,
                64,
                "Unknown element(s): Q",
            ),
        ] {
            let mut heap = SourceHeap::default();
            let mut atom = inp_ATOM {
                elname: source_name(alias),
                charge: initial_charge,
                radical: initial_radical,
                ..inp_ATOM::default()
            };
            atom.elname.copy_from_slice(&source_name(alias));
            let atoms = allocate_source_fixture(&mut heap, vec![atom]);
            let input = inchi_Atom {
                elname: source_name(alias),
                num_iso_H: [-1, 0, 0, 0],
                ..inchi_Atom::default()
            };
            let inputs = allocate_source_fixture(&mut heap, vec![input]);
            let errors = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
            let mut error = 0;
            SetAtomAndBondProperties(
                &mut heap,
                atoms,
                inputs,
                0,
                0,
                Some(errors),
                Some(&mut error),
            )
            .unwrap();
            let atom = &heap.slice(atoms.as_const()).unwrap()[0];
            assert_eq!(
                fixed_source_text(&atom.elname).unwrap(),
                expected_name,
                "{alias}"
            );
            assert_eq!(atom.el_number, expected_number, "{alias}");
            assert_eq!(atom.num_H, expected_h, "{alias}");
            assert_eq!(atom.iso_atw_diff, expected_isotope, "{alias}");
            assert_eq!(atom.charge, expected_charge, "{alias}");
            assert_eq!(atom.radical, expected_radical, "{alias}");
            assert_eq!(atom.at_type, expected_type, "{alias}");
            assert_eq!(error, expected_error, "{alias}");
            assert_eq!(source_error_text(&heap, errors), expected_text, "{alias}");
        }

        for (element, initial_isotope, isotopic_mass, valence, expected_name, expected_isotope) in [
            ("D", 0_i8, 0_i16, 0_i8, "H", 2_i8),
            ("T", 0, 0, 0, "H", 3),
            ("H", 0, 10_000, 0, "H", 1),
            ("H", 0, 9_999, 0, "H", 0),
            ("H", 0, 2, 0, "H", 2),
            ("H", 0, 5, 1, "H", 0),
            ("H", 0, 5, 2, "H", 5),
            ("H", 7, 2, 0, "H", 2),
            ("C", 0, 13, 0, "C", 2),
            ("C", 0, 9_999, 0, "C", -1),
        ] {
            let mut heap = SourceHeap::default();
            let mut atom = inp_ATOM {
                elname: source_name(element),
                iso_atw_diff: initial_isotope,
                valence,
                ..inp_ATOM::default()
            };
            if element == "D" || element == "T" {
                atom.elname[2..].fill(b'X' as i8);
            }
            let atoms = allocate_source_fixture(&mut heap, vec![atom]);
            let input = inchi_Atom {
                elname: source_name(element),
                isotopic_mass,
                ..inchi_Atom::default()
            };
            let inputs = allocate_source_fixture(&mut heap, vec![input]);
            let mut error = 0;
            SetAtomAndBondProperties(&mut heap, atoms, inputs, 0, 0, None, Some(&mut error))
                .unwrap();
            let atom = &heap.slice(atoms.as_const()).unwrap()[0];
            assert_eq!(
                fixed_source_text(&atom.elname).unwrap(),
                expected_name,
                "{element}"
            );
            assert_eq!(
                atom.iso_atw_diff, expected_isotope,
                "{element}/{isotopic_mass}"
            );
            if element == "D" || element == "T" {
                assert_eq!(atom.elname, source_name("H"));
            }
            assert_eq!(error, 0);
        }

        for (
            input_hydrogens,
            do_not_add,
            initial_h,
            expected_h,
            expected_isotopes,
            expected_type,
        ) in [
            ([-1_i8, 1, 2, 3], 0, 7, 7, [1, 2, 3], 2_u16),
            ([-1, 1, 2, 3], 1, 7, 7, [1, 2, 3], 0),
            ([4, -1, -2, -3], 0, 7, 4, [-1, -2, -3], 0),
        ] {
            let mut heap = SourceHeap::default();
            let atoms = allocate_source_fixture(
                &mut heap,
                vec![inp_ATOM {
                    elname: source_name("C"),
                    num_H: initial_h,
                    ..inp_ATOM::default()
                }],
            );
            let inputs = allocate_source_fixture(
                &mut heap,
                vec![inchi_Atom {
                    elname: source_name("C"),
                    num_iso_H: input_hydrogens,
                    ..inchi_Atom::default()
                }],
            );
            SetAtomAndBondProperties(&mut heap, atoms, inputs, 0, do_not_add, None, None).unwrap();
            let atom = &heap.slice(atoms.as_const()).unwrap()[0];
            assert_eq!(atom.num_H, expected_h);
            assert_eq!(atom.num_iso_H, expected_isotopes);
            assert_eq!(atom.at_type, expected_type);
        }

        for metal_neighbor in [false, true] {
            let mut heap = SourceHeap::default();
            let mut nitrogen = named_atom(b'N');
            nitrogen.valence = 2;
            nitrogen.charge = 1;
            nitrogen.bond_type[..2].fill(BOND_TYPE_ALTERN as u8);
            nitrogen.neighbor.copy_from_slice(&[1_u16; MAXVAL as usize]);
            let mut neighbor = named_atom(if metal_neighbor { b'F' } else { b'C' });
            neighbor.el_number = if metal_neighbor { 26 } else { 6 };
            let atoms = allocate_source_fixture(&mut heap, vec![nitrogen, neighbor]);
            let inputs = allocate_source_fixture(
                &mut heap,
                vec![
                    inchi_Atom {
                        elname: source_name("N"),
                        num_iso_H: [1, 0, 0, 0],
                        ..inchi_Atom::default()
                    },
                    inchi_Atom::default(),
                ],
            );
            SetAtomAndBondProperties(&mut heap, atoms, inputs, 0, 0, None, None).unwrap();
            assert_eq!(
                heap.slice(atoms.as_const()).unwrap()[0].chem_bonds_valence,
                if metal_neighbor { 3 } else { 2 }
            );
        }

        let mut heap = SourceHeap::default();
        assert_eq!(
            SetAtomAndBondProperties(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                -1,
                0,
                None,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let atoms = allocate_source_fixture(
            &mut heap,
            vec![inp_ATOM {
                elname: [b'C' as i8; 6],
                ..inp_ATOM::default()
            }],
        );
        let inputs = allocate_source_fixture(&mut heap, vec![inchi_Atom::default()]);
        assert_eq!(
            SetAtomAndBondProperties(&mut heap, atoms, inputs, 0, 0, None, None),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__inchi_dll__parse_options_string__line_1037() {
        let (count, values, arguments, offsets, mutated) = parse_case(b"\0", 4);
        assert_eq!((count, values), (1, vec![String::new()]));
        assert_eq!(arguments[0], SourceArgvPointer::EmptyLiteral);
        assert_eq!(arguments[1], SourceArgvPointer::Null);
        assert_eq!(offsets, [-2, -1]);
        assert_eq!(mutated, b"\0");

        let (count, values, arguments, offsets, mutated) =
            parse_case(b"  alpha\t\"beta gamma\" delta\0", 8);
        assert_eq!(count, 4);
        assert_eq!(values, ["", "alpha", "beta gamma", "delta"]);
        assert_eq!(arguments[4], SourceArgvPointer::Null);
        assert_eq!(offsets, [-2, 2, 8, 21, -1]);
        assert_eq!(mutated, b"  alpha\0beta gamma\0\" delta\0");

        let (_, values, _, _, _) = parse_case(b"ab\"cd\"ef \"a\"\"b x\\\\\\\"y z\\\\\"q\0", 8);
        assert_eq!(values, ["", "abcdef", "a\"b", "x\\\"y", "z\\q"]);

        let (count, values, arguments, offsets, mutated) = parse_case(b"one two three\0", 3);
        assert_eq!((count, values), (2, vec![String::new(), "one".to_owned()]));
        assert_eq!(arguments[2], SourceArgvPointer::Null);
        assert_eq!(offsets, [-2, 0, -1]);
        assert_eq!(mutated, b"one\0two three\0");

        for (text, max_arguments, expected) in [
            (&b"  \t \t\0"[..], 8, vec![""]),
            (&b"\"\" tail\0"[..], 8, vec!["", "", "tail"]),
            (&b"ab\"cd ef\"gh\0"[..], 8, vec!["", "abcd efgh"]),
            (&b"\"a\"\"b\" \"\"\"\"\0"[..], 8, vec!["", "a\"b \""]),
            (&b"a\\\"b c\0"[..], 8, vec!["", "a\"b", "c"]),
            (&b"a\\\\\"b c\"\0"[..], 8, vec!["", "a\\b c"]),
            (&b"a\\\\\\\"b c\0"[..], 8, vec!["", "a\\\"b", "c"]),
            (&b"a\\\\\\\\\"b c\"\0"[..], 8, vec!["", "a\\\\b c"]),
            (&b"alpha\\\0"[..], 8, vec!["", "alpha\\"]),
            (&b"\"alpha beta\0"[..], 8, vec!["", "alpha beta"]),
        ] {
            let (_, values, arguments, _, _) = parse_case(text, max_arguments);
            assert_eq!(values, expected, "{text:?}");
            assert_eq!(arguments[expected.len()], SourceArgvPointer::Null);
        }

        let (count, values, arguments, offsets, mutated) = parse_case(
            b"00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 \
20 21 22 23 24 25 26 27 28 29 30 31 32 33 34\0",
            32,
        );
        assert_eq!(count, 31);
        assert_eq!(values.len(), 31);
        assert_eq!(values[1], "00");
        assert_eq!(values[30], "29");
        assert_eq!(arguments[31], SourceArgvPointer::Null);
        assert_eq!(offsets[0], -2);
        assert_eq!(offsets[1..31], (0_i64..=87).step_by(3).collect::<Vec<_>>());
        assert_eq!(offsets[31], -1);
        assert_eq!(&mutated[87..], b"29\030 31 32 33 34\0");

        assert_eq!(
            parse_options_string(
                &mut SourceHeap::default(),
                SourceMutPointer::null(),
                &mut [SourceArgvPointer::Null; 1],
                1,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn official_c_oracle__inchi_dll__parse_options_string__line_1037() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--parse-options-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );

        let records = String::from_utf8(oracle.stdout)
            .expect("official C oracle output must be UTF-8")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("oracle record must be JSON"))
            .collect::<Vec<_>>();
        assert_eq!(records.len(), 18);

        for official in records {
            assert_eq!(official["operation"], "parse_options_string");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let max_arguments = official["input"]["maxargs"]
                .as_i64()
                .and_then(|value| i32::try_from(value).ok())
                .expect("maxargs must be i32");
            let mut input = official["input"]["bytes"]
                .as_array()
                .expect("input bytes must be an array")
                .iter()
                .map(|value| value.as_u64().expect("input byte must be u8") as u8)
                .collect::<Vec<_>>();
            input.extend([0xa5; 16]);

            let mut heap = SourceHeap::default();
            let command = allocate_source_fixture(
                &mut heap,
                input.iter().map(|byte| *byte as i8).collect::<Vec<_>>(),
            );
            let mut arguments = vec![SourceArgvPointer::Null; max_arguments as usize];
            let count = parse_options_string(&mut heap, command, &mut arguments, max_arguments)
                .unwrap_or_else(|error| panic!("{case_id}: {error:?}"));

            let rust_offsets = arguments[..=count as usize]
                .iter()
                .map(|argument| match argument {
                    SourceArgvPointer::Null => -1,
                    SourceArgvPointer::EmptyLiteral => -2,
                    SourceArgvPointer::Command(pointer) => pointer.difference(command).unwrap(),
                })
                .collect::<Vec<_>>();
            let rust_argv_bytes = arguments[..count as usize]
                .iter()
                .map(|argument| match argument {
                    SourceArgvPointer::EmptyLiteral => Vec::new(),
                    SourceArgvPointer::Command(pointer) => {
                        let bytes = heap.slice(pointer.as_const()).unwrap();
                        let length = bytes.iter().position(|byte| *byte == 0).unwrap();
                        bytes[..length].iter().map(|byte| *byte as u8).collect()
                    }
                    SourceArgvPointer::Null => panic!("{case_id}: NULL inside argv"),
                })
                .collect::<Vec<Vec<u8>>>();
            let rust_buffer = heap
                .slice(command.as_const())
                .unwrap()
                .iter()
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>();

            assert_eq!(official["output"]["count"], count, "{case_id}: count");
            assert_eq!(
                serde_json::to_value(rust_offsets).unwrap(),
                official["output"]["argv_offsets"],
                "{case_id}: argv offsets"
            );
            assert_eq!(
                serde_json::to_value(rust_argv_bytes).unwrap(),
                official["output"]["argv_bytes"],
                "{case_id}: argv bytes"
            );
            assert_eq!(
                serde_json::to_value(rust_buffer).unwrap(),
                official["output"]["buffer_bytes"],
                "{case_id}: mutated buffer"
            );
            assert_eq!(official["output"]["terminal_null"], true);
            assert_eq!(arguments[count as usize], SourceArgvPointer::Null);
        }
    }

    #[test]
    fn source_port__inchi_dll__setatomproperties__line_1139() {
        let mut heap = SourceHeap::default();
        let mut target = inp_ATOM {
            elname: [
                b'X' as i8, b'Y' as i8, b'Z' as i8, b'W' as i8, b'Q' as i8, b'R' as i8,
            ],
            ..inp_ATOM::default()
        };
        let input = inchi_Atom {
            elname: [b'C' as i8, 0, b'A' as i8, b'B' as i8, b'C' as i8, 0],
            charge: -2,
            radical: tagINCHIRadical_INCHI_RADICAL_SINGLET as i8,
            x: 1.25,
            y: -2.5,
            z: 3.0e-6,
            ..inchi_Atom::default()
        };
        let targets = allocate_source_fixture(&mut heap, vec![target.clone(), target.clone()]);
        let inputs = allocate_source_fixture(&mut heap, vec![inchi_Atom::default(), input]);
        let coordinates = allocate_source_fixture(&mut heap, vec![[b'#' as i8; 32]; 2]);
        let mut dimensions = 4;
        let mut error = 0;
        assert_eq!(
            SetAtomProperties(
                &mut heap,
                targets,
                coordinates,
                inputs,
                1,
                Some(&mut dimensions),
                None,
                Some(&mut error),
            ),
            Ok(0)
        );
        target = heap.slice(targets.as_const()).unwrap()[1].clone();
        assert_eq!(
            target.elname,
            [
                b'C' as i8, 0, b'Z' as i8, b'W' as i8, b'Q' as i8, b'R' as i8
            ]
        );
        assert_eq!((target.charge, target.radical), (-2, RADICAL_TRIPLET as i8));
        assert_eq!((target.x, target.y, target.z), (1.25, -2.5, 3.0e-6));
        assert_eq!(target.orig_at_number, 2);
        assert_eq!((dimensions, error), (7, 0));
        let coordinate = &heap.slice(coordinates.as_const()).unwrap()[1];
        assert_eq!(
            &coordinate[..30],
            b"    1.2500   -2.5000    0.0000"
                .iter()
                .map(|byte| *byte as i8)
                .collect::<Vec<_>>()
        );
        assert_eq!(&coordinate[30..], &[b'#' as i8; 2]);

        for (source, expected) in [(0_i8, 0_i8), (1, 3), (2, 2), (3, 3)] {
            let known_input = allocate_source_fixture(
                &mut heap,
                vec![inchi_Atom {
                    elname: [b'N' as i8, 0, 0, 0, 0, 0],
                    radical: source,
                    ..inchi_Atom::default()
                }],
            );
            let known_target = allocate_source_fixture(&mut heap, vec![inp_ATOM::default()]);
            assert_eq!(
                SetAtomProperties(
                    &mut heap,
                    known_target,
                    SourceMutPointer::null(),
                    known_input,
                    0,
                    None,
                    None,
                    None,
                ),
                Ok(0)
            );
            assert_eq!(
                heap.slice(known_target.as_const()).unwrap()[0].radical,
                expected
            );
            inchi_free(&mut heap, known_input).unwrap();
            inchi_free(&mut heap, known_target).unwrap();
        }

        for (source, initial_error, expected_radical, expected_error, expected_text) in [
            (8_i8, 4, 2_i8, 4, "Radical center type replaced: 8->2"),
            (-1_i8, 1, -1_i8, 9, "Radical center type replaced: -1->-1"),
        ] {
            let unknown_input = allocate_source_fixture(
                &mut heap,
                vec![inchi_Atom {
                    elname: [b'B' as i8, 0, 0, 0, 0, 0],
                    radical: source,
                    ..inchi_Atom::default()
                }],
            );
            let unknown_target = allocate_source_fixture(&mut heap, vec![inp_ATOM::default()]);
            let error_buffer = allocate_source_fixture(&mut heap, vec![0_i8; 256]);
            let mut case_error = initial_error;
            assert_eq!(
                SetAtomProperties(
                    &mut heap,
                    unknown_target,
                    SourceMutPointer::null(),
                    unknown_input,
                    0,
                    None,
                    Some(error_buffer),
                    Some(&mut case_error),
                ),
                Ok(0)
            );
            assert_eq!(
                heap.slice(unknown_target.as_const()).unwrap()[0].radical,
                expected_radical
            );
            assert_eq!(case_error, expected_error);
            let bytes = heap.slice(error_buffer.as_const()).unwrap();
            let length = bytes.iter().position(|byte| *byte == 0).unwrap();
            assert_eq!(
                String::from_utf8(bytes[..length].iter().map(|byte| *byte as u8).collect())
                    .unwrap(),
                expected_text
            );
            inchi_free(&mut heap, unknown_input).unwrap();
            inchi_free(&mut heap, unknown_target).unwrap();
            inchi_free(&mut heap, error_buffer).unwrap();
        }

        let threshold_inputs = allocate_source_fixture(
            &mut heap,
            vec![
                inchi_Atom {
                    elname: [b'H' as i8, 0, 0, 0, 0, 0],
                    x: 1.0e-6,
                    ..inchi_Atom::default()
                },
                inchi_Atom {
                    elname: [b'H' as i8, 0, 0, 0, 0, 0],
                    y: 1.000_001e-6,
                    ..inchi_Atom::default()
                },
            ],
        );
        let threshold_targets = allocate_source_fixture(&mut heap, vec![inp_ATOM::default(); 2]);
        let mut threshold_dimensions = 4;
        for index in 0..2 {
            SetAtomProperties(
                &mut heap,
                threshold_targets,
                SourceMutPointer::null(),
                threshold_inputs,
                index,
                Some(&mut threshold_dimensions),
                None,
                None,
            )
            .unwrap();
            assert_eq!(threshold_dimensions, if index == 0 { 4 } else { 6 });
        }

        inchi_free(&mut heap, targets).unwrap();
        inchi_free(&mut heap, inputs).unwrap();
        inchi_free(&mut heap, coordinates).unwrap();
        inchi_free(&mut heap, threshold_inputs).unwrap();
        inchi_free(&mut heap, threshold_targets).unwrap();
    }

    #[test]
    fn official_c_oracle__inchi_dll__setatomproperties__line_1139() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        fn atom_integer_fields(atom: &inp_ATOM) -> Vec<i64> {
            let mut fields = Vec::new();
            fields.extend(atom.elname.iter().map(|value| i64::from(*value)));
            fields.push(i64::from(atom.el_number));
            fields.extend(atom.neighbor.iter().map(|value| i64::from(*value)));
            fields.push(i64::from(atom.orig_at_number));
            fields.push(i64::from(atom.orig_compt_at_numb));
            fields.extend(atom.bond_stereo.iter().map(|value| i64::from(*value)));
            fields.extend(atom.bond_type.iter().map(|value| i64::from(*value)));
            fields.push(i64::from(atom.valence));
            fields.push(i64::from(atom.chem_bonds_valence));
            fields.push(i64::from(atom.num_H));
            fields.extend(atom.num_iso_H.iter().map(|value| i64::from(*value)));
            fields.push(i64::from(atom.iso_atw_diff));
            fields.push(i64::from(atom.charge));
            fields.push(i64::from(atom.radical));
            fields.push(i64::from(atom.bAmbiguousStereo));
            fields.push(i64::from(atom.cFlags));
            fields.push(i64::from(atom.at_type));
            fields.push(i64::from(atom.component));
            fields.push(i64::from(atom.endpoint));
            fields.push(i64::from(atom.c_point));
            fields.push(i64::from(atom.bUsed0DParity));
            fields.push(i64::from(atom.p_parity));
            fields.extend(atom.p_orig_at_num.iter().map(|value| i64::from(*value)));
            fields.extend(atom.sb_ord.iter().map(|value| i64::from(*value)));
            fields.extend(atom.sn_ord.iter().map(|value| i64::from(*value)));
            fields.extend(atom.sb_parity.iter().map(|value| i64::from(*value)));
            fields.extend(atom.sn_orig_at_num.iter().map(|value| i64::from(*value)));
            fields.push(i64::from(atom.bCutVertex));
            fields.push(i64::from(atom.nRingSystem));
            fields.push(i64::from(atom.nNumAtInRingSystem));
            fields.push(i64::from(atom.nBlockSystem));
            fields
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--set-atom-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["operation"], "SetAtomProperties");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let input = &official["input"];
            let atom_index = input["atom_index"].as_i64().unwrap() as usize;
            let mut input_atoms = vec![inchi_Atom::default(); 2];
            let element = input["element"].as_str().unwrap().as_bytes();
            for (target, source) in input_atoms[atom_index].elname.iter_mut().zip(element) {
                *target = *source as i8;
            }
            input_atoms[atom_index].radical = input["radical"].as_i64().unwrap() as i8;
            input_atoms[atom_index].charge = input["charge"].as_i64().unwrap() as i8;
            let coordinate_bits = input["coordinate_bits"].as_array().unwrap();
            input_atoms[atom_index].x = f64::from_bits(coordinate_bits[0].as_u64().unwrap());
            input_atoms[atom_index].y = f64::from_bits(coordinate_bits[1].as_u64().unwrap());
            input_atoms[atom_index].z = f64::from_bits(coordinate_bits[2].as_u64().unwrap());
            let mut target_atoms = vec![inp_ATOM::default(); 2];
            target_atoms[atom_index].elname.fill(b'Z' as i8);
            let mut heap = SourceHeap::default();
            let targets = allocate_source_fixture(&mut heap, target_atoms);
            let inputs = allocate_source_fixture(&mut heap, input_atoms);
            let coordinates = allocate_source_fixture(&mut heap, vec![[b'#' as i8; 32]; 2]);
            let mut error_bytes = vec![0_i8; STR_ERR_LEN as usize];
            let initial_error_text = input["initial_error_text"].as_str().unwrap().as_bytes();
            for (target, source) in error_bytes.iter_mut().zip(initial_error_text) {
                *target = *source as i8;
            }
            let error_buffer = allocate_source_fixture(&mut heap, error_bytes);
            let mut dimensions = input["initial_dimensions"].as_i64().unwrap() as i32;
            let mut error = input["initial_error"].as_i64().unwrap() as i32;
            let use_coordinates = input["use_coordinates"].as_bool().unwrap();
            let status = SetAtomProperties(
                &mut heap,
                targets,
                if use_coordinates {
                    coordinates
                } else {
                    SourceMutPointer::null()
                },
                inputs,
                atom_index as i32,
                Some(&mut dimensions),
                Some(error_buffer),
                Some(&mut error),
            )
            .unwrap_or_else(|failure| panic!("{case_id}: {failure:?}"));
            let expected = &official["output"];
            assert_eq!(
                status,
                expected["status"].as_i64().unwrap() as i32,
                "{case_id}"
            );
            let atom = &heap.slice(targets.as_const()).unwrap()[atom_index];
            let expected_fields = expected["atom"]["integer_fields"]
                .as_array()
                .unwrap()
                .iter()
                .map(|value| value.as_i64().unwrap())
                .collect::<Vec<_>>();
            assert_eq!(atom_integer_fields(atom), expected_fields, "{case_id}");
            assert_eq!(
                [atom.x.to_bits(), atom.y.to_bits(), atom.z.to_bits()],
                expected["atom"]["coordinate_bits"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|value| value.as_u64().unwrap())
                    .collect::<Vec<_>>()
                    .as_slice(),
                "{case_id}"
            );
            let coordinate = &heap.slice(coordinates.as_const()).unwrap()[atom_index];
            let expected_coordinate = expected["coordinate_buffer"]
                .as_array()
                .unwrap()
                .iter()
                .map(|value| value.as_i64().unwrap() as i8)
                .collect::<Vec<_>>();
            assert_eq!(coordinate.as_slice(), expected_coordinate, "{case_id}");
            assert_eq!(
                dimensions,
                expected["dimensions"].as_i64().unwrap() as i32,
                "{case_id}"
            );
            assert_eq!(
                error,
                expected["error"].as_i64().unwrap() as i32,
                "{case_id}"
            );
            assert_eq!(
                source_error_text(&heap, error_buffer),
                expected["error_text"],
                "{case_id}"
            );
            record_count += 1;
        }
        assert_eq!(record_count, 264);
    }

    #[test]
    fn source_port__inchi_dll__setbondproperties__line_1235() {
        let mapping_cases = [
            (1_i8, 0_i8, 1_u8, 0_i8, 0_i8),
            (2, 1, 2, 1, -1),
            (3, 4, 3, 4, -4),
            (4, 6, 4, 6, -6),
            (1, -1, 1, -1, 1),
            (2, -4, 2, -4, 4),
            (3, -6, 3, -6, 6),
            (4, 3, 4, 3, 3),
            (1, -3, 1, 3, 3),
        ];
        for (source_type, source_stereo, expected_type, expected1, expected2) in mapping_cases {
            let mut heap = SourceHeap::default();
            let atoms =
                allocate_source_fixture(&mut heap, vec![named_atom(b'C'), named_atom(b'O')]);
            let mut input = inchi_Atom::default();
            input.neighbor[0] = 1;
            input.bond_type[0] = source_type;
            input.bond_stereo[0] = source_stereo;
            let inputs = allocate_source_fixture(&mut heap, vec![input, inchi_Atom::default()]);
            let error_buffer = allocate_source_fixture(&mut heap, vec![0_i8; 256]);
            let mut num_bonds = 0;
            let mut error = 0;
            assert_eq!(
                SetBondProperties(
                    &mut heap,
                    atoms,
                    inputs,
                    0,
                    0,
                    2,
                    Some(&mut num_bonds),
                    Some(error_buffer),
                    Some(&mut error),
                ),
                Ok(0)
            );
            let result = heap.slice(atoms.as_const()).unwrap();
            assert_eq!((result[0].valence, result[1].valence), (1, 1));
            assert_eq!((result[0].neighbor[0], result[1].neighbor[0]), (1, 0));
            assert_eq!(
                (result[0].bond_type[0], result[1].bond_type[0]),
                (expected_type, expected_type)
            );
            assert_eq!(
                (result[0].bond_stereo[0], result[1].bond_stereo[0]),
                (expected1, expected2)
            );
            assert_eq!(
                (num_bonds, error, source_error_text(&heap, error_buffer)),
                (1, 0, String::new())
            );
        }

        {
            let mut heap = SourceHeap::default();
            let atoms =
                allocate_source_fixture(&mut heap, vec![named_atom(b'C'), named_atom(b'N')]);
            let mut input = inchi_Atom::default();
            input.neighbor[0] = 1;
            input.bond_type[0] = 9;
            input.bond_stereo[0] = 9;
            let inputs = allocate_source_fixture(&mut heap, vec![input, inchi_Atom::default()]);
            let error_buffer = allocate_source_fixture(&mut heap, vec![0_i8; 256]);
            let mut num_bonds = 5;
            let mut error = 16;
            assert_eq!(
                SetBondProperties(
                    &mut heap,
                    atoms,
                    inputs,
                    0,
                    0,
                    2,
                    Some(&mut num_bonds),
                    Some(error_buffer),
                    Some(&mut error),
                ),
                Ok(0)
            );
            let result = heap.slice(atoms.as_const()).unwrap();
            assert_eq!((result[0].bond_type[0], result[0].bond_stereo[0]), (1, 0));
            assert_eq!((num_bonds, error), (6, 24));
            assert_eq!(
                source_error_text(&heap, error_buffer),
                "Unrecognized bond type: 9; Unrecognized bond stereo:"
            );
        }

        for (neighbor, atom_index, expected_text) in [
            (-1_i16, 0_i32, "Bond to nonexistent atom"),
            (2, 0, "Bond to nonexistent atom"),
            (0, 0, "Atom has a bond to itself"),
        ] {
            let mut heap = SourceHeap::default();
            let atoms =
                allocate_source_fixture(&mut heap, vec![named_atom(b'C'), named_atom(b'O')]);
            let mut input = inchi_Atom::default();
            input.neighbor[0] = neighbor;
            input.bond_type[0] = 1;
            let inputs = allocate_source_fixture(&mut heap, vec![input, inchi_Atom::default()]);
            let error_buffer = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
            let mut num_bonds = 3;
            let mut error = 8;
            assert_eq!(
                SetBondProperties(
                    &mut heap,
                    atoms,
                    inputs,
                    atom_index,
                    0,
                    2,
                    Some(&mut num_bonds),
                    Some(error_buffer),
                    Some(&mut error),
                ),
                Ok(1)
            );
            assert_eq!((num_bonds, error), (3, 9));
            assert_eq!(source_error_text(&heap, error_buffer), expected_text);
            assert_eq!(heap.slice(atoms.as_const()).unwrap()[0].valence, 0);
        }

        {
            let mut heap = SourceHeap::default();
            let mut atom1 = named_atom(b'C');
            atom1.valence = 1;
            atom1.neighbor[0] = 1;
            atom1.bond_type[0] = 2;
            atom1.bond_stereo[0] = 1;
            let mut atom2 = named_atom(b'O');
            atom2.valence = 1;
            atom2.neighbor[0] = 0;
            atom2.bond_type[0] = 2;
            atom2.bond_stereo[0] = -1;
            let atoms = allocate_source_fixture(&mut heap, vec![atom1, atom2]);
            let mut input = inchi_Atom::default();
            input.neighbor[0] = 1;
            input.bond_type[0] = 2;
            input.bond_stereo[0] = 1;
            let inputs = allocate_source_fixture(&mut heap, vec![input, inchi_Atom::default()]);
            let mut num_bonds = 7;
            let mut error = 0;
            assert_eq!(
                SetBondProperties(
                    &mut heap,
                    atoms,
                    inputs,
                    0,
                    0,
                    2,
                    Some(&mut num_bonds),
                    None,
                    Some(&mut error),
                ),
                Ok(0)
            );
            assert_eq!((num_bonds, error), (7, 0));

            heap.slice_mut(atoms).unwrap()[0].bond_type[0] = 1;
            let error_buffer = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
            assert_eq!(
                SetBondProperties(
                    &mut heap,
                    atoms,
                    inputs,
                    0,
                    0,
                    2,
                    Some(&mut num_bonds),
                    Some(error_buffer),
                    Some(&mut error),
                ),
                Ok(0)
            );
            assert_eq!(error, 2);
            assert_eq!(
                source_error_text(&heap, error_buffer),
                "Multiple bonds between two atoms"
            );
            assert_eq!(heap.slice(atoms.as_const()).unwrap()[0].bond_type[0], 2);
        }

        for mismatch in [false, true] {
            let mut heap = SourceHeap::default();
            let mut atom1 = named_atom(b'C');
            atom1.valence = 1;
            atom1.neighbor[0] = 1;
            atom1.bond_type[0] = 1;
            atom1.bond_stereo[0] = 1;
            let atoms = allocate_source_fixture(&mut heap, vec![atom1, named_atom(b'O')]);
            let mut input = inchi_Atom::default();
            input.neighbor[0] = 1;
            input.bond_type[0] = if mismatch { 2 } else { 1 };
            input.bond_stereo[0] = 1;
            let inputs = allocate_source_fixture(&mut heap, vec![input, inchi_Atom::default()]);
            let error_buffer = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
            let mut num_bonds = 9;
            let mut error = 0;
            assert_eq!(
                SetBondProperties(
                    &mut heap,
                    atoms,
                    inputs,
                    0,
                    0,
                    2,
                    Some(&mut num_bonds),
                    Some(error_buffer),
                    Some(&mut error),
                ),
                Ok(0)
            );
            let result = heap.slice(atoms.as_const()).unwrap();
            assert_eq!((result[0].valence, result[1].valence), (1, 1));
            assert_eq!((result[1].neighbor[0], result[1].bond_stereo[0]), (0, -1));
            assert_eq!(num_bonds, 9);
            if mismatch {
                assert_eq!(
                    (error, result[0].bond_type[0], result[1].bond_type[0]),
                    (2, 2, 2)
                );
                assert_eq!(
                    source_error_text(&heap, error_buffer),
                    "Multiple bonds between two atoms"
                );
            } else {
                assert_eq!(
                    (error, result[0].bond_type[0], result[1].bond_type[0]),
                    (0, 1, 1)
                );
                assert_eq!(
                    source_error_text(&heap, error_buffer),
                    "Duplicated bond(s) between two atoms"
                );
            }
        }

        {
            let mut heap = SourceHeap::default();
            let mut atom1 = named_atom(b'C');
            atom1.valence = 2;
            atom1.neighbor[..2].copy_from_slice(&[1, 1]);
            atom1.bond_type[..2].copy_from_slice(&[1, 1]);
            let mut atom2 = named_atom(b'O');
            atom2.valence = 1;
            atom2.neighbor[0] = 0;
            atom2.bond_type[0] = 1;
            let atoms = allocate_source_fixture(&mut heap, vec![atom1, atom2]);
            let mut input = inchi_Atom::default();
            input.neighbor[0] = 1;
            input.bond_type[0] = 1;
            let inputs = allocate_source_fixture(&mut heap, vec![input, inchi_Atom::default()]);
            let error_buffer = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
            let mut num_bonds = 1;
            let mut error = 0;
            assert_eq!(
                SetBondProperties(
                    &mut heap,
                    atoms,
                    inputs,
                    0,
                    0,
                    2,
                    Some(&mut num_bonds),
                    Some(error_buffer),
                    Some(&mut error),
                ),
                Ok(0)
            );
            assert_eq!((num_bonds, error), (1, 2));
            assert_eq!(
                source_error_text(&heap, error_buffer),
                "Multiple bonds between two atoms"
            );
        }

        for one_sided_full in [false, true] {
            let mut heap = SourceHeap::default();
            let mut atom1 = named_atom(b'C');
            let mut atom2 = named_atom(b'O');
            if one_sided_full {
                atom1.valence = 1;
                atom1.neighbor[0] = 1;
                atom1.bond_type[0] = 1;
                atom2.valence = MAXVAL as i8;
                atom2.neighbor.fill(2);
            } else {
                atom1.valence = MAXVAL as i8;
                atom1.neighbor.fill(0);
            }
            let atoms = allocate_source_fixture(&mut heap, vec![atom1, atom2]);
            let mut input = inchi_Atom::default();
            input.neighbor[0] = 1;
            input.bond_type[0] = 1;
            let inputs = allocate_source_fixture(&mut heap, vec![input, inchi_Atom::default()]);
            let error_buffer = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
            let mut num_bonds = 4;
            let mut error = 0;
            assert_eq!(
                SetBondProperties(
                    &mut heap,
                    atoms,
                    inputs,
                    0,
                    0,
                    2,
                    Some(&mut num_bonds),
                    Some(error_buffer),
                    Some(&mut error),
                ),
                Ok(1)
            );
            assert_eq!((num_bonds, error), (4, 4));
            assert_eq!(
                source_error_text(&heap, error_buffer),
                if one_sided_full {
                    "Atom 'O' has more than 20 bonds"
                } else {
                    "Atom 'C' has more than 20 bonds"
                }
            );
        }

        let mut heap = SourceHeap::default();
        assert_eq!(
            SetBondProperties(
                &mut heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                -1,
                0,
                0,
                None,
                None,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn official_c_oracle__inchi_dll__setbondproperties__line_1235() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        fn atom_fields(atom: &inp_ATOM) -> Vec<i64> {
            let mut out = Vec::new();
            out.extend(atom.elname.iter().map(|v| i64::from(*v)));
            out.push(i64::from(atom.el_number));
            out.extend(atom.neighbor.iter().map(|v| i64::from(*v)));
            out.extend([
                i64::from(atom.orig_at_number),
                i64::from(atom.orig_compt_at_numb),
            ]);
            out.extend(atom.bond_stereo.iter().map(|v| i64::from(*v)));
            out.extend(atom.bond_type.iter().map(|v| i64::from(*v)));
            out.extend([
                i64::from(atom.valence),
                i64::from(atom.chem_bonds_valence),
                i64::from(atom.num_H),
            ]);
            out.extend(atom.num_iso_H.iter().map(|v| i64::from(*v)));
            out.extend([
                i64::from(atom.iso_atw_diff),
                i64::from(atom.charge),
                i64::from(atom.radical),
                i64::from(atom.bAmbiguousStereo),
                i64::from(atom.cFlags),
                i64::from(atom.at_type),
                i64::from(atom.component),
                i64::from(atom.endpoint),
                i64::from(atom.c_point),
                i64::from(atom.bUsed0DParity),
                i64::from(atom.p_parity),
            ]);
            out.extend(atom.p_orig_at_num.iter().map(|v| i64::from(*v)));
            out.extend(atom.sb_ord.iter().map(|v| i64::from(*v)));
            out.extend(atom.sn_ord.iter().map(|v| i64::from(*v)));
            out.extend(atom.sb_parity.iter().map(|v| i64::from(*v)));
            out.extend(atom.sn_orig_at_num.iter().map(|v| i64::from(*v)));
            out.extend([
                i64::from(atom.bCutVertex),
                i64::from(atom.nRingSystem),
                i64::from(atom.nNumAtInRingSystem),
                i64::from(atom.nBlockSystem),
            ]);
            out
        }

        fn setup_atoms(mode: i32) -> Vec<inp_ATOM> {
            let mut atoms = vec![named_atom(b'C'), named_atom(b'O')];
            match mode {
                4 | 5 | 15 => {
                    atoms[0].valence = 1;
                    atoms[1].valence = 1;
                    atoms[0].neighbor[0] = 1;
                    atoms[1].neighbor[0] = 0;
                    atoms[0].bond_type[0] = if mode == 5 { 2 } else { 1 };
                    atoms[1].bond_type[0] = if mode == 5 { 2 } else { 1 };
                    atoms[0].bond_stereo[0] = if mode == 15 { 1 } else { 0 };
                    atoms[1].bond_stereo[0] = if mode == 15 { -1 } else { 0 };
                }
                6 | 13 => {
                    atoms[0].valence = 2;
                    atoms[0].neighbor[..2].fill(1);
                    atoms[0].bond_type[..2].fill(1);
                    atoms[1].valence = 1;
                    atoms[1].neighbor[0] = 0;
                    atoms[1].bond_type[0] = 1;
                }
                7 | 8 | 12 | 14 => {
                    atoms[0].valence = 1;
                    atoms[0].neighbor[0] = 1;
                    atoms[0].bond_type[0] = if mode == 8 || mode == 14 { 2 } else { 1 };
                    if mode == 12 {
                        atoms[1].valence = MAXVAL as i8;
                        atoms[1].neighbor.fill(99);
                    }
                }
                9 => {
                    atoms[1].valence = 1;
                    atoms[1].neighbor[0] = 0;
                    atoms[1].bond_type[0] = 1;
                }
                10 | 11 => {
                    let full = if mode == 10 {
                        &mut atoms[0]
                    } else {
                        &mut atoms[1]
                    };
                    full.valence = MAXVAL as i8;
                    full.neighbor.fill(99);
                }
                _ => {}
            }
            atoms
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--set-bond-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut count = 0;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            let case_id = official["case_id"].as_str().unwrap();
            let input = &official["input"];
            let mode = input["mode"].as_i64().unwrap() as i32;
            let atoms = setup_atoms(mode);
            let mut source = vec![inchi_Atom::default(); 2];
            source[0].neighbor[0] = if mode == 1 {
                -1
            } else if mode == 2 {
                2
            } else if mode == 3 {
                0
            } else {
                1
            };
            source[0].bond_type[0] = input["bond_type"].as_i64().unwrap() as i8;
            source[0].bond_stereo[0] = input["bond_stereo"].as_i64().unwrap() as i8;
            let mut heap = SourceHeap::default();
            let atom_pointer = allocate_source_fixture(&mut heap, atoms);
            let input_pointer = allocate_source_fixture(&mut heap, source);
            let error_buffer = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
            let mut num_bonds = 5;
            let mut error = 16;
            let status = SetBondProperties(
                &mut heap,
                atom_pointer,
                input_pointer,
                0,
                0,
                2,
                Some(&mut num_bonds),
                Some(error_buffer),
                Some(&mut error),
            )
            .unwrap_or_else(|failure| panic!("{case_id}: {failure:?}"));
            let expected = &official["output"];
            assert_eq!(
                status,
                expected["status"].as_i64().unwrap() as i32,
                "{case_id}"
            );
            let actual_atoms = heap.slice(atom_pointer.as_const()).unwrap();
            for index in 0..2 {
                let expected_atom = &expected["atoms"][index];
                let expected_fields = expected_atom["integer_fields"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_i64().unwrap())
                    .collect::<Vec<_>>();
                assert_eq!(
                    atom_fields(&actual_atoms[index]),
                    expected_fields,
                    "{case_id} atom {index}"
                );
                assert_eq!(
                    [
                        actual_atoms[index].x.to_bits(),
                        actual_atoms[index].y.to_bits(),
                        actual_atoms[index].z.to_bits()
                    ],
                    expected_atom["coordinate_bits"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|v| v.as_u64().unwrap())
                        .collect::<Vec<_>>()
                        .as_slice(),
                    "{case_id} atom {index}"
                );
            }
            assert_eq!(
                num_bonds,
                expected["num_bonds"].as_i64().unwrap() as i32,
                "{case_id}"
            );
            assert_eq!(
                error,
                expected["error"].as_i64().unwrap() as i32,
                "{case_id}"
            );
            assert_eq!(
                source_error_text(&heap, error_buffer),
                expected["error_text"].as_str().unwrap(),
                "{case_id}"
            );
            count += 1;
        }
        assert_eq!(count, 527);
    }
}
