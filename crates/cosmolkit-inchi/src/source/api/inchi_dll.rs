use crate::source::api::inchi_dll_b::add_source_error;
use crate::source::api::inchi_dll_ext::{
    SetExtOrigAtDataByInChIExtInput, SetInChIExtInputByExtOrigAtData,
};
use crate::source::base::ichi_io::{
    inchi_ios_close, inchi_ios_eprint, inchi_ios_init, inchi_ios_print_nodisplay, inchi_ios_reset,
    inchi_strbuf_close, inchi_strbuf_init,
};
use crate::source::base::ichican2::SetBitFree;
use crate::source::base::ichimak2::WriteCoord;
use crate::source::base::ichiparm::{
    HelpCommandLineParms, InchiBuildMetadata, PrintInputParms, ReadCommandLineParms,
};
use crate::source::base::ichiread::ReadWriteInChI;
use crate::source::base::mol2atom::{FreeExtOrigAtData, FreeOrigAtData};
use crate::source::base::readinch::Extract0DParities;
use crate::source::base::runichi::ProcessOneStructureEx;
use crate::source::base::runichi2::TreatErrorsInReadTheStructure;
use crate::source::base::runichi3::OAD_Polymer_SmartReopenCyclizedUnits;
use crate::source::base::runichi4::FreeAllINChIArrays;
use crate::source::base::util::{
    detect_unusual_el_valence, extract_charges_and_radicals, extract_h_atoms,
    extract_inchi_substring, get_atomic_mass_from_elnum, get_num_H, get_periodic_table_number,
    inchi_calloc, inchi_free, inchi_malloc, inchi_stricmp, is_el_a_metal, is_in_the_list,
    mystrncpy_slice, n_bonds_val_to_metal,
};
use crate::source_types::{
    _IS_EOF, _IS_ERROR, _IS_FATAL, AB_PARITY_UNDF, AB_PARITY_UNKN, BOND_TYPE_ALTERN,
    BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE, BOND_TYPE_TRIPLE, CANON_GLOBALS, FILE, FLAG_INP_AT_CHIRAL,
    FLAG_INP_AT_NONCHIRAL, FLAG_SET_INP_AT_CHIRAL, INCHI_CLOCK, INCHI_IOS_TYPE_STRING,
    INCHI_IOSTREAM, INCHI_MAX_NUM_ARG, INCHI_OUT_SDFILE_ONLY, INPUT_PARMS, ISOTOPIC_SHIFT_FLAG,
    ISOTOPIC_SHIFT_MAX, LOG_MASK_NO_WARN, MAX_ATOMS, MAX_NUM_PATHS, MAX_NUM_STEREO_ATOM_NEIGH,
    MAX_NUM_STEREO_BONDS, MAX_SDF_VALUE, MAXVAL, MOL_COORD, NO_ATOM,
    NORMALLY_ALLOWED_INP_MAX_ATOMS, NUM_H_ISOTOPES, OAD_Polymer, OAD_V3000, ORIG_ATOM_DATA,
    RADICAL_DOUBLET, RADICAL_TRIPLET, REQ_MODE_CHIR_FLG_STEREO, REQ_MODE_DIFF_UU_STEREO,
    REQ_MODE_RACEMIC_STEREO, REQ_MODE_RELATIVE_STEREO, REQ_MODE_STEREO, STEREO_DBLE_EITHER,
    STEREO_SNGL_DOWN, STEREO_SNGL_EITHER, STEREO_SNGL_UP, STRUCT_DATA, SourceArgvPointer,
    SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList, clock_t,
    inchi_Atom, inchi_Input, inchi_Input_Polymer, inchi_Input_V3000, inchi_InputEx,
    inchi_InputINCHI, inchi_Output, inchi_OutputStruct, inchi_OutputStructEx, inchi_Stereo0D,
    inp_ATOM, local_util::ERR_ELEM, tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER,
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
    tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral, tagInputType_INPUT_INCHI,
    tagRetValCheckINCHI_INCHI_FAIL_I2I, tagRetValCheckINCHI_INCHI_INVALID_LAYOUT,
    tagRetValCheckINCHI_INCHI_INVALID_PREFIX, tagRetValCheckINCHI_INCHI_INVALID_VERSION,
    tagRetValCheckINCHI_INCHI_VALID_BETA, tagRetValCheckINCHI_INCHI_VALID_NON_STANDARD,
    tagRetValCheckINCHI_INCHI_VALID_STANDARD, tagRetValGetINCHI_inchi_Ret_ERROR,
    tagRetValGetINCHI_inchi_Ret_FATAL,
};
use crate::source_types::{
    _IS_OKAY, _IS_SKIP, _IS_UNKNOWN, _IS_WARNING, INCHI_IOS_STRING, INCHI_OUT_SAVEOPT,
    INCHI_OUT_STDINCHI, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT, MAX_SDF_HEADER,
    REQ_MODE_BASIC, REQ_MODE_SB_IGN_ALL_UU, REQ_MODE_SC_IGN_ALL_UU, TG_FLAG_1_5_TAUT,
    TG_FLAG_KETO_ENOL_TAUT, TG_FLAG_RECONNECT_COORD, tagRetValGetINCHI_inchi_Ret_EOF,
    tagRetValGetINCHI_inchi_Ret_OKAY, tagRetValGetINCHI_inchi_Ret_SKIP,
    tagRetValGetINCHI_inchi_Ret_UNKNOWN, tagRetValGetINCHI_inchi_Ret_WARNING,
};

const SOURCE_SIZEOF_INCHI_ATOM: u64 = 3 * 8 + 20 * 2 + 20 + 20 + 6 + 2 + 4 + 2 + 1 + 1;
const SOURCE_SIZEOF_INCHI_STEREO0D: u64 = 5 * 2 + 2;

pub(crate) fn input_erroneously_contains_pseudoatoms(
    heap: &mut SourceHeap,
    input: Option<&inchi_Input>,
    output: Option<&mut inchi_Output>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:293 input_erroneously_contains_pseudoatoms
    // INCHI✔️❌: int input_erroneously_contains_pseudoatoms( inchi_Input *inp,
    // INCHI✔️❌:                                             inchi_Output *out)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char *str_noz = "Unsupported in this mode element \'*\'";
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     /* Supposed that no '*' or 'Zz' elements are allowed in the input. */
    // INCHI✔️❌:     for (i = 0; i < inp->num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!strcmp(inp->atom->elname, "Zz") || !strcmp(inp->atom->elname, "*"))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (out)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 memset(out, 0, sizeof(*out)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:                 if ((out->szMessage = (char *)inchi_malloc(strlen(str_noz) + 1))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     strcpy(out->szMessage, str_noz);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: input_erroneously_contains_pseudoatoms

    const MESSAGE: &[u8] = b"Unsupported in this mode element '*'\0";

    let input = input.ok_or(SourceHeapError::NullPointer)?;
    let mut output = output;
    for _ in 0..input.num_atoms {
        let first_atom = heap
            .slice(input.atom.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let element_length = first_atom
            .elname
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let element = &first_atom.elname[..element_length];
        if element == [b'Z' as i8, b'z' as i8] || element == [b'*' as i8] {
            if let Some(output) = output.as_deref_mut() {
                *output = inchi_Output::default();
                output.szMessage = match inchi_malloc(heap, MESSAGE.len() as u64) {
                    Ok(message) => {
                        for (target, source) in heap.slice_mut(message)?.iter_mut().zip(MESSAGE) {
                            *target = *source as i8;
                        }
                        message
                    }
                    Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                    Err(error) => return Err(error),
                };
            }
            return Ok(1);
        }
    }

    Ok(0)
}

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

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetStdINCHI(
    heap: &mut SourceHeap,
    input: Option<&inchi_Input>,
    mut output: Option<&mut inchi_Output>,
    interrupted: i32,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:242 GetStdINCHI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL GetStdINCHI( inchi_Input *inp, inchi_Output *out )
{
    inchi_InputEx extended_input;

    /* No '*' or 'Zz' elements are allowed in the input . */
    if (input_erroneously_contains_pseudoatoms(inp, out))
    {
        return _IS_ERROR;
    }

    extended_input.atom = inp->atom;
    extended_input.num_atoms = inp->num_atoms;
    extended_input.num_stereo0D = inp->num_stereo0D;
    extended_input.stereo0D = inp->stereo0D;
    extended_input.szOptions = inp->szOptions;
    extended_input.polymer = NULL;
    extended_input.v3000 = NULL;

    return GetINCHI1( &extended_input, out, 1 );
}
    */
    // END INCHI C FUNCTION: GetStdINCHI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetStdINCHI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local conditional compilation.
    // INCHI✔️❌: bInterrupted is modeled explicitly for the active GetINCHI1 callee; NULL inp is outside the C-defined domain.
    // INCHI✔️❌: GetINCHI1 uses SourceHeap records, map access, and cloning that are materially slower than native pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetStdINCHI

    if input_erroneously_contains_pseudoatoms(heap, input, output.as_deref_mut())? != 0 {
        return Ok(_IS_ERROR as i32);
    }
    let input = input.ok_or(SourceHeapError::NullPointer)?;
    let extended_input = inchi_InputEx {
        atom: input.atom,
        stereo0D: input.stereo0D,
        szOptions: input.szOptions,
        num_atoms: input.num_atoms,
        num_stereo0D: input.num_stereo0D,
        polymer: SourceMutPointer::null(),
        v3000: SourceMutPointer::null(),
    };
    GetINCHI1(
        heap,
        Some(&extended_input),
        output,
        1,
        interrupted,
        stdout,
        build,
        clock_result,
    )
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetINCHI(
    heap: &mut SourceHeap,
    input: Option<&inchi_Input>,
    mut output: Option<&mut inchi_Output>,
    interrupted: i32,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:270 GetINCHI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL GetINCHI( inchi_Input *inp, inchi_Output *out )
{
    inchi_InputEx extended_input;

    /* For back compatibility: no '*' or 'Zz' elements are allowed in the input to GetINCHI() ! */
    if ( input_erroneously_contains_pseudoatoms( inp, out) )
    {
        return _IS_ERROR;
    }

    extended_input.atom = inp->atom;
    extended_input.num_atoms = inp->num_atoms;
    extended_input.stereo0D = inp->stereo0D;
    extended_input.num_stereo0D = inp->num_stereo0D;
    extended_input.szOptions = inp->szOptions;
    extended_input.polymer = NULL;
    extended_input.v3000 = NULL;

    return GetINCHI1( &extended_input, out, 0 );
}
    */
    // END INCHI C FUNCTION: GetINCHI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetINCHI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local conditional compilation.
    // INCHI✔️❌: bInterrupted is modeled explicitly for the active GetINCHI1 callee; NULL inp is outside the C-defined domain.
    // INCHI✔️❌: GetINCHI1 uses SourceHeap records, map access, and cloning that are materially slower than native pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetINCHI

    if input_erroneously_contains_pseudoatoms(heap, input, output.as_deref_mut())? != 0 {
        return Ok(_IS_ERROR as i32);
    }
    let input = input.ok_or(SourceHeapError::NullPointer)?;
    let extended_input = inchi_InputEx {
        atom: input.atom,
        stereo0D: input.stereo0D,
        szOptions: input.szOptions,
        num_atoms: input.num_atoms,
        num_stereo0D: input.num_stereo0D,
        polymer: SourceMutPointer::null(),
        v3000: SourceMutPointer::null(),
    };
    GetINCHI1(
        heap,
        Some(&extended_input),
        output,
        0,
        interrupted,
        stdout,
        build,
        clock_result,
    )
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetINCHIEx(
    heap: &mut SourceHeap,
    input: Option<&mut inchi_InputEx>,
    output: Option<&mut inchi_Output>,
    interrupted: i32,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:325 GetINCHIEx
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL GetINCHIEx( inchi_InputEx *inp, inchi_Output *out )
{
    int i;

    /* Check for star atoms and replace them by Zz atoms */
    for (i = 0; i < inp->num_atoms; i++)
    {
        if (!strcmp( inp->atom[i].elname, "*" ))
        {
            strcpy( inp->atom[i].elname, "Zz" );
        }
    }

    return GetINCHI1( inp, out, 0 );
}
    */
    // END INCHI C FUNCTION: GetINCHIEx
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetINCHIEx
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local conditional compilation.
    // INCHI✔️❌: bInterrupted is modeled explicitly for the active GetINCHI1 callee; NULL inp is outside the C-defined domain.
    // INCHI✔️❌: GetINCHI1 uses SourceHeap records, map access, and cloning that are materially slower than native pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetINCHIEx

    let input = input.ok_or(SourceHeapError::NullPointer)?;
    for index in 0..input.num_atoms {
        let atom = heap
            .slice(input.atom.as_const())?
            .get(index as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let name_length = atom
            .elname
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if &atom.elname[..name_length] == [b'*' as i8] {
            let atom = heap
                .slice_mut(input.atom)?
                .get_mut(index as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            atom.elname[0] = b'Z' as i8;
            atom.elname[1] = b'z' as i8;
            atom.elname[2] = 0;
        }
    }

    GetINCHI1(
        heap,
        Some(input),
        output,
        0,
        interrupted,
        stdout,
        build,
        clock_result,
    )
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn FreeInChIExtInput(
    heap: &mut SourceHeap,
    polymer: SourceMutPointer<inchi_Input_Polymer>,
    v3000: SourceMutPointer<inchi_Input_V3000>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2923 FreeInChIExtInput
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void FreeInChIExtInput( inchi_Input_Polymer *polymer, inchi_Input_V3000 *v3000 )
{
    int k;
    if (polymer)
    {
        if (polymer->n && polymer->units)
        {
            for (k = 0; k < polymer->n; k++)
            {
                if (polymer->units[k])
                {
                    if (polymer->units[k]->alist)
                    {
                        inchi_free( polymer->units[k]->alist );  polymer->units[k]->alist = NULL;
                    }
                    if (polymer->units[k]->blist)
                    {
                        inchi_free( polymer->units[k]->blist );  polymer->units[k]->blist = NULL;
                    }
                }
                inchi_free( polymer->units[k] );
            }
            inchi_free( polymer->units );
            polymer->units = NULL;
            inchi_free( polymer );
        }
    }
    if (v3000)
    {
        if (v3000->atom_index_orig)
        {
            inchi_free( v3000->atom_index_orig );
            v3000->atom_index_orig = NULL;
        }
        if (v3000->atom_index_fin)
        {
            inchi_free( v3000->atom_index_fin );
            v3000->atom_index_fin = NULL;
        }
        if (v3000->n_haptic_bonds && v3000->lists_haptic_bonds)
        {
            for (k = 0; k < v3000->n_haptic_bonds; k++)
            {
                if (v3000->lists_haptic_bonds[k])
                {
                    inchi_free( v3000->lists_haptic_bonds[k] );
                    v3000->lists_haptic_bonds[k] = NULL;
                }
            }
            inchi_free( v3000->lists_haptic_bonds );
            v3000->lists_haptic_bonds = NULL;
        }
        if (v3000->n_steabs && v3000->lists_steabs)
        {
            for (k = 0; k < v3000->n_steabs; k++)
            {
                if (v3000->lists_steabs[k])
                {
                    inchi_free( v3000->lists_steabs[k] );
                    v3000->lists_steabs[k] = NULL;
                }
            }
            inchi_free( v3000->lists_steabs );
            v3000->lists_steabs = NULL;
        }
        if (v3000->n_sterel && v3000->lists_sterel)
        {
            for (k = 0; k < v3000->n_sterel; k++)
            {
                if (v3000->lists_sterel[k])
                {
                    inchi_free( v3000->lists_sterel[k] );
                    v3000->lists_sterel[k] = NULL;
                }
            }
            inchi_free( v3000->lists_sterel );
            v3000->lists_sterel = NULL;
        }
        if (v3000->n_sterac && v3000->lists_sterac)
        {
            for (k = 0; k < v3000->n_sterac; k++)
            {
                if (v3000->lists_sterac[k])
                {
                    inchi_free( v3000->lists_sterac[k] );
                    v3000->lists_sterac[k] = NULL;
                }
            }
            inchi_free( v3000->lists_sterac );
            v3000->lists_sterac = NULL;
        }
        inchi_free( v3000 );
        /*memset( v3000, 0, sizeof( *v3000 ) );*/
    }
}
    */
    // END INCHI C FUNCTION: FreeInChIExtInput
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FreeInChIExtInput
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; inchi_free is the active libc free macro.
    // INCHI✔️❌: Zero counts skip nested-list and polymer-container frees; negative nonzero counts free outer lists after zero loop iterations.
    // INCHI✔️❌: SourceHeap pointer validation and map-based allocation ownership are materially slower than native pointer frees.
    // END INCHI ACTIVE MACRO CONFIGURATION: FreeInChIExtInput

    if !polymer.is_null() {
        let (count, units) = {
            let value = heap
                .slice(polymer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (value.n, value.units)
        };
        if count != 0 && !units.is_null() {
            for index in 0..count {
                let unit = heap.slice(units.as_const())?[index as usize];
                if !unit.is_null() {
                    let (atom_list, bond_list) = {
                        let value = heap
                            .slice(unit.as_const())?
                            .first()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        (value.alist, value.blist)
                    };
                    if !atom_list.is_null() {
                        inchi_free(heap, atom_list)?;
                        heap.slice_mut(unit)?
                            .first_mut()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .alist = SourceMutPointer::null();
                    }
                    if !bond_list.is_null() {
                        inchi_free(heap, bond_list)?;
                        heap.slice_mut(unit)?
                            .first_mut()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .blist = SourceMutPointer::null();
                    }
                }
                inchi_free(heap, unit)?;
            }
            inchi_free(heap, units)?;
            heap.slice_mut(polymer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .units = SourceMutPointer::null();
            inchi_free(heap, polymer)?;
        }
    }

    if !v3000.is_null() {
        let (atom_index_orig, atom_index_fin) = {
            let value = heap
                .slice(v3000.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (value.atom_index_orig, value.atom_index_fin)
        };
        if !atom_index_orig.is_null() {
            inchi_free(heap, atom_index_orig)?;
            heap.slice_mut(v3000)?[0].atom_index_orig = SourceMutPointer::null();
        }
        if !atom_index_fin.is_null() {
            inchi_free(heap, atom_index_fin)?;
            heap.slice_mut(v3000)?[0].atom_index_fin = SourceMutPointer::null();
        }

        let (count, list) = {
            let value = &heap.slice(v3000.as_const())?[0];
            (value.n_haptic_bonds, value.lists_haptic_bonds)
        };
        if count != 0 && !list.is_null() {
            for index in 0..count {
                let row = heap.slice(list.as_const())?[index as usize];
                if !row.is_null() {
                    inchi_free(heap, row)?;
                    heap.slice_mut(list)?[index as usize] = SourceMutPointer::null();
                }
            }
            inchi_free(heap, list)?;
            heap.slice_mut(v3000)?[0].lists_haptic_bonds = SourceMutPointer::null();
        }

        let (count, list) = {
            let value = &heap.slice(v3000.as_const())?[0];
            (value.n_steabs, value.lists_steabs)
        };
        if count != 0 && !list.is_null() {
            for index in 0..count {
                let row = heap.slice(list.as_const())?[index as usize];
                if !row.is_null() {
                    inchi_free(heap, row)?;
                    heap.slice_mut(list)?[index as usize] = SourceMutPointer::null();
                }
            }
            inchi_free(heap, list)?;
            heap.slice_mut(v3000)?[0].lists_steabs = SourceMutPointer::null();
        }

        let (count, list) = {
            let value = &heap.slice(v3000.as_const())?[0];
            (value.n_sterel, value.lists_sterel)
        };
        if count != 0 && !list.is_null() {
            for index in 0..count {
                let row = heap.slice(list.as_const())?[index as usize];
                if !row.is_null() {
                    inchi_free(heap, row)?;
                    heap.slice_mut(list)?[index as usize] = SourceMutPointer::null();
                }
            }
            inchi_free(heap, list)?;
            heap.slice_mut(v3000)?[0].lists_sterel = SourceMutPointer::null();
        }

        let (count, list) = {
            let value = &heap.slice(v3000.as_const())?[0];
            (value.n_sterac, value.lists_sterac)
        };
        if count != 0 && !list.is_null() {
            for index in 0..count {
                let row = heap.slice(list.as_const())?[index as usize];
                if !row.is_null() {
                    inchi_free(heap, row)?;
                    heap.slice_mut(list)?[index as usize] = SourceMutPointer::null();
                }
            }
            inchi_free(heap, list)?;
            heap.slice_mut(v3000)?[0].lists_sterac = SourceMutPointer::null();
        }
        inchi_free(heap, v3000)?;
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn GetStructFromStdINCHI(
    heap: &mut SourceHeap,
    inpInChI: Option<&inchi_InputINCHI>,
    outStruct: &mut inchi_OutputStruct,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2461 GetStructFromStdINCHI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL GetStructFromStdINCHI( inchi_InputINCHI *inpInChI,
                                      inchi_OutputStruct *outStruct )
{
    if (( inpInChI ) &&
        ( inpInChI->szInChI ) &&
        ( strlen( inpInChI->szInChI ) >= LEN_INCHI_STRING_PREFIX + 3 ) &&
        ( inpInChI->szInChI[LEN_INCHI_STRING_PREFIX + 1] == 'S' ))
    {
         /* brief check indicated valid std input (more checks in GetStructFromINCHI) */
        return GetStructFromINCHI( inpInChI, outStruct );
    }
    else
    {
        /* non-std or just invalid input */
        return inchi_Ret_ERROR;
    }
}
    */
    // END INCHI C FUNCTION: GetStructFromStdINCHI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetStructFromStdINCHI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; LEN_INCHI_STRING_PREFIX=6.
    // INCHI✔️❌: SourceHeap ownership tracking and delegated GetStructFromINCHI execution are materially slower than native pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetStructFromStdINCHI

    let Some(input) = inpInChI else {
        return Ok(tagRetValGetINCHI_inchi_Ret_ERROR);
    };
    if input.szInChI.is_null() {
        return Ok(tagRetValGetINCHI_inchi_Ret_ERROR);
    }
    let length = source_strlen(heap, input.szInChI.as_const())?;
    if length < crate::source_types::LEN_INCHI_STRING_PREFIX as usize + 3
        || heap.slice(input.szInChI.as_const())?
            [crate::source_types::LEN_INCHI_STRING_PREFIX as usize + 1]
            != b'S' as i8
    {
        return Ok(tagRetValGetINCHI_inchi_Ret_ERROR);
    }
    GetStructFromINCHI(heap, Some(input), outStruct, stdout, build, clock_result)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetStructFromINCHIEx(
    heap: &mut SourceHeap,
    inpInChI: Option<&inchi_InputINCHI>,
    outStruct: &mut inchi_OutputStructEx,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2485 GetStructFromINCHIEx
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL GetStructFromINCHIEx( inchi_InputINCHI *inpInChI,
                                     inchi_OutputStructEx *outStruct )
{
    INCHI_CLOCK ic;
    CANON_GLOBALS CG;
    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;
    INPUT_PARMS inp_parms;
    INPUT_PARMS *ip = &inp_parms;
    INCHI_IOSTREAM inchi_file[3];
    INCHI_IOSTREAM *out_file = inchi_file, *log_file = inchi_file + 1, *input_file = inchi_file + 2;
    int    i, nRet = 0, nRet1;
    /* djb-rwth: removing redundant variables/code */
    int bReleaseVersion = bRELEASE_VERSION; /* djb-rwth: ignoring LLVM warning: variable used in function call */
    unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */
#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    int  num_repeat = REPEAT_ALL;
#endif
    static char szMainOption[] = " ?InChI2Struct";
    char szSdfDataValue[MAX_SDF_VALUE + 1];
    const char *argv[INCHI_MAX_NUM_ARG + 1];
    int   argc;
    char *szOptions = NULL;
    /* conversion result */
    inp_ATOM *at = NULL;
    int num_at = 0;
    OAD_Polymer *polymer = NULL;
    OAD_V3000    *v3000 = NULL;

#if( TRACE_MEMORY_LEAKS == 1 )
    _CrtSetDbgFlag( _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF );

/* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE
    _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode( _CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG );
#endif

    /* turn on floating point exceptions */
#if ( !defined(__STDC__) || __STDC__ != 1 )
    {
        /* Get the default control word. */
        int cw = _controlfp( 0, 0 );

        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &= ~( EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL );

        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif
#endif

    memset( outStruct, 0, sizeof( *outStruct ) ); /* djb-rwth: memset_s C11/Annex K variant? */

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    repeat:
    FreeStructFromINCHI( &outStruct );
    inchi_ios_reset( input_file );  /* do not close input_file - its string buffer may point to inpInChI->szInChI */
    inchi_ios_close( out_file );
    inchi_ios_close( log_file );
#endif

    sd->bUserQuit = 0;

    /* Initialize internal for this function I/O streams as string buffers */
    inchi_ios_init( input_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( out_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( log_file, INCHI_IOS_TYPE_STRING, NULL );

    /* clear original input structure */
    memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( ip, 0, sizeof( *ip ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( szSdfDataValue, 0, sizeof( szSdfDataValue ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    memset( &ic, 0, sizeof( ic ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( &CG, 0, sizeof( CG ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    szMainOption[1] = INCHI_OPTION_PREFX;

    if (!inpInChI)
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }

    /* options */
    if (inpInChI /*&& inpInChI->szOptions*/)
    {
        /* fix bug discovered by Burt Leland 2008-12-23 */
        int opt_len = ( inpInChI->szOptions ? strlen( inpInChI->szOptions ) : 0 ) + sizeof( szMainOption ) + 1;
        szOptions = (char*)inchi_calloc((long long)opt_len + 1, sizeof(szOptions[0])); /* djb-rwth: cast operator added */
        if (szOptions)
        {
            if (inpInChI->szOptions)
                /* fix bug discovered by Burt Leland 2008-12-23 */
                strcpy(szOptions, inpInChI->szOptions);
            strcat(szOptions, szMainOption);
            argc = parse_options_string( szOptions, argv, INCHI_MAX_NUM_ARG );
        }
        else
        {
            nRet = _IS_FATAL;
            goto translate_RetVal; /* emergency exit */
        }
    }
    else
    {
        argc = 1;
            argv[0] = "";
        argv[1] = NULL;
    }

    if ((argc == 1
#ifdef TARGET_API_LIB
        && ( !inpInChI || !inpInChI->szInChI ))
#endif
        || (argc == 2 && ( argv[1][0] == INCHI_OPTION_PREFX ) &&
        ( !strcmp( argv[1] + 1, "?" ) || !inchi_stricmp( argv[1] + 1, "help" ) ))) /* djb-rwth: addressing LLVM warnings */
    {
        HelpCommandLineParms( log_file );
        outStruct->szLog = log_file->s.pStr;
        nRet = _IS_EOF;
        goto translate_RetVal;
    }

    nRet1 = ReadCommandLineParms( argc, argv, ip, szSdfDataValue,
                                  &ulDisplTime, bReleaseVersion,
                                  log_file );

    if (szOptions)
    {
        /* argv pointed to strings in szOptions */
        inchi_free( szOptions );
        szOptions = NULL;
    }

    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if (0 > nRet1)
    {
        goto exit_function;
    }

    if (ip->bNoStructLabels)
    {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    }
    else if (ip->nInputType == INPUT_INCHI_XML ||
            ip->nInputType == INPUT_INCHI_PLAIN ||
            ip->nInputType == INPUT_CMLFILE ||
            ip->nInputType == INPUT_INCHI)
    {
        /* the input may contain both the header and the label of the structure */
        if (!ip->pSdfLabel)
            ip->pSdfLabel = ip->szSdfDataHeader;
        if (!ip->pSdfValue)
            ip->pSdfValue = szSdfDataValue;
    }

    if (ip->nInputType && ip->nInputType != INPUT_INCHI)
    {
        inchi_ios_eprint( log_file, "Input type set to INPUT_INCHI\n" );
        ip->nInputType = INPUT_INCHI;
    }

    if (!inpInChI->szInChI)
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }
    else
    {
        const int strict = 0;                     /* do not use strict mode, it may be too alarmous */
        nRet = CheckINCHI( inpInChI->szInChI, strict );
        if (nRet == INCHI_VALID_STANDARD || nRet == INCHI_VALID_NON_STANDARD || nRet == INCHI_VALID_BETA) /* djb-rwth: removing redundant code */
        {
            ;
        }
        else
        {
            nRet = _IS_ERROR;
            goto exit_function;
        }
    }

    PrintInputParms( log_file, ip );

    /*********************************/
    /* InChI -> Structure conversion */
    /*********************************/

    /* input_file simulation */

    /*
    that was incorrect simulation, and correct one is much simpler, see below
    input_file->s.pStr = inpInChI->szInChI;
    input_file->s.nUsedLength = (int) strlen(input_file->s.pStr)+1;
    input_file->s.nAllocatedLength = input_file->s.nUsedLength;
    input_file->s.nPtr = 0;
    */
    inchi_ios_print_nodisplay( input_file, inpInChI->szInChI );

    /* buffer for the message */
    /* outStruct->szMessage = (char *)inchi_calloc( MAX_MSG_LEN, sizeof(outStruct->szMessage[0])); */

    outStruct->szMessage = (char *) inchi_calloc( MAX_MSG_LEN, sizeof( char ) );
    if (!outStruct->szMessage)
    {
        inchi_ios_eprint( log_file, "Cannot allocate output message buffer.\n" );
        nRet = -1;
    }
    else
    {
        int num_bonds;
        nRet = ReadWriteInChI( &ic, &CG , input_file, out_file, log_file,
                                ip, sd, &at, &num_at, &num_bonds,
                                &polymer, &v3000,
                                outStruct->szMessage,
                                MAX_MSG_LEN, outStruct->WarningFlags );

        if (nRet >= 0 && polymer && at && (num_at > 0)) /* djb-rwth: fixing oss-fuzz issue #68329, #68286 */
        {
            OAD_Polymer_SmartReopenCyclizedUnits( polymer, at,
                                                 num_at, &num_bonds );
        }
    }

    if (nRet >= 0 && at && num_at)
    {
        /* success */
        nRet = InpAtom0DToInchiAtom( at, num_at,
                                    &outStruct->num_atoms,
                                    &outStruct->atom,
                                    &outStruct->num_stereo0D,
                                    &outStruct->stereo0D );

        if (at)
        {
            inchi_free( at );
            at = NULL;
        }

        if (nRet >= 0 && polymer)
        {
            /* Check for and then replace ZZ for star atoms if Polymer extension is supplied */
            for (i = 0; i < outStruct->num_atoms; i++)
            {
                if (!strcmp( outStruct->atom[i].elname, "Zz" ))
                {
                    strcpy( outStruct->atom[i].elname, "*" );
                }
            }
        }

        if (nRet >= 0)
        {
            if (polymer || v3000)
            {
                nRet = SetInChIExtInputByExtOrigAtData( polymer, v3000,
                                                        &outStruct->polymer,
                                                        &outStruct->v3000,
                                                        outStruct->num_atoms ); /* pair to SetExtOrigAtDataByInChIExtInput */
                FreeExtOrigAtData( polymer, v3000 );
                polymer = NULL;
                v3000 = NULL;
            }
        }
        if (nRet < 0)
        {
            inchi_ios_eprint( log_file, "Final structure conversion failed\n" );
        }
    }
    outStruct->szLog = log_file->s.pStr;

exit_function:;

    for (i = 0; i < MAX_NUM_PATHS; i++)
    {
        if (ip->path[i])
        {
            inchi_free( (char*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }

    SetBitFree( &CG );

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    if (num_repeat-- > 0)
    {
        goto repeat;
    }
#endif

#ifdef TARGET_API_LIB
    /* output */
    if (log_file->s.pStr && log_file->s.nUsedLength > 0)
    {
        while (log_file->s.nUsedLength &&
                '\n' == log_file->s.pStr[log_file->s.nUsedLength - 1])
        {
            log_file->s.pStr[--log_file->s.nUsedLength] = '\0'; /* remove last LF */
        }
        if (outStruct)
        {
            outStruct->szLog = log_file->s.pStr;
            log_file->s.pStr = NULL;
        }
    }
#endif

translate_RetVal:

    /* Close internal I/O streams */
    /* that was incorrect also
    inchi_ios_reset(input_file);  */    /* do not close input_file - its string buffer may point to inpInChI->szInChI */
    inchi_ios_close( input_file );
    inchi_ios_close( out_file );
    inchi_ios_close( log_file );

    switch (nRet)
    {
        case -3: nRet = inchi_Ret_ERROR; break; /* Error: no Structure has been created */
        case -2: nRet = inchi_Ret_ERROR; break; /* Error: no Structure has been created */
        case -1: nRet = inchi_Ret_FATAL; break; /* Severe error: no Structure has been created (typically; break; memory allocation failed) */
        default:
            if (outStruct) /* djb-rwth: fixing a NULL pointer dereference */
            {
                if (!outStruct->atom || !outStruct->num_atoms)
                {
                    nRet = inchi_Ret_EOF;
                }
                else
                {
                    int m, n, t = 0;
                    for (m = 0; m < 2; m++)
                    {
                        for (n = 0; n < 2; n++)
                        {
                            if (outStruct->WarningFlags[m][n])
                            {
                                t++;
                            }
                        }
                    }
                    nRet = t ? inchi_Ret_WARNING : inchi_Ret_OKAY;
                }
                break;
            }
    }

    return nRet;
}
    */
    // END INCHI C FUNCTION: GetStructFromINCHIEx
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetStructFromINCHIEx
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; bRELEASE_VERSION=1 and INCHI_OPTION_PREFX='-'.
    // INCHI✔️❌: TRACE_MEMORY_LEAKS=0, REPEAT_ALL undefined, and the non-ANSI floating-point branch is inactive.
    // INCHI✔️❌: The active TARGET_API_LIB branch strips trailing LF bytes and transfers the log allocation to outStruct.
    // INCHI✔️❌: SourceHeap stream records, checked pointer access, snapshots, and bytewise ownership modeling are materially slower than native pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetStructFromINCHIEx

    const MAX_MESSAGE_LENGTH: usize = 512;
    const MAIN_OPTION: &[u8] = b" -InChI2Struct\0";

    *outStruct = inchi_OutputStructEx::default();

    let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()])?;
    let canon = heap.allocate_model_storage(vec![CANON_GLOBALS::default()])?;
    let sdf_value = heap.allocate_model_storage(vec![0_i8; MAX_SDF_VALUE as usize + 1])?;
    let empty_argument = heap.allocate_model_storage(vec![0_i8])?;
    let print_format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;

    let mut streams: [INCHI_IOSTREAM; 3] = std::array::from_fn(|_| INCHI_IOSTREAM::default());
    for stream in &mut streams {
        inchi_ios_init(
            Some(stream),
            INCHI_IOS_TYPE_STRING as i32,
            SourceMutPointer::null(),
        )?;
    }

    let mut structure_data = STRUCT_DATA::default();
    let mut input_parameters = INPUT_PARMS::default();
    let mut display_time = 0_u64;
    let mut raw_return = 0_i32;
    let mut translate_directly = false;
    let mut atoms = SourceMutPointer::<inp_ATOM>::null();
    let mut atom_count = 0_i32;
    let mut polymer = SourceMutPointer::<OAD_Polymer>::null();
    let mut v3000 = SourceMutPointer::<OAD_V3000>::null();

    'source: {
        let Some(input) = inpInChI else {
            raw_return = _IS_ERROR as i32;
            break 'source;
        };

        let option_length = if input.szOptions.is_null() {
            0
        } else {
            source_strlen(heap, input.szOptions.as_const())?
        };
        let source_option_length = option_length
            .checked_add(MAIN_OPTION.len())
            .and_then(|value| value.checked_add(2))
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        let options = match inchi_calloc::<i8>(heap, source_option_length as u64, 1) {
            Ok(pointer) => pointer,
            Err(_) => {
                raw_return = _IS_FATAL as i32;
                translate_directly = true;
                break 'source;
            }
        };
        if option_length != 0 {
            let source = heap.slice(input.szOptions.as_const())?[..option_length].to_vec();
            heap.slice_mut(options)?[..option_length].copy_from_slice(&source);
        }
        let main = MAIN_OPTION
            .iter()
            .map(|byte| *byte as i8)
            .collect::<Vec<_>>();
        heap.slice_mut(options)?[option_length..option_length + MAIN_OPTION.len()]
            .copy_from_slice(&main);

        let mut parsed_arguments = vec![SourceArgvPointer::Null; INCHI_MAX_NUM_ARG as usize + 1];
        let argument_count = parse_options_string(
            heap,
            options,
            &mut parsed_arguments,
            INCHI_MAX_NUM_ARG as i32,
        )?;
        let mut argv = vec![SourceConstPointer::null(); argument_count as usize];
        for (index, argument) in parsed_arguments[..argument_count as usize]
            .iter()
            .copied()
            .enumerate()
        {
            argv[index] = match argument {
                SourceArgvPointer::EmptyLiteral => empty_argument.as_const(),
                SourceArgvPointer::Command(pointer) => pointer.as_const(),
                SourceArgvPointer::Null => return Err(SourceHeapError::NullPointer),
            };
        }

        let help_requested = if argument_count == 1 && input.szInChI.is_null() {
            true
        } else if argument_count == 2 {
            let argument = argv[1];
            let bytes = heap.slice(argument)?;
            bytes.first().copied() == Some(b'-' as i8)
                && (source_string_equals(heap, argument.offset(1)?, b"?")? || {
                    let help = heap.allocate_model_storage(
                        b"help\0".iter().map(|byte| *byte as i8).collect(),
                    )?;
                    let equal = inchi_stricmp(heap, argument.offset(1)?, help.as_const())? == 0;
                    heap.free(help)?;
                    equal
                })
        } else {
            false
        };
        if help_requested {
            HelpCommandLineParms(heap, Some(&mut streams[1]), stdout, build)?;
            outStruct.szLog = streams[1].s.pStr;
            raw_return = _IS_EOF as i32;
            translate_directly = true;
            break 'source;
        }

        let parameter_return = ReadCommandLineParms(
            heap,
            argument_count,
            &argv,
            &mut input_parameters,
            sdf_value,
            &mut display_time,
            1,
            Some(&mut streams[1]),
        )?;
        inchi_free(heap, options)?;
        input_parameters.bNoStructLabels = 1;

        if parameter_return < 0 {
            break 'source;
        }
        input_parameters.pSdfLabel = SourceMutPointer::null();
        input_parameters.pSdfValue = SourceMutPointer::null();

        if input_parameters.nInputType != 0
            && input_parameters.nInputType != tagInputType_INPUT_INCHI
        {
            append_api_log(heap, &mut streams[1], "Input type set to INPUT_INCHI\n")?;
            input_parameters.nInputType = tagInputType_INPUT_INCHI;
        }

        if input.szInChI.is_null() {
            raw_return = _IS_ERROR as i32;
            break 'source;
        }
        let validity = CheckINCHI(
            heap,
            input.szInChI.as_const(),
            0,
            stdout,
            build,
            clock_result,
        )?;
        if !matches!(
            validity,
            tagRetValCheckINCHI_INCHI_VALID_STANDARD
                | tagRetValCheckINCHI_INCHI_VALID_NON_STANDARD
                | tagRetValCheckINCHI_INCHI_VALID_BETA
        ) {
            raw_return = _IS_ERROR as i32;
            break 'source;
        }

        PrintInputParms(heap, Some(&mut streams[1]), &input_parameters)?;

        let arguments = SourceVaList {
            arguments: vec![crate::source_types::SourceFormatArgument::Bytes(
                input.szInChI.as_const(),
            )],
            position: 0,
        };
        let _ = inchi_ios_print_nodisplay(
            heap,
            Some(&mut streams[2]),
            stdout,
            print_format.as_const(),
            &arguments,
        )?;

        outStruct.szMessage = match inchi_calloc::<i8>(heap, MAX_MESSAGE_LENGTH as u64, 1) {
            Ok(pointer) => pointer,
            Err(_) => {
                append_api_log(
                    heap,
                    &mut streams[1],
                    "Cannot allocate output message buffer.\n",
                )?;
                raw_return = -1;
                SourceMutPointer::null()
            }
        };

        if !outStruct.szMessage.is_null() {
            let message = outStruct.szMessage;
            let mut bond_count = 0_i32;
            let (output_stream, remaining_streams) = streams.split_at_mut(1);
            let (log_stream, input_stream) = remaining_streams.split_at_mut(1);
            raw_return = heap.with_slice_mut_and_heap_mut(message, |message, heap| {
                ReadWriteInChI(
                    heap,
                    clock,
                    canon,
                    &mut input_stream[0],
                    &mut output_stream[0],
                    &mut log_stream[0],
                    &mut input_parameters,
                    &mut structure_data,
                    Some(&mut atoms),
                    Some(&mut atom_count),
                    Some(&mut bond_count),
                    Some(&mut polymer),
                    Some(&mut v3000),
                    Some(message),
                    MAX_MESSAGE_LENGTH as i32,
                    Some(&mut outStruct.WarningFlags),
                    stdout,
                    clock_result,
                )
            })?;

            if raw_return >= 0 && !polymer.is_null() && !atoms.is_null() && atom_count > 0 {
                OAD_Polymer_SmartReopenCyclizedUnits(
                    heap,
                    polymer,
                    atoms,
                    atom_count,
                    &mut bond_count,
                )?;
            }
        }

        if raw_return >= 0 && !atoms.is_null() && atom_count != 0 {
            raw_return = InpAtom0DToInchiAtom(
                heap,
                atoms,
                atom_count,
                Some(&mut outStruct.num_atoms),
                Some(&mut outStruct.atom),
                Some(&mut outStruct.num_stereo0D),
                Some(&mut outStruct.stereo0D),
            )?;

            if !atoms.is_null() {
                inchi_free(heap, atoms)?;
                atoms = SourceMutPointer::null();
            }

            if raw_return >= 0 && !polymer.is_null() {
                for index in 0..outStruct.num_atoms {
                    let atom = &mut heap.slice_mut(outStruct.atom)?[index as usize];
                    let nul = atom
                        .elname
                        .iter()
                        .position(|byte| *byte == 0)
                        .ok_or(SourceHeapError::MissingNulTerminator)?;
                    if atom.elname[..nul] == [b'Z' as i8, b'z' as i8] {
                        atom.elname[0] = b'*' as i8;
                        atom.elname[1] = 0;
                    }
                }
            }

            if raw_return >= 0 && (!polymer.is_null() || !v3000.is_null()) {
                raw_return = SetInChIExtInputByExtOrigAtData(
                    heap,
                    polymer,
                    v3000,
                    &mut outStruct.polymer,
                    &mut outStruct.v3000,
                    i32::from(outStruct.num_atoms),
                )?;
                FreeExtOrigAtData(heap, polymer, v3000)?;
                polymer = SourceMutPointer::null();
                v3000 = SourceMutPointer::null();
            }
            if raw_return < 0 {
                append_api_log(heap, &mut streams[1], "Final structure conversion failed\n")?;
            }
        }
        outStruct.szLog = streams[1].s.pStr;
    }

    if !translate_directly {
        for path in input_parameters.path.iter_mut().take(MAX_NUM_PATHS as usize) {
            if !path.is_null() {
                inchi_free(heap, path.as_mut())?;
                *path = SourceConstPointer::null();
            }
        }
        let mut globals = heap.slice(canon.as_const())?[0].clone();
        SetBitFree(heap, &mut globals)?;
        heap.slice_mut(canon)?[0] = globals;

        if !streams[1].s.pStr.is_null() && streams[1].s.nUsedLength > 0 {
            while streams[1].s.nUsedLength > 0 {
                let last = usize::try_from(streams[1].s.nUsedLength - 1)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if heap.slice(streams[1].s.pStr.as_const())?[last] != b'\n' as i8 {
                    break;
                }
                streams[1].s.nUsedLength -= 1;
                heap.slice_mut(streams[1].s.pStr)?[last] = 0;
            }
            outStruct.szLog = streams[1].s.pStr;
            streams[1].s.pStr = SourceMutPointer::null();
        }
    }

    inchi_ios_close(heap, Some(&mut streams[2]))?;
    inchi_ios_close(heap, Some(&mut streams[0]))?;
    inchi_ios_close(heap, Some(&mut streams[1]))?;

    raw_return = match raw_return {
        -3 | -2 => tagRetValGetINCHI_inchi_Ret_ERROR,
        -1 => tagRetValGetINCHI_inchi_Ret_FATAL,
        _ if outStruct.atom.is_null() || outStruct.num_atoms == 0 => {
            tagRetValGetINCHI_inchi_Ret_EOF
        }
        _ => {
            let warnings = outStruct
                .WarningFlags
                .iter()
                .flatten()
                .filter(|flag| **flag != 0)
                .count();
            if warnings != 0 {
                tagRetValGetINCHI_inchi_Ret_WARNING
            } else {
                tagRetValGetINCHI_inchi_Ret_OKAY
            }
        }
    };

    heap.free(sdf_value)?;
    heap.free(empty_argument)?;
    heap.free(print_format)?;
    heap.free(clock)?;
    heap.free(canon)?;
    Ok(raw_return)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn GetStructFromINCHI(
    heap: &mut SourceHeap,
    inpInChI: Option<&inchi_InputINCHI>,
    out: &mut inchi_OutputStruct,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2856 GetStructFromINCHI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL GetStructFromINCHI( inchi_InputINCHI *inpInChI,
                                   inchi_OutputStruct *out )
{
    int ret = 0;

    inchi_OutputStructEx outex;
    memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    ret = GetStructFromINCHIEx( inpInChI, &outex );

    out->szLog = outex.szLog;
    out->szMessage = outex.szMessage;
    out->WarningFlags[0][0] = outex.WarningFlags[0][0];
    out->WarningFlags[0][1] = outex.WarningFlags[0][1];
    out->WarningFlags[1][0] = outex.WarningFlags[1][0];
    out->WarningFlags[1][1] = outex.WarningFlags[1][1];

    if (ret == inchi_Ret_OKAY || ret == inchi_Ret_WARNING)
    {
        out->num_atoms = outex.num_atoms;
        out->atom = outex.atom;
        out->num_stereo0D = outex.num_stereo0D;
        out->stereo0D = outex.stereo0D;
    }

    return ret;
}
    */
    // END INCHI C FUNCTION: GetStructFromINCHI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetStructFromINCHI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no conditional branch occurs in this function.
    // INCHI✔️❌: SourceHeap ownership tracking and the complete GetStructFromINCHIEx Rust call graph are materially slower than native pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetStructFromINCHI

    *out = inchi_OutputStruct::default();
    let mut extended = inchi_OutputStructEx::default();
    let result = GetStructFromINCHIEx(heap, inpInChI, &mut extended, stdout, build, clock_result)?;

    out.szLog = extended.szLog;
    out.szMessage = extended.szMessage;
    out.WarningFlags[0][0] = extended.WarningFlags[0][0];
    out.WarningFlags[0][1] = extended.WarningFlags[0][1];
    out.WarningFlags[1][0] = extended.WarningFlags[1][0];
    out.WarningFlags[1][1] = extended.WarningFlags[1][1];

    if result == tagRetValGetINCHI_inchi_Ret_OKAY || result == tagRetValGetINCHI_inchi_Ret_WARNING {
        out.num_atoms = extended.num_atoms;
        out.atom = extended.atom;
        out.num_stereo0D = extended.num_stereo0D;
        out.stereo0D = extended.stereo0D;
    }

    Ok(result)
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

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn FreeStructFromINCHIEx(
    heap: &mut SourceHeap,
    output: Option<&mut inchi_OutputStructEx>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2892 FreeStructFromINCHIEx
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void INCHI_DECL FreeStructFromINCHIEx( inchi_OutputStructEx *out )
{
    if (!out)
        return;

    if (out->atom)
    {
        inchi_free( out->atom );
    }
    if (out->stereo0D)
    {
        inchi_free( out->stereo0D );
    }
    if (out->szLog)
    {
        inchi_free( out->szLog );
    }
    if (out->szMessage)
    {
        inchi_free( out->szMessage );
    }
    if (out->polymer || out->v3000)
    {
        FreeInChIExtInput( out->polymer, out->v3000 );
    }

    memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */
}
    */
    // END INCHI C FUNCTION: FreeStructFromINCHIEx
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FreeStructFromINCHIEx
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local conditional compilation.
    // INCHI✔️❌: SourceHeap allocation-map operations and nested pointer validation are materially slower than native inchi_free calls.
    // END INCHI ACTIVE MACRO CONFIGURATION: FreeStructFromINCHIEx

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
    if !output.polymer.is_null() || !output.v3000.is_null() {
        FreeInChIExtInput(heap, output.polymer, output.v3000)?;
    }
    *output = inchi_OutputStructEx::default();
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

fn source_strlen(
    heap: &SourceHeap,
    pointer: SourceConstPointer<i8>,
) -> Result<usize, SourceHeapError> {
    heap.slice(pointer)?
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)
}

fn source_string_equals(
    heap: &SourceHeap,
    pointer: SourceConstPointer<i8>,
    expected: &[u8],
) -> Result<bool, SourceHeapError> {
    let length = source_strlen(heap, pointer)?;
    Ok(heap.slice(pointer)?[..length]
        .iter()
        .map(|byte| *byte as u8)
        .eq(expected.iter().copied()))
}

fn append_api_log(
    heap: &mut SourceHeap,
    log: &mut INCHI_IOSTREAM,
    text: &str,
) -> Result<(), SourceHeapError> {
    let format = heap.allocate_model_storage(
        text.bytes()
            .chain(std::iter::once(0))
            .map(|byte| byte as i8)
            .collect(),
    )?;
    let result = inchi_ios_eprint(heap, Some(log), format.as_const(), &SourceVaList::default());
    heap.free(format)?;
    result.map(|_| ())
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetINCHI1(
    heap: &mut SourceHeap,
    extended_input: Option<&inchi_InputEx>,
    mut output: Option<&mut inchi_Output>,
    enforce_standard_format: i32,
    interrupted: i32,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:345 GetINCHI1
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
static int GetINCHI1( inchi_InputEx *extended_input,
                      inchi_Output *out,
                      int enforce_std_format )
{
    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;
    char szTitle[MAX_SDF_HEADER + MAX_SDF_VALUE + 256];

    int i;
    long num_inp, num_err; /* djb-rwth: ignoring LLVM warning: variable used */
    char      szSdfDataValue[MAX_SDF_VALUE + 1];
    PINChI2     *pINChI[INCHI_NUM];
    PINChI_Aux2 *pINChI_Aux[INCHI_NUM];

    unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */
    unsigned long  ulTotalProcessingTime = 0; /* djb-rwth: ignoring LLVM warning: variable used */

    INPUT_PARMS inp_parms;
    INPUT_PARMS *ip = &inp_parms;

    ORIG_ATOM_DATA OrigAtData; /* 0=> disconnected, 1=> original */
    ORIG_ATOM_DATA *orig_inp_data = &OrigAtData;
    ORIG_ATOM_DATA PrepAtData[2]; /* 0=> disconnected, 1=> original */
    ORIG_ATOM_DATA *prep_inp_data = PrepAtData;
    int             bReleaseVersion = bRELEASE_VERSION;
    int   nRet = 0, nRet1;

    CANON_GLOBALS CG;
    INCHI_CLOCK ic;

    STRUCT_FPTRS *pStructPtrs = NULL;

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    int  num_repeat = REPEAT_ALL;
#endif

    const char *argv[INCHI_MAX_NUM_ARG + 1];
    int   argc;
    char *szOptions = NULL;

    INCHI_IOSTREAM inchi_file[3], *out_file = inchi_file, *log_file = inchi_file + 1;
    INCHI_IOSTREAM prb_file0, *prb_file = &prb_file0;
    INCHI_IOS_STRING temp_string_container;
    INCHI_IOS_STRING *strbuf = &temp_string_container;

    inchi_Input prev_versions_input;
    inchi_Input *pvinp = &prev_versions_input;

#ifdef GHI100_FIX
#if ((SPRINTF_FLAG != 1) && (SPRINTF_FLAG != 2))
    setlocale(LC_ALL, "en-US"); /* djb-rwth: setting all locales to "en-US" */
#endif
#endif

    pvinp->atom = extended_input->atom;
    pvinp->num_atoms = extended_input->num_atoms;
    pvinp->num_stereo0D = extended_input->num_stereo0D;
    pvinp->stereo0D = extended_input->stereo0D;
    pvinp->szOptions = extended_input->szOptions;

#if( TRACE_MEMORY_LEAKS == 1 )
    _CrtSetDbgFlag( _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF );
/* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE
    _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode( _CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG );
#endif

#if ( !defined(__STDC__) || __STDC__ != 1 )
    /* turn on floating point exceptions */
    {
        /* Get the default control word. */
        int cw = _controlfp( 0, 0 );

        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &= ~( EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL );

        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif
#endif

    szTitle[0] = '\0';

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    repeat:
          inchi_ios_close( out_file );
          inchi_ios_close( log_file );
          inchi_ios_close( prb_file );
          pStr = NULL;
#endif

    /* Initialize internal for this function output streams as string buffers */
    inchi_ios_init( out_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( log_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( prb_file, INCHI_IOS_TYPE_STRING, NULL );

    num_inp = 0;
    num_err = 0;
    sd->bUserQuit = 0;

    /* clear original input structure */
    memset( pINChI, 0, sizeof( pINChI ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( pINChI_Aux, 0, sizeof( pINChI_Aux ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( ip, 0, sizeof( *ip ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( orig_inp_data, 0, sizeof( *orig_inp_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( prep_inp_data, 0, 2 * sizeof( *prep_inp_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( szSdfDataValue, 0, sizeof( szSdfDataValue ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    memset( &CG, 0, sizeof( CG ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( &ic, 0, sizeof( ic ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    if (!out)
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }
    memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* options */
    if (pvinp && pvinp->szOptions)
    {
        szOptions = (char*) inchi_malloc( strlen( pvinp->szOptions ) + 1 );
        if (szOptions)
        {
            strcpy( szOptions, pvinp->szOptions );
            argc = parse_options_string( szOptions, argv, INCHI_MAX_NUM_ARG );
        }
        else
        {
            nRet = _IS_FATAL;
            goto translate_RetVal; /* emergency exit */
        }
    }
    else
    {
        argc = 1;
        argv[0] = "";
        argv[1] = NULL;
    }

    if ((argc == 1
#ifdef TARGET_API_LIB
              && ( !pvinp || pvinp->num_atoms <= 0 || !pvinp->atom ))
#endif
              || (argc == 2 && ( argv[1][0] == INCHI_OPTION_PREFX ) &&
                    ( !strcmp( argv[1] + 1, "?" ) || !inchi_stricmp( argv[1] + 1, "help" )) )) /* djb-rwth: addressing LLVM warnings */
    {
        HelpCommandLineParms( log_file );
        out->szLog = log_file->s.pStr;
        memset( log_file, 0, sizeof( *log_file ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        nRet = _IS_EOF;
        goto translate_RetVal;
    }

    nRet1 = ReadCommandLineParms( argc, argv, ip, szSdfDataValue, &ulDisplTime, bReleaseVersion, log_file );
    if (szOptions)
    {
        inchi_free( szOptions );
        szOptions = NULL;
    }
    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if (0 > nRet1)
    {
        nRet = _IS_FATAL;
        goto exit_function;
    }
    if (ip->bNoStructLabels)
    {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    }
    else
    {
        if (ip->nInputType == INPUT_INCHI_XML || ip->nInputType == INPUT_INCHI_PLAIN || ip->nInputType == INPUT_CMLFILE)
        {
            /* the input may contain both the header and the label of the structure */
            if (!ip->pSdfLabel)
                ip->pSdfLabel = ip->szSdfDataHeader;
            if (!ip->pSdfValue)
                ip->pSdfValue = szSdfDataValue;
        }
    }

    /* Ensure standardness */
    if (enforce_std_format)
    {
        if (ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT)
        {
            ip->bINChIOutputOptions &= ~INCHI_OUT_SAVEOPT;
        }
        if (0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ))
        {
            ip->bTautFlags &= ~TG_FLAG_RECONNECT_COORD;
        }
        if (0 != ( ip->nMode & REQ_MODE_BASIC ))
        {
            ip->nMode &= ~REQ_MODE_BASIC;
        }
        if (0 != ( ip->nMode & REQ_MODE_RELATIVE_STEREO ))
        {
            ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO );
        }
        if (0 != ( ip->nMode & REQ_MODE_RACEMIC_STEREO ))
        {
            ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO );
        }
        if (0 != ( ip->nMode & REQ_MODE_CHIR_FLG_STEREO ))
        {
            ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO );
        }
        if (0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO ))
        {
            ip->nMode &= ~REQ_MODE_DIFF_UU_STEREO;
        }
        if (0 == ( ip->nMode & ( REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU ) ))
        {
            ip->nMode |= REQ_MODE_SB_IGN_ALL_UU;
            ip->nMode |= REQ_MODE_SC_IGN_ALL_UU;
        }
        if (0 != ( ip->bTautFlags & TG_FLAG_KETO_ENOL_TAUT ))
        {
            ip->bTautFlags &= ~TG_FLAG_KETO_ENOL_TAUT;
        }
        if (0 != ( ip->bTautFlags & TG_FLAG_1_5_TAUT ))
        {
            ip->bTautFlags &= ~TG_FLAG_1_5_TAUT;
        }
        /* And anyway... */
        ip->bINChIOutputOptions |= INCHI_OUT_STDINCHI;
        ip->bINChIOutputOptions &= ~INCHI_OUT_SAVEOPT;
    }
    /* */

    PrintInputParms( log_file, ip );

    if (0 >= inchi_strbuf_init( strbuf, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT ))
    {
        inchi_ios_eprint( log_file, "Cannot allocate internal string buffer. Terminating\n" );
        nRet = _IS_FATAL;
        goto exit_function;
    }

    /***************************************************
    Main cycle -- read input structures and create their INChI's */ /* djb-rwth: addressing LLVM warning */
    ulTotalProcessingTime = 0;

    if (pStructPtrs)
    {
        memset( pStructPtrs, 0, sizeof( pStructPtrs[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    }

    /* === possible improvement: convert inp to orig_inp_data ==== */
    if (!sd->bUserQuit && !bInterrupted)
    {
        if (ip->last_struct_number && num_inp >= ip->last_struct_number)
        {
            nRet = _IS_EOF; /*  simulate end of file */
            goto exit_function;
        }

        nRet = ExtractOneStructure( sd,ip, szTitle, extended_input,
                                    log_file, out_file, prb_file,
                                    orig_inp_data, &num_inp );

        if (pStructPtrs)
        {
            pStructPtrs->cur_fptr++;
        }

#ifndef TARGET_API_LIB
        if (sd->bUserQuit)
        {
            break;
        }
#endif
        switch (nRet)
        {
            case _IS_FATAL:
                num_err++;
                goto exit_function;
            case _IS_EOF:
                goto exit_function;
            case _IS_ERROR:
                num_err++;
                goto exit_function;
#ifndef TARGET_API_LIB
            case _IS_SKIP:
                continue;
#endif
        }

        /* Create INChI for each connected component of the structure and */
        /* optionally display them ; output INChI for the whole structure */

        nRet1 = ProcessOneStructureEx( &ic, &CG, sd, ip, szTitle,
                                        pINChI, pINChI_Aux,
                                        NULL, /* inp_file is not necessary as all input is already saved in 'ip' */
                                        log_file, out_file, prb_file,
                                        orig_inp_data, prep_inp_data,
                                        num_inp, strbuf, 0 /* save_opt_bits */ );

        /*  Free INChI memory */
        FreeAllINChIArrays( pINChI, pINChI_Aux, sd->num_components );

        /* Free structure data */
        FreeOrigAtData( orig_inp_data );
        FreeOrigAtData( prep_inp_data );
        FreeOrigAtData( prep_inp_data + 1 );

        ulTotalProcessingTime += sd->ulStructTime;
        nRet = inchi_max( nRet, nRet1 );
        switch (nRet)
        {
            case _IS_FATAL:
                /* num_err ++; */
                goto exit_function;
            case _IS_ERROR:
                ; /* num_err ++; */
#ifndef TARGET_API_LIB
                continue;
#endif
        }
    }

exit_function:
    /* Avoid memory leaks in case of fatal error */
    if (pStructPtrs && pStructPtrs->fptr)
    {
        inchi_free( pStructPtrs->fptr );
    }
    /* Free INChI memory */
    FreeAllINChIArrays( pINChI, pINChI_Aux, sd->num_components );
    /*    Free structure data */
    FreeOrigAtData( orig_inp_data );
    FreeOrigAtData( prep_inp_data );
    FreeOrigAtData( prep_inp_data + 1 );

    inchi_strbuf_close( strbuf );

    for (i = 0; i < MAX_NUM_PATHS; i++)
    {
        if (ip->path[i])
        {
            inchi_free( (char*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }

    SetBitFree( &CG );

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    if (num_repeat-- > 0)
    {
        goto repeat;
    }
#endif

    /* output */
    produce_generation_output( out, sd, ip, log_file, out_file );

translate_RetVal:

    /* Close inernal I/O streams */
    inchi_ios_close( log_file );
    inchi_ios_close( out_file );
    inchi_ios_close( prb_file );

    switch (nRet)
    {
        case _IS_SKIP: nRet = inchi_Ret_SKIP; break; /* not used in INChI dll */
        case _IS_EOF: nRet = inchi_Ret_EOF; break; /* no structural data has been provided */
        case _IS_OKAY: nRet = inchi_Ret_OKAY; break; /* Success; break; no errors or warnings */
        case _IS_WARNING: nRet = inchi_Ret_WARNING; break; /* Success; break; warning(s) issued */
        case _IS_ERROR: nRet = inchi_Ret_ERROR; break; /* Error: no INChI has been created */
        case _IS_FATAL: nRet = inchi_Ret_FATAL; break; /* Severe error: no INChI has been created (typically; break; memory allocation failed) */
        case _IS_UNKNOWN:
        default: nRet = inchi_Ret_UNKNOWN; break; /* Unlnown program error */
    }

    return nRet;
}
    */
    // END INCHI C FUNCTION: GetINCHI1
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetINCHI1
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; bRELEASE_VERSION=1.
    // INCHI✔️❌: TRACE_MEMORY_LEAKS=0, REPEAT_ALL undefined, GHI100_FIX undefined, and standalone branches inactive.
    // INCHI✔️❌: The mutable bInterrupted global is an explicit argument; long and unsigned long are modeled as i64 and u64.
    // INCHI✔️❌: SourceHeap stream/root records and ABI literal storage add map access and cloning absent from native stack pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetINCHI1

    let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()])?;
    let canonical = heap.allocate_model_storage(vec![CANON_GLOBALS::default()])?;
    let title = heap.allocate_model_storage(vec![
        0_i8;
        (MAX_SDF_HEADER + MAX_SDF_VALUE + 256) as usize
    ])?;
    let sdf_value = heap.allocate_model_storage(vec![0_i8; MAX_SDF_VALUE as usize + 1])?;
    let empty_argument = heap.allocate_model_storage(vec![0_i8])?;
    let original = heap.allocate_model_storage(vec![ORIG_ATOM_DATA::default()])?;
    let prepared = heap.allocate_model_storage(vec![ORIG_ATOM_DATA::default(); 2])?;

    let mut streams: [INCHI_IOSTREAM; 3] = std::array::from_fn(|_| INCHI_IOSTREAM::default());
    let mut problem_stream = INCHI_IOSTREAM::default();
    inchi_ios_init(Some(&mut streams[0]), INCHI_IOS_TYPE_STRING as i32, SourceMutPointer::null())?;
    inchi_ios_init(Some(&mut streams[1]), INCHI_IOS_TYPE_STRING as i32, SourceMutPointer::null())?;
    inchi_ios_init(Some(&mut problem_stream), INCHI_IOS_TYPE_STRING as i32, SourceMutPointer::null())?;

    let mut structure_data = STRUCT_DATA::default();
    let mut input_parameters = INPUT_PARMS::default();
    let mut string_buffer = INCHI_IOS_STRING::default();
    let mut inchi_components = [SourceMutPointer::null(); crate::source_types::INCHI_NUM as usize];
    let mut aux_components = [SourceMutPointer::null(); crate::source_types::INCHI_NUM as usize];
    let mut display_time = 0_u64;
    let mut input_number = 0_i64;
    let mut raw_return = 0_i32;
    let mut translate_directly = false;

    'source: {
        let Some(out) = output.as_deref_mut() else {
            raw_return = _IS_ERROR as i32;
            break 'source;
        };
        *out = inchi_Output::default();

        let option_pointer = extended_input
            .map(|input| input.szOptions)
            .unwrap_or_else(SourceMutPointer::null);
        let mut options = SourceMutPointer::null();
        let mut parsed_arguments = vec![SourceArgvPointer::Null; INCHI_MAX_NUM_ARG as usize + 1];
        let argument_count = if option_pointer.is_null() {
            parsed_arguments[0] = SourceArgvPointer::EmptyLiteral;
            1
        } else {
            let length = source_strlen(heap, option_pointer.as_const())?;
            options = match inchi_malloc(heap, (length + 1) as u64) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => {
                    raw_return = _IS_FATAL as i32;
                    translate_directly = true;
                    break 'source;
                }
                Err(error) => return Err(error),
            };
            let bytes = heap.slice(option_pointer.as_const())?[..=length].to_vec();
            heap.slice_mut(options)?[..=length].copy_from_slice(&bytes);
            parse_options_string(
                heap,
                options,
                &mut parsed_arguments,
                INCHI_MAX_NUM_ARG as i32,
            )?
        };

        let mut argv = Vec::with_capacity(argument_count as usize);
        for argument in parsed_arguments.iter().copied().take(argument_count as usize) {
            argv.push(match argument {
                SourceArgvPointer::EmptyLiteral => empty_argument.as_const(),
                SourceArgvPointer::Command(pointer) => pointer.as_const(),
                SourceArgvPointer::Null => return Err(SourceHeapError::NullPointer),
            });
        }

        let no_input = extended_input.is_none_or(|input| input.num_atoms <= 0 || input.atom.is_null());
        let help = if argument_count == 2 {
            let argument = argv[1];
            let first = heap.slice(argument)?.first().copied();
            first == Some(b'-' as i8)
                && (source_string_equals(heap, argument.offset(1)?, b"?")? || {
                    let help = heap.allocate_model_storage(b"help\0".iter().map(|byte| *byte as i8).collect())?;
                    let equal = inchi_stricmp(heap, argument.offset(1)?, help.as_const())? == 0;
                    heap.free(help)?;
                    equal
                })
        } else {
            false
        };
        if (argument_count == 1 && no_input) || help {
            HelpCommandLineParms(heap, Some(&mut streams[1]), stdout, build)?;
            out.szLog = streams[1].s.pStr;
            streams[1] = INCHI_IOSTREAM::default();
            raw_return = _IS_EOF as i32;
            translate_directly = true;
            break 'source;
        }

        let parameter_return = ReadCommandLineParms(
            heap,
            argument_count,
            &argv,
            &mut input_parameters,
            sdf_value,
            &mut display_time,
            1,
            Some(&mut streams[1]),
        )?;
        if !options.is_null() {
            inchi_free(heap, options)?;
        }
        input_parameters.bNoStructLabels = 1;
        if parameter_return < 0 {
            raw_return = _IS_FATAL as i32;
            break 'source;
        }
        input_parameters.pSdfLabel = SourceMutPointer::null();
        input_parameters.pSdfValue = SourceMutPointer::null();

        if enforce_standard_format != 0 {
            input_parameters.bINChIOutputOptions &= !(INCHI_OUT_SAVEOPT as i32);
            input_parameters.bTautFlags &= !(u64::from(TG_FLAG_RECONNECT_COORD));
            input_parameters.nMode &= !(u64::from(REQ_MODE_BASIC));
            if input_parameters.nMode
                & u64::from(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO | REQ_MODE_CHIR_FLG_STEREO)
                != 0
            {
                input_parameters.nMode &= !u64::from(
                    REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO,
                );
            }
            input_parameters.nMode &= !(u64::from(REQ_MODE_DIFF_UU_STEREO));
            if input_parameters.nMode & u64::from(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU) == 0 {
                input_parameters.nMode |= u64::from(REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU);
            }
            input_parameters.bTautFlags &= !u64::from(TG_FLAG_KETO_ENOL_TAUT | TG_FLAG_1_5_TAUT);
            input_parameters.bINChIOutputOptions |= INCHI_OUT_STDINCHI as i32;
            input_parameters.bINChIOutputOptions &= !(INCHI_OUT_SAVEOPT as i32);
        }

        PrintInputParms(heap, Some(&mut streams[1]), &input_parameters)?;
        if inchi_strbuf_init(
            heap,
            &mut string_buffer,
            INCHI_STRBUF_INITIAL_SIZE as i32,
            INCHI_STRBUF_SIZE_INCREMENT as i32,
        )? <= 0
        {
            append_api_log(heap, &mut streams[1], "Cannot allocate internal string buffer. Terminating\n")?;
            raw_return = _IS_FATAL as i32;
            break 'source;
        }

        if structure_data.bUserQuit == 0 && interrupted == 0 {
            if input_parameters.last_struct_number != 0
                && input_number >= input_parameters.last_struct_number
            {
                raw_return = _IS_EOF as i32;
                break 'source;
            }

            let (out_stream, remaining) = streams.split_at_mut(1);
            let log_stream = &mut remaining[0];
            raw_return = heap.with_slice_mut_and_heap_mut(original, |records, heap| {
                ExtractOneStructure(
                    heap,
                    &mut structure_data,
                    &mut input_parameters,
                    Some(title),
                    extended_input,
                    Some(log_stream),
                    Some(&mut out_stream[0]),
                    Some(&mut problem_stream),
                    records.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?,
                    &mut input_number,
                )
            })?;
            if matches!(raw_return, value if value == _IS_FATAL as i32 || value == _IS_EOF as i32 || value == _IS_ERROR as i32) {
                break 'source;
            }

            let stream_records = heap.allocate_model_storage(streams.to_vec())?;
            let problem_record = heap.allocate_model_storage(vec![problem_stream.clone()])?;
            let process_result = heap.with_slice_mut_and_heap_mut(title, |title, heap| {
                ProcessOneStructureEx(
                    heap,
                    clock,
                    canonical,
                    &mut structure_data,
                    &mut input_parameters,
                    Some(title),
                    &mut inchi_components,
                    &mut aux_components,
                    SourceMutPointer::null(),
                    stream_records.offset(1)?,
                    stream_records,
                    problem_record,
                    original,
                    prepared,
                    input_number,
                    Some(&mut string_buffer),
                    0,
                    stdout,
                    clock_result,
                )
            });
            streams.clone_from_slice(heap.slice(stream_records.as_const())?);
            problem_stream = heap.slice(problem_record.as_const())?[0].clone();
            heap.free(stream_records)?;
            heap.free(problem_record)?;
            let process_return = process_result?;

            FreeAllINChIArrays(
                heap,
                &mut inchi_components,
                &mut aux_components,
                &mut structure_data.num_components,
            )?;
            heap.with_slice_mut_and_heap_mut(original, |records, heap| {
                FreeOrigAtData(heap, records.first_mut())
            })?;
            heap.with_slice_mut_and_heap_mut(prepared, |records, heap| {
                for record in records.iter_mut().take(2) {
                    FreeOrigAtData(heap, Some(record))?;
                }
                Ok(())
            })?;
            raw_return = raw_return.max(process_return);
        }
    }

    if !translate_directly {
        FreeAllINChIArrays(
            heap,
            &mut inchi_components,
            &mut aux_components,
            &mut structure_data.num_components,
        )?;
        heap.with_slice_mut_and_heap_mut(original, |records, heap| {
            FreeOrigAtData(heap, records.first_mut())
        })?;
        heap.with_slice_mut_and_heap_mut(prepared, |records, heap| {
            for record in records.iter_mut().take(2) {
                FreeOrigAtData(heap, Some(record))?;
            }
            Ok(())
        })?;
        inchi_strbuf_close(heap, Some(&mut string_buffer))?;
        for path in input_parameters.path.iter_mut().take(MAX_NUM_PATHS as usize) {
            if !path.is_null() {
                inchi_free(heap, path.as_mut())?;
                *path = SourceConstPointer::null();
            }
        }
        let mut globals = heap.slice(canonical.as_const())?[0].clone();
        SetBitFree(heap, &mut globals)?;
        heap.slice_mut(canonical)?[0] = globals;
        let (out_stream, remaining_streams) = streams.split_at_mut(1);
        produce_generation_output(
            heap,
            output.as_deref_mut(),
            &structure_data,
            &input_parameters,
            &mut remaining_streams[0],
            &mut out_stream[0],
        )?;
    }

    inchi_ios_close(heap, Some(&mut streams[1]))?;
    inchi_ios_close(heap, Some(&mut streams[0]))?;
    inchi_ios_close(heap, Some(&mut problem_stream))?;
    heap.free(clock)?;
    heap.free(canonical)?;
    heap.free(title)?;
    heap.free(sdf_value)?;
    heap.free(empty_argument)?;
    heap.free(original)?;
    heap.free(prepared)?;

    Ok(match raw_return {
        value if value == _IS_SKIP as i32 => tagRetValGetINCHI_inchi_Ret_SKIP,
        value if value == _IS_EOF as i32 => tagRetValGetINCHI_inchi_Ret_EOF,
        value if value == _IS_OKAY as i32 => tagRetValGetINCHI_inchi_Ret_OKAY,
        value if value == _IS_WARNING as i32 => tagRetValGetINCHI_inchi_Ret_WARNING,
        value if value == _IS_ERROR as i32 => tagRetValGetINCHI_inchi_Ret_ERROR,
        value if value == _IS_FATAL as i32 => tagRetValGetINCHI_inchi_Ret_FATAL,
        value if value == _IS_UNKNOWN as i32 => tagRetValGetINCHI_inchi_Ret_UNKNOWN,
        _ => tagRetValGetINCHI_inchi_Ret_UNKNOWN,
    })
}

#[allow(non_snake_case)]
pub(crate) fn produce_generation_output(
    heap: &mut SourceHeap,
    output: Option<&mut inchi_Output>,
    structure_data: &STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    log_file: &mut INCHI_IOSTREAM,
    output_file: &mut INCHI_IOSTREAM,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:741 produce_generation_output
    // INCHI✔️❌: void produce_generation_output( inchi_Output *out,
    // INCHI✔️❌:                                 STRUCT_DATA *sd,
    // INCHI✔️❌:                                 INPUT_PARMS *ip,
    // INCHI✔️❌:                                 INCHI_IOSTREAM *log_file,
    // INCHI✔️❌:                                 INCHI_IOSTREAM *out_file )
    // INCHI✔️❌:
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (sd->pStrErrStruct[0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (out && ( out->szMessage = (char *) inchi_malloc( strlen( sd->pStrErrStruct ) + 1 ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             strcpy(out->szMessage, sd->pStrErrStruct);
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Make separate strings with InChI and AuxInfo */
    // INCHI✔️❌:     if (out_file->s.pStr && out_file->s.nUsedLength > 0 && out)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         char *p;
    // INCHI✔️❌:         out->szInChI = out_file->s.pStr;
    // INCHI✔️❌:         out->szAuxInfo = NULL;
    // INCHI✔️❌:         if (!( INCHI_OUT_SDFILE_ONLY & ip->bINChIOutputOptions )) /* do not remove last LF from SDF output - 2008-12-23 DT */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (p = strchr( out->szInChI, '\n' ); p; p = strchr( p + 1, '\n' ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!memcmp( p, "\nAuxInfo", 8 ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *p = '\0';            /* remove LF after INChI */
    // INCHI✔️❌:                     out->szAuxInfo = p + 1; /* save pointer to AuxInfo */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else if (out->szAuxInfo || !p[1])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* remove LF after aux info or from the last char */
    // INCHI✔️❌:                     *p = '\0';
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         out_file->s.pStr = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     copy_corrected_log_tail( out, log_file );
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: produce_generation_output

    let mut output = output;
    if structure_data.pStrErrStruct[0] != 0 {
        if let Some(output) = output.as_deref_mut() {
            let message_length = structure_data
                .pStrErrStruct
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            output.szMessage = match inchi_malloc(heap, (message_length + 1) as u64) {
                Ok(message) => {
                    heap.slice_mut(message)?[..=message_length]
                        .copy_from_slice(&structure_data.pStrErrStruct[..=message_length]);
                    message
                }
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                Err(error) => return Err(error),
            };
        }
    }

    if !output_file.s.pStr.is_null() && output_file.s.nUsedLength > 0 {
        if let Some(output) = output.as_deref_mut() {
            output.szInChI = output_file.s.pStr;
            output.szAuxInfo = SourceMutPointer::null();
            if input_parameters.bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY as i32 == 0 {
                let bytes = heap.slice_mut(output.szInChI)?;
                let nul = bytes
                    .iter()
                    .position(|byte| *byte == 0)
                    .ok_or(SourceHeapError::MissingNulTerminator)?;
                let mut search_start = 0_usize;
                while search_start < nul {
                    let Some(relative) = bytes[search_start..nul]
                        .iter()
                        .position(|byte| *byte == b'\n' as i8)
                    else {
                        break;
                    };
                    let position = search_start + relative;
                    let compared = bytes
                        .get(position..position + 8)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if compared
                        .iter()
                        .map(|byte| *byte as u8)
                        .eq(b"\nAuxInfo".iter().copied())
                    {
                        bytes[position] = 0;
                        output.szAuxInfo = output.szInChI.offset((position + 1) as i64)?;
                    } else if !output.szAuxInfo.is_null() || bytes[position + 1] == 0 {
                        bytes[position] = 0;
                        break;
                    }
                    search_start = position + 1;
                }
            }
            output_file.s.pStr = SourceMutPointer::null();
        }
    }

    copy_corrected_log_tail(heap, output, log_file)
}

pub(crate) fn copy_corrected_log_tail(
    heap: &mut SourceHeap,
    output: Option<&mut inchi_Output>,
    log_file: &mut INCHI_IOSTREAM,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:787 copy_corrected_log_tail
    // INCHI✔️❌: void copy_corrected_log_tail( inchi_Output *out, INCHI_IOSTREAM *log_file )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (log_file->s.pStr && log_file->s.nUsedLength > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         while (log_file->s.nUsedLength &&
    // INCHI✔️❌:                 '\n' == log_file->s.pStr[log_file->s.nUsedLength - 1])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             log_file->s.pStr[--log_file->s.nUsedLength] = '\0';
    // INCHI✔️❌:                                             /* remove last LF */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (out)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             char *p;
    // INCHI✔️❌:             out->szLog = log_file->s.pStr;
    // INCHI✔️❌:             log_file->s.pStr = NULL;
    // INCHI✔️❌:             for (p = strchr( out->szLog, ' ' ); p; p = strchr( p + 1, ' ' ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!memcmp( p, " structure #", 12 ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *p = '\0';
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: copy_corrected_log_tail

    if log_file.s.pStr.is_null() || log_file.s.nUsedLength <= 0 {
        return Ok(());
    }

    while log_file.s.nUsedLength != 0 {
        let last = usize::try_from(log_file.s.nUsedLength - 1)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let bytes = heap.slice(log_file.s.pStr.as_const())?;
        let byte = *bytes.get(last).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if byte != b'\n' as i8 {
            break;
        }
        log_file.s.nUsedLength -= 1;
        heap.slice_mut(log_file.s.pStr)?[last] = 0;
    }

    if let Some(output) = output {
        output.szLog = log_file.s.pStr;
        log_file.s.pStr = SourceMutPointer::null();
        let bytes = heap.slice_mut(output.szLog)?;
        let nul = bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let mut search_start = 0_usize;
        while search_start < nul {
            let Some(relative) = bytes[search_start..nul]
                .iter()
                .position(|byte| *byte == b' ' as i8)
            else {
                break;
            };
            let position = search_start + relative;
            let compared = bytes
                .get(position..position + 12)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if compared
                .iter()
                .map(|byte| *byte as u8)
                .eq(b" structure #".iter().copied())
            {
                bytes[position] = 0;
            }
            search_start = position + 1;
        }
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CheckINCHI(
    heap: &mut SourceHeap,
    szINCHI: SourceConstPointer<i8>,
    strict: i32,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
    clock_result: clock_t,
) -> Result<u32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:834 CheckINCHI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL CheckINCHI( const char *szINCHI, const int strict )
{
    int ret = INCHI_VALID_NON_STANDARD;
    int ret_i2i;
    inchi_InputINCHI    inchi_inp;
    inchi_Output        inchi_out;
    size_t slen, pos_slash1 = 0;
    char *str = NULL;
    size_t i;
    size_t slen0;
    char pp;

    /* .. non-empty */
    if (szINCHI == NULL)
    {
        return INCHI_INVALID_PREFIX;
    }

    slen = strlen( szINCHI );


    /* .. has valid prefix */
    if (slen < LEN_INCHI_STRING_PREFIX + 3)
    {
        return INCHI_INVALID_PREFIX;
    }
    if (memcmp( szINCHI, INCHI_STRING_PREFIX, LEN_INCHI_STRING_PREFIX ))
    {
        return INCHI_INVALID_PREFIX;
    }

    /* .. has InChI version 1 */
    /* if (!isdigit(szINCHI[LEN_INCHI_STRING_PREFIX]) )  */
    if (szINCHI[LEN_INCHI_STRING_PREFIX] != '1')
    {
        return INCHI_INVALID_VERSION;
    }

    /* .. optionally has a 'standard' flag character */
    pos_slash1 = LEN_INCHI_STRING_PREFIX + 1;
    if (szINCHI[pos_slash1] == 'S')
    {
        /* Standard InChI ==> standard InChIKey */
        ret = INCHI_VALID_STANDARD;
        pos_slash1++;
    }
    else if (szINCHI[pos_slash1] == 'B')
    {
        /* Beta version InChI ==> non-standard */
        ret = INCHI_VALID_BETA;
        pos_slash1++;
    }

    /* .. has trailing slash in the right place */
    if (szINCHI[pos_slash1] != '/')
    {
        return INCHI_INVALID_LAYOUT;
    }

    /* .. the rest of source string contains valid literals */


    /* adjust line len so we not check trailing whitespaces */
    i = slen - 1;
    while (isspace(UCINT szINCHI[i--])) slen--;

    /* Treat possible SaveOpt letters  */
    slen0 = slen;
    if (( szINCHI[slen - 3] == '\\' ) &&
        ( szINCHI[slen - 2] >= 'A' ) && ( szINCHI[slen - 2] <= 'Z' ) &&
        ( szINCHI[slen - 1] >= 'A' ) && ( szINCHI[slen - 1] <= 'Z' )
        )
    {
        slen0 = slen - 3;
    }

    int prev_is_slash = 1;
    for (i = pos_slash1 + 1; i < slen0; i++)
    {
        pp = szINCHI[i];
#if ( FIX_GAF_2020_GENERIC==1 )
        if (prev_is_slash)
        {
            /* After slash: */
            if (pp == '0')
            {
                /* '0' is never allowed */
                return INCHI_INVALID_LAYOUT;
            }
            if (i > pos_slash1 + 1)
            {
                /* Not in main formula layer... */ 
                if (!islower(pp))
                {
                    /* only lowercase letters are allowed */
                    return INCHI_INVALID_LAYOUT;
                }
            }
        }
        prev_is_slash = (pp != '/') ? 0 : 1;
#endif
        if (pp >= 'A' && pp <= 'Z')   continue;
        if (pp >= 'a' && pp <= 'z')   continue;
        if (pp >= '0' && pp <= '9')  continue;
        switch (pp)
        {
            case '(': case ')':
            case '*': case '+':
            case ',': case '-':
            case '.': case '/':
#if ( FIX_GAF_2020_GENERIC==1 )
            case ';': case '?':     continue;
#else
            case ';': case '=':
            case '?': case '@':     continue;
#endif
            default:            return INCHI_INVALID_LAYOUT;
        }
    }

    if (strict)
    {
        char opts[] = "?FixedH ?RecMet ?SUU ?SLUUD";
        extract_inchi_substring( &str, szINCHI, slen );
        if (NULL == str)
        {
            ret = INCHI_FAIL_I2I;
            goto fin;
        }

        inchi_inp.szInChI = str;
        opts[0] = opts[8] = opts[16] = opts[21] = INCHI_OPTION_PREFX;
        inchi_inp.szOptions = opts;

        ret_i2i = GetINCHIfromINCHI( &inchi_inp, &inchi_out );

        if (( ( ret_i2i != inchi_Ret_OKAY ) && ( ret_i2i != inchi_Ret_WARNING ) ) || !inchi_out.szInChI)
        {
            ret = INCHI_FAIL_I2I;
        }
        else
        {
            if (strcmp( inchi_inp.szInChI, inchi_out.szInChI ))
            {
                ret = INCHI_FAIL_I2I;
            }
        }
    }

fin:if (strict)
{
    if (NULL != str)
        inchi_free( str );
}

    return ret;
}
        */
    // END INCHI C FUNCTION: CheckINCHI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CheckINCHI
    // INCHI✔️❌: TARGET_API_LIB and COMPILE_ANSI_ONLY are defined for GCC/Linux.
    // INCHI✔️❌: FIX_GAF_2020_GENERIC=1 selects slash-layer validation and excludes '=' and '@' from accepted body literals.
    // INCHI✔️❌: UCINT casts through unsigned char and the selected C locale supplies ASCII isspace/islower behavior.
    // INCHI✔️❌: strict validation calls GetINCHIfromINCHI, which calls back only with strict=0; this is the complete active SCC.
    // INCHI✔️❌: SourceHeap snapshots and checked pointer access add work absent from C pointer arithmetic.
    // END INCHI ACTIVE MACRO CONFIGURATION: CheckINCHI

    const PREFIX: &[u8] = b"InChI=";
    if szINCHI.is_null() {
        return Ok(tagRetValCheckINCHI_INCHI_INVALID_PREFIX);
    }
    let bytes = heap.slice(szINCHI)?;
    let mut length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    if length < PREFIX.len() + 3 {
        return Ok(tagRetValCheckINCHI_INCHI_INVALID_PREFIX);
    }
    if !bytes[..PREFIX.len()]
        .iter()
        .map(|byte| *byte as u8)
        .eq(PREFIX.iter().copied())
    {
        return Ok(tagRetValCheckINCHI_INCHI_INVALID_PREFIX);
    }
    if bytes[PREFIX.len()] as u8 != b'1' {
        return Ok(tagRetValCheckINCHI_INCHI_INVALID_VERSION);
    }

    let mut result = tagRetValCheckINCHI_INCHI_VALID_NON_STANDARD;
    let mut first_slash = PREFIX.len() + 1;
    if bytes[first_slash] as u8 == b'S' {
        result = tagRetValCheckINCHI_INCHI_VALID_STANDARD;
        first_slash += 1;
    } else if bytes[first_slash] as u8 == b'B' {
        result = tagRetValCheckINCHI_INCHI_VALID_BETA;
        first_slash += 1;
    }
    if bytes[first_slash] as u8 != b'/' {
        return Ok(tagRetValCheckINCHI_INCHI_INVALID_LAYOUT);
    }

    while length > 0
        && matches!(
            bytes[length - 1] as u8,
            b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'
        )
    {
        length -= 1;
    }

    let mut body_length = length;
    if length >= 3
        && bytes[length - 3] as u8 == b'\\'
        && (bytes[length - 2] as u8).is_ascii_uppercase()
        && (bytes[length - 1] as u8).is_ascii_uppercase()
    {
        body_length = length - 3;
    }

    let mut previous_is_slash = true;
    for index in first_slash + 1..body_length {
        let byte = bytes[index] as u8;
        if previous_is_slash {
            if byte == b'0' {
                return Ok(tagRetValCheckINCHI_INCHI_INVALID_LAYOUT);
            }
            if index > first_slash + 1 && !byte.is_ascii_lowercase() {
                return Ok(tagRetValCheckINCHI_INCHI_INVALID_LAYOUT);
            }
        }
        previous_is_slash = byte == b'/';
        if byte.is_ascii_alphanumeric()
            || matches!(
                byte,
                b'(' | b')' | b'*' | b'+' | b',' | b'-' | b'.' | b'/' | b';' | b'?'
            )
        {
            continue;
        }
        return Ok(tagRetValCheckINCHI_INCHI_INVALID_LAYOUT);
    }
    if strict == 0 {
        return Ok(result);
    }

    let mut extracted = SourceMutPointer::null();
    match extract_inchi_substring(heap, &mut extracted, szINCHI, length as u64) {
        Ok(()) => {}
        Err(SourceHeapError::AllocationFailed) => {
            return Ok(tagRetValCheckINCHI_INCHI_FAIL_I2I);
        }
        Err(error) => return Err(error),
    }
    if extracted.is_null() {
        return Ok(tagRetValCheckINCHI_INCHI_FAIL_I2I);
    }

    let options = heap.allocate_model_storage(
        b"?FixedH ?RecMet ?SUU ?SLUUD\0"
            .iter()
            .map(|byte| *byte as i8)
            .collect(),
    )?;
    {
        let option_bytes = heap.slice_mut(options)?;
        for index in [0_usize, 8, 16, 21] {
            option_bytes[index] = b'-' as i8;
        }
    }
    let input = inchi_InputINCHI {
        szInChI: extracted,
        szOptions: options,
    };
    let mut output = inchi_Output::default();
    let conversion =
        GetINCHIfromINCHI(heap, Some(&input), &mut output, stdout, build, clock_result);
    heap.free(options)?;
    let conversion = match conversion {
        Ok(value) => value,
        Err(error) => {
            inchi_free(heap, extracted)?;
            return Err(error);
        }
    };
    let extracted_bytes = heap.slice(extracted.as_const())?
        [..source_strlen(heap, extracted.as_const())?]
        .iter()
        .map(|byte| *byte as u8)
        .collect::<Vec<_>>();
    if !matches!(
        conversion,
        crate::source_types::tagRetValGetINCHI_inchi_Ret_OKAY
            | crate::source_types::tagRetValGetINCHI_inchi_Ret_WARNING
    ) || output.szInChI.is_null()
        || !source_string_equals(heap, output.szInChI.as_const(), &extracted_bytes)?
    {
        result = tagRetValCheckINCHI_INCHI_FAIL_I2I;
    }

    inchi_free(heap, extracted)?;
    Ok(result)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetINCHIfromINCHI(
    heap: &mut SourceHeap,
    inpInChI: Option<&inchi_InputINCHI>,
    out: &mut inchi_Output,
    stdout: SourceMutPointer<FILE>,
    build: InchiBuildMetadata<'_>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2114 GetINCHIfromINCHI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int INCHI_DECL GetINCHIfromINCHI( inchi_InputINCHI *inpInChI,
                                  inchi_Output *out )
{
    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;

    static char szMainOption[] = " ?InChI2InChI";

    INCHI_CLOCK ic;
    CANON_GLOBALS CG;

    int i;
    char      szSdfDataValue[MAX_SDF_VALUE + 1];
    unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */

    INPUT_PARMS inp_parms;
    INPUT_PARMS *ip = &inp_parms;

    int             bReleaseVersion = bRELEASE_VERSION;
    int   nRet = 0, nRet1;

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    int  num_repeat = REPEAT_ALL;
#endif

    const char *argv[INCHI_MAX_NUM_ARG + 1];
    int   argc;
    char *szOptions = NULL;

    INCHI_IOSTREAM inchi_file[3], *out_file = inchi_file, *log_file = inchi_file + 1, *input_file = inchi_file + 2;

#if( TRACE_MEMORY_LEAKS == 1 )
    _CrtSetDbgFlag( _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF );
/* for execution outside the VC++ debugger uncomment one of the following two */

#ifdef MY_REPORT_FILE
    _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode( _CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG );
#endif

    /* turn on floating point exceptions */
#if ( !defined(__STDC__) || __STDC__ != 1 )
    {
        /* Get the default control word. */
        int cw = _controlfp( 0, 0 );

        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &= ~( EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL );

        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif
#endif

    memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */
#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
repeat:
    FreeINCHI( out );
    inchi_ios_close( out_file );
    inchi_ios_close( log_file );
    inchi_ios_reset( input_file );  /* do not close input_file - its string buffer may point to inpInChI->szInChI */
#endif

    /* Initialize internal for this function I/O streams as string buffers */
    inchi_ios_init( input_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( out_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( log_file, INCHI_IOS_TYPE_STRING, NULL );

    sd->bUserQuit = 0;

    /* clear original input structure */
    /* memset( inchi_file, 0, sizeof(inchi_file) ); */
    memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( ip, 0, sizeof( *ip ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( szSdfDataValue, 0, sizeof( szSdfDataValue ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    memset( &ic, 0, sizeof( ic ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( &CG, 0, sizeof( CG ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    szMainOption[1] = INCHI_OPTION_PREFX;

    if (!inpInChI)
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }

    /* options */
    if (inpInChI)
    {
        int opt_len = (int) ( ( inpInChI->szOptions ? strlen( inpInChI->szOptions ) : 0 ) + sizeof( szMainOption ) + 1 );
        szOptions = (char*) inchi_calloc( (long long)opt_len + 1, sizeof( szOptions[0] ) ); /* djb-rwth: cast operator added */
        if (szOptions)
        {
            if (inpInChI->szOptions)
            {
                strcpy( szOptions, inpInChI->szOptions );
            }
            strcat( szOptions, szMainOption );
            argc = parse_options_string( szOptions, argv, INCHI_MAX_NUM_ARG );
        }
        else
        {
            nRet = _IS_FATAL;
            goto translate_RetVal; /* emergency exit */
        }
    }
    else
    {
        argc = 1;
        argv[0] = "";
        argv[1] = NULL;
    }

    if ((argc == 1
#ifdef TARGET_API_LIB
        && ( !inpInChI || !inpInChI->szInChI ))
#endif

        || (argc == 2 && ( argv[1][0] == INCHI_OPTION_PREFX ) &&
        ( !strcmp( argv[1] + 1, "?" ) || !inchi_stricmp( argv[1] + 1, "help" ) ))) /* djb-rwth: addressing LLVM warnings */
    {
        HelpCommandLineParms( log_file );
        out->szLog = log_file->s.pStr;
        memset( log_file, 0, sizeof( *log_file ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        nRet = _IS_EOF;
        inchi_free(szOptions); /* djb-rwth: avoiding memory leak */
        goto translate_RetVal;
    }

    nRet1 = ReadCommandLineParms( argc, argv, ip, szSdfDataValue,
                                &ulDisplTime, bReleaseVersion, log_file );

    if (szOptions)
    {
        /* argv pointed to strings in szOptions */
        inchi_free( szOptions );
        szOptions = NULL;
    }
    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if (0 > nRet1)
    {
        goto exit_function;
    }

    if (ip->bNoStructLabels)
    {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    }
    else if (ip->nInputType == INPUT_INCHI_XML ||
            ip->nInputType == INPUT_INCHI_PLAIN ||
            ip->nInputType == INPUT_CMLFILE ||
            ip->nInputType == INPUT_INCHI)
    {
        /* the input may contain both the header and the label of the structure */
        if (!ip->pSdfLabel)
        {
            ip->pSdfLabel = ip->szSdfDataHeader;
        }
        if (!ip->pSdfValue)
        {
            ip->pSdfValue = szSdfDataValue;
        }
    }

    if (ip->nInputType && ip->nInputType != INPUT_INCHI)
    {
        inchi_ios_eprint( log_file, "Input type set to INPUT_INCHI\n" );
        ip->nInputType = INPUT_INCHI;
    }

    if (!inpInChI->szInChI)
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }
    else
    {
        const int strict = 0;
        nRet = CheckINCHI( inpInChI->szInChI, strict );
        if (nRet != INCHI_VALID_STANDARD     &&
            nRet != INCHI_VALID_NON_STANDARD &&
            nRet != INCHI_VALID_BETA)
        {
            nRet = _IS_ERROR;
            goto exit_function;
        }
    }


    PrintInputParms( log_file, ip );

    /********************************/
    /* InChI -> InChI               */
    /********************************/

    /* input_file simulation */
    input_file->s.pStr = inpInChI->szInChI;
    input_file->s.nUsedLength = (int) strlen( input_file->s.pStr ) + 1;
    input_file->s.nAllocatedLength = input_file->s.nUsedLength;
    input_file->s.nPtr = 0;

    /* buffer for the message */
    out->szMessage = (char *) inchi_calloc( MAX_MSG_LEN, sizeof( out->szMessage[0] ) );
    if (!out->szMessage)
    {
        inchi_ios_eprint( log_file, "Cannot allocate output message buffer.\n" );
        nRet = -1;
    }
    else
    {
        nRet = ReadWriteInChI( &ic, &CG, input_file, out_file, log_file,
                                ip, sd,
                                NULL, 0, NULL,
                                NULL, NULL,
                                out->szMessage, MAX_MSG_LEN,
                                NULL /*out->WarningFlags*/ );
    }

    if (nRet >= 0 && out_file->s.pStr)
    {
        /* success */
        char* p = NULL;
        /* djb-rwth: fixing oss-fuzz issue #40971 */
        int p_len, out_szinchi_len;
        out_szinchi_len = strlen(out_file->s.pStr);
        out->szInChI = out_file->s.pStr;
        out->szAuxInfo = NULL;

        for (p = strchr( out->szInChI, '\n' ); p; p = strchr( p + 1, '\n' ))
        {
            p_len = strlen(p);
            if ((p_len >= 8) && !memcmp( p, "\nAuxInfo", 8 ))
            {
                *p = '\0';            /* remove LF after INChI */
                out->szAuxInfo = p + 1; /* save pointer to AuxInfo */
            }
            else if (out->szAuxInfo || !p[1])
            {
                /* remove LF after aux info or from the last char */
                *p = '\0';
                break;
            }
        }
        out_file->s.pStr = NULL;
    }

    /*
    out->szLog = log_file->pStr;
    log_file->pStr   = NULL;
    */

exit_function:;

    for (i = 0; i < MAX_NUM_PATHS; i++)
    {
        if (ip->path[i])
        {
            inchi_free( (char*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }

    SetBitFree( &CG );

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    if (num_repeat-- > 0)
    {
        goto repeat;
    }
#endif

#ifdef TARGET_API_LIB
    /* output */

        if (log_file->s.pStr && log_file->s.nUsedLength > 0)
        {
            while (log_file->s.nUsedLength && '\n' == log_file->s.pStr[log_file->s.nUsedLength - 1])
            {
                log_file->s.pStr[--log_file->s.nUsedLength] = '\0'; /* remove last LF */
            }
            if (out)
            {
                out->szLog = log_file->s.pStr;
                log_file->s.pStr = NULL;
            }
        }

#endif

translate_RetVal:

    /* close internal output streams */
    inchi_ios_close( out_file );
    inchi_ios_close( log_file );
    inchi_ios_reset( input_file );  /* do not close input_file - its string buffer may point to inpInChI->szInChI */

    switch (nRet)
    {
        case -3: nRet = inchi_Ret_ERROR; break; /* Error: no Structure has been created */
        case -2: nRet = inchi_Ret_ERROR; break; /* Error: no Structure has been created */
        case -1: nRet = inchi_Ret_FATAL; break; /* Severe error: no Structure has been created (typically; break; memory allocation failed) */
        default:
            /*
            if ( !outStruct->atom || !outStruct->num_atoms )
            {
                nRet = inchi_Ret_EOF;
            }
            else
            {
                int m,n,t=0;
                for ( m=0; m < 2; m ++ )
                {
                    for ( n=0; n < 2; n ++ )
                    {
                        if ( outStruct->WarningFlags[m][n] ) {
                            t ++;
                        }
                    }
                }
                nRet = t? inchi_Ret_WARNING : inchi_Ret_OKAY;
            }
            */
            break;
    }

    return nRet;
}
        */
    // END INCHI C FUNCTION: GetINCHIfromINCHI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetINCHIfromINCHI
    // INCHI✔️❌: TARGET_API_LIB and COMPILE_ANSI_ONLY are defined for GCC/Linux; TARGET_EXE_STANDALONE is not defined.
    // INCHI✔️❌: TRACE_MEMORY_LEAKS=0 and REPEAT_ALL is undefined, excluding debugger and repeat branches.
    // INCHI✔️❌: bRELEASE_VERSION=1, INCHI_OPTION_PREFX='-', INCHI_MAX_NUM_ARG=32, MAX_MSG_LEN=512.
    // INCHI✔️❌: TARGET_API_LIB transfers a nonempty log buffer after stripping all trailing LF bytes.
    // INCHI✔️❌: SourceHeap stream buffers and checked temporary snapshots add allocation and copying overhead.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetINCHIfromINCHI

    const MAX_MESSAGE_LENGTH: usize = 512;
    const MAIN_OPTION: &[u8] = b" -InChI2InChI\0";

    *out = inchi_Output::default();

    let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()])?;
    let canon = heap.allocate_model_storage(vec![CANON_GLOBALS::default()])?;
    let sdf_value = heap.allocate_model_storage(vec![0_i8; MAX_SDF_VALUE as usize + 1])?;
    let empty_argument = heap.allocate_model_storage(vec![0_i8])?;

    let mut streams: [INCHI_IOSTREAM; 3] = std::array::from_fn(|_| INCHI_IOSTREAM::default());
    for stream in &mut streams {
        inchi_ios_init(
            Some(stream),
            INCHI_IOS_TYPE_STRING as i32,
            SourceMutPointer::null(),
        )?;
    }

    let mut structure_data = STRUCT_DATA::default();
    let mut input_parameters = INPUT_PARMS::default();
    let mut display_time = 0_u64;
    let mut raw_return = 0_i32;
    let mut translate_directly = false;

    'source: {
        let Some(input) = inpInChI else {
            raw_return = _IS_ERROR as i32;
            break 'source;
        };

        let option_length = if input.szOptions.is_null() {
            0
        } else {
            source_strlen(heap, input.szOptions.as_const())?
        };
        let source_option_length = option_length
            .checked_add(MAIN_OPTION.len())
            .and_then(|value| value.checked_add(2))
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        let options = match inchi_calloc::<i8>(heap, source_option_length as u64, 1) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                raw_return = _IS_FATAL as i32;
                translate_directly = true;
                break 'source;
            }
            Err(error) => return Err(error),
        };
        if option_length != 0 {
            let source = heap.slice(input.szOptions.as_const())?[..option_length].to_vec();
            heap.slice_mut(options)?[..option_length].copy_from_slice(&source);
        }
        let main = MAIN_OPTION
            .iter()
            .map(|byte| *byte as i8)
            .collect::<Vec<_>>();
        heap.slice_mut(options)?[option_length..option_length + MAIN_OPTION.len()]
            .copy_from_slice(&main);

        let mut parsed_arguments = vec![SourceArgvPointer::Null; INCHI_MAX_NUM_ARG as usize + 1];
        let argument_count = parse_options_string(
            heap,
            options,
            &mut parsed_arguments,
            INCHI_MAX_NUM_ARG as i32,
        )?;
        let mut argv = vec![SourceConstPointer::null(); argument_count as usize];
        for (index, argument) in parsed_arguments[..argument_count as usize]
            .iter()
            .copied()
            .enumerate()
        {
            argv[index] = match argument {
                SourceArgvPointer::EmptyLiteral => empty_argument.as_const(),
                SourceArgvPointer::Command(pointer) => pointer.as_const(),
                SourceArgvPointer::Null => return Err(SourceHeapError::NullPointer),
            };
        }

        let help_requested = if argument_count == 1 && input.szInChI.is_null() {
            true
        } else if argument_count == 2 {
            let argument = argv[1];
            let bytes = heap.slice(argument)?;
            bytes.first().copied() == Some(b'-' as i8)
                && (source_string_equals(heap, argument.offset(1)?, b"?")? || {
                    let help = heap.allocate_model_storage(
                        b"help\0".iter().map(|byte| *byte as i8).collect(),
                    )?;
                    let equal = inchi_stricmp(heap, argument.offset(1)?, help.as_const())? == 0;
                    heap.free(help)?;
                    equal
                })
        } else {
            false
        };
        if help_requested {
            HelpCommandLineParms(heap, Some(&mut streams[1]), stdout, build)?;
            out.szLog = streams[1].s.pStr;
            streams[1] = INCHI_IOSTREAM::default();
            raw_return = _IS_EOF;
            inchi_free(heap, options)?;
            translate_directly = true;
            break 'source;
        }

        let parameter_return = ReadCommandLineParms(
            heap,
            argument_count,
            &argv,
            &mut input_parameters,
            sdf_value,
            &mut display_time,
            1,
            Some(&mut streams[1]),
        )?;
        inchi_free(heap, options)?;

        input_parameters.bNoStructLabels = 1;
        if parameter_return < 0 {
            break 'source;
        }
        input_parameters.pSdfLabel = SourceMutPointer::null();
        input_parameters.pSdfValue = SourceMutPointer::null();

        if input_parameters.nInputType != 0
            && input_parameters.nInputType != tagInputType_INPUT_INCHI
        {
            append_api_log(heap, &mut streams[1], "Input type set to INPUT_INCHI\n")?;
            input_parameters.nInputType = tagInputType_INPUT_INCHI;
        }

        if input.szInChI.is_null() {
            raw_return = _IS_ERROR as i32;
            break 'source;
        }
        let validity = CheckINCHI(
            heap,
            input.szInChI.as_const(),
            0,
            stdout,
            build,
            clock_result,
        )?;
        if !matches!(
            validity,
            tagRetValCheckINCHI_INCHI_VALID_STANDARD
                | tagRetValCheckINCHI_INCHI_VALID_NON_STANDARD
                | tagRetValCheckINCHI_INCHI_VALID_BETA
        ) {
            raw_return = _IS_ERROR as i32;
            break 'source;
        }

        PrintInputParms(heap, Some(&mut streams[1]), &input_parameters)?;

        let input_length = source_strlen(heap, input.szInChI.as_const())?;
        streams[2].s.pStr = input.szInChI;
        streams[2].s.nUsedLength =
            i32::try_from(input_length + 1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        streams[2].s.nAllocatedLength = streams[2].s.nUsedLength;
        streams[2].s.nPtr = 0;

        out.szMessage = match inchi_calloc::<i8>(heap, MAX_MESSAGE_LENGTH as u64, 1) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                append_api_log(
                    heap,
                    &mut streams[1],
                    "Cannot allocate output message buffer.\n",
                )?;
                raw_return = -1;
                SourceMutPointer::null()
            }
            Err(error) => return Err(error),
        };

        if !out.szMessage.is_null() {
            let message = out.szMessage;
            let (output_stream, remaining_streams) = streams.split_at_mut(1);
            let (log_stream, input_stream) = remaining_streams.split_at_mut(1);
            raw_return = heap.with_slice_mut_and_heap_mut(message, |message, heap| {
                ReadWriteInChI(
                    heap,
                    clock,
                    canon,
                    &mut input_stream[0],
                    &mut output_stream[0],
                    &mut log_stream[0],
                    &mut input_parameters,
                    &mut structure_data,
                    None,
                    None,
                    None,
                    None,
                    None,
                    Some(message),
                    MAX_MESSAGE_LENGTH as i32,
                    None,
                    stdout,
                    clock_result,
                )
            })?;
        }

        if raw_return >= 0 && !streams[0].s.pStr.is_null() {
            let output = streams[0].s.pStr;
            let mut aux_info = SourceMutPointer::null();
            let mut search = 0_usize;
            loop {
                let nul = source_strlen(heap, output.as_const())?;
                let next = heap.slice(output.as_const())?[search..nul]
                    .iter()
                    .position(|byte| *byte == b'\n' as i8)
                    .map(|offset| search + offset);
                let Some(position) = next else {
                    break;
                };
                let suffix_length = nul - position;
                let is_aux = suffix_length >= 8
                    && heap.slice(output.as_const())?[position..position + 8]
                        == [
                            b'\n' as i8,
                            b'A' as i8,
                            b'u' as i8,
                            b'x' as i8,
                            b'I' as i8,
                            b'n' as i8,
                            b'f' as i8,
                            b'o' as i8,
                        ];
                if is_aux {
                    heap.slice_mut(output)?[position] = 0;
                    aux_info = output.offset((position + 1) as i64)?;
                    search = position + 1;
                } else {
                    let next_is_nul = heap.slice(output.as_const())?.get(position + 1) == Some(&0);
                    if !aux_info.is_null() || next_is_nul {
                        heap.slice_mut(output)?[position] = 0;
                        break;
                    }
                    search = position + 1;
                }
            }
            out.szInChI = output;
            out.szAuxInfo = aux_info;
            streams[0].s.pStr = SourceMutPointer::null();
        }
    }

    if !translate_directly {
        for index in 0..MAX_NUM_PATHS as usize {
            if !input_parameters.path[index].is_null() {
                inchi_free(heap, input_parameters.path[index].as_mut())?;
                input_parameters.path[index] = SourceConstPointer::null();
            }
        }
        let mut globals = heap.slice(canon.as_const())?[0].clone();
        SetBitFree(heap, &mut globals)?;
        heap.slice_mut(canon)?[0] = globals;

        if !streams[1].s.pStr.is_null() && streams[1].s.nUsedLength > 0 {
            while streams[1].s.nUsedLength > 0 {
                let last = usize::try_from(streams[1].s.nUsedLength - 1)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if heap.slice(streams[1].s.pStr.as_const())?[last] != b'\n' as i8 {
                    break;
                }
                streams[1].s.nUsedLength -= 1;
                heap.slice_mut(streams[1].s.pStr)?[last] = 0;
            }
            out.szLog = streams[1].s.pStr;
            streams[1].s.pStr = SourceMutPointer::null();
        }
    }

    inchi_ios_close(heap, Some(&mut streams[0]))?;
    inchi_ios_close(heap, Some(&mut streams[1]))?;
    inchi_ios_reset(heap, &mut streams[2])?;

    raw_return = match raw_return {
        -3 | -2 => tagRetValGetINCHI_inchi_Ret_ERROR,
        -1 => tagRetValGetINCHI_inchi_Ret_FATAL,
        value => value,
    };

    heap.free(sdf_value)?;
    heap.free(empty_argument)?;
    heap.free(clock)?;
    heap.free(canon)?;
    Ok(raw_return)
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
    fn source_port__inchi_dll__getinchi1__line_345() {
        fn c_string(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .bytes()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect(),
            )
            .unwrap()
        }

        fn text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Option<String> {
            if pointer.is_null() {
                return None;
            }
            let bytes = heap.slice(pointer.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            Some(String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap())
        }

        fn methane_input(heap: &mut SourceHeap, options: SourceMutPointer<i8>) -> inchi_InputEx {
            let atoms = heap
                .allocate_model_storage(vec![inchi_Atom {
                    elname: source_name("C"),
                    num_iso_H: [-1, 0, 0, 0],
                    ..inchi_Atom::default()
                }])
                .unwrap();
            inchi_InputEx {
                atom: atoms,
                num_atoms: 1,
                szOptions: options,
                ..inchi_InputEx::default()
            }
        }

        fn invoke(
            heap: &mut SourceHeap,
            input: &inchi_InputEx,
            output: Option<&mut inchi_Output>,
            enforce_standard: i32,
            interrupted: i32,
        ) -> Result<i32, SourceHeapError> {
            GetINCHI1(
                heap,
                Some(input),
                output,
                enforce_standard,
                interrupted,
                SourceMutPointer::null(),
                InchiBuildMetadata {
                    compiler: "gcc",
                    date: "Jan 01 1970",
                    time: "00:00:00",
                },
                1_000,
            )
        }

        let mut null_output_heap = SourceHeap::default();
        let null_output_input = methane_input(&mut null_output_heap, SourceMutPointer::null());
        let null_output_snapshot = null_output_input.clone();
        assert_eq!(
            invoke(&mut null_output_heap, &null_output_input, None, 0, 0),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        assert_eq!(null_output_input, null_output_snapshot);
        assert_eq!(null_output_heap.live_source_allocation_count(), 0);

        let mut help_heap = SourceHeap::default();
        let empty_input = inchi_InputEx::default();
        let old_inchi = c_string(&mut help_heap, "old inchi");
        let old_aux = c_string(&mut help_heap, "old aux");
        let old_message = c_string(&mut help_heap, "old message");
        let old_log = c_string(&mut help_heap, "old log");
        let mut help_output = inchi_Output {
            szInChI: old_inchi,
            szAuxInfo: old_aux,
            szMessage: old_message,
            szLog: old_log,
        };
        assert_eq!(
            invoke(&mut help_heap, &empty_input, Some(&mut help_output), 0, 0),
            Ok(tagRetValGetINCHI_inchi_Ret_EOF)
        );
        assert!(help_output.szInChI.is_null());
        assert!(help_output.szAuxInfo.is_null());
        assert!(help_output.szMessage.is_null());
        let help_log = text(&help_heap, help_output.szLog).unwrap();
        assert!(help_log.starts_with("InChI version 1"));
        assert!(help_log.contains("Usage:"));
        assert!(help_log.contains("SaveOpt"));
        for old in [old_inchi, old_aux, old_message, old_log] {
            assert!(help_heap.slice(old.as_const()).is_ok());
        }
        FreeINCHI(&mut help_heap, Some(&mut help_output)).unwrap();
        assert_eq!(help_heap.live_source_allocation_count(), 0);

        let mut help_option_heap = SourceHeap::default();
        let help_option = c_string(&mut help_option_heap, "-help");
        let help_option_input = methane_input(&mut help_option_heap, help_option);
        let mut help_option_output = inchi_Output::default();
        assert_eq!(
            invoke(
                &mut help_option_heap,
                &help_option_input,
                Some(&mut help_option_output),
                0,
                0,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_EOF)
        );
        FreeINCHI(&mut help_option_heap, Some(&mut help_option_output)).unwrap();
        assert_eq!(help_option_heap.live_source_allocation_count(), 1);
        assert_eq!(help_option_heap.live_source_allocations_of::<i8>(), 1);

        let mut standard_heap = SourceHeap::default();
        let standard_input = methane_input(&mut standard_heap, SourceMutPointer::null());
        let standard_atoms = standard_heap
            .slice(standard_input.atom.as_const())
            .unwrap()
            .to_vec();
        let mut standard_output = inchi_Output::default();
        assert_eq!(
            invoke(
                &mut standard_heap,
                &standard_input,
                Some(&mut standard_output),
                0,
                0,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert_eq!(
            text(&standard_heap, standard_output.szInChI).as_deref(),
            Some("InChI=1S/CH4/h1H4")
        );
        assert_eq!(
            text(&standard_heap, standard_output.szAuxInfo).as_deref(),
            Some("AuxInfo=1/0/N:1/rA:1C/rB:/rC:;")
        );
        assert!(standard_output.szMessage.is_null());
        assert_eq!(
            text(&standard_heap, standard_output.szLog).as_deref(),
            Some(concat!(
                "Generating standard InChI\n",
                "Input format: MOLfile\n",
                "Output format: Plain text\n",
                "Full Aux. info\n",
                "No timeout\n",
                "Up to 1024 atoms per structure"
            ))
        );
        assert_eq!(
            standard_heap.slice(standard_input.atom.as_const()).unwrap(),
            standard_atoms
        );
        FreeINCHI(&mut standard_heap, Some(&mut standard_output)).unwrap();
        assert_eq!(standard_heap.live_source_allocation_count(), 1);
        assert_eq!(
            standard_heap.live_source_allocations_of::<crate::source_types::T_GROUP>(),
            1
        );

        for enforce_standard in [0, 1] {
            let mut heap = SourceHeap::default();
            let options = c_string(&mut heap, "-FixedH -RecMet -SUU -SLUUD -SaveOpt -unknown");
            let input = methane_input(&mut heap, options);
            let input_snapshot = input.clone();
            let atom_snapshot = heap.slice(input.atom.as_const()).unwrap().to_vec();
            let mut output = inchi_Output::default();
            assert_eq!(
                invoke(&mut heap, &input, Some(&mut output), enforce_standard, 0,),
                Ok(tagRetValGetINCHI_inchi_Ret_OKAY),
                "enforce_standard={enforce_standard}"
            );
            assert_eq!(
                text(&heap, output.szInChI).as_deref(),
                Some(if enforce_standard == 0 {
                    "InChI=1/CH4/h1H4\\PA"
                } else {
                    "InChI=1S/CH4/h1H4"
                }),
                "enforce_standard={enforce_standard}"
            );
            assert!(
                text(&heap, output.szLog)
                    .unwrap()
                    .contains("Unrecognized optionQ3: \"unknown\"."),
                "enforce_standard={enforce_standard}"
            );
            assert_eq!(input, input_snapshot);
            assert_eq!(heap.slice(input.atom.as_const()).unwrap(), atom_snapshot);
            FreeINCHI(&mut heap, Some(&mut output)).unwrap();
            assert_eq!(heap.live_source_allocation_count(), 1);
            assert_eq!(
                heap.live_source_allocations_of::<crate::source_types::T_GROUP>(),
                1
            );
        }

        let mut interrupted_heap = SourceHeap::default();
        let interrupted_input = methane_input(&mut interrupted_heap, SourceMutPointer::null());
        let mut interrupted_output = inchi_Output::default();
        assert_eq!(
            invoke(
                &mut interrupted_heap,
                &interrupted_input,
                Some(&mut interrupted_output),
                0,
                1,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert!(interrupted_output.szInChI.is_null());
        assert!(interrupted_output.szAuxInfo.is_null());
        assert!(interrupted_output.szMessage.is_null());
        assert_eq!(
            text(&interrupted_heap, interrupted_output.szLog).as_deref(),
            Some(concat!(
                "Generating standard InChI\n",
                "Input format: MOLfile\n",
                "Output format: Plain text\n",
                "Full Aux. info\n",
                "No timeout\n",
                "Up to 1024 atoms per structure"
            ))
        );
        FreeINCHI(&mut interrupted_heap, Some(&mut interrupted_output)).unwrap();
        assert_eq!(interrupted_heap.live_source_allocation_count(), 0);

        let mut option_failure_heap = SourceHeap::default();
        let failed_options = c_string(&mut option_failure_heap, "-FixedH");
        let failed_input = methane_input(&mut option_failure_heap, failed_options);
        let mut failed_output = inchi_Output::default();
        option_failure_heap.fail_after_allocations(0);
        assert_eq!(
            invoke(
                &mut option_failure_heap,
                &failed_input,
                Some(&mut failed_output),
                0,
                0,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_FATAL)
        );
        assert_eq!(option_failure_heap.source_allocation_calls(), 1);
        assert_eq!(failed_output, inchi_Output::default());
        assert_eq!(option_failure_heap.live_source_allocation_count(), 0);

        let mut buffer_failure_heap = SourceHeap::default();
        let buffer_failure_input =
            methane_input(&mut buffer_failure_heap, SourceMutPointer::null());
        let mut buffer_failure_output = inchi_Output::default();
        buffer_failure_heap.fail_after_allocations(1);
        assert_eq!(
            invoke(
                &mut buffer_failure_heap,
                &buffer_failure_input,
                Some(&mut buffer_failure_output),
                0,
                0,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_FATAL)
        );
        assert_eq!(buffer_failure_heap.source_allocation_calls(), 2);
        assert!(buffer_failure_output.szInChI.is_null());
        assert!(buffer_failure_output.szAuxInfo.is_null());
        assert!(buffer_failure_output.szMessage.is_null());
        assert!(
            text(&buffer_failure_heap, buffer_failure_output.szLog)
                .unwrap()
                .ends_with("Cannot allocate internal string buffer. Terminating")
        );
        FreeINCHI(&mut buffer_failure_heap, Some(&mut buffer_failure_output)).unwrap();
        assert_eq!(buffer_failure_heap.live_source_allocation_count(), 0);
    }

    #[test]
    fn source_port__inchi_dll__getinchi__line_270() {
        fn c_string(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .bytes()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect(),
            )
            .unwrap()
        }

        fn text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Option<String> {
            if pointer.is_null() {
                return None;
            }
            let bytes = heap.slice(pointer.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            Some(String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap())
        }

        fn invoke(
            heap: &mut SourceHeap,
            input: Option<&inchi_Input>,
            output: Option<&mut inchi_Output>,
            interrupted: i32,
        ) -> Result<i32, SourceHeapError> {
            GetINCHI(
                heap,
                input,
                output,
                interrupted,
                SourceMutPointer::null(),
                InchiBuildMetadata {
                    compiler: "gcc",
                    date: "Jan 01 1970",
                    time: "00:00:00",
                },
                1_000,
            )
        }

        let mut null_heap = SourceHeap::default();
        let mut null_output = inchi_Output::default();
        assert_eq!(
            invoke(&mut null_heap, None, Some(&mut null_output), 0),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(null_output, inchi_Output::default());

        for element in ["*", "Zz"] {
            let mut heap = SourceHeap::default();
            let atom = heap
                .allocate_model_storage(vec![inchi_Atom {
                    elname: source_name(element),
                    ..inchi_Atom::default()
                }])
                .unwrap();
            let input = inchi_Input {
                atom,
                num_atoms: 1,
                ..inchi_Input::default()
            };
            let old = c_string(&mut heap, "old");
            let mut output = inchi_Output {
                szInChI: old,
                szAuxInfo: old,
                szMessage: old,
                szLog: old,
            };
            assert_eq!(
                invoke(&mut heap, Some(&input), Some(&mut output), 0),
                Ok(_IS_ERROR as i32)
            );
            assert!(output.szInChI.is_null());
            assert!(output.szAuxInfo.is_null());
            assert!(output.szLog.is_null());
            assert_eq!(
                text(&heap, output.szMessage).as_deref(),
                Some("Unsupported in this mode element '*'")
            );
            assert!(heap.slice(old.as_const()).is_ok());
            FreeINCHI(&mut heap, Some(&mut output)).unwrap();
        }

        let mut methane_heap = SourceHeap::default();
        let methane_options = c_string(&mut methane_heap, "-FixedH -RecMet -SUU -SLUUD");
        let methane_atom = methane_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("C"),
                num_iso_H: [-1, 0, 0, 0],
                ..inchi_Atom::default()
            }])
            .unwrap();
        let methane_input = inchi_Input {
            atom: methane_atom,
            stereo0D: SourceMutPointer::null(),
            szOptions: methane_options,
            num_atoms: 1,
            num_stereo0D: 0,
        };
        let methane_snapshot = methane_input.clone();
        let methane_atom_snapshot = methane_heap
            .slice(methane_atom.as_const())
            .unwrap()
            .to_vec();
        let methane_option_snapshot = methane_heap
            .slice(methane_options.as_const())
            .unwrap()
            .to_vec();
        let mut methane_output = inchi_Output::default();
        assert_eq!(
            invoke(
                &mut methane_heap,
                Some(&methane_input),
                Some(&mut methane_output),
                0,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert_eq!(
            text(&methane_heap, methane_output.szInChI).as_deref(),
            Some("InChI=1/CH4/h1H4")
        );
        assert_eq!(
            text(&methane_heap, methane_output.szAuxInfo).as_deref(),
            Some("AuxInfo=1/0/N:1/rA:1C/rB:/rC:;")
        );
        assert!(methane_output.szMessage.is_null());
        assert_eq!(
            text(&methane_heap, methane_output.szLog).as_deref(),
            Some(concat!(
                "Generating non-standard InChI with the options: \n",
                "  Mobile H Perception OFF (include FixedH layer)\n",
                "  Include bonds to metals\n",
                "  Absolute Stereo ON\n",
                "  Include undefined/unknown stereogenic centers and bonds\n",
                "  Make labels for unknown and undefined stereo different\n",
                "  Do not account for keto-enol tautomerism\n",
                "  Do not account for 1,5-tautomerism\n",
                "  Do not account for PT_22_00 tautomerism\n",
                "  Do not account for PT_16_00 tautomerism\n",
                "  Do not account for PT_06_00 tautomerism\n",
                "  Do not account for PT_39_00 tautomerism\n",
                "  Do not account for PT_13_00 tautomerism\n",
                "  Do not account for PT_18_00 tautomerism\n",
                "Input format: MOLfile\n",
                "Output format: Plain text\n",
                "Full Aux. info\n",
                "No timeout\n",
                "Up to 1024 atoms per structure"
            ))
        );
        assert_eq!(methane_input, methane_snapshot);
        assert_eq!(
            methane_heap.slice(methane_atom.as_const()).unwrap(),
            methane_atom_snapshot
        );
        assert_eq!(
            methane_heap.slice(methane_options.as_const()).unwrap(),
            methane_option_snapshot
        );
        FreeINCHI(&mut methane_heap, Some(&mut methane_output)).unwrap();

        let mut stereo_heap = SourceHeap::default();
        let stereo_atoms = stereo_heap
            .allocate_model_storage(vec![
                inchi_Atom {
                    neighbor: [1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    bond_type: [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    elname: source_name("C"),
                    num_bonds: 4,
                    ..inchi_Atom::default()
                },
                inchi_Atom {
                    neighbor: [0; 20],
                    bond_type: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    elname: source_name("F"),
                    num_bonds: 1,
                    ..inchi_Atom::default()
                },
                inchi_Atom {
                    neighbor: [0; 20],
                    bond_type: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    elname: source_name("Cl"),
                    num_bonds: 1,
                    ..inchi_Atom::default()
                },
                inchi_Atom {
                    neighbor: [0; 20],
                    bond_type: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    elname: source_name("Br"),
                    num_bonds: 1,
                    ..inchi_Atom::default()
                },
                inchi_Atom {
                    neighbor: [0; 20],
                    bond_type: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    elname: source_name("I"),
                    num_bonds: 1,
                    ..inchi_Atom::default()
                },
            ])
            .unwrap();
        let stereo = stereo_heap
            .allocate_model_storage(vec![inchi_Stereo0D {
                neighbor: [1, 2, 3, 4],
                central_atom: 0,
                type_: 2,
                parity: 1,
            }])
            .unwrap();
        let stereo_input = inchi_Input {
            atom: stereo_atoms,
            stereo0D: stereo,
            num_atoms: 5,
            num_stereo0D: 1,
            ..inchi_Input::default()
        };
        let mut stereo_output = inchi_Output::default();
        assert_eq!(
            invoke(
                &mut stereo_heap,
                Some(&stereo_input),
                Some(&mut stereo_output),
                0,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert_eq!(
            text(&stereo_heap, stereo_output.szInChI).as_deref(),
            Some("InChI=1S/CBrClFI/c2-1(3,4)5/t1-/m1/s1")
        );
        assert_eq!(
            text(&stereo_heap, stereo_output.szAuxInfo).as_deref(),
            Some("AuxInfo=1/0/N:1,4,3,2,5/it:im/rA:5C.oFClBrI/rB:s1;s1;s1;s1;/rC:;;;;;")
        );
        assert!(stereo_output.szMessage.is_null());
        assert_eq!(
            text(&stereo_heap, stereo_output.szLog).as_deref(),
            Some(concat!(
                "Generating standard InChI\n",
                "Input format: MOLfile\n",
                "Output format: Plain text\n",
                "Full Aux. info\n",
                "No timeout\n",
                "Up to 1024 atoms per structure"
            ))
        );
        FreeINCHI(&mut stereo_heap, Some(&mut stereo_output)).unwrap();

        let mut output_null_heap = SourceHeap::default();
        let output_null_atom = output_null_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("C"),
                ..inchi_Atom::default()
            }])
            .unwrap();
        let output_null_input = inchi_Input {
            atom: output_null_atom,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        assert_eq!(
            invoke(&mut output_null_heap, Some(&output_null_input), None, 0),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        assert_eq!(output_null_heap.live_source_allocation_count(), 0);

        let mut interrupted_heap = SourceHeap::default();
        let interrupted_atom = interrupted_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("C"),
                ..inchi_Atom::default()
            }])
            .unwrap();
        let interrupted_input = inchi_Input {
            atom: interrupted_atom,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        let mut interrupted_output = inchi_Output::default();
        assert_eq!(
            invoke(
                &mut interrupted_heap,
                Some(&interrupted_input),
                Some(&mut interrupted_output),
                1,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert!(interrupted_output.szInChI.is_null());
        assert!(interrupted_output.szAuxInfo.is_null());
        assert!(interrupted_output.szMessage.is_null());
        assert_eq!(
            text(&interrupted_heap, interrupted_output.szLog).as_deref(),
            Some(concat!(
                "Generating standard InChI\n",
                "Input format: MOLfile\n",
                "Output format: Plain text\n",
                "Full Aux. info\n",
                "No timeout\n",
                "Up to 1024 atoms per structure"
            ))
        );
        FreeINCHI(&mut interrupted_heap, Some(&mut interrupted_output)).unwrap();

        let mut failure_heap = SourceHeap::default();
        let failure_options = c_string(&mut failure_heap, "-FixedH");
        let failure_atom = failure_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("C"),
                ..inchi_Atom::default()
            }])
            .unwrap();
        let failure_input = inchi_Input {
            atom: failure_atom,
            szOptions: failure_options,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        let mut failure_output = inchi_Output::default();
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            invoke(
                &mut failure_heap,
                Some(&failure_input),
                Some(&mut failure_output),
                0,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_FATAL)
        );
        assert_eq!(failure_heap.source_allocation_calls(), 1);
        assert_eq!(failure_output, inchi_Output::default());
    }

    #[test]
    fn source_port__inchi_dll__getinchiex__line_325() {
        fn c_string(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .bytes()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect(),
            )
            .unwrap()
        }

        fn text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Option<String> {
            if pointer.is_null() {
                return None;
            }
            let bytes = heap.slice(pointer.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            Some(String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap())
        }

        fn invoke(
            heap: &mut SourceHeap,
            input: Option<&mut inchi_InputEx>,
            output: Option<&mut inchi_Output>,
        ) -> Result<i32, SourceHeapError> {
            GetINCHIEx(
                heap,
                input,
                output,
                0,
                SourceMutPointer::null(),
                InchiBuildMetadata {
                    compiler: "gcc",
                    date: "Jan 01 1970",
                    time: "00:00:00",
                },
                1_000,
            )
        }

        let mut null_heap = SourceHeap::default();
        let mut null_output = inchi_Output::default();
        assert_eq!(
            invoke(&mut null_heap, None, Some(&mut null_output)),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(null_output, inchi_Output::default());

        for count in [i16::MIN, -1, 0] {
            let mut heap = SourceHeap::default();
            let polymer = heap
                .allocate_model_storage(vec![inchi_Input_Polymer::default()])
                .unwrap();
            let v3000 = heap
                .allocate_model_storage(vec![inchi_Input_V3000::default()])
                .unwrap();
            let mut input = inchi_InputEx {
                atom: SourceMutPointer::null(),
                num_atoms: count,
                polymer,
                v3000,
                ..inchi_InputEx::default()
            };
            let snapshot = input.clone();
            assert_eq!(
                invoke(&mut heap, Some(&mut input), None),
                Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
            );
            assert_eq!(input, snapshot);
        }

        let mut mutation_heap = SourceHeap::default();
        let mut first_star = inchi_Atom::default();
        first_star.elname = [b'*' as i8, 0, 11, 12, 13, 14];
        let carbon = inchi_Atom {
            elname: source_name("C"),
            ..inchi_Atom::default()
        };
        let mut second_star = inchi_Atom::default();
        second_star.elname = [b'*' as i8, 0, 21, 22, 23, 24];
        let atoms = mutation_heap
            .allocate_model_storage(vec![first_star, carbon.clone(), second_star])
            .unwrap();
        let polymer = mutation_heap
            .allocate_model_storage(vec![inchi_Input_Polymer::default()])
            .unwrap();
        let v3000 = mutation_heap
            .allocate_model_storage(vec![inchi_Input_V3000::default()])
            .unwrap();
        let mut mutation_input = inchi_InputEx {
            atom: atoms,
            num_atoms: 3,
            polymer,
            v3000,
            ..inchi_InputEx::default()
        };
        assert_eq!(
            invoke(&mut mutation_heap, Some(&mut mutation_input), None),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        assert_eq!(mutation_input.atom, atoms);
        assert_eq!(mutation_input.polymer, polymer);
        assert_eq!(mutation_input.v3000, v3000);
        let mutated = mutation_heap.slice(atoms.as_const()).unwrap();
        assert_eq!(mutated[0].elname, [b'Z' as i8, b'z' as i8, 0, 12, 13, 14]);
        assert_eq!(mutated[1], carbon);
        assert_eq!(mutated[2].elname, [b'Z' as i8, b'z' as i8, 0, 22, 23, 24]);

        let mut methane_heap = SourceHeap::default();
        let options = c_string(&mut methane_heap, "-FixedH -RecMet -SUU -SLUUD");
        let methane_atom = methane_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("C"),
                num_iso_H: [-1, 0, 0, 0],
                ..inchi_Atom::default()
            }])
            .unwrap();
        let mut methane_input = inchi_InputEx {
            atom: methane_atom,
            szOptions: options,
            num_atoms: 1,
            ..inchi_InputEx::default()
        };
        let methane_snapshot = methane_input.clone();
        let mut methane_output = inchi_Output::default();
        assert_eq!(
            invoke(
                &mut methane_heap,
                Some(&mut methane_input),
                Some(&mut methane_output),
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert_eq!(methane_input, methane_snapshot);
        assert_eq!(
            text(&methane_heap, methane_output.szInChI).as_deref(),
            Some("InChI=1/CH4/h1H4")
        );
        assert_eq!(
            text(&methane_heap, methane_output.szAuxInfo).as_deref(),
            Some("AuxInfo=1/0/N:1/rA:1C/rB:/rC:;")
        );
        assert!(methane_output.szMessage.is_null());
        assert_eq!(
            text(&methane_heap, methane_output.szLog).as_deref(),
            Some(concat!(
                "Generating non-standard InChI with the options: \n",
                "  Mobile H Perception OFF (include FixedH layer)\n",
                "  Include bonds to metals\n",
                "  Absolute Stereo ON\n",
                "  Include undefined/unknown stereogenic centers and bonds\n",
                "  Make labels for unknown and undefined stereo different\n",
                "  Do not account for keto-enol tautomerism\n",
                "  Do not account for 1,5-tautomerism\n",
                "  Do not account for PT_22_00 tautomerism\n",
                "  Do not account for PT_16_00 tautomerism\n",
                "  Do not account for PT_06_00 tautomerism\n",
                "  Do not account for PT_39_00 tautomerism\n",
                "  Do not account for PT_13_00 tautomerism\n",
                "  Do not account for PT_18_00 tautomerism\n",
                "Input format: MOLfile\n",
                "Output format: Plain text\n",
                "Full Aux. info\n",
                "No timeout\n",
                "Up to 1024 atoms per structure"
            ))
        );
        FreeINCHI(&mut methane_heap, Some(&mut methane_output)).unwrap();

        let mut failure_heap = SourceHeap::default();
        let failure_options = c_string(&mut failure_heap, "-FixedH");
        let failure_atom = failure_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("C"),
                ..inchi_Atom::default()
            }])
            .unwrap();
        let mut failure_input = inchi_InputEx {
            atom: failure_atom,
            szOptions: failure_options,
            num_atoms: 1,
            ..inchi_InputEx::default()
        };
        let mut failure_output = inchi_Output::default();
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            invoke(
                &mut failure_heap,
                Some(&mut failure_input),
                Some(&mut failure_output),
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_FATAL)
        );
        assert_eq!(failure_heap.source_allocation_calls(), 1);
        assert_eq!(failure_output, inchi_Output::default());
    }

    #[test]
    fn source_port__inchi_dll__getstdinchi__line_242() {
        fn c_string(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .bytes()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect(),
            )
            .unwrap()
        }

        fn text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Option<String> {
            if pointer.is_null() {
                return None;
            }
            let bytes = heap.slice(pointer.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            Some(String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap())
        }

        fn invoke(
            heap: &mut SourceHeap,
            input: Option<&inchi_Input>,
            output: Option<&mut inchi_Output>,
            interrupted: i32,
        ) -> Result<i32, SourceHeapError> {
            GetStdINCHI(
                heap,
                input,
                output,
                interrupted,
                SourceMutPointer::null(),
                InchiBuildMetadata {
                    compiler: "gcc",
                    date: "Jan 01 1970",
                    time: "00:00:00",
                },
                1_000,
            )
        }

        let mut null_heap = SourceHeap::default();
        assert_eq!(
            invoke(&mut null_heap, None, None, 0),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(null_heap.live_source_allocation_count(), 0);

        let mut pseudo_heap = SourceHeap::default();
        let pseudo_atoms = pseudo_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("*"),
                ..inchi_Atom::default()
            }])
            .unwrap();
        let pseudo_input = inchi_Input {
            atom: pseudo_atoms,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        let old_inchi = c_string(&mut pseudo_heap, "old inchi");
        let old_log = c_string(&mut pseudo_heap, "old log");
        let mut pseudo_output = inchi_Output {
            szInChI: old_inchi,
            szLog: old_log,
            ..inchi_Output::default()
        };
        assert_eq!(
            invoke(
                &mut pseudo_heap,
                Some(&pseudo_input),
                Some(&mut pseudo_output),
                0,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert!(pseudo_output.szInChI.is_null());
        assert!(pseudo_output.szAuxInfo.is_null());
        assert!(pseudo_output.szLog.is_null());
        assert_eq!(
            text(&pseudo_heap, pseudo_output.szMessage).as_deref(),
            Some("Unsupported in this mode element '*'")
        );
        assert!(pseudo_heap.slice(old_inchi.as_const()).is_ok());
        assert!(pseudo_heap.slice(old_log.as_const()).is_ok());
        FreeStdINCHI(&mut pseudo_heap, Some(&mut pseudo_output)).unwrap();
        assert_eq!(pseudo_heap.live_source_allocation_count(), 0);

        let mut pseudo_null_output_heap = SourceHeap::default();
        let pseudo_null_atoms = pseudo_null_output_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("Zz"),
                ..inchi_Atom::default()
            }])
            .unwrap();
        let pseudo_null_input = inchi_Input {
            atom: pseudo_null_atoms,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        assert_eq!(
            invoke(
                &mut pseudo_null_output_heap,
                Some(&pseudo_null_input),
                None,
                0,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(pseudo_null_output_heap.source_allocation_calls(), 0);

        let mut pseudo_failure_heap = SourceHeap::default();
        let pseudo_failure_atoms = pseudo_failure_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("*"),
                ..inchi_Atom::default()
            }])
            .unwrap();
        let pseudo_failure_input = inchi_Input {
            atom: pseudo_failure_atoms,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        let mut pseudo_failure_output = inchi_Output::default();
        pseudo_failure_heap.fail_after_allocations(0);
        assert_eq!(
            invoke(
                &mut pseudo_failure_heap,
                Some(&pseudo_failure_input),
                Some(&mut pseudo_failure_output),
                0,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(pseudo_failure_heap.source_allocation_calls(), 1);
        assert_eq!(pseudo_failure_output, inchi_Output::default());

        let mut output_null_heap = SourceHeap::default();
        let output_null_atoms = output_null_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("C"),
                num_iso_H: [-1, 0, 0, 0],
                ..inchi_Atom::default()
            }])
            .unwrap();
        let output_null_input = inchi_Input {
            atom: output_null_atoms,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        assert_eq!(
            invoke(&mut output_null_heap, Some(&output_null_input), None, 0),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        assert_eq!(output_null_heap.live_source_allocation_count(), 0);

        let mut standard_heap = SourceHeap::default();
        let options = c_string(&mut standard_heap, "-FixedH -RecMet -SUU -SLUUD -SaveOpt");
        let atoms = standard_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("C"),
                num_iso_H: [-1, 0, 0, 0],
                ..inchi_Atom::default()
            }])
            .unwrap();
        let input = inchi_Input {
            atom: atoms,
            szOptions: options,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        let input_snapshot = input.clone();
        let atom_snapshot = standard_heap.slice(atoms.as_const()).unwrap().to_vec();
        let option_snapshot = standard_heap.slice(options.as_const()).unwrap().to_vec();
        let mut output = inchi_Output::default();
        assert_eq!(
            invoke(&mut standard_heap, Some(&input), Some(&mut output), 0),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert_eq!(
            text(&standard_heap, output.szInChI).as_deref(),
            Some("InChI=1S/CH4/h1H4")
        );
        assert_eq!(
            text(&standard_heap, output.szAuxInfo).as_deref(),
            Some("AuxInfo=1/0/N:1/rA:1C/rB:/rC:;")
        );
        assert!(output.szMessage.is_null());
        assert_eq!(input, input_snapshot);
        assert_eq!(
            standard_heap.slice(atoms.as_const()).unwrap(),
            atom_snapshot
        );
        assert_eq!(
            standard_heap.slice(options.as_const()).unwrap(),
            option_snapshot
        );
        FreeStdINCHI(&mut standard_heap, Some(&mut output)).unwrap();
        assert_eq!(standard_heap.live_source_allocation_count(), 1);
        assert_eq!(
            standard_heap.live_source_allocations_of::<crate::source_types::T_GROUP>(),
            1
        );

        let mut interrupted_heap = SourceHeap::default();
        let interrupted_atoms = interrupted_heap
            .allocate_model_storage(vec![inchi_Atom {
                elname: source_name("C"),
                ..inchi_Atom::default()
            }])
            .unwrap();
        let interrupted_input = inchi_Input {
            atom: interrupted_atoms,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        let mut interrupted_output = inchi_Output::default();
        assert_eq!(
            invoke(
                &mut interrupted_heap,
                Some(&interrupted_input),
                Some(&mut interrupted_output),
                1,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert!(interrupted_output.szInChI.is_null());
        assert!(interrupted_output.szAuxInfo.is_null());
        assert!(interrupted_output.szMessage.is_null());
        assert!(text(&interrupted_heap, interrupted_output.szLog).is_some());
        FreeStdINCHI(&mut interrupted_heap, Some(&mut interrupted_output)).unwrap();
        assert_eq!(interrupted_heap.live_source_allocation_count(), 0);
    }

    #[test]
    fn source_port__scc__getinchifrominchi__line_2114() {
        fn c_string(heap: &mut SourceHeap, text: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                text.bytes()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect(),
            )
            .unwrap()
        }

        fn text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Option<String> {
            if pointer.is_null() {
                return None;
            }
            let bytes = heap.slice(pointer.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            Some(String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap())
        }

        let build = InchiBuildMetadata {
            compiler: "gcc",
            date: "Jan 01 1970",
            time: "00:00:00",
        };
        let stdout = SourceMutPointer::null();

        let mut heap = SourceHeap::default();
        assert_eq!(
            CheckINCHI(
                &mut heap,
                SourceConstPointer::null(),
                0,
                stdout,
                build,
                1_000,
            ),
            Ok(tagRetValCheckINCHI_INCHI_INVALID_PREFIX)
        );
        for (input, expected) in [
            ("", tagRetValCheckINCHI_INCHI_INVALID_PREFIX),
            ("InChI=1", tagRetValCheckINCHI_INCHI_INVALID_PREFIX),
            ("inchi=1S/CH4", tagRetValCheckINCHI_INCHI_INVALID_PREFIX),
            ("InChI=2S/CH4", tagRetValCheckINCHI_INCHI_INVALID_VERSION),
            ("InChI=1SCH4", tagRetValCheckINCHI_INCHI_INVALID_LAYOUT),
            ("InChI=1S/0", tagRetValCheckINCHI_INCHI_INVALID_LAYOUT),
            ("InChI=1S/CH4/C1", tagRetValCheckINCHI_INCHI_INVALID_LAYOUT),
            (
                "InChI=1S/CH4/c1=2",
                tagRetValCheckINCHI_INCHI_INVALID_LAYOUT,
            ),
            (
                "InChI=1S/CH4/c1@2",
                tagRetValCheckINCHI_INCHI_INVALID_LAYOUT,
            ),
            ("InChI=1S/CH4", tagRetValCheckINCHI_INCHI_VALID_STANDARD),
            ("InChI=1/CH4", tagRetValCheckINCHI_INCHI_VALID_NON_STANDARD),
            ("InChI=1B/CH4", tagRetValCheckINCHI_INCHI_VALID_BETA),
            (
                "InChI=1S/CH4/c1?2;3",
                tagRetValCheckINCHI_INCHI_VALID_STANDARD,
            ),
            ("InChI=1S/CH4\\PB", tagRetValCheckINCHI_INCHI_VALID_STANDARD),
            (
                "InChI=1S/CH4  \t\r\n",
                tagRetValCheckINCHI_INCHI_VALID_STANDARD,
            ),
        ] {
            let pointer = c_string(&mut heap, input);
            assert_eq!(
                CheckINCHI(&mut heap, pointer.as_const(), 0, stdout, build, 1_000),
                Ok(expected),
                "{input:?}"
            );
            heap.free(pointer).unwrap();
        }

        let old_inchi = c_string(&mut heap, "old inchi");
        let old_aux = c_string(&mut heap, "old aux");
        let old_message = c_string(&mut heap, "old message");
        let old_log = c_string(&mut heap, "old log");
        let mut null_output = inchi_Output {
            szInChI: old_inchi,
            szAuxInfo: old_aux,
            szMessage: old_message,
            szLog: old_log,
        };
        assert_eq!(
            GetINCHIfromINCHI(&mut heap, None, &mut null_output, stdout, build, 1_000,),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        assert_eq!(null_output, inchi_Output::default());
        for pointer in [old_inchi, old_aux, old_message, old_log] {
            assert!(heap.slice(pointer.as_const()).is_ok());
        }

        const BASE_LOG: &str = concat!(
            "Generating standard InChI\n",
            "Input format: InChI (plain identifier)\n",
            "Output format: Plain text\n",
            "Aux. info suppressed\n",
            "No timeout\n",
            "Up to 1024 atoms per structure",
        );

        let methane = c_string(&mut heap, "InChI=1S/CH4/h1H4");
        let input = inchi_InputINCHI {
            szInChI: methane,
            szOptions: SourceMutPointer::null(),
        };
        let mut output = inchi_Output::default();
        let result = GetINCHIfromINCHI(&mut heap, Some(&input), &mut output, stdout, build, 1_000);
        assert_eq!(result, Ok(0));
        assert_eq!(
            text(&heap, output.szInChI).as_deref(),
            Some("InChI=1S/CH4/h1H4")
        );
        assert!(output.szAuxInfo.is_null());
        assert_eq!(text(&heap, output.szMessage).as_deref(), Some(""));
        assert_eq!(text(&heap, output.szLog).as_deref(), Some(BASE_LOG));

        let nonstandard = c_string(&mut heap, "InChI=1/CH4/h1H4");
        let nonstandard_input = inchi_InputINCHI {
            szInChI: nonstandard,
            szOptions: SourceMutPointer::null(),
        };
        let mut nonstandard_output = inchi_Output::default();
        let nonstandard_result = GetINCHIfromINCHI(
            &mut heap,
            Some(&nonstandard_input),
            &mut nonstandard_output,
            stdout,
            build,
            1_000,
        );
        assert_eq!(nonstandard_result, Ok(0));
        assert_eq!(
            text(&heap, nonstandard_output.szInChI).as_deref(),
            Some("InChI=1/CH4/h1H4")
        );
        assert!(nonstandard_output.szAuxInfo.is_null());
        assert_eq!(
            text(&heap, nonstandard_output.szMessage).as_deref(),
            Some("")
        );
        assert_eq!(
            text(&heap, nonstandard_output.szLog).as_deref(),
            Some(BASE_LOG)
        );

        let strict_options = c_string(&mut heap, "-FixedH -RecMet -SUU -SLUUD");
        let strict_input = inchi_InputINCHI {
            szInChI: methane,
            szOptions: strict_options,
        };
        let mut strict_output = inchi_Output::default();
        let strict_get = GetINCHIfromINCHI(
            &mut heap,
            Some(&strict_input),
            &mut strict_output,
            stdout,
            build,
            1_000,
        );
        assert_eq!(strict_get, Ok(0));
        assert_eq!(
            text(&heap, strict_output.szInChI).as_deref(),
            Some("InChI=1/CH4/h1H4")
        );
        assert!(strict_output.szAuxInfo.is_null());
        assert_eq!(text(&heap, strict_output.szMessage).as_deref(), Some(""));
        assert_eq!(
            text(&heap, strict_output.szLog).as_deref(),
            Some(concat!(
                "Generating non-standard InChI with the options: \n",
                "  Mobile H Perception OFF (include FixedH layer)\n",
                "  Include bonds to metals\n",
                "  Absolute Stereo ON\n",
                "  Include undefined/unknown stereogenic centers and bonds\n",
                "  Make labels for unknown and undefined stereo different\n",
                "  Do not account for keto-enol tautomerism\n",
                "  Do not account for 1,5-tautomerism\n",
                "  Do not account for PT_22_00 tautomerism\n",
                "  Do not account for PT_16_00 tautomerism\n",
                "  Do not account for PT_06_00 tautomerism\n",
                "  Do not account for PT_39_00 tautomerism\n",
                "  Do not account for PT_13_00 tautomerism\n",
                "  Do not account for PT_18_00 tautomerism\n",
                "Input format: InChI (plain identifier)\n",
                "Output format: Plain text\n",
                "Aux. info suppressed\n",
                "No timeout\n",
                "Up to 1024 atoms per structure",
            ))
        );

        let strict = c_string(&mut heap, "InChI=1S/CH4/h1H4");
        let strict_result = CheckINCHI(&mut heap, strict.as_const(), 1, stdout, build, 1_000);
        assert_eq!(strict_result, Ok(tagRetValCheckINCHI_INCHI_FAIL_I2I));

        let strict_nonstandard = c_string(&mut heap, "InChI=1/CH4/h1H4");
        assert_eq!(
            CheckINCHI(
                &mut heap,
                strict_nonstandard.as_const(),
                1,
                stdout,
                build,
                1_000,
            ),
            Ok(tagRetValCheckINCHI_INCHI_VALID_NON_STANDARD)
        );

        let stereo = c_string(
            &mut heap,
            "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1",
        );
        let filter_options = c_string(&mut heap, "-SNon");
        let filter_input = inchi_InputINCHI {
            szInChI: stereo,
            szOptions: filter_options,
        };
        let mut filter_output = inchi_Output::default();
        let filter_result = GetINCHIfromINCHI(
            &mut heap,
            Some(&filter_input),
            &mut filter_output,
            stdout,
            build,
            1_000,
        );
        assert_eq!(filter_result, Ok(0));
        assert_eq!(
            text(&heap, filter_output.szInChI).as_deref(),
            Some("InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)")
        );
        assert_eq!(text(&heap, filter_output.szMessage).as_deref(), Some(""));
        assert_eq!(
            text(&heap, filter_output.szLog).as_deref(),
            Some(concat!(
                "Using specific structure perception features:\n",
                "  Stereo OFF\n",
                "Generating standard InChI\n",
                "Input format: InChI (plain identifier)\n",
                "Output format: Plain text\n",
                "Aux. info suppressed\n",
                "No timeout\n",
                "Up to 1024 atoms per structure",
            ))
        );

        let invalid_get_text = c_string(&mut heap, "InChI=2S/CH4");
        let invalid_get_input = inchi_InputINCHI {
            szInChI: invalid_get_text,
            szOptions: SourceMutPointer::null(),
        };
        let mut invalid_get_output = inchi_Output::default();
        let invalid_get_result = GetINCHIfromINCHI(
            &mut heap,
            Some(&invalid_get_input),
            &mut invalid_get_output,
            stdout,
            build,
            1_000,
        );
        assert_eq!(invalid_get_result, Ok(tagRetValGetINCHI_inchi_Ret_ERROR));
        assert_eq!(invalid_get_output, inchi_Output::default());

        let mut trace_heap = SourceHeap::default();
        let trace_text = c_string(&mut trace_heap, "InChI=1S/CH4/h1H4");
        let trace_input = inchi_InputINCHI {
            szInChI: trace_text,
            szOptions: SourceMutPointer::null(),
        };
        let mut trace_output = inchi_Output::default();
        trace_heap.trace_source_allocations();
        let trace_result = GetINCHIfromINCHI(
            &mut trace_heap,
            Some(&trace_input),
            &mut trace_output,
            stdout,
            build,
            1_000,
        );
        assert_eq!(trace_result, Ok(0));
        assert_eq!(trace_heap.source_allocation_calls(), 20);
        assert_eq!(trace_heap.live_allocation_count(), 92);
        assert_eq!(
            text(&trace_heap, trace_output.szInChI).as_deref(),
            Some("InChI=1S/CH4/h1H4")
        );
        assert_eq!(
            text(&trace_heap, trace_output.szMessage).as_deref(),
            Some("")
        );
        assert_eq!(
            text(&trace_heap, trace_output.szLog).as_deref(),
            Some(BASE_LOG)
        );

        let expected_results = [3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 0, 0, 0, 0];
        let expected_calls = [
            1, 21, 3, 4, 5, 6, 7, 8, 9, 10, 20, 12, 14, 14, 15, 16, 20, 19, 19, 21,
        ];
        let expected_live = [
            1, 92, 2, 3, 3, 3, 3, 3, 3, 3, 91, 3, 3, 3, 3, 4, 92, 3, 3, 92,
        ];
        let expected_inchi = [
            None,
            Some("InChI=1S/CH4/h1H4"),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some("InChI=1S/h1H4"),
            None,
            None,
            None,
            None,
            None,
            Some("InChI=1S/CH4/h1H4"),
            None,
            None,
            Some("1S/CH4/h1H4"),
        ];
        for failure_after in 0..trace_heap.source_allocation_calls() {
            let mut failed_heap = SourceHeap::default();
            let failed_text = c_string(&mut failed_heap, "InChI=1S/CH4/h1H4");
            let failed_input = inchi_InputINCHI {
                szInChI: failed_text,
                szOptions: SourceMutPointer::null(),
            };
            let mut failed_output = inchi_Output::default();
            failed_heap.fail_after_allocations(failure_after);
            let failed_result = GetINCHIfromINCHI(
                &mut failed_heap,
                Some(&failed_input),
                &mut failed_output,
                stdout,
                build,
                1_000,
            );
            let index = failure_after as usize;
            assert_eq!(
                failed_result,
                Ok(expected_results[index]),
                "{failure_after}"
            );
            assert_eq!(
                failed_heap.source_allocation_calls(),
                expected_calls[index],
                "{failure_after}"
            );
            assert_eq!(
                failed_heap.live_allocation_count(),
                expected_live[index],
                "{failure_after}"
            );
            assert_eq!(
                text(&failed_heap, failed_output.szInChI).as_deref(),
                expected_inchi[index],
                "{failure_after}"
            );
            assert!(failed_output.szAuxInfo.is_null(), "{failure_after}");
            assert_eq!(
                text(&failed_heap, failed_output.szMessage).as_deref(),
                if matches!(index, 0 | 2) {
                    None
                } else {
                    Some("")
                },
                "{failure_after}"
            );
            let expected_log = match index {
                0 => None,
                1 => Some(
                    concat!(
                        "Input format: InChI (plain identifier)\n",
                        "Output format: Plain text\n",
                        "Aux. info suppressed\n",
                        "No timeout\n",
                        "Up to 1024 atoms per structure",
                    )
                    .to_owned(),
                ),
                2 => Some(format!(
                    "{BASE_LOG}\n\n\nCannot allocate output message buffer."
                )),
                4..=9 => {
                    let suffix = match index {
                        4..=6 => "MOBILE_H_FORMULA (0)",
                        7 => "MOBILE_H_CONNECTIONS (1)",
                        8 => "MOBILE_H (2)",
                        9 => "MATERIAL_BALANCE (35)",
                        _ => unreachable!(),
                    };
                    Some(format!(
                        "{BASE_LOG}\n\n\n\nStructure: 1 Allocation failed (-1) in {suffix}"
                    ))
                }
                _ => Some(BASE_LOG.to_owned()),
            };
            assert_eq!(
                text(&failed_heap, failed_output.szLog),
                expected_log,
                "{failure_after}"
            );
        }

        let mut failure_heap = SourceHeap::default();
        let failure_input = c_string(&mut failure_heap, "InChI=1S/CH4/h1H4");
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            CheckINCHI(
                &mut failure_heap,
                failure_input.as_const(),
                1,
                stdout,
                build,
                1_000,
            ),
            Ok(tagRetValCheckINCHI_INCHI_FAIL_I2I)
        );

        let mut get_failure_heap = SourceHeap::default();
        let get_failure_text = c_string(&mut get_failure_heap, "InChI=1S/CH4/h1H4");
        let get_failure_input = inchi_InputINCHI {
            szInChI: get_failure_text,
            szOptions: SourceMutPointer::null(),
        };
        let mut get_failure_output = inchi_Output::default();
        get_failure_heap.fail_after_allocations(0);
        assert_eq!(
            GetINCHIfromINCHI(
                &mut get_failure_heap,
                Some(&get_failure_input),
                &mut get_failure_output,
                stdout,
                build,
                1_000,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_FATAL)
        );
        assert_eq!(get_failure_output, inchi_Output::default());
    }

    #[test]
    fn source_port__inchi_dll__freeinchiextinput__line_2923() {
        use crate::source_types::inchi_Input_PolymerUnit;

        fn assert_freed<T: 'static>(heap: &SourceHeap, pointer: SourceMutPointer<T>) {
            assert_eq!(
                heap.slice(pointer.as_const()).map(|_| ()),
                Err(SourceHeapError::MissingAllocation)
            );
        }

        fn assert_live<T: 'static>(heap: &SourceHeap, pointer: SourceMutPointer<T>) {
            assert!(heap.slice(pointer.as_const()).is_ok());
        }

        fn row(heap: &mut SourceHeap, value: i32) -> SourceMutPointer<i32> {
            heap.allocate(vec![value]).unwrap()
        }

        fn outer(
            heap: &mut SourceHeap,
            row: SourceMutPointer<i32>,
        ) -> SourceMutPointer<SourceMutPointer<i32>> {
            heap.allocate(vec![row, SourceMutPointer::null()]).unwrap()
        }

        let mut null_heap = SourceHeap::default();
        assert_eq!(
            FreeInChIExtInput(
                &mut null_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(())
        );
        assert_eq!(null_heap.live_source_allocation_count(), 0);

        let mut polymer_heap = SourceHeap::default();
        let alist = polymer_heap.allocate(vec![1_i32, 2]).unwrap();
        let blist = polymer_heap.allocate(vec![3_i32, 4, 5, 6]).unwrap();
        let unit = polymer_heap
            .allocate(vec![inchi_Input_PolymerUnit {
                alist,
                blist,
                ..inchi_Input_PolymerUnit::default()
            }])
            .unwrap();
        let units = polymer_heap
            .allocate(vec![unit, SourceMutPointer::null()])
            .unwrap();
        let polymer = polymer_heap
            .allocate(vec![inchi_Input_Polymer { units, n: 2 }])
            .unwrap();
        assert_eq!(polymer_heap.live_source_allocation_count(), 5);
        FreeInChIExtInput(&mut polymer_heap, polymer, SourceMutPointer::null()).unwrap();
        for pointer in [alist, blist] {
            assert_freed(&polymer_heap, pointer);
        }
        assert_freed(&polymer_heap, unit);
        assert_freed(&polymer_heap, units);
        assert_freed(&polymer_heap, polymer);
        assert_eq!(polymer_heap.live_source_allocation_count(), 0);

        let mut zero_polymer_heap = SourceHeap::default();
        let zero_alist = row(&mut zero_polymer_heap, 7);
        let zero_unit = zero_polymer_heap
            .allocate(vec![inchi_Input_PolymerUnit {
                alist: zero_alist,
                ..inchi_Input_PolymerUnit::default()
            }])
            .unwrap();
        let zero_units = zero_polymer_heap.allocate(vec![zero_unit]).unwrap();
        let zero_polymer = zero_polymer_heap
            .allocate(vec![inchi_Input_Polymer {
                units: zero_units,
                n: 0,
            }])
            .unwrap();
        FreeInChIExtInput(
            &mut zero_polymer_heap,
            zero_polymer,
            SourceMutPointer::null(),
        )
        .unwrap();
        assert_live(&zero_polymer_heap, zero_alist);
        assert_live(&zero_polymer_heap, zero_unit);
        assert_live(&zero_polymer_heap, zero_units);
        assert_live(&zero_polymer_heap, zero_polymer);
        for pointer in [zero_alist] {
            inchi_free(&mut zero_polymer_heap, pointer).unwrap();
        }
        inchi_free(&mut zero_polymer_heap, zero_unit).unwrap();
        inchi_free(&mut zero_polymer_heap, zero_units).unwrap();
        inchi_free(&mut zero_polymer_heap, zero_polymer).unwrap();

        let mut null_units_heap = SourceHeap::default();
        let null_units_polymer = null_units_heap
            .allocate(vec![inchi_Input_Polymer {
                units: SourceMutPointer::null(),
                n: 1,
            }])
            .unwrap();
        FreeInChIExtInput(
            &mut null_units_heap,
            null_units_polymer,
            SourceMutPointer::null(),
        )
        .unwrap();
        assert_live(&null_units_heap, null_units_polymer);
        inchi_free(&mut null_units_heap, null_units_polymer).unwrap();

        let mut negative_polymer_heap = SourceHeap::default();
        let negative_alist = row(&mut negative_polymer_heap, 8);
        let negative_unit = negative_polymer_heap
            .allocate(vec![inchi_Input_PolymerUnit {
                alist: negative_alist,
                ..inchi_Input_PolymerUnit::default()
            }])
            .unwrap();
        let negative_units = negative_polymer_heap.allocate(vec![negative_unit]).unwrap();
        let negative_polymer = negative_polymer_heap
            .allocate(vec![inchi_Input_Polymer {
                units: negative_units,
                n: -1,
            }])
            .unwrap();
        FreeInChIExtInput(
            &mut negative_polymer_heap,
            negative_polymer,
            SourceMutPointer::null(),
        )
        .unwrap();
        assert_live(&negative_polymer_heap, negative_alist);
        assert_live(&negative_polymer_heap, negative_unit);
        assert_freed(&negative_polymer_heap, negative_units);
        assert_freed(&negative_polymer_heap, negative_polymer);
        inchi_free(&mut negative_polymer_heap, negative_alist).unwrap();
        inchi_free(&mut negative_polymer_heap, negative_unit).unwrap();

        let mut full_v3000_heap = SourceHeap::default();
        let atom_index_orig = full_v3000_heap.allocate(vec![0_i32, 1]).unwrap();
        let atom_index_fin = full_v3000_heap.allocate(vec![1_i32, 0]).unwrap();
        let haptic_row = row(&mut full_v3000_heap, 10);
        let steabs_row = row(&mut full_v3000_heap, 11);
        let sterel_row = row(&mut full_v3000_heap, 12);
        let sterac_row = row(&mut full_v3000_heap, 13);
        let haptic = outer(&mut full_v3000_heap, haptic_row);
        let steabs = outer(&mut full_v3000_heap, steabs_row);
        let sterel = outer(&mut full_v3000_heap, sterel_row);
        let sterac = outer(&mut full_v3000_heap, sterac_row);
        let full_v3000 = full_v3000_heap
            .allocate(vec![inchi_Input_V3000 {
                atom_index_orig,
                atom_index_fin,
                n_haptic_bonds: 2,
                lists_haptic_bonds: haptic,
                n_steabs: 2,
                lists_steabs: steabs,
                n_sterel: 2,
                lists_sterel: sterel,
                n_sterac: 2,
                lists_sterac: sterac,
                ..inchi_Input_V3000::default()
            }])
            .unwrap();
        FreeInChIExtInput(&mut full_v3000_heap, SourceMutPointer::null(), full_v3000).unwrap();
        for pointer in [atom_index_orig, atom_index_fin] {
            assert_freed(&full_v3000_heap, pointer);
        }
        for pointer in [haptic_row, steabs_row, sterel_row, sterac_row] {
            assert_freed(&full_v3000_heap, pointer);
        }
        for pointer in [haptic, steabs, sterel, sterac] {
            assert_freed(&full_v3000_heap, pointer);
        }
        assert_freed(&full_v3000_heap, full_v3000);
        assert_eq!(full_v3000_heap.live_source_allocation_count(), 0);

        let mut zero_v3000_heap = SourceHeap::default();
        let zero_rows = [
            row(&mut zero_v3000_heap, 20),
            row(&mut zero_v3000_heap, 21),
            row(&mut zero_v3000_heap, 22),
            row(&mut zero_v3000_heap, 23),
        ];
        let zero_lists = [
            outer(&mut zero_v3000_heap, zero_rows[0]),
            outer(&mut zero_v3000_heap, zero_rows[1]),
            outer(&mut zero_v3000_heap, zero_rows[2]),
            outer(&mut zero_v3000_heap, zero_rows[3]),
        ];
        let zero_v3000 = zero_v3000_heap
            .allocate(vec![inchi_Input_V3000 {
                lists_haptic_bonds: zero_lists[0],
                lists_steabs: zero_lists[1],
                lists_sterel: zero_lists[2],
                lists_sterac: zero_lists[3],
                ..inchi_Input_V3000::default()
            }])
            .unwrap();
        FreeInChIExtInput(&mut zero_v3000_heap, SourceMutPointer::null(), zero_v3000).unwrap();
        assert_freed(&zero_v3000_heap, zero_v3000);
        for pointer in zero_rows {
            assert_live(&zero_v3000_heap, pointer);
        }
        for pointer in zero_lists {
            assert_live(&zero_v3000_heap, pointer);
        }
        for pointer in zero_rows {
            inchi_free(&mut zero_v3000_heap, pointer).unwrap();
        }
        for pointer in zero_lists {
            inchi_free(&mut zero_v3000_heap, pointer).unwrap();
        }

        let mut negative_v3000_heap = SourceHeap::default();
        let negative_rows = [
            row(&mut negative_v3000_heap, 30),
            row(&mut negative_v3000_heap, 31),
            row(&mut negative_v3000_heap, 32),
            row(&mut negative_v3000_heap, 33),
        ];
        let negative_lists = [
            outer(&mut negative_v3000_heap, negative_rows[0]),
            outer(&mut negative_v3000_heap, negative_rows[1]),
            outer(&mut negative_v3000_heap, negative_rows[2]),
            outer(&mut negative_v3000_heap, negative_rows[3]),
        ];
        let negative_v3000 = negative_v3000_heap
            .allocate(vec![inchi_Input_V3000 {
                n_haptic_bonds: -1,
                lists_haptic_bonds: negative_lists[0],
                n_steabs: -1,
                lists_steabs: negative_lists[1],
                n_sterel: -1,
                lists_sterel: negative_lists[2],
                n_sterac: -1,
                lists_sterac: negative_lists[3],
                ..inchi_Input_V3000::default()
            }])
            .unwrap();
        FreeInChIExtInput(
            &mut negative_v3000_heap,
            SourceMutPointer::null(),
            negative_v3000,
        )
        .unwrap();
        assert_freed(&negative_v3000_heap, negative_v3000);
        for pointer in negative_rows {
            assert_live(&negative_v3000_heap, pointer);
            inchi_free(&mut negative_v3000_heap, pointer).unwrap();
        }
        for pointer in negative_lists {
            assert_freed(&negative_v3000_heap, pointer);
        }

        let mut null_lists_heap = SourceHeap::default();
        let null_lists_v3000 = null_lists_heap
            .allocate(vec![inchi_Input_V3000 {
                n_haptic_bonds: i32::MAX,
                n_steabs: 1,
                n_sterel: -1,
                n_sterac: i32::MIN,
                ..inchi_Input_V3000::default()
            }])
            .unwrap();
        FreeInChIExtInput(
            &mut null_lists_heap,
            SourceMutPointer::null(),
            null_lists_v3000,
        )
        .unwrap();
        assert_freed(&null_lists_heap, null_lists_v3000);
        assert_eq!(null_lists_heap.live_source_allocation_count(), 0);
    }

    #[test]
    fn source_port__inchi_dll__getstructfromstdinchi__line_2461() {
        fn c_string(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .bytes()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect(),
            )
            .unwrap()
        }

        fn invoke(
            heap: &mut SourceHeap,
            input: Option<&inchi_InputINCHI>,
            output: &mut inchi_OutputStruct,
        ) -> Result<i32, SourceHeapError> {
            GetStructFromStdINCHI(
                heap,
                input,
                output,
                SourceMutPointer::null(),
                InchiBuildMetadata {
                    compiler: "gcc",
                    date: "Jan 01 1970",
                    time: "00:00:00",
                },
                1_000,
            )
        }

        fn sentinel_output(heap: &mut SourceHeap) -> inchi_OutputStruct {
            inchi_OutputStruct {
                atom: heap
                    .allocate_model_storage(vec![inchi_Atom::default()])
                    .unwrap(),
                stereo0D: heap
                    .allocate_model_storage(vec![inchi_Stereo0D::default()])
                    .unwrap(),
                num_atoms: i16::MAX,
                num_stereo0D: i16::MIN,
                szMessage: heap.allocate_model_storage(vec![b'M' as i8, 0]).unwrap(),
                szLog: heap.allocate_model_storage(vec![b'L' as i8, 0]).unwrap(),
                WarningFlags: [[u64::MAX, 17], [23, 42]],
            }
        }

        let mut null_heap = SourceHeap::default();
        let unchanged = sentinel_output(&mut null_heap);
        let mut null_output = unchanged.clone();
        assert_eq!(
            invoke(&mut null_heap, None, &mut null_output),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        assert_eq!(null_output, unchanged);

        let null_input = inchi_InputINCHI::default();
        assert_eq!(
            invoke(&mut null_heap, Some(&null_input), &mut null_output),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        assert_eq!(null_output, unchanged);

        for value in ["", "1234567S", "InChI=1/CH4/h1H4", "InChI=1s/CH4/h1H4"] {
            let mut heap = SourceHeap::default();
            let unchanged = sentinel_output(&mut heap);
            let pointer = c_string(&mut heap, value);
            let input = inchi_InputINCHI {
                szInChI: pointer,
                szOptions: SourceMutPointer::null(),
            };
            let mut output = unchanged.clone();
            assert_eq!(
                invoke(&mut heap, Some(&input), &mut output),
                Ok(tagRetValGetINCHI_inchi_Ret_ERROR),
                "{value:?}"
            );
            assert_eq!(output, unchanged, "{value:?}");
        }

        let mut false_positive_heap = SourceHeap::default();
        let unchanged = sentinel_output(&mut false_positive_heap);
        let false_positive_text = c_string(&mut false_positive_heap, "xxxxxxxSx");
        let false_positive_input = inchi_InputINCHI {
            szInChI: false_positive_text,
            szOptions: SourceMutPointer::null(),
        };
        let mut false_positive_output = unchanged.clone();
        assert_eq!(
            invoke(
                &mut false_positive_heap,
                Some(&false_positive_input),
                &mut false_positive_output,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_EOF)
        );
        assert_eq!(false_positive_output, inchi_OutputStruct::default());

        let mut standard_heap = SourceHeap::default();
        let unchanged = sentinel_output(&mut standard_heap);
        let standard_text = c_string(&mut standard_heap, "InChI=1S/CH4/h1H4");
        let standard_input = inchi_InputINCHI {
            szInChI: standard_text,
            szOptions: SourceMutPointer::null(),
        };
        let mut standard_output = unchanged.clone();
        assert_eq!(
            invoke(
                &mut standard_heap,
                Some(&standard_input),
                &mut standard_output,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert_eq!(standard_output.num_atoms, 1);
        assert_eq!(standard_output.num_stereo0D, 0);
        assert_eq!(standard_output.WarningFlags, [[0, 0], [0, 0]]);
        let mut expected_atom = inchi_Atom::default();
        expected_atom.elname[0] = b'C' as i8;
        expected_atom.num_iso_H = [4, 0, 0, 0];
        assert_eq!(
            standard_heap
                .slice(standard_output.atom.as_const())
                .unwrap(),
            &[expected_atom]
        );
        assert!(standard_output.stereo0D.is_null());
        FreeStructFromINCHI(&mut standard_heap, Some(&mut standard_output)).unwrap();

        let mut failure_heap = SourceHeap::default();
        let unchanged = sentinel_output(&mut failure_heap);
        let failure_text = c_string(&mut failure_heap, "InChI=1S/CH4/h1H4");
        let failure_input = inchi_InputINCHI {
            szInChI: failure_text,
            szOptions: SourceMutPointer::null(),
        };
        let mut failure_output = unchanged;
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            invoke(&mut failure_heap, Some(&failure_input), &mut failure_output,),
            Ok(tagRetValGetINCHI_inchi_Ret_EOF)
        );
        assert_eq!(failure_heap.source_allocation_calls(), 1);
        assert_eq!(failure_output, inchi_OutputStruct::default());
    }

    #[test]
    fn source_port__inchi_dll__getstructfrominchi__line_2856() {
        fn c_string(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .bytes()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect(),
            )
            .unwrap()
        }

        fn text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Option<String> {
            if pointer.is_null() {
                return None;
            }
            let bytes = heap.slice(pointer.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            Some(String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap())
        }

        fn invoke(
            heap: &mut SourceHeap,
            input: Option<&inchi_InputINCHI>,
            output: &mut inchi_OutputStruct,
        ) -> Result<i32, SourceHeapError> {
            GetStructFromINCHI(
                heap,
                input,
                output,
                SourceMutPointer::null(),
                InchiBuildMetadata {
                    compiler: "gcc",
                    date: "Jan 01 1970",
                    time: "00:00:00",
                },
                1_000,
            )
        }

        const COMMON_LOG: &str = "Input format: InChI (plain identifier)\n\
Output format: Plain text\n\
Aux. info suppressed\n\
No timeout\n\
Up to 1024 atoms per structure";

        let mut okay_heap = SourceHeap::default();
        let old_atom = okay_heap.allocate(vec![inchi_Atom::default()]).unwrap();
        let old_stereo = okay_heap.allocate(vec![inchi_Stereo0D::default()]).unwrap();
        let old_message = okay_heap.allocate(vec![b'M' as i8, 0]).unwrap();
        let old_log = okay_heap.allocate(vec![b'L' as i8, 0]).unwrap();
        let methane = c_string(&mut okay_heap, "InChI=1S/CH4/h1H4");
        let okay_input = inchi_InputINCHI {
            szInChI: methane,
            szOptions: SourceMutPointer::null(),
        };
        let mut okay_output = inchi_OutputStruct {
            atom: old_atom,
            stereo0D: old_stereo,
            num_atoms: i16::MAX,
            num_stereo0D: i16::MIN,
            szMessage: old_message,
            szLog: old_log,
            WarningFlags: [[u64::MAX; 2]; 2],
        };
        assert_eq!(
            invoke(&mut okay_heap, Some(&okay_input), &mut okay_output),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert!(okay_heap.slice(old_atom.as_const()).is_ok());
        assert!(okay_heap.slice(old_stereo.as_const()).is_ok());
        for pointer in [old_message, old_log] {
            assert!(okay_heap.slice(pointer.as_const()).is_ok());
        }
        assert_eq!(okay_output.num_atoms, 1);
        assert_eq!(okay_output.num_stereo0D, 0);
        assert_eq!(okay_output.WarningFlags, [[0, 0], [0, 0]]);
        assert_eq!(text(&okay_heap, okay_output.szMessage).as_deref(), Some(""));
        assert_eq!(
            text(&okay_heap, okay_output.szLog).as_deref(),
            Some(COMMON_LOG)
        );
        let mut expected_methane = inchi_Atom::default();
        expected_methane.elname[0] = b'C' as i8;
        expected_methane.num_iso_H = [4, 0, 0, 0];
        assert_eq!(
            okay_heap.slice(okay_output.atom.as_const()).unwrap(),
            &[expected_methane]
        );
        assert!(okay_output.stereo0D.is_null());
        FreeStructFromINCHI(&mut okay_heap, Some(&mut okay_output)).unwrap();
        inchi_free(&mut okay_heap, old_atom).unwrap();
        inchi_free(&mut okay_heap, old_stereo).unwrap();
        for pointer in [old_message, old_log] {
            inchi_free(&mut okay_heap, pointer).unwrap();
        }

        let mut warning_heap = SourceHeap::default();
        let warning_text = c_string(
            &mut warning_heap,
            "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1",
        );
        let warning_input = inchi_InputINCHI {
            szInChI: warning_text,
            szOptions: SourceMutPointer::null(),
        };
        let mut warning_output = inchi_OutputStruct::default();
        assert_eq!(
            invoke(&mut warning_heap, Some(&warning_input), &mut warning_output,),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert_eq!(warning_output.num_atoms, 7);
        assert_eq!(warning_output.num_stereo0D, 1);
        assert_eq!(warning_output.WarningFlags, [[0, 0], [0, 0]]);
        assert_eq!(
            text(&warning_heap, warning_output.szMessage).as_deref(),
            Some("")
        );
        assert_eq!(
            text(&warning_heap, warning_output.szLog).as_deref(),
            Some(COMMON_LOG)
        );
        assert_eq!(
            warning_heap
                .slice(warning_output.stereo0D.as_const())
                .unwrap()[0],
            inchi_Stereo0D {
                neighbor: [6, 0, 2, 3],
                central_atom: 1,
                type_: 2,
                parity: 1,
            }
        );
        FreeStructFromINCHI(&mut warning_heap, Some(&mut warning_output)).unwrap();

        for invalid in [None, Some("bad"), Some("InChI=2S/CH4")] {
            let mut heap = SourceHeap::default();
            let input_pointer = invalid
                .map(|value| c_string(&mut heap, value))
                .unwrap_or_else(SourceMutPointer::null);
            let input = inchi_InputINCHI {
                szInChI: input_pointer,
                szOptions: SourceMutPointer::null(),
            };
            let mut output = inchi_OutputStruct {
                num_atoms: 17,
                WarningFlags: [[9; 2]; 2],
                ..inchi_OutputStruct::default()
            };
            assert_eq!(
                invoke(&mut heap, Some(&input), &mut output),
                Ok(tagRetValGetINCHI_inchi_Ret_EOF)
            );
            assert_eq!(output, inchi_OutputStruct::default());
        }

        let mut null_heap = SourceHeap::default();
        let mut null_output = inchi_OutputStruct {
            num_atoms: 17,
            WarningFlags: [[9; 2]; 2],
            ..inchi_OutputStruct::default()
        };
        assert_eq!(
            invoke(&mut null_heap, None, &mut null_output),
            Ok(tagRetValGetINCHI_inchi_Ret_EOF)
        );
        assert_eq!(null_output, inchi_OutputStruct::default());

        let mut failure_heap = SourceHeap::default();
        let failure_text = c_string(&mut failure_heap, "InChI=1S/CH4/h1H4");
        let failure_input = inchi_InputINCHI {
            szInChI: failure_text,
            szOptions: SourceMutPointer::null(),
        };
        let mut failure_output = inchi_OutputStruct::default();
        failure_heap.fail_after_allocations(3);
        assert_eq!(
            invoke(&mut failure_heap, Some(&failure_input), &mut failure_output,),
            Ok(tagRetValGetINCHI_inchi_Ret_FATAL)
        );
        assert_eq!(failure_heap.source_allocation_calls(), 4);
        assert!(failure_output.atom.is_null());
        assert!(failure_output.stereo0D.is_null());
        assert_eq!(failure_output.num_atoms, 0);
        assert_eq!(failure_output.num_stereo0D, 0);
        assert!(failure_output.szMessage.is_null());
        assert_eq!(failure_output.WarningFlags, [[0, 0], [0, 0]]);
        assert_eq!(
            text(&failure_heap, failure_output.szLog),
            Some(format!(
                "{COMMON_LOG}\n\n\nCannot allocate output message buffer."
            ))
        );
        FreeStructFromINCHI(&mut failure_heap, Some(&mut failure_output)).unwrap();

        let mut boundary_heap = SourceHeap::default();
        let boundary_text = c_string(&mut boundary_heap, "InChI=1S/CH4/h1H4");
        let boundary_input = inchi_InputINCHI {
            szInChI: boundary_text,
            szOptions: SourceMutPointer::null(),
        };
        let mut boundary_output = inchi_OutputStruct {
            num_atoms: 17,
            ..inchi_OutputStruct::default()
        };
        boundary_heap.fail_after_allocations(2);
        assert_eq!(
            invoke(
                &mut boundary_heap,
                Some(&boundary_input),
                &mut boundary_output,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(boundary_heap.source_allocation_calls(), 5);
        assert_eq!(boundary_output, inchi_OutputStruct::default());
    }

    #[test]
    fn source_port__inchi_dll__getstructfrominchiex__line_2485() {
        fn c_string(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .bytes()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect(),
            )
            .unwrap()
        }

        fn text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Option<String> {
            if pointer.is_null() {
                return None;
            }
            let bytes = heap.slice(pointer.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            Some(String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap())
        }

        fn free_output(heap: &mut SourceHeap, output: &mut inchi_OutputStructEx) {
            for pointer in [output.szMessage, output.szLog] {
                inchi_free(heap, pointer).unwrap();
            }
            inchi_free(heap, output.atom).unwrap();
            inchi_free(heap, output.stereo0D).unwrap();
            FreeInChIExtInput(heap, output.polymer, output.v3000).unwrap();
            *output = inchi_OutputStructEx::default();
        }

        fn invoke(
            heap: &mut SourceHeap,
            input: Option<&inchi_InputINCHI>,
            output: &mut inchi_OutputStructEx,
        ) -> Result<i32, SourceHeapError> {
            GetStructFromINCHIEx(
                heap,
                input,
                output,
                SourceMutPointer::null(),
                InchiBuildMetadata {
                    compiler: "gcc",
                    date: "Jan 01 1970",
                    time: "00:00:00",
                },
                1_000,
            )
        }

        fn expected_atom(
            element: &str,
            neighbors: &[i16],
            bond_types: &[i8],
            num_iso_h: [i8; 4],
            charge: i8,
        ) -> inchi_Atom {
            assert_eq!(neighbors.len(), bond_types.len());
            let mut atom = inchi_Atom::default();
            for (destination, source) in atom.elname.iter_mut().zip(element.bytes()) {
                *destination = source as i8;
            }
            atom.neighbor[..neighbors.len()].copy_from_slice(neighbors);
            atom.bond_type[..bond_types.len()].copy_from_slice(bond_types);
            atom.num_bonds = neighbors.len() as i16;
            atom.num_iso_H = num_iso_h;
            atom.charge = charge;
            atom
        }

        fn assert_case(
            value: &str,
            expected_result: i32,
            expected_atoms: Vec<inchi_Atom>,
            expected_stereo: Vec<inchi_Stereo0D>,
            expected_message: &str,
            expected_log: &str,
            expected_flags: [[u64; 2]; 2],
        ) {
            let mut heap = SourceHeap::default();
            let pointer = c_string(&mut heap, value);
            let input = inchi_InputINCHI {
                szInChI: pointer,
                szOptions: SourceMutPointer::null(),
            };
            let mut output = inchi_OutputStructEx::default();
            assert_eq!(
                invoke(&mut heap, Some(&input), &mut output),
                Ok(expected_result)
            );
            assert_eq!(output.num_atoms as usize, expected_atoms.len());
            assert_eq!(output.num_stereo0D as usize, expected_stereo.len());
            assert_eq!(heap.slice(output.atom.as_const()).unwrap(), expected_atoms);
            let allocated_stereo = heap.slice(output.stereo0D.as_const()).unwrap_or(&[]);
            assert_eq!(&allocated_stereo[..expected_stereo.len()], expected_stereo);
            assert!(
                allocated_stereo[expected_stereo.len()..]
                    .iter()
                    .all(|stereo| *stereo == inchi_Stereo0D::default())
            );
            assert_eq!(
                text(&heap, output.szMessage).as_deref(),
                Some(expected_message)
            );
            assert_eq!(text(&heap, output.szLog).as_deref(), Some(expected_log));
            assert!(output.polymer.is_null());
            assert!(output.v3000.is_null());
            assert_eq!(output.WarningFlags, expected_flags);
            free_output(&mut heap, &mut output);
        }

        fn assert_allocation_case(
            failure_after: u64,
            expected_result: Result<i32, SourceHeapError>,
            expected_calls: u64,
            expected_atoms: i16,
            expected_message: Option<&str>,
            expected_log: Option<&str>,
        ) {
            let mut heap = SourceHeap::default();
            let pointer = c_string(&mut heap, "InChI=1S/CH4/h1H4");
            let input = inchi_InputINCHI {
                szInChI: pointer,
                szOptions: SourceMutPointer::null(),
            };
            let mut output = inchi_OutputStructEx::default();
            heap.fail_after_allocations(failure_after);
            assert_eq!(
                invoke(&mut heap, Some(&input), &mut output),
                expected_result,
                "allocation ordinal {failure_after}"
            );
            assert_eq!(
                heap.source_allocation_calls(),
                expected_calls,
                "allocation ordinal {failure_after}"
            );
            assert_eq!(
                output.num_atoms, expected_atoms,
                "allocation ordinal {failure_after}"
            );
            assert_eq!(
                text(&heap, output.szMessage).as_deref(),
                expected_message,
                "allocation ordinal {failure_after}"
            );
            assert_eq!(
                text(&heap, output.szLog).as_deref(),
                expected_log,
                "allocation ordinal {failure_after}"
            );
            assert_eq!(output.num_stereo0D, 0);
            assert!(output.polymer.is_null());
            assert!(output.v3000.is_null());
            if expected_atoms == 0 {
                assert!(output.atom.is_null());
                assert!(output.stereo0D.is_null());
            }
            free_output(&mut heap, &mut output);
        }

        const COMMON_LOG: &str = "Input format: InChI (plain identifier)\n\
Output format: Plain text\n\
Aux. info suppressed\n\
No timeout\n\
Up to 1024 atoms per structure";

        let mut null_heap = SourceHeap::default();
        let old_atom = null_heap.allocate(vec![inchi_Atom::default()]).unwrap();
        let old_message = null_heap.allocate(vec![0_i8]).unwrap();
        let mut null_output = inchi_OutputStructEx {
            atom: old_atom,
            num_atoms: 1,
            szMessage: old_message,
            ..inchi_OutputStructEx::default()
        };
        assert_eq!(
            invoke(&mut null_heap, None, &mut null_output),
            Ok(tagRetValGetINCHI_inchi_Ret_EOF)
        );
        assert_eq!(null_output, inchi_OutputStructEx::default());
        assert!(null_heap.slice(old_atom.as_const()).is_ok());
        assert!(null_heap.slice(old_message.as_const()).is_ok());
        inchi_free(&mut null_heap, old_atom).unwrap();
        inchi_free(&mut null_heap, old_message).unwrap();

        for value in [None, Some("bad"), Some("InChI=2S/CH4")] {
            let mut heap = SourceHeap::default();
            let pointer = value
                .map(|value| c_string(&mut heap, value))
                .unwrap_or_else(SourceMutPointer::null);
            let input = inchi_InputINCHI {
                szInChI: pointer,
                szOptions: SourceMutPointer::null(),
            };
            let mut output = inchi_OutputStructEx::default();
            assert_eq!(
                invoke(&mut heap, Some(&input), &mut output),
                Ok(tagRetValGetINCHI_inchi_Ret_EOF),
                "{value:?}"
            );
            assert!(output.atom.is_null());
            assert_eq!(output.num_atoms, 0);
            free_output(&mut heap, &mut output);
        }

        let methane_atom = expected_atom("C", &[], &[], [4, 0, 0, 0], 0);
        for value in ["InChI=1S/CH4/h1H4", "InChI=1/CH4/h1H4"] {
            assert_case(
                value,
                tagRetValGetINCHI_inchi_Ret_OKAY,
                vec![methane_atom.clone()],
                vec![],
                "",
                COMMON_LOG,
                [[0, 0], [0, 0]],
            );
        }

        assert_case(
            "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1",
            tagRetValGetINCHI_inchi_Ret_OKAY,
            vec![
                expected_atom("C", &[1], &[1], [3, 0, 0, 0], 0),
                expected_atom("C", &[6, 0, 2, 3], &[1, 1, 1, 1], [0; 4], 0),
                expected_atom("C", &[1, 4, 5], &[1, 2, 1], [0; 4], 0),
                expected_atom("O", &[1], &[1], [1, 0, 0, 0], 0),
                expected_atom("O", &[2], &[2], [0; 4], 0),
                expected_atom("O", &[2], &[1], [1, 0, 0, 0], 0),
                expected_atom("H", &[1], &[1], [0; 4], 0),
            ],
            vec![inchi_Stereo0D {
                neighbor: [6, 0, 2, 3],
                central_atom: 1,
                type_: 2,
                parity: 1,
            }],
            "",
            COMMON_LOG,
            [[0, 0], [0, 0]],
        );

        assert_case(
            "InChI=1B/C4H4N4.2Zz/c1-5-2-7-4-8-3-6-1;;/h1-4H;;/z101-1-8(9,10-8,3,1,6,2,5,2,7,3,6,1,5,4,7,4,8)/b5-1-,5-2+,6-1+,6-3-,7-2+,7-4+,8-3+,8-4+;;",
            tagRetValGetINCHI_inchi_Ret_OKAY,
            vec![
                expected_atom("C", &[10, 4, 5], &[1, 2, 1], [0; 4], 0),
                expected_atom("C", &[11, 4, 6], &[1, 1, 2], [0; 4], 0),
                expected_atom("C", &[12, 5, 9], &[1, 2, 1], [0; 4], 0),
                expected_atom("C", &[13, 6, 7], &[1, 1, 2], [0; 4], 0),
                expected_atom("N", &[0, 1], &[2, 1], [0; 4], 0),
                expected_atom("N", &[0, 2], &[1, 2], [0; 4], 0),
                expected_atom("N", &[1, 3], &[2, 1], [0; 4], 0),
                expected_atom("N", &[3, 8], &[2, 1], [0; 4], 0),
                expected_atom("*", &[7], &[1], [0; 4], 0),
                expected_atom("*", &[2], &[1], [0; 4], 0),
                expected_atom("H", &[0], &[1], [0; 4], 0),
                expected_atom("H", &[1], &[1], [0; 4], 0),
                expected_atom("H", &[2], &[1], [0; 4], 0),
                expected_atom("H", &[3], &[1], [0; 4], 0),
            ],
            vec![
                inchi_Stereo0D {
                    neighbor: [10, 0, 4, 1],
                    central_atom: -1,
                    type_: 1,
                    parity: 2,
                },
                inchi_Stereo0D {
                    neighbor: [10, 0, 5, 2],
                    central_atom: -1,
                    type_: 1,
                    parity: 1,
                },
                inchi_Stereo0D {
                    neighbor: [11, 1, 4, 0],
                    central_atom: -1,
                    type_: 1,
                    parity: 1,
                },
                inchi_Stereo0D {
                    neighbor: [11, 1, 6, 3],
                    central_atom: -1,
                    type_: 1,
                    parity: 1,
                },
                inchi_Stereo0D {
                    neighbor: [12, 2, 5, 0],
                    central_atom: -1,
                    type_: 1,
                    parity: 2,
                },
                inchi_Stereo0D {
                    neighbor: [13, 3, 6, 1],
                    central_atom: -1,
                    type_: 1,
                    parity: 1,
                },
                inchi_Stereo0D {
                    neighbor: [13, 3, 7, 3],
                    central_atom: -1,
                    type_: 1,
                    parity: 2,
                },
            ],
            "",
            COMMON_LOG,
            [[0, 0], [0, 0]],
        );

        let mut option_failure_heap = SourceHeap::default();
        let option_failure_text = c_string(&mut option_failure_heap, "InChI=1S/CH4/h1H4");
        let option_failure_input = inchi_InputINCHI {
            szInChI: option_failure_text,
            szOptions: SourceMutPointer::null(),
        };
        let mut option_failure_output = inchi_OutputStructEx::default();
        option_failure_heap.fail_after_allocations(0);
        assert_eq!(
            invoke(
                &mut option_failure_heap,
                Some(&option_failure_input),
                &mut option_failure_output,
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_EOF)
        );
        assert_eq!(option_failure_heap.source_allocation_calls(), 1);
        assert_eq!(option_failure_output, inchi_OutputStructEx::default());

        assert_allocation_case(
            1,
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY),
            305,
            1,
            Some(""),
            Some(
                "\nOutput format: Plain text\nAux. info suppressed\nNo timeout\nUp to 1024 atoms per structure",
            ),
        );
        assert_allocation_case(2, Err(SourceHeapError::NullPointer), 5, 0, Some(""), None);
        assert_allocation_case(
            3,
            Ok(tagRetValGetINCHI_inchi_Ret_FATAL),
            4,
            0,
            None,
            Some(&format!(
                "{COMMON_LOG}\n\n\nCannot allocate output message buffer."
            )),
        );
        assert_allocation_case(
            4,
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY),
            304,
            1,
            Some(""),
            Some(COMMON_LOG),
        );
        assert_allocation_case(
            5,
            Ok(tagRetValGetINCHI_inchi_Ret_EOF),
            6,
            0,
            Some(""),
            Some(COMMON_LOG),
        );
        for (ordinal, calls, layer) in [
            (6, 7, "MOBILE_H_FORMULA (0)"),
            (9, 10, "MOBILE_H_CONNECTIONS (1)"),
            (10, 11, "MOBILE_H (2)"),
            (11, 12, "MATERIAL_BALANCE (35)"),
        ] {
            assert_allocation_case(
                ordinal,
                Ok(tagRetValGetINCHI_inchi_Ret_EOF),
                calls,
                0,
                Some(""),
                Some(&format!(
                    "{COMMON_LOG}\n\n\n\nStructure: 1 Allocation failed (-1) in {layer}"
                )),
            );
        }
        for (ordinal, expected_result, calls, message) in [
            (
                18,
                tagRetValGetINCHI_inchi_Ret_FATAL,
                19,
                "*Conversion failed*",
            ),
            (
                19,
                tagRetValGetINCHI_inchi_Ret_ERROR,
                20,
                "*Conversion failed*",
            ),
            (
                300,
                tagRetValGetINCHI_inchi_Ret_FATAL,
                302,
                "Merge Components error; *Conversion failed*",
            ),
            (
                302,
                tagRetValGetINCHI_inchi_Ret_FATAL,
                303,
                "Merge Components error; *Conversion failed*",
            ),
        ] {
            assert_allocation_case(
                ordinal,
                Ok(expected_result),
                calls,
                0,
                Some(message),
                Some(&format!("{COMMON_LOG}\n\n\n{message}")),
            );
        }
        assert_allocation_case(
            303,
            Ok(tagRetValGetINCHI_inchi_Ret_FATAL),
            304,
            0,
            Some(""),
            Some(&format!(
                "{COMMON_LOG}\n\n\nFinal structure conversion failed"
            )),
        );
    }

    #[test]
    fn source_port__inchi_dll__input_erroneously_contains_pseudoatoms__line_293() {
        fn bytes(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Option<Vec<u8>> {
            if pointer.is_null() {
                return None;
            }
            Some(
                heap.slice(pointer.as_const())
                    .unwrap()
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
        }

        let mut heap = SourceHeap::default();
        assert_eq!(
            input_erroneously_contains_pseudoatoms(&mut heap, None, None),
            Err(SourceHeapError::NullPointer)
        );

        let old_inchi = allocate_source_fixture(&mut heap, vec![b'I' as i8, 0]);
        let old_aux = allocate_source_fixture(&mut heap, vec![b'A' as i8, 0]);
        let old_message = allocate_source_fixture(&mut heap, vec![b'M' as i8, 0]);
        let old_log = allocate_source_fixture(&mut heap, vec![b'L' as i8, 0]);
        let unchanged = inchi_Output {
            szInChI: old_inchi,
            szAuxInfo: old_aux,
            szMessage: old_message,
            szLog: old_log,
        };

        for atom_count in [i16::MIN, -1, 0] {
            let input = inchi_Input {
                atom: SourceMutPointer::null(),
                num_atoms: atom_count,
                ..inchi_Input::default()
            };
            let mut output = unchanged.clone();
            assert_eq!(
                input_erroneously_contains_pseudoatoms(&mut heap, Some(&input), Some(&mut output),),
                Ok(0)
            );
            assert_eq!(output, unchanged);
        }

        let atoms = allocate_source_fixture(
            &mut heap,
            vec![
                inchi_Atom {
                    elname: source_name("C"),
                    ..inchi_Atom::default()
                },
                inchi_Atom {
                    elname: source_name("Zz"),
                    ..inchi_Atom::default()
                },
            ],
        );
        let input = inchi_Input {
            atom: atoms,
            num_atoms: 2,
            ..inchi_Input::default()
        };
        let mut output = unchanged.clone();
        assert_eq!(
            input_erroneously_contains_pseudoatoms(&mut heap, Some(&input), Some(&mut output),),
            Ok(0)
        );
        assert_eq!(output, unchanged);

        let star_atoms = allocate_source_fixture(
            &mut heap,
            vec![inchi_Atom {
                elname: source_name("*"),
                ..inchi_Atom::default()
            }],
        );
        let star_input = inchi_Input {
            atom: star_atoms,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        assert_eq!(
            input_erroneously_contains_pseudoatoms(&mut heap, Some(&star_input), None),
            Ok(1)
        );

        let pseudo_atoms = allocate_source_fixture(
            &mut heap,
            vec![inchi_Atom {
                elname: source_name("Zz"),
                ..inchi_Atom::default()
            }],
        );
        let pseudo_input = inchi_Input {
            atom: pseudo_atoms,
            num_atoms: i16::MAX,
            ..inchi_Input::default()
        };
        let live_before_message = heap.live_allocation_count();
        heap.trace_source_allocations();
        let mut success_output = unchanged.clone();
        assert_eq!(
            input_erroneously_contains_pseudoatoms(
                &mut heap,
                Some(&pseudo_input),
                Some(&mut success_output),
            ),
            Ok(1)
        );
        assert_eq!(heap.source_allocation_calls(), 1);
        assert_eq!(heap.live_allocation_count(), live_before_message + 1);
        assert!(success_output.szInChI.is_null());
        assert!(success_output.szAuxInfo.is_null());
        assert!(success_output.szLog.is_null());
        assert_eq!(
            bytes(&heap, success_output.szMessage).as_deref(),
            Some(b"Unsupported in this mode element '*'\0".as_slice())
        );
        for pointer in [old_inchi, old_aux, old_message, old_log] {
            assert!(heap.slice(pointer.as_const()).is_ok());
        }

        let mut failure_heap = SourceHeap::default();
        let failure_atoms = allocate_source_fixture(
            &mut failure_heap,
            vec![inchi_Atom {
                elname: source_name("*"),
                ..inchi_Atom::default()
            }],
        );
        let failure_old = allocate_source_fixture(&mut failure_heap, vec![b'X' as i8, 0]);
        let failure_input = inchi_Input {
            atom: failure_atoms,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        let mut failure_output = inchi_Output {
            szInChI: failure_old,
            szAuxInfo: failure_old,
            szMessage: failure_old,
            szLog: failure_old,
        };
        let failure_live = failure_heap.live_allocation_count();
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            input_erroneously_contains_pseudoatoms(
                &mut failure_heap,
                Some(&failure_input),
                Some(&mut failure_output),
            ),
            Ok(1)
        );
        assert_eq!(failure_output, inchi_Output::default());
        assert_eq!(failure_heap.source_allocation_calls(), 1);
        assert_eq!(failure_heap.live_allocation_count(), failure_live);
        assert_eq!(
            failure_heap.slice(failure_old.as_const()).unwrap(),
            &[b'X' as i8, 0]
        );

        let invalid_input = inchi_Input {
            atom: SourceMutPointer::null(),
            num_atoms: 1,
            ..inchi_Input::default()
        };
        let mut invalid_output = unchanged.clone();
        assert_eq!(
            input_erroneously_contains_pseudoatoms(
                &mut heap,
                Some(&invalid_input),
                Some(&mut invalid_output),
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(invalid_output, unchanged);

        let unterminated = allocate_source_fixture(
            &mut heap,
            vec![inchi_Atom {
                elname: [b'Z' as i8; 6],
                ..inchi_Atom::default()
            }],
        );
        let unterminated_input = inchi_Input {
            atom: unterminated,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        assert_eq!(
            input_erroneously_contains_pseudoatoms(&mut heap, Some(&unterminated_input), None,),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__inchi_dll__produce_generation_output__line_741() {
        fn structure_data(message: &str) -> STRUCT_DATA {
            let mut data = STRUCT_DATA::default();
            for (target, source) in data
                .pStrErrStruct
                .iter_mut()
                .zip(message.bytes().chain(std::iter::once(0)))
            {
                *target = source as i8;
            }
            data
        }

        fn stream(
            heap: &mut SourceHeap,
            text: &[u8],
            padding: usize,
        ) -> (INCHI_IOSTREAM, SourceMutPointer<i8>) {
            let mut bytes = text.iter().map(|byte| *byte as i8).collect::<Vec<_>>();
            bytes.push(0);
            bytes.extend(std::iter::repeat_n(b'#' as i8, padding));
            let pointer = allocate_source_fixture(heap, bytes.clone());
            let mut stream = INCHI_IOSTREAM::default();
            stream.s.pStr = pointer;
            stream.s.nAllocatedLength = bytes.len() as i32;
            stream.s.nUsedLength = text.len() as i32;
            stream.s.nPtr = 64;
            (stream, pointer)
        }

        fn text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Option<String> {
            if pointer.is_null() {
                return None;
            }
            let bytes = heap.slice(pointer.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            Some(String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect()).unwrap())
        }

        let mut heap = SourceHeap::default();
        let mut empty_log = INCHI_IOSTREAM::default();
        let mut empty_output = INCHI_IOSTREAM::default();
        assert_eq!(
            produce_generation_output(
                &mut heap,
                None,
                &STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut empty_log,
                &mut empty_output,
            ),
            Ok(())
        );

        let (mut no_owner_output, no_owner_output_pointer) =
            stream(&mut heap, b"InChI=1S/CH4/h1H4\n", 16);
        let (mut no_owner_log, no_owner_log_pointer) = stream(&mut heap, b"log\n\n", 16);
        assert_eq!(
            produce_generation_output(
                &mut heap,
                None,
                &structure_data("ignored without output"),
                &INPUT_PARMS::default(),
                &mut no_owner_log,
                &mut no_owner_output,
            ),
            Ok(())
        );
        assert_eq!(no_owner_output.s.pStr, no_owner_output_pointer);
        assert_eq!(no_owner_log.s.pStr, no_owner_log_pointer);
        assert_eq!(no_owner_log.s.nUsedLength, 3);
        assert_eq!(text(&heap, no_owner_log_pointer).as_deref(), Some("log"));

        let old_inchi = allocate_source_fixture(&mut heap, vec![b'I' as i8, 0]);
        let old_aux = allocate_source_fixture(&mut heap, vec![b'A' as i8, 0]);
        let old_message = allocate_source_fixture(&mut heap, vec![b'M' as i8, 0]);
        let old_log = allocate_source_fixture(&mut heap, vec![b'L' as i8, 0]);
        let mut output = inchi_Output {
            szInChI: old_inchi,
            szAuxInfo: old_aux,
            szMessage: old_message,
            szLog: old_log,
        };
        let (mut generated, generated_pointer) = stream(
            &mut heap,
            b"InChI=1S/CH4/h1H4\nAuxInfo=1/0/N:1/rA:1C/rB:/rC:;\n",
            16,
        );
        let generated_used = generated.s.nUsedLength;
        let (mut log, log_pointer) = stream(&mut heap, b"prefix structure #9\n\n", 16);
        let mut parameters = INPUT_PARMS::default();
        parameters.bINChIOutputOptions = 0;
        assert_eq!(
            produce_generation_output(
                &mut heap,
                Some(&mut output),
                &structure_data("source warning"),
                &parameters,
                &mut log,
                &mut generated,
            ),
            Ok(())
        );
        assert_eq!(
            text(&heap, output.szMessage).as_deref(),
            Some("source warning")
        );
        assert_eq!(
            text(&heap, output.szInChI).as_deref(),
            Some("InChI=1S/CH4/h1H4")
        );
        assert_eq!(
            text(&heap, output.szAuxInfo).as_deref(),
            Some("AuxInfo=1/0/N:1/rA:1C/rB:/rC:;")
        );
        assert_eq!(output.szInChI, generated_pointer);
        assert_eq!(
            output.szAuxInfo.difference(generated_pointer).unwrap(),
            "InChI=1S/CH4/h1H4\n".len() as i64
        );
        assert!(generated.s.pStr.is_null());
        assert_eq!(generated.s.nUsedLength, generated_used);
        assert_eq!(generated.s.nPtr, 64);
        assert_eq!(output.szLog, log_pointer);
        assert_eq!(text(&heap, output.szLog).as_deref(), Some("prefix"));
        assert!(log.s.pStr.is_null());
        assert_eq!(log.s.nUsedLength, b"prefix structure #9".len() as i32);
        for old in [old_inchi, old_aux, old_message, old_log] {
            assert!(heap.slice(old.as_const()).is_ok());
        }

        let (mut sdf_stream, sdf_pointer) = stream(&mut heap, b"line1\nline2\n", 16);
        let mut sdf_log = INCHI_IOSTREAM::default();
        let mut sdf_output = inchi_Output::default();
        let mut sdf_parameters = INPUT_PARMS::default();
        sdf_parameters.bINChIOutputOptions = INCHI_OUT_SDFILE_ONLY as i32 | 0x4000;
        assert_eq!(
            produce_generation_output(
                &mut heap,
                Some(&mut sdf_output),
                &STRUCT_DATA::default(),
                &sdf_parameters,
                &mut sdf_log,
                &mut sdf_stream,
            ),
            Ok(())
        );
        assert_eq!(sdf_output.szInChI, sdf_pointer);
        assert!(sdf_output.szAuxInfo.is_null());
        assert_eq!(
            text(&heap, sdf_output.szInChI).as_deref(),
            Some("line1\nline2\n")
        );
        assert!(sdf_stream.s.pStr.is_null());

        let (mut multiline_stream, multiline_pointer) = stream(&mut heap, b"one\ntwo\n", 16);
        let mut multiline_log = INCHI_IOSTREAM::default();
        let mut multiline_output = inchi_Output::default();
        assert_eq!(
            produce_generation_output(
                &mut heap,
                Some(&mut multiline_output),
                &STRUCT_DATA::default(),
                &INPUT_PARMS::default(),
                &mut multiline_log,
                &mut multiline_stream,
            ),
            Ok(())
        );
        assert_eq!(multiline_output.szInChI, multiline_pointer);
        assert!(multiline_output.szAuxInfo.is_null());
        assert_eq!(
            text(&heap, multiline_output.szInChI).as_deref(),
            Some("one\ntwo")
        );

        let mut failure_heap = SourceHeap::default();
        let failure_old_message = allocate_source_fixture(&mut failure_heap, vec![b'X' as i8, 0]);
        let (mut failure_stream, failure_pointer) =
            stream(&mut failure_heap, b"InChI=1S/CH4\n", 16);
        let mut failure_log = INCHI_IOSTREAM::default();
        let mut failure_output = inchi_Output {
            szMessage: failure_old_message,
            ..inchi_Output::default()
        };
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            produce_generation_output(
                &mut failure_heap,
                Some(&mut failure_output),
                &structure_data("allocation failure"),
                &INPUT_PARMS::default(),
                &mut failure_log,
                &mut failure_stream,
            ),
            Ok(())
        );
        assert_eq!(failure_heap.source_allocation_calls(), 1);
        assert!(failure_output.szMessage.is_null());
        assert_eq!(failure_output.szInChI, failure_pointer);
        assert_eq!(
            text(&failure_heap, failure_pointer).as_deref(),
            Some("InChI=1S/CH4")
        );
        assert!(failure_stream.s.pStr.is_null());
        assert!(failure_heap.slice(failure_old_message.as_const()).is_ok());

        let mut invalid_data = STRUCT_DATA::default();
        invalid_data.pStrErrStruct.fill(b'x' as i8);
        let mut invalid_log = INCHI_IOSTREAM::default();
        let mut invalid_stream = INCHI_IOSTREAM::default();
        let mut invalid_output = inchi_Output::default();
        assert_eq!(
            produce_generation_output(
                &mut heap,
                Some(&mut invalid_output),
                &invalid_data,
                &INPUT_PARMS::default(),
                &mut invalid_log,
                &mut invalid_stream,
            ),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(invalid_output, inchi_Output::default());
    }

    #[test]
    fn source_port__inchi_dll__copy_corrected_log_tail__line_787() {
        let mut heap = SourceHeap::default();
        let old_log = allocate_source_fixture(&mut heap, vec![b'O' as i8, 0]);
        let old_inchi = allocate_source_fixture(&mut heap, vec![b'I' as i8, 0]);
        let mut output = inchi_Output {
            szInChI: old_inchi,
            szAuxInfo: old_inchi,
            szMessage: old_inchi,
            szLog: old_log,
        };

        let mut null_stream = INCHI_IOSTREAM::default();
        null_stream.s.nUsedLength = 1;
        null_stream.s.nAllocatedLength = 91;
        null_stream.s.nPtr = 17;
        assert_eq!(
            copy_corrected_log_tail(&mut heap, Some(&mut output), &mut null_stream),
            Ok(())
        );
        assert_eq!(output.szLog, old_log);
        assert_eq!(null_stream.s.nUsedLength, 1);

        let dormant = allocate_source_fixture(&mut heap, vec![b'x' as i8, b'\n' as i8, 0]);
        for used in [i32::MIN, -1, 0] {
            let mut stream = INCHI_IOSTREAM::default();
            stream.s.pStr = dormant;
            stream.s.nUsedLength = used;
            let before = heap.slice(dormant.as_const()).unwrap().to_vec();
            let mut untouched_output = output.clone();
            assert_eq!(
                copy_corrected_log_tail(&mut heap, Some(&mut untouched_output), &mut stream,),
                Ok(())
            );
            assert_eq!(untouched_output, output);
            assert_eq!(stream.s.pStr, dormant);
            assert_eq!(stream.s.nUsedLength, used);
            assert_eq!(heap.slice(dormant.as_const()).unwrap(), before);
        }

        let trim_only =
            allocate_source_fixture(&mut heap, vec![b'a' as i8, b'\n' as i8, b'\n' as i8, 0]);
        let mut trim_stream = INCHI_IOSTREAM::default();
        trim_stream.s.pStr = trim_only;
        trim_stream.s.nAllocatedLength = 4;
        trim_stream.s.nUsedLength = 3;
        trim_stream.s.nPtr = 5;
        assert_eq!(
            copy_corrected_log_tail(&mut heap, None, &mut trim_stream),
            Ok(())
        );
        assert_eq!(trim_stream.s.pStr, trim_only);
        assert_eq!(trim_stream.s.nAllocatedLength, 4);
        assert_eq!(trim_stream.s.nUsedLength, 1);
        assert_eq!(trim_stream.s.nPtr, 5);
        assert_eq!(
            heap.slice(trim_only.as_const()).unwrap(),
            &[b'a' as i8, 0, 0, 0]
        );

        let all_lf = allocate_source_fixture(&mut heap, vec![b'\n' as i8, b'\n' as i8, 0]);
        let mut all_lf_stream = INCHI_IOSTREAM::default();
        all_lf_stream.s.pStr = all_lf;
        all_lf_stream.s.nAllocatedLength = 3;
        all_lf_stream.s.nUsedLength = 2;
        let mut all_lf_output = output.clone();
        assert_eq!(
            copy_corrected_log_tail(&mut heap, Some(&mut all_lf_output), &mut all_lf_stream,),
            Ok(())
        );
        assert_eq!(all_lf_stream.s.nUsedLength, 0);
        assert!(all_lf_stream.s.pStr.is_null());
        assert_eq!(all_lf_output.szLog, all_lf);
        assert_eq!(heap.slice(all_lf.as_const()).unwrap(), &[0, 0, 0]);

        let source = b"alpha structure #1 beta Structure #2 gamma structure #3\n\n";
        let mut storage = source.iter().map(|byte| *byte as i8).collect::<Vec<_>>();
        storage.push(0);
        storage.extend(std::iter::repeat_n(b'#' as i8, 16));
        let log = allocate_source_fixture(&mut heap, storage.clone());
        let mut stream = INCHI_IOSTREAM::default();
        stream.s.pStr = log;
        stream.s.nAllocatedLength = storage.len() as i32;
        stream.s.nUsedLength = source.len() as i32;
        stream.s.nPtr = 512;
        let previous_output_log = output.szLog;
        assert_eq!(
            copy_corrected_log_tail(&mut heap, Some(&mut output), &mut stream),
            Ok(())
        );
        assert_eq!(output.szInChI, old_inchi);
        assert_eq!(output.szAuxInfo, old_inchi);
        assert_eq!(output.szMessage, old_inchi);
        assert_eq!(output.szLog, log);
        assert!(stream.s.pStr.is_null());
        assert_eq!(stream.s.nAllocatedLength, storage.len() as i32);
        assert_eq!(stream.s.nUsedLength, source.len() as i32 - 2);
        assert_eq!(stream.s.nPtr, 512);
        assert!(heap.slice(previous_output_log.as_const()).is_ok());
        let bytes = heap.slice(log.as_const()).unwrap();
        let first = b"alpha".len();
        let third = b"alpha structure #1 beta Structure #2 gamma".len();
        assert_eq!(bytes[first], 0);
        assert_eq!(bytes[third], 0);
        assert_eq!(bytes[source.len() - 2], 0);
        assert_eq!(bytes[source.len() - 1], 0);
        let preserved = b" Structure #2"
            .iter()
            .map(|byte| *byte as i8)
            .collect::<Vec<_>>();
        assert_eq!(
            &bytes[b"alpha structure #1 beta".len()..b"alpha structure #1 beta Structure #2".len()],
            preserved
        );

        let short = allocate_source_fixture(&mut heap, vec![b'x' as i8, 0]);
        let mut short_stream = INCHI_IOSTREAM::default();
        short_stream.s.pStr = short;
        short_stream.s.nUsedLength = 5;
        let mut short_output = inchi_Output::default();
        assert_eq!(
            copy_corrected_log_tail(&mut heap, Some(&mut short_output), &mut short_stream,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(short_stream.s.pStr, short);
        assert!(short_output.szLog.is_null());

        let unterminated =
            allocate_source_fixture(&mut heap, vec![b'a' as i8, b'b' as i8, b'c' as i8]);
        let mut unterminated_stream = INCHI_IOSTREAM::default();
        unterminated_stream.s.pStr = unterminated;
        unterminated_stream.s.nUsedLength = 3;
        let mut unterminated_output = inchi_Output::default();
        assert_eq!(
            copy_corrected_log_tail(
                &mut heap,
                Some(&mut unterminated_output),
                &mut unterminated_stream,
            ),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert!(unterminated_stream.s.pStr.is_null());
        assert_eq!(unterminated_output.szLog, unterminated);
    }

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
    fn source_port__inchi_dll__freestructfrominchiex__line_2892() {
        use crate::source_types::inchi_Input_PolymerUnit;

        let mut null_heap = SourceHeap::default();
        assert_eq!(FreeStructFromINCHIEx(&mut null_heap, None), Ok(()));
        let mut empty = inchi_OutputStructEx::default();
        assert_eq!(
            FreeStructFromINCHIEx(&mut null_heap, Some(&mut empty)),
            Ok(())
        );
        assert_eq!(empty, inchi_OutputStructEx::default());

        let mut heap = SourceHeap::default();
        let atom = heap
            .allocate_model_storage(vec![inchi_Atom::default()])
            .unwrap();
        let stereo = heap
            .allocate_model_storage(vec![inchi_Stereo0D::default()])
            .unwrap();
        let log = heap.allocate_model_storage(vec![b'L' as i8, 0]).unwrap();
        let message = heap.allocate_model_storage(vec![b'M' as i8, 0]).unwrap();
        let polymer_unit = heap
            .allocate_model_storage(vec![inchi_Input_PolymerUnit::default()])
            .unwrap();
        let polymer_units = heap.allocate_model_storage(vec![polymer_unit]).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![inchi_Input_Polymer {
                units: polymer_units,
                n: 1,
            }])
            .unwrap();
        let v3000 = heap
            .allocate_model_storage(vec![inchi_Input_V3000::default()])
            .unwrap();
        let mut output = inchi_OutputStructEx {
            atom,
            stereo0D: stereo,
            num_atoms: i16::MAX,
            num_stereo0D: i16::MIN,
            szMessage: message,
            szLog: log,
            WarningFlags: [[u64::MAX, 17], [23, 42]],
            polymer,
            v3000,
        };
        assert_eq!(FreeStructFromINCHIEx(&mut heap, Some(&mut output)), Ok(()));
        assert_eq!(output, inchi_OutputStructEx::default());
        assert_eq!(
            heap.slice(atom.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(stereo.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(log.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(message.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(polymer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(polymer_units.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(polymer_unit.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(v3000.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut polymer_heap = SourceHeap::default();
        let polymer_only = polymer_heap
            .allocate_model_storage(vec![inchi_Input_Polymer::default()])
            .unwrap();
        let mut polymer_output = inchi_OutputStructEx {
            polymer: polymer_only,
            ..inchi_OutputStructEx::default()
        };
        FreeStructFromINCHIEx(&mut polymer_heap, Some(&mut polymer_output)).unwrap();
        assert_eq!(polymer_output, inchi_OutputStructEx::default());
        assert_eq!(
            polymer_heap.slice(polymer_only.as_const()),
            Ok(&[inchi_Input_Polymer::default()][..])
        );
        inchi_free(&mut polymer_heap, polymer_only).unwrap();

        let mut v3000_heap = SourceHeap::default();
        let v3000_only = v3000_heap
            .allocate_model_storage(vec![inchi_Input_V3000::default()])
            .unwrap();
        let mut v3000_output = inchi_OutputStructEx {
            v3000: v3000_only,
            ..inchi_OutputStructEx::default()
        };
        FreeStructFromINCHIEx(&mut v3000_heap, Some(&mut v3000_output)).unwrap();
        assert_eq!(v3000_output, inchi_OutputStructEx::default());
        assert_eq!(
            v3000_heap.slice(v3000_only.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
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

    #[test]
    fn official_c_oracle__rdkit_core_roots__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        crate::source::base::ikey_dll::tests::official_c_oracle__getinchikeyfrominchi__exact();

        fn c_string(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .bytes()
                    .chain(std::iter::once(0))
                    .map(|byte| byte as i8)
                    .collect(),
            )
            .expect("fixture string allocation")
        }

        fn text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Option<String> {
            if pointer.is_null() {
                return None;
            }
            let bytes = heap.slice(pointer.as_const()).expect("source string");
            let end = bytes
                .iter()
                .position(|byte| *byte == 0)
                .expect("source string terminator");
            Some(
                String::from_utf8(bytes[..end].iter().map(|byte| *byte as u8).collect())
                    .expect("source string must be UTF-8"),
            )
        }

        fn expected_text(value: &Value) -> Option<&str> {
            if value.is_null() {
                None
            } else {
                Some(value.as_str().expect("oracle text field"))
            }
        }

        fn expected_i64_array(value: &Value) -> Vec<i64> {
            value
                .as_array()
                .expect("oracle integer array")
                .iter()
                .map(|item| item.as_i64().expect("oracle signed integer"))
                .collect()
        }

        fn assert_atom(case_id: &str, atom: &inchi_Atom, expected: &Value) {
            assert_eq!(
                [atom.x.to_bits(), atom.y.to_bits(), atom.z.to_bits()],
                expected["coordinate_bits"]
                    .as_array()
                    .expect("coordinate bits")
                    .iter()
                    .map(|item| item.as_u64().expect("coordinate bit value"))
                    .collect::<Vec<_>>()
                    .as_slice(),
                "{case_id}: coordinates"
            );
            for (field, actual) in [
                (
                    "neighbor",
                    atom.neighbor
                        .iter()
                        .map(|value| i64::from(*value))
                        .collect::<Vec<_>>(),
                ),
                (
                    "bond_type",
                    atom.bond_type
                        .iter()
                        .map(|value| i64::from(*value))
                        .collect::<Vec<_>>(),
                ),
                (
                    "bond_stereo",
                    atom.bond_stereo
                        .iter()
                        .map(|value| i64::from(*value))
                        .collect::<Vec<_>>(),
                ),
                (
                    "elname",
                    atom.elname
                        .iter()
                        .map(|value| i64::from(*value))
                        .collect::<Vec<_>>(),
                ),
                (
                    "num_iso_h",
                    atom.num_iso_H
                        .iter()
                        .map(|value| i64::from(*value))
                        .collect::<Vec<_>>(),
                ),
            ] {
                assert_eq!(
                    actual,
                    expected_i64_array(&expected[field]),
                    "{case_id}: atom {field}"
                );
            }
            for (field, actual) in [
                ("num_bonds", i64::from(atom.num_bonds)),
                ("isotopic_mass", i64::from(atom.isotopic_mass)),
                ("radical", i64::from(atom.radical)),
                ("charge", i64::from(atom.charge)),
            ] {
                assert_eq!(
                    actual,
                    expected[field].as_i64().expect("atom scalar"),
                    "{case_id}: atom {field}"
                );
            }
        }

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--rdkit-core-root-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let records =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in records.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            let case_id = official["case_id"].as_str().expect("case_id");
            match official["operation"].as_str().expect("operation") {
                "get_inchi_root" => {
                    let kind = official["input"]["kind"].as_i64().expect("fixture kind");
                    let options = expected_text(&official["input"]["options"]);
                    let mut heap = SourceHeap::default();
                    let atom_count = match kind {
                        1 => 5,
                        4 => 2,
                        _ => 1,
                    };
                    let mut atoms = vec![inchi_Atom::default(); atom_count];
                    let mut stereo = Vec::new();
                    match kind {
                        0 => {
                            atoms[0].elname = source_name("C");
                            atoms[0].num_iso_H[0] = -1;
                        }
                        1 => {
                            for (index, element) in ["C", "F", "Cl", "Br", "I"].iter().enumerate() {
                                atoms[index].elname = source_name(element);
                            }
                            for index in 0..4 {
                                atoms[0].neighbor[index] = (index + 1) as i16;
                                atoms[0].bond_type[index] = 1;
                                atoms[index + 1].neighbor[0] = 0;
                                atoms[index + 1].bond_type[0] = 1;
                                atoms[index + 1].num_bonds = 1;
                            }
                            atoms[0].num_bonds = 4;
                            stereo.push(inchi_Stereo0D {
                                neighbor: [1, 2, 3, 4],
                                central_atom: 0,
                                type_: 2,
                                parity: 1,
                            });
                        }
                        2 => {
                            atoms[0].elname = source_name("C");
                            atoms[0].num_iso_H[0] = -1;
                            atoms[0].x = 1.25;
                            atoms[0].y = -2.5;
                        }
                        3 => {
                            atoms[0].elname = source_name("C");
                            atoms[0].num_iso_H[0] = -1;
                            atoms[0].x = 1.25;
                            atoms[0].y = -2.5;
                            atoms[0].z = 3.75;
                            atoms[0].isotopic_mass = 13;
                            atoms[0].charge = 1;
                            atoms[0].radical = 2;
                        }
                        4 => {
                            atoms[0].elname = source_name("C");
                            atoms[1].elname = source_name("Fe");
                            atoms[0].num_iso_H[0] = -1;
                            atoms[0].neighbor[0] = 1;
                            atoms[0].bond_type[0] = 1;
                            atoms[0].num_bonds = 1;
                            atoms[1].neighbor[0] = 0;
                            atoms[1].bond_type[0] = 1;
                            atoms[1].num_bonds = 1;
                        }
                        5 => atoms[0].elname = source_name("*"),
                        _ => panic!("{case_id}: unknown fixture kind {kind}"),
                    }
                    let atoms_snapshot = atoms.clone();
                    let stereo_snapshot = stereo.clone();
                    let atom_pointer = heap
                        .allocate_model_storage(atoms)
                        .expect("atom fixture allocation");
                    let stereo_pointer = if stereo.is_empty() {
                        SourceMutPointer::null()
                    } else {
                        heap.allocate_model_storage(stereo)
                            .expect("stereo fixture allocation")
                    };
                    let option_pointer = options
                        .map(|value| c_string(&mut heap, value))
                        .unwrap_or_else(SourceMutPointer::null);
                    let option_snapshot = if option_pointer.is_null() {
                        None
                    } else {
                        Some(
                            heap.slice(option_pointer.as_const())
                                .expect("option bytes")
                                .to_vec(),
                        )
                    };
                    let input = inchi_Input {
                        atom: atom_pointer,
                        stereo0D: stereo_pointer,
                        szOptions: option_pointer,
                        num_atoms: atom_count as i16,
                        num_stereo0D: stereo_snapshot.len() as i16,
                    };
                    let mut output = inchi_Output::default();
                    let status = GetINCHI(
                        &mut heap,
                        Some(&input),
                        Some(&mut output),
                        0,
                        SourceMutPointer::null(),
                        InchiBuildMetadata {
                            compiler: "gcc",
                            date: "Jan 01 1970",
                            time: "00:00:00",
                        },
                        1_000,
                    )
                    .unwrap_or_else(|error| panic!("{case_id}: GetINCHI failed: {error:?}"));
                    let expected = &official["output"];
                    assert_eq!(
                        i64::from(status),
                        expected["status"].as_i64().expect("status"),
                        "{case_id}: status"
                    );
                    for (field, actual) in [
                        ("inchi", text(&heap, output.szInChI)),
                        ("aux", text(&heap, output.szAuxInfo)),
                        ("message", text(&heap, output.szMessage)),
                        ("log", text(&heap, output.szLog)),
                    ] {
                        assert_eq!(
                            actual.as_deref(),
                            expected_text(&expected[field]),
                            "{case_id}: {field}"
                        );
                    }
                    assert_eq!(
                        heap.slice(atom_pointer.as_const()).expect("input atoms"),
                        atoms_snapshot,
                        "{case_id}: atom input mutation"
                    );
                    if !stereo_pointer.is_null() {
                        assert_eq!(
                            heap.slice(stereo_pointer.as_const()).expect("input stereo"),
                            stereo_snapshot,
                            "{case_id}: stereo input mutation"
                        );
                    }
                    if let Some(snapshot) = option_snapshot {
                        assert_eq!(
                            heap.slice(option_pointer.as_const())
                                .expect("input options"),
                            snapshot,
                            "{case_id}: option input mutation"
                        );
                    }
                    assert_eq!(expected["input_unchanged"], true, "{case_id}");
                    FreeINCHI(&mut heap, Some(&mut output))
                        .unwrap_or_else(|error| panic!("{case_id}: FreeINCHI failed: {error:?}"));
                    assert_eq!(output, inchi_Output::default(), "{case_id}: free reset");
                    for field in ["inchi_null", "aux_null", "message_null", "log_null"] {
                        assert_eq!(expected["after_free"][field], true, "{case_id}: {field}");
                    }
                }
                "get_struct_from_inchi_root" => {
                    let input_text = official["input"]["inchi"].as_str().expect("input InChI");
                    let options = expected_text(&official["input"]["options"]);
                    let mut heap = SourceHeap::default();
                    let input_pointer = c_string(&mut heap, input_text);
                    let option_pointer = options
                        .map(|value| c_string(&mut heap, value))
                        .unwrap_or_else(SourceMutPointer::null);
                    let input_snapshot = heap
                        .slice(input_pointer.as_const())
                        .expect("input bytes")
                        .to_vec();
                    let option_snapshot = if option_pointer.is_null() {
                        None
                    } else {
                        Some(
                            heap.slice(option_pointer.as_const())
                                .expect("option bytes")
                                .to_vec(),
                        )
                    };
                    let input = inchi_InputINCHI {
                        szInChI: input_pointer,
                        szOptions: option_pointer,
                    };
                    let mut output = inchi_OutputStruct {
                        num_atoms: 17,
                        num_stereo0D: -17,
                        WarningFlags: [[u64::MAX; 2]; 2],
                        ..inchi_OutputStruct::default()
                    };
                    let status = GetStructFromINCHI(
                        &mut heap,
                        Some(&input),
                        &mut output,
                        SourceMutPointer::null(),
                        InchiBuildMetadata {
                            compiler: "gcc",
                            date: "Jan 01 1970",
                            time: "00:00:00",
                        },
                        1_000,
                    )
                    .unwrap_or_else(|error| {
                        panic!("{case_id}: GetStructFromINCHI failed: {error:?}")
                    });
                    let expected = &official["output"];
                    assert_eq!(
                        i64::from(status),
                        expected["status"].as_i64().expect("status"),
                        "{case_id}: status"
                    );
                    assert_eq!(
                        i64::from(output.num_atoms),
                        expected["num_atoms"].as_i64().expect("atom count"),
                        "{case_id}: atom count"
                    );
                    assert_eq!(
                        i64::from(output.num_stereo0D),
                        expected["num_stereo"].as_i64().expect("stereo count"),
                        "{case_id}: stereo count"
                    );
                    assert_eq!(
                        text(&heap, output.szMessage).as_deref(),
                        expected_text(&expected["message"]),
                        "{case_id}: message"
                    );
                    assert_eq!(
                        text(&heap, output.szLog).as_deref(),
                        expected_text(&expected["log"]),
                        "{case_id}: log"
                    );
                    assert_eq!(
                        output
                            .WarningFlags
                            .iter()
                            .flatten()
                            .copied()
                            .collect::<Vec<_>>(),
                        expected["warning_flags"]
                            .as_array()
                            .expect("warning flags")
                            .iter()
                            .map(|value| value.as_u64().expect("warning flag"))
                            .collect::<Vec<_>>(),
                        "{case_id}: warning flags"
                    );
                    let expected_atoms = expected["atoms"].as_array().expect("atoms");
                    if output.num_atoms > 0 {
                        let actual_atoms = heap
                            .slice(output.atom.as_const())
                            .expect("output atom storage");
                        assert_eq!(actual_atoms.len(), expected_atoms.len(), "{case_id}");
                        for (atom, expected_atom) in actual_atoms.iter().zip(expected_atoms) {
                            assert_atom(case_id, atom, expected_atom);
                        }
                    } else {
                        assert!(output.atom.is_null(), "{case_id}: atom pointer");
                        assert!(expected_atoms.is_empty(), "{case_id}: expected atoms");
                    }
                    let expected_stereo = expected["stereo"].as_array().expect("stereo");
                    if output.num_stereo0D > 0 {
                        let actual_stereo = heap
                            .slice(output.stereo0D.as_const())
                            .expect("output stereo storage");
                        assert!(actual_stereo.len() >= expected_stereo.len(), "{case_id}");
                        for (actual, expected) in actual_stereo
                            .iter()
                            .take(expected_stereo.len())
                            .zip(expected_stereo)
                        {
                            assert_eq!(
                                actual
                                    .neighbor
                                    .iter()
                                    .map(|value| i64::from(*value))
                                    .collect::<Vec<_>>(),
                                expected_i64_array(&expected["neighbor"]),
                                "{case_id}: stereo neighbors"
                            );
                            for (field, value) in [
                                ("central_atom", i64::from(actual.central_atom)),
                                ("type", i64::from(actual.type_)),
                                ("parity", i64::from(actual.parity)),
                            ] {
                                assert_eq!(
                                    value,
                                    expected[field].as_i64().expect("stereo scalar"),
                                    "{case_id}: stereo {field}"
                                );
                            }
                        }
                        assert!(
                            actual_stereo[expected_stereo.len()..]
                                .iter()
                                .all(|stereo| *stereo == inchi_Stereo0D::default()),
                            "{case_id}: unused Rust stereo allocation slots"
                        );
                    } else {
                        assert!(output.stereo0D.is_null(), "{case_id}: stereo pointer");
                        assert!(expected_stereo.is_empty(), "{case_id}: expected stereo");
                    }
                    assert_eq!(
                        heap.slice(input_pointer.as_const()).expect("input bytes"),
                        input_snapshot,
                        "{case_id}: input mutation"
                    );
                    if let Some(snapshot) = option_snapshot {
                        assert_eq!(
                            heap.slice(option_pointer.as_const()).expect("option bytes"),
                            snapshot,
                            "{case_id}: option mutation"
                        );
                    }
                    assert_eq!(expected["input_unchanged"], true, "{case_id}");
                    FreeStructFromINCHI(&mut heap, Some(&mut output)).unwrap_or_else(|error| {
                        panic!("{case_id}: FreeStructFromINCHI failed: {error:?}")
                    });
                    assert_eq!(
                        output,
                        inchi_OutputStruct::default(),
                        "{case_id}: free reset"
                    );
                    let after_free = &expected["after_free"];
                    for field in ["atom_null", "stereo_null", "message_null", "log_null"] {
                        assert_eq!(after_free[field], true, "{case_id}: {field}");
                    }
                    assert_eq!(after_free["num_atoms"], 0, "{case_id}");
                    assert_eq!(after_free["num_stereo"], 0, "{case_id}");
                    assert_eq!(
                        after_free["warning_flags"],
                        serde_json::json!([0, 0, 0, 0]),
                        "{case_id}"
                    );
                }
                operation => panic!("{case_id}: unexpected operation {operation}"),
            }
            record_count += 1;
        }
        assert_eq!(record_count, 15);
    }
}
