use crate::source::api::inchi_dll_b::{
    CreateInchiAtom, FreeInchi_Atom, FreeInchi_Input, InchiToInchiAtom, add_source_error,
};
use crate::source::base::ichi_io::inchi_ios_init;
use crate::source::base::readinch::{CreateInchi_Stereo0D, FreeInchi_Stereo0D};
use crate::source::base::util::inchi_free;
use crate::source_types::{
    AB_PARITY_UNDF, AB_PARITY_UNKN, INCHI_IOS_TYPE_STRING, INCHI_IOSTREAM, INCHI_MODE, INPUT_TYPE,
    InchiInpData, MAX_ATOMS, MAX_SDF_HEADER, MAX_SDF_VALUE, STR_ERR_LEN, SourceHeap,
    SourceHeapError, SourceMutPointer, inchi_Input, tagInputType_INPUT_INCHI_PLAIN,
    tagRetValGetINCHI_inchi_Ret_EOF, tagRetValGetINCHI_inchi_Ret_ERROR,
    tagRetValGetINCHI_inchi_Ret_FATAL, tagRetValGetINCHI_inchi_Ret_OKAY,
    tagRetValGetINCHI_inchi_Ret_WARNING,
};

pub(crate) fn Get_std_inchi_Input_FromAuxInfo(
    heap: &mut SourceHeap,
    aux_info: SourceMutPointer<i8>,
    do_not_add_h: i32,
    output: Option<&mut InchiInpData>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:79 Get_std_inchi_Input_FromAuxInfo
    // INCHI✔❌: EXPIMP_TEMPLATE INCHI_API int INCHI_DECL Get_std_inchi_Input_FromAuxInfo
    // INCHI✔❌: ( char *szInchiAuxInfo,
    // INCHI✔❌: int bDoNotAddH,
    // INCHI✔❌: InchiInpData *pInchiInp )
    // INCHI✔❌: {
    // INCHI✔❌:     int bDiffUnkUndfStereo = 0;
    // INCHI✔❌:     return Get_inchi_Input_FromAuxInfo( szInchiAuxInfo, bDoNotAddH, bDiffUnkUndfStereo,
    // INCHI✔❌:                                         pInchiInp );
    // INCHI✔❌: }
    // END INCHI C FUNCTION: Get_std_inchi_Input_FromAuxInfo

    let distinguish_unknown_undefined_stereo = 0;
    Get_inchi_Input_FromAuxInfo(
        heap,
        aux_info,
        do_not_add_h,
        distinguish_unknown_undefined_stereo,
        output,
    )
}

pub(crate) fn Get_inchi_Input_FromAuxInfo(
    heap: &mut SourceHeap,
    aux_info: SourceMutPointer<i8>,
    do_not_add_h: i32,
    distinguish_unknown_undefined_stereo: i32,
    output: Option<&mut InchiInpData>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:89 Get_inchi_Input_FromAuxInfo
    // INCHI✔❌: EXPIMP_TEMPLATE INCHI_API int INCHI_DECL Get_inchi_Input_FromAuxInfo( char *szInchiAuxInfo,
    // INCHI✔❌:                                                   int bDoNotAddH, int bDiffUnkUndfStereo,
    // INCHI✔❌:                                                   InchiInpData *pInchiInp )
    // INCHI✔❌: {
    // INCHI✔❌:     INCHI_IOSTREAM inp;
    // INCHI✔❌:     int num_at, nRet = inchi_Ret_OKAY, err = 0;
    // INCHI✔❌:     INCHI_MODE bChiral = 0;
    // INCHI✔❌:     /* the input string may contain the following line: "Structure NNN. HHH=VVV" */
    // INCHI✔❌:     long         lNumber;                   /* structure number NNN from the input */
    // INCHI✔❌:     char         szHeader[MAX_SDF_HEADER];  /* stucture label header HHH from the input */
    // INCHI✔❌:     char         szLabel[MAX_SDF_VALUE];    /* stucture label VVV from the input */
    // INCHI✔❌:
    // INCHI✔❌:     /* vABParityUnknown holds actual value of an internal constant signifying       */
    // INCHI✔❌:     /* unknown parity: either the same as for undefined parity (default==standard)  */
    // INCHI✔❌:     /*  or a specific one (non-std; requested by SLUUD switch).                     */
    // INCHI✔❌:     int vABParityUnknown = AB_PARITY_UNDF;
    // INCHI✔❌:     if (0 != bDiffUnkUndfStereo)
    // INCHI✔❌:     {
    // INCHI✔❌:         /* Make labels for unknown and undefined stereo different */
    // INCHI✔❌:         vABParityUnknown = AB_PARITY_UNKN;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:     if (pInchiInp && pInchiInp->pInp)
    // INCHI✔❌:     {
    // INCHI✔❌: /* clear output fields */
    // INCHI✔❌:         inchi_Input *pInp = pInchiInp->pInp;
    // INCHI✔❌:         char        *szOptions = pInp->szOptions;
    // INCHI✔❌:         memset( pInchiInp, 0, sizeof( *pInchiInp ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:         memset( pInp, 0, sizeof( *pInp ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔❌:         pInp->szOptions = szOptions;
    // INCHI✔❌:         pInchiInp->pInp = pInp;
    // INCHI✔❌:     }
    // INCHI✔❌:     else
    // INCHI✔❌:     {
    // INCHI✔❌:         return inchi_Ret_ERROR;
    // INCHI✔❌:     }
    // INCHI✔❌:     szHeader[0] = '\0';
    // INCHI✔❌:     szLabel[0] = '\0';
    // INCHI✔❌:     lNumber = 0;
    // INCHI✔❌:     /* prepare input string pointers */
    // INCHI✔❌:     inchi_ios_init( &inp, INCHI_IOS_TYPE_STRING, NULL ); /* fix bug discovered by Burt Leland 2008-12-23 */
    // INCHI✔❌:     inp.s.pStr = szInchiAuxInfo;
    // INCHI✔❌:     inp.s.nUsedLength = (int) strlen( szInchiAuxInfo );
    // INCHI✔❌:     inp.s.nAllocatedLength = inp.s.nUsedLength + 1;
    // INCHI✔❌:     inp.s.nPtr = 0;
    // INCHI✔❌:
    // INCHI✔❌:     num_at = InchiToInchi_Input( &inp, pInchiInp->pInp, 1, bDoNotAddH, vABParityUnknown,
    // INCHI✔❌:                                  INPUT_INCHI_PLAIN, szHeader, szLabel,
    // INCHI✔❌:                                  &lNumber, &bChiral, &err, pInchiInp->szErrMsg );
    // INCHI✔❌:     pInchiInp->bChiral = bChiral;
    // INCHI✔❌:     if (num_at <= 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (10 < err && err < 20)
    // INCHI✔❌:         {
    // INCHI✔❌:             nRet = inchi_Ret_EOF;
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:             if (err == 9)
    // INCHI✔❌:             {
    // INCHI✔❌:                 nRet = inchi_Ret_ERROR; /*  sdfile bypassed to $$$$ */
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:                 if (err && err < 30)
    // INCHI✔❌:                 {
    // INCHI✔❌:                     nRet = inchi_Ret_FATAL;
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                     if (98 == err)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         nRet = inchi_Ret_WARNING; /* empty AuxInfo */
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else
    // INCHI✔❌:                         if (err)
    // INCHI✔❌:                         {
    // INCHI✔❌:                             nRet = inchi_Ret_ERROR;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                             if (pInchiInp->szErrMsg[0])
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 nRet = inchi_Ret_WARNING;
    // INCHI✔❌:                             }
    // INCHI✔❌:     }
    // INCHI✔❌:     if (nRet != inchi_Ret_OKAY && nRet != inchi_Ret_WARNING)
    // INCHI✔❌:     {
    // INCHI✔❌:         Free_inchi_Input( pInchiInp->pInp );
    // INCHI✔❌:         pInchiInp->bChiral = 0;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return nRet;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: Get_inchi_Input_FromAuxInfo

    let parity_unknown = if distinguish_unknown_undefined_stereo != 0 {
        AB_PARITY_UNKN as i32
    } else {
        AB_PARITY_UNDF as i32
    };
    let Some(output) = output else {
        return Ok(tagRetValGetINCHI_inchi_Ret_ERROR);
    };
    if output.pInp.is_null() {
        return Ok(tagRetValGetINCHI_inchi_Ret_ERROR);
    }

    let input_pointer = output.pInp;
    let options = heap
        .slice(input_pointer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .szOptions;
    *output = InchiInpData::default();
    output.pInp = input_pointer;
    let mut input_data = inchi_Input {
        szOptions: options,
        ..inchi_Input::default()
    };
    *heap
        .slice_mut(input_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = input_data.clone();

    let aux_length = heap
        .slice(aux_info.as_const())?
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let aux_length_i32 =
        i32::try_from(aux_length).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let allocated_length = aux_length_i32
        .checked_add(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;

    let header = heap.allocate_model_storage(vec![0_i8; MAX_SDF_HEADER as usize])?;
    let label = heap.allocate_model_storage(vec![0_i8; MAX_SDF_VALUE as usize])?;
    let error_buffer = heap.allocate_model_storage(vec![0_i8; STR_ERR_LEN as usize])?;
    let mut stream = INCHI_IOSTREAM::default();
    inchi_ios_init(
        Some(&mut stream),
        INCHI_IOS_TYPE_STRING as i32,
        SourceMutPointer::null(),
    )?;
    stream.s.pStr = aux_info;
    stream.s.nUsedLength = aux_length_i32;
    stream.s.nAllocatedLength = allocated_length;
    stream.s.nPtr = 0;

    let mut structure_number = 0_i64;
    let mut chiral = 0_u64;
    let mut error = 0_i32;
    let parsed = InchiToInchi_Input(
        heap,
        Some(&mut stream),
        Some(&mut input_data),
        1,
        do_not_add_h,
        parity_unknown,
        tagInputType_INPUT_INCHI_PLAIN,
        Some(header),
        Some(label),
        Some(&mut structure_number),
        Some(&mut chiral),
        Some(&mut error),
        Some(error_buffer),
    );

    let operation = parsed.and_then(|atom_count| {
        output.bChiral =
            i32::try_from(chiral).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        output.szErrMsg.copy_from_slice(
            heap.slice(error_buffer.as_const())?
                .get(..STR_ERR_LEN as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );

        let mut return_value = tagRetValGetINCHI_inchi_Ret_OKAY;
        if atom_count <= 0 {
            return_value = if 10 < error && error < 20 {
                tagRetValGetINCHI_inchi_Ret_EOF
            } else if error == 9 {
                tagRetValGetINCHI_inchi_Ret_ERROR
            } else if error != 0 && error < 30 {
                tagRetValGetINCHI_inchi_Ret_FATAL
            } else if error == 98 {
                tagRetValGetINCHI_inchi_Ret_WARNING
            } else if error != 0 {
                tagRetValGetINCHI_inchi_Ret_ERROR
            } else if output.szErrMsg[0] != 0 {
                tagRetValGetINCHI_inchi_Ret_WARNING
            } else {
                tagRetValGetINCHI_inchi_Ret_OKAY
            };
        }

        if return_value != tagRetValGetINCHI_inchi_Ret_OKAY
            && return_value != tagRetValGetINCHI_inchi_Ret_WARNING
        {
            Free_inchi_Input(heap, &mut input_data)?;
            output.bChiral = 0;
        }
        Ok(return_value)
    });

    *heap
        .slice_mut(input_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = input_data;
    let header_cleanup = inchi_free(heap, header);
    let label_cleanup = inchi_free(heap, label);
    let error_cleanup = inchi_free(heap, error_buffer);
    let cleanup = header_cleanup
        .and(label_cleanup)
        .and(error_cleanup)
        .map(|_| ());
    match operation {
        Ok(return_value) => cleanup.map(|_| return_value),
        Err(source_error) => Err(source_error),
    }
}

pub(crate) fn Free_inchi_Input(
    heap: &mut SourceHeap,
    input: &mut inchi_Input,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:187 Free_inchi_Input
    // INCHI✔️❌: void INCHI_DECL Free_inchi_Input( inchi_Input *pInp )
    // INCHI✔️❌: {
    // INCHI✔️❌:     FreeInchi_Atom( &pInp->atom );
    // INCHI✔️❌:     FreeInchi_Stereo0D( &pInp->stereo0D );
    // INCHI✔️❌:     pInp->num_atoms = 0;
    // INCHI✔️❌:     pInp->num_stereo0D = 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: Free_inchi_Input

    FreeInchi_Atom(heap, Some(&mut input.atom))?;
    FreeInchi_Stereo0D(heap, Some(&mut input.stereo0D))?;
    input.num_atoms = 0;
    input.num_stereo0D = 0;
    Ok(())
}

pub(crate) fn Free_std_inchi_Input(
    heap: &mut SourceHeap,
    input: &mut inchi_Input,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:182 Free_std_inchi_Input
    // INCHI✔️✔️: void INCHI_DECL Free_std_inchi_Input( inchi_Input *pInp )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     Free_inchi_Input( pInp );
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: Free_std_inchi_Input

    Free_inchi_Input(heap, input)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn InchiToInchi_Input(
    heap: &mut SourceHeap,
    input: Option<&mut INCHI_IOSTREAM>,
    mut orig_at_data: Option<&mut inchi_Input>,
    merge_all_input_structures: i32,
    do_not_add_h: i32,
    ab_parity_unknown: i32,
    input_type: INPUT_TYPE,
    sdf_label: Option<SourceMutPointer<i8>>,
    sdf_value: Option<SourceMutPointer<i8>>,
    mut sdf_id: Option<&mut i64>,
    mut input_atom_flags: Option<&mut INCHI_MODE>,
    error: Option<&mut i32>,
    error_text: Option<SourceMutPointer<i8>>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:204 InchiToInchi_Input
    // INCHI✔❌: int InchiToInchi_Input( INCHI_IOSTREAM *inp_molfile,
    // INCHI✔❌:                         inchi_Input *orig_at_data,
    // INCHI✔❌:                         int bMergeAllInputStructures,
    // INCHI✔❌:                         int bDoNotAddH,
    // INCHI✔❌:                         int vABParityUnknown,
    // INCHI✔❌:                         INPUT_TYPE nInputType,
    // INCHI✔❌:                         char *pSdfLabel,
    // INCHI✔❌:                         char *pSdfValue,
    // INCHI✔❌:                         long *lSdfId,
    // INCHI✔❌:                         INCHI_MODE *pInpAtomFlags,
    // INCHI✔❌:                         int *err,
    // INCHI✔❌:                         char *pStrErr )
    // INCHI✔❌: {
    // INCHI✔❌:     /* inp_ATOM       *at = NULL; */
    // INCHI✔❌:     int             num_dimensions_new;
    // INCHI✔❌:     int             num_inp_bonds_new;
    // INCHI✔❌:     int             num_inp_atoms_new;
    // INCHI✔❌:     int             num_inp_0D_new;
    // INCHI✔❌:     inchi_Atom     *at_new = NULL;
    // INCHI✔❌:     inchi_Atom     *at_old = NULL;
    // INCHI✔❌:     inchi_Stereo0D *stereo0D_new = NULL;
    // INCHI✔❌:     inchi_Stereo0D *stereo0D_old = NULL;
    // INCHI✔❌:     int             nNumAtoms = 0, nNumStereo0D = 0;
    // INCHI✔❌:     MOL_COORD      *szCoordNew = NULL;
    // INCHI✔❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔❌:     int            i, j;
    // INCHI✔❌:
    // INCHI✔❌:     if (pStrErr)
    // INCHI✔❌:     {
    // INCHI✔❌:         pStrErr[0] = '\0';
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     /*FreeOrigAtData( orig_at_data );*/
    // INCHI✔❌:     if (lSdfId)
    // INCHI✔❌:         *lSdfId = 0;
    // INCHI✔❌:     do
    // INCHI✔❌:     {
    // INCHI✔❌:
    // INCHI✔❌:         at_old = orig_at_data ? orig_at_data->atom : NULL; /*  save pointer to the previous allocation */
    // INCHI✔❌:         stereo0D_old = orig_at_data ? orig_at_data->stereo0D : NULL;
    // INCHI✔❌:         /* djb-rwth: removing redundant code */
    // INCHI✔❌:         num_inp_atoms_new =
    // INCHI✔❌:             InchiToInchiAtom( inp_molfile, orig_at_data ? &stereo0D_new : NULL, &num_inp_0D_new,
    // INCHI✔❌:                           bDoNotAddH, vABParityUnknown, nInputType,
    // INCHI✔❌:                           orig_at_data ? &at_new : NULL, MAX_ATOMS,
    // INCHI✔❌:                           &num_dimensions_new, &num_inp_bonds_new,
    // INCHI✔❌:                           pSdfLabel, pSdfValue, lSdfId, pInpAtomFlags, err, pStrErr );
    // INCHI✔❌:         if (num_inp_atoms_new <= 0 && !*err)
    // INCHI✔❌:         {
    // INCHI✔❌:             TREAT_ERR( *err, 0, "Empty structure" );
    // INCHI✔❌:             *err = 98;
    // INCHI✔❌:         }
    // INCHI✔❌:         else
    // INCHI✔❌:             if (orig_at_data && !num_inp_atoms_new && 10 < *err && *err < 20 && orig_at_data->num_atoms > 0 && bMergeAllInputStructures)
    // INCHI✔❌:             {
    // INCHI✔❌:                 *err = 0; /* end of file */
    // INCHI✔❌:                 break;
    // INCHI✔❌:             }
    // INCHI✔❌:             else
    // INCHI✔❌:                 if (num_inp_atoms_new > 0 && orig_at_data)
    // INCHI✔❌:                 {
    // INCHI✔❌: /*  merge pOrigDataTmp + orig_at_data => pOrigDataTmp; */
    // INCHI✔❌:                     nNumAtoms = num_inp_atoms_new + orig_at_data->num_atoms;
    // INCHI✔❌:                     nNumStereo0D = num_inp_0D_new + orig_at_data->num_stereo0D;
    // INCHI✔❌:                     if (nNumAtoms >= MAX_ATOMS)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         TREAT_ERR( *err, 0, "Too many atoms [did you forget 'LargeMolecules' switch?]" );
    // INCHI✔❌:                         *err = 70;
    // INCHI✔❌:                         orig_at_data->num_atoms = -1;
    // INCHI✔❌:                     }
    // INCHI✔❌:                     else
    // INCHI✔❌:                         if (!at_old)
    // INCHI✔❌:                         {
    // INCHI✔❌:             /* the first structure */
    // INCHI✔❌:                             orig_at_data->atom = at_new;
    // INCHI✔❌:                             at_new = NULL;
    // INCHI✔❌:                             orig_at_data->num_atoms = num_inp_atoms_new;
    // INCHI✔❌:                             num_inp_atoms_new = 0; /* djb-rwth: ignoring LLVM warning: variable value used */
    // INCHI✔❌:                             orig_at_data->stereo0D = stereo0D_new;
    // INCHI✔❌:                             stereo0D_new = NULL;
    // INCHI✔❌:                             orig_at_data->num_stereo0D = num_inp_0D_new;
    // INCHI✔❌:                             num_inp_0D_new = 0;
    // INCHI✔❌:                         }
    // INCHI✔❌:                         else
    // INCHI✔❌:                             if ((orig_at_data->atom = CreateInchiAtom( nNumAtoms ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                             {
    // INCHI✔❌: /*  switch at_new <--> orig_at_data->at; */
    // INCHI✔❌:                                 if (orig_at_data->num_atoms)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     memcpy( orig_at_data->atom, at_old, orig_at_data->num_atoms * sizeof( orig_at_data->atom[0] ) );
    // INCHI✔❌:                                     /*  adjust numbering in the newly read structure */
    // INCHI✔❌:                                     for (i = 0; i < num_inp_atoms_new; i++)
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         for (j = 0; j < at_new[i].num_bonds; j++)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             at_new[i].neighbor[j] += orig_at_data->num_atoms;
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 FreeInchi_Atom( &at_old );
    // INCHI✔❌:                                 /*  copy newly read structure */
    // INCHI✔❌:                                 if (at_new) /* djb-rwth: fixing a NULL pointer dereference */
    // INCHI✔❌:                                     memcpy(orig_at_data->atom + orig_at_data->num_atoms, at_new, num_inp_atoms_new * sizeof(orig_at_data->atom[0]));
    // INCHI✔❌:                                 /*  copy newly read 0D stereo */
    // INCHI✔❌:                                 if (num_inp_0D_new > 0 && stereo0D_new)
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     if ((orig_at_data->stereo0D = CreateInchi_Stereo0D( nNumStereo0D ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         memcpy(orig_at_data->stereo0D, stereo0D_old, orig_at_data->num_stereo0D * sizeof(orig_at_data->stereo0D[0]));
    // INCHI✔❌:                                         /*  adjust numbering in the newly read structure */
    // INCHI✔❌:                                         for (i = 0; i < num_inp_0D_new; i++)
    // INCHI✔❌:                                         {
    // INCHI✔❌:                                             if (stereo0D_new[i].central_atom >= 0)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 stereo0D_new[i].central_atom += orig_at_data->num_atoms;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                             for (j = 0; j < 4; j++)
    // INCHI✔❌:                                             {
    // INCHI✔❌:                                                 stereo0D_new[i].neighbor[j] += orig_at_data->num_atoms;
    // INCHI✔❌:                                             }
    // INCHI✔❌:                                         }
    // INCHI✔❌:                                         FreeInchi_Stereo0D( &stereo0D_old );
    // INCHI✔❌:                                         memcpy(orig_at_data->stereo0D + orig_at_data->num_stereo0D,
    // INCHI✔❌:                                             stereo0D_new,
    // INCHI✔❌:                                             num_inp_0D_new * sizeof(orig_at_data->stereo0D[0]));
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                     else
    // INCHI✔❌:                                     {
    // INCHI✔❌:                                         num_inp_0D_new = 0;
    // INCHI✔❌:                                         TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔❌:                                         *err = -1;
    // INCHI✔❌:                                     }
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 else
    // INCHI✔❌:                                 {
    // INCHI✔❌:                                     num_inp_0D_new = 0;
    // INCHI✔❌:                                 }
    // INCHI✔❌:                                 /* update lengths */
    // INCHI✔❌:                                 orig_at_data->num_atoms += num_inp_atoms_new;
    // INCHI✔❌:                                 orig_at_data->num_stereo0D += num_inp_0D_new;
    // INCHI✔❌:                             }
    // INCHI✔❌:                             else
    // INCHI✔❌:                             {
    // INCHI✔❌:                                 TREAT_ERR( *err, 0, "Out of RAM" );
    // INCHI✔❌:                                 *err = -1;
    // INCHI✔❌:                             }
    // INCHI✔❌:                 }
    // INCHI✔❌:                 else
    // INCHI✔❌:                     if (num_inp_atoms_new > 0)
    // INCHI✔❌:                     {
    // INCHI✔❌:                         nNumAtoms += num_inp_atoms_new;
    // INCHI✔❌:                     }
    // INCHI✔❌:         FreeInchi_Atom( &at_new );
    // INCHI✔❌:         /* djb-rwth: removing redundant code */
    // INCHI✔❌:         FreeInchi_Stereo0D( &stereo0D_new );
    // INCHI✔❌:         num_inp_0D_new = 0;
    // INCHI✔❌:     }
    // INCHI✔❌:     while (!*err && bMergeAllInputStructures);
    // INCHI✔❌:  /*
    // INCHI✔❌:  if ( !*err ) {
    // INCHI✔❌:      orig_at_data->num_components =
    // INCHI✔❌:          MarkDisconnectedComponents( orig_at_data );
    // INCHI✔❌:      if ( orig_at_data->num_components == 0 ) {
    // INCHI✔❌:          TREAT_ERR (*err, 0, "No components found");
    // INCHI✔❌:          *err = 99;
    // INCHI✔❌:      }
    // INCHI✔❌:      if ( orig_at_data->num_components < 0 ) {
    // INCHI✔❌:          TREAT_ERR (*err, 0, "Too many components");
    // INCHI✔❌:          *err = 99;
    // INCHI✔❌:      }
    // INCHI✔❌:  }
    // INCHI✔❌:  */
    // INCHI✔❌:     if (szCoordNew)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_free( szCoordNew );
    // INCHI✔❌:     }
    // INCHI✔❌:     if (at_new)
    // INCHI✔❌:     {
    // INCHI✔❌:         inchi_free( at_new );
    // INCHI✔❌:     }
    // INCHI✔❌:     /*
    // INCHI✔❌:     if ( !*err ) {
    // INCHI✔❌:         if ( ReconcileAllCmlBondParities( orig_at_data->atom, orig_at_data->num_atoms ) ) {
    // INCHI✔❌:             TREAT_ERR (*err, 0, "Cannot reconcile stereobond parities");
    // INCHI✔❌:             if (!orig_at_data->num_atoms) {
    // INCHI✔❌:                 *err = 1;
    // INCHI✔❌:             }
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:     */
    // INCHI✔❌:     if (*err)
    // INCHI✔❌:     {
    // INCHI✔❌:         FreeInchi_Input( orig_at_data );
    // INCHI✔❌:     }
    // INCHI✔❌:     if (*err && !( 10 < *err && *err < 20 ) && pStrErr && !pStrErr[0])
    // INCHI✔❌:     {
    // INCHI✔❌:         TREAT_ERR( *err, 0, "Unknown error" );  /*   <BRKPT> */
    // INCHI✔❌:     }
    // INCHI✔❌:     return orig_at_data ? orig_at_data->num_atoms : nNumAtoms;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: InchiToInchi_Input

    if let Some(error_text) = error_text {
        *heap
            .slice_mut(error_text)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    }
    if let Some(sdf_id) = sdf_id.as_deref_mut() {
        *sdf_id = 0;
    }

    let input = input.ok_or(SourceHeapError::NullPointer)?;
    let error = error.ok_or(SourceHeapError::NullPointer)?;
    let mut atom_new = SourceMutPointer::null();
    let mut stereo_new = SourceMutPointer::null();
    let mut accumulated_atom_count = 0_i32;

    loop {
        let atom_old = orig_at_data
            .as_deref()
            .map_or(SourceMutPointer::null(), |input| input.atom);
        let stereo_old = orig_at_data
            .as_deref()
            .map_or(SourceMutPointer::null(), |input| input.stereo0D);
        let mut new_stereo_count = 0_i32;
        let mut new_dimensions = 0_i32;
        let mut new_bond_count = 0_i32;
        let output_requested = orig_at_data.is_some();

        let new_atom_count = InchiToInchiAtom(
            heap,
            Some(&mut *input),
            output_requested.then_some(&mut stereo_new),
            Some(&mut new_stereo_count),
            do_not_add_h,
            ab_parity_unknown,
            input_type,
            output_requested.then_some(&mut atom_new),
            MAX_ATOMS as i32,
            Some(&mut new_dimensions),
            Some(&mut new_bond_count),
            sdf_label,
            sdf_value,
            sdf_id.as_deref_mut(),
            input_atom_flags.as_deref_mut(),
            Some(error),
            error_text,
        )?;

        if new_atom_count <= 0 && *error == 0 {
            add_source_error(heap, error_text, "Empty structure")?;
            *error = 98;
        } else if let Some(input_data) = orig_at_data.as_deref_mut()
            && new_atom_count == 0
            && 10 < *error
            && *error < 20
            && input_data.num_atoms > 0
            && merge_all_input_structures != 0
        {
            *error = 0;
            break;
        } else if new_atom_count > 0 {
            if let Some(input_data) = orig_at_data.as_deref_mut() {
                let total_atom_count = new_atom_count
                    .checked_add(i32::from(input_data.num_atoms))
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                let total_stereo_count = new_stereo_count
                    .checked_add(i32::from(input_data.num_stereo0D))
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;

                if total_atom_count >= MAX_ATOMS as i32 {
                    add_source_error(
                        heap,
                        error_text,
                        "Too many atoms [did you forget 'LargeMolecules' switch?]",
                    )?;
                    *error = 70;
                    input_data.num_atoms = -1;
                } else if atom_old.is_null() {
                    input_data.atom = atom_new;
                    atom_new = SourceMutPointer::null();
                    input_data.num_atoms = i16::try_from(new_atom_count)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    input_data.stereo0D = stereo_new;
                    stereo_new = SourceMutPointer::null();
                    input_data.num_stereo0D = i16::try_from(new_stereo_count)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    new_stereo_count = 0;
                    let _ = new_stereo_count;
                } else {
                    input_data.atom = match CreateInchiAtom(heap, total_atom_count) {
                        Ok(atoms) => atoms,
                        Err(SourceHeapError::AllocationFailed) => {
                            add_source_error(heap, error_text, "Out of RAM")?;
                            *error = -1;
                            SourceMutPointer::null()
                        }
                        Err(source_error) => return Err(source_error),
                    };

                    if !input_data.atom.is_null() {
                        let old_atom_count = usize::try_from(input_data.num_atoms)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                        for index in 0..old_atom_count {
                            let atom = heap
                                .slice(atom_old.as_const())?
                                .get(index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            *heap
                                .slice_mut(input_data.atom)?
                                .get_mut(index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)? = atom;
                        }

                        if input_data.num_atoms != 0 {
                            let offset = input_data.num_atoms;
                            let new_atom_count_usize = usize::try_from(new_atom_count)
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                            let atoms = heap.slice_mut(atom_new)?;
                            for atom in atoms
                                .get_mut(..new_atom_count_usize)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                            {
                                let bond_count = usize::try_from(atom.num_bonds)
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                                for neighbor in atom
                                    .neighbor
                                    .get_mut(..bond_count)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                                {
                                    *neighbor = neighbor
                                        .checked_add(offset)
                                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                                }
                            }
                        }

                        let mut atom_old_slot = atom_old;
                        FreeInchi_Atom(heap, Some(&mut atom_old_slot))?;

                        if !atom_new.is_null() {
                            let destination_start = usize::try_from(input_data.num_atoms)
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                            let new_atom_count_usize = usize::try_from(new_atom_count)
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                            for index in 0..new_atom_count_usize {
                                let atom = heap
                                    .slice(atom_new.as_const())?
                                    .get(index)
                                    .cloned()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                *heap
                                    .slice_mut(input_data.atom)?
                                    .get_mut(destination_start + index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)? = atom;
                            }
                        }

                        if new_stereo_count > 0 && !stereo_new.is_null() {
                            input_data.stereo0D =
                                match CreateInchi_Stereo0D(heap, total_stereo_count) {
                                    Ok(stereo) => stereo,
                                    Err(SourceHeapError::AllocationFailed) => {
                                        new_stereo_count = 0;
                                        add_source_error(heap, error_text, "Out of RAM")?;
                                        *error = -1;
                                        SourceMutPointer::null()
                                    }
                                    Err(source_error) => return Err(source_error),
                                };

                            if !input_data.stereo0D.is_null() {
                                let old_stereo_count = usize::try_from(input_data.num_stereo0D)
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                                for index in 0..old_stereo_count {
                                    let stereo = heap
                                        .slice(stereo_old.as_const())?
                                        .get(index)
                                        .cloned()
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    *heap
                                        .slice_mut(input_data.stereo0D)?
                                        .get_mut(index)
                                        .ok_or(SourceHeapError::PointerOutOfBounds)? = stereo;
                                }

                                let atom_offset = input_data.num_atoms;
                                let new_stereo_count_usize = usize::try_from(new_stereo_count)
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                                let new_stereo = heap
                                    .slice_mut(stereo_new)?
                                    .get_mut(..new_stereo_count_usize)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                for stereo in new_stereo {
                                    if stereo.central_atom >= 0 {
                                        stereo.central_atom = stereo
                                            .central_atom
                                            .checked_add(atom_offset)
                                            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                                    }
                                    for neighbor in &mut stereo.neighbor {
                                        *neighbor = neighbor
                                            .checked_add(atom_offset)
                                            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                                    }
                                }

                                let mut stereo_old_slot = stereo_old;
                                FreeInchi_Stereo0D(heap, Some(&mut stereo_old_slot))?;

                                let destination_start = usize::try_from(input_data.num_stereo0D)
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                                for index in 0..new_stereo_count_usize {
                                    let stereo = heap
                                        .slice(stereo_new.as_const())?
                                        .get(index)
                                        .cloned()
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    *heap
                                        .slice_mut(input_data.stereo0D)?
                                        .get_mut(destination_start + index)
                                        .ok_or(SourceHeapError::PointerOutOfBounds)? = stereo;
                                }
                            }
                        } else {
                            new_stereo_count = 0;
                        }

                        input_data.num_atoms = input_data
                            .num_atoms
                            .checked_add(
                                i16::try_from(new_atom_count)
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                            )
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                        input_data.num_stereo0D = input_data
                            .num_stereo0D
                            .checked_add(
                                i16::try_from(new_stereo_count)
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                            )
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    }
                }
            } else {
                accumulated_atom_count = accumulated_atom_count
                    .checked_add(new_atom_count)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
        }

        FreeInchi_Atom(heap, Some(&mut atom_new))?;
        FreeInchi_Stereo0D(heap, Some(&mut stereo_new))?;
        let _ = (new_dimensions, new_bond_count);

        if *error != 0 || merge_all_input_structures == 0 {
            break;
        }
    }

    if !atom_new.is_null() {
        inchi_free(heap, atom_new)?;
    }

    if *error != 0 {
        let input_data = orig_at_data
            .as_deref_mut()
            .ok_or(SourceHeapError::NullPointer)?;
        FreeInchi_Input(heap, input_data)?;
    }
    if *error != 0
        && !(10 < *error && *error < 20)
        && let Some(error_text) = error_text
        && heap
            .slice(error_text.as_const())?
            .first()
            .is_some_and(|byte| *byte == 0)
    {
        add_source_error(heap, Some(error_text), "Unknown error")?;
    }

    Ok(orig_at_data
        .as_deref()
        .map_or(accumulated_atom_count, |input| i32::from(input.num_atoms)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        FLAG_INP_AT_CHIRAL, INCHI_IOS_STRING, INCHI_IOS_TYPE_FILE, INCHI_IOS_TYPE_STRING,
        MAX_SDF_HEADER, MAX_SDF_VALUE, STR_ERR_LEN, SourceFile, SourceMutPointer, inchi_Atom,
        inchi_Stereo0D, tagINCHIStereoParity0D_INCHI_PARITY_ODD,
        tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral, tagInputType_INPUT_INCHI_PLAIN,
    };
    use crate::test_support::allocate_source_fixture;
    use serde_json::{Value, json};
    use std::path::Path;
    use std::process::Command;

    fn string_stream(heap: &mut SourceHeap, text: &[u8]) -> (INCHI_IOSTREAM, SourceMutPointer<i8>) {
        let allocation = allocate_source_fixture(
            heap,
            text.iter().map(|byte| *byte as i8).collect::<Vec<_>>(),
        );
        (
            INCHI_IOSTREAM {
                s: INCHI_IOS_STRING {
                    pStr: allocation,
                    nAllocatedLength: text.len() as i32,
                    nUsedLength: text.len() as i32,
                    nPtr: 0,
                },
                type_: INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            },
            allocation,
        )
    }

    fn c_text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> String {
        let bytes = heap.slice(pointer.as_const()).unwrap();
        let length = bytes.iter().position(|byte| *byte == 0).unwrap();
        String::from_utf8(bytes[..length].iter().map(|byte| *byte as u8).collect()).unwrap()
    }

    fn inline_c_text(bytes: &[i8]) -> String {
        let length = bytes.iter().position(|byte| *byte == 0).unwrap();
        String::from_utf8(bytes[..length].iter().map(|byte| *byte as u8).collect()).unwrap()
    }

    fn focused_aux_case(
        text: &str,
        do_not_add_h: i32,
        distinguish_unknown_undefined_stereo: i32,
    ) -> (
        SourceHeap,
        SourceMutPointer<i8>,
        SourceMutPointer<inchi_Input>,
        InchiInpData,
        i32,
    ) {
        let mut heap = SourceHeap::default();
        let mut bytes = text.bytes().map(|byte| byte as i8).collect::<Vec<_>>();
        bytes.push(0);
        let aux = allocate_source_fixture(&mut heap, bytes);
        let input = allocate_source_fixture(&mut heap, vec![inchi_Input::default()]);
        let mut output = InchiInpData {
            pInp: input,
            bChiral: 99,
            szErrMsg: [b'X' as i8; STR_ERR_LEN as usize],
        };
        let status = Get_inchi_Input_FromAuxInfo(
            &mut heap,
            aux,
            do_not_add_h,
            distinguish_unknown_undefined_stereo,
            Some(&mut output),
        )
        .unwrap();
        (heap, aux, input, output, status)
    }

    fn cleanup_focused_aux_case(
        mut heap: SourceHeap,
        aux: SourceMutPointer<i8>,
        input: SourceMutPointer<inchi_Input>,
    ) {
        let mut state = heap.slice(input.as_const()).unwrap()[0].clone();
        Free_inchi_Input(&mut heap, &mut state).unwrap();
        inchi_free(&mut heap, aux).unwrap();
        inchi_free(&mut heap, input).unwrap();
    }

    #[test]
    fn source_port__ichilnct__get_std_inchi_input_fromauxinfo__line_79() {
        let mut null_heap = SourceHeap::default();
        let null_aux = allocate_source_fixture(
            &mut null_heap,
            b"AuxInfo=1/rA:1C/rB:/rC:;\n\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect(),
        );
        assert_eq!(
            Get_std_inchi_Input_FromAuxInfo(&mut null_heap, null_aux, 0, None),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        inchi_free(&mut null_heap, null_aux).unwrap();

        for (do_not_add_h, expected_ordinary_h) in [(0, -1_i8), (1, 0_i8)] {
            let mut heap = SourceHeap::default();
            let aux = allocate_source_fixture(
                &mut heap,
                b"AuxInfo=1/rA:1C/rB:/rC:;\n\0"
                    .iter()
                    .map(|byte| *byte as i8)
                    .collect(),
            );
            let input_pointer = allocate_source_fixture(&mut heap, vec![inchi_Input::default()]);
            let mut output = InchiInpData {
                pInp: input_pointer,
                ..InchiInpData::default()
            };

            assert_eq!(
                Get_std_inchi_Input_FromAuxInfo(&mut heap, aux, do_not_add_h, Some(&mut output),),
                Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
            );
            assert_eq!(output.bChiral, 0);
            assert_eq!(inline_c_text(&output.szErrMsg), "");
            let mut parsed = heap.slice(input_pointer.as_const()).unwrap()[0].clone();
            assert_eq!((parsed.num_atoms, parsed.num_stereo0D), (1, 0));
            assert_eq!(
                heap.slice(parsed.atom.as_const()).unwrap()[0].num_iso_H[0],
                expected_ordinary_h
            );

            Free_inchi_Input(&mut heap, &mut parsed).unwrap();
            inchi_free(&mut heap, aux).unwrap();
            inchi_free(&mut heap, input_pointer).unwrap();
        }
    }

    #[test]
    fn source_port__ichilnct__get_inchi_input_fromauxinfo__line_89() {
        let mut heap = SourceHeap::default();
        let aux = allocate_source_fixture(
            &mut heap,
            b"AuxInfo=1/rA:1C/rB:/rC:;\n\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect(),
        );
        assert_eq!(
            Get_inchi_Input_FromAuxInfo(&mut heap, aux, 0, 0, None),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        let mut missing_input = InchiInpData {
            bChiral: 7,
            szErrMsg: [9; STR_ERR_LEN as usize],
            ..InchiInpData::default()
        };
        let missing_before = missing_input.clone();
        assert_eq!(
            Get_inchi_Input_FromAuxInfo(&mut heap, aux, 0, 0, Some(&mut missing_input)),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        assert_eq!(missing_input, missing_before);
        inchi_free(&mut heap, aux).unwrap();

        let options = allocate_source_fixture(&mut heap, vec![b'-' as i8, b'S' as i8, 0]);
        let old_atom = allocate_source_fixture(&mut heap, vec![inchi_Atom::default()]);
        let old_atom_alias = old_atom.as_const();
        let old_stereo = allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default()]);
        let old_stereo_alias = old_stereo.as_const();
        let input_pointer = allocate_source_fixture(
            &mut heap,
            vec![inchi_Input {
                atom: old_atom,
                stereo0D: old_stereo,
                szOptions: options,
                num_atoms: 1,
                num_stereo0D: 1,
            }],
        );
        let labeled_aux = allocate_source_fixture(
            &mut heap,
            b"Structure: 42. LABEL =value\nAuxInfo=1/0/N:2/rA:2cCO/rB:s1;/rC:;;\n\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect(),
        );
        let mut output = InchiInpData {
            pInp: input_pointer,
            bChiral: 99,
            szErrMsg: [b'X' as i8; STR_ERR_LEN as usize],
        };
        assert_eq!(
            Get_inchi_Input_FromAuxInfo(&mut heap, labeled_aux, 0, 0, Some(&mut output)),
            Ok(tagRetValGetINCHI_inchi_Ret_OKAY)
        );
        assert_eq!(output.pInp, input_pointer);
        assert_eq!(output.bChiral, FLAG_INP_AT_CHIRAL as i32);
        assert_eq!(inline_c_text(&output.szErrMsg), "");
        let parsed = heap.slice(input_pointer.as_const()).unwrap()[0].clone();
        assert_eq!(parsed.szOptions, options);
        assert_eq!((parsed.num_atoms, parsed.num_stereo0D), (2, 0));
        assert_eq!(
            heap.slice(parsed.atom.as_const()).unwrap()[0].neighbor[0],
            1
        );
        assert_eq!(heap.slice(old_atom_alias).unwrap().len(), 1);
        assert_eq!(heap.slice(old_stereo_alias).unwrap().len(), 1);
        let mut parsed_cleanup = parsed;
        Free_inchi_Input(&mut heap, &mut parsed_cleanup).unwrap();
        heap.slice_mut(input_pointer).unwrap()[0] = parsed_cleanup;
        inchi_free(&mut heap, old_atom).unwrap();
        inchi_free(&mut heap, old_stereo).unwrap();
        inchi_free(&mut heap, options).unwrap();
        inchi_free(&mut heap, labeled_aux).unwrap();
        inchi_free(&mut heap, input_pointer).unwrap();

        for (text, expected_return, expected_error, options_survive) in [
            (
                b"AuxInfo=1//\n\0".as_slice(),
                tagRetValGetINCHI_inchi_Ret_WARNING,
                "Empty structure",
                false,
            ),
            (b"\0".as_slice(), tagRetValGetINCHI_inchi_Ret_EOF, "", false),
            (
                b"AuxInfo=1/rA:2CO/rB:x/rC:;;\n\0".as_slice(),
                tagRetValGetINCHI_inchi_Ret_ERROR,
                "Wrong bonds data",
                false,
            ),
        ] {
            let mut case_heap = SourceHeap::default();
            let case_options =
                allocate_source_fixture(&mut case_heap, vec![b'-' as i8, b'X' as i8, 0]);
            let case_input = allocate_source_fixture(
                &mut case_heap,
                vec![inchi_Input {
                    szOptions: case_options,
                    ..inchi_Input::default()
                }],
            );
            let case_aux = allocate_source_fixture(
                &mut case_heap,
                text.iter().map(|byte| *byte as i8).collect(),
            );
            let mut case_output = InchiInpData {
                pInp: case_input,
                bChiral: 8,
                szErrMsg: [b'Y' as i8; STR_ERR_LEN as usize],
            };
            assert_eq!(
                Get_inchi_Input_FromAuxInfo(&mut case_heap, case_aux, 0, 1, Some(&mut case_output),),
                Ok(expected_return),
                "{}",
                String::from_utf8_lossy(text)
            );
            assert_eq!(inline_c_text(&case_output.szErrMsg), expected_error);
            assert_eq!(case_output.bChiral, 0);
            let case_state = &case_heap.slice(case_input.as_const()).unwrap()[0];
            assert_eq!(case_state.szOptions == case_options, options_survive);
            assert_eq!((case_state.num_atoms, case_state.num_stereo0D), (0, 0));
            inchi_free(&mut case_heap, case_options).unwrap();
            inchi_free(&mut case_heap, case_aux).unwrap();
            inchi_free(&mut case_heap, case_input).unwrap();
        }

        let mut warning_heap = SourceHeap::default();
        let warning_input =
            allocate_source_fixture(&mut warning_heap, vec![inchi_Input::default()]);
        let warning_aux = allocate_source_fixture(
            &mut warning_heap,
            b"AuxInfo=1/rA:2CO/rB:a1;/rC:;;\n\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect(),
        );
        let mut warning_output = InchiInpData {
            pInp: warning_input,
            ..InchiInpData::default()
        };
        assert_eq!(
            Get_inchi_Input_FromAuxInfo(
                &mut warning_heap,
                warning_aux,
                0,
                0,
                Some(&mut warning_output),
            ),
            Ok(tagRetValGetINCHI_inchi_Ret_ERROR)
        );
        assert_eq!(
            inline_c_text(&warning_output.szErrMsg),
            "Atom has 1 or more than 3 aromatic bonds"
        );
        assert_eq!(warning_output.bChiral, 0);
        let mut warning_state = warning_heap.slice(warning_input.as_const()).unwrap()[0].clone();
        assert_eq!(
            (warning_state.num_atoms, warning_state.num_stereo0D),
            (0, 0)
        );
        Free_inchi_Input(&mut warning_heap, &mut warning_state).unwrap();
        inchi_free(&mut warning_heap, warning_aux).unwrap();
        inchi_free(&mut warning_heap, warning_input).unwrap();

        let mut fatal_heap = SourceHeap::default();
        let fatal_input = allocate_source_fixture(&mut fatal_heap, vec![inchi_Input::default()]);
        let fatal_aux = allocate_source_fixture(
            &mut fatal_heap,
            b"AuxInfo=1/rA:1C/rB:/rC:;\n\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect(),
        );
        let mut fatal_output = InchiInpData {
            pInp: fatal_input,
            ..InchiInpData::default()
        };
        fatal_heap.fail_after_allocations(0);
        assert_eq!(
            Get_inchi_Input_FromAuxInfo(&mut fatal_heap, fatal_aux, 0, 0, Some(&mut fatal_output),),
            Ok(tagRetValGetINCHI_inchi_Ret_FATAL)
        );
        assert_eq!(inline_c_text(&fatal_output.szErrMsg), "Out of RAM");
        assert_eq!(fatal_output.bChiral, 0);
        assert_eq!(
            fatal_heap.slice(fatal_input.as_const()).unwrap()[0],
            inchi_Input::default()
        );
        inchi_free(&mut fatal_heap, fatal_aux).unwrap();
        inchi_free(&mut fatal_heap, fatal_input).unwrap();

        let (no_h_heap, no_h_aux, no_h_input, _, no_h_status) =
            focused_aux_case("AuxInfo=1/rA:1C/rB:/rC:;\n", 1, 0);
        assert_eq!(no_h_status, tagRetValGetINCHI_inchi_Ret_OKAY);
        let no_h_state = &no_h_heap.slice(no_h_input.as_const()).unwrap()[0];
        assert_eq!(
            no_h_heap.slice(no_h_state.atom.as_const()).unwrap()[0].num_iso_H[0],
            0
        );
        cleanup_focused_aux_case(no_h_heap, no_h_aux, no_h_input);

        let merge_text = "AuxInfo=1/rA:1C/rB:/rC:;\n\
AuxInfo=1/rA:2NO0o/rB:s1;/rC:;;\n";
        let (merge_heap, merge_aux, merge_input, _, merge_status) =
            focused_aux_case(merge_text, 0, 0);
        assert_eq!(merge_status, tagRetValGetINCHI_inchi_Ret_OKAY);
        let merge_state = &merge_heap.slice(merge_input.as_const()).unwrap()[0];
        assert_eq!((merge_state.num_atoms, merge_state.num_stereo0D), (3, 1));
        let merge_atoms = merge_heap.slice(merge_state.atom.as_const()).unwrap();
        assert_eq!(
            (merge_atoms[1].neighbor[0], merge_atoms[2].neighbor[0]),
            (2, 1)
        );
        assert_eq!(
            merge_heap.slice(merge_state.stereo0D.as_const()).unwrap()[0].central_atom,
            2
        );
        cleanup_focused_aux_case(merge_heap, merge_aux, merge_input);

        let unknown_text = "AuxInfo=1/rA:4CCCClu/rB:;;s1s2s3;/rC:;;;;\n";
        for distinguish in [0, 1] {
            let (unknown_heap, unknown_aux, unknown_input, _, unknown_status) =
                focused_aux_case(unknown_text, 0, distinguish);
            assert_eq!(unknown_status, tagRetValGetINCHI_inchi_Ret_OKAY);
            let unknown_state = &unknown_heap.slice(unknown_input.as_const()).unwrap()[0];
            assert_eq!(unknown_state.num_stereo0D, 1);
            assert_eq!(
                unknown_heap
                    .slice(unknown_state.stereo0D.as_const())
                    .unwrap()[0]
                    .parity,
                3
            );
            cleanup_focused_aux_case(unknown_heap, unknown_aux, unknown_input);
        }

        let max_text = generated_aux_input_text(1);
        let (max_heap, max_aux, max_input, max_output, max_status) =
            focused_aux_case(&max_text, 0, 0);
        assert_eq!(max_status, tagRetValGetINCHI_inchi_Ret_ERROR);
        assert_eq!(
            inline_c_text(&max_output.szErrMsg),
            "Too many atoms [did you forget 'LargeMolecules' switch?]"
        );
        assert_eq!(
            max_heap.slice(max_input.as_const()).unwrap()[0],
            inchi_Input::default()
        );
        cleanup_focused_aux_case(max_heap, max_aux, max_input);
    }

    #[test]
    fn source_port__ichilnct__inchitoinchi_input__line_204() {
        let mut heap = SourceHeap::default();
        let error_text = allocate_source_fixture(&mut heap, vec![b'X' as i8, 0, 0]);
        let mut id = 91_i64;
        let mut error = 0_i32;
        let mut output = inchi_Input::default();
        assert_eq!(
            InchiToInchi_Input(
                &mut heap,
                None,
                Some(&mut output),
                0,
                0,
                0,
                tagInputType_INPUT_INCHI_PLAIN,
                None,
                None,
                Some(&mut id),
                None,
                Some(&mut error),
                Some(error_text),
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(id, 0);
        assert_eq!(c_text(&heap, error_text), "");
        inchi_free(&mut heap, error_text).unwrap();

        let count_text = b"AuxInfo=1/rA:2CO/rB:s1;/rC:;;\n";
        let (mut count_stream, count_allocation) = string_stream(&mut heap, count_text);
        let mut count_error = 0;
        assert_eq!(
            InchiToInchi_Input(
                &mut heap,
                Some(&mut count_stream),
                None,
                0,
                0,
                0,
                tagInputType_INPUT_INCHI_PLAIN,
                None,
                None,
                None,
                None,
                Some(&mut count_error),
                None,
            ),
            Ok(2)
        );
        assert_eq!(count_error, 0);
        inchi_free(&mut heap, count_allocation).unwrap();

        let labeled = b"Structure: 42. LABEL =value\n\
AuxInfo=1/0/N:2/rA:2cCO/rB:s1;/rC:0,0,0;1,1,1;\n";
        let (mut labeled_stream, labeled_allocation) = string_stream(&mut heap, labeled);
        let label = allocate_source_fixture(&mut heap, vec![0_i8; MAX_SDF_HEADER as usize]);
        let value = allocate_source_fixture(&mut heap, vec![0_i8; MAX_SDF_VALUE as usize]);
        let labeled_error_text =
            allocate_source_fixture(&mut heap, vec![b'X' as i8; STR_ERR_LEN as usize]);
        heap.slice_mut(labeled_error_text).unwrap()[1] = 0;
        let mut labeled_output = inchi_Input::default();
        let mut labeled_id = -1;
        let mut flags = 0;
        let mut labeled_error = 0;
        assert_eq!(
            InchiToInchi_Input(
                &mut heap,
                Some(&mut labeled_stream),
                Some(&mut labeled_output),
                0,
                0,
                0,
                tagInputType_INPUT_INCHI_PLAIN,
                Some(label),
                Some(value),
                Some(&mut labeled_id),
                Some(&mut flags),
                Some(&mut labeled_error),
                Some(labeled_error_text),
            ),
            Ok(2)
        );
        assert_eq!(labeled_error, 0);
        assert_eq!(c_text(&heap, labeled_error_text), "");
        assert_eq!((labeled_id, flags), (42, u64::from(FLAG_INP_AT_CHIRAL)));
        assert_eq!(
            (c_text(&heap, label), c_text(&heap, value)),
            ("LABEL".into(), "value".into())
        );
        assert_eq!(labeled_output.num_atoms, 2);
        assert!(!labeled_output.atom.is_null());
        assert_eq!(
            heap.slice(labeled_output.atom.as_const()).unwrap()[0].neighbor[0],
            1
        );
        FreeInchi_Input(&mut heap, &mut labeled_output).unwrap();
        for allocation in [labeled_allocation, label, value, labeled_error_text] {
            inchi_free(&mut heap, allocation).unwrap();
        }

        let merged = b"Structure: 1. FIRST =one\n\
AuxInfo=1/rA:1C/rB:/rC:;\n\
Structure: 2. SECOND =two\n\
AuxInfo=1/rA:2NO/rB:s1;/rC:;;\n\
Structure: 3. THIRD =three\n\
AuxInfo=1/rA:1Coe/rB:/rC:;\n";
        let (mut merge_stream, merge_allocation) = string_stream(&mut heap, merged);
        let merge_label = allocate_source_fixture(&mut heap, vec![0_i8; MAX_SDF_HEADER as usize]);
        let merge_value = allocate_source_fixture(&mut heap, vec![0_i8; MAX_SDF_VALUE as usize]);
        let merge_error_text = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
        let mut merged_output = inchi_Input::default();
        let mut merge_id = -1;
        let mut merge_error = 0;
        assert_eq!(
            InchiToInchi_Input(
                &mut heap,
                Some(&mut merge_stream),
                Some(&mut merged_output),
                1,
                0,
                0,
                tagInputType_INPUT_INCHI_PLAIN,
                Some(merge_label),
                Some(merge_value),
                Some(&mut merge_id),
                None,
                Some(&mut merge_error),
                Some(merge_error_text),
            ),
            Ok(4)
        );
        assert_eq!(merge_error, 0, "merge EOF must be reset to success");
        assert_eq!(merged_output.num_atoms, 4);
        assert_eq!(merged_output.num_stereo0D, 1);
        assert_eq!(
            (
                merge_id,
                c_text(&heap, merge_label),
                c_text(&heap, merge_value)
            ),
            (3, "THIRD".into(), "three".into())
        );
        let merged_atoms = heap.slice(merged_output.atom.as_const()).unwrap();
        assert_eq!(
            (merged_atoms[1].neighbor[0], merged_atoms[2].neighbor[0]),
            (2, 1)
        );
        let merged_stereo = &heap.slice(merged_output.stereo0D.as_const()).unwrap()[0];
        assert_eq!(merged_stereo.central_atom, 3);
        assert_eq!(merged_stereo.neighbor, [3; 4]);
        FreeInchi_Input(&mut heap, &mut merged_output).unwrap();
        for allocation in [merge_allocation, merge_label, merge_value, merge_error_text] {
            inchi_free(&mut heap, allocation).unwrap();
        }

        let (mut empty_stream, empty_allocation) = string_stream(&mut heap, b"AuxInfo=1//\n");
        let empty_error_text = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
        let mut empty_output = inchi_Input::default();
        let mut empty_error = 0;
        assert_eq!(
            InchiToInchi_Input(
                &mut heap,
                Some(&mut empty_stream),
                Some(&mut empty_output),
                0,
                0,
                0,
                tagInputType_INPUT_INCHI_PLAIN,
                None,
                None,
                None,
                None,
                Some(&mut empty_error),
                Some(empty_error_text),
            ),
            Ok(0)
        );
        assert_eq!(empty_error, 98);
        assert_eq!(c_text(&heap, empty_error_text), "Empty structure");
        inchi_free(&mut heap, empty_allocation).unwrap();
        inchi_free(&mut heap, empty_error_text).unwrap();

        let (mut boundary_stream, boundary_allocation) =
            string_stream(&mut heap, b"AuxInfo=1/rA:1C/rB:/rC:;\n");
        let boundary_atom = allocate_source_fixture(&mut heap, vec![inchi_Atom::default()]);
        let boundary_alias = boundary_atom.as_const();
        let boundary_error_text =
            allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
        let mut boundary_output = inchi_Input {
            atom: boundary_atom,
            num_atoms: MAX_ATOMS as i16 - 1,
            ..inchi_Input::default()
        };
        let mut boundary_error = 0;
        assert_eq!(
            InchiToInchi_Input(
                &mut heap,
                Some(&mut boundary_stream),
                Some(&mut boundary_output),
                0,
                0,
                0,
                tagInputType_INPUT_INCHI_PLAIN,
                None,
                None,
                None,
                None,
                Some(&mut boundary_error),
                Some(boundary_error_text),
            ),
            Ok(0)
        );
        assert_eq!(boundary_error, 70);
        assert_eq!(
            c_text(&heap, boundary_error_text),
            "Too many atoms [did you forget 'LargeMolecules' switch?]"
        );
        assert_eq!(
            heap.slice(boundary_alias),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(boundary_output, inchi_Input::default());
        inchi_free(&mut heap, boundary_allocation).unwrap();
        inchi_free(&mut heap, boundary_error_text).unwrap();

        let (mut malformed_stream, malformed_allocation) =
            string_stream(&mut heap, b"AuxInfo=1/rA:2CO/rB:x/rC:;;\n");
        let malformed_atom = allocate_source_fixture(&mut heap, vec![inchi_Atom::default()]);
        let malformed_alias = malformed_atom.as_const();
        let malformed_error_text =
            allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
        let mut malformed_output = inchi_Input {
            atom: malformed_atom,
            num_atoms: 1,
            ..inchi_Input::default()
        };
        let mut malformed_error = 0;
        assert_eq!(
            InchiToInchi_Input(
                &mut heap,
                Some(&mut malformed_stream),
                Some(&mut malformed_output),
                0,
                0,
                0,
                tagInputType_INPUT_INCHI_PLAIN,
                None,
                None,
                None,
                None,
                Some(&mut malformed_error),
                Some(malformed_error_text),
            ),
            Ok(0)
        );
        assert_eq!(malformed_error, 40);
        assert_eq!(c_text(&heap, malformed_error_text), "Wrong bonds data");
        assert_eq!(
            heap.slice(malformed_alias),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(malformed_output, inchi_Input::default());
        inchi_free(&mut heap, malformed_allocation).unwrap();
        inchi_free(&mut heap, malformed_error_text).unwrap();

        let (mut unknown_stream, unknown_allocation) =
            string_stream(&mut heap, b"AuxInfo=1/rA:1C/rB:/rC:;\n");
        let unknown_error_text =
            allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);
        let mut unknown_output = inchi_Input::default();
        let mut unknown_error = 77;
        assert_eq!(
            InchiToInchi_Input(
                &mut heap,
                Some(&mut unknown_stream),
                Some(&mut unknown_output),
                0,
                0,
                0,
                tagInputType_INPUT_INCHI_PLAIN,
                None,
                None,
                None,
                None,
                Some(&mut unknown_error),
                Some(unknown_error_text),
            ),
            Ok(0)
        );
        assert_eq!(unknown_error, 77);
        assert_eq!(c_text(&heap, unknown_error_text), "Unknown error");
        assert_eq!(unknown_output, inchi_Input::default());
        inchi_free(&mut heap, unknown_allocation).unwrap();
        inchi_free(&mut heap, unknown_error_text).unwrap();

        for (text, allocations_before_failure) in [
            (
                b"AuxInfo=1/rA:1C/rB:/rC:;\nAuxInfo=1/rA:1N/rB:/rC:;\n".as_slice(),
                6_u64,
            ),
            (
                b"AuxInfo=1/rA:1Coe/rB:/rC:;\nAuxInfo=1/rA:1C0o/rB:/rC:;\n".as_slice(),
                7_u64,
            ),
        ] {
            let mut failure_heap = SourceHeap::default();
            let (mut failure_stream, failure_allocation) = string_stream(&mut failure_heap, text);
            let failure_error_text =
                allocate_source_fixture(&mut failure_heap, vec![0_i8; STR_ERR_LEN as usize]);
            let mut failure_output = inchi_Input::default();
            let mut failure_error = 0;
            failure_heap.fail_after_allocations(allocations_before_failure);
            assert_eq!(
                InchiToInchi_Input(
                    &mut failure_heap,
                    Some(&mut failure_stream),
                    Some(&mut failure_output),
                    1,
                    0,
                    0,
                    tagInputType_INPUT_INCHI_PLAIN,
                    None,
                    None,
                    None,
                    None,
                    Some(&mut failure_error),
                    Some(failure_error_text),
                ),
                Ok(0),
                "merge allocation failure after {allocations_before_failure} successful allocations"
            );
            assert_eq!(failure_error, -1);
            assert_eq!(c_text(&failure_heap, failure_error_text), "Out of RAM");
            assert_eq!(failure_output, inchi_Input::default());
            inchi_free(&mut failure_heap, failure_allocation).unwrap();
            inchi_free(&mut failure_heap, failure_error_text).unwrap();
        }
    }

    fn inchi_input_atom_json(atom: &inchi_Atom) -> Value {
        json!({
            "x_bits": atom.x.to_bits(),
            "y_bits": atom.y.to_bits(),
            "z_bits": atom.z.to_bits(),
            "neighbor": atom.neighbor,
            "bond_type": atom.bond_type,
            "bond_stereo": atom.bond_stereo,
            "elname_bytes": atom.elname,
            "num_bonds": atom.num_bonds,
            "num_iso_h": atom.num_iso_H,
            "isotopic_mass": atom.isotopic_mass,
            "radical": atom.radical,
            "charge": atom.charge,
        })
    }

    fn inchi_input_stereo_json(stereo: &inchi_Stereo0D) -> Value {
        json!({
            "central_atom": stereo.central_atom,
            "neighbor": stereo.neighbor,
            "type": stereo.type_,
            "parity": stereo.parity,
        })
    }

    fn official_c_inchi_to_input_rust_record(official: &Value) -> Value {
        let case_id = official["case_id"].as_str().expect("case_id must be text");
        let input = &official["input"];
        let text = input["text"].as_str().expect("input text must be text");
        let file_mode = input["mode"].as_str() == Some("file");
        let input_type = input["input_type"]
            .as_i64()
            .expect("input_type must be integer") as INPUT_TYPE;
        let merge_all = input["merge_all"]
            .as_i64()
            .expect("merge_all must be integer") as i32;
        let do_not_add_h = input["do_not_add_h"]
            .as_i64()
            .expect("do_not_add_h must be integer") as i32;
        let ab_parity_unknown = input["ab_parity_unknown"]
            .as_i64()
            .expect("ab_parity_unknown must be integer") as i32;
        let initial_atom_count = input["initial_atom_count"]
            .as_i64()
            .expect("initial_atom_count must be integer") as i16;
        let initial_atom_storage = input["initial_atom_storage"]
            .as_u64()
            .expect("initial_atom_storage must be integer")
            as usize;
        let initial_stereo_count = input["initial_stereo_count"]
            .as_i64()
            .expect("initial_stereo_count must be integer")
            as i16;
        let initial_stereo_storage = input["initial_stereo_storage"]
            .as_u64()
            .expect("initial_stereo_storage must be integer")
            as usize;
        let allocation_failure_ordinal = input["allocation_failure_ordinal"]
            .as_u64()
            .expect("allocation_failure_ordinal must be integer");
        let initial_error = input["initial_error"]
            .as_i64()
            .expect("initial_error must be integer") as i32;
        let omit_output = input["omit_output"]
            .as_i64()
            .expect("omit_output must be integer")
            != 0;
        let omit_optional_outputs = input["omit_optional_outputs"]
            .as_i64()
            .expect("omit_optional_outputs must be integer")
            != 0;

        let mut heap = SourceHeap::default();
        let mut string_allocation = SourceMutPointer::null();
        let mut file_allocation = SourceMutPointer::null();
        let mut stream = if file_mode {
            file_allocation = allocate_source_fixture(
                &mut heap,
                vec![SourceFile {
                    bytes: text.as_bytes().to_vec(),
                    ..SourceFile::default()
                }],
            );
            INCHI_IOSTREAM {
                f: file_allocation,
                type_: INCHI_IOS_TYPE_FILE as i32,
                ..INCHI_IOSTREAM::default()
            }
        } else {
            let mut bytes = text.bytes().map(|byte| byte as i8).collect::<Vec<_>>();
            bytes.push(0);
            string_allocation = allocate_source_fixture(&mut heap, bytes);
            INCHI_IOSTREAM {
                s: INCHI_IOS_STRING {
                    pStr: string_allocation,
                    nAllocatedLength: text.len() as i32 + 1,
                    nUsedLength: text.len() as i32,
                    nPtr: 0,
                },
                type_: INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            }
        };
        let label = allocate_source_fixture(&mut heap, vec![0_i8; MAX_SDF_HEADER as usize]);
        let value = allocate_source_fixture(&mut heap, vec![0_i8; MAX_SDF_VALUE as usize]);
        let error_text = allocate_source_fixture(&mut heap, vec![0_i8; STR_ERR_LEN as usize]);

        let mut initial_atom = inchi_Atom::default();
        initial_atom.x = -0.0;
        initial_atom.y = 2.5;
        initial_atom.z = -3.5;
        initial_atom.elname[0] = b'F' as i8;
        initial_atom.num_iso_H[0] = 4;
        initial_atom.charge = -1;
        let original_atoms = if initial_atom_storage > 0 {
            allocate_source_fixture(&mut heap, vec![initial_atom; initial_atom_storage])
        } else {
            SourceMutPointer::null()
        };

        let mut initial_stereo = inchi_Stereo0D::default();
        initial_stereo.type_ = tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8;
        initial_stereo.parity = tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8;
        let original_stereo = if initial_stereo_storage > 0 {
            allocate_source_fixture(&mut heap, vec![initial_stereo; initial_stereo_storage])
        } else {
            SourceMutPointer::null()
        };
        let original_options = if omit_output {
            SourceMutPointer::null()
        } else {
            allocate_source_fixture(&mut heap, vec![b'-' as i8, b'X' as i8, 0, 0x5a_i8])
        };
        let mut output = inchi_Input {
            atom: original_atoms,
            stereo0D: original_stereo,
            szOptions: original_options,
            num_atoms: initial_atom_count,
            num_stereo0D: initial_stereo_count,
        };
        let mut sdf_id = -1_i64;
        let mut atom_flags = 0_u64;
        let mut error = initial_error;

        if allocation_failure_ordinal == 0 {
            heap.trace_source_allocations();
        } else {
            heap.fail_after_allocations(allocation_failure_ordinal - 1);
        }
        let status = InchiToInchi_Input(
            &mut heap,
            Some(&mut stream),
            (!omit_output).then_some(&mut output),
            merge_all,
            do_not_add_h,
            ab_parity_unknown,
            input_type,
            (!omit_optional_outputs).then_some(label),
            (!omit_optional_outputs).then_some(value),
            (!omit_optional_outputs).then_some(&mut sdf_id),
            (!omit_optional_outputs).then_some(&mut atom_flags),
            Some(&mut error),
            (!omit_optional_outputs).then_some(error_text),
        )
        .unwrap_or_else(|source_error| {
            panic!("Rust source model returned {source_error:?} for {case_id}")
        });
        let allocation_calls = heap.source_allocation_calls();
        let input_position = if file_mode {
            heap.slice(file_allocation.as_const()).unwrap()[0].position as i64
        } else {
            i64::from(stream.s.nPtr)
        };
        let original_atom_freed = !original_atoms.is_null()
            && heap.slice(original_atoms.as_const()) == Err(SourceHeapError::MissingAllocation);
        let original_stereo_freed = !original_stereo.is_null()
            && heap.slice(original_stereo.as_const()) == Err(SourceHeapError::MissingAllocation);

        let (atoms, bond_fields, stereo0d) = if omit_output {
            (Vec::new(), Vec::new(), Vec::new())
        } else {
            let atom_count = usize::try_from(output.num_atoms.max(0)).unwrap();
            let atoms = if atom_count == 0 || output.atom.is_null() {
                Vec::new()
            } else {
                heap.slice(output.atom.as_const()).unwrap()[..atom_count]
                    .iter()
                    .map(inchi_input_atom_json)
                    .collect::<Vec<_>>()
            };
            let bond_fields = if atom_count == 0 || output.atom.is_null() {
                Vec::new()
            } else {
                heap.slice(output.atom.as_const()).unwrap()[..atom_count]
                    .iter()
                    .enumerate()
                    .flat_map(|(atom_index, atom)| {
                        (0..usize::try_from(atom.num_bonds).unwrap()).map(move |slot| {
                            json!({
                                "atom_index": atom_index,
                                "slot": slot,
                                "neighbor": atom.neighbor[slot],
                                "bond_type": atom.bond_type[slot],
                                "bond_stereo": atom.bond_stereo[slot],
                            })
                        })
                    })
                    .collect::<Vec<_>>()
            };
            let stereo_count = usize::try_from(output.num_stereo0D.max(0)).unwrap();
            let stereo = if stereo_count == 0 || output.stereo0D.is_null() {
                Vec::new()
            } else {
                heap.slice(output.stereo0D.as_const()).unwrap()[..stereo_count]
                    .iter()
                    .map(inchi_input_stereo_json)
                    .collect::<Vec<_>>()
            };
            (atoms, bond_fields, stereo)
        };
        let pointer_state = |pointer: SourceMutPointer<inchi_Atom>,
                             original: SourceMutPointer<inchi_Atom>| {
            if pointer.is_null() {
                "null"
            } else if !original.is_null() && pointer == original {
                "reused"
            } else {
                "allocated"
            }
        };
        let stereo_pointer_state = if omit_output {
            Value::Null
        } else if output.stereo0D.is_null() {
            json!("null")
        } else if !original_stereo.is_null() && output.stereo0D == original_stereo {
            json!("reused")
        } else {
            json!("allocated")
        };
        let options_pointer_state = if omit_output {
            Value::Null
        } else if output.szOptions.is_null() {
            json!("null")
        } else if output.szOptions == original_options {
            json!("reused")
        } else {
            json!("allocated")
        };
        let record = json!({
            "schema_version": "cosmolkit-inchi-official-c-v1",
            "oracle": {
                "project": "IUPAC InChI",
                "tag": "v1.07.5",
                "commit": "11a87982bb518f57ac013f0b258c283655e1ea1d",
                "tree_sha256": "4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd",
                "api_version": "1.07.5",
            },
            "case_id": case_id,
            "operation": "inchi_to_inchi_input",
            "input": official["input"].clone(),
            "output": {
                "status": status,
                "error_code": error,
                "error_text": (!omit_optional_outputs).then(|| c_text(&heap, error_text)),
                "input_position": input_position,
                "sdf_label": (!omit_optional_outputs).then(|| c_text(&heap, label)),
                "sdf_value": (!omit_optional_outputs).then(|| c_text(&heap, value)),
                "sdf_id": (!omit_optional_outputs).then_some(sdf_id),
                "atom_flags": (!omit_optional_outputs).then_some(atom_flags),
                "allocation_calls": allocation_calls,
                "num_atoms": (!omit_output).then_some(output.num_atoms),
                "num_stereo0d": (!omit_output).then_some(output.num_stereo0D),
                "atoms": atoms,
                "bond_fields": bond_fields,
                "stereo0d": stereo0d,
                "atom_pointer_state": (!omit_output)
                    .then(|| pointer_state(output.atom, original_atoms)),
                "stereo_pointer_state": stereo_pointer_state,
                "options_pointer_state": options_pointer_state,
                "original_atom_freed": original_atom_freed,
                "original_stereo_freed": original_stereo_freed,
            },
        });

        if !omit_output {
            FreeInchi_Input(&mut heap, &mut output).unwrap();
        }
        for pointer in [original_atoms, output.atom] {
            if !pointer.is_null() && heap.slice(pointer.as_const()).is_ok() {
                inchi_free(&mut heap, pointer).unwrap();
            }
        }
        for pointer in [original_stereo, output.stereo0D] {
            if !pointer.is_null() && heap.slice(pointer.as_const()).is_ok() {
                inchi_free(&mut heap, pointer).unwrap();
            }
        }
        if !original_options.is_null() {
            inchi_free(&mut heap, original_options).unwrap();
        }
        for pointer in [label, value, error_text] {
            inchi_free(&mut heap, pointer).unwrap();
        }
        if file_mode {
            inchi_free(&mut heap, file_allocation).unwrap();
        } else {
            inchi_free(&mut heap, string_allocation).unwrap();
        }
        record
    }

    #[test]
    fn official_c_oracle__ichilnct__inchitoinchi_input__line_204() {
        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--inchi-to-inchi-input-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let official_records = String::from_utf8(oracle.stdout)
            .expect("official C oracle output must be UTF-8")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("oracle record must be JSON"))
            .collect::<Vec<_>>();
        let mut expected_case_ids = [
            "single-labeled-string",
            "single-labeled-file",
            "single-stops-before-second",
            "merge-three-eof",
            "merge-preexisting-remap",
            "count-only",
            "empty-structure",
            "malformed-cleans-old",
            "max-atoms",
            "seeded-error-unknown-text",
            "xml-single",
            "standalone-eof",
            "warning-text-success",
            "xml-fatal-cleanup",
            "merge-double-bond-stereo",
            "non-null-zero-count-old",
            "null-optional-eof",
            "null-optional-malformed",
        ]
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
        expected_case_ids.extend((1..=16).map(|ordinal| format!("allocation-atom-{ordinal}")));
        expected_case_ids.extend((1..=16).map(|ordinal| format!("allocation-stereo-{ordinal}")));
        assert_eq!(
            official_records
                .iter()
                .map(|record| record["case_id"].as_str().unwrap().to_owned())
                .collect::<Vec<_>>(),
            expected_case_ids
        );

        for official in &official_records {
            let case_id = official["case_id"].as_str().unwrap();
            let rust = official_c_inchi_to_input_rust_record(official);
            assert_eq!(official, &rust, "exact field mismatch for {case_id}");
        }
    }

    fn generated_aux_input_text(kind: i64) -> String {
        assert_eq!(kind, 1, "unknown generated AuxInfo fixture");
        let atom_count = MAX_ATOMS as usize;
        let mut text = String::from("AuxInfo=1/rA:32766");
        text.push_str(&"C".repeat(atom_count));
        text.push_str("/rB:");
        text.push_str(&";".repeat(atom_count - 1));
        text.push_str("/rC:");
        text.push_str(&";".repeat(atom_count));
        text.push('\n');
        text
    }

    fn official_c_aux_input_rust_record(official: &Value, standard: bool) -> Value {
        let case_id = official["case_id"].as_str().expect("case_id must be text");
        let input = &official["input"];
        let generated_text = input["generated_text"]
            .as_i64()
            .expect("generated_text must be integer");
        let text = if generated_text == 0 {
            input["text"]
                .as_str()
                .expect("literal AuxInfo must be text")
                .to_owned()
        } else {
            generated_aux_input_text(generated_text)
        };
        let do_not_add_h = input["do_not_add_h"]
            .as_i64()
            .expect("do_not_add_h must be integer") as i32;
        let distinguish_unknown_undefined_stereo = input["distinguish_unknown_undefined_stereo"]
            .as_i64()
            .expect("distinguish flag must be integer")
            as i32;
        let allocation_failure_ordinal = input["allocation_failure_ordinal"]
            .as_u64()
            .expect("allocation failure ordinal must be integer");
        let output_mode = input["output_mode"]
            .as_i64()
            .expect("output mode must be integer") as i32;
        let seed_old_allocations = input["seed_old_allocations"]
            .as_i64()
            .expect("seed_old_allocations must be integer")
            != 0;

        let mut heap = SourceHeap::default();
        let mut aux_bytes = text.bytes().map(|byte| byte as i8).collect::<Vec<_>>();
        aux_bytes.push(0);
        let aux_info = allocate_source_fixture(&mut heap, aux_bytes);
        let original_options =
            allocate_source_fixture(&mut heap, vec![b'-' as i8, b'X' as i8, 0, b'Z' as i8]);

        let mut original_atom_value = inchi_Atom::default();
        original_atom_value.x = -0.0;
        original_atom_value.y = 2.5;
        original_atom_value.z = -3.5;
        original_atom_value.elname[..2].copy_from_slice(&[b'F' as i8, 0]);
        original_atom_value.num_iso_H[0] = 4;
        original_atom_value.charge = -1;
        let original_atoms = if seed_old_allocations {
            allocate_source_fixture(&mut heap, vec![original_atom_value.clone()])
        } else {
            SourceMutPointer::null()
        };

        let original_stereo_value = inchi_Stereo0D {
            central_atom: 0,
            neighbor: [0; 4],
            type_: tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral as i8,
            parity: tagINCHIStereoParity0D_INCHI_PARITY_ODD as i8,
        };
        let original_stereo = if seed_old_allocations {
            allocate_source_fixture(&mut heap, vec![original_stereo_value.clone()])
        } else {
            SourceMutPointer::null()
        };

        let input_pointer = if output_mode == 2 {
            SourceMutPointer::null()
        } else {
            allocate_source_fixture(
                &mut heap,
                vec![inchi_Input {
                    atom: original_atoms,
                    stereo0D: original_stereo,
                    szOptions: original_options,
                    num_atoms: i16::from(seed_old_allocations),
                    num_stereo0D: i16::from(seed_old_allocations),
                }],
            )
        };
        let mut output = InchiInpData {
            pInp: input_pointer,
            bChiral: 0x5959_5959,
            szErrMsg: [0x59; STR_ERR_LEN as usize],
        };

        if allocation_failure_ordinal == 0 {
            heap.trace_source_allocations();
        } else {
            heap.fail_after_allocations(allocation_failure_ordinal - 1);
        }
        let output_argument = (output_mode != 1).then_some(&mut output);
        let status = if standard {
            Get_std_inchi_Input_FromAuxInfo(&mut heap, aux_info, do_not_add_h, output_argument)
        } else {
            Get_inchi_Input_FromAuxInfo(
                &mut heap,
                aux_info,
                do_not_add_h,
                distinguish_unknown_undefined_stereo,
                output_argument,
            )
        }
        .unwrap_or_else(|source_error| {
            panic!("Rust source model returned {source_error:?} for {case_id}")
        });
        let allocation_calls = heap.source_allocation_calls();
        let input_bytes = heap.slice(aux_info.as_const()).unwrap()[..=text.len()]
            .iter()
            .map(|byte| *byte as u8)
            .collect::<Vec<_>>();
        let original_atom_freed = !original_atoms.is_null()
            && heap.slice(original_atoms.as_const()) == Err(SourceHeapError::MissingAllocation);
        let original_stereo_freed = !original_stereo.is_null()
            && heap.slice(original_stereo.as_const()) == Err(SourceHeapError::MissingAllocation);
        let original_atom_unchanged = original_atoms.is_null()
            || heap
                .slice(original_atoms.as_const())
                .is_ok_and(|atoms| atoms.first() == Some(&original_atom_value));
        let original_stereo_unchanged = original_stereo.is_null()
            || heap
                .slice(original_stereo.as_const())
                .is_ok_and(|stereo| stereo.first() == Some(&original_stereo_value));

        let mut input_state = None;
        let mut atoms = Vec::new();
        let mut bond_fields = Vec::new();
        let mut stereo0d = Vec::new();
        let mut atom_pointer_state = Value::Null;
        let mut stereo_pointer_state = Value::Null;
        let mut options_pointer_state = Value::Null;
        if output_mode != 1 && !output.pInp.is_null() {
            let state = heap.slice(output.pInp.as_const()).unwrap()[0].clone();
            let atom_count = usize::try_from(state.num_atoms.max(0)).unwrap();
            if atom_count > 0 && !state.atom.is_null() {
                let atom_slice = &heap.slice(state.atom.as_const()).unwrap()[..atom_count];
                atoms = atom_slice
                    .iter()
                    .map(inchi_input_atom_json)
                    .collect::<Vec<_>>();
                bond_fields = atom_slice
                    .iter()
                    .enumerate()
                    .flat_map(|(atom_index, atom)| {
                        (0..usize::try_from(atom.num_bonds).unwrap()).map(move |slot| {
                            json!({
                                "atom_index": atom_index,
                                "slot": slot,
                                "neighbor": atom.neighbor[slot],
                                "bond_type": atom.bond_type[slot],
                                "bond_stereo": atom.bond_stereo[slot],
                            })
                        })
                    })
                    .collect::<Vec<_>>();
            }
            let stereo_count = usize::try_from(state.num_stereo0D.max(0)).unwrap();
            if stereo_count > 0 && !state.stereo0D.is_null() {
                stereo0d = heap.slice(state.stereo0D.as_const()).unwrap()[..stereo_count]
                    .iter()
                    .map(inchi_input_stereo_json)
                    .collect::<Vec<_>>();
            }
            atom_pointer_state = json!(if state.atom.is_null() {
                "null"
            } else if !original_atoms.is_null() && state.atom == original_atoms {
                "reused"
            } else {
                "allocated"
            });
            stereo_pointer_state = json!(if state.stereo0D.is_null() {
                "null"
            } else if !original_stereo.is_null() && state.stereo0D == original_stereo {
                "reused"
            } else {
                "allocated"
            });
            options_pointer_state = json!(if state.szOptions.is_null() {
                "null"
            } else if state.szOptions == original_options {
                "reused"
            } else {
                "allocated"
            });
            input_state = Some(state);
        }

        let output_present = output_mode != 1;
        let input_present = output_present && !output.pInp.is_null();
        let error_terminator = output_present
            .then(|| output.szErrMsg.iter().position(|byte| *byte == 0))
            .flatten();
        let error_text = error_terminator.map(|length| {
            String::from_utf8(
                output.szErrMsg[..length]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap()
        });
        let record = json!({
            "schema_version": "cosmolkit-inchi-official-c-v1",
            "oracle": {
                "project": "IUPAC InChI",
                "tag": "v1.07.5",
                "commit": "11a87982bb518f57ac013f0b258c283655e1ea1d",
                "tree_sha256": "4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd",
                "api_version": "1.07.5",
            },
            "case_id": case_id,
            "operation": if standard {
                "get_std_inchi_input_from_aux_info"
            } else {
                "get_inchi_input_from_aux_info"
            },
            "input": official["input"].clone(),
            "output": {
                "status": status,
                "allocation_calls": allocation_calls,
                "input_bytes": input_bytes,
                "input_pointer_state": output_present.then(|| if output.pInp.is_null() {
                    "null"
                } else if output.pInp == input_pointer {
                    "reused"
                } else {
                    "other"
                }),
                "b_chiral": output_present.then_some(output.bChiral),
                "error_text": if output_present { error_text } else { None },
                "error_text_terminated": output_present.then_some(error_terminator.is_some()),
                "error_bytes": output_present.then(|| output.szErrMsg.iter()
                    .map(|byte| *byte as u8).collect::<Vec<_>>()),
                "num_atoms": input_state.as_ref().map(|state| state.num_atoms),
                "num_stereo0d": input_state.as_ref().map(|state| state.num_stereo0D),
                "atoms": atoms,
                "bond_fields": bond_fields,
                "stereo0d": stereo0d,
                "atom_pointer_state": if input_present { atom_pointer_state } else { Value::Null },
                "stereo_pointer_state": if input_present { stereo_pointer_state } else { Value::Null },
                "options_pointer_state": if input_present { options_pointer_state } else { Value::Null },
                "original_atom_freed": original_atom_freed,
                "original_stereo_freed": original_stereo_freed,
                "original_atom_unchanged": original_atom_unchanged,
                "original_stereo_unchanged": original_stereo_unchanged,
            },
        });

        if let Some(mut state) = input_state {
            Free_inchi_Input(&mut heap, &mut state).unwrap();
            for pointer in [original_atoms, state.atom] {
                if !pointer.is_null() && heap.slice(pointer.as_const()).is_ok() {
                    inchi_free(&mut heap, pointer).unwrap();
                }
            }
            for pointer in [original_stereo, state.stereo0D] {
                if !pointer.is_null() && heap.slice(pointer.as_const()).is_ok() {
                    inchi_free(&mut heap, pointer).unwrap();
                }
            }
            for pointer in [original_options, state.szOptions] {
                if !pointer.is_null() && heap.slice(pointer.as_const()).is_ok() {
                    inchi_free(&mut heap, pointer).unwrap();
                }
            }
        } else {
            for pointer in [original_atoms] {
                if !pointer.is_null() && heap.slice(pointer.as_const()).is_ok() {
                    inchi_free(&mut heap, pointer).unwrap();
                }
            }
            for pointer in [original_stereo] {
                if !pointer.is_null() && heap.slice(pointer.as_const()).is_ok() {
                    inchi_free(&mut heap, pointer).unwrap();
                }
            }
            if heap.slice(original_options.as_const()).is_ok() {
                inchi_free(&mut heap, original_options).unwrap();
            }
        }
        if !input_pointer.is_null() {
            inchi_free(&mut heap, input_pointer).unwrap();
        }
        inchi_free(&mut heap, aux_info).unwrap();
        record
    }

    fn assert_aux_input_oracle(argument: &str, standard: bool, expected_fixed: &[&str]) {
        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg(argument)
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let official_records = String::from_utf8(oracle.stdout)
            .expect("official C oracle output must be UTF-8")
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).expect("oracle record must be JSON"))
            .collect::<Vec<_>>();
        let mut expected_case_ids = expected_fixed
            .iter()
            .map(|case_id| (*case_id).to_owned())
            .collect::<Vec<_>>();
        expected_case_ids.extend((1..=12).map(|ordinal| format!("allocation-merge-{ordinal}")));
        assert_eq!(
            official_records
                .iter()
                .map(|record| record["case_id"].as_str().unwrap().to_owned())
                .collect::<Vec<_>>(),
            expected_case_ids
        );
        for official in &official_records {
            let case_id = official["case_id"].as_str().unwrap();
            let rust = official_c_aux_input_rust_record(official, standard);
            assert_eq!(official, &rust, "exact field mismatch for {case_id}");
        }
    }

    #[test]
    fn official_c_oracle__ichilnct__get_inchi_input_fromauxinfo__line_89() {
        assert_aux_input_oracle(
            "--get-inchi-input-from-aux-info-records",
            false,
            &[
                "null-outer",
                "null-inner",
                "single-default-h",
                "single-no-add-h",
                "single-chiral",
                "merge-remap-chiral",
                "tetrahedral-unknown-default",
                "tetrahedral-unknown-distinct",
                "empty-structure",
                "eof",
                "missing-atom",
                "wrong-bonds",
                "aromatic-error",
                "old-reset-success",
                "old-reset-error",
                "max-atoms",
            ],
        );
    }

    #[test]
    fn official_c_oracle__ichilnct__get_std_inchi_input_fromauxinfo__line_79() {
        assert_aux_input_oracle(
            "--get-std-inchi-input-from-aux-info-records",
            true,
            &[
                "null-outer",
                "null-inner",
                "single-default-h",
                "single-no-add-h",
                "single-chiral",
                "merge-remap-chiral",
                "tetrahedral-unknown-default",
                "empty-structure",
                "eof",
                "missing-atom",
                "wrong-bonds",
                "aromatic-error",
                "old-reset-success",
                "old-reset-error",
                "max-atoms",
            ],
        );
    }

    #[test]
    fn source_port__ichilnct__free_inchi_input__line_187() {
        let mut heap = SourceHeap::default();
        let options = allocate_source_fixture(&mut heap, vec![b'-' as i8, b'S' as i8, 0]);
        let atom = allocate_source_fixture(&mut heap, vec![inchi_Atom::default(); 2]);
        let atom_alias = atom.as_const();
        let stereo = allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default()]);
        let stereo_alias = stereo.as_const();
        let mut input = inchi_Input {
            atom,
            stereo0D: stereo,
            szOptions: options,
            num_atoms: 2,
            num_stereo0D: 1,
        };

        assert_eq!(Free_inchi_Input(&mut heap, &mut input), Ok(()));
        assert!(input.atom.is_null());
        assert!(input.stereo0D.is_null());
        assert_eq!(input.szOptions, options);
        assert_eq!(input.num_atoms, 0);
        assert_eq!(input.num_stereo0D, 0);
        assert_eq!(
            heap.slice(atom_alias),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(stereo_alias),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(options.as_const()).unwrap(),
            &[b'-' as i8, b'S' as i8, 0]
        );
        assert_eq!(inchi_free(&mut heap, options), Ok(()));

        let null_options = SourceMutPointer::null();
        let mut empty_input = inchi_Input {
            atom: SourceMutPointer::null(),
            stereo0D: SourceMutPointer::null(),
            szOptions: null_options,
            num_atoms: 9,
            num_stereo0D: 7,
        };
        assert_eq!(Free_inchi_Input(&mut heap, &mut empty_input), Ok(()));
        assert!(empty_input.atom.is_null());
        assert!(empty_input.stereo0D.is_null());
        assert_eq!(empty_input.szOptions, null_options);
        assert_eq!(empty_input.num_atoms, 0);
        assert_eq!(empty_input.num_stereo0D, 0);

        let atom_allocation = allocate_source_fixture(&mut heap, vec![inchi_Atom::default(); 2]);
        let stereo_allocation = allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default()]);
        let mut invalid_atom_input = inchi_Input {
            atom: atom_allocation.offset(1).unwrap(),
            stereo0D: stereo_allocation,
            szOptions: SourceMutPointer::null(),
            num_atoms: 2,
            num_stereo0D: 1,
        };
        assert_eq!(
            Free_inchi_Input(&mut heap, &mut invalid_atom_input),
            Err(SourceHeapError::FreeOfInteriorPointer)
        );
        assert_eq!(invalid_atom_input.atom, atom_allocation.offset(1).unwrap());
        assert_eq!(invalid_atom_input.stereo0D, stereo_allocation);
        assert_eq!(invalid_atom_input.num_atoms, 2);
        assert_eq!(invalid_atom_input.num_stereo0D, 1);
        assert_eq!(heap.slice(atom_allocation.as_const()).unwrap().len(), 2);
        assert_eq!(heap.slice(stereo_allocation.as_const()).unwrap().len(), 1);
        assert_eq!(inchi_free(&mut heap, atom_allocation), Ok(()));
        assert_eq!(inchi_free(&mut heap, stereo_allocation), Ok(()));

        let atom_allocation = allocate_source_fixture(&mut heap, vec![inchi_Atom::default()]);
        let atom_alias = atom_allocation.as_const();
        let stereo_allocation =
            allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default(); 2]);
        let mut invalid_stereo_input = inchi_Input {
            atom: atom_allocation,
            stereo0D: stereo_allocation.offset(1).unwrap(),
            szOptions: SourceMutPointer::null(),
            num_atoms: 1,
            num_stereo0D: 2,
        };
        assert_eq!(
            Free_inchi_Input(&mut heap, &mut invalid_stereo_input),
            Err(SourceHeapError::FreeOfInteriorPointer)
        );
        assert!(invalid_stereo_input.atom.is_null());
        assert_eq!(
            invalid_stereo_input.stereo0D,
            stereo_allocation.offset(1).unwrap()
        );
        assert_eq!(invalid_stereo_input.num_atoms, 1);
        assert_eq!(invalid_stereo_input.num_stereo0D, 2);
        assert_eq!(
            heap.slice(atom_alias),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(heap.slice(stereo_allocation.as_const()).unwrap().len(), 2);
        assert_eq!(inchi_free(&mut heap, stereo_allocation), Ok(()));
    }

    #[test]
    fn source_port__ichilnct__free_std_inchi_input__line_182() {
        let mut heap = SourceHeap::default();
        let options = allocate_source_fixture(&mut heap, vec![b'-' as i8, b'S' as i8, 0]);
        let atom = allocate_source_fixture(&mut heap, vec![inchi_Atom::default(); 2]);
        let atom_alias = atom.as_const();
        let stereo = allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default()]);
        let stereo_alias = stereo.as_const();
        let mut input = inchi_Input {
            atom,
            stereo0D: stereo,
            szOptions: options,
            num_atoms: 2,
            num_stereo0D: 1,
        };

        assert_eq!(Free_std_inchi_Input(&mut heap, &mut input), Ok(()));
        assert!(input.atom.is_null());
        assert!(input.stereo0D.is_null());
        assert_eq!(input.szOptions, options);
        assert_eq!((input.num_atoms, input.num_stereo0D), (0, 0));
        assert_eq!(
            heap.slice(atom_alias),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(stereo_alias),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(options.as_const()).unwrap(),
            &[b'-' as i8, b'S' as i8, 0]
        );
        assert_eq!(inchi_free(&mut heap, options), Ok(()));

        let mut empty = inchi_Input {
            num_atoms: 3,
            num_stereo0D: 4,
            ..inchi_Input::default()
        };
        assert_eq!(Free_std_inchi_Input(&mut heap, &mut empty), Ok(()));
        assert_eq!(empty, inchi_Input::default());

        let atom = allocate_source_fixture(&mut heap, vec![inchi_Atom::default(); 2]);
        let stereo = allocate_source_fixture(&mut heap, vec![inchi_Stereo0D::default()]);
        let mut invalid = inchi_Input {
            atom: atom.offset(1).unwrap(),
            stereo0D: stereo,
            num_atoms: 1,
            num_stereo0D: 1,
            ..inchi_Input::default()
        };
        assert_eq!(
            Free_std_inchi_Input(&mut heap, &mut invalid),
            Err(SourceHeapError::FreeOfInteriorPointer)
        );
        assert_eq!(invalid.atom, atom.offset(1).unwrap());
        assert_eq!(invalid.stereo0D, stereo);
        assert_eq!((invalid.num_atoms, invalid.num_stereo0D), (1, 1));
        assert_eq!(inchi_free(&mut heap, atom), Ok(()));
        assert_eq!(inchi_free(&mut heap, stereo), Ok(()));
    }
}
