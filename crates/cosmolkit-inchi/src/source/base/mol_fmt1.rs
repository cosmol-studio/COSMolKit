use crate::source::base::ichi_io::inchi_fgetsLf;
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichiprt2::inchi_strtol;
use crate::source::base::mol_fmt2::{FreeMolfileData, MolfileFieldData, MolfileReadField};
use crate::source::base::mol_fmt3::{
    MolfileV3000Init, MolfileV3000ReadAtomsBlock, MolfileV3000ReadBondsBlock,
    MolfileV3000ReadCTABBeginAndCountsLine, MolfileV3000ReadTailOfCTAB,
};
use crate::source::base::mol_fmt4::{
    IntArray_Append, MolFmtSgroups_Alloc, MolFmtSgroups_Append, MolFmtSgroups_Free,
    MolFmtSgroups_GetIndexBySgroupId, SDFileSkipExtraData,
};
use crate::source::base::util::{
    dotify_non_printable_chars, extract_charges_and_radicals, get_atomic_mass, inchi_calloc,
    lrtrim, mystrncpy_slice, normalize_string, remove_one_lf,
};
use crate::source_types::{
    INCHI_IOSTREAM, MAXVAL, MOL_COORD, MOL_FMT_ATOM, MOL_FMT_BOND, MOL_FMT_CHAR_INT_DATA,
    MOL_FMT_CTAB, MOL_FMT_DATA, MOL_FMT_DOUBLE_DATA, MOL_FMT_HEADER_BLOCK, MOL_FMT_INPLINELEN,
    MOL_FMT_JUMP_TO_RIGHT, MOL_FMT_LONG_INT_DATA, MOL_FMT_M_CONN_EU, MOL_FMT_M_CONN_HH,
    MOL_FMT_M_CONN_HT, MOL_FMT_M_CONN_NON, MOL_FMT_M_SST_ALT, MOL_FMT_M_SST_BLK, MOL_FMT_M_SST_NON,
    MOL_FMT_M_SST_RAN, MOL_FMT_M_STY_COP, MOL_FMT_M_STY_CRO, MOL_FMT_M_STY_MER, MOL_FMT_M_STY_MOD,
    MOL_FMT_M_STY_MON, MOL_FMT_M_STY_NON, MOL_FMT_M_STY_SRU, MOL_FMT_MAX_VALUE_LEN,
    MOL_FMT_MAXLINELEN, MOL_FMT_SGROUPS, MOL_FMT_SHORT_INT_DATA, MOL_FMT_STRING_DATA,
    MOL_FMT_v3000, POLYMERS_NO, RADICAL_DOUBLET, SD_FMT_END_OF_DATA, SourceConstPointer,
    SourceHeap, SourceHeapError, SourceMutPointer, ZERO_ATW_DIFF,
};

fn mol_fmt1_c_string_eq(value: &[i8], expected: &[u8]) -> bool {
    let length = value
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(value.len());
    length == expected.len()
        && value[..length]
            .iter()
            .zip(expected)
            .all(|(left, right)| *left as u8 == *right)
}

fn mol_fmt1_add_ascii_message(
    target: Option<&mut [i8]>,
    message: &[u8],
) -> Result<i32, SourceHeapError> {
    let terminated = message
        .iter()
        .map(|byte| *byte as i8)
        .chain(std::iter::once(0))
        .collect::<Vec<_>>();
    AddErrorMessage(target, Some(&terminated))
}

fn mol_fmt1_treat_error(
    error: &mut i32,
    new_error: i32,
    target: Option<&mut [i8]>,
    message: &[u8],
) -> Result<(), SourceHeapError> {
    if *error == 0 && new_error != 0 {
        *error = new_error;
    }
    let _ = mol_fmt1_add_ascii_message(target, message)?;
    Ok(())
}

fn mol_fmt1_trim_array<const N: usize>(
    heap: &mut SourceHeap,
    value: &mut [i8; N],
) -> Result<i32, SourceHeapError> {
    let pointer = heap.allocate_model_storage(value.to_vec())?;
    let mut length = 0_i32;
    lrtrim(heap, pointer, Some(&mut length))?;
    value.copy_from_slice(heap.slice(pointer.as_const())?);
    heap.free(pointer)?;
    Ok(length)
}

fn mol_fmt1_sgroup_pointer(
    heap: &SourceHeap,
    sgroups: &MOL_FMT_SGROUPS,
    index: i32,
) -> Result<SourceMutPointer<crate::source_types::MOL_FMT_SGROUP>, SourceHeapError> {
    let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    heap.slice(sgroups.group.as_const())?
        .get(index)
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn ReadMolfile(
    heap: &mut SourceHeap,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    OnlyHeaderBlock: SourceMutPointer<MOL_FMT_HEADER_BLOCK>,
    OnlyCTab: SourceMutPointer<MOL_FMT_CTAB>,
    bGetOrigCoord: i32,
    treat_polymers: i32,
    treat_NPZz: i32,
    pname: SourceMutPointer<i8>,
    lname: i32,
    mut Id: Option<&mut u64>,
    pSdfLabel: SourceConstPointer<i8>,
    pSdfValue: SourceMutPointer<i8>,
    err: &mut i32,
    mut pStrErr: Option<&mut [i8]>,
    bNoWarnings: i32,
) -> Result<SourceMutPointer<MOL_FMT_DATA>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:92 ReadMolfile
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap pointer ownership adds allocation-map overhead.
    /*
MOL_FMT_DATA* ReadMolfile( INCHI_IOSTREAM *inp_file,
                           MOL_FMT_HEADER_BLOCK *OnlyHeaderBlock,
                           MOL_FMT_CTAB *OnlyCTab,
                           int bGetOrigCoord,
                           int treat_polymers,
                           int treat_NPZz,
                           char *pname,
                           int lname,
                           unsigned long *Id,
                           const char *pSdfLabel,
                           char *pSdfValue,
                           int *err,
                           char *pStrErr,
                           int bNoWarnings )
{
    MOL_FMT_DATA* mfdata;

    if (pname && lname)
    {
        pname[0] = '\0';
    }
    if (Id)
    {
        *Id = 0LU;  /* ignore for now */
    }

    mfdata = MolfileReadDataLines( inp_file, OnlyHeaderBlock, OnlyCTab,
                                   bGetOrigCoord, treat_polymers,
                                   err, pStrErr, bNoWarnings );

    if (*err < 0)
    {
        /* read OK, and end of data encountered */
        *err = -*err;
    }
    else
    {
        /* unnecessary extra data may have present in SDF; skip them for now */
        int ret_skip_extras = SDFileSkipExtraData( inp_file, Id, NULL, 0,
                                                   pname, lname, *err,
                                                   pSdfLabel, pSdfValue,
                                                   pStrErr, bNoWarnings);


        if (ret_skip_extras)
        {
            /* important to continue to the next good structure */
            *err = ret_skip_extras;
        }
    }

    /*  Treat star/Zz atoms if present.
        This first-line treating affects only Molfile input.
        No direct input in API calls, no InChI strings.
        These will be processed on checking internal data structs.
    */
    if (mfdata)
    {
        int nzz = 0; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        int pseudos_allowed = (treat_NPZz == 1) || (treat_polymers != POLYMERS_NO);		
        nzz = MolfileTreatPseudoElementAtoms( &mfdata->ctab,
                                              pseudos_allowed,
                                              err,
                                              pStrErr ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    }

    return mfdata;
}
    */
    // END INCHI C FUNCTION: ReadMolfile
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: ReadMolfile
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; unsigned long is u64 and POLYMERS_NO is zero.
    // INCHI✔️❌: MolfileReadDataLines, SDFileSkipExtraData, and MolfileTreatPseudoElementAtoms are completed source ports and are called in source order with the same error mutation rules.
    // INCHI✔️❌: The source's OnlyHeaderBlock cast is defined only when its pointer denotes the hdr field at the start of storage also containing a following ctab; the Rust model uses the explicit OnlyCTab descriptor for that modeled path and does not emulate undefined adjacent-stack access.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: ReadMolfile

    if !pname.is_null() && lname != 0 {
        *heap
            .slice_mut(pname)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    }
    if let Some(value) = Id.as_deref_mut() {
        *value = 0;
    }

    let mfdata = MolfileReadDataLines(
        heap,
        inp_file.as_deref_mut(),
        OnlyHeaderBlock,
        OnlyCTab,
        bGetOrigCoord,
        treat_polymers,
        err,
        pStrErr.as_deref_mut(),
        bNoWarnings,
    )?;

    if *err < 0 {
        *err = err.wrapping_neg();
    } else {
        let ret_skip_extras = SDFileSkipExtraData(
            heap,
            inp_file.as_deref_mut(),
            Id.as_deref_mut(),
            SourceMutPointer::null(),
            0,
            pname,
            lname,
            *err,
            pSdfLabel,
            pSdfValue,
            pStrErr.as_deref_mut(),
            bNoWarnings,
        )?;
        if ret_skip_extras != 0 {
            *err = ret_skip_extras;
        }
    }

    if !mfdata.is_null() {
        let mut ctab = if OnlyHeaderBlock.is_null() {
            heap.slice(mfdata.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .ctab
                .clone()
        } else if !OnlyCTab.is_null() {
            heap.slice(OnlyCTab.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone()
        } else {
            return Err(SourceHeapError::PointerOutOfBounds);
        };
        let pseudos_allowed = i32::from(treat_NPZz == 1 || treat_polymers != POLYMERS_NO as i32);
        let _ = MolfileTreatPseudoElementAtoms(
            heap,
            &mut ctab,
            pseudos_allowed,
            err,
            pStrErr.as_deref_mut(),
        )?;
    }

    Ok(mfdata)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MolfileReadDataLines(
    heap: &mut SourceHeap, mut inp_file: Option<&mut INCHI_IOSTREAM>,
    only_header: SourceMutPointer<MOL_FMT_HEADER_BLOCK>,
    only_ctab: SourceMutPointer<MOL_FMT_CTAB>, get_coords: i32, treat_polymers: i32,
    err: &mut i32, mut errors: Option<&mut [i8]>, no_warnings: i32,
) -> Result<SourceMutPointer<MOL_FMT_DATA>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:165 MolfileReadDataLines
    // INCHI✔️❌: complete source frame follows verbatim; checked heap ownership and local state synchronization add overhead.
    /*
MOL_FMT_DATA * MolfileReadDataLines( INCHI_IOSTREAM *inp_file,
                                     MOL_FMT_HEADER_BLOCK *OnlyHeaderBlock,
                                     MOL_FMT_CTAB *OnlyCTab,
                                     int bGetOrigCoord,
                                     int treat_polymers,
                                     int *err,
                                     char *pStrErr,
                                     int bNoWarnings )
{
    int n_alloc_atoms;
    MOL_FMT_CTAB  ctab, *pCtab = NULL;
    MOL_FMT_HEADER_BLOCK *pHdr = NULL;
    MOL_FMT_DATA *mfdata = NULL;
    int retcode, prevcode;
    int data_ended = 0, should_read_all = 0;

    if (!OnlyHeaderBlock)
    {
        should_read_all = 1;
    }

    /* djb-rwth: removing redundant code */
    /* djb-rwth: addressing coverity ID #499502 -- TREAT_ERR_AND_FIN properly used in all cases */
    *err = 0;

    if (should_read_all)
    {

        mfdata = (MOL_FMT_DATA*) inchi_calloc( 1, sizeof( MOL_FMT_DATA ) );
        if (!mfdata)
        {
            retcode = 1;
            AddErrorMessage( pStrErr, "Out of RAM" );
            goto err_fin;
        }

        pHdr = &mfdata->hdr;
        pCtab = &mfdata->ctab;
    }
    else
    {
        pHdr = OnlyHeaderBlock;
        pCtab = OnlyCTab ? OnlyCTab : &ctab;
        memset( pHdr, 0, sizeof( MOL_FMT_HEADER_BLOCK ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        memset( pCtab, 0, sizeof( MOL_FMT_CTAB ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    }

    pCtab->bonds = NULL;
    pCtab->atoms = NULL;
    pCtab->coords = NULL;
    pCtab->v3000 = NULL;

    /* Header lines */
    retcode = MolfileReadHeaderLines( pHdr, inp_file, pStrErr );
    if (retcode)
    {
        /*  most likely end of file */
        retcode += 10;
        goto err_fin;
    }

    /* Read counts line and also check if we deal with V3000 Molfile */
    retcode = MolfileReadCountsLine( pCtab, inp_file, pStrErr );
    if (retcode)
    {
        retcode += 20;
        goto err_fin;
    }

    if (pCtab->v3000)
    {

        /*    In V3000, the previously read counts should be neglected
            and new counts line (1st in Ctab) should be read preceded
            by "M  V30 BEGIN CTAB"                                    */

        retcode = MolfileV3000ReadCTABBeginAndCountsLine( pCtab, inp_file, pStrErr );
        if (retcode)
        {
            retcode += 70;
            TREAT_ERR_AND_FIN( retcode, 1, err_fin, pStrErr );
        }

        retcode = MolfileV3000Init( pCtab, pStrErr );
        if (retcode)
        {
            retcode += 70;
            TREAT_ERR_AND_FIN( retcode, 1, err_fin, pStrErr );
        }
    }

    /* Atomic block */
    if (should_read_all)
    {
        n_alloc_atoms = inchi_max( mfdata->ctab.n_atoms, 1 );
        mfdata->ctab.atoms = (MOL_FMT_ATOM*)
            inchi_calloc( n_alloc_atoms, sizeof( MOL_FMT_ATOM ) );
        if (!mfdata->ctab.atoms)
        {
            retcode = 2;
            TREAT_ERR_AND_FIN( retcode, 2, err_fin, "Out of RAM" );
        }

        if (bGetOrigCoord)
        {
            mfdata->ctab.coords = (MOL_COORD*)
                inchi_calloc( n_alloc_atoms, sizeof( MOL_COORD ) );
            if (!mfdata->ctab.coords)
            {
                retcode = 2;
                TREAT_ERR_AND_FIN( retcode, 2, err_fin, "Out of RAM" );
            }
        }
    }

    if (!pCtab->v3000)
    {
        retcode = MolfileReadAtomsBlock( pCtab, inp_file, retcode, pStrErr );
    }
    else
    {
        retcode = MolfileV3000ReadAtomsBlock( pCtab, inp_file, retcode, pStrErr );
    }
    if (retcode)
    {
        if (retcode < 0)
        {
            retcode = -retcode;
            data_ended = 1;
        }
        retcode += 30;
        /* goto err_fin; */
    }

    /* Bonds block */
    if (should_read_all && retcode < 30)
    {
        if (!data_ended)
        {
            int n_alloc_bonds = inchi_max( mfdata->ctab.n_bonds, 1 );
            mfdata->ctab.bonds = (MOL_FMT_BOND *)
                inchi_calloc( n_alloc_bonds, sizeof( MOL_FMT_BOND ) );
            if (!mfdata->ctab.bonds)
            {
                /* can't allocate bonds structure */
                retcode = 3;
                TREAT_ERR_AND_FIN( retcode, 3, err_fin, "Out of RAM" );
            }
        }
    }

    prevcode = retcode;
    if (!data_ended)
    {
        if (!pCtab->v3000)
        {
            retcode = MolfileReadBondsBlock( pCtab, inp_file, retcode, pStrErr );
        }
        else
        {
            retcode = MolfileV3000ReadBondsBlock( pCtab, inp_file, retcode, pStrErr );
        }
        if (retcode)
        {
            if (retcode < 0)
            {
                retcode = -retcode;
                data_ended = 1;
            }
            retcode = prevcode ? prevcode
                : retcode + 40;
        }


        /* SGroup, 3D, link line(s), collections, END CTAB */
        if (pCtab->v3000)
        {
            retcode = MolfileV3000ReadTailOfCTAB( pCtab, inp_file, retcode, pStrErr );
            if (retcode)
            {
                if (retcode < 0)
                {
                    retcode = -retcode;
                    data_ended = 1;
                }
                retcode = prevcode ? prevcode : retcode + 70;
            }
        }
    }
    prevcode = retcode;

    /* SText */
    if (!data_ended)
    {
        retcode = MolfileReadSTextBlock( pCtab, inp_file, retcode, pStrErr );
        if (retcode)
        {
            retcode = prevcode ? prevcode : retcode + 50;
        }
    }
    prevcode = retcode;

    /* Prop block */
    if (!data_ended)
    {
        retcode = MolfileReadPropBlock( pCtab, pHdr, inp_file, treat_polymers, retcode, pStrErr, bNoWarnings );

        if (retcode)
        {
            if (retcode < 0)
            {
                retcode = -retcode;
                data_ended = 1;
            }
            retcode = prevcode ? prevcode : retcode + 60;
        }
    }

    /* Check that all valences are in allowed range ( <=MAXVAL, currently 20 )          */
    if (1)
#ifdef TARGET_LIB_FOR_WINCHI
        if (pCtab && pCtab->atoms)
#endif
    {

        int i;
        for (i = 0; i < pCtab->n_atoms; i++)
        {
            if (pCtab->atoms) /* djb-rwth: fixing a NULL pointer dereference */
            {
                if (pCtab->atoms[i].valence > MAXVAL)
                {
                    retcode = 70 + 9;
                    TREAT_ERR( retcode, 0, "Too large input atomic valence" );
                    break;
                }
            }
        }
#if ( FIX_CURE53_ISSUE_NULL_DEREFERENCE_MAKE_A_COPY_OF_T_GROUP_INFO==1 || defined(FIX_IMPOSSIBLE_H_ISOTOPE_BUG) )
        /* Do not eat H isotopes other than [H,D,T] */
        for (i = 0; i < pCtab->n_atoms; i++)
        {
            if (pCtab->atoms) /* djb-rwth: fixing a NULL pointer dereference */
            {
                int dmass = pCtab->atoms[i].mass_difference;
                if ((!strcmp(pCtab->atoms[i].symbol, "H") && dmass != 0 && dmass != 1 && dmass != 2 && dmass != 127) ||
                    (!strcmp(pCtab->atoms[i].symbol, "D") && dmass != 0 && dmass != 1 && dmass != -1) ||
                    (!strcmp(pCtab->atoms[i].symbol, "T") && dmass != 0 && dmass != -1 && dmass != -2))
                {
                    retcode = 70 + 8;
                    TREAT_ERR(retcode, 0, "Unacceptable isotope of hydrogen");
                    break;
                }
            }
        }
#endif
    }

err_fin:
    *err = data_ended ? -retcode : retcode;

    if (should_read_all)
    {
        if (retcode)
        {
            mfdata = FreeMolfileData( mfdata );        /* delete all results */
        }
        return mfdata;
    }
    else
    {
        if (retcode)
        {
            return NULL;
        }
        else
        {
            return (MOL_FMT_DATA*) OnlyHeaderBlock;
        }
    }
}
*/
    // END INCHI C FUNCTION: MolfileReadDataLines
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadDataLines
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; TARGET_LIB_FOR_WINCHI is undefined.
    // INCHI✔️❌: FIX_CURE53_ISSUE_NULL_DEREFERENCE_MAKE_A_COPY_OF_T_GROUP_INFO == 1 and FIX_IMPOSSIBLE_H_ISOTOPE_BUG == 1.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadDataLines
    let all = only_header.is_null();
    let mut data = SourceMutPointer::null();
    let mut hdr = MOL_FMT_HEADER_BLOCK::default();
    let mut ctab = MOL_FMT_CTAB::default();
    let mut rc = 0_i32;
    let mut ended = false;
    *err = 0;
    'read: {
        if all {
            data = match inchi_calloc::<MOL_FMT_DATA>(heap, 1, std::mem::size_of::<MOL_FMT_DATA>() as u64) {
                Ok(p) => p,
                Err(SourceHeapError::AllocationFailed) => {
                    rc = 1; let _ = mol_fmt1_add_ascii_message(errors.as_deref_mut(), b"Out of RAM")?; break 'read;
                }
                Err(e) => return Err(e),
            };
        }
        ctab.bonds = SourceMutPointer::null();
        ctab.atoms = SourceMutPointer::null();
        ctab.coords = SourceMutPointer::null();
        ctab.v3000 = SourceMutPointer::null();

        rc = MolfileReadHeaderLines(heap, &mut hdr, inp_file.as_deref_mut(), errors.as_deref_mut())?;
        if rc != 0 { rc = rc.wrapping_add(10); break 'read; }
        rc = MolfileReadCountsLine(heap, &mut ctab, inp_file.as_deref_mut(), errors.as_deref_mut())?;
        if rc != 0 { rc = rc.wrapping_add(20); break 'read; }
        if !ctab.v3000.is_null() {
            rc = MolfileV3000ReadCTABBeginAndCountsLine(heap, &mut ctab, inp_file.as_deref_mut(), errors.as_deref_mut())?;
            if rc != 0 { rc = rc.wrapping_add(70); break 'read; }
            rc = MolfileV3000Init(heap, &mut ctab, errors.as_deref_mut())?;
            if rc != 0 { rc = rc.wrapping_add(70); break 'read; }
        }
        if all {
            let n = ctab.n_atoms.max(1) as u64;
            ctab.atoms = match inchi_calloc::<MOL_FMT_ATOM>(heap, n, std::mem::size_of::<MOL_FMT_ATOM>() as u64) {
                Ok(p) => p,
                Err(SourceHeapError::AllocationFailed) => {
                    rc = 2; let _ = mol_fmt1_add_ascii_message(errors.as_deref_mut(), b"Out of RAM")?; break 'read;
                }
                Err(e) => return Err(e),
            };
            if get_coords != 0 {
                ctab.coords = match inchi_calloc::<MOL_COORD>(heap, n, std::mem::size_of::<MOL_COORD>() as u64) {
                    Ok(p) => p,
                    Err(SourceHeapError::AllocationFailed) => {
                        rc = 2; let _ = mol_fmt1_add_ascii_message(errors.as_deref_mut(), b"Out of RAM")?; break 'read;
                    }
                    Err(e) => return Err(e),
                };
            }
        }
        rc = if ctab.v3000.is_null() {
            MolfileReadAtomsBlock(heap, &mut ctab, inp_file.as_deref_mut(), rc, errors.as_deref_mut())?
        } else {
            MolfileV3000ReadAtomsBlock(heap, &mut ctab, inp_file.as_deref_mut(), rc, errors.as_deref_mut())?
        };
        if rc != 0 {
            if rc < 0 { rc = rc.wrapping_abs(); ended = true; }
            rc = rc.wrapping_add(30);
        }
        if all && rc < 30 && !ended {
            let n = ctab.n_bonds.max(1) as u64;
            ctab.bonds = match inchi_calloc::<MOL_FMT_BOND>(heap, n, std::mem::size_of::<MOL_FMT_BOND>() as u64) {
                Ok(p) => p,
                Err(SourceHeapError::AllocationFailed) => {
                    rc = 3; let _ = mol_fmt1_add_ascii_message(errors.as_deref_mut(), b"Out of RAM")?; break 'read;
                }
                Err(e) => return Err(e),
            };
        }
        let mut previous = rc;
        if !ended {
            rc = if ctab.v3000.is_null() {
                MolfileReadBondsBlock(heap, &mut ctab, inp_file.as_deref_mut(), rc, errors.as_deref_mut())?
            } else {
                MolfileV3000ReadBondsBlock(heap, &mut ctab, inp_file.as_deref_mut(), rc, errors.as_deref_mut())?
            };
            if rc != 0 {
                if rc < 0 { rc = rc.wrapping_abs(); ended = true; }
                rc = if previous != 0 { previous } else { rc.wrapping_add(40) };
            }
            if !ctab.v3000.is_null() {
                rc = MolfileV3000ReadTailOfCTAB(heap, &mut ctab, inp_file.as_deref_mut(), rc, errors.as_deref_mut())?;
                if rc != 0 {
                    if rc < 0 { rc = rc.wrapping_abs(); ended = true; }
                    rc = if previous != 0 { previous } else { rc.wrapping_add(70) };
                }
            }
        }
        previous = rc;
        if !ended {
            rc = MolfileReadSTextBlock(heap, &ctab, inp_file.as_deref_mut(), rc, errors.as_deref_mut())?;
            if rc != 0 { rc = if previous != 0 { previous } else { rc.wrapping_add(50) }; }
        }
        previous = rc;
        if !ended {
            rc = MolfileReadPropBlock(heap, &mut ctab, &mut hdr, inp_file.as_deref_mut(), treat_polymers, rc, errors.as_deref_mut(), no_warnings)?;
            if rc != 0 {
                if rc < 0 { rc = rc.wrapping_abs(); ended = true; }
                rc = if previous != 0 { previous } else { rc.wrapping_add(60) };
            }
        }
        if !ctab.atoms.is_null() {
            let n = usize::try_from(ctab.n_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atoms = heap.slice(ctab.atoms.as_const())?.get(..n).ok_or(SourceHeapError::PointerOutOfBounds)?;
            for atom in atoms {
                if atom.valence > MAXVAL as i8 {
                    rc = 79; let _ = mol_fmt1_add_ascii_message(errors.as_deref_mut(), b"Too large input atomic valence")?; break;
                }
            }
            for atom in atoms {
                let m = atom.mass_difference;
                let bad = (mol_fmt1_c_string_eq(&atom.symbol, b"H") && !matches!(m, 0 | 1 | 2 | 127))
                    || (mol_fmt1_c_string_eq(&atom.symbol, b"D") && !matches!(m, -1 | 0 | 1))
                    || (mol_fmt1_c_string_eq(&atom.symbol, b"T") && !matches!(m, -2 | -1 | 0));
                if bad {
                    rc = 78; let _ = mol_fmt1_add_ascii_message(errors.as_deref_mut(), b"Unacceptable isotope of hydrogen")?; break;
                }
            }
        }
    }
    *err = if ended { rc.wrapping_neg() } else { rc };
    if all {
        if !data.is_null() {
            let out = heap.slice_mut(data)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?;
            out.hdr = hdr; out.ctab = ctab;
        }
        if rc != 0 { data = FreeMolfileData(heap, data)?; }
        Ok(data)
    } else {
        *heap.slice_mut(only_header)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)? = hdr;
        if !only_ctab.is_null() {
            *heap.slice_mut(only_ctab)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)? = ctab;
        }
        if rc != 0 { Ok(SourceMutPointer::null()) } else { Ok(only_header.cast()) }
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileReadHeaderLines(
    heap: &mut SourceHeap,
    hdr: &mut MOL_FMT_HEADER_BLOCK,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    _p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:451 MolfileReadHeaderLines
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap stack-buffer modeling adds overhead.
    /*
int MolfileReadHeaderLines( MOL_FMT_HEADER_BLOCK *hdr,
                            INCHI_IOSTREAM *inp_file,
                            char *pStrErr )
{
/* All input lines can have up to 80 characters */
/* Header Block */

    char line[MOL_FMT_INPLINELEN]; /* + cr +lf +zero termination + reserve */
    int  err = 0, len; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    const int  line_len = sizeof( line );
    char *p;


        /* Header line #1: name */
    p = inchi_fgetsLf( line, line_len, inp_file );
    if (!p)
    {
        /* can't read the input file line */
        err = 1;
        /* AddErrorMessage( pStrErr, "Can't read header block name line" ); */
        goto err_fin;
    }

    remove_one_lf( line );

    /* -- Disabled to relax strictness: allow > 80 chars names. */
    /*
    if ( line[MOL_FMT_MAXLINELEN] )
    {
        //err = 2;             too long line
        goto err_fin;
    }
    */

    len = MolfileReadField( hdr->molname,
                            sizeof( hdr->molname ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* Header line #2 */
    p = inchi_fgetsLf( line, line_len, inp_file );
    if (!p)
    {
        /* can't read the input file line */
        err = 3;
        /* AddErrorMessage( pStrErr, "Can't read header block line 2" ); */
        goto err_fin;
    }

    remove_one_lf( line );

    /* -- Disabled to relax strictness: allow > 80 chars names. */
    /*
    if ( line[MOL_FMT_MAXLINELEN] )
    {
        err = 4;             too long input file line
        goto err_fin;
    }
    */

    len = MolfileReadField( hdr->user_initls,
                            sizeof( hdr->user_initls ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( hdr->prog_name,
                            sizeof( hdr->prog_name ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /*------------ Relax strictness -----------------------*/
    len = MolfileReadField( &hdr->month, 2, MOL_FMT_CHAR_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->day, 2, MOL_FMT_CHAR_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->year, 2, MOL_FMT_CHAR_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->hour, 2, MOL_FMT_CHAR_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->minute, 2, MOL_FMT_CHAR_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( hdr->dim_code, sizeof( hdr->dim_code ) - 1, MOL_FMT_STRING_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->scaling_factor1, 2, MOL_FMT_SHORT_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->scaling_factor2, 10, MOL_FMT_DOUBLE_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->energy, 12, MOL_FMT_DOUBLE_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    len = MolfileReadField( &hdr->internal_regno, 6, MOL_FMT_LONG_INT_DATA, &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* Save the whole line 2 */
    p = line;
    len = MolfileReadField( hdr->line2,
                            sizeof( hdr->line2 ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* Header line #3: comment */
    p = inchi_fgetsLf( line, line_len, inp_file );

    if (!p)
    {
        err = 7;             /* can't read the line */
        /* AddErrorMessage( pStrErr, "Can't read header block comment line" ); */
        goto err_fin;
    }
    remove_one_lf( line );
    /* -- Disabled to relax strictness: allow > 80 chars comments.
    if ( line[MOL_FMT_MAXLINELEN] ){
        err = 8;             too long line
        goto err_fin;
    }
    */
    len = MolfileReadField( hdr->comment,
                            sizeof( hdr->comment ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

err_fin:

    return err;
}

    */
    // END INCHI C FUNCTION: MolfileReadHeaderLines
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadHeaderLines
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; both long-line rejection blocks are comments and inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadHeaderLines

    let line = heap.allocate_model_storage(vec![0_i8; MOL_FMT_INPLINELEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let mut p = inchi_fgetsLf(heap, line, MOL_FMT_INPLINELEN as i32, inp_file.as_deref_mut())?;
        if p.is_null() { return Ok(1); }
        remove_one_lf(heap, line)?;
        let _ = MolfileReadField(heap, MolfileFieldData::String(&mut hdr.molname), 200, i32::from(MOL_FMT_STRING_DATA), &mut p)?;

        p = inchi_fgetsLf(heap, line, MOL_FMT_INPLINELEN as i32, inp_file.as_deref_mut())?;
        if p.is_null() { return Ok(3); }
        remove_one_lf(heap, line)?;
        let _ = MolfileReadField(heap, MolfileFieldData::String(&mut hdr.user_initls), 2, i32::from(MOL_FMT_STRING_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::String(&mut hdr.prog_name), 8, i32::from(MOL_FMT_STRING_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::Char(&mut hdr.month), 2, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::Char(&mut hdr.day), 2, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::Char(&mut hdr.year), 2, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::Char(&mut hdr.hour), 2, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::Char(&mut hdr.minute), 2, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::String(&mut hdr.dim_code), 2, i32::from(MOL_FMT_STRING_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::Short(&mut hdr.scaling_factor1), 2, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::Double(&mut hdr.scaling_factor2), 10, i32::from(MOL_FMT_DOUBLE_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::Double(&mut hdr.energy), 12, i32::from(MOL_FMT_DOUBLE_DATA), &mut p)?;
        let _ = MolfileReadField(heap, MolfileFieldData::Long(&mut hdr.internal_regno), 6, i32::from(MOL_FMT_LONG_INT_DATA), &mut p)?;
        p = line;
        let _ = MolfileReadField(heap, MolfileFieldData::String(&mut hdr.line2), 200, i32::from(MOL_FMT_STRING_DATA), &mut p)?;

        p = inchi_fgetsLf(heap, line, MOL_FMT_INPLINELEN as i32, inp_file.as_deref_mut())?;
        if p.is_null() { return Ok(7); }
        remove_one_lf(heap, line)?;
        let _ = MolfileReadField(heap, MolfileFieldData::String(&mut hdr.comment), 80, i32::from(MOL_FMT_STRING_DATA), &mut p)?;
        Ok(0)
    })();
    let cleanup = heap.free(line);
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileReadCountsLine(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    inp_file: Option<&mut INCHI_IOSTREAM>,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:569 MolfileReadCountsLine
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap stack modeling adds overhead.
    /*
int MolfileReadCountsLine( MOL_FMT_CTAB* ctab,
                           INCHI_IOSTREAM *inp_file,
                           char *pStrErr )
{
    char *p;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    int   err = 0, len; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* djb-rwth: addressing coverity ID #499502 -- TREAT_ERR properly used in all cases */

    p = inchi_fgetsLf( line, line_len, inp_file );

    if (!p)
    {
        TREAT_ERR_AND_FIN( err, 1, err_fin, "Cannot read counts line" );
        /* can't read the input file line */
    }

    remove_one_lf( line );

    if (line[MOL_FMT_MAXLINELEN])
    {
        TREAT_ERR( err, 0, "Too long counts line" );  /* too long input file line */
    }

    if (0 > MolfileReadField( &ctab->n_atoms, 3, MOL_FMT_SHORT_INT_DATA, &p ) /* V2000 only: short int */
         || 0 > MolfileReadField( &ctab->n_bonds, 3, MOL_FMT_SHORT_INT_DATA, &p ) /* V2000 only: short int */

#if ( MOL_FMT_QUERY == MOL_FMT_PRESENT )
         || 0 > MolfileReadField( &ctab->n_atom_lists, 3, MOL_FMT_SHORT_INT_DATA, &p )
#else
         || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

         || 0 > MolfileReadField( NULL, /*obsolete*/      3, MOL_FMT_JUMP_TO_RIGHT, &p )
         || 0 > MolfileReadField( &ctab->chiral_flag, 3, MOL_FMT_CHAR_INT_DATA, &p )
         || 0 > MolfileReadField( &ctab->n_stext_entries, 3, MOL_FMT_SHORT_INT_DATA, &p )

#if ( MOL_FMT_CPSS == MOL_FMT_PRESENT )
         || 0 > MolfileReadField( &ctab->n_reaction_components_plus_1, 3, MOL_FMT_SHORT_INT_DATA, &p )
         || 0 > MolfileReadField( &ctab->n_reactants, 3, MOL_FMT_SHORT_INT_DATA, &p )
         || 0 > MolfileReadField( &ctab->n_products, 3, MOL_FMT_SHORT_INT_DATA, &p )
         || 0 > MolfileReadField( &ctab->n_intermediates, 3, MOL_FMT_SHORT_INT_DATA, &p )
#else
         || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
         || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
         || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
         || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

         || 0 > MolfileReadField( &ctab->n_property_lines, 3, MOL_FMT_SHORT_INT_DATA, &p ))
    {
        err = 3;  /* can't interpret counts line */
        TREAT_ERR( err, 3, "Cannot interpret counts line:" );  /* too long input file line */
        dotify_non_printable_chars( line );
        AddErrorMessage( pStrErr, line );
        goto err_fin;
    }

    /* Get CTFile version (V2000 or other) */
    len = MolfileReadField( ctab->version_string,
                            sizeof( ctab->version_string ) - 1,
                            MOL_FMT_STRING_DATA,
                            &p ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* Allocate additional space if V3000 */
    if (!strcmp( ctab->version_string, "V3000" ))
    {
        ctab->v3000 = (MOL_FMT_v3000*) inchi_calloc( 1, sizeof( MOL_FMT_v3000 ) );

        if (!ctab->v3000)
        {
            AddErrorMessage( pStrErr, "Out of RAM" );
            return -1;
        }
    }
    else
        ctab->v3000 = NULL; /* paranoia */

    /* Polymer Sgroups */
    MolFmtSgroups_Alloc( &( ctab->sgroups ), 1 );

err_fin:

    return err;
}
    */
    // END INCHI C FUNCTION: MolfileReadCountsLine
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadCountsLine
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; MOL_FMT_QUERY and MOL_FMT_CPSS are MOL_FMT_ABSENT.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadCountsLine

    const READ_ERROR: &[i8] = &[b'C' as i8,b'a' as i8,b'n' as i8,b'n' as i8,b'o' as i8,b't' as i8,b' ' as i8,b'r' as i8,b'e' as i8,b'a' as i8,b'd' as i8,b' ' as i8,b'c' as i8,b'o' as i8,b'u' as i8,b'n' as i8,b't' as i8,b's' as i8,b' ' as i8,b'l' as i8,b'i' as i8,b'n' as i8,b'e' as i8,0];
    const LONG_ERROR: &[i8] = &[b'T' as i8,b'o' as i8,b'o' as i8,b' ' as i8,b'l' as i8,b'o' as i8,b'n' as i8,b'g' as i8,b' ' as i8,b'c' as i8,b'o' as i8,b'u' as i8,b'n' as i8,b't' as i8,b's' as i8,b' ' as i8,b'l' as i8,b'i' as i8,b'n' as i8,b'e' as i8,0];
    const INTERPRET_ERROR: &[i8] = &[b'C' as i8,b'a' as i8,b'n' as i8,b'n' as i8,b'o' as i8,b't' as i8,b' ' as i8,b'i' as i8,b'n' as i8,b't' as i8,b'e' as i8,b'r' as i8,b'p' as i8,b'r' as i8,b'e' as i8,b't' as i8,b' ' as i8,b'c' as i8,b'o' as i8,b'u' as i8,b'n' as i8,b't' as i8,b's' as i8,b' ' as i8,b'l' as i8,b'i' as i8,b'n' as i8,b'e' as i8,b':' as i8,0];
    const RAM_ERROR: &[i8] = &[b'O' as i8,b'u' as i8,b't' as i8,b' ' as i8,b'o' as i8,b'f' as i8,b' ' as i8,b'R' as i8,b'A' as i8,b'M' as i8,0];

    let line = heap.allocate_model_storage(vec![0_i8; MOL_FMT_INPLINELEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let mut p = inchi_fgetsLf(heap, line, MOL_FMT_INPLINELEN as i32, inp_file)?;
        if p.is_null() {
            let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(READ_ERROR))?;
            return Ok(1);
        }
        remove_one_lf(heap, line)?;
        if heap.slice(line.as_const())?[MOL_FMT_MAXLINELEN as usize] != 0 {
            let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(LONG_ERROR))?;
        }
        let invalid = MolfileReadField(heap, MolfileFieldData::ShortIntoInt(&mut ctab.n_atoms), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
            || MolfileReadField(heap, MolfileFieldData::ShortIntoInt(&mut ctab.n_bonds), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
            || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
            || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
            || MolfileReadField(heap, MolfileFieldData::Char(&mut ctab.chiral_flag), 3, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)? < 0
            || MolfileReadField(heap, MolfileFieldData::Short(&mut ctab.n_stext_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
            || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
            || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
            || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
            || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
            || MolfileReadField(heap, MolfileFieldData::Short(&mut ctab.n_property_lines), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0;
        if invalid {
            let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(INTERPRET_ERROR))?;
            let _ = dotify_non_printable_chars(heap, line)?;
            let message = heap.slice(line.as_const())?;
            let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(message))?;
            return Ok(3);
        }
        let _ = MolfileReadField(heap, MolfileFieldData::String(&mut ctab.version_string), 6, i32::from(MOL_FMT_STRING_DATA), &mut p)?;
        let version_length = ctab.version_string.iter().position(|byte| *byte == 0).ok_or(SourceHeapError::MissingNulTerminator)?;
        if &ctab.version_string[..version_length] == [b'V' as i8,b'3' as i8,b'0' as i8,b'0' as i8,b'0' as i8] {
            ctab.v3000 = match inchi_calloc::<MOL_FMT_v3000>(heap, 1, std::mem::size_of::<MOL_FMT_v3000>() as u64) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => {
                    ctab.v3000 = SourceMutPointer::null();
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(RAM_ERROR))?;
                    return Ok(-1);
                }
                Err(error) => return Err(error),
            };
        } else {
            ctab.v3000 = SourceMutPointer::null();
        }
        let _ = MolFmtSgroups_Alloc(heap, Some(&mut ctab.sgroups), 1)?;
        Ok(0)
    })();
    let cleanup = heap.free(line);
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileReadAtomsBlock(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    mut err: i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:661 MolfileReadAtomsBlock
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap stack and object access add overhead.
    /*
int MolfileReadAtomsBlock( MOL_FMT_CTAB* ctab,
                           INCHI_IOSTREAM *inp_file,
                           int err,
                           char *pStrErr )
{
    char *p;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    int i;
    S_SHORT chg;
    static const S_SHORT charge_val[] = { 0, 3, 2, 1, 'R', -1, -2, -3 };

    /* djb-rwth: addressing coverity ID #499580 -- TREAT_ERR properly used in all cases */
    for (i = 0; i < ctab->n_atoms; i++)
    {
        p = inchi_fgetsLf( line, line_len, inp_file );

        if (!p)
        {
            if (!err)
            {
                TREAT_ERR( err, 2, "Cannot read atom block line" );
            }
            break;
        }

        remove_one_lf( line );


        if (line[MOL_FMT_MAXLINELEN])
        {
            TREAT_ERR( err, 0, "Too long atom block line" );
        }
        if (err)
        {
            if (!strcmp( line, SD_FMT_END_OF_DATA ))
            {
                err = -abs( err );
                break;
            }
            continue; /* bypass the rest of the Atom block */
        }

        if (NULL != ctab->coords)
        {
            mystrncpy( ctab->coords[i], p, 31 ); /* original coordinates */
        }

        if (NULL != ctab->atoms)
        {
            if (0 > MolfileReadField( &ctab->atoms[i].fx, 10, MOL_FMT_DOUBLE_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].fy, 10, MOL_FMT_DOUBLE_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].fz, 10, MOL_FMT_DOUBLE_DATA, &p )
                || 0 > MolfileReadField( NULL, /* undescribed in article*/    1, MOL_FMT_JUMP_TO_RIGHT, &p )
                || 0 == MolfileReadField( &ctab->atoms[i].symbol, 3, MOL_FMT_STRING_DATA, &p ) /* was sizeof(ctab->atoms[0].symbol)-1 */
                || 0 > MolfileReadField( &ctab->atoms[i].mass_difference, 2, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].charge, 3, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].stereo_parity, 3, MOL_FMT_CHAR_INT_DATA, &p )

#if ( MOL_FMT_QUERY == MOL_FMT_PRESENT )
                || 0 > MolfileReadField( &ctab->atoms[i].H_count_plus_1, 3, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].stereo_care, 3, MOL_FMT_CHAR_INT_DATA, &p )
#else
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

                || 0 > MolfileReadField( &ctab->atoms[i].valence, 3, MOL_FMT_CHAR_INT_DATA, &p ))
            {

                err = 4;
                TREAT_ERR( err, 4, "Cannot interpret atom block line:" );
                dotify_non_printable_chars( line );
                AddErrorMessage( pStrErr, line );

                if (!strcmp( line, SD_FMT_END_OF_DATA ))
                {
                    err = -abs( err );
                    break;
                }
                continue; /* can't interpret a first half of atom block line */
            }


            if (2 == strlen( ctab->atoms[i].symbol ) && isupper( UCINT ctab->atoms[i].symbol[1] ))
            {
                ctab->atoms[i].symbol[1] = (char) tolower( UCINT ctab->atoms[i].symbol[1] ); /* 5-4-99 DCh*/
            }

            if (( chg = (S_SHORT) ctab->atoms[i].charge ) < 0 || chg >= (int) ( sizeof( charge_val ) / sizeof( charge_val[0] ) ))
            {
                /* ctab->atoms[i].charge = 0; */ /* error; ignore for now */
                ctab->atoms[i].charge = (S_CHAR) ( 4 - chg ); /*  allow greater charges to accommodate NCI structures. 8-20-2002 */
                ctab->atoms[i].radical = 0;
            }
            else if ('R' == ( chg = charge_val[chg] ))
            {
                ctab->atoms[i].charge = 0;
                ctab->atoms[i].radical = RADICAL_DOUBLET;
            }
            else
            {
                ctab->atoms[i].charge = (S_CHAR) chg; /* actual charge value */
                ctab->atoms[i].radical = 0;
            }

            if (

#if ( MOL_FMT_CPSS == MOL_FMT_PRESENT )
                   0 > MolfileReadField( &ctab->atoms[i].H0_designator, 3, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].reaction_component_type, 3, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].reaction_component_num, 3, MOL_FMT_CHAR_INT_DATA, &p )
#else
                   0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

#if ( MOL_FMT_REACT == MOL_FMT_PRESENT )
                || 0 > MolfileReadField( &ctab->atoms[i].atom_atom_mapping_num, 3, MOL_FMT_SHORT_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->atoms[i].reaction_component_type, 3, MOL_FMT_CHAR_INT_DATA, &p )
#else
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

#if ( MOL_FMT_REACT == MOL_FMT_PRESENT || MOL_FMT_QUERY == MOL_FMT_PRESENT )
                || 0 > MolfileReadField( &ctab->atoms[i].exact_change_flag, 3, MOL_FMT_CHAR_INT_DATA, &p )
#else
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

            )
            {
                err = 5; /* can't interpret a second half of atom block line */

                TREAT_ERR( err, 5, "Cannot interpret atom block line:" );
                dotify_non_printable_chars( line );
                AddErrorMessage( pStrErr, line );

                if (!strcmp( line, SD_FMT_END_OF_DATA ))
                {
                    err = -abs( err );
                    break;
                }
                continue;
            }
        }
    }

/* err_fin: */

    return err;
}
    */
    // END INCHI C FUNCTION: MolfileReadAtomsBlock
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadAtomsBlock
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; MOL_FMT_QUERY, MOL_FMT_CPSS, and MOL_FMT_REACT are MOL_FMT_ABSENT.
    // INCHI✔️❌: The selected C locale uses ASCII isupper/tolower behavior and GCC signed narrowing wraps modulo 2^8.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadAtomsBlock

    const READ_ERROR: &[i8] = &[b'C' as i8,b'a' as i8,b'n' as i8,b'n' as i8,b'o' as i8,b't' as i8,b' ' as i8,b'r' as i8,b'e' as i8,b'a' as i8,b'd' as i8,b' ' as i8,b'a' as i8,b't' as i8,b'o' as i8,b'm' as i8,b' ' as i8,b'b' as i8,b'l' as i8,b'o' as i8,b'c' as i8,b'k' as i8,b' ' as i8,b'l' as i8,b'i' as i8,b'n' as i8,b'e' as i8,0];
    const LONG_ERROR: &[i8] = &[b'T' as i8,b'o' as i8,b'o' as i8,b' ' as i8,b'l' as i8,b'o' as i8,b'n' as i8,b'g' as i8,b' ' as i8,b'a' as i8,b't' as i8,b'o' as i8,b'm' as i8,b' ' as i8,b'b' as i8,b'l' as i8,b'o' as i8,b'c' as i8,b'k' as i8,b' ' as i8,b'l' as i8,b'i' as i8,b'n' as i8,b'e' as i8,0];
    const INTERPRET_ERROR: &[i8] = &[b'C' as i8,b'a' as i8,b'n' as i8,b'n' as i8,b'o' as i8,b't' as i8,b' ' as i8,b'i' as i8,b'n' as i8,b't' as i8,b'e' as i8,b'r' as i8,b'p' as i8,b'r' as i8,b'e' as i8,b't' as i8,b' ' as i8,b'a' as i8,b't' as i8,b'o' as i8,b'm' as i8,b' ' as i8,b'b' as i8,b'l' as i8,b'o' as i8,b'c' as i8,b'k' as i8,b' ' as i8,b'l' as i8,b'i' as i8,b'n' as i8,b'e' as i8,b':' as i8,0];
    let line = heap.allocate_model_storage(vec![0_i8; MOL_FMT_INPLINELEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let mut i = 0_i32;
        while i < ctab.n_atoms {
            let mut p = inchi_fgetsLf(heap, line, MOL_FMT_INPLINELEN as i32, inp_file.as_deref_mut())?;
            if p.is_null() {
                if err == 0 {
                    err = 2;
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(READ_ERROR))?;
                }
                break;
            }
            remove_one_lf(heap, line)?;
            if heap.slice(line.as_const())?[MOL_FMT_MAXLINELEN as usize] != 0 {
                let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(LONG_ERROR))?;
            }
            let is_end = {
                let bytes = heap.slice(line.as_const())?;
                bytes.get(..SD_FMT_END_OF_DATA.len()).is_some_and(|value| value.iter().map(|byte| *byte as u8).eq(SD_FMT_END_OF_DATA.iter().copied()))
            };
            if err != 0 {
                if is_end {
                    err = err.wrapping_abs().wrapping_neg();
                    break;
                }
                i = i.wrapping_add(1);
                continue;
            }

            if !ctab.coords.is_null() {
                let source = heap.slice(p.as_const())?.to_vec();
                let index = usize::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let coords = heap.slice_mut(ctab.coords)?;
                let target = coords.get_mut(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
                let _ = mystrncpy_slice(Some(target), Some(&source), 31)?;
            }

            if !ctab.atoms.is_null() {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let mut atom = heap.slice(ctab.atoms.as_const())?.get(index).cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
                let invalid_first = MolfileReadField(heap, MolfileFieldData::Double(&mut atom.fx), 10, i32::from(MOL_FMT_DOUBLE_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::Double(&mut atom.fy), 10, i32::from(MOL_FMT_DOUBLE_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::Double(&mut atom.fz), 10, i32::from(MOL_FMT_DOUBLE_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::None, 1, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::String(&mut atom.symbol), 3, i32::from(MOL_FMT_STRING_DATA), &mut p)? == 0
                    || MolfileReadField(heap, MolfileFieldData::Char(&mut atom.mass_difference), 2, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::Char(&mut atom.charge), 3, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::Char(&mut atom.stereo_parity), 3, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::Char(&mut atom.valence), 3, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)? < 0;
                if invalid_first {
                    heap.slice_mut(ctab.atoms)?[index] = atom;
                    err = 4;
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(INTERPRET_ERROR))?;
                    let _ = dotify_non_printable_chars(heap, line)?;
                    let source_line = heap.slice(line.as_const())?;
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(source_line))?;
                    if is_end {
                        err = err.wrapping_abs().wrapping_neg();
                        break;
                    }
                    i = i.wrapping_add(1);
                    continue;
                }

                let symbol_length = atom.symbol.iter().position(|byte| *byte == 0).ok_or(SourceHeapError::MissingNulTerminator)?;
                if symbol_length == 2 && (atom.symbol[1] as u8).is_ascii_uppercase() {
                    atom.symbol[1] = (atom.symbol[1] as u8).to_ascii_lowercase() as i8;
                }
                const CHARGE_VALUE: [i16; 8] = [0, 3, 2, 1, b'R' as i16, -1, -2, -3];
                let mut charge = i16::from(atom.charge);
                if charge < 0 || charge >= CHARGE_VALUE.len() as i16 {
                    atom.charge = 4_i16.wrapping_sub(charge) as i8;
                    atom.radical = 0;
                } else {
                    charge = CHARGE_VALUE[charge as usize];
                    if charge == b'R' as i16 {
                        atom.charge = 0;
                        atom.radical = RADICAL_DOUBLET as i8;
                    } else {
                        atom.charge = charge as i8;
                        atom.radical = 0;
                    }
                }

                let invalid_second = MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0;
                heap.slice_mut(ctab.atoms)?[index] = atom;
                if invalid_second {
                    err = 5;
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(INTERPRET_ERROR))?;
                    let _ = dotify_non_printable_chars(heap, line)?;
                    let source_line = heap.slice(line.as_const())?;
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(source_line))?;
                    if is_end {
                        err = err.wrapping_abs().wrapping_neg();
                        break;
                    }
                }
            }
            i = i.wrapping_add(1);
        }
        Ok(err)
    })();
    let cleanup = heap.free(line);
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileReadBondsBlock(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    mut err: i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:820 MolfileReadBondsBlock
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap stack and object access add overhead.
    /*
int MolfileReadBondsBlock( MOL_FMT_CTAB* ctab,
                           INCHI_IOSTREAM *inp_file,
                           int err,
                           char *pStrErr )
{
    char *p;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    int i;

#if 0
    if (NULL == ctab->bonds)
    {
        err = 1;
        goto err_fin;    /*internal error: memory has not been allocated for bonds structure*/
    }
#endif 

    /* djb-rwth: addressing coverity ID #499538 -- TREAT_ERR properly used in all cases */
    for (i = 0; i < ctab->n_bonds; i++)
    {
        p = inchi_fgetsLf( line, line_len, inp_file );

        if (!p)
        {
            if (!err)
            {
                TREAT_ERR( err, 2, "Cannot read bond block line" );
            }
            break;
        }

        remove_one_lf( line );

        if (line[MOL_FMT_MAXLINELEN])
        {
            err = err ? err : 3;             /* too long input file line */
        }

        if (err)
        {
            if (!strcmp( line, SD_FMT_END_OF_DATA ))
            {
                err = -abs( err );
                break;
            }
            continue;
        }

        if (ctab->bonds)
        {

            if (0 > MolfileReadField( &ctab->bonds[i].atnum1, 3, MOL_FMT_SHORT_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->bonds[i].atnum2, 3, MOL_FMT_SHORT_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->bonds[i].bond_type, 3, MOL_FMT_CHAR_INT_DATA, &p )
                || 0 > MolfileReadField( &ctab->bonds[i].bond_stereo, 3, MOL_FMT_CHAR_INT_DATA, &p )

#if ( MOL_FMT_QUERY == MOL_FMT_PRESENT )
                || 0 > MolfileReadField( &ctab->bonds[i].cBondTopology, 3, MOL_FMT_CHAR_INT_DATA, &p ) /* ring/chain */
#else
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

#if ( MOL_FMT_REACT == MOL_FMT_PRESENT )
                || 0 > MolfileReadField( &ctab->bonds[i].cReactingCenterStatus, 3, MOL_FMT_CHAR_INT_DATA, &p )
#else
                || 0 > MolfileReadField( NULL, 3, MOL_FMT_JUMP_TO_RIGHT, &p )
#endif

                )
            {

                if (!err)
                {
                    /* can't interpret bonds block line */
                    TREAT_ERR( err, 4, "Cannot interpret bond block line:" );
                    dotify_non_printable_chars( line );
                    AddErrorMessage( pStrErr, line );
                }
                if (!strcmp( line, SD_FMT_END_OF_DATA ))
                {
                    err = -abs( err );
                    break;
                }
            }
        }
    }

/* err_fin: */

    return err;
}
    */
    // END INCHI C FUNCTION: MolfileReadBondsBlock
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadBondsBlock
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; the #if 0 null-bond guard is inactive and MOL_FMT_QUERY/MOL_FMT_REACT are MOL_FMT_ABSENT.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadBondsBlock

    const READ_ERROR: &[i8] = &[b'C' as i8,b'a' as i8,b'n' as i8,b'n' as i8,b'o' as i8,b't' as i8,b' ' as i8,b'r' as i8,b'e' as i8,b'a' as i8,b'd' as i8,b' ' as i8,b'b' as i8,b'o' as i8,b'n' as i8,b'd' as i8,b' ' as i8,b'b' as i8,b'l' as i8,b'o' as i8,b'c' as i8,b'k' as i8,b' ' as i8,b'l' as i8,b'i' as i8,b'n' as i8,b'e' as i8,0];
    const INTERPRET_ERROR: &[i8] = &[b'C' as i8,b'a' as i8,b'n' as i8,b'n' as i8,b'o' as i8,b't' as i8,b' ' as i8,b'i' as i8,b'n' as i8,b't' as i8,b'e' as i8,b'r' as i8,b'p' as i8,b'r' as i8,b'e' as i8,b't' as i8,b' ' as i8,b'b' as i8,b'o' as i8,b'n' as i8,b'd' as i8,b' ' as i8,b'b' as i8,b'l' as i8,b'o' as i8,b'c' as i8,b'k' as i8,b' ' as i8,b'l' as i8,b'i' as i8,b'n' as i8,b'e' as i8,b':' as i8,0];
    let line = heap.allocate_model_storage(vec![0_i8; MOL_FMT_INPLINELEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let mut i = 0_i32;
        while i < ctab.n_bonds {
            let mut p = inchi_fgetsLf(heap, line, MOL_FMT_INPLINELEN as i32, inp_file.as_deref_mut())?;
            if p.is_null() {
                if err == 0 {
                    err = 2;
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(READ_ERROR))?;
                }
                break;
            }
            remove_one_lf(heap, line)?;
            if heap.slice(line.as_const())?[MOL_FMT_MAXLINELEN as usize] != 0 && err == 0 {
                err = 3;
            }
            let is_end = {
                let bytes = heap.slice(line.as_const())?;
                bytes.get(..SD_FMT_END_OF_DATA.len()).is_some_and(|value| value.iter().map(|byte| *byte as u8).eq(SD_FMT_END_OF_DATA.iter().copied()))
            };
            if err != 0 {
                if is_end {
                    err = err.wrapping_abs().wrapping_neg();
                    break;
                }
                i = i.wrapping_add(1);
                continue;
            }
            if !ctab.bonds.is_null() {
                let index = usize::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let mut bond = heap.slice(ctab.bonds.as_const())?.get(index).cloned().ok_or(SourceHeapError::PointerOutOfBounds)?;
                let invalid = MolfileReadField(heap, MolfileFieldData::Short(&mut bond.atnum1), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::Short(&mut bond.atnum2), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::Char(&mut bond.bond_type), 3, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::Char(&mut bond.bond_stereo), 3, i32::from(MOL_FMT_CHAR_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::None, 3, i32::from(MOL_FMT_JUMP_TO_RIGHT), &mut p)? < 0;
                heap.slice_mut(ctab.bonds)?[index] = bond;
                if invalid {
                    err = 4;
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(INTERPRET_ERROR))?;
                    let _ = dotify_non_printable_chars(heap, line)?;
                    let source_line = heap.slice(line.as_const())?;
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(source_line))?;
                    if is_end {
                        err = err.wrapping_abs().wrapping_neg();
                        break;
                    }
                }
            }
            i = i.wrapping_add(1);
        }
        Ok(err)
    })();
    let cleanup = heap.free(line);
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileReadSTextBlock(
    heap: &mut SourceHeap,
    ctab: &MOL_FMT_CTAB,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    mut err: i32,
    mut p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:917 MolfileReadSTextBlock
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap stack-buffer modeling adds overhead.
    /*
int MolfileReadSTextBlock( MOL_FMT_CTAB* ctab,
                           INCHI_IOSTREAM *inp_file,
                           int err,
                           char *pStrErr )
{
    /* just pass by all stext enties without attemp to interpret them */
    char *p;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    S_SHORT i;

    for (i = 0; i < 2 * ctab->n_stext_entries; i++)
    {
        p = inchi_fgetsLf( line, line_len, inp_file );
        if (!p)
        {
            if (!err)
            {
                TREAT_ERR_AND_FIN( err, 2, err_fin, "Cannot read STEXT block line" ); /* djb-rwth: addressing coverity ID #499517 -- TREAT_ERR_AND_FIN properly used */
            }
            break;
            /* can't read the input file line */
        }

        /*
        remove_one_lf( line );
        if ( line[MOL_FMT_MAXLINELEN] ){
            TREAT_ERR( err, 2, "Warning: Too long STEXT block line");
            too long input file line
        }
        */
    }

err_fin:

    return err;
}
    */
    // END INCHI C FUNCTION: MolfileReadSTextBlock
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadSTextBlock
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; line trimming and long-line checking are inside an inactive block comment.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadSTextBlock

    const READ_ERROR: &[i8] = &[b'C' as i8,b'a' as i8,b'n' as i8,b'n' as i8,b'o' as i8,b't' as i8,b' ' as i8,b'r' as i8,b'e' as i8,b'a' as i8,b'd' as i8,b' ' as i8,b'S' as i8,b'T' as i8,b'E' as i8,b'X' as i8,b'T' as i8,b' ' as i8,b'b' as i8,b'l' as i8,b'o' as i8,b'c' as i8,b'k' as i8,b' ' as i8,b'l' as i8,b'i' as i8,b'n' as i8,b'e' as i8,0];
    let line = heap.allocate_model_storage(vec![0_i8; MOL_FMT_INPLINELEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let limit = 2_i32.wrapping_mul(i32::from(ctab.n_stext_entries));
        let mut i = 0_i16;
        while i32::from(i) < limit {
            let p = inchi_fgetsLf(heap, line, MOL_FMT_INPLINELEN as i32, inp_file.as_deref_mut())?;
            if p.is_null() {
                if err == 0 {
                    err = 2;
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(READ_ERROR))?;
                }
                break;
            }
            i = i.wrapping_add(1);
        }
        Ok(err)
    })();
    let cleanup = heap.free(line);
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MolfileReadPropBlock(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    p_hdr: &mut MOL_FMT_HEADER_BLOCK,
    mut inp_file: Option<&mut INCHI_IOSTREAM>,
    treat_polymers: i32,
    mut err: i32,
    mut p_str_err: Option<&mut [i8]>,
    b_no_warnings: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:960 MolfileReadPropBlock
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap stack storage and typed field adapters add overhead.
    /*
int MolfileReadPropBlock( MOL_FMT_CTAB* ctab,
                          MOL_FMT_HEADER_BLOCK *pHdr,
                          INCHI_IOSTREAM *inp_file,
                          int treat_polymers,
                          int err,
                          char *pStrErr,
                          int bNoWarnings )
{
    enum { MULTI_LINE_MODE_NO_MODE, MULTI_LINE_MODE_ISIS_ALIAS };
    char *p;
    char line[MOL_FMT_INPLINELEN];
    const int line_len = sizeof( line );
    int   nMultiLineMode = MULTI_LINE_MODE_NO_MODE, nAtomNumber = 0;
    S_SHORT i, j;
    char  charM[2];
    char  szBlank[3];
    char  szType[4];
    S_SHORT  skip_lines = 0;
    S_SHORT  num_entries;
    S_SHORT  num_atoms = ctab->n_atoms;

    int  charge_encountered = 0;
    int  radical_encountered = 0;
    int  isotope_encountered = 0;
    int     polymer_occurred = 0;

    szType[0] = '\0'; /* djb-rwth: adding zero termination */
    int debug_polymers = 0;
#if ( DEBUG_POLYMERS == 1 )
    debug_polymers = 1;
#elif  ( DEBUG_POLYMERS == 2 )
    debug_polymers = 2;
#endif


    for (i = 0; ctab->version_string[0] ? 1 : ( i < ctab->n_property_lines ); i++)
    {
        /* the last line should be M END */

        /* ctab->version_string[0] == 0:
              exactly ctab->n_property_lines lines including M END */
        /* ctab->version_string[0] != 0:
              read until M END line was encountered */

        p = inchi_fgetsLf( line, line_len, inp_file );
        /* djb-rwth: addressing coverity ID #499577 -- TREAT_ERR properly used in all cases */
        if (!p)
        {
            if (!err)
            {
                TREAT_ERR( err, 2, "Cannot read properties block line" );
            }
            goto err_fin;
        }

        remove_one_lf( line );

        if (line[MOL_FMT_MAXLINELEN])
        {
            TREAT_ERR( err, 3, "Too long properties block line" );
            continue;
        }

        if (skip_lines > 0)
        {
            skip_lines--;
            continue;
        }

        /* alias */
        if (nMultiLineMode == MULTI_LINE_MODE_ISIS_ALIAS && nAtomNumber)
        {
            int  len;

            nMultiLineMode = MULTI_LINE_MODE_NO_MODE;
            if (0 >= ( len = normalize_string( p ) ))
            {
                nAtomNumber = 0;
                continue;
            }

            if (0 < len && len < (int) ( sizeof( ctab->atoms->symbol ) ))
            {
                int  nCharge, nRad;

                MOL_FMT_ATOM*  atom = ctab->atoms + nAtomNumber - 1;
                /* ctab->atoms[nAtomNumber-1].atom_aliased_flag = 1; */
                /*  extract radicals & charges */

                extract_charges_and_radicals( p, &nRad, &nCharge );

                /*  Aliased atom cannot have charge, radical & mass difference */
                /*  in the atom table or "M  CHG", "M  RAD", "M  ISO" */
                /* if ( nCharge ) */
                atom->charge = (S_CHAR) nCharge;
                /* if ( nRad ) */
                atom->radical = (char) nRad;

                if (1 == len && 'D' == p[0])
                {
                    /*  H isotope */
                    p[0] = 'H';
                    atom->mass_difference = 1;
                }
                else
                {
                    if (1 == len && 'T' == p[0])
                    {
                        /*  H isotope */

                        p[0] = 'H';
                        atom->mass_difference = 2;
                    }
                    else
                    {
                        atom->mass_difference = 0;
                    }
                }
                if (strlen( p ) < sizeof( ctab->atoms[0].symbol ))
                {
                    strcpy(atom->symbol, p);
                }
                else
                {
                    strcpy(atom->symbol, "???");
                }
                atom->atom_aliased_flag++;
            }
            /* else if( 0 < len  )
            {
               ^^^ Just too long alias name.
                     For consistency with parsing {H,D,T}-containing alias names, this
                     would result in issuing error rather than ignoring... However,
                     for compatibility reasons, the 'ignore' behavior remained intact.
            }
            */

            skip_lines = 0;
            nAtomNumber = 0;
            continue;
        }

        if (1 != MolfileReadField( charM, sizeof( charM ) - 1, MOL_FMT_STRING_DATA, &p )
            || 0 != MolfileReadField( szBlank, sizeof( szBlank ) - 1, MOL_FMT_STRING_DATA, &p ) /* must contain 0 bytes */
            || 0 >= MolfileReadField( szType, sizeof( szType ) - 1, MOL_FMT_STRING_DATA, &p ) /* must contain 3 bytes */
        )
        {
            if (!strcmp( line, SD_FMT_END_OF_DATA ))
            {
                err = err ? -abs( err ) : -4;
                break;
            }
            continue;  /* ignore because cannot recognize */
        }

        if (charM[0] == 'V')
        {
            skip_lines = 0;   /* ISIS/Desktop Atom Value: one-line property */
            continue;
        }

        if (charM[0] == 'G')
        {
            skip_lines = 1;   /* ISIS/Desktop Group abbreviation: two-line property */
            continue;
        }

        if (charM[0] == 'A')
        {
            if (NULL != ctab->atoms &&
                 0 < ( nAtomNumber = (int) strtol( szType, NULL, 10 ) ) &&
                 nAtomNumber <= ctab->n_atoms)
            {
                /* Atom Alias [ISIS/Desktop] two-line property */
                nMultiLineMode = MULTI_LINE_MODE_ISIS_ALIAS;
                continue;
            }
            else
            {
                nAtomNumber = 0;
                skip_lines = 1;
                continue;
            }
        }

        if (charM[0] == 'S' && !strcmp( szType, "SKP" ))
        {  /* skip lines */
            if (0 >= MolfileReadField( &skip_lines, 3, MOL_FMT_SHORT_INT_DATA, &p ))
            {
                skip_lines = 0;
            }
            continue;
        }

        if (charM[0] != 'M')
        {
            /* cannot recognize a line */
            continue;
        }

        if (!strcmp( szType, "REG" ))
        {
            int len;
            p = p + strspn( p, " " );
            len = strcspn( p, " " );
            len = inchi_min( len, MOL_FMT_MAX_VALUE_LEN );
            MolfileReadField( &pHdr->internal_regno, len, MOL_FMT_LONG_INT_DATA, &p );
            continue;
        }

        if (!strcmp( szType, "END" ))
        {
            if (ctab->version_string[0])
            {
                break;  /* end of property lines */
            }
            continue;
        }

        if (NULL == ctab->atoms)
        {
            continue; /* ignore because the user requested to bypass all this stuff */
        }

        /* Charge: Generic */
        if (!strcmp( szType, "CHG" ) &&
             0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
             1 <= num_entries && num_entries <= 8)
        {

            S_SHORT atoms[8];
            S_SHORT charges[8];

            if (!charge_encountered && !radical_encountered)
            {
                /* first charge or radical record clears all Atom Block */
                /* entered charge and radical data to zeroes            */
                charge_encountered = -1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (0 > MolfileReadField( &atoms[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     0 > MolfileReadField( &charges[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     atoms[j] <= 0 || atoms[j] > num_atoms ||
                     charges[j] < -15 || charges[j]  > 15)
                {
                    goto charge_error;
                }
            }
            if (charge_encountered == -1)
            {
                for (j = 0; j < num_atoms; j++)
                {
                    if (!ctab->atoms[j].atom_aliased_flag) /* do not clear aliased atoms.*/
                    {
                        ctab->atoms[j].charge = ctab->atoms[j].radical = '\0';
                    }
                }
                charge_encountered = 1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (!ctab->atoms[atoms[j] - 1].atom_aliased_flag) /* do not change aliased atoms.*/
                {
                    ctab->atoms[atoms[j] - 1].charge = (S_CHAR) charges[j];
                }
            }
            continue;
        charge_error:
            TREAT_ERR( err, 0, "Charge not recognized:" );
            dotify_non_printable_chars( line );
            AddErrorMessage( pStrErr, line );
            continue; /* ignore for now */
        }

        /* Radical: Generic */
        if (!strcmp( szType, "RAD" ) &&
             0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
             1 <= num_entries && num_entries <= 8)
        {

            S_SHORT atoms[8];
            S_SHORT radicals[8];

            if (!charge_encountered && !radical_encountered)
            {
                /* first charge or radical record clears all Atom Block */
                /* entered charge and radical data to zeroes            */
                radical_encountered = -1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (0 > MolfileReadField( &atoms[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     0 > MolfileReadField( &radicals[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     atoms[j] <= 0 || atoms[j] > num_atoms ||
                     radicals[j] < 0 || radicals[j]  > 3)
                {
                    goto radical_error;
                }
            }
            if (radical_encountered == -1)
            {
                for (j = 0; j < num_atoms; j++)
                {
                    if (!ctab->atoms[j].atom_aliased_flag)  /* do not clear aliased atoms. 5-3-99 DCh */
                        ctab->atoms[j].charge = ctab->atoms[j].radical = '\0';
                }
                radical_encountered = 1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (!ctab->atoms[atoms[j] - 1].atom_aliased_flag)
                {
                    /* do not change aliased atoms. 5-3-99 DCh */
                    ctab->atoms[atoms[j] - 1].radical = (S_CHAR) radicals[j];
                }
            }
            continue;
        radical_error:
            TREAT_ERR( err, 0, "Radical not recognized:" );
            dotify_non_printable_chars( line );
            AddErrorMessage( pStrErr, line );
            continue; /* ignore error for now */
        }

        /* Isotope: Generic */
        if (!strcmp( szType, "ISO" ) &&
             0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
             1 <= num_entries && num_entries <= 8)
        {

            S_SHORT atoms[8];
            S_SHORT iso_mass[8]; /*  contains istotope mass number, not difference. 7-14-00 DCh. */

            if (!isotope_encountered)
            {
                /* first charge or radical record clears all Atom Block */
                /* entered charge and radical data to zeroes            */
                isotope_encountered = -1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (0 > MolfileReadField( &atoms[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     0 > MolfileReadField( &iso_mass[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                     atoms[j] <= 0 || atoms[j] > num_atoms
                     /*|| iso_mass[j] < -18 || iso_mass[j]  > 12*/)
                {
                    /* goto isotope_error; */
                    atoms[j] = -1; /*  flag error */
                    TREAT_ERR( err, 0, "Isotopic data not recognized:" );
                    dotify_non_printable_chars( line );
                    AddErrorMessage( pStrErr, line );
                    continue; /* ignore isotopic error for now */
                }
            }

            if (isotope_encountered == -1)
            {
                for (j = 0; j < num_atoms; j++)
                {
                    /*if ( !ctab->atoms[j].atom_aliased_flag )*/  /* clear even aliased atoms */
                    ctab->atoms[j].mass_difference = 0;
                }
                isotope_encountered = 1;
            }
            for (j = 0; j < num_entries; j++)
            {
                if (atoms[j] <= 0)
                {
                    continue; /* ignore isotopic error for now */
                }

                if (1 /* !ctab->atoms[atoms[j]-1].atom_aliased_flag */)
                {
                    char *at = ctab->atoms[atoms[j] - 1].symbol;
                    if (at[1] || (at[0] != 'D' && at[0] != 'T')) /* djb-rwth: addressing LLVM warning */
                    {  /*  D & T cannot have ISO */
                        /*  need atomic weight to calculate isotope difference. 7-14-00 DCh. */

                        int  atw, atw_diff;
                        /*
                        NB: According to CTFile specification, difference should be in
                        [-18; +12] range, not in [-19; +19] as is checked below. */
                        if (( atw = get_atomic_mass( at ) ) && abs( atw_diff = (int) iso_mass[j] - atw ) < 20)
                        {
                            ctab->atoms[atoms[j] - 1].mass_difference = (char) ( atw_diff ? atw_diff : ZERO_ATW_DIFF );
                        }
                    }
                }
            }
            continue;
        }

        /* Sgroup, polymeric */
        if ( !strcmp( szType, "STY" ) ||
             !strcmp( szType, "SST" ) ||
             !strcmp( szType, "SLB" ) ||
             !strcmp( szType, "SCN" ) ||
             !strcmp( szType, "SAL" ) ||
             !strcmp( szType, "SBL" ) ||
             !strcmp( szType, "SDI" ) ||
             !strcmp( szType, "SMT" ) ||
             !strcmp( szType, "SBT" ))
        {
            int result;
            if (!treat_polymers)
            {
                polymer_occurred = 1;
                continue;
            }
            result = MolfileReadSgroupOfPolymer( ctab, pHdr, inp_file, line, szType, p, err, pStrErr );
            if (result != 0)
            {
                TREAT_ERR( err, result, "Could not interpret Molfile polymer data:" );
                dotify_non_printable_chars( line );
                AddErrorMessage( pStrErr, line );
                continue;
            }
        }
    }

err_fin:
    if (!treat_polymers && polymer_occurred)
    {
        /* for compatibility reasons, inchi-1 by default
        ignores polymer related lines (as v. 1.04 did)    */
        if (!bNoWarnings)
        {
            WarningMessage( pStrErr, "Ignore polymer data" );
        }
    }

    if (( ctab->sgroups.used > 0 ) && ( debug_polymers > 1 ))
    {
        ITRACE_( "\n* THE FOLLOWING %-d POLYMER SGROUP(S) WERE RECOGNISED *\n", ctab->sgroups.used );
        for (i = 0; i < ctab->sgroups.used; i++)
        {
            char *sty[] = { "NON", "SRU", "MON", "COP", "MOD", "XL", "MER" };
            char *sst[] = { "NON", "ALT", "RAN", "BLK" };
            char *con[] = { "NON", "HT", "HH", "EU" };

            ITRACE_( "\n* GROUP %-d\n", i );
            ITRACE_( "* \tindex=%-d\n", ctab->sgroups.group[i]->id );
            ITRACE_( "* \ttype=%-d %-s\n", ctab->sgroups.group[i]->type, sty[ctab->sgroups.group[i]->type] );
            if (ctab->sgroups.group[i]->subtype > -1)
                ITRACE_( "* \tsubtype=%-d %-s\n", ctab->sgroups.group[i]->subtype, sst[ctab->sgroups.group[i]->subtype] );
            if (ctab->sgroups.group[i]->conn)
                ITRACE_( "* \tconnection_type=%-d %-s\n", ctab->sgroups.group[i]->conn,
                                                                 con[ctab->sgroups.group[i]->conn] );
            ITRACE_( "* \tlabel=%-d\n", ctab->sgroups.group[i]->label );
            ITRACE_( "* \t%-d atoms:\t", ctab->sgroups.group[i]->alist.used );
            IntArray_DebugPrint( &( ctab->sgroups.group[i]->alist ) );
            ITRACE_( "* \t%-d bonds:\t", ctab->sgroups.group[i]->blist.used );
            IntArray_DebugPrint( &( ctab->sgroups.group[i]->blist ) );
            ITRACE_( "\n" );
        }
    }

    return err;
}
    */
    // END INCHI C FUNCTION: MolfileReadPropBlock
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadPropBlock
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; DEBUG_POLYMERS is zero, so both debug assignment and ITRACE branches are inactive.
    // INCHI✔️❌: WarningMessage is the active AddErrorMessage macro alias; ZERO_ATW_DIFF is 127.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadPropBlock

    const MULTI_LINE_MODE_NO_MODE: i32 = 0;
    const MULTI_LINE_MODE_ISIS_ALIAS: i32 = 1;
    let line = heap.allocate_model_storage(vec![0_i8; MOL_FMT_INPLINELEN as usize])?;
    let result = (|| -> Result<i32, SourceHeapError> {
        let mut multi_line_mode = MULTI_LINE_MODE_NO_MODE;
        let mut atom_number = 0_i32;
        let mut skip_lines = 0_i16;
        let num_atoms = ctab.n_atoms;
        let mut charge_encountered = 0_i32;
        let mut radical_encountered = 0_i32;
        let mut isotope_encountered = 0_i32;
        let mut polymer_occurred = 0_i32;
        let mut i = 0_i16;

        while ctab.version_string[0] != 0 || i < ctab.n_property_lines {
            let mut p = inchi_fgetsLf(heap, line, MOL_FMT_INPLINELEN as i32, inp_file.as_deref_mut())?;
            if p.is_null() {
                if err == 0 {
                    mol_fmt1_treat_error(&mut err, 2, p_str_err.as_deref_mut(), b"Cannot read properties block line")?;
                }
                break;
            }
            remove_one_lf(heap, line)?;
            if heap.slice(line.as_const())?[MOL_FMT_MAXLINELEN as usize] != 0 {
                mol_fmt1_treat_error(&mut err, 3, p_str_err.as_deref_mut(), b"Too long properties block line")?;
                i = i.wrapping_add(1);
                continue;
            }
            if skip_lines > 0 {
                skip_lines = skip_lines.wrapping_sub(1);
                i = i.wrapping_add(1);
                continue;
            }

            if multi_line_mode == MULTI_LINE_MODE_ISIS_ALIAS && atom_number != 0 {
                multi_line_mode = MULTI_LINE_MODE_NO_MODE;
                let len = normalize_string(heap, p)?;
                if len <= 0 {
                    atom_number = 0;
                    i = i.wrapping_add(1);
                    continue;
                }
                if len < 6 {
                    let mut charge = 0_i32;
                    let mut radical = 0_i32;
                    extract_charges_and_radicals(Some(heap.slice_mut(p)?), Some(&mut radical), Some(&mut charge))?;
                    let symbol = {
                        let value = heap.slice(p.as_const())?;
                        let length = value.iter().position(|byte| *byte == 0).ok_or(SourceHeapError::MissingNulTerminator)?;
                        value[..=length].to_vec()
                    };
                    let atom_index = usize::try_from(atom_number.wrapping_sub(1)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom = heap.slice_mut(ctab.atoms)?.get_mut(atom_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
                    atom.charge = charge as i8;
                    atom.radical = radical as i8;
                    if len == 1 && symbol[0] == b'D' as i8 {
                        atom.symbol[0] = b'H' as i8;
                        atom.symbol[1] = 0;
                        atom.mass_difference = 1;
                    } else if len == 1 && symbol[0] == b'T' as i8 {
                        atom.symbol[0] = b'H' as i8;
                        atom.symbol[1] = 0;
                        atom.mass_difference = 2;
                    } else {
                        atom.mass_difference = 0;
                        if symbol.len() < atom.symbol.len() {
                            atom.symbol.fill(0);
                            atom.symbol[..symbol.len()].copy_from_slice(&symbol);
                        } else {
                            atom.symbol = [b'?' as i8, b'?' as i8, b'?' as i8, 0, 0, 0];
                        }
                    }
                    atom.atom_aliased_flag = atom.atom_aliased_flag.wrapping_add(1);
                }
                skip_lines = 0;
                atom_number = 0;
                i = i.wrapping_add(1);
                continue;
            }

            let mut char_m = [0_i8; 2];
            let mut blank = [0_i8; 3];
            let mut type_ = [0_i8; 4];
            let recognized = MolfileReadField(heap, MolfileFieldData::String(&mut char_m), 1, i32::from(MOL_FMT_STRING_DATA), &mut p)? == 1
                && MolfileReadField(heap, MolfileFieldData::String(&mut blank), 2, i32::from(MOL_FMT_STRING_DATA), &mut p)? == 0
                && MolfileReadField(heap, MolfileFieldData::String(&mut type_), 3, i32::from(MOL_FMT_STRING_DATA), &mut p)? > 0;
            if !recognized {
                if mol_fmt1_c_string_eq(heap.slice(line.as_const())?, &SD_FMT_END_OF_DATA[..SD_FMT_END_OF_DATA.len() - 1]) {
                    err = if err != 0 { err.wrapping_abs().wrapping_neg() } else { -4 };
                    break;
                }
                i = i.wrapping_add(1);
                continue;
            }
            match char_m[0] as u8 {
                b'V' => {
                    skip_lines = 0;
                    i = i.wrapping_add(1);
                    continue;
                }
                b'G' => {
                    skip_lines = 1;
                    i = i.wrapping_add(1);
                    continue;
                }
                b'A' => {
                    let temporary = heap.allocate_model_storage(type_.to_vec())?;
                    atom_number = inchi_strtol(heap, temporary.as_const(), None, 10)? as i32;
                    heap.free(temporary)?;
                    if !ctab.atoms.is_null() && atom_number > 0 && atom_number <= ctab.n_atoms {
                        multi_line_mode = MULTI_LINE_MODE_ISIS_ALIAS;
                    } else {
                        atom_number = 0;
                        skip_lines = 1;
                    }
                    i = i.wrapping_add(1);
                    continue;
                }
                b'S' if mol_fmt1_c_string_eq(&type_, b"SKP") => {
                    if MolfileReadField(heap, MolfileFieldData::Short(&mut skip_lines), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? <= 0 {
                        skip_lines = 0;
                    }
                    i = i.wrapping_add(1);
                    continue;
                }
                b'M' => {}
                _ => {
                    i = i.wrapping_add(1);
                    continue;
                }
            }

            if mol_fmt1_c_string_eq(&type_, b"REG") {
                let remaining = heap.slice(p.as_const())?;
                let nul = remaining.iter().position(|byte| *byte == 0).ok_or(SourceHeapError::MissingNulTerminator)?;
                let leading = remaining[..nul].iter().take_while(|byte| **byte == b' ' as i8).count();
                p = p.offset(i64::try_from(leading).map_err(|_| SourceHeapError::PointerOffsetOverflow)?)?;
                let remaining = heap.slice(p.as_const())?;
                let nul = remaining.iter().position(|byte| *byte == 0).ok_or(SourceHeapError::MissingNulTerminator)?;
                let width = remaining[..nul].iter().position(|byte| *byte == b' ' as i8).unwrap_or(nul).min(MOL_FMT_MAX_VALUE_LEN as usize);
                let _ = MolfileReadField(heap, MolfileFieldData::Long(&mut p_hdr.internal_regno), i32::try_from(width).map_err(|_| SourceHeapError::SourceIntegerOverflow)?, i32::from(MOL_FMT_LONG_INT_DATA), &mut p)?;
                i = i.wrapping_add(1);
                continue;
            }
            if mol_fmt1_c_string_eq(&type_, b"END") {
                if ctab.version_string[0] != 0 {
                    break;
                }
                i = i.wrapping_add(1);
                continue;
            }
            if ctab.atoms.is_null() {
                i = i.wrapping_add(1);
                continue;
            }

            let mut num_entries = 0_i16;
            if mol_fmt1_c_string_eq(&type_, b"CHG")
                && MolfileReadField(heap, MolfileFieldData::Short(&mut num_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? > 0
                && (1..=8).contains(&num_entries)
            {
                let mut atom_numbers = [0_i16; 8];
                let mut charges = [0_i16; 8];
                if charge_encountered == 0 && radical_encountered == 0 { charge_encountered = -1; }
                let mut invalid = false;
                for j in 0..usize::try_from(num_entries).map_err(|_| SourceHeapError::SourceIntegerOverflow)? {
                    if MolfileReadField(heap, MolfileFieldData::Short(&mut atom_numbers[j]), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                        || MolfileReadField(heap, MolfileFieldData::Short(&mut charges[j]), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                        || atom_numbers[j] <= 0 || i32::from(atom_numbers[j]) > num_atoms || charges[j] < -15 || charges[j] > 15 {
                        invalid = true;
                        break;
                    }
                }
                if invalid {
                    mol_fmt1_treat_error(&mut err, 0, p_str_err.as_deref_mut(), b"Charge not recognized:")?;
                    let _ = dotify_non_printable_chars(heap, line)?;
                    let source_line = heap.slice(line.as_const())?.to_vec();
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
                    i = i.wrapping_add(1);
                    continue;
                }
                let atoms = heap.slice_mut(ctab.atoms)?;
                if charge_encountered == -1 {
                    for atom in atoms.iter_mut().take(usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?) {
                        if atom.atom_aliased_flag == 0 { atom.charge = 0; atom.radical = 0; }
                    }
                    charge_encountered = 1;
                }
                for j in 0..usize::try_from(num_entries).map_err(|_| SourceHeapError::SourceIntegerOverflow)? {
                    let index = usize::try_from(atom_numbers[j] - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    if atoms[index].atom_aliased_flag == 0 { atoms[index].charge = charges[j] as i8; }
                }
                i = i.wrapping_add(1);
                continue;
            }

            if mol_fmt1_c_string_eq(&type_, b"RAD")
                && MolfileReadField(heap, MolfileFieldData::Short(&mut num_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? > 0
                && (1..=8).contains(&num_entries)
            {
                let mut atom_numbers = [0_i16; 8];
                let mut radicals = [0_i16; 8];
                if charge_encountered == 0 && radical_encountered == 0 { radical_encountered = -1; }
                let mut invalid = false;
                for j in 0..usize::try_from(num_entries).map_err(|_| SourceHeapError::SourceIntegerOverflow)? {
                    if MolfileReadField(heap, MolfileFieldData::Short(&mut atom_numbers[j]), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                        || MolfileReadField(heap, MolfileFieldData::Short(&mut radicals[j]), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                        || atom_numbers[j] <= 0 || i32::from(atom_numbers[j]) > num_atoms || radicals[j] < 0 || radicals[j] > 3 {
                        invalid = true;
                        break;
                    }
                }
                if invalid {
                    mol_fmt1_treat_error(&mut err, 0, p_str_err.as_deref_mut(), b"Radical not recognized:")?;
                    let _ = dotify_non_printable_chars(heap, line)?;
                    let source_line = heap.slice(line.as_const())?.to_vec();
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
                    i = i.wrapping_add(1);
                    continue;
                }
                let atoms = heap.slice_mut(ctab.atoms)?;
                if radical_encountered == -1 {
                    for atom in atoms.iter_mut().take(usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?) {
                        if atom.atom_aliased_flag == 0 { atom.charge = 0; atom.radical = 0; }
                    }
                    radical_encountered = 1;
                }
                for j in 0..usize::try_from(num_entries).map_err(|_| SourceHeapError::SourceIntegerOverflow)? {
                    let index = usize::try_from(atom_numbers[j] - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    if atoms[index].atom_aliased_flag == 0 { atoms[index].radical = radicals[j] as i8; }
                }
                i = i.wrapping_add(1);
                continue;
            }

            if mol_fmt1_c_string_eq(&type_, b"ISO")
                && MolfileReadField(heap, MolfileFieldData::Short(&mut num_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? > 0
                && (1..=8).contains(&num_entries)
            {
                let mut atom_numbers = [0_i16; 8];
                let mut isotope_mass = [0_i16; 8];
                if isotope_encountered == 0 { isotope_encountered = -1; }
                for j in 0..usize::try_from(num_entries).map_err(|_| SourceHeapError::SourceIntegerOverflow)? {
                    if MolfileReadField(heap, MolfileFieldData::Short(&mut atom_numbers[j]), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                        || MolfileReadField(heap, MolfileFieldData::Short(&mut isotope_mass[j]), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                        || atom_numbers[j] <= 0 || i32::from(atom_numbers[j]) > num_atoms {
                        atom_numbers[j] = -1;
                        mol_fmt1_treat_error(&mut err, 0, p_str_err.as_deref_mut(), b"Isotopic data not recognized:")?;
                        let _ = dotify_non_printable_chars(heap, line)?;
                        let source_line = heap.slice(line.as_const())?.to_vec();
                        let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
                    }
                }
                if isotope_encountered == -1 {
                    for atom in heap.slice_mut(ctab.atoms)?.iter_mut().take(usize::try_from(num_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?) { atom.mass_difference = 0; }
                    isotope_encountered = 1;
                }
                for j in 0..usize::try_from(num_entries).map_err(|_| SourceHeapError::SourceIntegerOverflow)? {
                    if atom_numbers[j] <= 0 { continue; }
                    let index = usize::try_from(atom_numbers[j] - 1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let symbol = heap.slice(ctab.atoms.as_const())?[index].symbol;
                    if symbol[1] != 0 || (symbol[0] != b'D' as i8 && symbol[0] != b'T' as i8) {
                        let atw = get_atomic_mass(Some(&symbol))?;
                        let difference = i32::from(isotope_mass[j]).wrapping_sub(atw);
                        if atw != 0 && difference.wrapping_abs() < 20 {
                            heap.slice_mut(ctab.atoms)?[index].mass_difference = if difference != 0 { difference as i8 } else { ZERO_ATW_DIFF as i8 };
                        }
                    }
                }
                i = i.wrapping_add(1);
                continue;
            }

            if [b"STY", b"SST", b"SLB", b"SCN", b"SAL", b"SBL", b"SDI", b"SMT", b"SBT"]
                .iter().any(|candidate| mol_fmt1_c_string_eq(&type_, *candidate))
            {
                if treat_polymers == 0 {
                    polymer_occurred = 1;
                    i = i.wrapping_add(1);
                    continue;
                }
                let polymer_result = MolfileReadSgroupOfPolymer(heap, ctab, Some(p_hdr), inp_file.as_deref_mut(), line, &type_, p, err, p_str_err.as_deref_mut())?;
                if polymer_result != 0 {
                    mol_fmt1_treat_error(&mut err, polymer_result, p_str_err.as_deref_mut(), b"Could not interpret Molfile polymer data:")?;
                    let _ = dotify_non_printable_chars(heap, line)?;
                    let source_line = heap.slice(line.as_const())?.to_vec();
                    let _ = AddErrorMessage(p_str_err.as_deref_mut(), Some(&source_line))?;
                }
            }
            i = i.wrapping_add(1);
        }

        if treat_polymers == 0 && polymer_occurred != 0 && b_no_warnings == 0 {
            let _ = mol_fmt1_add_ascii_message(p_str_err.as_deref_mut(), b"Ignore polymer data")?;
        }
        Ok(err)
    })();
    let cleanup = heap.free(line);
    match (result, cleanup) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(value), Ok(())) => Ok(value),
    }
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn MolfileReadSgroupOfPolymer(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    _p_hdr: Option<&mut MOL_FMT_HEADER_BLOCK>,
    _inp_file: Option<&mut INCHI_IOSTREAM>,
    _line: SourceMutPointer<i8>,
    sz_type: &[i8],
    mut p: SourceMutPointer<i8>,
    mut err: i32,
    _p_str_err: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:1425 MolfileReadSgroupOfPolymer
    // INCHI✔️❌: complete source frame follows verbatim; checked source-heap object access and stack-buffer modeling add overhead.
    /*
+int MolfileReadSgroupOfPolymer( MOL_FMT_CTAB* ctab,
                                MOL_FMT_HEADER_BLOCK *pHdr,
                                INCHI_IOSTREAM *inp_file,
                                char line[MOL_FMT_INPLINELEN],
                                char *szType,
                                char *p,
                                int err,
                                char *pStrErr )
{
    S_SHORT  num_entries;
    S_SHORT  num_atoms = ctab->n_atoms;
    S_SHORT  num_bonds = ctab->n_bonds;
    int j, index = -1, ret;
    char  stmp[4], stmplong[81];
    S_SHORT sg_nums[8], sg_num = -1, sg_atoms[15], sg_bonds[15], tmp;
    int q, fail = 0, len;
    /* djb-rwth: removing redundant variables */
#if ( DEBUG_POLYMERS == 1 )
    debug_polymers = 1;
#elif  ( DEBUG_POLYMERS == 2 )
    debug_polymers = 2;
#endif
    /* djb-rwth: removing redundant code */

    /* Check for possible lead codes */

    /*    STY - Sgroup type
                Polymer-related recognized types are:
                SRU = SRU type,
                MON = monomer,
                COP = copolymer,
                MER = Mer type,
                MOD
                CRO
    */
    if ( !strcmp( szType, "STY" )
         && 0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p )
         && 1 <= num_entries && num_entries <= 8 )
    {
        for (j = 0; j < num_entries; j++)
        {
            fail = 0 > MolfileReadField( &sg_nums[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                0 > MolfileReadField( stmp, 4, MOL_FMT_STRING_DATA, &p );
            if (!fail)
            {
                int type = MOL_FMT_M_STY_NON;

                lrtrim( stmp, &len );

                if (!strcmp( stmp, "SRU" ))
                {
                    type = MOL_FMT_M_STY_SRU;
                }
                else if (!strcmp( stmp, "MON" ))
                {
                    type = MOL_FMT_M_STY_MON;
                }
                else if (!strcmp( stmp, "COP" ))
                {
                    type = MOL_FMT_M_STY_COP;
                }
                else if (!strcmp( stmp, "MOD" ))
                {
                    type = MOL_FMT_M_STY_MOD;
                }
                else if (!strcmp( stmp, "CRO" ))
                {
                    type = MOL_FMT_M_STY_CRO;
                }
                else if (!strcmp( stmp, "MER" ))
                {
                    type = MOL_FMT_M_STY_MER;
                }
                else
                {
                    fail = 1;
                }
                if (!fail)
                {
                    index = MolFmtSgroups_GetIndexBySgroupId( sg_nums[j], &( ctab->sgroups ) );
                    if (-1 == index)
                    {
                        ret = MolFmtSgroups_Append( &ctab->sgroups, sg_nums[j], type );
                        if (0 != ret)
                            fail = 1;
                        else
                            index = ctab->sgroups.used - 1;
                    }
                    else
                    {
                        ctab->sgroups.group[index]->type = type;
                    }
                    if (!fail)
                    {
                        for (q = 0; q < 4; q++)
                        {
                            ctab->sgroups.group[index]->xbr1[q] = -777777.777;
                            ctab->sgroups.group[index]->xbr2[q] = -777777.777;
                        }
                        ctab->sgroups.group[index]->smt[0] = '\0';
                    }
                }
            }
            if (fail)
            {
                err = 5;
                goto err_exit;
            }
        }
    }
    /*
    SST - Polymer Sgroup subtypes:
                ALT = alternating,
                RAN = random,
                BLK = block
    */
    else if (!strcmp( szType, "SST" ) &&
              0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
              1 <= num_entries && num_entries <= 8)
    {
        for (j = 0; j < num_entries; j++)
        {
            fail = 0 > MolfileReadField( &sg_nums[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                0 > MolfileReadField( stmp, 4, MOL_FMT_STRING_DATA, &p );
            if (!fail)
            {
                index = MolFmtSgroups_GetIndexBySgroupId( sg_nums[j], &( ctab->sgroups ) );
                if (-1 == index)
                {
                    fail = 1;
                }
            }
            if (!fail)
            {
                ctab->sgroups.group[index]->subtype = MOL_FMT_M_SST_NON;
                lrtrim( stmp, &len );
                if (!strcmp( stmp, "ALT" ))
                {
                    ctab->sgroups.group[index]->subtype = MOL_FMT_M_SST_ALT;
                }
                else if (!strcmp( stmp, "RAN" ))
                {
                    ctab->sgroups.group[index]->subtype = MOL_FMT_M_SST_RAN;
                }
                else if (!strcmp( stmp, "BLO" ))
                {
                    ctab->sgroups.group[index]->subtype = MOL_FMT_M_SST_BLK;
                }
                else if (!strcmp( stmp, "BLK" ))
                {
                    ctab->sgroups.group[index]->subtype = MOL_FMT_M_SST_BLK;
                }
                else
                {
                    fail = 1;
                    break;
                }
            }
        }
        if (fail)
        {
            err = 6;
            goto err_exit;
        }
    }
    /*    SLB - Sgroup Labels */
    else if (!strcmp( szType, "SLB" ) &&
              0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
              1 <= num_entries && num_entries <= 8)
    {
        for (j = 0; j < num_entries; j++)
        {
            fail = 0 > MolfileReadField( &sg_nums[j], 0, MOL_FMT_SHORT_INT_DATA, &p )
                   || 0 > MolfileReadField( &tmp, 0, MOL_FMT_SHORT_INT_DATA, &p );
            if (!fail)
            {
                index = MolFmtSgroups_GetIndexBySgroupId( sg_nums[j], &( ctab->sgroups ) );
                if (-1 == index)
                {
                    fail = 1;
                }
            }
            if (!fail)
            {
                ctab->sgroups.group[index]->label = tmp;
            }
        }
        if (fail)
        {
            err = 7;
            goto err_exit;
        }
    }
    /*    SCN - Sgroup Connectivity
                    HH = head-to-head,
                    HT = head-to-tail,
                    EU = either unknown.
                    Left justified.
    */
    else if (!strcmp( szType, "SCN" ) &&
              0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ) &&
              1 <= num_entries && num_entries <= 8)
    {
        for (j = 0; j < num_entries; j++)
        {
            fail = 0 > MolfileReadField( &sg_nums[j], 0, MOL_FMT_SHORT_INT_DATA, &p )
                   || 0 > MolfileReadField( stmp, 4, MOL_FMT_STRING_DATA, &p );
            if (!fail)
            {
                index = MolFmtSgroups_GetIndexBySgroupId( sg_nums[j], &( ctab->sgroups ) );
                if (-1 == index)
                {
                    fail = 1;
                }
            }
            if (!fail)
            {
                ctab->sgroups.group[index]->conn = MOL_FMT_M_CONN_NON;
                lrtrim( stmp, &len );
                if (!strcmp( stmp, "HT" ))
                {
                    ctab->sgroups.group[index]->conn = MOL_FMT_M_CONN_HT;
                }
                else if (!strcmp( stmp, "HH" ))
                {
                    ctab->sgroups.group[index]->conn = MOL_FMT_M_CONN_HH;
                }
                else if (!strcmp( stmp, "EU" ))
                {
                    ctab->sgroups.group[index]->conn = MOL_FMT_M_CONN_EU;
                }
                else
                {
                    fail = 1;       /* NB: we do not allow explicit different abbreviation - but note that  */
                                    /* totally skipping SCN line is allowed ("EU" will be inserted further) */
                }
            }
            if (fail)
            {
                err = 8;
                goto err_exit;
            }
        }
    }
    /* SAL - Sgroup atoms list  */
    else if (!strcmp( szType, "SAL" ))
    {
        if (0 < MolfileReadField( &sg_num, 4, MOL_FMT_SHORT_INT_DATA, &p ) &&
              0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ))
        {
            index = MolFmtSgroups_GetIndexBySgroupId( sg_num, &( ctab->sgroups ) );
            if (-1 == index)
            {
                fail = 1;
            }
            if (!fail)
            {
                for (j = 0; j < num_entries; j++)
                {
                    if (0 > MolfileReadField( &sg_atoms[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                        sg_atoms[j] <= 0 || sg_atoms[j] > num_atoms)
                    {
                        fail = 1;
                        break;
                    }
                }
            }
            if (!fail)
            {
                for (j = 0; j < num_entries; j++)
                {
                    if (0 != IntArray_Append( &( ctab->sgroups.group[index]->alist ), sg_atoms[j] ))
                    {
                        fail = 1;
                        break;
                    }
                }
            }
        }
        else
        {
            fail = 1;
        }
        if (fail)
        {
            err = 9;
            goto err_exit;
        }
    }
    /*    SBL - Sgroup bonds list */
    else if (!strcmp( szType, "SBL" ))
    {
        if (0 < MolfileReadField( &sg_num, 4, MOL_FMT_SHORT_INT_DATA, &p ) &&
            0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ))
        {
            index = MolFmtSgroups_GetIndexBySgroupId( sg_num, &( ctab->sgroups ) );
            if (-1 == index)
            {
                fail = 1;
            }
            if (!fail)
            {
                for (j = 0; j < num_entries; j++)
                {
                    if (0 > MolfileReadField( &sg_bonds[j], 0, MOL_FMT_SHORT_INT_DATA, &p ) ||
                         sg_bonds[j] <= 0 || sg_bonds[j] > num_bonds)
                    {
                        fail = 1;
                        break;
                    }
                }
            }
            if (!fail)
            {
                for (j = 0; j < num_entries; j++)
                {
                    if (0 != IntArray_Append( &( ctab->sgroups.group[index]->blist ), sg_bonds[j] ))
                    {
                        fail = 1;
                        break;
                    }
                }
            }
        }
        else
        {
            fail = 1;
        }
        if (fail)
        {
            err = 10;
            goto err_exit;
        }
    }
    /* SDI */
    else if (!strcmp( szType, "SDI" ))
    {
        double x[4];

        if (0 < MolfileReadField( &sg_num, 4, MOL_FMT_SHORT_INT_DATA, &p ) &&
            0 < MolfileReadField( &num_entries, 3, MOL_FMT_SHORT_INT_DATA, &p ))
        {
            index = MolFmtSgroups_GetIndexBySgroupId( sg_num, &( ctab->sgroups ) );
            if (-1 == index)
            {
                fail = 1;
            }
            else if (num_entries != 4)
            {
                fail = 1;
            }
            if (!fail)
            {
                for (j = 0; j < num_entries; j++)
                {
                    if (0 > MolfileReadField( &x[j], 0, MOL_FMT_DOUBLE_DATA, &p ))
                    {
                        fail = 1;
                        break;
                    }
                }
            }
            if (!fail)
            {
                if (fabs( -fabs( ctab->sgroups.group[index]->xbr1[0] ) + 777777.777 ) < 1.e-7)    /* brkt1 coords not yet here */
                {
                    for (q = 0; q < 4; q++)
                    {
                        ctab->sgroups.group[index]->xbr1[q] = x[q];
                    }
                }
                else
                {
                    for (q = 0; q < 4; q++)
                    {
                        ctab->sgroups.group[index]->xbr2[q] = x[q];
                    }
                }
            }
        }
        else
        {
            fail = 1;
        }
        if (fail)
        {
            err = 11;
            goto err_exit;
        }
    }
    /* SMT - Sgroup Subscript */
    else if (!strcmp( szType, "SMT" ))
    {
        index = -1;
        if (0 < MolfileReadField( &sg_num, 4, MOL_FMT_SHORT_INT_DATA, &p ) &&
            0 < MolfileReadField( stmplong, 80, MOL_FMT_STRING_DATA, &p ))
        {
            index = MolFmtSgroups_GetIndexBySgroupId( sg_num, &( ctab->sgroups ) );
        }
        if (-1 == index)
        {
            fail = 1;
        }
        if (!fail)
        {
            lrtrim( stmplong, &len );
            strcpy(ctab->sgroups.group[index]->smt, stmplong);
        }
        if (fail)
        {
            err = 11;
            goto err_exit;
        }
    }


    ITRACE_( "\n" );
    return 0;

err_exit:
    MolFmtSgroups_Free( &ctab->sgroups );

    return err; /* ignore polymeric error for now */
}

    */
    // END INCHI C FUNCTION: MolfileReadSgroupOfPolymer
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadSgroupOfPolymer
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; DEBUG_POLYMERS is 0 and ITRACE_ short-circuits without evaluating diagnostics.
    // INCHI✔️❌: S_SHORT is signed 16-bit; ctab atom and bond counts narrow with GCC low-16-bit semantics.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileReadSgroupOfPolymer

    let num_atoms = ctab.n_atoms as i16;
    let num_bonds = ctab.n_bonds as i16;
    let mut num_entries = 0_i16;
    let mut sg_num = -1_i16;
    let mut sg_nums = [0_i16; 8];
    let mut sg_atoms = [0_i16; 15];
    let mut sg_bonds = [0_i16; 15];
    let mut stmp = [0_i8; 5];
    let mut stmplong = [0_i8; 81];

    if mol_fmt1_c_string_eq(sz_type, b"STY") {
        let count_result = MolfileReadField(heap, MolfileFieldData::Short(&mut num_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)?;
        if count_result > 0 && (1..=8).contains(&num_entries) {
            for j in 0..i32::from(num_entries) {
                let slot = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let mut fail = MolfileReadField(heap, MolfileFieldData::Short(&mut sg_nums[slot]), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::String(&mut stmp), 4, i32::from(MOL_FMT_STRING_DATA), &mut p)? < 0;
                let mut type_ = MOL_FMT_M_STY_NON as i32;
                if !fail {
                    let _ = mol_fmt1_trim_array(heap, &mut stmp)?;
                    type_ = if mol_fmt1_c_string_eq(&stmp, b"SRU") {
                        MOL_FMT_M_STY_SRU as i32
                    } else if mol_fmt1_c_string_eq(&stmp, b"MON") {
                        MOL_FMT_M_STY_MON as i32
                    } else if mol_fmt1_c_string_eq(&stmp, b"COP") {
                        MOL_FMT_M_STY_COP as i32
                    } else if mol_fmt1_c_string_eq(&stmp, b"MOD") {
                        MOL_FMT_M_STY_MOD as i32
                    } else if mol_fmt1_c_string_eq(&stmp, b"CRO") {
                        MOL_FMT_M_STY_CRO as i32
                    } else if mol_fmt1_c_string_eq(&stmp, b"MER") {
                        MOL_FMT_M_STY_MER as i32
                    } else {
                        fail = true;
                        MOL_FMT_M_STY_NON as i32
                    };
                }
                if !fail {
                    let mut index = MolFmtSgroups_GetIndexBySgroupId(heap, i32::from(sg_nums[slot]), &ctab.sgroups)?;
                    if index == -1 {
                        if MolFmtSgroups_Append(heap, Some(&mut ctab.sgroups), i32::from(sg_nums[slot]), type_)? != 0 {
                            fail = true;
                        } else {
                            index = ctab.sgroups.used.wrapping_sub(1);
                        }
                    } else {
                        let pointer = mol_fmt1_sgroup_pointer(heap, &ctab.sgroups, index)?;
                        heap.slice_mut(pointer)?[0].type_ = type_;
                    }
                    if !fail {
                        let pointer = mol_fmt1_sgroup_pointer(heap, &ctab.sgroups, index)?;
                        let group = &mut heap.slice_mut(pointer)?[0];
                        group.xbr1 = [-777777.777; 4];
                        group.xbr2 = [-777777.777; 4];
                        group.smt[0] = 0;
                    }
                }
                if fail {
                    err = 5;
                    MolFmtSgroups_Free(heap, Some(&mut ctab.sgroups))?;
                    return Ok(err);
                }
            }
        }
    } else if mol_fmt1_c_string_eq(sz_type, b"SST") {
        let count_result = MolfileReadField(heap, MolfileFieldData::Short(&mut num_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)?;
        if count_result > 0 && (1..=8).contains(&num_entries) {
            let mut fail = false;
            for j in 0..i32::from(num_entries) {
                let slot = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                fail = MolfileReadField(heap, MolfileFieldData::Short(&mut sg_nums[slot]), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::String(&mut stmp), 4, i32::from(MOL_FMT_STRING_DATA), &mut p)? < 0;
                let index = if fail { -1 } else { MolFmtSgroups_GetIndexBySgroupId(heap, i32::from(sg_nums[slot]), &ctab.sgroups)? };
                if index == -1 {
                    fail = true;
                }
                if !fail {
                    let pointer = mol_fmt1_sgroup_pointer(heap, &ctab.sgroups, index)?;
                    heap.slice_mut(pointer)?[0].subtype = MOL_FMT_M_SST_NON as i32;
                    let _ = mol_fmt1_trim_array(heap, &mut stmp)?;
                    let subtype = if mol_fmt1_c_string_eq(&stmp, b"ALT") {
                        Some(MOL_FMT_M_SST_ALT as i32)
                    } else if mol_fmt1_c_string_eq(&stmp, b"RAN") {
                        Some(MOL_FMT_M_SST_RAN as i32)
                    } else if mol_fmt1_c_string_eq(&stmp, b"BLO") || mol_fmt1_c_string_eq(&stmp, b"BLK") {
                        Some(MOL_FMT_M_SST_BLK as i32)
                    } else {
                        None
                    };
                    if let Some(subtype) = subtype {
                        heap.slice_mut(pointer)?[0].subtype = subtype;
                    } else {
                        fail = true;
                        break;
                    }
                }
            }
            if fail {
                err = 6;
                MolFmtSgroups_Free(heap, Some(&mut ctab.sgroups))?;
                return Ok(err);
            }
        }
    } else if mol_fmt1_c_string_eq(sz_type, b"SLB") {
        let count_result = MolfileReadField(heap, MolfileFieldData::Short(&mut num_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)?;
        if count_result > 0 && (1..=8).contains(&num_entries) {
            let mut fail = false;
            for j in 0..i32::from(num_entries) {
                let slot = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let mut tmp = 0_i16;
                fail = MolfileReadField(heap, MolfileFieldData::Short(&mut sg_nums[slot]), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::Short(&mut tmp), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0;
                let index = if fail { -1 } else { MolFmtSgroups_GetIndexBySgroupId(heap, i32::from(sg_nums[slot]), &ctab.sgroups)? };
                if index == -1 {
                    fail = true;
                } else if !fail {
                    let pointer = mol_fmt1_sgroup_pointer(heap, &ctab.sgroups, index)?;
                    heap.slice_mut(pointer)?[0].label = i32::from(tmp);
                }
            }
            if fail {
                err = 7;
                MolFmtSgroups_Free(heap, Some(&mut ctab.sgroups))?;
                return Ok(err);
            }
        }
    } else if mol_fmt1_c_string_eq(sz_type, b"SCN") {
        let count_result = MolfileReadField(heap, MolfileFieldData::Short(&mut num_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)?;
        if count_result > 0 && (1..=8).contains(&num_entries) {
            for j in 0..i32::from(num_entries) {
                let slot = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let mut fail = MolfileReadField(heap, MolfileFieldData::Short(&mut sg_nums[slot]), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                    || MolfileReadField(heap, MolfileFieldData::String(&mut stmp), 4, i32::from(MOL_FMT_STRING_DATA), &mut p)? < 0;
                let index = if fail { -1 } else { MolFmtSgroups_GetIndexBySgroupId(heap, i32::from(sg_nums[slot]), &ctab.sgroups)? };
                if index == -1 {
                    fail = true;
                }
                if !fail {
                    let pointer = mol_fmt1_sgroup_pointer(heap, &ctab.sgroups, index)?;
                    heap.slice_mut(pointer)?[0].conn = MOL_FMT_M_CONN_NON as i32;
                    let _ = mol_fmt1_trim_array(heap, &mut stmp)?;
                    let conn = if mol_fmt1_c_string_eq(&stmp, b"HT") {
                        Some(MOL_FMT_M_CONN_HT as i32)
                    } else if mol_fmt1_c_string_eq(&stmp, b"HH") {
                        Some(MOL_FMT_M_CONN_HH as i32)
                    } else if mol_fmt1_c_string_eq(&stmp, b"EU") {
                        Some(MOL_FMT_M_CONN_EU as i32)
                    } else {
                        None
                    };
                    if let Some(conn) = conn {
                        heap.slice_mut(pointer)?[0].conn = conn;
                    } else {
                        fail = true;
                    }
                }
                if fail {
                    err = 8;
                    MolFmtSgroups_Free(heap, Some(&mut ctab.sgroups))?;
                    return Ok(err);
                }
            }
        }
    }

    else if mol_fmt1_c_string_eq(sz_type, b"SAL") {
        let mut fail = !(MolfileReadField(heap, MolfileFieldData::Short(&mut sg_num), 4, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? > 0
            && MolfileReadField(heap, MolfileFieldData::Short(&mut num_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? > 0);
        let index = if fail { -1 } else { MolFmtSgroups_GetIndexBySgroupId(heap, i32::from(sg_num), &ctab.sgroups)? };
        if index == -1 {
            fail = true;
        }
        if !fail {
            for j in 0..i32::from(num_entries) {
                let slot = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let Some(atom) = sg_atoms.get_mut(slot) else { return Err(SourceHeapError::PointerOutOfBounds); };
                if MolfileReadField(heap, MolfileFieldData::Short(atom), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                    || *atom <= 0 || *atom > num_atoms
                {
                    fail = true;
                    break;
                }
            }
        }
        if !fail {
            for j in 0..i32::from(num_entries) {
                let slot = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let pointer = mol_fmt1_sgroup_pointer(heap, &ctab.sgroups, index)?;
                let mut group = heap.slice(pointer.as_const())?[0].clone();
                let append = IntArray_Append(heap, Some(&mut group.alist), i32::from(sg_atoms[slot]));
                heap.slice_mut(pointer)?[0] = group;
                if append? != 0 {
                    fail = true;
                    break;
                }
            }
        }
        if fail {
            err = 9;
            MolFmtSgroups_Free(heap, Some(&mut ctab.sgroups))?;
            return Ok(err);
        }
    } else if mol_fmt1_c_string_eq(sz_type, b"SBL") {
        let mut fail = !(MolfileReadField(heap, MolfileFieldData::Short(&mut sg_num), 4, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? > 0
            && MolfileReadField(heap, MolfileFieldData::Short(&mut num_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? > 0);
        let index = if fail { -1 } else { MolFmtSgroups_GetIndexBySgroupId(heap, i32::from(sg_num), &ctab.sgroups)? };
        if index == -1 {
            fail = true;
        }
        if !fail {
            for j in 0..i32::from(num_entries) {
                let slot = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let Some(bond) = sg_bonds.get_mut(slot) else { return Err(SourceHeapError::PointerOutOfBounds); };
                if MolfileReadField(heap, MolfileFieldData::Short(bond), 0, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? < 0
                    || *bond <= 0 || *bond > num_bonds
                {
                    fail = true;
                    break;
                }
            }
        }
        if !fail {
            for j in 0..i32::from(num_entries) {
                let slot = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let pointer = mol_fmt1_sgroup_pointer(heap, &ctab.sgroups, index)?;
                let mut group = heap.slice(pointer.as_const())?[0].clone();
                let append = IntArray_Append(heap, Some(&mut group.blist), i32::from(sg_bonds[slot]));
                heap.slice_mut(pointer)?[0] = group;
                if append? != 0 {
                    fail = true;
                    break;
                }
            }
        }
        if fail {
            err = 10;
            MolFmtSgroups_Free(heap, Some(&mut ctab.sgroups))?;
            return Ok(err);
        }
    }

    else if mol_fmt1_c_string_eq(sz_type, b"SDI") {
        let mut x = [0.0_f64; 4];
        let mut fail = !(MolfileReadField(heap, MolfileFieldData::Short(&mut sg_num), 4, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? > 0
            && MolfileReadField(heap, MolfileFieldData::Short(&mut num_entries), 3, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? > 0);
        let index = if fail { -1 } else { MolFmtSgroups_GetIndexBySgroupId(heap, i32::from(sg_num), &ctab.sgroups)? };
        if index == -1 || num_entries != 4 {
            fail = true;
        }
        if !fail {
            for j in 0..i32::from(num_entries) {
                let slot = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                if MolfileReadField(heap, MolfileFieldData::Double(&mut x[slot]), 0, i32::from(MOL_FMT_DOUBLE_DATA), &mut p)? < 0 {
                    fail = true;
                    break;
                }
            }
        }
        if !fail {
            let pointer = mol_fmt1_sgroup_pointer(heap, &ctab.sgroups, index)?;
            let group = &mut heap.slice_mut(pointer)?[0];
            if (-group.xbr1[0].abs() + 777777.777).abs() < 1.0e-7 {
                group.xbr1 = x;
            } else {
                group.xbr2 = x;
            }
        }
        if fail {
            err = 11;
            MolFmtSgroups_Free(heap, Some(&mut ctab.sgroups))?;
            return Ok(err);
        }
    } else if mol_fmt1_c_string_eq(sz_type, b"SMT") {
        let mut index = -1_i32;
        if MolfileReadField(heap, MolfileFieldData::Short(&mut sg_num), 4, i32::from(MOL_FMT_SHORT_INT_DATA), &mut p)? > 0
            && MolfileReadField(heap, MolfileFieldData::String(&mut stmplong), 80, i32::from(MOL_FMT_STRING_DATA), &mut p)? > 0
        {
            index = MolFmtSgroups_GetIndexBySgroupId(heap, i32::from(sg_num), &ctab.sgroups)?;
        }
        if index == -1 {
            err = 11;
            MolFmtSgroups_Free(heap, Some(&mut ctab.sgroups))?;
            return Ok(err);
        }
        let _ = mol_fmt1_trim_array(heap, &mut stmplong)?;
        let length = stmplong.iter().position(|byte| *byte == 0).ok_or(SourceHeapError::MissingNulTerminator)?;
        let pointer = mol_fmt1_sgroup_pointer(heap, &ctab.sgroups, index)?;
        let group = &mut heap.slice_mut(pointer)?[0];
        let destination = group.smt.get_mut(..=length).ok_or(SourceHeapError::PointerOutOfBounds)?;
        destination.copy_from_slice(&stmplong[..=length]);
    }

    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn MolfileTreatPseudoElementAtoms(
    heap: &mut SourceHeap,
    ctab: &mut MOL_FMT_CTAB,
    pseudos_allowed: i32,
    err: &mut i32,
    mut errors: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:1860 MolfileTreatPseudoElementAtoms
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
static int MolfileTreatPseudoElementAtoms( MOL_FMT_CTAB* ctab,
                                           int pseudos_allowed,
                                           int *err,
                                           char *pStrErr )
{
    int i, nzz = 0;

    /* djb-rwth: addressing coverity ID #499499 -- TREAT_ERR properly used in all cases */

    for (i = 0; i < ctab->n_atoms; i++)
    {
        int is_zz = 0, is_star = 0;

        /* Zy is specifically disabled */
        if (!strcmp(ctab->atoms[i].symbol, "Zy"))
        {
            TREAT_ERR( *err, ( 70 + 6 ), "Invalid element(s):" );
            TREAT_ERR( *err, ( 70 + 6 ), ctab->atoms[i].symbol );
        }

        is_star = !strcmp( ctab->atoms[i].symbol, "*" );
        if (!is_star)
        {
            is_zz = !strcmp( ctab->atoms[i].symbol, "Zz" );
        }

        if (is_star || is_zz)
        {
            nzz++;
            if (0== pseudos_allowed)
            {
                /* Pseudoelements totally disabled */
                TREAT_ERR( *err, ( 70 + 6 ), "Invalid element(s):" );
                TREAT_ERR( *err, ( 70 + 6 ), ctab->atoms[i].symbol );
            }
            else
            {
                /* That's allowed, if it's star then convert to Zz */
                if (is_star)
                {
                    mystrncpy( ctab->atoms[i].symbol, "Zz", sizeof( "Zz" ) );
                }
            }
        }
    }

    return nzz;
}
    */
    // END INCHI C FUNCTION: MolfileTreatPseudoElementAtoms
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileTreatPseudoElementAtoms
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; strcmp is the active libc behavior.
    // INCHI✔️✔️: TREAT_ERR preserves a pre-existing nonzero error and always calls AddErrorMessage; mystrncpy is the completed util.c behavior.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: MolfileTreatPseudoElementAtoms

    fn treat_symbol_error(
        err: &mut i32,
        errors: &mut Option<&mut [i8]>,
        symbol: &[i8],
    ) -> Result<(), SourceHeapError> {
        mol_fmt1_treat_error(err, 76, errors.as_deref_mut(), b"Invalid element(s):")?;
        let _ = AddErrorMessage(errors.as_deref_mut(), Some(symbol))?;
        Ok(())
    }

    let mut nzz = 0_i32;
    let count = if ctab.n_atoms > 0 {
        usize::try_from(ctab.n_atoms).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    if count == 0 {
        return Ok(0);
    }
    let atoms = heap.slice_mut(ctab.atoms)?;
    for atom in atoms.get_mut(..count).ok_or(SourceHeapError::PointerOutOfBounds)? {
        if mol_fmt1_c_string_eq(&atom.symbol, b"Zy") {
            treat_symbol_error(err, &mut errors, &atom.symbol)?;
        }
        let is_star = mol_fmt1_c_string_eq(&atom.symbol, b"*");
        let is_zz = !is_star && mol_fmt1_c_string_eq(&atom.symbol, b"Zz");
        if is_star || is_zz {
            nzz = nzz.wrapping_add(1);
            if pseudos_allowed == 0 {
                treat_symbol_error(err, &mut errors, &atom.symbol)?;
            } else if is_star {
                let _ = mystrncpy_slice(
                    Some(&mut atom.symbol),
                    Some(&[b'Z' as i8, b'z' as i8, 0]),
                    3,
                )?;
            }
        }
    }
    Ok(nzz)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        FILE, INCHI_IOS_STRING, INCHI_IOS_TYPE_FILE, INCHI_IOS_TYPE_STRING, MOL_COORD,
        MOL_FMT_ATOM, MOL_FMT_BOND,
    };

    fn bytes(value: &str) -> Vec<i8> {
        value
            .as_bytes()
            .iter()
            .map(|byte| *byte as i8)
            .chain(std::iter::once(0))
            .collect()
    }

    fn text(value: &[i8]) -> String {
        String::from_utf8(
            value
                .iter()
                .take_while(|byte| **byte != 0)
                .map(|byte| *byte as u8)
                .collect(),
        )
        .unwrap()
    }

    fn input() -> String {
        let line2 = format!(
            "{}{}{}{}{}{}{}{}{}{}{}{}",
            "AB",
            "PROGNAME",
            "01",
            "02",
            "24",
            "12",
            "34",
            "2D",
            " 7",
            "      1.25",
            "        -2.5",
            "    42"
        );
        assert_eq!(line2.len(), 52);
        format!("  Molecule Name  \r\n{line2}\n  Comment text  \n")
    }

    fn string_stream(heap: &mut SourceHeap, value: &str) -> INCHI_IOSTREAM {
        let data = heap.allocate(bytes(value)).unwrap();
        INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING {
                pStr: data,
                nAllocatedLength: value.len() as i32 + 1,
                nUsedLength: value.len() as i32,
                nPtr: 0,
            },
            ..INCHI_IOSTREAM::default()
        }
    }

    fn polymer_ctab(heap: &mut SourceHeap) -> MOL_FMT_CTAB {
        let mut ctab = MOL_FMT_CTAB {
            n_atoms: 10,
            n_bonds: 8,
            ..MOL_FMT_CTAB::default()
        };
        assert_eq!(MolFmtSgroups_Alloc(heap, Some(&mut ctab.sgroups), 8), Ok(0));
        ctab
    }

    fn polymer_parse(
        heap: &mut SourceHeap,
        ctab: &mut MOL_FMT_CTAB,
        type_: &str,
        fields: &str,
        err: i32,
    ) -> i32 {
        let mut input = bytes(fields);
        input.resize(256, 0);
        let pointer = heap.allocate_model_storage(input).unwrap();
        let type_bytes = bytes(type_);
        let result = MolfileReadSgroupOfPolymer(
            heap,
            ctab,
            None,
            None,
            SourceMutPointer::null(),
            &type_bytes,
            pointer,
            err,
            None,
        )
        .unwrap();
        heap.free(pointer).unwrap();
        result
    }

    fn polymer_group(
        heap: &SourceHeap,
        ctab: &MOL_FMT_CTAB,
        id: i32,
    ) -> crate::source_types::MOL_FMT_SGROUP {
        let index = MolFmtSgroups_GetIndexBySgroupId(heap, id, &ctab.sgroups).unwrap();
        assert_ne!(index, -1);
        let pointer = mol_fmt1_sgroup_pointer(heap, &ctab.sgroups, index).unwrap();
        heap.slice(pointer.as_const()).unwrap()[0].clone()
    }

    #[test]
    fn source_port__mol_fmt1__molfilereadpropblock__line_960() {
        fn atom(symbol: &str, charge: i8, radical: i8, mass: i8) -> MOL_FMT_ATOM {
            let mut atom = MOL_FMT_ATOM {
                charge,
                radical,
                mass_difference: mass,
                ..MOL_FMT_ATOM::default()
            };
            for (destination, source) in atom.symbol.iter_mut().zip(symbol.bytes()) {
                *destination = source as i8;
            }
            atom
        }
        fn versioned(ctab: &mut MOL_FMT_CTAB) {
            ctab.version_string[..6].copy_from_slice(&[
                b'V' as i8, b'2' as i8, b'0' as i8, b'0' as i8, b'0' as i8, 0,
            ]);
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                atom("C", 5, 3, 9),
                atom("O", 6, 1, 8),
                atom("D", 7, 2, 7),
            ])
            .unwrap();
        let mut ctab = MOL_FMT_CTAB {
            n_atoms: 3,
            atoms,
            ..MOL_FMT_CTAB::default()
        };
        versioned(&mut ctab);
        let source = concat!(
            "M  CHG  1   1  -1\n",
            "M  RAD  1   2   2\n",
            "M  ISO  2   1  13   2  16\n",
            "M  REG 12345\n",
            "M  END\n"
        );
        let mut stream = string_stream(&mut heap, source);
        let mut header = MOL_FMT_HEADER_BLOCK::default();
        let mut errors = [0_i8; 1024];
        assert_eq!(
            MolfileReadPropBlock(
                &mut heap,
                &mut ctab,
                &mut header,
                Some(&mut stream),
                1,
                0,
                Some(&mut errors),
                0,
            ),
            Ok(0)
        );
        assert_eq!(header.internal_regno, 12345);
        assert_eq!(text(&errors), "");
        let parsed = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(
            (
                parsed[0].charge,
                parsed[0].radical,
                parsed[0].mass_difference
            ),
            (-1, 0, 1)
        );
        assert_eq!(
            (
                parsed[1].charge,
                parsed[1].radical,
                parsed[1].mass_difference
            ),
            (0, 2, ZERO_ATW_DIFF as i8)
        );
        assert_eq!(
            (
                parsed[2].charge,
                parsed[2].radical,
                parsed[2].mass_difference
            ),
            (0, 0, 0)
        );

        let alias_atoms = heap
            .allocate_model_storage(vec![atom("C", 4, 3, 5), atom("C", 5, 2, 6)])
            .unwrap();
        let mut alias_ctab = MOL_FMT_CTAB {
            n_atoms: 2,
            atoms: alias_atoms,
            ..MOL_FMT_CTAB::default()
        };
        versioned(&mut alias_ctab);
        let mut alias_stream = string_stream(
            &mut heap,
            concat!(
                "A    1\n",
                " \tD  \r\n",
                "A    2\n",
                "D+\n",
                "M  CHG  1   1  -7\n",
                "M  END\n"
            ),
        );
        errors.fill(0);
        assert_eq!(
            MolfileReadPropBlock(
                &mut heap,
                &mut alias_ctab,
                &mut header,
                Some(&mut alias_stream),
                1,
                0,
                Some(&mut errors),
                0,
            ),
            Ok(0)
        );
        let aliases = heap.slice(alias_atoms.as_const()).unwrap();
        assert_eq!(&aliases[0].symbol[..2], &[b'H' as i8, 0]);
        assert_eq!(
            (
                aliases[0].mass_difference,
                aliases[0].charge,
                aliases[0].radical,
                aliases[0].atom_aliased_flag
            ),
            (1, 0, 0, 1)
        );
        assert_eq!(&aliases[1].symbol[..2], &[b'D' as i8, 0]);
        assert_eq!(
            (
                aliases[1].mass_difference,
                aliases[1].charge,
                aliases[1].radical,
                aliases[1].atom_aliased_flag
            ),
            (0, 1, 0, 1)
        );

        let control_atoms = heap
            .allocate_model_storage(vec![atom("N", 2, 1, 3)])
            .unwrap();
        let mut control_ctab = MOL_FMT_CTAB {
            n_atoms: 1,
            atoms: control_atoms,
            ..MOL_FMT_CTAB::default()
        };
        versioned(&mut control_ctab);
        let mut control_stream = string_stream(
            &mut heap,
            concat!(
                "V  FOO\n",
                "G  FOO\n",
                "M  CHG  1   1   9\n",
                "S  SKP  1\n",
                "M  RAD  1   1   3\n",
                "A    9\n",
                "M  ISO  1   1  15\n",
                "Q  FOO\n",
                "M  END\n"
            ),
        );
        assert_eq!(
            MolfileReadPropBlock(
                &mut heap,
                &mut control_ctab,
                &mut header,
                Some(&mut control_stream),
                1,
                0,
                None,
                0
            ),
            Ok(0)
        );
        assert_eq!(
            (
                heap.slice(control_atoms.as_const()).unwrap()[0].charge,
                heap.slice(control_atoms.as_const()).unwrap()[0].radical,
                heap.slice(control_atoms.as_const()).unwrap()[0].mass_difference
            ),
            (2, 1, 3)
        );

        for (record, message) in [
            ("M  CHG  1   9   1\nM  END\n", "Charge not recognized:"),
            ("M  RAD  1   1   9\nM  END\n", "Radical not recognized:"),
            (
                "M  ISO  1   9  13\nM  END\n",
                "Isotopic data not recognized:",
            ),
        ] {
            let one_atom = heap
                .allocate_model_storage(vec![atom("C", 1, 2, 3)])
                .unwrap();
            let mut invalid_ctab = MOL_FMT_CTAB {
                n_atoms: 1,
                atoms: one_atom,
                ..MOL_FMT_CTAB::default()
            };
            versioned(&mut invalid_ctab);
            let mut invalid_stream = string_stream(&mut heap, record);
            errors.fill(0);
            assert_eq!(
                MolfileReadPropBlock(
                    &mut heap,
                    &mut invalid_ctab,
                    &mut header,
                    Some(&mut invalid_stream),
                    1,
                    0,
                    Some(&mut errors),
                    0
                ),
                Ok(0)
            );
            assert!(text(&errors).contains(message), "{}", text(&errors));
        }

        let mut no_polymer_ctab = MOL_FMT_CTAB {
            n_atoms: 1,
            atoms: control_atoms,
            ..MOL_FMT_CTAB::default()
        };
        versioned(&mut no_polymer_ctab);
        let mut no_polymer_stream = string_stream(&mut heap, "M  STY  1   1 SRU\nM  END\n");
        errors.fill(0);
        assert_eq!(
            MolfileReadPropBlock(
                &mut heap,
                &mut no_polymer_ctab,
                &mut header,
                Some(&mut no_polymer_stream),
                0,
                0,
                Some(&mut errors),
                0
            ),
            Ok(0)
        );
        assert_eq!(text(&errors), "Ignore polymer data");

        let mut polymer_ctab = polymer_ctab(&mut heap);
        polymer_ctab.atoms = heap
            .allocate_model_storage(vec![MOL_FMT_ATOM::default(); 10])
            .unwrap();
        versioned(&mut polymer_ctab);
        let mut polymer_stream = string_stream(&mut heap, "M  STY  1   1 SRU\nM  END\n");
        errors.fill(0);
        assert_eq!(
            MolfileReadPropBlock(
                &mut heap,
                &mut polymer_ctab,
                &mut header,
                Some(&mut polymer_stream),
                1,
                0,
                Some(&mut errors),
                0
            ),
            Ok(0)
        );
        assert_eq!(polymer_ctab.sgroups.used, 1, "{}", text(&errors));
        assert_eq!(
            polymer_group(&heap, &polymer_ctab, 1).type_,
            MOL_FMT_M_STY_SRU as i32
        );

        let mut fixed_ctab = MOL_FMT_CTAB {
            n_property_lines: 2,
            n_atoms: 1,
            atoms: control_atoms,
            ..MOL_FMT_CTAB::default()
        };
        let mut fixed_stream = string_stream(&mut heap, "M  END\nM  REG 77\ntrailing\n");
        header.internal_regno = 0;
        assert_eq!(
            MolfileReadPropBlock(
                &mut heap,
                &mut fixed_ctab,
                &mut header,
                Some(&mut fixed_stream),
                1,
                0,
                None,
                0
            ),
            Ok(0)
        );
        assert_eq!(header.internal_regno, 77);
        assert!(stream.s.nPtr > 0);

        let mut sdf_ctab = MOL_FMT_CTAB::default();
        versioned(&mut sdf_ctab);
        let mut sdf_stream = string_stream(&mut heap, "$$$$\n");
        assert_eq!(
            MolfileReadPropBlock(
                &mut heap,
                &mut sdf_ctab,
                &mut header,
                Some(&mut sdf_stream),
                1,
                0,
                None,
                0
            ),
            Ok(-4)
        );

        let mut long_ctab = MOL_FMT_CTAB::default();
        versioned(&mut long_ctab);
        let mut long_stream = string_stream(&mut heap, &("X".repeat(201) + "\nM  END\n"));
        errors.fill(0);
        assert_eq!(
            MolfileReadPropBlock(
                &mut heap,
                &mut long_ctab,
                &mut header,
                Some(&mut long_stream),
                1,
                0,
                Some(&mut errors),
                0
            ),
            Ok(3)
        );
        assert_eq!(text(&errors), "Too long properties block line");

        let mut eof_ctab = MOL_FMT_CTAB::default();
        versioned(&mut eof_ctab);
        let mut eof_stream = string_stream(&mut heap, "");
        errors.fill(0);
        assert_eq!(
            MolfileReadPropBlock(
                &mut heap,
                &mut eof_ctab,
                &mut header,
                Some(&mut eof_stream),
                1,
                0,
                Some(&mut errors),
                0
            ),
            Ok(2)
        );
        assert_eq!(text(&errors), "Cannot read properties block line");
        errors.fill(0);
        let mut prior_stream = string_stream(&mut heap, "");
        assert_eq!(
            MolfileReadPropBlock(
                &mut heap,
                &mut eof_ctab,
                &mut header,
                Some(&mut prior_stream),
                1,
                9,
                Some(&mut errors),
                0
            ),
            Ok(9)
        );
        assert_eq!(text(&errors), "");
    }

    #[test]
    fn source_port__mol_fmt1__molfilereadsgroupofpolymer__line_1425() {
        let mut heap = SourceHeap::default();
        let mut ctab = polymer_ctab(&mut heap);
        assert_eq!(
            polymer_parse(
                &mut heap,
                &mut ctab,
                "STY",
                "  6 1 SRU 2 MON 3 COP 4 MOD 5 CRO 6 MER",
                41,
            ),
            0
        );
        assert_eq!(ctab.sgroups.used, 6);
        for (id, expected_type) in [
            (1, MOL_FMT_M_STY_SRU as i32),
            (2, MOL_FMT_M_STY_MON as i32),
            (3, MOL_FMT_M_STY_COP as i32),
            (4, MOL_FMT_M_STY_MOD as i32),
            (5, MOL_FMT_M_STY_CRO as i32),
            (6, MOL_FMT_M_STY_MER as i32),
        ] {
            let group = polymer_group(&heap, &ctab, id);
            assert_eq!(group.type_, expected_type);
            assert_eq!(group.xbr1, [-777777.777; 4]);
            assert_eq!(group.xbr2, [-777777.777; 4]);
            assert_eq!(group.smt[0], 0);
        }
        assert_eq!(
            polymer_parse(&mut heap, &mut ctab, "STY", "  1 1 MON", 0),
            0
        );
        assert_eq!(ctab.sgroups.used, 6);
        assert_eq!(
            polymer_group(&heap, &ctab, 1).type_,
            MOL_FMT_M_STY_MON as i32
        );

        assert_eq!(
            polymer_parse(
                &mut heap,
                &mut ctab,
                "SST",
                "  4 1 ALT 2 RAN 3 BLO 4 BLK",
                0,
            ),
            0
        );
        assert_eq!(
            polymer_group(&heap, &ctab, 1).subtype,
            MOL_FMT_M_SST_ALT as i32
        );
        assert_eq!(
            polymer_group(&heap, &ctab, 2).subtype,
            MOL_FMT_M_SST_RAN as i32
        );
        assert_eq!(
            polymer_group(&heap, &ctab, 3).subtype,
            MOL_FMT_M_SST_BLK as i32
        );
        assert_eq!(
            polymer_group(&heap, &ctab, 4).subtype,
            MOL_FMT_M_SST_BLK as i32
        );

        assert_eq!(
            polymer_parse(&mut heap, &mut ctab, "SLB", "  2 1 -12 2 32767", 0),
            0
        );
        assert_eq!(polymer_group(&heap, &ctab, 1).label, -12);
        assert_eq!(polymer_group(&heap, &ctab, 2).label, 32767);

        assert_eq!(
            polymer_parse(&mut heap, &mut ctab, "SCN", "  3 1 HT  2 HH  3 EU ", 0),
            0
        );
        assert_eq!(
            polymer_group(&heap, &ctab, 1).conn,
            MOL_FMT_M_CONN_HT as i32
        );
        assert_eq!(
            polymer_group(&heap, &ctab, 2).conn,
            MOL_FMT_M_CONN_HH as i32
        );
        assert_eq!(
            polymer_group(&heap, &ctab, 3).conn,
            MOL_FMT_M_CONN_EU as i32
        );

        assert_eq!(
            polymer_parse(&mut heap, &mut ctab, "SAL", "   1  3 1 2 10", 0),
            0
        );
        assert_eq!(
            polymer_parse(&mut heap, &mut ctab, "SAL", "   1  2 3 4", 0),
            0
        );
        let group = polymer_group(&heap, &ctab, 1);
        assert_eq!(group.alist.used, 5);
        assert_eq!(
            &heap.slice(group.alist.item.as_const()).unwrap()[..5],
            &[1, 2, 10, 3, 4]
        );

        assert_eq!(
            polymer_parse(&mut heap, &mut ctab, "SBL", "   1  3 1 4 8", 0),
            0
        );
        let group = polymer_group(&heap, &ctab, 1);
        assert_eq!(group.blist.used, 3);
        assert_eq!(
            &heap.slice(group.blist.item.as_const()).unwrap()[..3],
            &[1, 4, 8]
        );

        assert_eq!(
            polymer_parse(&mut heap, &mut ctab, "SDI", "   1  4 1.25 -2.5 3 4", 0),
            0
        );
        assert_eq!(polymer_group(&heap, &ctab, 1).xbr1, [1.25, -2.5, 3.0, 4.0]);
        assert_eq!(
            polymer_parse(&mut heap, &mut ctab, "SDI", "   1  4 5 6 7 8", 0),
            0
        );
        assert_eq!(polymer_group(&heap, &ctab, 1).xbr2, [5.0, 6.0, 7.0, 8.0]);

        assert_eq!(
            polymer_parse(&mut heap, &mut ctab, "SMT", "   1   n value   ", 0),
            0
        );
        assert_eq!(text(&polymer_group(&heap, &ctab, 1).smt), "n value");
        assert_eq!(
            MolFmtSgroups_Free(&mut heap, Some(&mut ctab.sgroups)),
            Ok(())
        );
        assert_eq!(heap.live_source_allocation_count(), 0);

        for (type_, fields, expected) in [
            ("STY", "  1 1 BAD", 5),
            ("SST", "  1 9 ALT", 6),
            ("SLB", "  1 9 1", 7),
            ("SCN", "  1 1 ZZ ", 8),
            ("SAL", "   1  1 11", 9),
            ("SBL", "   1  1 0", 10),
            ("SDI", "   1  3 1 2 3", 11),
            ("SMT", "   9 value", 11),
        ] {
            let mut case_heap = SourceHeap::default();
            let mut case_ctab = polymer_ctab(&mut case_heap);
            assert_eq!(
                polymer_parse(&mut case_heap, &mut case_ctab, "STY", "  1 1 SRU", 0),
                0
            );
            assert_eq!(
                polymer_parse(&mut case_heap, &mut case_ctab, type_, fields, 99),
                expected,
                "{type_} {fields}"
            );
            assert_eq!(case_ctab.sgroups, MOL_FMT_SGROUPS::default());
            assert_eq!(case_heap.live_source_allocation_count(), 0);
        }

        let mut ignored_heap = SourceHeap::default();
        let mut ignored = polymer_ctab(&mut ignored_heap);
        let before = ignored.sgroups.clone();
        assert_eq!(
            polymer_parse(&mut ignored_heap, &mut ignored, "STY", "  0", 77),
            0
        );
        assert_eq!(ignored.sgroups, before);
        assert_eq!(
            polymer_parse(&mut ignored_heap, &mut ignored, "XYZ", "bad", 77),
            0
        );
        assert_eq!(ignored.sgroups, before);
        assert_eq!(
            MolFmtSgroups_Free(&mut ignored_heap, Some(&mut ignored.sgroups)),
            Ok(())
        );

        let mut failure_heap = SourceHeap::default();
        let mut failure = polymer_ctab(&mut failure_heap);
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            polymer_parse(&mut failure_heap, &mut failure, "STY", "  1 1 SRU", 0),
            5
        );
        assert_eq!(failure.sgroups, MOL_FMT_SGROUPS::default());
        assert_eq!(failure_heap.live_source_allocation_count(), 0);
    }

    #[test]
    fn source_port__mol_fmt1__molfilereadheaderlines__line_451() {
        let mut heap = SourceHeap::default();
        let source = input();
        let mut stream = string_stream(&mut heap, &source);
        let mut header = MOL_FMT_HEADER_BLOCK::default();
        assert_eq!(
            MolfileReadHeaderLines(&mut heap, &mut header, Some(&mut stream), None),
            Ok(0)
        );
        assert_eq!(text(&header.molname), "Molecule Name");
        assert_eq!(text(&header.user_initls), "AB");
        assert_eq!(text(&header.prog_name), "PROGNAME");
        assert_eq!(
            (
                header.month,
                header.day,
                header.year,
                header.hour,
                header.minute
            ),
            (1, 2, 24, 12, 34)
        );
        assert_eq!(text(&header.dim_code), "2D");
        assert_eq!(
            (
                header.scaling_factor1,
                header.scaling_factor2,
                header.energy,
                header.internal_regno
            ),
            (7, 1.25, -2.5, 42)
        );
        assert_eq!(text(&header.line2).len(), 52);
        assert_eq!(text(&header.comment), "Comment text");

        let file = heap
            .allocate(vec![FILE {
                bytes: source.as_bytes().to_vec(),
                ..FILE::default()
            }])
            .unwrap();
        let mut file_stream = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_FILE as i32,
            f: file,
            ..INCHI_IOSTREAM::default()
        };
        let mut file_header = MOL_FMT_HEADER_BLOCK::default();
        assert_eq!(
            MolfileReadHeaderLines(&mut heap, &mut file_header, Some(&mut file_stream), None),
            Ok(0)
        );
        assert_eq!(file_header, header);

        for (partial, expected) in [("", 1), ("name\n", 3), ("name\nline2\n", 7)] {
            let mut partial_stream = string_stream(&mut heap, partial);
            let mut partial_header = MOL_FMT_HEADER_BLOCK::default();
            assert_eq!(
                MolfileReadHeaderLines(
                    &mut heap,
                    &mut partial_header,
                    Some(&mut partial_stream),
                    None
                ),
                Ok(expected)
            );
        }
    }

    #[test]
    fn source_port__mol_fmt1__molfilereadcountsline__line_569() {
        fn counts(version: &str) -> String {
            let fields = [2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 999];
            let fixed = fields
                .iter()
                .map(|value| format!("{value:>3}"))
                .collect::<String>();
            assert_eq!(fixed.len(), 33);
            format!("{fixed} {version}\n")
        }
        fn error_text(buffer: &[i8]) -> String {
            text(buffer)
        }

        let mut heap = SourceHeap::default();
        let mut stream = string_stream(&mut heap, &counts("V2000"));
        let mut ctab = MOL_FMT_CTAB {
            n_atoms: 0x1234_0000,
            n_bonds: 0x5678_0000,
            ..MOL_FMT_CTAB::default()
        };
        let mut errors = [0_i8; 256];
        assert_eq!(
            MolfileReadCountsLine(&mut heap, &mut ctab, Some(&mut stream), Some(&mut errors)),
            Ok(0)
        );
        assert_eq!(ctab.n_atoms, i32::from_ne_bytes([2, 0, 0x34, 0x12]));
        assert_eq!(ctab.n_bonds, i32::from_ne_bytes([1, 0, 0x78, 0x56]));
        assert_eq!(
            (
                ctab.chiral_flag,
                ctab.n_stext_entries,
                ctab.n_property_lines
            ),
            (1, 0, 999)
        );
        assert_eq!(text(&ctab.version_string), "V2000");
        assert!(ctab.v3000.is_null());
        assert_eq!(
            (
                ctab.sgroups.allocated,
                ctab.sgroups.used,
                ctab.sgroups.increment
            ),
            (1, 0, 1)
        );
        assert_eq!(
            heap.slice(ctab.sgroups.group.as_const()).unwrap(),
            &[SourceMutPointer::null()]
        );
        assert_eq!(error_text(&errors), "");

        let mut v3000_stream = string_stream(&mut heap, &counts("V3000"));
        let mut v3000 = MOL_FMT_CTAB::default();
        assert_eq!(
            MolfileReadCountsLine(&mut heap, &mut v3000, Some(&mut v3000_stream), None),
            Ok(0)
        );
        assert!(!v3000.v3000.is_null());
        assert_eq!(
            heap.slice(v3000.v3000.as_const()).unwrap(),
            &[MOL_FMT_v3000::default()]
        );

        let mut empty = string_stream(&mut heap, "");
        let mut empty_ctab = MOL_FMT_CTAB::default();
        errors.fill(0);
        assert_eq!(
            MolfileReadCountsLine(
                &mut heap,
                &mut empty_ctab,
                Some(&mut empty),
                Some(&mut errors)
            ),
            Ok(1)
        );
        assert_eq!(error_text(&errors), "Cannot read counts line");

        let mut malformed = string_stream(&mut heap, " xx\u{1}broken\n");
        let mut malformed_ctab = MOL_FMT_CTAB::default();
        errors.fill(0);
        assert_eq!(
            MolfileReadCountsLine(
                &mut heap,
                &mut malformed_ctab,
                Some(&mut malformed),
                Some(&mut errors)
            ),
            Ok(3)
        );
        assert_eq!(
            error_text(&errors),
            "Cannot interpret counts line:  xx.broken"
        );

        let base = counts("V2000").trim_end().to_owned();
        let long_line = format!("{}{}\n", base, " ".repeat(201 - base.len()));
        assert_eq!(long_line.as_bytes()[200], b' ');
        let mut long_stream = string_stream(&mut heap, &long_line);
        let mut long_ctab = MOL_FMT_CTAB::default();
        errors.fill(0);
        assert_eq!(
            MolfileReadCountsLine(
                &mut heap,
                &mut long_ctab,
                Some(&mut long_stream),
                Some(&mut errors)
            ),
            Ok(0)
        );
        assert_eq!(error_text(&errors), "Too long counts line");

        let mut failure_heap = SourceHeap::default();
        let mut failure_stream = string_stream(&mut failure_heap, &counts("V3000"));
        let mut failure_ctab = MOL_FMT_CTAB::default();
        let mut failure_errors = [0_i8; 256];
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            MolfileReadCountsLine(
                &mut failure_heap,
                &mut failure_ctab,
                Some(&mut failure_stream),
                Some(&mut failure_errors)
            ),
            Ok(-1)
        );
        assert!(failure_ctab.v3000.is_null());
        assert_eq!(error_text(&failure_errors), "Out of RAM");

        let mut ignored_heap = SourceHeap::default();
        let mut ignored_stream = string_stream(&mut ignored_heap, &counts("V2000"));
        let mut ignored_ctab = MOL_FMT_CTAB::default();
        ignored_heap.fail_after_allocations(0);
        assert_eq!(
            MolfileReadCountsLine(
                &mut ignored_heap,
                &mut ignored_ctab,
                Some(&mut ignored_stream),
                None
            ),
            Ok(0)
        );
        assert!(ignored_ctab.sgroups.group.is_null());
    }

    #[test]
    fn source_port__mol_fmt1__molfilereadatomsblock__line_661() {
        fn atom_line(symbol: &str, mass: i8, charge: i8, parity: i8, valence: i8) -> String {
            format!(
                "{:>10.4}{:>10.4}{:>10.4} {:<3}{:>2}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}\n",
                1.25, -2.5, 3.75, symbol, mass, charge, parity, 0, 0, valence, 0, 0, 0, 0, 0, 0
            )
        }

        let mut heap = SourceHeap::default();
        let source = format!(
            "{}{}{}",
            atom_line("CL", 2, 1, 3, 6),
            atom_line("N", -1, 4, 1, 3),
            atom_line("O", 0, 9, 0, 2)
        );
        let mut stream = string_stream(&mut heap, &source);
        let atoms = heap.allocate(vec![MOL_FMT_ATOM::default(); 3]).unwrap();
        let mut coordinate_seed = [[b'?' as i8; 32]; 3];
        coordinate_seed[1][31] = b'!' as i8;
        let coords = heap.allocate(coordinate_seed.to_vec()).unwrap();
        let mut ctab = MOL_FMT_CTAB {
            n_atoms: 3,
            atoms,
            coords,
            ..MOL_FMT_CTAB::default()
        };
        let mut errors = [0_i8; 256];
        assert_eq!(
            MolfileReadAtomsBlock(
                &mut heap,
                &mut ctab,
                Some(&mut stream),
                0,
                Some(&mut errors)
            ),
            Ok(0)
        );
        let parsed = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(
            (parsed[0].fx, parsed[0].fy, parsed[0].fz),
            (1.25, -2.5, 3.75)
        );
        assert_eq!(&parsed[0].symbol[..3], &[b'C' as i8, b'l' as i8, 0]);
        assert_eq!(
            (
                parsed[0].mass_difference,
                parsed[0].charge,
                parsed[0].radical,
                parsed[0].stereo_parity,
                parsed[0].valence
            ),
            (2, 3, 0, 3, 6)
        );
        assert_eq!(
            (parsed[1].charge, parsed[1].radical),
            (0, RADICAL_DOUBLET as i8)
        );
        assert_eq!((parsed[2].charge, parsed[2].radical), (-5, 0));
        let stored_coords = heap.slice(coords.as_const()).unwrap();
        let expected_prefix = format!("{:>10.4}{:>10.4}{:>10.4}", 1.25, -2.5, 3.75);
        assert_eq!(
            stored_coords[0][..30]
                .iter()
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>(),
            expected_prefix.as_bytes()
        );
        assert_eq!(stored_coords[0][30], 0);
        assert_eq!(stored_coords[0][31], b'?' as i8);
        assert_eq!(stored_coords[1][31], b'!' as i8);
        assert_eq!(text(&errors), "");

        let mut eof_stream = string_stream(&mut heap, "");
        let mut eof_ctab = MOL_FMT_CTAB {
            n_atoms: 1,
            ..MOL_FMT_CTAB::default()
        };
        errors.fill(0);
        assert_eq!(
            MolfileReadAtomsBlock(
                &mut heap,
                &mut eof_ctab,
                Some(&mut eof_stream),
                0,
                Some(&mut errors)
            ),
            Ok(2)
        );
        assert_eq!(text(&errors), "Cannot read atom block line");

        let mut malformed_stream = string_stream(&mut heap, "bad\u{1}line\n");
        let malformed_atoms = heap.allocate(vec![MOL_FMT_ATOM::default()]).unwrap();
        let mut malformed_ctab = MOL_FMT_CTAB {
            n_atoms: 1,
            atoms: malformed_atoms,
            ..MOL_FMT_CTAB::default()
        };
        errors.fill(0);
        assert_eq!(
            MolfileReadAtomsBlock(
                &mut heap,
                &mut malformed_ctab,
                Some(&mut malformed_stream),
                0,
                Some(&mut errors)
            ),
            Ok(4)
        );
        assert_eq!(text(&errors), "Cannot interpret atom block line: bad.line");

        let mut end_stream = string_stream(&mut heap, "$$$$\n");
        let end_atoms = heap.allocate(vec![MOL_FMT_ATOM::default()]).unwrap();
        let mut end_ctab = MOL_FMT_CTAB {
            n_atoms: 1,
            atoms: end_atoms,
            ..MOL_FMT_CTAB::default()
        };
        errors.fill(0);
        assert_eq!(
            MolfileReadAtomsBlock(
                &mut heap,
                &mut end_ctab,
                Some(&mut end_stream),
                0,
                Some(&mut errors)
            ),
            Ok(-4)
        );

        let mut prior_stream = string_stream(&mut heap, "ignored\n$$$$\n");
        let prior_atoms = heap.allocate(vec![MOL_FMT_ATOM::default(); 2]).unwrap();
        let prior_before = heap.slice(prior_atoms.as_const()).unwrap().to_vec();
        let mut prior_ctab = MOL_FMT_CTAB {
            n_atoms: 2,
            atoms: prior_atoms,
            ..MOL_FMT_CTAB::default()
        };
        errors.fill(0);
        assert_eq!(
            MolfileReadAtomsBlock(
                &mut heap,
                &mut prior_ctab,
                Some(&mut prior_stream),
                7,
                Some(&mut errors)
            ),
            Ok(-7)
        );
        assert_eq!(heap.slice(prior_atoms.as_const()).unwrap(), prior_before);
        assert_eq!(text(&errors), "");

        let valid = atom_line("C", 0, 0, 0, 4);
        let long_line = format!("{}{}\n", valid.trim_end(), "X".repeat(220));
        assert_ne!(long_line.as_bytes()[MOL_FMT_MAXLINELEN as usize], 0);
        let mut long_stream = string_stream(&mut heap, &long_line);
        let long_atoms = heap.allocate(vec![MOL_FMT_ATOM::default()]).unwrap();
        let mut long_ctab = MOL_FMT_CTAB {
            n_atoms: 1,
            atoms: long_atoms,
            ..MOL_FMT_CTAB::default()
        };
        errors.fill(0);
        assert_eq!(
            MolfileReadAtomsBlock(
                &mut heap,
                &mut long_ctab,
                Some(&mut long_stream),
                0,
                Some(&mut errors)
            ),
            Ok(0)
        );
        assert_eq!(text(&errors), "Too long atom block line");

        let mut null_stream = string_stream(&mut heap, &atom_line("C", 0, 0, 0, 4));
        let mut null_ctab = MOL_FMT_CTAB {
            n_atoms: 1,
            ..MOL_FMT_CTAB::default()
        };
        assert_eq!(
            MolfileReadAtomsBlock(&mut heap, &mut null_ctab, Some(&mut null_stream), 0, None),
            Ok(0)
        );
        let _coordinate_type_check: MOL_COORD = [0; 32];
    }

    #[test]
    fn source_port__mol_fmt1__molfilereadbondsblock__line_820() {
        fn bond_line(first: i16, second: i16, bond_type: i8, stereo: i8) -> String {
            format!(
                "{first:>3}{second:>3}{bond_type:>3}{stereo:>3}{:>3}{:>3}\n",
                0, 0
            )
        }

        let mut heap = SourceHeap::default();
        let source = format!("{}{}", bond_line(1, 2, 1, 0), bond_line(2, 3, 4, 6));
        let mut stream = string_stream(&mut heap, &source);
        let bonds = heap.allocate(vec![MOL_FMT_BOND::default(); 2]).unwrap();
        let mut ctab = MOL_FMT_CTAB {
            n_bonds: 2,
            bonds,
            ..MOL_FMT_CTAB::default()
        };
        let mut errors = [0_i8; 256];
        assert_eq!(
            MolfileReadBondsBlock(
                &mut heap,
                &mut ctab,
                Some(&mut stream),
                0,
                Some(&mut errors)
            ),
            Ok(0)
        );
        assert_eq!(
            heap.slice(bonds.as_const()).unwrap(),
            &[
                MOL_FMT_BOND {
                    atnum1: 1,
                    atnum2: 2,
                    bond_type: 1,
                    bond_stereo: 0
                },
                MOL_FMT_BOND {
                    atnum1: 2,
                    atnum2: 3,
                    bond_type: 4,
                    bond_stereo: 6
                },
            ]
        );
        assert_eq!(text(&errors), "");

        let mut eof_stream = string_stream(&mut heap, "");
        let mut eof_ctab = MOL_FMT_CTAB {
            n_bonds: 1,
            ..MOL_FMT_CTAB::default()
        };
        errors.fill(0);
        assert_eq!(
            MolfileReadBondsBlock(
                &mut heap,
                &mut eof_ctab,
                Some(&mut eof_stream),
                0,
                Some(&mut errors)
            ),
            Ok(2)
        );
        assert_eq!(text(&errors), "Cannot read bond block line");

        let malformed_bonds = heap.allocate(vec![MOL_FMT_BOND::default()]).unwrap();
        let mut malformed_ctab = MOL_FMT_CTAB {
            n_bonds: 1,
            bonds: malformed_bonds,
            ..MOL_FMT_CTAB::default()
        };
        let mut malformed_stream = string_stream(&mut heap, "  1 xx\u{1}bad\n");
        errors.fill(0);
        assert_eq!(
            MolfileReadBondsBlock(
                &mut heap,
                &mut malformed_ctab,
                Some(&mut malformed_stream),
                0,
                Some(&mut errors)
            ),
            Ok(4)
        );
        assert_eq!(
            text(&errors),
            "Cannot interpret bond block line:   1 xx.bad"
        );
        assert_eq!(heap.slice(malformed_bonds.as_const()).unwrap()[0].atnum1, 1);

        let valid = bond_line(1, 2, 1, 0);
        let long_line = format!("{}{}\n", valid.trim_end(), "X".repeat(220));
        let mut long_stream = string_stream(&mut heap, &long_line);
        let long_bonds = heap.allocate(vec![MOL_FMT_BOND::default()]).unwrap();
        let mut long_ctab = MOL_FMT_CTAB {
            n_bonds: 1,
            bonds: long_bonds,
            ..MOL_FMT_CTAB::default()
        };
        errors.fill(0);
        assert_eq!(
            MolfileReadBondsBlock(
                &mut heap,
                &mut long_ctab,
                Some(&mut long_stream),
                0,
                Some(&mut errors)
            ),
            Ok(3)
        );
        assert_eq!(
            heap.slice(long_bonds.as_const()).unwrap()[0],
            MOL_FMT_BOND::default()
        );
        assert_eq!(text(&errors), "");

        let prior_bonds = heap.allocate(vec![MOL_FMT_BOND::default(); 2]).unwrap();
        let mut prior_ctab = MOL_FMT_CTAB {
            n_bonds: 2,
            bonds: prior_bonds,
            ..MOL_FMT_CTAB::default()
        };
        let mut prior_stream = string_stream(&mut heap, "ignored\n$$$$\n");
        assert_eq!(
            MolfileReadBondsBlock(&mut heap, &mut prior_ctab, Some(&mut prior_stream), 9, None),
            Ok(-9)
        );

        let mut null_stream = string_stream(&mut heap, &bond_line(1, 2, 1, 0));
        let mut null_ctab = MOL_FMT_CTAB {
            n_bonds: 1,
            ..MOL_FMT_CTAB::default()
        };
        assert_eq!(
            MolfileReadBondsBlock(&mut heap, &mut null_ctab, Some(&mut null_stream), 0, None),
            Ok(0)
        );
    }

    #[test]
    fn source_port__mol_fmt1__molfilereadstextblock__line_917() {
        let mut heap = SourceHeap::default();
        let consumed = "first\n$$$$\n";
        let source = format!("{consumed}{}tail\n", "X".repeat(240) + "\nshort\n");
        let mut stream = string_stream(&mut heap, &source);
        let ctab = MOL_FMT_CTAB {
            n_stext_entries: 2,
            ..MOL_FMT_CTAB::default()
        };
        let mut errors = [0_i8; 256];
        assert_eq!(
            MolfileReadSTextBlock(&mut heap, &ctab, Some(&mut stream), 0, Some(&mut errors)),
            Ok(0)
        );
        let expected_consumed = consumed.len() + 241 + "short\n".len();
        assert_eq!(stream.s.nPtr, expected_consumed as i32);
        assert_eq!(text(&errors), "");

        let mut eof_stream = string_stream(&mut heap, "one\n");
        let eof_ctab = MOL_FMT_CTAB {
            n_stext_entries: 1,
            ..MOL_FMT_CTAB::default()
        };
        assert_eq!(
            MolfileReadSTextBlock(
                &mut heap,
                &eof_ctab,
                Some(&mut eof_stream),
                0,
                Some(&mut errors)
            ),
            Ok(2)
        );
        assert_eq!(text(&errors), "Cannot read STEXT block line");

        let mut prior_stream = string_stream(&mut heap, "");
        errors.fill(0);
        assert_eq!(
            MolfileReadSTextBlock(
                &mut heap,
                &eof_ctab,
                Some(&mut prior_stream),
                7,
                Some(&mut errors)
            ),
            Ok(7)
        );
        assert_eq!(text(&errors), "");

        let mut negative_stream = string_stream(&mut heap, "untouched\n");
        let negative_ctab = MOL_FMT_CTAB {
            n_stext_entries: -1,
            ..MOL_FMT_CTAB::default()
        };
        assert_eq!(
            MolfileReadSTextBlock(
                &mut heap,
                &negative_ctab,
                Some(&mut negative_stream),
                0,
                None
            ),
            Ok(0)
        );
        assert_eq!(negative_stream.s.nPtr, 0);
    }

    #[test]
    fn source_port__mol_fmt1__molfilereaddatalines__line_165() {
        fn counts(atoms: i32, bonds: i32, version: &str) -> String {
            let fields = [atoms, bonds, 0, 0, 1, 0, 0, 0, 0, 0, 999];
            let fixed = fields
                .iter()
                .map(|value| format!("{value:>3}"))
                .collect::<String>();
            assert_eq!(fixed.len(), 33);
            format!("{fixed} {version}\n")
        }

        fn atom_line(symbol: &str, mass: i8, valence: i8, x: f64) -> String {
            format!(
                "{:>10.4}{:>10.4}{:>10.4} {:<3}{:>2}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}\n",
                x, -2.5, 3.75, symbol, mass, 0, 0, 0, 0, valence, 0, 0, 0, 0, 0, 0
            )
        }

        fn bond_line(first: i16, second: i16, bond_type: i8, stereo: i8) -> String {
            format!(
                "{first:>3}{second:>3}{bond_type:>3}{stereo:>3}{:>3}{:>3}\n",
                0, 0
            )
        }

        fn v2000(atom_one: (&str, i8, i8), atom_two: Option<(&str, i8, i8)>) -> String {
            let atom_count = if atom_two.is_some() { 2 } else { 1 };
            let bond_count = if atom_two.is_some() { 1 } else { 0 };
            let mut source = format!(
                "{}{}{}",
                input(),
                counts(atom_count, bond_count, "V2000"),
                atom_line(atom_one.0, atom_one.1, atom_one.2, 1.25)
            );
            if let Some(atom) = atom_two {
                source.push_str(&atom_line(atom.0, atom.1, atom.2, 4.5));
                source.push_str(&bond_line(1, 2, 2, 1));
            }
            source.push_str("M  END\n");
            source
        }

        let mut heap = SourceHeap::default();
        let source = v2000(("C", 0, 4), Some(("O", 0, 2)));
        let mut stream = string_stream(&mut heap, &source);
        let mut err = -1;
        let mut errors = [0_i8; 256];
        let data = MolfileReadDataLines(
            &mut heap,
            Some(&mut stream),
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            1,
            0,
            &mut err,
            Some(&mut errors),
            0,
        )
        .unwrap();
        assert!(!data.is_null());
        assert_eq!(err, 0);
        assert_eq!(text(&errors), "");
        assert_eq!(stream.s.nPtr, source.len() as i32);
        {
            let parsed = &heap.slice(data.as_const()).unwrap()[0];
            assert_eq!(text(&parsed.hdr.molname), "Molecule Name");
            assert_eq!((parsed.ctab.n_atoms, parsed.ctab.n_bonds), (2, 1));
            assert!(!parsed.ctab.atoms.is_null());
            assert!(!parsed.ctab.bonds.is_null());
            assert!(!parsed.ctab.coords.is_null());
            let atoms = heap.slice(parsed.ctab.atoms.as_const()).unwrap();
            assert_eq!(&atoms[0].symbol[..2], &[b'C' as i8, 0]);
            assert_eq!((atoms[0].fx, atoms[0].fy, atoms[0].fz), (1.25, -2.5, 3.75));
            assert_eq!(&atoms[1].symbol[..2], &[b'O' as i8, 0]);
            let bonds = heap.slice(parsed.ctab.bonds.as_const()).unwrap();
            assert_eq!((bonds[0].atnum1, bonds[0].atnum2), (1, 2));
            assert_eq!((bonds[0].bond_type, bonds[0].bond_stereo), (2, 1));
            let coords = heap.slice(parsed.ctab.coords.as_const()).unwrap();
            assert_eq!(coords.len(), 2);
            assert_eq!(coords[0][30], 0);
        }
        assert!(FreeMolfileData(&mut heap, data).unwrap().is_null());
        assert_eq!(heap.live_source_allocation_count(), 1);

        let mut header_heap = SourceHeap::default();
        let header_source = v2000(("C", 0, 4), Some(("N", 0, 3)));
        let mut header_stream = string_stream(&mut header_heap, &header_source);
        let header = header_heap
            .allocate_model_storage(vec![MOL_FMT_HEADER_BLOCK::default()])
            .unwrap();
        let ctab = header_heap
            .allocate_model_storage(vec![MOL_FMT_CTAB::default()])
            .unwrap();
        let mut header_err = -1;
        let result = MolfileReadDataLines(
            &mut header_heap,
            Some(&mut header_stream),
            header,
            ctab,
            1,
            0,
            &mut header_err,
            None,
            1,
        )
        .unwrap();
        assert_eq!(result, header.cast::<MOL_FMT_DATA>());
        assert_eq!(header_err, 0);
        assert_eq!(
            text(&header_heap.slice(header.as_const()).unwrap()[0].molname),
            "Molecule Name"
        );
        let parsed_ctab = header_heap.slice(ctab.as_const()).unwrap()[0].clone();
        assert_eq!((parsed_ctab.n_atoms, parsed_ctab.n_bonds), (2, 1));
        assert!(parsed_ctab.atoms.is_null());
        assert!(parsed_ctab.bonds.is_null());
        assert!(parsed_ctab.coords.is_null());
        let mut sgroups = parsed_ctab.sgroups;
        MolFmtSgroups_Free(&mut header_heap, Some(&mut sgroups)).unwrap();

        let mut no_coords_heap = SourceHeap::default();
        let no_coords_source = v2000(("C", 0, 4), None);
        let mut no_coords_stream = string_stream(&mut no_coords_heap, &no_coords_source);
        let mut no_coords_err = -1;
        let no_coords = MolfileReadDataLines(
            &mut no_coords_heap,
            Some(&mut no_coords_stream),
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            0,
            0,
            &mut no_coords_err,
            None,
            0,
        )
        .unwrap();
        assert_eq!(no_coords_err, 0);
        assert!(
            no_coords_heap.slice(no_coords.as_const()).unwrap()[0]
                .ctab
                .coords
                .is_null()
        );
        assert!(
            FreeMolfileData(&mut no_coords_heap, no_coords)
                .unwrap()
                .is_null()
        );

        let v3000_source = format!(
            "{}{}M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 1 -2.5 3.75 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n",
            input(),
            counts(0, 0, "V3000")
        );
        let mut v3000_heap = SourceHeap::default();
        let mut v3000_stream = string_stream(&mut v3000_heap, &v3000_source);
        let mut v3000_err = -1;
        let v3000 = MolfileReadDataLines(
            &mut v3000_heap,
            Some(&mut v3000_stream),
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            1,
            0,
            &mut v3000_err,
            None,
            0,
        )
        .unwrap();
        assert_eq!(v3000_err, 0);
        let v3000_ctab = &v3000_heap.slice(v3000.as_const()).unwrap()[0].ctab;
        assert_eq!((v3000_ctab.n_atoms, v3000_ctab.n_bonds), (1, 0));
        assert!(!v3000_ctab.v3000.is_null());
        assert_eq!(v3000_stream.s.nPtr, v3000_source.len() as i32);
        assert!(FreeMolfileData(&mut v3000_heap, v3000).unwrap().is_null());

        let mut eof_heap = SourceHeap::default();
        let mut eof_stream = string_stream(&mut eof_heap, "");
        let mut eof_err = -1;
        let eof = MolfileReadDataLines(
            &mut eof_heap,
            Some(&mut eof_stream),
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            1,
            0,
            &mut eof_err,
            None,
            0,
        )
        .unwrap();
        assert!(eof.is_null());
        assert_eq!(eof_err, 11);
        assert_eq!(eof_heap.live_source_allocation_count(), 1);

        let ended_source = format!("{}{}$$$$\n", input(), counts(1, 0, "V2000"));
        let mut ended_heap = SourceHeap::default();
        let mut ended_stream = string_stream(&mut ended_heap, &ended_source);
        let mut ended_err = -1;
        let mut ended_errors = [0_i8; 256];
        let ended = MolfileReadDataLines(
            &mut ended_heap,
            Some(&mut ended_stream),
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            1,
            0,
            &mut ended_err,
            Some(&mut ended_errors),
            0,
        )
        .unwrap();
        assert!(ended.is_null());
        assert_eq!(ended_err, -34);
        assert!(text(&ended_errors).contains("Cannot interpret atom block line:"));
        assert_eq!(ended_heap.live_source_allocation_count(), 1);

        for (source, expected_err, expected_text) in [
            (
                v2000(("C", 0, 21), None),
                79,
                "Too large input atomic valence",
            ),
            (
                v2000(("H", 3, 1), None),
                78,
                "Unacceptable isotope of hydrogen",
            ),
        ] {
            let mut case_heap = SourceHeap::default();
            let mut case_stream = string_stream(&mut case_heap, &source);
            let mut case_err = -1;
            let mut case_errors = [0_i8; 256];
            let result = MolfileReadDataLines(
                &mut case_heap,
                Some(&mut case_stream),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
                0,
                &mut case_err,
                Some(&mut case_errors),
                0,
            )
            .unwrap();
            assert!(result.is_null());
            assert_eq!(case_err, expected_err);
            assert_eq!(text(&case_errors), expected_text);
            assert_eq!(case_heap.live_source_allocation_count(), 1);
        }

        let mut allocation_heap = SourceHeap::default();
        let allocation_source = v2000(("C", 0, 4), None);
        let mut allocation_stream = string_stream(&mut allocation_heap, &allocation_source);
        allocation_heap.fail_after_allocations(0);
        let mut allocation_err = -1;
        let mut allocation_errors = [0_i8; 256];
        let allocation = MolfileReadDataLines(
            &mut allocation_heap,
            Some(&mut allocation_stream),
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            1,
            0,
            &mut allocation_err,
            Some(&mut allocation_errors),
            0,
        )
        .unwrap();
        assert!(allocation.is_null());
        assert_eq!(allocation_err, 1);
        assert_eq!(text(&allocation_errors), "Out of RAM");
        assert_eq!(allocation_heap.live_source_allocation_count(), 1);
    }

    #[test]
    fn source_port__mol_fmt1__readmolfile__line_92() {
        fn counts(atoms: i32, bonds: i32) -> String {
            let fields = [atoms, bonds, 0, 0, 1, 0, 0, 0, 0, 0, 999];
            let fixed = fields
                .iter()
                .map(|value| format!("{value:>3}"))
                .collect::<String>();
            format!("{fixed} V2000\n")
        }

        fn atom_line(symbol: &str) -> String {
            format!(
                "{:>10.4}{:>10.4}{:>10.4} {:<3}{:>2}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}\n",
                1.25, -2.5, 3.75, symbol, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0
            )
        }

        fn molfile(symbol: &str) -> String {
            format!("{}{}{}M  END\n", input(), counts(1, 0), atom_line(symbol))
        }

        fn buffer(heap: &mut SourceHeap, capacity: usize, initial: &[u8]) -> SourceMutPointer<i8> {
            let mut bytes = vec![0_i8; capacity];
            for (target, source) in bytes.iter_mut().zip(initial) {
                *target = *source as i8;
            }
            heap.allocate_model_storage(bytes).unwrap()
        }

        let mut heap = SourceHeap::default();
        let source = format!(
            "{}>  <CAS>\n12-34\n\n>  <X>\n  Alpha\t\tBeta  \n\n$$$$\nNEXT\n",
            molfile("C")
        );
        let mut stream = string_stream(&mut heap, &source);
        let name = buffer(&mut heap, 16, b"old\0");
        let label = buffer(&mut heap, 2, b"x\0").as_const();
        let value = buffer(&mut heap, 256, b"old\0");
        let mut id = 99_u64;
        let mut err = -1_i32;
        let mut errors = [0_i8; 256];
        let data = ReadMolfile(
            &mut heap,
            Some(&mut stream),
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            1,
            POLYMERS_NO as i32,
            0,
            name,
            16,
            Some(&mut id),
            label,
            value,
            &mut err,
            Some(&mut errors),
            0,
        )
        .unwrap();
        assert!(!data.is_null());
        assert_eq!(err, 0);
        assert_eq!(id, 1_234);
        assert_eq!(text(heap.slice(name.as_const()).unwrap()), "");
        assert_eq!(text(heap.slice(value.as_const()).unwrap()), "Alpha Beta");
        assert_eq!(text(&errors), "");
        assert_eq!(
            &heap.slice(stream.s.pStr.as_const()).unwrap()[stream.s.nPtr as usize..][..5],
            &[b'N' as i8, b'E' as i8, b'X' as i8, b'T' as i8, b'\n' as i8]
        );
        let parsed_ctab = &heap.slice(data.as_const()).unwrap()[0].ctab;
        assert_eq!(
            text(&heap.slice(parsed_ctab.atoms.as_const()).unwrap()[0].symbol),
            "C"
        );
        assert!(FreeMolfileData(&mut heap, data).unwrap().is_null());

        for (symbol, treat_polymers, treat_npzz, expected_symbol, expected_err, expected_text) in [
            ("*", POLYMERS_NO as i32, 1, "Zz", 0, ""),
            ("*", 1, 0, "Zz", 0, ""),
            ("*", POLYMERS_NO as i32, 0, "*", 76, "Invalid element(s): *"),
            (
                "Zy",
                POLYMERS_NO as i32,
                1,
                "Zy",
                76,
                "Invalid element(s): Zy",
            ),
        ] {
            let mut case_heap = SourceHeap::default();
            let case_source = format!("{}$$$$\n", molfile(symbol));
            let mut case_stream = string_stream(&mut case_heap, &case_source);
            let mut case_err = -1;
            let mut case_errors = [0_i8; 256];
            let parsed = ReadMolfile(
                &mut case_heap,
                Some(&mut case_stream),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
                treat_polymers,
                treat_npzz,
                SourceMutPointer::null(),
                0,
                None,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                &mut case_err,
                Some(&mut case_errors),
                0,
            )
            .unwrap();
            assert!(!parsed.is_null(), "{symbol}");
            let ctab = &case_heap.slice(parsed.as_const()).unwrap()[0].ctab;
            let parsed_symbol = text(&case_heap.slice(ctab.atoms.as_const()).unwrap()[0].symbol);
            assert_eq!(parsed_symbol, expected_symbol, "{symbol}");
            assert_eq!(case_err, expected_err, "{symbol}");
            assert_eq!(text(&case_errors), expected_text, "{symbol}");
            assert!(FreeMolfileData(&mut case_heap, parsed).unwrap().is_null());
        }

        let ended_source = format!("{}{}$$$$\nNEXT\n", input(), counts(1, 0));
        let mut ended_heap = SourceHeap::default();
        let mut ended_stream = string_stream(&mut ended_heap, &ended_source);
        let mut ended_err = -1;
        let ended = ReadMolfile(
            &mut ended_heap,
            Some(&mut ended_stream),
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            1,
            0,
            0,
            SourceMutPointer::null(),
            0,
            None,
            SourceConstPointer::null(),
            SourceMutPointer::null(),
            &mut ended_err,
            None,
            0,
        )
        .unwrap();
        assert!(ended.is_null());
        assert_eq!(ended_err, 34);
        assert_eq!(
            &ended_heap.slice(ended_stream.s.pStr.as_const()).unwrap()
                [ended_stream.s.nPtr as usize..][..5],
            &[b'N' as i8, b'E' as i8, b'X' as i8, b'T' as i8, b'\n' as i8]
        );

        let mut eof_heap = SourceHeap::default();
        let mut eof_stream = string_stream(&mut eof_heap, "");
        let eof_name = buffer(&mut eof_heap, 8, b"old\0");
        let mut eof_id = 7_u64;
        let mut eof_err = -1;
        let eof = ReadMolfile(
            &mut eof_heap,
            Some(&mut eof_stream),
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            1,
            0,
            0,
            eof_name,
            8,
            Some(&mut eof_id),
            SourceConstPointer::null(),
            SourceMutPointer::null(),
            &mut eof_err,
            None,
            1,
        )
        .unwrap();
        assert!(eof.is_null());
        assert_eq!(eof_err, 11);
        assert_eq!(eof_id, 0);
        assert_eq!(text(eof_heap.slice(eof_name.as_const()).unwrap()), "");
    }

    #[test]
    fn source_port__mol_fmt1__molfiletreatpseudoelementatoms__line_1860() {
        fn atom(symbol: &[u8], tail: i8) -> MOL_FMT_ATOM {
            let mut atom = MOL_FMT_ATOM::default();
            for (destination, source) in atom.symbol.iter_mut().zip(symbol) {
                *destination = *source as i8;
            }
            atom.symbol[5] = tail;
            atom
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate(vec![
                atom(b"C\0", 11),
                atom(b"*\0", 12),
                atom(b"Zz\0", 13),
                atom(b"ZZ\0", 14),
            ])
            .unwrap();
        let mut ctab = MOL_FMT_CTAB {
            n_atoms: 4,
            atoms,
            ..MOL_FMT_CTAB::default()
        };
        let mut err = 0;
        let mut errors = [0_i8; 256];
        assert_eq!(
            MolfileTreatPseudoElementAtoms(&mut heap, &mut ctab, 1, &mut err, Some(&mut errors),),
            Ok(2)
        );
        assert_eq!(err, 0);
        assert_eq!(text(&errors), "");
        let converted = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(&converted[1].symbol[..3], &[b'Z' as i8, b'z' as i8, 0]);
        assert_eq!(converted[1].symbol[5], 12);
        assert_eq!(&converted[2].symbol[..3], &[b'Z' as i8, b'z' as i8, 0]);
        assert_eq!(converted[2].symbol[5], 13);
        assert_eq!(&converted[3].symbol[..3], &[b'Z' as i8, b'Z' as i8, 0]);

        let prohibited_atoms = heap
            .allocate(vec![atom(b"*\0", 21), atom(b"Zz\0", 22), atom(b"Zy\0", 23)])
            .unwrap();
        let mut prohibited = MOL_FMT_CTAB {
            n_atoms: 3,
            atoms: prohibited_atoms,
            ..MOL_FMT_CTAB::default()
        };
        err = 0;
        errors.fill(0);
        assert_eq!(
            MolfileTreatPseudoElementAtoms(
                &mut heap,
                &mut prohibited,
                0,
                &mut err,
                Some(&mut errors),
            ),
            Ok(2)
        );
        assert_eq!(err, 76);
        assert_eq!(text(&errors), "Invalid element(s): *; Zz; Zy");
        let unchanged = heap.slice(prohibited_atoms.as_const()).unwrap();
        assert_eq!(&unchanged[0].symbol[..2], &[b'*' as i8, 0]);
        assert_eq!(&unchanged[1].symbol[..3], &[b'Z' as i8, b'z' as i8, 0]);
        assert_eq!(&unchanged[2].symbol[..3], &[b'Z' as i8, b'y' as i8, 0]);

        let preexisting_atoms = heap.allocate(vec![atom(b"Zy\0", 31)]).unwrap();
        let mut preexisting = MOL_FMT_CTAB {
            n_atoms: 1,
            atoms: preexisting_atoms,
            ..MOL_FMT_CTAB::default()
        };
        err = 91;
        errors.fill(0);
        assert_eq!(
            MolfileTreatPseudoElementAtoms(
                &mut heap,
                &mut preexisting,
                -7,
                &mut err,
                Some(&mut errors),
            ),
            Ok(0)
        );
        assert_eq!(err, 91);
        assert_eq!(text(&errors), "Invalid element(s): Zy");

        let mut empty = MOL_FMT_CTAB {
            n_atoms: 0,
            atoms: SourceMutPointer::null(),
            ..MOL_FMT_CTAB::default()
        };
        err = 0;
        assert_eq!(
            MolfileTreatPseudoElementAtoms(&mut heap, &mut empty, 0, &mut err, None),
            Ok(0)
        );
        assert_eq!(err, 0);

        let mut negative = MOL_FMT_CTAB {
            n_atoms: -1,
            atoms: SourceMutPointer::null(),
            ..MOL_FMT_CTAB::default()
        };
        assert_eq!(
            MolfileTreatPseudoElementAtoms(&mut heap, &mut negative, 0, &mut err, None),
            Ok(0)
        );
    }
}
