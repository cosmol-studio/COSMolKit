use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichinorm::MarkRingSystemsInp;
use crate::source::base::ichisort::iisort;
use crate::source::base::ichister::ReconcileAllCmlBondParities;
use crate::source::base::strutil::{
    CompAtomData_GetNumMapping, DisconnectMetals, DisconnectSalts, MarkDisconnectedComponents,
    bMayDisconnectMetals, bNumHeterAtomHasIsotopicH, fix_odd_things, imat_free, imat_new,
    post_fix_odd_things, remove_ion_pairs, subgraf_free, subgraf_new,
    subgraf_pathfinder_collect_all, subgraf_pathfinder_free, subgraf_pathfinder_new,
    subgraf_pathfinder_run,
};
use crate::source::base::util::{
    detect_unusual_el_valence, inchi_calloc, inchi_free, is_el_a_metal, is_ilist_inside,
    is_in_the_ilist, mystrncpy_slice,
};
use crate::source_types::{
    _IS_ERROR, _IS_FATAL, _IS_OKAY, _IS_WARNING, AT_NUMB, BOND_TAUTOM, CLOSING_SRU_DIRADICAL,
    CLOSING_SRU_HIGHER_ORDER_BOND, CLOSING_SRU_NOT_APPLICABLE, CLOSING_SRU_RING, COMP_ATOM_DATA,
    INCHI_BAS, INCHI_CLOCK, INCHI_MODE, INCHI_REC, INPUT_PARMS, MAX_NUM_STEREO_BONDS, MAXVAL,
    MOL_COORD, NO_POLYMER, OAD_AtProps, OAD_Polymer, OAD_PolymerUnit, OAD_V3000, ORIG_ATOM_DATA,
    POLYMER_CONN_EU, POLYMER_CONN_HT, POLYMER_CONN_NON, POLYMER_REPRESENTATION_MIXED,
    POLYMER_REPRESENTATION_SOURCE_BASED, POLYMER_REPRESENTATION_STRUCTURE_BASED,
    POLYMER_REPRESENTATION_UNRECOGNIZED, POLYMER_SST_ALT, POLYMER_SST_BLK, POLYMER_SST_NON,
    POLYMER_SST_RAN, POLYMER_STY_COP, POLYMER_STY_CRO, POLYMER_STY_MER, POLYMER_STY_MOD,
    POLYMER_STY_MON, POLYMER_STY_SRU, POLYMERS_NO, RADICAL_DOUBLET, RADICAL_SINGLET,
    RADICAL_TRIPLET, SB_PARITY_FLAG, SB_PARITY_MASK, SB_PARITY_SHFT, STRUCT_DATA,
    SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer, TG_FLAG_CHECK_VALENCE_COORD,
    TG_FLAG_DISCONNECT_COORD, TG_FLAG_DISCONNECT_COORD_DONE, TG_FLAG_DISCONNECT_SALTS,
    TG_FLAG_DISCONNECT_SALTS_DONE, TG_FLAG_FIX_ODD_THINGS_DONE, TG_FLAG_FIX_SP3_BUG,
    TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE, TG_FLAG_FOUND_ISOTOPIC_H_DONE, TG_FLAG_RECONNECT_COORD,
    inp_ATOM, tagFrameShifScheme_FSS_NONE, tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE,
    tagINCHIBondType_INCHI_BOND_TYPE_SINGLE, tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE,
};

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OrigAtData_DebugTrace(
    heap: &SourceHeap,
    data: &ORIG_ATOM_DATA,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1108 OrigAtData_DebugTrace
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void OrigAtData_DebugTrace( ORIG_ATOM_DATA* d )
{
    int i, k;

    ITRACE_( "\n\n*********************************************************************\n* ORIG_ATOM_DATA @ 0x%p", d );
    ITRACE_( "\n*  num_inp_atoms = %-d\n*  num_inp_bonds = %-d\n*  num_dimensions = %-d\n*  num_components = %-d",
            d->num_inp_atoms, d->num_inp_bonds, d->num_dimensions, d->num_components );
    ITRACE_( "\n*  ATOMS" );
    for (i = 0; i < d->num_inp_atoms; i++)
    {
        ITRACE_( "\n*    #%-5d %s%-d ( charge %-d, rad %-d nH %-d val %-d) [%-f %-f %-f]",
            i, d->at[i].elname, d->at[i].orig_at_number, d->at[i].charge, d->at[i].radical, d->at[i].num_H, d->at[i].valence,
            d->at[i].x, d->at[i].y, d->at[i].z );
        if (d->at[i].valence > 0)
        {
            ITRACE_( "\n            bonds to     " );
            for (k = 0; k < d->at[i].valence; k++)
            {
                int nbr = d->at[i].neighbor[k]; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
                ITRACE_( "%s%-3d ", d->at[nbr].elname, nbr + 1 );
            }
        }
        if (d->at[i].valence > 0)
        {
            ITRACE_( "\n            bond types   " );
            for (k = 0; k < d->at[i].valence; k++)
                ITRACE_( "%-3d ", d->at[i].bond_type[k] );
        }
    }
    /*OAD_Polymer_DebugTrace( d->polymer );*/
    ITRACE_( "\n* V3000 INFO @ 0x%-p", d->v3000 );
    ITRACE_( "\n*\n" );
    if (d->v3000)
    {
        ITRACE_( "\n*  n_star_atoms = %-d\n*  n_haptic_bonds = %-d\n*  n_collections = %-d",
            d->v3000->n_star_atoms, d->v3000->n_haptic_bonds, d->v3000->n_collections );
    }
    ITRACE_( "\n*\n* End ORIG_ATOM_DATA\n*********************************************************************\n" );

    return;
}
    */
    // END INCHI C FUNCTION: OrigAtData_DebugTrace
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OrigAtData_DebugTrace
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️✔️: ITRACE_ expands to 0 && _inchi_trace, so trace arguments are not evaluated.
    // INCHI✔️✔️: Loop conditions, valence tests, and the nbr initializer remain active C behavior.
    // END INCHI ACTIVE MACRO CONFIGURATION: OrigAtData_DebugTrace

    let atom_count = usize::try_from(data.num_inp_atoms.max(0))
        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    for atom_index in 0..atom_count {
        let atom = heap
            .slice(data.at.as_const())?
            .get(atom_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atom.valence > 0 {
            let valence = usize::try_from(atom.valence)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            for bond_index in 0..valence {
                let _neighbor = atom
                    .neighbor
                    .get(bond_index)
                    .copied()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
            }
        }
        if atom.valence > 0 {
            for _ in 0..atom.valence {
                // The active ITRACE_ right operand is unevaluated.
            }
        }
    }
    if !data.v3000.is_null() {
        // The active ITRACE_ right operand is unevaluated.
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_DebugTrace(
    heap: &SourceHeap,
    polymer: SourceMutPointer<OAD_Polymer>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3757 OAD_Polymer_DebugTrace
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void OAD_Polymer_DebugTrace( OAD_Polymer *p )
{
    int i;

    if (!p)
    {
        return;
    }

    ITRACE_( "\n\n* POLYMER INFO @ %-p (%-d group(s))", p, p->n );
    ITRACE_( "\n\n* %-d star atoms: ", p->n_pzz );
    for (i = 0; i < p->n_pzz; i++)
    {
        ITRACE_( " %-d", p->pzz[i] );
    }

    for (i = 0; i < p->n; i++)
    {
        ITRACE_( "\n* Polymer unit %-d", i );
        OAD_PolymerUnit_DebugTrace( p->units[i] );
    }
    ITRACE_( "\n* Really-do-PS = %-d", p->really_do_frame_shift );
    ITRACE_( "\n* Frame_shift_scheme = %-d", p->frame_shift_scheme );
    ITRACE_( "\n* Edit-repeats = %-d", p->edit_repeats );
    ITRACE_( "\n* End POLYMER INFO\n" );
    return;

}
    */
    // END INCHI C FUNCTION: OAD_Polymer_DebugTrace
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_Polymer_DebugTrace
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️✔️: ITRACE_ expands to 0 && _inchi_trace, so pzz and formatting arguments are not evaluated.
    // INCHI✔️✔️: The p->n loop and p->units[i] evaluation remain active and invoke the source-backed unit helper.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_Polymer_DebugTrace

    if polymer.is_null() {
        return Ok(());
    }
    let polymer = heap
        .slice(polymer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    // The n_pzz loop body is the unevaluated right operand of 0 && _inchi_trace.
    if polymer.n > 0 {
        let unit_count = usize::try_from(polymer.n)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let unit_pointers = heap.slice(polymer.units.as_const())?;
        for unit_index in 0..unit_count {
            let unit_pointer = *unit_pointers
                .get(unit_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if unit_pointer.is_null() {
                OAD_PolymerUnit_DebugTrace(None);
            } else {
                let unit = heap
                    .slice(unit_pointer.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if unit.na == i32::MIN {
                    return Err(SourceHeapError::UnsupportedSourceBehavior);
                }
                OAD_PolymerUnit_DebugTrace(Some(unit));
            }
        }
    }
    Ok(())
}

const SOURCE_SIZEOF_INP_ATOM: u64 = 176;
const SOURCE_SIZEOF_OAD_POLYMER: u64 = 48;
const SOURCE_SIZEOF_OAD_V3000: u64 = 104;

fn add_preprocess_message(
    structure_data: &mut STRUCT_DATA,
    message: &[u8],
) -> Result<(), SourceHeapError> {
    let mut c_message = message.iter().map(|byte| *byte as i8).collect::<Vec<_>>();
    c_message.push(0);
    AddErrorMessage(Some(&mut structure_data.pStrErrStruct), Some(&c_message))?;
    Ok(())
}

fn mask_connected_parities(
    heap: &mut SourceHeap,
    data: &ORIG_ATOM_DATA,
) -> Result<(), SourceHeapError> {
    if data.at.is_null() {
        if data.num_inp_atoms <= 0 {
            return Ok(());
        }
        return Err(SourceHeapError::MissingAllocation);
    }
    let atom_count =
        usize::try_from(data.num_inp_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = heap
        .slice_mut(data.at)?
        .get_mut(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for atom in atoms {
        for parity in atom
            .sb_parity
            .iter_mut()
            .take(MAX_NUM_STEREO_BONDS as usize)
            .take_while(|parity| **parity != 0)
        {
            *parity &= SB_PARITY_MASK as i8;
        }
    }
    Ok(())
}

fn restore_disconnected_parities(
    heap: &mut SourceHeap,
    data: &ORIG_ATOM_DATA,
) -> Result<(), SourceHeapError> {
    let atom_count =
        usize::try_from(data.num_inp_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = heap
        .slice_mut(data.at)?
        .get_mut(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for atom in atoms {
        for parity in atom
            .sb_parity
            .iter_mut()
            .take(MAX_NUM_STEREO_BONDS as usize)
            .take_while(|parity| **parity != 0)
        {
            let value = *parity;
            if value & SB_PARITY_FLAG as i8 != 0 {
                *parity = (value >> SB_PARITY_SHFT) & SB_PARITY_MASK as i8;
            }
        }
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn PreprocessOneStructure(
    heap: &mut SourceHeap,
    _clock: Option<&mut INCHI_CLOCK>,
    structure_data: &mut STRUCT_DATA,
    input_parameters: &INPUT_PARMS,
    original: &mut ORIG_ATOM_DATA,
    prepared: &mut [ORIG_ATOM_DATA],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:488 PreprocessOneStructure
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: PreprocessOneStructure
    // INCHI✔️❌: int PreprocessOneStructure( struct tagINCHI_CLOCK *ic,
    // INCHI✔️❌:                             STRUCT_DATA           *sd,
    // INCHI✔️❌:                             INPUT_PARMS           *ip,
    // INCHI✔️❌:                             ORIG_ATOM_DATA        *orig_inp_data,
    // INCHI✔️❌:                             ORIG_ATOM_DATA        *prep_inp_data )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     INCHI_MODE bTautFlags = 0;
    // INCHI✔️❌:     INCHI_MODE bTautFlagsDone = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* 1. Copy orig_inp_data --> prep_inp_data */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 > OrigAtData_Duplicate( prep_inp_data, orig_inp_data ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         AddErrorMessage( sd->pStrErrStruct, "Out of RAM" );
    // INCHI✔️❌:         sd->nStructReadError = 99;
    // INCHI✔️❌:         sd->nErrorType = _IS_FATAL;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 && (EXTR_HAS_METAL_ATOM & (EXTR_MASK | EXTR_FLAG) ) )
    // INCHI✔️❌:     if (bHasMetalAtom( orig_inp_data ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sd->bExtract |= EXTR_HAS_METAL_ATOM;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /* 2. Fix odd things in prep_inp_data            */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (0 < fix_odd_things( prep_inp_data->num_inp_atoms, prep_inp_data->at, /*0*/ip->bTautFlags & TG_FLAG_FIX_SP3_BUG, ip->bFixNonUniformDraw ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* changed 2010-03-17 DT */
    // INCHI✔️❌:         if (!ip->bNoWarnings)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             WarningMessage( sd->pStrErrStruct, "Charges were rearranged" );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (sd->nErrorType < _IS_WARNING)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->nErrorType = _IS_WARNING;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( FIX_ADJ_RAD == 1 )
    // INCHI✔️❌:     if (ip->bTautFlags & TG_FLAG_FIX_ADJ_RADICALS)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (0 < FixAdjacentRadicals( prep_inp_data->num_inp_atoms, prep_inp_data->at ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ADJ_RADICALS_DONE;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 && (EXTR_FLAGS & EXTR_HAS_FEATURE) )
    // INCHI✔️❌:     if (bFoundFeature( prep_inp_data->at, prep_inp_data->num_inp_atoms ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sd->bExtract |= EXTR_HAS_FEATURE;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Find whether the structure can be disconnected or is a salt */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Needs salt disconnection? */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bTautFlags & TG_FLAG_DISCONNECT_SALTS)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         prep_inp_data->bDisconnectSalts = ( 0 < DisconnectSalts( prep_inp_data, 0 ) );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         prep_inp_data->bDisconnectSalts = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Needs metal disconnection? */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ip->bTautFlags & TG_FLAG_DISCONNECT_COORD)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         i = ( 0 != ( ip->bTautFlags & TG_FLAG_CHECK_VALENCE_COORD ) );
    // INCHI✔️❌:         bMayDisconnectMetals( prep_inp_data, i, &bTautFlagsDone ); /* changes prep_inp_data->bDisconnectCoord */
    // INCHI✔️❌:         sd->bTautFlagsDone[INCHI_BAS] |= bTautFlagsDone; /* whether any disconnection has been rejected because of the metal proper valence */
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 )
    // INCHI✔️❌:         if (i && ( bTautFlagsDone & TG_FLAG_CHECK_VALENCE_COORD_DONE ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->bExtract |= EXTR_METAL_WAS_NOT_DISCONNECTED;
    // INCHI✔️❌:         }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         prep_inp_data->bDisconnectCoord = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     orig_inp_data->bDisconnectSalts = prep_inp_data->bDisconnectSalts;
    // INCHI✔️❌:     orig_inp_data->bDisconnectCoord = prep_inp_data->bDisconnectCoord;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* 3. if( orig_inp_data->bDisconnectSalts ) then
    // INCHI✔️❌:           disconnect salts in prep_inp_data    */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (( ip->bTautFlags & TG_FLAG_DISCONNECT_SALTS ) && prep_inp_data->bDisconnectSalts &&
    // INCHI✔️❌:          0 < ( i = DisconnectSalts( prep_inp_data, 1 ) )) /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!ip->bNoWarnings)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             WarningMessage( sd->pStrErrStruct, "Salt was disconnected" );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_DISCONNECT_SALTS_DONE;
    // INCHI✔️❌:         if (sd->nErrorType < _IS_WARNING)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->nErrorType = _IS_WARNING;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if ((i = ReconcileAllCmlBondParities( prep_inp_data->at, prep_inp_data->num_inp_atoms, 0 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             char szErrCode[16];
    // INCHI✔️❌:             sprintf(szErrCode, "%d", i);
    // INCHI✔️❌:             AddErrorMessage( sd->pStrErrStruct, "0D Parities Reconciliation failed:" );
    // INCHI✔️❌:             AddErrorMessage( sd->pStrErrStruct, szErrCode );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 )
    // INCHI✔️❌:         sd->bExtract |= EXTR_SALT_WAS_DISCONNECTED;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         prep_inp_data->bDisconnectSalts = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Mark the (disconnected) components in prep_inp_data    */
    // INCHI✔️❌:
    // INCHI✔️❌:     prep_inp_data->num_components = MarkDisconnectedComponents( prep_inp_data, 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (prep_inp_data->num_components < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         AddErrorMessage( sd->pStrErrStruct, "Out of RAM" );
    // INCHI✔️❌:         sd->nStructReadError = 99;
    // INCHI✔️❌:         sd->nErrorType = _IS_FATAL;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Detect isotopic H on heteroatoms -- necessary condition
    // INCHI✔️❌:        for global isotopic tautomerism */
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((i = bNumHeterAtomHasIsotopicH( prep_inp_data->at, prep_inp_data->num_inp_atoms ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i & 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FOUND_ISOTOPIC_H_DONE;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i & 2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* 4a. Detect unusual valences                                              */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (OrigAtData_bCheckUnusualValences( prep_inp_data, 1, sd->pStrErrStruct, ip->bNoWarnings ))
    // INCHI✔️❌:     {
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 )
    // INCHI✔️❌:         sd->bExtract |= EXTR_UNUSUAL_VALENCES;
    // INCHI✔️❌: #else
    // INCHI✔️❌:         ;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*    5. if( orig_inp_data->bDisconnectCoord ) then
    // INCHI✔️❌:               -- copy prep_inp_data --> prep_inp_data+1
    // INCHI✔️❌:               -- disconnect metals in prep_inp_data            */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (prep_inp_data->bDisconnectCoord)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         prep_inp_data->num_components = MarkDisconnectedComponents( prep_inp_data, 0 );
    // INCHI✔️❌:         if (prep_inp_data->num_components < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             AddErrorMessage( sd->pStrErrStruct, "Out of RAM" );
    // INCHI✔️❌:             sd->nStructReadError = 99;
    // INCHI✔️❌:             sd->nErrorType = _IS_FATAL;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Save reconnected structure in prep_inp_data+1 if requested */
    // INCHI✔️❌:         if (0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (0 > OrigAtData_Duplicate( prep_inp_data + 1, prep_inp_data ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 AddErrorMessage( sd->pStrErrStruct, "Out of RAM" );
    // INCHI✔️❌:                 sd->nStructReadError = 99;
    // INCHI✔️❌:                 sd->nErrorType = _IS_FATAL;
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             sd->bTautFlags[INCHI_REC] = sd->bTautFlags[INCHI_BAS];
    // INCHI✔️❌:             sd->bTautFlagsDone[INCHI_REC] = sd->bTautFlagsDone[INCHI_BAS];
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Remove "parity undefined in disconnected structure" flag from reconnected structure */
    // INCHI✔️❌:                 int k, m; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:                 inp_ATOM *at = ( prep_inp_data + 1 )->at;
    // INCHI✔️❌:                 int       num_at = ( prep_inp_data + 1 )->num_inp_atoms;
    // INCHI✔️❌:                 for (k = 0; k < num_at; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (m = 0; m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m]; m++) /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at[k].sb_parity[m] &= SB_PARITY_MASK;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Make disconnected structure in prep_inp_data */
    // INCHI✔️❌:         i = ( 0 != ( ip->bTautFlags & TG_FLAG_CHECK_VALENCE_COORD ) );
    // INCHI✔️❌:
    // INCHI✔️❌:         /*    prep_inp_data->bDisconnectCoord > 1 means add
    // INCHI✔️❌:                 prep_inp_data->bDisconnectCoord-1 explicit H atoms    */
    // INCHI✔️❌:         if (0 < ( i = DisconnectMetals( prep_inp_data, i, &bTautFlagsDone ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!ip->bNoWarnings)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 WarningMessage( sd->pStrErrStruct, "Metal was disconnected" );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_DISCONNECT_COORD_DONE;
    // INCHI✔️❌:             if (sd->nErrorType < _IS_WARNING)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sd->nErrorType = _IS_WARNING;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( bRELEASE_VERSION == 0 )
    // INCHI✔️❌:             sd->bExtract |= EXTR_METAL_WAS_DISCONNECTED;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             /* last parm=1 means find link between unchanged by Metal Disconnection components */
    // INCHI✔️❌:             prep_inp_data->num_components = MarkDisconnectedComponents( prep_inp_data, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:             if (prep_inp_data->num_components < 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 AddErrorMessage( sd->pStrErrStruct, "Out of RAM" );
    // INCHI✔️❌:                 sd->nStructReadError = 99;
    // INCHI✔️❌:                 sd->nErrorType = _IS_FATAL;
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Set parities for the disconnected structure */
    // INCHI✔️❌:                 int k, m, p;
    // INCHI✔️❌:                 inp_ATOM *at = ( prep_inp_data )->at;
    // INCHI✔️❌:                 int       num_at = ( prep_inp_data )->num_inp_atoms;
    // INCHI✔️❌:                 for (k = 0; k < num_at; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (m = 0; m < MAX_NUM_STEREO_BONDS && ( p = at[k].sb_parity[m] ); m++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (p & SB_PARITY_FLAG)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             at[k].sb_parity[m] = ( p >> SB_PARITY_SHFT ) & SB_PARITY_MASK;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if ((i = ReconcileAllCmlBondParities( prep_inp_data->at, prep_inp_data->num_inp_atoms, 1 ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 char szErrCode[16];
    // INCHI✔️❌:                 sprintf(szErrCode, "%d", i);
    // INCHI✔️❌:                 AddErrorMessage( sd->pStrErrStruct, "0D Parities Reconciliation failed:" );
    // INCHI✔️❌:                 AddErrorMessage( sd->pStrErrStruct, szErrCode );
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( REMOVE_ION_PAIRS_DISC_STRU == 1 )
    // INCHI✔️❌:             if (0 < remove_ion_pairs( prep_inp_data->num_inp_atoms, prep_inp_data->at ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!ip->bNoWarnings)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     WarningMessage( sd->pStrErrStruct, "Charges were rearranged" );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (sd->nErrorType < _IS_WARNING)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     sd->nErrorType = _IS_WARNING;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 sd->bTautFlagsDone[INCHI_REC] |= TG_FLAG_FIX_ODD_THINGS_DONE;
    // INCHI✔️❌:                 sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
    // INCHI✔️❌:             }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:             /*
    // INCHI✔️❌:               if prep_inp_data->nOldCompNumber[i] = iINChI+1 > 0 then
    // INCHI✔️❌:               component #(i+1) in prep_inp_data is identical to component #(iINChI+1) in prep_inp_data+1
    // INCHI✔️❌:             */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (i < 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             AddErrorMessage( sd->pStrErrStruct, "Cannot disconnect metal error" );
    // INCHI✔️❌:             sd->nStructReadError = i;
    // INCHI✔️❌:             sd->nErrorType = _IS_ERROR;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Remove "disconnected structure parities" from the structure */
    // INCHI✔️❌:         int k, m; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:         inp_ATOM *at = ( prep_inp_data )->at;
    // INCHI✔️❌:         int       num_at = ( prep_inp_data )->num_inp_atoms;
    // INCHI✔️❌:         for (k = 0; k < num_at; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (m = 0; m < MAX_NUM_STEREO_BONDS && at[k].sb_parity[m]; m++) /* djb-rwth: removing redundant code */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at[k].sb_parity[m] &= SB_PARITY_MASK;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌: if (sd->nErrorType < _IS_ERROR && prep_inp_data)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (0 < post_fix_odd_things( prep_inp_data->num_inp_atoms, prep_inp_data->at ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!ip->bNoWarnings)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 WarningMessage( sd->pStrErrStruct, "Charges were rearranged" );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (sd->nErrorType < _IS_WARNING)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sd->nErrorType = _IS_WARNING;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
    // INCHI✔️❌:             ( prep_inp_data + 1 )->at && ( prep_inp_data + 1 )->num_inp_atoms > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (0 < post_fix_odd_things( ( prep_inp_data + 1 )->num_inp_atoms, ( prep_inp_data + 1 )->at ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!ip->bNoWarnings)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     WarningMessage( sd->pStrErrStruct, "Charges were rearranged" );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (sd->nErrorType < _IS_WARNING)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     sd->nErrorType = _IS_WARNING;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 sd->bTautFlagsDone[INCHI_REC] |= TG_FLAG_FIX_ODD_THINGS_DONE;
    // INCHI✔️❌:                 sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     sd->bTautFlags[INCHI_BAS] |= bTautFlags;  /* TG_FLAG_CHECK_VALENCE_COORD_DONE, TG_FLAG_MOVE_CHARGE_COORD_DONE */
    // INCHI✔️❌:     sd->bTautFlagsDone[INCHI_BAS] |= bTautFlagsDone;  /* TG_FLAG_CHECK_VALENCE_COORD_DONE, TG_FLAG_MOVE_CHARGE_COORD_DONE */
    // INCHI✔️❌:
    // INCHI✔️❌:     return sd->nErrorType;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: PreprocessOneStructure
    // END INCHI C FUNCTION: PreprocessOneStructure
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: PreprocessOneStructure
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; bRELEASE_VERSION=1.
    // INCHI✔️❌: #define FIX_ADJ_RAD 0
    // INCHI✔️❌: #define REMOVE_ION_PAIRS_DISC_STRU 1
    // INCHI✔️❌: #define WarningMessage AddErrorMessage
    // INCHI✔️❌: Debug-only EXTR_* blocks and FixAdjacentRadicals are inactive.
    // INCHI✔️❌: SourceHeap deep-copy/allocation modeling retains known overhead in completed callees.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: PreprocessOneStructure
    let base = INCHI_BAS as usize;
    let reconnected = INCHI_REC as usize;
    let Some(_) = prepared.first() else {
        return Err(SourceHeapError::PointerOutOfBounds);
    };
    let taut_flags: INCHI_MODE = 0;
    let mut taut_flags_done: INCHI_MODE = 0;

    'processing: {
        if OrigAtData_Duplicate(heap, &mut prepared[0], original)? < 0 {
            add_preprocess_message(structure_data, b"Out of RAM")?;
            structure_data.nStructReadError = 99;
            structure_data.nErrorType = _IS_FATAL as i32;
            break 'processing;
        }

        let odd_changes = {
            let atom_count = prepared[0].num_inp_atoms;
            let atoms = if prepared[0].at.is_null() && atom_count <= 0 {
                &mut [][..]
            } else {
                heap.slice_mut(prepared[0].at)?
            };
            fix_odd_things(
                atom_count,
                atoms,
                (input_parameters.bTautFlags & TG_FLAG_FIX_SP3_BUG as INCHI_MODE) as i32,
                input_parameters.bFixNonUniformDraw,
            )?
        };
        if odd_changes > 0 {
            if input_parameters.bNoWarnings == 0 {
                add_preprocess_message(structure_data, b"Charges were rearranged")?;
            }
            if structure_data.nErrorType < _IS_WARNING as i32 {
                structure_data.nErrorType = _IS_WARNING as i32;
            }
            structure_data.bTautFlagsDone[base] |= TG_FLAG_FIX_ODD_THINGS_DONE as INCHI_MODE;
        }

        if input_parameters.bTautFlags & TG_FLAG_DISCONNECT_SALTS as INCHI_MODE != 0 {
            prepared[0].bDisconnectSalts =
                i32::from(DisconnectSalts(heap, &mut prepared[0], 0)? > 0);
        } else {
            prepared[0].bDisconnectSalts = 0;
        }

        if input_parameters.bTautFlags & TG_FLAG_DISCONNECT_COORD as INCHI_MODE != 0 {
            let check_valence = i32::from(
                input_parameters.bTautFlags & TG_FLAG_CHECK_VALENCE_COORD as INCHI_MODE != 0,
            );
            let _ = bMayDisconnectMetals(
                heap,
                &mut prepared[0],
                check_valence,
                Some(&mut taut_flags_done),
            )?;
            structure_data.bTautFlagsDone[base] |= taut_flags_done;
        } else {
            prepared[0].bDisconnectCoord = 0;
        }
        original.bDisconnectSalts = prepared[0].bDisconnectSalts;
        original.bDisconnectCoord = prepared[0].bDisconnectCoord;

        let salt_changes = if input_parameters.bTautFlags & TG_FLAG_DISCONNECT_SALTS as INCHI_MODE
            != 0
            && prepared[0].bDisconnectSalts != 0
        {
            DisconnectSalts(heap, &mut prepared[0], 1)?
        } else {
            0
        };
        if salt_changes > 0 {
            if input_parameters.bNoWarnings == 0 {
                add_preprocess_message(structure_data, b"Salt was disconnected")?;
            }
            structure_data.bTautFlagsDone[base] |= TG_FLAG_DISCONNECT_SALTS_DONE as INCHI_MODE;
            if structure_data.nErrorType < _IS_WARNING as i32 {
                structure_data.nErrorType = _IS_WARNING as i32;
            }
            let parity_status =
                ReconcileAllCmlBondParities(heap, prepared[0].at, prepared[0].num_inp_atoms, 0)?;
            if parity_status != 0 {
                add_preprocess_message(structure_data, b"0D Parities Reconciliation failed:")?;
                add_preprocess_message(structure_data, parity_status.to_string().as_bytes())?;
            }
        } else {
            prepared[0].bDisconnectSalts = 0;
        }

        prepared[0].num_components = MarkDisconnectedComponents(heap, &mut prepared[0], 0)?;
        if prepared[0].num_components < 0 {
            add_preprocess_message(structure_data, b"Out of RAM")?;
            structure_data.nStructReadError = 99;
            structure_data.nErrorType = _IS_FATAL as i32;
            break 'processing;
        }

        let isotope_flags = {
            let atom_count = prepared[0].num_inp_atoms;
            let atoms = if prepared[0].at.is_null() && atom_count <= 0 {
                &[][..]
            } else {
                heap.slice(prepared[0].at.as_const())?
            };
            bNumHeterAtomHasIsotopicH(atoms, atom_count)?
        };
        if isotope_flags & 1 != 0 {
            structure_data.bTautFlagsDone[base] |= TG_FLAG_FOUND_ISOTOPIC_H_DONE as INCHI_MODE;
        }
        if isotope_flags & 2 != 0 {
            structure_data.bTautFlagsDone[base] |= TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE as INCHI_MODE;
        }

        let _ = OrigAtData_bCheckUnusualValences(
            heap,
            Some(&prepared[0]),
            1,
            Some(&mut structure_data.pStrErrStruct),
            input_parameters.bNoWarnings,
        )?;

        if prepared[0].bDisconnectCoord != 0 {
            prepared[0].num_components = MarkDisconnectedComponents(heap, &mut prepared[0], 0)?;
            if prepared[0].num_components < 0 {
                add_preprocess_message(structure_data, b"Out of RAM")?;
                structure_data.nStructReadError = 99;
                structure_data.nErrorType = _IS_FATAL as i32;
                break 'processing;
            }

            if input_parameters.bTautFlags & TG_FLAG_RECONNECT_COORD as INCHI_MODE != 0 {
                let (base_data, saved_data) = prepared.split_at_mut(1);
                let saved = saved_data
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if OrigAtData_Duplicate(heap, saved, &base_data[0])? < 0 {
                    add_preprocess_message(structure_data, b"Out of RAM")?;
                    structure_data.nStructReadError = 99;
                    structure_data.nErrorType = _IS_FATAL as i32;
                    break 'processing;
                }
                structure_data.bTautFlags[reconnected] = structure_data.bTautFlags[base];
                structure_data.bTautFlagsDone[reconnected] = structure_data.bTautFlagsDone[base];
                mask_connected_parities(heap, saved)?;
            }

            let check_valence = i32::from(
                input_parameters.bTautFlags & TG_FLAG_CHECK_VALENCE_COORD as INCHI_MODE != 0,
            );
            let metal_changes = DisconnectMetals(
                heap,
                &mut prepared[0],
                check_valence,
                Some(&mut taut_flags_done),
            )?;
            if metal_changes > 0 {
                if input_parameters.bNoWarnings == 0 {
                    add_preprocess_message(structure_data, b"Metal was disconnected")?;
                }
                structure_data.bTautFlagsDone[base] |= TG_FLAG_DISCONNECT_COORD_DONE as INCHI_MODE;
                if structure_data.nErrorType < _IS_WARNING as i32 {
                    structure_data.nErrorType = _IS_WARNING as i32;
                }

                prepared[0].num_components = MarkDisconnectedComponents(heap, &mut prepared[0], 1)?;
                if prepared[0].num_components < 0 {
                    add_preprocess_message(structure_data, b"Out of RAM")?;
                    structure_data.nStructReadError = 99;
                    structure_data.nErrorType = _IS_FATAL as i32;
                    break 'processing;
                }

                restore_disconnected_parities(heap, &prepared[0])?;
                let parity_status = ReconcileAllCmlBondParities(
                    heap,
                    prepared[0].at,
                    prepared[0].num_inp_atoms,
                    1,
                )?;
                if parity_status != 0 {
                    add_preprocess_message(structure_data, b"0D Parities Reconciliation failed:")?;
                    add_preprocess_message(structure_data, parity_status.to_string().as_bytes())?;
                }

                let rearranged = {
                    let atom_count = prepared[0].num_inp_atoms;
                    let atoms = heap.slice_mut(prepared[0].at)?;
                    remove_ion_pairs(atom_count, atoms)?
                };
                if rearranged > 0 {
                    if input_parameters.bNoWarnings == 0 {
                        add_preprocess_message(structure_data, b"Charges were rearranged")?;
                    }
                    if structure_data.nErrorType < _IS_WARNING as i32 {
                        structure_data.nErrorType = _IS_WARNING as i32;
                    }
                    structure_data.bTautFlagsDone[reconnected] |=
                        TG_FLAG_FIX_ODD_THINGS_DONE as INCHI_MODE;
                    structure_data.bTautFlagsDone[base] |=
                        TG_FLAG_FIX_ODD_THINGS_DONE as INCHI_MODE;
                }
            } else if metal_changes < 0 {
                add_preprocess_message(structure_data, b"Cannot disconnect metal error")?;
                structure_data.nStructReadError = metal_changes;
                structure_data.nErrorType = _IS_ERROR as i32;
                break 'processing;
            }
        } else {
            mask_connected_parities(heap, &prepared[0])?;
        }
    }

    if structure_data.nErrorType < _IS_ERROR as i32 {
        let base_atom_count = prepared[0].num_inp_atoms;
        let base_atoms = if prepared[0].at.is_null() && base_atom_count <= 0 {
            &mut [][..]
        } else {
            heap.slice_mut(prepared[0].at)?
        };
        let _ = post_fix_odd_things(base_atom_count, base_atoms);
        if structure_data.bTautFlagsDone[base] & TG_FLAG_DISCONNECT_COORD_DONE as INCHI_MODE != 0 {
            if let Some(saved) = prepared.get(1) {
                if !saved.at.is_null() && saved.num_inp_atoms > 0 {
                    let saved_atom_count = saved.num_inp_atoms;
                    let saved_atoms = heap.slice_mut(saved.at)?;
                    let _ = post_fix_odd_things(saved_atom_count, saved_atoms);
                }
            }
        }
    }

    structure_data.bTautFlags[base] |= taut_flags;
    structure_data.bTautFlagsDone[base] |= taut_flags_done;
    Ok(structure_data.nErrorType)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_bCheckUnusualValences(
    heap: &SourceHeap,
    original: Option<&ORIG_ATOM_DATA>,
    add_isotopic_hydrogen: i32,
    mut error_text: Option<&mut [i8]>,
    no_warnings: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:89 OrigAtData_bCheckUnusualValences
    // INCHI✔️✔️: direct indexed scan, fixed-size message construction, and completed callees match source complexity.
    /*
    int OrigAtData_bCheckUnusualValences( ORIG_ATOM_DATA *orig_at_data,
                                          int bAddIsoH,
                                          char *pStrErrStruct,
                                          int bNoWarnings)
    {
        int i, val, num_found = 0;
        char msg[32];
        int len, num_H;

        int already_here = ( orig_at_data && orig_at_data->num_inp_atoms > 0 );

        inp_ATOM *at = already_here ? orig_at_data->at : NULL;

        if (at)
        {
            for (i = 0, num_found = 0; i < orig_at_data->num_inp_atoms; i++)
            {
                num_H = bAddIsoH ? NUMH( at, i ) : at[i].num_H;

                val = detect_unusual_el_valence( at[i].el_number,
                                                 at[i].charge,
                                                 at[i].radical,
                                                 at[i].chem_bonds_valence,
                                                 num_H,
                                                 at[i].valence );
                if (val)
                {
                    num_found++;
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErrStruct, "Accepted unusual valence(s):" );
                    }
                    len = sprintf(msg, "%s", at[i].elname);
                    if (at[i].charge)
                    {
                        len += sprintf(msg + len, "%+d", at[i].charge);
                    }
                    if (at[i].radical)
                    {
                        len += sprintf(msg + len, ",%s", at[i].radical == RADICAL_SINGLET ? "s" :
                            at[i].radical == RADICAL_DOUBLET ? "d" :
                            at[i].radical == RADICAL_TRIPLET ? "t" : "?");
                    }
                    len += sprintf(msg + len, "(%d)", val);
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErrStruct, msg );
                    }
                }
            }
        }

        return num_found;
    }
    */
    // END INCHI C FUNCTION: OrigAtData_bCheckUnusualValences
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_bCheckUnusualValences
    // INCHI✔️✔️: #define WarningMessage AddErrorMessage
    // INCHI✔️✔️: #define NUM_ISO_H(AT,N) (AT[N].num_iso_H[0]+AT[N].num_iso_H[1]+AT[N].num_iso_H[2])
    // INCHI✔️✔️: #define NUMH(AT,N) (AT[N].num_H + NUM_ISO_H(AT,N))
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_bCheckUnusualValences

    let Some(original) = original.filter(|value| value.num_inp_atoms > 0) else {
        return Ok(0);
    };
    if original.at.is_null() {
        return Ok(0);
    }
    let atoms = heap
        .slice(original.at.as_const())?
        .get(
            ..usize::try_from(original.num_inp_atoms)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
        )
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut found = 0_i32;
    for atom in atoms {
        let hydrogen_count = if add_isotopic_hydrogen != 0 {
            atom.num_iso_H
                .iter()
                .fold(i32::from(atom.num_H), |sum, count| {
                    sum.wrapping_add(i32::from(*count))
                })
        } else {
            i32::from(atom.num_H)
        };
        let valence = detect_unusual_el_valence(
            i32::from(atom.el_number),
            i32::from(atom.charge),
            i32::from(atom.radical),
            i32::from(atom.chem_bonds_valence),
            hydrogen_count,
            i32::from(atom.valence),
        )?;
        if valence == 0 {
            continue;
        }
        found = found.wrapping_add(1);
        if no_warnings != 0 {
            continue;
        }
        let heading = (*b"Accepted unusual valence(s):\0").map(|byte| byte as i8);
        AddErrorMessage(error_text.as_deref_mut(), Some(&heading))?;

        let symbol_length = atom
            .elname
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let mut message = atom.elname[..symbol_length].to_vec();
        if atom.charge != 0 {
            let charge = i32::from(atom.charge);
            if charge > 0 {
                message.push(b'+' as i8);
            }
            message.extend(charge.to_string().bytes().map(|byte| byte as i8));
        }
        if atom.radical != 0 {
            message.push(b',' as i8);
            message.push(match u32::try_from(atom.radical).ok() {
                Some(RADICAL_SINGLET) => b's' as i8,
                Some(RADICAL_DOUBLET) => b'd' as i8,
                Some(RADICAL_TRIPLET) => b't' as i8,
                _ => b'?' as i8,
            });
        }
        message.push(b'(' as i8);
        message.extend(valence.to_string().bytes().map(|byte| byte as i8));
        message.extend([b')' as i8, 0]);
        AddErrorMessage(error_text.as_deref_mut(), Some(&message))?;
    }
    Ok(found)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_Duplicate(
    heap: &mut SourceHeap,
    destination: &mut ORIG_ATOM_DATA,
    source: &ORIG_ATOM_DATA,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:148 OrigAtData_Duplicate
    // INCHI✔️❌: complete verbatim source frame is inserted below; SourceHeap deep copies add model overhead.
    /*
    int OrigAtData_Duplicate( ORIG_ATOM_DATA *new_orig_atom,
                              ORIG_ATOM_DATA *orig_atom )
    {
        inp_ATOM  *at = NULL;
        AT_NUMB   *nCurAtLen = NULL;
        AT_NUMB   *nOldCompNumber = NULL;
        int k, m, nn;
        int orig_nat = orig_atom->num_inp_atoms;

        int ret = -1; /* fail; 0 - OK */

        if (new_orig_atom->at &&
             new_orig_atom->num_inp_atoms >= orig_nat)
        {
            at = new_orig_atom->at;
        }
        else
        {
            at = (inp_ATOM *) inchi_calloc( (long long)orig_nat + 1, sizeof( at[0] ) ); /* djb-rwth: cast operator added */
            if (!at)
            {
                goto exit_function;
            }
        }

        if (new_orig_atom->nOldCompNumber &&
             new_orig_atom->num_components >= orig_atom->num_components)
        {
            nCurAtLen = new_orig_atom->nCurAtLen;
        }
        else
        {
            nCurAtLen = (AT_NUMB *) inchi_calloc( (long long)orig_atom->num_components + 1, sizeof( nCurAtLen[0] ) ); /* djb-rwth: cast operator added */
            if (!nCurAtLen)
            {
                goto exit_function;
            }
        }

        if (new_orig_atom->nCurAtLen &&
             new_orig_atom->num_components >= orig_atom->num_components)
        {
            nOldCompNumber = new_orig_atom->nOldCompNumber;
        }
        else
        {
            nOldCompNumber = (AT_NUMB *) inchi_calloc( (long long)orig_atom->num_components + 1,
                                                      sizeof( nOldCompNumber[0] ) ); /* djb-rwth: cast operator added */
            if (!nOldCompNumber)
            {
                goto exit_function;
            }
        }

        if (at && nCurAtLen && nOldCompNumber)
        {
            /* Copy */
            if (orig_atom->at)
            {
                memcpy(at, orig_atom->at,
                    orig_nat * sizeof(new_orig_atom->at[0]));
            }
            if (orig_atom->nCurAtLen)
            {
                memcpy(nCurAtLen, orig_atom->nCurAtLen,
                    orig_atom->num_components * sizeof(nCurAtLen[0]));
            }
            if (orig_atom->nOldCompNumber)
            {
                memcpy(nOldCompNumber, orig_atom->nOldCompNumber,
                    orig_atom->num_components * sizeof(nOldCompNumber[0]));
            }

            /* Deallocate */
            if (new_orig_atom->at && new_orig_atom->at != at)
            {
                inchi_free( new_orig_atom->at );
            }
            if (new_orig_atom->nCurAtLen && new_orig_atom->nCurAtLen != nCurAtLen)
            {
                inchi_free( new_orig_atom->nCurAtLen );
            }
            if (new_orig_atom->nOldCompNumber &&
                 new_orig_atom->nOldCompNumber != nOldCompNumber)
            {
                inchi_free( new_orig_atom->nOldCompNumber );
            }

            *new_orig_atom = *orig_atom;
            new_orig_atom->at = at;
            new_orig_atom->nCurAtLen = nCurAtLen;
            new_orig_atom->nOldCompNumber = nOldCompNumber;

            /* Data that are not to be copied */
            new_orig_atom->nNumEquSets = 0;
            memset( new_orig_atom->bSavedInINCHI_LIB, 0, sizeof( new_orig_atom->bSavedInINCHI_LIB ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( new_orig_atom->bPreprocessed, 0, sizeof( new_orig_atom->bPreprocessed ) ); /* djb-rwth: memset_s C11/Annex K variant? */



            new_orig_atom->szCoord = NULL;
            if (orig_atom->szCoord)
            {
                new_orig_atom->szCoord = (MOL_COORD *) inchi_calloc(orig_nat, sizeof(new_orig_atom->szCoord[0]));
                if (!new_orig_atom->szCoord)
                {
                    goto exit_function;
                }
                memcpy(new_orig_atom->szCoord, orig_atom->szCoord, orig_nat * sizeof(new_orig_atom->szCoord[0]));
            }


            /* Arrays that are not to be copied */

            new_orig_atom->nEquLabels = NULL;
            new_orig_atom->nSortedOrder = NULL;

            new_orig_atom->polymer = NULL;
            if (orig_atom->polymer)
            {
                /* Polymer stuff -- deep copy */
                OAD_Polymer *oldp = orig_atom->polymer;
                OAD_Polymer *newp = NULL;

                newp = (OAD_Polymer *) inchi_calloc( 1, sizeof( OAD_Polymer ) );
                if (!newp)
                {
                    inchi_free(newp); /* djb-rwth: avoiding memory leak */
                    goto exit_function;
                }
                memcpy(newp, orig_atom->polymer, sizeof(OAD_Polymer));
                newp->units = (OAD_PolymerUnit**) inchi_calloc( newp->n, sizeof(OAD_PolymerUnit*) ); /* djb-rwth: inchi_calloc must return OAD_PolymerUnit** */
                if (!newp->units)
                {
                    inchi_free(newp); /* djb-rwth: avoiding memory leak */
                    goto exit_function;
                }
                for (k = 0; k < orig_atom->polymer->n; k++)
                {
                    newp->units[k] = OAD_PolymerUnit_CreateCopy( orig_atom->polymer->units[k] );
                }
                if (oldp->n_pzz > 0)
                {
                    newp->n_pzz = oldp->n_pzz;
                    newp->pzz = (int *) inchi_calloc( newp->n_pzz, sizeof( int ) );
                    if (!newp->pzz)
                    {
                        inchi_free(newp->units); /* djb-rwth: fixing coverity ID #499546 */
                        inchi_free(newp); /* djb-rwth: avoiding memory leak */
                        goto exit_function;
                    }
                    memcpy(newp->pzz, oldp->pzz, newp->n_pzz * sizeof(oldp->pzz[0]));
                }
                new_orig_atom->polymer = newp;
            }

            new_orig_atom->v3000 = NULL;
            if (orig_atom->v3000)
            {
                /* V3000 features -- deep copy */
                OAD_V3000 *new_v3000 = NULL;
                new_v3000 = (OAD_V3000 *) inchi_calloc( 1, sizeof( OAD_V3000 ) );
                if (!new_v3000)
                {
                    inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                    goto exit_function;
                }
                memcpy(new_v3000, orig_atom->v3000, sizeof(OAD_V3000));
                if (orig_atom->v3000->atom_index_orig)
                {
                    new_v3000->atom_index_orig = (int *) inchi_calloc( orig_nat, sizeof( int ) );
                    /* if ( NULL==new_v3000->atom_index_orig ) {TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; } */
                    if (!new_v3000->atom_index_orig)
                    {
                        inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                        goto exit_function;
                    }
                    memcpy(new_v3000->atom_index_orig, orig_atom->v3000->atom_index_orig, orig_nat * sizeof(int));
                }
                if (orig_atom->v3000->atom_index_fin)
                {
                    new_v3000->atom_index_fin = (int *) inchi_calloc( orig_nat, sizeof( int ) );
                    /* if ( NULL==new_v3000->atom_index_fin ) {TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; } */
                    if (!new_v3000->atom_index_fin)
                    {
                        inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                        goto exit_function;
                    }
                    memcpy(new_v3000->atom_index_fin, orig_atom->v3000->atom_index_fin, orig_nat * sizeof(int));
                }
                if (orig_atom->v3000->n_haptic_bonds && orig_atom->v3000->lists_haptic_bonds)
                {
                    new_v3000->lists_haptic_bonds = (int **) inchi_calloc( orig_atom->v3000->n_haptic_bonds, sizeof( int* ) );
                    /* if ( NULL==new_v3000->lists_haptic_bonds ) { TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; }*/
                    for (m = 0; m < orig_atom->v3000->n_haptic_bonds; m++)
                    {
                        int *lst = NULL;
                        int *old_lst = orig_atom->v3000->lists_haptic_bonds[m];
                        nn = old_lst[2] + 3;
                        lst = new_v3000->lists_haptic_bonds[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                        if (!lst)
                        {
                            inchi_free(new_v3000->lists_haptic_bonds); /* djb-rwth: fixing coverity ID #499504 */
                            inchi_free(new_v3000->atom_index_orig); /* djb-rwth: fixing coverity ID #499540 */
                            inchi_free(new_v3000->atom_index_fin); /* djb-rwth: fixing coverity ID #499613 */
                            inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                            goto exit_function;
                        }
                        memcpy(lst, old_lst, nn * sizeof(int));
                    }
                }
                if (orig_atom->v3000->n_steabs && orig_atom->v3000->lists_steabs)
                {
                    new_v3000->lists_steabs = (int **) inchi_calloc( orig_atom->v3000->n_steabs, sizeof( int* ) );
                    /* if ( NULL==new_v3000->lists_steabs ) { TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; }*/
                    for (m = 0; m < orig_atom->v3000->n_steabs; m++)
                    {
                        int *lst = NULL;
                        int *old_lst = orig_atom->v3000->lists_steabs[m];
                        nn = old_lst[1] + 2;
                        lst = new_v3000->lists_steabs[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                        if (!lst)
                        {
                            inchi_free(new_v3000->lists_haptic_bonds); /* djb-rwth: fixing coverity ID #499504 */
                            inchi_free(new_v3000->lists_steabs); /* djb-rwth: fixing coverity ID #499543 */
                            inchi_free(new_v3000->atom_index_orig); /* djb-rwth: fixing coverity ID #499540 */
                            inchi_free(new_v3000->atom_index_fin); /* djb-rwth: fixing coverity ID #499613 */
                            inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                            goto exit_function;
                        }
                        memcpy(lst, old_lst, nn * sizeof(int));
                    }
                }
                if (orig_atom->v3000->n_sterel && orig_atom->v3000->lists_sterel)
                {
                    new_v3000->lists_sterel = (int **) inchi_calloc( orig_atom->v3000->n_sterel, sizeof( int* ) );
                    if (!new_v3000)
                    {
                        inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                        goto exit_function;
                    }
                    /* if ( NULL==new_v3000->lists_sterel ) { TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; }*/
                    for (m = 0; m < orig_atom->v3000->n_sterel; m++)
                    {
                        int *lst = NULL;
                        int *old_lst = orig_atom->v3000->lists_sterel[m];
                        nn = old_lst[1] + 2;
                        if (new_v3000->lists_sterel) /* djb-rwth: fixing a NULL pointer dereference */
                            lst = new_v3000->lists_sterel[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                        if (!lst)
                        {
                            inchi_free(new_v3000->lists_haptic_bonds); /* djb-rwth: fixing coverity ID #499504 */
                            inchi_free(new_v3000->lists_steabs); /* djb-rwth: fixing coverity ID #499543 */
                            inchi_free(new_v3000->lists_sterel); /* djb-rwth: fixing coverity ID #499504 */
                            inchi_free(new_v3000->atom_index_orig); /* djb-rwth: fixing coverity ID #499540 */
                            inchi_free(new_v3000->atom_index_fin); /* djb-rwth: fixing coverity ID #499613 */
                            inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                            goto exit_function;
                        }
                        memcpy(lst, old_lst, nn * sizeof(int));
                    }
                }
                if (orig_atom->v3000->n_sterac && orig_atom->v3000->lists_sterac)
                {
                    new_v3000->lists_sterac = (int **) inchi_calloc( orig_atom->v3000->n_sterac, sizeof( int* ) );
                    /* if ( NULL==new_v3000->lists_sterac ) { TREAT_ERR( err, 9001, "Out of RAM"); goto exit_function; }*/
                    if (new_v3000->lists_sterac) /* djb-rwth: fixing a NULL pointer dereference */
                    {
                        for (m = 0; m < orig_atom->v3000->n_sterac; m++)
                        {
                            int* lst = NULL;
                            int* old_lst = orig_atom->v3000->lists_sterac[m];
                            nn = old_lst[1] + 2;
                            lst = new_v3000->lists_sterac[m] = (int*)inchi_calloc(nn, sizeof(int));
                            if (!lst)
                            {
                                inchi_free(new_v3000->lists_haptic_bonds); /* djb-rwth: fixing coverity ID #499504 */
                                inchi_free(new_v3000->lists_steabs); /* djb-rwth: fixing coverity ID #499543 */
                                inchi_free(new_v3000->lists_sterel); /* djb-rwth: fixing coverity ID #499504 */
                                inchi_free(new_v3000->lists_sterac); /* djb-rwth: fixing coverity ID #499575 */
                                inchi_free(new_v3000->atom_index_orig); /* djb-rwth: fixing coverity ID #499540 */
                                inchi_free(new_v3000->atom_index_fin); /* djb-rwth: fixing coverity ID #499613 */
                                inchi_free(new_v3000); /* djb-rwth: avoiding memory leak */
                                goto exit_function;
                            }
                            memcpy(lst, old_lst, nn * sizeof(int));
                        }
                    }
                }

                new_orig_atom->v3000 = new_v3000;
            }

            /* Success */
            ret = 0;
        }

    exit_function:

        if (ret != 0)
        {
            if (at && new_orig_atom->at != at)
                inchi_free( at );
            if (nCurAtLen && new_orig_atom->nCurAtLen != nCurAtLen)
                inchi_free( nCurAtLen );
            if (nOldCompNumber && new_orig_atom->nOldCompNumber != nOldCompNumber)
                inchi_free( nOldCompNumber );
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: OrigAtData_Duplicate
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_Duplicate
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; inchi_calloc=calloc; inchi_free=free.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_Duplicate

    let original_atom_count = source.num_inp_atoms;
    let original_component_count = source.num_components;
    let mut atoms = SourceMutPointer::<inp_ATOM>::null();
    let mut current_lengths = SourceMutPointer::<AT_NUMB>::null();
    let mut old_component_numbers = SourceMutPointer::<AT_NUMB>::null();
    let mut result = -1_i32;

    'duplicate: {
        if !destination.at.is_null() && destination.num_inp_atoms >= original_atom_count {
            atoms = destination.at;
        } else {
            let count = match i64::from(original_atom_count)
                .checked_add(1)
                .and_then(|value| u64::try_from(value).ok())
            {
                Some(count) => count,
                None => break 'duplicate,
            };
            atoms = match inchi_calloc(heap, count, SOURCE_SIZEOF_INP_ATOM) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => break 'duplicate,
                Err(error) => return Err(error),
            };
        }

        if !destination.nOldCompNumber.is_null()
            && destination.num_components >= original_component_count
        {
            current_lengths = destination.nCurAtLen;
        } else {
            let count = match i64::from(original_component_count)
                .checked_add(1)
                .and_then(|value| u64::try_from(value).ok())
            {
                Some(count) => count,
                None => break 'duplicate,
            };
            current_lengths = match inchi_calloc(heap, count, 2) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => break 'duplicate,
                Err(error) => return Err(error),
            };
        }

        if !destination.nCurAtLen.is_null()
            && destination.num_components >= original_component_count
        {
            old_component_numbers = destination.nOldCompNumber;
        } else {
            let count = match i64::from(original_component_count)
                .checked_add(1)
                .and_then(|value| u64::try_from(value).ok())
            {
                Some(count) => count,
                None => break 'duplicate,
            };
            old_component_numbers = match inchi_calloc(heap, count, 2) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => break 'duplicate,
                Err(error) => return Err(error),
            };
        }

        if atoms.is_null() || current_lengths.is_null() || old_component_numbers.is_null() {
            break 'duplicate;
        }

        let atom_count = usize::try_from(original_atom_count)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let component_count = usize::try_from(original_component_count)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if !source.at.is_null() {
            let values = heap
                .slice(source.at.as_const())?
                .get(..atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            heap.slice_mut(atoms)?
                .get_mut(..atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone_from_slice(&values);
        }
        if !source.nCurAtLen.is_null() {
            let values = heap
                .slice(source.nCurAtLen.as_const())?
                .get(..component_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            heap.slice_mut(current_lengths)?
                .get_mut(..component_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .copy_from_slice(&values);
        }
        if !source.nOldCompNumber.is_null() {
            let values = heap
                .slice(source.nOldCompNumber.as_const())?
                .get(..component_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            heap.slice_mut(old_component_numbers)?
                .get_mut(..component_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .copy_from_slice(&values);
        }

        if !destination.at.is_null() && destination.at != atoms {
            inchi_free(heap, destination.at)?;
        }
        if !destination.nCurAtLen.is_null() && destination.nCurAtLen != current_lengths {
            inchi_free(heap, destination.nCurAtLen)?;
        }
        if !destination.nOldCompNumber.is_null()
            && destination.nOldCompNumber != old_component_numbers
        {
            inchi_free(heap, destination.nOldCompNumber)?;
        }

        *destination = source.clone();
        destination.at = atoms;
        destination.nCurAtLen = current_lengths;
        destination.nOldCompNumber = old_component_numbers;
        destination.nNumEquSets = 0;
        destination.bSavedInINCHI_LIB = [0; 2];
        destination.bPreprocessed = [0; 2];

        destination.szCoord = SourceMutPointer::null();
        if !source.szCoord.is_null() {
            let count = match u64::try_from(original_atom_count) {
                Ok(count) => count,
                Err(_) => break 'duplicate,
            };
            let coordinates = match inchi_calloc::<MOL_COORD>(heap, count, 32) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => break 'duplicate,
                Err(error) => return Err(error),
            };
            let values = heap
                .slice(source.szCoord.as_const())?
                .get(..atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec();
            heap.slice_mut(coordinates)?
                .get_mut(..atom_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .copy_from_slice(&values);
            destination.szCoord = coordinates;
        }

        destination.nEquLabels = SourceMutPointer::null();
        destination.nSortedOrder = SourceMutPointer::null();

        destination.polymer = SourceMutPointer::null();
        if !source.polymer.is_null() {
            let old_polymer = heap
                .slice(source.polymer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let new_polymer = match inchi_calloc::<OAD_Polymer>(heap, 1, SOURCE_SIZEOF_OAD_POLYMER)
            {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => break 'duplicate,
                Err(error) => return Err(error),
            };
            heap.slice_mut(new_polymer)?[0] = old_polymer.clone();
            let unit_count = match u64::try_from(old_polymer.n) {
                Ok(count) => count,
                Err(_) => {
                    inchi_free(heap, new_polymer)?;
                    break 'duplicate;
                }
            };
            let units = match inchi_calloc::<SourceMutPointer<OAD_PolymerUnit>>(heap, unit_count, 8)
            {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => {
                    inchi_free(heap, new_polymer)?;
                    break 'duplicate;
                }
                Err(error) => return Err(error),
            };
            heap.slice_mut(new_polymer)?[0].units = units;
            for index in 0..old_polymer.n {
                let old_unit = *heap
                    .slice(old_polymer.units.as_const())?
                    .get(index as usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let new_unit = OAD_PolymerUnit_CreateCopy(heap, old_unit)?;
                heap.slice_mut(units)?[index as usize] = new_unit;
            }
            if old_polymer.n_pzz > 0 {
                let pzz_count = u64::try_from(old_polymer.n_pzz)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let pzz = match inchi_calloc::<i32>(heap, pzz_count, 4) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => {
                        inchi_free(heap, units)?;
                        inchi_free(heap, new_polymer)?;
                        break 'duplicate;
                    }
                    Err(error) => return Err(error),
                };
                let values = heap
                    .slice(old_polymer.pzz.as_const())?
                    .get(..old_polymer.n_pzz as usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                heap.slice_mut(pzz)?.copy_from_slice(&values);
                let copied_polymer = &mut heap.slice_mut(new_polymer)?[0];
                copied_polymer.n_pzz = old_polymer.n_pzz;
                copied_polymer.pzz = pzz;
            }
            destination.polymer = new_polymer;
        }

        destination.v3000 = SourceMutPointer::null();
        if !source.v3000.is_null() {
            let old_v3000 = heap
                .slice(source.v3000.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let new_v3000 = match inchi_calloc::<OAD_V3000>(heap, 1, SOURCE_SIZEOF_OAD_V3000) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => break 'duplicate,
                Err(error) => return Err(error),
            };
            heap.slice_mut(new_v3000)?[0] = old_v3000.clone();

            if !old_v3000.atom_index_orig.is_null() {
                let indices = match inchi_calloc::<i32>(heap, atom_count as u64, 4) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => {
                        inchi_free(heap, new_v3000)?;
                        break 'duplicate;
                    }
                    Err(error) => return Err(error),
                };
                let values = heap
                    .slice(old_v3000.atom_index_orig.as_const())?
                    .get(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                heap.slice_mut(indices)?.copy_from_slice(&values);
                heap.slice_mut(new_v3000)?[0].atom_index_orig = indices;
            }
            if !old_v3000.atom_index_fin.is_null() {
                let indices = match inchi_calloc::<i32>(heap, atom_count as u64, 4) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => {
                        inchi_free(heap, new_v3000)?;
                        break 'duplicate;
                    }
                    Err(error) => return Err(error),
                };
                let values = heap
                    .slice(old_v3000.atom_index_fin.as_const())?
                    .get(..atom_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                heap.slice_mut(indices)?.copy_from_slice(&values);
                heap.slice_mut(new_v3000)?[0].atom_index_fin = indices;
            }

            if old_v3000.n_haptic_bonds != 0 && !old_v3000.lists_haptic_bonds.is_null() {
                let count = u64::try_from(old_v3000.n_haptic_bonds)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let lists = inchi_calloc::<SourceMutPointer<i32>>(heap, count, 8)?;
                heap.slice_mut(new_v3000)?[0].lists_haptic_bonds = lists;
                for index in 0..old_v3000.n_haptic_bonds as usize {
                    let old_list = heap.slice(old_v3000.lists_haptic_bonds.as_const())?[index];
                    let list_count = i64::from(heap.slice(old_list.as_const())?[2])
                        .checked_add(3)
                        .and_then(|value| u64::try_from(value).ok())
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    let list = match inchi_calloc::<i32>(heap, list_count, 4) {
                        Ok(pointer) => pointer,
                        Err(SourceHeapError::AllocationFailed) => {
                            inchi_free(heap, lists)?;
                            let v = heap.slice(new_v3000.as_const())?[0].clone();
                            inchi_free(heap, v.atom_index_orig)?;
                            inchi_free(heap, v.atom_index_fin)?;
                            inchi_free(heap, new_v3000)?;
                            break 'duplicate;
                        }
                        Err(error) => return Err(error),
                    };
                    let values = heap.slice(old_list.as_const())?[..list_count as usize].to_vec();
                    heap.slice_mut(list)?.copy_from_slice(&values);
                    heap.slice_mut(lists)?[index] = list;
                }
            }

            if old_v3000.n_steabs != 0 && !old_v3000.lists_steabs.is_null() {
                let count = u64::try_from(old_v3000.n_steabs)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let lists = inchi_calloc::<SourceMutPointer<i32>>(heap, count, 8)?;
                heap.slice_mut(new_v3000)?[0].lists_steabs = lists;
                for index in 0..old_v3000.n_steabs as usize {
                    let old_list = heap.slice(old_v3000.lists_steabs.as_const())?[index];
                    let list_count = i64::from(heap.slice(old_list.as_const())?[1])
                        .checked_add(2)
                        .and_then(|value| u64::try_from(value).ok())
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    let list = match inchi_calloc::<i32>(heap, list_count, 4) {
                        Ok(pointer) => pointer,
                        Err(SourceHeapError::AllocationFailed) => {
                            let v = heap.slice(new_v3000.as_const())?[0].clone();
                            inchi_free(heap, v.lists_haptic_bonds)?;
                            inchi_free(heap, lists)?;
                            inchi_free(heap, v.atom_index_orig)?;
                            inchi_free(heap, v.atom_index_fin)?;
                            inchi_free(heap, new_v3000)?;
                            break 'duplicate;
                        }
                        Err(error) => return Err(error),
                    };
                    let values = heap.slice(old_list.as_const())?[..list_count as usize].to_vec();
                    heap.slice_mut(list)?.copy_from_slice(&values);
                    heap.slice_mut(lists)?[index] = list;
                }
            }

            if old_v3000.n_sterel != 0 && !old_v3000.lists_sterel.is_null() {
                let count = u64::try_from(old_v3000.n_sterel)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let lists = match inchi_calloc::<SourceMutPointer<i32>>(heap, count, 8) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                    Err(error) => return Err(error),
                };
                heap.slice_mut(new_v3000)?[0].lists_sterel = lists;
                for index in 0..old_v3000.n_sterel as usize {
                    let old_list = heap.slice(old_v3000.lists_sterel.as_const())?[index];
                    let list_count = i64::from(heap.slice(old_list.as_const())?[1])
                        .checked_add(2)
                        .and_then(|value| u64::try_from(value).ok())
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    let list = if lists.is_null() {
                        SourceMutPointer::null()
                    } else {
                        match inchi_calloc::<i32>(heap, list_count, 4) {
                            Ok(pointer) => pointer,
                            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                            Err(error) => return Err(error),
                        }
                    };
                    if list.is_null() {
                        let v = heap.slice(new_v3000.as_const())?[0].clone();
                        inchi_free(heap, v.lists_haptic_bonds)?;
                        inchi_free(heap, v.lists_steabs)?;
                        inchi_free(heap, lists)?;
                        inchi_free(heap, v.atom_index_orig)?;
                        inchi_free(heap, v.atom_index_fin)?;
                        inchi_free(heap, new_v3000)?;
                        break 'duplicate;
                    }
                    let values = heap.slice(old_list.as_const())?[..list_count as usize].to_vec();
                    heap.slice_mut(list)?.copy_from_slice(&values);
                    heap.slice_mut(lists)?[index] = list;
                }
            }

            if old_v3000.n_sterac != 0 && !old_v3000.lists_sterac.is_null() {
                let count = u64::try_from(old_v3000.n_sterac)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let lists = match inchi_calloc::<SourceMutPointer<i32>>(heap, count, 8) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                    Err(error) => return Err(error),
                };
                heap.slice_mut(new_v3000)?[0].lists_sterac = lists;
                if !lists.is_null() {
                    for index in 0..old_v3000.n_sterac as usize {
                        let old_list = heap.slice(old_v3000.lists_sterac.as_const())?[index];
                        let list_count = i64::from(heap.slice(old_list.as_const())?[1])
                            .checked_add(2)
                            .and_then(|value| u64::try_from(value).ok())
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                        let list = match inchi_calloc::<i32>(heap, list_count, 4) {
                            Ok(pointer) => pointer,
                            Err(SourceHeapError::AllocationFailed) => {
                                let v = heap.slice(new_v3000.as_const())?[0].clone();
                                inchi_free(heap, v.lists_haptic_bonds)?;
                                inchi_free(heap, v.lists_steabs)?;
                                inchi_free(heap, v.lists_sterel)?;
                                inchi_free(heap, lists)?;
                                inchi_free(heap, v.atom_index_orig)?;
                                inchi_free(heap, v.atom_index_fin)?;
                                inchi_free(heap, new_v3000)?;
                                break 'duplicate;
                            }
                            Err(error) => return Err(error),
                        };
                        let values =
                            heap.slice(old_list.as_const())?[..list_count as usize].to_vec();
                        heap.slice_mut(list)?.copy_from_slice(&values);
                        heap.slice_mut(lists)?[index] = list;
                    }
                }
            }
            destination.v3000 = new_v3000;
        }

        result = 0;
    }

    if result != 0 {
        if !atoms.is_null() && destination.at != atoms {
            inchi_free(heap, atoms)?;
        }
        if !current_lengths.is_null() && destination.nCurAtLen != current_lengths {
            inchi_free(heap, current_lengths)?;
        }
        if !old_component_numbers.is_null() && destination.nOldCompNumber != old_component_numbers {
            inchi_free(heap, old_component_numbers)?;
        }
    }
    Ok(result)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OAD_PolymerUnit_New(
    heap: &mut SourceHeap,
    maxatoms: i32,
    maxbonds: i32,
    id: i32,
    label: i32,
    type_: i32,
    subtype: i32,
    conn: i32,
    smt: SourceConstPointer<i8>,
    na: i32,
    alist: Option<&crate::source_types::INT_ARRAY>,
    nb: i32,
    blist: Option<&crate::source_types::INT_ARRAY>,
    nbkbonds: i32,
    _bkbonds: SourceMutPointer<SourceMutPointer<i32>>,
) -> Result<SourceMutPointer<OAD_PolymerUnit>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1162 OAD_PolymerUnit_New
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
OAD_PolymerUnit* OAD_PolymerUnit_New( int       maxatoms,
                                      int       maxbonds,
                                      int       id,
                                      int       label,
                                      int       type,
                                      int       subtype,
                                      int       conn,
                                      char      *smt,
                                      int       na,
                                      INT_ARRAY *alist,
                                      int       nb,
                                      INT_ARRAY *blist,
                                      int       nbkbonds,
                                      int       **bkbonds )
{
    int k, err = 0;
    OAD_PolymerUnit *u2 = NULL;

    u2 = (OAD_PolymerUnit*) inchi_calloc( 1, sizeof( OAD_PolymerUnit ) );
    if (NULL == u2)
    {
        err = 1;
        goto exit_function;
    }

    u2->id = id;
    u2->label = label;
    u2->type = type;
    u2->subtype = subtype;
    u2->conn = conn;
    u2->na = na;
    u2->nb = nb;
    u2->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
    u2->cyclized = 0;
    for (k = 0; k < 4; k++)
    {
        u2->xbr1[k] = 0.0;
        u2->xbr2[k] = 0.0;
    }
    strcpy(u2->smt, smt);
    u2->cap1 = -1;
    u2->end_atom1 = -1;
    u2->cap2 = -1;
    u2->end_atom2 = -1;
    u2->maxbkbonds = maxbonds;
    u2->nbkbonds = nbkbonds;
    u2->cap1_is_undef = 0;
    u2->cap2_is_undef = 0;

    u2->alist = NULL;
    if (na > 0 || maxatoms > 0)
    {
        u2->alist	= (int *) inchi_calloc( na > 0 ? na : maxatoms, sizeof( int ) );
        if (!u2->alist )
        {
            err = 2;
            goto exit_function;
        }
        for (k = 0; k < na; k++)
        {
            u2->alist[k]	= alist->item[k];
        }
    }
    u2->blist = NULL;
    if (nb > 0 || maxbonds > 0)
    {
        u2->blist = (int *) inchi_calloc( nb > 0 ? 2 * nb : 2 * maxbonds, sizeof( int ) );
        if (!u2->blist)
        {
            err = 3;
            goto exit_function;
        }
        if (blist)
        {
            for (k = 0; k < 2 * nb; k++)
            {
                u2->blist[k]	= blist->item[k];
            }
        }
        
    }
    u2->bkbonds = NULL;
    
exit_function:

    if (err)
    {
        OAD_PolymerUnit_Free( u2 );
        return NULL;
    }

    return u2;
}
    */
    // END INCHI C FUNCTION: OAD_PolymerUnit_New
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_PolymerUnit_New
    // INCHI✔️❌: TARGET_API_LIB and COMPILE_ANSI_ONLY include this polymer constructor unchanged.
    // INCHI✔️❌: inchi_calloc maps to GCC/Linux calloc and CLOSING_SRU_NOT_APPLICABLE is the active zero-valued enum constant.
    // INCHI✔️❌: The unit/alist/blist allocation order, source-list copy bounds, ignored bkbonds input, and failure cleanup are preserved.
    // INCHI✔️❌: Rust checked heap lookups and temporary source slices add work absent from direct C pointer access.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_PolymerUnit_New

    fn calloc_or_null<T: Clone + Default + 'static>(
        heap: &mut SourceHeap,
        count: u64,
    ) -> Result<SourceMutPointer<T>, SourceHeapError> {
        match inchi_calloc::<T>(heap, count, std::mem::size_of::<T>() as u64) {
            Ok(pointer) => Ok(pointer),
            Err(
                SourceHeapError::AllocationFailed
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationSizeOverflow,
            ) => Ok(SourceMutPointer::null()),
            Err(error) => Err(error),
        }
    }

    let unit = calloc_or_null::<OAD_PolymerUnit>(heap, 1)?;
    if unit.is_null() {
        return Ok(SourceMutPointer::null());
    }
    let smt_bytes = heap.slice(smt)?;
    let smt_len = smt_bytes.iter().position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    if smt_len >= 80 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let smt_bytes = smt_bytes[..=smt_len].to_vec();
    {
        let value = heap.slice_mut(unit)?.first_mut().ok_or(SourceHeapError::PointerOutOfBounds)?;
        value.id = id;
        value.label = label;
        value.type_ = type_;
        value.subtype = subtype;
        value.conn = conn;
        value.na = na;
        value.nb = nb;
        value.cyclizable = CLOSING_SRU_NOT_APPLICABLE as i32;
        value.cyclized = 0;
        value.xbr1 = [0.0; 4];
        value.xbr2 = [0.0; 4];
        value.smt[..=smt_len].copy_from_slice(&smt_bytes);
        value.cap1 = -1;
        value.end_atom1 = -1;
        value.cap2 = -1;
        value.end_atom2 = -1;
        value.maxbkbonds = maxbonds;
        value.nbkbonds = nbkbonds;
        value.cap1_is_undef = 0;
        value.cap2_is_undef = 0;
        value.alist = SourceMutPointer::null();
        value.blist = SourceMutPointer::null();
        value.bkbonds = SourceMutPointer::null();
    }

    if na > 0 || maxatoms > 0 {
        let count = if na > 0 { na } else { maxatoms };
        let pointer = calloc_or_null::<i32>(heap, count as u64)?;
        heap.slice_mut(unit)?[0].alist = pointer;
        if pointer.is_null() {
            OAD_PolymerUnit_Free(heap, unit)?;
            return Ok(SourceMutPointer::null());
        }
        if na > 0 {
            let source = alist.ok_or(SourceHeapError::MissingAllocation)?;
            let count = usize::try_from(na).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let values = heap.slice(source.item.as_const())?.get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?.to_vec();
            heap.slice_mut(pointer)?.get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?.copy_from_slice(&values);
        }
    }
    if nb > 0 || maxbonds > 0 {
        let count = if nb > 0 { nb.wrapping_mul(2) } else { maxbonds.wrapping_mul(2) };
        let pointer = calloc_or_null::<i32>(heap, count as u64)?;
        heap.slice_mut(unit)?[0].blist = pointer;
        if pointer.is_null() {
            OAD_PolymerUnit_Free(heap, unit)?;
            return Ok(SourceMutPointer::null());
        }
        if let Some(source) = blist {
            let copy_count = nb.wrapping_mul(2);
            if copy_count > 0 {
                let copy_count = usize::try_from(copy_count)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let values = heap.slice(source.item.as_const())?.get(..copy_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?.to_vec();
                heap.slice_mut(pointer)?.get_mut(..copy_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?.copy_from_slice(&values);
            }
        }
    }
    heap.slice_mut(unit)?[0].bkbonds = SourceMutPointer::null();
    Ok(unit)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_Free(
    heap: &mut SourceHeap,
    unit: SourceMutPointer<OAD_PolymerUnit>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1342 OAD_PolymerUnit_Free
    // INCHI✔❌: void OAD_PolymerUnit_Free( OAD_PolymerUnit *unit )
    // INCHI✔❌: {
    // INCHI✔❌:
    // INCHI✔❌:     ITRACE_( "\n************** About to free OAD_PolymerUnit @ %-p\n", unit );
    // INCHI✔❌:     OAD_PolymerUnit_DebugTrace( unit );
    // INCHI✔❌:
    // INCHI✔❌:     if (unit)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (unit->alist)
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_free( unit->alist );
    // INCHI✔❌:             unit->alist = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (unit->blist)
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_free( unit->blist );
    // INCHI✔❌:             unit->blist = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (unit->bkbonds)
    // INCHI✔❌:         {
    // INCHI✔❌:             imat_free( unit->maxbkbonds, unit->bkbonds );
    // INCHI✔❌:             unit->bkbonds = NULL;
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     inchi_free( unit );
    // INCHI✔❌:
    // INCHI✔❌:     return;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_Free
    // BEGIN INCHI ACTIVE HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.h:133 ITRACE_
    // INCHI✔❌: #define ITRACE_ 0 && _inchi_trace
    // END INCHI ACTIVE HEADER MACRO: ITRACE_

    if unit.is_null() {
        OAD_PolymerUnit_DebugTrace(None);
    } else {
        let unit_ref = heap
            .slice(unit.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        OAD_PolymerUnit_DebugTrace(Some(unit_ref));
    }

    if !unit.is_null() {
        let (alist, blist, maxbkbonds, bkbonds) = {
            let unit_ref = heap
                .slice(unit.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (
                unit_ref.alist,
                unit_ref.blist,
                unit_ref.maxbkbonds,
                unit_ref.bkbonds,
            )
        };
        if !alist.is_null() {
            inchi_free(heap, alist)?;
            heap.slice_mut(unit)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .alist = SourceMutPointer::null();
        }
        if !blist.is_null() {
            inchi_free(heap, blist)?;
            heap.slice_mut(unit)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .blist = SourceMutPointer::null();
        }
        if !bkbonds.is_null() {
            imat_free(heap, maxbkbonds, bkbonds)?;
            heap.slice_mut(unit)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .bkbonds = SourceMutPointer::null();
        }
    }

    inchi_free(heap, unit)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_CreateCopy(
    heap: &mut SourceHeap,
    unit: SourceMutPointer<OAD_PolymerUnit>,
) -> Result<SourceMutPointer<OAD_PolymerUnit>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1260 OAD_PolymerUnit_CreateCopy
    // INCHI✔️❌: OAD_PolymerUnit* OAD_PolymerUnit_CreateCopy( OAD_PolymerUnit *u )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int k, err = 0;
    // INCHI✔️❌:     OAD_PolymerUnit *u2 = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     u2 = (OAD_PolymerUnit*) inchi_calloc( 1, sizeof( OAD_PolymerUnit ) );
    // INCHI✔️❌:     if (NULL == u2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         err = 1;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     u2->id = u->id;
    // INCHI✔️❌:     u2->type = u->type;
    // INCHI✔️❌:     u2->subtype = u->subtype;
    // INCHI✔️❌:     u2->conn = u->conn;
    // INCHI✔️❌:     u2->label = u->label;
    // INCHI✔️❌:     u2->na = u->na;
    // INCHI✔️❌:     u2->nb = u->nb;
    // INCHI✔️❌:     u2->cyclizable = u->cyclizable;
    // INCHI✔️❌:     u2->cyclized = u->cyclized;
    // INCHI✔️❌:     u2->cap1_is_undef = u->cap1_is_undef;
    // INCHI✔️❌:     u2->cap2_is_undef = u->cap2_is_undef;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (k = 0; k < 4; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u2->xbr1[k] = u->xbr1[k];
    // INCHI✔️❌:         u2->xbr2[k] = u->xbr2[k];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     strcpy(u2->smt, u->smt);
    // INCHI✔️❌:
    // INCHI✔️❌:     u2->cap1 = u->cap1;
    // INCHI✔️❌:     u2->end_atom1 = u->end_atom1;
    // INCHI✔️❌:     u2->cap2 = u->cap2;
    // INCHI✔️❌:     u2->end_atom2 = u->end_atom2;
    // INCHI✔️❌:     u2->nbkbonds = u->nbkbonds;
    // INCHI✔️❌:     u2->maxbkbonds = inchi_max( u->maxbkbonds, u->nbkbonds );
    // INCHI✔️❌:
    // INCHI✔️❌:     u2->alist = (int *) inchi_calloc( u2->na, sizeof( int ) );
    // INCHI✔️❌:     if (!u2->alist)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         err = 2;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (k = 0; k < u2->na; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u2->alist[k] = u->alist[k];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     u2->blist = (int *) inchi_calloc( 2 * (long long)u2->nb, sizeof( int ) );
    // INCHI✔️❌:     if (!u2->blist)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         err = 2;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (k = 0; k < 2 * u2->nb; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u2->blist[k] = u->blist[k];
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     err = imat_new( u2->maxbkbonds, 2, &( u2->bkbonds ) );
    // INCHI✔️❌:     if (!err)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (k = 0; k < u2->nbkbonds; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             u2->bkbonds[k][0] = u->bkbonds[k][0];
    // INCHI✔️❌:             u2->bkbonds[k][1] = u->bkbonds[k][1];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     if (err)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         OAD_PolymerUnit_Free( u2 );
    // INCHI✔️❌:         return NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return u2;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_CreateCopy

    let source = heap
        .slice(unit.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let copy = match inchi_calloc::<OAD_PolymerUnit>(
        heap,
        1,
        std::mem::size_of::<OAD_PolymerUnit>() as u64,
    ) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
        Err(error) => return Err(error),
    };
    let mut copied = source.clone();
    copied.alist = SourceMutPointer::null();
    copied.blist = SourceMutPointer::null();
    copied.bkbonds = SourceMutPointer::null();
    copied.representation = 0;
    copied.maxbkbonds = source.maxbkbonds.max(source.nbkbonds);
    copied.smt = [0; 80];
    let smt_length = source
        .smt
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(source.smt.len());
    copied.smt[..smt_length].copy_from_slice(&source.smt[..smt_length]);
    if smt_length < copied.smt.len() {
        copied.smt[smt_length] = 0;
    }
    heap.slice_mut(copy)?[0] = copied;

    let fail_copy = |heap: &mut SourceHeap,
                     copy: SourceMutPointer<OAD_PolymerUnit>|
     -> Result<SourceMutPointer<OAD_PolymerUnit>, SourceHeapError> {
        OAD_PolymerUnit_Free(heap, copy)?;
        Ok(SourceMutPointer::null())
    };

    let alist_count = match u64::try_from(source.na) {
        Ok(value) => value,
        Err(_) => return fail_copy(heap, copy),
    };
    let alist = match inchi_calloc::<i32>(heap, alist_count, 4) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return fail_copy(heap, copy),
        Err(error) => return Err(error),
    };
    if !source.alist.is_null() {
        let source_values = heap.slice(source.alist.as_const())?.to_vec();
        let target = heap.slice_mut(alist)?;
        for (index, value) in target.iter_mut().enumerate() {
            *value = *source_values
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
        }
    }
    heap.slice_mut(copy)?[0].alist = alist;

    let blist_count = i64::from(source.nb)
        .checked_mul(2)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let blist_count = match u64::try_from(blist_count) {
        Ok(value) => value,
        Err(_) => return fail_copy(heap, copy),
    };
    let blist = match inchi_calloc::<i32>(heap, blist_count, 4) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return fail_copy(heap, copy),
        Err(error) => return Err(error),
    };
    if !source.blist.is_null() {
        let source_values = heap.slice(source.blist.as_const())?.to_vec();
        let target = heap.slice_mut(blist)?;
        for (index, value) in target.iter_mut().enumerate() {
            *value = *source_values
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
        }
    }
    heap.slice_mut(copy)?[0].blist = blist;

    let max_bonds = heap.slice(copy.as_const())?[0].maxbkbonds;
    let mut backbone_bonds = heap.slice(copy.as_const())?[0].bkbonds;
    let bond_error = imat_new(heap, max_bonds, 2, &mut backbone_bonds)?;
    heap.slice_mut(copy)?[0].bkbonds = backbone_bonds;
    if bond_error != 0 {
        return fail_copy(heap, copy);
    }
    if !source.bkbonds.is_null() {
        for index in 0..source.nbkbonds {
            let source_row = heap
                .slice(source.bkbonds.as_const())?
                .get(index as usize)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let target_rows = heap.slice(copy.as_const())?[0].bkbonds;
            let target_row = heap
                .slice(target_rows.as_const())?
                .get(index as usize)
                .copied()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let values = heap
                .slice(source_row.as_const())?
                .get(..2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let values = values.to_vec();
            heap.slice_mut(target_row)?
                .get_mut(..2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .copy_from_slice(&values);
        }
    }
    Ok(copy)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_CompareAtomListsMod(
    heap: &SourceHeap,
    first: SourceConstPointer<OAD_PolymerUnit>,
    second: SourceConstPointer<OAD_PolymerUnit>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1377 OAD_PolymerUnit_CompareAtomListsMod
    // INCHI❌❌: int  OAD_PolymerUnit_CompareAtomListsMod( OAD_PolymerUnit *u1,
    // INCHI❌❌:                                           OAD_PolymerUnit *u2 )
    // INCHI❌❌: {
    // INCHI❌❌:     int i;
    // INCHI❌❌:     int n1 = u1->na;
    // INCHI❌❌:     int n2 = u2->na;
    // INCHI❌❌:     int n = n1;
    // INCHI❌❌:     if (n1 < n2)    return -1;
    // INCHI❌❌:     if (n1 > n2)    return 1;
    // INCHI❌❌:     /* n1 == n2 == n */
    // INCHI❌❌:     for (i = 0; i < n; i++)
    // INCHI❌❌:     {
    // INCHI❌❌:         if (u1->alist[i] < u2->alist[i])    return -1;
    // INCHI❌❌:         if (u1->alist[i] > u2->alist[i])    return    1;
    // INCHI❌❌:     }
    // INCHI❌❌:
    // INCHI❌❌:     return 0;
    // INCHI❌❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_CompareAtomListsMod

    let first = heap
        .slice(first)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second = heap
        .slice(second)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if first.na < second.na {
        return Ok(-1);
    }
    if first.na > second.na {
        return Ok(1);
    }
    if first.na <= 0 {
        return Ok(0);
    }
    let count = usize::try_from(first.na).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let first_values = heap.slice(first.alist.as_const())?;
    let second_values = heap.slice(second.alist.as_const())?;
    for index in 0..count {
        if first_values[index] < second_values[index] {
            return Ok(-1);
        }
        if first_values[index] > second_values[index] {
            return Ok(1);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
    heap: &mut SourceHeap,
    unit: &mut OAD_PolymerUnit,
    number_of_star_atoms: i32,
    star_atoms: SourceConstPointer<i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1437 OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves
    // INCHI❌❌: int  OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves( OAD_PolymerUnit  *u,
    // INCHI❌❌:                                                        int n_star_atoms,
    // INCHI❌❌:                                                        int *star_atoms )
    // INCHI❌❌: {
    // INCHI❌❌:     int k;
    // INCHI❌❌:
    // INCHI❌❌:     /* Sort bond atoms */
    // INCHI❌❌:     for (k = 0; k < u->nb; k++)
    // INCHI❌❌:     {
    // INCHI❌❌:         /* Place not-in-unit bond end to first place */
    // INCHI❌❌:         int a1 = u->blist[2 * k];
    // INCHI❌❌:         int a2 = u->blist[2 * k + 1];
    // INCHI❌❌:         int a1_is_not_in_alist = 0;
    // INCHI❌❌:         int a1_is_star_atom = 0;
    // INCHI❌❌:         int a2_is_not_in_alist = 0;
    // INCHI❌❌:         int a2_is_star_atom = 0;
    // INCHI❌❌:
    // INCHI❌❌:         if (!is_in_the_ilist( u->alist, a1, u->na ))
    // INCHI❌❌:         {
    // INCHI❌❌:             a1_is_not_in_alist = 1;
    // INCHI❌❌:         }
    // INCHI❌❌:         if (is_in_the_ilist( star_atoms, a1, n_star_atoms ))
    // INCHI❌❌:         {
    // INCHI❌❌:             a1_is_star_atom = 1;
    // INCHI❌❌:         }
    // INCHI❌❌:
    // INCHI❌❌:         if (!is_in_the_ilist( u->alist, a2, u->na ))
    // INCHI❌❌:         {
    // INCHI❌❌:             a2_is_not_in_alist = 1;
    // INCHI❌❌:         }
    // INCHI❌❌:         if (is_in_the_ilist( star_atoms, a2, n_star_atoms ))
    // INCHI❌❌:         {
    // INCHI❌❌:             a2_is_star_atom = 1;
    // INCHI❌❌:         }
    // INCHI❌❌:
    // INCHI❌❌:         if (( a1_is_not_in_alist || a1_is_star_atom ) &&
    // INCHI❌❌:             ( a2_is_not_in_alist || a2_is_star_atom ))
    // INCHI❌❌:         {
    // INCHI❌❌:             /* Both the ends are out of unit: the crossing bond is invalid */
    // INCHI❌❌:             return 1;
    // INCHI❌❌:         }
    // INCHI❌❌:         /* If a2 is star atom or non-star external to the current unit, swap(a2,a1) */
    // INCHI❌❌:         if (a2_is_star_atom || a2_is_not_in_alist)
    // INCHI❌❌:         {
    // INCHI❌❌:             u->blist[2 * k] = a2;
    // INCHI❌❌:             u->blist[2 * k + 1] = a1;
    // INCHI❌❌:         }
    // INCHI❌❌:     }
    // INCHI❌❌:
    // INCHI❌❌:     /* Sort bond themselves
    // INCHI❌❌:         for now, consider only the simplest cases of 2 bonds
    // INCHI❌❌:     */
    // INCHI❌❌:     if (u->nb == 2)            /* two bonds in SBL */
    // INCHI❌❌:     {
    // INCHI❌❌:         int b1a1 = u->blist[0];
    // INCHI❌❌:         int b1a2 = u->blist[1];
    // INCHI❌❌:         int b2a1 = u->blist[2];
    // INCHI❌❌:         int b2a2 = u->blist[3];
    // INCHI❌❌:         if (b1a1 > b2a1)
    // INCHI❌❌:         {
    // INCHI❌❌:             /* swap */
    // INCHI❌❌:             u->blist[0] = b2a1; u->blist[1] = b2a2;
    // INCHI❌❌:             u->blist[2] = b1a1; u->blist[3] = b1a2;
    // INCHI❌❌:         }
    // INCHI❌❌:     }
    // INCHI❌❌:
    // INCHI❌❌:     /* for single or no bonds, do nothing
    // INCHI❌❌:     else
    // INCHI❌❌:         ;
    // INCHI❌❌:     */
    // INCHI❌❌:
    // INCHI❌❌:     return 0;
    // INCHI❌❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves

    if unit.nb <= 0 {
        return Ok(0);
    }
    let alist = if unit.na == 0 {
        None
    } else {
        Some(heap.slice(unit.alist.as_const())?.to_vec())
    };
    let stars = if number_of_star_atoms == 0 {
        None
    } else {
        Some(heap.slice(star_atoms)?.to_vec())
    };
    for bond_index in 0..unit.nb {
        let offset = usize::try_from(i64::from(bond_index) * 2)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let (first, second) = {
            let bonds = heap.slice(unit.blist.as_const())?;
            (
                *bonds
                    .get(offset)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                *bonds
                    .get(offset + 1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
        };
        let first_external = is_in_the_ilist(alist.as_deref(), first, unit.na)?.is_none();
        let first_star = is_in_the_ilist(stars.as_deref(), first, number_of_star_atoms)?.is_some();
        let second_external = is_in_the_ilist(alist.as_deref(), second, unit.na)?.is_none();
        let second_star =
            is_in_the_ilist(stars.as_deref(), second, number_of_star_atoms)?.is_some();
        if (first_external || first_star) && (second_external || second_star) {
            return Ok(1);
        }
        if second_star || second_external {
            let bonds = heap.slice_mut(unit.blist)?;
            bonds[offset] = second;
            bonds[offset + 1] = first;
        }
    }
    if unit.nb == 2 {
        let bonds = heap.slice_mut(unit.blist)?;
        let [first_first, first_second, second_first, second_second] =
            [bonds[0], bonds[1], bonds[2], bonds[3]];
        if first_first > second_first {
            bonds[0] = second_first;
            bonds[1] = second_second;
            bonds[2] = first_first;
            bonds[3] = first_second;
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_PrepareWorkingSet(
    heap: &mut SourceHeap,
    polymer: &mut OAD_Polymer,
    canonical_numbers: SourceConstPointer<i32>,
    _component_numbers: SourceConstPointer<i32>,
    units2: SourceMutPointer<SourceMutPointer<OAD_PolymerUnit>>,
    unit_numbers: SourceMutPointer<i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2342 OAD_Polymer_PrepareWorkingSet
    // INCHI✔️❌: int OAD_Polymer_PrepareWorkingSet( OAD_Polymer     *p,
    // INCHI✔️❌:                                    int             *cano_nums,
    // INCHI✔️❌:                                    int             *compnt_nums,
    // INCHI✔️❌:                                    OAD_PolymerUnit **units2,
    // INCHI✔️❌:                                    int             *unum )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, k, err = 0, cano_num1 = -1, cano_num2 = -1;
    // INCHI✔️❌:     OAD_PolymerUnit *u;
    // INCHI✔️❌:
    // INCHI✔️❌:     /*OAD_Polymer_DebugTrace( p );*/
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Replace original atom numbers in polymer data with canonical plus 1.                    */
    // INCHI✔️❌:     /*  Note that we use 'cano1 nums', that is, 1-based (InChI internal 'cano nums' are 0-based)*/
    // INCHI✔️❌:     /*  Also remove from the list atoms who mapped to cano number 0  ( == -1 + 1_offset ),      */
    // INCHI✔️❌:     /*  they are explicit H's which have already been deleted.                                  */
    // INCHI✔️❌:     for (k = 0; k < p->n_pzz; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         cano_num1 = cano_nums[p->pzz[k]] + 1;
    // INCHI✔️❌:         if (cano_num1 == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* we shouldn't arrive here */
    // INCHI✔️❌:             err = 10;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         p->pzz[k] = cano_num1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < p->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int na_new = -1;
    // INCHI✔️❌:         u = units2[i];
    // INCHI✔️❌:
    // INCHI✔️❌:         for (k = 0; k < u->na; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cano_num1 = cano_nums[u->alist[k]] + 1;
    // INCHI✔️❌:             if (cano_num1 == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             u->alist[++na_new] = cano_num1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         u->na = na_new + 1;
    // INCHI✔️❌:         for (k = 0; k < 2 * u->nb; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cano_num1 = cano_nums[u->blist[k]] + 1;
    // INCHI✔️❌:             if (cano_num1 == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Can not proceed further as one of PU crossing bond ends
    // INCHI✔️❌:                    leads to explicit H which has been removed already       */
    // INCHI✔️❌:                 err = 11;
    // INCHI✔️❌:                 goto exit_function;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             u->blist[k] = cano_num1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         cano_num1 = cano_nums[u->cap1] + 1;
    // INCHI✔️❌:         if (cano_num1 == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             err = 11;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         u->cap1 = cano_num1;
    // INCHI✔️❌:
    // INCHI✔️❌:         cano_num1 = cano_nums[u->cap2] + 1;
    // INCHI✔️❌:         if (cano_num1 == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             err = 11;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         u->cap2 = cano_num1;
    // INCHI✔️❌:
    // INCHI✔️❌:         cano_num1 = cano_nums[u->end_atom1] + 1;
    // INCHI✔️❌:         if (cano_num1 == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             err = 11;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         u->end_atom1 = cano_num1;
    // INCHI✔️❌:
    // INCHI✔️❌:         cano_num1 = cano_nums[u->end_atom2] + 1;
    // INCHI✔️❌:         if (cano_num1 == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             err = 11;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         u->end_atom2 = cano_num1;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (k = 0; k < u->nbkbonds; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cano_num1 = cano_nums[u->bkbonds[k][0]] + 1;
    // INCHI✔️❌:             if (cano_num1 == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             cano_num2 = cano_nums[u->bkbonds[k][1]] + 1;
    // INCHI✔️❌:             if (cano_num2 == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             u->bkbonds[k][0] = inchi_min( cano_num1, cano_num2 );
    // INCHI✔️❌:             u->bkbonds[k][1] = inchi_max( cano_num1, cano_num2 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Sort the atoms and the bonds in all units */
    // INCHI✔️❌:     for (i = 0; i < p->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u = units2[i];
    // INCHI✔️❌:
    // INCHI✔️❌:         /* sort atoms (alist) */
    // INCHI✔️❌:         iisort( u->alist, u->na );
    // INCHI✔️❌:
    // INCHI✔️❌:         /*ITRACE_( "\n*** Polymer unit %-d : ( ", i );
    // INCHI✔️❌:         for (k = 0; k < u->na - 1; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ITRACE_( "%-d-", u->alist[k] );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ITRACE_( "%-d )\n", u->alist[u->na - 1] );*/
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Sort bonds (blist) */
    // INCHI✔️❌:         err = OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves( u, p->n_pzz, p->pzz );
    // INCHI✔️❌:         if (err)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* crossing bonds in blist are invalid */
    // INCHI✔️❌:             err = 12;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Check each unit for >1 connected components */
    // INCHI✔️❌: #if 0
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int icompnt;
    // INCHI✔️❌:             icompnt = compnt_nums[u->alist[0] - 1];
    // INCHI✔️❌:             for (k = 1; k < u->na; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (compnt_nums[u->alist[k] - 1] != icompnt)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     u->disjoint = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Sort all units in modified alist's lexicographic order
    // INCHI✔️❌:     (modification is: longer list always go first )\t\t\t*/
    // INCHI✔️❌:     for (i = 0; i < p->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         unum[i] = i;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 1; i < p->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int tmp = unum[i];
    // INCHI✔️❌:         int j = i - 1;
    // INCHI✔️❌:         while (j >= 0 && OAD_PolymerUnit_CompareAtomListsMod( units2[unum[j]], units2[tmp] ) > 0)
    // INCHI✔️❌:         /*while ( j >= 0 &&    OAD_PolymerUnit_CompareAtomLists( units2[ unum[j] ], units2[ tmp ] ) > 0  )*/
    // INCHI✔️❌:         {
    // INCHI✔️❌:             unum[j + 1] = unum[j];
    // INCHI✔️❌:             j--;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         unum[j + 1] = tmp;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     return err;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_Polymer_PrepareWorkingSet

    let canonical_values = if polymer.n_pzz > 0 || polymer.n > 0 {
        heap.slice(canonical_numbers)?.to_vec()
    } else {
        Vec::new()
    };
    let map = |number: i32| -> Result<i32, SourceHeapError> {
        let index = usize::try_from(number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        Ok(canonical_values
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .wrapping_add(1))
    };
    for index in 0..polymer.n_pzz {
        let pointer = polymer.pzz.offset(i64::from(index))?;
        let original = *heap
            .slice(pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mapped = map(original)?;
        if mapped == 0 {
            return Ok(10);
        }
        heap.slice_mut(pointer)?[0] = mapped;
    }

    for index in 0..polymer.n {
        let unit_pointer = *heap
            .slice(units2.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let mut kept = 0_i32;
        for atom_index in 0..unit.na {
            let atom_pointer = unit.alist.offset(i64::from(atom_index))?;
            let original = *heap
                .slice(atom_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mapped = map(original)?;
            if mapped != 0 {
                heap.slice_mut(unit.alist.offset(i64::from(kept))?)?[0] = mapped;
                kept = kept.wrapping_add(1);
            }
        }
        heap.slice_mut(unit_pointer)?[0].na = kept;
        for bond_index in 0..unit.nb.wrapping_mul(2) {
            let pointer = unit.blist.offset(i64::from(bond_index))?;
            let original = *heap
                .slice(pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mapped = map(original)?;
            if mapped == 0 {
                return Ok(11);
            }
            heap.slice_mut(pointer)?[0] = mapped;
        }
        for field_index in 0..4 {
            let original = {
                let current = &heap.slice(unit_pointer.as_const())?[0];
                match field_index {
                    0 => current.cap1,
                    1 => current.cap2,
                    2 => current.end_atom1,
                    _ => current.end_atom2,
                }
            };
            let mapped = map(original)?;
            if mapped == 0 {
                return Ok(11);
            }
            let current = &mut heap.slice_mut(unit_pointer)?[0];
            match field_index {
                0 => current.cap1 = mapped,
                1 => current.cap2 = mapped,
                2 => current.end_atom1 = mapped,
                _ => current.end_atom2 = mapped,
            }
        }
        for bond_index in 0..unit.nbkbonds {
            let row_pointer = *heap
                .slice(unit.bkbonds.as_const().offset(i64::from(bond_index))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let first = *heap
                .slice(row_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let first_mapped = map(first)?;
            if first_mapped == 0 {
                continue;
            }
            let second = *heap
                .slice(row_pointer.as_const().offset(1)?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let second_mapped = map(second)?;
            if second_mapped == 0 {
                continue;
            }
            let row = heap.slice_mut(row_pointer)?;
            row[0] = first_mapped.min(second_mapped);
            row[1] = first_mapped.max(second_mapped);
        }
    }

    for index in 0..polymer.n {
        let unit_pointer = *heap
            .slice(units2.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        iisort(heap, unit.alist, unit.na)?;
        let mut sorted_unit = unit;
        let order_error = OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
            heap,
            &mut sorted_unit,
            polymer.n_pzz,
            polymer.pzz.as_const(),
        )?;
        if order_error != 0 {
            return Ok(12);
        }
        heap.slice_mut(unit_pointer)?[0] = sorted_unit;
    }

    for index in 0..polymer.n {
        heap.slice_mut(unit_numbers.offset(i64::from(index))?)?[0] = index;
    }
    for index in 1..polymer.n {
        let temporary = *heap
            .slice(unit_numbers.offset(i64::from(index))?.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut previous = index - 1;
        while previous >= 0 {
            let left_number = *heap
                .slice(unit_numbers.offset(i64::from(previous))?.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let left_pointer = *heap
                .slice(units2.offset(i64::from(left_number))?.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let right_pointer = *heap
                .slice(units2.offset(i64::from(temporary))?.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if OAD_PolymerUnit_CompareAtomListsMod(
                heap,
                left_pointer.as_const(),
                right_pointer.as_const(),
            )? <= 0
            {
                break;
            }
            heap.slice_mut(unit_numbers.offset(i64::from(previous + 1))?)?[0] = left_number;
            previous -= 1;
        }
        heap.slice_mut(unit_numbers.offset(i64::from(previous + 1))?)?[0] = temporary;
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_RemoveHalfBond(
    heap: &mut SourceHeap,
    this_atom: i32,
    other_atom: i32,
    atoms: SourceMutPointer<inp_ATOM>,
    bond_type: &mut i32,
    bond_stereo: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2517 OrigAtData_RemoveHalfBond
    // INCHI✔️✔️: int  OrigAtData_RemoveHalfBond( int      this_atom,
    // INCHI✔️✔️:                                 int      other_atom,
    // INCHI✔️✔️:                                 inp_ATOM *at,
    // INCHI✔️✔️:                                 int      *bond_type,
    // INCHI✔️✔️:                                 int      *bond_stereo )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int k, kk;
    // INCHI✔️✔️:     /* djb-rwth: fixing oss-fuzz issues #68286, #30342 */
    // INCHI✔️✔️:     if (at && (this_atom >= 0) && (other_atom >= 0))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inp_ATOM* a = &(at[this_atom]);
    // INCHI✔️✔️:         if (a)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             for (k = 0; k < a->valence; k++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (a->neighbor[k] != other_atom)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     continue;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 *bond_type = a->bond_type[k];
    // INCHI✔️✔️:                 *bond_stereo = a->bond_stereo[k];
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 a->neighbor[k] = a->bond_type[k] = a->bond_stereo[k] = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 for (kk = k + 1; kk < a->valence; kk++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     a->neighbor[kk - 1] = a->neighbor[kk];
    // INCHI✔️✔️:                     a->bond_type[kk - 1] = a->bond_type[kk];
    // INCHI✔️✔️:                     a->bond_stereo[kk - 1] = a->bond_stereo[kk];
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 for (kk = a->valence - 1; kk < MAXVAL; kk++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     a->neighbor[kk] = 0;
    // INCHI✔️✔️:                     a->bond_type[kk] = (U_CHAR)0;
    // INCHI✔️✔️:                     a->bond_stereo[kk] = (S_CHAR)0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 return 1;
    // INCHI✔️✔️:             } /* k */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OrigAtData_RemoveHalfBond

    if atoms.is_null() || this_atom < 0 || other_atom < 0 {
        return Ok(0);
    }
    let atom_pointer = atoms.offset(i64::from(this_atom))?;
    let atom = heap
        .slice_mut(atom_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence = i32::from(atom.valence);
    for k in 0..valence {
        let index = usize::try_from(k).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if i32::from(
            *atom
                .neighbor
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ) != other_atom
        {
            continue;
        }
        *bond_type = i32::from(atom.bond_type[index]);
        *bond_stereo = i32::from(atom.bond_stereo[index]);
        atom.neighbor[index] = 0;
        atom.bond_type[index] = 0;
        atom.bond_stereo[index] = 0;
        for kk in (k + 1)..valence {
            let source = usize::try_from(kk).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let target = source - 1;
            atom.neighbor[target] = atom.neighbor[source];
            atom.bond_type[target] = atom.bond_type[source];
            atom.bond_stereo[target] = atom.bond_stereo[source];
        }
        for kk in (valence - 1)..20 {
            let index = usize::try_from(kk).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            atom.neighbor[index] = 0;
            atom.bond_type[index] = 0;
            atom.bond_stereo[index] = 0;
        }
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_RemoveBond(
    heap: &mut SourceHeap,
    this_atom: i32,
    other_atom: i32,
    atoms: SourceMutPointer<inp_ATOM>,
    bond_type: &mut i32,
    bond_stereo: &mut i32,
    num_inp_bonds: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2577 OrigAtData_RemoveBond
    // INCHI✔️✔️: int  OrigAtData_RemoveBond( int      this_atom,
    // INCHI✔️✔️:                             int      other_atom,
    // INCHI✔️✔️:                             inp_ATOM *at,
    // INCHI✔️✔️:                             int      *bond_type,
    // INCHI✔️✔️:                             int      *bond_stereo,
    // INCHI✔️✔️:                             int      *num_inp_bonds )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int del = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (at && (this_atom >= 0) && (other_atom >= 0)) /* djb-rwth: fixing oss-fuzz issue #68329, #68286 */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         del = OrigAtData_RemoveHalfBond(this_atom, other_atom, at, bond_type, bond_stereo);
    // INCHI✔️✔️:         del += OrigAtData_RemoveHalfBond(other_atom, this_atom, at, bond_type, bond_stereo);
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (del == 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             (*num_inp_bonds)--;
    // INCHI✔️✔️:             at[this_atom].valence--;
    // INCHI✔️✔️:             at[this_atom].chem_bonds_valence -= *bond_type;
    // INCHI✔️✔️:             at[other_atom].valence--;
    // INCHI✔️✔️:             at[other_atom].chem_bonds_valence -= *bond_type;
    // INCHI✔️✔️:             return 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OrigAtData_RemoveBond

    if atoms.is_null() || this_atom < 0 || other_atom < 0 {
        return Ok(0);
    }
    let first =
        OrigAtData_RemoveHalfBond(heap, this_atom, other_atom, atoms, bond_type, bond_stereo)?;
    let second =
        OrigAtData_RemoveHalfBond(heap, other_atom, this_atom, atoms, bond_type, bond_stereo)?;
    if first + second == 2 {
        *num_inp_bonds = num_inp_bonds.wrapping_sub(1);
        let first_pointer = atoms.offset(i64::from(this_atom))?;
        let second_pointer = atoms.offset(i64::from(other_atom))?;
        {
            let first_atom = heap
                .slice_mut(first_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            first_atom.valence = first_atom.valence.wrapping_sub(1);
            first_atom.chem_bonds_valence =
                first_atom.chem_bonds_valence.wrapping_sub(*bond_type as i8);
        }
        {
            let second_atom = heap
                .slice_mut(second_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            second_atom.valence = second_atom.valence.wrapping_sub(1);
            second_atom.chem_bonds_valence = second_atom
                .chem_bonds_valence
                .wrapping_sub(*bond_type as i8);
        }
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_AddBond(
    heap: &mut SourceHeap,
    this_atom: i32,
    other_atom: i32,
    atoms: SourceMutPointer<inp_ATOM>,
    mut bond_type: i32,
    bond_stereo: i32,
    num_bonds: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2607 OrigAtData_AddBond
    // INCHI✔️✔️: int  OrigAtData_AddBond( int        this_atom,
    // INCHI✔️✔️:                          int        other_atom,
    // INCHI✔️✔️:                          inp_ATOM   *at,
    // INCHI✔️✔️:                          int        bond_type,
    // INCHI✔️✔️:                          int        bond_stereo,
    // INCHI✔️✔️:                          int        *num_bonds )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (at)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* djb-rwth: fixing oss-fuzz issue #68286 */
    // INCHI✔️✔️:         int i, k, already_here;
    // INCHI✔️✔️:         inp_ATOM* a = &(at[this_atom]);
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (at[this_atom].valence >= MAXVAL ||
    // INCHI✔️✔️:             at[other_atom].valence >= MAXVAL)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (bond_type != INCHI_BOND_TYPE_DOUBLE && bond_type != INCHI_BOND_TYPE_TRIPLE)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             bond_type = INCHI_BOND_TYPE_SINGLE;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         k = a->valence;
    // INCHI✔️✔️:         already_here = 0;
    // INCHI✔️✔️:         for (i = 0; i < k; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (a->neighbor[i] == other_atom)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 already_here = 1; break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (!already_here)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a->neighbor[k] = other_atom;
    // INCHI✔️✔️:             a->bond_type[k] = (U_CHAR)bond_type;
    // INCHI✔️✔️:             a->bond_stereo[k] = (S_CHAR)bond_stereo;
    // INCHI✔️✔️:             a->chem_bonds_valence += bond_type;
    // INCHI✔️✔️:             a->valence++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         a = &(at[other_atom]);
    // INCHI✔️✔️:         k = a->valence;
    // INCHI✔️✔️:         already_here = 0;
    // INCHI✔️✔️:         for (i = 0; i < k; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (a->neighbor[i] == this_atom)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 already_here = 1; break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (!already_here && (k < MAXVAL)) /* djb-rwth: condition added to prevent buffer overrun */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a->neighbor[k] = this_atom;
    // INCHI✔️✔️:             a->bond_type[k] = (U_CHAR)bond_type;
    // INCHI✔️✔️:             a->bond_stereo[k] = (S_CHAR)bond_stereo;
    // INCHI✔️✔️:             a->chem_bonds_valence += bond_type;
    // INCHI✔️✔️:             a->valence++;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         (*num_bonds)++;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OrigAtData_AddBond

    if atoms.is_null() {
        return Ok(0);
    }
    let first_pointer = atoms.offset(i64::from(this_atom))?;
    let second_pointer = atoms.offset(i64::from(other_atom))?;
    let first_valence = heap
        .slice(first_pointer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .valence;
    let second_valence = heap
        .slice(second_pointer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .valence;
    if first_valence >= 20 || second_valence >= 20 {
        return Ok(0);
    }
    if bond_type != tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE as i32
        && bond_type != tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE as i32
    {
        bond_type = tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i32;
    }
    {
        let atom = heap
            .slice_mut(first_pointer)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let k =
            usize::try_from(atom.valence).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let already_here = atom.neighbor[..k]
            .iter()
            .any(|neighbor| i32::from(*neighbor) == other_atom);
        if !already_here {
            atom.neighbor[k] = other_atom as u16;
            atom.bond_type[k] = bond_type as u8;
            atom.bond_stereo[k] = bond_stereo as i8;
            atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_add(bond_type as i8);
            atom.valence = atom.valence.wrapping_add(1);
        }
    }
    {
        let atom = heap
            .slice_mut(second_pointer)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let k =
            usize::try_from(atom.valence).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let already_here = atom.neighbor[..k]
            .iter()
            .any(|neighbor| i32::from(*neighbor) == this_atom);
        if !already_here && k < 20 {
            atom.neighbor[k] = this_atom as u16;
            atom.bond_type[k] = bond_type as u8;
            atom.bond_stereo[k] = bond_stereo as i8;
            atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_add(bond_type as i8);
            atom.valence = atom.valence.wrapping_add(1);
        }
    }
    *num_bonds = num_bonds.wrapping_add(1);
    Ok(1)
}

#[allow(non_snake_case)]
pub(crate) fn UnMarkRingSystemsInp(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1961 UnMarkRingSystemsInp
    // INCHI✔️✔️: int UnMarkRingSystemsInp( inp_ATOM *at, int num_atoms )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int i;
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         at[i].bCutVertex = 0;
    // INCHI✔️✔️:         at[i].nRingSystem = 0;
    // INCHI✔️✔️:         at[i].nNumAtInRingSystem = 0;
    // INCHI✔️✔️:         at[i].nBlockSystem = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: UnMarkRingSystemsInp

    for index in 0..num_atoms {
        let atom = heap
            .slice_mut(atoms.offset(i64::from(index))?)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        atom.bCutVertex = 0;
        atom.nRingSystem = 0;
        atom.nNumAtInRingSystem = 0;
        atom.nBlockSystem = 0;
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_AddSingleStereolessBond(
    heap: &mut SourceHeap,
    this_atom: i32,
    other_atom: i32,
    atoms: SourceMutPointer<inp_ATOM>,
    num_bonds: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2682 OrigAtData_AddSingleStereolessBond
    // INCHI✔️✔️: int  OrigAtData_AddSingleStereolessBond( int      this_atom,
    // INCHI✔️✔️:                                          int      other_atom,
    // INCHI✔️✔️:                                          inp_ATOM *at,
    // INCHI✔️✔️:                                          int      *num_bonds )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     return OrigAtData_AddBond( this_atom, other_atom, at, INCHI_BOND_TYPE_SINGLE, 0, num_bonds );
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OrigAtData_AddSingleStereolessBond

    OrigAtData_AddBond(
        heap,
        this_atom,
        other_atom,
        atoms,
        tagINCHIBondType_INCHI_BOND_TYPE_SINGLE as i32,
        0,
        num_bonds,
    )
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_FindRingSystems(
    heap: &mut SourceHeap,
    polymer: SourceConstPointer<OAD_Polymer>,
    atoms: SourceMutPointer<inp_ATOM>,
    nat: i32,
    num_inp_bonds: &mut i32,
    num_ring_sys: SourceMutPointer<i32>,
    size_ring_sys: SourceMutPointer<i32>,
    start: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3236 OAD_Polymer_FindRingSystems
    // INCHI✔️❌: int  OAD_Polymer_FindRingSystems( OAD_Polymer *pd,
    // INCHI✔️❌:                                   inp_ATOM    *at,
    // INCHI✔️❌:                                   int         nat,
    // INCHI✔️❌:                                   int         *num_inp_bonds,
    // INCHI✔️❌:                                   int         *num_ring_sys,
    // INCHI✔️❌:                                   int         *size_ring_sys,
    // INCHI✔️❌:                                   int         start )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, nrings = 0, bond_type, bond_stereo;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (NULL == num_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Remove polymer SRU 'cyclizing' bonds if any */
    // INCHI✔️❌:     for (j = 0; j < pd->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pd->units[j]->cyclized)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             OrigAtData_RemoveBond( pd->units[j]->end_atom1 - 1,
    // INCHI✔️❌:                                    pd->units[j]->end_atom2 - 1,
    // INCHI✔️❌:                                    at, &bond_type, &bond_stereo,
    // INCHI✔️❌:                                    num_inp_bonds );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     MarkRingSystemsInp( at, nat, start ); /*0 );*/
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i <= nat; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_ring_sys[i] = -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 0; i < nat; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[i].nNumAtInRingSystem > 2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int atnum = at[i].orig_at_number;
    // INCHI✔️❌:             num_ring_sys[atnum] = at[i].nRingSystem;
    // INCHI✔️❌:             if (NULL != size_ring_sys)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 size_ring_sys[atnum] = at[i].nNumAtInRingSystem;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     UnMarkRingSystemsInp( at, nat );
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < nat; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (num_ring_sys[i] > -1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nrings++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Restore polymer SRU 'cyclizing' bonds if applicable */
    // INCHI✔️❌:     for (j = 0; j < pd->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pd->units[j]->cyclized)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             OrigAtData_AddSingleStereolessBond( pd->units[j]->end_atom1 - 1,
    // INCHI✔️❌:                                                 pd->units[j]->end_atom2 - 1,
    // INCHI✔️❌:                                                 at,
    // INCHI✔️❌:                                                 num_inp_bonds );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return nrings;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_Polymer_FindRingSystems

    if num_ring_sys.is_null() {
        return Ok(0);
    }
    let polymer_value = heap
        .slice(polymer)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let mut bond_type = 0;
    let mut bond_stereo = 0;
    for index in 0..polymer_value.n {
        let unit_pointer = *heap
            .slice(polymer_value.units.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        if unit.cyclized != 0 {
            let _ = OrigAtData_RemoveBond(
                heap,
                unit.end_atom1.wrapping_sub(1),
                unit.end_atom2.wrapping_sub(1),
                atoms,
                &mut bond_type,
                &mut bond_stereo,
                num_inp_bonds,
            )?;
        }
    }
    let _ = MarkRingSystemsInp(heap, atoms, nat, start)?;
    for index in 0..=nat {
        heap.slice_mut(num_ring_sys.offset(i64::from(index))?)?[0] = -1;
    }
    for index in 0..nat {
        let atom = heap
            .slice(atoms.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        if atom.nNumAtInRingSystem > 2 {
            let atom_index = i32::from(atom.orig_at_number);
            heap.slice_mut(num_ring_sys.offset(i64::from(atom_index))?)?[0] =
                i32::from(atom.nRingSystem);
            if !size_ring_sys.is_null() {
                heap.slice_mut(size_ring_sys.offset(i64::from(atom_index))?)?[0] =
                    i32::from(atom.nNumAtInRingSystem);
            }
        }
    }
    UnMarkRingSystemsInp(heap, atoms, nat)?;
    let mut rings = 0;
    for index in 0..nat {
        if heap
            .slice(num_ring_sys.as_const().offset(i64::from(index))?)?
            .first()
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            > -1
        {
            rings += 1;
        }
    }
    for index in 0..polymer_value.n {
        let unit_pointer = *heap
            .slice(polymer_value.units.as_const().offset(i64::from(index))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        if unit.cyclized != 0 {
            let _ = OrigAtData_AddSingleStereolessBond(
                heap,
                unit.end_atom1.wrapping_sub(1),
                unit.end_atom2.wrapping_sub(1),
                atoms,
                num_inp_bonds,
            )?;
        }
    }
    Ok(rings)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_SetAtProps(
    heap: &mut SourceHeap,
    polymer: SourceConstPointer<OAD_Polymer>,
    atoms: SourceMutPointer<inp_ATOM>,
    nat: i32,
    num_inp_bonds: &mut i32,
    atom_properties: SourceMutPointer<OAD_AtProps>,
    canonical_numbers: SourceConstPointer<i32>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3312 OAD_Polymer_SetAtProps
    // INCHI✔️❌: void OAD_Polymer_SetAtProps( OAD_Polymer *pd,
    // INCHI✔️❌:                              inp_ATOM    *at,
    // INCHI✔️❌:                              int         nat,
    // INCHI✔️❌:                              int         *num_inp_bonds,
    // INCHI✔️❌:                              OAD_AtProps *aprops,
    // INCHI✔️❌:                              int         *cano_nums )
    // INCHI✔️❌: {
    // INCHI✔️❌: /*  Max rank for in-ring atom is 216 which is achieved for N (element number 7 in Periodic system & erank_rule2[] ),*/
    // INCHI✔️❌: /*    then goes O with rank 215 (element number 8), and so on... lowest rank is 1 for H .                           */
    // INCHI✔️❌: /*                                                                                                                  */
    // INCHI✔️❌: /*  This follows to IUPAC rule 2 [Pure Appl. Chem., Vol. 74, No. 10, 2002, p. 1926] which states:                   */
    // INCHI✔️❌: /*  a. a ring or ring system containing nitrogen;                                                                   */
    // INCHI✔️❌: /*  b. a ring or ring system containing the heteroatom occurring earliest in the order given in Rule 4;             */
    // INCHI✔️❌: /*  ( which is     O > S > Se > Te > N > P > As > Sb > Bi > Si > Ge > Sn > Pb > B > Hg )                            */
    // INCHI✔️❌: /*  ...                                                                                                             */
    // INCHI✔️❌:
    // INCHI✔️❌:     int erank_rule2[] = { 0,1,198,197,196,202,2,216,215,191,190,189,188,187,206,210,214,183,182,181,180,179,178,177,176,
    // INCHI✔️❌:                           175,174,173,172,171,170,169,205,209,213,165,164,163,162,161,160,159,158,157,156,155,154,153,152,
    // INCHI✔️❌:                           151,204,208,212,147,146,145,144,143,142,141,140,139,138,137,136,135,134,133,132,131,130,129,128,
    // INCHI✔️❌:                           127,126,125,124,123,122,121,201,119,203,207,116,115,114,113,112,111,110,109,108,107,106,105,104,
    // INCHI✔️❌:                           103,102,101,100,99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81 };
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Max rank for chain atom is 215 which is achieved for O (element number 8 in Periodic system & erank_rule4[] ),  */
    // INCHI✔️❌:     /*  then goes N with rank 212 (element number 8), and so on... lowest rank is 1 for H .                             */
    // INCHI✔️❌:     /*                                                                                                                  */
    // INCHI✔️❌:     /*  This follows to IUPAC rule 4 [Pure Appl. Chem., Vol. 74, No. 10, 2002, p. 1927] which states:                   */
    // INCHI✔️❌:     /*  O > S > Se > Te > N > P > As > Sb > Bi > Si > Ge > Sn > Pb > B > Hg                                             */
    // INCHI✔️❌:     /*  Note: Other heteroatoms may be placed within this order as indicated by their positions in the                  */
    // INCHI✔️❌:     /*  periodic table [5].                                                                                             */
    // INCHI✔️❌:
    // INCHI✔️❌:     int erank_rule4[] = { 0,1,198,197,196,202,2,211,215,191,190,189,188,187,206,210,214,183,182,181,180,179,178,177,176,
    // INCHI✔️❌:                           175,174,173,172,171,170,169,205,209,213,165,164,163,162,161,160,159,158,157,156,155,154,153,152,
    // INCHI✔️❌:                           151,204,208,212,147,146,145,144,143,142,141,140,139,138,137,136,135,134,133,132,131,130,129,128,
    // INCHI✔️❌:                           127,126,125,124,123,122,121,201,119,203,207,116,115,114,113,112,111,110,109,108,107,106,105,104,
    // INCHI✔️❌:                           103,102,101,100,99,98,97,96,95,94,93,92,91,90,89,88,87,86,85,84,83,82,81 };
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     int i, j, k, nrings = 0;
    // INCHI✔️❌:     int a1, a2, dummy = 0, bond_type = -1, bond_stereo = -1;
    // INCHI✔️❌:     int *num_ring_sys = NULL, *size_ring_sys = NULL;
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #68112 */
    // INCHI✔️❌:     int err2_len = sizeof(erank_rule2) / sizeof(erank_rule2[0]);
    // INCHI✔️❌:     int err4_len = sizeof(erank_rule4) / sizeof(erank_rule4[0]);
    // INCHI✔️❌:
    // INCHI✔️❌:     if ((NULL == aprops) || !at || !pd) /* djb-rwth: fixing oss-fuzz issue #68329, #68286 */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Establish element ranks for atoms */
    // INCHI✔️❌:     for (k = 0; k < nat; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int atnum = at[k].orig_at_number, index = k;
    // INCHI✔️❌:         U_CHAR err4_ind = at[k].el_number;
    // INCHI✔️❌:         if (cano_nums)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             index = cano_nums[atnum];
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (index >= 0 && err4_ind < err4_len)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             aprops[index].erank = erank_rule4[err4_ind];
    // INCHI✔️❌:             aprops[index].ring_erank = 0;
    // INCHI✔️❌:             aprops[index].ring_size = 0;
    // INCHI✔️❌:             aprops[index].ring_num = -1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* deleted H's go here */
    // INCHI✔️❌:             ;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Establish ring systems assignments for atoms */
    // INCHI✔️❌:     num_ring_sys = (int *) inchi_calloc( (long long)nat + 1, sizeof( int ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     if (NULL == num_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     size_ring_sys = (int *) inchi_calloc( (long long)nat + 1, sizeof( int ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     if (NULL == size_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Note that we get here on the way of InChI2Struct conversion.            */
    // INCHI✔️❌:     /* Break temporarily any of (actually, the first) SRU 'cyclizing' bonds    */
    // INCHI✔️❌:     for (j = 0; j < pd->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pd->units[j]->na > 2 && pd->units[j]->nbkbonds > 0 &&
    // INCHI✔️❌:              pd->units[j]->cyclized == 0 &&
    // INCHI✔️❌:              pd->units[j]->cyclizable == CLOSING_SRU_RING)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             a1 = pd->units[j]->bkbonds[0][0] - 1;
    // INCHI✔️❌:             a2 = pd->units[j]->bkbonds[0][1] - 1;
    // INCHI✔️❌:             OrigAtData_RemoveBond( a1, a2, at, &bond_type, &bond_stereo, &dummy );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     nrings = OAD_Polymer_FindRingSystems( pd, at, nat, num_inp_bonds, num_ring_sys, size_ring_sys, 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Immediately restore just broken bond(s) */
    // INCHI✔️❌:     for (j = 0; j < pd->n; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pd->units[j]->na > 2 &&
    // INCHI✔️❌:              pd->units[j]->nbkbonds > 0 &&
    // INCHI✔️❌:              pd->units[j]->cyclized == 0 &&
    // INCHI✔️❌:              pd->units[j]->cyclizable == CLOSING_SRU_RING)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             a1 = pd->units[j]->bkbonds[0][0] - 1;
    // INCHI✔️❌:             a2 = pd->units[j]->bkbonds[0][1] - 1;
    // INCHI✔️❌:             /* OrigAtData_AddSingleStereolessBond( a1, a2, at, &dummy ); */
    // INCHI✔️❌:             OrigAtData_AddBond( a1, a2, at, bond_type, bond_stereo, &dummy );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nrings)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int max_ring_num = 0;
    // INCHI✔️❌:         /* SRU contains ring[s], proceed with them following (not totally) the IUPAC guidelines */
    // INCHI✔️❌:         for (k = 0; k < nat; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* Browse 0-based original atoms, go to 1-based cano nums domain if cano_nums mapping is suppied */
    // INCHI✔️❌:             int atnum = at[k].orig_at_number, index = k;
    // INCHI✔️❌:             if (cano_nums)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 index = cano_nums[atnum] + 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (num_ring_sys[atnum] >= 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 aprops[index].ring_num = num_ring_sys[atnum];  /* temporarily */
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (max_ring_num < aprops[index].ring_num)
    // INCHI✔️❌:                     max_ring_num = aprops[index].ring_num;          /* NB: OAD_Polymer_FindRingSystems may return num_ring_sys[]  */
    // INCHI✔️❌:                                                                     /* which is not a list of consecutive numbers                 */
    // INCHI✔️❌:
    // INCHI✔️❌:                 aprops[index].ring_size = size_ring_sys[atnum];     /* Size of ring system which includes the atom k .            */
    // INCHI✔️❌:                                                                     /* It is used as an additional score for in-ring              */
    // INCHI✔️❌:                                                                     /* atoms' prioritizing (instead of criteria in                */
    // INCHI✔️❌:                                                                     /* 2c-2h of IUPAC rule 2 which deal with ring sizes).         */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0; i <= max_ring_num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int erank, max_erank = 0;
    // INCHI✔️❌:             for (k = 0; k < nat; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int atnum = at[k].orig_at_number, index = k;
    // INCHI✔️❌:                 if (cano_nums)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     index = cano_nums[atnum] + 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (aprops[index].ring_num == i)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     erank = erank_rule2[at[k].el_number];
    // INCHI✔️❌:                     if (erank > max_erank)
    // INCHI✔️❌:                         max_erank = erank;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             for (k = 0; k < nat; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int atnum = at[k].orig_at_number, index = k;
    // INCHI✔️❌:                 if (cano_nums)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     index = cano_nums[atnum] + 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (aprops[index].ring_num == i)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (aprops[index].ring_size > 2)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         aprops[index].ring_erank = max_erank;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     if (num_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( num_ring_sys );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (size_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( size_ring_sys );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_Polymer_SetAtProps

    const ERANK_RULE2: &[i32] = &[
        0, 1, 198, 197, 196, 202, 2, 216, 215, 191, 190, 189, 188, 187, 206, 210, 214, 183, 182,
        181, 180, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169, 205, 209, 213, 165, 164,
        163, 162, 161, 160, 159, 158, 157, 156, 155, 154, 153, 152, 151, 204, 208, 212, 147, 146,
        145, 144, 143, 142, 141, 140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 130, 129, 128,
        127, 126, 125, 124, 123, 122, 121, 201, 119, 203, 207, 116, 115, 114, 113, 112, 111, 110,
        109, 108, 107, 106, 105, 104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90,
        89, 88, 87, 86, 85, 84, 83, 82, 81,
    ];
    const ERANK_RULE4: &[i32] = &[
        0, 1, 198, 197, 196, 202, 2, 211, 215, 191, 190, 189, 188, 187, 206, 210, 214, 183, 182,
        181, 180, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169, 205, 209, 213, 165, 164,
        163, 162, 161, 160, 159, 158, 157, 156, 155, 154, 153, 152, 151, 204, 208, 212, 147, 146,
        145, 144, 143, 142, 141, 140, 139, 138, 137, 136, 135, 134, 133, 132, 131, 130, 129, 128,
        127, 126, 125, 124, 123, 122, 121, 201, 119, 203, 207, 116, 115, 114, 113, 112, 111, 110,
        109, 108, 107, 106, 105, 104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90,
        89, 88, 87, 86, 85, 84, 83, 82, 81,
    ];

    if atom_properties.is_null() || atoms.is_null() || polymer.is_null() {
        return Ok(());
    }
    let canonical_index = |heap: &SourceHeap, atom_number: i32| -> Result<i32, SourceHeapError> {
        if canonical_numbers.is_null() {
            return Ok(atom_number);
        }
        heap.slice(canonical_numbers.offset(i64::from(atom_number))?)?
            .first()
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    };
    for k in 0..nat {
        let atom = heap
            .slice(atoms.as_const().offset(i64::from(k))?)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let atom_number = i32::from(atom.orig_at_number);
        let index = canonical_index(heap, atom_number)?;
        let element = usize::from(atom.el_number);
        if index >= 0 && element < ERANK_RULE4.len() {
            let properties = heap
                .slice_mut(atom_properties.offset(i64::from(index))?)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            properties.erank = ERANK_RULE4[element];
            properties.ring_erank = 0;
            properties.ring_size = 0;
            properties.ring_num = -1;
        }
    }

    let allocation_count = match usize::try_from(i64::from(nat) + 1) {
        Ok(count) => count,
        Err(_) => return Ok(()),
    };
    let ring_numbers = match heap.allocate(vec![0_i32; allocation_count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(()),
        Err(error) => return Err(error),
    };
    let ring_sizes = match heap.allocate(vec![0_i32; allocation_count]) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => {
            inchi_free(heap, ring_numbers)?;
            return Ok(());
        }
        Err(error) => {
            inchi_free(heap, ring_numbers)?;
            return Err(error);
        }
    };

    let result = (|| -> Result<(), SourceHeapError> {
        let polymer_value = heap
            .slice(polymer)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let mut dummy = 0_i32;
        let mut bond_type = -1_i32;
        let mut bond_stereo = -1_i32;
        for j in 0..polymer_value.n {
            let unit_pointer = *heap
                .slice(polymer_value.units.as_const().offset(i64::from(j))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let unit = heap
                .slice(unit_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if unit.na > 2
                && unit.nbkbonds > 0
                && unit.cyclized == 0
                && unit.cyclizable == CLOSING_SRU_RING as i32
            {
                let row = *heap
                    .slice(unit.bkbonds.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let bond = heap.slice(row.as_const())?;
                let a1 = bond[0].wrapping_sub(1);
                let a2 = bond[1].wrapping_sub(1);
                let _ = OrigAtData_RemoveBond(
                    heap,
                    a1,
                    a2,
                    atoms,
                    &mut bond_type,
                    &mut bond_stereo,
                    &mut dummy,
                )?;
            }
        }

        let rings = OAD_Polymer_FindRingSystems(
            heap,
            polymer,
            atoms,
            nat,
            num_inp_bonds,
            ring_numbers,
            ring_sizes,
            0,
        )?;

        for j in 0..polymer_value.n {
            let unit_pointer = *heap
                .slice(polymer_value.units.as_const().offset(i64::from(j))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let unit = heap
                .slice(unit_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if unit.na > 2
                && unit.nbkbonds > 0
                && unit.cyclized == 0
                && unit.cyclizable == CLOSING_SRU_RING as i32
            {
                let row = *heap
                    .slice(unit.bkbonds.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let bond = heap.slice(row.as_const())?;
                let a1 = bond[0].wrapping_sub(1);
                let a2 = bond[1].wrapping_sub(1);
                let _ =
                    OrigAtData_AddBond(heap, a1, a2, atoms, bond_type, bond_stereo, &mut dummy)?;
            }
        }

        if rings != 0 {
            let mut max_ring_number = 0_i32;
            for k in 0..nat {
                let atom = heap
                    .slice(atoms.as_const().offset(i64::from(k))?)?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                let atom_number = i32::from(atom.orig_at_number);
                let index = if canonical_numbers.is_null() {
                    k
                } else {
                    canonical_index(heap, atom_number)?.wrapping_add(1)
                };
                let ring_number =
                    heap.slice(ring_numbers.as_const())?[usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                if ring_number >= 0 {
                    let ring_size =
                        heap.slice(ring_sizes.as_const())?[usize::try_from(atom_number)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                    let properties = heap
                        .slice_mut(atom_properties.offset(i64::from(index))?)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    properties.ring_num = ring_number;
                    if max_ring_number < properties.ring_num {
                        max_ring_number = properties.ring_num;
                    }
                    properties.ring_size = ring_size;
                }
            }
            for ring_number in 0..=max_ring_number {
                let mut max_erank = 0_i32;
                for k in 0..nat {
                    let atom = heap
                        .slice(atoms.as_const().offset(i64::from(k))?)?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let atom_number = i32::from(atom.orig_at_number);
                    let index = if canonical_numbers.is_null() {
                        k
                    } else {
                        canonical_index(heap, atom_number)?.wrapping_add(1)
                    };
                    let properties = heap
                        .slice(atom_properties.as_const().offset(i64::from(index))?)?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if properties.ring_num == ring_number {
                        max_erank = max_erank.max(ERANK_RULE2[usize::from(atom.el_number)]);
                    }
                }
                for k in 0..nat {
                    let atom = heap
                        .slice(atoms.as_const().offset(i64::from(k))?)?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone();
                    let atom_number = i32::from(atom.orig_at_number);
                    let index = if canonical_numbers.is_null() {
                        k
                    } else {
                        canonical_index(heap, atom_number)?.wrapping_add(1)
                    };
                    let properties = heap
                        .slice_mut(atom_properties.offset(i64::from(index))?)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if properties.ring_num == ring_number && properties.ring_size > 2 {
                        properties.ring_erank = max_erank;
                    }
                }
            }
        }
        Ok(())
    })();

    inchi_free(heap, ring_numbers)?;
    inchi_free(heap, ring_sizes)?;
    result
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_CyclizeCloseableUnits(
    heap: &mut SourceHeap,
    original_atom_data: &mut ORIG_ATOM_DATA,
    _treat_polymers: i32,
    mut error_text: Option<&mut [i8]>,
    no_warnings: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1979 OAD_Polymer_CyclizeCloseableUnits
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int OAD_Polymer_CyclizeCloseableUnits( ORIG_ATOM_DATA *orig_at_data,
                                       int treat_polymers,
                                       char *pStrErr,
                                       int bNoWarnings )
{
    int i, err = 0; /* djb-rwth: removing redundant variables */

    for (i = 0; i < orig_at_data->polymer->n; i++)
    {
        OAD_PolymerUnit *unit = orig_at_data->polymer->units[i];

        if (!unit->cyclizable)
        {
            continue;
        }

        /* Find stars and their partners */
        OAD_PolymerUnit_SetEndsAndCaps( unit, orig_at_data, &err, pStrErr );
            /*	Reveal and store CRU caps and ends('stars and partners')
                Also set `unit->cap1_is_undef`, `unit->cap2_is_undef`, `unit->cyclizable` 
            */
        if (err)
        {
            break;
        }
        if (!unit->cyclizable)
        {
            continue;
        }


        if (OAD_PolymerUnit_HasMetal( unit, orig_at_data->at ))
        {
            /*unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;*/
            if (unit->cyclizable == CLOSING_SRU_RING)
            {
                /*unit->cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;*/
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErr, "Frame shift in metallated polymer unit may be missed" );
                }
            }
        }

        /* Now remove bonds to cap ("star atoms and cyclize a SRU */
        OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms( unit, orig_at_data, &err, pStrErr );

        if (err)
        {
            break;
        }
        if (!unit->cyclizable)
        {
            continue;
        }

        /* djb-rwth: removing redundant code */
    }

    /*
    if ( ncyclized )
    {
        if (!bNoWarnings)
        {
            WarningMessage( pStrErr, "Made provision for frame shift in polymer unit(s)" );
        }
    }
    */

    return err;
}
    */
    // END INCHI C FUNCTION: OAD_Polymer_CyclizeCloseableUnits
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_Polymer_CyclizeCloseableUnits
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; WarningMessage expands to AddErrorMessage.
    // INCHI✔️❌: The final ncyclized warning block is a C comment and inactive; treat_polymers is unused in the active body.
    // INCHI✔️❌: Performance is materially worse because modeled pointer access and unit writeback require SourceHeap lookups and cloning.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_Polymer_CyclizeCloseableUnits

    let polymer = heap
        .slice(original_atom_data.polymer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut error = 0_i32;
    for index in 0..polymer.n {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let unit_pointer = *heap
            .slice(polymer.units.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if unit.cyclizable == 0 {
            continue;
        }

        let set_result = OAD_PolymerUnit_SetEndsAndCaps(
            heap,
            &mut unit,
            original_atom_data,
            &mut error,
            error_text.as_deref_mut(),
        );
        heap.slice_mut(unit_pointer)?[0] = unit.clone();
        set_result?;
        if error != 0 {
            break;
        }
        if unit.cyclizable == 0 {
            continue;
        }

        if OAD_PolymerUnit_HasMetal(heap, &unit, original_atom_data.at.as_const())? != 0
            && unit.cyclizable == CLOSING_SRU_RING as i32
            && no_warnings == 0
        {
            let _ = AddErrorMessage(
                error_text.as_deref_mut(),
                Some(
                    &b"Frame shift in metallated polymer unit may be missed\0"
                        .map(|byte| byte as i8),
                ),
            )?;
        }

        OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(
            heap,
            unit_pointer,
            original_atom_data,
            &mut error,
            error_text.as_deref_mut(),
        )?;
        if error != 0 {
            break;
        }
        if heap
            .slice(unit_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .cyclizable
            == 0
        {
            continue;
        }
    }
    Ok(error)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_HasMetal(
    heap: &SourceHeap,
    unit: &OAD_PolymerUnit,
    atoms: SourceConstPointer<inp_ATOM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2055 OAD_PolymerUnit_HasMetal
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int OAD_PolymerUnit_HasMetal( OAD_PolymerUnit *u, inp_ATOM *at )
{
    int i;
    for (i = 0; i < u->na; i++)
    {
        if (is_el_a_metal( at[u->alist[i] - 1].el_number ))
        {
            return 1;
        }
    }

    return 0;
}
    */
    // END INCHI C FUNCTION: OAD_PolymerUnit_HasMetal
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_PolymerUnit_HasMetal
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no conditional branches.
    // INCHI✔️❌: Performance is materially worse because each modeled pointer access performs a SourceHeap allocation-table lookup.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_PolymerUnit_HasMetal

    if unit.na <= 0 {
        return Ok(0);
    }
    let atom_list = heap.slice(unit.alist.as_const())?;
    let atoms = heap.slice(atoms)?;
    for index in 0..unit.na {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atom_number = *atom_list
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atom_index = atom_number
            .checked_sub(1)
            .and_then(|value| usize::try_from(value).ok())
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let periodic_number = i32::from(
            atoms
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .el_number,
        );
        if is_el_a_metal(periodic_number)? != 0 {
            return Ok(1);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_Free(
    heap: &mut SourceHeap,
    mut pd: SourceMutPointer<OAD_Polymer>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2071 OAD_Polymer_Free
    // INCHI✔❌: void OAD_Polymer_Free( OAD_Polymer *pd )
    // INCHI✔❌: {
    // INCHI✔❌:     if (pd)
    // INCHI✔❌:     {
    // INCHI✔❌:         if (pd->pzz)
    // INCHI✔❌:         {
    // INCHI✔❌:             inchi_free( pd->pzz );
    // INCHI✔❌:             pd->pzz = NULL;
    // INCHI✔❌:             pd->n_pzz = 0;
    // INCHI✔❌:         }
    // INCHI✔❌:         if (pd->n && pd->units)
    // INCHI✔❌:         {
    // INCHI✔❌:             int k;
    // INCHI✔❌:             for (k = 0; k < pd->n; k++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 OAD_PolymerUnit_Free( pd->units[k] );
    // INCHI✔❌:             }
    // INCHI✔❌:             inchi_free( pd->units );
    // INCHI✔❌:             pd->units = NULL;
    // INCHI✔❌:             pd->n = 0;
    // INCHI✔❌:         }
    // INCHI✔❌:         inchi_free( pd );
    // INCHI✔❌:         pd = NULL;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     return;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: OAD_Polymer_Free

    if !pd.is_null() {
        let pzz = heap
            .slice(pd.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .pzz;
        if !pzz.is_null() {
            inchi_free(heap, pzz)?;
            let polymer = heap
                .slice_mut(pd)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            polymer.pzz = SourceMutPointer::null();
            polymer.n_pzz = 0;
        }

        let (n, units) = {
            let polymer = heap
                .slice(pd.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            (polymer.n, polymer.units)
        };
        if n != 0 && !units.is_null() {
            for k in 0..n {
                let unit = heap.slice(units.as_const())?[k as usize];
                OAD_PolymerUnit_Free(heap, unit)?;
            }
            inchi_free(heap, units)?;
            let polymer = heap
                .slice_mut(pd)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            polymer.units = SourceMutPointer::null();
            polymer.n = 0;
        }
        inchi_free(heap, pd)?;
        pd = SourceMutPointer::null();
    }
    let _ = pd;
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(
    heap: &mut SourceHeap,
    unit: SourceMutPointer<OAD_PolymerUnit>,
    original_atom_data: &mut ORIG_ATOM_DATA,
    error: &mut i32,
    _error_text: Option<&mut [i8]>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2101 OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms( OAD_PolymerUnit *unit,
                                                   ORIG_ATOM_DATA  *orig_inp_data,
                                                   int             *err,
                                                   char            *pStrErr )
{
    int bond_type, bond_stereo;

    *err = 0;
    if (!unit->cyclizable)
    {
        return;
    }

    if (unit->cyclizable == CLOSING_SRU_RING)
    {
        /* Disconnect both star atoms */
        OrigAtData_RemoveBond( unit->cap1 - 1, unit->end_atom1 - 1, orig_inp_data->at,
                               &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );

        OrigAtData_RemoveBond( unit->cap2 - 1, unit->end_atom2 - 1, orig_inp_data->at,
                               &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );

        OrigAtData_AddSingleStereolessBond( unit->end_atom1 - 1, unit->end_atom2 - 1,
                                            orig_inp_data->at, &orig_inp_data->num_inp_bonds );
    }

    else if (unit->cyclizable == CLOSING_SRU_HIGHER_ORDER_BOND)
    {
        int elevated; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        elevated = OrigAtData_IncreaseBondOrder( unit->end_atom1 - 1, unit->end_atom2 - 1, orig_inp_data->at ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
#if 0
/* the bond may already be broken at metal disconnection, so ignore the result here */
        if (!elevated)
        {
            /* *err = 1; */
            WarningMessage( pStrErr, "SRU closure via higher order bond failed" );
            unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
            return;
        }
#endif
        OrigAtData_RemoveBond( unit->cap1 - 1, unit->end_atom1 - 1, orig_inp_data->at,
                               &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );
        OrigAtData_RemoveBond( unit->cap2 - 1, unit->end_atom2 - 1, orig_inp_data->at,
                               &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );
    }

    else if (unit->cyclizable == CLOSING_SRU_DIRADICAL)
    {
        orig_inp_data->at[unit->end_atom1 - 1].radical = RADICAL_TRIPLET;
        OrigAtData_RemoveBond( unit->cap1 - 1, unit->end_atom1 - 1, orig_inp_data->at,
                               &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );
        OrigAtData_RemoveBond( unit->cap2 - 1, unit->end_atom2 - 1, orig_inp_data->at,
                                &bond_type, &bond_stereo, &orig_inp_data->num_inp_bonds );
    }

    if (!*err)
    {
        unit->cyclized = 1;
    }
    
    return;
}
    */
    // END INCHI C FUNCTION: OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; the #if 0 elevated-result warning branch is inactive.
    // INCHI✔️❌: Performance is materially worse because each modeled pointer access performs a SourceHeap allocation-table lookup.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms

    *error = 0;
    let initial = heap
        .slice(unit.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if initial.cyclizable == 0 {
        return Ok(());
    }

    let mut bond_type = 0_i32;
    let mut bond_stereo = 0_i32;
    if initial.cyclizable == CLOSING_SRU_RING as i32 {
        let _ = OrigAtData_RemoveBond(
            heap,
            initial.cap1.wrapping_sub(1),
            initial.end_atom1.wrapping_sub(1),
            original_atom_data.at,
            &mut bond_type,
            &mut bond_stereo,
            &mut original_atom_data.num_inp_bonds,
        )?;
        let _ = OrigAtData_RemoveBond(
            heap,
            initial.cap2.wrapping_sub(1),
            initial.end_atom2.wrapping_sub(1),
            original_atom_data.at,
            &mut bond_type,
            &mut bond_stereo,
            &mut original_atom_data.num_inp_bonds,
        )?;
        let _ = OrigAtData_AddSingleStereolessBond(
            heap,
            initial.end_atom1.wrapping_sub(1),
            initial.end_atom2.wrapping_sub(1),
            original_atom_data.at,
            &mut original_atom_data.num_inp_bonds,
        )?;
    } else if initial.cyclizable == CLOSING_SRU_HIGHER_ORDER_BOND as i32 {
        let _elevated = OrigAtData_IncreaseBondOrder(
            heap,
            initial.end_atom1.wrapping_sub(1),
            initial.end_atom2.wrapping_sub(1),
            original_atom_data.at,
        )?;
        let _ = OrigAtData_RemoveBond(
            heap,
            initial.cap1.wrapping_sub(1),
            initial.end_atom1.wrapping_sub(1),
            original_atom_data.at,
            &mut bond_type,
            &mut bond_stereo,
            &mut original_atom_data.num_inp_bonds,
        )?;
        let _ = OrigAtData_RemoveBond(
            heap,
            initial.cap2.wrapping_sub(1),
            initial.end_atom2.wrapping_sub(1),
            original_atom_data.at,
            &mut bond_type,
            &mut bond_stereo,
            &mut original_atom_data.num_inp_bonds,
        )?;
    } else if initial.cyclizable == CLOSING_SRU_DIRADICAL as i32 {
        let end_index = initial
            .end_atom1
            .checked_sub(1)
            .and_then(|value| usize::try_from(value).ok())
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        heap.slice_mut(original_atom_data.at)?
            .get_mut(end_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .radical = RADICAL_TRIPLET as i8;
        let _ = OrigAtData_RemoveBond(
            heap,
            initial.cap1.wrapping_sub(1),
            initial.end_atom1.wrapping_sub(1),
            original_atom_data.at,
            &mut bond_type,
            &mut bond_stereo,
            &mut original_atom_data.num_inp_bonds,
        )?;
        let _ = OrigAtData_RemoveBond(
            heap,
            initial.cap2.wrapping_sub(1),
            initial.end_atom2.wrapping_sub(1),
            original_atom_data.at,
            &mut bond_type,
            &mut bond_stereo,
            &mut original_atom_data.num_inp_bonds,
        )?;
    }

    if *error == 0 {
        heap.slice_mut(unit)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .cyclized = 1;
    }
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OAD_PolymerUnit_FindEndsAndCaps(
    heap: &SourceHeap,
    unit: &mut OAD_PolymerUnit,
    original_atom_data: &ORIG_ATOM_DATA,
    end1: &mut i32,
    cap1: &mut i32,
    cap1_is_star: &mut i32,
    end2: &mut i32,
    cap2: &mut i32,
    cap2_is_star: &mut i32,
    error: &mut i32,
    mut error_text: Option<&mut [i8]>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2166 OAD_PolymerUnit_FindEndsAndCaps
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void OAD_PolymerUnit_FindEndsAndCaps( OAD_PolymerUnit *unit,
                                      ORIG_ATOM_DATA *orig_at_data,
                                      int *end1,
                                      int *cap1,
                                      int *cap1_is_star,
                                      int *end2,
                                      int *cap2,
                                      int *cap2_is_star,
                                      int *err,
                                      char *pStrErr )
{
    int i, j, i_inside = 0, j_inside = 0;
    int num_atoms = orig_at_data->num_inp_atoms;

    *end1 = *end2 = *cap1 = *cap2 = 0;
    *cap1_is_star = *cap2_is_star = 0;
    *err = 0;

    if (!unit->blist || unit->nb < 1)
    {
        return;
    }
    /* Left crossing bond */
    i = unit->blist[0];
    j = unit->blist[1];
    i_inside = (NULL != is_in_the_ilist(unit->alist, i, unit->na));
    j_inside = (NULL != is_in_the_ilist(unit->alist, j, unit->na));
    if (i_inside && j_inside)
    {
        TREAT_ERR(*err, 9032, "Polymer CRU cap(s) lie inside CRU");
        return;
    }
    if (i_inside)
    {
        *end1 = i;
        *cap1 = j;
    }
    else
    {
        *end1 = j;
        *cap1 = i;
    }
    if (!strcmp(orig_at_data->at[*cap1 - 1].elname, "Zz"))
    {
        *cap1_is_star = 1;
    }
    /* Right crossing bond */
    i = unit->blist[2];
    j = unit->blist[3];
    i_inside = NULL != is_in_the_ilist(unit->alist, i, unit->na);
    j_inside = NULL != is_in_the_ilist(unit->alist, j, unit->na);
    if (i_inside && j_inside)
    {
        TREAT_ERR(*err, 9032, "Polymer CRU cap(s) lie inside CRU");
    }
    if (i_inside)
    {
        *end2 = i;
        *cap2 = j;
    }
    else
    {
        *end2 = j;
        *cap2 = i;
    }
    if (!strcmp(orig_at_data->at[*cap2 - 1].elname, "Zz"))
    {
        *cap2_is_star = 1;
    }
    /* Checks */
    if (*end1 <= 0 || *end1 > num_atoms || *cap1 <= 0 || *cap1 > num_atoms)
    {
        TREAT_ERR(*err, 9090, "Invalid polymer CRU crossing bond");
        return;
    }
    if (*end2 <= 0 || *end2 > num_atoms || *cap2 <= 0 || *cap2 > num_atoms)
    {
        TREAT_ERR(*err, 9091, "Invalid polymer CRU crossing bond");
        return;
    }
    if  ( *cap1 == *cap2 ) /* || (*end1 == *end2 && unit->na>1)  */
    {
        TREAT_ERR(*err, 9090, "Invalid polymer CRU surrounding");
        return;
    }

    /* Paranoia, 2020-05-22 */
    unit->end_atom1 = *end1;
    unit->end_atom2 = *end2;
    unit->cap1 = *cap1;
    unit->cap2 = *cap2;

    *err = 0;
    return;
}
    */
    // END INCHI C FUNCTION: OAD_PolymerUnit_FindEndsAndCaps
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_PolymerUnit_FindEndsAndCaps
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; TREAT_ERR preserves an
    // INCHI✔️✔️: existing error and always invokes AddErrorMessage; strcmp and list lookup are active.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_PolymerUnit_FindEndsAndCaps

    fn treat_error(
        error: &mut i32,
        code: i32,
        error_text: &mut Option<&mut [i8]>,
        message: &[i8],
    ) -> Result<(), SourceHeapError> {
        if *error == 0 && code != 0 {
            *error = code;
        }
        AddErrorMessage(error_text.as_deref_mut(), Some(message))?;
        Ok(())
    }

    fn is_star(atom: &inp_ATOM) -> Result<bool, SourceHeapError> {
        let nul = atom
            .elname
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
        Ok(atom.elname[..nul] == [b'Z' as i8, b'z' as i8])
    }

    *end1 = 0;
    *end2 = 0;
    *cap1 = 0;
    *cap2 = 0;
    *cap1_is_star = 0;
    *cap2_is_star = 0;
    *error = 0;

    if unit.blist.is_null() || unit.nb < 1 {
        return Ok(());
    }
    let crossing_bonds = heap.slice(unit.blist.as_const())?;
    let left_i = *crossing_bonds
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let left_j = *crossing_bonds
        .get(1)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let atom_list = if unit.na == 0 && unit.alist.is_null() {
        None
    } else {
        Some(heap.slice(unit.alist.as_const())?)
    };
    let left_i_inside = is_in_the_ilist(atom_list, left_i, unit.na)?.is_some();
    let left_j_inside = is_in_the_ilist(atom_list, left_j, unit.na)?.is_some();
    if left_i_inside && left_j_inside {
        treat_error(
            error,
            9032,
            &mut error_text,
            &b"Polymer CRU cap(s) lie inside CRU\0".map(|byte| byte as i8),
        )?;
        return Ok(());
    }
    if left_i_inside {
        *end1 = left_i;
        *cap1 = left_j;
    } else {
        *end1 = left_j;
        *cap1 = left_i;
    }
    let cap1_index = cap1
        .checked_sub(1)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
    let atoms = heap.slice(original_atom_data.at.as_const())?;
    let first_cap_atom = atoms
        .get(cap1_index)
        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
    if is_star(first_cap_atom)? {
        *cap1_is_star = 1;
    }

    let right_i = *crossing_bonds
        .get(2)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let right_j = *crossing_bonds
        .get(3)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let right_i_inside = is_in_the_ilist(atom_list, right_i, unit.na)?.is_some();
    let right_j_inside = is_in_the_ilist(atom_list, right_j, unit.na)?.is_some();
    if right_i_inside && right_j_inside {
        treat_error(
            error,
            9032,
            &mut error_text,
            &b"Polymer CRU cap(s) lie inside CRU\0".map(|byte| byte as i8),
        )?;
    }
    if right_i_inside {
        *end2 = right_i;
        *cap2 = right_j;
    } else {
        *end2 = right_j;
        *cap2 = right_i;
    }
    let cap2_index = cap2
        .checked_sub(1)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
    let second_cap_atom = atoms
        .get(cap2_index)
        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
    if is_star(second_cap_atom)? {
        *cap2_is_star = 1;
    }

    let number_of_atoms = original_atom_data.num_inp_atoms;
    if *end1 <= 0 || *end1 > number_of_atoms || *cap1 <= 0 || *cap1 > number_of_atoms {
        treat_error(
            error,
            9090,
            &mut error_text,
            &b"Invalid polymer CRU crossing bond\0".map(|byte| byte as i8),
        )?;
        return Ok(());
    }
    if *end2 <= 0 || *end2 > number_of_atoms || *cap2 <= 0 || *cap2 > number_of_atoms {
        treat_error(
            error,
            9091,
            &mut error_text,
            &b"Invalid polymer CRU crossing bond\0".map(|byte| byte as i8),
        )?;
        return Ok(());
    }
    if *cap1 == *cap2 {
        treat_error(
            error,
            9090,
            &mut error_text,
            &b"Invalid polymer CRU surrounding\0".map(|byte| byte as i8),
        )?;
        return Ok(());
    }

    unit.end_atom1 = *end1;
    unit.end_atom2 = *end2;
    unit.cap1 = *cap1;
    unit.cap2 = *cap2;
    *error = 0;
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_SetEndsAndCaps(
    heap: &SourceHeap,
    unit: &mut OAD_PolymerUnit,
    original_atom_data: &ORIG_ATOM_DATA,
    error: &mut i32,
    error_text: Option<&mut [i8]>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2266 OAD_PolymerUnit_SetEndsAndCaps
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void OAD_PolymerUnit_SetEndsAndCaps( OAD_PolymerUnit *unit,
                                     ORIG_ATOM_DATA  *orig_at_data,
                                     int             *err,
                                     char            *pStrErr )
{
    int k;

    unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
    unit->end_atom1 = unit->end_atom2 = unit->cap1 = unit->cap2 = -1;
    unit->cap1_is_undef = unit->cap2_is_undef = 0;

    OAD_PolymerUnit_FindEndsAndCaps( unit, orig_at_data,
                                     &unit->end_atom1, &unit->cap1, &unit->cap1_is_undef,
                                     &unit->end_atom2, &unit->cap2, &unit->cap2_is_undef,
                                     err, pStrErr);

    if (*err)
    {
        goto exit_function;
    }

#if ( defined(DEBUG_POLYMERS) && ( DEBUG_POLYMERS != 0 ) )
    ITRACE_( "Cap-end_atom pairs (numbers are from 1) are: %-d-%-d and %-d-%-d\n",
        unit->cap1, unit->end_atom1, unit->cap2, unit->end_atom2 );
#endif

    if (!unit->cap1_is_undef && !unit->cap2_is_undef)
    {
        goto exit_function;
    }

    /* The rest is applicable only to *---SRU---* case */

    /* Stars are separated by one atom - that's not error but do nothing */
    if (unit->end_atom1 == unit->end_atom2)
    {
#ifdef ALLOW_CLOSING_SRU_VIA_DIRADICAL
        unit->cyclizable = CLOSING_SRU_DIRADICAL;
#else
        unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
#endif
        goto exit_function;
    }

    /* Stars are separated by two atoms - that's not error but do nothing */
    for (k = 0; k < orig_at_data->at[unit->end_atom1 - 1].valence; k++)
    {
        if (orig_at_data->at[unit->end_atom1 - 1].neighbor[k] == unit->end_atom2 - 1)
        {
#ifdef ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND
            unit->cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;
#else
            unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
#endif
            goto exit_function;
        }
    }

    unit->cyclizable = CLOSING_SRU_RING;

exit_function:

    return;
}
    */
    // END INCHI C FUNCTION: OAD_PolymerUnit_SetEndsAndCaps
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_PolymerUnit_SetEndsAndCaps
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64;
    // INCHI✔️✔️: ALLOW_CLOSING_SRU_VIA_DIRADICAL=1 and ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND=1;
    // INCHI✔️✔️: DEBUG_POLYMERS is undefined, so its ITRACE_ branch is inactive.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_PolymerUnit_SetEndsAndCaps

    unit.cyclizable = CLOSING_SRU_NOT_APPLICABLE as i32;
    unit.end_atom1 = -1;
    unit.end_atom2 = -1;
    unit.cap1 = -1;
    unit.cap2 = -1;
    unit.cap1_is_undef = 0;
    unit.cap2_is_undef = 0;

    let mut end1 = unit.end_atom1;
    let mut cap1 = unit.cap1;
    let mut cap1_is_star = unit.cap1_is_undef;
    let mut end2 = unit.end_atom2;
    let mut cap2 = unit.cap2;
    let mut cap2_is_star = unit.cap2_is_undef;
    let find_result = OAD_PolymerUnit_FindEndsAndCaps(
        heap,
        unit,
        original_atom_data,
        &mut end1,
        &mut cap1,
        &mut cap1_is_star,
        &mut end2,
        &mut cap2,
        &mut cap2_is_star,
        error,
        error_text,
    );
    unit.end_atom1 = end1;
    unit.cap1 = cap1;
    unit.cap1_is_undef = cap1_is_star;
    unit.end_atom2 = end2;
    unit.cap2 = cap2;
    unit.cap2_is_undef = cap2_is_star;
    find_result?;

    if *error != 0 {
        return Ok(());
    }
    if unit.cap1_is_undef == 0 && unit.cap2_is_undef == 0 {
        return Ok(());
    }
    if unit.end_atom1 == unit.end_atom2 {
        unit.cyclizable = CLOSING_SRU_DIRADICAL as i32;
        return Ok(());
    }

    let first_end_index = unit
        .end_atom1
        .checked_sub(1)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
    let atoms = heap.slice(original_atom_data.at.as_const())?;
    let first_end = atoms
        .get(first_end_index)
        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
    let valence = usize::try_from(first_end.valence)
        .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
    if valence > first_end.neighbor.len() {
        return Err(SourceHeapError::UnsupportedSourceBehavior);
    }
    let second_end_index = unit
        .end_atom2
        .checked_sub(1)
        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
    for neighbor in &first_end.neighbor[..valence] {
        if i32::from(*neighbor) == second_end_index {
            unit.cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND as i32;
            return Ok(());
        }
    }
    unit.cyclizable = CLOSING_SRU_RING as i32;
    Ok(())
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OAD_CollectReachableAtoms(
    heap: &mut SourceHeap,
    original_atom_data: &ORIG_ATOM_DATA,
    start_atom: i32,
    forbidden_bond_count: i32,
    forbidden_bond_atoms: SourceMutPointer<i32>,
    reachable_count: &mut i32,
    reachable: SourceMutPointer<i32>,
    _error: &mut i32,
    _error_text: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3019 OAD_CollectReachableAtoms
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int OAD_CollectReachableAtoms( ORIG_ATOM_DATA  *orig_at_data,
                               int start_atom,
                               int nforbidden_bonds,
                               int *forbidden_bond_atoms,
                               int *n_reachable,
                               int *reachable,
                               int *err,
                               char *pStrErr)
{
    int iatom, natnums, max_atoms, j, ret = _IS_OKAY;
    int max_n_reachable = *n_reachable;
    int *atnums = NULL;
    subgraf *sg = NULL;
    subgraf_pathfinder *spf = NULL;

    /* djb-rwth: removing redundant code */
    max_atoms = orig_at_data->num_inp_atoms;
    iatom = start_atom - 1;
    *n_reachable = 0;

    atnums = (int *)inchi_calloc(max_atoms, sizeof(int));
    if (!atnums)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    for (j = 0; j < max_atoms; j++)
    {
        atnums[j] = orig_at_data->at[j].orig_at_number;
    }
    sg = subgraf_new(orig_at_data, max_atoms, atnums);
    if (!sg)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    spf = subgraf_pathfinder_new(sg, orig_at_data, iatom, iatom);
    if (!spf)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    /* move from orig# to node# */
    spf->start = iatom;
    for (j = 0; j < nforbidden_bonds; j++)
    {
        forbidden_bond_atoms[2 * j]      = sg->orig2node[forbidden_bond_atoms[2 * j] ];
        forbidden_bond_atoms[2 * j + 1]  = sg->orig2node[forbidden_bond_atoms[2 * j + 1] ];
    }

    /*memset(atnums, -1, max_atoms * sizeof(int));*/
    for (j = 0; j < max_atoms; j++)
    {
        atnums[j] = -1;
    }

    spf->nseen = 0;
    natnums = subgraf_pathfinder_collect_all(spf, nforbidden_bonds, forbidden_bond_atoms, atnums);

    if (natnums)
    {
        if (natnums > max_n_reachable)
        {
            ret = _IS_ERROR;
            goto exit_function;
        }

        for (j = 0; j < natnums && j < max_atoms; j++) /* djb-rwth: fixing buffer overruns */
        {
            reachable[(*n_reachable)++ ] = atnums[j];
        }
    }

exit_function:
    subgraf_free(sg);
    subgraf_pathfinder_free(spf);
    if (atnums)
    {
        inchi_free(atnums);
    }
    
    return ret;
}
    */
    // END INCHI C FUNCTION: OAD_CollectReachableAtoms
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_CollectReachableAtoms
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; sizeof(int)=4.
    // INCHI✔️❌: inchi_calloc and inchi_free expand to calloc and free; _IS_OKAY=0 and _IS_ERROR=2.
    // INCHI✔️❌: err and pStrErr are unused; forbidden endpoints are mutated in place only after all graph/pathfinder allocations succeed.
    // INCHI✔️❌: Rust preserves DFS order and cleanup order, but SourceHeap lookups and the node-list clone are materially more expensive than native pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_CollectReachableAtoms

    fn is_allocation_error(error: SourceHeapError) -> bool {
        matches!(
            error,
            SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed
        )
    }

    let max_reachable = *reachable_count;
    let max_atoms = original_atom_data.num_inp_atoms;
    let start = start_atom
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    *reachable_count = 0;

    let mut atom_numbers = SourceMutPointer::null();
    let mut graph = SourceMutPointer::null();
    let mut pathfinder = SourceMutPointer::null();

    let operation = (|| -> Result<i32, SourceHeapError> {
        atom_numbers = match inchi_calloc::<i32>(heap, max_atoms as u64, 4) {
            Ok(pointer) => pointer,
            Err(error) if is_allocation_error(error) => return Ok(_IS_ERROR as i32),
            Err(error) => return Err(error),
        };

        let atom_count = usize::try_from(max_atoms.max(0))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let original_numbers = heap
            .slice(original_atom_data.at.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .iter()
            .map(|atom| i32::from(atom.orig_at_number))
            .collect::<Vec<_>>();
        for index in 0..atom_count {
            heap.slice_mut(atom_numbers)?[index] = original_numbers[index];
        }

        let nodes = heap.slice(atom_numbers.as_const())?.to_vec();
        graph = subgraf_new(heap, original_atom_data, max_atoms, &nodes)?;
        if graph.is_null() {
            return Ok(_IS_ERROR as i32);
        }
        pathfinder = subgraf_pathfinder_new(heap, graph, original_atom_data, start, start)?;
        if pathfinder.is_null() {
            return Ok(_IS_ERROR as i32);
        }

        heap.slice_mut(pathfinder)?[0].start = start;
        let map_pointer = heap.slice(graph.as_const())?[0].orig2node;
        for forbidden_index in 0..forbidden_bond_count.max(0) {
            let offset = usize::try_from(forbidden_index)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                .checked_mul(2)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let first_original = *heap
                .slice(forbidden_bond_atoms.as_const())?
                .get(offset)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let first_node = *heap
                .slice(map_pointer.as_const())?
                .get(usize::try_from(first_original).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            *heap
                .slice_mut(forbidden_bond_atoms)?
                .get_mut(offset)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = first_node;

            let second_original = *heap
                .slice(forbidden_bond_atoms.as_const())?
                .get(offset + 1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let second_node = *heap
                .slice(map_pointer.as_const())?
                .get(usize::try_from(second_original).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            *heap
                .slice_mut(forbidden_bond_atoms)?
                .get_mut(offset + 1)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = second_node;
        }

        heap.slice_mut(atom_numbers)?.fill(-1);
        heap.slice_mut(pathfinder)?[0].nseen = 0;
        let collected = subgraf_pathfinder_collect_all(
            heap,
            pathfinder,
            forbidden_bond_count,
            forbidden_bond_atoms,
            atom_numbers,
        )?;

        if collected != 0 {
            if collected > max_reachable {
                return Ok(_IS_ERROR as i32);
            }
            let copy_count = usize::try_from(collected.max(0))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                .min(atom_count);
            for index in 0..copy_count {
                let atom_number = heap.slice(atom_numbers.as_const())?[index];
                let output_index = usize::try_from(*reachable_count)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                *heap
                    .slice_mut(reachable)?
                    .get_mut(output_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = atom_number;
                *reachable_count = reachable_count
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
        }
        Ok(_IS_OKAY as i32)
    })();

    let graph_cleanup = subgraf_free(heap, graph);
    let pathfinder_cleanup = subgraf_pathfinder_free(heap, pathfinder);
    let atom_numbers_cleanup = inchi_free(heap, atom_numbers);
    match operation {
        Err(error) => Err(error),
        Ok(status) => {
            graph_cleanup?;
            pathfinder_cleanup?;
            atom_numbers_cleanup?;
            Ok(status)
        }
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OAD_CollectBackboneBonds(
    heap: &mut SourceHeap,
    atom_data: &ORIG_ATOM_DATA,
    atom_count: i32,
    atom_list: &[i32],
    end_atom1: i32,
    end_atom2: i32,
    backbone_bond_count: &mut i32,
    backbone_bonds: SourceMutPointer<SourceMutPointer<i32>>,
    error: &mut i32,
    error_text: Option<&mut [i8]>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3109 OAD_CollectBackboneBonds
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: OAD_CollectBackboneBonds
    // INCHI✔️❌: void OAD_CollectBackboneBonds(ORIG_ATOM_DATA  *at_data,
    // INCHI✔️❌:                               int na,
    // INCHI✔️❌:                               int *alist,
    // INCHI✔️❌:                               int end_atom1,
    // INCHI✔️❌:                               int end_atom2,
    // INCHI✔️❌:                               int *nbkbonds,
    // INCHI✔️❌:                               int **bkbonds,
    // INCHI✔️❌:                               int *err,
    // INCHI✔️❌:                               char *pStrErr )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int start = 0, end = 0, dummy;
    // INCHI✔️❌:     subgraf *sg = NULL;
    // INCHI✔️❌:     subgraf_pathfinder *spf = NULL;
    // INCHI✔️❌:     /* Establish subgraph for na atoms of the alist */
    // INCHI✔️❌:     *nbkbonds = 0;
    // INCHI✔️❌:     sg = subgraf_new( at_data, na, alist );
    // INCHI✔️❌:     if (!sg)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         TREAT_ERR( *err, 9037, "Not enough memory (polymers)" );
    // INCHI✔️❌:         /* unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE; */
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     start = sg->orig2node[end_atom1];
    // INCHI✔️❌:     end = sg->orig2node[end_atom2];
    // INCHI✔️❌: #if 0
    // INCHI✔️❌:     if (start > end)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int tmp = end;
    // INCHI✔️❌:         end = start;
    // INCHI✔️❌:         start = tmp;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     spf = subgraf_pathfinder_new( sg, at_data, start, end );
    // INCHI✔️❌:     if (!spf)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         TREAT_ERR( *err, 9039, "Not enough memory (polymers)" ); /* djb-rwth: addressing coverity ID #499539 -- TREAT_ERR properly used */
    // INCHI✔️❌:         /*unit->cyclizable = CLOSING_SRU_NOT_APPLICABLE;*/
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     spf->seen[0] = spf->start;
    // INCHI✔️❌:     spf->nseen = 1;
    // INCHI✔️❌:     *nbkbonds = 0;
    // INCHI✔️❌:     subgraf_pathfinder_run( spf, 0, NULL,
    // INCHI✔️❌:                             nbkbonds,
    // INCHI✔️❌:                             bkbonds, /* we will collect backbone CRU bonds here */
    // INCHI✔️❌:                             &dummy,
    // INCHI✔️❌:                             NULL );
    // INCHI✔️❌:
    // INCHI✔️❌:     subgraf_free( sg );
    // INCHI✔️❌:     subgraf_pathfinder_free( spf );
    // INCHI✔️❌:     *err = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: OAD_CollectBackboneBonds
    // END INCHI C FUNCTION: OAD_CollectBackboneBonds
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_CollectBackboneBonds
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; #if 0 endpoint sorting is inactive; TREAT_ERR retains a pre-existing nonzero error and always calls AddErrorMessage.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_CollectBackboneBonds

    fn treat_error(
        error: &mut i32,
        code: i32,
        error_text: Option<&mut [i8]>,
    ) -> Result<(), SourceHeapError> {
        if *error == 0 && code != 0 {
            *error = code;
        }
        AddErrorMessage(
            error_text,
            Some(
                b"Not enough memory (polymers)\0"
                    .map(|byte| byte as i8)
                    .as_slice(),
            ),
        )?;
        Ok(())
    }

    *backbone_bond_count = 0;
    let graph = subgraf_new(heap, atom_data, atom_count, atom_list)?;
    if graph.is_null() {
        treat_error(error, 9037, error_text)?;
        return Ok(());
    }
    let graph_value = heap
        .slice(graph.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let start_index =
        usize::try_from(end_atom1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let end_index = usize::try_from(end_atom2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let mapping = heap.slice(graph_value.orig2node.as_const())?;
    let start = *mapping
        .get(start_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let end = *mapping
        .get(end_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let pathfinder = subgraf_pathfinder_new(heap, graph, atom_data, start, end)?;
    if pathfinder.is_null() {
        treat_error(error, 9039, error_text)?;
        return Ok(());
    }
    let seen = heap
        .slice(pathfinder.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .seen;
    *heap
        .slice_mut(seen)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = start;
    heap.slice_mut(pathfinder)?[0].nseen = 1;
    *backbone_bond_count = 0;
    let mut dummy = 0_i32;
    subgraf_pathfinder_run(
        heap,
        pathfinder,
        0,
        SourceMutPointer::null(),
        backbone_bond_count,
        backbone_bonds,
        &mut dummy,
        SourceMutPointer::null(),
    )?;
    subgraf_free(heap, graph)?;
    subgraf_pathfinder_free(heap, pathfinder)?;
    *error = 0;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_RemoveLinkFromCRUChain(
    heap: &mut SourceHeap,
    atom1: i32,
    atom2: i32,
    bond_count: &mut i32,
    bonds: SourceMutPointer<SourceMutPointer<i32>>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3607 OAD_PolymerUnit_RemoveLinkFromCRUChain
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: OAD_PolymerUnit_RemoveLinkFromCRUChain
    // INCHI✔️✔️: void OAD_PolymerUnit_RemoveLinkFromCRUChain( int at1,
    // INCHI✔️✔️:                                              int at2,
    // INCHI✔️✔️:                                              int *nbonds,
    // INCHI✔️✔️:                                              int **bonds )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int p, q;
    // INCHI✔️✔️: #if 0
    // INCHI✔️✔️:     if (at1 > at2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         int tmp = at1;
    // INCHI✔️✔️:         at1 = at2;
    // INCHI✔️✔️:         at2 = tmp;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     for (p = 0; p < *nbonds; p++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (bonds[p][0] == at1 && bonds[p][1] == at2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             for (q = p + 1; q < *nbonds; q++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 bonds[q - 1][0] = bonds[q][0];
    // INCHI✔️✔️:                 bonds[q - 1][1] = bonds[q][1];
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             ( *nbonds )--;
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return;
    // INCHI✔️✔️: }
    // END COMPLETE VERBATIM SOURCE FRAME: OAD_PolymerUnit_RemoveLinkFromCRUChain
    // END INCHI C FUNCTION: OAD_PolymerUnit_RemoveLinkFromCRUChain
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_PolymerUnit_RemoveLinkFromCRUChain
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; the #if 0 endpoint-sorting block is inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_PolymerUnit_RemoveLinkFromCRUChain

    let count =
        usize::try_from((*bond_count).max(0)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    for position in 0..count {
        let row = *heap
            .slice(bonds.as_const())?
            .get(position)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let values = heap.slice(row.as_const())?;
        let first = *values.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        let second = *values.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if first == atom1 && second == atom2 {
            for source_position in position + 1..count {
                let source_row = *heap
                    .slice(bonds.as_const())?
                    .get(source_position)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let source_values = heap.slice(source_row.as_const())?;
                let source_first = *source_values
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let source_second = *source_values
                    .get(1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let target_row = *heap
                    .slice(bonds.as_const())?
                    .get(source_position - 1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                *heap
                    .slice_mut(target_row)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = source_first;
                *heap
                    .slice_mut(target_row)?
                    .get_mut(1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = source_second;
            }
            *bond_count = bond_count.wrapping_sub(1);
            break;
        }
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_DelistIntraRingBackboneBonds(
    heap: &mut SourceHeap,
    unit: SourceMutPointer<OAD_PolymerUnit>,
    atom_data: &mut ORIG_ATOM_DATA,
    error: &mut i32,
    _error_text: Option<&mut [i8]>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3168 OAD_PolymerUnit_DelistIntraRingBackboneBonds
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: OAD_PolymerUnit_DelistIntraRingBackboneBonds
    // INCHI✔️❌: void OAD_PolymerUnit_DelistIntraRingBackboneBonds( OAD_PolymerUnit *unit,
    // INCHI✔️❌:                                                    ORIG_ATOM_DATA  *at_data,
    // INCHI✔️❌:                                                    int             *err,
    // INCHI✔️❌:                                                    char            *pStrErr )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nrings = 0;
    // INCHI✔️❌:     int *num_ring_sys = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!unit)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (unit->nbkbonds < 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Establish ring systems assignments for all related atoms */
    // INCHI✔️❌:
    // INCHI✔️❌:     *err = 1;
    // INCHI✔️❌:     num_ring_sys = (int *) inchi_calloc( (long long)at_data->num_inp_atoms + 1, sizeof( int ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     if (!num_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     *err = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     nrings = OAD_Polymer_FindRingSystems( at_data->polymer, at_data->at, at_data->num_inp_atoms, &at_data->num_inp_bonds,
    // INCHI✔️❌:                                           num_ring_sys, NULL, unit->end_atom1 - 1 ); /* NB: start dfs within connected compt! */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (nrings == 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int at1, at2, j = 0;
    // INCHI✔️❌: repeatj:
    // INCHI✔️❌:         at1 = unit->bkbonds[j][0];
    // INCHI✔️❌:         at2 = unit->bkbonds[j][1];
    // INCHI✔️❌:         if (( num_ring_sys[at1] == num_ring_sys[at2] ) && ( num_ring_sys[at1] != -1 ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             OAD_PolymerUnit_RemoveLinkFromCRUChain( at1, at2, &unit->nbkbonds, unit->bkbonds );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ++j;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (j < unit->nbkbonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto repeatj;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     if (num_ring_sys)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( num_ring_sys );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: OAD_PolymerUnit_DelistIntraRingBackboneBonds
    // END INCHI C FUNCTION: OAD_PolymerUnit_DelistIntraRingBackboneBonds
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_PolymerUnit_DelistIntraRingBackboneBonds
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; inchi_calloc is calloc and inchi_free is free; pStrErr is inactive data in this function.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_PolymerUnit_DelistIntraRingBackboneBonds

    if unit.is_null() {
        return Ok(());
    }
    let initial = heap
        .slice(unit.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if initial.nbkbonds < 1 {
        return Ok(());
    }
    *error = 1;
    let allocation_count = i64::from(atom_data.num_inp_atoms).wrapping_add(1) as u64;
    let ring_systems = match inchi_calloc::<i32>(heap, allocation_count, 4) {
        Ok(pointer) => pointer,
        Err(
            SourceHeapError::AllocationSizeOverflow
            | SourceHeapError::AllocationElementCountOutOfRange
            | SourceHeapError::AllocationFailed,
        ) => return Ok(()),
        Err(other) => return Err(other),
    };
    *error = 0;
    let rings = OAD_Polymer_FindRingSystems(
        heap,
        atom_data.polymer.as_const(),
        atom_data.at,
        atom_data.num_inp_atoms,
        &mut atom_data.num_inp_bonds,
        ring_systems,
        SourceMutPointer::null(),
        initial.end_atom1.wrapping_sub(1),
    )?;
    if rings != 0 {
        let mut position = 0_i32;
        loop {
            let current = heap
                .slice(unit.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if position >= current.nbkbonds {
                break;
            }
            let position_index =
                usize::try_from(position).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let row = *heap
                .slice(current.bkbonds.as_const())?
                .get(position_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let values = heap.slice(row.as_const())?;
            let atom1 = *values.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom2 = *values.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom1_index =
                usize::try_from(atom1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom2_index =
                usize::try_from(atom2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let systems = heap.slice(ring_systems.as_const())?;
            let first_system = *systems
                .get(atom1_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let second_system = *systems
                .get(atom2_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if first_system == second_system && first_system != -1 {
                let mut count = current.nbkbonds;
                OAD_PolymerUnit_RemoveLinkFromCRUChain(
                    heap,
                    atom1,
                    atom2,
                    &mut count,
                    current.bkbonds,
                )?;
                heap.slice_mut(unit)?[0].nbkbonds = count;
            } else {
                position = position.wrapping_add(1);
            }
        }
    }
    inchi_free(heap, ring_systems)?;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_DelistHighOrderBackboneBonds(
    heap: &mut SourceHeap,
    unit: SourceMutPointer<OAD_PolymerUnit>,
    original_atom_data: &ORIG_ATOM_DATA,
    composite_normalized_data: Option<&COMP_ATOM_DATA>,
    _error: &mut i32,
    _error_text: Option<&mut [i8]>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3509 OAD_PolymerUnit_DelistHighOrderBackboneBonds
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: OAD_PolymerUnit_DelistHighOrderBackboneBonds
    // INCHI✔️❌: void OAD_PolymerUnit_DelistHighOrderBackboneBonds( OAD_PolymerUnit *unit,
    // INCHI✔️❌:                                                    ORIG_ATOM_DATA  *orig_at_data,
    // INCHI✔️❌:                                                    COMP_ATOM_DATA  *composite_norm_data,
    // INCHI✔️❌:                                                    int             *err,
    // INCHI✔️❌:                                                    char            *pStrErr )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int at1, at2, j = 0, k, check_taut = 0, remove; /* djb-rwth: removing redundant variables/code */
    // INCHI✔️❌:
    // INCHI✔️❌:     int *orig_num = NULL, *curr_num = NULL;
    // INCHI✔️❌:     int bond_is_untouchable = 0, btype;  /* DT: moved from below 2024-09-01 DT */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (unit->na < 2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (unit->nb < 2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (unit->nbkbonds < 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Care on tautomeric bonds */
    // INCHI✔️❌:     if (composite_norm_data)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         check_taut = 1;
    // INCHI✔️❌:         orig_num = (int *) inchi_calloc( (long long)orig_at_data->num_inp_atoms + 2, sizeof( int ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:         curr_num = (int *) inchi_calloc( (long long)orig_at_data->num_inp_atoms + 2, sizeof( int ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:         if (orig_num && curr_num)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             check_taut = 1;
    // INCHI✔️❌:             CompAtomData_GetNumMapping( composite_norm_data, orig_num, curr_num );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: repeatj:
    // INCHI✔️❌:     remove = 0;
    // INCHI✔️❌:     at1 = unit->bkbonds[j][0];
    // INCHI✔️❌:     at2 = unit->bkbonds[j][1];
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     for (k = 0; k < orig_at_data->at[at1 - 1].valence; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (orig_at_data->at[at1 - 1].neighbor[k] != at2 - 1)
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*if ( border > 1 ) */
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     bond_is_untouchable = 0;
    // INCHI✔️❌:     if (check_taut && composite_norm_data && composite_norm_data->at && curr_num) /* djb-rwth: fixing a NULL pointer dereference */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (k = 0; k < composite_norm_data->at[curr_num[at1]].valence; k++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (composite_norm_data->at[curr_num[at1]].neighbor[k] != curr_num[at2])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             btype = composite_norm_data->at[curr_num[at1]].bond_type[k];
    // INCHI✔️❌:             bond_is_untouchable = ( btype == BOND_TAUTOM ); /*|| btype == BOND_ALTERN );*/
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (bond_is_untouchable)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         remove = 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (remove)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         OAD_PolymerUnit_RemoveLinkFromCRUChain( at1, at2, &unit->nbkbonds, unit->bkbonds );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ++j;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (j < unit->nbkbonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto repeatj;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (orig_num)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( orig_num );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (curr_num)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free( curr_num );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: OAD_PolymerUnit_DelistHighOrderBackboneBonds
    // END INCHI C FUNCTION: OAD_PolymerUnit_DelistHighOrderBackboneBonds
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_PolymerUnit_DelistHighOrderBackboneBonds
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; BOND_TAUTOM is 8; calloc zero-initializes each mapping independently; err and pStrErr are not accessed.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_PolymerUnit_DelistHighOrderBackboneBonds

    let initial = heap
        .slice(unit.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if initial.na < 2 || initial.nb < 2 || initial.nbkbonds < 1 {
        return Ok(());
    }
    let mut check_taut = false;
    let mut original_numbers = SourceMutPointer::null();
    let mut current_numbers = SourceMutPointer::null();
    if let Some(composite) = composite_normalized_data {
        check_taut = true;
        let count = i64::from(original_atom_data.num_inp_atoms).wrapping_add(2) as u64;
        original_numbers = match inchi_calloc::<i32>(heap, count, 4) {
            Ok(pointer) => pointer,
            Err(
                SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed,
            ) => SourceMutPointer::null(),
            Err(other) => return Err(other),
        };
        current_numbers = match inchi_calloc::<i32>(heap, count, 4) {
            Ok(pointer) => pointer,
            Err(
                SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed,
            ) => SourceMutPointer::null(),
            Err(other) => return Err(other),
        };
        if !original_numbers.is_null() && !current_numbers.is_null() {
            CompAtomData_GetNumMapping(heap, composite, original_numbers, current_numbers)?;
        }
    }

    let mut position = 0_i32;
    loop {
        let current = heap
            .slice(unit.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if position >= current.nbkbonds {
            break;
        }
        let position_index =
            usize::try_from(position).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let row = *heap
            .slice(current.bkbonds.as_const())?
            .get(position_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let row_values = heap.slice(row.as_const())?;
        let atom1 = *row_values
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atom2 = *row_values
            .get(1)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut untouchable = false;
        if check_taut {
            if let Some(composite) = composite_normalized_data {
                if !composite.at.is_null() && !current_numbers.is_null() {
                    let atom1_index =
                        usize::try_from(atom1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom2_index =
                        usize::try_from(atom2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let mapping = heap.slice(current_numbers.as_const())?;
                    let current1 = *mapping
                        .get(atom1_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let current2 = *mapping
                        .get(atom2_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let current1_index = usize::try_from(current1)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let normalized_atom = heap
                        .slice(composite.at.as_const())?
                        .get(current1_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    for neighbor_index in 0..i32::from(normalized_atom.valence).max(0) {
                        let neighbor_index = usize::try_from(neighbor_index)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        if i32::from(normalized_atom.neighbor[neighbor_index]) != current2 {
                            continue;
                        }
                        untouchable =
                            u32::from(normalized_atom.bond_type[neighbor_index]) == BOND_TAUTOM;
                        break;
                    }
                }
            }
        }
        if untouchable {
            let mut count = current.nbkbonds;
            OAD_PolymerUnit_RemoveLinkFromCRUChain(
                heap,
                atom1,
                atom2,
                &mut count,
                current.bkbonds,
            )?;
            heap.slice_mut(unit)?[0].nbkbonds = count;
        } else {
            position = position.wrapping_add(1);
        }
    }
    if !original_numbers.is_null() {
        inchi_free(heap, original_numbers)?;
    }
    if !current_numbers.is_null() {
        inchi_free(heap, current_numbers)?;
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_FindBackbones(
    heap: &mut SourceHeap,
    atom_data: &mut ORIG_ATOM_DATA,
    composite_normalized_data: Option<&COMP_ATOM_DATA>,
    error: &mut i32,
    mut error_text: Option<&mut [i8]>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2872 OAD_Polymer_FindBackbones
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: OAD_Polymer_FindBackbones
    // INCHI✔️❌: void OAD_Polymer_FindBackbones( ORIG_ATOM_DATA *at_data,
    // INCHI✔️❌:                                 COMP_ATOM_DATA *composite_norm_data,
    // INCHI✔️❌:                                 int            *err,
    // INCHI✔️❌:                                 char           *pStrErr )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:
    // INCHI✔️❌:     *err = 0;
    // INCHI✔️❌:     for (i = 0; i < at_data->polymer->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!at_data->polymer->units[i]->cyclizable)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         OAD_CollectBackboneBonds( at_data,
    // INCHI✔️❌:                                   at_data->polymer->units[i]->na,
    // INCHI✔️❌:                                   at_data->polymer->units[i]->alist,
    // INCHI✔️❌:                                   at_data->polymer->units[i]->end_atom1,
    // INCHI✔️❌:                                   at_data->polymer->units[i]->end_atom2,
    // INCHI✔️❌:                                   &(at_data->polymer->units[i]->nbkbonds),
    // INCHI✔️❌:                                   at_data->polymer->units[i]->bkbonds,
    // INCHI✔️❌:                                   err, pStrErr );
    // INCHI✔️❌:         if (*err)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at_data->polymer->units[i]->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at_data->polymer->units[i]->nbkbonds < 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at_data->polymer->units[i]->nbkbonds == 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  Special case: we got only one bond between end_atom1 and   */
    // INCHI✔️❌:             /*  end_atom2 (this may be the result of metal disconnection)  */
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         OAD_PolymerUnit_DelistIntraRingBackboneBonds( at_data->polymer->units[i], at_data, err, pStrErr );
    // INCHI✔️❌:         if (*err)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         OAD_PolymerUnit_DelistHighOrderBackboneBonds( at_data->polymer->units[i],
    // INCHI✔️❌:                                                       at_data, composite_norm_data,
    // INCHI✔️❌:                                                       err, pStrErr );
    // INCHI✔️❌:         if (*err)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at_data->polymer->units[i]->nbkbonds == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* We already cyclized frame-shiftable unit and preprocessed it (in 'prep_inp_data').                       */
    // INCHI✔️❌:             /* Despite that, now we discovered that there are no bonds eligible for frame shift                         */
    // INCHI✔️❌:             /* (as all potentially eligible in-unit bonds are either in-ring or alternate ones).                        */
    // INCHI✔️❌:             /* We can not simply restore original connections as the structure may have been already heavily touched.   */
    // INCHI✔️❌:             /* The most viable action is to hold a single frame-shift bond (between original caps of CRU ends).   */
    // INCHI✔️❌:             /* It is for sure will be converted to original bonds to star atoms on possible inchi2struct.               */
    // INCHI✔️❌:             at_data->polymer->units[i]->cyclizable = 1;
    // INCHI✔️❌:             at_data->polymer->units[i]->nbkbonds = 1;
    // INCHI✔️❌:             at_data->polymer->units[i]->bkbonds[0][0] = at_data->polymer->units[i]->end_atom1;
    // INCHI✔️❌:             at_data->polymer->units[i]->bkbonds[0][1] = at_data->polymer->units[i]->end_atom2;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: OAD_Polymer_FindBackbones
    // END INCHI C FUNCTION: OAD_Polymer_FindBackbones
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_Polymer_FindBackbones
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; CLOSING_SRU_NOT_APPLICABLE is 0; caller owns each preallocated bkbonds matrix.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_Polymer_FindBackbones

    *error = 0;
    let polymer = heap
        .slice(atom_data.polymer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for index in 0..polymer.n.max(0) {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let unit_pointer = *heap
            .slice(polymer.units.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if unit.cyclizable == 0 {
            continue;
        }
        let listed_count =
            usize::try_from(unit.na.max(0)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom_list = if listed_count == 0 {
            Vec::new()
        } else {
            heap.slice(unit.alist.as_const())?
                .get(..listed_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .to_vec()
        };
        let mut backbone_count = unit.nbkbonds;
        OAD_CollectBackboneBonds(
            heap,
            atom_data,
            unit.na,
            &atom_list,
            unit.end_atom1,
            unit.end_atom2,
            &mut backbone_count,
            unit.bkbonds,
            error,
            error_text.as_deref_mut(),
        )?;
        heap.slice_mut(unit_pointer)?[0].nbkbonds = backbone_count;
        if *error != 0 {
            heap.slice_mut(unit_pointer)?[0].cyclizable = CLOSING_SRU_NOT_APPLICABLE as i32;
            continue;
        }
        if backbone_count < 1 || backbone_count == 1 {
            continue;
        }
        OAD_PolymerUnit_DelistIntraRingBackboneBonds(
            heap,
            unit_pointer,
            atom_data,
            error,
            error_text.as_deref_mut(),
        )?;
        if *error != 0 {
            continue;
        }
        OAD_PolymerUnit_DelistHighOrderBackboneBonds(
            heap,
            unit_pointer,
            atom_data,
            composite_normalized_data,
            error,
            error_text.as_deref_mut(),
        )?;
        if *error != 0 {
            continue;
        }
        let filtered = heap
            .slice(unit_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if filtered.nbkbonds == 0 {
            {
                let unit = &mut heap.slice_mut(unit_pointer)?[0];
                unit.cyclizable = 1;
                unit.nbkbonds = 1;
            }
            let first_row = *heap
                .slice(filtered.bkbonds.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            *heap
                .slice_mut(first_row)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = filtered.end_atom1;
            *heap
                .slice_mut(first_row)?
                .get_mut(1)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = filtered.end_atom2;
        }
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_GetRepresentation(
    heap: &mut SourceHeap,
    polymer: SourceMutPointer<OAD_Polymer>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3788 OAD_Polymer_GetRepresentation
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: OAD_Polymer_GetRepresentation
    // INCHI✔️❌: int  OAD_Polymer_GetRepresentation( OAD_Polymer *p )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, n_source_based_units = 0, n_structure_based_units = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return NO_POLYMER;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < p->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (p->units[i]->nb == 2 || p->units[i]->nbkbonds > 0 ||
    // INCHI✔️❌:             ( ( p->units[i]->cap1 > 0 ) && ( p->units[i]->cap2 > 0 ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             p->units[i]->representation = POLYMER_REPRESENTATION_STRUCTURE_BASED;
    // INCHI✔️❌:             n_structure_based_units++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (p->units[i]->nb == 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             p->units[i]->representation = POLYMER_REPRESENTATION_SOURCE_BASED;
    // INCHI✔️❌:             n_source_based_units++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (p->n == n_source_based_units)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return POLYMER_REPRESENTATION_SOURCE_BASED;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (p->n == n_structure_based_units)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return POLYMER_REPRESENTATION_STRUCTURE_BASED;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (n_source_based_units        &&
    // INCHI✔️❌:               n_structure_based_units &&
    // INCHI✔️❌:               ( n_source_based_units + n_structure_based_units ) == p->n)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* TODO: check if SRU/MON are embedded in a single COP (is this check really necessary? ??) */
    // INCHI✔️❌:         return POLYMER_REPRESENTATION_MIXED;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #if 0
    // INCHI✔️❌:     else if (p->n == ( n_source_based_units + n_structure_based_units ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*
    // INCHI✔️❌:             Structure based presentation may include no-crossing bond units
    // INCHI✔️❌:             which only serve as embedding for ( >1 ) structure-based SRU's.
    // INCHI✔️❌:             The code below accounts for this.
    // INCHI✔️❌:         */
    // INCHI✔️❌:         if (n_source_based_units < n_structure_based_units)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int j, atom, atom_is_shared_with_struct_based_unit = 0;
    // INCHI✔️❌:             for (i = 0; i < p->n; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int k;
    // INCHI✔️❌:                 if (p->units[i]->representation != POLYMER_REPRESENTATION_SOURCE_BASED)
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 for (k = 0; k < p->units[i]->na; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     atom = p->units[i]->alist[k];
    // INCHI✔️❌:                     if (is_in_the_ilist( p->pzz, atom, p->n_pzz ))
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     atom_is_shared_with_struct_based_unit = 0;
    // INCHI✔️❌:                     for (j = 0; j < p->n; j++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (p->units[j]->representation != POLYMER_REPRESENTATION_STRUCTURE_BASED)
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         if (is_in_the_ilist( p->units[j]->alist, atom, p->units[j]->na ))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             atom_is_shared_with_struct_based_unit = 1;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (!atom_is_shared_with_struct_based_unit)
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (!atom_is_shared_with_struct_based_unit)
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (atom_is_shared_with_struct_based_unit)
    // INCHI✔️❌:                 return POLYMER_REPRESENTATION_STRUCTURE_BASED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return POLYMER_REPRESENTATION_MIXED;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     return POLYMER_REPRESENTATION_UNRECOGNIZED;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: OAD_Polymer_GetRepresentation
    // END INCHI C FUNCTION: OAD_Polymer_GetRepresentation
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_Polymer_GetRepresentation
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; the #if 0 embedding analysis is inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_Polymer_GetRepresentation

    if polymer.is_null() {
        return Ok(NO_POLYMER);
    }
    let value = heap
        .slice(polymer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut source_units = 0_i32;
    let mut structure_units = 0_i32;
    for index in 0..value.n.max(0) {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let unit_pointer = *heap
            .slice(value.units.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if unit.nb == 2 || unit.nbkbonds > 0 || (unit.cap1 > 0 && unit.cap2 > 0) {
            heap.slice_mut(unit_pointer)?[0].representation =
                POLYMER_REPRESENTATION_STRUCTURE_BASED as i32;
            structure_units = structure_units.wrapping_add(1);
        } else if unit.nb == 0 {
            heap.slice_mut(unit_pointer)?[0].representation =
                POLYMER_REPRESENTATION_SOURCE_BASED as i32;
            source_units = source_units.wrapping_add(1);
        }
    }
    if value.n == source_units {
        Ok(POLYMER_REPRESENTATION_SOURCE_BASED as i32)
    } else if value.n == structure_units {
        Ok(POLYMER_REPRESENTATION_STRUCTURE_BASED as i32)
    } else if source_units != 0
        && structure_units != 0
        && source_units.wrapping_add(structure_units) == value.n
    {
        Ok(POLYMER_REPRESENTATION_MIXED as i32)
    } else {
        Ok(POLYMER_REPRESENTATION_UNRECOGNIZED as i32)
    }
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OrigAtData_IncreaseBondOrder(
    heap: &mut SourceHeap,
    this_atom: i32,
    other_atom: i32,
    atoms: SourceMutPointer<inp_ATOM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2692 OrigAtData_IncreaseBondOrder
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int  OrigAtData_IncreaseBondOrder( int this_atom, int other_atom, inp_ATOM *at )
{
    int i, k, n_up = 0;
    inp_ATOM *a;

    if (at[this_atom].valence >= MAXVAL ||
         at[other_atom].valence >= MAXVAL)
    {
        return 0;
    }

    a = &( at[this_atom] );
    if (a->chem_bonds_valence > MAXVAL - 1)
    {
        return 0;
    }

    k = a->valence;
    for (i = 0; i < k; i++)
    {
        if (a->neighbor[i] != other_atom)
            continue;
        if (a->bond_type[i] > 3)
            return 0;
        a->bond_type[i]++;
        a->chem_bonds_valence++;
        n_up++;
        break;
    }

    a = &( at[other_atom] );
    if (a->chem_bonds_valence > MAXVAL - 1)
    {
        return 0;
    }
    k = a->valence;

    for (i = 0; i < k; i++)
    {
        if (a->neighbor[i] != this_atom)
        {
            continue;
        }
        if (a->bond_type[i] > 3)
        {
            return 0;
        }
        a->bond_type[i]++;
        a->chem_bonds_valence++;
        n_up++;
        break;
    }

    return n_up;
}
    */
    // END INCHI C FUNCTION: OrigAtData_IncreaseBondOrder
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OrigAtData_IncreaseBondOrder
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; MAXVAL is the generated source constant.
    // INCHI✔️❌: Performance is materially worse because modeled pointer dereference requires SourceHeap allocation-table lookup.
    // END INCHI ACTIVE MACRO CONFIGURATION: OrigAtData_IncreaseBondOrder

    let this_index =
        usize::try_from(this_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let other_index =
        usize::try_from(other_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let values = heap.slice_mut(atoms)?;
    let this_valence = values
        .get(this_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .valence;
    let other_valence = values
        .get(other_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .valence;
    if i32::from(this_valence) >= MAXVAL as i32
        || i32::from(other_valence) >= MAXVAL as i32
    {
        return Ok(0);
    }

    let mut increased = 0_i32;
    {
        let atom = values
            .get_mut(this_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if i32::from(atom.chem_bonds_valence) > MAXVAL as i32 - 1 {
            return Ok(0);
        }
        for bond_index in 0..i32::from(atom.valence) {
            let bond_index = usize::try_from(bond_index)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            if i32::from(atom.neighbor[bond_index]) != other_atom {
                continue;
            }
            if atom.bond_type[bond_index] > 3 {
                return Ok(0);
            }
            atom.bond_type[bond_index] = atom.bond_type[bond_index].wrapping_add(1);
            atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_add(1);
            increased = increased.wrapping_add(1);
            break;
        }
    }

    {
        let atom = values
            .get_mut(other_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if i32::from(atom.chem_bonds_valence) > MAXVAL as i32 - 1 {
            return Ok(0);
        }
        for bond_index in 0..i32::from(atom.valence) {
            let bond_index = usize::try_from(bond_index)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            if i32::from(atom.neighbor[bond_index]) != this_atom {
                continue;
            }
            if atom.bond_type[bond_index] > 3 {
                return Ok(0);
            }
            atom.bond_type[bond_index] = atom.bond_type[bond_index].wrapping_add(1);
            atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_add(1);
            increased = increased.wrapping_add(1);
            break;
        }
    }
    Ok(increased)
}

#[allow(non_snake_case)]
pub(crate) fn OrigAtData_DecreaseBondOrder(
    heap: &mut SourceHeap,
    this_atom: i32,
    other_atom: i32,
    atoms: SourceMutPointer<inp_ATOM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2750 OrigAtData_DecreaseBondOrder
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: OrigAtData_DecreaseBondOrder
    // INCHI✔️❌: int  OrigAtData_DecreaseBondOrder( int      this_atom,
    // INCHI✔️❌:                                    int      other_atom,
    // INCHI✔️❌:                                    inp_ATOM *at )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, k, n_dn = 0;
    // INCHI✔️❌:     inp_ATOM *a;
    // INCHI✔️❌:
    // INCHI✔️❌:     a = &( at[this_atom] );
    // INCHI✔️❌:     if (a->chem_bonds_valence > MAXVAL - 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     k = a->valence;
    // INCHI✔️❌:     for (i = 0; i < k; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (a->neighbor[i] != other_atom)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (a->bond_type[i] < 2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         a->bond_type[i]--;
    // INCHI✔️❌:         a->chem_bonds_valence--;
    // INCHI✔️❌:         n_dn++;
    // INCHI✔️❌:         break;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     a = &( at[other_atom] );
    // INCHI✔️❌:     k = a->valence;
    // INCHI✔️❌:     for (i = 0; i < k; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (a->neighbor[i] != this_atom)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (a->bond_type[i] < 2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         a->bond_type[i]--;
    // INCHI✔️❌:         a->chem_bonds_valence--;
    // INCHI✔️❌:         n_dn++;
    // INCHI✔️❌:         break;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return n_dn;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: OrigAtData_DecreaseBondOrder
    // END INCHI C FUNCTION: OrigAtData_DecreaseBondOrder
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_DecreaseBondOrder
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; MAXVAL is source-defined; neighbor is AT_NUMB/unsigned short and bond_type is U_CHAR.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigAtData_DecreaseBondOrder

    let first_index =
        usize::try_from(this_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let second_index =
        usize::try_from(other_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let first = heap
        .slice(atoms.as_const())?
        .get(first_index)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if i32::from(first.chem_bonds_valence) > MAXVAL as i32 - 1 {
        return Ok(0);
    }
    let mut decreased = 0_i32;
    for neighbor_index in 0..i32::from(first.valence).max(0) {
        let neighbor_index =
            usize::try_from(neighbor_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if i32::from(first.neighbor[neighbor_index]) != other_atom {
            continue;
        }
        if first.bond_type[neighbor_index] < 2 {
            return Ok(0);
        }
        let atom = heap
            .slice_mut(atoms)?
            .get_mut(first_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        atom.bond_type[neighbor_index] = atom.bond_type[neighbor_index].wrapping_sub(1);
        atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_sub(1);
        decreased = decreased.wrapping_add(1);
        break;
    }

    let second = heap
        .slice(atoms.as_const())?
        .get(second_index)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for neighbor_index in 0..i32::from(second.valence).max(0) {
        let neighbor_index =
            usize::try_from(neighbor_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if i32::from(second.neighbor[neighbor_index]) != this_atom {
            continue;
        }
        if second.bond_type[neighbor_index] < 2 {
            return Ok(0);
        }
        let atom = heap
            .slice_mut(atoms)?
            .get_mut(second_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        atom.bond_type[neighbor_index] = atom.bond_type[neighbor_index].wrapping_sub(1);
        atom.chem_bonds_valence = atom.chem_bonds_valence.wrapping_sub(1);
        decreased = decreased.wrapping_add(1);
        break;
    }
    Ok(decreased)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OAD_PolymerUnit_ReopenCyclized(
    heap: &mut SourceHeap,
    unit: SourceMutPointer<OAD_PolymerUnit>,
    atoms: SourceMutPointer<inp_ATOM>,
    _atom_properties: SourceMutPointer<OAD_AtProps>,
    _atom_count: i32,
    input_bond_count: &mut i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3953 OAD_PolymerUnit_ReopenCyclized
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: OAD_PolymerUnit_ReopenCyclized
    // INCHI✔️❌: void OAD_PolymerUnit_ReopenCyclized( OAD_PolymerUnit *u,
    // INCHI✔️❌:                                      inp_ATOM        *at,
    // INCHI✔️❌:                                      OAD_AtProps     *aprops,
    // INCHI✔️❌:                                      int             nat,
    // INCHI✔️❌:                                      int             *num_inp_bonds )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int bond_type, bond_stereo;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (u->cyclizable == CLOSING_SRU_RING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Decyclize artificially introduced bond */
    // INCHI✔️❌:         OrigAtData_RemoveBond( u->end_atom1 - 1, u->end_atom2 - 1,
    // INCHI✔️❌:                                at, &bond_type, &bond_stereo, num_inp_bonds );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (u->cyclizable == CLOSING_SRU_HIGHER_ORDER_BOND)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         OrigAtData_DecreaseBondOrder( u->end_atom1 - 1, u->end_atom2 - 1, at );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (u->cyclizable == CLOSING_SRU_DIRADICAL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[u->end_atom1 - 1].radical == RADICAL_TRIPLET)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             at[u->end_atom1 - 1].radical = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Add explicitly connections to star atoms */
    // INCHI✔️❌:     OrigAtData_AddSingleStereolessBond( u->cap1 - 1, u->end_atom1 - 1,
    // INCHI✔️❌:                                         at, num_inp_bonds );
    // INCHI✔️❌:     OrigAtData_AddSingleStereolessBond( u->cap2 - 1, u->end_atom2 - 1,
    // INCHI✔️❌:                                         at, num_inp_bonds );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Create crossing bonds */
    // INCHI✔️❌:     u->nb = 2;
    // INCHI✔️❌:     u->nbkbonds = 0;
    // INCHI✔️❌:     if (!u->blist)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         u->blist = (int *) inchi_calloc( 2 * (long long)u->nb, sizeof( int ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!u->blist)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     u->blist[0] = u->cap1;
    // INCHI✔️❌:     u->blist[1] = u->end_atom1;
    // INCHI✔️❌:     u->blist[2] = u->cap2;
    // INCHI✔️❌:     u->blist[3] = u->end_atom2;
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: OAD_PolymerUnit_ReopenCyclized
    // END INCHI C FUNCTION: OAD_PolymerUnit_ReopenCyclized
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_PolymerUnit_ReopenCyclized
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; sizeof(int)=4; RADICAL_TRIPLET and CLOSING_SRU_* use generated source constants.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_PolymerUnit_ReopenCyclized

    let initial = heap
        .slice(unit.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if initial.cyclizable == CLOSING_SRU_RING as i32 {
        let mut bond_type = 0_i32;
        let mut bond_stereo = 0_i32;
        let _ = OrigAtData_RemoveBond(
            heap,
            initial.end_atom1.wrapping_sub(1),
            initial.end_atom2.wrapping_sub(1),
            atoms,
            &mut bond_type,
            &mut bond_stereo,
            input_bond_count,
        )?;
    } else if initial.cyclizable == CLOSING_SRU_HIGHER_ORDER_BOND as i32 {
        let _ = OrigAtData_DecreaseBondOrder(
            heap,
            initial.end_atom1.wrapping_sub(1),
            initial.end_atom2.wrapping_sub(1),
            atoms,
        )?;
    } else if initial.cyclizable == CLOSING_SRU_DIRADICAL as i32 {
        let end = usize::try_from(initial.end_atom1.wrapping_sub(1))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom = heap
            .slice_mut(atoms)?
            .get_mut(end)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if i32::from(atom.radical) == RADICAL_TRIPLET as i32 {
            atom.radical = 0;
        }
    }
    let _ = OrigAtData_AddSingleStereolessBond(
        heap,
        initial.cap1.wrapping_sub(1),
        initial.end_atom1.wrapping_sub(1),
        atoms,
        input_bond_count,
    )?;
    let _ = OrigAtData_AddSingleStereolessBond(
        heap,
        initial.cap2.wrapping_sub(1),
        initial.end_atom2.wrapping_sub(1),
        atoms,
        input_bond_count,
    )?;
    {
        let unit = &mut heap.slice_mut(unit)?[0];
        unit.nb = 2;
        unit.nbkbonds = 0;
    }
    let mut bond_list = heap.slice(unit.as_const())?[0].blist;
    if bond_list.is_null() {
        bond_list = match inchi_calloc::<i32>(heap, 4, 4) {
            Ok(pointer) => pointer,
            Err(
                SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
                | SourceHeapError::AllocationFailed,
            ) => SourceMutPointer::null(),
            Err(other) => return Err(other),
        };
        heap.slice_mut(unit)?[0].blist = bond_list;
    }
    if bond_list.is_null() {
        return Ok(());
    }
    let values = heap.slice_mut(bond_list)?;
    *values
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = initial.cap1;
    *values
        .get_mut(1)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = initial.end_atom1;
    *values
        .get_mut(2)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = initial.cap2;
    *values
        .get_mut(3)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = initial.end_atom2;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_SmartReopenCyclizedUnits(
    heap: &mut SourceHeap,
    polymer: SourceMutPointer<OAD_Polymer>,
    atoms: SourceMutPointer<inp_ATOM>,
    atom_count: i32,
    input_bond_count: &mut i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3878 OAD_Polymer_SmartReopenCyclizedUnits
    // BEGIN COMPLETE VERBATIM SOURCE FRAME: OAD_Polymer_SmartReopenCyclizedUnits
    // INCHI✔️❌: void OAD_Polymer_SmartReopenCyclizedUnits( OAD_Polymer *p,
    // INCHI✔️❌:                                            inp_ATOM    *at,
    // INCHI✔️❌:                                            int         nat,
    // INCHI✔️❌:                                            int         *num_inp_bonds )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i;
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #68329 */
    // INCHI✔️❌:     OAD_AtProps *aprops = (OAD_AtProps*)inchi_calloc((long long)nat + 1, sizeof(OAD_AtProps)); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     /* nat + 1: add extra element for possibe 1-based indexing */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(aprops); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (p->n < 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(aprops); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!p->really_do_frame_shift)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(aprops); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #68329 */
    // INCHI✔️❌:     if (nat <= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(aprops); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*ITRACE_( "\n\n*********************************************************************\n* ENTERING OAD_Polymer_SmartReopenCyclizedUnits()" );
    // INCHI✔️❌:     OAD_Polymer_DebugTrace( p );*/
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Set atom properties for sorting */
    // INCHI✔️❌:     if (!aprops || !at) /* djb-rwth: fixing oss-fuzz issue #68329, #68286 */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(aprops); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     OAD_Polymer_SetAtProps( p, at, nat, num_inp_bonds, aprops, NULL ); /* NULL as we alredy are in 1-based cano_nums while at i2s/i2i*/
    // INCHI✔️❌:     for (i = 0; i < p->n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (p->units[i]) /* djb-rwth: fixing oss-fuzz issue #68329 */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             OAD_PolymerUnit *u = p->units[i];
    // INCHI✔️❌:             if (p->frame_shift_scheme == FSS_NONE)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ( /* !u->cyclizable || u->cyclized  || */
    // INCHI✔️❌:                 u->nbkbonds < 1 ||
    // INCHI✔️❌:                 u->cap1 < 1 || u->cap2 < 1 ||
    // INCHI✔️❌:                 u->cap1 > nat || u->cap2 > nat)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (OAD_PolymerUnit_SetReopeningDetails(u, at))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int senior_bond;
    // INCHI✔️❌:                 OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(u, at, aprops, &senior_bond);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             OAD_PolymerUnit_ReopenCyclized(u, at, aprops, nat, num_inp_bonds);
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     p->really_do_frame_shift = 0;
    // INCHI✔️❌:     inchi_free( aprops );
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END COMPLETE VERBATIM SOURCE FRAME: OAD_Polymer_SmartReopenCyclizedUnits
    // END INCHI C FUNCTION: OAD_Polymer_SmartReopenCyclizedUnits
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_Polymer_SmartReopenCyclizedUnits
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; sizeof(OAD_AtProps)=16; FSS_NONE is generated source enum value 1.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OAD_Polymer_SmartReopenCyclizedUnits

    let allocation_count = i64::from(atom_count).wrapping_add(1) as u64;
    let atom_properties = match inchi_calloc::<OAD_AtProps>(heap, allocation_count, 16) {
        Ok(pointer) => pointer,
        Err(
            SourceHeapError::AllocationSizeOverflow
            | SourceHeapError::AllocationElementCountOutOfRange
            | SourceHeapError::AllocationFailed,
        ) => SourceMutPointer::null(),
        Err(other) => return Err(other),
    };
    if polymer.is_null() {
        inchi_free(heap, atom_properties)?;
        return Ok(());
    }
    let polymer_value = heap
        .slice(polymer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if polymer_value.n < 1
        || polymer_value.really_do_frame_shift == 0
        || atom_count <= 0
        || atom_properties.is_null()
        || atoms.is_null()
    {
        inchi_free(heap, atom_properties)?;
        return Ok(());
    }
    OAD_Polymer_SetAtProps(
        heap,
        polymer.as_const(),
        atoms,
        atom_count,
        input_bond_count,
        atom_properties,
        SourceConstPointer::null(),
    )?;
    for index in 0..polymer_value.n.max(0) {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let unit_pointer = *heap
            .slice(polymer_value.units.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if unit_pointer.is_null() {
            continue;
        }
        if polymer_value.frame_shift_scheme == tagFrameShifScheme_FSS_NONE as i32 {
            continue;
        }
        let mut unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if unit.nbkbonds < 1
            || unit.cap1 < 1
            || unit.cap2 < 1
            || unit.cap1 > atom_count
            || unit.cap2 > atom_count
        {
            continue;
        }
        let reopening =
            OAD_PolymerUnit_SetReopeningDetails(heap, &mut unit, heap.slice(atoms.as_const())?)?;
        if reopening != 0 {
            let properties = heap.slice(atom_properties.as_const())?.to_vec();
            let mut senior_bond = 0_i32;
            OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(
                heap,
                &mut unit,
                atoms,
                &properties,
                &mut senior_bond,
            )?;
        }
        heap.slice_mut(unit_pointer)?[0] = unit;
        OAD_PolymerUnit_ReopenCyclized(
            heap,
            unit_pointer,
            atoms,
            atom_properties,
            atom_count,
            input_bond_count,
        )?;
    }
    heap.slice_mut(polymer)?[0].really_do_frame_shift = 0;
    inchi_free(heap, atom_properties)?;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_SetReopeningDetails(
    heap: &SourceHeap,
    unit: &mut OAD_PolymerUnit,
    atoms: &[inp_ATOM],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4007 OAD_PolymerUnit_SetReopeningDetails
    // INCHI✔️✔️: int OAD_PolymerUnit_SetReopeningDetails( OAD_PolymerUnit *u, inp_ATOM *at )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int k;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Check reopening  type */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Caps are separated by one atom - that's not error but do nothing */
    // INCHI✔️✔️:     if (u->nbkbonds == 0)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (u->nbkbonds == 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         u->end_atom1 = u->bkbonds[0][0];
    // INCHI✔️✔️:         u->end_atom2 = u->bkbonds[0][1];
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if (u->end_atom1 == u->end_atom2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️: #ifdef ALLOW_CLOSING_SRU_VIA_DIRADICAL
    // INCHI✔️✔️:             u->cyclizable = CLOSING_SRU_DIRADICAL;
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:             u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             /* If caps are separated by two atoms - that's not error but do nothing */
    // INCHI✔️✔️:             for (k = 0; k < at[u->end_atom1 - 1].valence; k++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (at[u->end_atom1 - 1].neighbor[k] == u->end_atom2 - 1)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     if (at[u->end_atom1 - 1].bond_type[k] > 1)
    // INCHI✔️✔️: #ifdef ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND
    // INCHI✔️✔️:                         u->cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;
    // INCHI✔️✔️: #else
    // INCHI✔️✔️: /*                  u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;*/
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️: break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     } /*    */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return u->nbkbonds;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_SetReopeningDetails

    if unit.nbkbonds == 0 {
        return Ok(0);
    }
    if unit.nbkbonds == 1 {
        let first_row = *heap
            .slice(unit.bkbonds.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let row = heap.slice(first_row.as_const())?;
        unit.end_atom1 = *row.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
        unit.end_atom2 = *row.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if unit.end_atom1 == unit.end_atom2 {
            unit.cyclizable = CLOSING_SRU_DIRADICAL as i32;
        } else {
            let atom_index = usize::try_from(
                unit.end_atom1
                    .checked_sub(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?,
            )
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom = atoms
                .get(atom_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let target = unit
                .end_atom2
                .checked_sub(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let valence = usize::try_from(i32::from(atom.valence).max(0))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            for index in 0..valence {
                if i32::from(atom.neighbor[index]) == target {
                    if atom.bond_type[index] > 1 {
                        unit.cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND as i32;
                    }
                    break;
                }
            }
        }
    }
    Ok(unit.nbkbonds)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(
    heap: &mut SourceHeap,
    unit: &mut OAD_PolymerUnit,
    _atoms: SourceMutPointer<inp_ATOM>,
    atom_properties: &[OAD_AtProps],
    senior_bond: &mut i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4056 OAD_PolymerUnit_SortBackboneBondsAndSetSeniors
    // INCHI✔️❌: void OAD_PolymerUnit_SortBackboneBondsAndSetSeniors( OAD_PolymerUnit *u,
    // INCHI✔️❌:                                                      inp_ATOM        *at,
    // INCHI✔️❌:                                                      OAD_AtProps     *aprops,
    // INCHI✔️❌:                                                      int             *senior_bond )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int j, *bnum = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     *senior_bond = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Sort backbone (== frame shiftable) bonds if necessary */
    // INCHI✔️❌:     if (u->nbkbonds > 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bnum = (int *) inchi_calloc( u->nbkbonds, sizeof( int ) );
    // INCHI✔️❌:         if (bnum)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (j = 0; j < u->nbkbonds; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bnum[j] = j;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             OAD_PolymerUnit_SortBackboneBonds( u, aprops, bnum );
    // INCHI✔️❌:             *senior_bond = bnum[0];
    // INCHI✔️❌:             inchi_free( bnum );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* v. 1.05+ : place senior atom the first ("left") in the senior bond */
    // INCHI✔️❌:     if (OAD_Polymer_IsFirstAtomRankLower( u->bkbonds[*senior_bond][0], u->bkbonds[*senior_bond][1], aprops ) == 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int tmp = u->bkbonds[*senior_bond][0];
    // INCHI✔️❌:         u->bkbonds[*senior_bond][0] = u->bkbonds[*senior_bond][1];
    // INCHI✔️❌:         u->bkbonds[*senior_bond][1] = tmp;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     u->end_atom1 = u->bkbonds[*senior_bond][0];
    // INCHI✔️❌:     u->end_atom2 = u->bkbonds[*senior_bond][1];
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_SortBackboneBondsAndSetSeniors

    *senior_bond = 0;
    if unit.nbkbonds > 1 {
        match inchi_calloc::<i32>(heap, unit.nbkbonds as u64, 4) {
            Ok(bond_number_pointer) => {
                for (index, value) in heap.slice_mut(bond_number_pointer)?.iter_mut().enumerate() {
                    *value =
                        i32::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                }
                let sort_result = OAD_PolymerUnit_SortBackboneBonds(
                    heap,
                    unit,
                    atom_properties,
                    bond_number_pointer,
                );
                if sort_result.is_ok() {
                    *senior_bond = *heap
                        .slice(bond_number_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                }
                inchi_free(heap, bond_number_pointer)?;
                sort_result?;
            }
            Err(SourceHeapError::AllocationFailed) => {}
            Err(error) => return Err(error),
        }
    }

    let senior_index =
        usize::try_from(*senior_bond).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let row_pointer = *heap
        .slice(unit.bkbonds.as_const())?
        .get(senior_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let (first, second) = {
        let row = heap.slice(row_pointer.as_const())?;
        (
            *row.first().ok_or(SourceHeapError::PointerOutOfBounds)?,
            *row.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    if OAD_Polymer_IsFirstAtomRankLower(first, second, atom_properties)? == 1 {
        let row = heap.slice_mut(row_pointer)?;
        row[0] = second;
        row[1] = first;
    }
    let row = heap.slice(row_pointer.as_const())?;
    unit.end_atom1 = row[0];
    unit.end_atom2 = row[1];
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_SortBackboneBonds(
    heap: &mut SourceHeap,
    unit: &OAD_PolymerUnit,
    atom_properties: &[OAD_AtProps],
    bond_numbers: SourceMutPointer<i32>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4097 OAD_PolymerUnit_SortBackboneBonds
    // INCHI✔️❌: void OAD_PolymerUnit_SortBackboneBonds( OAD_PolymerUnit *u,
    // INCHI✔️❌:                                         OAD_AtProps     *aprops,
    // INCHI✔️❌:                                         int             *bnum )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, tmp;
    // INCHI✔️❌:     int n = u->nbkbonds;
    // INCHI✔️❌:     if (NULL == bnum)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 1; i < n; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         tmp = bnum[i];
    // INCHI✔️❌:         j = i - 1;
    // INCHI✔️❌:         while (j >= 0 && OAD_Polymer_CompareBackboneBondsSeniority( u->bkbonds[bnum[j]], u->bkbonds[tmp], aprops ) > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             bnum[j + 1] = bnum[j];
    // INCHI✔️❌:             j--;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         bnum[j + 1] = tmp;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_SortBackboneBonds

    if bond_numbers.is_null() {
        return Ok(());
    }
    let bond_at = |heap: &SourceHeap, number: i32| -> Result<[i32; 2], SourceHeapError> {
        let index = usize::try_from(number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let row_pointer = *heap
            .slice(unit.bkbonds.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let row = heap.slice(row_pointer.as_const())?;
        Ok([
            *row.first().ok_or(SourceHeapError::PointerOutOfBounds)?,
            *row.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?,
        ])
    };

    for i in 1..unit.nbkbonds {
        let i = usize::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let temporary = *heap
            .slice(bond_numbers.as_const())?
            .get(i)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut j = i32::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)? - 1;
        while j >= 0 {
            let j_index = usize::try_from(j).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let current = *heap
                .slice(bond_numbers.as_const())?
                .get(j_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if OAD_Polymer_CompareBackboneBondsSeniority(
                &bond_at(heap, current)?,
                &bond_at(heap, temporary)?,
                atom_properties,
            )? <= 0
            {
                break;
            }
            heap.slice_mut(bond_numbers)?[j_index + 1] = current;
            j -= 1;
        }
        let destination =
            usize::try_from(j + 1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        heap.slice_mut(bond_numbers)?[destination] = temporary;
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_CompareBackboneBondsSeniority(
    bond1: &[i32; 2],
    bond2: &[i32; 2],
    atom_properties: &[OAD_AtProps],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4128 OAD_Polymer_CompareBackboneBondsSeniority
    // INCHI✔️✔️: int  OAD_Polymer_CompareBackboneBondsSeniority( int* b1, int* b2, OAD_AtProps *aprops )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int b1min, b1max, b2min, b2max, tmp, cmp = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Find min and max ext-ranked ends of the both bonds */
    // INCHI✔️✔️:     b1max = b1[0]; b1min = b1[1];
    // INCHI✔️✔️:     b2max = b2[0]; b2min = b2[1];
    // INCHI✔️✔️:     if (OAD_Polymer_IsFirstAtomRankLower( b1min, b1max, aprops ) == -1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         tmp = b1max;
    // INCHI✔️✔️:         b1max = b1min;
    // INCHI✔️✔️:         b1min = tmp;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (OAD_Polymer_IsFirstAtomRankLower( b2min, b2max, aprops ) == -1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         tmp = b2max;
    // INCHI✔️✔️:         b2max = b2min;
    // INCHI✔️✔️:         b2min = tmp;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Compare bonds' seniority */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* First, favor the bond which has greater ext-rank end
    // INCHI✔️✔️:        NB: the result may be 0, that is, equal max ext. ranks */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     cmp = OAD_Polymer_CompareRanksOfTwoAtoms( b1max, b2max, aprops );
    // INCHI✔️✔️:     if (cmp == 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return   1;        /* rank(b1max) < rank(b2max), so bond2 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (cmp == -1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  -1;        /* rank(b1max) > rank(b2max), so bond1 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Max ends are of the same rank, so favor the bond with lesser min-rank end
    // INCHI✔️✔️:        NB: the result may NOT be 0, that is, the case is always resolved    */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     cmp = OAD_Polymer_CompareRanksOfTwoAtoms( b1min, b2min, aprops ); /*OAD_Polymer_IsFirstAtomRankLower( b1min, b2min, aprops );*/
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (cmp == 1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  -1;         /* rank(b1min) < rank(b2min), so bond1 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (cmp == -1)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return   1;         /* rank(b1min) > rank(b2min), so bond2 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Min ends are of the same rank. Here is the time to compare directly
    // INCHI✔️✔️:        which canonical number is larger of max-ends ... */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (b1max < b2max)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (b1max > b2max)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* ... they are the same, so compare which canonical number is larger for min-ends ... */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (b1min < b2min)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;          /* b1min < b2min, so bond1 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (b1min > b2min)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;          /* b1min > b2min, so bond2 is senior */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;    /* we should not reach there */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OAD_Polymer_CompareBackboneBondsSeniority

    let (mut bond1_max, mut bond1_min) = (bond1[0], bond1[1]);
    let (mut bond2_max, mut bond2_min) = (bond2[0], bond2[1]);
    if OAD_Polymer_IsFirstAtomRankLower(bond1_min, bond1_max, atom_properties)? == -1 {
        std::mem::swap(&mut bond1_max, &mut bond1_min);
    }
    if OAD_Polymer_IsFirstAtomRankLower(bond2_min, bond2_max, atom_properties)? == -1 {
        std::mem::swap(&mut bond2_max, &mut bond2_min);
    }

    let comparison = OAD_Polymer_CompareRanksOfTwoAtoms(bond1_max, bond2_max, atom_properties)?;
    if comparison == 1 {
        return Ok(1);
    } else if comparison == -1 {
        return Ok(-1);
    }

    let comparison = OAD_Polymer_CompareRanksOfTwoAtoms(bond1_min, bond2_min, atom_properties)?;
    if comparison == 1 {
        return Ok(-1);
    } else if comparison == -1 {
        return Ok(1);
    }

    if bond1_max < bond2_max {
        return Ok(1);
    }
    if bond1_max > bond2_max {
        return Ok(-1);
    }
    if bond1_min < bond2_min {
        return Ok(-1);
    }
    if bond1_min > bond2_min {
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_CompareRanksOfTwoAtoms(
    atom1: i32,
    atom2: i32,
    atom_properties: &[OAD_AtProps],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4209 OAD_Polymer_CompareRanksOfTwoAtoms
    // INCHI✔️✔️: int OAD_Polymer_CompareRanksOfTwoAtoms( int atom1, int atom2, OAD_AtProps *aprops )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     const int HETEROCYC = 3, HETEROAT = 2, CARBOCYC = 1, CARBOAT = 0;
    // INCHI✔️✔️:         /* NB: Carbon's rank is always 2, next to the lowest */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     int a1 = atom1 - 1;
    // INCHI✔️✔️:     int a2 = atom2 - 1;
    // INCHI✔️✔️:     int a1typ = CARBOAT;
    // INCHI✔️✔️:     int a2typ = CARBOAT;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* djb-rwth: fixing oss-fuzz issue #69501, #68277 */
    // INCHI✔️✔️:     if ((a1 < 0) || (a2 < 0))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (aprops[a1].ring_size > 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (aprops[a1].ring_erank <= 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a1typ = CARBOCYC;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a1typ = HETEROCYC;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (aprops[a1].erank == 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a1typ = CARBOAT;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a1typ = HETEROAT;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (aprops[a2].ring_size > 2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (aprops[a2].ring_erank <= 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a2typ = CARBOCYC;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a2typ = HETEROCYC;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (aprops[a2].erank == 2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a2typ = CARBOAT;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         else
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             a2typ = HETEROAT;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Compare */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:         Follow IUPAC Rule 1
    // INCHI✔️✔️:             'The basic order of seniority of subunits is:
    // INCHI✔️✔️:                 heterocyclic rings and ring systems > heteroatom chains >
    // INCHI✔️✔️:                     > carbocyclic rings and ring systems > acyclic carbon chains'
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (a1typ == HETEROCYC && a2typ == HETEROCYC)   /* a1 and a2 are HETEROCYC */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* Try resolving by senior-heteroatom ring */
    // INCHI✔️✔️:         if (aprops[a1].ring_erank < aprops[a2].ring_erank)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return  1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (aprops[a1].ring_erank > aprops[a2].ring_erank)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* Same senior-heteroatom rings, try resolving by total ring size */
    // INCHI✔️✔️:         if (aprops[a1].ring_size < aprops[a2].ring_size)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return  1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (aprops[a1].ring_size > aprops[a2].ring_size)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* Could not resolve... */
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a1typ == HETEROCYC)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;  /* a1 is HETEROCYC, a2 is any other (==junior) */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a2typ == HETEROCYC)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;  /* a2 is HETEROCYC, a1 is any other (==junior) */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* HETEROCYC left out */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (a1typ == HETEROAT && a2typ == HETEROAT)  /* a1 and a2 are HETEROAT */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (aprops[a1].erank < aprops[a2].erank)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return  1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (aprops[a1].erank > aprops[a2].erank)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* Could not resolve... */
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a1typ == HETEROAT)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;  /* a1 is HETEROAT, a2 is any other (==junior) */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a2typ == HETEROAT)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;  /* a2 is HETEROAT, a1 is any other (==junior) */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* HETEROAT left out */
    // INCHI✔️✔️:     if (a1typ == CARBOCYC && a2typ == CARBOCYC) /* a1 and a2 are CARBOCYC */
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* Same senior-atom (C) ring, try resolving by total ring size */
    // INCHI✔️✔️:         if (aprops[a1].ring_size < aprops[a2].ring_size)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return  1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (aprops[a1].ring_size > aprops[a2].ring_size)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* Could not resolve... */
    // INCHI✔️✔️:         return 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a1typ == CARBOCYC)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (a2typ == CARBOCYC)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;        /* 0 means unresolved. It is legal here */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OAD_Polymer_CompareRanksOfTwoAtoms

    const HETEROCYC: i32 = 3;
    const HETEROAT: i32 = 2;
    const CARBOCYC: i32 = 1;
    const CARBOAT: i32 = 0;

    let index1 = atom1
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let index2 = atom2
        .checked_sub(1)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if index1 < 0 || index2 < 0 {
        return Ok(0);
    }
    let property1 = atom_properties
        .get(usize::try_from(index1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let property2 = atom_properties
        .get(usize::try_from(index2).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let classify = |property: &OAD_AtProps| {
        if property.ring_size > 2 {
            if property.ring_erank <= 2 {
                CARBOCYC
            } else {
                HETEROCYC
            }
        } else if property.erank == 2 {
            CARBOAT
        } else {
            HETEROAT
        }
    };
    let type1 = classify(property1);
    let type2 = classify(property2);
    if type1 == HETEROCYC && type2 == HETEROCYC {
        if property1.ring_erank < property2.ring_erank {
            return Ok(1);
        }
        if property1.ring_erank > property2.ring_erank {
            return Ok(-1);
        }
        if property1.ring_size < property2.ring_size {
            return Ok(1);
        }
        if property1.ring_size > property2.ring_size {
            return Ok(-1);
        }
        return Ok(0);
    } else if type1 == HETEROCYC {
        return Ok(-1);
    } else if type2 == HETEROCYC {
        return Ok(1);
    }
    if type1 == HETEROAT && type2 == HETEROAT {
        if property1.erank < property2.erank {
            return Ok(1);
        }
        if property1.erank > property2.erank {
            return Ok(-1);
        }
        return Ok(0);
    } else if type1 == HETEROAT {
        return Ok(-1);
    } else if type2 == HETEROAT {
        return Ok(1);
    }
    if type1 == CARBOCYC && type2 == CARBOCYC {
        if property1.ring_size < property2.ring_size {
            return Ok(1);
        }
        if property1.ring_size > property2.ring_size {
            return Ok(-1);
        }
        return Ok(0);
    } else if type1 == CARBOCYC {
        return Ok(-1);
    } else if type2 == CARBOCYC {
        return Ok(1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn OAD_Polymer_IsFirstAtomRankLower(
    atom1: i32,
    atom2: i32,
    atom_properties: &[OAD_AtProps],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4369 OAD_Polymer_IsFirstAtomRankLower
    // INCHI✔️✔️: int OAD_Polymer_IsFirstAtomRankLower( int atom1, int atom2, OAD_AtProps *aprops )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     /* Compare ext-ranks */
    // INCHI✔️✔️:     int result = OAD_Polymer_CompareRanksOfTwoAtoms( atom1, atom2, aprops );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (result)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return result;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Could not resolve who is junior by extended-ranks...             */
    // INCHI✔️✔️:     /* As a last resort, simply check which canonical number is lesser  */
    // INCHI✔️✔️:     if (atom1 < atom2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return  1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (atom1 > atom2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* should not reach there */
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: OAD_Polymer_IsFirstAtomRankLower

    let result = OAD_Polymer_CompareRanksOfTwoAtoms(atom1, atom2, atom_properties)?;
    if result != 0 {
        return Ok(result);
    }
    if atom1 < atom2 {
        return Ok(1);
    }
    if atom1 > atom2 {
        return Ok(-1);
    }
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_ValidatePolymerAndPseudoElementData(
    heap: &mut SourceHeap,
    original_atom_data: &mut ORIG_ATOM_DATA,
    treat_polymers: i32,
    allow_non_polymer_zz: i32,
    mut error_text: Option<&mut [i8]>,
    no_warnings: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1516 OAD_ValidatePolymerAndPseudoElementData
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int OAD_ValidatePolymerAndPseudoElementData( ORIG_ATOM_DATA *orig_at_data,
                                             int treat_polymers,
                                             int bNPZz,
                                             char *pStrErr,
                                             int bNoWarnings)
{
    int i, k, kk, type, subtype, representation, err = 0;
    int nat = orig_at_data->num_inp_atoms;
    int nsgroups = 0;
    OAD_PolymerUnit* u = NULL;
    OAD_Polymer *pd = orig_at_data->polymer;


    /* Assign polymer type and subunits type and check polymer data for consistency */
    /* djb-rwth: addressing coverity ID #499497 -- TREAT_ERR properly used in all cases */
    
    orig_at_data->valid_polymer = 0;
    if (treat_polymers && pd)
    {
        orig_at_data->valid_polymer = 1;
    }
    if (orig_at_data->valid_polymer)
    {
        nsgroups = pd->n;
    }
    if (nsgroups == 1)
    {
        /* Check if copolymer */
        type = pd->units[0]->type;
        if (type == POLYMER_STY_COP)
        {
            TREAT_ERR( err, 9001, "Copolymer must contain more than one unit" );
            goto exit_function;
        }
        /* Check if copolymer subtype */
        subtype = pd->units[0]->subtype;
        if (subtype == POLYMER_SST_RAN || subtype == POLYMER_SST_ALT || subtype == POLYMER_SST_BLK)
        {
            TREAT_ERR( err, 9002, "Single polymer unit may not be RAN/ALT/BLO" );
            goto exit_function;
        }
    }

    /* For each CRU */
    for (i = 0; i < nsgroups; i++)
    {
        /* Check if unit data makes sense */
        u = pd->units[i];
        if (u->nb != 0 && u->nb != 2)
        {
            TREAT_ERR( err, 9003, "Number of crossing bonds in polymer unit is not 0 or 2" );
            goto exit_function;
        }
        if (u->na < 1)
        {
            TREAT_ERR( err, 9004, "Empty polymer unit" );
            goto exit_function;
        }
        if (u->na > nat)
        {
            TREAT_ERR( err, 9005, "Too large polymer unit" );
            goto exit_function;
        }
        for (k = 0; k < u->na; k++)
        {
            int atom = u->alist[k];
            if (atom < 1 || atom > nat)
            {
                TREAT_ERR( err, 9006, "Invalid atom number in polymer unit" );
                goto exit_function;
            }
            /* was not accounting for COP ...
            if (is_in_the_ilist( pd->pzz, atom, pd->n_pzz ))
            {
                TREAT_ERR( err, 9007, "Star atom inside polymer unit" );
                goto exit_function;
            }
            */
        }

        OAD_PolymerUnit_SetEndsAndCaps( u, orig_at_data, &err, pStrErr );
            /*	Reveal and store CRU caps and ends('stars and partners')
                Also set `unit->cap1_is_undef`, `unit->cap2_is_undef`, `unit->cyclizable` 
            */
        if (err)
        {
            goto exit_function;
        }

        
        /* Set possibly missing unit parameters */
        u->nbkbonds = 0;
        u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
        u->cyclized = 0;
    }
    

    OAD_ValidateAndSortOutPseudoElementAtoms( orig_at_data, treat_polymers, bNPZz, &err, pStrErr );
    /* Here we:
                Make more polymer and pseudoatom data checks
                Convert both "*" and "Zz" temporarily to "Zy" (polymer-unrelated interal pseudoatoms)
                If applicable, check each CRU and back-convert "Zy" to "Zz" (polymer-related 
                pseudoelement atoms) if they are for valid bi-undef-end CRU
    */
    
    if (err)
    {
        /* already treated TREAT_ERR( err, 9040, "Improper pseudoelement atoms" ); */
        goto exit_function;
    }


    /* Make more polymer and pseudoatom data checks */

    /* Check if non-polymer-related Zz/star atoms enabled */
    if (orig_at_data->n_zy > 0 && bNPZz==0 )
    {
        TREAT_ERR( err, 9, "Non-polymer-related Zz/star atoms are not allowed" );
        goto exit_function;
    }

    if (!pd || !orig_at_data->valid_polymer )
    {
        goto exit_function;
    }

    if (pd->n_pzz > 0)
    {
        /* Allocate memory for polymer-related pseudoatoms */
        if ( pd->treat == POLYMERS_NO)
        {
            TREAT_ERR( err, 9, "Pseudoelement endgroups are not allowed" );
            goto exit_function;
        }
        if (pd->pzz)
        {
            inchi_free(pd->pzz);
            pd->pzz = NULL;
        }
        pd->pzz = (int *) inchi_calloc( pd->n_pzz, sizeof( int ) );
        if (!pd->pzz)
        {
            TREAT_ERR( err, 9010, "Not enough memory" );
            goto exit_function;
        }
        kk = 0;
        for (k = 0; k < nat; k++)
        {
            if (!strcmp( orig_at_data->at[k].elname, "Zz" ))
            {
                pd->pzz[kk++] = k + 1; /* djb-rwth: buffer overrun avoided implicitly */
            }
        }
    }

    /* Check copolymers and ensure that COP includes > 1 SRU */
    for (i = 0; i < pd->n; i++)
    {
        u = pd->units[i];

        if ( u->type == POLYMER_STY_COP ||
             u->type == POLYMER_STY_SRU     /* what drawn as 'SRU' [xyz]n may be actually copolymer [xyz]co */
           )
        {
            int j, in_units = 0;

            if (u->nb > 0)
            {
                /* crossing bonds present, either valid SRU or invalid copolymer */
                if (u->type == POLYMER_STY_COP)
                {
                    TREAT_ERR( err, 9026, "Polymer COP unit contains bracket-crossing bonds, not supported" );
                    goto exit_function;
                }
                else
                {
                    continue;
                }
            }
            /* now we have no crossing bonds units */
            for (j = 0; j < pd->n; j++)
            {
                if (pd->units[j]->type == POLYMER_STY_COP)
                {
                    continue;
                }
                if (is_ilist_inside( pd->units[j]->alist, pd->units[j]->na, pd->units[i]->alist, pd->units[i]->na ))
                {
                    in_units++;
                    if (in_units == 2)
                    {
                        break;
                    }
                }
            }
            if (in_units > 1)
            {
                if (u->type != POLYMER_STY_COP)
                {
                    u->type = POLYMER_STY_COP;
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErr, "Convert multiple-subunits unit to copolymer" );
                    }
                }
            }
            else    /* in_units <= 1)*/
            {
                if (u->type == POLYMER_STY_COP)
                {
                    TREAT_ERR( err, 9027, "Polymer COP unit contains a single SRU instead of multiple" );
                    goto exit_function;
                }
            }
        }
    }

    representation = OAD_Polymer_GetRepresentation( pd );

    /* Make more polymer data checks and perform some corrections*/
    if (representation == POLYMER_REPRESENTATION_SOURCE_BASED)
    {
        for (i = 0; i < nsgroups; i++)
        {
            /* Replace source-based 'SRU' with 'MON' */
            if (pd->units[i]->type == POLYMER_STY_SRU)
            {
                pd->units[i]->type = POLYMER_STY_MON;
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErr, "Converted src-based polymer unit type to MON" );
                }
            }
            if (pd->units[i]->type == POLYMER_STY_COP)
            {
                /* Set missing copolymer subtype to RAN */
                if (pd->units[i]->subtype == POLYMER_SST_NON)
                {
                    pd->units[i]->subtype = POLYMER_SST_RAN;
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErr, "Set missing copolymer subtype to RAN" );
                    }
                }
            }
            /* Suppress connectivity ("HH", "HT", "EU") */
            if (pd->units[i]->conn != POLYMER_CONN_NON)
            {
                pd->units[i]->conn = POLYMER_CONN_NON;
                if (!bNoWarnings)
                {
                    WarningMessage( pStrErr, "Ignore connection pattern for src-based polymer unit" );
                }
            }
        }
    }

#ifdef ALLOW_MIXED_SRU_AND_MON
    else if (representation == POLYMER_REPRESENTATION_STRUCTURE_BASED ||
              representation == POLYMER_REPRESENTATION_MIXED)
#else
    else if (representation == POLYMER_REPRESENTATION_STRUCTURE_BASED)
#endif
    {
        for (i = 0; i < nsgroups; i++)
        {
            int a1, a2, a1_is_not_in_alist, a1_is_star_atom, a2_is_not_in_alist, a2_is_star_atom;

            u = pd->units[i];

            /*    SRU that is copolymer unit embedding other SRU's */
            if (u->nb == 0)
            {
                if (u->type == POLYMER_STY_COP)
                {
                    ;
                }
                else if (u->type == POLYMER_STY_SRU)
                {
                    u->type = POLYMER_STY_COP;
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErr, "Set copolymer embedding unit mark to COP" );
                    }
                }
            }
            if (u->type == POLYMER_STY_COP)
            {
                u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
                /* Set possibly missing copolymer subtype to RAN */
                if (u->subtype == POLYMER_SST_NON)
                {
                    u->subtype = POLYMER_SST_RAN;
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErr, "Set missing copolymer subtype to RAN" );
                    }
                }
                continue;
            }

#ifdef ALLOW_MIXED_SRU_AND_MON
            if (u->type == POLYMER_STY_MON)
            {
                continue;
            }
#endif
            /*    SRU with endgroups or stars. Check it. */
            for (k = 0; k < u->nb; k++)
            {
                /* Check that there are no H end groups */
                a1 = u->blist[2 * k]; a2 = u->blist[2 * k + 1];
                if (!strcmp( orig_at_data->at[a1 - 1].elname, "H" ) ||
                     !strcmp( orig_at_data->at[a1 - 1].elname, "D" ) ||
                     !strcmp( orig_at_data->at[a1 - 1].elname, "T" ))
                {
                    TREAT_ERR( err, 9030, "H as polymer end group is not supported" );
                    goto exit_function;
                }
                if ( !strcmp( orig_at_data->at[a2 - 1].elname, "H" ) ||
                     !strcmp( orig_at_data->at[a2 - 1].elname, "D" ) ||
                     !strcmp( orig_at_data->at[a2 - 1].elname, "T" ))
                {
                    TREAT_ERR( err, 9031, "H as polymer end group is not supported" );
                    goto exit_function;
                }
                /* Ensure that caps of polymer unit lie outside it */
                a1_is_not_in_alist = a1_is_star_atom = 0;
                a2_is_not_in_alist = a2_is_star_atom = 0;
                if (!is_in_the_ilist( u->alist, a1, u->na ))
                {
                    a1_is_not_in_alist = 1;
                }
                if (is_in_the_ilist( pd->pzz, a1, pd->n_pzz ))
                {
                    a1_is_star_atom = 1;
                }
                if (!is_in_the_ilist( u->alist, a2, u->na ))
                {
                    a2_is_not_in_alist = 1;
                }
                if (is_in_the_ilist( pd->pzz, a2, pd->n_pzz ))
                {
                    a2_is_star_atom = 1;
                }
                if (( a1_is_not_in_alist || a1_is_star_atom ) &&
                    ( a2_is_not_in_alist || a2_is_star_atom ))
                {
                    TREAT_ERR( err, 9032, "Caps of polymer unit lie inside it" );
                    goto exit_function;
                }
            }

            if (u->type == POLYMER_STY_SRU || u->type == POLYMER_STY_MOD ||
                 u->type == POLYMER_STY_CRO || u->type == POLYMER_STY_MER)
            {
                /* If SRU connection is missing, set to default ('either') */
                if (u->conn == POLYMER_CONN_NON)
                {
                    if (!bNoWarnings)
                    {
                        WarningMessage( pStrErr, "Set missing copolymer unit connection to EU" );
                    }
                    u->conn = POLYMER_CONN_EU;
                }

                if (u->cap1 && u->cap2)
                {
                    /* Set SRU closure type */
                    if (u->na == 1)
                    {
#ifdef ALLOW_CLOSING_SRU_VIA_DIRADICAL
                        u->cyclizable = CLOSING_SRU_DIRADICAL;
#else
                        u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
#ifdef  CLOSING_STARRED_SRU_IS_A_MUST
                        TREAT_ERR( err, 9029, "Could not perform SRU closure" );
                        goto exit_function;
#endif
#endif
                    }
                    else if (u->na == 2)
                    {

#ifdef ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND
                        u->cyclizable = CLOSING_SRU_HIGHER_ORDER_BOND;
#else
                        u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
#ifdef  CLOSING_STARRED_SRU_IS_A_MUST
                        TREAT_ERR( err, 9029, "Could not perform SRU closure" );
                        goto exit_function;
#endif
#endif
                    }
                    else
                    {
                        u->cyclizable = CLOSING_SRU_RING;
                    }
                }

                if (u->conn != POLYMER_CONN_HT)
                {
                    /* frame shift/SRU cyclization is for head-to-tail connections only */
                    u->cyclizable = CLOSING_SRU_NOT_APPLICABLE;
                }

                if (u->cyclizable != CLOSING_SRU_NOT_APPLICABLE)
                {
                    /* Allocate PS (frame-shiftable) bonds */
                    if (u->bkbonds)
                    {
                        imat_free(u->maxbkbonds, u->bkbonds);
                        u->bkbonds = NULL;
                    }
                    u->maxbkbonds = orig_at_data->num_inp_bonds + 2;
                    err = imat_new( u->maxbkbonds, 2, &( u->bkbonds ) );
                    if (err)
                    {
                        TREAT_ERR( err, 9034, "Not enough memory (polymers)" );
                        goto exit_function;
                    }
                }
            }

        }
    }
    else
    {
       TREAT_ERR( err, 9035, "Invalid kind of polymer representation" );
       goto exit_function;
    }

    orig_at_data->valid_polymer = 1;

exit_function:
    if (err)
    {
        orig_at_data->valid_polymer = 0;
    }

    return err;
}
    */
    // END INCHI C FUNCTION: OAD_ValidatePolymerAndPseudoElementData
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_ValidatePolymerAndPseudoElementData
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; ALLOW_MIXED_SRU_AND_MON=1.
    // INCHI✔️❌: ALLOW_CLOSING_SRU_VIA_DIRADICAL=1 and ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND=1.
    // INCHI✔️❌: The nested CLOSING_STARRED_SRU_IS_A_MUST error branches are therefore inactive.
    // INCHI✔️❌: WarningMessage expands to AddErrorMessage; TREAT_ERR preserves an existing nonzero error and appends its message.
    // INCHI✔️❌: Performance remains materially worse because OAD_Polymer_GetRepresentation clones SourceHeap-backed model values.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_ValidatePolymerAndPseudoElementData

    fn treat_error(
        error: &mut i32,
        code: i32,
        error_text: &mut Option<&mut [i8]>,
        message: &'static [u8],
    ) -> Result<(), SourceHeapError> {
        if *error == 0 && code != 0 {
            *error = code;
        }
        let message = message.iter().map(|byte| *byte as i8).collect::<Vec<_>>();
        AddErrorMessage(error_text.as_deref_mut(), Some(&message))?;
        Ok(())
    }

    fn warning(
        error_text: &mut Option<&mut [i8]>,
        no_warnings: i32,
        message: &'static [u8],
    ) -> Result<(), SourceHeapError> {
        if no_warnings == 0 {
            let message = message.iter().map(|byte| *byte as i8).collect::<Vec<_>>();
            AddErrorMessage(error_text.as_deref_mut(), Some(&message))?;
        }
        Ok(())
    }

    fn unit_pointer(
        heap: &SourceHeap,
        polymer: &OAD_Polymer,
        index: i32,
    ) -> Result<SourceMutPointer<OAD_PolymerUnit>, SourceHeapError> {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let pointer = *heap
            .slice(polymer.units.as_const())?
            .get(index)
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
        if pointer.is_null() {
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        }
        Ok(pointer)
    }

    fn atom_name_is(atom: &inp_ATOM, expected: &[u8]) -> Result<bool, SourceHeapError> {
        let length = atom
            .elname
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
        Ok(atom.elname[..length]
            .iter()
            .map(|byte| *byte as u8)
            .eq(expected.iter().copied()))
    }

    let mut error = 0_i32;
    let polymer_pointer = original_atom_data.polymer;
    let atom_count = original_atom_data.num_inp_atoms;
    original_atom_data.valid_polymer = 0;
    if treat_polymers != 0 && !polymer_pointer.is_null() {
        original_atom_data.valid_polymer = 1;
    }

    let execution = (|| -> Result<(), SourceHeapError> {
        let mut number_of_groups = 0_i32;
        let mut polymer = None;
        if original_atom_data.valid_polymer != 0 {
            let value = heap
                .slice(polymer_pointer.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            number_of_groups = value.n;
            polymer = Some(value);
        }

        if number_of_groups == 1 {
            let polymer = polymer.as_ref().ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            let pointer = unit_pointer(heap, polymer, 0)?;
            let unit = heap
                .slice(pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            if unit.type_ == POLYMER_STY_COP as i32 {
                treat_error(
                    &mut error,
                    9001,
                    &mut error_text,
                    b"Copolymer must contain more than one unit\0",
                )?;
                return Ok(());
            }
            if matches!(
                unit.subtype,
                value if value == POLYMER_SST_RAN as i32
                    || value == POLYMER_SST_ALT as i32
                    || value == POLYMER_SST_BLK as i32
            ) {
                treat_error(
                    &mut error,
                    9002,
                    &mut error_text,
                    b"Single polymer unit may not be RAN/ALT/BLO\0",
                )?;
                return Ok(());
            }
        }

        for unit_index in 0..number_of_groups {
            let polymer = polymer.as_ref().ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            let pointer = unit_pointer(heap, polymer, unit_index)?;
            let unit = heap
                .slice(pointer.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            if unit.nb != 0 && unit.nb != 2 {
                treat_error(
                    &mut error,
                    9003,
                    &mut error_text,
                    b"Number of crossing bonds in polymer unit is not 0 or 2\0",
                )?;
                return Ok(());
            }
            if unit.na < 1 {
                treat_error(&mut error, 9004, &mut error_text, b"Empty polymer unit\0")?;
                return Ok(());
            }
            if unit.na > atom_count {
                treat_error(&mut error, 9005, &mut error_text, b"Too large polymer unit\0")?;
                return Ok(());
            }
            let atom_list = heap.slice(unit.alist.as_const())?;
            for atom_index in 0..unit.na {
                let atom_index = usize::try_from(atom_index)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let atom = *atom_list
                    .get(atom_index)
                    .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                if atom < 1 || atom > atom_count {
                    treat_error(
                        &mut error,
                        9006,
                        &mut error_text,
                        b"Invalid atom number in polymer unit\0",
                    )?;
                    return Ok(());
                }
            }

            heap.with_slice_mut_and_heap(pointer, |units, heap| {
                let unit = units
                    .first_mut()
                    .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                OAD_PolymerUnit_SetEndsAndCaps(
                    heap,
                    unit,
                    original_atom_data,
                    &mut error,
                    error_text.as_deref_mut(),
                )
            })?;
            if error != 0 {
                return Ok(());
            }
            let unit = heap
                .slice_mut(pointer)?
                .first_mut()
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            unit.nbkbonds = 0;
            unit.cyclizable = CLOSING_SRU_NOT_APPLICABLE as i32;
            unit.cyclized = 0;
        }

        OAD_ValidateAndSortOutPseudoElementAtoms(
            heap,
            original_atom_data,
            treat_polymers,
            allow_non_polymer_zz,
            &mut error,
            error_text.as_deref_mut(),
        )?;
        if error != 0 {
            return Ok(());
        }

        if original_atom_data.n_zy > 0 && allow_non_polymer_zz == 0 {
            treat_error(
                &mut error,
                9,
                &mut error_text,
                b"Non-polymer-related Zz/star atoms are not allowed\0",
            )?;
            return Ok(());
        }
        if polymer_pointer.is_null() || original_atom_data.valid_polymer == 0 {
            return Ok(());
        }

        let polymer_value = heap
            .slice(polymer_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
        if polymer_value.n_pzz > 0 {
            if polymer_value.treat == POLYMERS_NO as i32 {
                treat_error(
                    &mut error,
                    9,
                    &mut error_text,
                    b"Pseudoelement endgroups are not allowed\0",
                )?;
                return Ok(());
            }
            if !polymer_value.pzz.is_null() {
                inchi_free(heap, polymer_value.pzz)?;
                heap.slice_mut(polymer_pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::UnsupportedSourceBehavior)?
                    .pzz = SourceMutPointer::null();
            }
            let count = u64::try_from(polymer_value.n_pzz)
                .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
            let pseudoatoms = match inchi_calloc::<i32>(heap, count, 4) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => {
                    treat_error(&mut error, 9010, &mut error_text, b"Not enough memory\0")?;
                    return Ok(());
                }
                Err(other) => return Err(other),
            };
            heap.slice_mut(polymer_pointer)?
                .first_mut()
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?
                .pzz = pseudoatoms;
            let mut output_index = 0_usize;
            for atom_index in 0..atom_count {
                let atom_index = usize::try_from(atom_index)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let is_zz = {
                    let atoms = heap.slice(original_atom_data.at.as_const())?;
                    let atom = atoms
                        .get(atom_index)
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    atom_name_is(atom, b"Zz")?
                };
                if is_zz {
                    *heap
                        .slice_mut(pseudoatoms)?
                        .get_mut(output_index)
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)? =
                        i32::try_from(atom_index)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                            .checked_add(1)
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    output_index = output_index
                        .checked_add(1)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                }
            }
        }

        let polymer_value = heap
            .slice(polymer_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
        for unit_index in 0..polymer_value.n {
            let pointer = unit_pointer(heap, &polymer_value, unit_index)?;
            let unit = heap
                .slice(pointer.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            if unit.type_ == POLYMER_STY_COP as i32 || unit.type_ == POLYMER_STY_SRU as i32 {
                if unit.nb > 0 {
                    if unit.type_ == POLYMER_STY_COP as i32 {
                        treat_error(
                            &mut error,
                            9026,
                            &mut error_text,
                            b"Polymer COP unit contains bracket-crossing bonds, not supported\0",
                        )?;
                        return Ok(());
                    }
                    continue;
                }
                let mut included_units = 0_i32;
                for candidate_index in 0..polymer_value.n {
                    let candidate_pointer = unit_pointer(heap, &polymer_value, candidate_index)?;
                    let candidate = heap
                        .slice(candidate_pointer.as_const())?
                        .first()
                        .cloned()
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    if candidate.type_ == POLYMER_STY_COP as i32 {
                        continue;
                    }
                    let candidate_list = if candidate.alist.is_null() {
                        None
                    } else {
                        Some(heap.slice(candidate.alist.as_const())?)
                    };
                    let embedding_list = if unit.alist.is_null() {
                        None
                    } else {
                        Some(heap.slice(unit.alist.as_const())?)
                    };
                    if is_ilist_inside(
                        candidate_list,
                        candidate.na,
                        embedding_list,
                        unit.na,
                    )? != 0 {
                        included_units = included_units.wrapping_add(1);
                        if included_units == 2 {
                            break;
                        }
                    }
                }
                if included_units > 1 {
                    if unit.type_ != POLYMER_STY_COP as i32 {
                        heap.slice_mut(pointer)?
                            .first_mut()
                            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?
                            .type_ = POLYMER_STY_COP as i32;
                        warning(
                            &mut error_text,
                            no_warnings,
                            b"Convert multiple-subunits unit to copolymer\0",
                        )?;
                    }
                } else if unit.type_ == POLYMER_STY_COP as i32 {
                    treat_error(
                        &mut error,
                        9027,
                        &mut error_text,
                        b"Polymer COP unit contains a single SRU instead of multiple\0",
                    )?;
                    return Ok(());
                }
            }
        }

        let representation = OAD_Polymer_GetRepresentation(heap, polymer_pointer)?;
        if representation == POLYMER_REPRESENTATION_SOURCE_BASED as i32 {
            for unit_index in 0..number_of_groups {
                let polymer_value = heap
                    .slice(polymer_pointer.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                let pointer = unit_pointer(heap, &polymer_value, unit_index)?;
                let unit = heap
                    .slice_mut(pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                if unit.type_ == POLYMER_STY_SRU as i32 {
                    unit.type_ = POLYMER_STY_MON as i32;
                    warning(
                        &mut error_text,
                        no_warnings,
                        b"Converted src-based polymer unit type to MON\0",
                    )?;
                }
                if unit.type_ == POLYMER_STY_COP as i32
                    && unit.subtype == POLYMER_SST_NON as i32
                {
                    unit.subtype = POLYMER_SST_RAN as i32;
                    warning(
                        &mut error_text,
                        no_warnings,
                        b"Set missing copolymer subtype to RAN\0",
                    )?;
                }
                if unit.conn != POLYMER_CONN_NON as i32 {
                    unit.conn = POLYMER_CONN_NON as i32;
                    warning(
                        &mut error_text,
                        no_warnings,
                        b"Ignore connection pattern for src-based polymer unit\0",
                    )?;
                }
            }
        } else if representation == POLYMER_REPRESENTATION_STRUCTURE_BASED as i32
            || representation == POLYMER_REPRESENTATION_MIXED as i32
        {
            for unit_index in 0..number_of_groups {
                let polymer_value = heap
                    .slice(polymer_pointer.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                let pointer = unit_pointer(heap, &polymer_value, unit_index)?;
                {
                    let unit = heap
                        .slice_mut(pointer)?
                        .first_mut()
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    if unit.nb == 0 && unit.type_ == POLYMER_STY_SRU as i32 {
                        unit.type_ = POLYMER_STY_COP as i32;
                        warning(
                            &mut error_text,
                            no_warnings,
                            b"Set copolymer embedding unit mark to COP\0",
                        )?;
                    }
                    if unit.type_ == POLYMER_STY_COP as i32 {
                        unit.cyclizable = CLOSING_SRU_NOT_APPLICABLE as i32;
                        if unit.subtype == POLYMER_SST_NON as i32 {
                            unit.subtype = POLYMER_SST_RAN as i32;
                            warning(
                                &mut error_text,
                                no_warnings,
                                b"Set missing copolymer subtype to RAN\0",
                            )?;
                        }
                        continue;
                    }
                    if unit.type_ == POLYMER_STY_MON as i32 {
                        continue;
                    }
                }

                let unit = heap
                    .slice(pointer.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                let bond_list = if unit.blist.is_null() {
                    None
                } else {
                    Some(heap.slice(unit.blist.as_const())?)
                };
                for bond_index in 0..unit.nb {
                    let offset = bond_index
                        .checked_mul(2)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    let offset = usize::try_from(offset)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    let bonds = bond_list.ok_or(SourceHeapError::NullPointer)?;
                    let atom1 = *bonds
                        .get(offset)
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    let atom2 = *bonds
                        .get(offset + 1)
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    let atom1_index = atom1
                        .checked_sub(1)
                        .and_then(|value| usize::try_from(value).ok())
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    let atom2_index = atom2
                        .checked_sub(1)
                        .and_then(|value| usize::try_from(value).ok())
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    let atoms = heap.slice(original_atom_data.at.as_const())?;
                    let first = atoms
                        .get(atom1_index)
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    if atom_name_is(first, b"H")?
                        || atom_name_is(first, b"D")?
                        || atom_name_is(first, b"T")?
                    {
                        treat_error(
                            &mut error,
                            9030,
                            &mut error_text,
                            b"H as polymer end group is not supported\0",
                        )?;
                        return Ok(());
                    }
                    let second = atoms
                        .get(atom2_index)
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    if atom_name_is(second, b"H")?
                        || atom_name_is(second, b"D")?
                        || atom_name_is(second, b"T")?
                    {
                        treat_error(
                            &mut error,
                            9031,
                            &mut error_text,
                            b"H as polymer end group is not supported\0",
                        )?;
                        return Ok(());
                    }
                    let atom_list = if unit.alist.is_null() {
                        None
                    } else {
                        Some(heap.slice(unit.alist.as_const())?)
                    };
                    let pseudo_list = if polymer_value.pzz.is_null() {
                        None
                    } else {
                        Some(heap.slice(polymer_value.pzz.as_const())?)
                    };
                    let atom1_outside =
                        is_in_the_ilist(atom_list, atom1, unit.na)?.is_none();
                    let atom1_star =
                        is_in_the_ilist(pseudo_list, atom1, polymer_value.n_pzz)?.is_some();
                    let atom2_outside =
                        is_in_the_ilist(atom_list, atom2, unit.na)?.is_none();
                    let atom2_star =
                        is_in_the_ilist(pseudo_list, atom2, polymer_value.n_pzz)?.is_some();
                    if (atom1_outside || atom1_star) && (atom2_outside || atom2_star) {
                        treat_error(
                            &mut error,
                            9032,
                            &mut error_text,
                            b"Caps of polymer unit lie inside it\0",
                        )?;
                        return Ok(());
                    }
                }

                let unit = heap
                    .slice_mut(pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                if matches!(
                    unit.type_,
                    value if value == POLYMER_STY_SRU as i32
                        || value == POLYMER_STY_MOD as i32
                        || value == POLYMER_STY_CRO as i32
                        || value == POLYMER_STY_MER as i32
                ) {
                    if unit.conn == POLYMER_CONN_NON as i32 {
                        warning(
                            &mut error_text,
                            no_warnings,
                            b"Set missing copolymer unit connection to EU\0",
                        )?;
                        unit.conn = POLYMER_CONN_EU as i32;
                    }
                    if unit.cap1 != 0 && unit.cap2 != 0 {
                        unit.cyclizable = if unit.na == 1 {
                            CLOSING_SRU_DIRADICAL as i32
                        } else if unit.na == 2 {
                            CLOSING_SRU_HIGHER_ORDER_BOND as i32
                        } else {
                            CLOSING_SRU_RING as i32
                        };
                    }
                    if unit.conn != POLYMER_CONN_HT as i32 {
                        unit.cyclizable = CLOSING_SRU_NOT_APPLICABLE as i32;
                    }
                }
                let (cyclizable, old_bonds, old_maximum) = {
                    let unit = heap
                        .slice(pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    (unit.cyclizable, unit.bkbonds, unit.maxbkbonds)
                };
                if cyclizable != CLOSING_SRU_NOT_APPLICABLE as i32 {
                    if !old_bonds.is_null() {
                        imat_free(heap, old_maximum, old_bonds)?;
                        heap.slice_mut(pointer)?
                            .first_mut()
                            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?
                            .bkbonds = SourceMutPointer::null();
                    }
                    let maximum = original_atom_data
                        .num_inp_bonds
                        .checked_add(2)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    heap.slice_mut(pointer)?
                        .first_mut()
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?
                        .maxbkbonds = maximum;
                    let mut bonds = heap
                        .slice(pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?
                        .bkbonds;
                    let allocation_error = imat_new(heap, maximum, 2, &mut bonds)?;
                    heap.slice_mut(pointer)?
                        .first_mut()
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?
                        .bkbonds = bonds;
                    if allocation_error != 0 {
                        error = allocation_error;
                        treat_error(
                            &mut error,
                            9034,
                            &mut error_text,
                            b"Not enough memory (polymers)\0",
                        )?;
                        return Ok(());
                    }
                }
            }
        } else {
            treat_error(
                &mut error,
                9035,
                &mut error_text,
                b"Invalid kind of polymer representation\0",
            )?;
            return Ok(());
        }

        original_atom_data.valid_polymer = 1;
        Ok(())
    })();
    execution?;
    if error != 0 {
        original_atom_data.valid_polymer = 0;
    }
    Ok(error)
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn OAD_ValidateAndSortOutPseudoElementAtoms(
    heap: &mut SourceHeap,
    original_atom_data: &mut ORIG_ATOM_DATA,
    treat_polymers: i32,
    allow_non_polymer_zz: i32,
    error: &mut i32,
    mut error_text: Option<&mut [i8]>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4396 OAD_ValidateAndSortOutPseudoElementAtoms
    // INCHI✔️✔️: complete source frame follows verbatim.
    /*
void OAD_ValidateAndSortOutPseudoElementAtoms( ORIG_ATOM_DATA *orig_at_data,
                                               int treat_polymers,
                                               int bNPZz,
                                               int *err,
                                               char *pStrErr )
{
    int i, k, n_pseudo = 0;
    int nsgroups = 0, nzz = 0;
    int nat = orig_at_data->num_inp_atoms;
    OAD_Polymer *pd = orig_at_data->polymer;
    OAD_PolymerUnit* u = NULL;

    int pseudos_allowed = (bNPZz == 1) || (treat_polymers != POLYMERS_NO);

    for (k = 0; k < nat; k++)
    {
        int is_zz = 0, is_star = 0, is_zy=0;


        /* Though "Zy" is present in our internal Periodic Table,
           _input_ "Zy" atoms are prohibited.
        */
        if (!strcmp( orig_at_data->at[k].elname, "Zy" ))
        {
            is_zy = 1;
#if 0
            Disabled 2020-04-07
            TREAT_ERR(*err, (70 + 5), "Invalid element(s):");
            TREAT_ERR(*err, (70 + 5), orig_at_data->at[k].elname);
            continue;
#endif 
        }
        is_star = !strcmp( orig_at_data->at[k].elname, "*" );
        if (!is_star)
        {
            is_zz = !strcmp( orig_at_data->at[k].elname, "Zz" );
        }

        if (is_star || is_zz || is_zy)
        {
            n_pseudo++;
            if (0==pseudos_allowed)
            {
                /* That's an error */
                TREAT_ERR( *err, ( 70 + 5 ), "Invalid element(s):" );
                TREAT_ERR( *err, ( 70 + 5 ), orig_at_data->at[k].elname );
                continue;
            }

            /* Now check if valid pseudoelement atom */

            /* Should be strictly univalent and single-bonded */
            /* Should not have isotopic enrichment */
            if ( orig_at_data->at[k].valence > 1          ||
                 orig_at_data->at[k].chem_bonds_valence>1 /* || orig_at_data->at[k].iso_atw_diff != 0    */
                                                           )
            {
                TREAT_ERR( *err, ( 70 + 7 ), "Invalid pseudo element(s) bonding" );
                /*TREAT_ERR( *err, ( 70 + 7 ), orig_at_data->at[k].elname );*/
                continue;
            }


            /* Now convert both "*" and "Zz" temporarily to "Zy" */
            mystrncpy( orig_at_data->at[k].elname, "Zy", sizeof( "Zy" ) );
        }
    }

    orig_at_data->n_zy = 0;
    nzz = 0;
    if (orig_at_data->valid_polymer)
    {
        nsgroups = pd->n;
        /* If applicable, check each CRU and back-convert "Zy" to "Zz" (polymer-related pseudoelement atoms)
           if they are from valid paired CRU crossing bond out-of-bracket caps */
        for (i = 0; i < nsgroups; i++)
        {
            u = pd->units[i];
            if (u)
            {
                if (u->cap1_is_undef && u->cap2_is_undef)
                {
                    /* valid pair: CRU is capped with two undefined-nature atoms, call them "Zz", finally */
                    strcpy( orig_at_data->at[u->cap1 - 1].elname, "Zz" );
                    strcpy( orig_at_data->at[u->cap2 - 1].elname, "Zz" );
                    nzz+= 2;
                }
            }
        }
        orig_at_data->polymer->n_pzz = nzz;
    }

    orig_at_data->n_zy = n_pseudo - nzz;
    if (orig_at_data->n_zy)
    {
        /* Have non-polymer-related pseudoelement atoms */
        if (0==bNPZz)
        {
            TREAT_ERR( *err, ( 70 + 4 ), "Polymer-unrelated pseudoatoms are not allowed" );
        }
    }

    if (*err)
    {
        orig_at_data->valid_polymer = 0;
    }

}
    */
    // END INCHI C FUNCTION: OAD_ValidateAndSortOutPseudoElementAtoms
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: OAD_ValidateAndSortOutPseudoElementAtoms
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; POLYMERS_NO=0.
    // INCHI✔️✔️: The #if 0 rejection of input "Zy" is inactive; Zy remains a counted pseudoatom.
    // INCHI✔️✔️: TREAT_ERR preserves a pre-existing nonzero error and always invokes AddErrorMessage.
    // INCHI✔️✔️: mystrncpy is the active util.c function; strcmp/strcpy are active libc behavior.
    // END INCHI ACTIVE MACRO CONFIGURATION: OAD_ValidateAndSortOutPseudoElementAtoms

    fn treat_error(
        error: &mut i32,
        code: i32,
        error_text: &mut Option<&mut [i8]>,
        message: &[i8],
    ) -> Result<(), SourceHeapError> {
        if *error == 0 && code != 0 {
            *error = code;
        }
        AddErrorMessage(error_text.as_deref_mut(), Some(message))?;
        Ok(())
    }

    let atom_count = original_atom_data.num_inp_atoms;
    let polymer_pointer = original_atom_data.polymer;
    let pseudos_allowed = allow_non_polymer_zz == 1 || treat_polymers != POLYMERS_NO as i32;
    let mut pseudo_count = 0_i32;

    if atom_count > 0 {
        let atom_count = usize::try_from(atom_count)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atoms = heap.slice_mut(original_atom_data.at)?;
        for atom_index in 0..atom_count {
            let atom = atoms
                .get_mut(atom_index)
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            let name_length = atom
                .elname
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            let name = &atom.elname[..name_length];
            let is_zy = name == [b'Z' as i8, b'y' as i8];
            let is_star = name == [b'*' as i8];
            let is_zz = !is_star && name == [b'Z' as i8, b'z' as i8];

            if is_star || is_zz || is_zy {
                pseudo_count = pseudo_count
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                if !pseudos_allowed {
                    treat_error(
                        error,
                        75,
                        &mut error_text,
                        &b"Invalid element(s):\0".map(|byte| byte as i8),
                    )?;
                    treat_error(error, 75, &mut error_text, &atom.elname)?;
                    continue;
                }
                if atom.valence > 1 || atom.chem_bonds_valence > 1 {
                    treat_error(
                        error,
                        77,
                        &mut error_text,
                        &b"Invalid pseudo element(s) bonding\0".map(|byte| byte as i8),
                    )?;
                    continue;
                }
                mystrncpy_slice(
                    Some(&mut atom.elname),
                    Some(&[b'Z' as i8, b'y' as i8, 0]),
                    3,
                )?;
            }
        }
    }

    original_atom_data.n_zy = 0;
    let mut zz_count = 0_i32;
    if original_atom_data.valid_polymer != 0 {
        if polymer_pointer.is_null() {
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        }
        let (number_of_groups, unit_pointers) = {
            let polymer = heap
                .slice(polymer_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
            (polymer.n, polymer.units)
        };
        if number_of_groups > 0 {
            let number_of_groups = usize::try_from(number_of_groups)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            for unit_index in 0..number_of_groups {
                let unit_pointer = *heap
                    .slice(unit_pointers.as_const())?
                    .get(unit_index)
                    .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                if unit_pointer.is_null() {
                    continue;
                }
                let (both_caps_undefined, cap1, cap2) = {
                    let unit = heap
                        .slice(unit_pointer.as_const())?
                        .first()
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    (
                        unit.cap1_is_undef != 0 && unit.cap2_is_undef != 0,
                        unit.cap1,
                        unit.cap2,
                    )
                };
                if both_caps_undefined {
                    let cap1_index = cap1
                        .checked_sub(1)
                        .and_then(|value| usize::try_from(value).ok())
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    let first_cap = heap
                        .slice_mut(original_atom_data.at)?
                        .get_mut(cap1_index)
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    first_cap.elname[..3]
                        .copy_from_slice(&[b'Z' as i8, b'z' as i8, 0]);

                    let cap2_index = cap2
                        .checked_sub(1)
                        .and_then(|value| usize::try_from(value).ok())
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    let second_cap = heap
                        .slice_mut(original_atom_data.at)?
                        .get_mut(cap2_index)
                        .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
                    second_cap.elname[..3]
                        .copy_from_slice(&[b'Z' as i8, b'z' as i8, 0]);
                    zz_count = zz_count
                        .checked_add(2)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                }
            }
        }
        heap.slice_mut(polymer_pointer)?
            .first_mut()
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?
            .n_pzz = zz_count;
    }

    original_atom_data.n_zy = pseudo_count
        .checked_sub(zz_count)
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    if original_atom_data.n_zy != 0 && allow_non_polymer_zz == 0 {
        treat_error(
            error,
            74,
            &mut error_text,
            &b"Polymer-unrelated pseudoatoms are not allowed\0".map(|byte| byte as i8),
        )?;
    }
    if *error != 0 {
        original_atom_data.valid_polymer = 0;
    }
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn OAD_PolymerUnit_DebugTrace(unit: Option<&OAD_PolymerUnit>) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3642 OAD_PolymerUnit_DebugTrace
    // INCHI✔❌: void OAD_PolymerUnit_DebugTrace( OAD_PolymerUnit *u )
    // INCHI✔❌: {
    // INCHI✔❌:     char *conn = "ABSENT", *typ = "ABSENT", *styp = "ABSENT";
    // INCHI✔❌:
    // INCHI✔❌:     if (!u)
    // INCHI✔❌:     {
    // INCHI✔❌:         return;
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (u->conn == 1)
    // INCHI✔❌:     {
    // INCHI✔❌:         conn = "HT"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->conn == 2)
    // INCHI✔❌:     {
    // INCHI✔❌:         conn = "HH"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->conn == 3)
    // INCHI✔❌:     {
    // INCHI✔❌:         conn = "EU"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (u->type == 0)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "NONE"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->type == 1)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "SRU"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->type == 2)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "MON"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->type == 3)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "COP"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->type == 4)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "MOD"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->type == 5)
    // INCHI✔❌:     {
    // INCHI✔❌:         typ = "MER"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     if (u->subtype == 1)
    // INCHI✔❌:     {
    // INCHI✔❌:         styp = "ALT"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->subtype == 2)
    // INCHI✔❌:     {
    // INCHI✔❌:         styp = "RAN"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:     else if (u->subtype == 3)
    // INCHI✔❌:     {
    // INCHI✔❌:         styp = "BLK"; /* djb-rwth: ignoring LLVM warning: possible presence of global variables */
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:     {
    // INCHI✔❌:         int i, k;
    // INCHI✔❌:         int na, nb;
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_("\n\nPOLYMER UNIT @ %-p", u);
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_( "\n\tid=%-d   label=%-d   type=%-s   subtype=%-s   conn=%-s   subscr='%-s'\n",
    // INCHI✔❌:                 u->id, u->label, typ, styp, conn, u->smt );
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_( "\tBracket1 coords: %-f, %-f, %-f, %-f\n", u->xbr1[0], u->xbr1[1], u->xbr1[2], u->xbr1[3] );
    // INCHI✔❌:         ITRACE_( "\tBracket2 coords: %-f, %-f, %-f, %-f\n", u->xbr2[0], u->xbr2[1], u->xbr2[2], u->xbr2[3] );
    // INCHI✔❌:
    // INCHI✔❌:         na = u->na;
    // INCHI✔❌:         ITRACE_( "\t%-d atoms { ", na );
    // INCHI✔❌:         for (k = 0; k < na - 1; k++)
    // INCHI✔❌:         {
    // INCHI✔❌:             ITRACE_( " %-d, ", u->alist[k] );
    // INCHI✔❌:         }
    // INCHI✔❌:         ITRACE_( " %-d }\n", u->alist[na - 1] );
    // INCHI✔❌:
    // INCHI✔❌:         nb = u->nb;
    // INCHI✔❌:         ITRACE_( "\t%-d bonds crossing unit borders { ", nb );
    // INCHI✔❌:
    // INCHI✔❌:         for (k = 0; k < nb; k++)
    // INCHI✔❌:         {
    // INCHI✔❌:             ITRACE_( " %-d-%-d ", u->blist[2 * k], u->blist[2 * k + 1] );
    // INCHI✔❌:         }
    // INCHI✔❌:         ITRACE_( "}\n" );
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_( "\tCRU caps and end atoms { " );
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_( "*%-d-[-%-d(end1) ... ", u->cap1, u->end_atom1 );
    // INCHI✔❌:         ITRACE_( "%-d(end2)-]-*%-d", u->end_atom2, u->cap2 );
    // INCHI✔❌:         ITRACE_( " }\n" );
    // INCHI✔❌:
    // INCHI✔❌:         ITRACE_( "\tBackbone bonds (may include 'artificially cyclizing' one) : %-d bonds ", u->nbkbonds );
    // INCHI✔❌:         if (u->nbkbonds)
    // INCHI✔❌:         {
    // INCHI✔❌:             ITRACE_(" { ");
    // INCHI✔❌:             for (i = 0; i < u->nbkbonds; i++)
    // INCHI✔❌:             {
    // INCHI✔❌:                 ITRACE_( "(%-d, %-d)  ", u->bkbonds[i][0], u->bkbonds[i][1] );
    // INCHI✔❌:             }
    // INCHI✔❌:             ITRACE_("}\n");
    // INCHI✔❌:         }
    // INCHI✔❌:     }
    // INCHI✔❌:
    // INCHI✔❌:
    // INCHI✔❌:     return;
    // INCHI✔❌: }
    // END INCHI C FUNCTION: OAD_PolymerUnit_DebugTrace
    // BEGIN INCHI ACTIVE HEADER MACRO: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.h:133 ITRACE_
    // INCHI✔❌: #define ITRACE_ 0 && _inchi_trace
    // END INCHI ACTIVE HEADER MACRO: ITRACE_

    let Some(unit) = unit else {
        return;
    };
    let connection = match unit.conn {
        1 => "HT",
        2 => "HH",
        3 => "EU",
        _ => "ABSENT",
    };
    let unit_type = match unit.type_ {
        0 => "NONE",
        1 => "SRU",
        2 => "MON",
        3 => "COP",
        4 => "MOD",
        5 => "MER",
        _ => "ABSENT",
    };
    let subtype = match unit.subtype {
        1 => "ALT",
        2 => "RAN",
        3 => "BLK",
        _ => "ABSENT",
    };
    let _ = (connection, unit_type, subtype);
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn Inp_Atom_GetBondType(
    heap: &SourceHeap,
    at: SourceConstPointer<inp_ATOM>,
    iatom1: i32,
    iatom2: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4507 Inp_Atom_GetBondType
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int Inp_Atom_GetBondType(inp_ATOM *at, int iatom1, int iatom2)
{
    int i;

    for (i = 0; i < at[iatom1].valence; i++)
    {
        if (at[iatom1].neighbor[i] == iatom2)
        {
                return at[iatom1].bond_type[i];
        }
    }

    return -1;
}
    */
    // END INCHI C FUNCTION: Inp_Atom_GetBondType
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: Inp_Atom_GetBondType
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no conditional source branches or active macros.
    // INCHI✔️❌: `S_CHAR valence`, `AT_NUMB neighbor`, and `U_CHAR bond_type` undergo C integer promotion before comparison or return.
    // INCHI✔️❌: Performance is materially worse because the modeled pointer access performs a SourceHeap allocation-table lookup.
    // END INCHI ACTIVE MACRO CONFIGURATION: Inp_Atom_GetBondType

    let atom_index = usize::try_from(iatom1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = heap
        .slice(at)?
        .get(atom_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.valence <= 0 {
        return Ok(-1);
    }
    let valence = usize::try_from(atom.valence)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if valence > atom.neighbor.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    for i in 0..valence {
        if i32::from(atom.neighbor[i]) == iatom2 {
            return Ok(i32::from(atom.bond_type[i]));
        }
    }
    Ok(-1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn preprocess_atom(element: u8, charge: i8) -> inp_ATOM {
        inp_ATOM {
            el_number: element,
            charge,
            ..inp_ATOM::default()
        }
    }

    fn preprocess_bond(atoms: &mut [inp_ATOM], first: usize, second: usize, bond_type: u8) {
        let first_ordinal = atoms[first].valence as usize;
        let second_ordinal = atoms[second].valence as usize;
        atoms[first].neighbor[first_ordinal] = second as AT_NUMB;
        atoms[first].bond_type[first_ordinal] = bond_type;
        atoms[first].valence += 1;
        atoms[first].chem_bonds_valence += bond_type as i8;
        atoms[second].neighbor[second_ordinal] = first as AT_NUMB;
        atoms[second].bond_type[second_ordinal] = bond_type;
        atoms[second].valence += 1;
        atoms[second].chem_bonds_valence += bond_type as i8;
    }

    fn preprocess_original(heap: &mut SourceHeap, atoms: Vec<inp_ATOM>) -> ORIG_ATOM_DATA {
        ORIG_ATOM_DATA {
            num_inp_atoms: atoms.len() as i32,
            at: heap.allocate_model_storage(atoms).unwrap(),
            ..ORIG_ATOM_DATA::default()
        }
    }

    #[test]
    fn source_port__runichi3__inp_atom_getbondtype__line_4507() {
        let mut heap = SourceHeap::default();

        let mut populated = inp_ATOM {
            valence: 5,
            ..inp_ATOM::default()
        };
        populated.neighbor[..5].copy_from_slice(&[0, AT_NUMB::MAX, 7, 7, 8]);
        populated.bond_type[..5].copy_from_slice(&[0, u8::MAX, 3, 4, 5]);

        let mut zero = inp_ATOM {
            valence: 0,
            ..inp_ATOM::default()
        };
        zero.neighbor[0] = 11;
        zero.bond_type[0] = 91;

        let mut negative = inp_ATOM {
            valence: -1,
            ..inp_ATOM::default()
        };
        negative.neighbor[0] = 12;
        negative.bond_type[0] = 92;

        let atoms = heap
            .allocate_model_storage(vec![populated, zero, negative])
            .unwrap();
        let before = heap.slice(atoms.as_const()).unwrap().to_vec();

        assert_eq!(Inp_Atom_GetBondType(&heap, atoms.as_const(), 0, 0), Ok(0));
        assert_eq!(
            Inp_Atom_GetBondType(&heap, atoms.as_const(), 0, i32::from(AT_NUMB::MAX)),
            Ok(i32::from(u8::MAX))
        );
        assert_eq!(Inp_Atom_GetBondType(&heap, atoms.as_const(), 0, 7), Ok(3));
        assert_eq!(Inp_Atom_GetBondType(&heap, atoms.as_const(), 0, 8), Ok(5));
        assert_eq!(Inp_Atom_GetBondType(&heap, atoms.as_const(), 0, 9), Ok(-1));
        assert_eq!(Inp_Atom_GetBondType(&heap, atoms.as_const(), 0, -1), Ok(-1));
        assert_eq!(
            Inp_Atom_GetBondType(&heap, atoms.as_const(), 0, i32::MIN),
            Ok(-1)
        );
        assert_eq!(
            Inp_Atom_GetBondType(&heap, atoms.as_const(), 0, i32::MAX),
            Ok(-1)
        );
        assert_eq!(Inp_Atom_GetBondType(&heap, atoms.as_const(), 1, 11), Ok(-1));
        assert_eq!(Inp_Atom_GetBondType(&heap, atoms.as_const(), 2, 12), Ok(-1));
        assert_eq!(heap.slice(atoms.as_const()).unwrap(), before);

        assert_eq!(
            Inp_Atom_GetBondType(&heap, atoms.as_const().offset(1).unwrap(), 0, 11),
            Ok(-1)
        );

        let mut full = inp_ATOM {
            valence: MAXVAL as i8,
            ..inp_ATOM::default()
        };
        for index in 0..MAXVAL as usize {
            full.neighbor[index] = index as AT_NUMB;
            full.bond_type[index] = (index as u8).wrapping_add(101);
        }
        let full = heap.allocate_model_storage(vec![full]).unwrap();
        assert_eq!(
            Inp_Atom_GetBondType(
                &heap,
                full.as_const(),
                0,
                i32::try_from(MAXVAL - 1).unwrap(),
            ),
            Ok((MAXVAL as i32 - 1).wrapping_add(101))
        );
        assert_eq!(Inp_Atom_GetBondType(&heap, full.as_const(), 0, 20), Ok(-1));

        let overfull = heap
            .allocate_model_storage(vec![inp_ATOM {
                valence: (MAXVAL + 1) as i8,
                ..inp_ATOM::default()
            }])
            .unwrap();
        assert_eq!(
            Inp_Atom_GetBondType(&heap, overfull.as_const(), 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            Inp_Atom_GetBondType(&heap, atoms.as_const(), -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            Inp_Atom_GetBondType(&heap, atoms.as_const(), 3, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            Inp_Atom_GetBondType(&heap, SourceConstPointer::null(), 0, 0),
            Err(SourceHeapError::NullPointer)
        );
    }

    fn preprocess_error_text(data: &STRUCT_DATA) -> Vec<u8> {
        let length = data
            .pStrErrStruct
            .iter()
            .position(|byte| *byte == 0)
            .unwrap();
        data.pStrErrStruct[..length]
            .iter()
            .map(|byte| *byte as u8)
            .collect()
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_findendsandcaps__line_2166() {
        fn atom(name: &[u8]) -> inp_ATOM {
            let mut atom = inp_ATOM::default();
            for (target, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
                *target = source as i8;
            }
            atom
        }

        fn unit(
            heap: &mut SourceHeap,
            atom_list: Vec<i32>,
            bonds: Option<Vec<i32>>,
            number_of_bonds: i32,
        ) -> OAD_PolymerUnit {
            OAD_PolymerUnit {
                na: atom_list.len() as i32,
                alist: heap.allocate_model_storage(atom_list).unwrap(),
                blist: bonds
                    .map(|values| heap.allocate_model_storage(values).unwrap())
                    .unwrap_or_default(),
                nb: number_of_bonds,
                end_atom1: 71,
                cap1: 72,
                end_atom2: 73,
                cap2: 74,
                ..OAD_PolymerUnit::default()
            }
        }

        fn invoke(
            heap: &SourceHeap,
            unit: &mut OAD_PolymerUnit,
            original: &ORIG_ATOM_DATA,
            error_text: &mut [i8; 256],
        ) -> (Result<(), SourceHeapError>, [i32; 7]) {
            let mut end1 = 81;
            let mut cap1 = 82;
            let mut cap1_star = 83;
            let mut end2 = 84;
            let mut cap2 = 85;
            let mut cap2_star = 86;
            let mut error = 87;
            let result = OAD_PolymerUnit_FindEndsAndCaps(
                heap,
                unit,
                original,
                &mut end1,
                &mut cap1,
                &mut cap1_star,
                &mut end2,
                &mut cap2,
                &mut cap2_star,
                &mut error,
                Some(error_text),
            );
            (
                result,
                [end1, cap1, cap1_star, end2, cap2, cap2_star, error],
            )
        }

        fn error_bytes(buffer: &[i8; 256]) -> Vec<u8> {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            buffer[..end].iter().map(|byte| *byte as u8).collect()
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                atom(b"Zz\0"),
                atom(b"C\0"),
                atom(b"N\0"),
                atom(b"O\0"),
                atom(b"Zz\0"),
                atom(b"H\0"),
            ])
            .unwrap();
        let original = ORIG_ATOM_DATA {
            at: atoms,
            num_inp_atoms: 6,
            ..ORIG_ATOM_DATA::default()
        };

        for mut early in [
            unit(&mut heap, vec![2, 3], None, 1),
            unit(&mut heap, vec![2, 3], Some(vec![1, 2, 3, 4]), 0),
        ] {
            let before = early.clone();
            let mut message = [0_i8; 256];
            let (result, outputs) = invoke(&heap, &mut early, &original, &mut message);
            assert_eq!(result, Ok(()));
            assert_eq!(outputs, [0; 7]);
            assert_eq!(early, before);
            assert!(error_bytes(&message).is_empty());
        }

        for bonds in [vec![1, 2, 3, 4], vec![2, 1, 4, 3]] {
            let mut success = unit(&mut heap, vec![2, 3], Some(bonds), 2);
            let mut message = [0_i8; 256];
            let (result, outputs) = invoke(&heap, &mut success, &original, &mut message);
            assert_eq!(result, Ok(()));
            assert_eq!(outputs, [2, 1, 1, 3, 4, 0, 0]);
            assert_eq!(
                (
                    success.end_atom1,
                    success.cap1,
                    success.end_atom2,
                    success.cap2,
                ),
                (2, 1, 3, 4)
            );
            assert!(error_bytes(&message).is_empty());
        }

        let mut left_inside = unit(&mut heap, vec![2, 3], Some(vec![2, 3, 1, 4]), 2);
        let left_before = left_inside.clone();
        let mut message = [0_i8; 256];
        let (result, outputs) = invoke(&heap, &mut left_inside, &original, &mut message);
        assert_eq!(result, Ok(()));
        assert_eq!(outputs, [0, 0, 0, 0, 0, 0, 9032]);
        assert_eq!(left_inside, left_before);
        assert_eq!(error_bytes(&message), b"Polymer CRU cap(s) lie inside CRU");

        let mut right_inside = unit(&mut heap, vec![2, 3], Some(vec![1, 2, 2, 3]), 2);
        let mut message = [0_i8; 256];
        let (result, outputs) = invoke(&heap, &mut right_inside, &original, &mut message);
        assert_eq!(result, Ok(()));
        assert_eq!(outputs, [2, 1, 1, 2, 3, 0, 0]);
        assert_eq!(
            (
                right_inside.end_atom1,
                right_inside.cap1,
                right_inside.end_atom2,
                right_inside.cap2,
            ),
            (2, 1, 2, 3)
        );
        assert_eq!(error_bytes(&message), b"Polymer CRU cap(s) lie inside CRU");

        for (bonds, expected_error, expected_message) in [
            (
                vec![1, 9, 3, 4],
                9090,
                b"Invalid polymer CRU crossing bond".as_slice(),
            ),
            (
                vec![1, 2, 4, 9],
                9091,
                b"Invalid polymer CRU crossing bond".as_slice(),
            ),
            (
                vec![1, 2, 1, 3],
                9090,
                b"Invalid polymer CRU surrounding".as_slice(),
            ),
        ] {
            let mut invalid = unit(&mut heap, vec![2, 3], Some(bonds), 2);
            let before = invalid.clone();
            let mut message = [0_i8; 256];
            let (result, outputs) = invoke(&heap, &mut invalid, &original, &mut message);
            assert_eq!(result, Ok(()));
            assert_eq!(outputs[6], expected_error);
            assert_eq!(invalid, before);
            assert_eq!(error_bytes(&message), expected_message);
        }

        let mut undefined_cap = unit(&mut heap, vec![2, 3], Some(vec![0, 2, 3, 4]), 2);
        let undefined_before = undefined_cap.clone();
        let mut message = [0_i8; 256];
        let (result, outputs) = invoke(&heap, &mut undefined_cap, &original, &mut message);
        assert_eq!(result, Err(SourceHeapError::UnsupportedSourceBehavior));
        assert_eq!(&outputs[..3], &[2, 0, 0]);
        assert_eq!(outputs[6], 0);
        assert_eq!(undefined_cap, undefined_before);
        assert!(error_bytes(&message).is_empty());

        let malformed_atoms = heap
            .allocate_model_storage(vec![inp_ATOM {
                elname: [b'Z' as i8; 6],
                ..inp_ATOM::default()
            }])
            .unwrap();
        let malformed_original = ORIG_ATOM_DATA {
            at: malformed_atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        let mut malformed = unit(&mut heap, Vec::new(), Some(vec![1, 1, 1, 1]), 2);
        let mut message = [0_i8; 256];
        assert_eq!(
            invoke(&heap, &mut malformed, &malformed_original, &mut message,).0,
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_setendsandcaps__line_2266() {
        fn atom(name: &[u8]) -> inp_ATOM {
            let mut atom = inp_ATOM::default();
            for (target, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
                *target = source as i8;
            }
            atom
        }

        fn unit(
            heap: &mut SourceHeap,
            atom_list: Vec<i32>,
            crossing_bonds: Option<Vec<i32>>,
        ) -> OAD_PolymerUnit {
            OAD_PolymerUnit {
                na: atom_list.len() as i32,
                alist: heap.allocate_model_storage(atom_list).unwrap(),
                nb: i32::from(crossing_bonds.is_some()) * 2,
                blist: crossing_bonds
                    .map(|values| heap.allocate_model_storage(values).unwrap())
                    .unwrap_or_default(),
                cyclizable: 91,
                end_atom1: 92,
                cap1: 93,
                end_atom2: 94,
                cap2: 95,
                cap1_is_undef: 96,
                cap2_is_undef: 97,
                ..OAD_PolymerUnit::default()
            }
        }

        fn original(heap: &mut SourceHeap, atoms: Vec<inp_ATOM>) -> ORIG_ATOM_DATA {
            ORIG_ATOM_DATA {
                num_inp_atoms: atoms.len() as i32,
                at: heap.allocate_model_storage(atoms).unwrap(),
                ..ORIG_ATOM_DATA::default()
            }
        }

        let mut heap = SourceHeap::default();
        let ordinary = original(
            &mut heap,
            vec![atom(b"C\0"), atom(b"N\0"), atom(b"O\0"), atom(b"H\0")],
        );
        let mut empty = unit(&mut heap, vec![2, 3], None);
        let mut error = 88;
        let mut message = [0_i8; 256];
        assert_eq!(
            OAD_PolymerUnit_SetEndsAndCaps(
                &heap,
                &mut empty,
                &ordinary,
                &mut error,
                Some(&mut message),
            ),
            Ok(())
        );
        assert_eq!(error, 0);
        assert_eq!(
            (
                empty.cyclizable,
                empty.end_atom1,
                empty.cap1,
                empty.end_atom2,
                empty.cap2,
                empty.cap1_is_undef,
                empty.cap2_is_undef,
            ),
            (CLOSING_SRU_NOT_APPLICABLE as i32, 0, 0, 0, 0, 0, 0)
        );

        let mut non_star = unit(&mut heap, vec![2, 3], Some(vec![1, 2, 3, 4]));
        assert_eq!(
            OAD_PolymerUnit_SetEndsAndCaps(
                &heap,
                &mut non_star,
                &ordinary,
                &mut error,
                Some(&mut message),
            ),
            Ok(())
        );
        assert_eq!(non_star.cyclizable, CLOSING_SRU_NOT_APPLICABLE as i32);
        assert_eq!(
            (
                non_star.end_atom1,
                non_star.cap1,
                non_star.end_atom2,
                non_star.cap2,
            ),
            (2, 1, 3, 4)
        );

        let stars = original(&mut heap, vec![atom(b"Zz\0"), atom(b"C\0"), atom(b"Zz\0")]);
        let mut diradical = unit(&mut heap, vec![2], Some(vec![1, 2, 2, 3]));
        assert_eq!(
            OAD_PolymerUnit_SetEndsAndCaps(
                &heap,
                &mut diradical,
                &stars,
                &mut error,
                Some(&mut message),
            ),
            Ok(())
        );
        assert_eq!(diradical.cyclizable, CLOSING_SRU_DIRADICAL as i32);
        assert_eq!((diradical.cap1_is_undef, diradical.cap2_is_undef), (1, 1));

        let mut adjacent_atoms = vec![atom(b"Zz\0"), atom(b"C\0"), atom(b"N\0"), atom(b"Zz\0")];
        adjacent_atoms[1].valence = 1;
        adjacent_atoms[1].neighbor[0] = 2;
        adjacent_atoms[1].bond_type[0] = 1;
        let adjacent = original(&mut heap, adjacent_atoms);
        let mut higher_order = unit(&mut heap, vec![2, 3], Some(vec![1, 2, 3, 4]));
        assert_eq!(
            OAD_PolymerUnit_SetEndsAndCaps(
                &heap,
                &mut higher_order,
                &adjacent,
                &mut error,
                Some(&mut message),
            ),
            Ok(())
        );
        assert_eq!(
            higher_order.cyclizable,
            CLOSING_SRU_HIGHER_ORDER_BOND as i32
        );

        let ring_atoms = original(
            &mut heap,
            vec![atom(b"Zz\0"), atom(b"C\0"), atom(b"N\0"), atom(b"Zz\0")],
        );
        let mut ring = unit(&mut heap, vec![2, 3], Some(vec![1, 2, 3, 4]));
        assert_eq!(
            OAD_PolymerUnit_SetEndsAndCaps(
                &heap,
                &mut ring,
                &ring_atoms,
                &mut error,
                Some(&mut message),
            ),
            Ok(())
        );
        assert_eq!(ring.cyclizable, CLOSING_SRU_RING as i32);

        let one_star_atoms = original(
            &mut heap,
            vec![atom(b"Zz\0"), atom(b"C\0"), atom(b"N\0"), atom(b"O\0")],
        );
        let mut one_star = unit(&mut heap, vec![2, 3], Some(vec![1, 2, 3, 4]));
        assert_eq!(
            OAD_PolymerUnit_SetEndsAndCaps(
                &heap,
                &mut one_star,
                &one_star_atoms,
                &mut error,
                Some(&mut message),
            ),
            Ok(())
        );
        assert_eq!(one_star.cyclizable, CLOSING_SRU_RING as i32);

        let mut invalid = unit(&mut heap, vec![2, 3], Some(vec![2, 3, 1, 4]));
        message.fill(0);
        assert_eq!(
            OAD_PolymerUnit_SetEndsAndCaps(
                &heap,
                &mut invalid,
                &ordinary,
                &mut error,
                Some(&mut message),
            ),
            Ok(())
        );
        assert_eq!(error, 9032);
        assert_eq!(
            (
                invalid.cyclizable,
                invalid.end_atom1,
                invalid.cap1,
                invalid.end_atom2,
                invalid.cap2,
            ),
            (CLOSING_SRU_NOT_APPLICABLE as i32, 0, 0, 0, 0)
        );

        let mut undefined = unit(&mut heap, vec![2, 3], Some(vec![0, 2, 3, 4]));
        assert_eq!(
            OAD_PolymerUnit_SetEndsAndCaps(
                &heap,
                &mut undefined,
                &ordinary,
                &mut error,
                Some(&mut message),
            ),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(
            (
                undefined.cyclizable,
                undefined.end_atom1,
                undefined.cap1,
                undefined.cap1_is_undef,
            ),
            (CLOSING_SRU_NOT_APPLICABLE as i32, 2, 0, 0)
        );
    }

    #[test]
    fn source_port__runichi3__origatdata_debugtrace__line_1108() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            OrigAtData_DebugTrace(&heap, &ORIG_ATOM_DATA::default()),
            Ok(())
        );
        assert_eq!(
            OrigAtData_DebugTrace(
                &heap,
                &ORIG_ATOM_DATA {
                    num_inp_atoms: -1,
                    at: SourceMutPointer::null(),
                    ..ORIG_ATOM_DATA::default()
                },
            ),
            Ok(())
        );

        let dangling_v3000 = heap
            .allocate_model_storage(vec![OAD_V3000 {
                n_star_atoms: i32::MAX,
                n_haptic_bonds: i32::MIN,
                n_collections: -7,
                ..OAD_V3000::default()
            }])
            .unwrap();
        inchi_free(&mut heap, dangling_v3000).unwrap();
        assert_eq!(
            OrigAtData_DebugTrace(
                &heap,
                &ORIG_ATOM_DATA {
                    v3000: dangling_v3000,
                    ..ORIG_ATOM_DATA::default()
                },
            ),
            Ok(())
        );

        let atom = inp_ATOM {
            valence: 2,
            neighbor: {
                let mut values = [0_u16; 20];
                values[0] = u16::MAX;
                values[1] = 19;
                values
            },
            charge: i8::MIN,
            radical: i8::MAX,
            x: f64::from_bits(0x7ff8_0000_0000_0042),
            y: f64::INFINITY,
            z: -0.0,
            ..inp_ATOM::default()
        };
        let atoms = heap.allocate_model_storage(vec![atom.clone()]).unwrap();
        let data = ORIG_ATOM_DATA {
            at: atoms,
            num_inp_atoms: 1,
            v3000: dangling_v3000,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(OrigAtData_DebugTrace(&heap, &data), Ok(()));
        let after = &heap.slice(atoms.as_const()).unwrap()[0];
        let mut after_without_coordinates = after.clone();
        after_without_coordinates.x = 0.0;
        after_without_coordinates.y = 0.0;
        after_without_coordinates.z = 0.0;
        let mut expected_without_coordinates = atom.clone();
        expected_without_coordinates.x = 0.0;
        expected_without_coordinates.y = 0.0;
        expected_without_coordinates.z = 0.0;
        assert_eq!(after_without_coordinates, expected_without_coordinates);
        assert_eq!(after.x.to_bits(), 0x7ff8_0000_0000_0042);
        assert_eq!(after.y.to_bits(), f64::INFINITY.to_bits());
        assert_eq!(after.z.to_bits(), (-0.0_f64).to_bits());

        assert_eq!(
            OrigAtData_DebugTrace(
                &heap,
                &ORIG_ATOM_DATA {
                    at: atoms,
                    num_inp_atoms: 2,
                    ..ORIG_ATOM_DATA::default()
                },
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let oversized = heap
            .allocate_model_storage(vec![inp_ATOM {
                valence: 21,
                ..inp_ATOM::default()
            }])
            .unwrap();
        assert_eq!(
            OrigAtData_DebugTrace(
                &heap,
                &ORIG_ATOM_DATA {
                    at: oversized,
                    num_inp_atoms: 1,
                    ..ORIG_ATOM_DATA::default()
                },
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let negative_valence = heap
            .allocate_model_storage(vec![inp_ATOM {
                valence: -1,
                ..inp_ATOM::default()
            }])
            .unwrap();
        assert_eq!(
            OrigAtData_DebugTrace(
                &heap,
                &ORIG_ATOM_DATA {
                    at: negative_valence,
                    num_inp_atoms: 1,
                    ..ORIG_ATOM_DATA::default()
                },
            ),
            Ok(())
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymer_debugtrace__line_3757() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            OAD_Polymer_DebugTrace(&heap, SourceMutPointer::null()),
            Ok(())
        );

        let dangling_pzz = heap.allocate_model_storage(vec![71_i32]).unwrap();
        inchi_free(&mut heap, dangling_pzz).unwrap();
        let no_units = heap
            .allocate_model_storage(vec![OAD_Polymer {
                n: 0,
                n_pzz: i32::MAX,
                pzz: dangling_pzz,
                units: SourceMutPointer::null(),
                ..OAD_Polymer::default()
            }])
            .unwrap();
        assert_eq!(OAD_Polymer_DebugTrace(&heap, no_units), Ok(()));

        let dangling_alist = heap.allocate_model_storage(vec![81_i32]).unwrap();
        let dangling_blist = heap.allocate_model_storage(vec![82_i32]).unwrap();
        let dangling_bond_row = heap.allocate_model_storage(vec![83_i32]).unwrap();
        let dangling_bonds = heap
            .allocate_model_storage(vec![dangling_bond_row])
            .unwrap();
        inchi_free(&mut heap, dangling_alist).unwrap();
        inchi_free(&mut heap, dangling_blist).unwrap();
        inchi_free(&mut heap, dangling_bond_row).unwrap();
        inchi_free(&mut heap, dangling_bonds).unwrap();
        let unit_value = OAD_PolymerUnit {
            conn: 3,
            type_: 5,
            subtype: 2,
            na: 3,
            nb: 2,
            nbkbonds: 1,
            alist: dangling_alist,
            blist: dangling_blist,
            bkbonds: dangling_bonds,
            ..OAD_PolymerUnit::default()
        };
        let unit = heap
            .allocate_model_storage(vec![unit_value.clone()])
            .unwrap();
        let units = heap
            .allocate_model_storage(vec![SourceMutPointer::null(), unit])
            .unwrap();
        let polymer_value = OAD_Polymer {
            n: 2,
            units,
            n_pzz: i32::MAX,
            pzz: dangling_pzz,
            really_do_frame_shift: i32::MIN,
            frame_shift_scheme: i32::MAX,
            edit_repeats: -7,
            ..OAD_Polymer::default()
        };
        let polymer = heap
            .allocate_model_storage(vec![polymer_value.clone()])
            .unwrap();
        assert_eq!(OAD_Polymer_DebugTrace(&heap, polymer), Ok(()));
        assert_eq!(heap.slice(polymer.as_const()).unwrap()[0], polymer_value);
        assert_eq!(heap.slice(unit.as_const()).unwrap()[0], unit_value);

        let dangling_units = heap
            .allocate_model_storage(vec![SourceMutPointer::<OAD_PolymerUnit>::null()])
            .unwrap();
        inchi_free(&mut heap, dangling_units).unwrap();
        let negative = heap
            .allocate_model_storage(vec![OAD_Polymer {
                n: -1,
                units: dangling_units,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        assert_eq!(OAD_Polymer_DebugTrace(&heap, negative), Ok(()));

        let short_units = heap.allocate_model_storage(vec![unit]).unwrap();
        let short = heap
            .allocate_model_storage(vec![OAD_Polymer {
                n: 2,
                units: short_units,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_Polymer_DebugTrace(&heap, short),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let dangling_unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit::default()])
            .unwrap();
        inchi_free(&mut heap, dangling_unit).unwrap();
        let dangling_unit_list = heap.allocate_model_storage(vec![dangling_unit]).unwrap();
        let dangling = heap
            .allocate_model_storage(vec![OAD_Polymer {
                n: 1,
                units: dangling_unit_list,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_Polymer_DebugTrace(&heap, dangling),
            Err(SourceHeapError::MissingAllocation)
        );

        let overflow_unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: i32::MIN,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let overflow_units = heap.allocate_model_storage(vec![overflow_unit]).unwrap();
        let overflow = heap
            .allocate_model_storage(vec![OAD_Polymer {
                n: 1,
                units: overflow_units,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_Polymer_DebugTrace(&heap, overflow),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
    }

    #[test]
    fn source_port__runichi3__preprocessonestructure__line_488() {
        let mut missing_heap = SourceHeap::default();
        let mut missing_data = STRUCT_DATA::default();
        let mut missing_original = ORIG_ATOM_DATA::default();
        assert_eq!(
            PreprocessOneStructure(
                &mut missing_heap,
                None,
                &mut missing_data,
                &INPUT_PARMS::default(),
                &mut missing_original,
                &mut [],
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut failure_heap = SourceHeap::default();
        let mut failure_original =
            preprocess_original(&mut failure_heap, vec![preprocess_atom(6, 0)]);
        let mut failure_prepared = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut failure_data = STRUCT_DATA::default();
        failure_data.bTautFlags[INCHI_BAS as usize] = 0x1000_0000;
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            PreprocessOneStructure(
                &mut failure_heap,
                None,
                &mut failure_data,
                &INPUT_PARMS::default(),
                &mut failure_original,
                &mut failure_prepared,
            ),
            Ok(_IS_FATAL as i32)
        );
        assert_eq!(failure_data.nStructReadError, 99);
        assert_eq!(preprocess_error_text(&failure_data), b"Out of RAM");
        assert!(failure_prepared[0].at.is_null());
        assert_eq!(failure_data.bTautFlags[INCHI_BAS as usize], 0x1000_0000);

        let mut component_failure_heap = SourceHeap::default();
        let mut component_failure_original =
            preprocess_original(&mut component_failure_heap, vec![preprocess_atom(6, 0)]);
        let mut component_failure_prepared = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut component_failure_data = STRUCT_DATA::default();
        // OrigAtData_Duplicate performs three source allocations; the next
        // allocation is MarkDisconnectedComponents' first work buffer.
        component_failure_heap.fail_after_allocations(3);
        assert_eq!(
            PreprocessOneStructure(
                &mut component_failure_heap,
                None,
                &mut component_failure_data,
                &INPUT_PARMS::default(),
                &mut component_failure_original,
                &mut component_failure_prepared,
            ),
            Ok(_IS_FATAL as i32)
        );
        assert_eq!(component_failure_data.nStructReadError, 99);
        assert_eq!(
            preprocess_error_text(&component_failure_data),
            b"Out of RAM"
        );
        assert!(!component_failure_prepared[0].at.is_null());
        assert_eq!(component_failure_prepared[0].num_components, -1);

        let mut empty_heap = SourceHeap::default();
        let mut empty_original = ORIG_ATOM_DATA::default();
        let mut empty_prepared = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut empty_data = STRUCT_DATA::default();
        assert_eq!(
            PreprocessOneStructure(
                &mut empty_heap,
                Some(&mut INCHI_CLOCK::default()),
                &mut empty_data,
                &INPUT_PARMS::default(),
                &mut empty_original,
                &mut empty_prepared,
            ),
            Ok(0)
        );
        assert_eq!(empty_prepared[0].num_inp_atoms, 0);
        assert_eq!(empty_prepared[0].num_components, 0);
        assert!(preprocess_error_text(&empty_data).is_empty());

        let mut parity_heap = SourceHeap::default();
        let mut parity_atom = preprocess_atom(6, 0);
        parity_atom.sb_parity = [0x39, 0, 0];
        let mut parity_original = preprocess_original(&mut parity_heap, vec![parity_atom]);
        let mut parity_prepared = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut parity_data = STRUCT_DATA::default();
        assert_eq!(
            PreprocessOneStructure(
                &mut parity_heap,
                None,
                &mut parity_data,
                &INPUT_PARMS::default(),
                &mut parity_original,
                &mut parity_prepared,
            ),
            Ok(0)
        );
        assert_eq!(
            parity_heap.slice(parity_prepared[0].at.as_const()).unwrap()[0].sb_parity,
            [1, 0, 0]
        );
        assert_eq!(parity_prepared[0].num_components, 1);

        let mut isotope_heap = SourceHeap::default();
        let mut isotope_atom = preprocess_atom(8, 0);
        isotope_atom.chem_bonds_valence = 1;
        isotope_atom.num_iso_H[0] = 1;
        let mut isotope_original = preprocess_original(&mut isotope_heap, vec![isotope_atom]);
        let mut isotope_prepared = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut isotope_data = STRUCT_DATA::default();
        let isotope_parameters = INPUT_PARMS {
            bNoWarnings: 1,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            PreprocessOneStructure(
                &mut isotope_heap,
                None,
                &mut isotope_data,
                &isotope_parameters,
                &mut isotope_original,
                &mut isotope_prepared,
            ),
            Ok(0)
        );
        assert_ne!(
            isotope_data.bTautFlagsDone[INCHI_BAS as usize]
                & TG_FLAG_FOUND_ISOTOPIC_H_DONE as INCHI_MODE,
            0
        );
        assert_ne!(
            isotope_data.bTautFlagsDone[INCHI_BAS as usize]
                & TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE as INCHI_MODE,
            0
        );

        let mut oxo_heap = SourceHeap::default();
        let mut oxo_atoms = vec![preprocess_atom(17, -1), preprocess_atom(8, 0)];
        preprocess_bond(&mut oxo_atoms, 0, 1, 2);
        let mut oxo_original = preprocess_original(&mut oxo_heap, oxo_atoms);
        let mut oxo_prepared = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut oxo_data = STRUCT_DATA::default();
        let oxo_parameters = INPUT_PARMS {
            bFixNonUniformDraw: 1,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            PreprocessOneStructure(
                &mut oxo_heap,
                None,
                &mut oxo_data,
                &oxo_parameters,
                &mut oxo_original,
                &mut oxo_prepared,
            ),
            Ok(_IS_WARNING as i32)
        );
        assert_eq!(preprocess_error_text(&oxo_data), b"Charges were rearranged");
        assert_ne!(
            oxo_data.bTautFlagsDone[INCHI_BAS as usize] & TG_FLAG_FIX_ODD_THINGS_DONE as INCHI_MODE,
            0
        );
        let oxo_output = oxo_heap.slice(oxo_prepared[0].at.as_const()).unwrap();
        assert_eq!((oxo_output[0].charge, oxo_output[1].charge), (0, -1));

        let mut quiet_heap = SourceHeap::default();
        let mut quiet_atoms = vec![preprocess_atom(17, -1), preprocess_atom(8, 0)];
        preprocess_bond(&mut quiet_atoms, 0, 1, 2);
        let mut quiet_original = preprocess_original(&mut quiet_heap, quiet_atoms);
        let mut quiet_prepared = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut quiet_data = STRUCT_DATA::default();
        let quiet_parameters = INPUT_PARMS {
            bFixNonUniformDraw: 1,
            bNoWarnings: 1,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            PreprocessOneStructure(
                &mut quiet_heap,
                None,
                &mut quiet_data,
                &quiet_parameters,
                &mut quiet_original,
                &mut quiet_prepared,
            ),
            Ok(_IS_WARNING as i32)
        );
        assert!(preprocess_error_text(&quiet_data).is_empty());

        let mut salt_heap = SourceHeap::default();
        let mut salt_atoms = vec![
            preprocess_atom(7, 1),
            preprocess_atom(8, -1),
            preprocess_atom(6, 0),
        ];
        salt_atoms[0].num_H = 4;
        preprocess_bond(&mut salt_atoms, 0, 1, 1);
        preprocess_bond(&mut salt_atoms, 1, 2, 1);
        let mut salt_original = preprocess_original(&mut salt_heap, salt_atoms);
        salt_original.num_inp_bonds = 2;
        let mut salt_prepared = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut salt_data = STRUCT_DATA::default();
        let salt_parameters = INPUT_PARMS {
            bTautFlags: TG_FLAG_DISCONNECT_SALTS as INCHI_MODE,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            PreprocessOneStructure(
                &mut salt_heap,
                None,
                &mut salt_data,
                &salt_parameters,
                &mut salt_original,
                &mut salt_prepared,
            ),
            Ok(_IS_WARNING as i32)
        );
        assert_eq!(salt_original.bDisconnectSalts, 1);
        assert_eq!(salt_prepared[0].bDisconnectSalts, 1);
        assert_eq!(salt_prepared[0].num_inp_bonds, 1);
        assert_ne!(
            salt_data.bTautFlagsDone[INCHI_BAS as usize]
                & TG_FLAG_DISCONNECT_SALTS_DONE as INCHI_MODE,
            0
        );
        assert_eq!(
            preprocess_error_text(&salt_data),
            b"Salt was disconnected; Accepted unusual valence(s): (1)"
        );

        let mut metal_heap = SourceHeap::default();
        let mut metal_atoms = vec![preprocess_atom(11, 0), preprocess_atom(6, 0)];
        preprocess_bond(&mut metal_atoms, 0, 1, 1);
        let mut metal_original = preprocess_original(&mut metal_heap, metal_atoms);
        metal_original.num_inp_bonds = 1;
        let mut metal_prepared = [ORIG_ATOM_DATA::default(), ORIG_ATOM_DATA::default()];
        let mut metal_data = STRUCT_DATA::default();
        metal_data.bTautFlags[INCHI_BAS as usize] = 0x8000_0000;
        let metal_parameters = INPUT_PARMS {
            bTautFlags: (TG_FLAG_DISCONNECT_COORD | TG_FLAG_RECONNECT_COORD) as INCHI_MODE,
            ..INPUT_PARMS::default()
        };
        assert_eq!(
            PreprocessOneStructure(
                &mut metal_heap,
                None,
                &mut metal_data,
                &metal_parameters,
                &mut metal_original,
                &mut metal_prepared,
            ),
            Ok(_IS_WARNING as i32)
        );
        assert_eq!(metal_original.bDisconnectCoord, 1);
        assert_ne!(
            metal_data.bTautFlagsDone[INCHI_BAS as usize]
                & TG_FLAG_DISCONNECT_COORD_DONE as INCHI_MODE,
            0
        );
        assert_eq!(metal_data.bTautFlags[INCHI_REC as usize], 0x8000_0000);
        assert_eq!(
            preprocess_error_text(&metal_data),
            b"Accepted unusual valence(s): (1); Metal was disconnected"
        );
        let disconnected = metal_heap.slice(metal_prepared[0].at.as_const()).unwrap();
        let saved = metal_heap.slice(metal_prepared[1].at.as_const()).unwrap();
        assert_eq!((disconnected[0].valence, disconnected[1].valence), (0, 0));
        assert_eq!((saved[0].valence, saved[1].valence), (1, 1));
        assert_ne!(metal_prepared[0].at, metal_prepared[1].at);
    }

    fn duplicate_source(heap: &mut SourceHeap) -> ORIG_ATOM_DATA {
        let mut atoms = vec![inp_ATOM::default(); 2];
        atoms[0].orig_at_number = 11;
        atoms[1].orig_at_number = 22;
        let at = heap.allocate_model_storage(atoms).unwrap();
        let current_lengths = heap.allocate_model_storage(vec![2_u16]).unwrap();
        let old_components = heap.allocate_model_storage(vec![7_u16]).unwrap();
        let coordinates = heap
            .allocate_model_storage(vec![[1_i8; 32], [2_i8; 32]])
            .unwrap();

        let alist = heap.allocate_model_storage(vec![1_i32]).unwrap();
        let blist = heap.allocate_model_storage(vec![1_i32, 2]).unwrap();
        let backbone_row = heap.allocate_model_storage(vec![3_i32, 4]).unwrap();
        let backbone = heap.allocate_model_storage(vec![backbone_row]).unwrap();
        let unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                id: 19,
                na: 1,
                nb: 1,
                alist,
                blist,
                maxbkbonds: 1,
                nbkbonds: 1,
                bkbonds: backbone,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = heap.allocate_model_storage(vec![unit]).unwrap();
        let pzz = heap.allocate_model_storage(vec![31_i32, 32]).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units,
                n: 1,
                n_pzz: 2,
                pzz,
                representation: 9,
                ..OAD_Polymer::default()
            }])
            .unwrap();

        let haptic = heap
            .allocate_model_storage(vec![10_i32, 20, 2, 30, 40])
            .unwrap();
        let steabs = heap
            .allocate_model_storage(vec![50_i32, 2, 60, 70])
            .unwrap();
        let sterel = heap.allocate_model_storage(vec![80_i32, 1, 90]).unwrap();
        let sterac = heap.allocate_model_storage(vec![100_i32, 1, 110]).unwrap();
        let haptic_lists = heap.allocate_model_storage(vec![haptic]).unwrap();
        let steabs_lists = heap.allocate_model_storage(vec![steabs]).unwrap();
        let sterel_lists = heap.allocate_model_storage(vec![sterel]).unwrap();
        let sterac_lists = heap.allocate_model_storage(vec![sterac]).unwrap();
        let atom_index_orig = heap.allocate_model_storage(vec![4_i32, 5]).unwrap();
        let atom_index_fin = heap.allocate_model_storage(vec![6_i32, 7]).unwrap();
        let v3000 = heap
            .allocate_model_storage(vec![OAD_V3000 {
                n_non_star_atoms: 2,
                atom_index_orig,
                atom_index_fin,
                n_haptic_bonds: 1,
                lists_haptic_bonds: haptic_lists,
                n_steabs: 1,
                lists_steabs: steabs_lists,
                n_sterel: 1,
                lists_sterel: sterel_lists,
                n_sterac: 1,
                lists_sterac: sterac_lists,
                ..OAD_V3000::default()
            }])
            .unwrap();

        ORIG_ATOM_DATA {
            at,
            num_dimensions: 3,
            num_inp_atoms: 2,
            num_inp_bonds: 1,
            num_components: 1,
            nCurAtLen: current_lengths,
            nOldCompNumber: old_components,
            nNumEquSets: 5,
            nEquLabels: heap.allocate_model_storage(vec![12_u16]).unwrap(),
            nSortedOrder: heap.allocate_model_storage(vec![13_u16]).unwrap(),
            bSavedInINCHI_LIB: [1, 2],
            bPreprocessed: [3, 4],
            szCoord: coordinates,
            polymer,
            v3000,
            valid_polymer: 1,
            ..ORIG_ATOM_DATA::default()
        }
    }

    #[test]
    fn source_port__runichi3__origatdata_bcheckunusualvalences__line_89() {
        fn named_atom(name: u8, charge: i8, radical: i8) -> inp_ATOM {
            let mut atom = inp_ATOM {
                el_number: 6,
                charge,
                radical,
                chem_bonds_valence: 2,
                valence: 1,
                ..inp_ATOM::default()
            };
            atom.elname[..2].copy_from_slice(&[name as i8, 0]);
            atom
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                named_atom(b'C', 3, RADICAL_SINGLET as i8),
                named_atom(b'N', -3, RADICAL_DOUBLET as i8),
                named_atom(b'O', 3, RADICAL_TRIPLET as i8),
                named_atom(b'F', -3, 4),
            ])
            .unwrap();
        let original = ORIG_ATOM_DATA {
            at: atoms,
            num_inp_atoms: 4,
            ..ORIG_ATOM_DATA::default()
        };
        let mut errors = [0_i8; 256];
        assert_eq!(
            OrigAtData_bCheckUnusualValences(&heap, Some(&original), 0, Some(&mut errors), 0,),
            Ok(4)
        );
        let error_length = errors.iter().position(|byte| *byte == 0).unwrap();
        let error_bytes = errors[..error_length]
            .iter()
            .map(|byte| *byte as u8)
            .collect::<Vec<_>>();
        assert_eq!(
            error_bytes,
            b"Accepted unusual valence(s): C+3,s(2); N-3,d(2); O+3,t(2); F-3,?(2)"
        );

        let mut isotope_atom = inp_ATOM {
            el_number: 6,
            chem_bonds_valence: 3,
            valence: 3,
            num_iso_H: [1, 0, 0],
            ..inp_ATOM::default()
        };
        isotope_atom.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
        let isotope_atoms = heap.allocate_model_storage(vec![isotope_atom]).unwrap();
        let isotope_original = ORIG_ATOM_DATA {
            at: isotope_atoms,
            num_inp_atoms: 1,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            OrigAtData_bCheckUnusualValences(&heap, Some(&isotope_original), 0, None, 1,),
            Ok(1)
        );
        assert_eq!(
            OrigAtData_bCheckUnusualValences(&heap, Some(&isotope_original), 1, None, 1,),
            Ok(0)
        );

        let mut unchanged = [0_i8; 256];
        unchanged[..2].copy_from_slice(&[b'x' as i8, 0]);
        assert_eq!(
            OrigAtData_bCheckUnusualValences(&heap, Some(&original), 0, Some(&mut unchanged), 1,),
            Ok(4)
        );
        assert_eq!(&unchanged[..2], &[b'x' as i8, 0]);
        assert_eq!(
            OrigAtData_bCheckUnusualValences(&heap, None, 1, None, 0),
            Ok(0)
        );
        assert_eq!(
            OrigAtData_bCheckUnusualValences(&heap, Some(&ORIG_ATOM_DATA::default()), 1, None, 0,),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__origatdata_duplicate__line_148() {
        let mut heap = SourceHeap::default();
        let source = duplicate_source(&mut heap);
        let mut destination = ORIG_ATOM_DATA::default();
        assert_eq!(
            OrigAtData_Duplicate(&mut heap, &mut destination, &source),
            Ok(0)
        );
        assert_eq!(destination.num_inp_atoms, 2);
        assert_eq!(destination.num_components, 1);
        assert_ne!(destination.at, source.at);
        assert_ne!(destination.nCurAtLen, source.nCurAtLen);
        assert_ne!(destination.nOldCompNumber, source.nOldCompNumber);
        assert_ne!(destination.szCoord, source.szCoord);
        assert_eq!(
            heap.slice(destination.at.as_const()).unwrap()[0].orig_at_number,
            11
        );
        assert_eq!(
            heap.slice(destination.at.as_const()).unwrap()[1].orig_at_number,
            22
        );
        assert_eq!(heap.slice(destination.nCurAtLen.as_const()).unwrap()[0], 2);
        assert_eq!(
            heap.slice(destination.nOldCompNumber.as_const()).unwrap()[0],
            7
        );
        assert_eq!(
            heap.slice(destination.szCoord.as_const()).unwrap(),
            &[[1_i8; 32], [2_i8; 32]]
        );
        assert_eq!(destination.nNumEquSets, 0);
        assert_eq!(destination.bSavedInINCHI_LIB, [0, 0]);
        assert_eq!(destination.bPreprocessed, [0, 0]);
        assert!(destination.nEquLabels.is_null());
        assert!(destination.nSortedOrder.is_null());

        assert_ne!(destination.polymer, source.polymer);
        let source_polymer = heap.slice(source.polymer.as_const()).unwrap()[0].clone();
        let copied_polymer = heap.slice(destination.polymer.as_const()).unwrap()[0].clone();
        assert_eq!(copied_polymer.n, 1);
        assert_eq!(copied_polymer.n_pzz, 2);
        assert_ne!(copied_polymer.units, source_polymer.units);
        assert_ne!(copied_polymer.pzz, source_polymer.pzz);
        assert_eq!(
            heap.slice(copied_polymer.pzz.as_const()).unwrap(),
            &[31, 32]
        );
        let source_unit = heap.slice(source_polymer.units.as_const()).unwrap()[0];
        let copied_unit = heap.slice(copied_polymer.units.as_const()).unwrap()[0];
        assert_ne!(copied_unit, source_unit);
        assert_eq!(heap.slice(copied_unit.as_const()).unwrap()[0].id, 19);

        assert_ne!(destination.v3000, source.v3000);
        let source_v3000 = heap.slice(source.v3000.as_const()).unwrap()[0].clone();
        let copied_v3000 = heap.slice(destination.v3000.as_const()).unwrap()[0].clone();
        assert_ne!(copied_v3000.atom_index_orig, source_v3000.atom_index_orig);
        assert_ne!(copied_v3000.atom_index_fin, source_v3000.atom_index_fin);
        assert_eq!(
            heap.slice(copied_v3000.atom_index_orig.as_const()).unwrap(),
            &[4, 5]
        );
        assert_eq!(
            heap.slice(copied_v3000.atom_index_fin.as_const()).unwrap(),
            &[6, 7]
        );
        for (copied_lists, source_lists, expected) in [
            (
                copied_v3000.lists_haptic_bonds,
                source_v3000.lists_haptic_bonds,
                vec![10, 20, 2, 30, 40],
            ),
            (
                copied_v3000.lists_steabs,
                source_v3000.lists_steabs,
                vec![50, 2, 60, 70],
            ),
            (
                copied_v3000.lists_sterel,
                source_v3000.lists_sterel,
                vec![80, 1, 90],
            ),
            (
                copied_v3000.lists_sterac,
                source_v3000.lists_sterac,
                vec![100, 1, 110],
            ),
        ] {
            assert_ne!(copied_lists, source_lists);
            let copied_list = heap.slice(copied_lists.as_const()).unwrap()[0];
            let source_list = heap.slice(source_lists.as_const()).unwrap()[0];
            assert_ne!(copied_list, source_list);
            assert_eq!(heap.slice(copied_list.as_const()).unwrap(), expected);
        }
        heap.slice_mut(destination.at).unwrap()[0].orig_at_number = 99;
        assert_eq!(
            heap.slice(source.at.as_const()).unwrap()[0].orig_at_number,
            11
        );

        let mut reuse_heap = SourceHeap::default();
        let mut reuse_source = duplicate_source(&mut reuse_heap);
        reuse_source.szCoord = SourceMutPointer::null();
        reuse_source.polymer = SourceMutPointer::null();
        reuse_source.v3000 = SourceMutPointer::null();
        let reused_atoms = reuse_heap
            .allocate_model_storage(vec![inp_ATOM::default(); 4])
            .unwrap();
        let reused_lengths = reuse_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let reused_old = reuse_heap.allocate_model_storage(vec![0_u16; 2]).unwrap();
        let mut reused = ORIG_ATOM_DATA {
            at: reused_atoms,
            num_inp_atoms: 4,
            nCurAtLen: reused_lengths,
            nOldCompNumber: reused_old,
            num_components: 2,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            OrigAtData_Duplicate(&mut reuse_heap, &mut reused, &reuse_source),
            Ok(0)
        );
        assert_eq!(reused.at, reused_atoms);
        assert_eq!(reused.nCurAtLen, reused_lengths);
        assert_eq!(reused.nOldCompNumber, reused_old);

        let mut early_failure_heap = SourceHeap::default();
        let early_source = duplicate_source(&mut early_failure_heap);
        let mut unchanged = ORIG_ATOM_DATA {
            num_inp_atoms: 77,
            ..ORIG_ATOM_DATA::default()
        };
        early_failure_heap.fail_after_allocations(0);
        assert_eq!(
            OrigAtData_Duplicate(&mut early_failure_heap, &mut unchanged, &early_source),
            Ok(-1)
        );
        assert_eq!(unchanged.num_inp_atoms, 77);
        assert!(unchanged.at.is_null());

        let mut coordinate_failure_heap = SourceHeap::default();
        let mut coordinate_source = duplicate_source(&mut coordinate_failure_heap);
        coordinate_source.polymer = SourceMutPointer::null();
        coordinate_source.v3000 = SourceMutPointer::null();
        let mut partial = ORIG_ATOM_DATA::default();
        coordinate_failure_heap.fail_after_allocations(3);
        assert_eq!(
            OrigAtData_Duplicate(
                &mut coordinate_failure_heap,
                &mut partial,
                &coordinate_source
            ),
            Ok(-1)
        );
        assert_eq!(partial.num_inp_atoms, coordinate_source.num_inp_atoms);
        assert!(!partial.at.is_null());
        assert!(partial.szCoord.is_null());
        assert_eq!(partial.nNumEquSets, 0);

        let mut unit_failure_heap = SourceHeap::default();
        let mut unit_source = duplicate_source(&mut unit_failure_heap);
        unit_source.v3000 = SourceMutPointer::null();
        let mut unit_destination = ORIG_ATOM_DATA::default();
        unit_failure_heap.fail_after_allocations(6);
        assert_eq!(
            OrigAtData_Duplicate(&mut unit_failure_heap, &mut unit_destination, &unit_source),
            Ok(0)
        );
        let copied_polymer = unit_failure_heap
            .slice(unit_destination.polymer.as_const())
            .unwrap()[0]
            .clone();
        assert!(
            unit_failure_heap
                .slice(copied_polymer.units.as_const())
                .unwrap()[0]
                .is_null()
        );

        let mut count_heap = SourceHeap::default();
        let count_source = duplicate_source(&mut count_heap);
        let mut counted = ORIG_ATOM_DATA::default();
        count_heap.trace_source_allocations();
        assert_eq!(
            OrigAtData_Duplicate(&mut count_heap, &mut counted, &count_source),
            Ok(0)
        );
        let allocation_count = count_heap.source_allocation_calls();
        assert!(allocation_count >= 2);

        let mut sterac_failure_heap = SourceHeap::default();
        let sterac_source = duplicate_source(&mut sterac_failure_heap);
        let mut sterac_destination = ORIG_ATOM_DATA::default();
        sterac_failure_heap.fail_after_allocations(allocation_count - 2);
        assert_eq!(
            OrigAtData_Duplicate(
                &mut sterac_failure_heap,
                &mut sterac_destination,
                &sterac_source
            ),
            Ok(0)
        );
        let copied_v3000 = sterac_failure_heap
            .slice(sterac_destination.v3000.as_const())
            .unwrap()[0]
            .clone();
        assert_eq!(copied_v3000.n_sterac, 1);
        assert!(copied_v3000.lists_sterac.is_null());
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_orderbondatomsandbondsthemselves__line_1437() {
        let mut heap = SourceHeap::default();
        let alist = heap.allocate_model_storage(vec![1_i32, 2]).unwrap();
        let stars = heap.allocate_model_storage(vec![9_i32, 8, 7]).unwrap();
        let bonds = heap.allocate_model_storage(vec![1_i32, 9, 2, 8]).unwrap();
        let mut unit = OAD_PolymerUnit {
            na: 2,
            alist,
            nb: 2,
            blist: bonds,
            ..OAD_PolymerUnit::default()
        };
        assert_eq!(
            OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
                &mut heap,
                &mut unit,
                3,
                stars.as_const(),
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(bonds.as_const()).unwrap(), &[8, 2, 9, 1]);

        let partial_bonds = heap.allocate_model_storage(vec![1_i32, 9, 8, 7]).unwrap();
        let mut partial = OAD_PolymerUnit {
            na: 2,
            alist,
            nb: 2,
            blist: partial_bonds,
            ..OAD_PolymerUnit::default()
        };
        assert_eq!(
            OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
                &mut heap,
                &mut partial,
                3,
                stars.as_const(),
            ),
            Ok(1)
        );
        assert_eq!(heap.slice(partial_bonds.as_const()).unwrap(), &[9, 1, 8, 7]);

        let external = heap.allocate_model_storage(vec![6_i32, 1]).unwrap();
        let mut single = OAD_PolymerUnit {
            na: 2,
            alist,
            nb: 1,
            blist: external,
            ..OAD_PolymerUnit::default()
        };
        assert_eq!(
            OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
                &mut heap,
                &mut single,
                0,
                SourceConstPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(external.as_const()).unwrap(), &[6, 1]);

        let mut negative = OAD_PolymerUnit {
            na: i32::MIN,
            nb: i32::MIN,
            alist: SourceMutPointer::null(),
            blist: SourceMutPointer::null(),
            ..OAD_PolymerUnit::default()
        };
        assert_eq!(
            OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves(
                &mut heap,
                &mut negative,
                i32::MIN,
                SourceConstPointer::null(),
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymer_prepareworkingset__line_2342() {
        fn allocate_unit(
            heap: &mut SourceHeap,
            atom_list: Vec<i32>,
            bond_list: Vec<i32>,
            fields: [i32; 4],
            backbone_bonds: &[[i32; 2]],
        ) -> (
            SourceMutPointer<OAD_PolymerUnit>,
            Vec<SourceMutPointer<i32>>,
        ) {
            let na = i32::try_from(atom_list.len()).unwrap();
            let nb = i32::try_from(bond_list.len() / 2).unwrap();
            let alist = if atom_list.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(atom_list).unwrap()
            };
            let blist = if bond_list.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(bond_list).unwrap()
            };
            let rows = backbone_bonds
                .iter()
                .map(|row| heap.allocate_model_storage(row.to_vec()).unwrap())
                .collect::<Vec<_>>();
            let bkbonds = if rows.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(rows.clone()).unwrap()
            };
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    na,
                    nb,
                    cap1: fields[0],
                    cap2: fields[1],
                    end_atom1: fields[2],
                    end_atom2: fields[3],
                    alist,
                    blist,
                    nbkbonds: i32::try_from(rows.len()).unwrap(),
                    maxbkbonds: i32::try_from(rows.len()).unwrap(),
                    bkbonds,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            (unit, rows)
        }

        let mut heap = SourceHeap::default();
        let canonical = heap
            .allocate_model_storage(vec![4_i32, -1, 2, 0, 5, 1, 3, 6])
            .unwrap();
        let pzz = heap.allocate_model_storage(vec![6_i32]).unwrap();
        let (first, first_rows) = allocate_unit(
            &mut heap,
            vec![3, 1, 0, 2],
            vec![2, 4, 3, 0],
            [0, 2, 3, 5],
            &[[4, 0], [1, 2], [0, 1]],
        );
        let (second, _) = allocate_unit(&mut heap, vec![7, 5], vec![], [0, 0, 0, 0], &[]);
        let units = heap.allocate_model_storage(vec![first, second]).unwrap();
        let unit_numbers = heap.allocate_model_storage(vec![-9_i32, -9]).unwrap();
        let mut polymer = OAD_Polymer {
            n: 2,
            n_pzz: 1,
            pzz,
            ..OAD_Polymer::default()
        };
        assert_eq!(
            OAD_Polymer_PrepareWorkingSet(
                &mut heap,
                &mut polymer,
                canonical.as_const(),
                SourceConstPointer::null(),
                units,
                unit_numbers,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(pzz.as_const()).unwrap(), &[4]);
        assert_eq!(heap.slice(unit_numbers.as_const()).unwrap(), &[1, 0]);
        let first_value = &heap.slice(first.as_const()).unwrap()[0];
        assert_eq!(first_value.na, 3);
        assert_eq!(
            (
                first_value.cap1,
                first_value.cap2,
                first_value.end_atom1,
                first_value.end_atom2
            ),
            (5, 3, 1, 2)
        );
        assert_eq!(
            heap.slice(first_value.alist.as_const()).unwrap(),
            &[1, 3, 5, 2]
        );
        assert_eq!(
            heap.slice(first_value.blist.as_const()).unwrap(),
            &[1, 5, 6, 3]
        );
        assert_eq!(heap.slice(first_rows[0].as_const()).unwrap(), &[5, 6]);
        assert_eq!(heap.slice(first_rows[1].as_const()).unwrap(), &[1, 2]);
        assert_eq!(heap.slice(first_rows[2].as_const()).unwrap(), &[0, 1]);
        let second_value = &heap.slice(second.as_const()).unwrap()[0];
        assert_eq!(heap.slice(second_value.alist.as_const()).unwrap(), &[2, 7]);

        let mut error10_heap = SourceHeap::default();
        let error10_canonical = error10_heap
            .allocate_model_storage(vec![0_i32, -1])
            .unwrap();
        let error10_pzz = error10_heap.allocate_model_storage(vec![0_i32, 1]).unwrap();
        let mut error10_polymer = OAD_Polymer {
            n_pzz: 2,
            pzz: error10_pzz,
            ..OAD_Polymer::default()
        };
        assert_eq!(
            OAD_Polymer_PrepareWorkingSet(
                &mut error10_heap,
                &mut error10_polymer,
                error10_canonical.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(10)
        );
        assert_eq!(error10_heap.slice(error10_pzz.as_const()).unwrap(), &[1, 1]);

        let mut bond_error_heap = SourceHeap::default();
        let bond_error_canonical = bond_error_heap
            .allocate_model_storage(vec![4_i32, 3, -1])
            .unwrap();
        let (bond_error_unit, _) =
            allocate_unit(&mut bond_error_heap, vec![0], vec![1, 2], [0, 0, 0, 0], &[]);
        let bond_error_units = bond_error_heap
            .allocate_model_storage(vec![bond_error_unit])
            .unwrap();
        let mut bond_error_polymer = OAD_Polymer {
            n: 1,
            ..OAD_Polymer::default()
        };
        assert_eq!(
            OAD_Polymer_PrepareWorkingSet(
                &mut bond_error_heap,
                &mut bond_error_polymer,
                bond_error_canonical.as_const(),
                SourceConstPointer::null(),
                bond_error_units,
                SourceMutPointer::null(),
            ),
            Ok(11)
        );
        let bond_error_value = &bond_error_heap.slice(bond_error_unit.as_const()).unwrap()[0];
        assert_eq!(bond_error_value.na, 1);
        assert_eq!(
            bond_error_heap
                .slice(bond_error_value.blist.as_const())
                .unwrap(),
            &[4, 2]
        );
        assert_eq!(bond_error_value.cap1, 0);

        for failed_field in 0..4 {
            let mut error_heap = SourceHeap::default();
            let mut canonical_values = vec![9_i32, 8, 7, 6, 5, 4];
            canonical_values[failed_field + 1] = -1;
            let canonical = error_heap.allocate_model_storage(canonical_values).unwrap();
            let (first_unit, _) =
                allocate_unit(&mut error_heap, vec![5, 0], vec![], [0, 0, 0, 0], &[]);
            let (failing_unit, _) =
                allocate_unit(&mut error_heap, vec![0], vec![], [1, 2, 3, 4], &[]);
            let units = error_heap
                .allocate_model_storage(vec![first_unit, failing_unit])
                .unwrap();
            let mut polymer = OAD_Polymer {
                n: 2,
                ..OAD_Polymer::default()
            };
            assert_eq!(
                OAD_Polymer_PrepareWorkingSet(
                    &mut error_heap,
                    &mut polymer,
                    canonical.as_const(),
                    SourceConstPointer::null(),
                    units,
                    SourceMutPointer::null(),
                ),
                Ok(11)
            );
            let first_value = &error_heap.slice(first_unit.as_const()).unwrap()[0];
            assert_eq!(
                error_heap.slice(first_value.alist.as_const()).unwrap(),
                &[5, 10]
            );
            let value = &error_heap.slice(failing_unit.as_const()).unwrap()[0];
            let actual = [value.cap1, value.cap2, value.end_atom1, value.end_atom2];
            for index in 0..4 {
                let expected = if index < failed_field {
                    9 - i32::try_from(index).unwrap()
                } else {
                    i32::try_from(index + 1).unwrap()
                };
                assert_eq!(
                    actual[index], expected,
                    "failed field {failed_field}, field {index}"
                );
            }
        }

        let mut error12_heap = SourceHeap::default();
        let error12_canonical = error12_heap
            .allocate_model_storage(vec![0_i32, 1, 2])
            .unwrap();
        let (error12_unit, _) =
            allocate_unit(&mut error12_heap, vec![0], vec![1, 2], [0, 0, 0, 0], &[]);
        let error12_units = error12_heap
            .allocate_model_storage(vec![error12_unit])
            .unwrap();
        let error12_numbers = error12_heap.allocate_model_storage(vec![-7_i32]).unwrap();
        let mut error12_polymer = OAD_Polymer {
            n: 1,
            ..OAD_Polymer::default()
        };
        assert_eq!(
            OAD_Polymer_PrepareWorkingSet(
                &mut error12_heap,
                &mut error12_polymer,
                error12_canonical.as_const(),
                SourceConstPointer::null(),
                error12_units,
                error12_numbers,
            ),
            Ok(12)
        );
        assert_eq!(
            error12_heap.slice(error12_numbers.as_const()).unwrap(),
            &[-7]
        );

        let mut empty_heap = SourceHeap::default();
        let mut empty_polymer = OAD_Polymer {
            n: i32::MIN,
            n_pzz: i32::MIN,
            pzz: SourceMutPointer::null(),
            ..OAD_Polymer::default()
        };
        assert_eq!(
            OAD_Polymer_PrepareWorkingSet(
                &mut empty_heap,
                &mut empty_polymer,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__origatdata_removehalfbond__line_2517() {
        let mut atom = inp_ATOM {
            valence: 3,
            neighbor: {
                let mut values = [0; 20];
                values[..3].copy_from_slice(&[4, 7, 9]);
                values
            },
            bond_type: {
                let mut values = [0; 20];
                values[..3].copy_from_slice(&[1, 2, 3]);
                values
            },
            bond_stereo: {
                let mut values = [0; 20];
                values[..3].copy_from_slice(&[-1, 0, 1]);
                values
            },
            ..inp_ATOM::default()
        };
        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(vec![atom.clone()]).unwrap();
        let mut bond_type = -8;
        let mut bond_stereo = -9;
        assert_eq!(
            OrigAtData_RemoveHalfBond(&mut heap, 0, 7, atoms, &mut bond_type, &mut bond_stereo,),
            Ok(1)
        );
        assert_eq!((bond_type, bond_stereo), (2, 0));
        atom = heap.slice(atoms.as_const()).unwrap()[0].clone();
        assert_eq!(&atom.neighbor[..4], &[4, 9, 0, 0]);
        assert_eq!(&atom.bond_type[..4], &[1, 3, 0, 0]);
        assert_eq!(&atom.bond_stereo[..4], &[-1, 1, 0, 0]);

        let before = atom.clone();
        assert_eq!(
            OrigAtData_RemoveHalfBond(&mut heap, 0, 99, atoms, &mut bond_type, &mut bond_stereo,),
            Ok(0)
        );
        assert_eq!(heap.slice(atoms.as_const()).unwrap()[0], before);
        assert_eq!((bond_type, bond_stereo), (2, 0));
        assert_eq!(
            OrigAtData_RemoveHalfBond(&mut heap, -1, 4, atoms, &mut bond_type, &mut bond_stereo,),
            Ok(0)
        );
        assert_eq!(
            OrigAtData_RemoveHalfBond(
                &mut heap,
                0,
                4,
                SourceMutPointer::null(),
                &mut bond_type,
                &mut bond_stereo,
            ),
            Ok(0)
        );

        let mut full = inp_ATOM {
            valence: 20,
            ..inp_ATOM::default()
        };
        for index in 0..20 {
            full.neighbor[index] = u16::try_from(index + 1).unwrap();
            full.bond_type[index] = u8::try_from(index + 1).unwrap();
            full.bond_stereo[index] = i8::try_from(index as i32 - 10).unwrap();
        }
        let full_atoms = heap.allocate_model_storage(vec![full]).unwrap();
        assert_eq!(
            OrigAtData_RemoveHalfBond(
                &mut heap,
                0,
                1,
                full_atoms,
                &mut bond_type,
                &mut bond_stereo,
            ),
            Ok(1)
        );
        assert_eq!((bond_type, bond_stereo), (1, -10));
        let full = &heap.slice(full_atoms.as_const()).unwrap()[0];
        assert_eq!(full.neighbor[0], 2);
        assert_eq!(full.bond_type[0], 2);
        assert_eq!(full.bond_stereo[0], -9);
        assert_eq!(full.neighbor[19], 0);
        assert_eq!(full.bond_type[19], 0);
        assert_eq!(full.bond_stereo[19], 0);
    }

    #[test]
    fn source_port__runichi3__origatdata_removebond__line_2577() {
        fn atom(neighbors: &[u16], types: &[u8], stereos: &[i8], chem: i8) -> inp_ATOM {
            let mut atom = inp_ATOM {
                valence: i8::try_from(neighbors.len()).unwrap(),
                chem_bonds_valence: chem,
                ..inp_ATOM::default()
            };
            atom.neighbor[..neighbors.len()].copy_from_slice(neighbors);
            atom.bond_type[..types.len()].copy_from_slice(types);
            atom.bond_stereo[..stereos.len()].copy_from_slice(stereos);
            atom
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                atom(&[1, 2], &[2, 1], &[-1, 0], 3),
                atom(&[0, 2], &[2, 3], &[1, 0], 5),
            ])
            .unwrap();
        let mut bond_type = -1;
        let mut bond_stereo = -2;
        let mut num_bonds = 2;
        assert_eq!(
            OrigAtData_RemoveBond(
                &mut heap,
                0,
                1,
                atoms,
                &mut bond_type,
                &mut bond_stereo,
                &mut num_bonds,
            ),
            Ok(1)
        );
        assert_eq!((bond_type, bond_stereo, num_bonds), (2, 1, 1));
        let values = heap.slice(atoms.as_const()).unwrap();
        assert_eq!((values[0].valence, values[0].chem_bonds_valence), (1, 1));
        assert_eq!((values[1].valence, values[1].chem_bonds_valence), (1, 3));
        assert_eq!(&values[0].neighbor[..3], &[2, 0, 0]);
        assert_eq!(&values[1].neighbor[..3], &[2, 0, 0]);

        let partial = heap
            .allocate_model_storage(vec![atom(&[1], &[3], &[-1], 3), atom(&[2], &[4], &[1], 4)])
            .unwrap();
        bond_type = -1;
        bond_stereo = -2;
        num_bonds = 7;
        assert_eq!(
            OrigAtData_RemoveBond(
                &mut heap,
                0,
                1,
                partial,
                &mut bond_type,
                &mut bond_stereo,
                &mut num_bonds,
            ),
            Ok(0)
        );
        assert_eq!((bond_type, bond_stereo, num_bonds), (3, -1, 7));
        let values = heap.slice(partial.as_const()).unwrap();
        assert_eq!(values[0].neighbor[0], 0);
        assert_eq!((values[0].valence, values[0].chem_bonds_valence), (1, 3));
        assert_eq!(values[1].neighbor[0], 2);

        let before = heap.slice(partial.as_const()).unwrap().to_vec();
        assert_eq!(
            OrigAtData_RemoveBond(
                &mut heap,
                -1,
                1,
                partial,
                &mut bond_type,
                &mut bond_stereo,
                &mut num_bonds,
            ),
            Ok(0)
        );
        assert_eq!(heap.slice(partial.as_const()).unwrap(), before.as_slice());
        assert_eq!(
            OrigAtData_RemoveBond(
                &mut heap,
                0,
                1,
                SourceMutPointer::null(),
                &mut bond_type,
                &mut bond_stereo,
                &mut num_bonds,
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__origatdata_addbond__line_2607() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default(), inp_ATOM::default()])
            .unwrap();
        let mut num_bonds = 4;
        assert_eq!(
            OrigAtData_AddBond(&mut heap, 0, 1, atoms, 99, 200, &mut num_bonds),
            Ok(1)
        );
        assert_eq!(num_bonds, 5);
        let values = heap.slice(atoms.as_const()).unwrap();
        for (atom, neighbor) in [(&values[0], 1_u16), (&values[1], 0_u16)] {
            assert_eq!(atom.valence, 1);
            assert_eq!(atom.chem_bonds_valence, 1);
            assert_eq!(atom.neighbor[0], neighbor);
            assert_eq!(atom.bond_type[0], 1);
            assert_eq!(atom.bond_stereo[0], -56);
        }

        let before = heap.slice(atoms.as_const()).unwrap().to_vec();
        assert_eq!(
            OrigAtData_AddBond(&mut heap, 0, 1, atoms, 3, -7, &mut num_bonds),
            Ok(1)
        );
        assert_eq!(num_bonds, 6);
        assert_eq!(heap.slice(atoms.as_const()).unwrap(), before.as_slice());

        let asymmetric = heap
            .allocate_model_storage(vec![
                {
                    let mut atom = inp_ATOM {
                        valence: 1,
                        chem_bonds_valence: 2,
                        ..inp_ATOM::default()
                    };
                    atom.neighbor[0] = 1;
                    atom.bond_type[0] = 2;
                    atom
                },
                inp_ATOM::default(),
            ])
            .unwrap();
        num_bonds = 0;
        assert_eq!(
            OrigAtData_AddBond(&mut heap, 0, 1, asymmetric, 2, -3, &mut num_bonds),
            Ok(1)
        );
        let values = heap.slice(asymmetric.as_const()).unwrap();
        assert_eq!(values[0].valence, 1);
        assert_eq!(values[1].valence, 1);
        assert_eq!(values[1].bond_type[0], 2);
        assert_eq!(values[1].bond_stereo[0], -3);

        let mut full = inp_ATOM {
            valence: 20,
            chem_bonds_valence: 77,
            ..inp_ATOM::default()
        };
        full.neighbor = [8; 20];
        let full_atoms = heap
            .allocate_model_storage(vec![full, inp_ATOM::default()])
            .unwrap();
        let before = heap.slice(full_atoms.as_const()).unwrap().to_vec();
        num_bonds = 10;
        assert_eq!(
            OrigAtData_AddBond(&mut heap, 0, 1, full_atoms, 2, 3, &mut num_bonds),
            Ok(0)
        );
        assert_eq!(num_bonds, 10);
        assert_eq!(
            heap.slice(full_atoms.as_const()).unwrap(),
            before.as_slice()
        );

        let self_atom = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        num_bonds = -1;
        assert_eq!(
            OrigAtData_AddBond(&mut heap, 0, 0, self_atom, 3, 1, &mut num_bonds),
            Ok(1)
        );
        let value = &heap.slice(self_atom.as_const()).unwrap()[0];
        assert_eq!((value.valence, value.chem_bonds_valence), (1, 3));
        assert_eq!((value.neighbor[0], value.bond_type[0]), (0, 3));
        assert_eq!(num_bonds, 0);

        assert_eq!(
            OrigAtData_AddBond(
                &mut heap,
                0,
                1,
                SourceMutPointer::null(),
                2,
                0,
                &mut num_bonds,
            ),
            Ok(0)
        );
        assert_eq!(num_bonds, 0);
    }

    #[test]
    fn source_port__runichi3__unmarkringsystemsinp__line_1961() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                inp_ATOM {
                    bCutVertex: -7,
                    nRingSystem: 11,
                    nNumAtInRingSystem: 12,
                    nBlockSystem: 13,
                    charge: -2,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    bCutVertex: 8,
                    nRingSystem: u16::MAX,
                    nNumAtInRingSystem: u16::MAX,
                    nBlockSystem: u16::MAX,
                    charge: 3,
                    ..inp_ATOM::default()
                },
            ])
            .unwrap();
        assert_eq!(UnMarkRingSystemsInp(&mut heap, atoms, 2), Ok(0));
        let values = heap.slice(atoms.as_const()).unwrap();
        for atom in values {
            assert_eq!(
                (
                    atom.bCutVertex,
                    atom.nRingSystem,
                    atom.nNumAtInRingSystem,
                    atom.nBlockSystem,
                ),
                (0, 0, 0, 0)
            );
        }
        assert_eq!((values[0].charge, values[1].charge), (-2, 3));
        assert_eq!(
            UnMarkRingSystemsInp(&mut heap, SourceMutPointer::null(), 0),
            Ok(0)
        );
        assert_eq!(
            UnMarkRingSystemsInp(&mut heap, SourceMutPointer::null(), i32::MIN),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__origatdata_addsinglestereolessbond__line_2682() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default(), inp_ATOM::default()])
            .unwrap();
        let mut num_bonds = 0;
        assert_eq!(
            OrigAtData_AddSingleStereolessBond(&mut heap, 0, 1, atoms, &mut num_bonds,),
            Ok(1)
        );
        assert_eq!(num_bonds, 1);
        let values = heap.slice(atoms.as_const()).unwrap();
        assert_eq!((values[0].bond_type[0], values[0].bond_stereo[0]), (1, 0));
        assert_eq!((values[1].bond_type[0], values[1].bond_stereo[0]), (1, 0));
        assert_eq!(
            OrigAtData_AddSingleStereolessBond(
                &mut heap,
                0,
                1,
                SourceMutPointer::null(),
                &mut num_bonds,
            ),
            Ok(0)
        );
        assert_eq!(num_bonds, 1);
    }

    #[test]
    fn source_port__runichi3__oad_polymer_findringsystems__line_3236() {
        let mut heap = SourceHeap::default();
        let mut triangle_atoms = vec![inp_ATOM::default(); 3];
        triangle_atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
        triangle_atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
        triangle_atoms[2].neighbor[..2].copy_from_slice(&[0, 1]);
        for (index, atom) in triangle_atoms.iter_mut().enumerate() {
            atom.valence = 2;
            atom.orig_at_number = u16::try_from(index).unwrap();
        }
        let atoms = heap.allocate_model_storage(triangle_atoms).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        let ring_numbers = heap.allocate_model_storage(vec![-8_i32; 4]).unwrap();
        let ring_sizes = heap.allocate_model_storage(vec![-9_i32; 4]).unwrap();
        let mut num_bonds = 3;
        assert_eq!(
            OAD_Polymer_FindRingSystems(
                &mut heap,
                polymer.as_const(),
                atoms,
                3,
                &mut num_bonds,
                ring_numbers,
                ring_sizes,
                0,
            ),
            Ok(3)
        );
        assert_eq!(num_bonds, 3);
        assert_eq!(heap.slice(ring_numbers.as_const()).unwrap(), &[1, 1, 1, -1]);
        assert_eq!(heap.slice(ring_sizes.as_const()).unwrap(), &[3, 3, 3, -9]);
        assert!(
            heap.slice(atoms.as_const())
                .unwrap()
                .iter()
                .all(|atom| atom.nRingSystem == 0 && atom.nBlockSystem == 0)
        );

        assert_eq!(
            OAD_Polymer_FindRingSystems(
                &mut heap,
                SourceConstPointer::null(),
                SourceMutPointer::null(),
                3,
                &mut num_bonds,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
            ),
            Ok(0)
        );

        let mut path_atoms = vec![inp_ATOM::default(); 3];
        path_atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
        path_atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
        path_atoms[2].neighbor[..2].copy_from_slice(&[0, 1]);
        for atom in &mut path_atoms {
            atom.valence = 2;
            atom.chem_bonds_valence = 2;
            atom.bond_type[..2].copy_from_slice(&[1, 1]);
        }
        let path_atoms = heap.allocate_model_storage(path_atoms).unwrap();
        let unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                cyclized: 1,
                end_atom1: 1,
                end_atom2: 3,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = heap.allocate_model_storage(vec![unit]).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                n: 1,
                units,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let numbers = heap.allocate_model_storage(vec![-4_i32; 4]).unwrap();
        let mut num_bonds = 3;
        assert_eq!(
            OAD_Polymer_FindRingSystems(
                &mut heap,
                polymer.as_const(),
                path_atoms,
                3,
                &mut num_bonds,
                numbers,
                SourceMutPointer::null(),
                0,
            ),
            Ok(0)
        );
        assert_eq!(num_bonds, 3);
        assert_eq!(heap.slice(numbers.as_const()).unwrap(), &[-1, -1, -1, -1]);
    }

    #[test]
    fn source_port__runichi3__oad_polymer_setatprops__line_3312() {
        let mut heap = SourceHeap::default();
        let mut atoms_values = vec![inp_ATOM::default(); 3];
        atoms_values[0].neighbor[..2].copy_from_slice(&[1, 2]);
        atoms_values[1].neighbor[..2].copy_from_slice(&[0, 2]);
        atoms_values[2].neighbor[..2].copy_from_slice(&[0, 1]);
        for (index, atom) in atoms_values.iter_mut().enumerate() {
            atom.valence = 2;
            atom.orig_at_number = u16::try_from(index).unwrap();
            atom.el_number = [6, 7, 8][index];
        }
        let atoms = heap.allocate_model_storage(atoms_values).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        let properties = heap
            .allocate_model_storage(vec![OAD_AtProps::default(); 3])
            .unwrap();
        let mut num_bonds = 3;
        assert_eq!(
            OAD_Polymer_SetAtProps(
                &mut heap,
                polymer.as_const(),
                atoms,
                3,
                &mut num_bonds,
                properties,
                SourceConstPointer::null(),
            ),
            Ok(())
        );
        let values = heap.slice(properties.as_const()).unwrap();
        assert_eq!(
            values,
            &[
                OAD_AtProps {
                    erank: 2,
                    ring_erank: 216,
                    ring_num: 1,
                    ring_size: 3,
                },
                OAD_AtProps {
                    erank: 211,
                    ring_erank: 216,
                    ring_num: 1,
                    ring_size: 3,
                },
                OAD_AtProps {
                    erank: 215,
                    ring_erank: 216,
                    ring_num: 1,
                    ring_size: 3,
                },
            ]
        );
        assert_eq!(num_bonds, 3);

        let canonical = heap.allocate_model_storage(vec![2_i32, 0, 1]).unwrap();
        let mapped_properties = heap
            .allocate_model_storage(vec![OAD_AtProps::default(); 4])
            .unwrap();
        assert_eq!(
            OAD_Polymer_SetAtProps(
                &mut heap,
                polymer.as_const(),
                atoms,
                3,
                &mut num_bonds,
                mapped_properties,
                canonical.as_const(),
            ),
            Ok(())
        );
        let values = heap.slice(mapped_properties.as_const()).unwrap();
        assert_eq!(
            (values[1].ring_num, values[2].ring_num, values[3].ring_num),
            (1, 1, 1)
        );
        assert_eq!(
            (
                values[1].ring_size,
                values[2].ring_size,
                values[3].ring_size
            ),
            (3, 3, 3)
        );

        let failure_properties = heap
            .allocate_model_storage(vec![OAD_AtProps::default(); 3])
            .unwrap();
        heap.trace_source_allocations();
        heap.fail_after_allocations(0);
        assert_eq!(
            OAD_Polymer_SetAtProps(
                &mut heap,
                polymer.as_const(),
                atoms,
                3,
                &mut num_bonds,
                failure_properties,
                SourceConstPointer::null(),
            ),
            Ok(())
        );
        assert_eq!(
            heap.slice(failure_properties.as_const()).unwrap()[0].erank,
            2
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_compareatomlistsmod__line_1377() {
        let mut heap = SourceHeap::default();
        let first_values = heap.allocate_model_storage(vec![1_i32, 4]).unwrap();
        let second_values = heap.allocate_model_storage(vec![1_i32, 5]).unwrap();
        let first = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 2,
                alist: first_values,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let second = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 2,
                alist: second_values,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(&heap, first.as_const(), second.as_const()),
            Ok(-1)
        );
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(&heap, second.as_const(), first.as_const()),
            Ok(1)
        );

        let same_values = heap.allocate_model_storage(vec![1_i32, 4]).unwrap();
        let same = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 2,
                alist: same_values,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(&heap, first.as_const(), same.as_const()),
            Ok(0)
        );

        let short = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 1,
                alist: first_values,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(&heap, short.as_const(), first.as_const()),
            Ok(-1)
        );
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(&heap, first.as_const(), short.as_const()),
            Ok(1)
        );

        let negative = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: -1,
                alist: SourceMutPointer::null(),
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let negative_other = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: -1,
                alist: SourceMutPointer::null(),
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_PolymerUnit_CompareAtomListsMod(
                &heap,
                negative.as_const(),
                negative_other.as_const()
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_new__line_1162() {
        use crate::source_types::INT_ARRAY;

        fn c_string(heap: &mut SourceHeap, text: &str) -> SourceConstPointer<i8> {
            let mut bytes: Vec<i8> = text.bytes().map(|byte| byte as i8).collect();
            bytes.push(0);
            heap.allocate_model_storage(bytes).unwrap().as_const()
        }

        fn list(heap: &mut SourceHeap, values: &[i32]) -> INT_ARRAY {
            INT_ARRAY {
                item: heap.allocate_model_storage(values.to_vec()).unwrap(),
                allocated: values.len() as i32,
                used: values.len() as i32,
                increment: values.len() as i32,
            }
        }

        let mut heap = SourceHeap::default();
        let smt = c_string(&mut heap, "HT");
        let atoms = list(&mut heap, &[4, 7, 9]);
        let bonds = list(&mut heap, &[1, 2, 3, 4]);
        let ignored_row = heap.allocate_model_storage(vec![31_i32, 32]).unwrap();
        let ignored_bkbonds = heap.allocate_model_storage(vec![ignored_row]).unwrap();
        let unit = OAD_PolymerUnit_New(
            &mut heap,
            20,
            30,
            11,
            12,
            13,
            14,
            15,
            smt,
            3,
            Some(&atoms),
            2,
            Some(&bonds),
            17,
            ignored_bkbonds,
        )
        .unwrap();
        assert!(!unit.is_null());
        let value = &heap.slice(unit.as_const()).unwrap()[0];
        assert_eq!(
            (
                value.id,
                value.label,
                value.type_,
                value.subtype,
                value.conn,
                value.na,
                value.nb,
            ),
            (11, 12, 13, 14, 15, 3, 2)
        );
        assert_eq!(value.cyclizable, CLOSING_SRU_NOT_APPLICABLE as i32);
        assert_eq!(value.cyclized, 0);
        assert_eq!(value.xbr1, [0.0; 4]);
        assert_eq!(value.xbr2, [0.0; 4]);
        assert_eq!(&value.smt[..3], &[b'H' as i8, b'T' as i8, 0]);
        assert!(value.smt[3..].iter().all(|byte| *byte == 0));
        assert_eq!(
            (
                value.representation,
                value.cap1,
                value.end_atom1,
                value.cap2,
                value.end_atom2,
                value.maxbkbonds,
                value.nbkbonds,
                value.cap1_is_undef,
                value.cap2_is_undef,
            ),
            (0, -1, -1, -1, -1, 30, 17, 0, 0)
        );
        assert_eq!(heap.slice(value.alist.as_const()).unwrap(), &[4, 7, 9]);
        assert_eq!(heap.slice(value.blist.as_const()).unwrap(), &[1, 2, 3, 4]);
        assert!(value.bkbonds.is_null());

        let mut capacity_heap = SourceHeap::default();
        let smt = c_string(&mut capacity_heap, "");
        let ignored_empty = INT_ARRAY {
            item: SourceMutPointer::null(),
            allocated: 99,
            used: 99,
            increment: 99,
        };
        let unit = OAD_PolymerUnit_New(
            &mut capacity_heap,
            3,
            2,
            1,
            2,
            3,
            4,
            5,
            smt,
            0,
            None,
            0,
            Some(&ignored_empty),
            -8,
            SourceMutPointer::null(),
        )
        .unwrap();
        let value = &capacity_heap.slice(unit.as_const()).unwrap()[0];
        assert_eq!(
            capacity_heap.slice(value.alist.as_const()).unwrap(),
            &[0, 0, 0]
        );
        assert_eq!(
            capacity_heap.slice(value.blist.as_const()).unwrap(),
            &[0, 0, 0, 0]
        );
        assert_eq!((value.na, value.nb, value.nbkbonds), (0, 0, -8));

        let mut null_heap = SourceHeap::default();
        let smt = c_string(&mut null_heap, "S");
        let unit = OAD_PolymerUnit_New(
            &mut null_heap,
            0,
            0,
            i32::MIN,
            i32::MAX,
            -1,
            -2,
            -3,
            smt,
            0,
            None,
            0,
            None,
            0,
            SourceMutPointer::null(),
        )
        .unwrap();
        let value = &null_heap.slice(unit.as_const()).unwrap()[0];
        assert!(value.alist.is_null());
        assert!(value.blist.is_null());
        assert_eq!((value.id, value.label), (i32::MIN, i32::MAX));

        let mut zero_blist_heap = SourceHeap::default();
        let smt = c_string(&mut zero_blist_heap, "SRU");
        let unit = OAD_PolymerUnit_New(
            &mut zero_blist_heap,
            0,
            2,
            0,
            0,
            0,
            0,
            0,
            smt,
            0,
            None,
            1,
            None,
            0,
            SourceMutPointer::null(),
        )
        .unwrap();
        let value = &zero_blist_heap.slice(unit.as_const()).unwrap()[0];
        assert_eq!(
            zero_blist_heap.slice(value.blist.as_const()).unwrap(),
            &[0, 0]
        );

        for failure_after in 0_u64..=2 {
            let mut failure_heap = SourceHeap::default();
            let smt = c_string(&mut failure_heap, "X");
            let atoms = list(&mut failure_heap, &[1]);
            let bonds = list(&mut failure_heap, &[1, 2]);
            failure_heap.fail_after_allocations(failure_after);
            assert_eq!(
                OAD_PolymerUnit_New(
                    &mut failure_heap,
                    1,
                    1,
                    1,
                    2,
                    3,
                    4,
                    5,
                    smt,
                    1,
                    Some(&atoms),
                    1,
                    Some(&bonds),
                    1,
                    SourceMutPointer::null(),
                ),
                Ok(SourceMutPointer::null()),
                "failure_after={failure_after}"
            );
        }
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_createcopy__line_1260() {
        fn source_unit(heap: &mut SourceHeap) -> SourceMutPointer<OAD_PolymerUnit> {
            let alist = heap.allocate_model_storage(vec![4_i32, 5]).unwrap();
            let blist = heap.allocate_model_storage(vec![1_i32, 2, 3, 4]).unwrap();
            let first_row = heap.allocate_model_storage(vec![6_i32, 7]).unwrap();
            let second_row = heap.allocate_model_storage(vec![8_i32, 9]).unwrap();
            let bkbonds = heap
                .allocate_model_storage(vec![first_row, second_row])
                .unwrap();
            let mut smt = [99_i8; 80];
            smt[..4].copy_from_slice(&[b'A' as i8, b'B' as i8, b'C' as i8, 0]);
            heap.allocate_model_storage(vec![OAD_PolymerUnit {
                id: 11,
                type_: 12,
                subtype: 13,
                conn: 14,
                label: 15,
                na: 2,
                nb: 2,
                cyclizable: 16,
                cyclized: 17,
                xbr1: [1.0, -0.0, f64::INFINITY, f64::NAN],
                xbr2: [-1.0, 0.0, f64::NEG_INFINITY, f64::MIN],
                smt,
                representation: 91,
                cap1: 18,
                end_atom1: 19,
                cap2: 20,
                end_atom2: 21,
                cap1_is_undef: 22,
                cap2_is_undef: 23,
                alist,
                blist,
                maxbkbonds: 1,
                nbkbonds: 2,
                bkbonds,
            }])
            .unwrap()
        }

        let mut heap = SourceHeap::default();
        let source = source_unit(&mut heap);
        let copy = OAD_PolymerUnit_CreateCopy(&mut heap, source).unwrap();
        assert!(!copy.is_null());
        let source_value = heap.slice(source.as_const()).unwrap()[0].clone();
        let copy_value = heap.slice(copy.as_const()).unwrap()[0].clone();
        assert_eq!(copy_value.id, source_value.id);
        assert_eq!(copy_value.type_, source_value.type_);
        assert_eq!(copy_value.subtype, source_value.subtype);
        assert_eq!(copy_value.conn, source_value.conn);
        assert_eq!(copy_value.label, source_value.label);
        assert_eq!(copy_value.na, source_value.na);
        assert_eq!(copy_value.nb, source_value.nb);
        assert_eq!(copy_value.cyclizable, source_value.cyclizable);
        assert_eq!(copy_value.cyclized, source_value.cyclized);
        for (actual, expected) in copy_value.xbr1.iter().zip(source_value.xbr1.iter()) {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
        for (actual, expected) in copy_value.xbr2.iter().zip(source_value.xbr2.iter()) {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
        assert_eq!(
            &copy_value.smt[..4],
            &[b'A' as i8, b'B' as i8, b'C' as i8, 0]
        );
        assert!(copy_value.smt[4..].iter().all(|byte| *byte == 0));
        assert_eq!(copy_value.representation, 0);
        assert_eq!(copy_value.maxbkbonds, 2);
        assert_ne!(copy_value.alist, source_value.alist);
        assert_ne!(copy_value.blist, source_value.blist);
        assert_ne!(copy_value.bkbonds, source_value.bkbonds);
        assert_eq!(heap.slice(copy_value.alist.as_const()).unwrap(), &[4, 5]);
        assert_eq!(
            heap.slice(copy_value.blist.as_const()).unwrap(),
            &[1, 2, 3, 4]
        );
        let copy_rows = heap.slice(copy_value.bkbonds.as_const()).unwrap();
        let source_rows = heap.slice(source_value.bkbonds.as_const()).unwrap();
        assert_ne!(copy_rows[0], source_rows[0]);
        assert_ne!(copy_rows[1], source_rows[1]);
        assert_eq!(heap.slice(copy_rows[0].as_const()).unwrap(), &[6, 7]);
        assert_eq!(heap.slice(copy_rows[1].as_const()).unwrap(), &[8, 9]);
        OAD_PolymerUnit_Free(&mut heap, copy).unwrap();
        assert!(heap.slice(source.as_const()).is_ok());
        OAD_PolymerUnit_Free(&mut heap, source).unwrap();

        for successful_allocations in 0..=5 {
            let mut failure_heap = SourceHeap::default();
            let source = source_unit(&mut failure_heap);
            let source_value = failure_heap.slice(source.as_const()).unwrap()[0].clone();
            failure_heap.fail_after_allocations(successful_allocations);
            assert!(
                OAD_PolymerUnit_CreateCopy(&mut failure_heap, source)
                    .unwrap()
                    .is_null(),
                "successful_allocations={successful_allocations}"
            );
            assert!(failure_heap.slice(source.as_const()).is_ok());
            assert_eq!(
                failure_heap.slice(source_value.alist.as_const()).unwrap(),
                &[4, 5]
            );
            OAD_PolymerUnit_Free(&mut failure_heap, source).unwrap();
        }
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_free__line_1342() {
        let mut heap = SourceHeap::default();

        OAD_PolymerUnit_Free(&mut heap, SourceMutPointer::null()).unwrap();

        let empty_unit = heap.allocate(vec![OAD_PolymerUnit::default()]).unwrap();
        OAD_PolymerUnit_Free(&mut heap, empty_unit).unwrap();
        assert_eq!(
            heap.slice(empty_unit.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let alist = heap.allocate(vec![1_i32, 2, 3]).unwrap();
        let blist = heap.allocate(vec![1_i32, 2, 2, 3]).unwrap();
        let row_zero = heap.allocate(vec![1_i32, 2]).unwrap();
        let row_one = heap.allocate(vec![2_i32, 3]).unwrap();
        let bkbonds = heap.allocate(vec![row_zero, row_one]).unwrap();
        let complete_unit = heap
            .allocate(vec![OAD_PolymerUnit {
                alist,
                blist,
                maxbkbonds: 2,
                bkbonds,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        OAD_PolymerUnit_Free(&mut heap, complete_unit).unwrap();
        for result in [
            heap.slice(alist.as_const()),
            heap.slice(blist.as_const()),
            heap.slice(row_zero.as_const()),
            heap.slice(row_one.as_const()),
        ] {
            assert_eq!(result, Err(SourceHeapError::MissingAllocation));
        }
        assert_eq!(
            heap.slice(bkbonds.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(complete_unit.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let only_blist = heap.allocate(vec![4_i32, 5]).unwrap();
        let mixed_unit = heap
            .allocate(vec![OAD_PolymerUnit {
                na: i32::MAX,
                nb: i32::MIN,
                nbkbonds: i32::MAX,
                alist: SourceMutPointer::null(),
                blist: only_blist,
                bkbonds: SourceMutPointer::null(),
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        OAD_PolymerUnit_Free(&mut heap, mixed_unit).unwrap();
        assert_eq!(
            heap.slice(only_blist.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let unvisited_row = heap.allocate(vec![6_i32, 7]).unwrap();
        let negative_bkbonds = heap.allocate(vec![unvisited_row]).unwrap();
        let negative_count_unit = heap
            .allocate(vec![OAD_PolymerUnit {
                maxbkbonds: i32::MIN,
                bkbonds: negative_bkbonds,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        OAD_PolymerUnit_Free(&mut heap, negative_count_unit).unwrap();
        assert_eq!(heap.slice(unvisited_row.as_const()).unwrap(), &[6, 7]);
        assert_eq!(
            heap.slice(negative_bkbonds.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        heap.free(unvisited_row).unwrap();
    }

    #[test]
    fn source_port__runichi3__oad_polymer_free__line_2071() {
        let mut heap = SourceHeap::default();

        OAD_Polymer_Free(&mut heap, SourceMutPointer::null()).unwrap();

        let empty_polymer = heap.allocate(vec![OAD_Polymer::default()]).unwrap();
        OAD_Polymer_Free(&mut heap, empty_polymer).unwrap();
        assert_eq!(
            heap.slice(empty_polymer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let first_alist = heap.allocate(vec![1_i32]).unwrap();
        let first_unit = heap
            .allocate(vec![OAD_PolymerUnit {
                alist: first_alist,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let second_blist = heap.allocate(vec![1_i32, 2]).unwrap();
        let second_unit = heap
            .allocate(vec![OAD_PolymerUnit {
                blist: second_blist,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = heap.allocate(vec![first_unit, second_unit]).unwrap();
        let pzz = heap.allocate(vec![7_i32, 8]).unwrap();
        let complete_polymer = heap
            .allocate(vec![OAD_Polymer {
                units,
                n: 2,
                n_pzz: 2,
                pzz,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        OAD_Polymer_Free(&mut heap, complete_polymer).unwrap();
        assert_eq!(
            heap.slice(first_alist.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(second_blist.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(first_unit.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(second_unit.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(units.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(pzz.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(complete_polymer.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let retained_zero_unit = heap.allocate(vec![OAD_PolymerUnit::default()]).unwrap();
        let retained_zero_units = heap.allocate(vec![retained_zero_unit]).unwrap();
        let zero_count_polymer = heap
            .allocate(vec![OAD_Polymer {
                units: retained_zero_units,
                n: 0,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        OAD_Polymer_Free(&mut heap, zero_count_polymer).unwrap();
        assert_eq!(
            heap.slice(retained_zero_units.as_const()).unwrap(),
            &[retained_zero_unit]
        );
        assert!(heap.slice(retained_zero_unit.as_const()).is_ok());
        OAD_PolymerUnit_Free(&mut heap, retained_zero_unit).unwrap();
        heap.free(retained_zero_units).unwrap();

        let retained_negative_unit = heap.allocate(vec![OAD_PolymerUnit::default()]).unwrap();
        let negative_units = heap.allocate(vec![retained_negative_unit]).unwrap();
        let negative_count_polymer = heap
            .allocate(vec![OAD_Polymer {
                units: negative_units,
                n: i32::MIN,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        OAD_Polymer_Free(&mut heap, negative_count_polymer).unwrap();
        assert_eq!(
            heap.slice(negative_units.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert!(heap.slice(retained_negative_unit.as_const()).is_ok());
        OAD_PolymerUnit_Free(&mut heap, retained_negative_unit).unwrap();

        let no_units_pzz = heap.allocate(vec![9_i32]).unwrap();
        let no_units_polymer = heap
            .allocate(vec![OAD_Polymer {
                units: SourceMutPointer::null(),
                n: i32::MAX,
                n_pzz: i32::MIN,
                pzz: no_units_pzz,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        OAD_Polymer_Free(&mut heap, no_units_polymer).unwrap();
        assert_eq!(
            heap.slice(no_units_pzz.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_debugtrace__line_3642() {
        OAD_PolymerUnit_DebugTrace(None);

        for connection in [i32::MIN, 0, 1, 2, 3, 4, i32::MAX] {
            for unit_type in [i32::MIN, -1, 0, 1, 2, 3, 4, 5, 6, i32::MAX] {
                for subtype in [i32::MIN, 0, 1, 2, 3, 4, i32::MAX] {
                    let unit = OAD_PolymerUnit {
                        conn: connection,
                        type_: unit_type,
                        subtype,
                        na: i32::MIN,
                        nb: i32::MAX,
                        nbkbonds: i32::MAX,
                        alist: SourceMutPointer::null(),
                        blist: SourceMutPointer::null(),
                        bkbonds: SourceMutPointer::null(),
                        ..OAD_PolymerUnit::default()
                    };
                    let before = unit.clone();
                    OAD_PolymerUnit_DebugTrace(Some(&unit));
                    assert_eq!(unit, before);
                }
            }
        }
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_setreopeningdetails__line_4007() {
        fn run(
            nbkbonds: i32,
            endpoints: [i32; 2],
            atom: inp_ATOM,
            initial_cyclizable: i32,
        ) -> (i32, OAD_PolymerUnit) {
            let mut heap = SourceHeap::default();
            let row = heap.allocate_model_storage(endpoints.to_vec()).unwrap();
            let rows = heap.allocate_model_storage(vec![row]).unwrap();
            let mut unit = OAD_PolymerUnit {
                nbkbonds,
                bkbonds: rows,
                end_atom1: 91,
                end_atom2: 92,
                cyclizable: initial_cyclizable,
                ..OAD_PolymerUnit::default()
            };
            let mut atoms = vec![inp_ATOM::default(); 4];
            atoms[0] = atom;
            let status = OAD_PolymerUnit_SetReopeningDetails(&heap, &mut unit, &atoms).unwrap();
            (status, unit)
        }

        for count in [i32::MIN, -1, 0, 2, i32::MAX] {
            let (status, unit) = run(count, [1, 2], inp_ATOM::default(), 7);
            assert_eq!(status, count);
            assert_eq!(
                (unit.end_atom1, unit.end_atom2, unit.cyclizable),
                (91, 92, 7)
            );
        }

        let (_, diradical) = run(1, [1, 1], inp_ATOM::default(), 7);
        assert_eq!((diradical.end_atom1, diradical.end_atom2), (1, 1));
        assert_eq!(diradical.cyclizable, CLOSING_SRU_DIRADICAL as i32);

        let mut higher_order = inp_ATOM {
            valence: 2,
            ..inp_ATOM::default()
        };
        higher_order.neighbor[..2].copy_from_slice(&[3, 1]);
        higher_order.bond_type[..2].copy_from_slice(&[1, 2]);
        let (_, higher) = run(1, [1, 2], higher_order, 7);
        assert_eq!(higher.cyclizable, CLOSING_SRU_HIGHER_ORDER_BOND as i32);

        let mut first_match_is_single = inp_ATOM {
            valence: 3,
            ..inp_ATOM::default()
        };
        first_match_is_single.neighbor[..3].copy_from_slice(&[1, 1, 2]);
        first_match_is_single.bond_type[..3].copy_from_slice(&[1, 3, 3]);
        let (_, single) = run(1, [1, 2], first_match_is_single, 7);
        assert_eq!(single.cyclizable, 7);

        let mut no_match = inp_ATOM {
            valence: -1,
            ..inp_ATOM::default()
        };
        no_match.neighbor[0] = 1;
        no_match.bond_type[0] = 3;
        let (_, unchanged) = run(1, [1, 2], no_match, 7);
        assert_eq!(unchanged.cyclizable, 7);
    }

    #[test]
    fn source_port__runichi3__oad_polymer_compareranksoftwoatoms__line_4209() {
        let categories = [
            OAD_AtProps {
                erank: 2,
                ring_size: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 2,
                ring_erank: 2,
                ring_size: 3,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 3,
                ring_size: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 3,
                ring_erank: 3,
                ring_size: 3,
                ..OAD_AtProps::default()
            },
        ];
        for first in 0..categories.len() {
            for second in 0..categories.len() {
                let expected = if first > second {
                    -1
                } else if first < second {
                    1
                } else {
                    0
                };
                assert_eq!(
                    OAD_Polymer_CompareRanksOfTwoAtoms(
                        i32::try_from(first + 1).unwrap(),
                        i32::try_from(second + 1).unwrap(),
                        &categories,
                    )
                    .unwrap(),
                    expected
                );
            }
        }

        for (first, second, expected) in [
            (
                OAD_AtProps {
                    ring_erank: 3,
                    ring_size: 4,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    ring_erank: 4,
                    ring_size: 3,
                    ..OAD_AtProps::default()
                },
                1,
            ),
            (
                OAD_AtProps {
                    ring_erank: 4,
                    ring_size: 3,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    ring_erank: 3,
                    ring_size: 4,
                    ..OAD_AtProps::default()
                },
                -1,
            ),
            (
                OAD_AtProps {
                    ring_erank: 3,
                    ring_size: 5,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    ring_erank: 3,
                    ring_size: 4,
                    ..OAD_AtProps::default()
                },
                -1,
            ),
            (
                OAD_AtProps {
                    erank: 8,
                    ring_size: 2,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 7,
                    ring_size: 2,
                    ..OAD_AtProps::default()
                },
                -1,
            ),
            (
                OAD_AtProps {
                    erank: 2,
                    ring_erank: 2,
                    ring_size: 8,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 2,
                    ring_erank: 2,
                    ring_size: 7,
                    ..OAD_AtProps::default()
                },
                -1,
            ),
        ] {
            let values = [first, second];
            assert_eq!(
                OAD_Polymer_CompareRanksOfTwoAtoms(1, 2, &values).unwrap(),
                expected
            );
            assert_eq!(
                OAD_Polymer_CompareRanksOfTwoAtoms(2, 1, &values).unwrap(),
                -expected
            );
        }

        assert_eq!(OAD_Polymer_CompareRanksOfTwoAtoms(0, 1, &[]).unwrap(), 0);
        assert_eq!(OAD_Polymer_CompareRanksOfTwoAtoms(-1, 1, &[]).unwrap(), 0);
        assert!(matches!(
            OAD_Polymer_CompareRanksOfTwoAtoms(i32::MIN, 1, &[]),
            Err(SourceHeapError::SourceIntegerOverflow)
        ));
    }

    #[test]
    fn source_port__runichi3__oad_polymer_isfirstatomranklower__line_4369() {
        let properties = [
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 8,
                ..OAD_AtProps::default()
            },
        ];
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(1, 2, &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(2, 1, &properties).unwrap(),
            -1
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(1, 1, &properties).unwrap(),
            0
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(1, 3, &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(3, 1, &properties).unwrap(),
            -1
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(0, 1, &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_IsFirstAtomRankLower(-1, -2, &properties).unwrap(),
            -1
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymer_comparebackbonebondsseniority__line_4128() {
        let carbon = OAD_AtProps {
            erank: 2,
            ..OAD_AtProps::default()
        };
        let heteroatom = OAD_AtProps {
            erank: 8,
            ..OAD_AtProps::default()
        };
        let properties = [
            carbon.clone(),
            carbon.clone(),
            carbon.clone(),
            carbon,
            heteroatom,
        ];

        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[1, 5], &[2, 3], &properties).unwrap(),
            -1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[2, 3], &[5, 1], &properties).unwrap(),
            1
        );

        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[1, 4], &[2, 5], &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[2, 5], &[4, 1], &properties).unwrap(),
            -1
        );

        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[1, 3], &[2, 4], &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[2, 4], &[3, 1], &properties).unwrap(),
            -1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[4, 1], &[4, 2], &properties).unwrap(),
            -1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[4, 2], &[1, 4], &properties).unwrap(),
            1
        );
        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[3, 1], &[1, 3], &properties).unwrap(),
            0
        );

        assert_eq!(
            OAD_Polymer_CompareBackboneBondsSeniority(&[i32::MIN, 1], &[2, 3], &properties),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_sortbackbonebonds__line_4097() {
        let mut heap = SourceHeap::default();
        let row0 = heap.allocate(vec![1_i32, 2]).unwrap();
        let row1 = heap.allocate(vec![3_i32, 5]).unwrap();
        let row2 = heap.allocate(vec![1_i32, 4]).unwrap();
        let row3 = heap.allocate(vec![2_i32, 1]).unwrap();
        let rows = heap.allocate(vec![row0, row1, row2, row3]).unwrap();
        let unit = OAD_PolymerUnit {
            nbkbonds: 4,
            bkbonds: rows,
            ..OAD_PolymerUnit::default()
        };
        let properties = [
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 2,
                ..OAD_AtProps::default()
            },
            OAD_AtProps {
                erank: 8,
                ..OAD_AtProps::default()
            },
        ];

        let bond_numbers = heap.allocate(vec![0, 1, 2, 3]).unwrap();
        OAD_PolymerUnit_SortBackboneBonds(&mut heap, &unit, &properties, bond_numbers).unwrap();
        assert_eq!(heap.slice(bond_numbers.as_const()).unwrap(), &[1, 2, 0, 3]);

        let unchanged = heap.allocate(vec![3, 2, 1, 0]).unwrap();
        let no_bonds = OAD_PolymerUnit {
            nbkbonds: 0,
            ..unit.clone()
        };
        OAD_PolymerUnit_SortBackboneBonds(&mut heap, &no_bonds, &properties, unchanged).unwrap();
        assert_eq!(heap.slice(unchanged.as_const()).unwrap(), &[3, 2, 1, 0]);
        OAD_PolymerUnit_SortBackboneBonds(&mut heap, &unit, &properties, SourceMutPointer::null())
            .unwrap();

        let invalid = heap.allocate(vec![0, -1, 2, 3]).unwrap();
        assert_eq!(
            OAD_PolymerUnit_SortBackboneBonds(&mut heap, &unit, &properties, invalid,),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_sortbackbonebondsandsetseniors__line_4056() {
        fn properties() -> Vec<OAD_AtProps> {
            vec![
                OAD_AtProps {
                    erank: 2,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 2,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 2,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 2,
                    ..OAD_AtProps::default()
                },
                OAD_AtProps {
                    erank: 8,
                    ..OAD_AtProps::default()
                },
            ]
        }

        let mut heap = SourceHeap::default();
        let row0 = heap.allocate(vec![1_i32, 2]).unwrap();
        let row1 = heap.allocate(vec![3_i32, 5]).unwrap();
        let rows = heap.allocate(vec![row0, row1]).unwrap();
        let mut unit = OAD_PolymerUnit {
            nbkbonds: 2,
            bkbonds: rows,
            ..OAD_PolymerUnit::default()
        };
        let mut senior = -1;
        OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(
            &mut heap,
            &mut unit,
            SourceMutPointer::null(),
            &properties(),
            &mut senior,
        )
        .unwrap();
        assert_eq!(senior, 1);
        assert_eq!(&heap.slice(row1.as_const()).unwrap()[..2], &[5, 3]);
        assert_eq!((unit.end_atom1, unit.end_atom2), (5, 3));

        let mut failure_heap = SourceHeap::default();
        let failure_row0 = failure_heap.allocate(vec![1_i32, 2]).unwrap();
        let failure_row1 = failure_heap.allocate(vec![3_i32, 5]).unwrap();
        let failure_rows = failure_heap
            .allocate(vec![failure_row0, failure_row1])
            .unwrap();
        let mut failure_unit = OAD_PolymerUnit {
            nbkbonds: 2,
            bkbonds: failure_rows,
            ..OAD_PolymerUnit::default()
        };
        failure_heap.fail_after_allocations(0);
        senior = -1;
        OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(
            &mut failure_heap,
            &mut failure_unit,
            SourceMutPointer::null(),
            &properties(),
            &mut senior,
        )
        .unwrap();
        assert_eq!(senior, 0);
        assert_eq!(
            &failure_heap.slice(failure_row0.as_const()).unwrap()[..2],
            &[2, 1]
        );
        assert_eq!((failure_unit.end_atom1, failure_unit.end_atom2), (2, 1));

        let mut single_heap = SourceHeap::default();
        let single_row = single_heap.allocate(vec![2_i32, 1]).unwrap();
        let single_rows = single_heap.allocate(vec![single_row]).unwrap();
        let mut single_unit = OAD_PolymerUnit {
            nbkbonds: 1,
            bkbonds: single_rows,
            ..OAD_PolymerUnit::default()
        };
        senior = -1;
        OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(
            &mut single_heap,
            &mut single_unit,
            SourceMutPointer::null(),
            &properties(),
            &mut senior,
        )
        .unwrap();
        assert_eq!(senior, 0);
        assert_eq!(
            &single_heap.slice(single_row.as_const()).unwrap()[..2],
            &[2, 1]
        );
        assert_eq!((single_unit.end_atom1, single_unit.end_atom2), (2, 1));
    }

    #[test]
    fn source_port__runichi3__oad_collectreachableatoms__line_3019() {
        fn chain_atoms(count: usize) -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); count];
            for index in 0..count {
                atoms[index].orig_at_number = (index + 1) as AT_NUMB;
                let mut degree = 0_usize;
                if index != 0 {
                    atoms[index].neighbor[degree] = (index - 1) as AT_NUMB;
                    atoms[index].bond_type[degree] = 1;
                    degree += 1;
                }
                if index + 1 < count {
                    atoms[index].neighbor[degree] = (index + 1) as AT_NUMB;
                    atoms[index].bond_type[degree] = 1;
                    degree += 1;
                }
                atoms[index].valence = degree as i8;
            }
            atoms
        }

        fn atom_data(heap: &mut SourceHeap, count: usize) -> ORIG_ATOM_DATA {
            ORIG_ATOM_DATA {
                at: heap.allocate_model_storage(chain_atoms(count)).unwrap(),
                num_inp_atoms: count as i32,
                num_inp_bonds: count.saturating_sub(1) as i32,
                ..ORIG_ATOM_DATA::default()
            }
        }

        let mut heap = SourceHeap::default();
        let data = atom_data(&mut heap, 5);
        let reachable = heap.allocate_model_storage(vec![-9_i32; 7]).unwrap();
        let baseline = heap.live_allocation_count();
        let mut reachable_count = 5_i32;
        let mut error = 77_i32;
        let mut error_text = [b'o' as i8, b'l' as i8, b'd' as i8, 0, -1];
        heap.trace_source_allocations();
        assert_eq!(
            OAD_CollectReachableAtoms(
                &mut heap,
                &data,
                1,
                0,
                SourceMutPointer::null(),
                &mut reachable_count,
                reachable,
                &mut error,
                Some(&mut error_text),
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(heap.source_allocation_calls(), 13);
        assert_eq!(heap.live_allocation_count(), baseline);
        assert_eq!(reachable_count, 5);
        assert_eq!(
            &heap.slice(reachable.as_const()).unwrap()[..5],
            &[1, 2, 3, 4, 5]
        );
        assert_eq!(&heap.slice(reachable.as_const()).unwrap()[5..], &[-9, -9]);
        assert_eq!(error, 77);
        assert_eq!(error_text, [b'o' as i8, b'l' as i8, b'd' as i8, 0, -1]);

        let forbidden = heap.allocate_model_storage(vec![2_i32, 3, 3, 4]).unwrap();
        let isolated = heap.allocate_model_storage(vec![-7_i32; 3]).unwrap();
        let isolated_baseline = heap.live_allocation_count();
        reachable_count = 3;
        error = -42;
        assert_eq!(
            OAD_CollectReachableAtoms(
                &mut heap,
                &data,
                3,
                2,
                forbidden,
                &mut reachable_count,
                isolated,
                &mut error,
                None,
            ),
            Ok(_IS_OKAY as i32)
        );
        assert_eq!(heap.live_allocation_count(), isolated_baseline);
        assert_eq!(heap.slice(forbidden.as_const()).unwrap(), &[1, 2, 2, 3]);
        assert_eq!(reachable_count, 1);
        assert_eq!(heap.slice(isolated.as_const()).unwrap(), &[3, -7, -7]);
        assert_eq!(error, -42);

        reachable_count = 4;
        error = 91;
        assert_eq!(
            OAD_CollectReachableAtoms(
                &mut heap,
                &data,
                1,
                0,
                SourceMutPointer::null(),
                &mut reachable_count,
                SourceMutPointer::null(),
                &mut error,
                None,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(reachable_count, 0);
        assert_eq!(error, 91);

        let short_output = heap.allocate_model_storage(vec![-3_i32; 2]).unwrap();
        reachable_count = 5;
        assert_eq!(
            OAD_CollectReachableAtoms(
                &mut heap,
                &data,
                1,
                -1,
                SourceMutPointer::null(),
                &mut reachable_count,
                short_output,
                &mut error,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(reachable_count, 2);
        assert_eq!(heap.slice(short_output.as_const()).unwrap(), &[1, 2]);

        let null_forbidden_baseline = heap.live_allocation_count();
        reachable_count = 5;
        assert_eq!(
            OAD_CollectReachableAtoms(
                &mut heap,
                &data,
                1,
                1,
                SourceMutPointer::null(),
                &mut reachable_count,
                reachable,
                &mut error,
                None,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(reachable_count, 0);
        assert_eq!(heap.live_allocation_count(), null_forbidden_baseline);

        let invalid_start_baseline = heap.live_allocation_count();
        reachable_count = 5;
        assert_eq!(
            OAD_CollectReachableAtoms(
                &mut heap,
                &data,
                0,
                0,
                SourceMutPointer::null(),
                &mut reachable_count,
                reachable,
                &mut error,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(reachable_count, 0);
        assert_eq!(heap.live_allocation_count(), invalid_start_baseline);

        reachable_count = 123;
        assert_eq!(
            OAD_CollectReachableAtoms(
                &mut heap,
                &data,
                i32::MIN,
                0,
                SourceMutPointer::null(),
                &mut reachable_count,
                reachable,
                &mut error,
                None,
            ),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
        assert_eq!(reachable_count, 123);

        let negative_data = ORIG_ATOM_DATA {
            num_inp_atoms: -1,
            ..ORIG_ATOM_DATA::default()
        };
        reachable_count = 8;
        assert_eq!(
            OAD_CollectReachableAtoms(
                &mut heap,
                &negative_data,
                1,
                0,
                SourceMutPointer::null(),
                &mut reachable_count,
                SourceMutPointer::null(),
                &mut error,
                None,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(reachable_count, 0);

        for successful_allocations in 0_u64..12 {
            let mut failure_heap = SourceHeap::default();
            let failure_data = atom_data(&mut failure_heap, 4);
            let failure_forbidden = failure_heap.allocate_model_storage(vec![2_i32, 3]).unwrap();
            let failure_reachable = failure_heap
                .allocate_model_storage(vec![-11_i32; 4])
                .unwrap();
            let failure_baseline = failure_heap.live_allocation_count();
            let mut failure_count = 4_i32;
            let mut failure_error = 314_i32;
            let mut failure_text = [9_i8, 8, 7, 0];
            failure_heap.fail_after_allocations(successful_allocations);
            assert_eq!(
                OAD_CollectReachableAtoms(
                    &mut failure_heap,
                    &failure_data,
                    1,
                    1,
                    failure_forbidden,
                    &mut failure_count,
                    failure_reachable,
                    &mut failure_error,
                    Some(&mut failure_text),
                ),
                Ok(_IS_ERROR as i32),
                "allocation ordinal {successful_allocations}"
            );
            assert_eq!(
                failure_count, 0,
                "allocation ordinal {successful_allocations}"
            );
            assert_eq!(
                failure_heap.slice(failure_forbidden.as_const()).unwrap(),
                &[2, 3],
                "allocation ordinal {successful_allocations}"
            );
            assert_eq!(
                failure_heap.slice(failure_reachable.as_const()).unwrap(),
                &[-11, -11, -11, -11],
                "allocation ordinal {successful_allocations}"
            );
            assert_eq!(failure_error, 314);
            assert_eq!(failure_text, [9, 8, 7, 0]);
            assert_eq!(
                failure_heap.live_allocation_count(),
                failure_baseline,
                "allocation ordinal {successful_allocations}"
            );
        }
    }

    #[test]
    fn source_port__runichi3__oad_collectbackbonebonds__line_3109() {
        fn diamond_atoms() -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 4];
            atoms[0].valence = 2;
            atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
            atoms[0].bond_type[..2].copy_from_slice(&[1, 1]);
            atoms[1].valence = 2;
            atoms[1].neighbor[..2].copy_from_slice(&[0, 3]);
            atoms[1].bond_type[..2].copy_from_slice(&[1, 1]);
            atoms[2].valence = 2;
            atoms[2].neighbor[..2].copy_from_slice(&[0, 3]);
            atoms[2].bond_type[..2].copy_from_slice(&[1, 1]);
            atoms[3].valence = 2;
            atoms[3].neighbor[..2].copy_from_slice(&[1, 2]);
            atoms[3].bond_type[..2].copy_from_slice(&[1, 1]);
            atoms
        }

        fn output_rows(
            heap: &mut SourceHeap,
        ) -> (
            SourceMutPointer<SourceMutPointer<i32>>,
            Vec<SourceMutPointer<i32>>,
        ) {
            let rows = (0..6)
                .map(|_| heap.allocate(vec![-1_i32, -1]).unwrap())
                .collect::<Vec<_>>();
            (heap.allocate(rows.clone()).unwrap(), rows)
        }

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(diamond_atoms()).unwrap();
        let atom_data = ORIG_ATOM_DATA {
            at: atoms,
            num_inp_atoms: 4,
            ..ORIG_ATOM_DATA::default()
        };
        let (bonds, rows) = output_rows(&mut heap);
        let baseline = heap.live_allocation_count();
        let mut count = 99_i32;
        let mut error = 77_i32;
        let mut error_text = [0_i8; 256];
        OAD_CollectBackboneBonds(
            &mut heap,
            &atom_data,
            4,
            &[1, 2, 3, 4],
            1,
            4,
            &mut count,
            bonds,
            &mut error,
            Some(&mut error_text),
        )
        .unwrap();
        assert_eq!(count, 4);
        assert_eq!(error, 0);
        assert_eq!(error_text[0], 0);
        assert_eq!(heap.live_allocation_count(), baseline);
        let actual = rows[..4]
            .iter()
            .map(|row| heap.slice(row.as_const()).unwrap().to_vec())
            .collect::<Vec<_>>();
        assert_eq!(actual, vec![vec![1, 2], vec![2, 4], vec![1, 3], vec![3, 4]]);

        let (reverse_bonds, reverse_rows) = output_rows(&mut heap);
        let reverse_baseline = heap.live_allocation_count();
        count = -8;
        error = -9;
        OAD_CollectBackboneBonds(
            &mut heap,
            &atom_data,
            4,
            &[1, 2, 3, 4],
            4,
            1,
            &mut count,
            reverse_bonds,
            &mut error,
            None,
        )
        .unwrap();
        assert_eq!(count, 4);
        assert_eq!(error, 0);
        assert_eq!(heap.live_allocation_count(), reverse_baseline);
        assert_eq!(heap.slice(reverse_rows[0].as_const()).unwrap(), &[4, 2]);
        assert_eq!(heap.slice(reverse_rows[1].as_const()).unwrap(), &[2, 1]);

        let mut graph_failure_heap = SourceHeap::default();
        let failure_atoms = graph_failure_heap
            .allocate_model_storage(diamond_atoms())
            .unwrap();
        let failure_data = ORIG_ATOM_DATA {
            at: failure_atoms,
            num_inp_atoms: 4,
            ..ORIG_ATOM_DATA::default()
        };
        let (failure_bonds, _) = output_rows(&mut graph_failure_heap);
        let failure_baseline = graph_failure_heap.live_allocation_count();
        graph_failure_heap.fail_after_allocations(0);
        let mut failure_count = 123_i32;
        let mut failure_error = 0_i32;
        let mut failure_text = [0_i8; 256];
        OAD_CollectBackboneBonds(
            &mut graph_failure_heap,
            &failure_data,
            4,
            &[1, 2, 3, 4],
            1,
            4,
            &mut failure_count,
            failure_bonds,
            &mut failure_error,
            Some(&mut failure_text),
        )
        .unwrap();
        assert_eq!(failure_count, 0);
        assert_eq!(failure_error, 9037);
        let expected = b"Not enough memory (polymers)\0";
        assert_eq!(
            &failure_text[..expected.len()],
            &expected.map(|byte| byte as i8)
        );
        assert_eq!(graph_failure_heap.live_allocation_count(), failure_baseline);

        for successful_allocations in [9_u64, 10] {
            let mut pathfinder_failure_heap = SourceHeap::default();
            let pathfinder_failure_atoms = pathfinder_failure_heap
                .allocate_model_storage(diamond_atoms())
                .unwrap();
            let pathfinder_failure_data = ORIG_ATOM_DATA {
                at: pathfinder_failure_atoms,
                num_inp_atoms: 4,
                ..ORIG_ATOM_DATA::default()
            };
            let (pathfinder_failure_bonds, _) = output_rows(&mut pathfinder_failure_heap);
            let pathfinder_failure_baseline = pathfinder_failure_heap.live_allocation_count();
            pathfinder_failure_heap.fail_after_allocations(successful_allocations);
            let mut pathfinder_failure_count = 55_i32;
            let mut pathfinder_failure_error = 41_i32;
            let mut pathfinder_failure_text = [0_i8; 256];
            pathfinder_failure_text[..4].copy_from_slice(&[b'o' as i8, b'l' as i8, b'd' as i8, 0]);
            OAD_CollectBackboneBonds(
                &mut pathfinder_failure_heap,
                &pathfinder_failure_data,
                4,
                &[1, 2, 3, 4],
                1,
                4,
                &mut pathfinder_failure_count,
                pathfinder_failure_bonds,
                &mut pathfinder_failure_error,
                Some(&mut pathfinder_failure_text),
            )
            .unwrap();
            assert_eq!(pathfinder_failure_count, 0);
            assert_eq!(pathfinder_failure_error, 41);
            let expected = b"old; Not enough memory (polymers)\0";
            assert_eq!(
                &pathfinder_failure_text[..expected.len()],
                &expected.map(|byte| byte as i8),
                "allocation ordinal {successful_allocations}"
            );
            assert_eq!(
                pathfinder_failure_heap.live_allocation_count(),
                pathfinder_failure_baseline + 9,
                "the source leaks sg when pathfinder allocation fails at ordinal {successful_allocations}"
            );
        }
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_removelinkfromcruchain__line_3607() {
        fn make_bonds(
            heap: &mut SourceHeap,
            values: &[[i32; 2]],
        ) -> (
            SourceMutPointer<SourceMutPointer<i32>>,
            Vec<SourceMutPointer<i32>>,
        ) {
            let rows = values
                .iter()
                .map(|value| heap.allocate(value.to_vec()).unwrap())
                .collect::<Vec<_>>();
            (heap.allocate(rows.clone()).unwrap(), rows)
        }

        let mut heap = SourceHeap::default();
        let (bonds, rows) = make_bonds(&mut heap, &[[1, 2], [3, 4], [5, 6]]);
        let mut count = 3_i32;
        OAD_PolymerUnit_RemoveLinkFromCRUChain(&mut heap, 3, 4, &mut count, bonds).unwrap();
        assert_eq!(count, 2);
        assert_eq!(heap.slice(rows[0].as_const()).unwrap(), &[1, 2]);
        assert_eq!(heap.slice(rows[1].as_const()).unwrap(), &[5, 6]);
        assert_eq!(heap.slice(rows[2].as_const()).unwrap(), &[5, 6]);

        let (first_bonds, first_rows) =
            make_bonds(&mut heap, &[[i32::MIN, i32::MAX], [7, 8], [9, 10]]);
        let mut first_count = 3_i32;
        OAD_PolymerUnit_RemoveLinkFromCRUChain(
            &mut heap,
            i32::MIN,
            i32::MAX,
            &mut first_count,
            first_bonds,
        )
        .unwrap();
        assert_eq!(first_count, 2);
        assert_eq!(heap.slice(first_rows[0].as_const()).unwrap(), &[7, 8]);
        assert_eq!(heap.slice(first_rows[1].as_const()).unwrap(), &[9, 10]);
        assert_eq!(heap.slice(first_rows[2].as_const()).unwrap(), &[9, 10]);

        let (last_bonds, last_rows) = make_bonds(&mut heap, &[[1, 2], [3, 4]]);
        let mut last_count = 2_i32;
        OAD_PolymerUnit_RemoveLinkFromCRUChain(&mut heap, 3, 4, &mut last_count, last_bonds)
            .unwrap();
        assert_eq!(last_count, 1);
        assert_eq!(heap.slice(last_rows[1].as_const()).unwrap(), &[3, 4]);

        let (reverse_bonds, reverse_rows) = make_bonds(&mut heap, &[[2, 1], [3, 4]]);
        let mut reverse_count = 2_i32;
        OAD_PolymerUnit_RemoveLinkFromCRUChain(&mut heap, 1, 2, &mut reverse_count, reverse_bonds)
            .unwrap();
        assert_eq!(reverse_count, 2);
        assert_eq!(heap.slice(reverse_rows[0].as_const()).unwrap(), &[2, 1]);

        let (duplicate_bonds, duplicate_rows) = make_bonds(&mut heap, &[[4, 5], [4, 5], [6, 7]]);
        let mut duplicate_count = 3_i32;
        OAD_PolymerUnit_RemoveLinkFromCRUChain(
            &mut heap,
            4,
            5,
            &mut duplicate_count,
            duplicate_bonds,
        )
        .unwrap();
        assert_eq!(duplicate_count, 2);
        assert_eq!(heap.slice(duplicate_rows[0].as_const()).unwrap(), &[4, 5]);
        assert_eq!(heap.slice(duplicate_rows[1].as_const()).unwrap(), &[6, 7]);

        let mut negative_count = i32::MIN;
        OAD_PolymerUnit_RemoveLinkFromCRUChain(
            &mut heap,
            1,
            2,
            &mut negative_count,
            SourceMutPointer::null(),
        )
        .unwrap();
        assert_eq!(negative_count, i32::MIN);
        let mut zero_count = 0_i32;
        OAD_PolymerUnit_RemoveLinkFromCRUChain(
            &mut heap,
            1,
            2,
            &mut zero_count,
            SourceMutPointer::null(),
        )
        .unwrap();
        assert_eq!(zero_count, 0);
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_delistintraringbackbonebonds__line_3168() {
        fn ring_atom_data(heap: &mut SourceHeap) -> ORIG_ATOM_DATA {
            let mut atoms = vec![inp_ATOM::default(); 4];
            atoms[0].valence = 2;
            atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
            atoms[1].valence = 2;
            atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
            atoms[2].valence = 2;
            atoms[2].neighbor[..2].copy_from_slice(&[0, 1]);
            for (index, atom) in atoms.iter_mut().enumerate() {
                atom.orig_at_number = u16::try_from(index + 1).unwrap();
                atom.bond_type[..usize::try_from(atom.valence).unwrap()].fill(1);
                atom.chem_bonds_valence = atom.valence;
            }
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer::default()])
                .unwrap();
            ORIG_ATOM_DATA {
                at: atoms,
                polymer,
                num_inp_atoms: 4,
                num_inp_bonds: 3,
                ..ORIG_ATOM_DATA::default()
            }
        }

        fn make_unit(
            heap: &mut SourceHeap,
            values: &[[i32; 2]],
        ) -> (
            SourceMutPointer<OAD_PolymerUnit>,
            Vec<SourceMutPointer<i32>>,
        ) {
            let rows = values
                .iter()
                .map(|value| heap.allocate(value.to_vec()).unwrap())
                .collect::<Vec<_>>();
            let bonds = heap.allocate(rows.clone()).unwrap();
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    end_atom1: 1,
                    end_atom2: 4,
                    nbkbonds: i32::try_from(values.len()).unwrap(),
                    bkbonds: bonds,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            (unit, rows)
        }

        let mut heap = SourceHeap::default();
        let mut atom_data = ring_atom_data(&mut heap);
        let (unit, rows) = make_unit(&mut heap, &[[1, 2], [3, 4], [2, 3]]);
        let baseline = heap.live_allocation_count();
        let mut error = 88_i32;
        let mut text = [0_i8; 256];
        OAD_PolymerUnit_DelistIntraRingBackboneBonds(
            &mut heap,
            unit,
            &mut atom_data,
            &mut error,
            Some(&mut text),
        )
        .unwrap();
        assert_eq!(error, 0);
        assert_eq!(text[0], 0);
        assert_eq!(heap.live_allocation_count(), baseline);
        assert_eq!(heap.slice(unit.as_const()).unwrap()[0].nbkbonds, 1);
        assert_eq!(heap.slice(rows[0].as_const()).unwrap(), &[3, 4]);
        assert_eq!(heap.slice(rows[1].as_const()).unwrap(), &[2, 3]);
        assert_eq!(heap.slice(rows[2].as_const()).unwrap(), &[2, 3]);

        let mut chain_atoms = vec![inp_ATOM::default(); 3];
        chain_atoms[0].valence = 1;
        chain_atoms[0].neighbor[0] = 1;
        chain_atoms[1].valence = 2;
        chain_atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
        chain_atoms[2].valence = 1;
        chain_atoms[2].neighbor[0] = 1;
        for (index, atom) in chain_atoms.iter_mut().enumerate() {
            atom.orig_at_number = u16::try_from(index + 1).unwrap();
            atom.bond_type[..usize::try_from(atom.valence).unwrap()].fill(1);
            atom.chem_bonds_valence = atom.valence;
        }
        let chain_atom_pointer = heap.allocate_model_storage(chain_atoms).unwrap();
        let chain_polymer = heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        let mut chain_data = ORIG_ATOM_DATA {
            at: chain_atom_pointer,
            polymer: chain_polymer,
            num_inp_atoms: 3,
            num_inp_bonds: 2,
            ..ORIG_ATOM_DATA::default()
        };
        let (chain_unit, chain_rows) = make_unit(&mut heap, &[[1, 2], [2, 3]]);
        error = 66;
        OAD_PolymerUnit_DelistIntraRingBackboneBonds(
            &mut heap,
            chain_unit,
            &mut chain_data,
            &mut error,
            None,
        )
        .unwrap();
        assert_eq!(error, 0);
        assert_eq!(heap.slice(chain_unit.as_const()).unwrap()[0].nbkbonds, 2);
        assert_eq!(heap.slice(chain_rows[0].as_const()).unwrap(), &[1, 2]);

        let mut early_data = ORIG_ATOM_DATA::default();
        error = 19;
        OAD_PolymerUnit_DelistIntraRingBackboneBonds(
            &mut heap,
            SourceMutPointer::null(),
            &mut early_data,
            &mut error,
            None,
        )
        .unwrap();
        assert_eq!(error, 19);
        let empty_unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                nbkbonds: i32::MIN,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        OAD_PolymerUnit_DelistIntraRingBackboneBonds(
            &mut heap,
            empty_unit,
            &mut early_data,
            &mut error,
            None,
        )
        .unwrap();
        assert_eq!(error, 19);

        let mut failure_heap = SourceHeap::default();
        let mut failure_data = ring_atom_data(&mut failure_heap);
        let (failure_unit, failure_rows) = make_unit(&mut failure_heap, &[[1, 2], [3, 4]]);
        let failure_baseline = failure_heap.live_allocation_count();
        failure_heap.fail_after_allocations(0);
        let mut failure_error = 72_i32;
        let mut failure_text = [b'x' as i8, 0_i8]
            .into_iter()
            .chain(std::iter::repeat_n(0_i8, 254))
            .collect::<Vec<_>>();
        OAD_PolymerUnit_DelistIntraRingBackboneBonds(
            &mut failure_heap,
            failure_unit,
            &mut failure_data,
            &mut failure_error,
            Some(&mut failure_text),
        )
        .unwrap();
        assert_eq!(failure_error, 1);
        assert_eq!(&failure_text[..2], &[b'x' as i8, 0]);
        assert_eq!(failure_heap.live_allocation_count(), failure_baseline);
        assert_eq!(
            failure_heap.slice(failure_unit.as_const()).unwrap()[0].nbkbonds,
            2
        );
        assert_eq!(
            failure_heap.slice(failure_rows[0].as_const()).unwrap(),
            &[1, 2]
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_delisthighorderbackbonebonds__line_3509() {
        fn make_unit(
            heap: &mut SourceHeap,
            values: &[[i32; 2]],
        ) -> (
            SourceMutPointer<OAD_PolymerUnit>,
            Vec<SourceMutPointer<i32>>,
        ) {
            let rows = values
                .iter()
                .map(|value| heap.allocate(value.to_vec()).unwrap())
                .collect::<Vec<_>>();
            let bonds = heap.allocate(rows.clone()).unwrap();
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    na: 3,
                    nb: 2,
                    nbkbonds: i32::try_from(values.len()).unwrap(),
                    bkbonds: bonds,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            (unit, rows)
        }

        fn original_data(heap: &mut SourceHeap) -> ORIG_ATOM_DATA {
            ORIG_ATOM_DATA {
                at: heap
                    .allocate_model_storage(vec![inp_ATOM::default(); 3])
                    .unwrap(),
                num_inp_atoms: 3,
                ..ORIG_ATOM_DATA::default()
            }
        }

        fn mapped_composite(heap: &mut SourceHeap) -> COMP_ATOM_DATA {
            let mut atoms = vec![inp_ATOM::default(); 3];
            atoms[0].orig_at_number = 1;
            atoms[0].valence = 1;
            atoms[0].neighbor[0] = 1;
            atoms[0].bond_type[0] = BOND_TAUTOM as u8;
            atoms[1].orig_at_number = 2;
            atoms[1].valence = 2;
            atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
            atoms[1].bond_type[..2].copy_from_slice(&[BOND_TAUTOM as u8, 1]);
            atoms[2].orig_at_number = 3;
            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 1;
            atoms[2].bond_type[0] = 1;
            COMP_ATOM_DATA {
                at: heap.allocate_model_storage(atoms).unwrap(),
                num_at: 3,
                ..COMP_ATOM_DATA::default()
            }
        }

        let mut heap = SourceHeap::default();
        let original = original_data(&mut heap);
        let composite = mapped_composite(&mut heap);
        let (unit, rows) = make_unit(&mut heap, &[[1, 2], [2, 3]]);
        let baseline = heap.live_allocation_count();
        let mut error = 71_i32;
        let mut text = [0_i8; 256];
        text[..2].copy_from_slice(&[b'x' as i8, 0]);
        OAD_PolymerUnit_DelistHighOrderBackboneBonds(
            &mut heap,
            unit,
            &original,
            Some(&composite),
            &mut error,
            Some(&mut text),
        )
        .unwrap();
        assert_eq!(heap.slice(unit.as_const()).unwrap()[0].nbkbonds, 1);
        assert_eq!(heap.slice(rows[0].as_const()).unwrap(), &[2, 3]);
        assert_eq!(heap.slice(rows[1].as_const()).unwrap(), &[2, 3]);
        assert_eq!(heap.live_allocation_count(), baseline);
        assert_eq!(error, 71);
        assert_eq!(&text[..2], &[b'x' as i8, 0]);

        let (no_composite_unit, no_composite_rows) = make_unit(&mut heap, &[[1, 2]]);
        let original_atoms = heap.slice_mut(original.at).unwrap();
        original_atoms[0].valence = 1;
        original_atoms[0].neighbor[0] = 1;
        original_atoms[0].bond_type[0] = 3;
        OAD_PolymerUnit_DelistHighOrderBackboneBonds(
            &mut heap,
            no_composite_unit,
            &original,
            None,
            &mut error,
            None,
        )
        .unwrap();
        assert_eq!(
            heap.slice(no_composite_unit.as_const()).unwrap()[0].nbkbonds,
            1
        );
        assert_eq!(
            heap.slice(no_composite_rows[0].as_const()).unwrap(),
            &[1, 2]
        );

        let mut second_failure_heap = SourceHeap::default();
        let second_failure_original = original_data(&mut second_failure_heap);
        let second_failure_composite = mapped_composite(&mut second_failure_heap);
        let (second_failure_unit, _) = make_unit(&mut second_failure_heap, &[[1, 2]]);
        let second_failure_baseline = second_failure_heap.live_allocation_count();
        second_failure_heap.fail_after_allocations(1);
        let mut second_failure_error = 82_i32;
        OAD_PolymerUnit_DelistHighOrderBackboneBonds(
            &mut second_failure_heap,
            second_failure_unit,
            &second_failure_original,
            Some(&second_failure_composite),
            &mut second_failure_error,
            None,
        )
        .unwrap();
        assert_eq!(
            second_failure_heap
                .slice(second_failure_unit.as_const())
                .unwrap()[0]
                .nbkbonds,
            1
        );
        assert_eq!(
            second_failure_heap.live_allocation_count(),
            second_failure_baseline
        );
        assert_eq!(second_failure_error, 82);

        let mut first_failure_heap = SourceHeap::default();
        let first_failure_original = original_data(&mut first_failure_heap);
        let mut zero_mapped_atoms = vec![inp_ATOM::default(); 2];
        zero_mapped_atoms[0].orig_at_number = 1;
        zero_mapped_atoms[0].valence = 1;
        zero_mapped_atoms[0].neighbor[0] = 0;
        zero_mapped_atoms[0].bond_type[0] = BOND_TAUTOM as u8;
        zero_mapped_atoms[1].orig_at_number = 2;
        let first_failure_composite = COMP_ATOM_DATA {
            at: first_failure_heap
                .allocate_model_storage(zero_mapped_atoms)
                .unwrap(),
            num_at: 2,
            ..COMP_ATOM_DATA::default()
        };
        let (first_failure_unit, _) = make_unit(&mut first_failure_heap, &[[1, 2]]);
        let first_failure_baseline = first_failure_heap.live_allocation_count();
        first_failure_heap.fail_after_allocations(0);
        let mut first_failure_error = 93_i32;
        OAD_PolymerUnit_DelistHighOrderBackboneBonds(
            &mut first_failure_heap,
            first_failure_unit,
            &first_failure_original,
            Some(&first_failure_composite),
            &mut first_failure_error,
            None,
        )
        .unwrap();
        assert_eq!(
            first_failure_heap
                .slice(first_failure_unit.as_const())
                .unwrap()[0]
                .nbkbonds,
            0
        );
        assert_eq!(
            first_failure_heap.live_allocation_count(),
            first_failure_baseline
        );
        assert_eq!(first_failure_error, 93);

        for (na, nb, nbkbonds) in [(1, 9, 1), (9, 1, 1), (9, 9, 0)] {
            let early = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    na,
                    nb,
                    nbkbonds,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            OAD_PolymerUnit_DelistHighOrderBackboneBonds(
                &mut heap,
                early,
                &original,
                Some(&composite),
                &mut error,
                None,
            )
            .unwrap();
            assert_eq!(heap.slice(early.as_const()).unwrap()[0].nbkbonds, nbkbonds);
        }
    }

    #[test]
    fn source_port__runichi3__oad_polymer_findbackbones__line_2872() {
        fn triangle_data(
            heap: &mut SourceHeap,
            include_skipped: bool,
        ) -> (
            ORIG_ATOM_DATA,
            SourceMutPointer<OAD_PolymerUnit>,
            Vec<SourceMutPointer<i32>>,
            Option<(SourceMutPointer<OAD_PolymerUnit>, SourceMutPointer<i32>)>,
        ) {
            let mut atoms = vec![inp_ATOM::default(); 3];
            atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
            atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
            atoms[2].neighbor[..2].copy_from_slice(&[0, 1]);
            for (index, atom) in atoms.iter_mut().enumerate() {
                atom.valence = 2;
                atom.chem_bonds_valence = 2;
                atom.bond_type[..2].copy_from_slice(&[1, 1]);
                atom.orig_at_number = u16::try_from(index + 1).unwrap();
            }
            let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
            let atom_list = heap.allocate(vec![1_i32, 2, 3]).unwrap();
            let rows = (0..6)
                .map(|_| heap.allocate(vec![-1_i32, -1]).unwrap())
                .collect::<Vec<_>>();
            let bonds = heap.allocate(rows.clone()).unwrap();
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    na: 3,
                    nb: 3,
                    cyclizable: 1,
                    end_atom1: 1,
                    end_atom2: 3,
                    alist: atom_list,
                    maxbkbonds: 6,
                    nbkbonds: 55,
                    bkbonds: bonds,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            let mut unit_pointers = vec![unit];
            let skipped = if include_skipped {
                let skipped_row = heap.allocate(vec![77_i32, 88]).unwrap();
                let skipped_rows = heap.allocate(vec![skipped_row]).unwrap();
                let skipped_unit = heap
                    .allocate_model_storage(vec![OAD_PolymerUnit {
                        cyclizable: 0,
                        nbkbonds: 1,
                        bkbonds: skipped_rows,
                        ..OAD_PolymerUnit::default()
                    }])
                    .unwrap();
                unit_pointers.push(skipped_unit);
                Some((skipped_unit, skipped_row))
            } else {
                None
            };
            let units = heap.allocate_model_storage(unit_pointers).unwrap();
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer {
                    units,
                    n: if include_skipped { 2 } else { 1 },
                    ..OAD_Polymer::default()
                }])
                .unwrap();
            (
                ORIG_ATOM_DATA {
                    at: atom_pointer,
                    num_inp_atoms: 3,
                    num_inp_bonds: 3,
                    polymer,
                    ..ORIG_ATOM_DATA::default()
                },
                unit,
                rows,
                skipped,
            )
        }

        let mut heap = SourceHeap::default();
        let (mut atom_data, unit, rows, skipped) = triangle_data(&mut heap, true);
        let baseline = heap.live_allocation_count();
        let mut error = 91_i32;
        let mut text = [0_i8; 256];
        OAD_Polymer_FindBackbones(&mut heap, &mut atom_data, None, &mut error, Some(&mut text))
            .unwrap();
        let result = heap.slice(unit.as_const()).unwrap()[0].clone();
        assert_eq!(error, 0);
        assert_eq!(result.cyclizable, 1);
        assert_eq!(result.nbkbonds, 1);
        assert_eq!(heap.slice(rows[0].as_const()).unwrap(), &[1, 3]);
        assert_eq!(heap.live_allocation_count(), baseline);
        let (skipped_unit, skipped_row) = skipped.unwrap();
        assert_eq!(heap.slice(skipped_unit.as_const()).unwrap()[0].nbkbonds, 1);
        assert_eq!(heap.slice(skipped_row.as_const()).unwrap(), &[77, 88]);

        let mut collect_failure_heap = SourceHeap::default();
        let (mut collect_failure_data, collect_failure_unit, collect_failure_rows, _) =
            triangle_data(&mut collect_failure_heap, false);
        let collect_failure_baseline = collect_failure_heap.live_allocation_count();
        collect_failure_heap.fail_after_allocations(0);
        let mut collect_failure_error = 17_i32;
        let mut collect_failure_text = [0_i8; 256];
        OAD_Polymer_FindBackbones(
            &mut collect_failure_heap,
            &mut collect_failure_data,
            None,
            &mut collect_failure_error,
            Some(&mut collect_failure_text),
        )
        .unwrap();
        let failed = collect_failure_heap
            .slice(collect_failure_unit.as_const())
            .unwrap()[0]
            .clone();
        assert_eq!(failed.cyclizable, CLOSING_SRU_NOT_APPLICABLE as i32);
        assert_eq!(failed.nbkbonds, 0);
        assert_eq!(collect_failure_error, 9037);
        let expected = b"Not enough memory (polymers)\0";
        assert_eq!(
            &collect_failure_text[..expected.len()],
            &expected.map(|byte| byte as i8)
        );
        assert_eq!(
            collect_failure_heap.live_allocation_count(),
            collect_failure_baseline
        );
        assert_eq!(
            collect_failure_heap
                .slice(collect_failure_rows[0].as_const())
                .unwrap(),
            &[-1, -1]
        );

        let mut intra_failure_heap = SourceHeap::default();
        let (mut intra_failure_data, intra_failure_unit, _, _) =
            triangle_data(&mut intra_failure_heap, false);
        let intra_failure_baseline = intra_failure_heap.live_allocation_count();
        intra_failure_heap.fail_after_allocations(10);
        let mut intra_failure_error = 33_i32;
        OAD_Polymer_FindBackbones(
            &mut intra_failure_heap,
            &mut intra_failure_data,
            None,
            &mut intra_failure_error,
            None,
        )
        .unwrap();
        let partial = intra_failure_heap
            .slice(intra_failure_unit.as_const())
            .unwrap()[0]
            .clone();
        assert_eq!(intra_failure_error, 1);
        assert_eq!(partial.cyclizable, 1);
        assert!(partial.nbkbonds > 1);
        assert_eq!(
            intra_failure_heap.live_allocation_count(),
            intra_failure_baseline
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymer_getrepresentation__line_3788() {
        fn polymer_with_units(
            heap: &mut SourceHeap,
            units: Vec<OAD_PolymerUnit>,
        ) -> (
            SourceMutPointer<OAD_Polymer>,
            Vec<SourceMutPointer<OAD_PolymerUnit>>,
        ) {
            let pointers = units
                .into_iter()
                .map(|unit| heap.allocate_model_storage(vec![unit]).unwrap())
                .collect::<Vec<_>>();
            let unit_array = heap.allocate_model_storage(pointers.clone()).unwrap();
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer {
                    units: unit_array,
                    n: i32::try_from(pointers.len()).unwrap(),
                    ..OAD_Polymer::default()
                }])
                .unwrap();
            (polymer, pointers)
        }

        let mut heap = SourceHeap::default();
        assert_eq!(
            OAD_Polymer_GetRepresentation(&mut heap, SourceMutPointer::null()),
            Ok(NO_POLYMER)
        );
        let empty = heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        assert_eq!(
            OAD_Polymer_GetRepresentation(&mut heap, empty),
            Ok(POLYMER_REPRESENTATION_SOURCE_BASED as i32)
        );
        let negative = heap
            .allocate_model_storage(vec![OAD_Polymer {
                n: i32::MIN,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        assert_eq!(
            OAD_Polymer_GetRepresentation(&mut heap, negative),
            Ok(POLYMER_REPRESENTATION_UNRECOGNIZED as i32)
        );

        let (source, source_units) = polymer_with_units(
            &mut heap,
            vec![OAD_PolymerUnit {
                nb: 0,
                representation: 77,
                ..OAD_PolymerUnit::default()
            }],
        );
        assert_eq!(
            OAD_Polymer_GetRepresentation(&mut heap, source),
            Ok(POLYMER_REPRESENTATION_SOURCE_BASED as i32)
        );
        assert_eq!(
            heap.slice(source_units[0].as_const()).unwrap()[0].representation,
            POLYMER_REPRESENTATION_SOURCE_BASED as i32
        );

        let (structure, structure_units) = polymer_with_units(
            &mut heap,
            vec![
                OAD_PolymerUnit {
                    nb: 2,
                    ..OAD_PolymerUnit::default()
                },
                OAD_PolymerUnit {
                    nb: 0,
                    nbkbonds: 1,
                    ..OAD_PolymerUnit::default()
                },
                OAD_PolymerUnit {
                    nb: 0,
                    cap1: 1,
                    cap2: i32::MAX,
                    ..OAD_PolymerUnit::default()
                },
            ],
        );
        assert_eq!(
            OAD_Polymer_GetRepresentation(&mut heap, structure),
            Ok(POLYMER_REPRESENTATION_STRUCTURE_BASED as i32)
        );
        for unit in structure_units {
            assert_eq!(
                heap.slice(unit.as_const()).unwrap()[0].representation,
                POLYMER_REPRESENTATION_STRUCTURE_BASED as i32
            );
        }

        let (mixed, _) = polymer_with_units(
            &mut heap,
            vec![
                OAD_PolymerUnit::default(),
                OAD_PolymerUnit {
                    nb: 2,
                    ..OAD_PolymerUnit::default()
                },
            ],
        );
        assert_eq!(
            OAD_Polymer_GetRepresentation(&mut heap, mixed),
            Ok(POLYMER_REPRESENTATION_MIXED as i32)
        );

        let (unrecognized, unrecognized_units) = polymer_with_units(
            &mut heap,
            vec![OAD_PolymerUnit {
                nb: 1,
                cap1: 1,
                cap2: 0,
                representation: 93,
                ..OAD_PolymerUnit::default()
            }],
        );
        assert_eq!(
            OAD_Polymer_GetRepresentation(&mut heap, unrecognized),
            Ok(POLYMER_REPRESENTATION_UNRECOGNIZED as i32)
        );
        assert_eq!(
            heap.slice(unrecognized_units[0].as_const()).unwrap()[0].representation,
            93
        );
    }

    #[test]
    fn source_port__runichi3__origatdata_increasebondorder__line_2692() {
        fn pair(
            heap: &mut SourceHeap,
            first_type: u8,
            second_type: u8,
        ) -> SourceMutPointer<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 2];
            atoms[0].valence = 1;
            atoms[0].neighbor[0] = 1;
            atoms[0].bond_type[0] = first_type;
            atoms[0].chem_bonds_valence = first_type as i8;
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].bond_type[0] = second_type;
            atoms[1].chem_bonds_valence = second_type as i8;
            heap.allocate_model_storage(atoms).unwrap()
        }

        let mut heap = SourceHeap::default();
        for initial in 0_u8..=3 {
            let normal = pair(&mut heap, initial, initial);
            assert_eq!(OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, normal), Ok(2));
            let values = heap.slice(normal.as_const()).unwrap();
            assert_eq!(
                (values[0].bond_type[0], values[0].chem_bonds_valence),
                (initial + 1, initial as i8 + 1)
            );
            assert_eq!(
                (values[1].bond_type[0], values[1].chem_bonds_valence),
                (initial + 1, initial as i8 + 1)
            );
        }

        for limited_index in [0_usize, 1] {
            let limited = pair(&mut heap, 1, 1);
            heap.slice_mut(limited).unwrap()[limited_index].valence = MAXVAL as i8;
            let before = heap.slice(limited.as_const()).unwrap().to_vec();
            assert_eq!(
                OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, limited),
                Ok(0)
            );
            assert_eq!(heap.slice(limited.as_const()).unwrap(), before);
        }

        for negative_index in [0_usize, 1] {
            let negative = pair(&mut heap, 1, 1);
            heap.slice_mut(negative).unwrap()[negative_index].valence = -1;
            assert_eq!(
                OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, negative),
                Ok(1)
            );
            let values = heap.slice(negative.as_const()).unwrap();
            if negative_index == 0 {
                assert_eq!((values[0].bond_type[0], values[1].bond_type[0]), (1, 2));
            } else {
                assert_eq!((values[0].bond_type[0], values[1].bond_type[0]), (2, 1));
            }
        }

        let first_valence_max = pair(&mut heap, 1, 1);
        heap.slice_mut(first_valence_max).unwrap()[0].chem_bonds_valence = MAXVAL as i8;
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, first_valence_max),
            Ok(0)
        );
        assert_eq!(
            heap.slice(first_valence_max.as_const()).unwrap()[0].bond_type[0],
            1
        );

        let second_valence_max = pair(&mut heap, 1, 1);
        heap.slice_mut(second_valence_max).unwrap()[1].chem_bonds_valence = MAXVAL as i8;
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, second_valence_max),
            Ok(0)
        );
        let values = heap.slice(second_valence_max.as_const()).unwrap();
        assert_eq!(
            (values[0].bond_type[0], values[0].chem_bonds_valence),
            (2, 2)
        );
        assert_eq!(
            (values[1].bond_type[0], values[1].chem_bonds_valence),
            (1, MAXVAL as i8)
        );

        let first_unknown = pair(&mut heap, 4, 1);
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, first_unknown),
            Ok(0)
        );
        let values = heap.slice(first_unknown.as_const()).unwrap();
        assert_eq!((values[0].bond_type[0], values[1].bond_type[0]), (4, 1));

        let second_unknown = pair(&mut heap, 1, 4);
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, second_unknown),
            Ok(0)
        );
        let values = heap.slice(second_unknown.as_const()).unwrap();
        assert_eq!((values[0].bond_type[0], values[1].bond_type[0]), (2, 4));
        assert_eq!(
            (values[0].chem_bonds_valence, values[1].chem_bonds_valence),
            (2, 4)
        );

        let first_missing = pair(&mut heap, 1, 1);
        heap.slice_mut(first_missing).unwrap()[0].neighbor[0] = 0;
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, first_missing),
            Ok(1)
        );
        let values = heap.slice(first_missing.as_const()).unwrap();
        assert_eq!((values[0].bond_type[0], values[1].bond_type[0]), (1, 2));

        let second_missing = pair(&mut heap, 1, 1);
        heap.slice_mut(second_missing).unwrap()[1].neighbor[0] = 1;
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, second_missing),
            Ok(1)
        );
        let values = heap.slice(second_missing.as_const()).unwrap();
        assert_eq!((values[0].bond_type[0], values[1].bond_type[0]), (2, 1));

        let both_missing = pair(&mut heap, 1, 1);
        {
            let values = heap.slice_mut(both_missing).unwrap();
            values[0].neighbor[0] = 0;
            values[1].neighbor[0] = 1;
        }
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, both_missing),
            Ok(0)
        );

        let mut scanned_atoms = vec![inp_ATOM::default(); 3];
        scanned_atoms[0].valence = 2;
        scanned_atoms[0].neighbor[..2].copy_from_slice(&[2, 1]);
        scanned_atoms[0].bond_type[..2].copy_from_slice(&[3, 1]);
        scanned_atoms[0].chem_bonds_valence = 4;
        scanned_atoms[1].valence = 2;
        scanned_atoms[1].neighbor[..2].copy_from_slice(&[2, 0]);
        scanned_atoms[1].bond_type[..2].copy_from_slice(&[3, 1]);
        scanned_atoms[1].chem_bonds_valence = 4;
        let scanned = heap.allocate_model_storage(scanned_atoms).unwrap();
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, 0, 1, scanned),
            Ok(2)
        );
        let values = heap.slice(scanned.as_const()).unwrap();
        assert_eq!(values[0].bond_type[..2], [3, 2]);
        assert_eq!(values[1].bond_type[..2], [3, 2]);
        assert_eq!(
            (values[0].chem_bonds_valence, values[1].chem_bonds_valence),
            (5, 5)
        );

        let mut self_atom = inp_ATOM {
            valence: 1,
            chem_bonds_valence: 1,
            ..inp_ATOM::default()
        };
        self_atom.neighbor[0] = 0;
        self_atom.bond_type[0] = 1;
        let self_pointer = heap.allocate_model_storage(vec![self_atom]).unwrap();
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, 0, 0, self_pointer),
            Ok(2)
        );
        let self_value = &heap.slice(self_pointer.as_const()).unwrap()[0];
        assert_eq!(
            (self_value.bond_type[0], self_value.chem_bonds_valence),
            (3, 3)
        );

        let normal = pair(&mut heap, 1, 1);
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, -1, 1, normal),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            OrigAtData_IncreaseBondOrder(&mut heap, 0, 2, normal),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__runichi3__origatdata_decreasebondorder__line_2750() {
        fn pair(
            heap: &mut SourceHeap,
            first_type: u8,
            second_type: u8,
        ) -> SourceMutPointer<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 2];
            atoms[0].valence = 1;
            atoms[0].neighbor[0] = 1;
            atoms[0].bond_type[0] = first_type;
            atoms[0].chem_bonds_valence = first_type as i8;
            atoms[1].valence = 1;
            atoms[1].neighbor[0] = 0;
            atoms[1].bond_type[0] = second_type;
            atoms[1].chem_bonds_valence = second_type as i8;
            heap.allocate_model_storage(atoms).unwrap()
        }

        let mut heap = SourceHeap::default();
        let normal = pair(&mut heap, 3, 3);
        assert_eq!(OrigAtData_DecreaseBondOrder(&mut heap, 0, 1, normal), Ok(2));
        let values = heap.slice(normal.as_const()).unwrap();
        assert_eq!(
            (values[0].bond_type[0], values[0].chem_bonds_valence),
            (2, 2)
        );
        assert_eq!(
            (values[1].bond_type[0], values[1].chem_bonds_valence),
            (2, 2)
        );

        let maxed = pair(&mut heap, 3, 3);
        heap.slice_mut(maxed).unwrap()[0].chem_bonds_valence = MAXVAL as i8;
        assert_eq!(OrigAtData_DecreaseBondOrder(&mut heap, 0, 1, maxed), Ok(0));
        assert_eq!(heap.slice(maxed.as_const()).unwrap()[0].bond_type[0], 3);

        let first_low = pair(&mut heap, 1, 3);
        assert_eq!(
            OrigAtData_DecreaseBondOrder(&mut heap, 0, 1, first_low),
            Ok(0)
        );
        assert_eq!(heap.slice(first_low.as_const()).unwrap()[1].bond_type[0], 3);

        let second_low = pair(&mut heap, 2, 1);
        assert_eq!(
            OrigAtData_DecreaseBondOrder(&mut heap, 0, 1, second_low),
            Ok(0)
        );
        let values = heap.slice(second_low.as_const()).unwrap();
        assert_eq!(
            (values[0].bond_type[0], values[0].chem_bonds_valence),
            (1, 1)
        );
        assert_eq!(
            (values[1].bond_type[0], values[1].chem_bonds_valence),
            (1, 1)
        );

        let first_missing = pair(&mut heap, 3, 3);
        heap.slice_mut(first_missing).unwrap()[0].neighbor[0] = 0;
        assert_eq!(
            OrigAtData_DecreaseBondOrder(&mut heap, 0, 1, first_missing),
            Ok(1)
        );
        let values = heap.slice(first_missing.as_const()).unwrap();
        assert_eq!(values[0].bond_type[0], 3);
        assert_eq!(values[1].bond_type[0], 2);

        let second_missing = pair(&mut heap, 3, 3);
        heap.slice_mut(second_missing).unwrap()[1].neighbor[0] = 1;
        assert_eq!(
            OrigAtData_DecreaseBondOrder(&mut heap, 0, 1, second_missing),
            Ok(1)
        );
        let values = heap.slice(second_missing.as_const()).unwrap();
        assert_eq!(values[0].bond_type[0], 2);
        assert_eq!(values[1].bond_type[0], 3);

        let self_bond = heap
            .allocate_model_storage(vec![inp_ATOM {
                valence: 1,
                chem_bonds_valence: 3,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[0] = 0;
                    neighbors
                },
                bond_type: {
                    let mut types = [0; 20];
                    types[0] = 3;
                    types
                },
                ..inp_ATOM::default()
            }])
            .unwrap();
        assert_eq!(
            OrigAtData_DecreaseBondOrder(&mut heap, 0, 0, self_bond),
            Ok(2)
        );
        let self_value = &heap.slice(self_bond.as_const()).unwrap()[0];
        assert_eq!(
            (self_value.bond_type[0], self_value.chem_bonds_valence),
            (1, 1)
        );

        let self_partial = heap
            .allocate_model_storage(vec![inp_ATOM {
                valence: 1,
                chem_bonds_valence: 2,
                neighbor: [0; 20],
                bond_type: {
                    let mut types = [0; 20];
                    types[0] = 2;
                    types
                },
                ..inp_ATOM::default()
            }])
            .unwrap();
        assert_eq!(
            OrigAtData_DecreaseBondOrder(&mut heap, 0, 0, self_partial),
            Ok(0)
        );
        let self_value = &heap.slice(self_partial.as_const()).unwrap()[0];
        assert_eq!(
            (self_value.bond_type[0], self_value.chem_bonds_valence),
            (1, 1)
        );

        assert_eq!(
            OrigAtData_DecreaseBondOrder(&mut heap, -1, 0, normal),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_reopencyclized__line_3953() {
        fn fixture(
            heap: &mut SourceHeap,
            cyclizable: i32,
            end_bond_type: Option<i32>,
            blist: SourceMutPointer<i32>,
        ) -> (
            SourceMutPointer<OAD_PolymerUnit>,
            SourceMutPointer<inp_ATOM>,
            i32,
        ) {
            let atoms = heap
                .allocate_model_storage(vec![inp_ATOM::default(); 4])
                .unwrap();
            let mut bond_count = 0_i32;
            if let Some(bond_type) = end_bond_type {
                OrigAtData_AddBond(heap, 1, 2, atoms, bond_type, 0, &mut bond_count).unwrap();
            }
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    cyclizable,
                    cap1: 1,
                    end_atom1: 2,
                    end_atom2: 3,
                    cap2: 4,
                    nb: 91,
                    nbkbonds: 92,
                    blist,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            (unit, atoms, bond_count)
        }

        let mut heap = SourceHeap::default();
        let (ring_unit, ring_atoms, mut ring_bonds) = fixture(
            &mut heap,
            CLOSING_SRU_RING as i32,
            Some(1),
            SourceMutPointer::null(),
        );
        OAD_PolymerUnit_ReopenCyclized(
            &mut heap,
            ring_unit,
            ring_atoms,
            SourceMutPointer::null(),
            i32::MIN,
            &mut ring_bonds,
        )
        .unwrap();
        let ring = heap.slice(ring_unit.as_const()).unwrap()[0].clone();
        assert_eq!(ring_bonds, 2);
        assert_eq!((ring.nb, ring.nbkbonds), (2, 0));
        assert_eq!(heap.slice(ring.blist.as_const()).unwrap(), &[1, 2, 4, 3]);
        let values = heap.slice(ring_atoms.as_const()).unwrap();
        assert_eq!((values[0].neighbor[0], values[0].bond_type[0]), (1, 1));
        assert_eq!((values[3].neighbor[0], values[3].bond_type[0]), (2, 1));
        assert_eq!(values[1].valence, 1);
        assert_eq!(values[2].valence, 1);

        let (higher_unit, higher_atoms, mut higher_bonds) = fixture(
            &mut heap,
            CLOSING_SRU_HIGHER_ORDER_BOND as i32,
            Some(2),
            SourceMutPointer::null(),
        );
        OAD_PolymerUnit_ReopenCyclized(
            &mut heap,
            higher_unit,
            higher_atoms,
            SourceMutPointer::null(),
            0,
            &mut higher_bonds,
        )
        .unwrap();
        assert_eq!(higher_bonds, 3);
        let values = heap.slice(higher_atoms.as_const()).unwrap();
        assert_eq!(values[1].bond_type[0], 1);
        assert_eq!(values[2].bond_type[0], 1);

        let existing = heap.allocate(vec![9_i32, 8, 7, 6, 5]).unwrap();
        let (radical_unit, radical_atoms, mut radical_bonds) =
            fixture(&mut heap, CLOSING_SRU_DIRADICAL as i32, None, existing);
        heap.slice_mut(radical_atoms).unwrap()[1].radical = RADICAL_TRIPLET as i8;
        OAD_PolymerUnit_ReopenCyclized(
            &mut heap,
            radical_unit,
            radical_atoms,
            SourceMutPointer::null(),
            i32::MAX,
            &mut radical_bonds,
        )
        .unwrap();
        assert_eq!(heap.slice(radical_atoms.as_const()).unwrap()[1].radical, 0);
        assert_eq!(heap.slice(existing.as_const()).unwrap(), &[1, 2, 4, 3, 5]);
        assert_eq!(
            heap.slice(radical_unit.as_const()).unwrap()[0].blist,
            existing
        );

        let (unknown_unit, unknown_atoms, mut unknown_bonds) =
            fixture(&mut heap, i32::MAX, None, SourceMutPointer::null());
        OAD_PolymerUnit_ReopenCyclized(
            &mut heap,
            unknown_unit,
            unknown_atoms,
            SourceMutPointer::null(),
            0,
            &mut unknown_bonds,
        )
        .unwrap();
        assert_eq!(unknown_bonds, 2);

        let mut failure_heap = SourceHeap::default();
        let (failure_unit, failure_atoms, mut failure_bonds) = fixture(
            &mut failure_heap,
            CLOSING_SRU_RING as i32,
            Some(1),
            SourceMutPointer::null(),
        );
        let baseline = failure_heap.live_allocation_count();
        failure_heap.fail_after_allocations(0);
        OAD_PolymerUnit_ReopenCyclized(
            &mut failure_heap,
            failure_unit,
            failure_atoms,
            SourceMutPointer::null(),
            0,
            &mut failure_bonds,
        )
        .unwrap();
        let failed = failure_heap.slice(failure_unit.as_const()).unwrap()[0].clone();
        assert_eq!(failure_bonds, 2);
        assert_eq!((failed.nb, failed.nbkbonds), (2, 0));
        assert!(failed.blist.is_null());
        assert_eq!(failure_heap.live_allocation_count(), baseline);
        let values = failure_heap.slice(failure_atoms.as_const()).unwrap();
        assert_eq!(values[0].neighbor[0], 1);
        assert_eq!(values[3].neighbor[0], 2);
    }

    #[test]
    fn source_port__runichi3__oad_polymer_smartreopencyclizedunits__line_3878() {
        fn fixture(
            heap: &mut SourceHeap,
            frame_scheme: i32,
        ) -> (
            SourceMutPointer<OAD_Polymer>,
            SourceMutPointer<OAD_PolymerUnit>,
            SourceMutPointer<inp_ATOM>,
        ) {
            let mut atoms = vec![inp_ATOM::default(); 3];
            for (index, atom) in atoms.iter_mut().enumerate() {
                atom.orig_at_number = u16::try_from(index + 1).unwrap();
                atom.el_number = 6;
            }
            atoms[1].radical = RADICAL_TRIPLET as i8;
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let row = heap.allocate(vec![2_i32, 2]).unwrap();
            let rows = heap.allocate(vec![row]).unwrap();
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    cap1: 1,
                    cap2: 3,
                    nbkbonds: 1,
                    bkbonds: rows,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            let units = heap.allocate_model_storage(vec![unit]).unwrap();
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer {
                    units,
                    n: 1,
                    really_do_frame_shift: 1,
                    frame_shift_scheme: frame_scheme,
                    ..OAD_Polymer::default()
                }])
                .unwrap();
            (polymer, unit, atoms)
        }

        let mut heap = SourceHeap::default();
        let (polymer, unit, atoms) = fixture(&mut heap, 0);
        let baseline = heap.live_allocation_count();
        let mut bond_count = 0_i32;
        OAD_Polymer_SmartReopenCyclizedUnits(&mut heap, polymer, atoms, 3, &mut bond_count)
            .unwrap();
        assert_eq!(
            heap.slice(polymer.as_const()).unwrap()[0].really_do_frame_shift,
            0
        );
        let reopened = heap.slice(unit.as_const()).unwrap()[0].clone();
        assert_eq!(reopened.cyclizable, CLOSING_SRU_DIRADICAL as i32);
        assert_eq!((reopened.nb, reopened.nbkbonds), (2, 0));
        assert_eq!(
            heap.slice(reopened.blist.as_const()).unwrap(),
            &[1, 2, 3, 2]
        );
        assert_eq!(heap.slice(atoms.as_const()).unwrap()[1].radical, 0);
        assert_eq!(bond_count, 2);
        assert_eq!(heap.live_allocation_count(), baseline + 1);

        let null_baseline = heap.live_allocation_count();
        OAD_Polymer_SmartReopenCyclizedUnits(
            &mut heap,
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            3,
            &mut bond_count,
        )
        .unwrap();
        assert_eq!(heap.live_allocation_count(), null_baseline);

        let empty = heap
            .allocate_model_storage(vec![OAD_Polymer {
                really_do_frame_shift: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let empty_baseline = heap.live_allocation_count();
        OAD_Polymer_SmartReopenCyclizedUnits(
            &mut heap,
            empty,
            SourceMutPointer::null(),
            3,
            &mut bond_count,
        )
        .unwrap();
        assert_eq!(
            heap.slice(empty.as_const()).unwrap()[0].really_do_frame_shift,
            1
        );
        assert_eq!(heap.live_allocation_count(), empty_baseline);

        let (disabled, disabled_unit, disabled_atoms) = fixture(&mut heap, 0);
        heap.slice_mut(disabled).unwrap()[0].really_do_frame_shift = 0;
        let disabled_baseline = heap.live_allocation_count();
        OAD_Polymer_SmartReopenCyclizedUnits(
            &mut heap,
            disabled,
            disabled_atoms,
            3,
            &mut bond_count,
        )
        .unwrap();
        assert_eq!(heap.slice(disabled_unit.as_const()).unwrap()[0].nbkbonds, 1);
        assert_eq!(heap.live_allocation_count(), disabled_baseline);

        let (zero_nat, zero_nat_unit, zero_nat_atoms) = fixture(&mut heap, 0);
        let zero_baseline = heap.live_allocation_count();
        OAD_Polymer_SmartReopenCyclizedUnits(
            &mut heap,
            zero_nat,
            zero_nat_atoms,
            0,
            &mut bond_count,
        )
        .unwrap();
        assert_eq!(heap.slice(zero_nat_unit.as_const()).unwrap()[0].nbkbonds, 1);
        assert_eq!(
            heap.slice(zero_nat.as_const()).unwrap()[0].really_do_frame_shift,
            1
        );
        assert_eq!(heap.live_allocation_count(), zero_baseline);

        let (null_atoms, null_atoms_unit, _) = fixture(&mut heap, 0);
        let null_atoms_baseline = heap.live_allocation_count();
        OAD_Polymer_SmartReopenCyclizedUnits(
            &mut heap,
            null_atoms,
            SourceMutPointer::null(),
            3,
            &mut bond_count,
        )
        .unwrap();
        assert_eq!(
            heap.slice(null_atoms_unit.as_const()).unwrap()[0].nbkbonds,
            1
        );
        assert_eq!(
            heap.slice(null_atoms.as_const()).unwrap()[0].really_do_frame_shift,
            1
        );
        assert_eq!(heap.live_allocation_count(), null_atoms_baseline);

        let mut failure_heap = SourceHeap::default();
        let (failure_polymer, failure_unit, failure_atoms) = fixture(&mut failure_heap, 0);
        let failure_baseline = failure_heap.live_allocation_count();
        failure_heap.fail_after_allocations(0);
        let mut failure_bonds = 0_i32;
        OAD_Polymer_SmartReopenCyclizedUnits(
            &mut failure_heap,
            failure_polymer,
            failure_atoms,
            3,
            &mut failure_bonds,
        )
        .unwrap();
        assert_eq!(
            failure_heap.slice(failure_polymer.as_const()).unwrap()[0].really_do_frame_shift,
            1
        );
        assert_eq!(
            failure_heap.slice(failure_unit.as_const()).unwrap()[0].nbkbonds,
            1
        );
        assert_eq!(failure_heap.live_allocation_count(), failure_baseline);

        let (none_scheme, none_unit, none_atoms) =
            fixture(&mut heap, tagFrameShifScheme_FSS_NONE as i32);
        let none_baseline = heap.live_allocation_count();
        let mut none_bonds = 0_i32;
        OAD_Polymer_SmartReopenCyclizedUnits(
            &mut heap,
            none_scheme,
            none_atoms,
            3,
            &mut none_bonds,
        )
        .unwrap();
        assert_eq!(
            heap.slice(none_scheme.as_const()).unwrap()[0].really_do_frame_shift,
            0
        );
        assert_eq!(heap.slice(none_unit.as_const()).unwrap()[0].nbkbonds, 1);
        assert_eq!(none_bonds, 0);
        assert_eq!(heap.live_allocation_count(), none_baseline);

        let (invalid_cap, invalid_unit, invalid_atoms) = fixture(&mut heap, 0);
        heap.slice_mut(invalid_unit).unwrap()[0].cap2 = 4;
        let invalid_baseline = heap.live_allocation_count();
        OAD_Polymer_SmartReopenCyclizedUnits(
            &mut heap,
            invalid_cap,
            invalid_atoms,
            3,
            &mut none_bonds,
        )
        .unwrap();
        assert_eq!(heap.slice(invalid_unit.as_const()).unwrap()[0].nbkbonds, 1);
        assert_eq!(heap.live_allocation_count(), invalid_baseline);
    }

    #[test]
    fn source_port__runichi3__oad_validateandsortoutpseudoelementatoms__line_4396() {
        fn atom(name: &[u8], tail: i8) -> inp_ATOM {
            let mut atom = inp_ATOM {
                elname: [tail; 6],
                ..inp_ATOM::default()
            };
            for (target, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
                *target = source as i8;
            }
            atom
        }

        fn original(heap: &mut SourceHeap, atoms: Vec<inp_ATOM>) -> ORIG_ATOM_DATA {
            ORIG_ATOM_DATA {
                num_inp_atoms: atoms.len() as i32,
                at: heap.allocate_model_storage(atoms).unwrap(),
                n_zy: 91,
                ..ORIG_ATOM_DATA::default()
            }
        }

        fn polymer(
            heap: &mut SourceHeap,
            units: Vec<Option<OAD_PolymerUnit>>,
            initial_pzz: i32,
        ) -> SourceMutPointer<OAD_Polymer> {
            let mut unit_pointers = Vec::with_capacity(units.len());
            for unit in units {
                unit_pointers.push(
                    unit.map(|value| heap.allocate_model_storage(vec![value]).unwrap())
                        .unwrap_or_default(),
                );
            }
            let number_of_units = unit_pointers.len() as i32;
            let unit_pointers = if unit_pointers.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(unit_pointers).unwrap()
            };
            heap.allocate_model_storage(vec![OAD_Polymer {
                n: number_of_units,
                units: unit_pointers,
                n_pzz: initial_pzz,
                ..OAD_Polymer::default()
            }])
            .unwrap()
        }

        fn error_bytes(buffer: &[i8; 256]) -> Vec<u8> {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            buffer[..end].iter().map(|byte| *byte as u8).collect()
        }

        let mut heap = SourceHeap::default();

        let mut ordinary = original(&mut heap, vec![atom(b"C\0", 11), atom(b"N\0", 12)]);
        let ordinary_before = heap.slice(ordinary.at.as_const()).unwrap().to_vec();
        let mut error = 0;
        let mut message = [0_i8; 256];
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut ordinary,
            POLYMERS_NO as i32,
            0,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 0);
        assert_eq!(ordinary.n_zy, 0);
        assert_eq!(heap.slice(ordinary.at.as_const()).unwrap(), ordinary_before);
        assert!(error_bytes(&message).is_empty());

        let mut allowed = original(
            &mut heap,
            vec![
                atom(b"*\0", 21),
                atom(b"Zz\0", 22),
                atom(b"Zy\0", 23),
                atom(b"O\0", 24),
            ],
        );
        heap.slice_mut(allowed.at).unwrap()[2].iso_atw_diff = i8::MAX;
        heap.slice_mut(allowed.at).unwrap()[2].valence = -1;
        heap.slice_mut(allowed.at).unwrap()[2].chem_bonds_valence = -1;
        message.fill(0);
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut allowed,
            POLYMERS_NO as i32,
            1,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 0);
        assert_eq!(allowed.n_zy, 3);
        let values = heap.slice(allowed.at.as_const()).unwrap();
        for (index, tail) in [(0, 21), (1, 22), (2, 23)] {
            assert_eq!(&values[index].elname[..3], &[b'Z' as i8, b'y' as i8, 0]);
            assert_eq!(&values[index].elname[3..], &[tail; 3]);
        }
        assert_eq!(&values[3].elname[..2], &[b'O' as i8, 0]);

        let mut bonding = original(&mut heap, vec![atom(b"*\0", 31), atom(b"Zz\0", 32)]);
        heap.slice_mut(bonding.at).unwrap()[0].valence = 2;
        heap.slice_mut(bonding.at).unwrap()[1].chem_bonds_valence = 2;
        error = 0;
        message.fill(0);
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut bonding,
            POLYMERS_NO as i32,
            1,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 77);
        assert_eq!(bonding.n_zy, 2);
        let values = heap.slice(bonding.at.as_const()).unwrap();
        assert_eq!(&values[0].elname[..2], &[b'*' as i8, 0]);
        assert_eq!(&values[1].elname[..3], &[b'Z' as i8, b'z' as i8, 0]);
        assert_eq!(error_bytes(&message), b"Invalid pseudo element(s) bonding");

        let mut prohibited = original(
            &mut heap,
            vec![atom(b"*\0", 41), atom(b"Zz\0", 42), atom(b"Zy\0", 43)],
        );
        error = 0;
        message.fill(0);
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut prohibited,
            POLYMERS_NO as i32,
            0,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 75);
        assert_eq!(prohibited.n_zy, 3);
        assert_eq!(
            error_bytes(&message),
            b"Invalid element(s): *; Zz; Zy; Polymer-unrelated pseudoatoms are not allowed"
        );
        let values = heap.slice(prohibited.at.as_const()).unwrap();
        assert_eq!(&values[0].elname[..2], &[b'*' as i8, 0]);
        assert_eq!(&values[1].elname[..3], &[b'Z' as i8, b'z' as i8, 0]);
        assert_eq!(&values[2].elname[..3], &[b'Z' as i8, b'y' as i8, 0]);

        let mut treatment_only = original(&mut heap, vec![atom(b"*\0", 51)]);
        treatment_only.valid_polymer = 1;
        treatment_only.polymer = polymer(&mut heap, Vec::new(), 83);
        error = 0;
        message.fill(0);
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut treatment_only,
            1,
            0,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 74);
        assert_eq!(treatment_only.valid_polymer, 0);
        assert_eq!(treatment_only.n_zy, 1);
        assert_eq!(
            heap.slice(treatment_only.polymer.as_const()).unwrap()[0].n_pzz,
            0
        );
        assert_eq!(
            error_bytes(&message),
            b"Polymer-unrelated pseudoatoms are not allowed"
        );

        let mut non_boolean_disallowed = original(&mut heap, vec![atom(b"Zz\0", 61)]);
        error = 0;
        message.fill(0);
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut non_boolean_disallowed,
            POLYMERS_NO as i32,
            2,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 75);
        assert_eq!(error_bytes(&message), b"Invalid element(s): Zz");

        let mut non_boolean_allowed = original(&mut heap, vec![atom(b"Zz\0", 62)]);
        error = 0;
        message.fill(0);
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut non_boolean_allowed,
            1,
            2,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 0);
        assert_eq!(non_boolean_allowed.n_zy, 1);
        assert_eq!(
            &heap.slice(non_boolean_allowed.at.as_const()).unwrap()[0].elname[..3],
            &[b'Z' as i8, b'y' as i8, 0]
        );

        let units = vec![
            Some(OAD_PolymerUnit {
                cap1: 1,
                cap2: 3,
                cap1_is_undef: 1,
                cap2_is_undef: -1,
                ..OAD_PolymerUnit::default()
            }),
            None,
            Some(OAD_PolymerUnit {
                cap1: 3,
                cap2: 5,
                cap1_is_undef: 1,
                cap2_is_undef: 0,
                ..OAD_PolymerUnit::default()
            }),
        ];
        let polymer_pointer = polymer(&mut heap, units, 97);
        let mut polymer_data = original(
            &mut heap,
            vec![
                atom(b"Zz\0", 71),
                atom(b"C\0", 72),
                atom(b"*\0", 73),
                atom(b"N\0", 74),
                atom(b"Zy\0", 75),
            ],
        );
        polymer_data.polymer = polymer_pointer;
        polymer_data.valid_polymer = 1;
        error = 0;
        message.fill(0);
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut polymer_data,
            1,
            1,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 0);
        assert_eq!(polymer_data.valid_polymer, 1);
        assert_eq!(polymer_data.n_zy, 1);
        assert_eq!(heap.slice(polymer_pointer.as_const()).unwrap()[0].n_pzz, 2);
        let values = heap.slice(polymer_data.at.as_const()).unwrap();
        assert_eq!(&values[0].elname[..3], &[b'Z' as i8, b'z' as i8, 0]);
        assert_eq!(&values[0].elname[3..], &[71; 3]);
        assert_eq!(&values[2].elname[..3], &[b'Z' as i8, b'z' as i8, 0]);
        assert_eq!(&values[2].elname[3..], &[73; 3]);
        assert_eq!(&values[4].elname[..3], &[b'Z' as i8, b'y' as i8, 0]);

        let repeated_polymer = polymer(
            &mut heap,
            vec![Some(OAD_PolymerUnit {
                cap1: 1,
                cap2: 1,
                cap1_is_undef: 1,
                cap2_is_undef: 1,
                ..OAD_PolymerUnit::default()
            })],
            98,
        );
        let mut repeated = original(&mut heap, vec![atom(b"Zz\0", 81)]);
        repeated.polymer = repeated_polymer;
        repeated.valid_polymer = 1;
        error = 0;
        message.fill(0);
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut repeated,
            1,
            0,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 74);
        assert_eq!(repeated.n_zy, -1);
        assert_eq!(repeated.valid_polymer, 0);
        assert_eq!(heap.slice(repeated_polymer.as_const()).unwrap()[0].n_pzz, 2);

        let empty_polymer = polymer(&mut heap, Vec::new(), 99);
        let mut preexisting = original(&mut heap, vec![atom(b"C\0", 91)]);
        preexisting.polymer = empty_polymer;
        preexisting.valid_polymer = 1;
        error = 1234;
        message.fill(0);
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut preexisting,
            POLYMERS_NO as i32,
            1,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 1234);
        assert_eq!(preexisting.valid_polymer, 0);
        assert_eq!(preexisting.n_zy, 0);
        assert_eq!(heap.slice(empty_polymer.as_const()).unwrap()[0].n_pzz, 0);
        assert!(error_bytes(&message).is_empty());

        let mut null_error_text = original(&mut heap, vec![atom(b"*\0", 92)]);
        error = 0;
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut null_error_text,
            POLYMERS_NO as i32,
            0,
            &mut error,
            None,
        )
        .unwrap();
        assert_eq!(error, 75);

        let mut unterminated = original(
            &mut heap,
            vec![inp_ATOM {
                elname: [b'Z' as i8; 6],
                ..inp_ATOM::default()
            }],
        );
        unterminated.n_zy = 101;
        error = 0;
        assert_eq!(
            OAD_ValidateAndSortOutPseudoElementAtoms(
                &mut heap,
                &mut unterminated,
                POLYMERS_NO as i32,
                1,
                &mut error,
                Some(&mut message),
            ),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(unterminated.n_zy, 101);

        let partial_atoms = heap
            .allocate_model_storage(vec![atom(b"*\0", 102)])
            .unwrap();
        let mut short_atoms = ORIG_ATOM_DATA {
            at: partial_atoms,
            num_inp_atoms: 2,
            n_zy: 103,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            OAD_ValidateAndSortOutPseudoElementAtoms(
                &mut heap,
                &mut short_atoms,
                POLYMERS_NO as i32,
                1,
                &mut error,
                Some(&mut message),
            ),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(short_atoms.n_zy, 103);
        assert_eq!(
            &heap.slice(partial_atoms.as_const()).unwrap()[0].elname[..3],
            &[b'Z' as i8, b'y' as i8, 0]
        );

        let mut null_polymer = original(&mut heap, vec![atom(b"*\0", 104)]);
        null_polymer.valid_polymer = 1;
        null_polymer.n_zy = 105;
        error = 0;
        assert_eq!(
            OAD_ValidateAndSortOutPseudoElementAtoms(
                &mut heap,
                &mut null_polymer,
                1,
                1,
                &mut error,
                Some(&mut message),
            ),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(null_polymer.n_zy, 0);
        assert_eq!(null_polymer.valid_polymer, 1);
        assert_eq!(
            &heap.slice(null_polymer.at.as_const()).unwrap()[0].elname[..3],
            &[b'Z' as i8, b'y' as i8, 0]
        );

        let invalid_cap_polymer = polymer(
            &mut heap,
            vec![Some(OAD_PolymerUnit {
                cap1: 1,
                cap2: 9,
                cap1_is_undef: 1,
                cap2_is_undef: 1,
                ..OAD_PolymerUnit::default()
            })],
            106,
        );
        let mut invalid_cap = original(&mut heap, vec![atom(b"Zz\0", 107), atom(b"Zz\0", 108)]);
        invalid_cap.polymer = invalid_cap_polymer;
        invalid_cap.valid_polymer = 1;
        invalid_cap.n_zy = 109;
        error = 0;
        assert_eq!(
            OAD_ValidateAndSortOutPseudoElementAtoms(
                &mut heap,
                &mut invalid_cap,
                1,
                1,
                &mut error,
                Some(&mut message),
            ),
            Err(SourceHeapError::UnsupportedSourceBehavior)
        );
        assert_eq!(invalid_cap.n_zy, 0);
        assert_eq!(invalid_cap.valid_polymer, 1);
        assert_eq!(
            heap.slice(invalid_cap_polymer.as_const()).unwrap()[0].n_pzz,
            106
        );
        let values = heap.slice(invalid_cap.at.as_const()).unwrap();
        assert_eq!(&values[0].elname[..3], &[b'Z' as i8, b'z' as i8, 0]);
        assert_eq!(&values[1].elname[..3], &[b'Z' as i8, b'y' as i8, 0]);

        let mut negative_count = ORIG_ATOM_DATA {
            num_inp_atoms: -1,
            n_zy: 110,
            ..ORIG_ATOM_DATA::default()
        };
        error = 0;
        OAD_ValidateAndSortOutPseudoElementAtoms(
            &mut heap,
            &mut negative_count,
            POLYMERS_NO as i32,
            1,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(negative_count.n_zy, 0);
        assert_eq!(error, 0);
    }

    #[test]
    fn source_port__runichi3__oad_validatepolymerandpseudoelementdata__line_1516() {
        #[derive(Clone)]
        struct UnitSpec {
            type_: i32,
            subtype: i32,
            conn: i32,
            atoms: Vec<i32>,
            bonds: Option<Vec<i32>>,
            nb: i32,
        }

        fn atom(name: &[u8]) -> inp_ATOM {
            let mut atom = inp_ATOM::default();
            for (target, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
                *target = source as i8;
            }
            atom
        }

        fn spec(
            type_: u32,
            subtype: u32,
            conn: u32,
            atoms: &[i32],
            bonds: Option<&[i32]>,
            nb: i32,
        ) -> UnitSpec {
            UnitSpec {
                type_: type_ as i32,
                subtype: subtype as i32,
                conn: conn as i32,
                atoms: atoms.to_vec(),
                bonds: bonds.map(<[i32]>::to_vec),
                nb,
            }
        }

        fn fixture(
            heap: &mut SourceHeap,
            names: &[&[u8]],
            specs: &[UnitSpec],
            polymer_treat: i32,
            number_of_bonds: i32,
        ) -> (
            ORIG_ATOM_DATA,
            SourceMutPointer<OAD_Polymer>,
            Vec<SourceMutPointer<OAD_PolymerUnit>>,
        ) {
            let atoms = heap
                .allocate_model_storage(names.iter().map(|name| atom(name)).collect())
                .unwrap();
            let mut unit_pointers = Vec::new();
            for value in specs {
                let atom_list = heap.allocate_model_storage(value.atoms.clone()).unwrap();
                let bond_list = value
                    .bonds
                    .as_ref()
                    .map(|bonds| heap.allocate_model_storage(bonds.clone()).unwrap())
                    .unwrap_or_default();
                unit_pointers.push(
                    heap.allocate_model_storage(vec![OAD_PolymerUnit {
                        type_: value.type_,
                        subtype: value.subtype,
                        conn: value.conn,
                        na: i32::try_from(value.atoms.len()).unwrap(),
                        nb: value.nb,
                        alist: atom_list,
                        blist: bond_list,
                        cyclizable: 81,
                        cyclized: 82,
                        nbkbonds: 83,
                        ..OAD_PolymerUnit::default()
                    }])
                    .unwrap(),
                );
            }
            let units = heap.allocate_model_storage(unit_pointers.clone()).unwrap();
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer {
                    units,
                    n: i32::try_from(unit_pointers.len()).unwrap(),
                    treat: polymer_treat,
                    ..OAD_Polymer::default()
                }])
                .unwrap();
            (
                ORIG_ATOM_DATA {
                    at: atoms,
                    num_inp_atoms: i32::try_from(names.len()).unwrap(),
                    num_inp_bonds: number_of_bonds,
                    polymer,
                    valid_polymer: 91,
                    n_zy: 92,
                    ..ORIG_ATOM_DATA::default()
                },
                polymer,
                unit_pointers,
            )
        }

        fn invoke(
            heap: &mut SourceHeap,
            original: &mut ORIG_ATOM_DATA,
            treat_polymers: i32,
            allow_non_polymer_zz: i32,
            no_warnings: i32,
        ) -> (Result<i32, SourceHeapError>, [i8; 256]) {
            let mut message = [0_i8; 256];
            let result = OAD_ValidatePolymerAndPseudoElementData(
                heap,
                original,
                treat_polymers,
                allow_non_polymer_zz,
                Some(&mut message),
                no_warnings,
            );
            (result, message)
        }

        fn message(buffer: &[i8; 256]) -> Vec<u8> {
            let end = buffer.iter().position(|byte| *byte == 0).unwrap();
            buffer[..end].iter().map(|byte| *byte as u8).collect()
        }

        fn assert_error(
            names: &[&[u8]],
            specs: &[UnitSpec],
            expected_error: i32,
            expected_message: &[u8],
        ) {
            let mut heap = SourceHeap::default();
            let (mut original, _, _) = fixture(&mut heap, names, specs, 1, 0);
            let (result, text) = invoke(&mut heap, &mut original, 1, 1, 0);
            assert_eq!(result, Ok(expected_error));
            assert_eq!(message(&text), expected_message);
            assert_eq!(original.valid_polymer, 0);
        }

        let mut no_polymer_heap = SourceHeap::default();
        let mut no_polymer = ORIG_ATOM_DATA::default();
        let (result, text) = invoke(&mut no_polymer_heap, &mut no_polymer, 0, 1, 0);
        assert_eq!(result, Ok(0));
        assert_eq!(no_polymer.valid_polymer, 0);
        assert!(message(&text).is_empty());

        let mut enabled_without_polymer = ORIG_ATOM_DATA::default();
        let (result, _) = invoke(&mut no_polymer_heap, &mut enabled_without_polymer, 1, 1, 0);
        assert_eq!(result, Ok(0));
        assert_eq!(enabled_without_polymer.valid_polymer, 0);

        let mut empty_polymer_heap = SourceHeap::default();
        let (mut empty_polymer, _, _) = fixture(&mut empty_polymer_heap, &[], &[], 1, 0);
        let (result, text) = invoke(&mut empty_polymer_heap, &mut empty_polymer, 1, 1, 0);
        assert_eq!(result, Ok(0));
        assert_eq!(empty_polymer.valid_polymer, 1);
        assert!(message(&text).is_empty());

        for (allow, expected, expected_message, expected_zy) in [
            (
                0,
                75,
                b"Invalid element(s): Zz; Polymer-unrelated pseudoatoms are not allowed".as_slice(),
                1,
            ),
            (1, 0, b"".as_slice(), 1),
        ] {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(vec![atom(b"Zz\0")]).unwrap();
            let mut original = ORIG_ATOM_DATA {
                at: atoms,
                num_inp_atoms: 1,
                ..ORIG_ATOM_DATA::default()
            };
            let (result, text) = invoke(&mut heap, &mut original, 0, allow, 0);
            assert_eq!(result, Ok(expected));
            assert_eq!(message(&text), expected_message);
            assert_eq!(original.n_zy, expected_zy);
            assert_eq!(original.valid_polymer, 0);
        }

        assert_error(
            &[b"C\0"],
            &[spec(
                POLYMER_STY_COP,
                POLYMER_SST_NON,
                POLYMER_CONN_NON,
                &[1],
                None,
                0,
            )],
            9001,
            b"Copolymer must contain more than one unit",
        );
        for subtype in [POLYMER_SST_RAN, POLYMER_SST_ALT, POLYMER_SST_BLK] {
            assert_error(
                &[b"C\0"],
                &[spec(
                    POLYMER_STY_MON,
                    subtype,
                    POLYMER_CONN_NON,
                    &[1],
                    None,
                    0,
                )],
                9002,
                b"Single polymer unit may not be RAN/ALT/BLO",
            );
        }
        assert_error(
            &[b"C\0"],
            &[spec(
                POLYMER_STY_MON,
                POLYMER_SST_NON,
                POLYMER_CONN_NON,
                &[1],
                Some(&[1, 1]),
                1,
            )],
            9003,
            b"Number of crossing bonds in polymer unit is not 0 or 2",
        );
        assert_error(
            &[b"C\0"],
            &[spec(
                POLYMER_STY_MON,
                POLYMER_SST_NON,
                POLYMER_CONN_NON,
                &[],
                None,
                0,
            )],
            9004,
            b"Empty polymer unit",
        );
        assert_error(
            &[b"C\0"],
            &[spec(
                POLYMER_STY_MON,
                POLYMER_SST_NON,
                POLYMER_CONN_NON,
                &[1, 1],
                None,
                0,
            )],
            9005,
            b"Too large polymer unit",
        );
        for invalid_atom in [0, 2] {
            assert_error(
                &[b"C\0"],
                &[spec(
                    POLYMER_STY_MON,
                    POLYMER_SST_NON,
                    POLYMER_CONN_NON,
                    &[invalid_atom],
                    None,
                    0,
                )],
                9006,
                b"Invalid atom number in polymer unit",
            );
        }

        let source_spec = spec(
            POLYMER_STY_SRU,
            POLYMER_SST_NON,
            POLYMER_CONN_HT,
            &[1],
            None,
            0,
        );
        let mut source_heap = SourceHeap::default();
        let (mut source, _, source_units) =
            fixture(&mut source_heap, &[b"C\0"], &[source_spec.clone()], 1, 0);
        let (result, text) = invoke(&mut source_heap, &mut source, 1, 1, 0);
        assert_eq!(result, Ok(0));
        assert_eq!(source.valid_polymer, 1);
        let source_unit = &source_heap.slice(source_units[0].as_const()).unwrap()[0];
        assert_eq!(source_unit.type_, POLYMER_STY_MON as i32);
        assert_eq!(source_unit.conn, POLYMER_CONN_NON as i32);
        assert_eq!(source_unit.nbkbonds, 0);
        assert_eq!(source_unit.cyclizable, CLOSING_SRU_NOT_APPLICABLE as i32);
        assert_eq!(source_unit.cyclized, 0);
        assert_eq!(
            message(&text),
            b"Converted src-based polymer unit type to MON; Ignore connection pattern for src-based polymer unit"
        );

        let mut quiet_heap = SourceHeap::default();
        let (mut quiet, _, quiet_units) = fixture(&mut quiet_heap, &[b"C\0"], &[source_spec], 1, 0);
        let (result, text) = invoke(&mut quiet_heap, &mut quiet, 1, 1, 1);
        assert_eq!(result, Ok(0));
        assert!(message(&text).is_empty());
        assert_eq!(
            quiet_heap.slice(quiet_units[0].as_const()).unwrap()[0].type_,
            POLYMER_STY_MON as i32
        );

        let embedding_specs = [
            spec(
                POLYMER_STY_SRU,
                POLYMER_SST_NON,
                POLYMER_CONN_NON,
                &[1, 2],
                None,
                0,
            ),
            spec(
                POLYMER_STY_MON,
                POLYMER_SST_NON,
                POLYMER_CONN_NON,
                &[1],
                None,
                0,
            ),
            spec(
                POLYMER_STY_MON,
                POLYMER_SST_NON,
                POLYMER_CONN_NON,
                &[2],
                None,
                0,
            ),
        ];
        let mut embedding_heap = SourceHeap::default();
        let (mut embedding, _, embedding_units) = fixture(
            &mut embedding_heap,
            &[b"C\0", b"N\0"],
            &embedding_specs,
            1,
            0,
        );
        let (result, text) = invoke(&mut embedding_heap, &mut embedding, 1, 1, 0);
        assert_eq!(result, Ok(0));
        let outer = &embedding_heap.slice(embedding_units[0].as_const()).unwrap()[0];
        assert_eq!(outer.type_, POLYMER_STY_COP as i32);
        assert_eq!(outer.subtype, POLYMER_SST_RAN as i32);
        assert_eq!(
            message(&text),
            b"Convert multiple-subunits unit to copolymer; Set missing copolymer subtype to RAN"
        );

        assert_error(
            &[b"C\0", b"N\0"],
            &[
                spec(
                    POLYMER_STY_COP,
                    POLYMER_SST_NON,
                    POLYMER_CONN_NON,
                    &[1, 2],
                    None,
                    0,
                ),
                spec(
                    POLYMER_STY_MON,
                    POLYMER_SST_NON,
                    POLYMER_CONN_NON,
                    &[1],
                    None,
                    0,
                ),
            ],
            9027,
            b"Polymer COP unit contains a single SRU instead of multiple",
        );
        assert_error(
            &[b"C\0", b"N\0", b"O\0"],
            &[
                spec(
                    POLYMER_STY_COP,
                    POLYMER_SST_NON,
                    POLYMER_CONN_NON,
                    &[2],
                    Some(&[1, 2, 3, 2]),
                    2,
                ),
                spec(
                    POLYMER_STY_MON,
                    POLYMER_SST_NON,
                    POLYMER_CONN_NON,
                    &[1],
                    None,
                    0,
                ),
            ],
            9026,
            b"Polymer COP unit contains bracket-crossing bonds, not supported",
        );

        fn structure_case(
            names: &[&[u8]],
            atoms: &[i32],
            bonds: &[i32],
            type_: u32,
            conn: u32,
            failure_after: Option<u64>,
        ) -> (
            SourceHeap,
            ORIG_ATOM_DATA,
            SourceMutPointer<OAD_PolymerUnit>,
            Result<i32, SourceHeapError>,
            [i8; 256],
        ) {
            let mut heap = SourceHeap::default();
            let (mut original, _, units) = fixture(
                &mut heap,
                names,
                &[spec(type_, POLYMER_SST_NON, conn, atoms, Some(bonds), 2)],
                1,
                2,
            );
            if let Some(successful) = failure_after {
                heap.fail_after_allocations(successful);
            }
            let (result, text) = invoke(&mut heap, &mut original, 1, 1, 0);
            (heap, original, units[0], result, text)
        }

        for (atoms, bonds, expected) in [
            (vec![2], vec![1, 2, 3, 2], CLOSING_SRU_DIRADICAL),
            (vec![2, 3], vec![1, 2, 4, 3], CLOSING_SRU_HIGHER_ORDER_BOND),
            (vec![2, 3, 4], vec![1, 2, 5, 4], CLOSING_SRU_RING),
        ] {
            let names = vec![b"C\0".as_slice(); atoms.len() + 2];
            let (heap, original, unit_pointer, result, text) = structure_case(
                &names,
                &atoms,
                &bonds,
                POLYMER_STY_SRU,
                POLYMER_CONN_HT,
                None,
            );
            assert_eq!(result, Ok(0));
            assert_eq!(original.valid_polymer, 1);
            assert!(message(&text).is_empty());
            let unit = &heap.slice(unit_pointer.as_const()).unwrap()[0];
            assert_eq!(unit.cyclizable, expected as i32);
            assert_eq!(unit.maxbkbonds, 4);
            assert!(!unit.bkbonds.is_null());
            assert_eq!(heap.slice(unit.bkbonds.as_const()).unwrap().len(), 4);
        }

        for type_ in [POLYMER_STY_MOD, POLYMER_STY_CRO, POLYMER_STY_MER] {
            let (heap, original, unit_pointer, result, _) = structure_case(
                &[b"C\0", b"N\0", b"O\0"],
                &[2],
                &[1, 2, 3, 2],
                type_,
                POLYMER_CONN_HT,
                None,
            );
            assert_eq!(result, Ok(0));
            assert_eq!(original.valid_polymer, 1);
            assert_eq!(
                heap.slice(unit_pointer.as_const()).unwrap()[0].cyclizable,
                CLOSING_SRU_DIRADICAL as i32
            );
        }

        let (heap, original, unit_pointer, result, text) = structure_case(
            &[b"C\0", b"N\0", b"O\0"],
            &[2],
            &[1, 2, 3, 2],
            POLYMER_STY_SRU,
            POLYMER_CONN_NON,
            None,
        );
        assert_eq!(result, Ok(0));
        assert_eq!(original.valid_polymer, 1);
        let unit = &heap.slice(unit_pointer.as_const()).unwrap()[0];
        assert_eq!(unit.conn, POLYMER_CONN_EU as i32);
        assert_eq!(unit.cyclizable, CLOSING_SRU_NOT_APPLICABLE as i32);
        assert!(unit.bkbonds.is_null());
        assert_eq!(
            message(&text),
            b"Set missing copolymer unit connection to EU"
        );

        for (names, expected_error) in [
            (vec![b"H\0".as_slice(), b"C\0", b"O\0"], 9030),
            (vec![b"C\0".as_slice(), b"H\0", b"O\0"], 9031),
            (vec![b"D\0".as_slice(), b"C\0", b"O\0"], 9030),
            (vec![b"T\0".as_slice(), b"C\0", b"O\0"], 9030),
        ] {
            let (_, original, _, result, text) = structure_case(
                &names,
                &[2],
                &[1, 2, 3, 2],
                POLYMER_STY_SRU,
                POLYMER_CONN_HT,
                None,
            );
            assert_eq!(result, Ok(expected_error));
            assert_eq!(original.valid_polymer, 0);
            assert_eq!(message(&text), b"H as polymer end group is not supported");
        }

        for successful_allocations in [0_u64, 1] {
            let (heap, original, unit_pointer, result, text) = structure_case(
                &[b"C\0", b"N\0", b"O\0"],
                &[2],
                &[1, 2, 3, 2],
                POLYMER_STY_SRU,
                POLYMER_CONN_HT,
                Some(successful_allocations),
            );
            assert_eq!(result, Ok(1));
            assert_eq!(original.valid_polymer, 0);
            assert_eq!(message(&text), b"Not enough memory (polymers)");
            let unit = &heap.slice(unit_pointer.as_const()).unwrap()[0];
            assert_eq!(unit.maxbkbonds, 4);
            if successful_allocations == 0 {
                assert!(unit.bkbonds.is_null());
            } else {
                assert!(!unit.bkbonds.is_null());
                assert_eq!(heap.slice(unit.bkbonds.as_const()).unwrap().len(), 4);
                assert!(heap.slice(unit.bkbonds.as_const()).unwrap()[0].is_null());
            }
        }

        let mut replacement_heap = SourceHeap::default();
        let (mut replacement, _, replacement_units) = fixture(
            &mut replacement_heap,
            &[b"C\0", b"N\0", b"O\0"],
            &[spec(
                POLYMER_STY_SRU,
                POLYMER_SST_NON,
                POLYMER_CONN_HT,
                &[2],
                Some(&[1, 2, 3, 2]),
                2,
            )],
            1,
            2,
        );
        let old_row = replacement_heap
            .allocate_model_storage(vec![71_i32, 72])
            .unwrap();
        let old_matrix = replacement_heap
            .allocate_model_storage(vec![old_row])
            .unwrap();
        {
            let unit = &mut replacement_heap.slice_mut(replacement_units[0]).unwrap()[0];
            unit.maxbkbonds = 1;
            unit.bkbonds = old_matrix;
        }
        let (result, _) = invoke(&mut replacement_heap, &mut replacement, 1, 1, 0);
        assert_eq!(result, Ok(0));
        assert!(matches!(
            replacement_heap.slice(old_matrix.as_const()),
            Err(SourceHeapError::MissingAllocation)
        ));
        assert!(matches!(
            replacement_heap.slice(old_row.as_const()),
            Err(SourceHeapError::MissingAllocation)
        ));
        let replacement_unit = &replacement_heap
            .slice(replacement_units[0].as_const())
            .unwrap()[0];
        assert!(!replacement_unit.bkbonds.is_null());
        assert_ne!(replacement_unit.bkbonds, old_matrix);

        let (_, original, _, result, text) = structure_case(
            &[b"C\0", b"N\0", b"O\0", b"P\0"],
            &[2],
            &[1, 3, 4, 3],
            POLYMER_STY_SRU,
            POLYMER_CONN_HT,
            None,
        );
        assert_eq!(result, Ok(9032));
        assert_eq!(original.valid_polymer, 0);
        assert_eq!(message(&text), b"Caps of polymer unit lie inside it");

        let mixed_specs = [
            spec(
                POLYMER_STY_SRU,
                POLYMER_SST_NON,
                POLYMER_CONN_NON,
                &[1],
                None,
                0,
            ),
            spec(
                POLYMER_STY_MON,
                POLYMER_SST_NON,
                POLYMER_CONN_NON,
                &[3],
                Some(&[2, 3, 4, 3]),
                2,
            ),
        ];
        let mut mixed_heap = SourceHeap::default();
        let (mut mixed, _, mixed_units) = fixture(
            &mut mixed_heap,
            &[b"C\0", b"C\0", b"N\0", b"O\0"],
            &mixed_specs,
            1,
            2,
        );
        let (result, text) = invoke(&mut mixed_heap, &mut mixed, 1, 1, 0);
        assert_eq!(result, Ok(0));
        let source_unit = &mixed_heap.slice(mixed_units[0].as_const()).unwrap()[0];
        assert_eq!(source_unit.type_, POLYMER_STY_COP as i32);
        assert_eq!(source_unit.subtype, POLYMER_SST_RAN as i32);
        assert_eq!(
            message(&text),
            b"Set copolymer embedding unit mark to COP; Set missing copolymer subtype to RAN"
        );
        assert_eq!(
            mixed_heap.slice(mixed_units[1].as_const()).unwrap()[0].type_,
            POLYMER_STY_MON as i32
        );

        let mut invalid_representation_heap = SourceHeap::default();
        let (mut invalid_representation, polymer, _) =
            fixture(&mut invalid_representation_heap, &[], &[], 1, 0);
        invalid_representation_heap.slice_mut(polymer).unwrap()[0].n = -1;
        let (result, text) = invoke(
            &mut invalid_representation_heap,
            &mut invalid_representation,
            1,
            1,
            0,
        );
        assert_eq!(result, Ok(9035));
        assert_eq!(invalid_representation.valid_polymer, 0);
        assert_eq!(message(&text), b"Invalid kind of polymer representation");

        let mut pseudo_error_heap = SourceHeap::default();
        let (mut pseudo_error, _, _) = fixture(
            &mut pseudo_error_heap,
            &[b"Zz\0"],
            &[spec(
                POLYMER_STY_MON,
                POLYMER_SST_NON,
                POLYMER_CONN_NON,
                &[1],
                None,
                0,
            )],
            1,
            0,
        );
        pseudo_error_heap.slice_mut(pseudo_error.at).unwrap()[0].valence = 2;
        let (result, text) = invoke(&mut pseudo_error_heap, &mut pseudo_error, 1, 1, 0);
        assert_eq!(result, Ok(77));
        assert_eq!(pseudo_error.valid_polymer, 0);
        assert_eq!(message(&text), b"Invalid pseudo element(s) bonding");

        fn pseudo_cap_case(
            polymer_treat: i32,
            fail_after: Option<u64>,
        ) -> (
            SourceHeap,
            ORIG_ATOM_DATA,
            SourceMutPointer<OAD_Polymer>,
            Result<i32, SourceHeapError>,
            [i8; 256],
        ) {
            let mut heap = SourceHeap::default();
            let (mut original, polymer, _) = fixture(
                &mut heap,
                &[b"Zz\0", b"C\0", b"Zz\0"],
                &[spec(
                    POLYMER_STY_SRU,
                    POLYMER_SST_NON,
                    POLYMER_CONN_HT,
                    &[2],
                    Some(&[1, 2, 3, 2]),
                    2,
                )],
                polymer_treat,
                2,
            );
            if let Some(successful) = fail_after {
                heap.fail_after_allocations(successful);
            }
            let (result, text) = invoke(&mut heap, &mut original, 1, 1, 0);
            (heap, original, polymer, result, text)
        }

        let (heap, original, polymer, result, _) = pseudo_cap_case(1, None);
        assert_eq!(result, Ok(0));
        assert_eq!(original.valid_polymer, 1);
        let polymer = &heap.slice(polymer.as_const()).unwrap()[0];
        assert_eq!(polymer.n_pzz, 2);
        assert_eq!(heap.slice(polymer.pzz.as_const()).unwrap(), &[1, 3]);

        let mut old_pzz_heap = SourceHeap::default();
        let (mut old_pzz_original, old_pzz_polymer, _) = fixture(
            &mut old_pzz_heap,
            &[b"Zz\0", b"C\0", b"Zz\0"],
            &[spec(
                POLYMER_STY_SRU,
                POLYMER_SST_NON,
                POLYMER_CONN_HT,
                &[2],
                Some(&[1, 2, 3, 2]),
                2,
            )],
            1,
            2,
        );
        let stale_pzz = old_pzz_heap
            .allocate_model_storage(vec![91_i32, 92])
            .unwrap();
        old_pzz_heap.slice_mut(old_pzz_polymer).unwrap()[0].pzz = stale_pzz;
        let (result, _) = invoke(&mut old_pzz_heap, &mut old_pzz_original, 1, 1, 0);
        assert_eq!(result, Ok(0));
        assert!(matches!(
            old_pzz_heap.slice(stale_pzz.as_const()),
            Err(SourceHeapError::MissingAllocation)
        ));
        let replacement_pzz = old_pzz_heap.slice(old_pzz_polymer.as_const()).unwrap()[0].pzz;
        assert_ne!(replacement_pzz, stale_pzz);
        assert_eq!(
            old_pzz_heap.slice(replacement_pzz.as_const()).unwrap(),
            &[1, 3]
        );

        let (_, original, _, result, text) = pseudo_cap_case(POLYMERS_NO as i32, None);
        assert_eq!(result, Ok(9));
        assert_eq!(original.valid_polymer, 0);
        assert_eq!(message(&text), b"Pseudoelement endgroups are not allowed");

        let (heap, original, polymer, result, text) = pseudo_cap_case(1, Some(0));
        assert_eq!(result, Ok(9010));
        assert_eq!(original.valid_polymer, 0);
        assert_eq!(message(&text), b"Not enough memory");
        assert!(heap.slice(polymer.as_const()).unwrap()[0].pzz.is_null());
    }

    #[test]
    fn source_port__runichi3__oad_polymer_cyclizecloseableunits__line_1979() {
        fn atom(name: &[u8], element: u8) -> inp_ATOM {
            let mut atom = inp_ATOM {
                el_number: element,
                ..inp_ATOM::default()
            };
            for (target, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
                *target = source as i8;
            }
            atom
        }

        fn ring_case(
            heap: &mut SourceHeap,
            metal: bool,
            end_bond: bool,
        ) -> (ORIG_ATOM_DATA, SourceMutPointer<OAD_PolymerUnit>) {
            let atoms = heap
                .allocate_model_storage(vec![
                    atom(b"Zz\0", 0),
                    atom(b"C\0", if metal { 26 } else { 6 }),
                    atom(b"N\0", 7),
                    atom(b"Zz\0", 0),
                ])
                .unwrap();
            let mut bond_count = 0_i32;
            OrigAtData_AddBond(heap, 0, 1, atoms, 1, 0, &mut bond_count).unwrap();
            OrigAtData_AddBond(heap, 3, 2, atoms, 1, 0, &mut bond_count).unwrap();
            if end_bond {
                OrigAtData_AddBond(heap, 1, 2, atoms, 1, 0, &mut bond_count).unwrap();
            }
            let atom_list = heap.allocate_model_storage(vec![2_i32, 3]).unwrap();
            let bond_list = heap.allocate_model_storage(vec![1_i32, 2, 3, 4]).unwrap();
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    na: 2,
                    nb: 2,
                    cyclizable: 91,
                    cyclized: 92,
                    alist: atom_list,
                    blist: bond_list,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            let units = heap.allocate_model_storage(vec![unit]).unwrap();
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer {
                    units,
                    n: 1,
                    ..OAD_Polymer::default()
                }])
                .unwrap();
            (
                ORIG_ATOM_DATA {
                    at: atoms,
                    polymer,
                    num_inp_atoms: 4,
                    num_inp_bonds: bond_count,
                    ..ORIG_ATOM_DATA::default()
                },
                unit,
            )
        }

        fn text(buffer: &[i8]) -> Vec<u8> {
            buffer
                .iter()
                .take_while(|byte| **byte != 0)
                .map(|byte| *byte as u8)
                .collect()
        }

        let mut heap = SourceHeap::default();
        let (mut ring, ring_unit) = ring_case(&mut heap, true, false);
        let mut message = [0_i8; 256];
        assert_eq!(
            OAD_Polymer_CyclizeCloseableUnits(
                &mut heap,
                &mut ring,
                i32::MIN,
                Some(&mut message),
                0,
            ),
            Ok(0)
        );
        assert_eq!(
            text(&message),
            b"Frame shift in metallated polymer unit may be missed"
        );
        let unit = heap.slice(ring_unit.as_const()).unwrap()[0].clone();
        assert_eq!(unit.cyclizable, CLOSING_SRU_RING as i32);
        assert_eq!(unit.cyclized, 1);
        assert_eq!(
            (unit.cap1, unit.end_atom1, unit.end_atom2, unit.cap2),
            (1, 2, 3, 4)
        );
        assert_eq!(ring.num_inp_bonds, 1);
        let atoms = heap.slice(ring.at.as_const()).unwrap();
        assert_eq!((atoms[1].neighbor[0], atoms[2].neighbor[0]), (2, 1));

        let (mut quiet, quiet_unit) = ring_case(&mut heap, true, false);
        message.fill(0);
        assert_eq!(
            OAD_Polymer_CyclizeCloseableUnits(
                &mut heap,
                &mut quiet,
                i32::MAX,
                Some(&mut message),
                1,
            ),
            Ok(0)
        );
        assert!(text(&message).is_empty());
        assert_eq!(heap.slice(quiet_unit.as_const()).unwrap()[0].cyclized, 1);

        let (mut nonmetal, nonmetal_unit) = ring_case(&mut heap, false, false);
        assert_eq!(
            OAD_Polymer_CyclizeCloseableUnits(&mut heap, &mut nonmetal, 0, Some(&mut message), 0,),
            Ok(0)
        );
        assert!(text(&message).is_empty());
        assert_eq!(heap.slice(nonmetal_unit.as_const()).unwrap()[0].cyclized, 1);

        let (mut higher, higher_unit) = ring_case(&mut heap, true, true);
        assert_eq!(
            OAD_Polymer_CyclizeCloseableUnits(&mut heap, &mut higher, 71, None, 0),
            Ok(0)
        );
        let unit = &heap.slice(higher_unit.as_const()).unwrap()[0];
        assert_eq!(unit.cyclizable, CLOSING_SRU_HIGHER_ORDER_BOND as i32);
        assert_eq!(unit.cyclized, 1);
        assert_eq!(higher.num_inp_bonds, 1);
        let atoms = heap.slice(higher.at.as_const()).unwrap();
        assert_eq!((atoms[1].bond_type[0], atoms[2].bond_type[0]), (2, 2));

        let inactive_unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                cyclizable: 0,
                cyclized: 81,
                cap1: 82,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let inactive_units = heap.allocate_model_storage(vec![inactive_unit]).unwrap();
        let inactive_polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units: inactive_units,
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let mut inactive = ORIG_ATOM_DATA {
            polymer: inactive_polymer,
            ..ORIG_ATOM_DATA::default()
        };
        assert_eq!(
            OAD_Polymer_CyclizeCloseableUnits(&mut heap, &mut inactive, -9, None, -1),
            Ok(0)
        );
        let unit = &heap.slice(inactive_unit.as_const()).unwrap()[0];
        assert_eq!((unit.cyclized, unit.cap1), (81, 82));

        for count in [0_i32, -1, i32::MIN] {
            let polymer = heap
                .allocate_model_storage(vec![OAD_Polymer {
                    n: count,
                    ..OAD_Polymer::default()
                }])
                .unwrap();
            let mut empty = ORIG_ATOM_DATA {
                polymer,
                ..ORIG_ATOM_DATA::default()
            };
            assert_eq!(
                OAD_Polymer_CyclizeCloseableUnits(&mut heap, &mut empty, 0, None, 0),
                Ok(0)
            );
        }

        let ordinary_atoms = heap
            .allocate_model_storage(vec![
                atom(b"C\0", 6),
                atom(b"N\0", 7),
                atom(b"O\0", 8),
                atom(b"F\0", 9),
            ])
            .unwrap();
        let invalid_alist = heap.allocate_model_storage(vec![2_i32, 3]).unwrap();
        let invalid_blist = heap.allocate_model_storage(vec![2_i32, 3, 1, 4]).unwrap();
        let invalid_unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 2,
                nb: 2,
                cyclizable: 1,
                cyclized: 33,
                alist: invalid_alist,
                blist: invalid_blist,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let untouched_unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                cyclizable: 1,
                cyclized: 44,
                cap1: 45,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = heap
            .allocate_model_storage(vec![invalid_unit, untouched_unit])
            .unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units,
                n: 2,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let mut invalid = ORIG_ATOM_DATA {
            at: ordinary_atoms,
            polymer,
            num_inp_atoms: 4,
            ..ORIG_ATOM_DATA::default()
        };
        message.fill(0);
        assert_eq!(
            OAD_Polymer_CyclizeCloseableUnits(&mut heap, &mut invalid, 0, Some(&mut message), 0,),
            Ok(9032)
        );
        assert_eq!(text(&message), b"Polymer CRU cap(s) lie inside CRU");
        let invalid_result = &heap.slice(invalid_unit.as_const()).unwrap()[0];
        assert_eq!(invalid_result.cyclizable, CLOSING_SRU_NOT_APPLICABLE as i32);
        assert_eq!(invalid_result.cyclized, 33);
        let untouched = &heap.slice(untouched_unit.as_const()).unwrap()[0];
        assert_eq!((untouched.cyclized, untouched.cap1), (44, 45));

        let mut null_polymer = ORIG_ATOM_DATA::default();
        assert_eq!(
            OAD_Polymer_CyclizeCloseableUnits(&mut heap, &mut null_polymer, 0, None, 0),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_unlinkcapsandconnectendatoms__line_2101() {
        fn fixture(
            heap: &mut SourceHeap,
            cyclizable: i32,
            end_bond_type: Option<i32>,
        ) -> (SourceMutPointer<OAD_PolymerUnit>, ORIG_ATOM_DATA) {
            let atoms = heap
                .allocate_model_storage(vec![inp_ATOM::default(); 4])
                .unwrap();
            let mut bond_count = 0_i32;
            OrigAtData_AddBond(heap, 0, 1, atoms, 1, -1, &mut bond_count).unwrap();
            OrigAtData_AddBond(heap, 3, 2, atoms, 1, 1, &mut bond_count).unwrap();
            if let Some(bond_type) = end_bond_type {
                OrigAtData_AddBond(heap, 1, 2, atoms, bond_type, 0, &mut bond_count).unwrap();
            }
            let unit = heap
                .allocate_model_storage(vec![OAD_PolymerUnit {
                    cyclizable,
                    cyclized: 73,
                    cap1: 1,
                    end_atom1: 2,
                    end_atom2: 3,
                    cap2: 4,
                    ..OAD_PolymerUnit::default()
                }])
                .unwrap();
            (
                unit,
                ORIG_ATOM_DATA {
                    at: atoms,
                    num_inp_atoms: 4,
                    num_inp_bonds: bond_count,
                    ..ORIG_ATOM_DATA::default()
                },
            )
        }

        let mut heap = SourceHeap::default();
        let (inactive_unit, mut inactive_data) = fixture(&mut heap, 0, None);
        let inactive_atoms = heap.slice(inactive_data.at.as_const()).unwrap().to_vec();
        let mut error = 91_i32;
        let mut message = [0_i8; 64];
        message[..5].copy_from_slice(&[b'k' as i8, b'e' as i8, b'e' as i8, b'p' as i8, 0]);
        OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(
            &mut heap,
            inactive_unit,
            &mut inactive_data,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!(error, 0);
        assert_eq!(inactive_data.num_inp_bonds, 2);
        assert_eq!(
            heap.slice(inactive_unit.as_const()).unwrap()[0].cyclized,
            73
        );
        assert_eq!(
            heap.slice(inactive_data.at.as_const()).unwrap(),
            inactive_atoms.as_slice()
        );
        assert_eq!(
            &message[..5],
            &[b'k' as i8, b'e' as i8, b'e' as i8, b'p' as i8, 0]
        );

        let (ring_unit, mut ring_data) = fixture(&mut heap, CLOSING_SRU_RING as i32, None);
        error = i32::MIN;
        OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(
            &mut heap,
            ring_unit,
            &mut ring_data,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!((error, ring_data.num_inp_bonds), (0, 1));
        assert_eq!(heap.slice(ring_unit.as_const()).unwrap()[0].cyclized, 1);
        let atoms = heap.slice(ring_data.at.as_const()).unwrap();
        assert_eq!((atoms[0].valence, atoms[3].valence), (0, 0));
        assert_eq!((atoms[1].valence, atoms[2].valence), (1, 1));
        assert_eq!((atoms[1].neighbor[0], atoms[2].neighbor[0]), (2, 1));
        assert_eq!((atoms[1].bond_type[0], atoms[2].bond_type[0]), (1, 1));
        assert_eq!((atoms[1].bond_stereo[0], atoms[2].bond_stereo[0]), (0, 0));

        let (higher_unit, mut higher_data) =
            fixture(&mut heap, CLOSING_SRU_HIGHER_ORDER_BOND as i32, Some(1));
        OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(
            &mut heap,
            higher_unit,
            &mut higher_data,
            &mut error,
            None,
        )
        .unwrap();
        assert_eq!((error, higher_data.num_inp_bonds), (0, 1));
        let atoms = heap.slice(higher_data.at.as_const()).unwrap();
        assert_eq!((atoms[1].bond_type[0], atoms[2].bond_type[0]), (2, 2));
        assert_eq!(heap.slice(higher_unit.as_const()).unwrap()[0].cyclized, 1);

        let (missing_end_unit, mut missing_end_data) =
            fixture(&mut heap, CLOSING_SRU_HIGHER_ORDER_BOND as i32, None);
        OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(
            &mut heap,
            missing_end_unit,
            &mut missing_end_data,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!((error, missing_end_data.num_inp_bonds), (0, 0));
        assert_eq!(
            heap.slice(missing_end_unit.as_const()).unwrap()[0].cyclized,
            1
        );
        assert_eq!(
            &message[..5],
            &[b'k' as i8, b'e' as i8, b'e' as i8, b'p' as i8, 0]
        );

        let (diradical_unit, mut diradical_data) =
            fixture(&mut heap, CLOSING_SRU_DIRADICAL as i32, None);
        OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(
            &mut heap,
            diradical_unit,
            &mut diradical_data,
            &mut error,
            None,
        )
        .unwrap();
        assert_eq!((error, diradical_data.num_inp_bonds), (0, 0));
        let atoms = heap.slice(diradical_data.at.as_const()).unwrap();
        assert_eq!(atoms[1].radical, RADICAL_TRIPLET as i8);
        assert_eq!(
            (
                atoms[0].valence,
                atoms[1].valence,
                atoms[2].valence,
                atoms[3].valence
            ),
            (0, 0, 0, 0)
        );
        assert_eq!(
            heap.slice(diradical_unit.as_const()).unwrap()[0].cyclized,
            1
        );

        let (unknown_unit, mut unknown_data) = fixture(&mut heap, i32::MAX, Some(3));
        let unknown_atoms = heap.slice(unknown_data.at.as_const()).unwrap().to_vec();
        OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(
            &mut heap,
            unknown_unit,
            &mut unknown_data,
            &mut error,
            Some(&mut message),
        )
        .unwrap();
        assert_eq!((error, unknown_data.num_inp_bonds), (0, 3));
        assert_eq!(heap.slice(unknown_unit.as_const()).unwrap()[0].cyclized, 1);
        assert_eq!(
            heap.slice(unknown_data.at.as_const()).unwrap(),
            unknown_atoms.as_slice()
        );

        let (partial_unit, mut partial_data) = fixture(&mut heap, CLOSING_SRU_RING as i32, None);
        heap.slice_mut(partial_data.at).unwrap()[3].neighbor[0] = 1;
        OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(
            &mut heap,
            partial_unit,
            &mut partial_data,
            &mut error,
            None,
        )
        .unwrap();
        assert_eq!((error, partial_data.num_inp_bonds), (0, 2));
        let atoms = heap.slice(partial_data.at.as_const()).unwrap();
        assert_eq!(atoms[3].valence, 1);
        assert_eq!(atoms[2].neighbor[0], 0);
        assert_eq!(atoms[2].valence, 2);
        assert_eq!((atoms[2].neighbor[1], atoms[2].bond_type[1]), (1, 1));
        assert_eq!((atoms[1].neighbor[0], atoms[1].bond_type[0]), (2, 1));

        let (invalid_unit, mut invalid_data) =
            fixture(&mut heap, CLOSING_SRU_DIRADICAL as i32, None);
        heap.slice_mut(invalid_unit).unwrap()[0].end_atom1 = 0;
        error = 44;
        assert_eq!(
            OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms(
                &mut heap,
                invalid_unit,
                &mut invalid_data,
                &mut error,
                None,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(error, 0);
        assert_eq!(heap.slice(invalid_unit.as_const()).unwrap()[0].cyclized, 73);
    }

    #[test]
    fn source_port__runichi3__oad_polymerunit_hasmetal__line_2055() {
        let mut heap = SourceHeap::default();
        let empty = OAD_PolymerUnit::default();
        assert_eq!(
            OAD_PolymerUnit_HasMetal(&heap, &empty, SourceConstPointer::null()),
            Ok(0)
        );
        let negative = OAD_PolymerUnit {
            na: i32::MIN,
            ..OAD_PolymerUnit::default()
        };
        assert_eq!(
            OAD_PolymerUnit_HasMetal(&heap, &negative, SourceConstPointer::null()),
            Ok(0)
        );

        for periodic_number in 0_u8..=121 {
            let atoms = heap
                .allocate_model_storage(vec![inp_ATOM {
                    el_number: periodic_number,
                    ..inp_ATOM::default()
                }])
                .unwrap();
            let atom_list = heap.allocate_model_storage(vec![1_i32]).unwrap();
            let unit = OAD_PolymerUnit {
                na: 1,
                alist: atom_list,
                ..OAD_PolymerUnit::default()
            };
            let expected = i32::from(matches!(
                periodic_number,
                3..=4 | 11..=13 | 19..=31 | 37..=51 | 55..=84 | 87..=118
            ));
            assert_eq!(
                OAD_PolymerUnit_HasMetal(&heap, &unit, atoms.as_const()),
                Ok(expected),
                "periodic number {periodic_number}"
            );
        }

        let atoms = heap
            .allocate_model_storage(
                [6_u8, 7, 26, 8, 29]
                    .map(|el_number| inp_ATOM {
                        el_number,
                        ..inp_ATOM::default()
                    })
                    .to_vec(),
            )
            .unwrap();
        for atom_list in [vec![3, 1, 2], vec![1, 3, 2], vec![1, 2, 5]] {
            let list = heap.allocate_model_storage(atom_list).unwrap();
            let unit = OAD_PolymerUnit {
                na: 3,
                alist: list,
                ..OAD_PolymerUnit::default()
            };
            assert_eq!(
                OAD_PolymerUnit_HasMetal(&heap, &unit, atoms.as_const()),
                Ok(1)
            );
        }
        let nonmetals = heap.allocate_model_storage(vec![1, 2, 4]).unwrap();
        assert_eq!(
            OAD_PolymerUnit_HasMetal(
                &heap,
                &OAD_PolymerUnit {
                    na: 3,
                    alist: nonmetals,
                    ..OAD_PolymerUnit::default()
                },
                atoms.as_const(),
            ),
            Ok(0)
        );

        let early_metal = heap.allocate_model_storage(vec![3, 0]).unwrap();
        assert_eq!(
            OAD_PolymerUnit_HasMetal(
                &heap,
                &OAD_PolymerUnit {
                    na: 2,
                    alist: early_metal,
                    ..OAD_PolymerUnit::default()
                },
                atoms.as_const(),
            ),
            Ok(1)
        );
        let invalid_index = heap.allocate_model_storage(vec![0]).unwrap();
        assert_eq!(
            OAD_PolymerUnit_HasMetal(
                &heap,
                &OAD_PolymerUnit {
                    na: 1,
                    alist: invalid_index,
                    ..OAD_PolymerUnit::default()
                },
                atoms.as_const(),
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let short_list = heap.allocate_model_storage(vec![1]).unwrap();
        assert_eq!(
            OAD_PolymerUnit_HasMetal(
                &heap,
                &OAD_PolymerUnit {
                    na: 2,
                    alist: short_list,
                    ..OAD_PolymerUnit::default()
                },
                atoms.as_const(),
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }
}
