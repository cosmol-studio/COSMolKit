use crate::source::base::ichi_bns::mark_alt_bonds_and_taut_groups;
use crate::source::base::ichi_io::inchi_ios_init;
use crate::source::base::ichican2::{DeAllocBCN, GetBaseCanonRanking};
use crate::source::base::ichicano::{
    AllocateCS, Canon_INChI, DeAllocateCS, GetCanonLengths, InchiTimeAddMsec, InchiTimeElapsed,
    InchiTimeGet,
};
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichiisot::set_atom_iso_sort_keys;
use crate::source::base::ichimak2::{FillOutINChIBehavior, FillOutINChIWithBehavior};
use crate::source::base::ichimake::{CheckCanonNumberingCorrectness, inp2spATOM};
use crate::source::base::ichinorm::{FreeInpAtomData, MarkRingSystemsInp};
use crate::source::base::ichister::set_stereo_parity;
use crate::source::base::ichitaut::{
    CountTautomerGroups, free_t_group_info, make_a_copy_of_t_group_info, set_tautomer_iso_sort_keys,
};
use crate::source::base::mol2atom::{
    CreateInpAtom, CreateInpAtomData, FreeCompAtomData, FreeInpAtom,
};
use crate::source::base::runichi2::{GetOneComponent, TreatErrorsInReadTheStructure};
use crate::source::base::runichi3::PreprocessOneStructure;
use crate::source::base::runichi4::{
    GetProcessingWarningsOneComponentInChI, TreatErrorsInCreateOneComponentINChI,
};
use crate::source::base::strutil::{
    Alloc_INChI, Alloc_INChI_Aux, SetConnectedComponentNumber, add_DT_to_num_H, remove_terminal_HDT,
};
use crate::source::base::util::{inchi_calloc, inchi_free};
use crate::source_types::{
    _IS_ERROR, _IS_FATAL, AB_PARITY_UNDF, AB_PARITY_UNKN, AT_NUMB, ATT_PROTON, BCN, CANON_GLOBALS,
    CANON_MODE_CT, CANON_MODE_ISO, CANON_MODE_ISO_STEREO, CANON_MODE_STEREO, CANON_MODE_TAUT,
    CANON_STAT, CMODE_NO_ALT_SBONDS, CMODE_NOEQ_STEREO, CMODE_RACEMIC_STEREO, CMODE_REDNDNT_STEREO,
    CMODE_RELATIVE_STEREO, CMODE_SB_IGN_ALL_UU, CMODE_SC_IGN_ALL_UU, COMP_ATOM_DATA,
    COMPONENT_TREAT_INFO, CT_ERR_MAX, CT_ERR_MIN, CT_OUT_OF_RAM, CT_TAUCOUNT_ERR, CT_USER_QUIT_ERR,
    FIX_ISO_FIXEDH_BUG, FIX_TERM_H_CHRG_BUG, FLAG_NORM_CONSIDER_TAUT, INCHI_BAS, INCHI_CLOCK,
    INCHI_IOS_TYPE_FILE, INCHI_MODE, INCHI_NUM, INCHI_OUT_NO_AUX_INFO, INCHI_OUT_SHORT_AUX_INFO,
    INCHI_REC, INCHIGEN_CONTROL, INCHIGEN_DATA, INChI, INChI_Aux, INP_ATOM_DATA, INP_ATOM_DATA2,
    LOG_MASK_ALL, NORM_ATOMS, NUM_H_ISOTOPES, PES_BIT_ARSINE_STEREO, PES_BIT_FIX_SP3_BUG,
    PES_BIT_PHOSPHINE_STEREO, PES_BIT_POINT_EDGE_STEREO, REQ_MODE_BASIC, REQ_MODE_DEFAULT,
    REQ_MODE_DIFF_UU_STEREO, REQ_MODE_ISO, REQ_MODE_ISO_STEREO, REQ_MODE_NO_ALT_SBONDS,
    REQ_MODE_NOEQ_STEREO, REQ_MODE_NON_ISO, REQ_MODE_RACEMIC_STEREO, REQ_MODE_REDNDNT_STEREO,
    REQ_MODE_RELATIVE_STEREO, REQ_MODE_SB_IGN_ALL_UU, REQ_MODE_SC_IGN_ALL_UU, REQ_MODE_STEREO,
    REQ_MODE_TAUT, SourceHeap, SourceHeapError, SourceMutPointer, TAUT_INI, TAUT_NON, TAUT_NUM,
    TAUT_YES, TG_FLAG_ALL_TAUTOMERIC, TG_FLAG_ARSINE_STEREO, TG_FLAG_FIX_ISO_FIXEDH_BUG,
    TG_FLAG_FIX_SP3_BUG, TG_FLAG_FIX_TERM_H_CHRG_BUG, TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE,
    TG_FLAG_FOUND_ISOTOPIC_H_DONE, TG_FLAG_H_ALREADY_REMOVED, TG_FLAG_PHOSPHINE_STEREO,
    TG_FLAG_POINTED_EDGE_STEREO, clock_t, inchiTime, inp_ATOM, sp_ATOM,
};

const SOURCE_SIZEOF_AT_NUMB: u64 = 2;
const SOURCE_SIZEOF_INP_ATOM_DATA2: u64 = 192;
const SOURCE_SIZEOF_PINCHI2: u64 = 16;
const SOURCE_SIZEOF_NORM_ATOMS: u64 = 96;
const SOURCE_SIZEOF_COMPONENT_TREAT_INFO: u64 = 848;

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Normalization_step(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    ic: SourceMutPointer<INCHI_CLOCK>,
    pp_inchi: &[SourceMutPointer<INChI>; TAUT_NUM as usize],
    pp_inchi_aux: &[SourceMutPointer<INChI_Aux>; TAUT_NUM as usize],
    inp_at: SourceMutPointer<inp_ATOM>,
    out_norm_data: &mut [INP_ATOM_DATA; TAUT_NUM as usize],
    num_inp_at: i32,
    max_time: SourceMutPointer<inchiTime>,
    taut_flags: &mut INCHI_MODE,
    taut_flags_done: &mut INCHI_MODE,
    z: &mut COMPONENT_TREAT_INFO,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1138 Normalization_step
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int  Normalization_step( CANON_GLOBALS *pCG,
                         INCHI_CLOCK *ic,
                         INChI **ppINChI,
                         INChI_Aux **ppINChI_Aux,
                         inp_ATOM *inp_at,
                         INP_ATOM_DATA *out_norm_data[2],
                         int num_inp_at,
                         struct tagInchiTime *ulMaxTime,
                         INCHI_MODE *pbTautFlags,
                         INCHI_MODE *pbTautFlagsDone,
                         COMPONENT_TREAT_INFO *z )
{

    int i, ret = 0;

    T_GROUP_INFO * /*const*/  t_group_info = &( z->vt_group_info );
    T_GROUP_INFO * /*const*/  t_group_info_orig = &( z->vt_group_info_orig );

    BCN *pBCN = &( z->Bcn );

    /* */
    z->fix_isofixedh = 0;
    z->fix_termhchrg = 0;
    /* */
#if( FIX_ISO_FIXEDH_BUG == 1 )
    if (TG_FLAG_FIX_ISO_FIXEDH_BUG & *pbTautFlags)
        z->fix_isofixedh = 1;
#endif
#if( FIX_TERM_H_CHRG_BUG == 1 )
    if (TG_FLAG_FIX_TERM_H_CHRG_BUG & *pbTautFlags)
        z->fix_termhchrg = 1;
#endif


    z->bPointedEdgeStereo = ( ( TG_FLAG_POINTED_EDGE_STEREO & *pbTautFlags ) ? PES_BIT_POINT_EDGE_STEREO : 0 )
        | ( ( TG_FLAG_PHOSPHINE_STEREO    & *pbTautFlags ) ? PES_BIT_PHOSPHINE_STEREO : 0 )
        | ( ( TG_FLAG_ARSINE_STEREO       & *pbTautFlags ) ? PES_BIT_ARSINE_STEREO : 0 )
        | ( ( TG_FLAG_FIX_SP3_BUG         & *pbTautFlags ) ? PES_BIT_FIX_SP3_BUG : 0 );
    z->bTautFlags = ( *pbTautFlags     & ( ~(INCHI_MODE) TG_FLAG_ALL_TAUTOMERIC ) );
    z->bTautFlagsDone = ( *pbTautFlagsDone /*& (~(INCHI_MODE)TG_FLAG_ALL_TAUTOMERIC) */ );

    z->out_at = NULL;       /*, *norm_at_fixed_bonds[TAUT_NUM]; */ /*  = {out_norm_nontaut_at, out_norm_taut_at} ; */

    /* Init: internal structs */

    memset( z->s, 0, sizeof( z->s ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    if (pBCN) 
        memset( pBCN, 0, sizeof( pBCN[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    memset( t_group_info, 0, sizeof( *t_group_info ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( t_group_info_orig, 0, sizeof( *t_group_info_orig ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* Allocate: at[] */

    for (i = 0; i < TAUT_NUM; i++)
    {
        if (out_norm_data[i]->at)
        {
            z->at[i] =
                (sp_ATOM  *) inchi_malloc( num_inp_at * sizeof( *( z->at[0] ) ) );

            if (!z->at[i])
            {
                ret = -1;
            }
        }
        else
        {
            z->at[i] = NULL;
        }
    }

    if ((!out_norm_data[TAUT_NON]->at && !out_norm_data[TAUT_YES]->at) || !inp_at || ret) /* djb-rwth: addressing LLVM warning */
    {
        ret = -1;
        goto exit_function;
    }

    /* the first struct to process: tautomeric if exists else non-tautomeric */

    z->out_at = out_norm_data[TAUT_YES]->at ? out_norm_data[TAUT_YES]->at : out_norm_data[TAUT_NON]->at;

    /* copy the input structure to be normalized to the buffer for the normalization data */

    memcpy( z->out_at, inp_at, num_inp_at * sizeof( z->out_at[0] ) );

    /*  tautomeric groups setting */

    t_group_info->bIgnoreIsotopic = 0;   /*  include tautomeric group isotopic info in MarkTautomerGroups() */
    t_group_info->bTautFlags = *pbTautFlags;
    t_group_info->bTautFlagsDone = *pbTautFlagsDone;


    /*  Preprocess the structure; here THE NUMBER OF ATOMS MAY BE REDUCED */
    /*  ??? Ambiguity: H-D may become HD or DH (that is, H+implicit D or D+implicit H) */

    if (TG_FLAG_H_ALREADY_REMOVED & z->bTautFlags)
    {
        INP_ATOM_DATA *out_norm_data1 = out_norm_data[TAUT_YES]->at ? out_norm_data[TAUT_YES] :
            out_norm_data[TAUT_NON]->at ? out_norm_data[TAUT_NON] : NULL;
        if (out_norm_data1)
        {
            z->num_at_tg =
                z->num_atoms = out_norm_data1->num_at - out_norm_data1->num_removed_H;
            z->num_deleted_H = out_norm_data1->num_removed_H;
            t_group_info->tni.nNumRemovedExplicitH = z->num_deleted_H;
        }
        else
        {
            ret = -1;
            goto exit_function;
        }
    }
    else
    {
        z->num_at_tg = z->num_atoms = remove_terminal_HDT( num_inp_at, z->out_at, z->fix_termhchrg );

        z->num_deleted_H = num_inp_at - z->num_atoms;
        t_group_info->tni.nNumRemovedExplicitH = z->num_deleted_H;

        add_DT_to_num_H( z->num_atoms, z->out_at );
    }

    /* fix_odd_things( z->num_atoms, z->out_at );*/

#if( FIND_RING_SYSTEMS == 1 )
    MarkRingSystemsInp( z->out_at, z->num_atoms, 0 );
#endif

    /*  duplicate the preprocessed structure so that all supplied out_norm_data[]->at buffers are filled */
    if (z->out_at != out_norm_data[TAUT_YES]->at && out_norm_data[TAUT_YES]->at)
    {
        memcpy( out_norm_data[TAUT_YES]->at, z->out_at, num_inp_at * sizeof( z->out_at[0] ) );
    }

    if (out_norm_data[TAUT_YES]->at_fixed_bonds && out_norm_data[TAUT_YES]->at)
    {
        memcpy( out_norm_data[TAUT_YES]->at_fixed_bonds, z->out_at, num_inp_at * sizeof( z->out_at[0] ) );
    }

    if (z->out_at != out_norm_data[TAUT_NON]->at && out_norm_data[TAUT_NON]->at)
    {
        memcpy( out_norm_data[TAUT_NON]->at, z->out_at, num_inp_at * sizeof( z->out_at[0] ) );
    }

    /*******************************************************************************
     * ??? not true ??? duplicate inp_at and keep inp_at[] unchanged after terminal hydrogens removal
     * set stereo parities in taut_at[], non_taut_at[]
     * obtain max. lenghts of the name stereo parts
     * Ignore absence/presence of isotopic stereo for now
     * mark isotopic atoms
     *******************************************************************************/

    if (out_norm_data[TAUT_YES]->at && z->at[TAUT_YES])
    {

        /* final normalization of possibly tautomeric structure */

        ret = mark_alt_bonds_and_taut_groups( ic, pCG,
                                               out_norm_data[TAUT_YES]->at,
                                               out_norm_data[TAUT_YES]->at_fixed_bonds,
                                               z->num_atoms, ulMaxTime, t_group_info,
                                               NULL, NULL, 0, NULL );

        if (ret < 0)
        {
            goto exit_function;/*  out of RAM or other normalization problem */
        }

        z->num_taut_at = ret; /* number of atoms without removed H? */
        z->num_deleted_H_taut = t_group_info->tni.nNumRemovedExplicitH;
        out_norm_data[TAUT_YES]->num_at = z->num_atoms + z->num_deleted_H_taut; /* protons might have been removed */
        out_norm_data[TAUT_YES]->num_removed_H = z->num_deleted_H_taut;
        out_norm_data[TAUT_YES]->nNumRemovedProtons += t_group_info->tni.nNumRemovedProtons;

        for (i = 0; i < NUM_H_ISOTOPES; i++)
        {
            out_norm_data[TAUT_YES]->nNumRemovedProtonsIsotopic[i] += t_group_info->tni.nNumRemovedProtonsIsotopic[i] /*+ t_group_info->num_iso_H[i]*/;
            out_norm_data[TAUT_YES]->num_iso_H[i] += t_group_info->num_iso_H[i];
        }

        /* mark deleted isolated tautomeric H(+) */

        if (z->num_taut_at == 1 && out_norm_data[TAUT_YES]->at[0].at_type == ATT_PROTON &&
             t_group_info && t_group_info->tni.nNumRemovedProtons == 1)
        {
            out_norm_data[TAUT_YES]->bDeleted = 1;

            FreeInpAtom( &out_norm_data[TAUT_YES]->at_fixed_bonds );
        }
        else if (( t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT ) &&
                   out_norm_data[TAUT_YES]->at_fixed_bonds)
        {
            out_norm_data[TAUT_YES]->bTautPreprocessed = 1;
        }

        out_norm_data[TAUT_YES]->bTautFlags = *pbTautFlags = t_group_info->bTautFlags;
        out_norm_data[TAUT_YES]->bTautFlagsDone = *pbTautFlagsDone = t_group_info->bTautFlagsDone;
        out_norm_data[TAUT_YES]->bNormalizationFlags = t_group_info->tni.bNormalizationFlags;

        /* create internal sp_ATOM at[] out of out_norm_data[]->at */

        inp2spATOM( out_norm_data[TAUT_YES]->at, num_inp_at, z->at[TAUT_YES] );

        /* set stereo parities to at[]; nUserMode: accept alt. stereo bonds, min ring size */

        ret = set_stereo_parity( pCG,
                                 out_norm_data[TAUT_YES]->at,
                                 z->at[TAUT_YES],
                                 z->num_taut_at,
                                 z->num_deleted_H_taut,
                                 &( z->s[TAUT_YES].nMaxNumStereoAtoms ),
                                 &( z->s[TAUT_YES].nMaxNumStereoBonds ),
                                 z->nUserMode,
                                 z->bPointedEdgeStereo,
                                 z->vABParityUnknown,
                                 z->bLooseTSACheck,
                                 z->bStereoAtZz );

        if (RETURNED_ERROR( ret ))
        {
            goto exit_function; /*  stereo bond error */
        }

        z->s[TAUT_YES].bMayHaveStereo = ( z->s[TAUT_YES].nMaxNumStereoAtoms || z->s[TAUT_YES].nMaxNumStereoBonds );

        /*
         * mark isotopic atoms and atoms that have non-tautomeric
         * isotopic terminal hydrogen atoms 1H, 2H(D), 3H(T)
         */

        z->s[TAUT_YES].num_isotopic_atoms =

            set_atom_iso_sort_keys( z->num_taut_at, z->at[TAUT_YES], t_group_info,
                                    &( z->s[TAUT_YES].bHasIsotopicTautGroups ) );


        /**************************************************************************
         *  prepare tautomeric (if no tautomerism found then prepare non-tautomeric)
         *  structure for canonicalizaton:
         **************************************************************************
         *   remove t-groups that have no H,
         *   remove charges from t-groups if requested
         *   renumber t-groups and find final t_group_info->num_t_groups
         *   add to t-groups lists of endpoints tgroup->nEndpointAtomNumber[]
         *   calculate length of the t-group part of the connection table
         **************************************************************************/

        z->s[TAUT_YES].nLenLinearCTTautomer =

            CountTautomerGroups( z->at[TAUT_YES], z->num_taut_at, t_group_info );


        if (RETURNED_ERROR( z->s[TAUT_YES].nLenLinearCTTautomer ))
        {
            /* added error treatment 9-11-2003 */
            ret = z->s[TAUT_YES].nLenLinearCTTautomer;
            goto exit_function;
            /*  error has happened; no breakpoint here
            z->s[TAUT_YES].nLenLinearCTTautomer = 0;
            */
        }
        else if (z->s[TAUT_YES].nLenLinearCTTautomer > 0)
        {
            z->num_at_tg = z->num_taut_at + t_group_info->num_t_groups;

            /*  ??? -not true- create t_group_info_orig for multiple calls with atom renumbering */

            make_a_copy_of_t_group_info( t_group_info_orig /* dest*/, t_group_info /* source*/ );

            /*  mark isotopic tautomer groups: calculate t_group->iWeight */
            z->s[TAUT_YES].nLenLinearCTIsotopicTautomer = set_tautomer_iso_sort_keys( t_group_info );
            if (z->s[TAUT_YES].nLenLinearCTIsotopicTautomer < 0)
            {
                /* ??? -error cannot happen- error has happened; no breakpoint here */
                z->s[TAUT_YES].nLenLinearCTIsotopicTautomer = 0;
            }
            out_norm_data[TAUT_YES]->bTautomeric = z->s[TAUT_YES].nLenLinearCTTautomer;
        }

        /*  new variable: z->s[TAUT_YES].nLenCT introduced 7-22-2002 */

        GetCanonLengths( z->num_taut_at, z->at[TAUT_YES], &( z->s[TAUT_YES] ), t_group_info );
    } /* end of: final normalization of possibly tautomeric structure */



    if (out_norm_data[TAUT_NON]->at && out_norm_data[TAUT_YES]->at && z->at[TAUT_NON] && !z->s[TAUT_YES].nLenLinearCTTautomer)
    {
        /* the structure is non-tautomeric: use tautomeric treatment results only for it */

        inchi_free( z->at[TAUT_NON] );

        z->at[TAUT_NON] = NULL;
    }

    else if (!out_norm_data[TAUT_NON]->at && out_norm_data[TAUT_YES]->at &&
         !z->at[TAUT_NON] && z->at[TAUT_YES] && !z->s[TAUT_YES].nLenLinearCTTautomer)
    {
        /* requested tautomeric; found non-tautomeric; it is located in out_norm_data[TAUT_YES]->at */

        out_norm_data[TAUT_YES]->bTautomeric = 0;
    }

    else if (out_norm_data[TAUT_NON]->at && z->at[TAUT_NON])
    {
        /* the structure needs non-tautomeric treatment: final normalization of non-tautomeric structure */

        ret = mark_alt_bonds_and_taut_groups( ic, pCG,
                                              out_norm_data[TAUT_NON]->at,
                                              NULL, z->num_atoms, ulMaxTime, NULL,
                                              &( z->bTautFlags ), &( z->bTautFlagsDone ),
                                              0, NULL );

        if (ret < 0)
        {
            goto exit_function;  /*  out of RAM or other normalization problem */
        }

        out_norm_data[TAUT_NON]->num_at = z->num_atoms + z->num_deleted_H;
        out_norm_data[TAUT_NON]->num_removed_H = z->num_deleted_H;
        out_norm_data[TAUT_NON]->bTautFlags = *pbTautFlags;
        out_norm_data[TAUT_NON]->bTautFlagsDone = *pbTautFlagsDone;
        out_norm_data[TAUT_NON]->bNormalizationFlags = 0;

        /* create internal sp_ATOM at[] out of out_norm_data[]->at */

        inp2spATOM( out_norm_data[TAUT_NON]->at, num_inp_at, z->at[TAUT_NON] );

        /* set stereo parities to at[]; nUserMode: accept alt. stereo bonds, min ring size */

        ret = set_stereo_parity( pCG,
                                 out_norm_data[TAUT_NON]->at,
                                 z->at[TAUT_NON],
                                 z->num_atoms, z->num_deleted_H,
                                 &( z->s[TAUT_NON].nMaxNumStereoAtoms ),
                                 &( z->s[TAUT_NON].nMaxNumStereoBonds ),
                                 z->nUserMode,
                                 z->bPointedEdgeStereo, z->vABParityUnknown,
                                 z->bLooseTSACheck,
                                 z->bStereoAtZz );

        if (RETURNED_ERROR( ret ))
        {
            goto exit_function; /*  stereo bond error */
        }


        z->s[TAUT_NON].bMayHaveStereo = ( z->s[TAUT_NON].nMaxNumStereoAtoms || z->s[TAUT_NON].nMaxNumStereoBonds );

        /*
         * mark isotopic atoms and atoms that have non-tautomeric
         * isotopic terminal hydrogen atoms 1H, 2H(D), 3H(T)
         */

        z->s[TAUT_NON].num_isotopic_atoms =

            set_atom_iso_sort_keys( z->num_atoms, z->at[TAUT_NON], NULL, NULL );


        GetCanonLengths( z->num_atoms, z->at[TAUT_NON], &( z->s[TAUT_NON] ), NULL );


        out_norm_data[TAUT_NON]->bTautomeric = 0;

    } /* the structure needs non-tautomeric treatment: final normalization of non-tautomeric structure */


    /**********************************************************/

    /*  common */
    z->bMayHaveStereo = z->s[TAUT_YES].bMayHaveStereo || z->s[TAUT_NON].bMayHaveStereo;
    z->bHasIsotopicAtoms = z->s[TAUT_NON].num_isotopic_atoms > 0 || z->s[TAUT_NON].bHasIsotopicTautGroups > 0 ||
        z->s[TAUT_YES].num_isotopic_atoms > 0 || z->s[TAUT_YES].bHasIsotopicTautGroups > 0;
    /* */
    if (z->fix_isofixedh)
    {
        /* 2008-03-21 DT */
        z->bHasIsotopicAtoms = z->bHasIsotopicAtoms ||
            (z->s[TAUT_YES].nLenLinearCTTautomer > 0 && t_group_info &&
            ( (0 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[0]) ||
             (1 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[1]) ||
             (2 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[2]) )); /* djb-rwth: addressing LLVM warnings */
    }
    /* */
    z->bHasIsotopicAtoms = z->bHasIsotopicAtoms ||
        (z->s[TAUT_YES].nLenIsotopicEndpoints > 1 && t_group_info &&
        ( t_group_info->bTautFlagsDone & ( TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE ) )); /* djb-rwth: addressing LLVM warning */


    /* Set mode */

    /*  default mode */

    if (!( z->nUserMode & REQ_MODE_DEFAULT ))
    {
        z->nUserMode |= REQ_MODE_DEFAULT;
    }

    /*  adjust the mode to the reality */

    if (( z->nUserMode & REQ_MODE_ISO ) && !z->bHasIsotopicAtoms)
    {
        z->nUserMode ^= REQ_MODE_ISO;
        z->nUserMode |= REQ_MODE_NON_ISO;  /*  at least one is needed */
    }

    if (( z->nUserMode & REQ_MODE_STEREO ) && ( z->nUserMode & REQ_MODE_ISO ))
    {
        z->nUserMode |= REQ_MODE_ISO_STEREO;
    }

    if (( z->nUserMode & REQ_MODE_STEREO ) && !( z->nUserMode & REQ_MODE_NON_ISO ))
    {
        z->nUserMode ^= REQ_MODE_STEREO;
    }

    if (!z->bMayHaveStereo)
    {
        if (z->nUserMode & REQ_MODE_STEREO)
        {
            z->nUserMode ^= REQ_MODE_STEREO;
        }
        if (z->nUserMode & REQ_MODE_ISO_STEREO)
        {
            z->nUserMode ^= REQ_MODE_ISO_STEREO;
        }
    }

    if (( z->nUserMode & REQ_MODE_BASIC ) && ( !out_norm_data[TAUT_NON]->at || !ppINChI[TAUT_NON] || !ppINChI_Aux[TAUT_NON] || !z->at[TAUT_NON] ))
    {
        z->nUserMode ^= REQ_MODE_BASIC;
    }
    if (( z->nUserMode & REQ_MODE_TAUT ) && ( !out_norm_data[TAUT_YES]->at || !ppINChI[TAUT_YES] || !ppINChI_Aux[TAUT_YES] || !z->at[TAUT_YES] ))
    {
        z->nUserMode ^= REQ_MODE_TAUT;
    }

    /* Set n1, n2 according to the mode */

    switch ((int) z->nUserMode & ( REQ_MODE_BASIC | REQ_MODE_TAUT ))
    {
        case REQ_MODE_BASIC:
            z->n1 = TAUT_NON;
            z->n2 = TAUT_NON;
            break;
        case REQ_MODE_TAUT:
            z->n1 = TAUT_YES;
            z->n2 = TAUT_YES;
            break;
        case ( REQ_MODE_BASIC | REQ_MODE_TAUT ):
            z->n1 = TAUT_NON;
            z->n2 = TAUT_YES;
            break;
        default:
            /*  program error: inconsistent nUserMode or missing taut/non-taut allocation */ /*   <BRKPT> */
            ret = -3;
            goto exit_function;
    }

    if (ret == 0)
    {
        ret = z->num_atoms;
    }
    /*  treat the results later */

exit_function:

    return ret;
}
    */
    // END INCHI C FUNCTION: Normalization_step

    let mut ret = 0_i32;
    z.fix_isofixedh = 0;
    z.fix_termhchrg = 0;
    if FIX_ISO_FIXEDH_BUG == 1 && *taut_flags & u64::from(TG_FLAG_FIX_ISO_FIXEDH_BUG) != 0 {
        z.fix_isofixedh = 1;
    }
    if FIX_TERM_H_CHRG_BUG == 1 && *taut_flags & u64::from(TG_FLAG_FIX_TERM_H_CHRG_BUG) != 0 {
        z.fix_termhchrg = 1;
    }
    z.bPointedEdgeStereo = i32::from(*taut_flags & u64::from(TG_FLAG_POINTED_EDGE_STEREO) != 0)
        * PES_BIT_POINT_EDGE_STEREO as i32
        | i32::from(*taut_flags & u64::from(TG_FLAG_PHOSPHINE_STEREO) != 0)
            * PES_BIT_PHOSPHINE_STEREO as i32
        | i32::from(*taut_flags & u64::from(TG_FLAG_ARSINE_STEREO) != 0)
            * PES_BIT_ARSINE_STEREO as i32
        | i32::from(*taut_flags & u64::from(TG_FLAG_FIX_SP3_BUG) != 0) * PES_BIT_FIX_SP3_BUG as i32;
    z.bTautFlags = *taut_flags & !u64::from(TG_FLAG_ALL_TAUTOMERIC);
    z.bTautFlagsDone = *taut_flags_done;
    z.out_at = SourceMutPointer::null();
    z.s = std::array::from_fn(|_| Default::default());
    z.Bcn = Default::default();
    z.vt_group_info = Default::default();
    z.vt_group_info_orig = Default::default();

    let count = usize::try_from(num_inp_at).ok();
    for index in 0..TAUT_NUM as usize {
        z.at[index] = if !out_norm_data[index].at.is_null() {
            match count {
                Some(count) => {
                    let mut atoms = Vec::new();
                    if atoms.try_reserve_exact(count).is_err() {
                        ret = -1;
                        SourceMutPointer::null()
                    } else {
                        atoms.resize(count, sp_ATOM::default());
                        match heap.allocate(atoms) {
                            Ok(pointer) => pointer,
                            Err(SourceHeapError::AllocationFailed) => {
                                ret = -1;
                                SourceMutPointer::null()
                            }
                            Err(error) => return Err(error),
                        }
                    }
                }
                None => {
                    ret = -1;
                    SourceMutPointer::null()
                }
            }
        } else {
            SourceMutPointer::null()
        };
    }
    if (out_norm_data[TAUT_NON as usize].at.is_null()
        && out_norm_data[TAUT_YES as usize].at.is_null())
        || inp_at.is_null()
        || ret != 0
    {
        return Ok(-1);
    }

    z.out_at = if !out_norm_data[TAUT_YES as usize].at.is_null() {
        out_norm_data[TAUT_YES as usize].at
    } else {
        out_norm_data[TAUT_NON as usize].at
    };
    let count = count.ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let input_atoms = heap
        .slice(inp_at.as_const())?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    heap.slice_mut(z.out_at)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone_from_slice(&input_atoms);
    z.vt_group_info.bIgnoreIsotopic = 0;
    z.vt_group_info.bTautFlags = *taut_flags;
    z.vt_group_info.bTautFlagsDone = *taut_flags_done;

    if z.bTautFlags & u64::from(TG_FLAG_H_ALREADY_REMOVED) != 0 {
        let source = if !out_norm_data[TAUT_YES as usize].at.is_null() {
            &out_norm_data[TAUT_YES as usize]
        } else if !out_norm_data[TAUT_NON as usize].at.is_null() {
            &out_norm_data[TAUT_NON as usize]
        } else {
            return Ok(-1);
        };
        z.num_atoms = source.num_at.wrapping_sub(source.num_removed_H);
        z.num_at_tg = z.num_atoms;
        z.num_deleted_H = source.num_removed_H;
        z.vt_group_info.tni.nNumRemovedExplicitH = z.num_deleted_H as i16;
    } else {
        z.num_atoms = remove_terminal_HDT(heap, num_inp_at, z.out_at, z.fix_termhchrg)?;
        z.num_at_tg = z.num_atoms;
        z.num_deleted_H = num_inp_at.wrapping_sub(z.num_atoms);
        z.vt_group_info.tni.nNumRemovedExplicitH = z.num_deleted_H as i16;
        add_DT_to_num_H(z.num_atoms, heap.slice_mut(z.out_at)?)?;
    }
    let _ = MarkRingSystemsInp(heap, z.out_at, z.num_atoms, 0)?;

    for index in [TAUT_YES as usize, TAUT_NON as usize] {
        if z.out_at != out_norm_data[index].at && !out_norm_data[index].at.is_null() {
            let source = heap.slice(z.out_at.as_const())?[..count].to_vec();
            heap.slice_mut(out_norm_data[index].at)?[..count].clone_from_slice(&source);
        }
    }
    if !out_norm_data[TAUT_YES as usize].at_fixed_bonds.is_null()
        && !out_norm_data[TAUT_YES as usize].at.is_null()
    {
        let source = heap.slice(z.out_at.as_const())?[..count].to_vec();
        heap.slice_mut(out_norm_data[TAUT_YES as usize].at_fixed_bonds)?[..count]
            .clone_from_slice(&source);
    }

    if !out_norm_data[TAUT_YES as usize].at.is_null() && !z.at[TAUT_YES as usize].is_null() {
        let t_group_pointer = heap.allocate_model_storage(vec![z.vt_group_info.clone()])?;
        let result = mark_alt_bonds_and_taut_groups(
            heap,
            ic,
            pCG,
            out_norm_data[TAUT_YES as usize].at,
            out_norm_data[TAUT_YES as usize].at_fixed_bonds,
            z.num_atoms,
            max_time,
            t_group_pointer,
            SourceMutPointer::null(),
            SourceMutPointer::null(),
            0,
            SourceMutPointer::null(),
            clock_result,
        );
        z.vt_group_info = heap.slice(t_group_pointer.as_const())?[0].clone();
        heap.free(t_group_pointer)?;
        ret = result?;
        if ret < 0 {
            return Ok(ret);
        }
        z.num_taut_at = ret;
        z.num_deleted_H_taut = i32::from(z.vt_group_info.tni.nNumRemovedExplicitH);
        let taut = &mut out_norm_data[TAUT_YES as usize];
        taut.num_at = z.num_atoms.wrapping_add(z.num_deleted_H_taut);
        taut.num_removed_H = z.num_deleted_H_taut;
        taut.nNumRemovedProtons = taut
            .nNumRemovedProtons
            .wrapping_add(i32::from(z.vt_group_info.tni.nNumRemovedProtons));
        for index in 0..NUM_H_ISOTOPES as usize {
            taut.nNumRemovedProtonsIsotopic[index] = taut.nNumRemovedProtonsIsotopic[index]
                .wrapping_add(z.vt_group_info.tni.nNumRemovedProtonsIsotopic[index]);
            taut.num_iso_H[index] =
                taut.num_iso_H[index].wrapping_add(z.vt_group_info.num_iso_H[index]);
        }
        if z.num_taut_at == 1
            && heap.slice(taut.at.as_const())?[0].at_type == ATT_PROTON as u16
            && z.vt_group_info.tni.nNumRemovedProtons == 1
        {
            taut.bDeleted = 1;
            FreeInpAtom(heap, Some(&mut taut.at_fixed_bonds))?;
        } else if z.vt_group_info.tni.bNormalizationFlags & u64::from(FLAG_NORM_CONSIDER_TAUT) != 0
            && !taut.at_fixed_bonds.is_null()
        {
            taut.bTautPreprocessed = 1;
        }
        taut.bTautFlags = z.vt_group_info.bTautFlags;
        *taut_flags = z.vt_group_info.bTautFlags;
        taut.bTautFlagsDone = z.vt_group_info.bTautFlagsDone;
        *taut_flags_done = z.vt_group_info.bTautFlagsDone;
        taut.bNormalizationFlags = z.vt_group_info.tni.bNormalizationFlags;
        inp2spATOM(
            heap,
            taut.at.as_const(),
            num_inp_at,
            z.at[TAUT_YES as usize],
        )?;
        ret = set_stereo_parity(
            pCG,
            heap,
            taut.at,
            z.at[TAUT_YES as usize],
            z.num_taut_at,
            z.num_deleted_H_taut,
            Some(&mut z.s[TAUT_YES as usize].nMaxNumStereoAtoms),
            Some(&mut z.s[TAUT_YES as usize].nMaxNumStereoBonds),
            z.nUserMode,
            z.bPointedEdgeStereo,
            z.vABParityUnknown,
            z.bLooseTSACheck,
            z.bStereoAtZz,
        )?;
        if ret < 0 {
            return Ok(ret);
        }
        z.s[TAUT_YES as usize].bMayHaveStereo = i32::from(
            z.s[TAUT_YES as usize].nMaxNumStereoAtoms != 0
                || z.s[TAUT_YES as usize].nMaxNumStereoBonds != 0,
        );
        z.s[TAUT_YES as usize].num_isotopic_atoms = set_atom_iso_sort_keys(
            heap,
            z.num_taut_at,
            z.at[TAUT_YES as usize],
            Some(&z.vt_group_info),
            Some(&mut z.s[TAUT_YES as usize].bHasIsotopicTautGroups),
        )?;
        z.s[TAUT_YES as usize].nLenLinearCTTautomer = CountTautomerGroups(
            heap,
            z.at[TAUT_YES as usize],
            z.num_taut_at,
            Some(&mut z.vt_group_info),
        )?;
        if z.s[TAUT_YES as usize].nLenLinearCTTautomer < 0 {
            return Ok(z.s[TAUT_YES as usize].nLenLinearCTTautomer);
        } else if z.s[TAUT_YES as usize].nLenLinearCTTautomer > 0 {
            z.num_at_tg = z.num_taut_at.wrapping_add(z.vt_group_info.num_t_groups);
            let _ = make_a_copy_of_t_group_info(
                heap,
                Some(&mut z.vt_group_info_orig),
                Some(&mut z.vt_group_info),
            )?;
            z.s[TAUT_YES as usize].nLenLinearCTIsotopicTautomer =
                set_tautomer_iso_sort_keys(heap, Some(&mut z.vt_group_info))?;
            if z.s[TAUT_YES as usize].nLenLinearCTIsotopicTautomer < 0 {
                z.s[TAUT_YES as usize].nLenLinearCTIsotopicTautomer = 0;
            }
            taut.bTautomeric = z.s[TAUT_YES as usize].nLenLinearCTTautomer;
        }
        let _ = GetCanonLengths(
            heap,
            z.num_taut_at,
            z.at[TAUT_YES as usize].as_const(),
            &mut z.s[TAUT_YES as usize],
            Some(&z.vt_group_info),
        )?;
    }

    if !out_norm_data[TAUT_NON as usize].at.is_null()
        && !out_norm_data[TAUT_YES as usize].at.is_null()
        && !z.at[TAUT_NON as usize].is_null()
        && z.s[TAUT_YES as usize].nLenLinearCTTautomer == 0
    {
        inchi_free(heap, z.at[TAUT_NON as usize])?;
        z.at[TAUT_NON as usize] = SourceMutPointer::null();
    } else if out_norm_data[TAUT_NON as usize].at.is_null()
        && !out_norm_data[TAUT_YES as usize].at.is_null()
        && z.at[TAUT_NON as usize].is_null()
        && !z.at[TAUT_YES as usize].is_null()
        && z.s[TAUT_YES as usize].nLenLinearCTTautomer == 0
    {
        out_norm_data[TAUT_YES as usize].bTautomeric = 0;
    } else if !out_norm_data[TAUT_NON as usize].at.is_null() && !z.at[TAUT_NON as usize].is_null() {
        let flags_pointer = heap.allocate_model_storage(vec![z.bTautFlags])?;
        let done_pointer = heap.allocate_model_storage(vec![z.bTautFlagsDone])?;
        let result = mark_alt_bonds_and_taut_groups(
            heap,
            ic,
            pCG,
            out_norm_data[TAUT_NON as usize].at,
            SourceMutPointer::null(),
            z.num_atoms,
            max_time,
            SourceMutPointer::null(),
            flags_pointer,
            done_pointer,
            0,
            SourceMutPointer::null(),
            clock_result,
        );
        z.bTautFlags = heap.slice(flags_pointer.as_const())?[0];
        z.bTautFlagsDone = heap.slice(done_pointer.as_const())?[0];
        heap.free(flags_pointer)?;
        heap.free(done_pointer)?;
        ret = result?;
        if ret < 0 {
            return Ok(ret);
        }
        let non = &mut out_norm_data[TAUT_NON as usize];
        non.num_at = z.num_atoms.wrapping_add(z.num_deleted_H);
        non.num_removed_H = z.num_deleted_H;
        non.bTautFlags = *taut_flags;
        non.bTautFlagsDone = *taut_flags_done;
        non.bNormalizationFlags = 0;
        inp2spATOM(heap, non.at.as_const(), num_inp_at, z.at[TAUT_NON as usize])?;
        ret = set_stereo_parity(
            pCG,
            heap,
            non.at,
            z.at[TAUT_NON as usize],
            z.num_atoms,
            z.num_deleted_H,
            Some(&mut z.s[TAUT_NON as usize].nMaxNumStereoAtoms),
            Some(&mut z.s[TAUT_NON as usize].nMaxNumStereoBonds),
            z.nUserMode,
            z.bPointedEdgeStereo,
            z.vABParityUnknown,
            z.bLooseTSACheck,
            z.bStereoAtZz,
        )?;
        if ret < 0 {
            return Ok(ret);
        }
        z.s[TAUT_NON as usize].bMayHaveStereo = i32::from(
            z.s[TAUT_NON as usize].nMaxNumStereoAtoms != 0
                || z.s[TAUT_NON as usize].nMaxNumStereoBonds != 0,
        );
        z.s[TAUT_NON as usize].num_isotopic_atoms =
            set_atom_iso_sort_keys(heap, z.num_atoms, z.at[TAUT_NON as usize], None, None)?;
        let _ = GetCanonLengths(
            heap,
            z.num_atoms,
            z.at[TAUT_NON as usize].as_const(),
            &mut z.s[TAUT_NON as usize],
            None,
        )?;
        non.bTautomeric = 0;
    }

    z.bMayHaveStereo = i32::from(
        z.s[TAUT_YES as usize].bMayHaveStereo != 0 || z.s[TAUT_NON as usize].bMayHaveStereo != 0,
    );
    z.bHasIsotopicAtoms = i32::from(
        z.s[TAUT_NON as usize].num_isotopic_atoms > 0
            || z.s[TAUT_NON as usize].bHasIsotopicTautGroups > 0
            || z.s[TAUT_YES as usize].num_isotopic_atoms > 0
            || z.s[TAUT_YES as usize].bHasIsotopicTautGroups > 0,
    );
    if z.fix_isofixedh != 0 {
        z.bHasIsotopicAtoms = i32::from(
            z.bHasIsotopicAtoms != 0
                || (z.s[TAUT_YES as usize].nLenLinearCTTautomer > 0
                    && z.vt_group_info
                        .tni
                        .nNumRemovedProtonsIsotopic
                        .iter()
                        .any(|value| *value != 0)),
        );
    }
    z.bHasIsotopicAtoms = i32::from(
        z.bHasIsotopicAtoms != 0
            || (z.s[TAUT_YES as usize].nLenIsotopicEndpoints > 1
                && z.vt_group_info.bTautFlagsDone
                    & u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)
                    != 0),
    );
    if z.nUserMode & u64::from(REQ_MODE_DEFAULT) == 0 {
        z.nUserMode |= u64::from(REQ_MODE_DEFAULT);
    }
    if z.nUserMode & u64::from(REQ_MODE_ISO) != 0 && z.bHasIsotopicAtoms == 0 {
        z.nUserMode ^= u64::from(REQ_MODE_ISO);
        z.nUserMode |= u64::from(REQ_MODE_NON_ISO);
    }
    if z.nUserMode & u64::from(REQ_MODE_STEREO) != 0 && z.nUserMode & u64::from(REQ_MODE_ISO) != 0 {
        z.nUserMode |= u64::from(REQ_MODE_ISO_STEREO);
    }
    if z.nUserMode & u64::from(REQ_MODE_STEREO) != 0
        && z.nUserMode & u64::from(REQ_MODE_NON_ISO) == 0
    {
        z.nUserMode ^= u64::from(REQ_MODE_STEREO);
    }
    if z.bMayHaveStereo == 0 {
        if z.nUserMode & u64::from(REQ_MODE_STEREO) != 0 {
            z.nUserMode ^= u64::from(REQ_MODE_STEREO);
        }
        if z.nUserMode & u64::from(REQ_MODE_ISO_STEREO) != 0 {
            z.nUserMode ^= u64::from(REQ_MODE_ISO_STEREO);
        }
    }
    if z.nUserMode & u64::from(REQ_MODE_BASIC) != 0
        && (out_norm_data[TAUT_NON as usize].at.is_null()
            || pp_inchi[TAUT_NON as usize].is_null()
            || pp_inchi_aux[TAUT_NON as usize].is_null()
            || z.at[TAUT_NON as usize].is_null())
    {
        z.nUserMode ^= u64::from(REQ_MODE_BASIC);
    }
    if z.nUserMode & u64::from(REQ_MODE_TAUT) != 0
        && (out_norm_data[TAUT_YES as usize].at.is_null()
            || pp_inchi[TAUT_YES as usize].is_null()
            || pp_inchi_aux[TAUT_YES as usize].is_null()
            || z.at[TAUT_YES as usize].is_null())
    {
        z.nUserMode ^= u64::from(REQ_MODE_TAUT);
    }
    match z.nUserMode & u64::from(REQ_MODE_BASIC | REQ_MODE_TAUT) {
        value if value == u64::from(REQ_MODE_BASIC) => {
            z.n1 = TAUT_NON as i32;
            z.n2 = TAUT_NON as i32;
        }
        value if value == u64::from(REQ_MODE_TAUT) => {
            z.n1 = TAUT_YES as i32;
            z.n2 = TAUT_YES as i32;
        }
        value if value == u64::from(REQ_MODE_BASIC | REQ_MODE_TAUT) => {
            z.n1 = TAUT_NON as i32;
            z.n2 = TAUT_YES as i32;
        }
        _ => return Ok(-3),
    }
    if ret == 0 {
        ret = z.num_atoms;
    }
    Ok(ret)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn NormOneStructureINChI(
    heap: &mut SourceHeap,
    canonical_globals: &mut CANON_GLOBALS,
    clock: SourceMutPointer<INCHI_CLOCK>,
    generation_data: &mut INCHIGEN_DATA,
    control: &mut INCHIGEN_CONTROL,
    inchi_kind: i32,
    mut input_file: Option<&mut crate::source_types::INCHI_IOSTREAM>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:188 NormOneStructureINChI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int NormOneStructureINChI( CANON_GLOBALS *pCG,
                           INCHI_CLOCK *ic,
                           INCHIGEN_DATA *gendata,
                           INCHIGEN_CONTROL *genctl,
                           int iINChI,
                           INCHI_IOSTREAM *inp_file )
{
    int k, i, j, nRet = 0;


    STRUCT_DATA *sd = &( genctl->StructData );

    INPUT_PARMS *ip = &( genctl->InpParms );
    ORIG_ATOM_DATA *prep_inp_data = &( genctl->PrepInpData[0] );
    ORIG_ATOM_DATA *orig_inp_data = &( genctl->OrigInpData );
    INCHI_IOSTREAM *out_file = genctl->inchi_file, *log_file = genctl->inchi_file + 1;
    INCHI_IOSTREAM prbstr, *prb_file = &prbstr;

    PINChI2 **pINChI2 = genctl->pINChI;
    PINChI_Aux2 **pINChI_Aux2 = genctl->pINChI_Aux;
    NORM_CANON_FLAGS *pncFlags = &( genctl->ncFlags );


    INP_ATOM_DATA *inp_cur_data = NULL;

    long num_inp = genctl->num_inp;

    INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
    ORIG_ATOM_DATA *cur_prep_inp_data = prep_inp_data + iINChI;

    inchiTime      ulTStart;

    /* To save intermediate data... */

    COMP_ATOM_DATA *composite_norm_data = genctl->composite_norm_data[iINChI];
    INP_ATOM_DATA2 *all_inp_norm_data = NULL;
    memset( composite_norm_data + TAUT_NON, 0, sizeof( composite_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( composite_norm_data + TAUT_YES, 0, sizeof( composite_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( composite_norm_data + TAUT_INI, 0, sizeof( composite_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    inchi_ios_init( prb_file, INCHI_IOS_TYPE_FILE, NULL );

    /*
        if ( orig_inp_data is NOT empty AND
             prep_inp_data[0] IS empty ) then:

            1. copy orig_inp_data --> prep_inp_data[0]
            2. fix odd things in prep_inp_data[0]
            3. if( orig_inp_data->bDisconnectSalts ) then
                  -- disconnect salts in prep_inp_data[0]
            4. move protons to neutralize charges on heteroatoms
            5. if( orig_inp_data->bDisconnectCoord ) then
                  -- copy prep_inp_data[0] --> prep_inp_data[1]
                  -- disconnect metals in prep_inp_data[0]

            [ This all is done in PreprocessOneStructure() ]

        iINChI = 0
        =========
        (normal/disconnected layer)

            1. normalize prep_inp_data[0] in inp_norm_data[0,1]
            2. create INChI[ iINChI ] out of inp_norm_data[0,1]


        iINChI = 1 AND orig_inp_data->bDisconnectCoord > 0
        =================================================
        (reconnected layer)

            1. normalize prep_inp_data[1] in inp_norm_data[0,1]
            2. create INChI[ iINChI ] out of inp_norm_data[0,1]

    */


    ip->msec_LeftTime = ip->msec_MaxTime; /* start timeout countdown for each component */

    if (ip->bAllowEmptyStructure && !orig_inp_data->at && !orig_inp_data->num_inp_atoms)
    {
        ;
    }
    else
    {
        if (!orig_inp_data->at || !orig_inp_data->num_inp_atoms)
        {
            return 0; /* nothing to do */
        }
    }
    if (iINChI == 1 && orig_inp_data->bDisconnectCoord <= 0)
    {
        return 0;
    }

   /* m = iINChI; */ /* orig_inp_data index */

    if (iINChI != INCHI_BAS && iINChI != INCHI_REC)
    {
        AddErrorMessage( sd->pStrErrStruct, "Fatal undetermined program error" );
        sd->nStructReadError = 97;
        nRet = sd->nErrorType = _IS_FATAL;
        goto exit_function;
    }

    /*******************************************************************
     *                                                                 *
     *                                                                 *
     *  Whole structure preprocessing: 1st step of the normalization   *
     *                                                                 *
     *                                                                 *
     *                                                                 *
     *******************************************************************/

    if (( !prep_inp_data->at || !prep_inp_data->num_inp_atoms ) && orig_inp_data->num_inp_atoms > 0)
    {
        /* the structure has not been preprocessed */
        if (ip->msec_MaxTime)
        {
            InchiTimeGet( &ulTStart );
        }

        PreprocessOneStructure( ic, sd, ip, orig_inp_data, prep_inp_data );
        pncFlags->bTautFlags[iINChI][TAUT_YES] =
            pncFlags->bTautFlags[iINChI][TAUT_NON] = sd->bTautFlags[INCHI_BAS] | ip->bTautFlags;
        pncFlags->bTautFlagsDone[iINChI][TAUT_YES] =
            pncFlags->bTautFlagsDone[iINChI][TAUT_NON] = sd->bTautFlagsDone[INCHI_BAS] | ip->bTautFlagsDone;

        switch (sd->nErrorType)
        {
            case _IS_ERROR:
            case _IS_FATAL:
                /* error message */
                nRet = TreatErrorsInReadTheStructure( sd,
                                                      ip,
                                                      LOG_MASK_ALL,
                                                      inp_file,
                                                      log_file,
                                                      out_file,
                                                      prb_file,
                                                      prep_inp_data,
                                                      &num_inp );
                goto exit_function;
        }
    }

    /* To save intermediate data... */
    if (prep_inp_data[iINChI].num_components > 1)
    {
        all_inp_norm_data = (INP_ATOM_DATA2 *) inchi_calloc( prep_inp_data[iINChI].num_components, sizeof( all_inp_norm_data[0] ) );
    }


    /* allocate pINChI[iINChI] and pINChI_Aux2[iINChI] -- arrays of pointers to INChI and INChI_Aux */
    /* assign values to sd->num_components[]                                                  */

    /* djb-rwth: MYREALLOC2( PINChI2, PINChI_Aux2, pINChI2[iINChI], pINChI_Aux2[iINChI], sd->num_components[iINChI], cur_prep_inp_data->num_components, k ) has been replaced and the whole block rewritten to address memory leaks and reading from freed memory locations */
    do 
    { 
        if( (sd->num_components[iINChI]) <= (cur_prep_inp_data->num_components) ) 
        {
            PINChI2* newPTR1 = (PINChI2 *)inchi_calloc( (long long)cur_prep_inp_data->num_components + 1, sizeof(PINChI2) );
            PINChI_Aux2* newPTR2 = (PINChI_Aux2*)inchi_calloc( (long long)cur_prep_inp_data->num_components + 1, sizeof(PINChI_Aux2) );
            if ( newPTR1 && newPTR2 )
            { 
                if (pINChI2[iINChI] && sd->num_components[iINChI] > 0)
                    memcpy( newPTR1, pINChI2[iINChI], (sd->num_components[iINChI]) * sizeof(PINChI2) );
                if (pINChI_Aux2[iINChI] && sd->num_components[iINChI] > 0)
                    memcpy( newPTR2, pINChI_Aux2[iINChI], (sd->num_components[iINChI]) * sizeof(PINChI_Aux2) );
                if (pINChI2[iINChI]) 
                    inchi_free(pINChI2[iINChI]);
                if (pINChI_Aux2[iINChI])
                    inchi_free(pINChI_Aux2[iINChI]);
                pINChI2[iINChI] = newPTR1;
                pINChI_Aux2[iINChI] = newPTR2;
                sd->num_components[iINChI] = cur_prep_inp_data->num_components;
                k  = 0;
            } 
            else 
            {        
                inchi_free(newPTR1); 
                inchi_free(newPTR2); 
                k = 1;
            }             
        } 
        else 
        { 
            k = 0; 
        }
    } while (0);

    if (k)
    {
        AddErrorMessage( sd->pStrErrStruct, "Cannot allocate output data. Terminating" );
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_FATAL;
        inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
        goto exit_function;
    }

    /* Allocate */

    /* visible */
    gendata->NormAtomsNontaut[iINChI] = (NORM_ATOMS *) inchi_calloc( sd->num_components[iINChI], sizeof( NORM_ATOMS ) );
    gendata->NormAtomsTaut[iINChI] = (NORM_ATOMS *) inchi_calloc( sd->num_components[iINChI], sizeof( NORM_ATOMS ) );
    /* invisible */
    genctl->InpNormAtData[iINChI] = (INP_ATOM_DATA *) inchi_calloc( sd->num_components[iINChI], sizeof( INP_ATOM_DATA ) );
    genctl->InpNormTautData[iINChI] = (INP_ATOM_DATA *) inchi_calloc( sd->num_components[iINChI], sizeof( INP_ATOM_DATA ) );
    genctl->InpCurAtData[iINChI] = (INP_ATOM_DATA *) inchi_calloc( sd->num_components[iINChI], sizeof( INP_ATOM_DATA ) );
    genctl->cti[iINChI] = (COMPONENT_TREAT_INFO *) inchi_calloc( sd->num_components[iINChI], sizeof( COMPONENT_TREAT_INFO ) );
    if (genctl->cti[iINChI]) /* djb-rwth: fixing a NULL pointer dereference */
        memset( genctl->cti[iINChI], 0, sd->num_components[iINChI] * sizeof( COMPONENT_TREAT_INFO ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* Second normalization step - component by component */

    for (i = 0, nRet = 0; !sd->bUserQuitComponent && i < cur_prep_inp_data->num_components; i++)
    {

        if (ip->msec_MaxTime)
        {
            InchiTimeGet(&ulTStart);
        }

        if (genctl->InpCurAtData[iINChI]) /* djb-rwth: fixing a NULL pointer dereference */
        {
            inp_cur_data = &(genctl->InpCurAtData[iINChI][i]);

            /*  a) allocate memory and extract current component */
            nRet = GetOneComponent(ic, sd, ip, log_file, out_file, inp_cur_data,
                cur_prep_inp_data, i, num_inp);
        }

        if (ip->msec_MaxTime)
        {
            ip->msec_LeftTime -= InchiTimeElapsed( ic, &ulTStart );
        }

        switch (nRet)
        {
            case _IS_ERROR:
            case _IS_FATAL:
                goto exit_cycle;
        }


        /*  c) Create the component's INChI ( copies ip->bTautFlags into sd->bTautFlags)*/

        inp_norm_data[TAUT_NON] = &( genctl->InpNormAtData[iINChI][i] );
        memset( inp_norm_data[TAUT_NON], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        inp_norm_data[TAUT_YES] = &( genctl->InpNormTautData[iINChI][i] );
        memset( inp_norm_data[TAUT_YES], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        nRet = NormOneComponentINChI( pCG, ic, genctl, iINChI, i );

        /* To save intermediate data... */
        if (all_inp_norm_data)
        {
            for (j = 0; j < TAUT_NUM; j++)
            {
                if (inp_norm_data[j]->bExists)
                {
                    all_inp_norm_data[i][j] = *inp_norm_data[j];
                    memset( inp_norm_data[j], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                }
            }
        }

        if (nRet)
        {
            nRet = TreatErrorsInCreateOneComponentINChI( sd, ip,
                                                         cur_prep_inp_data,
                                                         i, num_inp, inp_file,
                                                         log_file, out_file,
                                                         prb_file );
            break;
        }
    }

exit_cycle:

    switch (nRet)
    {
        case _IS_FATAL:
        case _IS_ERROR: break;
        default:

            /* To save intermediate data... */
            if (all_inp_norm_data)
            {
                CreateCompositeNormAtom( composite_norm_data,
                                         all_inp_norm_data,
                                         prep_inp_data[iINChI].num_components );
            }
            break;
    }

    /* When saving intermediate data - avoid memory leaks in case of error */
    if (all_inp_norm_data)
    {
        for (i = 0; i < prep_inp_data[iINChI].num_components; i++)
        {
            for (k = 0; k < TAUT_NUM; k++)
            {
                FreeInpAtomData( &all_inp_norm_data[i][k] );
            }
        }
        inchi_free( all_inp_norm_data );
        all_inp_norm_data = NULL;
    }

exit_function:

    return nRet;
}
    */
    // END INCHI C FUNCTION: NormOneStructureINChI
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: NormOneStructureINChI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; no function-local display branch is active.
    // INCHI✔️❌: SourceHeap map access and record cloning are materially slower than native pointer access.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: NormOneStructureINChI

    let allocation_is_null = |error: SourceHeapError| {
        matches!(
            error,
            SourceHeapError::AllocationFailed
                | SourceHeapError::AllocationSizeOverflow
                | SourceHeapError::AllocationElementCountOutOfRange
        )
    };
    let mut result = 0_i32;
    let mut problem_file = crate::source_types::INCHI_IOSTREAM::default();
    inchi_ios_init(
        Some(&mut problem_file),
        INCHI_IOS_TYPE_FILE as i32,
        SourceMutPointer::null(),
    )?;

    if let Ok(kind) = usize::try_from(inchi_kind) {
        if let Some(composite) = control.composite_norm_data.get_mut(kind) {
            *composite = std::array::from_fn(|_| COMP_ATOM_DATA::default());
        }
    }
    control.InpParms.msec_LeftTime = control.InpParms.msec_MaxTime;

    if control.InpParms.bAllowEmptyStructure != 0
        && control.OrigInpData.at.is_null()
        && control.OrigInpData.num_inp_atoms == 0
    {
    } else if control.OrigInpData.at.is_null() || control.OrigInpData.num_inp_atoms == 0 {
        return Ok(0);
    }
    if inchi_kind == INCHI_REC as i32 && control.OrigInpData.bDisconnectCoord <= 0 {
        return Ok(0);
    }
    let kind = match inchi_kind {
        value if value == INCHI_BAS as i32 => INCHI_BAS as usize,
        value if value == INCHI_REC as i32 => INCHI_REC as usize,
        _ => {
            let message = b"Fatal undetermined program error\0".map(|byte| byte as i8);
            AddErrorMessage(Some(&mut control.StructData.pStrErrStruct), Some(&message))?;
            control.StructData.nStructReadError = 97;
            control.StructData.nErrorType = _IS_FATAL as i32;
            return Ok(_IS_FATAL as i32);
        }
    };
    control.composite_norm_data[kind] =
        std::array::from_fn(|_| COMP_ATOM_DATA::default());

    let mut input_number = control.num_inp;
    if (control.PrepInpData[0].at.is_null()
        || control.PrepInpData[0].num_inp_atoms == 0)
        && control.OrigInpData.num_inp_atoms > 0
    {
        let mut start = inchiTime::default();
        if control.InpParms.msec_MaxTime != 0 {
            InchiTimeGet(&mut start, clock_result);
        }
        let mut clock_value = heap
            .slice(clock.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let _ = PreprocessOneStructure(
            heap,
            Some(&mut clock_value),
            &mut control.StructData,
            &control.InpParms,
            &mut control.OrigInpData,
            &mut control.PrepInpData,
        )?;
        heap.slice_mut(clock)?[0] = clock_value;
        let flags = control.StructData.bTautFlags[INCHI_BAS as usize]
            | control.InpParms.bTautFlags;
        control.ncFlags.bTautFlags[kind] = [flags; TAUT_NUM as usize];
        let flags_done = control.StructData.bTautFlagsDone[INCHI_BAS as usize]
            | control.InpParms.bTautFlagsDone;
        control.ncFlags.bTautFlagsDone[kind] = [flags_done; TAUT_NUM as usize];
        if matches!(
            control.StructData.nErrorType,
            value if value == _IS_ERROR as i32 || value == _IS_FATAL as i32
        ) {
            let (output_stream, remainder) = control.inchi_file.split_at_mut(1);
            let output_stream = &mut output_stream[0];
            let log_stream = &mut remainder[0];
            result = TreatErrorsInReadTheStructure(
                heap,
                &mut control.StructData,
                &control.InpParms,
                LOG_MASK_ALL as i32,
                input_file.as_deref_mut(),
                Some(log_stream),
                Some(output_stream),
                Some(&mut problem_file),
                &control.PrepInpData[0],
                &mut input_number,
            )?;
            return Ok(result);
        }
    }

    let component_count_i32 = control.PrepInpData[kind].num_components;
    let component_count = usize::try_from(component_count_i32.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut all_normalized = SourceMutPointer::<INP_ATOM_DATA2>::null();
    if component_count_i32 > 1 {
        all_normalized = match inchi_calloc::<INP_ATOM_DATA2>(
            heap,
            component_count_i32 as u64,
            SOURCE_SIZEOF_INP_ATOM_DATA2,
        ) {
            Ok(pointer) => pointer,
            Err(error) if allocation_is_null(error) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    }

    if control.StructData.num_components[kind] <= component_count_i32 {
        let row_count = i64::from(component_count_i32).wrapping_add(1) as u64;
        let new_inchi = match inchi_calloc::<crate::source_types::PINChI2>(
            heap,
            row_count,
            SOURCE_SIZEOF_PINCHI2,
        ) {
            Ok(pointer) => pointer,
            Err(error) if allocation_is_null(error) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        let new_aux = match inchi_calloc::<crate::source_types::PINChI_Aux2>(
            heap,
            row_count,
            SOURCE_SIZEOF_PINCHI2,
        ) {
            Ok(pointer) => pointer,
            Err(error) if allocation_is_null(error) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        if !new_inchi.is_null() && !new_aux.is_null() {
            let old_count = usize::try_from(control.StructData.num_components[kind].max(0))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            if !control.pINChI[kind].is_null() && old_count > 0 {
                let old = heap
                    .slice(control.pINChI[kind].as_const())?
                    .get(..old_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                heap.slice_mut(new_inchi)?[..old_count].copy_from_slice(&old);
            }
            if !control.pINChI_Aux[kind].is_null() && old_count > 0 {
                let old = heap
                    .slice(control.pINChI_Aux[kind].as_const())?
                    .get(..old_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                heap.slice_mut(new_aux)?[..old_count].copy_from_slice(&old);
            }
            inchi_free(heap, control.pINChI[kind])?;
            inchi_free(heap, control.pINChI_Aux[kind])?;
            control.pINChI[kind] = new_inchi;
            control.pINChI_Aux[kind] = new_aux;
            control.StructData.num_components[kind] = component_count_i32;
        } else {
            inchi_free(heap, new_inchi)?;
            inchi_free(heap, new_aux)?;
            let message = b"Cannot allocate output data. Terminating\0".map(|byte| byte as i8);
            AddErrorMessage(Some(&mut control.StructData.pStrErrStruct), Some(&message))?;
            control.StructData.nStructReadError = 99;
            control.StructData.nErrorType = _IS_FATAL as i32;
            if !all_normalized.is_null() {
                inchi_free(heap, all_normalized)?;
            }
            return Ok(result);
        }
    }

    let allocation_count = component_count_i32 as u64;
    generation_data.NormAtomsNontaut[kind] =
        match inchi_calloc::<NORM_ATOMS>(heap, allocation_count, SOURCE_SIZEOF_NORM_ATOMS) {
            Ok(pointer) => pointer,
            Err(error) if allocation_is_null(error) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    generation_data.NormAtomsTaut[kind] =
        match inchi_calloc::<NORM_ATOMS>(heap, allocation_count, SOURCE_SIZEOF_NORM_ATOMS) {
            Ok(pointer) => pointer,
            Err(error) if allocation_is_null(error) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    control.InpNormAtData[kind] =
        match inchi_calloc::<INP_ATOM_DATA>(heap, allocation_count, SOURCE_SIZEOF_NORM_ATOMS) {
            Ok(pointer) => pointer,
            Err(error) if allocation_is_null(error) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    control.InpNormTautData[kind] =
        match inchi_calloc::<INP_ATOM_DATA>(heap, allocation_count, SOURCE_SIZEOF_NORM_ATOMS) {
            Ok(pointer) => pointer,
            Err(error) if allocation_is_null(error) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    control.InpCurAtData[kind] =
        match inchi_calloc::<INP_ATOM_DATA>(heap, allocation_count, SOURCE_SIZEOF_NORM_ATOMS) {
            Ok(pointer) => pointer,
            Err(error) if allocation_is_null(error) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
    control.cti[kind] = match inchi_calloc::<COMPONENT_TREAT_INFO>(
        heap,
        allocation_count,
        SOURCE_SIZEOF_COMPONENT_TREAT_INFO,
    ) {
        Ok(pointer) => pointer,
        Err(error) if allocation_is_null(error) => SourceMutPointer::null(),
        Err(error) => return Err(error),
    };

    result = 0;
    for component in 0..component_count {
        if control.StructData.bUserQuitComponent != 0 {
            break;
        }
        let mut start = inchiTime::default();
        if control.InpParms.msec_MaxTime != 0 {
            InchiTimeGet(&mut start, clock_result);
        }
        if !control.InpCurAtData[kind].is_null() {
            let current_pointer = control.InpCurAtData[kind].offset(component as i64)?;
            let mut current_input = heap.slice(current_pointer.as_const())?[0].clone();
            let mut clock_value = heap.slice(clock.as_const())?[0].clone();
            let (output_stream, remainder) = control.inchi_file.split_at_mut(1);
            result = GetOneComponent(
                heap,
                &mut clock_value,
                &mut control.StructData,
                &control.InpParms,
                Some(&mut remainder[0]),
                Some(&mut output_stream[0]),
                &mut current_input,
                &control.PrepInpData[kind],
                component as i32,
                control.num_inp,
                clock_result,
                clock_result,
            )?;
            heap.slice_mut(clock)?[0] = clock_value;
            heap.slice_mut(current_pointer)?[0] = current_input;
        }
        if control.InpParms.msec_MaxTime != 0 {
            let elapsed = normalization_component_elapsed(heap, clock, &start, clock_result)?;
            control.InpParms.msec_LeftTime =
                control.InpParms.msec_LeftTime.wrapping_sub(elapsed);
        }
        if result == _IS_ERROR as i32 || result == _IS_FATAL as i32 {
            break;
        }

        let non_pointer = control.InpNormAtData[kind].offset(component as i64)?;
        let taut_pointer = control.InpNormTautData[kind].offset(component as i64)?;
        heap.slice_mut(non_pointer)?[0] = INP_ATOM_DATA::default();
        heap.slice_mut(taut_pointer)?[0] = INP_ATOM_DATA::default();
        result = NormOneComponentINChI(
            heap,
            canonical_globals,
            clock,
            control,
            inchi_kind,
            component as i32,
            clock_result,
        )?;

        if !all_normalized.is_null() {
            for (representation, pointer) in [non_pointer, taut_pointer].into_iter().enumerate() {
                if heap.slice(pointer.as_const())?[0].bExists != 0 {
                    let moved = std::mem::take(&mut heap.slice_mut(pointer)?[0]);
                    heap.slice_mut(all_normalized)?[component][representation] = moved;
                }
            }
        }
        if result != 0 {
            let (output_stream, remainder) = control.inchi_file.split_at_mut(1);
            result = TreatErrorsInCreateOneComponentINChI(
                heap,
                &mut control.StructData,
                &control.InpParms,
                &control.PrepInpData[kind],
                component as i32,
                control.num_inp,
                input_file.as_deref_mut(),
                Some(&mut remainder[0]),
                Some(&mut output_stream[0]),
                Some(&mut problem_file),
            )?;
            break;
        }
    }

    if result != _IS_FATAL as i32 && result != _IS_ERROR as i32 && !all_normalized.is_null() {
        let records = heap.slice(all_normalized.as_const())?.to_vec();
        let _ = CreateCompositeNormAtom(
            heap,
            &mut control.composite_norm_data[kind],
            &records,
            component_count_i32,
        )?;
    }
    if !all_normalized.is_null() {
        for component in 0..component_count {
            for representation in 0..TAUT_NUM as usize {
                let mut normalized = {
                    let records = heap.slice_mut(all_normalized)?;
                    std::mem::take(&mut records[component][representation])
                };
                FreeInpAtomData(heap, &mut normalized)?;
            }
        }
        inchi_free(heap, all_normalized)?;
    }
    Ok(result)
}

fn normalization_component_elapsed(
    heap: &mut SourceHeap,
    clock: SourceMutPointer<INCHI_CLOCK>,
    start: &inchiTime,
    clock_result: clock_t,
) -> Result<i64, SourceHeapError> {
    let mut clock_value = heap
        .slice(clock.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let elapsed = InchiTimeElapsed(&mut clock_value, Some(start), clock_result);
    *heap
        .slice_mut(clock)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = clock_value;
    Ok(elapsed)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn NormOneComponentINChI(
    heap: &mut SourceHeap,
    canonical_globals: &mut CANON_GLOBALS,
    clock: SourceMutPointer<INCHI_CLOCK>,
    control: &mut INCHIGEN_CONTROL,
    inchi_kind: i32,
    component_index: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:654 NormOneComponentINChI
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int NormOneComponentINChI( CANON_GLOBALS *pCG,
                           INCHI_CLOCK *ic,
                           INCHIGEN_CONTROL * genctl,
                           int iINChI,
                           int i )
{
    STRUCT_DATA *sd = &( genctl->StructData );
    INPUT_PARMS *ip = &( genctl->InpParms );
    PINChI2 **pINChI2 = genctl->pINChI;
    PINChI_Aux2 **pINChI_Aux2 = genctl->pINChI_Aux;
    NORM_CANON_FLAGS *pncFlags = &( genctl->ncFlags );

    inchiTime     ulTStart, ulTEnd, *pulTEnd = NULL;
    int           k, num_at, ret = 0;
    int           bOrigCoord;
    INCHI_MODE     bTautFlags = ip->bTautFlags;
    INCHI_MODE     bTautFlagsDone = ( ip->bTautFlagsDone | sd->bTautFlagsDone[INCHI_BAS] );
    long          lElapsedTime;
    /*
    PINChI2     *pINChI     = pINChI2[iINChI];
    PINChI_Aux2 *pINChI_Aux = pINChI_Aux2[iINChI];
    */

    PINChI2     *pINChI = NULL;
    PINChI_Aux2 *pINChI_Aux = NULL;

    INChI       *cur_INChI[TAUT_NUM];
    INChI_Aux   *cur_INChI_Aux[TAUT_NUM];


    /* pINChI2[m=iINChI-1][j< prep_inp_data[m].num_components][TAUT_NON] */

    INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
    INP_ATOM_DATA *inp_cur_data = NULL;

    COMPONENT_TREAT_INFO *cti = NULL;

    inp_cur_data = &( genctl->InpCurAtData[iINChI][i] );

    cti = &( genctl->cti[iINChI][i] );

    inp_norm_data[TAUT_NON] = &( genctl->InpNormAtData[iINChI][i] );
    inp_norm_data[TAUT_YES] = &( genctl->InpNormTautData[iINChI][i] );


    pINChI = pINChI2[iINChI];
    pINChI_Aux = pINChI_Aux2[iINChI];

    InchiTimeGet( &ulTStart );
    bOrigCoord = !( ip->bINChIOutputOptions & ( INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO ) );


    for (k = 0; k < TAUT_NUM; k++)
    {
        cur_INChI[k] = pINChI[i][k];
        cur_INChI_Aux[k] = pINChI_Aux[i][k];
    }

    /*  allocate memory for non-tautimeric (k=0) and tautomeric (k=1) results */
    for (k = 0; k < TAUT_NUM; k++)
    {
        int nAllocMode = ( k == TAUT_YES ? REQ_MODE_TAUT : 0 ) || /* djb-rwth: bitwise operator | changed with logical operator ||*/
            ( bTautFlagsDone & ( TG_FLAG_FOUND_ISOTOPIC_H_DONE |
                TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE ) ) ?
                ( ip->nMode & REQ_MODE_ISO ) : 0;

        if (inp_cur_data && ((k == TAUT_NON && ( ip->nMode & REQ_MODE_BASIC )) ||
             (k == TAUT_YES && ( ip->nMode & REQ_MODE_TAUT )))) /* djb-rwth: addressing LLVM warnings; fixing a NULL pointer dereference */
        {
            /*  alloc INChI and INChI_Aux */
            cur_INChI[k] = Alloc_INChI( inp_cur_data->at, inp_cur_data->num_at, &inp_cur_data->num_bonds,
                                          &inp_cur_data->num_isotopic, nAllocMode );
            cur_INChI_Aux[k] = Alloc_INChI_Aux( inp_cur_data->num_at,
                                          inp_cur_data->num_isotopic, nAllocMode, bOrigCoord );
            if (cur_INChI_Aux[k])
            {
                cur_INChI_Aux[k]->bIsIsotopic = inp_cur_data->num_isotopic;
            }
            /*  alloc memory for the output structure: non-tautomeric and tautomeric (for displaying) */
            CreateInpAtomData( inp_norm_data[k], inp_cur_data->num_at, k );
        }
        else
        {
            FreeInpAtomData( inp_norm_data[k] );
        }
    }

    lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
    if (ip->msec_MaxTime)
    {
        ip->msec_LeftTime -= lElapsedTime;
    }
    sd->ulStructTime += lElapsedTime;


    /******************************************************
     *
     *  Get one component canonical numberings, etc.
     *
     ******************************************************/

    /*
     * Create_INChI() return value:
     * num_at <= 0: error code
     * num_at >  0: number of atoms (excluding terminal hydrogen atoms)
     * inp_norm_data[0] => non-tautomeric, inp_norm_data[1] => tautomeric
     */
    InchiTimeGet( &ulTStart );
    if (ip->msec_MaxTime)
    {
        ulTEnd = ulTStart;
        pulTEnd = &ulTEnd;
        if (ip->msec_LeftTime > 0)
        {
            InchiTimeAddMsec( ic, pulTEnd, ip->msec_LeftTime );
        }
    }

    if (cti) /* djb-rwth: fixing a NULL pointer dereference */
    {
        cti->nUserMode = ip->nMode;
        cti->bLooseTSACheck = ip->bLooseTSACheck;
        cti->bStereoAtZz = ip->bStereoAtZz;

        /* vABParityUnknown holds actual value of an internal constant signifying       */
        /* unknown parity: either the same as for undefined parity (default==standard)  */
        /*  or a specific one (non-std; requested by SLUUD switch).                     */
        cti->vABParityUnknown = AB_PARITY_UNDF;
        if (0 != (ip->nMode & REQ_MODE_DIFF_UU_STEREO))
        {
            /* Make labels for unknown and undefined stereo different */
            cti->vABParityUnknown = AB_PARITY_UNKN;
        }
    }

    if (inp_cur_data) /* djb-rwth: fixing a NULL pointer dereference */
    {
        num_at = Normalization_step(pCG, ic,
            cur_INChI, cur_INChI_Aux,
            inp_cur_data->at, inp_norm_data, inp_cur_data->num_at,
            pulTEnd, &bTautFlags, &bTautFlagsDone, cti);

        SetConnectedComponentNumber(inp_cur_data->at, inp_cur_data->num_at, i + 1); /*  normalization alters structure component number */

        for (k = 0; k < TAUT_NUM; k++)
        {
            if (cur_INChI_Aux[k] && cur_INChI_Aux[k]->nNumberOfAtoms > 0)
            {
                pncFlags->bNormalizationFlags[iINChI][k] |= cur_INChI_Aux[k]->bNormalizationFlags;
                pncFlags->bTautFlags[iINChI][k] |= cur_INChI_Aux[k]->bTautFlags;
                pncFlags->bTautFlagsDone[iINChI][k] |= cur_INChI_Aux[k]->bTautFlagsDone;
                pncFlags->nCanonFlags[iINChI][k] |= cur_INChI_Aux[k]->nCanonFlags;
            }
        }

        /*  Detect errors */
        if (num_at < 0)
        {
            sd->nErrorCode = num_at;
        }
        else if (num_at == 0)
        {
            sd->nErrorCode = -1;
        }
        else if (cur_INChI[TAUT_NON] && cur_INChI[TAUT_NON]->nErrorCode)
        {
            /*  non-tautomeric error */
            sd->nErrorCode = cur_INChI[TAUT_NON]->nErrorCode;
        }
        else if (cur_INChI[TAUT_YES] && cur_INChI[TAUT_YES]->nErrorCode)
        {
            /*  tautomeric error */
            sd->nErrorCode = cur_INChI[TAUT_YES]->nErrorCode;
        }

        /*  detect and store stereo warnings */
        if (!sd->nErrorCode)
        {
            GetProcessingWarningsOneComponentInChI( cur_INChI,
                                                    inp_norm_data,
                                                    sd,
                                                    0 /*bNoWarnings */ );
        }

        lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
        if (ip->msec_MaxTime)
        {
            ip->msec_LeftTime -= lElapsedTime;
        }

        sd->ulStructTime += lElapsedTime;
    #ifndef TARGET_API_LIB
        /*  Display the results */
        if (ip->bDisplay)
            eat_keyboard_input();
    #endif
        /*  a) No matter what happened save the allocated INChI pointers */
        /*  save the INChI of the current component */

        InchiTimeGet(&ulTStart);
        for (k = 0; k < TAUT_NUM; k++)
        {
            pINChI[i][k] = cur_INChI[k];
            pINChI_Aux[i][k] = cur_INChI_Aux[k];

            cur_INChI[k] = NULL;
            cur_INChI_Aux[k] = NULL;
        }

        /*  b) Count one component structure and/or INChI results only if there was no error */
        /*     Set inp_norm_data[j]->num_removed_H = number of removed explicit H           */

        if (!sd->nErrorCode)
        {

            /*  find where the current processed structure is located */
            int cur_is_in_non_taut = (pINChI[i][TAUT_NON] && pINChI[i][TAUT_NON]->nNumberOfAtoms > 0);
            int cur_is_in_taut = (pINChI[i][TAUT_YES] && pINChI[i][TAUT_YES]->nNumberOfAtoms > 0);
            int cur_is_non_taut = (cur_is_in_non_taut && 0 == pINChI[i][TAUT_NON]->lenTautomer) ||
                (cur_is_in_taut && 0 == pINChI[i][TAUT_YES]->lenTautomer); /* djb-rwth: addressing LLVM warnings */
            int cur_is_taut = cur_is_in_taut && 0 < pINChI[i][TAUT_YES]->lenTautomer;

            if (cur_is_non_taut + cur_is_taut)
            {
                /*  count tautomeric and non-tautomeric components of the structures */
                int j1 = cur_is_in_non_taut ? TAUT_NON : TAUT_YES;
                int j2 = cur_is_in_taut ? TAUT_YES : TAUT_NON;
                int j;
                sd->num_non_taut[iINChI] += cur_is_non_taut;
                sd->num_taut[iINChI] += cur_is_taut;
                for (j = j1; j <= j2; j++)
                {
                    int bIsotopic = (pINChI[i][j]->nNumberOfIsotopicAtoms ||
                        pINChI[i][j]->nNumberOfIsotopicTGroups ||
                        (pINChI[i][j]->nPossibleLocationsOfIsotopicH && pINChI[i][j]->nPossibleLocationsOfIsotopicH[0] > 1)); /* djb-rwth: addressing LLVM warning */
                    if (j == TAUT_YES && pINChI_Aux[i][j]) /* djb-rwth: fixing a NULL pointer dereference */
                    {
                        bIsotopic |= (0 < pINChI_Aux[i][j]->nNumRemovedIsotopicH[0] +
                            pINChI_Aux[i][j]->nNumRemovedIsotopicH[1] +
                            pINChI_Aux[i][j]->nNumRemovedIsotopicH[2]);
                    }
                    inp_norm_data[j]->bExists = 1; /*  j=0: non-taut exists, j=1: taut exists */
                    inp_norm_data[j]->bHasIsotopicLayer = bIsotopic;
                }
            }
        }
    }

    if (sd->nErrorCode == CT_OUT_OF_RAM || sd->nErrorCode == CT_USER_QUIT_ERR)
    {
        ret = _IS_FATAL;
    }
    else if (sd->nErrorCode)
    {
        ret = _IS_ERROR;
    }

    lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
    if (ip->msec_MaxTime)
    {
        ip->msec_LeftTime -= lElapsedTime;
    }

    sd->ulStructTime += lElapsedTime;
    return ret;
}
    */
    // END INCHI C FUNCTION: NormOneComponentINChI
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: NormOneComponentINChI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; the !TARGET_API_LIB display branch is inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: NormOneComponentINChI

    let kind = usize::try_from(inchi_kind)
        .ok()
        .filter(|index| *index < control.InpCurAtData.len())
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let row = usize::try_from(component_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current_pointer = control.InpCurAtData[kind].offset(i64::from(component_index))?;
    let non_pointer = control.InpNormAtData[kind].offset(i64::from(component_index))?;
    let taut_pointer = control.InpNormTautData[kind].offset(i64::from(component_index))?;
    let cti_pointer = control.cti[kind].offset(i64::from(component_index))?;
    let inchi_row_pointer = control.pINChI[kind].offset(i64::from(component_index))?;
    let aux_row_pointer = control.pINChI_Aux[kind].offset(i64::from(component_index))?;

    let mut current_input = heap
        .slice(current_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut normalized_data = [
        heap.slice(non_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
        heap.slice(taut_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    ];
    let mut cti = heap
        .slice(cti_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut current_inchi = heap
        .slice(inchi_row_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut current_aux = heap
        .slice(aux_row_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut input_parameters = control.InpParms.clone();
    let mut structure_data = control.StructData.clone();
    let mut normalization_flags = control.ncFlags.clone();

    let mut start = inchiTime::default();
    InchiTimeGet(&mut start, clock_result);
    let original_coordinates = i32::from(
        input_parameters.bINChIOutputOptions
            & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO) as i32
            == 0,
    );
    let mut taut_flags = input_parameters.bTautFlags;
    let mut taut_flags_done =
        input_parameters.bTautFlagsDone | structure_data.bTautFlagsDone[INCHI_BAS as usize];
    let atom_count = usize::try_from(current_input.num_at.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let atoms = if current_input.at.is_null() {
        Vec::new()
    } else {
        heap.slice(current_input.at.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };

    for representation in 0..TAUT_NUM as usize {
        let allocation_condition = representation == TAUT_YES as usize
            || taut_flags_done
                & u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)
                != 0;
        let allocation_mode = if allocation_condition {
            (input_parameters.nMode & u64::from(REQ_MODE_ISO)) as i32
        } else {
            0
        };
        let enabled = (representation == TAUT_NON as usize
            && input_parameters.nMode & u64::from(REQ_MODE_BASIC) != 0)
            || (representation == TAUT_YES as usize
                && input_parameters.nMode & u64::from(REQ_MODE_TAUT) != 0);
        if enabled {
            current_inchi[representation] = Alloc_INChI(
                heap,
                &atoms,
                current_input.num_at,
                &mut current_input.num_bonds,
                &mut current_input.num_isotopic,
                allocation_mode,
            )?;
            current_aux[representation] = Alloc_INChI_Aux(
                heap,
                current_input.num_at,
                current_input.num_isotopic,
                allocation_mode,
                original_coordinates,
            )?;
            if !current_aux[representation].is_null() {
                heap.slice_mut(current_aux[representation])?[0].bIsIsotopic =
                    current_input.num_isotopic;
            }
            let _ = CreateInpAtomData(
                heap,
                &mut normalized_data[representation],
                current_input.num_at,
                representation as i32,
            )?;
        } else {
            FreeInpAtomData(heap, &mut normalized_data[representation])?;
        }
    }

    let mut elapsed = normalization_component_elapsed(heap, clock, &start, clock_result)?;
    if input_parameters.msec_MaxTime != 0 {
        input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
    }
    structure_data.ulStructTime = structure_data.ulStructTime.wrapping_add(elapsed as u64);

    InchiTimeGet(&mut start, clock_result);
    let maximum_time = if input_parameters.msec_MaxTime != 0 {
        let pointer = heap.allocate_model_storage(vec![start.clone()])?;
        if input_parameters.msec_LeftTime > 0 {
            let mut clock_value = heap.slice(clock.as_const())?[0].clone();
            InchiTimeAddMsec(
                &mut clock_value,
                heap.slice_mut(pointer)?.first_mut(),
                input_parameters.msec_LeftTime as u64,
            );
            heap.slice_mut(clock)?[0] = clock_value;
        }
        pointer
    } else {
        SourceMutPointer::null()
    };

    cti.nUserMode = input_parameters.nMode;
    cti.bLooseTSACheck = input_parameters.bLooseTSACheck;
    cti.bStereoAtZz = input_parameters.bStereoAtZz;
    cti.vABParityUnknown = AB_PARITY_UNDF as i32;
    if input_parameters.nMode & u64::from(REQ_MODE_DIFF_UU_STEREO) != 0 {
        cti.vABParityUnknown = AB_PARITY_UNKN as i32;
    }
    let normalization_result = Normalization_step(
        heap,
        canonical_globals,
        clock,
        &current_inchi,
        &current_aux,
        current_input.at,
        &mut normalized_data,
        current_input.num_at,
        maximum_time,
        &mut taut_flags,
        &mut taut_flags_done,
        &mut cti,
        clock_result,
    );
    if !maximum_time.is_null() {
        heap.free(maximum_time)?;
    }
    let num_atoms = normalization_result?;

    if !current_input.at.is_null() {
        SetConnectedComponentNumber(
            heap.slice_mut(current_input.at)?,
            current_input.num_at,
            component_index.wrapping_add(1),
        )?;
    }
    for representation in 0..TAUT_NUM as usize {
        if !current_aux[representation].is_null() {
            let aux = &heap.slice(current_aux[representation].as_const())?[0];
            if aux.nNumberOfAtoms > 0 {
                normalization_flags.bNormalizationFlags[kind][representation] |=
                    aux.bNormalizationFlags;
                normalization_flags.bTautFlags[kind][representation] |= aux.bTautFlags;
                normalization_flags.bTautFlagsDone[kind][representation] |= aux.bTautFlagsDone;
                normalization_flags.nCanonFlags[kind][representation] |= aux.nCanonFlags;
            }
        }
    }
    if num_atoms < 0 {
        structure_data.nErrorCode = num_atoms;
    } else if num_atoms == 0 {
        structure_data.nErrorCode = -1;
    } else if !current_inchi[TAUT_NON as usize].is_null()
        && heap.slice(current_inchi[TAUT_NON as usize].as_const())?[0].nErrorCode != 0
    {
        structure_data.nErrorCode =
            heap.slice(current_inchi[TAUT_NON as usize].as_const())?[0].nErrorCode;
    } else if !current_inchi[TAUT_YES as usize].is_null()
        && heap.slice(current_inchi[TAUT_YES as usize].as_const())?[0].nErrorCode != 0
    {
        structure_data.nErrorCode =
            heap.slice(current_inchi[TAUT_YES as usize].as_const())?[0].nErrorCode;
    }
    if structure_data.nErrorCode == 0 {
        let _ = GetProcessingWarningsOneComponentInChI(
            heap,
            &current_inchi,
            &normalized_data,
            &mut structure_data,
            0,
        )?;
    }
    elapsed = normalization_component_elapsed(heap, clock, &start, clock_result)?;
    if input_parameters.msec_MaxTime != 0 {
        input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
    }
    structure_data.ulStructTime = structure_data.ulStructTime.wrapping_add(elapsed as u64);

    InchiTimeGet(&mut start, clock_result);
    *heap
        .slice_mut(inchi_row_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = current_inchi;
    *heap
        .slice_mut(aux_row_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = current_aux;

    if structure_data.nErrorCode == 0 {
        let inchi_row = heap.slice(inchi_row_pointer.as_const())?[0].clone();
        let aux_row = heap.slice(aux_row_pointer.as_const())?[0].clone();
        let non_taut_in = !inchi_row[TAUT_NON as usize].is_null()
            && heap.slice(inchi_row[TAUT_NON as usize].as_const())?[0].nNumberOfAtoms > 0;
        let taut_in = !inchi_row[TAUT_YES as usize].is_null()
            && heap.slice(inchi_row[TAUT_YES as usize].as_const())?[0].nNumberOfAtoms > 0;
        let non_taut = (non_taut_in
            && heap.slice(inchi_row[TAUT_NON as usize].as_const())?[0].lenTautomer == 0)
            || (taut_in
                && heap.slice(inchi_row[TAUT_YES as usize].as_const())?[0].lenTautomer == 0);
        let taut =
            taut_in && heap.slice(inchi_row[TAUT_YES as usize].as_const())?[0].lenTautomer > 0;
        if non_taut || taut {
            structure_data.num_non_taut[kind] =
                structure_data.num_non_taut[kind].wrapping_add(i32::from(non_taut));
            structure_data.num_taut[kind] =
                structure_data.num_taut[kind].wrapping_add(i32::from(taut));
            let first = if non_taut_in { TAUT_NON } else { TAUT_YES } as usize;
            let last = if taut_in { TAUT_YES } else { TAUT_NON } as usize;
            for representation in first..=last {
                let inchi = &heap.slice(inchi_row[representation].as_const())?[0];
                let possible_isotope_h = !inchi.nPossibleLocationsOfIsotopicH.is_null()
                    && heap.slice(inchi.nPossibleLocationsOfIsotopicH.as_const())?[0] > 1;
                let mut isotopic = inchi.nNumberOfIsotopicAtoms != 0
                    || inchi.nNumberOfIsotopicTGroups != 0
                    || possible_isotope_h;
                if representation == TAUT_YES as usize && !aux_row[representation].is_null() {
                    let aux = &heap.slice(aux_row[representation].as_const())?[0];
                    isotopic |= aux
                        .nNumRemovedIsotopicH
                        .iter()
                        .fold(0_i32, |sum, value| sum.wrapping_add(i32::from(*value)))
                        > 0;
                }
                normalized_data[representation].bExists = 1;
                normalized_data[representation].bHasIsotopicLayer = i32::from(isotopic);
            }
        }
    }

    let ret = if structure_data.nErrorCode == CT_OUT_OF_RAM
        || structure_data.nErrorCode == CT_USER_QUIT_ERR
    {
        _IS_FATAL as i32
    } else if structure_data.nErrorCode != 0 {
        _IS_ERROR as i32
    } else {
        0
    };
    elapsed = normalization_component_elapsed(heap, clock, &start, clock_result)?;
    if input_parameters.msec_MaxTime != 0 {
        input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
    }
    structure_data.ulStructTime = structure_data.ulStructTime.wrapping_add(elapsed as u64);

    *heap
        .slice_mut(current_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = current_input;
    *heap
        .slice_mut(non_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = normalized_data[TAUT_NON as usize].clone();
    *heap
        .slice_mut(taut_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = normalized_data[TAUT_YES as usize].clone();
    *heap
        .slice_mut(cti_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = cti;
    control.InpParms = input_parameters;
    control.StructData = structure_data;
    control.ncFlags = normalization_flags;
    let _ = row;
    Ok(ret)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CanonOneStructureINChI(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    ic: SourceMutPointer<INCHI_CLOCK>,
    genctl: &mut INCHIGEN_CONTROL,
    iINChI: i32,
    mut inp_file: Option<&mut crate::source_types::INCHI_IOSTREAM>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:503 CanonOneStructureINChI
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap row access adds overhead.
    /*
int CanonOneStructureINChI( CANON_GLOBALS *pCG,
                            INCHI_CLOCK *ic,
                            INCHIGEN_CONTROL *genctl,
                            int iINChI,
                            INCHI_IOSTREAM *inp_file )
{
    int i, /*m,*/ nRet = 0;

    STRUCT_DATA *sd = &( genctl->StructData );

    INPUT_PARMS *ip = &( genctl->InpParms );
    INCHI_IOSTREAM *out_file = genctl->inchi_file, *log_file = genctl->inchi_file + 1;
    INCHI_IOSTREAM prbstr, *prb_file = &prbstr;

    ORIG_ATOM_DATA *prep_inp_data = &( genctl->PrepInpData[0] );

    long num_inp = genctl->num_inp;



    INP_ATOM_DATA *inp_cur_data = NULL;

    INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */

    ORIG_ATOM_DATA *cur_prep_inp_data = prep_inp_data + iINChI;

    inchiTime      ulTStart;

    inchi_ios_init( prb_file, INCHI_IOS_TYPE_FILE, NULL );

    for (i = 0; i < TAUT_NUM; i++) /* initialize in case no InChI to generate 2008-12-23 DT */
    {
        inp_norm_data[i] = NULL;
    }

    /**************************************************************************/
    /*                                                                        */
    /*                                                                        */
    /*   M A I N   C Y C L E:   P R O C E S S    C O M P O N E N T S          */
    /*                                                                        */
    /*                                                                        */
    /*                     O N E   B Y   O N E                                */
    /*                                                                        */
    /*                                                                        */
    /**************************************************************************/

    for (i = 0, nRet = 0; !sd->bUserQuitComponent && i < cur_prep_inp_data->num_components; i++)
    {
        if (ip->msec_MaxTime)
        {
            InchiTimeGet( &ulTStart );
        }


        /*****************************************************/
        /*  a) allocate memory and extract current component */
        /*****************************************************/

        inp_cur_data = &( genctl->InpCurAtData[iINChI][i] );


        nRet = GetOneComponent( ic, sd, ip, log_file, out_file,
                                inp_cur_data, cur_prep_inp_data,
                                i, num_inp );

        if (ip->msec_MaxTime)
        {
            ip->msec_LeftTime -= InchiTimeElapsed( ic, &ulTStart );
        }

        switch (nRet) { case _IS_ERROR: case _IS_FATAL: goto exit_cycle; }

#ifndef TARGET_API_LIB
        /*  console request: Display the component? */
        if (ip->bDisplay && inp_file != stdin)
        {
            if (user_quit( "Enter=Display Component, Esc=Stop ?", ip->ulDisplTime ))
            {
                sd->bUserQuitComponent = 1;
                break;
            }
        }
#endif




        /*******************************************************************************/
        /*                                                                             */
        /*                      C A N O N I C A L I Z A T I O N                        */
        /*                                                                             */
        /*         (both tautomeric and non-tautomeric if requested)                   */
        /*                                                                             */
        /*******************************************************************************/
        /*  c) Create the component's INChI ( copies ip->bTautFlags into sd->bTautFlags)*/
        /*******************************************************************************/

        inp_norm_data[TAUT_NON] = &( genctl->InpNormAtData[iINChI][i] );
        inp_norm_data[TAUT_YES] = &( genctl->InpNormTautData[iINChI][i] );

        nRet = CanonOneComponentINChI( pCG, ic, genctl, iINChI, i );




        if (nRet)
        {
            nRet = TreatErrorsInCreateOneComponentINChI( sd,
                                                         ip,
                                                         cur_prep_inp_data,
                                                         i,
                                                         num_inp,
                                                         inp_file,
                                                         log_file,
                                                         out_file,
                                                         prb_file );
            break;
        }
    }

    /**************************************************************************/
    /*                                                                        */
    /*                                                                        */
    /*   E N D   O F   T H E    M A I N   C Y C L E   P R O C E S S I N G     */
    /*                                                                        */
    /*          C O M P O N E N T S    O N E   B Y   O N E                    */
    /*                                                                        */
    /*                                                                        */
    /**************************************************************************/

exit_cycle:

    switch (nRet)
    {
        case _IS_FATAL:
        case _IS_ERROR: break;
        default:

            break;
    }

    for (i = 0; i < TAUT_NUM; i++)
    {
        FreeInpAtomData( inp_norm_data[i] );
    }

    return nRet;
}
    */
    // END INCHI C FUNCTION: CanonOneStructureINChI
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CanonOneStructureINChI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: The !TARGET_API_LIB console display/user_quit branch is inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CanonOneStructureINChI

    let kind = usize::try_from(iINChI)
        .ok()
        .filter(|index| *index < genctl.PrepInpData.len())
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut nRet = 0_i32;
    let mut problem_file = crate::source_types::INCHI_IOSTREAM::default();
    inchi_ios_init(
        Some(&mut problem_file),
        INCHI_IOS_TYPE_FILE as i32,
        SourceMutPointer::null(),
    )?;
    let component_count = usize::try_from(genctl.PrepInpData[kind].num_components.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let mut last_normalized = [
        SourceMutPointer::<INP_ATOM_DATA>::null();
        TAUT_NUM as usize
    ];

    for component in 0..component_count {
        if genctl.StructData.bUserQuitComponent != 0 {
            break;
        }
        let mut start = inchiTime::default();
        if genctl.InpParms.msec_MaxTime != 0 {
            InchiTimeGet(&mut start, clock_result);
        }

        let current_pointer =
            genctl.InpCurAtData[kind].offset(i64::try_from(component).unwrap_or(i64::MAX))?;
        let mut current_input = heap
            .slice(current_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut clock_value = heap
            .slice(ic.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        {
            let (out_file, remainder) = genctl.inchi_file.split_at_mut(1);
            nRet = GetOneComponent(
                heap,
                &mut clock_value,
                &mut genctl.StructData,
                &genctl.InpParms,
                Some(&mut remainder[0]),
                Some(&mut out_file[0]),
                &mut current_input,
                &genctl.PrepInpData[kind],
                component as i32,
                genctl.num_inp,
                clock_result,
                clock_result,
            )?;
        }
        *heap
            .slice_mut(ic)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = clock_value;
        *heap
            .slice_mut(current_pointer)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = current_input;

        if genctl.InpParms.msec_MaxTime != 0 {
            let elapsed = normalization_component_elapsed(heap, ic, &start, clock_result)?;
            genctl.InpParms.msec_LeftTime =
                genctl.InpParms.msec_LeftTime.wrapping_sub(elapsed);
        }
        if nRet == _IS_ERROR as i32 || nRet == _IS_FATAL as i32 {
            break;
        }

        last_normalized[TAUT_NON as usize] =
            genctl.InpNormAtData[kind].offset(component as i64)?;
        last_normalized[TAUT_YES as usize] =
            genctl.InpNormTautData[kind].offset(component as i64)?;
        nRet = CanonOneComponentINChI(
            heap,
            pCG,
            ic,
            genctl,
            iINChI,
            component as i32,
            clock_result,
        )?;

        if nRet != 0 {
            let input_parameters = genctl.InpParms.clone();
            let prepared = genctl.PrepInpData[kind].clone();
            let input_number = genctl.num_inp;
            let (out_file, remainder) = genctl.inchi_file.split_at_mut(1);
            nRet = TreatErrorsInCreateOneComponentINChI(
                heap,
                &mut genctl.StructData,
                &input_parameters,
                &prepared,
                component as i32,
                input_number,
                inp_file.as_deref_mut(),
                Some(&mut remainder[0]),
                Some(&mut out_file[0]),
                Some(&mut problem_file),
            )?;
            break;
        }
    }

    for pointer in last_normalized {
        if !pointer.is_null() {
            let mut normalized = std::mem::take(
                heap.slice_mut(pointer)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            FreeInpAtomData(heap, &mut normalized)?;
        }
    }
    Ok(nRet)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CanonOneComponentINChI(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    ic: SourceMutPointer<INCHI_CLOCK>,
    genctl: &mut INCHIGEN_CONTROL,
    iINChI: i32,
    i: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:923 CanonOneComponentINChI
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap row access and cloned control state add overhead.
    /*
int CanonOneComponentINChI( CANON_GLOBALS *pCG,
                           INCHI_CLOCK *ic,
                            INCHIGEN_CONTROL *genctl,
                            int iINChI,
                            int i )
{
    STRUCT_DATA *sd = &( genctl->StructData );
    INPUT_PARMS *ip = &( genctl->InpParms );
    PINChI2 **pINChI2 = genctl->pINChI;
    PINChI_Aux2 **pINChI_Aux2 = genctl->pINChI_Aux;
    NORM_CANON_FLAGS *pncFlags = &( genctl->ncFlags );

    inchiTime     ulTStart, ulTEnd, *pulTEnd = NULL;
    int           k, num_at, ret = 0;
    INChI       *cur_INChI[TAUT_NUM];
    INChI_Aux   *cur_INChI_Aux[TAUT_NUM];
    long          lElapsedTime;
    /*
    PINChI2     *pINChI     = pINChI2[iINChI];
    PINChI_Aux2 *pINChI_Aux = pINChI_Aux2[iINChI];
    */

    PINChI2     *pINChI = NULL;
    PINChI_Aux2 *pINChI_Aux = NULL;

    INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
    INP_ATOM_DATA *inp_cur_data = NULL;

    COMPONENT_TREAT_INFO *cti = NULL;

    inp_cur_data = &( genctl->InpCurAtData[iINChI][i] );

    cti = &( genctl->cti[iINChI][i] );

    inp_norm_data[TAUT_NON] = &( genctl->InpNormAtData[iINChI][i] );
    inp_norm_data[TAUT_YES] = &( genctl->InpNormTautData[iINChI][i] );

    pINChI = pINChI2[iINChI];
    pINChI_Aux = pINChI_Aux2[iINChI];

    InchiTimeGet( &ulTStart );

    for (k = 0; k < TAUT_NUM; k++)
    {
        cur_INChI[k] = pINChI[i][k];
        cur_INChI_Aux[k] = pINChI_Aux[i][k];
    }


    lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
    if (ip->msec_MaxTime)
    {
        ip->msec_LeftTime -= lElapsedTime;
    }
    sd->ulStructTime += lElapsedTime;


    /******************************************************
     *
     *  Get one component canonical numberings, etc.
     *
     ******************************************************/

    /*
     * Create_INChI() return value:
     * num_at <= 0: error code
     * num_at >  0: number of atoms (excluding terminal hydrogen atoms)
     * inp_norm_data[0] => non-tautomeric, inp_norm_data[1] => tautomeric
     */
    InchiTimeGet( &ulTStart );
    if (ip->msec_MaxTime)
    {
        ulTEnd = ulTStart;
        pulTEnd = &ulTEnd;
        if (ip->msec_LeftTime > 0)
        {
            InchiTimeAddMsec( ic, pulTEnd, ip->msec_LeftTime );
        }
    }
    num_at = Canonicalization_step( pCG, ic, cur_INChI, cur_INChI_Aux,
                                    inp_norm_data, pulTEnd, NULL, sd->pStrErrStruct,
                                    cti, ip->bLargeMolecules );

#ifndef FIX_DOCANON_RETCODE_RESET_BUG
/* The next line erroneously resets num_at which just became the return value */
/* of canon procedure. Thanks to David Foss of Accelrys for finding this. 2013-01-24 */
    num_at = cti->num_atoms;
#endif

    SetConnectedComponentNumber( inp_cur_data->at, inp_cur_data->num_at, i + 1 ); /*  normalization alters structure component number */
    for (k = 0; k < TAUT_NUM; k++)
    {
        if (cur_INChI_Aux[k] && cur_INChI_Aux[k]->nNumberOfAtoms > 0)
        {
            pncFlags->bNormalizationFlags[iINChI][k] |= cur_INChI_Aux[k]->bNormalizationFlags;
            pncFlags->bTautFlags[iINChI][k] |= cur_INChI_Aux[k]->bTautFlags;
            pncFlags->bTautFlagsDone[iINChI][k] |= cur_INChI_Aux[k]->bTautFlagsDone;
            pncFlags->nCanonFlags[iINChI][k] |= cur_INChI_Aux[k]->nCanonFlags;
        }
    }

    /*  Detect errors */
    if (num_at < 0)
    {
        sd->nErrorCode = num_at;
    }
    else if (num_at == 0)
    {
        sd->nErrorCode = -1;
    }
    else if (cur_INChI[TAUT_NON] && cur_INChI[TAUT_NON]->nErrorCode)
    {
        /*  non-tautomeric error */
        sd->nErrorCode = cur_INChI[TAUT_NON]->nErrorCode;
    }
    else if (cur_INChI[TAUT_YES] && cur_INChI[TAUT_YES]->nErrorCode)
    {
        /*  tautomeric error */
        sd->nErrorCode = cur_INChI[TAUT_YES]->nErrorCode;
    }

    /*  detect and store stereo warnings */
    if (!sd->nErrorCode)
    {
        GetProcessingWarningsOneComponentInChI( cur_INChI,
                                                inp_norm_data,
                                                sd,
                                                0 /*bNoWarnings */ );
    }

    lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
    if (ip->msec_MaxTime)
    {
        ip->msec_LeftTime -= lElapsedTime;
    }

    sd->ulStructTime += lElapsedTime;
#ifndef TARGET_API_LIB
    /*  Display the results */
    if (ip->bDisplay)
        eat_keyboard_input( );
#endif
    /*  a) No matter what happened save the allocated INChI pointers */
    /*  save the INChI of the current component */

    InchiTimeGet( &ulTStart );
    for (k = 0; k < TAUT_NUM; k++)
    {
        pINChI[i][k] = cur_INChI[k];
        pINChI_Aux[i][k] = cur_INChI_Aux[k];

        cur_INChI[k] = NULL;
        cur_INChI_Aux[k] = NULL;
    }

    /*  b) Count one component structure and/or INChI results only if there was no error */
    /*     Set inp_norm_data[j]->num_removed_H = number of removed explicit H           */

    if (!sd->nErrorCode)
    {

        /*  find where the current processed structure is located */
        int cur_is_in_non_taut = ( pINChI[i][TAUT_NON] && pINChI[i][TAUT_NON]->nNumberOfAtoms > 0 );
        int cur_is_in_taut = ( pINChI[i][TAUT_YES] && pINChI[i][TAUT_YES]->nNumberOfAtoms > 0 );
        int cur_is_non_taut = (cur_is_in_non_taut && 0 == pINChI[i][TAUT_NON]->lenTautomer) ||
            (cur_is_in_taut && 0 == pINChI[i][TAUT_YES]->lenTautomer); /* djb-rwth: addressing LLVM warnings */
        int cur_is_taut = cur_is_in_taut && 0 < pINChI[i][TAUT_YES]->lenTautomer;

        if (cur_is_non_taut + cur_is_taut)
        {
            /*  count tautomeric and non-tautomeric components of the structures */
            int j1 = cur_is_in_non_taut ? TAUT_NON : TAUT_YES;
            int j2 = cur_is_in_taut ? TAUT_YES : TAUT_NON;
            int j;
            sd->num_non_taut[iINChI] += cur_is_non_taut;
            sd->num_taut[iINChI] += cur_is_taut;
            for (j = j1; j <= j2; j++)
            {
                int bIsotopic = ( pINChI[i][j]->nNumberOfIsotopicAtoms ||
                                 pINChI[i][j]->nNumberOfIsotopicTGroups ||
                                 (pINChI[i][j]->nPossibleLocationsOfIsotopicH && pINChI[i][j]->nPossibleLocationsOfIsotopicH[0] > 1) ); /* djb-rwth: addressing LLVM warning */
                if (j == TAUT_YES && pINChI_Aux[i][j]) /* djb-rwth: fixing a NULL pointer dereference */
                {
                    bIsotopic |= ( 0 < pINChI_Aux[i][j]->nNumRemovedIsotopicH[0] +
                                      pINChI_Aux[i][j]->nNumRemovedIsotopicH[1] +
                                      pINChI_Aux[i][j]->nNumRemovedIsotopicH[2] );
                }
                inp_norm_data[j]->bExists = 1; /*  j=0: non-taut exists, j=1: taut exists */
                inp_norm_data[j]->bHasIsotopicLayer = bIsotopic;
            }
        }
    }

    if (sd->nErrorCode == CT_OUT_OF_RAM || sd->nErrorCode == CT_USER_QUIT_ERR)
    {
        ret = _IS_FATAL;
    }
    else if (sd->nErrorCode)
    {
        ret = _IS_ERROR;
    }

    lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
    if (ip->msec_MaxTime)
    {
        ip->msec_LeftTime -= lElapsedTime;
    }

    sd->ulStructTime += lElapsedTime;

    return ret;
}
    */
    // END INCHI C FUNCTION: CanonOneComponentINChI
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CanonOneComponentINChI
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: FIX_DOCANON_RETCODE_RESET_BUG == 1, so the obsolete num_at reset is inactive.
    // INCHI✔️❌: The !TARGET_API_LIB display/eat_keyboard_input branch is inactive.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CanonOneComponentINChI

    let kind = usize::try_from(iINChI)
        .ok()
        .filter(|index| *index < genctl.InpCurAtData.len())
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let component = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current_pointer = genctl.InpCurAtData[kind].offset(i64::from(i))?;
    let non_pointer = genctl.InpNormAtData[kind].offset(i64::from(i))?;
    let taut_pointer = genctl.InpNormTautData[kind].offset(i64::from(i))?;
    let cti_pointer = genctl.cti[kind].offset(i64::from(i))?;
    let inchi_row_pointer = genctl.pINChI[kind].offset(i64::from(i))?;
    let aux_row_pointer = genctl.pINChI_Aux[kind].offset(i64::from(i))?;

    let current_input = heap
        .slice(current_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut normalized_data = [
        heap.slice(non_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
        heap.slice(taut_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    ];
    let mut cti = heap
        .slice(cti_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let current_inchi = heap
        .slice(inchi_row_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let current_aux = heap
        .slice(aux_row_pointer.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut input_parameters = genctl.InpParms.clone();
    let mut structure_data = genctl.StructData.clone();
    let mut canonical_flags = genctl.ncFlags.clone();

    let mut start = inchiTime::default();
    InchiTimeGet(&mut start, clock_result);
    let mut elapsed = normalization_component_elapsed(heap, ic, &start, clock_result)?;
    if input_parameters.msec_MaxTime != 0 {
        input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
    }
    structure_data.ulStructTime = structure_data.ulStructTime.wrapping_add(elapsed as u64);

    InchiTimeGet(&mut start, clock_result);
    let maximum_time = if input_parameters.msec_MaxTime != 0 {
        let pointer = heap.allocate_model_storage(vec![start.clone()])?;
        if input_parameters.msec_LeftTime > 0 {
            let mut clock_value = heap
                .slice(ic.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            InchiTimeAddMsec(
                &mut clock_value,
                Some(
                    heap.slice_mut(pointer)?
                        .first_mut()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                ),
                input_parameters.msec_LeftTime as u64,
            );
            *heap
                .slice_mut(ic)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = clock_value;
        }
        pointer
    } else {
        SourceMutPointer::null()
    };

    let canonical_result = Canonicalization_step(
        heap,
        pCG,
        ic,
        &current_inchi,
        &current_aux,
        &normalized_data,
        maximum_time,
        None,
        Some(&mut structure_data.pStrErrStruct),
        &mut cti,
        input_parameters.bLargeMolecules,
        clock_result,
    );
    if !maximum_time.is_null() {
        heap.free(maximum_time)?;
    }
    let num_at = canonical_result?;

    if !current_input.at.is_null() {
        SetConnectedComponentNumber(
            heap.slice_mut(current_input.at)?,
            current_input.num_at,
            i.wrapping_add(1),
        )?;
    }
    for representation in 0..TAUT_NUM as usize {
        if !current_aux[representation].is_null() {
            let aux = heap
                .slice(current_aux[representation].as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if aux.nNumberOfAtoms > 0 {
                canonical_flags.bNormalizationFlags[kind][representation] |=
                    aux.bNormalizationFlags;
                canonical_flags.bTautFlags[kind][representation] |= aux.bTautFlags;
                canonical_flags.bTautFlagsDone[kind][representation] |= aux.bTautFlagsDone;
                canonical_flags.nCanonFlags[kind][representation] |= aux.nCanonFlags;
            }
        }
    }

    if num_at < 0 {
        structure_data.nErrorCode = num_at;
    } else if num_at == 0 {
        structure_data.nErrorCode = -1;
    } else if !current_inchi[TAUT_NON as usize].is_null()
        && heap
            .slice(current_inchi[TAUT_NON as usize].as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nErrorCode
            != 0
    {
        structure_data.nErrorCode = heap
            .slice(current_inchi[TAUT_NON as usize].as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nErrorCode;
    } else if !current_inchi[TAUT_YES as usize].is_null()
        && heap
            .slice(current_inchi[TAUT_YES as usize].as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nErrorCode
            != 0
    {
        structure_data.nErrorCode = heap
            .slice(current_inchi[TAUT_YES as usize].as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nErrorCode;
    }

    if structure_data.nErrorCode == 0 {
        let _ = GetProcessingWarningsOneComponentInChI(
            heap,
            &current_inchi,
            &normalized_data,
            &mut structure_data,
            0,
        )?;
    }

    elapsed = normalization_component_elapsed(heap, ic, &start, clock_result)?;
    if input_parameters.msec_MaxTime != 0 {
        input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
    }
    structure_data.ulStructTime = structure_data.ulStructTime.wrapping_add(elapsed as u64);

    InchiTimeGet(&mut start, clock_result);
    *heap
        .slice_mut(inchi_row_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = current_inchi;
    *heap
        .slice_mut(aux_row_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = current_aux;

    if structure_data.nErrorCode == 0 {
        let inchi_row = heap
            .slice(inchi_row_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let aux_row = heap
            .slice(aux_row_pointer.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let non_taut_in = !inchi_row[TAUT_NON as usize].is_null()
            && heap.slice(inchi_row[TAUT_NON as usize].as_const())?[0].nNumberOfAtoms > 0;
        let taut_in = !inchi_row[TAUT_YES as usize].is_null()
            && heap.slice(inchi_row[TAUT_YES as usize].as_const())?[0].nNumberOfAtoms > 0;
        let non_taut = (non_taut_in
            && heap.slice(inchi_row[TAUT_NON as usize].as_const())?[0].lenTautomer == 0)
            || (taut_in
                && heap.slice(inchi_row[TAUT_YES as usize].as_const())?[0].lenTautomer == 0);
        let taut =
            taut_in && heap.slice(inchi_row[TAUT_YES as usize].as_const())?[0].lenTautomer > 0;

        if non_taut || taut {
            structure_data.num_non_taut[kind] =
                structure_data.num_non_taut[kind].wrapping_add(i32::from(non_taut));
            structure_data.num_taut[kind] =
                structure_data.num_taut[kind].wrapping_add(i32::from(taut));
            let first = if non_taut_in { TAUT_NON } else { TAUT_YES } as usize;
            let last = if taut_in { TAUT_YES } else { TAUT_NON } as usize;
            for representation in first..=last {
                let inchi = &heap.slice(inchi_row[representation].as_const())?[0];
                let possible_isotopic_h = !inchi.nPossibleLocationsOfIsotopicH.is_null()
                    && heap.slice(inchi.nPossibleLocationsOfIsotopicH.as_const())?[0] > 1;
                let mut isotopic = inchi.nNumberOfIsotopicAtoms != 0
                    || inchi.nNumberOfIsotopicTGroups != 0
                    || possible_isotopic_h;
                if representation == TAUT_YES as usize && !aux_row[representation].is_null() {
                    let aux = &heap.slice(aux_row[representation].as_const())?[0];
                    isotopic |= aux
                        .nNumRemovedIsotopicH
                        .iter()
                        .fold(0_i32, |sum, value| sum.wrapping_add(i32::from(*value)))
                        > 0;
                }
                normalized_data[representation].bExists = 1;
                normalized_data[representation].bHasIsotopicLayer = i32::from(isotopic);
            }
        }
    }

    let ret = if structure_data.nErrorCode == CT_OUT_OF_RAM
        || structure_data.nErrorCode == CT_USER_QUIT_ERR
    {
        _IS_FATAL as i32
    } else if structure_data.nErrorCode != 0 {
        _IS_ERROR as i32
    } else {
        0
    };

    elapsed = normalization_component_elapsed(heap, ic, &start, clock_result)?;
    if input_parameters.msec_MaxTime != 0 {
        input_parameters.msec_LeftTime = input_parameters.msec_LeftTime.wrapping_sub(elapsed);
    }
    structure_data.ulStructTime = structure_data.ulStructTime.wrapping_add(elapsed as u64);

    *heap
        .slice_mut(non_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? =
        normalized_data[TAUT_NON as usize].clone();
    *heap
        .slice_mut(taut_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? =
        normalized_data[TAUT_YES as usize].clone();
    *heap
        .slice_mut(cti_pointer)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = cti;
    genctl.InpParms = input_parameters;
    genctl.StructData = structure_data;
    genctl.ncFlags = canonical_flags;
    let _ = component;
    Ok(ret)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Canonicalization_step(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    ic: SourceMutPointer<INCHI_CLOCK>,
    ppINChI: &[SourceMutPointer<INChI>; TAUT_NUM as usize],
    ppINChI_Aux: &[SourceMutPointer<INChI_Aux>; TAUT_NUM as usize],
    out_norm_data: &[INP_ATOM_DATA; TAUT_NUM as usize],
    ulMaxTime: SourceMutPointer<inchiTime>,
    mut ti_out: Option<&mut crate::source_types::T_GROUP_INFO>,
    mut pStrErrStruct: Option<&mut [i8]>,
    z: &mut COMPONENT_TREAT_INFO,
    LargeMolecules: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1612 Canonicalization_step
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access and stack-object pointer modeling add overhead.
    /*
int  Canonicalization_step( CANON_GLOBALS *pCG,
                            INCHI_CLOCK *ic,
                            INChI **ppINChI,
                            INChI_Aux **ppINChI_Aux,
                            INP_ATOM_DATA *out_norm_data[2],
                            struct tagInchiTime *ulMaxTime,
                            T_GROUP_INFO *ti_out,
                            char *pStrErrStruct,
                            COMPONENT_TREAT_INFO *z,
                            int LargeMolecules )
{
    int i, ret = 0, ret2 = 0;



    T_GROUP_INFO * /*const*/  t_group_info = &( z->vt_group_info );
    T_GROUP_INFO * /*const*/  t_group_info_orig = &( z->vt_group_info_orig );

    CANON_STAT  CS, CS2;
    CANON_STAT *pCS = &CS;
    CANON_STAT *pCS2 = &CS2;  /*  save all allocations to avoid memory leaks in case Canon_INChI() removes the pointer */

    BCN *pBCN = &( z->Bcn );


    INChI     *pINChI = NULL;         /* added initialization 2006-03 */
    INChI_Aux *pINChI_Aux = NULL;     /* added initialization 2006-03 */

    /************************************************************
     *                                                          *
     *       Obtain all non-stereo canonical numberings         *
     *                                                          *
     ************************************************************/

    if (( z->nUserMode & REQ_MODE_NON_ISO ) && !( z->nUserMode & REQ_MODE_ISO ))
    {

        /* added for special non-isotopic test mode 2004-10-04 */
        if (t_group_info)
        {
            t_group_info->bIgnoreIsotopic = 1;
            if (t_group_info->nIsotopicEndpointAtomNumber)
            {
                t_group_info->nIsotopicEndpointAtomNumber[0] = inchi_min( 1, t_group_info->nIsotopicEndpointAtomNumber[0] );
            }
            memset( t_group_info->num_iso_H, 0, sizeof( t_group_info->num_iso_H ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( t_group_info->tni.nNumRemovedProtonsIsotopic, 0, sizeof( t_group_info->tni.nNumRemovedProtonsIsotopic ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            t_group_info->bTautFlagsDone &= ~( TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE );
        }

        for (i = 0; i < TAUT_NUM; i++)
        {
            z->s[i].bHasIsotopicTautGroups = 0;
            z->s[i].bIgnoreIsotopic = 1;
            z->s[i].nLenIsotopic = 0;
            z->s[i].nLenIsotopicEndpoints = 0;
            z->s[i].nLenLinearCTIsotopicTautomer = 0;
            z->s[i].num_isotopic_atoms = 0;
        }
        z->bHasIsotopicAtoms = 0;
    }

    ret = GetBaseCanonRanking( ic,
                               z->num_atoms, z->num_at_tg, z->at,
                               t_group_info, z->s,
                               pBCN,
                               ulMaxTime,
                               pCG,
                               z->fix_isofixedh, LargeMolecules );

    if (ret < 0)
    {
        goto exit_function; /*  program error */
    }

    /* added for special non-isotopic test mode 2004-10-04 */
    if (!pBCN->ftcn[z->n1].PartitionCt.Rank)
    {
        z->n1 = ALT_TAUT( z->n1 );
    }

    if (!pBCN->ftcn[z->n2].PartitionCt.Rank)
    {
        z->n2 = ALT_TAUT( z->n2 );
    }

    if (z->n1 > z->n2)
    {
        ret = CT_TAUCOUNT_ERR;
        goto exit_function; /*  program error */
    }


    /************************************************************
     *                                                          *
     *       Obtain stereo canonical numberings                 *
     *                                                          *
     ************************************************************/

    for (i = z->n2; i >= z->n1 && !RETURNED_ERROR( ret ); i--)
    {

        memset( pCS, 0, sizeof( *pCS ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        switch (i)
        {

            case TAUT_NON: /*  non-tautomeric */

                z->nMode = 0;
                z->nMode = ( z->s[i].nLenLinearCTTautomer == 0 ) ? CANON_MODE_CT : CANON_MODE_TAUT;
                z->nMode |= ( z->bHasIsotopicAtoms && ( z->nUserMode & REQ_MODE_ISO ) ) ? CANON_MODE_ISO : 0;
                z->nMode |= ( z->s[TAUT_NON].bMayHaveStereo && ( z->nUserMode & REQ_MODE_STEREO ) ) ? CANON_MODE_STEREO : 0;
                z->nMode |= ( z->bHasIsotopicAtoms && z->s[TAUT_NON].bMayHaveStereo && ( z->nUserMode & REQ_MODE_ISO_STEREO ) ) ? CANON_MODE_ISO_STEREO : 0;
                z->nMode |= ( z->nUserMode & REQ_MODE_NOEQ_STEREO ) ? CMODE_NOEQ_STEREO : 0;
                z->nMode |= ( z->nUserMode & REQ_MODE_REDNDNT_STEREO ) ? CMODE_REDNDNT_STEREO : 0;
                z->nMode |= ( z->nUserMode & REQ_MODE_NO_ALT_SBONDS ) ? CMODE_NO_ALT_SBONDS : 0;

                /* 2010-01-12 */
                z->nMode |= ( z->vABParityUnknown == AB_PARITY_UNDF ) ? 0 : REQ_MODE_DIFF_UU_STEREO;

                if (( z->nMode & CANON_MODE_STEREO ) == CANON_MODE_STEREO ||
                    ( z->nMode & CANON_MODE_ISO_STEREO ) == CANON_MODE_ISO_STEREO)
                {
                    z->nMode |= ( z->nUserMode & REQ_MODE_RELATIVE_STEREO ) ? CMODE_RELATIVE_STEREO : 0;
                    z->nMode |= ( z->nUserMode & REQ_MODE_RACEMIC_STEREO ) ? CMODE_RACEMIC_STEREO : 0;
                    z->nMode |= ( z->nUserMode & REQ_MODE_SC_IGN_ALL_UU ) ? CMODE_SC_IGN_ALL_UU : 0;
                    z->nMode |= ( z->nUserMode & REQ_MODE_SB_IGN_ALL_UU ) ? CMODE_SB_IGN_ALL_UU : 0;
                }

                if ((ret = AllocateCS( pCS, z->num_atoms, z->num_atoms, z->s[TAUT_NON].nLenCT, z->s[TAUT_NON].nLenCTAtOnly,
                    z->s[TAUT_NON].nLenLinearCTStereoDble, z->s[TAUT_NON].nMaxNumStereoBonds,
                    z->s[TAUT_NON].nLenLinearCTStereoCarb, z->s[TAUT_NON].nMaxNumStereoAtoms,
                    0, 0, z->s[TAUT_NON].nLenIsotopic, z->nMode, pBCN ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }

                *pCS2 = *pCS;
                break;


            case TAUT_YES: /*  tautomeric */

                z->nMode = 0;
                z->nMode = ( z->s[i].nLenLinearCTTautomer == 0 ) ? CANON_MODE_CT : CANON_MODE_TAUT;
                z->nMode |= ( z->bHasIsotopicAtoms && ( z->nUserMode & REQ_MODE_ISO ) ) ? CANON_MODE_ISO : 0;
                z->nMode |= ( z->s[TAUT_YES].bMayHaveStereo && ( z->nUserMode & REQ_MODE_STEREO ) ) ? CANON_MODE_STEREO : 0;
                z->nMode |= ( z->bHasIsotopicAtoms && z->s[TAUT_YES].bMayHaveStereo && ( z->nUserMode & REQ_MODE_ISO_STEREO ) ) ? CANON_MODE_ISO_STEREO : 0;
                z->nMode |= ( z->nUserMode & REQ_MODE_NOEQ_STEREO ) ? CMODE_NOEQ_STEREO : 0;
                z->nMode |= ( z->nUserMode & REQ_MODE_REDNDNT_STEREO ) ? CMODE_REDNDNT_STEREO : 0;
                z->nMode |= ( z->nUserMode & REQ_MODE_NO_ALT_SBONDS ) ? CMODE_NO_ALT_SBONDS : 0;

                /* 2010-01-12 */
                z->nMode |= ( z->vABParityUnknown == AB_PARITY_UNDF ) ? 0 : REQ_MODE_DIFF_UU_STEREO;

                if (( z->nMode & CANON_MODE_STEREO ) == CANON_MODE_STEREO ||
                    ( z->nMode & CANON_MODE_ISO_STEREO ) == CANON_MODE_ISO_STEREO)
                {
                    z->nMode |= ( z->nUserMode & REQ_MODE_RELATIVE_STEREO ) ? CMODE_RELATIVE_STEREO : 0;
                    z->nMode |= ( z->nUserMode & REQ_MODE_RACEMIC_STEREO ) ? CMODE_RACEMIC_STEREO : 0;
                    z->nMode |= ( z->nUserMode & REQ_MODE_SC_IGN_ALL_UU ) ? CMODE_SC_IGN_ALL_UU : 0;
                    z->nMode |= ( z->nUserMode & REQ_MODE_SB_IGN_ALL_UU ) ? CMODE_SB_IGN_ALL_UU : 0;
                }

                if ((ret = AllocateCS( pCS, z->num_atoms, z->num_at_tg, z->s[TAUT_YES].nLenCT, z->s[TAUT_YES].nLenCTAtOnly,
                    z->s[TAUT_YES].nLenLinearCTStereoDble, z->s[TAUT_YES].nMaxNumStereoBonds,
                    z->s[TAUT_YES].nLenLinearCTStereoCarb, z->s[TAUT_YES].nMaxNumStereoAtoms,
                    z->s[TAUT_YES].nLenLinearCTTautomer, z->s[TAUT_YES].nLenLinearCTIsotopicTautomer,
                    z->s[TAUT_YES].nLenIsotopic, z->nMode, pBCN ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }

                *pCS2 = *pCS;
                break;
        } /* switch () */


        /*  settings */
        pCS->lNumDecreasedCT = -1;
        pCS->bDoubleBondSquare = DOUBLE_BOND_NEIGH_LIST ? 2 : 0;  /*  2 => special mode */
        pCS->bIgnoreIsotopic = !( ( z->s[TAUT_NON].num_isotopic_atoms ||
            z->s[TAUT_YES].num_isotopic_atoms ||
            z->s[TAUT_YES].bHasIsotopicTautGroups ) ||
            ( z->nUserMode & REQ_MODE_NON_ISO ) ||
                                               !( z->nUserMode & REQ_MODE_ISO ) );

        if (( z->nUserMode & REQ_MODE_NON_ISO ) && !( z->nUserMode & REQ_MODE_ISO ))
        {
            pCS->bIgnoreIsotopic = 1; /* 10-04-2004 */
        }

        if (i == TAUT_YES)
        {
            /* tautomeric */
            pCS->t_group_info = t_group_info; /*  ??? make a copy or reuse ???  */
            pCS->t_group_info->bIgnoreIsotopic = !( z->s[TAUT_YES].bHasIsotopicTautGroups ||
                ( z->nUserMode & REQ_MODE_NON_ISO ) ||
                                                   !( z->nUserMode & REQ_MODE_ISO ) );
            if (( z->nUserMode & REQ_MODE_NON_ISO ) && !( z->nUserMode & REQ_MODE_ISO ))
            {
                pCS->t_group_info->bIgnoreIsotopic = 1; /* 10-04-2004 */
            }
        }
        pCS->ulTimeOutTime = pBCN->ulTimeOutTime;
        /*=========== Obsolete Mode Bits (bit 0 is Least Significant Bit) ===========
         *
         *  Mode      Bits       Description
         *   '0' c    0          Only one connection table canonicalization
         *   '1' C    1          Recalculate CT using fixed nSymmRank
         *   '2' i    1|2        Isotopic canonicalization (internal)
         *   '3' I    1|2|4      Isotopic canonicalization (output)
         *   '4' s    1|8        Stereo canonicalization
         *   '5' S    1|2|4|16   Stereo isotopic canonicalization
         *   '6' A    1|2|4|8|16 Output All
         */

        /***************************************
           The last canonicalization step
         ***************************************/
        if ((i >= 0) && (i < 2)) /* djb-rwth: fixing buffer overruns */
        {
            if (pBCN)
            {
                /* USE_CANON2 == 1 */
                pCS->NeighList = NULL;
                pCS->pBCN = pBCN;
                ret = Canon_INChI(ic, z->num_atoms,
                    i ? z->num_at_tg : z->num_atoms,
                    z->at[i], pCS, pCG, z->nMode, i);
            }
            else
            {
                /* old way */
                pCS->NeighList = CreateNeighList(z->num_atoms, i ? z->num_at_tg : z->num_atoms, z->at[i], pCS->bDoubleBondSquare, pCS->t_group_info);
                pCS->pBCN = NULL;
                ret = Canon_INChI(ic, z->num_atoms,
                    i ? z->num_at_tg : z->num_atoms,
                    z->at[i], pCS, pCG, z->nMode, i);
            }
        }

        pINChI = ppINChI[i];      /* pointers to already allocated still empty InChI */
        pINChI_Aux = ppINChI_Aux[i];

        if (ret <= 0)
        {
            /***************************************/
            /*  failure in Canon_INChI()            */
            /***************************************/
            pINChI->nErrorCode = ret;
            pINChI_Aux->nErrorCode = ret;
        }
        else
        {
            /***************************************/
            /*  success Canon_INChI()               */
            /*  save canonicalization results in   */
            /*  pINChI and pINChI_Aux                */
            /***************************************/
            pINChI->nErrorCode = 0;
            pINChI_Aux->nErrorCode = 0;
            pINChI->bDeleted = pINChI_Aux->bDeleted = out_norm_data[i]->bDeleted;
            pINChI_Aux->nCanonFlags = pCS->nCanonFlags;
            pINChI_Aux->bTautFlags = out_norm_data[i]->bTautFlags;
            pINChI_Aux->bTautFlagsDone = out_norm_data[i]->bTautFlagsDone;
            pINChI_Aux->bNormalizationFlags = out_norm_data[i]->bNormalizationFlags;

            /*  may return an error or a warning */
            ret = FillOutINChIReducedWarn( pINChI, pINChI_Aux,
                                           z->num_atoms,
                                           i ? z->num_at_tg : z->num_atoms,
                                           i ? z->num_deleted_H_taut : z->num_deleted_H,
                                           z->at[i],
                                           out_norm_data[i]->at,
                                           pCS,
                                           pCG,
                                           i,
                                           z->nUserMode,
                                           pStrErrStruct );

            if (RETURNED_ERROR( ret ))
            {
                /* failure in FillOutINChI() */
                pINChI->nErrorCode = ret;
                pINChI_Aux->nErrorCode = ret;
            }
            else
            {
                /****************************/
                /* success in FillOutINChI() */
                /****************************/

                /* mark non-tautomeric representation as having another, tautomeric representation */
                if (pINChI_Aux && z->s[TAUT_YES].nLenLinearCTTautomer)
                    pINChI_Aux->bIsTautomeric = z->s[TAUT_YES].nLenLinearCTTautomer;


                ret2 = CheckCanonNumberingCorrectness( z->num_atoms,
                                                       i ? z->num_at_tg : z->num_atoms,
                                                       z->at[i],
                                                       pCS,
                                                       pCG,
                                                       i,
                                                       pStrErrStruct );

                if (ret2 && pINChI_Aux && pINChI) /* djb-rwth: fixing a NULL pointer dereference */
                {
                    pINChI->nErrorCode = ret2;
                    pINChI_Aux->nErrorCode = ret2;
                    ret = ret2;
                }
            } /* success in FillOutINChI */
        } /* success Canon_INChI */

        FreeNeighList( pCS->NeighList );
        DeAllocateCS( pCS2 );

        pINChI = NULL;      /* avoid dangling pointers */
        pINChI_Aux = NULL;  /* avoid dangling pointers */
    } /* for ( i = z->n2; i >= z->n1 && !RETURNED_ERROR( ret ); i -- )  */

    if (ret == 0)
    {
        ret = z->num_atoms;
    }

exit_function:

    DeAllocBCN( pBCN );

    if (z->at[TAUT_YES])
    {
        inchi_free( z->at[TAUT_YES] );
        z->at[TAUT_YES] = NULL;
    }

    if (z->at[TAUT_NON])
    {
        inchi_free( z->at[TAUT_NON] );
        z->at[TAUT_NON] = NULL;
    }

    if (ti_out)
    {
        *ti_out = *t_group_info;
    }
    else
    {
        free_t_group_info( t_group_info );
        t_group_info = NULL;
    }

    free_t_group_info( t_group_info_orig );

    return ret;
}
    */
    // END INCHI C FUNCTION: Canonicalization_step
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: Canonicalization_step
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; USE_CANON2 is active.
    // INCHI✔️❌: DOUBLE_BOND_NEIGH_LIST == 0; pBCN is always the address of z->Bcn.
    // INCHI✔️❌: ALT_TAUT(X) selects 1-X for the active TAUT_NON/TAUT_YES indices.
    // INCHI✔️❌: RETURNED_ERROR(X) is CT_ERR_MIN <= X && X <= CT_ERR_MAX.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: Canonicalization_step

    let returned_error = |value: i32| (CT_ERR_MIN..=CT_ERR_MAX).contains(&value);
    let mut clock = heap
        .slice(ic.as_const())?
        .first()
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let mut cs = CANON_STAT::default();
    let mut cs_allocations = CANON_STAT::default();
    let mut bcn_storage = SourceMutPointer::<BCN>::null();
    let mut tgi_storage = SourceMutPointer::<crate::source_types::T_GROUP_INFO>::null();
    let mut ret = 0_i32;

    let computation = (|| -> Result<i32, SourceHeapError> {
        if z.nUserMode & u64::from(REQ_MODE_NON_ISO) != 0
            && z.nUserMode & u64::from(REQ_MODE_ISO) == 0
        {
            z.vt_group_info.bIgnoreIsotopic = 1;
            if !z.vt_group_info.nIsotopicEndpointAtomNumber.is_null() {
                let endpoint = heap
                    .slice_mut(z.vt_group_info.nIsotopicEndpointAtomNumber)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                *endpoint = (*endpoint).min(1);
            }
            z.vt_group_info.num_iso_H.fill(0);
            z.vt_group_info
                .tni
                .nNumRemovedProtonsIsotopic
                .fill(0);
            z.vt_group_info.bTautFlagsDone &=
                !u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE);
            for size in &mut z.s {
                size.bHasIsotopicTautGroups = 0;
                size.bIgnoreIsotopic = 1;
                size.nLenIsotopic = 0;
                size.nLenIsotopicEndpoints = 0;
                size.nLenLinearCTIsotopicTautomer = 0;
                size.num_isotopic_atoms = 0;
            }
            z.bHasIsotopicAtoms = 0;
        }

        ret = GetBaseCanonRanking(
            heap,
            &mut clock,
            z.num_atoms,
            z.num_at_tg,
            z.at,
            Some(&mut z.vt_group_info),
            &z.s,
            &mut z.Bcn,
            ulMaxTime,
            pCG,
            z.fix_isofixedh,
            LargeMolecules,
            None,
            None,
            clock_result,
        )?;
        if ret < 0 {
            return Ok(ret);
        }

        if z.Bcn.ftcn[z.n1 as usize].PartitionCt.Rank.is_null() {
            z.n1 = 1_i32.wrapping_sub(z.n1);
        }
        if z.Bcn.ftcn[z.n2 as usize].PartitionCt.Rank.is_null() {
            z.n2 = 1_i32.wrapping_sub(z.n2);
        }
        if z.n1 > z.n2 {
            return Ok(CT_TAUCOUNT_ERR);
        }

        bcn_storage = heap.allocate_model_storage(vec![z.Bcn.clone()])?;

        let mut i = z.n2;
        while i >= z.n1 && !returned_error(ret) {
            cs = CANON_STAT::default();
            let index = i as usize;
            let mut canon_mode = if z.s[index].nLenLinearCTTautomer == 0 {
                CANON_MODE_CT
            } else {
                CANON_MODE_TAUT
            };
            if z.bHasIsotopicAtoms != 0 && z.nUserMode & u64::from(REQ_MODE_ISO) != 0 {
                canon_mode |= CANON_MODE_ISO;
            }
            if z.s[index].bMayHaveStereo != 0
                && z.nUserMode & u64::from(REQ_MODE_STEREO) != 0
            {
                canon_mode |= CANON_MODE_STEREO;
            }
            if z.bHasIsotopicAtoms != 0
                && z.s[index].bMayHaveStereo != 0
                && z.nUserMode & u64::from(REQ_MODE_ISO_STEREO) != 0
            {
                canon_mode |= CANON_MODE_ISO_STEREO;
            }
            if z.nUserMode & u64::from(REQ_MODE_NOEQ_STEREO) != 0 {
                canon_mode |= CMODE_NOEQ_STEREO;
            }
            if z.nUserMode & u64::from(REQ_MODE_REDNDNT_STEREO) != 0 {
                canon_mode |= CMODE_REDNDNT_STEREO;
            }
            if z.nUserMode & u64::from(REQ_MODE_NO_ALT_SBONDS) != 0 {
                canon_mode |= CMODE_NO_ALT_SBONDS;
            }
            if z.vABParityUnknown != AB_PARITY_UNDF as i32 {
                canon_mode |= REQ_MODE_DIFF_UU_STEREO;
            }
            if canon_mode & CANON_MODE_STEREO == CANON_MODE_STEREO
                || canon_mode & CANON_MODE_ISO_STEREO == CANON_MODE_ISO_STEREO
            {
                if z.nUserMode & u64::from(REQ_MODE_RELATIVE_STEREO) != 0 {
                    canon_mode |= CMODE_RELATIVE_STEREO;
                }
                if z.nUserMode & u64::from(REQ_MODE_RACEMIC_STEREO) != 0 {
                    canon_mode |= CMODE_RACEMIC_STEREO;
                }
                if z.nUserMode & u64::from(REQ_MODE_SC_IGN_ALL_UU) != 0 {
                    canon_mode |= CMODE_SC_IGN_ALL_UU;
                }
                if z.nUserMode & u64::from(REQ_MODE_SB_IGN_ALL_UU) != 0 {
                    canon_mode |= CMODE_SB_IGN_ALL_UU;
                }
            }
            z.nMode = u64::from(canon_mode);

            ret = AllocateCS(
                heap,
                &mut cs,
                z.num_atoms,
                if i == TAUT_YES as i32 {
                    z.num_at_tg
                } else {
                    z.num_atoms
                },
                z.s[index].nLenCT,
                z.s[index].nLenCTAtOnly,
                z.s[index].nLenLinearCTStereoDble,
                z.s[index].nMaxNumStereoBonds,
                z.s[index].nLenLinearCTStereoCarb,
                z.s[index].nMaxNumStereoAtoms,
                if i == TAUT_YES as i32 {
                    z.s[index].nLenLinearCTTautomer
                } else {
                    0
                },
                if i == TAUT_YES as i32 {
                    z.s[index].nLenLinearCTIsotopicTautomer
                } else {
                    0
                },
                z.s[index].nLenIsotopic,
                z.nMode,
                bcn_storage,
            )?;
            if ret != 0 {
                return Ok(ret);
            }
            cs_allocations = cs.clone();

            cs.lNumDecreasedCT = -1;
            cs.bDoubleBondSquare = 0;
            cs.bIgnoreIsotopic = i32::from(
                !((z.s[TAUT_NON as usize].num_isotopic_atoms != 0
                    || z.s[TAUT_YES as usize].num_isotopic_atoms != 0
                    || z.s[TAUT_YES as usize].bHasIsotopicTautGroups != 0)
                    || z.nUserMode & u64::from(REQ_MODE_NON_ISO) != 0
                    || z.nUserMode & u64::from(REQ_MODE_ISO) == 0),
            );
            if z.nUserMode & u64::from(REQ_MODE_NON_ISO) != 0
                && z.nUserMode & u64::from(REQ_MODE_ISO) == 0
            {
                cs.bIgnoreIsotopic = 1;
            }

            if i == TAUT_YES as i32 {
                if tgi_storage.is_null() {
                    tgi_storage =
                        heap.allocate_model_storage(vec![z.vt_group_info.clone()])?;
                } else {
                    heap.slice_mut(tgi_storage)?[0] = z.vt_group_info.clone();
                }
                cs.t_group_info = tgi_storage;
                z.vt_group_info.bIgnoreIsotopic = i32::from(
                    !(z.s[TAUT_YES as usize].bHasIsotopicTautGroups != 0
                        || z.nUserMode & u64::from(REQ_MODE_NON_ISO) != 0
                        || z.nUserMode & u64::from(REQ_MODE_ISO) == 0),
                );
                if z.nUserMode & u64::from(REQ_MODE_NON_ISO) != 0
                    && z.nUserMode & u64::from(REQ_MODE_ISO) == 0
                {
                    z.vt_group_info.bIgnoreIsotopic = 1;
                }
                heap.slice_mut(tgi_storage)?[0] = z.vt_group_info.clone();
            }

            cs.ulTimeOutTime = heap.slice(bcn_storage.as_const())?[0].ulTimeOutTime;
            cs.NeighList = SourceMutPointer::null();
            cs.pBCN = bcn_storage;
            ret = Canon_INChI(
                heap,
                &mut clock,
                clock_result,
                None,
                None,
                z.num_atoms,
                if i == TAUT_YES as i32 {
                    z.num_at_tg
                } else {
                    z.num_atoms
                },
                z.at[index],
                &mut cs,
                pCG,
                canon_mode,
                i,
            )?;
            if i == TAUT_YES as i32 {
                z.vt_group_info = heap.slice(tgi_storage.as_const())?[0].clone();
            }

            let mut inchi = heap
                .slice(ppINChI[index].as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut aux = heap
                .slice(ppINChI_Aux[index].as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;

            if ret <= 0 {
                inchi.nErrorCode = ret;
                aux.nErrorCode = ret;
            } else {
                inchi.nErrorCode = 0;
                aux.nErrorCode = 0;
                inchi.bDeleted = out_norm_data[index].bDeleted;
                aux.bDeleted = out_norm_data[index].bDeleted;
                aux.nCanonFlags = cs.nCanonFlags;
                aux.bTautFlags = out_norm_data[index].bTautFlags;
                aux.bTautFlagsDone = out_norm_data[index].bTautFlagsDone;
                aux.bNormalizationFlags = out_norm_data[index].bNormalizationFlags;

                ret = FillOutINChIReducedWarn(
                    heap,
                    &mut inchi,
                    &mut aux,
                    z.num_atoms,
                    if i == TAUT_YES as i32 {
                        z.num_at_tg
                    } else {
                        z.num_atoms
                    },
                    if i == TAUT_YES as i32 {
                        z.num_deleted_H_taut
                    } else {
                        z.num_deleted_H
                    },
                    z.at[index],
                    out_norm_data[index].at,
                    &mut cs,
                    pCG,
                    i,
                    z.nUserMode,
                    pStrErrStruct.as_deref_mut(),
                )?;
                if returned_error(ret) {
                    inchi.nErrorCode = ret;
                    aux.nErrorCode = ret;
                } else {
                    if z.s[TAUT_YES as usize].nLenLinearCTTautomer != 0 {
                        aux.bIsTautomeric =
                            z.s[TAUT_YES as usize].nLenLinearCTTautomer;
                    }
                    let ret2 = CheckCanonNumberingCorrectness(
                        heap,
                        z.num_atoms,
                        if i == TAUT_YES as i32 {
                            z.num_at_tg
                        } else {
                            z.num_atoms
                        },
                        z.at[index].as_const(),
                        &mut cs,
                        pCG,
                        i,
                        pStrErrStruct.as_deref_mut(),
                    )?;
                    if ret2 != 0 {
                        inchi.nErrorCode = ret2;
                        aux.nErrorCode = ret2;
                        ret = ret2;
                    }
                }
            }

            *heap
                .slice_mut(ppINChI[index])?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = inchi;
            *heap
                .slice_mut(ppINChI_Aux[index])?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = aux;

            DeAllocateCS(heap, &mut cs_allocations)?;
            cs = CANON_STAT::default();
            cs_allocations = CANON_STAT::default();
            i = i.wrapping_sub(1);
        }

        if ret == 0 {
            ret = z.num_atoms;
        }
        Ok(ret)
    })();

    DeAllocateCS(heap, &mut cs_allocations)?;
    if !tgi_storage.is_null() {
        z.vt_group_info = heap.slice(tgi_storage.as_const())?[0].clone();
    }
    if !bcn_storage.is_null() {
        z.Bcn = heap.slice(bcn_storage.as_const())?[0].clone();
    }
    DeAllocBCN(heap, Some(&mut z.Bcn))?;

    for index in [TAUT_YES as usize, TAUT_NON as usize] {
        if !z.at[index].is_null() {
            inchi_free(heap, z.at[index])?;
            z.at[index] = SourceMutPointer::null();
        }
    }

    if let Some(output) = ti_out.as_deref_mut() {
        *output = z.vt_group_info.clone();
    } else {
        free_t_group_info(heap, Some(&mut z.vt_group_info))?;
    }
    free_t_group_info(heap, Some(&mut z.vt_group_info_orig))?;

    if !tgi_storage.is_null() {
        heap.free(tgi_storage)?;
    }
    if !bcn_storage.is_null() {
        heap.free(bcn_storage)?;
    }
    *heap
        .slice_mut(ic)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)? = clock;

    computation
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FillOutINChIReducedWarn(
    heap: &mut SourceHeap,
    pINChI: &mut INChI,
    pINChI_Aux: &mut INChI_Aux,
    num_atoms: i32,
    num_at_tg: i32,
    num_removed_H: i32,
    at: SourceMutPointer<sp_ATOM>,
    norm_at: SourceMutPointer<inp_ATOM>,
    pCS: &mut CANON_STAT,
    pCG: &mut CANON_GLOBALS,
    bTautomeric: i32,
    nUserMode: INCHI_MODE,
    pStrErrStruct: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2254 FillOutINChIReducedWarn
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
int FillOutINChIReducedWarn( INChI *pINChI,
                             INChI_Aux *pINChI_Aux,
                             int num_atoms,
                             int num_at_tg,
                             int num_removed_H,
                             sp_ATOM *at,
                             inp_ATOM *norm_at,
                             CANON_STAT *pCS,
                             CANON_GLOBALS *pCG,
                             int bTautomeric,
                             INCHI_MODE nUserMode,
                             char *pStrErrStruct )
{
    int i, j, m, n, g, len, ii, ret = 0;

    AT_NUMB   *pSymmRank, *pOrigNosInCanonOrd, *pConstitEquNumb, *pCanonOrd = NULL, *pCanonOrdInv = NULL, *pCanonOrdTaut;
    T_GROUP_INFO     *t_group_info = pCS->t_group_info;
    T_GROUP *t_group;
    int nErrorCode = 0;
    AT_NUMB *pCanonRank, *pCanonRankInv; /* canonical ranks of the atoms or tautomeric groups */
    AT_NUMB *pCanonRankAtoms = NULL, *pSortOrd = NULL;
    AT_RANK nMinOrd;
    INChI_Stereo *Stereo;
    int          bUseNumberingInv = 0, bUseIsotopicNumberingInv = 0;
    INCHI_MODE    nStereoUnmarkMode;

    /*AT_NUMB  *pCanonOrdNonIso = NULL, *pCanonOrdIso = NULL;*/
    /*AT_NUMB  *nOrigAtNosInCanonOrdNonIso = NULL, *nOrigAtNosInCanonOrdIso = NULL;*/

    /*  Check for warnings */
    if (pCS->nLenLinearCTStereoCarb < 0 || pCS->nLenLinearCTStereoDble < 0 ||
         pCS->nLenCanonOrdStereo < 0 || pCS->nLenCanonOrdStereoTaut < 0)
    {
        nErrorCode |= WARN_FAILED_STEREO;
    }
    if (pCS->nLenLinearCTIsotopic < 0 || pCS->nLenLinearCTIsotopicTautomer < 0 ||
         pCS->nLenCanonOrdIsotopic < 0 || pCS->nLenCanonOrdIsotopicTaut < 0)
    {
        nErrorCode |= WARN_FAILED_ISOTOPIC;
    }
    if (pCS->nLenLinearCTIsotopicStereoCarb < 0 || pCS->nLenLinearCTIsotopicStereoDble < 0 ||
         pCS->nLenCanonOrdIsotopicStereo < 0 || pCS->nLenCanonOrdIsotopicStereoTaut < 0)
    {
        nErrorCode |= WARN_FAILED_ISOTOPIC_STEREO;
    }
    pCanonRankAtoms = (AT_NUMB *) inchi_calloc( (long long)num_at_tg + 1, sizeof( pCanonRankAtoms[0] ) ); /* djb-rwth: cast operator added */
    pSortOrd = (AT_NUMB *) inchi_calloc( (long long)num_at_tg + 1, sizeof( pSortOrd[0] ) ); /*  must have more than num_atoms */ /* djb-rwth: cast operator added */

    if (!pCanonRankAtoms || !pSortOrd)
    {
        nErrorCode = 0;
        ret = CT_OUT_OF_RAM;  /*   <BRKPT> */
        pINChI->nErrorCode = pINChI_Aux->nErrorCode = CT_OUT_OF_RAM;
        goto exit_function;
    }

    /*  total charge */
    for (i = 0, n = 0; i < num_atoms + num_removed_H; i++)
    {
        n += at[i].charge;
    }
    pINChI->nTotalCharge = n;

    /*  number of atoms */
    pINChI->nNumberOfAtoms = num_atoms;
    pINChI_Aux->nNumberOfAtoms = num_atoms;

    /* removed protons and detachable isotopic H */
    if (bTautomeric && t_group_info)
    {
        pINChI_Aux->nNumRemovedProtons = t_group_info->tni.nNumRemovedProtons;
        for (i = 0; i < NUM_H_ISOTOPES; i++)
        {
            pINChI_Aux->nNumRemovedIsotopicH[i] = t_group_info->num_iso_H[i]
                + t_group_info->tni.nNumRemovedProtonsIsotopic[i];
        }
        if (pINChI_Aux->bNormalizationFlags & FLAG_FORCE_SALT_TAUT)
        {
            pINChI->nFlags |= INCHI_FLAG_HARD_ADD_REM_PROTON;
        }
        /*^^^
        if ( pINChI_Aux->bNormalizationFlags & (FLAG_NORM_CONSIDER_TAUT &~FLAG_PROTON_CHARGE_CANCEL) ) {
            WarningMessage(pStrErrStruct, "Proton(s) added/removed");
        }

        if ( pINChI_Aux->bNormalizationFlags & FLAG_PROTON_CHARGE_CANCEL ) {
            WarningMessage(pStrErrStruct, "Charges neutralized");
        }
        ^^^*/
    }

    /* abs or rel stereo may establish one of two canonical numberings */
    if (( pCS->nLenLinearCTStereoCarb > 0 || pCS->nLenLinearCTStereoDble > 0 ) &&
          pCS->nLenCanonOrdStereo > 0 &&
          ( (pCS->LinearCTStereoCarb && pCS->LinearCTStereoCarbInv) ||
              (pCS->LinearCTStereoDble && pCS->LinearCTStereoDbleInv) ) &&
          pCS->nCanonOrdStereo && pCS->nCanonOrdStereoInv
       ) /* djb-rwth: addressing LLVM warnings */
    {

        pCanonRank = pCanonRankAtoms;
        pCanonOrd = pCS->nCanonOrdStereo;
        pCanonRankInv = pSortOrd;
        pCanonOrdInv = pCS->nCanonOrdStereoInv;
        Stereo = pINChI->Stereo;
        for (i = 0; i < num_at_tg; i++)
        {
            pCanonRankInv[pCanonOrdInv[i]] =
                pCanonRank[pCanonOrd[i]] = (AT_NUMB) ( i + 1 );
        }
        /********************************************************************/
        /* copy stereo bonds and stereo centers; compare Inv and Abs stereo */
        /********************************************************************/
        nErrorCode = CopyLinearCTStereoToINChIStereo( Stereo,
                           pCS->LinearCTStereoCarb, pCS->nLenLinearCTStereoCarb,
                           pCS->LinearCTStereoDble, pCS->nLenLinearCTStereoDble
                           , pCanonOrd, pCanonRank, at, 0 /* non-isotopic */
                           , pCS->LinearCTStereoCarbInv
                           , pCS->LinearCTStereoDbleInv
                           , pCanonOrdInv, pCanonRankInv );

        if (Stereo->t_parityInv && Stereo->nNumberInv)
        {
            if (nUserMode & REQ_MODE_RELATIVE_STEREO)
            {
                pINChI->nFlags |= INCHI_FLAG_REL_STEREO;
            }
            if (nUserMode & REQ_MODE_RACEMIC_STEREO)
            {
                pINChI->nFlags |= INCHI_FLAG_RAC_STEREO;
            }
            if (Stereo->nCompInv2Abs)
            {
                if (Stereo->nCompInv2Abs == -1)
                {
                    /* switch pointers in Stereo so that the stereo becomes the smallest (relative)  */
                    /* flag Stereo->nCompInv2Abs == -1 will keep track of this exchange */
                    AT_NUMB    *nNumberInv = Stereo->nNumberInv;
                    S_CHAR     *t_parityInv = Stereo->t_parityInv;
                    Stereo->nNumberInv = Stereo->nNumber;
                    Stereo->t_parityInv = Stereo->t_parity;
                    Stereo->nNumber = nNumberInv;
                    Stereo->t_parity = t_parityInv;
                    /* switch pointers to set rel. stereo to pINChI_Aux->nOrigAtNosInCanonOrd
                                       and inv. stereo to pINChI_Aux->nOrigAtNosInCanonOrdInv */
                    switch_ptrs( &pCanonRank, &pCanonRankInv );
                    switch_ptrs( &pCanonOrd, &pCanonOrdInv );
                    bUseNumberingInv = 1; /* use inverted stereo numbering instead of normal */
                }
            }
        }

        for (i = 0; i < num_atoms; i++)
        {
            pINChI_Aux->nOrigAtNosInCanonOrdInv[i] = at[pCanonOrdInv[i]].orig_at_number;
            pINChI_Aux->nOrigAtNosInCanonOrd[i] = at[pCanonOrd[i]].orig_at_number;
        }
        if (bUseNumberingInv)
        {
            /* switch ptrs back to avoid confusion */
            switch_ptrs( &pCanonRank, &pCanonRankInv );
            switch_ptrs( &pCanonOrd, &pCanonOrdInv );
            /* save inverted stereo ranks & order because it represents the smallest (relative) */
            memcpy(pCanonRank, pCanonRankInv, num_at_tg * sizeof(pCanonRank[0]));
            /* change pCS->nCanonOrdStereo[] to inverted: */
            memcpy(pCanonOrd, pCanonOrdInv, num_at_tg * sizeof(pCanonOrd[0]));
        }
        pCanonRankInv = NULL;
        pCanonOrdInv = NULL;
        pOrigNosInCanonOrd = NULL;
    }
    else
    { /*------------------------------ no stereo */

        pCanonOrd = pCS->nLenCanonOrdStereo > 0 ? pCS->nCanonOrdStereo :
            pCS->nLenCanonOrd > 0 ? pCS->nCanonOrd : NULL;
        pCanonRank = pCanonRankAtoms;
        pOrigNosInCanonOrd = pINChI_Aux->nOrigAtNosInCanonOrd;
        if (pCanonOrd && pCanonRank)
        {
            for (i = 0; i < num_atoms; i++)
            {
                pCanonRank[pCanonOrd[i]] = (AT_NUMB) ( i + 1 );
                pOrigNosInCanonOrd[i] = at[pCanonOrd[i]].orig_at_number;
            }
            for (; i < num_at_tg; i++)
            {
                pCanonRank[pCanonOrd[i]] = (AT_NUMB) ( i + 1 );
            }
        }
    }
    /*pCanonOrdNonIso = pCanonOrd;*/  /* save for aux info */


    if (pINChI_Aux->OrigInfo)
    {
        /* charges, radicals, valences */
        for (i = 0; i < num_atoms; i++)
        {
            if (pCanonOrd) /* djb-rwth: fixing a NULL pointer dereference */
            {
                ii = pCanonOrd[i];
                if (norm_at[ii].valence || norm_at[ii].num_H)
                {
                    pINChI_Aux->OrigInfo[i].cCharge = norm_at[ii].charge;
                    pINChI_Aux->OrigInfo[i].cRadical = (norm_at[ii].radical == RADICAL_SINGLET) ? 0 :
                        (norm_at[ii].radical == RADICAL_DOUBLET) ? 1 :
                        (norm_at[ii].radical == RADICAL_TRIPLET) ? 2 :
                        norm_at[ii].radical ? 3 : 0;
                    pINChI_Aux->OrigInfo[i].cUnusualValence =
                        get_unusual_el_valence(norm_at[ii].el_number, norm_at[ii].charge, norm_at[ii].radical,
                            norm_at[ii].chem_bonds_valence, norm_at[ii].num_H, norm_at[ii].valence);
                }
                else
                {
                    /* charge of a single atom component is in the INChI; valence = 0 is standard */
                    pINChI_Aux->OrigInfo[i].cRadical = (norm_at[ii].radical == RADICAL_SINGLET) ? 0 :
                        (norm_at[ii].radical == RADICAL_DOUBLET) ? 1 :
                        (norm_at[ii].radical == RADICAL_TRIPLET) ? 2 :
                        norm_at[ii].radical ? 3 : 0;
                }
            }
        }
    }

    /* non-isotopic canonical numbers and equivalence of atoms (Aux) */
    pConstitEquNumb = pINChI_Aux->nConstitEquNumbers;  /*  contitutional equivalence */
    pSymmRank = pCS->nSymmRank;
    if (pCanonOrd && pCanonRank && pSymmRank && pConstitEquNumb)
    {
        for (i = 0; i < num_atoms; i++)
        {
            pConstitEquNumb[i] = pSymmRank[pCanonOrd[i]]; /*  constit. equ. ranks in order of canonical numbers */
            pSortOrd[i] = i;
        }
        for (; i < num_at_tg; i++)
        {
            pSortOrd[i] = MAX_ATOMS; /* for debugging only */
        }

        pCG->m_pn_RankForSort = pConstitEquNumb;
        inchi_qsort( pCG, pSortOrd, num_atoms, sizeof( pSortOrd[0] ), CompRanksOrd );

        for (i = 0, nMinOrd = pSortOrd[0], j = 1; j <= num_atoms; j++)
        {
            if (j == num_atoms || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]])
            {
                nMinOrd++;
                if (j - i > 1)
                {
                    /*  found a sequence of equivalent atoms: i..j-1 */
                    while (i < j)
                    {
                        pConstitEquNumb[pSortOrd[i++]] = nMinOrd; /*  = min. canon. rank in the group of equ. atoms */
                    }
                    /*  at this point j == i */
                }
                else
                {
                    pConstitEquNumb[pSortOrd[i++]] = 0; /*  means the atom is not equivalent to any other */
                }
                nMinOrd = pSortOrd[j]; /*  at the end j = num_atoms */
            }
        }
    }
    else
    {
        nErrorCode |= ERR_NO_CANON_RESULTS;
        ret = -1;  /*  program error; no breakpoint here */
        goto exit_function;
    }
    /*  atomic numbers from the Periodic Table */
    for (i = 0; i < num_atoms; i++)
    {
        pINChI->nAtom[i] = (int) at[pCanonOrd[i]].el_number;
    }
    /*  connection table: atoms only (before 7-29-2003 pCS->LinearCT2 contained non-isotopic CT) */
    if (pCS->nLenLinearCTAtOnly <= 0 || !pCS->LinearCT || !pINChI->nConnTable)
    {
        nErrorCode |= ERR_NO_CANON_RESULTS;
        ret = -2;
        goto exit_function;
    }
    memcpy( pINChI->nConnTable, pCS->LinearCT, sizeof( pINChI->nConnTable[0] )*pCS->nLenLinearCTAtOnly );
    pINChI->lenConnTable = pCS->nLenLinearCTAtOnly;

    /*  tautomeric group(s) canonical representation */
    len = 0;
    if (bTautomeric &&
         0 < ( n = SortTautomerGroupsAndEndpoints( pCG, t_group_info,
             num_atoms, num_at_tg, pCanonRank ) ))
    {
        /* SortTautomerGroupsAndEndpoints() produces canonically ordered t-groups */
        pINChI->nFlags |= ( t_group_info->bTautFlagsDone & TG_FLAG_ALL_SALT_DONE ) ? INCHI_FLAG_ACID_TAUT : 0;
        /*  number of tautomeric groups */
        pINChI->nTautomer[len++] = (AT_NUMB) n;
        /* store each tautomeric group, one by one */
        for (i = 0; i < n; i++)
        {
            g = (int) t_group_info->tGroupNumber[i]; /* original group numbers in sorted order */
            t_group = t_group_info->t_group + g;    /* pointer to the tautomeric group */
            /*  NumAt+INCHI_T_NUM_MOVABLE (group length excluding this number) */
            pINChI->nTautomer[len++] = t_group->nNumEndpoints + INCHI_T_NUM_MOVABLE;
            /*  Num(H), Num(-) */
            for (j = 0; j < INCHI_T_NUM_MOVABLE; j++) /* djb-rwth: removing redundant code */
                pINChI->nTautomer[len++] = t_group->num[j];
            for (j = T_NUM_NO_ISOTOPIC; j < INCHI_T_NUM_MOVABLE; j++) /* djb-rwth: redundant code as the loop is never executed -- discussion required */ /* djb-rwth: unresolved issue -- revision required */
                pINChI->nTautomer[len++] = 0; /* should not happen */
            /* tautomeric group endpoint canonical numbers, pre-sorted in ascending order */
            for (j = (int) t_group->nFirstEndpointAtNoPos,
                  m = j + (int) t_group->nNumEndpoints; j < m; j++)
            {
                pINChI->nTautomer[len++] = pCanonRank[(int) t_group_info->nEndpointAtomNumber[j]]; /*  At[j] */
            }
        }
        pINChI->lenTautomer = len;
        pINChI_Aux->nNumberOfTGroups = n;
    }
    else
    {
        pINChI->lenTautomer = 0;
        pINChI_Aux->nNumberOfTGroups = 0;
        if (t_group_info && ( ( t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT ) ||
            (t_group_info->nNumIsotopicEndpoints > 1 &&
            ( t_group_info->bTautFlagsDone & ( TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE ) )) )
           ) /* djb-rwth: addressing LLVM warning */
        {
            /* only protons (re)moved or added */
            pINChI->lenTautomer = 1;
            pINChI->nTautomer[0] = 0;
        }
    }

    /*  number of H (excluding tautomeric) */
    if (pCS->nNum_H)
    {
        for (i = 0; i < num_atoms; i++)
        {
            pINChI->nNum_H[i] = pCS->nNum_H[i];
        }
    }

    /*  number of fixed H (tautomeric H in non-tautomeric representation) */
    if (pCS->nNum_H_fixed && !pINChI->lenTautomer)
    {
        for (i = 0; i < num_atoms; i++)
        {
            pINChI->nNum_H_fixed[i] = pCS->nNum_H_fixed[i];
            pINChI->nNum_H[i] += pCS->nNum_H_fixed[i];
        }
    }

    /***********************************************************
     *  tautomeric group(s) numbering and symmetry;
     *  should not depend on switching to rel. stereo numbering
     */
    if (pINChI->lenTautomer && ( n = pINChI_Aux->nNumberOfTGroups ))
    {
        pCanonOrdTaut = pCS->nLenCanonOrdStereoTaut > 0 ? pCS->nCanonOrdStereoTaut :
            pCS->nLenCanonOrdTaut > 0 ? pCS->nCanonOrdTaut : NULL;
        pConstitEquNumb = pINChI_Aux->nConstitEquTGroupNumbers;
        pSymmRank = pCS->nSymmRankTaut;
        if (pCanonOrdTaut && pSymmRank && pConstitEquNumb)
        {
            for (i = 0; i < n; i++)
            {
                pConstitEquNumb[i] = pSymmRank[pCanonOrdTaut[i]];
                pSortOrd[i] = i;
            }

            pCG->m_pn_RankForSort = pConstitEquNumb;
            inchi_qsort( pCG, pSortOrd, n, sizeof( pSortOrd[0] ), CompRanksOrd );

            for (i = 0, nMinOrd = pSortOrd[0], j = 1; j <= n; j++)
            {
                if (j == n || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]])
                {
                    nMinOrd++; /* make is start from 1, not from zero */
                    if (j - i > 1)
                    {
                        /*  found a sequence of more than one equivalent t-groups: i..j-1 */
                        while (i < j)
                        {
                            pConstitEquNumb[pSortOrd[i++]] = nMinOrd;
                        }
                    }
                    else
                    {
                        pConstitEquNumb[pSortOrd[i++]] = 0;
                    }
                    nMinOrd = pSortOrd[j]; /*  at the end j == n */
                }
            }
        }
    }

    /*  Allocate and fill Hill formula */
    if (!( pINChI->szHillFormula = AllocateAndFillHillFormula( pINChI ) ))
    {
        nErrorCode = 0;
        ret = CT_WRONG_FORMULA; /* CT_OUT_OF_RAM;*/  /*   <BRKPT> */
        pINChI->nErrorCode = pINChI_Aux->nErrorCode = ret;
        goto exit_function;
    }

    if ((nStereoUnmarkMode = UnmarkAllUndefinedUnknownStereo( pINChI->Stereo, nUserMode ))) /* djb-rwth: addressing LLVM warning */
    {
        pINChI->nFlags |= ( nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU ) ? INCHI_FLAG_SC_IGN_ALL_UU : 0;
        pINChI->nFlags |= ( nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU ) ? INCHI_FLAG_SB_IGN_ALL_UU : 0;
        if (( nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU ) ||
            ( nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU ))
        {
            WarningMessage( pStrErrStruct, "Omitted undefined stereo" );
        }
    }

    /*************************/
    /* mark ambiguous stereo */
    /*************************/
    MarkAmbiguousStereo( at, norm_at, 0 /* non-isotopic */, pCanonOrd,
           pCS->LinearCTStereoCarb, pCS->nLenLinearCTStereoCarb,
           pCS->LinearCTStereoDble, pCS->nLenLinearCTStereoDble );


    /************************************************************************
     *
     *  isotopic part
     */
    /* abs or rel stereo may establish one of two canonical numberings */

    if (( pCS->nLenLinearCTIsotopicStereoCarb > 0 || pCS->nLenLinearCTIsotopicStereoDble > 0 ) &&
          pCS->nLenCanonOrdIsotopicStereo > 0 &&
          ( (pCS->LinearCTIsotopicStereoCarb && pCS->LinearCTIsotopicStereoCarbInv) ||
              (pCS->LinearCTIsotopicStereoDble && pCS->LinearCTIsotopicStereoDbleInv) ) &&
          pCS->nCanonOrdIsotopicStereo    && pCS->nCanonOrdIsotopicStereoInv
          ) /* djb-rwth: addressing LLVM warnings */
    {
        /* found isotopic stereo */
        pCanonRank = pCanonRankAtoms;
        pCanonOrd = pCS->nCanonOrdIsotopicStereo;
        pCanonRankInv = pSortOrd;
        pCanonOrdInv = pCS->nCanonOrdIsotopicStereoInv;
        Stereo = pINChI->StereoIsotopic;
        for (i = 0; i < num_at_tg; i++)
        {
            pCanonRankInv[pCanonOrdInv[i]] =
                pCanonRank[pCanonOrd[i]] = (AT_NUMB) ( i + 1 );
        }
        /********************************************************************/
        /* copy stereo bonds and stereo centers; compare Inv and Abs stereo */
        /********************************************************************/
        nErrorCode = CopyLinearCTStereoToINChIStereo( Stereo,
                           pCS->LinearCTIsotopicStereoCarb, pCS->nLenLinearCTIsotopicStereoCarb,
                           pCS->LinearCTIsotopicStereoDble, pCS->nLenLinearCTIsotopicStereoDble
                           , pCanonOrd, pCanonRank, at, 1 /* isotopic */
                           , pCS->LinearCTIsotopicStereoCarbInv
                           , pCS->LinearCTIsotopicStereoDbleInv
                           , pCanonOrdInv, pCanonRankInv );

        if (Stereo->t_parityInv && Stereo->nNumberInv)
        {
            if (nUserMode & REQ_MODE_RELATIVE_STEREO)
            {
                pINChI->nFlags |= INCHI_FLAG_REL_STEREO;
            }
            if (nUserMode & REQ_MODE_RACEMIC_STEREO)
            {
                pINChI->nFlags |= INCHI_FLAG_RAC_STEREO;
            }
            if (Stereo->nCompInv2Abs)
            {
                if (Stereo->nCompInv2Abs == -1)
                {
                    /* switch pointers so that the stereo becomes the smallest (relative)  */
                    /* flag Stereo->nCompInv2Abs == -1 will keep track of this exchange */
                    AT_NUMB    *nNumberInv = Stereo->nNumberInv;
                    S_CHAR     *t_parityInv = Stereo->t_parityInv;
                    Stereo->nNumberInv = Stereo->nNumber;
                    Stereo->t_parityInv = Stereo->t_parity;
                    Stereo->nNumber = nNumberInv;
                    Stereo->t_parity = t_parityInv;
                    switch_ptrs( &pCanonRank, &pCanonRankInv );
                    switch_ptrs( &pCanonOrd, &pCanonOrdInv );
                    bUseIsotopicNumberingInv = 1;
                }
            }
        }

        for (i = 0; i < num_atoms; i++)
        {
            pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv[i] = at[pCanonOrdInv[i]].orig_at_number;
            pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[i] = at[pCanonOrd[i]].orig_at_number;
        }
        if (bUseIsotopicNumberingInv)
        {
            switch_ptrs( &pCanonRank, &pCanonRankInv );
            switch_ptrs( &pCanonOrd, &pCanonOrdInv );
            memcpy(pCanonRank, pCanonRankInv, num_at_tg * sizeof(pCanonRank[0]));
            memcpy(pCanonOrd, pCanonOrdInv, num_at_tg * sizeof(pCanonOrd[0]));
        }
        pCanonRankInv = NULL;
        pCanonOrdInv = NULL;
        pOrigNosInCanonOrd = NULL;
    }
    else
    {
        /* no isotopic stereo */
        pCanonOrd = pCS->nLenCanonOrdIsotopicStereo > 0 ? pCS->nCanonOrdIsotopicStereo :
            pCS->nLenCanonOrdIsotopic > 0 ? pCS->nCanonOrdIsotopic : NULL;
        pCanonRank = pCanonRankAtoms;
        pOrigNosInCanonOrd = pINChI_Aux->nIsotopicOrigAtNosInCanonOrd;
        if (pCanonOrd && pCanonRank)
        {
            for (i = 0; i < num_atoms; i++)
            { /* Fix13 -- out of bounds */
                pCanonRank[pCanonOrd[i]] = (AT_NUMB) ( i + 1 );
                pOrigNosInCanonOrd[i] = at[pCanonOrd[i]].orig_at_number;
            }
            for (; i < num_at_tg; i++)
            { /* Fix13 -- out of bounds */
                pCanonRank[pCanonOrd[i]] = (AT_NUMB) ( i + 1 );
            }
        }
    }
    /*pCanonOrdIso = pCanonOrd;*/

    pConstitEquNumb = pINChI_Aux->nConstitEquIsotopicNumbers;
    pSymmRank = pCS->nSymmRankIsotopic;
    if (pCanonOrd && pCanonRank && pConstitEquNumb && pSymmRank)
    {
        for (i = 0; i < num_atoms; i++)
        {
            pConstitEquNumb[i] = pSymmRank[pCanonOrd[i]];
            pSortOrd[i] = i;
        }
        for (; i < num_at_tg; i++)
        {
            pSortOrd[i] = i;
        }

        pCG->m_pn_RankForSort = pConstitEquNumb;
        inchi_qsort( pCG, pSortOrd, num_atoms, sizeof( pSortOrd[0] ), CompRanksOrd );

        for (i = 0, nMinOrd = pSortOrd[0], j = 1; j <= num_atoms; j++)
        {
            if (j == num_atoms || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]])
            {
                nMinOrd++;
                if (j - i > 1)
                {
                    /*  found a sequence of equivalent atoms: i..j-1 */
                    while (i < j)
                    {
                        pConstitEquNumb[pSortOrd[i++]] = nMinOrd;
                    }
                }
                else
                {
                    pConstitEquNumb[pSortOrd[i++]] = 0; /* nMinOrd; */
                }
                nMinOrd = pSortOrd[j];
            }
        }
    }
    else
    {
        goto exit_function; /*  no isotopic info available */
    }
    /*  isotopic atoms */
    n = pINChI->nNumberOfIsotopicAtoms = pCS->nLenLinearCTIsotopic;
    for (i = 0; i < n; i++)
    {
        pINChI->IsotopicAtom[i].nAtomNumber = pCS->LinearCTIsotopic[i].at_num;
        pINChI->IsotopicAtom[i].nIsoDifference = pCS->LinearCTIsotopic[i].iso_atw_diff;
        pINChI->IsotopicAtom[i].nNum_H = pCS->LinearCTIsotopic[i].num_1H;
        pINChI->IsotopicAtom[i].nNum_D = pCS->LinearCTIsotopic[i].num_D;
        pINChI->IsotopicAtom[i].nNum_T = pCS->LinearCTIsotopic[i].num_T;
    }
    /*  isotopic tautomeric groups */
    n = pINChI->nNumberOfIsotopicTGroups = pCS->nLenLinearCTIsotopicTautomer;
    for (i = 0; i < n; i++)
    {
        pINChI->IsotopicTGroup[i].nTGroupNumber = pCS->LinearCTIsotopicTautomer[i].tgroup_num;
        pINChI->IsotopicTGroup[i].nNum_H = pCS->LinearCTIsotopicTautomer[i].num[2];
        pINChI->IsotopicTGroup[i].nNum_D = pCS->LinearCTIsotopicTautomer[i].num[1];
        pINChI->IsotopicTGroup[i].nNum_T = pCS->LinearCTIsotopicTautomer[i].num[0];
    }
    /* atoms that may exchange isotopic H-atoms */
    if (pCS->nExchgIsoH && pINChI->nPossibleLocationsOfIsotopicH)
    {
        for (i = 0, j = 1; i < num_atoms; i++)
        {
            if (pCS->nExchgIsoH[i])
            {
                pINChI->nPossibleLocationsOfIsotopicH[j++] = (AT_NUMB) ( i + 1 ); /* canonical number */
            }
        }
        pINChI->nPossibleLocationsOfIsotopicH[0] = (AT_NUMB) j; /* length including the 0th element */
    }

    if ((nStereoUnmarkMode = UnmarkAllUndefinedUnknownStereo( pINChI->StereoIsotopic, nUserMode ))) /* djb-rwth: addressing LLVM warning */
    {
        pINChI->nFlags |= ( nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU ) ? INCHI_FLAG_SC_IGN_ALL_ISO_UU : 0;
        pINChI->nFlags |= ( nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU ) ? INCHI_FLAG_SC_IGN_ALL_ISO_UU : 0;
        if (( nStereoUnmarkMode & REQ_MODE_SC_IGN_ALL_UU ) ||
            ( nStereoUnmarkMode & REQ_MODE_SB_IGN_ALL_UU ))
        {
            WarningMessage( pStrErrStruct, "Omitted undefined stereo" );
        }
    }
    /* mark ambiguous stereo */
    MarkAmbiguousStereo( at, norm_at, 1 /* isotopic */, pCanonOrd,
           pCS->LinearCTIsotopicStereoCarb, pCS->nLenLinearCTIsotopicStereoCarb,
           pCS->LinearCTIsotopicStereoDble, pCS->nLenLinearCTIsotopicStereoDble );

    /***********************************************************
     *  isotopic tautomeric group(s) numbering and symmetry;
     *  should not depend on switching to rel. stereo numbering
     */
    if (pINChI->lenTautomer && pINChI_Aux->nConstitEquIsotopicTGroupNumbers && pCS->nSymmRankIsotopicTaut &&
        ( pCS->nLenLinearCTIsotopic || pCS->nLenLinearCTIsotopicTautomer ) &&
          t_group_info && t_group_info->num_t_groups > 0)
    {
        n = t_group_info->num_t_groups; /* djb-rwth: ignoring LLVM warning: value used */
        pCanonOrdTaut = pCS->nLenCanonOrdIsotopicStereoTaut > 0 ?
            ( n = pCS->nLenCanonOrdIsotopicStereoTaut, pCS->nCanonOrdIsotopicStereoTaut ) :
            pCS->nLenCanonOrdIsotopicTaut > 0 ?
            ( n = pCS->nLenCanonOrdIsotopicTaut, pCS->nCanonOrdIsotopicTaut ) : ( n = 0, (AT_RANK*) NULL );
        pConstitEquNumb = pINChI_Aux->nConstitEquIsotopicTGroupNumbers;
        pSymmRank = pCS->nSymmRankIsotopicTaut;
        if (pCanonOrdTaut && pSymmRank && pConstitEquNumb && n > 0)
        {
            for (i = 0; i < n; i++)
            {
                pConstitEquNumb[i] = pSymmRank[pCanonOrdTaut[i]];
                pSortOrd[i] = i;
            }

            pCG->m_pn_RankForSort = pConstitEquNumb;
            inchi_qsort( pCG, pSortOrd, n, sizeof( pSortOrd[0] ), CompRanksOrd );

            for (i = 0, nMinOrd = pSortOrd[0], j = 1; j <= n; j++)
            {
                if (j == n || pConstitEquNumb[pSortOrd[i]] != pConstitEquNumb[pSortOrd[j]])
                {
                    nMinOrd++;
                    if (j - i > 1)
                    {
                        /*  found a sequence of equivalent t-groups: i..j-1 */
                        while (i < j)
                        {
                            pConstitEquNumb[pSortOrd[i++]] = nMinOrd;
                        }
                    }
                    else
                    {
                        pConstitEquNumb[pSortOrd[i++]] = 0; /*  nMinOrd; */
                    }
                    nMinOrd = pSortOrd[j]; /*  at the end j = n */
                }
            }
        }
    }


exit_function:
    if (pCanonRankAtoms)
    {
        inchi_free( pCanonRankAtoms );
    }
    if (pSortOrd)
    {
        inchi_free( pSortOrd );
    }

    pINChI->nErrorCode |= nErrorCode;
    pINChI_Aux->nErrorCode |= nErrorCode;

    return ret;
}
    */
    // END INCHI C FUNCTION: FillOutINChIReducedWarn
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: FillOutINChIReducedWarn
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64.
    // INCHI✔️❌: WarningMessage is the active AddErrorMessage macro alias.
    // INCHI✔️❌: ReducedWarn suppresses normalization warnings, emits omitted-stereo warnings, and does not promote CopyLinearCTStereoToINChIStereo return 1 to CT_OUT_OF_RAM.
    // INCHI✔️❌: Shared typed SourceHeap implementation retains known allocation and bounds-check overhead.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: FillOutINChIReducedWarn

    FillOutINChIWithBehavior(
        heap,
        pINChI,
        pINChI_Aux,
        num_atoms,
        num_at_tg,
        num_removed_H,
        at,
        norm_at,
        pCS,
        pCG,
        bTautomeric,
        nUserMode,
        pStrErrStruct,
        FillOutINChIBehavior::ReducedWarn,
    )
}

#[allow(non_snake_case)]
pub(crate) fn CreateCompAtomData(
    heap: &mut SourceHeap,
    input_atom_data: &mut COMP_ATOM_DATA,
    num_atoms: i32,
    num_components: i32,
    intermediate_tautomer: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2232 CreateCompAtomData
    // INCHI✔️❌: int CreateCompAtomData( COMP_ATOM_DATA *inp_at_data,
    // INCHI✔️❌:                         int num_atoms,
    // INCHI✔️❌:                         int num_components,
    // INCHI✔️❌:                         int bIntermediateTaut )
    // INCHI✔️❌: {
    // INCHI✔️❌:     FreeCompAtomData( inp_at_data );
    // INCHI✔️❌:     if (( inp_at_data->at = CreateInpAtom( num_atoms ) ) &&
    // INCHI✔️❌:         ( num_components <= 1 || bIntermediateTaut ||
    // INCHI✔️❌:         ( inp_at_data->nOffsetAtAndH = (AT_NUMB*) inchi_calloc( 2 * ( (long long)num_components + 1 ), sizeof( inp_at_data->nOffsetAtAndH[0] ) ) ) )) /* djb-rwth: cast operator added */
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         inp_at_data->num_at = num_atoms;
    // INCHI✔️❌:         inp_at_data->num_components = ( num_components > 1 ) ? num_components : 0;
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     FreeCompAtomData( inp_at_data );
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CreateCompAtomData
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateCompAtomData
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; sizeof(long long) == 8; sizeof(AT_NUMB) == 2
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateCompAtomData

    FreeCompAtomData(heap, input_atom_data)?;
    input_atom_data.at = CreateInpAtom(heap, num_atoms)?;
    if !input_atom_data.at.is_null() {
        if num_components <= 1 || intermediate_tautomer != 0 {
            input_atom_data.num_at = num_atoms;
            input_atom_data.num_components = if num_components > 1 {
                num_components
            } else {
                0
            };
            return Ok(1);
        }

        let offset_count = 2_u64
            .checked_mul(
                u64::try_from(i64::from(num_components) + 1)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
            )
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        input_atom_data.nOffsetAtAndH =
            match inchi_calloc::<AT_NUMB>(heap, offset_count, SOURCE_SIZEOF_AT_NUMB) {
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                result => result?,
            };
        if !input_atom_data.nOffsetAtAndH.is_null() {
            input_atom_data.num_at = num_atoms;
            input_atom_data.num_components = num_components;
            return Ok(1);
        }
    }

    FreeCompAtomData(heap, input_atom_data)?;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CreateCompositeNormAtom(
    heap: &mut SourceHeap,
    composite_norm_data: &mut [COMP_ATOM_DATA; TAUT_NUM as usize + 1],
    all_inp_norm_data: &[INP_ATOM_DATA2],
    num_components: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1973 CreateCompositeNormAtom
    // INCHI✔️❌: int CreateCompositeNormAtom( COMP_ATOM_DATA *composite_norm_data,
    // INCHI✔️❌:                              INP_ATOM_DATA2 *all_inp_norm_data,
    // INCHI✔️❌:                              int num_components )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, jj, k, n, m, tot_num_at, tot_num_H, cur_num_at, cur_num_H; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int num_comp[TAUT_NUM + 1], num_taut[TAUT_NUM + 1], num_del[TAUT_NUM + 1], num_at[TAUT_NUM + 1], num_inp_at[TAUT_NUM + 1];
    // INCHI✔️❌:     int ret = 0, indicator = 1;
    // INCHI✔️❌:     inp_ATOM *at, *at_from;
    // INCHI✔️❌:     memset( num_comp, 0, sizeof( num_comp ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( num_taut, 0, sizeof( num_taut ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     memset( num_del, 0, sizeof( num_taut ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     /* count taut and non-taut components */
    // INCHI✔️❌:     for (j = 0; j < TAUT_NUM; j++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_comp[j] = num_taut[j] = 0;
    // INCHI✔️❌:         for (i = 0; i < num_components; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (all_inp_norm_data[i][j].bExists)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num_del[j] += ( 0 != all_inp_norm_data[i][j].bDeleted );
    // INCHI✔️❌:                 num_comp[j] ++;
    // INCHI✔️❌:                 num_taut[j] += ( 0 != all_inp_norm_data[i][j].bTautomeric );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count intermediate taut structure components */
    // INCHI✔️❌:     if (num_comp[TAUT_YES] > num_del[TAUT_YES] && num_taut[TAUT_YES])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         num_comp[TAUT_INI] = num_comp[TAUT_YES] - num_del[TAUT_YES];
    // INCHI✔️❌:         */
    // INCHI✔️❌:
    // INCHI✔️❌:         for (i = 0, j = TAUT_YES; i < num_components; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (all_inp_norm_data[i][j].bExists &&
    // INCHI✔️❌:                 ( all_inp_norm_data[i][j].bDeleted ||
    // INCHI✔️❌:                     (all_inp_norm_data[i][j].bTautomeric &&
    // INCHI✔️❌:                     all_inp_norm_data[i][j].at_fixed_bonds &&
    // INCHI✔️❌:                     all_inp_norm_data[i][j].bTautPreprocessed) )) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num_comp[TAUT_INI] ++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* count atoms and allocate composite atom data */
    // INCHI✔️❌:     for (jj = 0; jj <= TAUT_INI; jj++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_at[jj] = num_inp_at[jj] = 0;
    // INCHI✔️❌:         j = inchi_min( jj, TAUT_YES );
    // INCHI✔️❌:         if (num_comp[jj])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < num_components; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (all_inp_norm_data[i][j].bDeleted)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* find k = the normaized structure index */
    // INCHI✔️❌:                 if (jj == TAUT_INI)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (all_inp_norm_data[i][j].bExists &&
    // INCHI✔️❌:                          all_inp_norm_data[i][j].at_fixed_bonds)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k = j;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (all_inp_norm_data[i][ALT_TAUT( j )].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted &&
    // INCHI✔️❌:                                 !all_inp_norm_data[i][j].bDeleted)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k = ALT_TAUT( j );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (all_inp_norm_data[i][j].bExists)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 k = j;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (all_inp_norm_data[i][j].bExists)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k = j;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (all_inp_norm_data[i][ALT_TAUT( j )].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k = ALT_TAUT( j );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 num_inp_at[jj] += all_inp_norm_data[i][k].num_at; /* all atoms including terminal H */
    // INCHI✔️❌:                 num_at[jj] += all_inp_norm_data[i][k].num_at - all_inp_norm_data[i][k].num_removed_H;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (num_inp_at[jj])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!CreateCompAtomData( composite_norm_data + jj, num_inp_at[jj], num_components, jj == TAUT_INI ))
    // INCHI✔️❌:                     goto exit_error;
    // INCHI✔️❌:                 composite_norm_data[jj].num_removed_H = num_inp_at[jj] - num_at[jj];
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* fill out composite atom */
    // INCHI✔️❌:     for (jj = 0; jj <= TAUT_INI; jj++, indicator <<= 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         j = inchi_min( jj, TAUT_YES );
    // INCHI✔️❌:         if (num_comp[jj])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             tot_num_at = 0;
    // INCHI✔️❌:             tot_num_H = 0;
    // INCHI✔️❌:             for (i = 0; i < num_components; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (all_inp_norm_data[i][j].bDeleted)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     composite_norm_data[jj].nNumRemovedProtons += all_inp_norm_data[i][j].nNumRemovedProtons;
    // INCHI✔️❌:                     for (n = 0; n < NUM_H_ISOTOPES; n++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] += all_inp_norm_data[i][j].nNumRemovedProtonsIsotopic[n];
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                 k = TAUT_NUM; /* djb-rwth: ignoring LLVM warning: value used */
    // INCHI✔️❌:                 /* find k = the normaized structure index */
    // INCHI✔️❌:                 if (jj == TAUT_INI)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (all_inp_norm_data[i][j].bExists && all_inp_norm_data[i][j].at_fixed_bonds)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k = j;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (all_inp_norm_data[i][ALT_TAUT( j )].bExists)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k = ALT_TAUT( j );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (all_inp_norm_data[i][j].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 k = j;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (all_inp_norm_data[i][j].bExists)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         k = j;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (all_inp_norm_data[i][ALT_TAUT( j )].bExists && !all_inp_norm_data[i][ALT_TAUT( j )].bDeleted)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             k = ALT_TAUT( j );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* copy main atoms */
    // INCHI✔️❌:                 cur_num_H = all_inp_norm_data[i][k].num_removed_H;       /* number of terminal H atoms */
    // INCHI✔️❌:                 cur_num_at = all_inp_norm_data[i][k].num_at - cur_num_H;  /* number of all but explicit terminal H atoms */
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (( tot_num_at + cur_num_at ) > num_at[jj] ||
    // INCHI✔️❌:                     ( num_at[jj] + tot_num_H + cur_num_H ) > num_inp_at[jj])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     goto exit_error; /* miscount */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 at = composite_norm_data[jj].at + tot_num_at; /* points to the 1st destination atom */
    // INCHI✔️❌:                 at_from = ( jj == TAUT_INI && k == TAUT_YES && all_inp_norm_data[i][k].at_fixed_bonds ) ?
    // INCHI✔️❌:                     all_inp_norm_data[i][k].at_fixed_bonds : all_inp_norm_data[i][k].at;
    // INCHI✔️❌:                 memcpy( at, at_from, sizeof( composite_norm_data[0].at[0] ) * cur_num_at ); /* copy atoms except terminal H */
    // INCHI✔️❌:                 /* shift neighbors of main atoms */
    // INCHI✔️❌:                 for (n = 0; n < cur_num_at; n++, at++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (m = 0; m < at->valence; m++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         at->neighbor[m] += tot_num_at;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* copy explicit H */
    // INCHI✔️❌:                 if (cur_num_H)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     at = composite_norm_data[jj].at + num_at[jj] + tot_num_H; /* points to the 1st destination atom */
    // INCHI✔️❌:                     memcpy( at, at_from + cur_num_at, sizeof( composite_norm_data[0].at[0] ) * cur_num_H );
    // INCHI✔️❌:                     /* shift neighbors of explicit H atoms */
    // INCHI✔️❌:                     for (n = 0; n < cur_num_H; n++, at++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         for (m = 0; m < at->valence; m++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             at->neighbor[m] += tot_num_at;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* composite counts */
    // INCHI✔️❌:                 composite_norm_data[jj].bHasIsotopicLayer |= all_inp_norm_data[i][k].bHasIsotopicLayer;
    // INCHI✔️❌:                 composite_norm_data[jj].num_isotopic += all_inp_norm_data[i][k].num_isotopic;
    // INCHI✔️❌:                 composite_norm_data[jj].num_bonds += all_inp_norm_data[i][k].num_bonds;
    // INCHI✔️❌:                 composite_norm_data[jj].bTautomeric += ( j == jj ) && all_inp_norm_data[i][k].bTautomeric;
    // INCHI✔️❌:                 composite_norm_data[jj].nNumRemovedProtons += all_inp_norm_data[i][k].nNumRemovedProtons;
    // INCHI✔️❌:                 for (n = 0; n < NUM_H_ISOTOPES; n++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] += all_inp_norm_data[i][k].nNumRemovedProtonsIsotopic[n];
    // INCHI✔️❌:                     composite_norm_data[jj].num_iso_H[n] += all_inp_norm_data[i][k].num_iso_H[n];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /*
    // INCHI✔️❌:                 composite_norm_data[j].num_at            += cur_num_at + cur_num_H;
    // INCHI✔️❌:                 composite_norm_data[j].num_removed_H     += cur_num_H;
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 /* total count */
    // INCHI✔️❌:                 tot_num_at += cur_num_at;
    // INCHI✔️❌:                 tot_num_H += cur_num_H;
    // INCHI✔️❌:                 /* offset for the next component */
    // INCHI✔️❌:                 if (composite_norm_data[jj].nOffsetAtAndH)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     composite_norm_data[jj].nOffsetAtAndH[2 * i] = tot_num_at;
    // INCHI✔️❌:                     composite_norm_data[jj].nOffsetAtAndH[2 * i + 1] = num_at[jj] + tot_num_H;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (tot_num_at != num_at[jj] ||
    // INCHI✔️❌:                  num_at[jj] + tot_num_H != num_inp_at[jj])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_error; /* miscount */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             composite_norm_data[jj].bExists = ( tot_num_at > 0 );
    // INCHI✔️❌:             ret |= indicator;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌:
    // INCHI✔️❌: exit_error:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CreateCompositeNormAtom
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateCompositeNormAtom
    // INCHI✔️❌: #define ALT_TAUT(X) ((X)>TAUT_YES? TAUT_YES : 1-(X))
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; TAUT_NUM == 2; TAUT_YES == 1; TAUT_INI == 2
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CreateCompositeNormAtom

    let component_count = if num_components > 0 {
        usize::try_from(num_components).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
    } else {
        0
    };
    if component_count > all_inp_norm_data.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    let taut_num = TAUT_NUM as usize;
    let taut_yes = TAUT_YES as usize;
    let taut_ini = TAUT_INI as usize;
    let mut num_comp = [0_i32; TAUT_NUM as usize + 1];
    let mut num_taut = [0_i32; TAUT_NUM as usize + 1];
    let mut num_del = [0_i32; TAUT_NUM as usize + 1];
    let mut num_at = [0_i32; TAUT_NUM as usize + 1];
    let mut num_inp_at = [0_i32; TAUT_NUM as usize + 1];
    let mut ret = 0_i32;

    for j in 0..taut_num {
        num_comp[j] = 0;
        num_taut[j] = 0;
        for component in &all_inp_norm_data[..component_count] {
            if component[j].bExists != 0 {
                num_del[j] = num_del[j].wrapping_add(i32::from(component[j].bDeleted != 0));
                num_comp[j] = num_comp[j].wrapping_add(1);
                num_taut[j] = num_taut[j].wrapping_add(i32::from(component[j].bTautomeric != 0));
            }
        }
    }

    if num_comp[taut_yes] > num_del[taut_yes] && num_taut[taut_yes] != 0 {
        let j = taut_yes;
        for component in &all_inp_norm_data[..component_count] {
            if component[j].bExists != 0
                && (component[j].bDeleted != 0
                    || (component[j].bTautomeric != 0
                        && !component[j].at_fixed_bonds.is_null()
                        && component[j].bTautPreprocessed != 0))
            {
                num_comp[taut_ini] = num_comp[taut_ini].wrapping_add(1);
            }
        }
    }

    for jj in 0..=taut_ini {
        num_at[jj] = 0;
        num_inp_at[jj] = 0;
        let j = jj.min(taut_yes);
        let alternate = if j > taut_yes { taut_yes } else { 1 - j };
        if num_comp[jj] != 0 {
            for component in &all_inp_norm_data[..component_count] {
                if component[j].bDeleted != 0 {
                    continue;
                }
                let k;
                if jj == taut_ini {
                    if component[j].bExists != 0 && !component[j].at_fixed_bonds.is_null() {
                        k = j;
                    } else if component[alternate].bExists != 0
                        && component[alternate].bDeleted == 0
                        && component[j].bDeleted == 0
                    {
                        k = alternate;
                    } else if component[j].bExists != 0 {
                        k = j;
                    } else {
                        continue;
                    }
                } else if component[j].bExists != 0 {
                    k = j;
                } else if component[alternate].bExists != 0 && component[alternate].bDeleted == 0 {
                    k = alternate;
                } else {
                    continue;
                }
                num_inp_at[jj] = num_inp_at[jj].wrapping_add(component[k].num_at);
                num_at[jj] = num_at[jj]
                    .wrapping_add(component[k].num_at.wrapping_sub(component[k].num_removed_H));
            }
            if num_inp_at[jj] != 0 {
                if CreateCompAtomData(
                    heap,
                    &mut composite_norm_data[jj],
                    num_inp_at[jj],
                    num_components,
                    i32::from(jj == taut_ini),
                )? == 0
                {
                    return Ok(ret);
                }
                composite_norm_data[jj].num_removed_H = num_inp_at[jj].wrapping_sub(num_at[jj]);
            }
        }
    }

    let mut indicator = 1_i32;
    for jj in 0..=taut_ini {
        let j = jj.min(taut_yes);
        let alternate = if j > taut_yes { taut_yes } else { 1 - j };
        if num_comp[jj] != 0 {
            let mut total_atoms = 0_i32;
            let mut total_hydrogens = 0_i32;
            for (i, component) in all_inp_norm_data[..component_count].iter().enumerate() {
                if component[j].bDeleted != 0 {
                    composite_norm_data[jj].nNumRemovedProtons = composite_norm_data[jj]
                        .nNumRemovedProtons
                        .wrapping_add(component[j].nNumRemovedProtons);
                    for n in 0..composite_norm_data[jj].nNumRemovedProtonsIsotopic.len() {
                        composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] =
                            composite_norm_data[jj].nNumRemovedProtonsIsotopic[n]
                                .wrapping_add(component[j].nNumRemovedProtonsIsotopic[n]);
                    }
                    continue;
                }

                let k;
                if jj == taut_ini {
                    if component[j].bExists != 0 && !component[j].at_fixed_bonds.is_null() {
                        k = j;
                    } else if component[alternate].bExists != 0 {
                        k = alternate;
                    } else if component[j].bExists != 0 && component[alternate].bDeleted == 0 {
                        k = j;
                    } else {
                        continue;
                    }
                } else if component[j].bExists != 0 {
                    k = j;
                } else if component[alternate].bExists != 0 && component[alternate].bDeleted == 0 {
                    k = alternate;
                } else {
                    continue;
                }

                let current_hydrogens = component[k].num_removed_H;
                let current_atoms = component[k].num_at.wrapping_sub(current_hydrogens);
                if total_atoms.wrapping_add(current_atoms) > num_at[jj]
                    || num_at[jj]
                        .wrapping_add(total_hydrogens)
                        .wrapping_add(current_hydrogens)
                        > num_inp_at[jj]
                {
                    return Ok(ret);
                }

                let source_pointer =
                    if jj == taut_ini && k == taut_yes && !component[k].at_fixed_bonds.is_null() {
                        component[k].at_fixed_bonds
                    } else {
                        component[k].at
                    };
                let source_count = usize::try_from(component[k].num_at)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let source_atoms = heap
                    .slice(source_pointer.as_const())?
                    .get(..source_count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                let current_atom_count = usize::try_from(current_atoms)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let current_hydrogen_count = usize::try_from(current_hydrogens)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let main_start = usize::try_from(total_atoms)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let hydrogen_start = usize::try_from(num_at[jj].wrapping_add(total_hydrogens))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;

                {
                    let destination = heap.slice_mut(composite_norm_data[jj].at)?;
                    let main_end = main_start
                        .checked_add(current_atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let main = destination
                        .get_mut(main_start..main_end)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    main.clone_from_slice(&source_atoms[..current_atom_count]);
                    for atom in main {
                        let valence = usize::try_from(atom.valence.max(0))
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        if valence > atom.neighbor.len() {
                            return Err(SourceHeapError::PointerOutOfBounds);
                        }
                        for neighbor in &mut atom.neighbor[..valence] {
                            *neighbor = i32::from(*neighbor).wrapping_add(total_atoms) as AT_NUMB;
                        }
                    }

                    if current_hydrogens != 0 {
                        let hydrogen_end = hydrogen_start
                            .checked_add(current_hydrogen_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let hydrogens = destination
                            .get_mut(hydrogen_start..hydrogen_end)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        hydrogens.clone_from_slice(
                            &source_atoms
                                [current_atom_count..current_atom_count + current_hydrogen_count],
                        );
                        for atom in hydrogens {
                            let valence = usize::try_from(atom.valence.max(0))
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            if valence > atom.neighbor.len() {
                                return Err(SourceHeapError::PointerOutOfBounds);
                            }
                            for neighbor in &mut atom.neighbor[..valence] {
                                *neighbor =
                                    i32::from(*neighbor).wrapping_add(total_atoms) as AT_NUMB;
                            }
                        }
                    }
                }

                composite_norm_data[jj].bHasIsotopicLayer |= component[k].bHasIsotopicLayer;
                composite_norm_data[jj].num_isotopic = composite_norm_data[jj]
                    .num_isotopic
                    .wrapping_add(component[k].num_isotopic);
                composite_norm_data[jj].num_bonds = composite_norm_data[jj]
                    .num_bonds
                    .wrapping_add(component[k].num_bonds);
                composite_norm_data[jj].bTautomeric = composite_norm_data[jj]
                    .bTautomeric
                    .wrapping_add(i32::from(j == jj && component[k].bTautomeric != 0));
                composite_norm_data[jj].nNumRemovedProtons = composite_norm_data[jj]
                    .nNumRemovedProtons
                    .wrapping_add(component[k].nNumRemovedProtons);
                for n in 0..composite_norm_data[jj].nNumRemovedProtonsIsotopic.len() {
                    composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] = composite_norm_data[jj]
                        .nNumRemovedProtonsIsotopic[n]
                        .wrapping_add(component[k].nNumRemovedProtonsIsotopic[n]);
                    composite_norm_data[jj].num_iso_H[n] = composite_norm_data[jj].num_iso_H[n]
                        .wrapping_add(component[k].num_iso_H[n]);
                }
                total_atoms = total_atoms.wrapping_add(current_atoms);
                total_hydrogens = total_hydrogens.wrapping_add(current_hydrogens);
                if !composite_norm_data[jj].nOffsetAtAndH.is_null() {
                    let offsets = heap.slice_mut(composite_norm_data[jj].nOffsetAtAndH)?;
                    offsets[2 * i] = total_atoms as AT_NUMB;
                    offsets[2 * i + 1] = num_at[jj].wrapping_add(total_hydrogens) as AT_NUMB;
                }
            }
            if total_atoms != num_at[jj]
                || num_at[jj].wrapping_add(total_hydrogens) != num_inp_at[jj]
            {
                return Ok(ret);
            }
            composite_norm_data[jj].bExists = i32::from(total_atoms > 0);
            ret |= indicator;
        }
        indicator <<= 1;
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn make_norm_atoms_from_inp_atoms(
    heap: &mut SourceHeap,
    generation_data: &mut INCHIGEN_DATA,
    control: &INCHIGEN_CONTROL,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2937 make_norm_atoms_from_inp_atoms
    // INCHI✔️❌: void make_norm_atoms_from_inp_atoms( INCHIGEN_DATA *gendata,
    // INCHI✔️❌:                                      INCHIGEN_CONTROL *genctl )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* TODO: make a full copy (with allocs) of atom arrays */
    // INCHI✔️❌:     size_t t1;
    // INCHI✔️❌:     int k;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (k = 0; k < INCHI_NUM; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (NULL != genctl->InpNormAtData[k])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             t1 = genctl->StructData.num_components[k] * sizeof( NORM_ATOMS );
    // INCHI✔️❌:             memcpy( gendata->NormAtomsNontaut[k], genctl->InpNormAtData[k], t1 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (NULL != genctl->InpNormTautData[k])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             t1 = genctl->StructData.num_components[k] * sizeof( NORM_ATOMS );
    // INCHI✔️❌:             memcpy( gendata->NormAtomsTaut[k], genctl->InpNormTautData[k], t1 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: make_norm_atoms_from_inp_atoms
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: make_norm_atoms_from_inp_atoms
    // INCHI✔️❌: INCHI_NUM is 2 in the selected TARGET_API_LIB configuration.
    // INCHI✔️❌: NORM_ATOMS and INP_ATOM_DATA have the same active C layout; memcpy is shallow.
    // INCHI✔️❌: Typed SourceHeap conversion allocates a temporary Vec, retaining known overhead.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: make_norm_atoms_from_inp_atoms

    for k in 0..INCHI_NUM as usize {
        let mut copy_shallow = |destination: SourceMutPointer<NORM_ATOMS>,
                                source: SourceMutPointer<INP_ATOM_DATA>|
         -> Result<(), SourceHeapError> {
            if source.is_null() {
                return Ok(());
            }
            let count = usize::try_from(control.StructData.num_components[k])
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let source = heap
                .slice(source.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .iter()
                .map(|value| NORM_ATOMS {
                    at: value.at.cast(),
                    at_fixed_bonds: value.at_fixed_bonds.cast(),
                    num_at: value.num_at,
                    num_removed_H: value.num_removed_H,
                    num_bonds: value.num_bonds,
                    num_isotopic: value.num_isotopic,
                    bExists: value.bExists,
                    bDeleted: value.bDeleted,
                    bHasIsotopicLayer: value.bHasIsotopicLayer,
                    bTautomeric: value.bTautomeric,
                    bTautPreprocessed: value.bTautPreprocessed,
                    nNumRemovedProtons: value.nNumRemovedProtons,
                    nNumRemovedProtonsIsotopic: value.nNumRemovedProtonsIsotopic,
                    num_iso_H: value.num_iso_H,
                    bTautFlags: value.bTautFlags,
                    bTautFlagsDone: value.bTautFlagsDone,
                    bNormalizationFlags: value.bNormalizationFlags,
                })
                .collect::<Vec<_>>();
            heap.slice_mut(destination)?
                .get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone_from_slice(&source);
            Ok(())
        };

        copy_shallow(
            generation_data.NormAtomsNontaut[k],
            control.InpNormAtData[k],
        )?;
        copy_shallow(generation_data.NormAtomsTaut[k], control.InpNormTautData[k])?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::base::mol2atom::FreeCompAtomData;
    use crate::source_types::{
        BOND_DOUBLE, BOND_SINGLE, CT_ATOMCOUNT_ERR, FLAG_PROTON_SINGLE_REMOVED,
        TG_FLAG_TEST_TAUT__ATOMS, TG_FLAG_TEST_TAUT__ATOMS_DONE, TG_FLAG_VARIABLE_PROTONS,
        inp_ATOM,
    };

    fn atom(number: AT_NUMB, neighbor: AT_NUMB) -> inp_ATOM {
        inp_ATOM {
            orig_at_number: number,
            neighbor: {
                let mut neighbors = [0; 20];
                neighbors[0] = neighbor;
                neighbors
            },
            valence: 1,
            ..inp_ATOM::default()
        }
    }

    #[test]
    fn source_port__inchi_dll_a2__make_norm_atoms_from_inp_atoms__line_2937() {
        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let fixed_atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let first = INP_ATOM_DATA {
            at: atoms,
            at_fixed_bonds: fixed_atoms,
            num_at: 1,
            num_removed_H: 2,
            num_bonds: 3,
            num_isotopic: 4,
            bExists: 5,
            bDeleted: 6,
            bHasIsotopicLayer: 7,
            bTautomeric: 8,
            bTautPreprocessed: 9,
            nNumRemovedProtons: 10,
            nNumRemovedProtonsIsotopic: [11, 12, 13],
            num_iso_H: [14, 15, 16],
            bTautFlags: 17,
            bTautFlagsDone: 18,
            bNormalizationFlags: 19,
        };
        let second = INP_ATOM_DATA {
            num_at: -1,
            num_removed_H: -2,
            num_bonds: -3,
            num_isotopic: -4,
            bExists: -5,
            bDeleted: -6,
            bHasIsotopicLayer: -7,
            bTautomeric: -8,
            bTautPreprocessed: -9,
            nNumRemovedProtons: -10,
            nNumRemovedProtonsIsotopic: [-11, -12, -13],
            num_iso_H: [-14, -15, -16],
            bTautFlags: u64::MAX - 2,
            bTautFlagsDone: u64::MAX - 1,
            bNormalizationFlags: u64::MAX,
            ..INP_ATOM_DATA::default()
        };
        let expected_first = NORM_ATOMS {
            at: atoms.cast(),
            at_fixed_bonds: fixed_atoms.cast(),
            num_at: 1,
            num_removed_H: 2,
            num_bonds: 3,
            num_isotopic: 4,
            bExists: 5,
            bDeleted: 6,
            bHasIsotopicLayer: 7,
            bTautomeric: 8,
            bTautPreprocessed: 9,
            nNumRemovedProtons: 10,
            nNumRemovedProtonsIsotopic: [11, 12, 13],
            num_iso_H: [14, 15, 16],
            bTautFlags: 17,
            bTautFlagsDone: 18,
            bNormalizationFlags: 19,
        };
        let expected_second = NORM_ATOMS {
            num_at: -1,
            num_removed_H: -2,
            num_bonds: -3,
            num_isotopic: -4,
            bExists: -5,
            bDeleted: -6,
            bHasIsotopicLayer: -7,
            bTautomeric: -8,
            bTautPreprocessed: -9,
            nNumRemovedProtons: -10,
            nNumRemovedProtonsIsotopic: [-11, -12, -13],
            num_iso_H: [-14, -15, -16],
            bTautFlags: u64::MAX - 2,
            bTautFlagsDone: u64::MAX - 1,
            bNormalizationFlags: u64::MAX,
            ..NORM_ATOMS::default()
        };
        let sentinel = NORM_ATOMS {
            num_at: 777,
            bNormalizationFlags: 888,
            ..NORM_ATOMS::default()
        };

        let mut control = INCHIGEN_CONTROL::default();
        control.StructData.num_components = [2, 1];
        control.InpNormAtData[0] = heap
            .allocate_model_storage(vec![first.clone(), second.clone()])
            .unwrap();
        control.InpNormTautData[0] = heap
            .allocate_model_storage(vec![second.clone(), first.clone()])
            .unwrap();
        control.InpNormTautData[1] = heap.allocate_model_storage(vec![first]).unwrap();

        let mut generation = INCHIGEN_DATA {
            num_components: [91, 92],
            NormAtomsNontaut: [
                heap.allocate_model_storage(vec![sentinel.clone(); 3])
                    .unwrap(),
                heap.allocate_model_storage(vec![sentinel.clone()]).unwrap(),
            ],
            NormAtomsTaut: [
                heap.allocate_model_storage(vec![sentinel.clone(); 3])
                    .unwrap(),
                heap.allocate_model_storage(vec![sentinel.clone(); 2])
                    .unwrap(),
            ],
            ..INCHIGEN_DATA::default()
        };

        assert_eq!(
            make_norm_atoms_from_inp_atoms(&mut heap, &mut generation, &control),
            Ok(())
        );
        assert_eq!(generation.num_components, [91, 92]);
        assert_eq!(
            heap.slice(generation.NormAtomsNontaut[0].as_const())
                .unwrap(),
            &[
                expected_first.clone(),
                expected_second.clone(),
                sentinel.clone()
            ]
        );
        assert_eq!(
            heap.slice(generation.NormAtomsTaut[0].as_const()).unwrap(),
            &[expected_second, expected_first.clone(), sentinel.clone()]
        );
        assert_eq!(
            heap.slice(generation.NormAtomsNontaut[1].as_const())
                .unwrap(),
            &[sentinel.clone()]
        );
        assert_eq!(
            heap.slice(generation.NormAtomsTaut[1].as_const()).unwrap(),
            &[expected_first, sentinel]
        );
    }

    fn input_data(
        heap: &mut SourceHeap,
        atoms: Vec<inp_ATOM>,
        removed_hydrogens: i32,
    ) -> crate::source_types::INP_ATOM_DATA {
        let atom_count = i32::try_from(atoms.len()).unwrap();
        crate::source_types::INP_ATOM_DATA {
            at: heap.allocate_model_storage(atoms).unwrap(),
            num_at: atom_count,
            num_removed_H: removed_hydrogens,
            bExists: 1,
            ..crate::source_types::INP_ATOM_DATA::default()
        }
    }

    fn free_composites(
        heap: &mut SourceHeap,
        composites: &mut [COMP_ATOM_DATA; TAUT_NUM as usize + 1],
    ) {
        for composite in composites {
            FreeCompAtomData(heap, composite).unwrap();
        }
    }

    fn normalization_atoms() -> Vec<inp_ATOM> {
        let mut atoms = vec![
            inp_ATOM {
                el_number: 6,
                orig_at_number: 1,
                chem_bonds_valence: 1,
                num_H: 3,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 6,
                orig_at_number: 2,
                chem_bonds_valence: 1,
                num_H: 3,
                ..inp_ATOM::default()
            },
        ];
        atoms[0].neighbor[0] = 1;
        atoms[0].bond_type[0] = BOND_SINGLE as u8;
        atoms[0].valence = 1;
        atoms[1].neighbor[0] = 0;
        atoms[1].bond_type[0] = BOND_SINGLE as u8;
        atoms[1].valence = 1;
        atoms
    }

    fn normalization_outputs(
        heap: &mut SourceHeap,
        non_taut: bool,
        taut: bool,
        fixed_bonds: bool,
        atom_count: usize,
    ) -> [INP_ATOM_DATA; TAUT_NUM as usize] {
        let mut outputs = std::array::from_fn(|_| INP_ATOM_DATA::default());
        for (index, enabled) in [non_taut, taut].into_iter().enumerate() {
            if enabled {
                outputs[index].at = heap
                    .allocate_model_storage(vec![inp_ATOM::default(); atom_count])
                    .unwrap();
                outputs[index].num_at = atom_count as i32;
                outputs[index].bExists = 1;
            }
        }
        if taut && fixed_bonds {
            outputs[TAUT_YES as usize].at_fixed_bonds = heap
                .allocate_model_storage(vec![inp_ATOM::default(); atom_count])
                .unwrap();
        }
        outputs
    }

    fn normalization_result(
        heap: &mut SourceHeap,
        outputs: &mut [INP_ATOM_DATA; TAUT_NUM as usize],
        input: SourceMutPointer<inp_ATOM>,
        mode: INCHI_MODE,
        flags: INCHI_MODE,
        z: &mut COMPONENT_TREAT_INFO,
        num_inp_at: i32,
    ) -> Result<i32, SourceHeapError> {
        let inchi = [
            heap.allocate_model_storage(vec![INChI::default()])?,
            heap.allocate_model_storage(vec![INChI::default()])?,
        ];
        let aux = [
            heap.allocate_model_storage(vec![INChI_Aux::default()])?,
            heap.allocate_model_storage(vec![INChI_Aux::default()])?,
        ];
        let clock = heap.allocate_model_storage(vec![INCHI_CLOCK::default()])?;
        z.nUserMode = mode;
        let mut taut_flags = flags;
        let mut taut_flags_done = 0;
        Normalization_step(
            heap,
            &mut CANON_GLOBALS::default(),
            clock,
            &inchi,
            &aux,
            input,
            outputs,
            num_inp_at,
            SourceMutPointer::null(),
            &mut taut_flags,
            &mut taut_flags_done,
            z,
            0,
        )
    }

    #[test]
    fn source_port__inchi_dll_a2__normalization_step__line_1138() {
        let all_flags = u64::from(
            TG_FLAG_H_ALREADY_REMOVED
                | TG_FLAG_FIX_ISO_FIXEDH_BUG
                | TG_FLAG_FIX_TERM_H_CHRG_BUG
                | TG_FLAG_POINTED_EDGE_STEREO
                | TG_FLAG_PHOSPHINE_STEREO
                | TG_FLAG_ARSINE_STEREO
                | TG_FLAG_FIX_SP3_BUG,
        );
        let all_stereo_bits = (PES_BIT_POINT_EDGE_STEREO
            | PES_BIT_PHOSPHINE_STEREO
            | PES_BIT_ARSINE_STEREO
            | PES_BIT_FIX_SP3_BUG) as i32;

        let mut no_output_heap = SourceHeap::default();
        let no_output_input = no_output_heap
            .allocate_model_storage(normalization_atoms())
            .unwrap();
        let mut no_outputs = std::array::from_fn(|_| INP_ATOM_DATA::default());
        let mut no_output_z = COMPONENT_TREAT_INFO::default();
        assert_eq!(
            normalization_result(
                &mut no_output_heap,
                &mut no_outputs,
                no_output_input,
                u64::from(REQ_MODE_BASIC),
                all_flags,
                &mut no_output_z,
                2,
            ),
            Ok(-1)
        );
        assert!(no_output_z.at.iter().all(|pointer| pointer.is_null()));
        assert_eq!(no_output_z.fix_isofixedh, 1);
        assert_eq!(no_output_z.fix_termhchrg, 1);
        assert_eq!(no_output_z.bPointedEdgeStereo, all_stereo_bits);

        let mut null_input_heap = SourceHeap::default();
        let mut null_input_outputs =
            normalization_outputs(&mut null_input_heap, true, false, false, 2);
        let mut null_input_z = COMPONENT_TREAT_INFO::default();
        assert_eq!(
            normalization_result(
                &mut null_input_heap,
                &mut null_input_outputs,
                SourceMutPointer::null(),
                u64::from(REQ_MODE_BASIC),
                u64::from(TG_FLAG_H_ALREADY_REMOVED),
                &mut null_input_z,
                2,
            ),
            Ok(-1)
        );
        assert!(!null_input_z.at[TAUT_NON as usize].is_null());

        let mut negative_count_heap = SourceHeap::default();
        let negative_count_input = negative_count_heap
            .allocate_model_storage(normalization_atoms())
            .unwrap();
        let mut negative_count_outputs =
            normalization_outputs(&mut negative_count_heap, true, false, false, 2);
        let mut negative_count_z = COMPONENT_TREAT_INFO::default();
        assert_eq!(
            normalization_result(
                &mut negative_count_heap,
                &mut negative_count_outputs,
                negative_count_input,
                u64::from(REQ_MODE_BASIC),
                u64::from(TG_FLAG_H_ALREADY_REMOVED),
                &mut negative_count_z,
                -1,
            ),
            Ok(-1)
        );
        assert!(negative_count_z.at[TAUT_NON as usize].is_null());

        for allocations_before_failure in [0_u64, 1] {
            let mut heap = SourceHeap::default();
            let input = heap.allocate_model_storage(normalization_atoms()).unwrap();
            let mut outputs = normalization_outputs(&mut heap, true, true, false, 2);
            heap.fail_after_allocations(allocations_before_failure);
            let mut z = COMPONENT_TREAT_INFO::default();
            assert_eq!(
                normalization_result(
                    &mut heap,
                    &mut outputs,
                    input,
                    u64::from(REQ_MODE_BASIC | REQ_MODE_TAUT),
                    u64::from(TG_FLAG_H_ALREADY_REMOVED),
                    &mut z,
                    2,
                ),
                Ok(-1)
            );
            assert_eq!(heap.source_allocation_calls(), 2);
            assert_eq!(
                z.at[TAUT_NON as usize].is_null(),
                allocations_before_failure == 0
            );
            assert_eq!(
                z.at[TAUT_YES as usize].is_null(),
                allocations_before_failure == 1
            );
        }

        for (non_taut, taut, mode, expected_n1, expected_n2) in [
            (true, false, REQ_MODE_BASIC, TAUT_NON, TAUT_NON),
            (false, true, REQ_MODE_TAUT, TAUT_YES, TAUT_YES),
            (
                true,
                true,
                REQ_MODE_BASIC | REQ_MODE_TAUT,
                TAUT_YES,
                TAUT_YES,
            ),
        ] {
            let mut heap = SourceHeap::default();
            let atoms = normalization_atoms();
            let input = heap.allocate_model_storage(atoms.clone()).unwrap();
            let mut outputs = normalization_outputs(&mut heap, non_taut, taut, taut, 2);
            let mut z = COMPONENT_TREAT_INFO::default();
            let run_flags = if non_taut && !taut {
                all_flags & !u64::from(TG_FLAG_H_ALREADY_REMOVED)
            } else {
                all_flags
            };
            assert_eq!(
                normalization_result(
                    &mut heap,
                    &mut outputs,
                    input,
                    u64::from(mode),
                    run_flags,
                    &mut z,
                    2,
                ),
                Ok(2),
                "non_taut={non_taut} taut={taut}"
            );
            assert_eq!(z.num_atoms, 2);
            assert_eq!(z.num_deleted_H, 0);
            assert_eq!(z.n1, expected_n1 as i32);
            assert_eq!(z.n2, expected_n2 as i32);
            assert_eq!(z.fix_isofixedh, 1);
            assert_eq!(z.fix_termhchrg, 1);
            assert_eq!(z.bPointedEdgeStereo, all_stereo_bits);
            assert_ne!(z.nUserMode & u64::from(REQ_MODE_DEFAULT), 0);
            let mut normalized_atoms = atoms;
            normalized_atoms[0].nRingSystem = 2;
            normalized_atoms[0].nNumAtInRingSystem = 1;
            normalized_atoms[0].nBlockSystem = 1;
            normalized_atoms[1].nRingSystem = 1;
            normalized_atoms[1].nNumAtInRingSystem = 1;
            normalized_atoms[1].nBlockSystem = 1;
            for output in &outputs {
                if !output.at.is_null() {
                    let mut expected = normalized_atoms.clone();
                    if !taut {
                        expected[1].component = 1;
                    }
                    assert_eq!(heap.slice(output.at.as_const()).unwrap(), expected);
                }
            }
            if taut {
                assert_eq!(outputs[TAUT_YES as usize].num_at, 2);
                assert_eq!(outputs[TAUT_YES as usize].num_removed_H, 0);
                if !outputs[TAUT_YES as usize].at_fixed_bonds.is_null() {
                    assert_eq!(
                        heap.slice(outputs[TAUT_YES as usize].at_fixed_bonds.as_const())
                            .unwrap(),
                        normalized_atoms
                    );
                }
            }
        }

        let mut terminal_h_heap = SourceHeap::default();
        let mut terminal_h_atoms = vec![
            inp_ATOM {
                elname: [b'C' as i8, 0, 0, 0, 0, 0],
                el_number: 6,
                orig_at_number: 1,
                neighbor: {
                    let mut neighbors = [0; 20];
                    neighbors[0] = 1;
                    neighbors
                },
                bond_type: {
                    let mut bonds = [0; 20];
                    bonds[0] = BOND_SINGLE as u8;
                    bonds
                },
                valence: 1,
                chem_bonds_valence: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                elname: [b'H' as i8, 0, 0, 0, 0, 0],
                el_number: 1,
                orig_at_number: 2,
                neighbor: [0; 20],
                bond_type: {
                    let mut bonds = [0; 20];
                    bonds[0] = BOND_SINGLE as u8;
                    bonds
                },
                valence: 1,
                chem_bonds_valence: 1,
                ..inp_ATOM::default()
            },
        ];
        terminal_h_atoms[1].neighbor[0] = 0;
        let terminal_h_input = terminal_h_heap
            .allocate_model_storage(terminal_h_atoms)
            .unwrap();
        let mut terminal_h_outputs =
            normalization_outputs(&mut terminal_h_heap, true, false, false, 2);
        let mut terminal_h_z = COMPONENT_TREAT_INFO::default();
        assert_eq!(
            normalization_result(
                &mut terminal_h_heap,
                &mut terminal_h_outputs,
                terminal_h_input,
                u64::from(REQ_MODE_BASIC),
                u64::from(TG_FLAG_FIX_TERM_H_CHRG_BUG),
                &mut terminal_h_z,
                2,
            ),
            Ok(1)
        );
        assert_eq!(
            (
                terminal_h_z.num_atoms,
                terminal_h_z.num_deleted_H,
                terminal_h_outputs[TAUT_NON as usize].num_at,
                terminal_h_outputs[TAUT_NON as usize].num_removed_H,
            ),
            (1, 1, 2, 1)
        );
        let terminal_h_output = terminal_h_heap
            .slice(terminal_h_outputs[TAUT_NON as usize].at.as_const())
            .unwrap();
        assert_eq!(
            (
                terminal_h_output[0].valence,
                terminal_h_output[0].chem_bonds_valence,
                terminal_h_output[0].num_H,
                terminal_h_output[1].valence,
                terminal_h_output[1].neighbor[0],
            ),
            (0, 0, 1, 1, 0)
        );

        let mut proton_heap = SourceHeap::default();
        let proton_input = proton_heap
            .allocate_model_storage(vec![inp_ATOM {
                elname: [b'H' as i8, 0, 0, 0, 0, 0],
                el_number: 1,
                orig_at_number: 1,
                charge: 1,
                at_type: ATT_PROTON as AT_NUMB,
                ..inp_ATOM::default()
            }])
            .unwrap();
        let mut proton_outputs = normalization_outputs(&mut proton_heap, false, true, true, 1);
        let proton_fixed = proton_outputs[TAUT_YES as usize].at_fixed_bonds;
        let mut proton_z = COMPONENT_TREAT_INFO::default();
        assert_eq!(
            normalization_result(
                &mut proton_heap,
                &mut proton_outputs,
                proton_input,
                u64::from(REQ_MODE_TAUT),
                u64::from(TG_FLAG_H_ALREADY_REMOVED),
                &mut proton_z,
                1,
            ),
            Ok(1)
        );
        assert_eq!(proton_outputs[TAUT_YES as usize].bDeleted, 1);
        assert_eq!(proton_outputs[TAUT_YES as usize].nNumRemovedProtons, 1);
        assert_eq!(
            proton_outputs[TAUT_YES as usize].bNormalizationFlags,
            u64::from(FLAG_PROTON_SINGLE_REMOVED)
        );
        assert!(proton_outputs[TAUT_YES as usize].at_fixed_bonds.is_null());
        assert_eq!(
            proton_heap.slice(proton_fixed.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );

        let mut taut_heap = SourceHeap::default();
        let mut taut_atoms = vec![
            inp_ATOM {
                el_number: 8,
                orig_at_number: 1,
                chem_bonds_valence: 1,
                num_H: 1,
                num_iso_H: [1, 0, 0],
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 6,
                orig_at_number: 2,
                chem_bonds_valence: 3,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 8,
                orig_at_number: 3,
                chem_bonds_valence: 2,
                ..inp_ATOM::default()
            },
        ];
        taut_atoms[0].neighbor[0] = 1;
        taut_atoms[0].bond_type[0] = BOND_SINGLE as u8;
        taut_atoms[0].valence = 1;
        taut_atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
        taut_atoms[1].bond_type[..2].copy_from_slice(&[BOND_SINGLE as u8, BOND_DOUBLE as u8]);
        taut_atoms[1].valence = 2;
        taut_atoms[2].neighbor[0] = 1;
        taut_atoms[2].bond_type[0] = BOND_DOUBLE as u8;
        taut_atoms[2].valence = 1;
        let taut_input = taut_heap.allocate_model_storage(taut_atoms).unwrap();
        let mut taut_outputs = normalization_outputs(&mut taut_heap, false, true, true, 3);
        let mut taut_z = COMPONENT_TREAT_INFO::default();
        assert_eq!(
            normalization_result(
                &mut taut_heap,
                &mut taut_outputs,
                taut_input,
                u64::from(REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_NON_ISO | REQ_MODE_STEREO),
                u64::from(TG_FLAG_H_ALREADY_REMOVED | TG_FLAG_TEST_TAUT__ATOMS),
                &mut taut_z,
                3,
            ),
            Ok(3)
        );
        assert_eq!(taut_outputs[TAUT_YES as usize].bTautPreprocessed, 0);
        assert!(taut_outputs[TAUT_YES as usize].bTautomeric > 0);
        assert_ne!(
            taut_outputs[TAUT_YES as usize].bTautFlagsDone
                & u64::from(TG_FLAG_TEST_TAUT__ATOMS_DONE),
            0
        );
        assert_eq!(taut_z.bHasIsotopicAtoms, 1);
        assert_ne!(taut_z.nUserMode & u64::from(REQ_MODE_ISO), 0);
        assert_eq!(
            taut_z.nUserMode & u64::from(REQ_MODE_STEREO | REQ_MODE_ISO_STEREO),
            0
        );

        let mut preprocessed_heap = SourceHeap::default();
        let mut preprocessed_atoms = vec![
            inp_ATOM {
                el_number: 7,
                orig_at_number: 1,
                num_H: 4,
                charge: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 8,
                orig_at_number: 2,
                chem_bonds_valence: 2,
                num_H: 1,
                charge: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 6,
                orig_at_number: 3,
                chem_bonds_valence: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 6,
                orig_at_number: 4,
                chem_bonds_valence: 1,
                ..inp_ATOM::default()
            },
        ];
        preprocessed_atoms[1].neighbor[..2].copy_from_slice(&[2, 3]);
        preprocessed_atoms[1].bond_type[..2].fill(BOND_SINGLE as u8);
        preprocessed_atoms[1].valence = 2;
        for index in [2_usize, 3] {
            preprocessed_atoms[index].neighbor[0] = 1;
            preprocessed_atoms[index].bond_type[0] = BOND_SINGLE as u8;
            preprocessed_atoms[index].valence = 1;
        }
        let preprocessed_input = preprocessed_heap
            .allocate_model_storage(preprocessed_atoms)
            .unwrap();
        let mut preprocessed_outputs =
            normalization_outputs(&mut preprocessed_heap, false, true, true, 4);
        let preprocessed_fixed = preprocessed_outputs[TAUT_YES as usize].at_fixed_bonds;
        let mut preprocessed_z = COMPONENT_TREAT_INFO::default();
        assert_eq!(
            normalization_result(
                &mut preprocessed_heap,
                &mut preprocessed_outputs,
                preprocessed_input,
                u64::from(REQ_MODE_TAUT),
                u64::from(TG_FLAG_H_ALREADY_REMOVED | TG_FLAG_VARIABLE_PROTONS),
                &mut preprocessed_z,
                4,
            ),
            Ok(4)
        );
        assert_eq!(preprocessed_outputs[TAUT_YES as usize].bTautPreprocessed, 1);
        assert_ne!(
            preprocessed_outputs[TAUT_YES as usize].bNormalizationFlags
                & u64::from(FLAG_NORM_CONSIDER_TAUT),
            0
        );
        assert!(
            !preprocessed_outputs[TAUT_YES as usize]
                .at_fixed_bonds
                .is_null()
        );
        assert!(
            preprocessed_heap
                .slice(preprocessed_fixed.as_const())
                .is_ok()
        );

        let mut invalid_mode_heap = SourceHeap::default();
        let invalid_input = invalid_mode_heap
            .allocate_model_storage(normalization_atoms())
            .unwrap();
        let mut invalid_outputs =
            normalization_outputs(&mut invalid_mode_heap, true, false, false, 2);
        let mut invalid_z = COMPONENT_TREAT_INFO::default();
        assert_eq!(
            normalization_result(
                &mut invalid_mode_heap,
                &mut invalid_outputs,
                invalid_input,
                u64::from(REQ_MODE_ISO),
                u64::from(TG_FLAG_H_ALREADY_REMOVED),
                &mut invalid_z,
                2,
            ),
            Ok(-3)
        );
        assert_eq!(
            invalid_z.nUserMode & u64::from(REQ_MODE_BASIC | REQ_MODE_TAUT),
            0
        );
    }

    #[test]
    fn source_port__inchi_dll_a2__canonicalization_step__line_1612() {
        fn prepared_basic_case(
            heap: &mut SourceHeap,
        ) -> (
            [SourceMutPointer<INChI>; TAUT_NUM as usize],
            [SourceMutPointer<INChI_Aux>; TAUT_NUM as usize],
            [INP_ATOM_DATA; TAUT_NUM as usize],
            SourceMutPointer<INCHI_CLOCK>,
            COMPONENT_TREAT_INFO,
        ) {
            let mut atoms = normalization_atoms();
            for atom in &mut atoms {
                atom.elname[0] = b'C' as i8;
            }
            let input = heap.allocate_model_storage(atoms.clone()).unwrap();
            let mut outputs = normalization_outputs(heap, true, false, false, atoms.len());
            let mode = u64::from(REQ_MODE_BASIC | REQ_MODE_NON_ISO);
            let mut found_bonds = 0;
            let mut found_isotopes = 0;
            let inchi_non = Alloc_INChI(
                heap,
                &atoms,
                atoms.len() as i32,
                &mut found_bonds,
                &mut found_isotopes,
                mode as i32,
            )
            .unwrap();
            let aux_non =
                Alloc_INChI_Aux(heap, atoms.len() as i32, found_isotopes, mode as i32, 1).unwrap();
            let inchi = [inchi_non, SourceMutPointer::null()];
            let aux = [aux_non, SourceMutPointer::null()];
            let clock = heap
                .allocate_model_storage(vec![INCHI_CLOCK::default()])
                .unwrap();
            let mut z = COMPONENT_TREAT_INFO {
                nUserMode: mode,
                vABParityUnknown: AB_PARITY_UNDF as i32,
                ..COMPONENT_TREAT_INFO::default()
            };
            let mut taut_flags = u64::from(TG_FLAG_H_ALREADY_REMOVED);
            let mut taut_flags_done = 0;
            assert_eq!(
                Normalization_step(
                    heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &inchi,
                    &aux,
                    input,
                    &mut outputs,
                    atoms.len() as i32,
                    SourceMutPointer::null(),
                    &mut taut_flags,
                    &mut taut_flags_done,
                    &mut z,
                    0,
                ),
                Ok(atoms.len() as i32)
            );
            (inchi, aux, outputs, clock, z)
        }

        let mut heap = SourceHeap::default();
        let (inchi, aux, mut outputs, clock, mut z) = prepared_basic_case(&mut heap);
        outputs[TAUT_NON as usize].bDeleted = 7;
        outputs[TAUT_NON as usize].bTautFlags = 11;
        outputs[TAUT_NON as usize].bTautFlagsDone = 13;
        outputs[TAUT_NON as usize].bNormalizationFlags = 17;
        z.vt_group_info.num_iso_H = [3, 4, 5];
        z.vt_group_info.tni.nNumRemovedProtonsIsotopic = [6, 7, 8];
        z.vt_group_info.bTautFlagsDone =
            u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE);
        z.s[TAUT_NON as usize].bHasIsotopicTautGroups = 1;
        z.s[TAUT_NON as usize].nLenIsotopic = 19;
        z.s[TAUT_NON as usize].nLenIsotopicEndpoints = 23;
        z.s[TAUT_NON as usize].nLenLinearCTIsotopicTautomer = 29;
        z.s[TAUT_NON as usize].num_isotopic_atoms = 31;
        z.bHasIsotopicAtoms = 1;
        let live_before = heap.live_source_allocation_count();
        let types_before = heap.live_source_allocation_types();
        let atom_allocations_before = types_before
            .iter()
            .filter(|name| **name == std::any::type_name::<sp_ATOM>())
            .count();
        let mut errors = [0_i8; 256];
        assert_eq!(
            Canonicalization_step(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                clock,
                &inchi,
                &aux,
                &outputs,
                SourceMutPointer::null(),
                None,
                Some(&mut errors),
                &mut z,
                0,
                0,
            ),
            Ok(2)
        );
        assert!(z.at.iter().all(|pointer| pointer.is_null()));
        assert_eq!(
            z.vt_group_info,
            crate::source_types::T_GROUP_INFO::default()
        );
        assert_eq!(
            z.vt_group_info_orig,
            crate::source_types::T_GROUP_INFO::default()
        );
        assert_eq!(z.bHasIsotopicAtoms, 0);
        assert_eq!(
            (
                z.s[TAUT_NON as usize].bHasIsotopicTautGroups,
                z.s[TAUT_NON as usize].bIgnoreIsotopic,
                z.s[TAUT_NON as usize].nLenIsotopic,
                z.s[TAUT_NON as usize].nLenIsotopicEndpoints,
                z.s[TAUT_NON as usize].nLenLinearCTIsotopicTautomer,
                z.s[TAUT_NON as usize].num_isotopic_atoms,
            ),
            (0, 1, 0, 0, 0, 0)
        );
        let types_after = heap.live_source_allocation_types();
        let atom_allocations_after = types_after
            .iter()
            .filter(|name| **name == std::any::type_name::<sp_ATOM>())
            .count();
        assert_eq!(atom_allocations_before, atom_allocations_after + 1);
        assert_eq!(heap.live_source_allocation_count(), live_before + 1);
        let generated = &heap.slice(inchi[TAUT_NON as usize].as_const()).unwrap()[0];
        assert_eq!(
            (
                generated.nErrorCode,
                generated.nNumberOfAtoms,
                generated.nTotalCharge,
                generated.bDeleted,
            ),
            (0, 2, 0, 7)
        );
        assert_eq!(
            heap.slice(generated.nAtom.as_const()).unwrap(),
            &[6_u8, 6_u8]
        );
        let generated_aux = &heap.slice(aux[TAUT_NON as usize].as_const()).unwrap()[0];
        assert_eq!(
            (
                generated_aux.nErrorCode,
                generated_aux.nNumberOfAtoms,
                generated_aux.bDeleted,
                generated_aux.bTautFlags,
                generated_aux.bTautFlagsDone,
                generated_aux.bNormalizationFlags,
            ),
            (0, 2, 7, 11, 13, 17)
        );

        let mut transfer_heap = SourceHeap::default();
        let (transfer_inchi, transfer_aux, transfer_outputs, transfer_clock, mut transfer_z) =
            prepared_basic_case(&mut transfer_heap);
        let transferred_endpoint =
            inchi_calloc::<AT_NUMB>(&mut transfer_heap, 1, std::mem::size_of::<AT_NUMB>() as u64)
                .unwrap();
        transfer_heap.slice_mut(transferred_endpoint).unwrap()[0] = 1;
        transfer_z.vt_group_info.nEndpointAtomNumber = transferred_endpoint;
        let mut transferred = crate::source_types::T_GROUP_INFO::default();
        assert_eq!(
            Canonicalization_step(
                &mut transfer_heap,
                &mut CANON_GLOBALS::default(),
                transfer_clock,
                &transfer_inchi,
                &transfer_aux,
                &transfer_outputs,
                SourceMutPointer::null(),
                Some(&mut transferred),
                None,
                &mut transfer_z,
                0,
                0,
            ),
            Ok(2)
        );
        assert_eq!(
            transferred.nEndpointAtomNumber,
            transfer_z.vt_group_info.nEndpointAtomNumber
        );
        assert_eq!(
            transfer_heap
                .slice(transferred.nEndpointAtomNumber.as_const())
                .unwrap(),
            &[1]
        );
        inchi_free(&mut transfer_heap, transferred.nEndpointAtomNumber).unwrap();

        let mut failure_heap = SourceHeap::default();
        let (failure_inchi, failure_aux, failure_outputs, failure_clock, mut failure_z) =
            prepared_basic_case(&mut failure_heap);
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            Canonicalization_step(
                &mut failure_heap,
                &mut CANON_GLOBALS::default(),
                failure_clock,
                &failure_inchi,
                &failure_aux,
                &failure_outputs,
                SourceMutPointer::null(),
                None,
                None,
                &mut failure_z,
                0,
                0,
            ),
            Ok(CT_OUT_OF_RAM)
        );
        assert!(failure_z.at.iter().all(|pointer| pointer.is_null()));
        assert_eq!(
            failure_z.vt_group_info,
            crate::source_types::T_GROUP_INFO::default()
        );
        assert_eq!(
            failure_z.vt_group_info_orig,
            crate::source_types::T_GROUP_INFO::default()
        );
    }

    fn norm_component_fixture(
        heap: &mut SourceHeap,
        atoms: Vec<inp_ATOM>,
        mode: u64,
        output_options: i32,
        taut_flags_done: u64,
        existing_inchi: [SourceMutPointer<INChI>; TAUT_NUM as usize],
        existing_aux: [SourceMutPointer<INChI_Aux>; TAUT_NUM as usize],
    ) -> (INCHIGEN_CONTROL, SourceMutPointer<INCHI_CLOCK>) {
        let atom_count = i32::try_from(atoms.len()).unwrap();
        let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
        let mut control = INCHIGEN_CONTROL::default();
        control.InpParms.nMode = mode;
        control.InpParms.bINChIOutputOptions = output_options;
        control.InpParms.bTautFlags = u64::from(TG_FLAG_H_ALREADY_REMOVED);
        control.InpParms.bTautFlagsDone = taut_flags_done;
        control.InpParms.bLooseTSACheck = 17;
        control.InpParms.bStereoAtZz = 19;
        control.InpCurAtData[INCHI_BAS as usize] = heap
            .allocate_model_storage(vec![INP_ATOM_DATA {
                at: atom_pointer,
                num_at: atom_count,
                ..INP_ATOM_DATA::default()
            }])
            .unwrap();
        control.InpNormAtData[INCHI_BAS as usize] = heap
            .allocate_model_storage(vec![INP_ATOM_DATA::default()])
            .unwrap();
        control.InpNormTautData[INCHI_BAS as usize] = heap
            .allocate_model_storage(vec![INP_ATOM_DATA::default()])
            .unwrap();
        control.cti[INCHI_BAS as usize] = heap
            .allocate_model_storage(vec![COMPONENT_TREAT_INFO::default()])
            .unwrap();
        control.pINChI[INCHI_BAS as usize] =
            heap.allocate_model_storage(vec![existing_inchi]).unwrap();
        control.pINChI_Aux[INCHI_BAS as usize] =
            heap.allocate_model_storage(vec![existing_aux]).unwrap();
        let clock = heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        (control, clock)
    }

    fn norm_component_atom() -> inp_ATOM {
        inp_ATOM {
            elname: [b'C' as i8, 0, 0, 0, 0, 0],
            el_number: 6,
            orig_at_number: 1,
            iso_atw_diff: 1,
            ..inp_ATOM::default()
        }
    }

    fn norm_structure_fixture(
        heap: &mut SourceHeap,
        component_count: usize,
    ) -> (
        INCHIGEN_DATA,
        INCHIGEN_CONTROL,
        SourceMutPointer<INCHI_CLOCK>,
    ) {
        let atoms = (0..component_count)
            .map(|component| inp_ATOM {
                elname: [b'C' as i8, 0, 0, 0, 0, 0],
                el_number: 6,
                orig_at_number: u16::try_from(component + 1).unwrap(),
                component: u16::try_from(component + 1).unwrap(),
                ..inp_ATOM::default()
            })
            .collect::<Vec<_>>();
        let atom_pointer = heap.allocate_model_storage(atoms).unwrap();
        let lengths = heap
            .allocate_model_storage(vec![1_u16; component_count])
            .unwrap();
        let atom_count = i32::try_from(component_count).unwrap();
        let mut control = INCHIGEN_CONTROL::default();
        control.num_inp = 41;
        control.InpParms.nMode = u64::from(REQ_MODE_BASIC);
        control.InpParms.bTautFlags = u64::from(TG_FLAG_H_ALREADY_REMOVED);
        control.OrigInpData = crate::source_types::ORIG_ATOM_DATA {
            at: atom_pointer,
            num_inp_atoms: atom_count,
            num_components: atom_count,
            ..crate::source_types::ORIG_ATOM_DATA::default()
        };
        control.PrepInpData[INCHI_BAS as usize] = crate::source_types::ORIG_ATOM_DATA {
            at: atom_pointer,
            num_inp_atoms: atom_count,
            num_components: atom_count,
            nCurAtLen: lengths,
            ..crate::source_types::ORIG_ATOM_DATA::default()
        };
        let clock = heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        (INCHIGEN_DATA::default(), control, clock)
    }

    #[test]
    fn source_port__inchi_dll_a2__normonecomponentinchi__line_654() {
        let null_inchi = [SourceMutPointer::null(); TAUT_NUM as usize];
        let null_aux = [SourceMutPointer::null(); TAUT_NUM as usize];

        let mut basic_heap = SourceHeap::default();
        let (mut basic, basic_clock) = norm_component_fixture(
            &mut basic_heap,
            vec![norm_component_atom()],
            u64::from(REQ_MODE_BASIC | REQ_MODE_ISO | REQ_MODE_DIFF_UU_STEREO),
            0,
            0,
            null_inchi,
            null_aux,
        );
        basic.InpParms.msec_MaxTime = 1;
        basic.InpParms.msec_LeftTime = 37;
        basic.StructData.ulStructTime = u64::MAX;
        assert_eq!(
            NormOneComponentINChI(
                &mut basic_heap,
                &mut CANON_GLOBALS::default(),
                basic_clock,
                &mut basic,
                INCHI_BAS as i32,
                0,
                125_000,
            ),
            Ok(0)
        );
        assert_eq!(basic.StructData.nErrorCode, 0);
        assert_eq!(basic.InpParms.msec_LeftTime, 37);
        assert_eq!(basic.StructData.ulStructTime, u64::MAX);
        let basic_row = basic_heap
            .slice(basic.pINChI[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        let basic_aux_row = basic_heap
            .slice(basic.pINChI_Aux[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        assert!(!basic_row[TAUT_NON as usize].is_null());
        assert!(basic_row[TAUT_YES as usize].is_null());
        let basic_inchi = &basic_heap
            .slice(basic_row[TAUT_NON as usize].as_const())
            .unwrap()[0];
        assert!(basic_inchi.StereoIsotopic.is_null());
        assert!(basic_inchi.nPossibleLocationsOfIsotopicH.is_null());
        let basic_aux_value = &basic_heap
            .slice(basic_aux_row[TAUT_NON as usize].as_const())
            .unwrap()[0];
        assert_eq!(basic_aux_value.bIsIsotopic, 1);
        assert!(!basic_aux_value.szOrigCoord.is_null());
        assert!(basic_aux_value.nIsotopicOrigAtNosInCanonOrd.is_null());
        let current = &basic_heap
            .slice(basic.InpCurAtData[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        assert_eq!((current.num_bonds, current.num_isotopic), (0, 1));
        assert_eq!(
            basic_heap.slice(current.at.as_const()).unwrap()[0].component,
            1
        );
        let basic_cti = &basic_heap
            .slice(basic.cti[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        assert_eq!(basic_cti.nUserMode, basic.InpParms.nMode);
        assert_eq!(basic_cti.bLooseTSACheck, 17);
        assert_eq!(basic_cti.bStereoAtZz, 19);
        assert_eq!(basic_cti.vABParityUnknown, AB_PARITY_UNKN as i32);
        assert_ne!(
            basic_heap.slice(basic_clock.as_const()).unwrap()[0].m_MaxPositiveClock,
            0
        );
        let basic_norm = &basic_heap
            .slice(basic.InpNormAtData[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        let basic_taut_norm = &basic_heap
            .slice(basic.InpNormTautData[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        assert!(!basic_norm.at.is_null());
        assert!(basic_taut_norm.at.is_null());

        let mut old_heap = SourceHeap::default();
        let old_basic = old_heap
            .allocate_model_storage(vec![INChI::default()])
            .unwrap();
        let old_basic_aux = old_heap
            .allocate_model_storage(vec![INChI_Aux::default()])
            .unwrap();
        let (mut taut, taut_clock) = norm_component_fixture(
            &mut old_heap,
            vec![norm_component_atom()],
            u64::from(REQ_MODE_TAUT | REQ_MODE_ISO),
            INCHI_OUT_NO_AUX_INFO as i32,
            0,
            [old_basic, SourceMutPointer::null()],
            [old_basic_aux, SourceMutPointer::null()],
        );
        assert_eq!(
            NormOneComponentINChI(
                &mut old_heap,
                &mut CANON_GLOBALS::default(),
                taut_clock,
                &mut taut,
                INCHI_BAS as i32,
                0,
                -1,
            ),
            Ok(0)
        );
        let taut_row = old_heap
            .slice(taut.pINChI[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        let taut_aux_row = old_heap
            .slice(taut.pINChI_Aux[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        assert_eq!(taut_row[TAUT_NON as usize], old_basic);
        assert_eq!(taut_aux_row[TAUT_NON as usize], old_basic_aux);
        assert!(!taut_row[TAUT_YES as usize].is_null());
        let taut_inchi = &old_heap
            .slice(taut_row[TAUT_YES as usize].as_const())
            .unwrap()[0];
        assert!(!taut_inchi.StereoIsotopic.is_null());
        assert!(!taut_inchi.nPossibleLocationsOfIsotopicH.is_null());
        let taut_aux_value = &old_heap
            .slice(taut_aux_row[TAUT_YES as usize].as_const())
            .unwrap()[0];
        assert!(!taut_aux_value.nIsotopicOrigAtNosInCanonOrd.is_null());
        assert!(taut_aux_value.szOrigCoord.is_null());
        assert_eq!(
            old_heap
                .slice(taut.cti[INCHI_BAS as usize].as_const())
                .unwrap()[0]
                .vABParityUnknown,
            AB_PARITY_UNDF as i32
        );

        for (flags_done, expect_basic_iso) in [
            (0_u64, false),
            (u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE), true),
            (u64::from(TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE), true),
        ] {
            let mut heap = SourceHeap::default();
            let (mut both, clock) = norm_component_fixture(
                &mut heap,
                vec![norm_component_atom()],
                u64::from(REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO),
                INCHI_OUT_SHORT_AUX_INFO as i32,
                flags_done,
                null_inchi,
                null_aux,
            );
            assert_eq!(
                NormOneComponentINChI(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &mut both,
                    INCHI_BAS as i32,
                    0,
                    7,
                ),
                Ok(0)
            );
            let row = heap
                .slice(both.pINChI[INCHI_BAS as usize].as_const())
                .unwrap()[0];
            let aux_row = heap
                .slice(both.pINChI_Aux[INCHI_BAS as usize].as_const())
                .unwrap()[0];
            assert_eq!(
                !heap.slice(row[TAUT_NON as usize].as_const()).unwrap()[0]
                    .StereoIsotopic
                    .is_null(),
                expect_basic_iso
            );
            assert!(
                !heap.slice(row[TAUT_YES as usize].as_const()).unwrap()[0]
                    .StereoIsotopic
                    .is_null()
            );
            assert!(
                heap.slice(aux_row[TAUT_NON as usize].as_const()).unwrap()[0]
                    .szOrigCoord
                    .is_null()
            );
            assert!(
                heap.slice(aux_row[TAUT_YES as usize].as_const()).unwrap()[0]
                    .szOrigCoord
                    .is_null()
            );
        }

        for (old_error, expected_return) in [
            (91, _IS_ERROR as i32),
            (CT_OUT_OF_RAM, _IS_FATAL as i32),
            (CT_USER_QUIT_ERR, _IS_FATAL as i32),
        ] {
            let mut heap = SourceHeap::default();
            let old_taut = heap
                .allocate_model_storage(vec![INChI {
                    nErrorCode: old_error,
                    ..INChI::default()
                }])
                .unwrap();
            let (mut control, clock) = norm_component_fixture(
                &mut heap,
                vec![norm_component_atom()],
                u64::from(REQ_MODE_BASIC),
                0,
                0,
                [SourceMutPointer::null(), old_taut],
                null_aux,
            );
            assert_eq!(
                NormOneComponentINChI(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &mut control,
                    INCHI_BAS as i32,
                    0,
                    11,
                ),
                Ok(expected_return)
            );
            assert_eq!(control.StructData.nErrorCode, old_error);
            assert_eq!(
                heap.slice(control.pINChI[INCHI_BAS as usize].as_const())
                    .unwrap()[0][TAUT_YES as usize],
                old_taut
            );
        }

        let mut flags_heap = SourceHeap::default();
        let old_taut = flags_heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                lenTautomer: 0,
                ..INChI::default()
            }])
            .unwrap();
        let old_taut_aux = flags_heap
            .allocate_model_storage(vec![INChI_Aux {
                nNumberOfAtoms: 1,
                nNumRemovedIsotopicH: [1, 0, 0],
                bNormalizationFlags: 1 << 8,
                bTautFlags: 1 << 9,
                bTautFlagsDone: 1 << 10,
                nCanonFlags: 1 << 11,
                ..INChI_Aux::default()
            }])
            .unwrap();
        let (mut flags, flags_clock) = norm_component_fixture(
            &mut flags_heap,
            vec![norm_component_atom()],
            u64::from(REQ_MODE_BASIC),
            0,
            0,
            [SourceMutPointer::null(), old_taut],
            [SourceMutPointer::null(), old_taut_aux],
        );
        assert_eq!(
            NormOneComponentINChI(
                &mut flags_heap,
                &mut CANON_GLOBALS::default(),
                flags_clock,
                &mut flags,
                INCHI_BAS as i32,
                0,
                3,
            ),
            Ok(0)
        );
        assert_eq!(
            flags.ncFlags.bNormalizationFlags[INCHI_BAS as usize][TAUT_YES as usize],
            1 << 8
        );
        assert_eq!(
            flags.ncFlags.bTautFlags[INCHI_BAS as usize][TAUT_YES as usize],
            1 << 9
        );
        assert_eq!(
            flags.ncFlags.bTautFlagsDone[INCHI_BAS as usize][TAUT_YES as usize],
            1 << 10
        );
        assert_eq!(
            flags.ncFlags.nCanonFlags[INCHI_BAS as usize][TAUT_YES as usize],
            1 << 11
        );
        assert_eq!(flags.StructData.num_non_taut[INCHI_BAS as usize], 1);
        assert_eq!(flags.StructData.num_taut[INCHI_BAS as usize], 0);
        let retained_norm = &flags_heap
            .slice(flags.InpNormTautData[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        assert_eq!(
            (retained_norm.bExists, retained_norm.bHasIsotopicLayer),
            (1, 1)
        );

        let mut observed_inchi_failure = false;
        let mut observed_aux_failure = false;
        let mut observed_output_failure = false;
        for failure_index in 0_u64..96 {
            let mut heap = SourceHeap::default();
            let (mut control, clock) = norm_component_fixture(
                &mut heap,
                vec![norm_component_atom()],
                u64::from(REQ_MODE_BASIC),
                0,
                0,
                null_inchi,
                null_aux,
            );
            heap.fail_after_allocations(failure_index);
            let _ = NormOneComponentINChI(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                clock,
                &mut control,
                INCHI_BAS as i32,
                0,
                5,
            )
            .unwrap();
            let row = heap
                .slice(control.pINChI[INCHI_BAS as usize].as_const())
                .unwrap()[0];
            let aux_row = heap
                .slice(control.pINChI_Aux[INCHI_BAS as usize].as_const())
                .unwrap()[0];
            let output = &heap
                .slice(control.InpNormAtData[INCHI_BAS as usize].as_const())
                .unwrap()[0];
            observed_inchi_failure |= row[TAUT_NON as usize].is_null();
            observed_aux_failure |=
                !row[TAUT_NON as usize].is_null() && aux_row[TAUT_NON as usize].is_null();
            observed_output_failure |= !row[TAUT_NON as usize].is_null()
                && !aux_row[TAUT_NON as usize].is_null()
                && output.at.is_null();
            if observed_inchi_failure && observed_aux_failure && observed_output_failure {
                break;
            }
        }
        assert!(observed_inchi_failure);
        assert!(observed_aux_failure);
        assert!(observed_output_failure);

        for (count, expected_error) in [(0_i32, -3_i32), (-1_i32, -1_i32)] {
            let mut heap = SourceHeap::default();
            let (mut control, clock) = norm_component_fixture(
                &mut heap,
                Vec::new(),
                u64::from(REQ_MODE_BASIC),
                0,
                0,
                null_inchi,
                null_aux,
            );
            heap.slice_mut(control.InpCurAtData[INCHI_BAS as usize])
                .unwrap()[0]
                .num_at = count;
            assert_eq!(
                NormOneComponentINChI(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &mut control,
                    INCHI_BAS as i32,
                    0,
                    0,
                ),
                Ok(_IS_ERROR as i32)
            );
            assert_eq!(control.StructData.nErrorCode, expected_error);
        }

        for left_time in [-9_i64, 0, 9] {
            let mut heap = SourceHeap::default();
            let (mut control, clock) = norm_component_fixture(
                &mut heap,
                vec![norm_component_atom()],
                u64::from(REQ_MODE_BASIC),
                0,
                0,
                null_inchi,
                null_aux,
            );
            control.InpParms.msec_MaxTime = 1;
            control.InpParms.msec_LeftTime = left_time;
            assert_eq!(
                NormOneComponentINChI(
                    &mut heap,
                    &mut CANON_GLOBALS::default(),
                    clock,
                    &mut control,
                    INCHI_BAS as i32,
                    0,
                    i64::MAX,
                ),
                Ok(0)
            );
            assert_eq!(control.InpParms.msec_LeftTime, left_time);
        }

        let mut invalid_heap = SourceHeap::default();
        let invalid_clock = invalid_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut invalid = INCHIGEN_CONTROL::default();
        assert_eq!(
            NormOneComponentINChI(
                &mut invalid_heap,
                &mut CANON_GLOBALS::default(),
                invalid_clock,
                &mut invalid,
                -1,
                0,
                0,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            NormOneComponentINChI(
                &mut invalid_heap,
                &mut CANON_GLOBALS::default(),
                invalid_clock,
                &mut invalid,
                INCHI_BAS as i32,
                i32::MAX,
                0,
            ),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__inchi_dll_a2__canononecomponentinchi__line_923() {
        let null_inchi = [SourceMutPointer::<INChI>::null(); TAUT_NUM as usize];
        let null_aux = [SourceMutPointer::<INChI_Aux>::null(); TAUT_NUM as usize];
        let mut heap = SourceHeap::default();
        let (mut control, clock) = norm_component_fixture(
            &mut heap,
            vec![norm_component_atom()],
            u64::from(REQ_MODE_BASIC | REQ_MODE_NON_ISO),
            0,
            0,
            null_inchi,
            null_aux,
        );
        control.InpParms.msec_MaxTime = 1;
        control.InpParms.msec_LeftTime = 23;
        control.StructData.ulStructTime = u64::MAX;
        assert_eq!(
            NormOneComponentINChI(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                clock,
                &mut control,
                INCHI_BAS as i32,
                0,
                125_000,
            ),
            Ok(0)
        );
        assert_eq!(
            CanonOneComponentINChI(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                clock,
                &mut control,
                INCHI_BAS as i32,
                0,
                125_000,
            ),
            Ok(0)
        );
        assert_eq!(control.StructData.nErrorCode, 0);
        assert_eq!(control.InpParms.msec_LeftTime, 23);
        assert_eq!(control.StructData.ulStructTime, u64::MAX);
        assert_eq!(control.StructData.num_non_taut[INCHI_BAS as usize], 1);
        assert_eq!(control.StructData.num_taut[INCHI_BAS as usize], 0);
        let row = heap
            .slice(control.pINChI[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        let generated = &heap.slice(row[TAUT_NON as usize].as_const()).unwrap()[0];
        assert_eq!(
            (
                generated.nErrorCode,
                generated.nNumberOfAtoms,
                generated.nTotalCharge,
                generated.lenTautomer,
            ),
            (0, 1, 0, 0)
        );
        assert_eq!(heap.slice(generated.nAtom.as_const()).unwrap(), &[6]);
        assert_eq!(
            &heap.slice(generated.szHillFormula.as_const()).unwrap()[..2],
            &[b'C' as i8, 0]
        );
        let normalized = &heap
            .slice(control.InpNormAtData[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        assert_eq!((normalized.bExists, normalized.bHasIsotopicLayer), (1, 0));
        let cti = &heap
            .slice(control.cti[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        assert!(cti.at.iter().all(|pointer| pointer.is_null()));
        let current = &heap
            .slice(control.InpCurAtData[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        assert_eq!(heap.slice(current.at.as_const()).unwrap()[0].component, 1);

        let mut failure_heap = SourceHeap::default();
        let (mut failure, failure_clock) = norm_component_fixture(
            &mut failure_heap,
            vec![norm_component_atom()],
            u64::from(REQ_MODE_BASIC | REQ_MODE_NON_ISO),
            0,
            0,
            null_inchi,
            null_aux,
        );
        assert_eq!(
            NormOneComponentINChI(
                &mut failure_heap,
                &mut CANON_GLOBALS::default(),
                failure_clock,
                &mut failure,
                INCHI_BAS as i32,
                0,
                0,
            ),
            Ok(0)
        );
        failure_heap.fail_after_allocations(0);
        assert_eq!(
            CanonOneComponentINChI(
                &mut failure_heap,
                &mut CANON_GLOBALS::default(),
                failure_clock,
                &mut failure,
                INCHI_BAS as i32,
                0,
                0,
            ),
            Ok(_IS_FATAL as i32)
        );
        assert_eq!(failure.StructData.nErrorCode, CT_OUT_OF_RAM);
        assert!(
            failure_heap
                .slice(failure.cti[INCHI_BAS as usize].as_const())
                .unwrap()[0]
                .at
                .iter()
                .all(|pointer| pointer.is_null())
        );
    }

    #[test]
    fn source_port__inchi_dll_a2__canononestructureinchi__line_503() {
        let mut heap = SourceHeap::default();
        let (mut generation, mut control, clock) = norm_structure_fixture(&mut heap, 1);
        control.InpParms.nMode = u64::from(REQ_MODE_BASIC | REQ_MODE_NON_ISO);
        control.InpParms.bTautFlags = u64::from(TG_FLAG_H_ALREADY_REMOVED);
        control.InpParms.msec_MaxTime = 1;
        control.InpParms.msec_LeftTime = 31;
        assert_eq!(
            NormOneStructureINChI(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                clock,
                &mut generation,
                &mut control,
                INCHI_BAS as i32,
                None,
                125_000,
            ),
            Ok(0)
        );
        assert_eq!(
            CanonOneStructureINChI(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                clock,
                &mut control,
                INCHI_BAS as i32,
                None,
                125_000,
            ),
            Ok(0)
        );
        assert_eq!(control.StructData.nErrorCode, 0);
        assert_eq!(control.StructData.num_non_taut[INCHI_BAS as usize], 1);
        let row = heap
            .slice(control.pINChI[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        let generated = &heap.slice(row[TAUT_NON as usize].as_const()).unwrap()[0];
        assert_eq!((generated.nErrorCode, generated.nNumberOfAtoms), (0, 1));
        assert_eq!(heap.slice(generated.nAtom.as_const()).unwrap(), &[6]);
        assert_eq!(
            heap.slice(control.InpCurAtData[INCHI_BAS as usize].as_const())
                .unwrap()[0]
                .num_at,
            1
        );
        let normalized = &heap
            .slice(control.InpNormAtData[INCHI_BAS as usize].as_const())
            .unwrap()[0];
        assert_eq!(normalized, &INP_ATOM_DATA::default());

        let mut quit_heap = SourceHeap::default();
        let (_quit_generation, mut quit_control, quit_clock) =
            norm_structure_fixture(&mut quit_heap, 1);
        quit_control.StructData.bUserQuitComponent = 1;
        assert_eq!(
            CanonOneStructureINChI(
                &mut quit_heap,
                &mut CANON_GLOBALS::default(),
                quit_clock,
                &mut quit_control,
                INCHI_BAS as i32,
                None,
                0,
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__inchi_dll_a2__normonestructureinchi__line_188() {
        let mut heap = SourceHeap::default();
        let (mut generation, mut control, clock) = norm_structure_fixture(&mut heap, 2);
        control.InpParms.msec_MaxTime = 29;
        assert_eq!(
            NormOneStructureINChI(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                clock,
                &mut generation,
                &mut control,
                INCHI_BAS as i32,
                None,
                73,
            ),
            Ok(0)
        );
        assert_eq!(control.InpParms.msec_LeftTime, 29);
        assert_eq!(control.StructData.num_components[INCHI_BAS as usize], 2);
        assert_eq!(
            heap.slice(control.pINChI[INCHI_BAS as usize].as_const())
                .unwrap()
                .len(),
            3
        );
        assert_eq!(
            heap.slice(control.pINChI_Aux[INCHI_BAS as usize].as_const())
                .unwrap()
                .len(),
            3
        );
        assert!(!generation.NormAtomsNontaut[INCHI_BAS as usize].is_null());
        assert!(!generation.NormAtomsTaut[INCHI_BAS as usize].is_null());
        assert_eq!(
            heap.slice(control.InpCurAtData[INCHI_BAS as usize].as_const())
                .unwrap()
                .len(),
            2
        );
        for component in 0..2_usize {
            let rows = heap
                .slice(control.pINChI[INCHI_BAS as usize].as_const())
                .unwrap();
            assert!(!rows[component][TAUT_NON as usize].is_null());
            assert!(rows[component][TAUT_YES as usize].is_null());
            let current = &heap
                .slice(control.InpCurAtData[INCHI_BAS as usize].as_const())
                .unwrap()[component];
            assert_eq!(current.num_at, 1);
            assert_eq!(
                heap.slice(current.at.as_const()).unwrap()[0].component,
                u16::try_from(component + 1).unwrap()
            );
        }
        assert_eq!(
            control.composite_norm_data[INCHI_BAS as usize],
            std::array::from_fn(|_| COMP_ATOM_DATA::default())
        );
        assert_eq!(heap.live_source_allocations_of::<INP_ATOM_DATA2>(), 0);

        let mut empty_heap = SourceHeap::default();
        let empty_clock = empty_heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut empty_generation = INCHIGEN_DATA::default();
        let mut empty = INCHIGEN_CONTROL::default();
        empty.InpParms.msec_MaxTime = 17;
        empty.composite_norm_data[INCHI_BAS as usize][TAUT_NON as usize].num_at = 91;
        assert_eq!(
            NormOneStructureINChI(
                &mut empty_heap,
                &mut CANON_GLOBALS::default(),
                empty_clock,
                &mut empty_generation,
                &mut empty,
                INCHI_BAS as i32,
                None,
                0,
            ),
            Ok(0)
        );
        assert_eq!(empty.InpParms.msec_LeftTime, 17);
        assert_eq!(
            empty.composite_norm_data[INCHI_BAS as usize],
            std::array::from_fn(|_| COMP_ATOM_DATA::default())
        );
        empty.InpParms.bAllowEmptyStructure = 1;
        assert_eq!(
            NormOneStructureINChI(
                &mut empty_heap,
                &mut CANON_GLOBALS::default(),
                empty_clock,
                &mut empty_generation,
                &mut empty,
                INCHI_BAS as i32,
                None,
                0,
            ),
            Ok(0)
        );
        assert!(!empty.pINChI[INCHI_BAS as usize].is_null());
        assert!(!empty.pINChI_Aux[INCHI_BAS as usize].is_null());

        let mut reconnect_heap = SourceHeap::default();
        let (mut reconnect_generation, mut reconnect, reconnect_clock) =
            norm_structure_fixture(&mut reconnect_heap, 1);
        assert_eq!(
            NormOneStructureINChI(
                &mut reconnect_heap,
                &mut CANON_GLOBALS::default(),
                reconnect_clock,
                &mut reconnect_generation,
                &mut reconnect,
                INCHI_REC as i32,
                None,
                0,
            ),
            Ok(0)
        );
        assert!(reconnect.pINChI[INCHI_REC as usize].is_null());

        let mut invalid_heap = SourceHeap::default();
        let (mut invalid_generation, mut invalid, invalid_clock) =
            norm_structure_fixture(&mut invalid_heap, 1);
        assert_eq!(
            NormOneStructureINChI(
                &mut invalid_heap,
                &mut CANON_GLOBALS::default(),
                invalid_clock,
                &mut invalid_generation,
                &mut invalid,
                2,
                None,
                0,
            ),
            Ok(_IS_FATAL as i32)
        );
        assert_eq!(invalid.StructData.nStructReadError, 97);
        assert_eq!(invalid.StructData.nErrorType, _IS_FATAL as i32);
        let invalid_message = invalid
            .StructData
            .pStrErrStruct
            .iter()
            .take_while(|byte| **byte != 0)
            .map(|byte| *byte as u8)
            .collect::<Vec<_>>();
        assert_eq!(invalid_message, b"Fatal undetermined program error");

        for failure_index in [0_u64, 1] {
            let mut failure_heap = SourceHeap::default();
            let (mut failure_generation, mut failure, failure_clock) =
                norm_structure_fixture(&mut failure_heap, 1);
            let retained_inchi = failure_heap
                .allocate_model_storage(vec![[SourceMutPointer::null(); TAUT_NUM as usize]])
                .unwrap();
            let retained_aux = failure_heap
                .allocate_model_storage(vec![[SourceMutPointer::null(); TAUT_NUM as usize]])
                .unwrap();
            failure.pINChI[INCHI_BAS as usize] = retained_inchi;
            failure.pINChI_Aux[INCHI_BAS as usize] = retained_aux;
            failure_heap.fail_after_allocations(failure_index);
            assert_eq!(
                NormOneStructureINChI(
                    &mut failure_heap,
                    &mut CANON_GLOBALS::default(),
                    failure_clock,
                    &mut failure_generation,
                    &mut failure,
                    INCHI_BAS as i32,
                    None,
                    0,
                ),
                Ok(0)
            );
            assert_eq!(failure.StructData.nStructReadError, 99);
            assert_eq!(failure.StructData.nErrorType, _IS_FATAL as i32);
            assert_eq!(failure.pINChI[INCHI_BAS as usize], retained_inchi);
            assert_eq!(failure.pINChI_Aux[INCHI_BAS as usize], retained_aux);
            assert_eq!(failure_heap.source_allocation_calls(), 2);
        }

        let mut visible_heap = SourceHeap::default();
        let (mut visible_generation, mut visible, visible_clock) =
            norm_structure_fixture(&mut visible_heap, 1);
        visible_heap.fail_after_allocations(2);
        assert_eq!(
            NormOneStructureINChI(
                &mut visible_heap,
                &mut CANON_GLOBALS::default(),
                visible_clock,
                &mut visible_generation,
                &mut visible,
                INCHI_BAS as i32,
                None,
                0,
            ),
            Ok(0)
        );
        assert!(visible_generation.NormAtomsNontaut[INCHI_BAS as usize].is_null());
        assert!(!visible_generation.NormAtomsTaut[INCHI_BAS as usize].is_null());

        let mut intermediate_heap = SourceHeap::default();
        let (mut intermediate_generation, mut intermediate, intermediate_clock) =
            norm_structure_fixture(&mut intermediate_heap, 2);
        intermediate_heap.fail_after_allocations(0);
        assert_eq!(
            NormOneStructureINChI(
                &mut intermediate_heap,
                &mut CANON_GLOBALS::default(),
                intermediate_clock,
                &mut intermediate_generation,
                &mut intermediate,
                INCHI_BAS as i32,
                None,
                0,
            ),
            Ok(0)
        );
        assert_eq!(intermediate.StructData.nErrorType, 0);

        let mut preprocessing_heap = SourceHeap::default();
        let (mut preprocessing_generation, mut preprocessing, preprocessing_clock) =
            norm_structure_fixture(&mut preprocessing_heap, 1);
        preprocessing.PrepInpData = std::array::from_fn(|_| Default::default());
        assert_eq!(
            NormOneStructureINChI(
                &mut preprocessing_heap,
                &mut CANON_GLOBALS::default(),
                preprocessing_clock,
                &mut preprocessing_generation,
                &mut preprocessing,
                INCHI_BAS as i32,
                None,
                0,
            ),
            Ok(0)
        );
        assert_eq!(
            preprocessing.PrepInpData[INCHI_BAS as usize].num_components,
            1
        );
        assert_eq!(
            preprocessing.ncFlags.bTautFlags[INCHI_BAS as usize],
            [preprocessing.StructData.bTautFlags[INCHI_BAS as usize]
                | preprocessing.InpParms.bTautFlags; TAUT_NUM as usize]
        );
        assert_eq!(
            preprocessing.ncFlags.bTautFlagsDone[INCHI_BAS as usize],
            [preprocessing.StructData.bTautFlagsDone[INCHI_BAS as usize]
                | preprocessing.InpParms.bTautFlagsDone; TAUT_NUM as usize]
        );

        let mut preprocessing_failure_heap = SourceHeap::default();
        let (
            mut preprocessing_failure_generation,
            mut preprocessing_failure,
            preprocessing_failure_clock,
        ) = norm_structure_fixture(&mut preprocessing_failure_heap, 1);
        preprocessing_failure.PrepInpData = std::array::from_fn(|_| Default::default());
        preprocessing_failure_heap.fail_after_allocations(0);
        assert_eq!(
            NormOneStructureINChI(
                &mut preprocessing_failure_heap,
                &mut CANON_GLOBALS::default(),
                preprocessing_failure_clock,
                &mut preprocessing_failure_generation,
                &mut preprocessing_failure,
                INCHI_BAS as i32,
                None,
                0,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(preprocessing_failure.StructData.nStructReadError, 99);
        assert_eq!(
            preprocessing_failure.StructData.nErrorType,
            _IS_ERROR as i32
        );

        let mut extraction_heap = SourceHeap::default();
        let (mut extraction_generation, mut extraction, extraction_clock) =
            norm_structure_fixture(&mut extraction_heap, 1);
        let lengths = extraction.PrepInpData[INCHI_BAS as usize].nCurAtLen;
        extraction_heap.slice_mut(lengths).unwrap()[0] = 2;
        assert_eq!(
            NormOneStructureINChI(
                &mut extraction_heap,
                &mut CANON_GLOBALS::default(),
                extraction_clock,
                &mut extraction_generation,
                &mut extraction,
                INCHI_BAS as i32,
                None,
                0,
            ),
            Ok(_IS_ERROR as i32)
        );
        assert_eq!(extraction.StructData.nErrorCode, CT_ATOMCOUNT_ERR);

        let mut quit_heap = SourceHeap::default();
        let (mut quit_generation, mut quit, quit_clock) = norm_structure_fixture(&mut quit_heap, 1);
        quit.StructData.bUserQuitComponent = 1;
        assert_eq!(
            NormOneStructureINChI(
                &mut quit_heap,
                &mut CANON_GLOBALS::default(),
                quit_clock,
                &mut quit_generation,
                &mut quit,
                INCHI_BAS as i32,
                None,
                0,
            ),
            Ok(0)
        );
        assert!(
            quit_heap
                .slice(quit.pINChI[INCHI_BAS as usize].as_const())
                .unwrap()[0][TAUT_NON as usize]
                .is_null()
        );
    }

    #[test]
    fn source_port__inchi_dll_a2__createcompatomdata__line_2232() {
        let mut simple_heap = SourceHeap::default();
        let old_atoms = simple_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let old_offsets = simple_heap
            .allocate_model_storage(vec![91_u16, 92_u16])
            .unwrap();
        let mut simple = COMP_ATOM_DATA {
            at: old_atoms,
            num_at: 77,
            num_removed_H: 76,
            nOffsetAtAndH: old_offsets,
            num_components: 75,
            ..COMP_ATOM_DATA::default()
        };
        assert_eq!(
            CreateCompAtomData(&mut simple_heap, &mut simple, 2, 1, 0),
            Ok(1)
        );
        assert_eq!(simple.num_at, 2);
        assert_eq!(simple.num_components, 0);
        assert_eq!(simple.num_removed_H, 0);
        assert!(!simple.at.is_null());
        assert!(simple.nOffsetAtAndH.is_null());
        assert_eq!(simple_heap.slice(simple.at.as_const()).unwrap().len(), 2);
        assert_eq!(simple_heap.live_allocation_count(), 1);
        FreeCompAtomData(&mut simple_heap, &mut simple).unwrap();

        for (components, intermediate, expected_components) in [(-7, 0, 0), (2, 1, 2), (2, -1, 2)] {
            let mut heap = SourceHeap::default();
            let mut data = COMP_ATOM_DATA::default();
            assert_eq!(
                CreateCompAtomData(&mut heap, &mut data, 1, components, intermediate),
                Ok(1)
            );
            assert_eq!(data.num_at, 1);
            assert_eq!(data.num_components, expected_components);
            assert!(!data.at.is_null());
            assert!(data.nOffsetAtAndH.is_null());
            assert_eq!(heap.live_allocation_count(), 1);
            FreeCompAtomData(&mut heap, &mut data).unwrap();
        }

        let mut composite_heap = SourceHeap::default();
        let mut composite = COMP_ATOM_DATA::default();
        assert_eq!(
            CreateCompAtomData(&mut composite_heap, &mut composite, 3, 4, 0),
            Ok(1)
        );
        assert_eq!(composite.num_at, 3);
        assert_eq!(composite.num_components, 4);
        assert_eq!(
            composite_heap
                .slice(composite.nOffsetAtAndH.as_const())
                .unwrap(),
            &[0_u16; 10]
        );
        assert_eq!(composite_heap.live_allocation_count(), 2);
        FreeCompAtomData(&mut composite_heap, &mut composite).unwrap();

        let mut zero_heap = SourceHeap::default();
        let mut zero = COMP_ATOM_DATA::default();
        assert_eq!(
            CreateCompAtomData(&mut zero_heap, &mut zero, 0, 0, 0),
            Ok(1)
        );
        assert_eq!(zero.num_at, 0);
        assert!(!zero.at.is_null());
        assert_eq!(zero_heap.slice(zero.at.as_const()).unwrap().len(), 0);
        FreeCompAtomData(&mut zero_heap, &mut zero).unwrap();

        let mut negative_heap = SourceHeap::default();
        let mut negative = COMP_ATOM_DATA {
            num_at: 19,
            ..COMP_ATOM_DATA::default()
        };
        assert_eq!(
            CreateCompAtomData(&mut negative_heap, &mut negative, -1, 2, 0),
            Ok(0)
        );
        assert_eq!(negative, COMP_ATOM_DATA::default());
        assert_eq!(negative_heap.live_allocation_count(), 0);

        let mut first_failure_heap = SourceHeap::default();
        first_failure_heap.fail_after_allocations(0);
        let mut first_failure = COMP_ATOM_DATA::default();
        assert_eq!(
            CreateCompAtomData(&mut first_failure_heap, &mut first_failure, 1, 2, 0),
            Ok(0)
        );
        assert_eq!(first_failure, COMP_ATOM_DATA::default());
        assert_eq!(first_failure_heap.source_allocation_calls(), 1);
        assert_eq!(first_failure_heap.live_allocation_count(), 0);

        let mut second_failure_heap = SourceHeap::default();
        second_failure_heap.fail_after_allocations(1);
        let mut second_failure = COMP_ATOM_DATA::default();
        assert_eq!(
            CreateCompAtomData(&mut second_failure_heap, &mut second_failure, 1, 2, 0),
            Ok(0)
        );
        assert_eq!(second_failure, COMP_ATOM_DATA::default());
        assert_eq!(second_failure_heap.source_allocation_calls(), 2);
        assert_eq!(second_failure_heap.live_allocation_count(), 0);
    }

    #[test]
    fn source_port__inchi_dll_a2__createcompositenormatom__line_1973() {
        let mut heap = SourceHeap::default();
        let mut component0_non = input_data(&mut heap, vec![atom(10, 0), atom(11, 0)], 1);
        component0_non.num_bonds = 1;
        component0_non.num_isotopic = 2;
        component0_non.bHasIsotopicLayer = 1;
        component0_non.nNumRemovedProtons = 3;
        component0_non.nNumRemovedProtonsIsotopic = [1, 2, 3];
        component0_non.num_iso_H = [4, 5, 6];

        let mut component0_taut = input_data(&mut heap, vec![atom(20, 0), atom(21, 0)], 1);
        component0_taut.at_fixed_bonds = heap
            .allocate_model_storage(vec![atom(30, 0), atom(31, 0)])
            .unwrap();
        component0_taut.bTautomeric = 1;
        component0_taut.bTautPreprocessed = 1;
        component0_taut.num_bonds = 2;
        component0_taut.num_isotopic = 3;
        component0_taut.bHasIsotopicLayer = 2;
        component0_taut.nNumRemovedProtons = 4;
        component0_taut.nNumRemovedProtonsIsotopic = [2, 3, 4];
        component0_taut.num_iso_H = [5, 6, 7];

        let mut component1_non = input_data(&mut heap, vec![atom(40, 0)], 0);
        component1_non.num_bonds = 3;
        let component1_taut = crate::source_types::INP_ATOM_DATA {
            bExists: 1,
            bDeleted: 1,
            nNumRemovedProtons: 7,
            nNumRemovedProtonsIsotopic: [1, 1, 1],
            ..crate::source_types::INP_ATOM_DATA::default()
        };

        let mut component2_non = input_data(&mut heap, vec![atom(50, 0)], 0);
        component2_non.num_bonds = 4;
        component2_non.num_isotopic = 5;
        component2_non.bHasIsotopicLayer = 4;
        component2_non.nNumRemovedProtons = 6;
        component2_non.nNumRemovedProtonsIsotopic = [3, 4, 5];
        component2_non.num_iso_H = [6, 7, 8];

        let inputs = [
            [component0_non, component0_taut],
            [component1_non, component1_taut],
            [
                component2_non,
                crate::source_types::INP_ATOM_DATA::default(),
            ],
        ];
        let mut composites = std::array::from_fn(|_| COMP_ATOM_DATA::default());
        assert_eq!(
            CreateCompositeNormAtom(&mut heap, &mut composites, &inputs, 3),
            Ok(7)
        );

        let non_atoms = heap.slice(composites[0].at.as_const()).unwrap();
        assert_eq!(
            non_atoms
                .iter()
                .map(|entry| entry.orig_at_number)
                .collect::<Vec<_>>(),
            vec![10, 40, 50, 11]
        );
        assert_eq!(non_atoms[1].neighbor[0], 1);
        assert_eq!(non_atoms[2].neighbor[0], 2);
        assert_eq!(non_atoms[3].neighbor[0], 0);
        assert_eq!(composites[0].num_at, 4);
        assert_eq!(composites[0].num_removed_H, 1);
        assert_eq!(composites[0].num_bonds, 8);
        assert_eq!(composites[0].bHasIsotopicLayer, 5);
        assert_eq!(
            heap.slice(composites[0].nOffsetAtAndH.as_const()).unwrap(),
            &[1, 4, 2, 4, 3, 4, 0, 0]
        );

        let taut_atoms = heap.slice(composites[1].at.as_const()).unwrap();
        assert_eq!(
            taut_atoms
                .iter()
                .map(|entry| entry.orig_at_number)
                .collect::<Vec<_>>(),
            vec![20, 50, 21]
        );
        assert_eq!(taut_atoms[1].neighbor[0], 1);
        assert_eq!(composites[1].bTautomeric, 1);
        assert_eq!(composites[1].nNumRemovedProtons, 17);
        assert_eq!(composites[1].nNumRemovedProtonsIsotopic, [6, 8, 10]);
        assert_eq!(composites[1].num_iso_H, [11, 13, 15]);
        assert_eq!(
            heap.slice(composites[1].nOffsetAtAndH.as_const()).unwrap(),
            &[1, 3, 0, 0, 2, 3, 0, 0]
        );

        let intermediate_atoms = heap.slice(composites[2].at.as_const()).unwrap();
        assert_eq!(
            intermediate_atoms
                .iter()
                .map(|entry| entry.orig_at_number)
                .collect::<Vec<_>>(),
            vec![30, 50, 31]
        );
        assert_eq!(intermediate_atoms[1].neighbor[0], 1);
        assert_eq!(composites[2].bTautomeric, 0);
        assert_eq!(composites[2].nNumRemovedProtons, 17);
        assert!(composites[2].nOffsetAtAndH.is_null());
        free_composites(&mut heap, &mut composites);

        let mut no_component_heap = SourceHeap::default();
        let mut no_component = std::array::from_fn(|index| COMP_ATOM_DATA {
            num_at: index as i32 + 70,
            ..COMP_ATOM_DATA::default()
        });
        let before = no_component.clone();
        assert_eq!(
            CreateCompositeNormAtom(&mut no_component_heap, &mut no_component, &[], -1),
            Ok(0)
        );
        assert_eq!(no_component, before);

        let mut first_failure_heap = SourceHeap::default();
        let first_failure_input = [[
            input_data(&mut first_failure_heap, vec![atom(1, 0)], 0),
            crate::source_types::INP_ATOM_DATA::default(),
        ]];
        first_failure_heap.fail_after_allocations(0);
        let mut first_failure = std::array::from_fn(|_| COMP_ATOM_DATA::default());
        assert_eq!(
            CreateCompositeNormAtom(
                &mut first_failure_heap,
                &mut first_failure,
                &first_failure_input,
                1,
            ),
            Ok(0)
        );
        assert_eq!(
            first_failure,
            std::array::from_fn(|_| COMP_ATOM_DATA::default())
        );

        let mut later_failure_heap = SourceHeap::default();
        let later_failure_input = [[input_data(&mut later_failure_heap, vec![atom(1, 0)], 0), {
            let mut data = input_data(&mut later_failure_heap, vec![atom(2, 0)], 0);
            data.bTautomeric = 0;
            data
        }]];
        later_failure_heap.fail_after_allocations(1);
        let mut later_failure = std::array::from_fn(|_| COMP_ATOM_DATA::default());
        assert_eq!(
            CreateCompositeNormAtom(
                &mut later_failure_heap,
                &mut later_failure,
                &later_failure_input,
                1,
            ),
            Ok(0)
        );
        assert!(!later_failure[0].at.is_null());
        assert!(later_failure[1].at.is_null());
        free_composites(&mut later_failure_heap, &mut later_failure);

        let mut miscount_heap = SourceHeap::default();
        let mut miscount_mobile = input_data(&mut miscount_heap, vec![atom(10, 0), atom(11, 0)], 0);
        miscount_mobile.bTautomeric = 1;
        let mut miscount_deleted_alt = input_data(
            &mut miscount_heap,
            vec![atom(20, 0), atom(21, 0), atom(22, 0), atom(23, 0)],
            0,
        );
        miscount_deleted_alt.bDeleted = 1;
        let trigger_non = crate::source_types::INP_ATOM_DATA::default();
        let mut trigger_mobile = input_data(&mut miscount_heap, vec![atom(30, 0)], 0);
        trigger_mobile.bTautomeric = 1;
        trigger_mobile.bTautPreprocessed = 1;
        trigger_mobile.at_fixed_bonds = miscount_heap
            .allocate_model_storage(vec![atom(31, 0)])
            .unwrap();
        let miscount_inputs = [
            [miscount_deleted_alt, miscount_mobile],
            [trigger_non, trigger_mobile],
        ];
        let mut miscount = std::array::from_fn(|_| COMP_ATOM_DATA::default());
        assert_eq!(
            CreateCompositeNormAtom(&mut miscount_heap, &mut miscount, &miscount_inputs, 2),
            Ok(3)
        );
        assert_eq!(miscount[0].bExists, 1);
        assert_eq!(
            miscount_heap.slice(miscount[0].at.as_const()).unwrap()[0].orig_at_number,
            30
        );
        assert_eq!(miscount[1].bExists, 1);
        assert_eq!(miscount[2].bExists, 0);
        free_composites(&mut miscount_heap, &mut miscount);
    }
}
