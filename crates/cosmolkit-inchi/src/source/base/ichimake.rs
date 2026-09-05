use crate::source::base::ichi_bns::mark_alt_bonds_and_taut_groups;
use crate::source::base::ichican2::{DeAllocBCN, GetBaseCanonRanking};
use crate::source::base::ichicano::{
    AllocateCS, Canon_INChI, DeAllocateCS, GetCanonLengths, UpdateFullLinearCT,
};
use crate::source::base::ichiisot::set_atom_iso_sort_keys;
use crate::source::base::ichimak2::FillOutINChI;
use crate::source::base::ichinorm::MarkRingSystemsInp;
use crate::source::base::ichiprt2::{Eql_INChI_Stereo, inchi_strtol};
use crate::source::base::ichisort::{
    CreateNeighListFromLinearCT, FreeNeighList, insertions_sort, insertions_sort_AT_RANK,
};
use crate::source::base::ichister::set_stereo_parity;
use crate::source::base::ichitaut::{
    CountTautomerGroups, make_a_copy_of_t_group_info, set_tautomer_iso_sort_keys,
};
use crate::source::base::mol2atom::FreeInpAtom;
use crate::source::base::strutil::{add_DT_to_num_H, remove_terminal_HDT};
use crate::source::base::util::{
    get_periodic_table_number, inchi_calloc, inchi_free, inchi_malloc,
};
use crate::source_types::{
    _IS_ERROR, _IS_FATAL, _IS_OKAY, _IS_WARNING, AB_PARITY_UNDF, AB_PARITY_UNKN, AT_NUMB, AT_RANK,
    ATOM_SIZES, ATT_PROTON, BCN, CANON_GLOBALS, CANON_MODE_CT, CANON_MODE_ISO,
    CANON_MODE_ISO_STEREO, CANON_MODE_STEREO, CANON_MODE_TAUT, CANON_STAT, CMODE_NO_ALT_SBONDS,
    CMODE_NOEQ_STEREO, CMODE_RACEMIC_STEREO, CMODE_REDNDNT_STEREO, CMODE_RELATIVE_STEREO,
    CMODE_SB_IGN_ALL_UU, CMODE_SC_IGN_ALL_UU, CT_CANON_ERR, CT_ERR_MAX, CT_ERR_MIN,
    CT_TAUCOUNT_ERR, EQL_SP2, EQL_SP3, ERR_NO_CANON_RESULTS, FLAG_NORM_CONSIDER_TAUT, ICR,
    ICR_MAX_DIFF_FIXED_H, ICR_MAX_ENDP_IN1_ONLY, ICR_MAX_ENDP_IN2_ONLY, ICR_MAX_SB_IN1_ONLY,
    ICR_MAX_SB_IN2_ONLY, ICR_MAX_SB_UNDF, ICR_MAX_SC_IN1_ONLY, ICR_MAX_SC_IN2_ONLY,
    ICR_MAX_SC_UNDF, INCHI_CLOCK, INCHI_FLAG_RAC_STEREO, INCHI_FLAG_REL_STEREO, INCHI_MODE,
    INCHI_SORT, INChI, INChI_Aux, INChI_Stereo, INP_ATOM_DATA, INPUT_PARMS, MAX_ATOMS,
    NUM_H_ISOTOPES, ORIG_ATOM_DATA, PES_BIT_ARSINE_STEREO, PES_BIT_FIX_SP3_BUG,
    PES_BIT_PHOSPHINE_STEREO, PES_BIT_POINT_EDGE_STEREO, REQ_MODE_BASIC, REQ_MODE_DEFAULT,
    REQ_MODE_DIFF_UU_STEREO, REQ_MODE_ISO, REQ_MODE_ISO_STEREO, REQ_MODE_NO_ALT_SBONDS,
    REQ_MODE_NOEQ_STEREO, REQ_MODE_NON_ISO, REQ_MODE_RACEMIC_STEREO, REQ_MODE_REDNDNT_STEREO,
    REQ_MODE_RELATIVE_STEREO, REQ_MODE_SB_IGN_ALL_UU, REQ_MODE_SC_IGN_ALL_UU, REQ_MODE_STEREO,
    REQ_MODE_TAUT, SourceConstPointer, SourceHeap, SourceHeapError, SourceMutPointer,
    StableSourceConstSlice, StableSourceSlice, T_GROUP_INFO, TAUT_NON, TAUT_NUM, TAUT_YES,
    TG_FLAG_ALL_TAUTOMERIC, TG_FLAG_ARSINE_STEREO, TG_FLAG_FIX_ISO_FIXEDH_BUG, TG_FLAG_FIX_SP3_BUG,
    TG_FLAG_FIX_TERM_H_CHRG_BUG, TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE, TG_FLAG_FOUND_ISOTOPIC_H_DONE,
    TG_FLAG_H_ALREADY_REMOVED, TG_FLAG_PHOSPHINE_STEREO, TG_FLAG_POINTED_EDGE_STEREO,
    WARN_FAILED_ISOTOPIC, WARN_FAILED_ISOTOPIC_STEREO, WARN_FAILED_STEREO, clock_t, inchiTime,
    inp_ATOM, sp_ATOM, tagDiffINChILayers_DIFL_F, tagDiffINChILayers_DIFL_FI,
    tagDiffINChILayers_DIFL_M, tagDiffINChILayers_DIFL_MI, tagDiffINChISegments_DIFS_b_SBONDS,
    tagDiffINChISegments_DIFS_f_FORMULA, tagDiffINChISegments_DIFS_h_H_ATOMS,
    tagDiffINChISegments_DIFS_i_IATOMS, tagDiffINChISegments_DIFS_m_SP3INV,
    tagDiffINChISegments_DIFS_o_TRANSP, tagDiffINChISegments_DIFS_q_CHARGE,
    tagDiffINChISegments_DIFS_s_STYPE, tagDiffINChISegments_DIFS_t_SATOMS,
    tagInchiDiffBits_IDIF_ATOMS, tagInchiDiffBits_IDIF_CHARGE, tagInchiDiffBits_IDIF_CON_LEN,
    tagInchiDiffBits_IDIF_CON_TBL, tagInchiDiffBits_IDIF_DIFF_TG_ENDP,
    tagInchiDiffBits_IDIF_EXTRA_TG_ENDP, tagInchiDiffBits_IDIF_ISO_AT,
    tagInchiDiffBits_IDIF_LESS_FH, tagInchiDiffBits_IDIF_LESS_H,
    tagInchiDiffBits_IDIF_MISS_TG_ENDP, tagInchiDiffBits_IDIF_MORE_FH,
    tagInchiDiffBits_IDIF_MORE_H, tagInchiDiffBits_IDIF_MULTIPLE_TG, tagInchiDiffBits_IDIF_NO_TAUT,
    tagInchiDiffBits_IDIF_NUM_AT, tagInchiDiffBits_IDIF_NUM_EL, tagInchiDiffBits_IDIF_NUM_ISO_AT,
    tagInchiDiffBits_IDIF_NUM_TG, tagInchiDiffBits_IDIF_POSITION_H, tagInchiDiffBits_IDIF_PROBLEM,
    tagInchiDiffBits_IDIF_REM_ISO_H, tagInchiDiffBits_IDIF_REM_PROT,
    tagInchiDiffBits_IDIF_SB_EXTRA, tagInchiDiffBits_IDIF_SB_EXTRA_UNDF,
    tagInchiDiffBits_IDIF_SB_MISS, tagInchiDiffBits_IDIF_SB_MISS_UNDF,
    tagInchiDiffBits_IDIF_SB_PARITY, tagInchiDiffBits_IDIF_SC_EXTRA,
    tagInchiDiffBits_IDIF_SC_EXTRA_UNDF, tagInchiDiffBits_IDIF_SC_INV,
    tagInchiDiffBits_IDIF_SC_MISS, tagInchiDiffBits_IDIF_SC_MISS_UNDF,
    tagInchiDiffBits_IDIF_SC_PARITY, tagInchiDiffBits_IDIF_SINGLE_TG, tagInchiDiffBits_IDIF_TG,
    tagInchiDiffBits_IDIF_WRONG_TAUT, tagMarkDiff_DIFV_BOTH_EMPTY, tagMarkDiff_DIFV_EQL2PRECED,
    tagMarkDiff_DIFV_FI_EQ_MI, tagMarkDiff_DIFV_IS_EMPTY, tagMarkDiff_DIFV_NEQ2PRECED,
    tagMarkDiff_DIFV_OUTPUT_OMIT_F,
};

fn copy_inp_atom_prefix(
    heap: &mut SourceHeap,
    destination: SourceMutPointer<inp_ATOM>,
    source: SourceConstPointer<inp_ATOM>,
    count: usize,
) -> Result<(), SourceHeapError> {
    if destination.as_const() == source {
        heap.slice(source)?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        return Ok(());
    }
    if destination.allocation_identity() != source.as_mut().allocation_identity() {
        return heap.with_slice_mut_and_heap(destination, |destination, heap| {
            let source = heap
                .slice(source)?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            destination
                .get_mut(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .copy_from_slice(source);
            Ok(())
        });
    }

    // Keep the established checked behavior for overlapping modeled pointers.
    // Normal source calls always pass distinct allocations to memcpy.
    let source = heap
        .slice(source)?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    heap.slice_mut(destination)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .copy_from_slice(&source);
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn inp2spATOM(
    heap: &mut SourceHeap,
    inp_at: SourceConstPointer<inp_ATOM>,
    num_inp_at: i32,
    at: SourceMutPointer<sp_ATOM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:119 inp2spATOM
    // INCHI✔️❌: complete source frame follows verbatim; input prefix cloning adds an allocation.
    /*
    int inp2spATOM(inp_ATOM* inp_at, int num_inp_at, sp_ATOM* at)
    {
        int i, j, val, elname_len;

        memset(at, 0, sizeof(at[0]) * num_inp_at); /* djb-rwth: memset_s C11/Annex K variant? */

        for (i = 0; i < num_inp_at; i++)
        {
            elname_len = sizeof(at[0].elname) - 1; /* djb-rwth: fixing coverity ID #499609 */
            strncpy(at[i].elname, inp_at[i].elname, elname_len);
            at[i].elname[elname_len] = '\0';
            at[i].el_number = (U_CHAR)get_periodic_table_number(at[i].elname);
            val = at[i].valence = inp_at[i].valence;
            for (j = 0; j < val; j++)
            {
                at[i].neighbor[j] = inp_at[i].neighbor[j];
                at[i].bond_type[j] = inp_at[i].bond_type[j];
            }
            at[i].chem_bonds_valence = inp_at[i].chem_bonds_valence;
            at[i].orig_at_number = inp_at[i].orig_at_number;
            at[i].orig_compt_at_numb = inp_at[i].orig_compt_at_numb;
            at[i].endpoint = inp_at[i].endpoint;
            at[i].iso_atw_diff = inp_at[i].iso_atw_diff;
            at[i].num_H = inp_at[i].num_H;
            at[i].cFlags = inp_at[i].cFlags;
            for (j = 0; j < NUM_H_ISOTOPES; j++)
            {
                at[i].num_iso_H[j] = inp_at[i].num_iso_H[j];
            }
            at[i].charge = inp_at[i].charge;
            at[i].radical = inp_at[i].radical;

    #if ( FIND_RING_SYSTEMS == 1 )
            at[i].nBlockSystem = inp_at[i].nBlockSystem;
            at[i].bCutVertex = inp_at[i].bCutVertex;
            at[i].nRingSystem = inp_at[i].nRingSystem;
            at[i].nNumAtInRingSystem = inp_at[i].nNumAtInRingSystem;
    #if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
            at[i].nDistanceFromTerminal = inp_at[i].nDistanceFromTerminal;
    #endif
    #endif

            /*
                    at[i].x                  = inp_at[i].x;
                    at[i].y                  = inp_at[i].y;
                    at[i].z                  = inp_at[i].z;
            */
        }

        return 0;
    }


    */
    // END INCHI C FUNCTION: inp2spATOM
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: inp2spATOM
    // INCHI✔️❌: #define FIND_RING_SYSTEMS 1
    // INCHI✔️❌: #define FIND_RINS_SYSTEMS_DISTANCES 0
    // END INCHI ACTIVE MACRO CONFIGURATION: inp2spATOM

    let count = usize::try_from(num_inp_at).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if count == 0 {
        return Ok(0);
    }
    let copy_atoms = |input: &[inp_ATOM], output: &mut [sp_ATOM]| {
        output.fill(sp_ATOM::default());
        for (source, destination) in input.iter().zip(output.iter_mut()) {
            let nul = source
                .elname
                .iter()
                .position(|byte| *byte == 0)
                .unwrap_or(source.elname.len());
            for index in 0..5 {
                destination.elname[index] = if index < nul { source.elname[index] } else { 0 };
            }
            destination.elname[5] = 0;
            destination.el_number = get_periodic_table_number(Some(&destination.elname))? as u8;
            destination.valence = source.valence;
            let valence = i32::from(source.valence);
            if valence > 20 {
                return Err(SourceHeapError::PointerOutOfBounds);
            }
            for index in 0..usize::try_from(valence.max(0))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
            {
                destination.neighbor[index] = source.neighbor[index];
                destination.bond_type[index] = source.bond_type[index];
            }
            destination.chem_bonds_valence = source.chem_bonds_valence;
            destination.orig_at_number = source.orig_at_number;
            destination.orig_compt_at_numb = source.orig_compt_at_numb;
            destination.endpoint = source.endpoint;
            destination.iso_atw_diff = source.iso_atw_diff;
            destination.num_H = source.num_H;
            destination.cFlags = source.cFlags;
            for index in 0..NUM_H_ISOTOPES as usize {
                destination.num_iso_H[index] = source.num_iso_H[index];
            }
            destination.charge = source.charge;
            destination.radical = source.radical;
            destination.nBlockSystem = source.nBlockSystem;
            destination.bCutVertex = source.bCutVertex;
            destination.nRingSystem = source.nRingSystem;
            destination.nNumAtInRingSystem = source.nNumAtInRingSystem;
        }
        Ok(())
    };

    if inp_at.as_mut().allocation_identity() != at.allocation_identity() {
        // SAFETY: Create_INChI allocates the input and output as distinct,
        // fixed buffers. inp2spATOM reads the first and writes the second.
        let input: StableSourceConstSlice<inp_ATOM> = unsafe { heap.stable_slice(inp_at)? };
        let mut output: StableSourceSlice<sp_ATOM> = unsafe { heap.stable_slice_mut(at)? };
        copy_atoms(input.prefix(count)?, output.prefix_mut(count)?)?;
        return Ok(0);
    }

    let input = heap
        .slice(inp_at)?
        .get(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    let output = heap
        .slice_mut(at)?
        .get_mut(..count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    copy_atoms(&input, output)?;
    Ok(0)
}

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn Create_INChI(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    ic: SourceMutPointer<INCHI_CLOCK>,
    ip: &INPUT_PARMS,
    ppINChI: [SourceMutPointer<INChI>; TAUT_NUM as usize],
    ppINChI_Aux: [SourceMutPointer<INChI_Aux>; TAUT_NUM as usize],
    _orig_inp_data: Option<&ORIG_ATOM_DATA>,
    inp_at: SourceMutPointer<inp_ATOM>,
    out_norm_data: &mut [INP_ATOM_DATA; TAUT_NUM as usize],
    num_inp_at: i32,
    mut nUserMode: INCHI_MODE,
    pbTautFlags: &mut INCHI_MODE,
    pbTautFlagsDone: &mut INCHI_MODE,
    ulMaxTime: SourceMutPointer<inchiTime>,
    mut ti_out: Option<&mut T_GROUP_INFO>,
    mut pStrErrStruct: Option<&mut [i8]>,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:3707 Create_INChI
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access and stack-object modeling add overhead.
    /*
    int  Create_INChI(CANON_GLOBALS* pCG,
        INCHI_CLOCK* ic,
        INPUT_PARMS* ip,
        INChI** ppINChI,
        INChI_Aux** ppINChI_Aux,
        ORIG_ATOM_DATA* orig_inp_data,
        inp_ATOM* inp_at,
        INP_ATOM_DATA* out_norm_data[2],
        int num_inp_at,
        INCHI_MODE nUserMode,
        INCHI_MODE* pbTautFlags,
        INCHI_MODE* pbTautFlagsDone,
        struct tagInchiTime* ulMaxTime,
        T_GROUP_INFO* ti_out,
        char* pStrErrStruct)
    {
        /*
        #define NON_TAUT 0
        #define TAUT     1
        */
        int nebend = 0, * ebend = NULL;
    
        sp_ATOM* at[TAUT_NUM]; /* at[0]=>non-tautomeric, at[1]=>tautomeric */
        int                       i, n1, n2, num_atoms, num_at_tg, num_removed_H, num_removed_H_taut = 0, ret = 0, ret2 = 0;
        INCHI_MODE                 nMode = 0;
        T_GROUP_INFO              vt_group_info;
        T_GROUP_INFO              vt_group_info_orig;
        T_GROUP_INFO* /*const*/  t_group_info = &vt_group_info;
        T_GROUP_INFO* /*const*/  t_group_info_orig = &vt_group_info_orig;
    
        CANON_STAT  CS, CS2;
        CANON_STAT* pCS = &CS;
        CANON_STAT* pCS2 = &CS2;  /*  save all allocations to avoid memory leaks in case Canon_INChI() removes the pointer */
    
        ATOM_SIZES  s[TAUT_NUM];
    
        BCN Bcn;
        BCN* pBCN = &Bcn;
    
        int bHasIsotopicAtoms = 0;
        int bMayHaveStereo = 0;
        int num_taut_at = 0;
    
        inp_ATOM* out_at = NULL;     /*, *norm_at_fixed_bonds[TAUT_NUM]; */ /*  = {out_norm_nontaut_at, out_norm_taut_at} ; */
        INChI* pINChI = NULL;      /* added initialization 2006-03 */
        INChI_Aux* pINChI_Aux = NULL;  /* added initialization 2006-03 */
        int        bPointedEdgeStereo = ((TG_FLAG_POINTED_EDGE_STEREO & *pbTautFlags) ? PES_BIT_POINT_EDGE_STEREO : 0)
            | ((TG_FLAG_PHOSPHINE_STEREO & *pbTautFlags) ? PES_BIT_PHOSPHINE_STEREO : 0)
            | ((TG_FLAG_ARSINE_STEREO & *pbTautFlags) ? PES_BIT_ARSINE_STEREO : 0)
            | ((TG_FLAG_FIX_SP3_BUG & *pbTautFlags) ? PES_BIT_FIX_SP3_BUG : 0);
        INCHI_MODE bTautFlags = (*pbTautFlags & (~(INCHI_MODE)TG_FLAG_ALL_TAUTOMERIC));
        INCHI_MODE bTautFlagsDone = (*pbTautFlagsDone /*& (~(INCHI_MODE)TG_FLAG_ALL_TAUTOMERIC) */);
    #if ( bRELEASE_VERSION == 0 )
        int bExtract = 0; /*  EXTR_HAS_ATOM_WITH_DEFINED_PARITY; */
    #endif
    
    #ifdef GHI100_FIX
    #if ((SPRINTF_FLAG != 1) && (SPRINTF_FLAG != 2))
        setlocale(LC_ALL, "en-US"); /* djb-rwth: setting all locales to "en-US" */
    #endif
    #endif
    
        /* */
        int bFixIsoFixedH = 0;
        int bFixTermHChrg = 0;
    
        int LargeMolecules = ip->bLargeMolecules;
        /* djb-rwth: removing redundant variables */
    
        /*    vABParityUnknown holds actual value of an internal constant signifying
            unknown parity: either the same as for undefined parity (default==standard)
            or a specific one (non-std; requested by SLUUD switch).                 */
        int vABParityUnknown = AB_PARITY_UNDF;
        if (0 != (nUserMode & REQ_MODE_DIFF_UU_STEREO))
        {
            /* Make labels for unknown and undefined stereo different */
            vABParityUnknown = AB_PARITY_UNKN;
        }
    
        /* djb-rwth: removing redundant code */
    
    #if ( FIX_ISO_FIXEDH_BUG == 1 )
        if (TG_FLAG_FIX_ISO_FIXEDH_BUG & *pbTautFlags)
            bFixIsoFixedH = 1;
    #endif
    #if ( FIX_TERM_H_CHRG_BUG == 1 )
        if (TG_FLAG_FIX_TERM_H_CHRG_BUG & *pbTautFlags)
            bFixTermHChrg = 1;
    #endif
    
    #ifdef FIX_SRU_CYCLIZING_PS_BONDS_IN_BNS
        /* Polymer related */
        if (orig_inp_data && orig_inp_data->polymer && orig_inp_data->polymer->n > 0 && orig_inp_data->polymer->valid)
        {
            int j, jj;
            nebend = 0;
            for (j = 0; j < orig_inp_data->polymer->n; j++)
                if (orig_inp_data->polymer->units[j]->cyclized == 1)
                    nebend++;
            if (nebend)
            {
                nebend *= 2;
                ebend = inchi_calloc(2 * nebend, sizeof(int));
                if (!ebend)
                {
                    ret = CT_OUT_OF_RAM; goto exit_function;
                }
                jj = 0;
                for (j = 0; j < orig_inp_data->polymer->n; j++)
                {
                    if (orig_inp_data->polymer->units[j]->cyclized == 1)
                    {
                        ebend[jj] = orig_inp_data->polymer->units[j]->end_atom1;
                        ebend[jj + 1] = orig_inp_data->polymer->units[j]->end_atom2;
                        jj += 2;
                    }
                }
            }
        }
    #endif
    
        memset(s, 0, sizeof(s)); /* djb-rwth: memset_s C11/Annex K variant? */
        if (pBCN)
        {
            memset(pBCN, 0, sizeof(pBCN[0])); /* djb-rwth: memset_s C11/Annex K variant? */
        }
        memset(t_group_info, 0, sizeof(*t_group_info)); /* djb-rwth: memset_s C11/Annex K variant? */
        memset(t_group_info_orig, 0, sizeof(*t_group_info_orig)); /* djb-rwth: memset_s C11/Annex K variant? */
        /*norm_at[TAUT_NON] = out_norm_data[TAUT_NON]->at; *//* output normalized non-tautomeric component */
        /*norm_at[TAUT_YES] = out_norm_data[TAUT_YES]->at; *//* output normalized tautomeric component */
        /*norm_at_fixed_bonds[TAUT_NON] = NULL;*/
        /*norm_at_fixed_bonds[TAUT_YES] = out_norm_data[TAUT_YES]->at_fixed_bonds;*/
        for (i = 0; i < TAUT_NUM; i++)
        {
            if (out_norm_data[i]->at)
            {
                if (!(at[i] = (sp_ATOM*)inchi_malloc(num_inp_at * sizeof(*at[0]))))
                {
                    ret = -1;
                }
            }
            else
            {
                at[i] = NULL;
            }
        }
    
        if ((!out_norm_data[TAUT_NON]->at && !out_norm_data[TAUT_YES]->at)
            || !inp_at || ret) /* djb-rwth: addressing LLVM warning */
        {
            ret = -1;
            goto exit_function;
        }
    
        /* the first struct to process: tautomeric if exists else non-tautomeric */
        out_at = out_norm_data[TAUT_YES]->at ? out_norm_data[TAUT_YES]->at : out_norm_data[TAUT_NON]->at;
        /* copy the input structure to be normalized to the buffer for the normalization data */
        memcpy(out_at, inp_at, num_inp_at * sizeof(out_at[0]));
        /*  tautomeric groups setting */
        t_group_info->bIgnoreIsotopic = 0;   /*  include tautomeric group isotopic info in MarkTautomerGroups() */
        t_group_info->bTautFlags = *pbTautFlags;
        t_group_info->bTautFlagsDone = *pbTautFlagsDone;
        t_group_info->t_group = NULL; /* djb-rwth: fixing oss-fuzz issue #70475 */
    
        /*
            Preprocess the structure
            (here THE NUMBER OF ATOMS MAY BE REDUCED)
        */
    
        /*  ??? Ambiguity: H-D may become HD or DH
            (that is, H+implicit D or D+implicit H)    */
        if (TG_FLAG_H_ALREADY_REMOVED & bTautFlags)
        {
            INP_ATOM_DATA* out_norm_data1 = out_norm_data[TAUT_YES]->at ? out_norm_data[TAUT_YES] :
                out_norm_data[TAUT_NON]->at ? out_norm_data[TAUT_NON] : NULL;
            if (out_norm_data1)
            {
                num_at_tg =
                    num_atoms = out_norm_data1->num_at - out_norm_data1->num_removed_H;
                num_removed_H = out_norm_data1->num_removed_H;
                t_group_info->tni.nNumRemovedExplicitH = num_removed_H;
            }
            else
            {
                ret = -1;
                goto exit_function;
            }
        }
        else
        {
            num_at_tg = num_atoms = remove_terminal_HDT(num_inp_at, out_at, bFixTermHChrg);
            num_removed_H = num_inp_at - num_atoms;
            t_group_info->tni.nNumRemovedExplicitH = num_removed_H;
            add_DT_to_num_H(num_atoms, out_at);
        }
        /*fix_odd_things( num_atoms, out_at );*/
    #if ( FIND_RING_SYSTEMS == 1 )
        MarkRingSystemsInp(out_at, num_atoms, 0);
    #endif
        /*  duplicate the preprocessed structure so that all supplied out_norm_data[]->at buffers are filled */
        if (out_at != out_norm_data[TAUT_YES]->at && out_norm_data[TAUT_YES]->at)
        {
            memcpy(out_norm_data[TAUT_YES]->at, out_at, num_inp_at * sizeof(out_at[0]));
        }
        if (out_norm_data[TAUT_YES]->at_fixed_bonds && out_norm_data[TAUT_YES]->at)
        {
            memcpy(out_norm_data[TAUT_YES]->at_fixed_bonds, out_at, num_inp_at * sizeof(out_at[0]));
        }
        if (out_at != out_norm_data[TAUT_NON]->at && out_norm_data[TAUT_NON]->at)
        {
            memcpy(out_norm_data[TAUT_NON]->at, out_at, num_inp_at * sizeof(out_at[0]));
        }
    
        /*
          ??? not true ??? duplicate inp_at and keep inp_at[] unchanged after terminal hydrogens removal
          set stereo parities in taut_at[], non_taut_at[]
          obtain max. lenghts of the name stereo parts
          Ignore absence/presence of isotopic stereo for now
          mark isotopic atoms
        */
        if (out_norm_data[TAUT_YES]->at && at[TAUT_YES])
        {
            /* final normalization of possibly tautomeric structure */
            ret = mark_alt_bonds_and_taut_groups(ic, pCG,
                out_norm_data[TAUT_YES]->at,
                out_norm_data[TAUT_YES]->at_fixed_bonds,
                num_atoms, ulMaxTime, t_group_info,
                NULL, NULL,
                nebend, ebend);
            if (ret < 0)
            {
                goto exit_function;/*  out of RAM or other normalization problem */
            }
            num_taut_at = ret; /* number of atoms without removed H? */
            num_removed_H_taut = t_group_info->tni.nNumRemovedExplicitH;
            out_norm_data[TAUT_YES]->num_at = num_atoms + num_removed_H_taut; /* protons might have been removed */
            out_norm_data[TAUT_YES]->num_removed_H = num_removed_H_taut;
            out_norm_data[TAUT_YES]->nNumRemovedProtons += t_group_info->tni.nNumRemovedProtons;
            for (i = 0; i < NUM_H_ISOTOPES; i++)
            {
                out_norm_data[TAUT_YES]->nNumRemovedProtonsIsotopic[i] +=
                    t_group_info->tni.nNumRemovedProtonsIsotopic[i] /*+ t_group_info->num_iso_H[i]*/;
                out_norm_data[TAUT_YES]->num_iso_H[i] +=
                    t_group_info->num_iso_H[i];
            }
            /* mark deleted isolated tautomeric H(+) */
            if (num_taut_at == 1 &&
                out_norm_data[TAUT_YES]->at[0].at_type == ATT_PROTON &&
                t_group_info && t_group_info->tni.nNumRemovedProtons == 1)
            {
                out_norm_data[TAUT_YES]->bDeleted = 1;
                FreeInpAtom(&out_norm_data[TAUT_YES]->at_fixed_bonds);
            }
            else
            {
                if ((t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT) &&
                    out_norm_data[TAUT_YES]->at_fixed_bonds)
                {
                    out_norm_data[TAUT_YES]->bTautPreprocessed = 1;
                }
            }
    
            /*
                if ( !(t_group_info->tni.bNormalizationFlags & (FLAG_NORM_CONSIDER_TAUT & ~FLAG_PROTON_SINGLE_REMOVED)) &&
                     out_norm_data[TAUT_YES]->at_fixed_bonds) {
                     FreeInpAtom( &out_norm_data[TAUT_YES]->at_fixed_bonds );
                }
            */
    
            /*out_norm_data[TAUT_YES]->num_removed_H = num_removed_H_taut;*/
    
            out_norm_data[TAUT_YES]->bTautFlags = *pbTautFlags = t_group_info->bTautFlags;
            out_norm_data[TAUT_YES]->bTautFlagsDone = *pbTautFlagsDone = t_group_info->bTautFlagsDone;
            out_norm_data[TAUT_YES]->bNormalizationFlags = t_group_info->tni.bNormalizationFlags;
    
            /* create internal sp_ATOM at[] out of out_norm_data[]->at */
            inp2spATOM(out_norm_data[TAUT_YES]->at, num_inp_at, at[TAUT_YES]);
    
            /* set stereo parities to at[]; nUserMode: accept alt. stereo bonds, min ring size */
            ret = set_stereo_parity(pCG, out_norm_data[TAUT_YES]->at, at[TAUT_YES],
                num_taut_at, num_removed_H_taut,
                &s[TAUT_YES].nMaxNumStereoAtoms,
                &s[TAUT_YES].nMaxNumStereoBonds,
                nUserMode, bPointedEdgeStereo,
                vABParityUnknown, ip->bLooseTSACheck, ip->bStereoAtZz);
    
    #if ( bRELEASE_VERSION == 0 )
            if (0 < ret)
            {
                bExtract |= EXTR_HAS_ATOM_WITH_DEFINED_PARITY;
            }
            if (t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT)
            {
                bExtract |= EXTR_TAUT_TREATMENT_CHARGES;
            }
    #endif
    
            if (RETURNED_ERROR(ret))
            {
                goto exit_function; /*  stereo bond error */
            }
    
            s[TAUT_YES].bMayHaveStereo =
                (s[TAUT_YES].nMaxNumStereoAtoms ||
                    s[TAUT_YES].nMaxNumStereoBonds);
    
            /*
                mark isotopic atoms and atoms that have non-tautomeric
                isotopic terminal hydrogen atoms 1H, 2H(D), 3H(T)
            */
    
            s[TAUT_YES].num_isotopic_atoms =
                set_atom_iso_sort_keys(num_taut_at,
                    at[TAUT_YES],
                    t_group_info,
                    &s[TAUT_YES].bHasIsotopicTautGroups);
    
            /*
                Prepare tautomeric (if no tautomerism found then prepare non-tautomeric)
                structure for canonicalizaton:
    
                    remove t-groups that have no H,
                    remove charges from t-groups if requested
                    renumber t-groups and find final t_group_info->num_t_groups
                    add to t-groups lists of endpoints tgroup->nEndpointAtomNumber[]
                    calculate length of the t-group part of the connection table
             */
    
            s[TAUT_YES].nLenLinearCTTautomer = CountTautomerGroups(at[TAUT_YES], num_taut_at, t_group_info);
    
            if (RETURNED_ERROR(s[TAUT_YES].nLenLinearCTTautomer))
            {
                /* added error treatment 9-11-2003 */
                ret = s[TAUT_YES].nLenLinearCTTautomer;
                goto exit_function;
                /*  error has happened; no breakpoint here
                s[TAUT_YES].nLenLinearCTTautomer = 0;
                */
            }
    
            else if (s[TAUT_YES].nLenLinearCTTautomer > 0)
            {
                num_at_tg = num_taut_at + t_group_info->num_t_groups;
                /*  ??? -not true- create t_group_info_orig for multiple calls with atom renumbering */
                make_a_copy_of_t_group_info(t_group_info_orig /* dest*/, t_group_info /* source*/); /* djb-rwth: addressing coverity ID #499544 -- properly used sequence of arguments according to the comment in previous line */
                /*  mark isotopic tautomer groups: calculate t_group->iWeight */
                s[TAUT_YES].nLenLinearCTIsotopicTautomer = set_tautomer_iso_sort_keys(t_group_info);
                if (s[TAUT_YES].nLenLinearCTIsotopicTautomer < 0)
                {
                    /* ??? -error cannot happen- error has happened; no breakpoint here */
                    s[TAUT_YES].nLenLinearCTIsotopicTautomer = 0;
                }
                out_norm_data[TAUT_YES]->bTautomeric = s[TAUT_YES].nLenLinearCTTautomer;
            }
    
            /*  new variable: s[TAUT_YES].nLenCT introduced 7-22-2002 */
            GetCanonLengths(num_taut_at, at[TAUT_YES], &s[TAUT_YES], t_group_info);
        }
    
        if (out_norm_data[TAUT_NON]->at &&
            out_norm_data[TAUT_YES]->at &&
            at[TAUT_NON] &&
            !s[TAUT_YES].nLenLinearCTTautomer)
        {
            /* the structure is non-tautomeric: use tautomeric treatment results only for it */
            inchi_free(at[TAUT_NON]);
            at[TAUT_NON] = NULL;
        }
    
        else if (!out_norm_data[TAUT_NON]->at &&
            out_norm_data[TAUT_YES]->at &&
            !at[TAUT_NON] &&
            at[TAUT_YES] &&
            !s[TAUT_YES].nLenLinearCTTautomer)
        {
            /* requested tautomeric; found non-tautomeric; it is located in out_norm_data[TAUT_YES]->at */
            out_norm_data[TAUT_YES]->bTautomeric = 0;
        }
    
        else if (out_norm_data[TAUT_NON]->at && at[TAUT_NON])
        {
            /* the structure needs non-tautomeric treatment:
            final normalization of non-tautomeric structure */
            ret = mark_alt_bonds_and_taut_groups(ic, pCG,
                out_norm_data[TAUT_NON]->at,
                NULL,
                num_atoms,
                ulMaxTime,
                NULL,
                &bTautFlags,
                &bTautFlagsDone,
                nebend, ebend);
            if (ret < 0)
            {
                goto exit_function;  /*  out of RAM or other normalization problem */
            }
            out_norm_data[TAUT_NON]->num_at = num_atoms + num_removed_H;
            out_norm_data[TAUT_NON]->num_removed_H = num_removed_H;
            out_norm_data[TAUT_NON]->bTautFlags = *pbTautFlags;
            out_norm_data[TAUT_NON]->bTautFlagsDone = *pbTautFlagsDone;
            out_norm_data[TAUT_NON]->bNormalizationFlags = 0;
    
            /* create internal sp_ATOM at[] out of out_norm_data[]->at */
            inp2spATOM(out_norm_data[TAUT_NON]->at, num_inp_at, at[TAUT_NON]);
    
            /* set stereo parities to at[]; nUserMode: accept alt. stereo bonds, min ring size */
            ret = set_stereo_parity(pCG, out_norm_data[TAUT_NON]->at,
                at[TAUT_NON], num_atoms, num_removed_H,
                &s[TAUT_NON].nMaxNumStereoAtoms,
                &s[TAUT_NON].nMaxNumStereoBonds, nUserMode,
                bPointedEdgeStereo, vABParityUnknown,
                ip->bLooseTSACheck, ip->bStereoAtZz);
    #if ( bRELEASE_VERSION == 0 )
            if (0 < ret)
            {
                bExtract |= EXTR_HAS_ATOM_WITH_DEFINED_PARITY;
            }
    #endif
            if (RETURNED_ERROR(ret))
            {
                goto exit_function; /*  stereo bond error */
            }
            s[TAUT_NON].bMayHaveStereo = (s[TAUT_NON].nMaxNumStereoAtoms ||
                s[TAUT_NON].nMaxNumStereoBonds);
    
            /*
             * mark isotopic atoms and atoms that have non-tautomeric
             * isotopic terminal hydrogen atoms 1H, 2H(D), 3H(T)
             */
            s[TAUT_NON].num_isotopic_atoms = set_atom_iso_sort_keys(num_atoms, at[TAUT_NON], NULL, NULL);
            GetCanonLengths(num_atoms, at[TAUT_NON], &s[TAUT_NON], NULL);
            out_norm_data[TAUT_NON]->bTautomeric = 0;
        }
    
        /* common  */
        bMayHaveStereo = s[TAUT_YES].bMayHaveStereo || s[TAUT_NON].bMayHaveStereo;
        bHasIsotopicAtoms = s[TAUT_NON].num_isotopic_atoms > 0 || s[TAUT_NON].bHasIsotopicTautGroups > 0 ||
            s[TAUT_YES].num_isotopic_atoms > 0 || s[TAUT_YES].bHasIsotopicTautGroups > 0;
        if (bFixIsoFixedH)
        {
            /* 2008-03-21 DT */
            bHasIsotopicAtoms = bHasIsotopicAtoms
                ||
                (s[TAUT_YES].nLenLinearCTTautomer > 0 && t_group_info &&
                    ((0 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[0]) ||
                        (1 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[1]) ||
                        (2 < NUM_H_ISOTOPES && t_group_info->tni.nNumRemovedProtonsIsotopic[2])
                        )); /* djb-rwth: addressing LLVM warning */
        }
        bHasIsotopicAtoms = bHasIsotopicAtoms
            ||
            (s[TAUT_YES].nLenIsotopicEndpoints > 1 && t_group_info &&
                (t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))); /* djb-rwth: addressing LLVM warning */
    
        /* default mode */
        if (!(nUserMode & REQ_MODE_DEFAULT))
        {
            /*  default */
            nUserMode |= REQ_MODE_DEFAULT;
        }
    
        /* adjust the mode to the reality */
        if ((nUserMode & REQ_MODE_ISO) && !bHasIsotopicAtoms)
        {
            nUserMode ^= REQ_MODE_ISO;
            nUserMode |= REQ_MODE_NON_ISO;  /*  at least one is needed */
        }
        if ((nUserMode & REQ_MODE_STEREO) && (nUserMode & REQ_MODE_ISO))
        {
            nUserMode |= REQ_MODE_ISO_STEREO;
        }
        if ((nUserMode & REQ_MODE_STEREO) && !(nUserMode & REQ_MODE_NON_ISO))
        {
            nUserMode ^= REQ_MODE_STEREO;
        }
        if (!bMayHaveStereo)
        {
            if (nUserMode & REQ_MODE_STEREO)
                nUserMode ^= REQ_MODE_STEREO;
            if (nUserMode & REQ_MODE_ISO_STEREO)
                nUserMode ^= REQ_MODE_ISO_STEREO;
        }
    
        if ((nUserMode & REQ_MODE_BASIC) &&
            (!out_norm_data[TAUT_NON]->at || !ppINChI[TAUT_NON] ||
                !ppINChI_Aux[TAUT_NON] || !at[TAUT_NON]))
        {
            nUserMode ^= REQ_MODE_BASIC;
        }
        if ((nUserMode & REQ_MODE_TAUT) &&
            (!out_norm_data[TAUT_YES]->at || !ppINChI[TAUT_YES] ||
                !ppINChI_Aux[TAUT_YES] || !at[TAUT_YES]))
        {
            nUserMode ^= REQ_MODE_TAUT;
        }
    
        switch ((int)nUserMode & (REQ_MODE_BASIC | REQ_MODE_TAUT))
        {
        case REQ_MODE_BASIC:
            n1 = TAUT_NON;
            n2 = TAUT_NON;
            break;
        case REQ_MODE_TAUT:
            n1 = TAUT_YES;
            n2 = TAUT_YES;
            break;
        case (REQ_MODE_BASIC | REQ_MODE_TAUT):
            n1 = TAUT_NON;
            n2 = TAUT_YES;
            break;
        default:
            ret = -3;
            goto exit_function; /*  program error: inconsistent nUserMode or missing taut/non-taut allocation */ /*   <BRKPT> */
        }
    
        /*
            Obtain all non-stereo canonical numberings
        */
    
        if ((nUserMode & REQ_MODE_NON_ISO) && !(nUserMode & REQ_MODE_ISO))
        {
            /* added for special non-isotopic test mode 2004-10-04 */
            if (t_group_info)
            {
                t_group_info->bIgnoreIsotopic = 1;
                if (t_group_info->nIsotopicEndpointAtomNumber)
                {
                    t_group_info->nIsotopicEndpointAtomNumber[0] = inchi_min(1, t_group_info->nIsotopicEndpointAtomNumber[0]);
                }
                memset(t_group_info->num_iso_H, 0, sizeof(t_group_info->num_iso_H)); /* djb-rwth: memset_s C11/Annex K variant? */
                memset(t_group_info->tni.nNumRemovedProtonsIsotopic, 0, sizeof(t_group_info->tni.nNumRemovedProtonsIsotopic)); /* djb-rwth: memset_s C11/Annex K variant? */
                t_group_info->bTautFlagsDone &= ~(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE);
            }
            for (i = 0; i < TAUT_NUM; i++)
            {
                s[i].bHasIsotopicTautGroups = 0;
                s[i].bIgnoreIsotopic = 1;
                s[i].nLenIsotopic = 0;
                s[i].nLenIsotopicEndpoints = 0;
                s[i].nLenLinearCTIsotopicTautomer = 0;
                s[i].num_isotopic_atoms = 0;
            }
            bHasIsotopicAtoms = 0;
        }
    
        ret = GetBaseCanonRanking(ic, num_atoms, num_at_tg, at,
            t_group_info, s, pBCN, ulMaxTime,
            pCG, bFixIsoFixedH, LargeMolecules);
    
        if (ret < 0)
        {
            goto exit_function; /*  program error */
        }
    #if ( bRELEASE_VERSION == 0 && FIND_CANON_NE_EQUITABLE == 1 )
        /* Debug only: find whether canonical equivalence is different from equitable partition */
        if (bCanonIsFinerThanEquitablePartition(num_atoms, at[n1], pBCN->ftcn[TAUT_NON].nSymmRankCt))
        {
            bExtract |= EXTR_CANON_NE_EQUITABLE;
        }
    #endif
    
        /* added for special non-isotopic test mode 2004-10-04 */
        if (!pBCN->ftcn[n1].PartitionCt.Rank)
        {
            n1 = ALT_TAUT(n1);
        }
        if (!pBCN->ftcn[n2].PartitionCt.Rank)
        {
            n2 = ALT_TAUT(n2);
        }
        if (n1 > n2)
        {
            ret = CT_TAUCOUNT_ERR;
            goto exit_function; /*  program error */
        }
    
        /*
            Obtain stereo canonical numberings
        */
    
        for (i = n2; i >= n1 && !RETURNED_ERROR(ret); i--)
        {
            memset(pCS, 0, sizeof(*pCS)); /* djb-rwth: memset_s C11/Annex K variant? */
    
            switch (i)
            {
            case TAUT_NON:
                /*  non-tautomeric */
                /* djb-rwth: removing redundant code */
                nMode = (s[i].nLenLinearCTTautomer == 0) ? CANON_MODE_CT : CANON_MODE_TAUT;
                nMode |= (bHasIsotopicAtoms && (nUserMode & REQ_MODE_ISO)) ? CANON_MODE_ISO : 0;
                nMode |= (s[TAUT_NON].bMayHaveStereo && (nUserMode & REQ_MODE_STEREO)) ? CANON_MODE_STEREO : 0;
                nMode |= (bHasIsotopicAtoms && s[TAUT_NON].bMayHaveStereo && (nUserMode & REQ_MODE_ISO_STEREO)) ? CANON_MODE_ISO_STEREO : 0;
                nMode |= (nUserMode & REQ_MODE_NOEQ_STEREO) ? CMODE_NOEQ_STEREO : 0;
                nMode |= (nUserMode & REQ_MODE_REDNDNT_STEREO) ? CMODE_REDNDNT_STEREO : 0;
                nMode |= (nUserMode & REQ_MODE_NO_ALT_SBONDS) ? CMODE_NO_ALT_SBONDS : 0;
                if ((nMode & CANON_MODE_STEREO) == CANON_MODE_STEREO ||
                    (nMode & CANON_MODE_ISO_STEREO) == CANON_MODE_ISO_STEREO)
                {
                    nMode |= (nUserMode & REQ_MODE_RELATIVE_STEREO) ? CMODE_RELATIVE_STEREO : 0;
                    nMode |= (nUserMode & REQ_MODE_RACEMIC_STEREO) ? CMODE_RACEMIC_STEREO : 0;
                    nMode |= (nUserMode & REQ_MODE_SC_IGN_ALL_UU) ? CMODE_SC_IGN_ALL_UU : 0;
                    nMode |= (nUserMode & REQ_MODE_SB_IGN_ALL_UU) ? CMODE_SB_IGN_ALL_UU : 0;
                }
                if ((ret = AllocateCS(pCS, num_atoms, num_atoms, s[TAUT_NON].nLenCT, s[TAUT_NON].nLenCTAtOnly,
                    s[TAUT_NON].nLenLinearCTStereoDble, s[TAUT_NON].nMaxNumStereoBonds,
                    s[TAUT_NON].nLenLinearCTStereoCarb, s[TAUT_NON].nMaxNumStereoAtoms,
                    0, 0, s[TAUT_NON].nLenIsotopic, nMode, pBCN))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                *pCS2 = *pCS;
                break;
            case TAUT_YES: /*  tautomeric */
                /* djb-rwth: removing redundant code */
                nMode = (s[i].nLenLinearCTTautomer == 0) ? CANON_MODE_CT : CANON_MODE_TAUT;
                nMode |= (bHasIsotopicAtoms && (nUserMode & REQ_MODE_ISO)) ? CANON_MODE_ISO : 0;
                nMode |= (s[TAUT_YES].bMayHaveStereo && (nUserMode & REQ_MODE_STEREO)) ? CANON_MODE_STEREO : 0;
                nMode |= (bHasIsotopicAtoms && s[TAUT_YES].bMayHaveStereo && (nUserMode & REQ_MODE_ISO_STEREO)) ? CANON_MODE_ISO_STEREO : 0;
                nMode |= (nUserMode & REQ_MODE_NOEQ_STEREO) ? CMODE_NOEQ_STEREO : 0;
                nMode |= (nUserMode & REQ_MODE_REDNDNT_STEREO) ? CMODE_REDNDNT_STEREO : 0;
                nMode |= (nUserMode & REQ_MODE_NO_ALT_SBONDS) ? CMODE_NO_ALT_SBONDS : 0;
                if ((nMode & CANON_MODE_STEREO) == CANON_MODE_STEREO ||
                    (nMode & CANON_MODE_ISO_STEREO) == CANON_MODE_ISO_STEREO)
                {
                    nMode |= (nUserMode & REQ_MODE_RELATIVE_STEREO) ? CMODE_RELATIVE_STEREO : 0;
                    nMode |= (nUserMode & REQ_MODE_RACEMIC_STEREO) ? CMODE_RACEMIC_STEREO : 0;
                    nMode |= (nUserMode & REQ_MODE_SC_IGN_ALL_UU) ? CMODE_SC_IGN_ALL_UU : 0;
                    nMode |= (nUserMode & REQ_MODE_SB_IGN_ALL_UU) ? CMODE_SB_IGN_ALL_UU : 0;
                }
                if ((ret = AllocateCS(pCS, num_atoms, num_at_tg, s[TAUT_YES].nLenCT, s[TAUT_YES].nLenCTAtOnly,
                    s[TAUT_YES].nLenLinearCTStereoDble, s[TAUT_YES].nMaxNumStereoBonds,
                    s[TAUT_YES].nLenLinearCTStereoCarb, s[TAUT_YES].nMaxNumStereoAtoms,
                    s[TAUT_YES].nLenLinearCTTautomer, s[TAUT_YES].nLenLinearCTIsotopicTautomer,
                    s[TAUT_YES].nLenIsotopic, nMode, pBCN))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                *pCS2 = *pCS;
                break;
            }
    
            /* 2009-12-05 */
            nMode |= (nUserMode & REQ_MODE_DIFF_UU_STEREO) ? REQ_MODE_DIFF_UU_STEREO : 0;
            /* 2009-12-05 */
    
            /*  settings */
            pCS->lNumDecreasedCT = -1;
            pCS->bDoubleBondSquare = DOUBLE_BOND_NEIGH_LIST ? 2 : 0;  /*  2 => special mode */
            pCS->bIgnoreIsotopic = !((s[TAUT_NON].num_isotopic_atoms ||
                s[TAUT_YES].num_isotopic_atoms ||
                s[TAUT_YES].bHasIsotopicTautGroups) ||
                (nUserMode & REQ_MODE_NON_ISO) ||
                !(nUserMode & REQ_MODE_ISO));
    
            if ((nUserMode & REQ_MODE_NON_ISO) && !(nUserMode & REQ_MODE_ISO))
            {
                pCS->bIgnoreIsotopic = 1; /* 10-04-2004 */
            }
    
            if (i == TAUT_YES)
            {
                /* tautomeric */
                pCS->t_group_info = t_group_info; /*  ??? make a copy or reuse ???  */
                pCS->t_group_info->bIgnoreIsotopic = !(s[TAUT_YES].bHasIsotopicTautGroups ||
                    (nUserMode & REQ_MODE_NON_ISO) ||
                    !(nUserMode & REQ_MODE_ISO));
                if ((nUserMode & REQ_MODE_NON_ISO) && !(nUserMode & REQ_MODE_ISO))
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
    
             /*
                 The last canonicalization step
             */
    
            if (pBCN)
            {
                /* USE_CANON2 == 1 */
                pCS->NeighList = NULL;
                pCS->pBCN = pBCN;
    
                ret = Canon_INChI(ic,
                    num_atoms,
                    i ? num_at_tg : num_atoms,
                    at[i], pCS,
                    pCG,
                    nMode, i);
            }
            else
            {
                /* old way */
                pCS->NeighList = CreateNeighList(num_atoms,
                    i ? num_at_tg : num_atoms,
                    at[i],
                    pCS->bDoubleBondSquare,
                    pCS->t_group_info);
                pCS->pBCN = NULL;
    
                ret = Canon_INChI(ic,
                    num_atoms,
                    i ? num_at_tg : num_atoms,
                    at[i], pCS,
                    pCG,
                    nMode, i);
            }
    
            pINChI = ppINChI[i];      /* pointers to already allocated still empty InChI */
            pINChI_Aux = ppINChI_Aux[i];
    
            if (ret <= 0)
            {
                /*
                    Failure in Canon_INChI()
                */
                pINChI->nErrorCode = ret;
                pINChI_Aux->nErrorCode = ret;
            }
            else
            {
                /*
                    Success Canon_INChI()
    
                    save canonicalization results in
                    pINChI and pINChI_Aux
                */
                pINChI->nErrorCode = 0;
                pINChI_Aux->nErrorCode = 0;
                pINChI->bDeleted = pINChI_Aux->bDeleted = out_norm_data[i]->bDeleted;
                pINChI_Aux->nCanonFlags = pCS->nCanonFlags;
                pINChI_Aux->bTautFlags = out_norm_data[i]->bTautFlags;
                pINChI_Aux->bTautFlagsDone = out_norm_data[i]->bTautFlagsDone;
                pINChI_Aux->bNormalizationFlags = out_norm_data[i]->bNormalizationFlags;
    
                /*  may return an error or a warning */
                ret = FillOutINChI(pINChI, pINChI_Aux,
                    num_atoms, i ? num_at_tg : num_atoms,
                    i ? num_removed_H_taut : num_removed_H, at[i],
                    out_norm_data[i]->at, pCS,
                    pCG,
                    i, nUserMode,
                    pStrErrStruct, ip->bNoWarnings);
    
                if (RETURNED_ERROR(ret))
                {
                    /* Failure in FillOutINChI() */
                    pINChI->nErrorCode = ret;
                    pINChI_Aux->nErrorCode = ret;
                }
                else
                {
                    /* Success in FillOutINChI() */
    
    #if ( bRELEASE_VERSION == 0 )
                    if (pINChI->Stereo &&
                        (pINChI->Stereo->nCompInv2Abs && !pINChI->Stereo->bTrivialInv) ||
                        pINChI->StereoIsotopic &&
                        (pINChI->StereoIsotopic->nCompInv2Abs && !pINChI->StereoIsotopic->bTrivialInv))
                    {
                        bExtract |= EXTR_NON_TRIVIAL_STEREO;
                    }
    #endif
                    /*    Mark non-tautomeric representation as having
                        another, tautomeric representation */
                    if (pINChI_Aux && s[TAUT_YES].nLenLinearCTTautomer)
                    {
                        pINChI_Aux->bIsTautomeric = s[TAUT_YES].nLenLinearCTTautomer;
                    }
    #if ( bRELEASE_VERSION == 0 )
                    pCS->bExtract |= bExtract;
                    pINChI->bExtract |= pCS->bExtract;
    #endif
    
                    ret2 = CheckCanonNumberingCorrectness(num_atoms,
                        i ? num_at_tg : num_atoms,
                        at[i], pCS,
                        pCG,
                        i, pStrErrStruct);
                    if (ret2 && pINChI_Aux) /* djb-rwth: fixing a NULL pointer dereference */
                    {
                        pINChI->nErrorCode = ret2;
                        pINChI_Aux->nErrorCode = ret2;
                        ret = ret2;
                    }
                }
            }
    
            FreeNeighList(pCS->NeighList);
            DeAllocateCS(pCS2);
    
            pINChI = NULL;      /* avoid dangling pointers */
            pINChI_Aux = NULL;  /* avoid dangling pointers */
        }
    
        if (ret == 0)
        {
            ret = num_atoms;
        }
        /*  treat the results later */
    
    exit_function:
    
        DeAllocBCN(pBCN);
        if (at[TAUT_YES])
        {
            inchi_free(at[TAUT_YES]);
        }
        if (at[TAUT_NON])
        {
            inchi_free(at[TAUT_NON]);
        }
        if (ti_out)
        {
            *ti_out = *t_group_info;
        }
        else
        {
            /* free_t_group_info(t_group_info); */
            if (t_group_info) /* djb-rwth: fixing oss-fuzz issue #42537161/70475 */
            {
                if (t_group_info->nEndpointAtomNumber)
                {
                    inchi_free(t_group_info->nEndpointAtomNumber);
                }
                if (t_group_info->tGroupNumber)
                {
                    inchi_free(t_group_info->tGroupNumber);
                }
                if (t_group_info->nIsotopicEndpointAtomNumber)
                {
                    inchi_free(t_group_info->nIsotopicEndpointAtomNumber);
                }
                memset(t_group_info, 0, sizeof(*t_group_info)); /* djb-rwth: memset_s C11/Annex K variant? */
            }
        }
        free_t_group_info(t_group_info_orig);
    
        if (ebend)
        {
            inchi_free(ebend);
        }
    
        return ret;
    }
    */
    // END INCHI C FUNCTION: Create_INChI
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: Create_INChI
    // INCHI✔️❌: GCC/Linux: bRELEASE_VERSION == 1, SPRINTF_FLAG == 2.
    // INCHI✔️❌: #define FIX_ISO_FIXEDH_BUG 1
    // INCHI✔️❌: #define FIX_TERM_H_CHRG_BUG 1
    // INCHI✔️❌: #define FIND_RING_SYSTEMS 1
    // INCHI✔️❌: #define FIND_CANON_NE_EQUITABLE 0
    // INCHI✔️❌: #define DOUBLE_BOND_NEIGH_LIST 0
    // INCHI✔️❌: GHI100_FIX and FIX_SRU_CYCLIZING_PS_BONDS_IN_BNS are undefined.
    // END INCHI ACTIVE MACRO CONFIGURATION: Create_INChI
    let returned_error = |value: i32| (CT_ERR_MIN..=CT_ERR_MAX).contains(&value);
    let mut at = [SourceMutPointer::<sp_ATOM>::null(); TAUT_NUM as usize];
    let mut sizes = std::array::from_fn(|_| ATOM_SIZES::default());
    let mut t_group_info = T_GROUP_INFO::default();
    let mut t_group_info_orig = T_GROUP_INFO::default();
    let mut bcn = BCN::default();
    let mut cs = CANON_STAT::default();
    let mut cs_allocations = CANON_STAT::default();
    let mut bcn_storage = SourceMutPointer::<BCN>::null();
    let mut tgi_storage = SourceMutPointer::<T_GROUP_INFO>::null();
    let mut ret = 0_i32;

    let computation = (|| -> Result<i32, SourceHeapError> {
        if num_inp_at < 0 {
            return Ok(-1);
        }
        let count = usize::try_from(num_inp_at)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let mut clock = heap
            .slice(ic.as_const())?
            .first()
            .cloned()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let b_pointed_edge_stereo =
            i32::from((*pbTautFlags & u64::from(TG_FLAG_POINTED_EDGE_STEREO)) != 0)
                * PES_BIT_POINT_EDGE_STEREO as i32
                | i32::from((*pbTautFlags & u64::from(TG_FLAG_PHOSPHINE_STEREO)) != 0)
                    * PES_BIT_PHOSPHINE_STEREO as i32
                | i32::from((*pbTautFlags & u64::from(TG_FLAG_ARSINE_STEREO)) != 0)
                    * PES_BIT_ARSINE_STEREO as i32
                | i32::from((*pbTautFlags & u64::from(TG_FLAG_FIX_SP3_BUG)) != 0)
                    * PES_BIT_FIX_SP3_BUG as i32;
        let mut b_taut_flags = *pbTautFlags & !u64::from(TG_FLAG_ALL_TAUTOMERIC);
        let mut b_taut_flags_done = *pbTautFlagsDone;
        let b_fix_iso_fixed_h =
            i32::from((*pbTautFlags & u64::from(TG_FLAG_FIX_ISO_FIXEDH_BUG)) != 0);
        let b_fix_term_h_charge =
            i32::from((*pbTautFlags & u64::from(TG_FLAG_FIX_TERM_H_CHRG_BUG)) != 0);
        let large_molecules = ip.bLargeMolecules;
        let v_ab_parity_unknown = if nUserMode & u64::from(REQ_MODE_DIFF_UU_STEREO) != 0 {
            AB_PARITY_UNKN as i32
        } else {
            AB_PARITY_UNDF as i32
        };

        for index in 0..TAUT_NUM as usize {
            if !out_norm_data[index].at.is_null() {
                at[index] = match inchi_calloc::<sp_ATOM>(
                    heap,
                    u64::try_from(num_inp_at)
                        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?,
                    std::mem::size_of::<sp_ATOM>() as u64,
                ) {
                    Ok(pointer) => pointer,
                    Err(SourceHeapError::AllocationFailed) => {
                        ret = -1;
                        SourceMutPointer::null()
                    }
                    Err(error) => return Err(error),
                };
            }
        }
        if (out_norm_data[TAUT_NON as usize].at.is_null()
            && out_norm_data[TAUT_YES as usize].at.is_null())
            || inp_at.is_null()
            || ret != 0
        {
            return Ok(-1);
        }

        let out_at = if !out_norm_data[TAUT_YES as usize].at.is_null() {
            out_norm_data[TAUT_YES as usize].at
        } else {
            out_norm_data[TAUT_NON as usize].at
        };
        // INCHI✔️✔️:     memcpy(out_at, inp_at, num_inp_at * sizeof(out_at[0]));
        copy_inp_atom_prefix(heap, out_at, inp_at.as_const(), count)?;

        t_group_info.bIgnoreIsotopic = 0;
        t_group_info.bTautFlags = *pbTautFlags;
        t_group_info.bTautFlagsDone = *pbTautFlagsDone;
        t_group_info.t_group = SourceMutPointer::null();

        let mut num_atoms;
        let mut num_at_tg;
        let num_removed_h;
        if b_taut_flags & u64::from(TG_FLAG_H_ALREADY_REMOVED) != 0 {
            let selected = if !out_norm_data[TAUT_YES as usize].at.is_null() {
                &out_norm_data[TAUT_YES as usize]
            } else {
                &out_norm_data[TAUT_NON as usize]
            };
            num_atoms = selected.num_at.wrapping_sub(selected.num_removed_H);
            num_at_tg = num_atoms;
            num_removed_h = selected.num_removed_H;
            t_group_info.tni.nNumRemovedExplicitH = num_removed_h as i16;
        } else {
            num_atoms = remove_terminal_HDT(heap, num_inp_at, out_at, b_fix_term_h_charge)?;
            num_at_tg = num_atoms;
            num_removed_h = num_inp_at.wrapping_sub(num_atoms);
            t_group_info.tni.nNumRemovedExplicitH = num_removed_h as i16;
            let atoms = heap.slice_mut(out_at)?;
            add_DT_to_num_H(num_atoms, atoms)?;
        }
        MarkRingSystemsInp(heap, out_at, num_atoms, 0)?;

        for index in 0..TAUT_NUM as usize {
            let destination = out_norm_data[index].at;
            if !destination.is_null() && destination != out_at {
                // INCHI✔️✔️:         memcpy(out_norm_data[TAUT_YES]->at, out_at, num_inp_at * sizeof(out_at[0]));
                // INCHI✔️✔️:         memcpy(out_norm_data[TAUT_NON]->at, out_at, num_inp_at * sizeof(out_at[0]));
                copy_inp_atom_prefix(heap, destination, out_at.as_const(), count)?;
            }
        }
        if !out_norm_data[TAUT_YES as usize].at_fixed_bonds.is_null()
            && !out_norm_data[TAUT_YES as usize].at.is_null()
        {
            // INCHI✔️✔️:         memcpy(out_norm_data[TAUT_YES]->at_fixed_bonds, out_at, num_inp_at * sizeof(out_at[0]));
            copy_inp_atom_prefix(
                heap,
                out_norm_data[TAUT_YES as usize].at_fixed_bonds,
                out_at.as_const(),
                count,
            )?;
        }

        let mut num_removed_h_taut = 0_i32;
        let mut num_taut_at = 0_i32;
        if !out_norm_data[TAUT_YES as usize].at.is_null() && !at[TAUT_YES as usize].is_null() {
            tgi_storage = heap.allocate_model_storage(vec![t_group_info.clone()])?;
            let normalization_result = mark_alt_bonds_and_taut_groups(
                heap,
                ic,
                pCG,
                out_norm_data[TAUT_YES as usize].at,
                out_norm_data[TAUT_YES as usize].at_fixed_bonds,
                num_atoms,
                ulMaxTime,
                tgi_storage,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
                SourceMutPointer::null(),
                clock_result,
            );
            ret = normalization_result?;
            t_group_info = heap.slice(tgi_storage.as_const())?[0].clone();
            if ret < 0 {
                return Ok(ret);
            }
            num_taut_at = ret;
            num_removed_h_taut = i32::from(t_group_info.tni.nNumRemovedExplicitH);
            let taut = &mut out_norm_data[TAUT_YES as usize];
            taut.num_at = num_atoms.wrapping_add(num_removed_h_taut);
            taut.num_removed_H = num_removed_h_taut;
            taut.nNumRemovedProtons = taut
                .nNumRemovedProtons
                .wrapping_add(i32::from(t_group_info.tni.nNumRemovedProtons));
            for isotope in 0..NUM_H_ISOTOPES as usize {
                taut.nNumRemovedProtonsIsotopic[isotope] = taut.nNumRemovedProtonsIsotopic
                    [isotope]
                    .wrapping_add(t_group_info.tni.nNumRemovedProtonsIsotopic[isotope]);
                taut.num_iso_H[isotope] =
                    taut.num_iso_H[isotope].wrapping_add(t_group_info.num_iso_H[isotope]);
            }
            if num_taut_at == 1
                && heap.slice(taut.at.as_const())?[0].at_type == ATT_PROTON as AT_NUMB
                && t_group_info.tni.nNumRemovedProtons == 1
            {
                taut.bDeleted = 1;
                FreeInpAtom(heap, Some(&mut taut.at_fixed_bonds))?;
            } else if t_group_info.tni.bNormalizationFlags
                & u64::from(FLAG_NORM_CONSIDER_TAUT)
                != 0
                && !taut.at_fixed_bonds.is_null()
            {
                taut.bTautPreprocessed = 1;
            }
            taut.bTautFlags = t_group_info.bTautFlags;
            *pbTautFlags = t_group_info.bTautFlags;
            taut.bTautFlagsDone = t_group_info.bTautFlagsDone;
            *pbTautFlagsDone = t_group_info.bTautFlagsDone;
            taut.bNormalizationFlags = t_group_info.tni.bNormalizationFlags;

            inp2spATOM(
                    heap,
                    taut.at.as_const(),
                    num_inp_at,
                    at[TAUT_YES as usize],
                )?;
            ret = set_stereo_parity(
                pCG,
                heap,
                taut.at,
                at[TAUT_YES as usize],
                num_taut_at,
                num_removed_h_taut,
                Some(&mut sizes[TAUT_YES as usize].nMaxNumStereoAtoms),
                Some(&mut sizes[TAUT_YES as usize].nMaxNumStereoBonds),
                nUserMode,
                b_pointed_edge_stereo,
                v_ab_parity_unknown,
                ip.bLooseTSACheck,
                ip.bStereoAtZz,
            )?;
            if returned_error(ret) {
                return Ok(ret);
            }
            sizes[TAUT_YES as usize].bMayHaveStereo = i32::from(
                sizes[TAUT_YES as usize].nMaxNumStereoAtoms != 0
                    || sizes[TAUT_YES as usize].nMaxNumStereoBonds != 0,
            );
            sizes[TAUT_YES as usize].num_isotopic_atoms = set_atom_iso_sort_keys(
                heap,
                num_taut_at,
                at[TAUT_YES as usize],
                Some(&t_group_info),
                Some(&mut sizes[TAUT_YES as usize].bHasIsotopicTautGroups),
            )?;
            sizes[TAUT_YES as usize].nLenLinearCTTautomer = CountTautomerGroups(
                heap,
                at[TAUT_YES as usize],
                num_taut_at,
                Some(&mut t_group_info),
            )?;
            if returned_error(sizes[TAUT_YES as usize].nLenLinearCTTautomer) {
                return Ok(sizes[TAUT_YES as usize].nLenLinearCTTautomer);
            } else if sizes[TAUT_YES as usize].nLenLinearCTTautomer > 0 {
                num_at_tg = num_taut_at.wrapping_add(t_group_info.num_t_groups);
                let _ = make_a_copy_of_t_group_info(
                    heap,
                    Some(&mut t_group_info_orig),
                    Some(&mut t_group_info),
                )?;
                sizes[TAUT_YES as usize].nLenLinearCTIsotopicTautomer =
                    set_tautomer_iso_sort_keys(heap, Some(&mut t_group_info))?;
                if sizes[TAUT_YES as usize].nLenLinearCTIsotopicTautomer < 0 {
                    sizes[TAUT_YES as usize].nLenLinearCTIsotopicTautomer = 0;
                }
                out_norm_data[TAUT_YES as usize].bTautomeric =
                    sizes[TAUT_YES as usize].nLenLinearCTTautomer;
            }
            GetCanonLengths(
                heap,
                num_taut_at,
                at[TAUT_YES as usize].as_const(),
                &mut sizes[TAUT_YES as usize],
                Some(&t_group_info),
            )?;
        }

        if !out_norm_data[TAUT_NON as usize].at.is_null()
            && !out_norm_data[TAUT_YES as usize].at.is_null()
            && !at[TAUT_NON as usize].is_null()
            && sizes[TAUT_YES as usize].nLenLinearCTTautomer == 0
        {
            inchi_free(heap, at[TAUT_NON as usize])?;
            at[TAUT_NON as usize] = SourceMutPointer::null();
        } else if out_norm_data[TAUT_NON as usize].at.is_null()
            && !out_norm_data[TAUT_YES as usize].at.is_null()
            && at[TAUT_NON as usize].is_null()
            && !at[TAUT_YES as usize].is_null()
            && sizes[TAUT_YES as usize].nLenLinearCTTautomer == 0
        {
            out_norm_data[TAUT_YES as usize].bTautomeric = 0;
        } else if !out_norm_data[TAUT_NON as usize].at.is_null()
            && !at[TAUT_NON as usize].is_null()
        {
            let external_flags = heap.allocate_model_storage(vec![b_taut_flags])?;
            let external_flags_done = heap.allocate_model_storage(vec![b_taut_flags_done])?;
            ret = mark_alt_bonds_and_taut_groups(
                heap,
                ic,
                pCG,
                out_norm_data[TAUT_NON as usize].at,
                SourceMutPointer::null(),
                num_atoms,
                ulMaxTime,
                SourceMutPointer::null(),
                external_flags,
                external_flags_done,
                0,
                SourceMutPointer::null(),
                clock_result,
            )?;
            b_taut_flags = heap.slice(external_flags.as_const())?[0];
            b_taut_flags_done = heap.slice(external_flags_done.as_const())?[0];
            inchi_free(heap, external_flags)?;
            inchi_free(heap, external_flags_done)?;
            if ret < 0 {
                return Ok(ret);
            }
            let non = &mut out_norm_data[TAUT_NON as usize];
            non.num_at = num_atoms.wrapping_add(num_removed_h);
            non.num_removed_H = num_removed_h;
            non.bTautFlags = *pbTautFlags;
            non.bTautFlagsDone = *pbTautFlagsDone;
            non.bNormalizationFlags = 0;
            inp2spATOM(heap, non.at.as_const(), num_inp_at, at[TAUT_NON as usize])?;
            ret = set_stereo_parity(
                pCG,
                heap,
                non.at,
                at[TAUT_NON as usize],
                num_atoms,
                num_removed_h,
                Some(&mut sizes[TAUT_NON as usize].nMaxNumStereoAtoms),
                Some(&mut sizes[TAUT_NON as usize].nMaxNumStereoBonds),
                nUserMode,
                b_pointed_edge_stereo,
                v_ab_parity_unknown,
                ip.bLooseTSACheck,
                ip.bStereoAtZz,
            )?;
            if returned_error(ret) {
                return Ok(ret);
            }
            sizes[TAUT_NON as usize].bMayHaveStereo = i32::from(
                sizes[TAUT_NON as usize].nMaxNumStereoAtoms != 0
                    || sizes[TAUT_NON as usize].nMaxNumStereoBonds != 0,
            );
            sizes[TAUT_NON as usize].num_isotopic_atoms = set_atom_iso_sort_keys(
                heap,
                num_atoms,
                at[TAUT_NON as usize],
                None,
                None,
            )?;
            GetCanonLengths(
                heap,
                num_atoms,
                at[TAUT_NON as usize].as_const(),
                &mut sizes[TAUT_NON as usize],
                None,
            )?;
            non.bTautomeric = 0;
        }

        let b_may_have_stereo = sizes[TAUT_YES as usize].bMayHaveStereo != 0
            || sizes[TAUT_NON as usize].bMayHaveStereo != 0;
        let mut b_has_isotopic_atoms = sizes[TAUT_NON as usize].num_isotopic_atoms > 0
            || sizes[TAUT_NON as usize].bHasIsotopicTautGroups > 0
            || sizes[TAUT_YES as usize].num_isotopic_atoms > 0
            || sizes[TAUT_YES as usize].bHasIsotopicTautGroups > 0;
        if b_fix_iso_fixed_h != 0 {
            b_has_isotopic_atoms |= sizes[TAUT_YES as usize].nLenLinearCTTautomer > 0
                && t_group_info
                    .tni
                    .nNumRemovedProtonsIsotopic
                    .iter()
                    .any(|value| *value != 0);
        }
        b_has_isotopic_atoms |= sizes[TAUT_YES as usize].nLenIsotopicEndpoints > 1
            && t_group_info.bTautFlagsDone
                & u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE)
                != 0;

        if nUserMode & u64::from(REQ_MODE_DEFAULT) == 0 {
            nUserMode |= u64::from(REQ_MODE_DEFAULT);
        }
        if nUserMode & u64::from(REQ_MODE_ISO) != 0 && !b_has_isotopic_atoms {
            nUserMode ^= u64::from(REQ_MODE_ISO);
            nUserMode |= u64::from(REQ_MODE_NON_ISO);
        }
        if nUserMode & u64::from(REQ_MODE_STEREO) != 0
            && nUserMode & u64::from(REQ_MODE_ISO) != 0
        {
            nUserMode |= u64::from(REQ_MODE_ISO_STEREO);
        }
        if nUserMode & u64::from(REQ_MODE_STEREO) != 0
            && nUserMode & u64::from(REQ_MODE_NON_ISO) == 0
        {
            nUserMode ^= u64::from(REQ_MODE_STEREO);
        }
        if !b_may_have_stereo {
            nUserMode &= !u64::from(REQ_MODE_STEREO | REQ_MODE_ISO_STEREO);
        }
        if nUserMode & u64::from(REQ_MODE_BASIC) != 0
            && (out_norm_data[TAUT_NON as usize].at.is_null()
                || ppINChI[TAUT_NON as usize].is_null()
                || ppINChI_Aux[TAUT_NON as usize].is_null()
                || at[TAUT_NON as usize].is_null())
        {
            nUserMode ^= u64::from(REQ_MODE_BASIC);
        }
        if nUserMode & u64::from(REQ_MODE_TAUT) != 0
            && (out_norm_data[TAUT_YES as usize].at.is_null()
                || ppINChI[TAUT_YES as usize].is_null()
                || ppINChI_Aux[TAUT_YES as usize].is_null()
                || at[TAUT_YES as usize].is_null())
        {
            nUserMode ^= u64::from(REQ_MODE_TAUT);
        }
        let requested = nUserMode & u64::from(REQ_MODE_BASIC | REQ_MODE_TAUT);
        let (mut n1, mut n2) = match requested {
            value if value == u64::from(REQ_MODE_BASIC) => (TAUT_NON, TAUT_NON),
            value if value == u64::from(REQ_MODE_TAUT) => (TAUT_YES, TAUT_YES),
            value if value == u64::from(REQ_MODE_BASIC | REQ_MODE_TAUT) => {
                (TAUT_NON, TAUT_YES)
            }
            _ => return Ok(-3),
        };

        if nUserMode & u64::from(REQ_MODE_NON_ISO) != 0
            && nUserMode & u64::from(REQ_MODE_ISO) == 0
        {
            t_group_info.bIgnoreIsotopic = 1;
            if !t_group_info.nIsotopicEndpointAtomNumber.is_null() {
                let first = &mut heap.slice_mut(t_group_info.nIsotopicEndpointAtomNumber)?[0];
                *first = (*first).min(1);
            }
            t_group_info.num_iso_H.fill(0);
            t_group_info.tni.nNumRemovedProtonsIsotopic.fill(0);
            t_group_info.bTautFlagsDone &=
                !u64::from(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE);
            for size in &mut sizes {
                size.bHasIsotopicTautGroups = 0;
                size.bIgnoreIsotopic = 1;
                size.nLenIsotopic = 0;
                size.nLenIsotopicEndpoints = 0;
                size.nLenLinearCTIsotopicTautomer = 0;
                size.num_isotopic_atoms = 0;
            }
            b_has_isotopic_atoms = false;
        }

        ret = GetBaseCanonRanking(
            heap,
            &mut clock,
            num_atoms,
            num_at_tg,
            at,
            Some(&mut t_group_info),
            &sizes,
            &mut bcn,
            ulMaxTime,
            pCG,
            b_fix_iso_fixed_h,
            large_molecules,
            None,
            None,
            clock_result,
        )?;
        if ret < 0 {
            return Ok(ret);
        }
        if bcn.ftcn[n1 as usize].PartitionCt.Rank.is_null() {
            n1 = 1 - n1;
        }
        if bcn.ftcn[n2 as usize].PartitionCt.Rank.is_null() {
            n2 = 1 - n2;
        }
        if n1 > n2 {
            return Ok(CT_TAUCOUNT_ERR);
        }
        bcn_storage = heap.allocate_model_storage(vec![bcn.clone()])?;

        for index in (n1..=n2).rev() {
            if returned_error(ret) {
                break;
            }
            cs = CANON_STAT::default();
            let i = index as usize;
            let mut n_mode = if sizes[i].nLenLinearCTTautomer == 0 {
                CANON_MODE_CT
            } else {
                CANON_MODE_TAUT
            };
            if b_has_isotopic_atoms && nUserMode & u64::from(REQ_MODE_ISO) != 0 {
                n_mode |= CANON_MODE_ISO;
            }
            if sizes[i].bMayHaveStereo != 0 && nUserMode & u64::from(REQ_MODE_STEREO) != 0 {
                n_mode |= CANON_MODE_STEREO;
            }
            if b_has_isotopic_atoms
                && sizes[i].bMayHaveStereo != 0
                && nUserMode & u64::from(REQ_MODE_ISO_STEREO) != 0
            {
                n_mode |= CANON_MODE_ISO_STEREO;
            }
            n_mode |= if nUserMode & u64::from(REQ_MODE_NOEQ_STEREO) != 0 {
                CMODE_NOEQ_STEREO
            } else {
                0
            };
            n_mode |= if nUserMode & u64::from(REQ_MODE_REDNDNT_STEREO) != 0 {
                CMODE_REDNDNT_STEREO
            } else {
                0
            };
            n_mode |= if nUserMode & u64::from(REQ_MODE_NO_ALT_SBONDS) != 0 {
                CMODE_NO_ALT_SBONDS
            } else {
                0
            };
            if n_mode & CANON_MODE_STEREO == CANON_MODE_STEREO
                || n_mode & CANON_MODE_ISO_STEREO == CANON_MODE_ISO_STEREO
            {
                n_mode |= if nUserMode & u64::from(REQ_MODE_RELATIVE_STEREO) != 0 {
                    CMODE_RELATIVE_STEREO
                } else {
                    0
                };
                n_mode |= if nUserMode & u64::from(REQ_MODE_RACEMIC_STEREO) != 0 {
                    CMODE_RACEMIC_STEREO
                } else {
                    0
                };
                n_mode |= if nUserMode & u64::from(REQ_MODE_SC_IGN_ALL_UU) != 0 {
                    CMODE_SC_IGN_ALL_UU
                } else {
                    0
                };
                n_mode |= if nUserMode & u64::from(REQ_MODE_SB_IGN_ALL_UU) != 0 {
                    CMODE_SB_IGN_ALL_UU
                } else {
                    0
                };
            }
            let allocation_result = AllocateCS(
                heap,
                &mut cs,
                num_atoms,
                if index == TAUT_YES { num_at_tg } else { num_atoms },
                sizes[i].nLenCT,
                sizes[i].nLenCTAtOnly,
                sizes[i].nLenLinearCTStereoDble,
                sizes[i].nMaxNumStereoBonds,
                sizes[i].nLenLinearCTStereoCarb,
                sizes[i].nMaxNumStereoAtoms,
                if index == TAUT_YES {
                    sizes[i].nLenLinearCTTautomer
                } else {
                    0
                },
                if index == TAUT_YES {
                    sizes[i].nLenLinearCTIsotopicTautomer
                } else {
                    0
                },
                sizes[i].nLenIsotopic,
                u64::from(n_mode),
                bcn_storage,
            )?;
            cs_allocations = cs.clone();
            if allocation_result != 0 {
                ret = allocation_result;
                break;
            }
            n_mode |= if nUserMode & u64::from(REQ_MODE_DIFF_UU_STEREO) != 0 {
                REQ_MODE_DIFF_UU_STEREO
            } else {
                0
            };
            cs.lNumDecreasedCT = -1;
            cs.bDoubleBondSquare = 0;
            cs.bIgnoreIsotopic = i32::from(
                !((sizes[TAUT_NON as usize].num_isotopic_atoms != 0
                    || sizes[TAUT_YES as usize].num_isotopic_atoms != 0
                    || sizes[TAUT_YES as usize].bHasIsotopicTautGroups != 0)
                    || nUserMode & u64::from(REQ_MODE_NON_ISO) != 0
                    || nUserMode & u64::from(REQ_MODE_ISO) == 0),
            );
            if nUserMode & u64::from(REQ_MODE_NON_ISO) != 0
                && nUserMode & u64::from(REQ_MODE_ISO) == 0
            {
                cs.bIgnoreIsotopic = 1;
            }
            if index == TAUT_YES {
                heap.slice_mut(tgi_storage)?[0] = t_group_info.clone();
                cs.t_group_info = tgi_storage;
                t_group_info.bIgnoreIsotopic = i32::from(
                    !(sizes[TAUT_YES as usize].bHasIsotopicTautGroups != 0
                        || nUserMode & u64::from(REQ_MODE_NON_ISO) != 0
                        || nUserMode & u64::from(REQ_MODE_ISO) == 0),
                );
                if nUserMode & u64::from(REQ_MODE_NON_ISO) != 0
                    && nUserMode & u64::from(REQ_MODE_ISO) == 0
                {
                    t_group_info.bIgnoreIsotopic = 1;
                }
                heap.slice_mut(tgi_storage)?[0] = t_group_info.clone();
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
                num_atoms,
                if index == TAUT_YES { num_at_tg } else { num_atoms },
                at[i],
                &mut cs,
                pCG,
                n_mode,
                index as i32,
            )?;
            if index == TAUT_YES {
                t_group_info = heap.slice(tgi_storage.as_const())?[0].clone();
            }

            let mut inchi = heap.slice(ppINChI[i].as_const())?[0].clone();
            let mut aux = heap.slice(ppINChI_Aux[i].as_const())?[0].clone();
            if ret <= 0 {
                inchi.nErrorCode = ret;
                aux.nErrorCode = ret;
            } else {
                inchi.nErrorCode = 0;
                aux.nErrorCode = 0;
                inchi.bDeleted = out_norm_data[i].bDeleted;
                aux.bDeleted = out_norm_data[i].bDeleted;
                aux.nCanonFlags = cs.nCanonFlags;
                aux.bTautFlags = out_norm_data[i].bTautFlags;
                aux.bTautFlagsDone = out_norm_data[i].bTautFlagsDone;
                aux.bNormalizationFlags = out_norm_data[i].bNormalizationFlags;
                ret = FillOutINChI(
                    heap,
                    &mut inchi,
                    &mut aux,
                    num_atoms,
                    if index == TAUT_YES { num_at_tg } else { num_atoms },
                    if index == TAUT_YES {
                        num_removed_h_taut
                    } else {
                        num_removed_h
                    },
                    at[i],
                    out_norm_data[i].at,
                    &mut cs,
                    pCG,
                    index as i32,
                    nUserMode,
                    pStrErrStruct.as_deref_mut(),
                    ip.bNoWarnings,
                )?;
                if returned_error(ret) {
                    inchi.nErrorCode = ret;
                    aux.nErrorCode = ret;
                } else {
                    if sizes[TAUT_YES as usize].nLenLinearCTTautomer != 0 {
                        aux.bIsTautomeric = sizes[TAUT_YES as usize].nLenLinearCTTautomer;
                    }
                    let ret2 = CheckCanonNumberingCorrectness(
                        heap,
                        num_atoms,
                        if index == TAUT_YES { num_at_tg } else { num_atoms },
                        at[i].as_const(),
                        &mut cs,
                        pCG,
                        index as i32,
                        pStrErrStruct.as_deref_mut(),
                    )?;
                    if ret2 != 0 {
                        inchi.nErrorCode = ret2;
                        aux.nErrorCode = ret2;
                        ret = ret2;
                    }
                }
            }
            heap.slice_mut(ppINChI[i])?[0] = inchi;
            heap.slice_mut(ppINChI_Aux[i])?[0] = aux;
            DeAllocateCS(heap, &mut cs_allocations)?;
            cs = CANON_STAT::default();
            cs_allocations = CANON_STAT::default();
        }
        if ret == 0 {
            ret = num_atoms;
        }
        Ok(ret)
    })();

    DeAllocateCS(heap, &mut cs_allocations)?;
    if !bcn_storage.is_null() {
        bcn = heap.slice(bcn_storage.as_const())?[0].clone();
    }
    DeAllocBCN(heap, Some(&mut bcn))?;
    if !bcn_storage.is_null() {
        inchi_free(heap, bcn_storage)?;
    }
    for pointer in at {
        if !pointer.is_null() {
            inchi_free(heap, pointer)?;
        }
    }
    if let Some(output) = ti_out.as_deref_mut() {
        *output = t_group_info;
    } else {
        if !t_group_info.nEndpointAtomNumber.is_null() {
            inchi_free(heap, t_group_info.nEndpointAtomNumber)?;
        }
        if !t_group_info.tGroupNumber.is_null() {
            inchi_free(heap, t_group_info.tGroupNumber)?;
        }
        if !t_group_info.nIsotopicEndpointAtomNumber.is_null() {
            inchi_free(heap, t_group_info.nIsotopicEndpointAtomNumber)?;
        }
    }
    if !t_group_info_orig.t_group.is_null() {
        inchi_free(heap, t_group_info_orig.t_group)?;
    }
    if !t_group_info_orig.nEndpointAtomNumber.is_null() {
        inchi_free(heap, t_group_info_orig.nEndpointAtomNumber)?;
    }
    if !t_group_info_orig.tGroupNumber.is_null() {
        inchi_free(heap, t_group_info_orig.tGroupNumber)?;
    }
    if !t_group_info_orig.nIsotopicEndpointAtomNumber.is_null() {
        inchi_free(heap, t_group_info_orig.nIsotopicEndpointAtomNumber)?;
    }
    if !tgi_storage.is_null() {
        inchi_free(heap, tgi_storage)?;
    }
    computation
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn CheckCanonNumberingCorrectness(
    heap: &mut SourceHeap,
    num_atoms: i32,
    num_at_tg: i32,
    at: SourceConstPointer<sp_ATOM>,
    pCS: &mut CANON_STAT,
    pCG: &mut CANON_GLOBALS,
    _bTautomeric: i32,
    _pStrErrStruct: Option<&mut [i8]>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:6230 CheckCanonNumberingCorrectness
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int CheckCanonNumberingCorrectness(int num_atoms,
        int num_at_tg,
        sp_ATOM* at,
        CANON_STAT* pCS,
        CANON_GLOBALS* pCG,
        int bTautomeric,
        char* pStrErrStruct)
    {
        int i, ret = 0;
        AT_NUMB* pCanonOrd = NULL;
        int nErrorCode = 0;
        AT_NUMB* pCanonRank; /* canonical ranks of the atoms or tautomeric groups */
        AT_NUMB* pCanonRankAtoms = NULL;

        /* djb-rwth: removing redundant code */

        pCanonRankAtoms = (AT_NUMB*)inchi_calloc((long long)num_at_tg + 1, sizeof(pCanonRankAtoms[0])); /* djb-rwth: cast operator added */

        /*
            Non-isotopic part
        */

        pCanonOrd = pCS->nLenCanonOrdStereo > 0 ? pCS->nCanonOrdStereo :
            pCS->nLenCanonOrd > 0 ? pCS->nCanonOrd : NULL;
        pCanonRank = pCanonRankAtoms;
        if (pCanonOrd && pCanonRank)
        {
            for (i = 0; i < num_at_tg; i++)
            {
                pCanonRank[pCanonOrd[i]] = (AT_NUMB)(i + 1);
            }
            ret = UpdateFullLinearCT(num_atoms, num_at_tg, at, pCanonRank, pCanonOrd, pCS,
                pCG,
                0);
            if (ret /*|| memcmp(pCS->LinearCT, pCS->LinearCT2, sizeof(AT_RANK) * pCS->nLenLinearCT )*/)
            {
                nErrorCode |= WARN_FAILED_STEREO;
            }
        }
        else
        {
            nErrorCode |= ERR_NO_CANON_RESULTS;
            goto exit_function;
        }

        /*
            Isotopic part
        */

        pCanonOrd = pCS->nLenCanonOrdIsotopicStereo > 0 ? pCS->nCanonOrdIsotopicStereo :
            pCS->nLenCanonOrdIsotopic > 0 ? pCS->nCanonOrdIsotopic : NULL;
        pCanonRank = pCanonRankAtoms;

        if (pCanonOrd && pCanonRank)
        {
            for (i = 0; i < num_at_tg; i++)
            {
                pCanonRank[pCanonOrd[i]] = (AT_NUMB)(i + 1);
            }
            ret = UpdateFullLinearCT(num_atoms, num_at_tg, at, pCanonRank, pCanonOrd, pCS,
                pCG,
                0);
            if (ret /*|| memcmp(pCS->LinearCT, pCS->LinearCT2, sizeof(AT_RANK) * pCS->nLenLinearCT )*/)
            {
                nErrorCode |= (pCS->nLenCanonOrdIsotopicStereo ? WARN_FAILED_ISOTOPIC_STEREO : WARN_FAILED_ISOTOPIC);
            }
        }

    exit_function:
        if (pCanonRankAtoms)
        {
            inchi_free(pCanonRankAtoms);
        }

        if (nErrorCode)
        {
            return CT_CANON_ERR;
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: CheckCanonNumberingCorrectness

    let count = i64::from(num_at_tg)
        .checked_add(1)
        .and_then(|value| u64::try_from(value).ok())
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let ranks = match inchi_calloc::<AT_NUMB>(heap, count, 2) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
        Err(error) => return Err(error),
    };
    let computation = (|| -> Result<i32, SourceHeapError> {
        let mut error_code = 0_i32;
        let total =
            usize::try_from(num_at_tg).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;

        let nonisotopic_order = if pCS.nLenCanonOrdStereo > 0 {
            pCS.nCanonOrdStereo
        } else if pCS.nLenCanonOrd > 0 {
            pCS.nCanonOrd
        } else {
            SourceMutPointer::null()
        };
        if nonisotopic_order.is_null() || ranks.is_null() {
            error_code |= ERR_NO_CANON_RESULTS as i32;
        } else {
            for i in 0..total {
                let order = heap.slice(nonisotopic_order.as_const())?[i];
                heap.slice_mut(ranks)?[usize::from(order)] = (i as AT_NUMB).wrapping_add(1);
            }
            if UpdateFullLinearCT(
                heap,
                num_atoms,
                num_at_tg,
                at,
                ranks,
                nonisotopic_order,
                pCS,
                pCG,
                0,
            )? != 0
            {
                error_code |= WARN_FAILED_STEREO as i32;
            }

            let isotopic_order = if pCS.nLenCanonOrdIsotopicStereo > 0 {
                pCS.nCanonOrdIsotopicStereo
            } else if pCS.nLenCanonOrdIsotopic > 0 {
                pCS.nCanonOrdIsotopic
            } else {
                SourceMutPointer::null()
            };
            if !isotopic_order.is_null() {
                for i in 0..total {
                    let order = heap.slice(isotopic_order.as_const())?[i];
                    heap.slice_mut(ranks)?[usize::from(order)] = (i as AT_NUMB).wrapping_add(1);
                }
                if UpdateFullLinearCT(
                    heap,
                    num_atoms,
                    num_at_tg,
                    at,
                    ranks,
                    isotopic_order,
                    pCS,
                    pCG,
                    0,
                )? != 0
                {
                    error_code |= if pCS.nLenCanonOrdIsotopicStereo != 0 {
                        WARN_FAILED_ISOTOPIC_STEREO as i32
                    } else {
                        WARN_FAILED_ISOTOPIC as i32
                    };
                }
            }
        }
        Ok(if error_code != 0 { CT_CANON_ERR } else { 0 })
    })();
    inchi_free(heap, ranks)?;
    computation
}

#[allow(non_snake_case)]
pub(crate) fn GetElementAndCount(
    heap: &mut SourceHeap,
    formula: &mut SourceConstPointer<i8>,
    element: &mut [i8],
    count: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:173 GetElementAndCount
    // INCHI✔️❌: int GetElementAndCount(const char** f, char* szEl, int* count)
    // INCHI✔️❌: {
    // INCHI✔️❌:     const char* p = *f;
    // INCHI✔️❌:     char* q;
    // INCHI✔️❌:     int   i = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: fixing oss-fuzz issue #37224 */
    // INCHI✔️❌:     if (p && *p)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (isupper(UCINT * p))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szEl[i++] = *p++;
    // INCHI✔️❌:             if (*p && islower(UCINT * p))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 szEl[i++] = *p++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             szEl[i] = '\0';
    // INCHI✔️❌:             if (1 == i && szEl[0] == 'C')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 szEl[0] = 'A'; /*  less than any element: */
    // INCHI✔️❌:                 /*  carbon-containing compounds should be first */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (*p && isdigit(UCINT * p))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *count = strtol(p, &q, 10);
    // INCHI✔️❌:                 p = q;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *count = 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             *f = p; /*  next element; */
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return -1; /*  not a chemical formula */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* v. 1.06 Changed "Zz" to "Zzz" as "Zz" is valid symbol now */
    // INCHI✔️❌:     strcpy(szEl, "Zzz");
    // INCHI✔️❌:     /*strcpy( szEl, "Zz" );*/
    // INCHI✔️❌:     /*  zero termination 'element' is larger than any other element */
    // INCHI✔️❌:     *count = 99999;         /* zero termination 'element count' is larger than any other count */
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetElementAndCount

    if !formula.is_null() {
        let bytes = heap.slice(*formula)?;
        let length = bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if length > 0 {
            let first = bytes[0] as u8;
            if first.is_ascii_uppercase() {
                let mut consumed = 1_usize;
                let mut symbol = [0_i8; 3];
                symbol[0] = bytes[0];
                if consumed < length && (bytes[consumed] as u8).is_ascii_lowercase() {
                    symbol[1] = bytes[consumed];
                    consumed += 1;
                }
                if consumed == 1 && symbol[0] as u8 == b'C' {
                    symbol[0] = b'A' as i8;
                }
                let has_digits = consumed < length && (bytes[consumed] as u8).is_ascii_digit();
                let symbol_length = consumed + 1;
                element
                    .get_mut(..symbol_length)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .copy_from_slice(&symbol[..symbol_length]);
                let number = (*formula).offset(
                    i64::try_from(consumed).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                )?;
                if has_digits {
                    let mut end = SourceConstPointer::null();
                    *count = inchi_strtol(heap, number, Some(&mut end), 10)? as i32;
                    *formula = end;
                } else {
                    *count = 1;
                    *formula = number;
                }
                return Ok(1);
            }
            return Ok(-1);
        }
    }
    element
        .get_mut(..4)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .copy_from_slice(&[b'Z' as i8, b'z' as i8, b'z' as i8, 0]);
    *count = 99999;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompareHillFormulas(
    heap: &mut SourceHeap,
    mut formula1: SourceConstPointer<i8>,
    mut formula2: SourceConstPointer<i8>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:241 CompareHillFormulas
    // INCHI✔️❌: int CompareHillFormulas(const char* f1, const char* f2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szEl1[4], szEl2[4];
    // INCHI✔️❌:     int  count1, count2, ret1, ret2, ret;
    // INCHI✔️❌:
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret1 = GetElementAndCount(&f1, szEl1, &count1);
    // INCHI✔️❌:         ret2 = GetElementAndCount(&f2, szEl2, &count2);
    // INCHI✔️❌:         if (0 <= ret1 && 0 <= ret2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = strcmp(szEl1, szEl2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /*  lexicographic order, string termination > any character */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = count2 - count1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /*  inverse atom count order */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0; /*  program error <BRKPT> */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } while (0 < ret1 && 0 < ret2);
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareHillFormulas

    loop {
        let mut element1 = [0_i8; 4];
        let mut element2 = [0_i8; 4];
        let mut count1 = 0;
        let mut count2 = 0;
        let ret1 = GetElementAndCount(heap, &mut formula1, &mut element1, &mut count1)?;
        let ret2 = GetElementAndCount(heap, &mut formula2, &mut element2, &mut count2)?;
        if ret1 < 0 || ret2 < 0 {
            return Ok(0);
        }
        for index in 0..element1.len() {
            let ret = i32::from(element1[index] as u8) - i32::from(element2[index] as u8);
            if ret != 0 {
                return Ok(ret);
            }
            if element1[index] == 0 {
                break;
            }
        }
        let ret = count2.wrapping_sub(count1);
        if ret != 0 {
            return Ok(ret);
        }
        if ret1 <= 0 || ret2 <= 0 {
            return Ok(0);
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn CompareHillFormulasNoH(
    heap: &mut SourceHeap,
    mut formula1: SourceConstPointer<i8>,
    mut formula2: SourceConstPointer<i8>,
    num_h1: &mut i32,
    num_h2: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:272 CompareHillFormulasNoH
    // INCHI✔️❌: int CompareHillFormulasNoH(const char* f1,
    // INCHI✔️❌:     const char* f2,
    // INCHI✔️❌:     int* num_H1,
    // INCHI✔️❌:     int* num_H2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szEl1[4], szEl2[4];
    // INCHI✔️❌:     int  count1, count2, ret1, ret2, ret;
    // INCHI✔️❌:
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret1 = GetElementAndCount(&f1, szEl1, &count1);
    // INCHI✔️❌:         if (0 < ret1 && szEl1[0] == 'H' && !szEl1[1])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *num_H1 += count1;
    // INCHI✔️❌:             ret1 = GetElementAndCount(&f1, szEl1, &count1);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         ret2 = GetElementAndCount(&f2, szEl2, &count2);
    // INCHI✔️❌:         if (0 < ret2 && szEl2[0] == 'H' && !szEl2[1])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *num_H2 += count2;
    // INCHI✔️❌:             ret2 = GetElementAndCount(&f2, szEl2, &count2);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (0 <= ret1 && 0 <= ret2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = strcmp(szEl1, szEl2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /*  lexicographic order, string termination > any character */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = count2 - count1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /*  inverse atom count order */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0; /*  program error <BRKPT> */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     } while (0 < ret1 && 0 < ret2);
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareHillFormulasNoH

    loop {
        let mut element1 = [0_i8; 4];
        let mut element2 = [0_i8; 4];
        let mut count1 = 0;
        let mut count2 = 0;
        let mut ret1 = GetElementAndCount(heap, &mut formula1, &mut element1, &mut count1)?;
        if ret1 > 0 && element1[0] as u8 == b'H' && element1[1] == 0 {
            *num_h1 = num_h1
                .checked_add(count1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            ret1 = GetElementAndCount(heap, &mut formula1, &mut element1, &mut count1)?;
        }
        let mut ret2 = GetElementAndCount(heap, &mut formula2, &mut element2, &mut count2)?;
        if ret2 > 0 && element2[0] as u8 == b'H' && element2[1] == 0 {
            *num_h2 = num_h2
                .checked_add(count2)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            ret2 = GetElementAndCount(heap, &mut formula2, &mut element2, &mut count2)?;
        }
        if ret1 < 0 || ret2 < 0 {
            return Ok(0);
        }
        let mut ret = 0;
        for index in 0..4 {
            ret = i32::from(element1[index] as u8) - i32::from(element2[index] as u8);
            if ret != 0 || element1[index] == 0 {
                break;
            }
        }
        if ret != 0 {
            return Ok(ret);
        }
        let ret = count2 - count1;
        if ret != 0 {
            return Ok(ret);
        }
        if ret1 <= 0 || ret2 <= 0 {
            return Ok(0);
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn CompareTautNonIsoPartOfINChI(
    heap: &SourceHeap,
    inchi1: &INChI,
    inchi2: &INChI,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:316 CompareTautNonIsoPartOfINChI
    // INCHI✔️❌: int CompareTautNonIsoPartOfINChI(const INChI* i1, const INChI* i2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int len1, len2, ret, i;
    // INCHI✔️❌:
    // INCHI✔️❌:     len1 = i1->lenTautomer > 0 && i1->nTautomer[0] ? i1->lenTautomer : 0;
    // INCHI✔️❌:     len2 = i2->lenTautomer > 0 && i2->nTautomer[0] ? i2->lenTautomer : 0;
    // INCHI✔️❌:     if ((ret = len2 - len1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 0; i < len1; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = (int)i2->nTautomer[i] - (int)i1->nTautomer[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareTautNonIsoPartOfINChI

    let len1 = if inchi1.lenTautomer > 0 && heap.slice(inchi1.nTautomer.as_const())?[0] != 0 {
        inchi1.lenTautomer
    } else {
        0
    };
    let len2 = if inchi2.lenTautomer > 0 && heap.slice(inchi2.nTautomer.as_const())?[0] != 0 {
        inchi2.lenTautomer
    } else {
        0
    };
    let ret = len2 - len1;
    if ret != 0 {
        return Ok(ret);
    }
    if len1 > 0 {
        let count = usize::try_from(len1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let tautomer1 = heap
            .slice(inchi1.nTautomer.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let tautomer2 = heap
            .slice(inchi2.nTautomer.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for index in 0..count {
            let ret = i32::from(tautomer2[index]) - i32::from(tautomer1[index]);
            if ret != 0 {
                return Ok(ret);
            }
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompINChITautVsNonTaut(
    heap: &mut SourceHeap,
    p1: &INCHI_SORT,
    p2: &INCHI_SORT,
    compare_isotopic: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:341 CompINChITautVsNonTaut
    // INCHI✔️❌: int CompINChITautVsNonTaut(const INCHI_SORT* p1,
    // INCHI✔️❌:     const INCHI_SORT* p2,
    // INCHI✔️❌:     int bCompareIsotopic)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret, num, i, num_H1, num_H2;
    // INCHI✔️❌:
    // INCHI✔️❌:     const INChI* i1 = NULL; /* Mobile-H layers in Mobile-H sorting order */
    // INCHI✔️❌:     const INChI* i2 = NULL; /* Fixed-H  layers in Fixed-H  sorting order */
    // INCHI✔️❌:
    // INCHI✔️❌:     int   n1;               /* TAUT_YES if tautomeric i1 exists, otherwise TAUT_NON */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* INChI_Stereo *Stereo1, *Stereo2; */
    // INCHI✔️❌:
    // INCHI✔️❌:     n1 = (p1->pINChI[TAUT_YES] && p1->pINChI[TAUT_YES]->nNumberOfAtoms) ? TAUT_YES : TAUT_NON;
    // INCHI✔️❌:
    // INCHI✔️❌:     i1 = p1->pINChI[n1];
    // INCHI✔️❌:     i2 = (n1 == TAUT_YES && p2->pINChI[TAUT_NON] &&
    // INCHI✔️❌:         p2->pINChI[TAUT_NON]->nNumberOfAtoms) ? p2->pINChI[TAUT_NON] : (const INChI*)NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* non-deleted-non-empty < deleted < empty */
    // INCHI✔️❌:     if (i1 && !i2)
    // INCHI✔️❌:         return 0;   /* non-empty is the smallest (first) */
    // INCHI✔️❌:     if (!i1 && i2)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     if (!i1 && !i2)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     if (i1->bDeleted)
    // INCHI✔️❌:         return 1;    /* deleted is the largest (last) among non-empty */
    // INCHI✔️❌:     if (i2->bDeleted)
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1->nNumberOfAtoms > 0 && !i2->nNumberOfAtoms)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* i2 = i2;         djb-rwth: an obviously useless statement */
    // INCHI✔️❌:
    // INCHI✔️❌:     num_H1 = num_H2 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* do not compare terminal H */
    // INCHI✔️❌:     if ((ret = CompareHillFormulasNoH(i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;  /* lexicographic order except the shorter one is greater (last): CH2O < CH2; C3XX < C2XX */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         compare non-isotopic non-tautomeric part
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* compare number of atoms (excluding terminal H) */
    // INCHI✔️❌:     if ((ret = i2->nNumberOfAtoms - i1->nNumberOfAtoms)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return ret; /*  more atoms first */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  compare elements  (excluding terminal H) */
    // INCHI✔️❌:     num = i1->nNumberOfAtoms;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     { /* should always be equal if Hill formulas are same */
    // INCHI✔️❌:         if ((ret = (int)i2->nAtom[i] - (int)i1->nAtom[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret; /* greater periodic number first */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         compare connection tables
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if ((ret = i2->lenConnTable - i1->lenConnTable)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return ret; /* longer connection table first */
    // INCHI✔️❌:     num = i2->lenConnTable;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = (int)i2->nConnTable[i] - (int)i1->nConnTable[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret; /* greater connection table first */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:       compare total number of H (inverse: H3 < H2 )
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if ((ret = num_H2 - num_H1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     /*
    // INCHI✔️❌:       compare non-tautomeric num_H: N < NH3 < NH2 < NH
    // INCHI✔️❌:     */
    // INCHI✔️❌:     num = i1->nNumberOfAtoms;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i2->nNum_H[i] != i1->nNum_H[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return !i2->nNum_H[i] ? 1 :  /* no H first */
    // INCHI✔️❌:                 !i1->nNum_H[i] ? -1 :
    // INCHI✔️❌:                 (int)i2->nNum_H[i] - (int)i1->nNum_H[i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         compare non-isotopic tautomeric part
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if ((ret = CompareTautNonIsoPartOfINChI(i1, i2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( ret = i2->lenTautomer - i1->lenTautomer )
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     num = inchi_min( i2->lenTautomer, i1->lenTautomer );
    // INCHI✔️❌:     for ( i = 0; i < num; i ++ ) {
    // INCHI✔️❌:         if ( ret = (int)i2->nTautomer[i] - (int)i1->nTautomer[i] )
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         at this point both components are either tautomeric
    // INCHI✔️❌:         or non-tautomeric
    // INCHI✔️❌:      */
    // INCHI✔️❌:
    // INCHI✔️❌:      /*
    // INCHI✔️❌:          non-tautomeric "fixed H" specific
    // INCHI✔️❌:      */
    // INCHI✔️❌:     if ( /*TAUT_NON == bTaut && (i2 &&*/ i2->nNum_H_fixed) /* djb-rwth: fixing coverity ID #499493 */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* first, compare non-tautomeric chem. formulas -- they may be different */
    // INCHI✔️❌:         /* secondly, compare fixed-H distribution */
    // INCHI✔️❌:         if (i2->nNum_H_fixed)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num = i2->nNumberOfAtoms;
    // INCHI✔️❌:             for (i = 0; i < num; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i2->nNum_H_fixed[i] != 0)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 1;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         compare non-isotopic stereo
    // INCHI✔️❌:     */
    // INCHI✔️❌:     ret = CompareInchiStereo(i1->Stereo, i1->nFlags, i2->Stereo, i2->nFlags);
    // INCHI✔️❌:     if (ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         do not switch back to tautomeric i1, i2
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* -- how to switch back --
    // INCHI✔️❌:     if ( i1t ) {
    // INCHI✔️❌:         i1  = i1t;
    // INCHI✔️❌:         i1t = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ( i2t ) {
    // INCHI✔️❌:         i2  = i2t;
    // INCHI✔️❌:         i2t = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:          compare isotopic non-tautomeric part
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (bCompareIsotopic)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = i2->nNumberOfIsotopicAtoms - i1->nNumberOfIsotopicAtoms)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         num = i1->nNumberOfIsotopicAtoms;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  compare isotopic atoms */
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nAtomNumber - (int)i1->IsotopicAtom[i].nAtomNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nIsoDifference - (int)i1->IsotopicAtom[i].nIsoDifference)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* compare isotopic H */
    // INCHI✔️❌:         /* if tautomeric comparison mode then here are compared only non-tautomeric H */
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_T - (int)i1->IsotopicAtom[i].nNum_T)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_D - (int)i1->IsotopicAtom[i].nNum_D)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_H - (int)i1->IsotopicAtom[i].nNum_H)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* compare isotopic tautomeric part */
    // INCHI✔️❌:         if ((ret = i2->nNumberOfIsotopicTGroups || i1->nNumberOfIsotopicTGroups)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         num = i1->nNumberOfIsotopicTGroups;
    // INCHI✔️❌:         for ( i = 0; i < num; i ++ ) {
    // INCHI✔️❌:             if ( ret = (int)i2->IsotopicTGroup[i].nTGroupNumber - (int)i1->IsotopicTGroup[i].nTGroupNumber )
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ( ret = (int)i2->IsotopicTGroup[i].nNum_T - (int)i1->IsotopicTGroup[i].nNum_T )
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ( ret = (int)i2->IsotopicTGroup[i].nNum_D - (int)i1->IsotopicTGroup[i].nNum_D )
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             if ( ret = (int)i2->IsotopicTGroup[i].nNum_H - (int)i1->IsotopicTGroup[i].nNum_H )
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* compare isotopic stereo */
    // INCHI✔️❌:         ret = CompareInchiStereo(i1->StereoIsotopic, i1->nFlags,
    // INCHI✔️❌:             i2->StereoIsotopic, i2->nFlags);
    // INCHI✔️❌:         if (ret)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         compare charges: non-charged first, then in order of
    // INCHI✔️❌:         ascending charges (negative first)
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i2->nTotalCharge && i1->nTotalCharge)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  both are charged; smaller charges first */
    // INCHI✔️❌:         ret = (int)i1->nTotalCharge - (int)i2->nTotalCharge;
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((ret = (i1->nTotalCharge ? 1 : 0) - (i2->nTotalCharge ? 1 : 0))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  only one is charged; uncharged first */
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* stable sort */
    // INCHI✔️❌:     /*ret = p1->ord_number - p2->ord_number;*/
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompINChITautVsNonTaut

    let taut_yes = TAUT_YES as usize;
    let taut_non = TAUT_NON as usize;
    let p1_taut = p1.pINChI[taut_yes];
    let n1 = if !p1_taut.is_null() && heap.slice(p1_taut.as_const())?[0].nNumberOfAtoms != 0 {
        taut_yes
    } else {
        taut_non
    };
    let i1_pointer = p1.pINChI[n1];
    let p2_non = p2.pINChI[taut_non];
    let i2_pointer = if n1 == taut_yes
        && !p2_non.is_null()
        && heap.slice(p2_non.as_const())?[0].nNumberOfAtoms != 0
    {
        p2_non
    } else {
        crate::source_types::SourceMutPointer::null()
    };
    if i1_pointer.is_null() || i2_pointer.is_null() {
        return Ok(0);
    }
    let i1 = heap.slice(i1_pointer.as_const())?[0].clone();
    let i2 = heap.slice(i2_pointer.as_const())?[0].clone();
    if i1.bDeleted != 0 {
        return Ok(1);
    }
    if i2.bDeleted != 0 {
        return Ok(-1);
    }
    if i1.nNumberOfAtoms > 0 && i2.nNumberOfAtoms == 0 {
        return Ok(0);
    }

    let mut num_h1 = 0;
    let mut num_h2 = 0;
    let ret = CompareHillFormulasNoH(
        heap,
        i1.szHillFormula.as_const(),
        i2.szHillFormula.as_const(),
        &mut num_h1,
        &mut num_h2,
    )?;
    if ret != 0 {
        return Ok(ret);
    }
    let ret = i2.nNumberOfAtoms.wrapping_sub(i1.nNumberOfAtoms);
    if ret != 0 {
        return Ok(ret);
    }
    let atom_count = usize::try_from(i1.nNumberOfAtoms.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let atoms1 = heap
        .slice(i1.nAtom.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let atoms2 = heap
        .slice(i2.nAtom.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for index in 0..atom_count {
        let ret = i32::from(atoms2[index]) - i32::from(atoms1[index]);
        if ret != 0 {
            return Ok(ret);
        }
    }
    let ret = i2.lenConnTable.wrapping_sub(i1.lenConnTable);
    if ret != 0 {
        return Ok(ret);
    }
    let connection_count = usize::try_from(i2.lenConnTable.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if connection_count > 0 {
        let connections1 = heap
            .slice(i1.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let connections2 = heap
            .slice(i2.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for index in 0..connection_count {
            let ret = i32::from(connections2[index]) - i32::from(connections1[index]);
            if ret != 0 {
                return Ok(ret);
            }
        }
    }
    let ret = num_h2.wrapping_sub(num_h1);
    if ret != 0 {
        return Ok(ret);
    }
    let hydrogens1 = heap
        .slice(i1.nNum_H.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let hydrogens2 = heap
        .slice(i2.nNum_H.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for index in 0..atom_count {
        if hydrogens2[index] != hydrogens1[index] {
            return Ok(if hydrogens2[index] == 0 {
                1
            } else if hydrogens1[index] == 0 {
                -1
            } else {
                i32::from(hydrogens2[index]) - i32::from(hydrogens1[index])
            });
        }
    }
    let ret = CompareTautNonIsoPartOfINChI(heap, &i1, &i2)?;
    if ret != 0 {
        return Ok(ret);
    }
    if !i2.nNum_H_fixed.is_null() {
        let fixed = heap
            .slice(i2.nNum_H_fixed.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if fixed.iter().any(|value| *value != 0) {
            return Ok(1);
        }
    }
    let stereo1 = if i1.Stereo.is_null() {
        None
    } else {
        Some(heap.slice(i1.Stereo.as_const())?[0].clone())
    };
    let stereo2 = if i2.Stereo.is_null() {
        None
    } else {
        Some(heap.slice(i2.Stereo.as_const())?[0].clone())
    };
    let ret = CompareInchiStereo(
        heap,
        stereo1.as_ref(),
        i1.nFlags,
        stereo2.as_ref(),
        i2.nFlags,
    )?;
    if ret != 0 {
        return Ok(ret);
    }

    if compare_isotopic != 0 {
        let ret = i2
            .nNumberOfIsotopicAtoms
            .wrapping_sub(i1.nNumberOfIsotopicAtoms);
        if ret != 0 {
            return Ok(ret);
        }
        let isotope_count = usize::try_from(i1.nNumberOfIsotopicAtoms.max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if isotope_count > 0 {
            let isotopes1 = heap
                .slice(i1.IsotopicAtom.as_const())?
                .get(..isotope_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let isotopes2 = heap
                .slice(i2.IsotopicAtom.as_const())?
                .get(..isotope_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in 0..isotope_count {
                for ret in [
                    i32::from(isotopes2[index].nAtomNumber)
                        - i32::from(isotopes1[index].nAtomNumber),
                    i32::from(isotopes2[index].nIsoDifference)
                        - i32::from(isotopes1[index].nIsoDifference),
                ] {
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
            for index in 0..isotope_count {
                for ret in [
                    i32::from(isotopes2[index].nNum_T) - i32::from(isotopes1[index].nNum_T),
                    i32::from(isotopes2[index].nNum_D) - i32::from(isotopes1[index].nNum_D),
                    i32::from(isotopes2[index].nNum_H) - i32::from(isotopes1[index].nNum_H),
                ] {
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
        }
        if i2.nNumberOfIsotopicTGroups != 0 || i1.nNumberOfIsotopicTGroups != 0 {
            return Ok(1);
        }
        let stereo1 = if i1.StereoIsotopic.is_null() {
            None
        } else {
            Some(heap.slice(i1.StereoIsotopic.as_const())?[0].clone())
        };
        let stereo2 = if i2.StereoIsotopic.is_null() {
            None
        } else {
            Some(heap.slice(i2.StereoIsotopic.as_const())?[0].clone())
        };
        let ret = CompareInchiStereo(
            heap,
            stereo1.as_ref(),
            i1.nFlags,
            stereo2.as_ref(),
            i2.nFlags,
        )?;
        if ret != 0 {
            return Ok(ret);
        }
    }
    if i2.nTotalCharge != 0 && i1.nTotalCharge != 0 {
        return Ok(i1.nTotalCharge.wrapping_sub(i2.nTotalCharge));
    }
    Ok(i32::from(i1.nTotalCharge != 0) - i32::from(i2.nTotalCharge != 0))
}

#[allow(non_snake_case)]
pub(crate) fn GetSp3RelRacAbs(inchi: Option<&INChI>, stereo: Option<&INChI_Stereo>) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:593 GetSp3RelRacAbs
    // INCHI✔️❌: int GetSp3RelRacAbs(const INChI* pINChI, INChI_Stereo* Stereo)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int nRet = SP3_NONE;
    // INCHI✔️❌:     if (pINChI && !pINChI->bDeleted && Stereo && 0 < Stereo->nNumberOfStereoCenters)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (0 != Stereo->nCompInv2Abs)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (pINChI->nFlags & INCHI_FLAG_REL_STEREO)
    // INCHI✔️❌:             {
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:                 if (1 < Stereo->nNumberOfStereoCenters)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nRet = SP3_REL;
    // INCHI✔️❌:                 }
    // INCHI✔️❌: #else
    // INCHI✔️❌:                 nRet = SP3_REL;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:                 if (pINChI->nFlags & INCHI_FLAG_RAC_STEREO)
    // INCHI✔️❌:                 {
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:                     if (1 < Stereo->nNumberOfStereoCenters)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         nRet = SP3_REL;
    // INCHI✔️❌:                     }
    // INCHI✔️❌: #else
    // INCHI✔️❌:                     nRet = SP3_RAC;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     nRet = SP3_ABS;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:             if (!((pINChI->nFlags & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) && 1 == Stereo->nNumberOfStereoCenters))
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nRet = SP3_ONLY; /*  SP3_NONE if relative stereo and 1 stereocenter */
    // INCHI✔️❌:             }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return nRet;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetSp3RelRacAbs

    const SP3_NONE: i32 = 0;
    const SP3_ONLY: i32 = 1;
    const SP3_ABS: i32 = 2;
    const SP3_REL: i32 = 4;
    const SP3_RAC: i32 = 8;

    let (Some(inchi), Some(stereo)) = (inchi, stereo) else {
        return SP3_NONE;
    };
    if inchi.bDeleted != 0 || stereo.nNumberOfStereoCenters <= 0 {
        return SP3_NONE;
    }
    if stereo.nCompInv2Abs != 0 {
        if inchi.nFlags & INCHI_MODE::from(INCHI_FLAG_REL_STEREO) != 0 {
            return SP3_REL;
        }
        if inchi.nFlags & INCHI_MODE::from(INCHI_FLAG_RAC_STEREO) != 0 {
            return SP3_RAC;
        }
        return SP3_ABS;
    }
    SP3_ONLY
}

#[allow(non_snake_case)]
pub(crate) fn CompINChILayers(
    heap: &mut SourceHeap,
    p1: &INCHI_SORT,
    p2: &INCHI_SORT,
    difference_segments: &mut [[i8; 11]; 4],
    fix_transposed_charge_bug: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:645 CompINChILayers
    // INCHI✔️❌: int CompINChILayers(const INCHI_SORT* p1,
    // INCHI✔️❌:     const INCHI_SORT* p2,
    // INCHI✔️❌:     char sDifSegs[][DIFS_LENGTH],
    // INCHI✔️❌:     int bFixTranspChargeBug)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret = 0, num, i, num_H1, num_H2;
    // INCHI✔️❌:
    // INCHI✔️❌:     const INChI* i1 = NULL; /* Mobile-H layers in Mobile-H sorting order */
    // INCHI✔️❌:     const INChI* i2 = NULL; /* Fixed-H  layers in Fixed-H  sorting order */
    // INCHI✔️❌:
    // INCHI✔️❌:     int   n1;               /* TAUT_YES if tautomeric i1 exists, otherwise TAUT_NON */
    // INCHI✔️❌:
    // INCHI✔️❌:     INChI_Stereo* Stereo1, * Stereo2;
    // INCHI✔️❌:     INChI_Stereo* IsoStereo1, * IsoStereo2;
    // INCHI✔️❌:     int bRelRac[DIFL_LENGTH];
    // INCHI✔️❌:     char* psDifSegs;
    // INCHI✔️❌:
    // INCHI✔️❌:     n1 = (p1->pINChI[TAUT_YES] && p1->pINChI[TAUT_YES]->nNumberOfAtoms) ? TAUT_YES : TAUT_NON;
    // INCHI✔️❌:
    // INCHI✔️❌:     i1 = p1->pINChI[n1];
    // INCHI✔️❌:     i2 = (n1 == TAUT_YES && p2->pINChI[TAUT_NON] &&
    // INCHI✔️❌:         p2->pINChI[TAUT_NON]->nNumberOfAtoms) ? p2->pINChI[TAUT_NON] : (const INChI*)NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     num_H1 = num_H2 = 0;
    // INCHI✔️❌:     memset(bRelRac, DIFV_BOTH_EMPTY, sizeof(bRelRac)); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /f    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     if (i1 && !i1->bDeleted && i1->szHillFormula && i1->szHillFormula[0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         if (i2 && !i2->bDeleted && i2->szHillFormula && i2->szHillFormula[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!CompareHillFormulasNoH(i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2) &&
    // INCHI✔️❌:                 num_H1 == num_H2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[DIFL_F][DIFS_f_FORMULA] |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         if (i2 && !i2->bDeleted && i2->szHillFormula && i2->szHillFormula[0])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[DIFL_F][DIFS_f_FORMULA] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /c    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     if (i1 && !i1->bDeleted && i1->lenConnTable > 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_f_FORMULA] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /h    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /* M: H atoms */
    // INCHI✔️❌:     if (i1 && !i1->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_H1 = (i1->lenTautomer > 0 && i1->nTautomer && i1->nTautomer[0]) ? 1 : 0; /* number of t-groups */
    // INCHI✔️❌:         if (!num_H1 && i1->nNum_H)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:             { /* immobile H */
    // INCHI✔️❌:                 if (i1->nNum_H[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num_H1 = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_h_H_ATOMS] |= num_H1 ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_h_H_ATOMS] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* F: fixed mobile H */
    // INCHI✔️❌:     if (i2 && !i2->bDeleted && i2->nNum_H_fixed)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         num_H2 = 0;
    // INCHI✔️❌:         if (i1 && !i1->bDeleted)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (i = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i2->nNum_H_fixed[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num_H2 = 1;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         sDifSegs[DIFL_F][DIFS_h_H_ATOMS] |= num_H2 ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_F][DIFS_h_H_ATOMS] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* MI: exchangable isotopic H: see OutputINChI1(), num_iso_H[] */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /q    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_F][DIFS_q_CHARGE];
    // INCHI✔️❌:     if (i1 && !i1->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i1->nTotalCharge)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[DIFL_M][DIFS_q_CHARGE] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[DIFL_M][DIFS_q_CHARGE] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i2 && !i2->bDeleted)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i1->nTotalCharge)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i1->nTotalCharge == i2->nTotalCharge)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (i2->nTotalCharge)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         *psDifSegs |= DIFV_IS_EMPTY;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i2->nTotalCharge)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!i2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bFixTranspChargeBug == 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* bug explanation:
    // INCHI✔️❌:
    // INCHI✔️❌:                     component #1 is tautomeric, component #2 is not
    // INCHI✔️❌:                     Mobile-H(#2) > Mobile-H(#1)
    // INCHI✔️❌:                     Fixed-H(#2) = Mobile-H(#2) < Fixed-H(#1)
    // INCHI✔️❌:
    // INCHI✔️❌:                     Layer       first_charge   second_charge
    // INCHI✔️❌:
    // INCHI✔️❌:                     Mobile-H    0    (comp#1)  -1 (comp#2)
    // INCHI✔️❌:                     Fixed-H     none (comp#2)  -1 (comp#1)
    // INCHI✔️❌:
    // INCHI✔️❌:                     v1.01 charge compared decided that charge layers are same and omitted Fixed-H /q layer
    // INCHI✔️❌:
    // INCHI✔️❌:                     Solution: when component permutation is detected AND fixed-H component does not exist,
    // INCHI✔️❌:                     compare Mobile-H charge [0 (comp#1) in the example] to the charge of Mobile-H [-1 (comp#2)]
    // INCHI✔️❌:                     of the component that has none Fixed-H charge
    // INCHI✔️❌:                     */
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* Fixed-H i2 is empty because Fixed-H struct is same as Mobile-H */
    // INCHI✔️❌:                     if (p1->ord_number != p2->ord_number && /* component order in Fixed-H is different from Mobile-H */
    // INCHI✔️❌:                         n1 == TAUT_YES && p2->pINChI[TAUT_YES] && !p2->pINChI[TAUT_YES]->bDeleted &&
    // INCHI✔️❌:                         p2->pINChI[TAUT_YES]->nNumberOfAtoms)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         int i2_nTotalCharge = p2->pINChI[TAUT_YES]->nTotalCharge;
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (i1->nTotalCharge)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (i1->nTotalCharge == i2_nTotalCharge)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (i2_nTotalCharge)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                                 else
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     *psDifSegs |= DIFV_IS_EMPTY;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (i2_nTotalCharge)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         *psDifSegs |= i1->nTotalCharge ? DIFV_EQL2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else /* if (bFixTranspChargeBug==1) */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     *psDifSegs |= i1->nTotalCharge ? DIFV_EQL2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else /* if ( !i2 ) { */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* i2 && i2->bDeleted */
    // INCHI✔️❌:                 *psDifSegs |= i1->nTotalCharge ? DIFV_IS_EMPTY : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_q_CHARGE] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         if (i2 && !i2->bDeleted)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i2->nTotalCharge)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sDifSegs[DIFL_F][DIFS_q_CHARGE] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 sDifSegs[DIFL_F][DIFS_q_CHARGE] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************** stereo *****************/
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1 && !i1->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Stereo1 = i1->Stereo;
    // INCHI✔️❌:         IsoStereo1 = i1->StereoIsotopic;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Stereo1 = NULL;
    // INCHI✔️❌:         IsoStereo1 = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i2 && !i2->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Stereo2 = i2->Stereo;
    // INCHI✔️❌:         IsoStereo2 = i2->StereoIsotopic;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         Stereo2 = NULL;
    // INCHI✔️❌:         IsoStereo2 = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /b    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /* M double bond stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_M][DIFS_b_SBONDS];
    // INCHI✔️❌:     if (Stereo1 && Stereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* F double bond stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_F][DIFS_b_SBONDS];
    // INCHI✔️❌:     if (Stereo2 && Stereo2->nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Stereo1 && Stereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (Eql_INChI_Stereo(Stereo1, EQL_SP2, Stereo2, EQL_SP2, 0))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Stereo1 && Stereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* MI double bond stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_MI][DIFS_b_SBONDS];
    // INCHI✔️❌:     if (IsoStereo1 && IsoStereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Eql_INChI_Stereo(IsoStereo1, EQL_SP2, Stereo1, EQL_SP2, 0))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Stereo1 && Stereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* FI double bond stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_FI][DIFS_b_SBONDS];
    // INCHI✔️❌:     if (IsoStereo2 && IsoStereo2->nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Eql_INChI_Stereo(IsoStereo2, EQL_SP2, Stereo2, EQL_SP2, 0))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!(Stereo1 && Stereo1->nNumberOfStereoBonds) &&
    // INCHI✔️❌:                 !(Stereo2 && Stereo2->nNumberOfStereoBonds) &&
    // INCHI✔️❌:                 Eql_INChI_Stereo(IsoStereo2, EQL_SP2, IsoStereo1, EQL_SP2, 0))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_FI_EQ_MI;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* the solution table for FI stereo,
    // INCHI✔️❌:            in case of FI stereo is empty
    // INCHI✔️❌:            E = segment is empty, NE = not empty
    // INCHI✔️❌:            +==============================+
    // INCHI✔️❌:            | M   | MI  | F   |  result    |
    // INCHI✔️❌:            +=====+=====+=====+============+
    // INCHI✔️❌:            | E   | E   | E   | both empty |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | NE  | E   | E   | both empty |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | E   | NE  | E   | is empty   |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | NE  | NE  | E   | both empty |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | E   | E   | NE  | is empty   |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | NE  | E   | NE  | is empty   |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | E   | NE  | NE  | is empty   |
    // INCHI✔️❌:            +-----+-----+-----+------------+
    // INCHI✔️❌:            | NE  | NE  | ME  | is empty   |
    // INCHI✔️❌:            +==============================+
    // INCHI✔️❌:         */
    // INCHI✔️❌:         if (Stereo2 && Stereo2->nNumberOfStereoBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:             if (IsoStereo1 && IsoStereo1->nNumberOfStereoBonds &&
    // INCHI✔️❌:                 !(Stereo1 && Stereo1->nNumberOfStereoBonds)
    // INCHI✔️❌:                 )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*==================================*/
    // INCHI✔️❌:     /*====     /t, /m, /s for M   ======*/
    // INCHI✔️❌:     /*==================================*/
    // INCHI✔️❌:     /* M sp3 stereo */
    // INCHI✔️❌:     bRelRac[DIFL_M] = GetSp3RelRacAbs(i1, Stereo1);       /* Mobile-H */
    // INCHI✔️❌:     bRelRac[DIFL_MI] = GetSp3RelRacAbs(i1, IsoStereo1);
    // INCHI✔️❌:     bRelRac[DIFL_F] = GetSp3RelRacAbs(i2, Stereo2);       /* Fixed-H */
    // INCHI✔️❌:     bRelRac[DIFL_FI] = GetSp3RelRacAbs(i2, IsoStereo2);
    // INCHI✔️❌:     if (SP3_NONE != bRelRac[DIFL_M])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_t_SATOMS] |= (bRelRac[DIFL_M] & SP3_ANY) ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_m_SP3INV] |= (bRelRac[DIFL_M] & SP3_ABS) ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_s_STYPE] |= (bRelRac[DIFL_M] & SP3_TYPE) ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_t_SATOMS] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_m_SP3INV] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         sDifSegs[DIFL_M][DIFS_s_STYPE] |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /t    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /* F sp3 stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_F][DIFS_t_SATOMS];
    // INCHI✔️❌:     if (SP3_ANY & bRelRac[DIFL_F])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Eql_INChI_Stereo(Stereo2, EQL_SP3, Stereo1, EQL_SP3, 0))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (SP3_ANY & bRelRac[DIFL_M])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* MI sp3 stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_MI][DIFS_t_SATOMS];
    // INCHI✔️❌:     if (SP3_ANY & bRelRac[DIFL_MI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Eql_INChI_Stereo(IsoStereo1, EQL_SP3, Stereo1, EQL_SP3, 0))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (SP3_ANY & bRelRac[DIFL_M])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* FI sp3 stereo */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_FI][DIFS_t_SATOMS];
    // INCHI✔️❌:     if (SP3_ANY & bRelRac[DIFL_FI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Eql_INChI_Stereo(IsoStereo2, EQL_SP3, Stereo2, EQL_SP3, 0))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!(SP3_ANY & bRelRac[DIFL_M]) &&
    // INCHI✔️❌:                 !(SP3_ANY & bRelRac[DIFL_F]) &&
    // INCHI✔️❌:                 Eql_INChI_Stereo(IsoStereo2, EQL_SP3, IsoStereo1, EQL_SP3, 0))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_FI_EQ_MI;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else /* similar to /b */
    // INCHI✔️❌:         if ((SP3_ANY & bRelRac[DIFL_F]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((SP3_ANY & bRelRac[DIFL_MI]) && !(SP3_ANY & bRelRac[DIFL_M]))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /m    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /* F sp3 abs stereo inversion */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_F][DIFS_m_SP3INV];
    // INCHI✔️❌:     if (bRelRac[DIFL_F] & SP3_ABS)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* the order of || operands below is critically important: || is not a commutative operation */
    // INCHI✔️❌:         if (!(bRelRac[DIFL_M] & SP3_ABS) || Stereo2->nCompInv2Abs != Stereo1->nCompInv2Abs)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (bRelRac[DIFL_M] & SP3_ABS)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* MI sp3 abs stereo inversion */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_MI][DIFS_m_SP3INV];
    // INCHI✔️❌:     if (SP3_ABS & bRelRac[DIFL_MI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((SP3_ABS & bRelRac[DIFL_M]) && IsoStereo1->nCompInv2Abs == Stereo1->nCompInv2Abs)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (SP3_ABS & bRelRac[DIFL_M])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* FI sp3 abs stereo inversion */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_FI][DIFS_m_SP3INV];
    // INCHI✔️❌:     if (SP3_ABS & bRelRac[DIFL_FI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((SP3_ABS & bRelRac[DIFL_F]) && IsoStereo2->nCompInv2Abs == Stereo2->nCompInv2Abs)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!(SP3_ABS & bRelRac[DIFL_M]) &&
    // INCHI✔️❌:                 !(SP3_ABS & bRelRac[DIFL_F]) &&
    // INCHI✔️❌:                 (SP3_ABS & bRelRac[DIFL_MI]) && /* make sure IsoStereo1 != NULL */
    // INCHI✔️❌:                 IsoStereo2->nCompInv2Abs == IsoStereo1->nCompInv2Abs)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_FI_EQ_MI;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* similar to /b */
    // INCHI✔️❌:         /* the order of || operands below is critically important: || is no a commutative operation */
    // INCHI✔️❌:         if ((SP3_ABS & bRelRac[DIFL_F]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((SP3_ABS & bRelRac[DIFL_MI]) && !(SP3_ABS & bRelRac[DIFL_M]))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /s    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /* F sp3 stereo type */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_F][DIFS_s_STYPE];
    // INCHI✔️❌:     if (bRelRac[DIFL_F] & SP3_TYPE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((bRelRac[DIFL_F] & SP3_TYPE) == (bRelRac[DIFL_M] & SP3_TYPE))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (bRelRac[DIFL_M] & SP3_TYPE)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* MI sp3 stereo type */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_MI][DIFS_s_STYPE];
    // INCHI✔️❌:     if (SP3_TYPE & bRelRac[DIFL_MI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((SP3_TYPE & bRelRac[DIFL_MI]) == (SP3_TYPE & bRelRac[DIFL_M]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (SP3_TYPE & bRelRac[DIFL_M])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* FI sp3 stereo type */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_FI][DIFS_s_STYPE];
    // INCHI✔️❌:     if (SP3_TYPE & bRelRac[DIFL_FI])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((SP3_TYPE & bRelRac[DIFL_FI]) == (SP3_TYPE & bRelRac[DIFL_F]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!(SP3_TYPE & bRelRac[DIFL_M]) &&
    // INCHI✔️❌:                 !(SP3_TYPE & bRelRac[DIFL_F]) &&
    // INCHI✔️❌:                 (SP3_TYPE & bRelRac[DIFL_MI]))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_FI_EQ_MI;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* similar to /b */
    // INCHI✔️❌:         /* the order of || operands below is critically important: || is not a commutative operation */
    // INCHI✔️❌:         if ((SP3_TYPE & bRelRac[DIFL_F]))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             *psDifSegs |= DIFV_EQL2PRECED; /* isotopic is missing because there is no isotopes */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((SP3_TYPE & bRelRac[DIFL_MI]) && !(SP3_TYPE & bRelRac[DIFL_M]))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= i2 ? DIFV_IS_EMPTY : DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /o    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     if (p1 && p2 && p1->ord_number != p2->ord_number)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sDifSegs[DIFL_F][DIFS_o_TRANSP] |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:     /*====     /i    ======*/
    // INCHI✔️❌:     /*=====================*/
    // INCHI✔️❌:
    // INCHI✔️❌:     /* M isotopic atoms */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_MI][DIFS_i_IATOMS];
    // INCHI✔️❌:     if (i1 && !i1->bDeleted && (i1->nNumberOfIsotopicAtoms || i1->nNumberOfIsotopicTGroups))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* F isotopic atoms */
    // INCHI✔️❌:     psDifSegs = &sDifSegs[DIFL_FI][DIFS_i_IATOMS];
    // INCHI✔️❌:     if (i2 && !i2->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i2->nNumberOfIsotopicAtoms || i2->nNumberOfIsotopicTGroups)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (!i1 || i1->bDeleted ||
    // INCHI✔️❌:                 i2->nNumberOfIsotopicAtoms != i1->nNumberOfIsotopicAtoms ||
    // INCHI✔️❌:                 i2->nNumberOfIsotopicTGroups != i1->nNumberOfIsotopicTGroups)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_NEQ2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int diff;
    // INCHI✔️❌:                 num = i1->nNumberOfIsotopicAtoms;
    // INCHI✔️❌:                 diff = 0;
    // INCHI✔️❌:                 for (i = 0; i < num; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* compare isotopic atoms */
    // INCHI✔️❌:                     if ((diff = (int)i2->IsotopicAtom[i].nAtomNumber - (int)i1->IsotopicAtom[i].nAtomNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if ((diff = (int)i2->IsotopicAtom[i].nIsoDifference - (int)i1->IsotopicAtom[i].nIsoDifference)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* compare isotopic H */
    // INCHI✔️❌:                     if ((diff = (int)i2->IsotopicAtom[i].nNum_T - (int)i1->IsotopicAtom[i].nNum_T)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if ((diff = (int)i2->IsotopicAtom[i].nNum_D - (int)i1->IsotopicAtom[i].nNum_D)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if ((diff = (int)i2->IsotopicAtom[i].nNum_H - (int)i1->IsotopicAtom[i].nNum_H)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (!diff)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     num = i1->nNumberOfIsotopicTGroups;
    // INCHI✔️❌:                     for (i = 0; i < num; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if ((diff = (int)i2->IsotopicTGroup[i].nTGroupNumber - (int)i1->IsotopicTGroup[i].nTGroupNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if ((diff = (int)i2->IsotopicTGroup[i].nNum_T - (int)i1->IsotopicTGroup[i].nNum_T)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if ((diff = (int)i2->IsotopicTGroup[i].nNum_D - (int)i1->IsotopicTGroup[i].nNum_D)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             return diff;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if ((diff = (int)i2->IsotopicTGroup[i].nNum_H - (int)i1->IsotopicTGroup[i].nNum_H)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 *psDifSegs |= diff ? DIFV_NEQ2PRECED : DIFV_FI_EQ_MI;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i1 && !i1->bDeleted && (i1->nNumberOfIsotopicAtoms || i1->nNumberOfIsotopicTGroups))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_IS_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (!i2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i1 && !i1->bDeleted && (i1->nNumberOfIsotopicAtoms || i1->nNumberOfIsotopicTGroups))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_EQL2PRECED;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *psDifSegs |= DIFV_BOTH_EMPTY;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // INCHI✔️❌:
    // END INCHI C FUNCTION: CompINChILayers

    const M: usize = tagDiffINChILayers_DIFL_M as usize;
    const MI: usize = tagDiffINChILayers_DIFL_MI as usize;
    const F: usize = tagDiffINChILayers_DIFL_F as usize;
    const FI: usize = tagDiffINChILayers_DIFL_FI as usize;
    const FORMULA: usize = tagDiffINChISegments_DIFS_f_FORMULA as usize;
    const H_ATOMS: usize = tagDiffINChISegments_DIFS_h_H_ATOMS as usize;
    const CHARGE: usize = tagDiffINChISegments_DIFS_q_CHARGE as usize;
    const SBONDS: usize = tagDiffINChISegments_DIFS_b_SBONDS as usize;
    const SATOMS: usize = tagDiffINChISegments_DIFS_t_SATOMS as usize;
    const SP3INV: usize = tagDiffINChISegments_DIFS_m_SP3INV as usize;
    const STYPE: usize = tagDiffINChISegments_DIFS_s_STYPE as usize;
    const IATOMS: usize = tagDiffINChISegments_DIFS_i_IATOMS as usize;
    const TRANSP: usize = tagDiffINChISegments_DIFS_o_TRANSP as usize;
    const BOTH_EMPTY: i8 = tagMarkDiff_DIFV_BOTH_EMPTY as i8;
    const EQUAL: i8 = tagMarkDiff_DIFV_EQL2PRECED as i8;
    const NOT_EQUAL: i8 = tagMarkDiff_DIFV_NEQ2PRECED as i8;
    const IS_EMPTY: i8 = tagMarkDiff_DIFV_IS_EMPTY as i8;
    const FI_EQ_MI: i8 = tagMarkDiff_DIFV_FI_EQ_MI as i8;
    const SP3_NONE: i32 = 0;
    const SP3_ONLY: i32 = 1;
    const SP3_ABS: i32 = 2;
    const SP3_REL: i32 = 4;
    const SP3_RAC: i32 = 8;
    const SP3_TYPE: i32 = SP3_ABS | SP3_REL | SP3_RAC;
    const SP3_ANY: i32 = SP3_TYPE | SP3_ONLY;

    let taut_yes = TAUT_YES as usize;
    let taut_non = TAUT_NON as usize;
    let p1_taut = p1.pINChI[taut_yes];
    let n1 = if !p1_taut.is_null() && heap.slice(p1_taut.as_const())?[0].nNumberOfAtoms != 0 {
        taut_yes
    } else {
        taut_non
    };
    let i1_pointer = p1.pINChI[n1];
    let p2_non = p2.pINChI[taut_non];
    let i2_pointer = if n1 == taut_yes
        && !p2_non.is_null()
        && heap.slice(p2_non.as_const())?[0].nNumberOfAtoms != 0
    {
        p2_non
    } else {
        Default::default()
    };
    let i1 = if i1_pointer.is_null() {
        None
    } else {
        Some(heap.slice(i1_pointer.as_const())?[0].clone())
    };
    let i2 = if i2_pointer.is_null() {
        None
    } else {
        Some(heap.slice(i2_pointer.as_const())?[0].clone())
    };
    let i1_live = i1.as_ref().filter(|inchi| inchi.bDeleted == 0);
    let i2_live = i2.as_ref().filter(|inchi| inchi.bDeleted == 0);

    let formula1_present = if let Some(inchi) = i1_live {
        !inchi.szHillFormula.is_null() && heap.slice(inchi.szHillFormula.as_const())?[0] != 0
    } else {
        false
    };
    let formula2_present = if let Some(inchi) = i2_live {
        !inchi.szHillFormula.is_null() && heap.slice(inchi.szHillFormula.as_const())?[0] != 0
    } else {
        false
    };
    let mut num_h1 = 0;
    let mut num_h2 = 0;
    if formula1_present {
        difference_segments[M][FORMULA] |= NOT_EQUAL;
        if formula2_present {
            let first = i1_live.ok_or(SourceHeapError::NullPointer)?;
            let second = i2_live.ok_or(SourceHeapError::NullPointer)?;
            if CompareHillFormulasNoH(
                heap,
                first.szHillFormula.as_const(),
                second.szHillFormula.as_const(),
                &mut num_h1,
                &mut num_h2,
            )? == 0
                && num_h1 == num_h2
            {
                difference_segments[F][FORMULA] |= EQUAL;
            } else {
                difference_segments[F][FORMULA] |= NOT_EQUAL;
            }
        } else {
            difference_segments[F][FORMULA] |= if i2.is_some() { IS_EMPTY } else { EQUAL };
        }
    } else {
        difference_segments[M][FORMULA] |= BOTH_EMPTY;
        difference_segments[F][FORMULA] |= if formula2_present {
            NOT_EQUAL
        } else {
            BOTH_EMPTY
        };
    }

    difference_segments[M][FORMULA] |= if i1_live.is_some_and(|inchi| inchi.lenConnTable > 1) {
        NOT_EQUAL
    } else {
        BOTH_EMPTY
    };

    let mobile_h_present = if let Some(inchi) = i1_live {
        let tautomeric = inchi.lenTautomer > 0
            && !inchi.nTautomer.is_null()
            && heap.slice(inchi.nTautomer.as_const())?[0] != 0;
        if tautomeric {
            true
        } else if !inchi.nNum_H.is_null() {
            let count = usize::try_from(inchi.nNumberOfAtoms.max(0))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            heap.slice(inchi.nNum_H.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .iter()
                .any(|value| *value != 0)
        } else {
            false
        }
    } else {
        false
    };
    difference_segments[M][H_ATOMS] |= if mobile_h_present {
        NOT_EQUAL
    } else {
        BOTH_EMPTY
    };
    let fixed_h_present = if let (Some(first), Some(second)) = (i1_live, i2_live) {
        if second.nNum_H_fixed.is_null() {
            false
        } else {
            let count = usize::try_from(first.nNumberOfAtoms.max(0))
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            heap.slice(second.nNum_H_fixed.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .iter()
                .any(|value| *value != 0)
        }
    } else {
        false
    };
    difference_segments[F][H_ATOMS] |= if fixed_h_present {
        NOT_EQUAL
    } else {
        BOTH_EMPTY
    };

    if let Some(first) = i1_live {
        difference_segments[M][CHARGE] |= if first.nTotalCharge != 0 {
            NOT_EQUAL
        } else {
            BOTH_EMPTY
        };
        if let Some(second) = i2_live {
            difference_segments[F][CHARGE] |= if first.nTotalCharge != 0 {
                if first.nTotalCharge == second.nTotalCharge {
                    EQUAL
                } else if second.nTotalCharge != 0 {
                    NOT_EQUAL
                } else {
                    IS_EMPTY
                }
            } else if second.nTotalCharge != 0 {
                NOT_EQUAL
            } else {
                BOTH_EMPTY
            };
        } else if i2.is_none() {
            let fixed_mark = if fix_transposed_charge_bug == 1
                && p1.ord_number != p2.ord_number
                && n1 == taut_yes
            {
                let p2_taut = p2.pINChI[taut_yes];
                if !p2_taut.is_null() {
                    let second_taut = &heap.slice(p2_taut.as_const())?[0];
                    if second_taut.bDeleted == 0 && second_taut.nNumberOfAtoms != 0 {
                        if first.nTotalCharge != 0 {
                            if first.nTotalCharge == second_taut.nTotalCharge {
                                EQUAL
                            } else if second_taut.nTotalCharge != 0 {
                                NOT_EQUAL
                            } else {
                                IS_EMPTY
                            }
                        } else if second_taut.nTotalCharge != 0 {
                            NOT_EQUAL
                        } else {
                            BOTH_EMPTY
                        }
                    } else if first.nTotalCharge != 0 {
                        EQUAL
                    } else {
                        BOTH_EMPTY
                    }
                } else if first.nTotalCharge != 0 {
                    EQUAL
                } else {
                    BOTH_EMPTY
                }
            } else if first.nTotalCharge != 0 {
                EQUAL
            } else {
                BOTH_EMPTY
            };
            difference_segments[F][CHARGE] |= fixed_mark;
        } else {
            difference_segments[F][CHARGE] |= if first.nTotalCharge != 0 {
                IS_EMPTY
            } else {
                BOTH_EMPTY
            };
        }
    } else {
        difference_segments[M][CHARGE] |= BOTH_EMPTY;
        if let Some(second) = i2_live {
            difference_segments[F][CHARGE] |= if second.nTotalCharge != 0 {
                NOT_EQUAL
            } else {
                BOTH_EMPTY
            };
        }
    }

    let stereo1 = if let Some(inchi) = i1_live.filter(|inchi| !inchi.Stereo.is_null()) {
        Some(heap.slice(inchi.Stereo.as_const())?[0].clone())
    } else {
        None
    };
    let iso_stereo1 = if let Some(inchi) = i1_live.filter(|inchi| !inchi.StereoIsotopic.is_null()) {
        Some(heap.slice(inchi.StereoIsotopic.as_const())?[0].clone())
    } else {
        None
    };
    let stereo2 = if let Some(inchi) = i2_live.filter(|inchi| !inchi.Stereo.is_null()) {
        Some(heap.slice(inchi.Stereo.as_const())?[0].clone())
    } else {
        None
    };
    let iso_stereo2 = if let Some(inchi) = i2_live.filter(|inchi| !inchi.StereoIsotopic.is_null()) {
        Some(heap.slice(inchi.StereoIsotopic.as_const())?[0].clone())
    } else {
        None
    };
    let has_bonds = |stereo: Option<&INChI_Stereo>| {
        stereo.is_some_and(|stereo| stereo.nNumberOfStereoBonds != 0)
    };
    let s1 = stereo1.as_ref();
    let is1 = iso_stereo1.as_ref();
    let s2 = stereo2.as_ref();
    let is2 = iso_stereo2.as_ref();

    difference_segments[M][SBONDS] |= if has_bonds(s1) { NOT_EQUAL } else { BOTH_EMPTY };
    difference_segments[F][SBONDS] |= if has_bonds(s2) {
        if has_bonds(s1) {
            if Eql_INChI_Stereo(heap, s1, EQL_SP2 as i32, s2, EQL_SP2 as i32, 0)? != 0 {
                EQUAL
            } else {
                NOT_EQUAL
            }
        } else {
            NOT_EQUAL
        }
    } else if has_bonds(s1) {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };
    difference_segments[MI][SBONDS] |= if has_bonds(is1) {
        if Eql_INChI_Stereo(heap, is1, EQL_SP2 as i32, s1, EQL_SP2 as i32, 0)? != 0 {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if has_bonds(s1) {
        EQUAL
    } else {
        BOTH_EMPTY
    };
    difference_segments[FI][SBONDS] |= if has_bonds(is2) {
        if Eql_INChI_Stereo(heap, is2, EQL_SP2 as i32, s2, EQL_SP2 as i32, 0)? != 0 {
            EQUAL
        } else if !has_bonds(s1)
            && !has_bonds(s2)
            && Eql_INChI_Stereo(heap, is2, EQL_SP2 as i32, is1, EQL_SP2 as i32, 0)? != 0
        {
            FI_EQ_MI
        } else {
            NOT_EQUAL
        }
    } else if has_bonds(s2) {
        EQUAL
    } else if has_bonds(is1) && !has_bonds(s1) {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };

    let rel_rac = [
        GetSp3RelRacAbs(i1_live, s1),
        GetSp3RelRacAbs(i1_live, is1),
        GetSp3RelRacAbs(i2_live, s2),
        GetSp3RelRacAbs(i2_live, is2),
    ];
    if rel_rac[M] != SP3_NONE {
        difference_segments[M][SATOMS] |= if rel_rac[M] & SP3_ANY != 0 {
            NOT_EQUAL
        } else {
            BOTH_EMPTY
        };
        difference_segments[M][SP3INV] |= if rel_rac[M] & SP3_ABS != 0 {
            NOT_EQUAL
        } else {
            BOTH_EMPTY
        };
        difference_segments[M][STYPE] |= if rel_rac[M] & SP3_TYPE != 0 {
            NOT_EQUAL
        } else {
            BOTH_EMPTY
        };
    } else {
        difference_segments[M][SATOMS] |= BOTH_EMPTY;
        difference_segments[M][SP3INV] |= BOTH_EMPTY;
        difference_segments[M][STYPE] |= BOTH_EMPTY;
    }

    difference_segments[F][SATOMS] |= if rel_rac[F] & SP3_ANY != 0 {
        if Eql_INChI_Stereo(heap, s2, EQL_SP3 as i32, s1, EQL_SP3 as i32, 0)? != 0 {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[M] & SP3_ANY != 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };
    difference_segments[MI][SATOMS] |= if rel_rac[MI] & SP3_ANY != 0 {
        if Eql_INChI_Stereo(heap, is1, EQL_SP3 as i32, s1, EQL_SP3 as i32, 0)? != 0 {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[M] & SP3_ANY != 0 {
        EQUAL
    } else {
        BOTH_EMPTY
    };
    difference_segments[FI][SATOMS] |= if rel_rac[FI] & SP3_ANY != 0 {
        if Eql_INChI_Stereo(heap, is2, EQL_SP3 as i32, s2, EQL_SP3 as i32, 0)? != 0 {
            EQUAL
        } else if rel_rac[M] & SP3_ANY == 0
            && rel_rac[F] & SP3_ANY == 0
            && Eql_INChI_Stereo(heap, is2, EQL_SP3 as i32, is1, EQL_SP3 as i32, 0)? != 0
        {
            FI_EQ_MI
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[F] & SP3_ANY != 0 {
        EQUAL
    } else if rel_rac[MI] & SP3_ANY != 0 && rel_rac[M] & SP3_ANY == 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };

    difference_segments[F][SP3INV] |= if rel_rac[F] & SP3_ABS != 0 {
        if rel_rac[M] & SP3_ABS == 0
            || s2.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
                != s1.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
        {
            NOT_EQUAL
        } else {
            EQUAL
        }
    } else if rel_rac[M] & SP3_ABS != 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };
    difference_segments[MI][SP3INV] |= if rel_rac[MI] & SP3_ABS != 0 {
        if rel_rac[M] & SP3_ABS != 0
            && is1.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
                == s1.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
        {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[M] & SP3_ABS != 0 {
        EQUAL
    } else {
        BOTH_EMPTY
    };
    difference_segments[FI][SP3INV] |= if rel_rac[FI] & SP3_ABS != 0 {
        if rel_rac[F] & SP3_ABS != 0
            && is2.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
                == s2.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
        {
            EQUAL
        } else if rel_rac[M] & SP3_ABS == 0
            && rel_rac[F] & SP3_ABS == 0
            && rel_rac[MI] & SP3_ABS != 0
            && is2.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
                == is1.ok_or(SourceHeapError::NullPointer)?.nCompInv2Abs
        {
            FI_EQ_MI
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[F] & SP3_ABS != 0 {
        EQUAL
    } else if rel_rac[MI] & SP3_ABS != 0 && rel_rac[M] & SP3_ABS == 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };

    difference_segments[F][STYPE] |= if rel_rac[F] & SP3_TYPE != 0 {
        if rel_rac[F] & SP3_TYPE == rel_rac[M] & SP3_TYPE {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[M] & SP3_TYPE != 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };
    difference_segments[MI][STYPE] |= if rel_rac[MI] & SP3_TYPE != 0 {
        if rel_rac[MI] & SP3_TYPE == rel_rac[M] & SP3_TYPE {
            EQUAL
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[M] & SP3_TYPE != 0 {
        EQUAL
    } else {
        BOTH_EMPTY
    };
    difference_segments[FI][STYPE] |= if rel_rac[FI] & SP3_TYPE != 0 {
        if rel_rac[FI] & SP3_TYPE == rel_rac[F] & SP3_TYPE {
            EQUAL
        } else if rel_rac[M] & SP3_TYPE == 0
            && rel_rac[F] & SP3_TYPE == 0
            && rel_rac[MI] & SP3_TYPE != 0
        {
            FI_EQ_MI
        } else {
            NOT_EQUAL
        }
    } else if rel_rac[F] & SP3_TYPE != 0 {
        EQUAL
    } else if rel_rac[MI] & SP3_TYPE != 0 && rel_rac[M] & SP3_TYPE == 0 {
        if i2.is_some() { IS_EMPTY } else { EQUAL }
    } else {
        BOTH_EMPTY
    };

    if p1.ord_number != p2.ord_number {
        difference_segments[F][TRANSP] |= NOT_EQUAL;
    }

    let isotopic1_present = i1_live.is_some_and(|inchi| {
        inchi.nNumberOfIsotopicAtoms != 0 || inchi.nNumberOfIsotopicTGroups != 0
    });
    difference_segments[MI][IATOMS] |= if isotopic1_present {
        NOT_EQUAL
    } else {
        BOTH_EMPTY
    };
    if let Some(second) = i2_live {
        let isotopic2_present =
            second.nNumberOfIsotopicAtoms != 0 || second.nNumberOfIsotopicTGroups != 0;
        if isotopic2_present {
            let Some(first) = i1_live else {
                difference_segments[FI][IATOMS] |= NOT_EQUAL;
                return Ok(0);
            };
            if second.nNumberOfIsotopicAtoms != first.nNumberOfIsotopicAtoms
                || second.nNumberOfIsotopicTGroups != first.nNumberOfIsotopicTGroups
            {
                difference_segments[FI][IATOMS] |= NOT_EQUAL;
            } else {
                let atom_count = usize::try_from(first.nNumberOfIsotopicAtoms.max(0))
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let mut different = false;
                if atom_count > 0 {
                    let atoms1 = heap
                        .slice(first.IsotopicAtom.as_const())?
                        .get(..atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let atoms2 = heap
                        .slice(second.IsotopicAtom.as_const())?
                        .get(..atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    different = atoms1 != atoms2;
                }
                if !different {
                    let group_count = usize::try_from(first.nNumberOfIsotopicTGroups.max(0))
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    if group_count > 0 {
                        let groups1 = heap
                            .slice(first.IsotopicTGroup.as_const())?
                            .get(..group_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let groups2 = heap
                            .slice(second.IsotopicTGroup.as_const())?
                            .get(..group_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        for index in 0..group_count {
                            if groups2[index].nTGroupNumber != groups1[index].nTGroupNumber
                                || groups2[index].nNum_T != groups1[index].nNum_T
                            {
                                different = true;
                                break;
                            }
                            let d_difference =
                                i32::from(groups2[index].nNum_D) - i32::from(groups1[index].nNum_D);
                            if d_difference != 0 {
                                return Ok(d_difference);
                            }
                            if groups2[index].nNum_H != groups1[index].nNum_H {
                                different = true;
                                break;
                            }
                        }
                    }
                }
                difference_segments[FI][IATOMS] |= if different { NOT_EQUAL } else { FI_EQ_MI };
            }
        } else if isotopic1_present {
            difference_segments[FI][IATOMS] |= IS_EMPTY;
        }
    } else if i2.is_none() {
        difference_segments[FI][IATOMS] |= if isotopic1_present { EQUAL } else { BOTH_EMPTY };
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompareInchiStereo(
    heap: &SourceHeap,
    stereo1: Option<&INChI_Stereo>,
    flags1: INCHI_MODE,
    stereo2: Option<&INChI_Stereo>,
    flags2: INCHI_MODE,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1607 CompareInchiStereo
    // INCHI✔️❌: int CompareInchiStereo(INChI_Stereo* Stereo1,
    // INCHI✔️❌:     INCHI_MODE nFlags1,
    // INCHI✔️❌:     INChI_Stereo* Stereo2,
    // INCHI✔️❌:     INCHI_MODE nFlags2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, num, ret;
    // INCHI✔️❌:     if (Stereo2 && Stereo1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  compare stereogenic bonds */
    // INCHI✔️❌:         num = inchi_min(Stereo2->nNumberOfStereoBonds, Stereo1->nNumberOfStereoBonds);
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)Stereo2->nBondAtom1[i] - (int)Stereo1->nBondAtom1[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)Stereo2->nBondAtom2[i] - (int)Stereo1->nBondAtom2[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)Stereo2->b_parity[i] - (int)Stereo1->b_parity[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if ((ret = (int)Stereo2->nNumberOfStereoBonds - (int)Stereo1->nNumberOfStereoBonds)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  compare stereogenic atoms */
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:         if (((nFlags1 | nFlags2) & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
    // INCHI✔️❌:             1 == Stereo2->nNumberOfStereoCenters &&
    // INCHI✔️❌:             1 == Stereo1->nNumberOfStereoCenters)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ; /*  do not compare single stereocenters in case of relative stereo */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌: #endif
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num = inchi_min(Stereo2->nNumberOfStereoCenters, Stereo1->nNumberOfStereoCenters);
    // INCHI✔️❌:             for (i = 0; i < num; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if ((ret = (int)Stereo2->nNumber[i] - (int)Stereo1->nNumber[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return ret;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if ((ret = (int)Stereo2->t_parity[i] - (int)Stereo1->t_parity[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return ret;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)Stereo2->nNumberOfStereoCenters - (int)Stereo1->nNumberOfStereoCenters)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             /*  compare stereo-abs-is-inverted flags  for non-relative, non-racemic */
    // INCHI✔️❌:             if (!((nFlags1 | nFlags2) & (INCHI_FLAG_RAC_STEREO | INCHI_FLAG_REL_STEREO)))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if ((ret = (Stereo2->nCompInv2Abs < 0) - (Stereo1->nCompInv2Abs < 0))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return ret;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (Stereo2 && (Stereo2->nNumberOfStereoBonds > 0 || Stereo2->nNumberOfStereoCenters > 0
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:             && /*  do not compare single stereocenters in case of relative stereo */
    // INCHI✔️❌:             !((nFlags2 & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
    // INCHI✔️❌:                 1 == Stereo2->nNumberOfStereoCenters
    // INCHI✔️❌:                 )
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (Stereo1 && (Stereo1->nNumberOfStereoBonds > 0 || Stereo1->nNumberOfStereoCenters > 0
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:                 && /*  do not compare single stereocenters in case of relative stereo */
    // INCHI✔️❌:                 !((nFlags1 & (INCHI_FLAG_REL_STEREO | INCHI_FLAG_RAC_STEREO)) &&
    // INCHI✔️❌:                     1 == Stereo1->nNumberOfStereoCenters
    // INCHI✔️❌:                     )
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareInchiStereo

    if let (Some(stereo1), Some(stereo2)) = (stereo1, stereo2) {
        let bond_count = stereo1
            .nNumberOfStereoBonds
            .min(stereo2.nNumberOfStereoBonds);
        if bond_count > 0 {
            let count =
                usize::try_from(bond_count).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let atom11 = heap
                .slice(stereo1.nBondAtom1.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom12 = heap
                .slice(stereo2.nBondAtom1.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom21 = heap
                .slice(stereo1.nBondAtom2.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let atom22 = heap
                .slice(stereo2.nBondAtom2.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let parity1 = heap
                .slice(stereo1.b_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let parity2 = heap
                .slice(stereo2.b_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in 0..count {
                let ret = i32::from(atom12[index]) - i32::from(atom11[index]);
                if ret != 0 {
                    return Ok(ret);
                }
                let ret = i32::from(atom22[index]) - i32::from(atom21[index]);
                if ret != 0 {
                    return Ok(ret);
                }
                let ret = i32::from(parity2[index]) - i32::from(parity1[index]);
                if ret != 0 {
                    return Ok(ret);
                }
            }
        }
        let ret = stereo2.nNumberOfStereoBonds - stereo1.nNumberOfStereoBonds;
        if ret != 0 {
            return Ok(ret);
        }
        let center_count = stereo1
            .nNumberOfStereoCenters
            .min(stereo2.nNumberOfStereoCenters);
        if center_count > 0 {
            let count = usize::try_from(center_count)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let number1 = heap
                .slice(stereo1.nNumber.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let number2 = heap
                .slice(stereo2.nNumber.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let parity1 = heap
                .slice(stereo1.t_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let parity2 = heap
                .slice(stereo2.t_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in 0..count {
                let ret = i32::from(number2[index]) - i32::from(number1[index]);
                if ret != 0 {
                    return Ok(ret);
                }
                let ret = i32::from(parity2[index]) - i32::from(parity1[index]);
                if ret != 0 {
                    return Ok(ret);
                }
            }
        }
        let ret = stereo2.nNumberOfStereoCenters - stereo1.nNumberOfStereoCenters;
        if ret != 0 {
            return Ok(ret);
        }
        if (flags1 | flags2) & INCHI_MODE::from(INCHI_FLAG_RAC_STEREO | INCHI_FLAG_REL_STEREO) == 0
        {
            let ret = i32::from(stereo2.nCompInv2Abs < 0) - i32::from(stereo1.nCompInv2Abs < 0);
            if ret != 0 {
                return Ok(ret);
            }
        }
    } else {
        if stereo2.is_some_and(|stereo| {
            stereo.nNumberOfStereoBonds > 0 || stereo.nNumberOfStereoCenters > 0
        }) {
            return Ok(1);
        }
        if stereo1.is_some_and(|stereo| {
            stereo.nNumberOfStereoBonds > 0 || stereo.nNumberOfStereoCenters > 0
        }) {
            return Ok(-1);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompINChI2(
    heap: &mut SourceHeap,
    p1: &INCHI_SORT,
    p2: &INCHI_SORT,
    b_taut: u32,
    compare_isotopic: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1712 CompINChI2
    // INCHI✔️❌: int CompINChI2(const INCHI_SORT* p1,
    // INCHI✔️❌:     const INCHI_SORT* p2,
    // INCHI✔️❌:     int bTaut,
    // INCHI✔️❌:     int bCompareIsotopic)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret, num, i, num_H1, num_H2;
    // INCHI✔️❌:
    // INCHI✔️❌:     const INChI* i1 = NULL; /* tautomeric if exists, otherwise non-tautomeric */
    // INCHI✔️❌:     const INChI* i2 = NULL; /* tautomeric if exists, otherwise non-tautomeric */
    // INCHI✔️❌:
    // INCHI✔️❌:     int   n1;               /* TAUT_YES if tautomeric i1 exists, otherwise TAUT_NON */
    // INCHI✔️❌:     int   n2;               /* TAUT_YES if tautomeric i2 exists, otherwise TAUT_NON */
    // INCHI✔️❌:
    // INCHI✔️❌:     const INChI* i1n = NULL; /* non-tautomeric if both tautomeric AND non-tautomeric exist */
    // INCHI✔️❌:     const INChI* i2n = NULL; /* non-tautomeric if both tautomeric AND non-tautomeric exist */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*const INChI *i1t = NULL;*/ /* temp for i1 if both tautomeric AND non-tautomeric exist */
    // INCHI✔️❌:     /*const INChI *i2t = NULL;*/ /* temp for i2 if both tautomeric AND non-tautomeric exist */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* INChI_Stereo *Stereo1, *Stereo2; */
    // INCHI✔️❌:
    // INCHI✔️❌:     n1 = (p1->pINChI[TAUT_YES] && p1->pINChI[TAUT_YES]->nNumberOfAtoms) ? TAUT_YES : TAUT_NON;
    // INCHI✔️❌:     n2 = (p2->pINChI[TAUT_YES] && p2->pINChI[TAUT_YES]->nNumberOfAtoms) ? TAUT_YES : TAUT_NON;
    // INCHI✔️❌:
    // INCHI✔️❌:     i1 = p1->pINChI[n1];
    // INCHI✔️❌:     i1n = (n1 == TAUT_YES && p1->pINChI[TAUT_NON] &&
    // INCHI✔️❌:         p1->pINChI[TAUT_NON]->nNumberOfAtoms) ? p1->pINChI[TAUT_NON] : (const INChI*)NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     i2 = p2->pINChI[n2];
    // INCHI✔️❌:     i2n = (n2 == TAUT_YES && p2->pINChI[TAUT_NON] &&
    // INCHI✔️❌:         p2->pINChI[TAUT_NON]->nNumberOfAtoms) ? p2->pINChI[TAUT_NON] : (const INChI*)NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* non-deleted-non-empty < deleted < empty */
    // INCHI✔️❌:     if (i1 && !i2)
    // INCHI✔️❌:         return -1;   /* non-empty is the smallest (first) */
    // INCHI✔️❌:     if (!i1 && i2)
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     if (!i1 && !i2)
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     if (i1->bDeleted && !i2->bDeleted)
    // INCHI✔️❌:         return 1;    /* deleted is the largest (last) among non-empty */
    // INCHI✔️❌:     if (!i1->bDeleted && i2->bDeleted)
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:
    // INCHI✔️❌:     num_H1 = num_H2 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* do not compare terminal H */
    // INCHI✔️❌:     if ((ret = CompareHillFormulasNoH(i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;  /* lexicographic order except the shorter one is greater (last): CH2O < CH2; C3XX < C2XX */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:             compare non-isotopic non-tautomeric part
    // INCHI✔️❌:      *********************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:      /* compare number of atoms (excluding terminal H) */
    // INCHI✔️❌:     if ((ret = i2->nNumberOfAtoms - i1->nNumberOfAtoms)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return ret; /*  more atoms first */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  compare elements  (excluding terminal H) */
    // INCHI✔️❌:     num = i1->nNumberOfAtoms;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* should always be equal if Hill formulas are same */
    // INCHI✔️❌:         if ((ret = (int)i2->nAtom[i] - (int)i1->nAtom[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret; /* greater periodic number first */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /**********************************************************
    // INCHI✔️❌:         compare connection tables
    // INCHI✔️❌:     ***********************************************************/
    // INCHI✔️❌:     if ((ret = i2->lenConnTable - i1->lenConnTable)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         return ret; /* longer connection table first */
    // INCHI✔️❌:     num = i2->lenConnTable;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = (int)i2->nConnTable[i] - (int)i1->nConnTable[i])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret; /* greater connection table first */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:       compare compare total number of H (inverse: H3 < H2 )
    // INCHI✔️❌:     **********************************************************/
    // INCHI✔️❌:     if ((ret = num_H2 - num_H1)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:       compare non-tautomeric num_H: N < NH3 < NH2 < NH
    // INCHI✔️❌:     **********************************************************/
    // INCHI✔️❌:     num = i1->nNumberOfAtoms;
    // INCHI✔️❌:     for (i = 0; i < num; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i2->nNum_H[i] != i1->nNum_H[i])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return !i2->nNum_H[i] ? 1 :  /* no H first */
    // INCHI✔️❌:                 !i1->nNum_H[i] ? -1 :
    // INCHI✔️❌:                 (int)i2->nNum_H[i] - (int)i1->nNum_H[i];
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:          compare non-isotopic tautomeric part
    // INCHI✔️❌:      *********************************************************/
    // INCHI✔️❌:     if ((ret = CompareTautNonIsoPartOfINChI(i1, i2))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*
    // INCHI✔️❌:     if ( ret = i2->lenTautomer - i1->lenTautomer )
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     num = inchi_min( i2->lenTautomer, i1->lenTautomer );
    // INCHI✔️❌:     for ( i = 0; i < num; i ++ ) {
    // INCHI✔️❌:         if ( ret = (int)i2->nTautomer[i] - (int)i1->nTautomer[i] )
    // INCHI✔️❌:            return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:      *                                                       *
    // INCHI✔️❌:      *  at this point both components are either tautomeric  *
    // INCHI✔️❌:      *  or non-tautomeric                                    *
    // INCHI✔️❌:      *                                                       *
    // INCHI✔️❌:      *********************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:      /*********************************************************
    // INCHI✔️❌:         non-tautomeric "fixed H" specific
    // INCHI✔️❌:       *********************************************************/
    // INCHI✔️❌:     if (TAUT_NON == bTaut && ((i1n && i1n->nNum_H_fixed) || (i2n && i2n->nNum_H_fixed))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* first, compare non-tautomeric chem. formulas -- they may be different */
    // INCHI✔️❌:         const char* f1 = (i1n /*&& i1n->nNum_H_fixed*/) ? i1n->szHillFormula : i1->szHillFormula;
    // INCHI✔️❌:         const char* f2 = (i2n /*&& i2n->nNum_H_fixed*/) ? i2n->szHillFormula : i2->szHillFormula;
    // INCHI✔️❌:         if (f1 && f2 && (ret = CompareHillFormulas(f1, f2)))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* secondly, compare fixed-H distribution */
    // INCHI✔️❌:         if (i1n && i1n->nNum_H_fixed && i2n && i2n->nNum_H_fixed)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num = inchi_min(i1n->nNumberOfAtoms, i2n->nNumberOfAtoms);
    // INCHI✔️❌:             for (i = 0; i < num; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (i2n->nNum_H_fixed[i] != i1n->nNum_H_fixed[i])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return !i2n->nNum_H_fixed[i] ? 1 : /* no fixed H first */
    // INCHI✔️❌:                         !i1n->nNum_H_fixed[i] ? -1 :
    // INCHI✔️❌:                         (int)i2n->nNum_H_fixed[i] - (int)i1n->nNum_H_fixed[i];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2n->nNumberOfAtoms - (int)i1n->nNumberOfAtoms)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /* should not happen <BRKPT> */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i1n && i1n->nNum_H_fixed)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num = i1n->nNumberOfAtoms;
    // INCHI✔️❌:                 for (i = 0; i < num; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* added 2004-05-04 */
    // INCHI✔️❌:                     if (i1n->nNum_H_fixed[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return -1; /* i1n->nNum_H_fixed[i] > 0? -1:1;*/
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* p1 is tautomeric, p2 is not tautomeric; this must have been detected earlier */
    // INCHI✔️❌:                 /*return -1;*/ /* has fixed H first *//* <BRKPT> */ /* removed 2004-05-04 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num = i2n->nNumberOfAtoms;
    // INCHI✔️❌:                 for (i = 0; i < num; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* added 2004-05-04 */
    // INCHI✔️❌:                     if (i2n->nNum_H_fixed[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         return 1; /* i2n->nNum_H_fixed[i] > 0? 1:-1;*/
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* p2 is tautomeric, p1 is not tautomeric; this must have been detected earlier */
    // INCHI✔️❌:                 /*return 1; */ /* has fixed H first *//* <BRKPT> */ /* removed 2004-05-04 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*************************************************************************
    // INCHI✔️❌:               if requested non-tautomeric comparison then
    // INCHI✔️❌:               prepare to compare non-taut non-isotopic stereo, etc.
    // INCHI✔️❌:      *************************************************************************/
    // INCHI✔️❌:     if (TAUT_NON == bTaut)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i1n)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*i1t = i1;*/
    // INCHI✔️❌:             i1 = i1n;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (i2n)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*i2t = i2;*/
    // INCHI✔️❌:             i2 = i2n;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*********************************************************
    // INCHI✔️❌:         compare non-isotopic stereo
    // INCHI✔️❌:      *********************************************************/
    // INCHI✔️❌:     ret = CompareInchiStereo(i1->Stereo, i1->nFlags, i2->Stereo, i2->nFlags);
    // INCHI✔️❌:     if (ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /*******************************************************
    // INCHI✔️❌:         do not switch back to tautomeric i1, i2
    // INCHI✔️❌:      *******************************************************/
    // INCHI✔️❌:      /* -- how to switch back --
    // INCHI✔️❌:      if ( i1t ) {
    // INCHI✔️❌:          i1  = i1t;
    // INCHI✔️❌:          i1t = NULL;
    // INCHI✔️❌:      }
    // INCHI✔️❌:      if ( i2t ) {
    // INCHI✔️❌:          i2  = i2t;
    // INCHI✔️❌:          i2t = NULL;
    // INCHI✔️❌:      }
    // INCHI✔️❌:      */
    // INCHI✔️❌:
    // INCHI✔️❌:      /******************************************************
    // INCHI✔️❌:           compare isotopic non-tautomeric part
    // INCHI✔️❌:       ******************************************************/
    // INCHI✔️❌:     if (bCompareIsotopic)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = i2->nNumberOfIsotopicAtoms - i1->nNumberOfIsotopicAtoms)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         num = i1->nNumberOfIsotopicAtoms;
    // INCHI✔️❌:         /*  compare isotopic atoms */
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nAtomNumber - (int)i1->IsotopicAtom[i].nAtomNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nIsoDifference - (int)i1->IsotopicAtom[i].nIsoDifference)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* compare isotopic H */
    // INCHI✔️❌:         /* if tautomeric comparison mode then here are compared only non-tautomeric H */
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_T - (int)i1->IsotopicAtom[i].nNum_T)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_D - (int)i1->IsotopicAtom[i].nNum_D)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicAtom[i].nNum_H - (int)i1->IsotopicAtom[i].nNum_H)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*****************************************************
    // INCHI✔️❌:              compare isotopic tautomeric part
    // INCHI✔️❌:          *****************************************************/
    // INCHI✔️❌:         if ((ret = i2->nNumberOfIsotopicTGroups - i1->nNumberOfIsotopicTGroups)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         num = i1->nNumberOfIsotopicTGroups;
    // INCHI✔️❌:         for (i = 0; i < num; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicTGroup[i].nTGroupNumber - (int)i1->IsotopicTGroup[i].nTGroupNumber)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicTGroup[i].nNum_T - (int)i1->IsotopicTGroup[i].nNum_T)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicTGroup[i].nNum_D - (int)i1->IsotopicTGroup[i].nNum_D)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((ret = (int)i2->IsotopicTGroup[i].nNum_H - (int)i1->IsotopicTGroup[i].nNum_H)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /****************************************************
    // INCHI✔️❌:             compare isotopic stereo
    // INCHI✔️❌:          ****************************************************/
    // INCHI✔️❌:         ret = CompareInchiStereo(i1->StereoIsotopic, i1->nFlags, i2->StereoIsotopic, i2->nFlags);
    // INCHI✔️❌:         if (ret)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /**********************************************************
    // INCHI✔️❌:         compare charges: non-charged first, then in order of
    // INCHI✔️❌:         ascending charges (negative first)
    // INCHI✔️❌:     ***********************************************************/
    // INCHI✔️❌:     if (i2->nTotalCharge && i1->nTotalCharge)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  both are charged; smaller charges first */
    // INCHI✔️❌:         ret = (int)i1->nTotalCharge - (int)i2->nTotalCharge;
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((ret = (i1->nTotalCharge ? 1 : 0) - (i2->nTotalCharge ? 1 : 0))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  only one is charged; uncharged first */
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* stable sort */
    // INCHI✔️❌:     /*ret = p1->ord_number - p2->ord_number;*/
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompINChI2

    let taut_yes = TAUT_YES as usize;
    let taut_non = TAUT_NON as usize;
    let p1_taut = p1.pINChI[taut_yes];
    let p2_taut = p2.pINChI[taut_yes];
    let n1 = if !p1_taut.is_null() && heap.slice(p1_taut.as_const())?[0].nNumberOfAtoms != 0 {
        taut_yes
    } else {
        taut_non
    };
    let n2 = if !p2_taut.is_null() && heap.slice(p2_taut.as_const())?[0].nNumberOfAtoms != 0 {
        taut_yes
    } else {
        taut_non
    };
    let i1_pointer = p1.pINChI[n1];
    let i2_pointer = p2.pINChI[n2];
    let mut i1 = if i1_pointer.is_null() {
        None
    } else {
        Some(heap.slice(i1_pointer.as_const())?[0].clone())
    };
    let mut i2 = if i2_pointer.is_null() {
        None
    } else {
        Some(heap.slice(i2_pointer.as_const())?[0].clone())
    };
    let p1_non = p1.pINChI[taut_non];
    let p2_non = p2.pINChI[taut_non];
    let i1n = if n1 == taut_yes
        && !p1_non.is_null()
        && heap.slice(p1_non.as_const())?[0].nNumberOfAtoms != 0
    {
        Some(heap.slice(p1_non.as_const())?[0].clone())
    } else {
        None
    };
    let i2n = if n2 == taut_yes
        && !p2_non.is_null()
        && heap.slice(p2_non.as_const())?[0].nNumberOfAtoms != 0
    {
        Some(heap.slice(p2_non.as_const())?[0].clone())
    } else {
        None
    };

    match (&i1, &i2) {
        (Some(_), None) => return Ok(-1),
        (None, Some(_)) => return Ok(1),
        (None, None) => return Ok(0),
        (Some(first), Some(second)) if first.bDeleted != 0 && second.bDeleted == 0 => return Ok(1),
        (Some(first), Some(second)) if first.bDeleted == 0 && second.bDeleted != 0 => {
            return Ok(-1);
        }
        _ => {}
    }
    let first = i1.as_ref().ok_or(SourceHeapError::NullPointer)?;
    let second = i2.as_ref().ok_or(SourceHeapError::NullPointer)?;
    let mut num_h1 = 0;
    let mut num_h2 = 0;
    let ret = CompareHillFormulasNoH(
        heap,
        first.szHillFormula.as_const(),
        second.szHillFormula.as_const(),
        &mut num_h1,
        &mut num_h2,
    )?;
    if ret != 0 {
        return Ok(ret);
    }
    let ret = second.nNumberOfAtoms.wrapping_sub(first.nNumberOfAtoms);
    if ret != 0 {
        return Ok(ret);
    }
    let atom_count = usize::try_from(first.nNumberOfAtoms.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if atom_count > 0 {
        let atoms1 = heap
            .slice(first.nAtom.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atoms2 = heap
            .slice(second.nAtom.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for index in 0..atom_count {
            let ret = i32::from(atoms2[index]) - i32::from(atoms1[index]);
            if ret != 0 {
                return Ok(ret);
            }
        }
    }
    let ret = second.lenConnTable.wrapping_sub(first.lenConnTable);
    if ret != 0 {
        return Ok(ret);
    }
    let connection_count = usize::try_from(second.lenConnTable.max(0))
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if connection_count > 0 {
        let connections1 = heap
            .slice(first.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let connections2 = heap
            .slice(second.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for index in 0..connection_count {
            let ret = i32::from(connections2[index]) - i32::from(connections1[index]);
            if ret != 0 {
                return Ok(ret);
            }
        }
    }
    let ret = num_h2.wrapping_sub(num_h1);
    if ret != 0 {
        return Ok(ret);
    }
    if atom_count > 0 {
        let hydrogens1 = heap
            .slice(first.nNum_H.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let hydrogens2 = heap
            .slice(second.nNum_H.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for index in 0..atom_count {
            if hydrogens2[index] != hydrogens1[index] {
                return Ok(if hydrogens2[index] == 0 {
                    1
                } else if hydrogens1[index] == 0 {
                    -1
                } else {
                    i32::from(hydrogens2[index]) - i32::from(hydrogens1[index])
                });
            }
        }
    }
    let ret = CompareTautNonIsoPartOfINChI(heap, first, second)?;
    if ret != 0 {
        return Ok(ret);
    }

    if b_taut == TAUT_NON
        && (i1n
            .as_ref()
            .is_some_and(|value| !value.nNum_H_fixed.is_null())
            || i2n
                .as_ref()
                .is_some_and(|value| !value.nNum_H_fixed.is_null()))
    {
        let formula1 = i1n.as_ref().unwrap_or(first).szHillFormula;
        let formula2 = i2n.as_ref().unwrap_or(second).szHillFormula;
        if !formula1.is_null() && !formula2.is_null() {
            let ret = CompareHillFormulas(heap, formula1.as_const(), formula2.as_const())?;
            if ret != 0 {
                return Ok(ret);
            }
        }
        let fixed1 = i1n.as_ref().filter(|value| !value.nNum_H_fixed.is_null());
        let fixed2 = i2n.as_ref().filter(|value| !value.nNum_H_fixed.is_null());
        match (fixed1, fixed2) {
            (Some(left), Some(right)) => {
                let count = usize::try_from(left.nNumberOfAtoms.min(right.nNumberOfAtoms).max(0))
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let values1 = heap
                    .slice(left.nNum_H_fixed.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let values2 = heap
                    .slice(right.nNum_H_fixed.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                for index in 0..count {
                    if values2[index] != values1[index] {
                        return Ok(if values2[index] == 0 {
                            1
                        } else if values1[index] == 0 {
                            -1
                        } else {
                            i32::from(values2[index]) - i32::from(values1[index])
                        });
                    }
                }
                let ret = right.nNumberOfAtoms.wrapping_sub(left.nNumberOfAtoms);
                if ret != 0 {
                    return Ok(ret);
                }
            }
            (Some(left), None) => {
                let count = usize::try_from(left.nNumberOfAtoms.max(0))
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                if heap
                    .slice(left.nNum_H_fixed.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .iter()
                    .any(|value| *value != 0)
                {
                    return Ok(-1);
                }
            }
            (None, Some(right)) => {
                let count = usize::try_from(right.nNumberOfAtoms.max(0))
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                if heap
                    .slice(right.nNum_H_fixed.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .iter()
                    .any(|value| *value != 0)
                {
                    return Ok(1);
                }
            }
            (None, None) => {}
        }
    }

    if b_taut == TAUT_NON {
        if let Some(fixed) = i1n {
            i1 = Some(fixed);
        }
        if let Some(fixed) = i2n {
            i2 = Some(fixed);
        }
    }
    let first = i1.as_ref().ok_or(SourceHeapError::NullPointer)?;
    let second = i2.as_ref().ok_or(SourceHeapError::NullPointer)?;
    let stereo1 = if first.Stereo.is_null() {
        None
    } else {
        Some(heap.slice(first.Stereo.as_const())?[0].clone())
    };
    let stereo2 = if second.Stereo.is_null() {
        None
    } else {
        Some(heap.slice(second.Stereo.as_const())?[0].clone())
    };
    let ret = CompareInchiStereo(
        heap,
        stereo1.as_ref(),
        first.nFlags,
        stereo2.as_ref(),
        second.nFlags,
    )?;
    if ret != 0 {
        return Ok(ret);
    }

    if compare_isotopic != 0 {
        let ret = second
            .nNumberOfIsotopicAtoms
            .wrapping_sub(first.nNumberOfIsotopicAtoms);
        if ret != 0 {
            return Ok(ret);
        }
        let isotope_count = usize::try_from(first.nNumberOfIsotopicAtoms.max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if isotope_count > 0 {
            let isotopes1 = heap
                .slice(first.IsotopicAtom.as_const())?
                .get(..isotope_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let isotopes2 = heap
                .slice(second.IsotopicAtom.as_const())?
                .get(..isotope_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in 0..isotope_count {
                for ret in [
                    i32::from(isotopes2[index].nAtomNumber)
                        - i32::from(isotopes1[index].nAtomNumber),
                    i32::from(isotopes2[index].nIsoDifference)
                        - i32::from(isotopes1[index].nIsoDifference),
                ] {
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
            for index in 0..isotope_count {
                for ret in [
                    i32::from(isotopes2[index].nNum_T) - i32::from(isotopes1[index].nNum_T),
                    i32::from(isotopes2[index].nNum_D) - i32::from(isotopes1[index].nNum_D),
                    i32::from(isotopes2[index].nNum_H) - i32::from(isotopes1[index].nNum_H),
                ] {
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
        }
        let ret = second
            .nNumberOfIsotopicTGroups
            .wrapping_sub(first.nNumberOfIsotopicTGroups);
        if ret != 0 {
            return Ok(ret);
        }
        let group_count = usize::try_from(first.nNumberOfIsotopicTGroups.max(0))
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if group_count > 0 {
            let groups1 = heap
                .slice(first.IsotopicTGroup.as_const())?
                .get(..group_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let groups2 = heap
                .slice(second.IsotopicTGroup.as_const())?
                .get(..group_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            for index in 0..group_count {
                for ret in [
                    i32::from(groups2[index].nTGroupNumber)
                        - i32::from(groups1[index].nTGroupNumber),
                    i32::from(groups2[index].nNum_T) - i32::from(groups1[index].nNum_T),
                    i32::from(groups2[index].nNum_D) - i32::from(groups1[index].nNum_D),
                    i32::from(groups2[index].nNum_H) - i32::from(groups1[index].nNum_H),
                ] {
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
        }
        let stereo1 = if first.StereoIsotopic.is_null() {
            None
        } else {
            Some(heap.slice(first.StereoIsotopic.as_const())?[0].clone())
        };
        let stereo2 = if second.StereoIsotopic.is_null() {
            None
        } else {
            Some(heap.slice(second.StereoIsotopic.as_const())?[0].clone())
        };
        let ret = CompareInchiStereo(
            heap,
            stereo1.as_ref(),
            first.nFlags,
            stereo2.as_ref(),
            second.nFlags,
        )?;
        if ret != 0 {
            return Ok(ret);
        }
    }

    if second.nTotalCharge != 0 && first.nTotalCharge != 0 {
        return Ok(first.nTotalCharge.wrapping_sub(second.nTotalCharge));
    }
    Ok(i32::from(first.nTotalCharge != 0) - i32::from(second.nTotalCharge != 0))
}

#[allow(non_snake_case)]
pub(crate) fn CompINChINonTaut2(
    heap: &mut SourceHeap,
    p1: &INCHI_SORT,
    p2: &INCHI_SORT,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2042 CompINChINonTaut2
    // INCHI✔️❌: int CompINChINonTaut2(const void* p1, const void* p2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:     ret = CompINChI2((const INCHI_SORT*)p1, (const INCHI_SORT*)p2, TAUT_NON, 1);
    // INCHI✔️❌: #if ( CANON_FIXH_TRANS == 1 )
    // INCHI✔️❌:     if (!ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* to obtain canonical transposition 2004-05-10 */
    // INCHI✔️❌:         ret = CompINChI2((const INCHI_SORT*)p1, (const INCHI_SORT*)p2, TAUT_YES, 1);
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     if (!ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* stable sort */
    // INCHI✔️❌:         ret = ((const INCHI_SORT*)p1)->ord_number - ((const INCHI_SORT*)p2)->ord_number;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompINChINonTaut2

    let mut ret = CompINChI2(heap, p1, p2, TAUT_NON, 1)?;
    if ret == 0 {
        ret = CompINChI2(heap, p1, p2, TAUT_YES, 1)?;
    }
    if ret == 0 {
        ret = i32::from(p1.ord_number) - i32::from(p2.ord_number);
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn CompINChITaut2(
    heap: &mut SourceHeap,
    p1: &INCHI_SORT,
    p2: &INCHI_SORT,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2064 CompINChITaut2
    // INCHI✔️❌: int CompINChITaut2(const void* p1, const void* p2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:     ret = CompINChI2((const INCHI_SORT*)p1, (const INCHI_SORT*)p2, TAUT_YES, 1);
    // INCHI✔️❌: #if ( CANON_FIXH_TRANS == 1 )
    // INCHI✔️❌:     if (!ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* to obtain canonical transposition 2004-05-10 */
    // INCHI✔️❌:         ret = CompINChI2((const INCHI_SORT*)p1, (const INCHI_SORT*)p2, TAUT_NON, 1);
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     if (!ret)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* stable sort */
    // INCHI✔️❌:         ret = ((const INCHI_SORT*)p1)->ord_number - ((const INCHI_SORT*)p2)->ord_number;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompINChITaut2

    let mut ret = CompINChI2(heap, p1, p2, TAUT_YES, 1)?;
    if ret == 0 {
        ret = CompINChI2(heap, p1, p2, TAUT_NON, 1)?;
    }
    if ret == 0 {
        ret = i32::from(p1.ord_number) - i32::from(p2.ord_number);
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn INChI_SegmentAction(c_dif_segs: i8) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1486 INChI_SegmentAction
    // INCHI✔️✔️: int INChI_SegmentAction(char cDifSegs)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     if (!(cDifSegs & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return INCHI_SEGM_OMIT;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if ((cDifSegs & DIFV_OUTPUT_EMPTY_T) && !(cDifSegs & DIFV_OUTPUT_EMPTY_F))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return INCHI_SEGM_EMPTY;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if ((cDifSegs & DIFV_OUTPUT_FILL_T))
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return INCHI_SEGM_FILL;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return INCHI_SEGM_OMIT; /* the control flow shoul never reach this point */
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: INChI_SegmentAction

    let value = i32::from(c_dif_segs);
    if value & crate::source_types::tagMarkDiff_DIFV_OUTPUT_OMIT_F as i32 == 0 {
        return crate::source_types::tagINChISegmAction_INCHI_SEGM_OMIT as i32;
    }
    if value & crate::source_types::tagMarkDiff_DIFV_OUTPUT_EMPTY_T as i32 != 0
        && value & crate::source_types::tagMarkDiff_DIFV_OUTPUT_EMPTY_F as i32 == 0
    {
        return crate::source_types::tagINChISegmAction_INCHI_SEGM_EMPTY as i32;
    }
    if value & crate::source_types::tagMarkDiff_DIFV_OUTPUT_FILL_T as i32 != 0 {
        return crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32;
    }
    crate::source_types::tagINChISegmAction_INCHI_SEGM_OMIT as i32
}

#[allow(non_snake_case)]
pub(crate) fn MarkUnusedAndEmptyLayers(s_dif_segs: &mut [[i8; 11]; 4]) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1525 MarkUnusedAndEmptyLayers
    // INCHI✔️❌: int MarkUnusedAndEmptyLayers(char sDifSegs[][DIFS_LENGTH])
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, nLayer, sBits, nFirstSegm;
    // INCHI✔️❌: #define nFirstFmlSegm   DIFS_f_FORMULA
    // INCHI✔️❌: #define nFirstIsoSegm   DIFS_i_IATOMS
    // INCHI✔️❌:     /* FI */
    // INCHI✔️❌:     nLayer = DIFL_FI;
    // INCHI✔️❌:     nFirstSegm = nFirstIsoSegm;
    // INCHI✔️❌:     sBits = 0;
    // INCHI✔️❌:     for (i = 0; i < DIFS_idf_LENGTH; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sBits |= sDifSegs[nLayer][i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!(sBits & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Omit the FI layer */
    // INCHI✔️❌:         memset(sDifSegs[nLayer], DIFV_BOTH_EMPTY, DIFS_idf_LENGTH); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (sDifSegs[nLayer][nFirstSegm] == DIFV_BOTH_EMPTY ||
    // INCHI✔️❌:             !(sDifSegs[nLayer][nFirstSegm] & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[nLayer][nFirstSegm] = DIFV_IS_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* MI */
    // INCHI✔️❌:     nLayer = DIFL_MI;
    // INCHI✔️❌:     nFirstSegm = nFirstIsoSegm;
    // INCHI✔️❌:     sBits = 0;
    // INCHI✔️❌:     for (i = 0; i < DIFS_idf_LENGTH; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sBits |= sDifSegs[nLayer][i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!(sBits & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Omit the MI layer */
    // INCHI✔️❌:         memset(sDifSegs[nLayer], DIFV_BOTH_EMPTY, DIFS_idf_LENGTH); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (sDifSegs[nLayer][nFirstSegm] == DIFV_BOTH_EMPTY ||
    // INCHI✔️❌:             !(sDifSegs[nLayer][nFirstSegm] & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[nLayer][nFirstSegm] = DIFV_IS_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* F */
    // INCHI✔️❌:     nLayer = DIFL_F;
    // INCHI✔️❌:     nFirstSegm = nFirstFmlSegm;
    // INCHI✔️❌:     sBits = 0;
    // INCHI✔️❌:     for (i = 0; i < DIFS_idf_LENGTH; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sBits |= sDifSegs[nLayer][i];
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!(sBits & DIFV_OUTPUT_OMIT_F) &&
    // INCHI✔️❌:         sDifSegs[DIFL_FI][nFirstIsoSegm] == DIFV_BOTH_EMPTY)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Omit the F layer: no non-iotopic and no isotopic segments */
    // INCHI✔️❌:         memset(sDifSegs[nLayer], DIFV_BOTH_EMPTY, DIFS_idf_LENGTH); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {    /* do not omit fixed-H layer */
    // INCHI✔️❌:         if (sDifSegs[nLayer][nFirstSegm] == DIFV_BOTH_EMPTY ||
    // INCHI✔️❌:             !(sDifSegs[nLayer][nFirstSegm] & DIFV_OUTPUT_OMIT_F))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             sDifSegs[nLayer][nFirstSegm] = DIFV_IS_EMPTY;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* M -- leave as it is */
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌:
    // INCHI✔️❌: #undef nFirstFmlSegm
    // INCHI✔️❌: #undef nFirstIsoSegm
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MarkUnusedAndEmptyLayers

    let first_iso = tagDiffINChISegments_DIFS_i_IATOMS as usize;
    let output_omit = tagMarkDiff_DIFV_OUTPUT_OMIT_F as i8;
    for layer in [
        tagDiffINChILayers_DIFL_FI as usize,
        tagDiffINChILayers_DIFL_MI as usize,
    ] {
        let bits = s_dif_segs[layer]
            .iter()
            .fold(0_i8, |bits, value| bits | value);
        if bits & output_omit == 0 {
            s_dif_segs[layer].fill(tagMarkDiff_DIFV_BOTH_EMPTY as i8);
        } else if s_dif_segs[layer][first_iso] == tagMarkDiff_DIFV_BOTH_EMPTY as i8
            || s_dif_segs[layer][first_iso] & output_omit == 0
        {
            s_dif_segs[layer][first_iso] = tagMarkDiff_DIFV_IS_EMPTY as i8;
        }
    }
    let fixed = tagDiffINChILayers_DIFL_F as usize;
    let formula = tagDiffINChISegments_DIFS_f_FORMULA as usize;
    let bits = s_dif_segs[fixed]
        .iter()
        .fold(0_i8, |bits, value| bits | value);
    if bits & output_omit == 0
        && s_dif_segs[tagDiffINChILayers_DIFL_FI as usize][first_iso]
            == tagMarkDiff_DIFV_BOTH_EMPTY as i8
    {
        s_dif_segs[fixed].fill(tagMarkDiff_DIFV_BOTH_EMPTY as i8);
    } else if s_dif_segs[fixed][formula] == tagMarkDiff_DIFV_BOTH_EMPTY as i8
        || s_dif_segs[fixed][formula] & output_omit == 0
    {
        s_dif_segs[fixed][formula] = tagMarkDiff_DIFV_IS_EMPTY as i8;
    }
    0
}

#[allow(non_snake_case)]
pub(crate) fn CompareReversedStereoINChI(
    heap: &SourceHeap,
    s1: Option<&INChI_Stereo>,
    s2: Option<&INChI_Stereo>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2555 CompareReversedStereoINChI
    // INCHI✔️❌: int CompareReversedStereoINChI(INChI_Stereo* s1/* InChI from reversed struct */,
    // INCHI✔️❌:     INChI_Stereo* s2 /* input InChI */
    // INCHI✔️❌: )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (s1 == NULL && s2 == NULL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((s1 == NULL) ^ (s2 == NULL))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         INChI_Stereo* s = s1 ? s1 : s2;
    // INCHI✔️❌:         if (s->nNumberOfStereoCenters || s->nNumberOfStereoBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 20; /* Diff: Missing Stereo */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (s1->nNumberOfStereoCenters != s2->nNumberOfStereoCenters)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 21;      /* Diff: Number of sp3 stereocenters */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (s1->nNumberOfStereoCenters > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (memcmp(s1->nNumber, s2->nNumber, s1->nNumberOfStereoCenters * sizeof(s1->nNumber[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 22;  /* Diff: sp3 stereocenter locations */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (memcmp(s1->t_parity, s2->t_parity, s1->nNumberOfStereoCenters * sizeof(s1->t_parity[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 23;  /* Diff: sp3 stereocenter parities */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (s1->nCompInv2Abs != s2->nCompInv2Abs && s1->nCompInv2Abs && s2->nCompInv2Abs)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 24;  /* Diff: sp3 inversion */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /*
    // INCHI✔️❌:         if ( s1->nNumberInv && s2->nNumberInv ) {
    // INCHI✔️❌:             if ( memcmp( s1->nNumberInv, s2->nNumberInv, s1->nNumberOfStereoCenters*sizeof(s1->nNumber[0]) ) )
    // INCHI✔️❌:                 return 25;
    // INCHI✔️❌:             if ( memcmp( s1->t_parityInv, s2->t_parityInv, s1->nNumberOfStereoCenters*sizeof(s1->t_parity[0]) ) )
    // INCHI✔️❌:                 return 26;
    // INCHI✔️❌:             if ( s1->nCompInv2Abs != s2->nCompInv2Abs ||
    // INCHI✔️❌:                  s1->bTrivialInv  != s2->bTrivialInv ) {
    // INCHI✔️❌:                 return 27;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         } else
    // INCHI✔️❌:         if ( s1->nNumberInv || s2->nNumberInv ) {
    // INCHI✔️❌:             return 28;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (s1->nNumberOfStereoBonds != s2->nNumberOfStereoBonds)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 25;      /* Diff: Number of stereobonds */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (s1->nNumberOfStereoBonds > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (memcmp(s1->nBondAtom1, s2->nBondAtom1, s1->nNumberOfStereoBonds * sizeof(s1->nBondAtom1[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 26; /* Diff: Stereobond 1st atom locations */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (memcmp(s1->nBondAtom2, s2->nBondAtom2, s1->nNumberOfStereoBonds * sizeof(s1->nBondAtom2[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 27; /* Diff: Stereobond 2nd atom locations */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (memcmp(s1->b_parity, s2->b_parity, s1->nNumberOfStereoBonds * sizeof(s1->b_parity[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 28; /* Diff: Stereobond parities */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareReversedStereoINChI

    let (Some(s1), Some(s2)) = (s1, s2) else {
        let Some(stereo) = s1.or(s2) else {
            return Ok(0);
        };
        return Ok(
            if stereo.nNumberOfStereoCenters != 0 || stereo.nNumberOfStereoBonds != 0 {
                20
            } else {
                0
            },
        );
    };
    if s1.nNumberOfStereoCenters != s2.nNumberOfStereoCenters {
        return Ok(21);
    }
    if s1.nNumberOfStereoCenters > 0 {
        let count = usize::try_from(s1.nNumberOfStereoCenters)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if heap
            .slice(s1.nNumber.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(s2.nNumber.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            return Ok(22);
        }
        if heap
            .slice(s1.t_parity.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(s2.t_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            return Ok(23);
        }
        if s1.nCompInv2Abs != s2.nCompInv2Abs && s1.nCompInv2Abs != 0 && s2.nCompInv2Abs != 0 {
            return Ok(24);
        }
    }
    if s1.nNumberOfStereoBonds != s2.nNumberOfStereoBonds {
        return Ok(25);
    }
    if s1.nNumberOfStereoBonds > 0 {
        let count = usize::try_from(s1.nNumberOfStereoBonds)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if heap
            .slice(s1.nBondAtom1.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(s2.nBondAtom1.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            return Ok(26);
        }
        if heap
            .slice(s1.nBondAtom2.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(s2.nBondAtom2.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            return Ok(27);
        }
        if heap
            .slice(s1.b_parity.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(s2.b_parity.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            return Ok(28);
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn CompareReversedStereoINChI2(
    heap: &SourceHeap,
    s1: Option<&INChI_Stereo>,
    s2: Option<&INChI_Stereo>,
    picr: &mut ICR,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2635 CompareReversedStereoINChI2
    /*
    int CompareReversedStereoINChI2(INChI_Stereo* s1/* InChI from reversed struct */,
        INChI_Stereo* s2 /* input InChI */,
        ICR* picr)
    {
        int ret = 0;
        int j1, j2, num_dif, num_extra_undf, num_miss_undf, num_in1_only, num_in2_only; /* djb-rwth: removing redundant variables */
        int bAddSb = !(picr->num_sb_undef_in1_only + picr->num_sb_in1_only + picr->num_sb_in2_only);
        int bAddSc = !(picr->num_sc_undef_in1_only + picr->num_sc_in1_only + picr->num_sc_in2_only);

        int nNumSc1 = s1 ? s1->nNumberOfStereoCenters : 0;
        int nNumSc2 = s2 ? s2->nNumberOfStereoCenters : 0;
        int nNumSb1 = s1 ? s1->nNumberOfStereoBonds : 0;
        int nNumSb2 = s2 ? s2->nNumberOfStereoBonds : 0;

        /* djb-rwth: redundant part of the condition corrected; fixing a NULL pointer dereference */
        if (s1 && s2 && (nNumSc1 || nNumSc2) && (nNumSc1 != nNumSc2 ||
            memcmp(s1->nNumber, s2->nNumber, nNumSc1 * sizeof(s1->nNumber[0])) ||
            memcmp(s1->t_parity, s2->t_parity, nNumSc1 * sizeof(s1->t_parity[0]))))
        {
            num_dif = num_extra_undf = num_miss_undf = num_in1_only = num_in2_only = 0;
            for (j1 = j2 = 0; j1 < nNumSc1 && j2 < nNumSc2; )
            {
                if (s1->nNumber[j1] == s2->nNumber[j2])
                {
                    if (s1->t_parity[j1] != s2->t_parity[j2]) /* djb-rwth: removing redundant code */
                    {
                        num_dif++;
                    }
                    j1++;
                    j2++;
                }
                else
                {
                    if (s1->nNumber[j1] < s2->nNumber[j2])
                    {
                        num_in1_only++;
                        if (s1->t_parity[j1] == AB_PARITY_UNDF)
                        {
                            num_extra_undf++;
                        }
                        if (bAddSc)
                        {
                            if (picr->num_sc_in1_only < ICR_MAX_SC_IN1_ONLY)
                            {
                                picr->sc_in1_only[picr->num_sc_in1_only++] = j1;
                            }
                            if (s1->t_parity[j1] == AB_PARITY_UNDF)
                            {
                                if (picr->num_sc_undef_in1_only < ICR_MAX_SC_UNDF)
                                {
                                    picr->sc_undef_in1_only[picr->num_sc_undef_in1_only++] = j1;
                                }
                            }
                        }
                        j1++;
                    }
                    else
                    {
                        num_in2_only++;
                        if (s2->t_parity[j2] == AB_PARITY_UNDF)
                        {
                            num_miss_undf++;
                        }
                        if (bAddSc)
                        {
                            if (picr->num_sc_in2_only < ICR_MAX_SC_IN2_ONLY)
                            {
                                picr->sc_in2_only[picr->num_sc_in2_only++] = j2;
                            }
                            if (s2->t_parity[j2] == AB_PARITY_UNDF)
                            {
                                if (picr->num_sc_undef_in2_only < ICR_MAX_SC_UNDF)
                                {
                                    picr->sc_undef_in2_only[picr->num_sc_undef_in2_only++] = j1;
                                }
                            }
                        }
                        j2++;
                    }
                }
            }
            while (j1 < nNumSc1)
            {
                if (s1->t_parity[j1] == AB_PARITY_UNDF)
                {
                    num_extra_undf++;
                }
                num_in1_only++;
                if (bAddSc)
                {
                    if (picr->num_sc_in1_only < ICR_MAX_SC_IN1_ONLY)
                    {
                        picr->sc_in1_only[picr->num_sc_in1_only++] = j1;
                    }
                    if (s1->t_parity[j1] == AB_PARITY_UNDF)
                    {
                        if (picr->num_sc_undef_in1_only < ICR_MAX_SC_UNDF)
                        {
                            picr->sc_undef_in1_only[picr->num_sc_undef_in1_only++] = j1;
                        }
                    }
                }
                j1++;
            }
            while (j2 < nNumSc2)
            {
                if (s2->t_parity[j2] == AB_PARITY_UNDF)
                {
                    num_miss_undf++;
                }
                num_in2_only++;
                if (bAddSc)
                {
                    if (picr->num_sc_in2_only < ICR_MAX_SC_IN2_ONLY)
                    {
                        picr->sc_in2_only[picr->num_sc_in2_only++] = j2;
                    }
                }
                j2++;
            }
            if (num_dif)
            {
                ret |= IDIF_SC_PARITY;
            }
            if (num_in1_only)
            {
                if (num_extra_undf)
                {
                    ret |= IDIF_SC_EXTRA_UNDF;
                }
                if (num_in1_only != num_extra_undf)
                {
                    ret |= IDIF_SC_EXTRA;
                }
            }
            if (num_in2_only)
            {
                if (num_miss_undf)
                {
                    ret |= IDIF_SC_MISS_UNDF;
                }
                if (num_in2_only != num_miss_undf)
                {
                    ret |= IDIF_SC_MISS;
                }
            }
        }
        if (s1 && s2 && s1->nCompInv2Abs != s2->nCompInv2Abs && s1->nCompInv2Abs && s2->nCompInv2Abs)
        {
            ret |= IDIF_SC_INV;
        }

        if ((nNumSb1 || nNumSb2) &&
            (nNumSb1 != nNumSb2 ||
                memcmp(s1->nBondAtom1, s2->nBondAtom1, nNumSb1 * sizeof(s1->nBondAtom1[0])) ||
                memcmp(s1->nBondAtom2, s2->nBondAtom2, nNumSb1 * sizeof(s1->nBondAtom2[0])) ||
                memcmp(s1->b_parity, s2->b_parity, nNumSb1 * sizeof(s1->b_parity[0]))))
        {
            num_dif = num_extra_undf = num_miss_undf = num_in1_only = num_in2_only = 0; /* djb-rwth: removing redundant code */
            for (j1 = j2 = 0; j1 < nNumSb1 && j2 < nNumSb2; )
            {
                if (s1->nBondAtom1[j1] == s2->nBondAtom1[j2] &&
                    s1->nBondAtom2[j1] == s2->nBondAtom2[j2])
                {
                    if (s1->b_parity[j1] != s2->b_parity[j2]) /* djb-rwth: removing redundant code */
                    {
                        num_dif++;
                    }
                    j1++;
                    j2++;
                }
                else
                {
                    if (s1->nBondAtom1[j1] < s2->nBondAtom1[j2] ||
                        (s1->nBondAtom1[j1] == s2->nBondAtom1[j2] && s1->nBondAtom2[j1] < s2->nBondAtom2[j2])) /* djb-rwth: addressing LLVM warning */
                    {
                        num_in1_only++;
                        if (s1->b_parity[j1] == AB_PARITY_UNDF)
                        {
                            num_extra_undf++;
                        }
                        if (bAddSb)
                        {
                            if (picr->num_sb_in1_only < ICR_MAX_SB_IN1_ONLY)
                            {
                                picr->sb_in1_only[picr->num_sb_in1_only++] = j1;
                            }
                            if (s1->b_parity[j1] == AB_PARITY_UNDF)
                            {
                                if (picr->num_sb_undef_in1_only < ICR_MAX_SB_UNDF)
                                {
                                    picr->sb_undef_in1_only[picr->num_sb_undef_in1_only++] = j1;
                                }
                            }
                        }
                        j1++;
                    }
                    else
                    {
                        num_in2_only++;
                        if (s2->b_parity[j2] == AB_PARITY_UNDF)
                        {
                            num_miss_undf++;
                        }
                        if (bAddSb)
                        {
                            if (picr->num_sb_in2_only < ICR_MAX_SB_IN2_ONLY)
                            {
                                picr->sb_in2_only[picr->num_sb_in2_only++] = j2;
                            }
                            if (s2->b_parity[j2] == AB_PARITY_UNDF)
                            {
                                if (picr->num_sb_undef_in2_only < ICR_MAX_SB_UNDF)
                                {
                                    picr->sb_undef_in2_only[picr->num_sb_undef_in2_only++] = j1;
                                }
                            }
                        }
                        j2++;
                    }
                }
            }
            while (j1 < nNumSb1)
            {
                num_in1_only++;
                if (s1->b_parity[j1] == AB_PARITY_UNDF)
                {
                    num_extra_undf++;
                }
                if (bAddSb)
                {
                    if (picr->num_sb_in1_only < ICR_MAX_SB_IN1_ONLY)
                    {
                        picr->sb_in1_only[picr->num_sb_in1_only++] = j1;
                    }
                    if (s1->b_parity[j1] == AB_PARITY_UNDF)
                    {
                        if (picr->num_sb_undef_in1_only < ICR_MAX_SB_UNDF)
                        {
                            picr->sb_undef_in1_only[picr->num_sb_undef_in1_only++] = j1;
                        }
                    }
                }
                j1++;
            }
            while (j2 < nNumSb2)
            {
                num_in2_only++;
                if (s2->b_parity[j2] == AB_PARITY_UNDF)
                {
                    num_miss_undf++;
                }
                if (bAddSb)
                {
                    if (picr->num_sb_in2_only < ICR_MAX_SB_IN2_ONLY)
                    {
                        picr->sb_in2_only[picr->num_sb_in2_only++] = j2;
                    }
                    if (s2->b_parity[j2] == AB_PARITY_UNDF)
                    {
                        if (picr->num_sb_undef_in2_only < ICR_MAX_SB_UNDF)
                        {
                            picr->sb_undef_in2_only[picr->num_sb_undef_in2_only++] = j1;
                        }
                    }
                }
                j2++;
            }
            if (num_dif)
            {
                ret |= IDIF_SB_PARITY;
            }
            if (num_in1_only)
            {
                if (num_extra_undf)
                {
                    ret |= IDIF_SB_EXTRA_UNDF;
                }
                if (num_in1_only != num_extra_undf)
                {
                    ret |= IDIF_SB_EXTRA;
                }
            }
            if (num_in2_only)
            {
                if (num_miss_undf)
                {
                    ret |= IDIF_SB_MISS_UNDF;
                }
                if (num_in2_only != num_miss_undf)
                {
                    ret |= IDIF_SB_MISS;
                }
            }
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: CompareReversedStereoINChI2
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CompareReversedStereoINChI2
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: ICR_MAX_SC_*=32, ICR_MAX_SB_*=32, and AB_PARITY_UNDF=3.
    // INCHI✔️❌: Checked SourceHeap slice lookup adds overhead versus memcmp/direct pointers.
    // END INCHI ACTIVE MACRO CONFIGURATION: CompareReversedStereoINChI2

    let add_sb = picr
        .num_sb_undef_in1_only
        .wrapping_add(picr.num_sb_in1_only)
        .wrapping_add(picr.num_sb_in2_only)
        == 0;
    let add_sc = picr
        .num_sc_undef_in1_only
        .wrapping_add(picr.num_sc_in1_only)
        .wrapping_add(picr.num_sc_in2_only)
        == 0;
    let num_sc1 = s1.map_or(0, |stereo| stereo.nNumberOfStereoCenters);
    let num_sc2 = s2.map_or(0, |stereo| stereo.nNumberOfStereoCenters);
    let num_sb1 = s1.map_or(0, |stereo| stereo.nNumberOfStereoBonds);
    let num_sb2 = s2.map_or(0, |stereo| stereo.nNumberOfStereoBonds);
    let mut ret = 0_i32;

    if let (Some(s1), Some(s2)) = (s1, s2) {
        if num_sc1 != 0 || num_sc2 != 0 {
            let count1 =
                usize::try_from(num_sc1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let count2 =
                usize::try_from(num_sc2).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            let numbers1 = if count1 == 0 {
                &[][..]
            } else {
                heap.slice(s1.nNumber.as_const())?
                    .get(..count1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            };
            let numbers2 = if count2 == 0 {
                &[][..]
            } else {
                heap.slice(s2.nNumber.as_const())?
                    .get(..count2)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            };
            let parities1 = if count1 == 0 {
                &[][..]
            } else {
                heap.slice(s1.t_parity.as_const())?
                    .get(..count1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            };
            let parities2 = if count2 == 0 {
                &[][..]
            } else {
                heap.slice(s2.t_parity.as_const())?
                    .get(..count2)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            };
            if num_sc1 != num_sc2 || numbers1 != numbers2 || parities1 != parities2 {
                let mut j1 = 0_usize;
                let mut j2 = 0_usize;
                let mut num_dif = 0_i32;
                let mut num_extra_undf = 0_i32;
                let mut num_miss_undf = 0_i32;
                let mut num_in1_only = 0_i32;
                let mut num_in2_only = 0_i32;
                while j1 < count1 && j2 < count2 {
                    if numbers1[j1] == numbers2[j2] {
                        num_dif = num_dif.wrapping_add(i32::from(parities1[j1] != parities2[j2]));
                        j1 += 1;
                        j2 += 1;
                    } else if numbers1[j1] < numbers2[j2] {
                        num_in1_only = num_in1_only.wrapping_add(1);
                        let undefined = parities1[j1] == AB_PARITY_UNDF as i8;
                        num_extra_undf = num_extra_undf.wrapping_add(i32::from(undefined));
                        if add_sc {
                            if picr.num_sc_in1_only < ICR_MAX_SC_IN1_ONLY as i32 {
                                picr.sc_in1_only[picr.num_sc_in1_only as usize] = j1 as AT_NUMB;
                                picr.num_sc_in1_only += 1;
                            }
                            if undefined && picr.num_sc_undef_in1_only < ICR_MAX_SC_UNDF as i32 {
                                picr.sc_undef_in1_only[picr.num_sc_undef_in1_only as usize] =
                                    j1 as AT_NUMB;
                                picr.num_sc_undef_in1_only += 1;
                            }
                        }
                        j1 += 1;
                    } else {
                        num_in2_only = num_in2_only.wrapping_add(1);
                        let undefined = parities2[j2] == AB_PARITY_UNDF as i8;
                        num_miss_undf = num_miss_undf.wrapping_add(i32::from(undefined));
                        if add_sc {
                            if picr.num_sc_in2_only < ICR_MAX_SC_IN2_ONLY as i32 {
                                picr.sc_in2_only[picr.num_sc_in2_only as usize] = j2 as AT_NUMB;
                                picr.num_sc_in2_only += 1;
                            }
                            if undefined && picr.num_sc_undef_in2_only < ICR_MAX_SC_UNDF as i32 {
                                picr.sc_undef_in2_only[picr.num_sc_undef_in2_only as usize] =
                                    j1 as AT_NUMB;
                                picr.num_sc_undef_in2_only += 1;
                            }
                        }
                        j2 += 1;
                    }
                }
                while j1 < count1 {
                    let undefined = parities1[j1] == AB_PARITY_UNDF as i8;
                    num_extra_undf = num_extra_undf.wrapping_add(i32::from(undefined));
                    num_in1_only = num_in1_only.wrapping_add(1);
                    if add_sc {
                        if picr.num_sc_in1_only < ICR_MAX_SC_IN1_ONLY as i32 {
                            picr.sc_in1_only[picr.num_sc_in1_only as usize] = j1 as AT_NUMB;
                            picr.num_sc_in1_only += 1;
                        }
                        if undefined && picr.num_sc_undef_in1_only < ICR_MAX_SC_UNDF as i32 {
                            picr.sc_undef_in1_only[picr.num_sc_undef_in1_only as usize] =
                                j1 as AT_NUMB;
                            picr.num_sc_undef_in1_only += 1;
                        }
                    }
                    j1 += 1;
                }
                while j2 < count2 {
                    let undefined = parities2[j2] == AB_PARITY_UNDF as i8;
                    num_miss_undf = num_miss_undf.wrapping_add(i32::from(undefined));
                    num_in2_only = num_in2_only.wrapping_add(1);
                    if add_sc && picr.num_sc_in2_only < ICR_MAX_SC_IN2_ONLY as i32 {
                        picr.sc_in2_only[picr.num_sc_in2_only as usize] = j2 as AT_NUMB;
                        picr.num_sc_in2_only += 1;
                    }
                    j2 += 1;
                }
                if num_dif != 0 {
                    ret |= tagInchiDiffBits_IDIF_SC_PARITY as i32;
                }
                if num_in1_only != 0 {
                    if num_extra_undf != 0 {
                        ret |= tagInchiDiffBits_IDIF_SC_EXTRA_UNDF as i32;
                    }
                    if num_in1_only != num_extra_undf {
                        ret |= tagInchiDiffBits_IDIF_SC_EXTRA as i32;
                    }
                }
                if num_in2_only != 0 {
                    if num_miss_undf != 0 {
                        ret |= tagInchiDiffBits_IDIF_SC_MISS_UNDF as i32;
                    }
                    if num_in2_only != num_miss_undf {
                        ret |= tagInchiDiffBits_IDIF_SC_MISS as i32;
                    }
                }
            }
        }
        if s1.nCompInv2Abs != s2.nCompInv2Abs && s1.nCompInv2Abs != 0 && s2.nCompInv2Abs != 0 {
            ret |= tagInchiDiffBits_IDIF_SC_INV as i32;
        }
    }

    if num_sb1 != 0 || num_sb2 != 0 {
        let s1 = s1.ok_or(SourceHeapError::PointerOutOfBounds)?;
        let s2 = s2.ok_or(SourceHeapError::PointerOutOfBounds)?;
        let count1 =
            usize::try_from(num_sb1).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let count2 =
            usize::try_from(num_sb2).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atom11 = if count1 == 0 {
            &[][..]
        } else {
            heap.slice(s1.nBondAtom1.as_const())?
                .get(..count1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        let atom12 = if count1 == 0 {
            &[][..]
        } else {
            heap.slice(s1.nBondAtom2.as_const())?
                .get(..count1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        let parity1 = if count1 == 0 {
            &[][..]
        } else {
            heap.slice(s1.b_parity.as_const())?
                .get(..count1)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        let atom21 = if count2 == 0 {
            &[][..]
        } else {
            heap.slice(s2.nBondAtom1.as_const())?
                .get(..count2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        let atom22 = if count2 == 0 {
            &[][..]
        } else {
            heap.slice(s2.nBondAtom2.as_const())?
                .get(..count2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        let parity2 = if count2 == 0 {
            &[][..]
        } else {
            heap.slice(s2.b_parity.as_const())?
                .get(..count2)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        if num_sb1 != num_sb2 || atom11 != atom21 || atom12 != atom22 || parity1 != parity2 {
            let mut j1 = 0_usize;
            let mut j2 = 0_usize;
            let mut num_dif = 0_i32;
            let mut num_extra_undf = 0_i32;
            let mut num_miss_undf = 0_i32;
            let mut num_in1_only = 0_i32;
            let mut num_in2_only = 0_i32;
            while j1 < count1 && j2 < count2 {
                if atom11[j1] == atom21[j2] && atom12[j1] == atom22[j2] {
                    num_dif = num_dif.wrapping_add(i32::from(parity1[j1] != parity2[j2]));
                    j1 += 1;
                    j2 += 1;
                } else if atom11[j1] < atom21[j2]
                    || (atom11[j1] == atom21[j2] && atom12[j1] < atom22[j2])
                {
                    num_in1_only += 1;
                    let undefined = parity1[j1] == AB_PARITY_UNDF as i8;
                    num_extra_undf += i32::from(undefined);
                    if add_sb {
                        if picr.num_sb_in1_only < ICR_MAX_SB_IN1_ONLY as i32 {
                            picr.sb_in1_only[picr.num_sb_in1_only as usize] = j1 as AT_NUMB;
                            picr.num_sb_in1_only += 1;
                        }
                        if undefined && picr.num_sb_undef_in1_only < ICR_MAX_SB_UNDF as i32 {
                            picr.sb_undef_in1_only[picr.num_sb_undef_in1_only as usize] =
                                j1 as AT_NUMB;
                            picr.num_sb_undef_in1_only += 1;
                        }
                    }
                    j1 += 1;
                } else {
                    num_in2_only += 1;
                    let undefined = parity2[j2] == AB_PARITY_UNDF as i8;
                    num_miss_undf += i32::from(undefined);
                    if add_sb {
                        if picr.num_sb_in2_only < ICR_MAX_SB_IN2_ONLY as i32 {
                            picr.sb_in2_only[picr.num_sb_in2_only as usize] = j2 as AT_NUMB;
                            picr.num_sb_in2_only += 1;
                        }
                        if undefined && picr.num_sb_undef_in2_only < ICR_MAX_SB_UNDF as i32 {
                            picr.sb_undef_in2_only[picr.num_sb_undef_in2_only as usize] =
                                j1 as AT_NUMB;
                            picr.num_sb_undef_in2_only += 1;
                        }
                    }
                    j2 += 1;
                }
            }
            while j1 < count1 {
                num_in1_only += 1;
                let undefined = parity1[j1] == AB_PARITY_UNDF as i8;
                num_extra_undf += i32::from(undefined);
                if add_sb {
                    if picr.num_sb_in1_only < ICR_MAX_SB_IN1_ONLY as i32 {
                        picr.sb_in1_only[picr.num_sb_in1_only as usize] = j1 as AT_NUMB;
                        picr.num_sb_in1_only += 1;
                    }
                    if undefined && picr.num_sb_undef_in1_only < ICR_MAX_SB_UNDF as i32 {
                        picr.sb_undef_in1_only[picr.num_sb_undef_in1_only as usize] = j1 as AT_NUMB;
                        picr.num_sb_undef_in1_only += 1;
                    }
                }
                j1 += 1;
            }
            while j2 < count2 {
                num_in2_only += 1;
                let undefined = parity2[j2] == AB_PARITY_UNDF as i8;
                num_miss_undf += i32::from(undefined);
                if add_sb {
                    if picr.num_sb_in2_only < ICR_MAX_SB_IN2_ONLY as i32 {
                        picr.sb_in2_only[picr.num_sb_in2_only as usize] = j2 as AT_NUMB;
                        picr.num_sb_in2_only += 1;
                    }
                    if undefined && picr.num_sb_undef_in2_only < ICR_MAX_SB_UNDF as i32 {
                        picr.sb_undef_in2_only[picr.num_sb_undef_in2_only as usize] = j1 as AT_NUMB;
                        picr.num_sb_undef_in2_only += 1;
                    }
                }
                j2 += 1;
            }
            if num_dif != 0 {
                ret |= tagInchiDiffBits_IDIF_SB_PARITY as i32;
            }
            if num_in1_only != 0 {
                if num_extra_undf != 0 {
                    ret |= tagInchiDiffBits_IDIF_SB_EXTRA_UNDF as i32;
                }
                if num_in1_only != num_extra_undf {
                    ret |= tagInchiDiffBits_IDIF_SB_EXTRA as i32;
                }
            }
            if num_in2_only != 0 {
                if num_miss_undf != 0 {
                    ret |= tagInchiDiffBits_IDIF_SB_MISS_UNDF as i32;
                }
                if num_in2_only != num_miss_undf {
                    ret |= tagInchiDiffBits_IDIF_SB_MISS as i32;
                }
            }
        }
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn CompareIcr(
    picr1: &ICR,
    picr2: &ICR,
    pin1: Option<&mut INCHI_MODE>,
    pin2: Option<&mut INCHI_MODE>,
    mask: INCHI_MODE,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:3189 CompareIcr
    // INCHI✔️✔️: complete active source frame follows verbatim; active IDIFF
    // INCHI✔️✔️: masks end at bit 30 and both implementations use one bitwise pass.
    /*
    int CompareIcr(ICR* picr1,
        ICR* picr2,
        INCHI_MODE* pin1,
        INCHI_MODE* pin2,
        INCHI_MODE mask)
    {
        int nNumExtraBits1 = 0, nNumExtraBits2 = 0, bit1, bit2;
        INCHI_MODE Flg1 = picr1->flags, Flg2 = picr2->flags, cur_bit = 1, in1, in2;
        int i, ret;

        /* compare flags */
        in1 = in2 = 0;
        for (i = 0; Flg1 || Flg2; i++, Flg1 >>= 1, Flg2 >>= 1, cur_bit <<= 1)
        {
            if (!(mask & cur_bit))
            {
                continue;
            }
            bit1 = Flg1 & 1;
            bit2 = Flg2 & 1;
            if (bit1 && !bit2)
            {
                in1 |= 1 << i;
                nNumExtraBits1++;
            }
            else
            {
                if (!bit1 && bit2)
                {
                    in2 |= 1 << i;
                    nNumExtraBits2++;
                }
            }
        }
        if (nNumExtraBits1 && !nNumExtraBits2)
        {
            ret = 1;
        }
        else
        {
            if (!nNumExtraBits1 && nNumExtraBits2)
            {
                ret = -1;
            }
            else
            {
                if (!in1 && !in2)
                {
                    ret = 0;
                }
                else
                {
                    ret = 2; /* compare produced undefined results */
                }
            }
        }

        if (pin1)
        {
            *pin1 = in1;
        }
        if (pin2)
        {
            *pin2 = in2;
        }

        /* more detailed compare not implemented */

        return ret;
    }
    */
    // END INCHI C FUNCTION: CompareIcr
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CompareIcr
    // INCHI✔️✔️: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; INCHI_MODE=unsigned long=u64.
    // END INCHI ACTIVE MACRO CONFIGURATION: CompareIcr
    let mut n_num_extra_bits1 = 0_i32;
    let mut n_num_extra_bits2 = 0_i32;
    let mut flg1 = picr1.flags;
    let mut flg2 = picr2.flags;
    let mut cur_bit = 1_u64;
    let mut in1 = 0_u64;
    let mut in2 = 0_u64;
    let mut i = 0_u32;

    while flg1 != 0 || flg2 != 0 {
        if mask & cur_bit != 0 {
            let bit1 = flg1 & 1;
            let bit2 = flg2 & 1;
            if bit1 != 0 && bit2 == 0 {
                in1 |= 1_u64 << i;
                n_num_extra_bits1 += 1;
            } else if bit1 == 0 && bit2 != 0 {
                in2 |= 1_u64 << i;
                n_num_extra_bits2 += 1;
            }
        }
        i += 1;
        flg1 >>= 1;
        flg2 >>= 1;
        cur_bit = cur_bit.wrapping_shl(1);
    }

    let ret = if n_num_extra_bits1 != 0 && n_num_extra_bits2 == 0 {
        1
    } else if n_num_extra_bits1 == 0 && n_num_extra_bits2 != 0 {
        -1
    } else if in1 == 0 && in2 == 0 {
        0
    } else {
        2
    };

    if let Some(pin1) = pin1 {
        *pin1 = in1;
    }
    if let Some(pin2) = pin2 {
        *pin2 = in2;
    }

    ret
}

#[allow(non_snake_case)]
pub(crate) fn CompareReversedINChI2(
    heap: &mut SourceHeap,
    i1: Option<&INChI>,
    i2: Option<&INChI>,
    a1: Option<&INChI_Aux>,
    a2: Option<&INChI_Aux>,
    picr: &mut ICR,
    err: &mut i32,
) -> Result<INCHI_MODE, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:3262 CompareReversedINChI2
    // INCHI✔️❌: complete source frame follows verbatim; endpoint staging adds Rust-owned vectors.
    /*
    INCHI_MODE CompareReversedINChI2(INChI* i1 /* InChI from reversed struct */,
        INChI* i2 /* input InChI */,
        INChI_Aux* a1,
        INChI_Aux* a2,
        ICR* picr,
        int* err)
    {
        INCHI_MODE ret = 0;
        INChI_Stereo* Stereo1 = NULL, * Stereo2 = NULL;
        int  n1, n2, m, j, j1, j2, ret2, num_H1, num_H2;

        *err = 0;

        memset(picr, 0, sizeof(*picr)); /* djb-rwth: memset_s C11/Annex K variant? */

        if (i1 == NULL && i2 == NULL)
        {
            return 0;
        }
        if ((i1 == NULL) ^ (i2 == NULL))
        {
            ret |= IDIF_PROBLEM; /* one InChI exists while another doesn't */
            goto exit_function;
        }

        if (i1->nErrorCode == i2->nErrorCode)
        {
            if (i1->nErrorCode)
            {
                ret |= IDIF_PROBLEM; /* both InChI have same error codes */
                goto exit_function;
            }
        }
        else
        {
            ret |= IDIF_PROBLEM; /* at least one InChI has an error code */
            goto exit_function;
        }

        if (i1->nNumberOfAtoms != i2->nNumberOfAtoms)
        {
            ret |= IDIF_NUM_AT;
            goto exit_function;
        }
        if (i1->nNumberOfAtoms > 0)
        {
            if (memcmp(i1->nAtom, i2->nAtom, i1->nNumberOfAtoms * sizeof(i1->nAtom[0])))
            {
                ret |= IDIF_ATOMS;
                goto exit_function;
            }
            /* IDIF_NON_TAUT_H,  IDIF_MORE_FH, IDIF_LESS_FH */
            if (memcmp(i1->nNum_H, i2->nNum_H, i1->nNumberOfAtoms * sizeof(i1->nNum_H[0])))
            {
                ret |= IDIF_POSITION_H;
                for (j1 = 0; j1 < i1->nNumberOfAtoms; j1++)
                {
                    if (i1->nNum_H[j1] != i2->nNum_H[j1] && picr->num_diff_pos_H < ICR_MAX_DIFF_FIXED_H)
                    {
                        picr->diff_pos_H_at[picr->num_diff_pos_H] = j1;
                        picr->diff_pos_H_nH[picr->num_diff_pos_H] = i1->nNum_H[j1] - i2->nNum_H[j1];
                        picr->num_diff_pos_H++;
                    }
                }
            }
            /* fixed H */
            if (i1->nNum_H_fixed || i2->nNum_H_fixed)
            {
                int bHasFixedH1 = 0, bHasFixedH2 = 0, i;
                if (i1->nNum_H_fixed)
                {
                    for (i = 0; i < i1->nNumberOfAtoms; i++)
                    {
                        if (i1->nNum_H_fixed[i])
                        {
                            bHasFixedH1++;
                        }
                    }
                }
                if (i2->nNum_H_fixed)
                {
                    for (i = 0; i < i2->nNumberOfAtoms; i++)
                    {
                        if (i2->nNum_H_fixed[i])
                        {
                            bHasFixedH2++;
                        }
                    }
                }
                if (bHasFixedH1 && !bHasFixedH2)
                {
                    for (i = j = 0; i < i1->nNumberOfAtoms; i++)
                    {
                        if (i1->nNum_H_fixed[i])
                        {
                            if (j < ICR_MAX_DIFF_FIXED_H)
                            {
                                picr->fixed_H_at1_more[j] = i;
                                picr->fixed_H_nH1_more[j] = i1->nNum_H_fixed[i];
                                j++;
                            }
                        }
                    }
                    picr->num_fixed_H1_more = j;
                    ret |= IDIF_MORE_FH; /* Extra Fixed-H */
                }
                else
                {
                    if (!bHasFixedH1 && bHasFixedH2)
                    {
                        for (i = j = 0; i < i2->nNumberOfAtoms; i++)
                        {
                            if (i2->nNum_H_fixed[i])
                            {
                                if (j < ICR_MAX_DIFF_FIXED_H)
                                {
                                    picr->fixed_H_at2_more[j] = i;
                                    picr->fixed_H_nH2_more[j] = i2->nNum_H_fixed[i];
                                    j++;
                                }
                            }
                        }
                        picr->num_fixed_H2_more = j;
                        ret |= IDIF_LESS_FH; /* Missed Fixed-H */
                    }
                    else
                    {
                        if (bHasFixedH1 && bHasFixedH2 &&
                            memcmp(i1->nNum_H_fixed, i2->nNum_H_fixed, i1->nNumberOfAtoms * sizeof(i1->nNum_H_fixed[0])))
                        {
                            for (i = j1 = j2 = 0; i < i1->nNumberOfAtoms; i++)
                            {
                                if (i1->nNum_H_fixed[i] > i2->nNum_H_fixed[i])
                                {
                                    if (j1 < ICR_MAX_DIFF_FIXED_H)
                                    {
                                        picr->fixed_H_at1_more[j1] = i;
                                        picr->fixed_H_nH1_more[j1] = i1->nNum_H_fixed[i] - i2->nNum_H_fixed[i];
                                        j1++;
                                    }
                                }
                                else
                                {
                                    if (i1->nNum_H_fixed[i] < i2->nNum_H_fixed[i])
                                    {
                                        if (j2 < ICR_MAX_DIFF_FIXED_H)
                                        {
                                            picr->fixed_H_at2_more[j2] = i;
                                            picr->fixed_H_nH2_more[j2] = i2->nNum_H_fixed[i] - i1->nNum_H_fixed[i];
                                            j2++;
                                        }
                                    }
                                }
                            }
                            ret |= (j1 ? IDIF_MORE_FH : 0) | (j2 ? IDIF_LESS_FH : 0);
                            picr->num_fixed_H1_more = j1;
                            picr->num_fixed_H2_more = j2;
                        }
                    }
                }
            }
        }

        /* compare formulas and H */
        num_H1 = 0;
        num_H2 = 0;
        ret2 = CompareHillFormulasNoH(i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2);
        picr->tot_num_H1 = num_H1;
        picr->tot_num_H2 = num_H2;
        if (ret2)
        {
            ret |= IDIF_NUM_EL;
            goto exit_function;
        }
        if (num_H1 > num_H2)
        {
            ret |= IDIF_MORE_H;
        }
        if (num_H1 < num_H2)
        {
            ret |= IDIF_LESS_H;
        }
        if (i1->lenConnTable != i2->lenConnTable)
        {
            ret |= IDIF_CON_LEN;
            goto exit_function;
        }
        else
        {
            if (i1->lenConnTable > 0 && memcmp(i1->nConnTable, i2->nConnTable, i1->lenConnTable * sizeof(i1->nConnTable[0])))
            {
                ret |= IDIF_CON_TBL;
                goto exit_function;
            }
        }

        /* output special cases: different number of t-groups, different sizes of t-groups, different endpoints */
        /* in isotopic or deprotonated cases i1->lenTautomer == 1 && i1->nTautomer[0] = 0 */

        /*
            if ( i1->lenTautomer != i2->lenTautomer && (i1->lenTautomer > 1 || i2->lenTautomer > 1) ) {
                ret |=  IDIF_TAUT_LEN;
            }
        */

        /* compare number of t-groups */
        n1 = i1->lenTautomer ? i1->nTautomer[0] : 0;
        n2 = i2->lenTautomer ? i2->nTautomer[0] : 0;
        if (!n1 && n2)
        {
            ret |= IDIF_NO_TAUT;
        }
        else
        {
            if (n1 && !n2)
            {
                ret |= IDIF_WRONG_TAUT;
            }
            else
            {
                if (n1 == 1 && n2 > 1)
                {
                    ret |= IDIF_SINGLE_TG;
                }
                else
                {
                    if (n1 > 1 && n2 == 1)
                    {
                        ret |= IDIF_MULTIPLE_TG;
                    }
                    else
                    {
                        if (n1 != n2)
                        {
                            ret |= IDIF_NUM_TG;
                        }
                    }
                }
            }
        }

        if (n1 || n2)
        {
            /* number of endpoints */
            int num1 = 0, num2 = 0, num_M1 = 0, num_M2 = 0;
            int len, num_eq, num_in1_only, num_in2_only;
            AT_NUMB* pe1 = (AT_NUMB*)inchi_malloc(((long long)i1->lenTautomer + 1) * sizeof(pe1[0])); /* djb-rwth: cast operator added */
            AT_NUMB* pe2 = (AT_NUMB*)inchi_malloc(((long long)i2->lenTautomer + 1) * sizeof(pe2[0])); /* djb-rwth: cast operator added */
            num_H1 = num_H2 = 0;
            /* collect endpoints, H, (-) */
            if (!pe1 || !pe2)
            {
                if (pe1)
                {
                    inchi_free(pe1);
                }
                if (pe2)
                {
                    inchi_free(pe2);
                }
                *err = -1; /* allocation error */
                goto exit_function;
            }
            for (m = 1; m < i1->lenTautomer; m += len)
            {
                len = i1->nTautomer[m++];
                num_H1 += i1->nTautomer[m];
                num_M1 += i1->nTautomer[m + 1];
                for (j = 2; j < len; j++)
                {
                    pe1[num1++] = i1->nTautomer[m + j];
                }
            }
            for (m = 1; m < i2->lenTautomer; m += len)
            {
                len = i2->nTautomer[m++];
                num_H2 += i2->nTautomer[m];
                num_M2 += i2->nTautomer[m + 1];
                for (j = 2; j < len; j++)
                {
                    pe2[num2++] = i2->nTautomer[m + j];
                }
            }
            picr->num_taut_H1 = num_H1;
            picr->num_taut_H2 = num_H2;
            picr->num_taut_M1 = num_M1;
            picr->num_taut_M2 = num_M2;
            /* sort endpoints */
            insertions_sort_AT_RANK(pe1, num1);
            insertions_sort_AT_RANK(pe2, num2);
            /* compare */
            /*
            if ( num1 < num2 ) {
                ret |= IDIF_LESS_TG_ENDP;
            } else
            if ( num1 > num2 ) {
                ret |= IDIF_MORE_TG_ENDP;
            }
            */
            /* compare all */
            num_eq = num_in1_only = num_in2_only = 0;
            for (j1 = j2 = 0; j1 < num1 && j2 < num2; )
            {
                if (pe1[j1] == pe2[j2])
                {
                    j1++;
                    j2++;
                    num_eq++;
                }
                else
                {
                    if (pe1[j1] < pe2[j2])
                    { /* BC: fixed, was pe2[j1] 2006-03-27 */
                        if (picr->num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY)
                        {
                            picr->endp_in1_only[picr->num_endp_in1_only++] = pe1[j1];
                        }
                        j1++;
                        num_in1_only++;
                    }
                    else
                    {
                        if (picr->num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY)
                        {
                            picr->endp_in2_only[picr->num_endp_in2_only++] = pe2[j2];
                        }
                        j2++;
                        num_in2_only++;
                    }
                }
            }
            while (j1 < num1)
            {
                if (picr->num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY)
                {
                    picr->endp_in1_only[picr->num_endp_in1_only++] = pe1[j1];
                }
                j1++;
                num_in1_only++;
            }
            while (j2 < num2)
            {
                if (picr->num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY)
                {
                    picr->endp_in2_only[picr->num_endp_in2_only++] = pe2[j2];
                }
                j2++;
                num_in2_only++;
            }
            if (num_in1_only)
            {
                ret |= IDIF_EXTRA_TG_ENDP;
            }
            if (num_in2_only)
            {
                ret |= IDIF_MISS_TG_ENDP;
            }
            if (!num_in1_only && !num_in2_only && num_eq)
            {
                ; /* same t-groups endpoints */
            }
            else
            {
                ret |= IDIF_DIFF_TG_ENDP;
            }
            inchi_free(pe1);
            inchi_free(pe2);
        }

        if ((i1->lenTautomer > 1 && i2->lenTautomer > 1) &&
            (i1->lenTautomer != i2->lenTautomer ||
                memcmp(i1->nTautomer, i2->nTautomer, i1->lenTautomer * sizeof(i1->nTautomer[0]))))
        {
            ret |= IDIF_TG;
        }

        if (i1->nNumberOfIsotopicAtoms != i2->nNumberOfIsotopicAtoms)
        {
            ret |= IDIF_NUM_ISO_AT;
        }
        else
        {
            if (i1->nNumberOfIsotopicAtoms > 0 && memcmp(i1->IsotopicAtom, i2->IsotopicAtom, i1->nNumberOfIsotopicAtoms * sizeof(i1->IsotopicAtom[0])))
            {
                ret |= IDIF_ISO_AT;
            }
        }
        if (i1->nTotalCharge != i2->nTotalCharge)
        {
            ret |= IDIF_CHARGE;
        }
        if (a1 && a1->nNumRemovedProtons && (!a2 || a2->nNumRemovedProtons != a1->nNumRemovedProtons))
        {
            ret |= IDIF_REM_PROT;
        }
        if (a1 && (!a2 ||
            a2->nNumRemovedIsotopicH[0] != a1->nNumRemovedIsotopicH[0] ||
            a2->nNumRemovedIsotopicH[1] != a1->nNumRemovedIsotopicH[1] ||
            a2->nNumRemovedIsotopicH[2] != a1->nNumRemovedIsotopicH[2]))
        {
            ret |= IDIF_REM_ISO_H;
        }

        /*
        if ( i1->nPossibleLocationsOfIsotopicH && i2->nPossibleLocationsOfIsotopicH ) {
            if ( i1->nPossibleLocationsOfIsotopicH[0] != i2->nPossibleLocationsOfIsotopicH[0] ||
                 memcmp(i1->nPossibleLocationsOfIsotopicH, i2->nPossibleLocationsOfIsotopicH,
                        sizeof(i1->nPossibleLocationsOfIsotopicH[0])*i1->nPossibleLocationsOfIsotopicH[0]) )
                return 18;
        } else
        if ( !i1->nPossibleLocationsOfIsotopicH != !i2->nPossibleLocationsOfIsotopicH ) {
            return 19;
        }
        */

        if (i1->StereoIsotopic &&
            i1->StereoIsotopic->nNumberOfStereoBonds + i1->StereoIsotopic->nNumberOfStereoCenters)
        {
            Stereo1 = i1->StereoIsotopic;
        }
        else
        {
            Stereo1 = i1->Stereo;
        }
        if (i2->StereoIsotopic &&
            i2->StereoIsotopic->nNumberOfStereoBonds + i2->StereoIsotopic->nNumberOfStereoCenters)
        {
            Stereo2 = i2->StereoIsotopic;
        }
        else
        {
            Stereo2 = i2->Stereo;
        }
        ret |= CompareReversedStereoINChI2(Stereo1, Stereo2, picr);

    exit_function:

        picr->flags = ret;

        return ret;
    }
    */
    // END INCHI C FUNCTION: CompareReversedINChI2
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: CompareReversedINChI2
    // INCHI✔️❌: READ_INCHI_STRING=1; COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux.
    // INCHI✔️❌: inchi_malloc/free resolve to the active mode.h malloc/free macros.
    // INCHI✔️❌: ICR_MAX_DIFF_FIXED_H=32 and ICR_MAX_ENDP_IN{1,2}_ONLY=32.
    // END INCHI ACTIVE MACRO CONFIGURATION: CompareReversedINChI2

    *err = 0;
    *picr = ICR::default();
    let mut ret = 0_u32;

    macro_rules! finish {
        () => {{
            picr.flags = INCHI_MODE::from(ret);
            return Ok(INCHI_MODE::from(ret));
        }};
    }

    let (i1, i2) = match (i1, i2) {
        (None, None) => return Ok(0),
        (Some(i1), Some(i2)) => (i1, i2),
        _ => {
            ret |= tagInchiDiffBits_IDIF_PROBLEM;
            finish!();
        }
    };

    if i1.nErrorCode == i2.nErrorCode {
        if i1.nErrorCode != 0 {
            ret |= tagInchiDiffBits_IDIF_PROBLEM;
            finish!();
        }
    } else {
        ret |= tagInchiDiffBits_IDIF_PROBLEM;
        finish!();
    }

    if i1.nNumberOfAtoms != i2.nNumberOfAtoms {
        ret |= tagInchiDiffBits_IDIF_NUM_AT;
        finish!();
    }
    if i1.nNumberOfAtoms > 0 {
        let atom_count = usize::try_from(i1.nNumberOfAtoms)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atoms1 = heap
            .slice(i1.nAtom.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atoms2 = heap
            .slice(i2.nAtom.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atoms1 != atoms2 {
            ret |= tagInchiDiffBits_IDIF_ATOMS;
            finish!();
        }

        let hydrogens1 = heap
            .slice(i1.nNum_H.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let hydrogens2 = heap
            .slice(i2.nNum_H.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if hydrogens1 != hydrogens2 {
            ret |= tagInchiDiffBits_IDIF_POSITION_H;
            for index in 0..atom_count {
                if hydrogens1[index] != hydrogens2[index]
                    && picr.num_diff_pos_H < ICR_MAX_DIFF_FIXED_H as i32
                {
                    let output = picr.num_diff_pos_H as usize;
                    picr.diff_pos_H_at[output] = index as AT_NUMB;
                    picr.diff_pos_H_nH[output] = hydrogens1[index].wrapping_sub(hydrogens2[index]);
                    picr.num_diff_pos_H = picr.num_diff_pos_H.wrapping_add(1);
                }
            }
        }

        if !i1.nNum_H_fixed.is_null() || !i2.nNum_H_fixed.is_null() {
            let fixed1 = if i1.nNum_H_fixed.is_null() {
                None
            } else {
                Some(
                    heap.slice(i1.nNum_H_fixed.as_const())?
                        .get(..atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let fixed2 = if i2.nNum_H_fixed.is_null() {
                None
            } else {
                Some(
                    heap.slice(i2.nNum_H_fixed.as_const())?
                        .get(..atom_count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let has_fixed1 = fixed1.map_or(0_i32, |values| {
                values.iter().fold(0_i32, |count, value| {
                    count.wrapping_add(i32::from(*value != 0))
                })
            });
            let has_fixed2 = fixed2.map_or(0_i32, |values| {
                values.iter().fold(0_i32, |count, value| {
                    count.wrapping_add(i32::from(*value != 0))
                })
            });

            if has_fixed1 != 0 && has_fixed2 == 0 {
                let values = fixed1.expect("nonzero count requires a source array");
                let mut count = 0_i32;
                for (index, &value) in values.iter().enumerate() {
                    if value != 0 && count < ICR_MAX_DIFF_FIXED_H as i32 {
                        picr.fixed_H_at1_more[count as usize] = index as AT_NUMB;
                        picr.fixed_H_nH1_more[count as usize] = value;
                        count = count.wrapping_add(1);
                    }
                }
                picr.num_fixed_H1_more = count;
                ret |= tagInchiDiffBits_IDIF_MORE_FH;
            } else if has_fixed1 == 0 && has_fixed2 != 0 {
                let values = fixed2.expect("nonzero count requires a source array");
                let mut count = 0_i32;
                for (index, &value) in values.iter().enumerate() {
                    if value != 0 && count < ICR_MAX_DIFF_FIXED_H as i32 {
                        picr.fixed_H_at2_more[count as usize] = index as AT_NUMB;
                        picr.fixed_H_nH2_more[count as usize] = value;
                        count = count.wrapping_add(1);
                    }
                }
                picr.num_fixed_H2_more = count;
                ret |= tagInchiDiffBits_IDIF_LESS_FH;
            } else if has_fixed1 != 0 && has_fixed2 != 0 {
                let values1 = fixed1.expect("nonzero count requires a source array");
                let values2 = fixed2.expect("nonzero count requires a source array");
                if values1 != values2 {
                    let mut count1 = 0_i32;
                    let mut count2 = 0_i32;
                    for index in 0..atom_count {
                        if values1[index] > values2[index] {
                            if count1 < ICR_MAX_DIFF_FIXED_H as i32 {
                                picr.fixed_H_at1_more[count1 as usize] = index as AT_NUMB;
                                picr.fixed_H_nH1_more[count1 as usize] =
                                    values1[index].wrapping_sub(values2[index]);
                                count1 = count1.wrapping_add(1);
                            }
                        } else if values1[index] < values2[index]
                            && count2 < ICR_MAX_DIFF_FIXED_H as i32
                        {
                            picr.fixed_H_at2_more[count2 as usize] = index as AT_NUMB;
                            picr.fixed_H_nH2_more[count2 as usize] =
                                values2[index].wrapping_sub(values1[index]);
                            count2 = count2.wrapping_add(1);
                        }
                    }
                    if count1 != 0 {
                        ret |= tagInchiDiffBits_IDIF_MORE_FH;
                    }
                    if count2 != 0 {
                        ret |= tagInchiDiffBits_IDIF_LESS_FH;
                    }
                    picr.num_fixed_H1_more = count1;
                    picr.num_fixed_H2_more = count2;
                }
            }
        }
    }

    let mut total_h1 = 0_i32;
    let mut total_h2 = 0_i32;
    let formula_difference = CompareHillFormulasNoH(
        heap,
        i1.szHillFormula.as_const(),
        i2.szHillFormula.as_const(),
        &mut total_h1,
        &mut total_h2,
    )?;
    picr.tot_num_H1 = total_h1;
    picr.tot_num_H2 = total_h2;
    if formula_difference != 0 {
        ret |= tagInchiDiffBits_IDIF_NUM_EL;
        finish!();
    }
    if total_h1 > total_h2 {
        ret |= tagInchiDiffBits_IDIF_MORE_H;
    }
    if total_h1 < total_h2 {
        ret |= tagInchiDiffBits_IDIF_LESS_H;
    }

    if i1.lenConnTable != i2.lenConnTable {
        ret |= tagInchiDiffBits_IDIF_CON_LEN;
        finish!();
    }
    if i1.lenConnTable > 0 {
        let connection_count =
            usize::try_from(i1.lenConnTable).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let connections1 = heap
            .slice(i1.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let connections2 = heap
            .slice(i2.nConnTable.as_const())?
            .get(..connection_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if connections1 != connections2 {
            ret |= tagInchiDiffBits_IDIF_CON_TBL;
            finish!();
        }
    }

    let num_groups1 = if i1.lenTautomer != 0 {
        *heap
            .slice(i1.nTautomer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)? as i32
    } else {
        0
    };
    let num_groups2 = if i2.lenTautomer != 0 {
        *heap
            .slice(i2.nTautomer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)? as i32
    } else {
        0
    };
    if num_groups1 == 0 && num_groups2 != 0 {
        ret |= tagInchiDiffBits_IDIF_NO_TAUT;
    } else if num_groups1 != 0 && num_groups2 == 0 {
        ret |= tagInchiDiffBits_IDIF_WRONG_TAUT;
    } else if num_groups1 == 1 && num_groups2 > 1 {
        ret |= tagInchiDiffBits_IDIF_SINGLE_TG;
    } else if num_groups1 > 1 && num_groups2 == 1 {
        ret |= tagInchiDiffBits_IDIF_MULTIPLE_TG;
    } else if num_groups1 != num_groups2 {
        ret |= tagInchiDiffBits_IDIF_NUM_TG;
    }

    if num_groups1 != 0 || num_groups2 != 0 {
        let byte_count1 = u64::try_from(
            (i64::from(i1.lenTautomer) + 1)
                .checked_mul(std::mem::size_of::<AT_NUMB>() as i64)
                .ok_or(SourceHeapError::AllocationSizeOverflow)?,
        )
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        let byte_count2 = u64::try_from(
            (i64::from(i2.lenTautomer) + 1)
                .checked_mul(std::mem::size_of::<AT_NUMB>() as i64)
                .ok_or(SourceHeapError::AllocationSizeOverflow)?,
        )
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;

        let endpoint_allocation1 = match inchi_malloc(heap, byte_count1) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                *err = -1;
                finish!();
            }
            Err(error) => return Err(error),
        };
        let endpoint_allocation2 = match inchi_malloc(heap, byte_count2) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                inchi_free(heap, endpoint_allocation1)?;
                *err = -1;
                finish!();
            }
            Err(error) => {
                inchi_free(heap, endpoint_allocation1)?;
                return Err(error);
            }
        };

        let endpoint_result =
            (|| -> Result<(Vec<AT_NUMB>, Vec<AT_NUMB>, i32, i32, i32, i32), SourceHeapError> {
                let tautomer_count1 = usize::try_from(i1.lenTautomer)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let tautomer_count2 = usize::try_from(i2.lenTautomer)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let tautomer1 = if tautomer_count1 == 0 {
                    Vec::new()
                } else {
                    heap.slice(i1.nTautomer.as_const())?
                        .get(..tautomer_count1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec()
                };
                let tautomer2 = if tautomer_count2 == 0 {
                    Vec::new()
                } else {
                    heap.slice(i2.nTautomer.as_const())?
                        .get(..tautomer_count2)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .to_vec()
                };

                let mut endpoints1 = Vec::new();
                let mut endpoints2 = Vec::new();
                let mut taut_h1 = 0_i32;
                let mut taut_h2 = 0_i32;
                let mut taut_m1 = 0_i32;
                let mut taut_m2 = 0_i32;

                let mut m = 1_usize;
                while m < tautomer_count1 {
                    let length = usize::from(tautomer1[m]);
                    m += 1;
                    let record = tautomer1
                        .get(m..m.saturating_add(length))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    taut_h1 = taut_h1.wrapping_add(i32::from(
                        *record.first().ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ));
                    taut_m1 = taut_m1.wrapping_add(i32::from(
                        *record.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ));
                    endpoints1.extend_from_slice(
                        record
                            .get(2..length)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    m = m
                        .checked_add(length)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                }

                m = 1;
                while m < tautomer_count2 {
                    let length = usize::from(tautomer2[m]);
                    m += 1;
                    let record = tautomer2
                        .get(m..m.saturating_add(length))
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    taut_h2 = taut_h2.wrapping_add(i32::from(
                        *record.first().ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ));
                    taut_m2 = taut_m2.wrapping_add(i32::from(
                        *record.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?,
                    ));
                    endpoints2.extend_from_slice(
                        record
                            .get(2..length)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    m = m
                        .checked_add(length)
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                }

                let count1 = i32::try_from(endpoints1.len())
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let count2 = i32::try_from(endpoints2.len())
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                insertions_sort_AT_RANK(&mut endpoints1, count1)?;
                insertions_sort_AT_RANK(&mut endpoints2, count2)?;
                Ok((endpoints1, endpoints2, taut_h1, taut_h2, taut_m1, taut_m2))
            })();

        let free1 = inchi_free(heap, endpoint_allocation1);
        let free2 = inchi_free(heap, endpoint_allocation2);
        free1?;
        free2?;
        let (endpoints1, endpoints2, taut_h1, taut_h2, taut_m1, taut_m2) = endpoint_result?;

        picr.num_taut_H1 = taut_h1;
        picr.num_taut_H2 = taut_h2;
        picr.num_taut_M1 = taut_m1;
        picr.num_taut_M2 = taut_m2;

        let mut index1 = 0_usize;
        let mut index2 = 0_usize;
        let mut num_equal = 0_i32;
        let mut num_in1_only = 0_i32;
        let mut num_in2_only = 0_i32;
        while index1 < endpoints1.len() && index2 < endpoints2.len() {
            if endpoints1[index1] == endpoints2[index2] {
                index1 += 1;
                index2 += 1;
                num_equal = num_equal.wrapping_add(1);
            } else if endpoints1[index1] < endpoints2[index2] {
                if picr.num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY as i32 {
                    picr.endp_in1_only[picr.num_endp_in1_only as usize] = endpoints1[index1];
                    picr.num_endp_in1_only = picr.num_endp_in1_only.wrapping_add(1);
                }
                index1 += 1;
                num_in1_only = num_in1_only.wrapping_add(1);
            } else {
                if picr.num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY as i32 {
                    picr.endp_in2_only[picr.num_endp_in2_only as usize] = endpoints2[index2];
                    picr.num_endp_in2_only = picr.num_endp_in2_only.wrapping_add(1);
                }
                index2 += 1;
                num_in2_only = num_in2_only.wrapping_add(1);
            }
        }
        while index1 < endpoints1.len() {
            if picr.num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY as i32 {
                picr.endp_in1_only[picr.num_endp_in1_only as usize] = endpoints1[index1];
                picr.num_endp_in1_only = picr.num_endp_in1_only.wrapping_add(1);
            }
            index1 += 1;
            num_in1_only = num_in1_only.wrapping_add(1);
        }
        while index2 < endpoints2.len() {
            if picr.num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY as i32 {
                picr.endp_in2_only[picr.num_endp_in2_only as usize] = endpoints2[index2];
                picr.num_endp_in2_only = picr.num_endp_in2_only.wrapping_add(1);
            }
            index2 += 1;
            num_in2_only = num_in2_only.wrapping_add(1);
        }
        if num_in1_only != 0 {
            ret |= tagInchiDiffBits_IDIF_EXTRA_TG_ENDP;
        }
        if num_in2_only != 0 {
            ret |= tagInchiDiffBits_IDIF_MISS_TG_ENDP;
        }
        if num_in1_only != 0 || num_in2_only != 0 || num_equal == 0 {
            ret |= tagInchiDiffBits_IDIF_DIFF_TG_ENDP;
        }
    }

    if i1.lenTautomer > 1 && i2.lenTautomer > 1 {
        let raw_tautomer_difference = if i1.lenTautomer != i2.lenTautomer {
            true
        } else {
            let count = usize::try_from(i1.lenTautomer)
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            heap.slice(i1.nTautomer.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                != heap
                    .slice(i2.nTautomer.as_const())?
                    .get(..count)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
        };
        if raw_tautomer_difference {
            ret |= tagInchiDiffBits_IDIF_TG;
        }
    }

    if i1.nNumberOfIsotopicAtoms != i2.nNumberOfIsotopicAtoms {
        ret |= tagInchiDiffBits_IDIF_NUM_ISO_AT;
    } else if i1.nNumberOfIsotopicAtoms > 0 {
        let count = usize::try_from(i1.nNumberOfIsotopicAtoms)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if heap
            .slice(i1.IsotopicAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != heap
                .slice(i2.IsotopicAtom.as_const())?
                .get(..count)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
        {
            ret |= tagInchiDiffBits_IDIF_ISO_AT;
        }
    }
    if i1.nTotalCharge != i2.nTotalCharge {
        ret |= tagInchiDiffBits_IDIF_CHARGE;
    }
    if let Some(aux1) = a1 {
        if aux1.nNumRemovedProtons != 0
            && a2.is_none_or(|aux2| aux2.nNumRemovedProtons != aux1.nNumRemovedProtons)
        {
            ret |= tagInchiDiffBits_IDIF_REM_PROT;
        }
        if a2.is_none_or(|aux2| {
            aux2.nNumRemovedIsotopicH[0] != aux1.nNumRemovedIsotopicH[0]
                || aux2.nNumRemovedIsotopicH[1] != aux1.nNumRemovedIsotopicH[1]
                || aux2.nNumRemovedIsotopicH[2] != aux1.nNumRemovedIsotopicH[2]
        }) {
            ret |= tagInchiDiffBits_IDIF_REM_ISO_H;
        }
    }

    let isotopic_stereo1 = if i1.StereoIsotopic.is_null() {
        None
    } else {
        heap.slice(i1.StereoIsotopic.as_const())?.first()
    };
    let regular_stereo1 = if i1.Stereo.is_null() {
        None
    } else {
        heap.slice(i1.Stereo.as_const())?.first()
    };
    let stereo1 = if isotopic_stereo1.is_some_and(|stereo| {
        stereo
            .nNumberOfStereoBonds
            .wrapping_add(stereo.nNumberOfStereoCenters)
            != 0
    }) {
        isotopic_stereo1
    } else {
        regular_stereo1
    };

    let isotopic_stereo2 = if i2.StereoIsotopic.is_null() {
        None
    } else {
        heap.slice(i2.StereoIsotopic.as_const())?.first()
    };
    let regular_stereo2 = if i2.Stereo.is_null() {
        None
    } else {
        heap.slice(i2.Stereo.as_const())?.first()
    };
    let stereo2 = if isotopic_stereo2.is_some_and(|stereo| {
        stereo
            .nNumberOfStereoBonds
            .wrapping_add(stereo.nNumberOfStereoCenters)
            != 0
    }) {
        isotopic_stereo2
    } else {
        regular_stereo2
    };

    ret |= CompareReversedStereoINChI2(heap, stereo1, stereo2, picr)? as u32;
    picr.flags = INCHI_MODE::from(ret);
    Ok(INCHI_MODE::from(ret))
}

#[allow(non_snake_case)]
pub(crate) fn CompareReversedINChI(
    heap: &SourceHeap,
    i1: Option<&INChI>,
    i2: Option<&INChI>,
    a1: Option<&INChI_Aux>,
    a2: Option<&INChI_Aux>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2936 CompareReversedINChI
    // INCHI✔️❌: int CompareReversedINChI(INChI* i1 /* InChI from reversed struct */,
    // INCHI✔️❌:     INChI* i2 /* input InChI */,
    // INCHI✔️❌:     INChI_Aux* a1,
    // INCHI✔️❌:     INChI_Aux* a2)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int ret;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1 == NULL && i2 == NULL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((i1 == NULL) ^ (i2 == NULL))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1; /* Diff: Missing InChI */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1->nErrorCode == i2->nErrorCode)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (i1->nErrorCode)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 2; /* Diff: Error codes */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i1->bDeleted != i2->bDeleted)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1; /* Diff: Missing InChI */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i1->nNumberOfAtoms != i2->nNumberOfAtoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 3;  /* Diff: Num. atoms */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1->nNumberOfAtoms > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (memcmp(i1->nAtom, i2->nAtom, i1->nNumberOfAtoms * sizeof(i1->nAtom[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 4; /* Diff: Elements */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (strcmp(i1->szHillFormula, i2->szHillFormula))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 7; /* Diff: Hill Formulas */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (memcmp(i1->nNum_H, i2->nNum_H, i1->nNumberOfAtoms * sizeof(i1->nNum_H[0])))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (i1->lenConnTable > 1 || i2->lenConnTable > 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 5; /* Diff: H Locations (mobile H present) */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 6; /* Diff: H Locations (no mobile H) */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* fixed H */
    // INCHI✔️❌:         if (i1->nNum_H_fixed || i2->nNum_H_fixed)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int bHasFixedH1 = 0, bHasFixedH2 = 0, i, j1, j2;
    // INCHI✔️❌:             if (i1->nNum_H_fixed)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (i1->nNum_H_fixed[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bHasFixedH1++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (i2->nNum_H_fixed)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = 0; i < i2->nNumberOfAtoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (i2->nNum_H_fixed[i])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bHasFixedH2++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* count the differences */
    // INCHI✔️❌:             j1 = j2 = 0;
    // INCHI✔️❌:             if (bHasFixedH1 && !bHasFixedH2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (i = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (i1->nNum_H_fixed[i] > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         j1++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                         if (i1->nNum_H_fixed[i] < 0)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             j2++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 return 18; /* Diff: Extra Fixed-H */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (!bHasFixedH1 && bHasFixedH2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (i = j1 = j2 = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (0 > i2->nNum_H_fixed[i])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             j1++;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                             if (0 < i2->nNum_H_fixed[i])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 j2++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     return 19; /* Diff: Missed Fixed-H */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (bHasFixedH1 && bHasFixedH2 &&
    // INCHI✔️❌:                         memcmp(i1->nNum_H_fixed, i2->nNum_H_fixed, i1->nNumberOfAtoms * sizeof(i1->nNum_H_fixed[0])))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         for (i = j1 = j2 = 0; i < i1->nNumberOfAtoms; i++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (i1->nNum_H_fixed[i] > i2->nNum_H_fixed[i])
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 j1++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (i1->nNum_H_fixed[i] < i2->nNum_H_fixed[i])
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     j2++;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             ret = (j1 && j2) ? 20 : j1 ? 18 : j2 ? 19 : 0;
    // INCHI✔️❌:             if (ret)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return ret; /* 20 => Diff: NotEql Fixed-H */
    // INCHI✔️❌:                 /* 19 => Diff: Missed Fixed-H (i1 has less) */
    // INCHI✔️❌:                 /* 18 => Diff: Extra Fixed-H  (i1 has more) */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1->lenConnTable != i2->lenConnTable)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 8; /* Diff: Connections length */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i1->lenConnTable > 0 && memcmp(i1->nConnTable, i2->nConnTable, i1->lenConnTable * sizeof(i1->nConnTable[0])))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 9; /* Diff: Connections */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* output special cases: different number of t-groups, different sizes of t-groups, different endpoints */
    // INCHI✔️❌:     if (i1->lenTautomer != i2->lenTautomer && (i1->lenTautomer > 1 || i2->lenTautomer > 1))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 10; /* Diff: Mobile groups length */ /* in isotopic or deprotonated cases i1->lenTautomer == 1 && i1->nTautomer[0] = 0 */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if ((i1->lenTautomer > 1 && i2->lenTautomer > 1) &&
    // INCHI✔️❌:         memcmp(i1->nTautomer, i2->nTautomer, i1->lenTautomer * sizeof(i1->nTautomer[0])))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 11; /* Diff: Mobile groups */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (i1->nNumberOfIsotopicAtoms != i2->nNumberOfIsotopicAtoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 12; /* Diff: Isotopic atoms number */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i1->nNumberOfIsotopicAtoms > 0 && memcmp(i1->IsotopicAtom, i2->IsotopicAtom, i1->nNumberOfIsotopicAtoms * sizeof(i1->IsotopicAtom[0])))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 13; /* Diff: Isotopic atoms */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (i1->nTotalCharge != i2->nTotalCharge)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 14; /* Diff: Charge */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         if ( i1->nNumberOfIsotopicTGroups != i2->nNumberOfIsotopicTGroups )
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 14;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if ( i1->nNumberOfIsotopicTGroups > 0 && memcmp( i1->IsotopicTGroup, i2->IsotopicTGroup, i1->nNumberOfIsotopicTGroups*sizeof(i1->IsotopicTGroup[0]) ) )
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 15;
    // INCHI✔️❌:             }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (a1 && a2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (a1->nNumRemovedProtons != a2->nNumRemovedProtons)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 16; /* Diff: Number of removed protons */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (memcmp(a1->nNumRemovedIsotopicH, a2->nNumRemovedIsotopicH, sizeof(a1->nNumRemovedIsotopicH)))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 17; /* Diff: Removed isotopic H */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         if ( i1->nPossibleLocationsOfIsotopicH && i2->nPossibleLocationsOfIsotopicH ) {
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ( i1->nPossibleLocationsOfIsotopicH[0] != i2->nPossibleLocationsOfIsotopicH[0] ||
    // INCHI✔️❌:                  memcmp(i1->nPossibleLocationsOfIsotopicH, i2->nPossibleLocationsOfIsotopicH,
    // INCHI✔️❌:                         sizeof(i1->nPossibleLocationsOfIsotopicH[0])*i1->nPossibleLocationsOfIsotopicH[0]) )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 18;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         } else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ( !i1->nPossibleLocationsOfIsotopicH != !i2->nPossibleLocationsOfIsotopicH ) {
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 19;}
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     /* ret = 20..31 => 40..51 */
    // INCHI✔️❌:     if ((ret = CompareReversedStereoINChI(i1->Stereo, i2->Stereo))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return ret + 20;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* ret = 40..51 => 60..71 */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!i2->StereoIsotopic && i2->Stereo && i1->StereoIsotopic &&
    // INCHI✔️❌:         0 < (i1->StereoIsotopic->nNumberOfStereoBonds + i1->StereoIsotopic->nNumberOfStereoCenters) &&
    // INCHI✔️❌:         0 == CompareReversedStereoINChI(i1->StereoIsotopic, i2->Stereo))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* InChI from reversed structure does not contain fully duplicated isotopic stereo */
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if ((ret = CompareReversedStereoINChI(i1->StereoIsotopic, i2->StereoIsotopic))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return ret + 40;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CompareReversedINChI

    let (Some(i1), Some(i2)) = (i1, i2) else {
        return Ok(if i1.is_none() && i2.is_none() { 0 } else { 1 });
    };
    if i1.nErrorCode == i2.nErrorCode {
        if i1.nErrorCode != 0 {
            return Ok(0);
        }
    } else {
        return Ok(2);
    }
    if i1.bDeleted != i2.bDeleted {
        return Ok(1);
    }
    if i1.nNumberOfAtoms != i2.nNumberOfAtoms {
        return Ok(3);
    }
    if i1.nNumberOfAtoms > 0 {
        let count = usize::try_from(i1.nNumberOfAtoms)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let atoms1 = heap
            .slice(i1.nAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let atoms2 = heap
            .slice(i2.nAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if atoms1 != atoms2 {
            return Ok(4);
        }
        let formula1 = heap.slice(i1.szHillFormula.as_const())?;
        let formula2 = heap.slice(i2.szHillFormula.as_const())?;
        let end1 = formula1
            .iter()
            .position(|value| *value == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let end2 = formula2
            .iter()
            .position(|value| *value == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        if formula1[..end1] != formula2[..end2] {
            return Ok(7);
        }
        let hydrogens1 = heap
            .slice(i1.nNum_H.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let hydrogens2 = heap
            .slice(i2.nNum_H.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if hydrogens1 != hydrogens2 {
            return Ok(if i1.lenConnTable > 1 || i2.lenConnTable > 1 {
                5
            } else {
                6
            });
        }
        if !i1.nNum_H_fixed.is_null() || !i2.nNum_H_fixed.is_null() {
            let fixed1 = if i1.nNum_H_fixed.is_null() {
                None
            } else {
                Some(
                    heap.slice(i1.nNum_H_fixed.as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let fixed2 = if i2.nNum_H_fixed.is_null() {
                None
            } else {
                Some(
                    heap.slice(i2.nNum_H_fixed.as_const())?
                        .get(..count)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
            };
            let has1 = fixed1.map_or(0, |values| {
                values.iter().filter(|value| **value != 0).count()
            });
            let has2 = fixed2.map_or(0, |values| {
                values.iter().filter(|value| **value != 0).count()
            });
            if has1 != 0 && has2 == 0 {
                return Ok(18);
            }
            if has1 == 0 && has2 != 0 {
                return Ok(19);
            }
            if has1 != 0 && has2 != 0 {
                let fixed1 = fixed1.ok_or(SourceHeapError::NullPointer)?;
                let fixed2 = fixed2.ok_or(SourceHeapError::NullPointer)?;
                if fixed1 != fixed2 {
                    let mut greater = false;
                    let mut less = false;
                    for (&left, &right) in fixed1.iter().zip(fixed2) {
                        greater |= left > right;
                        less |= left < right;
                    }
                    let ret = if greater && less {
                        20
                    } else if greater {
                        18
                    } else if less {
                        19
                    } else {
                        0
                    };
                    if ret != 0 {
                        return Ok(ret);
                    }
                }
            }
        }
    }
    if i1.lenConnTable != i2.lenConnTable {
        return Ok(8);
    }
    if i1.lenConnTable > 0 {
        let count =
            usize::try_from(i1.lenConnTable).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let connections1 = heap
            .slice(i1.nConnTable.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let connections2 = heap
            .slice(i2.nConnTable.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if connections1 != connections2 {
            return Ok(9);
        }
    }
    if i1.lenTautomer != i2.lenTautomer && (i1.lenTautomer > 1 || i2.lenTautomer > 1) {
        return Ok(10);
    }
    if i1.lenTautomer > 1 && i2.lenTautomer > 1 {
        let count =
            usize::try_from(i1.lenTautomer).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let tautomer1 = heap
            .slice(i1.nTautomer.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let tautomer2 = heap
            .slice(i2.nTautomer.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if tautomer1 != tautomer2 {
            return Ok(11);
        }
    }
    if i1.nNumberOfIsotopicAtoms != i2.nNumberOfIsotopicAtoms {
        return Ok(12);
    }
    if i1.nNumberOfIsotopicAtoms > 0 {
        let count = usize::try_from(i1.nNumberOfIsotopicAtoms)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let isotopic_atoms1 = heap
            .slice(i1.IsotopicAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let isotopic_atoms2 = heap
            .slice(i2.IsotopicAtom.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if isotopic_atoms1 != isotopic_atoms2 {
            return Ok(13);
        }
    }
    if i1.nTotalCharge != i2.nTotalCharge {
        return Ok(14);
    }
    if let (Some(a1), Some(a2)) = (a1, a2) {
        if a1.nNumRemovedProtons != a2.nNumRemovedProtons {
            return Ok(16);
        }
        if a1.nNumRemovedIsotopicH != a2.nNumRemovedIsotopicH {
            return Ok(17);
        }
    }
    let stereo1 = if i1.Stereo.is_null() {
        None
    } else {
        Some(&heap.slice(i1.Stereo.as_const())?[0])
    };
    let stereo2 = if i2.Stereo.is_null() {
        None
    } else {
        Some(&heap.slice(i2.Stereo.as_const())?[0])
    };
    let ret = CompareReversedStereoINChI(heap, stereo1, stereo2)?;
    if ret != 0 {
        return Ok(ret + 20);
    }
    let isotopic1 = if i1.StereoIsotopic.is_null() {
        None
    } else {
        Some(&heap.slice(i1.StereoIsotopic.as_const())?[0])
    };
    let isotopic2 = if i2.StereoIsotopic.is_null() {
        None
    } else {
        Some(&heap.slice(i2.StereoIsotopic.as_const())?[0])
    };
    if isotopic2.is_none()
        && stereo2.is_some()
        && isotopic1.is_some()
        && isotopic1
            .is_some_and(|stereo| stereo.nNumberOfStereoBonds + stereo.nNumberOfStereoCenters > 0)
        && CompareReversedStereoINChI(heap, isotopic1, stereo2)? == 0
    {
        return Ok(0);
    }
    let ret = CompareReversedStereoINChI(heap, isotopic1, isotopic2)?;
    Ok(if ret != 0 { ret + 40 } else { 0 })
}

#[allow(non_snake_case)]
pub(crate) fn GetInpStructErrorType(
    input_parameters: Option<&INPUT_PARMS>,
    error: i32,
    structure_error: Option<&[i8]>,
    num_input_atoms: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2480 GetInpStructErrorType
    // INCHI✔️❌: int GetInpStructErrorType(INPUT_PARMS* ip,
    // INCHI✔️❌:     int err,
    // INCHI✔️❌:     char* pStrErrStruct,
    // INCHI✔️❌:     int num_inp_atoms)
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (err && err == 9)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return _IS_ERROR; /*  sdfile bypassed to $$$$ */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (err && err < 30)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return _IS_FATAL;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (num_inp_atoms <= 0 || err)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (98 == err && 0 == num_inp_atoms && ip->bAllowEmptyStructure)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return _IS_WARNING;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         return _IS_ERROR;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (pStrErrStruct[0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return _IS_WARNING;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return _IS_OKAY;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: GetInpStructErrorType

    if error != 0 && error == 9 {
        return Ok(_IS_ERROR as i32);
    }
    if error != 0 && error < 30 {
        return Ok(_IS_FATAL as i32);
    }
    if num_input_atoms <= 0 || error != 0 {
        if error == 98
            && num_input_atoms == 0
            && input_parameters
                .ok_or(SourceHeapError::NullPointer)?
                .bAllowEmptyStructure
                != 0
        {
            return Ok(_IS_WARNING as i32);
        }
        return Ok(_IS_ERROR as i32);
    }
    let first = structure_error
        .and_then(|error| error.first())
        .ok_or(SourceHeapError::NullPointer)?;
    if *first != 0 {
        return Ok(_IS_WARNING as i32);
    }
    Ok(_IS_OKAY as i32)
}

pub(crate) fn mystrrev(
    heap: &mut SourceHeap,
    string: SourceMutPointer<i8>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2090 mystrrev
    // INCHI✔️✔️: void mystrrev(char* p)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     char c, * q = p;
    // INCHI✔️✔️:     while (*q++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         ;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     q -= 2; /*  pointer to the last character */
    // INCHI✔️✔️:     while (p < q)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         c = *q;  /*  swap */
    // INCHI✔️✔️:         *q-- = *p;
    // INCHI✔️✔️:         *p++ = c;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: mystrrev

    let bytes = heap.slice_mut(string)?;
    let length = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    bytes[..length].reverse();
    Ok(())
}

struct OrderStruct<'a> {
    dfs_number: &'a [AT_NUMB],
    number_of_descendants: &'a [AT_NUMB],
    current_atom: i32,
}

#[allow(non_snake_case)]
fn CompareDfsDescendants4CT(
    neighbor1: AT_RANK,
    neighbor2: AT_RANK,
    order: &OrderStruct<'_>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2119 CompareDfsDescendants4CT
    // INCHI✔️✔️: static int CompareDfsDescendants4CT(const void* a1, const void* a2, void* p)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     OrderStruct* os = (OrderStruct*)p;
    // INCHI✔️✔️:     int neigh1 = (int)*(const AT_RANK*)a1;
    // INCHI✔️✔️:     int neigh2 = (int)*(const AT_RANK*)a2;
    // INCHI✔️✔️:     if (neigh1 > MAX_ATOMS)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (neigh2 > MAX_ATOMS)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 0;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else if (neigh2 > MAX_ATOMS)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         AT_RANK nCurDfsNumber = os->m_gDfs4CT_nDfsNumber[os->m_gDfs4CT_nCurrentAtom];
    // INCHI✔️✔️:         int nDesc1 = nCurDfsNumber > os->m_gDfs4CT_nDfsNumber[neigh1] ?
    // INCHI✔️✔️:             0 : (int)os->m_gDfs4CT_nNumDescendants[neigh1];
    // INCHI✔️✔️:         int nDesc2 = nCurDfsNumber > os->m_gDfs4CT_nDfsNumber[neigh2] ?
    // INCHI✔️✔️:             0 : (int)os->m_gDfs4CT_nNumDescendants[neigh2];
    // INCHI✔️✔️:         int ret;
    // INCHI✔️✔️:         if ((ret = nDesc1 - nDesc2)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return ret;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         return  (int)neigh1 - (int)neigh2; /*  canon. numbers difference */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: CompareDfsDescendants4CT

    if u32::from(neighbor1) > MAX_ATOMS {
        return Ok(if u32::from(neighbor2) > MAX_ATOMS {
            0
        } else {
            1
        });
    }
    if u32::from(neighbor2) > MAX_ATOMS {
        return Ok(-1);
    }
    let current =
        usize::try_from(order.current_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current_dfs = *order
        .dfs_number
        .get(current)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let first = usize::from(neighbor1);
    let second = usize::from(neighbor2);
    let descendants1 = if current_dfs
        > *order
            .dfs_number
            .get(first)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
    {
        0
    } else {
        i32::from(
            *order
                .number_of_descendants
                .get(first)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let descendants2 = if current_dfs
        > *order
            .dfs_number
            .get(second)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
    {
        0
    } else {
        i32::from(
            *order
                .number_of_descendants
                .get(second)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        )
    };
    let result = descendants1 - descendants2;
    if result != 0 {
        return Ok(result);
    }
    Ok(i32::from(neighbor1) - i32::from(neighbor2))
}

#[allow(non_snake_case)]
pub(crate) fn GetDfsOrder4CT(
    heap: &mut SourceHeap,
    linear_ct: SourceMutPointer<AT_NUMB>,
    length_ct: i32,
    number_hydrogens: SourceConstPointer<i8>,
    number_of_atoms: i32,
    ct_mode: i32,
) -> Result<SourceMutPointer<AT_NUMB>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2154 GetDfsOrder4CT
    // INCHI✔️✔️: AT_NUMB* GetDfsOrder4CT(CANON_GLOBALS* pCG,
    // INCHI✔️✔️:     AT_NUMB* LinearCT,
    // INCHI✔️✔️:     int nLenCT,
    // INCHI✔️✔️:     S_CHAR* nNum_H,
    // INCHI✔️✔️:     int num_atoms,
    // INCHI✔️✔️:     int nCtMode)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     AT_NUMB* nStackAtom = NULL;
    // INCHI✔️✔️:     int         nTopStackAtom = -1;
    // INCHI✔️✔️:     AT_NUMB* nNumDescendants = NULL; /*  number of descendants incl. closures and the atom itself */
    // INCHI✔️✔️:     AT_NUMB* nDfsNumber = NULL;
    // INCHI✔️✔️:     S_CHAR* cNeighNumb = NULL;
    // INCHI✔️✔️:     NEIGH_LIST* nl = NULL;
    // INCHI✔️✔️:     AT_NUMB     nDfs;
    // INCHI✔️✔️:     int         i, j, u, k, start, num_rings, nTotOutputStringLen;
    // INCHI✔️✔️:     AT_NUMB* nOutputString = NULL, cDelim;
    // INCHI✔️✔️:     int         bCtPredecessors = (nCtMode & CT_MODE_PREDECESSORS);
    // INCHI✔️✔️:     OrderStruct os;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*  allocate arrays */
    // INCHI✔️✔️:     nStackAtom = (AT_NUMB*)inchi_malloc(num_atoms * sizeof(nStackAtom[0]));
    // INCHI✔️✔️:     nNumDescendants = (AT_NUMB*)inchi_malloc(num_atoms * sizeof(nNumDescendants[0]));
    // INCHI✔️✔️:     nDfsNumber = (AT_NUMB*)inchi_malloc(num_atoms * sizeof(nDfsNumber[0]));
    // INCHI✔️✔️:     cNeighNumb = (S_CHAR*)inchi_malloc(num_atoms * sizeof(cNeighNumb[0]));
    // INCHI✔️✔️:     nl = CreateNeighListFromLinearCT(LinearCT, nLenCT, num_atoms);
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*  check allocation */
    // INCHI✔️✔️:     if (!nStackAtom || !nNumDescendants || !nDfsNumber || !cNeighNumb || !nl)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* ret = CT_OUT_OF_RAM; */ /*  program error */ /*   <BRKPT> */
    // INCHI✔️✔️:         goto exit_function;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (bCtPredecessors)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         start = 0;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /*  find DFS start vertex (atom) */
    // INCHI✔️✔️:         for (i = 1, start = 0; i < num_atoms; i++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             if (nl[i][0] < nl[start][0])
    // INCHI✔️✔️:             { /*  index = nRank-1 */
    // INCHI✔️✔️:                 start = i;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /*
    // INCHI✔️✔️:       vertex information:
    // INCHI✔️✔️:         1. Number of (forward edges) + (back edges, first visit -- ring closures): nl[i][0]
    // INCHI✔️✔️:         2. Number of vertices traversed from this vertex, including the vertex:    nNumDescendants[i]
    // INCHI✔️✔️:         3. Each edge information:
    // INCHI✔️✔️:            a. forward edge (0) or back edge (1) indicator: nDfsNumber[i] > nDfsNumber[neigh]
    // INCHI✔️✔️:            b. neighbor at another end of the edge neigh = nl[i][k+1], k < i
    // INCHI✔️✔️:
    // INCHI✔️✔️:         Total per edge: 2 + 2*(number of edges)
    // INCHI✔️✔️:     */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* DFS initiation */
    // INCHI✔️✔️:     u = start; /* start atom */
    // INCHI✔️✔️:     nDfs = 0;
    // INCHI✔️✔️:     nTopStackAtom = -1;
    // INCHI✔️✔️:     memset(nDfsNumber, 0, num_atoms * sizeof(nDfsNumber[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     memset(nNumDescendants, 0, num_atoms * sizeof(nNumDescendants[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     memset(cNeighNumb, 0, num_atoms * sizeof(cNeighNumb[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:     /*  push the start atom on the stack */
    // INCHI✔️✔️:     nDfsNumber[u] = ++nDfs;
    // INCHI✔️✔️:     if (bCtPredecessors)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nNumDescendants[u] = 0; /* atom #1 has no predecessor */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         nNumDescendants[u] = 1; /* count itself as a descendant */
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    // INCHI✔️✔️:     /* nNumStartChildren = 0; */
    // INCHI✔️✔️:     num_rings = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* DFS */
    // INCHI✔️✔️:     do
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         /* advance */
    // INCHI✔️✔️:         while (i = (int)nStackAtom[nTopStackAtom], j = (int)cNeighNumb[i] + 1, (int)nl[i][0] >= j)
    // INCHI✔️✔️:             /*while ( (int)nl[i=nStackAtom[nTopStackAtom]][0] >= (j = (int)cNeighNumb[i]+1) )*/
    // INCHI✔️✔️:
    // INCHI✔️✔️:             /* replaced due to missing sequence point; undefined behavior, pointed by Geoffrey Hutchison */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             cNeighNumb[i]++;
    // INCHI✔️✔️:             u = (int)nl[i][j]; /*  jth neighbor of the vertex i */
    // INCHI✔️✔️:             if (!nDfsNumber[u])
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* tree edge, 1st visit -- advance */
    // INCHI✔️✔️:                 /* put unexplored vertex u on the stack for further examination */
    // INCHI✔️✔️:                 nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    // INCHI✔️✔️:                 nDfsNumber[u] = ++nDfs;
    // INCHI✔️✔️:                 if (bCtPredecessors)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     nNumDescendants[u] = i + 1; /* predecessor's rank */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     nNumDescendants[u]++; /* count atom u as its descendant */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 if (nTopStackAtom && u != (int)nStackAtom[nTopStackAtom - 1] &&
    // INCHI✔️✔️:                     /* back edge: u is not a predecessor of i */
    // INCHI✔️✔️:                     nDfsNumber[u] < nDfsNumber[i])
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* Back edge, 1st visit: u is an ancestor of i (ring closure) */
    // INCHI✔️✔️:                     if (!bCtPredecessors)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         nNumDescendants[i]++; /* count closures as descendants */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     num_rings++;          /* count ring closures */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 else
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     nl[i][j] = MAX_ATOMS + 1; /* back edge, 2nd visit: mark as deleted */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         cNeighNumb[i] = 0; /* all neighbors of the ith atom have been
    // INCHI✔️✔️:                               traversed; resore the neighbor counter */
    // INCHI✔️✔️:                               /* back up */
    // INCHI✔️✔️:         if (!bCtPredecessors && nTopStackAtom /* that is, i != start */)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             u = (int)nStackAtom[nTopStackAtom - 1]; /* predecessor of i */
    // INCHI✔️✔️:             nNumDescendants[u] += nNumDescendants[i]; /* add descendants */
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     } while (--nTopStackAtom >= 0);
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* Sort the neighbors in ascending order so that:
    // INCHI✔️✔️:        primary key   = number of descendants in the DFS tree; closure neighbor is 0
    // INCHI✔️✔️:        secondary key = canonical number (here vertex number = canonical number - 1)
    // INCHI✔️✔️:      */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     os.m_gDfs4CT_nDfsNumber = nDfsNumber;
    // INCHI✔️✔️:     os.m_gDfs4CT_nNumDescendants = nNumDescendants;
    // INCHI✔️✔️:     os.m_gDfs4CT_nCurrentAtom = -1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     /* sorting; deleted will be the last neighbors */
    // INCHI✔️✔️:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (nl[i][0] > 1)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             os.m_gDfs4CT_nCurrentAtom = i;
    // INCHI✔️✔️:             insertions_sort(&os, &nl[i][1], nl[i][0], sizeof(nl[i][1]), CompareDfsDescendants4CT);
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         /* reduce number of neighbors to exclude deleted */
    // INCHI✔️✔️:         for (k = 0; k < nl[i][0] && nl[i][k + 1] <= MAX_ATOMS; k++)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             ;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         nl[i][0] = k;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     nTotOutputStringLen = 3 * (num_atoms + num_rings + 1); /*  last 3 elements are a 'zero termination' */
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (bCtPredecessors)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if ((nOutputString = (AT_RANK*)inchi_calloc(nTotOutputStringLen, sizeof(nOutputString[0])))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             cDelim = '-';
    // INCHI✔️✔️:             for (u = 0, k = -3; u < num_atoms; u++)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 k += 3;
    // INCHI✔️✔️:                 if (k + 6 > nTotOutputStringLen)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     goto exit_error;  /* program error */
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:                 nOutputString[k] = nNumDescendants[u] ? nNumDescendants[u] : MAX_ATOMS + 1;
    // INCHI✔️✔️:                 nOutputString[k + 1] = nNum_H ? 16 + nNum_H[u] : 0;
    // INCHI✔️✔️:                 nOutputString[k + 2] = k ? ',' : '\0';
    // INCHI✔️✔️:                 for (j = 1; j <= nl[u][0] && nDfsNumber[u] > nDfsNumber[i = nl[u][j]]; j++)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     /* closures */
    // INCHI✔️✔️:                     k += 3;
    // INCHI✔️✔️:                     if (k + 6 > nTotOutputStringLen)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         goto exit_error;  /* program error */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     nOutputString[k] = i + 1;  /* closure */
    // INCHI✔️✔️:                     nOutputString[k + 1] = 0;
    // INCHI✔️✔️:                     nOutputString[k + 2] = cDelim;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     else
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (nNumDescendants)
    // INCHI✔️✔️:         {  /* do not need anymore */
    // INCHI✔️✔️:             inchi_free(nNumDescendants);
    // INCHI✔️✔️:             nNumDescendants = NULL;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:
    // INCHI✔️✔️:         /*
    // INCHI✔️✔️:             the output string contains:
    // INCHI✔️✔️:               (num_atoms) atoms for the DFS (spanning) tree
    // INCHI✔️✔️:               (num_atoms-1) delimiters for the DFS (spanning) tree
    // INCHI✔️✔️:               1 character for each atom that has 1 terminal hydrogen atoms
    // INCHI✔️✔️:               2 characters  for each atom that has 2-9 terminal hydrogen atoms
    // INCHI✔️✔️:               3 characters  for each atom that has 10-99 terminal hydrogen atoms, etc.
    // INCHI✔️✔️:               (num_rings) atoms for the ring closures
    // INCHI✔️✔️:               (num_rings) delimiters for the ring closures
    // INCHI✔️✔️:         */
    // INCHI✔️✔️:
    // INCHI✔️✔️:         if ((nOutputString = (AT_RANK*)inchi_calloc(nTotOutputStringLen, sizeof(nOutputString[0])))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             u = start; /*  start atom */
    // INCHI✔️✔️:             nTopStackAtom = -1;
    // INCHI✔️✔️:             memset(cNeighNumb, 0, num_atoms * sizeof(cNeighNumb[0])); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️✔️:             /*  push the start atom on the stack */
    // INCHI✔️✔️:             nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    // INCHI✔️✔️:             /*  output the starting atom */
    // INCHI✔️✔️:             k = 0;
    // INCHI✔️✔️:             nOutputString[k] = u + 1;
    // INCHI✔️✔️:             nOutputString[k + 1] = nNum_H ? 16 + nNum_H[u] : 0;
    // INCHI✔️✔️:             nOutputString[k + 2] = '\0';
    // INCHI✔️✔️:
    // INCHI✔️✔️:             do
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 /* advance */
    // INCHI✔️✔️:                 while (i = (int)nStackAtom[nTopStackAtom], j = (int)cNeighNumb[i] + 1, i < num_atoms && (int)nl[i][0] >= j) /* djb-rwth: additional condition to avoid buffer overruns */
    // INCHI✔️✔️:                     /*while ( (int)nl[i=nStackAtom[nTopStackAtom]][0] >= (j = (int)cNeighNumb[i]+1) )*/
    // INCHI✔️✔️:                     /* replaced due to missing sequence point; undefined behavior, reported by Geoffrey Hutchison */
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     k += 3;
    // INCHI✔️✔️:                     if (k + 6 > nTotOutputStringLen)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         goto exit_error;  /* program error */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     cNeighNumb[i]++;
    // INCHI✔️✔️:                     u = (int)nl[i][j]; /* neighbor */
    // INCHI✔️✔️:
    // INCHI✔️✔️:                     /* output neighbor's canonical number */
    // INCHI✔️✔️:                     nOutputString[k] = u + 1;
    // INCHI✔️✔️:
    // INCHI✔️✔️:                     if (nDfsNumber[u] > nDfsNumber[i])
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         /* tree edge, 1st visit -- advance */
    // INCHI✔️✔️:                         /* put 'unexplored' vertex u on the stack */
    // INCHI✔️✔️:                         nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
    // INCHI✔️✔️:
    // INCHI✔️✔️:                         /* output neighbor's number of H */
    // INCHI✔️✔️:                         nOutputString[k + 1] = nNum_H ? 16 + nNum_H[u] : 0;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         nOutputString[k + 1] = 0;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     /* output a delimiter preceding the neighbor */
    // INCHI✔️✔️:                     if (1 < nl[i][0])
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         if (j == 1)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             cDelim = '(';
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         else
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             if (j == nl[i][0])
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 cDelim = ')';
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                             else
    // INCHI✔️✔️:                             {
    // INCHI✔️✔️:                                 cDelim = ',';
    // INCHI✔️✔️:                             }
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     else
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         cDelim = '-';
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     nOutputString[k + 2] = cDelim;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 if ((i >= 0) && (i < num_atoms)) /* djb-rwth: fixing coverity ID #499483 */
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     cNeighNumb[i] = 0;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 /* back up: nothing else to do */
    // INCHI✔️✔️:             } while (--nTopStackAtom >= 0);
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     goto exit_function;
    // INCHI✔️✔️:
    // INCHI✔️✔️: exit_error:
    // INCHI✔️✔️:     if (nOutputString)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free(nOutputString);
    // INCHI✔️✔️:         nOutputString = NULL;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️: exit_function:
    // INCHI✔️✔️:     if (nStackAtom)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free(nStackAtom);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nNumDescendants)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free(nNumDescendants);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nDfsNumber)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free(nDfsNumber);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (cNeighNumb)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         inchi_free(cNeighNumb);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (nl)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         FreeNeighList(nl);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return nOutputString;
    // INCHI✔️✔️: }
    // INCHI✔️✔️:
    // END INCHI C FUNCTION: GetDfsOrder4CT

    let atom_count =
        usize::try_from(number_of_atoms).map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
    if atom_count == 0 {
        return Err(SourceHeapError::UnsupportedSourceBehavior);
    }
    let allocate_u16 = |heap: &mut SourceHeap| match inchi_calloc::<u16>(
        heap,
        atom_count as u64,
        std::mem::size_of::<u16>() as u64,
    ) {
        Ok(pointer) => Ok(pointer),
        Err(SourceHeapError::AllocationFailed) => Ok(SourceMutPointer::null()),
        Err(error) => Err(error),
    };
    let stack = allocate_u16(heap)?;
    let mut number_of_descendants = allocate_u16(heap)?;
    let dfs_number = allocate_u16(heap)?;
    let neighbor_number = match inchi_calloc::<i8>(heap, atom_count as u64, 1) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
        Err(error) => return Err(error),
    };
    let neighbor_lists =
        CreateNeighListFromLinearCT(heap, linear_ct.as_const(), length_ct, number_of_atoms)?;
    let predecessors = ct_mode & crate::source_types::CT_MODE_PREDECESSORS as i32;

    let computation = (|| -> Result<SourceMutPointer<AT_NUMB>, SourceHeapError> {
        if stack.is_null()
            || number_of_descendants.is_null()
            || dfs_number.is_null()
            || neighbor_number.is_null()
            || neighbor_lists.is_null()
        {
            return Ok(SourceMutPointer::null());
        }

        let mut start = 0_usize;
        if predecessors == 0 {
            for atom in 1..atom_count {
                let pointers = heap.slice(neighbor_lists.as_const())?;
                let valence = heap.slice(pointers[atom].as_const())?[0];
                let start_valence = heap.slice(pointers[start].as_const())?[0];
                if valence < start_valence {
                    start = atom;
                }
            }
        }

        heap.slice_mut(dfs_number)?.fill(0);
        heap.slice_mut(number_of_descendants)?.fill(0);
        heap.slice_mut(neighbor_number)?.fill(0);
        let mut dfs = 1_u16;
        heap.slice_mut(dfs_number)?[start] = dfs;
        heap.slice_mut(number_of_descendants)?[start] = if predecessors != 0 { 0 } else { 1 };
        heap.slice_mut(stack)?[0] = start as u16;
        let mut stack_top = 0_i32;
        let mut number_of_rings = 0_i32;

        loop {
            let mut atom = usize::from(heap.slice(stack.as_const())?[stack_top as usize]);
            loop {
                let neighbor_index = i32::from(heap.slice(neighbor_number.as_const())?[atom]) + 1;
                let list_pointer = heap.slice(neighbor_lists.as_const())?[atom];
                if i32::from(heap.slice(list_pointer.as_const())?[0]) < neighbor_index {
                    break;
                }
                heap.slice_mut(neighbor_number)?[atom] =
                    heap.slice(neighbor_number.as_const())?[atom].wrapping_add(1);
                let adjacent =
                    usize::from(heap.slice(list_pointer.as_const())?[neighbor_index as usize]);
                if heap.slice(dfs_number.as_const())?[adjacent] == 0 {
                    stack_top += 1;
                    heap.slice_mut(stack)?[stack_top as usize] = adjacent as u16;
                    dfs = dfs.wrapping_add(1);
                    heap.slice_mut(dfs_number)?[adjacent] = dfs;
                    if predecessors != 0 {
                        heap.slice_mut(number_of_descendants)?[adjacent] = (atom + 1) as u16;
                    } else {
                        heap.slice_mut(number_of_descendants)?[adjacent] =
                            heap.slice(number_of_descendants.as_const())?[adjacent].wrapping_add(1);
                    }
                } else if stack_top != 0
                    && adjacent
                        != usize::from(heap.slice(stack.as_const())?[(stack_top - 1) as usize])
                    && heap.slice(dfs_number.as_const())?[adjacent]
                        < heap.slice(dfs_number.as_const())?[atom]
                {
                    if predecessors == 0 {
                        heap.slice_mut(number_of_descendants)?[atom] =
                            heap.slice(number_of_descendants.as_const())?[atom].wrapping_add(1);
                    }
                    number_of_rings = number_of_rings.wrapping_add(1);
                } else {
                    heap.slice_mut(list_pointer)?[neighbor_index as usize] = (MAX_ATOMS + 1) as u16;
                }
                atom = usize::from(heap.slice(stack.as_const())?[stack_top as usize]);
            }
            heap.slice_mut(neighbor_number)?[atom] = 0;
            if predecessors == 0 && stack_top != 0 {
                let parent = usize::from(heap.slice(stack.as_const())?[(stack_top - 1) as usize]);
                heap.slice_mut(number_of_descendants)?[parent] = heap
                    .slice(number_of_descendants.as_const())?[parent]
                    .wrapping_add(heap.slice(number_of_descendants.as_const())?[atom]);
            }
            stack_top -= 1;
            if stack_top < 0 {
                break;
            }
        }

        for atom in 0..atom_count {
            let list_pointer = heap.slice(neighbor_lists.as_const())?[atom];
            let neighbor_count = usize::from(heap.slice(list_pointer.as_const())?[0]);
            if neighbor_count > 1 {
                heap.with_slice_mut_and_heap(list_pointer.offset(1)?, |neighbors, read_heap| {
                    let records = bytemuck::cast_slice_mut::<u16, u8>(
                        neighbors
                            .get_mut(..neighbor_count)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let order = OrderStruct {
                        dfs_number: read_heap.slice(dfs_number.as_const())?,
                        number_of_descendants: read_heap.slice(number_of_descendants.as_const())?,
                        current_atom: atom as i32,
                    };
                    insertions_sort(records, neighbor_count, 2, &mut |first, second| {
                        CompareDfsDescendants4CT(
                            u16::from_ne_bytes([first[0], first[1]]),
                            u16::from_ne_bytes([second[0], second[1]]),
                            &order,
                        )
                    })?;
                    Ok(())
                })?;
            }
            let list = heap.slice_mut(list_pointer)?;
            let mut retained = 0_usize;
            while retained < neighbor_count && u32::from(list[retained + 1]) <= MAX_ATOMS {
                retained += 1;
            }
            list[0] = retained as u16;
        }

        let total_length = number_of_atoms
            .checked_add(number_of_rings)
            .and_then(|value| value.checked_add(1))
            .and_then(|value| value.checked_mul(3))
            .ok_or(SourceHeapError::UnsupportedSourceBehavior)?;
        let output_length = usize::try_from(total_length)
            .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
        let output = match inchi_calloc::<u16>(heap, output_length as u64, 2) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => return Ok(SourceMutPointer::null()),
            Err(error) => return Err(error),
        };

        if predecessors != 0 {
            let mut position = 0_usize;
            for atom in 0..atom_count {
                if position + 6 > output_length {
                    inchi_free(heap, output)?;
                    return Ok(SourceMutPointer::null());
                }
                let predecessor = heap.slice(number_of_descendants.as_const())?[atom];
                heap.slice_mut(output)?[position] = if predecessor != 0 {
                    predecessor
                } else {
                    (MAX_ATOMS + 1) as u16
                };
                heap.slice_mut(output)?[position + 1] = if number_hydrogens.is_null() {
                    0
                } else {
                    (16_i32 + i32::from(heap.slice(number_hydrogens)?[atom])) as u16
                };
                heap.slice_mut(output)?[position + 2] = if position != 0 { b',' as u16 } else { 0 };
                let list_pointer = heap.slice(neighbor_lists.as_const())?[atom];
                let list_count = usize::from(heap.slice(list_pointer.as_const())?[0]);
                let mut list_index = 1_usize;
                while list_index <= list_count {
                    let adjacent = heap.slice(list_pointer.as_const())?[list_index];
                    if heap.slice(dfs_number.as_const())?[atom]
                        <= heap.slice(dfs_number.as_const())?[usize::from(adjacent)]
                    {
                        break;
                    }
                    position += 3;
                    if position + 6 > output_length {
                        inchi_free(heap, output)?;
                        return Ok(SourceMutPointer::null());
                    }
                    heap.slice_mut(output)?[position] = adjacent.wrapping_add(1);
                    heap.slice_mut(output)?[position + 1] = 0;
                    heap.slice_mut(output)?[position + 2] = b'-' as u16;
                    list_index += 1;
                }
                position += 3;
            }
        } else {
            inchi_free(heap, number_of_descendants)?;
            number_of_descendants = SourceMutPointer::null();
            heap.slice_mut(neighbor_number)?.fill(0);
            heap.slice_mut(stack)?[0] = start as u16;
            stack_top = 0;
            let mut position = 0_usize;
            heap.slice_mut(output)?[0] = (start + 1) as u16;
            heap.slice_mut(output)?[1] = if number_hydrogens.is_null() {
                0
            } else {
                (16_i32 + i32::from(heap.slice(number_hydrogens)?[start])) as u16
            };
            heap.slice_mut(output)?[2] = 0;
            loop {
                let atom = usize::from(heap.slice(stack.as_const())?[stack_top as usize]);
                let list_pointer = heap.slice(neighbor_lists.as_const())?[atom];
                let list_count = usize::from(heap.slice(list_pointer.as_const())?[0]);
                let neighbor_index =
                    usize::try_from(i32::from(heap.slice(neighbor_number.as_const())?[atom]))
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                        + 1;
                if atom < atom_count && neighbor_index <= list_count {
                    position += 3;
                    if position + 6 > output_length {
                        inchi_free(heap, output)?;
                        return Ok(SourceMutPointer::null());
                    }
                    heap.slice_mut(neighbor_number)?[atom] =
                        heap.slice(neighbor_number.as_const())?[atom].wrapping_add(1);
                    let adjacent =
                        usize::from(heap.slice(list_pointer.as_const())?[neighbor_index]);
                    heap.slice_mut(output)?[position] = (adjacent + 1) as u16;
                    if heap.slice(dfs_number.as_const())?[adjacent]
                        > heap.slice(dfs_number.as_const())?[atom]
                    {
                        stack_top += 1;
                        heap.slice_mut(stack)?[stack_top as usize] = adjacent as u16;
                        heap.slice_mut(output)?[position + 1] = if number_hydrogens.is_null() {
                            0
                        } else {
                            (16_i32 + i32::from(heap.slice(number_hydrogens)?[adjacent])) as u16
                        };
                    } else {
                        heap.slice_mut(output)?[position + 1] = 0;
                    }
                    heap.slice_mut(output)?[position + 2] = if list_count > 1 {
                        if neighbor_index == 1 {
                            b'(' as u16
                        } else if neighbor_index == list_count {
                            b')' as u16
                        } else {
                            b',' as u16
                        }
                    } else {
                        b'-' as u16
                    };
                } else {
                    if atom < atom_count {
                        heap.slice_mut(neighbor_number)?[atom] = 0;
                    }
                    stack_top -= 1;
                    if stack_top < 0 {
                        break;
                    }
                }
            }
        }
        Ok(output)
    })();

    inchi_free(heap, stack)?;
    inchi_free(heap, number_of_descendants)?;
    inchi_free(heap, dfs_number)?;
    inchi_free(heap, neighbor_number)?;
    FreeNeighList(heap, neighbor_lists)?;
    computation
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichimake__compareicr__line_3189() {
        fn compare(
            flags1: INCHI_MODE,
            flags2: INCHI_MODE,
            mask: INCHI_MODE,
        ) -> (i32, INCHI_MODE, INCHI_MODE) {
            let left = ICR {
                flags: flags1,
                tot_num_H1: i32::MIN,
                ..ICR::default()
            };
            let right = ICR {
                flags: flags2,
                tot_num_H2: i32::MAX,
                ..ICR::default()
            };
            let left_before = left.clone();
            let right_before = right.clone();
            let mut in1 = INCHI_MODE::MAX;
            let mut in2 = INCHI_MODE::MAX;
            let ret = CompareIcr(&left, &right, Some(&mut in1), Some(&mut in2), mask);
            assert_eq!(left, left_before);
            assert_eq!(right, right_before);
            (ret, in1, in2)
        }

        assert_eq!(compare(0, 0, 0), (0, 0, 0));
        for bit_index in 0..=30 {
            let bit = 1_u64 << bit_index;
            assert_eq!(compare(bit, 0, bit), (1, bit, 0), "left bit {bit_index}");
            assert_eq!(compare(0, bit, bit), (-1, 0, bit), "right bit {bit_index}");
            assert_eq!(compare(bit, bit, bit), (0, 0, 0), "shared bit {bit_index}");
            assert_eq!(compare(bit, 0, 0), (0, 0, 0), "masked bit {bit_index}");
        }

        let low = 1_u64;
        let middle = 1_u64 << 17;
        let high = 1_u64 << 30;
        assert_eq!(
            compare(low | middle, high, INCHI_MODE::MAX),
            (2, low | middle, high)
        );
        assert_eq!(
            compare(low | middle, middle | high, INCHI_MODE::MAX),
            (2, low, high)
        );
        assert_eq!(
            compare(low | middle, high, low | middle),
            (1, low | middle, 0)
        );
        assert_eq!(compare(low | middle, high, high), (-1, 0, high));

        let left = ICR {
            flags: low,
            ..ICR::default()
        };
        let right = ICR {
            flags: high,
            ..ICR::default()
        };
        assert_eq!(CompareIcr(&left, &right, None, None, low | high), 2);

        let outside_active_masks = 1_u64 << 63;
        assert_eq!(compare(outside_active_masks, 0, 0), (0, 0, 0));
    }

    #[test]
    fn source_port__ichimake__create_inchi__line_3707() {
        let mut heap = SourceHeap::default();
        let clock = heap
            .allocate_model_storage(vec![INCHI_CLOCK::default()])
            .unwrap();
        let mut empty_outputs = std::array::from_fn(|_| INP_ATOM_DATA::default());
        assert_eq!(
            Create_INChI(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                clock,
                &INPUT_PARMS::default(),
                [SourceMutPointer::null(); TAUT_NUM as usize],
                [SourceMutPointer::null(); TAUT_NUM as usize],
                None,
                SourceMutPointer::null(),
                &mut empty_outputs,
                0,
                u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
                &mut 0,
                &mut 0,
                SourceMutPointer::null(),
                None,
                None,
                0,
            ),
            Ok(-1)
        );

        let mut carbon = inp_ATOM {
            el_number: 6,
            num_H: 4,
            orig_at_number: 1,
            ..inp_ATOM::default()
        };
        carbon.elname[0] = b'C' as i8;
        let input = heap.allocate_model_storage(vec![carbon]).unwrap();
        let normalized = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let inchi_atoms = heap.allocate_model_storage(vec![0_u8]).unwrap();
        let inchi_connections = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let inchi_hydrogens = heap.allocate_model_storage(vec![0_i8]).unwrap();
        let inchi_fixed_hydrogens = heap.allocate_model_storage(vec![0_i8]).unwrap();
        let inchi = heap
            .allocate_model_storage(vec![INChI {
                nAtom: inchi_atoms,
                nConnTable: inchi_connections,
                nNum_H: inchi_hydrogens,
                nNum_H_fixed: inchi_fixed_hydrogens,
                ..INChI::default()
            }])
            .unwrap();
        let aux_order = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let aux_equivalence = heap.allocate_model_storage(vec![0_u16]).unwrap();
        let aux = heap
            .allocate_model_storage(vec![INChI_Aux {
                nOrigAtNosInCanonOrd: aux_order,
                nConstitEquNumbers: aux_equivalence,
                ..INChI_Aux::default()
            }])
            .unwrap();
        let mut outputs = std::array::from_fn(|_| INP_ATOM_DATA::default());
        outputs[TAUT_YES as usize] = INP_ATOM_DATA {
            at: normalized,
            num_at: 1,
            ..INP_ATOM_DATA::default()
        };
        let mut flags = u64::from(TG_FLAG_FIX_ISO_FIXEDH_BUG | TG_FLAG_FIX_TERM_H_CHRG_BUG);
        let mut flags_done = 0_u64;
        let mut globals = CANON_GLOBALS::default();

        assert_eq!(
            Create_INChI(
                &mut heap,
                &mut globals,
                clock,
                &INPUT_PARMS::default(),
                [SourceMutPointer::null(), inchi],
                [SourceMutPointer::null(), aux],
                None,
                input,
                &mut outputs,
                1,
                u64::from(REQ_MODE_TAUT | REQ_MODE_NON_ISO),
                &mut flags,
                &mut flags_done,
                SourceMutPointer::null(),
                None,
                None,
                0,
            ),
            Ok(1)
        );
        let generated = &heap.slice(inchi.as_const()).unwrap()[0];
        assert_eq!(generated.nErrorCode, 0);
        assert_eq!(generated.nNumberOfAtoms, 1);
        assert_eq!(heap.slice(generated.nAtom.as_const()).unwrap(), &[6]);
        assert_eq!(heap.slice(generated.nNum_H.as_const()).unwrap(), &[4]);
        assert_eq!(
            &heap.slice(generated.szHillFormula.as_const()).unwrap()[..4],
            &[b'C' as i8, b'H' as i8, b'4' as i8, 0]
        );
        assert_eq!(outputs[TAUT_YES as usize].num_at, 1);
        assert_eq!(outputs[TAUT_YES as usize].num_removed_H, 0);
        assert_eq!(heap.slice(aux.as_const()).unwrap()[0].nNumberOfAtoms, 1);
    }

    #[test]
    fn source_port__ichimake__checkcanonnumberingcorrectness__line_6230() {
        fn fixture(
            with_isotope: bool,
        ) -> (
            SourceHeap,
            SourceMutPointer<sp_ATOM>,
            CANON_STAT,
            CANON_GLOBALS,
        ) {
            let mut heap = SourceHeap::default();
            let atoms = heap
                .allocate_model_storage(vec![sp_ATOM::default()])
                .unwrap();
            let order = heap.allocate_model_storage(vec![0_u16]).unwrap();
            let isotope = if with_isotope {
                heap.allocate_model_storage(vec![0_u16]).unwrap()
            } else {
                SourceMutPointer::null()
            };
            let linear = heap.allocate_model_storage(vec![1_u16]).unwrap();
            let canon = CANON_STAT {
                LinearCT: linear,
                nMaxLenLinearCT: 1,
                nLenLinearCT: 1,
                nLenLinearCTAtOnly: 1,
                nCanonOrd: order,
                nLenCanonOrd: 1,
                nCanonOrdIsotopic: isotope,
                nLenCanonOrdIsotopic: i32::from(with_isotope),
                ..CANON_STAT::default()
            };
            (heap, atoms, canon, CANON_GLOBALS::default())
        }

        let (mut heap, atoms, mut canon, mut globals) = fixture(true);
        let live = heap.live_allocation_count();
        assert_eq!(
            CheckCanonNumberingCorrectness(
                &mut heap,
                1,
                1,
                atoms.as_const(),
                &mut canon,
                &mut globals,
                0,
                None,
            ),
            Ok(0)
        );
        assert_eq!(heap.live_allocation_count(), live);

        let (mut heap, atoms, mut canon, mut globals) = fixture(false);
        assert_eq!(
            CheckCanonNumberingCorrectness(
                &mut heap,
                1,
                1,
                atoms.as_const(),
                &mut canon,
                &mut globals,
                0,
                None,
            ),
            Ok(0)
        );

        let (mut heap, atoms, mut missing, mut globals) = fixture(false);
        missing.nLenCanonOrd = 0;
        missing.nCanonOrd = SourceMutPointer::null();
        assert_eq!(
            CheckCanonNumberingCorrectness(
                &mut heap,
                1,
                1,
                atoms.as_const(),
                &mut missing,
                &mut globals,
                0,
                None,
            ),
            Ok(CT_CANON_ERR)
        );

        let (mut heap, atoms, mut mismatch, mut globals) = fixture(false);
        heap.slice_mut(mismatch.LinearCT).unwrap()[0] = 0;
        assert_eq!(
            CheckCanonNumberingCorrectness(
                &mut heap,
                1,
                1,
                atoms.as_const(),
                &mut mismatch,
                &mut globals,
                0,
                None,
            ),
            Ok(CT_CANON_ERR)
        );
        assert_eq!(heap.slice(mismatch.LinearCT.as_const()).unwrap(), &[0]);

        let (mut heap, atoms, mut allocation_failure, mut globals) = fixture(false);
        let live = heap.live_allocation_count();
        heap.fail_after_allocations(0);
        assert_eq!(
            CheckCanonNumberingCorrectness(
                &mut heap,
                1,
                1,
                atoms.as_const(),
                &mut allocation_failure,
                &mut globals,
                0,
                None,
            ),
            Ok(CT_CANON_ERR)
        );
        assert_eq!(heap.live_allocation_count(), live);

        let (mut heap, _atoms, mut invalid_access, mut globals) = fixture(false);
        let live = heap.live_allocation_count();
        assert_eq!(
            CheckCanonNumberingCorrectness(
                &mut heap,
                1,
                1,
                SourceConstPointer::null(),
                &mut invalid_access,
                &mut globals,
                0,
                None,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(heap.live_allocation_count(), live);
    }

    #[test]
    fn source_port__ichimake__copy_inp_atom_prefix__source_memcpy_paths()
    -> Result<(), SourceHeapError> {
        fn atom(number: AT_NUMB) -> inp_ATOM {
            inp_ATOM {
                orig_at_number: number,
                ..inp_ATOM::default()
            }
        }

        let mut heap = SourceHeap::default();
        let source = heap
            .allocate_model_storage(vec![atom(1), atom(2), atom(3)])
            .unwrap();
        let destination = heap
            .allocate_model_storage(vec![atom(8), atom(8), atom(8)])
            .unwrap();

        copy_inp_atom_prefix(&mut heap, destination, source.as_const(), 2).unwrap();
        assert_eq!(
            heap.slice(destination.as_const())?
                .iter()
                .map(|atom| atom.orig_at_number)
                .collect::<Vec<_>>(),
            vec![1, 2, 8]
        );

        copy_inp_atom_prefix(&mut heap, source, source.as_const(), 3).unwrap();
        copy_inp_atom_prefix(&mut heap, source.offset(1).unwrap(), source.as_const(), 2).unwrap();
        assert_eq!(
            heap.slice(source.as_const())?
                .iter()
                .map(|atom| atom.orig_at_number)
                .collect::<Vec<_>>(),
            vec![1, 1, 2]
        );

        Ok::<(), SourceHeapError>(())
    }

    #[test]
    fn source_port__ichimake__inp2spatom__line_119() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            inp2spATOM(
                &mut heap,
                SourceConstPointer::null(),
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );

        let mut carbon = inp_ATOM::default();
        carbon.elname = [b'C' as i8, 0, 91, 92, 93, 94];
        carbon.valence = 2;
        carbon.neighbor[..3].copy_from_slice(&[7, 8, 99]);
        carbon.bond_type[..3].copy_from_slice(&[1, 2, 99]);
        carbon.chem_bonds_valence = 3;
        carbon.orig_at_number = 101;
        carbon.orig_compt_at_numb = 202;
        carbon.endpoint = 303;
        carbon.iso_atw_diff = -4;
        carbon.num_H = 5;
        carbon.cFlags = -6;
        carbon.num_iso_H = [7, 8, 9];
        carbon.charge = -2;
        carbon.radical = 3;
        carbon.nBlockSystem = 404;
        carbon.bCutVertex = -1;
        carbon.nRingSystem = 505;
        carbon.nNumAtInRingSystem = 606;
        carbon.bAmbiguousStereo = 44;
        carbon.x = 1.25;
        carbon.y = -2.5;
        carbon.z = 3.75;

        let mut truncated = inp_ATOM::default();
        truncated.elname = [
            b'A' as i8, b'b' as i8, b'c' as i8, b'd' as i8, b'e' as i8, b'f' as i8,
        ];
        truncated.valence = -1;
        truncated.neighbor[0] = 88;
        truncated.bond_type[0] = 77;

        let input = heap
            .allocate_model_storage(vec![carbon.clone(), truncated])
            .unwrap();
        let mut sentinel = sp_ATOM::default();
        sentinel.elname = [9; 6];
        sentinel.orig_at_number = 999;
        let output = heap
            .allocate_model_storage(vec![sentinel.clone(), sentinel.clone(), sentinel.clone()])
            .unwrap();
        assert_eq!(inp2spATOM(&mut heap, input.as_const(), 2, output), Ok(0));
        let converted = heap.slice(output.as_const()).unwrap();
        assert_eq!(converted[0].elname, [b'C' as i8, 0, 0, 0, 0, 0]);
        assert_eq!(converted[0].el_number, 6);
        assert_eq!(converted[0].valence, 2);
        assert_eq!(&converted[0].neighbor[..3], &[7, 8, 0]);
        assert_eq!(&converted[0].bond_type[..3], &[1, 2, 0]);
        assert_eq!(converted[0].chem_bonds_valence, 3);
        assert_eq!(converted[0].orig_at_number, 101);
        assert_eq!(converted[0].orig_compt_at_numb, 202);
        assert_eq!(converted[0].endpoint, 303);
        assert_eq!(converted[0].iso_atw_diff, -4);
        assert_eq!(converted[0].num_H, 5);
        assert_eq!(converted[0].cFlags, -6);
        assert_eq!(converted[0].num_iso_H, [7, 8, 9]);
        assert_eq!(converted[0].charge, -2);
        assert_eq!(converted[0].radical, 3);
        assert_eq!(converted[0].nBlockSystem, 404);
        assert_eq!(converted[0].bCutVertex, -1);
        assert_eq!(converted[0].nRingSystem, 505);
        assert_eq!(converted[0].nNumAtInRingSystem, 606);
        assert_eq!(converted[0].bAmbiguousStereo, 0);
        assert_eq!(
            converted[1].elname,
            [
                b'A' as i8, b'b' as i8, b'c' as i8, b'd' as i8, b'e' as i8, 0
            ]
        );
        assert_eq!(converted[1].el_number, u8::MAX);
        assert_eq!(converted[1].valence, -1);
        assert_eq!(converted[1].neighbor[0], 0);
        assert_eq!(converted[1].bond_type[0], 0);
        assert_eq!(converted[2], sentinel);

        let mut maximum = inp_ATOM::default();
        maximum.elname = [b'O' as i8, 0, 0, 0, 0, 0];
        maximum.valence = 20;
        for index in 0..20 {
            maximum.neighbor[index] = index as u16 + 100;
            maximum.bond_type[index] = index as u8 + 1;
        }
        let maximum = heap.allocate_model_storage(vec![maximum]).unwrap();
        let maximum_output = heap
            .allocate_model_storage(vec![sp_ATOM::default()])
            .unwrap();
        assert_eq!(
            inp2spATOM(&mut heap, maximum.as_const(), 1, maximum_output),
            Ok(0)
        );
        assert_eq!(
            heap.slice(maximum_output.as_const()).unwrap()[0].neighbor,
            std::array::from_fn(|index| index as u16 + 100)
        );

        let mut invalid = carbon;
        invalid.valence = 21;
        let invalid = heap.allocate_model_storage(vec![invalid]).unwrap();
        let invalid_output = heap.allocate_model_storage(vec![sentinel.clone()]).unwrap();
        assert_eq!(
            inp2spATOM(&mut heap, invalid.as_const(), 1, invalid_output),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            heap.slice(invalid_output.as_const()).unwrap()[0].valence,
            21
        );
        assert_eq!(
            inp2spATOM(&mut heap, invalid.as_const(), -1, invalid_output),
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__ichimake__mystrrev__line_2090() {
        for (input, expected) in [
            ("", ""),
            ("a", "a"),
            ("ab", "ba"),
            ("abc", "cba"),
            ("abcdef", "fedcba"),
        ] {
            let mut heap = SourceHeap::default();
            let mut bytes: Vec<i8> = input.bytes().map(|byte| byte as i8).collect();
            bytes.push(0);
            bytes.push(79);
            let string = heap.allocate_model_storage(bytes).unwrap();
            assert_eq!(mystrrev(&mut heap, string), Ok(()));
            let bytes = heap.slice(string.as_const()).unwrap();
            assert_eq!(
                &bytes[..expected.len()],
                &expected.bytes().map(|byte| byte as i8).collect::<Vec<_>>()
            );
            assert_eq!(bytes[expected.len()], 0);
            assert_eq!(bytes[expected.len() + 1], 79);
        }

        let mut heap = SourceHeap::default();
        let unterminated = heap.allocate_model_storage(vec![b'a' as i8]).unwrap();
        assert_eq!(
            mystrrev(&mut heap, unterminated),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(
            mystrrev(&mut heap, SourceMutPointer::null()),
            Err(SourceHeapError::NullPointer)
        );
    }
    use crate::source_types::{
        INChI_IsotopicAtom, INChI_IsotopicTGroup, tagMarkDiff_DIFV_NEQ2PRECED,
    };

    fn inchi_fixture(heap: &mut SourceHeap) -> INChI {
        INChI {
            nNumberOfAtoms: 2,
            szHillFormula: heap.allocate(vec![b'C' as i8, b'2' as i8, 0]).unwrap(),
            nAtom: heap.allocate(vec![6_u8, 6]).unwrap(),
            lenConnTable: 2,
            nConnTable: heap.allocate(vec![1_u16, 2]).unwrap(),
            lenTautomer: 2,
            nTautomer: heap.allocate(vec![2_u16, 1]).unwrap(),
            nNum_H: heap.allocate(vec![3_i8, 3]).unwrap(),
            nNumberOfIsotopicAtoms: 1,
            IsotopicAtom: heap
                .allocate(vec![INChI_IsotopicAtom {
                    nAtomNumber: 1,
                    ..INChI_IsotopicAtom::default()
                }])
                .unwrap(),
            ..INChI::default()
        }
    }

    fn reversed_inchi2_fixture(heap: &mut SourceHeap) -> INChI {
        INChI {
            nNumberOfAtoms: 2,
            szHillFormula: heap
                .allocate(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'6' as i8, 0])
                .unwrap(),
            nAtom: heap.allocate(vec![6_u8, 6]).unwrap(),
            lenConnTable: 2,
            nConnTable: heap.allocate(vec![1_u16, 2]).unwrap(),
            nNum_H: heap.allocate(vec![3_i8, 3]).unwrap(),
            ..INChI::default()
        }
    }

    #[test]
    fn source_port__ichimake__markunusedandemptylayers__line_1525() {
        let main = [9_i8; 11];
        let mut omitted = [main, [1_i8; 11], [1_i8; 11], [1_i8; 11]];
        assert_eq!(MarkUnusedAndEmptyLayers(&mut omitted), 0);
        assert_eq!(omitted[0], main);
        assert_eq!(omitted[1], [0; 11]);
        assert_eq!(omitted[2], [0; 11]);
        assert_eq!(omitted[3], [0; 11]);

        let mut exposed = [[0_i8; 11]; 4];
        exposed[0] = main;
        exposed[tagDiffINChILayers_DIFL_FI as usize][1] = tagMarkDiff_DIFV_NEQ2PRECED as i8;
        exposed[tagDiffINChILayers_DIFL_MI as usize][2] = tagMarkDiff_DIFV_IS_EMPTY as i8;
        assert_eq!(MarkUnusedAndEmptyLayers(&mut exposed), 0);
        assert_eq!(exposed[0], main);
        assert_eq!(
            exposed[tagDiffINChILayers_DIFL_FI as usize]
                [tagDiffINChISegments_DIFS_i_IATOMS as usize],
            tagMarkDiff_DIFV_IS_EMPTY as i8
        );
        assert_eq!(
            exposed[tagDiffINChILayers_DIFL_MI as usize]
                [tagDiffINChISegments_DIFS_i_IATOMS as usize],
            tagMarkDiff_DIFV_IS_EMPTY as i8
        );
        assert_eq!(
            exposed[tagDiffINChILayers_DIFL_F as usize]
                [tagDiffINChISegments_DIFS_f_FORMULA as usize],
            tagMarkDiff_DIFV_IS_EMPTY as i8
        );

        let mut preserved = [[0_i8; 11]; 4];
        preserved[0] = main;
        preserved[tagDiffINChILayers_DIFL_FI as usize]
            [tagDiffINChISegments_DIFS_i_IATOMS as usize] = tagMarkDiff_DIFV_OUTPUT_OMIT_F as i8;
        preserved[tagDiffINChILayers_DIFL_MI as usize]
            [tagDiffINChISegments_DIFS_i_IATOMS as usize] = tagMarkDiff_DIFV_OUTPUT_OMIT_F as i8;
        preserved[tagDiffINChILayers_DIFL_F as usize]
            [tagDiffINChISegments_DIFS_f_FORMULA as usize] = tagMarkDiff_DIFV_OUTPUT_OMIT_F as i8;
        let before = preserved;
        assert_eq!(MarkUnusedAndEmptyLayers(&mut preserved), 0);
        assert_eq!(preserved, before);
    }

    #[test]
    fn source_port__ichimake__compareinchistereo__line_1607() {
        let mut heap = SourceHeap::default();
        assert_eq!(CompareInchiStereo(&heap, None, 0, None, 0), Ok(0));
        let empty = INChI_Stereo::default();
        assert_eq!(CompareInchiStereo(&heap, Some(&empty), 0, None, 0), Ok(0));
        let nonempty = stereo_fixture(&mut heap, &[2], &[1], 1, &[3], &[4], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, None, 0, Some(&nonempty), 0),
            Ok(1)
        );
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, None, 0),
            Ok(-1)
        );

        let equal = stereo_fixture(&mut heap, &[2], &[1], 1, &[3], &[4], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&equal), 0),
            Ok(0)
        );
        let mut changed = stereo_fixture(&mut heap, &[2], &[1], 1, &[5], &[4], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(2)
        );
        changed = stereo_fixture(&mut heap, &[2], &[1], 1, &[3], &[6], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(2)
        );
        changed = stereo_fixture(&mut heap, &[2], &[1], 1, &[3], &[4], &[3]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(2)
        );
        changed = equal.clone();
        changed.nNumberOfStereoBonds = 2;
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(1)
        );

        changed = stereo_fixture(&mut heap, &[5], &[1], 1, &[3], &[4], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(3)
        );
        changed = stereo_fixture(&mut heap, &[2], &[3], 1, &[3], &[4], &[1]);
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(2)
        );
        changed = equal.clone();
        changed.nNumberOfStereoCenters = 2;
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(1)
        );

        changed = equal.clone();
        changed.nCompInv2Abs = -1;
        assert_eq!(
            CompareInchiStereo(&heap, Some(&nonempty), 0, Some(&changed), 0),
            Ok(1)
        );
        assert_eq!(
            CompareInchiStereo(
                &heap,
                Some(&nonempty),
                INCHI_MODE::from(INCHI_FLAG_REL_STEREO),
                Some(&changed),
                0,
            ),
            Ok(0)
        );
        assert_eq!(
            CompareInchiStereo(
                &heap,
                Some(&nonempty),
                0,
                Some(&changed),
                INCHI_MODE::from(INCHI_FLAG_RAC_STEREO),
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichimake__getelementandcount__line_173() {
        let mut heap = SourceHeap::default();
        let formula = heap
            .allocate(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'6' as i8, 0])
            .unwrap();
        let mut output = [99_i8; 6];
        let mut cursor = formula.as_const();
        let mut count = -7;
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(1)
        );
        assert_eq!(count, 2);
        assert_eq!(cursor, formula.as_const().offset(2).unwrap());
        assert_eq!(output, [65, 0, 99, 99, 99, 99]);
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(1)
        );
        assert_eq!(count, 6);
        assert_eq!(cursor, formula.as_const().offset(4).unwrap());
        assert_eq!(output, [72, 0, 99, 99, 99, 99]);
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(0)
        );
        assert_eq!(count, 99999);
        assert_eq!(cursor, formula.as_const().offset(4).unwrap());
        assert_eq!(output, [90, 122, 122, 0, 99, 99]);

        let chlorine = heap
            .allocate(vec![
                b'C' as i8, b'l' as i8, b'1' as i8, b'2' as i8, b'X' as i8, 0,
            ])
            .unwrap();
        cursor = chlorine.as_const();
        count = 0;
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(1)
        );
        assert_eq!(count, 12);
        assert_eq!(cursor, chlorine.as_const().offset(4).unwrap());
        assert_eq!(&output[..4], &[67, 108, 0, 0]);

        let nitrogen = heap.allocate(vec![b'N' as i8, 0]).unwrap();
        cursor = nitrogen.as_const();
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(1)
        );
        assert_eq!(count, 1);
        assert_eq!(cursor, nitrogen.as_const().offset(1).unwrap());

        let invalid = heap.allocate(vec![b'c' as i8, b'2' as i8, 0]).unwrap();
        cursor = invalid.as_const();
        count = 77;
        output.fill(88);
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(-1)
        );
        assert_eq!(cursor, invalid.as_const());
        assert_eq!(count, 77);
        assert_eq!(output, [88; 6]);

        cursor = SourceConstPointer::null();
        count = 0;
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(0)
        );
        assert!(cursor.is_null());
        assert_eq!(count, 99999);
        assert_eq!(&output[..5], &[90, 122, 122, 0, 88]);

        let overflow = heap
            .allocate(
                b"N9223372036854775808X\0"
                    .iter()
                    .map(|byte| *byte as i8)
                    .collect(),
            )
            .unwrap();
        cursor = overflow.as_const();
        assert_eq!(
            GetElementAndCount(&mut heap, &mut cursor, &mut output, &mut count),
            Ok(1)
        );
        assert_eq!(count, -1);
        assert_eq!(cursor, overflow.as_const().offset(20).unwrap());
    }

    #[test]
    fn source_port__ichimake__comparehillformulas__line_241() {
        let mut heap = SourceHeap::default();
        let mut formula = |bytes: &[u8]| {
            let mut bytes: Vec<i8> = bytes.iter().map(|byte| *byte as i8).collect();
            bytes.push(0);
            heap.allocate(bytes).unwrap().as_const()
        };
        let c2h6 = formula(b"C2H6");
        let c3h6 = formula(b"C3H6");
        let c2n = formula(b"C2N");
        let c2o = formula(b"C2O");
        let oxygen = formula(b"O");
        let empty = formula(b"");
        let invalid = formula(b"c2");
        let max_count = formula(b"N2147483647");
        let wrapped_count = formula(b"N2147483648");

        assert_eq!(CompareHillFormulas(&mut heap, c2h6, c2h6), Ok(0));
        assert_eq!(CompareHillFormulas(&mut heap, c2h6, c3h6), Ok(1));
        assert_eq!(CompareHillFormulas(&mut heap, c3h6, c2h6), Ok(-1));
        assert_eq!(CompareHillFormulas(&mut heap, c2n, c2o), Ok(-1));
        assert_eq!(CompareHillFormulas(&mut heap, c2o, c2n), Ok(1));
        assert_eq!(CompareHillFormulas(&mut heap, oxygen, empty), Ok(-11));
        assert_eq!(CompareHillFormulas(&mut heap, empty, oxygen), Ok(11));
        assert_eq!(CompareHillFormulas(&mut heap, invalid, oxygen), Ok(0));
        assert_eq!(
            CompareHillFormulas(&mut heap, SourceConstPointer::null(), empty),
            Ok(0)
        );
        assert_eq!(
            CompareHillFormulas(&mut heap, max_count, wrapped_count),
            Ok(1)
        );
    }

    #[test]
    fn source_port__ichimake__compinchi2__line_1712() {
        fn inchi(heap: &mut SourceHeap, formula: &[u8], atoms: &[u8], hydrogens: &[i8]) -> INChI {
            let mut formula_bytes: Vec<i8> = formula.iter().map(|byte| *byte as i8).collect();
            formula_bytes.push(0);
            INChI {
                nNumberOfAtoms: atoms.len() as i32,
                szHillFormula: heap.allocate(formula_bytes).unwrap(),
                nAtom: heap.allocate(atoms.to_vec()).unwrap(),
                nNum_H: heap.allocate(hydrogens.to_vec()).unwrap(),
                ..INChI::default()
            }
        }

        fn compare(
            heap: &mut SourceHeap,
            left: INChI,
            right: INChI,
            taut: u32,
            isotopic: i32,
        ) -> i32 {
            let mut first = INCHI_SORT::default();
            let mut second = INCHI_SORT::default();
            first.pINChI[TAUT_NON as usize] = heap.allocate(vec![left]).unwrap();
            second.pINChI[TAUT_NON as usize] = heap.allocate(vec![right]).unwrap();
            CompINChI2(heap, &first, &second, taut, isotopic).unwrap()
        }

        let mut heap = SourceHeap::default();
        let empty1 = INCHI_SORT::default();
        let empty2 = INCHI_SORT::default();
        assert_eq!(CompINChI2(&mut heap, &empty1, &empty2, TAUT_NON, 1), Ok(0));

        let live = inchi(&mut heap, b"C", &[6], &[0]);
        let mut only_left = INCHI_SORT::default();
        only_left.pINChI[TAUT_NON as usize] = heap.allocate(vec![live.clone()]).unwrap();
        assert_eq!(
            CompINChI2(&mut heap, &only_left, &empty2, TAUT_NON, 1),
            Ok(-1)
        );
        assert_eq!(
            CompINChI2(&mut heap, &empty1, &only_left, TAUT_NON, 1),
            Ok(1)
        );

        let mut deleted = live.clone();
        deleted.bDeleted = 1;
        assert_eq!(
            compare(&mut heap, deleted.clone(), live.clone(), TAUT_NON, 1),
            1
        );
        assert_eq!(compare(&mut heap, live.clone(), deleted, TAUT_NON, 1), -1);

        let nitrogen = inchi(&mut heap, b"N", &[7], &[0]);
        assert!(compare(&mut heap, live.clone(), nitrogen, TAUT_NON, 1) < 0);
        let c2 = inchi(&mut heap, b"C2", &[6, 6], &[0, 0]);
        assert_eq!(compare(&mut heap, live.clone(), c2, TAUT_NON, 1), 1);

        let mut atom2 = live.clone();
        atom2.nAtom = heap.allocate(vec![7]).unwrap();
        assert_eq!(compare(&mut heap, live.clone(), atom2, TAUT_NON, 1), 1);
        let mut conn1 = live.clone();
        conn1.lenConnTable = 1;
        conn1.nConnTable = heap.allocate(vec![1_u16]).unwrap();
        let mut conn2 = live.clone();
        conn2.lenConnTable = 1;
        conn2.nConnTable = heap.allocate(vec![2_u16]).unwrap();
        assert_eq!(
            compare(&mut heap, conn1.clone(), conn2.clone(), TAUT_NON, 1),
            1
        );
        conn2.lenConnTable = 2;
        conn2.nConnTable = heap.allocate(vec![1_u16, 2]).unwrap();
        assert_eq!(compare(&mut heap, conn1, conn2, TAUT_NON, 1), 1);

        let ch2 = inchi(&mut heap, b"CH2", &[6], &[0]);
        let ch3 = inchi(&mut heap, b"CH3", &[6], &[0]);
        assert_eq!(compare(&mut heap, ch2, ch3, TAUT_NON, 1), 1);
        let h1 = inchi(&mut heap, b"C", &[6], &[1]);
        let h2 = inchi(&mut heap, b"C", &[6], &[2]);
        assert_eq!(compare(&mut heap, h1, h2, TAUT_NON, 1), 1);

        let mut taut1 = live.clone();
        taut1.lenTautomer = 1;
        taut1.nTautomer = heap.allocate(vec![1_u16]).unwrap();
        let mut taut2 = live.clone();
        taut2.lenTautomer = 1;
        taut2.nTautomer = heap.allocate(vec![2_u16]).unwrap();
        assert_eq!(compare(&mut heap, taut1, taut2, TAUT_YES, 1), 1);

        let mut fixed1 = live.clone();
        fixed1.nNum_H_fixed = heap.allocate(vec![1_i8]).unwrap();
        let mut fixed2 = live.clone();
        fixed2.nNum_H_fixed = heap.allocate(vec![2_i8]).unwrap();
        let mut fixed_sort1 = INCHI_SORT::default();
        fixed_sort1.pINChI[TAUT_YES as usize] = heap.allocate(vec![live.clone()]).unwrap();
        fixed_sort1.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed1]).unwrap();
        let mut fixed_sort2 = INCHI_SORT::default();
        fixed_sort2.pINChI[TAUT_YES as usize] = heap.allocate(vec![live.clone()]).unwrap();
        fixed_sort2.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed2]).unwrap();
        assert_eq!(
            CompINChI2(&mut heap, &fixed_sort1, &fixed_sort2, TAUT_NON, 1),
            Ok(1)
        );

        let mut stereo_inchi = live.clone();
        let stereo = INChI_Stereo {
            nNumberOfStereoBonds: 1,
            nBondAtom1: heap.allocate(vec![1_u16]).unwrap(),
            nBondAtom2: heap.allocate(vec![2_u16]).unwrap(),
            b_parity: heap.allocate(vec![1_i8]).unwrap(),
            ..INChI_Stereo::default()
        };
        stereo_inchi.Stereo = heap.allocate(vec![stereo]).unwrap();
        assert_eq!(
            compare(&mut heap, live.clone(), stereo_inchi, TAUT_NON, 1),
            1
        );

        let mut isotope1 = live.clone();
        isotope1.nNumberOfIsotopicAtoms = 1;
        isotope1.IsotopicAtom = heap
            .allocate(vec![crate::source_types::INChI_IsotopicAtom {
                nAtomNumber: 1,
                nIsoDifference: 1,
                nNum_H: 1,
                nNum_D: 1,
                nNum_T: 1,
            }])
            .unwrap();
        let mut isotope2 = isotope1.clone();
        isotope2.IsotopicAtom = heap
            .allocate(vec![crate::source_types::INChI_IsotopicAtom {
                nAtomNumber: 1,
                nIsoDifference: 1,
                nNum_H: 1,
                nNum_D: 1,
                nNum_T: 2,
            }])
            .unwrap();
        assert_eq!(
            compare(&mut heap, isotope1.clone(), isotope2.clone(), TAUT_NON, 0),
            0
        );
        assert_eq!(compare(&mut heap, isotope1, isotope2, TAUT_NON, 1), 1);

        let mut group1 = live.clone();
        group1.nNumberOfIsotopicTGroups = 1;
        group1.IsotopicTGroup = heap
            .allocate(vec![crate::source_types::INChI_IsotopicTGroup {
                nTGroupNumber: 1,
                nNum_H: 1,
                nNum_D: 1,
                nNum_T: 1,
            }])
            .unwrap();
        let mut group2 = group1.clone();
        group2.IsotopicTGroup = heap
            .allocate(vec![crate::source_types::INChI_IsotopicTGroup {
                nTGroupNumber: 2,
                nNum_H: 1,
                nNum_D: 1,
                nNum_T: 1,
            }])
            .unwrap();
        assert_eq!(compare(&mut heap, group1, group2, TAUT_NON, 1), 1);

        let mut charged1 = live.clone();
        charged1.nTotalCharge = -1;
        let mut charged2 = live.clone();
        charged2.nTotalCharge = 1;
        assert_eq!(compare(&mut heap, charged1, charged2, TAUT_NON, 1), -2);
        let mut charged = live.clone();
        charged.nTotalCharge = 1;
        assert_eq!(compare(&mut heap, charged, live, TAUT_NON, 1), 1);
    }

    #[test]
    fn source_port__ichimake__compinchinontaut2__line_2042() {
        fn inchi(heap: &mut SourceHeap, charge: i32) -> INChI {
            INChI {
                nTotalCharge: charge,
                nNumberOfAtoms: 1,
                szHillFormula: heap.allocate(vec![b'C' as i8, 0]).unwrap(),
                nAtom: heap.allocate(vec![6_u8]).unwrap(),
                nNum_H: heap.allocate(vec![0_i8]).unwrap(),
                ..INChI::default()
            }
        }

        let mut heap = SourceHeap::default();
        let mut first = INCHI_SORT::default();
        let mut second = INCHI_SORT::default();
        let charged_left = inchi(&mut heap, -1);
        let charged_right = inchi(&mut heap, 1);
        first.pINChI[TAUT_NON as usize] = heap.allocate(vec![charged_left]).unwrap();
        second.pINChI[TAUT_NON as usize] = heap.allocate(vec![charged_right]).unwrap();
        assert_eq!(CompINChINonTaut2(&mut heap, &first, &second), Ok(-2));

        let mut stable1 = INCHI_SORT {
            ord_number: i16::MAX,
            ..INCHI_SORT::default()
        };
        let mut stable2 = INCHI_SORT {
            ord_number: i16::MIN,
            ..INCHI_SORT::default()
        };
        let equal1 = inchi(&mut heap, 0);
        let equal2 = inchi(&mut heap, 0);
        stable1.pINChI[TAUT_NON as usize] = heap.allocate(vec![equal1]).unwrap();
        stable2.pINChI[TAUT_NON as usize] = heap.allocate(vec![equal2]).unwrap();
        assert_eq!(CompINChINonTaut2(&mut heap, &stable1, &stable2), Ok(65535));

        let mut mobile_stereo = inchi(&mut heap, 0);
        let bond_atom1 = heap.allocate(vec![1_u16]).unwrap();
        let bond_atom2 = heap.allocate(vec![2_u16]).unwrap();
        let bond_parity = heap.allocate(vec![1_i8]).unwrap();
        mobile_stereo.Stereo = heap
            .allocate(vec![INChI_Stereo {
                nNumberOfStereoBonds: 1,
                nBondAtom1: bond_atom1,
                nBondAtom2: bond_atom2,
                b_parity: bond_parity,
                ..INChI_Stereo::default()
            }])
            .unwrap();
        let mut fallback1 = INCHI_SORT::default();
        let mut fallback2 = INCHI_SORT::default();
        let mobile_plain = inchi(&mut heap, 0);
        let fixed_plain1 = inchi(&mut heap, 0);
        let fixed_plain2 = inchi(&mut heap, 0);
        fallback1.pINChI[TAUT_YES as usize] = heap.allocate(vec![mobile_plain]).unwrap();
        fallback2.pINChI[TAUT_YES as usize] = heap.allocate(vec![mobile_stereo]).unwrap();
        fallback1.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed_plain1]).unwrap();
        fallback2.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed_plain2]).unwrap();
        assert_eq!(CompINChINonTaut2(&mut heap, &fallback1, &fallback2), Ok(1));
    }

    #[test]
    fn source_port__ichimake__compinchitaut2__line_2064() {
        fn inchi(heap: &mut SourceHeap, charge: i32) -> INChI {
            INChI {
                nTotalCharge: charge,
                nNumberOfAtoms: 1,
                szHillFormula: heap.allocate(vec![b'C' as i8, 0]).unwrap(),
                nAtom: heap.allocate(vec![6_u8]).unwrap(),
                nNum_H: heap.allocate(vec![0_i8]).unwrap(),
                ..INChI::default()
            }
        }

        let mut heap = SourceHeap::default();
        let mut first = INCHI_SORT::default();
        let mut second = INCHI_SORT::default();
        let charged_left = inchi(&mut heap, -1);
        let charged_right = inchi(&mut heap, 1);
        first.pINChI[TAUT_YES as usize] = heap.allocate(vec![charged_left]).unwrap();
        second.pINChI[TAUT_YES as usize] = heap.allocate(vec![charged_right]).unwrap();
        assert_eq!(CompINChITaut2(&mut heap, &first, &second), Ok(-2));

        let mut stable1 = INCHI_SORT {
            ord_number: i16::MAX,
            ..INCHI_SORT::default()
        };
        let mut stable2 = INCHI_SORT {
            ord_number: i16::MIN,
            ..INCHI_SORT::default()
        };
        let equal1 = inchi(&mut heap, 0);
        let equal2 = inchi(&mut heap, 0);
        stable1.pINChI[TAUT_YES as usize] = heap.allocate(vec![equal1]).unwrap();
        stable2.pINChI[TAUT_YES as usize] = heap.allocate(vec![equal2]).unwrap();
        assert_eq!(CompINChITaut2(&mut heap, &stable1, &stable2), Ok(65535));

        let mut fixed_stereo = inchi(&mut heap, 0);
        let bond_atom1 = heap.allocate(vec![1_u16]).unwrap();
        let bond_atom2 = heap.allocate(vec![2_u16]).unwrap();
        let bond_parity = heap.allocate(vec![1_i8]).unwrap();
        fixed_stereo.Stereo = heap
            .allocate(vec![INChI_Stereo {
                nNumberOfStereoBonds: 1,
                nBondAtom1: bond_atom1,
                nBondAtom2: bond_atom2,
                b_parity: bond_parity,
                ..INChI_Stereo::default()
            }])
            .unwrap();
        let mut fallback1 = INCHI_SORT::default();
        let mut fallback2 = INCHI_SORT::default();
        let mobile_plain1 = inchi(&mut heap, 0);
        let mobile_plain2 = inchi(&mut heap, 0);
        let fixed_plain = inchi(&mut heap, 0);
        fallback1.pINChI[TAUT_YES as usize] = heap.allocate(vec![mobile_plain1]).unwrap();
        fallback2.pINChI[TAUT_YES as usize] = heap.allocate(vec![mobile_plain2]).unwrap();
        fallback1.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed_plain]).unwrap();
        fallback2.pINChI[TAUT_NON as usize] = heap.allocate(vec![fixed_stereo]).unwrap();
        assert_eq!(CompINChITaut2(&mut heap, &fallback1, &fallback2), Ok(1));
    }

    #[test]
    fn source_port__ichimake__comparehillformulasnoh__line_272() {
        let mut heap = SourceHeap::default();
        let mut formula = |bytes: &[u8]| {
            let mut bytes: Vec<i8> = bytes.iter().map(|byte| *byte as i8).collect();
            bytes.push(0);
            heap.allocate(bytes).unwrap().as_const()
        };
        let c2h6 = formula(b"C2H6");
        let c2h4 = formula(b"C2H4");
        let c2n = formula(b"C2N");
        let c3n = formula(b"C3N");
        let c2o = formula(b"C2O");
        let h2o = formula(b"H2O");
        let oxygen = formula(b"O");
        let hydrogen = formula(b"H2");
        let empty = formula(b"");
        let invalid = formula(b"c2");
        let mut h1 = 10;
        let mut h2 = 20;
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, c2h6, c2h4, &mut h1, &mut h2),
            Ok(0)
        );
        assert_eq!((h1, h2), (16, 24));

        h1 = 0;
        h2 = 0;
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, c2n, c3n, &mut h1, &mut h2),
            Ok(1)
        );
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, c2n, c2o, &mut h1, &mut h2),
            Ok(-1)
        );
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, h2o, oxygen, &mut h1, &mut h2),
            Ok(0)
        );
        assert_eq!((h1, h2), (2, 0));
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, hydrogen, empty, &mut h1, &mut h2),
            Ok(0)
        );
        assert_eq!(h1, 4);
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, empty, empty, &mut h1, &mut h2),
            Ok(0)
        );
        assert_eq!(
            CompareHillFormulasNoH(&mut heap, invalid, oxygen, &mut h1, &mut h2),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichimake__comparetautnonisopartofinchi__line_316() {
        let mut heap = SourceHeap::default();
        let empty1 = INChI::default();
        let empty2 = INChI {
            lenTautomer: -1,
            ..INChI::default()
        };
        assert_eq!(CompareTautNonIsoPartOfINChI(&heap, &empty1, &empty2), Ok(0));

        let zero_first1 = heap.allocate(vec![0_u16, 9]).unwrap();
        let zero_first2 = heap.allocate(vec![0_u16, 7]).unwrap();
        let ignored1 = INChI {
            lenTautomer: 2,
            nTautomer: zero_first1,
            ..INChI::default()
        };
        let ignored2 = INChI {
            lenTautomer: 2,
            nTautomer: zero_first2,
            ..INChI::default()
        };
        assert_eq!(
            CompareTautNonIsoPartOfINChI(&heap, &ignored1, &ignored2),
            Ok(0)
        );

        let values1 = heap.allocate(vec![2_u16, 3]).unwrap();
        let values2 = heap.allocate(vec![2_u16, 3]).unwrap();
        let mut left = INChI {
            lenTautomer: 2,
            nTautomer: values1,
            ..INChI::default()
        };
        let mut right = INChI {
            lenTautomer: 2,
            nTautomer: values2,
            ..INChI::default()
        };
        assert_eq!(CompareTautNonIsoPartOfINChI(&heap, &left, &right), Ok(0));
        right.lenTautomer = 3;
        assert_eq!(CompareTautNonIsoPartOfINChI(&heap, &left, &right), Ok(1));
        right.lenTautomer = 2;
        right.nTautomer = heap.allocate(vec![2_u16, 5]).unwrap();
        assert_eq!(CompareTautNonIsoPartOfINChI(&heap, &left, &right), Ok(2));
        left.nTautomer = heap.allocate(vec![4_u16, 3]).unwrap();
        right.nTautomer = heap.allocate(vec![2_u16, 3]).unwrap();
        assert_eq!(CompareTautNonIsoPartOfINChI(&heap, &left, &right), Ok(-2));
    }

    fn comp_taut_fixture(heap: &mut SourceHeap) -> (INChI, INChI) {
        let make = |heap: &mut SourceHeap| INChI {
            nNumberOfAtoms: 2,
            szHillFormula: heap.allocate(vec![b'C' as i8, b'2' as i8, 0]).unwrap(),
            nAtom: heap.allocate(vec![6_u8, 6]).unwrap(),
            lenConnTable: 2,
            nConnTable: heap.allocate(vec![1_u16, 2]).unwrap(),
            lenTautomer: 2,
            nTautomer: heap.allocate(vec![2_u16, 1]).unwrap(),
            nNum_H: heap.allocate(vec![1_i8, 1]).unwrap(),
            ..INChI::default()
        };
        (make(heap), make(heap))
    }

    fn comp_taut_pair(
        heap: &mut SourceHeap,
        i1: INChI,
        i2: INChI,
        compare_isotopic: i32,
    ) -> Result<i32, SourceHeapError> {
        let i1 = heap.allocate(vec![i1]).unwrap();
        let i2 = heap.allocate(vec![i2]).unwrap();
        let mut p1 = INCHI_SORT::default();
        let mut p2 = INCHI_SORT::default();
        p1.pINChI[TAUT_YES as usize] = i1;
        p2.pINChI[TAUT_NON as usize] = i2;
        CompINChITautVsNonTaut(heap, &p1, &p2, compare_isotopic)
    }

    fn assert_comp_taut<F>(expected: i32, compare_isotopic: i32, update: F)
    where
        F: FnOnce(&mut SourceHeap, &mut INChI, &mut INChI),
    {
        let mut heap = SourceHeap::default();
        let (mut i1, mut i2) = comp_taut_fixture(&mut heap);
        update(&mut heap, &mut i1, &mut i2);
        assert_eq!(
            comp_taut_pair(&mut heap, i1, i2, compare_isotopic),
            Ok(expected)
        );
    }

    #[test]
    fn source_port__ichimake__compinchitautvsnontaut__line_341() {
        let mut heap = SourceHeap::default();
        let (i1, i2) = comp_taut_fixture(&mut heap);
        let i1_pointer = heap.allocate(vec![i1.clone()]).unwrap();
        let i2_pointer = heap.allocate(vec![i2.clone()]).unwrap();
        let mut p1 = INCHI_SORT::default();
        let mut p2 = INCHI_SORT::default();
        p1.pINChI[TAUT_NON as usize] = i1_pointer;
        p2.pINChI[TAUT_NON as usize] = i2_pointer;
        assert_eq!(CompINChITautVsNonTaut(&mut heap, &p1, &p2, 1), Ok(0));

        p1 = INCHI_SORT::default();
        p1.pINChI[TAUT_YES as usize] = i1_pointer;
        p2 = INCHI_SORT::default();
        assert_eq!(CompINChITautVsNonTaut(&mut heap, &p1, &p2, 1), Ok(0));
        assert_eq!(
            CompINChITautVsNonTaut(&mut heap, &INCHI_SORT::default(), &INCHI_SORT::default(), 1,),
            Ok(0)
        );

        assert_comp_taut(1, 0, |_, i1, _| i1.bDeleted = 1);
        assert_comp_taut(-1, 0, |_, _, i2| i2.bDeleted = 1);
        assert_comp_taut(1, 0, |heap, i1, i2| {
            i1.szHillFormula = heap
                .allocate(vec![b'C' as i8, b'2' as i8, b'N' as i8, 0])
                .unwrap();
            i2.szHillFormula = heap
                .allocate(vec![b'C' as i8, b'3' as i8, b'N' as i8, 0])
                .unwrap();
        });
        assert_comp_taut(1, 0, |heap, i1, i2| {
            i1.szHillFormula = heap
                .allocate(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'2' as i8, 0])
                .unwrap();
            i2.szHillFormula = heap
                .allocate(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'3' as i8, 0])
                .unwrap();
        });
        assert_comp_taut(1, 0, |_, _, i2| i2.nNumberOfAtoms = 3);
        assert_comp_taut(1, 0, |heap, _, i2| {
            i2.nAtom = heap.allocate(vec![6_u8, 7]).unwrap();
        });
        assert_comp_taut(1, 0, |_, _, i2| i2.lenConnTable = 3);
        assert_comp_taut(2, 0, |heap, _, i2| {
            i2.nConnTable = heap.allocate(vec![1_u16, 4]).unwrap();
        });
        assert_comp_taut(0, 0, |_, i1, i2| {
            i1.lenConnTable = 0;
            i1.nConnTable = Default::default();
            i2.lenConnTable = 0;
            i2.nConnTable = Default::default();
        });

        assert_comp_taut(1, 0, |heap, _, i2| {
            i2.nNum_H = heap.allocate(vec![0_i8, 1]).unwrap();
        });
        assert_comp_taut(-1, 0, |heap, i1, _| {
            i1.nNum_H = heap.allocate(vec![0_i8, 1]).unwrap();
        });
        assert_comp_taut(2, 0, |heap, _, i2| {
            i2.nNum_H = heap.allocate(vec![3_i8, 1]).unwrap();
        });
        assert_comp_taut(1, 0, |_, _, i2| i2.lenTautomer = 3);
        assert_comp_taut(2, 0, |heap, _, i2| {
            i2.nTautomer = heap.allocate(vec![2_u16, 3]).unwrap();
        });
        assert_comp_taut(0, 0, |heap, _, i2| {
            i2.nNum_H_fixed = heap.allocate(vec![0_i8, 0]).unwrap();
        });
        assert_comp_taut(1, 0, |heap, _, i2| {
            i2.nNum_H_fixed = heap.allocate(vec![0_i8, 1]).unwrap();
        });
        assert_comp_taut(1, 0, |heap, _, i2| {
            let stereo = stereo_fixture(heap, &[1], &[1], 0, &[], &[], &[]);
            i2.Stereo = heap.allocate(vec![stereo]).unwrap();
        });

        assert_comp_taut(0, 0, |_, _, i2| i2.nNumberOfIsotopicAtoms = 1);
        assert_comp_taut(1, 1, |_, _, i2| i2.nNumberOfIsotopicAtoms = 1);
        for (field, expected) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)] {
            assert_comp_taut(expected, 1, |heap, i1, i2| {
                let left = INChI_IsotopicAtom {
                    nAtomNumber: 1,
                    ..INChI_IsotopicAtom::default()
                };
                let mut right = left.clone();
                match field {
                    0 => right.nAtomNumber += 1,
                    1 => right.nIsoDifference += 2,
                    2 => right.nNum_T += 3,
                    3 => right.nNum_D += 4,
                    _ => right.nNum_H += 5,
                }
                i1.nNumberOfIsotopicAtoms = 1;
                i2.nNumberOfIsotopicAtoms = 1;
                i1.IsotopicAtom = heap.allocate(vec![left]).unwrap();
                i2.IsotopicAtom = heap.allocate(vec![right]).unwrap();
            });
        }
        assert_comp_taut(1, 1, |_, i1, _| i1.nNumberOfIsotopicTGroups = -7);
        assert_comp_taut(1, 1, |heap, _, i2| {
            let stereo = stereo_fixture(heap, &[1], &[1], 0, &[], &[], &[]);
            i2.StereoIsotopic = heap.allocate(vec![stereo]).unwrap();
        });

        assert_comp_taut(-3, 1, |_, i1, i2| {
            i1.nTotalCharge = -2;
            i2.nTotalCharge = 1;
        });
        assert_comp_taut(1, 1, |_, i1, _| i1.nTotalCharge = 1);
        assert_comp_taut(-1, 1, |_, _, i2| i2.nTotalCharge = 1);
        assert_comp_taut(0, 1, |_, _, _| {});
    }

    #[test]
    fn source_port__ichimake__getsp3relracabs__line_593() {
        const SP3_NONE: i32 = 0;
        const SP3_ONLY: i32 = 1;
        const SP3_ABS: i32 = 2;
        const SP3_REL: i32 = 4;
        const SP3_RAC: i32 = 8;

        let plain = INChI::default();
        let mut stereo = INChI_Stereo::default();
        assert_eq!(GetSp3RelRacAbs(None, None), SP3_NONE);
        assert_eq!(GetSp3RelRacAbs(Some(&plain), None), SP3_NONE);
        assert_eq!(GetSp3RelRacAbs(None, Some(&stereo)), SP3_NONE);
        assert_eq!(GetSp3RelRacAbs(Some(&plain), Some(&stereo)), SP3_NONE);

        stereo.nNumberOfStereoCenters = 1;
        stereo.nCompInv2Abs = 1;
        let deleted = INChI {
            bDeleted: 1,
            ..INChI::default()
        };
        assert_eq!(GetSp3RelRacAbs(Some(&deleted), Some(&stereo)), SP3_NONE);
        assert_eq!(GetSp3RelRacAbs(Some(&plain), Some(&stereo)), SP3_ABS);

        let relative = INChI {
            nFlags: INCHI_MODE::from(INCHI_FLAG_REL_STEREO),
            ..INChI::default()
        };
        let racemic = INChI {
            nFlags: INCHI_MODE::from(INCHI_FLAG_RAC_STEREO),
            ..INChI::default()
        };
        assert_eq!(GetSp3RelRacAbs(Some(&relative), Some(&stereo)), SP3_REL);
        assert_eq!(GetSp3RelRacAbs(Some(&racemic), Some(&stereo)), SP3_RAC);
        stereo.nNumberOfStereoCenters = 2;
        assert_eq!(GetSp3RelRacAbs(Some(&relative), Some(&stereo)), SP3_REL);
        assert_eq!(GetSp3RelRacAbs(Some(&racemic), Some(&stereo)), SP3_RAC);

        stereo.nCompInv2Abs = 0;
        stereo.nNumberOfStereoCenters = 1;
        assert_eq!(GetSp3RelRacAbs(Some(&plain), Some(&stereo)), SP3_ONLY);
        assert_eq!(GetSp3RelRacAbs(Some(&relative), Some(&stereo)), SP3_ONLY);
        assert_eq!(GetSp3RelRacAbs(Some(&racemic), Some(&stereo)), SP3_ONLY);
        stereo.nNumberOfStereoCenters = 2;
        assert_eq!(GetSp3RelRacAbs(Some(&relative), Some(&stereo)), SP3_ONLY);
        assert_eq!(GetSp3RelRacAbs(Some(&racemic), Some(&stereo)), SP3_ONLY);
    }

    fn layers_inchi(heap: &mut SourceHeap) -> INChI {
        let stereo = stereo_fixture(heap, &[1, 2], &[1, 2], 1, &[1], &[2], &[1]);
        let isotopic_stereo = stereo_fixture(heap, &[1, 2], &[1, 2], 1, &[1], &[2], &[1]);
        INChI {
            nTotalCharge: 1,
            nNumberOfAtoms: 2,
            szHillFormula: heap.allocate(vec![b'C' as i8, b'2' as i8, 0]).unwrap(),
            nAtom: heap.allocate(vec![6_u8, 6]).unwrap(),
            lenConnTable: 2,
            nConnTable: heap.allocate(vec![1_u16, 2]).unwrap(),
            lenTautomer: 1,
            nTautomer: heap.allocate(vec![0_u16]).unwrap(),
            nNum_H: heap.allocate(vec![1_i8, 0]).unwrap(),
            nNum_H_fixed: heap.allocate(vec![1_i8, 0]).unwrap(),
            nNumberOfIsotopicAtoms: 1,
            IsotopicAtom: heap
                .allocate(vec![INChI_IsotopicAtom {
                    nAtomNumber: 1,
                    nIsoDifference: 1,
                    nNum_H: 1,
                    nNum_D: 1,
                    nNum_T: 1,
                }])
                .unwrap(),
            nNumberOfIsotopicTGroups: 1,
            IsotopicTGroup: heap
                .allocate(vec![INChI_IsotopicTGroup {
                    nTGroupNumber: 1,
                    nNum_H: 1,
                    nNum_D: 1,
                    nNum_T: 1,
                }])
                .unwrap(),
            Stereo: heap.allocate(vec![stereo]).unwrap(),
            StereoIsotopic: heap.allocate(vec![isotopic_stereo]).unwrap(),
            ..INChI::default()
        }
    }

    fn run_layers(
        heap: &mut SourceHeap,
        first: Option<INChI>,
        second: Option<INChI>,
        fix_charge: i32,
    ) -> (i32, [[i8; 11]; 4]) {
        let mut p1 = INCHI_SORT::default();
        let mut p2 = INCHI_SORT::default();
        if let Some(first) = first {
            p1.pINChI[TAUT_YES as usize] = heap.allocate(vec![first]).unwrap();
        }
        if let Some(second) = second {
            p2.pINChI[TAUT_NON as usize] = heap.allocate(vec![second]).unwrap();
        }
        let mut segments = [[0_i8; 11]; 4];
        let result = CompINChILayers(heap, &p1, &p2, &mut segments, fix_charge).unwrap();
        (result, segments)
    }

    #[test]
    fn source_port__ichimake__compinchilayers__line_645() {
        const EQL: i8 = tagMarkDiff_DIFV_EQL2PRECED as i8;
        const NEQ: i8 = tagMarkDiff_DIFV_NEQ2PRECED as i8;
        const EMPTY: i8 = tagMarkDiff_DIFV_IS_EMPTY as i8;
        const FI_MI: i8 = tagMarkDiff_DIFV_FI_EQ_MI as i8;

        let mut heap = SourceHeap::default();
        assert_eq!(run_layers(&mut heap, None, None, 0), (0, [[0; 11]; 4]));

        let first = layers_inchi(&mut heap);
        let second = layers_inchi(&mut heap);
        let (result, segments) = run_layers(&mut heap, Some(first), Some(second), 0);
        assert_eq!(result, 0);
        let mut expected = [[0_i8; 11]; 4];
        expected[tagDiffINChILayers_DIFL_M as usize]
            [tagDiffINChISegments_DIFS_f_FORMULA as usize] = NEQ;
        expected[tagDiffINChILayers_DIFL_F as usize]
            [tagDiffINChISegments_DIFS_f_FORMULA as usize] = EQL;
        expected[tagDiffINChILayers_DIFL_M as usize]
            [tagDiffINChISegments_DIFS_h_H_ATOMS as usize] = NEQ;
        expected[tagDiffINChILayers_DIFL_F as usize]
            [tagDiffINChISegments_DIFS_h_H_ATOMS as usize] = NEQ;
        expected[tagDiffINChILayers_DIFL_M as usize][tagDiffINChISegments_DIFS_q_CHARGE as usize] =
            NEQ;
        expected[tagDiffINChILayers_DIFL_F as usize][tagDiffINChISegments_DIFS_q_CHARGE as usize] =
            EQL;
        expected[tagDiffINChILayers_DIFL_M as usize][tagDiffINChISegments_DIFS_b_SBONDS as usize] =
            NEQ;
        expected[tagDiffINChILayers_DIFL_F as usize][tagDiffINChISegments_DIFS_b_SBONDS as usize] =
            EQL;
        expected[tagDiffINChILayers_DIFL_MI as usize]
            [tagDiffINChISegments_DIFS_b_SBONDS as usize] = EQL;
        expected[tagDiffINChILayers_DIFL_FI as usize]
            [tagDiffINChISegments_DIFS_b_SBONDS as usize] = EQL;
        for segment in [
            tagDiffINChISegments_DIFS_t_SATOMS,
            tagDiffINChISegments_DIFS_m_SP3INV,
            tagDiffINChISegments_DIFS_s_STYPE,
        ] {
            expected[tagDiffINChILayers_DIFL_M as usize][segment as usize] = NEQ;
            expected[tagDiffINChILayers_DIFL_F as usize][segment as usize] = EQL;
            expected[tagDiffINChILayers_DIFL_MI as usize][segment as usize] = EQL;
            expected[tagDiffINChILayers_DIFL_FI as usize][segment as usize] = EQL;
        }
        expected[tagDiffINChILayers_DIFL_MI as usize]
            [tagDiffINChISegments_DIFS_i_IATOMS as usize] = NEQ;
        expected[tagDiffINChILayers_DIFL_FI as usize]
            [tagDiffINChISegments_DIFS_i_IATOMS as usize] = FI_MI;
        assert_eq!(segments, expected);

        let first = layers_inchi(&mut heap);
        let mut second = layers_inchi(&mut heap);
        second.szHillFormula = heap.allocate(vec![b'C' as i8, b'3' as i8, 0]).unwrap();
        second.nTotalCharge = 0;
        second.nNum_H_fixed = heap.allocate(vec![0_i8, 0]).unwrap();
        second.Stereo = Default::default();
        second.StereoIsotopic = Default::default();
        second.nNumberOfIsotopicAtoms = 0;
        second.nNumberOfIsotopicTGroups = 0;
        let (_, changed) = run_layers(&mut heap, Some(first), Some(second), 0);
        assert_eq!(changed[tagDiffINChILayers_DIFL_F as usize][0], NEQ);
        assert_eq!(changed[tagDiffINChILayers_DIFL_F as usize][2], 0);
        assert_eq!(changed[tagDiffINChILayers_DIFL_F as usize][3], EMPTY);
        assert_eq!(changed[tagDiffINChILayers_DIFL_F as usize][5], EMPTY);
        assert_eq!(changed[tagDiffINChILayers_DIFL_FI as usize][9], EMPTY);

        let first = layers_inchi(&mut heap);
        let mut deleted = layers_inchi(&mut heap);
        deleted.bDeleted = 1;
        let (_, missing) = run_layers(&mut heap, Some(first), Some(deleted), 0);
        assert_eq!(missing[tagDiffINChILayers_DIFL_F as usize][0], EMPTY);
        assert_eq!(missing[tagDiffINChILayers_DIFL_F as usize][3], EMPTY);

        let first = layers_inchi(&mut heap);
        let mut second = layers_inchi(&mut heap);
        second.IsotopicTGroup = heap
            .allocate(vec![INChI_IsotopicTGroup {
                nTGroupNumber: 1,
                nNum_H: 1,
                nNum_D: 4,
                nNum_T: 1,
            }])
            .unwrap();
        let (difference, _) = run_layers(&mut heap, Some(first), Some(second), 0);
        assert_eq!(difference, 3);

        let mut first = layers_inchi(&mut heap);
        first.nTotalCharge = 0;
        let first_pointer = heap.allocate(vec![first]).unwrap();
        let mut p1 = INCHI_SORT::default();
        p1.pINChI[TAUT_YES as usize] = first_pointer;
        p1.ord_number = 1;
        let mut p2_taut = layers_inchi(&mut heap);
        p2_taut.nTotalCharge = -1;
        let mut p2 = INCHI_SORT::default();
        p2.pINChI[TAUT_YES as usize] = heap.allocate(vec![p2_taut]).unwrap();
        p2.ord_number = 2;
        let mut fixed = [[0_i8; 11]; 4];
        assert_eq!(CompINChILayers(&mut heap, &p1, &p2, &mut fixed, 1), Ok(0));
        assert_eq!(fixed[tagDiffINChILayers_DIFL_F as usize][3], NEQ);
        assert_eq!(fixed[tagDiffINChILayers_DIFL_F as usize][10], NEQ);
    }

    fn stereo_fixture(
        heap: &mut SourceHeap,
        centers: &[u16],
        center_parity: &[i8],
        inversion: i32,
        bond_atom_1: &[u16],
        bond_atom_2: &[u16],
        bond_parity: &[i8],
    ) -> INChI_Stereo {
        INChI_Stereo {
            nNumberOfStereoCenters: centers.len() as i32,
            nNumber: heap.allocate(centers.to_vec()).unwrap(),
            t_parity: heap.allocate(center_parity.to_vec()).unwrap(),
            nCompInv2Abs: inversion,
            nNumberOfStereoBonds: bond_atom_1.len() as i32,
            nBondAtom1: heap.allocate(bond_atom_1.to_vec()).unwrap(),
            nBondAtom2: heap.allocate(bond_atom_2.to_vec()).unwrap(),
            b_parity: heap.allocate(bond_parity.to_vec()).unwrap(),
            ..INChI_Stereo::default()
        }
    }

    #[test]
    fn source_port__ichimake__comparereversedstereoinchi__line_2555() {
        let mut heap = SourceHeap::default();
        assert_eq!(CompareReversedStereoINChI(&heap, None, None), Ok(0));
        let empty = INChI_Stereo::default();
        assert_eq!(CompareReversedStereoINChI(&heap, Some(&empty), None), Ok(0));
        let center_only = INChI_Stereo {
            nNumberOfStereoCenters: 1,
            ..INChI_Stereo::default()
        };
        assert_eq!(
            CompareReversedStereoINChI(&heap, None, Some(&center_only)),
            Ok(20)
        );
        let bond_only = INChI_Stereo {
            nNumberOfStereoBonds: 1,
            ..INChI_Stereo::default()
        };
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&bond_only), None),
            Ok(20)
        );

        let base = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 1, &[2], &[4], &[2]);
        let equal = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 1, &[2], &[4], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&equal)),
            Ok(0)
        );

        let mut changed = equal.clone();
        changed.nNumberOfStereoCenters = 1;
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(21)
        );
        changed = stereo_fixture(&mut heap, &[1, 4], &[1, 2], 1, &[2], &[4], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(22)
        );
        changed = stereo_fixture(&mut heap, &[1, 3], &[1, 1], 1, &[2], &[4], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(23)
        );
        changed = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 2, &[2], &[4], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(24)
        );
        changed.nCompInv2Abs = 0;
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(0)
        );

        changed = equal.clone();
        changed.nNumberOfStereoBonds = 0;
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(25)
        );
        changed = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 1, &[3], &[4], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(26)
        );
        changed = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 1, &[2], &[5], &[2]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(27)
        );
        changed = stereo_fixture(&mut heap, &[1, 3], &[1, 2], 1, &[2], &[4], &[1]);
        assert_eq!(
            CompareReversedStereoINChI(&heap, Some(&base), Some(&changed)),
            Ok(28)
        );
    }

    #[test]
    fn source_port__ichimake__comparereversedstereoinchi2__line_2635() {
        let mut heap = SourceHeap::default();
        let mut empty_results = ICR::default();
        assert_eq!(
            CompareReversedStereoINChI2(&heap, None, None, &mut empty_results),
            Ok(0)
        );
        let one_sided = stereo_fixture(&mut heap, &[1], &[1], 0, &[], &[], &[]);
        assert_eq!(
            CompareReversedStereoINChI2(&heap, Some(&one_sided), None, &mut empty_results,),
            Ok(0)
        );

        let zero_count = INChI_Stereo::default();
        let one_center = stereo_fixture(&mut heap, &[1], &[1], 0, &[], &[], &[]);
        let mut zero_count_results = ICR::default();
        assert_eq!(
            CompareReversedStereoINChI2(
                &heap,
                Some(&zero_count),
                Some(&one_center),
                &mut zero_count_results,
            ),
            Ok(tagInchiDiffBits_IDIF_SC_MISS as i32)
        );
        assert_eq!(
            (
                zero_count_results.num_sc_in2_only,
                zero_count_results.sc_in2_only[0]
            ),
            (1, 0)
        );
        zero_count_results = ICR::default();
        assert_eq!(
            CompareReversedStereoINChI2(
                &heap,
                Some(&one_center),
                Some(&zero_count),
                &mut zero_count_results,
            ),
            Ok(tagInchiDiffBits_IDIF_SC_EXTRA as i32)
        );
        assert_eq!(
            (
                zero_count_results.num_sc_in1_only,
                zero_count_results.sc_in1_only[0]
            ),
            (1, 0)
        );

        let one_bond = stereo_fixture(&mut heap, &[], &[], 0, &[1], &[2], &[1]);
        zero_count_results = ICR::default();
        assert_eq!(
            CompareReversedStereoINChI2(
                &heap,
                Some(&zero_count),
                Some(&one_bond),
                &mut zero_count_results,
            ),
            Ok(tagInchiDiffBits_IDIF_SB_MISS as i32)
        );
        assert_eq!(
            (
                zero_count_results.num_sb_in2_only,
                zero_count_results.sb_in2_only[0]
            ),
            (1, 0)
        );
        zero_count_results = ICR::default();
        assert_eq!(
            CompareReversedStereoINChI2(
                &heap,
                Some(&one_bond),
                Some(&zero_count),
                &mut zero_count_results,
            ),
            Ok(tagInchiDiffBits_IDIF_SB_EXTRA as i32)
        );
        assert_eq!(
            (
                zero_count_results.num_sb_in1_only,
                zero_count_results.sb_in1_only[0]
            ),
            (1, 0)
        );

        let center1 = stereo_fixture(
            &mut heap,
            &[1, 3, 5],
            &[AB_PARITY_UNDF as i8, 1, 2],
            1,
            &[],
            &[],
            &[],
        );
        let center2 = stereo_fixture(
            &mut heap,
            &[2, 3, 6],
            &[AB_PARITY_UNDF as i8, 2, AB_PARITY_UNDF as i8],
            2,
            &[],
            &[],
            &[],
        );
        let mut center_results = ICR::default();
        let center_bits = tagInchiDiffBits_IDIF_SC_PARITY as i32
            | tagInchiDiffBits_IDIF_SC_EXTRA_UNDF as i32
            | tagInchiDiffBits_IDIF_SC_EXTRA as i32
            | tagInchiDiffBits_IDIF_SC_MISS_UNDF as i32
            | tagInchiDiffBits_IDIF_SC_INV as i32;
        assert_eq!(
            CompareReversedStereoINChI2(&heap, Some(&center1), Some(&center2), &mut center_results,),
            Ok(center_bits)
        );
        assert_eq!(
            (
                center_results.num_sc_in1_only,
                &center_results.sc_in1_only[..2],
                center_results.num_sc_undef_in1_only,
                center_results.sc_undef_in1_only[0],
            ),
            (2, &[0, 2][..], 1, 0)
        );
        assert_eq!(
            (
                center_results.num_sc_in2_only,
                &center_results.sc_in2_only[..2],
                center_results.num_sc_undef_in2_only,
                center_results.sc_undef_in2_only[0],
            ),
            (2, &[0, 2][..], 1, 1)
        );
        let recorded_centers = center_results.clone();
        assert_eq!(
            CompareReversedStereoINChI2(&heap, Some(&center1), Some(&center2), &mut center_results,),
            Ok(center_bits)
        );
        assert_eq!(center_results, recorded_centers);

        let bond1 = stereo_fixture(
            &mut heap,
            &[],
            &[],
            0,
            &[1, 3, 5],
            &[2, 4, 6],
            &[AB_PARITY_UNDF as i8, 1, 2],
        );
        let bond2 = stereo_fixture(
            &mut heap,
            &[],
            &[],
            0,
            &[1, 3, 7],
            &[3, 4, 8],
            &[AB_PARITY_UNDF as i8, 2, AB_PARITY_UNDF as i8],
        );
        let mut bond_results = ICR::default();
        let bond_bits = tagInchiDiffBits_IDIF_SB_PARITY as i32
            | tagInchiDiffBits_IDIF_SB_EXTRA_UNDF as i32
            | tagInchiDiffBits_IDIF_SB_EXTRA as i32
            | tagInchiDiffBits_IDIF_SB_MISS_UNDF as i32;
        assert_eq!(
            CompareReversedStereoINChI2(&heap, Some(&bond1), Some(&bond2), &mut bond_results,),
            Ok(bond_bits)
        );
        assert_eq!(
            (
                bond_results.num_sb_in1_only,
                &bond_results.sb_in1_only[..2],
                bond_results.num_sb_undef_in1_only,
                bond_results.sb_undef_in1_only[0],
            ),
            (2, &[0, 2][..], 1, 0)
        );
        assert_eq!(
            (
                bond_results.num_sb_in2_only,
                &bond_results.sb_in2_only[..2],
                bond_results.num_sb_undef_in2_only,
                &bond_results.sb_undef_in2_only[..2],
            ),
            (2, &[0, 2][..], 2, &[1, 3][..])
        );

        let equal1 = stereo_fixture(&mut heap, &[1], &[1], 0, &[2], &[3], &[1]);
        let equal2 = stereo_fixture(&mut heap, &[1], &[1], 0, &[2], &[3], &[1]);
        assert_eq!(
            CompareReversedStereoINChI2(&heap, Some(&equal1), Some(&equal2), &mut ICR::default(),),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichimake__comparereversedinchi2__line_3262() {
        fn compare(
            heap: &mut SourceHeap,
            left: Option<&INChI>,
            right: Option<&INChI>,
            aux_left: Option<&INChI_Aux>,
            aux_right: Option<&INChI_Aux>,
        ) -> (u32, ICR, i32) {
            let mut results = ICR {
                flags: INCHI_MODE::MAX,
                tot_num_H1: -77,
                ..ICR::default()
            };
            let mut error = 77;
            let difference = CompareReversedINChI2(
                heap,
                left,
                right,
                aux_left,
                aux_right,
                &mut results,
                &mut error,
            )
            .unwrap();
            (difference as u32, results, error)
        }

        let mut heap = SourceHeap::default();
        let (difference, results, error) = compare(&mut heap, None, None, None, None);
        assert_eq!((difference, results, error), (0, ICR::default(), 0));

        let one_sided = reversed_inchi2_fixture(&mut heap);
        let (difference, results, error) = compare(&mut heap, Some(&one_sided), None, None, None);
        assert_eq!(difference, tagInchiDiffBits_IDIF_PROBLEM);
        assert_eq!(results.flags as u32, difference);
        assert_eq!(error, 0);

        let left = reversed_inchi2_fixture(&mut heap);
        let right = reversed_inchi2_fixture(&mut heap);
        let (difference, results, error) =
            compare(&mut heap, Some(&left), Some(&right), None, None);
        assert_eq!((difference, results.flags as u32, error), (0, 0, 0));
        assert_eq!((results.tot_num_H1, results.tot_num_H2), (6, 6));

        let mut error_left = left.clone();
        let mut error_right = right.clone();
        error_left.nErrorCode = 4;
        error_right.nErrorCode = 4;
        assert_eq!(
            compare(&mut heap, Some(&error_left), Some(&error_right), None, None).0,
            tagInchiDiffBits_IDIF_PROBLEM
        );
        error_right.nErrorCode = 5;
        assert_eq!(
            compare(&mut heap, Some(&error_left), Some(&error_right), None, None).0,
            tagInchiDiffBits_IDIF_PROBLEM
        );

        let mut changed = reversed_inchi2_fixture(&mut heap);
        changed.nNumberOfAtoms = 1;
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&changed), None, None).0,
            tagInchiDiffBits_IDIF_NUM_AT
        );
        changed = reversed_inchi2_fixture(&mut heap);
        changed.nAtom = heap.allocate(vec![6_u8, 7]).unwrap();
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&changed), None, None).0,
            tagInchiDiffBits_IDIF_ATOMS
        );

        let mut hydrogen_left = reversed_inchi2_fixture(&mut heap);
        let mut hydrogen_right = reversed_inchi2_fixture(&mut heap);
        hydrogen_left.nNum_H = heap.allocate(vec![i8::MIN, 4]).unwrap();
        hydrogen_right.nNum_H = heap.allocate(vec![1_i8, 2]).unwrap();
        hydrogen_left.nNum_H_fixed = heap.allocate(vec![2_i8, -1]).unwrap();
        hydrogen_right.nNum_H_fixed = heap.allocate(vec![1_i8, 3]).unwrap();
        let (difference, results, _) = compare(
            &mut heap,
            Some(&hydrogen_left),
            Some(&hydrogen_right),
            None,
            None,
        );
        assert_eq!(
            difference,
            tagInchiDiffBits_IDIF_POSITION_H
                | tagInchiDiffBits_IDIF_MORE_FH
                | tagInchiDiffBits_IDIF_LESS_FH
        );
        assert_eq!(results.num_diff_pos_H, 2);
        assert_eq!(&results.diff_pos_H_at[..2], &[0, 1]);
        assert_eq!(&results.diff_pos_H_nH[..2], &[127, 2]);
        assert_eq!(
            (
                results.num_fixed_H1_more,
                results.fixed_H_at1_more[0],
                results.fixed_H_nH1_more[0],
            ),
            (1, 0, 1)
        );
        assert_eq!(
            (
                results.num_fixed_H2_more,
                results.fixed_H_at2_more[0],
                results.fixed_H_nH2_more[0],
            ),
            (1, 1, 4)
        );

        let mut fixed_only = reversed_inchi2_fixture(&mut heap);
        fixed_only.nNum_H_fixed = heap.allocate(vec![0_i8, -3]).unwrap();
        let (difference, results, _) =
            compare(&mut heap, Some(&fixed_only), Some(&right), None, None);
        assert_eq!(difference, tagInchiDiffBits_IDIF_MORE_FH);
        assert_eq!(
            (
                results.num_fixed_H1_more,
                results.fixed_H_at1_more[0],
                results.fixed_H_nH1_more[0],
            ),
            (1, 1, -3)
        );
        let (difference, results, _) =
            compare(&mut heap, Some(&right), Some(&fixed_only), None, None);
        assert_eq!(difference, tagInchiDiffBits_IDIF_LESS_FH);
        assert_eq!(results.fixed_H_nH2_more[0], -3);

        let mut fewer_formula_h = reversed_inchi2_fixture(&mut heap);
        fewer_formula_h.szHillFormula = heap
            .allocate(vec![b'C' as i8, b'2' as i8, b'H' as i8, b'4' as i8, 0])
            .unwrap();
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&fewer_formula_h), None, None).0,
            tagInchiDiffBits_IDIF_MORE_H
        );
        assert_eq!(
            compare(&mut heap, Some(&fewer_formula_h), Some(&left), None, None).0,
            tagInchiDiffBits_IDIF_LESS_H
        );
        let mut different_element = reversed_inchi2_fixture(&mut heap);
        different_element.szHillFormula = heap
            .allocate(vec![b'N' as i8, b'2' as i8, b'H' as i8, b'6' as i8, 0])
            .unwrap();
        let (difference, results, _) =
            compare(&mut heap, Some(&left), Some(&different_element), None, None);
        assert_eq!(difference, tagInchiDiffBits_IDIF_NUM_EL);
        assert_eq!((results.tot_num_H1, results.tot_num_H2), (0, 0));

        let mut connection = reversed_inchi2_fixture(&mut heap);
        connection.lenConnTable = 1;
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&connection), None, None).0,
            tagInchiDiffBits_IDIF_CON_LEN
        );
        connection = reversed_inchi2_fixture(&mut heap);
        connection.nConnTable = heap.allocate(vec![1_u16, 3]).unwrap();
        assert_eq!(
            compare(&mut heap, Some(&left), Some(&connection), None, None).0,
            tagInchiDiffBits_IDIF_CON_TBL
        );

        let mut tautomer_left = reversed_inchi2_fixture(&mut heap);
        tautomer_left.nTautomer = heap.allocate(vec![1_u16, 4, 2, 1, 5, 1]).unwrap();
        tautomer_left.lenTautomer = 6;
        let mut tautomer_right = reversed_inchi2_fixture(&mut heap);
        tautomer_right.nTautomer = heap.allocate(vec![1_u16, 4, 1, 0, 2, 5]).unwrap();
        tautomer_right.lenTautomer = 6;
        let (difference, results, _) = compare(
            &mut heap,
            Some(&tautomer_left),
            Some(&tautomer_right),
            None,
            None,
        );
        assert_eq!(
            difference,
            tagInchiDiffBits_IDIF_EXTRA_TG_ENDP
                | tagInchiDiffBits_IDIF_MISS_TG_ENDP
                | tagInchiDiffBits_IDIF_DIFF_TG_ENDP
                | tagInchiDiffBits_IDIF_TG
        );
        assert_eq!(
            (
                results.num_taut_H1,
                results.num_taut_H2,
                results.num_taut_M1,
                results.num_taut_M2,
            ),
            (2, 1, 1, 0)
        );
        assert_eq!(
            (
                results.num_endp_in1_only,
                results.endp_in1_only[0],
                results.num_endp_in2_only,
                results.endp_in2_only[0],
            ),
            (1, 1, 1, 2)
        );

        let mut two_groups = reversed_inchi2_fixture(&mut heap);
        two_groups.nTautomer = heap.allocate(vec![2_u16, 3, 1, 0, 2, 3, 1, 0, 4]).unwrap();
        two_groups.lenTautomer = 9;
        let mut three_groups = reversed_inchi2_fixture(&mut heap);
        three_groups.nTautomer = heap
            .allocate(vec![3_u16, 3, 1, 0, 2, 3, 1, 0, 4, 3, 1, 0, 6])
            .unwrap();
        three_groups.lenTautomer = 13;
        assert_ne!(
            compare(
                &mut heap,
                Some(&tautomer_left),
                Some(&two_groups),
                None,
                None
            )
            .0 & tagInchiDiffBits_IDIF_SINGLE_TG,
            0
        );
        assert_ne!(
            compare(
                &mut heap,
                Some(&two_groups),
                Some(&tautomer_left),
                None,
                None
            )
            .0 & tagInchiDiffBits_IDIF_MULTIPLE_TG,
            0
        );
        assert_ne!(
            compare(
                &mut heap,
                Some(&two_groups),
                Some(&three_groups),
                None,
                None
            )
            .0 & tagInchiDiffBits_IDIF_NUM_TG,
            0
        );
        assert_ne!(
            compare(&mut heap, Some(&right), Some(&tautomer_left), None, None).0
                & tagInchiDiffBits_IDIF_NO_TAUT,
            0
        );
        assert_ne!(
            compare(&mut heap, Some(&tautomer_left), Some(&right), None, None).0
                & tagInchiDiffBits_IDIF_WRONG_TAUT,
            0
        );

        let baseline_allocations = heap.live_allocation_count();
        heap.fail_after_allocations(0);
        let (difference, results, error) = compare(
            &mut heap,
            Some(&tautomer_left),
            Some(&tautomer_right),
            None,
            None,
        );
        assert_eq!(error, -1);
        assert_eq!(results.flags as u32, difference);
        assert_eq!(heap.live_allocation_count(), baseline_allocations);
        heap.fail_after_allocations(1);
        let (difference, results, error) = compare(
            &mut heap,
            Some(&tautomer_left),
            Some(&tautomer_right),
            None,
            None,
        );
        assert_eq!(error, -1);
        assert_eq!(results.flags as u32, difference);
        assert_eq!(heap.live_allocation_count(), baseline_allocations);

        let mut isotope_left = reversed_inchi2_fixture(&mut heap);
        isotope_left.nNumberOfIsotopicAtoms = 1;
        isotope_left.IsotopicAtom = heap
            .allocate(vec![INChI_IsotopicAtom {
                nAtomNumber: 1,
                ..INChI_IsotopicAtom::default()
            }])
            .unwrap();
        let mut isotope_right = reversed_inchi2_fixture(&mut heap);
        assert_eq!(
            compare(
                &mut heap,
                Some(&isotope_left),
                Some(&isotope_right),
                None,
                None
            )
            .0,
            tagInchiDiffBits_IDIF_NUM_ISO_AT
        );
        isotope_right.nNumberOfIsotopicAtoms = 1;
        isotope_right.IsotopicAtom = heap
            .allocate(vec![INChI_IsotopicAtom {
                nAtomNumber: 2,
                ..INChI_IsotopicAtom::default()
            }])
            .unwrap();
        assert_eq!(
            compare(
                &mut heap,
                Some(&isotope_left),
                Some(&isotope_right),
                None,
                None
            )
            .0,
            tagInchiDiffBits_IDIF_ISO_AT
        );

        let mut charged = reversed_inchi2_fixture(&mut heap);
        charged.nTotalCharge = -1;
        assert_eq!(
            compare(&mut heap, Some(&charged), Some(&right), None, None).0,
            tagInchiDiffBits_IDIF_CHARGE
        );
        let aux = INChI_Aux {
            nNumRemovedProtons: 2,
            nNumRemovedIsotopicH: [1, 2, 3],
            ..INChI_Aux::default()
        };
        let (difference, _, _) = compare(&mut heap, Some(&left), Some(&right), Some(&aux), None);
        assert_eq!(
            difference,
            tagInchiDiffBits_IDIF_REM_PROT | tagInchiDiffBits_IDIF_REM_ISO_H
        );

        let regular_left = stereo_fixture(&mut heap, &[1], &[1], 0, &[], &[], &[]);
        let regular_right = stereo_fixture(&mut heap, &[1], &[2], 0, &[], &[], &[]);
        let empty_isotopic_left = INChI_Stereo::default();
        let empty_isotopic_right = INChI_Stereo::default();
        let mut stereo_left = reversed_inchi2_fixture(&mut heap);
        stereo_left.Stereo = heap.allocate(vec![regular_left]).unwrap();
        stereo_left.StereoIsotopic = heap.allocate(vec![empty_isotopic_left]).unwrap();
        let mut stereo_right = reversed_inchi2_fixture(&mut heap);
        stereo_right.Stereo = heap.allocate(vec![regular_right]).unwrap();
        stereo_right.StereoIsotopic = heap.allocate(vec![empty_isotopic_right]).unwrap();
        assert_ne!(
            compare(
                &mut heap,
                Some(&stereo_left),
                Some(&stereo_right),
                None,
                None
            )
            .0 & tagInchiDiffBits_IDIF_SC_PARITY,
            0
        );

        let isotopic_left = stereo_fixture(&mut heap, &[2], &[1], 0, &[], &[], &[]);
        let isotopic_right = stereo_fixture(&mut heap, &[2], &[2], 0, &[], &[], &[]);
        stereo_left.StereoIsotopic = heap.allocate(vec![isotopic_left]).unwrap();
        stereo_right.StereoIsotopic = heap.allocate(vec![isotopic_right]).unwrap();
        assert_ne!(
            compare(
                &mut heap,
                Some(&stereo_left),
                Some(&stereo_right),
                None,
                None
            )
            .0 & tagInchiDiffBits_IDIF_SC_PARITY,
            0
        );
    }

    #[test]
    fn source_port__ichimake__comparereversedinchi__line_2936() {
        let mut heap = SourceHeap::default();
        assert_eq!(CompareReversedINChI(&heap, None, None, None, None), Ok(0));
        let base = inchi_fixture(&mut heap);
        let equal = inchi_fixture(&mut heap);
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), None, None, None),
            Ok(1)
        );
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&equal), None, None),
            Ok(0)
        );

        let mut left = base.clone();
        let mut right = equal.clone();
        left.nErrorCode = 7;
        right.nErrorCode = 7;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(0)
        );
        right.nErrorCode = 8;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(2)
        );

        right = equal.clone();
        right.bDeleted = 1;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(1)
        );
        right = equal.clone();
        right.nNumberOfAtoms = 1;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(3)
        );
        right = equal.clone();
        right.nAtom = heap.allocate(vec![6_u8, 8]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(4)
        );
        right = equal.clone();
        right.szHillFormula = heap.allocate(vec![b'C' as i8, b'H' as i8, 0]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(7)
        );
        right = equal.clone();
        right.nNum_H = heap.allocate(vec![3_i8, 2]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(5)
        );
        left = base.clone();
        left.lenConnTable = 1;
        right = equal.clone();
        right.lenConnTable = 1;
        right.nNum_H = heap.allocate(vec![2_i8, 3]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(6)
        );

        left = base.clone();
        right = equal.clone();
        left.nNum_H_fixed = heap.allocate(vec![1_i8, 0]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(18)
        );
        left = base.clone();
        right.nNum_H_fixed = heap.allocate(vec![0_i8, -1]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(19)
        );
        left.nNum_H_fixed = heap.allocate(vec![2_i8, 0]).unwrap();
        right.nNum_H_fixed = heap.allocate(vec![1_i8, 1]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(20)
        );
        left.nNum_H_fixed = heap.allocate(vec![2_i8, 1]).unwrap();
        right.nNum_H_fixed = heap.allocate(vec![1_i8, 1]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(18)
        );
        left.nNum_H_fixed = heap.allocate(vec![1_i8, 1]).unwrap();
        right.nNum_H_fixed = heap.allocate(vec![2_i8, 1]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(19)
        );

        right = equal.clone();
        right.lenConnTable = 1;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(8)
        );
        right = equal.clone();
        right.nConnTable = heap.allocate(vec![1_u16, 3]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(9)
        );
        right = equal.clone();
        right.lenTautomer = 3;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(10)
        );
        right = equal.clone();
        right.nTautomer = heap.allocate(vec![2_u16, 2]).unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(11)
        );
        right = equal.clone();
        right.nNumberOfIsotopicAtoms = 0;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(12)
        );
        right = equal.clone();
        right.IsotopicAtom = heap
            .allocate(vec![INChI_IsotopicAtom {
                nAtomNumber: 2,
                ..INChI_IsotopicAtom::default()
            }])
            .unwrap();
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(13)
        );
        right = equal.clone();
        right.nTotalCharge = -1;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&right), None, None),
            Ok(14)
        );

        let aux1 = INChI_Aux::default();
        let mut aux2 = INChI_Aux {
            nNumRemovedProtons: 1,
            ..INChI_Aux::default()
        };
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&equal), Some(&aux1), Some(&aux2)),
            Ok(16)
        );
        aux2.nNumRemovedProtons = 0;
        aux2.nNumRemovedIsotopicH = [0, 1, 0];
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&equal), Some(&aux1), Some(&aux2)),
            Ok(17)
        );
        assert_eq!(
            CompareReversedINChI(&heap, Some(&base), Some(&equal), Some(&aux1), None),
            Ok(0)
        );

        let stereo = stereo_fixture(&mut heap, &[1], &[1], 1, &[], &[], &[]);
        let stereo_pointer = heap.allocate(vec![stereo]).unwrap();
        left = base.clone();
        left.Stereo = stereo_pointer;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&equal), None, None),
            Ok(40)
        );
        left = base.clone();
        left.StereoIsotopic = stereo_pointer;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&equal), None, None),
            Ok(60)
        );
        left.Stereo = stereo_pointer;
        right = equal.clone();
        right.Stereo = stereo_pointer;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(0)
        );

        left = base.clone();
        right = equal.clone();
        left.lenTautomer = 1;
        right.lenTautomer = 0;
        assert_eq!(
            CompareReversedINChI(&heap, Some(&left), Some(&right), None, None),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichimake__getinpstructerrortype__line_2480() {
        for error in [i32::MIN, -1, 1, 8, 10, 29] {
            assert_eq!(
                GetInpStructErrorType(None, error, None, i32::MAX).unwrap(),
                _IS_FATAL as i32
            );
        }
        assert_eq!(
            GetInpStructErrorType(None, 9, None, 0).unwrap(),
            _IS_ERROR as i32
        );
        for error in [30, 31, 97, 98, i32::MAX] {
            assert_eq!(
                GetInpStructErrorType(None, error, None, 1).unwrap(),
                _IS_ERROR as i32
            );
        }
        assert_eq!(
            GetInpStructErrorType(None, 0, None, 0).unwrap(),
            _IS_ERROR as i32
        );
        assert_eq!(
            GetInpStructErrorType(None, 0, None, -1).unwrap(),
            _IS_ERROR as i32
        );

        let mut parameters = INPUT_PARMS::default();
        assert_eq!(
            GetInpStructErrorType(Some(&parameters), 98, None, 0).unwrap(),
            _IS_ERROR as i32
        );
        parameters.bAllowEmptyStructure = 1;
        assert_eq!(
            GetInpStructErrorType(Some(&parameters), 98, None, 0).unwrap(),
            _IS_WARNING as i32
        );
        assert_eq!(
            GetInpStructErrorType(None, 98, None, 0),
            Err(SourceHeapError::NullPointer)
        );

        assert_eq!(
            GetInpStructErrorType(None, 0, Some(&[0]), 1).unwrap(),
            _IS_OKAY as i32
        );
        assert_eq!(
            GetInpStructErrorType(None, 0, Some(&[b'x' as i8, 0]), 1).unwrap(),
            _IS_WARNING as i32
        );
        assert_eq!(
            GetInpStructErrorType(None, 0, None, 1),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichimake__comparedfsdescendants4ct__line_2119() {
        let dfs = [1_u16, 2, 3, 4];
        let descendants = [4_u16, 3, 1, 2];
        let deleted = (MAX_ATOMS + 1) as u16;
        let mut order = OrderStruct {
            dfs_number: &dfs,
            number_of_descendants: &descendants,
            current_atom: 0,
        };

        assert_eq!(
            CompareDfsDescendants4CT(deleted, deleted + 1, &order),
            Ok(0)
        );
        assert_eq!(CompareDfsDescendants4CT(deleted, 1, &order), Ok(1));
        assert_eq!(CompareDfsDescendants4CT(1, deleted, &order), Ok(-1));
        assert_eq!(CompareDfsDescendants4CT(1, 2, &order), Ok(2));
        assert_eq!(CompareDfsDescendants4CT(2, 3, &order), Ok(-1));

        order.current_atom = 3;
        assert_eq!(CompareDfsDescendants4CT(0, 1, &order), Ok(-1));
        assert_eq!(CompareDfsDescendants4CT(2, 1, &order), Ok(1));
        order.current_atom = -1;
        assert_eq!(
            CompareDfsDescendants4CT(0, 1, &order),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichimake__getdfsorder4ct__line_2154() {
        let mut heap = SourceHeap::default();
        let single = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let output =
            GetDfsOrder4CT(&mut heap, single, 1, SourceConstPointer::null(), 1, 0).unwrap();
        assert_eq!(heap.slice(output.as_const()).unwrap(), &[1, 0, 0, 0, 0, 0]);
        inchi_free(&mut heap, output).unwrap();

        let chain = heap.allocate_model_storage(vec![2_u16, 1, 3, 2]).unwrap();
        let hydrogens = heap.allocate_model_storage(vec![1_i8, 2, 3]).unwrap();
        let output = GetDfsOrder4CT(&mut heap, chain, 4, hydrogens.as_const(), 3, 0).unwrap();
        assert_eq!(
            heap.slice(output.as_const()).unwrap(),
            &[1, 17, 0, 2, 18, b'-' as u16, 3, 19, b'-' as u16, 0, 0, 0]
        );
        inchi_free(&mut heap, output).unwrap();

        let output = GetDfsOrder4CT(
            &mut heap,
            chain,
            4,
            hydrogens.as_const(),
            3,
            crate::source_types::CT_MODE_PREDECESSORS as i32,
        )
        .unwrap();
        assert_eq!(
            heap.slice(output.as_const()).unwrap(),
            &[
                (MAX_ATOMS + 1) as u16,
                17,
                0,
                1,
                18,
                b',' as u16,
                2,
                19,
                b',' as u16,
                0,
                0,
                0,
            ]
        );
        inchi_free(&mut heap, output).unwrap();

        for failure_after in 0..8 {
            let mut failing_heap = SourceHeap::default();
            let chain = failing_heap
                .allocate_model_storage(vec![2_u16, 1, 3, 2])
                .unwrap();
            failing_heap.fail_after_allocations(failure_after);
            let output = GetDfsOrder4CT(
                &mut failing_heap,
                chain,
                4,
                SourceConstPointer::null(),
                3,
                0,
            )
            .unwrap();
            assert!(output.is_null(), "allocation failure {failure_after}");
        }
    }

    #[test]
    fn source_port__ichimake__inchi_segmentaction__line_1486() {
        for value in i8::MIN..=i8::MAX {
            let promoted = i32::from(value);
            let expected = if promoted & crate::source_types::tagMarkDiff_DIFV_OUTPUT_OMIT_F as i32
                == 0
            {
                crate::source_types::tagINChISegmAction_INCHI_SEGM_OMIT as i32
            } else if promoted & crate::source_types::tagMarkDiff_DIFV_OUTPUT_EMPTY_T as i32 != 0
                && promoted & crate::source_types::tagMarkDiff_DIFV_OUTPUT_EMPTY_F as i32 == 0
            {
                crate::source_types::tagINChISegmAction_INCHI_SEGM_EMPTY as i32
            } else if promoted & crate::source_types::tagMarkDiff_DIFV_OUTPUT_FILL_T as i32 != 0 {
                crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32
            } else {
                crate::source_types::tagINChISegmAction_INCHI_SEGM_OMIT as i32
            };
            assert_eq!(INChI_SegmentAction(value), expected, "value={value}");
        }
    }
}
