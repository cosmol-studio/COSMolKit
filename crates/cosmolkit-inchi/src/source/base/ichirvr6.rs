use crate::source::base::ichimake::{CompareIcr, CompareReversedINChI2};
use crate::source::base::ichirvr1::{
    AddToEdgeList, AllocEdgeList, GetChargeFlowerUpperEdge, MakeOneInChIOutOfStrFromINChI2,
    RemoveForbiddenEdgeMask, RunBnsRestoreOnce, RunBnsTestOnce, SetForbiddenEdgeMask,
};
use crate::source::base::ichirvr4::FillOutExtraFixedHDataRestr;
use crate::source::base::ichirvr5::GetPlusMinusVertex;
use crate::source_types::{
    ALL_TC_GROUPS, BN_DATA, BN_STRUCT, CANON_GLOBALS, EDGE_LIST, EDGE_LIST_CLEAR, EDGE_LIST_FREE,
    ICR, INCHI_CLOCK, INCHI_MODE, INChI, INPUT_PARMS, NO_VERTEX, RI_ERR_ALLOC, RI_ERR_PROGR,
    RI_ERR_SYNTAX, STRUCT_DATA, SourceHeap, SourceHeapError, SourceMutPointer,
    SourceTGroupInfoPointer, StrFromINChI, TAUT_NON, TAUT_YES, VAL_AT, clock_t, inp_ATOM,
};

const INC_ADD_EDGE: i32 = 64;
const IDIF_PROBLEM: INCHI_MODE = 0x0000_0001;
const IDIF_SB_EXTRA_UNDF: INCHI_MODE = 0x0800_0000;
const IDIF_SB_MISS: INCHI_MODE = 0x4000_0000;
const IDIFF_CONSTIT: INCHI_MODE = 0x000f_fffe;
const IDIFF_STEREO: INCHI_MODE = 0x7ff0_0000;
const IDIFF_SB: INCHI_MODE = 0x7c00_0000;

#[rustfmt::skip]
#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FixRestoredStructureStereo(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    ic: SourceMutPointer<INCHI_CLOCK>,
    mut cmp_inchi: INCHI_MODE,
    icr: &mut ICR,
    mut cmp_inchi2: INCHI_MODE,
    icr2: &mut ICR,
    ip: &INPUT_PARMS,
    sd: &STRUCT_DATA,
    pBNS: &mut BN_STRUCT,
    pBD: &mut BN_DATA,
    pStruct: &mut StrFromINChI,
    at: SourceMutPointer<inp_ATOM>,
    at2: SourceMutPointer<inp_ATOM>,
    at3: SourceMutPointer<inp_ATOM>,
    pVA: &mut [VAL_AT],
    pTCGroups: &mut ALL_TC_GROUPS,
    mut ppt_group_info: Option<&mut SourceTGroupInfoPointer>,
    mut ppat_norm: Option<&mut SourceMutPointer<inp_ATOM>>,
    mut ppat_prep: Option<&mut SourceMutPointer<inp_ATOM>>,
    pInChI: [SourceMutPointer<INChI>; 2],
    _num_inp: i64,
    pn_num_run_bns: &mut i32,
    _pn_total_delta: &mut i32,
    forbidden_edge_mask: i32,
    forbidden_stereo_edge_mask: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr6.c:59 FixRestoredStructureStereo
    // INCHI✔❌: complete active source frame follows verbatim; focused branch tests
    // INCHI✔❌: confirm source behavior; SourceHeap allocation-map access adds overhead.
    /*
    int FixRestoredStructureStereo( struct tagCANON_GLOBALS *pCG,
                                    INCHI_CLOCK *ic,
                                    INCHI_MODE cmpInChI,
                                    ICR *icr,
                                    INCHI_MODE cmpInChI2,
                                    ICR *icr2,
                                    ICHICONST INPUT_PARMS *ip,
                                    STRUCT_DATA *sd,
                                    BN_STRUCT *pBNS,
                                    BN_DATA *pBD,
                                    StrFromINChI *pStruct,
                                    inp_ATOM *at,
                                    inp_ATOM *at2,
                                    inp_ATOM *at3,
                                    VAL_AT *pVA,
                                    ALL_TC_GROUPS *pTCGroups,
                                    T_GROUP_INFO **ppt_group_info,
                                    inp_ATOM **ppat_norm,
                                    inp_ATOM **ppat_prep,
                                    INChI *pInChI[],
                                    long num_inp,
                                    int *pnNumRunBNS,
                                    int *pnTotalDelta,
                                    int forbidden_edge_mask,
                                    int forbidden_stereo_edge_mask )
    {
        /*--------- process extra or missing Fixed-H on non-tautomeric atoms ------*/
        /* at2 should be the most recently restored atom, Fixed-H */
        int i, j, k, delta, max_success, cur_success, ret = 0; /* djb-rwth: removing redundant variables */
        int err, iOrigInChI, iRevrInChI;
        int j12, v1, v2, e, vRad;
        BNS_VERTEX *pv1, *pv2, *pvRad;
        BNS_EDGE   *pe, *peRad;
        EDGE_LIST AllChargeEdges, CurrEdges, NFlowerEdges, OtherNFlowerEdges, FixedStereoEdges, AllRadList;
        EDGE_LIST TautMinusEdges[2]; /* 0 -> O & O(+), 1=> N & N(+) */
    
        Vertex     vPathStart, vPathEnd;
        int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
        INChI_Stereo *pStereoInChI, *pStereo2InChI, *pStereoRevrs, *pStereo2Revrs;
    
        /* Stereo */
    
        /* currently being processed layer */
        pStereoInChI = ( pInChI[0]->StereoIsotopic &&
                         pInChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                         pInChI[0]->StereoIsotopic->nNumberOfStereoCenters )
            ? pInChI[0]->StereoIsotopic
            : pInChI[0]->Stereo;
    
        /* mobile-H layer in case of Fixed-H */
        pStereo2InChI = ( pStruct->bMobileH == TAUT_YES || !pInChI[1] ||
                          !pInChI[1]->nNumberOfAtoms || pInChI[1]->bDeleted )
            ? NULL
            : ( pInChI[1]->StereoIsotopic &&
                pInChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                pInChI[1]->StereoIsotopic->nNumberOfStereoCenters ) ?
                pInChI[1]->StereoIsotopic :
                pInChI[1]->Stereo;
    
        /* currently being processed layer */
        pStereoRevrs = ( pStruct->pOneINChI[0]->StereoIsotopic &&
                         pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                         pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters )
            ? pStruct->pOneINChI[0]->StereoIsotopic
            : pStruct->pOneINChI[0]->Stereo;
    
        /* mobile-H layer in case of Fixed-H */
        pStereo2Revrs = ( pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                          !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted )
            ? NULL
            : ( pStruct->pOneINChI[1]->StereoIsotopic &&
                pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters ) ?
                pStruct->pOneINChI[1]->StereoIsotopic :
                pStruct->pOneINChI[1]->Stereo;
    
        INCHI_HEAPCHK
    
        AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &CurrEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &NFlowerEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &FixedStereoEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &AllRadList, EDGE_LIST_CLEAR );
    
        AllocEdgeList( TautMinusEdges + 0, EDGE_LIST_CLEAR );
        AllocEdgeList( TautMinusEdges + 1, EDGE_LIST_CLEAR );
    
        cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
        if (cmpInChI & IDIF_PROBLEM)
        {
            ret = RI_ERR_PROGR; /* severe restore problem */
            goto exit_function;
        }
        if (err)
        {
            ret = RI_ERR_ALLOC;
            goto exit_function;
        }
    
        cmpInChI2 = 0;
    
        if (pStruct->bMobileH == TAUT_NON)
        {
            /* these indexes are used to compare Mobile-H InChI */
            iOrigInChI = ( pInChI[1] && pInChI[1]->nNumberOfAtoms && !pInChI[1]->bDeleted ) ? 1 : 0;
            iRevrInChI = ( pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNumberOfAtoms && !pStruct->pOneINChI[1]->bDeleted ) ? 1 : 0;
        }
        else
        {
            iOrigInChI = 0;
            iRevrInChI = 0;
        }
    
        memset( icr2, 0, sizeof( *icr2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        if (iRevrInChI || iOrigInChI)
        {
            /* additional mobile-H compare in case of Fixed-H */
            cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
            if (cmpInChI & IDIF_PROBLEM)
            {
                ret = RI_ERR_PROGR; /* severe restore problem */
                goto exit_function;
            }
            if (err)
            {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
        }
    
        if (!( cmpInChI & IDIFF_SB ) && !( cmpInChI2 & IDIFF_SB ))
        {
            goto exit_function;
        }
        /* need to temporarily remove fixing of stereogenic bonds */
        for (i = 0; i < pStruct->num_atoms; i++)
        {
            pv1 = pBNS->vert + i;
            for (j = 0; j < at2[i].valence; j++)
            {
                pe = pBNS->edge + ( e = pv1->iedge[j] );
                if (j == pe->neighbor1)
                {
                    /* do not store same bond 2 times */
                    if (( pe->forbidden & forbidden_stereo_edge_mask ) &&
                        ( ret = AddToEdgeList( &FixedStereoEdges, e, INC_ADD_EDGE ) )) /* djb-rwth: ignoring LLVM warning as there should be no memory leak */
                    {
                        /* djb-rwth: fixing coverity ID #499482 */
                        goto exit_function;
                    }
                }
            }
        }
    
        /* djb-rwth: removing redundant code */
        cur_success = 0;
        if (( cmpInChI & IDIF_SB_MISS ) && ( !cmpInChI2 || ( cmpInChI2 & IDIF_SB_MISS ) ) &&
             0 < ( max_success = pBNS->tot_st_cap - pBNS->tot_st_flow ))
        {
            /*----------------------------------------------------*/
            /* case 01: extra stereogenic bond, radical present   */
            /* X=N-O*  => X=N=O and eliminate radical             */
            /*----------------------------------------------------*/
            int aN;
            BNS_VERTEX *pvO, *pvN;
            BNS_EDGE   *peNO;
    
            RemoveForbiddenEdgeMask( pBNS, &FixedStereoEdges, forbidden_stereo_edge_mask );
    
            for (i = 0; i < icr->num_sb_in2_only && cur_success < max_success; i++)
            {
                j12 = icr->sb_in2_only[i];
                pv1 = pBNS->vert + ( v1 = pStereoInChI->nBondAtom1[j12] - 1 );
                pv2 = pBNS->vert + ( v2 = pStereoInChI->nBondAtom2[j12] - 1 );
                for (k = 0; k < at2[v1].valence; k++)
                {
                    pe = pBNS->edge + ( e = pv1->iedge[k] );
                    if (v2 == ( pe->neighbor12 ^ v1 ))
                        break; /* the edge has been found */
                }
                if (k >= at2[v1].valence) /* djb-rwth: addressing LLVM warning */
                {
                    ret = RI_ERR_SYNTAX;
                    goto exit_function;
                }
                /* check v1 */
                pv1->st_edge.cap--;
                pv1->st_edge.flow--;
                pv2->st_edge.flow--;
                pe->flow--; /* new radical on v2 */
                /* djb-rwth: removing redundant code */
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                pv1->st_edge.cap++;
                pv1->st_edge.flow++;
                pv2->st_edge.flow++;
                pe->flow++; /* remove new radical on v2 */
    
                if (ret == 1 /*&& !nDeltaH*/ && !nDeltaCharge && ( v2 == vPathStart || v2 == vPathEnd ))
                {
                    vRad = ( v2 == vPathStart ) ? vPathEnd : vPathStart;
                }
                else
                {
                    pv2->st_edge.cap--;
                    pv2->st_edge.flow--;
                    pv1->st_edge.flow--;
                    pe->flow--; /* new radical on v1 */
                    vRad = NO_VERTEX;
                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                    pv2->st_edge.cap++;
                    pv2->st_edge.flow++;
                    pv1->st_edge.flow++;
                    pe->flow++; /* remove new radical on v1 */
                    if (ret == 1 /*&& !nDeltaH*/ && !nDeltaCharge && ( v1 == vPathStart || v1 == vPathEnd ))
                    {
                        vRad = ( v1 == vPathStart ) ? vPathEnd : vPathStart;
                    }
                }
                if (vRad == NO_VERTEX)
                {
                    continue; /* radical did not affect this bond */
                }
                pvRad = pBNS->vert + vRad;
                /* detect =N-O*  */
                if (pVA[vRad].cNumValenceElectrons == 6 && at2[vRad].valence == 1 &&
                    ( peRad = pBNS->edge + pvRad->iedge[0] )->flow == 0 &&
                     pVA[aN = peRad->neighbor12 ^ vRad].cNumValenceElectrons == 5 &&
                     at2[aN].valence == 2)
                {
                    /*------------------------------------------------------------
                      Fix Metal disconnection/normalization inconsistency :
                                             disconnected  restored
                      R=N(+)-M     R=N--M     R=N  + M     R=N   + M
                        |       ->   ||    ->   ||     ->    |
                        O(-)         O          O            O* <- radical
    
                      The correct     R=N    + M(+)
                      disconnection     |
                      would be this:    O(-)
                    --------------------------------------------------------------*/
                    pvN = pBNS->vert + aN;
                    pvO = pvRad;
                    peNO = peRad;
    
                    /* N-O*  => N=O */
                    peNO->flow++;
                    pvO->st_edge.flow++;
                    pvN->st_edge.cap++;
                    pvN->st_edge.flow++;
                    pBNS->tot_st_cap += 1;
                    pBNS->tot_st_flow += 2;
                    cur_success++;
                }
                else
                {
                         /* all other radicals that affect stereo */
                    delta = pvRad->st_edge.cap - pvRad->st_edge.flow;
                    pvRad->st_edge.cap -= delta;
                    pBNS->tot_st_cap -= delta;
                }
            }
    /*exit_case_01:*/
            SetForbiddenEdgeMask( pBNS, &FixedStereoEdges, forbidden_stereo_edge_mask );
            if (cur_success)
            {
                /* djb-rwth: removing redundant code */
                /* recalculate InChI from the structure */
                if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep ) ))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                /*
                if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                    goto exit_function;
                }
                */
                cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
                if (cmpInChI & IDIF_PROBLEM)
                {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if (err)
                {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
                cmpInChI2 = 0;
                memset( icr2, 0, sizeof( *icr2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                if (iRevrInChI || iOrigInChI)
                {
                    /* additional mobile-H compare in case of Fixed-H */
                    cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                    if (cmpInChI & IDIF_PROBLEM)
                    {
                        ret = RI_ERR_PROGR; /* severe restore problem */
                        goto exit_function;
                    }
                    if (err)
                    {
                        ret = RI_ERR_ALLOC;
                        goto exit_function;
                    }
                }
    
                pStereoRevrs = ( pStruct->pOneINChI[0]->StereoIsotopic &&
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters )
                    ? pStruct->pOneINChI[0]->StereoIsotopic
                    : pStruct->pOneINChI[0]->Stereo;
    
    
                pStereo2Revrs = ( pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                           !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted )
                    ? NULL
                    : ( pStruct->pOneINChI[1]->StereoIsotopic &&
                        pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                        pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters ) ?
                        pStruct->pOneINChI[1]->StereoIsotopic :
                        pStruct->pOneINChI[1]->Stereo;
            }
        }
    
        /* djb-rwth: removing redundant code */
        if (!( cmpInChI & IDIF_SB_MISS ) && ( cmpInChI2 & IDIF_SB_MISS ) &&
             icr2->num_sb_in2_only &&
             0 < ( pBNS->tot_st_cap - pBNS->tot_st_flow )) /* djb-rwth: removing redundant code */
        {
            /*----------------------------------------------------*/
            /* case 02: missing stereogenic bond in Mobile-H only */
            /* X=N-O*  => X=N=O and eliminate radical             */
            /*----------------------------------------------------*/
            int retC, ret2C, retS, ret2S;
            /* djb-rwth: removing redundant variables */
            ICR  icr_Prev, icr2_Prev;
    
            /* blind attepmt */
            icr_Prev = *icr;
            icr2_Prev = *icr2;
            /* djb-rwth: removing redundant code */
            for (i = AllRadList.num_edges = 0; i < pStruct->num_atoms; i++)
            {
                if (pBNS->vert[i].st_edge.cap - pBNS->vert[i].st_edge.flow == 1 &&
                    ( ret = AddToEdgeList( &AllRadList, i, INC_ADD_EDGE ) ))
                {
                    goto exit_function;
                }
            }
            for (i = 0; i < AllRadList.num_edges; i++)
            {
                j = AllRadList.pnEdges[i];
                pBNS->vert[j].st_edge.cap -= 1;
                pBNS->tot_st_cap -= 1;
            }
            /*-------------------------------------------------*/
            /* re-create InChI and see whether it looks better */
            /*-------------------------------------------------*/
            if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                ppt_group_info, ppat_norm, ppat_prep ) ))
            {
                goto exit_function;
            }
            if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
            /*
            if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            */
            cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
            if (cmpInChI & IDIF_PROBLEM)
            {
                ret = RI_ERR_PROGR; /* severe restore problem */
                goto exit_function;
            }
            if (err)
            {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            cmpInChI2 = 0;
            memset( icr2, 0, sizeof( *icr2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            if (iRevrInChI || iOrigInChI)
            {
                /* additional mobile-H compare in case of Fixed-H */
                cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                if (cmpInChI & IDIF_PROBLEM)
                {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if (err)
                {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
            }
            retC = CompareIcr( icr, &icr_Prev, NULL, NULL, IDIFF_CONSTIT );
            retS = CompareIcr( icr, &icr_Prev, NULL, NULL, IDIFF_STEREO );
            ret2C = CompareIcr( icr2, &icr2_Prev, NULL, NULL, IDIFF_CONSTIT );
            ret2S = CompareIcr( icr2, &icr2_Prev, NULL, NULL, IDIFF_STEREO );
    
            if (0 >= retC &&
                 0 >= retS &&
                 0 >= ret2C &&
                 0 > ret2S)
            {
                ; /* accept */
            }
            else
            {
                 /* reject */
                for (i = 0; i < AllRadList.num_edges; i++)
                {
                    j = AllRadList.pnEdges[i];
                    pBNS->vert[j].st_edge.cap += 1;
                    pBNS->tot_st_cap += 1;
                }
    
                /*-------------------------------------------------*/
                /* re-create InChI-- return to previous state      */
                /*-------------------------------------------------*/
                if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep ) ))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                /*
                if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                    goto exit_function;
                }
                */
                cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
                if (cmpInChI & IDIF_PROBLEM)
                {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if (err)
                {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
                cmpInChI2 = 0;
                memset( icr2, 0, sizeof( *icr2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                if (iRevrInChI || iOrigInChI)
                {
                    /* additional mobile-H compare in case of Fixed-H */
                    cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                    if (cmpInChI & IDIF_PROBLEM)
                    {
                        ret = RI_ERR_PROGR; /* severe restore problem */
                        goto exit_function;
                    }
                    if (err)
                    {
                        ret = RI_ERR_ALLOC;
                        goto exit_function;
                    }
                }
                pStereoRevrs = ( pStruct->pOneINChI[0]->StereoIsotopic &&
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters ) ?
                    pStruct->pOneINChI[0]->StereoIsotopic : pStruct->pOneINChI[0]->Stereo;
    
    
                pStereo2Revrs = ( pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                           !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted ) ?
                    NULL :
                    ( pStruct->pOneINChI[1]->StereoIsotopic &&
                    pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                    pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters ) ?
                    pStruct->pOneINChI[1]->StereoIsotopic :
                    pStruct->pOneINChI[1]->Stereo;
            }
    /*exit_case_02:;*/
        }
    
        cur_success = 0;
        if (pStruct->bMobileH == TAUT_NON && ( cmpInChI & IDIF_SB_EXTRA_UNDF ) &&
             pStruct->endpoint)
        {
            /*------------------------------------------------------*/
            /* case 03: extra stereogenic bond in Fixed-H  only     */
            /* in Mobile-H this bond is not stereogenic.            */
            /* Since this bond parity is not known, it is UNDEFINED */
            /*------------------------------------------------------*/
            int bDone, num_endpoints;
    
            TautMinusEdges[0].num_edges = 0;
            TautMinusEdges[1].num_edges = 0;
            AllChargeEdges.num_edges = 0;
            /* in1 => in restored structure; in2 => in original InChI */
            for (i = 0; i < icr->num_sb_undef_in1_only; i++)
            {
                j12 = icr->sb_undef_in1_only[i];
                pv1 = pBNS->vert + ( v1 = pStereoRevrs->nBondAtom1[j12] - 1 );
                pv2 = pBNS->vert + ( v2 = pStereoRevrs->nBondAtom2[j12] - 1 ); /* djb-rwth: ignoring LLVM warning: variable used */
    
                if (pStereo2Revrs)
                {
                    /* reject if it is extra in Mobile-H also */
                    if (icr2->num_sb_undef_in1_only)
                    {
                        for (j = 0; j < icr2->num_sb_undef_in1_only; j++)
                        {
                            k = icr2->sb_undef_in1_only[j];
                            if (v1 == pStereo2Revrs->nBondAtom1[k] &&
                                 v2 == pStereo2Revrs->nBondAtom2[k])
                            {
                                break;
                            }
                        }
                        if (j < icr->num_sb_in1_only)
                        {
                            continue; /* extra stereobond in Mobile H also */
                        }
                    }
                }
                /* reject if it is a stereobond in Mobile-H also */
                if (pStereo2InChI && pStereo2InChI->nNumberOfStereoBonds)
                {
                    for (j = 0; j < pStereo2InChI->nNumberOfStereoBonds; j++)
                    {
                        if (v1 == pStereo2InChI->nBondAtom1[j] &&
                             v2 == pStereo2InChI->nBondAtom1[j])
                        {
                            break;
                        }
                    }
                    if (j < pStereo2InChI->nNumberOfStereoBonds)
                    {
                        continue; /* ignore this extra stereo bond: it is in Mobile-H */
                    }
                }
                /* find the edge between v1 and v2 */
                for (k = 0; k < at2[v1].valence; k++)
                {
                    pe = pBNS->edge + ( e = pv1->iedge[k] );
                    if (v2 == ( pe->neighbor12 ^ v1 ))
                        break; /* the edge has been found */
                }
                if (k >= at2[v1].valence) /* djb-rwth: addressing LLVM warning */
                {
                    ret = RI_ERR_SYNTAX;
                    goto exit_function;
                }
                /* Fix all charges except negative charges on tautomeric endpoints */
                if (!AllChargeEdges.num_edges && !TautMinusEdges[0].num_edges && !TautMinusEdges[1].num_edges)
                {
                    for (j = 0; j < pStruct->num_atoms; j++)
                    {
                        if (( k = pVA[j].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[k].forbidden)
                        {
                            if (!pStruct->endpoint[j])
                            {
                                if ((ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_function;
                                }
                            }
                            else
                                if (pVA[j].cNumValenceElectrons == 6)
                                {
                                    /* O */
                                    if ((ret = AddToEdgeList( TautMinusEdges + 0, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                    {
                                        goto exit_function;
                                    }
                                }
                                else
                                {
                                                     /* N */
                                    if ((ret = AddToEdgeList( TautMinusEdges + 1, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                    {
                                        goto exit_function;
                                    }
                                }
                        }
                        if (( k = pVA[j].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[k].forbidden)
                        {
                            if ((ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                            /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
                            if (pVA[j].cNumValenceElectrons == 5 && !pVA[j].cMetal && /* N, P, As */
                                 NO_VERTEX != ( k = GetChargeFlowerUpperEdge( pBNS, pVA, k ) ))
                            {
    
                                if (!pBNS->edge[j].forbidden && pBNS->edge[k].flow)
                                {
                                    if ((ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                    {
                                        goto exit_function;
                                    }
                                }
                            }
                        }
                    }
                }
                if (!pe->flow)
                    continue;
                /* fix all charges except tautomeric; first allow only O, then only N, finally both N and O */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                for (k = 1, bDone = 0; k < 4 && !bDone; k++)
                {
                    /* fix tautomeric charges */
                    num_endpoints = ( TautMinusEdges + 0 )->num_edges + ( TautMinusEdges + 1 )->num_edges;
                    if (k == 2)
                    {
                        /* fix charges on O */
                        SetForbiddenEdgeMask( pBNS, TautMinusEdges + 0, forbidden_edge_mask );
                        num_endpoints -= ( TautMinusEdges + 0 )->num_edges;
                    }
                    if (k == 1)
                    {
                        SetForbiddenEdgeMask( pBNS, TautMinusEdges + 1, forbidden_edge_mask );
                        num_endpoints -= ( TautMinusEdges + 1 )->num_edges;
                    }
                    if (num_endpoints >= 2)
                    {
                        delta = 1;
                        pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                        pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );
    
                        pe->forbidden |= forbidden_edge_mask; /* fix stereobond */
                        pe->flow -= delta;                    /* decrement stereobond order */
                        pv1->st_edge.flow -= delta;
                        pv2->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2 * delta;
    
                        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
    
                        if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                            (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 0) /* djb-rwth: addressing LLVM warnings */
                        {
                            /* Negative charge has been moved, no change in number of charges */
                            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                            if (ret > 0)
                            {
                                ( *pnNumRunBNS )++;
                                cur_success++; /* 01 */
                                bDone = 1;
                            }
                        }
                        else
                        {
                            pe->forbidden &= ~forbidden_edge_mask;
                            pe->flow += delta;
                            pv1->st_edge.flow += delta;
                            pv2->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2 * delta;
                        }
                    }
                    /* unfix tautomeric charges */
                    if (k == 2)
                        RemoveForbiddenEdgeMask( pBNS, TautMinusEdges + 0, forbidden_edge_mask );
                    if (k == 1)
                        RemoveForbiddenEdgeMask( pBNS, TautMinusEdges + 1, forbidden_edge_mask );
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            }
    /*exit_case_03:*/
            if (cur_success)
            {
                /* djb-rwth: removing redundant code */
                /* recalculate InChI from the structure */
                if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep ) ))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                /*
                if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                    goto exit_function;
                }
                */
                cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
                if (cmpInChI & IDIF_PROBLEM)
                {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if (err)
                {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
                cmpInChI2 = 0;
                memset( icr2, 0, sizeof( *icr2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                if (iRevrInChI || iOrigInChI)
                {
                    /* additional mobile-H compare in case of Fixed-H */
                    cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                    if (cmpInChI & IDIF_PROBLEM)
                    {
                        ret = RI_ERR_PROGR; /* severe restore problem */
                        goto exit_function;
                    }
                    if (err)
                    {
                        ret = RI_ERR_ALLOC;
                        goto exit_function;
                    }
                }
                pStereoRevrs = ( pStruct->pOneINChI[0]->StereoIsotopic &&
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters )
                    ? pStruct->pOneINChI[0]->StereoIsotopic
                    : pStruct->pOneINChI[0]->Stereo;
    
    
                pStereo2Revrs = ( pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                           !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted )
                    ? NULL
                    : ( pStruct->pOneINChI[1]->StereoIsotopic &&
                        pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                        pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters ) ?
                        pStruct->pOneINChI[1]->StereoIsotopic :
                        pStruct->pOneINChI[1]->Stereo;
            }
        }
    
        cur_success = 0;
        if (( cmpInChI & IDIF_SB_EXTRA_UNDF ))
        {
            /*------------------------------------------------------*/
            /* case 04: extra stereogenic bond                      */
            /* Since this bond parity is not known, it is UNDEFINED */
            /*------------------------------------------------------*/
            int bDone, num_endpoints;
    
            TautMinusEdges[0].num_edges = 0;
            TautMinusEdges[1].num_edges = 0;
            AllChargeEdges.num_edges = 0;
            /* in1 => in restored structure; in2 => in original InChI */
            for (i = 0; i < icr->num_sb_undef_in1_only; i++)
            {
                j12 = icr->sb_undef_in1_only[i];
                pv1 = pBNS->vert + ( v1 = pStereoRevrs->nBondAtom1[j12] - 1 );
                pv2 = pBNS->vert + ( v2 = pStereoRevrs->nBondAtom2[j12] - 1 ); /* djb-rwth: ignoring LLVM warning: variable used */
                
                /* djb-rwth: fixing oss-fuzz issue #67650 */
                pe = pBNS->edge + (e = pv1->iedge[0]); /* djb-rwth: proper initialisation required to avoid garbage values */
                /* find the edge between v1 and v2 */
                for (k = 0; k < at2[v1].valence; k++)
                {
                    pe = pBNS->edge + ( e = pv1->iedge[k] );
                    if (v2 == ( pe->neighbor12 ^ v1 ))
                        break; /* the edge has been found */
                }
                if (k == at2[v1].valence)
                {
                    ret = RI_ERR_SYNTAX;
                    goto exit_function;
                }
                if (pStereo2Revrs)
                {
                    /* reject if it is not extra in Mobile-H also */
                    if (icr2->num_sb_undef_in1_only)
                    {
                        for (j = 0; j < icr2->num_sb_undef_in1_only; j++)
                        {
                            k = icr2->sb_undef_in1_only[j];
                            if (v1 == pStereo2Revrs->nBondAtom1[k] &&
                                 v2 == pStereo2Revrs->nBondAtom2[k])
                            {
                                break;
                            }
                        }
                        if (j == icr->num_sb_in1_only)
                        {
                            continue; /* extra stereobond only in Fixed-H, not in Mobile H also */
                        }
                    }
                }
    
                /* Fix all charges except negative charges on tautomeric endpoints */
                if (!AllChargeEdges.num_edges && !TautMinusEdges[0].num_edges && !TautMinusEdges[1].num_edges)
                {
                    for (j = 0; j < pStruct->num_atoms; j++)
                    {
                        if (( k = pVA[j].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[k].forbidden)
                        {
                            if ((ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                        if (( k = pVA[j].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[k].forbidden)
                        {
                            int bMayBeUnfixed = !at2[j].num_H && !( pStruct->endpoint && pStruct->endpoint[j] );
                            if ((bMayBeUnfixed && pVA[j].cNumValenceElectrons == 6) ||
                                 (pVA[j].cNumValenceElectrons == 5 && pVA[j].cPeriodicRowNumber > 1)) /* djb-rwth: addressing LLVM warning */
                            {
                                /* O & P */
                                if ((ret = AddToEdgeList( TautMinusEdges + 0, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_function;
                                }
                            }
                            else
                            {
                                if (bMayBeUnfixed &&
                                        pVA[j].cNumValenceElectrons == 5 && pVA[j].cPeriodicRowNumber == 1)
                                {
                                    /* N */
                                    if ((ret = AddToEdgeList( TautMinusEdges + 1, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                    {
                                        goto exit_function;
                                    }
                                }
                                else
                                {
                                    if ((ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                    {
                                        goto exit_function;
                                    }
                                }
                            }
                            /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
                            if (pVA[j].cNumValenceElectrons == 5 && !pVA[j].cMetal && /* N, P, As */
                                 NO_VERTEX != ( k = GetChargeFlowerUpperEdge( pBNS, pVA, k ) ))
                            {
                                if (!pBNS->edge[j].forbidden && pBNS->edge[k].flow)
                                {
                                    if ((ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                    {
                                        goto exit_function;
                                    }
                                }
                            }
                        }
                    }
                }
                if (!pe->flow)
                    continue;
                /* fix all charges except tautomeric; first allow only O, then only N, finally both N and O */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                for (k = 1, bDone = 0; k < 4 && !bDone; k++)
                {
                    /* fix positive charges on heteroatoms */
                    num_endpoints = ( TautMinusEdges + 0 )->num_edges + ( TautMinusEdges + 1 )->num_edges;
                    if (k == 2)
                    {
                        /* fix charges on O */
                        SetForbiddenEdgeMask( pBNS, TautMinusEdges + 0, forbidden_edge_mask );
                        num_endpoints -= ( TautMinusEdges + 0 )->num_edges;
                    }
                    if (k == 1)
                    {
                        /* fix charges on N */
                        SetForbiddenEdgeMask( pBNS, TautMinusEdges + 1, forbidden_edge_mask );
                        num_endpoints -= ( TautMinusEdges + 1 )->num_edges;
                    }
                    if (num_endpoints >= 2)
                    {
                        delta = 1;
                        pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                        pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );
    
                        pe->forbidden |= forbidden_edge_mask; /* fix stereobond */
                        pe->flow -= delta;                    /* decrement stereobond order */
                        pv1->st_edge.flow -= delta;
                        pv2->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2 * delta;
    
                        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
    
                        if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                            (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 0) /* djb-rwth: addressing LLVM warnings */
                        {
                            /* Negative charge has been moved, no change in number of charges */
                            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                            if (ret > 0)
                            {
                                ( *pnNumRunBNS )++;
                                cur_success++; /* 01 */
                                bDone = 1;
                            }
                        }
                        else
                        {
                            pe->forbidden &= ~forbidden_edge_mask;
                            pe->flow += delta;
                            pv1->st_edge.flow += delta;
                            pv2->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2 * delta;
                        }
                    }
                    /* unfix tautomeric charges */
                    if (k == 2)
                        RemoveForbiddenEdgeMask( pBNS, TautMinusEdges + 0, forbidden_edge_mask );
                    if (k == 1)
                        RemoveForbiddenEdgeMask( pBNS, TautMinusEdges + 1, forbidden_edge_mask );
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            }
    /*exit_case_04:*/
            if (cur_success)
            {
                /* djb-rwth: removing redundant code */
                /* recalculate InChI from the structure */
                if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep ) ))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                /*
                if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                    goto exit_function;
                }
                */
                cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
                if (cmpInChI & IDIF_PROBLEM)
                {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if (err)
                {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
                cmpInChI2 = 0;
                memset( icr2, 0, sizeof( *icr2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                if (iRevrInChI || iOrigInChI)
                {
    /* additional mobile-H compare in case of Fixed-H */
                    cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                    if (cmpInChI & IDIF_PROBLEM)
                    {
                        ret = RI_ERR_PROGR; /* severe restore problem */
                        goto exit_function;
                    }
                    if (err)
                    {
                        ret = RI_ERR_ALLOC;
                        goto exit_function;
                    }
                }
                pStereoRevrs = ( pStruct->pOneINChI[0]->StereoIsotopic &&
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters )
                    ? pStruct->pOneINChI[0]->StereoIsotopic
                    : pStruct->pOneINChI[0]->Stereo;
    
    
                pStereo2Revrs = ( pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                           !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted )
                    ? NULL
                    : ( pStruct->pOneINChI[1]->StereoIsotopic &&
                        pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                        pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters ) ?
                        pStruct->pOneINChI[1]->StereoIsotopic :
                        pStruct->pOneINChI[1]->Stereo;
            }
        }
    
        cur_success = 0;
        if (pStruct->bMobileH == TAUT_YES &&
            ( cmpInChI & IDIF_SB_EXTRA_UNDF &&
                !pStruct->ti.num_t_groups )
             /*pStruct->bMobileH == TAUT_NON && (cmpInChI2 & IDIF_SB_EXTRA_UNDF)*/)
        {
            /*----------------------------------------------------------*/
            /* case 05: extra stereogenic bond on =NH2(+), (B, Mobile-H)*/
            /*                   H             H                        */
            /* original: N(+)=-N<   ->  N--==N/                         */
            /* (A)               H                                      */
            /*                             double bond is marked as     */
            /*                             not stereogenic due to       */
            /*                             its change during proton     */
            /*                             removal => No Stereo bond    */
            /*                             (=NH may be tautomeric)      */
            /*                                                          */
            /*                   H             H                        */
            /* original: N=-N(+)<   ->  N--==N/                         */
            /* (B)               H                                      */
            /*                             double bond was not          */
            /*                             changed during proton        */
            /* In Fixed-H this bond        removal => Undef Stereo      */
            /* may not be stereogenic      (=NH is not tautomeric)      */
            /* (a) due to (+) movement                                  */
            /* (b) due to symmetry (2H), even if isotopic               */
            /*                                                          */
            /* Fixed-H: move (+) to or from NH2 for Undef or No stereo  */
            /*          respectively                                    */
            /* Mobile-H: Add H(+) to =NH and move the charge to =N-     */
            /*           to eliminate Undef stereo                      */
            /*          Move charge from N to -NH2 to create            */
            /*           Undef Stereo                                   */
            /* Since this bond parity is not known, it is UNDEFINED     */
            /*                                                          */
            /* Solution: Add H(+) to =NH and move charge to -N=         */
            /*                                                          */
            /*----------------------------------------------------------*/
            int aN, aC, i1, i2, vPlusMinus;
            AllChargeEdges.num_edges = 0;
            /* in1 => in restored structure; in2 => in original InChI */
            for (i = 0; i < icr->num_sb_undef_in1_only; i++)
            {
                j12 = icr->sb_undef_in1_only[i];
                pv1 = pBNS->vert + ( v1 = pStereoRevrs->nBondAtom1[j12] - 1 );
                pv2 = pBNS->vert + ( v2 = pStereoRevrs->nBondAtom2[j12] - 1 ); /* djb-rwth: ignoring LLVM warning: variable used */
                /* indicators of -NH: */
                i1 = at2[v1].valence == 1 && at2[v1].num_H == 1 && !at2[v1].endpoint &&
                    pVA[v1].cNumValenceElectrons == 5 && pVA[v1].cPeriodicRowNumber == 1;
                i2 = at2[v2].valence == 1 && at2[v2].num_H == 1 && !at2[v2].endpoint &&
                    pVA[v2].cNumValenceElectrons == 5 && pVA[v2].cPeriodicRowNumber == 1;
                if ((!i1 && !i2) || (i1 && i2)) /* djb-rwth: addressing LLVM warnings */
                {
                    continue;
                }
                /* find the edge between v1 and v2 */
                for (k = 0; k < at2[v1].valence; k++)
                {
                    pe = pBNS->edge + ( e = pv1->iedge[k] );
                    if (v2 == ( pe->neighbor12 ^ v1 ))
                    {
                        break; /* the edge has been found */
                    }
                }
                if (k == at2[v1].valence)
                {
                    ret = RI_ERR_SYNTAX;
                    goto exit_function;
                }
                if (pe->flow != 1)
                {
                    continue; /* already charged */
                }
                aN = i1 ? v1 : v2; /* -NH atom */
                aC = i1 ? v2 : v1; /* neighbor */
                /* Replace =NH with -NH2
                   Create such a charge on some -N< that may be moved to NH2 to remove H(+):
                   transformation:
                   from:  HN=C-=-N=(+vert)-Y=(+super)-(+/-)
                   to:   2HN-C*-=-N=(+vert)-Y=(+super)-(+/-)*
                   Run BNS to obtain:
                         2HN-C=-=N(+)-(+vert)=Y-(+super)=(+/-)
                */
                vPlusMinus = GetPlusMinusVertex( pBNS, pTCGroups, 1, 0 );
                if (NO_VERTEX == vPlusMinus)
                {
                    break; /* cannot do anything */
                }
                /* increase edges to -Y-(+/-)-Y- capacities */
                delta = 1;
                for (i1 = 0; i1 < pBNS->vert[vPlusMinus].num_adj_edges; i1++)
                {
                    i2 = pBNS->edge[pBNS->vert[vPlusMinus].iedge[i1]].neighbor12 ^ vPlusMinus;
                    for (k = 0; k < pBNS->vert[i2].num_adj_edges; k++)
                    {
                        pBNS->edge[pBNS->vert[i2].iedge[k]].cap += delta;
                    }
                }
                /* Fix all charges except (+) on -N< */
                if (!AllChargeEdges.num_edges)
                {
                    for (j = 0; j < pStruct->num_atoms; j++)
                    {
                        if (( k = pVA[j].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[k].forbidden)
                        {
                            if ((ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                        if (( k = pVA[j].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[k].forbidden)
                        {
                            if (pVA[j].cNumValenceElectrons == 5 && pVA[j].cPeriodicRowNumber == 1 &&
                                 !at2[j].num_H && at2[j].valence == 3 &&
                                 !( at2[j].endpoint || (pStruct->endpoint && pStruct->endpoint[j]) )) /* djb-rwth: addressing LLVM warning */
                            {
                                ; /* do not fix -N< or =N(+)< */
                            }
                            else
                            {
                                                 /* all others */
                                if ((ret = AddToEdgeList( TautMinusEdges + 0, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_function;
                                }
                            }
                            /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
                            if (pVA[j].cNumValenceElectrons == 5 && !pVA[j].cMetal && /* N, P, As */
                                 NO_VERTEX != ( k = GetChargeFlowerUpperEdge( pBNS, pVA, k ) ))
                            {
                                if (!pBNS->edge[j].forbidden && pBNS->edge[k].flow)
                                {
                                    if ((ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                    {
                                        goto exit_function;
                                    }
                                }
                            }
                        }
                    }
                }
                /* Make bond to =NH single, add radical to aC */
                pe->flow -= delta; /* make single bond */
                pBNS->vert[aN].st_edge.flow -= delta;
                pBNS->vert[aN].st_edge.cap -= delta; /* avoid radical on N */
                pBNS->vert[aC].st_edge.flow -= delta; /* create radical on C */
                pBNS->vert[vPlusMinus].st_edge.cap += delta; /* create radical on (+/-) */
                pBNS->tot_st_flow -= 2 * delta;
                /* fix C-NH bond */
                if ((ret = AddToEdgeList( &AllChargeEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                /* pBNS->tot_st_cap is unchanged */
                /* find all aC edges except pe to fix them */
                /* 2. Check whether it would work and do if it would */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );/* fix aC edges */
                pe->cap++;
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
    
                if (ret == 1 && ( (vPathEnd == vPlusMinus && vPathStart == aC) ||
                    (vPathEnd == aC && vPathStart == vPlusMinus) ) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warnings */
                {
                    /* Negative charge has been moved, no change in number of charges */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if (ret > 0)
                    {
                        ( *pnNumRunBNS )++;
                        /* 3. Add  H to -NH and register increaded charge */
                        pStruct->at[aN].num_H++;
                        pTCGroups->total_charge++;
                        cur_success++; /* 01 */
                    }
                }
                else
                {
                    pe->flow += delta; /* make single bond */
                    pBNS->vert[aN].st_edge.flow += delta;
                    pBNS->vert[aN].st_edge.cap += delta; /* avoid radical on N */
                    pBNS->vert[aC].st_edge.flow += delta; /* create radical on C */
                    pBNS->vert[vPlusMinus].st_edge.cap -= delta; /* create radical on (+/-) */
                    pBNS->tot_st_flow += 2 * delta;
                    RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );/* fix aC edges */
                    AllChargeEdges.num_edges--; /* remove pe from the list */
                    CurrEdges.num_edges = 0;
                    continue; /* should not happen */
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );/* fix aC edges */
                AllChargeEdges.num_edges--; /* remove pe from the list */
                CurrEdges.num_edges = 0;
            }
    /*exit_case_05:*/
            if (cur_success)
            {
                /* djb-rwth: removing redundant code */
                /* recalculate InChI from the structure */
                if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep ) ))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                /*
                if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                    goto exit_function;
                }
                */
                cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
                if (cmpInChI & IDIF_PROBLEM)
                {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if (err)
                {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
                cmpInChI2 = 0;
                memset( icr2, 0, sizeof( *icr2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                if (iRevrInChI || iOrigInChI)
                {
                    /* additional mobile-H compare in case of Fixed-H */
                    cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                    if (cmpInChI & IDIF_PROBLEM)
                    {
                        ret = RI_ERR_PROGR; /* severe restore problem */
                        goto exit_function;
                    }
                    if (err)
                    {
                        ret = RI_ERR_ALLOC;
                        goto exit_function;
                    }
                }
                pStereoRevrs = ( pStruct->pOneINChI[0]->StereoIsotopic &&
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters )
                    ? pStruct->pOneINChI[0]->StereoIsotopic
                    : pStruct->pOneINChI[0]->Stereo; /* djb-rwth: ignoring LLVM warning: variable used */
    
    
                pStereo2Revrs = ( pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                           !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted )
                    ? NULL
                    : ( pStruct->pOneINChI[1]->StereoIsotopic &&
                        pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                        pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters ) ?
                        pStruct->pOneINChI[1]->StereoIsotopic :
                        pStruct->pOneINChI[1]->Stereo;
            }
        }
    
        cur_success = 0;
        if (pStruct->bMobileH == TAUT_NON && pStereo2Revrs /* added check 2006-04-05 */ &&
            ( cmpInChI2 & IDIF_SB_EXTRA_UNDF &&
                !pStruct->ti.num_t_groups )
             /*pStruct->bMobileH == TAUT_NON && (cmpInChI2 & IDIF_SB_EXTRA_UNDF)*/)
        {
            /*----------------------------------------------------------*/
            /* case 06: extra stereogenic bond on =NH2(+), (B, Fixed-H) */
            /*                   H                H        ===========  */
            /* original: N(+)=-N<   ->  N--==N(+)<                      */
            /* (A)               H                H                     */
            /*                             double bond in Mobile-H      */
            /*                             layer has Undef stereo       */
            /*                                                          */
            /*                                                          */
            /* Fixed-H: move (+) to or from NH2 for Undef or No stereo  */
            /*          respectively                                    */
            /* Mobile-H: Add H(+) to =NH and move the charge to =N-     */
            /*           to eliminate Undef stereo                      */
            /*          Move charge from N to -NH2 to create            */
            /*           Undef Stereo                                   */
            /* Since this bond parity is not known, it is UNDEFINED     */
            /*                                                          */
            /* Solution: Move (+) from -NH2(+) to othe -N<              */
            /*                                                          */
            /*----------------------------------------------------------*/
            int aN, i1, i2, ePlus; /* djb-rwth: removing redundant variables */
            BNS_EDGE   *pePlus;
            AllChargeEdges.num_edges = 0;
            /* in1 => in restored structure; in2 => in original InChI */
            for (i = 0; i < icr2->num_sb_undef_in1_only; i++)
            {
                j12 = icr2->sb_undef_in1_only[i];
                pv1 = pBNS->vert + ( v1 = pStereo2Revrs->nBondAtom1[j12] - 1 );
                pv2 = pBNS->vert + ( v2 = pStereo2Revrs->nBondAtom2[j12] - 1 ); /* djb-rwth: ignoring LLVM warning: variable used */
                /* indicators of -NH: */
                i1 = at2[v1].valence == 1 && at2[v1].num_H == 2 && !at2[v1].endpoint &&
                    pVA[v1].cNumValenceElectrons == 5 && pVA[v1].cPeriodicRowNumber == 1;
                i2 = at2[v2].valence == 1 && at2[v2].num_H == 2 && !at2[v2].endpoint &&
                    pVA[v2].cNumValenceElectrons == 5 && pVA[v2].cPeriodicRowNumber == 1;
                if ((!i1 && !i2) || (i1 && i2)) /* djb-rwth: addressing LLVM warnings */
                {
                    continue;
                }
                /* find the edge between v1 and v2 */
                for (k = 0; k < at2[v1].valence; k++)
                {
                    pe = pBNS->edge + ( e = pv1->iedge[k] ); /* djb-rwth: ignoring LLVM warning: variable used */
                    if (v2 == ( pe->neighbor12 ^ v1 ))
                        break; /* the edge has been found */
                }
                if (k >= at2[v1].valence) /* djb-rwth: addressing LLVM warning */
                {
                    ret = RI_ERR_SYNTAX;
                    goto exit_function;
                }
                if (pe->flow != 1)
                {
                    continue; /* already charged */
                }
                aN = i1 ? v1 : v2; /* -NH atom */
                /* djb-rwth: removing redundant code */
                if (0 > ( ePlus = pVA[aN].nCPlusGroupEdge - 1 ) ||
                    ( pePlus = pBNS->edge + ePlus )->flow ||  /* must be (+) charged */
                     pePlus->forbidden)
                {
                    continue;
                }
                /* Move (+) from =NH2(+) to some other -N<
                */
                /* Fix all charges except (+) on -N< */
                if (!AllChargeEdges.num_edges)
                {
                    for (j = 0; j < pStruct->num_atoms; j++)
                    {
                        if (( k = pVA[j].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[k].forbidden)
                        {
                            if ((ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                        if (( k = pVA[j].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[k].forbidden)
                        {
                            if (pVA[j].cNumValenceElectrons == 5 && pVA[j].cPeriodicRowNumber == 1 &&
                                 !at2[j].num_H && at2[j].valence == 3 &&
                                 !( at2[j].endpoint || (pStruct->endpoint && pStruct->endpoint[j]) )) /* djb-rwth: addressing LLVM warning */
                            {
                                ; /* do not fix -N< or =N(+)< */
                            }
                            else
                            {
                                                 /* all others */
                                if ((ret = AddToEdgeList( TautMinusEdges + 0, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_function;
                                }
                            }
                            /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
                            if (pVA[j].cNumValenceElectrons == 5 && !pVA[j].cMetal && /* N, P, As */
                                 NO_VERTEX != ( k = GetChargeFlowerUpperEdge( pBNS, pVA, k ) ))
                            {
                                if (!pBNS->edge[j].forbidden && pBNS->edge[k].flow)
                                {
                                    if ((ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                    {
                                        goto exit_function;
                                    }
                                }
                            }
                        }
                    }
                }
                /* pePlus edge is already fixed; unfix it */
                /* To decrement (+) on =NH2(+) decrement its double bond order */
                /* djb-rwth: removing redundant code */
                if (!pe->flow)
                    continue;
                pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );
    
                delta = 1;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2 * delta;
    
                pe->forbidden |= forbidden_edge_mask;
                pePlus->forbidden &= ~forbidden_edge_mask;
    
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
    
                if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                    (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 0) /* djb-rwth: addressing LLVM warnings */
                {
                    /* (+)charge was just moved, no change in number of charges */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if (ret > 0)
                    {
                        ( *pnNumRunBNS )++;
                        cur_success++; /* 01 */
                    }
                }
                else
                {
                    pe->flow += delta; /* roll back */
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2 * delta;
                }
                pe->forbidden &= ~forbidden_edge_mask;
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );/* fix aC edges */
            }
    /*exit_case_06:*/
            if (cur_success)
            {
                /* djb-rwth: removing redundant code */
                /* recalculate InChI from the structure */
                if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep ) ))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                /*
                if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                    goto exit_function;
                }
                */
                cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr2, &err );
                if (cmpInChI & IDIF_PROBLEM)
                {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if (err)
                {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
                /* djb-rwth: removing redundant code */
                memset( icr2, 0, sizeof( *icr2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                if (iRevrInChI || iOrigInChI)
                {
                    /* additional mobile-H compare in case of Fixed-H */
                    cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
                    if (cmpInChI & IDIF_PROBLEM)
                    {
                        ret = RI_ERR_PROGR; /* severe restore problem */
                        goto exit_function;
                    }
                    if (err)
                    {
                        ret = RI_ERR_ALLOC;
                        goto exit_function;
                    }
                }
                pStereoRevrs = ( pStruct->pOneINChI[0]->StereoIsotopic &&
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                           pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters )
                    ? pStruct->pOneINChI[0]->StereoIsotopic
                    : pStruct->pOneINChI[0]->Stereo; /* djb-rwth: ignoring LLVM warning: variable used */
    
    
                pStereo2Revrs = ( pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                           !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted )
                    ? NULL
                    : ( pStruct->pOneINChI[1]->StereoIsotopic &&
                        pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                        pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters ) ?
                        pStruct->pOneINChI[1]->StereoIsotopic :
                        pStruct->pOneINChI[1]->Stereo; /* djb-rwth: ignoring LLVM warning: variable used */
            }
        }
    
    
    exit_function:
        SetForbiddenEdgeMask( pBNS, &FixedStereoEdges, forbidden_stereo_edge_mask );
        AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );
        AllocEdgeList( &CurrEdges, EDGE_LIST_FREE );
        AllocEdgeList( &NFlowerEdges, EDGE_LIST_FREE );
        AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_FREE );
        AllocEdgeList( &FixedStereoEdges, EDGE_LIST_FREE );
        AllocEdgeList( &AllRadList, EDGE_LIST_FREE ); /* eliminate memory leak */
        AllocEdgeList( TautMinusEdges + 0, EDGE_LIST_FREE );
        AllocEdgeList( TautMinusEdges + 1, EDGE_LIST_FREE );
    
        return ret;
    }
    */
    // END INCHI C FUNCTION: FixRestoredStructureStereo
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FixRestoredStructureStereo
    // INCHI✔❌: READ_INCHI_STRING=1; COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux.
    // INCHI✔❌: INCHI_HEAPCHK is empty in the selected non-CHECK_WIN32_VC_HEAP build.
    // END INCHI ACTIVE MACRO CONFIGURATION: FixRestoredStructureStereo
    let mut ret = 0_i32;
    let mut all_charge_edges = EDGE_LIST::default();
    let mut current_edges = EDGE_LIST::default();
    let mut nitrogen_flower_edges = EDGE_LIST::default();
    let mut other_nitrogen_flower_edges = EDGE_LIST::default();
    let mut fixed_stereo_edges = EDGE_LIST::default();
    let mut all_radicals = EDGE_LIST::default();
    let mut taut_minus_edges = [EDGE_LIST::default(), EDGE_LIST::default()];
    
    let _ = AllocEdgeList(heap, &mut all_charge_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut current_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut nitrogen_flower_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut other_nitrogen_flower_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut fixed_stereo_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut all_radicals, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut taut_minus_edges[0], EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut taut_minus_edges[1], EDGE_LIST_CLEAR)?;
    
    let execution = (|| -> Result<(), SourceHeapError> {
        macro_rules! first_model {
            ($pointer:expr) => {{
                heap.slice($pointer.as_const())?
                    .first()
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            }};
        }
        macro_rules! selected_stereo {
            ($inchi:expr) => {{
                let inchi = first_model!($inchi);
                if !inchi.StereoIsotopic.is_null() {
                    let isotopic = first_model!(inchi.StereoIsotopic);
                    if isotopic
                        .nNumberOfStereoBonds
                        .wrapping_add(isotopic.nNumberOfStereoCenters)
                        != 0
                    {
                        inchi.StereoIsotopic
                    } else {
                        inchi.Stereo
                    }
                } else {
                    inchi.Stereo
                }
            }};
        }
        macro_rules! compare_pair {
            ($reversed_index:expr, $input_index:expr, $output:expr) => {{
                let reversed = first_model!(pStruct.pOneINChI[$reversed_index]);
                let input = first_model!(pInChI[$input_index]);
                let auxiliary = if pStruct.pOneINChI_Aux[$reversed_index].is_null() {
                    None
                } else {
                    Some(first_model!(pStruct.pOneINChI_Aux[$reversed_index]))
                };
                let mut comparison_error = 0_i32;
                let flags = CompareReversedINChI2(
                    heap,
                    Some(&reversed),
                    Some(&input),
                    auxiliary.as_ref(),
                    None,
                    $output,
                    &mut comparison_error,
                )?;
                (flags, comparison_error)
            }};
        }
        macro_rules! refresh_reversed_stereo {
            ($p_stereo_reversed:ident, $p_stereo2_reversed:ident) => {{
                $p_stereo_reversed = selected_stereo!(pStruct.pOneINChI[0]);
                $p_stereo2_reversed = if i32::from(pStruct.bMobileH) == TAUT_YES as i32
                    || pStruct.pOneINChI[1].is_null()
                {
                    SourceMutPointer::null()
                } else {
                    let reversed_mobile = first_model!(pStruct.pOneINChI[1]);
                    if reversed_mobile.nNumberOfAtoms == 0 || reversed_mobile.bDeleted != 0 {
                        SourceMutPointer::null()
                    } else {
                        selected_stereo!(pStruct.pOneINChI[1])
                    }
                };
            }};
        }
        macro_rules! rebuild_and_compare {
            ($primary_output:expr, $original_index:ident, $reversed_index:ident,
             $p_stereo_reversed:ident, $p_stereo2_reversed:ident, $refresh:expr) => {{
                ret = MakeOneInChIOutOfStrFromINChI2(
                    heap,
                    pCG,
                    ic,
                    ip,
                    sd,
                    pBNS,
                    pStruct,
                    at,
                    at2,
                    at3,
                    pVA,
                    pTCGroups,
                    ppt_group_info.as_deref_mut(),
                    ppat_norm.as_deref_mut(),
                    ppat_prep.as_deref_mut(),
                    clock_result,
                )?;
                if ret < 0 {
                    return Ok(());
                }
                ret = FillOutExtraFixedHDataRestr(heap, pStruct)?;
                if ret != 0 {
                    return Ok(());
                }
                let compared = compare_pair!(0, 0, $primary_output);
                cmp_inchi = compared.0;
                if cmp_inchi & IDIF_PROBLEM != 0 {
                    ret = RI_ERR_PROGR;
                    return Ok(());
                }
                if compared.1 != 0 {
                    ret = RI_ERR_ALLOC;
                    return Ok(());
                }
                cmp_inchi2 = 0;
                *icr2 = ICR::default();
                if $reversed_index != 0 || $original_index != 0 {
                    let compared = compare_pair!($reversed_index, $original_index, icr2);
                    cmp_inchi2 = compared.0;
                    if cmp_inchi & IDIF_PROBLEM != 0 {
                        ret = RI_ERR_PROGR;
                        return Ok(());
                    }
                    if compared.1 != 0 {
                        ret = RI_ERR_ALLOC;
                        return Ok(());
                    }
                }
                if $refresh {
                    refresh_reversed_stereo!($p_stereo_reversed, $p_stereo2_reversed);
                }
            }};
        }
        macro_rules! edge_snapshot {
            ($edge:expr) => {{
                heap.slice(pBNS.edge.as_const())?
                    .get(usize::try_from($edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            }};
        }
        macro_rules! vertex_snapshot {
            ($vertex:expr) => {{
                heap.slice(pBNS.vert.as_const())?
                    .get(usize::try_from($vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            }};
        }
        macro_rules! atom_snapshot {
            ($atom:expr) => {{
                heap.slice(at2.as_const())?
                    .get(usize::try_from($atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            }};
        }
        macro_rules! incident_edge {
            ($vertex:expr, $order:expr) => {{
                let vertex = vertex_snapshot!($vertex);
                *heap
                    .slice(vertex.iedge.as_const())?
                    .get(usize::try_from($order).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
            }};
        }
        macro_rules! find_bond_edge {
            ($first:expr, $second:expr) => {{
                let atom = atom_snapshot!($first);
                let mut order = 0_i32;
                let mut found = None;
                while order < i32::from(atom.valence) {
                    let edge_number = incident_edge!($first, order);
                    let edge = edge_snapshot!(edge_number);
                    if $second == i32::from(edge.neighbor12) ^ $first {
                        found = Some(edge_number);
                        break;
                    }
                    order = order.wrapping_add(1);
                }
                found
            }};
        }
    
        let p_stereo_input = selected_stereo!(pInChI[0]);
        let p_stereo2_input = if i32::from(pStruct.bMobileH) == TAUT_YES as i32
            || pInChI[1].is_null()
        {
            SourceMutPointer::null()
        } else {
            let input_mobile = first_model!(pInChI[1]);
            if input_mobile.nNumberOfAtoms == 0 || input_mobile.bDeleted != 0 {
                SourceMutPointer::null()
            } else {
                selected_stereo!(pInChI[1])
            }
        };
        let mut p_stereo_reversed = selected_stereo!(pStruct.pOneINChI[0]);
        let mut p_stereo2_reversed = if i32::from(pStruct.bMobileH) == TAUT_YES as i32
            || pStruct.pOneINChI[1].is_null()
        {
            SourceMutPointer::null()
        } else {
            let reversed_mobile = first_model!(pStruct.pOneINChI[1]);
            if reversed_mobile.nNumberOfAtoms == 0 || reversed_mobile.bDeleted != 0 {
                SourceMutPointer::null()
            } else {
                selected_stereo!(pStruct.pOneINChI[1])
            }
        };
    
        let compared = compare_pair!(0, 0, icr);
        cmp_inchi = compared.0;
        if cmp_inchi & IDIF_PROBLEM != 0 {
            ret = RI_ERR_PROGR;
            return Ok(());
        }
        if compared.1 != 0 {
            ret = RI_ERR_ALLOC;
            return Ok(());
        }
        cmp_inchi2 = 0;
    
        let (original_inchi_index, reversed_inchi_index) =
            if i32::from(pStruct.bMobileH) == TAUT_NON as i32 {
                let original = if !pInChI[1].is_null() {
                    let value = first_model!(pInChI[1]);
                    usize::from(value.nNumberOfAtoms != 0 && value.bDeleted == 0)
                } else {
                    0
                };
                let reversed = if !pStruct.pOneINChI[1].is_null() {
                    let value = first_model!(pStruct.pOneINChI[1]);
                    usize::from(value.nNumberOfAtoms != 0 && value.bDeleted == 0)
                } else {
                    0
                };
                (original, reversed)
            } else {
                (0, 0)
            };
    
        *icr2 = ICR::default();
        if reversed_inchi_index != 0 || original_inchi_index != 0 {
            let compared = compare_pair!(reversed_inchi_index, original_inchi_index, icr2);
            cmp_inchi2 = compared.0;
            if cmp_inchi & IDIF_PROBLEM != 0 {
                ret = RI_ERR_PROGR;
                return Ok(());
            }
            if compared.1 != 0 {
                ret = RI_ERR_ALLOC;
                return Ok(());
            }
        }
    
        if cmp_inchi & IDIFF_SB == 0 && cmp_inchi2 & IDIFF_SB == 0 {
            return Ok(());
        }
    
        let mut atom_number = 0_i32;
        while atom_number < pStruct.num_atoms {
            let atom = atom_snapshot!(atom_number);
            let mut order = 0_i32;
            while order < i32::from(atom.valence) {
                let edge_number = incident_edge!(atom_number, order);
                let edge = edge_snapshot!(edge_number);
                if atom_number == i32::from(edge.neighbor1)
                    && i32::from(edge.forbidden) & forbidden_stereo_edge_mask != 0
                {
                    ret = AddToEdgeList(heap, &mut fixed_stereo_edges, edge_number, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                order = order.wrapping_add(1);
            }
            atom_number = atom_number.wrapping_add(1);
        }
    
        let mut current_success = 0_i32;
        let maximum_success = pBNS.tot_st_cap.wrapping_sub(pBNS.tot_st_flow);
        if cmp_inchi & IDIF_SB_MISS != 0
            && (cmp_inchi2 == 0 || cmp_inchi2 & IDIF_SB_MISS != 0)
            && maximum_success > 0
        {
            RemoveForbiddenEdgeMask(heap, pBNS, &fixed_stereo_edges, forbidden_stereo_edge_mask)?;
            let stereo = first_model!(p_stereo_input);
            let bond_atoms1 = heap.slice(stereo.nBondAtom1.as_const())?.to_vec();
            let bond_atoms2 = heap.slice(stereo.nBondAtom2.as_const())?.to_vec();
            let count = usize::try_from(icr.num_sb_in2_only.max(0))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let missing_bonds = icr.sb_in2_only[..count.min(icr.sb_in2_only.len())].to_vec();
            for difference in missing_bonds {
                if current_success >= maximum_success {
                    break;
                }
                let difference = usize::from(difference);
                let vertex1 = i32::from(
                    *bond_atoms1
                        .get(difference)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
                .wrapping_sub(1);
                let vertex2 = i32::from(
                    *bond_atoms2
                        .get(difference)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                )
                .wrapping_sub(1);
                let Some(edge_number) = find_bond_edge!(vertex1, vertex2) else {
                    ret = RI_ERR_SYNTAX;
                    return Ok(());
                };
    
                let vertex1_index = usize::try_from(vertex1)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let vertex2_index = usize::try_from(vertex2)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let edge_index = usize::try_from(edge_number)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let vertices = heap.slice_mut(pBNS.vert)?;
                vertices[vertex1_index].st_edge.cap =
                    vertices[vertex1_index].st_edge.cap.wrapping_sub(1);
                vertices[vertex1_index].st_edge.flow =
                    vertices[vertex1_index].st_edge.flow.wrapping_sub(1);
                vertices[vertex2_index].st_edge.flow =
                    vertices[vertex2_index].st_edge.flow.wrapping_sub(1);
                heap.slice_mut(pBNS.edge)?[edge_index].flow =
                    heap.slice(pBNS.edge.as_const())?[edge_index].flow.wrapping_sub(1);
                let mut path_start = 0_i32;
                let mut path_end = 0_i32;
                let mut path_length = 0_i32;
                let mut delta_h = 0_i32;
                let mut delta_charge = 0_i32;
                let mut visited_atoms = 0_i32;
                ret = RunBnsTestOnce(
                    heap,
                    pBNS,
                    pBD,
                    pVA,
                    &mut path_start,
                    &mut path_end,
                    &mut path_length,
                    &mut delta_h,
                    &mut delta_charge,
                    &mut visited_atoms,
                )?;
                let vertices = heap.slice_mut(pBNS.vert)?;
                vertices[vertex1_index].st_edge.cap =
                    vertices[vertex1_index].st_edge.cap.wrapping_add(1);
                vertices[vertex1_index].st_edge.flow =
                    vertices[vertex1_index].st_edge.flow.wrapping_add(1);
                vertices[vertex2_index].st_edge.flow =
                    vertices[vertex2_index].st_edge.flow.wrapping_add(1);
                heap.slice_mut(pBNS.edge)?[edge_index].flow =
                    heap.slice(pBNS.edge.as_const())?[edge_index].flow.wrapping_add(1);
    
                let mut radical_vertex = if ret == 1
                    && delta_charge == 0
                    && (vertex2 == path_start || vertex2 == path_end)
                {
                    if vertex2 == path_start {
                        path_end
                    } else {
                        path_start
                    }
                } else {
                    NO_VERTEX
                };
                if radical_vertex == NO_VERTEX {
                    let vertices = heap.slice_mut(pBNS.vert)?;
                    vertices[vertex2_index].st_edge.cap =
                        vertices[vertex2_index].st_edge.cap.wrapping_sub(1);
                    vertices[vertex2_index].st_edge.flow =
                        vertices[vertex2_index].st_edge.flow.wrapping_sub(1);
                    vertices[vertex1_index].st_edge.flow =
                        vertices[vertex1_index].st_edge.flow.wrapping_sub(1);
                    heap.slice_mut(pBNS.edge)?[edge_index].flow =
                        heap.slice(pBNS.edge.as_const())?[edge_index].flow.wrapping_sub(1);
                    ret = RunBnsTestOnce(
                        heap,
                        pBNS,
                        pBD,
                        pVA,
                        &mut path_start,
                        &mut path_end,
                        &mut path_length,
                        &mut delta_h,
                        &mut delta_charge,
                        &mut visited_atoms,
                    )?;
                    let vertices = heap.slice_mut(pBNS.vert)?;
                    vertices[vertex2_index].st_edge.cap =
                        vertices[vertex2_index].st_edge.cap.wrapping_add(1);
                    vertices[vertex2_index].st_edge.flow =
                        vertices[vertex2_index].st_edge.flow.wrapping_add(1);
                    vertices[vertex1_index].st_edge.flow =
                        vertices[vertex1_index].st_edge.flow.wrapping_add(1);
                    heap.slice_mut(pBNS.edge)?[edge_index].flow =
                        heap.slice(pBNS.edge.as_const())?[edge_index].flow.wrapping_add(1);
                    if ret == 1
                        && delta_charge == 0
                        && (vertex1 == path_start || vertex1 == path_end)
                    {
                        radical_vertex = if vertex1 == path_start {
                            path_end
                        } else {
                            path_start
                        };
                    }
                }
                if radical_vertex == NO_VERTEX {
                    continue;
                }
    
                let radical_index = usize::try_from(radical_vertex)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let radical_atom = atom_snapshot!(radical_vertex);
                let radical_valence_electrons = pVA
                    .get(radical_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .cNumValenceElectrons;
                let radical = vertex_snapshot!(radical_vertex);
                let mut radical_edge_number = NO_VERTEX;
                let mut nitrogen_vertex = NO_VERTEX;
                let mut is_nitrogen_oxygen = false;
                if radical_valence_electrons == 6 && radical_atom.valence == 1 {
                    radical_edge_number = *heap
                        .slice(radical.iedge.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let radical_edge = edge_snapshot!(radical_edge_number);
                    if radical_edge.flow == 0 {
                        nitrogen_vertex = i32::from(radical_edge.neighbor12) ^ radical_vertex;
                        let nitrogen_index = usize::try_from(nitrogen_vertex)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        if pVA
                            .get(nitrogen_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .cNumValenceElectrons
                            == 5
                            && atom_snapshot!(nitrogen_vertex).valence == 2
                        {
                            is_nitrogen_oxygen = true;
                        }
                    }
                }
                if is_nitrogen_oxygen {
                    let nitrogen_index = usize::try_from(nitrogen_vertex)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let radical_edge_index = usize::try_from(radical_edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    heap.slice_mut(pBNS.edge)?[radical_edge_index].flow =
                        heap.slice(pBNS.edge.as_const())?[radical_edge_index]
                            .flow
                            .wrapping_add(1);
                    let vertices = heap.slice_mut(pBNS.vert)?;
                    vertices[radical_index].st_edge.flow =
                        vertices[radical_index].st_edge.flow.wrapping_add(1);
                    vertices[nitrogen_index].st_edge.cap =
                        vertices[nitrogen_index].st_edge.cap.wrapping_add(1);
                    vertices[nitrogen_index].st_edge.flow =
                        vertices[nitrogen_index].st_edge.flow.wrapping_add(1);
                    pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_add(1);
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                    current_success = current_success.wrapping_add(1);
                } else {
                    let delta = radical.st_edge.cap.wrapping_sub(radical.st_edge.flow);
                    heap.slice_mut(pBNS.vert)?[radical_index].st_edge.cap = radical
                        .st_edge
                        .cap
                        .wrapping_sub(delta);
                    pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_sub(delta);
                }
            }
            SetForbiddenEdgeMask(heap, pBNS, &fixed_stereo_edges, forbidden_stereo_edge_mask)?;
            if current_success != 0 {
                rebuild_and_compare!(
                    icr,
                    original_inchi_index,
                    reversed_inchi_index,
                    p_stereo_reversed,
                    p_stereo2_reversed,
                    true
                );
            }
        }
    
        if cmp_inchi & IDIF_SB_MISS == 0
            && cmp_inchi2 & IDIF_SB_MISS != 0
            && icr2.num_sb_in2_only != 0
            && pBNS.tot_st_cap.wrapping_sub(pBNS.tot_st_flow) > 0
        {
            let previous_icr = icr.clone();
            let previous_icr2 = icr2.clone();
            all_radicals.num_edges = 0;
            let mut atom_number = 0_i32;
            while atom_number < pStruct.num_atoms {
                let vertex = vertex_snapshot!(atom_number);
                if vertex.st_edge.cap.wrapping_sub(vertex.st_edge.flow) == 1 {
                    ret = AddToEdgeList(heap, &mut all_radicals, atom_number, INC_ADD_EDGE)?;
                    if ret != 0 {
                        return Ok(());
                    }
                }
                atom_number = atom_number.wrapping_add(1);
            }
            let radical_count = usize::try_from(all_radicals.num_edges.max(0))
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let radical_vertices = if radical_count == 0 {
                Vec::new()
            } else {
                heap.slice(all_radicals.pnEdges.as_const())?[..radical_count].to_vec()
            };
            for vertex in &radical_vertices {
                let index = usize::try_from(*vertex)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                heap.slice_mut(pBNS.vert)?[index].st_edge.cap =
                    heap.slice(pBNS.vert.as_const())?[index]
                        .st_edge
                        .cap
                        .wrapping_sub(1);
                pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_sub(1);
            }
            rebuild_and_compare!(
                icr,
                original_inchi_index,
                reversed_inchi_index,
                p_stereo_reversed,
                p_stereo2_reversed,
                false
            );
            let ret_c = CompareIcr(icr, &previous_icr, None, None, IDIFF_CONSTIT);
            let ret_s = CompareIcr(icr, &previous_icr, None, None, IDIFF_STEREO);
            let ret2_c = CompareIcr(icr2, &previous_icr2, None, None, IDIFF_CONSTIT);
            let ret2_s = CompareIcr(icr2, &previous_icr2, None, None, IDIFF_STEREO);
            if !(ret_c <= 0 && ret_s <= 0 && ret2_c <= 0 && ret2_s < 0) {
                for vertex in radical_vertices {
                    let index = usize::try_from(vertex)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    heap.slice_mut(pBNS.vert)?[index].st_edge.cap =
                        heap.slice(pBNS.vert.as_const())?[index]
                            .st_edge
                            .cap
                            .wrapping_add(1);
                    pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_add(1);
                }
                rebuild_and_compare!(
                    icr,
                    original_inchi_index,
                    reversed_inchi_index,
                    p_stereo_reversed,
                    p_stereo2_reversed,
                    true
                );
            }
        }
    
    current_success = 0;
    if i32::from(pStruct.bMobileH) == TAUT_NON as i32
        && cmp_inchi & IDIF_SB_EXTRA_UNDF != 0
        && !pStruct.endpoint.is_null()
    {
        taut_minus_edges[0].num_edges = 0;
        taut_minus_edges[1].num_edges = 0;
        all_charge_edges.num_edges = 0;
        let stereo_reversed = first_model!(p_stereo_reversed);
        let reversed_bond1 = heap.slice(stereo_reversed.nBondAtom1.as_const())?.to_vec();
        let reversed_bond2 = heap.slice(stereo_reversed.nBondAtom2.as_const())?.to_vec();
        let count = usize::try_from(icr.num_sb_undef_in1_only.max(0))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let differences = icr.sb_undef_in1_only[..count.min(icr.sb_undef_in1_only.len())].to_vec();
        for difference in differences {
            let difference = usize::from(difference);
            let vertex1 = i32::from(
                *reversed_bond1
                    .get(difference)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
            .wrapping_sub(1);
            let vertex2 = i32::from(
                *reversed_bond2
                    .get(difference)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
            .wrapping_sub(1);
    
            if !p_stereo2_reversed.is_null() && icr2.num_sb_undef_in1_only != 0 {
                let stereo2_reversed = first_model!(p_stereo2_reversed);
                let mobile_bond1 = heap.slice(stereo2_reversed.nBondAtom1.as_const())?;
                let mobile_bond2 = heap.slice(stereo2_reversed.nBondAtom2.as_const())?;
                let mobile_count = usize::try_from(icr2.num_sb_undef_in1_only.max(0))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let mut order = 0_usize;
                while order < mobile_count {
                    let mobile_difference = usize::from(
                        *icr2
                            .sb_undef_in1_only
                            .get(order)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    if vertex1
                        == i32::from(
                            *mobile_bond1
                                .get(mobile_difference)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        )
                        && vertex2
                            == i32::from(
                                *mobile_bond2
                                    .get(mobile_difference)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            )
                    {
                        break;
                    }
                    order = order.wrapping_add(1);
                }
                if i32::try_from(order).map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    < icr.num_sb_in1_only
                {
                    continue;
                }
            }
    
            if !p_stereo2_input.is_null() {
                let stereo2_input = first_model!(p_stereo2_input);
                if stereo2_input.nNumberOfStereoBonds != 0 {
                    let input_bond1 = heap.slice(stereo2_input.nBondAtom1.as_const())?;
                    let mut order = 0_i32;
                    while order < stereo2_input.nNumberOfStereoBonds {
                        let index = usize::try_from(order)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        if vertex1
                            == i32::from(
                                *input_bond1
                                    .get(index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            )
                            && vertex2
                                == i32::from(
                                    *input_bond1
                                        .get(index)
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                )
                        {
                            break;
                        }
                        order = order.wrapping_add(1);
                    }
                    if order < stereo2_input.nNumberOfStereoBonds {
                        continue;
                    }
                }
            }
    
            let Some(edge_number) = find_bond_edge!(vertex1, vertex2) else {
                ret = RI_ERR_SYNTAX;
                return Ok(());
            };
            if all_charge_edges.num_edges == 0
                && taut_minus_edges[0].num_edges == 0
                && taut_minus_edges[1].num_edges == 0
            {
                let endpoints = heap.slice(pStruct.endpoint.as_const())?.to_vec();
                let mut atom_number = 0_i32;
                while atom_number < pStruct.num_atoms {
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let valence = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    if minus_edge >= 0 && edge_snapshot!(minus_edge).forbidden == 0 {
                        let target = if *endpoints
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            == 0
                        {
                            &mut all_charge_edges
                        } else if valence.cNumValenceElectrons == 6 {
                            &mut taut_minus_edges[0]
                        } else {
                            &mut taut_minus_edges[1]
                        };
                        ret = AddToEdgeList(heap, target, minus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                    let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                    if plus_edge >= 0 && edge_snapshot!(plus_edge).forbidden == 0 {
                        ret = AddToEdgeList(heap, &mut all_charge_edges, plus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                        if valence.cNumValenceElectrons == 5 && valence.cMetal == 0 {
                            let upper = GetChargeFlowerUpperEdge(heap, pBNS, pVA, plus_edge)?;
                            if upper != NO_VERTEX
                                && edge_snapshot!(atom_number).forbidden == 0
                                && edge_snapshot!(upper).flow != 0
                            {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut all_charge_edges,
                                    upper,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                        }
                    }
                    atom_number = atom_number.wrapping_add(1);
                }
            }
            if edge_snapshot!(edge_number).flow == 0 {
                continue;
            }
            SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            let mut mode = 1_i32;
            let mut done = false;
            while mode < 4 && !done {
                let mut number_endpoints = taut_minus_edges[0]
                    .num_edges
                    .wrapping_add(taut_minus_edges[1].num_edges);
                if mode == 2 {
                    SetForbiddenEdgeMask(heap, pBNS, &taut_minus_edges[0], forbidden_edge_mask)?;
                    number_endpoints =
                        number_endpoints.wrapping_sub(taut_minus_edges[0].num_edges);
                }
                if mode == 1 {
                    SetForbiddenEdgeMask(heap, pBNS, &taut_minus_edges[1], forbidden_edge_mask)?;
                    number_endpoints =
                        number_endpoints.wrapping_sub(taut_minus_edges[1].num_edges);
                }
                if number_endpoints >= 2 {
                    let edge = edge_snapshot!(edge_number);
                    let first = i32::from(edge.neighbor1);
                    let second = i32::from(edge.neighbor12) ^ first;
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let first_index = usize::try_from(first)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let second_index = usize::try_from(second)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    heap.slice_mut(pBNS.edge)?[edge_index].forbidden =
                        (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                    heap.slice_mut(pBNS.edge)?[edge_index].flow =
                        edge.flow.wrapping_sub(1);
                    let vertices = heap.slice_mut(pBNS.vert)?;
                    vertices[first_index].st_edge.flow =
                        vertices[first_index].st_edge.flow.wrapping_sub(1);
                    vertices[second_index].st_edge.flow =
                        vertices[second_index].st_edge.flow.wrapping_sub(1);
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);
                    let mut path_start = 0_i32;
                    let mut path_end = 0_i32;
                    let mut path_length = 0_i32;
                    let mut delta_h = 0_i32;
                    let mut delta_charge = 0_i32;
                    let mut visited_atoms = 0_i32;
                    ret = RunBnsTestOnce(
                        heap,
                        pBNS,
                        pBD,
                        pVA,
                        &mut path_start,
                        &mut path_end,
                        &mut path_length,
                        &mut delta_h,
                        &mut delta_charge,
                        &mut visited_atoms,
                    )?;
                    if ret == 1
                        && ((path_end == first && path_start == second)
                            || (path_end == second && path_start == first))
                        && delta_charge == 0
                    {
                        ret = RunBnsRestoreOnce(
                            heap,
                            pBNS,
                            pBD,
                            pVA,
                            pTCGroups,
                            clock_result,
                        )?;
                        if ret > 0 {
                            *pn_num_run_bns = pn_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                            done = true;
                        }
                    } else {
                        heap.slice_mut(pBNS.edge)?[edge_index].forbidden =
                            (i32::from(heap.slice(pBNS.edge.as_const())?[edge_index].forbidden)
                                & !forbidden_edge_mask) as i8;
                        let current_flow = heap.slice(pBNS.edge.as_const())?[edge_index].flow;
                        heap.slice_mut(pBNS.edge)?[edge_index].flow = current_flow.wrapping_add(1);
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        vertices[first_index].st_edge.flow =
                            vertices[first_index].st_edge.flow.wrapping_add(1);
                        vertices[second_index].st_edge.flow =
                            vertices[second_index].st_edge.flow.wrapping_add(1);
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                    }
                }
                if mode == 2 {
                    RemoveForbiddenEdgeMask(
                        heap,
                        pBNS,
                        &taut_minus_edges[0],
                        forbidden_edge_mask,
                    )?;
                }
                if mode == 1 {
                    RemoveForbiddenEdgeMask(
                        heap,
                        pBNS,
                        &taut_minus_edges[1],
                        forbidden_edge_mask,
                    )?;
                }
                mode = mode.wrapping_add(1);
            }
            RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
        }
        if current_success != 0 {
            rebuild_and_compare!(
                icr,
                original_inchi_index,
                reversed_inchi_index,
                p_stereo_reversed,
                p_stereo2_reversed,
                true
            );
        }
    }
    
    current_success = 0;
    if cmp_inchi & IDIF_SB_EXTRA_UNDF != 0 {
        taut_minus_edges[0].num_edges = 0;
        taut_minus_edges[1].num_edges = 0;
        all_charge_edges.num_edges = 0;
        let stereo_reversed = first_model!(p_stereo_reversed);
        let reversed_bond1 = heap.slice(stereo_reversed.nBondAtom1.as_const())?.to_vec();
        let reversed_bond2 = heap.slice(stereo_reversed.nBondAtom2.as_const())?.to_vec();
        let count = usize::try_from(icr.num_sb_undef_in1_only.max(0))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let differences = icr.sb_undef_in1_only[..count.min(icr.sb_undef_in1_only.len())].to_vec();
        for difference in differences {
            let difference = usize::from(difference);
            let vertex1 = i32::from(
                *reversed_bond1
                    .get(difference)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
            .wrapping_sub(1);
            let vertex2 = i32::from(
                *reversed_bond2
                    .get(difference)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
            .wrapping_sub(1);
            let first_incident = incident_edge!(vertex1, 0);
            let mut edge_number = first_incident;
            let atom = atom_snapshot!(vertex1);
            let mut order = 0_i32;
            while order < i32::from(atom.valence) {
                edge_number = incident_edge!(vertex1, order);
                if vertex2 == i32::from(edge_snapshot!(edge_number).neighbor12) ^ vertex1 {
                    break;
                }
                order = order.wrapping_add(1);
            }
            if order == i32::from(atom.valence) {
                ret = RI_ERR_SYNTAX;
                return Ok(());
            }
            if !p_stereo2_reversed.is_null() && icr2.num_sb_undef_in1_only != 0 {
                let stereo2_reversed = first_model!(p_stereo2_reversed);
                let mobile_bond1 = heap.slice(stereo2_reversed.nBondAtom1.as_const())?;
                let mobile_bond2 = heap.slice(stereo2_reversed.nBondAtom2.as_const())?;
                let mobile_count = usize::try_from(icr2.num_sb_undef_in1_only.max(0))
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let mut mobile_order = 0_usize;
                while mobile_order < mobile_count {
                    let mobile_difference = usize::from(
                        *icr2
                            .sb_undef_in1_only
                            .get(mobile_order)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    if vertex1
                        == i32::from(
                            *mobile_bond1
                                .get(mobile_difference)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        )
                        && vertex2
                            == i32::from(
                                *mobile_bond2
                                    .get(mobile_difference)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            )
                    {
                        break;
                    }
                    mobile_order = mobile_order.wrapping_add(1);
                }
                if i32::try_from(mobile_order)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    == icr.num_sb_in1_only
                {
                    continue;
                }
            }
    
            if all_charge_edges.num_edges == 0
                && taut_minus_edges[0].num_edges == 0
                && taut_minus_edges[1].num_edges == 0
            {
                let endpoints = if pStruct.endpoint.is_null() {
                    None
                } else {
                    Some(heap.slice(pStruct.endpoint.as_const())?.to_vec())
                };
                let mut atom_number = 0_i32;
                while atom_number < pStruct.num_atoms {
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let valence = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    if minus_edge >= 0 && edge_snapshot!(minus_edge).forbidden == 0 {
                        ret = AddToEdgeList(
                            heap,
                            &mut all_charge_edges,
                            minus_edge,
                            INC_ADD_EDGE,
                        )?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                    let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                    if plus_edge >= 0 && edge_snapshot!(plus_edge).forbidden == 0 {
                        let atom = atom_snapshot!(atom_number);
                        let may_be_unfixed = atom.num_H == 0
                            && endpoints
                                .as_ref()
                                .and_then(|values| values.get(atom_index))
                                .is_none_or(|endpoint| *endpoint == 0);
                        let target = if (may_be_unfixed && valence.cNumValenceElectrons == 6)
                            || (valence.cNumValenceElectrons == 5
                                && valence.cPeriodicRowNumber > 1)
                        {
                            &mut taut_minus_edges[0]
                        } else if may_be_unfixed
                            && valence.cNumValenceElectrons == 5
                            && valence.cPeriodicRowNumber == 1
                        {
                            &mut taut_minus_edges[1]
                        } else {
                            &mut all_charge_edges
                        };
                        ret = AddToEdgeList(heap, target, plus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            return Ok(());
                        }
                        if valence.cNumValenceElectrons == 5 && valence.cMetal == 0 {
                            let upper = GetChargeFlowerUpperEdge(heap, pBNS, pVA, plus_edge)?;
                            if upper != NO_VERTEX
                                && edge_snapshot!(atom_number).forbidden == 0
                                && edge_snapshot!(upper).flow != 0
                            {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut all_charge_edges,
                                    upper,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                        }
                    }
                    atom_number = atom_number.wrapping_add(1);
                }
            }
            if edge_snapshot!(edge_number).flow == 0 {
                continue;
            }
            SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            let mut mode = 1_i32;
            let mut done = false;
            while mode < 4 && !done {
                let mut number_endpoints = taut_minus_edges[0]
                    .num_edges
                    .wrapping_add(taut_minus_edges[1].num_edges);
                if mode == 2 {
                    SetForbiddenEdgeMask(heap, pBNS, &taut_minus_edges[0], forbidden_edge_mask)?;
                    number_endpoints =
                        number_endpoints.wrapping_sub(taut_minus_edges[0].num_edges);
                }
                if mode == 1 {
                    SetForbiddenEdgeMask(heap, pBNS, &taut_minus_edges[1], forbidden_edge_mask)?;
                    number_endpoints =
                        number_endpoints.wrapping_sub(taut_minus_edges[1].num_edges);
                }
                if number_endpoints >= 2 {
                    let edge = edge_snapshot!(edge_number);
                    let first = i32::from(edge.neighbor1);
                    let second = i32::from(edge.neighbor12) ^ first;
                    let edge_index = usize::try_from(edge_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let first_index = usize::try_from(first)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let second_index = usize::try_from(second)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    heap.slice_mut(pBNS.edge)?[edge_index].forbidden =
                        (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                    heap.slice_mut(pBNS.edge)?[edge_index].flow =
                        edge.flow.wrapping_sub(1);
                    let vertices = heap.slice_mut(pBNS.vert)?;
                    vertices[first_index].st_edge.flow =
                        vertices[first_index].st_edge.flow.wrapping_sub(1);
                    vertices[second_index].st_edge.flow =
                        vertices[second_index].st_edge.flow.wrapping_sub(1);
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);
                    let mut path_start = 0_i32;
                    let mut path_end = 0_i32;
                    let mut path_length = 0_i32;
                    let mut delta_h = 0_i32;
                    let mut delta_charge = 0_i32;
                    let mut visited_atoms = 0_i32;
                    ret = RunBnsTestOnce(
                        heap,
                        pBNS,
                        pBD,
                        pVA,
                        &mut path_start,
                        &mut path_end,
                        &mut path_length,
                        &mut delta_h,
                        &mut delta_charge,
                        &mut visited_atoms,
                    )?;
                    if ret == 1
                        && ((path_end == first && path_start == second)
                            || (path_end == second && path_start == first))
                        && delta_charge == 0
                    {
                        ret = RunBnsRestoreOnce(
                            heap,
                            pBNS,
                            pBD,
                            pVA,
                            pTCGroups,
                            clock_result,
                        )?;
                        if ret > 0 {
                            *pn_num_run_bns = pn_num_run_bns.wrapping_add(1);
                            current_success = current_success.wrapping_add(1);
                            done = true;
                        }
                    } else {
                        heap.slice_mut(pBNS.edge)?[edge_index].forbidden =
                            (i32::from(heap.slice(pBNS.edge.as_const())?[edge_index].forbidden)
                                & !forbidden_edge_mask) as i8;
                        let current_flow = heap.slice(pBNS.edge.as_const())?[edge_index].flow;
                        heap.slice_mut(pBNS.edge)?[edge_index].flow = current_flow.wrapping_add(1);
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        vertices[first_index].st_edge.flow =
                            vertices[first_index].st_edge.flow.wrapping_add(1);
                        vertices[second_index].st_edge.flow =
                            vertices[second_index].st_edge.flow.wrapping_add(1);
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                    }
                }
                if mode == 2 {
                    RemoveForbiddenEdgeMask(
                        heap,
                        pBNS,
                        &taut_minus_edges[0],
                        forbidden_edge_mask,
                    )?;
                }
                if mode == 1 {
                    RemoveForbiddenEdgeMask(
                        heap,
                        pBNS,
                        &taut_minus_edges[1],
                        forbidden_edge_mask,
                    )?;
                }
                mode = mode.wrapping_add(1);
            }
            RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
        }
        if current_success != 0 {
            rebuild_and_compare!(
                icr,
                original_inchi_index,
                reversed_inchi_index,
                p_stereo_reversed,
                p_stereo2_reversed,
                true
            );
        }
    }
    
    current_success = 0;
    if i32::from(pStruct.bMobileH) == TAUT_YES as i32
        && cmp_inchi & IDIF_SB_EXTRA_UNDF != 0
        && pStruct.ti.num_t_groups == 0
    {
        all_charge_edges.num_edges = 0;
        let stereo_reversed = first_model!(p_stereo_reversed);
        let reversed_bond1 = heap.slice(stereo_reversed.nBondAtom1.as_const())?.to_vec();
        let reversed_bond2 = heap.slice(stereo_reversed.nBondAtom2.as_const())?.to_vec();
        let count = usize::try_from(icr.num_sb_undef_in1_only.max(0))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let differences = icr.sb_undef_in1_only[..count.min(icr.sb_undef_in1_only.len())].to_vec();
        for difference in differences {
            let difference = usize::from(difference);
            let vertex1 = i32::from(
                *reversed_bond1
                    .get(difference)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
            .wrapping_sub(1);
            let vertex2 = i32::from(
                *reversed_bond2
                    .get(difference)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
            .wrapping_sub(1);
            let atom1 = atom_snapshot!(vertex1);
            let atom2 = atom_snapshot!(vertex2);
            let index1 = usize::try_from(vertex1)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let index2 = usize::try_from(vertex2)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let first_is_nh = atom1.valence == 1
                && atom1.num_H == 1
                && atom1.endpoint == 0
                && pVA[index1].cNumValenceElectrons == 5
                && pVA[index1].cPeriodicRowNumber == 1;
            let second_is_nh = atom2.valence == 1
                && atom2.num_H == 1
                && atom2.endpoint == 0
                && pVA[index2].cNumValenceElectrons == 5
                && pVA[index2].cPeriodicRowNumber == 1;
            if first_is_nh == second_is_nh {
                continue;
            }
            let Some(edge_number) = find_bond_edge!(vertex1, vertex2) else {
                ret = RI_ERR_SYNTAX;
                return Ok(());
            };
            let edge = edge_snapshot!(edge_number);
            if edge.flow != 1 {
                continue;
            }
            let nitrogen = if first_is_nh { vertex1 } else { vertex2 };
            let carbon = if first_is_nh { vertex2 } else { vertex1 };
            let plus_minus = GetPlusMinusVertex(heap, pBNS, pTCGroups, 1, 0)?;
            if plus_minus == NO_VERTEX {
                break;
            }
    
            let plus_minus_vertex = vertex_snapshot!(plus_minus);
            let mut plus_order = 0_i32;
            while plus_order < i32::from(plus_minus_vertex.num_adj_edges) {
                let first_edge_number = incident_edge!(plus_minus, plus_order);
                let neighbor = i32::from(edge_snapshot!(first_edge_number).neighbor12) ^ plus_minus;
                let neighbor_vertex = vertex_snapshot!(neighbor);
                let mut neighbor_order = 0_i32;
                while neighbor_order < i32::from(neighbor_vertex.num_adj_edges) {
                    let capacity_edge = incident_edge!(neighbor, neighbor_order);
                    let capacity_index = usize::try_from(capacity_edge)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let capacity = heap.slice(pBNS.edge.as_const())?[capacity_index].cap;
                    heap.slice_mut(pBNS.edge)?[capacity_index].cap = capacity.wrapping_add(1);
                    neighbor_order = neighbor_order.wrapping_add(1);
                }
                plus_order = plus_order.wrapping_add(1);
            }
    
            if all_charge_edges.num_edges == 0 {
                let endpoints = if pStruct.endpoint.is_null() {
                    None
                } else {
                    Some(heap.slice(pStruct.endpoint.as_const())?.to_vec())
                };
                let mut atom_number = 0_i32;
                while atom_number < pStruct.num_atoms {
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let valence = &pVA[atom_index];
                    let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    if minus_edge >= 0 && edge_snapshot!(minus_edge).forbidden == 0 {
                        ret = AddToEdgeList(
                            heap,
                            &mut all_charge_edges,
                            minus_edge,
                            INC_ADD_EDGE,
                        )?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                    let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                    if plus_edge >= 0 && edge_snapshot!(plus_edge).forbidden == 0 {
                        let atom = atom_snapshot!(atom_number);
                        let movable_nitrogen = valence.cNumValenceElectrons == 5
                            && valence.cPeriodicRowNumber == 1
                            && atom.num_H == 0
                            && atom.valence == 3
                            && atom.endpoint == 0
                            && endpoints
                                .as_ref()
                                .and_then(|values| values.get(atom_index))
                                .is_none_or(|endpoint| *endpoint == 0);
                        if !movable_nitrogen {
                            ret = AddToEdgeList(
                                heap,
                                &mut taut_minus_edges[0],
                                plus_edge,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                return Ok(());
                            }
                        }
                        if valence.cNumValenceElectrons == 5 && valence.cMetal == 0 {
                            let upper = GetChargeFlowerUpperEdge(heap, pBNS, pVA, plus_edge)?;
                            if upper != NO_VERTEX
                                && edge_snapshot!(atom_number).forbidden == 0
                                && edge_snapshot!(upper).flow != 0
                            {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut all_charge_edges,
                                    upper,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                        }
                    }
                    atom_number = atom_number.wrapping_add(1);
                }
            }
    
            let edge_index = usize::try_from(edge_number)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let nitrogen_index = usize::try_from(nitrogen)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let carbon_index = usize::try_from(carbon)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let plus_minus_index = usize::try_from(plus_minus)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(pBNS.edge)?[edge_index].flow = edge.flow.wrapping_sub(1);
            let vertices = heap.slice_mut(pBNS.vert)?;
            vertices[nitrogen_index].st_edge.flow =
                vertices[nitrogen_index].st_edge.flow.wrapping_sub(1);
            vertices[nitrogen_index].st_edge.cap =
                vertices[nitrogen_index].st_edge.cap.wrapping_sub(1);
            vertices[carbon_index].st_edge.flow =
                vertices[carbon_index].st_edge.flow.wrapping_sub(1);
            vertices[plus_minus_index].st_edge.cap =
                vertices[plus_minus_index].st_edge.cap.wrapping_add(1);
            pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);
            ret = AddToEdgeList(heap, &mut all_charge_edges, edge_number, INC_ADD_EDGE)?;
            if ret != 0 {
                return Ok(());
            }
            SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            let current_cap = heap.slice(pBNS.edge.as_const())?[edge_index].cap;
            heap.slice_mut(pBNS.edge)?[edge_index].cap = current_cap.wrapping_add(1);
            let mut path_start = 0_i32;
            let mut path_end = 0_i32;
            let mut path_length = 0_i32;
            let mut delta_h = 0_i32;
            let mut delta_charge = 0_i32;
            let mut visited_atoms = 0_i32;
            ret = RunBnsTestOnce(
                heap,
                pBNS,
                pBD,
                pVA,
                &mut path_start,
                &mut path_end,
                &mut path_length,
                &mut delta_h,
                &mut delta_charge,
                &mut visited_atoms,
            )?;
            let test_matched = ret == 1
                && ((path_end == plus_minus && path_start == carbon)
                    || (path_end == carbon && path_start == plus_minus))
                && delta_charge == 1;
            if test_matched {
                ret = RunBnsRestoreOnce(
                    heap,
                    pBNS,
                    pBD,
                    pVA,
                    pTCGroups,
                    clock_result,
                )?;
                if ret > 0 {
                    *pn_num_run_bns = pn_num_run_bns.wrapping_add(1);
                    let structure_atom = heap
                        .slice_mut(pStruct.at)?
                        .get_mut(nitrogen_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    structure_atom.num_H = structure_atom.num_H.wrapping_add(1);
                    pTCGroups.total_charge = pTCGroups.total_charge.wrapping_add(1);
                    current_success = current_success.wrapping_add(1);
                }
            }
            if !test_matched {
                let current_flow = heap.slice(pBNS.edge.as_const())?[edge_index].flow;
                heap.slice_mut(pBNS.edge)?[edge_index].flow = current_flow.wrapping_add(1);
                let vertices = heap.slice_mut(pBNS.vert)?;
                vertices[nitrogen_index].st_edge.flow =
                    vertices[nitrogen_index].st_edge.flow.wrapping_add(1);
                vertices[nitrogen_index].st_edge.cap =
                    vertices[nitrogen_index].st_edge.cap.wrapping_add(1);
                vertices[carbon_index].st_edge.flow =
                    vertices[carbon_index].st_edge.flow.wrapping_add(1);
                vertices[plus_minus_index].st_edge.cap =
                    vertices[plus_minus_index].st_edge.cap.wrapping_sub(1);
                pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                all_charge_edges.num_edges = all_charge_edges.num_edges.wrapping_sub(1);
                current_edges.num_edges = 0;
                continue;
            }
            RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
            all_charge_edges.num_edges = all_charge_edges.num_edges.wrapping_sub(1);
            current_edges.num_edges = 0;
        }
        if current_success != 0 {
            rebuild_and_compare!(
                icr,
                original_inchi_index,
                reversed_inchi_index,
                p_stereo_reversed,
                p_stereo2_reversed,
                true
            );
        }
    }
    
    current_success = 0;
    if i32::from(pStruct.bMobileH) == TAUT_NON as i32
        && !p_stereo2_reversed.is_null()
        && cmp_inchi2 & IDIF_SB_EXTRA_UNDF != 0
        && pStruct.ti.num_t_groups == 0
    {
        all_charge_edges.num_edges = 0;
        let stereo2_reversed = first_model!(p_stereo2_reversed);
        let reversed_bond1 = heap.slice(stereo2_reversed.nBondAtom1.as_const())?.to_vec();
        let reversed_bond2 = heap.slice(stereo2_reversed.nBondAtom2.as_const())?.to_vec();
        let count = usize::try_from(icr2.num_sb_undef_in1_only.max(0))
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let differences = icr2.sb_undef_in1_only[..count.min(icr2.sb_undef_in1_only.len())].to_vec();
        for difference in differences {
            let difference = usize::from(difference);
            let vertex1 = i32::from(
                *reversed_bond1
                    .get(difference)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
            .wrapping_sub(1);
            let vertex2 = i32::from(
                *reversed_bond2
                    .get(difference)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            )
            .wrapping_sub(1);
            let atom1 = atom_snapshot!(vertex1);
            let atom2 = atom_snapshot!(vertex2);
            let index1 = usize::try_from(vertex1)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let index2 = usize::try_from(vertex2)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let first_is_nh2 = atom1.valence == 1
                && atom1.num_H == 2
                && atom1.endpoint == 0
                && pVA[index1].cNumValenceElectrons == 5
                && pVA[index1].cPeriodicRowNumber == 1;
            let second_is_nh2 = atom2.valence == 1
                && atom2.num_H == 2
                && atom2.endpoint == 0
                && pVA[index2].cNumValenceElectrons == 5
                && pVA[index2].cPeriodicRowNumber == 1;
            if first_is_nh2 == second_is_nh2 {
                continue;
            }
            let Some(edge_number) = find_bond_edge!(vertex1, vertex2) else {
                ret = RI_ERR_SYNTAX;
                return Ok(());
            };
            let edge = edge_snapshot!(edge_number);
            if edge.flow != 1 {
                continue;
            }
            let nitrogen = if first_is_nh2 { vertex1 } else { vertex2 };
            let nitrogen_index = usize::try_from(nitrogen)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let plus_edge_number = pVA[nitrogen_index].nCPlusGroupEdge.wrapping_sub(1);
            if plus_edge_number < 0 {
                continue;
            }
            let plus_edge = edge_snapshot!(plus_edge_number);
            if plus_edge.flow != 0 || plus_edge.forbidden != 0 {
                continue;
            }
    
            if all_charge_edges.num_edges == 0 {
                let endpoints = if pStruct.endpoint.is_null() {
                    None
                } else {
                    Some(heap.slice(pStruct.endpoint.as_const())?.to_vec())
                };
                let mut atom_number = 0_i32;
                while atom_number < pStruct.num_atoms {
                    let atom_index = usize::try_from(atom_number)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let valence = &pVA[atom_index];
                    let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    if minus_edge >= 0 && edge_snapshot!(minus_edge).forbidden == 0 {
                        ret = AddToEdgeList(
                            heap,
                            &mut all_charge_edges,
                            minus_edge,
                            INC_ADD_EDGE,
                        )?;
                        if ret != 0 {
                            return Ok(());
                        }
                    }
                    let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                    if plus_edge >= 0 && edge_snapshot!(plus_edge).forbidden == 0 {
                        let atom = atom_snapshot!(atom_number);
                        let movable_nitrogen = valence.cNumValenceElectrons == 5
                            && valence.cPeriodicRowNumber == 1
                            && atom.num_H == 0
                            && atom.valence == 3
                            && atom.endpoint == 0
                            && endpoints
                                .as_ref()
                                .and_then(|values| values.get(atom_index))
                                .is_none_or(|endpoint| *endpoint == 0);
                        if !movable_nitrogen {
                            ret = AddToEdgeList(
                                heap,
                                &mut taut_minus_edges[0],
                                plus_edge,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                return Ok(());
                            }
                        }
                        if valence.cNumValenceElectrons == 5 && valence.cMetal == 0 {
                            let upper = GetChargeFlowerUpperEdge(heap, pBNS, pVA, plus_edge)?;
                            if upper != NO_VERTEX
                                && edge_snapshot!(atom_number).forbidden == 0
                                && edge_snapshot!(upper).flow != 0
                            {
                                ret = AddToEdgeList(
                                    heap,
                                    &mut all_charge_edges,
                                    upper,
                                    INC_ADD_EDGE,
                                )?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                        }
                    }
                    atom_number = atom_number.wrapping_add(1);
                }
            }
            if edge_snapshot!(edge_number).flow == 0 {
                continue;
            }
            let edge = edge_snapshot!(edge_number);
            let first = i32::from(edge.neighbor1);
            let second = i32::from(edge.neighbor12) ^ first;
            let edge_index = usize::try_from(edge_number)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let first_index = usize::try_from(first)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let second_index = usize::try_from(second)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let plus_edge_index = usize::try_from(plus_edge_number)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(pBNS.edge)?[edge_index].flow = edge.flow.wrapping_sub(1);
            let vertices = heap.slice_mut(pBNS.vert)?;
            vertices[first_index].st_edge.flow =
                vertices[first_index].st_edge.flow.wrapping_sub(1);
            vertices[second_index].st_edge.flow =
                vertices[second_index].st_edge.flow.wrapping_sub(1);
            pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);
            heap.slice_mut(pBNS.edge)?[edge_index].forbidden =
                (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
            heap.slice_mut(pBNS.edge)?[plus_edge_index].forbidden =
                (i32::from(plus_edge.forbidden) & !forbidden_edge_mask) as i8;
            let mut path_start = 0_i32;
            let mut path_end = 0_i32;
            let mut path_length = 0_i32;
            let mut delta_h = 0_i32;
            let mut delta_charge = 0_i32;
            let mut visited_atoms = 0_i32;
            ret = RunBnsTestOnce(
                heap,
                pBNS,
                pBD,
                pVA,
                &mut path_start,
                &mut path_end,
                &mut path_length,
                &mut delta_h,
                &mut delta_charge,
                &mut visited_atoms,
            )?;
            if ret == 1
                && ((path_end == first && path_start == second)
                    || (path_end == second && path_start == first))
                && delta_charge == 0
            {
                ret = RunBnsRestoreOnce(
                    heap,
                    pBNS,
                    pBD,
                    pVA,
                    pTCGroups,
                    clock_result,
                )?;
                if ret > 0 {
                    *pn_num_run_bns = pn_num_run_bns.wrapping_add(1);
                    current_success = current_success.wrapping_add(1);
                }
            } else {
                let current_flow = heap.slice(pBNS.edge.as_const())?[edge_index].flow;
                heap.slice_mut(pBNS.edge)?[edge_index].flow = current_flow.wrapping_add(1);
                let vertices = heap.slice_mut(pBNS.vert)?;
                vertices[first_index].st_edge.flow =
                    vertices[first_index].st_edge.flow.wrapping_add(1);
                vertices[second_index].st_edge.flow =
                    vertices[second_index].st_edge.flow.wrapping_add(1);
                pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
            }
            let current_forbidden = heap.slice(pBNS.edge.as_const())?[edge_index].forbidden;
            heap.slice_mut(pBNS.edge)?[edge_index].forbidden =
                (i32::from(current_forbidden) & !forbidden_edge_mask) as i8;
            RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
        }
        if current_success != 0 {
            rebuild_and_compare!(
                icr2,
                original_inchi_index,
                reversed_inchi_index,
                p_stereo_reversed,
                p_stereo2_reversed,
                true
            );
        }
    }
        Ok(())
    })();


    let set_result =
        SetForbiddenEdgeMask(heap, pBNS, &fixed_stereo_edges, forbidden_stereo_edge_mask);
    let free_all_charge = AllocEdgeList(heap, &mut all_charge_edges, EDGE_LIST_FREE);
    let free_current = AllocEdgeList(heap, &mut current_edges, EDGE_LIST_FREE);
    let free_nitrogen = AllocEdgeList(heap, &mut nitrogen_flower_edges, EDGE_LIST_FREE);
    let free_other_nitrogen =
        AllocEdgeList(heap, &mut other_nitrogen_flower_edges, EDGE_LIST_FREE);
    let free_stereo = AllocEdgeList(heap, &mut fixed_stereo_edges, EDGE_LIST_FREE);
    let free_radicals = AllocEdgeList(heap, &mut all_radicals, EDGE_LIST_FREE);
    let free_taut_oxygen = AllocEdgeList(heap, &mut taut_minus_edges[0], EDGE_LIST_FREE);
    let free_taut_nitrogen = AllocEdgeList(heap, &mut taut_minus_edges[1], EDGE_LIST_FREE);

    execution?;
    set_result?;
    let _ = free_all_charge?;
    let _ = free_current?;
    let _ = free_nitrogen?;
    let _ = free_other_nitrogen?;
    let _ = free_stereo?;
    let _ = free_radicals?;
    let _ = free_taut_oxygen?;
    let _ = free_taut_nitrogen?;
    Ok(ret)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{AB_PARITY_UNDF, BNS_EDGE, BNS_ST_EDGE, BNS_VERTEX, INChI_Stereo};

    const FIXED_MASK: i32 = 8;
    const FORBIDDEN_MASK: i32 = 4;

    fn inchi(heap: &mut SourceHeap, parity: Option<i8>, error: i32) -> SourceMutPointer<INChI> {
        let stereo = if let Some(parity) = parity {
            let bond_atom1 = heap.allocate_model_storage(vec![1_u16]).unwrap();
            let bond_atom2 = heap.allocate_model_storage(vec![2_u16]).unwrap();
            let bond_parity = heap.allocate_model_storage(vec![parity]).unwrap();
            INChI_Stereo {
                nNumberOfStereoBonds: 1,
                nBondAtom1: bond_atom1,
                nBondAtom2: bond_atom2,
                b_parity: bond_parity,
                ..INChI_Stereo::default()
            }
        } else {
            INChI_Stereo::default()
        };
        let stereo = heap.allocate_model_storage(vec![stereo]).unwrap();
        let formula = heap
            .allocate_model_storage(vec![b'C' as i8, b'2' as i8, 0])
            .unwrap();
        let atoms = heap.allocate_model_storage(vec![6_u8, 6]).unwrap();
        let hydrogens = heap.allocate_model_storage(vec![0_i8, 0]).unwrap();
        heap.allocate_model_storage(vec![INChI {
            nErrorCode: error,
            nNumberOfAtoms: 2,
            szHillFormula: formula,
            nAtom: atoms,
            nNum_H: hydrogens,
            Stereo: stereo,
            ..INChI::default()
        }])
        .unwrap()
    }

    fn atoms(heap: &mut SourceHeap, values: [inp_ATOM; 2]) -> SourceMutPointer<inp_ATOM> {
        heap.allocate_model_storage(values.to_vec()).unwrap()
    }

    fn bond_network(
        heap: &mut SourceHeap,
        edge_flow: i32,
        forbidden: i8,
        total_capacity: i32,
    ) -> BN_STRUCT {
        let incident0 = heap.allocate_model_storage(vec![0_i32]).unwrap();
        let incident1 = heap.allocate_model_storage(vec![0_i32]).unwrap();
        let vertices = heap
            .allocate_model_storage(vec![
                BNS_VERTEX {
                    st_edge: BNS_ST_EDGE {
                        cap: total_capacity,
                        ..BNS_ST_EDGE::default()
                    },
                    num_adj_edges: 1,
                    max_adj_edges: 1,
                    iedge: incident0,
                    ..BNS_VERTEX::default()
                },
                BNS_VERTEX {
                    num_adj_edges: 1,
                    max_adj_edges: 1,
                    iedge: incident1,
                    ..BNS_VERTEX::default()
                },
            ])
            .unwrap();
        let edges = heap
            .allocate_model_storage(vec![BNS_EDGE {
                neighbor1: 1,
                neighbor12: 1,
                neigh_ord: [0, 0],
                cap: 1,
                flow: edge_flow,
                forbidden,
                ..BNS_EDGE::default()
            }])
            .unwrap();
        BN_STRUCT {
            num_atoms: 2,
            num_vertices: 2,
            num_bonds: 1,
            num_edges: 1,
            tot_st_cap: total_capacity,
            vert: vertices,
            edge: edges,
            ..BN_STRUCT::default()
        }
    }

    fn invoke(
        heap: &mut SourceHeap,
        structure: &mut StrFromINChI,
        network: &mut BN_STRUCT,
        at: SourceMutPointer<inp_ATOM>,
        valence: &mut [VAL_AT],
        input: [SourceMutPointer<INChI>; 2],
    ) -> (i32, ICR, ICR, i32) {
        let mut canonical = CANON_GLOBALS::default();
        let mut data = BN_DATA::default();
        let mut primary = ICR::default();
        let mut secondary = ICR::default();
        let mut groups = ALL_TC_GROUPS::default();
        let mut runs = 0;
        let mut delta = 0;
        let result = FixRestoredStructureStereo(
            heap,
            &mut canonical,
            SourceMutPointer::null(),
            INCHI_MODE::MAX,
            &mut primary,
            INCHI_MODE::MAX,
            &mut secondary,
            &INPUT_PARMS::default(),
            &STRUCT_DATA::default(),
            network,
            &mut data,
            structure,
            at,
            at,
            at,
            valence,
            &mut groups,
            None,
            None,
            None,
            input,
            i64::MAX,
            &mut runs,
            &mut delta,
            FORBIDDEN_MASK,
            FIXED_MASK,
            0,
        )
        .unwrap();
        (result, primary, secondary, runs)
    }

    fn primary_structure(
        reversed: SourceMutPointer<INChI>,
        mobile_h: i8,
        atom_count: i32,
    ) -> StrFromINChI {
        let mut structure = StrFromINChI {
            bMobileH: mobile_h,
            num_atoms: atom_count,
            ..StrFromINChI::default()
        };
        structure.pOneINChI[0] = reversed;
        structure
    }

    #[test]
    fn source_port__ichirvr6__fixrestoredstructurestereo__line_59() {
        // No stereo differences: both comparisons are refreshed, all local lists stay empty,
        // and the source exit cleanup performs no source allocation.
        let mut heap = SourceHeap::default();
        let reversed = inchi(&mut heap, None, 0);
        let input = inchi(&mut heap, None, 0);
        let at = atoms(&mut heap, [inp_ATOM::default(), inp_ATOM::default()]);
        let mut structure = primary_structure(reversed, TAUT_YES as i8, 0);
        let mut network = BN_STRUCT::default();
        let baseline = heap.live_allocation_count();
        let (result, primary, secondary, runs) = invoke(
            &mut heap,
            &mut structure,
            &mut network,
            at,
            &mut [],
            [input, SourceMutPointer::null()],
        );
        assert_eq!((result, primary.flags, secondary.flags, runs), (0, 0, 0, 0));
        assert_eq!(heap.live_allocation_count(), baseline);

        // Equal nonzero source error codes set IDIF_PROBLEM and map to RI_ERR_PROGR.
        let mut heap = SourceHeap::default();
        let reversed = inchi(&mut heap, None, 7);
        let input = inchi(&mut heap, None, 7);
        let at = atoms(&mut heap, [inp_ATOM::default(), inp_ATOM::default()]);
        let mut structure = primary_structure(reversed, TAUT_YES as i8, 0);
        let mut network = BN_STRUCT::default();
        let (result, primary, _, runs) = invoke(
            &mut heap,
            &mut structure,
            &mut network,
            at,
            &mut [],
            [input, SourceMutPointer::null()],
        );
        assert_eq!(result, RI_ERR_PROGR);
        assert_eq!(primary.flags & IDIF_PROBLEM, IDIF_PROBLEM);
        assert_eq!(runs, 0);

        // Case 01: a missing input stereobond without a graph bond is a syntax error.
        let mut heap = SourceHeap::default();
        let reversed = inchi(&mut heap, None, 0);
        let input = inchi(&mut heap, Some(1), 0);
        let at = atoms(&mut heap, [inp_ATOM::default(), inp_ATOM::default()]);
        let mut structure = primary_structure(reversed, TAUT_YES as i8, 2);
        let mut network = bond_network(&mut heap, 0, 0, 1);
        let (result, primary, _, _) = invoke(
            &mut heap,
            &mut structure,
            &mut network,
            at,
            &mut [VAL_AT::default(), VAL_AT::default()],
            [input, SourceMutPointer::null()],
        );
        assert_eq!(result, RI_ERR_SYNTAX);
        assert_eq!(primary.flags & IDIF_SB_MISS, IDIF_SB_MISS);
        assert_eq!((primary.num_sb_in2_only, primary.sb_in2_only[0]), (1, 0));

        // Fixed-edge collection uses atom index (not neighbor order), propagates calloc failure,
        // restores the source mask, and leaves no allocation behind.
        let mut heap = SourceHeap::default();
        let reversed = inchi(&mut heap, None, 0);
        let input = inchi(&mut heap, Some(1), 0);
        let at = atoms(
            &mut heap,
            [
                inp_ATOM {
                    valence: 1,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    valence: 1,
                    ..inp_ATOM::default()
                },
            ],
        );
        let mut structure = primary_structure(reversed, TAUT_YES as i8, 2);
        let mut network = bond_network(&mut heap, 1, FIXED_MASK as i8, 1);
        let baseline = heap.live_allocation_count();
        heap.fail_after_allocations(0);
        let (result, _, _, _) = invoke(
            &mut heap,
            &mut structure,
            &mut network,
            at,
            &mut [VAL_AT::default(), VAL_AT::default()],
            [input, SourceMutPointer::null()],
        );
        assert_eq!(result, RI_ERR_ALLOC);
        assert_eq!(heap.live_allocation_count(), baseline);
        assert_eq!(
            heap.slice(network.edge.as_const()).unwrap()[0].forbidden,
            FIXED_MASK as i8
        );

        // Case 02: only the secondary comparison misses a bond; the first radical-list
        // growth is the observable allocation and maps its failure to RI_ERR_ALLOC.
        let mut heap = SourceHeap::default();
        let reversed_primary = inchi(&mut heap, None, 0);
        let input_primary = inchi(&mut heap, None, 0);
        let reversed_secondary = inchi(&mut heap, None, 0);
        let input_secondary = inchi(&mut heap, Some(1), 0);
        let at = atoms(&mut heap, [inp_ATOM::default(), inp_ATOM::default()]);
        let mut structure = primary_structure(reversed_primary, TAUT_NON as i8, 2);
        structure.pOneINChI[1] = reversed_secondary;
        let mut network = bond_network(&mut heap, 0, 0, 1);
        let baseline = heap.live_allocation_count();
        heap.fail_after_allocations(0);
        let (result, primary, secondary, _) = invoke(
            &mut heap,
            &mut structure,
            &mut network,
            at,
            &mut [VAL_AT::default(), VAL_AT::default()],
            [input_primary, input_secondary],
        );
        assert_eq!((result, primary.flags & IDIFF_SB), (RI_ERR_ALLOC, 0));
        assert_eq!(secondary.flags & IDIF_SB_MISS, IDIF_SB_MISS);
        assert_eq!(heap.live_allocation_count(), baseline);

        // Case 03: fixed-H undefined reversed bond plus endpoint data requires a real
        // graph bond; absent adjacency returns the source syntax status before case 04.
        let mut heap = SourceHeap::default();
        let reversed = inchi(&mut heap, Some(AB_PARITY_UNDF as i8), 0);
        let input = inchi(&mut heap, None, 0);
        let at = atoms(&mut heap, [inp_ATOM::default(), inp_ATOM::default()]);
        let mut structure = primary_structure(reversed, TAUT_NON as i8, 2);
        structure.endpoint = heap.allocate_model_storage(vec![0_u16, 0]).unwrap();
        let mut network = BN_STRUCT::default();
        let (result, primary, _, _) = invoke(
            &mut heap,
            &mut structure,
            &mut network,
            at,
            &mut [VAL_AT::default(), VAL_AT::default()],
            [input, SourceMutPointer::null()],
        );
        assert_eq!(result, RI_ERR_SYNTAX);
        assert_eq!(primary.flags & IDIF_SB_EXTRA_UNDF, IDIF_SB_EXTRA_UNDF);
        assert_eq!(primary.num_sb_undef_in1_only, 1);

        // Cases 04 and 05: the undefined mobile-H bond is found, case 04 rejects zero
        // flow, and the NH-specific case 05 independently rejects the same zero flow.
        let mut heap = SourceHeap::default();
        let reversed = inchi(&mut heap, Some(AB_PARITY_UNDF as i8), 0);
        let input = inchi(&mut heap, None, 0);
        let at = atoms(
            &mut heap,
            [
                inp_ATOM {
                    valence: 1,
                    num_H: 1,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    valence: 1,
                    ..inp_ATOM::default()
                },
            ],
        );
        let mut structure = primary_structure(reversed, TAUT_YES as i8, 2);
        let mut network = bond_network(&mut heap, 0, 0, 0);
        let mut valence = [VAL_AT::default(), VAL_AT::default()];
        valence[0].cNumValenceElectrons = 5;
        valence[0].cPeriodicRowNumber = 1;
        let before = network.clone();
        let (result, primary, _, runs) = invoke(
            &mut heap,
            &mut structure,
            &mut network,
            at,
            &mut valence,
            [input, SourceMutPointer::null()],
        );
        assert_eq!(result, 0);
        assert_eq!(primary.flags & IDIF_SB_EXTRA_UNDF, IDIF_SB_EXTRA_UNDF);
        assert_eq!(
            (network.tot_st_cap, network.tot_st_flow, runs),
            (before.tot_st_cap, before.tot_st_flow, 0)
        );
        assert_eq!(heap.slice(network.edge.as_const()).unwrap()[0].flow, 0);

        // Case 06: the secondary fixed-H undefined bond reaches the NH2 predicate and
        // preserves all state when the source bond flow is zero.
        let mut heap = SourceHeap::default();
        let reversed_primary = inchi(&mut heap, None, 0);
        let input_primary = inchi(&mut heap, None, 0);
        let reversed_secondary = inchi(&mut heap, Some(AB_PARITY_UNDF as i8), 0);
        let input_secondary = inchi(&mut heap, None, 0);
        let at = atoms(
            &mut heap,
            [
                inp_ATOM {
                    valence: 1,
                    num_H: 2,
                    ..inp_ATOM::default()
                },
                inp_ATOM {
                    valence: 1,
                    ..inp_ATOM::default()
                },
            ],
        );
        let mut structure = primary_structure(reversed_primary, TAUT_NON as i8, 2);
        structure.pOneINChI[1] = reversed_secondary;
        let mut network = bond_network(&mut heap, 0, 0, 0);
        let mut valence = [VAL_AT::default(), VAL_AT::default()];
        valence[0].cNumValenceElectrons = 5;
        valence[0].cPeriodicRowNumber = 1;
        let (result, primary, secondary, runs) = invoke(
            &mut heap,
            &mut structure,
            &mut network,
            at,
            &mut valence,
            [input_primary, input_secondary],
        );
        assert_eq!((result, primary.flags & IDIFF_SB, runs), (0, 0, 0));
        assert_eq!(secondary.flags & IDIF_SB_EXTRA_UNDF, IDIF_SB_EXTRA_UNDF);
        assert_eq!(secondary.num_sb_undef_in1_only, 1);
        assert_eq!(heap.slice(network.edge.as_const()).unwrap()[0].flow, 0);
    }
}
