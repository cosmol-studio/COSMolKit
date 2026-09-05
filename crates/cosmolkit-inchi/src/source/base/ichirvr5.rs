use crate::source::base::ichiring::is_bond_in_Nmax_memb_ring;
use crate::source::base::ichirvr1::{
    AddToEdgeList, AllocEdgeList, FindInEdgeList, GetChargeFlowerUpperEdge, MakeOneInChIOutOfStrFromINChI2,
    RemoveForbiddenEdgeMask, RunBnsRestoreOnce, RunBnsTestOnce, SetForbiddenEdgeMask,
};
use crate::source::base::ichirvr4::{FillOutCMP2MHINCHI, FillOutExtraFixedHDataRestr};
use crate::source_types::{
    ALL_TC_GROUPS, BFS_Q, BN_DATA, BN_STRUCT, BOND_TYPE_SINGLE, CANON_GLOBALS, CMP2MHINCHI, EDGE_LIST, EDGE_LIST_CLEAR,
    EDGE_LIST_FREE, INCHI_CLOCK, INChI, INPUT_PARMS, MAX_DIFF_FIXH, MAX_DIFF_MOBH, NO_VERTEX, RI_ERR_PROGR,
    STRUCT_DATA, SourceHeap, SourceHeapError, SourceMutPointer, SourceTGroupInfoPointer, StrFromINChI, VAL_AT, clock_t,
    inp_ATOM, tagTCGroupTypes_TCG_Minus as TCG_Minus, tagTCGroupTypes_TCG_Plus as TCG_Plus,
};

const INC_ADD_EDGE: i32 = 64;
const TGRF_MINUS_FIRST: i16 = 1;

#[allow(non_snake_case)]
pub(crate) fn GetPlusMinusVertex(
    heap: &SourceHeap,
    pBNS: &BN_STRUCT,
    pTCGroups: &ALL_TC_GROUPS,
    bCheckForbiddenPlus: i32,
    bCheckForbiddenMinus: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:58 GetPlusMinusVertex
    /*
    int GetPlusMinusVertex( BN_STRUCT *pBNS,
                            ALL_TC_GROUPS *pTCGroups,
                            int bCheckForbiddenPlus,
                            int bCheckForbiddenMinus )
    {
        int k, ePlusSuper, eMinusSuper, vPlusSuper, vPlusMinus1 = NO_VERTEX, vPlusMinus2 = NO_VERTEX; /* djb-rwth: removing redundant variables */
        BNS_EDGE *pEdge;
        if (( k = pTCGroups->nGroup[TCG_Plus] ) >= 0 &&
            ( ePlusSuper = pTCGroups->pTCG[k].nForwardEdge ) > 0 &&
            ( vPlusSuper = pTCGroups->pTCG[k].nVertexNumber ) >= pBNS->num_atoms &&
             !( ( pEdge = pBNS->edge + ePlusSuper )->forbidden && bCheckForbiddenPlus ))
        {

            vPlusMinus1 = pEdge->neighbor12 ^ vPlusSuper;
        }
        if (( k = pTCGroups->nGroup[TCG_Minus] ) >= 0 &&
            ( eMinusSuper = pTCGroups->pTCG[k].nForwardEdge ) > 0 &&
            ( pTCGroups->pTCG[k].nVertexNumber ) >= pBNS->num_atoms && /* djb-rwth: removing redundant code */
             !( ( pEdge = pBNS->edge + eMinusSuper )->forbidden && bCheckForbiddenMinus ))
        {

            vPlusMinus2 = pEdge->neighbor12 ^ eMinusSuper;
        }
        if ((bCheckForbiddenPlus && NO_VERTEX == vPlusMinus1) ||
             (bCheckForbiddenMinus && NO_VERTEX == vPlusMinus2)) /* djb-rwth: addressing LLVM warnings */
        {
            return NO_VERTEX;
        }

        return ( NO_VERTEX != vPlusMinus1 ) ? vPlusMinus1 : vPlusMinus2;
    }
    */
    // END INCHI C FUNCTION: GetPlusMinusVertex
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetPlusMinusVertex
    // INCHI✔️❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: The minus branch intentionally XORs neighbor12 with eMinusSuper, exactly as the active source does.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetPlusMinusVertex

    let mut plus_minus_1 = NO_VERTEX;
    let mut plus_minus_2 = NO_VERTEX;

    let group_index = pTCGroups.nGroup[TCG_Plus as usize];
    if group_index >= 0 {
        let group = heap
            .slice(pTCGroups.pTCG.as_const())?
            .get(group_index as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let edge_index = group.nForwardEdge;
        if edge_index > 0 && group.nVertexNumber >= pBNS.num_atoms {
            let edge = heap
                .slice(pBNS.edge.as_const())?
                .get(edge_index as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if !(edge.forbidden != 0 && bCheckForbiddenPlus != 0) {
                plus_minus_1 = i32::from(edge.neighbor12) ^ group.nVertexNumber;
            }
        }
    }

    let group_index = pTCGroups.nGroup[TCG_Minus as usize];
    if group_index >= 0 {
        let group = heap
            .slice(pTCGroups.pTCG.as_const())?
            .get(group_index as usize)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let edge_index = group.nForwardEdge;
        if edge_index > 0 && group.nVertexNumber >= pBNS.num_atoms {
            let edge = heap
                .slice(pBNS.edge.as_const())?
                .get(edge_index as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if !(edge.forbidden != 0 && bCheckForbiddenMinus != 0) {
                plus_minus_2 = i32::from(edge.neighbor12) ^ edge_index;
            }
        }
    }

    if (bCheckForbiddenPlus != 0 && plus_minus_1 == NO_VERTEX)
        || (bCheckForbiddenMinus != 0 && plus_minus_2 == NO_VERTEX)
    {
        return Ok(NO_VERTEX);
    }
    Ok(if plus_minus_1 != NO_VERTEX {
        plus_minus_1
    } else {
        plus_minus_2
    })
}

#[allow(non_snake_case)]
pub(crate) fn bIsUnsatCarbonInASmallRing(
    heap: &mut SourceHeap,
    at2: SourceMutPointer<inp_ATOM>,
    pVA: &[VAL_AT],
    iat: i32,
    pbfsq: &BFS_Q,
    min_ring_size: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:92 bIsUnsatCarbonInASmallRing
    /*
    int bIsUnsatCarbonInASmallRing( inp_ATOM *at2,
                                    VAL_AT *pVA,
                                    int iat,
                                    BFS_Q *pbfsq,
                                    int min_ring_size )
    {
        int j, nCurRingSize, nMinRingSize;
        if (min_ring_size < 5)
        {
            /* =C= in a small ring  */
            if (at2[iat].valence == 2 &&
                 pVA[iat].cMinRingSize <= 5 &&
                 at2[iat].chem_bonds_valence == 4)
            {
                return 1;
            }
        }
        else
        {
            if (at2[iat].valence == 2 &&
                 pVA[iat].cMinRingSize &&
                 pVA[iat].cMinRingSize <= min_ring_size &&
                 at2[iat].chem_bonds_valence == 3)
            {
                return 1;
            }
            nCurRingSize = nMinRingSize = min_ring_size + 1;
            if (( at2[iat].valence == 2 || at2[iat].valence == 3 ) &&
                 at2[iat].chem_bonds_valence == at2[iat].valence + 1)
            {
                for (j = 0; j < at2[iat].valence; j++)
                {
                    nCurRingSize = is_bond_in_Nmax_memb_ring( at2, iat, j, pbfsq->q,
                                                 pbfsq->nAtomLevel,
                                                 pbfsq->cSource, (AT_RANK) nMinRingSize /* max ring size */ );
                    if (0 < nCurRingSize && nCurRingSize < nMinRingSize)
                    {
                        nMinRingSize = nCurRingSize;
                    }
                }
                return ( 0 <= nCurRingSize ) ? ( nMinRingSize <= min_ring_size ) : nCurRingSize;
            }
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: bIsUnsatCarbonInASmallRing
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: bIsUnsatCarbonInASmallRing
    // INCHI✔️❌: COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux; READ_INCHI_STRING=1.
    // INCHI✔️❌: AT_RANK is unsigned 16-bit; is_bond_in_Nmax_memb_ring is the existing complete source port.
    // END INCHI ACTIVE MACRO CONFIGURATION: bIsUnsatCarbonInASmallRing

    let index = usize::try_from(iat).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom = heap
        .slice(at2.as_const())?
        .get(index)
        .cloned()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let valence = pVA.get(index).ok_or(SourceHeapError::PointerOutOfBounds)?;
    if min_ring_size < 5 {
        if atom.valence == 2 && valence.cMinRingSize <= 5 && atom.chem_bonds_valence == 4 {
            return Ok(1);
        }
    } else {
        if atom.valence == 2
            && valence.cMinRingSize != 0
            && i32::from(valence.cMinRingSize) <= min_ring_size
            && atom.chem_bonds_valence == 3
        {
            return Ok(1);
        }
        let mut current_ring_size = min_ring_size.wrapping_add(1);
        let mut minimum_ring_size = current_ring_size;
        if (atom.valence == 2 || atom.valence == 3) && i32::from(atom.chem_bonds_valence) == i32::from(atom.valence) + 1
        {
            let mut neighbor_order = 0_i32;
            while neighbor_order < i32::from(atom.valence) {
                current_ring_size = is_bond_in_Nmax_memb_ring(
                    heap,
                    at2,
                    iat,
                    neighbor_order,
                    pbfsq.q,
                    pbfsq.nAtomLevel,
                    pbfsq.cSource,
                    minimum_ring_size as u16,
                )?;
                if current_ring_size > 0 && current_ring_size < minimum_ring_size {
                    minimum_ring_size = current_ring_size;
                }
                neighbor_order = neighbor_order.wrapping_add(1);
            }
            return Ok(if current_ring_size >= 0 {
                i32::from(minimum_ring_size <= min_ring_size)
            } else {
                current_ring_size
            });
        }
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FixMobileHRestoredStructure(
    heap: &mut SourceHeap,
    pCG: &mut CANON_GLOBALS,
    ic: SourceMutPointer<INCHI_CLOCK>,
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
    _bHasSomeFixedH: i32,
    _pnNumRunBNS: Option<&mut i32>,
    _pnTotalDelta: Option<&mut i32>,
    forbidden_edge_mask: i32,
    forbidden_stereo_edge_mask: i32,
    clock_result: clock_t,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:141 FixMobileHRestoredStructure
    // INCHI✔️❌: complete active source frame follows verbatim; implementation
    // INCHI✔️❌: preserves active behavior, with known SourceHeap overhead.
    /*
    int FixMobileHRestoredStructure( CANON_GLOBALS *pCG,
                                     INCHI_CLOCK *ic,
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
                                     int bHasSomeFixedH,
                                     int *pnNumRunBNS,
                                     int *pnTotalDelta,
                                     int forbidden_edge_mask,
                                     int forbidden_stereo_edge_mask )
    {
        /*--------- process extra or missing Fixed-H on non-tautomeric atoms ------*/
        /* at2 should be the most recently restored atom, Fixed-H */
        int i, j, k, iat, delta, cur_success, ret = 0; /* djb-rwth: removing redundant variables/code */
        CMP2MHINCHI c2i;
        CMP2MHINCHI *pc2i = &c2i;

        EDGE_LIST AllChargeEdges, CurrEdges, CurrEdges2, CurrEdges3, TautEdges, NFlowerEdges, OtherNFlowerEdges, FixedLargeRingStereoEdges;
        EDGE_LIST  *pEdgeList = NULL;

        EdgeIndex e;
        BNS_EDGE  *pe;
        Vertex v1, v2, vPlusMinus;
        BNS_VERTEX *pv1, *pv2;

        Vertex     vPathStart, vPathEnd;
        int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;

        int        forbidden_edge_mask_inv = ~forbidden_edge_mask; /* djb-rwth: removing redundant variables */

        INCHI_HEAPCHK

        AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &CurrEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &NFlowerEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &CurrEdges2, EDGE_LIST_CLEAR );
        AllocEdgeList( &CurrEdges3, EDGE_LIST_CLEAR );
        AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &FixedLargeRingStereoEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &TautEdges, EDGE_LIST_CLEAR );

        /* djb-rwth: removing redundant code */

        if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
        {
            goto exit_function;  /* no fixed-H found */
        }
        /* taut group edges */
        for (i = 0; i < pTCGroups->num_tgroups; i++)
        {
            pv1 = pBNS->vert + ( v1 = pTCGroups->pTCG[i].nVertexNumber ); /* t-group vertex */ /* djb-rwth: ignoring LLVM warning: see comments below */
            for (j = 0; j < pv1->num_adj_edges; j++)
            {
                /* e, pe - tautomeric atom edge; pv2 - endpoint vertex */
                /* Note: pe, pv2, v1 are not used here; they are to show how to traverse t-group */
                pv2 = pBNS->vert + ( pe = pBNS->edge + ( e = pv1->iedge[j] ) )->neighbor1; /* djb-rwth: ignoring LLVM warning: see comments above */
                if ((ret = AddToEdgeList( &TautEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
            }
        }
        /* charge and flower edges */
        for (i = 0; i < pStruct->num_atoms; i++)
        {
            if (( e = pVA[i].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                ( ret = AddToEdgeList( &AllChargeEdges, e, INC_ADD_EDGE ) ))
            {
                goto exit_function;
            }
            if (( e = pVA[i].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden)
            {
                if ((ret = AddToEdgeList( &AllChargeEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }

                /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
                if (pVA[i].cNumValenceElectrons == 5 && !pVA[i].cMetal && /* N, P, As */
                     NO_VERTEX != ( j = GetChargeFlowerUpperEdge( pBNS, pVA, e ) ))
                {

                    if (!pBNS->edge[j].forbidden && pBNS->edge[j].flow)
                    {
                        if ((ret = AddToEdgeList( &AllChargeEdges, j, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                        if ((ret = AddToEdgeList( &NFlowerEdges, j, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                    else
                    {
                        if ((ret = AddToEdgeList( &OtherNFlowerEdges, j, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                }
            }
        }
        if (forbidden_stereo_edge_mask)
        {
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                for (j = 0; j < at2[i].valence; j++)
                {
                    if (pBNS->edge[k = pBNS->vert[i].iedge[j]].forbidden == forbidden_stereo_edge_mask)
                    {
                        int nMinRingSize = is_bond_in_Nmax_memb_ring( at2, i, j, pStruct->pbfsq->q,
                                                                 pStruct->pbfsq->nAtomLevel,
                                                                 pStruct->pbfsq->cSource, 99 /* max ring size */ );
                        if (0 < nMinRingSize && ( ret = AddToEdgeList( &FixedLargeRingStereoEdges, k, INC_ADD_EDGE ) ))
                        {
                            goto exit_function;
                        }
                    }
                }
            }
        }

        INCHI_HEAPCHK

        if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
        {
            goto exit_function;
        }

        INCHI_HEAPCHK
        if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
        {
            goto exit_function;
        }

        INCHI_HEAPCHK




        if (pc2i->nNumTgInChI == 1 && ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1 ) &&
                pc2i->nNumTgDBNMinusRevrs + pc2i->nNumTgNHMinusRevrs == 0 && pc2i->nNumTgOMinusInChI &&
                !( pTCGroups->pTCG[0].tg_RestoreFlags & TGRF_MINUS_FIRST ))
        {            /*----------------------------------------------------*/
            /* case 01: restored has -O(-) and does not have N(-) */
            /*          endpoints defined by the original InChI   */
            /*          restored has single taut. group or more   */
            /*          tautomeric endpoints.                     */
            /* Solution: move (-) from endp. -O(-) to endpoints N */
            /*-------
    ---------------------------------------------*/
            pTCGroups->pTCG[0].tg_RestoreFlags |= TGRF_MINUS_FIRST;
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
            if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
            {
                goto exit_function;  /* no fixed-H found */
            }
            if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
            if (!pc2i->bHasDifference)
            {
                goto exit_function; /* nothing to do */
            }
        }
        if (pc2i->nNumTgInChI == 1 && ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1 ) &&
             pc2i->nNumTgDBNMinusRevrs + pc2i->nNumTgNHMinusRevrs == 0 && pc2i->nNumTgOMinusInChI == 0)
        {
            /*-------------------------------------------------------*/
            /* case 02: restored has no -O(-) and does not have N(-) */
            /*          restored has single taut. group or more      */
            /*          tautomeric endpoints.                        */
            /* Solution: >N-AB=N-  => >N(+)=AB-NH- (add H(+))        */
            /* Solution: >N-AB=NH  => >N(+)=AB-NH2 (add H(+))        */
            /*      SB_N_III  DB_N_III                               */
            /*-------------------------------------------------------*/
            int iat_SB_N_III[MAX_DIFF_MOBH], iat_DB_N_III[MAX_DIFF_MOBH];
            int num_SB_N_III = 0, num_DB_N_III = 0, k1, k2;
            CurrEdges.num_edges = 0;
            cur_success = 0;
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                iat = i;
                if (pVA[iat].cNumValenceElectrons == 5 && pVA[i].cPeriodicRowNumber == 1 &&
                     !at2[iat].endpoint && !at2[iat].charge && !at2[iat].radical)
                {
                    if (num_DB_N_III < MAX_DIFF_MOBH && !at2[iat].num_H &&
                         at2[iat].valence == 2 &&
                         at2[iat].chem_bonds_valence == 3 &&
                         !at2[iat].sb_parity[0] &&  /* do not eliminate stereobonds */
                         ( e = pVA[iat].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                         pBNS->edge[e].cap && !pBNS->edge[e].flow)
                    {
                        /* -N= */
                        iat_DB_N_III[num_DB_N_III++] = iat;
                    }
                    else
                    {
                        if (num_DB_N_III < MAX_DIFF_MOBH && 1 == at2[iat].num_H &&
                             at2[iat].valence == 1 &&
                             at2[iat].chem_bonds_valence == 2 &&
                             !at2[iat].sb_parity[0] &&  /* do not eliminate stereobonds */
                             ( e = pVA[iat].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                             pBNS->edge[e].cap && !pBNS->edge[e].flow)
                        {
                            /* -N= */
                            iat_DB_N_III[num_DB_N_III++] = iat;
                        }
                        else
                        {
                            if (num_SB_N_III < MAX_DIFF_MOBH && !at2[iat].num_H &&
                                    at2[iat].valence == 3 &&
                                    at2[iat].chem_bonds_valence == 3 &&
                                    ( e = pVA[iat].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                                    pBNS->edge[e].cap && pBNS->edge[e].flow)
                            {
                                /* -N< */
                                iat_SB_N_III[num_SB_N_III++] = iat;
                                if ((ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_function;
                                }
                            }
                        }
                    }
                }
            }
            if (num_DB_N_III && num_SB_N_III)
            {
                EdgeIndex ieMinus;
                BNS_EDGE  *peMinus;
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
                for (i = 0; i < num_DB_N_III && !cur_success; i++)
                {
                    iat = iat_DB_N_III[i];
                    e = pBNS->edge[k1 = pBNS->vert[iat].iedge[0]].flow ? k1 :
                        pBNS->edge[k2 = pBNS->vert[iat].iedge[1]].flow ? k2 : NO_VERTEX;
                    if (e == NO_VERTEX)
                    {
                        continue; /* should not happen */
                    }
                    ieMinus = pVA[iat].nCMinusGroupEdge - 1;
                    peMinus = pBNS->edge + ieMinus;
                    pe = pBNS->edge + e;
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                    pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                    pe->forbidden |= forbidden_edge_mask;     /* fix double bond */
                    peMinus->forbidden &= forbidden_edge_mask_inv; /* allow negative charge */
                    delta = 1;
                    pe->flow -= delta; /* remove (-) from AB-O(-) */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 2) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Added (-)charge -N= and (+) to -N< => nDeltaCharge == 2 */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if (ret > 0)
                        {
                            /* djb-rwth: removing redundant code */
                            cur_success++; /* 01 */

                            /* eliminate (-) charge and add H */
                            pv1 = pBNS->vert + ( v1 = peMinus->neighbor1 );      /* atom */
                            pv2 = pBNS->vert + ( v2 = peMinus->neighbor12 ^ v1 );/* (=) vertex */ /* djb-rwth: ignoring LLVM warning: consistency of the code */
                            /* effectively eliminate (-) edge by setting its cap=flow= 0 */
                            peMinus->cap--;
                            peMinus->flow--;
                            pv1->st_edge.cap--;
                            pv1->st_edge.flow--;
                            pv2->st_edge.cap--;
                            pv2->st_edge.flow--;
                            pBNS->tot_st_flow -= 2;
                            pBNS->tot_st_cap -= 2;
                            /* add H */
                            pStruct->at[iat].num_H++;
                            /* register total charge increase */
                            pTCGroups->total_charge++;
                            pStruct->nNumRemovedProtonsByRevrs -= 1;
                        }
                    }
                    else
                    {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        peMinus->forbidden |= forbidden_edge_mask;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                CurrEdges.num_edges = 0; /* clear current edge list */

                if (cur_success)
                {
                    /* djb-rwth: removing redundant code */
                    /* recalculate InChI from the structure */
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
                    if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                    if (!pc2i->bHasDifference)
                    {
                        goto exit_function; /* nothing to do */
                    }
                }
            }
        }
        if (pc2i->nNumTgInChI == 1 && ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1 ) && /* ADP */
            pc2i->nNumTgMInChI == 0 && pc2i->nNumTgNInChI && pc2i->nNumTgOInChI)
        {
            /*-------------------------------------------------------*/
            /* case 03: restored has N and O endpoints, no (-) endp  */
            /* case 04: original has single taut. group or more      */
            /*          tautomeric endpoints.                        */
            /* Solution: 1. Move taut attachment from O to N         */
            /* Solution: 2. Replace the attachment with (-)          */
            /*      SB_N_III  DB_N_III                               */
            /*-------------------------------------------------------*/
            /*
              int iat_SB_N_III[MAX_DIFF_MOBH], iat_DB_N_III[MAX_DIFF_MOBH];
              int num_SB_N_III = 0, num_DB_N_III = 0, k1, k2,
            */
            int itg, j1, j2, bAction = 0;
            BNS_VERTEX *pTg, *pvEndp, *pvEndp2, *pvCent; /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
            Vertex     vEndp, vEndp2, vCent;
            BNS_EDGE   *peTg, *peTg2, *peCent1, *peCent2;
            EdgeIndex  eTg, eTg2;

            CurrEdges.num_edges = 0;
            CurrEdges2.num_edges = 0;
            cur_success = 0;

            /* 1st attempt: -NH-=O => -N(-)-=O  or -N=-OH => -N(-)-=O */
            for (itg = 0; itg < pTCGroups->num_tgroups && !cur_success; itg++)
            {
                pTg = pBNS->vert + pTCGroups->pTCG[itg].nVertexNumber;
                for (i = 0; i < pTg->num_adj_edges && !cur_success; i++)
                {
                    pvEndp = pBNS->vert + ( vEndp = ( peTg = pBNS->edge + ( eTg = pTg->iedge[i] ) )->neighbor1 ); /* djb-rwth: ignoring LLVM warning: value used */
                    eTg2 = -1;
                    if (pVA[vEndp].cNumValenceElectrons == 6 && peTg->cap)
                    {
                        /* endpoint -OH or =O found; search for a possible centerpoint */
                        for (j1 = 0; j1 < at2[vEndp].valence && eTg2 < 0; j1++)
                        {
                            peCent1 = pBNS->edge + pvEndp->iedge[j1]; /* edge from O to a centerpoint */
                            pvCent = pBNS->vert + ( vCent = peCent1->neighbor12 ^ vEndp ); /* centerpoint */
                            if (at2[vCent].endpoint || !peCent1->cap ||
                                 peCent1->flow + ( peTg->cap == peTg->flow ) != 1)
                            {
                                continue;
                            }
                            /* search for another endpoint, N, around vCent */
                            for (j2 = 0; j2 < at2[vCent].valence; j2++)
                            {
                                peCent2 = pBNS->edge + pvCent->iedge[j2];
                                pvEndp2 = pBNS->vert + ( vEndp2 = peCent2->neighbor12 ^ vCent ); /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
                                if (!peCent2->cap || peCent2->flow + peCent1->flow != 1 ||
                                     at2[vEndp2].endpoint != itg + 1 ||
                                     pVA[vEndp2].cNumValenceElectrons != 5 ||
                                     0 > ( j = pVA[vEndp2].nTautGroupEdge - 1 ) ||
                                     ( peTg2 = pBNS->edge + j )->forbidden ||
                                     peCent2->flow + ( peTg2->cap == peTg2->flow ) != 1)
                                {
                                    continue;
                                }
                                eTg2 = j;
                                break; /* found OH-C=N- or O=C-NH- */
                            }
                        }
                    }
                    if (eTg2 >= 0)
                    {
                        /*--------------------------------------------
                                        tg                        tg
                                eTg //\ eTg2              eTg / \\eTg2
                                    //  \                     /   \\
                            vEndp HO--C==N vEndp2 -->  vEndp O==C--NH vEndp2
                                    ^ ^ ^                     ^ ^ ^
                                eCent1 | eCent2           eCent1 | eCent2
                                        vCent                     vCent

                        additional action: -OH-C=N- => O=C-NH-
                        -------------------------------------------*/
                        if (0 == peTg->cap - peTg->flow && 1 == peTg2->cap - peTg2->flow &&
                             0 == peCent1->flow && 1 == peCent2->flow)
                        {
                            peTg->flow--;          /* 03 prepare */
                            peTg2->flow++;
                            peCent2->flow--;
                            peCent1->flow++;
                            bAction |= 1; /* switched H position */
                        }
                        if (1 == peTg->cap - peTg->flow && 0 == peTg2->cap - peTg2->flow &&
                             1 == peCent1->flow && 0 == peCent2->flow)
                        {
                            /* replace -NH- with -N(-)- */
                            pTCGroups->pTCG[itg].tg_num_H--;
                            pTCGroups->pTCG[itg].tg_num_Minus++;
                            pTCGroups->pTCG[itg].tg_RestoreFlags |= TGRF_MINUS_FIRST;
                            pTCGroups->pTCG[itg].tg_set_Minus = vEndp2 + 1;
                            pStruct->ti.t_group[itg].num[1] ++; /* increment number of (-), keep number of taut attachments */
                            pTCGroups->total_charge--;
                            pTCGroups->tgroup_charge--;
                            pStruct->nNumRemovedProtonsByRevrs += 1;
                            bAction |= 2; /* single NH (at2[vEndp2]) replaced with N(-) */
                            cur_success++; /* 03/04 */
                        }
                    }
                }
            }

            if (0 == pc2i->nNumTgNHInChI + pc2i->nNumTgNH2InChI && pc2i->nNumTgOHInChI && !cur_success)
            {
                /* transfer an attachement to N */
                for (itg = 0; itg < pTCGroups->num_tgroups; itg++)
                {
                    pTg = pBNS->vert + pTCGroups->pTCG[itg].nVertexNumber;
                    for (i = 0; i < pTg->num_adj_edges; i++)
                    {
                        pvEndp = pBNS->vert + ( vEndp = ( peTg = pBNS->edge + ( eTg = pTg->iedge[i] ) )->neighbor1 );
                        if (pVA[vEndp].cNumValenceElectrons == 6 &&
                             at2[vEndp].valence == at2[vEndp].chem_bonds_valence &&
                             peTg->flow && peTg->flow == peTg->cap)
                        {
                            /* endpoint -OH found; save the tautomeric group edge */
                            if ((ret = AddToEdgeList( &CurrEdges, eTg, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                        else
                        {
                            if (pVA[vEndp].cNumValenceElectrons == 5 &&
                                 pVA[vEndp].cPeriodicRowNumber == 1 &&
                                 at2[vEndp].valence + 1 == at2[vEndp].chem_bonds_valence &&
                                 peTg->cap && peTg->flow + 1 == peTg->cap)
                            {
                                /* endpoint -N= or =NH found, check for -N=-OH */
                                e = -1;
                                for (j1 = 0; j1 < at2[vEndp].valence && e < 0; j1++)
                                {
                                    peCent1 = pBNS->edge + pvEndp->iedge[j1];
                                    if (peCent1->flow == 1)
                                    {
                                        /* double bond */
                                        pvCent = pBNS->vert + ( vCent = peCent1->neighbor12 ^ vEndp );
                                        if (at2[vCent].endpoint)
                                            continue;
                                        for (j2 = 0; j2 < at2[vCent].valence; j2++)
                                        {
                                            peCent2 = pBNS->edge + pvCent->iedge[j2];
                                            pvEndp2 = pBNS->vert + ( vEndp2 = peCent2->neighbor12 ^ vCent ); /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
                                            if (peCent2->flow || at2[vEndp2].endpoint != itg + 1 ||
                                                 pVA[vEndp2].cNumValenceElectrons != 6 ||
                                                 0 >= ( e = pVA[vEndp2].nTautGroupEdge - 1 ) ||
                                                 pBNS->edge[e].forbidden || !pBNS->edge[e].flow)
                                            {
                                                e = -1;
                                                continue;
                                            }
                                            /*********************/
                                            /* found -N=X-OH     */
                                            /*    vEndp ^ vEndp2 */
                                            /*          vCent    */
                                            /*********************/
                                            /* save this -OH taut edge */
                                            if ((ret = AddToEdgeList( &CurrEdges2, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                            {
                                                goto exit_function;
                                            }
                                            break;
                                        }
                                    }
                                }
                                if (e < 0 && ( ret = AddToEdgeList( &CurrEdges, eTg, INC_ADD_EDGE ) ))
                                {
                                    goto exit_function;
                                }
                            }
                        }
                    }
                }
                /* rearrange the flows */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                SetForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
                SetForbiddenEdgeMask( pBNS, &CurrEdges2, forbidden_edge_mask );
                pEdgeList = CurrEdges2.num_edges ? &CurrEdges2 : CurrEdges.num_edges ? &CurrEdges : NULL;

                for (i = 0; pEdgeList && i < pEdgeList->num_edges && !cur_success; i++)
                {
                    pe = pBNS->edge + pEdgeList->pnEdges[i]; /* pe->flow = 1 <=> -OH */
                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + ( v1 = pe->neighbor1 );       /* -OH atom */
                    pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 ); /* t-group vertex */
                    /* locate the t-group */
                    for (itg = 0; itg < pTCGroups->num_tgroups; itg++)
                    {
                        if (v2 == pTCGroups->pTCG[itg].nVertexNumber)
                        {
                            break;
                        }
                    }
                    if (itg == pTCGroups->num_tgroups)
                    {
                        /* tgroup not found -- should not happen */
                        continue;
                    }

                    delta = 1;
                    pe->flow -= delta; /* add one attachment to  */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 2) /* djb-rwth: addressing LLVM warning */
                    {
                        /* Added (-)charge -N= and (+) to -N< => nDeltaCharge == 2 */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if (ret > 0)
                        {
                            /* djb-rwth: removing redundant code */
                            cur_success++; /* 03 */
                            /* replace -NH- with -N(-)- */
                            pTCGroups->pTCG[itg].tg_num_H--;
                            pTCGroups->pTCG[itg].tg_num_Minus++;
                            pTCGroups->pTCG[itg].tg_RestoreFlags |= TGRF_MINUS_FIRST;
                            pStruct->ti.t_group[itg].num[1] ++;
                            pTCGroups->total_charge--;
                            pTCGroups->tgroup_charge--;
                            pStruct->nNumRemovedProtonsByRevrs += 1;
                            bAction |= 4; /* H in the 1st available NH was replaced with (-) */
                        }
                    }
                    else
                    {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
            }
            else
            {
                if (pc2i->nNumTgNHInChI + pc2i->nNumTgNH2InChI && pc2i->nNumTgOInChI && !cur_success)
                {
                    /* change an attachement to N from H to (-) */
                    for (itg = 0; itg < pTCGroups->num_tgroups && !cur_success; itg++)
                    {
                        pTg = pBNS->vert + pTCGroups->pTCG[itg].nVertexNumber;
                        for (i = 0; i < pTg->num_adj_edges && !cur_success; i++)
                        {
                            pvEndp2 = pBNS->vert + ( vEndp2 = ( peTg = pBNS->edge + pTg->iedge[i] )->neighbor1 ); /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
                            if (pVA[vEndp2].cNumValenceElectrons == 5 && pVA[vEndp2].cPeriodicRowNumber == 1 &&
                                 at2[vEndp2].valence == at2[vEndp2].chem_bonds_valence &&
                                 peTg->flow && peTg->flow == peTg->cap)
                            {
                                /* endpoint -NHn found; change its charge */
                                cur_success++; /* 04 */
                                /* replace -NH- with -N(-)- */
                                pTCGroups->pTCG[itg].tg_num_H--;
                                pTCGroups->pTCG[itg].tg_num_Minus++;
                                pTCGroups->pTCG[itg].tg_RestoreFlags |= TGRF_MINUS_FIRST;
                                pTCGroups->pTCG[itg].tg_set_Minus = vEndp2 + 1;
                                pStruct->ti.t_group[itg].num[1] ++;
                                pTCGroups->total_charge--;
                                pTCGroups->tgroup_charge--;
                                pStruct->nNumRemovedProtonsByRevrs += 1;
                                bAction |= 8; /* manually set (-) charge to NH atom, vEndp2 */
                            }
                        }
                    }
                }
            }
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
                if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (pStruct->One_ti.num_t_groups == 1 && pStruct->One_ti.t_group[0].num[1])
                {
                    /* this method did not work: no alt path from N(-) to =O */
                    itg = 0;
                    if (bAction & ( 8 | 2 ))
                    {
                        /* roll back NH -> N(-) replacement; H move from OH to N is not undone */
                        pTCGroups->pTCG[itg].tg_num_H++;
                        pTCGroups->pTCG[itg].tg_num_Minus--;
                        pTCGroups->pTCG[itg].tg_RestoreFlags &= ~TGRF_MINUS_FIRST;
                        pTCGroups->pTCG[itg].tg_set_Minus = 0;
                        pStruct->ti.t_group[itg].num[1] --;
                        pTCGroups->total_charge++;
                        pTCGroups->tgroup_charge++;
                        pStruct->nNumRemovedProtonsByRevrs -= 1;
                        cur_success--;
                    }
                    else
                    {
                        if (bAction & 4)
                        {
                            pTCGroups->pTCG[itg].tg_num_H++;
                            pTCGroups->pTCG[itg].tg_num_Minus--;
                            pTCGroups->pTCG[itg].tg_RestoreFlags &= ~TGRF_MINUS_FIRST;
                            pStruct->ti.t_group[itg].num[1] --;
                            pTCGroups->total_charge++;
                            pTCGroups->tgroup_charge++;
                            pStruct->nNumRemovedProtonsByRevrs -= 1;
                            cur_success--;
                        }
                        else
                        {
                            ret = RI_ERR_PROGR;
                            goto exit_function;
                        }
                    }
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
                    if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }

        if (pc2i->nNumTgInChI == 1 && ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1 ) && /* ADP */
            pc2i->nNumTgMInChI == 0 && ( pc2i->nNumTgNInChI || pc2i->nNumTgOInChI ) &&
            NO_VERTEX != ( vPlusMinus = GetPlusMinusVertex( pBNS, pTCGroups, 1, 1 ) ))
        {
            /*---------------------------------------------------------------------------*/
            /* case 05: restored has N endpoints, no (-) endpoints                       */
            /*          original has single taut. group or more                          */
            /*          tautomeric endpoints.                                            */
            /* Solution: Find -N< and allow (+) charge change                            */
            /*           Fix all charges and taut attachments exept                      */
            /*           =N- and =O (taut. endpoints)                                    */
            /*           Increment st_edge.cap on (+/-) vertex => add (+) charge to -N<  */
            /*           Increment tot. charge in other places                           */
            /*           Increment t-group st_edge.cap                                   */
            /*           Run BNS                                                         */
            /*                                                                           */
            /*      (+/-)*               (+/-)           Result:                         */
            /*        |                    ||                                            */
            /*        |                    ||            - Added (+) to -N<              */
            /*       (+)super             (+)super       - Added attachment point to O   */
            /*        ||                   |                                             */
            /*        ||          =>       |             To make this attachment H,      */
            /*       (Y)                  (Y)            increment                       */
            /*        |                    ||            pTCGroups->pTCG[itg].tg_num_H   */
            /*        |                    ||                                            */
            /*       (+)hetero            (+)hetero      Technical details:              */
            /*         \\                   \            increase capacities of          */
            /*           N                    N(+)       edges to (+/-) otherwise        */
            /*           |                    ||         flow may not be able to         */
            /*   *(t)--O=R.            (t)==O-R.         increase                        */
            /*                                                                           */
            /*                                                                           */
            /*---------------------------------------------------------------------------*/
            int itg;
            BNS_VERTEX *pTg, *pvEndp; /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
            Vertex     vEndp, vTg;
            BNS_EDGE   *peTg;
            EdgeIndex  eTg;
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];

            CurrEdges.num_edges = 0;
            CurrEdges2.num_edges = 0;
            cur_success = 0;
            /* find -N< and non-taut =N- or =O */
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                iat = nCanon2AtnoRevrs[i];
                /* -N< */
                if (!at2[iat].endpoint && !at2[iat].charge && !at2[iat].radical && !at2[iat].num_H &&
                     pVA[i].cNumValenceElectrons == 5 && pVA[i].cPeriodicRowNumber == 1 &&
                     0 <= ( e = pVA[iat].nCPlusGroupEdge - 1 ) && pBNS->edge[e].flow && !pBNS->edge[e].forbidden)
                {

                    if ((ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
            if (!CurrEdges.num_edges)
            {
                goto exit_case_05;
            }
            /* find taut -N= and =O */
            for (itg = 0; itg < pTCGroups->num_tgroups && !cur_success; itg++)
            {
                CurrEdges2.num_edges = 0;
                pTg = pBNS->vert + ( vTg = pTCGroups->pTCG[itg].nVertexNumber );
                for (i = 0; i < pTg->num_adj_edges; i++)
                {
                    pvEndp = pBNS->vert + ( vEndp = ( peTg = pBNS->edge + ( eTg = pTg->iedge[i] ) )->neighbor1 ); /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
                    if (at2[vEndp].charge || at2[vEndp].radical || peTg->cap - peTg->flow != 1)
                    {
                        continue;
                    }
                    /* t-group edges to -N= and =O */
                    if ((ret = AddToEdgeList( &CurrEdges2, eTg, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
                if (!CurrEdges2.num_edges)
                {
                    goto exit_case_05;
                }
                /* fix all charge edges except -N< and all taut. edges except =O and =N- */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                SetForbiddenEdgeMask( pBNS, &TautEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurrEdges2, forbidden_edge_mask );
                delta = 1;
                /* Increment st_edge.cap on (+/-) vertex */
                pBNS->vert[vPlusMinus].st_edge.cap += delta;
                /* Increment st_edge.cap on t-group */
                pTg->st_edge.cap += delta;
                /* total cap count */
                pBNS->tot_st_cap += 2 * delta;

                v1 = vPlusMinus;
                v2 = vTg;

                /* increase capacities of edges to Y  */
                for (i = 0; i < pBNS->vert[vPlusMinus].num_adj_edges; i++)
                {
                    j = pBNS->edge[pBNS->vert[vPlusMinus].iedge[i]].neighbor12 ^ vPlusMinus;
                    for (k = 0; k < pBNS->vert[j].num_adj_edges; k++)
                    {
                        pBNS->edge[pBNS->vert[j].iedge[k]].cap += delta;
                    }
                }

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                    (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warnings */
                {
                    /* Added (+)charge to -N< => nDeltaCharge == 1 */
                    /* Flow change on pe (-)charge edge (atom B-O(-)) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if (ret > 0)
                    {
                        /* djb-rwth: removing redundant code */
                        cur_success++; /* 01 */
                        /* update bookkeeping */
                        pTCGroups->total_charge += delta;
                        pTCGroups->pTCG[itg].edges_cap += delta;
                        pTCGroups->pTCG[itg].tg_num_H += delta;
                        pStruct->nNumRemovedProtonsByRevrs -= delta;
                    }
                }
                else
                {
                    pBNS->vert[vPlusMinus].st_edge.cap -= delta;
                    pTg->st_edge.cap -= delta;
                    /*pTCGroups->pTCG[itg].edges_cap     -= delta;*/ /* ???bug??? - commented out 2006-03-22 */
                    pBNS->tot_st_cap -= 2 * delta;
                    /* decrease capacities of edges to Y  */
                    for (i = 0; i < pBNS->vert[vPlusMinus].num_adj_edges; i++)
                    {
                        j = pBNS->edge[pBNS->vert[vPlusMinus].iedge[i]].neighbor12 ^ vPlusMinus;
                        for (k = 0; k < pBNS->vert[j].num_adj_edges; k++)
                        {
                            pBNS->edge[pBNS->vert[j].iedge[k]].cap -= delta;
                        }
                    }
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &TautEdges, forbidden_edge_mask );
            }
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
                if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }

        exit_case_05:;
        }

        while (pc2i->nNumDiffMobH && pc2i->nChargeMobHRevrs > pc2i->nChargeMobHInChI)
        {
            /*----------------------------------------------------*/
            /* case 06: restored has extra H attached to -O(-)    */
            /*          while the chrge should be on C, most pro- */
            /*          bably in a small ring.ut. group or more   */
            /*          tautomeric endpoints.                     */
            /* Solution: move (-) from O to C                     */
            /*----------------------------------------------------*/
            int iO, mode;
            EdgeIndex e2;
            BNS_EDGE  *pe2;
            cur_success = 0;
            for (i = 0; !cur_success && i < pc2i->len_c2at; i++)
            {

                if (pc2i->c2at[i].nMobHRevrs == pc2i->c2at[i].nMobHInChI + 1 &&
                     pc2i->c2at[i].nNumHRevrs == pc2i->c2at[i].nMobHInChI &&
                     !pc2i->c2at[i].endptInChI && !pc2i->c2at[i].endptRevrs &&
                     at2[iO = pc2i->c2at[i].atomNumber].charge == -1 &&
                     0 <= ( e = pVA[iO].nCMinusGroupEdge - 1 ) && ( pe = pBNS->edge + e )->flow)
                {
                    /* try suitable atoms C */
                    /* first look for =C= in a small ring */
                    for (mode = 4; !cur_success && mode <= 8; mode++)
                    {

                        if (mode == 8)
                            mode = 99;

                        for (iat = 0; !cur_success && iat < pStruct->num_atoms; iat++)
                        {

                            if (!at2[iat].charge && !at2[iat].radical &&
                                 pVA[iat].cNumValenceElectrons == 4 &&
                                 0 <= ( e2 = pVA[iat].nCMinusGroupEdge - 1 ) && !( pe2 = pBNS->edge + e2 )->flow &&
                                 0 < bIsUnsatCarbonInASmallRing( at2, pVA, iat, pStruct->pbfsq, mode ))
                            {

                                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                                /* allow negative charge on the chosen carbon */
                                pe2->forbidden &= forbidden_edge_mask_inv;

                                delta = 1;
                                if (!pe->flow)
                                    continue;
                                pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                                pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );
                                pe->flow -= delta;
                                pv1->st_edge.flow -= delta;
                                pv2->st_edge.flow -= delta;
                                pBNS->tot_st_flow -= 2 * delta;

                                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                                if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                    (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warning */
                                {
                                    /* Added (-)charge to unsaturated C => nDeltaCharge == 2 */
                                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                                    if (ret > 0)
                                    {
                                        /* djb-rwth: removing redundant code */
                                        cur_success++; /* 01 */
                                        /* djb-rwth: removing redundant code */
                                    }
                                }
                                else
                                {
                                    pe->forbidden |= forbidden_edge_mask;
                                    pe->flow += delta;
                                    pv1->st_edge.flow += delta;
                                    pv2->st_edge.flow += delta;
                                    pBNS->tot_st_flow += 2 * delta;
                                }
                                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                            }
                        }
                    }
                }
            }
            if (cur_success)
            {
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
                if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
            else
            {
                break;
            }
        }
        if (pc2i->len_c2at && pc2i->nChargeMobHRevrs > pc2i->nChargeMobHInChI)
        {
            /*------------------------------------------------------------------*/
            /* case 07: -NO2 are to be tautomeric but they are not AND          */
            /*          InChI has a SINGLE tautomeric group                     */
            /*                                                                  */
            /*                   (-)O                   (-)O                    */
            /* Solution: convert     \                      \                   */
            /*                        N-X=...-Z(-)   =>      N(+)=X- ...=Z      */
            /*                      //                      /                   */
            /*                     O                    (-)O                    */
            /*                                                                  */
            /*                     O                       O                    */
            /*        or            \\                      \\                  */
            /*                        N-X=...-Z(-)    =>      N=X-  ...=Z       */
            /*                      //                       /                  */
            /*                     O                     (-)O                   */
            /*                                                                  */
            /*                                                                  */
            /*  (a) move (-) from other tautomeric atom to O in O=N-X           */
            /*          or   from other atom that has to be tautomeric          */
            /*               but is not                                         */
            /*  (b) create (+) [ion pair creation] on N as in                   */
            /*                                                                  */
            /*       OH             OH                                          */
            /*      /              /                                            */
            /*  -C=N     =>  =C-N(+)                                            */
            /*     \\             \\                                            */
            /*       O              O                                           */
            /*                                                                  */
            /*------------------------------------------------------------------*/
            int num_DB_O = 0;
            short iat_DB_O[MAX_DIFF_FIXH], iat_NO2[MAX_DIFF_FIXH];
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            /*
            AT_NUMB  *nAtno2CanonRevrs = pStruct->nAtno2Canon[0];
            */
            inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[0] &&
                                         pStruct->pOne_norm_data[0]->at ) ? pStruct->pOne_norm_data[0]->at : NULL;

            int iN, one_success;
            BNS_EDGE *peDB_O_Minus;
            int neigh, nNumO, nNumOthers;
    #define CHG_SET_WRONG_TAUT_N   0
    #define CHG_SET_WRONG_TAUT_O   1
    #define CHG_SET_WRONG_TAUT_ALL 2
    #define CHG_LAST_SET           2 /* the last index in trying */
    #define CHG_SET_O_FIXED        3
    #define CHG_SET_NUM            4
            EDGE_LIST ChangeableEdges[CHG_SET_NUM];
            memset( ChangeableEdges, 0, sizeof( ChangeableEdges ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            /* equivalent to AllocEdgeList( &EdgeList, EDGE_LIST_CLEAR ); */
            /*
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                                   pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            */
            CurrEdges.num_edges = 0; /* clear current edge list */
            cur_success = 0;
            for (i = 0; i < pStruct->num_atoms; i++)
            {
                iat = nCanon2AtnoRevrs[i];
                if ( /* orig. InChI info: taut in orig. InChI =O located in -NO2 that is not taut in Reconstructed InChI */
                     num_DB_O < MAX_DIFF_FIXH &&
                     pVA[iat].cNumValenceElectrons == 6 /* O, S, Se, Te */ &&
                     ( !at2[iat].endpoint /*|| pc2i->c2at[i].nMobHInChI*/ ) &&
                     ( e = pVA[iat].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                     at2[iat].num_H == 0 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                     /* reversed structure info: */
                     !( at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint ) /*|| pc2i->c2at[i].nMobHRevrs*/ &&
                     !at2[iat].charge &&
                     at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 &&
                     /* find whether it belongs to NO2 */
                     pVA[iN = at2[iat].neighbor[0]].cNumValenceElectrons == 5 &&
                     at2[iN].valence == 3 && ( at2[iN].charge == 0 || at2[iN].charge == 1 ) &&
                     at2[iN].chem_bonds_valence == 5 - at2[iN].charge)
                {
                    /* find the second O */
                    nNumO = nNumOthers = 0;
                    for (k = 0; k < at2[iN].valence; k++)
                    {
                        neigh = at2[iN].neighbor[k];
                        if (neigh == iat)
                        {
                            continue;
                        }
                        if (pVA[neigh].cNumValenceElectrons == 6 &&
                             !at2[neigh].endpoint &&
                             !( at_Mobile_H_Revrs && at_Mobile_H_Revrs[neigh].endpoint ) &&
                             at2[neigh].valence == 1 && at2[neigh].num_H == 0 &&
                             at2[neigh].radical == 0 && ( at2[neigh].charge == 0 || at2[neigh].charge == -1 ) &&
                             at2[neigh].chem_bonds_valence - at2[neigh].charge == 2)
                        {
                            nNumO++;
                        }
                        else
                        {
                            if (at2[iN].bond_type[k] == BOND_TYPE_SINGLE &&
                                  at2[neigh].valence > 1 &&
                                  at2[neigh].valence < at2[neigh].chem_bonds_valence)
                            {
                                nNumOthers++;
                            }
                        }
                    }
                    if (nNumO != 1 || nNumOthers != 1)
                    {
                        continue;
                    }
                    for (k = 0; k < num_DB_O; k++)
                    {
                        if (iat_NO2[k] == iN)
                        {
                            break;
                        }
                    }
                    if (k == num_DB_O)
                    {
                        iat_NO2[num_DB_O] = iN;
                        iat_DB_O[num_DB_O++] = iat;
                    }
                    /* save the =O (-)-edge to avoid interference */
                    if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_case_07;
                    }
                }
            }
            if (num_DB_O)
            {
                /* search for falsely tautomeric negatively charged atoms N and O */
                for (i = 0; i < pc2i->len_c2at; i++)
                {
                    iat = pc2i->c2at[i].atomNumber;
                    if (pc2i->c2at[i].endptRevrs && !pc2i->c2at[i].endptInChI &&
                         pc2i->c2at[i].nAtChargeRevrs == -1 &&
                         0 <= ( e = pVA[iat].nCMinusGroupEdge - 1 ) && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                         0 > FindInEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e ))
                    {
                        if (pc2i->c2at[i].nValElectr == 6)
                        {
                            if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_WRONG_TAUT_O], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_case_07;
                            }
                        }
                        else
                            if (pc2i->c2at[i].nValElectr == 5)
                            {
                                if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_WRONG_TAUT_N], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                {
                                    goto exit_case_07;
                                }
                            }
                        if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_WRONG_TAUT_ALL], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_case_07;
                        }
                    }
                }
                /* ------- finally, try to move charges from O=N --------------*/
                for (i = 0; i < num_DB_O; i++)
                {
                    int nDeltaChargeExpected;
                    one_success = 0;
                    delta = 1;
                    iat = iat_DB_O[i];
                    peDB_O_Minus = pBNS->edge + ( (long long)pVA[iat].nCMinusGroupEdge - 1 ); /* djb-rwth: cast operator added */
                    pe = pBNS->edge + pBNS->vert[iat].iedge[0];

                    if (!pe->flow)
                        continue;
                    pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                    pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                    pe->forbidden |= forbidden_edge_mask;

                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2 * delta;

                    for (k = 0; !one_success && k <= CHG_LAST_SET; k++)
                    {
                        if (!ChangeableEdges[k].num_edges)
                        {
                            continue;
                        }
                        nDeltaChargeExpected = 0;

                        SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                        RemoveForbiddenEdgeMask( pBNS, &ChangeableEdges[k], forbidden_edge_mask );
                        /* allow (-) charge to move to N=O */
                        peDB_O_Minus->forbidden &= forbidden_edge_mask_inv;

                        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                        if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                            (vPathEnd == v2 && vPathStart == v1) ) &&
                                          nDeltaCharge == nDeltaChargeExpected) /* djb-rwth: addressing LLVM warnings */
                        {
                            /* Move (-) charge to =O and remove it an endpoint => nDeltaCharge == 0 */
                            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                            if (ret > 0)
                            {
                                /* djb-rwth: removing redundant code */
                                one_success++; /* 07 */
                            }
                        }
                        INCHI_HEAPCHK
                    }
                    cur_success += one_success;

                    RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                    pe->forbidden &= forbidden_edge_mask_inv;

                    if (!one_success)
                    {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2 * delta;
                    }
                }
            }
        exit_case_07:
            for (i = 0; i < CHG_SET_NUM; i++)
            {
                AllocEdgeList( &ChangeableEdges[i], EDGE_LIST_FREE );
            }

            CurrEdges.num_edges = 0; /* clear current edge list */
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
                if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
    #undef CHG_SET_NOOH
    #undef CHG_SET_WRONG_TAUT
    #undef CHG_SET_TAUT
    #undef CHG_LAST_SET
    #undef CHG_SET_O_FIXED
    #undef CHG_SET_NUM
        }

    exit_function:
        AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );
        AllocEdgeList( &CurrEdges, EDGE_LIST_FREE );
        AllocEdgeList( &CurrEdges2, EDGE_LIST_FREE );
        AllocEdgeList( &CurrEdges3, EDGE_LIST_FREE );
        AllocEdgeList( &NFlowerEdges, EDGE_LIST_FREE );
        AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_FREE );
        AllocEdgeList( &FixedLargeRingStereoEdges, EDGE_LIST_FREE );
        AllocEdgeList( &TautEdges, EDGE_LIST_FREE );

        return ret;
    }
        */
    // END INCHI C FUNCTION: FixMobileHRestoredStructure
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FixMobileHRestoredStructure
    // INCHI✔️❌: READ_INCHI_STRING=1; COMPILE_ANSI_ONLY=1; TARGET_API_LIB=1; GCC/Linux.
    // END INCHI ACTIVE MACRO CONFIGURATION: FixMobileHRestoredStructure
    let mut ret = 0_i32;
    let forbidden_edge_mask_inv = !forbidden_edge_mask;
    let mut comparison = CMP2MHINCHI::default();
    let mut all_charge_edges = EDGE_LIST::default();
    let mut current_edges = EDGE_LIST::default();
    let mut current_edges_2 = EDGE_LIST::default();
    let mut current_edges_3 = EDGE_LIST::default();
    let mut taut_edges = EDGE_LIST::default();
    let mut nitrogen_flower_edges = EDGE_LIST::default();
    let mut other_nitrogen_flower_edges = EDGE_LIST::default();
    let mut fixed_large_ring_stereo_edges = EDGE_LIST::default();

    let _ = AllocEdgeList(heap, &mut all_charge_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut current_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut nitrogen_flower_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut current_edges_2, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut current_edges_3, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut other_nitrogen_flower_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut fixed_large_ring_stereo_edges, EDGE_LIST_CLEAR)?;
    let _ = AllocEdgeList(heap, &mut taut_edges, EDGE_LIST_CLEAR)?;

    let execution = (|| -> Result<(), SourceHeapError> {
        macro_rules! rebuild {
            () => {{
                let rebuild_result = MakeOneInChIOutOfStrFromINChI2(
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
                );
                ret = rebuild_result?;
                if ret < 0 {
                    return Ok(());
                }
            }};
        }
        macro_rules! compare {
            () => {{
                ret = FillOutExtraFixedHDataRestr(heap, pStruct)?;
                if ret != 0 {
                    return Ok(());
                }
                let at2_snapshot = heap.slice(at2.as_const())?.to_vec();
                ret = FillOutCMP2MHINCHI(
                    heap,
                    pStruct,
                    pTCGroups,
                    &at2_snapshot,
                    pVA,
                    pInChI,
                    &mut comparison,
                )?;
                if ret != 0 {
                    return Ok(());
                }
            }};
        }

        let input_fixed = heap
            .slice(pInChI[0].as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nNum_H_fixed;
        let reversed_fixed = heap
            .slice(pStruct.pOneINChI[0].as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nNum_H_fixed;
        if input_fixed.is_null() && reversed_fixed.is_null() {
            return Ok(());
        }

        let mut group_index = 0_i32;
        while group_index < pTCGroups.num_tgroups {
            let group = heap
                .slice(pTCGroups.pTCG.as_const())?
                .get(usize::try_from(group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let vertex = heap
                .slice(pBNS.vert.as_const())?
                .get(usize::try_from(group.nVertexNumber).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut order = 0_i32;
            while order < i32::from(vertex.num_adj_edges) {
                let edge = i32::from(
                    *heap
                        .slice(vertex.iedge.as_const())?
                        .get(usize::try_from(order).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                ret = AddToEdgeList(heap, &mut taut_edges, edge, INC_ADD_EDGE)?;
                if ret != 0 {
                    return Ok(());
                }
                order = order.wrapping_add(1);
            }
            group_index = group_index.wrapping_add(1);
        }

        let mut atom_number = 0_i32;
        while atom_number < pStruct.num_atoms {
            let atom_index = usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let valence = pVA.get(atom_index).ok_or(SourceHeapError::PointerOutOfBounds)?.clone();
            let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
            if minus_edge >= 0
                && heap
                    .slice(pBNS.edge.as_const())?
                    .get(usize::try_from(minus_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .forbidden
                    == 0
            {
                ret = AddToEdgeList(heap, &mut all_charge_edges, minus_edge, INC_ADD_EDGE)?;
                if ret != 0 {
                    return Ok(());
                }
            }
            let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
            if plus_edge >= 0
                && heap
                    .slice(pBNS.edge.as_const())?
                    .get(usize::try_from(plus_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .forbidden
                    == 0
            {
                ret = AddToEdgeList(heap, &mut all_charge_edges, plus_edge, INC_ADD_EDGE)?;
                if ret != 0 {
                    return Ok(());
                }
                if valence.cNumValenceElectrons == 5 && valence.cMetal == 0 {
                    let upper = GetChargeFlowerUpperEdge(heap, pBNS, pVA, plus_edge)?;
                    if upper != NO_VERTEX {
                        let flower = heap
                            .slice(pBNS.edge.as_const())?
                            .get(usize::try_from(upper).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if flower.forbidden == 0 && flower.flow != 0 {
                            ret = AddToEdgeList(heap, &mut all_charge_edges, upper, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(());
                            }
                            ret = AddToEdgeList(heap, &mut nitrogen_flower_edges, upper, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(());
                            }
                        } else {
                            ret = AddToEdgeList(heap, &mut other_nitrogen_flower_edges, upper, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(());
                            }
                        }
                    }
                }
            }
            atom_number = atom_number.wrapping_add(1);
        }

        if forbidden_stereo_edge_mask != 0 {
            let bfs = heap
                .slice(pStruct.pbfsq.as_const())?
                .first()
                .cloned()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            atom_number = 0;
            while atom_number < pStruct.num_atoms {
                let atom_index = usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let vertex = heap
                    .slice(pBNS.vert.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let mut order = 0_i32;
                while order < i32::from(atom.valence) {
                    let edge_number = i32::from(
                        *heap
                            .slice(vertex.iedge.as_const())?
                            .get(usize::try_from(order).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    if i32::from(
                        heap.slice(pBNS.edge.as_const())?
                            .get(usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .forbidden,
                    ) == forbidden_stereo_edge_mask
                    {
                        let ring = is_bond_in_Nmax_memb_ring(
                            heap,
                            at2,
                            atom_number,
                            order,
                            bfs.q,
                            bfs.nAtomLevel,
                            bfs.cSource,
                            99,
                        )?;
                        if ring > 0 {
                            ret = AddToEdgeList(heap, &mut fixed_large_ring_stereo_edges, edge_number, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(());
                            }
                        }
                    }
                    order = order.wrapping_add(1);
                }
                atom_number = atom_number.wrapping_add(1);
            }
        }

        compare!();

        if comparison.nNumTgInChI == 1
            && (comparison.nNumEndpRevrs < comparison.nNumEndpInChI || comparison.nNumTgRevrs > 1)
            && i32::from(comparison.nNumTgDBNMinusRevrs) + i32::from(comparison.nNumTgNHMinusRevrs) == 0
            && comparison.nNumTgOMinusInChI != 0
        {
            let flags = heap
                .slice(pTCGroups.pTCG.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .tg_RestoreFlags;
            if flags & TGRF_MINUS_FIRST == 0 {
                heap.slice_mut(pTCGroups.pTCG)?
                    .first_mut()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .tg_RestoreFlags |= TGRF_MINUS_FIRST;
                rebuild!();
                ret = FillOutExtraFixedHDataRestr(heap, pStruct)?;
                if ret != 0 {
                    return Ok(());
                }
                let input_fixed = heap
                    .slice(pInChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H_fixed;
                let reversed_fixed = heap
                    .slice(pStruct.pOneINChI[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .nNum_H_fixed;
                if input_fixed.is_null() && reversed_fixed.is_null() {
                    return Ok(());
                }
                let at2_snapshot = heap.slice(at2.as_const())?.to_vec();
                ret = FillOutCMP2MHINCHI(heap, pStruct, pTCGroups, &at2_snapshot, pVA, pInChI, &mut comparison)?;
                if ret != 0 || comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        if comparison.nNumTgInChI == 1
            && (comparison.nNumEndpRevrs < comparison.nNumEndpInChI || comparison.nNumTgRevrs > 1)
            && i32::from(comparison.nNumTgDBNMinusRevrs) + i32::from(comparison.nNumTgNHMinusRevrs) == 0
            && comparison.nNumTgOMinusInChI == 0
        {
            let mut double_bond_n_iii = [0_i32; MAX_DIFF_MOBH as usize];
            let mut single_bond_n_iii = [0_i32; MAX_DIFF_MOBH as usize];
            let mut number_double_bond_n_iii = 0_i32;
            let mut number_single_bond_n_iii = 0_i32;
            current_edges.num_edges = 0;
            let mut current_success = 0_i32;

            let mut atom_number = 0_i32;
            while atom_number < pStruct.num_atoms {
                let atom_index = usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let valence = pVA.get(atom_index).ok_or(SourceHeapError::PointerOutOfBounds)?;
                if valence.cNumValenceElectrons == 5
                    && valence.cPeriodicRowNumber == 1
                    && atom.endpoint == 0
                    && atom.charge == 0
                    && atom.radical == 0
                {
                    let minus_edge = valence.nCMinusGroupEdge.wrapping_sub(1);
                    let first_double_shape = number_double_bond_n_iii < MAX_DIFF_MOBH as i32
                        && atom.num_H == 0
                        && atom.valence == 2
                        && atom.chem_bonds_valence == 3
                        && atom.sb_parity[0] == 0;
                    let second_double_shape = number_double_bond_n_iii < MAX_DIFF_MOBH as i32
                        && atom.num_H == 1
                        && atom.valence == 1
                        && atom.chem_bonds_valence == 2
                        && atom.sb_parity[0] == 0;
                    let mut classified_double = false;
                    if (first_double_shape || second_double_shape) && minus_edge >= 0 {
                        let minus_edge_available = {
                            let edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(usize::try_from(minus_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden == 0 && edge.cap != 0 && edge.flow == 0
                        };
                        if minus_edge_available {
                            double_bond_n_iii[number_double_bond_n_iii as usize] = atom_number;
                            number_double_bond_n_iii = number_double_bond_n_iii.wrapping_add(1);
                            classified_double = true;
                        }
                    }
                    if !classified_double
                        && number_single_bond_n_iii < MAX_DIFF_MOBH as i32
                        && atom.num_H == 0
                        && atom.valence == 3
                        && atom.chem_bonds_valence == 3
                    {
                        let plus_edge = valence.nCPlusGroupEdge.wrapping_sub(1);
                        if plus_edge >= 0 {
                            let plus_edge_available = {
                                let edge = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(usize::try_from(plus_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.forbidden == 0 && edge.cap != 0 && edge.flow != 0
                            };
                            if plus_edge_available {
                                single_bond_n_iii[number_single_bond_n_iii as usize] = atom_number;
                                number_single_bond_n_iii = number_single_bond_n_iii.wrapping_add(1);
                                ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                        }
                    }
                }
                atom_number = atom_number.wrapping_add(1);
            }

            if number_double_bond_n_iii != 0 && number_single_bond_n_iii != 0 {
                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                let mut candidate_index = 0_i32;
                while candidate_index < number_double_bond_n_iii && current_success == 0 {
                    let atom_number = double_bond_n_iii[candidate_index as usize];
                    let atom_index = usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let atom_vertex = heap
                        .slice(pBNS.vert.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let incident_edges = heap.slice(atom_vertex.iedge.as_const())?;
                    let first_edge_number = *incident_edges.first().ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let first_has_flow = heap
                        .slice(pBNS.edge.as_const())?
                        .get(usize::try_from(first_edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .flow
                        != 0;
                    let bond_edge_number = if first_has_flow {
                        first_edge_number
                    } else {
                        let second_edge_number = *incident_edges.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let second_has_flow = heap
                            .slice(pBNS.edge.as_const())?
                            .get(usize::try_from(second_edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .flow
                            != 0;
                        if second_has_flow { second_edge_number } else { NO_VERTEX }
                    };
                    if bond_edge_number == NO_VERTEX {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    }
                    let minus_edge_number = pVA[atom_index].nCMinusGroupEdge.wrapping_sub(1);
                    let bond_edge_index =
                        usize::try_from(bond_edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let minus_edge_index =
                        usize::try_from(minus_edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let bond_edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(bond_edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if bond_edge_before.flow == 0 {
                        candidate_index = candidate_index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(bond_edge_before.neighbor1);
                    let second_vertex = i32::from(bond_edge_before.neighbor12) ^ first_vertex;
                    {
                        let edges = heap.slice_mut(pBNS.edge)?;
                        edges[bond_edge_index].forbidden =
                            (i32::from(edges[bond_edge_index].forbidden) | forbidden_edge_mask) as i8;
                        edges[minus_edge_index].forbidden =
                            (i32::from(edges[minus_edge_index].forbidden) & forbidden_edge_mask_inv) as i8;
                        edges[bond_edge_index].flow = edges[bond_edge_index].flow.wrapping_sub(1);
                    }
                    {
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_vertex, second_vertex] {
                            let vertex = vertices
                                .get_mut(
                                    usize::try_from(vertex_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(1);
                        }
                    }
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);

                    let mut path_start = 0_i32;
                    let mut path_end = 0_i32;
                    let mut path_length = 0_i32;
                    let mut delta_h = 0_i32;
                    let mut delta_charge = 0_i32;
                    let mut number_visited_atoms = 0_i32;
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
                        &mut number_visited_atoms,
                    )?;
                    if ret == 1
                        && ((path_end == first_vertex && path_start == second_vertex)
                            || (path_end == second_vertex && path_start == first_vertex))
                        && delta_charge == 2
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        if ret > 0 {
                            current_success = current_success.wrapping_add(1);
                            let minus_edge_before = heap
                                .slice(pBNS.edge.as_const())?
                                .get(minus_edge_index)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let minus_first_vertex = i32::from(minus_edge_before.neighbor1);
                            let minus_second_vertex = i32::from(minus_edge_before.neighbor12) ^ minus_first_vertex;
                            {
                                let edges = heap.slice_mut(pBNS.edge)?;
                                edges[minus_edge_index].cap = edges[minus_edge_index].cap.wrapping_sub(1);
                                edges[minus_edge_index].flow = edges[minus_edge_index].flow.wrapping_sub(1);
                            }
                            {
                                let vertices = heap.slice_mut(pBNS.vert)?;
                                for vertex_number in [minus_first_vertex, minus_second_vertex] {
                                    let vertex = vertices
                                        .get_mut(
                                            usize::try_from(vertex_number)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    vertex.st_edge.cap = vertex.st_edge.cap.wrapping_sub(1);
                                    vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(1);
                                }
                            }
                            pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);
                            pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_sub(2);
                            let structure_atom = heap
                                .slice_mut(pStruct.at)?
                                .get_mut(atom_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            structure_atom.num_H = structure_atom.num_H.wrapping_add(1);
                            pTCGroups.total_charge = pTCGroups.total_charge.wrapping_add(1);
                            pStruct.nNumRemovedProtonsByRevrs = pStruct.nNumRemovedProtonsByRevrs.wrapping_sub(1);
                        }
                    } else {
                        {
                            let edges = heap.slice_mut(pBNS.edge)?;
                            edges[bond_edge_index].forbidden =
                                (i32::from(edges[bond_edge_index].forbidden) & forbidden_edge_mask_inv) as i8;
                            edges[minus_edge_index].forbidden =
                                (i32::from(edges[minus_edge_index].forbidden) | forbidden_edge_mask) as i8;
                            edges[bond_edge_index].flow = edges[bond_edge_index].flow.wrapping_add(1);
                        }
                        {
                            let vertices = heap.slice_mut(pBNS.vert)?;
                            for vertex_number in [first_vertex, second_vertex] {
                                let vertex = vertices
                                    .get_mut(
                                        usize::try_from(vertex_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(1);
                            }
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                    }
                    candidate_index = candidate_index.wrapping_add(1);
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                current_edges.num_edges = 0;

                if current_success != 0 {
                    rebuild!();
                    compare!();
                    if comparison.bHasDifference == 0 {
                        return Ok(());
                    }
                }
            }
        }

        if comparison.nNumTgInChI == 1
            && (comparison.nNumEndpRevrs < comparison.nNumEndpInChI || comparison.nNumTgRevrs > 1)
            && comparison.nNumTgMInChI == 0
            && comparison.nNumTgNInChI != 0
            && comparison.nNumTgOInChI != 0
        {
            let mut action = 0_i32;
            current_edges.num_edges = 0;
            current_edges_2.num_edges = 0;
            let mut current_success = 0_i32;
            let mut taut_group_index = 0_i32;
            while taut_group_index < pTCGroups.num_tgroups && current_success == 0 {
                let group = heap
                    .slice(pTCGroups.pTCG.as_const())?
                    .get(usize::try_from(taut_group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let group_vertex = heap
                    .slice(pBNS.vert.as_const())?
                    .get(usize::try_from(group.nVertexNumber).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let mut group_edge_order = 0_i32;
                while group_edge_order < i32::from(group_vertex.num_adj_edges) && current_success == 0 {
                    let taut_edge_number = i32::from(
                        *heap
                            .slice(group_vertex.iedge.as_const())?
                            .get(usize::try_from(group_edge_order).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                    );
                    let taut_edge = heap
                        .slice(pBNS.edge.as_const())?
                        .get(usize::try_from(taut_edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let endpoint = i32::from(taut_edge.neighbor1);
                    let mut second_taut_edge = -1_i32;
                    let mut center_edge_1 = -1_i32;
                    let mut center_edge_2 = -1_i32;
                    let mut second_endpoint = NO_VERTEX;
                    if pVA
                        .get(usize::try_from(endpoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .cNumValenceElectrons
                        == 6
                        && taut_edge.cap != 0
                    {
                        let endpoint_vertex = heap
                            .slice(pBNS.vert.as_const())?
                            .get(usize::try_from(endpoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let endpoint_atom = heap
                            .slice(at2.as_const())?
                            .get(usize::try_from(endpoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let mut neighbor_order = 0_i32;
                        while neighbor_order < i32::from(endpoint_atom.valence) && second_taut_edge < 0 {
                            center_edge_1 = i32::from(
                                *heap
                                    .slice(endpoint_vertex.iedge.as_const())?
                                    .get(
                                        usize::try_from(neighbor_order)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            );
                            let first_center_edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(usize::try_from(center_edge_1).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let center = i32::from(first_center_edge.neighbor12) ^ endpoint;
                            let center_atom = heap
                                .slice(at2.as_const())?
                                .get(usize::try_from(center).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if center_atom.endpoint != 0
                                || first_center_edge.cap == 0
                                || i32::from(first_center_edge.flow) + i32::from(taut_edge.cap == taut_edge.flow) != 1
                            {
                                neighbor_order = neighbor_order.wrapping_add(1);
                                continue;
                            }
                            let center_vertex = heap
                                .slice(pBNS.vert.as_const())?
                                .get(usize::try_from(center).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let mut center_order = 0_i32;
                            while center_order < i32::from(center_atom.valence) {
                                center_edge_2 = i32::from(
                                    *heap
                                        .slice(center_vertex.iedge.as_const())?
                                        .get(
                                            usize::try_from(center_order)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                );
                                let second_center_edge = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(center_edge_2)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .cloned()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                second_endpoint = i32::from(second_center_edge.neighbor12) ^ center;
                                let second_index = usize::try_from(second_endpoint)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                if second_center_edge.cap != 0
                                    && i32::from(second_center_edge.flow) + i32::from(first_center_edge.flow) == 1
                                    && i32::from(heap.slice(at2.as_const())?[second_index].endpoint)
                                        == taut_group_index.wrapping_add(1)
                                    && pVA[second_index].cNumValenceElectrons == 5
                                {
                                    let candidate_taut_edge = pVA[second_index].nTautGroupEdge.wrapping_sub(1);
                                    if candidate_taut_edge >= 0 {
                                        let edge = heap
                                            .slice(pBNS.edge.as_const())?
                                            .get(
                                                usize::try_from(candidate_taut_edge)
                                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                            )
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                        if edge.forbidden == 0
                                            && i32::from(second_center_edge.flow) + i32::from(edge.cap == edge.flow)
                                                == 1
                                        {
                                            second_taut_edge = candidate_taut_edge;
                                            break;
                                        }
                                    }
                                }
                                center_order = center_order.wrapping_add(1);
                            }
                            neighbor_order = neighbor_order.wrapping_add(1);
                        }
                    }
                    if second_taut_edge >= 0 {
                        let taut_index =
                            usize::try_from(taut_edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let second_taut_index =
                            usize::try_from(second_taut_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let center_1_index =
                            usize::try_from(center_edge_1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let center_2_index =
                            usize::try_from(center_edge_2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let before = heap.slice(pBNS.edge.as_const())?;
                        let first_taut = before[taut_index].clone();
                        let second_taut = before[second_taut_index].clone();
                        let first_center = before[center_1_index].clone();
                        let second_center = before[center_2_index].clone();
                        if i32::from(first_taut.cap) - i32::from(first_taut.flow) == 0
                            && i32::from(second_taut.cap) - i32::from(second_taut.flow) == 1
                            && first_center.flow == 0
                            && second_center.flow == 1
                        {
                            let edges = heap.slice_mut(pBNS.edge)?;
                            edges[taut_index].flow = edges[taut_index].flow.wrapping_sub(1);
                            edges[second_taut_index].flow = edges[second_taut_index].flow.wrapping_add(1);
                            edges[center_2_index].flow = edges[center_2_index].flow.wrapping_sub(1);
                            edges[center_1_index].flow = edges[center_1_index].flow.wrapping_add(1);
                            action |= 1;
                        }
                        let after = heap.slice(pBNS.edge.as_const())?;
                        if i32::from(after[taut_index].cap) - i32::from(after[taut_index].flow) == 1
                            && i32::from(after[second_taut_index].cap) - i32::from(after[second_taut_index].flow) == 0
                            && after[center_1_index].flow == 1
                            && after[center_2_index].flow == 0
                        {
                            let group_index =
                                usize::try_from(taut_group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let groups = heap.slice_mut(pTCGroups.pTCG)?;
                            groups[group_index].tg_num_H = groups[group_index].tg_num_H.wrapping_sub(1);
                            groups[group_index].tg_num_Minus = groups[group_index].tg_num_Minus.wrapping_add(1);
                            groups[group_index].tg_RestoreFlags |= TGRF_MINUS_FIRST;
                            groups[group_index].tg_set_Minus = second_endpoint.wrapping_add(1);
                            let taut_groups = heap.slice_mut(pStruct.ti.t_group)?;
                            taut_groups[group_index].num[1] = taut_groups[group_index].num[1].wrapping_add(1);
                            pTCGroups.total_charge = pTCGroups.total_charge.wrapping_sub(1);
                            pTCGroups.tgroup_charge = pTCGroups.tgroup_charge.wrapping_sub(1);
                            pStruct.nNumRemovedProtonsByRevrs = pStruct.nNumRemovedProtonsByRevrs.wrapping_add(1);
                            action |= 2;
                            current_success = current_success.wrapping_add(1);
                        }
                    }
                    group_edge_order = group_edge_order.wrapping_add(1);
                }
                taut_group_index = taut_group_index.wrapping_add(1);
            }

            if i32::from(comparison.nNumTgNHInChI) + i32::from(comparison.nNumTgNH2InChI) == 0
                && comparison.nNumTgOHInChI != 0
                && current_success == 0
            {
                taut_group_index = 0;
                while taut_group_index < pTCGroups.num_tgroups {
                    let group = heap
                        .slice(pTCGroups.pTCG.as_const())?
                        .get(usize::try_from(taut_group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let group_vertex = heap
                        .slice(pBNS.vert.as_const())?
                        .get(usize::try_from(group.nVertexNumber).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut order = 0_i32;
                    while order < i32::from(group_vertex.num_adj_edges) {
                        let taut_edge_number = i32::from(
                            *heap
                                .slice(group_vertex.iedge.as_const())?
                                .get(usize::try_from(order).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        );
                        let taut_edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(usize::try_from(taut_edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let endpoint = usize::from(taut_edge.neighbor1);
                        let endpoint_atom = heap
                            .slice(at2.as_const())?
                            .get(endpoint)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if pVA[endpoint].cNumValenceElectrons == 6
                            && endpoint_atom.valence == endpoint_atom.chem_bonds_valence
                            && taut_edge.flow != 0
                            && taut_edge.flow == taut_edge.cap
                        {
                            ret = AddToEdgeList(heap, &mut current_edges, taut_edge_number, INC_ADD_EDGE)?;
                            if ret != 0 {
                                return Ok(());
                            }
                        } else if pVA[endpoint].cNumValenceElectrons == 5
                            && pVA[endpoint].cPeriodicRowNumber == 1
                            && i32::from(endpoint_atom.valence).wrapping_add(1)
                                == i32::from(endpoint_atom.chem_bonds_valence)
                            && taut_edge.cap != 0
                            && taut_edge.flow.wrapping_add(1) == taut_edge.cap
                        {
                            let endpoint_vertex = heap
                                .slice(pBNS.vert.as_const())?
                                .get(endpoint)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let mut found_edge = -1_i32;
                            let mut neighbor_order = 0_i32;
                            while neighbor_order < i32::from(endpoint_atom.valence) && found_edge < 0 {
                                let first_edge_number = i32::from(
                                    *heap
                                        .slice(endpoint_vertex.iedge.as_const())?
                                        .get(
                                            usize::try_from(neighbor_order)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                );
                                let first_edge = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(first_edge_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .cloned()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if first_edge.flow == 1 {
                                    let center = i32::from(first_edge.neighbor12)
                                        ^ i32::try_from(endpoint).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                    let center_index =
                                        usize::try_from(center).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                    if heap.slice(at2.as_const())?[center_index].endpoint != 0 {
                                        neighbor_order = neighbor_order.wrapping_add(1);
                                        continue;
                                    }
                                    let center_atom = heap.slice(at2.as_const())?[center_index].clone();
                                    let center_vertex = heap.slice(pBNS.vert.as_const())?[center_index].clone();
                                    let mut center_order = 0_i32;
                                    while center_order < i32::from(center_atom.valence) {
                                        let second_edge_number = i32::from(
                                            *heap
                                                .slice(center_vertex.iedge.as_const())?
                                                .get(
                                                    usize::try_from(center_order)
                                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                                )
                                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                        );
                                        let second_edge = heap
                                            .slice(pBNS.edge.as_const())?
                                            .get(
                                                usize::try_from(second_edge_number)
                                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                            )
                                            .cloned()
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                        let second_endpoint = i32::from(second_edge.neighbor12) ^ center;
                                        let second_index = usize::try_from(second_endpoint)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                        if second_edge.flow == 0
                                            && i32::from(heap.slice(at2.as_const())?[second_index].endpoint)
                                                == taut_group_index.wrapping_add(1)
                                            && pVA[second_index].cNumValenceElectrons == 6
                                        {
                                            let candidate = pVA[second_index].nTautGroupEdge.wrapping_sub(1);
                                            if candidate > 0 {
                                                let edge = heap
                                                    .slice(pBNS.edge.as_const())?
                                                    .get(
                                                        usize::try_from(candidate)
                                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                                    )
                                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                                if edge.forbidden == 0 && edge.flow != 0 {
                                                    found_edge = candidate;
                                                    ret = AddToEdgeList(
                                                        heap,
                                                        &mut current_edges_2,
                                                        candidate,
                                                        INC_ADD_EDGE,
                                                    )?;
                                                    if ret != 0 {
                                                        return Ok(());
                                                    }
                                                    break;
                                                }
                                            }
                                        }
                                        center_order = center_order.wrapping_add(1);
                                    }
                                }
                                neighbor_order = neighbor_order.wrapping_add(1);
                            }
                            if found_edge < 0 {
                                ret = AddToEdgeList(heap, &mut current_edges, taut_edge_number, INC_ADD_EDGE)?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                        }
                        order = order.wrapping_add(1);
                    }
                    taut_group_index = taut_group_index.wrapping_add(1);
                }

                SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                SetForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                SetForbiddenEdgeMask(heap, pBNS, &current_edges_2, forbidden_edge_mask)?;
                let selected = if current_edges_2.num_edges != 0 {
                    &current_edges_2
                } else {
                    &current_edges
                };
                let selected_edges = if selected.num_edges > 0 {
                    heap.slice(selected.pnEdges.as_const())?
                        [..usize::try_from(selected.num_edges).map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                        .to_vec()
                } else {
                    Vec::new()
                };
                for edge_number in selected_edges {
                    if current_success != 0 {
                        break;
                    }
                    let edge_index = usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if edge_before.flow == 0 {
                        continue;
                    }
                    let first_vertex = i32::from(edge_before.neighbor1);
                    let second_vertex = i32::from(edge_before.neighbor12) ^ first_vertex;
                    let mut found_group = None;
                    let mut index = 0_i32;
                    while index < pTCGroups.num_tgroups {
                        if heap
                            .slice(pTCGroups.pTCG.as_const())?
                            .get(usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nVertexNumber
                            == second_vertex
                        {
                            found_group = Some(index);
                            break;
                        }
                        index = index.wrapping_add(1);
                    }
                    let Some(found_group) = found_group else {
                        continue;
                    };
                    heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow.wrapping_sub(1);
                    {
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_vertex, second_vertex] {
                            let vertex = vertices
                                .get_mut(
                                    usize::try_from(vertex_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(1);
                        }
                    }
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);
                    let mut path_start = 0_i32;
                    let mut path_end = 0_i32;
                    let mut path_length = 0_i32;
                    let mut delta_h = 0_i32;
                    let mut delta_charge = 0_i32;
                    let mut visited = 0_i32;
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
                        &mut visited,
                    )?;
                    if ret == 1
                        && ((path_end == first_vertex && path_start == second_vertex)
                            || (path_end == second_vertex && path_start == first_vertex))
                        && delta_charge == 2
                    {
                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                        if ret > 0 {
                            current_success = current_success.wrapping_add(1);
                            let group_index =
                                usize::try_from(found_group).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let groups = heap.slice_mut(pTCGroups.pTCG)?;
                            groups[group_index].tg_num_H = groups[group_index].tg_num_H.wrapping_sub(1);
                            groups[group_index].tg_num_Minus = groups[group_index].tg_num_Minus.wrapping_add(1);
                            groups[group_index].tg_RestoreFlags |= TGRF_MINUS_FIRST;
                            let taut_groups = heap.slice_mut(pStruct.ti.t_group)?;
                            let taut_group = taut_groups
                                .get_mut(group_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            taut_group.num[1] = taut_group.num[1].wrapping_add(1);
                            pTCGroups.total_charge = pTCGroups.total_charge.wrapping_sub(1);
                            pTCGroups.tgroup_charge = pTCGroups.tgroup_charge.wrapping_sub(1);
                            pStruct.nNumRemovedProtonsByRevrs = pStruct.nNumRemovedProtonsByRevrs.wrapping_add(1);
                            action |= 4;
                        }
                    } else {
                        heap.slice_mut(pBNS.edge)?[edge_index].flow = edge_before.flow;
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_vertex, second_vertex] {
                            let vertex = vertices
                                .get_mut(
                                    usize::try_from(vertex_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(1);
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                    }
                }
                RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
            } else if i32::from(comparison.nNumTgNHInChI) + i32::from(comparison.nNumTgNH2InChI) != 0
                && comparison.nNumTgOInChI != 0
                && current_success == 0
            {
                taut_group_index = 0;
                while taut_group_index < pTCGroups.num_tgroups && current_success == 0 {
                    let group_index =
                        usize::try_from(taut_group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let group = heap
                        .slice(pTCGroups.pTCG.as_const())?
                        .get(group_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let vertex = heap
                        .slice(pBNS.vert.as_const())?
                        .get(usize::try_from(group.nVertexNumber).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut order = 0_i32;
                    while order < i32::from(vertex.num_adj_edges) && current_success == 0 {
                        let edge_number = i32::from(
                            *heap
                                .slice(vertex.iedge.as_const())?
                                .get(usize::try_from(order).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        );
                        let edge = heap
                            .slice(pBNS.edge.as_const())?
                            .get(usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let endpoint = usize::from(edge.neighbor1);
                        let atom = heap
                            .slice(at2.as_const())?
                            .get(endpoint)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if pVA[endpoint].cNumValenceElectrons == 5
                            && pVA[endpoint].cPeriodicRowNumber == 1
                            && atom.valence == atom.chem_bonds_valence
                            && edge.flow != 0
                            && edge.flow == edge.cap
                        {
                            current_success = current_success.wrapping_add(1);
                            let groups = heap.slice_mut(pTCGroups.pTCG)?;
                            groups[group_index].tg_num_H = groups[group_index].tg_num_H.wrapping_sub(1);
                            groups[group_index].tg_num_Minus = groups[group_index].tg_num_Minus.wrapping_add(1);
                            groups[group_index].tg_RestoreFlags |= TGRF_MINUS_FIRST;
                            groups[group_index].tg_set_Minus = i32::try_from(endpoint)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?
                                .wrapping_add(1);
                            let taut_groups = heap.slice_mut(pStruct.ti.t_group)?;
                            let taut_group = taut_groups
                                .get_mut(group_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            taut_group.num[1] = taut_group.num[1].wrapping_add(1);
                            pTCGroups.total_charge = pTCGroups.total_charge.wrapping_sub(1);
                            pTCGroups.tgroup_charge = pTCGroups.tgroup_charge.wrapping_sub(1);
                            pStruct.nNumRemovedProtonsByRevrs = pStruct.nNumRemovedProtonsByRevrs.wrapping_add(1);
                            action |= 8;
                        }
                        order = order.wrapping_add(1);
                    }
                    taut_group_index = taut_group_index.wrapping_add(1);
                }
            }

            if current_success != 0 {
                rebuild!();
                compare!();
                let reversed_minus = if pStruct.One_ti.num_t_groups == 1 {
                    heap.slice(pStruct.One_ti.t_group.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .num[1]
                } else {
                    0
                };
                if pStruct.One_ti.num_t_groups == 1 && reversed_minus != 0 {
                    let group_index = 0_usize;
                    if action & (8 | 2) != 0 {
                        let groups = heap.slice_mut(pTCGroups.pTCG)?;
                        groups[group_index].tg_num_H = groups[group_index].tg_num_H.wrapping_add(1);
                        groups[group_index].tg_num_Minus = groups[group_index].tg_num_Minus.wrapping_sub(1);
                        groups[group_index].tg_RestoreFlags &= !TGRF_MINUS_FIRST;
                        groups[group_index].tg_set_Minus = 0;
                    } else if action & 4 != 0 {
                        let groups = heap.slice_mut(pTCGroups.pTCG)?;
                        groups[group_index].tg_num_H = groups[group_index].tg_num_H.wrapping_add(1);
                        groups[group_index].tg_num_Minus = groups[group_index].tg_num_Minus.wrapping_sub(1);
                        groups[group_index].tg_RestoreFlags &= !TGRF_MINUS_FIRST;
                    } else {
                        ret = RI_ERR_PROGR;
                        return Ok(());
                    }
                    let taut_group = heap
                        .slice_mut(pStruct.ti.t_group)?
                        .get_mut(group_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    taut_group.num[1] = taut_group.num[1].wrapping_sub(1);
                    pTCGroups.total_charge = pTCGroups.total_charge.wrapping_add(1);
                    pTCGroups.tgroup_charge = pTCGroups.tgroup_charge.wrapping_add(1);
                    pStruct.nNumRemovedProtonsByRevrs = pStruct.nNumRemovedProtonsByRevrs.wrapping_sub(1);
                    current_success = current_success.wrapping_sub(1);
                    rebuild!();
                    compare!();
                }
                if comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        if comparison.nNumTgInChI == 1
            && (comparison.nNumEndpRevrs < comparison.nNumEndpInChI || comparison.nNumTgRevrs > 1)
            && comparison.nNumTgMInChI == 0
            && (comparison.nNumTgNInChI != 0 || comparison.nNumTgOInChI != 0)
        {
            let plus_minus_vertex = GetPlusMinusVertex(heap, pBNS, pTCGroups, 1, 1)?;
            if plus_minus_vertex != NO_VERTEX {
                'case_05: {
                    current_edges.num_edges = 0;
                    current_edges_2.num_edges = 0;
                    let mut current_success = 0_i32;
                    let canonical_to_atom = pStruct.nCanon2Atno[0];
                    let mut canonical_number = 0_i32;
                    while canonical_number < pStruct.num_atoms {
                        let canonical_index =
                            usize::try_from(canonical_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let atom_number = i32::from(
                            *heap
                                .slice(canonical_to_atom.as_const())?
                                .get(canonical_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?,
                        );
                        let atom_index =
                            usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let atom = heap
                            .slice(at2.as_const())?
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if atom.endpoint == 0
                            && atom.charge == 0
                            && atom.radical == 0
                            && atom.num_H == 0
                            && pVA[canonical_index].cNumValenceElectrons == 5
                            && pVA[canonical_index].cPeriodicRowNumber == 1
                        {
                            let plus_edge = pVA[atom_index].nCPlusGroupEdge.wrapping_sub(1);
                            if plus_edge >= 0 {
                                let edge = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(usize::try_from(plus_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                if edge.flow != 0 && edge.forbidden == 0 {
                                    ret = AddToEdgeList(heap, &mut current_edges, plus_edge, INC_ADD_EDGE)?;
                                    if ret != 0 {
                                        return Ok(());
                                    }
                                }
                            }
                        }
                        canonical_number = canonical_number.wrapping_add(1);
                    }
                    if current_edges.num_edges == 0 {
                        break 'case_05;
                    }

                    let mut taut_group_index = 0_i32;
                    while taut_group_index < pTCGroups.num_tgroups && current_success == 0 {
                        current_edges_2.num_edges = 0;
                        let group_index =
                            usize::try_from(taut_group_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let group = heap
                            .slice(pTCGroups.pTCG.as_const())?
                            .get(group_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let group_vertex_number = group.nVertexNumber;
                        let group_vertex = heap
                            .slice(pBNS.vert.as_const())?
                            .get(
                                usize::try_from(group_vertex_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let mut order = 0_i32;
                        while order < i32::from(group_vertex.num_adj_edges) {
                            let edge_number = i32::from(
                                *heap
                                    .slice(group_vertex.iedge.as_const())?
                                    .get(usize::try_from(order).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            );
                            let edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let endpoint = usize::from(edge.neighbor1);
                            let atom = heap
                                .slice(at2.as_const())?
                                .get(endpoint)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if atom.charge == 0 && atom.radical == 0 && i32::from(edge.cap) - i32::from(edge.flow) == 1
                            {
                                ret = AddToEdgeList(heap, &mut current_edges_2, edge_number, INC_ADD_EDGE)?;
                                if ret != 0 {
                                    return Ok(());
                                }
                            }
                            order = order.wrapping_add(1);
                        }
                        if current_edges_2.num_edges == 0 {
                            break 'case_05;
                        }
                        SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        SetForbiddenEdgeMask(heap, pBNS, &taut_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(heap, pBNS, &current_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(heap, pBNS, &current_edges_2, forbidden_edge_mask)?;

                        {
                            let vertices = heap.slice_mut(pBNS.vert)?;
                            let plus = vertices
                                .get_mut(
                                    usize::try_from(plus_minus_vertex)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            plus.st_edge.cap = plus.st_edge.cap.wrapping_add(1);
                            let taut = vertices
                                .get_mut(
                                    usize::try_from(group_vertex_number)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            taut.st_edge.cap = taut.st_edge.cap.wrapping_add(1);
                        }
                        pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_add(2);

                        let plus_vertex = heap
                            .slice(pBNS.vert.as_const())?
                            .get(usize::try_from(plus_minus_vertex).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let mut plus_order = 0_i32;
                        while plus_order < i32::from(plus_vertex.num_adj_edges) {
                            let edge_number = i32::from(
                                *heap
                                    .slice(plus_vertex.iedge.as_const())?
                                    .get(usize::try_from(plus_order).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            );
                            let edge = heap
                                .slice(pBNS.edge.as_const())?
                                .get(usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let intermediate = i32::from(edge.neighbor12) ^ plus_minus_vertex;
                            let intermediate_vertex = heap
                                .slice(pBNS.vert.as_const())?
                                .get(usize::try_from(intermediate).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                .cloned()
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            let edge_numbers = heap.slice(intermediate_vertex.iedge.as_const())?
                                [..usize::from(intermediate_vertex.num_adj_edges)]
                                .to_vec();
                            for edge_number in edge_numbers {
                                let index =
                                    usize::try_from(edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                let edge = heap
                                    .slice_mut(pBNS.edge)?
                                    .get_mut(index)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.cap = edge.cap.wrapping_add(1);
                            }
                            plus_order = plus_order.wrapping_add(1);
                        }

                        let mut path_start = 0_i32;
                        let mut path_end = 0_i32;
                        let mut path_length = 0_i32;
                        let mut delta_h = 0_i32;
                        let mut delta_charge = 0_i32;
                        let mut visited = 0_i32;
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
                            &mut visited,
                        )?;
                        if ret == 1
                            && ((path_end == plus_minus_vertex && path_start == group_vertex_number)
                                || (path_end == group_vertex_number && path_start == plus_minus_vertex))
                            && delta_charge == 1
                        {
                            ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                            if ret > 0 {
                                current_success = current_success.wrapping_add(1);
                                pTCGroups.total_charge = pTCGroups.total_charge.wrapping_add(1);
                                let groups = heap.slice_mut(pTCGroups.pTCG)?;
                                groups[group_index].edges_cap = groups[group_index].edges_cap.wrapping_add(1);
                                groups[group_index].tg_num_H = groups[group_index].tg_num_H.wrapping_add(1);
                                pStruct.nNumRemovedProtonsByRevrs = pStruct.nNumRemovedProtonsByRevrs.wrapping_sub(1);
                            }
                        } else {
                            {
                                let vertices = heap.slice_mut(pBNS.vert)?;
                                let plus_index = usize::try_from(plus_minus_vertex)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                vertices[plus_index].st_edge.cap = vertices[plus_index].st_edge.cap.wrapping_sub(1);
                                let group_vertex_index = usize::try_from(group_vertex_number)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                vertices[group_vertex_index].st_edge.cap =
                                    vertices[group_vertex_index].st_edge.cap.wrapping_sub(1);
                            }
                            pBNS.tot_st_cap = pBNS.tot_st_cap.wrapping_sub(2);
                            plus_order = 0;
                            while plus_order < i32::from(plus_vertex.num_adj_edges) {
                                let edge_number = i32::from(
                                    *heap
                                        .slice(plus_vertex.iedge.as_const())?
                                        .get(
                                            usize::try_from(plus_order)
                                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                        )
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                );
                                let edge = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(
                                        usize::try_from(edge_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .cloned()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let intermediate = i32::from(edge.neighbor12) ^ plus_minus_vertex;
                                let vertex = heap
                                    .slice(pBNS.vert.as_const())?
                                    .get(
                                        usize::try_from(intermediate)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .cloned()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let edge_numbers =
                                    heap.slice(vertex.iedge.as_const())?[..usize::from(vertex.num_adj_edges)].to_vec();
                                for edge_number in edge_numbers {
                                    let index = usize::try_from(edge_number)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                    let edge = heap
                                        .slice_mut(pBNS.edge)?
                                        .get_mut(index)
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    edge.cap = edge.cap.wrapping_sub(1);
                                }
                                plus_order = plus_order.wrapping_add(1);
                            }
                        }
                        RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(heap, pBNS, &taut_edges, forbidden_edge_mask)?;
                        taut_group_index = taut_group_index.wrapping_add(1);
                    }
                    if current_success != 0 {
                        rebuild!();
                        compare!();
                        if comparison.bHasDifference == 0 {
                            return Ok(());
                        }
                    }
                }
            }
        }

        while comparison.nNumDiffMobH != 0 && comparison.nChargeMobHRevrs > comparison.nChargeMobHInChI {
            let mut current_success = 0_i32;
            let mut difference_index = 0_i32;
            while current_success == 0 && difference_index < i32::from(comparison.len_c2at) {
                let difference = comparison.c2at
                    [usize::try_from(difference_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                .clone();
                if i32::from(difference.nMobHRevrs) != i32::from(difference.nMobHInChI).wrapping_add(1)
                    || difference.nNumHRevrs != difference.nMobHInChI
                    || difference.endptInChI != 0
                    || difference.endptRevrs != 0
                {
                    difference_index = difference_index.wrapping_add(1);
                    continue;
                }
                let oxygen_number = i32::from(difference.atomNumber);
                let oxygen_index = usize::try_from(oxygen_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let oxygen = heap
                    .slice(at2.as_const())?
                    .get(oxygen_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if oxygen.charge == -1 {
                    let oxygen_minus_edge = pVA
                        .get(oxygen_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let oxygen_edge_has_flow = oxygen_minus_edge >= 0
                        && heap
                            .slice(pBNS.edge.as_const())?
                            .get(usize::try_from(oxygen_minus_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .flow
                            != 0;
                    if oxygen_edge_has_flow {
                        let mut mode = 4_i32;
                        while current_success == 0 && mode <= 8 {
                            if mode == 8 {
                                mode = 99;
                            }
                            let mut atom_number = 0_i32;
                            while current_success == 0 && atom_number < pStruct.num_atoms {
                                let atom_index =
                                    usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                let atom = heap
                                    .slice(at2.as_const())?
                                    .get(atom_index)
                                    .cloned()
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                let mut carbon_minus_edge = -1_i32;
                                let candidate_matches = if atom.charge == 0
                                    && atom.radical == 0
                                    && pVA[atom_index].cNumValenceElectrons == 4
                                {
                                    carbon_minus_edge = pVA[atom_index].nCMinusGroupEdge.wrapping_sub(1);
                                    if carbon_minus_edge >= 0 {
                                        let carbon_edge_no_flow = heap
                                            .slice(pBNS.edge.as_const())?
                                            .get(
                                                usize::try_from(carbon_minus_edge)
                                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                            )
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                                            .flow
                                            == 0;
                                        if carbon_edge_no_flow {
                                            let bfs = heap
                                                .slice(pStruct.pbfsq.as_const())?
                                                .first()
                                                .cloned()
                                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                            bIsUnsatCarbonInASmallRing(heap, at2, pVA, atom_number, &bfs, mode)? > 0
                                        } else {
                                            false
                                        }
                                    } else {
                                        false
                                    }
                                } else {
                                    false
                                };
                                if candidate_matches {
                                    SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                                    {
                                        let carbon_edge = heap
                                            .slice_mut(pBNS.edge)?
                                            .get_mut(
                                                usize::try_from(carbon_minus_edge)
                                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                            )
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                        carbon_edge.forbidden =
                                            (i32::from(carbon_edge.forbidden) & forbidden_edge_mask_inv) as i8;
                                    }
                                    let oxygen_edge_index = usize::try_from(oxygen_minus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                                    let oxygen_edge_before = heap
                                        .slice(pBNS.edge.as_const())?
                                        .get(oxygen_edge_index)
                                        .cloned()
                                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                    if oxygen_edge_before.flow == 0 {
                                        atom_number = atom_number.wrapping_add(1);
                                        continue;
                                    }
                                    let first_vertex = i32::from(oxygen_edge_before.neighbor1);
                                    let second_vertex = i32::from(oxygen_edge_before.neighbor12) ^ first_vertex;
                                    heap.slice_mut(pBNS.edge)?[oxygen_edge_index].flow =
                                        oxygen_edge_before.flow.wrapping_sub(1);
                                    {
                                        let vertices = heap.slice_mut(pBNS.vert)?;
                                        for vertex_number in [first_vertex, second_vertex] {
                                            let vertex = vertices
                                                .get_mut(
                                                    usize::try_from(vertex_number)
                                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                                )
                                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(1);
                                        }
                                    }
                                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);

                                    let mut path_start = 0_i32;
                                    let mut path_end = 0_i32;
                                    let mut path_length = 0_i32;
                                    let mut delta_h = 0_i32;
                                    let mut delta_charge = 0_i32;
                                    let mut number_visited_atoms = 0_i32;
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
                                        &mut number_visited_atoms,
                                    )?;
                                    if ret == 1
                                        && ((path_end == first_vertex && path_start == second_vertex)
                                            || (path_end == second_vertex && path_start == first_vertex))
                                        && delta_charge == 1
                                    {
                                        ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                                        if ret > 0 {
                                            current_success = current_success.wrapping_add(1);
                                        }
                                    } else {
                                        {
                                            let oxygen_edge = heap
                                                .slice_mut(pBNS.edge)?
                                                .get_mut(oxygen_edge_index)
                                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                            oxygen_edge.forbidden =
                                                (i32::from(oxygen_edge.forbidden) | forbidden_edge_mask) as i8;
                                            oxygen_edge.flow = oxygen_edge.flow.wrapping_add(1);
                                        }
                                        {
                                            let vertices = heap.slice_mut(pBNS.vert)?;
                                            for vertex_number in [first_vertex, second_vertex] {
                                                let vertex = vertices
                                                    .get_mut(
                                                        usize::try_from(vertex_number)
                                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                                    )
                                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(1);
                                            }
                                        }
                                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                                    }
                                    SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                                }
                                atom_number = atom_number.wrapping_add(1);
                            }
                            mode = mode.wrapping_add(1);
                        }
                    }
                }
                difference_index = difference_index.wrapping_add(1);
            }
            if current_success != 0 {
                rebuild!();
                compare!();
                if comparison.bHasDifference == 0 {
                    return Ok(());
                }
            } else {
                break;
            }
        }

        if comparison.len_c2at != 0 && comparison.nChargeMobHRevrs > comparison.nChargeMobHInChI {
            const CHG_SET_WRONG_TAUT_N: usize = 0;
            const CHG_SET_WRONG_TAUT_O: usize = 1;
            const CHG_SET_WRONG_TAUT_ALL: usize = 2;
            const CHG_LAST_SET: usize = 2;
            const CHG_SET_O_FIXED: usize = 3;
            const CHG_SET_NUM: usize = 4;

            let mut double_bond_oxygen = [0_i16; MAX_DIFF_FIXH as usize];
            let mut nitrogen_dioxide = [0_i16; MAX_DIFF_FIXH as usize];
            let mut number_double_bond_oxygen = 0_i32;
            let canonical_to_atom = pStruct.nCanon2Atno[0];
            let mobile_h_reversed = if pStruct.pOne_norm_data[0].is_null() {
                SourceMutPointer::null()
            } else {
                heap.slice(pStruct.pOne_norm_data[0].as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .at
            };
            let mut changeable_edges: [EDGE_LIST; CHG_SET_NUM] = std::array::from_fn(|_| EDGE_LIST::default());
            current_edges.num_edges = 0;
            let mut current_success = 0_i32;
            let mut leave_case = false;

            let mut canonical_number = 0_i32;
            while canonical_number < pStruct.num_atoms && !leave_case {
                let canonical_index =
                    usize::try_from(canonical_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom_number = i32::from(
                    *heap
                        .slice(canonical_to_atom.as_const())?
                        .get(canonical_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let atom_index = usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice(at2.as_const())?
                    .get(atom_index)
                    .cloned()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let mut minus_edge = -1_i32;
                let mut dioxide_center = None;
                if number_double_bond_oxygen < MAX_DIFF_FIXH as i32
                    && pVA[atom_index].cNumValenceElectrons == 6
                    && atom.endpoint == 0
                    && {
                        minus_edge = pVA[atom_index].nCMinusGroupEdge.wrapping_sub(1);
                        minus_edge >= 0
                            && heap
                                .slice(pBNS.edge.as_const())?
                                .get(usize::try_from(minus_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .forbidden
                                == 0
                    }
                    && atom.num_H == 0
                    && !(!mobile_h_reversed.is_null()
                        && heap
                            .slice(mobile_h_reversed.as_const())?
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .endpoint
                            != 0)
                    && atom.charge == 0
                    && atom.valence == 1
                    && atom.chem_bonds_valence == 2
                {
                    let nitrogen_index = usize::from(atom.neighbor[0]);
                    let nitrogen = heap
                        .slice(at2.as_const())?
                        .get(nitrogen_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if pVA
                        .get(nitrogen_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .cNumValenceElectrons
                        == 5
                        && nitrogen.valence == 3
                        && (nitrogen.charge == 0 || nitrogen.charge == 1)
                        && i32::from(nitrogen.chem_bonds_valence) == 5_i32.wrapping_sub(i32::from(nitrogen.charge))
                    {
                        dioxide_center = Some(nitrogen_index);
                    }
                }
                if let Some(nitrogen_index) = dioxide_center {
                    let nitrogen = heap
                        .slice(at2.as_const())?
                        .get(nitrogen_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut number_oxygen = 0_i32;
                    let mut number_others = 0_i32;
                    let mut neighbor_order = 0_i32;
                    while neighbor_order < i32::from(nitrogen.valence) {
                        let position =
                            usize::try_from(neighbor_order).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let neighbor_number = i32::from(nitrogen.neighbor[position]);
                        if neighbor_number == atom_number {
                            neighbor_order = neighbor_order.wrapping_add(1);
                            continue;
                        }
                        let neighbor_index =
                            usize::try_from(neighbor_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let neighbor = heap
                            .slice(at2.as_const())?
                            .get(neighbor_index)
                            .cloned()
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let neighbor_reversed_endpoint = !mobile_h_reversed.is_null()
                            && heap
                                .slice(mobile_h_reversed.as_const())?
                                .get(neighbor_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .endpoint
                                != 0;
                        if pVA
                            .get(neighbor_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .cNumValenceElectrons
                            == 6
                            && neighbor.endpoint == 0
                            && !neighbor_reversed_endpoint
                            && neighbor.valence == 1
                            && neighbor.num_H == 0
                            && neighbor.radical == 0
                            && (neighbor.charge == 0 || neighbor.charge == -1)
                            && i32::from(neighbor.chem_bonds_valence).wrapping_sub(i32::from(neighbor.charge)) == 2
                        {
                            number_oxygen = number_oxygen.wrapping_add(1);
                        } else if nitrogen.bond_type[position] == BOND_TYPE_SINGLE as u8
                            && neighbor.valence > 1
                            && neighbor.valence < neighbor.chem_bonds_valence
                        {
                            number_others = number_others.wrapping_add(1);
                        }
                        neighbor_order = neighbor_order.wrapping_add(1);
                    }
                    if number_oxygen == 1 && number_others == 1 {
                        let nitrogen_number =
                            i32::try_from(nitrogen_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let mut duplicate_index = 0_i32;
                        while duplicate_index < number_double_bond_oxygen
                            && i32::from(nitrogen_dioxide[duplicate_index as usize]) != nitrogen_number
                        {
                            duplicate_index = duplicate_index.wrapping_add(1);
                        }
                        if duplicate_index == number_double_bond_oxygen {
                            nitrogen_dioxide[number_double_bond_oxygen as usize] = nitrogen_number as i16;
                            double_bond_oxygen[number_double_bond_oxygen as usize] = atom_number as i16;
                            number_double_bond_oxygen = number_double_bond_oxygen.wrapping_add(1);
                        }
                        ret = AddToEdgeList(heap, &mut changeable_edges[CHG_SET_O_FIXED], minus_edge, INC_ADD_EDGE)?;
                        if ret != 0 {
                            leave_case = true;
                        }
                    }
                }
                canonical_number = canonical_number.wrapping_add(1);
            }

            if number_double_bond_oxygen != 0 && !leave_case {
                let mut difference_index = 0_i32;
                while difference_index < i32::from(comparison.len_c2at) && !leave_case {
                    let difference = comparison.c2at
                        [usize::try_from(difference_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?]
                    .clone();
                    if difference.endptRevrs != 0 && difference.endptInChI == 0 && difference.nAtChargeRevrs == -1 {
                        let atom_number = i32::from(difference.atomNumber);
                        let atom_index =
                            usize::try_from(atom_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let minus_edge = pVA
                            .get(atom_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .nCMinusGroupEdge
                            .wrapping_sub(1);
                        let charge_edge_ok = minus_edge >= 0
                            && {
                                let edge = heap
                                    .slice(pBNS.edge.as_const())?
                                    .get(usize::try_from(minus_edge).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                edge.forbidden == 0 && edge.flow != 0
                            }
                            && FindInEdgeList(heap, &changeable_edges[CHG_SET_O_FIXED], minus_edge)? < 0;
                        if charge_edge_ok && difference.nValElectr == 6 {
                            ret = AddToEdgeList(
                                heap,
                                &mut changeable_edges[CHG_SET_WRONG_TAUT_O],
                                minus_edge,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                leave_case = true;
                            }
                        } else if charge_edge_ok && difference.nValElectr == 5 {
                            ret = AddToEdgeList(
                                heap,
                                &mut changeable_edges[CHG_SET_WRONG_TAUT_N],
                                minus_edge,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                leave_case = true;
                            }
                        }
                        if charge_edge_ok && !leave_case {
                            ret = AddToEdgeList(
                                heap,
                                &mut changeable_edges[CHG_SET_WRONG_TAUT_ALL],
                                minus_edge,
                                INC_ADD_EDGE,
                            )?;
                            if ret != 0 {
                                leave_case = true;
                            }
                        }
                    }
                    difference_index = difference_index.wrapping_add(1);
                }

                let mut target_index = 0_i32;
                while target_index < number_double_bond_oxygen && !leave_case {
                    let atom_index = usize::try_from(i32::from(double_bond_oxygen[target_index as usize]))
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let target_minus_edge = pVA
                        .get(atom_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .nCMinusGroupEdge
                        .wrapping_sub(1);
                    let atom_vertex = heap
                        .slice(pBNS.vert.as_const())?
                        .get(atom_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let bond_edge_number = *heap
                        .slice(atom_vertex.iedge.as_const())?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let bond_edge_index =
                        usize::try_from(bond_edge_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let bond_edge_before = heap
                        .slice(pBNS.edge.as_const())?
                        .get(bond_edge_index)
                        .cloned()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if bond_edge_before.flow == 0 {
                        target_index = target_index.wrapping_add(1);
                        continue;
                    }
                    let first_vertex = i32::from(bond_edge_before.neighbor1);
                    let second_vertex = i32::from(bond_edge_before.neighbor12) ^ first_vertex;
                    {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(bond_edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) | forbidden_edge_mask) as i8;
                        edge.flow = edge.flow.wrapping_sub(1);
                    }
                    {
                        let vertices = heap.slice_mut(pBNS.vert)?;
                        for vertex_number in [first_vertex, second_vertex] {
                            let vertex = vertices
                                .get_mut(
                                    usize::try_from(vertex_number).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            vertex.st_edge.flow = vertex.st_edge.flow.wrapping_sub(1);
                        }
                    }
                    pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_sub(2);

                    let mut one_success = 0_i32;
                    let mut set_index = 0_usize;
                    while one_success == 0 && set_index <= CHG_LAST_SET {
                        if changeable_edges[set_index].num_edges == 0 {
                            set_index += 1;
                            continue;
                        }
                        let expected_delta_charge = 0_i32;
                        SetForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                        RemoveForbiddenEdgeMask(heap, pBNS, &changeable_edges[set_index], forbidden_edge_mask)?;
                        {
                            let edge = heap
                                .slice_mut(pBNS.edge)?
                                .get_mut(
                                    usize::try_from(target_minus_edge)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            edge.forbidden = (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                        }
                        let mut path_start = 0_i32;
                        let mut path_end = 0_i32;
                        let mut path_length = 0_i32;
                        let mut delta_h = 0_i32;
                        let mut delta_charge = 0_i32;
                        let mut number_visited_atoms = 0_i32;
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
                            &mut number_visited_atoms,
                        )?;
                        if ret == 1
                            && ((path_end == first_vertex && path_start == second_vertex)
                                || (path_end == second_vertex && path_start == first_vertex))
                            && delta_charge == expected_delta_charge
                        {
                            ret = RunBnsRestoreOnce(heap, pBNS, pBD, pVA, pTCGroups, clock_result)?;
                            if ret > 0 {
                                one_success = one_success.wrapping_add(1);
                            }
                        }
                        set_index += 1;
                    }
                    current_success = current_success.wrapping_add(one_success);
                    RemoveForbiddenEdgeMask(heap, pBNS, &all_charge_edges, forbidden_edge_mask)?;
                    {
                        let edge = heap
                            .slice_mut(pBNS.edge)?
                            .get_mut(bond_edge_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        edge.forbidden = (i32::from(edge.forbidden) & forbidden_edge_mask_inv) as i8;
                    }
                    if one_success == 0 {
                        heap.slice_mut(pBNS.edge)?[bond_edge_index].flow = bond_edge_before.flow;
                        {
                            let vertices = heap.slice_mut(pBNS.vert)?;
                            for vertex_number in [first_vertex, second_vertex] {
                                let vertex = vertices
                                    .get_mut(
                                        usize::try_from(vertex_number)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                                vertex.st_edge.flow = vertex.st_edge.flow.wrapping_add(1);
                            }
                        }
                        pBNS.tot_st_flow = pBNS.tot_st_flow.wrapping_add(2);
                    }
                    target_index = target_index.wrapping_add(1);
                }
            }

            for edges in &mut changeable_edges {
                let _ = AllocEdgeList(heap, edges, EDGE_LIST_FREE)?;
            }
            current_edges.num_edges = 0;
            if current_success != 0 {
                rebuild!();
                compare!();
                if comparison.bHasDifference == 0 {
                    return Ok(());
                }
            }
        }

        Ok(())
    })();

    let cleanup_result = (|| -> Result<(), SourceHeapError> {
        let _ = AllocEdgeList(heap, &mut all_charge_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut current_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut current_edges_2, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut current_edges_3, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut nitrogen_flower_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut other_nitrogen_flower_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut fixed_large_ring_stereo_edges, EDGE_LIST_FREE)?;
        let _ = AllocEdgeList(heap, &mut taut_edges, EDGE_LIST_FREE)?;
        Ok(())
    })();

    execution?;
    cleanup_result?;
    Ok(ret)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{
        BNS_EDGE, BNS_VERTEX, INChI_Aux, RI_ERR_ALLOC, SourceMutPointer, T_GROUP, T_GROUP_INFO, TC_GROUP,
    };

    fn call_fix_mobile_h_restored_structure(
        heap: &mut SourceHeap,
        structure: &mut StrFromINChI,
        atoms: SourceMutPointer<inp_ATOM>,
        bns: &mut BN_STRUCT,
        valence: &mut [VAL_AT],
        groups: &mut ALL_TC_GROUPS,
        input: [SourceMutPointer<INChI>; 2],
        runs: &mut i32,
        total_delta: &mut i32,
    ) -> Result<i32, SourceHeapError> {
        FixMobileHRestoredStructure(
            heap,
            &mut CANON_GLOBALS::default(),
            SourceMutPointer::null(),
            &INPUT_PARMS::default(),
            &STRUCT_DATA::default(),
            bns,
            &mut BN_DATA::default(),
            structure,
            atoms,
            atoms,
            atoms,
            valence,
            groups,
            None,
            None,
            None,
            input,
            i64::MIN,
            i32::MAX,
            Some(runs),
            Some(total_delta),
            0x20,
            0,
            clock_t::default(),
        )
    }

    fn fixed_h_inchi(fixed_h: SourceMutPointer<i8>, mobile_h: SourceMutPointer<i8>, total_charge: i32) -> INChI {
        INChI {
            nNum_H_fixed: fixed_h,
            nNum_H: mobile_h,
            nTotalCharge: total_charge,
            ..INChI::default()
        }
    }

    fn mapping_aux(heap: &mut SourceHeap, atom_numbers: Vec<u16>) -> SourceMutPointer<INChI_Aux> {
        let canonical = heap.allocate_model_storage(atom_numbers).unwrap();
        heap.allocate_model_storage(vec![INChI_Aux {
            nOrigAtNosInCanonOrd: canonical,
            ..INChI_Aux::default()
        }])
        .unwrap()
    }

    fn group_info(heap: &mut SourceHeap, groups: Vec<T_GROUP>, endpoints: Vec<u16>) -> T_GROUP_INFO {
        T_GROUP_INFO {
            num_t_groups: groups.len() as i32,
            max_num_t_groups: groups.len() as i32,
            nNumEndpoints: endpoints.len() as i32,
            t_group: heap.allocate_model_storage(groups).unwrap(),
            nEndpointAtomNumber: heap.allocate_model_storage(endpoints).unwrap(),
            ..T_GROUP_INFO::default()
        }
    }

    #[test]
    fn source_port__ichirvr5__getplusminusvertex__line_58() {
        let mut heap = SourceHeap::default();
        let edges = heap
            .allocate_model_storage(vec![
                BNS_EDGE::default(),
                BNS_EDGE {
                    neighbor12: 14,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    neighbor12: 5,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    neighbor12: 14,
                    forbidden: -1,
                    ..BNS_EDGE::default()
                },
                BNS_EDGE {
                    neighbor12: 5,
                    forbidden: 1,
                    ..BNS_EDGE::default()
                },
            ])
            .unwrap();
        let tcg = heap
            .allocate_model_storage(vec![
                TC_GROUP {
                    nForwardEdge: 1,
                    nVertexNumber: 10,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    nForwardEdge: 2,
                    nVertexNumber: 11,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    nForwardEdge: 0,
                    nVertexNumber: 10,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    nForwardEdge: 1,
                    nVertexNumber: 9,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    nForwardEdge: 3,
                    nVertexNumber: 10,
                    ..TC_GROUP::default()
                },
                TC_GROUP {
                    nForwardEdge: 4,
                    nVertexNumber: 11,
                    ..TC_GROUP::default()
                },
            ])
            .unwrap();
        let bns = BN_STRUCT {
            num_atoms: 10,
            edge: edges,
            ..BN_STRUCT::default()
        };

        let mut groups = ALL_TC_GROUPS {
            pTCG: tcg,
            nGroup: [-1; 18],
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &groups, 0, 0), Ok(NO_VERTEX));

        groups.nGroup[TCG_Plus as usize] = 0;
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &groups, 0, 0), Ok(4));
        groups.nGroup[TCG_Minus as usize] = 1;
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &groups, 0, 0), Ok(4));

        groups.nGroup[TCG_Plus as usize] = -1;
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &groups, 0, 0), Ok(7));
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &groups, 1, 0), Ok(NO_VERTEX));

        groups.nGroup[TCG_Plus as usize] = 0;
        groups.nGroup[TCG_Minus as usize] = -1;
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &groups, 0, 1), Ok(NO_VERTEX));

        groups.nGroup[TCG_Plus as usize] = 4;
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &groups, 0, 0), Ok(4));
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &groups, -1, 0), Ok(NO_VERTEX));

        groups.nGroup[TCG_Plus as usize] = -1;
        groups.nGroup[TCG_Minus as usize] = 5;
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &groups, 0, 0), Ok(1));
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &groups, 0, 2), Ok(NO_VERTEX));

        groups.nGroup[TCG_Minus as usize] = -1;
        for rejected_group in [2, 3] {
            groups.nGroup[TCG_Plus as usize] = rejected_group;
            assert_eq!(
                GetPlusMinusVertex(&heap, &bns, &groups, 0, 0),
                Ok(NO_VERTEX),
                "group={rejected_group}"
            );
        }

        let null_groups = ALL_TC_GROUPS {
            pTCG: SourceMutPointer::null(),
            nGroup: [-1; 18],
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(GetPlusMinusVertex(&heap, &bns, &null_groups, 0, 0), Ok(NO_VERTEX));
    }

    #[test]
    fn source_port__ichirvr5__bisunsatcarboninasmallring__line_92() {
        fn atom(valence: i8, chemical_valence: i8, neighbors: &[u16]) -> inp_ATOM {
            let mut atom = inp_ATOM {
                valence,
                chem_bonds_valence: chemical_valence,
                ..inp_ATOM::default()
            };
            atom.neighbor[..neighbors.len()].copy_from_slice(neighbors);
            atom
        }

        fn evaluate_fast(atom: inp_ATOM, ring_size: i8, minimum: i32) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(vec![atom]).unwrap();
            bIsUnsatCarbonInASmallRing(
                &mut heap,
                atoms,
                &[VAL_AT {
                    cMinRingSize: ring_size,
                    ..VAL_AT::default()
                }],
                0,
                &BFS_Q::default(),
                minimum,
            )
        }

        for minimum in [i32::MIN, 4] {
            assert_eq!(evaluate_fast(atom(2, 4, &[]), 5, minimum), Ok(1));
        }
        assert_eq!(evaluate_fast(atom(1, 4, &[]), 5, 4), Ok(0));
        assert_eq!(evaluate_fast(atom(2, 3, &[]), 5, 4), Ok(0));
        assert_eq!(evaluate_fast(atom(2, 4, &[]), 6, 4), Ok(0));
        assert_eq!(evaluate_fast(atom(2, 4, &[]), -1, 4), Ok(1));

        for ring_size in [-1, 5] {
            assert_eq!(evaluate_fast(atom(2, 3, &[]), ring_size, 5), Ok(1));
        }
        assert_eq!(evaluate_fast(atom(2, 2, &[]), 0, 5), Ok(0));
        assert_eq!(evaluate_fast(atom(2, 2, &[]), 6, 5), Ok(0));
        assert_eq!(evaluate_fast(atom(4, 5, &[]), 3, 5), Ok(0));
        assert_eq!(evaluate_fast(atom(2, 2, &[]), 3, 5), Ok(0));

        let mut ring_heap = SourceHeap::default();
        let ring_atoms = ring_heap
            .allocate_model_storage(vec![atom(2, 3, &[1, 2]), atom(2, 2, &[0, 2]), atom(2, 2, &[0, 1])])
            .unwrap();
        let levels = ring_heap.allocate_model_storage(vec![0_u16; 3]).unwrap();
        let sources = ring_heap.allocate_model_storage(vec![0_i8; 3]).unwrap();
        let queue_values = ring_heap.allocate_model_storage(vec![0_u16; 3]).unwrap();
        let queue = ring_heap
            .allocate_model_storage(vec![crate::source_types::QUEUE {
                Val: queue_values,
                nTotLength: 3,
                ..crate::source_types::QUEUE::default()
            }])
            .unwrap();
        let bfs = BFS_Q {
            q: queue,
            nAtomLevel: levels,
            cSource: sources,
            num_at: 3,
            ..BFS_Q::default()
        };
        assert_eq!(
            bIsUnsatCarbonInASmallRing(&mut ring_heap, ring_atoms, &vec![VAL_AT::default(); 3], 0, &bfs, 5,),
            Ok(1)
        );
        assert_eq!(ring_heap.slice(levels.as_const()).unwrap(), [0, 0, 0]);
        assert_eq!(ring_heap.slice(sources.as_const()).unwrap(), [0, 0, 0]);

        let mut open_heap = SourceHeap::default();
        let open_atoms = open_heap
            .allocate_model_storage(vec![atom(2, 3, &[1, 2]), atom(1, 1, &[0]), atom(1, 1, &[0])])
            .unwrap();
        let levels = open_heap.allocate_model_storage(vec![0_u16; 3]).unwrap();
        let sources = open_heap.allocate_model_storage(vec![0_i8; 3]).unwrap();
        let values = open_heap.allocate_model_storage(vec![0_u16; 3]).unwrap();
        let queue = open_heap
            .allocate_model_storage(vec![crate::source_types::QUEUE {
                Val: values,
                nTotLength: 3,
                ..crate::source_types::QUEUE::default()
            }])
            .unwrap();
        assert_eq!(
            bIsUnsatCarbonInASmallRing(
                &mut open_heap,
                open_atoms,
                &vec![VAL_AT::default(); 3],
                0,
                &BFS_Q {
                    q: queue,
                    nAtomLevel: levels,
                    cSource: sources,
                    num_at: 3,
                    ..BFS_Q::default()
                },
                5,
            ),
            Ok(0)
        );
    }

    #[test]
    fn source_port__ichirvr5__fixmobilehrestoredstructure__line_141() {
        let mut no_fixed_heap = SourceHeap::default();
        let input = no_fixed_heap.allocate_model_storage(vec![INChI::default()]).unwrap();
        let reversed = no_fixed_heap.allocate_model_storage(vec![INChI::default()]).unwrap();
        let mut no_fixed_structure = StrFromINChI {
            pOneINChI: [reversed, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        let no_fixed_allocations = no_fixed_heap.live_allocation_count();
        let mut runs = i32::MIN;
        let mut total_delta = i32::MAX;
        assert_eq!(
            call_fix_mobile_h_restored_structure(
                &mut no_fixed_heap,
                &mut no_fixed_structure,
                SourceMutPointer::null(),
                &mut BN_STRUCT::default(),
                &mut [],
                &mut ALL_TC_GROUPS::default(),
                [input, SourceMutPointer::null()],
                &mut runs,
                &mut total_delta,
            ),
            Ok(0)
        );
        assert_eq!((runs, total_delta), (i32::MIN, i32::MAX));
        assert_eq!(no_fixed_heap.live_allocation_count(), no_fixed_allocations);

        for successful_allocations in [0_u64, 1] {
            let mut heap = SourceHeap::default();
            let fixed_h = heap.allocate_model_storage(vec![0_i8]).unwrap();
            let input = heap
                .allocate_model_storage(vec![fixed_h_inchi(fixed_h, SourceMutPointer::null(), 0)])
                .unwrap();
            let reversed = heap
                .allocate_model_storage(vec![fixed_h_inchi(fixed_h, SourceMutPointer::null(), 0)])
                .unwrap();
            let auxiliary = mapping_aux(&mut heap, vec![1]);
            let atoms = heap.allocate_model_storage(Vec::<inp_ATOM>::new()).unwrap();
            let mut structure = StrFromINChI {
                num_atoms: 0,
                pOneINChI: [reversed, SourceMutPointer::null()],
                pOneINChI_Aux: [auxiliary, SourceMutPointer::null()],
                ..StrFromINChI::default()
            };
            let fixture_allocations = heap.live_allocation_count();
            heap.fail_after_allocations(successful_allocations);
            let mut runs = 17;
            let mut total_delta = -19;
            assert_eq!(
                call_fix_mobile_h_restored_structure(
                    &mut heap,
                    &mut structure,
                    atoms,
                    &mut BN_STRUCT::default(),
                    &mut [],
                    &mut ALL_TC_GROUPS::default(),
                    [input, SourceMutPointer::null()],
                    &mut runs,
                    &mut total_delta,
                ),
                Ok(RI_ERR_ALLOC),
                "allocation ordinal={successful_allocations}"
            );
            assert_eq!((runs, total_delta), (17, -19));
            assert_eq!(
                heap.live_allocation_count(),
                fixture_allocations + successful_allocations as usize
            );
            assert_eq!(structure.nCanon2Atno[0].is_null(), successful_allocations == 0);
            assert!(structure.nAtno2Canon[0].is_null());
        }

        for fail_first_edge_list_growth in [true, false] {
            let mut heap = SourceHeap::default();
            let fixed_h = heap.allocate_model_storage(vec![0_i8]).unwrap();
            let input = heap
                .allocate_model_storage(vec![fixed_h_inchi(fixed_h, SourceMutPointer::null(), 0)])
                .unwrap();
            let reversed = heap
                .allocate_model_storage(vec![fixed_h_inchi(fixed_h, SourceMutPointer::null(), 0)])
                .unwrap();
            let auxiliary = mapping_aux(&mut heap, vec![1]);
            let atoms = heap.allocate_model_storage(vec![inp_ATOM::default()]).unwrap();
            let incident = heap.allocate_model_storage(vec![0_i32]).unwrap();
            let vertices = heap
                .allocate_model_storage(vec![BNS_VERTEX {
                    num_adj_edges: 1,
                    iedge: incident,
                    ..BNS_VERTEX::default()
                }])
                .unwrap();
            let edges = heap.allocate_model_storage(vec![BNS_EDGE::default()]).unwrap();
            let tcg = heap
                .allocate_model_storage(vec![TC_GROUP {
                    nVertexNumber: 0,
                    ..TC_GROUP::default()
                }])
                .unwrap();
            let mut groups = ALL_TC_GROUPS {
                num_tgroups: 1,
                pTCG: tcg,
                ..ALL_TC_GROUPS::default()
            };
            let mut structure = StrFromINChI {
                num_atoms: 1,
                pOneINChI: [reversed, SourceMutPointer::null()],
                pOneINChI_Aux: [auxiliary, SourceMutPointer::null()],
                ..StrFromINChI::default()
            };
            let mut bns = BN_STRUCT {
                num_atoms: 1,
                num_vertices: 1,
                num_edges: 1,
                vert: vertices,
                edge: edges,
                ..BN_STRUCT::default()
            };
            let mut valence = vec![VAL_AT {
                nCMinusGroupEdge: 1,
                ..VAL_AT::default()
            }];
            let fixture_allocations = heap.live_allocation_count();
            if fail_first_edge_list_growth {
                heap.fail_after_allocations(0);
            } else {
                heap.trace_source_allocations();
            }
            let mut runs = 23;
            let mut total_delta = 29;
            assert_eq!(
                call_fix_mobile_h_restored_structure(
                    &mut heap,
                    &mut structure,
                    atoms,
                    &mut bns,
                    &mut valence,
                    &mut groups,
                    [input, SourceMutPointer::null()],
                    &mut runs,
                    &mut total_delta,
                ),
                Ok(if fail_first_edge_list_growth { RI_ERR_ALLOC } else { 0 })
            );
            assert_eq!((runs, total_delta), (23, 29));
            assert_eq!(heap.slice(edges.as_const()).unwrap(), &[BNS_EDGE::default()]);
            assert_eq!(
                heap.live_allocation_count(),
                fixture_allocations + usize::from(!fail_first_edge_list_growth) * 2
            );
            if !fail_first_edge_list_growth {
                assert_eq!(heap.source_allocation_calls(), 4);
                assert!(!structure.nCanon2Atno[0].is_null());
                assert!(!structure.nAtno2Canon[0].is_null());
            }
        }

        let mut branch_heap = SourceHeap::default();
        let fixed_h = branch_heap.allocate_model_storage(vec![0_i8; 2]).unwrap();
        let input_mobile = branch_heap.allocate_model_storage(vec![0_i8; 2]).unwrap();
        let reversed_mobile = branch_heap.allocate_model_storage(vec![2_i8; 2]).unwrap();
        let input = branch_heap
            .allocate_model_storage(vec![fixed_h_inchi(fixed_h, input_mobile, -1)])
            .unwrap();
        let reversed = branch_heap
            .allocate_model_storage(vec![fixed_h_inchi(fixed_h, reversed_mobile, 1)])
            .unwrap();
        let auxiliary = mapping_aux(&mut branch_heap, vec![1, 2]);
        let atoms = branch_heap
            .allocate_model_storage(vec![inp_ATOM::default(), inp_ATOM::default()])
            .unwrap();
        let original_groups = group_info(
            &mut branch_heap,
            vec![T_GROUP {
                nNumEndpoints: 2,
                nFirstEndpointAtNoPos: 0,
                ..T_GROUP::default()
            }],
            vec![0, 1],
        );
        let reversed_groups = group_info(
            &mut branch_heap,
            vec![
                T_GROUP {
                    nNumEndpoints: 2,
                    nFirstEndpointAtNoPos: 0,
                    ..T_GROUP::default()
                },
                T_GROUP::default(),
            ],
            vec![0, 1],
        );
        let mut structure = StrFromINChI {
            num_atoms: 2,
            ti: original_groups,
            One_ti: reversed_groups,
            pOneINChI: [reversed, SourceMutPointer::null()],
            pOneINChI_Aux: [auxiliary, SourceMutPointer::null()],
            ..StrFromINChI::default()
        };
        let mut valence = vec![
            VAL_AT {
                cNumValenceElectrons: 5,
                cPeriodicRowNumber: 1,
                ..VAL_AT::default()
            },
            VAL_AT {
                cNumValenceElectrons: 6,
                ..VAL_AT::default()
            },
        ];
        let mut runs = 31;
        let mut total_delta = 37;
        let mut branch_groups = ALL_TC_GROUPS {
            nGroup: [-1; 18],
            ..ALL_TC_GROUPS::default()
        };
        assert_eq!(
            call_fix_mobile_h_restored_structure(
                &mut branch_heap,
                &mut structure,
                atoms,
                &mut BN_STRUCT::default(),
                &mut valence,
                &mut branch_groups,
                [input, SourceMutPointer::null()],
                &mut runs,
                &mut total_delta,
            ),
            Ok(0)
        );
        assert_eq!((runs, total_delta), (31, 37));
        assert_eq!(
            branch_heap.slice(atoms.as_const()).unwrap(),
            &[inp_ATOM::default(), inp_ATOM::default()]
        );
        assert_eq!(structure.nNumRemovedProtonsByRevrs, 0);
    }
}
